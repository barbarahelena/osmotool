"""
quantifier.py — read/hit counting and RPKM/RPM normalisation for osmotool

Counting rules
--------------
* **Best-hit per read**: each query read is assigned to the single
  highest-scoring (best bitscore) subject.  Ties are broken by lowest
  evalue, then by lexicographic subject ID order for reproducibility.
* **Paired-end weighting**: each read in a pair contributes 0.5 counts,
  so one read *pair* = 1 fragment count.  Single-end reads contribute 1.0.
* **Gene family roll-up**: subjects are mapped to their gene-family label
  (first '|'-delimited field of the subject ID).  Final counts are
  summed per family.

Normalisation
-------------
* **RPM** (profile mode — no gene lengths available):
    RPM = (raw_count × 10⁶) / total_reads
* **RPKM** (annotate mode — protein lengths from Prodigal available):
    RPKM = (raw_count × 10⁹) / (gene_length_bp × total_reads)
  gene_length_bp = length_aa × 3  (coding sequence approximation)
"""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Iterable

from osmotool.utils import gene_family_from_header

log = logging.getLogger("osmotool")

# Families present in osmodiamond. Only cosmetic: controls which families
# are guaranteed to appear as an explicit 0-count row in gene_counts.tsv
# even with zero hits -- a family that gets >=1 hit shows up regardless,
# via count_hits()'s defaultdict. Keep in sync with osmo_refdb's
# families.yaml when adding a new gene family (that's the actual source
# of truth for what the reference database contains).
KNOWN_FAMILIES = (
    "ectA", "ectB", "ectC", "betL", "kdpA", "nhaA",
    "proX", "proP", "otsA", "otsB", "mscL", "mscS",
    "galE", "mazG", "murB",
)


# ---------------------------------------------------------------------------
# Best-hit selection
# ---------------------------------------------------------------------------

def select_best_hits(alignments: Iterable[dict]) -> dict[str, dict]:
    """
    Return {qseqid: best_alignment_row} keeping the top-scoring hit per read.

    Tie-breaking: higher bitscore → lower evalue → lexicographically smallest sseqid.
    """
    best: dict[str, dict] = {}
    for row in alignments:
        qid = row["qseqid"]
        if qid not in best:
            best[qid] = row
        else:
            prev = best[qid]
            if _is_better(row, prev):
                best[qid] = row
    return best


def _is_better(candidate: dict, current: dict) -> bool:
    """Return True if *candidate* is a better alignment than *current*."""
    if candidate["bitscore"] > current["bitscore"]:
        return True
    if candidate["bitscore"] == current["bitscore"]:
        if candidate["evalue"] < current["evalue"]:
            return True
        if candidate["evalue"] == current["evalue"]:
            return candidate["sseqid"] < current["sseqid"]
    return False


def filter_by_family_cutoff(
    alignments: Iterable[dict], cutoffs: dict[str, float]
) -> Iterable[dict]:
    """
    Drop alignment rows whose bitscore falls below *their own family's*
    calibrated cutoff (see functions.load_diamond_cutoffs), applied to raw
    alignments BEFORE select_best_hits(). This matters: a query's top hit
    by raw bitscore might be to family X but still below X's cutoff (e.g.
    betT scoring highest against betL, but not high enough to be
    confidently betL rather than a related-but-distinct paralog) while a
    weaker but still-legitimate hit to family Y clears Y's own cutoff --
    filtering the raw rows first lets that second-best hit become the
    query's best surviving call, instead of the query being dropped
    entirely just because its single top hit failed one family's cutoff.

    A family absent from *cutoffs* (no calibration data available, or
    *cutoffs* is empty/not supplied) passes through ungated.
    """
    for row in alignments:
        family = gene_family_from_header(row["sseqid"])
        if family in cutoffs and row["bitscore"] < cutoffs[family]:
            continue
        yield row


def apply_cascade_check(
    best_hits: dict[str, dict],
    cascade_scores: dict[str, dict[str, float]],
    cutoffs: dict[str, float],
) -> dict[str, dict]:
    """
    `profile` mode's DIAMOND+HMM cascade: drop a DIAMOND best-hit if its
    family was flagged in osmo_refdb's profile_cascade.tsv (real,
    demonstrated precision problems, e.g. betL vs betT/CaiT) AND the read
    fails a short-read-calibrated HMM check for that same family.

    Parameters
    ----------
    best_hits:
        DIAMOND's best hit per read (from select_best_hits()).
    cascade_scores:
        {read_id: {family: best_raw_hmm_bitscore}} from running hmmscan
        (use_cut_ga=False) on just the subset of reads whose DIAMOND
        family needed checking, then grouping parse_hmmscan_output() rows
        by read_id (via functions.strip_orf_suffix) and family.
    cutoffs:
        {family: short_read_bitscore_threshold} from
        functions.load_cascade_config.

    A family absent from *cutoffs* is never checked (DIAMOND's call is
    trusted as-is, no evidence it needs a second opinion). A read whose
    family IS in cutoffs but has no recorded HMM score at all against that
    family (i.e. absent from cascade_scores, or that family missing from
    its per-family score dict) is treated as scoring 0 -- hmmscan found
    no signal at all for it, which should never pass a real threshold.
    """
    kept: dict[str, dict] = {}
    for qid, row in best_hits.items():
        family = gene_family_from_header(row["sseqid"])
        if family not in cutoffs:
            kept[qid] = row
            continue
        score = cascade_scores.get(qid, {}).get(family, 0.0)
        if score >= cutoffs[family]:
            kept[qid] = row
    return kept


# ---------------------------------------------------------------------------
# Raw counting
# ---------------------------------------------------------------------------

def count_hits(
    best_hits: dict[str, dict],
    *,
    paired: bool = False,
) -> dict[str, float]:
    """
    Aggregate best hits into per-gene-family raw counts.

    Parameters
    ----------
    best_hits:
        Mapping of qseqid → alignment row (from :func:`select_best_hits`).
    paired:
        If True each read contributes 0.5 (one read-pair = one fragment).

    Returns
    -------
    dict[family_label, raw_count_float]
        All known families are present; families with no hits have count 0.
    """
    weight = 0.5 if paired else 1.0
    counts: dict[str, float] = defaultdict(float)

    # Pre-populate so all families appear in output
    for fam in KNOWN_FAMILIES:
        counts[fam] = 0.0

    for row in best_hits.values():
        family = gene_family_from_header(row["sseqid"])
        counts[family] += weight

    return dict(counts)


# ---------------------------------------------------------------------------
# RPM / RPKM normalisation
# ---------------------------------------------------------------------------

def compute_rpm(
    counts: dict[str, float],
    total_reads: int,
) -> dict[str, float]:
    """
    RPM = (raw_count × 10⁶) / total_reads

    Used in ``profile`` mode where gene lengths are unavailable.
    """
    if total_reads == 0:
        return {fam: 0.0 for fam in counts}
    return {fam: (cnt * 1_000_000) / total_reads for fam, cnt in counts.items()}


def compute_rpkm(
    counts: dict[str, float],
    total_reads: int,
    gene_lengths_aa: dict[str, float],
) -> dict[str, float]:
    """
    RPKM = (raw_count × 10⁹) / (gene_length_bp × total_reads)

    gene_length_bp is estimated as the mean protein length (aa) × 3 for
    each family, derived from the proteins that received hits.

    Parameters
    ----------
    counts:
        Per-family raw counts.
    total_reads:
        Total input reads (denominator for normalisation).
    gene_lengths_aa:
        Mapping of family_label → mean gene length in amino acids.
        Families absent from this dict fall back to RPM normalisation.
    """
    if total_reads == 0:
        return {fam: 0.0 for fam in counts}

    rpkm: dict[str, float] = {}
    for fam, cnt in counts.items():
        if fam in gene_lengths_aa and gene_lengths_aa[fam] > 0:
            gene_len_bp = gene_lengths_aa[fam] * 3
            rpkm[fam] = (cnt * 1_000_000_000) / (gene_len_bp * total_reads)
        else:
            # Fallback: report RPM when length is unknown
            rpkm[fam] = (cnt * 1_000_000) / total_reads
    return rpkm


# ---------------------------------------------------------------------------
# Gene-length estimation from hit subjects
# ---------------------------------------------------------------------------

def estimate_family_lengths(
    best_hits: dict[str, dict],
    protein_lengths: dict[str, int],
) -> dict[str, float]:
    """
    Compute the mean hit-subject protein length (aa) per family.

    Parameters
    ----------
    best_hits:
        Output of :func:`select_best_hits`.
    protein_lengths:
        {protein_id: length_aa} from :func:`~osmotool.functions.parse_protein_lengths`,
        keyed by the *query* protein ID (the assembly's own Prodigal-called
        ORF, not the reference database subject it aligned to) -- RPKM
        normalises by the length of the gene actually found in the genome
        being annotated, not by the length of whichever reference protein
        it happened to match.

    Returns
    -------
    dict[family_label, mean_length_aa]
    """
    family_lens: dict[str, list[int]] = defaultdict(list)
    for row in best_hits.values():
        query_id = row["qseqid"]
        family = gene_family_from_header(row["sseqid"])
        if query_id in protein_lengths:
            family_lens[family].append(protein_lengths[query_id])
    return {
        fam: sum(lens) / len(lens)
        for fam, lens in family_lens.items()
        if lens
    }
