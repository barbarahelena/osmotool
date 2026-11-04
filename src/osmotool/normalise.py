"""
normalise.py — copies-per-genome normalisation for ``profile`` mode

RPM (quantifier.compute_rpm) depends on everything else in the sample: host
or eukaryotic DNA, viral load, and the community's average genome size all
shift it. Dividing a gene family's length-normalised read count by that of
universal single-copy marker genes instead gives the average number of
copies of the family per genome in the community (equivalently, for a
family that is at most single-copy, the fraction of cells carrying it).

Method — the approach of MicrobeCensus (Nayfach & Pollard 2015, Genome
Biol 16:51), reimplemented here for DIAMOND rather than copied (MicrobeCensus
is GPL-3.0, and its trained cutoffs/coefficients are specific to RAPsearch2
v2.15 scores, so they would not transfer anyway):

1. Reads are searched against a separate marker DIAMOND db (osmo_refdb
   release's ``markers/``: the ten fetchMG/mOTUs universal single-copy COGs
   -- see osmo_refdb's markers.yaml) in a second pass, so markers never
   compete with osmo families for best hits.
2. For each family f (osmo or marker), reads per kilobase of *effective*
   gene length::

       RPK_f = count_f / (L_eff_f / 1000)
       L_eff = 3 * L_aa - (2c - 1) * R + 1

   L_aa is the hit-weighted mean length of the reference proteins f's reads
   aligned to, R the sample's median read length, c the --min_query_cover
   threshold. L_eff counts the read start positions that leave at least
   c * R of the read inside the gene (holds for any L >= c * R, including
   genes shorter than the read). Without it, short families are
   underestimated relative to long markers: cspA (~70 aa, L_eff 121 of
   210 bp) vs leuS (~860 aa, L_eff 2491 of 2580 bp) by ~1.7x at 150 bp
   reads and the default c = 0.8.
3. Genome coverage G: drop markers with zero reads (MicrobeCensus does the
   same), drop markers whose RPK is >= 1 scaled MAD (1.48 * MAD) from the
   median -- MicrobeCensus's outlier rule, which removes e.g. a marker
   duplicated in a dominant lineage -- and take the mean of the rest.
   Unweighted: MicrobeCensus's per-marker weights were trained against
   RAPsearch2 and would need recalibrating for DIAMOND (phase 2).
4. ``copies_per_genome_f = RPK_f / G``.

What this does NOT correct for (phase 2, via osmo_refdb's simulated-read
benchmark): per-family detection sensitivity. Both the osmo family and the
markers are counted with the same DIAMOND thresholds, so sensitivity
differences partly cancel, but not fully -- treat copies_per_genome as
comparable across samples, not as an exact absolute copy number.
"""

from __future__ import annotations

import gzip
import logging
import statistics
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from osmotool.utils import gene_family_from_header

log = logging.getLogger("osmotool")

#: MicrobeCensus's MAD scale constant (~1.4826 for normal data).
MAD_CONST = 1.48
#: Below this many markers with >=1 read, no genome-coverage estimate is made.
MIN_MARKERS = 5
#: Below this many marker fragments in total, the estimate is flagged as noisy.
LOW_MARKER_FRAGMENTS = 100


# ---------------------------------------------------------------------------
# Inputs: read length, reference lengths, marker manifest
# ---------------------------------------------------------------------------

def median_read_length(paths: list[str | Path], n_reads: int = 10_000) -> int:
    """Median read length (bp) of the first *n_reads* reads of the first
    readable FASTQ in *paths* (handles .gz). Same sampling idea as
    MicrobeCensus's auto_detect_read_length."""
    for path in paths:
        if path is None:
            continue
        opener = gzip.open if str(path).endswith(".gz") else open
        lengths: list[int] = []
        with opener(path, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    lengths.append(len(line.rstrip("\r\n")))
                    if len(lengths) >= n_reads:
                        break
        if lengths:
            return int(statistics.median(lengths))
    raise ValueError("no reads found to estimate read length from")


def reference_lengths(
    best_hits: dict[str, dict],
    ref_lengths_aa: dict[str, int],
) -> dict[str, float]:
    """Hit-weighted mean reference (subject) protein length, in aa, per
    family: each read contributes the length of the reference it aligned
    to. Unlike annotate mode's quantifier.estimate_family_lengths (keyed by
    the *query* ORF, since annotate has full-length query proteins), profile
    reads are fragments, so the reference is the only length available."""
    per_family: dict[str, list[int]] = defaultdict(list)
    for row in best_hits.values():
        length = ref_lengths_aa.get(row["sseqid"])
        if length:
            per_family[gene_family_from_header(row["sseqid"])].append(length)
    return {fam: sum(ls) / len(ls) for fam, ls in per_family.items()}


def load_marker_names(manifest: str | Path) -> list[str]:
    """Marker family names from osmo_refdb's markers/osmo_markers.manifest.tsv."""
    with open(manifest) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if not header or header[0] != "marker":
            raise ValueError(f"{manifest}: expected a 'marker' first column")
        return [line.split("\t", 1)[0] for line in fh if line.strip()]


# ---------------------------------------------------------------------------
# Effective length / RPK
# ---------------------------------------------------------------------------

def effective_length_bp(length_aa: float, read_length: int, min_query_cover: float) -> float:
    """Read start positions that leave >= min_query_cover of the read inside
    a gene of *length_aa*: ``3 * L_aa - (2c - 1) * R + 1`` (see module
    docstring). 0.0 when the gene is too short for any read to qualify."""
    gene_bp = 3.0 * length_aa
    if gene_bp < min_query_cover * read_length:
        return 0.0
    return max(gene_bp - (2.0 * min_query_cover - 1.0) * read_length + 1.0, 0.0)


def compute_rpk(
    counts: dict[str, float],
    lengths_aa: dict[str, float],
    read_length: int,
    min_query_cover: float,
) -> dict[str, float | None]:
    """Reads per kb of effective length. A family with zero reads is 0.0
    (no length needed); a family WITH reads but no usable length is None."""
    rpk: dict[str, float | None] = {}
    for fam, cnt in counts.items():
        if cnt == 0:
            rpk[fam] = 0.0
            continue
        eff = effective_length_bp(lengths_aa.get(fam, 0.0), read_length, min_query_cover)
        rpk[fam] = cnt / (eff / 1000.0) if eff > 0 else None
    return rpk


# ---------------------------------------------------------------------------
# Genome coverage from markers
# ---------------------------------------------------------------------------

@dataclass
class MarkerEstimate:
    """Result of :func:`estimate_genome_coverage`."""
    coverage_rpk: float | None          # G; None if not estimable
    used: list[str] = field(default_factory=list)
    outliers: list[str] = field(default_factory=list)
    undetected: list[str] = field(default_factory=list)
    total_fragments: float = 0.0
    note: str = ""


def estimate_genome_coverage(
    marker_rpk: dict[str, float | None],
    marker_counts: dict[str, float],
    *,
    min_markers: int = MIN_MARKERS,
    mad_const: float = MAD_CONST,
) -> MarkerEstimate:
    """Mean marker RPK after removing undetected markers and MAD outliers.

    Outlier rule (MicrobeCensus): drop |x - median| >= mad_const * MAD. When
    MAD is 0 (half or more of the markers share the median value exactly)
    that rule would drop everything, so keep only values equal to the median.
    """
    total = sum(marker_counts.values())
    detected = {m: v for m, v in marker_rpk.items() if v}   # drops 0.0 and None
    undetected = sorted(set(marker_rpk) - set(detected))
    est = MarkerEstimate(coverage_rpk=None, undetected=undetected, total_fragments=total)

    if len(detected) < min_markers:
        est.note = (f"only {len(detected)}/{len(marker_rpk)} markers detected "
                    f"(need >= {min_markers}); copies_per_genome not estimable")
        return est

    values = list(detected.values())
    med = statistics.median(values)
    mad = mad_const * statistics.median(abs(v - med) for v in values)
    for m, v in sorted(detected.items()):
        if (abs(v - med) < mad) if mad > 0 else (v == med):
            est.used.append(m)
        else:
            est.outliers.append(m)
    est.coverage_rpk = statistics.mean(detected[m] for m in est.used)
    if total < LOW_MARKER_FRAGMENTS:
        est.note = (f"only {total:.1f} marker fragments; copies_per_genome is noisy "
                    f"at this depth")
    return est


def copies_per_genome(
    rpk: dict[str, float | None],
    coverage_rpk: float | None,
) -> dict[str, float | None]:
    """RPK_f / G for every family; None where either side is unknown."""
    if not coverage_rpk:
        return {fam: None for fam in rpk}
    return {fam: (v / coverage_rpk if v is not None else None) for fam, v in rpk.items()}
