"""
runners.py — high-level orchestrators for osmotool subcommands

profile  :  FASTQ reads → diamond blastx → RPM counts
annotate :  Assembly FASTA (or proteins) → Prodigal → diamond blastp → RPKM counts
"""

from __future__ import annotations

import logging
from pathlib import Path

from osmotool.functions import (
    run_diamond_blastx,
    run_diamond_blastp,
    run_hmmscan,
    run_orfm,
    run_prodigal,
    parse_diamond_output,
    parse_hmmscan_output,
    parse_protein_lengths,
    load_diamond_cutoffs,
    load_cascade_config,
    load_excluded_families,
    extract_fastq_subset,
    strip_orf_suffix,
    count_reads,
)
from osmotool.quantifier import (
    select_best_hits,
    count_hits,
    compute_rpm,
    compute_rpkm,
    estimate_family_lengths,
    filter_by_family_cutoff,
    apply_cascade_check,
)
from osmotool.output import write_gene_counts, write_aln_stats, build_aln_stats
from osmotool.utils import resolve_out_prefix, check_tool

log = logging.getLogger("osmotool")


# ---------------------------------------------------------------------------
# profile subcommand
# ---------------------------------------------------------------------------

def run_profile(
    db: str,
    out_prefix: str,
    *,
    reads1: str | None = None,
    reads2: str | None = None,
    singles: str | None = None,
    min_identity: float = 0.80,
    min_query_cover: float = 0.80,
    min_subject_cover: float = 0.0,
    min_seqlen: int = 20,
    evalue: float = 1e-5,
    threads: int = 4,
    keep_aln: bool = False,
    tmpdir: str | None = None,
    total_reads: int | None = None,
    hmm_db: str | None = None,
    cascade_config: str | None = None,
    exclude_families: str | None = None,
) -> tuple[Path, Path]:
    """
    Profile FASTQ reads against the osmodiamond database.

    Parameters
    ----------
    tmpdir:
        Directory for temporary files (merged paired-end FASTQ).
        Defaults to ``$TMPDIR`` (node-local scratch on HPC).
    total_reads:
        Pre-computed total read count.  Supply to skip on-the-fly counting
        (useful when already known from FastQC / MultiQC).
    hmm_db, cascade_config:
        Optional two-stage DIAMOND+HMM cascade. DIAMOND stays the primary
        caller for every read (needed for speed at millions-of-reads
        scale, and HMM's GA cutoffs don't apply to read fragments
        anyway), but reads whose DIAMOND call lands on a family flagged
        in osmo_refdb's <release>.profile_cascade.tsv (real, demonstrated
        precision problems, e.g. betL vs betT/CaiT) get pulled out,
        6-frame translated, and re-scored with hmmscan using a
        short-read-calibrated threshold (NOT --cut_ga -- see
        run_hmmscan's docstring) as a second opinion; failing reads are
        dropped. Only that flagged subset pays the extra cost, not the
        whole dataset. Both must be supplied together, and require orfm +
        hmmscan on PATH.
    exclude_families:
        Optional path to osmo_refdb's <release>.profile_excluded_families.txt
        (families marked scope: annotate_only, e.g. murB -- a
        near-universal housekeeping gene with no osmoadaptation-specific
        signal at the read level). Those families are still searched
        normally; this only drops them from the reported gene_counts.tsv.
    """
    check_tool("diamond")
    if cascade_config is not None:
        if hmm_db is None:
            raise ValueError("hmm_db is required when cascade_config is given")
        check_tool("orfm")
        check_tool("hmmscan")

    paired = reads1 is not None and reads2 is not None
    resolve_out_prefix(out_prefix)

    # --- count input reads (skip if pre-supplied) ---
    if total_reads is None:
        log.info("Counting input reads…")
        total_reads = count_reads(reads1, reads2, singles)
    log.info("  total reads: %d", total_reads)

    # --- run DIAMOND blastx ---
    aln_file = Path(str(out_prefix) + ".blastx.tsv")
    run_diamond_blastx(
        db,
        aln_file,
        reads1=reads1,
        reads2=reads2,
        singles=singles,
        min_identity=min_identity,
        min_query_cover=min_query_cover,
        min_seqlen=min_seqlen,
        evalue=evalue,
        threads=threads,
        tmpdir=tmpdir,
    )

    # --- parse alignments & best-hit selection ---
    alignments = list(parse_diamond_output(aln_file, min_subject_cover=min_subject_cover))
    filtered_reads = len({row["qseqid"] for row in alignments})
    best_hits = select_best_hits(iter(alignments))

    # --- DIAMOND+HMM cascade: second opinion for reads on flagged families ---
    if cascade_config is not None:
        best_hits = _apply_profile_cascade(
            best_hits, reads1, reads2, singles, hmm_db, cascade_config,
            out_prefix, threads, keep_aln,
        )

    # --- counting & RPM ---
    counts = count_hits(best_hits, paired=paired)
    if exclude_families is not None:
        excluded = load_excluded_families(exclude_families)
        counts = {fam: cnt for fam, cnt in counts.items() if fam not in excluded}
    rpm = compute_rpm(counts, total_reads)

    # --- write outputs ---
    counts_path = write_gene_counts(
        out_prefix, counts, rpm,
        total_reads=total_reads,
        filtered_reads=filtered_reads,
        normalisation="rpm",
    )
    stats = build_aln_stats(
        mode="profile",
        total_reads=total_reads,
        filtered_reads=filtered_reads,
        counts=counts,
        db=str(db),
        min_identity=min_identity,
        min_query_cover=min_query_cover,
        min_subject_cover=min_subject_cover,
        evalue=evalue,
    )
    stats_path = write_aln_stats(out_prefix, stats)

    if not keep_aln:
        aln_file.unlink(missing_ok=True)
    else:
        log.info("Alignment file retained: %s", aln_file)

    log.info("profile complete: %s  %s", counts_path, stats_path)
    return counts_path, stats_path


def _apply_profile_cascade(
    best_hits: dict[str, dict],
    reads1: str | None,
    reads2: str | None,
    singles: str | None,
    hmm_db: str,
    cascade_config: str,
    out_prefix: str,
    threads: int,
    keep_aln: bool,
) -> dict[str, dict]:
    """Second-opinion check for `profile`'s DIAMOND+HMM cascade: pull just
    the reads whose DIAMOND call landed on a flagged family, re-score them
    with hmmscan (raw bit scores), and drop any that fail that family's
    short-read-calibrated threshold. See run_profile's docstring."""
    cutoffs = load_cascade_config(cascade_config)
    flagged_read_ids = {
        qid for qid, row in best_hits.items()
        if row["sseqid"].split("|")[0] in cutoffs
    }
    if not flagged_read_ids:
        log.info("Cascade: no reads called to a flagged family, nothing to re-check.")
        return best_hits

    log.info("Cascade: re-checking %d/%d reads called to a flagged family (%s)",
              len(flagged_read_ids), len(best_hits), ", ".join(sorted(cutoffs)))

    subset_fastq = Path(str(out_prefix) + ".cascade_subset.fastq")
    input_files = [f for f in (reads1, reads2, singles) if f is not None]
    extract_fastq_subset(input_files, flagged_read_ids, subset_fastq)

    orfs_faa = Path(str(out_prefix) + ".cascade.orfs.faa")
    run_orfm(subset_fastq, orfs_faa)

    hmm_tblout = Path(str(out_prefix) + ".cascade.hmmscan.tblout")
    run_hmmscan(hmm_db, orfs_faa, hmm_tblout, threads=threads, use_cut_ga=False)

    cascade_scores: dict[str, dict[str, float]] = {}
    for hit in parse_hmmscan_output(hmm_tblout):
        read_id = strip_orf_suffix(hit["qseqid"])
        family = hit["sseqid"]
        per_read = cascade_scores.setdefault(read_id, {})
        if hit["bitscore"] > per_read.get(family, 0.0):
            per_read[family] = hit["bitscore"]

    filtered = apply_cascade_check(best_hits, cascade_scores, cutoffs)
    log.info("Cascade: %d/%d flagged-family reads confirmed by HMM, rest dropped",
              sum(1 for qid in flagged_read_ids if qid in filtered), len(flagged_read_ids))

    if not keep_aln:
        subset_fastq.unlink(missing_ok=True)
        orfs_faa.unlink(missing_ok=True)
        hmm_tblout.unlink(missing_ok=True)
    else:
        log.info("Cascade intermediate files retained: %s %s %s",
                  subset_fastq, orfs_faa, hmm_tblout)

    return filtered


# ---------------------------------------------------------------------------
# annotate subcommand
# ---------------------------------------------------------------------------

def run_annotate(
    db: str,
    assembly: str,
    out_prefix: str,
    *,
    proteins: str | None = None,
    method: str = "hmm",
    hmm_db: str | None = None,
    diamond_cutoffs: str | None = None,
    min_identity: float = 0.80,
    min_query_cover: float = 0.80,
    min_subject_cover: float = 0.0,
    evalue: float = 1e-5,
    threads: int = 4,
    keep_aln: bool = False,
    keep_proteins: bool = False,
    tmpdir: str | None = None,
    exclude_families: str | None = None,
) -> dict[str, tuple[Path, Path]]:
    """
    Annotate an assembly (or pre-called proteins) against osmodiamond.

    If *proteins* is None, Prodigal is run first in metagenomic mode.
    *tmpdir* is currently unused in annotate mode but accepted for a
    consistent interface (e.g. future chunked Prodigal runs).

    method: "hmm" (default), "diamond", or "both".

    Unlike `profile` (short reads), annotate works on full-length
    Prodigal-called ORFs -- exactly the scale osmo_refdb's HMM GA cutoffs
    are calibrated against, so hmmscan --cut_ga is a meaningful, cheap
    check here (a bacterial genome is a few thousand proteins, not
    millions of reads) and gives a per-family specificity gate DIAMOND's
    flat --min_identity doesn't have. Real-genome testing (E. coli K-12)
    confirmed this in practice: HMM correctly rejected a betT/CaiT
    paralog DIAMOND called as betL, and correctly resolved proP to its
    true single copy where DIAMOND over-called 4. "hmm"/"both" require
    --hmm_db and hmmscan on PATH; pass method="diamond" to avoid the
    HMMER dependency, at the cost of reintroducing those failure modes.

    diamond_cutoffs: optional path to osmo_refdb's
    <release>.diamond_cutoffs.tsv, giving DIAMOND calls the same kind of
    per-family specificity gate HMM's GA cutoff already has. Deliberately
    NOT offered in `profile`: those cutoffs are calibrated against
    full-length protein-vs-protein DIAMOND alignments (same scale as
    Prodigal ORFs here), and applying them to short read fragments would
    hit the identical scale-mismatch problem HMM's GA cutoffs have on
    short reads (see run_hmmscan's docstring) -- a short read can never
    reach a bitscore calibrated against full-length sequences.

    exclude_families:
        Optional path to osmo_refdb's <release>.annotate_excluded_families.txt
        -- currently just decoy references (families.yaml:
        decoy_from_negatives, e.g. betL_decoy), which must never appear as
        a reported gene family in either mode. Those sequences are still
        searched normally (that's the whole point -- they need to be able
        to win DIAMOND's best-hit contest away from a mislabeled call);
        this only drops them from the reported gene_counts.tsv. Unlike
        `profile`'s exclude_families, annotate_only families (e.g. murB)
        are deliberately NOT in this list -- they should stay visible here.

    Returns
    -------
    {method_name: (gene_counts_path, aln_stats_path)}
    """
    if method not in ("diamond", "hmm", "both"):
        raise ValueError(f"Unknown method '{method}': expected diamond, hmm, or both")
    if method in ("hmm", "both") and hmm_db is None:
        raise ValueError("hmm_db is required when method is 'hmm' or 'both'")

    if method in ("diamond", "both"):
        check_tool("diamond")
    if method in ("hmm", "both"):
        check_tool("hmmscan")
    if proteins is None:
        check_tool("prodigal")

    resolve_out_prefix(out_prefix)

    # --- optional Prodigal step ---
    if proteins is None:
        log.info("Running Prodigal on assembly: %s", assembly)
        proteins_path = Path(str(out_prefix) + ".prodigal.faa")
        run_prodigal(assembly, proteins_path)
    else:
        proteins_path = Path(proteins)

    # --- protein lengths for RPKM ---
    log.info("Parsing protein lengths from %s", proteins_path)
    protein_lengths = parse_protein_lengths(proteins_path)
    total_proteins = len(protein_lengths)
    log.info("  %d proteins", total_proteins)

    results: dict[str, tuple[Path, Path]] = {}

    if method in ("diamond", "both"):
        aln_file = Path(str(out_prefix) + ".blastp.tsv")
        run_diamond_blastp(
            db,
            proteins_path,
            aln_file,
            min_identity=min_identity,
            min_query_cover=min_query_cover,
            evalue=evalue,
            threads=threads,
        )
        alignments = list(parse_diamond_output(aln_file, min_subject_cover=min_subject_cover))
        if diamond_cutoffs is not None:
            cutoffs = load_diamond_cutoffs(diamond_cutoffs)
            n_before = len(alignments)
            alignments = list(filter_by_family_cutoff(alignments, cutoffs))
            log.info("Applied per-family DIAMOND cutoffs (%s): %d -> %d alignment rows",
                     diamond_cutoffs, n_before, len(alignments))
        results["diamond"] = _finish_annotate_method(
            "diamond", out_prefix if method == "diamond" else f"{out_prefix}.diamond",
            alignments, protein_lengths, total_proteins, db=str(db),
            min_identity=min_identity, min_query_cover=min_query_cover,
            min_subject_cover=min_subject_cover, evalue=evalue,
            exclude_families=exclude_families,
        )
        if not keep_aln:
            aln_file.unlink(missing_ok=True)
        else:
            log.info("Alignment file retained: %s", aln_file)

    if method in ("hmm", "both"):
        hmm_aln_file = Path(str(out_prefix) + ".hmmscan.tblout")
        run_hmmscan(hmm_db, proteins_path, hmm_aln_file, threads=threads)
        hmm_alignments = list(parse_hmmscan_output(hmm_aln_file))
        results["hmm"] = _finish_annotate_method(
            "hmm", out_prefix if method == "hmm" else f"{out_prefix}.hmm",
            hmm_alignments, protein_lengths, total_proteins, db=str(hmm_db),
            min_identity=None, min_query_cover=None, min_subject_cover=None, evalue=None,
            exclude_families=exclude_families,
        )
        if not keep_aln:
            hmm_aln_file.unlink(missing_ok=True)
        else:
            log.info("Alignment file retained: %s", hmm_aln_file)

    if not keep_proteins and proteins is None:
        # remove temp Prodigal file we created
        proteins_path.unlink(missing_ok=True)
    else:
        log.info("Protein FASTA retained: %s", proteins_path)

    for method_name, (counts_path, stats_path) in results.items():
        log.info("annotate (%s) complete: %s  %s", method_name, counts_path, stats_path)
    return results


def _finish_annotate_method(
    method_name: str,
    out_prefix: str | Path,
    alignments: list[dict],
    protein_lengths: dict[str, int],
    total_proteins: int,
    *,
    db: str,
    min_identity: float | None,
    min_query_cover: float | None,
    min_subject_cover: float | None,
    evalue: float | None,
    exclude_families: str | None = None,
) -> tuple[Path, Path]:
    """Shared best-hit selection / counting / RPKM / output-writing tail for
    one annotate method (diamond or hmm), factored out since both methods
    produce identically-shaped alignment rows (see parse_hmmscan_output's
    docstring) and go through the same quantification and output steps."""
    filtered_reads = len({row["qseqid"] for row in alignments})
    best_hits = select_best_hits(iter(alignments))

    # annotate mode: single-end proteins (no paired weighting)
    counts = count_hits(best_hits, paired=False)
    if exclude_families is not None:
        excluded = load_excluded_families(exclude_families)
        counts = {fam: cnt for fam, cnt in counts.items() if fam not in excluded}

    # RPKM using mean hit-subject protein lengths
    family_lengths_aa = estimate_family_lengths(best_hits, protein_lengths)
    rpkm = compute_rpkm(counts, total_proteins, family_lengths_aa)

    counts_path = write_gene_counts(
        out_prefix, counts, rpkm,
        total_reads=total_proteins,
        filtered_reads=filtered_reads,
        normalisation="rpkm",
    )
    stats = build_aln_stats(
        mode=f"annotate:{method_name}",
        total_reads=total_proteins,
        filtered_reads=filtered_reads,
        counts=counts,
        db=db,
        min_identity=min_identity if min_identity is not None else "n/a",
        min_query_cover=min_query_cover if min_query_cover is not None else "n/a",
        min_subject_cover=min_subject_cover if min_subject_cover is not None else "n/a",
        evalue=evalue if evalue is not None else "n/a",
    )
    stats_path = write_aln_stats(out_prefix, stats)

    return counts_path, stats_path
