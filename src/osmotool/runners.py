"""
runners.py — high-level orchestrators for osmotool subcommands

profile  :  FASTQ reads → diamond blastx → RPM counts
annotate :  Assembly FASTA (or proteins) → Prodigal → diamond blastp → copies_per_kb counts
"""

from __future__ import annotations

import logging
from pathlib import Path

from osmotool.functions import (
    AssemblyStats,
    compute_assembly_stats,
    run_diamond_blastx,
    run_diamond_blastp,
    run_hmmsearch,
    run_orfm,
    run_prodigal,
    parse_diamond_output,
    parse_hmmsearch_output,
    parse_protein_lengths,
    parse_orf_coordinates,
    load_diamond_cutoffs,
    load_cascade_config,
    load_excluded_families,
    extract_fastq_subset,
    strip_orf_suffix,
    strip_mate_tag,
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
    extract_gene_coordinates,
)
from osmotool.normalise import (
    median_read_length,
    reference_lengths,
    load_marker_names,
    compute_rpk,
    estimate_genome_coverage,
    copies_per_genome,
)
from osmotool.output import (
    write_gene_counts, write_aln_stats, build_aln_stats, write_system_calls,
    write_gene_coordinates, write_marker_report, PROFILE_SYSTEMS_NOTE,
)
from osmotool.systems import detect_systems, DEFAULT_SYSTEMS, load_systems_config
from osmotool.utils import resolve_out_prefix, check_tool, MarkerDb

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
    markers: MarkerDb | None = None,
    train_refs: str | None = None,
) -> tuple[Path, Path]:
    """
    Profile FASTQ reads against the osmodiamond database.

    Parameters
    ----------
    markers, train_refs:
        Copies-per-genome normalisation (see normalise.py). *markers* is the
        release's markers/ subdirectory (universal single-copy genes, searched
        in a second, separate DIAMOND pass with the same thresholds);
        *train_refs* is osmo_refdb.train_refs.faa, for osmo reference
        lengths. Both required together; adds a copies_per_genome column to
        gene_counts.tsv and writes <out_prefix>.markers.tsv. None (the
        library default) skips it; the CLI enables it whenever DATABASE
        ships markers/, unless --no_markers.
    tmpdir:
        Directory for temporary files (merged paired-end FASTQ).
        Defaults to ``$TMPDIR`` (node-local scratch on HPC).
    total_reads:
        Pre-computed total read count.  Supply to skip on-the-fly counting
        (useful when already known from FastQC / MultiQC).
    hmm_db, cascade_config:
        Two-stage DIAMOND+HMM cascade. DIAMOND stays the primary caller for
        every read (needed for speed at millions-of-reads scale, and HMM's
        GA cutoffs don't apply to read fragments anyway), but reads whose
        DIAMOND call lands on a family flagged in osmo_refdb's
        <release>.profile_cascade.tsv (real, demonstrated precision
        problems -- e.g. betL vs betT/CaiT, and several families as low as
        5-25% DIAMOND-only precision in osmo_refdb's own benchmark) get
        pulled out, 6-frame translated, and re-scored with hmmsearch using a
        short-read-calibrated threshold (NOT --cut_ga -- see
        run_hmmsearch's docstring) as a second opinion; failing reads are
        dropped. Only that flagged subset pays the extra cost, not the
        whole dataset. Both must be supplied together, and require orfm +
        hmmsearch on PATH. At this function level these default to None
        (DIAMOND-only) for library callers, but osmotool's CLI (`profile`)
        enables the cascade by default -- given the precision numbers
        above, DIAMOND-only should be treated as the opt-out, not the
        default, for any caller reporting results.
    exclude_families:
        Optional path to osmo_refdb's <release>.profile_excluded_families.txt
        (families marked scope: annotate_only, e.g. murB -- a
        near-universal housekeeping gene with no osmoadaptation-specific
        signal at the read level). Those families are still searched
        normally; this only drops them from the reported gene_counts.tsv.

    No multi-subunit system detection here (see systems.py / `annotate`'s
    --skip_systems): short reads don't carry contig-level coordinates, so
    the proximity/colocalization check that makes a system-level "present"
    call meaningful (are these subunit genes actually next to each other,
    or just independently present somewhere in the community) simply isn't
    checkable at read length. This mode reports each gene family's
    abundance independently, same as before Task 2 added the systems
    layer -- a system-level claim (opuA/opuB/opuC/mrp/ect/ots/trk/ktr)
    requires `annotate` mode (or a dedicated assembly+binning step) on the
    same samples. See this project's README for the same note.
    """
    check_tool("diamond")
    if markers is not None and train_refs is None:
        raise ValueError("train_refs is required when markers is given")
    if cascade_config is not None:
        if hmm_db is None:
            raise ValueError("hmm_db is required when cascade_config is given")
        check_tool("orfm")
        check_tool("hmmsearch")

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

    # --- copies per genome, via universal single-copy markers ---
    copies, marker_stats = None, None
    if markers is not None:
        copies, marker_stats = _normalise_to_markers(
            markers, train_refs, best_hits, counts,
            reads1=reads1, reads2=reads2, singles=singles, paired=paired,
            out_prefix=out_prefix, keep_aln=keep_aln,
            min_identity=min_identity, min_query_cover=min_query_cover,
            min_subject_cover=min_subject_cover, min_seqlen=min_seqlen,
            evalue=evalue, threads=threads, tmpdir=tmpdir,
        )

    # --- write outputs ---
    counts_path = write_gene_counts(
        out_prefix, counts, rpm,
        total_reads=total_reads,
        filtered_reads=filtered_reads,
        normalisation="rpm",
        copies_per_genome=copies,
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
        note_systems=PROFILE_SYSTEMS_NOTE,
        extra=marker_stats,
    )
    stats_path = write_aln_stats(out_prefix, stats)

    if not keep_aln:
        aln_file.unlink(missing_ok=True)
    else:
        log.info("Alignment file retained: %s", aln_file)

    log.info("profile complete: %s  %s", counts_path, stats_path)
    return counts_path, stats_path


def _normalise_to_markers(
    markers: MarkerDb,
    train_refs: str,
    best_hits: dict[str, dict],
    counts: dict[str, float],
    *,
    reads1: str | None,
    reads2: str | None,
    singles: str | None,
    paired: bool,
    out_prefix: str,
    keep_aln: bool,
    min_identity: float,
    min_query_cover: float,
    min_subject_cover: float,
    min_seqlen: int,
    evalue: float,
    threads: int,
    tmpdir: str | None,
) -> tuple[dict[str, float | None], dict]:
    """Second DIAMOND pass against the marker db and the copies-per-genome
    computation (see normalise.py). Markers get exactly the same alignment
    thresholds and best-hit/paired counting as the osmo families, so
    detection-sensitivity differences between them are as small as this
    phase-1 method can make them. Returns (copies_per_genome, aln_stats rows)."""
    read_length = median_read_length([reads1, singles])
    log.info("Marker normalisation: median read length %d bp", read_length)

    marker_aln = Path(str(out_prefix) + ".markers.blastx.tsv")
    run_diamond_blastx(
        markers.dmnd, marker_aln,
        reads1=reads1, reads2=reads2, singles=singles,
        min_identity=min_identity, min_query_cover=min_query_cover,
        min_seqlen=min_seqlen, evalue=evalue, threads=threads, tmpdir=tmpdir,
    )
    marker_hits = select_best_hits(parse_diamond_output(marker_aln, min_subject_cover=min_subject_cover))

    names = load_marker_names(markers.manifest)
    decoys = load_excluded_families(markers.exclude_families) if markers.exclude_families else set()
    raw = count_hits(marker_hits, paired=paired, known_families=names)
    # Only real markers count: decoy hits (and any unexpected label) are
    # dropped here, after best-hit selection has already let them win reads
    # away from the markers they mimic.
    unexpected = sorted(set(raw) - set(names) - decoys)
    if unexpected:
        log.warning("Marker db has labels not in its manifest, ignored: %s", ", ".join(unexpected))
    marker_counts = {m: raw.get(m, 0.0) for m in names}

    marker_len = reference_lengths(marker_hits, parse_protein_lengths(markers.refs))
    marker_rpk = compute_rpk(marker_counts, marker_len, read_length, min_query_cover)
    estimate = estimate_genome_coverage(marker_rpk, marker_counts)
    if estimate.note:
        log.warning("Marker normalisation: %s", estimate.note)
    if estimate.outliers:
        log.info("Marker normalisation: outlier markers excluded: %s", ", ".join(estimate.outliers))

    osmo_len = reference_lengths(best_hits, parse_protein_lengths(train_refs))
    osmo_rpk = compute_rpk(counts, osmo_len, read_length, min_query_cover)
    copies = copies_per_genome(osmo_rpk, estimate.coverage_rpk)

    write_marker_report(out_prefix, marker_counts, marker_len, marker_rpk, estimate,
                        read_length=read_length)

    if not keep_aln:
        marker_aln.unlink(missing_ok=True)
    else:
        log.info("Marker alignment file retained: %s", marker_aln)

    cov = estimate.coverage_rpk
    stats = {
        "marker_db":            str(markers.dmnd),
        "read_length":          read_length,
        "marker_fragments":     f"{estimate.total_fragments:.1f}",
        "genome_coverage_rpk":  "NA" if cov is None else f"{cov:.6g}",
        "markers_used":         ",".join(estimate.used) or "none",
        "markers_outliers":     ",".join(estimate.outliers) or "none",
        "markers_undetected":   ",".join(estimate.undetected) or "none",
    }
    if estimate.note:
        stats["note_markers"] = estimate.note
    return copies, stats


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
    with hmmsearch (raw bit scores), and drop any that fail that family's
    short-read-calibrated threshold. See run_profile's docstring.

    best_hits' qseqids are _R1/_R2-tagged for paired input (see
    functions._merge_paired_with_mate_tags), but the actual R1/R2 FASTQ
    files on disk are untagged -- extract_fastq_subset() reads reads back
    out of *those* original files, so flagged_read_ids (and later,
    cascade_scores) must use the untagged id via strip_mate_tag() to find
    them. Single-end qseqids are never tagged, so strip_mate_tag() is a
    no-op there.
    """
    cutoffs = load_cascade_config(cascade_config)
    flagged_qids = {
        qid for qid, row in best_hits.items()
        if row["sseqid"].split("|")[0] in cutoffs
    }
    if not flagged_qids:
        log.info("Cascade: no reads called to a flagged family, nothing to re-check.")
        return best_hits

    log.info("Cascade: re-checking %d/%d reads called to a flagged family (%s)",
              len(flagged_qids), len(best_hits), ", ".join(sorted(cutoffs)))

    subset_fastq = Path(str(out_prefix) + ".cascade_subset.fastq")
    input_files = [f for f in (reads1, reads2, singles) if f is not None]
    # extract_fastq_subset reads back from the original, untagged R1/R2
    # files, so it needs untagged ids -- see this function's docstring.
    extract_fastq_subset(input_files, {strip_mate_tag(qid) for qid in flagged_qids}, subset_fastq)

    orfs_faa = Path(str(out_prefix) + ".cascade.orfs.faa")
    run_orfm(subset_fastq, orfs_faa)

    hmm_tblout = Path(str(out_prefix) + ".cascade.hmmsearch.tblout")
    run_hmmsearch(hmm_db, orfs_faa, hmm_tblout, threads=threads, use_cut_ga=False)

    cascade_scores: dict[str, dict[str, float]] = {}
    for hit in parse_hmmsearch_output(hmm_tblout):
        read_id = strip_orf_suffix(hit["qseqid"])
        family = hit["sseqid"]
        per_read = cascade_scores.setdefault(read_id, {})
        if hit["bitscore"] > per_read.get(family, 0.0):
            per_read[family] = hit["bitscore"]

    # cascade_scores is keyed by the untagged original read ID (recovered
    # via strip_orf_suffix from headers orfm/hmmsearch saw, which trace back
    # through extract_fastq_subset's verbatim copy of the original,
    # untagged FASTQ headers) -- re-key it to match best_hits' _R1/_R2-
    # tagged qseqids so apply_cascade_check's per-qid lookup lines up.
    # NOTE: if both mates of one fragment share an identical original
    # header (exactly the ambiguity the tagging exists to resolve) AND
    # both independently land on a flagged family, their re-check scores
    # get merged under that one shared original ID in cascade_scores
    # (orfm/hmmsearch have no mate awareness), so each mate's cascade
    # decision here uses the best-of-both-mates score rather than a
    # strictly mate-specific one. A narrow, compounding edge case
    # (--cascade + paired + a real header collision + both mates
    # independently flagged) -- fully closing it would require
    # extract_fastq_subset to also re-tag headers per source file.
    scores_by_qid = {qid: cascade_scores.get(strip_mate_tag(qid), {}) for qid in best_hits}
    filtered = apply_cascade_check(best_hits, scores_by_qid, cutoffs)
    log.info("Cascade: %d/%d flagged-family reads confirmed by HMM, rest dropped",
              sum(1 for qid in flagged_qids if qid in filtered), len(flagged_qids))

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
    skip_systems: bool = False,
    systems_config: str | None = None,
) -> dict[str, tuple[Path, Path, Path | None, Path | None]]:
    """
    Annotate an assembly (or pre-called proteins) against osmodiamond.

    If *proteins* is None, Prodigal is run first in metagenomic mode.
    *tmpdir* is currently unused in annotate mode but accepted for a
    consistent interface (e.g. future chunked Prodigal runs).

    method: "hmm" (default), "diamond", or "both".

    Unlike `profile` (short reads), annotate works on full-length
    Prodigal-called ORFs -- exactly the scale osmo_refdb's HMM GA cutoffs
    are calibrated against, so hmmsearch --cut_ga is a meaningful, cheap
    check here (a bacterial genome is a few thousand proteins, not
    millions of reads) and gives a per-family specificity gate DIAMOND's
    flat --min_identity doesn't have. Real-genome testing (E. coli K-12)
    confirmed this in practice: HMM correctly rejected a betT/CaiT
    paralog DIAMOND called as betL, and correctly resolved proP to its
    true single copy where DIAMOND over-called 4. "hmm"/"both" require
    --hmm_db and hmmsearch on PATH; pass method="diamond" to avoid the
    HMMER dependency, at the cost of reintroducing those failure modes.

    diamond_cutoffs: optional path to osmo_refdb's
    <release>.diamond_cutoffs.tsv, giving DIAMOND calls the same kind of
    per-family specificity gate HMM's GA cutoff already has. Deliberately
    NOT offered in `profile`: those cutoffs are calibrated against
    full-length protein-vs-protein DIAMOND alignments (same scale as
    Prodigal ORFs here), and applying them to short read fragments would
    hit the identical scale-mismatch problem HMM's GA cutoffs have on
    short reads (see run_hmmsearch's docstring) -- a short read can never
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

    skip_systems:
        Skip the Task 2 complex/system-aware detection layer (systems.py)
        entirely -- no <out_prefix>.systems.tsv, no Prodigal-header
        coordinate parsing, and consequently no <out_prefix>.gene_coordinates.tsv
        either (both depend on the same parsed ORF coordinates). Off by
        default; set True to save the (small) extra parsing cost, or if you
        just don't need it.

    systems_config:
        Optional path to a TSV overriding max_distance_genes/max_distance_bp
        per system (columns: system, max_distance_genes, max_distance_bp) --
        see systems.load_systems_config. Mandatory/accessory/fusion
        membership is not overridable this way (structural, not a
        calibration knob); only distance defaults are.

    Returns
    -------
    {method_name: (gene_counts_path, aln_stats_path, systems_path_or_None, coordinates_path_or_None)}
    systems_path/coordinates_path are None when skip_systems=True, or when
    no ORF had coordinates (e.g. --proteins supplied with non-Prodigal
    headers).
    """
    if method not in ("diamond", "hmm", "both"):
        raise ValueError(f"Unknown method '{method}': expected diamond, hmm, or both")
    if method in ("hmm", "both") and hmm_db is None:
        raise ValueError("hmm_db is required when method is 'hmm' or 'both'")

    if method in ("diamond", "both"):
        check_tool("diamond")
    if method in ("hmm", "both"):
        check_tool("hmmsearch")
    if proteins is None:
        check_tool("prodigal")

    resolve_out_prefix(out_prefix)

    # --- whole-assembly stats for the gene_counts.tsv header ---
    # Read from the assembly even when --proteins skips Prodigal; ASSEMBLY
    # is a required argument either way.
    assembly_stats = compute_assembly_stats(assembly)
    log.info(
        "Assembly: %d bp, GC %s, %d contigs, N50 %d",
        assembly_stats.genome_bp,
        "NA" if assembly_stats.gc is None else f"{assembly_stats.gc:.5f}",
        assembly_stats.n_contigs,
        assembly_stats.contig_n50,
    )

    # --- optional Prodigal step ---
    if proteins is None:
        log.info("Running Prodigal on assembly: %s", assembly)
        proteins_path = Path(str(out_prefix) + ".prodigal.faa")
        run_prodigal(assembly, proteins_path)
    else:
        proteins_path = Path(proteins)

    # --- protein lengths for copies_per_kb normalisation ---
    log.info("Parsing protein lengths from %s", proteins_path)
    protein_lengths = parse_protein_lengths(proteins_path)
    total_proteins = len(protein_lengths)
    log.info("  %d proteins", total_proteins)

    # --- ORF coordinates for Task 2's system/complex detection layer ---
    orf_coords: dict = {}
    if not skip_systems:
        orf_coords = parse_orf_coordinates(proteins_path)
        if not orf_coords:
            log.warning(
                "No ORF coordinates parsed from %s (expected Prodigal-style "
                "headers) -- skipping system/complex detection for this run",
                proteins_path,
            )
        systems_def = (
            load_systems_config(systems_config) if systems_config is not None
            else DEFAULT_SYSTEMS
        )
    else:
        systems_def = DEFAULT_SYSTEMS

    results: dict[str, tuple[Path, Path, Path | None, Path | None]] = {}

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
            orf_coords=orf_coords, systems_def=systems_def,
            assembly_stats=assembly_stats,
        )
        if not keep_aln:
            aln_file.unlink(missing_ok=True)
        else:
            log.info("Alignment file retained: %s", aln_file)

    if method in ("hmm", "both"):
        hmm_aln_file = Path(str(out_prefix) + ".hmmsearch.tblout")
        run_hmmsearch(hmm_db, proteins_path, hmm_aln_file, threads=threads)
        hmm_alignments = list(parse_hmmsearch_output(hmm_aln_file))
        results["hmm"] = _finish_annotate_method(
            "hmm", out_prefix if method == "hmm" else f"{out_prefix}.hmm",
            hmm_alignments, protein_lengths, total_proteins, db=str(hmm_db),
            min_identity=None, min_query_cover=None, min_subject_cover=None, evalue=None,
            exclude_families=exclude_families,
            orf_coords=orf_coords, systems_def=systems_def,
            assembly_stats=assembly_stats,
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

    for method_name, (counts_path, stats_path, systems_path, coords_path) in results.items():
        log.info("annotate (%s) complete: %s  %s%s%s", method_name, counts_path, stats_path,
                  f"  {systems_path}" if systems_path else "",
                  f"  {coords_path}" if coords_path else "")
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
    orf_coords: dict | None = None,
    systems_def: dict | None = None,
    assembly_stats: AssemblyStats | None = None,
) -> tuple[Path, Path, Path | None, Path | None]:
    """Shared best-hit selection / counting / copies_per_kb / output-writing tail for
    one annotate method (diamond or hmm), factored out since both methods
    produce identically-shaped alignment rows (see parse_hmmsearch_output's
    docstring) and go through the same quantification and output steps.

    orf_coords/systems_def: Task 2's system/complex detection (systems.py),
    run against this method's own best_hits -- diamond and hmm each get
    their own independent systems.tsv/gene_coordinates.tsv, same as they
    already get independent gene_counts.tsv/aln_stats.tsv. Both skipped
    (systems_path/coords_path=None) when orf_coords is empty/None -- see
    run_annotate's docstring for why that happens and how it degrades
    gracefully."""
    # Proteins with at least one hit that passed filtering, counted before
    # best-hit selection and decoy exclusion. Written under its legacy
    # "filtered_reads" name (shared with profile) and as proteins_with_hits.
    filtered_reads = len({row["qseqid"] for row in alignments})
    best_hits = select_best_hits(iter(alignments))

    excluded: set[str] = (
        load_excluded_families(exclude_families) if exclude_families is not None else set()
    )

    # annotate mode: single-end proteins (no paired weighting)
    counts = count_hits(best_hits, paired=False)
    if excluded:
        counts = {fam: cnt for fam, cnt in counts.items() if fam not in excluded}

    # copies_per_kb: gene-family copy number normalised by mean hit-subject
    # protein length and total predicted proteins in this genome. Same
    # formula shape as profile's RPM/RPKM but there's no sequencing depth
    # here to normalise against -- total_proteins is a per-genome ORF
    # count, not a read count -- so this is a compositional statistic
    # ("how much of this genome, length-corrected, is this family"), not
    # an abundance/expression estimate. Don't compare it directly against
    # profile mode's RPM values.
    family_lengths_aa = estimate_family_lengths(best_hits, protein_lengths)
    copies_per_kb = compute_rpkm(counts, total_proteins, family_lengths_aa)

    counts_path = write_gene_counts(
        out_prefix, counts, copies_per_kb,
        total_reads=total_proteins,
        filtered_reads=filtered_reads,
        normalisation="copies_per_kb",
        total_reads_label="total_proteins",
        assembly_stats=assembly_stats,
        proteins_with_hits=filtered_reads,
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

    systems_path: Path | None = None
    coords_path: Path | None = None
    if orf_coords:
        calls = detect_systems(best_hits, orf_coords, systems_def or DEFAULT_SYSTEMS)
        systems_path = write_system_calls(out_prefix, calls)

        coord_rows = extract_gene_coordinates(best_hits, orf_coords)
        if excluded:
            coord_rows = [row for row in coord_rows if row["gene"] not in excluded]
        coords_path = write_gene_coordinates(out_prefix, coord_rows)

    return counts_path, stats_path, systems_path, coords_path
