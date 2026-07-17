"""
handle_args.py — argument parsing for osmotool
"""

from __future__ import annotations

import argparse


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="osmotool",
        description=(
            "Screen osmoadaptation genes in metagenomic datasets using DIAMOND. "
            "Gene families: ectA, ectB, ectC (ectoine), betL (betaine transport), "
            "kdpA (K⁺ uptake), nhaA (Na⁺/H⁺ antiporter)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--version", action="version", version="%(prog)s (see pyproject.toml)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Enable debug logging."
    )

    sub = parser.add_subparsers(dest="subcommand", metavar="subcommand")
    sub.required = True

    _add_profile_parser(sub)
    _add_annotate_parser(sub)

    return parser


# ---------------------------------------------------------------------------
# Shared filter arguments
# ---------------------------------------------------------------------------

def _add_filter_args(p: argparse.ArgumentParser) -> None:
    filt = p.add_argument_group("DIAMOND filters")
    filt.add_argument(
        "--min_identity", type=float, default=0.80, metavar="FRAC",
        help="Minimum fractional sequence identity (0–1).",
    )
    filt.add_argument(
        "--min_query_cover", type=float, default=0.80, metavar="FRAC",
        help="Minimum query coverage (0–1).",
    )
    filt.add_argument(
        "--min_subject_cover", type=float, default=0.0, metavar="FRAC",
        help="Minimum subject/target coverage (0-1, default 0 = off). "
             "Filters out hits that only cover a small fraction of the "
             "matched reference protein -- catches reads that align well "
             "to just one domain of a multi-domain fusion protein (e.g. a "
             "shared substrate-binding domain fused onto an otherwise "
             "unrelated transporter). Unlike --min_query_cover, this is "
             "deliberately off by default: a short read is always much "
             "shorter than most full-length reference proteins, so a "
             "strict threshold here will also discard genuine short-read "
             "hits. Start low (0.1-0.3) and check the impact on recall "
             "before raising it.",
    )
    filt.add_argument(
        "--evalue", type=float, default=1e-5, metavar="FLOAT",
        help="Maximum e-value.",
    )
    filt.add_argument(
        "--threads", type=int, default=4, metavar="N",
        help="Number of threads for DIAMOND.",
    )


def _add_output_args(p: argparse.ArgumentParser) -> None:
    out = p.add_argument_group("output")
    out.add_argument(
        "--out_prefix", required=True, metavar="PREFIX",
        help="Prefix for output files (e.g. results/sample).",
    )
    out.add_argument(
        "--keep_aln", action="store_true", default=False,
        help="Retain the raw DIAMOND alignment TSV alongside outputs.",
    )


def _add_hpc_args(p: argparse.ArgumentParser) -> None:
    hpc = p.add_argument_group("HPC / performance")
    hpc.add_argument(
        "--tmpdir", default=None, metavar="DIR",
        help=(
            "Directory for temporary files (e.g. merged paired-end FASTQ). "
            "Defaults to $TMPDIR → $SCRATCH → /tmp.  "
            "Set to fast node-local scratch to avoid writing to a shared filesystem."
        ),
    )
    hpc.add_argument(
        "--total_reads", type=int, default=None, metavar="N",
        help=(
            "Total read count (skip on-the-fly counting). "
            "Provide when already known (e.g. from FastQC/MultiQC) "
            "to avoid iterating large files."
        ),
    )


# ---------------------------------------------------------------------------
# profile subcommand
# ---------------------------------------------------------------------------

def _add_profile_parser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "profile",
        help="Profile FASTQ reads against osmodiamond (diamond blastx).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "db", metavar="DIAMOND_DB",
        help="Path to the osmodiamond DIAMOND database (.dmnd).",
    )

    reads = p.add_argument_group("reads (paired-end or single-end)")
    me = reads.add_mutually_exclusive_group(required=True)
    me.add_argument(
        "-1", "--reads1", metavar="R1.FASTQ[.GZ]",
        help="Forward reads (R1) in FASTQ format (optionally gzip-compressed).",
    )
    me.add_argument(
        "--singles", metavar="READS.FASTQ[.GZ]",
        help="Single-end reads.",
    )
    reads.add_argument(
        "-2", "--reads2", metavar="R2.FASTQ[.GZ]", default=None,
        help="Reverse reads (R2).  Only valid with -1/--reads1.",
    )
    reads.add_argument(
        "--min_seqlen", type=int, default=20, metavar="AA",
        help="Minimum aligned ORF length in AMINO ACIDS, passed to diamond "
             "--min-orf (NOT bp/nucleotides). Default (20 aa = 60bp) is safe "
             "for typical 100-150bp short reads; a too-high value silently "
             "drops all short-read hits since a short read can only ever "
             "contain a short partial ORF.",
    )

    cascade = p.add_argument_group("DIAMOND+HMM cascade (optional)")
    cascade.add_argument(
        "--cascade_config", default=None, metavar="OSMO_REFDB.PROFILE_CASCADE.TSV",
        help="Path to osmo_refdb's <release>.profile_cascade.tsv. When "
             "given (with --hmm_db), reads whose DIAMOND call lands on a "
             "family flagged there (real, demonstrated precision "
             "problems, e.g. betL vs betT/CaiT) get a second opinion from "
             "hmmscan using a short-read-calibrated threshold before the "
             "call is kept -- only that flagged subset pays the extra "
             "cost, not the whole dataset. DIAMOND's GA-cutoff-equivalent "
             "specificity gate this project has otherwise lacked in "
             "`profile` mode. Requires orfm + hmmscan on PATH.",
    )
    cascade.add_argument(
        "--hmm_db", default=None, metavar="OSMO_REFDB.HMM",
        help="Path to the pressed osmo_refdb HMM database. Required when "
             "--cascade_config is given.",
    )
    cascade.add_argument(
        "--exclude_families", default=None, metavar="OSMO_REFDB.PROFILE_EXCLUDED_FAMILIES.TXT",
        help="Path to osmo_refdb's <release>.profile_excluded_families.txt "
             "(families marked scope: annotate_only, e.g. murB -- a "
             "near-universal housekeeping gene with no osmoadaptation-"
             "specific signal at the read level). Those families are "
             "still searched normally; this only drops them from the "
             "reported gene_counts.tsv.",
    )

    _add_filter_args(p)
    _add_output_args(p)
    _add_hpc_args(p)


# ---------------------------------------------------------------------------
# annotate subcommand
# ---------------------------------------------------------------------------

def _add_annotate_parser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "annotate",
        help="Annotate an assembly against osmodiamond (Prodigal + diamond blastp).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "db", metavar="DIAMOND_DB",
        help="Path to the osmodiamond DIAMOND database (.dmnd).",
    )
    p.add_argument(
        "assembly", metavar="ASSEMBLY.FASTA",
        help="Assembled contigs in FASTA format.",
    )
    p.add_argument(
        "--proteins", default=None, metavar="PROTEINS.FAA",
        help="Pre-called protein FASTA (skip Prodigal if provided).",
    )
    p.add_argument(
        "--method", choices=["diamond", "hmm", "both"], default="hmm",
        help="Detection method. 'hmm' (default) runs hmmscan --cut_ga "
             "against --hmm_db using each family's calibrated "
             "gathering-threshold cutoff -- meaningful here (unlike in "
             "`profile`) because Prodigal-called ORFs are full-length, the "
             "same scale those cutoffs were calibrated against. Validated "
             "against a real E. coli K-12 genome to be meaningfully more "
             "specific than DIAMOND for this reference database (fixed a "
             "false positive and an over-counted family DIAMOND got "
             "wrong). Requires --hmm_db and hmmscan on PATH -- pass "
             "'diamond' instead if you want to avoid the HMMER dependency, "
             "or 'both' to get DIAMOND's output alongside HMM's for "
             "comparison.",
    )
    p.add_argument(
        "--hmm_db", default=None, metavar="OSMO_REFDB.HMM",
        help="Path to the pressed osmo_refdb HMM database "
             "(hmms/osmo_refdb.hmm, alongside its .h3f/.h3i/.h3m/.h3p "
             "press indices). Required when --method is 'hmm' (the "
             "default) or 'both'; pass --method diamond if you don't have "
             "this file.",
    )
    p.add_argument(
        "--diamond_cutoffs", default=None, metavar="OSMO_REFDB.DIAMOND_CUTOFFS.TSV",
        help="Path to osmo_refdb's <release>.diamond_cutoffs.tsv (from "
             "08b_calibrate_diamond_cutoffs.py). Gives DIAMOND calls a "
             "per-family specificity gate analogous to HMM's GA cutoff -- "
             "without it, a single flat --min_identity applies to every "
             "family, which can't separate a true gene from a genuinely "
             "close paralog (e.g. betL vs betT/CaiT) above that threshold. "
             "Only offered here, not in `profile`: these cutoffs are "
             "calibrated against full-length protein-vs-protein "
             "alignments, the same scale as Prodigal ORFs -- applying them "
             "to short read fragments would hit the same scale mismatch "
             "problem HMM's GA cutoffs have on short reads.",
    )
    p.add_argument(
        "--keep_proteins", action="store_true", default=False,
        help="Keep the Prodigal protein FASTA (ignored when --proteins is given).",
    )
    p.add_argument(
        "--exclude_families", default=None, metavar="OSMO_REFDB.ANNOTATE_EXCLUDED_FAMILIES.TXT",
        help="Path to osmo_refdb's <release>.annotate_excluded_families.txt "
             "(currently just decoy references, e.g. betL_decoy -- "
             "families.yaml: decoy_from_negatives -- which exist purely to "
             "win DIAMOND's best-hit contest away from a mislabeled call "
             "and must never appear as a reported gene family). Those "
             "sequences are still searched normally; this only drops them "
             "from the reported gene_counts.tsv. Note this is a different "
             "file from `profile`'s --exclude_families: annotate_only "
             "families like murB stay visible here.",
    )

    _add_filter_args(p)
    _add_output_args(p)
    _add_hpc_args(p)
