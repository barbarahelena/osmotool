"""
handle_args.py — argument parsing for osmotool
"""

from __future__ import annotations

import argparse

from osmotool.download import RELEASES, LATEST_RELEASE


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="osmotool",
        description=(
            "Screen osmoadaptation genes in metagenomic datasets using DIAMOND "
            "and HMMER. "
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
    _add_download_parser(sub)

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
        "database", metavar="DATABASE",
        help="Path to an unpacked osmo_refdb release directory (e.g. "
             "releases/v5/). osmotool finds osmo_refdb.dmnd, "
             "hmms/osmo_refdb.hmm, osmo_refdb.profile_cascade.tsv, and "
             "osmo_refdb.profile_excluded_families.txt inside it by their "
             "fixed names -- no need to point at each file separately.",
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
        "--cascade", action="store_true", default=False,
        help="Give reads a second opinion from hmmscan when their DIAMOND "
             "call lands on a family with demonstrated precision problems "
             "(e.g. betL vs betT/CaiT), using DATABASE's "
             "osmo_refdb.profile_cascade.tsv and hmms/osmo_refdb.hmm -- "
             "only that flagged subset pays the extra cost, not the whole "
             "dataset. DIAMOND's GA-cutoff-equivalent specificity gate "
             "this project has otherwise lacked in `profile` mode. "
             "Requires orfm + hmmscan on PATH, and both files present in "
             "DATABASE.",
    )

    _add_filter_args(p)
    _add_output_args(p)
    _add_hpc_args(p)


# ---------------------------------------------------------------------------
# download-db subcommand
# ---------------------------------------------------------------------------

def _add_download_parser(sub: argparse._SubParsersAction) -> None:
    known = ", ".join(sorted(RELEASES))
    p = sub.add_parser(
        "download-db",
        help="Download and unpack an osmo_refdb release from Zenodo.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--release", default="latest", metavar="NAME",
        help=f"osmo_refdb release to download, or 'latest' (currently "
             f"{LATEST_RELEASE}). Known releases: {known}.",
    )
    p.add_argument(
        "--location", default=".", metavar="DIR",
        help="Directory to download and unpack into. The release unpacks "
             "into a subdirectory named after the release (e.g. "
             "<location>/v5/) -- pass that path as osmotool's DATABASE "
             "argument.",
    )
    p.add_argument(
        "--force", action="store_true", default=False,
        help="Re-download and overwrite even if the release already exists at --location.",
    )
    p.add_argument(
        "--keep_archive", action="store_true", default=False,
        help="Keep the downloaded .tar.gz alongside the extracted directory "
             "instead of deleting it after extraction.",
    )


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
        "database", metavar="DATABASE",
        help="Path to an unpacked osmo_refdb release directory (e.g. "
             "releases/v5/). osmotool finds osmo_refdb.dmnd, "
             "hmms/osmo_refdb.hmm, osmo_refdb.diamond_cutoffs.tsv, and "
             "osmo_refdb.annotate_excluded_families.txt inside it by "
             "their fixed names -- no need to point at each file "
             "separately.",
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
             "against DATABASE's hmms/osmo_refdb.hmm using each family's "
             "calibrated gathering-threshold cutoff -- meaningful here "
             "(unlike in `profile`) because Prodigal-called ORFs are "
             "full-length, the same scale those cutoffs were calibrated "
             "against. Validated against a real E. coli K-12 genome to be "
             "meaningfully more specific than DIAMOND for this reference "
             "database (fixed a false positive and an over-counted family "
             "DIAMOND got wrong). Requires hmms/osmo_refdb.hmm present in "
             "DATABASE and hmmscan on PATH -- pass 'diamond' instead if "
             "you want to avoid the HMMER dependency, or 'both' to get "
             "DIAMOND's output alongside HMM's for comparison.",
    )
    p.add_argument(
        "--keep_proteins", action="store_true", default=False,
        help="Keep the Prodigal protein FASTA (ignored when --proteins is given).",
    )

    _add_filter_args(p)
    _add_output_args(p)
    _add_hpc_args(p)
