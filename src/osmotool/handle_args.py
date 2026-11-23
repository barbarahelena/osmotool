"""
handle_args.py — argument parsing for osmotool
"""

from __future__ import annotations

import argparse
from importlib.metadata import PackageNotFoundError, version

from osmotool.download import ANNOTATE_ONLY_RELEASES, RELEASES, LATEST_RELEASE

try:
    __version__ = version("osmotool")
except PackageNotFoundError:
    __version__ = "unknown"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="osmotool",
        description=(
            "Screen osmoadaptation genes in metagenomic datasets using DIAMOND "
            "and HMMER. The gene panel itself (which families are searched) is "
            "defined entirely by the osmo_refdb release passed as DATABASE -- "
            "see that release's families.yaml / README for the current panel, "
            "rather than a fixed list baked into osmotool."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Enable debug logging."
    )

    sub = parser.add_subparsers(dest="subcommand", metavar="subcommand")
    sub.required = True

    _add_profile_parser(sub)
    _add_annotate_parser(sub)
    _add_predict_nacl_parser(sub)
    _add_build_catalog_parser(sub)
    _add_index_catalog_parser(sub)
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


# ---------------------------------------------------------------------------
# profile subcommand
# ---------------------------------------------------------------------------

def _add_profile_parser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "profile",
        help="Profile FASTQ reads: gene catalog first, translated search for the rest.",
        description="Profile osmoadaptation genes in FASTQ reads, in two steps. "
                    "1) Reads are assigned against the release's gene catalog, "
                    "annotated once at full length with the release's HMMs: a "
                    "read matching a catalog gene near-exactly takes that "
                    "gene's family. 2) Whatever the catalog cannot settle falls "
                    "through to translated search (diamond blastx, then the "
                    "DIAMOND+HMM cascade for families DIAMOND confuses). "
                    "Needs an osmo_refdb release that ships a catalog (v10 or "
                    "later).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "database", metavar="DATABASE",
        help="Path to an unpacked osmo_refdb release directory (v10 or "
             "later, e.g. releases/v10/). osmotool finds the gene catalog "
             "under catalogs/, osmo_refdb.dmnd, hmms/osmo_refdb.hmm, "
             "osmo_refdb.profile_cascade.tsv, "
             "osmo_refdb.profile_excluded_families.txt and "
             "osmo_refdb.profile_report_levels.tsv inside it by their "
             "fixed names -- no need to point at each file separately. With "
             "the report levels, <out_prefix>.unit_counts.tsv is written "
             "next to gene_counts.tsv: one row per reporting unit (a gene, "
             "or a paralog group short reads can't split), with units not "
             "reliable at any level listed as excluded. The levels were "
             "decided on translated-search calibration reads, so they "
             "describe step 2; catalog-settled reads are full-length calls.",
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
    total = reads.add_mutually_exclusive_group()
    total.add_argument(
        "--total_fragments", type=int, default=None, metavar="N",
        help="Number of sequenced fragments, the RPM denominator: read PAIRS "
             "for paired-end (the read count of R1 alone), reads for "
             "single-end. Skips counting the input; provide when already "
             "known (e.g. from FastQC/MultiQC) to avoid iterating large files.",
    )
    total.add_argument(
        "--total_reads", type=int, default=None, metavar="N",
        help="Like --total_fragments, but the reads summed over ALL input "
             "files (R1 + R2 for paired-end); halved for paired-end input.",
    )
    reads.add_argument(
        "--min_seqlen", type=int, default=20, metavar="AA",
        help="Minimum aligned ORF length in AMINO ACIDS, passed to diamond "
             "--min-orf (NOT bp/nucleotides). Default (20 aa = 60bp) is safe "
             "for typical 100-150bp short reads; a too-high value silently "
             "drops all short-read hits since a short read can only ever "
             "contain a short partial ORF.",
    )

    cat = p.add_argument_group(
        "gene catalog",
        "Reads are first assigned against the gene catalog DATABASE ships, "
        "annotated at full length with DATABASE's HMMs; whatever the catalog "
        "cannot settle falls through to translated search (DIAMOND + the "
        "cascade).")
    cat.add_argument(
        "--habitat", metavar="NAME", default=None,
        help="Which of DATABASE's catalogs (catalogs/NAME/) to assign reads "
             "against. Needed only when the release ships more than one; "
             "osmo_refdb v10 ships gmgc_gut alone.",
    )
    # A catalog built locally with `osmotool build-catalog` for DATABASE's
    # HMMs, used instead of the one DATABASE ships. Kept out of the main help
    # (documented in docs/profile.md): a shipped catalog is versioned with its
    # release, a local one is the user's to keep in step.
    cat.add_argument("--catalog_dir", metavar="DIR", default=None, help=argparse.SUPPRESS)
    cat.add_argument(
        "--catalog_only", action="store_true", default=False,
        help="Count catalog-settled reads only; skip translated search for "
             "the rest. Only sensible when the community is well covered by "
             "the catalog (e.g. human gut with a gut catalog).",
    )
    cat.add_argument(
        "--catalog_index", metavar="DIR", default=None,
        help="Map against this full index of the catalog's genes (built by "
             "`osmotool index-catalog`, distributed separately from the "
             "release) instead of the index the catalog ships. Must be built "
             "from the same gene file the catalog was annotated from.",
    )
    cat.add_argument(
        "--aligner_bin", metavar="EXE", default=None,
        help="Map with this executable instead of the one recorded in the "
             "catalog (same CLI and index format, e.g. a reimplementation).",
    )

    tune = p.add_argument_group("catalog tuning")
    tune.add_argument(
        "--catalog_min_identity", type=float, default=0.97, metavar="FRAC",
        help="Minimum nucleotide identity (1 - NM / alignment columns) for a "
             "read to be settled by the catalog.",
    )
    tune.add_argument(
        "--catalog_min_aligned", type=int, default=45, metavar="BP",
        help="Minimum aligned read bases for a read to be settled by the catalog "
             "(soft-clipped ends, e.g. a read running off a gene end, don't count).",
    )

    cascade = p.add_argument_group("translated search: DIAMOND+HMM cascade")
    cascade.add_argument(
        "--no_cascade", action="store_true", default=False,
        help="Disable the DIAMOND+HMM cascade (on by default). The cascade "
             "gives reads a second opinion from hmmsearch when their DIAMOND "
             "call lands on a family with demonstrated precision problems "
             "(e.g. betL vs betT/CaiT, or several others as low as "
             "5-25%% DIAMOND-only precision -- see osmo_refdb's "
             "profile_cascade.tsv), using DATABASE's "
             "osmo_refdb.profile_cascade.tsv and hmms/osmo_refdb.hmm -- "
             "only that flagged subset pays the extra cost, not the whole "
             "dataset. Requires orfm + hmmsearch on PATH, and both files "
             "present in DATABASE; pass this flag to run DIAMOND-only if "
             "you don't have them.",
    )

    norm = p.add_argument_group("copies-per-genome normalisation")
    norm.add_argument(
        "--no_markers", action="store_true", default=False,
        help="Skip copies-per-genome normalisation. By default, when "
             "DATABASE ships a markers/ subdirectory (universal single-copy "
             "marker genes), reads are also searched against it in a second "
             "DIAMOND pass and gene_counts.tsv gains a copies_per_genome "
             "column (family reads per kb / mean marker reads per kb, "
             "MicrobeCensus-style MAD outlier removal), plus "
             "<out_prefix>.markers.tsv with the per-marker evidence. "
             "Releases without markers/ are unaffected.",
    )

    _add_filter_args(p)
    _add_output_args(p)
    _add_hpc_args(p)

    bench = p.add_argument_group("benchmarking")
    bench.add_argument(
        "--translated_only", action="store_true", default=False,
        help="Skip the catalog and assign every read by translated search "
             "(DIAMOND + the cascade), for comparing the two. Also runs on a "
             "release that ships no catalog.",
    )


# ---------------------------------------------------------------------------
# download-db subcommand
# ---------------------------------------------------------------------------

def _add_download_parser(sub: argparse._SubParsersAction) -> None:
    known = ", ".join(f"{r} (annotate only)" if r in ANNOTATE_ONLY_RELEASES else r
                      for r in sorted(RELEASES, key=lambda r: int(r.lstrip("v"))))
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
             "<location>/v10/) -- pass that path as osmotool's DATABASE "
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

def _add_predict_nacl_parser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "predict-nacl",
        help="Predict maximum NaCl tolerance from existing `annotate` output.",
        description="Predict each genome's maximum NaCl tolerance (% w/v) from "
                    "the output of an earlier `osmotool annotate --method hmm` "
                    "run (osmotool >= 0.8.0) against the same DATABASE, with "
                    "the model in DATABASE's nacl_model/. Writes "
                    "<prefix>.nacl_prediction.tsv and "
                    "<prefix>.nacl_contributions.tsv next to each prefix's "
                    "annotate output, and prints one summary table. Use "
                    "`annotate --predict_nacl` to do both in one run.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "database", metavar="DATABASE",
        help="The osmo_refdb release directory the genomes were annotated with "
             "(it must contain nacl_model/).",
    )
    p.add_argument(
        "prefixes", metavar="PREFIX", nargs="+",
        help="annotate --out_prefix of each genome (so PREFIX.gene_counts.tsv "
             "exists). For --method both output, pass PREFIX.hmm.",
    )
    p.add_argument(
        "-o", "--output", default=None, metavar="TSV",
        help="Write the summary table here instead of stdout.",
    )


def _add_build_catalog_parser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "build-catalog",
        help="Annotate a gene catalog once at full length for `profile --catalog_dir`.",
        description="Translate a nucleotide gene catalog (e.g. a prevalence-"
                    "filtered GMGC subset, or UHGG genes), search it with "
                    "DATABASE's HMMs under their gathering cutoffs (the decision "
                    "`annotate --method hmm` makes), and write a catalog "
                    "directory with only the genes that hit any model: annotated "
                    "genes, plus sub-cutoff homologs and decoy-model genes kept so "
                    "that reads from them are recognised and dropped. The whole "
                    "catalog is indexed (bwa-mem3 by default), so every read maps "
                    "competitively to its own gene.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("database", metavar="DATABASE",
                   help="osmo_refdb release directory (needs hmms/osmo_refdb.hmm).")
    p.add_argument("--genes", required=True, metavar="FNA[.GZ]",
                   help="Nucleotide CDS FASTA, one in-frame gene per record.")
    p.add_argument("--proteins", default=None, metavar="FAA[.GZ]",
                   help="Protein FASTA with the same IDs (otherwise genes are "
                        "translated with table 11).")
    p.add_argument("-o", "--out_dir", required=True, metavar="DIR")
    p.add_argument("--aligner", default="bwa-mem3", choices=["bwa-mem3", "bwa-mem2"],
                   help="Aligner that indexes the catalog; profile maps with the same one. "
                        "bwa-mem3 needs AVX2 on x86-64 (or ARM NEON).")
    p.add_argument("--aligner_bin", default=None, metavar="EXE",
                   help="Executable for --aligner if not the default name on PATH.")
    p.add_argument("--index", default="all", choices=["all", "hits", "neighbours", "none"],
                   help="all: index every catalog gene, so reads map competitively. "
                        "hits: index only genes with an HMM hit -- a much "
                        "smaller index, but reads from other genes can only fall back. "
                        "neighbours: the hit genes plus every gene without a hit that a "
                        "read window aligns to a hit gene at >= 90%% over >= 45 bp, so "
                        "reads from those competitors map to their own gene and are "
                        "dropped; about as small as hits. "
                        "none: annotation only, used with a full index from "
                        "`osmotool index-catalog` (profile --catalog_index).")
    p.add_argument("--min_protein_length", type=int, default=30, metavar="AA")
    p.add_argument("--evalue", type=float, default=1e-5,
                   help="hmmsearch reporting threshold; reported genes below the "
                        "gathering cutoff are kept as homologs.")
    p.add_argument("-t", "--threads", type=int, default=4)


def _add_index_catalog_parser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "index-catalog",
        help="Build the full, competitive index of a gene catalog (release-independent).",
        description="Index every gene of a catalog with bwa-mem3 (default) for "
                    "`profile --catalog_index DIR`. The index depends "
                    "only on the gene file and the aligner, not on an osmo_refdb "
                    "release, so it is built and distributed once and reused across "
                    "releases; profile checks it was built from the same gene file "
                    "the release's catalog was annotated from.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--genes", required=True, metavar="FNA[.GZ]",
                   help="The nucleotide gene FASTA the catalog was annotated from.")
    p.add_argument("-o", "--out_dir", required=True, metavar="DIR")
    p.add_argument("--aligner", default="bwa-mem3", choices=["bwa-mem3", "bwa-mem2"])
    p.add_argument("--aligner_bin", default=None, metavar="EXE")


def _add_annotate_parser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "annotate",
        help="Annotate an assembly: Prodigal, then hmmsearch (default) or diamond blastp.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "database", metavar="DATABASE",
        help="Path to an unpacked osmo_refdb release directory (e.g. "
             "releases/v10/). osmotool finds osmo_refdb.dmnd, "
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
        help="Detection method. 'hmm' (default) runs hmmsearch --cut_ga "
             "against DATABASE's hmms/osmo_refdb.hmm using each family's "
             "calibrated gathering-threshold cutoff -- meaningful here "
             "(unlike in `profile`) because Prodigal-called ORFs are "
             "full-length, the same scale those cutoffs were calibrated "
             "against. Validated against a real E. coli K-12 genome to be "
             "meaningfully more specific than DIAMOND for this reference "
             "database (fixed a false positive and an over-counted family "
             "DIAMOND got wrong). Requires hmms/osmo_refdb.hmm present in "
             "DATABASE and hmmsearch on PATH -- pass 'diamond' instead if "
             "you want to avoid the HMMER dependency, or 'both' to get "
             "DIAMOND's output alongside HMM's for comparison.",
    )
    p.add_argument(
        "--keep_proteins", action="store_true", default=False,
        help="Keep the Prodigal protein FASTA (ignored when --proteins is given).",
    )

    systems = p.add_argument_group("multi-subunit system detection")
    systems.add_argument(
        "--skip_systems", action="store_true", default=False,
        help="Skip the complex/system-aware detection layer (e.g. opuA, mrp) "
             "-- no <out_prefix>.systems.tsv, no ORF-coordinate parsing. Off "
             "by default.",
    )
    systems.add_argument(
        "--systems_config", default=None, metavar="TSV",
        help="Optional TSV overriding each system's max_distance_genes/"
             "max_distance_bp proximity window (columns: system, "
             "max_distance_genes, max_distance_bp). The defaults are 3 genes / "
             "2000 bp between consecutive subunits (30 / 35000 for "
             "galE_murB_mazG), checked against operons in six species -- "
             "see docs/annotate.md. Mandatory/accessory/fusion membership "
             "is not overridable this way.",
    )

    nacl = p.add_argument_group("NaCl tolerance prediction")
    nacl.add_argument(
        "--predict_nacl", action="store_true", default=False,
        help="Also predict the genome's maximum NaCl tolerance (%% w/v) with the "
             "model in DATABASE's nacl_model/ (osmo_refdb v10 and later), from "
             "the HMM counts. Writes <out_prefix>.nacl_prediction.tsv (value, 80%% "
             "interval, flags, largest contributions) and "
             "<out_prefix>.nacl_contributions.tsv (every feature's SHAP "
             "contribution). Needs --method hmm or both, and the xgboost package "
             "(pip install 'osmotool[predict]'). Salt-out organisms only: genomes "
             "that look salt-in are flagged.",
    )

    _add_filter_args(p)
    _add_output_args(p)
    _add_hpc_args(p)
