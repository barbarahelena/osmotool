"""
functions.py — DIAMOND wrappers and alignment parsing for osmotool
"""

from __future__ import annotations

import csv
import logging
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from osmotool.utils import get_tmpdir

log = logging.getLogger("osmotool")

# DIAMOND tabular output columns used throughout
#
# NOTE: scovhsp (subject/target coverage) is requested here as an *output*
# column only -- it is NOT passed to diamond as a --subject-cover CLI
# filter alongside --query-cover. Combining --query-cover and
# --subject-cover in the same diamond invocation triggers an unresolved
# upstream crash (vector bounds exception) in some diamond builds/modes
# (https://github.com/bbuchfink/diamond/issues/812). Subject-coverage
# filtering is instead done in Python in parse_diamond_output() below.
DIAMOND_FIELDS = (
    "qseqid", "sseqid", "pident", "qcovhsp", "scovhsp", "length", "evalue", "bitscore"
)
DIAMOND_OUTFMT_ARGS = ["6", *DIAMOND_FIELDS]


# ---------------------------------------------------------------------------
# DIAMOND blastx  (reads → protein DB)
# ---------------------------------------------------------------------------

def run_diamond_blastx(
    db: str | Path,
    out_file: str | Path,
    *,
    reads1: str | Path | None = None,
    reads2: str | Path | None = None,
    singles: str | Path | None = None,
    min_identity: float = 0.80,
    min_query_cover: float = 0.80,
    min_seqlen: int = 20,
    evalue: float = 1e-5,
    threads: int = 4,
    tmpdir: str | Path | None = None,
    extra_args: list[str] | None = None,
) -> Path:
    """
    Run ``diamond blastx`` on paired-end or single-end FASTQ reads.

    For paired-end input both files are concatenated into a temp file
    in *tmpdir* (defaults to ``$TMPDIR`` / node-local scratch on HPC).

    NOTE: ``min_seqlen`` is passed to DIAMOND's ``--min-orf`` and is in
    AMINO ACIDS (not nucleotides/bp), despite what the option name might
    suggest. A too-high value silently drops every alignment for short
    reads: e.g. a 100-150bp Illumina read can contain at most ~33-50 aa
    of ORF even with a perfect in-frame match and no flanking sequence,
    so a default of 45 aa (used prior to this fix) rejected essentially
    all short-read hits. Default here (20 aa = 60bp) is safe for typical
    100-150bp short reads; raise it only if you are profiling long reads
    or assembled contigs where partial/spurious short ORFs are a bigger
    concern than losing true short-read hits.

    Returns the path to the tabular output file.
    """
    out_file = Path(out_file)
    db = Path(db)

    if not db.exists() and not Path(str(db) + ".dmnd").exists():
        raise FileNotFoundError(f"DIAMOND database not found: {db}")

    query_path, query_is_temp = _prepare_query(
        reads1=reads1, reads2=reads2, singles=singles, tmpdir=tmpdir
    )

    cmd = [
        "diamond", "blastx",
        "--db", str(db),
        "--query", str(query_path),
        "--out", str(out_file),
        "--outfmt", *DIAMOND_OUTFMT_ARGS,
        "--id", str(min_identity * 100),       # diamond takes percentage
        "--query-cover", str(min_query_cover * 100),
        "--min-score", "0",
        "--evalue", str(evalue),
        "--min-orf", str(min_seqlen),
        "--threads", str(threads),
        "--quiet",
    ]
    if extra_args:
        cmd.extend(extra_args)

    log.info("Running diamond blastx → %s", out_file)
    log.debug("Command: %s", " ".join(cmd))

    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise RuntimeError(f"diamond blastx failed (exit {result.returncode})")

    if query_is_temp:
        query_path.unlink(missing_ok=True)

    return out_file


def _prepare_query(
    reads1: str | Path | None,
    reads2: str | Path | None,
    singles: str | Path | None,
    tmpdir: str | Path | None = None,
) -> tuple[Path, bool]:
    """
    Return (query_path, is_temporary).

    If only one file is given, return it directly (is_temporary=False).
    For paired-end, merge R1 + R2 into one temp FASTQ under *tmpdir*
    (resolved via :func:`~osmotool.utils.get_tmpdir`, honouring ``$TMPDIR``
    for HPC node-local scratch) via _merge_paired_with_mate_tags, which
    tags each read ID so the two mates can never collide under DIAMOND's
    qseqid.
    """
    if singles is not None:
        return Path(singles), False
    if reads1 is not None and reads2 is None:
        return Path(reads1), False
    if reads1 is not None and reads2 is not None:
        scratch = get_tmpdir(tmpdir)
        tmp = tempfile.NamedTemporaryFile(
            suffix=".fastq.gz", delete=False, prefix="osmotool_",
            dir=scratch,
        )
        tmp.close()
        merged = Path(tmp.name)
        log.debug("Merging paired reads → %s (tmpdir=%s)", merged, scratch)
        _merge_paired_with_mate_tags(Path(reads1), Path(reads2), merged)
        return merged, True
    raise ValueError("Provide either --singles or -1/--reads1 (with optional -2/--reads2).")


def _merge_paired_with_mate_tags(reads1: Path, reads2: Path, dest: Path) -> None:
    """
    Merge R1 + R2 into one gzipped FASTQ for a single DIAMOND invocation,
    appending _R1/_R2 onto each read ID so the two mates of a pair can
    never collide under DIAMOND's qseqid (the header up to the first
    whitespace).

    This matters for modern Illumina headers, which put the mate marker
    AFTER that whitespace (e.g. "... 1:N:0:INDEX" / "... 2:N:0:INDEX") --
    without tagging, both mates of a fragment would produce an identical
    qseqid, and quantifier.select_best_hits()'s dict (keyed by qseqid)
    would silently keep only the higher-scoring mate's hit, discarding the
    other mate's 0.5 count entirely instead of reporting both (see
    quantifier.count_hits' paired-weighting docstring). Older
    "READID/1"/"READID/2"-style headers already avoid this by luck, not by
    design -- tagging here makes correctness independent of header
    convention rather than relying on it.
    """
    import gzip

    with gzip.open(dest, "wt", compresslevel=1) as out_fh:
        for path, tag in ((reads1, "_R1"), (reads2, "_R2")):
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rt") as in_fh:
                for i, line in enumerate(in_fh):
                    if i % 4 != 0:
                        out_fh.write(line)
                        continue
                    if not line.startswith("@"):
                        raise ValueError(
                            f"Malformed FASTQ header in {path} "
                            f"(record starting at line {i + 1}): {line!r}"
                        )
                    read_id, sep, rest = line[1:].rstrip("\n").partition(" ")
                    header = f"@{read_id}{tag}"
                    if sep:
                        header += f" {rest}"
                    out_fh.write(header + "\n")


# ---------------------------------------------------------------------------
# DIAMOND blastp  (proteins → protein DB)
# ---------------------------------------------------------------------------

def run_diamond_blastp(
    db: str | Path,
    proteins: str | Path,
    out_file: str | Path,
    *,
    min_identity: float = 0.80,
    min_query_cover: float = 0.80,
    evalue: float = 1e-5,
    threads: int = 4,
    extra_args: list[str] | None = None,
) -> Path:
    """
    Run ``diamond blastp`` on a protein FASTA (e.g. Prodigal output).

    Returns the path to the tabular output file.
    """
    out_file = Path(out_file)
    db = Path(db)

    cmd = [
        "diamond", "blastp",
        "--db", str(db),
        "--query", str(proteins),
        "--out", str(out_file),
        "--outfmt", *DIAMOND_OUTFMT_ARGS,
        "--id", str(min_identity * 100),
        "--query-cover", str(min_query_cover * 100),
        "--evalue", str(evalue),
        "--threads", str(threads),
        "--no-self-hits",
        "--quiet",
    ]
    if extra_args:
        cmd.extend(extra_args)

    log.info("Running diamond blastp → %s", out_file)
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise RuntimeError(f"diamond blastp failed (exit {result.returncode})")

    return out_file


# ---------------------------------------------------------------------------
# hmmsearch  (osmo_refdb HMM database → proteins)
# ---------------------------------------------------------------------------
#
# --cut_ga is only applied in `annotate` mode: osmo_refdb's HMM GA cutoffs
# are calibrated against full-length protein sequences, so applying them
# to short translated read fragments would silently produce zero hits for
# essentially every read (see osmo_refdb's
# pipeline/translate_and_hmmscan.sh for the same pitfall it works around
# by sweeping raw bit scores instead). Prodigal-called ORFs in `annotate`
# ARE full-length, so --cut_ga is exactly the right tool there: it gives
# HMM a per-family specificity gate that DIAMOND's flat --min_identity
# does not have, which matters for genuinely close paralogs within the
# same Pfam family (e.g. betL vs betT/CaiT). `profile`'s cascade runs
# without --cut_ga and applies its own short-read-calibrated threshold.
#
# hmmsearch rather than hmmscan: osmo_refdb is a few dozen models searched
# against thousands (annotate) to millions (profile cascade) of
# sequences. hmmscan re-reads the whole profile database for every query
# sequence; hmmsearch loads each model once and streams the sequences,
# which is ~5x faster on a single bacterial genome. Bit scores -- and so
# --cut_ga and the cascade's bit-score thresholds -- are identical between
# the two. E-values are not by default: hmmscan's search space Z is the
# number of models, hmmsearch's is the number of sequences. run_hmmsearch
# passes -Z/--domZ = number of models so E-values (used for best-hit
# tie-breaking, and for HMMER's default E<=10 reporting threshold, which
# decides whether weak raw-score cascade hits are reported at all) stay
# exactly what hmmscan produced, independent of input size.

def count_hmm_models(hmm_db: str | Path) -> int:
    """Number of profile HMMs in an HMMER3 .hmm file (one NAME line each)."""
    n = 0
    with open(hmm_db) as fh:
        for line in fh:
            if line.startswith("NAME "):
                n += 1
    if n == 0:
        raise ValueError(f"No HMMs found in {hmm_db} -- not an HMMER3 .hmm file?")
    return n


def run_hmmsearch(
    hmm_db: str | Path,
    proteins: str | Path,
    out_file: str | Path,
    *,
    threads: int = 4,
    use_cut_ga: bool = True,
) -> Path:
    """
    Run ``hmmsearch`` with the osmo_refdb HMM database (the .hmm file; the
    .h3* hmmpress indices are not needed) against a protein FASTA.

    use_cut_ga=True (default, used by `annotate`) applies each family's
    calibrated gathering-threshold cutoff, valid for full-length ORFs.
    use_cut_ga=False (used by `profile`'s DIAMOND+HMM cascade) reports raw
    bit scores instead -- short read fragments can never reach a cutoff
    calibrated against full-length sequences, so the cascade instead
    applies its own short-read-calibrated threshold (see
    quantifier.apply_cascade_check) to whatever parse_hmmsearch_output()
    reports here.

    -Z/--domZ are fixed to the number of models in *hmm_db* so E-values
    match hmmscan's (see the comment block above).

    Returns the path to the tblout output file.
    """
    out_file = Path(out_file)
    hmm_db = Path(hmm_db)

    if not hmm_db.exists():
        raise FileNotFoundError(f"HMM database not found: {hmm_db}")

    n_models = str(count_hmm_models(hmm_db))
    cmd = [
        "hmmsearch",
        *(["--cut_ga"] if use_cut_ga else []),
        "-Z", n_models,
        "--domZ", n_models,
        "--tblout", str(out_file),
        "--noali",
        "--cpu", str(threads),
        str(hmm_db),
        str(proteins),
    ]

    log.info("Running hmmsearch → %s", out_file)
    log.debug("Command: %s", " ".join(cmd))

    result = subprocess.run(cmd, stdout=subprocess.DEVNULL)
    if result.returncode != 0:
        raise RuntimeError(f"hmmsearch failed (exit {result.returncode})")

    return out_file


def parse_hmmsearch_output(tblout_file: str | Path) -> Iterator[dict]:
    """
    Yield one dict per hit row from an hmmsearch --tblout file, shaped to
    match parse_diamond_output()'s row format (qseqid/sseqid/bitscore/
    evalue) so the same quantifier.select_best_hits()/count_hits()
    functions work unchanged for either method.

    hmmsearch tblout columns (whitespace-delimited):
      target name  accession  query name  accession  E-value  score  bias  ...
    In hmmsearch the *target* is the protein sequence and the *query* is
    the HMM (the reverse of hmmscan), so qseqid = column 0 and sseqid =
    column 2. osmo_refdb builds one HMM per family, named after the family
    itself (e.g. "betL"), so the HMM name is used directly as sseqid --
    gene_family_from_header() splits on '|' and a bare family name with
    no '|' passes through unchanged.

    If hmmsearch was run with --cut_ga (run_hmmsearch's default), every
    row here already passed that family's calibrated gathering threshold
    -- no further score filtering is needed. If run with use_cut_ga=False
    (profile mode's cascade), these are raw scores and the caller is
    responsible for applying its own threshold.
    """
    with open(tblout_file) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 6:
                continue
            yield {
                "qseqid":   fields[0],
                "sseqid":   fields[2],
                "evalue":   float(fields[4]),
                "bitscore": float(fields[5]),
            }


# orfm's ORF headers are exactly "<original read ID>_<orf>_<frame>_<n>"
# (same convention osmo_refdb's own 11_compute_metrics.py relies on and
# has verified against real orfm output), so stripping the trailing three
# underscore-separated numbers recovers the original read ID an ORF was
# called from. The read ID may or may not carry a "/1"/"/2" mate marker:
# single-end reads and modern Illumina headers (mate marker after the
# first whitespace, which orfm drops) don't, so the marker must not be
# required here -- otherwise no ORF maps back to its read and the cascade
# silently drops every flagged-family read. The greedy (.*) strips exactly
# the last three numeric fields, so read IDs that themselves contain
# underscores are preserved.
ORF_SUFFIX_RE = re.compile(r"^(.*)_\d+_\d+_\d+$")


def strip_orf_suffix(orf_id: str) -> str:
    """Recover the original read ID from an orfm-called ORF ID (see
    ORF_SUFFIX_RE). Returns orf_id unchanged if it doesn't match the
    expected suffix pattern."""
    m = ORF_SUFFIX_RE.match(orf_id)
    return m.group(1) if m else orf_id


def strip_mate_tag(read_id: str) -> str:
    """Undo the _R1/_R2 suffix _merge_paired_with_mate_tags appends to
    disambiguate mate qseqids, recovering the original FASTQ read ID (as
    it appears, untagged, in the actual R1/R2 input files). Returns
    read_id unchanged if it doesn't end in _R1/_R2 -- e.g. single-end
    reads, which are never tagged since there's no mate to collide with."""
    if read_id.endswith("_R1") or read_id.endswith("_R2"):
        return read_id[:-3]
    return read_id


def run_orfm(fastq_path: str | Path, out_faa: str | Path) -> Path:
    """
    6-frame-translate *fastq_path* with ``orfm``, matching osmo_refdb's own
    pipeline/translate_and_hmmscan.sh usage. Used by `profile`'s
    DIAMOND+HMM cascade to translate just the subset of reads that need a
    second opinion (see extract_fastq_subset), not the whole input.
    """
    out_faa = Path(out_faa)
    fastq_path = Path(fastq_path)

    log.info("Running orfm → %s", out_faa)
    with open(out_faa, "w") as out_fh:
        result = subprocess.run(["orfm", str(fastq_path)], stdout=out_fh)
    if result.returncode != 0:
        raise RuntimeError(f"orfm failed (exit {result.returncode})")

    return out_faa


def extract_fastq_subset(
    fastq_paths: list[str | Path], wanted_ids: set[str], out_path: str | Path
) -> int:
    """
    Write only the FASTQ records whose read ID (first whitespace-delimited
    token of the header, matching how iter_fastq_ids_and_seqs works
    elsewhere in this project) is in *wanted_ids* to *out_path*, across
    one or more (optionally gzipped) input files -- e.g. both R1 and R2
    for paired-end input. Returns the number of records written.

    Used by `profile`'s DIAMOND+HMM cascade to pull just the reads that
    need a second HMM opinion out of what can be a very large input file,
    rather than re-processing everything.
    """
    import gzip

    out_path = Path(out_path)
    n_written = 0
    with open(out_path, "w") as out_fh:
        for fastq_path in fastq_paths:
            fastq_path = Path(fastq_path)
            opener = gzip.open if fastq_path.suffix == ".gz" else open
            with opener(fastq_path, "rt") as in_fh:
                while True:
                    header = in_fh.readline()
                    if not header:
                        break
                    seq = in_fh.readline()
                    plus = in_fh.readline()
                    qual = in_fh.readline()
                    read_id = header[1:].rstrip("\n").split()[0]
                    if read_id in wanted_ids:
                        out_fh.write(header)
                        out_fh.write(seq)
                        out_fh.write(plus)
                        out_fh.write(qual)
                        n_written += 1
    return n_written


# ---------------------------------------------------------------------------
# Prodigal ORF prediction
# ---------------------------------------------------------------------------

def run_prodigal(
    assembly: str | Path,
    proteins_out: str | Path,
    *,
    gff_out: str | Path | None = None,
) -> tuple[Path, Path | None]:
    """
    Run Prodigal in metagenomic mode (``-p meta``) on *assembly*.

    Returns (proteins_faa_path, gff_path_or_None).
    """
    proteins_out = Path(proteins_out)
    assembly = Path(assembly)

    cmd = [
        "prodigal",
        "-i", str(assembly),
        "-a", str(proteins_out),
        "-p", "meta",
        "-q",          # quiet
    ]
    if gff_out is not None:
        cmd += ["-f", "gff", "-o", str(gff_out)]
    else:
        cmd += ["-o", "/dev/null"]

    log.info("Running Prodigal → %s", proteins_out)
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise RuntimeError(f"prodigal failed (exit {result.returncode})")

    return proteins_out, (Path(gff_out) if gff_out else None)


# ---------------------------------------------------------------------------
# Alignment file parsing
# ---------------------------------------------------------------------------

def parse_diamond_output(
    aln_file: str | Path,
    *,
    min_identity: float = 0.0,
    min_query_cover: float = 0.0,
    min_subject_cover: float = 0.0,
    max_evalue: float = 1.0,
) -> Iterator[dict]:
    """
    Yield one dict per passing alignment row from a DIAMOND tabular file.

    Columns match DIAMOND_FIELDS:
      qseqid  sseqid  pident  qcovhsp  scovhsp  length  evalue  bitscore

    min_subject_cover filters on how much of the *target/reference*
    sequence the alignment covers, as opposed to min_query_cover (how much
    of the read/ORF aligns). This matters for multi-domain reference
    proteins: a short read can align at high identity and high query
    coverage to just one domain of a much longer fusion protein, passing
    every existing filter while only sampling a small fraction of that
    reference's actual length -- e.g. a read matching only the fused-in
    substrate-binding domain of an otherwise unrelated multi-domain
    transporter. min_subject_cover defaults to 0.0 (off) for backward
    compatibility; it's a deliberately blunt tool since a read is always
    much shorter than most full-length reference proteins, so a strict
    threshold will also discard many genuine short-read hits -- start low
    (e.g. 0.1-0.3) and check impact on recall before raising it.
    """
    with open(aln_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            row = {
                "qseqid":   parts[0],
                "sseqid":   parts[1],
                "pident":   float(parts[2]),
                "qcovhsp":  float(parts[3]),
                "scovhsp":  float(parts[4]),
                "length":   int(parts[5]),
                "evalue":   float(parts[6]),
                "bitscore": float(parts[7]),
            }
            if row["pident"] < min_identity * 100:
                continue
            if row["qcovhsp"] < min_query_cover * 100:
                continue
            if row["scovhsp"] < min_subject_cover * 100:
                continue
            if row["evalue"] > max_evalue:
                continue
            yield row


# ---------------------------------------------------------------------------
# Per-family DIAMOND cutoffs (osmo_refdb's 08b_calibrate_diamond_cutoffs.py)
# ---------------------------------------------------------------------------

def load_diamond_cutoffs(path: str | Path) -> dict[str, float]:
    """
    Load a per-family minimum-bitscore manifest produced by osmo_refdb's
    08b_calibrate_diamond_cutoffs.py (``<release>.diamond_cutoffs.tsv``).

    This exists because DIAMOND's --min_identity/--min_query_cover/--evalue
    filters are applied uniformly across every family, unlike HMM's
    per-family calibrated GA cutoff -- a flat identity threshold can't
    separate a true gene from a genuinely close paralog (e.g. betL vs
    betT/CaiT) whose sequence identity to it happens to exceed that flat
    threshold. Applying these cutoffs is optional (see
    quantifier.filter_by_family_cutoff): a family missing from this file
    (or the file not being supplied at all) is simply not gated any
    further than the existing global filters.

    Returns {family: cutoff_bitscore}.
    """
    cutoffs: dict[str, float] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cutoffs[row["family"]] = float(row["cutoff_bitscore"])
    return cutoffs


# ---------------------------------------------------------------------------
# Profile-mode DIAMOND+HMM cascade config (osmo_refdb's 11_compute_metrics.py)
# ---------------------------------------------------------------------------

def load_cascade_config(path: str | Path) -> dict[str, float]:
    """
    Load osmo_refdb's <release>.profile_cascade.tsv: families where
    DIAMOND's benchmark precision fell below a threshold, paired with
    HMM's best short-read bitscore cutoff for that family. Only rows with
    needs_cascade_check == "yes" are returned.

    Returns {family: hmm_short_read_threshold}.
    """
    cutoffs: dict[str, float] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["needs_cascade_check"] == "yes":
                cutoffs[row["family"]] = float(row["hmm_short_read_threshold"])
    return cutoffs


# ---------------------------------------------------------------------------
# Profile-mode scope exclusions (osmo_refdb's 08c_write_scope_manifest.py)
# ---------------------------------------------------------------------------

def load_excluded_families(path: str | Path) -> set[str]:
    """
    Load osmo_refdb's <release>.profile_excluded_families.txt: families
    marked scope: annotate_only in families.yaml (e.g. murB -- a
    near-universal housekeeping gene included so `annotate` can check
    genome-level co-occurrence, but carrying no osmoadaptation-specific
    signal at the read level). These families are still searched normally
    in `profile` (harmless, since they don't sequence-overlap with the
    osmoadaptation genes) -- this only controls what appears in the
    reported gene_counts.tsv.

    One family name per line; blank lines ignored.
    """
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


# ---------------------------------------------------------------------------
# Protein length extraction (from Prodigal FASTA headers)
# ---------------------------------------------------------------------------

def parse_protein_lengths(proteins_faa: str | Path) -> dict[str, int]:
    """
    Parse a Prodigal protein FASTA and return {protein_id: length_aa}.

    Prodigal header example:
      >NODE_1_length_50000_cov_5.0_1 # 1 # 1500 # 1 # ...
    Length in amino acids is derived from sequence characters.
    """
    lengths: dict[str, int] = {}
    current_id: str | None = None
    current_len = 0
    with open(proteins_faa) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_len
                current_id = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line.replace("*", ""))
    if current_id is not None:
        lengths[current_id] = current_len
    return lengths


@dataclass(frozen=True)
class OrfCoord:
    """One ORF's coordinates, parsed from a Prodigal-style protein header."""
    contig: str
    gene_index: int
    start: int
    end: int
    strand: int


# Prodigal's default protein FASTA header shape:
#   >NODE_1_length_50000_cov_5.0_7 # 1801 # 2400 # 1 # ID=1_7;partial=00;...
# i.e. "<contig>_<gene_index> # <start> # <end> # <strand> # ...". gene_index
# is Prodigal's own 1-based, per-contig sequential gene number, so
# abs(idx_a - idx_b) - 1 directly gives "number of intervening genes" with no
# extra bookkeeping. Matches the same trailing "_<index>" convention
# parse_protein_lengths() already relies on for the ID itself.
PRODIGAL_HEADER_RE = re.compile(
    r"^(?P<contig>.+)_(?P<gene_index>\d+)\s+#\s*(?P<start>\d+)\s*#\s*(?P<end>\d+)\s*#\s*(?P<strand>-?1)\s*#"
)


def parse_orf_coordinates(proteins_faa: str | Path) -> dict[str, OrfCoord]:
    """
    Parse a Prodigal protein FASTA and return {protein_id: OrfCoord}.

    Best-effort: a header that doesn't match Prodigal's own
    "<contig>_<gene_index> # start # end # strand # ..." shape is silently
    skipped rather than raising -- this is what lets the `annotate`
    --proteins path (pre-called, possibly non-Prodigal proteins) degrade
    gracefully to "systems detection skipped" instead of crashing (see
    systems.detect_systems callers in runners.run_annotate).

    Returns an empty dict if no header matched at all.
    """
    coords: dict[str, OrfCoord] = {}
    with open(proteins_faa) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].rstrip("\n")
            protein_id = header.split()[0]
            m = PRODIGAL_HEADER_RE.match(header)
            if not m:
                continue
            coords[protein_id] = OrfCoord(
                contig=m.group("contig"),
                gene_index=int(m.group("gene_index")),
                start=int(m.group("start")),
                end=int(m.group("end")),
                strand=int(m.group("strand")),
            )
    return coords


def count_reads(
    reads1: str | Path | None,
    reads2: str | Path | None,
    singles: str | Path | None,
) -> int:
    """
    Count total reads across input FASTQ files (handles .gz).

    Uses ``zcat | wc -l`` (gzipped) or ``wc -l`` (plain) via the shell —
    far faster than iterating in Python on the large files common on HPC.
    FASTQ has 4 lines per read, so read count = line count / 4.

    Falls back to pure-Python iteration if the shell command fails
    (e.g. ``zcat`` not available).
    """
    files: list[Path] = []
    if reads1 is not None:
        files.append(Path(reads1))
    if reads2 is not None:
        files.append(Path(reads2))
    if singles is not None:
        files.append(Path(singles))

    total = 0
    for path in files:
        total += _count_reads_file(path)
    return total


def _count_reads_file(path: Path) -> int:
    """Count reads in a single FASTQ file using the shell (fast path)."""
    is_gz = str(path).endswith(".gz")
    try:
        if is_gz:
            # zcat → wc -l; use pipefail so errors surface
            cmd = f"zcat {path} | wc -l"
        else:
            cmd = f"wc -l < {path}"
        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, check=True
        )
        line_count = int(result.stdout.strip().split()[0])
        return line_count // 4
    except Exception:
        log.debug("Shell read-count failed for %s, falling back to Python", path)
        return _count_reads_python(path)


def _count_reads_python(path: Path) -> int:
    """Pure-Python fallback read counter (slower, no shell dependency)."""
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    n = 0
    with opener(path, "rt") as fh:
        for i, _ in enumerate(fh):
            if i % 4 == 0:
                n += 1
    return n
