"""
catalog.py — catalog-based read assignment for `profile` (Cayman-style)

The hard step, deciding which family a sequence belongs to, is done once,
at full length, on a gene catalog (e.g. a prevalence-filtered GMGC subset,
or UHGG genes): every catalog gene is translated and searched with the
release's HMMs under their calibrated gathering cutoffs, exactly the
decision `annotate --method hmm` makes on Prodigal ORFs. Reads are then
assigned by near-exact nucleotide matching (bwa-mem3 mem by default) to those genes.

  build_catalog       genes.fna (+ optional proteins) -> catalog directory:
                      catalog.fna (only genes with an HMM hit), bwa index,
                      catalog.annotation.tsv, catalog.json
  map_reads           reads -> per-read catalog decisions:
                        call      primary alignment >= min_identity over
                                  >= min_aligned bp, on an annotated gene,
                                  inside its HMM envelope
                        drop      same alignment quality, but on a gene the
                                  full-length HMMs did NOT assign (a homolog
                                  below GA, or an HMM decoy), or outside the
                                  envelope: confidently not that family
                        fallback  unmapped or below the thresholds: the read
                                  goes on to the translated search (DIAMOND
                                  + cascade + report levels), unless
                                  profile runs catalog-only

Genes with only a sub-GA hit are kept in catalog.fna on purpose: a read from
a hard negative (same domain, different gene) then maps to its own gene and
is dropped, instead of being forced onto the nearest annotated family or
passed to translated search, where those reads are the main source of false
positives.

Catalog-bound by construction: above ~97% nucleotide identity only strains
whose genes are in the catalog are seen. Everything else falls back.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import logging
import re
import subprocess
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator

from osmotool.functions import count_hmm_models

log = logging.getLogger("osmotool")

CATALOG_FNA = "catalog.fna"
CATALOG_ANNOTATION = "catalog.annotation.tsv"
CATALOG_META = "catalog.json"

# Aligners with bwa's CLI ("<bin> index ref", "<bin> mem -t N ref r1 [r2]")
# and SAM output carrying NM. Index file suffixes differ, so the catalog
# records which one built it and profile maps with the same one.
# bwa-mem3 (Fulcrum Genomics; the engine behind the bwa-mem3-rs Rust crate)
# is the default: a drop-in for bwa mem / bwa-mem2, on bioconda, reading
# bwa-mem2 indexes. It needs AVX2 on x86-64 (or NEON / Apple M1+); on older
# CPUs use aligner="bwa". aligner_bin swaps in any executable with the same
# CLI and index format (e.g. a bwa-mem2 reimplementation).
DEFAULT_ALIGNER = "bwa-mem3"
ALIGNERS = {
    "bwa-mem3": {"bin": "bwa-mem3",
                 "index_suffixes": (".amb", ".ann", ".bwt.2bit.64", ".pac")},
    "bwa-mem2": {"bin": "bwa-mem2",
                 "index_suffixes": (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")},
    "bwa": {"bin": "bwa", "index_suffixes": (".amb", ".ann", ".bwt", ".pac", ".sa")},
}

ANNOTATED, HOMOLOG, DECOY = "annotated", "homolog", "decoy"

# Bacterial/archaeal code (NCBI table 11); only the amino acids matter here.
_BASES = "TCAG"
_AA11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: _AA11[16 * i + 4 * j + k]
               for i, a in enumerate(_BASES) for j, b in enumerate(_BASES)
               for k, c in enumerate(_BASES)}


def translate(nt: str) -> str:
    """Translate a CDS in frame 1 (table 11). Codons with ambiguous bases
    become X; a trailing partial codon and a final stop are dropped, an
    internal stop stays '*' (HMMER accepts it)."""
    nt = nt.upper().replace("U", "T")
    aa = "".join(CODON_TABLE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))
    return aa[:-1] if aa.endswith("*") else aa


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """(id, sequence) per record; id is the header up to the first whitespace,
    as bwa and HMMER read it. Gzip is detected from the suffix."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    name, chunks = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                parts = line[1:].split()
                name, chunks = (parts[0] if parts else ""), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def read_ga_cutoffs(hmm_db: str | Path) -> dict[str, tuple[float, float]]:
    """{model name: (GA sequence score, GA domain score)} from an HMMER3 file;
    models without a GA line are absent (they can never annotate a gene)."""
    cutoffs: dict[str, tuple[float, float]] = {}
    name = None
    with open(hmm_db) as fh:
        for line in fh:
            if line.startswith("NAME "):
                name = line.split()[1]
            elif line.startswith("GA ") and name is not None:
                vals = line.split()[1:3]
                cutoffs[name] = (float(vals[0]), float(vals[1].rstrip(";")))
    return cutoffs


def sha256_of(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Full-length annotation
# ---------------------------------------------------------------------------

def parse_domtblout(path: str | Path) -> dict[str, dict[str, dict]]:
    """{gene: {model: {"score", "evalue", "domains": [(env_from, env_to, score)]}}}
    from hmmsearch --domtblout (target = gene protein, query = model)."""
    hits: dict[str, dict[str, dict]] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            gene, model = f[0], f[3]
            h = hits.setdefault(gene, {}).setdefault(
                model, {"score": float(f[7]), "evalue": float(f[6]), "domains": []})
            h["domains"].append((int(f[19]), int(f[20]), float(f[13])))
    return hits


def assign_gene(hits: dict[str, dict], ga: dict[str, tuple[float, float]],
                decoys: set[str]) -> dict:
    """The full-length decision for one gene, as annotate --method hmm makes
    it: keep the models whose sequence score reaches GA1 with at least one
    domain reaching GA2 (= hmmsearch --cut_ga), then the best one by score
    (ties: lower E-value, then name, as quantifier.select_best_hits). A
    decoy model winning marks the gene a decoy; no model passing marks it a
    homolog (below GA). The envelope is the span of the winning model's
    domains that reach GA2 (nucleotide coordinates are derived later)."""
    passing = []
    for model, h in hits.items():
        if model not in ga:
            continue
        ga1, ga2 = ga[model]
        good = [d for d in h["domains"] if d[2] >= ga2]
        if h["score"] >= ga1 and good:
            passing.append((-h["score"], h["evalue"], model, good))
    best_any = min(hits.items(), key=lambda kv: (-kv[1]["score"], kv[1]["evalue"], kv[0]))
    if not passing:
        return {"label": "", "status": HOMOLOG, "score": "",
                "best_model": best_any[0], "best_score": best_any[1]["score"],
                "env_aa": None}
    _, _, model, good = min(passing)
    env = (min(d[0] for d in good), max(d[1] for d in good))
    return {"label": model, "status": DECOY if model in decoys else ANNOTATED,
            "score": hits[model]["score"], "best_model": best_any[0],
            "best_score": best_any[1]["score"], "env_aa": env}


def build_catalog(
    hmm_db: str | Path,
    genes: str | Path,
    out_dir: str | Path,
    *,
    proteins: str | Path | None = None,
    decoy_models: set[str] | None = None,
    release: str | Path | None = None,
    threads: int = 4,
    min_protein_length: int = 30,
    evalue: float = 1e-5,
    aligner: str = DEFAULT_ALIGNER,
    aligner_bin: str | None = None,
) -> Path:
    """
    Annotate a nucleotide gene catalog once with the release's HMMs and write
    a catalog directory `profile --catalog` maps reads to.

    genes:    nucleotide CDS FASTA (gz ok), one gene per record, in frame.
    proteins: optional protein FASTA with the same IDs (e.g. GMGC's .faa);
              otherwise genes are translated (table 11).
    decoy_models: HMM decoy names (osmo_refdb hmm_decoys, e.g. pepABC_decoy):
              a gene they win is kept as a decoy, never counted.
    evalue:   hmmsearch reporting threshold (per-sequence and per-domain, with
              -Z/--domZ fixed to the number of models as in run_hmmsearch).
              Genes reported but below GA are kept as homologs: see the
              module docstring for why.
    aligner, aligner_bin: which bwa-style aligner indexes the catalog
              (ALIGNERS); aligner_bin overrides the executable.
    """
    if aligner not in ALIGNERS:
        raise ValueError(f"Unknown aligner {aligner!r}: expected one of {', '.join(ALIGNERS)}")
    aligner_bin = aligner_bin or ALIGNERS[aligner]["bin"]
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    decoy_models = decoy_models or set()
    ga = read_ga_cutoffs(hmm_db)
    if not ga:
        raise ValueError(f"No GA cutoffs in {hmm_db}: cannot annotate a catalog without them")

    # --- proteins for hmmsearch ---
    faa = out_dir / "catalog.search.faa"
    n_genes = n_short = 0
    gene_len: dict[str, int] = {}
    seen: set[str] = set()
    source = proteins if proteins is not None else genes
    with open(faa, "w") as out:
        for gid, seq in iter_fasta(source):
            if gid in seen:
                raise ValueError(f"{source}: duplicate sequence ID {gid}")
            seen.add(gid)
            n_genes += 1
            aa = seq.rstrip("*") if proteins is not None else translate(seq)
            if len(aa) < min_protein_length:
                n_short += 1
                continue
            out.write(f">{gid}\n{aa}\n")
    seen.clear()
    log.info("Catalog: %d genes read, %d shorter than %d aa skipped", n_genes, n_short,
             min_protein_length)

    # --- full-length HMM search (no --cut_ga: sub-GA homologs are kept) ---
    domtbl = out_dir / "catalog.hmmsearch.domtblout"
    n_models = str(count_hmm_models(hmm_db))
    cmd = ["hmmsearch", "-Z", n_models, "--domZ", n_models, "-E", str(evalue),
           "--domE", str(evalue), "--domtblout", str(domtbl), "--noali",
           "--cpu", str(threads), str(hmm_db), str(faa)]
    log.info("Running hmmsearch on the catalog → %s", domtbl)
    if subprocess.run(cmd, stdout=subprocess.DEVNULL).returncode != 0:
        raise RuntimeError("hmmsearch failed on the catalog")
    hits = parse_domtblout(domtbl)
    decisions = {g: assign_gene(h, ga, decoy_models) for g, h in hits.items()}

    # --- nucleotide sequences of every gene with a hit ---
    fna = out_dir / CATALOG_FNA
    with open(fna, "w") as out:
        for gid, seq in iter_fasta(genes):
            if gid in decisions:
                gene_len[gid] = len(seq)
                out.write(f">{gid}\n{seq}\n")
    missing = set(decisions) - set(gene_len)
    if missing:
        raise ValueError(f"{len(missing)} protein IDs have no gene in {genes}, e.g. "
                         f"{sorted(missing)[:3]}")

    # --- annotation table ---
    by_label: dict[str, int] = {}
    with open(out_dir / CATALOG_ANNOTATION, "w") as out:
        out.write("gene\tlabel\tstatus\tgene_nt\tenv_nt_start\tenv_nt_end\tscore\t"
                  "best_model\tbest_score\n")
        for gid in sorted(decisions):
            d = decisions[gid]
            if d["env_aa"] is None:
                start = end = ""
            else:
                start, end = (d["env_aa"][0] - 1) * 3 + 1, min(d["env_aa"][1] * 3, gene_len[gid])
            key = d["label"] if d["status"] == ANNOTATED else d["status"]
            by_label[key] = by_label.get(key, 0) + 1
            out.write(f"{gid}\t{d['label']}\t{d['status']}\t{gene_len[gid]}\t{start}\t{end}\t"
                      f"{d['score']}\t{d['best_model']}\t{d['best_score']}\n")

    # --- bwa index ---
    log.info("Indexing %d catalog genes with %s", len(decisions), aligner_bin)
    if subprocess.run([aligner_bin, "index", str(fna)], stdout=subprocess.DEVNULL,
                      stderr=subprocess.DEVNULL).returncode != 0:
        raise RuntimeError(f"{aligner_bin} index failed")

    from osmotool import __version__
    meta = {
        "built": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "osmotool_version": __version__,
        "release": str(release) if release is not None else "",
        "hmm_db": str(hmm_db),
        "hmm_sha256": sha256_of(hmm_db),
        "genes": str(genes),
        "proteins": str(proteins) if proteins is not None else "",
        "n_genes_input": n_genes,
        "n_genes_kept": len(decisions),
        "min_protein_length": min_protein_length,
        "hmmsearch_evalue": evalue,
        "aligner": aligner,
        "aligner_bin": aligner_bin,
        "genes_per_label": dict(sorted(by_label.items())),
    }
    (out_dir / CATALOG_META).write_text(json.dumps(meta, indent=2) + "\n")
    faa.unlink(missing_ok=True)
    log.info("Catalog written: %s (%d annotated, %d homologs, %d decoys)", out_dir,
             sum(v for k, v in by_label.items() if k not in (HOMOLOG, DECOY)),
             by_label.get(HOMOLOG, 0), by_label.get(DECOY, 0))
    return out_dir


# ---------------------------------------------------------------------------
# Loading a catalog
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CatalogGene:
    label: str            # family for annotated genes, else ""
    status: str           # annotated | homolog | decoy
    length_nt: int
    env: tuple[int, int] | None   # 1-based inclusive nucleotide envelope


@dataclass
class Catalog:
    root: Path
    fna: Path
    genes: dict[str, CatalogGene]
    meta: dict = field(default_factory=dict)
    aligner_bin: str = ALIGNERS[DEFAULT_ALIGNER]["bin"]


def load_catalog(path: str | Path, *, hmm_db: str | Path | None = None,
                 aligner_bin: str | None = None) -> Catalog:
    """Load a build_catalog directory. With *hmm_db*, refuse a catalog that
    was annotated with different HMMs (another release): its family calls
    would not be this release's. Reads are mapped with the aligner that
    built the index (catalog.json), or *aligner_bin* (same CLI and index
    format, e.g. a reimplementation) if given."""
    root = Path(path)
    fna = root / CATALOG_FNA
    meta_path = root / CATALOG_META
    if not meta_path.exists():
        raise FileNotFoundError(f"{root} is not an osmotool catalog (no {CATALOG_META}); "
                                "build it with `osmotool build-catalog`")
    meta = json.loads(meta_path.read_text())
    aligner = meta.get("aligner", DEFAULT_ALIGNER)
    if aligner not in ALIGNERS:
        raise ValueError(f"{root}: built with unknown aligner {aligner!r}")
    missing = [p.name for p in [root / CATALOG_ANNOTATION, fna]
               + [Path(str(fna) + s) for s in ALIGNERS[aligner]["index_suffixes"]] if not p.exists()]
    if missing:
        raise FileNotFoundError(f"{root} is not a complete osmotool catalog (missing "
                                f"{', '.join(missing)}); rebuild it with `osmotool build-catalog`")
    if hmm_db is not None and meta.get("hmm_sha256") != sha256_of(hmm_db):
        raise ValueError(f"{root} was annotated with a different HMM database "
                         f"({meta.get('hmm_db')}) than {hmm_db}; rebuild it for this release")
    genes: dict[str, CatalogGene] = {}
    with open(root / CATALOG_ANNOTATION) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            env = ((int(f[idx["env_nt_start"]]), int(f[idx["env_nt_end"]]))
                   if f[idx["env_nt_start"]] else None)
            genes[f[idx["gene"]]] = CatalogGene(f[idx["label"]], f[idx["status"]],
                                                int(f[idx["gene_nt"]]), env)
    return Catalog(root=root, fna=fna, genes=genes, meta=meta,
                   aligner_bin=aligner_bin or meta.get("aligner_bin") or ALIGNERS[aligner]["bin"])


# ---------------------------------------------------------------------------
# Read mapping
# ---------------------------------------------------------------------------

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
MATE_SUFFIX_RE = re.compile(r"/[12]$")


def norm_read_id(read_id: str) -> str:
    """bwa-style aligners drop a trailing /1 or /2 from read names; FASTQ IDs are compared
    after the same normalisation."""
    return MATE_SUFFIX_RE.sub("", read_id)


def classify_alignment(fields: list[str], catalog: Catalog, *, min_identity: float,
                       min_aligned: int, min_envelope_overlap: float = 0.5
                       ) -> tuple[str, str, float, int]:
    """
    Decision for one primary SAM record: (decision, label, identity, aligned_bp)
    with decision call | drop | fallback (module docstring).

    identity = 1 - NM / alignment columns (M/=/X + I + D); aligned_bp = read
    bases in the alignment (M/=/X + I). Soft clips do not count against
    identity: a read hanging off a gene end aligns partially, and is judged
    on its aligned part as long as that is >= min_aligned bp.
    A call also needs >= min_envelope_overlap of the aligned reference span
    inside the gene's HMM envelope (reads on another domain of a multi-domain
    gene are not the family).
    """
    flag = int(fields[1])
    if flag & 0x4:
        return "fallback", "", 0.0, 0
    cols = aligned = ref_span = 0
    for n, op in CIGAR_RE.findall(fields[5]):
        n = int(n)
        if op in "M=X":
            cols += n; aligned += n; ref_span += n
        elif op == "I":
            cols += n; aligned += n
        elif op in "DN":
            cols += n; ref_span += n
    nm = next((int(t[5:]) for t in fields[11:] if t.startswith("NM:i:")), None)
    if nm is None or cols == 0:
        return "fallback", "", 0.0, aligned
    identity = 1.0 - nm / cols
    if identity < min_identity or aligned < min_aligned:
        return "fallback", "", identity, aligned
    gene = catalog.genes.get(fields[2])
    if gene is None or gene.status != ANNOTATED or gene.env is None:
        return "drop", gene.status if gene else "unknown_gene", identity, aligned
    start = int(fields[3])
    end = start + ref_span - 1
    overlap = max(0, min(end, gene.env[1]) - max(start, gene.env[0]) + 1)
    if overlap < min_envelope_overlap * ref_span:
        return "drop", "outside_envelope", identity, aligned
    return "call", gene.label, identity, aligned


@dataclass
class CatalogResult:
    counts: dict[str, float]                  # family -> weighted reads called
    lengths_aa: dict[str, float]              # family -> read-weighted mean gene length
    accepted: dict[int, set[str]]             # mate (0 single, 1, 2) -> normalised read IDs
    n_reads: int = 0
    n_called: int = 0
    n_dropped: int = 0
    n_fallback: int = 0
    dropped_by_reason: dict[str, int] = field(default_factory=dict)


def map_reads(
    catalog: Catalog,
    *,
    reads1: str | None = None,
    reads2: str | None = None,
    singles: str | None = None,
    threads: int = 4,
    min_identity: float = 0.97,
    min_aligned: int = 60,
    hits_out: str | Path | None = None,
) -> CatalogResult:
    """Run the catalog's aligner (bwa-mem3 by default) against the catalog and classify every read's primary
    alignment (classify_alignment). Paired mates are judged independently
    and weigh 0.5 each, as in count_hits. With *hits_out*, every call/drop
    is written as a TSV (read, mate, gene, decision, label, identity,
    aligned_bp) for auditing and benchmarking."""
    paired = reads1 is not None and reads2 is not None
    inputs = [reads1, reads2] if paired else [singles or reads1]
    weight = 0.5 if paired else 1.0
    cmd = [catalog.aligner_bin, "mem", "-v", "1", "-t", str(threads), str(catalog.fna),
           *map(str, inputs)]
    log.info("Mapping reads to the catalog with %s mem (min identity %.2f, >= %d bp)",
             catalog.aligner_bin, min_identity, min_aligned)
    res = CatalogResult(counts={}, lengths_aa={}, accepted={0: set(), 1: set(), 2: set()})
    len_sum: dict[str, float] = {}
    out_fh = open(hits_out, "w") if hits_out else None
    if out_fh:
        out_fh.write("read\tmate\tgene\tdecision\tlabel\tidentity\taligned_bp\n")
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    try:
        for line in proc.stdout:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            flag = int(fields[1])
            if flag & 0x900:          # secondary / supplementary
                continue
            res.n_reads += 1
            mate = 1 if flag & 0x40 else 2 if flag & 0x80 else 0
            decision, label, ident, aligned = classify_alignment(
                fields, catalog, min_identity=min_identity, min_aligned=min_aligned)
            if decision == "fallback":
                res.n_fallback += 1
                continue
            rid = norm_read_id(fields[0])
            res.accepted[mate].add(rid)
            if decision == "call":
                res.n_called += 1
                res.counts[label] = res.counts.get(label, 0.0) + weight
                len_sum[label] = len_sum.get(label, 0.0) + weight * catalog.genes[fields[2]].length_nt / 3
            else:
                res.n_dropped += 1
                res.dropped_by_reason[label] = res.dropped_by_reason.get(label, 0) + 1
            if out_fh:
                out_fh.write(f"{rid}\t{mate}\t{fields[2]}\t{decision}\t{label}\t"
                             f"{ident:.4f}\t{aligned}\n")
    finally:
        proc.stdout.close()
        rc = proc.wait()
        if out_fh:
            out_fh.close()
    if rc != 0:
        raise RuntimeError(f"{catalog.aligner_bin} mem failed (exit {rc})")
    res.lengths_aa = {f: len_sum[f] / c for f, c in res.counts.items() if c}
    log.info("Catalog: %d reads, %d called, %d dropped (%s), %d to fallback",
             res.n_reads, res.n_called, res.n_dropped,
             ", ".join(f"{k} {v}" for k, v in sorted(res.dropped_by_reason.items())) or "none",
             res.n_fallback)
    return res


def write_fallback_fastq(source: str | Path, accepted: set[str], dest: str | Path) -> int:
    """Copy every FASTQ record of *source* whose (normalised) read ID was not
    settled by the catalog to *dest*; returns the number written. Memory is
    proportional to the settled reads, which in a metagenome are few."""
    source = Path(source)
    opener = gzip.open if source.suffix == ".gz" else open
    n = 0
    with opener(source, "rt") as fin, open(dest, "w") as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            rec = header + fin.readline() + fin.readline() + fin.readline()
            if norm_read_id(header[1:].split()[0]) not in accepted:
                fout.write(rec)
                n += 1
    return n
