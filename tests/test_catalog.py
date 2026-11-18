"""
test_catalog.py — catalog-based read assignment (catalog.py) and its outputs
"""

from __future__ import annotations

import dataclasses
import json
import random
import shutil
from pathlib import Path

import pytest

from osmotool.catalog import (
    ALIGNERS, ANNOTATED, BACKGROUND, DECOY, HOMOLOG, INDEX_META, Catalog, CatalogGene,
    alternative_placements, assign_gene, classify_alignment, find_neighbours, load_catalog,
    neighbour_windows, norm_read_id, parse_domtblout, read_ga_cutoffs, select_neighbours,
    settle_records, translate, write_fallback_fastq,
)
from osmotool.functions import ReportUnit
from osmotool.output import write_gene_counts, write_unit_counts
from osmotool.quantifier import aggregate_report_units
from osmotool.utils import default_catalog, resolve_refdb, select_catalog


# ---------------------------------------------------------------------------
# translation, HMM file, domtblout
# ---------------------------------------------------------------------------

def test_translate_table11_drops_final_stop_and_partial_codon():
    assert translate("ATGGCTTGGTAA") == "MAW"
    assert translate("ATGGCTTGGTAAGC") == "MAW*"[:-1]   # partial codon dropped, then stop
    assert translate("ATGTAAGCT") == "M*A"              # internal stop kept
    assert translate("ATGNNNGCT") == "MXA"


def test_read_ga_cutoffs(tmp_path: Path):
    hmm = tmp_path / "x.hmm"
    hmm.write_text("HMMER3/f\nNAME  ectA\nGA    25.00 20.00;\n//\n"
                   "HMMER3/f\nNAME  noga\n//\n"
                   "HMMER3/f\nNAME  pep_decoy\nGA    30.5 30.5;\n//\n")
    assert read_ga_cutoffs(hmm) == {"ectA": (25.0, 20.0), "pep_decoy": (30.5, 30.5)}


DOMTBL = (
    "# comment\n"
    # target tacc tlen query qacc qlen E score bias # of c-E i-E dscore dbias hf ht af at ef et acc desc
    "g1 - 200 ectA - 180 1e-30 120.0 0.1 1 1 1e-31 1e-31 118.0 0.1 1 180 5 190 3 195 0.9 -\n"
    "g1 - 200 ectB - 300 1e-10 40.0 0.1 1 1 1e-11 1e-11 39.0 0.1 1 90 10 100 8 104 0.9 -\n"
    "g2 - 150 ectA - 180 1e-6 22.0 0.1 1 1 1e-7 1e-7 21.0 0.1 1 90 10 100 8 104 0.9 -\n"
    "g3 - 300 pep_decoy - 250 1e-40 200.0 0.1 1 1 1e-41 1e-41 199.0 0.1 1 250 2 260 1 262 0.9 -\n"
    "g3 - 300 ectA - 180 1e-20 90.0 0.1 1 1 1e-21 1e-21 89.0 0.1 1 180 2 182 1 185 0.9 -\n"
    # g4: ectA split into two domains, each below GA 25 but the sequence above it
    "g4 - 400 ectA - 180 1e-12 40.0 0.1 1 2 1e-8 1e-8 21.0 0.1 1 90 5 95 5 98 0.9 -\n"
    "g4 - 400 ectA - 180 1e-12 40.0 0.1 2 2 1e-8 1e-8 19.0 0.1 91 180 200 300 199 305 0.9 -\n"
    # g5: fused ORF, ectA and ectB both pass in different regions
    "g5 - 600 ectA - 180 1e-30 100.0 0.1 1 1 1e-31 1e-31 99.0 0.1 1 180 1 190 1 200 0.9 -\n"
    "g5 - 600 ectB - 300 1e-25 80.0 0.1 1 1 1e-26 1e-26 79.0 0.1 1 300 260 580 255 590 0.9 -\n"
)
GA = {"ectA": (25.0, 20.0), "ectB": (30.0, 30.0), "pep_decoy": (30.0, 30.0)}


def test_parse_and_assign(tmp_path: Path):
    p = tmp_path / "d.tbl"
    p.write_text(DOMTBL)
    hits = parse_domtblout(p)
    a1 = assign_gene(hits["g1"], GA, {"pep_decoy"})
    assert (a1["label"], a1["status"], a1["env_aa"]) == ("ectA", ANNOTATED, (3, 195))
    a2 = assign_gene(hits["g2"], GA, {"pep_decoy"})          # 22 < GA 25
    assert a2["status"] == HOMOLOG and a2["label"] == "" and a2["best_model"] == "ectA"
    a3 = assign_gene(hits["g3"], GA, {"pep_decoy"})          # decoy wins the gene
    assert (a3["label"], a3["status"]) == ("pep_decoy", DECOY)
    a4 = assign_gene(hits["g4"], GA, {"pep_decoy"})          # as hmmsearch --cut_ga's tblout
    assert (a4["label"], a4["status"], a4["env_aa"]) == ("ectA", ANNOTATED, (5, 305))
    a5 = assign_gene(hits["g5"], GA, {"pep_decoy"})          # fused: ectB kept as other region
    assert a5["label"] == "ectA" and a5["others"] == [("ectB", (255, 590))]


# ---------------------------------------------------------------------------
# read classification
# ---------------------------------------------------------------------------

def _catalog() -> Catalog:
    return Catalog(root=Path("."), fna=Path("catalog.fna"), genes={
        "geneA": CatalogGene("ectA", ANNOTATED, 600, (31, 540)),
        "geneA2": CatalogGene("ectA", ANNOTATED, 600, (31, 540)),
        "geneB": CatalogGene("ectB", ANNOTATED, 900, (1, 900)),
        "geneH": CatalogGene("", HOMOLOG, 600, None),
        "geneD": CatalogGene("pep_decoy", DECOY, 900, (1, 900)),
        "geneF": CatalogGene("ectA", ANNOTATED, 1800, (1, 600), (("ectB", 765, 1770),)),
    })


def _sam(rname="geneA", pos=100, cigar="150M", nm=1, flag=0, xa=None, qname="r1"):
    f = [qname, str(flag), rname, str(pos), "60", cigar, "*", "0", "0", "*", "*", f"NM:i:{nm}"]
    if xa:
        f.append("XA:Z:" + "".join(f"{g},+{p},{c},{n};" for g, p, c, n in xa))
    return f


@pytest.mark.parametrize("fields,decision,label", [
    (_sam(), "call", "ectA"),
    (_sam(nm=6), "fallback", ""),                              # 96% identity
    (_sam(cigar="50M100S"), "fallback", ""),                   # only 50 bp aligned
    (_sam(cigar="80M70S", nm=0), "call", "ectA"),              # runs off the gene end
    (_sam(rname="geneH"), "drop", HOMOLOG),
    (_sam(rname="geneD"), "drop", DECOY),
    (_sam(rname="geneX"), "drop", BACKGROUND),                 # indexed, no HMM hit
    (_sam(pos=500), "call", "ectA"),                           # past the envelope: whole gene counts
    (_sam(rname="geneF", pos=900), "call", "ectB"),            # fused gene: in ectB's region
    (_sam(rname="geneF", pos=100), "call", "ectA"),
    (_sam(rname="geneF", pos=620), "call", "ectA"),            # between regions: the gene's family
    (_sam(flag=4, rname="*"), "fallback", ""),
])
def test_classify_alignment(fields, decision, label):
    d, lab, _, _ = classify_alignment(fields, _catalog(), min_identity=0.97, min_aligned=60)
    assert (d, lab) == (decision, label)


def test_identity_counts_indels_as_columns():
    _, _, ident, aligned = classify_alignment(_sam(cigar="70M2I78M", nm=2), _catalog(),
                                              min_identity=0.97, min_aligned=60)
    assert aligned == 150 and ident == pytest.approx(1 - 2 / 150)


def test_norm_read_id():
    assert norm_read_id("abc/1") == "abc" and norm_read_id("abc_1") == "abc_1"


def test_write_fallback_fastq_skips_settled_reads(tmp_path: Path):
    fq = tmp_path / "r.fastq"
    fq.write_text("@a/1 x\nACGT\n+\nIIII\n@b/1\nACGT\n+\nIIII\n@c/1\nACGT\n+\nIIII\n")
    out = tmp_path / "o.fastq"
    assert write_fallback_fastq(fq, {"b"}, out) == 2
    assert [l for l in out.read_text().splitlines() if l.startswith("@")] == ["@a/1 x", "@c/1"]


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------

def _catalog_dir(tmp_path: Path, aligner="bwa-mem3", hmm_sha="abc", index="all",
                 genes_sha="g1") -> Path:
    (tmp_path / "catalog.fna").write_text(">geneA\nACGT\n")
    if index != "none":
        for suf in ALIGNERS[aligner]["index_suffixes"]:
            (tmp_path / f"catalog.fna{suf}").write_text("")
    (tmp_path / "catalog.annotation.tsv").write_text(
        "gene\tlabel\tstatus\tgene_nt\tenv_nt_start\tenv_nt_end\tscore\tbest_model\tbest_score\tother_envelopes\n"
        "geneA\tectA\tannotated\t600\t31\t540\t120.0\tectA\t120.0\tectB:301-540\n"
        "geneH\t\thomolog\t600\t\t\t\tectA\t22.0\t\n")
    (tmp_path / "catalog.json").write_text(json.dumps(
        {"aligner": aligner, "hmm_sha256": hmm_sha, "index": index, "genes_sha256": genes_sha,
         "genes": "genes.fna.gz"}))
    return tmp_path


def _index_dir(tmp_path: Path, genes_sha="g1", aligner="bwa-mem3") -> Path:
    tmp_path.mkdir(parents=True, exist_ok=True)
    for suf in ALIGNERS[aligner]["index_suffixes"]:
        (tmp_path / f"catalog.fna{suf}").write_text("")
    (tmp_path / INDEX_META).write_text(json.dumps(
        {"index": "all", "aligner": aligner, "genes_sha256": genes_sha, "genes": "genes.fna.gz"}))
    return tmp_path


def test_load_catalog(tmp_path: Path):
    cat = load_catalog(_catalog_dir(tmp_path))
    assert cat.aligner_bin == "bwa-mem3"
    assert cat.genes["geneA"] == CatalogGene("ectA", ANNOTATED, 600, (31, 540), (("ectB", 301, 540),))
    assert cat.genes["geneH"].env is None
    assert load_catalog(tmp_path, aligner_bin="/opt/x/bwa-mem3").aligner_bin == "/opt/x/bwa-mem3"


def test_load_catalog_needs_no_fasta(tmp_path: Path):
    d = _catalog_dir(tmp_path)
    (d / "catalog.fna").unlink()          # index="all": only the index is in the catalog dir
    assert load_catalog(d).fna == d / "catalog.fna"


def test_load_catalog_missing_index(tmp_path: Path):
    d = _catalog_dir(tmp_path)
    (d / "catalog.fna.bwt.2bit.64").unlink()
    with pytest.raises(FileNotFoundError, match="bwt.2bit.64"):
        load_catalog(d)


def test_annotation_only_catalog_needs_an_index(tmp_path: Path):
    (tmp_path / "cat").mkdir()
    d = _catalog_dir(tmp_path / "cat", index="none")
    with pytest.raises(ValueError, match="--catalog_index"):
        load_catalog(d)
    cat = load_catalog(d, index_dir=_index_dir(tmp_path / "idx"))
    assert cat.fna == tmp_path / "idx" / "catalog.fna" and cat.index_kind == "all"


def test_external_index_must_come_from_the_same_gene_file(tmp_path: Path):
    (tmp_path / "cat").mkdir()
    d = _catalog_dir(tmp_path / "cat", index="hits")
    assert load_catalog(d).index_kind == "hits"
    with pytest.raises(ValueError, match="different gene file"):
        load_catalog(d, index_dir=_index_dir(tmp_path / "idx", genes_sha="other"))
    with pytest.raises(FileNotFoundError, match="catalog_index.json"):
        load_catalog(d, index_dir=tmp_path / "cat")


def test_load_catalog_other_release(tmp_path: Path):
    hmm = tmp_path / "osmo_refdb.hmm"
    hmm.write_text("NAME x\n")
    cdir = tmp_path / "cat"
    cdir.mkdir()
    with pytest.raises(ValueError, match="different HMM"):
        load_catalog(_catalog_dir(cdir), hmm_db=hmm)


# ---------------------------------------------------------------------------
# outputs with count sources
# ---------------------------------------------------------------------------

UNITS = [ReportUnit("ectA", "gene", ("ectA",)),
         ReportUnit("opu_permease", "group", ("opuAB", "opuBB")),
         ReportUnit("opuAC", "none", ("opuAC",))]


def test_units_with_sources():
    counts = {"ectA": 3.0, "opuAB": 1.0, "opuBB": 2.0, "opuAC": 4.0}
    sources = {"ectA": (2.0, 0.5, 1.0), "opuAB": (1.0, 0.0, 0.0), "opuBB": (0.5, 0.25, 1.5),
               "opuAC": (1.5, 0.0, 2.5)}
    rows = {r["unit"]: r for r in aggregate_report_units(UNITS, counts, counts, sources=sources)}
    assert (rows["ectA"]["catalog_count"], rows["ectA"]["catalog_ambiguous"],
            rows["ectA"]["translated_count"]) == (2.0, 0.5, 1.0)
    assert (rows["opu_permease"]["catalog_count"], rows["opu_permease"]["catalog_ambiguous"],
            rows["opu_permease"]["translated_count"]) == (1.5, 0.25, 1.5)
    ac = rows["opuAC"]
    assert ac["status"] == "catalog_only" and ac["catalog_count"] == 1.5
    assert ac["raw_count"] is None and ac["translated_count"] is None


def test_unit_counts_file_with_sources(tmp_path: Path):
    counts = {"ectA": 3.0, "opuAC": 4.0}
    sources = {"ectA": (2.0, 0.5, 1.0), "opuAC": (1.5, 0.0, 2.5)}
    rows = aggregate_report_units(UNITS, counts, counts, sources=sources)
    lines = write_unit_counts(tmp_path / "s", rows, total_reads=10, filtered_reads=5,
                              report_levels="r.tsv", with_sources=True).read_text().splitlines()
    assert lines[3].endswith("RPM\tcatalog_count\tcatalog_ambiguous\ttranslated_count")
    assert lines[4] == "ectA\tgene\tectA\treported\t3.0\t3\t2.0\t0.50\t1.0"
    assert lines[6] == "opuAC\tnone\topuAC\tcatalog_only\tNA\tNA\t1.5\t0.00\tNA"


def test_gene_counts_with_sources(tmp_path: Path):
    path = write_gene_counts(tmp_path / "s", {"ectA": 3.0}, {"ectA": 30.0}, total_reads=10,
                             filtered_reads=3, count_sources={"ectA": (2.0, 0.25, 1.0)},
                             copies_per_genome={"ectA": 0.5})
    lines = path.read_text().splitlines()
    assert lines[2] == ("gene\traw_count\tRPM\tcatalog_count\tcatalog_ambiguous\ttranslated_count"
                        "\tcopies_per_genome")
    assert lines[3] == "ectA\t3.0\t30\t2.0\t0.25\t1.0\t0.5"


# ---------------------------------------------------------------------------
# settle_records: unique, ambiguous, fallback
# ---------------------------------------------------------------------------

def _line(**kw) -> str:
    return "\t".join(_sam(**kw)) + "\n"


def _settle(lines, paired=False):
    return settle_records(["@HD\tVN:1.6\n", *lines], _catalog(), paired=paired,
                          min_identity=0.97, min_aligned=45)


def test_alternative_placements_parses_xa():
    f = _sam(xa=[("geneB", 10, "150M", 0), ("geneH", 5, "100M50S", 2)])
    assert alternative_placements(f) == [("geneB", 10, "150M", 0), ("geneH", 5, "100M50S", 2)]
    assert alternative_placements(_sam()) == []


def test_settle_unique_call_and_same_family_alternative():
    # an alternative on another ectA gene does not make the read ambiguous
    res = _settle([_line(xa=[("geneA2", 100, "150M", 1)])])
    assert res.counts == {"ectA": 1.0} and res.n_ambiguous == 0 and res.ambiguous == {}


def test_settle_ambiguous_between_families_is_split_and_kept_out_of_counts():
    res = _settle([_line(xa=[("geneB", 10, "150M", 0)])])
    assert res.counts == {} and res.n_ambiguous == 1
    assert res.ambiguous == {"ectA": 0.5, "ectB": 0.5}
    assert res.accepted[0] == {"r1"}                    # settled: no fallback


def test_settle_ambiguous_call_vs_background_keeps_only_the_family_share():
    res = _settle([_line(xa=[("geneX", 10, "150M", 0)])])
    assert res.ambiguous == {"ectA": 0.5} and res.counts == {}


def test_settle_ignores_alternatives_below_thresholds():
    res = _settle([_line(xa=[("geneB", 10, "150M", 9)])])   # 94% identity
    assert res.counts == {"ectA": 1.0} and res.n_ambiguous == 0


def test_settle_fallback_secondary_and_mates():
    lines = [
        _line(qname="p/1", flag=0x41),                      # mate 1, call
        _line(qname="p/2", flag=0x81, nm=10),               # mate 2, 93%: fallback
        _line(qname="p/1", flag=0x141, rname="geneB"),      # secondary: skipped
    ]
    res = _settle(lines, paired=True)
    assert res.n_reads == 2 and res.n_fallback == 1
    assert res.counts == {"ectA": 0.5}
    assert res.accepted[1] == {"p"} and res.accepted[2] == set()


# ---------------------------------------------------------------------------
# catalogs shipped with a release
# ---------------------------------------------------------------------------

def _release_with_catalog(tmp_path: Path) -> Path:
    rel = tmp_path / "v99"
    rel.mkdir()
    (rel / "osmo_refdb.dmnd").write_text("")
    cdir = rel / "catalogs" / "gmgc_gut"
    cdir.mkdir(parents=True)
    _catalog_dir(cdir)
    (rel / "catalogs" / "not_a_catalog").mkdir()
    return rel


def test_resolve_refdb_finds_release_catalogs(tmp_path: Path):
    rel = _release_with_catalog(tmp_path)
    refdb = resolve_refdb(str(rel), mode="profile")
    assert refdb.catalogs == {"gmgc_gut": rel / "catalogs" / "gmgc_gut"}
    assert resolve_refdb(str(rel), mode="annotate").catalogs == {}


def test_select_catalog_by_name_or_path(tmp_path: Path):
    rel = _release_with_catalog(tmp_path)
    refdb = resolve_refdb(str(rel), mode="profile")
    assert select_catalog(refdb, "gmgc_gut") == rel / "catalogs" / "gmgc_gut"
    own = tmp_path / "own"
    own.mkdir()
    _catalog_dir(own)
    assert select_catalog(refdb, str(own)) == own
    with pytest.raises(ValueError, match="gmgc_gut"):
        select_catalog(refdb, "marine")


def test_default_catalog_is_the_one_the_release_ships(tmp_path: Path):
    rel = _release_with_catalog(tmp_path)
    refdb = resolve_refdb(str(rel), mode="profile")
    shipped = rel / "catalogs" / "gmgc_gut"
    assert default_catalog(refdb, None) == shipped
    assert default_catalog(refdb, "gmgc_gut") == shipped
    assert default_catalog(refdb, None, no_catalog=True) is None
    with pytest.raises(ValueError, match="exclude each other"):
        default_catalog(refdb, "gmgc_gut", no_catalog=True)


def test_default_catalog_none_shipped_or_several(tmp_path: Path):
    rel = _release_with_catalog(tmp_path)
    refdb = resolve_refdb(str(rel), mode="profile")
    assert default_catalog(dataclasses.replace(refdb, catalogs={}), None) is None
    several = dataclasses.replace(refdb, catalogs={"gmgc_gut": rel / "catalogs" / "gmgc_gut",
                                                   "marine": tmp_path / "marine"})
    with pytest.raises(ValueError, match="several gene catalogs"):
        default_catalog(several, None)
    assert default_catalog(several, "gmgc_gut") == rel / "catalogs" / "gmgc_gut"


# ---------------------------------------------------------------------------
# index="neighbours": the genes without a hit close to a hit gene
# ---------------------------------------------------------------------------

def test_neighbour_windows_cover_every_gene_without_a_hit(tmp_path: Path):
    fna = tmp_path / "g.fna"
    fna.write_text(">hit\n" + "A" * 400 + "\n>long\n" + "C" * 320 + "\n>short\n" + "G" * 90 + "\n")
    wins = list(neighbour_windows(fna, {"hit"}, window=150, step=50))
    assert [n for n, _ in wins] == ["long|0", "long|50", "long|100", "long|150", "long|170", "short|0"]
    assert all(len(w) == 150 for n, w in wins if n.startswith("long"))
    assert wins[-1][1] == "G" * 90


def test_select_neighbours_identity_length_and_alternatives():
    sam = ["@HD\tVN:1.6\n"] + ["\t".join(f) + "\n" for f in (
        _sam(qname="n1|0", nm=10),                          # 150M, 93.3%: a neighbour
        _sam(qname="n1|50", nm=2),                          # better window of the same gene
        _sam(qname="far|0", nm=20),                         # 86.7%: not
        _sam(qname="short|0", cigar="40M110S", nm=0),       # 40 bp aligned: not
        _sam(qname="alt|0", nm=30, xa=[("geneB", 10, "150M", 3)]),   # only its XA placement passes
        _sam(qname="unmapped|0", flag=4, rname="*"),
    )]
    best = select_neighbours(sam, min_identity=0.90, min_aligned=45)
    assert set(best) == {"n1", "alt"}
    assert best["n1"] == pytest.approx(1 - 2 / 150) and best["alt"] == pytest.approx(0.98)


@pytest.mark.skipif(shutil.which("bwa-mem3") is None, reason="needs bwa-mem3")
def test_find_neighbours_with_bwa_mem3(tmp_path: Path):
    rng = random.Random(5)
    hit = "".join(rng.choice("ACGT") for _ in range(1200))
    near = "".join(c if rng.random() > 0.05 else {"A": "C", "C": "G", "G": "T", "T": "A"}[c] for c in hit)
    other = "".join(rng.choice("ACGT") for _ in range(1200))
    genes = tmp_path / "genes.fna"
    genes.write_text(f">hit\n{hit}\n>near\n{near}\n>other\n{other}\n")
    hits_fna = tmp_path / "hits.fna"
    hits_fna.write_text(f">hit\n{hit}\n")
    work = tmp_path / "work"
    best = find_neighbours(genes, hits_fna, {"hit"}, work, aligner_bin="bwa-mem3", threads=1)
    assert set(best) == {"near"} and best["near"] >= 0.90
    assert not work.exists()
