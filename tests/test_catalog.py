"""
test_catalog.py — catalog-based read assignment (catalog.py) and its outputs
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from osmotool.catalog import (
    ALIGNERS, ANNOTATED, DECOY, HOMOLOG, Catalog, CatalogGene, assign_gene,
    classify_alignment, load_catalog, norm_read_id, parse_domtblout, read_ga_cutoffs,
    translate, write_fallback_fastq,
)
from osmotool.functions import ReportUnit
from osmotool.output import write_gene_counts, write_unit_counts
from osmotool.quantifier import aggregate_report_units


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


# ---------------------------------------------------------------------------
# read classification
# ---------------------------------------------------------------------------

def _catalog() -> Catalog:
    return Catalog(root=Path("."), fna=Path("catalog.fna"), genes={
        "geneA": CatalogGene("ectA", ANNOTATED, 600, (31, 540)),
        "geneH": CatalogGene("", HOMOLOG, 600, None),
        "geneD": CatalogGene("pep_decoy", DECOY, 900, (1, 900)),
    })


def _sam(rname="geneA", pos=100, cigar="150M", nm=1, flag=0):
    return ["r1", str(flag), rname, str(pos), "60", cigar, "*", "0", "0", "*", "*", f"NM:i:{nm}"]


@pytest.mark.parametrize("fields,decision,label", [
    (_sam(), "call", "ectA"),
    (_sam(nm=6), "fallback", ""),                              # 96% identity
    (_sam(cigar="50M100S"), "fallback", ""),                   # only 50 bp aligned
    (_sam(cigar="80M70S", nm=0), "call", "ectA"),              # runs off the gene end
    (_sam(rname="geneH"), "drop", HOMOLOG),
    (_sam(rname="geneD"), "drop", DECOY),
    (_sam(pos=500), "drop", "outside_envelope"),               # 500-649 vs env ..540
    (_sam(pos=430), "call", "ectA"),                           # 430-579: 111 of 150 inside
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

def _catalog_dir(tmp_path: Path, aligner="bwa-mem3", hmm_sha="abc") -> Path:
    (tmp_path / "catalog.fna").write_text(">geneA\nACGT\n")
    for suf in ALIGNERS[aligner]["index_suffixes"]:
        (tmp_path / f"catalog.fna{suf}").write_text("")
    (tmp_path / "catalog.annotation.tsv").write_text(
        "gene\tlabel\tstatus\tgene_nt\tenv_nt_start\tenv_nt_end\tscore\tbest_model\tbest_score\n"
        "geneA\tectA\tannotated\t600\t31\t540\t120.0\tectA\t120.0\n"
        "geneH\t\thomolog\t600\t\t\t\tectA\t22.0\n")
    (tmp_path / "catalog.json").write_text(json.dumps({"aligner": aligner, "hmm_sha256": hmm_sha}))
    return tmp_path


def test_load_catalog(tmp_path: Path):
    cat = load_catalog(_catalog_dir(tmp_path))
    assert cat.aligner_bin == "bwa-mem3"
    assert cat.genes["geneA"] == CatalogGene("ectA", ANNOTATED, 600, (31, 540))
    assert cat.genes["geneH"].env is None
    assert load_catalog(tmp_path, aligner_bin="/opt/x/bwa-mem3").aligner_bin == "/opt/x/bwa-mem3"


def test_load_catalog_missing_index(tmp_path: Path):
    d = _catalog_dir(tmp_path)
    (d / "catalog.fna.bwt.2bit.64").unlink()
    with pytest.raises(FileNotFoundError, match="bwt.2bit.64"):
        load_catalog(d)


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
    sources = {"ectA": (2.0, 1.0), "opuAB": (1.0, 0.0), "opuBB": (0.5, 1.5), "opuAC": (1.5, 2.5)}
    rows = {r["unit"]: r for r in aggregate_report_units(UNITS, counts, counts, sources=sources)}
    assert (rows["ectA"]["catalog_count"], rows["ectA"]["translated_count"]) == (2.0, 1.0)
    assert (rows["opu_permease"]["catalog_count"], rows["opu_permease"]["translated_count"]) == (1.5, 1.5)
    ac = rows["opuAC"]
    assert ac["status"] == "catalog_only" and ac["catalog_count"] == 1.5
    assert ac["raw_count"] is None and ac["translated_count"] is None


def test_unit_counts_file_with_sources(tmp_path: Path):
    counts = {"ectA": 3.0, "opuAC": 4.0}
    sources = {"ectA": (2.0, 1.0), "opuAC": (1.5, 2.5)}
    rows = aggregate_report_units(UNITS, counts, counts, sources=sources)
    lines = write_unit_counts(tmp_path / "s", rows, total_reads=10, filtered_reads=5,
                              report_levels="r.tsv", with_sources=True).read_text().splitlines()
    assert lines[3].endswith("RPM\tcatalog_count\ttranslated_count")
    assert lines[4] == "ectA\tgene\tectA\treported\t3.0\t3\t2.0\t1.0"
    assert lines[6] == "opuAC\tnone\topuAC\tcatalog_only\tNA\tNA\t1.5\tNA"


def test_gene_counts_with_sources(tmp_path: Path):
    path = write_gene_counts(tmp_path / "s", {"ectA": 3.0}, {"ectA": 30.0}, total_reads=10,
                             filtered_reads=3, count_sources={"ectA": (2.0, 1.0)},
                             copies_per_genome={"ectA": 0.5})
    lines = path.read_text().splitlines()
    assert lines[2] == "gene\traw_count\tRPM\tcatalog_count\ttranslated_count\tcopies_per_genome"
    assert lines[3] == "ectA\t3.0\t30\t2.0\t1.0\t0.5"
