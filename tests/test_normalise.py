"""
test_normalise.py — unit tests for osmotool.normalise (profile mode's
copies-per-genome normalisation) and the pieces of output/utils/quantifier
it relies on.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from osmotool.normalise import (
    median_read_length,
    reference_lengths,
    load_marker_names,
    effective_length_bp,
    compute_rpk,
    estimate_genome_coverage,
    copies_per_genome,
)
from osmotool.output import write_gene_counts, write_marker_report, build_aln_stats
from osmotool.quantifier import count_hits
from osmotool.utils import resolve_refdb


def _hit(qid, sid):
    return {"qseqid": qid, "sseqid": sid, "bitscore": 100.0, "evalue": 1e-10}


def _fastq(lengths):
    return "".join(f"@r{i}\n{'A' * n}\n+\n{'I' * n}\n" for i, n in enumerate(lengths))


# ---------------------------------------------------------------------------
# median_read_length
# ---------------------------------------------------------------------------

class TestMedianReadLength:
    def test_plain_fastq(self, tmp_path: Path):
        p = tmp_path / "r.fastq"
        p.write_text(_fastq([100, 150, 150, 151, 90]))
        assert median_read_length([p]) == 150

    def test_gzipped_fastq(self, tmp_path: Path):
        p = tmp_path / "r.fastq.gz"
        with gzip.open(p, "wt") as fh:
            fh.write(_fastq([250, 250, 249]))
        assert median_read_length([p]) == 250

    def test_skips_none_and_samples_first_n(self, tmp_path: Path):
        p = tmp_path / "r.fastq"
        p.write_text(_fastq([100, 100, 100, 300, 300, 300, 300]))
        assert median_read_length([None, p], n_reads=3) == 100

    def test_no_reads_raises(self, tmp_path: Path):
        p = tmp_path / "empty.fastq"
        p.write_text("")
        with pytest.raises(ValueError):
            median_read_length([p])


# ---------------------------------------------------------------------------
# reference_lengths / load_marker_names
# ---------------------------------------------------------------------------

class TestReferenceLengths:
    def test_hit_weighted_mean_of_subject_lengths(self):
        hits = {
            "r1": _hit("r1", "leuS|A|x"),
            "r2": _hit("r2", "leuS|A|x"),
            "r3": _hit("r3", "leuS|B|y"),
            "r4": _hit("r4", "ffh|C|z"),
        }
        lengths = {"leuS|A|x": 800, "leuS|B|y": 900, "ffh|C|z": 450}
        out = reference_lengths(hits, lengths)
        assert out["leuS"] == pytest.approx((800 + 800 + 900) / 3)
        assert out["ffh"] == 450

    def test_subject_missing_from_fasta_ignored(self):
        hits = {"r1": _hit("r1", "ffh|C|z")}
        assert reference_lengths(hits, {}) == {}


class TestLoadMarkerNames:
    def test_reads_first_column(self, tmp_path: Path):
        p = tmp_path / "m.tsv"
        p.write_text("marker\tn_train_refs\nleuS\t10\nffh\t5\n")
        assert load_marker_names(p) == ["leuS", "ffh"]

    def test_wrong_header_raises(self, tmp_path: Path):
        p = tmp_path / "m.tsv"
        p.write_text("family\tn\nleuS\t1\n")
        with pytest.raises(ValueError):
            load_marker_names(p)


# ---------------------------------------------------------------------------
# effective_length_bp / compute_rpk
# ---------------------------------------------------------------------------

class TestEffectiveLength:
    def test_formula(self):
        # 3*300 - (2*0.8 - 1)*150 + 1 = 900 - 90 + 1
        assert effective_length_bp(300, 150, 0.8) == pytest.approx(811.0)

    def test_matches_brute_force_window_count(self):
        """L_eff must equal the number of read start positions leaving
        >= c*R bases of the read inside the gene -- including a gene
        shorter than the read."""
        for gene_bp, R, c in [(900, 150, 0.8), (210, 150, 0.8), (120, 150, 0.8), (600, 100, 0.5)]:
            need = c * R
            n = sum(1 for s in range(-R, gene_bp + 1)
                    if min(s + R, gene_bp) - max(s, 0) >= need)
            assert effective_length_bp(gene_bp / 3, R, c) == pytest.approx(n, abs=1)

    def test_gene_too_short_for_any_read(self):
        assert effective_length_bp(30, 150, 0.8) == 0.0   # 90 bp < 120 bp needed


class TestComputeRpk:
    def test_basic(self):
        rpk = compute_rpk({"a": 81.1}, {"a": 300}, 150, 0.8)
        assert rpk["a"] == pytest.approx(100.0)   # 81.1 reads / 0.811 kb

    def test_zero_count_needs_no_length(self):
        assert compute_rpk({"a": 0.0}, {}, 150, 0.8) == {"a": 0.0}

    def test_reads_but_no_length_is_none(self):
        assert compute_rpk({"a": 3.0}, {}, 150, 0.8) == {"a": None}

    def test_short_gene_corrected_relative_to_long(self):
        """Same per-base coverage on a short and a long gene must give the
        same RPK once effective length is used (the reason it exists)."""
        R, c = 150, 0.8
        short_eff, long_eff = effective_length_bp(70, R, c), effective_length_bp(860, R, c)
        counts = {"short": short_eff / 100, "long": long_eff / 100}
        rpk = compute_rpk(counts, {"short": 70, "long": 860}, R, c)
        assert rpk["short"] == pytest.approx(rpk["long"])


# ---------------------------------------------------------------------------
# estimate_genome_coverage / copies_per_genome
# ---------------------------------------------------------------------------

def _markers(values):
    names = [f"m{i}" for i in range(len(values))]
    rpk = dict(zip(names, values))
    counts = {m: (200.0 if v else 0.0) for m, v in rpk.items()}
    return rpk, counts


class TestEstimateGenomeCoverage:
    def test_tight_values_all_used(self):
        rpk, counts = _markers([10.0, 10.5, 9.5, 10.2, 9.8, 10.1, 9.9, 10.0, 10.3, 9.7])
        est = estimate_genome_coverage(rpk, counts)
        assert est.coverage_rpk == pytest.approx(10.0, rel=0.01)
        assert len(est.used) >= 5

    def test_duplicated_marker_removed_as_outlier(self):
        """A marker at 2x (lineage-specific duplication) must not pull G up."""
        rpk, counts = _markers([10.0, 10.2, 9.8, 10.1, 9.9, 10.0, 20.0, 10.0, 9.9, 10.1])
        est = estimate_genome_coverage(rpk, counts)
        assert "m6" in est.outliers
        assert est.coverage_rpk == pytest.approx(10.0, rel=0.02)

    def test_mad_zero_keeps_median_values_only(self):
        rpk, counts = _markers([10.0] * 6 + [30.0])
        est = estimate_genome_coverage(rpk, counts)
        assert est.coverage_rpk == 10.0
        assert est.outliers == ["m6"]

    def test_undetected_markers_skipped(self):
        rpk, counts = _markers([10.0, 10.0, 10.0, 10.0, 10.0, 0.0, None])
        est = estimate_genome_coverage(rpk, counts)
        assert est.coverage_rpk == 10.0
        assert est.undetected == ["m5", "m6"]

    def test_too_few_markers_not_estimable(self):
        rpk, counts = _markers([10.0, 10.0, 10.0, 0.0, 0.0])
        est = estimate_genome_coverage(rpk, counts)
        assert est.coverage_rpk is None
        assert "not estimable" in est.note

    def test_low_depth_flagged(self):
        rpk = {f"m{i}": 1.0 for i in range(6)}
        counts = {m: 2.0 for m in rpk}
        est = estimate_genome_coverage(rpk, counts)
        assert est.coverage_rpk == 1.0
        assert "noisy" in est.note


class TestCopiesPerGenome:
    def test_ratio(self):
        out = copies_per_genome({"a": 20.0, "b": 0.0, "c": None}, 10.0)
        assert out == {"a": 2.0, "b": 0.0, "c": None}

    def test_no_coverage_all_none(self):
        assert copies_per_genome({"a": 20.0}, None) == {"a": None}


# ---------------------------------------------------------------------------
# count_hits(known_families=...)
# ---------------------------------------------------------------------------

def test_count_hits_marker_family_list():
    hits = {"r1": _hit("r1", "leuS|A|x"), "r2": _hit("r2", "leuS_decoy|B|y")}
    counts = count_hits(hits, paired=True, known_families=["leuS", "ffh"])
    assert counts == {"leuS": 0.5, "ffh": 0.0, "leuS_decoy": 0.5}
    assert "ectA" not in counts


# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------

class TestOutputs:
    def test_gene_counts_without_copies_unchanged(self, tmp_path: Path):
        out = write_gene_counts(str(tmp_path / "s"), {"ectA": 2.0}, {"ectA": 1.0},
                                total_reads=10, filtered_reads=5)
        assert out.read_text().splitlines()[2] == "gene\traw_count\tRPM"

    def test_gene_counts_with_copies_column(self, tmp_path: Path):
        out = write_gene_counts(str(tmp_path / "s"), {"ectA": 2.0, "ectB": 0.0, "ectC": 1.0},
                                {"ectA": 1.0, "ectB": 0.0, "ectC": 0.5},
                                total_reads=10, filtered_reads=5,
                                copies_per_genome={"ectA": 0.25, "ectB": 0.0, "ectC": None})
        lines = out.read_text().splitlines()
        assert lines[2] == "gene\traw_count\tRPM\tcopies_per_genome"
        assert lines[3].split("\t")[-1] == "0.25"
        assert lines[4].split("\t")[-1] == "0"
        assert lines[5].split("\t")[-1] == "NA"

    def test_marker_report(self, tmp_path: Path):
        rpk, counts = _markers([10.0, 10.0, 10.0, 10.0, 10.0, 30.0, 0.0])
        est = estimate_genome_coverage(rpk, counts)
        out = write_marker_report(str(tmp_path / "s"), counts, {"m0": 500.0}, rpk, est,
                                  read_length=150)
        text = out.read_text().splitlines()
        assert text[0] == "# read_length\t150"
        assert text[1] == "# genome_coverage_rpk\t10"
        rows = {l.split("\t")[0]: l.split("\t") for l in text if not l.startswith(("#", "marker"))}
        assert rows["m0"][2] == "500.0" and rows["m0"][-1] == "used"
        assert rows["m5"][-1] == "outlier"
        assert rows["m6"][-1] == "undetected"
        assert rows["m1"][2] == "NA"

    def test_aln_stats_extra_rows(self):
        stats = build_aln_stats(mode="profile", total_reads=1, filtered_reads=1, counts={},
                                db="x", min_identity=0.8, min_query_cover=0.8, evalue=1e-5,
                                extra={"genome_coverage_rpk": "10"})
        assert stats["genome_coverage_rpk"] == "10"


# ---------------------------------------------------------------------------
# resolve_refdb: markers/ discovery
# ---------------------------------------------------------------------------

class TestResolveMarkers:
    def _release(self, tmp_path: Path, with_markers: bool, complete: bool = True):
        (tmp_path / "osmo_refdb.dmnd").write_text("")
        (tmp_path / "osmo_refdb.train_refs.faa").write_text(">ectA|A|x\nM\n")
        if with_markers:
            d = tmp_path / "markers"
            d.mkdir()
            (d / "osmo_markers.dmnd").write_text("")
            (d / "osmo_markers.manifest.tsv").write_text("marker\nleuS\n")
            if complete:
                (d / "osmo_markers.train_refs.faa").write_text(">leuS|A|x\nM\n")
        return tmp_path

    def test_old_release_has_no_markers(self, tmp_path: Path):
        refdb = resolve_refdb(str(self._release(tmp_path, False)), mode="profile")
        assert refdb.markers is None
        assert refdb.train_refs is not None

    def test_markers_found_in_profile_mode(self, tmp_path: Path):
        refdb = resolve_refdb(str(self._release(tmp_path, True)), mode="profile")
        assert refdb.markers is not None
        assert refdb.markers.exclude_families is None

    def test_annotate_mode_ignores_markers(self, tmp_path: Path):
        refdb = resolve_refdb(str(self._release(tmp_path, True)), mode="annotate")
        assert refdb.markers is None and refdb.train_refs is None

    def test_incomplete_markers_dir_raises(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError):
            resolve_refdb(str(self._release(tmp_path, True, complete=False)), mode="profile")
