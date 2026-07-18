"""
test_output.py — unit tests for osmotool.output
"""

from __future__ import annotations

from pathlib import Path

import pytest

from osmotool.output import (
    write_gene_counts,
    write_aln_stats,
    build_aln_stats,
)


# ---------------------------------------------------------------------------
# write_gene_counts
# ---------------------------------------------------------------------------

class TestWriteGeneCounts:
    def test_creates_file(self, tmp_path: Path, sample_counts, sample_rpm):
        prefix = str(tmp_path / "sample")
        out = write_gene_counts(
            prefix, sample_counts, sample_rpm,
            total_reads=6, filtered_reads=6, normalisation="rpm",
        )
        assert out.exists()

    def test_header_rows(self, tmp_path: Path, sample_counts, sample_rpm):
        prefix = str(tmp_path / "sample")
        out = write_gene_counts(
            prefix, sample_counts, sample_rpm,
            total_reads=6, filtered_reads=6, normalisation="rpm",
        )
        lines = out.read_text().splitlines()
        assert lines[0].startswith("# total_reads\t6")
        assert lines[1].startswith("# filtered_reads\t6")

    def test_column_header(self, tmp_path: Path, sample_counts, sample_rpm):
        prefix = str(tmp_path / "sample")
        out = write_gene_counts(
            prefix, sample_counts, sample_rpm,
            total_reads=6, filtered_reads=6, normalisation="rpm",
        )
        lines = out.read_text().splitlines()
        assert lines[2] == "gene\traw_count\tRPM"

    def test_all_families_present(self, tmp_path: Path, sample_counts, sample_rpm):
        prefix = str(tmp_path / "sample")
        out = write_gene_counts(
            prefix, sample_counts, sample_rpm,
            total_reads=6, filtered_reads=6,
        )
        content = out.read_text()
        for fam in ("ectA", "ectB", "ectC", "betL", "kdpA", "nhaA"):
            assert fam in content

    def test_rpkm_column_label(self, tmp_path: Path, sample_counts, sample_rpm):
        prefix = str(tmp_path / "sample")
        out = write_gene_counts(
            prefix, sample_counts, sample_rpm,
            total_reads=6, filtered_reads=6, normalisation="rpkm",
        )
        lines = out.read_text().splitlines()
        assert "RPKM" in lines[2]

    def test_output_filename(self, tmp_path: Path, sample_counts, sample_rpm):
        prefix = str(tmp_path / "mysample")
        out = write_gene_counts(
            prefix, sample_counts, sample_rpm,
            total_reads=100, filtered_reads=10,
        )
        assert out.name == "mysample.gene_counts.tsv"

    def test_total_reads_label_override(self, tmp_path: Path, sample_counts, sample_rpm):
        prefix = str(tmp_path / "sample")
        out = write_gene_counts(
            prefix, sample_counts, sample_rpm,
            total_reads=4128, filtered_reads=187,
            normalisation="copies_per_kb", total_reads_label="total_proteins",
        )
        lines = out.read_text().splitlines()
        assert lines[0] == "# total_proteins\t4128"
        assert lines[2] == "gene\traw_count\tCOPIES_PER_KB"

    def test_raw_count_values(self, tmp_path: Path, sample_counts, sample_rpm):
        prefix = str(tmp_path / "sample")
        out = write_gene_counts(
            prefix, sample_counts, sample_rpm,
            total_reads=6, filtered_reads=6,
        )
        lines = [l for l in out.read_text().splitlines() if not l.startswith("#")]
        data = {parts[0]: float(parts[1]) for parts in (l.split("\t") for l in lines[1:])}
        for fam in sample_counts:
            assert data[fam] == pytest.approx(sample_counts[fam])


# ---------------------------------------------------------------------------
# write_aln_stats
# ---------------------------------------------------------------------------

class TestWriteAlnStats:
    def test_creates_file(self, tmp_path: Path):
        prefix = str(tmp_path / "sample")
        out = write_aln_stats(prefix, {"mode": "profile", "total_reads": 100})
        assert out.exists()

    def test_header_row(self, tmp_path: Path):
        prefix = str(tmp_path / "sample")
        out = write_aln_stats(prefix, {"mode": "profile"})
        lines = out.read_text().splitlines()
        assert lines[0] == "metric\tvalue"

    def test_key_value_written(self, tmp_path: Path):
        prefix = str(tmp_path / "sample")
        out = write_aln_stats(prefix, {"total_reads": 42, "mode": "annotate"})
        content = out.read_text()
        assert "total_reads\t42" in content
        assert "mode\tannotate" in content

    def test_output_filename(self, tmp_path: Path):
        prefix = str(tmp_path / "mysample")
        out = write_aln_stats(prefix, {})
        assert out.name == "mysample.aln_stats.tsv"


# ---------------------------------------------------------------------------
# build_aln_stats
# ---------------------------------------------------------------------------

class TestBuildAlnStats:
    def _stats(self, **kwargs):
        defaults = dict(
            mode="profile",
            total_reads=1000,
            filtered_reads=50,
            counts={"ectA": 20.0, "ectB": 10.0, "ectC": 5.0,
                    "betL": 5.0, "kdpA": 5.0, "nhaA": 5.0},
            db="/path/to/osmodiamond.dmnd",
            min_identity=0.80,
            min_query_cover=0.80,
            evalue=1e-5,
        )
        defaults.update(kwargs)
        return build_aln_stats(**defaults)

    def test_required_keys_present(self):
        stats = self._stats()
        for key in ("mode", "total_reads", "filtered_reads",
                    "reads_assigned", "pct_reads_assigned",
                    "min_identity", "min_query_cover", "evalue"):
            assert key in stats

    def test_per_family_keys(self):
        stats = self._stats()
        for fam in ("ectA", "ectB", "ectC", "betL", "kdpA", "nhaA"):
            assert f"count_{fam}" in stats

    def test_pct_assigned_calculation(self):
        counts = {"ectA": 50.0, "ectB": 0.0, "ectC": 0.0,
                  "betL": 0.0, "kdpA": 0.0, "nhaA": 0.0}
        stats = build_aln_stats(
            mode="profile", total_reads=1000, filtered_reads=100,
            counts=counts, db="db", min_identity=0.8,
            min_query_cover=0.8, evalue=1e-5,
        )
        assert float(stats["pct_reads_assigned"]) == pytest.approx(50.0)

    def test_zero_filtered_reads_no_division_error(self):
        counts = {"ectA": 0.0}
        stats = build_aln_stats(
            mode="profile", total_reads=0, filtered_reads=0,
            counts=counts, db="db", min_identity=0.8,
            min_query_cover=0.8, evalue=1e-5,
        )
        assert float(stats["pct_reads_assigned"]) == pytest.approx(0.0)
