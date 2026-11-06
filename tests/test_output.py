"""
test_output.py — unit tests for osmotool.output
"""

from __future__ import annotations

from pathlib import Path

import pytest

from osmotool.output import (
    write_gene_counts,
    write_gene_coordinates,
    write_aln_stats,
    build_aln_stats,
    write_system_calls,
    PROFILE_SYSTEMS_NOTE,
)
from osmotool.functions import AssemblyStats
from osmotool.systems import SystemCall


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

    def test_assembly_stats_lines(self, tmp_path: Path, sample_counts, sample_rpm):
        common = dict(total_reads=4128, filtered_reads=187,
                      normalisation="copies_per_kb", total_reads_label="total_proteins")
        plain = write_gene_counts(str(tmp_path / "plain"), sample_counts, sample_rpm, **common)
        out = write_gene_counts(
            str(tmp_path / "stats"), sample_counts, sample_rpm, **common,
            assembly_stats=AssemblyStats(genome_bp=4641652, gc=0.507907098,
                                         n_contigs=1, contig_n50=4641652),
        )
        lines = out.read_text().splitlines()
        assert lines[:6] == [
            "# total_proteins\t4128",
            "# genome_bp\t4641652",
            "# gc\t0.50791",
            "# n_contigs\t1",
            "# contig_n50\t4641652",
            "# filtered_reads\t187",
        ]
        # everything except the four new lines is byte-identical
        assert lines[:1] + lines[5:] == plain.read_text().splitlines()

    def test_assembly_stats_gc_na(self, tmp_path: Path, sample_counts, sample_rpm):
        out = write_gene_counts(
            str(tmp_path / "s"), sample_counts, sample_rpm,
            total_reads=1, filtered_reads=0,
            assembly_stats=AssemblyStats(genome_bp=4, gc=None, n_contigs=1, contig_n50=4),
        )
        assert "# gc\tNA" in out.read_text().splitlines()

    def test_no_assembly_stats_by_default(self, tmp_path: Path, sample_counts, sample_rpm):
        out = write_gene_counts(
            str(tmp_path / "s"), sample_counts, sample_rpm,
            total_reads=6, filtered_reads=6,
        )
        comments = [l for l in out.read_text().splitlines() if l.startswith("#")]
        assert comments == ["# total_reads\t6", "# filtered_reads\t6"]

    def test_proteins_with_hits_line(self, tmp_path: Path, sample_counts, sample_rpm):
        common = dict(total_reads=4128, filtered_reads=187,
                      normalisation="copies_per_kb", total_reads_label="total_proteins")
        plain = write_gene_counts(str(tmp_path / "plain"), sample_counts, sample_rpm, **common)
        out = write_gene_counts(str(tmp_path / "pwh"), sample_counts, sample_rpm, **common,
                                proteins_with_hits=187)
        lines = out.read_text().splitlines()
        assert lines[1:3] == ["# filtered_reads\t187", "# proteins_with_hits\t187"]
        assert lines[:2] + lines[3:] == plain.read_text().splitlines()

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
# write_gene_coordinates
# ---------------------------------------------------------------------------

class TestWriteGeneCoordinates:
    def _rows(self):
        return [
            {"gene": "ectB", "protein_id": "NODE_1_8", "contig": "NODE_1",
             "start": 2410, "end": 3200, "strand": 1, "gene_index": 8},
            {"gene": "ectA", "protein_id": "NODE_1_7", "contig": "NODE_1",
             "start": 1801, "end": 2400, "strand": 1, "gene_index": 7},
        ]

    def test_creates_file(self, tmp_path: Path):
        prefix = str(tmp_path / "sample")
        out = write_gene_coordinates(prefix, self._rows())
        assert out.exists()

    def test_output_filename(self, tmp_path: Path):
        prefix = str(tmp_path / "mysample")
        out = write_gene_coordinates(prefix, self._rows())
        assert out.name == "mysample.gene_coordinates.tsv"

    def test_header_row(self, tmp_path: Path):
        prefix = str(tmp_path / "sample")
        out = write_gene_coordinates(prefix, self._rows())
        lines = out.read_text().splitlines()
        assert lines[0] == "gene\tprotein_id\tcontig\tstart\tend\tstrand\tgene_index"

    def test_rows_sorted_by_gene_then_start(self, tmp_path: Path):
        prefix = str(tmp_path / "sample")
        out = write_gene_coordinates(prefix, self._rows())
        lines = out.read_text().splitlines()[1:]
        assert lines[0].startswith("ectA\tNODE_1_7\tNODE_1\t1801\t2400\t1\t7")
        assert lines[1].startswith("ectB\tNODE_1_8\tNODE_1\t2410\t3200\t1\t8")

    def test_empty_rows(self, tmp_path: Path):
        prefix = str(tmp_path / "sample")
        out = write_gene_coordinates(prefix, [])
        lines = out.read_text().splitlines()
        assert lines == ["gene\tprotein_id\tcontig\tstart\tend\tstrand\tgene_index"]


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

    def test_note_systems_present_when_passed(self):
        stats = self._stats(note_systems=PROFILE_SYSTEMS_NOTE)
        assert stats["note_systems"] == PROFILE_SYSTEMS_NOTE

    def test_note_systems_absent_by_default(self):
        stats = self._stats()
        assert "note_systems" not in stats


# ---------------------------------------------------------------------------
# write_system_calls
# ---------------------------------------------------------------------------

class TestWriteSystemCalls:
    def _calls(self) -> dict[str, SystemCall]:
        return {
            "opuA": SystemCall(
                system="opuA", present=True, subunit_count=3, colocalized=True,
                fragmented_assembly_flag=False, contig="NODE_1",
                mandatory_hits=("opuAA", "opuAB", "opuAC"), accessory_hits=(),
            ),
            "mrp": SystemCall(
                system="mrp", present=False, subunit_count=2, colocalized=False,
                fragmented_assembly_flag=True, contig=None,
                mandatory_hits=("mrpA", "mrpC"), accessory_hits=(),
            ),
            "opuB": SystemCall(
                system="opuB", present=False, subunit_count=2, colocalized=False,
                fragmented_assembly_flag=False, contig=None,
                mandatory_hits=("opuBB", "opuBC"), accessory_hits=(),
                rescued_hits=("opuBC<-opuCC",),
            ),
        }

    def test_creates_file(self, tmp_path: Path):
        prefix = str(tmp_path / "sample")
        out = write_system_calls(prefix, self._calls())
        assert out.exists()

    def test_output_filename(self, tmp_path: Path):
        prefix = str(tmp_path / "mysample")
        out = write_system_calls(prefix, self._calls())
        assert out.name == "mysample.systems.tsv"

    def test_header_row(self, tmp_path: Path):
        prefix = str(tmp_path / "sample")
        out = write_system_calls(prefix, self._calls())
        lines = out.read_text().splitlines()
        assert lines[0] == (
            "system\tpresent\tsubunit_count\tcolocalized\t"
            "fragmented_assembly_flag\tcontig\tmandatory_hits\taccessory_hits\trescued_hits"
        )

    def test_rows_sorted_and_populated(self, tmp_path: Path):
        prefix = str(tmp_path / "sample")
        out = write_system_calls(prefix, self._calls())
        lines = out.read_text().splitlines()[1:]
        assert lines[0].startswith("mrp\tFalse\t2\tFalse\tTrue\t-\tmrpA,mrpC\t-\t-")
        assert lines[1].startswith("opuA\tTrue\t3\tTrue\tFalse\tNODE_1\topuAA,opuAB,opuAC\t-\t-")

    def test_rescued_hits_column(self, tmp_path: Path):
        prefix = str(tmp_path / "sample")
        out = write_system_calls(prefix, self._calls())
        content = out.read_text()
        assert "opuBC<-opuCC" in content
