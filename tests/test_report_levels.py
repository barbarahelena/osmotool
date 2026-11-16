"""
test_report_levels.py — profile's unit table from osmo_refdb's
profile_report_levels.tsv: loading, aggregation, output, release discovery
"""

from __future__ import annotations

from pathlib import Path

import pytest

from osmotool.functions import ReportUnit, load_report_levels
from osmotool.output import write_unit_counts
from osmotool.quantifier import aggregate_report_units
from osmotool.utils import resolve_refdb

MANIFEST = (
    "# osmo_refdb profile report levels: one row per reporting unit\n"
    "# rule: ...\n"
    "unit\treport_level\tmembers\tparalog_group\tcascade_members\tf1_all\n"
    "ectA\tgene\tectA\t\t\t0.796\n"
    "opu_permease\tgroup\topuAB,opuBB,opuCB\topu_permease\topuAB,opuBB,opuCB\t0.527\n"
    "opuAC\tnone\topuAC\topu_substrate_binding\t\t\n"
)


def _manifest(tmp_path: Path, text: str = MANIFEST) -> Path:
    p = tmp_path / "osmo_refdb.profile_report_levels.tsv"
    p.write_text(text)
    return p


# ---------------------------------------------------------------------------
# load_report_levels
# ---------------------------------------------------------------------------

class TestLoadReportLevels:
    def test_reads_units_in_order_skipping_comments(self, tmp_path: Path):
        units = load_report_levels(_manifest(tmp_path))
        assert units == [
            ReportUnit("ectA", "gene", ("ectA",)),
            ReportUnit("opu_permease", "group", ("opuAB", "opuBB", "opuCB")),
            ReportUnit("opuAC", "none", ("opuAC",)),
        ]

    def test_unknown_level_raises(self, tmp_path: Path):
        text = "unit\treport_level\tmembers\nectA\tfamily\tectA\n"
        with pytest.raises(ValueError, match="report_level"):
            load_report_levels(_manifest(tmp_path, text))

    def test_family_in_two_units_raises(self, tmp_path: Path):
        text = "unit\treport_level\tmembers\nectA\tgene\tectA\ng\tgroup\tectA,ectB\n"
        with pytest.raises(ValueError, match="both"):
            load_report_levels(_manifest(tmp_path, text))

    def test_missing_column_raises(self, tmp_path: Path):
        with pytest.raises(ValueError, match="members"):
            load_report_levels(_manifest(tmp_path, "unit\treport_level\nectA\tgene\n"))


# ---------------------------------------------------------------------------
# aggregate_report_units
# ---------------------------------------------------------------------------

UNITS = [
    ReportUnit("ectA", "gene", ("ectA",)),
    ReportUnit("opu_permease", "group", ("opuAB", "opuBB", "opuCB")),
    ReportUnit("opuAC", "none", ("opuAC",)),
]


class TestAggregateReportUnits:
    def _rows(self, counts, copies=None):
        rpm = {k: v * 10 for k, v in counts.items()}
        return {r["unit"]: r for r in aggregate_report_units(UNITS, counts, rpm, copies)}

    def test_group_sums_members(self):
        rows = self._rows({"ectA": 2.0, "opuAB": 1.5, "opuBB": 0.5, "opuCB": 0.0, "opuAC": 3.0})
        assert rows["opu_permease"]["raw_count"] == 2.0
        assert rows["opu_permease"]["rpm"] == 20.0
        assert rows["opu_permease"]["status"] == "reported"
        assert rows["ectA"]["raw_count"] == 2.0

    def test_excluded_unit_has_no_count_even_with_reads(self):
        rows = self._rows({"ectA": 0.0, "opuAC": 3.0})
        assert rows["opuAC"]["status"] == "excluded"
        assert rows["opuAC"]["raw_count"] is None and rows["opuAC"]["rpm"] is None

    def test_every_unit_present_with_zero_reads(self):
        rows = self._rows({})
        assert list(rows) == ["ectA", "opu_permease", "opuAC"]
        assert rows["ectA"]["raw_count"] == 0.0

    def test_copies_sum_and_unknown_propagates(self):
        counts = {"ectA": 1.0, "opuAB": 1.0, "opuBB": 1.0, "opuCB": 0.0}
        rows = self._rows(counts, {"ectA": None, "opuAB": 0.25, "opuBB": 0.5, "opuCB": 0.0})
        assert rows["opu_permease"]["copies_per_genome"] == 0.75
        assert rows["ectA"]["copies_per_genome"] is None

    def test_zero_read_member_without_copies_entry_counts_zero(self):
        rows = self._rows({"opuAB": 1.0, "opuBB": 0.0}, {"opuAB": 0.25})
        assert rows["opu_permease"]["copies_per_genome"] == 0.25

    def test_label_outside_manifest_with_reads_is_listed(self):
        rows = self._rows({"ectA": 1.0, "mrpA_mrpB_fused": 2.0, "kdpA": 0.0})
        assert rows["mrpA_mrpB_fused"]["status"] == "not_in_manifest"
        assert rows["mrpA_mrpB_fused"]["raw_count"] == 2.0
        assert "kdpA" not in rows      # no reads: nothing to announce


# ---------------------------------------------------------------------------
# write_unit_counts
# ---------------------------------------------------------------------------

class TestWriteUnitCounts:
    def _write(self, tmp_path: Path, copies=None) -> list[str]:
        counts = {"ectA": 2.0, "opuAB": 1.5, "opuBB": 0.5, "opuAC": 3.0}
        rows = aggregate_report_units(UNITS, counts, {k: v * 10 for k, v in counts.items()}, copies)
        path = write_unit_counts(tmp_path / "s", rows, total_reads=100, filtered_reads=7,
                                 report_levels="rel/levels.tsv", with_copies=copies is not None)
        assert path == tmp_path / "s.unit_counts.tsv"
        return path.read_text().splitlines()

    def test_format(self, tmp_path: Path):
        lines = self._write(tmp_path)
        assert lines[:3] == ["# total_reads\t100", "# filtered_reads\t7",
                             "# report_levels\trel/levels.tsv"]
        assert lines[3] == "unit\treport_level\tmembers\tstatus\traw_count\tRPM"
        assert lines[4] == "ectA\tgene\tectA\treported\t2.0\t20"
        assert lines[5] == "opu_permease\tgroup\topuAB,opuBB,opuCB\treported\t2.0\t20"
        assert lines[6] == "opuAC\tnone\topuAC\texcluded\tNA\tNA"

    def test_copies_column(self, tmp_path: Path):
        lines = self._write(tmp_path, {"ectA": 0.5, "opuAB": 0.1, "opuBB": 0.2, "opuAC": 1.0})
        assert lines[3].endswith("\tcopies_per_genome")
        assert lines[5].endswith("\t0.3")
        assert lines[6].endswith("\tNA\tNA\tNA")


# ---------------------------------------------------------------------------
# resolve_refdb: report-level manifest discovery
# ---------------------------------------------------------------------------

class TestResolveReportLevels:
    def _release(self, tmp_path: Path, with_levels: bool) -> str:
        (tmp_path / "osmo_refdb.dmnd").write_text("")
        if with_levels:
            _manifest(tmp_path)
        return str(tmp_path)

    def test_found_in_profile_mode(self, tmp_path: Path):
        refdb = resolve_refdb(self._release(tmp_path, True), mode="profile")
        assert refdb.report_levels == tmp_path / "osmo_refdb.profile_report_levels.tsv"

    def test_absent_in_old_release(self, tmp_path: Path):
        assert resolve_refdb(self._release(tmp_path, False), mode="profile").report_levels is None

    def test_annotate_mode_ignores_it(self, tmp_path: Path):
        assert resolve_refdb(self._release(tmp_path, True), mode="annotate").report_levels is None
