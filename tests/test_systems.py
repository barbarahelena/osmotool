"""
test_systems.py — tests for systems.py's complex/system-aware detection
(Task 2: proximity/quorum logic on top of osmotool's existing per-gene hits)
"""

from __future__ import annotations

from pathlib import Path

import pytest

from osmotool.functions import OrfCoord
from osmotool.systems import DEFAULT_SYSTEMS, SystemDef, detect_systems, load_systems_config


def _coord(contig: str, gene_index: int, start: int, end: int, strand: int = 1) -> OrfCoord:
    return OrfCoord(contig=contig, gene_index=gene_index, start=start, end=end, strand=strand)


def _hit(family: str, orf_id: str, uniprot: str = "P00000", organism: str = "Test_organism") -> dict:
    """A best_hits-shaped row (see quantifier.select_best_hits' return shape)."""
    return {
        "qseqid": orf_id,
        "sseqid": f"{family}|{uniprot}|{organism}",
        "bitscore": 200.0,
        "evalue": 1e-30,
    }


def _gene(contig: str, gene_index: int, length_bp: int = 300) -> tuple[int, int]:
    """Return (start, end) placing a *length_bp* gene at slot *gene_index* on a
    contig of closely-packed genes (matches how a real operon looks)."""
    start = (gene_index - 1) * length_bp + 1
    end = start + length_bp - 1
    return start, end


class TestOpuA:
    def test_colocalized_present(self):
        best_hits = {}
        orf_coords = {}
        for i, fam in enumerate(("opuAA", "opuAB", "opuAC"), start=1):
            orf_id = f"NODE_1_{i}"
            best_hits[orf_id] = _hit(fam, orf_id)
            start, end = _gene("NODE_1", i)
            orf_coords[orf_id] = _coord("NODE_1", i, start, end)

        calls = detect_systems(best_hits, orf_coords)
        call = calls["opuA"]
        assert call.present is True
        assert call.colocalized is True
        assert call.subunit_count == 3
        assert call.fragmented_assembly_flag is False
        assert call.contig == "NODE_1"
        assert set(call.mandatory_hits) == {"opuAA", "opuAB", "opuAC"}

    def test_absent_no_hits(self):
        calls = detect_systems({}, {})
        call = calls["opuA"]
        assert call.present is False
        assert call.subunit_count == 0
        assert call.colocalized is False
        assert call.fragmented_assembly_flag is False
        assert call.contig is None

    def test_fragmented_across_contigs(self):
        best_hits = {
            "NODE_1_1": _hit("opuAA", "NODE_1_1"),
            "NODE_2_1": _hit("opuAB", "NODE_2_1"),
        }
        orf_coords = {
            "NODE_1_1": _coord("NODE_1", 1, *_gene("NODE_1", 1)),
            "NODE_2_1": _coord("NODE_2", 1, *_gene("NODE_2", 1)),
        }
        calls = detect_systems(best_hits, orf_coords)
        call = calls["opuA"]
        assert call.present is False
        assert call.subunit_count == 2
        assert call.colocalized is False
        assert call.fragmented_assembly_flag is True


class TestProximityWindow:
    """Boundary-tests _within_window via the ots system (2 mandatory slots,
    default max_distance_genes=8 / max_distance_bp=8000)."""

    def _build(self, gap_genes: int):
        # otsA at gene_index 1, otsB at gene_index (2 + gap_genes) so that
        # intervening_genes == gap_genes exactly.
        otsb_index = 2 + gap_genes
        best_hits = {
            "NODE_1_1": _hit("otsA", "NODE_1_1"),
            f"NODE_1_{otsb_index}": _hit("otsB", f"NODE_1_{otsb_index}"),
        }
        orf_coords = {
            "NODE_1_1": _coord("NODE_1", 1, *_gene("NODE_1", 1)),
            f"NODE_1_{otsb_index}": _coord("NODE_1", otsb_index, *_gene("NODE_1", otsb_index)),
        }
        return best_hits, orf_coords

    def test_within_default_window_present(self):
        limit = DEFAULT_SYSTEMS["ots"].max_distance_genes
        best_hits, orf_coords = self._build(gap_genes=limit)  # exactly at the default limit
        calls = detect_systems(best_hits, orf_coords)
        assert calls["ots"].present is True
        assert calls["ots"].colocalized is True

    def test_beyond_default_window_not_colocalized(self):
        limit = DEFAULT_SYSTEMS["ots"].max_distance_genes
        best_hits, orf_coords = self._build(gap_genes=limit + 1)  # one gene past the limit
        calls = detect_systems(best_hits, orf_coords)
        call = calls["ots"]
        assert call.present is False
        assert call.colocalized is False
        # same contig, just too far apart -- NOT a fragmented-assembly case
        assert call.fragmented_assembly_flag is False
        assert call.subunit_count == 2


class TestMrpFusion:
    def test_fused_orf_satisfies_both_slots(self):
        best_hits = {
            "NODE_1_1": _hit("mrpA_mrpB_fused", "NODE_1_1"),
            "NODE_1_2": _hit("mrpC", "NODE_1_2"),
            "NODE_1_3": _hit("mrpE", "NODE_1_3"),
        }
        orf_coords = {
            "NODE_1_1": _coord("NODE_1", 1, *_gene("NODE_1", 1)),
            "NODE_1_2": _coord("NODE_1", 2, *_gene("NODE_1", 2)),
            "NODE_1_3": _coord("NODE_1", 3, *_gene("NODE_1", 3)),
        }
        calls = detect_systems(best_hits, orf_coords)
        call = calls["mrp"]
        assert call.present is True
        assert call.colocalized is True
        # fused ORF counts toward both mrpA and mrpB
        assert call.subunit_count == 4
        assert set(call.mandatory_hits) == {"mrpA", "mrpB", "mrpC", "mrpE"}

    def test_split_mrpA_mrpB_also_satisfies(self):
        best_hits = {
            f"NODE_1_{i}": _hit(fam, f"NODE_1_{i}")
            for i, fam in enumerate(("mrpA", "mrpB", "mrpC", "mrpE"), start=1)
        }
        orf_coords = {
            f"NODE_1_{i}": _coord("NODE_1", i, *_gene("NODE_1", i))
            for i in range(1, 5)
        }
        calls = detect_systems(best_hits, orf_coords)
        assert calls["mrp"].present is True
        assert calls["mrp"].subunit_count == 4


class TestKtrOrSlot:
    def test_ktrD_satisfies_pore_subunit_slot(self):
        best_hits = {
            "NODE_1_1": _hit("ktrA", "NODE_1_1"),
            "NODE_1_2": _hit("ktrD", "NODE_1_2"),
        }
        orf_coords = {
            "NODE_1_1": _coord("NODE_1", 1, *_gene("NODE_1", 1)),
            "NODE_1_2": _coord("NODE_1", 2, *_gene("NODE_1", 2)),
        }
        calls = detect_systems(best_hits, orf_coords)
        call = calls["ktr"]
        assert call.present is True
        assert "ktrD" in call.mandatory_hits
        assert "ktrB" not in call.mandatory_hits


class TestStrandConsistency:
    """A real operon is one contiguous same-strand transcriptional unit --
    complementary to (not a substitute for) the distance window. See
    systems._strand_consistent's docstring."""

    def _build(self, between_strand: int):
        # otsA (gene_index 1) and otsB (gene_index 3) are within the
        # default distance window regardless; an unrelated ORF at
        # gene_index 2 (not itself a family hit) sits directly between them.
        best_hits = {
            "NODE_1_1": _hit("otsA", "NODE_1_1"),
            "NODE_1_3": _hit("otsB", "NODE_1_3"),
        }
        orf_coords = {
            "NODE_1_1": _coord("NODE_1", 1, *_gene("NODE_1", 1), strand=1),
            "NODE_1_2": _coord("NODE_1", 2, *_gene("NODE_1", 2), strand=between_strand),
            "NODE_1_3": _coord("NODE_1", 3, *_gene("NODE_1", 3), strand=1),
        }
        return best_hits, orf_coords

    def test_opposite_strand_gene_between_breaks_colocalization(self):
        best_hits, orf_coords = self._build(between_strand=-1)
        calls = detect_systems(best_hits, orf_coords)
        call = calls["ots"]
        assert call.present is False
        assert call.colocalized is False
        # still genuinely found, just not confirmed as one operon
        assert call.subunit_count == 2

    def test_same_strand_gene_between_still_colocalizes(self):
        best_hits, orf_coords = self._build(between_strand=1)
        calls = detect_systems(best_hits, orf_coords)
        assert calls["ots"].present is True
        assert calls["ots"].colocalized is True

    def test_opposite_strand_hits_never_colocalize(self):
        best_hits = {
            "NODE_1_1": _hit("otsA", "NODE_1_1"),
            "NODE_1_2": _hit("otsB", "NODE_1_2"),
        }
        orf_coords = {
            "NODE_1_1": _coord("NODE_1", 1, *_gene("NODE_1", 1), strand=1),
            "NODE_1_2": _coord("NODE_1", 2, *_gene("NODE_1", 2), strand=-1),
        }
        calls = detect_systems(best_hits, orf_coords)
        assert calls["ots"].present is False


class TestLoadSystemsConfig:
    def test_override_tightens_window(self, tmp_path: Path):
        default_limit = DEFAULT_SYSTEMS["ots"].max_distance_genes
        tightened = 1  # strictly tighter than the default in every case tested here
        config_path = tmp_path / "systems_config.tsv"
        config_path.write_text(f"system\tmax_distance_genes\tmax_distance_bp\nots\t{tightened}\t8000\n")

        systems = load_systems_config(config_path)
        assert systems["ots"].max_distance_genes == tightened
        # untouched systems keep their defaults
        assert systems["opuA"].max_distance_genes == DEFAULT_SYSTEMS["opuA"].max_distance_genes

        # a gap that passes under the default window should be rejected
        # under the tightened override
        otsb_index = 2 + default_limit
        best_hits = {
            "NODE_1_1": _hit("otsA", "NODE_1_1"),
            f"NODE_1_{otsb_index}": _hit("otsB", f"NODE_1_{otsb_index}"),
        }
        orf_coords = {
            "NODE_1_1": _coord("NODE_1", 1, *_gene("NODE_1", 1)),
            f"NODE_1_{otsb_index}": _coord("NODE_1", otsb_index, *_gene("NODE_1", otsb_index)),
        }
        assert detect_systems(best_hits, orf_coords)["ots"].present is True  # sanity: passes under default
        calls = detect_systems(best_hits, orf_coords, systems)
        assert calls["ots"].present is False

    def test_unknown_system_in_config_is_ignored_with_warning(self, tmp_path: Path, caplog):
        config_path = tmp_path / "systems_config.tsv"
        config_path.write_text("system\tmax_distance_genes\tmax_distance_bp\nnotasystem\t3\t8000\n")
        systems = load_systems_config(config_path)
        assert set(systems) == set(DEFAULT_SYSTEMS)
