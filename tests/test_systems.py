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


class TestParalogRescue:
    """Positional rescue for osmo_refdb's documented paralog-confusion
    cliques (systems.PARALOG_CONFUSION_CLIQUES) -- see
    systems._rescue_missing_mandatory_families's docstring. Mirrors the real
    finding on B. subtilis/S. aureus/etc.: opuBC gets systematically
    mislabeled opuCC (shared PF04069 substrate-binding domain), but only
    counts as opuBC when positionally confirmed adjacent to a real opuB
    subunit -- never on label alone."""

    def test_rescues_opuBC_via_mislabeled_opuCC(self):
        best_hits = {
            "NODE_1_1": _hit("opuBA", "NODE_1_1"),
            "NODE_1_2": _hit("opuBB", "NODE_1_2"),
            "NODE_1_3": _hit("opuCC", "NODE_1_3"),  # really opuBC, mislabeled
        }
        orf_coords = {
            f"NODE_1_{i}": _coord("NODE_1", i, *_gene("NODE_1", i)) for i in (1, 2, 3)
        }
        calls = detect_systems(best_hits, orf_coords)
        call = calls["opuB"]
        assert call.present is True
        assert call.colocalized is True
        assert "opuBC" in call.mandatory_hits
        assert call.rescued_hits == ("opuBC<-opuCC",)

    def test_rescue_alone_does_not_flip_present_when_another_slot_is_genuinely_missing(self):
        # opuBA has NO hit anywhere (a real recall gap, not a labeling
        # issue) -- rescuing opuBC still can't make the system present.
        best_hits = {
            "NODE_1_1": _hit("opuBB", "NODE_1_1"),
            "NODE_1_2": _hit("opuCC", "NODE_1_2"),  # really opuBC, mislabeled
        }
        orf_coords = {
            f"NODE_1_{i}": _coord("NODE_1", i, *_gene("NODE_1", i)) for i in (1, 2)
        }
        calls = detect_systems(best_hits, orf_coords)
        call = calls["opuB"]
        assert call.present is False
        # but the rescue is still surfaced, precisely the useful signal
        # even though the system overall isn't present
        assert call.rescued_hits == ("opuBC<-opuCC",)
        assert "opuBC" in call.mandatory_hits

    def test_rescue_requires_positional_confirmation_not_label_alone(self):
        # an "opuCC"-labeled hit far from any real opuB subunit must NOT
        # be rescued -- there's no anchor to justify the reattribution.
        best_hits = {
            "NODE_1_1": _hit("opuBA", "NODE_1_1"),
            "NODE_1_2": _hit("opuBB", "NODE_1_2"),
            "NODE_1_500": _hit("opuCC", "NODE_1_500"),  # far away, unrelated
        }
        orf_coords = {
            "NODE_1_1": _coord("NODE_1", 1, *_gene("NODE_1", 1)),
            "NODE_1_2": _coord("NODE_1", 2, *_gene("NODE_1", 2)),
            "NODE_1_500": _coord("NODE_1", 500, *_gene("NODE_1", 500)),
        }
        calls = detect_systems(best_hits, orf_coords)
        call = calls["opuB"]
        assert call.present is False
        assert call.rescued_hits == ()
        assert "opuBC" not in call.mandatory_hits

    def test_trk_rescues_trkH_via_mislabeled_ktrB(self):
        best_hits = {
            "NODE_1_1": _hit("trkA", "NODE_1_1"),
            "NODE_1_2": _hit("ktrB", "NODE_1_2"),  # really trkH, mislabeled
        }
        orf_coords = {
            "NODE_1_1": _coord("NODE_1", 1, *_gene("NODE_1", 1)),
            "NODE_1_2": _coord("NODE_1", 2, *_gene("NODE_1", 2)),
        }
        calls = detect_systems(best_hits, orf_coords)
        call = calls["trk"]
        assert call.present is True
        assert call.rescued_hits == ("trkH<-ktrB",)

    def test_proX_never_rescues_opuAC(self):
        # Regression guard: proX is a documented, deliberately SEPARATE,
        # non-orthologous system (ProU), not a duplicate copy of Opu the
        # way opuAC/opuBC/opuCC genuinely are -- rescuing opuAC from a real
        # proX hit would conflate two distinct systems (e.g. E. coli's own
        # real ProVWX operon getting misread as evidence of Firmicutes
        # OpuA), not correct a labeling ambiguity. proX must never appear
        # in PARALOG_CONFUSION_CLIQUES as a rescue source for opu*C.
        best_hits = {
            "NODE_1_1": _hit("opuAA", "NODE_1_1"),
            "NODE_1_2": _hit("opuAB", "NODE_1_2"),
            "NODE_1_3": _hit("proX", "NODE_1_3"),  # E. coli's real ProU subunit
        }
        orf_coords = {
            f"NODE_1_{i}": _coord("NODE_1", i, *_gene("NODE_1", i)) for i in (1, 2, 3)
        }
        calls = detect_systems(best_hits, orf_coords)
        call = calls["opuA"]
        assert call.present is False
        assert call.rescued_hits == ()
        assert "opuAC" not in call.mandatory_hits


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

    def test_override_preserves_require_strand_consistency(self, tmp_path: Path):
        # regression test: overriding distances must not silently reset
        # require_strand_consistency back to the dataclass default (True)
        config_path = tmp_path / "systems_config.tsv"
        config_path.write_text(
            "system\tmax_distance_genes\tmax_distance_bp\ngalE_murB_mazG\t40\t40000\n"
        )
        systems = load_systems_config(config_path)
        assert systems["galE_murB_mazG"].require_strand_consistency is False
        assert DEFAULT_SYSTEMS["galE_murB_mazG"].require_strand_consistency is False
        # a strict system's default (True) must also survive being untouched
        assert systems["ots"].require_strand_consistency is True


class TestGaleMurBMazG:
    """galE_murB_mazG (Culligan et al. 2012, ISME J -- see systems.py's
    DEFAULT_SYSTEMS entry): a wider genomic-neighborhood claim, not a
    transcriptional operon -- much looser distance window, and strand
    consistency deliberately NOT enforced (a real strand flip sits inside
    the one confirmed positive instance, both in the source fosmid clone
    and in C. aerofaciens ATCC 25986's own reference genome)."""

    def test_present_despite_wide_gap_and_strand_flip(self):
        # mirrors the real C. aerofaciens layout: galE/murB tight and
        # same-strand, mazG ~23 genes away on the OPPOSITE strand
        best_hits = {
            "NODE_1_1": _hit("mazG", "NODE_1_1"),
            "NODE_1_25": _hit("murB", "NODE_1_25"),
            "NODE_1_26": _hit("galE", "NODE_1_26"),
        }
        orf_coords = {
            "NODE_1_1": _coord("NODE_1", 1, *_gene("NODE_1", 1), strand=-1),
            "NODE_1_25": _coord("NODE_1", 25, *_gene("NODE_1", 25), strand=1),
            "NODE_1_26": _coord("NODE_1", 26, *_gene("NODE_1", 26), strand=1),
        }
        calls = detect_systems(best_hits, orf_coords)
        call = calls["galE_murB_mazG"]
        assert call.present is True
        assert call.colocalized is True

    def test_absent_when_beyond_the_wide_window(self):
        default = DEFAULT_SYSTEMS["galE_murB_mazG"]
        far_index = default.max_distance_genes + 5
        best_hits = {
            "NODE_1_1": _hit("galE", "NODE_1_1"),
            "NODE_1_2": _hit("murB", "NODE_1_2"),
            f"NODE_1_{far_index}": _hit("mazG", f"NODE_1_{far_index}"),
        }
        orf_coords = {
            "NODE_1_1": _coord("NODE_1", 1, *_gene("NODE_1", 1)),
            "NODE_1_2": _coord("NODE_1", 2, *_gene("NODE_1", 2)),
            f"NODE_1_{far_index}": _coord("NODE_1", far_index, *_gene("NODE_1", far_index)),
        }
        calls = detect_systems(best_hits, orf_coords)
        assert calls["galE_murB_mazG"].present is False

    def test_strict_systems_still_reject_strand_flips(self):
        # regression guard: require_strand_consistency=False on this ONE
        # system must not weaken the strict systems (e.g. ots)
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
