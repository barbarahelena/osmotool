"""
test_quantifier.py — unit tests for osmotool.quantifier
"""

from __future__ import annotations

import pytest

from osmotool.quantifier import (
    select_best_hits,
    count_hits,
    compute_rpm,
    compute_rpkm,
    estimate_family_lengths,
    KNOWN_FAMILIES,
)


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

def _row(qid, sid, bitscore=200.0, evalue=1e-20, pident=90.0, qcovhsp=85.0, length=150):
    return {
        "qseqid": qid,
        "sseqid": sid,
        "pident": pident,
        "qcovhsp": qcovhsp,
        "length": length,
        "evalue": evalue,
        "bitscore": bitscore,
    }


# ---------------------------------------------------------------------------
# select_best_hits
# ---------------------------------------------------------------------------

class TestSelectBestHits:
    def test_single_hit_returned(self):
        rows = [_row("r1", "ectA|P1|Org")]
        best = select_best_hits(iter(rows))
        assert best["r1"]["sseqid"] == "ectA|P1|Org"

    def test_best_bitscore_wins(self):
        rows = [
            _row("r1", "ectA|P1|Org", bitscore=100.0),
            _row("r1", "ectB|P2|Org", bitscore=200.0),
        ]
        best = select_best_hits(iter(rows))
        assert best["r1"]["sseqid"] == "ectB|P2|Org"

    def test_tied_bitscore_lower_evalue_wins(self):
        rows = [
            _row("r1", "ectA|P1|Org", bitscore=200.0, evalue=1e-10),
            _row("r1", "ectB|P2|Org", bitscore=200.0, evalue=1e-20),
        ]
        best = select_best_hits(iter(rows))
        assert best["r1"]["sseqid"] == "ectB|P2|Org"

    def test_tied_bitscore_evalue_lexicographic_tiebreak(self):
        rows = [
            _row("r1", "ectC|P3|Org", bitscore=200.0, evalue=1e-20),
            _row("r1", "ectA|P1|Org", bitscore=200.0, evalue=1e-20),
        ]
        best = select_best_hits(iter(rows))
        # 'ectA|P1|Org' < 'ectC|P3|Org' lexicographically
        assert best["r1"]["sseqid"] == "ectA|P1|Org"

    def test_different_reads_independent(self):
        rows = [
            _row("r1", "ectA|P1|Org", bitscore=200.0),
            _row("r2", "ectB|P2|Org", bitscore=150.0),
        ]
        best = select_best_hits(iter(rows))
        assert len(best) == 2
        assert best["r1"]["sseqid"] == "ectA|P1|Org"
        assert best["r2"]["sseqid"] == "ectB|P2|Org"

    def test_empty_input_returns_empty(self):
        assert select_best_hits(iter([])) == {}


# ---------------------------------------------------------------------------
# count_hits
# ---------------------------------------------------------------------------

class TestCountHits:
    def _best_hits(self, assignments: dict[str, str]) -> dict[str, dict]:
        """Build a minimal best_hits dict from {read_id: subject_id}."""
        return {qid: _row(qid, sid) for qid, sid in assignments.items()}

    def test_all_families_present_in_output(self):
        counts = count_hits({})
        assert set(counts.keys()) == set(KNOWN_FAMILIES)

    def test_single_end_weight_one(self):
        best = self._best_hits({"r1": "ectA|P1|Org"})
        counts = count_hits(best, paired=False)
        assert counts["ectA"] == 1.0

    def test_paired_weight_half(self):
        best = self._best_hits({"r1": "ectA|P1|Org", "r2": "ectA|P2|Org"})
        counts = count_hits(best, paired=True)
        assert counts["ectA"] == pytest.approx(1.0)  # 2 reads × 0.5

    def test_multi_family_counts(self):
        best = self._best_hits({
            "r1": "ectA|P1|Org",
            "r2": "ectB|P2|Org",
            "r3": "ectA|P3|Org",
        })
        counts = count_hits(best)
        assert counts["ectA"] == pytest.approx(2.0)
        assert counts["ectB"] == pytest.approx(1.0)
        assert counts["ectC"] == pytest.approx(0.0)

    def test_unknown_family_not_in_known_but_counted(self):
        """Unknown families beyond KNOWN_FAMILIES are still accumulated."""
        best = self._best_hits({"r1": "newGene|P9|Org"})
        counts = count_hits(best)
        assert counts.get("newGene", 0.0) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# compute_rpm
# ---------------------------------------------------------------------------

class TestComputeRpm:
    def test_basic_rpm(self):
        counts = {"ectA": 10.0, "ectB": 0.0}
        rpm = compute_rpm(counts, total_reads=1_000_000)
        assert rpm["ectA"] == pytest.approx(10.0)
        assert rpm["ectB"] == pytest.approx(0.0)

    def test_zero_total_reads_returns_zeros(self):
        counts = {"ectA": 5.0}
        rpm = compute_rpm(counts, total_reads=0)
        assert rpm["ectA"] == pytest.approx(0.0)

    def test_formula(self):
        # RPM = count * 1e6 / total
        counts = {"ectA": 50.0}
        rpm = compute_rpm(counts, total_reads=500_000)
        assert rpm["ectA"] == pytest.approx(100.0)


# ---------------------------------------------------------------------------
# compute_rpkm
# ---------------------------------------------------------------------------

class TestComputeRpkm:
    def test_basic_rpkm(self):
        counts = {"ectA": 10.0}
        lengths = {"ectA": 100.0}   # 100 aa → 300 bp
        rpkm = compute_rpkm(counts, total_reads=1_000, gene_lengths_aa=lengths)
        # RPKM = 10 * 1e9 / (300 * 1000) = 33333.3...
        expected = (10 * 1e9) / (300 * 1_000)
        assert rpkm["ectA"] == pytest.approx(expected)

    def test_missing_length_falls_back_to_rpm(self):
        counts = {"ectA": 10.0}
        rpkm = compute_rpkm(counts, total_reads=1_000_000, gene_lengths_aa={})
        # Falls back to RPM
        assert rpkm["ectA"] == pytest.approx(10.0)

    def test_zero_total_reads(self):
        counts = {"ectA": 5.0}
        rpkm = compute_rpkm(counts, total_reads=0, gene_lengths_aa={"ectA": 150.0})
        assert rpkm["ectA"] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# estimate_family_lengths
# ---------------------------------------------------------------------------

class TestEstimateFamilyLengths:
    # protein_lengths is keyed by the *query* protein ID (the assembly's own
    # Prodigal-called ORF), never by the reference database subject ID it
    # aligned to -- see runners.py: parse_protein_lengths() parses the query
    # proteins file, not the reference db. A dict keyed by sseqid here would
    # never match and RPKM would silently fall back to RPM every time.
    def test_single_hit(self):
        best = {"r1": _row("r1", "ectA|P1|Org")}
        prot_lens = {"r1": 200}
        lengths = estimate_family_lengths(best, prot_lens)
        assert lengths["ectA"] == pytest.approx(200.0)

    def test_mean_across_multiple_hits(self):
        best = {
            "r1": _row("r1", "ectA|P1|Org"),
            "r2": _row("r2", "ectA|P2|Org"),
        }
        prot_lens = {"r1": 100, "r2": 200}
        lengths = estimate_family_lengths(best, prot_lens)
        assert lengths["ectA"] == pytest.approx(150.0)

    def test_missing_protein_ignored(self):
        best = {"r1": _row("r1", "ectA|P1|Org")}
        lengths = estimate_family_lengths(best, {})
        assert "ectA" not in lengths

    def test_multi_family(self):
        best = {
            "r1": _row("r1", "ectA|P1|Org"),
            "r2": _row("r2", "ectB|P2|Org"),
        }
        prot_lens = {"r1": 120, "r2": 360}
        lengths = estimate_family_lengths(best, prot_lens)
        assert lengths["ectA"] == pytest.approx(120.0)
        assert lengths["ectB"] == pytest.approx(360.0)

    def test_reference_id_keyed_dict_never_matches(self):
        """Regression test for the RPKM bug: protein_lengths keyed by the
        reference/subject ID (the pre-fix, wrong assumption) must produce no
        length data at all, not a lucky match -- confirms the lookup really
        uses qseqid, not sseqid."""
        best = {"r1": _row("r1", "ectA|P1|Org")}
        prot_lens = {"ectA|P1|Org": 200}  # keyed by sseqid, not qseqid
        lengths = estimate_family_lengths(best, prot_lens)
        assert "ectA" not in lengths
