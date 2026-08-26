"""
test_quantifier.py — unit tests for osmotool.quantifier
"""

from __future__ import annotations

import pytest

from osmotool.functions import OrfCoord
from osmotool.quantifier import (
    select_best_hits,
    count_hits,
    compute_rpm,
    compute_rpkm,
    estimate_family_lengths,
    extract_gene_coordinates,
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
        """NOTE: "r1"/"r2" here are distinct synthetic qseqids -- this
        does NOT exercise the real-world case where both mates of a pair
        share an identical qseqid (see TestPairedQseqidCollision below for
        that; this test alone would still pass even if select_best_hits()
        silently collapsed colliding mates into one entry)."""
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
# Paired-end qseqid collision (see functions._merge_paired_with_mate_tags)
# ---------------------------------------------------------------------------
#
# Modern Illumina headers put the mate marker AFTER the first whitespace
# (e.g. "... 1:N:0:INDEX" / "... 2:N:0:INDEX"), so DIAMOND -- which
# truncates headers at the first whitespace to get qseqid -- assigns both
# mates of a fragment the SAME qseqid unless something disambiguates them
# first. select_best_hits() is a plain dict keyed by qseqid, so a
# collision doesn't average or combine the two mates' hits: it silently
# keeps only the higher-scoring one and drops the other, with no error.
# The tests below simulate DIAMOND output rows both *before* tagging (the
# bug) and *after* it (functions._merge_paired_with_mate_tags's fix, which
# appends _R1/_R2 to each read ID before the merged FASTQ ever reaches
# DIAMOND) to show the fix is exactly sufficient -- select_best_hits() and
# count_hits() need no changes of their own once qseqids are unique.

class TestPairedQseqidCollision:
    def test_colliding_qseqid_same_family_undercounts(self):
        """The bug: both mates hit the SAME family but share one qseqid
        (as they would under un-tagged modern Illumina headers) -- only
        one dict entry survives, so the pair is worth 0.5, not the
        expected 1.0. This is the case a naive fix (or no fix) would
        still pass "by luck" if it only checked totals, not per-family
        correctness -- see test_paired_weight_half's note."""
        rows = [
            _row("frag", "ectA|P1|Org", bitscore=200.0),
            _row("frag", "ectA|P2|Org", bitscore=150.0),
        ]
        best = select_best_hits(iter(rows))
        assert len(best) == 1  # the collision: one dict slot, not two
        counts = count_hits(best, paired=True)
        assert counts["ectA"] == pytest.approx(0.5)  # bug: should be 1.0

    def test_colliding_qseqid_different_family_loses_a_mate(self):
        """The bug, disagreeing-family variant: mate 1 hits ectA, mate 2
        hits ectB, but they collide on qseqid -- the loser (lower
        bitscore) vanishes entirely rather than contributing its own 0.5,
        instead of the expected 0.5/0.5 split."""
        rows = [
            _row("frag", "ectA|P1|Org", bitscore=200.0),
            _row("frag", "ectB|P2|Org", bitscore=150.0),
        ]
        best = select_best_hits(iter(rows))
        counts = count_hits(best, paired=True)
        assert counts["ectA"] == pytest.approx(0.5)
        assert counts["ectB"] == pytest.approx(0.0)  # bug: mate 2 vanished

    def test_tagged_qseqid_same_family_counts_correctly(self):
        """The fix applied: _merge_paired_with_mate_tags appends _R1/_R2
        before DIAMOND ever sees the reads, so what DIAMOND (and
        therefore select_best_hits/count_hits) actually receives are
        "frag_R1"/"frag_R2" -- distinct qseqids. Same scenario as the bug
        test above, but this is what would have caught it: expect the
        full 1.0, not 0.5."""
        rows = [
            _row("frag_R1", "ectA|P1|Org", bitscore=200.0),
            _row("frag_R2", "ectA|P2|Org", bitscore=150.0),
        ]
        best = select_best_hits(iter(rows))
        assert len(best) == 2
        counts = count_hits(best, paired=True)
        assert counts["ectA"] == pytest.approx(1.0)

    def test_tagged_qseqid_different_family_splits_correctly(self):
        """Fix applied, disagreeing-family variant: both mates now get
        their own dict entry, so mate 1's 0.5 lands on ectA and mate 2's
        0.5 lands on ectB independently -- the "independent by design"
        behavior (see runners.run_profile's docstring) working as
        intended, with no extra merge/reconciliation logic needed."""
        rows = [
            _row("frag_R1", "ectA|P1|Org", bitscore=200.0),
            _row("frag_R2", "ectB|P2|Org", bitscore=150.0),
        ]
        best = select_best_hits(iter(rows))
        counts = count_hits(best, paired=True)
        assert counts["ectA"] == pytest.approx(0.5)
        assert counts["ectB"] == pytest.approx(0.5)


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


# ---------------------------------------------------------------------------
# extract_gene_coordinates
# ---------------------------------------------------------------------------

def _coord(contig="NODE_1", gene_index=1, start=1, end=300, strand=1) -> OrfCoord:
    return OrfCoord(contig=contig, gene_index=gene_index, start=start, end=end, strand=strand)


class TestExtractGeneCoordinates:
    def test_single_hit(self):
        best = {"NODE_1_1": _row("NODE_1_1", "ectA|P1|Org")}
        coords = {"NODE_1_1": _coord(gene_index=1, start=1, end=300)}
        rows = extract_gene_coordinates(best, coords)
        assert rows == [{
            "gene": "ectA", "protein_id": "NODE_1_1", "contig": "NODE_1",
            "start": 1, "end": 300, "strand": 1, "gene_index": 1,
        }]

    def test_multiple_hits(self):
        best = {
            "NODE_1_1": _row("NODE_1_1", "ectA|P1|Org"),
            "NODE_1_2": _row("NODE_1_2", "ectB|P2|Org"),
        }
        coords = {
            "NODE_1_1": _coord(gene_index=1, start=1, end=300),
            "NODE_1_2": _coord(gene_index=2, start=310, end=600),
        }
        rows = extract_gene_coordinates(best, coords)
        assert len(rows) == 2
        genes = {row["gene"] for row in rows}
        assert genes == {"ectA", "ectB"}

    def test_hit_without_coordinates_skipped(self):
        """A best hit whose protein has no entry in orf_coords (e.g.
        non-Prodigal --proteins headers) is silently dropped, matching
        parse_orf_coordinates' own degrade-gracefully behaviour."""
        best = {"r1": _row("r1", "ectA|P1|Org")}
        rows = extract_gene_coordinates(best, {})
        assert rows == []

    def test_empty_best_hits(self):
        assert extract_gene_coordinates({}, {"NODE_1_1": _coord()}) == []
