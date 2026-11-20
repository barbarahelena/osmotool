"""
test_functions.py — unit tests for osmotool.functions

Focus: _merge_paired_with_mate_tags / strip_mate_tag, the fix for DIAMOND
qseqid collisions between paired-end mates (see runners.run_profile's
paired-end handling). Modern Illumina headers put the mate marker AFTER
the first whitespace (e.g. "... 1:N:0:INDEX" / "... 2:N:0:INDEX"), so
without tagging, R1 and R2 records for the same fragment produce an
identical qseqid under DIAMOND (which truncates headers at the first
whitespace) -- see test_quantifier.py's TestPairedQseqidCollision for
what that does downstream to select_best_hits()/count_hits().
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from osmotool.functions import _merge_paired_with_mate_tags, strip_mate_tag


def _read_lines(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return fh.readlines()


def _first_token(header_line: str) -> str:
    """First whitespace-delimited token of a FASTQ header, minus '@' --
    i.e. exactly what DIAMOND would use as qseqid."""
    return header_line[1:].rstrip("\n").split()[0]


class TestMergePairedWithMateTags:
    def test_modern_illumina_headers_collide_without_tag(self, tmp_path: Path):
        """Sanity check on the premise of the fix: confirm the realistic
        collision actually exists in the raw inputs before tagging --
        both mates' headers share the same first-whitespace token."""
        r1_header = "@M00123:45:000000000-A1B2C:1:1101:15000:1000 1:N:0:ATCGAT"
        r2_header = "@M00123:45:000000000-A1B2C:1:1101:15000:1000 2:N:0:ATCGAT"
        assert _first_token(r1_header) == _first_token(r2_header)

    def test_tags_disambiguate_colliding_headers(self, tmp_path: Path):
        r1 = tmp_path / "r1.fastq"
        r2 = tmp_path / "r2.fastq"
        r1.write_text(
            "@M00123:45:000000000-A1B2C:1:1101:15000:1000 1:N:0:ATCGAT\n"
            "ACGT\n+\nIIII\n"
        )
        r2.write_text(
            "@M00123:45:000000000-A1B2C:1:1101:15000:1000 2:N:0:ATCGAT\n"
            "TGCA\n+\nIIII\n"
        )
        dest = tmp_path / "merged.fastq.gz"
        _merge_paired_with_mate_tags(r1, r2, dest)

        lines = _read_lines(dest)
        assert len(lines) == 8  # 2 records x 4 lines
        header1, header2 = lines[0], lines[4]
        token1, token2 = _first_token(header1), _first_token(header2)

        assert token1 != token2, "mate qseqids must be unique after tagging"
        assert token1 == "M00123:45:000000000-A1B2C:1:1101:15000:1000_R1"
        assert token2 == "M00123:45:000000000-A1B2C:1:1101:15000:1000_R2"
        # descriptor after the original whitespace is preserved
        assert header1.rstrip("\n").endswith("1:N:0:ATCGAT")
        assert header2.rstrip("\n").endswith("2:N:0:ATCGAT")

    def test_sequence_and_quality_lines_untouched(self, tmp_path: Path):
        r1 = tmp_path / "r1.fastq"
        r2 = tmp_path / "r2.fastq"
        r1.write_text("@frag 1:N:0:X\nACGTACGT\n+\nIIIIIIII\n")
        r2.write_text("@frag 2:N:0:X\nTTTTGGGG\n+\nJJJJJJJJ\n")
        dest = tmp_path / "merged.fastq.gz"
        _merge_paired_with_mate_tags(r1, r2, dest)

        lines = _read_lines(dest)
        assert lines[1] == "ACGTACGT\n"
        assert lines[2] == "+\n"
        assert lines[3] == "IIIIIIII\n"
        assert lines[5] == "TTTTGGGG\n"
        assert lines[7] == "JJJJJJJJ\n"

    def test_old_style_already_unique_headers_still_get_tagged(self, tmp_path: Path):
        """Older READID/1 - READID/2 headers never collided in the first
        place, but the fix tags unconditionally rather than relying on
        header convention -- confirms that doesn't break anything."""
        r1 = tmp_path / "r1.fastq"
        r2 = tmp_path / "r2.fastq"
        r1.write_text("@READID/1\nACGT\n+\nIIII\n")
        r2.write_text("@READID/2\nTGCA\n+\nIIII\n")
        dest = tmp_path / "merged.fastq.gz"
        _merge_paired_with_mate_tags(r1, r2, dest)

        lines = _read_lines(dest)
        token1, token2 = _first_token(lines[0]), _first_token(lines[4])
        assert token1 == "READID/1_R1"
        assert token2 == "READID/2_R2"

    def test_handles_gzipped_input(self, tmp_path: Path):
        r1 = tmp_path / "r1.fastq.gz"
        r2 = tmp_path / "r2.fastq.gz"
        with gzip.open(r1, "wt") as fh:
            fh.write("@frag 1:N:0:X\nACGT\n+\nIIII\n")
        with gzip.open(r2, "wt") as fh:
            fh.write("@frag 2:N:0:X\nTGCA\n+\nIIII\n")
        dest = tmp_path / "merged.fastq.gz"
        _merge_paired_with_mate_tags(r1, r2, dest)

        lines = _read_lines(dest)
        assert _first_token(lines[0]) == "frag_R1"
        assert _first_token(lines[4]) == "frag_R2"

    def test_malformed_header_raises(self, tmp_path: Path):
        r1 = tmp_path / "r1.fastq"
        r2 = tmp_path / "r2.fastq"
        r1.write_text("NOT_A_HEADER\nACGT\n+\nIIII\n")
        r2.write_text("@frag 2:N:0:X\nTGCA\n+\nIIII\n")
        dest = tmp_path / "merged.fastq.gz"
        with pytest.raises(ValueError, match="Malformed FASTQ header"):
            _merge_paired_with_mate_tags(r1, r2, dest)


class TestStripMateTag:
    def test_strips_r1(self):
        assert strip_mate_tag("frag_R1") == "frag"

    def test_strips_r2(self):
        assert strip_mate_tag("frag_R2") == "frag"

    def test_unchanged_when_untagged(self):
        """Single-end qseqids are never tagged -- must round-trip as-is."""
        assert strip_mate_tag("frag") == "frag"

    def test_unchanged_for_unrelated_trailing_text(self):
        assert strip_mate_tag("frag_R12") == "frag_R12"


# ---------------------------------------------------------------------------
# hmmsearch tblout parsing
# ---------------------------------------------------------------------------
#
# hmmsearch's tblout puts the *sequence* in column 0 and the *HMM* in
# column 2 -- the reverse of hmmscan. A swap here would silently turn every
# protein ID into a "family" and vice versa, so pin the orientation down.

from osmotool.functions import count_hmm_models, parse_hmmsearch_output

_HMMSEARCH_TBLOUT = """\
#                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
NC_000913.3_890      -          opuCA                -              8e-50  160.9   0.1   1.1e-49  160.5   0.1   1.0   1   0   0   1   1   1   1 # 1 # 1200 # 1 # ID=1_890
NC_000913.3_12       -          betL                 -            3.2e-07   25.4   2.0   4.1e-07   25.0   2.0   1.1   1   0   0   1   1   1   1 # 5 # 900 # -1 # ID=1_12
#
# Program:         hmmsearch
"""


class TestParseHmmsearchOutput:
    def test_sequence_is_qseqid_and_hmm_is_sseqid(self, tmp_path):
        tbl = tmp_path / "x.tblout"
        tbl.write_text(_HMMSEARCH_TBLOUT)
        rows = list(parse_hmmsearch_output(tbl))
        assert rows == [
            {"qseqid": "NC_000913.3_890", "sseqid": "opuCA", "evalue": 8e-50, "bitscore": 160.9},
            {"qseqid": "NC_000913.3_12", "sseqid": "betL", "evalue": 3.2e-07, "bitscore": 25.4},
        ]

    def test_empty_tblout(self, tmp_path):
        tbl = tmp_path / "x.tblout"
        tbl.write_text("# only comments\n#\n")
        assert list(parse_hmmsearch_output(tbl)) == []


class TestCountHmmModels:
    def test_counts_name_lines(self, tmp_path):
        hmm = tmp_path / "db.hmm"
        hmm.write_text(
            "HMMER3/f [3.4 | Aug 2023]\nNAME  betL\nLENG  500\n//\n"
            "HMMER3/f [3.4 | Aug 2023]\nNAME  opuCA\nLENG  380\n//\n"
        )
        assert count_hmm_models(hmm) == 2

    def test_rejects_non_hmm_file(self, tmp_path):
        bad = tmp_path / "db.hmm"
        bad.write_text(">seq\nMKV\n")
        with pytest.raises(ValueError):
            count_hmm_models(bad)


# ---------------------------------------------------------------------------
# strip_orf_suffix (orfm ORF ID → original read ID, used by the cascade)
# ---------------------------------------------------------------------------

from osmotool.functions import strip_orf_suffix


class TestStripOrfSuffix:
    @pytest.mark.parametrize("orf_id, read_id", [
        # legacy /1 /2 mate markers
        ("read42/1_1_2_3", "read42/1"),
        ("read42/2_10_4_1", "read42/2"),
        # no mate marker: single-end reads, and modern Illumina headers
        # whose mate marker sits after the whitespace orfm drops
        ("r215486_3_3_4", "r215486"),
        ("A00123:8:H5:1:1101:1000:2000_1_1_4", "A00123:8:H5:1:1101:1000:2000"),
        # underscores inside the read ID itself are preserved
        ("SRR123_456.7_2_5_1", "SRR123_456.7"),
    ])
    def test_recovers_read_id(self, orf_id, read_id):
        assert strip_orf_suffix(orf_id) == read_id

    def test_non_orfm_id_passes_through(self):
        assert strip_orf_suffix("not_an_orf") == "not_an_orf"


# ---------------------------------------------------------------------------
# orfm minimum ORF length for short reads (cascade crash on 73-75 bp reads)
# ---------------------------------------------------------------------------

import shutil
from unittest import mock

from osmotool.functions import (
    ORFM_DEFAULT_MIN_LEN,
    fasta_has_sequences,
    orfm_min_length,
    run_orfm,
)


class TestOrfmMinLength:
    @pytest.mark.parametrize("read_length, expected", [
        (150, 96),   # long reads keep orfm's default
        (98, 96),    # shortest read that still fits 96 nt in every frame
        (97, 93),
        (75, 72),
        (73, 69),
        (4, 3),      # floored at a valid multiple of 3
        (1, 3),
    ])
    def test_values(self, read_length, expected):
        got = orfm_min_length(read_length)
        assert got == expected
        assert got % 3 == 0
        assert got <= ORFM_DEFAULT_MIN_LEN

    @pytest.mark.parametrize("bad", [0, -3, 70])
    def test_run_orfm_rejects_invalid_min_len(self, tmp_path, bad):
        with pytest.raises(ValueError):
            run_orfm(tmp_path / "in.fastq", tmp_path / "out.faa", min_orf_len=bad)

    def test_run_orfm_passes_min_len(self, tmp_path):
        with mock.patch("osmotool.functions.subprocess.run") as run:
            run.return_value.returncode = 0
            run_orfm(tmp_path / "in.fastq", tmp_path / "out.faa", min_orf_len=69)
        cmd = run.call_args.args[0]
        assert cmd[:3] == ["orfm", "-m", "69"]

    @pytest.mark.skipif(shutil.which("orfm") is None, reason="orfm not on PATH")
    def test_short_reads_yield_orfs_with_scaled_min_len(self, tmp_path):
        # 75 bp of stop-free sequence in all frames ("GCT" repeats)
        fq = tmp_path / "short.fastq"
        seq = ("GCT" * 25)
        fq.write_text(f"@r1\n{seq}\n+\n{'I' * len(seq)}\n")
        default_out = run_orfm(fq, tmp_path / "default.faa")
        assert not fasta_has_sequences(default_out)
        scaled = run_orfm(fq, tmp_path / "scaled.faa",
                          min_orf_len=orfm_min_length(len(seq)))
        assert fasta_has_sequences(scaled)


class TestFastaHasSequences:
    def test_empty(self, tmp_path):
        p = tmp_path / "e.faa"
        p.write_text("")
        assert not fasta_has_sequences(p)

    def test_with_record(self, tmp_path):
        p = tmp_path / "r.faa"
        p.write_text(">orf1\nMKV\n")
        assert fasta_has_sequences(p)


# ---------------------------------------------------------------------------
# count_fragments / fragments_from_reads: profile's RPM denominator
# ---------------------------------------------------------------------------

from osmotool.functions import count_fragments, fragments_from_reads
from osmotool.quantifier import compute_rpm, count_hits


def _write_fastq(path: Path, ids: list[str]) -> Path:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt") as fh:
        for rid in ids:
            fh.write(f"@{rid}\nACGT\n+\nIIII\n")
    return path


class TestCountFragments:
    def test_paired_counts_pairs(self, tmp_path):
        r1 = _write_fastq(tmp_path / "R1.fastq.gz", [f"r{i}" for i in range(6)])
        r2 = _write_fastq(tmp_path / "R2.fastq.gz", [f"r{i}" for i in range(6)])
        assert count_fragments(r1, r2, None) == 6

    def test_single_end_counts_reads(self, tmp_path):
        se = _write_fastq(tmp_path / "se.fastq", [f"r{i}" for i in range(5)])
        assert count_fragments(None, None, se) == 5
        assert count_fragments(se, None, None) == 5  # -1 without -2

    def test_out_of_sync_pairs_warn_and_halve_the_sum(self, tmp_path, caplog):
        r1 = _write_fastq(tmp_path / "R1.fastq", [f"r{i}" for i in range(6)])
        r2 = _write_fastq(tmp_path / "R2.fastq", [f"r{i}" for i in range(5)])
        assert count_fragments(r1, r2, None) == 5.5
        assert "out of sync" in caplog.text

    def test_no_input_raises(self):
        with pytest.raises(ValueError):
            count_fragments(None, None, None)

    @pytest.mark.parametrize("reads, paired, expected", [
        (1000, True, 500), (1001, True, 500.5), (1000, False, 1000),
    ])
    def test_fragments_from_reads(self, reads, paired, expected):
        assert fragments_from_reads(reads, paired=paired) == expected

    def test_paired_and_single_end_rpm_agree(self, tmp_path):
        """One library, run paired-end and as R1-only single-end, gives the
        same RPM: counts (0.5 per mate) and denominator (pairs) are both
        fragments. Before this fix paired-end RPM was half."""
        ids = [f"r{i}" for i in range(10)]
        r1 = _write_fastq(tmp_path / "R1.fastq", ids)
        r2 = _write_fastq(tmp_path / "R2.fastq", ids)
        # 4 R1 reads and 4 R2 reads hit ectA: same per-read hit rate in both mates
        pe_hits = {f"{rid}{tag}": {"sseqid": "ectA|ref"}
                   for rid in ids[:4] for tag in ("_R1", "_R2")}
        se_hits = {rid: {"sseqid": "ectA|ref"} for rid in ids[:4]}
        pe = compute_rpm(count_hits(pe_hits, paired=True), count_fragments(r1, r2, None))
        se = compute_rpm(count_hits(se_hits, paired=False), count_fragments(r1, None, None))
        assert pe["ectA"] == pytest.approx(se["ectA"]) == pytest.approx(400_000)
