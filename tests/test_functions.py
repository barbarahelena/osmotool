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
