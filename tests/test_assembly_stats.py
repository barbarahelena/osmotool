"""
Tests for annotate's whole-assembly stats (genome_bp / gc / n_contigs /
contig_n50): functions.compute_assembly_stats and the gene_counts.tsv
header lines run_annotate writes from it.
"""

from __future__ import annotations

from pathlib import Path
from unittest import mock

import pytest

from osmotool.functions import AssemblyStats, compute_assembly_stats, _n50
from osmotool.quantifier import KNOWN_FAMILIES
from osmotool.runners import run_annotate


def _write(path: Path, data: bytes) -> Path:
    path.write_bytes(data)
    return path


# ---------------------------------------------------------------------------
# compute_assembly_stats
# ---------------------------------------------------------------------------

class TestComputeAssemblyStats:
    def test_multi_record(self, tmp_path: Path):
        # ctg1 = 12 bp (wrapped over two lines), ctg2 = 4 bp, ctg3 = 8 bp
        fa = _write(tmp_path / "a.fna", (
            b">ctg1 some description\n"
            b"ACGTAC\n"
            b"GTACGT\n"
            b">ctg2\n"
            b"GGGG\n"
            b">ctg3\n"
            b"AATTAATT\n"
        ))
        st = compute_assembly_stats(fa)
        assert st.genome_bp == 24
        assert st.n_contigs == 3
        # G+C: ctg1 6, ctg2 4, ctg3 0 -> 10 / 24
        assert st.gc == pytest.approx(10 / 24)
        # sorted desc 12, 8, 4; 12*2 = 24 >= 24 -> N50 = 12
        assert st.contig_n50 == 12

    def test_lowercase_counts_same_as_uppercase(self, tmp_path: Path):
        upper = compute_assembly_stats(_write(tmp_path / "u.fna", b">c\nACGTGGCC\n"))
        lower = compute_assembly_stats(_write(tmp_path / "l.fna", b">c\nacgtggcc\n"))
        mixed = compute_assembly_stats(_write(tmp_path / "m.fna", b">c\nAcGtgGcC\n"))
        assert upper == lower == mixed
        assert upper.gc == pytest.approx(6 / 8)

    def test_n_runs_count_toward_length_not_gc(self, tmp_path: Path):
        fa = _write(tmp_path / "n.fna", b">c\nGCNNNNNNNNAT\nnnnn\n")
        st = compute_assembly_stats(fa)
        assert st.genome_bp == 16          # every sequence character
        assert st.gc == pytest.approx(0.5)  # GC / (GCAT) = 2 / 4

    def test_iupac_ambiguity_codes_excluded_from_gc(self, tmp_path: Path):
        # S (=G/C) and W (=A/T) etc. must not leak into either side of the ratio
        fa = _write(tmp_path / "i.fna", b">c\nGGAARYSWKMBDHVN\n")
        st = compute_assembly_stats(fa)
        assert st.genome_bp == 15
        assert st.gc == pytest.approx(2 / 4)

    def test_empty_record(self, tmp_path: Path):
        fa = _write(tmp_path / "e.fna", (
            b">ctg1\n"
            b"GGCC\n"
            b">empty\n"
            b">ctg3\n"
            b"AATT\n"
        ))
        st = compute_assembly_stats(fa)
        assert st.n_contigs == 3   # empty record still counts
        assert st.genome_bp == 8
        assert st.gc == pytest.approx(0.5)
        assert st.contig_n50 == 4

    def test_empty_record_at_end(self, tmp_path: Path):
        st = compute_assembly_stats(_write(tmp_path / "e.fna", b">a\nACGT\n>b\n"))
        assert (st.n_contigs, st.genome_bp, st.contig_n50) == (2, 4, 4)

    def test_crlf_line_endings(self, tmp_path: Path):
        lf = compute_assembly_stats(_write(
            tmp_path / "lf.fna", b">c1\nACGT\nGG\n>c2\nAT\n"))
        crlf = compute_assembly_stats(_write(
            tmp_path / "crlf.fna", b">c1\r\nACGT\r\nGG\r\n>c2\r\nAT\r\n"))
        assert crlf == lf
        assert crlf.genome_bp == 8   # \r not counted as sequence

    def test_blank_lines_and_missing_final_newline(self, tmp_path: Path):
        st = compute_assembly_stats(_write(tmp_path / "b.fna", b">c\nACGT\n\nGG"))
        assert st.genome_bp == 6
        assert st.gc == pytest.approx(4 / 6)

    def test_no_acgt_gives_gc_none(self, tmp_path: Path):
        st = compute_assembly_stats(_write(tmp_path / "n.fna", b">c\nNNNN\n"))
        assert st.genome_bp == 4
        assert st.gc is None

    def test_empty_file(self, tmp_path: Path):
        st = compute_assembly_stats(_write(tmp_path / "empty.fna", b""))
        assert st == AssemblyStats(genome_bp=0, gc=None, n_contigs=0, contig_n50=0)

    def test_sequence_before_header_raises(self, tmp_path: Path):
        with pytest.raises(ValueError, match="before the first '>'"):
            compute_assembly_stats(_write(tmp_path / "bad.fna", b"ACGT\n>c\nACGT\n"))


class TestN50:
    def test_hand_computed(self):
        # lengths 2, 3, 4, 5, 6 -> total 20, half 10.
        # desc cumulative: 6, 11 -> 11 >= 10 at length 5.
        assert _n50([2, 3, 4, 5, 6]) == 5

    def test_exactly_half_boundary(self):
        # 10, 5, 5 -> total 20; first record alone reaches exactly half
        assert _n50([5, 10, 5]) == 10

    def test_single_record(self):
        assert _n50([1234]) == 1234

    def test_zero_length_records_ignored(self):
        assert _n50([0, 0, 7, 3]) == 7

    def test_empty(self):
        assert _n50([]) == 0
        assert _n50([0, 0]) == 0


# ---------------------------------------------------------------------------
# run_annotate end to end (Prodigal and hmmsearch faked)
# ---------------------------------------------------------------------------

# ctg1 = 18 bp, ctg2 = 6 bp, ctg3 = 2 bp -> 26 bp, N50 18.
# G+C = 5 + 4 (ctg1) + 2 (ctg2) + 2 (ctg3) = 13; A+C+G+T = 14 + 4 + 2 = 20.
ANNOTATE_ASSEMBLY = (
    b">ctg1 test contig\n"
    b"ACGTACGTAC\n"
    b"GGCCNNNN\n"
    b">ctg2\n"
    b"acgtnn\n"
    b">ctg3\n"
    b"GC\n"
)

FAKE_PROTEINS = (
    ">ctg1_1 # 1 # 303 # 1 # ID=1_1\n" + "M" * 100 + "\n"
    ">ctg1_2 # 400 # 852 # -1 # ID=1_2\n" + "M" * 150 + "\n"
    ">ctg2_1 # 1 # 603 # 1 # ID=2_1\n" + "M" * 200 + "\n"
)

# hmmsearch --tblout: target (protein), accession, query (HMM), accession,
# E-value, score, ...
FAKE_TBLOUT = (
    "# target name  accession  query name  accession  E-value  score  bias\n"
    "ctg1_1  -  ectA  -  1e-50  200.0  0.0\n"
    "ctg1_2  -  ectB  -  1e-40  150.0  0.0\n"
    "ctg1_2  -  ectC  -  1e-10   50.0  0.0\n"   # loses best-hit to ectB
    "ctg2_1  -  kdpA  -  1e-60  300.0  0.0\n"
)

# The gene table main produced for this exact input before assembly stats
# were added (captured from commit 18f6b42). Families not listed here are
# all "<fam>\t0.0\t0".
EXPECTED_HIT_ROWS = {
    "ectA": "ectA\t1.0\t1.11111e+06",
    "ectB": "ectB\t1.0\t740741",
    "kdpA": "kdpA\t1.0\t555556",
}


def _run_fake_annotate(tmp_path: Path, **kwargs) -> dict:
    assembly = _write(tmp_path / "genome.fna", ANNOTATE_ASSEMBLY)

    def fake_prodigal(assembly_path, proteins_out, **_):
        Path(proteins_out).write_text(FAKE_PROTEINS)
        return Path(proteins_out), None

    def fake_hmmsearch(hmm_db, proteins, out_file, **_):
        Path(out_file).write_text(FAKE_TBLOUT)
        return Path(out_file)

    with mock.patch("osmotool.runners.check_tool"), \
         mock.patch("osmotool.runners.run_prodigal", side_effect=fake_prodigal), \
         mock.patch("osmotool.runners.run_hmmsearch", side_effect=fake_hmmsearch):
        return run_annotate(
            "unused.dmnd", str(assembly), str(tmp_path / "out"),
            method="hmm", hmm_db="unused.hmm", **kwargs,
        )


class TestAnnotateGeneCountsHeader:
    def test_header_lines(self, tmp_path: Path):
        counts_path = _run_fake_annotate(tmp_path)["hmm"][0]
        comments = [l for l in counts_path.read_text().splitlines() if l.startswith("#")]
        assert comments == [
            "# total_proteins\t3",
            "# genome_bp\t26",
            "# gc\t0.65000",
            "# n_contigs\t3",
            "# contig_n50\t18",
            "# filtered_reads\t3",
            "# proteins_with_hits\t3",
        ]

    def test_gene_table_unchanged(self, tmp_path: Path):
        counts_path = _run_fake_annotate(tmp_path)["hmm"][0]
        table = [l for l in counts_path.read_text().splitlines() if not l.startswith("#")]
        expected = ["gene\traw_count\tCOPIES_PER_KB"] + [
            EXPECTED_HIT_ROWS.get(fam, f"{fam}\t0.0\t0") for fam in sorted(KNOWN_FAMILIES)
        ]
        assert table == expected

    def test_stats_written_with_precalled_proteins(self, tmp_path: Path):
        proteins = tmp_path / "given.faa"
        proteins.write_text(FAKE_PROTEINS)
        counts_path = _run_fake_annotate(tmp_path, proteins=str(proteins))["hmm"][0]
        assert "# genome_bp\t26\n" in counts_path.read_text()
