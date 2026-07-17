"""
conftest.py — shared pytest fixtures for osmotool tests
"""

from __future__ import annotations

import gzip
import textwrap
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Minimal synthetic DIAMOND alignment output
# ---------------------------------------------------------------------------
#
# Six hits (one per osmoadaptation gene family), one read each.
# Format: qseqid sseqid pident qcovhsp length evalue bitscore
#
SAMPLE_ALN_LINES = [
    # read hits ectA
    "read_1\tectA|P12345|Halomonas_elongata\t95.0\t90.0\t150\t1e-30\t250.0",
    # read hits ectB
    "read_2\tectB|P23456|Chromohalobacter_salexigens\t88.0\t85.0\t140\t1e-25\t210.0",
    # read hits ectC (two reads, best selected by bitscore)
    "read_3\tectC|P34567|Virgibacillus_halodenitrificans\t92.0\t88.0\t145\t1e-28\t235.0",
    "read_3\tectC|P34568|Some_organism\t75.0\t80.0\t130\t1e-15\t180.0",  # worse hit
    # read hits betL
    "read_4\tbetL|A0A001|Bacillus_subtilis\t85.0\t82.0\t138\t1e-20\t200.0",
    # read hits kdpA
    "read_5\tkdpA|P56789|Escherichia_coli\t90.0\t87.0\t155\t1e-32\t260.0",
    # read hits nhaA
    "read_6\tnhaA|P67890|Escherichia_coli\t93.0\t91.0\t160\t1e-35\t270.0",
]


@pytest.fixture
def aln_file(tmp_path: Path) -> Path:
    """Write a synthetic DIAMOND tabular file and return its path."""
    p = tmp_path / "test.tsv"
    p.write_text("\n".join(SAMPLE_ALN_LINES) + "\n")
    return p


@pytest.fixture
def sample_counts() -> dict[str, float]:
    """Raw counts for one read per family (single-end)."""
    return {
        "ectA": 1.0,
        "ectB": 1.0,
        "ectC": 1.0,
        "betL": 1.0,
        "kdpA": 1.0,
        "nhaA": 1.0,
    }


@pytest.fixture
def sample_rpm(sample_counts) -> dict[str, float]:
    """Pre-computed RPM for sample_counts with total_reads=6."""
    total = 6
    return {fam: (cnt * 1_000_000) / total for fam, cnt in sample_counts.items()}


@pytest.fixture
def sample_protein_fasta(tmp_path: Path) -> Path:
    """Minimal Prodigal-style protein FASTA with known lengths."""
    content = textwrap.dedent("""\
        >NODE_1_1 # 1 # 450 # 1 # ID=1
        MAAAKLLVLSAAFAATAAQA
        MAAAKLLVLSAAFAATAAQAMAAAKLLVLSAAFAATAAQAMAAAKLLVLSAAFAATAAQAMAAAKLLVLSAAFAATAAQAEND
        >NODE_1_2 # 451 # 900 # -1 # ID=2
        MAAAKLLVLSAAFAATAAQA
        MAAAKLLVLSAAFAATAAQAMAAAKLLVLSAAFAATAAQAMAAAKLLVLSAAFAATAAQAEND
    """)
    p = tmp_path / "proteins.faa"
    p.write_text(content)
    return p
