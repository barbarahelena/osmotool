"""
output.py — write osmotool result tables
"""

from __future__ import annotations

import logging
from pathlib import Path

log = logging.getLogger("osmotool")


# ---------------------------------------------------------------------------
# gene_counts.tsv
# ---------------------------------------------------------------------------

def write_gene_counts(
    out_prefix: str | Path,
    counts: dict[str, float],
    normalised: dict[str, float],
    *,
    total_reads: int,
    filtered_reads: int,
    normalisation: str = "rpm",  # "rpm" or "copies_per_kb"
    total_reads_label: str = "total_reads",
) -> Path:
    """
    Write ``<out_prefix>.gene_counts.tsv``.

    Format::

        # <total_reads_label>    <N>
        # filtered_reads <N>
        gene    raw_count   <rpm|copies_per_kb>
        ectA    12.5        3.21
        ...

    total_reads_label:
        Header label for the first row. ``profile`` mode passes reads
        (the default, "total_reads"); ``annotate`` mode passes
        "total_proteins" since that count is Prodigal-called ORFs, not
        sequencing reads.

    Returns the output path.
    """
    out_path = Path(out_prefix).with_suffix("")
    # avoid double-suffix: strip any existing suffix only if it is not
    # a meaningful extension; just append directly
    out_path = Path(str(out_prefix) + ".gene_counts.tsv")

    norm_col = normalisation.upper()
    families = sorted(counts.keys())

    with open(out_path, "w") as fh:
        fh.write(f"# {total_reads_label}\t{total_reads}\n")
        fh.write(f"# filtered_reads\t{filtered_reads}\n")
        fh.write(f"gene\traw_count\t{norm_col}\n")
        for fam in families:
            raw = counts.get(fam, 0.0)
            norm = normalised.get(fam, 0.0)
            fh.write(f"{fam}\t{raw:.1f}\t{norm:.6g}\n")

    log.info("Wrote gene counts → %s", out_path)
    return out_path


# ---------------------------------------------------------------------------
# aln_stats.tsv
# ---------------------------------------------------------------------------

def write_aln_stats(
    out_prefix: str | Path,
    stats: dict,
) -> Path:
    """
    Write ``<out_prefix>.aln_stats.tsv``.

    *stats* is a flat dict of metric → value (strings or numbers).
    """
    out_path = Path(str(out_prefix) + ".aln_stats.tsv")

    with open(out_path, "w") as fh:
        fh.write("metric\tvalue\n")
        for key, val in stats.items():
            fh.write(f"{key}\t{val}\n")

    log.info("Wrote alignment stats → %s", out_path)
    return out_path


# ---------------------------------------------------------------------------
# Stats helpers
# ---------------------------------------------------------------------------

def build_aln_stats(
    *,
    mode: str,
    total_reads: int,
    filtered_reads: int,
    counts: dict[str, float],
    db: str,
    min_identity: float,
    min_query_cover: float,
    evalue: float,
    min_subject_cover: float = 0.0,
) -> dict:
    """Return a flat dict of summary statistics for :func:`write_aln_stats`."""
    total_assigned = sum(counts.values())
    pct_assigned = (
        (total_assigned / filtered_reads * 100) if filtered_reads > 0 else 0.0
    )
    stats = {
        "mode":                mode,
        "db":                  db,
        "total_reads":         total_reads,
        "filtered_reads":      filtered_reads,
        "reads_assigned":      f"{total_assigned:.1f}",
        "pct_reads_assigned":  f"{pct_assigned:.2f}",
        "min_identity":        min_identity,
        "min_query_cover":     min_query_cover,
        "min_subject_cover":   min_subject_cover,
        "evalue":              evalue,
    }
    # per-family breakdown
    for fam, cnt in sorted(counts.items()):
        stats[f"count_{fam}"] = f"{cnt:.1f}"
    return stats
