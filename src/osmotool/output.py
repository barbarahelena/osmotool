"""
output.py — write osmotool result tables
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from osmotool.functions import AssemblyStats
    from osmotool.systems import SystemCall

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
    copies_per_genome: dict[str, float | None] | None = None,
    assembly_stats: AssemblyStats | None = None,
    proteins_with_hits: int | None = None,
) -> Path:
    """
    Write ``<out_prefix>.gene_counts.tsv``.

    Format::

        # <total_reads_label>    <N>
        # filtered_reads <N>
        gene    raw_count   <rpm|copies_per_kb>   [copies_per_genome]
        ectA    12.5        3.21                  [0.42]
        ...

    total_reads_label:
        Header label for the first row. ``profile`` mode passes reads
        (the default, "total_reads"); ``annotate`` mode passes
        "total_proteins" since that count is Prodigal-called ORFs, not
        sequencing reads.
    copies_per_genome:
        ``profile`` mode with marker genes available (see normalise.py):
        adds a final ``copies_per_genome`` column, ``NA`` where not
        estimable. Omitted (None) leaves the file byte-identical to before.
    assembly_stats:
        ``annotate`` mode only. When given, four more comment lines follow
        the first one: ``# genome_bp``, ``# gc`` (5 decimals, or ``NA`` if
        the assembly has no A/C/G/T), ``# n_contigs`` and ``# contig_n50``.
        ``profile`` leaves this None, so its output is unchanged.
    proteins_with_hits:
        ``annotate`` mode only. When given, written as
        ``# proteins_with_hits`` directly after ``# filtered_reads``. It is
        the same number under a name that fits annotate; ``filtered_reads``
        stays for parsers that already read it.

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
        if assembly_stats is not None:
            gc = "NA" if assembly_stats.gc is None else f"{assembly_stats.gc:.5f}"
            fh.write(f"# genome_bp\t{assembly_stats.genome_bp}\n")
            fh.write(f"# gc\t{gc}\n")
            fh.write(f"# n_contigs\t{assembly_stats.n_contigs}\n")
            fh.write(f"# contig_n50\t{assembly_stats.contig_n50}\n")
        fh.write(f"# filtered_reads\t{filtered_reads}\n")
        if proteins_with_hits is not None:
            fh.write(f"# proteins_with_hits\t{proteins_with_hits}\n")
        extra_col = "\tcopies_per_genome" if copies_per_genome is not None else ""
        fh.write(f"gene\traw_count\t{norm_col}{extra_col}\n")
        for fam in families:
            raw = counts.get(fam, 0.0)
            norm = normalised.get(fam, 0.0)
            extra = ""
            if copies_per_genome is not None:
                extra = "\t" + _fmt_optional(copies_per_genome.get(fam))
            fh.write(f"{fam}\t{raw:.1f}\t{norm:.6g}{extra}\n")

    log.info("Wrote gene counts → %s", out_path)
    return out_path


def _fmt_optional(value: float | None) -> str:
    return "NA" if value is None else f"{value:.6g}"


# ---------------------------------------------------------------------------
# markers.tsv (profile mode, copies-per-genome normalisation -- see normalise.py)
# ---------------------------------------------------------------------------

def write_marker_report(
    out_prefix: str | Path,
    counts: dict[str, float],
    lengths_aa: dict[str, float],
    rpk: dict[str, float | None],
    estimate,  # normalise.MarkerEstimate
    *,
    read_length: int,
) -> Path:
    """
    Write ``<out_prefix>.markers.tsv``: per-marker evidence behind the
    genome-coverage estimate, so a copies_per_genome value can be audited.

    Format::

        # read_length         150
        # genome_coverage_rpk 12.3     (NA if not estimable)
        # markers_used        8
        # note                ...      (only when there is one)
        marker  raw_count  ref_length_aa  RPK  status
        argS    40.5       575.2          23.1 used | outlier | undetected
    """
    out_path = Path(str(out_prefix) + ".markers.tsv")
    status = {m: "used" for m in estimate.used}
    status.update({m: "outlier" for m in estimate.outliers})
    status.update({m: "undetected" for m in estimate.undetected})
    with open(out_path, "w") as fh:
        fh.write(f"# read_length\t{read_length}\n")
        fh.write(f"# genome_coverage_rpk\t{_fmt_optional(estimate.coverage_rpk)}\n")
        fh.write(f"# markers_used\t{len(estimate.used)}\n")
        if estimate.note:
            fh.write(f"# note\t{estimate.note}\n")
        fh.write("marker\traw_count\tref_length_aa\tRPK\tstatus\n")
        for m in sorted(counts):
            length = lengths_aa.get(m)
            fh.write(f"{m}\t{counts[m]:.1f}\t{'NA' if length is None else f'{length:.1f}'}\t"
                     f"{_fmt_optional(rpk.get(m))}\t{status.get(m, 'undetected')}\n")
    log.info("Wrote marker report → %s", out_path)
    return out_path


# ---------------------------------------------------------------------------
# gene_coordinates.tsv (annotate mode only -- see quantifier.extract_gene_coordinates)
# ---------------------------------------------------------------------------

def write_gene_coordinates(
    out_prefix: str | Path,
    rows: list[dict],
) -> Path:
    """
    Write ``<out_prefix>.gene_coordinates.tsv``: one row per called gene
    (best-hit protein with known Prodigal coordinates), unlike
    gene_counts.tsv which reports one row per *family* summed across all
    its hits.

    Format::

        gene    protein_id  contig  start   end     strand  gene_index
        ectA    NODE_1_7    NODE_1  1801    2400    1       7
        ectB    NODE_1_8    NODE_1  2410    3200    1       8

    Rows are sorted by (gene, contig, start) for stable, readable output.
    *rows* comes from :func:`osmotool.quantifier.extract_gene_coordinates`.
    """
    out_path = Path(str(out_prefix) + ".gene_coordinates.tsv")

    ordered = sorted(rows, key=lambda r: (r["gene"], r["contig"], r["start"]))

    with open(out_path, "w") as fh:
        fh.write("gene\tprotein_id\tcontig\tstart\tend\tstrand\tgene_index\n")
        for row in ordered:
            fh.write(
                f"{row['gene']}\t{row['protein_id']}\t{row['contig']}\t"
                f"{row['start']}\t{row['end']}\t{row['strand']}\t{row['gene_index']}\n"
            )

    log.info("Wrote gene coordinates → %s", out_path)
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

#: Task 3 -- explicit read-mode note carried into profile mode's aln_stats.tsv
#: (see run_profile's docstring in runners.py for the full rationale). Not a
#: log-only warning: this is written into the actual output file so it
#: survives past the terminal a run happened in.
PROFILE_SYSTEMS_NOTE = (
    "multi-subunit system calls (opuA/opuB/opuC/mrp/ect/ots/trk/ktr) are NOT "
    "confirmed co-located in read-profiling mode -- abundances below are per-gene "
    "and independent; a system-level 'present' claim requires annotate mode "
    "(or a dedicated assembly+binning step) on the same samples"
)


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
    note_systems: str | None = None,
    extra: dict | None = None,
) -> dict:
    """Return a flat dict of summary statistics for :func:`write_aln_stats`.

    note_systems: profile mode passes PROFILE_SYSTEMS_NOTE (Task 3); annotate
    mode leaves this None since it has its own systems.tsv output instead.
    extra: additional key/value rows inserted before the per-family counts
    (profile's marker-normalisation summary).
    """
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
    if note_systems is not None:
        stats["note_systems"] = note_systems
    if extra:
        stats.update(extra)
    # per-family breakdown
    for fam, cnt in sorted(counts.items()):
        stats[f"count_{fam}"] = f"{cnt:.1f}"
    return stats


# ---------------------------------------------------------------------------
# systems.tsv (annotate mode only -- see systems.detect_systems)
# ---------------------------------------------------------------------------

def write_system_calls(
    out_prefix: str | Path,
    calls: dict[str, "SystemCall"],
) -> Path:
    """
    Write ``<out_prefix>.systems.tsv``: one row per multi-subunit system
    (see systems.DEFAULT_SYSTEMS), reporting whether its mandatory subunits
    were found co-located within a proximity window on the same contig.

    Format::

        system  present subunit_count   colocalized fragmented_assembly_flag   contig  mandatory_hits  accessory_hits  rescued_hits
        opuA    True    3               True        False                      NODE_1  opuAA,opuAB,opuAC   -           -
        opuB    False   2               False       False                      -       opuBB,opuBC     -               opuBC<-opuCC
        mrp     False   2               False       True                       -       mrpA,mrpC       -               -

    present/colocalized/fragmented_assembly_flag: see systems.detect_systems's
    docstring for exactly what each means and how they're distinguished.

    rescued_hits: entries of the form "<family><-<original_label>" for any
    mandatory subunit that had zero hits under its own family label, but a
    hit under a documented paralog-confusion sibling (see
    systems.PARALOG_CONFUSION_CLIQUES) was positionally confirmed adjacent
    to another real subunit of the same system -- see
    systems._rescue_missing_mandatory_families. Populated whenever a rescue
    was confirmed, independent of whether the system as a whole ends up
    `present` (e.g. one subunit correctly reattributed via rescue doesn't
    help if a *different* mandatory subunit is genuinely absent elsewhere).
    """
    out_path = Path(str(out_prefix) + ".systems.tsv")

    with open(out_path, "w") as fh:
        fh.write(
            "system\tpresent\tsubunit_count\tcolocalized\t"
            "fragmented_assembly_flag\tcontig\tmandatory_hits\taccessory_hits\trescued_hits\n"
        )
        for name in sorted(calls):
            call = calls[name]
            fh.write(
                f"{call.system}\t{call.present}\t{call.subunit_count}\t"
                f"{call.colocalized}\t{call.fragmented_assembly_flag}\t"
                f"{call.contig or '-'}\t"
                f"{','.join(call.mandatory_hits) or '-'}\t"
                f"{','.join(call.accessory_hits) or '-'}\t"
                f"{','.join(call.rescued_hits) or '-'}\n"
            )

    log.info("Wrote system calls → %s", out_path)
    return out_path
