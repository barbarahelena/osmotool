"""
utils.py — shared helpers for osmotool
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

log = logging.getLogger("osmotool")


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    logging.getLogger("osmotool").setLevel(level)
    logging.getLogger("osmotool").addHandler(handler)


# ---------------------------------------------------------------------------
# External tool helpers
# ---------------------------------------------------------------------------

def check_tool(name: str) -> None:
    """Raise RuntimeError if *name* is not on PATH.

    Uses :func:`shutil.which` so that tools loaded via ``module load``
    or ``conda activate`` are found without spawning a subshell.
    """
    if shutil.which(name) is None:
        raise RuntimeError(
            f"Required tool '{name}' not found on PATH. "
            "Load the appropriate module or activate the conda environment."
        )


def run_command(cmd: list[str], *, desc: str = "", check: bool = True) -> subprocess.CompletedProcess:
    """Run *cmd*, stream stderr, raise on non-zero exit if check=True."""
    log.debug("Running: %s", " ".join(cmd))
    if desc:
        log.info(desc)
    result = subprocess.run(cmd, capture_output=False)
    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed (exit {result.returncode}): {' '.join(cmd)}")
    return result


# ---------------------------------------------------------------------------
# FASTA / path helpers
# ---------------------------------------------------------------------------

def gene_family_from_header(subject_id: str) -> str:
    """
    Extract the gene-family label from an osmodiamond subject ID.

    osmodiamond header format:  familyLabel|UniProtID|OrganismName
    e.g.  ectA|P0A393|Halomonas_elongata

    Returns the first pipe-delimited field, or the full string if no '|'.
    """
    return subject_id.split("|")[0]


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def resolve_out_prefix(out_prefix: str) -> Path:
    """Create parent directories for *out_prefix* and return a Path."""
    p = Path(out_prefix)
    ensure_dir(p.parent)
    return p


def get_tmpdir(override: str | Path | None = None) -> Path:
    """
    Return a writable temporary directory suitable for HPC scratch.

    Resolution order:
      1. *override* (explicit ``--tmpdir`` CLI argument)
      2. ``$TMPDIR`` environment variable  (set by SLURM/PBS on compute nodes)
      3. ``$SCRATCH``                      (common on some HPC systems)
      4. ``/tmp``                           (last-resort fallback)

    On most HPC schedulers the compute node sets ``$TMPDIR`` to fast
    local NVMe or SSD scratch — always prefer it over a shared filesystem.
    """
    import os
    if override is not None:
        p = Path(override)
        p.mkdir(parents=True, exist_ok=True)
        return p
    for env_var in ("TMPDIR", "SCRATCH"):
        val = os.environ.get(env_var)
        if val:
            return Path(val)
    return Path("/tmp")


# ---------------------------------------------------------------------------
# osmo_refdb release directory
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class RefDb:
    dmnd: Path
    hmm: Path | None
    cascade_config: Path | None
    diamond_cutoffs: Path | None
    exclude_families: Path | None


def resolve_refdb(database: str, mode: str) -> RefDb:
    """Resolve an osmo_refdb release directory into its fixed-name files.

    An unpacked osmo_refdb release (e.g. ``releases/v5/``) always uses the
    same filenames regardless of release version:
    ``osmo_refdb.dmnd``, ``hmms/osmo_refdb.hmm`` (+ pressed ``.h3*``
    indices), ``osmo_refdb.profile_cascade.tsv``,
    ``osmo_refdb.diamond_cutoffs.tsv``,
    ``osmo_refdb.profile_excluded_families.txt``, and
    ``osmo_refdb.annotate_excluded_families.txt``. This lets callers pass
    one directory instead of wiring up each file individually.

    *mode* is ``"profile"`` or ``"annotate"`` -- the only difference is
    which excluded-families list and cascade/cutoffs file apply. The
    DIAMOND database is required; every other file is optional and
    resolves to ``None`` when absent from *database*.
    """
    root = Path(database)
    dmnd = root / "osmo_refdb.dmnd"
    if not dmnd.exists():
        raise FileNotFoundError(
            f"DIAMOND database not found: {dmnd} -- DATABASE should be an "
            "unpacked osmo_refdb release directory containing osmo_refdb.dmnd"
        )

    hmm = root / "hmms" / "osmo_refdb.hmm"
    if not hmm.exists():
        hmm = None

    if mode == "profile":
        cascade_config = root / "osmo_refdb.profile_cascade.tsv"
        exclude_families = root / "osmo_refdb.profile_excluded_families.txt"
        diamond_cutoffs = None
    elif mode == "annotate":
        cascade_config = None
        diamond_cutoffs = root / "osmo_refdb.diamond_cutoffs.tsv"
        exclude_families = root / "osmo_refdb.annotate_excluded_families.txt"
    else:
        raise ValueError(f"Unknown mode '{mode}': expected 'profile' or 'annotate'")

    if cascade_config is not None and not cascade_config.exists():
        cascade_config = None
    if diamond_cutoffs is not None and not diamond_cutoffs.exists():
        diamond_cutoffs = None
    if exclude_families is not None and not exclude_families.exists():
        exclude_families = None

    return RefDb(
        dmnd=dmnd, hmm=hmm, cascade_config=cascade_config,
        diamond_cutoffs=diamond_cutoffs, exclude_families=exclude_families,
    )
