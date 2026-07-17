"""
utils.py — shared helpers for osmotool
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import sys
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
