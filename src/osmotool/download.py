"""
download.py — fetch and unpack osmo_refdb releases from Zenodo

Each osmo_refdb release is a separate Zenodo record with its own versioned
DOI (e.g. v5 -> 10.5281/zenodo.21420253, record id 21420253). RELEASES maps
release name -> record id; a new osmo_refdb release means adding an entry
here. The record's file list and checksums are fetched from Zenodo's REST
API at download time rather than hardcoding a download URL, so this stays
correct even if Zenodo's internal URL scheme changes.
"""

from __future__ import annotations

import hashlib
import json
import logging
import shutil
import tarfile
import urllib.request
from pathlib import Path

log = logging.getLogger("osmotool")

ZENODO_API = "https://zenodo.org/api/records/{record_id}"

RELEASES: dict[str, str] = {
    "v5": "21420253",
}
LATEST_RELEASE = "v5"


def _resolve_release(release: str) -> str:
    if release == "latest":
        release = LATEST_RELEASE
    if release not in RELEASES:
        known = ", ".join(sorted(RELEASES))
        raise ValueError(f"Unknown osmo_refdb release '{release}'. Known releases: {known}")
    return release


def _fetch_record_metadata(record_id: str) -> dict:
    url = ZENODO_API.format(record_id=record_id)
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            return json.load(resp)
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            f"Could not reach Zenodo ({url}): {exc}. Check your network "
            "connection, or download the release manually -- see the "
            "README's Reference database section."
        ) from exc


def _pick_archive_file(metadata: dict, record_id: str) -> dict:
    files = metadata.get("files", [])
    for f in files:
        if f.get("key", "").endswith(".tar.gz"):
            return f
    raise RuntimeError(
        f"No .tar.gz archive found in Zenodo record {record_id}'s file list "
        f"(files: {[f.get('key') for f in files]})"
    )


def _download_url_for(file_info: dict, record_id: str) -> str:
    url = file_info.get("links", {}).get("self")
    if url:
        return url
    # Fallback if Zenodo's response shape ever changes underneath us.
    return f"https://zenodo.org/records/{record_id}/files/{file_info['key']}?download=1"


def _download_file(url: str, dest: Path) -> None:
    log.info("Downloading %s -> %s", url, dest)
    try:
        with urllib.request.urlopen(url, timeout=30) as resp, open(dest, "wb") as fh:
            shutil.copyfileobj(resp, fh)
    except OSError as exc:
        dest.unlink(missing_ok=True)
        raise RuntimeError(f"Download failed ({url}): {exc}") from exc


def _verify_checksum(path: Path, checksum: str | None) -> None:
    if not checksum or ":" not in checksum:
        log.warning("No checksum reported by Zenodo for %s -- skipping verification", path.name)
        return
    algo, expected = checksum.split(":", 1)
    if algo != "md5":
        log.warning("Unsupported checksum algorithm '%s' -- skipping verification", algo)
        return
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    actual = digest.hexdigest()
    if actual != expected:
        raise RuntimeError(
            f"Checksum mismatch for {path.name}: expected md5:{expected}, "
            f"got md5:{actual}. The download may be corrupted or "
            "incomplete -- try again, or pass --force to retry from scratch."
        )
    log.info("Checksum OK (md5:%s)", actual)


def _safe_extractall(tar: tarfile.TarFile, dest: Path) -> None:
    """Extract *tar* into *dest*, refusing any member whose resolved path
    would land outside *dest* (path traversal via '../' or an absolute
    path in a corrupted or tampered archive)."""
    dest_resolved = dest.resolve()
    for member in tar.getmembers():
        member_path = (dest_resolved / member.name).resolve()
        if member_path != dest_resolved and dest_resolved not in member_path.parents:
            raise RuntimeError(
                f"Refusing to extract unsafe path outside destination: {member.name}"
            )
    tar.extractall(dest_resolved)


def download_refdb(
    release: str,
    location: str | Path,
    *,
    force: bool = False,
    keep_archive: bool = False,
) -> Path:
    """
    Download and unpack an osmo_refdb release from Zenodo into *location*.

    Parameters
    ----------
    release:
        Release name (e.g. "v5"), or "latest" for the newest known release.
    location:
        Directory to download/extract into. Created if it doesn't exist.
        The release unpacks into a subdirectory named after the release
        (e.g. ``<location>/v5/``), matching osmo_refdb's tarball layout --
        pass that path as osmotool's DATABASE argument.
    force:
        Re-download and overwrite even if ``<location>/<release>`` already
        exists.
    keep_archive:
        Keep the downloaded ``.tar.gz`` alongside the extracted directory
        instead of deleting it after a successful extraction.

    Returns
    -------
    Path to the extracted release directory.
    """
    release = _resolve_release(release)
    record_id = RELEASES[release]

    location = Path(location)
    location.mkdir(parents=True, exist_ok=True)
    out_dir = location / release

    if out_dir.exists() and not force:
        log.info("%s already exists, skipping download (use --force to re-download)", out_dir)
        return out_dir

    log.info("Fetching osmo_refdb %s metadata from Zenodo (record %s)", release, record_id)
    metadata = _fetch_record_metadata(record_id)
    file_info = _pick_archive_file(metadata, record_id)
    download_url = _download_url_for(file_info, record_id)
    archive_path = location / file_info["key"]

    _download_file(download_url, archive_path)
    _verify_checksum(archive_path, file_info.get("checksum"))

    log.info("Extracting %s -> %s", archive_path, location)
    with tarfile.open(archive_path) as tar:
        _safe_extractall(tar, location)

    if not out_dir.exists():
        raise RuntimeError(
            f"Extraction finished but expected directory {out_dir} was not "
            "created -- osmo_refdb's tarball layout may have changed."
        )

    if keep_archive:
        log.info("Archive retained: %s", archive_path)
    else:
        archive_path.unlink(missing_ok=True)

    log.info("osmo_refdb %s ready: %s", release, out_dir)
    return out_dir
