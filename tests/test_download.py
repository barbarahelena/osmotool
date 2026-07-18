"""
test_download.py — unit tests for osmotool.download

Network calls (_fetch_record_metadata, _download_file) are mocked
throughout -- these tests never hit Zenodo.
"""

from __future__ import annotations

import hashlib
import io
import tarfile
from pathlib import Path
from unittest import mock

import pytest

from osmotool.download import (
    RELEASES,
    LATEST_RELEASE,
    _resolve_release,
    _pick_archive_file,
    _download_url_for,
    _verify_checksum,
    download_refdb,
)


# ---------------------------------------------------------------------------
# _resolve_release
# ---------------------------------------------------------------------------

class TestResolveRelease:
    def test_latest_resolves_to_known_release(self):
        assert _resolve_release("latest") == LATEST_RELEASE
        assert LATEST_RELEASE in RELEASES

    def test_known_release_passthrough(self):
        assert _resolve_release("v5") == "v5"

    def test_unknown_release_raises(self):
        with pytest.raises(ValueError, match="Unknown osmo_refdb release"):
            _resolve_release("v99")


# ---------------------------------------------------------------------------
# _pick_archive_file / _download_url_for
# ---------------------------------------------------------------------------

class TestPickArchiveFile:
    def test_finds_tar_gz(self):
        metadata = {"files": [
            {"key": "qc_scorecard.tsv"},
            {"key": "v5.tar.gz", "links": {"self": "http://x/v5.tar.gz"}},
        ]}
        info = _pick_archive_file(metadata, "123")
        assert info["key"] == "v5.tar.gz"

    def test_no_archive_raises(self):
        metadata = {"files": [{"key": "readme.txt"}]}
        with pytest.raises(RuntimeError, match="No .tar.gz archive"):
            _pick_archive_file(metadata, "123")


class TestDownloadUrlFor:
    def test_uses_links_self_when_present(self):
        info = {"key": "v5.tar.gz", "links": {"self": "http://x/v5.tar.gz/content"}}
        assert _download_url_for(info, "123") == "http://x/v5.tar.gz/content"

    def test_falls_back_when_links_missing(self):
        info = {"key": "v5.tar.gz"}
        url = _download_url_for(info, "123")
        assert url == "https://zenodo.org/records/123/files/v5.tar.gz?download=1"


# ---------------------------------------------------------------------------
# _verify_checksum
# ---------------------------------------------------------------------------

class TestVerifyChecksum:
    def test_matching_md5_passes(self, tmp_path: Path):
        p = tmp_path / "f.tar.gz"
        p.write_bytes(b"hello world")
        digest = hashlib.md5(b"hello world").hexdigest()
        _verify_checksum(p, f"md5:{digest}")  # should not raise

    def test_mismatched_md5_raises(self, tmp_path: Path):
        p = tmp_path / "f.tar.gz"
        p.write_bytes(b"hello world")
        with pytest.raises(RuntimeError, match="Checksum mismatch"):
            _verify_checksum(p, "md5:deadbeef")

    def test_missing_checksum_skips_silently(self, tmp_path: Path):
        p = tmp_path / "f.tar.gz"
        p.write_bytes(b"hello world")
        _verify_checksum(p, None)  # should not raise


# ---------------------------------------------------------------------------
# download_refdb (network mocked)
# ---------------------------------------------------------------------------

def _make_fake_archive(path: Path, release: str = "v5", content: bytes = b"fake dmnd") -> str:
    """Write a tarball at *path* containing <release>/osmo_refdb.dmnd and
    return its md5 checksum."""
    with tarfile.open(path, "w:gz") as tar:
        info = tarfile.TarInfo(name=f"{release}/osmo_refdb.dmnd")
        info.size = len(content)
        tar.addfile(info, io.BytesIO(content))
    return hashlib.md5(path.read_bytes()).hexdigest()


class TestDownloadRefdb:
    def test_full_download_and_extract(self, tmp_path: Path):
        archive_src = tmp_path / "src.tar.gz"
        md5 = _make_fake_archive(archive_src)
        metadata = {"files": [
            {"key": "v5.tar.gz", "checksum": f"md5:{md5}", "links": {"self": "http://fake/v5.tar.gz"}}
        ]}

        location = tmp_path / "out"

        def fake_download(url, dest):
            dest.write_bytes(archive_src.read_bytes())

        with mock.patch("osmotool.download._fetch_record_metadata", return_value=metadata), \
             mock.patch("osmotool.download._download_file", side_effect=fake_download):
            out_dir = download_refdb("v5", location)

        assert out_dir == location / "v5"
        assert (out_dir / "osmo_refdb.dmnd").exists()
        # archive deleted by default
        assert not (location / "v5.tar.gz").exists()

    def test_keep_archive(self, tmp_path: Path):
        archive_src = tmp_path / "src.tar.gz"
        md5 = _make_fake_archive(archive_src)
        metadata = {"files": [
            {"key": "v5.tar.gz", "checksum": f"md5:{md5}", "links": {"self": "http://fake/v5.tar.gz"}}
        ]}
        location = tmp_path / "out"

        def fake_download(url, dest):
            dest.write_bytes(archive_src.read_bytes())

        with mock.patch("osmotool.download._fetch_record_metadata", return_value=metadata), \
             mock.patch("osmotool.download._download_file", side_effect=fake_download):
            download_refdb("v5", location, keep_archive=True)

        assert (location / "v5.tar.gz").exists()

    def test_skips_existing_without_force(self, tmp_path: Path):
        location = tmp_path / "out"
        existing = location / "v5"
        existing.mkdir(parents=True)
        (existing / "osmo_refdb.dmnd").write_bytes(b"already here")

        with mock.patch("osmotool.download._fetch_record_metadata") as fetch_mock:
            out_dir = download_refdb("v5", location)

        fetch_mock.assert_not_called()
        assert out_dir == existing

    def test_force_redownloads_existing(self, tmp_path: Path):
        archive_src = tmp_path / "src.tar.gz"
        md5 = _make_fake_archive(archive_src, content=b"new content")
        metadata = {"files": [
            {"key": "v5.tar.gz", "checksum": f"md5:{md5}", "links": {"self": "http://fake/v5.tar.gz"}}
        ]}
        location = tmp_path / "out"
        existing = location / "v5"
        existing.mkdir(parents=True)

        def fake_download(url, dest):
            dest.write_bytes(archive_src.read_bytes())

        with mock.patch("osmotool.download._fetch_record_metadata", return_value=metadata), \
             mock.patch("osmotool.download._download_file", side_effect=fake_download):
            download_refdb("v5", location, force=True)

        assert (existing / "osmo_refdb.dmnd").read_bytes() == b"new content"

    def test_checksum_mismatch_raises(self, tmp_path: Path):
        archive_src = tmp_path / "src.tar.gz"
        _make_fake_archive(archive_src)
        metadata = {"files": [
            {"key": "v5.tar.gz", "checksum": "md5:0000000000000000000000000000000",
             "links": {"self": "http://fake/v5.tar.gz"}}
        ]}
        location = tmp_path / "out"

        def fake_download(url, dest):
            dest.write_bytes(archive_src.read_bytes())

        with mock.patch("osmotool.download._fetch_record_metadata", return_value=metadata), \
             mock.patch("osmotool.download._download_file", side_effect=fake_download):
            with pytest.raises(RuntimeError, match="Checksum mismatch"):
                download_refdb("v5", location)

    def test_unknown_release_raises_before_any_network_call(self, tmp_path: Path):
        with mock.patch("osmotool.download._fetch_record_metadata") as fetch_mock:
            with pytest.raises(ValueError, match="Unknown osmo_refdb release"):
                download_refdb("v99", tmp_path / "out")
        fetch_mock.assert_not_called()


# ---------------------------------------------------------------------------
# _safe_extractall (path traversal guard)
# ---------------------------------------------------------------------------

class TestSafeExtractall:
    def test_rejects_path_traversal(self, tmp_path: Path):
        from osmotool.download import _safe_extractall

        archive_path = tmp_path / "evil.tar.gz"
        with tarfile.open(archive_path, "w:gz") as tar:
            data = b"pwned"
            info = tarfile.TarInfo(name="../../etc/evil")
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))

        dest = tmp_path / "dest"
        dest.mkdir()
        with tarfile.open(archive_path) as tar:
            with pytest.raises(RuntimeError, match="unsafe path"):
                _safe_extractall(tar, dest)
