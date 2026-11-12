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
# _download_file: resume after an early close (urlopen mocked)
# ---------------------------------------------------------------------------

class _FakeResponse(io.BytesIO):
    def __init__(self, body: bytes, status: int, content_length: int):
        super().__init__(body)
        self.status = status
        self.headers = {"Content-Length": str(content_length)}


def _fake_server(payload: bytes, cut_at: list[int], honour_range: bool = True):
    """urlopen stand-in serving *payload*; the n-th request's connection
    closes after cut_at[n] bytes of its body (if given)."""
    calls: list[str | None] = []

    def urlopen(req, timeout=None):
        rng = req.get_header("Range")
        calls.append(rng)
        start = int(rng.split("=")[1].rstrip("-")) if (rng and honour_range) else 0
        body = payload[start:]
        n = len(calls) - 1
        if n < len(cut_at):
            body = body[:cut_at[n]]
        status = 206 if (rng and honour_range) else 200
        return _FakeResponse(body, status, len(payload) - start)

    return urlopen, calls


class TestDownloadFile:
    PAYLOAD = bytes(range(256)) * 40  # 10240 bytes

    def test_resumes_after_early_close(self, tmp_path: Path):
        from osmotool.download import _download_file
        urlopen, calls = _fake_server(self.PAYLOAD, cut_at=[3000, 2000])
        dest = tmp_path / "a.tar.gz"
        with mock.patch("urllib.request.urlopen", side_effect=urlopen):
            _download_file("http://fake/a", dest)
        assert dest.read_bytes() == self.PAYLOAD
        assert calls == [None, "bytes=3000-", "bytes=5000-"]

    def test_restarts_when_range_ignored(self, tmp_path: Path):
        from osmotool.download import _download_file
        urlopen, calls = _fake_server(self.PAYLOAD, cut_at=[3000], honour_range=False)
        dest = tmp_path / "a.tar.gz"
        with mock.patch("urllib.request.urlopen", side_effect=urlopen):
            _download_file("http://fake/a", dest)
        assert dest.read_bytes() == self.PAYLOAD
        assert len(calls) == 2

    def test_gives_up_after_repeated_early_closes(self, tmp_path: Path):
        from osmotool.download import _download_file, DOWNLOAD_ATTEMPTS
        urlopen, _ = _fake_server(self.PAYLOAD, cut_at=[100] * DOWNLOAD_ATTEMPTS)
        dest = tmp_path / "a.tar.gz"
        with mock.patch("urllib.request.urlopen", side_effect=urlopen):
            with pytest.raises(RuntimeError, match="kept closing early"):
                _download_file("http://fake/a", dest)
        assert not dest.exists()


# ---------------------------------------------------------------------------
# download_refdb (network mocked)
# ---------------------------------------------------------------------------

def _make_fake_archive(path: Path, release: str = "v5", content: bytes = b"fake dmnd") -> str:
    """Write a tarball at *path* containing <release>/osmo_refdb.dmnd and
    return its md5 checksum. *release* is the archive's top-level directory
    name (``osmo_refdb_v9`` for v9's layout)."""
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

    def _v9_metadata(self, md5: str) -> dict:
        # v9's Zenodo record: archive named osmo_refdb_v9.tar.gz, listed
        # after a .sha256 sidecar that must not be picked.
        return {"files": [
            {"key": "osmo_refdb_v9.tar.gz.sha256", "checksum": "md5:x",
             "links": {"self": "http://fake/osmo_refdb_v9.tar.gz.sha256"}},
            {"key": "osmo_refdb_v9.tar.gz", "checksum": f"md5:{md5}",
             "links": {"self": "http://fake/osmo_refdb_v9.tar.gz"}},
        ]}

    def test_archive_root_renamed_to_release(self, tmp_path: Path):
        archive_src = tmp_path / "src.tar.gz"
        md5 = _make_fake_archive(archive_src, release="osmo_refdb_v9")
        location = tmp_path / "out"

        def fake_download(url, dest):
            dest.write_bytes(archive_src.read_bytes())

        with mock.patch("osmotool.download._fetch_record_metadata",
                        return_value=self._v9_metadata(md5)), \
             mock.patch("osmotool.download._download_file", side_effect=fake_download):
            out_dir = download_refdb("v9", location)

        assert out_dir == location / "v9"
        assert (out_dir / "osmo_refdb.dmnd").exists()
        assert not (location / "osmo_refdb_v9").exists()
        assert not (location / "osmo_refdb_v9.tar.gz").exists()

    def test_force_with_renamed_root_replaces_existing(self, tmp_path: Path):
        archive_src = tmp_path / "src.tar.gz"
        md5 = _make_fake_archive(archive_src, release="osmo_refdb_v9", content=b"new content")
        location = tmp_path / "out"
        existing = location / "v9"
        existing.mkdir(parents=True)
        (existing / "stale.txt").write_text("old")

        def fake_download(url, dest):
            dest.write_bytes(archive_src.read_bytes())

        with mock.patch("osmotool.download._fetch_record_metadata",
                        return_value=self._v9_metadata(md5)), \
             mock.patch("osmotool.download._download_file", side_effect=fake_download):
            download_refdb("v9", location, force=True)

        assert (existing / "osmo_refdb.dmnd").read_bytes() == b"new content"
        assert not (existing / "stale.txt").exists()

    def test_multiple_top_level_entries_raise(self, tmp_path: Path):
        archive_src = tmp_path / "src.tar.gz"
        with tarfile.open(archive_src, "w:gz") as tar:
            for name in ("a/osmo_refdb.dmnd", "b/osmo_refdb.dmnd"):
                info = tarfile.TarInfo(name=name)
                info.size = 1
                tar.addfile(info, io.BytesIO(b"x"))
        md5 = hashlib.md5(archive_src.read_bytes()).hexdigest()
        metadata = {"files": [
            {"key": "v5.tar.gz", "checksum": f"md5:{md5}", "links": {"self": "http://fake/v5.tar.gz"}}
        ]}

        def fake_download(url, dest):
            dest.write_bytes(archive_src.read_bytes())

        with mock.patch("osmotool.download._fetch_record_metadata", return_value=metadata), \
             mock.patch("osmotool.download._download_file", side_effect=fake_download):
            with pytest.raises(RuntimeError, match="single top-level directory"):
                download_refdb("v5", tmp_path / "out")

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
