"""Read single files out of a large zip archive on a web server, without
downloading the whole archive.

The ZTF SN Ia DR2 release is one 1.37 GB zip. Most of it is spectra and, by
accident, the release's own git history (850 MB). The light curves are 125 MB of
it. A zip file keeps its table of contents (the "central directory") at the very
end, and the server honours HTTP Range requests, so we can read the table of
contents first and then fetch only the byte ranges of the members we want.

Members stored next to each other are fetched in one request (up to
``chunk_bytes``), which keeps the request count small: ~20 requests for all 3,592
light curves instead of 3,592.
"""

from __future__ import annotations

import io
import struct
import time
import zipfile
import zlib
from dataclasses import dataclass

import requests

USER_AGENT = "astrophysics-atlas-D6/1.0 (python-requests; hobbyist research)"


class HTTPRangeFile(io.RawIOBase):
    """A read-only, seekable file object backed by HTTP Range requests."""

    def __init__(self, url: str, session: requests.Session | None = None,
                 retries: int = 5, timeout: float = 120.0):
        self.url = url
        self.session = session or requests.Session()
        self.session.headers.setdefault("User-Agent", USER_AGENT)
        self.retries = retries
        self.timeout = timeout
        self.pos = 0
        self.n_requests = 0
        r = self._request("HEAD")
        self.size = int(r.headers["Content-Length"])
        self.headers = dict(r.headers)
        if r.headers.get("Accept-Ranges", "").lower() != "bytes":
            # nginx only advertises it on 206 responses; check with a tiny GET.
            probe = self._request("GET", headers={"Range": "bytes=0-15"})
            if probe.status_code != 206:
                raise RuntimeError(f"{url} does not support byte ranges")

    def _request(self, method: str, headers: dict | None = None) -> requests.Response:
        delay = 2.0
        for attempt in range(self.retries):
            try:
                r = self.session.request(method, self.url, headers=headers,
                                         timeout=self.timeout, allow_redirects=True)
                r.raise_for_status()
                self.n_requests += 1
                return r
            except (requests.ConnectionError, requests.Timeout, requests.HTTPError):
                if attempt == self.retries - 1:
                    raise
                time.sleep(delay)
                delay *= 2
        raise AssertionError("unreachable")

    def get_range(self, start: int, stop: int) -> bytes:
        """Bytes [start, stop) of the remote file."""
        stop = min(stop, self.size)
        if stop <= start:
            return b""
        r = self._request("GET", headers={"Range": f"bytes={start}-{stop - 1}"})
        if r.status_code != 206 or len(r.content) != stop - start:
            raise RuntimeError(f"bad range response: HTTP {r.status_code}, "
                               f"{len(r.content)} of {stop - start} bytes")
        return r.content

    # --- io.RawIOBase interface, so zipfile can read the central directory ---
    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.pos

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        base = {io.SEEK_SET: 0, io.SEEK_CUR: self.pos, io.SEEK_END: self.size}[whence]
        self.pos = base + offset
        return self.pos

    def readinto(self, buffer) -> int:
        data = self.get_range(self.pos, self.pos + len(buffer))
        buffer[: len(data)] = data
        self.pos += len(data)
        return len(data)


@dataclass
class RemoteZip:
    """The table of contents of a remote zip, plus the means to fetch members."""

    url: str
    remote: HTTPRangeFile
    infos: dict[str, zipfile.ZipInfo]

    @classmethod
    def open(cls, url: str) -> "RemoteZip":
        remote = HTTPRangeFile(url)
        # A 1 MB read buffer pulls the whole central directory in a few requests.
        with zipfile.ZipFile(io.BufferedReader(remote, buffer_size=1 << 20)) as z:
            infos = {i.filename: i for i in z.infolist()}
        return cls(url=url, remote=remote, infos=infos)

    def names(self, prefix: str = "", suffix: str = "") -> list[str]:
        return sorted(n for n, i in self.infos.items()
                      if n.startswith(prefix) and n.endswith(suffix) and not i.is_dir())

    def fetch(self, names, chunk_bytes: int = 8 << 20, progress=None):
        """Yield (name, uncompressed bytes) for each requested member.

        Members are grouped into runs of neighbouring byte ranges and each run
        is one HTTP request. Every member is checked against the CRC-32 stored
        in the archive's table of contents.
        """
        infos = sorted((self.infos[n] for n in names), key=lambda i: i.header_offset)
        # Upper bound on a local header: 30 fixed bytes + name + extra field.
        # The local extra field can differ from the central one, so add slack.
        slack = 1024

        def span(i):
            return i.header_offset, i.header_offset + 30 + len(i.filename.encode()) + \
                len(i.extra) + i.compress_size + slack

        runs, current = [], []
        for info in infos:
            start, stop = span(info)
            if current and (stop - span(current[0])[0] > chunk_bytes
                            or start - span(current[-1])[1] > (1 << 20)):
                runs.append(current)
                current = []
            current.append(info)
        if current:
            runs.append(current)

        done = 0
        for run in runs:
            run_start = span(run[0])[0]
            run_stop = max(span(i)[1] for i in run)
            blob = self.remote.get_range(run_start, run_stop)
            for info in run:
                yield info.filename, _extract(blob, info.header_offset - run_start, info,
                                              self.remote)
                done += 1
            if progress:
                progress(done, len(infos))


def _extract(blob: bytes, offset: int, info: zipfile.ZipInfo, remote: HTTPRangeFile) -> bytes:
    """Decompress one member whose local header starts at blob[offset]."""
    header = blob[offset: offset + 30]
    if header[:4] != b"PK\x03\x04":
        raise RuntimeError(f"{info.filename}: no local file header at {info.header_offset}")
    (_, _, _, method, _, _, _, _, _, name_len, extra_len) = struct.unpack("<4sHHHHHIIIHH", header)
    data_start = offset + 30 + name_len + extra_len
    data = blob[data_start: data_start + info.compress_size]
    if len(data) < info.compress_size:  # the slack was not enough: fetch directly
        absolute = info.header_offset + 30 + name_len + extra_len
        data = remote.get_range(absolute, absolute + info.compress_size)
    if method == zipfile.ZIP_STORED:
        raw = data
    elif method == zipfile.ZIP_DEFLATED:
        raw = zlib.decompressobj(-zlib.MAX_WBITS).decompress(data)
    else:
        raise RuntimeError(f"{info.filename}: unsupported compression method {method}")
    if len(raw) != info.file_size or (zlib.crc32(raw) & 0xFFFFFFFF) != info.CRC:
        raise RuntimeError(f"{info.filename}: size or CRC-32 mismatch")
    return raw
