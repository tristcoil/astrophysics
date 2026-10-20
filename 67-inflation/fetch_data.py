"""Download and cache everything the analysis needs.

What this project reads, all public:

1. **BICEP/Keck BK18** (LAMBDA mirror): the B-mode bandpowers, BK18's own r
   likelihood and component decomposition, and the CosmoMC likelihood data,
   from which we rebuild the BK18 likelihood ourselves.
2. **ACT DR6** MCMC chains (LAMBDA): LCDM runs (ACT, P-ACT, P-ACT-LB) and the
   tensor-mode runs, every one of which includes BK18.
3. **SPT-3G D1** MCMC chains (LAMBDA): SPT alone and CMB-SPA, with and without
   DESI, read out of two large zip archives with HTTP range requests.
4. **Planck 2018** chains from the parameter grid (IRSA's mirror of the Planck
   Legacy Archive), likewise read member by member out of an 11.8 GB zip.
5. Planck's SMICA CMB maps for the horizon-problem sky (hard-linked from
   project D3 when it has them).

Everything lands in ``data/`` with each file's URL, size, SHA-256 and retrieval
time in ``data/manifest.json``; nothing cached is fetched twice, so after one
run the whole project works offline.

    python fetch_data.py            # everything (~0.2 GB new; 1.2 GB more without D3)
    python fetch_data.py --list     # show what would be fetched and what is cached
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import shutil
import struct
import time
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import config as cfg


# ---------------------------------------------------------------------------
# Manifest and plain downloads
# ---------------------------------------------------------------------------


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _load_manifest() -> dict:
    return json.loads(cfg.MANIFEST_PATH.read_text()) if cfg.MANIFEST_PATH.exists() else {}


def _save_manifest(manifest: dict) -> None:
    tmp = cfg.MANIFEST_PATH.with_suffix(".json.part")
    tmp.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    tmp.replace(cfg.MANIFEST_PATH)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _request(url: str, headers: dict | None = None, retries: int = 5,
             timeout: int = 300) -> bytes:
    hdr = {"User-Agent": cfg.USER_AGENT}
    hdr.update(headers or {})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=hdr),
                                        timeout=timeout) as r:
                return r.read()
        except (urllib.error.URLError, TimeoutError, ConnectionError) as err:
            if attempt == retries - 1:
                raise
            print(f"    retrying after {err!r}")
            time.sleep(2 + 4 * attempt)
    raise RuntimeError("unreachable")


def _key(dest: Path) -> str:
    return str(dest.relative_to(cfg.DATA_DIR))


def _adopt(dest: Path, manifest: dict, url: str, **extra) -> None:
    """Record a file that is already on disk but missing from the manifest."""
    if _key(dest) not in manifest:
        manifest[_key(dest)] = {
            "url": url, "bytes": dest.stat().st_size, "sha256": _sha256(dest),
            "retrieved_utc": datetime.fromtimestamp(dest.stat().st_mtime, timezone.utc)
            .isoformat(timespec="seconds"), **extra}


def _link_from_sibling(dest: Path, manifest: dict, url: str) -> bool:
    """Hard-link an identical file from an earlier project's cache, if it has one."""
    src = cfg.SIBLING_COPIES.get(dest)
    if src is None or not src.exists():
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(src, dest)
    except OSError:  # another filesystem: copy instead
        shutil.copy2(src, dest)
    manifest[_key(dest)] = {"url": url, "bytes": dest.stat().st_size, "sha256": _sha256(dest),
                            "linked_from": str(src.relative_to(cfg.ROOT_DIR)),
                            "retrieved_utc": _now()}
    print(f"  linked   {_key(dest)}  (from {src.relative_to(cfg.ROOT_DIR).parts[0]}, "
          "no download)")
    return True


def download(url: str, dest: Path, manifest: dict, *, refresh: bool = False) -> Path:
    """Fetch one file unless cached (or available from a sibling project)."""
    if dest.exists() and not refresh:
        _adopt(dest, manifest, url)
        return dest
    if not refresh and _link_from_sibling(dest, manifest, url):
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    tmp.write_bytes(_request(url))
    tmp.replace(dest)  # atomic: a half-written file is never mistaken for a cached one
    manifest[_key(dest)] = {"url": url, "bytes": dest.stat().st_size, "sha256": _sha256(dest),
                            "retrieved_utc": _now()}
    print(f"  fetched  {_key(dest)}  ({dest.stat().st_size / 1e6:.1f} MB)")
    return dest


# ---------------------------------------------------------------------------
# Reading single members out of large remote zip archives
# ---------------------------------------------------------------------------


class RangeFile(io.RawIOBase):
    """A read-only, seekable view of a remote file through HTTP range requests.

    Python's ``zipfile`` only needs ``seek``/``read``: it reads the central
    directory at the end of the archive, then each wanted member's bytes, so a
    few megabytes come over the wire instead of the whole archive.
    """

    def __init__(self, url: str, block: int = 1 << 20):
        self.url, self.block, self.pos, self._cache = url, block, 0, {}
        self.bytes_fetched = 0
        req = urllib.request.Request(url, headers={"User-Agent": cfg.USER_AGENT,
                                                   "Range": "bytes=0-0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            if r.status != 206:
                raise IOError(f"{url} does not support range requests")
            self.size = int(r.headers["Content-Range"].rsplit("/", 1)[1])

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.pos

    def seek(self, offset: int, whence: int = 0) -> int:
        self.pos = {0: offset, 1: self.pos + offset, 2: self.size + offset}[whence]
        return self.pos

    def _fetch_blocks(self, first: int, last: int) -> None:
        missing = [i for i in range(first, last + 1) if i not in self._cache]
        if not missing:
            return
        a = missing[0] * self.block
        b = min((missing[-1] + 1) * self.block, self.size) - 1
        blob = _request(self.url, {"Range": f"bytes={a}-{b}"})
        if len(blob) != b - a + 1:
            raise IOError(f"short range read from {self.url}")
        self.bytes_fetched += len(blob)
        for i in range(missing[0], missing[-1] + 1):
            off = i * self.block - a
            self._cache[i] = blob[off:off + self.block]
        if len(self._cache) > 256:  # keep memory bounded on long extractions
            for i in sorted(self._cache)[:-64]:
                del self._cache[i]

    def read(self, n: int = -1) -> bytes:
        if n < 0:
            n = self.size - self.pos
        n = max(0, min(n, self.size - self.pos))
        if n == 0:
            return b""
        first, last = self.pos // self.block, (self.pos + n - 1) // self.block
        self._fetch_blocks(first, last)
        out = bytearray()
        while n > 0:
            i, off = divmod(self.pos, self.block)
            chunk = self._cache[i][off:off + n]
            out += chunk
            self.pos += len(chunk)
            n -= len(chunk)
        return bytes(out)

    def readinto(self, b) -> int:
        data = self.read(len(b))
        b[:len(data)] = data
        return len(data)


class _SubFile(io.RawIOBase):
    """A window onto part of another seekable file (a stored zip inside a zip)."""

    def __init__(self, f, start: int, size: int):
        self.f, self.start, self.size, self.pos = f, start, size, 0

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.pos

    def seek(self, offset: int, whence: int = 0) -> int:
        self.pos = {0: offset, 1: self.pos + offset, 2: self.size + offset}[whence]
        return self.pos

    def read(self, n: int = -1) -> bytes:
        if n < 0:
            n = self.size - self.pos
        n = max(0, min(n, self.size - self.pos))
        self.f.seek(self.start + self.pos)
        data = self.f.read(n)
        self.pos += len(data)
        return data

    def readinto(self, b) -> int:
        data = self.read(len(b))
        b[:len(data)] = data
        return len(data)


def open_remote_zip(url: str, inner: str | None = None) -> tuple[zipfile.ZipFile, RangeFile]:
    """Open a remote zip (or a *stored* zip member inside it) without downloading it."""
    rf = RangeFile(url, block=1 << 20)
    outer = zipfile.ZipFile(rf)
    if inner is None:
        return outer, rf
    info = outer.getinfo(inner)
    if info.compress_type != zipfile.ZIP_STORED:
        raise IOError(f"{inner} is compressed inside {url}; cannot be read in place")
    rf.seek(info.header_offset)
    header = rf.read(30)
    name_len, extra_len = struct.unpack("<HH", header[26:30])
    start = info.header_offset + 30 + name_len + extra_len
    return zipfile.ZipFile(_SubFile(rf, start, info.compress_size)), rf


def fetch_zip_folder(url: str, folder: str, dest_dir: Path, manifest: dict, *,
                     inner: str | None = None, keep=lambda name: True,
                     refresh: bool = False, _open: dict = {}) -> list[Path]:
    """Extract the files of one folder of a remote zip that ``keep`` accepts."""
    marker = dest_dir / ".complete"
    if marker.exists() and not refresh:
        return sorted(p for p in dest_dir.iterdir() if p.name != ".complete")
    handle = _open.get((url, inner))
    if handle is None:
        print(f"  reading the central directory of {url.rsplit('/', 1)[-1]}"
              + (f" / {inner}" if inner else ""))
        handle = _open[(url, inner)] = open_remote_zip(url, inner)
    zf, rf = handle
    prefix = folder.rstrip("/") + "/"
    members = [i for i in zf.infolist()
               if i.filename.startswith(prefix) and "/" not in i.filename[len(prefix):]
               and not i.is_dir() and keep(i.filename[len(prefix):])]
    if not members:
        raise FileNotFoundError(f"nothing under {folder} in {url}")
    dest_dir.mkdir(parents=True, exist_ok=True)
    before = rf.bytes_fetched
    out = []
    for info in members:
        dest = dest_dir / info.filename[len(prefix):]
        if not dest.exists() or refresh:
            tmp = dest.with_suffix(dest.suffix + ".part")
            tmp.write_bytes(zf.read(info))
            tmp.replace(dest)
            manifest[_key(dest)] = {
                "url": url, "member": (inner + " :: " if inner else "") + info.filename,
                "bytes": dest.stat().st_size, "sha256": _sha256(dest),
                "crc32_verified": True, "retrieved_utc": _now()}
        out.append(dest)
    marker.write_text(_now())
    mb = (rf.bytes_fetched - before) / 1e6
    print(f"  fetched  {_key(dest_dir)}/  ({len(out)} files, {mb:.1f} MB over the wire)")
    return out


# ---------------------------------------------------------------------------
# Fetch everything
# ---------------------------------------------------------------------------


def _chain_file(name: str) -> bool:
    """Chain samples, parameter names and run settings; not GetDist plots or pickles."""
    low = name.lower()
    if low.endswith((".pkl", ".pdf", ".dill_pickle", ".locked", ".checkpoint", ".covmat")):
        return False
    return True


def fetch_all(*, refresh: bool = False) -> dict:
    manifest = _load_manifest()
    try:
        print("BICEP/Keck BK18 release (LAMBDA mirror)")
        for name in cfg.BK18_FILES:
            download(cfg.BK18_URL + name, cfg.BK18_DIR / name, manifest, refresh=refresh)
        print("ACT DR6 chains (LAMBDA)")
        for name in cfg.ACT_R_CHAINS:
            download(f"{cfg.ACT_URL}r/{name}.tar.gz", cfg.ACT_DIR / f"{name}.tar.gz",
                     manifest, refresh=refresh)
        for name in cfg.ACT_LCDM_CHAINS:
            download(f"{cfg.ACT_URL}lcdm/{name}.tar.gz", cfg.ACT_DIR / f"{name}.tar.gz",
                     manifest, refresh=refresh)
        print("SPT-3G D1 chains (LAMBDA, by byte range)")
        for archive, folder, _label in cfg.SPT_CHAINS.values():
            if archive == "main":
                url, inner = cfg.SPT_URL + cfg.SPT_MAIN_ZIP, cfg.SPT_MAIN_INNER
            else:
                url, inner = cfg.SPT_URL + cfg.SPT_FIXED_ZIP, None
            fetch_zip_folder(url, folder, cfg.SPT_DIR / folder.split("/")[-1], manifest,
                             inner=inner, keep=_chain_file, refresh=refresh)
        print("Planck 2018 parameter-grid chains (IRSA, by byte range)")
        for folder, _label in cfg.PLANCK_CHAINS.values():
            root = folder.replace("/", "_")

            def keep(name: str, root=root) -> bool:  # the main chains, not importance runs
                return name.startswith(root) and "_post_" not in name and (
                    name[len(root):] in (".paramnames", ".properties.ini", ".inputparams",
                                         ".likelihoods", ".ranges")
                    or (name.endswith(".txt") and name[len(root):-4].lstrip("_").isdigit()))

            fetch_zip_folder(cfg.PLANCK_GRID_URL, folder, cfg.PLANCK_DIR / root, manifest,
                             keep=keep, refresh=refresh)
            # Planck's own GetDist summaries of the same run: an answer key.
            fetch_zip_folder(cfg.PLANCK_GRID_URL, folder + "/dist", cfg.PLANCK_DIR / root / "dist",
                             manifest, keep=lambda n, root=root: n in (root + ".margestats",
                                                                       root + ".likestats"),
                             refresh=refresh)
        print("Planck SMICA CMB maps and mask (the horizon-problem sky)")
        for name in cfg.SMICA_FILES:
            download(cfg.SMICA_URL + name, cfg.SKY_DIR / name, manifest, refresh=refresh)
        download(cfg.MASK_URL + cfg.MASK_FILE, cfg.SKY_DIR / cfg.MASK_FILE, manifest,
                 refresh=refresh)
    finally:
        _save_manifest(manifest)
    print("Earlier atlas projects (results read in place)")
    for label, path in (("D3 CAMB fit (our own n_s)", cfg.D3_FIT),
                        ("O1 curvature", cfg.O1_RESULTS)):
        state = "found" if path.exists() else "MISSING (README fallback values will be used)"
        print(f"  {label:26s} {state}")
    return manifest


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--refresh", action="store_true", help="re-download everything")
    args = ap.parse_args(argv)
    fetch_all(refresh=args.refresh)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
