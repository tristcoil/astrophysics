"""Download and cache everything the analysis needs.

What this project reads, all public:

1. **CMS open data** (CERN Open Data portal, CC0): 2012 collisions at 8 TeV in
   reduced NanoAOD format -- four primary-dataset files for the four-lepton
   Higgs search (12.2 GB), the matching Higgs and ZZ simulation (0.56 GB), and
   61.5 million dimuon-triggered events for the particle-zoo spectrum (2.2 GB).
2. **The CMS outreach analysis code** (records 12360 and 12342), which we
   replicate before changing anything.
3. **Literature** (small): the arXiv source of papers whose formulas or tables
   we use (renormalization-group equations, matching conditions, lattice
   results), and HEPData tables of measured Higgs couplings.

The big files are fetched with parallel HTTPS byte ranges (the portal serves
them at ~1.2 MB/s per connection), written in place, resumable after any
interruption, and verified against the adler32 checksums the portal publishes.
Everything lands in ``data/`` with its URL, size, checksums and retrieval time
in ``data/manifest.json``; nothing cached is fetched twice, so after one run
the whole project works offline.

    python fetch_data.py                    # everything (~15 GB)
    python fetch_data.py --only simulation  # just the 0.56 GB of simulation
    python fetch_data.py --only dimuon      # just the dimuon file (2.2 GB)
    python fetch_data.py --only literature  # papers' sources, tables, SMDR (~5 MB)
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import sys
import threading
import time
import urllib.error
import urllib.request
import zipfile
import zlib
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import config as cfg

_manifest_lock = threading.Lock()


# ---------------------------------------------------------------------------
# Manifest and checksums
# ---------------------------------------------------------------------------


def _load_manifest() -> dict:
    return json.loads(cfg.MANIFEST_PATH.read_text()) if cfg.MANIFEST_PATH.exists() else {}


def _save_manifest(manifest: dict) -> None:
    tmp = cfg.MANIFEST_PATH.with_suffix(".json.part")
    tmp.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    tmp.replace(cfg.MANIFEST_PATH)


def _record(key: str, entry: dict) -> None:
    with _manifest_lock:
        manifest = _load_manifest()
        manifest[key] = entry
        _save_manifest(manifest)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def checksums(path: Path) -> tuple[str, str]:
    """(adler32 as 8 hex digits, the way CERN lists it; sha256)."""
    adler, sha = 1, hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            adler = zlib.adler32(block, adler)
            sha.update(block)
    return f"{adler & 0xffffffff:08x}", sha.hexdigest()


def _get(url: str, *, byte_range: tuple[int, int] | None = None, timeout: int = 120) -> bytes:
    headers = {"User-Agent": cfg.USER_AGENT}
    if byte_range is not None:
        headers["Range"] = f"bytes={byte_range[0]}-{byte_range[1]}"
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        if byte_range is not None and r.status != 206:
            raise ConnectionError(f"server ignored the byte range (HTTP {r.status})")
        return r.read()


def _get_retry(url: str, retries: int = 6, **kw) -> bytes:
    for attempt in range(retries):
        try:
            return _get(url, **kw)
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as err:
            if attempt == retries - 1:
                raise
            print(f"    retrying {url.rsplit('/', 1)[-1]} after {err!r}", flush=True)
            time.sleep(3 + 5 * attempt)
    raise RuntimeError("unreachable")


# ---------------------------------------------------------------------------
# Parallel, resumable byte-range downloads of the CMS files
# ---------------------------------------------------------------------------


class _Target:
    """One file being assembled in place from byte ranges."""

    def __init__(self, name: str):
        rec, sub, size, adler, _ = cfg.CMS_FILES[name]
        self.name, self.record, self.size, self.adler = name, rec, size, adler
        self.url = cfg.EOS_HTTPS + sub
        self.dest = cfg.RAW_DIR / f"{name}.root"
        self.part = self.dest.with_suffix(".root.part")
        self.done_path = self.dest.with_suffix(".root.chunks")
        self.n_chunks = -(-size // cfg.CHUNK_BYTES)
        self.lock = threading.Lock()
        self.done = set()
        if self.done_path.exists() and self.part.exists():
            self.done = {int(x) for x in self.done_path.read_text().split()}
        else:
            self.done_path.unlink(missing_ok=True)
        if not self.part.exists():
            with open(self.part, "wb") as f:
                f.truncate(size)          # sparse until filled
        self.fd = os.open(self.part, os.O_RDWR)

    def chunk_range(self, i: int) -> tuple[int, int]:
        start = i * cfg.CHUNK_BYTES
        return start, min(start + cfg.CHUNK_BYTES, self.size) - 1

    def fetch_chunk(self, i: int, progress: dict) -> None:
        """Stream one byte range straight into place (1 MB at a time, so twelve
        connections hold ~12 MB of memory, not twelve whole chunks)."""
        start, end = self.chunk_range(i)
        for attempt in range(8):
            try:
                req = urllib.request.Request(self.url, headers={
                    "User-Agent": cfg.USER_AGENT, "Range": f"bytes={start}-{end}"})
                pos = start
                with urllib.request.urlopen(req, timeout=180) as r:
                    if r.status != 206:
                        raise ConnectionError(f"server ignored the byte range (HTTP {r.status})")
                    while block := r.read(1 << 20):
                        os.pwrite(self.fd, block, pos)
                        pos += len(block)
                if pos != end + 1:
                    raise ConnectionError(f"short read {pos - start} of {end - start + 1}")
                os.fdatasync(self.fd)     # on disk before it is marked done (survives reboots)
                with self.lock:
                    self.done.add(i)
                    with open(self.done_path, "a") as f:
                        f.write(f"{i}\n")
                with _manifest_lock:
                    progress["bytes"] += end - start + 1
                return
            except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as err:
                print(f"    {self.name} chunk {i}: {err!r}; retry {attempt + 1}", flush=True)
                time.sleep(5 + 10 * attempt)
        raise RuntimeError(f"{self.name}: chunk {i} failed 8 times")

    def finish(self) -> bool:
        os.fsync(self.fd)
        os.close(self.fd)
        if len(self.done) != self.n_chunks:
            return False
        adler, sha = checksums(self.part)
        if adler != self.adler:
            print(f"  !! {self.name}: adler32 {adler} != published {self.adler}; "
                  "deleting the copy so the next run starts over", flush=True)
            self.part.unlink()
            self.done_path.unlink(missing_ok=True)
            return False
        self.part.replace(self.dest)
        self.done_path.unlink(missing_ok=True)
        _record(f"raw/{self.dest.name}", {
            "url": self.url, "record": f"https://opendata.cern.ch/record/{self.record}",
            "bytes": self.size, "adler32": adler, "sha256": sha, "retrieved_utc": _now()})
        print(f"  verified {self.dest.name} (adler32 {adler} as published)", flush=True)
        return True


def fetch_cms(names: list[str], workers: int = cfg.DOWNLOAD_WORKERS) -> None:
    manifest = _load_manifest()
    todo = []
    for name in sorted(names, key=lambda n: cfg.CMS_FILES[n][2]):   # small files first
        dest = cfg.RAW_DIR / f"{name}.root"
        if dest.exists() and f"raw/{dest.name}" in manifest:
            continue
        if dest.exists():        # present but unrecorded: verify before trusting it
            adler, sha = checksums(dest)
            if adler == cfg.CMS_FILES[name][3]:
                _record(f"raw/{dest.name}", {
                    "url": cfg.EOS_HTTPS + cfg.CMS_FILES[name][1], "bytes": dest.stat().st_size,
                    "adler32": adler, "sha256": sha, "retrieved_utc": _now()})
                continue
            dest.unlink()
        todo.append(_Target(name))
    if not todo:
        print("  CMS files: all cached and verified")
        return

    jobs = [(t, i) for t in todo for i in range(t.n_chunks) if i not in t.done]
    total = sum(t.chunk_range(i)[1] - t.chunk_range(i)[0] + 1 for t, i in jobs)
    print(f"  CMS files: {len(todo)} to fetch, {total / 1e9:.2f} GB in {len(jobs)} ranges "
          f"over {workers} connections", flush=True)
    progress = {"bytes": 0}
    t0 = time.time()
    remaining = {t.name: {i for i in range(t.n_chunks) if i not in t.done} for t in todo}
    finished = set()

    def work(job):
        t, i = job
        t.fetch_chunk(i, progress)
        return t, i

    stop = threading.Event()

    def report():
        while not stop.wait(60):
            done = progress["bytes"]
            rate = done / max(time.time() - t0, 1e-9)
            eta = (total - done) / rate / 60 if rate > 0 else float("nan")
            print(f"    {done / 1e9:6.2f} / {total / 1e9:.2f} GB  {rate / 1e6:5.1f} MB/s  "
                  f"~{eta:4.0f} min left", flush=True)

    reporter = threading.Thread(target=report, daemon=True)
    reporter.start()
    failed = False
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(work, j) for j in jobs]
        for fut in futures:
            try:
                t, i = fut.result()
            except RuntimeError as err:
                print(f"  !! {err}", flush=True)
                failed = True
                continue
            remaining[t.name].discard(i)
            if not remaining[t.name] and t.name not in finished:
                finished.add(t.name)
                t.finish()
    stop.set()
    for t in todo:
        if t.name not in finished:
            try:
                os.close(t.fd)
            except OSError:
                pass
    if failed:
        raise SystemExit("some byte ranges failed; run fetch_data.py again to resume")


# ---------------------------------------------------------------------------
# The outreach analysis code
# ---------------------------------------------------------------------------


def fetch_code() -> None:
    out = cfg.LITERATURE_DIR / "cms_outreach_code"
    out.mkdir(parents=True, exist_ok=True)
    manifest = _load_manifest()
    for fname, adler in cfg.OUTREACH_CODE_FILES.items():
        dest = out / fname
        key = str(dest.relative_to(cfg.DATA_DIR))
        if dest.exists() and key in manifest:
            continue
        blob = _get_retry(cfg.OUTREACH_CODE_URL + fname)
        got = f"{zlib.adler32(blob) & 0xffffffff:08x}"
        if got != adler:
            raise SystemExit(f"{fname}: adler32 {got} != published {adler}")
        dest.write_bytes(blob)
        _record(key, {"url": cfg.OUTREACH_CODE_URL + fname, "bytes": len(blob),
                      "adler32": got, "sha256": hashlib.sha256(blob).hexdigest(),
                      "retrieved_utc": _now()})
        print(f"  fetched  {key}")
    dest = out / "dimuonSpectrum.py"
    key = str(dest.relative_to(cfg.DATA_DIR))
    if not (dest.exists() and key in manifest):
        blob = _get_retry(cfg.DIMUON_CODE_URL)
        got = f"{zlib.adler32(blob) & 0xffffffff:08x}"
        if got != cfg.DIMUON_CODE_ADLER32:
            raise SystemExit(f"dimuon code zip: adler32 {got} != {cfg.DIMUON_CODE_ADLER32}")
        with zipfile.ZipFile(io.BytesIO(blob)) as z:
            member = next(n for n in z.namelist() if n.endswith("/dimuonSpectrum.py"))
            dest.write_bytes(z.read(member))
        _record(key, {"url": cfg.DIMUON_CODE_URL, "member": member, "bytes": dest.stat().st_size,
                      "zip_adler32": got, "retrieved_utc": _now()})
        print(f"  fetched  {key}")
    # the published plots of record 12360 (vector PDFs, read back bin by bin)
    plots = cfg.OUTREACH_PLOTS_DIR
    key = str(plots.relative_to(cfg.DATA_DIR))
    if not (plots.exists() and key in manifest):
        blob = _get_retry(cfg.OUTREACH_ZIP_URL)
        got = f"{zlib.adler32(blob) & 0xffffffff:08x}"
        if got != cfg.OUTREACH_ZIP_ADLER32:
            raise SystemExit(f"outreach zip: adler32 {got} != {cfg.OUTREACH_ZIP_ADLER32}")
        plots.mkdir(parents=True, exist_ok=True)
        names = []
        with zipfile.ZipFile(io.BytesIO(blob)) as z:
            for member in z.namelist():
                if "/plots/" in member and member.endswith(".pdf"):
                    (plots / Path(member).name).write_bytes(z.read(member))
                    names.append(Path(member).name)
        _record(key, {"url": cfg.OUTREACH_ZIP_URL, "zip_adler32": got,
                      "zip_sha256": hashlib.sha256(blob).hexdigest(),
                      "files": sorted(names), "retrieved_utc": _now()})
        print(f"  fetched  {key}/ ({len(names)} published plots)")


# ---------------------------------------------------------------------------
# Literature and the SMDR cross-check code
# ---------------------------------------------------------------------------


def fetch_literature() -> None:
    manifest = _load_manifest()
    for fname, (url, why) in cfg.LITERATURE.items():
        dest = cfg.LITERATURE_DIR / fname
        key = str(dest.relative_to(cfg.DATA_DIR))
        if dest.exists() and key in manifest:
            continue
        if "hepdata.net" in url:
            # HEPData answers 403 to Python's urllib (and to unknown user
            # agents) but serves the curl command-line tool
            import subprocess
            r = subprocess.run(["curl", "-s", "-f", "-L", "-m", "120", url], capture_output=True)
            if r.returncode != 0:
                raise SystemExit(f"{fname}: curl failed ({r.returncode})")
            blob = r.stdout
        else:
            req = urllib.request.Request(url, headers={"User-Agent": cfg.ARXIV_USER_AGENT})
            for attempt in range(5):
                try:
                    with urllib.request.urlopen(req, timeout=120) as r:
                        blob = r.read()
                    break
                except (urllib.error.URLError, TimeoutError, ConnectionError) as err:
                    if attempt == 4:
                        raise
                    print(f"    retrying {fname} after {err!r}")
                    time.sleep(5 + 5 * attempt)
        dest.write_bytes(blob)
        _record(key, {"url": url, "purpose": why, "bytes": len(blob),
                      "sha256": hashlib.sha256(blob).hexdigest(), "retrieved_utc": _now()})
        print(f"  fetched  {key} ({len(blob) / 1e3:.0f} kB)")


def build_smdr() -> bool:
    """Unpack and compile SMDR (C; gcc and make) into data/cache/smdr. GCC 14+
    defaults to C23, which rejects TSIL's old-style function pointers, so
    the build uses -std=gnu17. Returns False if the build is not possible."""
    import shutil
    import subprocess
    import tarfile
    root = cfg.CACHE_DIR / "smdr"
    src = root / "SMDR-1.3"
    if (src / "calc_fit").exists() and (src / "fig_RGrun_vs_Q").exists():
        return True
    if shutil.which("gcc") is None or shutil.which("make") is None:
        print("  SMDR: no gcc/make on this machine; the cross-check will be skipped")
        return False
    root.mkdir(parents=True, exist_ok=True)
    with tarfile.open(cfg.LITERATURE_DIR / "smdr-1.3.tar.gz") as t:
        t.extractall(root, filter="data")
    for _ in range(2):          # the parallel build occasionally races; a rerun completes it
        subprocess.run(["make", "-j4", "CC=gcc -std=gnu17", "SMDR_OPT=-O2"], cwd=src,
                       capture_output=True, text=True)
    ok = (src / "calc_fit").exists() and (src / "fig_RGrun_vs_Q").exists()
    print(f"  SMDR {'built' if ok else 'build FAILED (cross-check skipped)'} in {src}")
    return ok


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

GROUPS = {
    "simulation": list(cfg.SIMULATION),
    "dimuon": [cfg.DIMUON_FILE],
    "data": list(cfg.DATA_4L),
}


def fetch_all(only: str | None = None, workers: int = cfg.DOWNLOAD_WORKERS) -> None:
    if only in (None, "code"):
        print("Outreach analysis code:")
        fetch_code()
    if only in (None, "literature"):
        print("Literature:")
        fetch_literature()
        build_smdr()
    if only in ("code", "literature"):
        return
    names = GROUPS[only] if only in GROUPS else [n for g in GROUPS.values() for n in g]
    if only in (None, *GROUPS):
        print("CMS open data:")
        fetch_cms(names, workers)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--only", choices=["code", "literature", *GROUPS])
    ap.add_argument("--workers", type=int, default=cfg.DOWNLOAD_WORKERS)
    args = ap.parse_args(argv)
    fetch_all(args.only, args.workers)
    return 0


if __name__ == "__main__":
    sys.exit(main())
