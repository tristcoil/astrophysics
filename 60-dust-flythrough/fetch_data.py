"""Download and cache every input this project uses.

    python fetch_data.py              # everything (~8.1 GB)
    python fetch_data.py main 2kpc    # only some sources

1. **The Edenhofer et al. (2024) 3D dust map** from Zenodo (record 10658339):
   the posterior mean and standard deviation on HEALPix spheres, 3.25 GB for
   the main map (69-1250 pc) and 4.12 GB for the 2 kpc validation map. Zenodo
   serves one connection at ~1.3 MB/s from here, so the files come down as
   32 MB byte ranges over six connections (resumable; MD5-checked against
   Zenodo's checksum at the end).
2. **A second 3D map**: Vergely, Lallement & Cox (2022), the 10 pc-resolution
   cube (3 x 3 x 0.8 kpc, 5 pc sampling, A at 550 nm), from CDS (222 MB).
3. **Cloud distances**: Zucker et al. (2020) Table A.1, 326 sightlines toward
   local molecular clouds (the data Alves et al. 2020 fitted the Radcliffe
   Wave to), and Alves et al.'s 54 "bridging" sightlines (Harvard Dataverse).
4. **The Local Bubble**: O'Neill et al. (2024)'s shell model from the same
   dust map (138 MB FITS table, Harvard Dataverse) and Zucker et al. (2022)'s
   34 cluster tracebacks (Harvard Dataverse).
5. **Young star clusters**: Hunt & Reffert (2023) cluster catalogue (CDS).
6. **SFD (1998) total dust column**, through the dustmaps package's own fetch
   (two 64 MB polar maps), to check the map's column at high latitude.
7. **ZGR23 extinction curve** (Zhang, Green & Rix 2023, Zenodo 7811871), for
   converting the map's units to any band.
8. **Literature**: arXiv sources of the papers we check against.

Files land in ``data/``; nothing is re-downloaded once present, so after one
run the analysis works offline. ``data/manifest.json`` records each file's
URL, size, checksum and download date.
"""

from __future__ import annotations

import concurrent.futures as cf
import datetime as dt
import gzip
import hashlib
import json
import os
import shutil
import sys
import tarfile
import threading
import time
from pathlib import Path

import requests

import config as C

UA = "astro-atlas-fetch/1.0 (python-requests; hobby research project)"
UA_BROWSER = "Mozilla/5.0 (X11; Linux x86_64) astro-atlas hobby project"
UA_CDS = "Wget/1.21 (astro-atlas hobby project)"     # CDS FTP pages refuse browser UAs
CHUNK = 32 * 1024 * 1024
WORKERS = 6

_manifest_lock = threading.Lock()


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def record(dest: Path, url: str, checksum: str | None = None) -> None:
    """Add a file to data/manifest.json (path relative to data/)."""
    if checksum is None:
        checksum = "sha256:" + sha256sum(dest)
    with _manifest_lock:
        manifest = json.loads(C.MANIFEST_PATH.read_text()) if C.MANIFEST_PATH.exists() else {}
        manifest[str(dest.relative_to(C.DATA_DIR))] = {
            "url": url, "bytes": dest.stat().st_size, "checksum": checksum,
            "downloaded": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}
        C.MANIFEST_PATH.write_text(json.dumps(manifest, indent=1, sort_keys=True))


def _hash(path: Path, h, block: int = 16 * 1024 * 1024) -> str:
    with path.open("rb") as fh:
        while chunk := fh.read(block):
            h.update(chunk)
    return h.hexdigest()


def md5sum(path: Path) -> str:
    return _hash(path, hashlib.md5())


def sha256sum(path: Path) -> str:
    return _hash(path, hashlib.sha256())


def download(url: str, dest: Path, *, ua: str = UA, timeout: float = 120.0,
             retries: int = 6) -> Path:
    """Fetch a (small or medium) file once, atomically, with retries."""
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    last = None
    for attempt in range(retries):
        try:
            with requests.get(url, headers={"User-Agent": ua}, stream=True,
                              timeout=(30, timeout), allow_redirects=True) as resp:
                if resp.status_code == 429:
                    time.sleep(float(resp.headers.get("Retry-After", 60)) + 1)
                    continue
                resp.raise_for_status()
                with tmp.open("wb") as fh:
                    for block in resp.iter_content(1024 * 1024):
                        fh.write(block)
            if tmp.stat().st_size == 0:
                raise IOError("empty response")
            tmp.replace(dest)
            record(dest, url)
            print(f"   {dest.relative_to(C.DATA_DIR)} ({dest.stat().st_size / 1e6:.1f} MB)", flush=True)
            return dest
        except (requests.RequestException, IOError) as exc:
            last = exc
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"could not fetch {url}: {last}")


def _fetch_range(session: requests.Session, url: str, fd: int, start: int, stop: int) -> None:
    """Write bytes [start, stop] of ``url`` at the same offset of the open file."""
    for attempt in range(8):
        try:
            headers = {"User-Agent": UA, "Range": f"bytes={start}-{stop}"}
            with session.get(url, headers=headers, stream=True, timeout=(30, 90)) as resp:
                if resp.status_code == 429:
                    time.sleep(float(resp.headers.get("Retry-After", 60)) + 1)
                    continue
                if resp.status_code != 206:
                    raise IOError(f"HTTP {resp.status_code} for range {start}-{stop}")
                pos = start
                for block in resp.iter_content(1024 * 1024):
                    os.pwrite(fd, block, pos)
                    pos += len(block)
                if pos != stop + 1:
                    raise IOError(f"short read {pos - start} of {stop + 1 - start} bytes")
                return
        except (requests.RequestException, IOError) as exc:
            print(f"   range {start // CHUNK}: {exc} (attempt {attempt + 1})", flush=True)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"giving up on bytes {start}-{stop} of {url}")


def download_parallel(url: str, dest: Path, size: int, md5: str) -> Path:
    """Fetch a large file as parallel byte ranges, resumably, and verify its MD5.

    Progress is kept in ``<dest>.chunks.json`` (finished chunk indices), so an
    interrupted download restarts where it stopped.
    """
    if dest.exists() and dest.stat().st_size == size:
        return dest
    part = dest.with_suffix(dest.suffix + ".part")
    state_path = dest.with_suffix(dest.suffix + ".chunks.json")
    done = set(json.loads(state_path.read_text())) if state_path.exists() and part.exists() else set()
    n_chunks = (size + CHUNK - 1) // CHUNK
    todo = [i for i in range(n_chunks) if i not in done]
    print(f"-> {dest.name}: {size / 1e9:.2f} GB, {len(todo)} of {n_chunks} chunks to fetch", flush=True)
    fd = os.open(part, os.O_RDWR | os.O_CREAT, 0o644)
    os.ftruncate(fd, size)
    lock = threading.Lock()
    t0 = time.time()
    session = requests.Session()

    def work(i: int) -> int:
        start = i * CHUNK
        _fetch_range(session, url, fd, start, min(size, start + CHUNK) - 1)
        with lock:
            done.add(i)
            state_path.write_text(json.dumps(sorted(done)))
            rate = (len(done) - (n_chunks - len(todo))) * CHUNK / max(time.time() - t0, 1e-3)
            print(f"   {len(done)}/{n_chunks} chunks, {rate / 1e6:.1f} MB/s", flush=True)
        return i

    try:
        with cf.ThreadPoolExecutor(WORKERS) as pool:
            list(pool.map(work, todo))
    finally:
        os.fsync(fd)
        os.close(fd)
    got = md5sum(part)
    if got != md5:
        state_path.unlink(missing_ok=True)
        raise RuntimeError(f"MD5 mismatch for {dest.name}: {got} != {md5}; delete {part} and retry")
    part.replace(dest)
    state_path.unlink(missing_ok=True)
    record(dest, url, f"md5:{md5}")
    print(f"   {dest.name} OK ({(time.time() - t0) / 60:.1f} min)", flush=True)
    return dest


def dataverse_file(file_id: int, dest: Path) -> Path:
    return download(f"https://dataverse.harvard.edu/api/access/datafile/{file_id}", dest)


def dataverse_listing(doi: str) -> list[dict]:
    url = f"https://dataverse.harvard.edu/api/datasets/:persistentId/?persistentId=doi:{doi}"
    resp = requests.get(url, headers={"User-Agent": UA}, timeout=60)
    resp.raise_for_status()
    return [f["dataFile"] for f in resp.json()["data"]["latestVersion"]["files"]]


# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------

def fetch_edenhofer(which: str) -> Path:
    key, size, md5 = C.EDENHOFER_FILES[which]
    url = C.ZENODO_FILE_URL.format(record=C.ZENODO_RECORD, key=key)
    path = download_parallel(url, C.EDENHOFER_DIR / key, size, md5)
    for small in ("readme.md", "interp2box.py", "interp2lbd.py"):
        download(C.ZENODO_FILE_URL.format(record=C.ZENODO_RECORD, key=small), C.EDENHOFER_DIR / small)
    return path


def fetch_sample_shells() -> None:
    """A few shells of the 12 posterior samples, by HTTP byte range (~300 MB in all).

    Each block is one sample's one HEALPix sphere: 786432 big-endian float32
    (3.1 MB), saved raw as ``sample_shells/s<sample>_r<shell>.f32``.
    """
    url = C.ZENODO_FILE_URL.format(record=C.ZENODO_RECORD, key=C.SAMPLES_FILE)
    out = C.EDENHOFER_DIR / "sample_shells"
    out.mkdir(exist_ok=True)
    nbytes = 12 * 256 ** 2 * 4                        # 786432 pixels x 4 bytes
    jobs = []
    for i in C.SAMPLE_SHELL_PAIRS:
        for shell in (i, i + 1):
            for s in range(12):
                dest = out / f"s{s:02d}_r{shell:03d}.f32"
                if not (dest.exists() and dest.stat().st_size == nbytes):
                    jobs.append((s, shell, dest))
    session = requests.Session()

    def work(job):
        s, shell, dest = job
        start = C.SAMPLES_DATA_OFFSET + (s * 516 + shell) * nbytes
        tmp = dest.with_suffix(".part")
        fd = os.open(tmp, os.O_RDWR | os.O_CREAT | os.O_TRUNC, 0o644)
        try:
            _fetch_block(session, url, fd, start, nbytes)
        finally:
            os.close(fd)
        tmp.replace(dest)
        record(dest, f"{url} bytes={start}-{start + nbytes - 1}")
        return dest

    with cf.ThreadPoolExecutor(WORKERS) as pool:
        for dest in pool.map(work, jobs):
            print(f"   {dest.name}", flush=True)


def _fetch_block(session: requests.Session, url: str, fd: int, start: int, nbytes: int) -> None:
    """Write bytes [start, start + nbytes) of ``url`` to the start of an open file."""
    for attempt in range(8):
        try:
            headers = {"User-Agent": UA, "Range": f"bytes={start}-{start + nbytes - 1}"}
            with session.get(url, headers=headers, stream=True, timeout=(30, 90)) as resp:
                if resp.status_code == 429:
                    time.sleep(float(resp.headers.get("Retry-After", 60)) + 1)
                    continue
                if resp.status_code != 206:
                    raise IOError(f"HTTP {resp.status_code}")
                pos = 0
                for block in resp.iter_content(1024 * 1024):
                    os.pwrite(fd, block, pos)
                    pos += len(block)
                if pos != nbytes:
                    raise IOError(f"short read {pos} of {nbytes}")
                return
        except (requests.RequestException, IOError) as exc:
            print(f"   block at {start}: {exc} (attempt {attempt + 1})", flush=True)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"giving up on bytes {start}+{nbytes} of {url}")


def fetch_vergely() -> Path:
    base = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/664/A174/"
    download(base + "ReadMe", C.VERGELY_DIR / "ReadMe", ua=UA_CDS)
    download(base + "list.dat", C.VERGELY_DIR / "list.dat", ua=UA_CDS)
    return download(base + "fits/" + C.VERGELY_FILE, C.VERGELY_DIR / C.VERGELY_FILE, ua=UA_CDS,
                    timeout=600)


def fetch_clouds() -> None:
    dataverse_file(4361456, C.CATALOG_DIR / "Handbook_Distances_Zucker2020.dat")
    dataverse_file(4361453, C.CATALOG_DIR / "Bridging_Sightlines_Alves_2020.dat")


def fetch_bubble() -> None:
    dataverse_file(8943783, C.CATALOG_DIR / "ONeill2024_LocalBubble_ShellProperties_A0.5.fits")
    out = C.CATALOG_DIR / "zucker2022_tracebacks"
    for f in dataverse_listing("10.7910/DVN/E8PQOD"):
        dataverse_file(f["id"], out / f["filename"])
    # Zucker et al. (2022) Extended Data Tables 1 (clusters) and 2 (birth positions)
    for doi in ("10.7910/DVN/ZU97QD", "10.7910/DVN/1VT8BC"):
        for f in dataverse_listing(doi):
            dataverse_file(f["id"], C.CATALOG_DIR / "zucker2022" / f["filename"])


def fetch_clusters() -> None:
    base = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/673/A114/"
    out = C.CATALOG_DIR / "hunt_reffert_2023"
    download(base + "ReadMe", out / "ReadMe", ua=UA_CDS)
    gz = download(base + "clusters.dat.gz", out / "clusters.dat.gz", ua=UA_CDS)
    dat = out / "clusters.dat"
    if not dat.exists():
        with gzip.open(gz, "rb") as src, dat.open("wb") as dst:
            shutil.copyfileobj(src, dst)


def fetch_sfd() -> None:
    """SFD (1998) E(B-V) polar maps, laid out where the dustmaps package expects them.

    dustmaps' own ``sfd.fetch()`` gets HTTP 403 from Dataverse's search API
    (October 2026); the dataset listing below works, so we fetch the same files.
    """
    target = C.DUSTMAPS_DIR / "sfd"
    for f in dataverse_listing("10.7910/DVN/EWCNL5"):
        if f["filename"] in ("SFD_dust_4096_ngp.fits", "SFD_dust_4096_sgp.fits"):
            dataverse_file(f["id"], target / f["filename"])


def fetch_zgr23() -> None:
    url = "https://zenodo.org/api/records/7811871/files/{}/content"
    for key in ("extinction_curve.txt", "README.md"):
        download(url.format(key), C.CATALOG_DIR / "zgr23" / key)


#: arXiv sources of the papers whose numbers we check (TeX where available).
PAPERS = {
    "2308.01295": "Edenhofer et al. 2024, the 3D dust map",
    "2001.08748": "Alves et al. 2020, the Radcliffe Wave",
    "2402.12596": "Konietzka et al. 2024, the Radcliffe Wave is oscillating",
    "2201.05124": "Zucker et al. 2022, star formation and the Local Bubble",
    "2403.04961": "O'Neill et al. 2024, the Local Bubble is a Local Chimney",
    "2001.00591": "Zucker et al. 2020, distances to local molecular clouds",
    "2303.13424": "Hunt & Reffert 2023, Gaia DR3 star clusters",
    "2208.06469": "Li et al. 2022, kinematics of young stars in the Radcliffe Wave",
}


def fetch_papers() -> None:
    for arxiv_id in PAPERS:
        folder = C.LIT_DIR / arxiv_id
        src = download(f"https://arxiv.org/e-print/{arxiv_id}", folder / "src.tgz", ua=UA_BROWSER)
        with src.open("rb") as fh:
            magic = fh.read(4)
        if magic.startswith(b"%PDF"):                  # PDF-only submissions
            pdf = folder / "paper.pdf"
            if not pdf.exists():
                shutil.copyfile(src, pdf)
        elif not any(folder.glob("*.tex")):
            with tarfile.open(src, "r:gz") as tar:
                tar.extractall(folder, filter="data")


SOURCES = {
    "main": lambda: fetch_edenhofer("main"),
    "2kpc": lambda: fetch_edenhofer("2kpc"),
    "samples": fetch_sample_shells,
    "vergely": fetch_vergely,
    "clouds": fetch_clouds,
    "bubble": fetch_bubble,
    "clusters": fetch_clusters,
    "sfd": fetch_sfd,
    "zgr23": fetch_zgr23,
    "papers": fetch_papers,
}


def main(names: list[str]) -> None:
    for name in names or list(SOURCES):
        print(f"== {name}", flush=True)
        SOURCES[name]()


if __name__ == "__main__":
    main(sys.argv[1:])
