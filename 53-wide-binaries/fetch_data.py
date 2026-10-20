"""Download and cache every input of the project.

    python fetch_data.py            # everything not yet cached
    python fetch_data.py --catalog  # only the El-Badry et al. (2021) catalogues

Each file lands in data/ with its URL, retrieval date, size and SHA-256
recorded in data/provenance.json, and is never fetched twice.

The two big Zenodo files (1.42 GB + 0.35 GB) come down as parallel byte
ranges: on 9 Oct 2026 a single connection to Zenodo ran at ~110 kB/s (4.5 h),
while each ranged request ran at ~1.4 MB/s. Zenodo allows 133 requests a
minute per IP, so chunks are large (32 MB) and the worker count small.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import requests

from config import PROVENANCE, RAW_DIR, ZENODO_FILES, ZENODO_RECORD

USER_AGENT = "astrophysics-atlas-M3/1.0 (hobby project; python-requests)"
CHUNK = 32 * 2**20
WORKERS = 6
_prov_lock = threading.Lock()


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(2**22), b""):
            h.update(block)
    return h.hexdigest()


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(2**22), b""):
            h.update(block)
    return h.hexdigest()


def record_provenance(path: Path, url: str, **extra) -> None:
    """Add/replace this file's entry in data/provenance.json."""
    with _prov_lock:
        prov = json.loads(PROVENANCE.read_text()) if PROVENANCE.exists() else {}
        key = str(path.relative_to(PROVENANCE.parent))
        prov[key] = dict(url=url, retrieved=now_utc(), bytes=path.stat().st_size,
                         sha256=sha256(path), **extra)
        tmp = PROVENANCE.with_suffix(".part")
        tmp.write_text(json.dumps(prov, indent=1, sort_keys=True))
        tmp.replace(PROVENANCE)


def _get_with_retry(session: requests.Session, url: str, headers: dict, tries: int = 12,
                    timeout: float = 120.0) -> requests.Response:
    delay = 5.0
    for attempt in range(tries):
        try:
            r = session.get(url, headers=headers, stream=True, timeout=timeout)
        except requests.RequestException as exc:
            if attempt == tries - 1:
                raise
            print(f"    {type(exc).__name__}; retrying in {delay:.0f} s", flush=True)
            time.sleep(delay)
            delay = min(delay * 1.7, 120)
            continue
        if r.status_code == 429:  # rate limited: honour Retry-After
            wait = float(r.headers.get("Retry-After", 60))
            r.close()
            print(f"    429 from Zenodo; waiting {wait:.0f} s", flush=True)
            time.sleep(wait + 1)
            continue
        if r.status_code >= 500:
            r.close()
            time.sleep(delay)
            delay = min(delay * 1.7, 120)
            continue
        return r
    raise RuntimeError(f"giving up on {url}")


def ranged_download(url: str, dest: Path, size: int, *, chunk: int = CHUNK,
                    workers: int = WORKERS) -> None:
    """Fetch ``url`` (``size`` bytes) as parallel byte ranges into ``dest``.

    Each range is kept as its own file in ``dest.chunks/`` until the whole set
    is complete, so an interrupted download resumes where it stopped."""
    part_dir = dest.with_name(dest.name + ".chunks")
    part_dir.mkdir(exist_ok=True)
    ranges = [(i, lo, min(lo + chunk, size) - 1) for i, lo in enumerate(range(0, size, chunk))]
    session = requests.Session()
    session.headers["User-Agent"] = USER_AGENT
    t0 = time.time()
    done_bytes = [sum(p.stat().st_size for p in part_dir.glob("*.bin"))]
    lock = threading.Lock()

    def one(job):
        i, lo, hi = job
        out = part_dir / f"{i:05d}.bin"
        want = hi - lo + 1
        if out.exists() and out.stat().st_size == want:
            return
        for attempt in range(6):
            r = _get_with_retry(session, url, {"Range": f"bytes={lo}-{hi}"})
            if r.status_code != 206:
                r.close()
                raise RuntimeError(f"expected 206 for a range, got {r.status_code}")
            tmp = out.with_suffix(".tmp")
            got = 0
            try:
                with open(tmp, "wb") as f:
                    for block in r.iter_content(2**20):
                        f.write(block)
                        got += len(block)
            except requests.RequestException as exc:
                print(f"    chunk {i}: {type(exc).__name__} after {got:,} B; retrying", flush=True)
                time.sleep(5)
                continue
            if got != want:
                print(f"    chunk {i}: {got:,} of {want:,} B; retrying", flush=True)
                continue
            tmp.replace(out)
            with lock:
                done_bytes[0] += want
                rate = done_bytes[0] / max(time.time() - t0, 1e-3) / 2**20
                print(f"    {dest.name}: {done_bytes[0] / 2**20:7.0f} / {size / 2**20:.0f} MB "
                      f"({rate:.1f} MB/s incl. resumed)", flush=True)
            return
        raise RuntimeError(f"chunk {i} of {dest.name} failed repeatedly")

    with ThreadPoolExecutor(workers) as pool:
        list(pool.map(one, ranges))
    tmp = dest.with_suffix(dest.suffix + ".part")
    with open(tmp, "wb") as fout:
        for i, _, _ in ranges:
            with open(part_dir / f"{i:05d}.bin", "rb") as fin:
                shutil.copyfileobj(fin, fout, 2**22)
    if tmp.stat().st_size != size:
        raise RuntimeError(f"{dest.name}: assembled {tmp.stat().st_size} B, expected {size}")
    tmp.replace(dest)
    shutil.rmtree(part_dir)


def fetch_zenodo(name: str) -> Path:
    """One file of El-Badry et al.'s Zenodo record, MD5-checked."""
    info = ZENODO_FILES[name]
    dest = RAW_DIR / name
    url = f"https://zenodo.org/api/records/{ZENODO_RECORD}/files/{name}/content"
    if dest.exists() and dest.stat().st_size == info["size"]:
        return dest
    print(f"  downloading {name} ({info['size'] / 1e6:,.0f} MB) from Zenodo ...", flush=True)
    if info["size"] > 4 * CHUNK:
        ranged_download(url, dest, info["size"])
    else:
        r = _get_with_retry(requests.Session(), url, {"User-Agent": USER_AGENT})
        r.raise_for_status()
        dest.write_bytes(r.content)
    got = md5(dest)
    if got != info["md5"]:
        bad = dest.with_suffix(dest.suffix + ".badmd5")
        dest.rename(bad)
        raise RuntimeError(f"{name}: MD5 {got} != Zenodo's {info['md5']} (kept as {bad.name})")
    record_provenance(dest, f"https://zenodo.org/records/{ZENODO_RECORD}", md5=got,
                      citation="El-Badry, Rix & Heintz 2021, MNRAS 506, 2269")
    print(f"    {name}: MD5 verified", flush=True)
    return dest


def fetch_catalogs() -> None:
    # The two small code files document how the catalogue was built.
    for name in ("find_binaries_edr3.py", "num_neighbors_edr3.py",
                 "all_columns_catalog_shift.fits.gz", "all_columns_catalog.fits.gz"):
        fetch_zenodo(name)


# --------------------------------------------------------------------------
# Smaller inputs
# --------------------------------------------------------------------------
SFD_URL = "https://raw.githubusercontent.com/kbarbary/sfddata/master/{name}"
SFD_FILES = {"SFD_dust_4096_ngp.fits": 67_115_520, "SFD_dust_4096_sgp.fits": 67_115_520}
HWANG_GDRIVE_ID = "1h1pgexHUTpvE80PWCh6m1zY3QMMuYPnd"
PECAUT_URL = "https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt"


def _simple_get(url: str, dest: Path, *, expect_size: int | None = None, **prov) -> Path:
    if dest.exists() and (expect_size is None or dest.stat().st_size == expect_size):
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  downloading {dest.name} ...", flush=True)
    r = _get_with_retry(requests.Session(), url, {"User-Agent": USER_AGENT}, timeout=600)
    r.raise_for_status()
    tmp = dest.with_suffix(dest.suffix + ".part")
    with open(tmp, "wb") as f:
        for block in r.iter_content(2**20):
            f.write(block)
    if expect_size is not None and tmp.stat().st_size != expect_size:
        raise RuntimeError(f"{dest.name}: {tmp.stat().st_size} B, expected {expect_size}")
    tmp.replace(dest)
    record_provenance(dest, url, **prov)
    return dest


def fetch_sfd() -> None:
    """Schlegel, Finkbeiner & Davis (1998) reddening maps (Lambert projections
    of the two Galactic hemispheres), from the copy used by the sfdmap and
    sncosmo packages."""
    for name, size in SFD_FILES.items():
        _simple_get(SFD_URL.format(name=name), RAW_DIR / "sfd" / name, expect_size=size,
                    citation="Schlegel, Finkbeiner & Davis 1998, ApJ 500, 525")


def fetch_hwang() -> Path:
    """Hwang, Ting & Zakamska (2022) eccentricity estimates for every pair of
    El-Badry et al. (2021), linked from github.com/HC-Hwang/Eccentricity-of-
    wide-binaries (a Google Drive file; large files need the confirm form)."""
    dest = RAW_DIR / "hwang2022" / "hwang2022_eccentricities.fits"
    if dest.exists() and dest.stat().st_size > 2e8:
        return dest
    import re
    s = requests.Session()
    s.headers["User-Agent"] = USER_AGENT
    probe = s.get("https://drive.google.com/uc", params=dict(export="download", id=HWANG_GDRIVE_ID),
                  timeout=120)
    uuid = re.search(r'name="uuid" value="([^"]+)"', probe.text)
    params = dict(id=HWANG_GDRIVE_ID, export="download", confirm="t")
    if uuid:
        params["uuid"] = uuid.group(1)
    r = s.get("https://drive.usercontent.google.com/download", params=params, stream=True, timeout=600)
    r.raise_for_status()
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".part")
    with open(tmp, "wb") as f:
        for block in r.iter_content(2**20):
            f.write(block)
    tmp.replace(dest)
    record_provenance(dest, f"https://drive.google.com/file/d/{HWANG_GDRIVE_ID}",
                      citation="Hwang, Ting & Zakamska 2022, MNRAS 512, 3383")
    return dest


def fetch_pecaut() -> Path:
    return _simple_get(PECAUT_URL, RAW_DIR / "pecaut_mamajek" / "EEM_dwarf_UBVIJHK_colors_Teff.txt",
                       citation="Pecaut & Mamajek 2013, ApJS 208, 9 (table version 2022.04.16)")


#: Chae's public samples and code (Zenodo record, file) with MD5 where Zenodo lists one
CHAE_FILES = [
    (8416435, "README", None), (8416435, "gaia_dr3_MSMS_d200pc.csv", "8db3561bc1624d423c420efb96f93de6"),
    (8416435, "gaia_dr3_MSMS_purewb.csv", "2373d5132f1f5db2ff885818a86d636b"),
    (8416435, "M_MagG_polyfit_V_Ic.dat", None), (8416435, "M_MagG_polyfit_J_Bp_Rp.dat", None),
    (8416435, "M_MagG_polyfit_F.dat", None), (8416435, "MagG_M_polyfit_V_Ic.dat", None),
    (8416435, "MagG_M_polyfit_J_Bp_Rp.dat", None), (8416435, "MagG_M_polyfit_F.dat", None),
    (8416435, "velprofile_purewb.py", None),
    (10062232, "gaia_dr3_MSMS_purewb_revised.csv", None), (10062232, "velprofile_wb_v3.py", None),
    (10987140, "gaia_dr3_MSMS_purewb_ruwe.csv", None), (10987140, "erratum.pdf", None),
]
SAAD_TING_URL = "https://raw.githubusercontent.com/seratsaad/wb3d-gamma/main/{path}"
SAAD_TING_FILES = {"data/README.md": "README.md", "data/chae_2026_data.csv": "chae_2026_data.csv",
                   "data/chae_2026_gaia.csv": "chae_2026_gaia.csv",
                   "posterior_samples.npz": "posterior_samples_saad_ting.npz",
                   "run_gamma_analysis.py": "run_gamma_analysis.py"}


def fetch_chae() -> None:
    """Chae (2023, 2024a) samples, mass-magnitude polynomials and code."""
    for rec, name, want in CHAE_FILES:
        dest = RAW_DIR / "chae_zenodo" / str(rec) / name
        if dest.exists() and dest.stat().st_size > 0 and (want is None or md5(dest) == want):
            continue
        url = f"https://zenodo.org/api/records/{rec}/files/{name}/content"
        dest.unlink(missing_ok=True)
        _simple_get(url, dest, citation="Chae 2023, ApJ 952, 128; Chae 2024, ApJ 960, 114",
                    zenodo_record=rec)
        if want is not None and md5(dest) != want:
            raise RuntimeError(f"{name}: MD5 mismatch")


def fetch_chae2026_3d() -> None:
    """The 36 pairs with precise RVs of Chae et al. (2026, arXiv:2601.21728), as
    tabulated in Saad & Ting's (2026, arXiv:2603.11015) public repository."""
    for path, name in SAAD_TING_FILES.items():
        _simple_get(SAAD_TING_URL.format(path=path), RAW_DIR / "chae2026_3d" / name,
                    citation="Chae et al. 2026 (arXiv:2601.21728); Saad & Ting 2026 (arXiv:2603.11015)")


def ensure_data(quiet: bool = False) -> None:
    fetch_catalogs()
    fetch_hwang()
    fetch_pecaut()
    fetch_sfd()
    fetch_chae()
    fetch_chae2026_3d()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--catalog", action="store_true", help="only the Zenodo catalogues")
    args = ap.parse_args(argv)
    if args.catalog:
        fetch_catalogs()
    else:
        ensure_data()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
