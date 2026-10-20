"""Download and cache every input this project uses (~0.3 GB in ~6,000 files).

    python fetch_data.py                # everything
    python fetch_data.py --only ztf     # one part: ztf, ptf, literature, gaia

Files land in ``data/raw/`` and are never re-downloaded once complete, so after one
run the whole analysis works offline. ``data/manifest.json`` records each file's
source, size, SHA-256 and download date (the thousands of image cutouts are
summarised per directory). See the README for citations.

What comes from where:

* IRSA (ZTF): the list of every ZTF science exposure that covers ZTF J1539+5027,
  then for each one a 61x61-pixel cutout of the *difference image* (exposure minus
  a deep reference) and the PSF model of that difference image: ~5,900 requests.
  IRSA answers most in 2-5 s, but ~18% hang for 3-20 minutes (staged from slow
  storage); the first pass (~1 hour) sets those aside, and ``--patient`` fetches
  them with long timeouts (~40 minutes with 16 workers). ~145 exposures have no
  difference image at all (HTTP 404, recorded in data/raw/ztf/missing.json).
  Interrupt it at any time: a rerun continues where it stopped.
* IRSA (PTF): the Palomar Transient Factory R-band exposures of 2009-2012 that
  cover the target, as 8' cutouts of the calibrated science images.
* arXiv: LaTeX sources of Burdge et al. 2019 (the published ephemeris and
  masses), Kupfer et al. 2024 (the LISA binaries: Tables 1 and 2) and the papers
  whose equations we implement (Robson, Cornish & Liu 2019; Kupfer et al. 2018).
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import io
import json
import sys
import tarfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import config as C

UA_TOOL = "astro-atlas-fetch/1.0 (python-urllib; hobby research project)"
#: arXiv serves e-prints to browser-like agents more reliably.
UA_BROWSER = "Mozilla/5.0 (X11; Linux x86_64; rv:128.0) Gecko/20100101 Firefox/128.0"

_print_lock = threading.Lock()


def log(msg: str) -> None:
    with _print_lock:
        print(msg, flush=True)


def fetch_bytes(url: str, *, ua: str = UA_TOOL, retries: int = 5, timeout: float = 300.0) -> bytes:
    """GET ``url`` with retries; 403/404/410 are final and raise FileNotFoundError."""
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": ua})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as exc:
            if exc.code in (403, 404, 410):
                raise FileNotFoundError(f"{exc.code} for {url}") from exc
            last = exc
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as exc:
            last = exc
        if attempt + 1 < retries:
            time.sleep(min(60.0, 4.0 * 2**attempt))
    raise RuntimeError(f"could not fetch {url}: {last}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


_manifest_lock = threading.Lock()


def record(key: str, entry: dict) -> None:
    """Add or replace one entry in data/manifest.json."""
    with _manifest_lock:
        man = json.loads(C.MANIFEST_PATH.read_text()) if C.MANIFEST_PATH.exists() else {}
        entry = dict(entry, retrieved=dt.date.today().isoformat())
        man[key] = entry
        C.MANIFEST_PATH.write_text(json.dumps(man, indent=1, sort_keys=True))


def record_file(path: Path, url: str, note: str = "") -> None:
    record(str(path.relative_to(C.DATA_DIR)), {"url": url, "bytes": path.stat().st_size,
                                                "sha256": sha256(path), "note": note})


def record_directory(directory: Path, url: str, note: str) -> None:
    """Summarise a directory of many small files: count, bytes, hash of the listing."""
    files = sorted(p for p in directory.iterdir() if p.is_file())
    listing = "\n".join(f"{p.name} {p.stat().st_size}" for p in files).encode()
    record(str(directory.relative_to(C.DATA_DIR)) + "/",
           {"url": url, "files": len(files), "bytes": sum(p.stat().st_size for p in files),
            "sha256_of_listing": hashlib.sha256(listing).hexdigest(), "note": note})


def save(url: str, dest: Path, *, ua: str = UA_TOOL, note: str = "", force: bool = False) -> Path:
    if dest.exists() and dest.stat().st_size > 0 and not force:
        return dest
    data = fetch_bytes(url, ua=ua)
    tmp = dest.with_name(dest.name + ".part")
    tmp.write_bytes(data)
    tmp.replace(dest)
    record_file(dest, url, note)
    log(f"  {dest.relative_to(C.PROJECT_DIR)}: {len(data) / 1e3:.0f} kB")
    return dest


# ---------------------------------------------------------------------------
# ZTF
# ---------------------------------------------------------------------------

ZTF_EXPOSURES = C.ZTF_DIR / "exposures.csv"
ZTF_MISSING = C.ZTF_DIR / "missing.json"


def ztf_product_base(row) -> str:
    """IRSA path stem of one ZTF science exposure (one CCD quadrant)."""
    ff = str(row["filefracday"])
    return (f"{C.ZTF_DATA_URL}/{ff[:4]}/{ff[4:8]}/{ff[8:]}/"
            f"ztf_{ff}_{int(row['field']):06d}_{row['filtercode']}_c{int(row['ccdid']):02d}_"
            f"{row['imgtypecode']}_q{int(row['qid'])}")


def ztf_cutout_paths(row) -> tuple[Path, Path]:
    stem = (f"{row['filefracday']}_{int(row['field']):06d}_{row['filtercode']}_"
            f"c{int(row['ccdid']):02d}_q{int(row['qid'])}")
    return C.ZTF_CUTOUT_DIR / f"{stem}_diff.fits", C.ZTF_CUTOUT_DIR / f"{stem}_psf.fits"


def fetch_ztf_metadata(force: bool = False) -> Path:
    """Every ZTF science exposure whose footprint contains the target."""
    if ZTF_EXPOSURES.exists() and not force:
        return ZTF_EXPOSURES
    url = f"{C.ZTF_SEARCH_URL}?POS={C.TARGET_RA_2016:.6f},{C.TARGET_DEC_2016:.6f}&ct=csv"
    save(url, ZTF_EXPOSURES, note="IRSA IBE search: ZTF science exposures covering the target",
         force=force)
    return ZTF_EXPOSURES


ZTF_SLOW = C.ZTF_DIR / "slow.json"


def _ztf_jobs(rows, skip: set) -> list[tuple[Path, str]]:
    center = f"center={C.TARGET_RA_2016:.6f},{C.TARGET_DEC_2016:.6f}"
    size = f"size={C.ZTF_CUTOUT_PIX}pix"
    jobs = []
    for row in rows:
        diff, psf = ztf_cutout_paths(row)
        base = ztf_product_base(row)
        for dest, url in ((diff, f"{base}_scimrefdiffimg.fits.fz?{center}&{size}&gzip=false"),
                          (psf, f"{base}_diffimgpsf.fits")):
            if (dest.exists() and dest.stat().st_size > 0) or dest.name in skip:
                continue
            jobs.append((dest, url))
    return jobs


def fetch_ztf_cutouts(workers: int = C.ZTF_DOWNLOAD_WORKERS, patient: bool = False) -> None:
    """Difference-image cutout + difference-image PSF for every exposure (resumable).

    IRSA answers most requests in 2-5 s, but some exposures (mostly 2018-2019)
    hang for minutes before the first byte, cutout and PSF alike. A first pass
    therefore uses a short timeout and sets those aside; later passes give them
    up to 20 minutes each, with fewer parallel requests.
    """
    import pandas as pd

    rows = pd.read_csv(ZTF_EXPOSURES).to_dict("records")
    missing = json.loads(ZTF_MISSING.read_text()) if ZTF_MISSING.exists() else {}
    lock = threading.Lock()

    def run_pass(jobs, n_workers, timeout, retries, label):
        slow, t0, done = [], time.time(), 0

        def one(job):
            dest, url = job
            try:
                data = fetch_bytes(url, retries=retries, timeout=timeout)
            except FileNotFoundError:
                with lock:
                    missing[dest.name] = url
                    ZTF_MISSING.write_text(json.dumps(missing, indent=1))
                return "missing"
            except RuntimeError:
                return "slow"
            if not data.startswith(b"SIMPLE"):
                return "slow"  # an HTML error page; try again later
            tmp = dest.with_name(dest.name + ".part")
            tmp.write_bytes(data)
            tmp.replace(dest)
            return "ok"

        with ThreadPoolExecutor(n_workers) as ex:
            futures = {ex.submit(one, j): j for j in jobs}
            for fut in as_completed(futures):
                if fut.result() == "slow":
                    slow.append(futures[fut])
                done += 1
                if done % 200 == 0 or done == len(jobs):
                    rate = done / (time.time() - t0)
                    log(f"  ZTF {label}: {done}/{len(jobs)} ({rate:.2f} files/s, "
                        f"{len(slow)} set aside, ~{(len(jobs) - done) / max(rate, 1e-6) / 60:.0f} min left)")
        return slow

    jobs = _ztf_jobs(rows, set(missing))
    log(f"ZTF: {len(rows)} exposures; {len(jobs)} files still to fetch "
        f"({len(missing)} known to be missing at IRSA)")
    if jobs:
        # --patient: the files left over are the slow ones; skip the short-timeout pass
        slow = jobs if patient else run_pass(jobs, workers, timeout=25.0, retries=1, label="pass 1")
        for n_pass, (n_workers, timeout) in enumerate(((max(workers, 6), 1200.0), (3, 1200.0)), start=2):
            if not slow:
                break
            ZTF_SLOW.write_text(json.dumps([str(d.name) for d, _ in slow], indent=1))
            log(f"ZTF: {len(slow)} slow files; retrying with up to {timeout / 60:.0f} min each")
            slow = run_pass(slow, n_workers, timeout=timeout, retries=2, label=f"pass {n_pass}")
        ZTF_SLOW.write_text(json.dumps([str(d.name) for d, _ in slow], indent=1))
        if slow:
            log(f"ZTF: {len(slow)} files still not delivered; rerun fetch_data.py later "
                "(the analysis uses every exposure whose cutout and PSF are both present)")
    record_directory(C.ZTF_CUTOUT_DIR, C.ZTF_DATA_URL,
                     f"{C.ZTF_CUTOUT_PIX}x{C.ZTF_CUTOUT_PIX}-px difference-image cutouts "
                     "(*_diff.fits) and difference-image PSFs (*_psf.fits), one pair per exposure")


# ---------------------------------------------------------------------------
# PTF
# ---------------------------------------------------------------------------

PTF_EXPOSURES = C.PTF_DIR / "exposures.csv"


def fetch_ptf_metadata(force: bool = False) -> Path:
    url = f"{C.PTF_SEARCH_URL}?POS={C.TARGET_RA_2016:.6f},{C.TARGET_DEC_2016:.6f}&ct=csv"
    return save(url, PTF_EXPOSURES, force=force,
                note="IRSA IBE search: PTF Level 1 science images covering the target")


def ptf_cutout_path(row) -> Path:
    return C.PTF_CUTOUT_DIR / (Path(row["pfilename"]).stem + "_cut.fits")


def fetch_ptf_cutouts(workers: int = 4) -> None:
    import pandas as pd

    rows = pd.read_csv(PTF_EXPOSURES).to_dict("records")
    center = f"center={C.TARGET_RA_2016:.6f},{C.TARGET_DEC_2016:.6f}"
    size = f"size={C.PTF_CUTOUT_ARCMIN:g}arcmin"
    jobs = [(ptf_cutout_path(r), f"{C.PTF_DATA_URL}/{r['pfilename']}?{center}&{size}&gzip=false")
            for r in rows]
    jobs = [j for j in jobs if not (j[0].exists() and j[0].stat().st_size > 0)]
    log(f"PTF: {len(rows)} images; {len(jobs)} cutouts still to fetch")

    def one(job):
        dest, url = job
        data = fetch_bytes(url, retries=4, timeout=600)
        if not data.startswith(b"SIMPLE"):
            raise RuntimeError(f"not a FITS file: {url}")
        dest.write_bytes(data)

    with ThreadPoolExecutor(workers) as ex:
        for i, fut in enumerate(as_completed([ex.submit(one, j) for j in jobs])):
            fut.result()
            if (i + 1) % 20 == 0:
                log(f"  PTF cutouts: {i + 1}/{len(jobs)}")
    record_directory(C.PTF_CUTOUT_DIR, C.PTF_DATA_URL,
                     f"{C.PTF_CUTOUT_ARCMIN:g}' cutouts of the PTF Level 1 science images (R band)")


# ---------------------------------------------------------------------------
# Literature
# ---------------------------------------------------------------------------

#: arXiv e-prints whose LaTeX tables or equations we use (text only; nothing is run).
ARXIV = {
    "1907.11291": "Burdge et al. 2019, Nature 571, 528: discovery, ephemeris, masses",
    "2411.02238": "Shaw et al. 2024, AJ: JWST clock paper with Burdge's updated ephemeris",
    "2302.12719": "Kupfer et al. 2024, ApJ 963, 100: the LISA verification binaries",
    "1805.00482": "Kupfer et al. 2018, MNRAS 480, 302: GW amplitudes and LISA SNRs (check)",
    "1803.01944": "Robson, Cornish & Liu 2019, CQG 36, 105011: LISA sensitivity curves",
    "1801.04275": "Maoz, Hallakoun & Badenes 2018, MNRAS 476, 2584: DWD population",
    "astro-ph/0312577": "Marsh, Nelemans & Steeghs 2004, MNRAS 350, 113: mass-transfer stability",
}
OBRIEN_40PC = C.LIT_DIR / "obrien2024_40pc_white_dwarfs.tsv"


def arxiv_dir(arxiv_id: str) -> Path:
    return C.LIT_DIR / f"arxiv_{arxiv_id.replace('/', '_')}"


def fetch_literature() -> None:
    for arxiv_id, note in ARXIV.items():
        dest = arxiv_dir(arxiv_id)
        archive = C.LIT_DIR / f"arxiv_{arxiv_id.replace('/', '_')}.src"
        save(f"https://arxiv.org/e-print/{arxiv_id}", archive, ua=UA_BROWSER, note=note)
        if dest.exists() and any(dest.iterdir()):
            continue
        dest.mkdir(parents=True, exist_ok=True)
        raw = archive.read_bytes()
        try:
            with tarfile.open(fileobj=io.BytesIO(raw), mode="r:*") as tar:
                tar.extractall(dest, filter="data")
        except tarfile.ReadError:  # a single gzipped .tex file
            import gzip

            (dest / "main.tex").write_bytes(gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw)
        log(f"  extracted {archive.name} -> {dest.relative_to(C.PROJECT_DIR)}")
    # O'Brien et al. 2024: every confirmed white dwarf within 40 pc (updated table)
    url = ("https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/MNRAS/527/8687/tablea1"
           "&-out.max=unlimited&-out=Name,GaiaDR3,Plx,e_Plx,GMAG,BP-RP,SpType,Comp,Teff,Mass,e_Mass,Massc")
    save(url, OBRIEN_40PC, note="O'Brien et al. 2024 (MNRAS 527, 8687) 40 pc white dwarfs, VizieR tablea1")


PTF_GAIA_FIELD = C.PTF_DIR / "gaia_dr3_field.tsv"


def fetch_ptf_field_stars(force: bool = False) -> Path:
    """Gaia DR3 stars around the target (comparison stars for the PTF photometry)."""
    url = ("https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=I/355/gaiadr3"
           f"&-c={C.TARGET_RA_2016:.6f}%20{C.TARGET_DEC_2016:+.6f}&-c.rm=5.5&-out.max=5000"
           "&-out=Source,RA_ICRS,DE_ICRS,Gmag,BPmag,RPmag,BP-RP,Plx,pmRA,pmDE,RUWE,VarFlag"
           "&Gmag=%3C20.5")
    return save(url, PTF_GAIA_FIELD, force=force,
                note="VizieR I/355/gaiadr3: Gaia DR3 sources within 5.5' of the target, G < 20.5")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--only", choices=["ztf", "ptf", "literature", "gaia"], default=None)
    parser.add_argument("--workers", type=int, default=C.ZTF_DOWNLOAD_WORKERS)
    parser.add_argument("--patient", action="store_true",
                        help="ZTF: give every remaining file up to 20 minutes from the start (IRSA stages "
                             "some exposures from slow storage)")
    args = parser.parse_args()
    parts = [args.only] if args.only else ["literature", "gaia", "ptf", "ztf"]
    for part in parts:
        print(f"== {part}", flush=True)
        if part == "ztf":
            fetch_ztf_metadata()
            fetch_ztf_cutouts(args.workers, patient=args.patient)
        elif part == "ptf":
            fetch_ptf_metadata()
            fetch_ptf_field_stars()
            fetch_ptf_cutouts()
        elif part == "literature":
            fetch_literature()
    print("done.")


if __name__ == "__main__":
    sys.exit(main())
