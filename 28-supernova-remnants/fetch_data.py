"""Download and cache everything the project needs.

About 0.2 GB in total, all public, no accounts:

1. **Plates** (STScI): Digitized Sky Survey scans of the Palomar Schmidt plates
   of the Crab (1951, 1983, 1991, 1992, 1998) and Cas A (1954, 1991), and
   Pan-STARRS1 stacked images of both (2009-2014).
2. **Gaia DR3** (CDS VizieR): every star down to G = 20 around both remnants,
   with proper motions, to wind back to each plate's date and register it.
3. **Chandra** (Chandra Data Archive): the pipeline's 1024 x 1024 "centre
   image" of Cas A from every full-frame ACIS-S observation, 2000-2019, and
   the field-of-view outline of each.
4. **Catalogues**: Green's catalogue of 310 Galactic supernova remnants
   (2024 October), and the ATNF Pulsar Catalogue v2.8.1.

Everything lands in ``data/`` with ``manifest.json`` recording each file's
source URL, retrieval time and SHA-256 hash. Files already present are never
re-downloaded; ``--refresh`` forces it. Requests go one at a time and back off
when a server is busy, so this is gentle on the archives.

    python fetch_data.py                  # everything
    python fetch_data.py --only chandra   # plates | gaia | chandra | catalogs
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import http.client
import json
import math
import os
import re
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone

import config
from config import CATALOG_DIR, CHANDRA_DIR, DATA_DIR, MANIFEST_PATH, TARGETS, USER_AGENT, Source

GROUPS = ("plates", "gaia", "chandra", "catalogs")


class BadPayload(ValueError):
    """The server answered, but not with the file we asked for."""


def _check_payload(source: Source, head: bytes) -> None:
    """Cheap check that the bytes are the expected kind of file, not an HTML
    error page served with status 200."""
    text = head[:4000].decode("latin-1", errors="replace")
    if source.kind == "fits" and not head.startswith(b"SIMPLE"):
        raise BadPayload(f"{source.filename}: not a FITS file ({' '.join(text[:200].split())!r})")
    if source.kind == "gz" and not head.startswith(b"\x1f\x8b"):
        raise BadPayload(f"{source.filename}: not gzip ({text[:120]!r})")
    if source.kind == "tsv" and ("#INFO" not in text and "#RESOURCE" not in text):
        raise BadPayload(f"{source.filename}: not a VizieR TSV answer ({text[:200]!r})")
    if source.kind == "html" and "<" not in text:
        raise BadPayload(f"{source.filename}: not HTML ({text[:120]!r})")


def _get(url: str, *, timeout: float = 120.0) -> bytes:
    """A small GET (directory listings, file-name lookups), with retries."""
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(1, 6):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read()
        except (urllib.error.URLError, http.client.HTTPException, OSError) as err:
            if isinstance(err, urllib.error.HTTPError) and err.code < 500 and err.code != 429:
                raise
            if attempt == 5:
                raise
            time.sleep(3.0 * attempt)
    raise AssertionError("unreachable")


def _download(source: Source, dest, *, attempts: int = 6) -> tuple[int, str]:
    """Stream one URL to ``dest``; return (bytes, sha256). Retries with backoff.

    HTTP 429 and 5xx answers are retried after a growing pause; other 4xx
    answers mean the request itself is wrong and fail at once.
    """
    url = source.url
    if source.params:
        url += "?" + urllib.parse.urlencode(source.params)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            sha = hashlib.sha256()
            n = 0
            with urllib.request.urlopen(request, timeout=source.timeout) as response, open(dest, "wb") as out:
                first = True
                while chunk := response.read(1 << 20):
                    if first:
                        _check_payload(source, chunk)
                        first = False
                    out.write(chunk)
                    sha.update(chunk)
                    n += len(chunk)
            if n == 0:
                raise BadPayload(f"{source.filename}: empty response")
            return n, sha.hexdigest()
        except urllib.error.HTTPError as err:
            last_error = err
            if err.code != 429 and err.code < 500:
                break
            wait = float(err.headers.get("Retry-After") or 10 * attempt)
        except (urllib.error.URLError, http.client.HTTPException, OSError, BadPayload) as err:
            last_error = err
            if isinstance(err, BadPayload) and attempt >= 2:
                break  # the same wrong answer twice: the request is the problem
            wait = 3.0 * attempt
        if attempt < attempts:
            time.sleep(wait)
    raise RuntimeError(f"could not download {source.filename} from {url[:160]}: {last_error}")


def _load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text())
    return {"files": {}}


def _record(filename: str, entry: dict) -> None:
    """Add one manifest entry under a lock (re-read, update, write atomically)."""
    lock = MANIFEST_PATH.with_suffix(".lock")
    with open(lock, "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        manifest = _load_manifest()
        manifest["files"][filename] = entry
        tmp = MANIFEST_PATH.with_suffix(f".tmp{os.getpid()}")
        tmp.write_text(json.dumps(manifest, indent=1))
        tmp.replace(MANIFEST_PATH)
        fcntl.flock(fh, fcntl.LOCK_UN)


def fetch_sources(sources: list[Source], *, refresh: bool = False) -> list[str]:
    """Make sure every source file is cached. Returns the list of failures.

    Downloads go to a temporary name unique to this process and are renamed
    when complete, so an interrupted run never leaves a truncated file that
    later looks cached.
    """
    failures = []
    for source in sources:
        path = DATA_DIR / source.filename
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists() and not refresh:
            continue
        if source.pause:
            time.sleep(source.pause)
        print(f"  fetching   {source.filename} ...", end="", flush=True)
        tmp = path.with_name(f"{path.name}.part{os.getpid()}")
        try:
            nbytes, sha = _download(source, tmp)
        except RuntimeError as err:
            tmp.unlink(missing_ok=True)
            failures.append(str(err))
            print(" FAILED")
            continue
        tmp.replace(path)
        url = source.url + ("?" + urllib.parse.urlencode(source.params) if source.params else "")
        _record(source.filename, {
            "url": url,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "bytes": nbytes,
            "sha256": sha,
            "what": source.what,
            "reference": source.reference,
        })
        print(f" {nbytes / 1e6:.2f} MB" if nbytes > 1e6 else f" {nbytes / 1e3:.0f} kB")
    return failures


# --------------------------------------------------------------------------
# 1. Plates
# --------------------------------------------------------------------------

DSS_REFERENCE = ("Digitized Sky Survey, STScI (Lasker et al. 1990; POSS-I: National Geographic "
                 "Society-Palomar; POSS-II: Reid et al. 1991)")
PS1_REFERENCE = "Pan-STARRS1 3pi survey stacks, DR2 (Chambers et al. 2016; Flewelling et al. 2020)"


def dss_sources() -> list[Source]:
    out = []
    for plate in config.PLATES:
        if plate.survey.startswith("ps1-"):
            continue
        t = TARGETS[plate.target]
        out.append(Source(
            f"plates/{plate.key}.fits", config.DSS_URL,
            f"{plate.label}: {t.box_arcmin:.0f}' x {t.box_arcmin:.0f}' DSS scan around {t.name}",
            DSS_REFERENCE, kind="fits", pause=1.0, timeout=180,
            params={"v": plate.survey, "r": f"{t.ra:.6f}", "d": f"{t.dec:.6f}", "e": "J2000",
                    "h": f"{t.box_arcmin:.1f}", "w": f"{t.box_arcmin:.1f}", "f": "fits", "c": "none"}))
    return out


def _ps1_skycells(t: config.Target, band: str) -> list[str]:
    """Every Pan-STARRS stack file whose skycell overlaps the cutout box.

    ``ps1filenames.py`` returns the one skycell whose centre is nearest a
    position, so we ask at the box centre, its corners and edge midpoints."""
    files = []
    half = t.box_arcmin / 2 / 60
    for fx in (-1, 0, 1):
        for fy in (-1, 0, 1):
            dec = t.dec + fy * half
            ra = t.ra + fx * half / max(0.2, math.cos(math.radians(dec)))
            q = urllib.parse.urlencode({"ra": f"{ra:.6f}", "dec": f"{dec:.6f}", "filters": band})
            text = _get(f"{config.PS1_FILENAMES_URL}?{q}").decode()
            for line in text.splitlines()[1:]:
                parts = line.split()
                if len(parts) >= 8 and parts[7] not in files:
                    files.append(parts[7])
            time.sleep(0.3)
    return files


def ps1_sources() -> list[Source]:
    out = []
    for plate in config.PLATES:
        if not plate.survey.startswith("ps1-"):
            continue
        t = TARGETS[plate.target]
        band = plate.survey.split("-")[1]
        size_px = int(round(t.box_arcmin * 60 / 0.25))
        for fname in _ps1_skycells(t, band):
            cell = re.search(r"skycell\.(\d+)\.(\d+)", fname)
            tag = f"{cell.group(1)}.{cell.group(2)}" if cell else "x"
            out.append(Source(
                f"plates/{plate.key}_{tag}.fits", config.PS1_FITSCUT_URL,
                f"{plate.label}: Pan-STARRS stack {fname.rsplit('/', 1)[-1]}, cut to "
                f"{t.box_arcmin:.0f}' around {t.name}",
                PS1_REFERENCE, kind="fits", pause=0.5, timeout=600,
                params={"red": fname, "format": "fits", "ra": f"{t.ra:.6f}", "dec": f"{t.dec:.6f}",
                        "size": str(size_px)}))
    return out


# --------------------------------------------------------------------------
# 2. Gaia
# --------------------------------------------------------------------------

def gaia_sources() -> list[Source]:
    out = []
    for key, t in TARGETS.items():
        params = {
            "-source": "I/355/gaiadr3",
            "-c": f"{t.ra:.6f} {t.dec:+.6f}",
            "-c.rm": f"{config.GAIA_RADIUS_ARCMIN[key]:.1f}",
            "-out": config.GAIA_COLUMNS,
            "-out.max": "unlimited",
            "Gmag": f"<{config.GAIA_GMAX}",
        }
        out.append(Source(
            f"gaia/{key}_gaia_dr3.tsv", config.VIZIER_ASU_URL,
            f"Gaia DR3 stars with G < {config.GAIA_GMAX} within {config.GAIA_RADIUS_ARCMIN[key]}' of {t.name}",
            "Gaia Collaboration, Vallenari et al. 2023, A&A 674, A1 (via CDS VizieR I/355)",
            kind="tsv", timeout=300, params=params))
    return out


# --------------------------------------------------------------------------
# 3. Chandra
# --------------------------------------------------------------------------

def chandra_sources() -> list[Source]:
    """The centre image and field-of-view file of every observation in
    ``config.CHANDRA_EPOCHS``. File names carry a processing version
    (``N00x``), so each observation's directory listing is read first."""
    out = []
    for epoch, obsids in config.CHANDRA_EPOCHS.items():
        for obsid in obsids:
            have = sorted(CHANDRA_DIR.glob(f"{obsid:05d}_*"))
            if len(have) >= 2:  # both files cached: skip the listing
                out.extend(Source(f"chandra/{p.name}", "", "", "") for p in have)
                continue
            base = f"{config.CHANDRA_FTP}/{obsid % 10}/{obsid}/primary/"
            listing = _get(base).decode()
            for pattern, what in ((r'href="(acisf\d+N\d+_cntr_img2\.fits\.gz)"',
                                   "centre image: every event at 0.492\" pixels, all energies"),
                                  (r'href="(acisf\d+_\d+N\d+_fov1\.fits\.gz)"',
                                   "field of view: outline of each CCD on the sky")):
                m = re.search(pattern, listing)
                if not m:
                    raise RuntimeError(f"ObsID {obsid}: no file matching {pattern} in {base}")
                name = m.group(1)
                out.append(Source(f"chandra/{obsid:05d}_{name}", base + name,
                                  f"Chandra ACIS-S ObsID {obsid} (epoch {epoch}): {what}",
                                  "Chandra Data Archive (CXC), standard pipeline products",
                                  kind="gz", pause=0.5, timeout=300))
            time.sleep(0.3)
    return out


# --------------------------------------------------------------------------
# 4. Catalogues
# --------------------------------------------------------------------------

def catalog_sources() -> list[Source]:
    return [
        Source("catalogs/green_snrs_2024oct.html", config.GREEN_URL,
               "Green's catalogue of Galactic SNRs, summary data table (2024 October version)",
               config.GREEN_REFERENCE, kind="html"),
        Source(f"catalogs/psrcat_pkg.v{config.ATNF_VERSION}.tar.gz", config.ATNF_URL,
               f"ATNF Pulsar Catalogue v{config.ATNF_VERSION} (database + software)",
               "Manchester, Hobbs, Teoh & Hobbs 2005, AJ 129, 1993", kind="gz"),
    ]


def extract_atnf() -> None:
    """Pull ``psrcat.db`` out of the release tarball (it sits in a subdirectory)."""
    if config.ATNF_DB.exists():
        return
    tarball = CATALOG_DIR / f"psrcat_pkg.v{config.ATNF_VERSION}.tar.gz"
    with tarfile.open(tarball, mode="r:gz") as tar:
        members = [m for m in tar.getmembers() if m.name.endswith("/psrcat.db")]
        if len(members) != 1:
            raise RuntimeError(f"expected one psrcat.db in {tarball.name}, found {len(members)}")
        config.ATNF_DB.write_bytes(tar.extractfile(members[0]).read())


# --------------------------------------------------------------------------

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--only", choices=GROUPS, action="append",
                        help="fetch only this group (repeatable)")
    parser.add_argument("--refresh", action="store_true", help="re-download cached files")
    args = parser.parse_args(argv)
    groups = args.only or GROUPS

    failures: list[str] = []
    if "plates" in groups:
        print("Plates (DSS) ...")
        failures += fetch_sources(dss_sources(), refresh=args.refresh)
        print("Pan-STARRS stacks ...")
        failures += fetch_sources(ps1_sources(), refresh=args.refresh)
    if "gaia" in groups:
        print("Gaia DR3 ...")
        failures += fetch_sources(gaia_sources(), refresh=args.refresh)
    if "chandra" in groups:
        print("Chandra ...")
        failures += fetch_sources([s for s in chandra_sources() if s.url], refresh=args.refresh)
    if "catalogs" in groups:
        print("Catalogues ...")
        failures += fetch_sources(catalog_sources(), refresh=args.refresh)
        if not failures:
            extract_atnf()

    if failures:
        print("\nSome downloads failed (re-run to retry):")
        for f in failures:
            print("  " + f)
        return 1
    print("All data cached in", DATA_DIR)
    return 0


if __name__ == "__main__":
    sys.exit(main())
