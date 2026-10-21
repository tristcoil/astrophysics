"""Download and cache every input the project uses. Safe to re-run: files already
on disk are kept (pass --refresh to fetch again).

What is fetched
---------------
1. JPL's DE440 planetary ephemeris (de440.bsp, 120 MB) and its constants
   (header.440, 23 kB) -- initial conditions, masses and the answer key.
2. JPL Horizons: Mercury's heliocentric state from DE441 (DE440's long-span
   twin) at 241 epochs from 1000 BC to AD 5000, and its osculating elements at
   109 epochs from 1560 to 2640 -- the spec's "check against Horizons' own
   elements", and a six-millennium test of our own integration.
3. Baroch et al. (2021)'s LaTeX source (arXiv:2103.03140v2, 5 MB) and their
   TESS eclipse times from VizieR: eclipsing binaries whose orbits precess
   mostly through general relativity.
4. JPL's small-body database: orbits of asteroids 1-1467, for the 343 whose
   masses DE440 includes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from config import (DE440_BSP, DE440_BSP_URL, DE440_HEADER, DE440_HEADER_URL, HORIZONS_API,
                    J2000_JD, MANIFEST_PATH, RAW_DIR, YEAR_DAYS)

USER_AGENT = "astrophysics-atlas-N1/1.0 (hobby project; python urllib)"
HORIZONS_DIR = RAW_DIR / "horizons"
HORIZONS_DIR.mkdir(parents=True, exist_ok=True)

#: Horizons answers a TLIST of at most 80 epochs per request.
TLIST_CHUNK = 80

#: Mercury every 25 Julian years, 3,000 years either side of J2000 (DE441).
LONG_EPOCHS = [J2000_JD + 25.0 * k * YEAR_DAYS for k in range(-120, 121)]
#: Mercury every 10 Julian years across DE440's span, for the element check.
DE440_SPAN_EPOCHS = [J2000_JD + 10.0 * k * YEAR_DAYS for k in range(-44, 65)]


# --------------------------------------------------------------------------
# Bookkeeping
# --------------------------------------------------------------------------

def _manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text()) if MANIFEST_PATH.exists() else {}


def _record(key: str, path: Path, url: str, **extra) -> None:
    m = _manifest()
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    m[key] = {"file": str(path.relative_to(MANIFEST_PATH.parent)), "url": url,
              "bytes": path.stat().st_size, "sha256": h.hexdigest(),
              "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), **extra}
    MANIFEST_PATH.write_text(json.dumps(m, indent=1, sort_keys=True))


def download(url: str, path: Path, *, key: str, refresh: bool = False, retries: int = 4,
             timeout: float = 600, data: bytes | None = None, **extra) -> Path:
    """Fetch ``url`` to ``path`` atomically (via a .part file) and log it."""
    if path.exists() and not refresh:
        if key not in _manifest():
            _record(key, path, url, **extra)
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as r, open(tmp, "wb") as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
            tmp.replace(path)
            _record(key, path, url, **extra)
            return path
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last = exc
            time.sleep(3.0 * (attempt + 1))
    raise RuntimeError(f"could not download {url}: {last}")


# --------------------------------------------------------------------------
# DE440
# --------------------------------------------------------------------------

def fetch_de440(refresh: bool = False) -> None:
    download(DE440_HEADER_URL, DE440_HEADER, key="de440_header", refresh=refresh,
             purpose="DE440 constants: GM values, AU, c, solar J2 and spin")
    download(DE440_BSP_URL, DE440_BSP, key="de440_bsp", refresh=refresh,
             purpose="DE440 Chebyshev ephemeris 1550-2650 (SPICE SPK)")
    size = DE440_BSP.stat().st_size
    if size != 119_799_808:
        print(f"  warning: de440.bsp is {size} bytes (expected 119,799,808)")


# --------------------------------------------------------------------------
# Horizons (DE441)
# --------------------------------------------------------------------------

def _horizons_params(kind: str, jds: list[float]) -> dict:
    p = {"format": "text", "COMMAND": "'199'", "OBJ_DATA": "'NO'", "MAKE_EPHEM": "'YES'",
         "EPHEM_TYPE": f"'{kind}'", "CENTER": "'500@10'", "REF_PLANE": "'ECLIPTIC'",
         "REF_SYSTEM": "'ICRF'", "OUT_UNITS": "'AU-D'", "CSV_FORMAT": "'YES'",
         "TIME_TYPE": "'TDB'", "TLIST_TYPE": "'JD'",
         "TLIST": "'" + " ".join(f"{j:.6f}" for j in jds) + "'"}
    if kind == "VECTORS":
        p.update({"VEC_TABLE": "'2'", "VEC_LABELS": "'NO'"})
    else:
        p.update({"ELM_LABELS": "'NO'"})
    return p


def horizons_files(label: str) -> list[Path]:
    return sorted(HORIZONS_DIR.glob(f"{label}_*.txt"))


def fetch_horizons_series(label: str, kind: str, jds: list[float], refresh: bool = False) -> None:
    for i in range(0, len(jds), TLIST_CHUNK):
        chunk = jds[i:i + TLIST_CHUNK]
        path = HORIZONS_DIR / f"{label}_{i // TLIST_CHUNK:02d}.txt"
        if path.exists() and not refresh:
            continue
        body = urllib.parse.urlencode(_horizons_params(kind, chunk)).encode()
        download(HORIZONS_API, path, key=f"horizons_{label}_{i // TLIST_CHUNK:02d}",
                 refresh=True, data=body, timeout=180,
                 purpose=f"Mercury heliocentric {kind.lower()} (DE441), {len(chunk)} epochs")
        text = path.read_text()
        n = text.count("\n", text.index("$$SOE"), text.index("$$EOE")) - 1 if "$$SOE" in text else 0
        if n != len(chunk):
            path.unlink()
            raise RuntimeError(f"Horizons returned {n} rows for {len(chunk)} epochs ({path.name})")


def fetch_horizons(refresh: bool = False) -> None:
    fetch_horizons_series("mercury_vectors_long", "VECTORS", LONG_EPOCHS, refresh)
    fetch_horizons_series("mercury_vectors_span", "VECTORS", DE440_SPAN_EPOCHS, refresh)
    fetch_horizons_series("mercury_elements_span", "ELEMENTS", DE440_SPAN_EPOCHS, refresh)


# --------------------------------------------------------------------------
# Literature: eclipsing binaries
# --------------------------------------------------------------------------

ARXIV_DIR = RAW_DIR / "arxiv"
VIZIER_DIR = RAW_DIR / "vizier"


def fetch_literature(refresh: bool = False) -> None:
    """Baroch et al. 2021 (A&A 649, A64): the LaTeX source holds its tables of
    system properties, apsidal rates and eclipse-timing differences; VizieR
    J/A+A/649/A64 holds its TESS eclipse times."""
    path = ARXIV_DIR / "2103.03140v2.tar.gz"
    download("https://arxiv.org/e-print/2103.03140v2", path, key="arxiv_2103.03140v2",
             refresh=refresh, purpose="Baroch et al. 2021 LaTeX source (tables 1, 4-9)")
    import tarfile
    if not tarfile.is_tarfile(path):
        raise RuntimeError(f"{path} is not a tar archive")
    base = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/649/A64/"
    for name in ("ReadMe", "table3.dat", "stars.dat"):
        download(base + name, VIZIER_DIR / "J_A+A_649_A64" / name,
                 key=f"vizier_J_A+A_649_A64_{name}", refresh=refresh,
                 purpose="Baroch et al. 2021 TESS eclipse times (VizieR)")


SBDB_DIR = RAW_DIR / "sbdb"
SBDB_ASTEROIDS = SBDB_DIR / "numbered_1_1467.json"


def fetch_asteroid_orbits(refresh: bool = False) -> None:
    """Semi-major axes of asteroids 1-1467 (JPL SBDB), for the 343 whose masses
    DE440 includes: their combined pull on Mercury's perihelion."""
    q = urllib.parse.urlencode({"fields": "spkid,pdes,a,e,i", "sb-ns": "n",
                                "sb-cdata": '{"AND":["spkid|RG|20000001|20001467"]}'})
    download("https://ssd-api.jpl.nasa.gov/sbdb_query.api?" + q, SBDB_ASTEROIDS,
             key="sbdb_numbered_1_1467", refresh=refresh,
             purpose="asteroid orbits for DE440's 343 massive asteroids")
    rows = json.loads(SBDB_ASTEROIDS.read_text())["data"]
    if len(rows) != 1467:
        raise RuntimeError(f"SBDB returned {len(rows)} asteroids, expected 1467")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refresh", action="store_true", help="download again even if cached")
    args = ap.parse_args(argv)
    steps = [("DE440 ephemeris and constants", fetch_de440),
             ("JPL Horizons (DE441) checks", fetch_horizons),
             ("Literature tables (eclipsing binaries)", fetch_literature),
             ("Asteroid orbits (JPL SBDB)", fetch_asteroid_orbits)]
    for label, fn in steps:
        print(f"* {label}", flush=True)
        fn(refresh=args.refresh)
    print("All inputs cached under", RAW_DIR)
    return 0


if __name__ == "__main__":
    sys.exit(main())
