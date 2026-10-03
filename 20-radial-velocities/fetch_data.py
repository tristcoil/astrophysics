"""Download and cache everything the analysis needs.

Six public sources, all free:

1. **Keck/HIRES** velocities of the Lick-Carnegie Exoplanet Survey (Butler et
   al. 2017), plus Tal-Or et al.'s (2019) zero-point-corrected version, from
   the CDS (VizieR) FTP mirror. Fetched whole: ~3 MB compressed.
2. **Lick/Hamilton** velocities from the 25-year Lick planet search (Fischer et
   al. 2014), also from the CDS.
3. **HARPS RVBank** (Trifonov et al. 2020; Ver. 02, Perdelwitz et al. 2024):
   one CSV per star, with activity indices for every spectrum.
4. **ELODIE** (Observatoire de Haute-Provence), the spectrograph 51 Peg b was
   found with: the archive's cross-correlation velocities for 51 Peg, and the
   FITS header of each spectrum, which holds its barycentric time.
5. **NASA Exoplanet Archive**: the confirmed-planet population, every published
   solution for our own planets, and host-star parameters.
6. **TESS** light curves of 55 Cnc from MAST (via lightkurve), for the one
   transiting planet among our targets, and Claret's (2017) limb-darkening
   tables for its star.

Every file lands under ``data/raw/``; a manifest (``data/manifest.json``)
records where each came from, when, its size and SHA-256. Nothing already
cached is fetched again, so after one run the project works offline.

    python fetch_data.py            # fill the cache
    python fetch_data.py --refresh  # re-download everything
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import io
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from config import (
    CLARET, ELODIE_CCF_URL, ELODIE_CITATION, ELODIE_HEADER_URL,
    ELODIE_POLITE_DELAY_S, KECK, KECK_TALOR, LICK, MANIFEST_PATH,
    NASA_POPULATION_COLUMNS, NASA_PS_COLUMNS, NASA_TAP_URL, RAW_DIR,
    RVBANK_CITATION, RVBANK_URL, STARS, TESS_AUTHOR, TESS_DIR, TESS_EXPTIME_S,
    ZENG_BASE_URL, ZENG_MODELS, CdsCatalog,
)

USER_AGENT = "astrophysics-atlas-A3/1.0 (hobbyist research; radial velocities)"
TIMEOUT_S = 180
RETRIES = 5

ELODIE_DIR = RAW_DIR / "elodie"
RVBANK_DIR = RAW_DIR / "rvbank"
NASA_DIR = RAW_DIR / "nasa"
for _d in (ELODIE_DIR, RVBANK_DIR, NASA_DIR, TESS_DIR):
    _d.mkdir(parents=True, exist_ok=True)

NASA_POPULATION_PATH = NASA_DIR / "pscomppars_population.csv"
NASA_PS_PATH = NASA_DIR / "ps_our_planets.csv"
NASA_HOSTS_PATH = NASA_DIR / "stellarhosts_our_stars.csv"
ZENG_DIR = RAW_DIR / "zeng2019"
ZENG_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def _now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text())
    return {"files": {}}


def _record(path: Path, url: str, source: str, **extra) -> None:
    """Add one downloaded file to the provenance manifest."""
    manifest = _load_manifest()
    data = path.read_bytes()
    manifest["files"][str(path.relative_to(RAW_DIR.parent))] = {
        "url": url,
        "source": source,
        "retrieved_utc": _now_utc(),
        "bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        **extra,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, sort_keys=True))


def http_get(url: str, *, params: dict | None = None, user_agent: str = USER_AGENT) -> bytes:
    """GET with retries. DNS lookups on this machine fail intermittently, so a
    ``URLError`` is retried quickly; HTTP errors other than 5xx are not."""
    if params:
        url = f"{url}?{urllib.parse.urlencode(params)}"
    last = None
    for attempt in range(RETRIES):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": user_agent})
            with urllib.request.urlopen(req, timeout=TIMEOUT_S) as resp:
                return resp.read()
        except urllib.error.HTTPError as err:
            if err.code < 500:
                raise
            last = err
        except (urllib.error.URLError, TimeoutError, ConnectionError) as err:
            last = err
        time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"giving up on {url}: {last}")


def _write_atomic(path: Path, data: bytes) -> None:
    """Write via a temporary name so an interrupted download never leaves a
    truncated file where the analysis expects a complete one."""
    tmp = path.with_name(path.name + ".part")
    tmp.write_bytes(data)
    tmp.replace(path)


# ---------------------------------------------------------------------------
# 1-2. VizieR catalogues (Keck, Tal-Or, Lick, Claret)
# ---------------------------------------------------------------------------

def fetch_cds(cat: CdsCatalog, *, refresh: bool = False, quiet: bool = False) -> None:
    for filename in cat.files:
        path = cat.local(filename)
        if path.exists() and not refresh:
            continue
        url = cat.url(filename)
        data = http_get(url)
        if filename.endswith(".gz"):
            gzip.decompress(data)  # raises if the file is not a complete gzip stream
        if len(data) < 500:
            raise RuntimeError(f"{url} returned only {len(data)} bytes")
        _write_atomic(path, data)
        _record(path, url, f"CDS/VizieR {cat.catalog}", citation=cat.citation)
        if not quiet:
            print(f"  {cat.catalog}/{filename}: {len(data) / 1e6:.2f} MB")


# ---------------------------------------------------------------------------
# 3. HARPS RVBank
# ---------------------------------------------------------------------------

def rvbank_path(name: str) -> Path:
    return RVBANK_DIR / f"{name}_HARPS_RVBank_ver02.csv"


def fetch_rvbank(name: str, *, refresh: bool = False, quiet: bool = False) -> None:
    path = rvbank_path(name)
    if path.exists() and not refresh:
        return
    url = RVBANK_URL.format(name=name)
    data = http_get(url)
    header = next(csv.reader(io.StringIO(data[:4000].decode())))
    for column in ("BJD", "RV_mlc_nzp", "e_RV_mlc_nzp", "Halpha", "f_RV"):
        if column not in header:
            raise RuntimeError(f"{url} does not look like an RVBank file (no {column})")
    _write_atomic(path, data)
    _record(path, url, "HARPS RVBank Ver. 02", citation=RVBANK_CITATION)
    if not quiet:
        rows = data.count(b"\n") - 1
        print(f"  RVBank {name}: {rows} spectra")


# ---------------------------------------------------------------------------
# 4. ELODIE
# ---------------------------------------------------------------------------

#: FITS keywords kept from each ELODIE header.
ELODIE_KEYWORDS = ("JDB1", "JDB2", "BERV", "DATE-OBS", "MJD-OBS", "EXPTIME",
                   "IMATYP", "SN", "OBJECT")


def elodie_ccf_path(name: str) -> Path:
    return ELODIE_DIR / f"{name}_ccf.tsv"


def elodie_headers_path(name: str) -> Path:
    return ELODIE_DIR / f"{name}_headers.json.gz"


def elodie_times_path(name: str) -> Path:
    return ELODIE_DIR / f"{name}_header_keywords.csv"


def parse_elodie_ccf(text: str) -> pd.DataFrame:
    """The archive's tab-separated CCF table: comment lines start with ``#``,
    the column-name line with ``$``."""
    rows = [line.rstrip("\n").split("\t") for line in text.splitlines()
            if line.strip() and not line.startswith(("#", "$"))]
    cols = ["objname", "j2000", "simbad", "sophie", "datenuit", "imanum", "type",
            "masque", "imatyp", "exptim", "sn", "vfit", "sigfit", "ampfit"]
    df = pd.DataFrame([r[:len(cols)] for r in rows], columns=cols)
    for c in ("exptim", "sn", "vfit", "sigfit", "ampfit"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    for c in ("datenuit", "imanum", "type", "masque", "imatyp"):
        df[c] = df[c].str.strip()
    return df.drop(columns=["simbad", "sophie"])


def _parse_header_html(raw: bytes) -> tuple[str, dict]:
    """Strip the archive's HTML from a header listing; return text + keywords."""
    text = html.unescape(re.sub(r"<[^>]+>", "", raw.decode("latin-1")))
    lines = [ln for ln in text.splitlines() if re.match(r"^[A-Z][A-Z0-9_-]{0,7}\s", ln)]
    keys = {}
    for ln in lines:
        key = ln[:8].strip()
        if key in ELODIE_KEYWORDS and key not in keys:
            m = re.match(r"^\S+\s+('([^']*)'|(\S+))", ln)
            if m:
                keys[key] = (m.group(2) if m.group(2) is not None else m.group(3)).strip()
    return "\n".join(lines), keys


def fetch_elodie(name: str, *, refresh: bool = False, quiet: bool = False) -> None:
    ccf_path = elodie_ccf_path(name)
    if not ccf_path.exists() or refresh:
        url = ELODIE_CCF_URL.format(name=name)
        data = http_get(url)
        ccf = parse_elodie_ccf(data.decode("latin-1"))
        if len(ccf) < 10 or ccf["vfit"].isna().all():
            raise RuntimeError(f"ELODIE CCF table for {name} looks wrong ({len(ccf)} rows)")
        _write_atomic(ccf_path, data)
        _record(ccf_path, url, "ELODIE archive, table e501 (cross-correlation functions)",
                citation=ELODIE_CITATION)
        if not quiet:
            print(f"  ELODIE {name}: {len(ccf)} cross-correlation velocities")
    ccf = parse_elodie_ccf(ccf_path.read_text(encoding="latin-1"))

    # One FITS header per object spectrum. Cached incrementally so an
    # interrupted run resumes where it stopped.
    hdr_path = elodie_headers_path(name)
    headers = {} if refresh or not hdr_path.exists() else json.loads(gzip.decompress(hdr_path.read_bytes()))
    wanted = [f"{n}/{i}" for n, i in ccf.loc[ccf["type"] == "obj", ["datenuit", "imanum"]].itertuples(index=False)]
    missing = [w for w in dict.fromkeys(wanted) if w not in headers]
    if missing and not quiet:
        print(f"  ELODIE {name}: fetching {len(missing)} FITS headers "
              f"(~{len(missing) * (0.6 + ELODIE_POLITE_DELAY_S):.0f} s)")
    for k, spec in enumerate(missing):
        night, imanum = spec.split("/")
        url = ELODIE_HEADER_URL.format(night=night, imanum=imanum)
        text, keys = _parse_header_html(http_get(url))
        if "JDB1" not in keys:
            raise RuntimeError(f"no JDB in the header of {spec}")
        headers[spec] = text
        if (k + 1) % 25 == 0 or k + 1 == len(missing):
            _write_atomic(hdr_path, gzip.compress(json.dumps(headers).encode()))
        time.sleep(ELODIE_POLITE_DELAY_S)
    if missing:
        _record(hdr_path, ELODIE_HEADER_URL, "ELODIE archive, table e500 (spectra) FITS headers",
                citation=ELODIE_CITATION, n_headers=len(headers))

    # A tidy table of the keywords the analysis needs, rebuilt from the cache.
    rows = []
    for spec, text in headers.items():
        _, keys = _parse_header_html(text.encode("latin-1"))
        night, imanum = spec.split("/")
        rows.append({"datenuit": night, "imanum": imanum, **keys})
    pd.DataFrame(rows).to_csv(elodie_times_path(name), index=False)


# ---------------------------------------------------------------------------
# 5. NASA Exoplanet Archive
# ---------------------------------------------------------------------------

def _tap(query: str) -> tuple[pd.DataFrame, str]:
    data = http_get(NASA_TAP_URL, params={"query": query, "format": "csv"})
    text = data.decode()
    if text.lstrip().startswith(("<?xml", "<VOTABLE", "ERROR")):
        raise RuntimeError(f"archive returned an error document:\n{text[:300]}")
    return pd.read_csv(io.StringIO(text), low_memory=False), text


def fetch_nasa(*, refresh: bool = False, quiet: bool = False) -> None:
    hosts = ",".join(f"'{s.archive_host}'" for s in STARS.values())
    jobs = [
        (NASA_POPULATION_PATH, f"select {','.join(NASA_POPULATION_COLUMNS)} from pscomppars",
         "pscomppars (every confirmed planet, one row each)", 4000),
        (NASA_PS_PATH, f"select {','.join(NASA_PS_COLUMNS)} from ps where hostname in ({hosts})",
         "ps (every published solution for our planets)", 20),
        (NASA_HOSTS_PATH,
         "select hostname,st_refname,st_teff,st_tefferr1,st_logg,st_loggerr1,st_met,"
         "st_meterr1,st_mass,st_masserr1,st_masserr2,st_rad,st_raderr1,st_raderr2,"
         f"st_rotp,st_rotperr1 from stellarhosts where hostname in ({hosts})",
         "stellarhosts (every published parameter set for our stars)", 4),
    ]
    for path, query, label, min_rows in jobs:
        if path.exists() and not refresh:
            continue
        table, text = _tap(query)
        if len(table) < min_rows:
            raise RuntimeError(f"only {len(table)} rows for {label}; refusing to cache")
        _write_atomic(path, text.encode())
        _record(path, f"{NASA_TAP_URL}?query={query}", "NASA Exoplanet Archive: " + label,
                query=query, rows=int(len(table)))
        if not quiet:
            print(f"  NASA archive {label.split(' ')[0]}: {len(table):,} rows")


# ---------------------------------------------------------------------------
# 6. TESS light curves
# ---------------------------------------------------------------------------

def tess_files() -> list[Path]:
    return sorted(TESS_DIR.rglob("*_lc.fits"))


def fetch_tess(*, refresh: bool = False, quiet: bool = False) -> None:
    targets = [s for s in STARS.values() if s.tess_target]
    if tess_files() and not refresh:
        return
    import warnings

    import lightkurve as lk

    for star in targets:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            search = lk.search_lightcurve(star.tess_target, mission="TESS",
                                          author=TESS_AUTHOR, exptime=TESS_EXPTIME_S)
            if len(search) == 0:
                raise RuntimeError(f"no {TESS_AUTHOR} {TESS_EXPTIME_S}-s light curves for {star.name}")
            search.download_all(download_dir=str(TESS_DIR))
        for path in tess_files():
            _record(path, "https://mast.stsci.edu (via lightkurve)", "TESS SPOC light curve",
                    target=star.tess_target)
        if not quiet:
            sectors = ", ".join(str(m).split()[-1] for m in search.table["mission"])
            print(f"  TESS {star.name}: sectors {sectors}")


# ---------------------------------------------------------------------------
# 7. Planet-interior models
# ---------------------------------------------------------------------------

def fetch_zeng(*, refresh: bool = False, quiet: bool = False) -> None:
    import numpy as np
    for key, (filename, _) in ZENG_MODELS.items():
        path = ZENG_DIR / filename
        if path.exists() and not refresh:
            continue
        url = ZENG_BASE_URL + filename
        data = http_get(url)
        table = np.loadtxt(io.StringIO(data.decode()))
        if table.ndim != 2 or table.shape[1] != 2 or len(table) < 5:
            raise RuntimeError(f"{url} is not a two-column mass-radius table")
        _write_atomic(path, data)
        _record(path, url, "Zeng et al. 2019 mass-radius model", citation="Zeng et al. 2019, PNAS 116, 9723")
        if not quiet:
            print(f"  Zeng et al. model: {filename}")


def load_zeng(key: str):
    import numpy as np
    m, r = np.loadtxt(ZENG_DIR / ZENG_MODELS[key][0], unpack=True)
    order = np.argsort(m)
    return m[order], r[order]


# ---------------------------------------------------------------------------

def ensure_data(*, refresh: bool = False, quiet: bool = False) -> None:
    """Fetch anything missing from the cache; a no-op when it is complete."""
    for cat in (KECK, KECK_TALOR, LICK, CLARET):
        fetch_cds(cat, refresh=refresh, quiet=quiet)
    for star in STARS.values():
        if star.rvbank:
            fetch_rvbank(star.rvbank, refresh=refresh, quiet=quiet)
        if star.elodie:
            fetch_elodie(star.elodie, refresh=refresh, quiet=quiet)
    fetch_nasa(refresh=refresh, quiet=quiet)
    fetch_tess(refresh=refresh, quiet=quiet)
    fetch_zeng(refresh=refresh, quiet=quiet)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true", help="re-download everything")
    args = parser.parse_args(argv)
    ensure_data(refresh=args.refresh)
    manifest = _load_manifest()
    total = sum(f["bytes"] for f in manifest["files"].values())
    print(f"\nCache complete: {len(manifest['files'])} files, {total / 1e6:.1f} MB "
          f"(provenance in {MANIFEST_PATH.relative_to(RAW_DIR.parent.parent)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
