"""Download and cache every input file, and record where each one came from.

Three archives, all small (~13 MB in total):

1. **WDC-SILSO** -- the international sunspot number (version 2.0), its official
   table of cycle minima and maxima, and four reconstructions of the sunspot group
   number that reach back to 1610. SILSO publishes SHA-256 hashes of its sunspot
   number files; each download is checked against them, so a truncated or
   half-updated file cannot slip through.
2. **The RGO/USAF/NOAA sunspot-group database** (Hathaway, solarcyclescience.com):
   one text file per year since 1874, one line per sunspot group per day.
3. **NOAA SWPC** -- the 2019 NOAA/NASA panel forecast for Cycle 25, with its stated
   uncertainty range, and SWPC's current forecast.

Everything lands in ``data/`` and is never downloaded twice (``--refresh`` forces
it). ``data/manifest.json`` records the URL, size, SHA-256 and retrieval time of
each file, so every number in the README can be traced to the exact bytes it came
from. SILSO updates its files monthly (recent months are provisional), so a later
download can shift the last few months of the record slightly.

Run directly to populate the cache::

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

from config import (
    AR_DATABASE_URL,
    AR_DIR,
    AR_EXTRA_FILES,
    AR_FILE_PATTERN,
    AR_FIRST_YEAR,
    DATA_DIR,
    GROUP_NUMBER_SERIES,
    MANIFEST_PATH,
    SILSO_CYCLE_TABLE,
    SILSO_DATA_URL,
    SILSO_DOI,
    SILSO_HASH_MANIFEST,
    SILSO_SERIES,
    SWPC_FILES,
    SWPC_URL,
)

SILSO_DIR = DATA_DIR / "silso"
SWPC_DIR = DATA_DIR / "swpc"
#: A pause between the ~150 yearly group files: this is a one-person research site.
POLITE_DELAY_S = 0.2
USER_AGENT = "astrophysics-atlas-F1/1.0 (hobby research; python-requests)"


def _download(url: str, dest: Path) -> None:
    """Stream ``url`` to ``dest`` via a temporary file, so an interrupted download
    never leaves a truncated file that later runs would mistake for the cache."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=60, headers={"User-Agent": USER_AGENT}) as r:
        r.raise_for_status()
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 16):
                fh.write(chunk)
    tmp.replace(dest)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 16), b""):
            h.update(block)
    return h.hexdigest()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _record(path: Path, url: str, retrieved: str) -> dict:
    return {"url": url, "bytes": path.stat().st_size, "sha256": _sha256(path),
            "retrieved_utc": retrieved}


def _fetch(url: str, dest: Path, refresh: bool, quiet: bool) -> bool:
    """Download unless cached; return True if a download happened."""
    if dest.exists() and not refresh:
        return False
    if not quiet:
        print(f"  downloading {url}")
    _download(url, dest)
    return True


# --- WDC-SILSO --------------------------------------------------------------------------


def _silso_published_hashes(path: Path) -> dict[str, str]:
    """SILSO's hash list: '<sha256>  /their/archive/path/<file>' per line."""
    out = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) == 2:
            out[Path(parts[1]).name] = parts[0]
    return out


def fetch_silso(*, refresh: bool = False, quiet: bool = False) -> dict:
    """Sunspot numbers (hash-verified), the cycle table and the group-number series."""
    SILSO_DIR.mkdir(parents=True, exist_ok=True)
    entries = {}
    names = list(SILSO_SERIES.values())
    downloaded = False
    for name in names:
        downloaded |= _fetch(SILSO_DATA_URL + name, SILSO_DIR / name, refresh, quiet)
    # The hash list must be fetched in the same session as the files it vouches for.
    hash_path = SILSO_DIR / SILSO_HASH_MANIFEST
    _fetch(SILSO_DATA_URL + SILSO_HASH_MANIFEST, hash_path, refresh or downloaded, quiet)
    published = _silso_published_hashes(hash_path)

    retrieved = _now()
    for name in names:
        path = SILSO_DIR / name
        rec = _record(path, SILSO_DATA_URL + name, retrieved)
        expected = published.get(name)
        rec["silso_sha256_verified"] = expected == rec["sha256"]
        if expected is not None and expected != rec["sha256"]:
            raise RuntimeError(
                f"{name}: SHA-256 does not match SILSO's published hash -- the file "
                "may have been caught mid-update; re-run with --refresh")
        entries[name] = rec

    extra = [SILSO_CYCLE_TABLE] + [s["file"] for s in GROUP_NUMBER_SERIES.values()]
    for rel in extra:
        dest = SILSO_DIR / Path(rel).name
        _fetch(SILSO_DATA_URL + rel, dest, refresh, quiet)
        entries[dest.name] = _record(dest, SILSO_DATA_URL + rel, retrieved)
    return entries


# --- The RGO / USAF / NOAA sunspot-group database -----------------------------------------


def fetch_groups(*, refresh: bool = False, quiet: bool = False) -> dict:
    """One file per year from 1874 through the current year, plus documentation."""
    entries = {}
    this_year = datetime.now(timezone.utc).year
    for year in range(AR_FIRST_YEAR, this_year + 1):
        name = AR_FILE_PATTERN.format(year=year)
        dest = AR_DIR / name
        try:
            if _fetch(AR_DATABASE_URL + name, dest, refresh, quiet):
                time.sleep(POLITE_DELAY_S)
        except requests.HTTPError as err:
            # Early in January the new year's file may not exist yet.
            if year == this_year and err.response is not None and err.response.status_code == 404:
                continue
            raise
        entries[name] = _record(dest, AR_DATABASE_URL + name, _now())
    for name in AR_EXTRA_FILES:
        dest = AR_DIR / name
        _fetch(AR_DATABASE_URL + name, dest, refresh, quiet)
        entries[name] = _record(dest, AR_DATABASE_URL + name, _now())
    return entries


# --- NOAA SWPC -------------------------------------------------------------------------------


def fetch_swpc(*, refresh: bool = False, quiet: bool = False) -> dict:
    SWPC_DIR.mkdir(parents=True, exist_ok=True)
    entries = {}
    for name in SWPC_FILES.values():
        dest = SWPC_DIR / name
        _fetch(SWPC_URL + name, dest, refresh, quiet)
        json.loads(dest.read_text())  # fail now, not mid-analysis, if it is not JSON
        entries[name] = _record(dest, SWPC_URL + name, _now())
    return entries


# --- Everything ------------------------------------------------------------------------------


def _last_silso_day(path: Path) -> str:
    last = path.read_text().strip().splitlines()[-1].split(";")
    return f"{last[0]}-{last[1]}-{last[2]}"


def _last_group_day() -> str:
    files = sorted(AR_DIR.glob("g[0-9][0-9][0-9][0-9].txt"))
    lines = [ln for ln in files[-1].read_text().splitlines() if ln.strip()]
    ln = lines[-1]
    return f"{ln[0:4]}-{int(ln[4:6]):02d}-{int(ln[6:8]):02d}"


def fetch_all(*, refresh: bool = False, quiet: bool = False) -> dict:
    """Make sure every input is cached; write and return the manifest."""
    if MANIFEST_PATH.exists() and not refresh:
        manifest = load_manifest()
        missing = [n for sec in ("silso", "sunspot_groups", "swpc")
                   for n in manifest.get(sec, {})
                   if not (DATA_DIR / manifest["dirs"][sec] / n).exists()]
        if not missing:
            if not quiet:
                print(f"all inputs already cached (retrieved {manifest['retrieved_utc'][:10]})")
            return manifest

    if not quiet:
        print("WDC-SILSO sunspot numbers and group numbers")
    silso = fetch_silso(refresh=refresh, quiet=quiet)
    if not quiet:
        print("RGO/USAF/NOAA sunspot-group database")
    groups = fetch_groups(refresh=refresh, quiet=quiet)
    if not quiet:
        print("NOAA SWPC forecasts")
    swpc = fetch_swpc(refresh=refresh, quiet=quiet)

    manifest = {
        "retrieved_utc": _now(),
        "dirs": {"silso": "silso", "sunspot_groups": AR_DIR.name, "swpc": "swpc"},
        "credits": {
            "silso": f"WDC-SILSO, Royal Observatory of Belgium, Brussels ({SILSO_DOI})",
            "sunspot_groups": ("Royal Greenwich Observatory Photoheliographic Results 1874-1976 "
                               "and USAF/NOAA SOON data 1977-, compiled by D. H. Hathaway "
                               "(solarcyclescience.com)"),
            "swpc": "NOAA Space Weather Prediction Center",
        },
        "last_day": {
            "silso_daily": _last_silso_day(SILSO_DIR / SILSO_SERIES["daily"]),
            "sunspot_groups": _last_group_day(),
        },
        "silso": silso,
        "sunspot_groups": groups,
        "swpc": swpc,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n")
    if not quiet:
        total = sum(e["bytes"] for sec in (silso, groups, swpc) for e in sec.values())
        n_ok = sum(e.get("silso_sha256_verified", False) for e in silso.values())
        print(f"cached {len(silso) + len(groups) + len(swpc)} files, {total / 1e6:.1f} MB; "
              f"{n_ok} SILSO files match SILSO's published SHA-256; "
              f"sunspot record to {manifest['last_day']['silso_daily']}, "
              f"group record to {manifest['last_day']['sunspot_groups']}")
    return manifest


def load_manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text())


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true",
                        help="download everything again even if cached")
    args = parser.parse_args(argv)
    fetch_all(refresh=args.refresh)
    return 0


if __name__ == "__main__":
    sys.exit(main())
