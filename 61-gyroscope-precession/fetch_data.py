"""Download and cache every input this project uses.

    python fetch_data.py              # everything
    python fetch_data.py orbits eop   # only some sources

1. **ILRS weekly precise orbits** of LAGEOS 1, LARES 2 and LAGEOS 2 from the
   EUROLAS Data Center (EDC, anonymous FTP): every analysis centre (AC) fits a
   7-day orbit to the satellites' laser ranges each week and publishes it as an
   SP3 file (Earth-fixed SLRF2020 frame, 2-min positions and velocities). We keep
   the latest version of each AC's file for every week from July 2022 on
   (~0.24 MB each).
2. **IERS Earth orientation parameters** (EOP 20 C04) to rotate the orbits into
   the inertial (GCRS) frame.
3. **Earth's gravity field** (GOCO06s, GGM05S, EIGEN-6S4 v2 from ICGEM) and the
   **FES2004 ocean tides** in the IERS Conventions' geopotential form.
4. **JPL DE440** for the Sun and Moon (hard link to N1's copy if present).
5. **IERS Conventions 2010** chapters 5, 6, 7 (the model recipes we follow).
6. **CelesTrak's current two-line elements** of the satellites, to compare
   them with the precise orbits.

Files land in ``data/raw``; nothing is re-downloaded once present, so after one
run the analysis works offline. ``data/manifest.json`` records each file's URL,
size, SHA-256 and download date.
"""

from __future__ import annotations

import concurrent.futures as cf
import datetime as dt
import hashlib
import json
import os
import re
import sys
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path

import config as C

UA = "astro-atlas-fetch/1.0 (python-urllib; hobby research project)"
UA_BROWSER = "Mozilla/5.0 (X11; Linux x86_64) astro-atlas hobby project"


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

class NotFound(Exception):
    pass


def fetch(url: str, *, timeout: float = 120.0, retries: int = 8, ua: str = UA) -> bytes:
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": ua})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as exc:
            if exc.code in (403, 404, 410):
                raise NotFound(f"HTTP {exc.code} for {url}") from exc
            last = exc
        except urllib.error.URLError as exc:
            if "550" in str(exc.reason):          # FTP: no such file or directory
                raise NotFound(f"FTP 550 for {url}") from exc
            if "421" in str(exc.reason):          # EDC: too many connections from one address
                time.sleep(20.0 * (attempt + 1))
            last = exc
        except (TimeoutError, ConnectionError, OSError) as exc:
            last = exc
        time.sleep(3.0 * (attempt + 1))
    raise RuntimeError(f"could not fetch {url}: {last}")


def save(dest: Path, payload: bytes) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    tmp.write_bytes(payload)
    tmp.replace(dest)
    return dest


_manifest_lock = threading.Lock()


def _read_manifest() -> dict:
    for _ in range(20):                     # another process may be replacing it right now
        try:
            return json.loads(C.MANIFEST_PATH.read_text()) if C.MANIFEST_PATH.exists() else {"files": {}}
        except json.JSONDecodeError:
            time.sleep(0.2)
    raise RuntimeError("manifest.json unreadable")


def record(dest: Path, url: str, payload: bytes, note: str = "") -> None:
    with _manifest_lock:
        man = _read_manifest()
        man["files"][str(dest.relative_to(C.PROJECT_DIR))] = {
            "url": url, "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
            "downloaded": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            **({"note": note} if note else {}),
        }
        tmp = C.MANIFEST_PATH.with_suffix(f".tmp{os.getpid()}")
        tmp.write_text(json.dumps(man, indent=1, sort_keys=True))
        tmp.replace(C.MANIFEST_PATH)          # atomic: readers never see a half-written file


def fetch_file(url: str, dest: Path, *, note: str = "", ua: str = UA, quiet: bool = False) -> Path:
    if dest.exists():
        return dest
    payload = fetch(url, ua=ua)
    save(dest, payload)
    record(dest, url, payload, note)
    if not quiet:
        print(f"  {dest.relative_to(C.PROJECT_DIR)}  ({len(payload) / 1e6:.2f} MB)")
    return dest


def ftp_listing(url: str) -> list[str]:
    """Names in an FTP directory listing (urllib returns 'ls -l' style lines)."""
    text = fetch(url, timeout=90.0).decode("latin-1")
    return [line.split()[-1] for line in text.splitlines() if line.strip()]


# ---------------------------------------------------------------------------
# 1. ILRS weekly orbits (EDC)
# ---------------------------------------------------------------------------

_SP3_RE = re.compile(r"^(?P<ac>[a-z]+)\.orb\.(?P<sat>[a-z0-9]+)\.(?P<week>\d{6})\.v(?P<ver>\d+)\.sp3\.gz$")


def orbit_inventory(sat: str, refresh: bool = False) -> list[dict]:
    """Every weekly SP3 file on EDC for one satellite since 2022 (cached as JSON)."""
    path = C.SP3_DIR / f"inventory_{sat}.json"
    if path.exists() and not refresh:
        return json.loads(path.read_text())
    base = "ftp://edc.dgfi.tum.de/pub/slr/products/orbits/TRF/{sat}/"
    years = [y for y in ftp_listing(base.format(sat=sat)) if y.isdigit() and 2022 <= int(y) <= 2026]
    jobs = []
    for year in years:
        for week in ftp_listing(base.format(sat=sat) + f"{year}/"):
            if week.isdigit():
                jobs.append((year, week))

    def list_week(job):
        year, week = job
        url = C.EDC_ORBITS.format(sat=sat, year=year, week=week)
        out = []
        for name in ftp_listing(url):
            m = _SP3_RE.match(name)
            if m:
                out.append({"ac": m["ac"], "week": week, "version": int(m["ver"]),
                            "name": name, "url": url + name})
        return out

    files = []
    with cf.ThreadPoolExecutor(max_workers=2) as pool:      # EDC refuses many connections
        for chunk in pool.map(list_week, jobs):
            files.extend(chunk)
    files.sort(key=lambda f: (f["week"], f["ac"], f["version"]))
    path.write_text(json.dumps(files, indent=0))
    print(f"  {sat}: {len(jobs)} weeks, {len(files)} SP3 files listed")
    return files


def latest_versions(files: list[dict], acs=C.ACS) -> dict:
    """{(ac, week): file} keeping each AC's highest version number per week."""
    best = {}
    for f in files:
        if f["ac"] not in acs:
            continue
        key = (f["ac"], f["week"])
        if key not in best or f["version"] > best[key]["version"]:
            best[key] = f
    return best


def fetch_orbits(sats=("lageos1", "lares2"), acs=C.ACS, workers: int = 2) -> None:
    """Download centre by centre (in the order of ``acs``), both satellites each time."""
    best = {sat: latest_versions(orbit_inventory(sat), acs) for sat in sats}

    def get(job):
        url, dest = job
        fetch_file(url, dest, quiet=True)
        return dest

    for ac in acs:
        todo = []
        for sat in sats:
            for (a, week), f in sorted(best[sat].items()):
                d = dt.date(2000 + int(week[:2]), int(week[2:4]), int(week[4:]))
                if a != ac or not (C.ANALYSIS_START <= d <= C.ANALYSIS_END + dt.timedelta(days=7)):
                    continue
                dest = C.SP3_DIR / sat / ac / f["name"]
                if not dest.exists():
                    todo.append((f["url"], dest))
        print(f"  {ac}: {len(todo)} files to download")
        with cf.ThreadPoolExecutor(max_workers=workers) as pool:
            for i, _ in enumerate(pool.map(get, todo), 1):
                if i % 100 == 0:
                    print(f"    {i}/{len(todo)}")


# ---------------------------------------------------------------------------
# 2-6. Small files
# ---------------------------------------------------------------------------

def fetch_eop() -> None:
    fetch_file(C.EOP_C04_URL, C.EOP_DIR / "eopc04.1962-now", ua=UA_BROWSER,
               note="IERS EOP 20 C04 (polar motion, UT1-UTC, dX, dY), daily, Paris Observatory")
    fetch_file(C.EOP_FINALS_URL, C.EOP_DIR / "finals2000A.all", ua=UA_BROWSER,
               note="IERS rapid service EOP (Bulletin A), for dates after the end of C04")


def fetch_tides() -> None:
    fetch_file(C.FES2004_URL, C.TIDE_DIR / "fes2004_Cnm-Snm.dat", ua=UA_BROWSER,
               note="IERS Conventions 2010 ch.6: FES2004 ocean tide geopotential coefficients")


def fetch_gravity() -> None:
    for name, key in C.GRAVITY_MODELS.items():
        url = C.ICGEM_URL.format(key=key, name=urllib.request.quote(name))
        fetch_file(url, C.GRAVITY_DIR / f"{name.replace(' ', '_')}.gfc", ua=UA_BROWSER,
                   note=f"ICGEM static gravity field {name}")


def fetch_ephemeris() -> None:
    dest = C.EPHEM_DIR / "de440.bsp"
    if dest.exists():
        return
    if C.DE440_SHARED.exists():
        os.link(C.DE440_SHARED, dest)          # hard link: no second copy on disk
        print(f"  {dest.relative_to(C.PROJECT_DIR)}  (hard link to N1's copy)")
        return
    fetch_file(C.DE440S_URL, C.EPHEM_DIR / "de440s.bsp", note="JPL DE440 short (1849-2150)")


def fetch_conventions() -> None:
    for ch in (5, 6, 7):
        fetch_file(C.IERS_CONVENTIONS_URL.format(ch=ch), C.LIT_DIR / f"iers2010_ch{ch}.pdf",
                   ua=UA_BROWSER, note=f"IERS Conventions (2010), chapter {ch}")


def fetch_c20() -> None:
    fetch_file(C.CSR_C20_URL, C.GRAVITY_DIR / "C20_RL06.txt", ua=UA_BROWSER,
               note="UT/CSR monthly C20 from SLR, RL06")


def fetch_tle() -> None:
    """Current element sets; CelesTrak keeps no public history (Space-Track needs an account)."""
    stamp = dt.date.today().isoformat()
    for sat, info in C.SATELLITES.items():
        dest = C.TLE_DIR / f"{sat}_{stamp}.tle"
        if not list(C.TLE_DIR.glob(f"{sat}_*.tle")):
            fetch_file(C.CELESTRAK_GP.format(norad=info["norad"]), dest, ua=UA_BROWSER,
                       note=f"CelesTrak GP element set for {info['name']} on {stamp}")


def reconcile_orbits() -> int:
    """Record any cached SP3 file missing from the manifest (e.g. after an interrupted run)."""
    man = _read_manifest()
    known = set(man["files"])
    urls = {}
    for inv in C.SP3_DIR.glob("inventory_*.json"):
        for f in json.loads(inv.read_text()):
            urls[f["name"]] = f["url"]
    added = 0
    for path in sorted(C.SP3_DIR.glob("*/*/*.sp3.gz")):
        key = str(path.relative_to(C.PROJECT_DIR))
        if key not in known and path.name in urls:
            record(path, urls[path.name], path.read_bytes(), note="recorded after download")
            added += 1
    return added


SOURCES = {
    "eop": fetch_eop, "tides": fetch_tides, "gravity": fetch_gravity,
    "ephemeris": fetch_ephemeris, "conventions": fetch_conventions, "tle": fetch_tle,
    "c20": fetch_c20,
    "orbits": fetch_orbits,
}


def main(names: list[str]) -> None:
    for name in names or list(SOURCES):
        print(f"== {name}")
        SOURCES[name]()


if __name__ == "__main__":
    main(sys.argv[1:])
