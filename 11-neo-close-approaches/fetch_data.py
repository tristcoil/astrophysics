"""Download and cache everything the analysis needs.

Everything comes from NASA/JPL's Solar System Dynamics group and its Center for
Near-Earth Object Studies (CNEOS), through their free public APIs (no key):

1. **Close approaches** (``cad.api``): every approach of a known near-Earth
   asteroid to within 0.05 au of Earth, 1900-2200 -- the past computed
   backwards and the future forwards from each asteroid's fitted orbit.
2. **The near-Earth asteroid catalogue** (``sbdb_query.api``): one row per known
   NEA, with the date it was first observed.
3. **Sentry**, the impact-monitoring system (``sentry.api``): every object with
   a non-zero chance of hitting Earth in the next century or so, every one of
   its possible impact dates, and the list of objects already removed.
4. **Fireballs** (``fireball.api``): bright meteors -- small asteroids actually
   hitting the atmosphere -- detected by US government sensors since 1988.
5. **Apophis**: its database record (``sbdb.api``) and, from the Horizons
   ephemeris system, its path past Earth in April 2029, the uncertainty of that
   path, and its orbit around the Sun before and after the flyby.

Each response is checked before it is written (these services report many
errors as an ordinary response), written atomically, and recorded in
``data/manifest.json`` with its URL, parameters and retrieval time. Nothing
already cached is fetched again, so once this has run the project works offline.

    python fetch_data.py             # fetch whatever is missing
    python fetch_data.py --refresh   # re-download everything (JPL updates daily)
"""

from __future__ import annotations

import argparse
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

from config import (
    APOPHIS_CAD_DIST_MAX_AU, APOPHIS_DES, APOPHIS_ELEMENTS_WINDOW, APOPHIS_HORIZONS_ID,
    APOPHIS_VECTOR_WINDOWS, CAD_CHUNK_EDGES, CAD_PARAMS, CAD_URL, DATA_DIR, FIREBALL_URL,
    HORIZONS_URL, MANIFEST_PATH, MOON_VECTOR_WINDOW, SBDB_QUERY_FIELDS, SBDB_QUERY_URL,
    SBDB_URL, SENTRY_URL,
)

TIMEOUT_S = 300
#: Pause between requests. JPL asks users not to hammer its APIs.
PAUSE_S = 1.0
USER_AGENT = "astrophysics-atlas/E3 (hobby research project; python-requests)"

CAD_FIELDS = ["des", "orbit_id", "jd", "cd", "dist", "dist_min", "dist_max", "v_rel",
              "v_inf", "t_sigma_f", "h", "diameter", "diameter_sigma", "fullname"]
FIREBALL_FIELDS = ["date", "energy", "impact-e", "lat", "lat-dir", "lon", "lon-dir", "alt", "vel"]


# ---------------------------------------------------------------------------
# Plumbing: provenance, atomic writes, validated requests
# ---------------------------------------------------------------------------

def _now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def load_manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text()) if MANIFEST_PATH.exists() else {}


def _record(key: str, path: Path, url: str, params: dict, rows: int) -> dict:
    manifest = load_manifest()
    entry = {
        "file": str(path.relative_to(DATA_DIR)),
        "url": url,
        "params": params,
        "retrieved_utc": _now_utc(),
        "rows": int(rows),
        "bytes": path.stat().st_size,
    }
    manifest[key] = entry
    tmp = MANIFEST_PATH.with_name(MANIFEST_PATH.name + ".part")
    tmp.write_text(json.dumps(manifest, indent=2))
    tmp.replace(MANIFEST_PATH)
    return entry


def _write_atomic(path: Path, text: str) -> None:
    """Write via a temporary name, so an interrupted download never leaves a
    truncated file where the analysis expects a complete one."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    tmp.write_text(text)
    tmp.replace(path)


def _get(url: str, params: dict, *, retries: int = 3) -> requests.Response:
    last_error = None
    for attempt in range(retries):
        try:
            response = requests.get(url, params=params, timeout=TIMEOUT_S,
                                    headers={"User-Agent": USER_AGENT})
            if response.status_code in (429, 503):  # busy: back off and retry
                raise requests.HTTPError(f"HTTP {response.status_code}")
            response.raise_for_status()
            time.sleep(PAUSE_S)
            return response
        except (requests.ConnectionError, requests.Timeout, requests.HTTPError) as exc:
            last_error = exc
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"{url} failed after {retries} attempts: {last_error}")


def _fetch_json(key: str, path: Path, url: str, params: dict, validate, *,
                refresh: bool, quiet: bool) -> dict:
    """Fetch a JSON API response unless cached; ``validate`` returns a row count
    or raises. The raw response text is what gets cached."""
    if path.exists() and not refresh and key in load_manifest():
        return json.loads(path.read_text())
    if not quiet:
        shown = ", ".join(f"{k}={v}" for k, v in params.items())
        print(f"  {url.rsplit('/', 1)[-1]}  {shown}")
    response = _get(url, params)
    try:
        payload = response.json()
    except ValueError as exc:
        raise RuntimeError(f"{url} did not return JSON: {response.text[:200]!r}") from exc
    rows = validate(payload)
    _write_atomic(path, response.text)
    _record(key, path, response.url, params, rows)
    if not quiet:
        print(f"    -> {path.name}: {rows:,} rows, {path.stat().st_size / 1e6:.2f} MB")
    return payload


# ---------------------------------------------------------------------------
# Validators: each returns the number of rows, or raises
# ---------------------------------------------------------------------------

def _check_table(payload: dict, fields: list[str], *, min_rows: int, what: str) -> int:
    if "fields" not in payload and int(payload.get("count", 0)) == 0 and min_rows == 0:
        return 0
    got = payload.get("fields")
    if got != fields:
        raise RuntimeError(f"{what}: unexpected fields {got!r} (message: {payload.get('message')!r})")
    n = int(payload["count"])
    if n != len(payload["data"]):
        raise RuntimeError(f"{what}: count {n} but {len(payload['data'])} rows")
    if n < min_rows:
        raise RuntimeError(f"{what}: only {n} rows -- refusing to cache")
    return n


def _check_records(payload: dict, *, min_rows: int, what: str, required: set[str]) -> int:
    data = payload.get("data")
    if not isinstance(data, list) or len(data) < min_rows:
        raise RuntimeError(f"{what}: bad or short response ({str(payload)[:200]})")
    missing = required - set(data[0])
    if missing:
        raise RuntimeError(f"{what}: records lack {sorted(missing)}")
    return len(data)


# ---------------------------------------------------------------------------
# The individual datasets
# ---------------------------------------------------------------------------

def cad_chunk_path(start: str, stop: str) -> Path:
    return DATA_DIR / "cad" / f"cad_{start}_{stop}.json"


def fetch_close_approaches(*, refresh: bool = False, quiet: bool = False) -> None:
    """All NEA approaches within 0.05 au, 1900-2200, in time chunks."""
    for start, stop in zip(CAD_CHUNK_EDGES[:-1], CAD_CHUNK_EDGES[1:]):
        params = {**CAD_PARAMS, "date-min": start, "date-max": stop}
        _fetch_json(f"cad_{start}_{stop}", cad_chunk_path(start, stop), CAD_URL, params,
                    lambda p: _check_table(p, CAD_FIELDS, min_rows=100, what="close approaches"),
                    refresh=refresh, quiet=quiet)


def fetch_apophis_approaches(*, refresh: bool = False, quiet: bool = False) -> None:
    params = {"des": APOPHIS_DES, "date-min": CAD_CHUNK_EDGES[0], "date-max": CAD_CHUNK_EDGES[-1],
              "dist-max": str(APOPHIS_CAD_DIST_MAX_AU), "diameter": "true", "fullname": "true"}
    _fetch_json("cad_apophis", DATA_DIR / "cad_apophis.json", CAD_URL, params,
                lambda p: _check_table(p, CAD_FIELDS, min_rows=20, what="Apophis approaches"),
                refresh=refresh, quiet=quiet)


def fetch_nea_catalog(*, refresh: bool = False, quiet: bool = False) -> None:
    params = {"fields": ",".join(SBDB_QUERY_FIELDS), "sb-group": "neo", "sb-kind": "a"}
    _fetch_json("nea_catalog", DATA_DIR / "nea_catalog.json", SBDB_QUERY_URL, params,
                lambda p: _check_table(p, SBDB_QUERY_FIELDS, min_rows=30_000, what="NEA catalogue"),
                refresh=refresh, quiet=quiet)


def fetch_sentry(*, refresh: bool = False, quiet: bool = False) -> None:
    """The risk list three ways: per object, per possible impact, and removed."""
    _fetch_json("sentry_objects", DATA_DIR / "sentry_objects.json", SENTRY_URL, {},
                lambda p: _check_records(p, min_rows=500, what="Sentry objects",
                                         required={"des", "ip", "ps_cum", "ts_max", "diameter", "v_inf"}),
                refresh=refresh, quiet=quiet)
    _fetch_json("sentry_impacts", DATA_DIR / "sentry_impacts.json", SENTRY_URL, {"all": "1"},
                lambda p: _check_records(p, min_rows=5_000, what="Sentry virtual impactors",
                                         required={"des", "date", "ip", "energy", "ps", "ts"}),
                refresh=refresh, quiet=quiet)
    _fetch_json("sentry_removed", DATA_DIR / "sentry_removed.json", SENTRY_URL, {"removed": "1"},
                lambda p: _check_records(p, min_rows=500, what="Sentry removed objects",
                                         required={"des", "removed"}),
                refresh=refresh, quiet=quiet)


def fetch_fireballs(*, refresh: bool = False, quiet: bool = False) -> None:
    _fetch_json("fireballs", DATA_DIR / "fireballs.json", FIREBALL_URL, {},
                lambda p: _check_table(p, FIREBALL_FIELDS, min_rows=500, what="fireballs"),
                refresh=refresh, quiet=quiet)


def fetch_apophis_record(*, refresh: bool = False, quiet: bool = False) -> None:
    def validate(p):
        if "object" not in p or "orbit" not in p:
            raise RuntimeError(f"SBDB record malformed: {str(p)[:200]}")
        return 1
    _fetch_json("sbdb_apophis", DATA_DIR / "sbdb_apophis.json", SBDB_URL,
                {"sstr": APOPHIS_DES, "phys-par": "1", "discovery": "1"}, validate,
                refresh=refresh, quiet=quiet)


# --- Horizons ---------------------------------------------------------------

def _horizons_validate(payload: dict) -> int:
    text = payload.get("result", "")
    if "$$SOE" not in text or "$$EOE" not in text:
        raise RuntimeError(f"Horizons returned no ephemeris:\n{text[:600]}")
    body = text.split("$$SOE", 1)[1].split("$$EOE", 1)[0]
    return sum(1 for line in body.splitlines() if line.strip())


def _q(value: str) -> str:
    """Horizons wants most parameter values in single quotes."""
    return f"'{value}'"


def _horizons(key: str, params: dict, *, refresh: bool, quiet: bool) -> None:
    full = {"format": "json", "OBJ_DATA": _q("NO"), "MAKE_EPHEM": _q("YES"),
            "CSV_FORMAT": _q("YES"), **params}
    _fetch_json(key, DATA_DIR / "horizons" / f"{key}.json", HORIZONS_URL, full,
                _horizons_validate, refresh=refresh, quiet=quiet)


def _vector_params(command: str, start: str, stop: str, step: str) -> dict:
    """Geometric geocentric state vectors in Earth's equatorial (ICRF) frame, km
    and km/s. Horizons tags vector tables with TDB times."""
    return {"COMMAND": _q(command), "EPHEM_TYPE": _q("VECTORS"), "CENTER": _q("500@399"),
            "REF_PLANE": _q("FRAME"), "REF_SYSTEM": _q("ICRF"), "VEC_TABLE": _q("2"),
            "VEC_CORR": _q("NONE"), "OUT_UNITS": _q("KM-S"),
            "START_TIME": _q(start), "STOP_TIME": _q(stop), "STEP_SIZE": _q(step)}


def fetch_horizons(*, refresh: bool = False, quiet: bool = False) -> None:
    for i, (start, stop, step) in enumerate(APOPHIS_VECTOR_WINDOWS):
        _horizons(f"apophis_vectors_{i}", _vector_params(APOPHIS_HORIZONS_ID, start, stop, step),
                  refresh=refresh, quiet=quiet)
    start, stop, step = MOON_VECTOR_WINDOW
    _horizons("moon_vectors", _vector_params("301", start, stop, step), refresh=refresh, quiet=quiet)

    # Geocentric range and its formal 3-sigma uncertainty (quantities 20 and 39),
    # minute by minute through closest approach. Observer tables are in UTC.
    start, stop, _ = APOPHIS_VECTOR_WINDOWS[1]
    _horizons("apophis_range_sigma", {
        "COMMAND": _q(APOPHIS_HORIZONS_ID), "EPHEM_TYPE": _q("OBSERVER"),
        "CENTER": _q("500@399"), "QUANTITIES": _q("20,39"), "TIME_TYPE": _q("UT"),
        "START_TIME": _q(start), "STOP_TIME": _q(stop), "STEP_SIZE": _q("1 m"),
    }, refresh=refresh, quiet=quiet)

    # Heliocentric osculating elements, ecliptic frame, daily across the flyby.
    start, stop, step = APOPHIS_ELEMENTS_WINDOW
    _horizons("apophis_elements", {
        "COMMAND": _q(APOPHIS_HORIZONS_ID), "EPHEM_TYPE": _q("ELEMENTS"),
        "CENTER": _q("500@10"), "REF_PLANE": _q("ECLIPTIC"), "REF_SYSTEM": _q("ICRF"),
        "OUT_UNITS": _q("AU-D"), "START_TIME": _q(start), "STOP_TIME": _q(stop),
        "STEP_SIZE": _q(step),
    }, refresh=refresh, quiet=quiet)


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------

def ensure_data(*, refresh: bool = False, quiet: bool = False) -> dict:
    """Fetch anything missing from the cache; a no-op when it is complete."""
    steps = [
        ("Close approaches, 1900-2200 (JPL CAD)", fetch_close_approaches),
        ("Apophis: every approach within 0.5 au", fetch_apophis_approaches),
        ("Near-Earth asteroid catalogue (JPL SBDB)", fetch_nea_catalog),
        ("Sentry impact monitoring (CNEOS)", fetch_sentry),
        ("Fireballs (CNEOS)", fetch_fireballs),
        ("Apophis database record (JPL SBDB)", fetch_apophis_record),
        ("Apophis and the Moon in April 2029 (JPL Horizons)", fetch_horizons),
    ]
    for title, step in steps:
        if not quiet:
            print(title)
        step(refresh=refresh, quiet=quiet)
    return load_manifest()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true",
                        help="re-download everything (JPL's tables change daily)")
    args = parser.parse_args(argv)
    manifest = ensure_data(refresh=args.refresh)
    total = sum(e["bytes"] for e in manifest.values())
    print(f"\nDone: {len(manifest)} cached files, {total / 1e6:.1f} MB in {DATA_DIR.name}/. "
          f"The analysis can now run offline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
