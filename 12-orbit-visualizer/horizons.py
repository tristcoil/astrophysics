"""JPL Horizons: what we ask for, how it is cached, and how the answers are read.

Every query is sent through ``astroquery.jplhorizons``, and the **raw text**
Horizons returns is what gets cached, one file per query, under
``data/horizons/``. Keeping the raw response rather than a parsed table costs a
little disk and buys two things. Nothing is lost in translation, because the
parser below is the only reader. And each file carries its own provenance: the
header says which planetary ephemeris, which spacecraft reconstruction, or which
comet orbit solution (with its number of observations and data arc) produced
the numbers.

Conventions for everything in this module and downstream:

* positions in au, velocities in au/day, times as Julian dates in TDB;
* frame: ecliptic and equinox of J2000 (Horizons' ``ECLIPTIC`` reference plane
  with the ICRF x-axis), the standard frame for solar-system orbits;
* *geometric* states -- where the body actually is at that instant, with no
  light-travel-time or aberration corrections -- which is what orbital dynamics
  needs. (What you would *see* is a different question, handled by the
  observer tables used for the night-sky chart.)
"""

from __future__ import annotations

import json
import re
import time
import warnings
from dataclasses import dataclass
from datetime import datetime, timezone

import numpy as np

from config import (
    CENTER_GEO, CENTER_SSB, CENTER_SUN, COMETS, DWARFS, FLYBY_HALF_WINDOW_DAYS,
    FLYBY_STEP, FLYBYS, HALE_BOPP_SPAN, HALLEY_SPAN, HORIZONS_DIR, INTERSTELLAR,
    INTERSTELLAR_PERIHELION, ISO_BARY_STEP, ISO_LONG_SPAN, ISO_LONG_STEP,
    ISO_PASSAGE_HALF_WINDOW_DAYS, ISO_PASSAGE_STEP, LIGHT_DAY_SPAN, LIGHT_DAY_STEP,
    MANIFEST_PATH, OBSERVER, ORBIT_SPAN_CERES, ORBIT_SPAN_INNER, ORBIT_SPAN_OUTER,
    PLANET_CENTER, PLANETS, SAMPLES_PER_ORBIT, SNAPSHOT_DATE, SPACECRAFT,
    SPACECRAFT_STEP, Body, body,
)

warnings.filterwarnings("ignore", module="astroquery")

#: Horizons observer-table quantities for the sky chart: 1 astrometric RA/Dec,
#: 4 apparent azimuth/elevation, 9 visual magnitude, 20 range, 29 constellation.
SKY_QUANTITIES = "1,4,9,20,29"

#: Bodies shown in the night-sky chart. Horizons major-body IDs (centres).
SKY_TARGETS = [("sun", "Sun", "10"), ("moon", "Moon", "301"), ("mercury", "Mercury", "199"),
               ("venus", "Venus", "299"), ("mars", "Mars", "499"), ("jupiter", "Jupiter", "599"),
               ("saturn", "Saturn", "699"), ("uranus", "Uranus", "799"),
               ("neptune", "Neptune", "899")]


@dataclass(frozen=True)
class Query:
    """One request to Horizons, and the file its answer is cached in."""

    key: str
    horizons_id: str
    kind: str  # vectors | elements | observer
    center: str | dict
    start: str | None = None
    stop: str | None = None
    step: str | None = None
    id_type: str | None = None
    closest_apparition: bool = False
    purpose: str = ""

    @property
    def path(self):
        return HORIZONS_DIR / f"{self.key}.txt"


# --------------------------------------------------------------------------
# The query plan
# --------------------------------------------------------------------------

def _step_for(period_days: float) -> str:
    """A sampling step giving ``SAMPLES_PER_ORBIT`` samples per orbit."""
    step_days = period_days / SAMPLES_PER_ORBIT
    if step_days >= 2:
        return f"{round(step_days)}d"
    return f"{round(step_days * 1440)}m"


def _shift(date: str, days: float) -> str:
    from datetime import date as _date, timedelta

    d = _date.fromisoformat(date) + timedelta(days=days)
    return d.isoformat()


def _next_day(date: str) -> str:
    return _shift(date, 1)


def orbit_key(b: Body) -> str:
    return f"orbit_{b.key}"


def track_key(b: Body) -> str:
    return f"track_{b.key}"


# Queries that depend on the snapshot date carry it in their key, so moving the
# date fetches fresh answers instead of silently reusing the old day's.

def snapshot_key(b: Body) -> str:
    return f"snapshot_{SNAPSHOT_DATE}_{b.key}"


def snapshot_geo_key(b: Body) -> str:
    return f"snapshotgeo_{SNAPSHOT_DATE}_{b.key}"


def elements_key(b: Body) -> str:
    return f"elements_{SNAPSHOT_DATE}_{b.key}"


def sky_key(name: str) -> str:
    return f"sky_{SNAPSHOT_DATE}_{name}"


def flyby_key(craft: str, planet: str) -> str:
    return f"flyby_{craft}_{planet}"


def _vec(key, b: Body, center, start, stop, step, purpose) -> Query:
    return Query(key, b.horizons_id, "vectors", center, start, stop, step,
                 id_type=b.id_type, closest_apparition=b.closest_apparition, purpose=purpose)


def query_plan() -> list[Query]:
    """Every Horizons request the project makes, in a fixed order."""
    plan: list[Query] = []
    snap_stop = _next_day(SNAPSHOT_DATE)

    # 1. Orbits long enough to measure periods and sizes: Kepler's third law.
    for b in PLANETS + DWARFS:
        if b.key == "ceres":
            span = ORBIT_SPAN_CERES
        elif b.key in ("mercury", "venus", "earth", "mars"):
            span = ORBIT_SPAN_INNER
        else:
            span = ORBIT_SPAN_OUTER
        plan.append(_vec(orbit_key(b), b, CENTER_SUN, *span, _step_for(b.period_days),
                         "Kepler's third law: period and size of the orbit"))

    # 2. Osculating elements at the snapshot: the GM Horizons itself uses for
    #    each orbit (Sun plus the planet), and a cross-check of our own
    #    state-vector-to-elements conversion.
    for b in PLANETS + DWARFS:
        plan.append(Query(elements_key(b), b.horizons_id, "elements", CENTER_SUN,
                          SNAPSHOT_DATE, snap_stop, "1d", id_type=b.id_type,
                          purpose="osculating elements and Keplerian GM at the snapshot"))

    # 3. Where everything is today.
    for b in PLANETS + DWARFS + COMETS + INTERSTELLAR + SPACECRAFT:
        plan.append(_vec(snapshot_key(b), b, CENTER_SUN, SNAPSHOT_DATE, snap_stop, "1d",
                         "heliocentric position today"))
    for b in SPACECRAFT:
        plan.append(_vec(snapshot_geo_key(b), b, CENTER_GEO, SNAPSHOT_DATE, snap_stop,
                         "1d", "distance and direction from Earth today"))

    # 4. Comets.
    halley = body("halley")
    plan.append(_vec(track_key(halley), halley, CENTER_SUN, *HALLEY_SPAN, "1d",
                     "Kepler's second law: one full orbit, perihelion 1986 to 2061"))
    hb = body("hale_bopp")
    plan.append(_vec(track_key(hb), hb, CENTER_SUN, *HALE_BOPP_SPAN, "2d",
                     "a long-period comet's plunge"))

    # 5. Interstellar objects.
    for b in INTERSTELLAR:
        tp = INTERSTELLAR_PERIHELION[b.key]
        plan.append(_vec(track_key(b), b, CENTER_SUN, *ISO_LONG_SPAN, ISO_LONG_STEP,
                         "trajectory through the solar system"))
        plan.append(_vec(f"passage_{b.key}", b, CENTER_SUN,
                         _shift(tp, -ISO_PASSAGE_HALF_WINDOW_DAYS),
                         _shift(tp, ISO_PASSAGE_HALF_WINDOW_DAYS), ISO_PASSAGE_STEP,
                         "perihelion passage at high time resolution"))
        plan.append(_vec(f"bary_{b.key}", b, CENTER_SSB, *ISO_LONG_SPAN, ISO_BARY_STEP,
                         "barycentric track: arrival and departure speed at infinity"))

    # 6. Spacecraft: the whole trajectory, then each planetary flyby up close.
    for b in SPACECRAFT:
        plan.append(_vec(track_key(b), b, CENTER_SUN, b.start, b.stop, SPACECRAFT_STEP,
                         "trajectory from launch to the end of the prediction"))
    for craft, planet, date in FLYBYS:
        b = body(craft)
        plan.append(_vec(flyby_key(craft, planet), b, PLANET_CENTER[planet],
                         _shift(date, -FLYBY_HALF_WINDOW_DAYS),
                         _shift(date, FLYBY_HALF_WINDOW_DAYS), FLYBY_STEP,
                         f"gravity assist at {planet.capitalize()}, planet-centred"))

    # 7. Voyager 1 and the light-day: its distance from Earth, and both bodies
    #    relative to the barycentre so light-travel time can be solved exactly.
    v1 = body("voyager1")
    plan.append(_vec("lightday_v1_geo", v1, CENTER_GEO, *LIGHT_DAY_SPAN, LIGHT_DAY_STEP,
                     "Voyager 1 distance from Earth"))
    plan.append(_vec("lightday_v1_ssb", v1, CENTER_SSB, *LIGHT_DAY_SPAN, LIGHT_DAY_STEP,
                     "Voyager 1 relative to the solar-system barycentre"))
    plan.append(Query("lightday_earth_ssb", "399", "vectors", CENTER_SSB, *LIGHT_DAY_SPAN,
                      LIGHT_DAY_STEP, purpose="Earth relative to the solar-system barycentre"))

    # 8. The sky tonight, from the observing site: evening of the snapshot date
    #    to the following morning, in UT.
    site = {"lon": OBSERVER["lon"], "lat": OBSERVER["lat"],
            "elevation": OBSERVER["elevation_km"]}
    for key, _name, hid in SKY_TARGETS:
        plan.append(Query(sky_key(key), hid, "observer", site,
                          f"{SNAPSHOT_DATE} 12:00", f"{snap_stop} 12:00", "10m",
                          purpose=f"{_name} as seen from {OBSERVER['name']}"))
    return plan


# --------------------------------------------------------------------------
# Fetching
# --------------------------------------------------------------------------

class HorizonsError(RuntimeError):
    pass


def _send(q: Query):
    from astroquery.jplhorizons import Horizons

    epochs = {"start": q.start, "stop": q.stop, "step": q.step}
    obj = Horizons(id=q.horizons_id, location=q.center, epochs=epochs, id_type=q.id_type)
    comet_kw = {"closest_apparition": q.closest_apparition} if q.id_type else {}
    if q.kind == "vectors":
        return obj.vectors_async(refplane="ecliptic", aberrations="geometric", cache=False,
                                 **comet_kw)
    if q.kind == "elements":
        return obj.elements_async(refplane="ecliptic", cache=False, **comet_kw)
    if q.kind == "observer":
        return obj.ephemerides_async(quantities=SKY_QUANTITIES, refraction=True,
                                     cache=False, **comet_kw)
    raise ValueError(f"unknown query kind {q.kind!r}")


def fetch(q: Query, *, refresh: bool = False, retries: int = 3) -> tuple[str, bool]:
    """Return the raw Horizons text for a query, downloading it only if needed.

    Returns ``(text, downloaded)``. A response without a data block is an error
    message from Horizons (an unknown target, a date outside a trajectory's
    coverage), and is raised rather than cached.
    """
    if q.path.exists() and not refresh:
        return q.path.read_text(), False

    last_exc = None
    for attempt in range(retries):
        try:
            resp = _send(q)
            text = resp.text
            break
        except Exception as exc:  # network hiccups: back off and retry
            last_exc = exc
            time.sleep(2.0 * (attempt + 1))
    else:
        raise HorizonsError(f"{q.key}: request failed after {retries} tries: {last_exc}")

    if "$$SOE" not in text:
        msg = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("API ")]
        raise HorizonsError(f"{q.key}: Horizons returned no data:\n  " + "\n  ".join(msg[:12]))

    q.path.write_text(text)
    _record(q, resp.url, text)
    return text, True


def _record(q: Query, url: str, text: str) -> None:
    manifest = json.loads(MANIFEST_PATH.read_text()) if MANIFEST_PATH.exists() else {}
    meta = header_metadata(text)
    manifest[q.key] = {
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "purpose": q.purpose,
        "url": url,
        "target": meta.get("target"),
        "target_source": meta.get("target_source"),
        "center": meta.get("center"),
        "rows": text.count("\n", text.index("$$SOE"), text.index("$$EOE")) - 1,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=1, sort_keys=True))


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------

def _data_block(text: str) -> tuple[list[str], list[list[str]]]:
    """Column names and rows of the CSV block between $$SOE and $$EOE."""
    head, _, rest = text.partition("$$SOE")
    body_txt, _, _ = rest.partition("$$EOE")
    # The column-name line is the last non-rule line before $$SOE.
    header_lines = [ln for ln in head.splitlines() if ln.strip() and not ln.startswith("***")]
    names = [c.strip() for c in header_lines[-1].split(",")]
    rows = []
    for ln in body_txt.strip().splitlines():
        cells = [c.strip() for c in ln.split(",")]
        rows.append(cells)
    return names, rows


def _floats(rows, idx) -> np.ndarray:
    out = np.empty(len(rows))
    for i, r in enumerate(rows):
        try:
            out[i] = float(r[idx])
        except (ValueError, IndexError):
            out[i] = np.nan
    return out


def parse_vectors(text: str) -> dict:
    """State vectors as arrays: jd, r (N,3) au, v (N,3) au/day, lt, rg, rr."""
    names, rows = _data_block(text)
    col = {n: i for i, n in enumerate(names)}
    out = {"jd": _floats(rows, col["JDTDB"]),
           "date": [r[col["Calendar Date (TDB)"]] for r in rows]}
    out["r"] = np.column_stack([_floats(rows, col[k]) for k in ("X", "Y", "Z")])
    out["v"] = np.column_stack([_floats(rows, col[k]) for k in ("VX", "VY", "VZ")])
    for k, name in (("lt", "LT"), ("rg", "RG"), ("rr", "RR")):
        if name in col:
            out[k] = _floats(rows, col[name])
    out["meta"] = header_metadata(text)
    return out


def parse_elements(text: str) -> dict:
    """Osculating elements: e, q, i, Omega, w, Tp, n, M, nu, a, Q, P (au, deg, days)."""
    names, rows = _data_block(text)
    col = {n: i for i, n in enumerate(names)}
    keys = {"e": "EC", "q": "QR", "inc": "IN", "raan": "OM", "argp": "W", "tp": "Tp",
            "n": "N", "M": "MA", "nu": "TA", "a": "A", "Q": "AD", "P": "PR"}
    out = {"jd": _floats(rows, col["JDTDB"])}
    for k, name in keys.items():
        out[k] = _floats(rows, col[name])
    out["meta"] = header_metadata(text)
    return out


def parse_observer(text: str) -> dict:
    """Observer table: time, apparent azimuth/elevation, magnitude, constellation."""
    names, rows = _data_block(text)
    col = {n: i for i, n in enumerate(names)}

    def find(prefix):
        for n, i in col.items():
            if n.startswith(prefix):
                return i
        raise KeyError(f"no column starting {prefix!r} in {list(col)}")

    out = {"jd_utc": _floats(rows, find("Date_________JDUT")),
           "az": _floats(rows, find("Azi")), "el": _floats(rows, find("Elev")),
           "ra": _floats(rows, find("R.A.")), "dec": _floats(rows, find("DEC")),
           "mag": _floats(rows, find("APmag")),
           "delta_au": _floats(rows, find("delta")),
           "constellation": [r[find("Cnst")] for r in rows]}
    out["meta"] = header_metadata(text)
    return out


_META_PATTERNS = {
    "target": r"^Target body name:\s*(.*?)\s*\{source:\s*(.*?)\}",
    "center": r"^Center body name:\s*(.*?)\s*\{source",
    "keplerian_gm": r"^Keplerian GM\s*:\s*([-+0-9.Ee]+)",
    "revised": r"Revised:\s*(.*?)\s{2,}",
    "soln": r"Soln\.date:\s*(\S+)\s+# obs:\s*(\d+)\s*\(([^)]*)\)",
    "soln_ref": r"soln ref\.=\s*([^,]+),\s*data arc:\s*(\S+ to \S+)",
}


def header_metadata(text: str) -> dict:
    """Provenance pulled from a Horizons response header."""
    meta: dict = {}
    head = text.split("$$SOE")[0]
    for key, pat in _META_PATTERNS.items():
        m = re.search(pat, head, flags=re.MULTILINE)
        if not m:
            continue
        if key == "target":
            meta["target"], meta["target_source"] = m.group(1), m.group(2)
        elif key == "keplerian_gm":
            meta["keplerian_gm_au3_d2"] = float(m.group(1))
        elif key == "soln":
            meta["solution_date"], meta["n_obs"], meta["obs_span"] = (
                m.group(1), int(m.group(2)), m.group(3))
        elif key == "soln_ref":
            meta["solution"], meta["data_arc"] = m.group(1).strip(), m.group(2)
        else:
            meta[key] = m.group(1).strip()
    return meta


# --------------------------------------------------------------------------
# Loading from the cache
# --------------------------------------------------------------------------

_PARSERS = {"vectors": parse_vectors, "elements": parse_elements, "observer": parse_observer}
_PLAN_CACHE: dict[str, Query] | None = None


def plan_by_key() -> dict[str, Query]:
    global _PLAN_CACHE
    if _PLAN_CACHE is None:
        _PLAN_CACHE = {q.key: q for q in query_plan()}
    return _PLAN_CACHE


def load(key: str) -> dict:
    """Parse a cached Horizons response. Never touches the network."""
    q = plan_by_key()[key]
    if not q.path.exists():
        raise FileNotFoundError(f"{q.path} is not cached -- run `python fetch_data.py` first")
    return _PARSERS[q.kind](q.path.read_text())
