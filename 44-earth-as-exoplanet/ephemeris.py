"""JPL Horizons vectors for DSCOVR, the Moon and the Sun, and a transit finder.

Horizons (ssd.jpl.nasa.gov, no key) carries DSCOVR as body -78. We fetch
geocentric ICRF vectors (km, geometric, UT time tags) on a coarse grid for the
whole mission and on a fine grid around events, and cache each request as a
CSV under ``data/raw/horizons/``.

The finder predicts, without looking at any image, every time the Moon enters
EPIC's field of view and every time it crosses Earth's disc as seen from
DSCOVR: in front (a transit) or behind (an occultation).
"""

from __future__ import annotations

import io
import re

import numpy as np
import pandas as pd

import geometry as G
from config import EPIC_FOV_DEG, RAW_DIR

HORIZONS_API = "https://ssd.jpl.nasa.gov/api/horizons.api"
HORIZONS_DIR = RAW_DIR / "horizons"
HORIZONS_DIR.mkdir(parents=True, exist_ok=True)
BODIES = {"dscovr": "-78", "moon": "301", "sun": "10"}
#: Mission span covered by the coarse grid.
MISSION = ("2015-06-13", "2026-10-05")
COARSE_STEP = "2 h"
FINE_STEP = "2 m"
#: Half-width of EPIC's square field of view.
HALF_FOV_RAD = np.radians(EPIC_FOV_DEG / 2)


def _cache_path(body: str, start: str, stop: str, step: str):
    tag = re.sub(r"[^0-9A-Za-z]+", "", f"{start}_{stop}_{step}")
    return HORIZONS_DIR / f"{body}_{tag}.csv"


def fetch_vectors(body: str, start: str, stop: str, step: str) -> pd.DataFrame:
    """Geocentric ICRF position vectors (km) of one body, cached."""
    path = _cache_path(body, start, stop, step)
    if not path.exists():
        from fetch_data import _record, _write_atomic, http_get

        params = {
            "format": "text", "COMMAND": f"'{BODIES[body]}'", "OBJ_DATA": "NO",
            "MAKE_EPHEM": "YES", "EPHEM_TYPE": "VECTORS", "CENTER": "'500@399'",
            "START_TIME": f"'{start}'", "STOP_TIME": f"'{stop}'", "STEP_SIZE": f"'{step}'",
            "REF_PLANE": "FRAME", "REF_SYSTEM": "ICRF", "VEC_TABLE": "1",
            "VEC_CORR": "NONE", "OUT_UNITS": "KM-S", "CSV_FORMAT": "YES", "TIME_TYPE": "UT",
        }
        text = http_get(HORIZONS_API, params=params).decode()
        if "$$SOE" not in text:
            raise RuntimeError(f"Horizons gave no ephemeris for {body}: {text[:500]}")
        block = text.split("$$SOE")[1].split("$$EOE")[0].strip()
        _write_atomic(path, ("jd,calendar,x,y,z\n" + "\n".join(
            ",".join(p.strip() for p in line.split(",")[:5]) for line in block.splitlines())).encode())
        _record(path, f"{HORIZONS_API} ({body} {start}..{stop} step {step})",
                "JPL Horizons geocentric ICRF vectors")
    df = pd.read_csv(path)
    df["t"] = pd.to_datetime(df["jd"], unit="D", origin="julian").dt.round("s")
    return df[["t", "x", "y", "z"]]


def vectors(start: str, stop: str, step: str) -> pd.DataFrame:
    """DSCOVR, Moon and Sun vectors on one time grid."""
    out = None
    for body in BODIES:
        df = fetch_vectors(body, start, stop, step).rename(
            columns={"x": f"{body}_x", "y": f"{body}_y", "z": f"{body}_z"})
        out = df if out is None else out.merge(df, on="t")
    return out


def _arr(df: pd.DataFrame, body: str) -> np.ndarray:
    return df[[f"{body}_x", f"{body}_y", f"{body}_z"]].to_numpy(float)


def moon_view(df: pd.DataFrame) -> pd.DataFrame:
    """Where the Moon appears relative to Earth from DSCOVR, at each time."""
    mv = G.moon_in_view(_arr(df, "dscovr"), _arr(df, "moon"), _arr(df, "sun"))
    out = pd.DataFrame({"t": df["t"]})
    for k in ("x_rad", "y_rad", "sep_rad", "r_moon_rad", "r_earth_rad", "moon_phase_deg"):
        out[k] = mv[k]
    out["in_front"] = mv["in_front"]
    out["phase_deg"] = G.phase_angle_deg(_arr(df, "dscovr"), _arr(df, "sun"))
    # Partly inside the square field of view.
    out["in_fov"] = ((np.abs(out["x_rad"]) < HALF_FOV_RAD + out["r_moon_rad"])
                     & (np.abs(out["y_rad"]) < HALF_FOV_RAD + out["r_moon_rad"]))
    out["overlaps_earth"] = out["sep_rad"] < out["r_earth_rad"] + out["r_moon_rad"]
    return out


def find_events(view: pd.DataFrame, gap: pd.Timedelta = pd.Timedelta("6h")) -> pd.DataFrame:
    """Group consecutive in-field-of-view times into events."""
    v = view[view["in_fov"]].copy()
    if v.empty:
        return pd.DataFrame()
    v["event"] = (v["t"].diff() > gap).cumsum()
    rows = []
    for _, g in v.groupby("event"):
        ov = g[g["overlaps_earth"]]
        rows.append({
            "fov_start": g["t"].iloc[0], "fov_end": g["t"].iloc[-1],
            "kind": ("transit" if g["in_front"].iloc[0] else "occultation") if len(ov) else "near miss",
            "overlap_start": ov["t"].iloc[0] if len(ov) else pd.NaT,
            "overlap_end": ov["t"].iloc[-1] if len(ov) else pd.NaT,
            "min_sep_earth_radii": float((g["sep_rad"] / g["r_earth_rad"]).min()),
            "t_min_sep": g["t"].iloc[int(np.argmin(g["sep_rad"].to_numpy()))],
            "moon_phase_deg": float(g["moon_phase_deg"].median()),
            "moon_diam_over_earth": float((g["r_moon_rad"] / g["r_earth_rad"]).median()),
        })
    return pd.DataFrame(rows)


def mission_events() -> pd.DataFrame:
    """Every Moon-in-view event of the mission, refined on a fine grid."""
    coarse = moon_view(vectors(*MISSION, COARSE_STEP))
    events = find_events(coarse)
    refined = []
    for _, e in events.iterrows():
        t0 = (e["fov_start"] - pd.Timedelta("3h")).strftime("%Y-%m-%d %H:%M")
        t1 = (e["fov_end"] + pd.Timedelta("3h")).strftime("%Y-%m-%d %H:%M")
        fine = find_events(moon_view(vectors(t0, t1, FINE_STEP)))
        refined.append(fine if len(fine) else e.to_frame().T)
    return pd.concat(refined, ignore_index=True) if refined else events
