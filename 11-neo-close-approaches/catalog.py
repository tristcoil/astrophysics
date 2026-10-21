"""Load the cached JPL responses into clean pandas tables.

Every derived column is defined here, once:

* distances in km and lunar distances (LD), with JPL's 3-sigma minimum and
  maximum;
* a **best diameter** for every asteroid: the measured one where JPL has it
  (radar, infrared), otherwise estimated from brightness (H) with a typical
  albedo of 0.14 -- plus the range that albedos of 0.05-0.25 would allow;
* the size class;
* for each approach, how many days before or after it the asteroid was first
  observed (the **discovery lag**), from the catalogue's first-observation date.

Times from JPL's close-approach and Horizons vector tables are TDB (the
timescale of planetary motion), which runs 69 s ahead of UTC in the 2020s.
That is irrelevant at the resolution of any plot here, so TDB timestamps are
stored as ordinary datetimes; the Apophis flyby time is converted properly.
"""

from __future__ import annotations

import json
from datetime import datetime, timedelta
from functools import lru_cache

import numpy as np
import pandas as pd

import neo_physics as nphys
from config import (
    ALBEDO_RANGE, ALBEDO_TYPICAL, CAD_CHUNK_EDGES, DATA_DIR, SIZE_CLASS_EDGES_M,
    SIZE_CLASS_LABELS,
)
from fetch_data import cad_chunk_path, load_manifest

JD_UNIX_EPOCH = 2440587.5


def _read(name: str) -> dict:
    path = DATA_DIR / name
    if not path.exists():
        raise FileNotFoundError(f"{path} not found -- run `python fetch_data.py` first")
    return json.loads(path.read_text())


def _table(payload: dict) -> pd.DataFrame:
    return pd.DataFrame(payload.get("data", []), columns=payload.get("fields"))


def jd_to_datetime(jd) -> pd.Series:
    return pd.to_datetime(np.asarray(jd, float) - JD_UNIX_EPOCH, unit="D")


def parse_time_sigma_days(values) -> np.ndarray:
    """JPL's 3-sigma timing uncertainty, e.g. '< 00:01', '00:05', '3_19:55'
    (days_hours:minutes), as days. '< 00:01' (under a minute) becomes 0."""
    out = np.full(len(values), np.nan)
    for i, s in enumerate(values):
        if not isinstance(s, str):
            continue
        s = s.strip()
        if s.startswith("<"):
            out[i] = 0.0
            continue
        days, _, hm = s.rpartition("_")
        hours, minutes = hm.split(":")
        out[i] = (float(days) if days else 0.0) + float(hours) / 24 + float(minutes) / 1440
    return out


def size_class(d_m) -> pd.Categorical:
    """Label each diameter (m) with its class from ``config.SIZE_CLASS_EDGES_M``."""
    d_m = np.asarray(d_m, float)
    labels = np.full(d_m.shape, SIZE_CLASS_LABELS[-1], dtype=object)
    for edge, label in reversed(list(zip(SIZE_CLASS_EDGES_M, SIZE_CLASS_LABELS))):
        labels[d_m >= edge] = label
    labels[~np.isfinite(d_m)] = None
    return pd.Categorical(labels, categories=SIZE_CLASS_LABELS, ordered=True)


def add_sizes(df: pd.DataFrame, h_col: str = "h", measured_km_col: str = "diameter_km") -> pd.DataFrame:
    """Best diameter (m): measured if available, otherwise from H at p = 0.14."""
    h = pd.to_numeric(df[h_col], errors="coerce")
    df["d_est_m"] = nphys.diameter_km(h, ALBEDO_TYPICAL) * 1e3
    # albedo 0.25 (bright) gives the small end, 0.05 (dark) the large end
    df["d_est_lo_m"] = nphys.diameter_km(h, ALBEDO_RANGE[1]) * 1e3
    df["d_est_hi_m"] = nphys.diameter_km(h, ALBEDO_RANGE[0]) * 1e3
    measured = pd.to_numeric(df[measured_km_col], errors="coerce") * 1e3
    df["d_measured"] = measured.notna()
    df["d_best_m"] = measured.where(measured.notna(), df["d_est_m"])
    df["size_class"] = size_class(df["d_best_m"])
    return df


# ---------------------------------------------------------------------------
# Close approaches
# ---------------------------------------------------------------------------

def _clean_approaches(raw: pd.DataFrame) -> pd.DataFrame:
    df = raw.rename(columns={"dist": "dist_au", "dist_min": "dist_min_au", "dist_max": "dist_max_au",
                             "diameter": "diameter_km", "t_sigma_f": "t_sigma"})
    for col in ["jd", "dist_au", "dist_min_au", "dist_max_au", "v_rel", "v_inf", "h",
                "diameter_km", "diameter_sigma"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["fullname"] = df["fullname"].str.strip()
    df["time"] = jd_to_datetime(df["jd"])
    df["year"] = df["time"].dt.year
    df["dist_km"] = nphys.au_to_km(df["dist_au"])
    df["dist_ld"] = nphys.au_to_ld(df["dist_au"])
    df["dist_min_ld"] = nphys.au_to_ld(df["dist_min_au"])
    df["dist_max_ld"] = nphys.au_to_ld(df["dist_max_au"])
    df["t_sigma_days"] = parse_time_sigma_days(df["t_sigma"].tolist())
    return add_sizes(df)


@lru_cache(maxsize=1)
def load_approaches() -> pd.DataFrame:
    """Every cached NEA approach within 0.05 au, 1900-2200, one row each,
    joined to the catalogue's first-observation date."""
    parts = []
    for start, stop in zip(CAD_CHUNK_EDGES[:-1], CAD_CHUNK_EDGES[1:]):
        path = cad_chunk_path(start, stop)
        parts.append(_table(_read(str(path.relative_to(DATA_DIR)))))
    raw = pd.concat(parts, ignore_index=True)
    # Chunks meet at midnight on their shared date; never count an approach twice.
    raw = raw.drop_duplicates(subset=["des", "jd"]).reset_index(drop=True)
    df = _clean_approaches(raw)

    cat = load_nea_catalog()[["pdes", "first_obs", "last_obs", "pha", "class", "condition_code"]]
    df = df.merge(cat, how="left", left_on="des", right_on="pdes").drop(columns="pdes")
    # Whole days between first observation and closest approach: < 0 means
    # seen before, 0 on the day itself, > 0 only afterwards.
    approach_day = df["time"].dt.normalize()
    df["lag_days"] = (df["first_obs"] - approach_day).dt.days
    return df.sort_values("jd").reset_index(drop=True)


def load_apophis_approaches() -> pd.DataFrame:
    return _clean_approaches(_table(_read("cad_apophis.json")))


# ---------------------------------------------------------------------------
# The NEA catalogue
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def load_nea_catalog() -> pd.DataFrame:
    df = _table(_read("nea_catalog.json"))
    for col in ["H", "diameter", "albedo", "moid", "n_obs_used", "data_arc"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in ["first_obs", "last_obs"]:
        df[col] = pd.to_datetime(df[col], errors="coerce")
    df["full_name"] = df["full_name"].str.strip()
    df = df.rename(columns={"diameter": "diameter_km"})
    return add_sizes(df, h_col="H")


# ---------------------------------------------------------------------------
# Sentry
# ---------------------------------------------------------------------------

def decimal_year(value) -> float:
    """A date as a decimal year. Sentry's dates reach the 2880s, beyond the
    year 2262 at which pandas timestamps stop, so its dates are kept as
    floats. Accepts 'YYYY-MM-DD.dd' (fractional day), 'YYYY-MM-DD[ HH:MM:SS]'
    or a datetime."""
    if isinstance(value, str):
        value = value.strip()
        if len(value) > 10 and value[10] == ".":  # 'YYYY-MM-DD.dd'
            moment = datetime.fromisoformat(value[:10]) + timedelta(days=float("0" + value[10:]))
        else:
            moment = datetime.fromisoformat(value.replace(" ", "T"))
    else:
        moment = pd.Timestamp(value).to_pydatetime().replace(tzinfo=None)
    start = datetime(moment.year, 1, 1)
    length = (datetime(moment.year + 1, 1, 1) - start).total_seconds()
    return moment.year + (moment - start).total_seconds() / length


@lru_cache(maxsize=1)
def load_sentry_objects() -> pd.DataFrame:
    df = pd.DataFrame(_read("sentry_objects.json")["data"])
    for col in ["ip", "ps_cum", "ps_max", "ts_max", "diameter", "h", "v_inf", "n_imp"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["fullname"] = df["fullname"].str.strip()
    df["last_obs"] = pd.to_datetime(df["last_obs"].str.slice(0, 10), errors="coerce")
    years = df["range"].str.split("-", expand=True)
    df["first_year"] = pd.to_numeric(years[0], errors="coerce")
    df["last_year"] = pd.to_numeric(years[1], errors="coerce")
    df["diameter_m"] = df["diameter"] * 1e3
    return df


@lru_cache(maxsize=1)
def load_sentry_impacts() -> pd.DataFrame:
    """One row per possible impact date ("virtual impactor")."""
    df = pd.DataFrame(_read("sentry_impacts.json")["data"])
    for col in ["ip", "energy", "ps", "ts", "sigma_vi"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["fullname"] = df["fullname"].str.strip()
    df["impact_year"] = [decimal_year(d) for d in df["date"]]
    return df


@lru_cache(maxsize=1)
def load_sentry_removed() -> pd.DataFrame:
    df = pd.DataFrame(_read("sentry_removed.json")["data"])
    df["removed"] = pd.to_datetime(df["removed"], errors="coerce")
    return df.sort_values("removed").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Fireballs
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def load_fireballs() -> pd.DataFrame:
    df = _table(_read("fireballs.json"))
    df["time"] = pd.to_datetime(df["date"])
    df["radiated_j"] = pd.to_numeric(df["energy"], errors="coerce") * 1e10
    df["impact_kt"] = pd.to_numeric(df["impact-e"], errors="coerce")
    for col in ["lat", "lon", "alt", "vel"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["lat"] = np.where(df["lat-dir"] == "S", -df["lat"], df["lat"])
    df["lon"] = np.where(df["lon-dir"] == "W", -df["lon"], df["lon"])
    df["year"] = df["time"].dt.year
    return df.drop(columns=["energy", "impact-e", "lat-dir", "lon-dir"]).sort_values("time").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Apophis: database record and Horizons tables
# ---------------------------------------------------------------------------

def load_apophis_record() -> dict:
    return _read("sbdb_apophis.json")


def _horizons_rows(key: str) -> tuple[list[str], list[list[str]]]:
    text = _read(f"horizons/{key}.json")["result"]
    lines = text.splitlines()
    soe = next(i for i, line in enumerate(lines) if line.startswith("$$SOE"))
    eoe = next(i for i, line in enumerate(lines) if line.startswith("$$EOE"))
    header = [h.strip() for h in lines[soe - 2].split(",")]
    rows = [[c.strip() for c in line.split(",")] for line in lines[soe + 1:eoe] if line.strip()]
    return header, rows


def load_horizons_vectors(key: str) -> pd.DataFrame:
    """Geocentric state vectors: time (TDB), x, y, z (km), vx, vy, vz (km/s)."""
    header, rows = _horizons_rows(key)
    df = pd.DataFrame(rows, columns=header)[["JDTDB", "X", "Y", "Z", "VX", "VY", "VZ"]]
    df = df.astype(float).rename(columns=str.lower).rename(columns={"jdtdb": "jd_tdb"})
    df["time_tdb"] = jd_to_datetime(df["jd_tdb"])
    df["r_km"] = np.sqrt(df["x"] ** 2 + df["y"] ** 2 + df["z"] ** 2)
    return df


def load_apophis_range_sigma() -> pd.DataFrame:
    """Geocentric range and its formal 3-sigma uncertainty, minute by minute (UTC)."""
    header, rows = _horizons_rows("apophis_range_sigma")
    df = pd.DataFrame(rows, columns=header)
    out = pd.DataFrame({
        "time_utc": pd.to_datetime(df["Date__(UT)__HR:MN"], format="%Y-%b-%d %H:%M"),
        "range_km": nphys.au_to_km(df["delta"].astype(float)),
        "range_rate_kms": df["deldot"].astype(float),
        "range_3sigma_km": df["RNG_3sigma"].astype(float),
    })
    return out


def load_apophis_elements() -> pd.DataFrame:
    """Heliocentric osculating orbit (ecliptic J2000), daily: a, e, q, Q in au,
    i in degrees, period in days."""
    header, rows = _horizons_rows("apophis_elements")
    df = pd.DataFrame(rows, columns=header)
    out = pd.DataFrame({
        "jd_tdb": df["JDTDB"].astype(float),
        "a_au": df["A"].astype(float), "e": df["EC"].astype(float),
        "q_au": df["QR"].astype(float), "Q_au": df["AD"].astype(float),
        "i_deg": df["IN"].astype(float), "period_days": df["PR"].astype(float),
        "node_deg": df["OM"].astype(float), "peri_deg": df["W"].astype(float),
        "true_anomaly_deg": df["TA"].astype(float),
    })
    out["time_tdb"] = jd_to_datetime(out["jd_tdb"])
    return out


def snapshot_time() -> pd.Timestamp:
    """When the close-approach table was retrieved: the 'today' of this run."""
    manifest = load_manifest()
    stamps = [v["retrieved_utc"] for k, v in manifest.items() if k.startswith("cad_1") or k.startswith("cad_2")]
    return pd.Timestamp(max(stamps)).tz_convert(None)
