"""Readers for the WDC-SILSO files, the group-number reconstructions and the SWPC
forecasts. Each returns a tidy pandas DataFrame; missing values (SILSO's -1) become NaN.

File formats (semicolon-separated, documented on sidc.be/SILSO/infosn*):

* daily:    year; month; day; decimal year; SN; std; n_obs; provisional flag
* monthly:  year; month; decimal year; SN; std; n_obs; provisional flag
* yearly:   decimal year; SN; std; n_obs; provisional flag
* hemispheric monthly: year; month; decimal year; total; north; south;
            std total/north/south; n_obs total/north/south; provisional flag

The final column is 1 for definitive values and 0 for provisional ones (the last
few months, still open to revision).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from config import (
    DATA_DIR,
    GROUP_NUMBER_SERIES,
    SILSO_CYCLE_TABLE,
    SILSO_SERIES,
    SWPC_FILES,
)

SILSO_DIR = DATA_DIR / "silso"
SWPC_DIR = DATA_DIR / "swpc"


def _read(name: str, columns: list[str]) -> pd.DataFrame:
    df = pd.read_csv(SILSO_DIR / name, sep=";", header=None, names=columns,
                     skipinitialspace=True)
    for c in columns:
        if c.startswith(("sn", "std", "n_obs")):
            df[c] = df[c].astype(float).where(df[c] >= 0)
    return df


def _month_dates(df: pd.DataFrame) -> pd.Series:
    return pd.to_datetime(dict(year=df["year"], month=df["month"], day=15))


def daily() -> pd.DataFrame:
    df = _read(SILSO_SERIES["daily"],
               ["year", "month", "day", "t", "sn", "std", "n_obs", "definitive"])
    df["date"] = pd.to_datetime(df[["year", "month", "day"]])
    return df


def monthly() -> pd.DataFrame:
    df = _read(SILSO_SERIES["monthly"],
               ["year", "month", "t", "sn", "std", "n_obs", "definitive"])
    df["date"] = _month_dates(df)
    return df


def smoothed() -> pd.DataFrame:
    """SILSO's own 13-month smoothed series (our smoothing is checked against it)."""
    df = _read(SILSO_SERIES["smoothed"],
               ["year", "month", "t", "sn", "std", "n_obs", "definitive"])
    df["date"] = _month_dates(df)
    return df


def yearly() -> pd.DataFrame:
    df = _read(SILSO_SERIES["yearly"], ["t", "sn", "std", "n_obs", "definitive"])
    df["year"] = np.floor(df["t"]).astype(int)
    return df


def hemispheric(smoothed_series: bool = False) -> pd.DataFrame:
    name = SILSO_SERIES["hemispheric_smoothed" if smoothed_series else "hemispheric_monthly"]
    df = _read(name, ["year", "month", "t", "sn", "sn_north", "sn_south",
                      "std", "std_north", "std_south",
                      "n_obs", "n_obs_north", "n_obs_south", "definitive"])
    df["date"] = _month_dates(df)
    return df


def official_cycles() -> pd.DataFrame:
    """SILSO's table of cycle minima and maxima (13-month smoothed values)."""
    path = SILSO_DIR / Path(SILSO_CYCLE_TABLE).name
    raw = pd.read_csv(path, sep=";", header=None, skiprows=2, encoding="latin-1",
                      names=["cycle", "min_year", "min_month", "min_sn",
                             "max_year", "max_month", "max_sn"])
    raw = raw.dropna(subset=["cycle"])
    raw["cycle"] = raw["cycle"].astype(int)
    for c in ("min_year", "min_month", "max_year", "max_month"):
        raw[c] = pd.to_numeric(raw[c], errors="coerce").astype("Int64")
    for c in ("min_sn", "max_sn"):
        raw[c] = pd.to_numeric(raw[c].astype(str).str.replace(",", "."), errors="coerce")
    return raw.reset_index(drop=True)


def group_numbers() -> dict[str, pd.DataFrame]:
    """The four yearly group-number reconstructions, keyed as in config.

    Columns: t (decimal year, mid-year), gn, gn_err and n_obs (the number of
    observations or observing days behind the value; NaN where a file lacks them).
    """
    out = {}
    for key, spec in GROUP_NUMBER_SERIES.items():
        path = SILSO_DIR / Path(spec["file"]).name
        arr = np.genfromtxt(path)
        if arr.shape[1] < 4:
            arr = np.column_stack([arr[:, :3], np.full(len(arr), np.nan)])
        df = pd.DataFrame(arr[:, :4], columns=["t", "gn", "gn_err", "n_obs"])
        for c in ("gn", "gn_err", "n_obs"):
            df[c] = df[c].where(df[c] >= 0)
        if "own_until" in spec:
            df = df[df["t"] <= spec["own_until"]]
        out[key] = df.dropna(subset=["gn"]).reset_index(drop=True)
    return out


def _swpc(key: str) -> list[dict]:
    return json.loads((SWPC_DIR / SWPC_FILES[key]).read_text())


def panel_forecast() -> pd.DataFrame:
    """The 2019 NOAA/NASA panel's predicted smoothed sunspot number for Cycle 25,
    with its range: amplitude +/-10% (``low``/``high``) and timing +/-6 months."""
    base = pd.DataFrame(_swpc("panel_2019"))
    base["date"] = pd.to_datetime(base["time_tag"] + "-15")
    rng = pd.DataFrame(_swpc("panel_2019_range")).replace(-1.0, np.nan)
    rng["date"] = pd.to_datetime(rng["time-tag"] + "-15")
    df = base[["date", "smoothed_ssn"]].rename(columns={"smoothed_ssn": "ssn"})
    df = df.merge(rng, on="date", how="outer").sort_values("date")
    variants = [c for c in rng.columns if c.startswith("smoothed_ssn")]
    df["envelope_low"] = df[variants].min(axis=1)
    df["envelope_high"] = df[variants].max(axis=1)
    return df[["date", "ssn", "envelope_low", "envelope_high"]].reset_index(drop=True)


def current_forecast() -> pd.DataFrame:
    """SWPC's current forecast for the remainder of the cycle."""
    df = pd.DataFrame(_swpc("current_forecast"))
    df["date"] = pd.to_datetime(df["time-tag"] + "-15")
    return df[["date", "predicted_ssn", "low_ssn", "high_ssn"]]
