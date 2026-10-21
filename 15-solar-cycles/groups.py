"""Reader for the RGO/USAF/NOAA sunspot-group database (one line per group per day).

Fixed-width columns (``data/sunspot_groups/format.txt``), 1-indexed:

    1-4   year          5-6  month        7-8  day      9-12 time (fraction of day)
    13-20 group number (Greenwich to 1976; USAF/NOAA region number after)
    21-24 sub-group / group-type codes (layout changes in 1977 and 1982)
    26-29 umbral area to 1981; McIntosh class from 1982
    31-34 whole-spot area, millionths of the visible disc, as observed
    36-39 corrected umbral area to 1981; number of spots from 1982
    41-44 whole-spot area corrected for foreshortening, millionths of a hemisphere
    46-50 distance from disc centre (radii)     52-56 position angle (deg)
    58-62 Carrington longitude (deg)            64-68 latitude (deg, south negative)
    70-74 central-meridian distance (deg, east negative)

Three conventions matter and are handled here:

* A day with **no record at all** was observed and had no spots.
* A day whose record has group number 0 and every other field zero was **not
  observed** (Hathaway marks the 123 gaps in the RGO record this way; the first
  days of May 1874, before the series starts, look the same).
* From 1977 the USAF/NOAA positions are whole degrees and the areas run ~1.4x
  smaller than RGO's for the same spots -- see ``analysis.measure_area_seam``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import solar_physics as sp
from config import AR_DIR, AR_FILE_PATTERN, AR_FIRST_YEAR, RGO_LAST_YEAR

_PICKLE = AR_DIR.parent / "sunspot_groups.pkl"


def _f(field: str) -> float:
    field = field.strip()
    return float(field) if field else 0.0


def _parse_line(s: str):
    s = s.ljust(74)
    return (int(s[0:4]), int(s[4:6]), int(s[6:8]), _f(s[8:12]),
            s[12:20].strip(), s[20:24], s[25:29].strip(), int(_f(s[30:34])),
            s[35:39].strip(), int(_f(s[40:44])), _f(s[45:50]), _f(s[51:56]),
            _f(s[57:62]), _f(s[63:68]), _f(s[69:74]))


_COLUMNS = ["year", "month", "day", "frac", "group", "codes", "col26", "area_observed",
            "col36", "area", "distance", "position_angle", "carrington_lon", "lat", "cmd"]


def _files():
    files = sorted(AR_DIR.glob(AR_FILE_PATTERN.replace("{year}", "[0-9]" * 4)))
    years = [int(f.stem[1:]) for f in files]
    if not files or years[0] != AR_FIRST_YEAR or years != list(range(years[0], years[-1] + 1)):
        raise RuntimeError(f"expected yearly files g{AR_FIRST_YEAR}.txt onwards with no gaps "
                           f"in {AR_DIR}; run fetch_data.py")
    return files


def _parse_all() -> pd.DataFrame:
    rows = []
    for path in _files():
        with open(path, encoding="latin-1") as fh:  # universal newlines strip the CRs
            rows.extend(_parse_line(ln.rstrip("\n")) for ln in fh if ln.strip())
    df = pd.DataFrame(rows, columns=_COLUMNS)
    df["date"] = pd.to_datetime(df[["year", "month", "day"]])
    df["time"] = df["date"] + pd.to_timedelta(df["frac"], unit="D")
    df["t"] = sp.decimal_year(df["time"].to_numpy())
    df["jd"] = sp.julian_date(df["time"].to_numpy())
    df["era"] = np.where(df["year"] <= RGO_LAST_YEAR, "RGO",
                         np.where(df["year"] < 1982, "USAF", "NOAA"))
    # 1982 onwards: McIntosh class and spot count replace the umbral areas.
    modern = df["year"] >= 1982
    df["mcintosh"] = df["col26"].where(modern, "")
    df["n_spots"] = pd.to_numeric(df["col36"].where(modern), errors="coerce")
    # Sub-groups share a number: RGO marks them in the first code digit (group 89
    # of May 1874 is '00' and '10', ten degrees apart); USAF in the third code
    # character ('00A0', '00B0'); NOAA with a letter suffix in the first. The track
    # key keeps them apart.
    sub = np.select([df["era"] == "RGO", df["era"] == "USAF"],
                    [df["codes"].str[0], df["codes"].str[2]], df["codes"].str[0])
    df["track"] = df["group"] + ":" + pd.Series(sub, index=df.index).str.strip()
    # A group lasts one disc passage (<= 14 days); a handful of numbers are reused
    # months or years apart. Split a track wherever its sightings are > 30 days apart.
    df = df.sort_values(["track", "jd"], kind="stable")
    gap = df.groupby("track")["jd"].diff() > 30.0
    seg = gap.astype(int).groupby(df["track"]).cumsum()
    df["track"] = np.where(seg > 0, df["track"] + "#" + seg.astype(str), df["track"])
    df = df.sort_values(["jd"], kind="stable").reset_index(drop=True)
    df["missing_marker"] = df["group"].isin(["", "0"])
    return df.drop(columns=["col26", "col36"])


def load(use_cache: bool = True) -> pd.DataFrame:
    """Every record, parsed. Cached as a pickle beside the raw files."""
    newest = max(f.stat().st_mtime for f in _files())
    if use_cache and _PICKLE.exists() and _PICKLE.stat().st_mtime >= newest:
        return pd.read_pickle(_PICKLE)
    df = _parse_all()
    df.to_pickle(_PICKLE)
    return df


def split(records: pd.DataFrame) -> tuple[pd.DataFrame, pd.DatetimeIndex, pd.DatetimeIndex]:
    """(sunspot records, all days covered, days not observed).

    The record starts on the first day with a real group and ends on the last day
    of the last file's last month with data.
    """
    real = records[~records["missing_marker"]]
    first = real["date"].min()
    last = records["date"].max() + pd.offsets.MonthEnd(0)
    days = pd.date_range(first, last, freq="D")
    missing = pd.DatetimeIndex(records.loc[records["missing_marker"], "date"].unique())
    missing = missing[(missing >= first) & (missing <= last)].sort_values()
    return real.reset_index(drop=True), days, missing


def daily_area(spots: pd.DataFrame, days: pd.DatetimeIndex, missing: pd.DatetimeIndex,
               usaf_factor: float = 1.0) -> pd.DataFrame:
    """Total corrected spot area per day (millionths of a hemisphere), whole disc and
    each hemisphere. Unobserved days are NaN; USAF/NOAA areas are multiplied by
    ``usaf_factor``."""
    a = spots["area"].astype(float) * np.where(spots["era"] == "RGO", 1.0, usaf_factor)
    frame = pd.DataFrame({"date": spots["date"], "total": a,
                          "north": a.where(spots["lat"] > 0, 0.0),
                          "south": a.where(spots["lat"] < 0, 0.0)})
    out = frame.groupby("date")[["total", "north", "south"]].sum().reindex(days, fill_value=0.0)
    out.loc[out.index.isin(missing)] = np.nan
    out.index.name = "date"
    return out


def hathaway_daily_area() -> pd.DataFrame:
    """Hathaway's own derived daily areas (USAF/NOAA x1.4, -1 = no data): the answer
    key our ``daily_area`` is checked against."""
    df = pd.read_csv(AR_DIR / "daily_area.txt", sep=r"\s+", skiprows=1,
                     names=["year", "month", "day", "total", "north", "south"])
    df.index = pd.to_datetime(df[["year", "month", "day"]])
    return df[["total", "north", "south"]].where(df[["total", "north", "south"]] >= 0)


def listed_missing_days() -> pd.DatetimeIndex:
    """The days Hathaway lists as missing from the RGO record (MissingDays.txt)."""
    days = []
    for ln in (AR_DIR / "MissingDays.txt").read_text().splitlines():
        s = ln.rstrip()
        if len(s) >= 8 and s[:4].isdigit():
            try:
                days.append(pd.Timestamp(int(s[0:4]), int(s[4:6]), int(s[6:8])))
            except ValueError:
                continue
    return pd.DatetimeIndex(sorted(days))
