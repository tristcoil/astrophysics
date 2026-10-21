"""Reader for the Sixth Catalog of Orbits of Visual Binary Stars (ORB6).

The fixed-width layout is documented in ``orb6format.txt`` (fetched with the
catalogue). Units are normalised on reading: periods to Julian years, semi-major
axes to arcsec, periastron epochs to decimal (Julian) years.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

import config as C
import orbits

# Identification columns are fixed (1-indexed, inclusive, from orb6format.txt).
ID_FIELDS = (
    ("ra_h", 1, 2), ("ra_m", 3, 4), ("ra_s", 5, 9),
    ("de_sign", 10, 10), ("de_d", 11, 12), ("de_m", 13, 14), ("de_s", 15, 18),
    ("wds", 20, 29), ("disc", 31, 44), ("ads", 46, 50), ("hd", 52, 57), ("hip", 59, 64),
    ("grade", 234, 234), ("note", 236, 236), ("ref", 238, 245), ("png", 247, 264),
)
#: Numeric fields: (name, first column, last column, unit/flag field or None).
#: The documented columns are not always respected: long periods in days,
#: large errors and a few other values spill one or two columns into the
#: neighbouring field, left or right (34 periods, 54 period errors, 32 epoch
#: errors in the October 2026 file). Each field therefore takes the
#: whitespace-delimited token that overlaps its documented columns most, and a
#: letter or symbol stuck to the number is its unit or flag.
NUM_FIELDS = (
    ("v1", 67, 71, "v1_flag"), ("v2", 74, 78, "v2_flag"),
    ("P_raw", 82, 92, "P_unit"), ("P_err_raw", 95, 104, None),
    ("a_raw", 106, 114, "a_unit"), ("a_err_raw", 117, 124, None),
    ("i", 126, 133, None), ("i_err", 135, 142, None),
    ("Omega", 144, 151, "Omega_flag"), ("Omega_err", 154, 161, None),
    ("T_raw", 163, 174, "T_unit"), ("T_err_raw", 177, 186, None),
    ("e", 188, 195, None), ("e_err", 197, 204, None),
    ("omega", 206, 213, "omega_flag"), ("omega_err", 215, 222, None),
    ("equinox", 224, 227, None), ("last_obs", 229, 232, None),
)
#: Where a separate unit/flag character sits when it is not glued to the number.
FLAG_COLS = {"v1_flag": 72, "v2_flag": 79, "P_unit": 93, "a_unit": 115, "Omega_flag": 152,
             "T_unit": 175, "omega_flag": 214}
INT_FIELDS = ("ra_h", "ra_m", "de_d", "de_m", "ads", "hd", "hip", "grade")
FLOAT_ID = ("ra_s", "de_s")

#: The first data line follows a fixed 7-line header (title, blank, 3 column
#: rulers, column legend, blank).
HEADER_LINES = 7
_NUM = re.compile(r"^([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)(.*)$")


def _num(s: str) -> float:
    s = s.strip()
    if s in ("", "."):
        return np.nan
    try:
        return float(s)
    except ValueError:
        return np.nan


def parse_line(line: str) -> dict:
    line = line.ljust(264)
    rec = {name: line[a - 1:b].strip() for name, a, b in ID_FIELDS}
    tokens = [(m.start() + 1, m.end(), m.group()) for m in re.finditer(r"\S+", line[64:232])]
    tokens = [(s0 + 64, e0 + 64, t) for s0, e0, t in tokens]
    used = set()
    for name, a, b, flag in NUM_FIELDS:
        best, ov = None, 0
        for k, (s0, e0, t) in enumerate(tokens):
            o = min(b, e0) - max(a, s0) + 1
            if o > ov and k not in used:
                best, ov = k, o
        value, suffix = np.nan, ""
        if best is not None:
            used.add(best)
            tok = tokens[best][2]
            m = _NUM.match(tok)
            if m:
                value, suffix = float(m.group(1)), m.group(2)
            elif tok == ".":
                value = np.nan
            else:
                suffix = tok if flag else ""
        rec[name] = value
        if flag:
            col = FLAG_COLS[flag]
            ch = line[col - 1].strip()
            rec[flag] = suffix or ch
    for k in INT_FIELDS:
        rec[k] = _num(rec[k])
    for k in FLOAT_ID:
        rec[k] = _num(rec[k])
    return rec


def read_orb6(path=C.ORB6_DIR / "orb6orbits.txt") -> pd.DataFrame:
    lines = path.read_text(encoding="latin-1").splitlines()[HEADER_LINES:]
    df = pd.DataFrame([parse_line(ln) for ln in lines if ln.strip()])

    sign = np.where(df["de_sign"] == "-", -1.0, 1.0)
    df["ra_deg"] = 15.0 * (df["ra_h"] + df["ra_m"] / 60.0 + df["ra_s"] / 3600.0)
    df["dec_deg"] = sign * (df["de_d"] + df["de_m"] / 60.0 + df["de_s"] / 3600.0)

    # Period -> Julian years. Codes: d days, y years, c centuries, h hours, m minutes.
    pscale = df["P_unit"].map({"d": 1.0 / C.JULIAN_YEAR, "y": 1.0, "c": 100.0, "h": 1.0 / 24 / C.JULIAN_YEAR,
                               "m": 1.0 / 1440 / C.JULIAN_YEAR})
    df["P_yr"] = df["P_raw"] * pscale
    df["P_err_yr"] = df["P_err_raw"] * pscale
    # Semi-major axis -> arcsec. Codes: a arcsec, m mas, M arcmin, u uas.
    ascale = df["a_unit"].map({"a": 1.0, "m": 1e-3, "M": 60.0, "u": 1e-6})
    df["a_arcsec"] = df["a_raw"] * ascale
    df["a_err_arcsec"] = df["a_err_raw"] * ascale
    # Periastron epoch -> Julian decimal year. Codes: y Besselian year, d JD-2400000,
    # m MJD, c centuries (Besselian year / 100); 'j' (undocumented, six recent
    # orbits with values like 2022.471) is taken as a Julian year.
    T = np.full(len(df), np.nan)
    Terr = np.full(len(df), np.nan)
    u = df["T_unit"].to_numpy()
    raw, rawerr = df["T_raw"].to_numpy(), df["T_err_raw"].to_numpy()
    for code, conv, econv in (
            ("y", lambda v: orbits.besselian_to_year(v), lambda e: e),
            ("c", lambda v: orbits.besselian_to_year(100.0 * v), lambda e: 100.0 * e),
            ("d", lambda v: orbits.jd_to_year(v + 2400000.0), lambda e: e / C.JULIAN_YEAR),
            ("m", lambda v: orbits.jd_to_year(v + 2400000.5), lambda e: e / C.JULIAN_YEAR),
            ("j", lambda v: v, lambda e: e)):
        m = u == code
        T[m] = conv(raw[m])
        Terr[m] = econv(rawerr[m])
    df["T_yr"] = T
    df["T_err_yr"] = Terr
    df["grade"] = df["grade"].fillna(0).astype(int)
    df["hip"] = df["hip"].fillna(0).astype(int)
    df["hd"] = df["hd"].fillna(0).astype(int)
    # The ascending node is identified (by radial velocities) when flagged '*'.
    df["node_known"] = df["Omega_flag"].str.contains(r"\*", regex=True)
    df["name"] = df["wds"] + " " + df["disc"]
    return df


def find(df: pd.DataFrame, wds: str, disc: str) -> pd.Series:
    """The ORB6 row for a WDS id and discoverer designation (spaces ignored)."""
    key = disc.replace(" ", "")
    m = (df["wds"] == wds) & (df["disc"].str.replace(" ", "", regex=False) == key)
    if m.sum() != 1:
        raise KeyError(f"{wds} {disc}: {m.sum()} matches in ORB6")
    return df[m].iloc[0]


def elements(row) -> tuple:
    """(P, T, e, a, i, Omega, omega) of an ORB6 row, ready for orbits.*."""
    return (row["P_yr"], row["T_yr"], row["e"], row["a_arcsec"], row["i"], row["Omega"], row["omega"])
