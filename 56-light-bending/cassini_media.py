"""The Earth's atmosphere on the Cassini radio links: troposphere and ionosphere.

Each leg of the two-way link crosses the troposphere at the station end. Its
delay is the zenith delay times a mapping function of the elevation; for a
spacecraft rising or setting, the delay changes by ~1 mm/s at 10 degrees, a
fractional frequency shift of ~3e-12, which the model must remove. Three
calibrations are available, from best to worst:

* **AMC**, the Advanced Media Calibration system at DSS 25: two water vapour
  radiometers (AWVR 1 and 2) pointed along the antenna's line of sight
  retrieve the zenith wet delay every 30 s; the zenith dry delay comes from
  surface pressure (PDS ``pd1``/``pd2`` files, format DORS-002). This is the
  calibration Bertotti, Iess & Tortora (2003) used.
* the **DSN's own calibrations** (``tro`` files, format TRK-2-23): zenith dry
  and wet delays per complex as normalized power series *added to seasonal
  models* (Fourier series). The seasonal model is printed only for Goldstone
  (TRK-2-23 Figure 3-2), so this works for DSS 15/25 only.
* a **standard atmosphere** (Saastamoinen dry delay from the station height,
  10 cm wet) for the Canberra and Madrid stations.

The mapping functions are Niell's (1996, JGR 101, 3227; NMF).

The Earth's ionosphere: DSN line-of-sight calibrations (``ion`` files,
TRK-2-23 CHPART: one-way delay in metres at 2295 MHz per pass and complex,
from GPS). A plasma advances the carrier phase, so the phase delay at
frequency f is -(2295 MHz / f)^2 times it.
"""

from __future__ import annotations

import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd

import config as C

CASSINI_DIR = C.RAW_DIR / "cassini_sce1"
#: DSN complex of each station
COMPLEX = {14: 10, 15: 10, 24: 10, 25: 10, 26: 10, 34: 40, 43: 40, 45: 40, 54: 60, 55: 60, 63: 60, 65: 60}
#: reference frequency of the ionosphere calibrations (TRK-2-23 3.2.1)
F_ION_REF = 2295e6
_EPOCH_1950 = np.datetime64("1950-01-01T00:00:00", "s")
_EPOCH_1970_IN_1950 = float((np.datetime64("1970-01-01") - np.datetime64("1950-01-01")).astype("timedelta64[s]").astype(int))


# ---------------------------------------------------------------------------
# Niell mapping functions
# ---------------------------------------------------------------------------

_NMF_LAT = np.array([15.0, 30.0, 45.0, 60.0, 75.0])
_NMF_HYD_AVG = np.array([
    [1.2769934e-3, 2.9153695e-3, 62.610505e-3],
    [1.2683230e-3, 2.9152299e-3, 62.837393e-3],
    [1.2465397e-3, 2.9288445e-3, 63.721774e-3],
    [1.2196049e-3, 2.9022565e-3, 63.824265e-3],
    [1.2045996e-3, 2.9024912e-3, 64.258455e-3]])
_NMF_HYD_AMP = np.array([
    [0.0, 0.0, 0.0],
    [1.2709626e-5, 2.1414979e-5, 9.0128400e-5],
    [2.6523662e-5, 3.0160779e-5, 4.3497037e-5],
    [3.4000452e-5, 7.2562722e-5, 84.795348e-5],
    [4.1202191e-5, 11.723375e-5, 170.37206e-5]])
_NMF_HT = (2.53e-5, 5.49e-3, 1.14e-3)
_NMF_WET = np.array([
    [5.8021897e-4, 1.4275268e-3, 4.3472961e-2],
    [5.6794847e-4, 1.5138625e-3, 4.6729510e-2],
    [5.8118019e-4, 1.4572752e-3, 4.3908931e-2],
    [5.9727542e-4, 1.5007428e-3, 4.4626982e-2],
    [6.1641693e-4, 1.7599082e-3, 5.4736038e-2]])


def _marini(sin_e, a, b, c):
    return (1 + a / (1 + b / (1 + c))) / (sin_e + a / (sin_e + b / (sin_e + c)))


def _interp_lat(table, lat_deg):
    lat = np.clip(abs(lat_deg), 15.0, 75.0)
    return np.array([np.interp(lat, _NMF_LAT, table[:, k]) for k in range(3)])


def niell(el_rad, lat_deg: float, height_m: float, doy) -> tuple[np.ndarray, np.ndarray]:
    """Niell (1996) hydrostatic and wet mapping functions at elevation el."""
    sin_e = np.sin(el_rad)
    phase = 2 * np.pi * (np.asarray(doy, dtype=float) - 28.0) / 365.25
    if lat_deg < 0:
        phase = phase + np.pi
    avg, amp = _interp_lat(_NMF_HYD_AVG, lat_deg), _interp_lat(_NMF_HYD_AMP, lat_deg)
    a, b, c = (avg[k] - amp[k] * np.cos(phase) for k in range(3))
    m_h = _marini(sin_e, a, b, c) + (1 / sin_e - _marini(sin_e, *_NMF_HT)) * height_m / 1e3
    m_w = _marini(sin_e, *_interp_lat(_NMF_WET, lat_deg))
    return m_h, m_w


# ---------------------------------------------------------------------------
# CSP media calibrations (TRK-2-23)
# ---------------------------------------------------------------------------

def _csp_time(s: str) -> float:
    """'02/06/01,06:00:00.001' (UTC) -> seconds since 1950 (86400 per day)."""
    d, t = s.split(",") if "," in s else (s, "00:00")
    yy, mm, dd = (int(x) for x in d.split("/"))
    year = 1900 + yy if yy >= 69 else 2000 + yy
    hms = [float(x) for x in t.split(":")] + [0.0, 0.0]
    day = (np.datetime64(f"{year:04d}-{mm:02d}-{dd:02d}") - np.datetime64("1950-01-01")).astype("timedelta64[D]").astype(int)
    return day * 86400.0 + hms[0] * 3600 + hms[1] * 60 + hms[2]


def read_csp(paths) -> pd.DataFrame:
    """Parse ADJUST commands: one row per calibration (kind = WET/DRY/CHPART)."""
    rows = []
    for p in paths:
        text = re.sub(r"#[^\n]*", "", Path(p).read_text())        # strip comments
        # numbers contain '.', which also ends a command: split on the verb instead
        for cmd in re.split(r"(?=ADJUST)", " ".join(text.split())):
            if not cmd.startswith("ADJUST"):
                continue
            m_kind = re.search(r"MODEL\s*\(\s*(WET NUPART|DRY NUPART|CHPART)\s*\)", cmd)
            m_from = re.search(r"FROM\s*\(([^)]*)\)\s*TO\s*\(([^)]*)\)", cmd)
            m_dsn = re.search(r"DSN\s*\(\s*C?(\d+)\s*\)", cmd)
            m_by = re.search(r"BY\s+(D?NRMPOW|D?CONST|D?TRIG)\s*\(([^)]*)\)", cmd)
            if not (m_kind and m_from and m_dsn and m_by):
                continue
            coef = [float(x.replace("D", "E")) for x in m_by.group(2).split(",")]
            rows.append({"kind": m_kind.group(1).split()[0], "complex": int(m_dsn.group(1)),
                         "t0": _csp_time(m_from.group(1)), "t1": _csp_time(m_from.group(2)),
                         "form": m_by.group(1).lstrip("D"), "coef": coef})
    return pd.DataFrame(rows)


def eval_csp(cal: pd.DataFrame, kind: str, complex_id: int, t1950: np.ndarray) -> np.ndarray:
    """Sum of the calibrations of one kind covering each time (NaN if none)."""
    t = np.asarray(t1950, dtype=float)
    out = np.full(t.shape, np.nan)
    sub = cal[(cal["kind"] == kind) & (cal["complex"] == complex_id)]
    for _, r in sub.iterrows():
        m = (t >= r["t0"]) & (t <= r["t1"])
        if not m.any():
            continue
        if r["form"] == "NRMPOW":
            x = 2 * (t[m] - r["t0"]) / (r["t1"] - r["t0"]) - 1
            v = np.polynomial.polynomial.polyval(x, r["coef"])
        elif r["form"] == "CONST":
            v = np.full(m.sum(), r["coef"][0])
        else:
            P, a = r["coef"][0], r["coef"][1:]
            x = 2 * np.pi * (t[m] - r["t0"]) / P
            v = a[0] + sum(a[2 * k - 1] * np.cos(k * x) + a[2 * k] * np.sin(k * x) for k in range(1, (len(a) - 1) // 2 + 1))
        out[m] = np.where(np.isnan(out[m]), v, out[m] + v)
    return out


#: Goldstone's seasonal zenith delays (m), TRK-2-23 Figure 3-2 (TRIG series
#: from 1972-01-01, period 31,557,600 s); the DSN calibrations in the ``tro``
#: files are corrections to these.
SEASONAL_C10 = {
    "WET": (31557600.0, [0.0870, -0.0360, -0.0336, 0.0002, 0.0200, 0.0008, -0.0021, -0.0036, -0.0002]),
    "DRY": (31557600.0, [2.0521, 0.0082, -0.0005, -0.0004, 0.0033, -0.0015, 0.0005, -0.0011, 0.0036]),
}


def _seasonal_c10(kind: str, t1950: np.ndarray) -> np.ndarray:
    P, a = SEASONAL_C10[kind]
    x = 2 * np.pi * (np.asarray(t1950) - _csp_time("72/01/01,00:00")) / P
    return a[0] + sum(a[2 * k - 1] * np.cos(k * x) + a[2 * k] * np.sin(k * x) for k in range(1, (len(a) - 1) // 2 + 1))


def saastamoinen_dry(lat_deg: float, height_m: float) -> float:
    """Zenith hydrostatic delay (m) for a standard-atmosphere surface pressure."""
    p_hpa = 1013.25 * (1 - 2.2557e-5 * height_m) ** 5.2559
    return 0.0022768 * p_hpa / (1 - 0.00266 * np.cos(2 * np.deg2rad(lat_deg)) - 0.28e-6 * height_m)


# ---------------------------------------------------------------------------
# AMC (water vapour radiometer) path delays at DSS 25
# ---------------------------------------------------------------------------

def read_amc(unit: int = 1) -> pd.DataFrame:
    """All AWVR path-delay records of one unit (1 or 2): t1950, zenith wet and
    dry delays (m), the AWVR elevation (deg, folded to 0-90) and flags. Only
    records with a valid retrieval (path-delay flag 0, AWVR flag 0) are kept."""
    rows = []
    for p in sorted(glob.glob(str(CASSINI_DIR / f"pd{unit}" / f"*.pd{unit}"))):
        for ln in Path(p).read_text().splitlines():
            if not ln.startswith("P "):
                continue
            f = ln.split()
            rows.append((float(f[1]), float(f[2]), float(f[3]), float(f[4]), float(f[5]), float(f[6]),
                         int(f[7]), *(int(x) for x in f[8:14])))
    df = pd.DataFrame(rows, columns=["t1970", "az", "el", "wet_cm", "dry_cm", "liq", "alg",
                                     "cloud", "wvr", "mtp", "surf", "wind", "pd"])
    df["t1950"] = df["t1970"] + _EPOCH_1970_IN_1950
    df["el_true"] = np.where(df["el"] > 90, 180 - df["el"], df["el"])
    ok = (df["pd"] == 0) & (df["wvr"] == 0) & (df["wet_cm"] < 900)
    df = df[ok].drop_duplicates("t1950").sort_values("t1950").reset_index(drop=True)
    df["wet_m"], df["dry_m"] = df["wet_cm"] / 100, df["dry_cm"] / 100
    return df


class Troposphere:
    """Zenith dry and wet delays at a station as functions of UTC (s since
    1950). ``source``: 'amc' (DSS 25 only; falls back to 'dsn' where the AWVRs
    have no valid data), 'dsn' (Goldstone only) or 'standard'."""

    def __init__(self, station: int, lat_deg: float, height_m: float, source: str = "amc",
                 amc_unit: int = 1, amc_smooth_s: float = 120.0):
        self.station, self.lat, self.h = station, lat_deg, height_m
        self.complex = COMPLEX[station]
        src = source
        if src == "amc" and station != 25:
            src = "dsn"
        if src == "dsn" and self.complex != 10:
            src = "standard"
        self.source = src
        self.cal = read_csp(sorted(glob.glob(str(CASSINI_DIR / "tro" / "*.tro")))) if src in ("amc", "dsn") else None
        self.amc = None
        if src == "amc":
            a = read_amc(amc_unit)
            if amc_smooth_s > 30:
                # running mean over amc_smooth_s within contiguous stretches
                w = int(round(amc_smooth_s / 30.0)) | 1
                grp = (a["t1950"].diff() > 90).cumsum()
                a["wet_m"] = a.groupby(grp)["wet_m"].transform(lambda s: s.rolling(w, center=True, min_periods=1).mean())
                a["dry_m"] = a.groupby(grp)["dry_m"].transform(lambda s: s.rolling(w, center=True, min_periods=1).mean())
            self.amc = a

    def _dsn(self, kind: str, t: np.ndarray) -> np.ndarray:
        corr = eval_csp(self.cal, kind, 10, t)
        return _seasonal_c10(kind, t) + np.nan_to_num(corr)

    def zenith(self, t1950) -> tuple[np.ndarray, np.ndarray]:
        t = np.asarray(t1950, dtype=float)
        if self.source == "standard":
            return np.full(t.shape, saastamoinen_dry(self.lat, self.h)), np.full(t.shape, 0.10)
        dry, wet = self._dsn("DRY", t), self._dsn("WET", t)
        if self.source == "amc":
            a = self.amc
            ta = a["t1950"].to_numpy()
            i = np.clip(np.searchsorted(ta, t), 1, len(ta) - 1)
            near = np.minimum(np.abs(t - ta[i - 1]), np.abs(t - ta[i])) <= 90.0
            # a gap of more than 3 minutes in the AWVR record: keep the DSN value there
            gap_ok = (ta[i] - ta[i - 1]) <= 120.0
            use = near & gap_ok
            w_amc = np.interp(t, ta, a["wet_m"].to_numpy())
            d_amc = np.interp(t, ta, a["dry_m"].to_numpy())
            wet = np.where(use, w_amc, wet)
            dry = np.where(use, d_amc, dry)
        return dry, wet

    def delay(self, t1950, el_rad) -> np.ndarray:
        """One-way tropospheric delay (s) along a line of sight at elevation el."""
        t = np.asarray(t1950, dtype=float)
        doy = (t / 86400.0 - (np.datetime64("2002-01-01") - np.datetime64("1950-01-01")).astype(int)) + 1.0
        m_h, m_w = niell(np.clip(el_rad, np.deg2rad(3.0), None), self.lat, self.h, doy)
        dry, wet = self.zenith(t)
        return (dry * m_h + wet * m_w) / C.C_LIGHT


class Ionosphere:
    """DSN line-of-sight ionosphere calibrations for Cassini (SCID 82)."""

    def __init__(self):
        self.cal = read_csp(sorted(glob.glob(str(CASSINI_DIR / "ion" / "*.ion"))))

    def phase_delay(self, station: int, t1950, freq_hz) -> np.ndarray:
        """One-way phase delay (s; negative: the plasma advances the phase)."""
        s = eval_csp(self.cal, "CHPART", COMPLEX[station], t1950)
        return -np.nan_to_num(s) * (F_ION_REF / np.asarray(freq_hz)) ** 2 / C.C_LIGHT
