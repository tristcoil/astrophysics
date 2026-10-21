"""X-ray outbursts: daily MAXI and Swift/BAT light curves in Crab units,
hardness ratios, state transitions and Eddington ratios.

Every flux is first divided by the Crab Nebula's flux measured by the same
instrument in the same band over the same year (a running median), so that
"1 Crab" means the same thing in every band and the instruments' slow drifts
cancel. A Crab-unit flux becomes an energy flux by multiplying with the
Crab's own energy flux in that band, from its canonical spectrum
N(E) = 9.7 E^-2.1 photons cm^-2 s^-1 keV^-1 (Toor & Seward 1974). That
assumes the source looks like the Crab within each band: rough for very
soft spectra, adequate for an order-of-magnitude Eddington ratio.

The hardness ratio is (4-10 keV)/(2-4 keV) in Crab units: 1 for a
Crab-like spectrum, above 1 in the hard state, near 0.2 in the soft state.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from config import (
    CRAB_INDEX, CRAB_NORM, CRAB_WINDOW_D, GX339_MIN_OUTBURST_DAYS,
    GX339_OUTBURST_THRESHOLD_CRAB, J1820_DISTANCE_KPC, J1820_MASS, J1820_OUTBURST,
    KPC_CM, L_EDD_PER_MSUN,
)
from fetch_data import bat_path, maxi_path

KEV_ERG = 1.602176634e-9
MAXI_BANDS = {"2-4": (2.0, 4.0), "4-10": (4.0, 10.0), "10-20": (10.0, 20.0)}
BAT_BAND = (15.0, 50.0)


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_maxi(source: str) -> pd.DataFrame:
    cols = ["mjd", "f2_20", "e2_20", "f2_4", "e2_4", "f4_10", "e4_10", "f10_20", "e10_20"]
    df = pd.DataFrame(np.loadtxt(maxi_path(source)), columns=cols)
    # days without a usable scan carry zero or absurd errors
    ok = (df["e2_4"] > 0) & (df["e4_10"] > 0) & (df["e10_20"] > 0) & (df["e2_4"] < 1.0)
    return df[ok].reset_index(drop=True)


def load_bat(source: str) -> pd.DataFrame:
    raw = np.loadtxt(bat_path(source), comments="#")
    df = pd.DataFrame({"mjd": raw[:, 0] + 0.5, "rate": raw[:, 1], "err": raw[:, 2],
                       "flag": raw[:, 7].astype(int)})
    return df[(df["flag"] == 0) & (df["err"] > 0)].reset_index(drop=True)


# ---------------------------------------------------------------------------
# The Crab as the yardstick
# ---------------------------------------------------------------------------

def crab_energy_flux(e1: float, e2: float) -> float:
    """Energy flux (erg cm^-2 s^-1) of the canonical Crab spectrum in [e1, e2] keV."""
    a = 2.0 - CRAB_INDEX
    return CRAB_NORM * (e2**a - e1**a) / a * KEV_ERG


def _running_median(mjd, value, err, grid, window=CRAB_WINDOW_D):
    out = np.full(len(grid), np.nan)
    for k, t in enumerate(grid):
        m = np.abs(mjd - t) < window / 2
        good = m & np.isfinite(value) & (err > 0)
        if good.sum() > 20:
            out[k] = np.median(value[good])
    # fill edges by nearest valid value
    valid = np.isfinite(out)
    return np.interp(grid, grid[valid], out[valid])


class CrabYardstick:
    """Running yearly median of the Crab in every MAXI band and in BAT."""

    def __init__(self):
        maxi, bat = load_maxi("Crab"), load_bat("Crab")
        self.grid = np.arange(min(maxi["mjd"].min(), bat["mjd"].min()),
                              max(maxi["mjd"].max(), bat["mjd"].max()) + 30, 30.0)
        self.curves = {}
        for band in ("2-4", "4-10", "10-20"):
            key = "f" + band.replace("-", "_")
            self.curves[band] = _running_median(maxi["mjd"].values, maxi[key].values,
                                                maxi["e" + key[1:]].values, self.grid)
        self.curves["2-20"] = _running_median(maxi["mjd"].values, maxi["f2_20"].values,
                                              maxi["e2_20"].values, self.grid)
        self.curves["bat"] = _running_median(bat["mjd"].values, bat["rate"].values,
                                             bat["err"].values, self.grid)

    def __call__(self, band: str, mjd):
        return np.interp(mjd, self.grid, self.curves[band])

    def summary(self) -> dict:
        return {b: {"median": float(np.median(c)), "min": float(c.min()), "max": float(c.max())}
                for b, c in self.curves.items()}


def light_curves(source: str, crab: CrabYardstick, mjd_range=None) -> pd.DataFrame:
    """MAXI bands and BAT on a common daily grid, in Crab units and erg/cm^2/s."""
    maxi, bat = load_maxi(source), load_bat(source)
    if mjd_range is not None:
        maxi = maxi[(maxi["mjd"] >= mjd_range[0]) & (maxi["mjd"] <= mjd_range[1])]
        bat = bat[(bat["mjd"] >= mjd_range[0]) & (bat["mjd"] <= mjd_range[1])]
    df = pd.DataFrame({"mjd": np.round(maxi["mjd"].values - 0.5) + 0.5})
    for band, (e1, e2) in MAXI_BANDS.items():
        key = "f" + band.replace("-", "_")
        c = crab(band, maxi["mjd"].values)
        df[f"crab_{band}"] = maxi[key].values / c
        df[f"err_{band}"] = maxi["e" + key[1:]].values / c
        df[f"eflux_{band}"] = df[f"crab_{band}"] * crab_energy_flux(e1, e2)
    df["crab_2-20"] = maxi["f2_20"].values / crab("2-20", maxi["mjd"].values)
    df["eflux_2-20"] = df[[f"eflux_{b}" for b in MAXI_BANDS]].sum(axis=1)
    b = pd.DataFrame({"mjd": np.round(bat["mjd"].values - 0.5) + 0.5,
                      "crab_bat": bat["rate"].values / crab("bat", bat["mjd"].values),
                      "err_bat": bat["err"].values / crab("bat", bat["mjd"].values)})
    df = df.merge(b, on="mjd", how="outer").sort_values("mjd").reset_index(drop=True)
    df["eflux_bat"] = df["crab_bat"] * crab_energy_flux(*BAT_BAND)
    # 20-50 keV from BAT (the Crab's 20-50/15-50 share), so 2-50 keV has no overlap
    share = crab_energy_flux(20.0, 50.0) / crab_energy_flux(*BAT_BAND)
    df["eflux_2-50"] = df["eflux_2-20"] + share * df["eflux_bat"]
    with np.errstate(divide="ignore", invalid="ignore"):
        df["hr"] = df["crab_4-10"] / df["crab_2-4"]
        df["hr_err"] = df["hr"] * np.hypot(df["err_4-10"] / df["crab_4-10"],
                                           df["err_2-4"] / df["crab_2-4"])
    good_hr = (df["crab_2-4"] > 3 * df["err_2-4"]) & (df["crab_4-10"] > 3 * df["err_4-10"])
    df.loc[~good_hr, ["hr", "hr_err"]] = np.nan
    return df


# ---------------------------------------------------------------------------
# MAXI J1820+070: states, transitions, Eddington ratios
# ---------------------------------------------------------------------------

def _smooth(y, width=3):
    s = pd.Series(y).rolling(width, center=True, min_periods=1).median()
    return s.values


def hardness_levels(df: pd.DataFrame, min_crab: float = 0.05) -> tuple[float, float]:
    """Hard and soft hardness levels: medians of the 3-day running median of the
    hardness in the hardest and softest thirds of the bright days."""
    bright = df[np.isfinite(df["hr"]) & (df["crab_2-20"] > min_crab)]
    hr = _smooth(bright["hr"].values)
    lo_cut, hi_cut = np.quantile(hr, [1 / 3, 2 / 3])
    return float(np.median(hr[hr >= hi_cut])), float(np.median(hr[hr <= lo_cut]))


def find_transitions(df: pd.DataFrame, min_crab: float = 0.05,
                     levels: tuple[float, float] | None = None) -> dict | None:
    """Hard-to-soft and soft-to-hard transitions from the hardness ratio.

    A transition is the day the 3-day running median of the hardness crosses
    the geometric mean of the hard and soft levels (the ratio is
    multiplicative). The levels come from this outburst's own bright days
    (2-20 keV above ``min_crab``) unless given. The transition 'window' is
    where the hardness is between 80% and 20% of the way from hard to soft,
    in log space."""
    bright = df[np.isfinite(df["hr"]) & (df["crab_2-20"] > min_crab)].copy()
    if len(bright) < 15:
        return None
    hr = _smooth(bright["hr"].values)
    hard, soft = levels if levels is not None else hardness_levels(df, min_crab)
    mid = np.sqrt(hard * soft)
    t = bright["mjd"].values
    n = len(t)
    below = hr < mid
    if not below.any() or below.all():
        return None
    first_soft = int(np.argmax(below))
    last_soft = int(n - 1 - np.argmax(below[::-1]))
    lg = (np.log(hr) - np.log(soft)) / (np.log(hard) - np.log(soft))  # 1 hard, 0 soft

    def window(i0, step):
        i0 = min(max(i0, 0), n - 1)
        a = i0
        while 0 < a < n - 1 and lg[a] < 0.8:
            a -= step
        b = i0
        while 0 < b < n - 1 and lg[b] > 0.2:
            b += step
        a, b = min(max(a, 0), n - 1), min(max(b, 0), n - 1)
        return sorted([float(t[a]), float(t[b])])

    returned = last_soft + 1 < n  # did the source come back to the hard state?
    return {
        "hr_hard": hard, "hr_soft": soft, "hr_mid": float(mid),
        "hard_to_soft_mjd": float(t[first_soft]),
        "hard_to_soft_window": window(first_soft, 1),
        "soft_to_hard_mjd": float(t[last_soft + 1]) if returned else float("nan"),
        "soft_to_hard_window": window(last_soft + 1, -1) if returned else None,
    }


def eddington_ratio_samples(eflux_cgs, rng, n: int = 50000):
    """L/L_Edd for MAXI J1820+070 over distance and mass uncertainties."""
    d = rng.normal(*J1820_DISTANCE_KPC, n) * KPC_CM
    u = rng.standard_normal(n)
    m = J1820_MASS[0] + np.where(u > 0, J1820_MASS[1] * u, J1820_MASS[2] * u)
    return 4 * np.pi * d**2 * eflux_cgs / (L_EDD_PER_MSUN * m)


def j1820_analysis(crab: CrabYardstick, rng) -> tuple[dict, pd.DataFrame]:
    df = light_curves("MAXI J1820+070", crab, J1820_OUTBURST)
    tr = find_transitions(df)

    def flux_at(mjd, col, half=1.0):
        m = (np.abs(df["mjd"] - mjd) <= half) & np.isfinite(df[col])
        return float(df.loc[m, col].mean())

    i_peak = int(np.nanargmax(df["eflux_2-50"].values))
    i_bat = int(np.nanargmax(df["crab_bat"].values))
    i_maxi = int(np.nanargmax(df["crab_2-4"].values))
    peak = df["eflux_2-50"].iloc[i_peak]
    s2h = flux_at(tr["soft_to_hard_mjd"], "eflux_2-50")
    h2s = flux_at(tr["hard_to_soft_mjd"], "eflux_2-50")

    def q(x):
        lo, med, hi = np.quantile(x, [0.16, 0.5, 0.84])
        return {"median": float(med), "plus": float(hi - med), "minus": float(med - lo)}

    out = {
        "transitions": tr,
        "peak_2_50_kev": {"mjd": float(df["mjd"].iloc[i_peak]), "eflux_cgs": float(peak),
                          "l_ledd": q(eddington_ratio_samples(peak, rng))},
        "peak_bat": {"mjd": float(df["mjd"].iloc[i_bat]), "crab": float(df["crab_bat"].iloc[i_bat])},
        "peak_maxi_2_4": {"mjd": float(df["mjd"].iloc[i_maxi]),
                          "crab": float(df["crab_2-4"].iloc[i_maxi])},
        "soft_to_hard": {"eflux_2_50_cgs": s2h, "l_ledd": q(eddington_ratio_samples(s2h, rng))},
        "hard_to_soft": {"eflux_2_50_cgs": h2s, "l_ledd": q(eddington_ratio_samples(h2s, rng))},
    }
    return out, df


# ---------------------------------------------------------------------------
# GX 339-4: every outburst of the MAXI era
# ---------------------------------------------------------------------------

def gx339_outbursts(crab: CrabYardstick) -> tuple[list[dict], pd.DataFrame]:
    """Outbursts = runs of days where the 15-day running median of the MAXI
    2-20 keV flux exceeds the threshold, lasting at least the minimum length
    (gaps of up to 20 days bridged), padded by 30 days on each side. Each is
    classified 'full' if its hardness reaches the soft state (below 0.4, well
    into the soft branch), else 'hard-only' (a 'failed' outburst)."""
    df = light_curves("GX 339-4", crab)
    df = df[np.isfinite(df["crab_2-20"])].reset_index(drop=True)
    sm = pd.Series(df["crab_2-20"].values).rolling(15, center=True, min_periods=5).median().values
    on = sm > GX339_OUTBURST_THRESHOLD_CRAB
    t = df["mjd"].values
    runs, start = [], None
    for k in range(len(t)):
        if on[k] and start is None:
            start = k
        if start is not None and (not on[k] or k == len(t) - 1):
            runs.append([start, k if on[k] else k - 1])
            start = None
    merged = []
    for a, b in runs:
        if merged and t[a] - t[merged[-1][1]] <= 20:
            merged[-1][1] = b
        else:
            merged.append([a, b])
    out, segs = [], []
    for a, b in merged:
        if t[b] - t[a] < GX339_MIN_OUTBURST_DAYS:
            continue
        t0, t1 = t[a] - 30, t[b] + 30
        seg = df[(df["mjd"] >= t0) & (df["mjd"] <= t1)]
        hr = _smooth(seg["hr"].values, 5)
        softest = float(np.nanmin(hr)) if np.isfinite(hr).any() else np.nan
        out.append({
            "start_mjd": float(t[a]), "end_mjd": float(t[b]),
            "window_mjd": [float(t0), float(t1)],
            "peak_crab_2_20": float(seg["crab_2-20"].max()),
            "peak_crab_bat": float(np.nanmax(seg["crab_bat"].values)),
            "min_hr": softest,
            "kind": "full" if softest < 0.4 else "hard-only",
        })
        segs.append(seg)
    # One pair of hard/soft levels for the source, from all its full outbursts,
    # so that every outburst's transitions are found the same way (a faint
    # outburst's hard days are too faint for a hardness of their own).
    full = pd.concat([sg for sg, o in zip(segs, out) if o["kind"] == "full"])
    levels = hardness_levels(full, GX339_OUTBURST_THRESHOLD_CRAB)
    for o, seg in zip(out, segs):
        if o["kind"] == "full":
            tr = find_transitions(seg, min_crab=GX339_OUTBURST_THRESHOLD_CRAB / 2, levels=levels)
            o["transitions"] = tr
            for key in ("hard_to_soft", "soft_to_hard"):
                mjd = tr[f"{key}_mjd"] if tr else np.nan
                m = (np.abs(seg["mjd"] - mjd) <= 1) & np.isfinite(seg["crab_2-20"])
                o[f"{key}_crab_2_20"] = float(seg.loc[m, "crab_2-20"].mean()) if m.any() else np.nan
                o[f"{key}_eflux_2_50"] = float(seg.loc[m, "eflux_2-50"].mean()) if m.any() else np.nan
    return out, df
