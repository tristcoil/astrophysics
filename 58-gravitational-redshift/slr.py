"""Laser ranging as an independent check on the orbits (the main systematic of the test).

Galileo satellites carry laser retroreflectors. ILRS stations fire laser pulses
at them and time the echoes: a range measurement that needs no satellite clock.
ESA publishes, every week, per-pass mean residuals (observed minus computed
range, mm) of its own orbits. A radial orbit error dr leaks into the estimated
satellite clock as dr/c (the clock and the orbit are estimated together from
one-way signals), so a pass residual res says the clock estimate at that time
is off by about -res/c (Herrmann et al. 2018).

Fitting the same once-per-orbit template to these SLR-derived "clock errors"
as to the clocks gives the bias that orbit errors of this kind put into alpha.
"""

from __future__ import annotations

import datetime as dt
import gzip
import re

import numpy as np
import pandas as pd

import config as C
import prepare

LINE = re.compile(
    r"^\s*(?P<station>\d{4})\s+(?P<domes>\S+)\s+(?P<prn>[GREJC]\d\d)\s+"
    r"(?P<date>\d\d/\d\d/\d\d)\s+(?P<time>\d\d:\d\d)\s+(?P<dur>\d+)\s+"
    r"(?P<ngood>\d+)\s+(?P<mean>-?\d+)\s+(?P<std>\d+)\s")


def read_summaries() -> pd.DataFrame:
    rows = []
    for path in sorted(C.SLR_DIR.glob("*_SLR.SUM.gz")):
        with gzip.open(path, "rt", errors="replace") as fh:
            for line in fh:
                m = LINE.match(line)
                if not m:
                    continue
                g = m.groupdict()
                yy, mo, dd = (int(x) for x in g["date"].split("/"))
                hh, mi = (int(x) for x in g["time"].split(":"))
                start = dt.datetime(2000 + yy, mo, dd, hh, mi)
                rows.append({"file": path.name, "station": int(g["station"]), "prn": g["prn"],
                             "start": start, "dur_min": int(g["dur"]), "n_good": int(g["ngood"]),
                             "mean_mm": float(g["mean"]), "std_mm": float(g["std"])})
    df = pd.DataFrame(rows)
    # weekly files can repeat a pass at the week boundary: keep one copy
    df = df.drop_duplicates(subset=["station", "prn", "start", "dur_min", "n_good", "mean_mm"])
    return df.sort_values("start").reset_index(drop=True)


def gps_seconds(t: pd.Series) -> np.ndarray:
    """SLR epochs are UTC; the 17-18 s offset to GPS time is irrelevant at orbit-phase level."""
    epoch = pd.Timestamp(C.GPS_EPOCH)
    return ((pd.to_datetime(t) - epoch).dt.total_seconds()).to_numpy()


def pass_templates(prn: str, passes: pd.DataFrame, field: str = "kc") -> np.ndarray:
    """Average of a cached per-epoch quantity over each pass (1-min sampling, linear interpolation)."""
    d = prepare.load(prn)
    t_all = (d["gps0"][:, None] + np.arange(prepare.N_EPOCHS)[None, :] * C.EPOCH_STEP_S).ravel()
    v_all = d[field].astype(float).ravel()
    # drop the duplicated 24:00 epoch of each day (it is the next day's 00:00)
    keep = np.ones_like(t_all, bool)
    keep[prepare.N_EPOCHS - 1::prepare.N_EPOCHS] = False
    t_all, v_all = t_all[keep], v_all[keep]
    order = np.argsort(t_all)
    t_all, v_all = t_all[order], v_all[order]
    t0 = gps_seconds(passes["start"])
    out = np.full(len(passes), np.nan)
    for k, (ts, dur) in enumerate(zip(t0, passes["dur_min"].to_numpy())):
        ts_grid = ts + np.arange(0, max(dur, 1) * 60 + 1, 60.0)
        j = np.searchsorted(t_all, ts_grid)
        if j.min() == 0 or j.max() >= len(t_all):
            continue
        # require cached epochs around the whole pass (no day missing)
        if np.any(t_all[j] - t_all[j - 1] > C.EPOCH_STEP_S + 1):
            continue
        out[k] = np.mean(np.interp(ts_grid, t_all, v_all))
    return out


def orbit_bias(prn: str, passes: pd.DataFrame, *, max_dur_min: int = 120,
               floor_mm: float = 10.0, per_year_offset: bool = True) -> dict:
    """Regress SLR-implied clock errors (-res/c) on the alpha template (half the convention term).

    alpha's template is the gravitational half of the periodic proper time,
    T_alpha = 0.5 * (-2 r.v/c^2) for a Keplerian orbit, so the fitted slope is the
    alpha that orbit errors with this phase would fake.
    """
    p = passes[(passes["prn"] == prn) & (passes["dur_min"] <= max_dur_min)
               & (passes["n_good"] > 0)].copy()
    p["tpl"] = 0.5 * pass_templates(prn, p, "kc")
    p = p[np.isfinite(p["tpl"])]
    # robust: drop passes more than 5 robust sigma from the median (blunders, bad stations)
    med = p["mean_mm"].median()
    rs = 1.4826 * np.median(np.abs(p["mean_mm"] - med))
    p = p[np.abs(p["mean_mm"] - med) < 5 * rs]
    y = -p["mean_mm"].to_numpy() * 1e-3 / C.C_LIGHT          # s
    sig = np.sqrt((p["std_mm"].to_numpy() ** 2) / np.maximum(p["n_good"].to_numpy(), 1)
                  + floor_mm**2) * 1e-3 / C.C_LIGHT
    cols = [p["tpl"].to_numpy()]
    if per_year_offset:
        years = pd.to_datetime(p["start"]).dt.year.to_numpy()
        for yv in np.unique(years):
            cols.append((years == yv).astype(float))
    else:
        cols.append(np.ones(len(p)))
    X = np.vstack(cols).T
    W = 1 / sig**2
    A = X.T @ (X * W[:, None])
    b = X.T @ (W * y)
    cov = np.linalg.inv(A)
    coef = cov @ b
    res = y - X @ coef
    chi2 = float(np.sum(W * res**2) / max(len(y) - X.shape[1], 1))
    return {"prn": prn, "alpha_orbit": float(coef[0]), "err": float(np.sqrt(cov[0, 0])),
            "err_scaled": float(np.sqrt(cov[0, 0] * max(chi2, 1))), "chi2_nu": chi2,
            "n_passes": int(len(p)), "rms_mm": float(np.std(p["mean_mm"])),
            "first": str(pd.to_datetime(p["start"]).min().date()),
            "last": str(pd.to_datetime(p["start"]).max().date()), "passes": p}


# ---------------------------------------------------------------------------
# Pass-averaged geometry, and the orbit bias with the clock fit's nuisance model
# ---------------------------------------------------------------------------

def _stacked(d: dict, values: np.ndarray):
    """Continuous (t, value) series from the per-day cache, dropping each day's 24:00 duplicate."""
    t_all = (d["gps0"][:, None] + np.arange(prepare.N_EPOCHS)[None, :] * C.EPOCH_STEP_S).ravel()
    v_all = values.ravel()
    keep = np.ones_like(t_all, bool)
    keep[prepare.N_EPOCHS - 1::prepare.N_EPOCHS] = False
    order = np.argsort(t_all[keep])
    return t_all[keep][order], v_all[keep][order]


def pass_features(prn: str, passes: pd.DataFrame, harmonics=(1, 2)) -> pd.DataFrame:
    """Template and Sun geometry averaged over each pass (1-min sampling, linear interpolation)."""
    d = prepare.load(prn)
    mu = d["mu"].astype(float)
    series = {"tpl": 0.5 * d["kc"], "beta": d["beta"].astype(float)}
    for k in harmonics:
        series[f"cos{k}"] = np.cos(k * mu)
        series[f"sin{k}"] = np.sin(k * mu)
    stacked = {name: _stacked(d, v) for name, v in series.items()}
    t_all = stacked["tpl"][0]
    t0 = gps_seconds(passes["start"])
    out = {name: np.full(len(passes), np.nan) for name in series}
    for k, (ts, dur) in enumerate(zip(t0, passes["dur_min"].to_numpy())):
        grid = ts + np.arange(0, max(dur, 1) * 60 + 1, 60.0)
        j = np.searchsorted(t_all, grid)
        if j.min() == 0 or j.max() >= len(t_all):
            continue
        if np.any(t_all[j] - t_all[j - 1] > C.EPOCH_STEP_S + 1):
            continue
        for name, (tt, vv) in stacked.items():
            out[name][k] = np.mean(np.interp(grid, tt, vv))
    res = passes.copy()
    for name, v in out.items():
        res[name] = v
    return res[np.isfinite(res["tpl"])]


def clean_passes(p: pd.DataFrame, max_dur_min: int = 120) -> pd.DataFrame:
    p = p[(p["dur_min"] <= max_dur_min) & (p["n_good"] > 0)]
    med = p["mean_mm"].median()
    rs = 1.4826 * np.median(np.abs(p["mean_mm"] - med))
    return p[np.abs(p["mean_mm"] - med) < 5 * rs]


def orbit_bias_nuisance(p: pd.DataFrame, *, origin: str, window_days: float | None,
                        harmonics=(), beta_degree: int = 0, floor_mm: float = 10.0,
                        station_offsets: bool = True) -> dict:
    """alpha implied by SLR, with Sun-locked nuisance terms like the clock fit.

    Columns: the alpha template; per-station range offsets (SLR station biases are
    real and stations track at different orbit phases); optionally per-window
    Sun-locked harmonics times powers of beta.
    """
    y = -p["mean_mm"].to_numpy() * 1e-3 / C.C_LIGHT
    sig = np.sqrt(p["std_mm"].to_numpy() ** 2 / np.maximum(p["n_good"].to_numpy(), 1)
                  + floor_mm**2) * 1e-3 / C.C_LIGHT
    cols, names = [p["tpl"].to_numpy()], ["alpha"]
    if station_offsets:
        for st in np.unique(p["station"]):
            cols.append((p["station"].to_numpy() == st).astype(float))
            names.append(f"st{st}")
    else:
        cols.append(np.ones(len(p)))
        names.append("offset")
    if window_days and harmonics:
        days = (pd.to_datetime(p["start"]) - pd.Timestamp(origin)).dt.days.to_numpy()
        win = np.floor(days / window_days).astype(int)
        beta = np.radians(p["beta"].to_numpy())
        for k in harmonics:
            for trig in ("cos", "sin"):
                base = p[f"{trig}{k}"].to_numpy()
                for w in np.unique(win):
                    for b in range(beta_degree + 1):
                        cols.append(base * beta**b * (win == w))
                        names.append(f"{trig}{k}_w{w}_b{b}")
    X = np.vstack(cols).T
    W = 1 / sig**2
    A = X.T @ (X * W[:, None])
    cov = np.linalg.pinv(A)
    coef = cov @ (X.T @ (W * y))
    res = y - X @ coef
    chi2 = float(np.sum(W * res**2) / max(len(y) - X.shape[1], 1))
    return {"alpha_orbit": float(coef[0]), "err": float(np.sqrt(cov[0, 0])),
            "err_scaled": float(np.sqrt(cov[0, 0] * max(chi2, 1))), "chi2_nu": chi2,
            "n_passes": int(len(p)), "n_params": X.shape[1]}


def orbit_bias_bootstrap(p: pd.DataFrame, n_boot: int = 300, block_days: int = C.BOOTSTRAP_BLOCK_DAYS,
                         seed: int = C.RANDOM_SEED, **kw) -> float:
    """Block bootstrap (30-day blocks of passes) of orbit_bias_nuisance's alpha."""
    rng = np.random.default_rng(seed)
    days = (pd.to_datetime(p["start"]) - pd.to_datetime(p["start"]).min()).dt.days.to_numpy()
    blocks = days // block_days
    groups = [np.flatnonzero(blocks == b) for b in np.unique(blocks)]
    vals = []
    for _ in range(n_boot):
        pick = np.concatenate([groups[g] for g in rng.integers(0, len(groups), len(groups))])
        vals.append(orbit_bias_nuisance(p.iloc[pick], **kw)["alpha_orbit"])
    return float(np.std(vals))
