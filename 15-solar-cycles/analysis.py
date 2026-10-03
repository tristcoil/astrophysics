"""Every measurement the README quotes, computed from the cached data.

    tables, results = analysis.analyse()

``tables`` holds what the figures draw (the smoothed record, the cycle table, the
sunspot groups with their cycle assignments, the rotation-by-rotation centroids,
the rotation pairs ...). ``results`` holds every number the README quotes and is
saved to ``data/results.json``, so each one can be traced to the run that made it.

Seven parts, one function each:

1. record()       the sunspot number since 1749: our smoothing and cycle finder,
                  checked against SILSO's own; cycle lengths and amplitudes;
                  spotless days at each minimum
2. rules()        the empirical rules of the cycle (Waldmeier, Gnevyshev-Ohl, ...)
3. cycle25()      the current cycle against the official forecast, and what its
                  early rise could have told us
4. hemispheres()  the two hemispheres keep their own, out-of-step cycles
5. butterfly()    sunspot groups since 1874: the butterfly diagram, the 1977 data
                  seam, and the equatorward drift of the sunspot zones
6. rotation()     the Sun's differential rotation from day-to-day group motions
7. long_record()  before 1749: the group-number reconstructions and the Maunder
                  minimum
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import signal

import groups
import silso
import solar_physics as sp
from config import (
    AREA_SEAM_WINDOW_YR,
    BUTTERFLY_LAT_BIN_DEG,
    BUTTERFLY_LAT_MAX_DEG,
    CARRINGTON_EPOCH_JD,
    CARRINGTON_SYNODIC_DAYS,
    CENTROID_MIN_AREA,
    CYCLE_EARLY_MONTHS,
    CYCLE_PRIOR_FLOOR,
    DALTON_MINIMUM,
    FIRST_NUMBERED_MINIMUM_YEAR,
    GROUP_NUMBER_SERIES,
    LIT_AREA_FACTOR,
    LIT_BALTHASAR_1986,
    LIT_HATHAWAY_DRIFT,
    LIT_HATHAWAY_T0,
    LIT_HS98_SCALE_V2,
    LIT_PANEL_2019,
    LIT_SILSO_CYCLE25_MAX,
    MAUNDER_MINIMUM,
    MIN_CYCLE_SEPARATION_YR,
    PREDICTION_CALIBRATION_CYCLES,
    PREDICTION_LEAD_MONTHS,
    RESULTS_PATH,
    RGO_LAST_YEAR,
    ROTATION_LAT_BIN_DEG,
    ROTATION_MAX_CMD_DEG,
    ROTATION_MAX_DRIFT_DEG_PER_DAY,
    ROTATION_PAIR_DT_DAYS,
    ROTATION_YEARS,
    USAF_AREA_FACTOR_HATHAWAY,
)

MONTH0_YEAR = 1749  # SILSO's monthly record begins January 1749
#: Drift-law fits and cycle statistics from the group data use these cycles: 12 is
#: the first whose whole life the RGO record covers (it began in 1878).
FIRST_FULL_GROUP_CYCLE = 12
N_BOOT = 2000


# --- Small helpers ------------------------------------------------------------------------


def _ym(t: float) -> str:
    """Decimal year (any time within a month) -> 'YYYY-MM'."""
    k = int(np.floor((t - MONTH0_YEAR) * 12 + 1e-6))
    return f"{MONTH0_YEAR + k // 12}-{k % 12 + 1:02d}"


def _t_of_ym(s: str) -> float:
    """'YYYY-MM' -> decimal year at mid-month."""
    y, m = map(int, s.split("-"))
    return y + (m - 0.5) / 12.0


def _months(t) -> np.ndarray:
    """Decimal years -> months since 1749-01-01."""
    return (np.asarray(t, dtype=float) - MONTH0_YEAR) * 12.0


def _month_diff(a: str, b: str) -> int:
    ya, ma = map(int, a.split("-"))
    yb, mb = map(int, b.split("-"))
    return (ya - yb) * 12 + (ma - mb)


def _decimal_year_of_jd(jd) -> np.ndarray:
    return sp.decimal_year(sp.datetime_from_jd(jd))


def _clean(obj):
    """Make results JSON-serialisable (numpy scalars/arrays, NaN -> None)."""
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [_clean(v) for v in obj.tolist()]
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        return None if not np.isfinite(obj) else float(obj)
    return obj


# --- 1. The sunspot record since 1749 -------------------------------------------------------


def record() -> tuple[dict, dict]:
    m = silso.monthly()
    sn = m["sn"].to_numpy(float)
    t = m["t"].to_numpy(float)
    ours = sp.smooth_13_month(sn)
    m["smoothed"] = ours
    m["smoothed_silso"] = silso.smoothed()["sn"].to_numpy(float)

    diff = ours - m["smoothed_silso"].to_numpy()
    ok = np.isfinite(diff)
    smoothing_check = {
        "n_months": int(ok.sum()),
        "max_abs_diff": float(np.max(np.abs(diff[ok]))),
        "first_month": _ym(t[ok][0]), "last_month": _ym(t[ok][-1]),
    }

    ex = sp.find_cycle_extrema(ours, min_separation_yr=MIN_CYCLE_SEPARATION_YR)
    if int(np.floor(t[ex.minima[0]])) != FIRST_NUMBERED_MINIMUM_YEAR:
        raise RuntimeError("the first minimum found is not 1755's; cycle numbering would be off")
    n_cyc = len(ex.minima)
    tm = _months(t)

    rows = []
    for k in range(n_cyc):
        i, j = int(ex.minima[k]), int(ex.maxima[k])
        nxt = int(ex.minima[k + 1]) if k + 1 < n_cyc else None
        end = nxt if nxt is not None else len(sn)
        shape = sp.fit_cycle_shape(tm[i:end], sn[i:end], t0_guess=tm[i] - 3,
                                   rmax_guess=ours[j], t0_bounds=(tm[i] - 30, tm[i] + 18))
        rows.append(dict(
            cycle=k + 1, i_min=i, i_max=j, t_min=t[i], t_max=t[j],
            min_ym=_ym(t[i]), max_ym=_ym(t[j]), min_sn=ours[i], max_sn=ours[j],
            rise_yr=t[j] - t[i],
            length_yr=(t[nxt] - t[i]) if nxt is not None else np.nan,
            decline_yr=(t[nxt] - t[j]) if nxt is not None else np.nan,
            total=np.nansum(sn[i:end]), complete=nxt is not None,
            s24=ours[i + 24] if i + 24 < len(ours) else np.nan,
            rise_20_80_yr=_rise_between(t[i:j + 1], ours[i:j + 1], ours[j], 0.2, 0.8),
            t0=MONTH0_YEAR + shape["t0_months"] / 12.0,
            t0_err_months=shape["t0_err_months"], shape_rmax=shape["rmax"],
            shape_b_months=shape["b_months"], shape_rms=shape["rms_residual"],
        ))
    cycles = pd.DataFrame(rows)
    last = cycles.iloc[-1]
    last_smoothed = int(np.flatnonzero(np.isfinite(ours))[-1])
    # The current cycle's maximum is settled once the smoothed number has fallen
    # well below it (SILSO waits for the same before declaring a maximum).
    fall = 1.0 - ours[last_smoothed] / last["max_sn"]

    # Against SILSO's official table.
    off = silso.official_cycles()
    comparison = []
    for _, o in off.iterrows():
        c = cycles.loc[cycles["cycle"] == o["cycle"]].iloc[0]
        row = {"cycle": int(o["cycle"]),
               "min_official": f"{int(o['min_year'])}-{int(o['min_month']):02d}",
               "min_ours": c["min_ym"], "min_sn_official": float(o["min_sn"]),
               "min_sn_ours": round(float(c["min_sn"]), 2)}
        row["min_offset_months"] = _month_diff(c["min_ym"], row["min_official"])
        if pd.notna(o["max_year"]):
            row["max_official"] = f"{int(o['max_year'])}-{int(o['max_month']):02d}"
            row["max_sn_official"] = float(o["max_sn"])
        elif int(o["cycle"]) == int(last["cycle"]):
            row["max_official"] = LIT_SILSO_CYCLE25_MAX["date"]
            row["max_sn_official"] = LIT_SILSO_CYCLE25_MAX["value"]
        row["max_ours"] = c["max_ym"]
        row["max_sn_ours"] = round(float(c["max_sn"]), 2)
        row["max_offset_months"] = _month_diff(c["max_ym"], row["max_official"])
        comparison.append(row)
    comp = pd.DataFrame(comparison)

    done = cycles[cycles["complete"]]
    i_long, i_short = done["length_yr"].idxmax(), done["length_yr"].idxmin()
    i_big, i_small = cycles["max_sn"].idxmax(), cycles["max_sn"].idxmin()
    dalton = cycles["cycle"].between(5, 7)
    in_dalton = (t >= DALTON_MINIMUM[0]) & (t < DALTON_MINIMUM[1])

    results = {
        "months": len(m), "first_month": _ym(t[0]), "last_month": _ym(t[-1]),
        "provisional_months": m.loc[m["definitive"] == 0, "date"].dt.strftime("%Y-%m").tolist(),
        "smoothing_check": smoothing_check,
        "n_cycles": n_cyc,
        "official_comparison": {
            "rows": comp.to_dict(orient="records"),
            "minima_exact": int((comp["min_offset_months"] == 0).sum()),
            "maxima_exact": int((comp["max_offset_months"] == 0).sum()),
            "max_abs_min_offset_months": int(comp["min_offset_months"].abs().max()),
        },
        "lengths": {
            "n": len(done), "mean_yr": done["length_yr"].mean(), "sd_yr": done["length_yr"].std(),
            "longest": {"cycle": int(cycles.at[i_long, "cycle"]), "yr": cycles.at[i_long, "length_yr"]},
            "shortest": {"cycle": int(cycles.at[i_short, "cycle"]), "yr": cycles.at[i_short, "length_yr"]},
            "mean_rise_yr": cycles["rise_yr"].mean(), "mean_decline_yr": done["decline_yr"].mean(),
        },
        "amplitudes": {
            "strongest": {"cycle": int(cycles.at[i_big, "cycle"]), "max_sn": cycles.at[i_big, "max_sn"],
                          "date": cycles.at[i_big, "max_ym"]},
            "weakest": {"cycle": int(cycles.at[i_small, "cycle"]), "max_sn": cycles.at[i_small, "max_sn"],
                        "date": cycles.at[i_small, "max_ym"]},
            "mean_max_sn": cycles["max_sn"].mean(),
            "dalton_cycles_mean_max": cycles.loc[dalton, "max_sn"].mean(),
            "other_cycles_mean_max": cycles.loc[~dalton, "max_sn"].mean(),
            "dalton_mean_sn": float(np.mean(sn[in_dalton])),
            "overall_mean_sn": float(np.mean(sn)),
        },
        "current_cycle": {
            "cycle": int(last["cycle"]), "min": last["min_ym"], "max": last["max_ym"],
            "max_sn": last["max_sn"], "last_smoothed_month": _ym(t[last_smoothed]),
            "last_smoothed_sn": float(ours[last_smoothed]), "fall_from_max": float(fall),
        },
        "shape_start_times": {int(r.cycle): {"t0": _ym(r.t0), "offset_from_min_months":
                                              (r.t0 - r.t_min) * 12.0, "rmax": r.shape_rmax}
                              for r in cycles.itertuples()},
        "periodogram": _periodogram(sn),
        "spotless_days": _spotless_days(cycles),
    }
    # Our start times against Hathaway (2011), who used version-1 sunspot numbers.
    offs = [(_months(r.t0) - _months(_t_of_ym(LIT_HATHAWAY_T0[r.cycle])))[()]
            for r in cycles.itertuples() if r.cycle in LIT_HATHAWAY_T0]
    results["shape_start_times_vs_hathaway"] = {
        "n": len(offs), "mean_offset_months": float(np.mean(offs)),
        "sd_offset_months": float(np.std(offs, ddof=1)),
        "range_months": [float(np.min(offs)), float(np.max(offs))]}

    tables = {"monthly": m, "cycles": cycles, "official": off, "extrema": ex,
              "yearly": silso.yearly(), "daily": silso.daily()}
    return tables, results


def _rise_between(t, s, peak: float, lo: float, hi: float) -> float:
    """Years the smoothed number takes to climb from lo x peak to hi x peak."""
    t_lo = t[np.argmax(s >= lo * peak)]
    t_hi = t[np.argmax(s >= hi * peak)]
    return float(t_hi - t_lo)


def _periodogram(sn: np.ndarray) -> dict:
    """Power spectrum of the monthly record; the peak between 5 and 20 years."""
    x = np.asarray(sn, dtype=float)
    f, p = signal.periodogram(x - x.mean(), fs=12.0, window="hann", detrend="linear",
                              nfft=1 << 16)
    sel = (f > 1 / 20) & (f < 1 / 5)
    k = np.argmax(np.where(sel, p, -np.inf))
    half = p[k] / 2
    lo = k
    while lo > 0 and p[lo] > half:
        lo -= 1
    hi = k
    while hi < len(p) - 1 and p[hi] > half:
        hi += 1
    return {"peak_period_yr": float(1 / f[k]), "half_power_periods_yr": [float(1 / f[hi]), float(1 / f[lo])],
            "freq_per_yr": f.tolist()[::64], "power": (p / p[k]).tolist()[::64]}


def _spotless_days(cycles: pd.DataFrame, min_coverage: float = 0.99) -> list[dict]:
    """Days with no spots at all around each minimum: between the maxima either
    side of it. Only minima whose window was observed on >= 99% of days."""
    d = silso.daily()
    t = d["t"].to_numpy(float)
    s = d["sn"].to_numpy(float)
    out = []
    for k in range(1, len(cycles)):
        lo, hi = cycles.at[k - 1, "t_max"], cycles.at[k, "t_max"]
        w = (t >= lo) & (t < hi)
        seen = w & np.isfinite(s)
        if w.sum() == 0 or seen.sum() / w.sum() < min_coverage:
            continue
        zero = seen & (s == 0)
        # longest unbroken run of spotless days
        z = zero[w].astype(int)
        edges = np.diff(np.r_[0, z, 0])
        runs = np.flatnonzero(edges == -1) - np.flatnonzero(edges == 1)
        out.append({"minimum_of_cycle": int(cycles.at[k, "cycle"]), "minimum": cycles.at[k, "min_ym"],
                    "spotless_days": int(zero.sum()), "longest_run_days": int(runs.max(initial=0)),
                    "coverage": float(seen.sum() / w.sum())})
    return out


# --- 2. The rules of the cycle -----------------------------------------------------------------


def rules(cycles: pd.DataFrame) -> dict:
    c = cycles
    done = c[c["complete"]]
    modern = c["cycle"] >= 10
    out = {
        # Waldmeier (1935): strong cycles rise fast. Every cycle's rise is complete.
        "waldmeier": sp.correlation(c["rise_yr"], c["max_sn"]),
        "waldmeier_from_cycle10": sp.correlation(c.loc[modern, "rise_yr"], c.loc[modern, "max_sn"]),
        "waldmeier_fit": sp.line_fit(c["rise_yr"], c["max_sn"]),
        # The same effect with rise time measured from 20% to 80% of the peak, which
        # ignores how deep the minimum was and how flat the top is.
        "waldmeier_20_80": sp.correlation(c["rise_20_80_yr"], c["max_sn"]),
        # The rate of rise: activity two years after minimum against the peak to come.
        "early_rise": sp.correlation(done["s24"], done["max_sn"]),
        "early_rise_from_cycle10": sp.correlation(done.loc[done["cycle"] >= 10, "s24"],
                                                  done.loc[done["cycle"] >= 10, "max_sn"]),
        # Amplitude-period effect: a long cycle is followed by a weak one.
        "length_vs_next_amplitude": sp.correlation(c["length_yr"].to_numpy()[:-1],
                                                   c["max_sn"].to_numpy()[1:]),
        "length_vs_next_fit": sp.line_fit(c["length_yr"].to_numpy()[:-1], c["max_sn"].to_numpy()[1:]),
        # The activity left at minimum foreshadows the next cycle (weakly).
        "minimum_vs_amplitude": sp.correlation(c["min_sn"], c["max_sn"]),
    }
    # Gnevyshev-Ohl rule: summed over a whole cycle, an odd cycle outdoes the even
    # cycle before it.
    pairs = []
    for even in range(2, int(c["cycle"].max()), 2):
        e = c.loc[c["cycle"] == even].iloc[0]
        o = c.loc[c["cycle"] == even + 1].iloc[0]
        pairs.append({"even": even, "odd": even + 1, "even_total": e["total"],
                      "odd_total": o["total"], "odd_complete": bool(o["complete"]),
                      "holds": bool(o["total"] > e["total"]),
                      "even_max": e["max_sn"], "odd_max": o["max_sn"]})
    out["gnevyshev_ohl"] = {
        "pairs": pairs,
        "n_complete": sum(p["odd_complete"] for p in pairs),
        "n_hold": sum(p["holds"] for p in pairs if p["odd_complete"]),
        "violations": [f"{p['even']}-{p['odd']}" for p in pairs if p["odd_complete"] and not p["holds"]],
    }
    return out


# --- 3. Cycle 25 against the forecast ---------------------------------------------------------


def cycle25(monthly: pd.DataFrame, cycles: pd.DataFrame) -> tuple[dict, dict]:
    c25 = cycles.iloc[-1]
    ours = monthly["smoothed"].to_numpy(float)
    t = monthly["t"].to_numpy(float)
    panel = silso.panel_forecast()
    current = silso.current_forecast()
    lit = LIT_PANEL_2019

    panel_peak_row = panel.loc[panel["ssn"].idxmax()]
    res = {
        "observed_max": c25["max_sn"], "observed_max_date": c25["max_ym"],
        "observed_min": c25["min_sn"], "observed_min_date": c25["min_ym"],
        "panel_max": lit["max_ssn"], "panel_max_range": list(lit["max_range"]),
        "panel_max_date": lit["max_date"], "panel_curve_peak": float(panel_peak_row["ssn"]),
        "panel_curve_peak_date": panel_peak_row["date"].strftime("%Y-%m"),
        "max_ratio_to_panel": c25["max_sn"] / lit["max_ssn"],
        "max_months_early": -_month_diff(c25["max_ym"], lit["max_date"]),
        "max_timing_within_panel_range": abs(_month_diff(c25["max_ym"], lit["max_date"]))
        <= lit["max_date_pm_months"],
        "max_within_panel_range": lit["max_range"][0] <= c25["max_sn"] <= lit["max_range"][1],
        "min_months_early": -_month_diff(c25["min_ym"], lit["min_date"]),
        "min_timing_within_panel_range": abs(_month_diff(c25["min_ym"], lit["min_date"]))
        <= lit["min_date_pm_months"],
        "silso_announced": LIT_SILSO_CYCLE25_MAX,
        "vs_cycle24": c25["max_sn"] / cycles.iloc[-2]["max_sn"],
        "vs_cycle23": c25["max_sn"] / cycles.iloc[-3]["max_sn"],
        "rank": int((cycles["max_sn"] > c25["max_sn"]).sum() + 1),
        "shape_fit": {"t0": _ym(c25["t0"]), "rmax": c25["shape_rmax"], "b_months": c25["shape_b_months"]},
    }

    # The early-rise "hindcast": how big would the cycle be, judged from the
    # smoothed number some months after minimum? Calibrated on earlier cycles,
    # in log space because the scatter grows with amplitude.
    cal = cycles[cycles["cycle"].isin(PREDICTION_CALIBRATION_CYCLES)]
    i25 = int(c25["i_min"])
    hind = []
    for lead in PREDICTION_LEAD_MONTHS:
        x = np.array([ours[int(i) + lead] for i in cal["i_min"]])
        y = cal["max_sn"].to_numpy(float)
        fit = sp.line_fit(np.log(x), np.log(y))
        x25 = ours[i25 + lead]
        yhat, lo68, hi68 = sp.predict(fit, np.log(x25), 0.68)
        _, lo95, hi95 = sp.predict(fit, np.log(x25), 0.95)
        hind.append({
            "lead_months": lead,
            # The smoothed value at min+lead needs 6 more months of data; the minimum
            # itself was only certain some months after it happened, too.
            "available": _ym(t[i25 + lead + 6]), "t_available": t[i25 + lead + 6],
            "s_at_lead": float(x25), "r": float(np.corrcoef(np.log(x), np.log(y))[0, 1]),
            "prediction": float(np.exp(yhat)), "lo68": float(np.exp(lo68)), "hi68": float(np.exp(hi68)),
            "lo95": float(np.exp(lo95)), "hi95": float(np.exp(hi95)),
            "observed_within_68": bool(np.exp(lo68) <= c25["max_sn"] <= np.exp(hi68)),
            "panel_within_68": bool(np.exp(lo68) <= lit["max_ssn"] <= np.exp(hi68)),
            "panel_within_95": bool(np.exp(lo95) <= lit["max_ssn"] <= np.exp(hi95)),
        })
    res["hindcast"] = hind
    res["hindcast_calibration"] = [int(c) for c in cal["cycle"]]
    first_excl = next((h for h in hind if not h["panel_within_68"]), None)
    res["panel_first_excluded_68"] = None if first_excl is None else {
        "lead_months": first_excl["lead_months"], "available": first_excl["available"]}

    # The average cycle (superposed on its minimum) for context, cycles 10-24.
    span = 150
    stack = np.array([ours[int(i):int(i) + span] for i in cal["i_min"]])
    avg = np.nanmean(stack, axis=0)
    sd = np.nanstd(stack, axis=0, ddof=1)
    k25 = min(span, len(ours) - i25)
    res["mean_max_calibration_cycles"] = float(cal["max_sn"].mean())
    res["median_max_calibration_cycles"] = float(cal["max_sn"].median())
    res["vs_mean_max"] = float(c25["max_sn"] / cal["max_sn"].mean())
    tables = {"panel": panel, "current_forecast": current,
              "average_cycle": pd.DataFrame({"months": np.arange(span), "mean": avg, "sd": sd}),
              "cycle25_smoothed": ours[i25:i25 + k25]}
    return tables, res


# --- 4. Two hemispheres ------------------------------------------------------------------------


def hemispheres(cycles: pd.DataFrame, area: pd.DataFrame) -> tuple[dict, dict]:
    """When each hemisphere peaked, cycle by cycle: from SILSO's hemispheric sunspot
    numbers (1992 on) and from the group areas (1874 on)."""
    h = silso.hemispheric()
    h["north_s"] = sp.smooth_13_month(h["sn_north"].to_numpy(float))
    h["south_s"] = sp.smooth_13_month(h["sn_south"].to_numpy(float))
    th = h["t"].to_numpy(float)
    silso_rows = []
    for k in range(len(cycles)):
        lo = cycles.at[k, "t_min"]
        hi = cycles.at[k + 1, "t_min"] if k + 1 < len(cycles) else np.inf
        sel = (th >= lo) & (th < hi) & np.isfinite(h["north_s"].to_numpy())
        if sel.sum() < 24 or lo < th[0]:
            continue
        iN = int(np.argmax(np.where(sel, h["north_s"], -np.inf)))
        iS = int(np.argmax(np.where(sel, h["south_s"], -np.inf)))
        silso_rows.append({"cycle": int(cycles.at[k, "cycle"]), "north_peak": _ym(th[iN]),
                           "north_sn": float(h["north_s"].iat[iN]), "south_peak": _ym(th[iS]),
                           "south_sn": float(h["south_s"].iat[iS]),
                           "south_minus_north_yr": float(th[iS] - th[iN])})

    mon = area.resample("MS").mean()
    tm = sp.decimal_year((mon.index + pd.Timedelta(days=14)).to_numpy())
    N = sp.smooth_13_month(mon["north"].to_numpy(float))
    S = sp.smooth_13_month(mon["south"].to_numpy(float))
    rows = []
    for k in range(len(cycles)):
        c = int(cycles.at[k, "cycle"])
        if c < FIRST_FULL_GROUP_CYCLE:
            continue
        lo = cycles.at[k, "t_min"]
        hi = cycles.at[k + 1, "t_min"] if k + 1 < len(cycles) else np.inf
        sel = (tm >= lo) & (tm < hi) & np.isfinite(N)
        iN = int(np.argmax(np.where(sel, N, -np.inf)))
        iS = int(np.argmax(np.where(sel, S, -np.inf)))
        tot_n = np.nansum(mon["north"].to_numpy()[(tm >= lo) & (tm < hi)])
        tot_s = np.nansum(mon["south"].to_numpy()[(tm >= lo) & (tm < hi)])
        rows.append({"cycle": c, "north_peak": tm[iN], "south_peak": tm[iS],
                     "north_peak_area": N[iN], "south_peak_area": S[iS],
                     "lag_south_minus_north_yr": tm[iS] - tm[iN],
                     "asymmetry": (tot_n - tot_s) / (tot_n + tot_s)})
    lags = pd.DataFrame(rows)
    res = {
        "silso": silso_rows,
        "areas": lags.to_dict(orient="records"),
        "mean_abs_lag_yr": float(lags["lag_south_minus_north_yr"].abs().mean()),
        "max_abs_lag_yr": float(lags["lag_south_minus_north_yr"].abs().max()),
        "n_south_later": int((lags["lag_south_minus_north_yr"] > 0).sum()),
        "n_cycles": len(lags),
        "n_north_dominant": int((lags["asymmetry"] > 0).sum()),
    }
    tables = {"hemispheric": h, "hemispheric_area_monthly": mon.assign(t=tm, north_s=N, south_s=S),
              "hemisphere_lags": lags}
    return tables, res


# --- 5. The butterfly diagram and the drift of the sunspot zones ------------------------------


def measure_area_seam(spots: pd.DataFrame, days, missing) -> tuple[dict, pd.DataFrame]:
    """How much smaller are the USAF/NOAA spot areas than RGO's?

    Use the sunspot number as a transfer standard: in monthly means, fit
    log(area) = a + b log(SN) + d [after 1976] across the seam. exp(-d) is the
    factor by which USAF/NOAA areas fall short *at the same sunspot number*.
    Months with SN <= 20 are left out (both measures are then small and noisy).
    Errors: bootstrap over whole years.
    """
    area = groups.daily_area(spots, days, missing, usaf_factor=1.0)
    sn = silso.daily().set_index("date")["sn"]
    j = area[["total"]].join(sn, how="inner").dropna()
    mon = j.resample("MS").mean().dropna()
    mon = mon[(mon["sn"] > 20) & (mon["total"] > 0)]
    seam = RGO_LAST_YEAR + 1
    yr = mon.index.year.to_numpy()

    def fit(frame, dummies):
        X = np.column_stack([np.ones(len(frame)), np.log(frame["sn"])] + dummies)
        coef, *_ = np.linalg.lstsq(X, np.log(frame["total"]), rcond=None)
        return coef

    win = mon[(yr >= seam - AREA_SEAM_WINDOW_YR) & (yr < seam + AREA_SEAM_WINDOW_YR)]
    wy = win.index.year.to_numpy()
    after = (wy >= seam).astype(float)
    coef = fit(win, [after])
    rng = np.random.default_rng(7)
    years = np.unique(wy)
    boot = []
    for _ in range(N_BOOT):
        pick = rng.choice(years, size=years.size, replace=True)
        idx = np.concatenate([np.flatnonzero(wy == y) for y in pick])
        sub = win.iloc[idx]
        sy = sub.index.year.to_numpy()
        if (sy >= seam).all() or (sy < seam).all():
            continue
        boot.append(np.exp(-fit(sub, [(sy >= seam).astype(float)])[2]))
    boot = np.array(boot)
    # USAF (1977-81) and NOAA (1982 on) separately.
    coef2 = fit(win, [((wy >= seam) & (wy < 1982)).astype(float), (wy >= 1982).astype(float)])

    # Area per sunspot number, decade by decade (common slope over the whole record),
    # relative to the RGO average: the seam, and how flat each side of it is. The
    # periods break exactly where the source changes: RGO to 1976, USAF 1977-81,
    # NOAA from 1982.
    start = np.where(yr < seam, (yr // 10) * 10,
                     np.where(yr < 1982, seam, np.maximum((yr // 10) * 10, 1982)))
    periods = np.unique(start)
    X = np.column_stack([np.log(mon["sn"])] + [(start == p).astype(float) for p in periods])
    c3, *_ = np.linalg.lstsq(X, np.log(mon["total"]), rcond=None)
    level = np.exp(c3[1:])
    rgo_mean = np.exp(np.mean(np.log(level[periods < seam])))
    first_year = [int(yr[start == p].min()) for p in periods]
    last_year = [int(yr[start == p].max()) for p in periods]
    per_decade = pd.DataFrame({"first_year": first_year, "last_year": last_year,
                               "source": ["RGO" if p < seam else "USAF" if p < 1982 else "NOAA"
                                          for p in periods],
                               "relative_area_per_sn": level / rgo_mean,
                               "n_months": [int((start == p).sum()) for p in periods]})
    res = {
        "window": [int(seam - AREA_SEAM_WINDOW_YR), int(seam + AREA_SEAM_WINDOW_YR - 1)],
        "n_months": len(win), "slope": float(coef[1]),
        "factor": float(np.exp(-coef[2])),
        "factor_ci68": [float(np.percentile(boot, 16)), float(np.percentile(boot, 84))],
        "factor_usaf_1977_1981": float(np.exp(-coef2[2])),
        "factor_noaa_1982_on": float(np.exp(-coef2[3])),
        "literature": LIT_AREA_FACTOR,
    }
    return res, per_decade


def _assign_cycles(t_months, alat, cycles: pd.DataFrame, lam0: float, tau: float,
                   sigma: float, early_months: float = CYCLE_EARLY_MONTHS,
                   floor: float = CYCLE_PRIOR_FLOOR) -> np.ndarray:
    """Which cycle does each sunspot group belong to?

    Near a minimum the old cycle's spots (near the equator) and the new cycle's (at
    high latitude) are on the disc together. Each group goes to the cycle that makes
    it most probable: that cycle's activity at the time (its fitted cycle shape)
    times a Gaussian in latitude about that cycle's drift-law latitude. A small
    floor lets a new cycle's first high-latitude spots appear up to three years
    before its fitted start. Returns cycle numbers.
    """
    t_months = np.asarray(t_months, dtype=float)
    alat = np.asarray(alat, dtype=float)
    logp = np.full((t_months.size, len(cycles)), -np.inf)
    for k, c in enumerate(cycles.itertuples()):
        t0 = _months(c.t0)
        dt = t_months - t0
        prior = sp.cycle_shape(t_months, t0, c.shape_rmax, c.shape_b_months)
        prior = prior + np.where((dt > -early_months) & (dt <= 0), floor, 0.0)
        ok = prior > 0
        lam = sp.drift_law(dt[ok], lam0, tau)
        logp[ok, k] = np.log(prior[ok]) - 0.5 * ((alat[ok] - lam) / sigma) ** 2
    best = np.argmax(logp, axis=1)
    out = cycles["cycle"].to_numpy()[best]
    out[~np.isfinite(logp.max(axis=1))] = 0
    return out


def _centroids(s: pd.DataFrame) -> pd.DataFrame:
    """Area-weighted mean |latitude| per (cycle, hemisphere, Carrington rotation)."""
    g = s.assign(w=s["area"].astype(float), wl=s["area"] * s["alat"]).groupby(["cycle", "hem", "rot"])
    cen = g[["w", "wl"]].sum().reset_index()
    cen["lat"] = cen["wl"] / cen["w"]
    cen = cen.rename(columns={"w": "area"}).drop(columns="wl")
    jd_mid = CARRINGTON_EPOCH_JD + (cen["rot"] - 0.5) * CARRINGTON_SYNODIC_DAYS
    cen["t"] = _decimal_year_of_jd(jd_mid.to_numpy())
    return cen


def butterfly(cycles: pd.DataFrame) -> tuple[dict, dict]:
    rec = groups.load()
    spots, days, missing = groups.split(rec)

    # --- Our parse against Hathaway's own derived files.
    listed = groups.listed_missing_days()
    area_h = groups.daily_area(spots, days, missing, usaf_factor=USAF_AREA_FACTOR_HATHAWAY)
    hd = groups.hathaway_daily_area()
    j = area_h[["total"]].join(hd[["total"]], rsuffix="_h", how="inner")
    both = j["total"].notna() & j["total_h"].notna()
    validation = {
        "days_compared": int(both.sum()),
        "days_matching": int((np.abs(j["total"] - j["total_h"])[both] < 0.5).sum()),
        "missing_days_ours": len(missing), "missing_days_listed": len(listed),
        "missing_days_agree": bool(set(missing) == set(listed)),
        "missing_flags_agree": int(((j["total"].isna()) == (j["total_h"].isna())).sum()),
    }

    seam, per_decade = measure_area_seam(spots, days, missing)

    # --- Census of the records.
    has_area = spots["area"] > 0
    eras = spots.groupby("era").agg(records=("area", "size"), groups=("track", "nunique"),
                                    first=("date", "min"), last=("date", "max"))
    w = spots.loc[has_area, "area"].astype(float) * np.where(
        spots.loc[has_area, "era"] == "RGO", 1.0, USAF_AREA_FACTOR_HATHAWAY)
    alat = spots.loc[has_area, "lat"].abs()
    order = np.argsort(alat.to_numpy())
    cum = np.cumsum(w.to_numpy()[order]) / w.sum()
    census = {
        "records": len(spots), "records_with_area": int(has_area.sum()),
        "groups": int(spots["track"].nunique()),
        "first_day": days[0].strftime("%Y-%m-%d"), "last_day": days[-1].strftime("%Y-%m-%d"),
        "days": len(days), "days_observed": len(days) - len(missing),
        "spotless_days": int(len(days) - len(missing) - spots["date"].nunique()),
        "eras": {e: {"records": int(r["records"]), "groups": int(r["groups"]),
                     "first": r["first"].strftime("%Y-%m-%d"), "last": r["last"].strftime("%Y-%m-%d")}
                 for e, r in eras.iterrows()},
        "area_fraction_within_5deg": float(w[alat < 5].sum() / w.sum()),
        "area_fraction_beyond_30deg": float(w[alat > 30].sum() / w.sum()),
        "area_fraction_beyond_40deg": float(w[alat > 40].sum() / w.sum()),
        "lat_containing_99pct_area": float(alat.to_numpy()[order][np.searchsorted(cum, 0.99)]),
        "median_abs_lat_area_weighted": float(alat.to_numpy()[order][np.searchsorted(cum, 0.5)]),
        "records_beyond_50deg": int((spots["lat"].abs() > 50).sum()),
        "highest_latitude": spots.loc[spots["lat"].abs().idxmax(), ["date", "lat", "area", "era"]]
        .astype(str).to_dict(),
        "usaf_latitudes_whole_degrees": float(np.mean(
            np.abs(spots.loc[spots["year"] >= 1982, "lat"] % 1) < 1e-9)),
    }

    # --- Butterfly image: mean daily spot area per 1-degree strip per rotation.
    s = spots[has_area & (spots["lat"] != 0)].copy()
    s["area_adj"] = s["area"] * np.where(s["era"] == "RGO", 1.0, USAF_AREA_FACTOR_HATHAWAY)
    s["alat"] = s["lat"].abs()
    s["hem"] = np.sign(s["lat"]).astype(int)
    s["rot"] = np.floor(sp.carrington_rotation_number(s["jd"])).astype(int)
    obs_days = pd.Series(1.0, index=days.difference(missing))
    day_rot = np.floor(sp.carrington_rotation_number(
        sp.julian_date(obs_days.index.to_numpy()) + 0.5)).astype(int)
    n_obs_per_rot = pd.Series(1, index=day_rot).groupby(level=0).sum()
    lat_edges = np.arange(-BUTTERFLY_LAT_MAX_DEG, BUTTERFLY_LAT_MAX_DEG + 1e-9, BUTTERFLY_LAT_BIN_DEG)
    rots = np.arange(s["rot"].min(), s["rot"].max() + 1)
    H, _, _ = np.histogram2d(s["rot"], s["lat"], bins=[np.r_[rots, rots[-1] + 1] - 0.5, lat_edges],
                             weights=s["area_adj"])
    nobs = n_obs_per_rot.reindex(rots).fillna(0).to_numpy()
    with np.errstate(invalid="ignore", divide="ignore"):
        grid = H / nobs[:, None]
    grid[nobs == 0] = np.nan
    rot_t_edges = _decimal_year_of_jd(CARRINGTON_EPOCH_JD + (np.r_[rots, rots[-1] + 1] - 1)
                                      * CARRINGTON_SYNODIC_DAYS)

    # --- Assign every group to its cycle, then measure the drift of the zones.
    s["t_months"] = _months(s["t"])
    lam0, tau, sigma = LIT_HATHAWAY_DRIFT["lambda0_deg"], LIT_HATHAWAY_DRIFT["tau_months"], 5.0
    history = []
    for _ in range(6):
        s["cycle"] = _assign_cycles(s["t_months"], s["alat"], cycles, lam0, tau, sigma)
        cen = _centroids(s)
        cen = cen[(cen["area"] >= CENTROID_MIN_AREA) & (cen["cycle"] >= FIRST_FULL_GROUP_CYCLE)]
        t0 = cycles.set_index("cycle")["t0"]
        cen["dt_months"] = (cen["t"] - cen["cycle"].map(t0)) * 12.0
        fit = sp.fit_drift_law(cen["dt_months"], cen["lat"])
        new = (fit["lambda0_deg"], fit["tau_months"])
        sel = s["cycle"] >= FIRST_FULL_GROUP_CYCLE
        resid = s.loc[sel, "alat"] - sp.drift_law(
            (s.loc[sel, "t"] - s.loc[sel, "cycle"].map(t0)) * 12.0, *new)
        sigma = float(np.sqrt(np.average(resid**2, weights=s.loc[sel, "area_adj"])))
        history.append([*new, sigma])
        if abs(new[0] - lam0) < 1e-3 and abs(new[1] - tau) < 1e-2:
            break
        lam0, tau = new
    s["cycle"] = _assign_cycles(s["t_months"], s["alat"], cycles, lam0, tau, sigma)
    cen = _centroids(s)
    cen = cen[(cen["area"] >= CENTROID_MIN_AREA) & (cen["cycle"] >= FIRST_FULL_GROUP_CYCLE)].copy()
    cyc = cycles.set_index("cycle")
    cen["dt_months"] = (cen["t"] - cen["cycle"].map(cyc["t0"])) * 12.0
    cen["dt_min_months"] = (cen["t"] - cen["cycle"].map(cyc["t_min"])) * 12.0
    fit_t0 = sp.fit_drift_law(cen["dt_months"], cen["lat"])
    fit_min = sp.fit_drift_law(cen["dt_min_months"], cen["lat"])

    # Bootstrap over whole cycles: the honest error bar when cycles differ.
    rng = np.random.default_rng(11)
    ids = np.sort(cen["cycle"].unique())
    by = {c: cen[cen["cycle"] == c] for c in ids}
    boot = []
    for _ in range(N_BOOT // 4):
        pick = rng.choice(ids, size=ids.size, replace=True)
        sub = pd.concat([by[c] for c in pick])
        f = sp.fit_drift_law(sub["dt_months"], sub["lat"])
        boot.append((f["lambda0_deg"], f["tau_months"]))
    boot = np.array(boot)
    fit_t0["lambda0_ci68"] = np.percentile(boot[:, 0], [16, 84]).tolist()
    fit_t0["tau_ci68"] = np.percentile(boot[:, 1], [16, 84]).tolist()

    # Hathaway's own start times for cycles 12-23, same data otherwise.
    h_t0 = {c: _t_of_ym(v) for c, v in LIT_HATHAWAY_T0.items()}
    ch = cen[cen["cycle"].isin(list(h_t0))]
    fit_hathaway_t0 = sp.fit_drift_law((ch["t"] - ch["cycle"].map(h_t0)) * 12.0, ch["lat"])
    fit_our_t0_same = sp.fit_drift_law(ch["dt_months"], ch["lat"])

    # Does the path depend on cycle strength? Per-cycle latitude scale at fixed tau,
    # against amplitude, with the cycle clock started at t0 or at the minimum.
    strength = {}
    for clock, col, f in (("t0", "dt_months", fit_t0), ("minimum", "dt_min_months", fit_min)):
        per = []
        for c in ids:
            d = cen[cen["cycle"] == c]
            e = np.exp(-d[col] / f["tau_months"])
            per.append((c, float(np.sum(d["lat"] * e) / np.sum(e * e))))
        per = pd.DataFrame(per, columns=["cycle", "lambda0"])
        per["amplitude"] = per["cycle"].map(cyc["max_sn"])
        strength[clock] = {"per_cycle": per.to_dict(orient="records"),
                           "correlation": sp.correlation(per["amplitude"], per["lambda0"]),
                           "spread_deg": float(per["lambda0"].std(ddof=1))}

    # A model-free version of the same question: the area-weighted mean latitude of
    # each cycle's spots in its first three years.
    early = []
    for c in ids:
        d = s[(s["cycle"] == c) & (s["t"] >= cyc.at[c, "t0"]) & (s["t"] < cyc.at[c, "t0"] + 3)]
        early.append((int(c), float(np.average(d["alat"], weights=d["area_adj"])), cyc.at[c, "max_sn"]))
    early = pd.DataFrame(early, columns=["cycle", "mean_lat_first_3yr", "amplitude"])
    strength["first_3_years"] = {"per_cycle": early.to_dict(orient="records"),
                                 "correlation": sp.correlation(early["amplitude"],
                                                               early["mean_lat_first_3yr"]),
                                 "fit": sp.line_fit(early["amplitude"], early["mean_lat_first_3yr"])}
    # ... and with Hathaway's own start times (cycles 12-23).
    per = []
    for c in sorted(h_t0):
        d = cen[cen["cycle"] == c]
        e = np.exp(-(d["t"] - h_t0[c]) * 12.0 / fit_t0["tau_months"])
        per.append((c, float(np.sum(d["lat"] * e) / np.sum(e * e)), cyc.at[c, "max_sn"]))
    per = pd.DataFrame(per, columns=["cycle", "lambda0", "amplitude"])
    strength["hathaway_t0"] = {"correlation": sp.correlation(per["amplitude"], per["lambda0"])}

    # Leaving out the last years of each cycle, when the next one overlaps it.
    late = cen[cen["dt_months"] < 108]
    fit_early_part = sp.fit_drift_law(late["dt_months"], late["lat"])

    # Around each minimum the two cycles share the disc: when the new cycle's first
    # groups appear, when the old cycle's last ones fade, and where each sits.
    first = s.groupby("track").agg(cycle=("cycle", "first"), t_first=("t", "min"),
                                   t_last=("t", "max"), lat=("alat", "mean"))
    overlap = []
    for c in range(FIRST_FULL_GROUP_CYCLE + 1, int(cycles["cycle"].max()) + 1):
        tmin = cyc.at[c, "t_min"]
        new = np.sort(first.loc[first["cycle"] == c, "t_first"].to_numpy())
        old = np.sort(first.loc[first["cycle"] == c - 1, "t_last"].to_numpy())
        near = first[(first["t_first"] > tmin - 1) & (first["t_first"] < tmin + 1)]
        overlap.append({
            "minimum_of_cycle": c, "minimum": cyc.at[c, "min_ym"],
            "first_new_months_before_min": (tmin - new[0]) * 12.0,
            "third_new_months_before_min": (tmin - new[2]) * 12.0,
            "last_old_months_after_min": (old[-1] - tmin) * 12.0,
            "third_last_old_months_after_min": (old[-3] - tmin) * 12.0,
            "groups_within_1yr_new": int((near["cycle"] == c).sum()),
            "groups_within_1yr_old": int((near["cycle"] == c - 1).sum()),
            "median_lat_new": float(near.loc[near["cycle"] == c, "lat"].median()),
            "median_lat_old": float(near.loc[near["cycle"] == c - 1, "lat"].median()),
        })
    ov = pd.DataFrame(overlap)

    res = {
        "validation": validation, "seam": seam, "census": census,
        "drift_law": {
            "fit": fit_t0, "fit_minimum_clock": fit_min, "fit_first_108_months": fit_early_part,
            "iterations": history, "group_sigma_deg": sigma,
            "cycles": [int(c) for c in ids], "n_centroids": int(len(cen)),
            "hathaway_t0_cycles_12_23": fit_hathaway_t0, "our_t0_cycles_12_23": fit_our_t0_same,
            "literature": LIT_HATHAWAY_DRIFT,
            "strength_dependence": strength,
        },
        "overlap": {"per_minimum": ov.to_dict(orient="records"),
                    "median_first_new_months_before_min": float(ov["first_new_months_before_min"].median()),
                    "range_first_new_months_before_min": [float(ov["first_new_months_before_min"].min()),
                                                          float(ov["first_new_months_before_min"].max())],
                    "median_last_old_months_after_min": float(ov["last_old_months_after_min"].median()),
                    "range_last_old_months_after_min": [float(ov["last_old_months_after_min"].min()),
                                                        float(ov["last_old_months_after_min"].max())],
                    "median_third_new_months_before_min": float(ov["third_new_months_before_min"].median()),
                    "median_third_last_old_months_after_min": float(ov["third_last_old_months_after_min"].median()),
                    "median_lat_new_near_min": float(ov["median_lat_new"].median()),
                    "median_lat_old_near_min": float(ov["median_lat_old"].median())},
    }
    tables = {"spots": s, "all_spots": spots, "days": days, "missing": missing,
              "daily_area": area_h, "centroids": cen, "area_per_decade": per_decade,
              "butterfly_grid": grid, "butterfly_t_edges": rot_t_edges, "butterfly_lat_edges": lat_edges}
    return tables, res


# --- 6. Differential rotation ---------------------------------------------------------------------


def rotation(spots: pd.DataFrame) -> tuple[dict, dict]:
    """Track each RGO group from one day to the next: its drift in Carrington
    longitude gives its sidereal rotation rate."""
    r = spots[(spots["era"] == "RGO") & spots["year"].between(*ROTATION_YEARS)]
    r = r[~r.duplicated(["track", "date"], keep=False)]  # ambiguous same-day sub-groups

    def pairs_for(max_cmd: float) -> pd.DataFrame:
        q = r[r["cmd"].abs() <= max_cmd].sort_values(["track", "jd"])
        g = q.groupby("track", sort=False)
        dt = g["jd"].diff()
        dlon = (g["carrington_lon"].diff() + 180.0) % 360.0 - 180.0
        dlat = g["lat"].diff()
        p = pd.DataFrame({"track": q["track"], "dt": dt, "dlon": dlon,
                          "lat": q["lat"] - dlat / 2, "dlat": dlat, "t": q["t"]}).dropna()
        p = p[p["dt"].between(*ROTATION_PAIR_DT_DAYS)]
        p["drift"] = p["dlon"] / p["dt"]
        n_all = len(p)
        p = p[p["drift"].abs() <= ROTATION_MAX_DRIFT_DEG_PER_DAY].copy()
        p.attrs["rejected"] = n_all - len(p)
        p["omega"] = sp.sidereal_rate(p["drift"])
        return p

    p = pairs_for(ROTATION_MAX_CMD_DEG)
    fit = sp.fit_rotation_law(p["lat"], p["omega"])

    # Errors: bootstrap over groups (pairs from one group are not independent).
    x = np.sin(np.radians(p["lat"].to_numpy())) ** 2
    y = p["omega"].to_numpy()
    codes, uniq = pd.factorize(p["track"])
    sums = np.zeros((uniq.size, 5))
    for col, v in enumerate((np.ones_like(x), x, x * x, y, x * y)):
        sums[:, col] = np.bincount(codes, weights=v, minlength=uniq.size)
    rng = np.random.default_rng(5)
    boot = []
    for _ in range(N_BOOT // 200):
        wts = rng.poisson(1.0, size=(200, uniq.size)).astype(float)
        S = wts @ sums  # n, Sx, Sxx, Sy, Sxy per resample
        det = S[:, 0] * S[:, 2] - S[:, 1] ** 2
        B = (S[:, 0] * S[:, 4] - S[:, 1] * S[:, 3]) / det
        A = (S[:, 3] - B * S[:, 1]) / S[:, 0]
        boot.append(np.column_stack([A, B]))
    boot = np.vstack(boot)
    fit["A_err_boot"] = float(np.std(boot[:, 0], ddof=1))
    fit["B_err_boot"] = float(np.std(boot[:, 1], ddof=1))

    edges = np.arange(0, 40 + ROTATION_LAT_BIN_DEG, ROTATION_LAT_BIN_DEG)
    b = pd.cut(p["lat"].abs(), edges)
    binned = p.groupby(b, observed=True)["omega"].agg(["size", "mean", "median", "std"]).reset_index()
    binned["lat_mid"] = [iv.mid for iv in binned["lat"]]
    binned["sem"] = binned["std"] / np.sqrt(binned["size"])
    binned = binned.drop(columns="lat")

    robustness = {}
    for cmd in (30.0, 45.0, 75.0):
        q = pairs_for(cmd)
        f = sp.fit_rotation_law(q["lat"], q["omega"])
        robustness[f"cmd_{int(cmd)}"] = {"A": f["A"], "B": f["B"], "n": len(q)}
    med = binned[binned["size"] >= 100]
    fm = sp.fit_rotation_law(med["lat_mid"], med["median"])
    robustness["binned_medians"] = {"A": fm["A"], "B": fm["B"]}

    periods = {}
    for lat in (0.0, 15.0, 30.0, 40.0):
        om = float(sp.rotation_law(lat, fit["A"], fit["B"]))
        periods[f"{int(lat)}deg"] = {"omega": om, "sidereal_days": float(sp.period_days(om)),
                                     "synodic_days": float(sp.synodic_from_sidereal(sp.period_days(om)))}
    lit = LIT_BALTHASAR_1986
    res = {
        "n_pairs": len(p), "n_groups": int(uniq.size), "rejected_pairs": int(p.attrs["rejected"]),
        "years": list(ROTATION_YEARS), "max_cmd_deg": ROTATION_MAX_CMD_DEG,
        "fit": fit, "literature": lit,
        "A_diff": fit["A"] - lit["A"], "B_diff": fit["B"] - lit["B"],
        "periods": periods, "robustness": robustness,
        "equator_to_40deg_lap_days": float(360.0 / (fit["A"] - float(sp.rotation_law(40.0, fit["A"], fit["B"])))),
    }
    return {"rotation_pairs": p, "rotation_binned": binned}, res


# --- 7. Before 1749 --------------------------------------------------------------------------------


def long_record(yearly: pd.DataFrame) -> tuple[dict, dict]:
    """Put each group-number reconstruction on the version-2 sunspot-number scale
    (a single factor, fitted where both overlap within 1874-1976), then ask what it
    says about the 18th century and about the Maunder minimum."""
    y = yearly.set_index("t")["sn"]
    gns = silso.group_numbers()
    out, scaled = {}, {}
    full_years = y[(y.index >= 1749) & (y.index < 2026)]
    for key, d in gns.items():
        d = d.set_index("t")
        j = d[["gn"]].join(y, how="inner").dropna()
        cal = j[(j.index >= 1874) & (j.index <= 1976)]
        k = float((cal["gn"] * cal["sn"]).sum() / (cal["gn"] ** 2).sum())
        sc = d.assign(sn_equiv=k * d["gn"])
        scaled[key] = sc.reset_index()
        early = j[(j.index >= 1749) & (j.index < 1850)]
        info = {"label": GROUP_NUMBER_SERIES[key]["label"], "family": GROUP_NUMBER_SERIES[key]["family"],
                "first_year": float(d.index.min()), "last_year": float(d.index.max()),
                "scale_to_sn_v2": k, "calibration_years": [float(cal.index.min()), float(cal.index.max())],
                "calibration_n": len(cal)}
        if len(early) >= 30:
            info["early_ratio_1749_1849"] = float((k * early["gn"]).mean() / early["sn"].mean())
        mm = sc[(sc.index >= MAUNDER_MINIMUM[0]) & (sc.index < MAUNDER_MINIMUM[1] + 1)]
        if len(mm):
            info["maunder_mean_sn_equiv"] = float(mm["sn_equiv"].mean())
            info["maunder_fraction_of_modern_mean"] = float(mm["sn_equiv"].mean() / full_years.mean())
            info["maunder_years"] = len(mm)
            if mm["n_obs"].notna().any():
                # Only Hoyt & Schatten list how many observations stand behind each
                # year -- the evidence that a blank year was watched, not missed.
                info["maunder_years_with_no_spot_seen"] = int((mm["gn"] == 0).sum())
                info["maunder_median_obs_per_year"] = float(mm["n_obs"].median())
                blank = mm[mm["gn"] == 0]
                info["maunder_median_obs_in_blank_years"] = float(blank["n_obs"].median())
        out[key] = info
    res = {"series": out, "modern_mean_sn_1749_2025": float(full_years.mean()),
           "hs98_scale_expected": LIT_HS98_SCALE_V2,
           "sn_1700_1748_mean": float(y[(y.index >= 1700) & (y.index < 1749)].mean())}
    return {"group_numbers": scaled}, res


# --- Everything ----------------------------------------------------------------------------------------


def analyse() -> tuple[dict, dict]:
    tables, results = {}, {}
    t, r = record()
    tables.update(t)
    results["record"] = r
    results["rules"] = rules(tables["cycles"])
    t, r = cycle25(tables["monthly"], tables["cycles"])
    tables.update(t)
    results["cycle25"] = r
    t, r = butterfly(tables["cycles"])
    tables.update(t)
    results["butterfly"] = r
    t, r = hemispheres(tables["cycles"], tables["daily_area"])
    tables.update(t)
    results["hemispheres"] = r
    t, r = rotation(tables["all_spots"])
    tables.update(t)
    results["rotation"] = r
    t, r = long_record(tables["yearly"])
    tables.update(t)
    results["long_record"] = r
    results["cycles"] = tables["cycles"].drop(columns=["i_min", "i_max"]).to_dict(orient="records")
    return tables, _clean(results)


def save_results(results: dict) -> None:
    RESULTS_PATH.write_text(json.dumps(results, indent=2) + "\n")
