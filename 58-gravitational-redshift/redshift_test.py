"""The gravitational-redshift test with the eccentric Galileo clocks, end to end.

    python redshift_test.py          # ~1 h on first run (null tests), cached afterwards

Steps (each cached under data/cache/):
1. per-day fits for E14, E18 (alpha, alpha_tot) and all satellites (J2 and tides);
2. which clock-days to use (E14's maser until it degraded; E18 after its 2016-17
   non-nominal maser period);
3. global fits of alpha with Sun-locked nuisance signals, block-bootstrap errors;
4. null tests: the same fit with "fake perigee" templates, which feel the same
   systematics but carry no gravity signal; their scatter calibrates the error;
5. robustness splits; laser-ranging (SLR) orbit check with the same nuisance model;
6. the error budget and the combined result (Delva et al.'s alpha, and Herrmann
   et al.'s alpha_tot ~ alpha/2).
"""

from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd

import clockfit as F
import config as C
import globalfit as G
import prepare
import slr

CACHE = C.CACHE_DIR

#: Clock-days used for the test. E14 (GSAT-0202) ran a passive hydrogen maser that
#: was excellent until its noise began to grow in early 2024 (and the operators
#: re-tuned it repeatedly from 2025); E18 (GSAT-0201) ran the maser that ESA
#: called non-nominal in 2016-17 (Herrmann et al. 2018): its first months here are
#: ten times noisier, and its frequency was re-set in March 2017.
SEGMENTS = {
    "E14": ("2016-08-01", "2023-12-31"),
    "E18": ("2017-04-01", "2026-10-06"),
}
#: The span the published tests covered (both ended on 2017-12-16).
PUBLISHED_SPAN = ("2016-08-01", "2017-12-16")

PRIMARY_MODEL = "sun12_1yr_beta2"
MODELS = ("none", "sun1_1yr", "sun1_1yr_beta2", "sun12_1yr_beta2", "sun12_6mo_beta2")
NULL_DRIFTS_DEG_PER_YR = (-90, -60, -40, -25, 25, 40, 60, 90)
NULL_PHASES_DEG = (0, 90)

#: Systematics we cannot measure without onboard telemetry, adopted from Delva et
#: al. (2018, Table II): GSAT0202 = E14, GSAT0201 = E18. Both are upper-limit-style
#: estimates from ground characterisation of the masers (magnetic field along the
#: orbit from IGRF; temperature from the satellite's attitude to the Sun).
ADOPTED = {"E14": {"temperature": 0.13e-5, "magnetic": 1.92e-5},
           "E18": {"temperature": 0.59e-5, "magnetic": 1.93e-5}}

PRNS_ALL = ("E14", "E18") + C.CONTROLS


# ---------------------------------------------------------------------------
# 1. Per-day fits
# ---------------------------------------------------------------------------

def daily(prn: str, fit_set: str, force: bool = False) -> pd.DataFrame:
    path = CACHE / f"daily_{fit_set}_{prn}.parquet"
    if path.exists() and not force:
        return pd.read_parquet(path)
    df, _ = F.run_fits(prn, fit_set)
    df.to_parquet(path)
    return df


def segment_days(prn: str, df: pd.DataFrame, span: tuple | None = None) -> pd.DataFrame:
    lo, hi = span or SEGMENTS[prn]
    return df[df["ok"] & (df["date"] >= lo) & (df["date"] <= hi)]


# ---------------------------------------------------------------------------
# 2-5. Global fits, null tests, robustness, SLR
# ---------------------------------------------------------------------------

def _gfit(prn, days, model, d, **kw):
    r = G.fit(prn, days, G.NUISANCE_MODELS[model], d=d, **kw)
    return {k: v for k, v in r.items() if k not in ("coef", "names")}


def null_tests(prn: str, days: pd.DataFrame, d: dict, model: str = PRIMARY_MODEL) -> dict:
    idx = {s: k for k, s in enumerate(d["date"])}
    rows = np.array([idx[s] for s in days["date"]])
    om = np.unwrap(np.array([G.perigee_angle(d, i) for i in rows]))
    tday = pd.to_datetime(days["date"]).to_numpy().astype("datetime64[D]").astype(np.int64)
    amp = np.sqrt(C.GM_EARTH * np.nanmean(d["a"][rows])) * np.nanmean(d["e"][rows]) / C.C_LIGHT**2
    om_of, t_of = dict(zip(rows, om)), dict(zip(rows, tday - tday[0]))
    vals = []
    for drift in NULL_DRIFTS_DEG_PER_YR:
        for phase in NULL_PHASES_DEG:
            fn = (lambda i, dr=drift, ph=phase:
                  om_of[i] + np.radians(dr) * t_of[i] / 365.25 + np.radians(ph))
            r = G.fit(prn, days, G.NUISANCE_MODELS[model], d=d, n_boot=50,
                      template_fn=lambda dd, i, f=fn: G.fake_template(dd, i, f, amp))
            vals.append({"drift_deg_yr": drift, "phase_deg": phase, "alpha": r["alpha"]})
    v = np.array([x["alpha"] for x in vals])
    return {"model": model, "n": len(v), "rms": float(np.sqrt(np.mean(v**2))), "sd": float(np.std(v)),
            "mean": float(np.mean(v)), "max_abs": float(np.max(np.abs(v))), "values": vals,
            "perigee_rate_deg_yr": float(np.degrees(np.polyfit(tday - tday[0], om, 1)[0]) * 365.25)}


def slr_orbit(prn: str, passes: pd.DataFrame, span: tuple) -> dict:
    lo, hi = span
    p = passes[(passes["prn"] == prn) & (passes["start"] >= lo) & (passes["start"] <= hi + " 23:59")]
    p = slr.clean_passes(slr.pass_features(prn, p))
    model = G.NUISANCE_MODELS[PRIMARY_MODEL]
    kw = dict(origin=lo, window_days=model.window_days, harmonics=model.harmonics,
              beta_degree=model.beta_degree)
    r = slr.orbit_bias_nuisance(p, **kw)
    r["err_boot"] = slr.orbit_bias_bootstrap(p, n_boot=200, **kw)
    plain = slr.orbit_bias_nuisance(p, origin=lo, window_days=None)
    r["alpha_orbit_no_nuisance"] = plain["alpha_orbit"]
    # station-by-station (station offsets only): SLR systematics show up as disagreement
    st = []
    for station, g in p.groupby("station"):
        if len(g) < 150:
            continue
        rs = slr.orbit_bias_nuisance(g, origin=lo, window_days=None, station_offsets=False)
        st.append({"station": int(station), "n": int(len(g)), "alpha_orbit": rs["alpha_orbit"],
                   "err": rs["err_scaled"]})
    st = pd.DataFrame(st)
    w = 1 / st["err"] ** 2
    mean = float(np.sum(w * st["alpha_orbit"]) / np.sum(w))
    chi2 = float(np.sum(w * (st["alpha_orbit"] - mean) ** 2) / (len(st) - 1))
    r.update(stations=st.to_dict("records"), station_mean=mean, station_chi2_nu=chi2,
             station_err_scaled=float(np.sqrt(chi2 / np.sum(w))))
    r["n_passes"] = int(len(p))
    return r


def yearly(prn: str, df: pd.DataFrame, passes: pd.DataFrame) -> list[dict]:
    """Per calendar year: the simple daily-fit combination (no nuisance) and SLR's orbit term."""
    out = []
    ok = df[df["ok"]].assign(year=pd.to_datetime(df["date"]).dt.year)
    for year, g in ok.groupby("year"):
        r = F.combine(g, "alpha", use_gls=True)
        row = {"year": int(year), "alpha": r["value"], "err_boot": r["err_boot"], "n_days": r["n_days"],
               "rms_ps_median": float(g["rms_ps"].median()), "q_median": float(g["q"].median())}
        p = passes[(passes["prn"] == prn) & (pd.to_datetime(passes["start"]).dt.year == year)]
        if len(p) > 50:
            try:
                pf = slr.clean_passes(slr.pass_features(prn, p))
                rs = slr.orbit_bias_nuisance(pf, origin=f"{year}-01-01", window_days=None)
                row.update(alpha_orbit=rs["alpha_orbit"], alpha_orbit_err=rs["err_scaled"],
                           slr_passes=rs["n_passes"])
            except Exception:
                pass
        out.append(row)
    return out


# ---------------------------------------------------------------------------
# 6. Budget and combination
# ---------------------------------------------------------------------------

def budget(prn: str, res: dict) -> dict:
    prim = res["global"][PRIMARY_MODEL]
    stat = max(prim["err_boot"], res["null"]["rms"])     # null-calibrated
    vals = [res["global"][m]["alpha"] for m in MODELS if m != "none"]
    model_spread = 0.5 * (max(vals) - min(vals))
    orbit = float(np.hypot(res["slr"]["alpha_orbit"], res["slr"]["station_err_scaled"]))
    t, m = ADOPTED[prn]["temperature"], ADOPTED[prn]["magnetic"]
    total = float(np.sqrt(stat**2 + model_spread**2 + orbit**2 + t**2 + m**2))
    return {"alpha": prim["alpha"], "stat_null_calibrated": stat, "nuisance_model": model_spread,
            "orbit_slr": orbit, "temperature_adopted": t, "magnetic_adopted": m, "total": total}


def combine_satellites(b14: dict, b18: dict) -> dict:
    """Weighted by the clock-noise (statistical) errors; systematics combined as in Delva et al.:
    orbit terms weighted alike (partly independent), adopted terms taken as common."""
    w = np.array([1 / b14["stat_null_calibrated"] ** 2, 1 / b18["stat_null_calibrated"] ** 2])
    w /= w.sum()
    a = float(w @ [b14["alpha"], b18["alpha"]])
    stat = float(np.sqrt(np.sum((w * [b14["stat_null_calibrated"], b18["stat_null_calibrated"]]) ** 2)))
    out = {"alpha": a, "weights": w.tolist(), "stat": stat}
    for key in ("nuisance_model", "orbit_slr", "temperature_adopted", "magnetic_adopted"):
        out[key] = float(w @ [b14[key], b18[key]])
    out["total"] = float(np.sqrt(sum(out[k] ** 2 for k in
                                     ("stat", "nuisance_model", "orbit_slr", "temperature_adopted",
                                      "magnetic_adopted"))))
    return out


def run(force: bool = False) -> dict:
    path = CACHE / "redshift_results.json"
    if path.exists() and not force:
        return json.loads(path.read_text())
    t0 = time.time()
    passes = slr.read_summaries()
    passes.to_parquet(CACHE / "slr_passes.parquet")
    out = {"segments": SEGMENTS, "primary_model": PRIMARY_MODEL, "per_satellite": {}}
    # J2 and tides for every satellite (does each clock feel the Earth's flattening?)
    j2 = []
    for prn in PRNS_ALL:
        row = {"prn": prn}
        # (a) everything beyond the Keplerian convention, as one template (GR = 1)
        df = daily(prn, "beyond", force=force)
        sel = df[df["ok"] & (df["q"] < 1.5e-24)]
        if prn in SEGMENTS:
            sel = sel[(sel["date"] >= SEGMENTS[prn][0]) & (sel["date"] <= SEGMENTS[prn][1])]
        r = F.combine(sel, "beyond", use_gls=True)
        row.update(days=int(len(sel)), beyond=r["value"], beyond_err=r["err_boot"])
        # (b) the J2 potential and the tides on their own (orbit bending removed with the model)
        df = daily(prn, "j2_tide", force=force)
        sel = df[df["ok"] & (df["q"] < 1.5e-24)]
        if prn in SEGMENTS:
            sel = segment_days(prn, df)
        for col in ("j2", "tide"):
            r = F.combine(sel, col, use_gls=True)
            row[col] = r["value"]
            row[col + "_err"] = r["err_boot"]
        j2.append(row)
        print(f"  {prn}: beyond-Kepler {row['beyond']:.3f} +- {row['beyond_err']:.3f}; J2 potential "
              f"{row['j2']:.3f} +- {row['j2_err']:.3f}, tides {row['tide']:.2f} +- {row['tide_err']:.2f}")
    out["j2_tide"] = j2
    for prn in ("E14", "E18"):
        res = {}
        df = daily(prn, "alpha", force=force)
        daily(prn, "alpha_tot", force=force)
        d = prepare.load(prn)
        days = segment_days(prn, df)
        res["n_days"] = int(len(days))
        res["global"] = {m: _gfit(prn, days, m, d, n_boot=500) for m in MODELS}
        for m, r in res["global"].items():
            print(f"  {prn} {m:16s} alpha {r['alpha'] * 1e5:+6.2f} +- {r['err_boot'] * 1e5:.2f}")
        tot = daily(prn, "alpha_tot")
        res["global_alpha_tot"] = _gfit(prn, segment_days(prn, tot), PRIMARY_MODEL, d, n_boot=500,
                                        template="alpha_tot")
        res["published_span"] = {
            "none": _gfit(prn, segment_days(prn, df, PUBLISHED_SPAN), "none", d, n_boot=300),
            PRIMARY_MODEL: _gfit(prn, segment_days(prn, df, PUBLISHED_SPAN), PRIMARY_MODEL, d, n_boot=300),
        } if prn == "E14" else {}
        res["null"] = null_tests(prn, days, d)
        print(f"  {prn} null tests: rms {res['null']['rms'] * 1e5:.2f}")
        if prn == "E14":   # why the simpler Sun-locked model is not adopted
            res["null_sun1_1yr"] = null_tests(prn, days, d, model="sun1_1yr")
            print(f"  {prn} null tests, sun1_1yr: rms {res['null_sun1_1yr']['rms'] * 1e5:.2f}")
        splits = {"no eclipse days": days[days["shadow_epochs"] == 0],
                  "|beta| > 15 deg": days[days["beta_deg"].abs() > 15],
                  "first half": days.iloc[: len(days) // 2], "second half": days.iloc[len(days) // 2:]}
        res["splits"] = {k: _gfit(prn, v, PRIMARY_MODEL, d, n_boot=300) for k, v in splits.items()}
        res["slr"] = slr_orbit(prn, passes, SEGMENTS[prn])
        res["slr"].pop("passes", None)
        res["yearly"] = yearly(prn, df, passes)
        res["budget"] = budget(prn, res)
        out["per_satellite"][prn] = res
        print(f"  {prn}: alpha = {res['budget']['alpha'] * 1e5:+.2f} +- {res['budget']['total'] * 1e5:.2f} (x1e-5)")
    out["combined"] = combine_satellites(out["per_satellite"]["E14"]["budget"],
                                         out["per_satellite"]["E18"]["budget"])
    out["runtime_s"] = time.time() - t0
    path.write_text(json.dumps(out, indent=1, default=float))
    return out


if __name__ == "__main__":
    import sys
    r = run(force="--force" in sys.argv)
    c = r["combined"]
    print(f"combined alpha = ({c['alpha'] * 1e5:+.2f} +- {c['total'] * 1e5:.2f}) x 1e-5")


# ---------------------------------------------------------------------------
# Diagnostics for the figures (cached)
# ---------------------------------------------------------------------------

def once_per_rev(prn: str, force: bool = False) -> pd.DataFrame:
    """Per day: amplitude and phase of the once-per-orbit residual signal, measured from
    the Sun's direction and from the perigee. A redshift violation is perigee-locked;
    radiation-pressure orbit errors and thermal effects are Sun-locked."""
    path = CACHE / f"onceperrev_{prn}.parquet"
    if path.exists() and not force:
        return pd.read_parquet(path)
    d = prepare.load(prn)
    rows = []
    for i in range(len(d["date"])):
        y, _ = F.series_for(d, i, "full")
        u = d["u"][i].astype(float)
        f = F.fit_day(y, {"cu": np.cos(u), "su": np.sin(u)})
        if not f.ok:
            continue
        A, B = f.gls["cu"], f.gls["su"]
        phi_u = np.arctan2(B, A)
        usun = np.angle(np.nanmean(np.exp(1j * (u - d["mu"][i].astype(float)))))
        omega = G.perigee_angle(d, i)
        rows.append({"date": d["date"][i], "amp_ps": float(np.hypot(A, B) * 1e12),
                     "phase_sun_deg": float(np.degrees(np.angle(np.exp(1j * (phi_u - usun))))),
                     "phase_perigee_deg": float(np.degrees(np.angle(np.exp(1j * (phi_u - omega))))),
                     "beta_deg": float(np.degrees(np.nanmean(d["beta"][i]))), "q": f.q, "rms_ps": f.rms * 1e12})
    df = pd.DataFrame(rows)
    df.to_parquet(path)
    return df


def folded(prn: str, nbins: int = 72, force: bool = False) -> dict:
    """Clock residuals (per-day offset and rate removed, GR subtracted, no alpha fitted) folded on
    the mean anomaly over the whole segment, with the physical clock's swing for scale.

    Each day contributes one mean per phase bin; bins are combined with a 5-sigma (MAD) clip so
    that a handful of disturbed days cannot dominate a bin."""
    path = CACHE / f"folded_{prn}.npz"
    if path.exists() and not force:
        with np.load(path) as z:
            return {k: z[k] for k in z.files}
    d = prepare.load(prn)
    df = daily(prn, "alpha")
    days = set(segment_days(prn, df)["date"])
    edges = np.linspace(-np.pi, np.pi, nbins + 1)
    res_b, phys_b, tpl_b = [], [], []
    for i, date in enumerate(d["date"]):
        if date not in days:
            continue
        y, tpl = F.series_for(d, i, "full")
        f = F.fit_day(y, {}, keep_resid=True)
        if not f.ok:
            continue
        M = np.angle(np.exp(1j * d["M"][i].astype(float)))
        phys = d["clk"][i] + d["kc"][i]
        ok = f.mask & np.isfinite(phys)
        X = np.vstack([np.ones(ok.sum()), F.TC[ok]]).T
        pc = phys[ok] - X @ np.linalg.lstsq(X, phys[ok], rcond=None)[0]
        tc = tpl["alpha"][ok] - X @ np.linalg.lstsq(X, tpl["alpha"][ok], rcond=None)[0]
        k = np.clip(np.digitize(M[ok], edges) - 1, 0, nbins - 1)
        rows = [np.full(nbins, np.nan) for _ in range(3)]
        for b in np.unique(k):
            sel = k == b
            rows[0][b] = f.resid[ok][sel].mean()
            rows[1][b] = pc[sel].mean()
            rows[2][b] = tc[sel].mean()
        res_b.append(rows[0])
        phys_b.append(rows[1])
        tpl_b.append(rows[2])
    R, P, T = np.array(res_b), np.array(phys_b), np.array(tpl_b)
    med = np.nanmedian(R, axis=0)
    mad = 1.4826 * np.nanmedian(np.abs(R - med), axis=0)
    keep = np.abs(R - med) < 5 * mad
    Rk = np.where(keep, R, np.nan)
    n = np.sum(np.isfinite(Rk), axis=0)
    out = {"centers": 0.5 * (edges[1:] + edges[:-1]), "resid": np.nanmean(Rk, axis=0),
           "resid_sem": np.nanstd(Rk, axis=0) / np.sqrt(n), "phys": np.nanmean(P, axis=0),
           "tpl": np.nanmean(T, axis=0), "n": n, "clipped_fraction": float(1 - keep[np.isfinite(R)].mean())}
    np.savez(path, **out)
    return out


def clock_history(prn: str, force: bool = False) -> pd.DataFrame:
    """Daily clock rate (us/day) against the ground reference and noise level, for the history figure."""
    path = CACHE / f"clockrate_{prn}.parquet"
    if path.exists() and not force:
        return pd.read_parquet(path)
    d = prepare.load(prn)
    rows = []
    t = F.T_SEC
    for i, date in enumerate(d["date"]):
        y = d["clk"][i] + d["kc"][i] - (d["mono"][i] + d["j2"][i] + d["tide"][i] + d["kin"][i])
        ok = np.isfinite(y)
        if ok.sum() < 200:
            continue
        rows.append({"date": date, "rate_us_day": float(np.polyfit(t[ok], y[ok], 1)[0] * 86400e6)})
    df = pd.DataFrame(rows).merge(daily(prn, "alpha")[["date", "q", "rms_ps", "ok"]], on="date", how="left")
    df.to_parquet(path)
    return df


def add_sun1_null(force: bool = False) -> dict:
    """Add the sun1_1yr null test to an existing results file (same code path as run())."""
    path = CACHE / "redshift_results.json"
    out = json.loads(path.read_text())
    res = out["per_satellite"]["E14"]
    if "null_sun1_1yr" in res and not force:
        return res["null_sun1_1yr"]
    d = prepare.load("E14")
    days = segment_days("E14", daily("E14", "alpha"))
    res["null_sun1_1yr"] = null_tests("E14", days, d, model="sun1_1yr")
    path.write_text(json.dumps(out, indent=1, default=float))
    return res["null_sun1_1yr"]
