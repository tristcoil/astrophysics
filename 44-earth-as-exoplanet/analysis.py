"""Every number in the README, computed from the cached data.

Sections follow the questions of the project spec:

1. geometry   - where DSCOVR is; the phase angle EPIC sees Earth at
2. transform  - what the display images' brightness means (tests of gamma)
3. rotation   - the length of the day, from one year of the dot's flicker
4. sampling   - how rarely could we look and still find it (Jiang et al. 2018)
5. surface    - continents and oceans in the dot's colour; rotation vs weather
6. seasons    - the poles tipping into view; the yearly brightness cycle
7. long       - eleven years of thumbnails: solar vs sidereal day, phase angle
8. moon       - the transits, the Moon's brightness, every predicted crossing
9. versions   - whether the archive's two processing versions agree

``run.py`` calls ``run_all`` and saves the results to ``data/results.json``.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle

import ephemeris as E
import epicdata
import geometry as G
import lightcurve as LC
import moon as M
import periodicity as P
import photometry as ph
import surface
from config import CACHE_DIR, JPG_WINDOWS, RESULTS_PATH
from transform import (
    GAMMA_GRID, GAMMA_NOMINAL, GAMMA_RANGE, OCEAN_REFLECTANCE, YANG2018,
)

WINDOWS = {"2015-16": JPG_WINDOWS["yang2018_year"], "2024": JPG_WINDOWS["year2024"]}
JIANG2018_TABLE2 = {  # Jiang et al. 2018, AJ 156, 26, Table 2 (780 nm): interval (h) -> P(detect)
    1.8008: 0.99, 3.6017: 1.00, 5.4025: 1.00, 7.2033: 0.94, 9.0042: 0.93, 10.8050: 0.80,
    12.6058: 0.60, 14.4067: 0.65, 16.2075: 0.56, 18.0083: 0.51, 19.8092: 0.47, 21.6100: 0.23,
    23.4108: 0.32, 25.2117: 0.28, 27.0125: 0.25, 28.8133: 0.12, 30.6142: 0.24, 32.4150: 0.10,
    34.2158: 0.07, 36.0167: 0.07,
}
EPIC_BASE_INTERVAL_S = 6483.0  # Jiang et al.'s nominal EPIC cadence
#: Longitude ranges for "the Pacific-facing" and "the Africa-facing" Earth.
PACIFIC = (-180.0, -140.0)
AFRICA = (0.0, 40.0)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, (pd.Timestamp,)):
        return x.isoformat()
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist())
    return x


def save_results(results: dict) -> None:
    RESULTS_PATH.write_text(json.dumps(_jsonable(results), indent=1))


def _lon_in(lon, rng):
    return (lon >= rng[0]) & (lon < rng[1])


def _detrended(w: pd.DataFrame, col: str, days: int = 15) -> pd.Series:
    """Relative deviation from a running mean of daily means (removes seasons)."""
    day = w["t"].dt.floor("D")
    dm = w.groupby(day)[col].mean()
    trend = dm.rolling(days, center=True, min_periods=5).mean()
    return w[col] / day.map(trend).to_numpy() - 1.0


# ---------------------------------------------------------------------------
# 1. geometry
# ---------------------------------------------------------------------------

def geometry_summary(meta: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    dsc, moon, sun = epicdata.vectors(meta)
    alpha = G.phase_angle_deg(dsc, sun)
    d = G.norm(dsc)
    # DSCOVR in a frame that turns with the Sun: x toward the Sun, z toward
    # ecliptic-ish north (the J2000 pole projected), y completing the triad.
    xs = G.unit(sun)
    zs = G.unit(G.POLE - np.sum(G.POLE * xs, axis=1)[:, None] * xs)
    ys = np.cross(zs, xs)
    track = pd.DataFrame({
        "t": meta["t"], "phase_deg": alpha, "dist_km": d,
        "along_sun_km": np.sum(dsc * xs, axis=1), "y_km": np.sum(dsc * ys, axis=1),
        "z_km": np.sum(dsc * zs, axis=1), "lat": meta["lat"],
    })
    yearly = track.groupby(track["t"].dt.year).agg(
        phase_min=("phase_deg", "min"), phase_max=("phase_deg", "max"),
        dist_min=("dist_km", "min"), dist_max=("dist_km", "max"),
        lat_min=("lat", "min"), lat_max=("lat", "max"), n=("phase_deg", "size"))
    res = {
        "n_images": int(len(meta)), "first": meta["t"].min(), "last": meta["t"].max(),
        "n_days": int(meta["t"].dt.normalize().nunique()),
        "images_per_day_median_by_month": meta.groupby(meta["t"].dt.month)["image"].size().div(
            meta.groupby(meta["t"].dt.month)["t"].apply(lambda s: s.dt.normalize().nunique())).round(1).to_dict(),
        "phase_deg": {"min": float(alpha.min()), "max": float(alpha.max()),
                      "median": float(np.median(alpha))},
        "dist_km": {"min": float(d.min()), "max": float(d.max())},
        "earth_diameter_deg": {"min": float(2 * np.degrees(np.arcsin(6378.137 / d.max()))),
                               "max": float(2 * np.degrees(np.arcsin(6378.137 / d.min())))},
        "subobs_lat": {"min": float(meta["lat"].min()), "max": float(meta["lat"].max())},
        "phase_period_d": _dominant_period(track),
        "yearly": yearly.round(3).to_dict(orient="index"),
        "versions": meta["version"].value_counts().to_dict(),
    }
    return res, track


def _dominant_period(track: pd.DataFrame) -> float:
    """Main period (days) of the daily-mean phase angle: half of DSCOVR's loop."""
    daily = track.groupby(track["t"].dt.floor("D"))["phase_deg"].mean()
    t = (daily.index - daily.index[0]).days.to_numpy().astype(float)
    f = np.linspace(1 / 400, 1 / 20, 20000)
    p = LombScargle(t, daily.to_numpy()).power(f)
    return float(1 / f[np.argmax(p)])


# ---------------------------------------------------------------------------
# 2. transform tests
# ---------------------------------------------------------------------------

def _mixing_gamma(vr, vy, lo=45, hi=240):
    """Fit one exponent to the ocean-cloud mixing locus (see transform.py)."""
    from scipy.optimize import minimize_scalar

    bins = np.arange(16, 246, 6)
    ridge = []
    for b in bins:
        s = (vr >= b) & (vr < b + 6)
        if s.sum() > 200:
            ridge.append((b + 3, np.median(vy[s])))
    ridge = np.array(ridge, float)
    o, c = ridge[1], ridge[-1]  # clear-ocean mode and brightest cloud

    def locus(g, x):
        kr = (c[0] / o[0]) ** g
        ky = (c[1] / o[1]) ** g
        f = ((x / o[0]) ** g - 1) / (kr - 1)
        return o[1] * (1 + f * (ky - 1)) ** (1 / g)

    sel = (ridge[:, 0] >= lo) & (ridge[:, 0] <= hi)
    x, y = ridge[sel, 0], ridge[sel, 1]
    r = minimize_scalar(lambda g: np.sum((locus(g, x) - y) ** 2), bounds=(0.5, 5), method="bounded")
    g = float(r.x)
    return {"gamma": g, "rms_dn": float(np.sqrt(r.fun / len(x))),
            "implied_ocean_colour_vs_cloud": float((c[0] / o[0]) ** g / (c[1] / o[1]) ** g),
            "ridge": ridge.tolist()}


def transform_tests(jpg: pd.DataFrame) -> tuple[dict, dict]:
    """The five tests of the display exponent described in transform.py."""
    w = LC.window(jpg, *YANG2018["window"])
    day = w["t"].dt.date.to_numpy()
    out = {"yang_window_images": int(len(w))}
    # 3/4. day-to-day scatter and peak-to-peak vs Yang et al. 2018
    stats = {}
    for g in GAMMA_GRID:
        row = {}
        for ch in ("red", "green", "blue"):
            s = pd.Series(w[ph.moment_column("inner", ch, g)].to_numpy()).groupby(day)
            row[ch] = {"rel_std": float((s.std() / s.mean()).mean()),
                       "peak_to_peak": float(((s.max() - s.min()) / s.mean()).mean())}
        stats[f"{g:.1f}"] = row
    out["yang_comparison"] = stats

    def solve(stat, ch):
        gs = np.array(GAMMA_GRID)
        vals = np.array([stats[f"{g:.1f}"][ch][stat] for g in gs])
        target = YANG2018["daily_rel_std" if stat == "rel_std" else "daily_peak_to_peak"][ch]
        return float(np.interp(target, vals, gs)) if vals.min() <= target <= vals.max() else None

    out["gamma_from_yang_std"] = {ch: solve("rel_std", ch) for ch in ("red", "blue")}
    out["gamma_from_yang_p2p"] = {ch: solve("peak_to_peak", ch) for ch in ("red", "blue")}

    # 1/2/5: pixel-level tests on Pacific-facing images (mostly ocean and cloud)
    pac = w[(w["lon"] < -140) | (w["lon"] > 170)]
    vr, vg, vb = [], [], []
    ocean_ratio, cloud_ratio = [], []
    for im in pac["image"]:
        rgb = ph.load_rgb(epicdata.image_path("jpg", im)).astype(np.int16)
        m = ph.measure(rgb.astype(np.uint8))
        yy, xx = np.mgrid[:rgb.shape[0], :rgb.shape[1]]
        central = np.hypot(xx - m["cx"], yy - m["cy"]) < 0.5 * m["r_extent"]
        r, g_, b = rgb[..., 0][central], rgb[..., 1][central], rgb[..., 2][central]
        vr.append(r), vg.append(g_), vb.append(b)
        oc = r < np.percentile(r, 3)
        cl = (r > np.percentile(r, 99)) & (np.abs(r - g_) < 12)
        ocean_ratio.append(np.median(b[oc] / np.maximum(r[oc], 1)))
        cloud_ratio.append(np.median(b[cl] / np.maximum(r[cl], 1)))
    vr, vg, vb = (np.concatenate(a).astype(float) for a in (vr, vg, vb))
    out["mixing_blue_vs_red"] = _mixing_gamma(vr, vb)
    out["mixing_green_vs_red"] = _mixing_gamma(vr, vg)
    vo, vc = float(np.median(ocean_ratio)), float(np.median(cloud_ratio))
    out["ocean_cloud_colour"] = {
        "display_blue_red_ocean": vo, "display_blue_red_cloud": vc,
        # reflectance blue/red: ocean 4.5-7 (Rayleigh + water), cloud 1.0-1.12
        "gamma_range": [float(np.log(4.5 / 1.12) / np.log(vo / vc)), float(np.log(7.0 / 1.0) / np.log(vo / vc))],
        "gamma_nominal": float(np.log(5.7 / 1.05) / np.log(vo / vc)),
    }
    # 5: darkest ocean against the disc mean
    from scipy.optimize import brentq
    hist = _hist_window("jpg", w["image"])
    H = hist[:, 1].sum(axis=0).astype(float)  # inner disc, all images
    v = np.arange(256) / 255.0
    dark = {}
    for ci, ch in ((0, "red"), (2, "blue")):
        cdf = np.cumsum(H[ci]) / H[ci].sum()
        target = OCEAN_REFLECTANCE[ch] / YANG2018["mean"][ch]
        gs = []
        for p in (0.005, 0.01, 0.02):
            vp = np.searchsorted(cdf, p) / 255.0
            gs.append(brentq(lambda g: vp**g / ((H[ci] * v**g).sum() / H[ci].sum()) - target, 0.3, 6))
        dark[ch] = {"gamma_min": float(min(gs)), "gamma_max": float(max(gs))}
    out["dark_ocean"] = dark
    out["adopted"] = {"gamma": GAMMA_NOMINAL, "range": list(GAMMA_RANGE)}
    # Earth's July-2016 mean linear light per exponent: the reflectance scale
    earth_ref = {f"{g:.1f}": {ch: float(w[ph.moment_column("inner", ch, g)].mean())
                              for ch in ph.CHANNELS} for g in GAMMA_GRID}
    out["reflectance_scale"] = {f"{g:.1f}": {ch: YANG2018["mean"][ch] / earth_ref[f"{g:.1f}"][ch]
                                             for ch in YANG2018["mean"]} for g in GAMMA_GRID}
    art = {"mixing_vr": vr, "mixing_vg": vg, "mixing_vb": vb, "earth_ref": earth_ref}
    return out, art


def _hist_window(kind: str, images) -> np.ndarray:
    z = np.load(ph.table_paths(kind)[1], allow_pickle=False)
    idx = pd.Index(z["image"]).get_indexer(list(images))
    return z["hist"][idx[idx >= 0]]


# ---------------------------------------------------------------------------
# 3. rotation
# ---------------------------------------------------------------------------

def rotation(jpg: pd.DataFrame) -> tuple[dict, dict]:
    out, art = {}, {}
    for name, (a, b) in WINDOWS.items():
        w = LC.window(jpg, a, b)
        t = LC.days_since(w["t"])
        months = w["t"].dt.strftime("%Y-%m").to_numpy()
        res = {"n_images": int(len(w)), "n_days": int(w["t"].dt.normalize().nunique())}
        for ch in ("red", "green", "blue"):
            y = w[f"F_{ch}"].to_numpy()
            jk = P.jackknife_peak(t, y, months, 1.0, 0.05)
            f, p = P.periodogram(t, y, 0.02, 6.0, oversample=5)
            top = f[np.argmax(p)]
            grid_t, grid = P.hourly_grid(t, y)
            r = P.acf(grid, 24 * 12)
            lag = np.arange(len(r)).astype(float)
            acf_p, acf_r = P.acf_peak(lag, r)
            harmonics = {}
            for k in (1, 2, 3, 4):
                fk, pk = P.peak_near(f, p, float(k), 0.03)
                harmonics[f"{24 / k:.0f}h"] = {"f": fk, "power": pk}
            res[ch] = {"period_h": jk["period_h"], "period_err_h": jk["period_err_h"],
                       "offset_s": (jk["period_h"] - 24.0) * 3600, "offset_err_s": jk["period_err_h"] * 3600,
                       "strongest_peak_period_h": 24.0 / top, "acf_peak_h": acf_p, "acf_r": acf_r,
                       "harmonics": harmonics}
            art[(name, ch)] = {"f": f, "p": p, "acf_lag_h": lag, "acf": r}
        out[name] = res
    return out, art


# ---------------------------------------------------------------------------
# 4. sampling
# ---------------------------------------------------------------------------

def sampling(jpg: pd.DataFrame, rng: np.random.Generator, n_segments: int = 300) -> dict:
    out = {"jiang2018": {f"{k:.2f}": v for k, v in JIANG2018_TABLE2.items()}}
    band = (1 - 0.5 / 7, 1 + 0.5 / 7)  # half a resolution element of a week-long segment
    for name, (a, b) in WINDOWS.items():
        w = LC.window(jpg, a, b)
        t = LC.days_since(w["t"])
        for mode in ("random", "regular"):
            out[f"{name} {mode}"] = P.sampling_test(t, w["F_red"].to_numpy(), EPIC_BASE_INTERVAL_S,
                                                    range(1, 21), rng, mode=mode,
                                                    n_segments=n_segments, band=band)
    out["band"] = list(band)
    return out


# ---------------------------------------------------------------------------
# 5. surface: colour, land, weather
# ---------------------------------------------------------------------------

def surface_analysis(jpg: pd.DataFrame) -> tuple[dict, dict]:
    out, art = {}, {}
    for name, (a, b) in WINDOWS.items():
        w = LC.window(jpg, a, b).copy()
        res = {}
        # longitude profiles (10-degree bins of the sub-DSCOVR longitude)
        w["lon_bin"] = (np.floor(w["lon"] / 10) * 10 + 5).astype(int)
        cols = ["F_red", "F_green", "F_blue", "blue_red", "land_flux", "polar_flux"]
        prof = w.groupby("lon_bin")[cols].agg(["mean", "std", "count"])
        art[(name, "profile")] = prof
        pac, afr = w[_lon_in(w["lon"], PACIFIC)], w[_lon_in(w["lon"], AFRICA)]
        res["pacific_vs_africa"] = {
            c: {"pacific": float(pac[c].mean()), "africa": float(afr[c].mean()),
                "pacific_over_africa": float(pac[c].mean() / afr[c].mean())}
            for c in ["F_red", "F_green", "F_blue", "blue_red", "land_flux"]}
        # same contrast in display units and at the ends of the gamma range
        for g in (1.0, GAMMA_RANGE[0], GAMMA_RANGE[1]):
            res[f"pacific_over_africa_g{g:.1f}"] = {
                ch: float(pac[f"F_{ch}_g{g:.1f}"].mean() / afr[f"F_{ch}_g{g:.1f}"].mean())
                for ch in ("red", "green", "blue")}
            br = pac[f"F_blue_g{g:.1f}"] / pac[f"F_red_g{g:.1f}"]
            ba = afr[f"F_blue_g{g:.1f}"] / afr[f"F_red_g{g:.1f}"]
            res[f"pacific_over_africa_g{g:.1f}"]["blue_red"] = float(br.mean() / ba.mean())
        res["corr_with_land"] = {c: float(np.corrcoef(w[c], w["land_flux"])[0, 1])
                                 for c in ["F_red", "F_green", "F_blue", "blue_red"]}
        res["corr_with_polar"] = {c: float(np.corrcoef(w[c], w["polar_flux"])[0, 1])
                                  for c in ["F_red", "F_green", "F_blue", "blue_red"]}
        fit = np.polyfit(w["land_flux"], w["blue_red"], 1)
        res["blue_red_vs_land_fit"] = {"slope": float(fit[0]), "intercept": float(fit[1]),
                                       "rms": float(np.std(w["blue_red"] - np.polyval(fit, w["land_flux"])))}
        # rotation vs weather: after removing seasons, how much variance repeats
        # with longitude (per quarter, since the pattern drifts with season)
        w["quarter"] = w["t"].dt.quarter
        split = {}
        for c in ["F_red", "F_green", "F_blue", "blue_red"]:
            r = _detrended(w, c)
            prof_q = r.groupby([w["quarter"], w["lon_bin"]]).transform("mean")
            vt, vr = r.var(), prof_q.var()
            split[c] = {"rel_std": float(np.sqrt(vt)), "rotation_share": float(vr / vt),
                        "weather_share": float((r - prof_q).var() / vt)}
        res["rotation_vs_weather"] = split
        # principal components of the three colours (Fan et al. 2019)
        X = np.log(w[["F_red", "F_green", "F_blue"]].to_numpy())
        X = (X - X.mean(0)) / X.std(0)
        U, S, Vt = np.linalg.svd(X, full_matrices=False)
        pcs = U * S
        sign = np.sign(Vt[:, 0])  # make PC1 "brighter = positive"
        sign[1] = np.sign(np.corrcoef(pcs[:, 1], w["land_flux"])[0, 1]) or 1
        pcs, Vt = pcs * sign, Vt * sign[:, None]
        res["pca"] = {"explained": (S**2 / np.sum(S**2)).tolist(), "loadings": Vt.tolist(),
                      "corr_pc1_land": float(np.corrcoef(pcs[:, 0], w["land_flux"])[0, 1]),
                      "corr_pc1_polar": float(np.corrcoef(pcs[:, 0], w["polar_flux"])[0, 1]),
                      "corr_pc2_land": float(np.corrcoef(pcs[:, 1], w["land_flux"])[0, 1]),
                      "corr_pc2_polar": float(np.corrcoef(pcs[:, 1], w["polar_flux"])[0, 1])}
        art[(name, "pcs")] = pd.DataFrame({"t": w["t"].to_numpy(), "pc1": pcs[:, 0], "pc2": pcs[:, 1],
                                           "land": w["land_flux"].to_numpy(), "lon": w["lon"].to_numpy()})
        out[name] = res
    return out, art


# ---------------------------------------------------------------------------
# 6. seasons
# ---------------------------------------------------------------------------

def seasons(jpg: pd.DataFrame) -> tuple[dict, dict]:
    out, art = {}, {}
    for name, (a, b) in WINDOWS.items():
        w = LC.window(jpg, a, b)
        daily = w.groupby(w["t"].dt.floor("D")).agg(
            F_red=("F_red", "mean"), F_green=("F_green", "mean"), F_blue=("F_blue", "mean"),
            blue_red=("blue_red", "mean"), lat=("lat", "mean"), polar=("polar_flux", "mean"),
            land=("land_flux", "mean"), n=("F_red", "size"))
        daily = daily[daily["n"] >= 5]  # as Yang et al.: at least five images a day
        monthly = w.groupby(w["t"].dt.month)[["F_red", "F_green", "F_blue", "blue_red", "lat",
                                              "polar_flux", "land_flux"]].mean()
        sm = daily["F_red"].rolling(31, center=True, min_periods=15).mean()
        res = {"monthly": monthly.round(4).to_dict(orient="index"),
               "brightest_month": int(monthly["F_red"].idxmax()),
               "faintest_month": int(monthly["F_red"].idxmin()),
               "seasonal_range_red": float(sm.max() / sm.min() - 1),
               "corr_daily_red_polar": float(np.corrcoef(daily["F_red"], daily["polar"])[0, 1]),
               "subobs_lat_range": [float(w["lat"].min()), float(w["lat"].max())]}
        # Yang et al. 2018's seasonal peaks: December and early June
        jun = daily[(daily.index.month == 6)]["F_red"].mean()
        aprsep = daily[daily.index.month.isin([4, 8, 9])]["F_red"].mean()
        res["june_over_apr_aug_sep"] = float(jun / aprsep)
        out[name] = res
        art[name] = daily
    return out, art


# ---------------------------------------------------------------------------
# 7. eleven years of thumbnails
# ---------------------------------------------------------------------------

def long_baseline(thumbs: pd.DataFrame) -> tuple[dict, dict]:
    w = thumbs[thumbs["ok"]]
    t = LC.days_since(w["t"])
    out, art = {"n_images": int(len(w)), "span_days": float(t.max() - t.min()),
                "first": w["t"].min(), "last": w["t"].max()}, {}
    sid = 24.0 / P.SIDEREAL_DAY_H  # cycles per solar day
    for ch, col in (("red", "F_red_g1.0"), ("blue", "F_blue_g1.0")):
        y = w[col].to_numpy()
        f = np.linspace(0.995, 1.006, 20001)
        p = LombScargle(t, y).power(f)
        fpk, ppk = P.peak_near(f, p, 1.0, 0.004)
        # power at the sidereal frequency relative to the solar peak
        psid = float(np.interp(sid, f, p))
        # error from leaving out one year at a time
        years = w["t"].dt.year.to_numpy()
        jk = P.jackknife_peak(t, y, years, 1.0, 0.002, oversample=10)
        out[ch] = {"peak_f": fpk, "peak_period_h": 24.0 / fpk, "period_err_h": jk["period_err_h"],
                   "offset_s": (24.0 / fpk - 24.0) * 3600, "offset_err_s": jk["period_err_h"] * 3600,
                   "sidereal_period_h": P.SIDEREAL_DAY_H, "power_peak": ppk,
                   "power_at_sidereal": psid, "resolution_s": 24 * 3600 / (t.max() - t.min()) / 1.0}
        art[ch] = {"f": f, "p": p}
        # broad periodogram for the low-frequency (seasonal) structure
        fb = np.linspace(1 / 1000, 4.2, 400000)
        art[(ch, "broad")] = {"f": fb, "p": LombScargle(t, y).power(fb)}
    return out, art


def phase_dependence(thumbs: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    """How the dot's daily-mean brightness depends on the phase angle.

    The Sun-Earth-DSCOVR angle swings with DSCOVR's ~6-month loop around L1,
    while seasons follow the calendar; over eleven years the two drift apart,
    so a regression with seasonal harmonics separates them."""
    w = thumbs[thumbs["ok"]].copy()
    daily = w.groupby(w["t"].dt.floor("D")).agg(
        F_red=("F_red_g1.0", "mean"), F_green=("F_green_g1.0", "mean"), F_blue=("F_blue_g1.0", "mean"),
        F_red_lin=(f"F_red_g{GAMMA_NOMINAL:.1f}", "mean"), F_blue_lin=(f"F_blue_g{GAMMA_NOMINAL:.1f}", "mean"),
        phase=("phase_deg", "mean"), n=("phase_deg", "size"), version=("version", lambda s: (s == "04").mean()))
    daily = daily[daily["n"] >= 8].copy()
    doy = daily.index.dayofyear.to_numpy() / 365.25
    yrs = (daily.index - pd.Timestamp("2020-01-01")).days.to_numpy() / 365.25
    out = {"n_days": int(len(daily)), "phase_range": [float(daily["phase"].min()), float(daily["phase"].max())]}
    for col in ("F_red", "F_green", "F_blue", "F_red_lin", "F_blue_lin"):
        y = np.log(daily[col].to_numpy())
        cols = [np.ones_like(y), daily["phase"].to_numpy(), yrs, daily["version"].to_numpy()]
        for k in (1, 2, 3):
            cols += [np.cos(2 * np.pi * k * doy), np.sin(2 * np.pi * k * doy)]
        A = np.column_stack(cols)
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        resid = y - A @ coef
        # errors: residuals are correlated over weeks; inflate by sqrt(7)
        cov = np.linalg.inv(A.T @ A) * np.var(resid) * 7.0
        out[col] = {"slope_per_deg": float(coef[1]), "slope_err": float(np.sqrt(cov[1, 1])),
                    "trend_per_year": float(coef[2]), "trend_err": float(np.sqrt(cov[2, 2])),
                    "v04_offset": float(coef[3]), "v04_offset_err": float(np.sqrt(cov[3, 3]))}
        daily[col + "_season_free"] = y - A[:, 2:] @ coef[2:]
    # Lambert sphere for comparison
    out["lambert_slope_per_deg_at_median"] = float(np.gradient(
        np.log(G.lambert_phase(np.array([6.0, 7.0, 8.0]))))[1])
    return out, daily


# ---------------------------------------------------------------------------
# 8. the Moon
# ---------------------------------------------------------------------------

def moon_analysis(earth_ref: dict) -> tuple[dict, dict]:
    path = CACHE_DIR / "transits.parquet"
    if path.exists():
        tab = pd.read_parquet(path)
    else:
        tab = pd.concat([M.measure_transit(k) for k in M.TRANSITS], ignore_index=True)
        tab.to_parquet(path)
    summary = M.summarise(tab, earth_ref)
    out = {"transits": summary, "expected": M.expected_ratios(), "stray_light": M.stray_light(tab),
           "gamma_matching_moon": {k: {ch: M.gamma_matching_moon(summary, k, ch) for ch in ("red", "blue")}
                                   for k in summary}}
    ev_path = CACHE_DIR / "moon_events.parquet"
    events = pd.read_parquet(ev_path) if ev_path.exists() else E.mission_events()
    if not ev_path.exists():
        events.to_parquet(ev_path)
    out["events"] = {"n_transits": int((events["kind"] == "transit").sum()),
                     "n_occultations": int((events["kind"] == "occultation").sum()),
                     "n_near_misses": int((events["kind"] == "near miss").sum())}
    out["imaged"] = event_coverage()
    return out, {"transits": tab, "events": events}


def archive_gaps(meta: pd.DataFrame, min_days: int = 7) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
    """Stretches of more than ``min_days`` without any natural-colour image."""
    days = pd.Series(sorted(meta["t"].dt.normalize().unique()))
    d = days.diff().dt.days
    return [(days[i - 1], days[i]) for i in d.index[d > min_days]]


def event_coverage() -> dict:
    """Did EPIC image the predicted Moon events? Granule counts of the login-
    gated L1A 'miscellaneous' collection (NASA CMR, cached by fetch_data.py)
    in each event window, against random 4-hour windows away from events."""
    from fetch_data import CMR_COVERAGE_PATH

    if not CMR_COVERAGE_PATH.exists():
        return {}
    cov = json.loads(CMR_COVERAGE_PATH.read_text())
    gaps = archive_gaps(epicdata.load_metadata())
    ev = pd.DataFrame(cov["events"])
    start = pd.to_datetime(ev["start"]).dt.tz_localize(None)
    ev["in_gap"] = [any(a < s < b for a, b in gaps) for s in start]
    use = ev[~ev["in_gap"]]
    ctl = pd.DataFrame(cov["controls"])
    return {"n_events": int(len(ev)), "n_in_archive_gaps": int(ev["in_gap"].sum()),
            "n_outside_gaps": int(len(use)), "n_with_misc_granules": int((use["L1A_misc"] > 0).sum()),
            "median_misc_granules": float(use.loc[use["L1A_misc"] > 0, "L1A_misc"].median()),
            "max_misc_granules": int(use["L1A_misc"].max()),
            "n_control_windows": int(len(ctl)), "n_controls_with_granules": int((ctl["L1A_misc"] > 0).sum()),
            "max_control_granules": int(ctl["L1A_misc"].max()),
            "archive_gaps": [(a.date().isoformat(), b.date().isoformat()) for a, b in gaps]}


def version_offset(jpg: pd.DataFrame) -> dict:
    """v04 vs v03 images at matched longitude and date (2024 has both)."""
    w = jpg[jpg["ok"]].copy()
    w["lon_bin"] = np.floor(w["lon"] / 15)
    w["week"] = w["t"].dt.to_period("W")
    ref = w[w["version"] == "03"].groupby(["week", "lon_bin"])["F_red"].mean()
    v4 = w[w["version"] == "04"]
    key = list(zip(v4["week"], v4["lon_bin"]))
    r = v4["F_red"].to_numpy() / ref.reindex(key).to_numpy()
    r = r[np.isfinite(r)]
    return {"n_v04_matched": int(len(r)), "median_ratio": float(np.median(r)),
            "spread": float(np.std(r)), "sem": float(np.std(r) / np.sqrt(max(len(r), 1)))}


# ---------------------------------------------------------------------------

def run_all(rng: np.random.Generator, *, quiet: bool = False) -> tuple[dict, dict]:
    meta = epicdata.load_metadata()
    results, art = {}, {}
    results["geometry"], art["track"] = geometry_summary(meta)

    jpg = LC.load("jpg")
    jpg = surface.add_to(jpg)
    art["jpg"] = jpg
    results["transform"], art["transform"] = transform_tests(jpg)
    results["rotation"], art["rotation"] = rotation(jpg)
    results["sampling"] = sampling(jpg, rng)
    results["surface"], art["surface"] = surface_analysis(jpg)
    results["seasons"], art["seasons"] = seasons(jpg)
    results["versions"] = version_offset(jpg)

    thumbs = LC.load("thumbs")
    art["thumbs"] = thumbs
    results["long"], art["long"] = long_baseline(thumbs)
    results["phase"], art["phase"] = phase_dependence(thumbs)

    results["moon"], art["moon"] = moon_analysis(art["transform"]["earth_ref"])
    return results, art
