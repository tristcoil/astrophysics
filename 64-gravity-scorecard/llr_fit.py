"""Fit the APOLLO range residuals: nuisance parameters plus three relativity tests.

For each lunar ephemeris (DE440 and, as an independent check, EPM2021) the
model of ``llr_model`` is evaluated once; its residuals are then fitted by
weighted least squares with

* nuisance parameters: a range bias per APOLLO hardware group (a-f, as APOLLO
  recommends), the station position (east, north, up) and plate velocity, the
  reflector coordinates (Apollo 15's x held fixed, since a common x-shift equals
  a bias), and a slowly varying along-track shift of the Moon (a time offset
  quadratic in time) -- see the README for why DE440 needs it;
* relativity parameters, each zero (or one) in GR:
  - the Sun's and Earth's Shapiro delay scaled by (1 + gamma)/2: it adds 7.5 m
    to every pulse's one-way range and varies by +-40 cm with the distances;
  - the change of TDB - TT during the 2.5-s round trip (up to +-18 cm), scaled
    by a factor k (GR: k = 1);
  - the Nordtvedt term eta * C * cos D with our simulated C = 12.98 m (GR: eta = 0).

The ephemeris itself was integrated with GR (beta = gamma = 1) and fitted to
LLR data assuming GR, so these are consistency tests of the ranges against a GR
solar-system model, not independent measurements: a real test re-fits the whole
ephemeris with the parameter free (Williams et al.; Hofmann & Mueller 2018).

Errors: the post-fit scatter (2-4 cm) is far above APOLLO's millimetre
statistical errors and is correlated within a night, so a model-error floor is
added in quadratrure (chosen to make chi^2/dof = 1) and the quoted uncertainties
come from a bootstrap over nights.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

import apollo
import config as C
import ephem
import llr_model
import nordtvedt

EPHEMERIDES = {"DE440": ephem.default, "EPM2021": ephem.epm2021}
REFLECTORS = ["Apollo 11", "Apollo 14", "Apollo 15", "Lunokhod 1", "Lunokhod 2"]
DT_STEP = 1e-4          # s, finite-difference step for the along-track template
N_BOOT = 400
SEED = 20261019


def model_table(name: str, force: bool = False) -> pd.DataFrame:
    """Nominal model and the along-track template for one ephemeris (cached)."""
    path = C.DERIVED_DIR / f"llr_model_{name}.pkl"
    if path.exists() and not force:
        return pd.read_pickle(path)
    d = apollo.load()
    eph = EPHEMERIDES[name]()
    m = llr_model.compute(d, ephemeris=eph)
    m2 = llr_model.compute(d, params={"eph_dt": DT_STEP}, ephemeris=eph, components=False)
    m["tmpl_dt"] = (m2["resid_m"] - m["resid_m"]) / DT_STEP           # m per s of ephemeris time
    for col in ("group", "reflector", "year", "month", "day", "sigma_s", "rt_s"):
        m[col] = d[col].to_numpy()
    m.to_pickle(path)
    return m


def design(m: pd.DataFrame, *, along_track: bool = True, reflectors: bool = True,
           relativity: bool = True, nordtvedt_c: float | None = None) -> tuple[np.ndarray, list[str]]:
    t = m["years"].to_numpy() - C.LLR_REFERENCE_EPOCH
    cols, names = [], []
    for g in "abcdef":
        cols.append((m["group"] == g).to_numpy(float))
        names.append(f"bias_{g}")
    for k, lab in (("d_east", "E"), ("d_north", "N"), ("d_up", "U")):
        cols.append(m[k].to_numpy())
        names.append(f"sta_{lab}")
    for k, lab in (("d_east", "E"), ("d_north", "N"), ("d_up", "U")):
        cols.append(m[k].to_numpy() * t)
        names.append(f"vel_{lab}")
    if reflectors:
        for r in REFLECTORS:
            sel = (m["reflector"] == r).to_numpy(float)
            for ax in "xyz":
                if r == "Apollo 15" and ax == "x":
                    continue
                cols.append(m[f"d_refl_{ax}"].to_numpy() * sel)
                names.append(f"refl_{r.replace(' ', '')}_{ax}")
    if along_track:
        T = m["tmpl_dt"].to_numpy()
        for p in range(3):
            cols.append(T * t**p)
            names.append(f"dt_{p}")
    if relativity:
        cols.append(-0.5 * (m["shapiro_sun_m"] + m["shapiro_earth_m"]).to_numpy())
        names.append("gamma_minus_1")
        cols.append(-m["tdb_term_m"].to_numpy())
        names.append("tdb_scale_minus_1")
        c_n = nordtvedt_c if nordtvedt_c is not None else nordtvedt_coefficient()
        cols.append(c_n * np.cos(m["D_rad"].to_numpy()))
        names.append("eta_nordtvedt")
        cols.append(c_n * np.sin(m["D_rad"].to_numpy()))
        names.append("eta_sinD")
    return np.column_stack(cols), names


def nordtvedt_coefficient() -> float:
    return nordtvedt.amplitude()["cosD"]


def wls(A, y, sig, keep):
    w = 1.0 / sig[keep]
    x, *_ = np.linalg.lstsq(A[keep] * w[:, None], y[keep] * w, rcond=None)
    return x


def fit(m: pd.DataFrame, *, floor: float | None = None, n_boot: int = N_BOOT, **kw) -> dict:
    y = m["resid_m"].to_numpy()
    A, names = design(m, **kw)
    stat = apollo.scaled_sigma(m) * C.C_LIGHT / 2.0
    keep = np.ones(len(y), bool)
    fl = C.LLR_MODEL_FLOOR_M if floor is None else floor
    for _ in range(10):                          # outlier rejection + floor tuning
        sig = np.hypot(stat, fl)
        x = wls(A, y, sig, keep)
        r = y - A @ x
        s_rob = 1.4826 * np.median(np.abs(r[keep] - np.median(r[keep])))
        new_keep = np.abs(r) < C.LLR_OUTLIER_SIGMA * max(s_rob, 0.005)
        chi2 = np.sum((r[new_keep] / sig[new_keep]) ** 2) / (new_keep.sum() - A.shape[1])
        if floor is None:                        # scale the floor so chi2/dof -> 1
            fl = float(np.sqrt(max(fl**2 * chi2 + (chi2 - 1) * np.median(stat) ** 2, 1e-6)))
        if (new_keep == keep).all() and abs(chi2 - 1) < 0.02:
            break
        keep = new_keep
    sig = np.hypot(stat, fl)
    x = wls(A, y, sig, keep)
    r = y - A @ x
    Aw = A[keep] / sig[keep, None]
    cov = np.linalg.pinv(Aw.T @ Aw)
    formal = np.sqrt(np.diag(cov))
    # bootstrap over nights (residual errors are correlated within a night)
    night = (m["year"] * 10000 + m["month"] * 100 + m["day"]).to_numpy()
    nights = np.unique(night[keep])
    rng = np.random.default_rng(SEED)
    idx_by_night = {n: np.flatnonzero(keep & (night == n)) for n in nights}
    boots = []
    for _ in range(n_boot):
        pick = rng.choice(nights, size=len(nights), replace=True)
        rows = np.concatenate([idx_by_night[n] for n in pick])
        w = 1.0 / sig[rows]
        xb, *_ = np.linalg.lstsq(A[rows] * w[:, None], y[rows] * w, rcond=None)
        boots.append(xb)
    boots = np.array(boots)
    boot_sd = boots.std(axis=0, ddof=1)
    out = {
        "n_used": int(keep.sum()), "n_rejected": int((~keep).sum()), "n_nights": int(len(nights)),
        "floor_m": fl, "rms_m": float(np.sqrt(np.mean(r[keep] ** 2))),
        "wrms_m": float(np.sqrt(np.sum((r[keep] / sig[keep]) ** 2) / np.sum(1 / sig[keep] ** 2))),
        "params": {n: {"value": float(v), "formal": float(f), "boot": float(b)}
                   for n, v, f, b in zip(names, x, formal, boot_sd)},
    }
    t = m["years"].to_numpy() - C.LLR_REFERENCE_EPOCH
    if "dt_0" in names:
        c = [x[names.index(f"dt_{p}")] for p in range(3)]
        out["along_track_ms_by_year"] = {int(yy): float(1e3 * (c[0] + c[1] * (yy + 0.5 - C.LLR_REFERENCE_EPOCH)
                                                                 + c[2] * (yy + 0.5 - C.LLR_REFERENCE_EPOCH) ** 2))
                                          for yy in range(2006, 2021)}
    out["_resid"] = r
    out["_keep"] = keep
    return out


def along_track_by_year(m: pd.DataFrame, floor: float) -> dict:
    """Year-by-year along-track (time) offset of the ephemeris Moon preferred by the ranges (ms),
    with the other nuisance parameters fitted jointly; errors from a bootstrap over nights."""
    A, names = design(m, along_track=False, relativity=False)
    T = m["tmpl_dt"].to_numpy()
    yr = np.floor(m["years"].to_numpy()).astype(int)
    years = sorted(set(yr))
    cols = [T * (yr == y) for y in years]
    A = np.column_stack([A] + cols)
    y = m["resid_m"].to_numpy()
    keep = m["used"].to_numpy() if "used" in m else np.ones(len(y), bool)
    sig = np.hypot(apollo.scaled_sigma(m) * C.C_LIGHT / 2.0, floor)
    x = wls(A, y, sig, keep)
    night = (m["year"] * 10000 + m["month"] * 100 + m["day"]).to_numpy()
    nights = np.unique(night[keep])
    idx = {n: np.flatnonzero(keep & (night == n)) for n in nights}
    rng = np.random.default_rng(SEED + 1)
    boots = []
    for _ in range(200):
        rows = np.concatenate([idx[n] for n in rng.choice(nights, size=len(nights), replace=True)])
        w = 1 / sig[rows]
        xb, *_ = np.linalg.lstsq(A[rows] * w[:, None], y[rows] * w, rcond=None)
        boots.append(xb[-len(years):])
    sd = np.array(boots).std(axis=0, ddof=1)
    # the template is d(resid)/d(dt); the Moon that fits is the ephemeris evaluated at t - x
    return {int(yy): {"offset_ms": float(-1e3 * x[len(names) + k]), "err_ms": float(1e3 * sd[k]),
                      "n": int(((yr == yy) & keep).sum())} for k, yy in enumerate(years)}


RESULTS_CACHE = C.DERIVED_DIR / "llr_results.json"
GPS_FILE = C.RAW_DIR / "gps" / "midas.IGS14.txt"


def gps_velocity(station: str = "P027") -> dict | None:
    """East/north/up velocity (mm/yr) of a GPS station from the NGL MIDAS table (IGS14 frame)."""
    if not GPS_FILE.exists():
        return None
    for line in GPS_FILE.read_text().splitlines():
        f = line.split()
        if f and f[0] == station:
            return {"station": station, "years": [float(f[2]), float(f[3])],
                    "v_mm_yr": [1e3 * float(f[8]), 1e3 * float(f[9]), 1e3 * float(f[10])],
                    "sigma_mm_yr": [1e3 * float(f[11]), 1e3 * float(f[12]), 1e3 * float(f[13])],
                    "lat_lon_h": [float(f[24]), float(f[25]), float(f[26])]}
    return None


def plate_motion_check(res: dict) -> dict:
    """Our fitted APOLLO velocity vs GPS station P027 (2.8 km away), mm/yr."""
    gps = gps_velocity()
    out = {"gps": gps}
    for name in EPHEMERIDES:
        p = res[name]["params"]
        v = [1e3 * p[f"vel_{k}"]["value"] for k in "ENU"]
        s = [1e3 * p[f"vel_{k}"]["boot"] for k in "ENU"]
        entry = {"v_mm_yr": v, "sigma_mm_yr": s}
        if gps:
            entry["diff_sigma"] = [(a - b) / np.hypot(sa, sb) for a, b, sa, sb in
                                   zip(v, gps["v_mm_yr"], s, gps["sigma_mm_yr"])]
        out[name] = entry
    return out


def run_all(force: bool = False) -> dict:
    tables_path = C.DERIVED_DIR / "llr_fit_tables.pkl"
    if RESULTS_CACHE.exists() and tables_path.exists() and not force:
        res = json.loads(RESULTS_CACHE.read_text())
        for name in EPHEMERIDES:     # json turns integer year keys into strings
            res[name]["along_track_by_year"] = {int(k): v for k, v in res[name]["along_track_by_year"].items()}
        res["plate_motion_check"] = plate_motion_check(res)
        return res
    results, tables = {}, {}
    for name in EPHEMERIDES:
        m = model_table(name, force=force)
        tables[name] = m
        res = fit(m)
        m = m.copy()
        m["post_resid_m"] = res.pop("_resid")
        m["used"] = res.pop("_keep")
        tables[name] = m
        results[name] = res
        # the same fit without the along-track freedom (shows the DE440 drift)
        r0 = fit(m, along_track=False, n_boot=50)
        r0.pop("_resid"); r0.pop("_keep")
        results[name + "_no_along_track"] = {k: r0[k] for k in ("rms_m", "wrms_m", "floor_m", "n_used")}
        results[name]["along_track_by_year"] = along_track_by_year(m, res["floor_m"])
        results[name]["rms_by_group_m"] = {g: float(np.sqrt(np.mean(m.loc[m["used"] & (m["group"] == g), "post_resid_m"] ** 2)))
                                           for g in "abcdef"}
        res_n = m.loc[m["used"]]
        night_mean = res_n.groupby(res_n["year"] * 10000 + res_n["month"] * 100 + res_n["day"])["post_resid_m"]
        results[name]["night_mean_rms_m"] = float(np.sqrt(np.mean(night_mean.mean() ** 2)))
        results[name]["within_night_sd_median_m"] = float(night_mean.std().median())
    results["nordtvedt_coefficient_m"] = nordtvedt_coefficient()
    results["plate_motion_check"] = plate_motion_check(results)
    pd.to_pickle(tables, tables_path)
    RESULTS_CACHE.write_text(json.dumps(results, default=float))
    return results


if __name__ == "__main__":
    res = run_all()
    for name in EPHEMERIDES:
        r = res[name]
        print(f"{name}: n={r['n_used']} rms {100 * r['rms_m']:.2f} cm wrms {100 * r['wrms_m']:.2f} cm floor {100 * r['floor_m']:.2f} cm")
        for k in ("gamma_minus_1", "tdb_scale_minus_1", "eta_nordtvedt", "eta_sinD", "dt_0", "dt_1", "dt_2"):
            p = r["params"][k]
            print(f"   {k:18s} {p['value']:+.3e}  formal {p['formal']:.1e}  boot {p['boot']:.1e}")
        print("   no along-track:", res[name + "_no_along_track"])
    print(json.dumps({k: v for k, v in res.items() if k.endswith("coefficient_m")}))
