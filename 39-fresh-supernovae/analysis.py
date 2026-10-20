"""Every number in D6, from cached data and fits.

The steps, in order:
  1. Load the DR2 sample and our two sets of light-curve fits: SALT3 (the
     spec's model; all ZTF-internal results) and SALT2-B21 (D1's own model;
     the test of D1). Apply DR2's basic cuts (cut flow).
  2. Read D1's universe (flat LambdaCDM, H0, M) and D1's standardisation
     recipe (alpha, beta, bias correction) out of D1's own data, and check the
     one convention sncosmo and SNANA do not share (MAG_OFFSET) by refitting
     Pantheon+'s own Foundation light curves.
  3. The landing test: fresh supernovae standardised exactly as D1 would,
     against D1's curve. Mean offset, trend with redshift, scatter.
  4. Scatter before and after standardisation; the standardisation fitted on
     ZTF itself, compared with D1 and the ZTF team (Ginolin et al. 2025a).
  5. The scatter budget: peculiar velocities and redshift errors, measured.
  6. Selection (Malmquist) bias above z = 0.06: data vs a forward model.
  7. Peculiar sub-types, the 4 supernovae Pantheon+ also has, the ZTF
     team's own fits (validation), the "pocket effect" split.
  8. The live snapshot.
Results go to data/results.json; per-supernova tables to data/cache/ (the
figures and the explorer read both).
"""

from __future__ import annotations

import json
import re
import time

import numpy as np
import pandas as pd
from scipy import optimize

import config
import cosmology as cosmo
import sample
import saltfit
import standardize as st

C_KM_S = cosmo.C_KM_S
POCKET_EFFECT_MJD = 58788.0  # 2019-11-01: CCD read-out change (Lacroix et al. 2025)
#: D1's sample covered |c| < 0.28 and |x1| < 2.9 (Pantheon+ cuts |c| < 0.3,
#: |x1| < 3, sigma_x1 < 1.5); D1's recipe is only tested inside that range.
D1_RANGE_C, D1_RANGE_X1, D1_MAX_X1_ERR = 0.3, 3.0, 1.5
SIGMA_INT_GUESS = 0.13  # for weighting means; the fitted value is ~0.12-0.15


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def load_d1_universe() -> dict:
    """D1's best-fit flat LambdaCDM with the Cepheid calibration."""
    path = config.D1_DIR / "data" / "results.json"
    if path.exists():
        best = json.loads(path.read_text())["calibrated"]["best"]
        return {"Om": float(best["Om"]), "H0": float(best["H0"]), "M": float(best["M"]),
                "source": "05-supernova-cosmology/data/results.json, calibrated best fit"}
    return {**config.D1_FALLBACK, "source": "config.D1_FALLBACK (D1 not run)"}


def load_fits(df: pd.DataFrame, kind: str, n_workers: int = 6) -> pd.DataFrame:
    name = {"salt3": "salt3_fits.parquet", "salt2-b21": "salt2b21_fits.parquet"}[kind]
    fits = saltfit.fit_dr2_sample(sample.fit_tasks(df), n_workers=n_workers,
                                  cache_path=config.CACHE_DIR / name, kind=kind)
    fits = fits.set_index("name")
    fits["ok"] = fits["ok"].fillna(False).astype(bool)
    return fits


def pantheon_covariance_diag(n: int) -> np.ndarray:
    """Diagonal of Pantheon+'s STAT+SYS covariance (the honest per-SN errors;
    D1 found the table's error column is ~1.4x larger than this)."""
    cache = config.CACHE_DIR / "pantheon_cov_diag.npy"
    if cache.exists():
        return np.load(cache)
    vals = pd.read_csv(config.PANTHEON_DIR / "Pantheon+SH0ES_STAT+SYS.cov", header=None,
                       engine="c").iloc[:, 0].to_numpy(float)
    assert int(vals[0]) == n and len(vals) == n * n + 1, "unexpected covariance layout"
    diag = vals[1:].reshape(n, n).diagonal().copy()
    np.save(cache, diag)
    return diag


# ---------------------------------------------------------------------------
# The convention check: sncosmo vs SNANA on Pantheon+'s own light curves
# ---------------------------------------------------------------------------

def read_snana_lightcurve(path) -> tuple[dict, pd.DataFrame]:
    """A Pantheon+ photometry file (SNANA format, FLUXCAL on zero point 27.5 AB)."""
    text = path.read_text()
    head = {}
    for key, pat in (("snid", r"SNID:\s*(\S+)"), ("z", r"REDSHIFT_HELIO:\s*([0-9.]+)"),
                     ("mwebv", r"MWEBV:\s*([0-9.]+)")):
        m = re.search(pat, text)
        head[key] = (m.group(1) if key == "snid" else float(m.group(1))) if m else None
    rows = [l.split() for l in text.splitlines() if l.startswith("OBS:")]
    lc = pd.DataFrame({"mjd": [float(r[1]) for r in rows],
                       "band": ["ps1::" + r[2] for r in rows],
                       "flux": [float(r[4]) for r in rows],
                       "fluxerr": [float(r[5]) for r in rows]})
    lc["zp"], lc["good"], lc["flag"] = 27.5, lc["fluxerr"] > 0, 0
    return head, lc


def snana_convention_check(pantheon: pd.DataFrame) -> dict:
    """Refit Pantheon+'s own Foundation light curves with our pipeline.

    With D1's model and x1, c, t0 fixed at Pantheon+'s values, the only free
    number is x0, so -2.5 log10(x0_ours / x0_Pantheon+) isolates any
    convention difference between sncosmo and SNANA. SNANA adds the B21
    model's MAG_OFFSET (0.27 mag) to model magnitudes; sncosmo does not.
    Free fits with SALT3 give the bridge between our SALT3 parameters and
    Pantheon+'s (offsets in mB, x1, c for the same photometry).
    """
    import sncosmo
    from astropy.table import Table

    found = pantheon[pantheon["IDSURVEY"] == 150].set_index("CID")
    rows = []
    b21 = saltfit.get_model("salt2-b21")
    for path in sorted(config.FOUNDATION_LC_DIR.glob("*.txt")):
        head, lc = read_snana_lightcurve(path)
        if head["snid"] not in found.index:
            continue
        q = found.loc[head["snid"]]
        q = q.iloc[0] if isinstance(q, pd.DataFrame) else q
        z, ebv = head["z"], head["mwebv"]
        row = {"snid": head["snid"], "z": z}
        # (a) D1's model, x1, c, t0 fixed: only the amplitude is fitted.
        sel = saltfit._restframe_bands(lc, z, "salt2-b21")
        ph = (sel["mjd"] - q["PKMJD"]) / (1 + z)
        sel = sel[(ph > -10) & (ph < 40)]
        b21.set(z=z, mwebv=ebv, t0=q["PKMJD"], x1=q["x1"], c=q["c"], x0=1e-2)
        try:
            _, fm = sncosmo.fit_lc(saltfit._table(sel), b21, ["x0"], modelcov=False)
            row["d_fixed"] = float(-2.5 * np.log10(fm["x0"] / q["x0"]))
        except Exception:
            row["d_fixed"] = np.nan
        # (b) free fits with both models through the full two-pass pipeline.
        for kind in ("salt2-b21", "salt3"):
            f = saltfit.fit_lightcurve(lc, z, ebv, q["PKMJD"], name=head["snid"], kind=kind)
            tag = "b21" if kind == "salt2-b21" else "s3"
            if f["ok"]:
                row[f"dmB_{tag}"] = f["mB"] - q["mB"]
                row[f"dx1_{tag}"] = f["x1"] - q["x1"]
                row[f"dc_{tag}"] = f["c"] - q["c"]
        rows.append(row)
    d = pd.DataFrame(rows)

    def stat(col):
        v = d[col].dropna().to_numpy()
        return {"mean": float(np.mean(v)), "err": float(np.std(v, ddof=1) / np.sqrt(len(v))),
                "median": float(np.median(v)), "std": float(np.std(v, ddof=1)), "n": int(len(v))}

    out = {"mag_offset_in_SALT2_INFO": saltfit.b21_info()["MAG_OFFSET"],
           "fixed_x1_c_t0": stat("d_fixed")}
    for col in ("dmB_b21", "dx1_b21", "dc_b21", "dmB_s3", "dx1_s3", "dc_s3"):
        out[col] = stat(col)
    d.to_csv(config.CACHE_DIR / "snana_convention_check.csv", index=False)
    return out


# ---------------------------------------------------------------------------
# Per-supernova table
# ---------------------------------------------------------------------------

def cut_flow(t: pd.DataFrame) -> tuple[pd.Series, list[tuple[str, int]]]:
    """DR2's basic cuts (Rigault et al. 2025, Table 1) on our SALT3 fits."""
    steps = [
        ("released", pd.Series(True, index=t.index)),
        ("has a light curve and a converged fit", t["has_lightcurve"] & t["ok"]),
        ("good light-curve sampling", t["good_sampling"].fillna(False).astype(bool)),
        (f"x1 in [{config.CUT_X1[0]:+g}, {config.CUT_X1[1]:+g}]", t["x1"].between(*config.CUT_X1)),
        (f"c in [{config.CUT_C[0]:+g}, {config.CUT_C[1]:+g}]", t["c"].between(*config.CUT_C)),
        (f"sigma_t0 <= {config.CUT_T0_ERR:g} d", t["t0_err"] <= config.CUT_T0_ERR),
        (f"sigma_x1 <= {config.CUT_X1_ERR:g}", t["x1_err"] <= config.CUT_X1_ERR),
        (f"sigma_c <= {config.CUT_C_ERR:g}", t["c_err"] <= config.CUT_C_ERR),
        ("fit probability >= 1e-7", t["fitprob"] >= config.CUT_FITPROB),
    ]
    mask = pd.Series(True, index=t.index)
    flow = []
    for label, m in steps:
        mask &= m.fillna(False).astype(bool)
        flow.append((label, int(mask.sum())))
    return mask, flow


FIT_COLS = ["ok", "fail", "t0", "x0", "x1", "c", "t0_err", "x0_err", "x1_err", "c_err",
            "cov_x0_x1", "cov_x0_c", "cov_x1_c", "mB", "mB_err", "cov_mB_x1", "cov_mB_c",
            "chisq", "ndof", "fitprob", "n_points", "bands_fitted", "lc_n_phases",
            "lc_n_before", "lc_n_after", "lc_n_bands", "det_bands", "good_sampling", "at_bound"]


def d1_distance(f: pd.DataFrame, d1: dict, recipe: st.D1Recipe, mag_offset: float,
                host_mass, z_cmb) -> pd.DataFrame:
    """D1's recipe on fits made with D1's model (columns mB, x1, c, ...).

    mB on SNANA's scale is sncosmo's minus MAG_OFFSET. Then
        nominal: mu = mB + alpha x1 - beta c + const - M
        full:    nominal minus Pantheon+'s low-z bias correction (emulated).
    """
    mb = f["mB"] - mag_offset
    out = pd.DataFrame(index=f.index)
    out["mu_nominal"] = mb + recipe.alpha * f["x1"] - recipe.beta * f["c"] + recipe.const - d1["M"]
    out["bias"] = recipe.bias_correction(f["x1"], f["c"], host_mass, z_cmb)
    out["mu_full"] = out["mu_nominal"] - out["bias"]
    out["sigma_fit"] = st.tripp_error(f, recipe.alpha, recipe.beta)
    return out


def build_table(df, s3, b21, d1, recipe, mag_offset):
    t = df.join(s3[[k for k in FIT_COLS if k in s3.columns]])
    t["ok"] = t["ok"].fillna(False).astype(bool)
    t["basic_cuts"], flow = cut_flow(t)
    t["normal"] = ~t["peculiar"]
    t["fresh"] = ~t["in_pantheon"]
    t["cosmo"] = t["basic_cuts"] & t["normal"] & t["fresh"]
    t["mu_d1"] = cosmo.distance_modulus(t["z_cmb"], t["redshift"], d1["Om"], d1["H0"])
    t["resid_raw"] = t["mB"] - t["mu_d1"]

    # D1's model fits, and D1's recipe on them.
    bb = b21[[k for k in FIT_COLS if k in b21.columns]].add_prefix("b21_")
    t = t.join(bb)
    t["b21_ok"] = t["b21_ok"].fillna(False).astype(bool)
    f = t[["b21_mB", "b21_x1", "b21_c", "b21_mB_err", "b21_x1_err", "b21_c_err",
           "b21_cov_mB_x1", "b21_cov_mB_c", "b21_cov_x1_c"]].rename(
        columns=lambda c: c.removeprefix("b21_"))
    dd = d1_distance(f, d1, recipe, mag_offset, t["host_mass"], t["z_cmb"])
    t["mu_d1_nominal"], t["mu_d1_full"] = dd["mu_nominal"], dd["mu_full"]
    t["bias_emulated"], t["sigma_fit_d1"] = dd["bias"], dd["sigma_fit"]
    t["resid_d1_nominal"] = t["mu_d1_nominal"] - t["mu_d1"]
    t["resid_d1_full"] = t["mu_d1_full"] - t["mu_d1"]
    t["d1_test"] = (t["cosmo"] & t["b21_ok"] & (t["b21_c"].abs() < D1_RANGE_C)
                    & (t["b21_x1"].abs() < D1_RANGE_X1) & (t["b21_x1_err"] < D1_MAX_X1_ERR))
    # What one gets by (wrongly) treating sncosmo's SALT3 mB as if it were on
    # SNANA's scale, with D1's alpha and beta: the trap this project fell into.
    t["resid_naive_salt3"] = (t["mB"] + recipe.alpha * t["x1"] - recipe.beta * t["c"]
                              + recipe.const - d1["M"] - t["mu_d1"])

    # Error budget per supernova (fit errors from D1's fits for D1's test).
    t["sigma_vpec"] = cosmo.sigma_mu_from_velocity(t["z_cmb"], config.SIGMA_VPEC)
    zerr = t["redshift_err"].fillna(0.0).to_numpy(float)
    imprecise = ~t["z_group"].eq("galaxy").to_numpy()
    zerr = np.where(imprecise, np.maximum(zerr, config.SIGMA_Z_SN_FEATURES), zerr)
    t["sigma_zerr"] = cosmo.sigma_mu_from_velocity(t["z_cmb"], zerr * C_KM_S)
    t["sigma_meas_d1"] = np.sqrt(t["sigma_fit_d1"] ** 2 + t["sigma_vpec"] ** 2 + t["sigma_zerr"] ** 2)
    return t, flow


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def scatter_summary(r) -> dict:
    r = np.asarray(r, float)
    r = r[np.isfinite(r)]
    keep = st.chauvenet_mask(r)
    return {"n": int(len(r)), "median": float(np.median(r)), "mean": float(np.mean(r)),
            "std": float(np.std(r, ddof=1)), "nmad": st.nmad(r),
            "n_clipped": int((~keep).sum()), "std_clipped": float(np.std(r[keep], ddof=1)),
            "nmad_clipped": st.nmad(r[keep])}


def weighted_offset(r, sigma_meas, sigma_int) -> dict:
    """Weighted mean residual with an intrinsic scatter, and its error."""
    r = np.asarray(r, float)
    w = 1 / (np.asarray(sigma_meas, float) ** 2 + sigma_int ** 2)
    mean = float(np.sum(w * r) / np.sum(w))
    return {"mean": mean, "err": float(1 / np.sqrt(np.sum(w))), "n": int(len(r))}


def linear_trend(x, y, sigma) -> dict:
    """Weighted straight-line fit y = a + b x; returns b and its error."""
    w = 1 / np.asarray(sigma, float) ** 2
    A = np.column_stack([np.ones(len(x)), np.asarray(x, float)])
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    coef = cov @ (A.T @ (w * np.asarray(y, float)))
    return {"intercept": float(coef[0]), "slope": float(coef[1]),
            "slope_err": float(np.sqrt(cov[1, 1])),
            "slope_sigma": float(coef[1] / np.sqrt(cov[1, 1]))}


def binned(x, y, edges, sigma=None) -> list[dict]:
    """Mean (weighted if sigma given), error and nMAD of y in bins of x."""
    out = []
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (x >= lo) & (x < hi) & np.isfinite(y)
        n = int(m.sum())
        if n < 3:
            continue
        ym = y[m]
        if sigma is not None:
            w = 1 / np.asarray(sigma, float)[m] ** 2
            mean = float(np.sum(w * ym) / np.sum(w))
            err = float(1 / np.sqrt(np.sum(w)))
        else:
            mean = float(np.mean(ym))
            err = float(np.std(ym, ddof=1) / np.sqrt(n))
        out.append({"lo": float(lo), "hi": float(hi), "x": float(np.median(x[m])), "n": n,
                    "mean": mean, "err": err, "median": float(np.median(ym)),
                    "nmad": st.nmad(ym), "std": float(np.std(ym, ddof=1))})
    return out


def fit_velocity_dispersion(r, sigma_fit, z, outlier_sigma: float = 1.0) -> dict:
    """Intrinsic scatter and random peculiar-velocity dispersion, robustly.

    var_i = sigma_fit_i^2 + sigma_SN^2 + sigma_mu(z_i; sigma_v)^2. The two
    terms are separable because the velocity term grows as 1/z. Hubble
    residuals have heavier tails than a Gaussian (std ~1.3 x nMAD: odd
    supernovae, wrong host galaxies), and a pure Gaussian likelihood lets a
    handful of those masquerade as a large sigma_v. So each supernova is
    either "good" (the Gaussian above) or, with probability f, an outlier
    drawn from a broad Gaussian (sigma = 1 mag); f is fitted too.
    """
    r = np.asarray(r, float)
    unit = cosmo.sigma_mu_from_velocity(np.asarray(z, float), 1.0)  # mag per km/s
    sf2 = np.asarray(sigma_fit, float) ** 2

    def nll(p, robust=True):
        s_sn, s_v, mu, logit_f = p
        var = sf2 + s_sn ** 2 + (unit * s_v) ** 2
        good = np.exp(-0.5 * (r - mu) ** 2 / var) / np.sqrt(2 * np.pi * var)
        if not robust:
            return -np.sum(np.log(good + 1e-300))
        f = 1 / (1 + np.exp(-logit_f))
        bad = np.exp(-0.5 * (r - mu) ** 2 / outlier_sigma ** 2) / np.sqrt(2 * np.pi * outlier_sigma ** 2)
        return -np.sum(np.log((1 - f) * good + f * bad))

    out = {"n": int(len(r))}
    for robust in (False, True):
        x0 = [0.12, 300.0, float(np.median(r)), -3.0]
        res = optimize.minimize(lambda p: nll(p, robust), x0, method="Nelder-Mead",
                                options={"xatol": 1e-6, "fatol": 1e-8, "maxiter": 20000})
        grid = np.linspace(0, 900, 91)
        prof = np.array([optimize.minimize(lambda q: nll([q[0], sv, q[1], q[2]], robust),
                                           [res.x[0], res.x[2], res.x[3]], method="Nelder-Mead",
                                           options={"xatol": 1e-6, "fatol": 1e-8}).fun for sv in grid])
        prof = 2 * (prof - prof.min())
        inside = grid[prof <= 1.0]
        key = "robust" if robust else "gaussian"
        out[key] = {"sigma_sn": float(abs(res.x[0])), "sigma_v": float(abs(res.x[1])),
                    "sigma_v_68": [float(inside.min()), float(inside.max())],
                    "outlier_fraction": float(1 / (1 + np.exp(-res.x[3]))) if robust else 0.0,
                    "profile_grid": grid.tolist(), "profile_dchi2": prof.tolist()}
    # The robust fit is the one quoted.
    out.update({k: out["robust"][k] for k in ("sigma_sn", "sigma_v", "sigma_v_68", "outlier_fraction",
                                              "profile_grid", "profile_dchi2")})
    return out


def fit_redshift_scatter(r, sigma_fit, z, sigma_sn, sigma_v, outlier_fraction=0.0,
                         outlier_sigma: float = 1.0) -> dict:
    """Extra redshift scatter (in z) needed on top of a known sigma_SN, sigma_v
    (same outlier mixture as fit_velocity_dispersion, fraction fitted)."""
    r = np.asarray(r, float)
    unit = cosmo.sigma_mu_from_velocity(np.asarray(z, float), 1.0)
    sf2 = np.asarray(sigma_fit, float) ** 2

    def nll(p):
        sz, mu, logit_f = p
        var = sf2 + sigma_sn ** 2 + (unit * sigma_v) ** 2 + (unit * sz * C_KM_S) ** 2
        good = np.exp(-0.5 * (r - mu) ** 2 / var) / np.sqrt(2 * np.pi * var)
        f = 1 / (1 + np.exp(-logit_f))
        bad = np.exp(-0.5 * (r - mu) ** 2 / outlier_sigma ** 2) / np.sqrt(2 * np.pi * outlier_sigma ** 2)
        return -np.sum(np.log((1 - f) * good + f * bad))

    res = optimize.minimize(nll, [3e-3, float(np.median(r)), -3.0], method="Nelder-Mead",
                            options={"xatol": 1e-8, "fatol": 1e-9})
    grid = np.linspace(0, 8e-3, 81)
    prof = np.array([optimize.minimize(lambda q: nll([g, q[0], q[1]]), [res.x[1], res.x[2]],
                                       method="Nelder-Mead").fun for g in grid])
    prof = 2 * (prof - prof.min())
    inside = grid[prof <= 1.0]
    return {"sigma_z": float(abs(res.x[0])), "sigma_z_68": [float(inside.min()), float(inside.max())],
            "sigma_v_equiv": float(abs(res.x[0]) * C_KM_S), "n": int(len(r)),
            "outlier_fraction": float(1 / (1 + np.exp(-res.x[2])))}


# ---------------------------------------------------------------------------
# Coherent flows: what the low-redshift residuals look like on the sky
# ---------------------------------------------------------------------------

def _unit_vectors(ra, dec):
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    g = SkyCoord(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg).galactic
    return np.column_stack([np.cos(g.b.rad) * np.cos(g.l.rad), np.cos(g.b.rad) * np.sin(g.l.rad),
                            np.sin(g.b.rad)])


def fit_flows(d: pd.DataFrame, resid_col: str, sigma_fit_col: str,
              sigma_sn: float = 0.15, sigma_v: float = 300.0, seed: int = 9) -> dict:
    """Peculiar-velocity models for low-z residuals.

    A galaxy moving away from us faster than the Hubble flow (radial
    peculiar velocity v > 0) has a redshift too high for its distance, so its
    supernova looks closer than predicted: residual = -k(z) v, with
    k = 5/ln10 (1+z) / (c z (1+z/2)) magnitudes per km/s. Models:
      offset only; + monopole (everything moving in or out together, or
      anything else that grows as 1/z, like a brightness-dependent
      photometry bias); + dipole (a bulk flow towards one direction); + both.
    """
    zc0 = d["z_cmb"].to_numpy()
    k0 = cosmo.sigma_mu_from_velocity(zc0, 1.0)
    var0 = d[sigma_fit_col].to_numpy() ** 2 + sigma_sn ** 2 + (k0 * sigma_v) ** 2
    r0 = d[resid_col].to_numpy()
    # Drop outliers first (Chauvenet on the normalised residuals).
    keep = st.chauvenet_mask((r0 - np.median(r0)) / np.sqrt(var0))
    d = d[keep]
    n = _unit_vectors(d["ra"], d["dec"])
    zc, k, var, r = zc0[keep], k0[keep], var0[keep], r0[keep]

    def chi2(p, mono, dip):
        off, vm, vx, vy, vz = p
        v = (vm if mono else 0.0) + ((n @ np.array([vx, vy, vz])) if dip else 0.0)
        return float(np.sum((r - off + k * v) ** 2 / var))

    out = {"n": int(len(r)), "n_clipped": int((~keep).sum()), "z_range": [float(zc.min()), float(zc.max())],
           "assumed_sigma_sn": sigma_sn, "assumed_sigma_v": sigma_v}
    base = None
    for label, mono, dip in (("offset", False, False), ("monopole", True, False),
                             ("dipole", False, True), ("monopole+dipole", True, True)):
        free = [0] + ([1] if mono else []) + ([2, 3, 4] if dip else [])

        def f(q):
            p = np.zeros(5)
            p[free] = q
            return chi2(p, mono, dip)

        q = optimize.minimize(f, np.zeros(len(free)), method="Powell",
                              options={"maxiter": 40000, "xtol": 1e-7, "ftol": 1e-10}).x
        q = optimize.minimize(f, q, method="BFGS").x
        H = st._hessian(f, q, rel=1e-3)
        cov = 2 * np.linalg.inv(H)
        p = np.zeros(5)
        p[free] = q
        res = {"chi2": f(q), "n_params": len(free), "offset": float(p[0]),
               "offset_err": float(np.sqrt(cov[0, 0]))}
        if mono:
            i = free.index(1)
            res["monopole_kms"], res["monopole_err"] = float(p[1]), float(np.sqrt(cov[i, i]))
        if dip:
            idx = [free.index(j) for j in (2, 3, 4)]
            V = p[2:]
            cv = cov[np.ix_(idx, idx)]
            draws = np.random.default_rng(seed).multivariate_normal(V, cv, 4000)
            amp = np.linalg.norm(draws, axis=1)
            ls = np.degrees(np.arctan2(draws[:, 1], draws[:, 0])) % 360
            bs = np.degrees(np.arcsin(draws[:, 2] / amp))
            l0 = float(np.degrees(np.arctan2(V[1], V[0])) % 360)
            dl = (ls - l0 + 180) % 360 - 180
            res.update({"dipole_kms": float(np.linalg.norm(V)), "dipole_err": float(np.std(amp)),
                        "l": l0, "l_err": float(np.std(dl)),
                        "b": float(np.degrees(np.arcsin(V[2] / np.linalg.norm(V)))),
                        "b_err": float(np.std(bs)), "vector": V.tolist()})
        if base is None:
            base = res["chi2"]
        res["delta_chi2_vs_offset"] = float(base - res["chi2"])
        out[label] = res
    return out


# ---------------------------------------------------------------------------
# Selection bias: a forward model (a labelled simulation built from real data)
# ---------------------------------------------------------------------------

def peak_magnitudes(t: pd.DataFrame, band: str = "ztfr") -> np.ndarray:
    """Observed peak magnitude in a ZTF band from each SN's fitted SALT3 model."""
    model = saltfit.get_model("salt3")
    out = np.full(len(t), np.nan)
    for i, (_, r) in enumerate(t.iterrows()):
        if not r["ok"] or not np.isfinite(r["x0"]):
            continue
        model.set(z=r["redshift"], t0=r["t0"], x0=r["x0"], x1=r["x1"], c=r["c"], mwebv=r["mwebv"])
        out[i] = model.bandmag(band, "ab", r["t0"])
    return out


def selection_model(vl: pd.DataFrame, fit: st.TrippFit, z_grid, d1: dict,
                    n_draw: int = 40000, seed: int = 6) -> dict:
    """Expected Hubble-residual bias and mean (x1, c) of selected SNe vs z.

    Real ingredients: (x1, c) drawn from our own volume-limited sample (which
    represents nature at z < 0.06), the intrinsic scatter and alpha, beta, M
    fitted on ZTF, the observed relation between SALT's mB and the ZTF r-band
    peak, and the ZTF DR2 selection function P(classified | peak mag)
    = 1 / (1 + exp(4.5 (m - 18.8))) (Rigault et al. 2025, Fig. 4).
    Simulated: which of them a magnitude-limited survey would have classified.
    """
    rng = np.random.default_rng(seed)
    vl = vl[np.isfinite(vl["m_peak_r"])]
    A = np.column_stack([np.ones(len(vl)), vl["redshift"], vl["c"], vl["x1"]])
    coef, *_ = np.linalg.lstsq(A, vl["m_peak_r"] - vl["mB"], rcond=None)
    a, b, M = fit.params["alpha"], fit.params["beta"], fit.params["M"]
    m50, k = config.DR2_PUBLISHED["selection_m50"], config.DR2_PUBLISHED["selection_k"]
    idx = rng.integers(0, len(vl), n_draw)
    x1 = vl["x1"].to_numpy()[idx]
    c = vl["c"].to_numpy()[idx]
    eps = rng.normal(0, fit.sigma_int, n_draw)
    out = []
    for z in z_grid:
        # M was fitted relative to D1's distances, so use D1's universe here.
        mu = float(cosmo.distance_modulus(np.array([z]), np.array([z]), d1["Om"], d1["H0"])[0])
        mB = mu + M - a * x1 + b * c + eps
        m_r = mB + coef[0] + coef[1] * z + coef[2] * c + coef[3] * x1
        p = 1 / (1 + np.exp(k * (m_r - m50)))
        out.append({"z": float(z), "bias": float(np.sum(p * eps) / np.sum(p)),
                    "mean_c": float(np.sum(p * c) / np.sum(p)),
                    "mean_x1": float(np.sum(p * x1) / np.sum(p)),
                    "fraction_selected": float(np.mean(p))})
    return {"grid": out, "peak_relation": dict(zip(["const", "z", "c", "x1"], map(float, coef))),
            "n_parent": int(len(vl)), "sigma_int": float(fit.sigma_int)}


# ---------------------------------------------------------------------------
# The whole analysis
# ---------------------------------------------------------------------------

def run(n_workers: int = 6) -> dict:
    t_start = time.time()
    R: dict = {}

    df = sample.load_dr2()
    s3 = load_fits(df, "salt3", n_workers)
    b21 = load_fits(df, "salt2-b21", n_workers)
    pantheon = sample.load_pantheon()
    d1 = load_d1_universe()
    recipe = st.recover_d1_recipe(pantheon)
    conv = snana_convention_check(pantheon)
    mag_offset = saltfit.b21_info()["MAG_OFFSET"]
    R["d1_universe"] = d1
    R["d1_recipe"] = {"alpha": recipe.alpha, "beta": recipe.beta, "const": recipe.const,
                      "regression_rms": recipe.regression_rms, "bias_coef": recipe.bias_coef,
                      "bias_rms": recipe.bias_rms, "bias_n": recipe.bias_n,
                      "mag_offset": mag_offset, "range": {"abs_c": D1_RANGE_C, "abs_x1": D1_RANGE_X1},
                      "published": config.PANTHEON_STANDARDISATION}
    R["convention_check"] = conv

    t, flow = build_table(df, s3, b21, d1, recipe, mag_offset)
    R["cut_flow"] = flow
    theirs = (t["t21_lccoverage_flag"] == 1) & (t["t21_fitquality_flag"] == 1)
    R["sample"] = {
        "n_released": int(len(t)), "n_lightcurves": int(t["has_lightcurve"].sum()),
        "n_in_pantheon": int(t["in_pantheon"].sum()),
        "in_pantheon": t.index[t["in_pantheon"]].tolist(),
        "n_basic_cuts": int(t["basic_cuts"].sum()),
        "n_basic_normal": int((t["basic_cuts"] & t["normal"]).sum()),
        "n_cosmo": int(t["cosmo"].sum()),
        "n_cosmo_volume_limited": int((t["cosmo"] & (t["redshift"] <= config.Z_VOLUME_LIMIT)).sum()),
        "n_d1_test": int(t["d1_test"].sum()),
        "n_b21_ok": int(t["b21_ok"].sum()),
        "fit_failures": t.loc[t["has_lightcurve"] & ~t["ok"], "fail"]
        .astype(str).str.split(":").str[0].value_counts().to_dict(),
        "z_groups_cosmo": t.loc[t["cosmo"], "z_group"].value_counts().to_dict(),
        "n_dr2_flags_pass": int(theirs.sum()),
        "agreement_with_dr2_flags": {
            "both": int((theirs & t["basic_cuts"]).sum()),
            "ours_only": int((~theirs & t["basic_cuts"]).sum()),
            "theirs_only": int((theirs & ~t["basic_cuts"]).sum())},
        "published": config.DR2_PUBLISHED,
    }

    z = t["redshift"]
    windows = {
        "all": pd.Series(True, index=t.index),
        "dr2_clean_0.03-0.1": z.between(*config.Z_DR2_CLEAN),
        "volume_limited_z<=0.06": z <= config.Z_VOLUME_LIMIT,
        "hubble_flow_0.023-0.06": z.between(config.Z_HUBBLE_FLOW_MIN, config.Z_VOLUME_LIMIT),
    }

    # --- The landing test: D1's model, D1's recipe, D1's universe ----------
    d1s = t["d1_test"]
    hf = d1s & windows["hubble_flow_0.023-0.06"]
    vl = d1s & windows["volume_limited_z<=0.06"] & (z >= 0.015)
    tot = np.sqrt(t["sigma_meas_d1"] ** 2 + SIGMA_INT_GUESS ** 2)
    land = {}
    for flavour in ("nominal", "full"):
        r = t[f"resid_d1_{flavour}"]
        off = weighted_offset(r[hf], t.loc[hf, "sigma_meas_d1"], SIGMA_INT_GUESS)
        hm = t["host_mass"]
        hi, lo = vl & (hm > 10), vl & (hm <= 10)
        land[flavour] = {
            "offset_hubble_flow": off,
            # If the offset were real, D1's H0 would have to change by this much
            # (fainter supernovae = farther = slower expansion). It is not real:
            # ZTF DR2's flux scale is blinded.
            "H0_that_would_absorb_it": float(d1["H0"] * 10 ** (-off["mean"] / 5)),
            "trend_vs_z_0.015-0.06": linear_trend(z[vl], r[vl], tot[vl]),
            "trend_vs_c_volume_limited": linear_trend(t.loc[vl, "b21_c"], r[vl] - off["mean"], tot[vl]),
            "trend_vs_x1_volume_limited": linear_trend(t.loc[vl, "b21_x1"], r[vl] - off["mean"], tot[vl]),
            "mass_step": {
                "high_minus_low": float(np.median(r[hi]) - np.median(r[lo])),
                "err": float(1.2533 * np.hypot(st.nmad(r[hi]) / np.sqrt(hi.sum()),
                                                st.nmad(r[lo]) / np.sqrt(lo.sum()))),
                "n_high": int(hi.sum()), "n_low": int(lo.sum())},
            "scatter_hubble_flow": scatter_summary(r[hf]),
            "scatter_dr2_clean": scatter_summary(r[d1s & windows["dr2_clean_0.03-0.1"]]),
        }
    naive = t["resid_naive_salt3"]
    land["naive_salt3_offset"] = weighted_offset(
        naive[t["cosmo"] & windows["hubble_flow_0.023-0.06"]],
        t.loc[t["cosmo"] & windows["hubble_flow_0.023-0.06"], "sigma_vpec"], SIGMA_INT_GUESS)
    land["sigma_int_assumed"] = SIGMA_INT_GUESS
    land["lacroix_2025_dr2_vs_smp_offset"] = config.LACROIX_2025_OFFSET_MAG
    R["landing"] = land

    # Pantheon+ itself against D1's curve (the reference D1 was fitted to).
    pz = pantheon
    p_mu_d1 = cosmo.distance_modulus(pz["zHD"], pz["zHEL"], d1["Om"], d1["H0"])
    p_resid = pz["m_b_corr"] - d1["M"] - p_mu_d1
    p_err = np.sqrt(pantheon_covariance_diag(len(pz)))
    pm = (pz["zHD"] > 0.01) & (pz["IS_CALIBRATOR"] == 0)
    p_hf = pm & pz["zHD"].between(config.Z_HUBBLE_FLOW_MIN, config.Z_VOLUME_LIMIT)
    R["pantheon_reference"] = {
        "offset_hubble_flow": weighted_offset(p_resid[p_hf], p_err[p_hf], 0.0),
        "scatter_hubble_flow": scatter_summary(p_resid[p_hf]),
        "scatter_dr2_clean": scatter_summary(p_resid[pm & pz["zHD"].between(*config.Z_DR2_CLEAN)]),
        "n_hubble_flow": int(p_hf.sum()),
    }

    # --- Standardisation fitted on ZTF itself (SALT3) -------------------------
    cos = t["cosmo"]
    fit_sel = cos & windows["volume_limited_z<=0.06"] & np.isfinite(t["host_mass"]) & (z >= 0.015)
    sub = t[fit_sel]
    extra = sub["sigma_vpec"] ** 2 + sub["sigma_zerr"] ** 2
    fits = {
        "linear_all": st.fit_tripp(sub, sub["mu_d1"], extra),
        "with_mass_step_all": st.fit_tripp(sub, sub["mu_d1"], extra, host_mass=sub["host_mass"]),
        "linear": st.fit_tripp(sub, sub["mu_d1"], extra, clip=3.5),
        "with_mass_step": st.fit_tripp(sub, sub["mu_d1"], extra, host_mass=sub["host_mass"], clip=3.5),
    }
    fit_mass = fits["with_mass_step"]
    R["ztf_standardisation"] = {
        "sample": "DR2 cuts, normal, new to D1, 0.015 <= z <= 0.06, host mass known (SALT3); "
                  "'linear'/'with_mass_step' drop > 3.5 sigma outliers, '_all' keep them",
        **{k: {"params": f.params, "errors": f.errors, "sigma_int": f.sigma_int, "n": f.n,
               "n_clipped": f.n_clipped} for k, f in fits.items()},
        "ginolin_2025a": config.GINOLIN_2025A,
    }
    a, b, M, g = (fit_mass.params[k] for k in ("alpha", "beta", "M", "gamma"))
    step = np.where(t["host_mass"] > 10, 0.5, np.where(t["host_mass"] <= 10, -0.5, 0.0))
    t["resid_ztf"] = t["mB"] + a * t["x1"] - b * t["c"] - M + g * step - t["mu_d1"]
    t["sigma_fit_ztf"] = st.tripp_error(t, a, b)
    t["sigma_tot_ztf"] = np.sqrt(t["sigma_fit_ztf"] ** 2 + t["sigma_vpec"] ** 2
                                 + t["sigma_zerr"] ** 2 + fit_mass.sigma_int ** 2)

    # --- Scatter before and after standardisation ----------------------------
    raw = t["resid_raw"] - t.loc[cos, "resid_raw"].median()
    R["scatter"] = {}
    for wname, m in windows.items():
        R["scatter"][wname] = {
            "raw_salt3": scatter_summary(raw[cos & m]),
            "ztf_fit_salt3": scatter_summary(t.loc[cos & m, "resid_ztf"]),
            "d1_recipe_b21": scatter_summary(t.loc[d1s & m, "resid_d1_full"]),
        }
    for grp in ("galaxy", "SN features"):
        m = cos & t["z_group"].eq(grp) & windows["dr2_clean_0.03-0.1"]
        R["scatter"][f"dr2_clean_0.03-0.1, {grp} redshifts"] = {
            "ztf_fit_salt3": scatter_summary(t.loc[m, "resid_ztf"])}

    # --- Scatter budget: velocities and redshift errors ----------------------
    gal = cos & t["z_group"].eq("galaxy") & (z <= config.Z_VOLUME_LIMIT) & (z >= 0.005)
    budget = fit_velocity_dispersion(t.loc[gal, "resid_ztf"], t.loc[gal, "sigma_fit_ztf"],
                                     t.loc[gal, "z_cmb"])
    snf = cos & t["z_group"].eq("SN features") & (z <= config.Z_VOLUME_LIMIT)
    zs = fit_redshift_scatter(t.loc[snf, "resid_ztf"], t.loc[snf, "sigma_fit_ztf"],
                              t.loc[snf, "z_cmb"], budget["sigma_sn"], budget["sigma_v"])
    # Coherent flows in the low-z residuals (D1's recipe, host-galaxy redshifts).
    lowz = t["d1_test"] & t["z_group"].eq("galaxy") & z.between(0.01, config.Z_VOLUME_LIMIT)
    flows = fit_flows(t[lowz], "resid_d1_full", "sigma_fit_d1")
    flows["carrick_2015_external_bulk_flow"] = {"v": 159, "v_err": 23, "l": 304, "l_err": 11,
                                                "b": 6, "b_err": 13}
    # The random-velocity dispersion again, once the coherent flow is removed.
    best = flows["monopole+dipole"]
    nvec = _unit_vectors(t.loc[gal, "ra"], t.loc[gal, "dec"])
    v_coh = best["monopole_kms"] + nvec @ np.array(best["vector"])
    r_flow = t.loc[gal, "resid_ztf"] + cosmo.sigma_mu_from_velocity(t.loc[gal, "z_cmb"], 1.0) * v_coh
    budget_after = fit_velocity_dispersion(r_flow, t.loc[gal, "sigma_fit_ztf"], t.loc[gal, "z_cmb"])
    flows["sigma_v_after_removing_flow"] = {k: budget_after[k] for k in
                                            ("sigma_sn", "sigma_v", "sigma_v_68", "outlier_fraction", "n")}
    R["flows"] = flows
    edges_b = np.array([0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06])
    R["scatter_budget"] = {
        "galaxy_redshifts": budget, "sn_feature_redshifts": zs,
        "binned_galaxy": binned(t.loc[gal, "z_cmb"], t.loc[gal, "resid_ztf"], edges_b),
        "binned_sn_features": binned(t.loc[snf, "z_cmb"], t.loc[snf, "resid_ztf"], edges_b),
        "median_sigma_fit": float(t.loc[gal, "sigma_fit_ztf"].median()),
        "dr2_quoted": {"sigma_v_typical": 300, "sigma_z_sn": 3e-3},
    }

    # --- Selection bias above the volume limit -------------------------------
    t["m_peak_r"] = peak_magnitudes(t)
    z_grid = np.round(np.arange(0.01, 0.2001, 0.005), 4)
    sel = selection_model(t[cos & windows["volume_limited_z<=0.06"]], fit_mass, z_grid, d1)
    edges = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
                      0.11, 0.12, 0.14, 0.16, 0.20])
    R["selection"] = {
        "model": sel,
        "data_resid_ztf": binned(z[cos], t.loc[cos, "resid_ztf"], edges, t.loc[cos, "sigma_tot_ztf"]),
        "data_mean_c": binned(z[cos], t.loc[cos, "c"], edges),
        "data_mean_x1": binned(z[cos], t.loc[cos, "x1"], edges),
    }

    # --- Peculiar sub-types (SALT3, ZTF standardisation) ----------------------
    usable = t["ok"] & t["good_sampling"].fillna(False).astype(bool) & \
        (t["c_err"] < 0.15) & (t["x1_err"] < 1.5) & t["fresh"]
    t["subtype"] = t["sub_type"].fillna("unspecified")
    zero = float(np.median(t.loc[cos, "resid_ztf"]))
    pec = {}
    for name, grp in t[usable].groupby("subtype"):
        r = grp["resid_ztf"] - zero
        pec[name] = {"n": int(len(grp)), "n_peculiar_label": int(grp["peculiar"].sum()),
                     "median_resid": float(np.median(r)), "nmad": st.nmad(r),
                     "median_x1": float(np.median(grp["x1"])), "median_c": float(np.median(grp["c"])),
                     "pass_basic_cuts": int(grp["basic_cuts"].sum())}
    R["subtypes"] = pec

    # --- The 4 supernovae Pantheon+ also has ----------------------------------
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    cp = SkyCoord(pantheon["RA"].to_numpy() * u.deg, pantheon["DEC"].to_numpy() * u.deg)
    over = []
    for name, r in t[t["in_pantheon"]].iterrows():
        sep = SkyCoord(r["ra"] * u.deg, r["dec"] * u.deg).separation(cp).arcsec
        for _, q in pantheon[sep < 5].iterrows():
            over.append({
                "ztf": name, "iau": r["iau_name"], "pantheon_cid": q["CID"],
                "pantheon_survey": int(q["IDSURVEY"]), "calibrator": bool(q["IS_CALIBRATOR"]),
                "z": float(r["redshift"]),
                "mB_ztf_d1model": float(r["b21_mB"] - mag_offset) if r["b21_ok"] else None,
                "mB_pantheon": float(q["mB"]),
                "x1_ztf_d1model": float(r["b21_x1"]) if r["b21_ok"] else None, "x1_pantheon": float(q["x1"]),
                "c_ztf_d1model": float(r["b21_c"]) if r["b21_ok"] else None, "c_pantheon": float(q["c"]),
                "mu_ztf_d1": float(r["mu_d1_full"]) if r["b21_ok"] else None,
                "mu_pantheon": float(q["m_b_corr"] - d1["M"]),
                "mu_cepheid": float(q["CEPH_DIST"]) if q["IS_CALIBRATOR"] else None,
                "salt3_fit_ok": bool(r["ok"]), "fail": None if r["ok"] else str(r["fail"]),
            })
    R["overlap_with_pantheon"] = over

    # --- Validation against the ZTF team's own fits ---------------------------
    both = t["ok"] & np.isfinite(t["s3zt_x1"]) & t["basic_cuts"]
    val = {"n": int(both.sum())}
    for k in ("x1", "c", "mB", "t0"):
        d = (t.loc[both, k] - t.loc[both, "s3zt_" + k]).to_numpy()
        err = t.loc[both, k + "_err"].to_numpy()
        val[k] = {"median_diff": float(np.median(d)), "nmad_diff": st.nmad(d),
                  "p95_abs": float(np.percentile(np.abs(d), 95)), "nmad_pull": st.nmad(d / err),
                  "median_err_ours": float(np.median(err)),
                  "median_err_theirs": float(np.median(t.loc[both, "s3zt_" + k + "_err"]))
                  if k != "mB" else None}
    val["vs_their_salt2_t21"] = {
        k: {"median_diff": float(np.median(t.loc[both, k] - t.loc[both, "t21_" + k])),
            "nmad_diff": st.nmad(t.loc[both, k] - t.loc[both, "t21_" + k])}
        for k in ("x1", "c", "mB")}
    bb = both & t["b21_ok"]
    val["salt3_vs_d1_model_on_ztf"] = {
        k: {"median_diff": float(np.median(t.loc[bb, k] - t.loc[bb, "b21_" + k])),
            "nmad_diff": st.nmad(t.loc[bb, k] - t.loc[bb, "b21_" + k])}
        for k in ("x1", "c", "mB")}
    R["validation"] = val

    # --- Pocket effect: before vs after November 2019 -------------------------
    vlx = cos & windows["volume_limited_z<=0.06"]
    early, late = vlx & (t["t0"] < POCKET_EFFECT_MJD), vlx & (t["t0"] >= POCKET_EFFECT_MJD)
    R["pocket_effect_split"] = {
        "before_nov2019": weighted_offset(t.loc[early, "resid_ztf"], t.loc[early, "sigma_tot_ztf"], 0.0),
        "after_nov2019": weighted_offset(t.loc[late, "resid_ztf"], t.loc[late, "sigma_tot_ztf"], 0.0),
    }

    # --- Binned Hubble residuals for the main figure --------------------------
    zedges = np.array([0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08,
                       0.1, 0.12, 0.15, 0.2, 0.3])
    R["binned"] = {
        "ztf_d1": binned(z[d1s], t.loc[d1s, "resid_d1_full"], zedges, tot[d1s]),
        "pantheon": binned(pz.loc[pm, "zHD"], p_resid[pm], zedges, p_err[pm]),
    }

    # --- Live snapshot ---------------------------------------------------------
    R["live"] = run_live(d1, recipe, mag_offset, R)

    t.to_parquet(config.CACHE_DIR / "hubble_table.parquet")
    pantheon.assign(mu=pz["m_b_corr"] - d1["M"], mu_d1=p_mu_d1, resid=p_resid,
                    err=p_err).to_parquet(config.CACHE_DIR / "pantheon_table.parquet")
    R["runtime_s"] = round(time.time() - t_start, 1)
    config.RESULTS_PATH.write_text(json.dumps(_jsonable(R), indent=1))
    return R


def run_live(d1: dict, recipe: st.D1Recipe, mag_offset: float, R: dict) -> dict:
    """Fit the newest supernovae like the DR2 sample and test them against D1."""
    import live

    meta, lt, lcs = live.fit_snapshot()
    if lt.empty:
        return {"available": False}
    lt["ok"] = lt["ok"].fillna(False).astype(bool)
    # D1's model too, for D1's recipe.
    b21_rows = []
    for name, lc in lcs.items():
        r = lt.loc[name]
        f = saltfit.fit_lightcurve(lc, r["redshift"], r["mwebv"], r.get("t0_guess", np.nan),
                                   name=name, kind="salt2-b21")
        b21_rows.append({"tns_name": name, **{f"b21_{k}": f.get(k) for k in
                                               ("ok", "mB", "x1", "c", "mB_err", "x1_err", "c_err",
                                                "cov_mB_x1", "cov_mB_c", "cov_x1_c")}})
    lt = lt.join(pd.DataFrame(b21_rows).set_index("tns_name"))
    lt["b21_ok"] = lt["b21_ok"].fillna(False).astype(bool)
    lt["mu_d1"] = cosmo.distance_modulus(lt["z_cmb"], lt["redshift"], d1["Om"], d1["H0"])
    f = lt[["b21_mB", "b21_x1", "b21_c", "b21_mB_err", "b21_x1_err", "b21_c_err",
            "b21_cov_mB_x1", "b21_cov_mB_c", "b21_cov_x1_c"]].rename(
        columns=lambda c: c.removeprefix("b21_")).astype(float)
    dd = d1_distance(f, d1, recipe, mag_offset, np.full(len(lt), np.nan), lt["z_cmb"])
    lt["mu_d1_full"] = dd["mu_full"]
    lt["resid_d1_full"] = lt["mu_d1_full"] - lt["mu_d1"]
    zerr = lt["redshift_err"].to_numpy(float)
    lt["sigma_meas_d1"] = np.sqrt(dd["sigma_fit"] ** 2
                                  + cosmo.sigma_mu_from_velocity(lt["z_cmb"], config.SIGMA_VPEC) ** 2
                                  + cosmo.sigma_mu_from_velocity(lt["z_cmb"], zerr * C_KM_S) ** 2)
    cuts = (lt["ok"] & lt["good_sampling"].fillna(False).astype(bool)
            & lt["x1"].between(*config.CUT_X1) & lt["c"].between(*config.CUT_C)
            & (lt["t0_err"] <= config.CUT_T0_ERR) & (lt["x1_err"] <= config.CUT_X1_ERR)
            & (lt["c_err"] <= config.CUT_C_ERR) & (lt["fitprob"] >= config.CUT_FITPROB))
    lt["basic_cuts"] = cuts
    lt["d1_test"] = cuts & lt["b21_ok"] & (lt["b21_c"].abs() < D1_RANGE_C) \
        & (lt["b21_x1"].abs() < D1_RANGE_X1) & (lt["b21_x1_err"] < D1_MAX_X1_ERR)
    lt["days_since_peak"] = (lt["last_mjd"] - lt["t0"]) / (1 + lt["redshift"])
    lt.to_parquet(config.CACHE_DIR / "live_table.parquet")
    if lcs:
        pd.concat([lc.assign(tns_name=n) for n, lc in lcs.items()]).to_parquet(
            config.CACHE_DIR / "live_lightcurves.parquet")
    good = lt[lt["d1_test"]]
    out = {"available": True, "meta": meta, "n_objects": int(len(lt)), "n_fit_ok": int(lt["ok"].sum()),
           "n_basic_cuts": int(cuts.sum()), "n_d1_test": int(len(good)),
           "z_groups": lt["z_group"].value_counts().to_dict()}
    if len(good) >= 3:
        out["offset"] = weighted_offset(good["resid_d1_full"], good["sigma_meas_d1"], SIGMA_INT_GUESS)
        out["scatter"] = scatter_summary(good["resid_d1_full"])
        out["dr2_offset_same_recipe"] = R["landing"]["full"]["offset_hubble_flow"]["mean"]
    cols = ["tns_name", "oid", "redshift", "z_group", "discovery_date", "t0", "x1", "x1_err", "c",
            "c_err", "resid_d1_full", "sigma_meas_d1", "basic_cuts", "d1_test", "good_sampling",
            "days_since_peak", "fail", "classifying_groups"]
    out["objects"] = _jsonable(lt[cols].to_dict(orient="records"))
    return out


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if o is pd.NA:
        return None
    return o


if __name__ == "__main__":
    res = run()
    print(json.dumps({k: res[k] for k in ("sample", "landing")}, indent=1)[:5000])
