"""Every number in the README and the data behind every figure, from the cached data.

Steps, each cached under ``data/derived/`` (``--force`` recomputes):

1. ``sn``      Pantheon+ alone on D1's (Omega_m, Omega_Lambda) grid
2. ``bao``     DESI BAO alone, curvature free (checked against DESI's own chain)
3. ``cmb``     the CMB's band in the plane (compressed Planck PR4; D3's own fit)
4. ``combos``  every combination of the three, on fine (Omega_m, Omega_k) grids,
               with H0 and the age of the universe derived at each point
5. ``mcmc``    an independent MCMC of BAO + CMB (checks the grid + Laplace method)
6. ``ours``    our universe: its parameters, past, future and equality epochs
7. ``gallery`` the named universes: fates and a(t) histories
8. ``map``     the plane as a map of fates, ages and bounce sizes
9. ``ages``    ages across the plane for two values of H0; the oldest stars
10. ``checks`` closed forms, CAMB, DESI's best fits

    python analysis.py           # use caches where present
    python analysis.py --force   # recompute everything (~5-10 min)
"""

from __future__ import annotations

import argparse
import json
import time

import numpy as np
from scipy import optimize

import background as bg
import config as cfg
import fetch_data as fd
import friedmann as fr
import likelihoods as lk

H_FID = cfg.UNCALIBRATED_H0 / 100.0
OB_FID = cfg.UNCALIBRATED_OMBH2


def _np(obj):
    """Make numpy scalars/arrays JSON-friendly."""
    if isinstance(obj, dict):
        return {str(k): _np(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_np(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _np(obj.tolist())
    if isinstance(obj, (np.floating, float)):
        v = float(obj)
        return v if np.isfinite(v) else (None if np.isnan(v) else ("inf" if v > 0 else "-inf"))
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def _cached(name: str, force: bool, fn):
    """Run ``fn`` (returning a dict of arrays + a 'summary' dict) or load its cache."""
    path = cfg.DERIVED_DIR / f"{name}.npz"
    jpath = cfg.DERIVED_DIR / f"{name}.json"
    if path.exists() and jpath.exists() and not force:
        with np.load(path, allow_pickle=False) as z:
            arrays = {k: z[k] for k in z.files}
        return arrays, json.loads(jpath.read_text())
    t0 = time.monotonic()
    arrays, summary = fn()
    summary["seconds"] = round(time.monotonic() - t0, 1)
    np.savez_compressed(path, **arrays)
    jpath.write_text(json.dumps(_np(summary), indent=1))
    print(f"  {name}: {summary['seconds']:.0f} s")
    return arrays, json.loads(jpath.read_text())


class Data:
    """The likelihoods, built once."""

    def __init__(self):
        self.sn = lk.SupernovaLikelihood(fd.load_pantheon())
        self.bao = lk.BAOLikelihood(fd.load_bao())
        self.rec = bg.Recombination()
        self.cmb = {"pr4": lk.cmb_prior_from_chains(fd.load_planck_early_chains())}
        d3 = fd.load_d3_fit()
        if d3 is not None:
            self.cmb["d3"] = lk.cmb_prior_from_d3(d3)
        self.marg = {}
        for key, prior in self.cmb.items():
            self.marg[(key, False)] = lk.PhysicalDensityMarginal(prior, self.rec)
            self.marg[(key, True)] = lk.PhysicalDensityMarginal(prior, self.rec, self.bao)


def _summaries(p, fields: dict) -> dict:
    return {k: lk.weighted_summary(v, p) for k, v in fields.items()}


def _mixture(p, mean, var) -> dict:
    """Summary of a quantity whose conditional distribution at each grid point is
    Gaussian (mean, var): the mixture's mean and total variance."""
    p = p.ravel() / p.sum()
    mean, var = mean.ravel(), var.ravel()
    keep = p > 0
    m = float(np.sum(p[keep] * mean[keep]))
    v = float(np.sum(p[keep] * (var[keep] + (mean[keep] - m) ** 2)))
    return {"mean": m, "std": float(np.sqrt(v))}


# ---------------------------------------------------------------------------
# 1-2. Supernovae alone; BAO alone
# ---------------------------------------------------------------------------


def step_sn(D: Data):
    Om = np.linspace(*cfg.SN_OM[:2], cfg.SN_OM[2])
    OL = np.linspace(*cfg.SN_OL[:2], cfg.SN_OL[2])
    OMg, OLg = np.meshgrid(Om, OL, indexing="ij")
    m = bg.Model.from_plane(OMg.ravel(), OLg.ravel(), H_FID, OB_FID)
    chi2, off = D.sn.chi2(m)
    chi2 = chi2.reshape(OMg.shape)
    p = lk.grid_posterior(chi2, mask=OMg >= 0)  # D1's only prior: Omega_m >= 0
    Ok = m.omk.reshape(OMg.shape)
    q0 = 0.5 * OMg - OLg + m.Or.reshape(OMg.shape)
    summ = _summaries(p, {"Om": OMg, "OL": OLg, "Ok": Ok})
    summ["chi2_min"] = float(np.min(chi2))
    summ["n_sn"] = D.sn.n
    summ["p_accelerating"] = float(p[q0 < 0].sum())
    summ["levels"] = lk.credible_levels(p)
    return {"Om": Om, "OL": OL, "chi2": chi2, "p": p}, summ


def step_bao(D: Data):
    Om = np.linspace(0.20, 0.40, 401)
    Ok = np.linspace(-0.22, 0.25, 471)
    OMg, OKg = np.meshgrid(Om, Ok, indexing="ij")
    m = bg.Model.from_omk(OMg.ravel(), OKg.ravel(), H_FID, OB_FID)
    r = D.bao.m2lnl_uncalibrated(m)
    m2 = r["m2lnl"].reshape(OMg.shape)
    prior = (OKg >= cfg.PRIOR_OMK[0]) & (OKg <= cfg.PRIOR_OMK[1])
    p = lk.grid_posterior(m2, mask=prior)
    OL = m.OL.reshape(OMg.shape)
    summ = _summaries(p, {"Om": OMg, "OL": OL, "Ok": OKg, "hrd": r["hrd"].reshape(OMg.shape)})
    summ["chi2_min"] = float(np.min(r["chi2"]))
    summ["levels"] = lk.credible_levels(p)
    # DESI's own chain for exactly this run, and its best fit
    try:
        ch = fd.load_desi_chain("desi-bao-all")
        w = ch["weight"]
        summ["desi_chain"] = {k: lk.weighted_summary(ch[c], w) for k, c in
                              (("Om", "omegam"), ("Ok", "omk"), ("OL", "omegal"),
                               ("hrd", "hrdrag"))}
        summ["desi_chain"]["n_rows"] = int(w.size)
        summ["shift_sigma"] = {k: (summ[k]["mean"] - summ["desi_chain"][k]["mean"])
                               / summ["desi_chain"][k]["std"] for k in ("Om", "Ok", "OL")}
        summ["width_ratio"] = {k: summ[k]["std"] / summ["desi_chain"][k]["std"]
                               for k in ("Om", "Ok", "OL")}
    except FileNotFoundError:
        pass
    return {"Om": Om, "Ok": Ok, "OL": OL, "m2lnl": m2, "p": p}, summ


# ---------------------------------------------------------------------------
# 3. The CMB band
# ---------------------------------------------------------------------------


def _cmb_band(D: Data, key: str):
    """The CMB's likelihood on a grid that follows its ridge across the plane.

    Alone, the CMB measures the physical densities and the angle theta_*, but not
    h: for every Omega_m there is a curvature (and an h) that puts the
    last-scattering surface at the right distance. That is the famous
    "geometric degeneracy", a thin band. The grid is (Omega_m, offset from the
    ridge in Omega_k): a shear, so cells keep equal areas in the plane.
    """
    marg = D.marg[(key, False)]
    om_line = np.linspace(*cfg.CMB_BAND_OM[:2], cfg.CMB_BAND_OM[2])
    # 1. the ridge: coarse scan in Omega_k for every Omega_m, then refine
    ok_scan = np.linspace(-0.5, 0.25, 301)
    OM, OK = np.meshgrid(om_line, ok_scan, indexing="ij")
    r = marg(OM.ravel(), OK.ravel())
    m2 = r["m2lnl"].reshape(OM.shape)
    ridge = ok_scan[np.argmin(np.where(np.isfinite(m2), m2, np.inf), axis=1)]
    for _ in range(3):
        f = lambda dk: marg(om_line, ridge + dk)["m2lnl"]  # noqa: E731
        h = 2e-4
        f0, fp, fm = f(0.0), f(h), f(-h)
        with np.errstate(invalid="ignore"):
            curv = (fp - 2 * f0 + fm) / h**2
            step = -(fp - fm) / (2 * h) / np.where(curv > 0, curv, np.inf)
        ridge = ridge + np.clip(np.where(np.isfinite(step), step, 0.0), -0.01, 0.01)
    # 2. the band grid around the ridge
    dk = np.linspace(-cfg.CMB_BAND_HALFWIDTH, cfg.CMB_BAND_HALFWIDTH, cfg.CMB_BAND_N)
    OMg = np.repeat(om_line[:, None], dk.size, axis=1)
    OKg = ridge[:, None] + dk[None, :]
    r = marg(OMg.ravel(), OKg.ravel())
    shape = OMg.shape
    m2 = r["m2lnl"].reshape(shape)
    p = lk.grid_posterior(m2)
    OL = r["OL"].reshape(shape)
    H0 = 100 * r["h"].reshape(shape)
    summ = {"name": D.cmb[key].name, "chi2_min": float(np.nanmin(r["chi2"])),
            "levels": lk.credible_levels(p),
            "ridge": {"Om": om_line[::20].tolist(), "Ok": ridge[::20].tolist(),
                      "H0": (100 * marg(om_line[::20], ridge[::20])["h"]).tolist()},
            "edge_mass": float(p[:, [0, -1]].sum())}
    # where the ridge crosses the flat line
    i = np.argmin(np.abs(ridge))
    summ["flat_crossing"] = {"Om": float(om_line[i]), "Ok": float(ridge[i]),
                             "H0": float(100 * marg([om_line[i]], [ridge[i]])["h"][0])}
    return {"Om": om_line, "dk": dk, "ridge": ridge, "Ok": OKg, "OL": OL, "H0": H0,
            "m2lnl": m2, "p": p}, summ


# ---------------------------------------------------------------------------
# 4. Combinations
# ---------------------------------------------------------------------------

COMBOS = {
    "bao+cmb": {"label": "DESI BAO + CMB", "sn": False, "bao": "calibrated", "cmb": "pr4"},
    "sn+cmb": {"label": "Pantheon+ + CMB", "sn": True, "bao": None, "cmb": "pr4"},
    "sn+bao": {"label": "Pantheon+ + DESI BAO", "sn": True, "bao": "uncalibrated", "cmb": None},
    "sn+bao+cmb": {"label": "all three", "sn": True, "bao": "calibrated", "cmb": "pr4"},
    "bao+cmb_d3": {"label": "DESI BAO + CMB (D3's own fit)", "sn": False, "bao": "calibrated",
                   "cmb": "d3"},
    "sn+bao+cmb_d3": {"label": "all three, D3's CMB", "sn": True, "bao": "calibrated",
                      "cmb": "d3"},
}


def combo_eval(D: Data, spec: dict, Om, Ok) -> dict:
    Om = np.atleast_1d(np.asarray(Om, float)).ravel()
    Ok = np.atleast_1d(np.asarray(Ok, float)).ravel()
    out = {"m2lnl": np.zeros(Om.size)}
    if spec["cmb"]:
        r = D.marg[(spec["cmb"], spec["bao"] == "calibrated")](Om, Ok)
        out.update(r)
    else:
        m = bg.Model.from_omk(Om, Ok, H_FID, OB_FID)
        out["OL"] = m.OL
    if spec["bao"] == "uncalibrated":
        b = D.bao.m2lnl_uncalibrated(bg.Model.from_omk(Om, Ok, H_FID, OB_FID))
        out["m2lnl"] = out["m2lnl"] + b["m2lnl"]
        out["hrd"] = b["hrd"]
    if spec["sn"]:
        c2, off = D.sn.chi2(bg.Model.from_omk(Om, Ok, H_FID, OB_FID))
        out["m2lnl"] = out["m2lnl"] + c2
        out["sn_chi2"] = c2
    return out


def _locate(D: Data, spec: dict, start=(0.30, 0.0)) -> tuple[np.ndarray, np.ndarray]:
    """Best (Omega_m, Omega_k) and a local covariance from a quadratic fit."""
    f = lambda p: float(combo_eval(D, spec, [p[0]], [p[1]])["m2lnl"][0])  # noqa: E731
    best = optimize.minimize(f, start, method="Nelder-Mead",
                             options={"xatol": 1e-6, "fatol": 1e-7, "maxiter": 2000}).x
    steps = np.array([0.01, 0.01])
    for _ in range(3):
        u = np.linspace(-2, 2, 9)
        U1, U2 = np.meshgrid(u, u, indexing="ij")
        a = best[0] + steps[0] * U1.ravel()
        b = best[1] + steps[1] * U2.ravel()
        y = combo_eval(D, spec, a, b)["m2lnl"]
        ok = np.isfinite(y)
        X = np.column_stack([np.ones(ok.sum()), U1.ravel()[ok], U2.ravel()[ok],
                             U1.ravel()[ok] ** 2, U1.ravel()[ok] * U2.ravel()[ok],
                             U2.ravel()[ok] ** 2])
        c = np.linalg.lstsq(X, y[ok], rcond=None)[0]
        Hs = np.array([[2 * c[3], c[4]], [c[4], 2 * c[5]]])  # in step units, of -2 lnL
        cov_u = 2.0 * np.linalg.inv(Hs)
        cov = cov_u * np.outer(steps, steps)
        steps = 0.6 * np.sqrt(np.diag(cov))
    return best, cov


def _derived_age(Om, Ok, ob, obc, cov_x):
    """Age at the best physical densities, and its conditional variance from theirs."""
    base = bg.age_gyr(bg.Model.from_physical(Om, Ok, ob, obc))
    d1 = 1e-5
    g_ob = (bg.age_gyr(bg.Model.from_physical(Om, Ok, ob + d1, obc)) - base) / d1
    d2 = 1e-4
    g_obc = (bg.age_gyr(bg.Model.from_physical(Om, Ok, ob, obc + d2)) - base) / d2
    var = g_ob**2 * cov_x[0] + 2 * g_ob * g_obc * cov_x[1] + g_obc**2 * cov_x[2]
    return base, var


def step_combo(D: Data, name: str):
    spec = COMBOS[name]
    best, cov = _locate(D, spec)
    sig = np.sqrt(np.diag(cov))
    n = 161
    Om = best[0] + sig[0] * np.linspace(-6, 6, n)
    Ok = best[1] + sig[1] * np.linspace(-6, 6, n)
    OMg, OKg = np.meshgrid(Om, Ok, indexing="ij")
    r = combo_eval(D, spec, OMg.ravel(), OKg.ravel())
    shape = OMg.shape
    m2 = r["m2lnl"].reshape(shape)
    prior = (OMg >= cfg.PRIOR_OM[0]) & (OKg >= cfg.PRIOR_OMK[0]) & (OKg <= cfg.PRIOR_OMK[1])
    p = lk.grid_posterior(m2, mask=prior)
    OL = r["OL"].reshape(shape)
    fields = {"Om": OMg, "Ok": OKg, "OL": OL, "q0": 0.5 * OMg - OL}
    arrays = {"Om": Om, "Ok": Ok, "OL": OL, "m2lnl": m2, "p": p}
    summ = {"label": spec["label"], "best_fit": {"Om": float(best[0]), "Ok": float(best[1])},
            "local_cov": cov.tolist(), "levels": lk.credible_levels(p),
            "edge_mass": float(p[[0, -1], :].sum() + p[:, [0, -1]].sum())}
    if spec["cmb"]:
        h = r["h"].reshape(shape)
        cov_x = (r["sig_ombh2"] ** 2, r["corr"] * r["sig_ombh2"] * r["sig_omega_bc"],
                 r["sig_omega_bc"] ** 2)
        # conditional spread of h at fixed (Om, Ok): h = sqrt((omega_bc + nu) / Om)
        var_h = (0.5 * h.ravel() / (r["omega_bc"] + bg.OMEGA_NU_MAPPING)) ** 2 * cov_x[2]
        age, var_age = _derived_age(OMg.ravel(), OKg.ravel(), r["ombh2"], r["omega_bc"], cov_x)
        arrays.update({"H0": 100 * h, "age": age.reshape(shape),
                       "ombh2": r["ombh2"].reshape(shape),
                       "omega_bc": r["omega_bc"].reshape(shape)})
        summ["H0"] = _mixture(p, 100 * h, (100**2) * var_h)
        summ["age_gyr"] = _mixture(p, age, var_age)
        summ["ombh2"] = _mixture(p, r["ombh2"], cov_x[0])
        summ["omega_bc"] = _mixture(p, r["omega_bc"], cov_x[2])
        k = np.unravel_index(np.argmax(p), shape)
        summ["map_point"] = {"Om": float(OMg[k]), "Ok": float(OKg[k]), "OL": float(OL[k]),
                             "H0": float(100 * h[k]), "age_gyr": float(age.reshape(shape)[k]),
                             "ombh2": float(r["ombh2"].reshape(shape)[k]),
                             "omega_bc": float(r["omega_bc"].reshape(shape)[k]),
                             "chi2": float(r["chi2"].reshape(shape)[k])
                             + (float(r["sn_chi2"].reshape(shape)[k]) if spec["sn"] else 0.0)}
    summ.update(_summaries(p, fields))
    if spec["sn"] and not spec["cmb"]:
        t0H0 = fr.age_grid(OMg, OL, bg.OMEGA_R / H_FID**2)
        arrays["t0H0"] = t0H0
        summ["t0H0"] = lk.weighted_summary(t0H0, p)
    return arrays, summ


def step_flat(D: Data, names=("bao+cmb", "sn+bao+cmb", "bao+cmb_d3")):
    """The same combinations with space assumed flat (Omega_k = 0), for a
    like-for-like comparison with Planck's and DESI's flat-LCDM ages."""
    arrays, summ = {}, {}
    for name in names:
        spec = COMBOS[name]
        if spec["cmb"] and spec["cmb"] not in D.cmb:
            continue
        om = np.linspace(0.27, 0.34, 1401)
        r = combo_eval(D, spec, om, np.zeros(om.size))
        p = lk.grid_posterior(r["m2lnl"])
        cov_x = (r["sig_ombh2"] ** 2, r["corr"] * r["sig_ombh2"] * r["sig_omega_bc"],
                 r["sig_omega_bc"] ** 2)
        age, var_age = _derived_age(om, np.zeros(om.size), r["ombh2"], r["omega_bc"], cov_x)
        var_h = (0.5 * r["h"] / (r["omega_bc"] + bg.OMEGA_NU_MAPPING)) ** 2 * cov_x[2]
        arrays[f"{name}_Om"], arrays[f"{name}_p"] = om, p
        summ[name] = {"label": spec["label"] + " (flat)", "Om": lk.weighted_summary(om, p),
                      "H0": _mixture(p, 100 * r["h"], 1e4 * var_h),
                      "age_gyr": _mixture(p, age, var_age),
                      "edge_mass": float(p[[0, -1]].sum())}
    return arrays, summ


# ---------------------------------------------------------------------------
# 5. An independent MCMC for BAO + CMB
# ---------------------------------------------------------------------------


def step_mcmc(D: Data, nsteps: int = 4000, nwalkers: int = 48, seed: int = 20261019):
    """Sample (Omega_m, Omega_k, omega_b, omega_b + omega_c) directly, with no
    Laplace approximation, to check the grid method on the same likelihood."""
    import emcee

    marg = D.marg[("pr4", True)]

    def log_prob(theta):
        theta = np.atleast_2d(theta)
        Om, Ok, x = theta[:, 0], theta[:, 1], theta[:, 2:]
        out = np.full(theta.shape[0], -np.inf)
        ok = (Om > cfg.PRIOR_OM[0]) & (Om < cfg.PRIOR_OM[1]) & (Ok > cfg.PRIOR_OMK[0]) & \
             (Ok < cfg.PRIOR_OMK[1]) & (x[:, 0] > 0.005) & (x[:, 0] < 0.1) & (x[:, 1] > 0.05)
        if ok.any():
            r, h = marg.predictions(Om[ok], Ok[ok], x[ok])
            chi2 = np.einsum("nk,nk->n", r, r)
            good = (h > cfg.PRIOR_H0[0] / 100) & (h < cfg.PRIOR_H0[1] / 100) & np.isfinite(chi2)
            vals = np.where(good, -0.5 * chi2, -np.inf)
            out[ok] = vals
        return out

    rng = np.random.default_rng(seed)
    start = np.array([0.30, 0.002, 0.0222, 0.1414])
    scale = np.array([0.004, 0.0015, 0.0001, 0.0008])
    p0 = start + scale * rng.standard_normal((nwalkers, 4))
    sampler = emcee.EnsembleSampler(nwalkers, 4, log_prob, vectorize=True)
    sampler.random_state = np.random.RandomState(seed).get_state()  # emcee seeding trap
    sampler.run_mcmc(p0, nsteps, progress=False)
    burn = nsteps // 4
    tau = sampler.get_autocorr_time(discard=burn, tol=0)
    chain = sampler.get_chain(discard=burn, thin=5, flat=True)
    Om, Ok, ob, obc = chain.T
    m = bg.Model.from_physical(Om, Ok, ob, obc)
    age = bg.age_gyr(m)
    summ = {"nsteps": nsteps, "nwalkers": nwalkers, "burn": burn,
            "tau": tau.tolist(), "n_eff": float(chain.shape[0] * 5 / np.max(tau)),
            "acceptance": float(np.mean(sampler.acceptance_fraction))}
    w = np.ones(Om.size)
    summ.update({k: lk.weighted_summary(v, w) for k, v in
                 (("Om", Om), ("Ok", Ok), ("OL", m.OL), ("H0", 100 * m.h), ("age_gyr", age),
                  ("ombh2", ob), ("omega_bc", obc))})
    return {"Om": Om, "Ok": Ok, "OL": m.OL, "H0": 100 * m.h, "age": age}, summ


# ---------------------------------------------------------------------------
# 6. Our universe: parameters, epochs, past and future
# ---------------------------------------------------------------------------


def _solve_z(fn, lo, hi):
    return float(optimize.brentq(fn, lo, hi, xtol=1e-10))


def epochs(m: bg.Model) -> dict:
    """When each ingredient took over, and the calendar of our universe."""
    K = cfg.HUBBLE_TIME_GYR_X_H0 / 100.0
    from scipy import integrate

    def t_of_z(z):  # Gyr since the Big Bang
        a = 1.0 / (1.0 + z)
        f = lambda u: 2 * u**3 / np.sqrt(bg.p_phys(np.array([u * u]), m)[0, 0])  # noqa: E731
        return K * integrate.quad(f, 0.0, np.sqrt(a), epsabs=0, epsrel=1e-11, limit=200)[0]

    ocb, Ok = float(m.omega_cb[0]), float(m.omk[0])
    h2 = float(m.h[0]) ** 2
    ode = float(m.omega_de[0])
    # matter-radiation equality as CAMB defines it: all neutrinos counted as radiation
    om_rad = bg.OMEGA_GAMMA * (1.0 + cfg.N_EFF * bg.NU_PER_SPECIES)
    z_eq = ocb / om_rad - 1.0
    # matter (incl. the now-slow neutrino) = Lambda
    om = ocb + bg.OMEGA_NU
    z_mL = (ode / om) ** (1.0 / 3.0) - 1.0

    def q(z):  # sign of the deceleration, rho + 3p; curvature exerts no force
        a = 1.0 / (1.0 + z)
        return 0.5 * om / a**3 + bg.OMEGA_R / a**4 - ode  # the slow neutrino is matter here

    z_acc = _solve_z(q, 0.0, 5.0)
    rec = bg.Recombination()
    zstar = float(rec("zstar", m.ombh2, m.omch2)[0])
    zdrag = float(rec("zdrag", m.ombh2, m.omch2)[0])
    age = t_of_z(0.0)
    out = {"age_gyr": age, "z_eq": z_eq, "t_eq_yr": t_of_z(z_eq) * 1e9,
           "z_star": zstar, "t_star_yr": t_of_z(zstar) * 1e9,
           "z_drag": zdrag, "r_drag": float(rec.r_drag(m)[0]),
           "z_acc": z_acc, "t_acc_gyr": t_of_z(z_acc), "lookback_acc_gyr": age - t_of_z(z_acc),
           "z_matter_lambda": z_mL, "t_matter_lambda_gyr": t_of_z(z_mL),
           "lookback_matter_lambda_gyr": age - t_of_z(z_mL),
           "Om": float(m.Om[0]), "OL": float(m.OL[0]), "Ok": Ok, "H0": 100 * float(m.h[0]),
           "Or": float(m.Or[0]), "omega_cb": ocb}
    # curvature takes over (never, if Lambda > 0 dominates first) -- for the record
    out["curvature_radius_gly"] = (cfg.C_KM_S / (100 * m.h[0]) / np.sqrt(abs(Ok)) * 3.26156e-3
                                   if Ok != 0 else np.inf)
    # e-folding time of the far future, 1 / (H0 sqrt(OL)) = (1/100) / sqrt(omega_Lambda)
    out["efold_time_gyr"] = K / np.sqrt(ode)
    out["Omega_m_h2"] = om
    out["h2"] = h2
    return out


def density_history(m: bg.Model, n: int = 1200) -> dict:
    """Shares of the critical density over cosmic time, past and future."""
    a = np.geomspace(1e-9, 1e2, n)
    ocb, ob = float(m.omega_cb[0]), float(m.ombh2[0])
    parts = {"photons": bg.OMEGA_GAMMA / a**4,
             "neutrinos": (bg.OMEGA_R - bg.OMEGA_GAMMA) / a**4 + bg.nu_a4(np.minimum(a, 1.0))
             / np.minimum(a, 1.0) ** 4 * np.where(a > 1, a ** -3.0, 1.0),
             "dark matter": (ocb - ob) / a**3, "atoms": ob / a**3,
             "dark energy": np.full(a.size, float(m.omega_de[0])),
             "curvature": float(m.omega_k[0]) / a**2}
    total = sum(parts.values())
    K = cfg.HUBBLE_TIME_GYR_X_H0 / 100.0
    # time at each a: cumulative integral of a da / sqrt(P), P = a^4 * total
    x = np.log(a)
    f = a**2 / np.sqrt(a**4 * total)  # dt/dlna in 1/100 units = a / (a^2 H) * a
    t = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(x))])
    t += a[0] ** 2 / (2 * np.sqrt(a[0] ** 4 * total[0]))  # radiation era before a[0]
    return {"a": a, "t_gyr": t * K, **{k: v / total for k, v in parts.items()}}


def step_ours(D: Data, combos: dict):
    """Our universe = the mean of the all-three posterior (curvature free)."""
    s = combos["sn+bao+cmb"]
    Om, Ok = s["Om"]["mean"], s["Ok"]["mean"]
    ob, obc = s["ombh2"]["mean"], s["omega_bc"]["mean"]
    m = bg.Model.from_physical(Om, Ok, ob, obc)
    ep = epochs(m)
    hist = density_history(m)
    # past and future of 200 universes drawn from the posterior grid
    arrays, _ = _cached("combo_sn+bao+cmb", False, lambda: None)
    p = arrays["p"].ravel()
    rng = np.random.default_rng(7)
    idx = rng.choice(p.size, size=200, p=p / p.sum())
    OMg, OKg = np.meshgrid(arrays["Om"], arrays["Ok"], indexing="ij")
    t_grid = np.linspace(-14.5, 60.0, 600)
    curves = []
    for i in idx:
        H0 = float(arrays["H0"].ravel()[i])
        u = fr.Universe(Om=float(OMg.ravel()[i]), OL=float(arrays["OL"].ravel()[i]),
                        Or=fr.omega_r_today(H0), H0=H0)
        pieces = u.history(t_grid[0], t_grid[-1], n=500)
        tt, aa = pieces[0]
        curves.append(np.interp(t_grid, tt, aa, left=np.nan, right=np.nan))
    curves = np.array(curves)
    best = fr.Universe(Om=float(m.Om[0]), OL=float(m.OL[0]), Or=float(m.Or[0]),
                       H0=100 * float(m.h[0]), name="our universe")
    bt, ba = best.history(t_grid[0], t_grid[-1], n=800)[0]
    fate = best.fate()
    # the future, in numbers
    fut = {"a_in_10gyr": float(np.interp(10.0, bt, ba)), "a_in_50gyr": float(np.interp(50.0, bt, ba)),
           "doubling_time_gyr": float(np.interp(2.0, ba, bt))}
    arrays_out = {"t_grid": t_grid, "curves": curves, "best_t": bt, "best_a": ba,
                  **{f"hist_{k}": v for k, v in hist.items()}}
    summ = {"params": {"Om": float(m.Om[0]), "OL": float(m.OL[0]), "Ok": float(m.omk[0]),
                       "H0": 100 * float(m.h[0]), "ombh2": ob, "omega_bc": obc,
                       "Or": float(m.Or[0])},
            "epochs": ep, "fate": fate, "future": fut,
            "density_today": {k: float(np.interp(1.0, hist["a"], v)) for k, v in hist.items()
                              if k not in ("a", "t_gyr")}}
    return arrays_out, summ


# ---------------------------------------------------------------------------
# 7. The gallery of universes
# ---------------------------------------------------------------------------


def gallery_universes(H0: float, ours: dict) -> list[tuple[str, fr.Universe, str]]:
    """(key, universe, one-line story). All expand at today's H0 at a = 1."""
    Or = fr.omega_r_today(H0)
    nbb = lambda Om: float(fr.no_big_bang_boundary([Om], Or)[0])  # noqa: E731
    U = lambda **kw: fr.Universe(H0=H0, **kw)  # noqa: E731
    return [
        ("ours", U(Om=ours["Om"], OL=ours["OL"], Or=Or, name="Our universe"),
         "Big Bang, then a slowing expansion that began to speed up ~6 billion years ago"),
        ("eds", U(Om=1.0, OL=0.0, Or=0.0, name="Einstein–de Sitter"),
         "matter only, exactly flat: a ~ t^(2/3), slowing forever (the pre-1998 favourite)"),
        ("milne", U(Om=0.0, OL=0.0, Or=0.0, name="Empty (Milne)"),
         "nothing at all: a ~ t, coasting; age exactly 1/H0"),
        ("radiation", U(Om=0.0, OL=0.0, Or=1.0, name="Radiation only"),
         "light only: a ~ t^(1/2), the shortest-lived past"),
        ("desitter", U(Om=0.0, OL=1.0, Or=0.0, name="de Sitter"),
         "vacuum energy only: exponential growth, no beginning"),
        ("open", U(Om=0.3, OL=0.0, Or=Or, name="Open, matter only"),
         "too little matter to stop: expands forever, slowing"),
        ("closed", U(Om=3.0, OL=0.0, Or=Or, name="Closed, matter only"),
         "Friedmann's 1922 closed universe: a cycloid into a Big Crunch"),
        ("negative", U(Om=0.3, OL=-0.3, Or=Or, name="Negative Lambda"),
         "a negative cosmological constant always wins: recollapse"),
        ("closed_lambda", U(Om=3.0, OL=0.1, Or=Or, name="Closed, weak Lambda"),
         "Lambda too weak to save a heavy universe: a later Big Crunch"),
        ("loiter", U(Om=0.1, OL=nbb(0.1) - 0.005, Or=Or, name="Loitering (Lemaître)"),
         "hovers near Einstein's static size at z ~ 2 (Petrosian, Salpeter & Szekeres 1967)"),
        ("bounce", U(Om=0.1, OL=1.6, Or=Or, name="Bounce"),
         "no Big Bang: shrank from infinity to a minimum size and re-expanded"),
        ("phantom", U(Om=ours["Om"], OL=1.0 - ours["Om"], Or=0.0, w=-1.5, name="Phantom (Big Rip)"),
         "dark energy with w = -1.5 grows as space expands and tears it apart"),
    ]


def step_gallery(ours_summary: dict):
    H0 = ours_summary["params"]["H0"]
    ours = ours_summary["params"]
    arrays, items = {}, []
    for key, u, story in gallery_universes(H0, ours):
        f = u.fate()
        pieces = u.history(-40.0, 60.0, n=700)
        for j, (t, a) in enumerate(pieces):
            arrays[f"{key}_t{j}"] = t
            arrays[f"{key}_a{j}"] = a
        f.update({"key": key, "story": story, "n_pieces": len(pieces)})
        items.append(f)
    return arrays, {"H0": H0, "universes": items}


# ---------------------------------------------------------------------------
# 8. The map of the plane
# ---------------------------------------------------------------------------


def step_map(H0: float):
    Or = fr.omega_r_today(H0)
    Om = np.linspace(*cfg.MAP_OM[:2], cfg.MAP_OM[2])
    OL = np.linspace(*cfg.MAP_OL[:2], cfg.MAP_OL[2])
    OMg, OLg = np.meshgrid(Om, OL, indexing="ij")
    code = fr.classify_grid(OMg, OLg, Or)
    t0H0 = fr.age_grid(OMg, OLg, Or)
    a_min = fr.bounce_size_grid(OMg, OLg, Or)
    om_line = np.linspace(0.0, cfg.MAP_OM[1], 3001)
    arrays = {"Om": Om, "OL": OL, "code": code.astype(np.int8), "t0H0": t0H0, "a_min": a_min,
              "line_Om": om_line,
              "nbb": fr.no_big_bang_boundary(om_line, Or), "nbb_norad": fr.no_big_bang_boundary(om_line),
              "recollapse": fr.recollapse_boundary(om_line, Or)}
    # static sizes along the no-Big-Bang edge: the redshift a bouncing universe can reach
    a_star = np.linspace(fr.smallest_physical_static_a(Or), 0.9, 2000)
    om_c, ol_c = fr.double_root_curve(a_star, Or)
    arrays.update({"edge_a_star": a_star, "edge_Om": om_c, "edge_OL": ol_c})
    summ = {"H0": H0, "Or": Or,
            "max_bounce_z_radiation": fr.max_bounce_redshift(Or),
            "max_bounce_z_no_radiation": "unbounded (as Omega_m -> 0)",
            "a_star_min": fr.smallest_physical_static_a(Or),
            "nbb_at_Om0_with_radiation": float(fr.no_big_bang_boundary([0.0], Or)[0]),
            "max_bounce_z_atoms": fr.max_bounce_redshift(Or, 0.049),
            "max_bounce_z_measured_Om": fr.max_bounce_redshift(Or, 0.30),
            "fractions": {name: float(np.mean(code == c)) for c, name in
                          enumerate(("accelerating", "decelerating", "recollapse", "bounce"))}}
    # sparkline lattice: a(t) at a lattice of points in the plane (off the boundary lines)
    lat_Om = np.arange(0.25, 3.0, 0.5)
    lat_OL = np.arange(-0.75, 3.0, 0.5)
    for i, om in enumerate(lat_Om):
        for j, ol in enumerate(lat_OL):
            u = fr.Universe(Om=float(om), OL=float(ol), Or=Or, H0=H0)
            try:
                pieces = u.history(-30.0, 40.0, n=200)
            except (ValueError, RuntimeError, ZeroDivisionError):
                continue
            for k, (t, a) in enumerate(pieces):
                arrays[f"spark_{i}_{j}_t{k}"] = t
                arrays[f"spark_{i}_{j}_a{k}"] = a
    arrays["lat_Om"], arrays["lat_OL"] = lat_Om, lat_OL
    return arrays, summ


# ---------------------------------------------------------------------------
# 9. Ages
# ---------------------------------------------------------------------------


def step_ages(map_arrays: dict, H0_values: dict):
    t0H0 = map_arrays["t0H0"]
    b6 = fd.load_b6_age_floor()
    arrays = {f"age_{k}": t0H0 * fr.hubble_time_gyr(H0) for k, H0 in H0_values.items()}
    tmin = b6["t_min_gyr"]
    summ = {"oldest_stars": b6, "H0_values": H0_values}
    for k, H0 in H0_values.items():
        tH = fr.hubble_time_gyr(H0)
        summ[k] = {"H0": H0, "hubble_time_gyr": tH, "eds_age_gyr": tH * 2 / 3,
                   "milne_age_gyr": tH,
                   "eds_vs_stars_sigma": (tmin - tH * 2 / 3) / b6["t_min_err"]}
    # the 1990s: Freedman et al.'s H0 = 80 and Bolte & Hogan's 15.8 Gyr clusters
    h90 = cfg.AGE_CRISIS_1990S["H0_freedman94"][0]
    summ["crisis_1990s"] = {"H0": h90, "eds_age_gyr": fr.hubble_time_gyr(h90) * 2 / 3,
                            "gc_age": cfg.AGE_CRISIS_1990S["gc_age_bolte_hogan95"],
                            "chaboyer_lower95": cfg.AGE_CRISIS_1990S["gc_age_chaboyer96_lower95"]}
    # the largest H0 the oldest stars allow, along the flat line
    om_line = np.linspace(0.05, 1.0, 96)
    t_flat = fr.age_grid(om_line, 1.0 - om_line, 0.0)
    arrays["flat_Om"] = om_line
    arrays["flat_H0_max"] = t_flat * cfg.HUBBLE_TIME_GYR_X_H0 / tmin
    arrays["flat_H0_max_p5"] = t_flat * cfg.HUBBLE_TIME_GYR_X_H0 / b6["t_min_p5"]
    return arrays, summ


# ---------------------------------------------------------------------------
# 10. Checks
# ---------------------------------------------------------------------------


def step_checks(D: Data):
    import camb

    out = {}
    # closed forms
    U = fr.Universe
    out["closed_forms"] = {
        "eds": [U(1, 0)._t_from_zero(1.0), fr.age_einstein_de_sitter()],
        "milne": [U(0, 0)._t_from_zero(1.0), 1.0],
        "radiation": [U(0, 0, Or=1.0)._t_from_zero(1.0), 0.5],
        "flat_lcdm_0.3": [U(0.3, 0.7)._t_from_zero(1.0), fr.age_flat_lcdm(0.3)],
        "open_matter_0.3": [U(0.3, 0)._t_from_zero(1.0), fr.age_open_matter(0.3)],
        "closed_matter_3_lifetime": [U(3.0, 0).fate()["lifetime_gyr"] / U(3.0, 0).t_hubble_gyr,
                                     fr.lifetime_closed_matter(3.0)],
    }
    # boundaries against Carroll, Press & Turner (1992) eqs. 11-12
    om = np.array([0.1, 0.3, 0.5, 0.8, 1.5, 2.5])
    x = (1 - om) / om
    cpt = np.where(om < 0.5, 4 * om * np.cosh(np.arccosh(np.maximum(x, 1)) / 3) ** 3,
                   4 * om * np.cos(np.arccos(np.clip(x, -1, 1)) / 3) ** 3)
    om2 = np.array([1.2, 2.0, 3.0, 5.0])
    rc = 4 * om2 * np.cos(np.arccos((1 - om2) / om2) / 3 + 4 * np.pi / 3) ** 3
    out["cpt92"] = {"no_big_bang_max_abs_diff": float(np.max(np.abs(fr.no_big_bang_boundary(om) - cpt))),
                    "recollapse_max_abs_diff": float(np.max(np.abs(fr.recollapse_boundary(om2) - rc)))}
    # the acceleration equation, integrated directly, against the quadrature
    errs = {}
    for name, u in (("ours", U(0.31, 0.69, Or=fr.OR_GALLERY)), ("closed", U(3.0, 0.0)),
                    ("bounce", U(0.3, 1.8)), ("negative", U(0.3, -0.5))):
        e = []
        for span in ((0, -3), (0, 3)):
            sol = fr.integrate_ode(u, span)
            for tg, a in u.history(-3 * u.t_hubble_gyr, 3 * u.t_hubble_gyr, n=400):
                th = tg / u.t_hubble_gyr
                msk = (th >= sol.t.min()) & (th <= sol.t.max()) & (a > 2e-3)
                if msk.any():
                    e.append(float(np.max(np.abs(sol.sol(th[msk])[0] - a[msk]))))
        errs[name] = max(e)
    out["ode_vs_quadrature_max_da"] = errs
    # CAMB, curvature included
    rows = []
    zs = np.array([0.295, 0.51, 0.706, 0.934, 1.321, 1.484, 2.33])
    for H0, ob, oc, ok in ((67.36, 0.02237, 0.1200, 0.0), (67.36, 0.02237, 0.1200, -0.05),
                           (55.0, 0.0222, 0.118, -0.2), (75.0, 0.0225, 0.121, 0.03),
                           (70.0, 0.0221, 0.125, 0.2), (45.0, 0.0224, 0.120, -0.25)):
        res = camb.get_background(bg.camb_params(ob, oc, H0=H0, omk=ok))
        d = res.get_derived_params()
        m = bg.Model(h=H0 / 100, ombh2=ob, omch2=oc, omk=ok)
        dm = bg.transverse_distance(bg.comoving_distance(zs, m), m)[0]
        # CAMB's own physical_time is accurate to ~2e-5; integrate its H(z) tightly
        from scipy import integrate
        t_camb = integrate.quad(lambda xx: 1.0 / res.hubble_parameter(np.expm1(xx)), 0,
                                np.log(1e8), epsrel=1e-12, epsabs=0, limit=500)[0] * \
            cfg.HUBBLE_TIME_GYR_X_H0
        rows.append({"H0": H0, "ombh2": ob, "omch2": oc, "omk": ok,
                     "theta_star": float(D.rec.theta_star(m)[0] / (d["thetastar"] / 100) - 1),
                     "r_drag": float(D.rec.r_drag(m)[0] / d["rdrag"] - 1),
                     "D_M": float(np.max(np.abs(dm / (res.angular_diameter_distance(zs) * (1 + zs)) - 1))),
                     "H": float(np.max(np.abs(bg.hubble(zs, m)[0] / res.hubble_parameter(zs) - 1))),
                     "age_minus_camb_quad_gyr": float(bg.age_gyr(m)[0] - t_camb),
                     "age_minus_camb_physical_time_gyr": float(bg.age_gyr(m)[0] - res.physical_time(0))})
    out["camb"] = rows
    # Planck 2018's best-fit parameters (TT,TE,EE+lowE+lensing means): age
    mp = bg.Model(h=0.6736, ombh2=0.02237, omch2=0.1200, omk=0.0)
    out["planck18_means_age_gyr"] = float(bg.age_gyr(mp)[0])
    om_rad = bg.OMEGA_GAMMA * (1.0 + cfg.N_EFF * bg.NU_PER_SPECIES)
    out["planck18_means_z_eq"] = float((0.02237 + 0.1200) / om_rad - 1)
    # DESI's best fits: our chi^2 at their parameters
    desi = {}
    try:
        bf = fd.load_desi_bestfit("desi-bao-all")
        mb = bg.Model.from_omk(bf["omm"], bf["omk"], H_FID, OB_FID)
        desi["bao"] = {"ours": float(D.bao.chi2_calibrated(mb, [bf["hrdrag"] / H_FID])[0]),
                       "desi": bf["chi2__BAO"]}
        bf = fd.load_desi_bestfit("desi-bao-all_pantheonplus")
        mb = bg.Model.from_omk(bf["omm"], bf["omk"], H_FID, OB_FID)
        desi["bao+sn"] = {"ours_bao": float(D.bao.chi2_calibrated(mb, [bf["hrdrag"] / H_FID])[0]),
                          "desi_bao": bf["chi2__BAO"],
                          "ours_sn": float(D.sn.chi2(mb)[0][0]), "desi_sn": bf["chi2__SN"]}
        bf = fd.load_desi_bestfit(f"desi-bao-all_{cfg.CMB_FULL}")
        mb = bg.Model(h=bf["H0"] / 100, ombh2=bf["ombh2"], omch2=bf["omch2"], omk=bf["omk"])
        desi["bao+cmb"] = {"ours_bao": float(D.bao.chi2_calibrated(mb, D.rec.r_drag(mb))[0]),
                           "desi_bao": bf["chi2__BAO"],
                           "ours_rdrag": float(D.rec.r_drag(mb)[0]), "desi_rdrag": bf["rdrag"],
                           "ours_age": float(bg.age_gyr(mb)[0]), "desi_age": bf["age"],
                           "ours_Om": float(mb.Om[0]), "desi_Om": bf["omegam"],
                           "ours_OL": float(mb.OL[0]), "desi_OL": bf["omegal"]}
    except FileNotFoundError:
        pass
    out["desi_bestfits"] = desi
    # DESI's chains with the full CMB, for the record
    chains = {}
    for combo, label in cfg.DESI_OMK_CHAINS.items():
        try:
            ch = fd.load_desi_chain(combo)
        except FileNotFoundError:
            continue
        w = ch["weight"]
        entry = {k: lk.weighted_summary(ch[c], w) for k, c in
                 (("Om", "omegam"), ("Ok", "omk"), ("OL", "omegal"), ("age_gyr", "age"))
                 if c in ch}
        if "H0" in ch:
            entry["H0"] = lk.weighted_summary(ch["H0"], w)
        entry["n_rows"] = int(w.size)
        chains[label] = entry
    out["desi_chains"] = chains
    return {}, out


# ---------------------------------------------------------------------------


def run(force: bool = False) -> dict:
    t_start = time.monotonic()
    print("Loading data and likelihoods")
    D = Data()
    results = {"inputs": {"n_sn": D.sn.n, "n_bao": D.bao.n,
                          "cmb_priors": {k: v.as_dict() for k, v in D.cmb.items()},
                          "d5_cmb_priors": fd.load_d5_cmb_priors(),
                          "oldest_stars": fd.load_b6_age_floor()}}
    print("1-2. Supernovae alone; BAO alone")
    _, results["sn"] = _cached("sn", force, lambda: step_sn(D))
    _, results["bao"] = _cached("bao", force, lambda: step_bao(D))
    print("3. The CMB band")
    results["cmb"] = {}
    for key in D.cmb:
        _, results["cmb"][key] = _cached(f"cmb_band_{key}", force, lambda k=key: _cmb_band(D, k))
    print("4. Combinations")
    results["combos"] = {}
    for name, spec in COMBOS.items():
        if spec["cmb"] and spec["cmb"] not in D.cmb:
            continue
        _, results["combos"][name] = _cached(f"combo_{name}", force,
                                             lambda n=name: step_combo(D, n))
    _, results["flat"] = _cached("flat", force, lambda: step_flat(D))
    print("5. MCMC check of BAO + CMB")
    _, results["mcmc"] = _cached("mcmc_bao_cmb", force, lambda: step_mcmc(D))
    print("6. Our universe")
    _, results["ours"] = _cached("ours", force, lambda: step_ours(D, results["combos"]))
    H0 = results["ours"]["params"]["H0"]
    print("7. The gallery")
    _, results["gallery"] = _cached("gallery", force, lambda: step_gallery(results["ours"]))
    print("8. The map of the plane")
    map_arrays, results["map"] = _cached("map", force, lambda: step_map(H0))
    print("9. Ages")
    d1 = fd.load_d1_h0()
    H0_values = {"cmb": results["combos"]["bao+cmb"]["H0"]["mean"], "cepheids": d1["H0"]}
    results["inputs"]["d1_h0"] = d1
    _, results["ages"] = _cached("ages", force, lambda: step_ages(map_arrays, H0_values))
    print("10. Checks")
    _, results["checks"] = _cached("checks", force, lambda: step_checks(D))
    results["seconds_total"] = round(time.monotonic() - t_start, 1)
    cfg.RESULTS_PATH.write_text(json.dumps(_np(results), indent=1))
    print(f"results -> {cfg.RESULTS_PATH.relative_to(cfg.PROJECT_DIR)} "
          f"({results['seconds_total']:.0f} s)")
    return results


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--force", action="store_true", help="recompute every cached step")
    args = parser.parse_args(argv)
    run(force=args.force)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
