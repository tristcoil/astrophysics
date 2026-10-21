"""Every number in the README and the data behind every figure, from the cached data.

Steps, each cached under ``data/derived/`` (``--force`` recomputes):

1. ``bk18``         our own BK18 likelihood, sampled: r and the foregrounds
2. ``models``       the inflation models solved exactly (background + Mukhanov-Sasaki)
3. ``constraints``  n_s and r from every chain; our Planck + BK18 combination;
                    2D credible regions; where each model lands
4. ``horizons``     the horizon and flatness problems in numbers
5. ``bmodes``       the B-mode spectra: gravitational waves, lensing, dust
6. ``checks``       our numbers against the published ones

    python analysis.py                 # use caches where present
    python analysis.py --force         # recompute everything
    python analysis.py --only bk18     # one step (and what it needs)
"""

from __future__ import annotations

import argparse
import json
import math
import time

import numpy as np

import config as cfg


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
    """Run ``fn`` (returning a dict of arrays + a summary dict) or load its cache."""
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
    print(f"  {name}: {summary['seconds']:.0f} s", flush=True)
    return arrays, json.loads(jpath.read_text())


def _interval(x: np.ndarray) -> dict:
    q = np.percentile(x, [2.275, 15.865, 50, 84.135, 97.725])
    return {"median": q[2], "lo68": q[1], "hi68": q[3], "lo95": q[0], "hi95": q[4],
            "mean": float(np.mean(x)), "std": float(np.std(x))}


# ---------------------------------------------------------------------------
# 1. Our own BK18 likelihood
# ---------------------------------------------------------------------------

#: r grid for BK18's marginal posterior (Rao-Blackwell; linear interpolation between).
R_GRID = np.linspace(0.0, 0.20, 201)


def step_bk18_mcmc():
    """Sample r and the seven foreground parameters with our BK18 likelihood."""
    import bk18

    like = bk18.BK18()
    th = bk18.camb_bb(bk18.BK18_COSMOLOGY, r=1.0, nt=0.0, lens_accuracy=2)
    like.set_cmb(th["lensed_scalar_bb"], th["tensor_bb"])
    res = bk18.sample(like)
    chain = res["chain"]
    best = int(np.argmax(res["log_post"]))
    p_best = dict(zip(bk18.PARAMS, chain[best]))
    summary = {
        "n_samples": int(chain.shape[0]), "tau": res["tau"],
        "n_independent": chain.shape[0] * res["thin"] / float(np.max(res["tau"])),
        "acceptance": res["acceptance"], "nwalkers": res["nwalkers"],
        "nsteps": res["nsteps"], "burn": res["burn"], "thin": res["thin"],
        "best_fit": p_best, "chi2_best": like.chi2(p_best),
        "n_bandpowers": like.nbins * like.ncl,
        "window_weight_dropped": like.window_weight_dropped,
        "params": {k: _interval(chain[:, j]) for j, k in enumerate(bk18.PARAMS)},
    }
    arrays = {"chain": chain, "log_post": res["log_post"],
              "dl_lens": th["lensed_scalar_bb"], "dl_tens_r1": th["tensor_bb"]}
    return arrays, summary


def _validation_worker(args):
    import bk18

    cosmo, r, fg = args
    th = bk18.camb_bb(cosmo, r=max(r, 1e-9), nt=None, lens_accuracy=2)
    like = _VALIDATION_LIKE
    cmb = (like.bin(th["lensed_scalar_bb"][like.lmin:like.lmax + 1]),
           like.bin(th["tensor_bb"][like.lmin:like.lmax + 1]))
    return like.chi2({**fg, "r": 1.0}, cmb)  # the tensor spectrum already has the sample's r


_VALIDATION_LIKE = None


def step_bk18_validation(n: int = 120, seed: int = 11, processes: int = 6):
    """Our BK18 chi^2 against Cobaya's own, sample by sample, in ACT DR6's Planck + BK18
    run: each sample's cosmology (through CAMB), r and foreground parameters."""
    import multiprocessing as mp

    import bk18
    import chains as ch

    global _VALIDATION_LIKE
    _VALIDATION_LIKE = bk18.BK18()
    d = ch.load("act_planck-bk_r_camb")
    rng = np.random.default_rng(seed)
    idx = rng.choice(d["weight"].size, n, replace=False, p=d["weight"] / d["weight"].sum())
    jobs = []
    for i in idx:
        cosmo = {"ombh2": d["ombh2"][i], "omch2": d["omch2"][i], "tau": d["tau"][i],
                 "ns": d["ns"][i], "theta100": 100 * d["theta"][i], "logA": d["logA"][i]}
        fg = {k: float(d[k][i]) for k in bk18.PARAMS if k != "r"}
        jobs.append((cosmo, float(d["r"][i]), fg))
    with mp.get_context("fork").Pool(processes) as pool:
        ours = np.array(pool.map(_validation_worker, jobs, chunksize=2))
    theirs = d["chi2_bk18"][idx]
    diff = ours - theirs
    summary = {"n": n, "chain": "act_planck-bk_r_camb (Cobaya bicep_keck_2018, CAMB 1.5.0)",
               "mean_diff": float(diff.mean()), "rms_diff": float(np.sqrt(np.mean(diff ** 2))),
               "max_abs_diff": float(np.max(np.abs(diff))),
               "chi2_range": [float(theirs.min()), float(theirs.max())]}
    return {"ours": ours, "cobaya": theirs, "r": d["r"][idx]}, summary


def step_bk18(mc: dict):
    """BK18's marginal r posterior from our chain (Rao-Blackwell), against BK18's own."""
    import bk18

    like = bk18.BK18()
    like.set_cmb(mc["dl_lens"], mc["dl_tens_r1"])
    rb = bk18.rao_blackwell_r(like, mc["chain"], R_GRID, n_use=3000)
    post = rb["posterior"] / rb["posterior"].max()
    r_pub, l_pub = bk18.published_r_likelihood()
    ours_at_pub = np.interp(r_pub, R_GRID, post)
    summary = {
        "method": "Rao-Blackwell average of p(r | foregrounds) over 3000 posterior samples",
        "r95": bk18._limit(R_GRID, post), "r95_mc_error": float(np.std(rb["boot_r95"])),
        "r68": bk18._limit(R_GRID, post, 0.6827),
        "r95_published_curve": bk18._limit(r_pub, l_pub),
        "r_peak": float(R_GRID[np.argmax(post)]),
        "r_peak_published": float(r_pub[np.argmax(l_pub)]),
        "L0_over_Lpeak": float(post[0]), "L0_over_Lpeak_published": float(l_pub[0]),
        "max_abs_diff_vs_published": float(np.max(np.abs(ours_at_pub - l_pub))),
        "r95_samples": float(np.percentile(mc["chain"][:, 0], 95)),
    }
    arrays = {"r_grid": R_GRID, "r_like": post, "r_like_boot_std": rb["boot_std"] /
              rb["posterior"].max(), "r_pub": r_pub, "l_pub": l_pub}
    return arrays, summary


# ---------------------------------------------------------------------------
# 2. The inflation models, solved exactly
# ---------------------------------------------------------------------------

#: Parameter grids for the families (in M_Pl units; alpha is dimensionless).
FAMILY_GRID_N = 16
N_FINE = np.arange(50.0, 60.01, 1.0)
#: Reheating scenarios that set N_*: (label, w, T_reh in GeV or None = instant).
REHEATING = (("instant", 0.0, None), ("w=0, T_reh=1e9 GeV", 0.0, 1e9),
             ("w=0, T_reh=4 MeV", 0.0, 4e-3), ("w=1, T_reh=4 MeV", 1.0, 4e-3))
#: Starobinsky's n_s over a wide range of N_*, and N_* against T_reh for several
#: equations of state during reheating.
STARO_N_SCAN = np.arange(40.0, 84.01, 2.0)
T_REH_SCAN = np.logspace(math.log10(4e-3), 16.0, 120)
W_SCAN = (0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0)


def _reh_consts() -> dict:
    return dict(k_star_mpc=cfg.K_PIVOT, h0_km_s_mpc=cfg.H0_FID, t0_k=cfg.T_CMB_K,
                g_s0=cfg.G_S_TODAY, mpc_inv_gev=cfg.MPC_INV_GEV,
                kb_gev_per_k=cfg.K_B_GEV_PER_K, c_km_s=cfg.C_KM_S)


FIELDS = ("n_star", "ns", "r", "alpha_s", "nt", "ln_pr", "ns_sr1", "r_sr1", "ns_sr2", "r_sr2",
          "V_star_quarter_gev", "H_star_gev", "delta_phi", "phi_star", "phi_end", "eps1",
          "eps2", "r_at_-3.219", "V0")


def _one_model(spec: dict, param, n_stars) -> dict:
    """Predictions of one potential (one parameter value) on a grid of N_*,
    plus N_* and the predictions for each reheating scenario."""
    import inflation as inf

    m = inf.model_from_spec(spec, param)
    off = (round(math.log(cfg.K_PIVOT_PLANCK_R / cfg.K_PIVOT), 3),)  # r at 0.002 Mpc^-1
    preds = inf.predictions(m, n_stars, ln10_10_as=cfg.LN10_10_AS, m_pl_gev=cfg.M_PL_GEV,
                            extra_offsets=off)
    row = {k: np.array([p[k] for p in preds]) for k in FIELDS}
    bg = inf.background(m, n_target=max(85.0, float(np.max(n_stars)) + 15.0))
    reh = []
    for _, w, t in REHEATING:
        n = inf.n_star_reheating(bg, ln10_10_as=cfg.LN10_10_AS, m_pl_gev=cfg.M_PL_GEV, w=w,
                                 t_reh_gev=t, consts=_reh_consts())
        p = inf.spectra_at(bg, n) if np.isfinite(n) else {"ns": np.nan, "r": np.nan}
        reh.append((n, p["ns"], p["r"]))
    row["reheating"] = np.array(reh)
    return row


def _family_params(spec: dict) -> np.ndarray:
    lo, hi = spec["range"]
    return np.geomspace(lo, hi, FAMILY_GRID_N)


def _model_worker(args):
    name, param, n_stars = args
    return name, param, _one_model(cfg.MODELS[name], param, n_stars)


def step_models(processes: int = 4):
    import multiprocessing as mp

    import inflation as inf

    jobs = []
    for name, spec in cfg.MODELS.items():
        if "param" in spec:
            for x in _family_params(spec):
                jobs.append((name, float(x), np.array(cfg.N_STAR_GRID)))
        else:
            jobs.append((name, None, N_FINE))
    with mp.get_context("fork").Pool(processes) as pool:
        done = pool.map(_model_worker, jobs, chunksize=1)
    arrays, summary = {}, {"models": {}, "reheating_scenarios": [r[0] for r in REHEATING],
                           "n_star_grid_fine": N_FINE, "n_star_grid_family": cfg.N_STAR_GRID}
    by_name: dict[str, list] = {}
    for name, param, row in done:
        by_name.setdefault(name, []).append((param, row))
    for name, rows in by_name.items():
        spec = cfg.MODELS[name]
        rows.sort(key=lambda t: (t[0] if t[0] is not None else 0.0))
        params = np.array([t[0] if t[0] is not None else np.nan for t in rows])
        for k in FIELDS:
            arrays[f"{name}__{k}"] = np.array([t[1][k] for t in rows])
        arrays[f"{name}__param"] = params
        arrays[f"{name}__reheating"] = np.array([t[1]["reheating"] for t in rows])
        s = {"label": spec["label"], "kind": spec["kind"]}
        if "param" in spec:
            s["param"], s["range"] = spec["param"], spec["range"]
        else:  # single-trajectory models: a table at N_* = 50, 55, 60
            row = rows[0][1]
            s["table"] = {f"N{int(n)}": {k: float(row[k][i]) for k in
                                         ("ns", "r", "alpha_s", "ns_sr1", "r_sr1", "ns_sr2",
                                          "r_sr2", "V_star_quarter_gev", "delta_phi")}
                          for i, n in enumerate(row["n_star"]) if n in (50.0, 55.0, 60.0)}
            s["reheating"] = {lab: {"n_star": float(v[0]), "ns": float(v[1]), "r": float(v[2])}
                              for (lab, _, _), v in zip(REHEATING, row["reheating"])}
        summary["models"][name] = s

    # Starobinsky over a wide range of N_*, and N_* against the reheating temperature
    m = inf.Model("starobinsky")
    bg = inf.background(m, n_target=100.0)
    scan = [inf.spectra_at(bg, float(n)) for n in STARO_N_SCAN]
    arrays["staro_scan_n"] = STARO_N_SCAN
    arrays["staro_scan_ns"] = np.array([p["ns"] for p in scan])
    arrays["staro_scan_r"] = np.array([p["r"] for p in scan])
    grid = np.full((len(W_SCAN), T_REH_SCAN.size), np.nan)
    for i, w in enumerate(W_SCAN):
        for j, t in enumerate(T_REH_SCAN):
            grid[i, j] = inf.n_star_reheating(bg, ln10_10_as=cfg.LN10_10_AS,
                                              m_pl_gev=cfg.M_PL_GEV, w=w, t_reh_gev=float(t),
                                              consts=_reh_consts())
    arrays["staro_treh"] = T_REH_SCAN
    arrays["staro_nstar_w"] = grid
    arrays["staro_w"] = np.array(W_SCAN)
    summary["starobinsky_nstar_range"] = {
        f"w={w:.2f}": [float(np.nanmin(grid[i])), float(np.nanmax(grid[i]))]
        for i, w in enumerate(W_SCAN)}
    return arrays, summary


def step_modes():
    """Single fluctuations followed through horizon exit (for the figures), and the
    slow-roll accuracy table."""
    import inflation as inf

    arrays, summary = {}, {}
    for name, model, n_star in (("starobinsky", inf.Model("starobinsky"), 55.0),
                                ("phi2", inf.Model("power", 2.0), 55.0)):
        bg = inf.background(model)
        lnk = bg.ln_ah(bg.n_end - n_star)
        h = inf.mode_history(bg, lnk, 1500)
        n_exit = bg.n_end - n_star
        for k, v in h.items():
            arrays[f"{name}__{k}"] = v
        arrays[f"{name}__n_since_exit"] = h["N"] - n_exit
        summary[name] = {"n_star": n_star, "n_exit": n_exit, "n_end": bg.n_end,
                         "final_PR": float(h["sqrt_PR"][-1] ** 2),
                         "final_Pt": float(h["sqrt_Pt"][-1] ** 2)}
    return arrays, summary


# ---------------------------------------------------------------------------
# 3. Constraints: n_s and r from the chains, and where the models land
# ---------------------------------------------------------------------------

NS_RANGE = (0.930, 1.000)
R_RANGE = (0.0, 0.30)
GRID_N = 420

#: The planes drawn in the keeper figure and used for the verdicts:
#: key -> (label, how it is built).
PLANES = {
    "planck18": ("Planck 2018 (TT,TE,EE+lowE+lensing)", "chain:planck18_r"),
    "planck18_bk18_ours": ("Planck 2018 + BK18 (our BK18 likelihood)", "ours"),
    "planck18_bk15": ("Planck 2018 + BK15 (Planck's run)", "chain:planck18_r_bk15"),
    "planck_bk18_act": ("Planck + BK18 (ACT's run)", "chain:act_planck-bk_r_camb"),
    "planck_bk18_lb_act": ("Planck + BK18 + lensing + BAO (ACT's run)",
                           "chain:act_planck-bk-l-b_r_camb"),
    "pact_bk18": ("P-ACT + BK18", "chain:act_p-actlite-bk_r_camb"),
    "pactlb_bk18": ("P-ACT-LB + BK18 (ACT DR6)", "chain:act_p-actlite-bk-l-b_r_camb"),
    "pactlb_x_bk18": ("P-ACT-LB × our BK18 r (product; validation)",
                      "product:act_p-actlite-l-b_lcdm_camb"),
    "cmbspa_x_bk18": ("CMB-SPA × our BK18 r (product)", "product:cmbspa"),
    "cmbspa_desi_x_bk18": ("CMB-SPA + DESI × our BK18 r (product)", "product:cmbspa_desi"),
}
#: Planes in which every model is graded.
VERDICT_PLANES = ("planck18_bk18_ours", "pactlb_bk18", "cmbspa_desi_x_bk18")


def _r_like_interp(bk: dict):
    grid, like = bk["r_grid"], bk["r_like"]
    return lambda r: np.interp(r, grid, like, right=0.0)


def _product_density(ns: np.ndarray, w: np.ndarray, r_like, xr, yr, n) -> dict:
    """p(n_s) from a LCDM chain times our BK18 r likelihood: the joint posterior if n_s
    and r were uncorrelated (validated on ACT, where the true joint chain exists)."""
    import chains as ch

    xe = np.linspace(*xr, n + 1)
    xc = 0.5 * (xe[:-1] + xe[1:])
    yc = 0.5 * (np.linspace(*yr, n + 1)[:-1] + np.linspace(*yr, n + 1)[1:])
    h, _ = np.histogram(ns, xe, weights=w)
    from scipy.ndimage import gaussian_filter1d
    sx = max(0.6, ch.n_eff(w) ** -0.2 * 1.06 * ch.mean_std(ns, w)[1] / (xe[1] - xe[0]))
    px = gaussian_filter1d(h, sx, mode="constant")
    p = np.outer(px, r_like(yc))
    p /= p.sum()
    return {"x": xc, "y": yc, "p": p, "levels": ch.hpd_levels(p, (0.6827, 0.9545)),
            "kernel_bins": (sx, 0.0), "mass_in_grid": 1.0}


def _verdicts(planes: dict, models: dict) -> dict:
    """For every model, the credibility of its (n_s, r) in each verdict plane:
    at N_* = 50, 55, 60 and at its instant-reheating N_*; for a family, the
    best parameter value."""
    import chains as ch

    out = {}
    for name, spec in cfg.MODELS.items():
        ns, r = models[f"{name}__ns"], models[f"{name}__r"]   # (nparam, nN)
        nst = models[f"{name}__n_star"][0]
        reh = models[f"{name}__reheating"]                     # (nparam, nscen, 3)
        params = models[f"{name}__param"]
        res = {}
        for pk in VERDICT_PLANES:
            d = planes[pk]
            cred = np.array([[ch.credibility(d, ns[i, j], r[i, j]) for j in range(ns.shape[1])]
                             for i in range(ns.shape[0])])
            inst = np.array([ch.credibility(d, reh[i, 0, 1], reh[i, 0, 2])
                             if np.isfinite(reh[i, 0, 1]) else 1.0 for i in range(ns.shape[0])])
            row = {}
            for nv in (50.0, 55.0, 60.0):
                j = int(np.argmin(np.abs(nst - nv)))
                i = int(np.argmin(cred[:, j]))
                row[f"N{int(nv)}"] = {"cred": float(cred[i, j]),
                                      "sigma": ch.sigma_equivalent(cred[i, j]),
                                      "beyond_resolution": bool(cred[i, j] > ch.CRED_RESOLVED),
                                      "param": None if np.isnan(params[i]) else float(params[i])}
            i = int(np.argmin(cred.min(axis=1)))
            j = int(np.argmin(cred[i]))
            row["best_in_band"] = {"cred": float(cred[i, j]), "sigma": ch.sigma_equivalent(cred[i, j]),
                                   "beyond_resolution": bool(cred[i, j] > ch.CRED_RESOLVED),
                                   "n_star": float(nst[j]),
                                   "param": None if np.isnan(params[i]) else float(params[i])}
            k = int(np.argmin(inst))
            row["instant_reheating"] = {"cred": float(inst[k]), "sigma": ch.sigma_equivalent(inst[k]),
                                        "beyond_resolution": bool(inst[k] > ch.CRED_RESOLVED),
                                        "n_star": float(reh[k, 0, 0]),
                                        "param": None if np.isnan(params[k]) else float(params[k])}
            res[pk] = row
        out[name] = res
    return out


def step_constraints(bk: dict, models: dict):
    import chains as ch

    arrays, summary = {}, {"chains": {}, "planes": {}}
    loaded = {}
    for key, entry in ch.catalogue().items():
        d = ch.load(key)
        loaded[key] = d
        w = d["weight"]
        s = {"label": entry["label"], "source": entry["source"], "n_samples": int(w.size),
             "n_eff": ch.n_eff(w), "ns": ch.summary(d["ns"], w)}
        if "tau" in d:
            s["tau"] = ch.summary(d["tau"], w)
        if "r" in d:
            s["r95"] = ch.upper_limit(d["r"], w)
            s["r_mean"] = ch.mean_std(d["r"], w)[0]
        if "r002" in d:
            s["r002_95"] = ch.upper_limit(d["r002"], w)
        summary["chains"][key] = s

    r_like = _r_like_interp(bk)
    # our combination: Planck 2018's r-free chain, importance-weighted by our BK18 r likelihood
    pl = loaded["planck18_r"]
    w_ours = pl["weight"] * r_like(pl["r"])
    pub_like = lambda r: np.interp(r, bk["r_pub"], bk["l_pub"], right=0.0)
    w_pub = pl["weight"] * pub_like(pl["r"])
    summary["ours"] = {
        "label": PLANES["planck18_bk18_ours"][0],
        "ns": ch.summary(pl["ns"], w_ours), "r95": ch.upper_limit(pl["r"], w_ours),
        "tau": ch.summary(pl["tau"], w_ours),
        "r002_95": ch.upper_limit(pl["r002"], w_ours),
        "n_eff": ch.n_eff(w_ours), "n_eff_before": ch.n_eff(pl["weight"]),
        "with_bk18_published_curve": {"ns": ch.summary(pl["ns"], w_pub),
                                      "r95": ch.upper_limit(pl["r"], w_pub)},
    }
    # the product approximation, validated where the true joint chain exists (ACT)
    for name in ("act_p-actlite-l-b_lcdm_camb", "cmbspa", "cmbspa_desi", "act_planck_lcdm_camb"):
        d = loaded[name]
        summary["chains"][name]["product_with_bk18_r95"] = ch.upper_limit(
            R_GRID, np.interp(R_GRID, bk["r_grid"], bk["r_like"]))

    for pk, (label, how) in PLANES.items():
        if how == "ours":
            dens = ch.density_2d(pl["ns"], pl["r"], w_ours, NS_RANGE, R_RANGE, n=GRID_N)
        elif how.startswith("chain:"):
            d = loaded[how.split(":", 1)[1]]
            dens = ch.density_2d(d["ns"], d["r"], d["weight"], NS_RANGE, R_RANGE, n=GRID_N)
        else:
            d = loaded[how.split(":", 1)[1]]
            dens = _product_density(d["ns"], d["weight"], r_like, NS_RANGE, R_RANGE, GRID_N)
        arrays[f"plane__{pk}__p"] = dens["p"].astype(np.float32)
        arrays[f"plane__{pk}__levels"] = np.array(dens["levels"])
        arrays["plane_x"], arrays["plane_y"] = dens["x"], dens["y"]
        # marginals of the plane itself (checks the smoothing)
        pns = dens["p"].sum(axis=1)
        pr = dens["p"].sum(axis=0)
        cr = np.cumsum(pr) / pr.sum()
        summary["planes"][pk] = {"label": label, "source": how,
                                 "levels": dens["levels"], "kernel_bins": dens["kernel_bins"],
                                 "ns_mean": float(np.sum(pns * dens["x"]) / pns.sum()),
                                 "r95_from_plane": float(np.interp(0.95, cr, dens["y"]))}
    planes = {pk: {"x": arrays["plane_x"], "y": arrays["plane_y"],
                   "p": arrays[f"plane__{pk}__p"].astype(float)} for pk in PLANES}
    # 1D r marginals (reflection-smoothed at r = 0) for the BK18 figure
    from scipy.ndimage import gaussian_filter1d
    edges = np.linspace(0.0, 0.25, 501)
    centres = 0.5 * (edges[:-1] + edges[1:])
    for name, (rv, wv) in {"planck18": (pl["r"], pl["weight"]),
                           "planck18_bk18_ours": (pl["r"], w_ours),
                           "planck_bk18_act": (loaded["act_planck-bk_r_camb"]["r"],
                                               loaded["act_planck-bk_r_camb"]["weight"]),
                           "pactlb_bk18": (loaded["act_p-actlite-bk-l-b_r_camb"]["r"],
                                           loaded["act_p-actlite-bk-l-b_r_camb"]["weight"])}.items():
        h, _ = np.histogram(np.concatenate([rv, -rv]), np.concatenate([-edges[::-1], edges[1:]]),
                            weights=np.concatenate([wv, wv]))
        # display curves: 2.5x Silverman's width (chains are correlated, n_eff overstates)
        sig = max(1.0, 2.5 * 1.06 * ch.mean_std(rv, wv)[1] * ch.n_eff(wv) ** -0.2
                  / (edges[1] - edges[0]))
        sm = gaussian_filter1d(h, sig, mode="constant")[edges.size - 1:]
        arrays[f"r1d__{name}"] = sm / sm.max()
    arrays["r1d_grid"] = centres
    summary["verdicts"] = _verdicts(planes, models)

    # validation of the product approximation on ACT: the true P-ACT-LB + BK18 joint
    # chain against P-ACT-LB (LCDM) times our BK18 r likelihood
    a = summary["planes"]["pactlb_bk18"]
    b = summary["planes"]["pactlb_x_bk18"]
    pa, pb = planes["pactlb_bk18"]["p"], planes["pactlb_x_bk18"]["p"]
    summary["product_validation"] = {
        "ns_mean_joint": a["ns_mean"], "ns_mean_product": b["ns_mean"],
        "r95_joint": a["r95_from_plane"], "r95_product": b["r95_from_plane"],
        "overlap": float(np.sum(np.sqrt(pa / pa.sum() * pb / pb.sum()))),  # Bhattacharyya
    }

    # D3's own n_s (TT only, tau fixed, r = 0)
    try:
        d3 = json.loads(cfg.D3_FIT.read_text())["parameters"]["ns"]
        summary["d3"] = {"value": d3["value"], "error": d3["error"], "source": str(
            cfg.D3_FIT.relative_to(cfg.ROOT_DIR))}
    except (FileNotFoundError, KeyError):
        summary["d3"] = {**cfg.D3_NS_FALLBACK, "source": "D3 README (fallback)"}
    return arrays, summary


# ---------------------------------------------------------------------------
# 4. The horizon and flatness problems
# ---------------------------------------------------------------------------


def _o1_curvature() -> dict:
    try:
        r = json.loads(cfg.O1_RESULTS.read_text())["combos"]["bao+cmb"]["Ok"]
        return {"mean": r["mean"], "std": r["std"],
                "source": str(cfg.O1_RESULTS.relative_to(cfg.ROOT_DIR)) + " [combos][bao+cmb]"}
    except (FileNotFoundError, KeyError):
        return {**cfg.O1_OK_FALLBACK, "source": "O1 README (fallback)"}


def _compare_with_p2() -> dict | None:
    """Our free-particle g*(T) against project P2's (interactions, lattice QCD, a full
    hadron gas, real neutrino decoupling), where P2 has been run."""
    import horizons as hz

    if not cfg.P2_CURVES.exists():
        return None
    with np.load(cfg.P2_CURVES) as z:
        t, g_rho, g_s = z["T"], z["g_rho"], z["g_s"]
    t_gev = t if t.max() < 1e5 else t / 1e3   # P2 tabulates in GeV (checked by range below)
    rows = {}
    for T in (1e4, 100.0, 10.0, 1.0, 0.2, 0.1565, 0.14, 0.1, 0.01, 1e-3, 1e-4):
        if not (t_gev.min() <= T <= t_gev.max()):
            continue
        ours = hz.g_star(T)
        theirs = (float(np.interp(math.log(T), np.log(t_gev[::-1]) if t_gev[0] > t_gev[-1] else np.log(t_gev),
                                  g_rho[::-1] if t_gev[0] > t_gev[-1] else g_rho)),
                  float(np.interp(math.log(T), np.log(t_gev[::-1]) if t_gev[0] > t_gev[-1] else np.log(t_gev),
                                  g_s[::-1] if t_gev[0] > t_gev[-1] else g_s)))
        rows[f"{T:g} GeV"] = {"ours": ours, "p2": theirs}
    above = [v for k, v in rows.items() if float(k.split()[0]) >= 0.01]
    worst = max(abs(v["ours"][0] / v["p2"][0] - 1) for v in above) if above else None
    # what it does to the e-fold counts: N depends on ln a_end, i.e. on g_s(T_reh)^(-1/3)
    # and T_reh ~ (rho / g)^(1/4); at T_reh ~ 1e15 GeV the two g* differ by < 2%
    hi = rows.get("10000 GeV")
    dn = None if hi is None else abs(math.log(hi["ours"][1] / hi["p2"][1]) / 3
                                     - math.log(hi["ours"][0] / hi["p2"][0]) / 4)
    return {"source": str(cfg.P2_CURVES.relative_to(cfg.ROOT_DIR)), "rows": rows,
            "max_rel_diff_above_10MeV": worst, "efold_shift_estimate": dn}


def step_horizons(models_summary: dict):
    import horizons as hz
    import inflation as inf

    horizon = hz.cmb_horizon()
    ok = _o1_curvature()
    bound = abs(ok["mean"]) + 1.96 * ok["std"]          # 95% bound on |Omega_k|
    # the inflationary part: Starobinsky, solved exactly, instant reheating
    n_star = models_summary["models"]["starobinsky"]["reheating"]["instant"]["n_star"]
    bg = inf.background(inf.Model("starobinsky"), n_target=100.0)
    res = inf.spectra_at(bg, n_star)
    scale = inf.energy_scale(res, cfg.LN10_10_AS, cfg.M_PL_GEV)
    hist = hz.history(res, bg, scale["V0"], omega_k_bound=bound)
    summary = {"horizon": horizon, "curvature": {**ok, "abs_bound_95": bound},
               "history": hist["summary"],
               "omega_dev_planck_temperature": hz.omega_dev_planck_temperature(bound),
               "starobinsky_instant": {"n_star": n_star, "ns": res["ns"], "r": res["r"],
                                       "V_star_quarter_gev": scale["V_star_quarter_gev"],
                                       "H_star_gev": scale["H_star_gev"],
                                       "e_folds_in_solution": bg.n_end}}
    summary["p2_comparison"] = _compare_with_p2()
    sky = hz.cmb_sky()
    arrays = dict(hist["arrays"])
    arrays.update({"sky_mollweide": sky.pop("mollweide"), "sky_zoom": sky.pop("zoom"),
                   "sky_zoom_center": sky.pop("zoom_center"),
                   "sky_moll_outline": sky.pop("moll_outline"),
                   "sky_moll_extent": sky.pop("moll_extent")})
    summary["sky"] = sky
    return arrays, summary


# ---------------------------------------------------------------------------
# 5. B-mode spectra
# ---------------------------------------------------------------------------

R_SHOW = (0.036, 0.01, 0.003, 0.001)


def _bk18_display():
    """BK18's released display products: component decomposition and bandpowers."""
    comp = np.loadtxt(cfg.BK18_DIR / "BK18_components_20210607.txt", comments="#")
    bp = np.loadtxt(cfg.BK18_DIR / "BK18_bandpowers_20210607.txt", comments="#")
    # bandpower columns: bin, ell, then (value, error) per spectrum in the header order;
    # the auto spectra come first: K95, B95, 150, 220
    autos = {"K95": 0, "B95": 1, "150": 2, "220": 3}
    out = {"comp": comp, "bp_ell": bp[:, 1]}
    for name, i in autos.items():
        out[f"bp_{name}"] = bp[:, 2 + 2 * i]
        out[f"bp_{name}_err"] = bp[:, 3 + 2 * i]
    return out


def step_bmodes(mc: dict):
    import bk18
    import camb

    pars = camb.set_params(H0=cfg.H0_FID, ombh2=cfg.OMBH2_FID, omch2=cfg.OMCH2_FID,
                           mnu=cfg.MNU_FID, tau=cfg.TAU_FID, ns=cfg.NS_FID,
                           As=math.exp(cfg.LN10_10_AS) * 1e-10, lmax=2500,
                           lens_potential_accuracy=2, WantTensors=True, r=1.0)
    pars.NonLinear = camb.model.NonLinear_lens
    pars.max_l_tensor = 1500
    pars.max_eta_k_tensor = 3600.0
    res = camb.get_results(pars)
    cls = res.get_cmb_power_spectra(pars, CMB_unit="muK")
    ell = np.arange(cls["total"].shape[0])
    arrays = {"ell": ell, "lens_bb": cls["lensed_scalar"][:, 2], "ee": cls["lensed_scalar"][:, 1],
              "tt": cls["lensed_scalar"][:, 0]}
    for r in R_SHOW:  # tensors with the consistency relation n_t = -r/8
        res.Params.InitPower.set_params(As=math.exp(cfg.LN10_10_AS) * 1e-10, ns=cfg.NS_FID,
                                        r=r, nt=-r / 8, pivot_scalar=0.05, pivot_tensor=0.05)
        res.power_spectra_from_transfer()
        arrays[f"tens_bb_r{r}"] = res.get_tensor_cls(CMB_unit="muK")[:, 2]
    # our dust model at 150 GHz (BICEP2/Keck) and 95 GHz (BICEP3), from our chain
    like = bk18.BK18()
    chain = mc["chain"][::12]
    i150, i95 = like.maps.index("BK18_150_B"), like.maps.index("BK18_B95e_B")
    lgrid = np.arange(20, 351)
    d150, d95 = [], []
    for th in chain:
        p = dict(zip(bk18.PARAMS, th))
        f = like.dust_scaling(p["BBbetadust"])
        shape = p["BBdust"] * (lgrid / 80.0) ** p["BBalphadust"]
        d150.append(shape * f[i150] ** 2)
        d95.append(shape * f[i95] ** 2)
    for name, d in (("dust150", np.array(d150)), ("dust95", np.array(d95))):
        arrays[f"{name}_l"] = lgrid
        arrays[f"{name}_q"] = np.percentile(d, [2.275, 15.865, 50, 84.135, 97.725], axis=0)
    disp = _bk18_display()
    arrays.update({f"bk18_{k}": v for k, v in disp.items()})
    summary = {
        "r_show": R_SHOW,
        "tensor_peak_ell": int(ell[2 + np.argmax(arrays["tens_bb_r0.036"][2:300])]),
        "tensor_peak_dl_r0.036": float(np.max(arrays["tens_bb_r0.036"][2:300])),
        "lensing_dl_at_80": float(arrays["lens_bb"][80]),
        "dust150_at_80": {"median": float(np.interp(80, lgrid, arrays["dust150_q"][2])),
                          "lo68": float(np.interp(80, lgrid, arrays["dust150_q"][1])),
                          "hi68": float(np.interp(80, lgrid, arrays["dust150_q"][3]))},
        "dust95_at_80": float(np.interp(80, lgrid, arrays["dust95_q"][2])),
        "ratio_dust150_to_r0.036_at_80": float(np.interp(80, lgrid, arrays["dust150_q"][2])
                                              / arrays["tens_bb_r0.036"][80]),
        "ratio_lens_to_r0.036_at_80": float(arrays["lens_bb"][80] / arrays["tens_bb_r0.036"][80]),
    }
    return arrays, summary


# ---------------------------------------------------------------------------
# 5b. Dense model grids for the interactive explorer
# ---------------------------------------------------------------------------

EXPLORER_N = np.arange(40.0, 80.01, 2.5)
EXPLORER_NPARAM = 24
EXPLORER_FIELDS = ("ns", "r", "alpha_s", "ns_sr1", "r_sr1", "V_star_quarter_gev", "delta_phi")


def _explorer_worker(args):
    name, param = args
    return name, param, _one_model(cfg.MODELS[name], param, EXPLORER_N)


def step_explorer_grid(processes: int = 4):
    import multiprocessing as mp

    jobs = []
    for name, spec in cfg.MODELS.items():
        if "param" in spec:
            lo, hi = spec["range"]
            jobs += [(name, float(x)) for x in np.geomspace(lo, hi, EXPLORER_NPARAM)]
        else:
            jobs.append((name, None))
    with mp.get_context("fork").Pool(processes) as pool:
        done = pool.map(_explorer_worker, jobs, chunksize=1)
    arrays = {"n_star": EXPLORER_N}
    by: dict[str, list] = {}
    for name, param, row in done:
        by.setdefault(name, []).append((param, row))
    for name, rows in by.items():
        rows.sort(key=lambda t: t[0] if t[0] is not None else 0.0)
        arrays[f"{name}__param"] = np.array([t[0] if t[0] is not None else np.nan for t in rows])
        for k in EXPLORER_FIELDS:
            arrays[f"{name}__{k}"] = np.array([t[1][k] for t in rows])
        arrays[f"{name}__reheating"] = np.array([t[1]["reheating"] for t in rows])
    return arrays, {"n_star": EXPLORER_N, "n_param": EXPLORER_NPARAM,
                    "reheating_scenarios": [r[0] for r in REHEATING]}


# ---------------------------------------------------------------------------
# 6. Checks against the published numbers
# ---------------------------------------------------------------------------


def _hpd_interval(grid, like, mass=0.6827):
    """Peak and highest-density interval of a 1D curve on a grid."""
    p = like / np.trapezoid(like, grid)
    order = np.argsort(p)[::-1]
    dx = np.gradient(grid)
    c = np.cumsum(p[order] * dx[order])
    keep = order[:np.searchsorted(c, mass) + 1]
    return float(grid[np.argmax(p)]), float(grid[keep].min()), float(grid[keep].max())


def step_checks(results: dict, bk_arrays: dict, mc_arrays: dict):
    P = cfg.PUBLISHED
    rows = []

    def add(quantity, ours, published, ref, verdict, note=""):
        rows.append({"quantity": quantity, "ours": ours, "published": published, "ref": ref,
                     "verdict": verdict, "note": note})

    def gauss(q, s, pkey, quantity, note=""):
        pub = P[pkey]
        d = (s["mean"] - pub["value"]) / pub["error"]
        add(quantity, f"{s['mean']:.4f} ± {s['std']:.4f}", f"{pub['value']:.4f} ± {pub['error']:.4f}",
            pub["ref"], f"{d:+.2f}σ", note)

    ch = results["constraints"]["chains"]
    gauss("n_s", ch["planck18_lcdm"]["ns"], "planck18_ns", "Planck 2018 n_s (LCDM)")
    gauss("n_s", ch["planck18_tt_lcdm"]["ns"], "planck18_tt_ns", "Planck 2018 TT n_s")
    add("Planck 2018 r_0.002, 95%", f"< {ch['planck18_r']['r002_95']:.4f}",
        f"< {P['planck18_r002']['value']}", P["planck18_r002"]["ref"], "match",
        "our chain is TT,TE,EE+lowE+lensing")
    add("Planck 2018 + BK15 r_0.002, 95%", f"< {ch['planck18_r_bk15']['r002_95']:.4f}",
        f"< {P['planck18_bk15_r002']['value']}", P["planck18_bk15_r002"]["ref"], "match")
    b = results["bk18"]
    add("BK18 alone: r_0.05, 95% (our likelihood)", f"< {b['r95']:.4f} ± {b['r95_mc_error']:.4f}",
        f"< {P['bk18_r95']['value']}", P["bk18_r95"]["ref"], "match")
    pk, lo, hi = _hpd_interval(bk_arrays["r_grid"], bk_arrays["r_like"])
    add("BK18 alone: r peak and 68% interval", f"{pk:.3f} +{hi - pk:.3f} −{pk - lo:.3f}",
        "0.014 +0.010 −0.011", P["bk18_r_peak"]["ref"], "match")
    add("BK18 alone: L(r=0)/L(peak)", f"{b['L0_over_Lpeak']:.3f}", f"{P['bk18_L0']['value']}",
        P["bk18_L0"]["ref"], "match")
    mcs = results["bk18_mcmc"]["params"]
    for key, pkey, name in (("BBdust", "bk18_Ad", "dust amplitude A_d [uK^2]"),
                            ("BBbetadust", "bk18_betad", "dust spectral index beta_d")):
        s = mcs[key]
        pub = P[pkey]
        add(f"BK18 {name}", f"{s['median']:.2f} +{s['hi68'] - s['median']:.2f} "
            f"−{s['median'] - s['lo68']:.2f}",
            f"{pub['value']} +{pub['plus']} −{pub['minus']}", pub["ref"], "match")
    a_sync95 = float(np.percentile(mc_arrays["chain"][:, 2], 95))
    add("BK18 synchrotron A_sync, 95%", f"< {a_sync95:.2f}", f"< {P['bk18_Async95']['value']}",
        P["bk18_Async95"]["ref"], "match", "one-sided 95% quantile of our chain")
    add("BK18 dust alpha_d (best fit)", f"{results['bk18_mcmc']['best_fit']['BBalphadust']:.2f}",
        f"{P['bk18_alphad_ml']['value']}", P["bk18_alphad_ml"]["ref"], "match",
        "maximum-posterior sample of our chain")
    v = results["bk18_validation"]
    add("BK18 chi^2: ours vs Cobaya's, per sample", f"rms {v['rms_diff']:.3f}, max {v['max_abs_diff']:.3f}",
        "identical", "Cobaya's chi2__bicep_keck_2018 in ACT DR6's Planck + BK18 chain",
        "match", f"{v['n']} posterior samples, each at its own cosmology")
    o = results["constraints"]["ours"]
    add("Planck 2018 + BK18: r_0.05, 95% (ours)", f"< {o['r95']:.4f}",
        f"< {ch['act_planck-bk_r_camb']['r95']:.4f} (ACT's joint run)",
        "ACT DR6.02 chain planck-bk_r_camb (plik-lite, Sroll2)", "match")
    add("Planck + BK18 + lensing + BAO: r_0.05, 95%", f"< {ch['act_planck-bk-l-b_r_camb']['r95']:.4f}",
        f"< {P['bk18_planck_bao_r95']['value']}", P["bk18_planck_bao_r95"]["ref"], "match",
        "ACT's chain; BK18's own uses 2018-era BAO")
    gauss("n_s", ch["act_actlite_lcdm_camb"]["ns"], "act_ns", "ACT DR6 n_s",
          "ACT-lite chain; paper: full multi-frequency likelihood")
    gauss("n_s", ch["act_p-actlite_lcdm_camb"]["ns"], "pact_ns", "P-ACT n_s", "ACT-lite chain")
    gauss("n_s", ch["act_p-actlite-l-b_lcdm_camb"]["ns"], "pactlb_ns", "P-ACT-LB n_s",
          "ACT-lite chain")
    gauss("n_s", ch["act_planck_lcdm_camb"]["ns"], "act_planck_ns", "Planck n_s (ACT's run)")
    add("P-ACT-LB + BK18: r_0.05, 95%", f"< {ch['act_p-actlite-bk-l-b_r_camb']['r95']:.4f}",
        f"< {P['pactlb_bk18_r95']['value']}", P["pactlb_bk18_r95"]["ref"], "match")
    gauss("n_s", ch["spt3g"]["ns"], "spt_ns", "SPT-3G D1 n_s",
          "released chain uses the lite likelihood")
    gauss("n_s", ch["spt3g_desi"]["ns"], "spt_desi_ns", "SPT-3G D1 + DESI n_s")
    gauss("n_s", ch["cmbspa"]["ns"], "spt_cmbspa_ns", "CMB-SPA n_s", "2026 bug-fix chains")
    gauss("n_s", ch["cmbspa_desi"]["ns"], "spt_cmbspa_desi_ns", "CMB-SPA + DESI n_s",
          "2026 bug-fix chains")
    m = results["models"]["models"]
    for key, pkey, name in (("phi2", "phi2_n60", "m²φ²"), ("phi4", "phi4_n60", "λφ⁴"),
                            ("starobinsky", "starobinsky_n60", "Starobinsky")):
        t = m[key]["table"]["N60"]
        add(f"{name} at N* = 60: n_s, r (exact)", f"{t['ns']:.4f}, {t['r']:.4f}",
            f"{P[pkey]['ns']}, {P[pkey]['r']} (slow roll)", "atlas spec (slow-roll estimates)",
            "match" if abs(t["ns"] - P[pkey]["ns"]) < 0.0015 else "differs",
            f"first-order slow roll: {t['ns_sr1']:.4f}, {t['r_sr1']:.4f}")
    a_s = math.exp(cfg.LN10_10_AS) * 1e-10
    coef = (1.5 * math.pi ** 2 * a_s * 0.01) ** 0.25 * cfg.M_PL_GEV
    v_bk15 = (1.5 * math.pi ** 2 * a_s * P["planck18_bk15_r002"]["value"]) ** 0.25 * cfg.M_PL_GEV
    add("energy scale V^1/4 at r = 0.01", f"{coef:.3e} GeV", "≈ 1e16 GeV (atlas spec)",
        "(3 pi^2 A_s r / 2)^(1/4) M_Pl with Planck 2018 A_s", "match")
    add("V^1/4 bound from r_0.002 < 0.056", f"< {v_bk15:.2e} GeV", "< 1.6e16 GeV",
        P["planck18_v_quarter_gev"]["ref"], "match")
    v_ours = (1.5 * math.pi ** 2 * a_s * o["r95"]) ** 0.25 * cfg.M_PL_GEV
    summary = {"rows": rows, "v_quarter_bound_ours_gev": v_ours, "energy_coefficient_gev": coef}
    return {"dummy": np.zeros(1)}, summary


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

STEPS = ("bk18", "models", "explorer", "modes", "constraints", "horizons", "bmodes", "checks")


def run(force: bool = False, only: str | None = None) -> dict:
    results = {}
    t0 = time.monotonic()
    if only in (None, "bk18"):
        mc, results["bk18_mcmc"] = _cached("bk18_mcmc", force, step_bk18_mcmc)
        _, results["bk18"] = _cached("bk18", force, lambda: step_bk18(mc))
        _, results["bk18_validation"] = _cached("bk18_validation", force, step_bk18_validation)
    if only in (None, "models"):
        _, results["models"] = _cached("models", force, step_models)
    if only in (None, "explorer"):
        _, results["explorer_grid"] = _cached("explorer_grid", force, step_explorer_grid)
    if only in (None, "modes"):
        _, results["modes"] = _cached("modes", force, step_modes)
    if only in (None, "constraints"):
        mc, _ = _cached("bk18_mcmc", False, step_bk18_mcmc)
        bk, _ = _cached("bk18", False, lambda: step_bk18(mc))
        mod, _ = _cached("models", False, step_models)
        _, results["constraints"] = _cached("constraints", force,
                                            lambda: step_constraints(bk, mod))
    if only in (None, "horizons"):
        _, ms = _cached("models", False, step_models)
        _, results["horizons"] = _cached("horizons", force, lambda: step_horizons(ms))
    if only in (None, "bmodes"):
        mc, _ = _cached("bk18_mcmc", False, step_bk18_mcmc)
        _, results["bmodes"] = _cached("bmodes", force, lambda: step_bmodes(mc))
    if only is None or only == "checks":
        need = {}
        for name, fn in (("bk18_mcmc", step_bk18_mcmc), ("bk18_validation", step_bk18_validation),
                         ("models", step_models)):
            need[name] = results.get(name) or _cached(name, False, fn)[1]
        mc_arr, _ = _cached("bk18_mcmc", False, step_bk18_mcmc)
        bk_arr, need["bk18"] = _cached("bk18", False, lambda: step_bk18(mc_arr))
        need["constraints"] = results.get("constraints") or json.loads(
            (cfg.DERIVED_DIR / "constraints.json").read_text())
        _, results["checks"] = _cached("checks", force or only == "checks",
                                       lambda: step_checks(need, bk_arr, mc_arr))
    results["seconds_total"] = round(time.monotonic() - t0, 1)
    if only is None:
        cfg.RESULTS_PATH.write_text(json.dumps(_np(results), indent=1))
    return results


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--force", action="store_true", help="recompute every cached step")
    ap.add_argument("--only", choices=STEPS, help="run one step")
    args = ap.parse_args(argv)
    run(force=args.force, only=args.only)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
