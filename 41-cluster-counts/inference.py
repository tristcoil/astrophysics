"""Fitting the cluster counts: how clumpy is the universe?

Two fits, both Poisson likelihoods in 20 count-rate x 14 redshift bins:

1. **Consistency scan.** The count rate-mass relation fixed at Ghirardini
   et al.'s published values; (Omega_m, sigma_8) scanned on a grid with flat
   priors. This asks: given *their* mass calibration, do *our* selection
   model and counts land on *their* cosmology?
2. **Marginalised fit (MCMC).** The five relation parameters free, each with
   a Gaussian prior of the published width (median, larger 68% half-width).
   This propagates the mass-calibration uncertainty, which dominates the
   error budget. The priors come from the same paper's joint fit (counts +
   lensing), so they already carry some count information: the result is a
   consistency check with honest error bars, not an independent measurement.

Robustness variants rerun the scan with one ingredient changed at a time
(selection threshold, extrapolation, absorption, optical depth limit,
contamination weights, fixed cosmological parameters, a hydrostatic mass bias).
"""

from __future__ import annotations

import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import config as C
import cosmology as CO
import model as MD
import selection as S

_WORKER = {}


def variant_kernel(spec: dict) -> dict:
    """Build (or load) the selection kernel a variant needs."""
    tag = spec.get("kernel_tag", "default")
    sel = S.SelectionModel(S.fit_nodes(use_clerc=spec.get("use_clerc", True)),
                           extrapolation=spec.get("extrapolation", "linear"),
                           s50_scale=spec.get("s50_scale", 1.0))
    return MD.build_kernel(sel, tag=tag, use_nh=spec.get("use_nh", True), use_zvlim=spec.get("use_zvlim", True))


def _init_worker(spec: dict, data: np.ndarray):
    from threadpoolctl import threadpool_limits

    threadpool_limits(1)
    _WORKER["model"] = MD.CountsModel(variant_kernel(spec), area_deg2=spec.get("area_deg2", C.AREA_COSMO_DEG2))
    _WORKER["spec"] = spec
    _WORKER["data"] = data
    z_edges = MD.bins()[1]
    zmask = np.ones(data.shape[1], bool)
    if spec.get("z_max") is not None:
        zmask &= z_edges[1:] <= spec["z_max"] + 1e-9
    if spec.get("z_min") is not None:
        zmask &= z_edges[:-1] >= spec["z_min"] - 1e-9
    _WORKER["zmask"] = zmask[None, :]


def _cosmo(spec, om, s8):
    kw = {k: spec[k] for k in ("h", "Ob", "ns") if k in spec}
    return CO.Cosmology(Om=om, sigma8=s8, **kw)


def _scaling(spec, vals=None):
    sr = MD.Scaling(mass_scale=spec.get("mass_scale", 1.0))
    if vals is not None:
        sr.A_X, sr.B_X, sr.F_X, sr.G_X, sr.sigma_X = vals
    return sr


def _lnl(om, s8, scaling_vals=None):
    m, spec, d, zm = _WORKER["model"], _WORKER["spec"], _WORKER["data"], _WORKER["zmask"]
    N = m.expected(_cosmo(spec, om, s8), _scaling(spec, scaling_vals))
    return MD.loglike(N[:, zm[0]], d[:, zm[0]]), float(N[:, zm[0]].sum())


def _grid_row(args):
    om, s8s = args
    return [_lnl(om, s8) for s8 in s8s]


def grid_scan(spec: dict, data: np.ndarray, om_grid=None, s8_grid=None, processes: int = 6) -> dict:
    om_grid = np.linspace(*C.OM_GRID) if om_grid is None else om_grid
    s8_grid = np.linspace(*C.S8GRID_SIGMA8) if s8_grid is None else s8_grid
    import multiprocessing as mp

    variant_kernel(spec)  # build the cache once, before the workers start
    # "fork": workers inherit the imported modules (forkserver, Python 3.14's
    # default, would re-import the calling script, impossible from stdin)
    with ProcessPoolExecutor(processes, mp_context=mp.get_context("fork"),
                             initializer=_init_worker, initargs=(spec, data)) as ex:
        rows = list(ex.map(_grid_row, [(om, s8_grid) for om in om_grid]))
    lnl = np.array([[r[0] for r in row] for row in rows])
    ntot = np.array([[r[1] for r in row] for row in rows])
    return {"om": om_grid, "s8": s8_grid, "lnl": lnl, "ntot": ntot, **summarize_grid(om_grid, s8_grid, lnl)}


def _hpd_bounds(x, w, frac=0.68):
    """Central credible interval from weighted samples."""
    order = np.argsort(x)
    cw = np.cumsum(w[order]) / w.sum()
    lo = x[order][np.searchsorted(cw, (1 - frac) / 2)]
    hi = x[order][np.searchsorted(cw, 1 - (1 - frac) / 2)]
    med = x[order][np.searchsorted(cw, 0.5)]
    return float(med), float(lo), float(hi)


def summarize_grid(om, s8, lnl) -> dict:
    """Marginal constraints from a grid posterior (flat priors)."""
    OM, S8 = np.meshgrid(om, s8, indexing="ij")
    p = np.exp(lnl - lnl.max())
    # refine the grid by bilinear interpolation for smoother marginals
    from scipy.interpolate import RegularGridInterpolator

    fine_om = np.linspace(om[0], om[-1], 4 * om.size)
    fine_s8 = np.linspace(s8[0], s8[-1], 4 * s8.size)
    FO, FS = np.meshgrid(fine_om, fine_s8, indexing="ij")
    interp = RegularGridInterpolator((om, s8), lnl - lnl.max())
    fp = np.exp(interp(np.stack([FO.ravel(), FS.ravel()], -1)))
    S8v = FS.ravel() * np.sqrt(FO.ravel() / 0.3)
    i = np.unravel_index(np.argmax(lnl), lnl.shape)
    edge = bool(p[0, :].max() > 0.01 or p[-1, :].max() > 0.01 or p[:, 0].max() > 0.01 or p[:, -1].max() > 0.01)
    return {"best": {"Om": float(OM[i]), "sigma8": float(S8[i]), "lnl": float(lnl[i])},
            "Om": _hpd_bounds(FO.ravel(), fp), "sigma8": _hpd_bounds(FS.ravel(), fp),
            "S8": _hpd_bounds(S8v, fp), "touches_edge": edge}


# ---------------------------------------------------------------------------
# MCMC with the scaling relation free
# ---------------------------------------------------------------------------

PARAMS = ["Om", "sigma8", "A_X", "B_X", "F_X", "G_X", "sigma_X"]


def _log_post(theta):
    om, s8, *sv = theta
    if not (0.1 < om < 0.6 and 0.5 < s8 < 1.3):
        return -np.inf
    if sv[0] <= 0 or sv[4] <= 0.05:
        return -np.inf
    lp = 0.0
    for name, v in zip(PARAMS[2:], sv):
        mu, sd = C.SCALING_PUBLISHED[name]
        lp -= 0.5 * ((v - mu) / sd) ** 2
    return lp + _lnl(om, s8, sv)[0]


def run_mcmc(spec: dict, data: np.ndarray, nwalkers=C.MCMC_WALKERS, nsteps=C.MCMC_STEPS,
             burn=C.MCMC_BURN, seed=C.MCMC_SEED, processes: int = 6) -> dict:
    import multiprocessing as mp

    import emcee

    variant_kernel(spec)
    rng = np.random.default_rng(seed)
    start = np.array([0.29, 0.88] + [C.SCALING_PUBLISHED[p][0] for p in PARAMS[2:]])
    scale = np.array([0.01, 0.01] + [0.3 * C.SCALING_PUBLISHED[p][1] for p in PARAMS[2:]])
    p0 = start + scale * rng.standard_normal((nwalkers, len(PARAMS)))
    t0 = time.time()
    ctx = mp.get_context("fork")
    with ctx.Pool(processes, initializer=_init_worker, initargs=(spec, data)) as pool:
        sampler = emcee.EnsembleSampler(nwalkers, len(PARAMS), _log_post, pool=pool)
        sampler.random_state = np.random.RandomState(seed).get_state()  # see emcee seeding note
        sampler.run_mcmc(p0, nsteps, progress=False)
    chain = sampler.get_chain(discard=burn, flat=True)
    lnp = sampler.get_log_prob(discard=burn, flat=True)
    try:
        tau = sampler.get_autocorr_time(discard=burn, quiet=True)
    except Exception:  # noqa: BLE001
        tau = np.full(len(PARAMS), np.nan)
    S8v = chain[:, 1] * np.sqrt(chain[:, 0] / 0.3)
    w = np.ones(len(chain))
    summ = {p: _hpd_bounds(chain[:, i], w) for i, p in enumerate(PARAMS)}
    summ["S8"] = _hpd_bounds(S8v, w)
    corr = np.corrcoef(np.column_stack([chain, S8v]).T)
    return {"chain": chain[:: max(1, len(chain) // 20000)], "lnp": lnp[:: max(1, len(chain) // 20000)],
            "summary": summ, "acceptance": float(np.mean(sampler.acceptance_fraction)),
            "autocorr": tau.tolist(), "seconds": time.time() - t0,
            "corr": corr.tolist(), "corr_names": PARAMS + ["S8"]}


# ---------------------------------------------------------------------------
# Variants
# ---------------------------------------------------------------------------

BASELINE = {"label": "baseline", "kernel_tag": "default"}

VARIANTS = [
    {"label": "selection threshold 10% lower", "kernel_tag": "s50_0p9", "s50_scale": 0.9},
    {"label": "selection threshold 10% higher", "kernel_tag": "s50_1p1", "s50_scale": 1.1},
    {"label": "selection held flat beyond z = 0.5", "kernel_tag": "flat", "extrapolation": "flat"},
    {"label": "only z_hat < 0.5 (no extrapolation needed)", "kernel_tag": "default", "z_max": 0.5},
    {"label": "only z_hat > 0.15 (drops the lowest bin)", "kernel_tag": "default", "z_min": 0.15},
    {"label": "selection without the 2000 s anchor", "kernel_tag": "noclerc", "use_clerc": False},
    {"label": "no Galactic-absorption term", "kernel_tag": "nonh", "use_nh": False},
    {"label": "no optical depth limit (ZVLIM)", "kernel_tag": "nozv", "use_zvlim": False},
    {"label": "no contamination weights", "kernel_tag": "default", "weighted": False},
    {"label": "area of our own mask (12,695 deg^2)", "kernel_tag": "default", "area_deg2": None},
    {"label": "H0 = 67.17 (prior centre - 1 sigma)", "kernel_tag": "default", "h": 0.6717},
    {"label": "H0 = 73.04 (distance ladder)", "kernel_tag": "default", "h": 0.7304},
    {"label": "n_s = 0.92", "kernel_tag": "default", "ns": 0.92},
    {"label": "n_s = 1.00", "kernel_tag": "default", "ns": 1.00},
    {"label": "Omega_b = 0.046", "kernel_tag": "default", "Ob": 0.046},
    {"label": "Omega_b = 0.052", "kernel_tag": "default", "Ob": 0.052},
]


def run_variants(df, processes: int = 6, om_grid=None, s8_grid=None, verbose=True) -> list:
    om_grid = np.linspace(0.18, 0.42, 41) if om_grid is None else om_grid
    s8_grid = np.linspace(0.72, 1.04, 41) if s8_grid is None else s8_grid
    out = []
    for spec in VARIANTS:
        spec = dict(spec)
        if "area_deg2" in spec and spec["area_deg2"] is None:
            spec["area_deg2"] = float(variant_kernel(spec)["area_deg2_mask"])
        data = MD.data_counts(df) if spec.get("weighted", True) else MD.data_counts(df, weighted=False)
        t0 = time.time()
        res = grid_scan(spec, data, om_grid, s8_grid, processes)
        res["label"] = spec["label"]
        if verbose:
            print(f"  {spec['label']:<45s} Om {res['Om'][0]:.3f}  sigma8 {res['sigma8'][0]:.3f}  "
                  f"S8 {res['S8'][0]:.3f}  ({time.time() - t0:.0f} s)")
        out.append(res)
    return out


def default_processes() -> int:
    return max(1, min(6, (os.cpu_count() or 2) - 2))
