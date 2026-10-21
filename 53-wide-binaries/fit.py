"""Fit the forward model to an observed sample: the wide binary test itself.

Parameters (Banik et al. 2024, table 4, same priors):
    a_break [kAU] (1, 15), beta (-15, -1), gamma (0, 4), f_CB (0, 0.99),
    log10 k_CB (-3, -0.25), f_LOS (0, 0.6), alpha (-0.5, 1.5)
(our alpha range is narrower than their (-2, 3.6) because our orbit
libraries span -0.5 ... 1.5; their posterior and ours sit well inside it).

``max_like(alpha=...)`` maximises the binomial likelihood with alpha fixed
(their gradient-ascent step); ``run_mcmc`` samples all seven with emcee,
seeded through ``sampler.random_state`` (emcee copies numpy's global RNG at
construction, so a later np.random.seed would not reach it).
"""

from __future__ import annotations

import json
import time

import numpy as np
from scipy.optimize import minimize

from config import DERIVED_DIR
import forward as FW

NAMES = ["a_break", "beta", "gamma", "f_CB", "log_k_CB", "f_LOS", "alpha"]
PRIOR = np.array([[1.0, 15.0], [-15.0, -1.0], [0.0, 4.0], [0.0, 0.99], [-3.0, -0.25], [0.0, 0.6],
                  [-0.5, 1.5]])
START = np.array([5.6, -2.8, 2.0, 0.63, np.log10(2.5e-3), 0.02, 0.5])


def observed_counts(rsky_kau, vt, binning: FW.Binning):
    H, _, _ = np.histogram2d(rsky_kau, vt, bins=[binning.r_edges, binning.v_edges])
    return H


def library_for(tag: str, m1, m2, *, binning=None, n: int = 1_000_000, quiet: bool = False,
                rebuild: bool = False, **build_kw) -> FW.ModelLibrary:
    """FW.build_library cached as data/derived/modellib_<tag>.npz."""
    binning = binning or FW.Binning()
    path = DERIVED_DIR / f"modellib_{tag}.npz"
    if path.exists() and not rebuild:
        z = np.load(path, allow_pickle=True)
        b = FW.Binning(r_edges=z["r_edges"], v_edges=z["v_edges"], a_lo_kau=float(z["a_lo"]),
                       a_hi_kau=float(z["a_hi"]), n_a=int(z["n_a"]), n_e=int(z["n_e"]))
        return FW.ModelLibrary(binning=b, alphas=z["alphas"], log_k=z["log_k"], H0=z["H0"], H1=z["H1"],
                               N0=z["N0"], N1=z["N1"], n_samples=int(z["n_samples"]),
                               meta=json.loads(str(z["meta"])))
    t0 = time.time()
    lib = FW.build_library(m1, m2, binning=binning, n=n, quiet=quiet, mass_tag=tag, **build_kw)
    np.savez(path, r_edges=binning.r_edges, v_edges=binning.v_edges, a_lo=binning.a_lo_kau,
             a_hi=binning.a_hi_kau, n_a=binning.n_a, n_e=binning.n_e, alphas=lib.alphas, log_k=lib.log_k,
             H0=lib.H0, H1=lib.H1, N0=lib.N0, N1=lib.N1, n_samples=lib.n_samples,
             meta=json.dumps(lib.meta))
    if not quiet:
        print(f"    model library '{tag}' built in {time.time() - t0:.0f} s", flush=True)
    return lib


class Fitter:
    def __init__(self, lib: FW.ModelLibrary, counts, los, prior=PRIOR, fixed: dict | None = None,
                 acceptance=None, gauss: dict | None = None):
        """``gauss``: optional Gaussian priors {name: (mean, sigma)} on top of
        the uniform bounds (e.g. Hwang et al.'s gamma = 1.32 +- 0.09)."""
        self.model = FW.Model(lib, los, acceptance)
        self.counts = counts
        self.prior = np.array(prior, float)
        self.fixed = dict(fixed or {})
        self.gauss = dict(gauss or {})
        self.free = [i for i, n in enumerate(NAMES) if n not in self.fixed]

    def full(self, x):
        th = np.empty(len(NAMES))
        th[self.free] = x
        for n, v in self.fixed.items():
            th[NAMES.index(n)] = v
        return th

    def lnprior(self, th):
        lo, hi = self.prior[:, 0], self.prior[:, 1]
        if not np.all((th >= lo) & (th <= hi)):
            return -np.inf
        return float(sum(-0.5 * ((th[NAMES.index(n)] - m) / s) ** 2 for n, (m, s) in self.gauss.items()))

    def lnpost(self, x):
        th = self.full(x)
        lp = self.lnprior(th)
        if not np.isfinite(lp):
            return -np.inf
        p = self.model.pixel_probabilities(th)
        return FW.binomial_lnlike(p, self.counts) + lp

    def max_like(self, starts=None, seed: int = 0, nm_maxiter: int = 6000, powell_maxiter: int = 20000):
        rng = np.random.default_rng(seed)
        if starts is None:
            base = START.copy()
            lo_all, hi_all = self.prior[:, 0], self.prior[:, 1]
            outside = (base <= lo_all) | (base >= hi_all)
            base[outside] = 0.5 * (lo_all + hi_all)[outside]     # e.g. log A_max for absolute companions
            base = base[self.free]
            starts = [base] + [base + rng.normal(0, 1, len(base)) * np.array(
                [1.0, 0.3, 0.4, 0.08, 0.2, 0.01, 0.2])[self.free] for _ in range(5)]
        best = None
        for x0 in starts:
            lo, hi = self.prior[self.free, 0], self.prior[self.free, 1]
            x0 = np.clip(x0, lo + 1e-6, hi - 1e-6)
            r = minimize(lambda x: -self.lnpost(x) if np.isfinite(self.lnpost(x)) else 1e30, x0,
                         method="Nelder-Mead", options=dict(maxiter=nm_maxiter, maxfev=nm_maxiter, xatol=1e-4,
                                                            fatol=1e-4, adaptive=True))
            r = minimize(lambda x: -self.lnpost(x) if np.isfinite(self.lnpost(x)) else 1e30, r.x,
                         method="Powell", options=dict(maxiter=powell_maxiter, xtol=1e-4, ftol=1e-6))
            if best is None or r.fun < best.fun:
                best = r
        return self.full(best.x), -best.fun

    def run_mcmc(self, x0, *, n_walkers: int = 32, n_steps: int = 3000, seed: int = 1, burn: int = 1000):
        import emcee
        ndim = len(self.free)
        rng = np.random.default_rng(seed)
        scale = np.array([0.3, 0.1, 0.1, 0.01, 0.05, 0.003, 0.03])[self.free]
        p0 = np.array(x0)[self.free] + scale * rng.standard_normal((n_walkers, ndim))
        lo, hi = self.prior[self.free, 0], self.prior[self.free, 1]
        p0 = np.clip(p0, lo + 1e-6, hi - 1e-6)
        sampler = emcee.EnsembleSampler(n_walkers, ndim, self.lnpost)
        sampler.random_state = np.random.RandomState(seed).get_state()
        sampler.run_mcmc(p0, n_steps, progress=False)
        chain = sampler.get_chain(discard=burn, flat=True)
        lnp = sampler.get_log_prob(discard=burn, flat=True)
        try:
            tau = sampler.get_autocorr_time(discard=burn, quiet=True)
        except Exception:  # noqa: BLE001
            tau = np.full(ndim, np.nan)
        return chain, lnp, tau, float(np.mean(sampler.acceptance_fraction))


def summarize(chain, names):
    out = {}
    for j, n in enumerate(names):
        q = np.percentile(chain[:, j], [15.865, 50, 84.135])
        # mode from a histogram (Banik et al. quote modes)
        h, e = np.histogram(chain[:, j], bins=60)
        mode = 0.5 * (e[np.argmax(h)] + e[np.argmax(h) + 1])
        out[n] = dict(median=float(q[1]), lo=float(q[0]), hi=float(q[2]), mode=float(mode))
    return out
