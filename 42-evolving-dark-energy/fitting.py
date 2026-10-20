"""Combining likelihoods, finding best fits, and sampling posteriors.

A ``Problem`` is one model (LCDM, wCDM or w0waCDM) confronted with one set of
data (any of: DESI BAO, a supernova sample, a CMB prior). Its parameters follow
DESI's choices exactly:

* without the CMB, BAO cannot tell H0 from r_d and measure only their product,
  so the parameters are Omega_m and h*r_d (plus w0, wa);
* with the CMB prior, the sound horizon is computed from omega_b and omega_c,
  and the parameters are H0, omega_b, Omega_m (plus w0, wa).

The supernova offset is always removed analytically. Priors are DESI's flat
ranges plus w0 + wa < 0.

Best fits ("MAP", maximum a posteriori) are found by Nelder-Mead from several
starting points; the significance of evolving dark energy is the drop in
chi^2 from LCDM to w0waCDM, Delta chi^2_MAP, converted to Gaussian sigmas for
two extra parameters exactly as DESI does. Posteriors are sampled with emcee.
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
from scipy import optimize, stats

import background as bg
from config import PRIORS, REQUIRE_W0_PLUS_WA_NEGATIVE, UNCALIBRATED_H0, UNCALIBRATED_OMBH2
from early_universe import EarlyUniverse
from likelihoods import BAOLikelihood, CMBPrior, SNLikelihood

MODEL_EXTRA = {"lcdm": [], "wcdm": ["w0"], "w0wa": ["w0", "wa"]}


class Problem:
    def __init__(self, model: str, *, bao: BAOLikelihood | None = None,
                 sn: SNLikelihood | None = None, cmb: CMBPrior | None = None,
                 early: EarlyUniverse | None = None, label: str = ""):
        if bao is None and sn is None and cmb is None:
            raise ValueError("a problem needs some data")
        self.model, self.bao, self.sn, self.cmb = model, bao, sn, cmb
        self.calibrated = cmb is not None
        self.label = label
        if self.calibrated:
            base = ["H0", "ombh2", "Om"]
            self.early = early if early is not None else EarlyUniverse()
        else:
            base = ["Om", "hrd"] if bao is not None else ["Om"]
            self.early = None
        self.names = base + MODEL_EXTRA[model]
        targets = []
        if bao is not None:
            targets.append(bao.unique_z)
        if sn is not None:
            targets.append(sn.z_hd)
        z_t = np.unique(np.concatenate(targets)) if targets else np.array([1.0])
        self.grid = bg.DistanceGrid(z_t, z_max=1400.0 if self.calibrated else float(z_t.max()))
        self._pos = {}
        if bao is not None:
            self._pos["bao"] = np.searchsorted(z_t, bao.unique_z)
        if sn is not None:
            self._pos["sn"] = np.searchsorted(z_t, sn.z_hd)

    # --- parameters -----------------------------------------------------------

    @property
    def ndim(self) -> int:
        return len(self.names)

    def free_names(self, profile_hrd: bool) -> list[str]:
        return [n for n in self.names if not (profile_hrd and n == "hrd")]

    def _get(self, theta, name, default=None):
        if name in self.names:
            return theta[:, self.names.index(name)]
        return np.full(theta.shape[0], default, float)

    def cosmology(self, theta: np.ndarray) -> bg.Cosmo:
        theta = np.atleast_2d(theta)
        w0 = self._get(theta, "w0", -1.0)
        wa = self._get(theta, "wa", 0.0)
        Om = self._get(theta, "Om")
        if self.calibrated:
            return bg.Cosmo.from_Om(Om, self._get(theta, "H0") / 100.0, self._get(theta, "ombh2"),
                                    w0, wa)
        return bg.Cosmo.from_Om(Om, np.full(Om.size, UNCALIBRATED_H0 / 100.0),
                                np.full(Om.size, UNCALIBRATED_OMBH2), w0, wa)

    def in_prior(self, theta: np.ndarray) -> np.ndarray:
        theta = np.atleast_2d(theta)
        ok = np.ones(theta.shape[0], bool)
        for j, n in enumerate(self.names):
            lo, hi = PRIORS[n]
            ok &= (theta[:, j] >= lo) & (theta[:, j] <= hi)
        if REQUIRE_W0_PLUS_WA_NEGATIVE and "wa" in self.names:
            ok &= self._get(theta, "w0") + self._get(theta, "wa") < 0.0
        if self.calibrated:
            c = self.cosmology(theta)
            ok &= self.early.inside(c.ombh2, c.omch2) & (c.omega_de > 0)
        return ok

    # --- chi^2 ----------------------------------------------------------------

    def chi2_parts(self, theta: np.ndarray, *, profile_hrd: bool = False) -> dict:
        """chi^2 of each likelihood for a batch of parameter vectors.

        With ``profile_hrd`` the vectors omit h*r_d, which is solved for exactly.
        Points outside the prior get inf.
        """
        theta = np.atleast_2d(np.asarray(theta, float))
        if profile_hrd and "hrd" in self.names:
            full = np.insert(theta, self.names.index("hrd"), 100.0, axis=1)
        else:
            full = theta
        n = full.shape[0]
        out = {"total": np.full(n, np.inf)}
        ok = self.in_prior(full)
        if not ok.any():
            return out
        th = full[ok]
        cosmo = self.cosmology(th)
        f = self.grid.integrand(cosmo)
        dc = self.grid.comoving_nodes(cosmo, f)
        dct = dc[:, self.grid.index]
        total = np.zeros(th.shape[0])
        parts = {}
        if self.bao is not None:
            zu = self.bao.unique_z
            dm = dct[:, self._pos["bao"]]
            hz = bg.hubble(zu, cosmo)
            unscaled = self.bao.model_unscaled(dm, hz)
            if self.calibrated:
                rd = self.early.rd(cosmo)
                parts["bao"] = self.bao.chi2(unscaled / rd[:, None])
                parts["rd"] = rd
            elif profile_hrd:
                parts["bao"], rd = self.bao.chi2_profile_scale(unscaled)
                parts["rd"] = rd
                parts["hrd"] = rd * cosmo.h
            else:
                rd = self._get(th, "hrd") / cosmo.h
                parts["bao"] = self.bao.chi2(unscaled / rd[:, None])
                parts["rd"] = rd
            total += parts["bao"]
        if self.sn is not None:
            parts["sn"], parts["sn_offset"] = self.sn.chi2(dct[:, self._pos["sn"]])
            total += parts["sn"]
        if self.cmb is not None:
            theta_star = self.early.theta_star(cosmo, self.grid, dc_nodes=dc, integrand=f)
            parts["theta_star"] = theta_star
            parts["cmb"] = self.cmb.chi2(theta_star, cosmo.ombh2, cosmo.omega_cb)
            total += parts["cmb"]
            parts.setdefault("rd", self.early.rd(cosmo))
        parts["total"] = total
        for k, v in parts.items():
            arr = np.full(n, np.inf if k in ("total", "bao", "sn", "cmb") else np.nan)
            arr[ok] = v
            out[k] = arr
        return out

    def log_prob(self, theta: np.ndarray) -> np.ndarray:
        """Log posterior (flat priors) for emcee, vectorised over walkers."""
        return -0.5 * self.chi2_parts(theta)["total"]


# ---------------------------------------------------------------------------
# Best fits
# ---------------------------------------------------------------------------

#: Typical scales, to give the simplex sensible shapes in every direction.
SCALE = {"H0": 1.0, "ombh2": 0.0002, "Om": 0.01, "hrd": 1.0, "w0": 0.05, "wa": 0.2}
#: Generic starting points (no published best fit is used to start).
STARTS = {"H0": (67.0, 64.0, 70.0), "ombh2": (0.0222, 0.0222, 0.0224),
          "Om": (0.31, 0.34, 0.29), "hrd": (100.0, 98.0, 102.0),
          "w0": (-1.0, -0.6, -0.85), "wa": (0.0, -1.5, -0.6)}


@dataclass
class MAPResult:
    names: list
    x: np.ndarray
    chi2: float
    parts: dict
    n_eval: int
    seconds: float

    def as_dict(self) -> dict:
        d = {n: float(v) for n, v in zip(self.names, self.x)}
        d.update({f"chi2_{k}": float(v) for k, v in self.parts.items()
                  if k in ("bao", "sn", "cmb")})
        for k in ("rd", "hrd", "theta_star", "sn_offset"):
            if k in self.parts and np.isfinite(self.parts[k]):
                d[k] = float(self.parts[k])
        d["chi2"] = float(self.chi2)
        return d


def find_map(problem: Problem, extra_starts=(), *, tol: float = 1e-7) -> MAPResult:
    """Minimise chi^2 (= -2 ln posterior + const) from several starting points."""
    profile = "hrd" in problem.names
    names = problem.free_names(profile)
    scale = np.array([SCALE[n] for n in names])
    t0 = time.monotonic()
    neval = 0

    def f(u):
        nonlocal neval
        neval += 1
        return float(problem.chi2_parts((u * scale)[None, :], profile_hrd=profile)["total"][0])

    starts = [np.array([STARTS[n][i] for n in names]) for i in range(3)]
    starts += [np.asarray(s, float) for s in extra_starts]
    best = None
    for x0 in starts:
        if not np.isfinite(f(x0 / scale)):
            continue
        res = optimize.minimize(f, x0 / scale, method="Nelder-Mead",
                                options={"xatol": tol, "fatol": tol, "maxiter": 20000,
                                         "maxfev": 40000, "adaptive": True})
        # restart once from the end point: Nelder-Mead can stall on long ridges
        res = optimize.minimize(f, res.x, method="Nelder-Mead",
                                options={"xatol": tol, "fatol": tol, "maxiter": 20000,
                                         "maxfev": 40000, "adaptive": True})
        if best is None or res.fun < best.fun:
            best = res
    x = best.x * scale
    parts = {k: float(v[0]) for k, v in
             problem.chi2_parts(x[None, :], profile_hrd=profile).items()}
    if profile:  # report the full parameter vector, with the solved-for h*r_d
        x = np.insert(x, problem.names.index("hrd"), parts["hrd"])
    return MAPResult(problem.names, x, parts["total"], parts, neval, time.monotonic() - t0)


def significance(dchi2: float, extra_params: int = 2) -> float:
    """Gaussian-equivalent sigma of a chi^2 improvement with k extra parameters.

    DESI's convention: p = P(chi^2_k > dchi2), expressed as a two-sided sigma.
    """
    p = stats.chi2.sf(max(dchi2, 0.0), extra_params)
    return float(stats.norm.isf(p / 2.0))


# ---------------------------------------------------------------------------
# Posterior sampling
# ---------------------------------------------------------------------------


def run_mcmc(problem: Problem, start: np.ndarray, *, nwalkers: int = 48, nsteps: int = 6000,
             seed: int = 1, burn_fraction: float = 0.3, scatter=None, progress: bool = False) -> dict:
    """Sample the posterior with emcee (affine-invariant ensemble), all walkers at once."""
    import emcee

    rng = np.random.default_rng(seed)
    ndim = problem.ndim
    if scatter is None:
        scatter = np.array([SCALE[n] for n in problem.names]) * 0.3
    p0 = []
    while len(p0) < nwalkers:  # a small ball around the best fit, inside the prior
        trial = start + scatter * rng.standard_normal((nwalkers, ndim))
        lp = problem.log_prob(trial)
        p0.extend(trial[np.isfinite(lp)])
    p0 = np.array(p0[:nwalkers])
    sampler = emcee.EnsembleSampler(nwalkers, ndim, problem.log_prob, vectorize=True)
    # emcee snapshots numpy's global RNG when built; seed its own state explicitly
    sampler.random_state = np.random.RandomState(seed).get_state()
    t0 = time.monotonic()
    sampler.run_mcmc(p0, nsteps, progress=progress)
    burn = int(burn_fraction * nsteps)
    tau = sampler.get_autocorr_time(discard=burn, tol=0)
    chain = sampler.get_chain(discard=burn, flat=True)
    lp = sampler.get_log_prob(discard=burn, flat=True)
    return {"names": problem.names, "samples": chain, "log_prob": lp, "tau": tau,
            "n_eff": chain.shape[0] / float(np.max(tau)),
            "acceptance": float(np.mean(sampler.acceptance_fraction)),
            "nsteps": nsteps, "nwalkers": nwalkers, "burn": burn,
            "seconds": time.monotonic() - t0}


def summarize(x: np.ndarray, w: np.ndarray | None = None) -> dict:
    """Mean, standard deviation and equal-tailed 68% interval (GetDist's two-tail)."""
    x = np.asarray(x, float)
    w = np.ones_like(x) if w is None else np.asarray(w, float)
    order = np.argsort(x)
    cdf = np.cumsum(w[order]) / w.sum()
    lo, med, hi = np.interp([0.158655, 0.5, 0.841345], cdf, x[order])
    mean = float(np.average(x, weights=w))
    std = float(np.sqrt(np.average((x - mean) ** 2, weights=w)))
    return {"mean": mean, "std": std, "median": float(med),
            "lo68": float(lo), "hi68": float(hi),
            "upper68": float(np.interp(0.68, cdf, x[order])),
            "lower68": float(np.interp(0.32, cdf, x[order]))}
