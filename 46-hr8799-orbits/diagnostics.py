"""Sanity checks of the nested-sampling fits.

A nested-sampling run can lose the best-fitting region of a rugged, many-peaked
likelihood early on and then converge confidently on a worse one. Two checks
catch this:

1. Nesting. The models are nested: near1248_coplanar is a special case of both
   coplanar and near1248, which are special cases of unconstrained. Any point
   of a constrained model is a valid point of the general one, so the general
   model's best likelihood can never be lower. We map each fit's best samples
   into every model that contains it.
2. Local optimisation. From the best samples (own and mapped), a bounded
   Powell search finds the best likelihood nearby. A posterior whose best
   sample lies far below that maximum (more than ~10 in ln L; a healthy
   20-26-parameter posterior sits ~3-5 below it) has missed the peak.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

import config as C
import kepler as K
import models as M

#: model -> the models that contain it
CONTAINED_IN = {"near1248_coplanar": ("coplanar", "near1248", "unconstrained"),
                "coplanar": ("unconstrained",), "near1248": ("unconstrained",), "unconstrained": ()}


def theta_from_physical(model: M.OrbitModel, phys: dict) -> np.ndarray | None:
    """Inverse of model.physical for one parameter set (None if not representable).

    Elements -> each planet's Jacobi sky state at model.t_state. The mirror orbit
    (z -> -z) is used where needed: it has the same likelihood.
    """
    el = phys["elements"]
    pos = {n: k for k, n in enumerate(model.labels)}
    th = np.zeros(model.ndim)
    plx = float(phys["plx"])
    mstar = float(phys["mstar"])
    masses = {p: float(phys["masses"][p]) for p in C.INNER_TO_OUTER} if model.fit_masses else \
        {p: 0.0 for p in C.INNER_TO_OUTER}
    m_inside = mstar
    for p in C.INNER_TO_OUTER:
        a, e, inc, aop, pan, tau = map(float, el[p])
        m_k = m_inside + masses[p]
        (x, y, z), (vx, vy, vz) = K.orbit_xyz(model.t_state, a, e, inc, aop, pan, tau, m_k, velocity=True)
        st = np.array([x, y, z, vx, vy, vz]) * plx
        if model.coplanar:
            ae, ee, ie, we, Oe, te = map(float, el["e"])
            if abs(inc - ie) > 1e-9 or abs(np.mod(pan - Oe + np.pi / 2, np.pi) - np.pi / 2) > 1e-9:
                return None
            if np.cos(ie) <= 0:
                return None
            th[pos["cos_i"]] = np.cos(ie)
            th[pos["Omega"]] = np.mod(Oe, np.pi)
        else:
            if st[2] < 0:
                st[2], st[5] = -st[2], -st[5]
            th[pos[f"z_{p}"]] = st[2]
            th[pos[f"vz_{p}"]] = st[5]
        for q, v in zip(("x", "y", "vx", "vy"), st[[0, 1, 3, 4]]):
            th[pos[f"{q}_{p}"]] = v
        if model.fit_masses:
            th[pos[f"m_{p}"]] = masses[p] / C.MJUP_MSUN
        m_inside = m_k
    th[pos["plx"]] = plx
    th[pos["mstar"]] = mstar
    for k, s in enumerate(model.data.jitter_sources):
        th[pos[f"log10_jit_{s}"]] = np.log10(phys["jitter"][k]) if len(phys["jitter"]) else -2.0
    return th


def _single(phys: dict, k: int) -> dict:
    """The k-th parameter set of a vectorised physical dict."""
    out = {"elements": {p: tuple(np.atleast_1d(v)[k] for v in el) for p, el in phys["elements"].items()},
           "plx": np.atleast_1d(phys["plx"])[k], "mstar": np.atleast_1d(phys["mstar"])[k],
           "masses": {p: np.atleast_1d(v)[k] if np.ndim(v) else v for p, v in phys["masses"].items()},
           "jitter": np.atleast_2d(phys["jitter"])[:, k] if np.size(phys["jitter"]) else np.zeros(0)}
    return out


def polish(model: M.OrbitModel, theta0: np.ndarray, maxfev: int = 30000) -> tuple[np.ndarray, float]:
    """Local maximisation of the data likelihood (no prior terms) within the
    model's prior support, by adaptive Nelder-Mead (never returns worse than
    its start; the support's walls are steep jumps that confuse line searches)."""

    def f(th):
        if model.loglike(th) < -1e11:          # outside the prior's support
            return 1e15
        return -model.data_loglike(th)

    r = minimize(f, theta0, method="Nelder-Mead",
                 options=dict(maxfev=maxfev, xatol=1e-8, fatol=1e-6, adaptive=True))
    return r.x, -r.fun


def check(fit_objects: dict, n_top: int = 6, n_scan: int = 3000) -> dict:
    """fit_objects: key -> fits.Fit (same data set). Returns a report per key.

    Compares data log-likelihoods only (the models have different priors)."""
    report = {}
    best_points = {}
    rng = np.random.default_rng(C.RANDOM_SEED)
    for key, fit in fit_objects.items():
        m = fit.model
        cand = np.vstack([fit.raw["samples"][np.argsort(-fit.raw["logl"])[:200]],
                          fit.equal[rng.choice(len(fit.equal), min(n_scan, len(fit.equal)), replace=False)]])
        lnl = np.array([m.data_loglike(th) for th in cand])
        order = np.argsort(-lnl)
        best_points[key] = cand[order[:n_top]]
        report[key] = {"posterior_max_lnL": float(lnl[order[0]]), "polished_max_lnL": -np.inf,
                       "source_of_best": None}
    kinds = {key: fit.model.kind for key, fit in fit_objects.items()}
    for key, fit in fit_objects.items():
        candidates = [(key, th) for th in best_points[key]]
        for other, okind in kinds.items():
            if kinds[key] in CONTAINED_IN.get(okind, ()):
                ofit = fit_objects[other]
                for th in best_points[other]:
                    ph = _single(ofit.model.physical(th[None, :]), 0)
                    mapped = theta_from_physical(fit.model, ph)
                    if mapped is not None and fit.model.loglike(mapped) > -1e11:
                        candidates.append((other, mapped))
        for src, th in candidates:
            x, val = polish(fit.model, th)
            if val > report[key]["polished_max_lnL"]:
                report[key].update(polished_max_lnL=float(val), source_of_best=src, best_theta=x)
        r = report[key]
        r["shortfall"] = r["polished_max_lnL"] - r["posterior_max_lnL"]
        r["ok"] = r["shortfall"] < 10.0
    return report


def box_margins(fit) -> dict:
    """How close each posterior comes to the edges of the sampled box.

    For every sky-state coordinate (cos i stops at the physical face-on bound;
    Omega is periodic) the smallest distance from a posterior sample (0.1 and 99.9
    percentiles) to the box edge, in units of the posterior's own standard
    deviation. Values > 5 mean the box did not truncate the posterior.
    z and vz are drawn inside their physical (bound-orbit) limits, so only
    their lower bound z >= 0 is a sampling edge."""
    m, s = fit.model, fit.equal
    out = {}
    for k, name in enumerate(m.labels):
        q = name.split("_")[0]
        if q not in ("x", "y", "vx", "vy"):
            continue
        lo, hi = min(m._lo[k], m._hi[k]), max(m._lo[k], m._hi[k])
        p_lo, p_hi = np.percentile(s[:, k], [0.1, 99.9])
        sd = np.std(s[:, k])
        out[name] = float(min(p_lo - lo, hi - p_hi) / sd) if sd > 0 else np.inf
    return out
