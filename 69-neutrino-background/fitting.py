"""Fitting universes to data: Gauss–Newton fits and profile likelihoods.

Every likelihood here is (close to) Gaussian, so each one can be written as a
vector of *whitened residuals* r with χ² = |r|². A data combination simply
concatenates its likelihoods' residual vectors, and fitting is non-linear least
squares.

* :func:`fit` — Gauss–Newton: linearise the model around the current guess
  (the Jacobian J = ∂r/∂θ by central differences, all 2n CAMB runs in
  parallel), solve the linear least-squares step, repeat. The CMB is so nearly
  linear in (ω_b, ω_c, θ, ln A_s, n_s, τ) that it converges in two or three
  steps; the covariance (JᵀJ)⁻¹ gives Gaussian ("Laplace") error bars.

* :func:`profile` — the profile likelihood of one parameter (N_eff or Σm_ν):
  fix it on a grid and minimise χ² over all the others at each grid value.
  Δχ²(value) is the honest summary when the parameter sits against a physical
  boundary (Σm_ν ≥ 0). Each grid point starts from the linear prediction of the
  full fit and iterates with that fit's Jacobian, refined by Broyden secant
  updates, so a point costs a handful of CAMB runs instead of a new Jacobian
  each step. Points run in parallel, one per worker.

For a Gaussian likelihood the profile and the Bayesian marginal coincide; the
Bayesian limits we quote integrate exp(−Δχ²/2) over the prior range.
"""

from __future__ import annotations

import multiprocessing as mp
import time
from dataclasses import dataclass, field

import numpy as np

import cmb_model as cm
import config as cfg


@dataclass
class Problem:
    """A data combination and the parameters it varies.

    ``likelihoods``: objects with ``name``, ``residual(theory, p)`` and an
    optional ``nuisance`` tuple of parameter names they own (for example
    Planck's calibration A_planck), which never require a new CAMB run.
    """

    name: str
    likelihoods: list
    free: list[str]
    nu_model: str = "one"
    fixed: dict = field(default_factory=dict)

    @property
    def nuisance_names(self) -> set:
        out = set()
        for lk in self.likelihoods:
            out.update(getattr(lk, "nuisance", ()))
        return out

    def is_camb(self, name: str) -> bool:
        return name not in self.nuisance_names

    def params(self, x: dict) -> dict:
        return {**cm.DEFAULTS, **self.fixed, **x}

    def residual(self, theory: dict, p: dict) -> np.ndarray:
        return np.concatenate([lk.residual(theory, p) for lk in self.likelihoods])

    def chi2_parts(self, theory: dict, p: dict) -> dict:
        """χ² of each likelihood in its published convention (constants restored)."""
        return {lk.name: float(np.sum(lk.residual(theory, p) ** 2))
                + float(getattr(lk, "chi2_constant", 0.0)) for lk in self.likelihoods}


def _clip(name: str, value: float) -> float:
    lo, hi = cfg.PARAM_BOUNDS.get(name, (-np.inf, np.inf))
    return float(min(max(value, lo), hi))


def _camb_part(p: dict) -> dict:
    return {k: p[k] for k in cm.CAMB_KEYS if k in p}


# ---------------------------------------------------------------------------
# Gauss–Newton
# ---------------------------------------------------------------------------


def _evaluate(problem: Problem, plist: list[dict], pool) -> list[tuple[np.ndarray, dict]]:
    theories = pool.map([_camb_part(p) for p in plist], problem.nu_model)
    return [(problem.residual(t, p), t) for t, p in zip(theories, plist)]


def jacobian(problem: Problem, p: dict, free: list[str], pool) -> tuple[np.ndarray, np.ndarray, dict]:
    """Central-difference Jacobian of the residuals; returns (J, r0, theory0)."""
    plist = [p]
    index = []
    for k in free:
        h = cfg.STEPS[k]
        up, dn = dict(p), dict(p)
        up[k] = p[k] + h
        dn[k] = p[k] - h
        lo, hi = cfg.PARAM_BOUNDS.get(k, (-np.inf, np.inf))
        if dn[k] < lo:  # one-sided at a boundary (Σm_ν near 0)
            dn[k] = p[k]
        if up[k] > hi:
            up[k] = p[k]
        index.append((len(plist), len(plist) + 1, up[k] - dn[k]))
        plist += [up, dn]
    evals = _evaluate(problem, plist, pool)
    r0, t0 = evals[0]
    cols = [(evals[i][0] - evals[j][0]) / d for i, j, d in index]
    return np.column_stack(cols), r0, t0


def fit(problem: Problem, start: dict, pool, *, max_iter: int = 10, tol: float = 0.005,
        verbose: bool = True) -> dict:
    """Minimise χ² over ``problem.free`` starting from ``start``."""
    free = list(problem.free)
    p = problem.params(start)
    for k in free:
        p[k] = _clip(k, p[k])
    t_start = time.monotonic()
    J, r, theory = jacobian(problem, p, free, pool)
    chi2 = float(r @ r)
    history = [chi2]
    if verbose:
        print(f"    [{problem.name}] start chi2 = {chi2:.2f}")
    for it in range(max_iter):
        step = np.linalg.lstsq(J, -r, rcond=None)[0]
        lam, accepted = 1.0, False
        while lam >= 1 / 16:
            trial = dict(p)
            for k, s in zip(free, step):
                trial[k] = _clip(k, p[k] + lam * s)
            (r_t, t_t), = _evaluate(problem, [trial], pool)
            chi2_t = float(r_t @ r_t)
            if chi2_t <= chi2 + 1e-9:
                accepted = True
                break
            lam /= 2
        if not accepted:
            break
        improvement = chi2 - chi2_t
        p, r, theory, chi2 = trial, r_t, t_t, chi2_t
        history.append(chi2)
        if verbose:
            print(f"      iteration {it + 1}: chi2 = {chi2:.3f} (step x{lam:g}, "
                  f"{time.monotonic() - t_start:.0f} s)")
        if improvement < tol:
            break
        J, r, theory = jacobian(problem, p, free, pool)
    # Final Jacobian at the best fit: the covariance and the profile's linear model.
    J, r, theory = jacobian(problem, p, free, pool)
    chi2 = float(r @ r)
    cov = np.linalg.pinv(J.T @ J)
    err = np.sqrt(np.diag(cov))
    return {
        "problem": problem.name,
        "free": free,
        "best": {k: float(p[k]) for k in free},
        "params": {k: float(v) for k, v in p.items() if isinstance(v, (int, float))},
        "errors": {k: float(e) for k, e in zip(free, err)},
        "cov": cov.tolist(),
        "jacobian": J.tolist(),
        "chi2": chi2,
        "chi2_parts": problem.chi2_parts(theory, p),
        "n_data": int(r.size),
        "derived": cm.derived_dict(theory),
        "history": history,
        "seconds": time.monotonic() - t_start,
    }


# ---------------------------------------------------------------------------
# Profile likelihood
# ---------------------------------------------------------------------------


def _local_jacobian(problem: Problem, p: dict, names: list[str], resid) -> np.ndarray:
    """Central-difference Jacobian evaluated in this process (a profile worker)."""
    cols = []
    for k in names:
        h = cfg.STEPS[k]
        up, dn = dict(p), dict(p)
        up[k] = _clip(k, p[k] + h)
        dn[k] = _clip(k, p[k] - h)
        cols.append((resid(up)[0] - resid(dn)[0]) / (up[k] - dn[k]))
    return np.column_stack(cols)


def _profile_point(args):
    """One grid point: minimise over the other parameters.

    Gauss-Newton steps with the full fit's Jacobian, refined by Broyden secant
    updates. Converged when the linear model predicts less than ``tol`` of
    further improvement. If a step fails to lower χ² even when halved three
    times, the Jacobian is recomputed here (2n CAMB runs) and the step retried;
    a second failure means the point sits at the minimum to within CAMB's
    numerical noise.
    """
    problem, name, value, start, J0, others, max_iter, tol = args
    p = dict(start)
    p[name] = value
    calls = 0

    def resid(q):
        nonlocal calls
        calls += 1
        t = cm.load_cached(_camb_part(q), problem.nu_model)
        if t is None:
            t = cm.compute_theory(_camb_part(q), problem.nu_model)
            cm.save_cached(_camb_part(q), problem.nu_model, t)
        return problem.residual(t, q), t

    r, t = resid(p)
    chi2 = float(r @ r)
    J = np.array(J0, float)
    history = [chi2]
    status = "max_iter"
    refreshed = False
    for _ in range(max_iter):
        step = np.linalg.lstsq(J, -r, rcond=None)[0]
        predicted = chi2 - float(np.sum((r + J @ step) ** 2))
        if predicted < tol:
            status = "converged"
            break
        lam, accepted = 1.0, False
        while lam >= 1 / 8:
            q = dict(p)
            for k, s_ in zip(others, step):
                q[k] = _clip(k, p[k] + lam * s_)
            r_q, t_q = resid(q)
            chi2_q = float(r_q @ r_q)
            if chi2_q < chi2:
                accepted = True
                break
            lam /= 2
        if not accepted:
            if refreshed:
                status = "stalled"  # at the minimum to within numerical noise
                break
            J = _local_jacobian(problem, p, others, resid)
            refreshed = True
            continue
        dx = np.array([q[k] - p[k] for k in others])
        dr = r_q - r
        denom = dx @ dx
        if denom > 0:  # Broyden secant update of the Jacobian
            J = J + np.outer(dr - J @ dx, dx) / denom
        improvement = chi2 - chi2_q
        p, r, t, chi2 = q, r_q, t_q, chi2_q
        history.append(chi2)
        if improvement < tol:
            status = "converged"
            break
    return {"value": float(value), "chi2": chi2, "params": {k: float(p[k]) for k in others},
            "chi2_parts": problem.chi2_parts(t, p), "derived": cm.derived_dict(t),
            "calls": calls, "converged": status != "max_iter", "status": status,
            "refreshed_jacobian": refreshed, "history": history}


def profile(problem: Problem, best: dict, name: str, grid, pool, *, max_iter: int = 12,
            tol: float = 0.003, verbose: bool = True) -> dict:
    """Δχ² profile of ``name`` over ``grid``, starting from a full fit ``best``.

    ``best`` must come from :func:`fit` with ``name`` among the free
    parameters, so that its Jacobian and covariance give each grid point a
    linear starting guess (the conditional mean of the others given ``name``).
    """
    free = best["free"]
    if name not in free:
        raise ValueError(f"{name} was not free in the base fit")
    others = [k for k in free if k != name]
    i = free.index(name)
    cov = np.array(best["cov"])
    J = np.array(best["jacobian"])
    J_others = J[:, [free.index(k) for k in others]]
    slope = cov[:, i] / cov[i, i]
    p0 = problem.params(best["params"])
    t0 = time.monotonic()
    tasks = []
    for v in grid:
        start = dict(p0)
        for k, s in zip(free, slope):
            if k != name:
                start[k] = _clip(k, p0[k] + s * (v - p0[name]))
        tasks.append((problem, name, float(v), start, J_others, others, max_iter, tol))
    ctx = mp.get_context("fork")
    if pool.n_workers > 1 and len(tasks) > 1:
        with ctx.Pool(min(pool.n_workers, len(tasks)), initializer=cm._init_worker) as wp:
            points = wp.map(_profile_point, tasks, chunksize=1)
    else:
        points = [_profile_point(a) for a in tasks]
    chi2_min = min(min(pt["chi2"] for pt in points), best["chi2"])
    for pt in points:
        pt["dchi2"] = pt["chi2"] - chi2_min
    if verbose:
        calls = sum(pt["calls"] for pt in points)
        nc = sum(not pt["converged"] for pt in points)
        nst = sum(pt.get("status") == "stalled" for pt in points)
        print(f"    [{problem.name}] profile of {name}: {len(points)} points, {calls} CAMB "
              f"runs, {time.monotonic() - t0:.0f} s" + (f", {nc} hit max_iter" if nc else "")
              + (f", {nst} stalled at noise level" if nst else ""))
    return {"problem": problem.name, "name": name, "points": points, "best": best["best"],
            "chi2_best": best["chi2"], "seconds": time.monotonic() - t0}


# ---------------------------------------------------------------------------
# Turning a profile into numbers
# ---------------------------------------------------------------------------


def parabola_fit(values, dchi2, *, max_dchi2: float = 9.0) -> dict:
    """Fit Δχ² = ((x − μ)/σ)² + c to the profile (points with Δχ² < max_dchi2).

    Extrapolates the minimum beyond a physical boundary: for Σm_ν this is the
    "effective" mass the data would prefer if negative values were allowed.
    """
    x = np.asarray(values, float)
    y = np.asarray(dchi2, float)
    m = y < max_dchi2
    if m.sum() < 3:
        m = np.argsort(y)[:3]
    a, b, c = np.polyfit(x[m], y[m], 2)
    sigma = 1.0 / np.sqrt(a) if a > 0 else np.inf
    mu = -b / (2 * a) if a > 0 else np.nan
    resid = y[m] - np.polyval([a, b, c], x[m])
    return {"mu": float(mu), "sigma": float(sigma), "offset": float(c - b**2 / (4 * a)),
            "rms_resid": float(np.sqrt(np.mean(resid**2)))}


def bayes_upper_limit(values, dchi2, *, lower: float = 0.0, level: float = 0.95,
                      n: int = 20001) -> float:
    """Upper limit at ``level`` for a flat prior above ``lower``, from a profile.

    The likelihood exp(−Δχ²/2) is interpolated in Δχ² (cubic in the value,
    linear extrapolation past the last grid point is not allowed).
    """
    from scipy.interpolate import PchipInterpolator

    x = np.asarray(values, float)
    y = np.asarray(dchi2, float)
    order = np.argsort(x)
    x, y = x[order], y[order]
    f = PchipInterpolator(x, y, extrapolate=False)
    grid = np.linspace(max(lower, x[0]), x[-1], n)
    like = np.exp(-0.5 * (f(grid) - np.nanmin(f(grid))))
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (like[1:] + like[:-1]) * np.diff(grid))])
    if like[-1] > 1e-3 * like.max():
        raise ValueError("profile grid does not reach far enough to bound the integral")
    cdf /= cdf[-1]
    return float(np.interp(level, cdf, grid))


def truncated_gaussian_limit(mu: float, sigma: float, *, lower: float = 0.0,
                             level: float = 0.95) -> float:
    """Upper limit for a Gaussian likelihood N(μ, σ) and a flat prior x > lower."""
    from scipy.stats import norm

    a = norm.cdf((lower - mu) / sigma)
    return float(mu + sigma * norm.ppf(a + level * (1 - a)))


def gaussian_interval(values, dchi2) -> dict:
    """Best value and 68% interval (Δχ² = 1 crossings) of a profile."""
    from scipy.interpolate import PchipInterpolator

    x = np.asarray(values, float)
    y = np.asarray(dchi2, float)
    order = np.argsort(x)
    x, y = x[order], y[order]
    fine = np.linspace(x[0], x[-1], 20001)
    yy = PchipInterpolator(x, y)(fine)
    i = int(np.argmin(yy))
    best = fine[i]
    lo = np.interp(1.0, yy[:i + 1][::-1], fine[:i + 1][::-1]) if yy[0] > 1 else np.nan
    hi = np.interp(1.0, yy[i:], fine[i:]) if yy[-1] > 1 else np.nan
    return {"best": float(best), "lo": float(lo), "hi": float(hi),
            "minus": float(best - lo), "plus": float(hi - best)}
