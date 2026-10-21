"""Fit our orbit model to one week of an analysis centre's precise orbit.

The weekly SP3 orbits are our "observations": positions every 2 minutes, which
the analysis centres fitted to the laser ranges with their own models. We rotate
them into the GCRS and fit our own dynamical model by Gauss-Newton least squares:

    parameters: the state (position, velocity) at the start of the week,
                along-track accelerations (constant, cos u, sin u), a constant
                radial acceleration, and optionally the radiation-pressure scale
                C_R. The adopted runs hold C_R at one value per satellite
                (fix_cr): fitted every week, it absorbs part of frame dragging
                (in eclipse seasons radiation pressure can turn the orbital
                plane), see test_physics.test_jump_gain.

Frame dragging is switched *off* in our model (mu_LT = 0). Over one week it
moves the orbital plane by only 0.6 mas (4 cm at the satellite), so the fit
simply averages over it; what we later measure is how far each new week's
orbital plane has moved beyond our model's prediction from the week before.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

import arcmodel as AM
import config as C
import forcemodel as FM
import frames as F

NOMINAL_CR = {"lageos1": 1.13, "lageos2": 1.12, "lares2": 1.08}


@dataclass
class ArcFit:
    sat: str
    ac: str
    week: str
    mjd0: float               # arc epoch (MJD TT) = first SP3 epoch
    span_s: float             # length of the fitted data (s)
    y0: np.ndarray            # fitted GCRS state at mjd0
    prm: np.ndarray           # fitted C_R, T0, Tc, Ts, R0
    cov: np.ndarray           # formal covariance (unit weight 1 cm), scaled by the residual rms
    rms_rtn: np.ndarray       # residual rms radial, along, cross (m)
    n_obs: int
    n_iter: int
    converged: bool


def observations(orbit: F.Orbit, mjd0: float, step_s: float = C.FIT_SAMPLE_S):
    """GCRS positions at the fit epochs (every step_s), as (t_since_mjd0, r)."""
    t = (orbit.mjd_tt - mjd0) * C.SECONDS_PER_DAY
    keep = np.isclose(np.mod(t + 1e-6, step_s), 0.0, atol=1e-3)
    r, v = F.itrs_to_gcrs(orbit.mjd_utc[keep], orbit.r_itrs[keep], orbit.v_itrs[keep])
    return t[keep], r, v


def rtn_basis(r, v):
    R = r / np.linalg.norm(r, axis=1)[:, None]
    N = np.cross(r, v)
    N /= np.linalg.norm(N, axis=1)[:, None]
    T = np.cross(N, R)
    return R, T, N


def fit_arc(orbit: F.Orbit, model: AM.ArcModel | None = None, *, extra_days: float = 1.0,
            max_iter: int = 10, fit_mask=None, h: float = C.INTEGRATION_STEP_S,
            step_s: float = C.FIT_SAMPLE_S, model_kw: dict | None = None,
            fix_cr: bool = False, cr: float | None = None):
    """Fit one weekly orbit. Returns (ArcFit, model) - the model covers the arc + extra_days.

    fix_cr: hold the radiation-pressure scale at its nominal value (or ``cr``).
    """
    if fix_cr:
        fit_mask = np.ones(6 + FM.NPAR, bool)
        fit_mask[6] = False
    if isinstance(cr, dict):
        cr = cr.get(orbit.sat)
    mjd0 = float(orbit.mjd_tt[0])
    t_obs, r_obs, v_obs = observations(orbit, mjd0, step_s)
    span = float(t_obs[-1])
    if model is None:
        model = AM.build(mjd0, span / C.SECONDS_PER_DAY + extra_days, orbit.sat, **(model_kw or {}))
    stride = int(round(step_s / h))
    y0 = np.concatenate([r_obs[0], v_obs[0]])
    prm = np.array([NOMINAL_CR.get(orbit.sat, 1.1) if cr is None else cr, 0.0, 0.0, 0.0, 0.0])
    mask = np.ones(6 + FM.NPAR, bool) if fit_mask is None else np.asarray(fit_mask, bool)
    converged = False
    sigma = 0.01
    for it in range(1, max_iter + 1):
        T, Y = model.run(y0, prm, 0.0, span, h=h, stride=stride, with_var=True)
        idx = np.searchsorted(T, t_obs - 1e-6)
        if np.any(np.abs(T[idx] - t_obs) > 1e-3):
            raise RuntimeError("output epochs do not match the observation epochs")
        Ym = Y[idx]
        res = (r_obs - Ym[:, :3]).ravel()
        ncol = 6 + FM.NPAR
        H = np.empty((res.size, ncol))
        for col in range(ncol):
            H[:, col] = Ym[:, 6 + col * 6: 6 + col * 6 + 3].ravel()
        Hm = H[:, mask] / sigma
        # scale columns for conditioning
        scale = np.linalg.norm(Hm, axis=0)
        A = Hm / scale
        dx, *_ = np.linalg.lstsq(A, res / sigma, rcond=None)
        dx = dx / scale
        step = np.zeros(ncol)
        step[mask] = dx
        y0 = y0 + step[:6]
        prm = prm + step[6:]
        if np.all(np.abs(step[:3]) < 1e-4) and np.all(np.abs(step[3:6]) < 1e-7):
            converged = True
            break
    # final residuals and covariance
    T, Y = model.run(y0, prm, 0.0, span, h=h, stride=stride, with_var=True)
    idx = np.searchsorted(T, t_obs - 1e-6)
    Ym = Y[idx]
    d = r_obs - Ym[:, :3]
    R, Tt, N = rtn_basis(r_obs, v_obs)
    rms = np.array([np.sqrt(np.mean(np.einsum("ij,ij->i", d, B) ** 2)) for B in (R, Tt, N)])
    ncol = 6 + FM.NPAR
    H = np.empty((d.size, ncol))
    for col in range(ncol):
        H[:, col] = Ym[:, 6 + col * 6: 6 + col * 6 + 3].ravel()
    Hm = H[:, mask]
    s2 = np.mean(d.ravel() ** 2) * d.size / max(1, d.size - mask.sum())
    cov = np.zeros((ncol, ncol))
    cov[np.ix_(mask, mask)] = np.linalg.pinv(Hm.T @ Hm) * s2
    converged = converged and bool(np.all(rms < 1.0))     # metre-level residuals: not a solution
    fit = ArcFit(sat=orbit.sat, ac=orbit.ac, week=orbit.week, mjd0=mjd0, span_s=span, y0=y0,
                 prm=prm, cov=cov, rms_rtn=rms, n_obs=len(t_obs), n_iter=it, converged=converged)
    return fit, model
