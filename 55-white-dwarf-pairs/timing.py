"""The ephemeris of ZTF J1539+5027 from survey photometry: P, P-dot and the O-C diagram.

Model. The orbital phase (in cycles, 0 = mid primary eclipse) at barycentric time t is

    phi(t) = phi0 + f0 x + f1 x^2 / 2 [+ f2 x^3 / 6],    x = t - T_REF (days),

so the orbital frequency f(t) = f0 + f1 x drifts linearly when the period derivative
is constant (f1 = -P-dot / P^2). In each band the flux is

    F(t) = c_group + sum_k a_k cos(2 pi k phi) [+ b_k sin(2 pi k phi)]
           - d1 E1(phi) - d2 E2(phi),

a smooth part (irradiated face of the cool star, tidal distortion: cosines, symmetric
about conjunction) and the two eclipses. E1 is the fraction of the hot primary's
limb-darkened disc hidden by the secondary, E2 the (normalised) fraction of the
secondary's disc hidden when the primary passes in front, both computed from the
geometry that Burdge et al. (2019) fitted to high-speed photometry (R1/a, R2/a and
the inclination, held fixed: survey data cannot measure them). c_group is a constant
per band and observing season (the difference images are relative to a reference
image whose own flux of the star is unknown). Every basis function is averaged over
the exposure (30 s for ZTF = 7% of an orbit; 60 s for PTF), with 15-point
Gauss-Legendre quadrature.

The model is linear in everything except (phi0, f0, f1), so those three are fitted
with every linear coefficient solved exactly at each step (variable projection).

Note: phi(t) counts orbits forward, so a shrinking period (f1 > 0) makes eclipses
arrive earlier and earlier than a constant-period clock predicts: the O-C diagram
bends down, as in Burdge et al.'s figure 2 and C9's Hulse-Taylor plot.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import optimize

import config as C

DAY = 86400.0
#: Reference epoch for the phase polynomial (BJD_TDB), near the middle of ZTF's span.
T_REF = 2459500.0
#: Burdge et al. (2019) geometry, fixed: radii in units of the separation, inclination.
R1_A = 1.562e-2 / 11.218e-2
R2_A = 3.140e-2 / 11.218e-2
INCL_DEG = 84.15
#: Linear limb darkening of the hot primary (Burdge et al.'s prior centre for g').
U_LD_PRIMARY = 0.15
N_NODES = 15
GL_X, GL_W = np.polynomial.legendre.leggauss(N_NODES)


# ---------------------------------------------------------------------------
# Eclipse geometry
# ---------------------------------------------------------------------------


def blocked_fraction(delta, r_star, r_occ, u_ld=0.0, n_r=600):
    """Fraction of a linearly limb-darkened disc (radius r_star) hidden by an opaque
    disc (radius r_occ) whose centre is ``delta`` away (all in the same units)."""
    delta = np.atleast_1d(np.asarray(delta, float))
    r = (np.arange(n_r) + 0.5) / n_r * r_star
    mu = np.sqrt(np.clip(1 - (r / r_star) ** 2, 0, None))
    intensity = 1 - u_ld * (1 - mu)
    w = intensity * 2 * np.pi * r * (r_star / n_r)
    d = delta[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        cos_t = (r[None, :] ** 2 + d**2 - r_occ**2) / (2 * r[None, :] * d)
    frac = np.where(cos_t <= -1, 1.0, np.where(cos_t >= 1, 0.0, np.arccos(np.clip(cos_t, -1, 1)) / np.pi))
    frac = np.where(d == 0, (r[None, :] < r_occ).astype(float), frac)
    return (frac * w).sum(axis=1) / w.sum()


class EclipseTables:
    """E1, E2 tabulated against projected separation, for fast interpolation."""

    def __init__(self, r1=R1_A, r2=R2_A, incl_deg=INCL_DEG, u1=U_LD_PRIMARY):
        self.r1, self.r2 = r1, r2
        self.cos_i = np.cos(np.radians(incl_deg))
        self.grid = np.linspace(0, r1 + r2 + 1e-6, 3001)
        self.e1 = blocked_fraction(self.grid, r1, r2, u1)
        e2 = blocked_fraction(self.grid, r2, r1, 0.0)
        self.e2 = e2 / e2.max()

    def separation(self, phi):
        ang = 2 * np.pi * phi
        return np.sqrt(np.sin(ang) ** 2 + (self.cos_i * np.cos(ang)) ** 2), np.cos(ang)

    def e1_e2(self, phi):
        sep, c = self.separation(phi)
        inside = sep < self.grid[-1]
        e1 = np.zeros_like(phi)
        e2 = np.zeros_like(phi)
        front_sec = inside & (c > 0)  # secondary in front of the primary: primary eclipse
        front_pri = inside & (c < 0)
        e1[front_sec] = np.interp(sep[front_sec], self.grid, self.e1)
        e2[front_pri] = np.interp(sep[front_pri], self.grid, self.e2)
        return e1, e2


ECL = EclipseTables()


# ---------------------------------------------------------------------------
# Data container and basis
# ---------------------------------------------------------------------------


@dataclass
class BandData:
    """Photometry of one band, with the quadrature times of each exposure."""

    name: str
    t: np.ndarray  # BJD_TDB of mid-exposure
    exptime_s: np.ndarray
    y: np.ndarray  # flux (any linear unit)
    sig: np.ndarray
    groups: np.ndarray  # integer group index per point
    group_names: list
    n_harm: int = 4
    sines: bool = False
    eclipses: bool = True
    xnodes: np.ndarray = field(init=False)

    def __post_init__(self):
        half = self.exptime_s / 2 / DAY
        self.xnodes = (self.t - T_REF)[:, None] + half[:, None] * GL_X[None, :]
        self.G = np.zeros((self.t.size, len(self.group_names)))
        self.G[np.arange(self.t.size), self.groups] = 1.0

    def subset(self, mask):
        return BandData(self.name, self.t[mask], self.exptime_s[mask], self.y[mask], self.sig[mask],
                        self.groups[mask], self.group_names, self.n_harm, self.sines, self.eclipses)


def phase(x, theta):
    """Orbital phase (cycles) at x = t - T_REF days for theta = (phi0, f0, f1[, f2])."""
    phi = theta[0] + theta[1] * x + 0.5 * theta[2] * x**2
    if len(theta) > 3:
        phi = phi + theta[3] * x**3 / 6
    return phi


def shape_basis(phi, n_harm, sines, eclipses):
    cols = []
    ang = 2 * np.pi * phi
    for k in range(1, n_harm + 1):
        cols.append(np.cos(k * ang))
        if sines:
            cols.append(np.sin(k * ang))
    if eclipses:
        e1, e2 = ECL.e1_e2(np.mod(phi, 1.0))
        cols += [-e1, -e2]
    return np.stack(cols, axis=-1)


def shape_names(n_harm, sines, eclipses):
    names = []
    for k in range(1, n_harm + 1):
        names.append(f"cos{k}")
        if sines:
            names.append(f"sin{k}")
    if eclipses:
        names += ["d1_primary", "d2_secondary"]
    return names


def smeared_basis(band: BandData, theta):
    phi = phase(band.xnodes, theta)  # (N, nodes)
    b = shape_basis(phi, band.n_harm, band.sines, band.eclipses)  # (N, nodes, k)
    return np.tensordot(b, GL_W / 2, axes=([1], [0]))


def solve_band(band: BandData, theta, shape_coef=None):
    """Weighted linear least squares for one band at fixed theta.

    If ``shape_coef`` is given the shape is held fixed (only group offsets and an
    overall amplitude scale are fitted), as used for per-season timing."""
    B = smeared_basis(band, theta)
    w = 1.0 / band.sig
    if shape_coef is None:
        X = np.hstack([band.G, B])
    else:
        X = np.hstack([band.G, (B @ shape_coef)[:, None]])
    beta, *_ = np.linalg.lstsq(X * w[:, None], band.y * w, rcond=None)
    model = X @ beta
    return beta, model, (band.y - model) * w


def residuals(theta, bands, shape_coefs=None):
    out = []
    for i, b in enumerate(bands):
        sc = None if shape_coefs is None else shape_coefs[i]
        out.append(solve_band(b, theta, sc)[2])
    return np.concatenate(out)


def chi2(theta, bands, shape_coefs=None):
    r = residuals(theta, bands, shape_coefs)
    return float(r @ r)


# ---------------------------------------------------------------------------
# Fitting
# ---------------------------------------------------------------------------


#: Typical 1-sigma sizes of (phi0, f0, f1, f2) in cycles, 1/d, 1/d^2, 1/d^3. The fit works
#: in units of these, so finite-difference steps (relative to max(1, |parameter|)) suit
#: every parameter; in raw units a step in f1 or f2 would shift the phase by many cycles.
PARAM_SCALE = np.array([1e-3, 1e-6, 1e-9, 1e-12])


def fit_ephemeris(bands, theta0, x_scale=None, shape_coefs=None, residual_fn=None, n_lin=None):
    """Least-squares (phi0, f0, f1[, f2]) with every linear coefficient profiled out.

    Returns the best theta, its covariance (scaled by the reduced chi^2 if > 1), chi^2, dof."""
    theta0 = np.asarray(theta0, float)
    scale = PARAM_SCALE[: theta0.size] if x_scale is None else np.asarray(x_scale, float)
    rfun = residual_fn if residual_fn is not None else (lambda th: residuals(th, bands, shape_coefs))
    res = optimize.least_squares(lambda u: rfun(theta0 + u * scale), np.zeros_like(theta0), method="lm",
                                 xtol=1e-12, ftol=1e-12, max_nfev=4000)
    theta = theta0 + res.x * scale
    # covariance from a central-difference Jacobian at the solution (steps of 1% of each
    # parameter's typical uncertainty), independent of the optimiser's internal one
    r0 = rfun(theta)
    J = np.empty((r0.size, theta.size))
    for j in range(theta.size):
        h = 0.01 * scale[j]
        tp, tm = theta.copy(), theta.copy()
        tp[j] += h
        tm[j] -= h
        J[:, j] = (rfun(tp) - rfun(tm)) / (2 * h)
    if n_lin is None:
        n_lin = sum(b.G.shape[1] + (len(shape_names(b.n_harm, b.sines, b.eclipses)) if shape_coefs is None else 1)
                    for b in bands)
    n = r0.size
    dof = n - n_lin - theta0.size
    c2 = float(r0 @ r0)
    cov = np.linalg.pinv(J.T @ J) * max(c2 / dof, 1.0)
    return theta, cov, c2, dof


def theta_to_ephemeris(theta, cov=None):
    """Translate (phi0, f0, f1) into P, P-dot and the eclipse time nearest T_REF."""
    phi0, f0, f1 = theta[:3]
    P = DAY / f0  # s
    pdot = -f1 / f0**2  # (d/d^2) / (1/d^2) -> dimensionless P-dot: dP/dt = -f1/f0^2 (days/day)
    n0 = np.round(phi0)  # eclipse number near T_REF
    # solve phi(t) = n0 near T_REF: x ~ (n0 - phi0)/f0
    x = (n0 - phi0) / f0
    for _ in range(5):
        x -= (phase(x, theta) - n0) / (f0 + f1 * x)
    out = {"P_s": P, "Pdot": pdot, "f0_per_day": f0, "f1_per_day2": f1,
           "T_ecl_bjd": T_REF + x, "phi0": phi0}
    if len(theta) > 3:
        out["f2_per_day3"] = theta[3]
    if cov is not None:
        g_P = np.array([0, -DAY / f0**2, 0] + [0] * (len(theta) - 3))
        g_Pd = np.array([0, 2 * f1 / f0**3, -1 / f0**2] + [0] * (len(theta) - 3))
        g_T = np.zeros(len(theta))
        g_T[0] = -1 / f0  # d(T_ecl)/d(phi0) (days)
        out["P_err_s"] = float(np.sqrt(g_P @ cov @ g_P))
        out["Pdot_err"] = float(np.sqrt(g_Pd @ cov @ g_Pd))
        out["T_ecl_err_s"] = float(np.sqrt(g_T @ cov @ g_T) * DAY)
        if len(theta) > 3:
            out["f2_err_per_day3"] = float(np.sqrt(cov[3, 3]))
    return out


def predicted_eclipse_time(theta, t_near):
    """Mid-eclipse BJD_TDB of the eclipse nearest ``t_near`` under ephemeris theta."""
    x = t_near - T_REF
    n = np.round(phase(x, theta))
    for _ in range(8):
        x -= (phase(x, theta) - n) / (theta[1] + theta[2] * x)
    return T_REF + x, int(n)


def phase_offset_fit(bands, theta, shape_coefs, grid_half_width=0.06, n_grid=241):
    """Best phase shift (cycles) of a fixed-shape template, with 1-sigma error.

    chi^2 is mapped on a grid of shifts and a parabola fitted around its minimum;
    the error comes from Delta chi^2 = 1 after rescaling chi^2 to its reduced value."""
    shifts = np.linspace(-grid_half_width, grid_half_width, n_grid)
    vals = np.array([chi2(np.r_[theta[0] + s, theta[1:]], bands, shape_coefs) for s in shifts])
    i = int(np.argmin(vals))
    lo, hi = max(i - 6, 0), min(i + 7, n_grid)
    a, b, c = np.polyfit(shifts[lo:hi], vals[lo:hi], 2)
    best = -b / (2 * a)
    n = sum(bb.t.size for bb in bands)
    n_par = sum(bb.G.shape[1] + 1 for bb in bands) + 1
    red = max((c - b**2 / (4 * a)) / max(n - n_par, 1), 1.0)
    err = np.sqrt(red / a) if a > 0 else np.nan
    return best, err, shifts, vals


def make_band(df: pd.DataFrame, name: str, group_col: str = "group", n_harm=4, sines=False,
              eclipses=True, flux_col="flux", err_col="err") -> BandData:
    names = sorted(df[group_col].unique())
    idx = {g: i for i, g in enumerate(names)}
    return BandData(name, df["bjd_tdb"].to_numpy(float), df["exptime_s"].to_numpy(float),
                    df[flux_col].to_numpy(float), df[err_col].to_numpy(float),
                    df[group_col].map(idx).to_numpy(int), names, n_harm, sines, eclipses)
