"""Mass from motion: gas surface densities, the gas's own pull, and density profiles.

**How much gas, where.** The 21-cm line is optically thin, so the integrated
brightness counts hydrogen atoms: N_HI = 1.823e18 cm^-2 per K km/s. Averaging
the THINGS moment-0 map around each tilted ring (including the empty pixels
in it) and correcting to face-on gives Sigma_HI(R); helium and metals add 40%.

**The gas's own pull.** A flat disc does not pull like a point mass. Its
potential is the sum of thin rings, each with the potential
-(2 G M / pi) K(m) / sqrt((R + a)^2 + z^2), m = 4aR / ((R + a)^2 + z^2),
with K the complete elliptic integral. We sum rings on a fine grid, average
over a sech^2 vertical layer, and differentiate: V^2 = R dPhi/dR. Checked
against Freeman's (1970) exact result for an exponential disc in
``test_physics.py``.

**Density from rotation.** For a spherical mass distribution, the rotation
curve alone gives the density (de Blok et al. 2001, eq. 1):

    rho(R) = (1 / 4 pi G) [ 2 (V / R) dV/dR + (V / R)^2 ].

A cold-dark-matter cusp has rho ~ R^-1 near the centre; a core has rho ~ R^0.
Oh et al. (2011) measured the inner slope alpha = dlog rho / dlog R with a
least-squares line through the points inside the "break radius", where the
slope changes fastest; we follow their recipe exactly so their DDO 154 value
(from de Blok et al.'s curve) can be reproduced before ours is quoted.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares
from scipy.special import ellipk

from config import G_KPC, GAS_PER_HI, NHI_PER_K_KMS

#: Hydrogen atoms per cm^2 in a column of 1 Msun/pc^2.
HI_ATOMS_PER_MSUN_PC2 = 1.2489e20


def brightness_temperature_k(s_jy_beam, bmaj_arcsec, bmin_arcsec, freq_ghz=1.420405752):
    """Rayleigh-Jeans brightness temperature of a flux density per beam."""
    return 1.222e6 * s_jy_beam / (freq_ghz ** 2 * bmaj_arcsec * bmin_arcsec)


def hi_surface_density(mom0_jy_beam_kms, beam) -> np.ndarray:
    """Observed (projected) HI surface density in Msun/pc^2, without helium."""
    nhi = NHI_PER_K_KMS * brightness_temperature_k(mom0_jy_beam_kms, beam[0], beam[1])
    return nhi / HI_ATOMS_PER_MSUN_PC2


def hi_mass_msun(total_flux_jy_kms: float, distance_mpc: float) -> float:
    """M_HI = 2.356e5 D^2 S dv for optically thin gas."""
    return 2.356e5 * distance_mpc ** 2 * total_flux_jy_kms


def ring_profile(values: np.ndarray, x: np.ndarray, y: np.ndarray, tr, *, fill=0.0) -> np.ndarray:
    """Average of a map over each tilted ring, empty pixels included.

    ``values``, ``x``, ``y`` cover the whole map (not just detected pixels):
    an empty pixel in a ring is a measurement of "no gas", not missing data.
    """
    from tiltedring import _select

    vals = np.where(np.isfinite(values), values, fill)
    out = np.full(len(tr.radius), np.nan)
    for k in range(len(tr.radius)):
        sel = _select(x, y, vals, tr.ring_params(k), tr.radius[k], tr.width, "both")
        if sel.any():
            out[k] = vals[sel].mean()
    return out


def disc_rotation_curve(r_eval_kpc, r_kpc, sigma_msun_pc2, *, z0_kpc=0.0, n_grid=3000,
                        n_z=24) -> np.ndarray:
    """Signed V^2 (km/s)^2 of an axisymmetric disc with surface density Sigma(R).

    Sigma is interpolated linearly between the given radii and taken as zero
    beyond the last one. With ``z0_kpc`` > 0 the layer is sech^2(z / z0).
    Negative values mean an outward pull (e.g. a central hole in the gas).
    """
    r_max = r_kpc.max()
    edges = np.linspace(0.0, r_max, n_grid + 1)
    a = 0.5 * (edges[1:] + edges[:-1])
    da = edges[1] - edges[0]
    sig = np.interp(a, r_kpc, sigma_msun_pc2, left=sigma_msun_pc2[0], right=0.0) * 1e6  # Msun/kpc^2
    mass = 2 * np.pi * a * sig * da  # Msun per ring

    if z0_kpc > 0:
        zq, wq = np.polynomial.legendre.leggauss(n_z)
        z = 0.5 * (zq + 1) * 6 * z0_kpc  # 0 .. 6 z0
        w = wq * 3 * z0_kpc / np.cosh(z / z0_kpc) ** 2
        w = w / w.sum()
    else:
        z, w = np.array([0.0]), np.array([1.0])

    def phi(rr):
        out = np.zeros_like(rr)
        for zz, ww in zip(z, w):
            denom2 = (rr[:, None] + a[None, :]) ** 2 + zz ** 2
            m = np.clip(4 * a[None, :] * rr[:, None] / denom2, 0, 1 - 1e-15)
            out += ww * np.sum(-(2 * G_KPC / np.pi) * mass[None, :] * ellipk(m) / np.sqrt(denom2), axis=1)
        return out

    r_eval = np.atleast_1d(np.asarray(r_eval_kpc, dtype=float))
    # The potential is evaluated on the cell edges, always half a cell from
    # the nearest ring (whose kernel has a logarithmic singularity), then
    # differentiated there.
    grid = edges[1:]
    if r_eval.max() > r_max:
        grid = np.concatenate([grid, np.arange(r_max + da, 1.05 * r_eval.max() + 2 * da, da)])
    dphi = np.gradient(phi(grid), grid)
    return r_eval * np.interp(r_eval, grid, dphi)


def freeman_disc_v2(r_kpc, sigma0_msun_pc2, h_kpc) -> np.ndarray:
    """Exact V^2 of a thin exponential disc (Freeman 1970)."""
    from scipy.special import i0, i1, k0, k1

    y = np.asarray(r_kpc) / (2 * h_kpc)
    return 4 * np.pi * G_KPC * sigma0_msun_pc2 * 1e6 * h_kpc * y ** 2 * (i0(y) * k0(y) - i1(y) * k1(y))


# --------------------------------------------------------------------------
# Density profiles and the inner slope
# --------------------------------------------------------------------------

def density_from_curve(r_kpc, v_kms) -> np.ndarray:
    """rho(R) in Msun/pc^3 for a spherical distribution (de Blok et al. 2001)."""
    r, v = np.asarray(r_kpc, float), np.asarray(v_kms, float)
    dvdr = np.gradient(v, r)
    rho = (2 * (v / r) * dvdr + (v / r) ** 2) / (4 * np.pi * G_KPC)  # Msun/kpc^3
    return rho / 1e9


def inner_slope(r_kpc, rho, *, max_break_fraction=0.5) -> dict:
    """Oh et al. (2011) recipe for the logarithmic inner slope alpha.

    The break radius is where the log-log slope changes most rapidly (largest
    |second difference|) within the inner ``max_break_fraction`` of the
    curve; alpha is a least-squares line through the points inside it. The
    uncertainty is the mean change in alpha when the next point outside is
    added and when the point at the break is dropped.
    """
    ok = np.isfinite(rho) & (rho > 0) & (np.asarray(r_kpc) > 0)
    lr, lrho = np.log10(np.asarray(r_kpc)[ok]), np.log10(np.asarray(rho)[ok])
    n = len(lr)
    slope = np.diff(lrho) / np.diff(lr)
    curv = np.abs(np.diff(slope))  # at points 1 .. n-2
    lim = max(3, int(max_break_fraction * n))
    j = int(np.argmax(curv[1:lim - 1])) + 2  # index of the break point, >= 2 so >= 3 points inside

    def fit(k):
        return float(np.polyfit(lr[:k + 1], lrho[:k + 1], 1)[0])

    alpha = fit(j)
    alt = [fit(j + 1), fit(j - 1)] if j - 1 >= 1 else [fit(j + 1)]
    return {"alpha": alpha, "alpha_err": float(np.mean([abs(a - alpha) for a in alt])),
            "r_break_kpc": float(10 ** lr[j]), "n_points": j + 1,
            "r_innermost_kpc": float(10 ** lr[0])}


# --------------------------------------------------------------------------
# Halo models
# --------------------------------------------------------------------------

def v_iso(r_kpc, rho0_msun_pc3, rc_kpc):
    """Pseudo-isothermal sphere: a cored halo."""
    r = np.asarray(r_kpc, float)
    return np.sqrt(4 * np.pi * G_KPC * rho0_msun_pc3 * 1e9 * rc_kpc ** 2 * (1 - (rc_kpc / r) * np.arctan(r / rc_kpc)))


def v_nfw(r_kpc, c, v200, h0=70.0):
    """NFW halo with concentration c and circular speed V200 at R200."""
    r = np.asarray(r_kpc, float)
    r200 = v200 / (10 * h0 / 1e3)  # kpc: V200 = 10 H0 R200 (H0 in km/s/kpc)
    x = r / r200
    num = np.log(1 + c * x) - c * x / (1 + c * x)
    den = x * (np.log(1 + c) - c / (1 + c))
    return v200 * np.sqrt(num / den)


def fit_halo(kind, r, v, err, v_baryon2=None):
    """Least-squares halo fit to V_obs^2 - V_baryons^2 (minimum disc if none)."""
    vb2 = np.zeros_like(v) if v_baryon2 is None else v_baryon2

    def model(p):
        vh = v_iso(r, *p) if kind == "ISO" else v_nfw(r, *p)
        return np.sqrt(np.clip(vh ** 2 + vb2, 0, None))

    p0, lo, hi = ((0.03, 1.0), (1e-4, 0.01), (10, 50)) if kind == "ISO" else ((8, 60), (0.1, 5), (100, 500))
    sol = least_squares(lambda p: (v - model(p)) / err, p0, bounds=(lo, hi))
    dof = max(len(v) - 2, 1)
    try:
        cov = np.linalg.inv(sol.jac.T @ sol.jac)
        perr = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        perr = np.full(2, np.nan)
    return {"params": sol.x.tolist(), "errors": perr.tolist(), "chi2_red": float(np.sum(sol.fun ** 2) / dof),
            "rms_kms": float(np.sqrt(np.mean((v - model(sol.x)) ** 2)))}


#: Hubble parameter h of the Planck cosmology Dutton & Maccio (2014) simulated
#: (as in project D4).
DUTTON_MACCIO_H = 0.671


def nfw_concentration_lcdm(m200_msun):
    """Median NFW concentration of cold-dark-matter haloes at z = 0.

    Dutton & Maccio (2014), eq. 8: log10 c200 = 0.905 - 0.101 log10(M200 /
    [1e12 Msun/h]), scatter 0.11 dex. The same relation D4 used.
    """
    m_h = np.asarray(m200_msun, dtype=float) * DUTTON_MACCIO_H
    return 10.0 ** (0.905 - 0.101 * np.log10(m_h / 1e12))


def m200_from_v200(v200, h0=70.0):
    """M200 = V200^3 / (10 G H0), H0 in km/s/Mpc."""
    return v200 ** 3 / (10 * G_KPC * h0 / 1e3)


def fit_lcdm_nfw(r, v, err, r_min_kpc):
    """NFW halo with the LCDM concentration for its mass, matched to the outer curve.

    The one free number, V200, is fitted to the points beyond ``r_min_kpc``;
    the concentration follows from the halo's mass. The answer is the cusp a
    cold-dark-matter universe would put in this galaxy if it holds the
    observed outer rotation.
    """
    outer = r >= r_min_kpc

    def resid(p):
        c = nfw_concentration_lcdm(m200_from_v200(p[0]))
        return (v[outer] - v_nfw(r[outer], c, p[0])) / err[outer]

    sol = least_squares(resid, [60.0], bounds=([5.0], [500.0]))
    v200 = float(sol.x[0])
    return {"v200": v200, "c": float(nfw_concentration_lcdm(m200_from_v200(v200))),
            "m200": float(m200_from_v200(v200)),
            "chi2_red_outer": float(np.sum(sol.fun ** 2) / max(outer.sum() - 1, 1))}


def smooth_derivative(r, v, window=5):
    """dV/dR from local quadratic fits over ``window`` points, anchored at V(0) = 0.

    A plain finite difference of a rotation curve is too noisy for the
    density formula (Oh et al. 2011 must have smoothed too: their plotted
    densities imply derivatives no finite difference of the curve gives).
    """
    rr = np.concatenate([[0.0], np.asarray(r, float)])
    vv = np.concatenate([[0.0], np.asarray(v, float)])
    out = np.empty(len(r))
    h = window // 2
    for k in range(len(r)):
        j = k + 1
        lo, hi = max(0, j - h), min(len(rr), j + h + 1)
        deg = 2 if hi - lo >= 4 else 1
        out[k] = np.polyfit(rr[lo:hi] - rr[j], vv[lo:hi], deg)[-2]
    return out


def density_smoothed(r_kpc, v_kms, window=5) -> np.ndarray:
    """rho(R) in Msun/pc^3 using ``smooth_derivative``."""
    r, v = np.asarray(r_kpc, float), np.asarray(v_kms, float)
    rho = (2 * (v / r) * smooth_derivative(r, v, window) + (v / r) ** 2) / (4 * np.pi * G_KPC)
    return rho / 1e9


def slope_within(r_kpc, rho, r_break_kpc) -> dict:
    """Oh et al.'s least-squares inner slope for a given break radius, with
    their uncertainty (mean change when the next point is added and when the
    break point is dropped)."""
    r = np.asarray(r_kpc, float)
    ok = np.isfinite(rho) & (rho > 0)
    lr, lrho = np.log10(r[ok]), np.log10(np.asarray(rho)[ok])
    n_in = int(np.sum(r[ok] <= r_break_kpc * 1.0001))

    def fit(n):
        return float(np.polyfit(lr[:n], lrho[:n], 1)[0])

    alpha = fit(n_in)
    alt = [fit(n_in + 1)] + ([fit(n_in - 1)] if n_in - 1 >= 2 else [])
    return {"alpha": alpha, "alpha_err": float(np.mean([abs(a - alpha) for a in alt])),
            "n_points": n_in, "r_break_kpc": float(10 ** lr[n_in - 1])}


def velocity_power_slope(r_kpc, v_kms, r_max_kpc, r_min_kpc=0.0) -> dict:
    """Inner density slope from the rotation curve itself, no derivatives.

    If rho ~ R^alpha then M ~ R^(3 + alpha) and V ~ R^beta with
    beta = (2 + alpha) / 2. A line through log V vs log R over the inner
    points gives beta, hence alpha = 2 beta - 2: 0 for a core (solid-body
    rotation, V ~ R), -1 for an NFW cusp (V ~ R^0.5).
    """
    r, v = np.asarray(r_kpc, float), np.asarray(v_kms, float)
    sel = (r <= r_max_kpc) & (r >= r_min_kpc) & (v > 0) & np.isfinite(v)
    coef, cov = np.polyfit(np.log10(r[sel]), np.log10(v[sel]), 1, cov=True) if sel.sum() > 3 else \
        (np.polyfit(np.log10(r[sel]), np.log10(v[sel]), 1), np.full((2, 2), np.nan))
    beta = float(coef[0])
    return {"beta": beta, "beta_err": float(np.sqrt(cov[0, 0])), "alpha": 2 * beta - 2,
            "alpha_err": 2 * float(np.sqrt(cov[0, 0])), "n_points": int(sel.sum())}
