"""Relative velocities of wide binaries from Gaia astrometry.

The observable is the scaled velocity (Pittordis & Sutherland 2019)

    v~ = v_sky / sqrt(G M / r_sky),

the sky-projected relative velocity in units of the Newtonian circular speed
at the projected separation. Following Banik et al. (2024, Section 2.3):

* both stars are put at one distance, the inverse-variance weighted mean of
  their parallax distances (their true distance difference, ~r_sky, is far
  smaller than the parallax errors at 100 pc, ~80 kAU);
* each star's 3D velocity is its tangential velocity plus the pair's systemic
  radial velocity along *its own* line of sight, which removes the
  "perspective" effect: a receding pair with no internal motion would
  otherwise seem to shrink on the sky (Shaya & Olling 2011; El-Badry 2019);
* the sky plane is perpendicular to the bisector of the two lines of sight.

Uncertainties come from a Monte Carlo over Gaia's parallax/proper-motion
covariance, the systemic RV error, a 0.024 dex mass scatter at fixed M_G, and
the unknown separation along the line of sight (their deprojection prior
P(x) ~ x^-2.6 / sqrt(x^2 - 1), x = r / r_sky).
"""

from __future__ import annotations

import numpy as np

from config import GM_SUN, AU_M, K_PM

#: (km/s)^2 for G * 1 Msun / 1 AU
GM_KMS2_AU = GM_SUN / AU_M / 1e6


def unit(ra_deg, dec_deg):
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    n = np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=-1)
    ea = np.stack([-np.sin(ra), np.cos(ra), np.zeros_like(ra)], axis=-1)
    ed = np.stack([-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)], axis=-1)
    return n, ea, ed


def weighted_distance(plx1, e1, plx2, e2):
    """Inverse-variance weighted mean of the two parallax *distances* (pc)."""
    d1, d2 = 1000.0 / plx1, 1000.0 / plx2
    s1, s2 = d1 * e1 / plx1, d2 * e2 / plx2
    w1, w2 = 1.0 / s1**2, 1.0 / s2**2
    return (w1 * d1 + w2 * d2) / (w1 + w2)


def sky_velocity(ra1, dec1, pmra1, pmdec1, ra2, dec2, pmra2, pmdec2, d1, d2, rv_sys,
                 perspective: bool = True):
    """Sky-plane relative velocity vector (km/s, 3D vector lying in the plane
    perpendicular to the bisector) and projected separation (AU).
    d1, d2 in pc; rv_sys in km/s (ignored if ``perspective`` is False)."""
    n1, ea1, ed1 = unit(ra1, dec1)
    n2, ea2, ed2 = unit(ra2, dec2)
    vt1 = K_PM * (d1 / 1000.0)[..., None] * (pmra1[..., None] * ea1 + pmdec1[..., None] * ed1)
    vt2 = K_PM * (d2 / 1000.0)[..., None] * (pmra2[..., None] * ea2 + pmdec2[..., None] * ed2)
    if perspective:
        vt1 = vt1 + rv_sys[..., None] * n1
        vt2 = vt2 + rv_sys[..., None] * n2
    dv = vt1 - vt2
    nm = n1 + n2
    nm = nm / np.linalg.norm(nm, axis=-1, keepdims=True)
    v_sky = dv - np.sum(dv * nm, axis=-1, keepdims=True) * nm
    # separation vector projected on the same plane (distances in pc -> AU)
    dr = (d1[..., None] * n1 - d2[..., None] * n2) * 206264.806
    r_sky = np.linalg.norm(dr - np.sum(dr * nm, axis=-1, keepdims=True) * nm, axis=-1)
    return v_sky, r_sky


def vtilde(v_sky_kms, r_sky_au, mtot):
    return v_sky_kms / np.sqrt(GM_KMS2_AU * mtot / r_sky_au)


def systemic_rv(rv1, e1, rv2, e2):
    """Inverse-variance weighted mean RV (one star's if only one is known)."""
    have1, have2 = np.isfinite(rv1) & np.isfinite(e1), np.isfinite(rv2) & np.isfinite(e2)
    w1 = np.where(have1, 1.0 / np.where(have1, e1, 1.0) ** 2, 0.0)
    w2 = np.where(have2, 1.0 / np.where(have2, e2, 1.0) ** 2, 0.0)
    wsum = w1 + w2
    rv = (w1 * np.nan_to_num(rv1) + w2 * np.nan_to_num(rv2)) / np.where(wsum > 0, wsum, 1.0)
    err = 1.0 / np.sqrt(np.where(wsum > 0, wsum, np.nan))
    rv = np.where(wsum > 0, rv, np.nan)
    return rv, err


# --------------------------------------------------------------------------
# Deprojection prior (Banik et al. 2024 eq. 6; Banik & Kroupa 2019 eq. 18)
# --------------------------------------------------------------------------
_X = np.concatenate([1.0 + np.logspace(-8, 0, 2000)[:-1], np.linspace(2.0, 200.0, 4000)])
_PDF = _X ** -2.6 / np.sqrt(_X**2 - 1.0)
_CDF = np.concatenate([[0.0], np.cumsum(0.5 * (_PDF[1:] + _PDF[:-1]) * np.diff(_X))])
_CDF /= _CDF[-1]


def sample_deprojection(rng, size):
    """Draws of x = r / r_sky from P(x) ~ x^-2.6 / sqrt(x^2 - 1)."""
    return np.interp(rng.random(size), _CDF, _X)


# --------------------------------------------------------------------------
# Monte Carlo error propagation
# --------------------------------------------------------------------------

def _cov3(s, which):
    """3x3 covariance of (parallax, pmra, pmdec) for star ``which`` (1/2)."""
    g = lambda k: np.asarray(s[f"{k}{which}"], float)
    ep, ea, ed = g("parallax_error"), g("pmra_error"), g("pmdec_error")
    rpa, rpd, rad = g("parallax_pmra_corr"), g("parallax_pmdec_corr"), g("pmra_pmdec_corr")
    C = np.empty((len(ep), 3, 3))
    C[:, 0, 0], C[:, 1, 1], C[:, 2, 2] = ep**2, ea**2, ed**2
    C[:, 0, 1] = C[:, 1, 0] = rpa * ep * ea
    C[:, 0, 2] = C[:, 2, 0] = rpd * ep * ed
    C[:, 1, 2] = C[:, 2, 1] = rad * ea * ed
    return C


def mc_vtilde(s: dict, mass_fn, *, n_trials: int = 1024, seed: int = 0, sigma_logm: float = 0.024,
              plx_err_scale: float = 1.0, chunk: int = 2000):
    """Monte Carlo spread of v~ for each pair in ``s`` (a dict of columns with
    suffixes 1/2: ra, dec, parallax, pmra, pmdec, errors and correlations,
    phot_g_mean_mag, phot_g_mean_flux_over_error, plus rv_sys, rv_sys_err and
    the nominal common distance ``dist`` in pc).
    ``mass_fn(M_G)`` gives stellar masses. Returns (median, std) of v~."""
    rng = np.random.default_rng(seed)
    n = len(s["ra1"])
    med = np.empty(n)
    std = np.empty(n)
    for lo in range(0, n, chunk):
        sl = slice(lo, min(lo + chunk, n))
        sub = {k: np.asarray(v)[sl] for k, v in s.items()}
        m = len(sub["ra1"])
        draws = []
        for w in (1, 2):
            C = _cov3(sub, w)
            if plx_err_scale != 1.0:
                S = np.diag([plx_err_scale, 1.0, 1.0])
                C = S @ C @ S
            L = np.linalg.cholesky(C + 1e-18 * np.eye(3)[None])
            z = rng.standard_normal((m, n_trials, 3))
            mean = np.stack([sub[f"parallax{w}"], sub[f"pmra{w}"], sub[f"pmdec{w}"]], axis=-1)
            draws.append(mean[:, None, :] + np.einsum("mij,mtj->mti", L, z))
        (p1, a1, d1), (p2, a2, d2) = [np.moveaxis(x, -1, 0) for x in draws]
        e1 = sub["parallax_error1"][:, None] * plx_err_scale
        e2 = sub["parallax_error2"][:, None] * plx_err_scale
        dist = weighted_distance(p1, e1, p2, e2)
        # line-of-sight separation (Banik et al. 2024 Section 2.3.2)
        theta = np.arccos(np.clip(np.sum(unit(sub["ra1"], sub["dec1"])[0] * unit(sub["ra2"], sub["dec2"])[0],
                                         axis=-1), -1, 1))
        rsky_au = theta * sub["dist"] * 206264.806
        sig_dd = np.hypot(1000 / sub["parallax1"] ** 2 * e1[:, 0], 1000 / sub["parallax2"] ** 2 * e2[:, 0])
        use_prior = (rsky_au / 206264.806) < sig_dd
        x = sample_deprojection(rng, (m, n_trials))
        dlos = rsky_au[:, None] * np.sqrt(x**2 - 1.0) / 206264.806 * np.where(rng.random((m, n_trials)) < 0.5, 1, -1)
        dd1 = np.where(use_prior[:, None], dist + 0.5 * dlos, 1000.0 / p1)
        dd2 = np.where(use_prior[:, None], dist - 0.5 * dlos, 1000.0 / p2)
        rv = sub["rv_sys"][:, None] + sub["rv_sys_err"][:, None] * rng.standard_normal((m, n_trials))
        rv = np.where(np.isfinite(rv), rv, 0.0)
        bc = lambda k: np.broadcast_to(sub[k][:, None], (m, n_trials))
        vsky, rsky = sky_velocity(bc("ra1"), bc("dec1"), a1, d1, bc("ra2"), bc("dec2"), a2, d2, dd1, dd2, rv)
        # masses: photometric error + distance + 0.024 dex at fixed M_G
        mags = []
        for w, dd in ((1, dd1), (2, dd2)):
            g = sub[f"phot_g_mean_mag{w}"][:, None]
            eg = 1.0857 / sub[f"phot_g_mean_flux_over_error{w}"][:, None]
            mg = g + eg * rng.standard_normal((m, n_trials)) - 5.0 * np.log10(dd / 10.0)
            mags.append(mass_fn(mg) * 10.0 ** (sigma_logm * rng.standard_normal((m, n_trials))))
        vt = vtilde(np.linalg.norm(vsky, axis=-1), rsky, mags[0] + mags[1])
        med[sl] = np.median(vt, axis=1)
        std[sl] = np.std(vt, axis=1)
    return med, std
