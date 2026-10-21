"""The history of cosmic star formation, from Madau & Dickinson's compilation.

Madau & Dickinson (2014, "MD14") collected measurements of the universe's
star-formation-rate density psi(z): how many solar masses of new stars form
per year in a typical (comoving) cubic megaparsec, at each epoch. Ultraviolet
surveys count young massive stars directly but lose light to dust; infrared
surveys catch that dust-absorbed light re-emitted. MD14 put both on one scale
(Salpeter IMF, a common faint-galaxy cut-off of 0.03 L*) in their Table 1, and
fitted

    psi(z) = A (1+z)^B / (1 + ((1+z)/C)^D)       [M_sun / yr / Mpc^3],

which rises as (1+z)^(B-D) at early times and falls as (1+z)^B at late times.

This module refits that function to their Table 1 ourselves, and derives from
it everything the review quotes: when the peak was, how fast star formation
has declined since, and -- by integrating over cosmic time -- how today's stars
built up. Each derived number is compared with MD14 in run.py.

Integrating psi over time gives the mass of stars ever formed. Not all of it
stays in stars: massive stars return part of their mass to the gas in winds
and supernovae. The fraction returned, R, is 0.27 for a Salpeter IMF (MD14,
section 2.3), so the stellar mass density today is

    rho_*(z) = (1 - R) * integral_{z}^{inf} psi(z') |dt/dz'| dz'.
"""

from __future__ import annotations

import numpy as np
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM, Planck18
from scipy.optimize import brentq, least_squares

import literature as lit
from config import H0, OMEGA_M, RETURN_FRACTION, SEED, Z_REIONIZATION_END

#: MD14's cosmology (also MPA-JHU's and GSWLC's): flat, H0 = 70, Omega_m = 0.3.
COSMO = FlatLambdaCDM(H0=H0, Om0=OMEGA_M)
#: Planck 2018 (project D3's cosmology), for checking how much the choice matters.
COSMO_PLANCK = Planck18

#: Redshift grid for time integrals: fine at low z, reaching z = 30 where the
#: integrand is negligible (psi(30) / psi_peak ~ 1e-4, and little time passes).
_Z_GRID = np.concatenate([np.linspace(0.0, 10.0, 4001)[:-1], np.geomspace(10.0, 30.0, 400)])


def psi_md(z, A: float, B: float, C: float, D: float) -> np.ndarray:
    """MD14's fitting function (their eq. 15), M_sun / yr / Mpc^3."""
    zp1 = 1.0 + np.asarray(z, dtype=float)
    return A * zp1**B / (1.0 + (zp1 / C) ** D)


def z_peak(B: float, C: float, D: float) -> float:
    """Redshift of the maximum of psi_md, analytically.

    d ln psi / d ln(1+z) = B - D x^D / (1 + x^D) with x = (1+z)/C; setting it to
    zero gives x^D = B / (D - B).
    """
    return C * (B / (D - B)) ** (1.0 / D) - 1.0


# ---------------------------------------------------------------------------
# Fitting MD14's Table 1
# ---------------------------------------------------------------------------


def table1_points(sfrd) -> dict:
    """Redshift (bin midpoint), log psi and asymmetric errors from Table 1."""
    z = 0.5 * (sfrd["z_lo"].to_numpy() + sfrd["z_hi"].to_numpy())
    return {
        "z": z,
        "z_lo": sfrd["z_lo"].to_numpy(),
        "z_hi": sfrd["z_hi"].to_numpy(),
        "logpsi": sfrd["log_value"].to_numpy(),
        "err_plus": sfrd["err_plus"].to_numpy(),
        "err_minus": sfrd["err_minus"].to_numpy(),
        "band": sfrd["band"].to_numpy(),
        "reference": sfrd["reference"].to_numpy(),
    }


def _residuals(theta, z, y, ep, em):
    """Asymmetric-error residuals in log psi.

    A point whose quoted upper error is large is allowed to sit far below the
    model, and vice versa: the error used is the one on the side the model
    lies.
    """
    A, B, C, D = theta
    model = np.log10(psi_md(z, 10**A, B, C, D))
    err = np.where(model > y, ep, em)
    return (y - model) / err


def fit_md(points: dict, start=None) -> dict:
    """Least-squares fit of psi_md to log psi, with asymmetric errors.

    MD14 do not say how they fitted; this is the plain chi-square choice. The
    normalisation is fitted as log10(A) so all four parameters are of order one.
    """
    p0 = start or [np.log10(lit.MD14_FIT["A"]), lit.MD14_FIT["B"], lit.MD14_FIT["C"],
                   lit.MD14_FIT["D"]]
    res = least_squares(_residuals, p0,
                        args=(points["z"], points["logpsi"], points["err_plus"],
                              points["err_minus"]),
                        bounds=([-4, 0.5, 1.2, 2.0], [0, 6, 8, 12]))
    A, B, C, D = res.x
    chi2 = float(np.sum(res.fun**2))
    dof = len(points["z"]) - 4
    return {"A": float(10**A), "B": float(B), "C": float(C), "D": float(D),
            "chi2": chi2, "dof": dof, "z_peak": float(z_peak(B, C, D))}


def bootstrap_fit(points: dict, n: int = 1000, seed: int = SEED) -> np.ndarray:
    """Refit to resampled point sets: the spread is our parameter uncertainty.

    The four parameters are strongly correlated (the peak's height, position
    and the two slopes trade off), and MD14 note that each survey's error bars
    were computed differently. Resampling whole points, rather than trusting
    the formal least-squares covariance, lets the data say how well the curve
    and the quantities derived from it are pinned down.
    """
    rng = np.random.default_rng(seed)
    m = len(points["z"])
    out = []
    for _ in range(n):
        idx = rng.integers(0, m, m)
        sub = {k: v[idx] for k, v in points.items()}
        try:
            f = fit_md(sub)
        except ValueError:
            continue
        out.append([f["A"], f["B"], f["C"], f["D"]])
    return np.array(out)


# ---------------------------------------------------------------------------
# Time integrals: the build-up of stellar mass
# ---------------------------------------------------------------------------


def cosmic_time_gyr(z, cosmo=COSMO) -> np.ndarray:
    """Age of the universe at redshift z, Gyr."""
    return cosmo.age(np.asarray(z, dtype=float)).to_value(u.Gyr)


def history(params: dict, cosmo=COSMO, R: float = RETURN_FRACTION) -> dict:
    """psi, the mass formed and the stellar mass density along the z grid.

    Integrated in cosmic time with the trapezoid rule on a grid fine enough
    that halving the spacing changes rho_*(0) by < 1e-5 (checked in tests).
    Returns arrays ordered by increasing redshift (decreasing time).
    """
    z = _Z_GRID
    t = cosmic_time_gyr(z, cosmo)  # decreasing with z
    psi = psi_md(z, params["A"], params["B"], params["C"], params["D"])
    # mass formed after z: integral from t(z) to t0 of psi dt (yr)
    seg = 0.5 * (psi[1:] + psi[:-1]) * (t[:-1] - t[1:]) * 1e9
    formed_after = np.concatenate([[0.0], np.cumsum(seg)])
    total = formed_after[-1]  # stars formed before z = 30 are negligible
    formed_before = total - formed_after
    return {"z": z, "t_gyr": t, "psi": psi, "formed_before": formed_before,
            "rho_star": (1.0 - R) * formed_before, "total_formed": total,
            "t0_gyr": float(t[0])}


def z_when_fraction_formed(h: dict, fraction: float) -> float:
    """Redshift by which ``fraction`` of today's stars had formed.

    With R constant, a fraction of stellar mass equals the same fraction of
    mass formed (MD14 use the same instantaneous-recycling approximation).
    """
    f = h["formed_before"] / h["total_formed"]  # decreasing with z
    return float(np.interp(fraction, f[::-1], h["z"][::-1]))


def fraction_formed_before(h: dict, z: float) -> float:
    f = h["formed_before"] / h["total_formed"]
    return float(np.interp(z, h["z"], f))


def mean_stellar_age_gyr(h: dict) -> float:
    """MD14's mass-weighted mean age: t0 - int t psi dt / int psi dt."""
    t, psi = h["t_gyr"][::-1], h["psi"][::-1]  # increasing time
    return float(h["t0_gyr"] - np.trapezoid(t * psi, t) / np.trapezoid(psi, t))


def efold_time_gyr(z, logpsi, sigma=None, z_max: float = 1.0, cosmo=COSMO) -> dict:
    """Fit psi ∝ exp(-t / tau) in cosmic time to points with z < z_max.

    MD14: "the SFRD data at z < 1 can also be fit quite well by an exponential
    decline with cosmic time and an e-folding timescale of 3.9 Gyr".
    """
    z, logpsi = np.asarray(z), np.asarray(logpsi)
    sel = z < z_max
    t = cosmic_time_gyr(z[sel], cosmo)
    y = logpsi[sel] * np.log(10.0)  # natural log
    if sigma is None:
        coef, cov = np.polyfit(t, y, 1, cov=True)
    else:  # polyfit's weights multiply the residuals: w = 1 / sigma
        w = 1.0 / (np.asarray(sigma)[sel] * np.log(10.0))
        coef, cov = np.polyfit(t, y, 1, w=w, cov="unscaled")
    slope, err = coef[0], np.sqrt(cov[0, 0])
    return {"tau_gyr": float(-1.0 / slope), "tau_err_gyr": float(err / slope**2),
            "n_points": int(sel.sum())}


def efold_of_curve(params: dict, z_max: float = 1.0, cosmo=COSMO) -> float:
    """The same e-folding time, measured on the fitted curve itself."""
    zz = np.linspace(0.0, z_max, 200)
    t = cosmic_time_gyr(zz, cosmo)
    y = np.log(psi_md(zz, **{k: params[k] for k in "ABCD"}))
    return float(-1.0 / np.polyfit(t, y, 1)[0])


def high_z_slope(z, logpsi, z_min: float = 3.0) -> dict:
    """Power-law index of psi vs (1+z) for the points at z > z_min."""
    z, logpsi = np.asarray(z), np.asarray(logpsi)
    sel = z > z_min
    coef, cov = np.polyfit(np.log10(1.0 + z[sel]), logpsi[sel], 1, cov=True)
    return {"slope": float(coef[0]), "slope_err": float(np.sqrt(cov[0, 0])),
            "n_points": int(sel.sum())}


def baryon_density(omega_b_h2: float) -> float:
    """Mean comoving baryon density, M_sun / Mpc^3, from Omega_b h^2."""
    rho_crit_h2 = 2.775e11  # 3 (100 km/s/Mpc)^2 / (8 pi G), M_sun/Mpc^3 per h^2
    return omega_b_h2 * rho_crit_h2


def summarise(params: dict, cosmo=COSMO) -> dict:
    """Every number MD14 quote, recomputed from a set of fit parameters."""
    h = history(params, cosmo)
    zp = z_peak(params["B"], params["C"], params["D"])
    return {
        "z_peak": float(zp),
        "t_peak_gyr": float(cosmic_time_gyr(zp, cosmo)),
        "psi_peak": float(psi_md(zp, **{k: params[k] for k in "ABCD"})),
        "psi_z0": float(psi_md(0.0, **{k: params[k] for k in "ABCD"})),
        "psi_z7_over_z0": float(psi_md(7.0, **{k: params[k] for k in "ABCD"})
                                / psi_md(0.0, **{k: params[k] for k in "ABCD"})),
        "rho_star_z0": float(h["rho_star"][0]),
        "z_half_mass": z_when_fraction_formed(h, 0.5),
        "frac_before_peak": fraction_formed_before(h, zp),
        "z_last_quarter": z_when_fraction_formed(h, 0.75),  # 25% formed after it
        "frac_before_reionization_end": fraction_formed_before(h, Z_REIONIZATION_END),
        "mean_stellar_age_gyr": mean_stellar_age_gyr(h),
        "efold_curve_gyr": efold_of_curve(params, cosmo=cosmo),
        "t0_gyr": h["t0_gyr"],
    }


def find_redshift_at_time(t_gyr: float, cosmo=COSMO) -> float:
    """Redshift at which the universe was t_gyr old."""
    if t_gyr >= cosmic_time_gyr(0.0, cosmo):  # today (or a rounding step past it)
        return 0.0
    return brentq(lambda z: cosmic_time_gyr(z, cosmo) - t_gyr, 0.0, 50.0)
