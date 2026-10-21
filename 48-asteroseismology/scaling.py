"""From numax and Dnu to mass, radius, surface gravity and luminosity.

Two scaling relations, both anchored on the Sun:

* numax scales with the acoustic cut-off frequency, i.e. with surface gravity
  over the square root of temperature (Brown et al. 1991; Kjeldsen & Bedding
  1995):   numax / numax_sun = (g / g_sun) (Teff / Teff_sun)^-1/2
* Dnu is the inverse of the sound-crossing time of the star, so it scales with
  the square root of the mean density (Ulrich 1986):
           Dnu / Dnu_sun = (rho / rho_sun)^1/2

Solving the pair for M and R:

    R / R_sun = (numax/numax_sun) (Dnu/Dnu_sun)^-2 (Teff/Teff_sun)^1/2
    M / M_sun = (numax/numax_sun)^3 (Dnu/Dnu_sun)^-4 (Teff/Teff_sun)^3/2

Correction factors are written APOKASC-3's way (Pinsonneault et al. 2025,
their Eqs. 1-2): numax -> f_numax numax and Dnu -> f_Dnu Dnu, with f = 1 for
the uncorrected relations. f_Dnu comes from stellar models (red giants are
not scaled-up Suns); f_numax is APOKASC-3's empirical calibration against
radii from Gaia parallaxes.
"""

from __future__ import annotations

import numpy as np

import config as C

LOGG_SUN = 4.438  # log10 of the nominal solar surface gravity, cgs (IAU 2015 B3: 27,400 cm/s^2)


def radius(numax, dnu, teff, *, numax_sun, dnu_sun, teff_sun=C.TEFF_SUN_NOMINAL,
           f_numax=1.0, f_dnu=1.0):
    x = f_numax * np.asarray(numax) / numax_sun
    d = f_dnu * np.asarray(dnu) / dnu_sun
    return x * d ** -2 * (np.asarray(teff) / teff_sun) ** 0.5


def mass(numax, dnu, teff, *, numax_sun, dnu_sun, teff_sun=C.TEFF_SUN_NOMINAL,
         f_numax=1.0, f_dnu=1.0):
    x = f_numax * np.asarray(numax) / numax_sun
    d = f_dnu * np.asarray(dnu) / dnu_sun
    return x ** 3 * d ** -4 * (np.asarray(teff) / teff_sun) ** 1.5


def logg(numax, teff, *, numax_sun, teff_sun=C.TEFF_SUN_NOMINAL, f_numax=1.0):
    return LOGG_SUN + np.log10(f_numax * np.asarray(numax) / numax_sun
                               * (np.asarray(teff) / teff_sun) ** 0.5)


def luminosity(r, teff, *, teff_sun=C.TEFF_SUN_NOMINAL):
    """L / L_sun = R^2 (Teff/Teff_sun)^4 (Stefan-Boltzmann)."""
    return np.asarray(r) ** 2 * (np.asarray(teff) / teff_sun) ** 4


def mean_density(dnu, *, dnu_sun, f_dnu=1.0):
    """rho / rho_sun = (Dnu / Dnu_sun)^2; the Sun's mean density is 1.41 g/cm^3."""
    return (f_dnu * np.asarray(dnu) / dnu_sun) ** 2


def monte_carlo(numax_samples, dnu_samples, teff, teff_err, *, numax_sun, dnu_sun,
                teff_sun=C.TEFF_SUN_NOMINAL, f_numax=1.0, f_dnu=1.0, f_dnu_err=0.0,
                seed: int = 0) -> dict:
    """Mass, radius, log g and L with uncertainties from posterior samples.

    numax and Dnu enter as samples (their posteriors), Teff and f_Dnu as
    Gaussians. Returns median and 16th/84th-percentile half-widths.
    """
    rng = np.random.default_rng(seed)
    n = min(len(numax_samples), len(dnu_samples))
    nm = rng.choice(np.asarray(numax_samples), n)
    dn = rng.choice(np.asarray(dnu_samples), n)
    te = teff + teff_err * rng.standard_normal(n)
    fd = f_dnu + f_dnu_err * rng.standard_normal(n)
    kw = dict(numax_sun=numax_sun, dnu_sun=dnu_sun, teff_sun=teff_sun, f_numax=f_numax)
    out = {}
    r = radius(nm, dn, te, f_dnu=fd, **kw)
    m = mass(nm, dn, te, f_dnu=fd, **kw)
    g = logg(nm, te, numax_sun=numax_sun, teff_sun=teff_sun, f_numax=f_numax)
    lum = luminosity(r, te, teff_sun=teff_sun)
    for name, v in (("radius", r), ("mass", m), ("logg", g), ("luminosity", lum)):
        q = np.percentile(v, [16, 50, 84])
        out[name] = float(q[1])
        out[name + "_err"] = float(0.5 * (q[2] - q[0]))
    return out


# ---------------------------------------------------------------------------
# APOKASC-3's numax calibration (Pinsonneault et al. 2025, Section 3.2)
# ---------------------------------------------------------------------------

#: Their Table 2, GARSTEC models with Mosser weighting (the default scale):
#: zero points for the red clump and the red giant branch above 50 uHz, and
#: the cubic in ln(numax) used on the red giant branch for 1.1 < numax < 50.
APOKASC_FNUMAX = {"RC": 0.9937, "RGB": 0.9959,
                  "poly": (-0.00556347, 0.04268084, -0.07452754, -0.02442568)}


def apokasc_f_numax(numax, state: str) -> float:
    """f_numax = (1 + p)^-1 for luminous RGB stars; constant zero points otherwise.

    Pinsonneault et al. write the RC/RGB zero points as "total scale factors
    relative to 3076 uHz" (an effective solar numax of 3063 and 3057 uHz,
    i.e. 3076 * 0.9959 and 3076 * 0.9937). Used only to illustrate the size
    of the calibration; per-star values come from their table.
    """
    if state == "RC":
        return APOKASC_FNUMAX["RC"]
    if numax >= 50.0:
        return APOKASC_FNUMAX["RGB"]
    a, b, c, d = APOKASC_FNUMAX["poly"]
    x = np.log(numax)
    return float(1.0 / (1.0 + (a * x ** 3 + b * x ** 2 + c * x + d)))
