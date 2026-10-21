"""The formulas behind the HR diagram.

Everything the diagram shows is derived from four Gaia measurements per star:
apparent brightness G, the colour BP-RP, the parallax, and the proper motion.
This module holds every conversion between those and physical quantities, so
each one is written once and checked by ``test_physics.py``.

The central relation is the distance modulus. A star's apparent magnitude m
depends on how far away it is; its absolute magnitude M is the magnitude it
would have at a standard 10 parsecs. With distance d = 1/parallax,

    M = m - 5 log10(d / 10 pc) = m + 5 log10(parallax / mas) - 10

and plotting M rather than m is what turns a scatter of stars at all distances
into the HR diagram: it compares their intrinsic brightnesses.
"""

from __future__ import annotations

import numpy as np
from astropy import units as u

#: 1 au per Julian year in km/s. A proper motion of 1 mas/yr at a parallax of
#: 1 mas (a distance of 1 kpc) is a sideways speed of this many km/s.
KM_S_PER_AU_PER_YR = float((1 * u.au / u.yr).to(u.km / u.s).value)  # 4.74047

#: 5 / ln(10): converts a fractional distance error into a magnitude error.
MAG_PER_FRACTIONAL_DISTANCE = 5.0 / np.log(10.0)


# --------------------------------------------------------------------------
# Distance and brightness
# --------------------------------------------------------------------------

def distance_pc(parallax_mas):
    """Naive distance 1/parallax, in parsecs.

    Only trustworthy when the parallax is measured to much better than 10%;
    below that, 1/parallax is a biased estimator of distance (the error bar is
    symmetric in parallax, not in distance). The 100 pc sample enforces
    parallax_over_error > 10 for exactly this reason.
    """
    parallax_mas = np.asarray(parallax_mas, dtype=float)
    if np.any(parallax_mas <= 0):
        raise ValueError("a naive distance needs a positive parallax")
    return 1000.0 / parallax_mas


def absolute_magnitude(apparent_mag, parallax_mas):
    """M = m + 5 log10(parallax / mas) - 10.

    Equivalent to M = m - 5 log10(d / 10 pc) with d = 1000 / parallax pc. No
    extinction term: within 100 pc the Sun sits inside the Local Bubble, a
    cavity largely cleared of interstellar dust, so reddening is negligible for
    this sample. That would not be true much further out.
    """
    parallax_mas = np.asarray(parallax_mas, dtype=float)
    if np.any(parallax_mas <= 0):
        raise ValueError("absolute magnitude needs a positive parallax")
    return np.asarray(apparent_mag, dtype=float) + 5.0 * np.log10(parallax_mas) - 10.0


def apparent_magnitude(absolute_mag, distance_parsec):
    """Inverse of the distance modulus: m = M + 5 log10(d / pc) - 5."""
    return np.asarray(absolute_mag, dtype=float) + 5.0 * np.log10(distance_parsec) - 5.0


def absolute_magnitude_error(parallax_over_error):
    """1-sigma scatter in M from the parallax error alone, in magnitudes.

    Since M depends on 5 log10(parallax), a fractional parallax error f becomes
    a magnitude error of (5 / ln 10) * f = 2.17 f. At the sample's limit of
    parallax_over_error = 10 that is 0.22 mag; for most of the sample it is far
    smaller.
    """
    return MAG_PER_FRACTIONAL_DISTANCE / np.asarray(parallax_over_error, dtype=float)


def luminosity_relative_to_sun(absolute_mag, sun_absolute_mag):
    """How many times brighter than the Sun a star is, *in the same band*.

    This is a band luminosity, not a bolometric one: a cool red dwarf emits much
    of its light in the infrared, outside Gaia's G band, so its total luminosity
    relative to the Sun is higher than this ratio suggests.
    """
    return 10.0 ** (-0.4 * (np.asarray(absolute_mag, dtype=float) - sun_absolute_mag))


def combined_magnitude_offset(flux_ratio):
    """Magnitude change when a second star's light is added to the first.

    ``flux_ratio`` is secondary/primary flux. Two identical stars (ratio 1)
    appear 2.5 log10(2) = 0.753 mag brighter than one: that is why unresolved
    twin binaries form a second main sequence 0.75 mag above the first.
    """
    return -2.5 * np.log10(1.0 + np.asarray(flux_ratio, dtype=float))


def parallax_zero_point_shift(parallax_mas, zero_point_mas=-0.017):
    """Change in M if a parallax zero-point offset were corrected.

    Gaia DR3 parallaxes are on average slightly too small: quasars, which have
    no parallax, show a median of -0.017 mas (Lindegren et al. 2021, A&A 649,
    A4). Correcting it (parallax_true = parallax - zero_point) makes every star
    slightly closer and so slightly fainter. For this sample the shift is tiny,
    which is why it is quantified rather than applied.
    """
    parallax_mas = np.asarray(parallax_mas, dtype=float)
    return 5.0 * np.log10((parallax_mas - zero_point_mas) / parallax_mas)


# --------------------------------------------------------------------------
# Motion
# --------------------------------------------------------------------------

def tangential_velocity(pmra_mas_yr, pmdec_mas_yr, parallax_mas):
    """Speed across the sky, in km/s, relative to the Sun.

    v_tan = 4.74047 * mu / parallax, with mu the total proper motion in mas/yr
    and the parallax in mas. (``pmra`` in Gaia already includes the cos(dec)
    factor, so the two components add in quadrature directly.)

    This is two of the three components of the star's velocity; the third,
    along the line of sight, needs a spectrum. So v_tan is a lower limit on the
    star's true speed relative to the Sun: a halo star moving mostly towards or
    away from us will look slow.
    """
    mu = np.hypot(np.asarray(pmra_mas_yr, dtype=float), np.asarray(pmdec_mas_yr, dtype=float))
    return KM_S_PER_AU_PER_YR * mu / np.asarray(parallax_mas, dtype=float)


# --------------------------------------------------------------------------
# Photometric quality: the corrected BP/RP flux excess (Riello et al. 2021)
# --------------------------------------------------------------------------

def flux_excess_corrected(bp_rp, excess_factor):
    """C*: the BP/RP flux excess with its normal colour dependence removed.

    The BP and RP photometers measure light through a window about 3.5 x 2.1
    arcsec; G is measured more sharply. So a neighbour's light leaks into BP
    and RP more than into G, and (BP flux + RP flux) / G flux -- the "excess
    factor" C -- rises. Even clean stars have C > 1 because the passbands
    overlap, by an amount that depends on colour. Subtracting that expected
    value leaves C*, which is ~0 for clean isolated sources.

    Coefficients from Riello et al. (2021, A&A 649, A3), Table 2, as implemented
    in the Gaia collaboration's reference code.
    """
    x = np.asarray(bp_rp, dtype=float)
    c = np.asarray(excess_factor, dtype=float)
    expected = np.full_like(x, np.nan)
    blue = x < 0.5
    green = (x >= 0.5) & (x < 4.0)
    red = x >= 4.0
    expected[blue] = 1.154360 + 0.033772 * x[blue] + 0.032277 * x[blue] ** 2
    expected[green] = (1.162004 + 0.011464 * x[green] + 0.049255 * x[green] ** 2
                       - 0.005879 * x[green] ** 3)
    expected[red] = 1.057572 + 0.140537 * x[red]
    return c - expected


def flux_excess_sigma(g_mag):
    """1-sigma scatter of C* for clean stars at magnitude G (Riello+2021, eq. 18).

    Photometric noise grows towards faint magnitudes, so the tolerance does
    too: about 0.006 for bright stars, 0.014 at G = 15, 0.08 at G = 20.
    """
    g = np.asarray(g_mag, dtype=float)
    return 0.0059898 + 8.817481e-12 * g ** 7.618399
