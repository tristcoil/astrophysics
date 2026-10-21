"""Gaia photometry: dust, extinction and a quality filter.

Dust between us and a star dims it and reddens it. At GD-1's high Galactic
latitudes the effect is small (E(B-V) ~ 0.01-0.1 mag) but not negligible next
to the ~0.05 mag width of an old population's colour-magnitude sequence, so
every colour and magnitude in the stream analysis is dereddened first.
"""

from __future__ import annotations

import functools

import numpy as np

from config import R_V, SFD_RECALIBRATION

#: Gaia DR2 extinction coefficients as a function of intrinsic colour and
#: extinction (Gaia Collaboration, Babusiaux et al. 2018, A&A 616, A10, Eq. 1
#: and Table 1): k_X = A_X / A_0 with A_0 = 3.1 E(B-V).
_K_COEFFS = {
    "G": (0.9761, -0.1704, 0.0086, 0.0011, -0.0438, 0.0013, 0.0099),
    "BP": (1.1517, -0.0871, -0.0333, 0.0173, -0.0230, 0.0006, 0.0043),
    "RP": (0.6104, -0.0170, -0.0026, -0.0017, -0.0078, 0.00005, 0.0006),
}


def extinction_coefficient(band: str, bp_rp0, a0):
    c1, c2, c3, c4, c5, c6, c7 = _K_COEFFS[band]
    x = np.asarray(bp_rp0, float)
    return c1 + c2 * x + c3 * x**2 + c4 * x**3 + c5 * a0 + c6 * a0**2 + c7 * x * a0


def deredden(g, bp_rp, ebv, iterations: int = 3):
    """Return (G_0, (BP-RP)_0). The coefficients depend on the *intrinsic*
    colour, so the correction is iterated from the observed colour."""
    a0 = R_V * np.asarray(ebv, float)
    colour0 = np.asarray(bp_rp, float).copy()
    for _ in range(iterations):
        k_bp = extinction_coefficient("BP", colour0, a0)
        k_rp = extinction_coefficient("RP", colour0, a0)
        colour0 = bp_rp - (k_bp - k_rp) * a0
    g0 = g - extinction_coefficient("G", colour0, a0) * a0
    return g0, colour0


@functools.lru_cache(maxsize=1)
def _sfd_map():
    from astropy.io import fits

    from fetch_data import SFD_PATH

    with fits.open(SFD_PATH) as h:
        return np.asarray(h[1].data["TEMPERATURE"], float), int(h[1].header["NSIDE"])


def ebv_sfd(ra, dec):
    """E(B-V) from Schlegel, Finkbeiner & Davis (1998), recalibrated by 0.86
    (Schlafly & Finkbeiner 2011). The map is the total dust column along the
    line of sight -- right for halo stars, which lie beyond all of it."""
    import healpy as hp

    from galaxy import galactic_lb

    m, nside = _sfd_map()
    l, b = galactic_lb(ra, dec)
    pix = hp.ang2pix(nside, l, b, lonlat=True, nest=True)
    return SFD_RECALIBRATION * m[pix]


# --------------------------------------------------------------------------
# Photometric quality (Riello et al. 2021, A&A 649, A3)
# --------------------------------------------------------------------------


def expected_excess(bp_rp):
    """BP + RP flux over G flux expected for clean photometry (Table 2)."""
    x = np.asarray(bp_rp, float)
    blue = 1.154360 + 0.033772 * x + 0.032277 * x**2
    mid = 1.162004 + 0.011464 * x + 0.049255 * x**2 - 0.005879 * x**3
    red = 1.057572 + 0.140537 * x
    return np.where(x < 0.5, blue, np.where(x < 4.0, mid, red))


def clean_colour(excess_factor, bp_rp, g, nsigma: float = 3.0):
    """True where BP/RP agree with G: |C*| < nsigma sigma_C*(G) (Eqs. 6, 18).

    A source whose BP and RP windows caught light from a neighbour (or
    background) that the narrower G window did not has a colour that cannot
    be trusted -- common in crowded fields such as M92's."""
    cstar = np.asarray(excess_factor, float) - expected_excess(bp_rp)
    sigma = 0.0059898 + 8.817481e-12 * np.asarray(g, float) ** 7.618399
    return np.isfinite(cstar) & (np.abs(cstar) < nsigma * sigma)


def distance_modulus(distance_kpc):
    return 5 * np.log10(np.asarray(distance_kpc, float) * 1000.0) - 5


def distance_from_modulus(dm):
    return 10 ** ((np.asarray(dm, float) + 5) / 5) / 1000.0
