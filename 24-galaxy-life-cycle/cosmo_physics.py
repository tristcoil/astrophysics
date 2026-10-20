"""Distances, luminosities and survey volumes for SDSS galaxies.

Adapted from project D2 (07-cosmic-web/cosmo_physics.py), with one change:
everything here uses the "concordance" cosmology (H0 = 70 km/s/Mpc,
Omega_m = 0.3, flat) because that is what the MPA-JHU and GSWLC catalogues
assumed when they turned fluxes into stellar masses and SFRs. Mixing
cosmologies would shift every mass and volume by a few per cent in ways that
are easy to miss.

Three pieces of physics:

1. **Distance from redshift**, through the expansion history of the
   adopted cosmology.
2. **K-corrections**: redshift slides a galaxy's spectrum through the camera's
   filters, so a distant galaxy is seen in bluer rest-frame light. We use the
   published polynomials of Chilingarian et al. (2010, 2012), as D2 did.
3. **V_max**, the volume within which each galaxy would have made it into a
   flux-limited survey (Schmidt 1968). A galaxy that could only have been seen
   nearby stands for many like it that the survey missed further out; weighting
   it by 1/V_max turns a flux-limited catalogue into densities per Mpc^3.
"""

from __future__ import annotations

import numpy as np
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM

from config import H0, OMEGA_M

C_KMS = 299_792.458
COSMO = FlatLambdaCDM(H0=H0, Om0=OMEGA_M)

# Tabulate once: over 0 < z < 0.4 interpolation is good to < 1e-7.
_Z_TABLE = np.linspace(0.0, 0.4, 8001)
_DC_TABLE = COSMO.comoving_distance(_Z_TABLE).to_value(u.Mpc)


def comoving_distance(z) -> np.ndarray:
    """Comoving distance in Mpc."""
    z = np.asarray(z, dtype=float)
    if np.any(z > _Z_TABLE[-1]):
        raise ValueError(f"redshift beyond interpolation table (z > {_Z_TABLE[-1]})")
    return np.interp(z, _Z_TABLE, _DC_TABLE)


def redshift_at_distance(dc_mpc) -> np.ndarray:
    return np.interp(np.asarray(dc_mpc, dtype=float), _DC_TABLE, _Z_TABLE)


def luminosity_distance(z) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    return (1.0 + z) * comoving_distance(z)


def angular_diameter_distance(z) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    return comoving_distance(z) / (1.0 + z)


def distance_modulus(z) -> np.ndarray:
    return 5.0 * np.log10(luminosity_distance(z) * 1e6 / 10.0)


def comoving_volume(z, solid_angle_sr: float) -> np.ndarray:
    """Comoving volume inside redshift z over a solid angle, Mpc^3."""
    return solid_angle_sr / 3.0 * comoving_distance(z) ** 3


# ---------------------------------------------------------------------------
# K-corrections: Chilingarian, Melchior & Zolotukhin (2010, MNRAS 405, 1409),
# coefficients of Chilingarian & Zolotukhin (2012, MNRAS 419, 1727), from the
# authors' calc_kcor.py (http://kcor.sai.msu.ru). Copied from D2, which
# checked them against the authors' code.
# ---------------------------------------------------------------------------

_KCOR_COEFFS = {
    "r_gr": [
        [0, 0, 0, 0],
        [1.83285, -2.71446, 4.97336, -3.66864],
        [-19.7595, 10.5033, 18.8196, 6.07785],
        [33.6059, -120.713, -49.299, 0],
        [144.371, 216.453, 0, 0],
        [-295.39, 0, 0, 0],
    ],
    "g_gr": [
        [0, 0, 0, 0],
        [-2.45204, 4.10188, 10.5258, -13.5889],
        [56.7969, -140.913, 144.572, 57.2155],
        [-466.949, 222.789, -917.46, -78.0591],
        [2906.77, 1500.8, 1689.97, 30.889],
        [-10453.7, -4419.56, -1011.01, 0],
        [17568, 3236.68, 0, 0],
        [-10820.7, 0, 0, 0],
    ],
}


def k_correction(band: str, z, g_minus_r) -> np.ndarray:
    """K-correction (mag) in SDSS ``band`` ('g' or 'r'), from observed g - r."""
    coeffs = _KCOR_COEFFS[f"{band}_gr"]
    z = np.asarray(z, dtype=float)
    c = np.asarray(g_minus_r, dtype=float)
    k = np.zeros(np.broadcast(z, c).shape)
    for x, row in enumerate(coeffs):
        for y, a in enumerate(row):
            if a:
                k = k + a * z**x * c**y
    return k


def absolute_magnitude(apparent, z, kcorr) -> np.ndarray:
    """Rest-frame absolute magnitude for H0 = 70."""
    return np.asarray(apparent, dtype=float) - distance_modulus(z) - np.asarray(kcorr)


# ---------------------------------------------------------------------------
# V_max
# ---------------------------------------------------------------------------


def z_at_magnitude(M, gr_obs, m_lim: float, z_lo: float, z_hi: float) -> np.ndarray:
    """Redshift at which a galaxy of absolute magnitude M appears at m_lim.

    Solves m_lim = M + DM(z) + K_r(z) by bisection on [z_lo, z_hi], holding the
    observed colour fixed in the K-correction (a small approximation over the
    range that matters). Returns z_lo or z_hi when the solution lies outside.
    """
    M = np.asarray(M, dtype=float)
    gr = np.asarray(gr_obs, dtype=float)
    lo = np.full_like(M, z_lo)
    hi = np.full_like(M, z_hi)

    def mag(z):
        return M + distance_modulus(z) + k_correction("r", z, gr)

    below = mag(lo) > m_lim   # fainter than the limit even at z_lo
    above = mag(hi) < m_lim   # still brighter than the limit at z_hi
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        bright = mag(mid) < m_lim
        lo = np.where(bright, mid, lo)
        hi = np.where(bright, hi, mid)
    out = 0.5 * (lo + hi)
    out = np.where(below, z_lo, out)
    return np.where(above, z_hi, out)


def vmax_limits(M, gr_obs, z_range, r_faint: float, r_bright: float):
    """(z_min, z_max) within which each galaxy would pass both flux limits.

    z_max: where it would fade to the faint limit. z_min: where it would be
    brighter than the bright limit. Both clipped to the sample's z_range.
    """
    zlo, zhi = z_range
    z_max = z_at_magnitude(M, gr_obs, r_faint, zlo, zhi)
    z_min = z_at_magnitude(M, gr_obs, r_bright, zlo, zhi)
    return z_min, z_max


def vmax(z_min, z_max, solid_angle_sr: float) -> np.ndarray:
    """Comoving volume between z_min and z_max over the survey area, Mpc^3."""
    return comoving_volume(z_max, solid_angle_sr) - comoving_volume(z_min, solid_angle_sr)


def cartesian(ra_deg, dec_deg, distance) -> np.ndarray:
    ra = np.radians(np.asarray(ra_deg, dtype=float))
    dec = np.radians(np.asarray(dec_deg, dtype=float))
    d = np.asarray(distance, dtype=float)
    return np.column_stack([d * np.cos(dec) * np.cos(ra), d * np.cos(dec) * np.sin(ra),
                            d * np.sin(dec)])


def unit_vectors(ra_deg, dec_deg) -> np.ndarray:
    return cartesian(ra_deg, dec_deg, 1.0)


def los_window_mpc(dv_kms: float, z) -> np.ndarray:
    """Comoving line-of-sight depth (Mpc) of a +/- dv velocity window at z."""
    z = np.asarray(z, dtype=float)
    h_z = COSMO.H(z).to_value(u.km / u.s / u.Mpc)
    return (1.0 + z) * dv_kms / h_z
