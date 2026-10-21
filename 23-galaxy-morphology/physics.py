"""The physics: from a redshift and a magnitude to a luminosity and a colour.

The same formulas as project D2 (``07-cosmic-web/cosmo_physics.py``), so that
a galaxy's luminosity, colour and environment mean exactly the same thing in
both projects. Copied rather than imported, so each project folder stands on
its own; ``test_physics.py`` checks the two agree.

Three pieces of physics are needed:

1. **Distance from redshift.** Space expands, so light from a distant galaxy
   is stretched on the way to us; the stretch (the redshift z) grows with
   distance. Converting z to distance needs a cosmological model; we use
   Planck 2018's, as D2 did.

2. **Brightness to luminosity.** How bright a galaxy *looks* depends on how far
   away it is. The absolute magnitude M says how luminous it really is.

3. **K-corrections.** Redshift also slides a galaxy's spectrum through the
   camera's filters, so the r filter samples bluer rest-frame light for a
   distant galaxy than for a nearby one. K-corrections undo that, so colours
   and luminosities of galaxies at different distances compare like for like.
   This project needs four bands: u, g and r for the classic galaxy colours
   g - r and u - r, and i for the g - i colour that Hart et al. (2016) used.
"""

from __future__ import annotations

import numpy as np
from astropy import units as u
from astropy.cosmology import Planck18 as COSMO

#: Speed of light, km/s.
C_KMS = 299_792.458

#: Dimensionless Hubble constant of the adopted cosmology (Planck 2018).
H_LITTLE = float(COSMO.h)

# Astropy's Planck18 includes radiation and massive neutrinos, so every
# distance is a numerical integral. For hundreds of thousands of galaxies we
# tabulate once and interpolate; over 0 < z < 0.35 the error is < 1e-7.
_Z_TABLE = np.linspace(0.0, 0.35, 7001)
_DC_TABLE = COSMO.comoving_distance(_Z_TABLE).to_value(u.Mpc)


def comoving_distance(z) -> np.ndarray:
    """Comoving distance in Mpc (Planck 2018): the distance *today*."""
    z = np.asarray(z, dtype=float)
    if np.any(z > _Z_TABLE[-1]):
        raise ValueError(f"redshift beyond interpolation table (z > {_Z_TABLE[-1]})")
    return np.interp(z, _Z_TABLE, _DC_TABLE)


def comoving_distance_h(z) -> np.ndarray:
    """Comoving distance in h^-1 Mpc, the unit of D2's density measure."""
    return comoving_distance(z) * H_LITTLE


def luminosity_distance(z) -> np.ndarray:
    """D_L = (1 + z) D_C in a flat universe: the distance for the inverse-square law."""
    z = np.asarray(z, dtype=float)
    return (1.0 + z) * comoving_distance(z)


def distance_modulus(z) -> np.ndarray:
    """m - M for a source at redshift z: 5 log10(D_L / 10 pc)."""
    return 5.0 * np.log10(luminosity_distance(z) * 1e6 / 10.0)


def kpc_per_arcsec(z) -> np.ndarray:
    """Physical size of one arcsecond at redshift z (angular diameter distance)."""
    z = np.asarray(z, dtype=float)
    d_a_kpc = comoving_distance(z) / (1.0 + z) * 1e3
    return d_a_kpc * np.pi / (180.0 * 3600.0)


# ---------------------------------------------------------------------------
# K-corrections
# ---------------------------------------------------------------------------
#
# Polynomial approximations of Chilingarian, Melchior & Zolotukhin (2010,
# MNRAS 405, 1409), coefficients as updated in Chilingarian & Zolotukhin
# (2012, MNRAS 419, 1727), from the authors' calc_kcor.py (kcor.sai.msu.ru,
# simplified BSD licence). In each band the K-correction is a polynomial in
# redshift and one *observed* colour,
#
#     K(z, c) = sum_{x,y} a[x][y] * z**x * c**y,
#
# giving the rest-frame (z = 0) magnitude in the same filter.

_KCOR_COEFFS = {
    # r and g from observed g - r: identical to project D2.
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
    # u and r from observed u - r.
    "u_ur": [
        [0, 0, 0, 0],
        [10.3686, -6.12658, 2.58748, -0.299322],
        [-138.069, 45.0511, -10.8074, 0.95854],
        [540.494, -43.7644, 3.84259, 0],
        [-1005.28, 10.9763, 0, 0],
        [710.482, 0, 0, 0],
    ],
    "r_ur": [
        [0, 0, 0, 0],
        [3.03458, -1.50775, 0.576228, -0.0754155],
        [-47.8362, 19.0053, -3.15116, 0.286009],
        [154.986, -35.6633, 1.09562, 0],
        [-188.094, 28.1876, 0, 0],
        [68.9867, 0, 0, 0],
    ],
    # g and i from observed g - i.
    "g_gi": [
        [0, 0, 0, 0],
        [1.59269, -2.97991, 7.31089, -3.46913],
        [-27.5631, -9.89034, 15.4693, 6.53131],
        [161.969, -76.171, -56.1923, 0],
        [-204.457, 217.977, 0, 0],
        [-50.6269, 0, 0, 0],
    ],
    "i_gi": [
        [0, 0, 0, 0],
        [-2.21853, 3.94007, 0.678402, -1.24751],
        [-15.7929, -19.3587, 15.0137, 2.27779],
        [118.791, -40.0709, -30.6727, 0],
        [-134.571, 125.799, 0, 0],
        [-55.4483, 0, 0, 0],
    ],
}


def k_correction(band: str, z, colour, colour_name: str = "gr") -> np.ndarray:
    """K-correction (magnitudes) in SDSS ``band`` to the rest frame.

    ``colour`` is the *observed* colour named by ``colour_name`` ('gr', 'ur' or
    'gi'), corrected for Milky Way dust.
    """
    coeffs = _KCOR_COEFFS[f"{band}_{colour_name}"]
    z = np.asarray(z, dtype=float)
    c = np.asarray(colour, dtype=float)
    k = np.zeros(np.broadcast(z, c).shape)
    for x, row in enumerate(coeffs):
        for y, a in enumerate(row):
            if a:
                k = k + a * z**x * c**y
    return k


def rest_frame_colour(blue: str, red: str, colour_obs, z) -> np.ndarray:
    """Rest-frame colour (e.g. 'g', 'r'): the observed one minus both K-corrections."""
    name = blue + red
    return (np.asarray(colour_obs, dtype=float)
            - k_correction(blue, z, colour_obs, name)
            + k_correction(red, z, colour_obs, name))


def absolute_magnitude(apparent, z, kcorr) -> np.ndarray:
    """Rest-frame absolute magnitude for Planck 2018's H0: M = m - (m - M)(z) - K.

    Brighter is more negative: the Milky Way is about M_r = -21.
    """
    return np.asarray(apparent, dtype=float) - distance_modulus(z) - np.asarray(kcorr)


def to_literature_magnitude(m_abs, z, band_shift_offset: float, q_evolution: float):
    """D2's conversion of M_r to Zehavi et al. (2011)'s convention.

    Their magnitudes are K-corrected to the r band shifted to z = 0.1 (a mean
    offset, ``band_shift_offset``, measured from the data), evolved to their
    z = 0.1 brightness with M -> M + Q (z - 0.1), and quoted for h = 1.
    Needed only to pick D2's density tracers, which are defined in it.
    """
    z = np.asarray(z, dtype=float)
    return (np.asarray(m_abs, dtype=float) + band_shift_offset
            + q_evolution * (z - 0.1) - 5.0 * np.log10(H_LITTLE))


def h_magnitude_shift(h_other: float) -> float:
    """Add this to our M (Planck h) to quote it for another Hubble constant.

    Papers in this comparison mostly assume H0 = 70 km/s/Mpc. Distances scale
    as 1/H0, so at h = 0.70 a galaxy is 3% closer and 0.07 mag fainter.
    """
    return 5.0 * np.log10(h_other / H_LITTLE)


def bandshift_offset(kr, z) -> float:
    """D2's offset between the z = 0.1-shifted r band and ours.

    For a galaxy observed at exactly z = 0.1 the observed r filter *is* the
    shifted band, whose K-correction is -2.5 log10(1.1) for any spectrum, so the
    offset is our K_r(0.1) + 2.5 log10(1.1), averaged over galaxies near z = 0.1.
    """
    near = (np.asarray(z) > 0.095) & (np.asarray(z) < 0.105)
    return float(np.median(np.asarray(kr)[near]) + 2.5 * np.log10(1.1))


# ---------------------------------------------------------------------------
# Geometry for the density measure
# ---------------------------------------------------------------------------

def cartesian(ra_deg, dec_deg, distance) -> np.ndarray:
    """Equatorial (x, y, z): x toward RA 0, y toward RA 90, z toward the pole."""
    ra = np.radians(np.asarray(ra_deg, dtype=float))
    dec = np.radians(np.asarray(dec_deg, dtype=float))
    d = np.asarray(distance, dtype=float)
    return np.column_stack([d * np.cos(dec) * np.cos(ra),
                            d * np.cos(dec) * np.sin(ra),
                            d * np.sin(dec)])


def los_window(dv_kms: float, z) -> np.ndarray:
    """Comoving line-of-sight depth (h^-1 Mpc) of a +/-dv velocity window at z.

    A rest-frame velocity difference dv at redshift z is an observed redshift
    difference dz = (1 + z) dv / c, which spans c dz / H(z) of comoving distance.
    """
    z = np.asarray(z, dtype=float)
    h_z = COSMO.H(z).to_value(u.km / u.s / u.Mpc)
    return (1.0 + z) * dv_kms / h_z * H_LITTLE
