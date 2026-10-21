"""The physics: from a redshift and a magnitude to a position and a luminosity.

A spectroscopic survey measures, for each galaxy, where it is on the sky, how
much its light is stretched (the redshift z), and how bright it looks. Turning
that into a 3D map and a census needs four pieces of physics, all here:

1. **Distance from redshift.** Space expands, so the light of a distant
   galaxy is stretched on its way to us; the stretch grows with distance
   (Hubble's law). Converting z to distance needs a cosmological model -- how
   fast the universe expanded at each epoch -- and we use Planck 2018's.

2. **Brightness to luminosity.** How bright a galaxy looks depends on how far
   away it is; its *absolute* magnitude says how luminous it really is.

3. **K-corrections.** Redshift also slides a galaxy's spectrum through the
   camera's filters, so the r filter samples bluer rest-frame light for a
   distant galaxy than for a nearby one. The K-correction undoes that, so
   galaxies at different distances are compared like for like.

4. **Positions in 3D.** Sky coordinates plus distance give x, y, z.

Distances come in two unit systems. For maps we quote Mpc (and light-years)
using Planck 2018's Hubble constant, H0 = 67.7 km/s/Mpc. For clustering we
quote h^-1 Mpc, where h = H0 / 100 km/s/Mpc: the literature's convention, which
makes the numbers independent of the (contested -- see project D1) value of H0.
"""

from __future__ import annotations

import numpy as np
from astropy import units as u
from astropy.cosmology import Planck18 as COSMO

#: Speed of light, km/s.
C_KMS = 299_792.458

#: Dimensionless Hubble constant of the adopted cosmology.
H_LITTLE = float(COSMO.h)

#: Megaparsecs to light-years.
MPC_IN_LY = float((1 * u.Mpc).to(u.lyr).value)

# Astropy's Planck18 includes radiation and massive neutrinos, so every
# distance is a numerical integral. For 650,000 galaxies we tabulate once and
# interpolate; over 0 < z < 0.35 this is accurate to better than 1 part in 10^7.
_Z_TABLE = np.linspace(0.0, 0.35, 7001)
_DC_TABLE = COSMO.comoving_distance(_Z_TABLE).to_value(u.Mpc)


def comoving_distance(z) -> np.ndarray:
    """Comoving distance in Mpc (Planck 2018 cosmology).

    The distance between us and the galaxy *today*, with the expansion that
    happened while its light travelled scaled out. It is the natural distance for
    a map: two galaxies drifting apart only through cosmic expansion keep a fixed
    comoving separation, so structures in the map are not distorted by when we
    happen to see them.
    """
    z = np.asarray(z, dtype=float)
    if np.any(z > _Z_TABLE[-1]):
        raise ValueError(f"redshift beyond interpolation table (z > {_Z_TABLE[-1]})")
    return np.interp(z, _Z_TABLE, _DC_TABLE)


def comoving_distance_h(z) -> np.ndarray:
    """Comoving distance in h^-1 Mpc, the clustering literature's unit."""
    return comoving_distance(z) * H_LITTLE


def redshift_at_distance(dc_mpc) -> np.ndarray:
    """Inverse of :func:`comoving_distance`: the redshift at a distance in Mpc."""
    return np.interp(np.asarray(dc_mpc, dtype=float), _DC_TABLE, _Z_TABLE)


def luminosity_distance(z) -> np.ndarray:
    """Luminosity distance in Mpc: D_L = (1 + z) D_C for a flat universe.

    The distance that makes the inverse-square law work for an expanding
    universe. The (1 + z) factor is the energy each photon loses to redshift
    together with the stretched time between photons' arrivals.
    """
    z = np.asarray(z, dtype=float)
    return (1.0 + z) * comoving_distance(z)


def distance_modulus(z) -> np.ndarray:
    """m - M for a source at redshift z: 5 log10(D_L / 10 pc)."""
    return 5.0 * np.log10(luminosity_distance(z) * 1e6 / 10.0)


def lookback_time_gyr(z) -> np.ndarray:
    """How long ago the light we see now left the galaxy, in billions of years."""
    return COSMO.lookback_time(np.asarray(z, dtype=float)).to_value(u.Gyr)


# ---------------------------------------------------------------------------
# K-corrections
# ---------------------------------------------------------------------------
#
# Polynomial approximations of Chilingarian, Melchior & Zolotukhin (2010,
# MNRAS 405, 1409), coefficients updated in Chilingarian & Zolotukhin (2012,
# MNRAS 419, 1727), from the authors' calc_kcor.py (http://kcor.sai.msu.ru,
# simplified BSD licence). The K-correction in a band is a polynomial in the
# redshift and one *observed* colour:
#
#     K(z, c) = sum_{x,y} a[x][y] * z**x * c**y
#
# giving the rest-frame (z = 0) magnitude in the same filter. They are fits to
# the full kcorrect template machinery, accurate to a few hundredths of a
# magnitude for z < 0.5 -- far better than we need.

_KCOR_COEFFS = {
    # r band from observed g - r
    "r_gr": [
        [0, 0, 0, 0],
        [1.83285, -2.71446, 4.97336, -3.66864],
        [-19.7595, 10.5033, 18.8196, 6.07785],
        [33.6059, -120.713, -49.299, 0],
        [144.371, 216.453, 0, 0],
        [-295.39, 0, 0, 0],
    ],
    # g band from observed g - r
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
    """K-correction (magnitudes) in SDSS ``band`` ('g' or 'r') to rest frame.

    ``g_minus_r`` is the *observed* colour, corrected for Milky Way dust.
    """
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
    """Rest-frame absolute magnitude: M = m - (m - M)(z) - K.

    Quoted for Planck 2018's H0. Brighter is more negative: the Milky Way is
    about M_r = -21, a typical bright galaxy in this survey -21 to -22.
    """
    return np.asarray(apparent, dtype=float) - distance_modulus(z) - np.asarray(kcorr)


def rest_frame_color(g_minus_r_obs, z) -> np.ndarray:
    """Rest-frame g - r: the observed colour with both K-corrections removed."""
    return (np.asarray(g_minus_r_obs, dtype=float)
            - k_correction("g", z, g_minus_r_obs)
            + k_correction("r", z, g_minus_r_obs))


def to_literature_magnitude(m_abs, z, band_shift_offset: float, q_evolution: float):
    """Convert our rest-frame M_r into the convention of Zehavi et al. (2011).

    Three differences separate the two:

    * **Bandpass.** They K-correct to the r filter blueshifted to z = 0.1
      (written ^{0.1}r), we to the r filter itself. For a given galaxy the
      difference is a colour term, the same at every redshift;
      ``band_shift_offset`` is its average over our galaxies (computed in
      ``analysis.py``, ~0.2 mag).
    * **Evolution.** They evolve every galaxy to its z = 0.1 brightness with
      M -> M + Q (z - 0.1): stars age and fade, so a galaxy seen nearby was
      brighter at z = 0.1.
    * **Hubble constant.** They quote magnitudes for h = 1 ("M - 5 log h");
      with distances ~1/h, that shifts every magnitude by -5 log10(h).
    """
    z = np.asarray(z, dtype=float)
    return (np.asarray(m_abs, dtype=float) + band_shift_offset
            + q_evolution * (z - 0.1) - 5.0 * np.log10(H_LITTLE))


def apparent_magnitude_limit_curve(z, r_limit: float, kcorr: float = 0.0):
    """The faintest absolute magnitude visible at redshift z for a flux limit.

    Everything fainter than this curve at a given z drops out of the survey --
    the edge of the empty wedge in the magnitude-redshift diagram.
    """
    return r_limit - distance_modulus(z) - kcorr


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def cartesian(ra_deg, dec_deg, distance) -> np.ndarray:
    """Equatorial (x, y, z) with x toward RA 0, y toward RA 90, z toward the pole."""
    ra = np.radians(np.asarray(ra_deg, dtype=float))
    dec = np.radians(np.asarray(dec_deg, dtype=float))
    d = np.asarray(distance, dtype=float)
    return np.column_stack([d * np.cos(dec) * np.cos(ra),
                            d * np.cos(dec) * np.sin(ra),
                            d * np.sin(dec)])


def rp_pi(p1, p2) -> tuple[np.ndarray, np.ndarray]:
    """Split pair separations into across (r_p) and along (pi) the line of sight.

    The line of sight is the direction to the pair's midpoint (Fisher et al.
    1994, MNRAS 267, 927): with s = p2 - p1 and l = (p1 + p2) / 2,
    pi = |s . l| / |l| and r_p^2 = |s|^2 - pi^2.

    Why split at all: redshift measures velocity, not distance. A galaxy's own
    orbital motion adds to its cosmological redshift, displacing it *along the
    line of sight only*. So pi is contaminated by those motions while r_p --
    measured on the sky -- is not.
    """
    p1 = np.atleast_2d(p1)
    p2 = np.atleast_2d(p2)
    s = p2 - p1
    l = p1 + p2  # the factor 1/2 cancels in the projection
    s2 = np.einsum("ij,ij->i", s, s)
    pi = np.abs(np.einsum("ij,ij->i", s, l)) / np.sqrt(np.einsum("ij,ij->i", l, l))
    rp = np.sqrt(np.clip(s2 - pi**2, 0.0, None))
    return rp, pi


def shell_volume(solid_angle_sr: float, d_min, d_max):
    """Comoving volume of a shell between two distances over a solid angle."""
    return solid_angle_sr / 3.0 * (np.asarray(d_max) ** 3 - np.asarray(d_min) ** 3)


def los_window(dv_kms: float, z) -> np.ndarray:
    """Comoving line-of-sight depth (h^-1 Mpc) of a +/-dv velocity window at z.

    A rest-frame velocity difference dv at redshift z is an observed redshift
    difference dz = (1 + z) dv / c, which spans c dz / H(z) of comoving distance.
    """
    z = np.asarray(z, dtype=float)
    h_z = COSMO.H(z).to_value(u.km / u.s / u.Mpc)
    return (1.0 + z) * dv_kms / h_z * H_LITTLE


def velocity_dispersion(z) -> tuple[float, float]:
    """Line-of-sight velocity dispersion of a cluster from member redshifts.

    Returns (sigma, standard error) in km/s, in the cluster's rest frame:
    sigma_v = c * std(z) / (1 + mean z). The error is the large-N approximation
    sigma / sqrt(2 (N - 1)).
    """
    z = np.asarray(z, dtype=float)
    n = len(z)
    if n < 3:
        return float("nan"), float("nan")
    zbar = z.mean()
    sigma = C_KMS * z.std(ddof=1) / (1.0 + zbar)
    return float(sigma), float(sigma / np.sqrt(2.0 * (n - 1)))


def wp_power_law(rp, r0, gamma):
    """Projected correlation function of a power law xi(r) = (r / r0)^-gamma.

    w_p(r_p) = r_p (r_p / r0)^-gamma Gamma(1/2) Gamma((gamma - 1)/2) / Gamma(gamma/2)
    (Davis & Peebles 1983), the form Zehavi et al. (2011) fit.
    """
    from scipy.special import gamma as G

    rp = np.asarray(rp, dtype=float)
    return rp * (rp / r0) ** (-gamma) * G(0.5) * G((gamma - 1.0) / 2.0) / G(gamma / 2.0)
