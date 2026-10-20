"""From angles and redshifts to megaparsecs, and from flux to power.

A radio galaxy's size on the sky is an angle. Its physical (proper) length,
projected onto the sky, is that angle times the *angular diameter distance*:

    l_p = phi * D_A(z),    D_A = D_C / (1 + z),
    D_C = (c / H0) * integral_0^z dz' / E(z'),   E = sqrt(Om (1+z)^3 + 1 - Om)

in a flat universe. D_A has a famous quirk: it stops growing near z ~ 1.6 and
then shrinks, so very distant objects look *larger* than nearer ones of the
same size. Over the redshifts of the giants (z < ~2) a giant of fixed length
is smallest on the sky near that turnover -- 0.7 Mpc never subtends less than
~1.3' -- which is why every giant search can ignore smaller sources.

The radio luminosity (spectral power, W/Hz) from a flux density S, assuming a
power-law spectrum S ~ nu^alpha:

    L_nu = 4 pi D_L^2 S_nu (1 + z)^(-1 - alpha),   D_L = (1 + z) D_C.

Everything here is our own implementation; ``test_physics.py`` checks it
against astropy to better than 1e-5.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.integrate import cumulative_trapezoid

from config import Cosmology

C_KM_S = 299_792.458
MPC_M = 3.0856775814913673e22
ARCMIN_RAD = np.pi / (180 * 60)
#: Seconds per gigayear.
GYR_S = 3.15576e16

_Z_GRID = np.concatenate([np.linspace(0, 1e-3, 50, endpoint=False), np.geomspace(1e-3, 12, 6000)])


@lru_cache(maxsize=16)
def _tables(cosmo: Cosmology) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Comoving distance (Mpc) and lookback time (Gyr) on a fine redshift grid."""
    z = _Z_GRID
    om = cosmo.Om0
    e = np.sqrt(om * (1 + z) ** 3 + (1 - om))
    dh = C_KM_S / cosmo.H0  # Hubble distance, Mpc
    dc = dh * cumulative_trapezoid(1 / e, z, initial=0.0)
    th_gyr = (MPC_M / 1e3 / cosmo.H0) / GYR_S  # Hubble time, Gyr
    tl = th_gyr * cumulative_trapezoid(1 / ((1 + z) * e), z, initial=0.0)
    return z, dc, tl


def comoving_distance(z, cosmo: Cosmology) -> np.ndarray:
    zg, dc, _ = _tables(cosmo)
    return np.interp(np.asarray(z, float), zg, dc)


def angular_diameter_distance(z, cosmo: Cosmology) -> np.ndarray:
    z = np.asarray(z, float)
    return comoving_distance(z, cosmo) / (1 + z)


def luminosity_distance(z, cosmo: Cosmology) -> np.ndarray:
    z = np.asarray(z, float)
    return comoving_distance(z, cosmo) * (1 + z)


def lookback_time(z, cosmo: Cosmology) -> np.ndarray:
    zg, _, tl = _tables(cosmo)
    return np.interp(np.asarray(z, float), zg, tl)


def age_of_universe(cosmo: Cosmology) -> float:
    """Age today (Gyr): lookback time to z = infinity (closed form, flat LCDM)."""
    om, ol = cosmo.Om0, 1 - cosmo.Om0
    th_gyr = (MPC_M / 1e3 / cosmo.H0) / GYR_S
    return float(th_gyr * 2 / (3 * np.sqrt(ol)) * np.arcsinh(np.sqrt(ol / om)))


def age_at(z, cosmo: Cosmology) -> np.ndarray:
    return age_of_universe(cosmo) - lookback_time(z, cosmo)


def kpc_per_arcmin(z, cosmo: Cosmology) -> np.ndarray:
    return angular_diameter_distance(z, cosmo) * 1e3 * ARCMIN_RAD


def projected_length_mpc(las_arcmin, z, cosmo: Cosmology) -> np.ndarray:
    """Projected proper length (Mpc) of something ``las_arcmin`` long at redshift z."""
    las = np.asarray(las_arcmin, float)
    return las * ARCMIN_RAD * angular_diameter_distance(z, cosmo)


def angle_arcmin(length_mpc, z, cosmo: Cosmology) -> np.ndarray:
    """Angle (arcmin) subtended by a projected proper length at redshift z."""
    return np.asarray(length_mpc, float) / angular_diameter_distance(z, cosmo) / ARCMIN_RAD


def minimum_angle_arcmin(length_mpc: float, cosmo: Cosmology, z_max: float = 6.0) -> tuple[float, float]:
    """The smallest angle a given length can subtend, and the redshift where."""
    z = np.linspace(0.01, z_max, 6000)
    a = angle_arcmin(length_mpc, z, cosmo)
    i = int(np.argmin(a))
    return float(a[i]), float(z[i])


def radio_luminosity(flux_jy, z, cosmo: Cosmology, alpha: float = -0.7) -> np.ndarray:
    """Rest-frame spectral luminosity (W/Hz) at the observed frequency."""
    dl_m = luminosity_distance(z, cosmo) * MPC_M
    z = np.asarray(z, float)
    return 4 * np.pi * dl_m**2 * np.asarray(flux_jy, float) * 1e-26 * (1 + z) ** (-1 - alpha)


def comoving_volume_element(z, cosmo: Cosmology) -> np.ndarray:
    """dV_C/dz per steradian (Mpc^3): D_C^2 * D_H / E(z)."""
    z = np.asarray(z, float)
    e = np.sqrt(cosmo.Om0 * (1 + z) ** 3 + (1 - cosmo.Om0))
    return comoving_distance(z, cosmo) ** 2 * (C_KM_S / cosmo.H0) / e
