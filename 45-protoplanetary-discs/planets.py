"""From a gap's width to the mass of the planet that might have carved it (DSHARP VII).

A planet in a gas disc pushes gas away from its orbit (it gives angular momentum to the
gas outside and takes it from the gas inside), so a gap opens. Viscosity -- turbulent
mixing, parametrised by alpha -- fills it back in, and the gas pressure resists. Zhang et
al. (2018) ran 45 hydrodynamic simulations (planet/star mass ratio q, disc thickness h/r,
alpha) with dust grains of many sizes, made synthetic 1.27 mm images, and fitted the
dust gap width Delta = (r_out - r_in)/r_out as

    Delta = A K'^B,    K' = q (h/r)^-0.18 alpha^-0.31,

with A and B depending on how well the largest grains couple to the gas (their Stokes
number St_max). Inverting gives q = (Delta/A)^(1/B) (h/r)^0.18 alpha^0.31.

The disc thickness comes from its temperature: a passively heated, flared disc has
T(r) = (phi L* / (8 pi r^2 sigma_SB))^(1/4) with phi = 0.02, and h/r = c_s / v_K with
c_s = sqrt(k T / (mu m_H)), mu = 2.35.
"""

from __future__ import annotations

import numpy as np

from config import (AU_CGS, DUST_OPACITY_CM2_G, FLARING_PHI, G_CGS, K_B_CGS, L_SUN_CGS, M_H_CGS, M_JUP_CGS,
                    M_SUN_CGS, MEAN_MOLECULAR_WEIGHT, SIGMA_SB_CGS, ZHANG_ST_UNIT, ZHANG_WIDTH_FITS,
                    ZHANG_WIDTH_FITS_DSD2)

M_SUN_IN_MJUP = M_SUN_CGS / M_JUP_CGS


def midplane_temperature(r_au, lstar_lsun: float, phi: float = FLARING_PHI):
    r = np.asarray(r_au, dtype=float) * AU_CGS
    return (phi * lstar_lsun * L_SUN_CGS / (8.0 * np.pi * r**2 * SIGMA_SB_CGS)) ** 0.25


def sound_speed(t_k):
    return np.sqrt(K_B_CGS * np.asarray(t_k, dtype=float) / (MEAN_MOLECULAR_WEIGHT * M_H_CGS))


def aspect_ratio(r_au, lstar_lsun: float, mstar_msun: float, phi: float = FLARING_PHI):
    """h/r = c_s / v_K at radius r."""
    r = np.asarray(r_au, dtype=float) * AU_CGS
    v_k = np.sqrt(G_CGS * mstar_msun * M_SUN_CGS / r)
    return sound_speed(midplane_temperature(r_au, lstar_lsun, phi)) / v_k


def toomre_sigma(r_au, lstar_lsun: float, mstar_msun: float, q_toomre: float = 1.0):
    """Gas surface density at which the disc would become gravitationally unstable,
    Sigma_GI = c_s Omega / (pi G Q) (DSHARP VII uses Q = 1)."""
    r = np.asarray(r_au, dtype=float) * AU_CGS
    omega = np.sqrt(G_CGS * mstar_msun * M_SUN_CGS / r**3)
    return sound_speed(midplane_temperature(r_au, lstar_lsun)) * omega / (np.pi * G_CGS * q_toomre)


def _nearest_st(st: float, table: dict) -> float:
    keys = [k for k in table if isinstance(k, float) and table[k] is not None]
    return min(keys, key=lambda k: abs(np.log10(k) - np.log10(st)))


def width_fit(st_units: float, model: str) -> tuple[float, float, float, float, float]:
    """(St_max, A, B, -err, +err in log10 K') for a Stokes number given in DSHARP VII's units
    of 5.23e-4. The "1mm" models use the DSD1 fit with the same St_max (their Sec. 4)."""
    st = st_units * ZHANG_ST_UNIT
    table = ZHANG_WIDTH_FITS_DSD2 if model == "DSD2" else ZHANG_WIDTH_FITS
    key = _nearest_st(st, table)
    return (key,) + tuple(table[key])


def mass_ratio(delta, h_over_r: float, alpha: float, a: float, b: float):
    """Invert Delta = A K'^B for q."""
    k_prime = (np.asarray(delta, dtype=float) / a) ** (1.0 / b)
    return k_prime * h_over_r**0.18 * alpha**0.31


def planet_mass_mjup(delta, h_over_r: float, alpha: float, st_units: float, model: str, mstar_msun: float):
    """Planet mass in Jupiter masses with the fit's own 68% range (in log10)."""
    st, a, b, lo, hi = width_fit(st_units, model)
    q = mass_ratio(delta, h_over_r, alpha, a, b)
    m = q * mstar_msun * M_SUN_IN_MJUP
    return m, m * 10.0 ** (-lo), m * 10.0 ** hi


def sigma_dust_thin(i_cgs, t_k, nu_hz: float, kappa: float = DUST_OPACITY_CM2_G):
    """Dust surface density for optically thin emission, I = B_nu(T) kappa Sigma."""
    from alma_images import planck

    return np.asarray(i_cgs) / (planck(nu_hz, t_k) * kappa)


def optical_depth(i_cgs, t_k, nu_hz: float):
    """tau from I = B_nu(T) (1 - exp(-tau)) (DSHARP II eq. 2); NaN where I >= B."""
    from alma_images import planck

    x = np.asarray(i_cgs) / planck(nu_hz, t_k)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where((x > 0) & (x < 1), -np.log1p(-x), np.nan)
