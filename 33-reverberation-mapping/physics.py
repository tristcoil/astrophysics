"""Small physics helpers: light-travel radii, virial products, luminosities.

The reverberation mass is M = f c tau sigma^2 / G, where tau is the rest-frame
lag (so R = c tau is the BLR radius in light-days), sigma the line-of-sight
velocity dispersion of the gas, and f a dimensionless factor for the BLR's
unknown geometry and orientation. The part we measure, c tau sigma^2 / G, is
called the virial product (VP).
"""

from __future__ import annotations

import numpy as np

import config as C


def virial_product(tau_days, sigma_kms):
    """c tau sigma^2 / G in solar masses (tau rest-frame days, sigma km/s)."""
    tau_days = np.asarray(tau_days, float)
    sigma_kms = np.asarray(sigma_kms, float)
    return C.C_CGS * tau_days * C.DAY_S * (sigma_kms * 1e5) ** 2 / C.G_CGS / C.MSUN_G


def virial_product_err(tau, tau_err, sigma, sigma_err):
    """Error on the virial product from independent lag and width errors."""
    vp = virial_product(tau, sigma)
    return vp * np.hypot(np.asarray(tau_err) / tau, 2 * np.asarray(sigma_err) / sigma)


def galactic_extinction_5100(ebv: float) -> float:
    """A_lambda at 5100 A in magnitudes, Cardelli et al. (1989) law with R_V = 3.1.

    CCM89 at 1/lambda = 1.96 um^-1 gives A_5100/A_V = 1.094 (between the B
    band's 1.32 and the V band's 1.00), i.e. A_5100 = 3.39 E(B-V).
    """
    x = 1.0 / 0.51  # inverse microns
    y = x - 1.82
    a = 1 + 0.17699 * y - 0.50447 * y**2 - 0.02427 * y**3 + 0.72085 * y**4 + 0.01979 * y**5 \
        - 0.77530 * y**6 + 0.32999 * y**7
    b = 1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3 - 5.38434 * y**4 - 0.62251 * y**5 \
        + 5.30260 * y**6 - 2.09002 * y**7
    rv = 3.1
    return float(ebv * rv * (a + b / rv))


def log_lambda_l5100(flux_1e15, z: float = C.Z_NGC5548, dl_mpc: float = C.DL_NGC5548_MPC,
                     ebv: float = 0.0):
    """log10 of lambda L_lambda at rest-frame 5100 A from an observed flux density.

    ``flux_1e15`` is F_lambda at observed wavelength 5100 (1 + z) in
    1e-15 erg/s/cm^2/A (the convention of Bentz et al. 2013 and Peterson et al.
    2013): lambda L_lambda = 4 pi D_L^2 [5100 (1+z)] F_lambda. ``ebv`` > 0 also
    corrects for Milky Way dust.
    """
    f = np.asarray(flux_1e15, float) * 1e-15
    if ebv:
        f = f * 10 ** (0.4 * galactic_extinction_5100(ebv))
    d = dl_mpc * C.MPC_CM
    return np.log10(4 * np.pi * d**2 * 5100.0 * (1 + z) * f)


def eddington_luminosity(mbh_msun):
    return C.L_EDD_PER_MSUN * np.asarray(mbh_msun, float)


def schwarzschild_radius_cm(mbh_msun):
    return 2 * C.G_CGS * np.asarray(mbh_msun, float) * C.MSUN_G / C.C_CGS**2
