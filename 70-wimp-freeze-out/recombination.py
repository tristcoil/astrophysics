"""The third freeze-out: the electrons that never found a proton.

When the universe cooled to ~3000 K (z ~ 1100), protons and electrons combined into hydrogen. Had
recombination kept up with the cooling (the Saha equation), the free-electron fraction would have
fallen exponentially to ~1e-30 by z = 200. It did not: as the gas thinned, the recombination rate
n_e alpha_B dropped below the expansion rate H and the leftover fraction froze at ~3 x 10^-4,
just as WIMPs and neutrons froze. Those leftover electrons later matter for the CMB (they keep the
gas temperature coupled to the photons until z ~ 150).

This is a minimal model, enough to show the freeze-out: Peebles' (1968) three-level hydrogen atom
with the case-B recombination coefficient of Pequignot, Petitjean & Boisson (1991) multiplied by
RECFAST's fudge factor 1.14 (Seager, Sasselov & Scott 1999), the matter temperature equal to the
radiation temperature, helium already neutral (valid below z ~ 1600), and the expansion rate H(z)
taken from CAMB for the same cosmology. CAMB's full recombination (RECFAST with its corrections)
is the check. O4 will treat recombination properly.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import integrate

import config as C

K_B_EV = 8.617333262e-5           # eV/K
E_ION = 13.605693122994           # eV (Rydberg energy, CODATA)
E_2 = E_ION / 4.0                 # ionisation energy from n = 2
LAMBDA_2S1S = 8.2245809           # s^-1, two-photon decay 2s -> 1s
LAMBDA_ALPHA_CM = 121.5670e-7     # Lyman-alpha wavelength, cm
M_E_EV = C.M_E * 1e9
HBARC_EV_CM = C.HBARC_GEV_CM * 1e9
FUDGE = 1.14                      # RECFAST (Seager et al. 1999) for a three-level atom
Y_HE = 0.245                      # helium mass fraction (BBN; P4's own estimate is 0.2426-0.2469)
M_H_GEV = 0.93878307              # hydrogen atom mass in GeV


@dataclass
class Cosmology:
    H0: float = 67.36
    ombh2: float = C.OMEGA_B_H2
    omch2: float = C.OMEGA_C_H2
    T0_K: float = C.T0_K
    Y_he: float = Y_HE


def _camb_background(cos: Cosmology):
    import camb

    p = camb.set_params(H0=cos.H0, ombh2=cos.ombh2, omch2=cos.omch2, mnu=0.06, YHe=cos.Y_he,
                        TCMB=cos.T0_K)
    return camb.get_background(p)


def alpha_B(T_K):
    """Case-B recombination coefficient (cm^3/s), Pequignot et al. 1991, times the RECFAST fudge."""
    t = np.asarray(T_K) / 1e4
    return FUDGE * 1e-13 * 4.309 * t ** -0.6166 / (1.0 + 0.6703 * t ** 0.5300)


def saha_xe(T_K, n_H):
    """Hydrogen ionisation fraction in Saha equilibrium."""
    kT = K_B_EV * np.asarray(T_K)
    # (m_e kT / 2 pi)^(3/2) in cm^-3
    s = (M_E_EV * kT / (2 * math.pi)) ** 1.5 / HBARC_EV_CM ** 3 * np.exp(-E_ION / kT) / n_H
    return 0.5 * (np.sqrt(s * s + 4 * s) - s)


@dataclass
class RecResult:
    z: np.ndarray
    x_e: np.ndarray
    x_saha: np.ndarray
    gamma_over_H: np.ndarray       # C alpha_B n_e / H : the rate at which a free electron recombines
    t_s: np.ndarray                # cosmic time in seconds
    z_freeze: float                # where gamma/H = 1
    x_e_200: float
    x_e_camb: np.ndarray
    z_camb: np.ndarray


def run(cos: Cosmology | None = None, z_start: float = 1600.0, z_end: float = 20.0) -> RecResult:
    cos = cos or Cosmology()
    b = _camb_background(cos)
    h = cos.H0 / 100.0
    n_H0 = (1 - cos.Y_he) * cos.ombh2 * C.RHO_CRIT_H2_GEV_CM3 / M_H_GEV        # cm^-3
    zg = np.linspace(0.0, z_start * 1.2, 6001)
    Hg = b.hubble_parameter(zg) / (C.C_M_S / 1e3) * (C.C_M_S / 3.0856775814913673e22)  # 1/s
    lnH = np.log(Hg)

    def H_of(z):
        return math.exp(np.interp(z, zg, lnH))

    def rhs(z, y):
        x = y[0]
        T = cos.T0_K * (1 + z)
        nH = n_H0 * (1 + z) ** 3
        H = H_of(z)
        a = float(alpha_B(T))
        kT = K_B_EV * T
        beta = a * (M_E_EV * kT / (2 * math.pi)) ** 1.5 / HBARC_EV_CM ** 3 * math.exp(-E_2 / kT)
        K = LAMBDA_ALPHA_CM ** 3 / (8 * math.pi * H)
        Cp = (1 + K * LAMBDA_2S1S * nH * (1 - x)) / (1 + K * (LAMBDA_2S1S + beta) * nH * (1 - x))
        dxdt = -Cp * (a * nH * x * x - beta * (1 - x) * math.exp(-(E_ION - E_2) / kT))
        return [-dxdt / (H * (1 + z))]

    T0 = cos.T0_K * (1 + z_start)
    x0 = float(saha_xe(T0, n_H0 * (1 + z_start) ** 3))
    sol = integrate.solve_ivp(rhs, (z_start, z_end), [x0], method="Radau", rtol=1e-9, atol=1e-14,
                              dense_output=True)
    z = np.linspace(z_start, z_end, 2000)
    x = sol.sol(z)[0]
    T = cos.T0_K * (1 + z)
    nH = n_H0 * (1 + z) ** 3
    xs = saha_xe(T, nH)
    Hz = np.array([H_of(zz) for zz in z])
    kT = K_B_EV * T
    beta = alpha_B(T) * (M_E_EV * kT / (2 * math.pi)) ** 1.5 / HBARC_EV_CM ** 3 * np.exp(-E_2 / kT)
    K = LAMBDA_ALPHA_CM ** 3 / (8 * math.pi * Hz)
    Cp = (1 + K * LAMBDA_2S1S * nH * (1 - x)) / (1 + K * (LAMBDA_2S1S + beta) * nH * (1 - x))
    gam = Cp * alpha_B(T) * nH * x / Hz
    # Gamma/H first rises (the Peebles bottleneck relaxes as the gas cools), peaks near z ~ 900,
    # then falls through 1: the freeze-out is that downward crossing.
    k0 = int(np.argmax(gam))
    j = np.flatnonzero(gam[k0:] < 1.0)
    if len(j):
        i = k0 + j[0]
        r0, r1 = math.log(gam[i - 1]), math.log(gam[i])
        z_fr = float(z[i - 1] + (z[i] - z[i - 1]) * r0 / (r0 - r1))
    else:
        z_fr = float("nan")
    t_s = b.physical_time(z) * 1e9 * 365.25 * 86400.0
    zc = np.geomspace(z_end, 2500, 400)[::-1]
    xc = b.get_background_redshift_evolution(zc, ["x_e"], format="array")[:, 0]
    return RecResult(z=z, x_e=x, x_saha=xs, gamma_over_H=gam, t_s=t_s, z_freeze=z_fr,
                     x_e_200=float(np.interp(200.0, z[::-1], x[::-1])), x_e_camb=xc, z_camb=zc)
