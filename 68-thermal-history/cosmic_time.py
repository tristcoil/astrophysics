"""The cosmic calendar: time, horizon and horizon mass at every temperature, 10 TeV to today.

The Friedmann equation H^2 = 8 pi rho / (3 M_Pl^2) with rho(T) = rho_EoS(T) (everything in the
thermal bath, from ``eos.build``) + rho_matter + rho_Lambda (+ the massive neutrino once it slows
down). Entropy conservation fixes the scale factor, a ∝ 1/(g_s^(1/3) T), so

    dt = -(dT/T) [1 + (1/3) d ln g_s/d ln T] / H,

and the age at temperature T is the integral of that from infinity down to T. Above 10 TeV, g*
is taken constant (the radiation-era t = 1/(2H)); what came before (inflation, reheating) is
unknown and shifts nothing measurable here: the clock is the "hot Big Bang" clock.

Below 10 keV the universe is photons + neutrinos + matter + dark energy: Planck 2018's late
universe (TT,TE,EE+lowE+lensing+BAO; omega_b, omega_c, H0, one 0.06 eV neutrino), flat.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline

import config as C

PI2 = np.pi ** 2
KMSMPC_TO_GEV = 1e3 / C.MPC_M * C.HBAR_GEV_S          # H0 in km/s/Mpc -> GeV
RHO_CRIT_H2 = 3 * (100 * KMSMPC_TO_GEV) ** 2 * C.M_PL ** 2 / (8 * np.pi)  # GeV^4 per h^2


@lru_cache(maxsize=None)
def _frozen_fd_table():
    """rho/(g T_nu^4) of a fermion frozen out while relativistic, f = 1/(e^(p/T_nu) + 1), vs y = m/T_nu."""
    ys = np.concatenate([[0.0], np.logspace(-4, 4, 500)])

    def rho(y):
        f = lambda q: q * q * np.sqrt(q * q + y * y) * np.exp(-q) / (1 + np.exp(-q))
        return quad(f, 0, np.inf, epsabs=0, epsrel=1e-11, limit=300)[0] / (2 * PI2)

    vals = np.array([rho(y) for y in ys])
    return CubicSpline(np.log(ys[1:]), np.log(vals[1:])), vals[0]


def rho_frozen_neutrino(T_nu, m, g=2.0):
    """Energy density (GeV^4) of one massive neutrino flavour (nu + anti-nu) with a frozen FD spectrum."""
    sp, r0 = _frozen_fd_table()
    T_nu = np.asarray(T_nu, dtype=float)
    y = m / T_nu
    out = np.where(y < 1e-4, r0, 0.0)
    big = y >= 1e-4
    out = np.where(big, np.exp(sp(np.log(np.clip(y, 1e-4, 1e4)))), out)
    # beyond y = 1e4 the neutrino is pure matter: rho = m n, n = 3 zeta(3)/(4 pi^2) g T^3 per state
    nr = y > 1e4
    out = np.where(nr, y * 3 * C.ZETA3 / (4 * PI2), out)
    return g * T_nu ** 4 * out


def _cumint_from_top(u, y, p):
    """Int_{u_i}^{u_top} y(u) e^(p u) du for every grid point, summed from the top down.

    y is taken linear between grid points and the exponential is integrated exactly, so integrands
    that change by many orders of magnitude (1/H ~ T^-2 over 17 decades) lose no precision.
    """
    u0, u1, y0, y1 = u[:-1], u[1:], y[:-1], y[1:]
    h = u1 - u0
    e0, e1 = np.exp(p * u0), np.exp(p * u1)
    # Int (y0 + s (y1 - y0)/h) e^(p u) du over [u0, u1]
    i0 = (e1 - e0) / p
    i1 = ((e1 * (h - 1 / p) + e0 / p) / p)  # Int (u - u0) e^(p u) du
    pieces = y0 * i0 + (y1 - y0) / h * i1
    out = np.zeros_like(u)
    out[:-1] = np.cumsum(pieces[::-1])[::-1]
    return out


def late_universe(g_rho0: float, g_s0: float, nu_ratio: float, cosmo: dict = C.PLANCK18) -> dict:
    """Densities today (GeV^4) of matter, radiation, the massive neutrino and Lambda (flat)."""
    h = cosmo["H0"] / 100.0
    rho_crit = RHO_CRIT_H2 * h ** 2
    T0 = C.T0_GEV
    rho_m0 = (cosmo["ombh2"] + cosmo["omch2"]) * RHO_CRIT_H2
    m = cosmo["m_nu_eV"] * 1e-9
    T_nu0 = T0 * nu_ratio
    rho_nu_massive0 = rho_frozen_neutrino(T_nu0, m)
    # radiation today: photons + the two massless flavours (the massive one is counted separately)
    rho_massless_nu_per_flavour = 7 / 8 * 2 * PI2 / 30 * T_nu0 ** 4
    rho_r0 = PI2 / 30 * g_rho0 * T0 ** 4 - rho_massless_nu_per_flavour  # remove one flavour
    rho_L = rho_crit - rho_m0 - rho_r0 - float(rho_nu_massive0)
    return {"rho_crit": rho_crit, "rho_m0": rho_m0, "rho_r0": rho_r0, "rho_nu_massive0": float(rho_nu_massive0),
            "rho_L": rho_L, "Omega_m": rho_m0 / rho_crit, "Omega_r": rho_r0 / rho_crit,
            "Omega_nu": float(rho_nu_massive0) / rho_crit, "Omega_L": rho_L / rho_crit,
            "T_nu0": T_nu0, "m_nu": m, "g_s0": g_s0, "g_rho0": g_rho0, "h": h}


def calendar(E: dict, cosmo: dict = C.PLANCK18, per_decade: int = 300) -> dict:
    """Age, Hubble rate, horizons and horizon mass from today's temperature up to 10 TeV.

    ``E`` is an equation of state from ``eos.build`` (lattice model), whose grid must reach below
    ~10 keV (the end of e+- annihilation).
    """
    T_eos, g_rho, g_s = E["T"], E["g_rho"], E["g_s"]
    g_rho0, g_s0 = float(g_rho[0]), float(g_s[0])
    # neutrinos today (one effective ratio for the massive flavour: the flavour average)
    nu_ratio = float(((E["T_nu_e_over_T"][0] ** 4 + 2 * E["T_nu_mu_over_T"][0] ** 4) / 3) ** 0.25)
    L = late_universe(g_rho0, g_s0, nu_ratio, cosmo)

    T = np.logspace(np.log10(C.T0_GEV), np.log10(T_eos[-1]), int(np.log10(T_eos[-1] / C.T0_GEV) * per_decade) + 1)
    lnT = np.log(T)
    in_eos = T >= T_eos[0]
    gr = np.where(in_eos, np.interp(lnT, np.log(T_eos), g_rho), g_rho0)
    gs = np.where(in_eos, np.interp(lnT, np.log(T_eos), g_s), g_s0)
    # thermal bath: the EoS above 10 keV; below, photons + 2 massless flavours (+ the massive one below)
    rho_bath = PI2 / 30 * gr * T ** 4
    a_rel = (g_s0 / gs) ** (1 / 3) * C.T0_GEV / T                 # a / a0
    T_nu = nu_ratio * C.T0_GEV / a_rel                             # massive flavour's temperature
    rho_nu_m = rho_frozen_neutrino(T_nu, L["m_nu"])
    rho_nu_as_massless = 7 / 8 * 2 * PI2 / 30 * T_nu ** 4
    rho_bath = rho_bath - rho_nu_as_massless + rho_nu_m           # swap in the massive flavour exactly
    rho_m = L["rho_m0"] / a_rel ** 3
    rho_tot = rho_bath + rho_m + L["rho_L"]
    H = np.sqrt(8 * np.pi * rho_tot / 3) / C.M_PL                 # GeV
    dlngs = np.gradient(np.log(gs), lnT)
    integrand = (1 + dlngs / 3) / H                                # dt/d(-ln T), GeV^-1
    t_top = 1 / (2 * H[-1])                                        # constant g* above the grid
    # t ~ T^-2 and chi ~ T^-1: integrate y e^(p ln T) with the slowly varying y = integrand T^-p
    t = (t_top + _cumint_from_top(lnT, integrand * T ** 2, -2.0)) * C.HBAR_GEV_S   # seconds
    # comoving particle horizon chi = Int c dt / a (radiation era above the grid: chi = 2t/a)
    chi_top = 2 * t_top / a_rel[-1]
    chi = (chi_top + _cumint_from_top(lnT, integrand / a_rel * T, -1.0)) * C.HBAR_GEV_S * C.C_M_S
    H_s = H / C.HBAR_GEV_S
    hubble_radius_m = C.C_M_S / H_s
    particle_horizon_m = a_rel * chi
    # mass-energy inside the Hubble sphere, (4 pi/3) rho R_H^3 / c^2 = c^3 / (2 G H)
    m_hubble_kg = C.C_M_S ** 3 / (2 * C.G_SI * H_s)
    # frequency today of a wave that entered the horizon (k = aH) at T
    f_today_hz = a_rel * H_s / (2 * np.pi)
    eq = np.argmin(np.abs(np.log(rho_bath / rho_m)))
    return {"T": T, "t_s": t, "H_s": H_s, "a": a_rel, "g_rho": gr, "g_s": gs,
            "rho_bath": rho_bath, "rho_m": rho_m, "rho_L": L["rho_L"],
            "hubble_radius_m": hubble_radius_m, "particle_horizon_m": particle_horizon_m,
            "m_hubble_msun": m_hubble_kg / C.M_SUN_KG, "f_today_hz": f_today_hz,
            "late": L, "nu_ratio": nu_ratio, "i_equality": int(eq)}


def time_at(cal: dict, T) -> np.ndarray:
    """Age (s) at photon temperature T (GeV), by log-log interpolation."""
    return np.exp(np.interp(np.log(np.asarray(T, dtype=float)), np.log(cal["T"]), np.log(cal["t_s"])))


def temperature_at(cal: dict, t_s) -> np.ndarray:
    o = np.argsort(cal["t_s"])
    return np.exp(np.interp(np.log(np.asarray(t_s, dtype=float)), np.log(cal["t_s"][o]), np.log(cal["T"][o])))


def value_at(cal: dict, key: str, T) -> np.ndarray:
    return np.exp(np.interp(np.log(np.asarray(T, dtype=float)), np.log(cal["T"]), np.log(cal[key])))
