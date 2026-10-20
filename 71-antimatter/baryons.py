"""Counting baryons and photons: the one-in-a-billion number.

The baryon-to-photon ratio eta = n_b / n_gamma is the cleanest way to state
how much matter survived. Two measurements give it: the deuterium made in the
first minutes (``bbn.py``) and the CMB's acoustic peaks, which measure the
baryon density as Omega_b h^2. This module converts between the two and works
out what the numbers mean:

* n_gamma = (2 zeta(3) / pi^2) (k T0 / hbar c)^3, the CMB photons per m^3;
* n_b = Omega_b h^2 rho_crit,100 / m_b, where m_b is the mean mass per baryon
  *including electrons and nuclear binding* (atomic masses of 1H and 4He), so
  the conversion depends slightly on the helium fraction Y_p;
* eta_10 = 10^10 eta ~ 273.3 Omega_b h^2 / (1 - 0.0071 Y_p).

It also follows the excess back in time. What is conserved as the universe
expands (without baryon-number violation) is the baryon number per unit
entropy, Y_B = n_B / s. Early on there were nearly as many antiquarks as
quarks, and the excess per quark is Y_B s / n_q, which we compute with the
Standard Model's entropy degrees of freedom g*s(T) (Saikawa & Shirai 2018).
"""

from __future__ import annotations

import math

import numpy as np

import config as cfg

#: Critical density for H0 = 100 km/s/Mpc (kg/m^3).
RHO_CRIT_100 = 3 * (100e3 / cfg.MPC) ** 2 / (8 * math.pi * cfg.G_NEWTON)


def n_gamma(t_cmb: float = cfg.T_CMB) -> float:
    """CMB photon number density today (per m^3)."""
    return 2 * cfg.ZETA3 / math.pi**2 * (cfg.K_B * t_cmb / (cfg.HBAR * cfg.C_LIGHT)) ** 3


def mass_per_baryon(yp: float = cfg.YP) -> float:
    """Mean mass per baryon in kg, from atomic masses (electrons included)."""
    return cfg.M_U_KG * (cfg.M_H1_U * (1 - yp) + cfg.M_HE4_U / 4 * yp)


def n_baryon_from_ombh2(ombh2: float, yp: float = cfg.YP) -> float:
    """Baryons per m^3 today."""
    return ombh2 * RHO_CRIT_100 / mass_per_baryon(yp)


def eta_from_ombh2(ombh2, yp: float = cfg.YP, t_cmb: float = cfg.T_CMB):
    return np.asarray(ombh2) * RHO_CRIT_100 / mass_per_baryon(yp) / n_gamma(t_cmb)


def ombh2_from_eta(eta, yp: float = cfg.YP, t_cmb: float = cfg.T_CMB):
    return np.asarray(eta) * n_gamma(t_cmb) * mass_per_baryon(yp) / RHO_CRIT_100


def eta10_coefficient(yp: float = cfg.YP, t_cmb: float = cfg.T_CMB) -> float:
    """eta_10 per unit Omega_b h^2 (273.3 / (1 - 0.0071 Y_p) at T0 = 2.7255 K)."""
    return 1e10 * float(eta_from_ombh2(1.0, yp, t_cmb))


def entropy_per_photon(gstar_s: float = 3.931) -> float:
    """s / n_gamma = (2 pi^4 / (45 zeta(3))) g*s / 2: 7.04 today."""
    return math.pi**4 / (45 * cfg.ZETA3) * gstar_s


def numbers_today(ombh2: float, t_cmb: float = cfg.T_CMB, yp: float = cfg.YP) -> dict:
    """The asymmetry in everyday units."""
    ng = n_gamma(t_cmb)
    nb = n_baryon_from_ombh2(ombh2, yp)
    eta = nb / ng
    gs0 = 3.931                     # Saikawa & Shirai's g*s today (photons + neutrinos)
    return {
        "n_gamma_per_cm3": ng * 1e-6,
        "n_b_per_m3": nb,
        "protons_per_m3_if_all_hydrogen": nb,
        "eta": eta,
        "photons_per_baryon": 1 / eta,
        "Y_B": eta / entropy_per_photon(gs0),
        "gstar_s_today": gs0,
        "rho_b_kg_m3": nb * mass_per_baryon(yp),
        "mass_per_baryon_u": mass_per_baryon(yp) / cfg.M_U_KG,
    }


# ---------------------------------------------------------------------------
# The excess in the early universe
# ---------------------------------------------------------------------------

#: Quarks: (name, mass in GeV for the thermal count, g = 2 spins x 3 colours).
#: Masses: PDG MS-bar values for u, d, s (at 2 GeV), m_c(m_c), m_b(m_b), and
#: the top-quark mass -- only used to decide when each flavour is
#: relativistic, so scheme differences do not matter here.
QUARKS = [("u", 0.00216), ("d", 0.00470), ("s", 0.0935), ("c", 1.273), ("b", 4.183),
          ("t", 172.57)]


def fermion_number_density(T_gev, m_gev, g: float):
    """n for one fermion species (particles only, no chemical potential),
    in GeV^3: g/(2 pi^2) int p^2 dp / (exp(E/T) + 1)."""
    T = np.atleast_1d(np.asarray(T_gev, float))
    out = np.empty_like(T)
    x = np.linspace(0, 60, 6001)                         # p / T
    for i, t in enumerate(T):
        e = np.sqrt(x**2 + (m_gev / t) ** 2)
        out[i] = g / (2 * math.pi**2) * t**3 * np.trapezoid(x**2 / (np.exp(np.minimum(e, 700)) + 1), x)
    return out


def quark_excess(T_gev, gstar_s, y_b: float):
    """Fraction of quarks that are 'extra', (n_q - n_qbar) / (n_q + n_qbar),
    above the QCD crossover: 3 Y_B s / (2 sum_flavours n_q)."""
    T = np.asarray(T_gev, float)
    s = 2 * math.pi**2 / 45 * np.asarray(gstar_s) * T**3
    nq = sum(fermion_number_density(T, m, 6.0) for _, m in QUARKS)
    return 3 * y_b * s / (2 * nq)


# ---------------------------------------------------------------------------
# The annihilation catastrophe: what a perfectly symmetric universe leaves
# ---------------------------------------------------------------------------

M_N_GEV = 0.93892                     # mean nucleon mass
G_NUCLEON = 4.0                       # p and n, two spin states each
#: <sigma_ann v> for nucleon-antinucleon annihilation at T ~ 10-60 MeV:
#: beta * sigma_ann = 44 +- 4 mb from 40 to 800 MeV/c (OBELIX; Zhou &
#: Timmermans 2012), i.e. sigma_ann * v_rel ~ 1.32e-15 cm^3/s, almost constant
#: (the s-wave 1/v law). Isospin partners (n-bar p, n-bar n) assumed equal.
SIGMA_V_CM3S = 44e-27 * 2.99792458e10
GEV2_TO_CM3S = 0.389379e-27 * 2.99792458e10          # 1 GeV^-2 times c


def nucleon_freezeout(sigma_v_cm3s: float = SIGMA_V_CM3S, x_end: float = 400.0, eta: float = 0.0):
    """Nucleon and antinucleon abundances Y = n/s as the universe cools.

    eta = 0: perfectly symmetric (the "annihilation catastrophe");
    otherwise the real universe with net baryon number Y_B = eta n_gamma/s.
    Solved in x = m_N/T from T = 300 MeV with Saikawa & Shirai's g*(T).
    Returns x, T (GeV), Y_N, Y_Nbar, Y_eq."""
    from scipy.integrate import solve_ivp

    import fetch_data as fd

    tt, grho, gs = fd.load_gstar()
    lt = np.log(tt)

    def g_rho(T):
        return np.interp(np.log(T), lt, grho)

    def g_s(T):
        return np.interp(np.log(T), lt, gs)

    def dlngs_dlnT(T):
        h = 0.01
        return (math.log(g_s(T * math.exp(h))) - math.log(g_s(T * math.exp(-h)))) / (2 * h)

    sv = sigma_v_cm3s / GEV2_TO_CM3S                       # GeV^-2
    y_b = eta / entropy_per_photon() if eta > 0 else 0.0

    def y_eq(x):
        T = M_N_GEV / x
        n = G_NUCLEON * (M_N_GEV * T / (2 * math.pi)) ** 1.5 * math.exp(-x)
        return n / (2 * math.pi**2 / 45 * g_s(T) * T**3)

    def lam(x):
        T = M_N_GEV / x
        s = 2 * math.pi**2 / 45 * g_s(T) * T**3
        H = math.sqrt(8 * math.pi**3 * g_rho(T) / 90) * T * T / cfg.M_PLANCK_GEV
        return (s * sv / (H * x)) * (1 - dlngs_dlnT(T) / 3)

    def rhs(x, y):
        # Y_Nbar = y[0], Y_N = Y_Nbar + Y_B:
        # dY/dx = -(s <sv> / (H x)) (1 - dln g*s/dln T / 3) (Y_N Y_Nbar - Y_eq^2)
        ye = y_eq(x)
        return [-lam(x) * ((y[0] + y_b) * y[0] - ye * ye)]

    def jac(x, y):
        return [[-lam(x) * (2 * y[0] + y_b)]]

    x0 = M_N_GEV / 0.3
    ye0 = y_eq(x0)
    ybar0 = (-y_b + math.sqrt(y_b * y_b + 4 * ye0 * ye0)) / 2
    xs = np.geomspace(x0, x_end, 600)
    sol = solve_ivp(rhs, (x0, x_end), [ybar0], method="Radau", jac=jac, t_eval=xs,
                    rtol=1e-9, atol=1e-40)
    ybar = np.maximum(sol.y[0], 0.0)
    yeq = np.array([y_eq(x) for x in xs])
    return {"x": xs, "T": M_N_GEV / xs, "Y_N": ybar + y_b, "Y_Nbar": ybar, "Y_eq": yeq,
            "Y_B": y_b}
