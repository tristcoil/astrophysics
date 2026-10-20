"""Neutrino decoupling and electron-positron annihilation, solved as a calculation.

Textbook version: neutrinos stop interacting before the electrons and positrons annihilate, the
annihilation heats only the photons, and entropy conservation gives T_nu/T_gamma = (4/11)^(1/3),
i.e. N_eff = 3 exactly.

Our version follows Escudero (2020, JCAP 05, 048; "NUDEC_BSM"): photons + e+- share one temperature
T_gamma, each neutrino flavour keeps a thermal (Fermi-Dirac) spectrum at its own temperature, and
energy flows between them at the weak-interaction rate

    d rho_nu_e/dt  = G_F^2/pi^5 [4 (g_eL^2 + g_eR^2) F(T_gamma, T_nu_e) + 2 F(T_nu_mu, T_nu_e)]
    d rho_nu_mu/dt = G_F^2/pi^5 [4 (g_muL^2 + g_muR^2) F(T_gamma, T_nu_mu) - F(T_nu_mu, T_nu_e)]
    F(T1, T2) = 32 f_a (T1^9 - T2^9) + 56 f_s T1^4 T2^4 (T1 - T2),

(his eqs. 3.5-3.7, massless-electron rates with Fermi-Dirac factors f_a = 0.884, f_s = 0.829).
The photon-electron plasma is not an ideal gas: finite-temperature QED changes its pressure by

    P_int = -e^2 T^4 (K/6 + K^2/2) + e^3 T^4 sqrt(2)/(6 pi) (2K + tau^2 k)^(3/2),     tau = m_e/T,

the O(e^2) and O(e^3) terms of Bennett et al. (2020, JCAP 03, 003, eqs. 3.5 and 3.17-3.18, written
with their functions K and k). Everything else is energy conservation: no fit, no free parameter.
The run starts at 10 MeV with everything in equilibrium and ends at 5 keV, when the electrons
are gone.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import CubicSpline

import config as C
import thermo
from thermo import Species

PI2 = np.pi ** 2
E2 = 4 * np.pi * C.ALPHA_EM          # e^2 in natural (Heaviside-Lorentz) units
E3 = E2 ** 1.5
ELECTRON = Species("electron", "e", "lepton", C.M_E, 4, -1, 11)
RHO_NU1 = 7 * PI2 / 120              # one flavour (nu + anti-nu) energy density / T_nu^4


# --- Finite-temperature QED --------------------------------------------------------------


def _fd(E):
    """Fermi-Dirac occupation 1/(e^E + 1), without overflow."""
    x = np.exp(-E)
    return x / (1 + x)


@lru_cache(maxsize=None)
def _qed_table():
    """phi(tau) = P_int / T^4 on a grid in u = ln tau, as splines for the O(e^2) and O(e^3) parts."""
    taus = np.exp(np.linspace(np.log(1e-4), np.log(150.0), 700))
    opts = dict(epsabs=0.0, epsrel=1e-11, limit=400)

    def K(t):
        return quad(lambda w: w * w / np.hypot(w, t) * _fd(np.hypot(w, t)), 0, np.inf, **opts)[0] / PI2

    def k(t):
        # (1/pi^2) Int dw / (E (e^E + 1)), E = sqrt(w^2 + t^2); split at w = t for the 1/E peak
        f = lambda w: _fd(np.hypot(w, t)) / np.hypot(w, t)
        return (quad(f, 0, t, **opts)[0] + quad(f, t, np.inf, **opts)[0]) / PI2

    Ks = np.array([K(t) for t in taus])
    ks = np.array([k(t) for t in taus])
    phi2 = -E2 * (Ks / 6 + Ks ** 2 / 2)
    phi3 = E3 * np.sqrt(2) / (6 * np.pi) * (2 * Ks + taus ** 2 * ks) ** 1.5
    u = np.log(taus)
    return u, CubicSpline(u, phi2), CubicSpline(u, phi3)


def qed(T, order: int = 3) -> dict:
    """P_int, dP_int/dT, d2P_int/dT2 of the e+- gamma plasma (GeV units). order 2: O(e^2) only."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    u_grid, s2, s3 = _qed_table()
    u = np.clip(np.log(C.M_E / T), u_grid[0], u_grid[-1])
    phi = s2(u) + (s3(u) if order >= 3 else 0)
    phi_u = s2(u, 1) + (s3(u, 1) if order >= 3 else 0)
    phi_uu = s2(u, 2) + (s3(u, 2) if order >= 3 else 0)
    # P = T^4 phi(tau): dP/dT = T^3 (4 phi - phi_u); d2P/dT2 = T^2 (12 phi - 7 phi_u + phi_uu)
    P = T ** 4 * phi
    dP = T ** 3 * (4 * phi - phi_u)
    d2P = T ** 2 * (12 * phi - 7 * phi_u + phi_uu)
    return {"P": P, "dP": dP, "d2P": d2P, "rho": -P + T * dP, "s": dP}


# --- Energy transfer -----------------------------------------------------------------------


def _F(T1, T2, fa=C.F_A_FD, fs=C.F_S_FD):
    return 32 * fa * (T1 ** 9 - T2 ** 9) + 56 * fs * T1 ** 4 * T2 ** 4 * (T1 - T2)


@lru_cache(maxsize=None)
def _ejls_table():
    """Escudero, Jackson, Laine & Sandner (2025): f_a1, f_a2, f_s1, f_s2 against tau = m_e/T_gamma."""
    a = np.loadtxt(C.RAW_DIR / "ejls2025_rate_coefficients.dat", comments="#")
    a = a[a[:, 0] <= 50.0]
    tau = a[:, 0]
    return {name: CubicSpline(tau, a[:, col]) for name, col in (("fa1", 1), ("fa2", 2), ("fs1", 5), ("fs2", 6))}


def _coefficients(Tg, gL, gR, rates):
    """Effective (f_a, f_s) for neutrino-electron processes of one flavour."""
    if rates == "fd":
        return C.F_A_FD, C.F_S_FD
    tab = _ejls_table()
    tau = min(C.M_E / Tg, 50.0)
    r = gL * gR / (gL ** 2 + gR ** 2)  # Y/X: the V^2 - A^2 part of the matrix element
    return float(tab["fa1"](tau) + r * tab["fa2"](tau)), float(tab["fs1"](tau) + r * tab["fs2"](tau))


def transfer_rates(Tg, Tne, Tnm, rates: str = "fd"):
    """d rho/dt for one nu_e flavour and one nu_mu (= nu_tau) flavour (GeV^5).

    rates "fd": Escudero 2020's massless-electron Fermi-Dirac coefficients (f_a = 0.884,
    f_s = 0.829); "fd_me": the electron-mass dependent coefficients of Escudero, Jackson, Laine
    & Sandner 2025 (their eqs. 2.22-2.28 with chemical potentials set to zero).
    """
    pref = C.G_FERMI ** 2 / np.pi ** 5
    fa_e, fs_e = _coefficients(Tg, C.G_EL, C.G_ER, rates)
    fa_m, fs_m = _coefficients(Tg, C.G_MUL, C.G_MUR, rates)
    if rates == "fd":
        fa_nn, fs_nn = C.F_A_FD, C.F_S_FD
    else:
        tab = _ejls_table()
        fa_nn, fs_nn = float(tab["fa1"](0.0)), float(tab["fs1"](0.0))
    nn = _F(Tnm, Tne, fa_nn, fs_nn)
    de = pref * (4 * (C.G_EL ** 2 + C.G_ER ** 2) * _F(Tg, Tne, fa_e, fs_e) + 2 * nn)
    dm = pref * (4 * (C.G_MUL ** 2 + C.G_MUR ** 2) * _F(Tg, Tnm, fa_m, fs_m) - nn)
    return de, dm


def relaxation_rates(T) -> dict:
    """Rates at which nu_e and nu_mu,tau share energy with the plasma for a small T difference.

    Gamma = (d rho_nu/dt) / (d rho_nu/dT_nu x Delta T): linearising F gives
    (288 f_a + 56 f_s) T^8 Delta T, so Gamma = G_F^2/pi^5 4 (g_L^2 + g_R^2)(288 f_a + 56 f_s) T^5 / (4 RHO_NU1 T^3).
    """
    T = np.asarray(T, dtype=float)
    lin = (288 * C.F_A_FD + 56 * C.F_S_FD) * T ** 8
    pref = C.G_FERMI ** 2 / np.pi ** 5 / (4 * RHO_NU1 * T ** 3)
    return {"nu_e": pref * 4 * (C.G_EL ** 2 + C.G_ER ** 2) * lin,
            "nu_mu": pref * 4 * (C.G_MUL ** 2 + C.G_MUR ** 2) * lin}


def hubble_equilibrium(T) -> np.ndarray:
    """H (GeV) of photons + e+- + three neutrino flavours all at temperature T, QED included."""
    T = np.asarray(T, dtype=float)
    rho = PI2 / 15 * T ** 4 + thermo.gas(ELECTRON, T)["rho"] + 3 * RHO_NU1 * T ** 4 + qed(T)["rho"]
    return np.sqrt(8 * np.pi * rho / 3) / C.M_PL


def decoupling_temperatures() -> dict:
    """Where the relaxation rate equals the expansion rate, for nu_e and nu_mu,tau (GeV)."""
    from scipy.optimize import brentq

    out = {}
    for fl in ("nu_e", "nu_mu"):
        f = lambda lt: float(np.log(relaxation_rates(np.exp(lt))[fl] / hubble_equilibrium(np.exp([lt]))[0]))
        out[fl] = float(np.exp(brentq(f, np.log(1e-4), np.log(0.1))))
    return out


# --- The decoupling run -----------------------------------------------------------------------


def solve(collisions: bool = True, qed_order: int = 3, separate: bool = True, rates: str = "fd_me",
          T_start: float = C.T_NU_START, T_end: float = C.T_NU_END, rate_scale: float = 1.0,
          n_out: int = 1500) -> dict:
    """Integrate T_gamma, T_nu_e, T_nu_mu and time over ln a from T_start to T_end.

    ``collisions=False`` is instantaneous decoupling at T_start; ``qed_order`` 0 turns the QED
    corrections off, 2 keeps O(e^2), 3 adds O(e^3); ``separate=False`` gives all three flavours
    one common temperature (mimicking flavour oscillations, as Escudero's headline 3.045);
    ``rates`` "fd" (Escudero 2020) or "fd_me" (with the electron mass, EJLS 2025).
    """

    def plasma(Tg):
        el = thermo.gas(ELECTRON, Tg)
        q = qed(Tg, qed_order) if qed_order else {"P": 0 * Tg, "dP": 0 * Tg, "d2P": 0 * Tg, "rho": 0 * Tg}
        return el, q

    def rhs(lna, y):
        Tg, Tne, Tnm = y[0], y[1], y[2]
        el, q = plasma(np.array([Tg]))
        rho_g = PI2 / 15 * Tg ** 4
        rho_tot = rho_g + el["rho"][0] + q["rho"][0] + RHO_NU1 * (Tne ** 4 + 2 * Tnm ** 4)
        H = np.sqrt(8 * np.pi * rho_tot / 3) / C.M_PL
        if collisions:
            de, dm = transfer_rates(Tg, Tne, Tnm, rates)
            de, dm = rate_scale * de, rate_scale * dm
            if not separate:
                tot = (de + 2 * dm) / 3
                de = dm = tot
        else:
            de = dm = 0.0
        num = 4 * rho_g + 3 * (el["rho"][0] + el["p"][0]) + 3 * Tg * q["dP"][0] + (de + 2 * dm) / H
        den = 4 * PI2 / 15 * Tg ** 3 + el["drho_dT"][0] + Tg * q["d2P"][0]
        dTg = -num / den
        dTne = -Tne + de / (H * 4 * RHO_NU1 * Tne ** 3)
        dTnm = -Tnm + dm / (H * 4 * RHO_NU1 * Tnm ** 3)
        return [dTg, dTne, dTnm, 1.0 / H]

    # ln a runs from 0 (T = T_start) to about ln(T_start/T_end) + ln(11/4)/3 + margin
    span = np.log(T_start / T_end) + 0.5
    lna = np.linspace(0, span, n_out)
    T0 = T_start
    H0 = hubble_equilibrium(np.array([T0]))[0]
    t0 = 1 / (2 * float(H0))  # radiation-dominated age at the start (Escudero's choice)
    sol = solve_ivp(rhs, (0, span), [T0, T0, T0, t0], t_eval=lna, method="LSODA",
                    rtol=1e-10, atol=1e-16)
    if not sol.success:
        raise RuntimeError(sol.message)
    Tg, Tne, Tnm, t = sol.y
    keep = Tg >= T_end * 0.999
    Tg, Tne, Tnm, t, lna = Tg[keep], Tne[keep], Tnm[keep], t[keep], lna[keep]
    rg_e, rg_m = Tg[-1] / Tne[-1], Tg[-1] / Tnm[-1]
    n_eff = (11 / 4) ** (4 / 3) * ((Tne[-1] / Tg[-1]) ** 4 + 2 * (Tnm[-1] / Tg[-1]) ** 4)
    # comoving neutrino temperature in units of the initial one, as Escudero's z_nu
    return {"lna": lna, "T_gamma": Tg, "T_nu_e": Tne, "T_nu_mu": Tnm, "t_s": t * C.HBAR_GEV_S,
            "ratio_e": rg_e, "ratio_mu": rg_m, "N_eff": n_eff,
            "ratio_common": Tg[-1] / ((Tne[-1] ** 4 + 2 * Tnm[-1] ** 4) / 3) ** 0.25,
            "settings": {"collisions": collisions, "qed_order": qed_order, "separate": separate,
                         "rates": rates, "T_start_MeV": T_start * 1e3, "rate_scale": rate_scale}}


def low_T_gstar(run: dict) -> dict:
    """g*rho and g*s today (photons + neutrinos) from a decoupling run's final temperatures."""
    re_, rm_ = 1 / run["ratio_e"], 1 / run["ratio_mu"]
    g_rho = 2 + 7 / 8 * 2 * (re_ ** 4 + 2 * rm_ ** 4)
    g_s = 2 + 7 / 8 * 2 * (re_ ** 3 + 2 * rm_ ** 3)
    return {"g_rho": g_rho, "g_s": g_s}


def instantaneous_ratio(T, T_dec: float | None = None) -> np.ndarray:
    """T_nu/T_gamma after sudden decoupling, from entropy conservation in the ideal e+- gamma gas.

    T_dec None: decoupling before any annihilation (the textbook (4/11)^(1/3) today).
    """
    T = np.asarray(T, dtype=float)

    def g_s_eg(TT):
        return 2 + thermo.gstar(ELECTRON, TT)["s"]

    g_dec = 2 + 7 / 8 * 4 if T_dec is None else float(g_s_eg(np.array([T_dec]))[0])
    r = (g_s_eg(T) / g_dec) ** (1 / 3)
    if T_dec is not None:
        r = np.where(T >= T_dec, 1.0, r)
    return r
