"""Why N_eff = 3.044 and not 3: neutrino decoupling, computed.

The textbook entropy argument assumes the neutrinos decoupled completely
*before* electrons and positrons annihilated. They did not, quite: at
T ~ 1–2 MeV the weak interactions were still trading a little energy between
the e± plasma and the neutrinos when annihilation began heating the plasma, so
the neutrinos caught a small share of that heat. Measured in units of one
ideally decoupled species, the neutrinos' radiation density is N_eff ≈ 3.044.

We follow the momentum-averaged method of Escudero (2019, JCAP 02, 007; 2020,
JCAP 05, 048): every species is a thermal Fermi–Dirac or Bose–Einstein gas with
its own temperature, and energy conservation plus the weak energy-transfer
rates evolve the temperatures (Escudero 2020, eqs. 3.2–3.4):

    dT_να/dt = −H T_να + (δρ_να/δt) / (dρ_να/dT_να)                 (each flavour)
    dT_γ/dt = −[4H ρ_γ + 3H (ρ_e + P_e) + 3H T dP_int/dT + δρ_νe/δt + 2 δρ_νμ/δt]
              / [dρ_γ/dT + dρ_e/dT + T d²P_int/dT²]

ν_μ and ν_τ behave identically (neither has charged-current reactions with
electrons), hence the factor 2. The energy-transfer rates are

    δρ_να/δt = 32 G_F²/π⁵ × { (X f_a1 + Y f_a2)(T_γ⁹ − T_να⁹)
                             + 7/4 (X f_s1 + Y f_s2) T_γ⁴ T_να⁴ (T_γ − T_να) + [ν–ν terms] }

with X = 4(g_L² + g_R²) and Y = 4 g_L g_R. The f's convert the exact
Maxwell–Boltzmann, massless-electron result (f = 1, 0) into Fermi–Dirac
statistics (f_a1 = 0.884, f_s1 = 0.829) and, through their dependence on
m_e/T_γ, include the electron mass (Escudero, Jackson, Laine & Sandner 2026,
whose published table we read). P_int is the finite-temperature QED
correction to the pressure of the e±γ plasma (Bennett et al. 2020):

    P₂ = −(e²T²/12π²) ∫dp (p²/E) n_D − (e²/8π⁴) [∫dp (p²/E) n_D]²
    P₃ = (e³ T / 12π⁴) I^{3/2},   I = ∫dp ((p² + E²)/E) n_D,   n_D = 2/(e^{E/T} + 1).

N_eff follows from the final energy densities (eq. 3.7):

    N_eff = (8/7) (11/4)^{4/3} (ρ_νe + 2 ρ_νμ) / ρ_γ.

Natural units (MeV, ħ = c = k_B = 1); time in MeV⁻¹.
"""

from __future__ import annotations

import numpy as np
from scipy import integrate, interpolate

import config as cfg

ME = 0.51099895          # electron mass, MeV (CODATA 2018)
GF = 1.1663788e-11       # Fermi constant, MeV^-2 (PDG 2024)
MPL = 1.220890e22        # Planck mass, MeV (non-reduced, G = 1/MPL^2)
ALPHA = 1 / 137.035999084
E2 = 4 * np.pi * ALPHA   # e^2 in Heaviside-Lorentz natural units

#: Neutral-current couplings of electrons to ν_e (including the W exchange) and
#: to ν_μ,τ, with electroweak radiative corrections for E < 10 MeV (Escudero
#: 2020, eq. 3.5). Tree level with sin²θ_W = 0.223 would give 1 ± 4s² + 8s⁴.
G_EL, G_ER, G_MUL, G_MUR = 0.727, 0.233, -0.273, 0.233


# ---------------------------------------------------------------------------
# Thermodynamics of each component
# ---------------------------------------------------------------------------


def rho_gamma(T):
    return np.pi**2 / 15.0 * T**4


def drho_gamma_dT(T):
    return 4.0 * np.pi**2 / 15.0 * T**3


def rho_nu(T):
    """One flavour, ν + ν̄, massless Fermi–Dirac: 7/8 × 2 × π²/30 T⁴."""
    return 7.0 * np.pi**2 / 120.0 * T**4


def drho_nu_dT(T):
    return 7.0 * np.pi**2 / 30.0 * T**3


def _fd_integrals(z: float) -> dict:
    """Electron–positron (g = 4) integrals at z = m_e/T, per power of T.

    With x = p/T, ε = sqrt(x² + z²) and n(ε) = 1/(e^ε + 1):
        ρ_e/T⁴ = (4/2π²) ∫ x² ε n,   P_e/T⁴ = (4/6π²) ∫ x⁴/ε n,
        (dρ_e/dT)/T³ = (4/2π²) ∫ x² ε² e^ε n²,
        J = ∫ x²/ε n   and   I = 2 ∫ (x² + ε²)/ε n   (the QED integrals, /T²).
    """
    def q(f):
        return integrate.quad(f, 0, 80, epsabs=0, epsrel=1e-12, limit=400)[0]

    def eps(x):
        return np.sqrt(x * x + z * z)

    def n(x):
        return 1.0 / (np.exp(eps(x)) + 1.0)

    g = 4.0
    return {
        "rho": g / (2 * np.pi**2) * q(lambda x: x * x * eps(x) * n(x)),
        "p": g / (6 * np.pi**2) * q(lambda x: x**4 / eps(x) * n(x)),
        "drho": g / (2 * np.pi**2) * q(lambda x: x * x * eps(x) ** 2 * n(x) * (1.0 - n(x))),
        "J": q(lambda x: x * x / eps(x) * n(x)),
        "I": 2.0 * q(lambda x: (x * x + eps(x) ** 2) / eps(x) * n(x)),
    }


class Plasma:
    """Splines in u = ln(m_e/T) of the e± thermodynamics and the QED pressure.

    The QED pressure is stored as p(u) = P_int/T⁴, from which the derivatives
    follow exactly: dP/dT = T³ (4p − p_u) and d²P/dT² = T² (12p − 7p_u + p_uu).
    """

    def __init__(self, zmin: float = 1e-3, zmax: float = 60.0, n: int = 1500):
        u = np.linspace(np.log(zmin), np.log(zmax), n)
        tab = [_fd_integrals(float(np.exp(v))) for v in u]
        self.u = u
        self.zmin, self.zmax = zmin, zmax
        col = {k: np.array([t[k] for t in tab]) for k in ("rho", "p", "drho", "J", "I")}
        self._s = {k: interpolate.CubicSpline(u, np.log(col[k])) for k in ("rho", "p", "drho")}
        J, I = col["J"], col["I"]
        p2 = -E2 / (6 * np.pi**2) * J - E2 / (2 * np.pi**4) * J**2
        p3 = E2**1.5 / (12 * np.pi**4) * I**1.5
        self._p = {2: interpolate.CubicSpline(u, p2), 3: interpolate.CubicSpline(u, p2 + p3)}
        self.massless = {"J": col["J"][0], "I": col["I"][0]}

    def _u(self, T):
        return np.log(np.clip(ME / T, self.zmin, self.zmax))

    def rho_e(self, T):
        return np.exp(self._s["rho"](self._u(T))) * T**4

    def p_e(self, T):
        return np.exp(self._s["p"](self._u(T))) * T**4

    def drho_e_dT(self, T):
        return np.exp(self._s["drho"](self._u(T))) * T**3

    def qed(self, T, order: int):
        """(P_int, dP_int/dT, d²P_int/dT²) to O(e²) (order 2) or O(e³) (order 3)."""
        if order == 0:
            return 0.0, 0.0, 0.0
        s = self._p[order]
        u = self._u(T)
        p, pu, puu = s(u), s(u, 1), s(u, 2)
        return p * T**4, T**3 * (4 * p - pu), T**2 * (12 * p - 7 * pu + puu)


# ---------------------------------------------------------------------------
# Weak-interaction energy transfer
# ---------------------------------------------------------------------------


class Rates:
    """Energy-transfer rates δρ/δt for ν_e and ν_μ (MeV⁵).

    ``mode``: "off" (instantaneous decoupling), "MB" (Maxwell–Boltzmann,
    massless electrons: f = 1), "FD" (Fermi–Dirac factors f_a1 = 0.884,
    f_s1 = 0.829, massless electrons), or "FD+me" (the published m_e/T_γ
    dependence of all four coefficients, Escudero et al. 2026).
    """

    def __init__(self, mode: str, table: np.ndarray | None = None):
        self.mode = mode
        if mode == "FD+me":
            if table is None:
                raise ValueError("FD+me needs the rate-coefficient table")
            tau = table[:, 0]
            self._f = {k: interpolate.CubicSpline(tau, table[:, j])
                       for k, j in (("a1", 1), ("a2", 2), ("s1", 5), ("s2", 6))}
            self.fa0, self.fs0 = table[0, 1], table[0, 5]
            self.tau_max = tau[-1]
        elif mode == "FD":
            self.fa0, self.fs0 = 0.884, 0.829
        else:
            self.fa0, self.fs0 = 1.0, 1.0

    def coeffs(self, T):
        if self.mode == "FD+me":
            tau = min(ME / T, self.tau_max)
            return (float(self._f["a1"](tau)), float(self._f["a2"](tau)),
                    float(self._f["s1"](tau)), float(self._f["s2"](tau)))
        return self.fa0, 0.0, self.fs0, 0.0

    def __call__(self, T, Te, Tm):
        if self.mode == "off":
            return 0.0, 0.0
        fa1, fa2, fs1, fs2 = self.coeffs(T)
        xe, ye = 4 * (G_EL**2 + G_ER**2), 4 * G_EL * G_ER
        xm, ym = 4 * (G_MUL**2 + G_MUR**2), 4 * G_MUL * G_MUR

        def plasma(x, y, Tn):
            return ((x * fa1 + y * fa2) * (T**9 - Tn**9)
                    + 1.75 * (x * fs1 + y * fs2) * T**4 * Tn**4 * (T - Tn))

        nunu = self.fa0 * (Tm**9 - Te**9) + 1.75 * self.fs0 * Tm**4 * Te**4 * (Tm - Te)
        pref = 32.0 * GF**2 / np.pi**5
        return pref * (plasma(xe, ye, Te) + 2.0 * nunu), pref * (plasma(xm, ym, Tm) - nunu)

    def relaxation_over_H(self, T, H):
        """How fast each flavour's temperature relaxes toward T_γ, in units of H.

        Γ_α = −∂(δρ_να/δt)/∂T_να / (dρ_ν/dT) at T_ν = T_γ (linear response); the
        neutrinos decouple when Γ_α falls below H.
        """
        if self.mode == "off":
            return 0.0, 0.0
        h = 1e-6 * T
        de0, dm0 = self(T, T, T)
        de1, _ = self(T, T - h, T)
        _, dm1 = self(T, T, T - h)
        return ((de1 - de0) / h) / drho_nu_dT(T) / H, ((dm1 - dm0) / h) / drho_nu_dT(T) / H


# ---------------------------------------------------------------------------
# The evolution
# ---------------------------------------------------------------------------


def solve(*, rates: str = "FD+me", qed_order: int = 3, T_start: float | None = None,
          T_end: float = 0.009, plasma: Plasma | None = None, table: np.ndarray | None = None,
          n_out: int = 600) -> dict:
    """Integrate the three temperatures from T_start to T_end (MeV), in ln a.

    Starts at 10 MeV (Escudero's choice: the weak rates exceed the expansion
    rate there by a factor of ~10⁴); runs without interactions ("instantaneous
    decoupling") start at 300 MeV so that e± are fully relativistic at
    decoupling, as the textbook argument assumes.
    """
    plasma = plasma or Plasma()
    rate = Rates(rates, table)
    if T_start is None:
        T_start = 300.0 if rates == "off" else 10.0

    def hubble(T, Te, Tm):
        P, dP, _ = plasma.qed(T, qed_order)
        rho = rho_gamma(T) + plasma.rho_e(T) + (-P + T * dP) + rho_nu(Te) + 2 * rho_nu(Tm)
        return np.sqrt(8 * np.pi * rho / 3.0) / MPL

    def rhs(lna, y):
        T, Te, Tm = y
        P, dP, d2P = plasma.qed(T, qed_order)
        H = hubble(T, Te, Tm)
        de, dm = rate(T, Te, Tm)
        num = (4 * H * rho_gamma(T) + 3 * H * (plasma.rho_e(T) + plasma.p_e(T)) + 3 * H * T * dP
               + de + 2 * dm)
        den = drho_gamma_dT(T) + plasma.drho_e_dT(T) + T * d2P
        return [-num / den / H, (-H * Te + de / drho_nu_dT(Te)) / H,
                (-H * Tm + dm / drho_nu_dT(Tm)) / H]

    def done(lna, y):
        return y[0] - T_end
    done.terminal = True

    span = (0.0, np.log(T_start / T_end) + 3.0)
    sol = integrate.solve_ivp(rhs, span, [T_start] * 3, method="LSODA", rtol=1e-10, atol=1e-13,
                              events=done, dense_output=True)
    lna = np.linspace(sol.t[0], sol.t[-1], n_out)
    T, Te, Tm = sol.sol(lna)
    k = (8.0 / 7.0) * (11.0 / 4.0) ** (4.0 / 3.0)
    track = k * (rho_nu(Te) + 2 * rho_nu(Tm)) / rho_gamma(T)
    a_end = np.exp(lna[-1]) / T_start  # comoving temperatures a T with a T_start = 1
    ratio_inst = (4.0 / 11.0) ** (1.0 / 3.0)
    gam = np.array([rate.relaxation_over_H(t, hubble(t, a, b)) for t, a, b in zip(T, Te, Tm)])
    return {
        "rates": rates, "qed_order": qed_order, "T_start": T_start,
        "neff": float(track[-1]),
        "z_gamma": float(T[-1] * a_end), "z_nue": float(Te[-1] * a_end),
        "z_numu": float(Tm[-1] * a_end),
        "Tnu_over_Tgamma_e": float(Te[-1] / T[-1]),
        "Tnu_over_Tgamma_mu": float(Tm[-1] / T[-1]),
        "ratio_vs_instantaneous_e": float(Te[-1] / T[-1] / ratio_inst),
        "ratio_vs_instantaneous_mu": float(Tm[-1] / T[-1] / ratio_inst),
        "lna": lna, "T": T, "Tnue": Te, "Tnumu": Tm, "neff_track": track,
        "gamma_e_over_H": gam[:, 0], "gamma_mu_over_H": gam[:, 1],
        "success": bool(sol.success), "nfev": int(sol.nfev),
    }


def decoupling_temperature(res: dict, which: str = "e") -> float:
    """Photon temperature (MeV) at which Γ_α/H falls through 1 (log-interpolated)."""
    g = np.asarray(res["gamma_e_over_H"] if which == "e" else res["gamma_mu_over_H"])
    T = np.asarray(res["T"])
    below = np.nonzero(g < 1.0)[0]
    if below.size == 0 or below[0] == 0:
        return float("nan")
    i = below[0]
    lg0, lg1 = np.log(g[i - 1]), np.log(g[i])
    w = lg0 / (lg0 - lg1)
    return float(np.exp(np.log(T[i - 1]) + w * (np.log(T[i]) - np.log(T[i - 1]))))


def load_rate_table() -> np.ndarray:
    """Escudero, Jackson, Laine & Sandner (2026)'s m_e-dependent rate coefficients."""
    path = cfg.RAW_DIR / "nudec" / "rate_coefficients_neutrinos.dat"
    if not path.exists():
        raise FileNotFoundError(f"{path} is not cached; run `python fetch_data.py` first")
    return np.loadtxt(path, comments="#")
