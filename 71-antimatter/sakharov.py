"""Why the Standard Model cannot make the matter excess: Sakharov's checklist.

Sakharov (1967) showed that a universe that starts symmetric can only end up
with more matter than antimatter if three things happen:

1. **Baryon number is violated.** In the Standard Model it is, through
   "sphalerons": transitions over an energy barrier of the electroweak
   fields that change baryon and lepton number together. They are fast in
   the hot early universe and switch off as the Higgs field turns on. We use
   the lattice rate of D'Onofrio, Rummukainen & Tranberg (2014) and compare
   it with the expansion rate H(T) to find when they stop.
2. **C and CP are violated.** The Standard Model's only CP violation (in the
   quark sector) is the phase of the CKM matrix, measured by the Jarlskog
   invariant J ~ 3e-5. But it only acts if all quark masses of each charge
   differ; Jarlskog's dimensionless measure at the electroweak temperature,
   J times the product of squared-mass differences over T^12, is ~1e-20,
   against the ~1e-10 needed.
3. **Departure from thermal equilibrium.** The expansion at the electroweak
   epoch is 16 orders of magnitude slower than particle reactions, so the
   only chance is a violent (first-order) phase transition with bubbles of
   the new phase. In the Standard Model with a 125 GeV Higgs there is none:
   lattice simulations find a smooth crossover (first order only if m_H were
   below ~72 GeV). Even the one-loop estimate of the transition strength,
   v(T_c)/T_c = 2E/lambda, gives ~0.15 for 125 GeV, far below the ~1 needed
   to stop sphalerons from washing out whatever was made.
"""

from __future__ import annotations

import math

import numpy as np
from scipy import optimize

import config as cfg

# Electroweak inputs (GeV): PDG values, filled in config.PDG_EW.
GF = 1.1663788e-5
V_EW = (math.sqrt(2) * GF) ** -0.5          # 246.22 GeV


def hubble_rate_gev(T, gstar):
    """H = sqrt(8 pi^3 g* / 90) T^2 / M_Pl (radiation era), in GeV."""
    return math.sqrt(8 * math.pi**3 / 90) * np.sqrt(gstar) * np.asarray(T) ** 2 / cfg.M_PLANCK_GEV


# ---------------------------------------------------------------------------
# 3. The strength of the electroweak transition at one loop
# ---------------------------------------------------------------------------


def one_loop_coefficients(m_h, m_w, m_z, m_t, v=V_EW):
    """High-temperature one-loop potential V = D (T^2 - T0^2) phi^2 - E T phi^3
    + lambda_T phi^4 / 4 (Anderson & Hall 1992; Dine et al. 1992)."""
    D = (2 * m_w**2 + m_z**2 + 2 * m_t**2) / (8 * v * v)
    E = (2 * m_w**3 + m_z**3) / (4 * math.pi * v**3)
    B = 3.0 / (64 * math.pi**2 * v**4) * (2 * m_w**4 + m_z**4 - 4 * m_t**4)
    T0sq = (m_h**2 - 8 * B * v * v) / (4 * D)
    lam = m_h**2 / (2 * v * v)
    ln_ab = 2 * math.log(4 * math.pi) - 2 * 0.5772156649
    ln_af = 2 * math.log(math.pi) - 2 * 0.5772156649

    def lam_T(T):
        return lam - 3.0 / (16 * math.pi**2 * v**4) * (
            2 * m_w**4 * (math.log(m_w**2 / T**2) - ln_ab)
            + m_z**4 * (math.log(m_z**2 / T**2) - ln_ab)
            - 4 * m_t**4 * (math.log(m_t**2 / T**2) - ln_af))

    return {"D": D, "E": E, "B": B, "T0sq": T0sq, "lambda": lam, "lambda_T": lam_T}


def transition_strength(m_h, m_w=80.37, m_z=91.1876, m_t=172.57):
    """(T_c, v(T_c)/T_c) of the one-loop high-T potential: degenerate minima
    when D (T_c^2 - T0^2) = E^2 T_c^2 / lambda_Tc, and v_c = 2 E T_c / lambda_Tc."""
    c = one_loop_coefficients(m_h, m_w, m_z, m_t)
    if c["T0sq"] <= 0:
        return float("nan"), float("nan")

    def f(T):
        return c["D"] * (T * T - c["T0sq"]) - c["E"] ** 2 * T * T / c["lambda_T"](T)

    t0 = math.sqrt(c["T0sq"])
    grid = t0 * np.geomspace(1 + 1e-9, 4.0, 400)
    vals = np.array([f(t) for t in grid])
    ok = np.where(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
    if len(ok) == 0 or c["lambda_T"](grid[ok[0]]) <= 0:
        return float("nan"), float("nan")
    tc = optimize.brentq(f, grid[ok[0]], grid[ok[0] + 1])
    return tc, 2 * c["E"] / c["lambda_T"](tc)


# ---------------------------------------------------------------------------
# 2. CP violation: the Jarlskog invariant and Jarlskog's measure
# ---------------------------------------------------------------------------


def jarlskog_from_wolfenstein(lam, A, rhobar, etabar):
    """J = Im(V_us V_cb V_ub* V_cs*) from the Wolfenstein parameters, using the
    exact (all-order) standard parametrisation (PDG CKM review, eqs. 12.4-12.5)."""
    s12 = lam
    s23 = A * lam * lam
    # rho + i eta = s13 e^{i delta} / (s12 s23) up to the PDG all-order relation
    z = complex(rhobar, etabar)
    s13e = A * lam**3 * z * math.sqrt(1 - A * A * lam**4) / (
        math.sqrt(1 - lam * lam) * (1 - A * A * lam**4 * z))
    s13 = abs(s13e)
    delta = math.atan2(s13e.imag, s13e.real)
    c12, c23, c13 = (math.sqrt(1 - x * x) for x in (s12, s23, s13))
    return c12 * c23 * c13**2 * s12 * s23 * s13 * math.sin(delta), {
        "s12": s12, "s23": s23, "s13": s13, "delta": delta}


def ckm_matrix(s12, s23, s13, delta):
    c12, c23, c13 = (math.sqrt(1 - x * x) for x in (s12, s23, s13))
    e = complex(math.cos(delta), math.sin(delta))
    return np.array([
        [c12 * c13, s12 * c13, s13 / e],
        [-s12 * c23 - c12 * s23 * s13 * e, c12 * c23 - s12 * s23 * s13 * e, s23 * c13],
        [s12 * s23 - c12 * c23 * s13 * e, -c12 * s23 - s12 * c23 * s13 * e, c23 * c13]])


def jarlskog_measure(J, up_masses, down_masses, T):
    """Jarlskog's dimensionless CP measure for baryogenesis:
    J prod_{i<j} (m_ui^2 - m_uj^2) prod_{i<j} (m_di^2 - m_dj^2) / T^12."""
    mu, md = np.asarray(up_masses, float), np.asarray(down_masses, float)
    prod = 1.0
    for m in (mu, md):
        for i in range(3):
            for j in range(i + 1, 3):
                prod *= abs(m[i] ** 2 - m[j] ** 2)
    return J * prod / np.asarray(T, float) ** 12
