"""The universe's own passage through the QCD crossover.

1. **When.** The age of the universe at each temperature comes from project
   P2's calendar (its full Standard-Model equation of state, ``t(T)`` in
   ``68-thermal-history/data/derived/curves.npz``). As a cross-check we also
   integrate the Friedmann equation here with the lattice equation of state
   (HotQCD 2014 ansatz, our hadron gas below 130 MeV) plus photons, leptons
   and neutrinos: with entropy conserved, d ln a = -(1/3) d ln s, so
   t(T) = t(T_max) + int (1/(3H)) d ln s, with H^2 = 8 pi G rho / 3.

2. **How close to mu_B = 0.** The universe's matter excess is fixed: one
   baryon per ~1.2e10 units of entropy (Planck's omega_b). Near mu = 0 every
   charge density is linear in the chemical potentials,
   n_X / T^3 = sum_Y chi_XY mu_Y / T, so the conditions
   n_B = b s, n_S = 0, n_Q = 0 (leptons included) and n_L = l s per flavour
   give a 6x6 linear system for mu_B, mu_Q, mu_S and the three lepton
   chemical potentials. The strong-interaction susceptibilities chi_XY come
   from our hadron gas below 145 MeV and from lattice QCD (HotQCD 2012) above.
   The lepton asymmetry l is not measured; we show l = 0 and the value
   sphalerons would give, l = -(51/28) b shared by three flavours.
"""

from __future__ import annotations

import numpy as np
from numpy.polynomial.laguerre import laggauss
from scipy.special import zeta

import config as cfg
import hrg
import lattice

ZETA3 = float(zeta(3))
LEPTONS = {"e": 0.51099895e-3, "mu": 0.1056583755, "tau": 1.77686}   # GeV
_X, _W = laggauss(64)
_W = _W * np.exp(_X)


def eta_from_omega_b(omega_b: float = cfg.PLANCK_OMEGA_B[0]) -> float:
    """Baryon-to-photon ratio from the physical baryon density."""
    rho_c_h2 = 1.87847e-26            # kg/m^3, critical density for h = 1
    m_p = 1.67262192e-27              # kg
    n_b = omega_b * rho_c_h2 / m_p    # m^-3
    kT_over_hc = cfg.T_CMB_K * 8.617333262e-5 / 1.97326980e-7   # 1/m
    n_gamma = 2 * ZETA3 / np.pi ** 2 * kT_over_hc ** 3
    return n_b / n_gamma


def baryon_per_entropy(g_s_today: float = 3.9302) -> float:
    """n_B / s, conserved since baryogenesis (g*s today from P2: 3.9302)."""
    s_over_ngamma = (2 * np.pi ** 2 / 45 * g_s_today) / (2 * ZETA3 / np.pi ** 2)
    return eta_from_omega_b() / s_over_ngamma


def _fermion(T: float, m: float, g: float) -> dict:
    """Energy, pressure and net-charge susceptibility of a fermion species
    (particles + antiparticles, g internal states each), over T^4 / T^2."""
    if m == 0:
        return {"e": g * 2 * 7 / 8 * np.pi ** 2 / 30, "p": g * 2 * 7 / 8 * np.pi ** 2 / 90,
                "chi": g / 6}
    x = _X                            # p / T
    E = np.sqrt(x * x + (m / T) ** 2)
    f = 1 / (np.exp(E) + 1)
    pref = 2 * g / (2 * np.pi ** 2)   # particle + antiparticle
    e = pref * np.sum(_W * x * x * E * f)
    p = pref * np.sum(_W * x ** 4 / (3 * E) * f)
    chi = pref * np.sum(_W * x * x * f * (1 - f))
    return {"e": e, "p": p, "chi": chi}


def leptons_and_photons(T: float) -> dict:
    """e+p over T^4 for photons, charged leptons and three neutrino species
    (all in equilibrium: valid above ~2 MeV), and lepton susceptibilities."""
    e = 2 * np.pi ** 2 / 30
    p = 2 * np.pi ** 2 / 90
    chi = {}
    for name, m in LEPTONS.items():
        r = _fermion(T, m, 2)
        e, p = e + r["e"], p + r["p"]
        chi[name] = r["chi"]
    nu = _fermion(T, 0.0, 1)
    e, p = e + 3 * nu["e"], p + 3 * nu["p"]
    return {"e": e, "p": p, "chi_l": chi, "chi_nu": nu["chi"]}


# ---------------------------------------------------------------------------
# The strong sector's equation of state (mu = 0)
# ---------------------------------------------------------------------------

def qcd_pressure(T_gev: np.ndarray) -> np.ndarray:
    """p/T^4 of QCD: our hadron gas below 120 MeV, the HotQCD ansatz above
    140 MeV, blended smoothly in between (they agree to 1% at 130 MeV)."""
    gas = hrg.model()
    T = np.atleast_1d(np.asarray(T_gev, float))
    out = np.empty_like(T)
    for i, t in enumerate(T):
        Tm = t * 1e3
        lat = float(lattice.hotqcd_pressure(np.array([Tm]))[0]) if Tm > 100 else 0.0
        if Tm >= 140:
            out[i] = lat
            continue
        tot = gas.totals(gas.state(t))
        h = tot["p"] * cfg.HBARC ** 3 / t ** 4
        w = np.clip((Tm - 120) / 20, 0, 1)
        w = w * w * (3 - 2 * w)
        out[i] = (1 - w) * h + w * lat
    return out


def friedmann_age(T_eval_mev: np.ndarray, T_max_gev: float = 2.0) -> np.ndarray:
    """Age (s) at temperatures T_eval from our own integration (cross-check of P2)."""
    lnT = np.linspace(np.log(T_max_gev), np.log(0.010), 1500)
    T = np.exp(lnT)
    pq = qcd_pressure(T)
    # s/T^3 = 4 p/T^4 + T d(p/T^4)/dT, e/T^4 = 3 p/T^4 + T d(p/T^4)/dT
    dpq = np.gradient(pq, lnT)
    eq, sq = 3 * pq + dpq, 4 * pq + dpq
    lep = [leptons_and_photons(t) for t in T]
    e = (eq + np.array([x["e"] for x in lep])) * T ** 4          # GeV^4
    s = (sq + np.array([x["e"] + x["p"] for x in lep])) * T ** 3  # GeV^3
    H = np.sqrt(8 * np.pi / 3 * e) / cfg.M_PLANCK                  # GeV
    lns = np.log(s)
    t0 = 1 / (2 * H[0])                       # radiation era, slowly varying g*
    dt = -np.diff(lns) / 3 * 0.5 * (1 / H[1:] + 1 / H[:-1])
    t_gev = t0 + np.concatenate([[0.0], np.cumsum(dt)])
    t_s = t_gev * cfg.HBAR_GEV_S
    Te = np.asarray(T_eval_mev, float) / 1e3
    return np.exp(np.interp(-np.log(Te), -lnT, np.log(t_s)))


def p2_calendar() -> dict | None:
    """P2's t(T) and g*(T), if that project has been run."""
    path = cfg.P2_DIR / "data" / "derived" / "curves.npz"
    if not path.exists():
        return None
    d = np.load(path)
    return {"T_gev": d["cal_T"], "t_s": d["cal_t"], "g_T": d["T"], "g_rho": d["g_rho"],
            "g_s": d["g_s"]}


def age_at(T_mev, cal: dict | None = None) -> np.ndarray:
    """Age of the universe (s) at temperature(s) T, from P2's calendar
    (falls back to our own integration if P2 has not been run)."""
    cal = p2_calendar() if cal is None else cal
    T = np.atleast_1d(np.asarray(T_mev, float)) / 1e3
    if cal is None:
        return friedmann_age(T * 1e3)
    order = np.argsort(cal["T_gev"])
    return np.exp(np.interp(np.log(T), np.log(cal["T_gev"][order]),
                            np.log(cal["t_s"][order])))


# ---------------------------------------------------------------------------
# The universe's baryon chemical potential
# ---------------------------------------------------------------------------

def strong_susceptibilities(T_gev: float) -> dict:
    """chi_XY / T^2 of the strong sector: hadron gas below 145 MeV, HotQCD 2012
    (150-250 MeV) above, blended over 145-155 MeV."""
    Tm = T_gev * 1e3
    gas = hrg.model()
    lat = lattice.hotqcd_2012_susceptibilities()
    keys = ("BB", "QQ", "SS", "BS", "BQ", "QS")
    if Tm > lat["T"][-1]:
        raise ValueError("lattice table ends at 250 MeV")
    h = gas.susceptibilities(T_gev) if Tm < 155 else None
    l_ = {k: float(np.interp(Tm, lat["T"], lat[k])) for k in keys} if Tm > 145 else None
    if h is None:
        return l_
    h = {k: h[k] for k in keys}
    if l_ is None:
        return h
    w = (Tm - 145) / 10
    return {k: (1 - w) * h[k] + w * l_[k] for k in keys}


def chemical_potentials(T_gev: float, lepton: str = "zero") -> dict:
    """mu/T of baryon number, charge, strangeness and lepton flavours."""
    b = baryon_per_entropy()
    l_each = 0.0 if lepton == "zero" else -(51 / 28) * b / 3
    chi = strong_susceptibilities(T_gev)
    lp = leptons_and_photons(T_gev)
    s_qcd = 4 * qcd_pressure(np.array([T_gev]))[0]
    # s/T^3 including d(p/T^4)/dlnT for QCD
    h = 1e-3
    dp = (qcd_pressure(np.array([T_gev * (1 + h)]))[0]
          - qcd_pressure(np.array([T_gev * (1 - h)]))[0]) / (2 * h)
    s_hat = s_qcd + dp + lp["e"] + lp["p"]
    # unknowns: xB, xQ, xS, xLe, xLmu, xLtau ; charged lepton mu = mu_L - mu_Q
    A = np.zeros((6, 6))
    rhs = np.zeros(6)
    A[0, :3] = [chi["BB"], chi["BQ"], chi["BS"]]
    rhs[0] = b * s_hat
    A[1, :3] = [chi["BS"], chi["QS"], chi["SS"]]
    A[2, :3] = [chi["BQ"], chi["QQ"], chi["QS"]]
    for k, name in enumerate(LEPTONS):
        cl = lp["chi_l"][name]
        # charge: charged lepton (Q = -1) net density cl (xL - xQ)
        A[2, 1] += cl
        A[2, 3 + k] -= cl
        A[3 + k, 1] = -cl
        A[3 + k, 3 + k] = cl + lp["chi_nu"]
        rhs[3 + k] = l_each * s_hat
    x = np.linalg.solve(A, rhs)
    return {"T_MeV": T_gev * 1e3, "muB_over_T": x[0], "muQ_over_T": x[1], "muS_over_T": x[2],
            "muL_over_T": x[3:].tolist(), "s_over_T3": s_hat, "b": b,
            "proton_excess": float(np.tanh(x[0] + x[1]))}


def entropy_over_T3(T_gev: float) -> float:
    """Total entropy density / T^3 (QCD + photons + leptons + neutrinos)."""
    h = 1e-3
    p = qcd_pressure(np.array([T_gev * (1 - h), T_gev, T_gev * (1 + h)]))
    lp = leptons_and_photons(T_gev)
    return 4 * p[1] + (p[2] - p[0]) / (2 * h) + lp["e"] + lp["p"]


def mu_b_nonlinear(T_gev: float) -> float:
    """mu_B (GeV) solving n_B(hadron gas) = b s exactly, for the cool end of
    the path where antibaryons annihilate and mu_B/T is no longer small
    (mu_S from strangeness neutrality; mu_Q ~ 1e-9 T is neglected). Nuclei are
    left out: below ~1 MeV they would condense into helium (that is BBN)."""
    from scipy.optimize import brentq

    gas = hrg.model(with_nuclei=False)
    target = baryon_per_entropy() * entropy_over_T3(T_gev) * T_gev ** 3 / cfg.HBARC ** 3

    def f(mu):
        st = gas.solve_mu(T_gev, mu, q_over_b=None)
        return st.n @ gas.hl.B - target

    hi = 0.938272 - 1e-6              # below the proton mass
    return brentq(f, 0.0, hi, xtol=1e-15, rtol=1e-12)


def path(T_mev: np.ndarray, lepton: str = "zero") -> dict:
    """The universe's (T, mu_B) path: linear response with lattice
    susceptibilities at T >= 150 MeV, the exact hadron-gas solution below."""
    out_mu, out_x = [], []
    for Tm in np.asarray(T_mev, float):
        if Tm >= 150:
            x = chemical_potentials(Tm / 1e3, lepton)["muB_over_T"]
            mu = x * Tm
        else:
            mu = mu_b_nonlinear(Tm / 1e3) * 1e3
            x = mu / Tm
        out_mu.append(mu)
        out_x.append(x)
    return {"T_MeV": np.asarray(T_mev, float), "muB_MeV": np.array(out_mu),
            "muB_over_T": np.array(out_x)}
