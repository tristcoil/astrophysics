"""The relic neutrinos today: how many, how cold, how fast.

One second after the Big Bang the weak interactions froze out and the
neutrinos stopped trading energy with everything else. From then on each
neutrino's momentum simply stretched with the expanding universe, p ∝ 1/a, so
the population kept the Fermi–Dirac *shape* it had at decoupling, at a
temperature falling as 1/a. That is true even after the neutrinos became
non-relativistic: the occupation number f(p) = 1/(e^(p/T_ν) + 1) is frozen in
comoving momentum. What changed later is only what each momentum means: for a
massive neutrino, speed v = p c / sqrt(p² + m² c²).

The photons share the same 1/a cooling, except for one event: when the
temperature fell below the electron mass, electrons and positrons annihilated
and handed their entropy to the photons, not to the already-decoupled
neutrinos. Entropy conservation (g_s = 2 + 7/8 × 4 = 11/2 before, 2 after)
says the photons came out warmer by (11/4)^(1/3), so today

    T_ν = (4/11)^(1/3) T_γ = 1.945 K.

This module turns that into numbers (densities, mean momenta, speed
distributions, the redshift at which each mass became non-relativistic, the
free-streaming scale) and into the expected rate of a PTOLEMY-style tritium
capture experiment. All quantities use T_CMB = 2.7255 K (Fixsen 2009) unless
stated otherwise. The non-instantaneous-decoupling corrections come from
``decoupling.py``.
"""

from __future__ import annotations

import numpy as np
from scipy import integrate, special

# ---------------------------------------------------------------------------
# Constants (CODATA 2018 / SI 2019 exact values where defined)
# ---------------------------------------------------------------------------

K_B_EV = 8.617333262e-5          # Boltzmann constant, eV/K (exact)
HBARC_EV_CM = 1.973269804e-5     # hbar c, eV cm
C_KM_S = 299_792.458             # speed of light, km/s (exact)
C_CM_S = 2.99792458e10
ZETA3 = float(special.zeta(3.0))  # Riemann zeta(3) = 1.2020569...
G_SI = 6.67430e-11               # m^3 kg^-1 s^-2
MPC_M = 3.0856775814913673e22    # m
EV_J = 1.602176634e-19           # J (exact)
AVOGADRO = 6.02214076e23         # exact
YEAR_S = 365.25 * 86400.0
#: Mean momentum of a relativistic Fermi–Dirac gas, in units of k T:
#: <p> = 7 π^4 / (180 ζ(3)) k T = 3.1514 k T.
MEAN_P_FD = 7.0 * np.pi**4 / (180.0 * ZETA3)
#: Atomic mass of tritium in daltons (AME2020: 3.01604928 u).
TRITIUM_U = 3.01604928


def critical_density_ev_cm3_per_h2() -> float:
    """3 H0^2 / (8 pi G) c^2 for H0 = 100 km/s/Mpc, in eV per cm^3 (≈ 10,536.7)."""
    h0 = 100.0 * 1e3 / MPC_M  # s^-1
    rho_kg_m3 = 3.0 * h0**2 / (8.0 * np.pi * G_SI)
    return rho_kg_m3 * (C_CM_S / 100.0) ** 2 / EV_J / 1e6


RHO_CRIT_H2 = critical_density_ev_cm3_per_h2()


# ---------------------------------------------------------------------------
# Photons and the entropy argument
# ---------------------------------------------------------------------------


def kT_ev(T_K: float) -> float:
    return K_B_EV * T_K


def photon_number_density(T_K: float) -> float:
    """n_γ = 2 ζ(3)/π² (kT/ħc)³, per cm³ (410.7 at 2.7255 K)."""
    return 2.0 * ZETA3 / np.pi**2 * (kT_ev(T_K) / HBARC_EV_CM) ** 3


def photon_energy_density(T_K: float) -> float:
    """ρ_γ = π²/15 (kT)⁴/(ħc)³, in eV per cm³."""
    return np.pi**2 / 15.0 * kT_ev(T_K) ** 4 / HBARC_EV_CM**3


def entropy_ratio_tnu_tgamma() -> float:
    """(4/11)^(1/3): what e± annihilation leaves between the neutrino and photon temperatures."""
    return (4.0 / 11.0) ** (1.0 / 3.0)


def neutrino_temperature_K(T_cmb_K: float, neff: float = 3.0) -> float:
    """Relic neutrino temperature today.

    With ``neff = 3`` this is the instantaneous-decoupling value (4/11)^(1/3) T_γ.
    A different N_eff is interpreted as the same three flavours sharing a
    slightly higher (thermal-equivalent) temperature: ρ_ν ∝ N_eff T⁴ at fixed
    T_γ, so T_ν ∝ (N_eff / 3)^(1/4).
    """
    return entropy_ratio_tnu_tgamma() * T_cmb_K * (neff / 3.0) ** 0.25


def neutrino_number_density(T_nu_K: float) -> float:
    """Number density of ν + ν̄ of one flavour (two helicity states), per cm³.

    n = (3/4) × (ζ(3)/π²) × g × T³ with g = 2 (one ν and one ν̄ state);
    = (3/11) n_γ when T_ν = (4/11)^(1/3) T_γ, i.e. 112 per cm³.
    """
    return 1.5 * ZETA3 / np.pi**2 * (kT_ev(T_nu_K) / HBARC_EV_CM) ** 3


def massless_energy_density(T_nu_K: float) -> float:
    """ρ of one flavour (ν + ν̄) while relativistic: 7/8 × 2 × π²/30 (kT)⁴, eV/cm³."""
    return 7.0 / 8.0 * np.pi**2 / 15.0 * kT_ev(T_nu_K) ** 4 / HBARC_EV_CM**3


# ---------------------------------------------------------------------------
# Momenta, speeds and densities of massive relic neutrinos
# ---------------------------------------------------------------------------


def fd_pdf_momentum(x: np.ndarray) -> np.ndarray:
    """Probability density of x = p/(kT_ν) for a relativistic Fermi–Dirac gas.

    f(x) ∝ x² / (e^x + 1), normalised by ∫ = (3/2) ζ(3).
    """
    x = np.asarray(x, float)
    return x**2 / (np.exp(x) + 1.0) / (1.5 * ZETA3)


def mean_momentum_ev(T_nu_K: float) -> float:
    """<p> c = 3.1514 k T_ν (5.28e-4 eV today)."""
    return MEAN_P_FD * kT_ev(T_nu_K)


def speed_from_momentum(p_ev: np.ndarray, m_ev: float) -> np.ndarray:
    """v / c = p / sqrt(p² + m²) (p and m in eV, i.e. pc and mc²)."""
    p = np.asarray(p_ev, float)
    return p / np.sqrt(p**2 + m_ev**2)


def speed_distribution(m_ev: float, T_nu_K: float, n: int = 4000, vmax_kms: float | None = None):
    """Distribution of speeds today, dN/dv (normalised), for one mass eigenstate.

    The occupation number is frozen in momentum, so dN/dv = dN/dp × dp/dv with
    p = m v / sqrt(1 − v²). Returns (v [km/s], pdf [per km/s]).
    """
    kT = kT_ev(T_nu_K)
    if vmax_kms is None:
        vmax_kms = min(C_KM_S * 0.999999, C_KM_S * speed_from_momentum(15 * kT, m_ev))
    beta = np.linspace(0, vmax_kms / C_KM_S, n)
    gamma = 1.0 / np.sqrt(1.0 - beta**2)
    p = m_ev * beta * gamma
    dp_dbeta = m_ev * gamma**3
    pdf_beta = fd_pdf_momentum(p / kT) / kT * dp_dbeta
    return beta * C_KM_S, pdf_beta / C_KM_S


def speed_quantiles(m_ev: float, T_nu_K: float, q=(0.16, 0.5, 0.84)) -> dict:
    """Mean speed and speed quantiles today (km/s), computed exactly in momentum space."""
    kT = kT_ev(T_nu_K)
    x = np.linspace(0, 40, 400_001)
    cdf = integrate.cumulative_trapezoid(fd_pdf_momentum(x), x, initial=0.0)
    cdf /= cdf[-1]
    xq = np.interp(q, cdf, x)
    v = C_KM_S * speed_from_momentum(x * kT, m_ev)
    mean_v = float(integrate.trapezoid(v * fd_pdf_momentum(x), x))
    return {"mean_kms": mean_v,
            "quantiles_kms": {f"{qq:.2f}": float(C_KM_S * speed_from_momentum(xx * kT, m_ev))
                              for qq, xx in zip(q, xq)}}


def energy_density_massive(m_ev: float, T_nu_K: float) -> float:
    """Exact ρ (eV/cm³) of one flavour (ν + ν̄) of mass m, Fermi–Dirac in momentum."""
    kT = kT_ev(T_nu_K)
    pref = 2.0 / (2.0 * np.pi**2) * kT**3 / HBARC_EV_CM**3  # g = 2, per unit x³
    val, _ = integrate.quad(lambda x: x**2 * np.sqrt((x * kT) ** 2 + m_ev**2) / (np.exp(x) + 1.0),
                            0, 60, epsabs=0, epsrel=1e-12, limit=200)
    return pref * val


def omega_nu_h2(masses_ev, T_nu_K: float) -> float:
    """Ω_ν h² of a set of mass eigenstates today (exact Fermi–Dirac integral)."""
    return sum(energy_density_massive(m, T_nu_K) for m in masses_ev) / RHO_CRIT_H2


def mass_to_density_factor(T_nu_K: float) -> float:
    """The "93 eV": Σm / Ω_ν h² for fully non-relativistic neutrinos, eV."""
    return RHO_CRIT_H2 / neutrino_number_density(T_nu_K)


def z_nonrelativistic(m_ev: float, T_nu_K: float) -> float:
    """Redshift at which the mean momentum equalled m c (one common definition)."""
    return m_ev / mean_momentum_ev(T_nu_K) - 1.0


def free_streaming_k(m_ev: float, T_nu_K: float, omega_m: float, z: float) -> float:
    """Free-streaming wavenumber in h/Mpc during matter domination (non-relativistic ν).

    k_fs = sqrt(3/2) a H Ω_m(z)^(1/2) / v_th with v_th = <p>/m: the scale
    over which neutrino thermal motion beats gravitational infall (Lesgourgues &
    Pastor 2006, Sect. 4.5.2). In h/Mpc it is sqrt(3/2) × 100 × sqrt(Ω_m (1+z)) / v_th[km/s].
    """
    v_th = C_KM_S * mean_momentum_ev(T_nu_K) * (1.0 + z) / m_ev
    return np.sqrt(1.5) * 100.0 * np.sqrt(omega_m * (1.0 + z)) / v_th


def k_nonrelativistic(m_ev: float, T_nu_K: float, omega_m: float) -> float:
    """k_fs at the moment of the non-relativistic transition: the minimum k_fs,
    ≈ 0.018 Ω_m^(1/2) (m/1 eV)^(1/2) h/Mpc (Lesgourgues & Pastor 2006, eq. 97; m is
    the mass of one state). Scales smaller than 1/k_nr are suppressed. (Hu,
    Eisenstein & Tegmark 1998 quote 0.026 for a differently defined scale.)"""
    return free_streaming_k(m_ev, T_nu_K, omega_m, z_nonrelativistic(m_ev, T_nu_K))


# ---------------------------------------------------------------------------
# Direct detection: neutrino capture on tritium (PTOLEMY)
# ---------------------------------------------------------------------------


def tritium_capture_rate_per_year(sigma_v_cm2: float, n_per_state_cm3: float, mass_g: float,
                                  *, majorana: bool, clustering: float = 1.0) -> float:
    """Captures per year of relic ν_e on a tritium target, ν_e + ³H → ³He + e⁻.

    For non-relativistic neutrinos σ ∝ 1/v, so σv is a constant (``sigma_v_cm2``,
    in units of cm² × c). Every mass eigenstate carries the full electron
    flavour summed over i (Σ|U_ei|² = 1), and only left-helicity neutrinos are
    captured. Helicity is conserved as the neutrinos free-stream, so of the 112
    per cm³ of each mass eigenstate, 56 are left-helical neutrinos and 56
    right-helical antineutrinos. For Dirac neutrinos only the 56 can be
    captured; for Majorana neutrinos the right-helical "antineutrinos" are the
    same particles and, now that they are slow, are captured too, doubling the
    rate (Long, Lunardini & Sabancilar 2014). ``n_per_state_cm3`` is the density
    of one helicity state of one mass eigenstate (56 per cm³).
    """
    n_tritium = mass_g / TRITIUM_U * AVOGADRO
    n_captured = n_per_state_cm3 * (2.0 if majorana else 1.0) * clustering
    return sigma_v_cm2 * C_CM_S * n_captured * n_tritium * YEAR_S


def relic_summary(T_cmb_K: float, neff: float, masses_ev: dict, omega_m: float,
                  sigma_v_cm2: float, ptolemy_mass_g: float = 100.0,
                  heating: tuple[float, float] | None = None) -> dict:
    """Every relic-neutrino number the README quotes, for one T_CMB and N_eff.

    ``heating`` = (T_νe/T_γ, T_νμ/T_γ) relative to (4/11)^(1/3), from our
    decoupling calculation: the neutrinos' small share of the e± heat raises
    their temperatures, and so their numbers, by a few tenths of a percent.
    """
    t_inst = neutrino_temperature_K(T_cmb_K, 3.0)
    t_eff = neutrino_temperature_K(T_cmb_K, neff)
    n_gamma = photon_number_density(T_cmb_K)
    n_flav = neutrino_number_density(t_inst)
    out = {
        "T_cmb_K": T_cmb_K,
        "n_gamma_cm3": n_gamma,
        "rho_gamma_eV_cm3": photon_energy_density(T_cmb_K),
        "T_nu_K": t_inst,
        "T_nu_eV": kT_ev(t_inst),
        "T_nu_K_neff": t_eff,
        "n_nu_per_flavour_cm3": n_flav,
        "n_nu_total_cm3": 3.0 * n_flav,
        "n_nu_total_cm3_neff_scaled": 3.0 * neutrino_number_density(t_eff),
        "mean_p_eV": mean_momentum_ev(t_inst),
        "mass_to_density_eV": mass_to_density_factor(t_inst),
        "mass_to_density_eV_neff": mass_to_density_factor(t_eff),
        "rho_crit_h2_eV_cm3": RHO_CRIT_H2,
        "states": {},
    }
    if heating is not None:
        re, rm = heating
        n_dec = n_flav * (re**3 + 2 * rm**3)
        out["heating_ratio_e"], out["heating_ratio_mu"] = re, rm
        out["T_nue_K_dec"], out["T_numu_K_dec"] = t_inst * re, t_inst * rm
        out["n_nu_total_cm3_dec"] = n_dec
        out["mass_to_density_eV_dec"] = RHO_CRIT_H2 / (n_dec / 3.0)
    for label, m in masses_ev.items():
        s = {"m_eV": m}
        if m > 0:
            s.update(speed_quantiles(m, t_inst))
            s["z_nr"] = z_nonrelativistic(m, t_inst)
            s["k_nr_h_Mpc"] = k_nonrelativistic(m, t_inst, omega_m)
            s["omega_nu_h2"] = omega_nu_h2([m], t_inst)
        else:
            s.update({"mean_kms": C_KM_S, "z_nr": None})
        out["states"][label] = s
    n_state = n_flav / 2.0
    out["ptolemy"] = {
        "target_g": ptolemy_mass_g,
        "sigma_v_cm2": sigma_v_cm2,
        "n_per_helicity_state_cm3": n_state,
        "rate_dirac_per_yr": tritium_capture_rate_per_year(sigma_v_cm2, n_state, ptolemy_mass_g,
                                                            majorana=False),
        "rate_majorana_per_yr": tritium_capture_rate_per_year(sigma_v_cm2, n_state, ptolemy_mass_g,
                                                               majorana=True),
    }
    return out
