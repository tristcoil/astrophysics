"""The physics of a radioactively powered supernova.

Explosive burning makes 56Ni. It decays (half-life 6.1 d) to 56Co, which decays
(half-life 77.2 d) to stable 56Fe. Each decay releases gamma rays (and, for 19%
of cobalt decays, a positron). In the dense young ejecta the gamma rays are
trapped and thermalised, and the energy comes out as optical light. As the
ejecta thin out, more and more gamma rays escape, while the positrons, being
charged, keep depositing their energy locally for years.

Everything here is a function of time since explosion ``t`` in days and is
written per solar mass of 56Ni, so that ``M_Ni * power(t)`` is a luminosity in
erg/s.

References: Arnett (1982, ApJ 253, 785) for the diffusion model and "Arnett's
rule"; Nadyozhin (1994, ApJS 92, 527) for the energies per decay; Jeffery
(1999, astro-ph/9907015), Wygoda et al. (2019, MNRAS 484, 3941) and Sharon &
Kushnir (2020, MNRAS 496, 4517) for the gamma-ray escape time T0; Katz, Kushnir
& Dong (2013, arXiv:1301.6766) for the integral test.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize

import config as C

LN2 = np.log(2.0)

#: Decay rates in 1/day, and mean lifetimes tau = 1/lambda (8.76 d and 111.4 d).
LAMBDA_NI = LN2 / C.HALF_LIFE_NI56_D
LAMBDA_CO = LN2 / C.HALF_LIFE_CO56_D
TAU_NI = 1.0 / LAMBDA_NI
TAU_CO = 1.0 / LAMBDA_CO

#: Number of 56Ni nuclei in one solar mass.
N_PER_MSUN = C.MSUN_G / (C.MASS_NI56_AMU * C.AMU_G)

#: Cobalt's decay rate in magnitudes per 100 days: 2.5 log10(e) * 100 / tau_Co.
#: This is the slope a fully trapped radioactive tail must follow.
CO_SLOPE_MAG_PER_100D = 2.5 / np.log(10.0) * 100.0 / TAU_CO


def nuclei(t_days, m_ni: float = 1.0) -> dict:
    """Numbers of 56Ni, 56Co and 56Fe nuclei at time t for ``m_ni`` solar masses."""
    t = np.asarray(t_days, dtype=float)
    n0 = m_ni * N_PER_MSUN
    ni = n0 * np.exp(-LAMBDA_NI * t)
    co = n0 * LAMBDA_NI / (LAMBDA_CO - LAMBDA_NI) * (np.exp(-LAMBDA_NI * t) - np.exp(-LAMBDA_CO * t))
    return {"Ni": ni, "Co": co, "Fe": n0 - ni - co}


def power(t_days) -> dict:
    """Radioactive power per solar mass of initial 56Ni, in erg/s, by channel.

    Returns ``ni_gamma`` (gamma rays from 56Ni), ``co_gamma`` (gamma rays and
    annihilation photons from 56Co), ``co_positron`` (positron kinetic energy)
    and their ``total``. Neutrino energy is excluded: it escapes.
    """
    n = nuclei(t_days, 1.0)
    per_day = 1.0 / C.DAY_S  # decay rates are per day; luminosities per second
    ni_gamma = LAMBDA_NI * per_day * n["Ni"] * C.E_GAMMA_NI_MEV * C.MEV_ERG
    co_rate = LAMBDA_CO * per_day * n["Co"]
    co_gamma = co_rate * C.E_GAMMA_CO_MEV * C.MEV_ERG
    co_positron = co_rate * C.E_POSITRON_CO_MEV * C.MEV_ERG
    return {
        "ni_gamma": ni_gamma,
        "co_gamma": co_gamma,
        "co_positron": co_positron,
        "total": ni_gamma + co_gamma + co_positron,
    }


def two_exponential_coefficients() -> tuple[float, float]:
    """Coefficients (A, B) of the textbook form q(t) = A e^(-t/tau_Ni) + B e^(-t/tau_Co).

    Rearranging the Bateman solution gives B = N0 E_Co lambda_Co lambda_Ni /
    (lambda_Ni - lambda_Co) and A = N0 E_Ni lambda_Ni - B. Nadyozhin (1994)
    quotes A = 6.45e43 and B = 1.45e43 erg/s per solar mass.
    """
    per_s = 1.0 / C.DAY_S
    e_ni = C.E_GAMMA_NI_MEV * C.MEV_ERG
    e_co = (C.E_GAMMA_CO_MEV + C.E_POSITRON_CO_MEV) * C.MEV_ERG
    b = N_PER_MSUN * e_co * LAMBDA_CO * per_s * LAMBDA_NI / (LAMBDA_NI - LAMBDA_CO)
    a = N_PER_MSUN * e_ni * LAMBDA_NI * per_s - b
    return a, b


def gamma_trapping(t_days, t0_days: float | None):
    """Fraction of gamma-ray energy deposited in the ejecta, 1 - exp(-(T0/t)^2).

    The gamma-ray optical depth of homologously expanding ejecta falls as
    t^-2 (column density ~ M / (v t)^2); T0 is the time at which it reaches one.
    ``t0_days=None`` means complete trapping.
    """
    t = np.asarray(t_days, dtype=float)
    if t0_days is None or not np.isfinite(t0_days):
        return np.ones_like(t)
    with np.errstate(divide="ignore", over="ignore"):
        x = np.where(t > 0, (t0_days / np.where(t > 0, t, 1.0)) ** 2, np.inf)
    return -np.expm1(-x)


def deposition(t_days, t0_days: float | None = None, positrons_trapped: bool = True):
    """Energy deposited per solar mass of 56Ni (erg/s): leaky gammas + local positrons."""
    p = power(t_days)
    f = gamma_trapping(t_days, t0_days)
    pos = p["co_positron"] if positrons_trapped else 0.0
    return (p["ni_gamma"] + p["co_gamma"]) * f + pos


def arnett_light_curve(t_days, m_ni: float, tau_m_days: float, t0_days: float | None = None,
                       dt: float = 0.05):
    """Arnett (1982) light curve for M_Ni solar masses, in erg/s.

    L(t) = exp(-x^2) * Integral_0^t (2 t'/tau_m^2) exp((t'/tau_m)^2) D(t') dt',
    with x = t / tau_m and D the deposited radioactive power. tau_m is the
    effective diffusion time, (2 kappa M_ej / (beta c v))^(1/2) with beta ~ 13.8.
    The integral is advanced step by step so the large exponentials never form.
    """
    t_out = np.atleast_1d(np.asarray(t_days, dtype=float))
    t_max = max(float(np.nanmax(t_out)), 1.0)
    grid = np.arange(0.0, t_max + dt, dt)
    dep = m_ni * deposition(grid, t0_days)
    dep[0] = m_ni * deposition(1e-6, t0_days)
    tm2 = tau_m_days**2
    lum = np.zeros_like(grid)
    for i in range(len(grid) - 1):
        t1, t2 = grid[i], grid[i + 1]
        damp = np.exp(-(t2 * t2 - t1 * t1) / tm2)
        inc = 0.5 * (t2 - t1) * (2 * t1 / tm2 * damp * dep[i] + 2 * t2 / tm2 * dep[i + 1])
        lum[i + 1] = lum[i] * damp + inc
    return np.interp(t_out, grid, lum)


def arnett_rule_mass(l_peak: float, t_rise_days: float, alpha: float = 1.0) -> float:
    """Nickel mass (Msun) from Arnett's rule, L_peak = alpha * q(t_rise) * M_Ni.

    Near peak the gamma rays are fully trapped, so the total radioactive power
    at the time of peak sets the peak luminosity, to within the factor alpha
    (models give 0.8-1.2).
    """
    return float(l_peak / (alpha * power(t_rise_days)["total"]))


def katz_ratio(t_days, lum, t0_days: float | None) -> float:
    """Nickel mass from the Katz integral: Int t L dt / Int t D_1(t) dt.

    For homologous ejecta with negligible initial thermal energy, the
    time-weighted integral of the light curve equals that of the deposited
    radioactive power (Katz et al. 2013). It is independent of the diffusion
    physics that Arnett's rule approximates. ``t_days`` must start near
    explosion.
    """
    t = np.asarray(t_days, dtype=float)
    lum = np.asarray(lum, dtype=float)
    fine = np.linspace(0.0, t.max(), 4000)
    num = np.trapezoid(t * lum, t)
    # The data start a little after explosion. The missing sliver is tiny and
    # is added assuming the early "fireball" rise L ~ t^2: Int t L dt = L0 t0^2 / 4.
    if t[0] > 0:
        num += lum[0] * t[0] ** 2 / 4.0
    den = np.trapezoid(fine * deposition(fine, t0_days), fine)
    return float(num / den)


def fit_tail(t_days, lum, lum_err=None, fit_t0: bool = True, t0_guess: float = 40.0,
             full_trapping_t0: float | None = None) -> dict:
    """Fit L(t) = M_Ni * D(t; T0) to a radioactive tail (log-space least squares).

    With ``fit_t0=False`` the gamma rays are taken as fully trapped (or trapped
    with a fixed ``full_trapping_t0``). Returns the best-fit nickel mass, T0,
    their 1-sigma errors, and the residual scatter in dex.
    """
    t = np.asarray(t_days, float)
    y = np.log10(np.asarray(lum, float))
    if lum_err is None:
        sig = np.full_like(y, 0.02)
    else:
        sig = np.asarray(lum_err, float) / (np.asarray(lum, float) * np.log(10))
        sig = np.clip(sig, 0.005, None)

    if fit_t0:
        def model(tt, log_m, log_t0):
            return log_m + np.log10(deposition(tt, 10**log_t0))
        p0 = [np.log10(np.median(10**y / deposition(t, t0_guess))), np.log10(t0_guess)]
        popt, pcov = optimize.curve_fit(model, t, y, p0=p0, sigma=sig, absolute_sigma=False)
        perr = np.sqrt(np.diag(pcov))
        m, t0 = 10 ** popt[0], 10 ** popt[1]
        res = y - model(t, *popt)
        return {"m_ni": m, "m_ni_err": m * np.log(10) * perr[0], "t0": t0,
                "t0_err": t0 * np.log(10) * perr[1], "rms_dex": float(np.std(res)),
                "n": int(len(t))}

    def model1(tt, log_m):
        return log_m + np.log10(deposition(tt, full_trapping_t0))
    p0 = [np.log10(np.median(10**y / deposition(t, full_trapping_t0)))]
    popt, pcov = optimize.curve_fit(model1, t, y, p0=p0, sigma=sig, absolute_sigma=False)
    m = 10 ** popt[0]
    res = y - model1(t, *popt)
    return {"m_ni": m, "m_ni_err": m * np.log(10) * np.sqrt(pcov[0, 0]),
            "t0": full_trapping_t0, "t0_err": None, "rms_dex": float(np.std(res)),
            "n": int(len(t))}


def fit_arnett(t_days, lum, lum_err=None, t_max_fit: float = 150.0) -> dict:
    """Fit Arnett's model with gamma-ray leakage: M_Ni, tau_m and T0 together.

    Fitting the peak and the early tail at once breaks the degeneracy a tail
    alone suffers from (more nickel with faster gamma-ray escape looks much
    like less nickel with slower escape): the peak fixes the nickel, the
    tail's steepening fixes T0. Errors are the 1-sigma statistical ones,
    rescaled by the fit's residual scatter.
    """
    t = np.asarray(t_days, float)
    y = np.asarray(lum, float)
    sel = (t <= t_max_fit) & (t > 0)
    t, y = t[sel], y[sel]
    sig = np.full_like(y, 0.03) if lum_err is None else (np.asarray(lum_err, float)[sel] / (y * np.log(10)))
    sig = np.clip(sig, 0.01, None)

    def model(tt, m, tau, t0):
        return np.log10(np.clip(arnett_light_curve(tt, m, tau, t0, dt=0.1), 1e30, None))

    popt, pcov = optimize.curve_fit(model, t, np.log10(y), p0=[0.5, 15.0, 35.0], sigma=sig,
                                    bounds=([0.01, 2.0, 5.0], [3.0, 60.0, 300.0]))
    perr = np.sqrt(np.diag(pcov))
    res = np.log10(y) - model(t, *popt)
    return {"m_ni": float(popt[0]), "m_ni_err": float(perr[0]), "tau_m": float(popt[1]),
            "tau_m_err": float(perr[1]), "t0": float(popt[2]), "t0_err": float(perr[2]),
            "rms_dex": float(np.std(res)), "t_max_fit": t_max_fit, "n": int(len(t))}


def fit_arnett_fixed_t0(t_days, lum, lum_err, t0_days: float) -> dict:
    """Arnett's model with T0 held fixed: fit M_Ni and tau_m only."""
    t = np.asarray(t_days, float)
    y = np.asarray(lum, float)
    sig = np.clip(np.asarray(lum_err, float) / (y * np.log(10)), 0.01, None)

    def model(tt, m, tau):
        return np.log10(np.clip(arnett_light_curve(tt, m, tau, t0_days, dt=0.1), 1e30, None))

    popt, pcov = optimize.curve_fit(model, t, np.log10(y), p0=[0.5, 15.0], sigma=sig,
                                    bounds=([0.01, 2.0], [3.0, 60.0]))
    perr = np.sqrt(np.diag(pcov))
    return {"m_ni": float(popt[0]), "m_ni_err": float(perr[0]), "tau_m": float(popt[1]),
            "tau_m_err": float(perr[1]), "t_max_fit": float(t.max())}

