"""Binary pulsars: periastron advance as a test, not an input.

A pulsar in a binary is a clock on an orbit. Timing it gives the Keplerian
orbit plus "post-Keplerian" effects, each of which general relativity predicts
from the two masses alone (Damour & Deruelle 1986; Taylor & Weisberg 1989):

* periastron advance  omega-dot = 3 (Pb/2pi)^(-5/3) (T M)^(2/3) / (1 - e^2)
* time dilation + gravitational redshift
                      gamma = e (Pb/2pi)^(1/3) T^(2/3) M^(-4/3) m_c (m_p + 2 m_c)
* orbital decay       Pb-dot = -(192 pi/5) (Pb/2pi)^(-5/3) f(e) T^(5/3) m_p m_c M^(-1/3)
* Shapiro delay shape s = sin i, with the mass function giving m_c sin i
* mass ratio          R = m_p/m_c = x_c/x_p (theory-independent, double pulsar only)

with T = GM_sun/c^3 = 4.925 microseconds and masses in solar units.

Published masses use omega-dot itself, so a plot of omega-dot against GR is
consistent by construction. Here the masses come from two *other* effects,
and omega-dot is then predicted and compared with the measured one.
Uncertainties by Monte Carlo over the measured inputs.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq, fsolve

from config import C_SI, GM_SUN_SI, JULIAN_YEAR_S, PUBLISHED, T_SUN_S

DAY_S = 86400.0


def omdot_gr(pb_s: float, e: float, m_total: float) -> float:
    """Periastron advance, rad/s."""
    return 3.0 * (pb_s / (2 * np.pi)) ** (-5 / 3) * (T_SUN_S * m_total) ** (2 / 3) / (1 - e**2)


def gamma_gr(pb_s, e, m_p, m_c):
    M = m_p + m_c
    return e * (pb_s / (2 * np.pi)) ** (1 / 3) * T_SUN_S ** (2 / 3) * M ** (-4 / 3) * m_c * (m_p + 2 * m_c)


def pbdot_gr(pb_s, e, m_p, m_c):
    M = m_p + m_c
    fe = (1 + 73 / 24 * e**2 + 37 / 96 * e**4) * (1 - e**2) ** (-3.5)
    return -(192 * np.pi / 5) * (pb_s / (2 * np.pi)) ** (-5 / 3) * fe * T_SUN_S ** (5 / 3) * m_p * m_c * M ** (-1 / 3)


def mass_function(pb_s, x_s):
    """(m_c sin i)^3 / M^2 in solar masses, from the pulsar's projected orbit."""
    return (2 * np.pi / pb_s) ** 2 * x_s**3 / T_SUN_S


def deg_per_year(rad_per_s):
    return np.degrees(rad_per_s) * JULIAN_YEAR_S


def masses_gamma_pbdot(pb_s, e, gamma, pbdot):
    """Solve gamma and Pb-dot for (m_p, m_c)."""
    def eqs(m):
        mp, mc = m
        return [gamma_gr(pb_s, e, mp, mc) / gamma - 1, pbdot_gr(pb_s, e, mp, mc) / pbdot - 1]
    sol, info, ier, msg = fsolve(eqs, [1.4, 1.4], full_output=True, xtol=1e-13)
    if ier != 1:
        raise RuntimeError(msg)
    return sol


def masses_ratio_gamma(pb_s, e, gamma, R):
    """Double pulsar: m_A = R m_B, and gamma (of A) fixes the scale analytically."""
    k = e * (pb_s / (2 * np.pi)) ** (1 / 3) * T_SUN_S ** (2 / 3) * (R + 2) * (1 + R) ** (-4 / 3)
    m_b = (gamma / k) ** 1.5
    return R * m_b, m_b


def masses_ratio_s(pb_s, x_a, s, R):
    """Double pulsar: mass function of A, sin i = s and R give m_B = f (1+R)^2 / s^3."""
    m_b = mass_function(pb_s, x_a) * (1 + R) ** 2 / s**3
    return R * m_b, m_b


def masses_s_pbdot(pb_s, e, x_a, s, pbdot):
    """Double pulsar: Shapiro shape plus orbital decay (the paper's own route)."""
    f = mass_function(pb_s, x_a)

    def resid(m_b):
        # m_b^3 s^3 / M^2 = f  ->  M = sqrt(m_b^3 s^3 / f)
        M = np.sqrt(m_b**3 * s**3 / f)
        return pbdot_gr(pb_s, e, M - m_b, m_b) / pbdot - 1
    m_b = brentq(resid, 0.5, 3.0, xtol=1e-14)
    M = np.sqrt(m_b**3 * s**3 / f)
    return M - m_b, m_b


def _mc(fun, draws: int = 20000, seed: int = 7) -> dict:
    rng = np.random.default_rng(seed)
    vals = np.array([fun(rng) for _ in range(draws)])
    ratio = vals[:, 0]
    return {"ratio": float(np.median(ratio)), "ratio_err": float(np.std(ratio)),
            "m_p": float(np.median(vals[:, 1])), "m_p_err": float(np.std(vals[:, 1])),
            "m_c": float(np.median(vals[:, 2])), "m_c_err": float(np.std(vals[:, 2])),
            "omdot_pred": float(np.median(vals[:, 3])), "omdot_pred_err": float(np.std(vals[:, 3]))}


def hulse_taylor() -> dict:
    p = PUBLISHED["b1913"]
    pb = p["pb_days"] * DAY_S

    def draw(rng):
        g = rng.normal(p["gamma_s"], p["gamma_err"])
        pd = rng.normal(p["pbdot_intrinsic"], p["pbdot_intrinsic_err"])
        mp, mc = masses_gamma_pbdot(pb, p["e"], g, pd)
        pred = deg_per_year(omdot_gr(pb, p["e"], mp + mc))
        return [p["omdot_deg_yr"] / pred, mp, mc, pred]
    out = _mc(draw, draws=4000)
    mp, mc = p["m_p"], p["m_c"]
    out.update({"route": "masses from gamma and the intrinsic Pb-dot",
                "omdot_obs": p["omdot_deg_yr"], "pb_days": p["pb_days"], "e": p["e"],
                "m_total_published": mp + mc,
                "compactness": compactness(pb, mp + mc, p["e"]),
                "per_orbit_rad_obs": np.radians(p["omdot_deg_yr"]) * pb / JULIAN_YEAR_S,
                "pbdot_ratio_published": p["pbdot_ratio"], "pbdot_ratio_err": p["pbdot_ratio_err"]})
    return out


def double_pulsar() -> dict:
    p = PUBLISHED["j0737"]
    pb = p["pb_days"] * DAY_S
    s_err = 0.5 * (p["s_err_plus"] + p["s_err_minus"])
    res = {}

    def draw_rg(rng):
        R = rng.normal(p["mass_ratio"], p["mass_ratio_err"])
        g = rng.normal(p["gamma_s"], p["gamma_err"])
        ma, mb = masses_ratio_gamma(pb, p["e"], g, R)
        pred = deg_per_year(omdot_gr(pb, p["e"], ma + mb))
        return [p["omdot_deg_yr"] / pred, ma, mb, pred]

    def draw_rs(rng):
        R = rng.normal(p["mass_ratio"], p["mass_ratio_err"])
        s = min(rng.normal(p["s"], s_err), 1.0)
        ma, mb = masses_ratio_s(pb, p["x_a_s"], s, R)
        pred = deg_per_year(omdot_gr(pb, p["e"], ma + mb))
        return [p["omdot_deg_yr"] / pred, ma, mb, pred]

    def draw_spb(rng):
        s = min(rng.normal(p["s"], s_err), 1.0)
        pd = rng.normal(p["pbdot_gw"], p["pbdot_gw_err"])
        ma, mb = masses_s_pbdot(pb, p["e"], p["x_a_s"], s, pd)
        pred = deg_per_year(omdot_gr(pb, p["e"], ma + mb))
        return [p["omdot_deg_yr"] / pred, ma, mb, pred]

    res["ratio_gamma"] = {"route": "masses from the mass ratio R and gamma", **_mc(draw_rg)}
    res["ratio_s"] = {"route": "masses from R and the Shapiro shape s", **_mc(draw_rs)}
    res["s_pbdot"] = {"route": "masses from s and the GW part of Pb-dot", **_mc(draw_spb, draws=4000)}
    M = p["m_a"] + p["m_b"]
    res.update({"omdot_obs": p["omdot_deg_yr"], "pb_days": p["pb_days"], "e": p["e"],
                "m_total_published": M, "compactness": compactness(pb, M, p["e"]),
                "per_orbit_rad_obs": np.radians(p["omdot_deg_yr"]) * pb / JULIAN_YEAR_S,
                "k_test_published": p["k_test"], "k_test_err": p["k_test_err"]})
    return res


def compactness(pb_s: float, m_total: float, e: float) -> float:
    gm = GM_SUN_SI * m_total
    a = (gm * pb_s**2 / (4 * np.pi**2)) ** (1 / 3)
    return gm / (C_SI**2 * a * (1 - e**2))


def analyse() -> dict:
    return {"b1913": hulse_taylor(), "j0737": double_pulsar()}
