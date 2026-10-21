"""Why light redshifts in an expanding universe: no time-translation symmetry.

A flat FLRW universe, ds^2 = -dt^2 + a(t)^2 (dx^2 + dy^2 + dz^2), has six
Killing vectors (three translations, three rotations of space) but no
timelike one, because the metric changes with time. So nothing makes a
photon's energy constant. What it does have is a **conformal** Killing vector,
a(t) d_t, with L g = 2 (da/dt) g. For light (and only light, since only for
null geodesics does a conformal factor drop out), g(a d_t, k) is conserved:

    a(t) x (photon energy measured by a comoving observer) = constant,

so wavelengths stretch with the universe: 1 + z = a_now / a_then. Every photon
loses energy by the same factor, so a thermal (Planck) spectrum stays thermal
with T = T_0 (1 + z), which D7 measured at 46 redshifts.

For a massive particle the spatial Killing vectors conserve a^2 dx/dtau, so its
momentum relative to the comoving observers falls as 1/a: peculiar velocities
die away as the universe expands.

Here: D7's fit of T(z) = T0 (1 + z)^(1 - beta) redone from D7's own table, and
geodesics integrated numerically in Planck 2018 LCDM to show the two
conservation laws (and that a E is NOT conserved for massive particles).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import integrate

import config as C


def load_tcmb() -> pd.DataFrame:
    df = pd.read_csv(C.D7_TCMB_TABLE, comment="#")
    df["sets"] = df["sets"].str.split(";")
    return df


def tcmb_fit(df=None, t0=C.T_CMB) -> dict:
    """D7's statistic, recoded: split-normal chi^2 on the 'baseline' set, T0 fixed."""
    df = load_tcmb() if df is None else df
    d = df[df["sets"].apply(lambda s: "baseline" in s) & (df["kind"] != "upper")]
    z, T = d["z"].to_numpy(float), d["T"].to_numpy(float)
    ep, em, es = (d[c_].to_numpy(float) for c_ in ("err_plus", "err_minus", "err_sys"))

    def chi2(beta):
        model = t0 * (1 + z) ** (1 - beta)
        res = model - T
        err = np.hypot(np.where(res > 0, ep, em), es)
        return float(np.sum((res / err) ** 2))

    grid = np.linspace(-0.3, 0.3, 60001)
    c2 = np.array([chi2(b_) for b_ in grid])
    i = int(np.argmin(c2))
    inside = grid[c2 <= c2[i] + 1]
    return {"n": int(len(d)), "beta": float(grid[i]), "plus": float(inside.max() - grid[i]),
            "minus": float(grid[i] - inside.min()), "chi2_min": float(c2[i]),
            "chi2_expansion": chi2(0.0), "chi2_static": chi2(1.0), "dof": int(len(d) - 1),
            "points": {"z": z.tolist(), "T": T.tolist(), "ep": ep.tolist(), "em": em.tolist(), "es": es.tolist(),
                       "method": d["method"].tolist()},
            "upper_limits": df[df["sets"].apply(lambda s: "baseline" in s) & (df["kind"] == "upper")][["z", "T"]]
            .to_dict("list")}


def hubble_lcdm(a):
    """H(a) in 1/Gyr for Planck 2018 LCDM (radiation + matter + Lambda, flat)."""
    H0 = C.PLANCK18["H0"] * 1e3 / (1e6 * C.PC) * C.GYR
    Om, Or = C.PLANCK18["Omega_m"], C.PLANCK18["Omega_r"]
    return H0 * np.sqrt(Or / a**4 + Om / a**3 + (1 - Om - Or))


def age_gyr(a=1.0) -> float:
    val, _ = integrate.quad(lambda lna: 1 / hubble_lcdm(np.exp(lna)), -40, np.log(a), limit=500)
    return float(val)


def geodesics_demo() -> dict:
    """Integrate a photon and a massive particle through LCDM with the geodesic
    equations, which know nothing about Killing vectors:

        d2t/dl2 = -a a' (dx/dl)^2,   d2x/dl2 = -2 (a'/a) (dt/dl)(dx/dl),

    with the scale factor carried along in the same integration
    (da/dl = a H(a) dt/dl), so a and a' = a H stay exactly consistent. We then
    check the products a E (conformal Killing vector a d_t) and a^2 dx/dl
    (spatial Killing vector d_x) along the way."""

    def rhs(_, s):
        aa, x_, ut, ux = s
        Hh = hubble_lcdm(aa)
        return [aa * Hh * ut, ux, -aa**2 * Hh * ux**2, -2 * Hh * ut * ux]

    def stop(_, s):
        return s[0] - 1.0
    stop.terminal = True

    def run(a0, ut0, ux0, lam_max):
        sol = integrate.solve_ivp(rhs, (0, lam_max), [a0, 0.0, ut0, ux0], method="DOP853", rtol=1e-12,
                                  atol=1e-16, events=stop, dense_output=True)
        lam = np.linspace(0, sol.t[-1], 800)
        return sol.sol(lam)

    out = {}
    # a photon emitted at recombination (z = 1090), heading along +x: null means dt = a dx
    a0 = 1 / 1091
    aa, xx, ut, ux = run(a0, 1.0, 1.0 / a0, 1e7)
    out["photon"] = {"a": aa.tolist(), "E": ut.tolist(), "aE": (aa * ut).tolist(), "a2dx": (aa**2 * ux).tolist(),
                     "null_rel": float(np.max(np.abs((-ut**2 + aa**2 * ux**2) / ut**2))),
                     "aE_spread": float(np.ptp(aa * ut) / np.mean(aa * ut)),
                     "E_ratio_start_end": float(ut[0] / ut[-1]), "a_ratio": float(aa[-1] / aa[0])}
    # a galaxy given a peculiar speed of 0.6c at z = 10, coasting since
    a0, v0 = 1 / 11, 0.6
    g0 = 1 / np.sqrt(1 - v0**2)
    aa, xx, ut, ux = run(a0, g0, g0 * v0 / a0, 1e5)   # u^t = gamma, a u^x = gamma v
    p = aa * ux                                         # momentum per unit mass seen by comoving observers
    v = p / np.sqrt(1 + p**2)
    out["massive"] = {"a": aa.tolist(), "v": v.tolist(), "p": p.tolist(), "a2dx": (aa**2 * ux).tolist(),
                      "aE": (aa * ut).tolist(), "a2dx_spread": float(np.ptp(aa**2 * ux) / np.mean(aa**2 * ux)),
                      "aE_change": float(aa[-1] * ut[-1] / (aa[0] * ut[0])),
                      "v_today": float(v[-1]), "v_today_from_killing": float((g0 * v0 * a0) / np.sqrt(1 + (g0 * v0 * a0) ** 2)),
                      "norm_rel": float(np.max(np.abs(-ut**2 + aa**2 * ux**2 + 1)))}
    out["age_gyr"] = age_gyr()
    return out
