"""Spinning black holes: frame dragging, the ergosphere, and the plunge from the innermost orbit.

Kerr metric, equatorial plane, units G = c = M = 1, spin a = J c / (G M^2) in [0, 1).
Boyer-Lindquist (BL) functions: Delta = r^2 - 2r + a^2, horizons r_pm = 1 pm sqrt(1 - a^2),
equatorial ergosphere edge at r = 2 (where g_tt changes sign).

Equatorial timelike geodesics with energy E and angular momentum L per unit mass:

    r^2 dr/dtau  = -sqrt(R),    R = P^2 - Delta [r^2 + (L - aE)^2],   P = E (r^2 + a^2) - a L
    r^2 dphi/dtau = -(aE - L) + a P / Delta           (BL: diverges at the horizon)
    r^2 dt/dtau   = -a (aE - L) + (r^2 + a^2) P / Delta

To follow a body through the horizon we use ingoing Kerr coordinates,
dphi~ = dphi + a dr/Delta, dv = dt + (r^2 + a^2) dr/Delta, in which (for infall)
the 1/Delta terms combine into (P - sqrt R)/Delta = [r^2 + (L - aE)^2] / (P + sqrt R),
which is finite at the horizon. Plotted in Kerr-Schild Cartesian coordinates,
x + i y = (r + i a) exp(i phi~), the horizon is a circle of radius sqrt(r_+^2 + a^2).
"""

from __future__ import annotations

import math

import numpy as np
from scipy import integrate


def horizons(a):
    s = math.sqrt(1.0 - a * a)
    return 1.0 + s, 1.0 - s


def isco(a, prograde: bool = True):
    """Innermost stable circular orbit radius (Bardeen, Press & Teukolsky 1972)."""
    z1 = 1.0 + (1.0 - a * a) ** (1 / 3) * ((1 + a) ** (1 / 3) + (1 - a) ** (1 / 3))
    z2 = math.sqrt(3 * a * a + z1 * z1)
    sgn = -1.0 if prograde else 1.0
    return 3.0 + z2 + sgn * math.sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2))


def circular_EL(r, a, prograde: bool = True):
    """Energy and angular momentum per unit mass of a circular equatorial orbit (BPT 1972)."""
    s = 1.0 if prograde else -1.0
    sa = s * a
    den = r**0.75 * math.sqrt(r**1.5 - 3.0 * r**0.5 + 2.0 * sa)
    E = (r**1.5 - 2.0 * r**0.5 + sa) / den
    L = s * (r**2 - 2.0 * sa * r**0.5 + a * a) / den
    return E, L


def efficiency(a, prograde: bool = True):
    """Fraction of rest-mass energy released by gas spiralling to the ISCO: 1 - E_isco."""
    return 1.0 - circular_EL(isco(a, prograde), a, prograde)[0]


def frame_dragging_omega(r, a):
    """Angular velocity d phi/dt of zero-angular-momentum observers in the equatorial plane."""
    r = np.asarray(r, float)
    return 2.0 * a / (r**3 + a * a * r + 2.0 * a * a)


def allowed_omega(r, a):
    """Range (Omega_min, Omega_max) of d phi/dt for any timelike world line at fixed r (equator).

    Inside the ergosphere (r < 2) Omega_min > 0: everything must co-rotate.
    """
    r = np.asarray(r, float)
    gtt = -(1.0 - 2.0 / r)
    gtp = -2.0 * a / r
    gpp = r * r + a * a + 2.0 * a * a / r
    disc = np.sqrt(np.maximum(gtp * gtp - gtt * gpp, 0.0))
    return (-gtp - disc) / gpp, (-gtp + disc) / gpp


def _rhs(tau, y, a, E, L):
    """y = (r, phi~, v); infalling equatorial geodesic in ingoing Kerr coordinates."""
    r = y[0]
    delta = r * r - 2.0 * r + a * a
    P = E * (r * r + a * a) - a * L
    K = r * r + (L - a * E) ** 2
    R = max(P * P - delta * K, 0.0)
    sR = math.sqrt(R)
    reg = K / (P + sR)                     # = (P - sqrt R) / Delta, finite at the horizon
    drdt = -sR / (r * r)
    dphit = (-(a * E - L) + a * reg) / (r * r)
    dv = (-a * (a * E - L) + (r * r + a * a) * reg) / (r * r)
    return [drdt, dphit, dv]


def plunge(a, E, L, r_start, r_stop, phi0=0.0, max_tau=1e5, n=4000):
    """Integrate an infalling equatorial geodesic from r_start to r_stop.

    Returns dict with tau, r, phi~ (ingoing), v, Kerr-Schild x, y, and the proper and
    ingoing-time spent between r_start and the outer horizon.
    """
    def stop(tau, y, *args):
        return y[0] - r_stop
    stop.terminal = True
    # near a turning point (circular orbit) sqrt(R) ~ 0 and the start is slow; that is physical
    sol = integrate.solve_ivp(_rhs, (0.0, max_tau), [r_start, phi0, 0.0], args=(a, E, L), method="DOP853",
                              rtol=1e-10, atol=1e-12, events=stop, dense_output=True, max_step=0.5)
    tau = np.linspace(0.0, sol.t[-1], n)
    r, ph, v = sol.sol(tau)
    x = r * np.cos(ph) - a * np.sin(ph)
    yy = r * np.sin(ph) + a * np.cos(ph)
    rp, _ = horizons(a)
    cross = np.flatnonzero(r <= rp)
    tau_h = float(tau[cross[0]]) if cross.size else math.nan
    v_h = float(v[cross[0]]) if cross.size else math.nan
    return {"tau": tau, "r": r, "phi": ph, "v": v, "x": x, "y": yy, "tau_to_horizon": tau_h,
            "v_to_horizon": v_h, "turns_to_horizon": float(abs(ph[cross[0]] - phi0) / (2 * math.pi)) if cross.size else math.nan}


def _rhs_bl(tau, y, a, E, L):
    """y = (r, phi, t) in Boyer-Lindquist coordinates (outside the horizon only)."""
    r = y[0]
    delta = r * r - 2.0 * r + a * a
    P = E * (r * r + a * a) - a * L
    R = max(P * P - delta * (r * r + (L - a * E) ** 2), 0.0)
    return [-math.sqrt(R) / (r * r), (-(a * E - L) + a * P / delta) / (r * r),
            (-a * (a * E - L) + (r * r + a * a) * P / delta) / (r * r)]


def plunge_bl(a, E, L, r_start, eps=2e-3, n=6000):
    """The same infall in Boyer-Lindquist coordinates, whose angle phi is measured against
    the distant stars, stopped at r = r_+ + eps. Seen this way a falling body never quite
    reaches the horizon: it is swept round ever faster, approaching the horizon's own
    angular velocity Omega_H (the coordinate time t diverges)."""
    rp, _ = horizons(a)

    def stop(tau, y, *args):
        return y[0] - (rp + eps)
    stop.terminal = True
    sol = integrate.solve_ivp(_rhs_bl, (0.0, 1e5), [r_start, 0.0, 0.0], args=(a, E, L), method="DOP853",
                              rtol=1e-10, atol=1e-12, events=stop, dense_output=True, max_step=0.05)
    T = sol.t[-1]
    # phi runs away near the horizon (d phi/dr ~ 1/Delta): sample geometrically toward the end
    tau = np.unique(np.concatenate([np.linspace(0.0, T, n), T - np.geomspace(T, T * 1e-13, n)]))
    r, ph, t = sol.sol(tau)
    return {"tau": tau, "r": r, "phi": ph, "t": t, "x": r * np.cos(ph), "y": r * np.sin(ph)}


def zero_L_drop(a, r_start=20.0, r_stop=None):
    """A body dropped from rest at infinity with zero angular momentum (E = 1, L = 0)."""
    _, rm = horizons(a)
    return plunge(a, 1.0, 0.0, r_start, r_stop if r_stop is not None else 1.02 * rm)


def isco_plunge(a, dr: float = 1e-3, r_stop=None):
    """The plunge of gas that has reached the ISCO: E, L of the ISCO, starting just inside it."""
    ri = isco(a)
    E, L = circular_EL(ri, a)
    _, rm = horizons(a)
    out = plunge(a, E, L, ri - dr, r_stop if r_stop is not None else 1.02 * rm, max_tau=4e4)
    out.update({"r_isco": ri, "E": E, "L": L})
    return out


def penrose_max_efficiency():
    """Largest fraction of mass extractable from a maximally spinning hole: 1 - 1/sqrt(2)."""
    return 1.0 - 1.0 / math.sqrt(2.0)
