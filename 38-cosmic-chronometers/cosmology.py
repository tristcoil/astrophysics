"""The physics: expansion rate, its slowdown and speed-up, and the long-run future.

Notation. ``a`` is the scale factor (1 today), ``z = 1/a - 1`` the redshift, and
``H = (da/dt)/a`` the expansion rate, in km/s/Mpc. ``E = H/H0``. A universe holds
radiation (``Or``), matter (``Om``), curvature (``Ok = 1 - Or - Om - Ode``) and dark
energy (``Ode``) whose equation of state follows the CPL form
``w(a) = w0 + wa (1 - a)`` (``w0 = -1, wa = 0`` is a cosmological constant).

    E(a)^2 = Or a^-4 + Om a^-3 + Ok a^-2 + Ode f(a),
    f(a)   = a^(-3 (1 + w0 + wa)) exp(-3 wa (1 - a)).

Three derived quantities carry the story of this project:

* ``adot = H / (1 + z) = H a``: how fast two galaxies a fixed comoving distance
  apart (1 Mpc today) separate, in km/s. It falls while gravity decelerates the
  expansion and rises once dark energy takes over.
* ``q = -a (d^2a/dt^2) / (da/dt)^2 = (1 + z) H'/H - 1``: the deceleration
  parameter; negative means accelerating. In terms of the contents,
  ``q = (1/2) sum_i Omega_i(z) (1 + 3 w_i)``.
* the transition redshift ``z_t`` where ``q = 0``: the minimum of ``adot``.

Everything here is vectorised: arrays broadcast, so a grid of universes can be
evaluated against a vector of redshifts in one call.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize

import config as cfg

T_H = cfg.HUBBLE_TIME_GYR  # 1/H in Gyr for H = 1 km/s/Mpc


def omega_r(H0) -> float:
    """Radiation density today for a given H0 (photons + massless neutrinos)."""
    return cfg.OMEGA_R_H2 / (np.asarray(H0) / 100.0) ** 2


def w_of_a(a, w0=-1.0, wa=0.0):
    return w0 + wa * (1.0 - np.asarray(a, dtype=float))


def de_density(a, w0=-1.0, wa=0.0):
    """Dark-energy density relative to today, rho_de(a) / rho_de(1), for CPL."""
    a = np.asarray(a, dtype=float)
    return a ** (-3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * (1.0 - a))


def _ode(Om, Ode, Or):
    """Dark energy today: ``None`` means a flat universe."""
    return (1.0 - Om - Or) if Ode is None else Ode


def E2_of_a(a, Om, Ode=None, w0=-1.0, wa=0.0, Or=0.0):
    """(H/H0)^2 at scale factor a. Far-future code calls this directly: going
    through z = 1/a - 1 loses all precision once a exceeds ~1e16."""
    a = np.asarray(a, dtype=float)
    Ode = _ode(Om, Ode, Or)
    Ok = 1.0 - Om - Ode - Or
    return (Or * a**-4 + Om * a**-3 + Ok * a**-2 + Ode * de_density(a, w0, wa))


def E2(z, Om, Ode=None, w0=-1.0, wa=0.0, Or=0.0):
    """(H/H0)^2 at redshift z. ``Ode=None`` makes the universe flat."""
    return E2_of_a(1.0 / (1.0 + np.asarray(z, dtype=float)), Om, Ode, w0, wa, Or)


def hubble(z, H0, Om, Ode=None, w0=-1.0, wa=0.0, Or=0.0):
    """H(z) in km/s/Mpc. NaN where the model has no real expansion rate."""
    e2 = E2(z, Om, Ode, w0, wa, Or)
    return H0 * np.sqrt(np.where(e2 > 0, e2, np.nan))


def adot(z, H):
    """da/dt in km/s/Mpc: the separation speed of a pair 1 Mpc apart today."""
    return np.asarray(H) / (1.0 + np.asarray(z))


def deceleration_of_a(a, Om, Ode=None, w0=-1.0, wa=0.0, Or=0.0):
    """q = (1/2) sum_i Omega_i (1 + 3 w_i) at scale factor a. Curvature
    (w = -1/3) drops out."""
    a = np.asarray(a, dtype=float)
    Ode = _ode(Om, Ode, Or)
    num = (2.0 * Or * a**-4 + Om * a**-3
           + Ode * de_density(a, w0, wa) * (1.0 + 3.0 * w_of_a(a, w0, wa)))
    return num / (2.0 * E2_of_a(a, Om, Ode, w0, wa, Or))


def deceleration(z, Om, Ode=None, w0=-1.0, wa=0.0, Or=0.0):
    """q(z), the deceleration parameter at redshift z."""
    return deceleration_of_a(1.0 / (1.0 + np.asarray(z, dtype=float)), Om, Ode, w0, wa, Or)


def z_transition_lcdm(Om, OL):
    """Transition redshift of a matter + Lambda universe (any curvature).

    q = 0 where Om (1 + z)^3 = 2 OL. NaN if the universe never accelerates.
    """
    Om = np.asarray(Om, dtype=float)
    OL = np.asarray(OL, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        zt = np.cbrt(2.0 * OL / Om) - 1.0
    return np.where((OL > 0) & (Om > 0) & (zt >= 0), zt, np.nan)


def z_transition_flat(Om):
    """Flat LCDM: z_t = (2 (1 - Om) / Om)^(1/3) - 1, defined for Om < 2/3."""
    return z_transition_lcdm(Om, 1.0 - np.asarray(Om, dtype=float))


def z_transition(Om, Ode=None, w0=-1.0, wa=0.0, Or=0.0, z_max=20.0):
    """Highest redshift at which q changes sign (numerical; any model)."""
    zz = np.linspace(0.0, z_max, 20001)
    q = deceleration(zz, Om, Ode, w0, wa, Or)
    if q[0] >= 0:
        return np.nan
    flips = np.nonzero(np.diff(np.sign(q)) > 0)[0]
    if flips.size == 0:
        return np.nan
    i = flips[-1]
    return optimize.brentq(lambda x: deceleration(x, Om, Ode, w0, wa, Or), zz[i], zz[i + 1])


# ---------------------------------------------------------------------------
# Time: ages, lookback times, and the scale factor from Big Bang to far future
# ---------------------------------------------------------------------------


def expansion_history(H0, Om, Ode=None, w0=-1.0, wa=0.0, Or=None,
                      ln_a_min=-25.0, ln_a_max=12.0, n=74001):
    """Cosmic time t(a) since the Big Bang, in Gyr, on a grid in ln a.

    Integrates dt = d(ln a) / H(a) with the trapezoid rule (Simpson-level
    accuracy is not needed: the grid step is 5e-4 in ln a). The first instant,
    deep in the radiation era, is added analytically: t = a^2 / (2 H0 sqrt(Or)).
    If the expansion would turn around (H^2 < 0), the grid stops there.
    Radiation defaults to photons + neutrinos for the given H0, and a flat
    universe (``Ode=None``) takes its dark energy as 1 - Om - Or.
    """
    Or = float(omega_r(H0)) if Or is None else Or
    ln_a = np.linspace(ln_a_min, ln_a_max, n)
    a = np.exp(ln_a)
    e2 = E2_of_a(a, Om, Ode, w0, wa, Or)
    ok = np.cumprod(e2 > 0).astype(bool)  # stop at the first turnaround
    ln_a, a, e2 = ln_a[ok], a[ok], e2[ok]
    inv_h = T_H / (H0 * np.sqrt(e2))  # 1/H in Gyr
    t = np.concatenate([[0.0], np.cumsum(0.5 * (inv_h[1:] + inv_h[:-1]) * np.diff(ln_a))])
    t += a[0] ** 2 / (2.0 * np.sqrt(Or)) * T_H / H0 if Or > 0 else 0.0
    return {"ln_a": ln_a, "a": a, "t": t, "H": H0 * np.sqrt(e2)}


def age(H0, Om, Ode=None, w0=-1.0, wa=0.0, Or=None) -> float:
    """Age of the universe today in Gyr."""
    hist = expansion_history(H0, Om, Ode, w0, wa, Or, ln_a_max=0.0, n=50001)
    return float(hist["t"][-1])


def lookback_time(z, H0, Om, Ode=None, w0=-1.0, wa=0.0, Or=None):
    """Lookback time to redshift z in Gyr."""
    hist = expansion_history(H0, Om, Ode, w0, wa, Or, ln_a_max=0.0, n=50001)
    t_of_lna = np.interp(-np.log1p(np.asarray(z, dtype=float)), hist["ln_a"], hist["t"])
    return hist["t"][-1] - t_of_lna


def big_rip_time(H0, Om, w, Or=None) -> float:
    """Time from today to the Big Rip (Gyr) for constant phantom w < -1, flat.

    Integrates t = int_{a=1}^{inf} d(ln a)/H out to a = e^60 and adds the
    tail analytically: once phantom energy dominates, H grows as
    a^(3|1+w|/2), so the remaining time beyond a is 2 / (3 |1+w| H(a)).
    """
    if w >= -1:
        return np.inf
    Or = float(omega_r(H0)) if Or is None else Or
    ln_a = np.linspace(0.0, 60.0, 600001)
    a = np.exp(ln_a)
    e2 = E2_of_a(a, Om, None, w, 0.0, Or)
    inv_h = T_H / (H0 * np.sqrt(e2))
    t = np.sum(0.5 * (inv_h[1:] + inv_h[:-1]) * np.diff(ln_a))
    return float(t + 2.0 / (3.0 * abs(1.0 + w)) * inv_h[-1])


def big_rip_time_approx(H0, Om, w) -> float:
    """Caldwell, Kamionkowski & Weinberg (2003): (2/3)|1+w|^-1 H0^-1 (1-Om)^-1/2."""
    return (2.0 / 3.0) / abs(1.0 + w) * (T_H / H0) / np.sqrt(1.0 - Om)


def rip_lead_time(period, w):
    """How long before the Big Rip a bound system of orbital period P comes apart.

    Caldwell et al. (2003): t = P sqrt(2 |1 + 3w|) / (6 pi |1 + w|), independent
    of H0 and Omega_m. (Same units as ``period``.)
    """
    return np.asarray(period) * np.sqrt(2.0 * abs(1.0 + 3.0 * w)) / (6.0 * np.pi * abs(1.0 + w))


def acceleration_end(H0, Om, w0, wa, Or=None, ln_a_max=8.0):
    """For a CPL universe accelerating today: the future a and t at which q
    returns to positive (the acceleration ends). ``None`` if it never does
    before a = e^ln_a_max.
    """
    Or = float(omega_r(H0)) if Or is None else Or
    ln_a = np.linspace(0.0, ln_a_max, 80001)
    a = np.exp(ln_a)
    q = deceleration_of_a(a, Om, None, w0, wa, Or)
    if q[0] >= 0:
        return None
    up = np.nonzero((q[:-1] < 0) & (q[1:] >= 0))[0]
    if up.size == 0:
        return None
    i = up[0]
    f = lambda x: deceleration_of_a(np.exp(x), Om, None, w0, wa, Or)
    ln_a_end = optimize.brentq(f, ln_a[i], ln_a[i + 1])
    hist = expansion_history(H0, Om, None, w0, wa, Or, ln_a_max=ln_a_max + 0.5)
    t_end = np.interp(ln_a_end, hist["ln_a"], hist["t"])
    t_now = np.interp(0.0, hist["ln_a"], hist["t"])
    return {"a": float(np.exp(ln_a_end)), "t_from_now_gyr": float(t_end - t_now),
            "w_at_end": float(w_of_a(np.exp(ln_a_end), w0, wa))}
