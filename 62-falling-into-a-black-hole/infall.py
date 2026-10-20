"""Falling in: clocks, tides and what a far-away friend receives. Schwarzschild, radial.

Units G = c = M = 1 unless a function says otherwise (multiply times by GM/c^3,
lengths by GM/c^2). The horizon is at r = 2.

Two kinds of faller:

* the **raindrop**, falling from rest at infinity (energy per unit mass E = 1).
  Its speed relative to hovering observers is beta = sqrt(2/r), the escape
  speed, and its own clock is Gullstrand-Painleve time:
  tau_left(r) = (sqrt 2 / 3) r^(3/2) to reach r = 0, exactly as in Newton.
* a body **dropped from rest** at r0 (E = sqrt(1 - 2/r0)). Its fall is a
  cycloid, r = (r0/2)(1 + cos eta), tau = sqrt(r0^3/8)(eta + sin eta), and it
  reaches r = 0 after pi sqrt(r0^3/8), Newton's free-fall time again.

Inside the horizon the longest any body can last is pi (in units of GM/c^3),
reached by falling from rest just above the horizon; the raindrop gets 4/3.
"""

from __future__ import annotations

import math

import numpy as np
from scipy import integrate, optimize

SQRT2 = math.sqrt(2.0)


# ---------------------------------------------------------------------------
# The raindrop (from rest at infinity)
# ---------------------------------------------------------------------------

def rain_tau_left(r):
    """Proper time to the singularity from radius r."""
    return SQRT2 / 3.0 * np.asarray(r, float) ** 1.5


def rain_schwarzschild_t(r):
    """Schwarzschild (far-away clock) time along the raindrop, up to a constant; r > 2.

    dt/dr = -sqrt(r/2) r/(r - 2). With x = sqrt(r/2):
    t = -4 [x^3/3 + x + (1/2) ln((x - 1)/(x + 1))] + const. Diverges at the horizon.
    """
    x = np.sqrt(np.asarray(r, float) / 2.0)
    return -4.0 * (x**3 / 3.0 + x + 0.5 * np.log((x - 1.0) / (x + 1.0)))


def rain_ef_v(r, tau):
    """Ingoing Eddington-Finkelstein time v = t + r* along the raindrop (regular at r = 2).

    v = T + r - 2 sqrt(2 r) + 4 ln(sqrt(r/2) + 1), with GP time T = tau.
    """
    r = np.asarray(r, float)
    return tau + r - 2.0 * np.sqrt(2.0 * r) + 4.0 * np.log(np.sqrt(r / 2.0) + 1.0)


# ---------------------------------------------------------------------------
# Dropped from rest at r0
# ---------------------------------------------------------------------------

def drop_r_tau(eta, r0):
    """Cycloid: radius and proper time (from release) at parameter eta in [0, pi]."""
    eta = np.asarray(eta, float)
    return 0.5 * r0 * (1.0 + np.cos(eta)), math.sqrt(r0**3 / 8.0) * (eta + np.sin(eta))


def drop_tau_total(r0):
    """Proper time from release at r0 to the singularity: pi sqrt(r0^3 / 8)."""
    return math.pi * math.sqrt(r0**3 / 8.0)


def drop_tau_inside(r0):
    """Proper time from the horizon to the singularity for a body dropped from rest at r0 > 2."""
    eta_h = math.acos(4.0 / r0 - 1.0)
    return math.sqrt(r0**3 / 8.0) * (math.pi - eta_h - math.sin(eta_h))


def drop_schwarzschild_t(eta, r0):
    """Schwarzschild time since release (MTW eq. 31.10), valid outside the horizon."""
    eta = np.asarray(eta, float)
    k = math.sqrt(r0 / 2.0 - 1.0)
    tan = np.tan(eta / 2.0)
    return (2.0 * np.log(np.abs((k + tan) / (k - tan)))
            + 2.0 * k * (eta + r0 / 4.0 * (eta + np.sin(eta))))


def drop_schwarzschild_t_numeric(r, r0):
    """Same as drop_schwarzschild_t, by quadrature of dt/dr = -E / ((1 - 2/r) sqrt(E^2 - 1 + 2/r))."""
    e = math.sqrt(1.0 - 2.0 / r0)

    def f(x):
        return e / ((1.0 - 2.0 / x) * math.sqrt(max(e * e - 1.0 + 2.0 / x, 1e-300)))

    # substitute x = r0 - s^2 near the turning point to remove the 1/sqrt singularity
    def g(s):
        return 2.0 * s * f(r0 - s * s)

    val, _ = integrate.quad(g, 0.0, math.sqrt(r0 - r), epsabs=1e-12, epsrel=1e-11, limit=400)
    return val


# ---------------------------------------------------------------------------
# Tides and hovering
# ---------------------------------------------------------------------------

def tidal_stretch_si(r_m, mass_msun, length_m, gm_sun):
    """Head-to-foot acceleration difference 2GM L / r^3 (m/s^2) for a radial faller (exact in GR)."""
    return 2.0 * gm_sun * mass_msun * length_m / np.asarray(r_m, float) ** 3


def tidal_squeeze_si(r_m, mass_msun, width_m, gm_sun):
    """Sideways compression GM w / r^3 (m/s^2)."""
    return gm_sun * mass_msun * width_m / np.asarray(r_m, float) ** 3


def hover_acceleration(r):
    """Proper acceleration needed to hover at r: 1 / (r^2 sqrt(1 - 2/r)) in units of c^4/(GM)."""
    r = np.asarray(r, float)
    return 1.0 / (r**2 * np.sqrt(1.0 - 2.0 / r))


# ---------------------------------------------------------------------------
# What a friend hovering far above receives (radial light, same line as the fall)
# ---------------------------------------------------------------------------

def tortoise(r):
    """r* = r + 2 ln(r/2 - 1), outside the horizon."""
    r = np.asarray(r, float)
    return r + 2.0 * np.log(r / 2.0 - 1.0)


def friend_signal(r_emit, t_emit, r_friend, energy=1.0):
    """Arrival time (friend's own clock) and 1 + z of radial light from a radial faller.

    The faller has energy per unit mass ``energy`` (1 = raindrop) and is at
    (t_emit, r_emit), r_emit > 2, in Schwarzschild coordinates. Light climbs along an
    outgoing radial null line: t_arrive = t_emit + r*(r_friend) - r*(r_emit). The
    friend hovers at r_friend; their clock runs at sqrt(1 - 2/r_friend) of t.
    1 + z = sqrt(1 - 2/r_f) [E + sqrt(E^2 - (1 - 2/r_e))] / (1 - 2/r_e)
    (= sqrt(1 - 2/r_f) / (1 - beta) for the raindrop).
    """
    r_emit = np.asarray(r_emit, float)
    t_arr = np.asarray(t_emit, float) + tortoise(r_friend) - tortoise(r_emit)
    lapse_f = math.sqrt(1.0 - 2.0 / r_friend)
    f = 1.0 - 2.0 / r_emit
    one_plus_z = lapse_f * (energy + np.sqrt(energy**2 - f)) / f
    return t_arr * lapse_f, one_plus_z


def rain_received(r_friend: float, r_start: float, n: int = 4000):
    """The raindrop seen by a friend hovering at r_friend, from r_start down to the horizon.

    Returns arrays: faller's proper time tau (since r_start), friend's clock time
    (since the first signal), 1 + z, and the radius of emission. The radial
    redshift grows as exp(t / 4M) at late times.
    """
    # sample densely toward the horizon: r - 2 geometric
    dr = np.geomspace(r_start - 2.0, 1e-12, n)
    r = 2.0 + dr
    tau = rain_tau_left(r_start) - rain_tau_left(r)
    t = rain_schwarzschild_t(r) - rain_schwarzschild_t(r_start)
    t_f, z1 = friend_signal(r, t, r_friend)
    return tau, t_f - t_f[0], z1, r


def late_time_efold(t_friend, one_plus_z):
    """Fit ln(1 + z) = t / tau_e + c over the last part of the signal; returns tau_e (units of M)."""
    sel = one_plus_z > 50
    slope, _ = np.polyfit(t_friend[sel], np.log(one_plus_z[sel]), 1)
    return 1.0 / slope


# ---------------------------------------------------------------------------
# Kruskal-Szekeres and Penrose coordinates
# ---------------------------------------------------------------------------

def kruskal_from_v_r(v, r):
    """Kruskal (U, V) from ingoing EF time v and radius r (works inside and outside).

    V = exp(v/4), and U V = (1 - r/2) exp(r/2), so U = (1 - r/2) exp(r/2) / V.
    Outside the horizon U < 0, inside U > 0; r = 0 is U V = 1.
    """
    v = np.asarray(v, float)
    r = np.asarray(r, float)
    V = np.exp(v / 4.0)
    U = (1.0 - r / 2.0) * np.exp(r / 2.0) / V
    return U, V


def penrose(U, V, p: float = 1.0):
    """Compactified coordinates: time up, space to the right.

    Ut = arctan(sign(U)|U|^p), Vt likewise. p = 1 is the textbook choice; any p > 0
    is an equally valid compactification (each null coordinate is rescaled
    monotonically, so light still travels at 45 degrees and r = 0, where U V = 1,
    is still the horizontal line Ut + Vt = pi/2). A small p undoes the exponential
    crowding of the far exterior into the corners, so distant observers are visible.
    """
    U = np.asarray(U, float)
    V = np.asarray(V, float)
    pu = np.arctan(np.sign(U) * np.abs(U) ** p)
    pv = np.arctan(np.sign(V) * np.abs(V) ** p)
    return pv - pu, pv + pu     # x, y (both in units where null lines are at 45 degrees)


def static_world_line(r, t):
    """Kruskal (U, V) of a hovering observer at radius r > 2 at Schwarzschild times t."""
    t = np.asarray(t, float)
    rs = float(tortoise(r))
    return -np.exp(-(t - rs) / 4.0), np.exp((t + rs) / 4.0)


# ---------------------------------------------------------------------------
# Spaghettification: where the stretch reaches a limit
# ---------------------------------------------------------------------------

def radius_for_stretch(mass_msun, length_m, accel_si, gm_sun):
    """Radius (m) at which 2GM L / r^3 equals accel_si."""
    return (2.0 * gm_sun * mass_msun * length_m / accel_si) ** (1.0 / 3.0)


def critical_mass_for_horizon_stretch(length_m, accel_si, gm_sun, c):
    """Black hole mass (Msun) whose horizon stretch across length_m equals accel_si.

    At r = 2GM/c^2 the stretch is L c^6 / (4 (GM)^2): M = sqrt(L c^6 / (4 a)) / GM_sun.
    """
    return math.sqrt(length_m * c**6 / (4.0 * accel_si)) / gm_sun


def drop_tau_inside_max():
    """sup over r0 of drop_tau_inside: pi, approached as r0 -> 2 (released just above the horizon)."""
    return math.pi


def rain_tau_inside():
    return float(rain_tau_left(2.0))   # 4/3


def solve_r_at_tau(tau_left_target):
    """Raindrop radius when its clock shows tau_left_target to the singularity."""
    return (3.0 * np.asarray(tau_left_target, float) / SQRT2) ** (2.0 / 3.0)


def drop_r_of_tau(tau, r0):
    """Radius at proper time tau since release from r0 (inverts the cycloid)."""
    def f(eta, target):
        return math.sqrt(r0**3 / 8.0) * (eta + math.sin(eta)) - target
    out = []
    for tt in np.atleast_1d(tau):
        eta = optimize.brentq(f, 0.0, math.pi, args=(tt,))
        out.append(0.5 * r0 * (1 + math.cos(eta)))
    return np.array(out)
