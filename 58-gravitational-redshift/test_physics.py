"""Checks of the physics and numerics, on orbits whose answers we know.

    python -m pytest -q test_physics.py      (or: python test_physics.py)
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

import config as C
import gnss

H = C.EPOCH_STEP_S
C2 = C.C_LIGHT**2


def kepler_orbit(t, a=27_977e3, e=0.162, inc=np.radians(50.0), raan=0.4, argp=1.1, m0=0.3):
    """Analytic two-body position and velocity (GCRS-like frame)."""
    n = np.sqrt(C.GM_EARTH / a**3)
    M = m0 + n * t
    E = M.copy()
    for _ in range(50):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    xo = a * (np.cos(E) - e)
    yo = a * np.sqrt(1 - e**2) * np.sin(E)
    edot = n / (1 - e * np.cos(E))
    vxo = -a * np.sin(E) * edot
    vyo = a * np.sqrt(1 - e**2) * np.cos(E) * edot
    cO, sO, ci, si, cw, sw = np.cos(raan), np.sin(raan), np.cos(inc), np.sin(inc), np.cos(argp), np.sin(argp)
    Rm = np.array([[cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci],
                   [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci],
                   [sw * si, cw * si]])
    r = (Rm @ np.vstack([xo, yo])).T
    v = (Rm @ np.vstack([vxo, vyo])).T
    return r, v, E


def detrend(t, y):
    X = np.vstack([np.ones_like(t), t]).T
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ c


def _no_bodies(n):
    far = np.tile([1e30, 0.0, 0.0], (n, 1))   # tides switched off
    return far, far


def test_derivative_matches_analytic_velocity():
    t = np.arange(289) * H
    r, v, _ = kepler_orbit(t)
    vnum = gnss.derivative(r, H)
    assert np.max(np.abs(vnum - v)) < 1e-6        # m/s, including the one-sided ends


def test_kepler_convention_is_the_integrated_two_body_term():
    """-2 r.v/c^2 equals the integral of -(GM/r + v^2/2)/c^2 (minus its mean rate) for a Kepler orbit."""
    t = np.arange(289) * H
    r, v, E = kepler_orbit(t)
    sun, moon = _no_bodies(len(t))
    # J2 does not enter the comparison below: put the orbit on the equator-free axis is unnecessary,
    # we simply compare mono + kin with the convention term.
    vnum, terms = gnss.relativity_terms(H, r, r, sun, moon)
    diff = detrend(t, terms["kepler_convention"] - terms["mono"] - terms["kin"])
    assert np.max(np.abs(diff)) < 0.02e-12           # 0.02 ps, against a 760 ns swing
    a, e = 27_977e3, 0.162
    analytic = -2 * np.sqrt(C.GM_EARTH * a) * e * np.sin(E) / C2
    assert np.max(np.abs(terms["kepler_convention"] - analytic)) < 1e-15
    # half of the periodic swing is gravitational, half is speed
    g = detrend(t, terms["mono"])
    k = detrend(t, terms["kin"])
    assert abs(np.ptp(g) / np.ptp(k) - 1) < 1e-6


def test_numerical_proper_time_on_a_j2_orbit():
    """5-min samples -> velocities -> Simpson integral reproduce an ODE-integrated proper time."""
    mu, J2, R = C.GM_EARTH, C.J2_EARTH, C.R_EARTH_EQ

    def u_total(x):
        rr = np.linalg.norm(x)
        s = x[2] / rr
        return mu / rr - mu * J2 * R**2 / rr**3 * 0.5 * (3 * s**2 - 1)

    def accel(x):
        rr = np.linalg.norm(x)
        z2 = (x[2] / rr) ** 2
        k = 1.5 * J2 * mu * R**2 / rr**5
        return np.array([-mu * x[0] / rr**3 + k * x[0] * (5 * z2 - 1),
                         -mu * x[1] / rr**3 + k * x[1] * (5 * z2 - 1),
                         -mu * x[2] / rr**3 + k * x[2] * (5 * z2 - 3)])

    def rhs(_, s):
        x, vv = s[:3], s[3:6]
        rate = (u_total(x) + 0.5 * vv @ vv) / C2
        return np.concatenate([vv, accel(x), [-rate]])

    r0, v0, _ = kepler_orbit(np.array([0.0]))
    t = np.arange(289) * H
    sol = solve_ivp(rhs, (0, t[-1]), np.concatenate([r0[0], v0[0], [0.0]]), t_eval=t,
                    method="DOP853", rtol=1e-13, atol=1e-9)
    r = sol.y[:3].T
    tau_exact = detrend(t, sol.y[6])
    sun, moon = _no_bodies(len(t))
    _, terms = gnss.relativity_terms(H, r, r, sun, moon)
    tau_num = detrend(t, terms["mono"] + terms["j2"] + terms["kin"])
    assert np.max(np.abs(tau_num - tau_exact)) < 0.05e-12          # 0.05 ps
    # the IGS two-body convention misses the flattening's effect: tens of ps, twice per orbit
    miss = detrend(t, terms["kepler_convention"]) - tau_num
    assert 20e-12 < np.ptp(miss) < 300e-12


def test_gps_budget_numbers():
    """A circular GPS orbit gains ~38.6 us/day on the geoid: +45.7 (height) - 7.2 (speed)."""
    a = 26_561.75e3
    grav = (C.W0_GEOID - C.GM_EARTH / a) / C2 * 86_400e6
    kin = -0.5 * C.GM_EARTH / a / C2 * 86_400e6
    assert abs(grav - 45.7) < 0.1 and abs(kin + 7.2) < 0.1
    assert abs(grav + kin - 38.5) < 0.15
    # the two effects cancel at r = 3 GM / (2 W0), about 1.5 Earth radii (3,170 km up)
    r0 = 1.5 * C.GM_EARTH / C.W0_GEOID
    assert abs(r0 / 6_371e3 - 1.5) < 0.01


def test_rotation_preserves_lengths():
    t = 1_179_100_800.0 + np.arange(5) * 300.0   # 2017-05-21
    R = gnss.terrestrial_to_celestial(t)
    eye = np.einsum("nij,nkj->nik", R, R)
    assert np.allclose(eye, np.eye(3), atol=1e-14)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
