"""Self-checks for the orbital mechanics, against cases with known answers.

Run with ``python test_physics.py``. No test framework and no network: these
check the formulas, not the data.

Every number in this project -- an interstellar object's speed, the direction it
came from, the mass of the Sun -- passes through these functions, and a slip in
any of them (a missing factor of two, a rotation the wrong way round) would give
answers that still look perfectly plausible on a plot.
"""

from __future__ import annotations

import sys

import numpy as np

import orbit_physics as op

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float, unit: str = "") -> None:
    """Assert ``got`` is within ``tol`` (fractional, or absolute when expected is 0)."""
    err = abs(got - expected) / abs(expected) if expected else abs(got)
    ok = err <= tol
    if not ok:
        FAILURES.append(label)
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<50} {got:>14.8g} vs {expected:>14.8g} "
          f"{unit:<7} (err {err:.1e}, tol {tol:.0e})")


def check_true(label: str, condition: bool) -> None:
    if not condition:
        FAILURES.append(label)
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label}")


def synthetic_orbit(a, e, inc, raan, argp, n_orbits, samples_per_orbit, mu=op.GM_SUN):
    """A perfect Keplerian orbit sampled uniformly in time, via Kepler's equation."""
    period = op.kepler_period(a, mu)
    t = np.linspace(0.0, n_orbits * period, int(n_orbits * samples_per_orbit) + 1)
    E = op.solve_kepler(2 * np.pi * t / period, e)
    nu = op.true_anomaly_from_eccentric(E, e)
    states = [op.state_from_elements(a, e, inc, raan, argp, n_, mu) for n_ in nu]
    r = np.array([s[0] for s in states])
    v = np.array([s[1] for s in states])
    return t + 2451545.0, r, v, period


def main() -> int:
    print("Units and constants")
    check("1 au/day in km/s", op.KMS_PER_AU_DAY, 1731.456837, 1e-9, "km/s")
    # One light-day is 173.14 au -- the number in NASA's Voyager 1 announcement.
    check("one light-day in au", op.LIGHT_DAY_AU, 173.1446, 1e-6, "au")
    # GM_sun in au^3/day^2 is the square of Gauss's gravitational constant
    # k = 0.01720209895, to the 1e-10 level at which DE440 refines it.
    check("GM_sun vs Gauss's k^2", op.GM_SUN, 0.01720209895**2, 1e-9, "au3/d2")

    print("\nSpeeds at Earth's distance")
    # Earth's mean orbital speed is 29.78 km/s; escape from the Sun at 1 au, 42.1.
    check("circular speed at 1 au", op.to_kms(op.circular_speed(1.0)), 29.78, 1e-3, "km/s")
    check("solar escape speed at 1 au", op.to_kms(op.escape_speed(1.0)), 42.12, 1e-3, "km/s")
    check("escape = sqrt(2) x circular", op.escape_speed(3.7) / op.circular_speed(3.7),
          np.sqrt(2), 1e-12)

    print("\nKepler's third law")
    # A 1 au orbit around the Sun alone takes 365.2569 days; Earth's actual
    # sidereal year (365.2564 d) is shorter because the Earth-Moon mass adds to
    # the Sun's.
    check("period of a 1 au orbit", op.kepler_period(1.0), 365.2569, 1e-6, "d")
    # JPL's mean Earth-Moon orbit (a = 1.00000261 au) and the sidereal year
    # (365.256363 d) should imply GM = G(M_sun + M_earth+moon), i.e. GM_sun
    # times (1 + 1/328900.56).
    gm_earth = op.gm_from_orbit(1.00000261, 365.256363)
    check("GM from Earth's orbit / expected", gm_earth / (op.GM_SUN * (1 + 1 / 328900.56)),
          1.0, 2e-5)

    print("\nElements <-> state round trips")
    cases = {
        # Halley-like: highly eccentric and retrograde.
        "ellipse (Halley-like)": (17.83, 0.9671, np.radians(162.2), np.radians(59.1),
                                  np.radians(112.2), np.radians(35.0)),
        # 'Oumuamua-like: hyperbolic, a < 0.
        "hyperbola ('Oumuamua-like)": (-1.272, 1.2011, np.radians(122.7), np.radians(24.6),
                                       np.radians(241.8), np.radians(-60.0)),
    }
    for label, (a, e, inc, raan, argp, nu) in cases.items():
        r, v = op.state_from_elements(a, e, inc, raan, argp, nu)
        el = op.elements_from_state(r, v)
        check(f"{label}: a", float(el.a), a, 1e-10, "au")
        check(f"{label}: e", float(el.e), e, 1e-10)
        check(f"{label}: inclination", float(np.degrees(el.inc)), np.degrees(inc), 1e-10, "deg")
        check(f"{label}: true anomaly", float(el.nu), float(np.mod(nu, 2 * np.pi)), 1e-9, "rad")

    print("\nHyperbolic asymptotes")
    # e = 2 with perihelion along +x: the asymptotes sit at +-120 degrees.
    a, e = -1.0, 2.0
    r, v = op.state_from_elements(a, e, 0.0, 0.0, 0.0, np.radians(10.0))
    radiant, heading, vinf = op.asymptotes(r, v)
    check("heading angle, e=2", float(np.degrees(np.arctan2(heading[1], heading[0]))),
          120.0, 1e-10, "deg")
    check("radiant angle, e=2", float(np.degrees(np.arctan2(radiant[1], radiant[0]))),
          -120.0, 1e-10, "deg")
    check("v_inf = sqrt(GM / -a)", float(vinf), np.sqrt(op.GM_SUN / 1.0), 1e-10, "au/d")
    # Far out along the outbound leg the velocity must point along the heading.
    r_far, v_far = op.state_from_elements(a, e, 0.0, 0.0, 0.0, np.radians(119.999))
    check_true("velocity far out is parallel to the heading (<0.1 deg)",
               np.degrees(op.angle_between(v_far, heading)) < 0.1)
    # And v^2 - v_esc^2 is the same everywhere: energy conservation.
    check("v^2 - v_esc^2 conserved",
          float(op.hyperbolic_excess_speed(r_far, v_far)), float(vinf), 1e-10)

    print("\nGravity-assist bending")
    # r_p v_inf^2 / GM = 1 gives e = 2 and a 60 degree turn.
    check("turn angle when e = 2", np.degrees(op.flyby_turn_angle(1.0, 1.0, 1.0)), 60.0,
          1e-12, "deg")
    check_true("closer pass bends more",
               op.flyby_turn_angle(10.0, 1e5, 1.27e8) > op.flyby_turn_angle(10.0, 1e6, 1.27e8))

    print("\nReference frames")
    # The ecliptic pole sits at RA 18h (270 deg), Dec 90 - 23.44 = +66.56 deg.
    ra, dec = op.radec_from_ecliptic(np.array([0.0, 0.0, 1.0]))
    check("ecliptic pole RA", float(ra), 270.0, 1e-12, "deg")
    check("ecliptic pole Dec", float(dec), 66.560708, 1e-7, "deg")
    # The Galactic centre, RA 266.405 Dec -28.936, must come out at l = b = 0.
    gc = np.array([np.cos(np.radians(-28.93617)) * np.cos(np.radians(266.40499)),
                   np.cos(np.radians(-28.93617)) * np.sin(np.radians(266.40499)),
                   np.sin(np.radians(-28.93617))])
    l, b = op.spherical_deg(op.icrs_to_galactic(gc))
    check("Galactic centre longitude (deg from 0)", float(min(l, 360 - l)), 0.0, 2e-3, "deg")
    check("Galactic centre latitude", float(b), 0.0, 2e-3, "deg")
    # 10 km/s straight at the Galactic centre is pure U.
    v_ecl = op.icrs_to_ecliptic(gc * 10.0 / op.KMS_PER_AU_DAY)
    U, V, W = op.galactic_uvw_kms(v_ecl)
    check("U for motion toward the Galactic centre", float(U), 10.0, 1e-4, "km/s")
    check("V for motion toward the Galactic centre", float(V), 0.0, 1e-3, "km/s")
    back = op.icrs_to_ecliptic(op.ecliptic_to_icrs(np.array([0.3, -0.5, 0.8])))
    check("ecliptic -> ICRS -> ecliptic", float(op.norm(back - [0.3, -0.5, 0.8])), 0.0, 1e-15)

    print("\nInterpolation and event finding")
    t = np.linspace(0, 20, 41)  # step 0.5: coarse on purpose
    got = [c[0] for c in op.crossings(t, np.sin(t), np.cos(t), 0.0)]
    worst = max(abs(g - np.pi * round(g / np.pi)) for g in got)
    check("zeros of sin(t) from 0.5-step samples", worst, 0.0, 1e-5, "")
    mins = op.extrema(t, np.cos(t), -np.sin(t), kind="min")
    check("minimum of cos(t) at pi, 0.5 step", mins[0][0], np.pi, 1e-3)

    # Accuracy is not the point of a coarse grid -- the *order* is. Cubic
    # Hermite with exact derivatives must be fourth order: halving the step
    # should cut the position error ~16x.
    def hermite_error(step):
        tt = np.arange(0, 20 + step, step)
        track = op.Track(tt, np.column_stack([np.sin(tt), np.cos(tt), tt]),
                         np.column_stack([np.cos(tt), -np.sin(tt), np.ones_like(tt)]))
        probe = np.linspace(1.03, 18.97, 400)
        truth = np.column_stack([np.sin(probe), np.cos(probe), probe])
        return float(np.max(op.norm(track.position(probe) - truth)))

    e1, e2 = hermite_error(0.5), hermite_error(0.25)
    check("Hermite error on a 0.5 grid (~h^4/384)", e1, 0.5**4 / 384, 0.35)
    check_true(f"Hermite is fourth order: halving the step cuts error {e1 / e2:.1f}x (~16)",
               12 < e1 / e2 < 20)

    print("\nMeasuring a perfect synthetic orbit")
    # A Mars-like ellipse followed for 5 orbits, sampled 120 times per orbit:
    # the measurement chain must hand back the period and size we put in, and
    # so the GM that generated it.
    a, e = 1.5237, 0.0934
    jd, r, v, period = synthetic_orbit(a, e, np.radians(1.85), np.radians(49.6),
                                       np.radians(286.5), 5, 120)
    m = op.measure_orbit(jd, r, v)
    check("measured sidereal period", m["period_days"], period, 1e-9, "d")
    check("measured semi-major axis", m["a_au"], a, 1e-9, "au")
    check("measured eccentricity", m["e"], e, 1e-8)
    check("recovered GM / true GM", op.gm_from_orbit(m["a_au"], m["period_days"]) / op.GM_SUN,
          1.0, 1e-8)
    # Kepler's second law on the same orbit: the areal velocity is constant.
    av = op.areal_velocity(r, v)
    check("areal velocity spread (max/min - 1)", float(av.max() / av.min() - 1), 0.0, 1e-12)
    # A retrograde, highly eccentric orbit (Halley-like) must work the same way.
    jd, r, v, period = synthetic_orbit(17.83, 0.967, np.radians(162.2), np.radians(59.1),
                                       np.radians(112.2), 2, 4000)
    m = op.measure_orbit(jd, r, v)
    check("retrograde e=0.967: period", m["period_days"], period, 1e-7, "d")
    check("retrograde e=0.967: semi-major axis", m["a_au"], 17.83, 1e-7, "au")

    print("\nA ring of mass (the orbit-averaged pull of another planet)")
    # Far outside a ring, it pulls like a point mass at its centre.
    check("ring far inside the orbit acts like a point mass", op.ring_gm_factor(1e-4), 1.0, 1e-7)
    check("leading term inside: 1 + 3/4 x^2", op.ring_gm_factor(0.1) - 1, 0.75 * 0.01, 0.02)
    check("leading term outside: -(1/2) y^3", op.ring_gm_factor(10.0), -0.5e-3, 0.02)

    # Brute force: 20,000 point masses around a circle of radius b, summed.
    def brute_force_factor(b, r=1.0, n=20000):
        phi = np.linspace(0, 2 * np.pi, n, endpoint=False)
        pts = np.column_stack([b * np.cos(phi), b * np.sin(phi)])
        d = np.array([r, 0.0]) - pts
        pull = -np.sum(d[:, 0] / np.hypot(d[:, 0], d[:, 1]) ** 3) / n  # radial, per unit GM
        return -pull * r**2

    for ratio in (0.3, 0.55, 0.8, 1.6, 3.0):
        check(f"series vs 20,000-point ring, b/r = {ratio}", op.ring_gm_factor(ratio),
              brute_force_factor(ratio), 1e-6)

    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: {', '.join(FAILURES)}")
        return 1
    print("All physics checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
