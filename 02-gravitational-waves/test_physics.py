"""Self-checks for the physics module, against cases with known answers.

Run with ``python test_physics.py``. No test framework and no network: these
check the formulas, not the data.

The most satisfying block is the second one. The discovery paper's text walks
from "chirp mass ~30 Msun" to "total mass >~70 Msun", "Schwarzschild radii
>~210 km" and "only ~350 km apart" -- and these formulas have to reproduce
each of those numbers, because the black-hole argument this project makes is
exactly that chain.
"""

from __future__ import annotations

import sys

import numpy as np

import gw_physics as gp

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float, unit: str = "") -> None:
    """Assert ``got`` is within ``tol`` (fractional) of ``expected``."""
    rel = abs(got - expected) / abs(expected) if expected else abs(got)
    status = "ok  " if rel <= tol else "FAIL"
    if rel > tol:
        FAILURES.append(label)
    print(f"  [{status}] {label:<50} {got:>12.5g} vs {expected:>10.5g} {unit:<7} "
          f"({rel * 100:.2f}% off, tol {tol * 100:g}%)")


def truth(label: str, condition: bool) -> None:
    if not condition:
        FAILURES.append(label)
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label}")


def main() -> int:
    print("Units")
    # G Msun / c^3 and G Msun / c^2: the Sun's mass as a time and as a length.
    check("G Msun / c^3", gp.T_SUN * 1e6, 4.9255, 1e-3, "us")
    check("G Msun / c^2", gp.L_SUN, 1476.6, 1e-3, "m")
    check("Schwarzschild radius of the Sun", gp.schwarzschild_radius(1.0) / 1e3, 2.953, 1e-3, "km")

    print("\nThe discovery paper's own back-of-envelope chain (PRL 116, 061102)")
    # "chirp mass of ~30 Msun, implying that the total mass ... is >~70 Msun"
    check("minimum total mass for Mc = 30", gp.min_total_mass(30.0), 70.0, 0.02, "Msun")
    # "This bounds the sum of the Schwarzschild radii ... to 2GM/c^2 >~ 210 km"
    check("2GM/c^2 for M = 70 Msun", gp.schwarzschild_radius(70.0) / 1e3, 210.0, 0.02, "km")
    # "equal Newtonian point masses orbiting at this frequency [75 Hz orbital,
    # 150 Hz GW] would be only ~350 km apart"
    check("separation at f_GW = 150 Hz, M = 70", gp.keplerian_separation(150.0, 70.0) / 1e3,
          350.0, 0.02, "km")
    # "peak gravitational-wave luminosity of 3.6e56 erg/s, equivalent to 200 Msun c^2/s"
    check("3.6e56 erg/s in Msun c^2 / s", gp.erg_per_s_to_msun_c2_per_s(3.6e56), 200.0, 0.02)

    print("\nChirp mass")
    # The canonical binary neutron star: two 1.4 Msun stars.
    check("Mc(1.4, 1.4)", gp.chirp_mass(1.4, 1.4), 1.2188, 1e-4, "Msun")
    check("Mc is symmetric", gp.chirp_mass(36.0, 29.0), gp.chirp_mass(29.0, 36.0), 1e-12)
    check("eta for equal masses", gp.symmetric_mass_ratio(10.0, 10.0), 0.25, 1e-12)
    # For equal masses the minimum-total-mass bound is attained exactly.
    check("min total mass is tight for m1 = m2", gp.min_total_mass(gp.chirp_mass(35, 35)), 70.0, 1e-9)
    m = gp.partner_mass(30.0, 1.4)
    check("partner_mass round-trips", gp.chirp_mass(1.4, m), 30.0, 1e-8, "Msun")

    print("\nThe Newtonian chirp")
    for mc in (1.2, 30.0):
        f = 57.0
        fdot = gp.fdot_newtonian(f, mc)
        check(f"Mc from (f, fdot) round-trips, Mc={mc}", gp.chirp_mass_from_fdot(f, fdot), mc, 1e-10)
        check(f"Mc from slope round-trips, Mc={mc}", gp.chirp_mass_from_slope(gp.chirp_slope(mc)),
              mc, 1e-10)
        tau = gp.time_to_coalescence(f, mc)
        check(f"frequency_at inverts time_to_coalescence, Mc={mc}", gp.frequency_at(tau, mc), f, 1e-10)
        # tau ~ f^(-8/3) implies tau * fdot = (3/8) f exactly.
        check(f"tau * fdot = 3f/8, Mc={mc}", tau * fdot, 3.0 * f / 8.0, 1e-10)
    # f^(-8/3) must fall linearly in time: check against the integrated law.
    mc, tau = 30.0, np.array([0.05, 0.10, 0.20])
    lhs = gp.frequency_at(tau, mc) ** (-8.0 / 3.0)
    check("f^(-8/3) is linear in time, slope = chirp_slope",
          np.polyfit(-tau, lhs, 1)[0], gp.chirp_slope(mc), 1e-9)
    # A 1.4+1.4 neutron-star binary enters LIGO's band at 10 Hz roughly 17
    # minutes (~1000 s) before merger -- the number usually quoted.
    check("BNS time in band from 10 Hz", gp.time_to_coalescence(10.0, 1.2188), 1000.0, 0.02, "s")
    # GW150914 in the paper's words: "Over 0.2 s ... from 35 to 150 Hz". The
    # leading-order law with Mc ~ 30 must put 35 Hz roughly 0.2 s before merger.
    check("Mc=30: time from 35 Hz to merger", gp.time_to_coalescence(35.0, 30.0), 0.2, 0.25, "s")

    print("\nOrbits and compactness")
    # Innermost stable circular orbit: ~4400 Hz for a one-solar-mass system.
    check("f_ISCO for 1 Msun", gp.isco_frequency(1.0), 4397.0, 1e-3, "Hz")
    # Two Suns in contact orbit in about 6 hours -> f_GW ~ 1e-4 Hz.
    check("two Suns touching: f_GW", gp.contact_frequency(2.0, 2 * gp.R_SUN), 1.0e-4, 0.05, "Hz")
    # Sirius B: mass 1.018 Msun (Bond et al. 2017), radius 0.0081 Rsun (Joyce et al. 2018).
    check("white dwarf radius, Sirius B", gp.radius_white_dwarf(1.018) / gp.R_SUN, 0.0081, 0.10, "Rsun")
    truth("white dwarfs shrink as mass grows",
          gp.radius_white_dwarf(1.2) < gp.radius_white_dwarf(0.6) < gp.radius_white_dwarf(0.3))
    check("orbital speed at ISCO is c/sqrt(6)", gp.orbital_speed(gp.isco_frequency(20.0), 20.0),
          1.0 / np.sqrt(6.0), 1e-9, "c")

    print("\nDetector geometry (constants from LALSuite)")
    h1, l1 = gp.DETECTORS["H1"], gp.DETECTORS["L1"]
    # LAL also publishes each vertex in Earth-fixed metres; our WGS-84
    # conversion of the published latitude/longitude/elevation must agree.
    lal_vertex = {
        "H1": (-2.16141492636e06, -3.83469517889e06, 4.60035022664e06),
        "L1": (-7.42760447238e04, -5.49628371971e06, 3.22425701744e06),
    }
    for name, det in (("H1", h1), ("L1", l1)):
        offset = np.linalg.norm(det.vertex - np.array(lal_vertex[name]))
        truth(f"{name} vertex matches LAL's Earth-fixed position ({offset * 1e3:.1f} mm off)",
              offset < 0.01)
        for arm in "xy":
            check(f"{name} {arm} arm is a unit vector", np.linalg.norm(det.arm(arm)), 1.0, 1e-12)
        truth(f"{name} arms are perpendicular (x.y = {det.arm('x') @ det.arm('y'):.1e})",
              abs(det.arm("x") @ det.arm("y")) < 1e-3)
    # "The 15-ms window is determined by the 10-ms intersite propagation time"
    check("H1-L1 light travel time", gp.light_travel_time(h1, l1) * 1e3, 10.0, 0.01, "ms")
    check("a detector overlaps itself fully", gp.tensor_overlap(h1, h1), 1.0, 1e-12)
    # Rotate H1 by 90 degrees about the vertical: x and y swap, D flips sign.
    rotated = gp.Detector("H1 rotated", h1.latitude, h1.longitude, h1.elevation,
                          h1.x_azimuth + np.pi / 2, h1.y_azimuth + np.pi / 2, 0.0, 0.0)
    flat = gp.Detector("H1 flat", h1.latitude, h1.longitude, h1.elevation,
                       h1.x_azimuth, h1.y_azimuth, 0.0, 0.0)
    check("co-located detector rotated 90 deg: overlap", gp.tensor_overlap(flat, rotated), -1.0, 1e-9)
    ov = gp.tensor_overlap(h1, l1)
    truth(f"H1 and L1 are nearly anti-aligned (overlap {ov:+.3f})", -1.0 < ov < -0.7)

    rng = np.random.default_rng(1)
    ra, dec = rng.uniform(0, 2 * np.pi, 20000), np.arcsin(rng.uniform(-1, 1, 20000))
    dt = gp.arrival_delay(h1, l1, ra, dec, gmst=1.234)
    travel = gp.light_travel_time(h1, l1)
    truth("no sky position gives a delay beyond the light travel time",
          np.all(np.abs(dt) <= travel * (1 + 1e-12)))
    check("delays reach the light travel time", np.abs(dt).max(), travel, 0.01, "s")
    truth("delay(a, b) = -delay(b, a)",
          np.allclose(dt, -gp.arrival_delay(l1, h1, ra, dec, gmst=1.234)))

    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: {', '.join(FAILURES)}")
        return 1
    print("All physics checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
