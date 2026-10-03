"""Self-checks for the physics module, against cases with known answers.

Run with ``python test_physics.py``. No test framework needed, and no network:
these check the formulas, not the data.

The point of these is that every transit result downstream is only as good as
these conversions, and a unit slip in one of them would produce a wrong planet
radius that still looks perfectly plausible on a plot.
"""

from __future__ import annotations

import sys

from astropy import units as u

import transit_physics as tp

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float, unit: str = "") -> None:
    """Assert ``got`` is within ``tol`` (fractional) of ``expected``."""
    rel = abs(got - expected) / abs(expected) if expected else abs(got)
    status = "ok  " if rel <= tol else "FAIL"
    if rel > tol:
        FAILURES.append(label)
    print(f"  [{status}] {label:<46} {got:>12.5g} vs {expected:>10.5g} {unit:<8} "
          f"({rel * 100:.2f}% off, tol {tol * 100:.0f}%)")


def main() -> int:
    print("Unit conversions")
    # The Sun is about 109 Earths across, and about 9.7 Jupiters across.
    check("Rsun / Rearth", tp.RSUN_PER_REARTH, 109.076, 1e-3)
    check("Rsun / Rjup", tp.RSUN_PER_RJUP, 9.7311, 1e-3)

    print("\nTransit depth of the solar system, viewed from outside")
    # A distant observer watching Earth cross the Sun sees an 84 ppm dip. This is
    # the canonical number quoted for why Earth-analogue detection is hard.
    earth_depth = (1.0 / tp.RSUN_PER_REARTH) ** 2
    check("Earth across Sun (ppm)", earth_depth * 1e6, 84.0, 0.02, "ppm")
    jupiter_depth = (1.0 / tp.RSUN_PER_RJUP) ** 2
    check("Jupiter across Sun (%)", jupiter_depth * 100, 1.056, 0.02, "%")

    print("\nDepth -> radius inversion round-trips")
    # Feeding a known radius in should give the same radius back out.
    for r_earth in (1.0, 4.0, 11.2):
        depth = (r_earth / tp.RSUN_PER_REARTH) ** 2
        rp, _ = tp.planet_radius(depth, 1.0 * u.R_sun)
        check(f"round-trip {r_earth} Rearth", tp.to_earth_radii(rp), r_earth, 1e-6, "Re")

    print("\nKepler's third law")
    # Earth's sidereal year around a 1 Msun star must give 1 AU.
    a = tp.semi_major_axis(365.256 * u.day, 1.0)
    check("a(365.256 d, 1 Msun)", a.to(u.au).value, 1.0, 1e-4, "AU")
    # Jupiter: 11.862 yr -> 5.204 AU
    a_j = tp.semi_major_axis(11.862 * 365.25 * u.day, 1.0)
    check("a(11.862 yr, 1 Msun)", a_j.to(u.au).value, 5.204, 2e-3, "AU")
    # Mercury: 87.969 d -> 0.3871 AU
    a_m = tp.semi_major_axis(87.969 * u.day, 1.0)
    check("a(87.969 d, 1 Msun)", a_m.to(u.au).value, 0.38710, 1e-3, "AU")

    print("\nScaled semi-major axis")
    # The Sun-Earth distance is about 215 solar radii.
    check("a/Rstar for Earth", tp.scaled_semi_major_axis(365.256 * u.day, 1.0, 1.0),
          215.03, 1e-3)

    print("\nEquilibrium temperature")
    # Earth's equilibrium temperature with a 0.3 Bond albedo is the textbook
    # 255 K -- the number that makes the greenhouse effect's size obvious.
    teq = tp.equilibrium_temperature(5772 * u.K, 1.0, 365.256 * u.day, 1.0, albedo=0.3)
    check("Earth Teq (A=0.3)", teq.value, 255.0, 0.01, "K")
    # With zero albedo it rises to 278.6 K.
    teq0 = tp.equilibrium_temperature(5772 * u.K, 1.0, 365.256 * u.day, 1.0, albedo=0.0)
    check("Earth Teq (A=0)", teq0.value, 278.6, 0.01, "K")

    print("\nTransit duration")
    # Earth's central transit of the Sun lasts about 13 hours.
    dur = tp.transit_duration_central(365.256 * u.day, 1.0, 1.0)
    check("Earth central transit", dur.to(u.hour).value, 13.0, 0.05, "h")

    print("\nImpact parameter from duration")
    # A transit lasting 86.6% of the central duration implies b = 0.5.
    b = tp.impact_parameter_from_duration(8.66 * u.hour, 10.0 * u.hour)
    check("b from T/Tmax = 0.866", b, 0.5, 0.01)
    # An observed duration longer than the central maximum is unphysical.
    assert tp.impact_parameter_from_duration(11 * u.hour, 10 * u.hour) is None
    print("  [ok  ] duration > central maximum returns None")

    print("\nLimb darkening correction")
    ld = tp.limb_darkening_note(0.01)
    # With u1=0.4, u2=0.26 the disc-averaged intensity is 1 - 0.4/3 - 0.26/6.
    check("disc-averaged intensity", ld["disc_average_intensity"],
          1 - 0.4 / 3 - 0.26 / 6, 1e-9)
    # The naive radius ratio should come out larger than the corrected one.
    assert ld["ratio_naive"] > ld["ratio_corrected"], "correction has the wrong sign"
    print(f"  [ok  ] naive ratio exceeds corrected by "
          f"{ld['overestimate_percent']:.1f}% (correct sign)")

    print("\nSignal-to-noise")
    # Doubling the number of in-transit points raises SNR by sqrt(2).
    s1 = tp.transit_snr(1e-4, 1e-3, 100)
    s2 = tp.transit_snr(1e-4, 1e-3, 200)
    check("SNR scales as sqrt(N)", s2 / s1, 2 ** 0.5, 1e-9)
    check("SNR(1e-4, 1e-3, 100)", s1, 1.0, 1e-9)

    print("\nEnd-to-end: HAT-P-7 b from its published depth")
    # Published: Rstar 1.99 Rsun, Mstar 1.56 Msun, P 2.204735 d, Rp ~1.43 Rjup.
    # The naive (no limb darkening) inversion should land somewhat high, which is
    # exactly the systematic the limb-darkening note warns about.
    d = tp.derive_all(depth=0.0062, period_days=2.204735, stellar_radius_rsun=1.99,
                      stellar_mass_msun=1.56, stellar_teff_k=6441)
    check("HAT-P-7 b a", d.semi_major_axis_au, 0.03790, 0.03, "AU")
    check("HAT-P-7 b Rp (naive, expected high)", d.planet_radius_rjup, 1.43, 0.10, "Rj")

    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: {', '.join(FAILURES)}")
        return 1
    print("All physics checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
