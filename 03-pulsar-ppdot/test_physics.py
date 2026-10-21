"""Self-checks for the physics module, against cases with known answers.

Run with ``python test_physics.py``. No test framework needed, and no network:
these check the formulas, not the data.

Every number on the P-Pdot diagram's guide lines, and every age, field and
luminosity quoted in the README, passes through these functions. A slip of 2 in
the field coefficient or a year of 365 vs 365.25 days would not look wrong on a
log-log plot spanning twelve decades -- which is exactly why it is checked here.
"""

from __future__ import annotations

import sys

import numpy as np

import pulsar_physics as pp

FAILURES: list[str] = []
N_CHECKS = 0


def check(label: str, got: float, expected: float, tol: float, unit: str = "") -> None:
    """Assert ``got`` is within ``tol`` (fractional) of ``expected``."""
    global N_CHECKS
    N_CHECKS += 1
    rel = abs(got - expected) / abs(expected) if expected else abs(got)
    status = "ok  " if rel <= tol else "FAIL"
    if rel > tol:
        FAILURES.append(label)
    print(f"  [{status}] {label:<48} {got:>12.5g} vs {expected:>10.5g} {unit:<8} "
          f"({rel * 100:.2f}% off, tol {tol * 100:.1f}%)")


def check_true(label: str, condition: bool) -> None:
    global N_CHECKS
    N_CHECKS += 1
    if not condition:
        FAILURES.append(label)
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label}")


def main() -> int:
    print("Magnetic dipole radiation")
    # B = sqrt(3 c^3 I P Pdot / (8 pi^2 R^6)) with I = 1e45 g cm^2 and R = 10 km
    # must reproduce the 3.2e19 the catalogue (and every textbook) uses.
    check("B coefficient from first principles", pp.dipole_field_coefficient(), 3.2e19,
          1e-3, "G")

    # The Crab pulsar, from the catalogue's timing solution (Lyne et al. 2015):
    # f0 = 29.9469230 Hz, f1 = -3.77535e-10 Hz/s, f2 = 1.1147e-20 Hz/s^2.
    f0, f1, f2 = 29.9469230, -3.77535e-10, 1.1147e-20
    p_crab, pdot_crab = 1.0 / f0, -f1 / f0**2
    print("\nThe Crab pulsar -- against the ATNF catalogue's own derived values")
    check("Crab characteristic age", pp.characteristic_age_yr(p_crab, pdot_crab), 1.26e3,
          0.01, "yr")
    check("Crab surface field", pp.surface_field_gauss(p_crab, pdot_crab), 3.79e12, 0.01, "G")
    check("Crab spin-down luminosity", pp.spin_down_luminosity(p_crab, pdot_crab), 4.46e38,
          0.01, "erg/s")
    # 45 years of Jodrell Bank timing give a mean braking index of 2.342(1).
    check("Crab braking index (Lyne et al. 2015)", pp.braking_index(f0, f1, f2), 2.342,
          1e-3)
    # The light cylinder -- where corotation would reach c -- is ~1600 km out.
    check("Crab light cylinder", pp.light_cylinder_radius_cm(p_crab) / 1e5, 1593, 0.01, "km")

    print("\nThe Vela pulsar -- against the ATNF catalogue's own derived values")
    f0v, f1v = 11.1946499395, -1.5666e-11
    p_vela, pdot_vela = 1.0 / f0v, -f1v / f0v**2
    check("Vela characteristic age", pp.characteristic_age_yr(p_vela, pdot_vela), 1.13e4,
          0.01, "yr")
    check("Vela surface field", pp.surface_field_gauss(p_vela, pdot_vela), 3.38e12, 0.01, "G")
    check("Vela spin-down luminosity", pp.spin_down_luminosity(p_vela, pdot_vela), 6.92e36,
          0.01, "erg/s")

    print("\nAges with a birth period")
    # Born spinning infinitely fast, the true age is the characteristic age.
    check("true age -> char. age as P0 -> 0",
          pp.true_age_yr(p_crab, pdot_crab, 1e-9), pp.characteristic_age_yr(p_crab, pdot_crab),
          1e-9, "yr")
    # Round trip: an age computed from a birth period must give the period back.
    for n in (2.342, 3.0):
        age = pp.true_age_yr(p_crab, pdot_crab, 0.019, braking_index=n)
        check(f"birth-period round trip (n = {n})",
              pp.birth_period_s(p_crab, pdot_crab, age, braking_index=n) * 1e3, 19.0, 1e-9, "ms")
    # An age longer than the spin-down age is impossible for that law.
    check_true("age > spin-down age gives NaN birth period",
               np.isnan(pp.birth_period_s(p_crab, pdot_crab, 5000.0)))

    print("\nGuide lines on the diagram invert their own formulas")
    periods = np.array([1e-3, 0.1, 10.0])
    for age in (1e3, 1e6, 1e9):
        got = pp.characteristic_age_yr(periods, pp.pdot_for_age(periods, age))
        check(f"constant-age line, {age:.0e} yr", float(np.max(np.abs(got / age - 1))) + 1,
              1.0, 1e-12)
    for field in (1e8, 1e12, 1e15):
        got = pp.surface_field_gauss(periods, pp.pdot_for_field(periods, field))
        check(f"constant-field line, {field:.0e} G", float(np.max(np.abs(got / field - 1))) + 1,
              1.0, 1e-12)
    got = pp.spin_down_luminosity(periods, pp.pdot_for_edot(periods, 1e33))
    check("constant-Edot line, 1e33 erg/s", float(np.max(np.abs(got / 1e33 - 1))) + 1, 1.0, 1e-12)

    print("\nDeath lines (Zhang, Harding & Muslimov 2000)")
    # At P = 1 s the dipolar line sits at 10^-14.62.
    check("dipole death line at P = 1 s", pp.death_line_pdot(1.0, "dipole"), 10**-14.62, 1e-12)
    # ZHM2000 single out the 8.5 s PSR J2144-3933 as lying beyond the Chen &
    # Ruderman death valley: below even the lower, multipolar line.
    p_j2144, pdot_j2144 = 8.509827, 4.96e-16
    check_true("J2144-3933 lies below the multipolar death line",
               pdot_j2144 < pp.death_line_pdot(p_j2144, "multipole"))
    check_true("the Crab lies far above the dipolar death line",
               pdot_crab > 1e3 * pp.death_line_pdot(p_crab, "dipole"))

    print("\nShklovskii effect")
    # The standard coefficient: Pdot_shk / P = 2.43e-21 (mu / mas/yr)^2 (d / kpc) s^-1.
    check("Shklovskii coefficient", pp.shklovskii_pdot(1.0, 1.0, 1.0), 2.43e-21, 1e-3, "s^-1")

    print("\nWhat a spin period proves about density")
    # A satellite skimming Earth's surface orbits in 84.34 min. The same
    # formula that bounds a pulsar's density turns that back into Earth's mean
    # density, 5.514 g/cm^3.
    check("Earth's density from a grazing orbit", pp.minimum_mean_density(84.34 * 60), 5.514,
          2e-3, "g/cm3")
    # Hessels et al. (2006): the 716 Hz PSR J1748-2446ad limits a neutron star's
    # radius to below 14.4 km for 1.4 Msun, and below 16 km for up to 2 Msun.
    check("R_max at 716 Hz, 1.4 Msun (Hessels+06)", pp.max_radius_km(716.36, 1.4), 14.4, 0.005,
          "km")
    check("R_max at 716 Hz, 2.0 Msun (Hessels+06)", pp.max_radius_km(716.36, 2.0), 16.0, 0.02,
          "km")
    check("max spin <-> max radius round trip",
          pp.max_spin_frequency_hz(1.4, pp.max_radius_km(716.36, 1.4)), 716.36, 1e-12, "Hz")

    print("\nNon-physical inputs are refused, not silently converted")
    check_true("negative Pdot -> NaN age", np.isnan(pp.characteristic_age_yr(0.005, -1e-20)))
    check_true("zero Pdot -> NaN field", np.isnan(pp.surface_field_gauss(0.001, 0.0)))
    check_true("negative Pdot -> NaN Edot", np.isnan(pp.spin_down_luminosity(0.005, -1e-20)))

    print()
    if FAILURES:
        print(f"{len(FAILURES)} of {N_CHECKS} CHECKS FAILED: {', '.join(FAILURES)}")
        return 1
    print(f"All {N_CHECKS} physics checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
