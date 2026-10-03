"""Self-checks for the physics module, against cases with known answers.

Run with ``python test_physics.py``. No test framework and no network: these
check the formulas, not the data.

Every point on the HR diagram passes through ``absolute_magnitude``, and every
velocity through ``tangential_velocity``. A slip in either -- arcseconds for
milliarcseconds, a sign in the distance modulus -- would move every star in a
way that still looks like a perfectly plausible diagram. These catch that.
"""

from __future__ import annotations

import sys

import numpy as np

import hr_physics as hp
from config import SUN_ABS_G

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float, unit: str = "") -> None:
    """Assert ``got`` is within ``tol`` (absolute) of ``expected``."""
    ok = abs(got - expected) <= tol
    if not ok:
        FAILURES.append(label)
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<50} {got:>11.5g} vs {expected:>10.5g} "
          f"{unit:<6} (tol {tol:g})")


def check_raises(label: str, fn) -> None:
    try:
        fn()
    except ValueError:
        print(f"  [ok  ] {label}")
        return
    FAILURES.append(label)
    print(f"  [FAIL] {label}: no ValueError raised")


def main() -> int:
    print("Constants")
    # 1 au/yr expressed in km/s: the factor in every tangential velocity.
    check("1 au/yr in km/s", hp.KM_S_PER_AU_PER_YR, 4.74047, 1e-5, "km/s")

    print("\nDistance and the distance modulus")
    check("10 mas -> 100 pc", float(hp.distance_pc(10.0)), 100.0, 1e-9, "pc")
    # At exactly 10 pc (parallax 100 mas), absolute and apparent magnitude agree
    # by definition.
    check("M = m at 10 pc", float(hp.absolute_magnitude(7.3, 100.0)), 7.3, 1e-12, "mag")
    # Each factor of 10 in distance is 5 magnitudes.
    check("M at 100 pc is m - 5", float(hp.absolute_magnitude(12.0, 10.0)), 7.0, 1e-12, "mag")
    # The Sun, M_G = 4.67, seen from 1 au must come out near its apparent G of
    # about -26.9 (V = -26.76 and G - V = -0.17 for a G2 dwarf).
    au_in_pc = 1.0 / 206264.806
    check("Sun's apparent G from 1 au", float(hp.apparent_magnitude(SUN_ABS_G, au_in_pc)),
          -26.90, 0.05, "mag")
    # Round trip: M -> m at some distance -> M from the matching parallax.
    m = hp.apparent_magnitude(10.0, 37.0)
    check("round trip M -> m -> M", float(hp.absolute_magnitude(m, 1000.0 / 37.0)),
          10.0, 1e-12, "mag")
    check_raises("negative parallax rejected", lambda: hp.absolute_magnitude(10.0, -1.0))
    check_raises("zero parallax rejected", lambda: hp.distance_pc(0.0))

    print("\nHow good the magnitudes are")
    # A 10% parallax error is a (5 / ln 10) * 0.1 = 0.217 mag error in M.
    check("sigma_M at parallax/error = 10", float(hp.absolute_magnitude_error(10.0)),
          0.2171, 1e-4, "mag")
    # Correcting Gaia's -0.017 mas parallax zero-point shifts a 100 pc star by
    # only ~0.004 mag -- far below anything visible on the diagram.
    check("zero-point shift at 10 mas", float(hp.parallax_zero_point_shift(10.0)),
          0.00369, 2e-5, "mag")

    print("\nBrightness relative to the Sun")
    check("M_G = M_sun is 1 Sun", float(hp.luminosity_relative_to_sun(SUN_ABS_G, SUN_ABS_G)),
          1.0, 1e-12)
    # Five magnitudes is exactly a factor of 100 in brightness.
    check("5 mag brighter is 100 Suns",
          float(hp.luminosity_relative_to_sun(SUN_ABS_G - 5, SUN_ABS_G)), 100.0, 1e-9)

    print("\nUnresolved binaries")
    # Two identical stars together are twice as bright: 2.5 log10(2) = 0.753 mag.
    check("twin binary offset", float(hp.combined_magnitude_offset(1.0)), -0.75257, 1e-5, "mag")
    check("faint companion, 1% of the light",
          float(hp.combined_magnitude_offset(0.01)), -0.0108, 1e-4, "mag")

    print("\nTangential velocity")
    # Barnard's Star: 10.39 arcsec/yr of proper motion at 1.83 pc is the
    # textbook ~90 km/s. Gaia DR3 values.
    v = hp.tangential_velocity(-801.551, 10362.394, 546.9759)
    check("Barnard's Star", float(v), 90.07, 0.05, "km/s")
    # 1 mas/yr at 1 mas parallax (1 kpc) is 4.74 km/s, by definition of the constant.
    check("1 mas/yr at 1 kpc", float(hp.tangential_velocity(1.0, 0.0, 1.0)), 4.74047, 1e-5, "km/s")
    # The two proper-motion components add in quadrature (3-4-5 triangle).
    check("components add in quadrature", float(hp.tangential_velocity(3.0, 4.0, 1.0)),
          5 * 4.74047, 1e-4, "km/s")

    print("\nCorrected flux excess C* (Riello et al. 2021)")
    # A star sitting exactly on the expected colour relation has C* = 0.
    x = 1.0
    expected = 1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3
    check("C* = 0 on the relation", float(hp.flux_excess_corrected(x, expected)), 0.0, 1e-12)
    # The published piecewise polynomial is continuous at its two colour
    # boundaries. A mistyped coefficient would open a visible step there.
    for edge in (0.5, 4.0):
        below = float(hp.flux_excess_corrected(edge - 1e-9, 0.0))
        above = float(hp.flux_excess_corrected(edge + 1e-9, 0.0))
        check(f"continuous at BP-RP = {edge}", above - below, 0.0, 1e-4)
    # The tolerance grows towards faint magnitudes.
    check("sigma_C* floor at bright G", float(hp.flux_excess_sigma(5.0)), 0.0060, 1e-4)
    check("sigma_C* at G = 20", float(hp.flux_excess_sigma(20.0)), 0.078, 0.002)
    s = hp.flux_excess_sigma(np.linspace(4, 21, 200))
    ok = bool(np.all(np.diff(s) > 0))
    if not ok:
        FAILURES.append("sigma_C* monotonic")
    print(f"  [{'ok  ' if ok else 'FAIL'}] sigma_C* rises monotonically with G")

    print("\nVelocity dispersion from proper motions alone (Dehnen & Binney 1998)")
    # A known answer built from random numbers: stars in random directions,
    # with a known mean velocity and isotropic random motions of 20 km/s per
    # component, seen only through their projection onto the sky. The method
    # must recover the mean, and S^2 = <|p - A<v>|^2> must equal 2 sigma^2,
    # since two of the three components are visible. (A test of the method,
    # not data: nothing here reaches a figure.)
    from analysis import proper_motion_dispersion

    rng = np.random.default_rng(1)
    n = 200_000
    rhat = rng.normal(size=(n, 3))
    rhat /= np.linalg.norm(rhat, axis=1)[:, None]
    v_true = np.array([-10.0, -20.0, -7.0])
    v = v_true + rng.normal(scale=20.0, size=(n, 3))
    p = v - np.sum(v * rhat, axis=1)[:, None] * rhat
    fit = proper_motion_dispersion(p, rhat, clip=np.inf)
    for axis, name in enumerate("UVW"):
        check(f"mean velocity {name} recovered", fit["v_mean"][axis], v_true[axis], 0.3, "km/s")
    check("S = sqrt(2) x 20 km/s", fit["S"], 20.0 * np.sqrt(2.0), 0.2, "km/s")

    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: {', '.join(FAILURES)}")
        return 1
    print("All physics checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
