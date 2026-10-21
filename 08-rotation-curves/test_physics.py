"""Self-checks for the physics module, against cases with known answers.

Run with ``python test_physics.py``. No test framework and no network: these
check the formulas, not the data.

Every dark-matter claim this project makes rests on comparing two speeds, and
both come out of these conversions. A slip of a factor of sqrt(2) in how the
stellar contribution scales, or a kpc/km mix-up in an acceleration, would
produce a wrong answer that still looks entirely plausible on a plot -- so each
formula is pinned to something whose value is known independently.
"""

from __future__ import annotations

import sys

import numpy as np
from astropy import constants as const
from astropy import units as u

import rotation_physics as rp

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float, unit: str = "") -> None:
    """Assert ``got`` is within ``tol`` (fractional) of ``expected``."""
    got = float(got)
    rel = abs(got - expected) / abs(expected) if expected else abs(got)
    status = "ok  " if rel <= tol else "FAIL"
    if rel > tol:
        FAILURES.append(label)
    print(f"  [{status}] {label:<50} {got:>12.5g} vs {expected:>11.5g} {unit:<9} "
          f"({rel * 100:.2f}% off, tol {tol * 100:.2g}%)")


def check_true(label: str, condition: bool) -> None:
    status = "ok  " if condition else "FAIL"
    if not condition:
        FAILURES.append(label)
    print(f"  [{status}] {label}")


def main() -> int:
    print("Constants")
    # G in galaxy units is the one number every velocity-mass conversion uses.
    check("G (kpc km^2 s^-2 Msun^-1)", rp.G_KPC, 4.3009e-6, 1e-4)
    check("1 kpc in metres", rp.KPC_M, 3.0857e19, 1e-4, "m")
    check("GM_sun", rp.GM_SUN_KM3_S2, 1.32712440e11, 1e-8, "km^3/s^2")

    print("\nThe Solar System is Keplerian")
    # Earth's circular speed, and the round trip speed -> enclosed mass -> speed.
    v_earth = rp.keplerian_velocity(1.0, rp.AU_KPC)
    check("circular speed at 1 AU", v_earth, 29.785, 1e-3, "km/s")
    check("mass enclosed by Earth's orbit", rp.enclosed_mass(v_earth, rp.AU_KPC), 1.0, 1e-9, "Msun")
    # The acceleration Earth feels from the Sun, GM/r^2, computed via V^2/r.
    g_earth = float((const.GM_sun / (1 * u.au) ** 2).to(u.m / u.s**2).value)
    check("V^2/r at 1 AU equals GM/r^2", rp.acceleration(v_earth, rp.AU_KPC), g_earth, 1e-9, "m/s^2")
    # NASA Planetary Fact Sheet mean orbital speeds (rounded to 0.1 km/s there),
    # from its distances and eccentricities. The mean speed of an eccentric
    # orbit is slightly below the circular one; Mercury is the stringent case.
    for name, a_1e6km, e, v_nasa in (
        ("Mercury", 57.9, 0.206, 47.4),
        ("Earth", 149.6, 0.017, 29.8),
        ("Jupiter", 778.5, 0.049, 13.1),
        ("Neptune", 4515.0, 0.010, 5.4),
        ("Pluto", 5906.4, 0.244, 4.7),
    ):
        v = rp.mean_orbital_speed(a_1e6km * 1e6, e)
        check(f"{name} mean orbital speed", v, v_nasa, 0.012, "km/s")
    # Without the eccentricity correction Mercury would miss by ~1%.
    naive = rp.mean_orbital_speed(57.9e6, 0.0)
    check_true(f"eccentricity lowers Mercury's mean speed ({naive:.2f} -> "
               f"{rp.mean_orbital_speed(57.9e6, 0.206):.2f} km/s)",
               rp.mean_orbital_speed(57.9e6, 0.206) < naive)

    print("\nCombining the visible components")
    # With Y = 1 and no gas or bulge, the baryonic curve is the disc curve.
    v = np.array([50.0, 100.0, 150.0])
    check("Y=1 pure disc returns Vdisk", rp.signed_sqrt(
        rp.baryonic_velocity_squared(0 * v, v, 0 * v, 1.0, 1.0))[1], 100.0, 1e-12, "km/s")
    # Mass scales V^2, not V: Y = 0.5 must give Vdisk / sqrt(2), not Vdisk / 2.
    check("Y=0.5 scales V by 1/sqrt(2)", rp.signed_sqrt(
        rp.baryonic_velocity_squared(0 * v, v, 0 * v, 0.5, 1.0))[1], 100.0 / np.sqrt(2), 1e-12, "km/s")
    # Components add in quadrature: 30 and 40 make 50.
    check("30 km/s gas + 40 km/s disc", rp.signed_sqrt(
        rp.baryonic_velocity_squared(30.0, 40.0, 0.0, 1.0, 1.0)), 50.0, 1e-12, "km/s")
    # A negative (outward-pulling) gas term must subtract, not add.
    check("50 km/s disc with -30 km/s gas", rp.signed_sqrt(
        rp.baryonic_velocity_squared(-30.0, 50.0, 0.0, 1.0, 1.0)), 40.0, 1e-12, "km/s")
    # Total baryonic mass: 1e10 Lsun of disc at Y=0.5 plus 1e9 Msun of HI.
    check("M_bar = 0.5 L + 1.33 M_HI", rp.baryonic_mass(1e10, 0.0, 1e9, 0.5, 0.7), 6.33e9, 1e-12, "Msun")

    print("\nAn exponential disc (Freeman 1970)")
    r_d, m_d = 3.0, 5e10
    r = np.linspace(0.05, 60, 60000)
    v_disk = rp.freeman_disk_velocity(r, m_d, r_d)
    # The textbook result: the curve peaks at 2.15 disc scale lengths...
    check("peak radius / R_d", r[np.argmax(v_disk)] / r_d, 2.15, 0.01)
    # ...where V^2 = 0.3877 G M / R_d.
    check("peak V^2 / (G M / R_d)", v_disk.max() ** 2 / (rp.G_KPC * m_d / r_d), 0.3877, 0.002)
    # Far out it must become Keplerian: a disc is a point mass from far enough away.
    far = 20 * r_d
    ratio = rp.freeman_disk_velocity(far, m_d, r_d) / rp.keplerian_velocity(m_d, far)
    check("V_disc / V_Kepler at 20 R_d", ratio, 1.0, 0.01)
    check_true(f"disc approaches Kepler from above ({ratio:.4f} > 1)", ratio > 1.0)

    print("\nDark matter haloes")
    h0 = 73.0
    # V200 = 10 H0 R200: a 200 km/s halo has R200 = 200/730 Mpc.
    check("R200 for V200 = 200 km/s", rp.r200_kpc(200.0, h0), 200.0 / 730.0 * 1e3, 1e-12, "kpc")
    # M200 two ways: from V200^2 R200 / G, and from 200 x critical density x volume.
    r200 = rp.r200_kpc(200.0, h0) * u.kpc
    rho_crit = 3 * (h0 * u.km / u.s / u.Mpc) ** 2 / (8 * np.pi * const.G)
    m200_direct = (200 * rho_crit * 4 / 3 * np.pi * r200**3).to(u.M_sun).value
    check("M200 matches 200 rho_crit volume", rp.m200_msun(200.0, h0), m200_direct, 1e-6, "Msun")
    # An NFW halo's speed at R200 is V200 by definition, for any concentration.
    for c in (5.0, 10.0, 20.0):
        check(f"NFW V(R200) = V200, c = {c:g}", rp.v_nfw(rp.r200_kpc(150.0, h0), 150.0, c, h0),
              150.0, 1e-10, "km/s")
    # ...and it peaks at 2.163 scale radii.
    c = 10.0
    r_s = rp.r200_kpc(150.0, h0) / c
    rr = np.linspace(0.01, 10, 200000) * r_s
    check("NFW peak radius / r_s", rr[np.argmax(rp.v_nfw(rr, 150.0, c, h0))] / r_s, 2.163, 1e-3)
    # Dutton & Maccio (2014): c200 = 10^0.905 = 8.04 at 1e12 Msun/h, and lighter
    # haloes are more concentrated (slope -0.101 per dex).
    check("LCDM c200 at 1e12 Msun/h", rp.nfw_concentration_lcdm(1e12 / rp.DUTTON_MACCIO_H),
          10**0.905, 1e-12)
    check("LCDM c200 at 1e10 Msun/h", rp.nfw_concentration_lcdm(1e10 / rp.DUTTON_MACCIO_H),
          10 ** (0.905 + 0.202), 1e-12)
    # The pseudo-isothermal halo flattens to V_inf, and at r = r_c reaches
    # V_inf * sqrt(1 - pi/4).
    check("ISO V(r_c) / V_inf", rp.v_iso(2.0, 100.0, 2.0) / 100.0, np.sqrt(1 - np.pi / 4), 1e-12)
    check("ISO V(100 r_c) / V_inf", rp.v_iso(200.0, 100.0, 2.0) / 100.0, 0.99218, 1e-4)
    # Central density: V_inf = 100 km/s, r_c = 1 kpc -> 100^2 / (4 pi G) Msun/kpc^3.
    check("ISO central density", rp.iso_central_density(100.0, 1.0),
          1e4 / (4 * np.pi * rp.G_KPC) / 1e9, 1e-12, "Msun/pc^3")

    print("\nThe radial acceleration relation")
    gd = 1.2e-10
    # High accelerations: no mass discrepancy.
    check("g_bar = 100 g_dagger -> g_obs = g_bar", rp.rar_mcgaugh(100 * gd, gd) / (100 * gd), 1.0, 1e-4)
    # Low accelerations: g_obs -> sqrt(g_bar g_dagger).
    gb = 1e-4 * gd
    check("g_bar = 1e-4 g_dagger -> sqrt(g_bar g_dagger)", rp.rar_mcgaugh(gb, gd),
          np.sqrt(gb * gd), 0.006, "m/s^2")
    # At g_bar = g_dagger the boost is 1 / (1 - 1/e) = 1.582.
    check("boost at g_bar = g_dagger", rp.rar_mcgaugh(gd, gd) / gd, 1 / (1 - np.exp(-1)), 1e-12)
    # Its logarithmic slope: 1 where visible matter suffices, 1/2 in the deep
    # regime, and the analytic form must match a numerical derivative.
    check("log slope at g_bar = 1e3 g_dagger", rp.rar_log_slope(1e3 * gd, gd), 1.0, 1e-9)
    check("log slope at g_bar = 1e-6 g_dagger", rp.rar_log_slope(1e-6 * gd, gd), 0.5, 1e-3)
    lg, h = np.log10(gd), 1e-5
    numeric = (np.log10(rp.rar_mcgaugh(10 ** (lg + h), gd))
               - np.log10(rp.rar_mcgaugh(10 ** (lg - h), gd))) / (2 * h)
    check("log slope matches numerical derivative", rp.rar_log_slope(gd, gd), numeric, 1e-6)
    # Flatness falls out of the relation. Far from a mass M, g_bar = G M / r^2;
    # in the low-acceleration regime g_obs = sqrt(g_bar g_dagger), and V^2 = r g_obs
    # then gives V^4 = G M g_dagger -- the same speed at every radius. (This is
    # also the baryonic Tully-Fisher relation, V^4 proportional to M.)
    m = 1e10
    v_flat = ((const.G * m * u.M_sun * gd * u.m / u.s**2) ** 0.25).to(u.km / u.s).value
    for r_kpc in (300.0, 3000.0):
        g_bar = rp.acceleration(rp.keplerian_velocity(m, r_kpc), r_kpc)
        v = np.sqrt(rp.rar_mcgaugh(g_bar, gd) * r_kpc * rp.KPC_M) / 1e3
        check(f"1e10 Msun: V at {r_kpc:g} kpc = (G M g_dagger)^1/4", v, v_flat, 0.005, "km/s")

    print("\nError propagation")
    # A pure velocity error: 5% in V is 10% in g, i.e. 0.0434 dex.
    check("5% velocity error -> dex", rp.log_gobs_error(100.0, 5.0, 90.0, 0.0, 10.0, 0.0),
          0.1 / np.log(10), 1e-9, "dex")
    # Edge-on discs are immune to inclination errors; face-on ones are not.
    check_true("inclination error hurts more at 35 deg than 80 deg",
               rp.log_gobs_error(100, 0, 35, 5, 10, 0) > rp.log_gobs_error(100, 0, 80, 5, 10, 0))
    # g_bar of a pure-gas point knows nothing about the stellar mass-to-light ratio.
    check("all-gas point: only the HI flux error", rp.log_gbar_error(0.0), np.log10(1.1), 1e-12, "dex")
    check("all-star point: only the M/L scatter", rp.log_gbar_error(1.0), 0.11, 1e-12, "dex")

    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: {', '.join(FAILURES)}")
        return 1
    print("All physics checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
