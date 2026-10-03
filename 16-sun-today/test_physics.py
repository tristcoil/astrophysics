"""Self-checks for the physics, against cases with known answers.

Run with ``python test_physics.py``. No test framework and no network needed:
these check the formulas and the coordinate conventions, not the data.

The three quantities most easily got subtly wrong are all here: the direction
of Carrington longitude (a sign slip would turn "faster than the Carrington
frame" into "slower"), the factor 2 between density and emission scale heights,
and the GOES flare-class rounding.
"""

from __future__ import annotations

import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")

import config as C  # noqa: E402
import goes  # noqa: E402
import solar_physics as sp  # noqa: E402

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
    print(f"  [{status}] {label:<52} {got:>12.5g} vs {expected:>10.5g} {unit:<9} "
          f"({rel * 100:.2f}% off, tol {tol * 100:.1f}%)")


def check_true(label: str, condition: bool) -> None:
    global N_CHECKS
    N_CHECKS += 1
    if not condition:
        FAILURES.append(label)
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label}")


def main() -> int:
    print("The corona as a hydrostatic atmosphere")
    check("mean mass per particle, He/H = 0.085", sp.mean_particle_mass(0.085),
          (1 + 4 * 0.085) / (2 + 3 * 0.085), 1e-12, "m_H")
    check("mean mass per particle, pure hydrogen", sp.mean_particle_mass(0.0), 0.5, 1e-12, "m_H")
    check("solar surface gravity", sp.surface_gravity(), 274.2, 0.001, "m/s^2")
    h1 = sp.pressure_scale_height_m(1e6) / 1e6
    check("pressure scale height at 1 MK (~50 Mm)", h1, 50.6, 0.01, "Mm")
    check("scale height -> temperature round trip",
          sp.temperature_from_scale_height(sp.pressure_scale_height_m(1.7e6)), 1.7e6, 1e-12, "K")
    check("density at the base is 1", float(sp.hydrostatic_density(1.0, 1e6)), 1.0, 1e-12)
    h = 1e-3  # plane-parallel limit: ln n = -h/H for h << R
    lam = C.R_SUN_M / sp.pressure_scale_height_m(1e6)
    check("near the base, ln n falls as -h/H",
          -np.log(float(sp.hydrostatic_density(1 + h, 1e6))) / (lam * h), 1.0, 0.002)
    # Line-of-sight integral: along the ray 1 - 1/r ~ 1 - 1/p + s^2/(2p^3), so for
    # H << R the integral is a Gaussian, I(p) ~ n^2(p) sqrt(pi p^3 H/R), and the
    # log-slope is -2 (R/H) / p^2 + 3/(2p). Check the numerical integral against it.
    p = np.array([1.04, 1.05, 1.06])
    i = sp.offlimb_intensity(p, 1e6)
    slope = (np.log(i[2]) - np.log(i[0])) / (p[2] - p[0])
    check("off-limb log-slope vs asymptotic formula at 1.05 R", slope,
          -2 * lam / 1.05**2 + 3 / (2 * 1.05), 0.01)
    hh = np.linspace(12, 80, 40)
    truth = sp.offlimb_intensity(1 + hh * 1e6 / C.R_SUN_M, 1.4e6)
    rng = np.random.default_rng(1)
    noisy = truth * 3.7e4 * np.exp(rng.normal(0, 0.03, truth.size))
    t_fit, rms, _ = sp.fit_offlimb_temperature(hh, noisy)
    check("hydrostatic fit recovers 1.4 MK through 3% noise", t_fit, 1.4e6, 0.03, "K")
    check("exponential e-folding height of a pure exponential",
          sp.emission_scale_height_mm(hh, np.exp(-hh / 23.0)), 23.0, 1e-9, "Mm")
    t_simple = sp.temperature_from_scale_height(2 * sp.emission_scale_height_mm(hh, truth) * 1e6)
    check_true("plain exponential fit overestimates T of a spherical atmosphere "
               f"({t_simple / 1e6:.2f} vs 1.40 MK)", t_simple > 1.4e6)

    print("\nRotation")
    lat = np.arange(-45, 46, 5.0)
    law = C.PUBLISHED_ROTATION["komm1993"]
    omega = sp.rotation_law(lat, law["A"], law["B"], law["C"])
    noisy = omega + rng.normal(0, 0.03, lat.size)
    coef, cov, chi2, dof = sp.fit_rotation_law(lat, noisy, np.full(lat.size, 0.03))
    check("A recovered from a noisy Komm et al. profile", coef[0], law["A"], 0.003, "deg/day")
    check_true(f"chi^2/dof sensible ({chi2:.1f}/{dof})", 0.3 < chi2 / dof < 2.5)
    check("equator: sidereal 25.38 d at the Carrington rate",
          360 / C.CARRINGTON_RATE_DEG_PER_DAY, 25.38, 0.001, "days")
    w_peri = sp.earth_orbital_rate_deg_per_day("2026-01-03T12:00")
    w_aph = sp.earth_orbital_rate_deg_per_day("2026-07-06T12:00")
    check("Earth's orbital rate at perihelion (Kepler II)", w_peri, 1.0193, 0.002, "deg/day")
    check("Earth's orbital rate at aphelion", w_aph, 0.9534, 0.002, "deg/day")
    check("Carrington synodic period (27.2753 d)",
          sp.synodic_period_days(C.CARRINGTON_RATE_DEG_PER_DAY, 360 / 365.25636), 27.2753, 0.0005,
          "days")
    # The sign convention: a point fixed in Carrington coordinates must drift WEST
    # (toward larger Stonyhurst longitude) by 14.18 - 0.99 = 13.2 deg/day.
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from sunpy.coordinates import frames

    t0 = "2026-09-28T12:00"
    obs = frames.HeliographicStonyhurst(0 * u.deg, 7 * u.deg, 1 * u.AU, obstime=t0)
    fixed = SkyCoord(100 * u.deg, 10 * u.deg, frame=frames.HeliographicCarrington,
                     obstime=t0, observer="earth")
    lon0 = fixed.transform_to(frames.HeliographicStonyhurst(obstime=t0)).lon.to_value(u.deg)
    later = SkyCoord(100 * u.deg, 10 * u.deg, frame=frames.HeliographicCarrington,
                     obstime="2026-09-29T12:00", observer="earth")
    lon1 = later.transform_to(frames.HeliographicStonyhurst(obstime="2026-09-29T12:00")).lon.to_value(u.deg)
    del obs
    check("fixed Carrington point drifts west 13.2 deg/day",
          (lon1 - lon0 + 180) % 360 - 180,
          C.CARRINGTON_RATE_DEG_PER_DAY - sp.earth_orbital_rate_deg_per_day("2026-09-29"), 0.005,
          "deg")

    print("\nThe photosphere")
    check("Wien peak of B_lambda at T_eff (500 nm)",
          _wien_peak(C.T_EFF_K), 2.897771955e-3 / C.T_EFF_K, 1e-3, "m")
    b = sp.planck_lambda(617.33e-9, 6100.0)
    check("Planck radiance -> temperature round trip", float(sp.planck_temperature(617.33e-9, b)),
          6100.0, 1e-10, "K")
    u_, v_ = sp.allen_limb_coefficients(6173.3)
    check("Allen's law at mu = 1 is 1", float(sp.allen_limb_darkening(1.0, 6173.3)), 1.0, 1e-12)
    check_true(f"Allen's u, v at 6173 A plausible (u = {u_:.3f}, v = {v_:.3f})",
               0.7 < u_ < 1.0 and -0.4 < v_ < 0.0)
    limb = [float(sp.allen_limb_darkening(0.0, wl)) for wl in (4500, 6173.3, 8000, 12000)]
    check_true("limb/centre contrast weakens toward the red "
               f"({', '.join(f'{x:.2f}' for x in limb)})", all(np.diff(limb) > 0))
    check("Eddington grey law at the limb", float(sp.eddington_limb_darkening(0.0)), 0.4, 1e-12)
    check("grey atmosphere: T(tau = 2/3) = T_eff", float(sp.grey_temperature(2 / 3)), C.T_EFF_K,
          1e-12, "K")
    check("grey atmosphere: surface T = T_eff / 2^(1/4)", float(sp.grey_temperature(0.0)),
          C.T_EFF_K * 0.5**0.25, 1e-12, "K")
    mu = np.linspace(0.1, 1, 10)
    wl = 617.33e-9
    ratio = sp.planck_lambda(wl, sp.grey_temperature(mu)) / sp.planck_lambda(wl, sp.grey_temperature(1))
    anchor = float(sp.planck_lambda(wl, C.T_EFF_K) / sp.planck_lambda(wl, sp.grey_temperature(1)))
    back = sp.eddington_barbier_temperature(mu, ratio, wl, ratio_at_anchor=anchor)
    check("Eddington-Barbier inversion returns the input T(tau)",
          float(np.max(np.abs(back / sp.grey_temperature(mu) - 1))), 0.0, 1e-9)

    print("\nFlares")
    for flux, cls in ((1.3762636e-6, "C1.3"), (2.3e-6, "C2.3"), (9.99e-5, "M9.9"),
                      (1.2e-4, "X1.2"), (2.4969406e-6, "C2.4"), (5.5e-8, "A5.5")):
        check_true(f"GOES class of {flux:.4g} W/m^2 is {cls} (got {goes.flare_class(flux)})",
                   goes.flare_class(flux) == cls)
    for code, idx in ((5, 0), (6, 1), (9, 2), (10, 3), (7, 1)):
        check_true(f"NCEI primary-channel byte {code} -> response row offset {idx}",
                   goes.decode_primary_chan(code) == idx)
    check_true("GOES-18 A1+B1 response row is 23", goes.response_index(18, 0) == 23)
    check("thermal energy 3kT sqrt(EM V)", float(sp.thermal_energy_erg(1e7, 1e48, 1e27)),
          3 * 1.380649e-16 * 1e7 * np.sqrt(1e75), 1e-6, "erg")
    check("radiative cooling time, 10 MK at 1e10 cm^-3",
          float(sp.radiative_cooling_time_s(1e7, 1e10)),
          3 * 1.380649e-16 * 1e7 / (1e10 * 1.2e-19 * 1e7**-0.5), 1e-6, "s")

    print("\nGeometry")
    import sdo

    check("mu at r = 0.6 is 0.8", float(sdo.mu_from_radius(0.6)), 0.8, 1e-12)
    check("a 1.45 Mm pixel at disc centre", float(sp.msh_per_pixel(1.45, 1.0)),
          1.45**2 / (2 * np.pi * 695.7**2) * 1e6, 1e-9, "MSH")

    print()
    if FAILURES:
        print(f"{len(FAILURES)} of {N_CHECKS} CHECKS FAILED: {', '.join(FAILURES)}")
        return 1
    print(f"All {N_CHECKS} physics checks passed.")
    return 0


def _wien_peak(t: float) -> float:
    lam = np.linspace(200e-9, 2000e-9, 200001)
    return float(lam[np.argmax(sp.planck_lambda(lam, t))])


if __name__ == "__main__":
    sys.exit(main())
