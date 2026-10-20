"""Self-checks for the physics and the measuring machinery, against known answers.

Run with ``python test_physics.py``. No test framework and no network: these
check the formulas and the code paths, each against something known
independently (astropy's cosmology, closed-form integrals, Monte Carlo
simulations, synthetic images).
"""

from __future__ import annotations

import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")

import config
import cosmology
import hosts
import lengths
import radiomaps as rm
import spectra

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float, unit: str = "", *, absolute=False) -> None:
    """``got`` within ``tol`` of ``expected`` (fractional unless ``absolute``)."""
    got = float(got)
    err = abs(got - expected) if absolute or expected == 0 else abs(got - expected) / abs(expected)
    ok = err <= tol
    if not ok:
        FAILURES.append(label)
    how = f"|diff| {err:.2g}" if absolute or expected == 0 else f"{err * 100:.3f}% off"
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<62} {got:>12.6g} vs {expected:>11.6g} {unit:<7} ({how})")


def check_true(label: str, condition: bool) -> None:
    if not condition:
        FAILURES.append(label)
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label}")


def test_cosmology() -> None:
    print("Cosmology (against astropy)")
    from astropy.cosmology import FlatLambdaCDM

    for c in config.COSMOLOGIES:
        ap = FlatLambdaCDM(H0=c.H0, Om0=c.Om0, Tcmb0=0)
        for z in (0.02, 0.5, 0.896, 3.0):
            check(f"{c.name}: D_A(z = {z})", cosmology.angular_diameter_distance(z, c),
                  ap.angular_diameter_distance(z).value, 1e-5, "Mpc")
        check(f"{c.name}: age of the universe", cosmology.age_of_universe(c), ap.age(0).value, 1e-6, "Gyr")
        check(f"{c.name}: lookback time to z = 0.896", cosmology.lookback_time(0.896, c),
              ap.lookback_time(0.896).value, 1e-5, "Gyr")
    c = config.PLANCK18
    ap = FlatLambdaCDM(H0=c.H0, Om0=c.Om0, Tcmb0=0)
    import astropy.units as u

    s, z, a = 0.1, 0.3, -0.7  # 100 mJy at z = 0.3
    expected = (4 * np.pi * ap.luminosity_distance(z) ** 2 * s * u.Jy * (1 + z) ** (-1 - a)).to(u.W / u.Hz).value
    check("radio luminosity of 0.1 Jy at z = 0.3", cosmology.radio_luminosity(s, z, c, a), expected, 1e-5, "W/Hz")
    ang, zm = cosmology.minimum_angle_arcmin(0.7, c)
    check("0.7 Mpc never looks smaller than 1.3' (Oei et al. 2023)", ang, 1.34, 0.01, "arcmin")
    check("Porphyrion: age of the universe at z = 0.896", cosmology.age_at(0.896, c), 6.3, 0.01, "Gyr")


def test_pareto() -> None:
    print("Pareto lengths and projection (Monte Carlo)")
    rng = np.random.default_rng(1)
    xi, lmin = -3.5, 0.3
    L = lmin * rng.random(4_000_000) ** (1 / (xi + 1))
    lp = L * np.sqrt(1 - rng.random(L.size) ** 2)  # isotropic: cos(theta) uniform
    est, err = lengths.pareto_mle(lp, 0.7)
    check("projected lengths keep the intrinsic tail index (Oei et al. Eq. 5)", est, xi, 3 * err / abs(xi), "")
    # deprojection: mean true length of giants seen at l_p
    sel = (lp > 2.95) & (lp < 3.05)
    check("mean true length of l_p = 3 Mpc giants (analytic vs simulated)",
          lengths.deprojection(3.0, xi)["mean"], np.mean(L[sel]), 0.01, "Mpc")
    check("median true length of l_p = 3 Mpc giants", lengths.deprojection(3.0, xi)["median"], np.median(L[sel]),
          0.01, "Mpc")
    check("xi = -3: E[L | L_p] = l_p * Gamma(1.5) Gamma(2.5) / Gamma(2)^2", lengths.deprojection(1.0, -3.0)["mean"],
          0.886227 * 1.329340, 1e-4, "")


def test_selection() -> None:
    print("Selection model")
    for sig in (0.3, 1.2, 2.5):
        lns, g = lengths.sb_detection_table(sig)
        check_true(f"sigma_ref {sig}: detection probability falls monotonically", np.all(np.diff(g) <= 1e-9))
        check(f"sigma_ref {sig}: faint limit -> every giant detected", g[0], 1.0, 1e-3)
        check(f"sigma_ref {sig}: bright limit -> none detected", g[-1], 0.0, 1e-3, absolute=True)
    # without surface-brightness selection the completeness is the volume
    # fraction where the angle exceeds the threshold
    grid = lengths.selection_grid(0.5)
    C = lengths.completeness(grid, 0.0, 1e12, 0.5)
    i = np.searchsorted(grid.lp, 2.0)
    frac = (grid.w * (grid.phi[i] >= config.LENGTH_MODEL["phi_threshold_arcmin"])).sum()
    check("completeness with only a sharp 5' threshold = volume fraction", C[i], frac, 1e-6)
    # recover a known tail index from a simulated catalogue (coarse grid)
    sim = lengths.simulate(1500, -3.5, 1.9, 720.0, 1.2, 0.5, seed=7)
    axes = {"xi": np.arange(-6.0, -1.5 + 1e-9, 0.1), "half_width": np.arange(0, 4.01, 0.5),
            "b_ref": np.arange(100.0, 3001.0, 200.0), "sigma_ref": np.arange(0.4, 2.41, 0.4)}
    post = lengths.grid_posterior(sim["lp"], 0.5, axes=axes)
    s = post.summary("xi")
    check(f"simulated survey with xi = -3.5: recovered {s['mean']:.2f} +- {s['sd']:.2f}", s["mean"], -3.5,
          2 * s["sd"] / 3.5, "")


def test_synchrotron() -> None:
    print("Synchrotron spectra and ages")
    x, s = spectra._jp_grid(-0.6)
    i, j = np.searchsorted(x, 1e-4), np.searchsorted(x, 1e-3)
    check("JP spectrum far below the break follows injection", np.log(s[j] / s[i]) / np.log(x[j] / x[i]), -0.6, 0.02,
          absolute=True)
    x, s = spectra._ci_grid(-0.55)
    i, j = np.searchsorted(x, 30), np.searchsorted(x, 300)
    check("CI spectrum far above the break steepens by 0.5", np.log(s[j] / s[i]) / np.log(x[j] / x[i]), -1.05, 0.01,
          absolute=True)
    z, nb = 0.1, 1000.0
    bmax = spectra.max_age_field(z)
    t0 = spectra.spectral_age_myr(nb, bmax, z)
    check_true("age is largest at B = B_CMB / sqrt(3)",
               t0 > spectra.spectral_age_myr(nb, bmax * 1.05, z) and t0 > spectra.spectral_age_myr(nb, bmax * 0.95, z))
    # Miley's formula: B ~ flux^(2/7) and ~ (volume)^(-2/7)
    b1 = spectra.equipartition_field_ug(1.0, 144, -0.8, 0.1, 600, 200, 300)
    b2 = spectra.equipartition_field_ug(2.0 ** 3.5, 144, -0.8, 0.1, 600, 200, 300)
    check("equipartition field scales as flux^(2/7)", b2 / b1, 2.0, 1e-9)
    # Shulevski et al. 2019, 3C 236 NW lobe: their fluxes and field, their CI age
    a = np.log(0.73 / 5.03) / np.log(1400 / 143)
    age = spectra.spectral_age_myr(spectra.break_from_index(a, -0.55, "CI"), 1.28, 0.1005)
    check("3C 236 NW lobe CI age from Shulevski et al.'s fluxes (theirs: 129 Myr)", age, 129.1, 0.1, "Myr")
    check("black hole mass at sigma = 200 km/s (Kormendy & Ho 2013)", hosts.mbh_from_sigma(200.0), 3.09e8, 1e-9, "Msun")


def test_images() -> None:
    print("Measuring on synthetic images")
    from astropy.wcs import WCS

    n, pix = 400, 4.5
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---SIN", "DEC--SIN"]
    w.wcs.crval = [150.0, 40.0]
    w.wcs.crpix = [n / 2 + 0.5, n / 2 + 0.5]
    w.wcs.cdelt = [-pix / 3600, pix / 3600]
    rng = np.random.default_rng(3)
    sigma = 1e-4
    data = rng.normal(0, sigma, (n, n))
    # a uniform bar 200 pixels (15') long, east-west, at 10 sigma
    data[195:206, 100:300] += 10 * sigma
    img = rm.RadioImage(data, w, 20.0, pix, 144.0, "synthetic")
    check("noise estimate of pure Gaussian noise plus a bar", rm.noise(img), sigma, 0.02, "Jy/beam")
    x0, y0 = img.pixel_to_world(110, 200)
    x1, y1 = img.pixel_to_world(290, 200)
    g = config.Giant("test", "test", 150.0, 40.0, "-", 0.1, 0.0, "s", "-", "-", 30.0,
                     seeds=((float(x0[0]), float(y0[0])), (float(x1[0]), float(y1[0]))))
    ext = rm.measure_extent(img, g)
    check("LAS of a 200-pixel bar (15.0')", ext.las_arcmin, 200 * pix / 60, 0.01, "arcmin")
    # the largest extent of a filled 11 x 200 pixel bar runs corner to corner
    # (between pixel centres 10 rows and 199 columns apart): atan(10/199) = 2.9 deg
    check("its position angle (the bar's diagonal)", abs(ext.pa_deg - 90.0), np.degrees(np.arctan(10 / 199)), 0.3,
          "deg", absolute=True)
    check("its flux density (sum / beam area)", ext.flux_jy, 10 * sigma * 11 * 200 / img.beam_area_pix, 0.02, "Jy")
    # smoothing a point source keeps its flux and lowers its peak by the beam-area ratio
    pt = np.zeros((n, n))
    pt[200, 200] = 1.0
    src = rm.smooth_to(rm.RadioImage(pt, w, 1e-3, pix, 144.0, "pt"), 20.0)  # make it a 20" beam
    big = rm.smooth_to(src, 60.0)
    check("smoothing 20\" -> 60\" conserves flux density", np.nansum(big.data) / big.beam_area_pix,
          np.nansum(src.data) / src.beam_area_pix, 1e-3)
    for rot in (0.0, 27.0, 135.0):
        wr = rm.tangent_wcs(232.38, 60.26, 10.0, 201, rotation_deg=rot)
        c = wr.pixel_to_world_values(100, 100)
        p = wr.pixel_to_world_values(150, 100)
        check(f"rotated grid: +x points to position angle {rot:g}", rm.position_angle(c[0], c[1], p[0], p[1]), rot,
              1e-6, "deg", absolute=True)


def main() -> int:
    test_cosmology()
    test_pareto()
    test_selection()
    test_synchrotron()
    test_images()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} check(s) FAILED: {FAILURES}")
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
