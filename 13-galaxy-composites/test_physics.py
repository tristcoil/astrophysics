"""Self-checks for the colour physics and the compositing, against known answers.

Run with ``python test_physics.py``. No test framework and no network (the
reference data must be cached): these check the formulas, not the galaxies.

Every claim about "true colour" in this project goes through a chain of
conversions -- spectrum to filter flux, spectrum to CIE XYZ, XYZ to screen
RGB, flux per square arcsecond to candela per square metre. Each link is pinned
here to something whose value is known independently: the definition of the
AB system and of the candela, the chromaticity of the CIE's standard
illuminants, the sRGB standard's own matrix, the Sun's total irradiance.
"""

from __future__ import annotations

import sys

import numpy as np

import color_physics as cp
import composite as cm
import imaging
from config import GALAXIES, Galaxy, target_wcs

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float, unit: str = "", *, absolute=False) -> None:
    """``got`` within ``tol`` of ``expected`` (fractional unless ``absolute``)."""
    got = float(got)
    err = abs(got - expected) if absolute or expected == 0 else abs(got - expected) / abs(expected)
    ok = err <= tol
    if not ok:
        FAILURES.append(label)
    how = f"|diff| {err:.2g}" if absolute or expected == 0 else f"{err * 100:.3f}% off"
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<58} {got:>12.6g} vs {expected:>11.6g} {unit:<7} ({how})")


def check_true(label: str, condition: bool) -> None:
    if not condition:
        FAILURES.append(label)
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label}")


def flat_fnu_spectrum(lam, jy=3631.0):
    """f_lambda of a source with constant f_nu."""
    return jy * cp.JY_CGS * cp.C_ANGSTROM_S / lam**2


def main() -> int:
    lam = np.linspace(1000.0, 30000.0, 29001)

    print("Photometry: the AB system")
    # A source with f_nu = 3631 Jy at every frequency is magnitude 0 in any band.
    for b in "ugriz":
        check(f"flat 3631 Jy spectrum, SDSS {b}", cp.ab_mag(lam, flat_fnu_spectrum(lam), f"sdss_{b}"),
              0.0, 1e-6, "mag", absolute=True)
    # Ten times fainter is exactly 2.5 mag fainter.
    check("10x fainter is +2.5 mag", cp.ab_mag(lam, flat_fnu_spectrum(lam, 363.1), "sdss_r"),
          2.5, 1e-9, "mag", absolute=True)
    check("one nanomaggy in erg/s/cm^2/Hz", cp.nanomaggy_in_cgs(), 3.631e-29, 1e-9)
    check("22.5 mag is one nanomaggy", cp.mag_to_nmgy(22.5), 1.0, 1e-12)

    print("\nThe Sun")
    lam_s, f_s = cp.load_sun()
    # Total solar irradiance, Kopp & Lean 2011 (GRL 38, L01706): 1360.8 +/- 0.5
    # W/m^2. CALSPEC stops at 2.7 microns, and a 5800 K body emits ~2.7% of its
    # energy beyond that; add the tail with the Rayleigh-Jeans slope
    # f_lambda ~ lambda^-4, whose integral from lambda_max is f(lambda_max) lambda_max / 3.
    tail = f_s[-1] * lam_s[-1] / 3.0
    tsi = (np.trapezoid(f_s, lam_s) + tail) * cp.ERG_S_CM2_IN_W_M2
    check("total irradiance of the solar spectrum at 1 AU", tsi, 1360.8, 0.005, "W/m^2")
    # Willmer 2018 (ApJS 236, 47) computed the Sun's absolute magnitudes from the
    # same CALSPEC spectrum through the same filter curves.
    mags = cp.sun_absolute_magnitudes()
    for b, pub in {"g": 5.11, "r": 4.65, "i": 4.53}.items():
        check(f"Sun's absolute AB magnitude, SDSS {b} (Willmer 2018)", mags[b], pub, 0.02, "mag",
              absolute=True)

    print("\nColorimetry")
    cl, cmf = cp.load_cmf()
    # By construction of the 1931 system, the three functions have equal areas,
    # so a flat ("equal-energy") spectrum sits at x = y = 1/3.
    xy_e = cp.xyz_to_xy(cp.spectrum_xyz(cl, np.ones_like(cl)))
    check("equal-energy spectrum: x", xy_e[0], 1 / 3, 2e-4, absolute=True)
    check("equal-energy spectrum: y", xy_e[1], 1 / 3, 2e-4, absolute=True)
    # CIE D65 is defined to have chromaticity (0.31271, 0.32902) for the 2 deg observer.
    lam_d, d65 = cp.load_d65()
    xy_d = cp.xyz_to_xy(cp.spectrum_xyz(lam_d, d65))
    check("D65 spectrum: x (CIE: 0.31271)", xy_d[0], 0.31271, 2e-4, absolute=True)
    check("D65 spectrum: y (CIE: 0.32902)", xy_d[1], 0.32902, 2e-4, absolute=True)
    # The sRGB standard (IEC 61966-2-1) publishes its XYZ -> linear RGB matrix;
    # ours is derived from the primaries and white point, so should agree.
    iec = np.array([[3.2406, -1.5372, -0.4986], [-0.9689, 1.8758, 0.0415], [0.0557, -0.2040, 1.0570]])
    check("sRGB matrix vs IEC 61966-2-1 (max element diff)", np.abs(cp.M_XYZ_TO_SRGB - iec).max(),
          0.0, 6e-4, absolute=True)
    rgb_white = cp.xyz_to_display_rgb(cp.d65_white_xyz())
    check("D65 white -> sRGB (1,1,1): worst channel", np.abs(rgb_white - 1).max(), 0.0, 1e-9, absolute=True)
    v = np.linspace(0, 1, 1001)
    check("sRGB encode/decode round trip", np.abs(cp.srgb_decode(cp.srgb_encode(v)) - v).max(),
          0.0, 1e-12, absolute=True)
    # Chromatic adaptation must carry the source white exactly onto the target white.
    sun_w = cp.sun_white_xyz()
    ad = cp.chromatic_adaptation(sun_w, cp.d65_white_xyz()) @ sun_w
    check("Bradford: sunlight -> D65 white", np.abs(ad - cp.d65_white_xyz()).max(), 0.0, 1e-12, absolute=True)
    _, hexs, _ = cp.swatch(cp.spectrum_xyz(lam_s, f_s), adapt_from=sun_w)
    check_true(f"the Sun adapted to itself renders pure white ({hexs})", hexs == "#FFFFFF")

    print("\nBlackbodies and colour temperature")
    for T in (2856.0, 5000.0, 6500.0, 10000.0):
        cct, duv = cp.correlated_color_temperature(cp.spectrum_xyz(cl, cp.blackbody_flam(cl, T)))
        check(f"blackbody {T:.0f} K: recovered CCT", cct, T, 2e-3, "K")
        check(f"blackbody {T:.0f} K: distance from locus", duv, 0.0, 2e-4, absolute=True)
    # D65's CCT is 6504 K (CIE 15:2004), a little above the locus (Duv ~ +0.003).
    cct, duv = cp.correlated_color_temperature(cp.spectrum_xyz(lam_d, d65))
    check("D65: CCT (CIE: 6504 K)", cct, 6504.0, 2e-3, "K")
    check_true(f"D65 sits on the green side of the locus (Duv = {duv:+.4f})", 0.002 < duv < 0.005)
    lam_w = np.linspace(1000, 100000, 400001)
    peak = lam_w[np.argmax(cp.blackbody_flam(lam_w, 5000.0))]
    check("Wien's law: peak x T", peak * 5000.0, 2.897772e7, 1e-4, "A K")

    print("\nLuminance")
    # The candela's definition: 1 W of 540 THz (555.0 nm) light is 683 lm, and
    # the 1931 ybar is 1.0000 there. A narrow line carrying 1 W/m^2 -> 683 lux.
    line = np.exp(-0.5 * ((cl - 5550.0) / 0.5) ** 2)
    line *= 1.0 / cp.ERG_S_CM2_IN_W_M2 / np.trapezoid(line, cl)  # 1 W/m^2 in erg/s/cm^2/A
    check("1 W/m^2 at 555 nm is 683 lux", cp.illuminance_lux(cl, line), 683.0, 2e-3, "lx")
    # A matte white surface under 1 W/m^2 of 555 nm light: L = E / pi.
    check("Lambertian paper, reflectance 1: L = E/pi", cp.paper_luminance(683.0, 1.0), 683.0 / np.pi, 1e-12)
    # Radiance round trip: a source of flux F spread over 1 arcsec^2 has
    # luminance F_lum / (1 arcsec^2 in sr).
    Y = cp.spectrum_xyz(cl, line)[1]
    check("luminance of that light spread over 1 arcsec^2", cp.luminance_cd_m2(Y),
          683.0 / cp.ARCSEC2_SR, 2e-3, "cd/m^2")
    check("1 arcsec^2 in steradians", cp.ARCSEC2_SR, 2.350443e-11, 1e-6, "sr")

    print("\nCompositing")
    rng = np.random.default_rng(1)
    r, g, b = (rng.gamma(0.8, 20.0, (64, 64)) for _ in range(3))
    Q, sat = 8.0, 250.0
    ours = cm.lupton_rgb(r, g, b, softening=cm.softening_for(Q, sat), saturation=sat)
    # astropy's LuptonAsinhStretch is asinh(Q x / stretch) * 0.1 / asinh(0.1 Q):
    # the same curve as ours when Q_a / stretch = 1 / softening and the
    # normalisations agree, i.e. asinh(0.1 Q_a) = 0.1 Q.
    from astropy.visualization import make_lupton_rgb
    q_a = np.sinh(0.1 * Q) / 0.1
    stretch_a = q_a * cm.softening_for(Q, sat)
    theirs = make_lupton_rgb(r, g, b, minimum=0, stretch=stretch_a, Q=q_a, output_dtype=np.float64)
    check("our Lupton RGB vs astropy.make_lupton_rgb (max diff)", np.abs(ours - theirs).max(), 0.0, 1e-9,
          absolute=True)
    # Colour preservation: channel ratios survive the stretch, even when saturated.
    px = np.array([[[400.0, 200.0, 100.0]]])  # a bright orange pixel, far past saturation
    out = cm.lupton_rgb(px[..., 0], px[..., 1], px[..., 2], softening=1.0, saturation=250.0)[0, 0]
    check("Lupton keeps a saturated pixel's G/R ratio", out[1] / out[0], 0.5, 1e-12)
    check("...and its B/R ratio", out[2] / out[0], 0.25, 1e-12)
    pc = cm.per_channel_asinh(px[..., 0], px[..., 1], px[..., 2], softening=1.0, saturation=250.0)[0, 0]
    check_true(f"per-channel stretch washes it out (G/R = {pc[1] / pc[0]:.2f}, not 0.50)", pc[1] / pc[0] > 0.8)
    # The white balance: a source with these colours must come out grey.
    w = cm.weights_that_whiten({"g": 17.44, "r": 17.0, "i": 16.9})
    f = {bd: cp.mag_to_nmgy(m) for bd, m in {"g": 17.44, "r": 17.0, "i": 16.9}.items()}
    R, G, B = cm.channels({bd: np.array([[f[bd]]]) for bd in "gri"}, w)
    check("white-balanced source: R/G", R[0, 0] / G[0, 0], 1.0, 1e-12)
    check("white-balanced source: B/G", B[0, 0] / G[0, 0], 1.0, 1e-12)
    nc = cm.neutral_color_of(w)
    check("neutral colour recovered from the weights: g-r", nc["g-r"], 0.44, 1e-9, absolute=True)
    # The viewer's recipe: an empty (zero) image maps to the grey set by m.
    zero = {bd: np.zeros((2, 2)) for bd in "gri"}
    vr = cm.viewer_sdss_rgb(zero, scales={"g": (2, 2.5), "r": (1, 1.5), "i": (0, 1.0)}, m=0.02, Q=20.0)
    check("viewer recipe: zero flux -> asinh(Q m)/sqrt(Q)", vr[0, 0, 0], np.arcsinh(0.4) / np.sqrt(20), 1e-6)

    print("\nGeometry")
    fake = Galaxy("test", "Test", "", "NGC 0", 150.0, 30.0, 10.0, 1.0, "", "")
    empty = {bd: np.zeros((600, 600)) for bd in "gri"}
    frame = imaging.Frame(fake, target_wcs(fake), empty, empty)
    ell = imaging.ellipse_mask(frame, 150.0, 30.0, 200.0, 0.5, 30.0)
    check("ellipse mask area = pi a b", ell.sum(), np.pi * 200 * 100, 5e-3, "px")
    north = imaging.ellipse_mask(frame, 150.0, 30.0, 100.0, 0.2, 0.0)
    ys, xs = np.nonzero(north)
    check_true("PA = 0: the major axis runs north-south", np.ptp(ys) > 3 * np.ptp(xs))
    east = imaging.ellipse_mask(frame, 150.0, 30.0, 100.0, 0.2, 90.0)
    ys, xs = np.nonzero(east)
    check_true("PA = 90: the major axis runs east-west", np.ptp(xs) > 3 * np.ptp(ys))
    check("star mask radius, G = 13 (legacypipe formula)", imaging.star_mask_radius_arcsec(13.0),
          1630.0 * 1.396**-13, 1e-12, "arcsec")
    check_true("every gallery frame fits the viewer's 3000 px limit", all(g.size_px <= 3000 for g in GALAXIES))

    print()
    if FAILURES:
        print(f"{len(FAILURES)} check(s) FAILED: {FAILURES}")
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
