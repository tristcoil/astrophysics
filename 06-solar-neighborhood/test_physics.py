"""Self-checks for the physics module, against cases with known answers.

Run with ``python test_physics.py``. No test framework and no network needed:
these check the formulas, not the data. (The dwarf-sequence checks use the
cached Mamajek table if fetch_data.py has been run, and are skipped otherwise.)

Every star on the map is only as well placed as these conversions. A slip in
one of them -- degrees for radians, a sign in the Galactic rotation, a
magnitude off by five -- would still draw a beautiful, entirely wrong map.
"""

from __future__ import annotations

import numpy as np

import stellar_physics as sp
from config import DWARF_SEQUENCE_PATH

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float, unit: str = "", *, absolute=False) -> None:
    """Assert ``got`` is within ``tol`` of ``expected`` (fractional unless ``absolute``)."""
    err = abs(got - expected) if absolute or expected == 0 else abs(got - expected) / abs(expected)
    ok = err <= tol
    if not ok:
        FAILURES.append(label)
    shown = f"{err:.3g} {unit}".strip() if absolute or expected == 0 else f"{err * 100:.2f}%"
    tol_shown = f"{tol:g} {unit}".strip() if absolute or expected == 0 else f"{tol * 100:g}%"
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<50} {got:>11.5g} vs {expected:>10.5g} {unit:<5} "
          f"(off by {shown}, tol {tol_shown})")


def check_true(label: str, condition: bool, detail: str = "") -> None:
    if not condition:
        FAILURES.append(label)
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label:<50} {detail}")


def main() -> int:
    print("Distances")
    # Proxima Centauri: Gaia DR3 parallax 768.0665 mas -> 4.2465 light years.
    d = float(sp.distance_pc(768.0665))
    check("Proxima: parallax -> distance", d, 1.30197, 1e-4, "pc")
    check("Proxima: distance in light years", d * sp.LY_PER_PC, 4.2465, 1e-4, "ly")
    check("1 pc in light years", sp.LY_PER_PC, 3.26156, 1e-5, "ly")

    print("\nMagnitudes")
    # A star at exactly 10 pc has M = m, by definition; at 100 pc it is 5 mag fainter.
    check("M at 10 pc equals m", float(sp.absolute_magnitude(7.0, 100.0)), 7.0, 1e-9, "mag", absolute=True)
    check("M at 100 pc is m - 5", float(sp.absolute_magnitude(10.0, 10.0)), 5.0, 1e-9, "mag", absolute=True)
    # Sirius: V = -1.46, Hipparcos parallax 379.21 mas -> M_V = +1.43 (textbook ~1.42).
    check("Sirius absolute magnitude", float(sp.absolute_magnitude(-1.46, 379.21)), 1.43, 0.02, "mag", absolute=True)
    # A late-M dwarf (M_G = 15.2) is visible to G = 20.7 out to 10^2.1 = 125.9 pc.
    check("M8 dwarf limiting distance", float(sp.limiting_distance_pc(15.2, 20.7)), 125.89, 1e-3, "pc")
    check("limiting distance when M = m_lim", float(sp.limiting_distance_pc(20.7, 20.7)), 10.0, 1e-9, "pc")

    print("\nGalactic coordinates")
    # Sgr A*, the black hole at the Galactic centre: l = 359.944, b = -0.046 deg.
    # 8178 pc is the GRAVITY Collaboration (2019) distance.
    x, y, z, l, b = sp.galactic_xyz(266.41683, -29.00781, 8178.0)
    check("Sgr A*: galactic longitude", float(l), 359.944, 0.01, "deg", absolute=True)
    check("Sgr A*: galactic latitude", float(b), -0.046, 0.01, "deg", absolute=True)
    check("Sgr A*: lies along +X", float(x), 8178.0, 1e-5, "pc")
    # North Galactic Pole (ICRS 192.85948, +27.12825): straight up the +Z axis.
    x, y, z, l, b = sp.galactic_xyz(192.85948, 27.12825, 100.0)
    check("NGP: galactic latitude", float(b), 90.0, 1e-4, "deg", absolute=True)
    check("NGP: lies along +Z", float(z), 100.0, 1e-6, "pc")
    # l = 90 deg, b = 0 is the direction of Galactic rotation: must be +Y. This is
    # the check that catches a flipped handedness. (ICRS of l=90, b=0 via astropy
    # would be circular, so use the published IAU value: 21h12m01s, +48d19m47s.)
    x, y, z, l, b = sp.galactic_xyz(318.00435, 48.32975, 100.0)
    check("l = 90: galactic longitude", float(l), 90.0, 0.01, "deg", absolute=True)
    check("l = 90: lies along +Y", float(y), 100.0, 1e-4, "pc")

    print("\nGaia colour-contamination flag (Riello et al. 2021)")
    # The published piecewise polynomial must be continuous at its joins; a
    # mistyped coefficient would show up as a jump.
    c_lo, _ = sp.corrected_excess_factor(0.0, 0.5 - 1e-9, 15.0)
    c_hi, _ = sp.corrected_excess_factor(0.0, 0.5 + 1e-9, 15.0)
    check("C* continuous at BP-RP = 0.5", float(c_lo - c_hi), 0.0, 1e-5, absolute=True)
    c_lo, _ = sp.corrected_excess_factor(0.0, 4.0 - 1e-9, 15.0)
    c_hi, _ = sp.corrected_excess_factor(0.0, 4.0 + 1e-9, 15.0)
    check("C* continuous at BP-RP = 4.0", float(c_lo - c_hi), 0.0, 1e-4, absolute=True)
    _, sig = sp.corrected_excess_factor(1.2, 1.0, 20.0)
    check_true("C* scatter grows toward faint stars", float(sp.corrected_excess_factor(1.2, 1.0, 12.0)[1]) < float(sig),
               f"sigma(G=12) < sigma(G=20) = {float(sig):.3f}")

    print("\nProper-motion propagation")
    # 1000 mas/yr in pmra*cos(dec) at dec = 60 deg moves RA by 1/cos(60) = 2 arcsec per year.
    ra, dec = sp.propagate_position(10.0, 60.0, 1000.0, 0.0, 1.0)
    check("RA shift at dec 60 (arcsec/yr)", float((ra - 10.0) * 3600), 2.0, 1e-9, "arcsec")
    ra, dec = sp.propagate_position(10.0, 0.0, 0.0, 10362.394, 24.75)
    check("Barnard's star Dec shift, 1991.25 -> 2016", float((dec - 0.0) * 3600), 256.469, 1e-5, "arcsec")

    print("\nThe eye's colour response (CIE 1931, analytic fit)")
    xb, yb, zb = sp.cie_xyz_bar(np.array([555.0, 600.0, 450.0]))
    check("y-bar(555 nm): peak of eye sensitivity", float(yb[0]), 1.0000, 0.01)
    check("x-bar(600 nm)", float(xb[1]), 1.0622, 0.01)
    check("z-bar(450 nm)", float(zb[2]), 1.7471, 0.03)
    lam = np.arange(360.0, 831.0, 1.0)
    xb, yb, zb = sp.cie_xyz_bar(lam)
    # CIE normalises the three curves to equal area: flat white light is neutral.
    check("equal-energy white: x/y areas", xb.sum() / yb.sum(), 1.0, 0.02)
    check("equal-energy white: z/y areas", zb.sum() / yb.sum(), 1.0, 0.02)

    print("\nBlackbody physics")
    # Wien's law: the Sun (5772 K) peaks at 2.898e-3 / 5772 m = 502 nm.
    lam = np.arange(300.0, 900.0, 0.1)
    peak = float(lam[np.argmax(sp.planck_lambda(lam, 5772.0))])
    check("Wien peak of a 5772 K blackbody", peak, sp.WIEN_B / 5772.0 * 1e9, 0.002, "nm")

    print("\nStar colours")
    rgb = {t: sp.blackbody_srgb(t)[0] for t in (3000.0, 5772.0, 6500.0, 10000.0)}
    check_true("3000 K is orange: R > G > B", rgb[3000][0] > rgb[3000][1] > rgb[3000][2],
               f"rgb = {np.round(rgb[3000], 3)}")
    check_true("5772 K (the Sun) is warm white: R >= G >= B > 0.8",
               rgb[5772][0] >= rgb[5772][1] >= rgb[5772][2] > 0.8, f"rgb = {np.round(rgb[5772], 3)}")
    check_true("6500 K is close to neutral white", float(rgb[6500].min()) > 0.93,
               f"rgb = {np.round(rgb[6500], 3)}")
    check_true("10000 K is blue-white: B > G > R", rgb[10000][2] > rgb[10000][1] > rgb[10000][0],
               f"rgb = {np.round(rgb[10000], 3)}")
    lut = sp.StarColours()
    check("lookup table matches direct calculation", float(np.abs(lut(4321.0) - sp.blackbody_srgb(4321.0)[0]).max()),
          0.0, 0.005, absolute=True)

    print("\nNames")
    for raw, want in [("* alf CMa", "α CMa"), ("* ome02 Tau", "ω2 Tau"), ("NAME Barnard's star", "Barnard's star"),
                      ("Wolf  359", "Wolf 359"), ("* 61 Cyg A", "61 Cyg A"), ("V* V645 Cen", "V645 Cen")]:
        got = sp.tidy_simbad_name(raw)
        check_true(f"tidy {raw!r}", got == want, f"-> {got!r}")

    if DWARF_SEQUENCE_PATH.exists():
        print("\nMamajek dwarf sequence")
        seq = sp.DwarfSequence.from_text(DWARF_SEQUENCE_PATH.read_text())
        # The table's own G2V row: Bp-Rp 0.823, Teff 5770 K, M_G 4.635.
        check("G2V colour -> temperature", float(seq.teff_from_bp_rp(0.823)), 5770.0, 1e-3, "K")
        check("G2V colour -> main-sequence M_G", float(seq.main_sequence_mg(0.823)), 4.635, 0.01, "mag", absolute=True)
        check("Sun-like B-V 0.65 -> temperature", float(seq.teff_from_b_v(0.65)), 5770.0, 1e-3, "K")
        bounds = seq.class_boundaries()
        check("K/M boundary midway between K9V and M0V", bounds["M"][1], 3890.0, 1e-9, "K")
        check("A/B boundary midway between B9.5V and A0V", bounds["A"][1], 10050.0, 1e-9, "K")
        cls = seq.spectral_class(np.array([30000.0, 9000.0, 6500.0, 5770.0, 4500.0, 3200.0, 2000.0]))
        check_true("spectral letters across the sequence", "".join(cls) == "BAFGKML", f"-> {''.join(cls)}")
        check_true("nearest subtype of 5770 K", str(seq.spectral_subtype(np.array([5770.0]))[0]) == "G2V")
    else:
        print("\n(dwarf-sequence checks skipped: run fetch_data.py first)")

    print()
    if FAILURES:
        print(f"{len(FAILURES)} check(s) FAILED: {FAILURES}")
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
