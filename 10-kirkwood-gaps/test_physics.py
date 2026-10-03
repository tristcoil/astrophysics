"""Self-checks against cases with known answers.

Run with ``python test_physics.py``. No test framework, and no network: the
two checks that use real data read it from the cache.

Every gap location in this project comes from Kepler's third law and Jupiter's
period, and every position on the top-down map comes from a chain of angle
rotations. A slip in either would still produce plausible-looking plots, which
is exactly why they are checked here against numbers known in advance.
"""

from __future__ import annotations

import json
import sys

import numpy as np

import mpcorb
import orbit_physics as op
from config import HORIZONS_PATH, MPCORB_PATH

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float, unit: str = "") -> None:
    """Pass if ``got`` is within ``tol`` of ``expected`` (fractional when
    ``expected`` is nonzero, absolute otherwise)."""
    err = abs(got - expected) / abs(expected) if expected else abs(got)
    ok = err <= tol
    if not ok:
        FAILURES.append(label)
    shown = f"{err * 100:.3g}% off" if expected else f"off by {err:.2g}"
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<44} {got:>14.7g} vs {expected:>12.7g} "
          f"{unit:<6} ({shown})")


def check_equal(label: str, got, expected) -> None:
    ok = got == expected
    if not ok:
        FAILURES.append(label)
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<44} {got!r} vs {expected!r}")


def test_packed_formats() -> None:
    print("MPC packed dates (examples from the MPC's own format page)")
    for packed, expected in [("J9611", (1996, 1, 1.0)), ("J961A", (1996, 1, 10.0)),
                             ("J969U", (1996, 9, 30.0)), ("J96A1", (1996, 10, 1.0)),
                             ("K01AM", (2001, 10, 22.0)), ("K2669", (2026, 6, 9.0))]:
        check_equal(f"unpack {packed}", mpcorb.unpack_date(packed), expected)
    y, m, d = mpcorb.unpack_date("K01AM138303")
    check("fractional day K01AM138303", d, 22.138303, 1e-12, "d")

    print("\nJulian dates")
    # J2000.0 is 2000 January 1.5 = JD 2451545.0 by definition.
    check("JD of 2000-01-01.5", mpcorb.calendar_to_jd(2000, 1, 1.5), 2451545.0, 1e-12)
    check("JD of 1996-01-01.0", mpcorb.packed_date_to_jd("J9611"), 2450083.5, 1e-12)
    # MPCORB's standard epoch, 2026 June 9.0 TT, is 112 days before the
    # 2026 September 29.0 snapshot, which Horizons reports as JD 2461312.5.
    check("JD of K2669 (2026-06-09.0)", mpcorb.packed_date_to_jd("K2669"), 2461200.5, 1e-12)

    print("\nMPC packed numbers")
    for packed, expected in [("00001", 1), ("99999", 99999), ("A0000", 100000),
                             ("A0345", 100345), ("Z9999", 359999), ("a0000", 360000),
                             ("z9999", 619999), ("~0000", 620000), ("~000z", 620061),
                             ("~AZaz", 3140113), ("~zzzz", 15396335),
                             ("K07Tf8A", 0)]:
        check_equal(f"unpack {packed}", mpcorb.unpack_number(packed), expected)


def test_kepler() -> None:
    print("\nKepler's third law")
    # A massless body with Earth's sidereal year sits at 1 AU (the AU and the
    # Gaussian constant are defined that way, to the precision that matters).
    check("a(365.256363 d)", float(op.semi_major_axis(365.256363)), 1.0, 1e-5, "AU")
    check("a(Mars, 686.980 d)", float(op.semi_major_axis(686.980)), 1.52368, 1e-4, "AU")
    # Jupiter's mean distance (Standish's J2000 mean elements: 5.2029 AU).
    check("Jupiter a from its mean period", op.jupiter_semi_major_axis(), 5.2029, 2e-4, "AU")

    print("\nThe MPC's convention for mean motion")
    # MPCORB lists Ceres with a = 2.7655526 AU and n = 0.21430445 deg/day
    # (the file retrieved for this project). Matching it to 8 digits confirms
    # the MPC uses exactly k/a^1.5 for a massless asteroid, as we do.
    check("n(Ceres) from a", float(op.mean_motion(2.7655526)), 0.21430445, 1e-7, "deg/d")

    print("\nResonance locations (as tabulated on Wikipedia's 'Kirkwood gap' page)")
    # Published values differ in the 3rd decimal depending on which value of
    # Jupiter's semi-major axis is used, so the tolerance is 0.1%.
    for (p, q), lit in [((4, 1), 2.065), ((3, 1), 2.502), ((5, 2), 2.825),
                        ((7, 3), 2.958), ((2, 1), 3.279), ((3, 2), 3.972),
                        ((4, 3), 4.296)]:
        check(f"{p}:{q} resonance", op.resonance_a(p, q), lit, 1e-3, "AU")
    # A three-body resonance must satisfy its defining condition exactly:
    # 5 n_J - 2 n_S - 2 n = 0 for the "5J-2S-2" resonance that cuts through
    # the Veritas family near 3.17 AU (Nesvorny & Morbidelli 1998).
    a3 = op.three_body_resonance_a(5, -2, -2)
    n_j, n_s = 360 / op.JUPITER_PERIOD_DAYS, 360 / op.SATURN_PERIOD_DAYS
    check("5 n_J - 2 n_S - 2 n(a) at 5J-2S-2", 5 * n_j - 2 * n_s - 2 * float(op.mean_motion(a3)),
          0.0, 1e-12, "deg/d")
    check("5J-2S-2 location", a3, 3.17, 0.01, "AU")
    check("period ratio at the 3:1", float(op.period_ratio(op.resonance_a(3, 1))), 3.0, 1e-12)


def test_orbits() -> None:
    print("\nKepler's equation")
    rng = np.random.default_rng(1)
    M = rng.uniform(0, 2 * np.pi, 100_000)
    e = rng.uniform(0, 0.99, 100_000)
    E = op.solve_kepler(M, e)
    resid = np.abs(np.mod(E - e * np.sin(E) - M + np.pi, 2 * np.pi) - np.pi).max()
    check("max |E - e sinE - M|, 1e5 random orbits", resid, 0.0, 1e-10, "rad")

    print("\nElements -> position, simple geometry")
    x, y, z = op.elements_to_xyz(2.0, 0.5, 0.0, 0.0, 0.0, 0.0)
    check("perihelion distance a(1-e)", float(x), 1.0, 1e-12, "AU")
    x, y, z = op.elements_to_xyz(2.0, 0.5, 0.0, 0.0, 0.0, 180.0)
    check("aphelion distance a(1+e)", float(-x), 3.0, 1e-12, "AU")
    x, y, z = op.elements_to_xyz(3.0, 0.0, 30.0, 0.0, 90.0, 0.0)
    check("circular orbit, 30 deg tilt: max height", float(z), 1.5, 1e-12, "AU")

    print("\nElements -> position, against JPL Horizons (real planets)")
    # Horizons returns both osculating elements and the position vector at
    # the same instant. Turning the former into the latter with our code must
    # reproduce Horizons' own vector: this checks every angle convention.
    if HORIZONS_PATH.exists():
        bodies = json.loads(HORIZONS_PATH.read_text())["bodies"]
        for name in ("Mars", "Jupiter", "Saturn"):
            b = bodies[name]
            xyz = op.elements_to_xyz(b["a_au"], b["e"], b["incl_deg"], b["node_deg"],
                                     b["peri_deg"], b["mean_anomaly_deg"])
            err = float(np.linalg.norm(np.array(xyz, dtype=float) - np.array(b["xyz_au"])))
            check(f"{name}: |our xyz - Horizons xyz|", err, 0.0, 1e-6, "AU")
    else:
        print("  (skipped: run fetch_data.py first)")

    print("\nResonant angle")
    # A Hilda at perihelion in conjunction with Jupiter: sigma must be 0.
    check("3:2 angle at perihelic conjunction", float(op.resonant_angle(3, 2, 40.0, 40.0, 40.0)),
          0.0, 1e-12, "deg")
    # A Trojan 60 deg ahead of Jupiter sits at L4: longitude from Jupiter +60.
    check("L4 Trojan: longitude from Jupiter", float(op.longitude_from_jupiter(100.0, 40.0)),
          60.0, 1e-12, "deg")
    check("wrap-around: 10 deg vs 350 deg", float(op.longitude_from_jupiter(10.0, 350.0)),
          20.0, 1e-12, "deg")

    print("\nSizes")
    # Ceres: H = 3.34 (MPCORB); geometric albedo 0.09 and mean diameter
    # 939 km, both measured by the Dawn spacecraft.
    check("Ceres diameter from H", float(op.diameter_km(3.34, 0.09)), 939.4, 0.03, "km")


def test_mpcorb_file() -> None:
    print("\nParsing the real MPCORB file")
    if not MPCORB_PATH.exists():
        print("  (skipped: run fetch_data.py first)")
        return
    d = mpcorb.load_mpcorb(quiet=True)
    first = {k: d[k][:4] for k in ("number", "a", "e", "incl", "H")}
    check_equal("first four numbers", first["number"].tolist(), [1, 2, 3, 4])
    # Well-known values: Ceres a = 2.766, Vesta i = 7.14 deg, Pallas e = 0.231.
    check("Ceres a", float(first["a"][0]), 2.766, 2e-3, "AU")
    check("Pallas e", float(first["e"][1]), 0.231, 1e-2)
    check("Vesta i", float(first["incl"][3]), 7.14, 1e-2, "deg")
    # Every numbered object must decode to a unique number, 1..N.
    num = d["number"][d["section"] == 0]
    check_equal("numbered objects are exactly 1..N",
                bool(np.array_equal(np.sort(num), np.arange(1, len(num) + 1))), True)
    # The file's own mean motions must agree with its semi-major axes, to the
    # rounding of the printed digits (n to 1e-8 deg/day, a to 1e-7 AU).
    # The rounding bound is 5e-9 in n plus the effect of 5e-8 AU in a,
    # which is dn = 1.5 (n/a) da.
    sel = (d["a"] > 1.5) & (d["a"] < 6)
    a, n = d["a"][sel], d["n"][sel]
    bound = 5e-9 + 1.5 * n / a * 5e-8
    excess = np.abs(op.mean_motion(a) - n) / bound
    check_equal(f"objects outside the rounding bound (of {len(a):,})",
                int((excess > 1.0).sum()), 0)


def main() -> int:
    test_packed_formats()
    test_kepler()
    test_orbits()
    test_mpcorb_file()
    try:
        import test_algorithms
    except ImportError:
        pass
    else:
        FAILURES.extend(test_algorithms.run())
    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILED: {FAILURES}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
