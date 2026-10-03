"""Self-checks for the physics module, against cases with known answers.

Run with ``python test_physics.py``. No test framework and no network: these check
the formulas, not the data (the data-level checks -- our smoothing against SILSO's,
our cycle dates against SILSO's table, our daily spot areas against Hathaway's --
are part of the analysis itself and printed by ``run.py``).

Every headline number here rests on a few small formulas: a smoothing filter, a
rotation-rate conversion, a foreshortening correction. An off-by-one in the filter
or a sign slip in the Carrington frame would produce results that still look
entirely plausible on a plot, so each one is pinned to something whose value is
known independently.
"""

from __future__ import annotations

import sys

import numpy as np

import solar_physics as sp

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float, unit: str = "") -> None:
    """Assert ``got`` is within ``tol`` (absolute) of ``expected``."""
    got = float(got)
    ok = abs(got - expected) <= tol
    if not ok:
        FAILURES.append(label)
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<58} {got:>14.6g} vs {expected:>12.6g} {unit:<8} "
          f"(tol {tol:g})")


def check_true(label: str, condition: bool) -> None:
    if not condition:
        FAILURES.append(label)
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label}")


def main() -> int:
    print("Time")
    # J2000.0 is 2000-01-01 12:00 TT; as a UT calendar date its JD is 2451545.0 exactly.
    check("JD of 2000-01-01 12:00", sp.julian_date(np.datetime64("2000-01-01T12:00"))[()], 2451545.0, 1e-6)
    check("JD round trip", sp.julian_date(sp.datetime_from_jd(2460000.25))[()], 2460000.25, 1e-6)
    # 2000 is a leap year: 1 July is day 183, so 182 days have elapsed of 366.
    check("decimal year of 2000-07-01", sp.decimal_year(np.datetime64("2000-07-01"))[()], 2000 + 182 / 366,
          1e-9)

    print("\nThe 13-month smoothing (SILSO's tapered boxcar)")
    check("weights sum to 1", sp.SMOOTHING_WEIGHTS.sum(), 1.0, 1e-12)
    check("end weights are half the inner ones", sp.SMOOTHING_WEIGHTS[0] / sp.SMOOTHING_WEIGHTS[1], 0.5,
          1e-12)
    flat = sp.smooth_13_month(np.full(40, 7.0))
    check("a constant stays constant", np.nanmax(np.abs(flat - 7.0)), 0.0, 1e-12)
    check_true("first and last six months are left undefined",
               bool(np.isnan(flat[:6]).all() and np.isnan(flat[-6:]).all() and np.isfinite(flat[6:-6]).all()))
    ramp = sp.smooth_13_month(np.arange(40.0))
    check("a straight line passes unchanged (symmetric filter)", np.nanmax(np.abs(ramp - np.arange(40.0))),
          0.0, 1e-12)
    spike = np.zeros(40)
    spike[20] = 24.0
    s = sp.smooth_13_month(spike)
    check("a spike of 24 contributes 24/12 at its own month", s[20], 2.0, 1e-12)
    check("... and half that, 24/24, six months away", s[26], 1.0, 1e-12)
    check("... and nothing seven months away", s[27], 0.0, 1e-12)
    gap = np.full(40, 5.0)
    gap[20] = np.nan
    check_true("a missing month blanks every window that touches it",
               bool(np.isnan(sp.smooth_13_month(gap)[14:27]).all()))

    print("\nFinding cycles")
    t = np.arange(0, 60 * 12) / 12.0
    wave = 100 * (1 - np.cos(2 * np.pi * t / 11.0)) / 2
    ex = sp.find_cycle_extrema(wave)
    check("sinusoid: minima are 11 years apart", np.diff(ex.minima).mean() / 12.0, 11.0, 1 / 12)
    check("sinusoid: each maximum is 5.5 years after its minimum",
          np.mean(ex.maxima[:-1] - ex.minima[:-1]) / 12.0, 5.5, 1 / 12)
    flat_bottom = wave.copy()
    i0 = ex.minima[1]
    flat_bottom[i0 - 4:i0 + 5] = 0.0  # a nine-month run at exactly zero, as in 1810
    ex2 = sp.find_cycle_extrema(flat_bottom)
    check("a flat-bottomed minimum resolves to the middle of the run", ex2.minima[1], i0, 0)

    print("\nThe cycle shape (Hathaway, Wilson & Reichmann 1994)")
    tt = np.linspace(0, 200, 200001)
    b = 50.0
    f = sp.cycle_shape(tt, 0.0, 100.0, b_months=b)
    check("peak height, in units of Rmax", f.max() / 100.0, 1.007, 0.001)
    check("peak time, in units of the width b", tt[np.argmax(f)] / b, 1.09, 0.01)
    check_true("zero before the start time", bool(np.all(sp.cycle_shape(np.array([-5.0, -1.0]), 0.0, 100.0) == 0)))
    # With Rmax = 100 on the old scale (166.7 on the new), (100/Rmax)^(1/4) = 1.
    check("width relation at Rmax(v1) = 100", sp.shape_width_months(100 / 0.6), 27.12 + 25.15, 1e-9,
          "months")

    print("\nThe drift of the sunspot zones")
    check("latitude at the cycle's start", sp.drift_law(0.0, 28.0, 90.0), 28.0, 1e-12, "deg")
    check("one e-folding time later", sp.drift_law(90.0, 28.0, 90.0), 28.0 / np.e, 1e-12, "deg")
    tm = np.linspace(0, 140, 60)
    fit = sp.fit_drift_law(tm, sp.drift_law(tm, 27.0, 95.0))
    check("fit recovers lambda0 from exact data", fit["lambda0_deg"], 27.0, 1e-6, "deg")
    check("fit recovers tau from exact data", fit["tau_months"], 95.0, 1e-4, "months")

    print("\nSolar rotation")
    # Carrington chose a sidereal period of 25.38 days; seen from the orbiting Earth
    # the same rotation takes 27.2753 days -- the length of a Carrington rotation.
    check("synodic period of the Carrington frame", sp.synodic_from_sidereal(25.38), 27.2753, 0.001, "days")
    check("... and back again", sp.sidereal_from_synodic(27.2753), 25.38, 0.001, "days")
    check("Carrington frame rate", 360 / 25.38, 14.1844, 1e-4, "deg/day")
    check("a spot fixed in Carrington longitude turns with the frame", sp.sidereal_rate(0.0), 360 / 25.38,
          1e-12, "deg/day")
    check("a spot drifting +0.37 deg/day turns faster", sp.sidereal_rate(0.37), 360 / 25.38 + 0.37, 1e-12,
          "deg/day")
    check("rotation law at the equator", sp.rotation_law(0.0, 14.551, -2.87), 14.551, 1e-12, "deg/day")
    check("rotation law at 30 deg (sin^2 = 1/4)", sp.rotation_law(30.0, 14.551, -2.87), 14.551 - 2.87 / 4,
          1e-12, "deg/day")
    lat = np.linspace(-35, 35, 71)
    rf = sp.fit_rotation_law(lat, sp.rotation_law(lat, 14.5, -2.8))
    check("rotation fit recovers A", rf["A"], 14.5, 1e-9)
    check("rotation fit recovers B", rf["B"], -2.8, 1e-9)
    check("period at 14.4 deg/day", sp.period_days(14.4), 25.0, 1e-12, "days")
    un = sp.unwrap_degrees(np.array([355.0, 358.0, 1.0, 4.0]))
    check("longitude unwrapping across 360", un[-1] - un[0], 9.0, 1e-9, "deg")

    print("\nFore-shortening")
    check("at disc centre a hemisphere is twice the disc", sp.corrected_area(100.0, 0.0), 50.0, 1e-9)
    # A real RGO record (1976 Jan 12, group 23680): 19 millionths of the disc
    # observed at 0.593 radii from centre, listed as 12 millionths of a hemisphere.
    check("RGO record 1976-01-12, group 23680", sp.corrected_area(19.0, 0.593), 12.0, 0.5)
    # 60 deg from disc centre a spot looks half its true size (cos 60 = 1/2), which
    # exactly cancels the factor of two between hemisphere and disc.
    check("60 deg from disc centre", sp.corrected_area(100.0, np.sin(np.radians(60.0))), 100.0, 1e-6)
    check("Carrington rotation 1 starts at its epoch", sp.carrington_rotation_number(2398167.4), 1.0, 1e-12)
    check("one synodic rotation later", sp.carrington_rotation_number(2398167.4 + 27.2753), 2.0, 1e-9)

    print("\nStatistics")
    rng = np.random.default_rng(3)
    x = rng.normal(size=500)
    c = sp.correlation(x, 2 * x + 1, n_boot=500)
    check("perfect line: r = 1", c["r"], 1.0, 1e-12)
    y = x + rng.normal(size=500)
    c = sp.correlation(x, y, n_boot=2000)
    check("r for y = x + unit noise is 1/sqrt(2)", c["r"], 1 / np.sqrt(2), 0.05)
    check_true("bootstrap 95% interval brackets the estimate", c["r_ci95"][0] < c["r"] < c["r_ci95"][1])
    lf = sp.line_fit(x, 3 + 2 * x + rng.normal(scale=0.5, size=500))
    check("line fit slope", lf["b"], 2.0, 0.1)
    # The Cycle 25 hindcast calibrates on 15 cycles; its 68% prediction interval
    # must cover 68% of new cycles. Checked over many simulated calibrations of that
    # size (any single one can miss by a few percent through its estimated scatter).
    hits = []
    for _ in range(4000):
        xc = rng.normal(size=15)
        fc = sp.line_fit(xc, 3 + 2 * xc + rng.normal(scale=0.5, size=15))
        xn = rng.normal()
        yn = 3 + 2 * xn + rng.normal(scale=0.5)
        _, lo, hi = sp.predict(fc, xn, 0.68)
        hits.append(lo <= yn <= hi)
    check("68% prediction interval (n = 15) covers 68% of new points", np.mean(hits), 0.68, 0.02)

    print()
    if FAILURES:
        print(f"{len(FAILURES)} check(s) FAILED: {', '.join(FAILURES)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
