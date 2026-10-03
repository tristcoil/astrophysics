"""Run the whole project end to end.

    python run.py

Fetches whatever is not already cached (3.9 GB the first time), measures the
terminal velocities, builds the rotation curve and the outer-disc maps, writes
every figure, and prints a summary that grades our numbers against the
published ones. Once the cache is warm this runs offline, in about 6 minutes.

    python run.py --refresh-outer

rebuilds the outer-disc maps (the slow step) even if they are cached.
"""

from __future__ import annotations

import argparse
import sys
import warnings

warnings.filterwarnings("ignore")

import analysis
import figures
from fetch_data import fetch_all, missing_tiles


def _f(x, fmt="{:.1f}"):
    return "n/a" if x is None else fmt.format(x)


def summary(r: dict) -> None:
    t = r["terminal"]
    print("\nTerminal velocities (the edge of the hydrogen emission)")
    q = t["q1_vs_mgd16"]
    print(f"  first quadrant vs McClure-Griffiths & Dickey 2016 (VLA, 1'), {q['n']} longitudes:")
    print(f"    median offset {q['median']:+.2f} km/s, rms {q['rms']:.2f} km/s, "
          f"worst {q['max_abs']:.1f} km/s")
    qs = t["q1_vs_mgd16_smoothed"]
    print(f"    after smoothing theirs to our 16' beam: rms {qs['rms']:.2f} km/s")
    q4 = t["q4_vs_mgd07_table"]
    print(f"  fourth quadrant vs the 10 values printed by McClure-Griffiths & Dickey 2007: "
          f"mean offset {q4['mean']:+.2f} km/s")
    for quad in ("Q1", "Q4"):
        lf = t["linear_fit"][quad]
        print(f"  {quad} fit |V_t| = v_circ (1-|sin l|) + v': v_circ {lf['v_circ']:.1f} "
              f"(published {lf['published'][0]}), v' {lf['v_prime']:.1f} "
              f"(published {lf['published'][1]})")
    a = t["AR0"]
    print(f"  Oort A x R0 (fourth quadrant): {a['Q4_0.35-0.95']['AR0']:.1f} +- "
          f"{a['Q4_0.35-0.95']['err']:.1f} km/s over |sin l| 0.35-0.95 (published 89 +- 4); "
          f"{a['Q4_0.80-0.95']['AR0']:.0f} over 0.80-0.95 (published 110 +- 2; ours swings "
          f"{a['Q4_0.80-0.90']['AR0']:.0f}-{a['Q4_0.85-0.95']['AR0']:.0f} as the range shifts)")
    ns = t["ns_asymmetry"]
    print(f"  north - south: {ns['at_25']:+.1f} / {ns['at_40']:+.1f} / {ns['at_65']:+.1f} km/s at "
          f"|l| = 25 / 40 / 65 (published -15 / +10 / -10); mean at 55-67: "
          f"{ns['mean_55_67']:+.2f} (published -5.2); LSR-outflow fit Pi0 = "
          f"{ns['lsr_outflow_fit_Pi0']:.1f} km/s leaves rms {ns['rms_after_Pi0']:.1f} of "
          f"{ns['rms_before']:.1f}")
    ms = t["method_sensitivity"]
    print(f"  method: 12 variants (threshold, window, latitude statistic) move a longitude by "
          f"{ms['Q1']['median_spread']:.2f} (Q1) / {ms['Q4']['median_spread']:.2f} (Q4) km/s "
          "(median spread)")

    rot = r["rotation"]
    print("\nRotation curve from the tangent points (R0 = 8.15 kpc, Reid et al. 2019 constants)")
    for quad in ("Q1", "Q4"):
        b = rot["binned"][quad]
        pts = [(R, th) for R, th in zip(b["R"], b["theta"]) if th is not None]
        sample = [(R, th) for R, th in pts if abs(R - round(R)) < 0.13 and round(R) in (3, 4, 5, 6, 7)]
        print(f"  {quad}: " + ", ".join(f"R={R:.1f}: {th:.0f}" for R, th in sample) + " km/s")
    m = rot["mgd_style_fit"]
    print(f"  with the old IAU constants, Theta/Theta0 = {m['k0']:.3f} R/R0 + {m['k1']:.3f} "
          f"(McClure-Griffiths & Dickey 2016: 0.171, 0.889; 2007: 0.186, 0.887)")

    mz = r["masers"]
    print(f"\nMasers: {mz['n_selected']} of {mz['n_total']} pass Reid et al.'s cuts; their orbits "
          f"lag the rotation curve by {-mz['vrot_minus_reid_urc']['median']:.1f} km/s (median) and "
          f"fall inward at {-mz['vrad_median_outward']:.1f} km/s (published: 4.3 +- 5.6, "
          "6.0 +- 1.4)")
    p = mz["perseus_kinematic_over_parallax"]
    print(f"  Perseus arm: kinematic distance / parallax distance = {p['median']:.2f} (median of "
          f"{p['n']}), range {p['range'][0]:.2f}-{p['range'][1]:.2f}")

    c = r["context"]
    mm = c["enclosed_mass_Msun"]
    print(f"\nMass: {mm['4'] / 1e10:.1f}e10 Msun inside 4 kpc, {mm['R0_reid_theta0'] / 1e11:.2f}e11 "
          f"inside the Sun's orbit, {mm['eilers_outermost']['M'] / 1e11:.1f}e11 inside "
          f"{mm['eilers_outermost']['R']:.0f} kpc (spherical V^2 R / G)")
    print(f"  Among SPARC's {c['sparc_n_with_vflat']} galaxies with a flat speed, "
          f"{100 * c['sparc_fraction_faster_than_mw']:.0f}% rotate faster than the Milky Way")

    o = r["outer"]
    print(f"\nOuter disc: HI surface density scale length "
          f"{_f(o['default']['sigma_scale_length_kpc'], '{:.2f}')} kpc (12-25 kpc)")
    for ring in o["lbh06"]["rings"]:
        if ring["R"] in (15.5, 20.5, 25.5):
            print(f"  warp at R = {ring['R']} (LBH06 constants): W0 {ring['W0']:.0f}, "
                  f"W1 {ring['W1']:.0f}, W2 {ring['W2']:.0f} pc")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh-outer", action="store_true",
                        help="rebuild the outer-disc maps even if cached")
    args = parser.parse_args(argv)
    if missing_tiles():
        print("Fetching data (first run only)")
    try:
        fetch_all(quiet=True)
    except RuntimeError as err:
        print(err, file=sys.stderr)
        return 1
    print("Measuring (terminal velocities, rotation curve, masers, outer disc)")
    results = analysis.run(refresh_outer=args.refresh_outer)
    print("Drawing figures")
    for path in figures.make_all():
        print(f"  {path}")
    summary(results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
