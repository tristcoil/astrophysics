"""Run the whole project end to end.

    python run.py

Downloads whatever is not already cached, analyses every target, writes the
figures, and prints a summary comparing our measurements to the literature.
Once the cache is warm this runs offline.
"""

from __future__ import annotations

import argparse
import sys
import warnings

warnings.filterwarnings("ignore")

import analysis
import figures
from config import FIGURE_DIR, TARGETS, target_by_name
from fetch_data import fetch_light_curves, fetch_stellar_parameters, load_stellar_parameters


def _fmt(value, spec=".4f", dash="—"):
    return format(value, spec) if value is not None else dash


def print_summary(results) -> None:
    print("\n" + "=" * 78)
    print("RESULTS — our measurements vs published values")
    print("=" * 78)

    for r in results:
        print(f"\n{r.planet}   ({r.kic}, {r.difficulty} case)")
        print(f"  {r.n_cadences:,} cadences over {r.baseline_days:.0f} days, "
              f"{r.n_transits_observed} transits")
        print(f"  detection: SDE = {r.sde:.1f}, folded SNR = {_fmt(r.snr, '.0f')}")

        rows = [
            ("Orbital period (d)", r.comparison.get("period"), ".6f"),
            ("Transit depth (ppm)", r.comparison.get("depth_ppm"), ".1f"),
            ("Duration T14 (h)", r.comparison.get("duration_hours"), ".3f"),
            ("Radius, raw depth (R_e)", r.comparison.get("radius_rearth"), ".3f"),
            ("Radius, limb-dark corr.", r.comparison.get("radius_rearth_ld_corrected"), ".3f"),
        ]
        print(f"    {'quantity':<26}{'ours':>14}{'published':>14}{'diff':>10}")
        for label, c, spec in rows:
            if not c:
                continue
            print(f"    {label:<26}{format(c['ours'], spec):>14}"
                  f"{format(c['published'], spec):>14}"
                  f"{c['percent_diff']:>9.2f}%")

        v = r.validation
        oe = v.get("odd_even", {})
        se = v.get("secondary_eclipse", {})
        trap = v.get("trapezoid", {})
        print("    validation:")
        print(f"      [INFO] SDE = {r.sde:.1f}  (noise alone can reach ~7-9 on "
              f"Kepler's gappy sampling, so only a large margin counts)")
        pl = v.get("plausibility", {})
        if pl.get("ok"):
            ok = "PASS" if pl["passes"] else "FLAG"
            print(f"      [{ok}] size gate: implies {pl['implied_radius_rjup']:.3f} R_jup "
                  f"— {pl['verdict']}")
        ds = v.get("depth_stability", {})
        if ds.get("ok"):
            ok = "PASS" if ds["consistent"] else "FLAG"
            print(f"      [{ok}] depth stable over the mission: chunk scatter "
                  f"{ds['observed_scatter_ppm']:.1f} ppm vs {ds['expected_scatter_ppm']:.1f} "
                  f"expected (ratio {ds['ratio']:.2f})")
        if oe.get("significance_sigma") is not None:
            ok = "PASS" if oe["consistent"] else "FLAG"
            frac = oe.get("fractional_difference")
            extra = f", {frac:.2%} fractionally" if frac is not None else ""
            print(f"      [{ok}] odd/even depths differ by "
                  f"{oe['significance_sigma']:.1f} sigma{extra}")
        if se.get("significance_sigma") is not None:
            ok = "FLAG" if se.get("binary_like") else "PASS"
            print(f"      [{ok}] secondary eclipse {se['depth_ppm']:.1f} ppm "
                  f"({se['significance_sigma']:.1f} sigma)")
            if se.get("interpretation"):
                print(f"             {se['interpretation']}")
        if trap.get("ok"):
            print(f"      [INFO] trapezoid fit: T14 = {trap['t14_hours']:.3f} h "
                  f"as observed, flat T23 = {trap['t23_hours']:.3f} h")
            if trap.get("smearing_correction_valid"):
                print(f"             correcting for the {trap['cadence_hours'] * 60:.1f} min "
                      f"exposure: T14 = {trap['t14_corrected_hours']:.3f} h, "
                      f"ingress {trap['ingress_corrected_minutes']:.1f} min")
            print(f"             residual RMS {trap['rms_residual_ppm']:.0f} ppm")

        d = r.derived
        if d:
            ld = d.get("limb_darkening", {})
            print("    derived:")
            br = r.comparison.get("radius_bracket_rearth")
            if br:
                inside = "CONTAINS" if br["published_inside"] else "does NOT contain"
                print(f"      Rp between {br['low']:.3f} and {br['high']:.3f} R_earth "
                      f"— {inside} the published {br['published']:.3f}")
            print(f"      Rp = {d['planet_radius_rearth']:.2f} R_earth "
                  f"= {d['planet_radius_rjup']:.3f} R_jup (from raw depth)")
            print(f"      a  = {d['semi_major_axis_au']:.4f} AU "
                  f"({d['scaled_semi_major_axis']:.1f} stellar radii)")
            print(f"      Teq ~ {d['equilibrium_temperature_k']:.0f} K (albedo 0.3 assumed)")
            if ld:
                print(f"      limb-darkening: naive sqrt(depth) overestimates "
                      f"Rp/R* by {ld['overestimate_percent']:.1f}%")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--target", action="append", help="run only this target (repeatable)")
    parser.add_argument("--no-figures", action="store_true", help="skip figure generation")
    parser.add_argument("--refresh-archive", action="store_true",
                        help="re-query the NASA Exoplanet Archive instead of using the cache")
    args = parser.parse_args(argv)

    targets = TARGETS
    if args.target:
        targets = [target_by_name(t) for t in args.target]

    if args.refresh_archive:
        published_all = fetch_stellar_parameters(targets)
    else:
        try:
            published_all = load_stellar_parameters()
        except FileNotFoundError:
            print("no cached archive parameters; querying now")
            published_all = fetch_stellar_parameters(targets)

    results = []
    for target in targets:
        lc = fetch_light_curves(target, quiet=True)
        pub = published_all["planets"].get(target.planet, {})
        result, art = analysis.analyse(target, lc, published=pub)
        results.append(result)

        if not args.no_figures:
            written = figures.make_all(result, art, target)
            print(f"  wrote {len(written)} figure files")

    analysis.save_results(results)
    if not args.no_figures and len(results) > 1:
        figures.figure_comparison(results)

    print_summary(results)
    print(f"\nFigures in {FIGURE_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
