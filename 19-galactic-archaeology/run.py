"""Run the whole project end to end.

    python run.py                  # all three parts, every figure
    python run.py --part gd1       # one part: gd1, gse, spiral (repeatable)
    python run.py --no-figures

Reads the cached Gaia data (run ``fetch_data.py`` once first), writes
data/results.json and figures/, and prints the headline numbers beside the
published ones. Needs no network once the cache is warm (the McMillan 2017
potential table is built from galpy on the first run and cached too).
"""

from __future__ import annotations

import argparse
import gc
import time
import warnings

warnings.filterwarnings("ignore")

import analysis  # noqa: E402
from config import FIGURE_DIR, REFERENCE  # noqa: E402


def summary(results) -> None:
    print("\n" + "=" * 96)
    print("RESULTS -- what we found, beside the published values")
    print("=" * 96)
    if "gd1" in results:
        r = results["gd1"]
        k = REFERENCE["koposov2010"]
        print("\nGD-1")
        print(f"  field: {r['field']['n_stars']:,} stars; members {r['selection']['extent']:,} "
              f"({r['n_band']:,} in the stream band, ~{r['n_band_contamination']:.0f} of them field stars)")
        print(f"  detected from phi1 = {r['extent_phi1'][0]:.0f} to {r['extent_phi1'][1]:.0f} deg")
        lin = r["distance_model"]["linear"]
        print(f"  distance gradient: {lin[0]:.4f} +- {r['distance_model']['linear_err'][0]:.4f} "
              "mag/deg (Koposov 2010: ~0.012)")
        for d in r["distance"]:
            print(f"    phi1 {d['phi1']:6.1f}: DM {d['dm']:.3f} +- {d['dm_err']:.3f}")
        print(f"    Koposov 2010: {list(zip(k['distance_phi1'], k['distance_kpc']))}")
        ps_ = r["parallax_summary"]
        print(f"  parallax check: offset {ps_['mean_offset_mas']:+.4f} +- {ps_['mean_offset_err_mas']:.4f} mas, "
              f"chi2 {ps_['chi2_stat']:.1f}/{ps_['dof']}")
        for g in r["gaps"]:
            print(f"  gap at {g['phi1']:6.1f}: {100 * g['depth_fraction']:.0f}% deep, {g['sigma']:.1f} sigma")
        print(f"  spur {r['spur']['n_feature']} vs {r['spur']['n_mirror']} (p = {r['spur']['p_value']:.1e}); "
              f"blob {r['blob']['n_feature']} vs {r['blob']['n_mirror']}; "
              f"control {r['spur_control']['n_feature']} vs {r['spur_control']['n_mirror']}")
        print(f"  Huang et al. 2019 members: {r['huang2019']['fate']}")
        print(f"  Gaia RVs of band members: {r['rv_check']['n_consistent']} of {r['rv_check']['n']} "
              "on the published velocity curve")
    if "gse" in results:
        r = results["gse"]
        print("\nGaia-Sausage-Enceladus")
        print(f"  halo stars (|v - v_LSR| > 210) within 2.5 kpc: {r['n_halo_2p5kpc']:,} "
              f"({100 * r['halo_fraction_random']:.2f}% of all 6D stars)")
        kb = r["halo_kinematic_ellipsoid"]
        print(f"  velocity-cut halo: beta = {kb['beta']:.2f}, sigma = {[round(x) for x in kb['sigma']]} km/s")
        for b in r["beta_by_metallicity"]:
            if "beta" in b:
                print(f"    [M/H] {b['feh_lo']:5.2f}..{b['feh_hi']:5.2f}: n={b['n']:6d} beta = "
                      f"{b['beta']:.2f} +- {b['beta_err']:.2f} (halo share {b['weight']:.2f})")
        print(f"  substructure (5 kpc): {r['substructure_counts_5kpc']}")
        print(f"  GSE chemistry: {r['gse_chemistry']}")
        c = r["alpha_difference_combined"]
        print(f"  [alpha/Fe] GSE - thick disc at fixed [M/H]: {c['value']:+.3f} +- {c['err']:.3f}")
    if "spiral" in results:
        r = results["spiral"]
        a = r["age_summary"]
        print("\nPhase spiral")
        print(f"  ring: {r['n_ring_good']:,} stars (Antoja 2018: 935,590 in DR2)")
        print(f"  age of the disturbance: {a['min_myr']:.0f}-{a['max_myr']:.0f} Myr, median "
              f"{a['median_myr']:.0f} (Antoja 2018: 300-900, ~500)")
        print(f"  Sigma_1.1: {r['surface_density_11']} vs measured {r['surface_density_11_measured']}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--part", action="append", choices=["gd1", "gse", "spiral"])
    parser.add_argument("--no-figures", action="store_true")
    parser.add_argument("--n-boot", type=int, default=200)
    args = parser.parse_args(argv)
    t0 = time.time()
    parts = args.part or ("gd1", "gse", "spiral")
    # One part at a time, freeing its (large) star tables once its figures
    # are drawn: all three at once need ~4 GB. The overview figure keeps only
    # a few small arrays from each.
    results, keep = {}, {}
    n_written = 0
    for part in parts:
        res, art = analysis.run((part,), n_boot=args.n_boot)
        results.update(res)
        analysis.save_results(res)
        if not args.no_figures:
            import figures

            n_written += len(figures.make_all(res, art))
            keep[part] = figures.overview_arrays(part, res[part], art[part])
        del art
        gc.collect()
    if not args.no_figures and all(p in keep for p in ("gd1", "gse", "spiral")):
        import figures

        n_written += len(figures.overview(results, keep))
    if not args.no_figures and "gse" in parts:
        import interactive

        interactive.write_explorer()
        n_written += 1
    if not args.no_figures:
        print(f"\nwrote {n_written} figure files to {FIGURE_DIR}")
    summary(results)
    print(f"\n(total {time.time() - t0:.0f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
