"""Run the whole project end to end.

    python run.py

Downloads whatever is not already cached, runs the analysis, writes every
figure and the interactive 3D map, and prints a summary of the measurements
next to the published values they are graded against. Once the cache is warm
this runs offline. The first run spends a few minutes counting galaxy pairs;
later runs reuse the cached counts.
"""

from __future__ import annotations

import argparse
import sys
import time
import warnings

warnings.filterwarnings("ignore")

import analysis
import literature as lit
from config import CATALOG_PATH, FIGURE_DIR
from fetch_data import build_catalog, fetch_abell, fetch_survey, load_catalog


def print_summary(res: dict) -> None:
    line = "=" * 78
    print(f"\n{line}\nRESULTS\n{line}")
    vl = res["volume_limited"]
    print(f"\nCatalogue: {res['n_galaxies']:,} main-sample galaxies "
          f"(median z = {res['median_z']:.4f}; Strauss et al. quote {lit.STRAUSS_MGS['median_z']})")
    print(f"  survey footprint {res['mask']['surveyed_deg2']:,.0f} deg2 "
          f"({res['mask']['interior_deg2']:,.0f} deg2 away from the edges)")
    print(f"  equatorial slice: {res['slice']['n']:,} galaxies, "
          f"{res['slice']['n_volume_limited']:,} of them in the volume-limited sample")

    print(f"\nVolume-limited sample (Mr - 5 log h < {vl['mag_limit_lit']}, "
          f"{vl['z_range'][0]} < z < {vl['z_range'][1]:.4f})")
    print(f"  {vl['n']:,} galaxies ({vl['n_interior']:,} in the clustering footprint); "
          f"Zehavi et al.: {lit.ZEHAVI_MR20['n_gal']:,}")
    print(f"  number density {vl['number_density_h3'] * 1e3:.2f} vs "
          f"{vl['number_density_zehavi_h3'] * 1e3:.2f} x 1e-3 h^3 Mpc^-3 "
          f"(ratio {vl['number_density_h3'] / vl['number_density_zehavi_h3']:.2f})")
    c = vl["completeness"]
    print(f"  [{'PASS' if c['complete'] else 'FLAG'}] complete to z_max: a threshold galaxy "
          f"there appears at r <= {max(c['r_at_zmax']):.3f} (limit {c['r_limit']}; "
          f"margin {c['margin_mag']:.3f} mag)")
    print(f"  red fraction {vl['red_fraction']:.3f} vs {vl['red_fraction_zehavi']:.3f} (Zehavi, "
          f"their own colour cut)")

    d = res["color_divider"]
    sign = "-" if d["intercept"] < 0 else "+"
    print(f"\nColour divider (from {d['n_bins_used']} bimodal luminosity bins): "
          f"(g-r)_0 = {d['slope']:.4f} M_r {sign} {abs(d['intercept']):.3f}")

    e = res["environment"]
    lo, hi = e["all"]["red_fraction"][0], e["all"]["red_fraction"][-1]
    print(f"\nEnvironment ({e['n']:,} galaxies): red fraction rises from {lo:.2f} in the "
          f"sparsest eighth to {hi:.2f} in the densest")
    for key in ("bright", "faint"):
        f = e[key]["red_fraction"]
        print(f"  {key:>6} (M_lit {'<' if key == 'bright' else '>='} "
              f"{e['luminosity_split_M_lit']}): {f[0]:.2f} -> {f[-1]:.2f}")

    grp = res["groups"]
    print(f"\nGroups: {grp['n_groups_10plus']} with >= 10 bright members "
          f"(linking {grp['linking_lengths_h'][0]:.2f} x {grp['linking_lengths_h'][1]:.2f} h^-1 Mpc)")
    for gg in grp["richest"][:5]:
        name = gg["abell"] or "(no Abell match)"
        print(f"  {name:<18} N={gg['n']:3d}  z={gg['z']:.4f}  sigma_v={gg['sigma_v']:5.0f} km/s  "
              f"looks {gg['elongation']:.1f}x deeper than wide")
    coma, ref = grp.get("coma"), grp.get("coma_colless_dunn")
    if coma and ref:
        print(f"  [CHECK] Coma (Abell 1656): sigma_v = {coma['sigma_v']:.0f} +/- "
              f"{coma['sigma_v_err']:.0f} km/s from {coma['n_dispersion']} galaxies; "
              f"Colless & Dunn (1996): {ref['sigma_v']:.0f} +/- {ref['sigma_v_err']:.0f}")
    print("  richest in the slice:")
    for gg in grp["richest_in_slice"][:4]:
        name = gg["abell"] or "(no Abell match)"
        print(f"  {name:<18} N={gg['n']:3d}  RA={gg['ra']:6.2f}  z={gg['z']:.4f}  "
              f"sigma_v={gg['sigma_v']:5.0f} km/s  elongation {gg['elongation']:.1f}")

    cl = res["clustering"]
    z = cl["zehavi"]
    print(f"\nClustering: {cl['n_data']:,} galaxies, {cl['n_random']:,} random points, "
          f"{cl['area_deg2']:,.0f} deg2")
    for key, f in cl["fits"].items():
        ref = (f"Zehavi {z['r0_diag']:.2f} +/- {z['r0_diag_err']:.2f}, gamma {z['gamma_diag']:.2f}"
               if f["diagonal"] else
               f"Zehavi {z['r0']:.2f} +/- {z['r0_err']:.2f}, gamma {z['gamma']:.2f}")
        print(f"  {key:<22} r0 = {f['r0']:.2f} +/- {f['r0_err']:.2f} h^-1 Mpc, "
              f"gamma = {f['gamma']:.2f} +/- {f['gamma_err']:.2f}, "
              f"chi2/dof = {f['chi2']:.1f}/{f['dof']}  [{ref}]")
    print("  w_p ratio to Zehavi et al. at r_p =")
    print("   " + "  ".join(f"{rp:>5}" for rp in lit.ZEHAVI_WP_RP))
    print("   " + "  ".join(f"{q:5.2f}" for q in cl["ratio_to_zehavi"]))

    w = res["great_wall"]
    for s in w["structures"]:
        print(f"\nLargest structure in the slice at delta > {s['threshold_delta']:.0f}: "
              f"{s['length_mpc']:.0f} Mpc end to end ({s['length_h']:.0f} h^-1 Mpc), "
              f"RA {s['ra_range'][0]:.0f}-{s['ra_range'][1]:.0f}, median z {s['median_z']:.3f}")
    g = lit.GOTT_SGW
    print(f"  Gott et al.: RA {g['ra_range_hours'][0] * 15:.0f}-{g['ra_range_hours'][1] * 15:.0f}, "
          f"z = {g['redshift']}, {g['comoving_length_mpc']:.0f} Mpc along its curve (H0=71)")
    v = res["largest_void"]
    print(f"Largest empty circle in the slice: radius {v['radius_mpc']:.0f} Mpc "
          f"({v['radius_h']:.1f} h^-1 Mpc) at RA {v['ra_center']:.0f}, "
          f"{v['distance_h'] / 0.6766:.0f} Mpc away")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--no-figures", action="store_true", help="skip figure generation")
    parser.add_argument("--workers", type=int, default=None,
                        help="processes for pair counting (default: half the cores, max 8)")
    args = parser.parse_args(argv)

    if not CATALOG_PATH.exists():
        fetch_survey()
        build_catalog()
    arrays, meta = load_catalog()
    abell = fetch_abell(quiet=True)

    t0 = time.time()
    print("Analysis")
    results, art = analysis.analyse(arrays, meta, abell, workers=args.workers)
    analysis.save_results(results)
    print(f"  done in {time.time() - t0:.0f}s")

    if not args.no_figures:
        import figures

        print("\nFigures")
        written = figures.make_all(results, art)
        print(f"  wrote {len(written)} files to {FIGURE_DIR}")

    print_summary(results)
    return 0


if __name__ == "__main__":
    sys.exit(main())
