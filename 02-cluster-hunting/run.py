"""Run the whole project end to end.

    python run.py

Reads the cached Gaia data (run ``fetch_data.py`` once first), hunts each
cluster, writes every figure, and prints what was found beside the published
values. Needs no network once the cache is warm.
"""

from __future__ import annotations

import argparse
import warnings

warnings.filterwarnings("ignore")

import analysis
from config import CLUSTERS, FIGURE_DIR, cluster_by_key
from fetch_data import (
    load_gaia, load_literature, load_named_stars, load_published_members,
)


def _row(label, c, spec):
    if not c:
        return
    ours = format(c["ours"], spec)
    err = f" ± {format(c['ours_err'], spec)}" if c.get("ours_err") is not None else ""
    pub = format(c["published"], spec)
    perr = f" ± {format(c['published_err'], spec)}" if c.get("published_err") is not None else ""
    print(f"    {label:<28}{ours + err:>22}{pub + perr:>24}   {c['source']}")


def print_summary(results) -> None:
    print("\n" + "=" * 100)
    print("RESULTS -- what we found, beside the published values")
    print("=" * 100)
    for r in results:
        k, m, d, s, h = r.knot, r.motion, r.distance, r.structure, r.hr_diagram
        mem, kin = r.membership, r.kinematic_only
        print(f"\n{r.cluster}  ({r.difficulty} case, {r.n_stars:,} stars in the field)")
        print(f"  knot: ({k['pmra']:+.2f}, {k['pmdec']:+.2f}) mas/yr, {k['contrast']:.0f}x denser "
              f"than its surroundings; members' parallaxes agree to "
              f"{k['parallax_spread_fraction']:.0%}")
        ru = k.get("runner_up")
        if ru:
            print(f"        next-best clump: contrast {ru['contrast']:.1f}x, parallax spread "
                  f"{ru['parallax_spread_fraction']:.0%} -- the field's own pile-up, not a cluster")
        print(f"  motion alone: {kin['n']} candidates; adding parallax confirms "
              f"{kin['n_confirmed_by_parallax']} and rejects {kin['n_rejected_by_parallax']} "
              f"(field stars sharing the motion at another distance)")
        halo = (f" ({mem['n_core']:.0f} core + {mem['n_halo']:.0f} halo)"
                if mem.get("n_halo") else "")
        print(f"  motion + parallax: {mem['n']} members{halo}; expected contamination "
              f"{mem['expected_contamination']:.0f}, expected missed {mem['expected_missed']:.0f}")
        print(f"  distance: {d['distance_pc']:.1f} ± {d['distance_err_pc']:.1f} pc "
              f"(parallax {d['parallax_mas']:.4f} mas after a {d['zero_point_mean_mas'] * 1000:+.0f} "
              f"uas zero-point correction)")
        print(f"  space velocity: {m['speed_kms']:.2f} km/s relative to the Sun, heading for "
              f"RA {m['convergent_point_deg'][0]:.1f}, Dec {m['convergent_point_deg'][1]:+.1f}; "
              f"UVW = ({m['uvw_kms'][0]:+.1f}, {m['uvw_kms'][1]:+.1f}, {m['uvw_kms'][2]:+.1f}) km/s")
        ds = m["depth_slope"]
        if s["depth_resolvable"]:
            print(f"  nearer members move faster: d(pm)/d(parallax) = "
                  f"({ds['pmra_c']['slope']:+.2f} ± {ds['pmra_c']['err']:.2f}, "
                  f"{ds['pmdec_c']['slope']:+.2f} ± {ds['pmdec_c']['err']:.2f}) vs "
                  f"({m['predicted_slope'][0]:+.2f}, {m['predicted_slope'][1]:+.2f}) predicted")
            depth = (f"depth sigma {s['depth_sigma_pc']:.1f} pc vs width "
                     f"{s['width_sigma_pc']:.1f} pc across the sky")
        else:
            print(f"  depth-motion test: not possible here -- {s['depth_note']}")
            depth = "depth unresolvable at this distance"
        print(f"  size: half the members within {s['r50_deg']:.2f} deg = {s['r50_pc']:.1f} pc; "
              f"{depth}")
        print(f"  internal velocity spread: {m['velocity_dispersion_kms']:.2f} km/s "
              f"(an upper limit: binaries' orbital motion adds to it)")
        print(f"  main sequence: {h['width_members_mag']:.3f} mag wide vs "
              f"{h['width_field_mag']:.3f} for the field -- {h['width_ratio']:.1f}x thinner")

        print(f"    {'quantity':<28}{'ours':>22}{'published':>24}   source")
        cmp_ = r.comparison
        _row("pmRA (mas/yr)", cmp_.get("pmra"), ".3f")
        _row("pmDec (mas/yr)", cmp_.get("pmdec"), ".3f")
        _row("parallax, uncorrected (mas)", cmp_.get("parallax_uncorrected"), ".4f")
        _row("distance (pc)", cmp_.get("distance_pc"), ".1f")
        _row("radial velocity (km/s)", cmp_.get("rv_kms"), ".2f")
        _row("radial velocity (km/s)", cmp_.get("rv_kms_gc2018"), ".2f")
        _row("r50 (deg)", cmp_.get("r50_deg"), ".3f")
        _row("speed (km/s)", cmp_.get("speed_kms"), ".2f")
        if cmp_.get("convergent_point_offset_deg"):
            print(f"    convergent point differs from GC2018 by "
                  f"{cmp_['convergent_point_offset_deg']['ours']:.2f} deg")
        mv = r.membership_vs_published
        if mv:
            print(f"  vs Hunt & Reffert (2023) members in the same footprint: we recover "
                  f"{mv['recovered_fraction']:.0%} of their {mv['published_in_footprint']}; "
                  f"{mv['ours_confirmed_fraction']:.0%} of ours are in their list")
        for name, f in r.famous_stars.items():
            if f.get("in_sample"):
                verdict = "member" if f["member"] else "NOT selected"
                print(f"    {name:<9} G={f['g']:.2f}  p={f['p_member']:.2f}  RUWE={f['ruwe']:.1f}"
                      f"  -> {verdict}")
            else:
                print(f"    {name:<9} G={f['g']:.2f}  -> {f['reason']}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--cluster", action="append", help="run only this cluster (repeatable)")
    parser.add_argument("--no-figures", action="store_true", help="skip figure generation")
    args = parser.parse_args(argv)

    clusters = [cluster_by_key(k) for k in args.cluster] if args.cluster else CLUSTERS
    literature = load_literature()

    results, artefacts = [], []
    for cluster in clusters:
        result, art = analysis.analyse(
            cluster, load_gaia(cluster), literature,
            load_published_members(cluster), load_named_stars(cluster),
        )
        results.append(result)
        artefacts.append(art)

    analysis.save_results(results)

    if not args.no_figures:
        import figures

        written = figures.make_all(results, artefacts, clusters, literature)
        print(f"\nwrote {len(written)} figure files to {FIGURE_DIR}")

    print_summary(results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
