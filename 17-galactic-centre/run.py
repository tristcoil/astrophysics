"""Run the whole project end to end.

    python fetch_data.py     # once: ~18 MB (the S2 data plus Do et al.'s chains)
    python run.py            # every fit, measurement and figure

The first run samples ten posteriors and re-fits seven analysis variants,
about half an hour on four cores; afterwards the chains are cached in
``data/chains/`` and a re-run takes about two minutes, offline. ``--resample``
throws the cache away. Every headline number is printed next to the published
value it is graded against.
"""

from __future__ import annotations

import argparse
import warnings

warnings.filterwarnings("ignore")

import analysis  # noqa: E402
from config import FIGURE_DIR, PUBLISHED  # noqa: E402


def _rule(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def pm(s: dict, d: int = 3) -> str:
    return f"{s['median']:.{d}f} ± {0.5 * (s['hi68'] - s['lo68']):.{d}f}"


def print_summary(res: dict) -> None:
    dat = res["data"]
    _rule("THE DATA")
    for key, label in (("keck_astro", "Keck positions"), ("vlt_astro", "NTT/VLT positions")):
        d = dat[key]
        parts = ", ".join(f"{k} {v['n']} (median error {v['median_error_mas']:.2f} mas)"
                          for k, v in d["by_instrument"].items())
        print(f"  {label}: {d['n_used']} used of {d['n_total']}, {d['first']:.1f}-{d['last']:.1f}: {parts}")
    r = dat["keck_rv"]
    print(f"  radial velocities: {r['n']} ({r['first']:.1f}-{r['last']:.1f}), from "
          + ", ".join(f"{k} {v['n']}" for k, v in r["by_instrument"].items()))
    ov = dat["rv_overlap"]
    print(f"  ({ov['n_identical_in_ucla_file']} of the VLT team's {ov['n_vlt_velocities']} SINFONI velocities are in "
          f"the UCLA file unchanged, so each is used once; its {ov['n_early_keck_velocities']} early Keck ones appear "
          f"there re-reduced, by {', '.join(f'{c:+.0f}' for c in ov['early_keck_revisions_kms'])} km/s)")

    f = res["fits"]
    _rule("REPEATING EACH TEAM'S OWN FIT ON ITS OWN DATA")
    for name, which in (("keck_newton", "Newtonian"), ("keck_gr", "GR"), ("keck_redshift", "redshift free")):
        cmp = f[name]["vs_do2019_chain"]
        worst = max(abs(v["difference_in_sigma"]) for v in cmp.values() if v["difference_in_sigma"] is not None)
        print(f"  Keck, {which:<13}: M = {pm(f[name]['parameters']['gm6'])}, R0 = {pm(f[name]['parameters']['r0'])} kpc"
              f"  | Do et al.'s chain: {cmp['gm6']['theirs_mean']:.3f} ± {cmp['gm6']['theirs_sd']:.3f}, "
              f"{cmp['r0']['theirs_mean']:.3f} ± {cmp['r0']['theirs_sd']:.3f}; worst parameter {worst:.2f}σ")
    v = f["vlt_kepler"]["vs_published"]
    print(f"  VLT, Keplerian     : M = {v['gm6']['ours']:.2f} ± {v['gm6']['ours_err']:.2f}, "
          f"R0 = {v['r0']['ours']:.2f} ± {v['r0']['ours_err']:.2f} kpc, vz = {v['vz_kms']['ours']:.1f} km/s"
          f"  | Gillessen+ 2017: {v['gm6']['published']} ± {v['gm6']['published_err']}, "
          f"{v['r0']['published']} ± {v['r0']['published_err']}, {v['vz_kms']['published']}")

    _rule("MASS AND DISTANCE (general relativity, each team's frame free)")
    for name in ("keck_gr", "vlt_gr", "joint_gr"):
        p = f[name]["parameters"]
        print(f"  {f[name]['label']:<34} M = {pm(p['gm6'])} million suns   R0 = {pm(p['r0'])} kpc   "
              f"(M ∝ R0^{f[name]['mass_distance_exponent']:.2f}, correlation {f[name]['mass_distance_correlation']:.3f})")
    for key in ("gravity2022", "g17_multistar"):
        pub = PUBLISHED[key]
        print(f"  {pub['ref']:<60} M = {pub['gm6'][0]}, R0 = {pub['r0'][0]}")
    for name, m in res["headline"]["mass_at_gravity_distance"].items():
        print(f"  {name}: at GRAVITY's R0 = {m['r0']} kpc our banana gives M = {m['gm6']:.3f} ± {m['gm6_err']:.3f} "
              f"(GRAVITY: {m['gravity2022_gm6']} ± {m['gravity2022_gm6_err']})")
    sy = res["systematics"]
    g = sy["summary"]["joint_gr"]
    print(f"  systematics (joint GR fit re-done {len(sy['variants']['joint_gr']) - 1} other ways): R0 moves by up to "
          f"{g['r0_max_shift']:.3f} kpc (rms {g['r0_systematic_rms']:.3f}), M by up to {g['gm6_max_shift']:.3f} "
          f"(rms {g['gm6_systematic_rms']:.3f})")
    for k, row in sy["variants"]["joint_gr"].items():
        print(f"    {k:<26} M = {row['gm6']:.3f} ± {row['gm6_err']:.3f}   R0 = {row['r0']:.3f} ± {row['r0_err']:.3f}")
    nf = res["naive_frames"]
    print(f"  one shared frame for both teams: χ² worse by {nf['delta_chi2']:.0f}; M = "
          f"{nf['one_shared_frame']['gm6']:.3f} vs {nf['separate_frames']['gm6']:.3f}, R0 = "
          f"{nf['one_shared_frame']['r0']:.3f} vs {nf['separate_frames']['r0']:.3f} with separate frames")

    k = res["kepler_by_hand"]
    _rule("KEPLER'S THIRD LAW BY HAND")
    print(f"  a = {k['a_mas']:.2f} mas = {k['a_au']:.0f} AU at {k['r0_kpc']:.2f} kpc; P = {k['period_yr']:.3f} yr")
    print(f"  M = a³/P² = {k['mass_msun']:.3e} suns  (it would be {k['mass_if_r0_were_7_kpc']:.2e} at 7 kpc, "
          f"{k['mass_if_r0_were_9_kpc']:.2e} at 9 kpc: positions alone cannot weigh it)")

    _rule("TESTING EINSTEIN")
    rt = res["relativity"]
    for key, lab in (("keck_redshift", "Keck positions"), ("joint_redshift", "Keck + VLT positions")):
        r = rt[key]
        print(f"  redshift factor f ({lab}): {r['f_median']:.2f} ± {r['f_sd']:.2f}  -> "
              f"{r['sigma_from_newton']:.1f}σ from Newton (0), {r['sigma_from_einstein']:+.1f}σ from Einstein (1)")
    print(f"  published: Do+ 2019 f = 0.88 ± 0.17; GRAVITY 2018 f = 0.90 ± 0.09 ± 0.15")
    fr = res["systematics"]["summary"]["joint_redshift"]
    print(f"  f under the {len(res['systematics']['variants']['joint_redshift']) - 1} analysis variants: "
          f"{fr['f_range'][0]:.2f} to {fr['f_range'][1]:.2f} (rms shift {fr['f_systematic_rms']:.2f})")
    for k, row in res["systematics"]["variants"]["joint_redshift"].items():
        print(f"    {k:<26} f = {row['f']:.2f} ± {row['f_err']:.2f}")
    r = rt["joint_precession"]
    print(f"  precession factor f_SP: {r['f_sp_median']:.2f} ± {r['f_sp_sd']:.2f}  (GRAVITY 2020 with VLTI data: "
          "1.10 ± 0.19)")
    pp = res["precession_profile"]
    print(f"    profile: with all positions Einstein fits worse than the best f_SP by Δχ² = "
          f"{pp['all_positions']['delta_chi2_gr_vs_best']:.1f}; with adaptive-optics positions only by "
          f"{pp['ao_only']['delta_chi2_gr_vs_best']:.1f} (Newton vs Einstein: {pp['ao_only']['delta_chi2_newton_vs_gr']:.1f})"
          " -> a false signal from the speckle era")
    for key in ("keck_gr_vs_keck_newton", "joint_gr_vs_joint_newton"):
        d = rt[key]
        print(f"  {key}: Δχ² = {d['delta_chi2']:.1f}, ΔBIC = {d['delta_bic']:.1f} in favour of GR; "
              f"Newton's orbit shifts M by {d['mass_shift']:+.3f}, R0 by {d['r0_shift']:+.3f} kpc")
    ro = res["roemer_delay_days"]
    print(f"  Roemer delay: {ro['min_days']:+.1f} days (apocentre, in front) to {ro['max_days']:+.1f} "
          "days (pericentre, behind)")

    d = f["joint_gr"]["derived"]
    _rule("S2 AT PERICENTRE, AND WHY IT MUST BE A BLACK HOLE")
    print(f"  pericentre {pm(d['pericentre_au'], 1)} AU = {d['pericentre_light_hours']['median']:.1f} light-hours "
          f"= {d['pericentre_in_rs']['median']:.0f} Schwarzschild radii, at {pm(d['pericentre_speed_kms'], 0)} km/s "
          f"({100 * d['pericentre_speed_frac_c']['median']:.2f}% of c)")
    print(f"  relativistic redshift there: {d['redshift_at_pericentre_kms']['median']:.0f} km/s; precession "
          f"{d['precession_arcmin_per_orbit']['median']:.1f}′ per orbit")
    c = res["black_hole_case"]
    print(f"  mean density inside pericentre ≥ {c['rho_mean_inside_pericentre']:.2e} M☉/pc³ "
          f"({c['density_gain_since_1997']:.0f}× the 1997 limit)")
    print(f"  longest-lived dark cluster: 1997 limit {c['max_lifetime_1997']['lifetime_yr']:.1e} yr; with S2 "
          f"{c['max_lifetime_now']['lifetime_yr']:.1e} yr ({c['max_lifetime_now']['class']}, "
          f"{c['max_lifetime_now']['m_star']:.2g} M☉); 10 Gyr needs ρ < {c['density_needed_for_10_gyr']:.1e}")
    print(f"  exotic loophole left: primordial black holes lighter than {c['primordial_bh_mass_limit_msun']:.1e} M☉ "
          f"(was {c['primordial_bh_mass_limit_1997_msun']:.1e} in 1997), or a boson star")

    _rule("PREDICTION FOR THE EVENT HORIZON TELESCOPE (project I5)")
    for name in ("keck_gr", "vlt_gr", "joint_gr"):
        s = f[name]["derived"]["shadow_diameter_uas"]
        g = f[name]["derived"]["theta_g_uas"]
        print(f"  {f[name]['label']:<34} shadow {pm(s, 1)} µas (θg = {pm(g, 2)} µas)")
    e = PUBLISHED["eht2022"]
    print(f"  EHT 2022: ring {e['ring_uas'][0]} ± {e['ring_uas'][1]} µas; shadow {e['shadow_uas'][0]} ± "
          f"{e['shadow_uas'][1]} µas")

    _rule("MCMC HEALTH")
    for name, s in f.items():
        m = s["mcmc"]
        print(f"  {name:<18} {m['walkers']} walkers × {m['steps']} steps, acceptance {m['acceptance']:.2f}, "
              f"longest autocorrelation {m['max_autocorr_steps']:.0f} steps ({m['steps_per_autocorr']:.0f}× shorter "
              "than the run)")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--resample", action="store_true", help="ignore cached chains and sample again")
    parser.add_argument("--no-figures", action="store_true", help="skip the figures")
    parser.add_argument("--no-viewer", action="store_true", help="skip the interactive orbit explorer")
    args = parser.parse_args(argv)

    results, art = analysis.run_all(resample=args.resample)
    analysis.save_results(results)
    if not args.no_figures:
        import figures

        written = figures.make_all(results, art)
        print(f"wrote {len(written)} figure files")
    if not args.no_viewer:
        try:
            import viewer
        except ImportError:
            viewer = None
        if viewer is not None:
            path, n = viewer.build_viewer(results, art)
            print(f"wrote {path} ({n / 1e6:.2f} MB of embedded data)")
    print_summary(results)
    print(f"\nFigures in {FIGURE_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
