"""Run the whole project end to end.

    python run.py

Fetches whatever is not cached yet, measures, writes every figure and
``data/results.json``, and prints a summary. Once the cache is warm it runs
offline in a few minutes (the 4-D posterior grid is computed once and cached).

    python run.py --no-figures
"""

from __future__ import annotations

import argparse
import sys
import time
import warnings

warnings.filterwarnings("ignore")

import numpy as np

import analysis
import config
import environment
import fetch_data
import hosts


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--no-figures", action="store_true", help="skip figure generation")
    args = parser.parse_args(argv)
    t0 = time.time()

    for group in fetch_data.GROUPS:
        if fetch_data.missing(group):
            print(f"Fetching {group}")
            fetch_data.fetch_all(only=group)

    print("Census")
    census = analysis.census()
    census["redshift_check"] = analysis.redshift_check()
    print("Length distribution")
    lres, lart = analysis.length_distribution()
    sim = analysis.simulation_check()
    print("Our own lengths from the images")
    measured, mart = analysis.measure_all()
    print("DR2 against DR3")
    releases = analysis.compare_releases(measured)
    print("Lobe spectra and ages")
    ares, aart = analysis.ageing_all(mart)
    print("Host galaxies and black holes")
    hres = hosts.compare_black_holes()
    hres_oei = hosts.compare_with_oei2023()
    mres = hosts.compare_stellar_masses()
    print("Environments")
    bres = environment.borg_densities()
    lumchk = environment.borg_luminosity_check()
    eres = environment.compare_environments()

    written = []
    if not args.no_figures:
        import figures

        print("Figures")
        written += figures.figure_gallery(measured, mart)
        written += figures.figure_lineup(measured, mart)
        written += figures.figure_porphyrion(measured, mart, lres)
        written += figures.figure_records(measured, mart)
        written += figures.figure_census(census)
        written += figures.figure_lengths(lres, lart)
        written += figures.figure_ageing(ares, aart)
        written += figures.figure_hosts(hres)
        written += figures.figure_environment(bres, lumchk, eres)
        written += figures.figure_validation(census, measured, sim, ares, hres_oei, bres)
        written += figures.figure_gallery(measured, mart, group="dr3")
        written += figures.figure_releases(releases)
        import interactive

        written.append(interactive.write_explorer(measured, mart))
        print(f"  wrote {len(written)} files to {config.FIGURE_DIR.name}/")

    analysis.save_results(census=census, length_distribution=lres, simulation_check=sim, giants=measured,
                          dr2_vs_dr3=releases,
                          ageing=ares, black_holes=hres, black_holes_vs_oei2023=hres_oei, stellar_masses=mres,
                          borg=bres, borg_luminosity_check=lumchk, environment_sdss=eres)
    print_summary(census, lres, sim, measured, ares, hres, hres_oei, mres, bres, lumchk, eres, releases)
    print(f"\nDone in {time.time() - t0:.0f} s. Results in {config.RESULTS_PATH.relative_to(config.PROJECT_DIR)}")
    return 0


def print_summary(census, lres, sim, measured, ares, hres, hres_oei, mres, bres, lumchk, eres, releases) -> None:
    print("\n" + "=" * 96)
    print("I18 -- Giant radio galaxies")
    print("=" * 96)
    print(f"\nCensus: {census['n_giants']:,} known giants; median length {census['median_lp_mpc']:.2f} Mpc, median z "
          f"{census['median_z']:.2f}; {census['n_above']['3']} above 3 Mpc, {census['n_above']['5']} above 5 Mpc; "
          f"{census['n_larger_than_moon']} larger on the sky than the full Moon")
    cc = census["cosmology"]
    print(f"  catalogue lengths reproduced from angle and z: median ratio {cc['catalogue_vs_planck_median_ratio']:.5f}")
    for c in config.COSMOLOGIES:
        print(f"  {c.name:<24} lengths x{cc[c.name]['median_ratio_to_planck']:.3f}; giants still above 0.7 Mpc: "
              f"{cc[c.name]['n_giants']:,}")
    print(f"  a 0.7 Mpc giant never appears smaller than {census['min_angle_0p7mpc_arcmin']:.2f}' "
          f"(at z = {census['min_angle_redshift']:.2f})")
    rc = census.get("redshift_check")
    if rc:
        print(f"  catalogue redshift vs the SDSS spectrum at the host ({rc['n_with_spectrum']:,} giants): "
              f"{rc['frac_within_0p01'] * 100:.1f}% agree to 1%; {rc['n_off_by_10pct']} ({rc['frac_off_by_10pct'] * 100:.1f}%) "
              f"are off by >10% ({rc['off_spectral_classes']}); {rc['n_no_longer_giants']} of these fall below 0.7 Mpc "
              f"with the SDSS redshift")

    print("\nLength distribution (Oei et al. 2023 sample and model)")
    for zm, r in lres.items():
        if not zm.startswith("0."):
            continue
        p = r["posterior"]["xi"]
        pub = r["published"]
        print(f"  z < {zm}: {r['n']:,} giants; naive xi {r['naive_xi']:.2f} +- {r['naive_xi_err']:.2f}; with selection "
              f"{p['mean']:.2f} +- {p['sd']:.2f} (mode {p['mode']:.2f}; published {pub['xi'][0]} +- {pub['xi'][1]})")
        for k in ("half_width", "b_ref", "sigma_ref"):
            q = r["posterior"][k]
            print(f"      {k:<11} {q['mean']:8.2f} +- {q['sd']:.2f}   published {pub[k if k != 'half_width' else 'half_width_arcmin']}")
    print(f"  simulations: mode bias {sim['mode_bias']:+.2f}, mean bias {sim['mean_bias']:+.2f}, 68% coverage "
          f"{sim['coverage_68'] * 100:.0f}%")
    d = lres["porphyrion_deprojection"]
    print(f"  Porphyrion in 3D (xi = {d['xi_used']:.2f}): median {d['median']:.2f}, mean {d['mean']:.2f} Mpc "
          f"(68%: {d['p16']:.2f}-{d['p84']:.2f}); published 6.8 / 7.28")
    t = lres["tail"]
    print(f"  of all giants, {t['fraction_above_5mpc_of_giants'] * 100:.2f}% exceed 5 Mpc, "
          f"{t['fraction_above_porphyrion'] * 100:.2f}% exceed Porphyrion")

    print("\nOur lengths (3-sigma extent on LoTSS DR2)")
    print(f"  {'giant':<24}{'z':>7}{'6\"':>8}{'20\"':>8}{'60\"':>8}{'adopted':>9}{'Mpc':>7}   claims (arcmin / Mpc)")
    for r in measured:
        at = r["at"]
        f = lambda k: f"{at[k]['las_arcmin']:.1f}" if k in at and at[k]["complete"] else "--"
        cl = "; ".join(f"{c['las_arcmin'] or '-'}'/{c['lp_mpc']}" for c in r["claims"])
        tag = " [DR3]" if r.get("survey") == "dr3" else ""
        print(f"  {r['name'] + tag:<24}{r['z']:>7.3f}{f('6'):>8}{f('20'):>8}{f('60'):>8}{r['adopted']['las_arcmin']:>8.1f}'"
              f"{r['adopted']['lp_mpc']:>7.2f}   {cl}")
    print(f"\nDR2 against DR3 ({releases['n_complete']} giants, same seeds): lengths x{releases['las_ratio_median']:.3f} "
          f"(range {releases['las_ratio_range'][0]:.3f}-{releases['las_ratio_range'][1]:.3f}), flux x"
          f"{releases['flux_ratio_median']:.3f} (range {releases['flux_ratio_range'][0]:.2f}-"
          f"{releases['flux_ratio_range'][1]:.2f}), noise x{releases['noise_ratio_median']:.2f}")

    print("\nLobe ages (LoTSS 144 MHz vs NVSS 1.4 GHz at 45\")")
    for k, r in ares.items():
        names = r.get("lobe_names", {"A": "A", "B": "B"})
        for lobe, v in r["lobes"].items():
            a = v["alpha"]
            s = (f"alpha {a:.2f} +- {v['alpha_err']:.2f}" if a is not None and np.isfinite(a)
                 else f"alpha < {v['alpha_upper']:.2f}")
            age = v.get("ci_age_myr")
            ages = (f"; CI age {age:.0f} Myr (max {v['ci_age_max_myr']:.0f})" if age is not None and np.isfinite(age)
                    else (f"; {v['ci_note']}" if v.get("ci_note") else ""))
            print(f"  {k:<10} {names[lobe]:<3} S144 {v['flux_144_jy']:.3f} Jy, S1400 {v['flux_1400_jy']:.3f} Jy, {s}, "
                  f"B_eq {v['b_eq_ug']:.2f} uG{ages}")
    if "3c236" in ares:
        mc = ares["3c236"]["model_check_with_their_fluxes"]
        print(f"  3C 236 model check (Shulevski et al. fluxes): NW {mc['NW']['ci_age_myr']:.0f} Myr (129), "
              f"SE {mc['SE']['ci_age_myr']:.0f} Myr (117)")

    print("\nBlack holes (M-sigma, SDSS DR17)")
    print(f"  {hres['giants_all']['n_with_mbh']:,} giant hosts: median {hres['giants_all']['median_mbh_1e9']:.2f} x 1e9 Msun "
          f"(10-90%: {hres['giants_all']['p10_p90_1e9'][0]:.2f}-{hres['giants_all']['p10_p90_1e9'][1]:.2f})")
    print(f"  vs Oei et al. 2023 for {hres_oei['n_both']} hosts: {hres_oei['median_log_ratio']:+.3f} dex "
          f"(scatter {hres_oei['scatter_dex']:.3f})")
    m = hres["matched"]
    print(f"  giants - ordinary radio galaxies, matched in z and L144: {m['median_difference_dex']:+.3f} +- "
          f"{m['err_dex']:.3f} dex ({m['n_giants']} giants)")
    print(f"  stellar masses, same matching: {mres['median_difference_dex']:+.3f} +- {mres['err_dex']:.3f} dex")
    mv = hres["mbh_vs_length"]
    print(f"  black hole mass vs length: Spearman rho {mv['spearman_rho']:.3f} (p = {mv['p']:.2f})")

    print("\nEnvironments")
    for method in ("fixed", "flex"):
        b = bres[method]
        print(f"  BORG {method:<5}: giants Gamma(k={b['giants']['gamma_k']:.2f}, theta={b['giants']['gamma_theta']:.2f}) "
              f"mean {b['giants']['gamma_mean']:.2f}; radio galaxies mean {b['radiogalaxies']['gamma_mean']:.2f}; "
              f"KS p {b['ks_p']:.1e}")
        lc = lumchk[method]
        print(f"         vs radio galaxies > 1e24 W/Hz (n = {lc['n_bright']}): {lc['giants_mean']:.2f} vs "
              f"{lc['bright_rg_mean']:.2f}, KS p {lc['ks_p_giants_vs_bright']:.2f}")
    for k in ("vs_radiogalaxies_z", "vs_radiogalaxies_z_L", "vs_lookalikes"):
        v = eres.get(k)
        if v:
            print(f"  SDSS neighbours {k:<22}: ratio {v['ratio']:.2f} +- {v['err']:.2f} ({v['n_giants']} giants)")


if __name__ == "__main__":
    sys.exit(main())
