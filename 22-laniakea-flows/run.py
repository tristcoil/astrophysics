"""Run the whole project end to end.

    python run.py

Fetches whatever is not cached yet (about 33 MB), measures everything, writes
every figure and ``data/results.json``, and prints a summary. The first run
takes about 15 minutes (Lambda-CDM covariances of 2.7 billion galaxy pairs, the
2MRS reconstruction, the beta* profile likelihood); later runs reuse
``data/derived/`` and take a couple of minutes.

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
import fetch_data


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--no-figures", action="store_true", help="skip figure generation")
    args = parser.parse_args(argv)
    t0 = time.time()
    if fetch_data.missing():
        print("Fetching data")
        failures = fetch_data.fetch_all()
        if failures:
            for f in failures:
                print("ERROR:", f, file=sys.stderr)
            return 1
    res, art = analysis.run_all(verbose=True)
    written = []
    if not args.no_figures:
        import figures
        import interactive

        print("Figures")
        for fn in figures.ALL:
            written += fn(res, art)
        written.append(interactive.write_explorer(res, art))
        print(f"  wrote {len(written)} files to {config.FIGURE_DIR.name}/")
    analysis.save_results(res)
    print_summary(res)
    print(f"\nDone in {time.time() - t0:.0f} s. Results in {config.RESULTS_PATH.relative_to(config.PROJECT_DIR)}")
    return 0


def print_summary(res) -> None:
    line = "=" * 96
    print("\n" + line + "\nJ4 -- Laniakea: the flows of galaxies\n" + line)
    fr = res["frames"]
    print(f"\nLocal Group through the CMB (Planck Sun-CMB minus Sun-LG): {fr['lg_cmb_ours']['speed']:.1f} km/s toward "
          f"(l, b) = ({fr['lg_cmb_ours']['l']:.1f}, {fr['lg_cmb_ours']['b']:.1f})")
    print("\nHubble constant (CF4 groups beyond 4,000 km/s)")
    for v in ("published", "corrected"):
        b = res["hubble"][v]["baseline"]
        print(f"  {v:<10} H0 = {b['H0']:.2f} +- {b['H0_stat']:.2f} (stat) +- {res['hubble'][v]['zero_point']['H0_zp_err']:.2f} "
              f"(zero point); {b['n']:,} groups, rms {b['rms_logH']:.4f} dex   [Tully et al. 2023: 74.6]")
    zp = res["hubble"]["published"]["zero_point"]
    print(f"  table zero point {zp['offset_mag']:+.3f} mag vs paper's Eq. 8 {zp['paper_eq8_offset']:+.3f}: "
          f"H0 {zp['H0_if_paper_offset']:.2f} with the paper's; Planck's H0 needs {zp['shift_to_planck_mag']:+.3f} mag")
    e = res["estimators"]
    print("\nCF4's published peculiar velocities")
    print(f"  Vpds / Vpwf / Vpec reproduced to {e['Vpds_reproduced']['median_abs_diff']:.2f} / "
          f"{e['Vpwf_reproduced']['median_abs_diff']:.2f} / {e['Vpec_reproduced']['median_abs_diff']:.2f} km/s "
          "(Vpwf, Vpec need a base-10 log)")
    r = e["ratio_to_direct_small_v"]
    print(f"  ratio to the direct estimate: CF4 column {r['cf4_column']:.3f}, natural log {r['natural_log']:.3f}")
    w = e["well_measured_beyond_ramp"]
    print(f"  {w['n']} groups with e(DM) <= 0.10: rms CF4 column {w['cf4_Vpec_column']['std']:.0f}, natural log "
          f"{w['wf_natural_log']['std']:.0f}, distance errors alone {w['expected_error_median']:.0f} km/s")
    s = res["sdss_sign"]
    ev = s["evidence"]
    print("\nSDSS Fundamental Plane distances in CF4")
    print(f"  {ev['n_matched']:,} galaxies: DM - DM_z = {ev['cf4_vs_eta']['slope']:+.3f} eta + {ev['cf4_vs_eta']['intercept']:.3f} "
          f"(r = {ev['cf4_vs_eta']['r']:.4f}); eta vs FP size residual slope {ev['eta_vs_fp_residual']['slope']:.2f}")
    print(f"  corrected SDSS zero point {s['zero_point']['dzp']:+.3f} +- {s['zero_point']['dzp_err']:.3f} mag "
          f"({s['zero_point']['corrected_uncalibrated']['n_groups']} shared groups); {s['n_groups_affected']:,} groups rebuilt")
    nt = res["nusser_tully_comparison"]
    print(f"  Nusser & Tully 2026 components reproduced: rms {nt['corrected']['rms_diff_sigma']:.2f} sigma (corrected), "
          f"{nt['published']['rms_diff_sigma']:.2f} sigma (as published)")
    print("\nBulk flows (Nusser & Tully bins; joint monopole + dipole fits)")
    for version in ("corrected", "published"):
        b = res["bulk_flow"][version]
        print(f"  {version}:")
        for k, d in b["bins"].items():
            print(f"    {k:<8} n={d['n']:6,} |V| {d['Vamp']:4.0f} +- {d['Vamp_err_meas']:3.0f} km/s toward SG "
                  f"({d['dir_sg'][0]:5.1f},{d['dir_sg'][1]:5.1f}); LCDM rms {d['lcdm_rms_V']:4.0f}; "
                  f"p_V {d['pV']:.3g} ({d['pV_sigma']:.2f} s), p_mV {d['pmV']:.3g} ({d['pmV_sigma']:.2f} s)")
        for k, gl in b["global"].items():
            print(f"    all bins {k}: p_V {gl['pV']:.3g} ({gl['pV_sigma']:.2f} s), p_mV {gl['pmV']:.3g} ({gl['pmV_sigma']:.2f} s)")
    sp = {round(d["R"]): d for d in res["bulk_flow"]["corrected"]["spheres"]}
    for R in (50, 100, 150, 200, 300):
        if R in sp:
            d = sp[R]
            print(f"  sphere {R:3d}: |V| {d['Vamp']:4.0f} +- {d['Vamp_err_meas']:3.0f} km/s (Reff {d['reff']:.0f}); "
                  f"LCDM 68% band {d['lcdm_band_V'][0]:.0f}-{d['lcdm_band_V'][2]:.0f}; p_V {d['pV']:.3g}")
    lf = res["luminosity_function"]
    print(f"\n2MRS luminosity function: alpha {lf['alpha']:.3f} +- {lf['alpha_err']:.3f}, M* - 5 log h "
          f"{lf['Mstar']:.3f} +- {lf['Mstar_err']:.3f}")
    fd = res["flux_dipole"]
    for fr_ in ("lg", "cmb"):
        b60, b130 = fd[f"beta_{fr_}_60"], fd[f"beta_{fr_}_130"]
        print(f"  flux dipole ({fr_} frame): misalignment {b60['misalignment']:.1f} deg at 60, {b130['misalignment']:.1f} "
              f"at 130 h^-1 Mpc; beta_L {b60['beta_L']:.2f} (60), {b130['beta_L']:.2f} (130)")
    print(f"  luminosity density {fd['rho_L_isophotal']:.3g} h Lsun/Mpc^3 (isophotal); ~{fd['rho_L_total_equiv']:.3g} total "
          f"[Erdogdu 2006: 7.67e8]")
    vc = res["velocity_comparison"]
    print(f"\nGravity vs flows ({vc['n_groups']:,} CF4 groups): beta* = {vc['beta']:.3f} +- {vc['beta_err']:.3f} "
          f"(free amplitude {vc['beta_free_fit']:.3f}); V_ext {vc['Vext']:.0f} km/s toward ({vc['Vext_l']:.0f}, {vc['Vext_b']:.0f}); "
          f"sigma_v {vc['sigma_v']:.0f}; H0 {vc['H0']:.2f}   [Carrick+2015: 0.431 +- 0.021, 159 toward (304, 6)]")
    print(f"  predicted Local Group motion: {vc['v_lg_pred_speed']:.0f} km/s toward ({vc['v_lg_pred_lb'][0]:.1f}, "
          f"{vc['v_lg_pred_lb'][1]:.1f}), {vc['v_lg_misalignment']:.1f} deg from the CMB's (internal part "
          f"{vc['v_lg_internal_speed']:.0f} km/s)")
    for d in res["lg_decomposition"]:
        print(f"    {str(d['shell']):<10} {d['speed']:5.0f} km/s toward ({d['l']:5.1f}, {d['b']:5.1f})")
    print("\nBasins of attraction (2MRS internal field)")
    for k, rep in res["basins"].items():
        members = [m for m, inn in rep["members"].items() if inn]
        print(f"  smoothing {k:<6}: our basin {rep['ours_equiv_diameter']:5.0f} h^-1 Mpc across, "
              f"{rep['ours_mass_msun']:.2g} Msun; attractor at r = {rep['ours_attractor_r']:.0f} toward "
              f"(l, b) = ({rep['ours_attractor_lb'][0]:.0f}, {rep['ours_attractor_lb'][1]:.0f}); contains {', '.join(members)}")


if __name__ == "__main__":
    sys.exit(main())
