"""Run the whole project end to end.

    python run.py                  # download what is missing, analyse, draw, summarise
    python run.py --figures-only   # redraw figures and the explorer from the last analysis
    python run.py --no-fetch       # analyse cached data only (offline)

The first run downloads ~0.75 GB (the ZTF catalogue, 950,000 AAVSO
observations of T CrB, CDS cross-matches with Gaia DR3 and AllWISE, the
Gaia DR3 eclipsing-binary table, Kepler/TESS/NSVS/SuperWASP/ASAS-SN/AAVSO
photometry of KIC 9832227, OGLE photometry of V1309 Sco and nine arXiv
sources); that takes 1-2 hours, mostly the CDS cross-match and the AAVSO API.
Afterwards the analysis and figures run offline in ~10 minutes.
"""

from __future__ import annotations

import argparse
import json
import os
import warnings

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
warnings.filterwarnings("ignore")

import config as C


def print_summary(r: dict) -> None:
    ro, ce, pl, k, v, t = r["roche"], r["census"], r["pl"], r["kic"], r["v1309"], r["tcrb"]
    print("\n" + "=" * 78)
    print("RESULTS: stars that touch, trade mass and merge")
    print("=" * 78)
    km = ro["kic_model"]
    print(f"\nRoche model of KIC 9832227 with Molnar et al.'s parameters: amplitude {km['amp_model_mag']:.3f} mag "
          f"vs Kepler {km['amp_kepler_mag']:.3f} mag (rms {100 * km['rms_diff']:.1f}%)")
    e = ce["edge"]
    print(f"\nZTF contact binaries: {ce['n_ew']:,}; period edge {e['p_edge']:.4f} d (logistic midpoint, "
          f"width {e['width_dex']:.3f} dex); 1% below {ce['quantiles']['0.01']:.4f} d")
    print(f"  edge in r slices: " + ", ".join(f"{s['p_edge']:.4f}" for s in ce["mag_slices"]))
    print(f"  peak of log-P histogram: all {ce['peak_period_logp_hist']:.3f} d, within 1 kpc "
          f"{ce['peak_period_within_1kpc']:.3f} d")
    se = ce["sub_edge"]
    print(f"  P < 0.21 d with parallaxes: {se['n_with_parallax']}, median M_G {se['median_M_G']:.2f}, "
          f"(BP-RP)0 {se['median_bp_rp0']:.2f}; {100 * se['frac_M_G_lt_3']:.1f}% luminous")
    g = ce["gaia"]
    print(f"  Gaia DR3 EBs P<=2d: {g['n_p_le_2d']:,}, shortest {g['p_min']:.5f} d (search floor); contact-like "
          f"{g['n_contact_like']:,}, edge {g['edge_contact_like']['p_edge']:.4f} d")
    zg = ce["ztf_vs_gaia"]
    print(f"  ZTF EW in Gaia's EB table: {zg['n_in_gaia_eb']:,}; same period {100 * zg['same']:.1f}%")
    c18 = pl["chen2018"]
    print(f"\nPL (W1): Chen 2018 reproduced: {c18['reproduce_W1_mean']['slope']:.3f} log P "
          f"{c18['reproduce_W1_mean']['b0']:+.3f}; DR3 recalibration (RUWE<1.4, n={c18['W1_mean_dr3_ruwe']['n']}): "
          f"{c18['W1_mean_dr3_ruwe']['slope']:.3f} ± {c18['W1_mean_dr3_ruwe']['slope_err']:.3f} log P "
          f"{c18['W1_mean_dr3_ruwe']['b0']:+.3f}, scatter {c18['W1_mean_dr3_ruwe']['scatter']:.3f} "
          f"(TGAS, same stars {c18['W1_mean_tgas_same_stars']['scatter']:.3f})")
    zp = pl["ztf_precise"]["W1_ols"]
    print(f"  ZTF sample (n={zp['n']:,}): {zp['slope']:.3f} log P {zp['b0']:+.3f}, scatter {zp['scatter']:.3f}")
    print("  slope by distance: " + ", ".join(f"{b['d0']}-{b['d1']} pc {b['slope']:.2f}" for b in pl["distance_bins"]))
    pt = pl["parallax_test"]
    print(f"  parallax-space test on {pt['n']:,} ZTF EWs: measured/predicted = {pt['ratio']:.4f} ± {pt['ratio_err']:.4f}")
    wv = pl["wesenheit_vs_parallax"]
    print(f"  Wesenheit distances: {wv['n_mapped']:,} mapped, median {wv['median_d_pl'] / 1000:.1f} kpc; "
          f"vs Gaia (n={wv['n']:,}) ratio {wv['median_ratio']:.3f}, scatter {wv['scatter_percent']:.0f}%")
    print(f"\nKIC 9832227: {k['n_timings']} timings ({k['n_ours']} ours), {k['years'][0]:.1f}-{k['years'][1]:.1f}")
    print(f"  NSVS O-C {k['nsvs']['oc_min']:+.1f} ± {k['nsvs']['err_min']:.1f} min; with the label taken literally "
          f"{k['nsvs']['oc_as_labelled_min']:+.1f} min (shift {k['nsvs']['shift_min']:.2f}, expected "
          f"{k['nsvs']['expected_shift_min']:.2f})")
    print(f"  vs Kovacs+19 ({k['vs_kovacs']['n']} pairs): median {k['vs_kovacs']['median_diff_min']:+.2f} min, rms "
          f"{k['vs_kovacs']['rms_diff_min']:.2f} min")
    for key in ("molnar_like_labelled", "molnar_like_correct", "all_correct"):
        p = k["exponential_profiles"][key]
        print(f"  exponential law, {key}: t0 = {p['t0_best']:.2f} ({p['t0_lo']:.2f}-{p['t0_hi']:.2f}), tau "
              f"{p['tau_d']:.3f} d, jitter {p['jitter_min']:.2f} min (smooth period change "
              f"{p['ltte_jitter_min']:.2f} min)")
    ql = k["fits_all"]["quad+ltte"]
    print(f"  all data: steady Pdot = {ql['pdot']:.2e} + {ql['period_yr']:.1f}-yr wave of {ql['amp_min']:.1f} min, "
          f"jitter {ql['jitter_min']:.2f} min; quadratic alone Pdot {k['fits_all']['quadratic']['pdot']:.2e}")
    print(f"  period 2011 -> 2025: {k['period_change']['p_2011']:.8f} -> {k['period_change']['p_2025']:.8f} d "
          f"({k['period_change']['dp_s']:+.2f} s); Kovacs+19 model misses 2024 by {k['kovacs_prediction']['miss_2024_min']:+.0f} min")
    lp, lph, on = v["law_from_periods"], v["law_from_phases"], v["onset"]
    print(f"\nV1309 Sco: {v['n_I']} I-band points; period drop {v['period_drop_percent']:.2f}%")
    print(f"  law from periods: P0 {lp['p0']:.5f}, c {lp['c']:.2f}, t0 {lp['t0_year']:.2f} ± {lp['t0_year_err']:.2f} "
          f"(Tylenda 1.4456, 15.29, 2010.10); from eclipse phases t0 {lph['t0_year']:.2f} ± {lph['t0_year_err']:.2f}")
    print(f"  vs Tylenda's 17 periods: median diff {v['vs_tylenda']['median_diff'] * 86400:.2f} s, chi2 "
          f"{v['vs_tylenda']['chi2']:.1f}/{v['vs_tylenda']['n']}")
    print(f"  rise begins {on['t_on10_year']:.3f} (e-folding {on['tau_d']:.1f} d)")
    er = t["eruptions"]
    print(f"\nT CrB: {t['n_obs']:,} AAVSO observations, {t['first_year']:.2f}-{t['last_year']:.2f}; "
          f"{t['visual_rejected']:,} visual estimates rejected")
    print(f"  1946: t2 {er['1946']['t2']:.1f} d, t3 {er['1946']['t3']:.1f} d; secondary maximum day "
          f"{er['1946']['secondary_peak_day']:.0f} (1866: day {er['1866']['secondary_peak_day']:.0f})")
    print(f"  B: quiescence {t['B_levels']['quiescence_2005_14']:.2f}, high {t['B_levels']['high_2016_23']:.2f}, "
          f"2023-24 dip {t['B_levels']['dip_2023_24']:.2f}, last 90 d {t['B_levels']['last_90d']:.2f}")
    print(f"  last 30 days: V {t['last30']['V']['median']:.2f}, B {t['last30']['B']['median']:.2f}: no eruption")
    hz = t["hazard"]
    print(f"  recurrence hazard: {100 * hz['p_rest_of_2026']:.0f}% rest of 2026, {100 * hz['p_next_12_months']:.0f}% "
          f"next 12 months (check vs Pei+26: {100 * hz['check_rest_of_2026_from_jul11']:.1f}% / "
          f"{100 * hz['check_12_months_from_jul11']:.1f}%)")
    print()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--figures-only", action="store_true")
    ap.add_argument("--no-fetch", action="store_true")
    a = ap.parse_args()
    if not a.figures_only:
        if not a.no_fetch:
            import fetch_data
            fetch_data.fetch_all()
        import analysis
        analysis.run_all(("roche", "kic", "v1309", "census", "pl", "tcrb"))
    import figures
    import explorer
    figures.main()
    print("explorer:", explorer.build())
    print_summary(json.loads(C.RESULTS_PATH.read_text()))


if __name__ == "__main__":
    main()
