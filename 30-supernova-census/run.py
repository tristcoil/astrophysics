"""Run the whole project end to end.

    python run.py                 # fetch what is missing, analyse, draw, summarise
    python run.py --no-figures    # numbers only
    python run.py --figures-only  # redraw from the cached analysis
    python run.py --refresh-lc    # refit every Type II light curve

Downloads whatever is not cached (fetch_data.py), runs the analysis, writes
every figure and the offline explorer, and prints each measurement next to the
published value it is graded against ([CHECK] lines). With a warm cache it
runs offline in about two minutes.
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
import time
import warnings

warnings.filterwarnings("ignore")

import numpy as np

import config as C
import literature as lit

PICKLE = C.DERIVED_DIR / "results.pkl"


def _ci(d, scale=1e4, nd=2):
    return f"{d['rate'] / scale:.{nd}f} (95%: {d['lo'] / scale:.{nd}f}-{d['hi'] / scale:.{nd}f})"


def print_summary(res: dict) -> None:
    say = print
    line = "=" * 78
    say(f"\n{line}\nRESULTS\n{line}")
    b = res["bts_check"]
    cen = res["census"]
    say(f"\nBTS catalogue: {cen['n_candidates_all']:,} candidates, peaks {cen['peak_date_range'][0]} to "
        f"{cen['peak_date_range'][1]}; {cen['n_sn_like_classified_all']:,} classified supernovae/TDEs")
    say(f"  statistical sample (quality + purity, m <= 18.5, Jun 2018 - Oct 2025): {cen['n_stat_sample']:,} "
        f"events, {100 * cen['frac_stat_classified']:.1f}% classified")
    say(f"  [CHECK] Perley et al. 2020: {100 * lit.PERLEY_2020['completeness_m18p5']:.0f}% classified")
    say(f"  our absolute magnitudes in the BTS convention vs the explorer's column: median "
        f"{b['median']:+.4f} mag, 99% within {b['p99_abs']:.3f} mag ({b['n']:,} events)")
    d = res["distances"]
    say(f"  Cosmicflows-4 distances for {d['n_cf4']} of {d['n_below_cz']} supernovae below cz = 4,000 km/s; "
        f"median |change| {d['dM_abs_median']:.2f} mag; SN 2023ixf: M = {d['sn2023ixf'][0]['M_bts']:.2f} "
        f"(redshift) -> {d['sn2023ixf'][0]['M_peak']:.2f} (Cepheid/TRGB distance)")

    r = res["rates"]
    rep = r["perley_replication"]
    say("\nRates (supernovae per Gpc^3 per year, units of 10^4)")
    say(f"  Perley et al.'s method, their window, today's classifications: Ia {_ci(rep['Ia'])}, "
        f"CC {_ci(rep['CC'])}, SLSN (M < -21) {rep['SLSN']['rate']:.1f}")
    say(f"  [CHECK] Perley et al. 2020: Ia {lit.PERLEY_2020['rate_Ia'][0] / 1e4:.2f} +/- "
        f"{lit.PERLEY_2020['rate_Ia'][1] / 1e4:.2f}, CC {lit.PERLEY_2020['rate_CC'][0] / 1e4:.1f} "
        f"(+{lit.PERLEY_2020['rate_CC'][1] / 1e4:.1f}/-{lit.PERLEY_2020['rate_CC'][2] / 1e4:.1f}), "
        f"SLSN {lit.PERLEY_2020['rate_SLSN'][0]}")
    fr = rep["f_rec_today"]
    say(f"  their recovery fraction re-measured today: {fr['n_pass']}/{fr['n_easy']} = {fr['f_rec']:.3f} "
        f"[CHECK] {lit.PERLEY_2020['f_rec'][0]}/{lit.PERLEY_2020['f_rec'][1]}")
    fe = r["f_ext"]
    say(f"  dust loss from the SFD map: f_ext = {fe['f_ext_r']:.3f} (r), {fe['f_ext_g']:.3f} (g) over "
        f"{fe['area_deg2_low_extinction']:,.0f} deg2 [CHECK] Perley et al.: {C.PERLEY_F_EXT}")
    say("  season by season (f_sky fixed at Perley's 0.35):")
    for _, s in r["seasonal_check"].iterrows():
        say(f"    {s['season']}: T = {s['T_eff_yr']:.2f} yr, f_rec = {s['f_rec']:.2f}, f_cl(18.4) = "
            f"{s['f_cl_18p5']:.2f}; Ia {s['rate_Ia'] / 1e4:.2f}, CC {s['rate_CC'] / 1e4:.1f}")
    for key in ("main", "full"):
        w = r[key]
        rr = w["rates"]
        say(f"  {w['label']}: exposure {w['exposure_yr']:.2f} yr (T = {w['T_eff_yr']:.2f}, mean f_rec "
            f"{w['f_rec_mean']:.2f})")
        say(f"    Ia (M < -16.5) {_ci(rr['Ia'])} [{rr['Ia']['n']}]; CC (M < -14) {_ci(rr['CC'])} [{rr['CC']['n']}]")
        say(f"    II {_ci(rr['II'])}, IIb {_ci(rr['IIb'], nd=3)}, IIn {_ci(rr['IIn'], nd=3)}, "
            f"Ib/c {_ci(rr['Ibc'])}, Ic-BL {_ci(rr['IcBL'], nd=3)}")
        say(f"    SLSN (M < -21) {rr['SLSN']['rate']:.1f} ({rr['SLSN']['lo']:.1f}-{rr['SLSN']['hi']:.1f}); "
            f"SLSN-I {rr['SLSN-I']['rate']:.1f}; TDE {rr['TDE']['rate']:.0f}")
        rat = w["ratio_cc_ia"]
        say(f"    CC / Ia (coverage-free): {rat['ratio']:.2f} (68%: {rat['lo68']:.2f}-{rat['hi68']:.2f})")
        f14 = w["fractions"]["M<-14.0"]
        say(f"    volumetric shares (M < -14): Ia {100 * f14['Ia']['frac']:.1f}%, Ib/c {100 * f14['Ibc']['frac']:.1f}%, "
            f"II {100 * f14['II']['frac']:.1f}% (IIb {100 * f14['IIb_of_II']:.1f}% and IIn "
            f"{100 * f14['IIn_of_II']:.1f}% of II)")
    say(f"  [CHECK] Li et al. 2011 volume-limited: Ia {lit.LI_2011_VOLUME['Ia'][0]}%, Ib/c "
        f"{lit.LI_2011_VOLUME['Ibc'][0]}%, II {lit.LI_2011_VOLUME['II'][0]}%; IIb 11.9%, IIn 8.6% of II")
    ml = r["mag_limited"].set_index("m_lim").loc[18.5]
    tot = ml["Ia"] + ml["Ibc"] + ml["II"] + ml["IIn"]
    say(f"  magnitude-limited (m < 18.5): Ia {100 * ml['Ia'] / tot:.1f}%, Ib/c {100 * ml['Ibc'] / tot:.1f}%, "
        f"II+IIn {100 * (ml['II'] + ml['IIn']) / tot:.1f}%  [CHECK] Li et al. 2011 mag-limited: 79.2 / 4.1 / 16.6%")
    say(f"  [CHECK] other surveys, Ia: PTF {lit.FROHMAIER_2019_IA[0] / 1e4:.2f} +/- {lit.FROHMAIER_2019_IA[1] / 1e4:.2f}; "
        f"ASAS-SN {lit.DESAI_2026_IA[0] / 1e4:.2f} +/- {lit.DESAI_2026_IA[1] / 1e4:.2f}; "
        f"CC: PTF {lit.FROHMAIER_2021['rate_CC'][0] / 1e4:.1f}, ASAS-SN {lit.PESSI_2025['rate_CC'][0] / 1e4:.1f}")
    for k, v in r["cc_anchored"].items():
        say(f"  CC rate through our CC/Ia ratio and the {k} Ia rate: {v['rate'] / 1e4:.1f} +/- {v['err68'] / 1e4:.1f}")
    zr = r["main_redshift_distances"]
    say(f"  with redshift distances only (no Cosmicflows-4): Ia {zr['Ia']['rate'] / 1e4:.2f}, CC {zr['CC']['rate'] / 1e4:.2f}")
    cs = r["cc_vs_sfr"]
    say(f"  core-collapse supernovae per 1,000 M_sun of stars formed: {1e3 * cs['k_measured_j3']:.1f}-"
        f"{1e3 * cs['k_measured_md14']:.1f} (J3 / Madau & Dickinson SFR at z = {cs['z_median_cc']:.3f}); "
        f"Kroupa IMF predicts {1e3 * cs['k_kroupa']['8-18']:.1f} (8-18 M_sun), {1e3 * cs['k_kroupa']['8-25']:.1f} "
        f"(8-25), {1e3 * cs['k_kroupa']['8-50']:.1f} (8-50)")

    h = res["hosts"]
    say(f"\nHosts: {h['n_with_host']:,} supernovae with an SDSS host, {h['n_with_mpa']:,} with MPA-JHU "
        f"measurements ({100 * h['frac_spec_z_agrees']:.1f}% of host spectra agree in redshift); "
        f"{h['n_plane']:,} in the SFR-mass analysis")
    for cut in ("-11.0", "-12.0"):
        q = h["quenched"][cut]
        say(f"  hosts with log sSFR < {cut}: " + "; ".join(
            f"{fam} {100 * v['obs']:.1f}% (pred. by mass {100 * v['pred_mass']:.0f}%, by SFR "
            f"{100 * v['pred_sfr']:.1f}%)" for fam, v in q.items()))
    pm = h["ia_per_mass_sf_over_q"]
    say(f"  SNe Ia per unit stellar mass, star-forming / quenched galaxies: {pm['ratio']:.2f} +/- {pm['err']:.2f}")
    for key, lab in (("per_mass_Ia", "Type Ia"), ("per_mass_CC", "core collapse")):
        t = h[key]
        say(f"  {lab} per unit stellar mass vs host sSFR (relative to the average galaxy): " + ", ".join(
            f"{r.ssfr_mid:.1f}: {r.rel:.2f}" for r in t.itertuples() if r.n > 0))
    say("  [CHECK] Mannucci et al. 2005: Type Ia per unit mass ~20x higher in late-type galaxies than in E/S0")
    ab = h["a_plus_b"]
    say(f"  A + B fit: log(B/A) = {ab['log_B_over_A_yr']:.2f} ({ab['lo68']:.2f}-{ab['hi68']:.2f}) yr  "
        f"[CHECK] Sullivan et al. 2006: {np.log10(lit.SULLIVAN_2006_AB['B'][0] / lit.SULLIVAN_2006_AB['A'][0]):.2f}")
    t2 = h["perley_table2"]
    say("  Perley et al. Table 2, today's sample (Ia / CC): " + "; ".join(
        f"{k} {100 * t2['Ia'].get(k, {'frac': 0})['frac']:.1f} / {100 * t2['CC'].get(k, {'frac': 0})['frac']:.1f}%"
        for k in ("red sequence", "green valley", "blue cloud", "subdwarf")))
    say("  [CHECK] Perley et al. 2020: red sequence 31 / 11%, green valley 19 / 19%, blue cloud 47 / 64%, "
        "subdwarf 1.0 / 5.0% (+ hostless 2.2 / 1.5%)")

    lc = res["lightcurves"]["stats"]
    say(f"\nType II light curves: {lc['n_lc']:,}; explosion epoch for {lc['n_explosion_epoch']:,}; fall off the "
        f"plateau seen in {lc['n_fall_seen']}; complete for {lc['n_full']} ({lc['n_ni_tail']} with a tail)")
    say(f"  plateau end t_PT: median {lc['t_pt_median']:.1f} d, mean {lc['t_pt_mean']:.1f} +/- {lc['t_pt_std']:.1f}; "
        f"volume-weighted median {lc['t_pt_vweighted_median']:.1f} d")
    say(f"  start of the fall (t_PT - 2 w): mean {lc['t_end_mean']:.1f} +/- {lc['t_end_std']:.1f} d  "
        f"[CHECK] Anderson et al. 2014 OPTd: {lit.ANDERSON_2014_OPTD[0]} +/- {lit.ANDERSON_2014_OPTD[1]} d")
    say(f"  M_r at day 50: median {lc['Mr50_median']:.2f} (volume-weighted {lc['Mr50_vweighted_median']:.2f})")
    say(f"  nickel from the tail: median {lc['ni_tail_median']:.3f} M_sun [CHECK] Anderson 2019 median "
        f"{lit.ANDERSON_2019_NI_MEDIAN} (literature, not magnitude-limited); M_Ni ~ L50^{lc['ni_vs_L50_slope']:.2f}")
    say(f"  R = 500 R_sun: E = {lc['E51_median']:.2f} x 10^51 erg (68%: {lc['E51_q'][0]:.2f}-{lc['E51_q'][1]:.2f}), "
        f"M_ej = {lc['Mej_median']:.1f} M_sun ({lc['Mej_q'][0]:.1f}-{lc['Mej_q'][1]:.1f}); for R = 300 / 1000: "
        f"M_ej {lc['Mej_median_Rlo']:.1f} / {lc['Mej_median_Rhi']:.1f}, E {lc['E51_median_Rlo']:.2f} / "
        f"{lc['E51_median_Rhi']:.2f}")
    mz = lit.MARTINEZ_2022
    say(f"  [CHECK] Martinez et al. 2022 (hydrodynamics + velocities): M_ej {mz['M_ej_range'][0]}-{mz['M_ej_range'][1]}, "
        f"E {mz['E51_range'][0]}-{mz['E51_range'][1]}")

    p = res["progenitors"]
    c = p["counting"]
    say(f"\nProgenitors: Smartt's counting: {100 * c['f_8_18_of_8_100']:.1f}% of 8-100 M_sun stars below 18 M_sun "
        f"-> {c['n_expected_above']:.1f} expected above (found <= 1), Poisson p = {c['p_at_most_found']:.1e}")
    say(f"  [CHECK] Smartt 2015: 70%, 13 expected, p = {lit.SMARTT_2015['p']:.0e}")
    for key, lab, chk in (("ml_s15_stars", "Smartt 2015 masses (STARS)", lit.SMARTT_2015["M_hi_STARS"]),
                          ("ml_s15_kepler", "Smartt 2015 masses (KEPLER)", lit.SMARTT_2015["M_hi_KEPLER"]),
                          ("ml_db18", "Davies & Beasor 2018 masses", lit.DAVIES_BEASOR_2018["M_hi"])):
        f = p[key]
        say(f"  {lab}: M_lo {f['M_lo']:.1f}; M_hi from detections {f['M_hi_detections_only']:.1f}, with limits "
            f"{f['M_hi']:.1f} (68% {f['M_hi_68'][0]:.1f}-{f['M_hi_68'][1]:.1f}, 95% {f['M_hi_95'][0]:.1f}-"
            f"{f['M_hi_95'][1]:.1f}); Delta chi2 at 25 M_sun = {f['delta_chi2_Mhi_25']:.1f}  [CHECK] {chk}")
    lf, lfb = p["lum_fit"], p["lum_fit_beasor"]
    say(f"  luminosity space (slope {p['gamma_L']:.3f}): log L_hi = {lf['logL_hi']:.2f} (68% {lf['logL_hi_68'][0]:.2f}-"
        f"{lf['logL_hi_68'][1]:.2f}); log L = 5.5 disfavoured at {lf['sigma']:.1f} sigma "
        f"[CHECK] Davies & Beasor 2020: 5.28 (+0.12/-0.06), ~2 sigma")
    say(f"    with Beasor et al.'s extra 0.2 dex: log L_hi = {lfb['logL_hi']:.2f}; 5.5 at {lfb['sigma']:.1f} sigma")
    nb, nbb = p["nothing_brighter"], p["nothing_brighter_beasor"]
    say(f"  a sample of {nb['n']} shows nothing brighter than log L = {nb['brightest_seen']:.2f} if stars up to 5.5 "
        f"exploded: {100 * nb['p']:.1f}% ({nb['sigma']:.1f} sigma); with the extra scatter {100 * nbb['p']:.1f}%")


def scalars(obj):
    """The JSON-friendly part of the results (numbers, short lists, dicts)."""
    import pandas as pd
    if isinstance(obj, dict):
        return {str(k): scalars(v) for k, v in obj.items()
                if not isinstance(v, (pd.DataFrame, pd.Series, np.ndarray)) or k in ()}
    if isinstance(obj, (list, tuple)):
        return [scalars(v) for v in obj]
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, (float, int, str, bool)) or obj is None:
        return obj
    return str(obj)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument("--figures-only", action="store_true")
    ap.add_argument("--refresh-lc", action="store_true")
    ap.add_argument("--no-fetch", action="store_true", help="never touch the network")
    args = ap.parse_args(argv)
    t0 = time.time()
    import bts
    if args.figures_only:
        res = pickle.loads(PICKLE.read_bytes())
        df = bts.load()
    else:
        if not args.no_fetch:
            import fetch_data
            fetch_data.main(["--only", "bts", "alerce", "hosts", "background", "cf4", "dust", "papers"])
        import analysis
        res, df = analysis.run_all(refresh_lc=args.refresh_lc)
        PICKLE.write_bytes(pickle.dumps(res))
        C.RESULTS_PATH.write_text(json.dumps(scalars(res), indent=1, default=str))
        print_summary(res)
    if not args.no_figures:
        import figures
        import explorer
        print("\nFigures:")
        figures.draw_all(res, df)
        print(f"  {explorer.write(res, df)}")
    print(f"\nDone in {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
