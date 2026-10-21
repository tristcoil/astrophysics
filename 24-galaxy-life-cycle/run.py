"""Run the whole project end to end.

    python run.py              # fetch what is missing, analyse, draw, summarise
    python run.py --no-figures # numbers only
    python run.py --figures-only  # redraw from the cached analysis

Downloads whatever is not already cached (see fetch_data.py), runs the
analysis, writes every figure, and prints each measurement next to the
published value it is graded against. Once the cache is warm it runs offline
in a few minutes.
"""

from __future__ import annotations

import argparse
import pickle
import sys
import time
import warnings

warnings.filterwarnings("ignore")

import analysis
import literature as lit
from config import GALAXY_CATALOG, TARGET_CATALOG


def _pm(d, k, fmt="{:.3f}"):
    b = d[k]
    return (fmt.format(b["p50"]) + " (68%: " + fmt.format(b["p16"]) + "–"
            + fmt.format(b["p84"]) + ")")


def print_summary(res: dict) -> None:
    line = "=" * 78
    say = print
    say(f"\n{line}\nRESULTS\n{line}")

    geo = res["geometry"]
    say(f"\nCatalogue: {res['n_galaxies']:,} Main Galaxy Sample galaxies; MPA-JHU entry for "
        f"{100 * res['mpa_coverage']['has_mpa']:.2f}%, valid mass + SFR for "
        f"{100 * res['mpa_coverage']['sfr_ok']:.2f}% ({res['mpa_coverage']['n_ssfr_above_minus8']} "
        f"with log sSFR > -8 dropped as unphysical)")
    say(f"Targets: {geo['n_targets_unique']:,} ({geo['n_targets_unobserved']:,} never observed); "
        f"fibre completeness {100 * geo['fibre_completeness_all']:.1f}% "
        f"(Strauss et al.: ~{100 * lit.MGS_FIBRE_LOSS:.0f}% lost to collisions alone)")
    wv = geo["window_vmax"]
    say(f"  V_max window 14.5 <= r < 17.6: completeness {100 * wv['fibre_completeness']:.1f}%; "
        f"of observed targets {100 * wv['frac_observed_not_in_sample']:.1f}% are not in the sample: "
        f"stars {100 * wv['frac_stars']:.1f}%, quasars {100 * wv['frac_qso']:.1f}%, failed "
        f"redshifts {100 * wv['frac_bad_z']:.1f}%, other {100 * wv['frac_missing_other']:.1f}%")
    say(f"  completeness weight scale: {geo['completeness_weight_scale']:.4f}")
    say(f"  area: {geo['area_deg2']:,.1f} deg2 of tiled sectors (DR7 Legacy spectroscopy: 8,032); "
        f"{100 * geo['frac_targets_outside_footprint']:.2f}% of targets fall outside; "
        f"cross-check from target density in a well-filled interior: "
        f"{geo['area_from_target_density_deg2']:,.0f} deg2 (biased low, see README)")

    ms = res["ms"]
    r = ms["ridge"]
    rp = lit.RENZINI_PENG_2015
    say(f"\nMain sequence (ridge line, AGN excluded, {ms['n']:,} galaxies at 0.02 < z < 0.085):")
    say(f"  log SFR = {r['slope']:.3f} (+/- {r['slope_err']:.3f}) log M "
        f"{r['intercept']:+.2f} (+/- {r['intercept_err']:.2f})")
    say(f"  [CHECK] Renzini & Peng 2015: {rp['slope'][0]} +/- {rp['slope'][1]}, "
        f"{rp['intercept'][0]} +/- {rp['intercept'][1]}")
    say(f"  at log M = 10.5: ours {r['slope'] * 10.5 + r['intercept']:.3f}, theirs "
        f"{rp['slope'][0] * 10.5 + rp['intercept'][0]:.3f}")
    rl = ms["ridge_lowmass"]
    say(f"  clean part only (9 < log M < 10.5): slope {rl['slope']:.3f} +/- {rl['slope_err']:.3f}, "
        f"intercept {rl['intercept']:+.2f}")
    say("  ridge points (ours vs R&P): " + ", ".join(
        f"{m:.1f}: {p:+.2f}/{rp['slope'][0] * m + rp['intercept'][0]:+.2f}"
        for m, p in zip(r["logm"], r["peak"])))
    ra = ms["ridge_with_agn"]
    say(f"  with AGN kept: slope {ra['slope']:.3f}, intercept {ra['intercept']:+.2f}")
    say(f"  SF peak width at log M = 10.5: sigma = {ms['sf_peak_10p5']['sigma']:.3f} dex "
        f"(R&P: {rp['sigma_dex']}); over 9-11: {ms['sf_peak_all']['sigma']:.3f}")
    say(f"  number x SFR ridge sits {ms['sfr_weighted_offset_dex']:+.2f} dex above the number "
        f"ridge (R&P: ~{rp['number_times_sfr_offset_dex']})")
    say(f"  green valley floor: {ms['valley']['position']:+.2f} dex below the sequence "
        f"(density there {100 * ms['valley']['depth']:.1f}% of the SF peak)")
    say("  valley by mass: " + ", ".join(
        f"{v['logm']:.1f}: {v['position']:+.2f}" if v["found"] else f"{v['logm']:.1f}: none"
        for v in ms["valley_by_mass"]))
    say(f"  half of galaxies quenched above log M = {ms['logm_half_quenched']:.2f}")
    for key in ("ms_gswlc", "ms_gswlc_m2"):
        m2 = res[key]
        v2 = m2["valley"]
        say(f"  {m2['label']}: slope {m2['ridge']['slope']:.3f} (9-10.5: "
            f"{m2['ridge_lowmass']['slope']:.3f}), intercept {m2['ridge']['intercept']:+.2f}; "
            f"own valley {'%+.2f' % v2['position'] if v2['found'] else 'none (single peak + tail)'}; "
            f"with the MPA-JHU valley, half quenched at {m2['logm_half_quenched']:.2f}")
    cd = res["color_divider"]
    say(f"\nColour: red if (g-r)_0 > {cd['slope']:.3f} log M {cd['intercept']:+.3f}; "
        f"red {100 * res['red_fraction_local'][0]:.1f}%, quenched "
        f"{100 * res['quenched_fraction_local'][0]:.1f}% (weighted)")
    say(f"  red but still forming stars: {100 * res['red_but_starforming'][0]:.1f}%; "
        f"blue but quenched: {100 * res['blue_but_quenched'][0]:.1f}%")

    e = res["env"]
    t = e["tracer"]
    say(f"\nEnvironment: {t['n']:,} tracers (M_r + Q z < {t['M0_limit']:.2f}), "
        f"n = {t['nbar_mpc3'] * 1e3:.2f} x 1e-3 Mpc^-3; median d5 {e['median_d5_mpc']:.1f} Mpc")
    say(f"  {e['n_usable']:,} of {e['n_local']:,} galaxies usable "
        f"({100 * e['frac_edge_dropped']:.1f}% dropped at survey edges)")
    pg = lit.PENG_2010
    for key, label in (("quench", "sSFR-quenched"), ("red", "colour-red")):
        f = e[key]["fit"]
        say(f"  separable fit ({label}): log p1 = {f['log_p1']:.2f} +/- {f['err'][0]:.2f}, "
            f"p2 = {f['p2']:.2f} +/- {f['err'][1]:.2f}, log p3 = {f['log_p3']:.2f} +/- "
            f"{f['err'][2]:.2f}, p4 = {f['p4']:.2f} +/- {f['err'][3]:.2f}; chi2 {f['chi2']:.0f} "
            f"for {f['dof']} dof; rms residual {100 * e[key]['rms_resid']:.1f}%")
        md = e[key]["mass_dependent_env"]
        say(f"    environment threshold vs mass: d log p1 / d log M = "
            f"{md['slope_logp1_per_dex']:+.3f} +/- {md['slope_err']:.3f} (0 = separable)")
    say(f"  [CHECK] Peng et al. 2010 (colour): log p1 = {pg['log_p1'][0]}, p2 = {pg['p2'][0]}, "
        f"log p3 = {pg['log_p3'][0]}, p4 = {pg['p4'][0]}")

    a = res["agn"]
    say("\nAGN (composite + Seyfert + LINER), weighted fractions at 10 < log M < 11:")
    for lab, d in (("MPA-JHU SFRs", a), ("GSWLC SFRs", res["agn_gswlc"])):
        z = d["zones_logm10_11"]
        say(f"  {lab}: star-forming {100 * z['star_forming'][0]:.1f}%, green valley "
            f"{100 * z['green_valley'][0]:.1f}%, quenched {100 * z['quenched'][0]:.1f}%")
        for row in d["by_dms"]:
            say(f"    log M {row['logm'][0]}-{row['logm'][1]}: AGN fraction peaks at "
                f"{100 * row['peak']:.0f}%, {row['dms_at_peak']:+.1f} dex from the sequence")
    b = res["bpt_agreement"]
    say(f"  our BPT classes vs MPA-JHU's: SF {100 * b['agree_sf']:.1f}%, composite "
        f"{100 * b['agree_composite']:.1f}%, AGN {100 * b['agree_agn']:.1f}% agree")

    c = res["gswlc_cmp"]
    s16 = lit.SALIM_2016
    say(f"\nGSWLC-2 vs MPA-JHU ({c['n']:,} galaxies):")
    say(f"  masses: MPA-JHU - GSWLC = {c['dlogm_median']:+.3f} dex, robust scatter "
        f"{c['dlogm_robust_sigma']:.3f} (Salim: GSWLC heavier by "
        f"{s16['mass_offset_range_dex'][0]}-{s16['mass_offset_range_dex'][1]}, scatter "
        f"{s16['mass_scatter_dex']})")
    for k in ("star_forming", "agn", "weak_or_none"):
        d = c[f"dssfr_{k}"]
        say(f"  sSFR, {k:>12}: median {d['median']:+.3f} dex, robust scatter "
            f"{d['robust_sigma']:.3f}, std {d['std']:.3f} (n = {d['n']:,})")
    say(f"  [CHECK] Salim et al. 2016, star-forming: scatter {s16['ssfr_scatter_sf_dex']} dex")

    cs = res["census"]
    b04 = lit.BRINCHMANN_2004_SFRD
    say(f"\nCensus, 0.02 < z < 0.2 ({cs['n']:,} galaxies, {cs['n_jackknife']} jackknife regions):")
    say(f"  SFR density {cs['rho_sfr']:.5f} at z_eff = {cs['z_eff_sfr']:.3f}; evolved to z = 0.1: "
        f"{cs['rho_sfr_z0.1']:.5f} +/- {cs['rho_sfr_err']:.5f} M_sun/yr/Mpc^3 (Kroupa)")
    say(f"    beta = 0: {cs['rho_sfr_z0.1_beta0']:.5f}, beta = 3: {cs['rho_sfr_z0.1_beta3']:.5f}")
    cm = res["census_mean"]
    say(f"    likelihood mean instead of median: {cm['rho_sfr_z0.1']:.5f} +/- "
        f"{cm['rho_sfr_err']:.5f} (Salpeter {cm['rho_sfr_salpeter']:.4f})")
    zm = res["census_mean_local"]["zones"]
    say("    mean-based shares (0.02 < z < 0.085): " + "; ".join(
        f"{k}: {100 * v['sfr_share']:.0f}%" for k, v in zm.items()))
    say(f"  [CHECK] Brinchmann et al. 2004: {b04['value']:.5f} (+{b04['sys'][0]:.4f} "
        f"-{b04['sys'][1]:.4f} sys)")
    from sfh import psi_md
    psi01 = psi_md(0.1, **lit.MD14_FIT)
    say(f"  Salpeter: {cs['rho_sfr_salpeter']:.4f}; [CHECK] Madau & Dickinson curve at z = 0.1: "
        f"{psi01:.4f} (ratio {cs['rho_sfr_salpeter'] / psi01:.2f})")
    gsc = res["census_gswlc"]
    say(f"  GSWLC-X2 (Chabrier): {gsc['rho_sfr_z0.1']:.5f} -> Salpeter {gsc['rho_sfr_salpeter']:.4f} "
        f"(covers {100 * gsc['gswlc_coverage']:.0f}% of the galaxies; missing area is "
        f"represented by matched galaxies)")
    say(f"  stellar mass density {cs['rho_star'] / 1e8:.3f} +/- {cs['rho_star_err'] / 1e8:.3f} "
        f"x 1e8 (Kroupa) -> Salpeter {cs['rho_star_salpeter'] / 1e8:.3f} x 1e8 (log "
        f"{__import__('math').log10(cs['rho_star_salpeter']):.3f})")
    say(f"  [CHECK] Li & White 2009: log {lit.LOCAL_STELLAR_MASS_DENSITY['li_white_2009'][0]}; "
        f"MD14 integral: 5.8 x 1e8; GSWLC: {gsc['rho_star_salpeter'] / 1e8:.3f} x 1e8")
    say(f"  half of today's star formation in galaxies below 10^{cs['logm_half_sfr']:.2f}; "
        f"half of the stellar mass above 10^{cs['logm_half_mass']:.2f}")
    z = res["census_local"]["zones"]
    say("  shares (0.02 < z < 0.085): " + "; ".join(
        f"{k}: {100 * v['sfr_share']:.0f}% of SFR, {100 * v['mass_share']:.0f}% of mass"
        for k, v in z.items()))

    h = res["sfh"]
    pub, ours = h["published"], h["ours"]
    sp, so = pub["summary"], ours["summary"]
    say("\nCosmic star-formation history (Madau & Dickinson Table 1, 51 points):")
    say(f"  our fit: A = {ours['A']:.4f}, B = {ours['B']:.2f}, C = {ours['C']:.2f}, "
        f"D = {ours['D']:.2f}; chi2 = {ours['chi2']:.1f} for {ours['dof']} dof "
        f"(their curve: chi2 = {pub['chi2']:.1f})")
    bs = ours["boot"]
    say(f"  peak: z = {_pm(bs, 'z_peak', '{:.2f}')}, {_pm(bs, 't_peak_gyr', '{:.2f}')} Gyr "
        f"[their curve: z = {sp['z_peak']:.3f}, {sp['t_peak_gyr']:.2f} Gyr; text: z ~ 1.9, 3.5 Gyr]")
    say(f"  half of today's stars by z = {_pm(bs, 'z_half_mass', '{:.2f}')} [their curve "
        f"{sp['z_half_mass']:.2f}; text 1.3]")
    say(f"  formed before the peak: {_pm(bs, 'frac_before_peak', '{:.2f}')} [their curve "
        f"{sp['frac_before_peak']:.2f}; text ~0.25]")
    say(f"  last quarter after z = {_pm(bs, 'z_last_quarter', '{:.2f}')} [their curve "
        f"{sp['z_last_quarter']:.2f}; text 0.7]")
    say(f"  before reionization ended (z = 6): {_pm(bs, 'frac_before_reionization_end', '{:.4f}')}"
        f" [their curve {sp['frac_before_reionization_end']:.4f}; text < 0.01]")
    say(f"  rho_*(z=0) = {_pm(bs, 'rho_star_z0', '{:.3g}')} [their curve "
        f"{sp['rho_star_z0']:.3g}; text 5.8e8]")
    say(f"  mass-weighted mean stellar age {_pm(bs, 'mean_stellar_age_gyr', '{:.2f}')} Gyr "
        f"[their curve {sp['mean_stellar_age_gyr']:.2f}; Planck18 "
        f"{pub['summary_planck']['mean_stellar_age_gyr']:.2f}; text 8.3]")
    ew, eu = h["efold_data_weighted"], h["efold_data_unweighted"]
    say(f"  e-folding time since z = 1: data {ew['tau_gyr']:.2f} +/- {ew['tau_err_gyr']:.2f} Gyr "
        f"(weighted), {eu['tau_gyr']:.2f} (unweighted); their curve "
        f"{sp['efold_curve_gyr']:.2f} [text 3.9]")
    hz = h["high_z_slope"]
    say(f"  high-z slope (z > 3): (1+z)^{hz['slope']:.2f} +/- {hz['slope_err']:.2f} [text -2.9]")
    say(f"  psi(7) / psi(0) = {sp['psi_z7_over_z0']:.2f} (their curve) [text: 'about the same']")
    bf = h["baryon_fraction"]
    say(f"  stellar baryon fraction today: {100 * bf['md14_omega_b']:.1f}% (MD14 Omega_b), "
        f"{100 * bf['d3_planck']:.1f}% with D3's Omega_b h^2 = {bf['omega_b_h2_d3']} "
        f"[text 9%; Chabrier {100 * bf['d3_planck_chabrier']:.1f}%, text 5%]")
    say("  direct stellar mass counts minus integral of their curve: " + ", ".join(
        f"z {r['z'][0]}-{r['z'][1]}: {r['median_resid_dex']:+.2f} dex" for r in h["smd_by_z"]))


CACHE = GALAXY_CATALOG.parent / "analysis_cache.pkl"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--no-figures", action="store_true")
    p.add_argument("--figures-only", action="store_true",
                   help="redraw figures from the cached analysis (data/analysis_cache.pkl)")
    args = p.parse_args(argv)
    t0 = time.time()
    if args.figures_only and CACHE.exists():
        with open(CACHE, "rb") as f:
            res = pickle.load(f)
    else:
        if not (GALAXY_CATALOG.exists() and TARGET_CATALOG.exists()):
            import fetch_data
            fetch_data.main([])
        res = analysis.analyse()
        analysis.save_results(res)
        with open(CACHE, "wb") as f:
            pickle.dump(res, f, protocol=pickle.HIGHEST_PROTOCOL)
    print_summary(res)
    if not args.no_figures:
        import figures
        print("\nFigures")
        for fn in (figures.fig_sfr_mass_plane, figures.fig_3d_surface, figures.fig_quenching,
                   figures.fig_agn, figures.fig_methods, figures.fig_cosmic_sfh,
                   figures.fig_mass_buildup, figures.fig_census, figures.fig_gallery,
                   figures.fig_web_slice):
            for path in fn(res):
                if path.endswith(".png"):
                    print(f"  {path}")
        import explorer
        print(f"  {explorer.write(res)}")
    print(f"\nDone in {time.time() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
