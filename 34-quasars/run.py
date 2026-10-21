"""Run the whole project end to end.

    python run.py                 # fetch what is missing, analyse, draw, summarise
    python run.py --figures-only  # redraw from the cached analysis
    python run.py --no-figures    # numbers only

Downloads whatever is not cached (see fetch_data.py), runs the analysis, writes
every figure and the offline explorer, and prints each measurement next to the
published value it is graded against. Once the cache is warm it runs offline in
about two minutes.
"""

from __future__ import annotations

import argparse
import pickle
import sys
import time
import warnings

warnings.filterwarnings("ignore")

import literature as lit
from config import (DR16Q_TABLE, DR7_TABLE, RESULTS_PICKLE, SPEC_LIST, WUSHEN_TABLE)


def need_fetch() -> bool:
    return not all(p.exists() for p in (DR16Q_TABLE, DR7_TABLE, WUSHEN_TABLE, SPEC_LIST))


def print_summary(res: dict) -> None:
    say = print
    line = "=" * 78
    c, ph, q, m, g = res["composite"], res["photometry"], res["qlf"], res["masses"], res["growth"]
    say(f"\n{line}\nRESULTS\n{line}")
    cz = res["census"]
    say(f"\nDR16Q: {cz['n']:,} quasars; z > 5: {cz['n_z_gt'][5.0]:,}, z > 6: {cz['n_z_gt'][6.0]:,}; "
        f"masses for {m['n']:,}")
    for r in cz["by_programme"]:
        say(f"  {r['programme']:<17s} {r['n']:>8,}   median z {r['median']:.2f}")

    say(f"\nComposite spectrum: {c['n_spectra']:,} spectra, {c['wave_min']:.0f}-{c['wave_max']:.0f} A, "
        f"S/N {c['snr_2800']:.0f} per A at Mg II")
    say(f"  alpha_nu (geometric mean, 1350/4215 A): {c['alpha_gm']:.2f} all, {c['alpha_gm_luminous']:.2f} for "
        f"M_i < -25   [CHECK] Vanden Berk 2001: {lit.VB01['alpha_nu_gm']}")
    say(f"  alpha_nu (median): {c['alpha_median']:.2f} all, {c['alpha_median_luminous']:.2f} luminous   "
        f"[CHECK] VB01: {lit.VB01['alpha_nu_median']}")
    for k, v in c["vb_comparison"].items():
        say(f"  ours / VB01 over {k} A: {v['median_ratio']:.3f}")
    say("  C IV equivalent width by luminosity (bright -> faint): "
        + ", ".join(f"{b['civ_ew']:.0f} A" for b in c["luminosity_bins"] if b["civ_ew"] == b["civ_ew"]))

    say("\nColours and K-corrections:")
    for k, v in ph["agreement"].items():
        say(f"  |synthetic - observed| {k}: {v:.3f} mag")
    for k, v in ph["k_comparison"].items():
        say(f"  K(ours) - K(Richards+06), z {k}: median {v['median']:+.3f}, max {v['max_abs']:.3f}")

    r6, sk, sk0, d37 = (q["r06_reproduction"], q["sk12_comparison"], q["sk12_comparison_without_image_factor"],
                        q["dr3_vs_dr7"])
    say("\nLuminosity function:")
    say(f"  Richards+06 rebuilt from their own 15,343 quasars: {r6['n_matched']} bins, median "
        f"{r6['median_dlog']:+.4f} dex, {100 * r6['frac_within_0p05']:.0f}% within 0.05 dex")
    say(f"  DR7 uniform ({q['n_dr7_uniform']:,} quasars) vs Shen & Kelly 2012: median {sk['median_dlog']:+.4f} dex "
        f"(their 5% image factor removed: {sk0['median_dlog']:+.4f}, scatter {sk0['robust_std_dlog']:.4f})")
    say(f"  DR7 vs DR3 (same bins): median {d37['median_dlog']:+.3f} dex, scatter {d37['robust_std_dlog']:.3f}")
    pk = q["era_peak"]
    say(f"  luminous quasars (M_i < -27.6) peak at z = {pk['z_peak']:.2f} "
        f"(68%: {pk['z_peak_p16']:.2f}-{pk['z_peak_p84']:.2f}), {pk['age_at_peak_gyr']:.2f} Gyr after the Big Bang, "
        f"{pk['density_peak_per_gpc3']:.0f} per Gpc^3   [CHECK] Richards+06: 'between z = 2 and 3'")
    for k, v in pk["variants"].items():
        say(f"    variant {k:<34s} z_peak = {v:.2f}")

    say("\nBlack holes:")
    kol = m["kollmeier"]
    say(f"  Eddington ratios (0.3<z<4, 1e45<Lbol<1e47, {kol['n']:,}): median {10 ** kol['median']:.3f}, "
        f"spread {kol['robust_std']:.2f} dex, {100 * kol['frac_above_1']:.1f}% above 1   "
        f"[CHECK] Kollmeier+06: peak {lit.KOLLMEIER06['peak']}, 0.3 dex")
    say(f"  Mg II - H beta: {m['hb_mg']['median']:+.3f} dex (scatter {m['hb_mg']['robust_std']:.2f}, n = {m['hb_mg']['n']:,})")
    say(f"  C IV - Mg II:   {m['mg_c4']['median']:+.3f} dex (scatter {m['mg_c4']['robust_std']:.2f}, n = {m['mg_c4']['n']:,})")
    bs = m["blueshift"]
    say(f"  C IV bias vs blueshift: {bs['slope_raw']:+.3f} dex per 1000 km/s raw, {bs['slope_c17']:+.3f} with "
        f"Coatman+17, {bs['slope_ours']:+.3f} with our cross-validated recalibration")
    st = m["z2_step"]
    s1, s2 = st.iloc[3], st.iloc[4]
    say(f"  the z = 2 step: fiducial median {s1['fiducial']:.2f} -> {s2['fiducial']:.2f}; Mg II {s1['mgii']:.2f} -> "
        f"{s2['mgii']:.2f}; recalibrated C IV {s1['civ_recal']:.2f} -> {s2['civ_recal']:.2f}")

    say("\nThe growth problem:")
    j = g["j0313"]
    say(f"  J0313-1806: {j['age_myr']:.0f} Myr after the Big Bang [CHECK] Wang+21: ~670; mass with Wang's width "
        f"(VO09) 10^{j['log_mbh_vo09_from_wang_width']:.2f} [CHECK] 10^{j['log_mbh_published']:.2f}; "
        f"with Yang's width (S11) 10^{g['j0313_s11']:.2f}")
    say("  minimum seed (Wang's mass), growth from z = 30/20/15: "
        + " / ".join(f"{j['seeds'][k]['seed_published']:.1e}" for k in (30, 20, 15))
        + "   [CHECK] Wang+21: 1e4-1e5 for z = 15-30")
    f = g["frac_seed_above_100"]
    say(f"  {len(g['hz'])} z > 5.7 quasars needing a seed > 100 Msun: {100 * f[30]:.0f}% (from z = 30), "
        f"{100 * f[20]:.0f}% (from z = 20); e-folding time {g['efold_myr']:.1f} Myr; 100 -> 1e9 Msun takes "
        f"{g['t_100_to_1e9_myr']:.0f} Myr")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--figures-only", action="store_true")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args(argv)
    t0 = time.time()
    if args.figures_only:
        with open(RESULTS_PICKLE, "rb") as fh:
            res = pickle.load(fh)
    else:
        if need_fetch():
            import fetch_data
            fetch_data.main([])
        import analysis
        res = analysis.run()
    if not args.no_figures:
        import explorer
        import figures
        figures.main([])
        p = explorer.build(res)
        print(f"explorer: {p}")
    print_summary(res)
    print(f"\nDone in {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
