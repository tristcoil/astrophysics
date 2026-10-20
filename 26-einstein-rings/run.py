"""Run the whole project end to end.

    python run.py

Fetches whatever is not already cached, fits lens models to every Hubble image
that does not have one yet (cached in data/lens_models/), runs the population
analysis, writes every figure and the interactive explorer, saves
data/results.json, and prints a summary that grades our numbers against the
published ones. Once the cache is warm this runs offline.

    python run.py --refit          # redo all lens models (~1 h on 4 cores)
    python run.py --no-figures     # numbers only
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import warnings

warnings.filterwarnings("ignore")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

import analysis as A
import catalog
import literature
from config import FIGURE_DIR, MODEL_DIR, RESULTS_PATH, SHOWCASE, SPECTRUM_LENS

ANATOMY_LENS = "J1627-0053"
CONSPIRACY_LENS = "J0946+1006"


def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()
                if not isinstance(v, np.ndarray) or v.size < 50}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.bool_):
        return bool(obj)
    if hasattr(obj, "to_dict"):
        return obj.to_dict("records")
    return obj


def comparison(t, ix, iso, slopes, dm, dc, ours, auger, treu):
    """Our headline numbers next to the published ones, for the summary and figure 0.
    Each row: label, ours, published, published *fractional* uncertainty (or None)."""
    a = auger.set_index("name")
    s = t[t["fstar_chab"].notna() & t["theta_e"].notna()]
    re2 = dm["re2_table"].set_index("name")
    common_re2 = re2.index.intersection(a.index)
    rows = [
        ("Mass inside the ring (median, vs SLACS IX)", 10 ** ix["log_m_ein_median_diff"], 1.0, None),
        ("σ★ in the fibre / σ_SIE", iso["f_sdss"][0], 0.948, 0.008 / 0.948),
        ("σ★ at R_e/8 / σ_SIE", iso["f_e8"][0], 1.019, 0.008 / 1.019),
        ("intrinsic scatter of that ratio", iso["intrinsic_scatter_frac"], 0.075, None),
        ("stellar share inside the ring, Chabrier (mean)", s["fstar_ein_chab_ours"].mean(),
         s["fstar_chab"].mean(), None),
        ("stellar share inside the ring, Salpeter (mean)", s["fstar_ein_salp_ours"].mean(),
         s["fstar_salp"].mean(), None),
        ("dark share inside R_e/2, Chabrier (mean)", re2.loc[common_re2, "fdm_re2_chab"].mean(),
         a.loc[common_re2, "fdm_chab"].mean(), None),
        ("density slope γ′ (population mean)", slopes["population"]["mu"], 2.078, 0.027 / 2.078),
        ("intrinsic scatter of γ′", slopes["population"]["scatter"], 0.16, 0.02 / 0.16),
        ("stellar mass, lensing+dynamics ÷ Salpeter (α)", 10 ** dc["population"]["mu"], 1.0,
         10**0.03 - 1),
    ]
    if ours.get("n_compared"):
        rows.append(("our Einstein radii from the images (median)", 1 + ours["median_frac_diff"],
                     1.0, 0.02))
    return rows


def print_summary(res, rows):
    ix, iso, sl, dm, dc, disc, ours = (res[k] for k in ("slacs_ix", "isothermal", "slopes",
                                                         "dark", "decomposition", "discovery",
                                                         "our_einstein_radii"))
    print("\n=== Weighing galaxies with Einstein rings ===")
    print(f"\nSample: {res['sample']['n_candidates']} SLACS candidates imaged by Hubble; "
          f"{res['sample']['n_grade_a']} definite lenses, {ix['n']} with an Einstein radius.")
    print(f"  Einstein radii {ix['theta_e_range'][0]:.2f}-{ix['theta_e_range'][1]:.2f}\" = "
          f"{ix['r_ein_kpc_range'][0]:.1f}-{ix['r_ein_kpc_range'][1]:.1f} kpc; masses inside them "
          f"{ix['m_ein_range'][0]:.2e}-{ix['m_ein_range'][1]:.2e} Msun")
    print(f"  reproduce SLACS IX: R_E to {ix['r_ein_max_abs_diff_kpc']:.3f} kpc, log M_E median "
          f"{ix['log_m_ein_median_diff']:+.4f} dex (max {ix['log_m_ein_max_abs_diff']:.4f})")

    print(f"\nHow the lenses were found: second redshift recovered blind in {disc['n_recovered']} of "
          f"{disc['n']} SDSS spectra")
    for r_ in disc["confirmation_vs_theta"]:
        print(f"  predicted ring {r_['lo']:.1f}-{r_['hi']:.1f}\": {r_['n']:>3} candidates, "
              f"{100 * r_['frac_a']:.0f}% confirmed (A), {100 * r_['frac_ab']:.0f}% A or B")

    if ours.get("n_compared"):
        print(f"\nOur own lens models: {ours['n_fitted']} fitted, {ours['n_compared']} compared with "
              f"published Einstein radii")
        print(f"  median {100 * ours['median_frac_diff']:+.2f}%, robust scatter "
              f"{100 * ours['robust_scatter']:.2f}%, {ours['n_within_2pct']} within 2%, "
              f"{ours['n_within_5pct']} within 5%, {ours['n_off_by_10pct']} off by >10% "
              f"({', '.join(ours['names_off_by_10pct'])})")
        print(f"  flagged by their own diagnostics: {ours['n_flagged']} "
              f"({', '.join(f['name'] + ': ' + f['flag'] for f in ours['flagged'])})")
        if ours.get("n_unflagged"):
            print(f"  unflagged ({ours['n_unflagged']}): median {100 * ours['unflagged_median_frac_diff']:+.2f}%, "
                  f"robust scatter {100 * ours['unflagged_robust_scatter']:.2f}%, "
                  f"{ours['unflagged_n_within_5pct']} within 5%, {ours['unflagged_n_off_by_10pct']} off by >10%")
        comb = res.get("our_einstein_radii_combined")
        if comb and comb.get("hand_masked"):
            print(f"  with neighbours masked by hand in {', '.join(comb['hand_masked'])}: median "
                  f"{100 * comb['median_frac_diff']:+.2f}%, robust scatter {100 * comb['robust_scatter']:.2f}%, "
                  f"{comb['n_within_5pct']} within 5%, {comb['n_off_by_10pct']} off by >10%")
        p1 = res.get("our_einstein_radii_pass1")
        if p1:
            print(f"  first blind pass ({p1['n_compared']} compared): median {100 * p1['median_frac_diff']:+.2f}%, "
                  f"robust scatter {100 * p1['robust_scatter']:.2f}%, {p1['n_within_5pct']} within 5%, "
                  f"{p1['n_off_by_10pct']} off by >10%")
        sv = res.get("shear_vs_sie")
        if sv:
            print(f"  adding external shear ({sv['n']} ring-subset lenses): radii change by "
                  f"{100 * sv['median_change']:+.2f}% median, {100 * sv['robust_scatter']:.2f}% scatter, "
                  f"{sv['n_changed_over_5pct']} by >5%; median shear {sv['median_shear']:.3f}")
        iso2 = res.get("isothermal_with_our_radii")
        if iso2:
            print(f"  isothermal test with our own radii ({iso2['n']} lenses): sigma_e8 = "
                  f"{iso2['f_e8'][0]:.3f} +/- {iso2['f_e8'][1]:.3f} sigma_SIE, fibre {iso2['f_sdss'][0]:.3f}")

    f8, fs = iso["f_e8"], iso["f_sdss"]
    print(f"\nIsothermal test ({iso['n']} lenses): sigma_fibre = {fs[0]:.3f} +/- {fs[1]:.3f} sigma_SIE "
          f"(published 0.948 +/- 0.008); sigma_e8 = {f8[0]:.3f} +/- {f8[1]:.3f} sigma_SIE (published "
          f"1.019 +/- 0.008); intrinsic scatter {100 * iso['intrinsic_scatter_frac']:.1f}% (published 7.5%)")

    p = sl["population"]
    print(f"\nDensity slopes ({sl['n']} lenses): <gamma'> = {p['mu']:.3f} +/- {p['mu_err']:.3f}, "
          f"intrinsic scatter {p['scatter']:.3f} (+{p['scatter_hi']:.3f}/-{p['scatter_lo']:.3f}); "
          f"published 2.078 +/- 0.027, 0.16 +/- 0.02 (Auger et al. 2010)")
    for k, v in res["slope_variants"].items():
        print(f"  variant {k}: <gamma'> = {v['mu']:.3f} +/- {v['mu_err']:.3f}")
    c = res["vs_auger2010"]
    print(f"  lens by lens vs Auger et al. 2010 ({c['n']}): median diff {c['median_diff']:+.3f}, "
          f"rms {c['rms_diff']:.3f}")

    print(f"\nDark matter inside the ring ({dm['n']} lenses with stellar masses):")
    for imf, lab in (("chab", "Chabrier"), ("salp", "Salpeter")):
        d = dm[imf]
        print(f"  {lab}: stars {d['fstar_mean']:.2f} (scatter {d['fstar_std']:.2f}) of the lensing mass, "
              f"dark {d['fdm_mean']:.2f}; stars exceed the total in {d['n_fstar_above_1']} "
              f"(by 2 sigma in {d['n_fstar_above_1_by_2sigma']})")
    print(f"  heaviest stellar population allowed: median {dm['alpha_max_median']:.2f}x Chabrier "
          f"(min {dm['alpha_max_min']:.2f}); Salpeter = {dm['salp_over_chab']:.2f}x, allowed in "
          f"{100 * dm['frac_lenses_allowing_salpeter']:.0f}% of lenses")
    print(f"  dark share vs sigma: Spearman {dm['spearman_fdm_vs_sigma']:+.2f} (p = {dm['p_fdm_vs_sigma']:.2g}); "
          f"vs stellar mass {dm['spearman_fdm_vs_mstar']:+.2f} (p = {dm['p_fdm_vs_mstar']:.2g})")
    print(f"  dark share inside r_e/2: {dm['fdm_re2_chab_mean']:.2f} (Chabrier), "
          f"{dm['fdm_re2_salp_mean']:.2f} (Salpeter)")

    p = dc["population"]
    c = res["vs_treu2010"]
    print(f"\nStars + NFW halo ({dc['n']} lenses): <log alpha> = {p['mu']:+.3f} +/- {p['mu_err']:.3f} "
          f"relative to Salpeter (Chabrier would be {-dc['log_salp_over_chab']:+.2f}); intrinsic "
          f"scatter {p['scatter']:.3f}; published 0.00 +/- 0.03 (Treu et al. 2010)")
    print(f"  trend: d log alpha / d log sigma = {dc['trend_slope']:.2f} +/- {dc['trend_slope_err']:.2f}; "
          f"Treu et al.'s own table through our fit: {c['trend_slope_their_table']:.2f} +/- "
          f"{c['trend_slope_err_their_table']:.2f} (published 1.31 +/- 0.16)")
    print(f"  lens by lens vs Treu et al. 2010 ({c['n']}): median diff {c['median_diff']:+.3f}, "
          f"rms {c['rms_diff']:.3f}")
    for k, v in res["decomposition_variants"].items():
        print(f"  variant {k}: <log alpha> = {v['mu']:+.3f} +/- {v['mu_err']:.3f}")

    print("\nGraded against the literature (ours / published - 1):")
    for label, o, pub, u in rows:
        band = f"   (published uncertainty {100 * u:.1f}%)" if u else ""
        print(f"  {label:<50} {100 * (o / pub - 1):+6.1f}%{band}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refit", action="store_true", help="redo every lens model")
    parser.add_argument("--processes", type=int, default=4)
    parser.add_argument("--no-figures", action="store_true")
    parser.add_argument("--offline", action="store_true", help="never touch the network")
    args = parser.parse_args(argv)

    if not args.offline:
        from fetch_data import fetch_all, missing
        if any(missing(g) for g in ("tables", "literature", "sdss")):
            fetch_all()

    import lensmodel
    grade_a = catalog.read_grade_a()
    lensmodel.fit_all(grade_a, processes=args.processes, refit=args.refit)
    lensmodel.fit_all(grade_a, variant="hand_masked", processes=args.processes, refit=args.refit)

    t = A.lens_table()
    res = {"sample": {"n_candidates": len(catalog.read_candidates()), "n_grade_a": len(t)}}
    res["slacs_ix"] = A.reproduce_slacs_ix(t)
    res["cosmology"] = A.cosmology_dependence(t)
    res["isothermal"] = A.isothermal_test(t)
    slopes = A.density_slopes(t)
    res["slopes"] = slopes
    variants = {"beta=+0.25": A.density_slopes(t, beta=0.25),
                "beta=-0.25": A.density_slopes(t, beta=-0.25),
                "seeing=1.0": A.density_slopes(t, seeing=1.0),
                "seeing=2.0": A.density_slopes(t, seeing=2.0)}
    res["slope_variants"] = {k: v["population"] for k, v in variants.items()}
    res["dark"] = A.dark_fractions(t, slopes)
    dc = A.decompositions(t)
    res["decomposition"] = dc
    dvar = {"beta=+0.25": A.decompositions(t, beta=0.25),
            "r_break=15 kpc": A.decompositions(t, r_break_kpc=15.0),
            "r_break=60 kpc": A.decompositions(t, r_break_kpc=60.0)}
    res["decomposition_variants"] = {k: v["population"] for k, v in dvar.items()}
    res["discovery"] = A.discovery()
    res["our_einstein_radii"] = A.our_einstein_radii(t)
    res["isothermal_with_our_radii"] = (A.isothermal_test_with_ours(t, res["our_einstein_radii"])
                                        if res["our_einstein_radii"].get("n_compared") else None)
    pass1_path = MODEL_DIR / "pass1" / "summary.json"
    res["our_einstein_radii_pass1"] = json.loads(pass1_path.read_text()) if pass1_path.exists() else None
    has_shear = any((MODEL_DIR / f"{n}_shear.json").exists() for n in t["name"])
    res["our_einstein_radii_shear"] = A.our_einstein_radii(t, "shear") if has_shear else None
    if has_shear:
        # How much adding external shear moves our own radii (the ring-subset lenses it was run on)
        sh = res["our_einstein_radii_shear"]["table"].set_index("name")
        base = res["our_einstein_radii"]["table"].set_index("name")
        common = sh.index[(sh["failed"] == "")].intersection(base.index[base["failed"] == ""])
        ch = sh.loc[common, "b_ours"] / base.loc[common, "b_ours"] - 1
        res["shear_vs_sie"] = {
            "n": int(len(common)), "median_change": float(ch.median()),
            "robust_scatter": float(1.4826 * np.median(np.abs(ch - ch.median()))),
            "n_changed_over_5pct": int((ch.abs() > 0.05).sum()),
            "median_shear": float(sh.loc[common, "shear"].median())}
    res["our_einstein_radii_combined"] = A.our_einstein_radii(t, "combined")

    # Grading against the per-lens literature tables.
    auger, treu = literature.auger2010(), literature.treu2010()
    m = slopes["table"].merge(auger[["name", "gamma"]], on="name", suffixes=("", "_pub"))
    d = m["gamma"] - m["gamma_pub"]
    res["vs_auger2010"] = {"n": len(m), "median_diff": float(d.median()),
                           "rms_diff": float(np.sqrt(np.mean(d**2))),
                           "mean_of_their_table": float(auger["gamma"].mean())}
    m = dc["table"].merge(treu[["name", "log_alpha"]], on="name", suffixes=("", "_pub"))
    d = m["log_alpha"] - m["log_alpha_pub"]
    # Their table through our trend fit
    e = 0.5 * (treu["e_log_alpha_lo"] + treu["e_log_alpha_hi"])
    x = np.log10(treu["sigma"]) - np.log10(250)
    w = 1 / e**2
    X = np.vstack([x, np.ones_like(x)]).T
    cov = np.linalg.inv(X.T @ (X * w.to_numpy()[:, None]))
    coef = cov @ (X.T @ (w * treu["log_alpha"]).to_numpy())
    res["vs_treu2010"] = {"n": len(m), "median_diff": float(d.median()),
                          "rms_diff": float(np.sqrt(np.mean(d**2))),
                          "trend_slope_their_table": float(coef[0]),
                          "trend_slope_err_their_table": float(np.sqrt(cov[0, 0]))}

    rows = comparison(t, res["slacs_ix"], res["isothermal"], slopes, res["dark"], dc,
                      res["our_einstein_radii"], auger, treu)
    res["comparison"] = [{"label": r[0], "ours": float(r[1]), "published": float(r[2]),
                          "published_frac_err": r[3]} for r in rows]

    if not args.no_figures:
        import figures as F
        import spectra
        import interactive

        written = []
        written += F.figure_comparison(rows)
        names = [n for n in SHOWCASE if (MODEL_DIR / f"{n}.npz").exists()]
        written += F.figure_gallery(t, names, res["our_einstein_radii"])
        if (MODEL_DIR / f"{ANATOMY_LENS}.npz").exists():
            written += F.figure_anatomy(t, ANATOMY_LENS)
        cand = catalog.read_candidates().set_index("name")
        sp = spectra.read_spectrum(SPECTRUM_LENS)
        search = spectra.search_background(sp, cand.loc[SPECTRUM_LENS, "z_lens"])
        written += F.figure_two_redshifts(SPECTRUM_LENS, cand.loc[SPECTRUM_LENS, "z_lens"],
                                          cand.loc[SPECTRUM_LENS, "z_source"], search)
        written += F.figure_discovery(res["discovery"])
        if res["our_einstein_radii"].get("n_compared"):
            written += F.figure_einstein_radii(res["our_einstein_radii"],
                                               res["our_einstein_radii_pass1"],
                                               res["our_einstein_radii_shear"])
        written += F.figure_isothermal(t, res["isothermal"])
        written += F.figure_density_slopes(
            t, slopes, auger[["name", "gamma"]],
            [("radial orbits, β = +0.25", variants["beta=+0.25"]),
             ("tangential orbits, β = −0.25", variants["beta=-0.25"]),
             ("seeing 1.0″", variants["seeing=1.0"]), ("seeing 2.0″", variants["seeing=2.0"])])
        written += F.figure_dark_matter(t, res["dark"])
        row = next(r for r in A.slope_sample(t).itertuples() if r.name == CONSPIRACY_LENS)
        fs = dc["table"].set_index("name").loc[CONSPIRACY_LENS, "fstar"]
        written += F.figure_conspiracy(t, dc, row, A.conspiracy_profiles(row, fs),
                                       treu[["name", "log_alpha"]], slopes)
        written += interactive.build(t, res["our_einstein_radii"])
        print(f"\nWrote {len(written)} files to {FIGURE_DIR.name}/")

    RESULTS_PATH.write_text(json.dumps(_jsonable(res), indent=1, default=str))
    print_summary(res, rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
