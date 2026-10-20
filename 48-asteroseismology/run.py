"""Run the whole project end to end.

    python run.py                  # download what is missing, analyse, draw, summarise
    python run.py --figures-only   # redraw figures and the explorer from the last analysis

The first run downloads ~1.2 GB (104 stars x 18 Kepler quarters, one year of
16 Cyg A and B at 1-minute cadence, 25 years of SOHO/VIRGO) and fits every
star (~10 minutes on 6 cores for the red giants, ~15 for 27 years of the Sun;
longer on a busy machine). Everything is cached in data/; afterwards the
analysis and figures run offline in a few minutes.
"""

from __future__ import annotations

import argparse
import json
import warnings

warnings.filterwarnings("ignore")

import config as C


def _pc(x, nd=2):
    return f"{100 * x:+.{nd}f}%"


def print_summary(res: dict) -> None:
    g = res["global"]
    print("\n" + "=" * 78)
    print("RESULTS: asteroseismology of 104 Kepler red giants")
    print("=" * 78)
    print(f"\nGlobal oscillation parameters, this work vs the catalogues (median, robust scatter):")
    for q, name in (("numax", "numax"), ("dnu", "Dnu  ")):
        for cat, cname in (("ak", "APOKASC-3"), ("yu", "Yu+2018 ")):
            s = g[f"{q}_vs_{cat}"]
            print(f"  {name} vs {cname}: {_pc(s['median'])} +- {100 * s['median_err']:.2f}%, scatter "
                  f"{100 * s['scatter']:.2f}% (N={s['n']})")
    s = g["dnu_acf_vs_ak"]
    print(f"  (autocorrelation Dnu alone vs APOKASC-3: {_pc(s['median'])}, scatter {100 * s['scatter']:.2f}%; "
          f"red clump {_pc(g['dnu_acf_vs_ak_RC']['median'])})")
    print(f"  amplitude per radial mode vs Yu: {_pc(g['amp_vs_yu']['median'], 1)}; granulation "
          f"{_pc(g['gran_vs_yu']['median'], 1)}; envelope width {_pc(g['width_vs_yu']['median'], 1)}")
    d = res["dnu_numax"]
    for k in ("this work", "Yu et al. 2018 catalogue (our fit)"):
        f = d[k]
        print(f"  Dnu = {f['alpha']:.3f} numax^{f['beta']:.3f} (+-{f['beta_err']:.3f}), {f['n']} stars: {k}")
    sun = res["sun"]
    print(f"\nThe Sun (VIRGO, {len(sun['years_used'])} years): numax = {sun['numax_sun']:.1f} +- "
          f"{sun['numax_sun_err']:.1f} uHz (year-to-year scatter {sun['numax_sun_scatter']:.1f}); "
          f"Dnu = {sun['dnu_sun']:.3f} +- {sun['dnu_sun_err']:.3f} uHz")
    mr = res["masses_radii"]
    print("\nMasses and radii:")
    for key, label in (("ak_scale_vs_apokasc", "APOKASC-3 scale vs APOKASC-3"),
                       ("yu_scale_vs_yu", "Yu scale vs Yu (uncorrected)"),
                       ("self_vs_apokasc", "self-calibrated vs APOKASC-3")):
        print(f"  {label}: mass {_pc(mr[key]['mass']['median'], 1)} (scatter "
              f"{100 * mr[key]['mass']['scatter']:.1f}%), radius {_pc(mr[key]['radius']['median'], 1)} "
              f"(scatter {100 * mr[key]['radius']['scatter']:.1f}%)")
    gr = res["gaia_radii"]
    print("\nSeismic / Gaia radius (median):")
    for s, label in (("self", "self-calibrated"), ("yu", "Yu solar values"), ("ak", "APOKASC-3 scale")):
        print(f"  {label:16s} RGB {gr[s]['RGB']['median']:.3f} +- {gr[s]['RGB']['median_err']:.3f}   "
              f"RC {gr[s]['RC']['median']:.3f} +- {gr[s]['RC']['median_err']:.3f}")
    st = res["states"]
    print(f"\nPeriod spacings: attempted {st['n_attempted']}, detected {st['n_detected']}")
    dv = st["dpi1_vs_vrard"]
    print(f"  vs Vrard+2016: {dv['n']} stars, median {_pc(dv['median'])}, scatter {100 * dv['scatter']:.2f}%, "
          f"{dv['n_within_2pc']} within 2%, outliers {dv['outliers']}")
    for k in ("APOKASC-3", "Yu et al. 2018", "Vrard et al. 2016"):
        print(f"  shell-H vs core-He classification agrees with {k}: {st[k]['n_agree']}/{st[k]['n_compared']}")
    print("\n16 Cyg A and B (self-calibrated scaling relations vs interferometry / models):")
    for name, s in res["sixteen_cyg"].items():
        sc = s["scaling"]["self"]
        lit = s["literature"]
        print(f"  {name}: numax {s['numax']:.0f} (Lund {s['numax_lund']:.0f}), Dnu {s['dnu']:.2f} (Lund "
              f"{s['dnu_lund']:.1f}); R = {sc['radius']:.3f} +- {sc['radius_err']:.3f} (CHARA "
              f"{lit['radius_interf']} +- {lit['radius_interf_err']}); M = {sc['mass']:.2f} +- "
              f"{sc['mass_err']:.2f} (models {lit['mass_model']} +- {lit['mass_model_err']})")
    print()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--figures-only", action="store_true")
    args = ap.parse_args()

    if not args.figures_only:
        import analysis
        import fetch_data
        import mainseq
        import sample
        import stars
        import sun

        print("1. Data")
        fetch_data.fetch_all()
        print("2. Red giants (per-star pipeline, cached)")
        stars.analyse_all(sample.load_or_build_sample(quiet=True)["kic"])
        print("3. The Sun (SOHO/VIRGO, one fit per year)")
        sun.analyse_sun()
        print("4. 16 Cyg A and B")
        mainseq.analyse_16cyg()
        print("5. Ensemble analysis")
        res = analysis.run_analysis()
    else:
        res = json.loads(C.RESULTS_PATH.read_text())

    import explorer
    import figures
    import sonify

    print("6. Sounds")
    picks = figures._ladder_stars(__import__("pandas").read_csv(C.DATA_DIR / "results_stars.csv",
                                                                index_col="kic"))
    extra = [k for k in (figures.HERO, figures.RGB_STAR, figures.RC_STAR) if k not in picks]
    kics = [int(k) for k in picks] + extra
    have = json.loads((C.AUDIO_DIR / "audio.json").read_text()) if (C.AUDIO_DIR / "audio.json").exists() else {}
    if not args.figures_only or any(str(k) not in have for k in kics):
        sonify.make_audio(kics)
    print("7. Figures")
    figures.make_all()
    print("8. Explorer")
    explorer.build(kics, res)
    print_summary(res)


if __name__ == "__main__":
    main()
