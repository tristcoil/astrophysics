"""Run L4 end to end.

    python fetch_data.py      # once: ~1.1 GB of inputs
    python run.py             # every measurement, figure and the explorer
    python run.py --refit     # redo the slow fits instead of using data/cache/
    python run.py --no-figures

Writes data/results.json and figures/, and prints the headline numbers next
to the published ones. Needs no network once the inputs are cached. The slow
steps (disc mixture fits, ~5 min each; nested sampling of the kilonova,
~10-20 min) are cached after the first run.
"""

from __future__ import annotations

import argparse
import os
import time
import warnings

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")   # be polite on a shared machine
warnings.filterwarnings("ignore")

import analysis  # noqa: E402
from config import FIGURES, REFERENCE  # noqa: E402


def _q(q, fmt="{:.3f}"):
    lo, mid, hi = q
    return f"{fmt.format(mid)} (+{fmt.format(hi - mid)} -{fmt.format(mid - lo)})"


def summary(r) -> None:
    print("\n" + "=" * 96)
    print("RESULTS -- what we found, beside the published values")
    print("=" * 96)
    s = r["sample"]
    print(f"\nGALAH DR4: {s['n_catalogue']:,} stars; {s['n_clean']:,} clean field stars "
          f"({s['n_giants']:,} giants, {s['n_dwarfs']:,} dwarfs; {s['n_disc']:,} on disc orbits)")
    k = r["knee"]
    print("\n1. The alpha-knee")
    for c in ("giants", "dwarfs"):
        f = k[f"disc_{c}"]
        print(f"  disc {c:6s}: knee [Fe/H] = {f['knee']:+.3f} ({f['knee_lo']:+.3f}..{f['knee_hi']:+.3f}), "
              f"plateau [Mg/Fe] = {f['plateau']:.3f}, ridge method {f['ridge_method']['knee']:+.3f}")
    if isinstance(k.get("apogee_disc_knee"), dict):
        a = k["apogee_disc_knee"]
        print(f"  disc APOGEE giants ({k['apogee_disc_giants']:,}): knee {a['knee']:+.3f} "
              f"({a['knee_lo']:+.3f}..{a['knee_hi']:+.3f}), plateau {a['plateau']:.3f}")
    for c in ("giants", "dwarfs"):
        if f"apogee_vs_galah_{c}" in k:
            a = k[f"apogee_vs_galah_{c}"]
            print(f"  GALAH - APOGEE, {a['n']:,} common {c}: [Fe/H] {a['d_feh_median']:+.3f} "
                  f"(sigma {a['d_feh_sigma']:.3f}), [Mg/Fe] {a['d_mgfe_median']:+.3f} (sigma {a['d_mgfe_sigma']:.3f})")
    for key, f in k.items():
        if key.startswith("gse[") and key.endswith("[giants]"):
            print(f"  {key[4:-9]:30s}: knee {f['knee']:+.2f} ({f['knee_lo']:+.2f}..{f['knee_hi']:+.2f}), "
                  f"plateau {f['plateau']:.3f}, n = {f['n']}, median [Fe/H] {f['median_feh']:+.2f}, "
                  f"knee vs straight line: dchi2 = {f['delta_chi2_vs_line']:.1f}")
    print(f"  published GSE knee: ~{REFERENCE['gse_knee']['mackereth2019']} (Mackereth+19), "
          f"{REFERENCE['gse_knee']['literature_range']} (others)")
    print(f"  H3 boxes reproduced: Helmi {k['h3_counts']['helmi']:,} (H3: 77,205), "
          f"Koppelman {k['h3_counts']['koppelman']:,} (H3: 12,343)")

    p = r["process"]
    print("\n2. Prompt shares at solar metallicity (f_cc)")
    print("  element   giants(ours)  dwarfs(G22 recipe)  dwarfs(own plateau)   Griffith+22 (DR3)")
    gr = REFERENCE["griffith2022_fcc"]
    for el in ["fe"] + [e.lower() for e in ("C", "N", "O", "Na", "Al", "Si", "K", "Ca", "Sc", "Ti", "V", "Cr",
                                            "Mn", "Co", "Ni", "Cu", "Zn", "Rb", "Sr", "Y", "Zr", "Mo", "Ru",
                                            "Ba", "La", "Ce", "Nd", "Sm")] + ["eu", "eu_dr3", "eu_dr4"]:
        cells = []
        for key in ("giants_ours", "dwarfs_griffith_recipe", "dwarfs_ours"):
            v = p[key].get(el, {})
            cells.append(f"{v['f_cc']:.2f}±{v['f_cc_err'] or 0:.2f}" if v.get("f_cc") is not None else "   --    ")
        ref = gr.get(el.capitalize())
        print(f"  {el:8s}  {cells[0]:>11s}   {cells[1]:>11s}        {cells[2]:>11s}        "
              f"{'' if ref is None else f'{ref:.2f}'}")
    print(f"  plateaus: giants {p['giants_ours']['_plateau_zp']:.3f}, dwarfs {p['dwarfs_ours']['_plateau_zp']:.3f} "
          f"(Griffith+22 assumed 0.30)")

    n = r["ncapture"]
    print("\n3. Neutron-capture elements")
    for q, off in n["gse_minus_disc"].items():
        o2 = n["gse_minus_disc_cnn_regime"][q]
        print(f"  GSE - disc [{q}]: {off['offset']:+.3f} +- {off['err']:.3f} (n={off['n']}); "
              f"[Fe/H] >= -1 only: {o2['offset']:+.3f} +- {o2['err']:.3f} (n={o2['n']})")
    gm = n["gse_median"]
    b22 = n["buder2022"]
    print(f"  GSE medians: [Fe/H] {gm['fe_h']:+.2f} (Buder+22: {b22['fe_h'][0]}), [Mg/Fe] {gm['mg_fe']:+.2f} "
          f"({b22['mg_fe'][0]}), [Eu/Fe] {gm['eu_fe']:+.2f} ({b22['eu_fe'][0]}), [Ba/Fe] {gm['ba_fe']:+.2f} "
          f"({b22['ba_fe'][0]})")
    print(f"  pure-r [Ba/Eu] = {n['ba_eu_pure_r']:+.2f}, pure-s = {n['ba_eu_pure_s']:+.2f}")
    for q in ("mg_fe", "ba_fe", "eu_fe"):
        j = n["scatter_jinabase"][q]
        print(f"  JINAbase spread of [{q}]: " + ", ".join(f"{x['x']:+.1f}: {x['robust_sigma']:.2f}" for x in j))
    print(f"  DR4's own Eu vs DR3: r = {n['eu_dr4_vs_dr3']['r']:.2f}, spread {n['eu_dr4_vs_dr3']['sigma']:.2f}; "
          f"DR4's own vs CNN (stars DR3 lacks): r = {n['eu_dr4_vs_cnn']['r']:.2f}, spread {n['eu_dr4_vs_cnn']['sigma']:.2f}")

    kn = r["kilonova"]
    v = REFERENCE["villar2017"]
    print(f"\n4. AT2017gfo ({kn['n_points']} points, {kn['n_upper']} upper limits, {kn['n_filters']} filters; "
          f"log Z = {kn['logz']:.1f})")
    for name, m in kn["modes"].items():
        q = m["quantiles"]
        print(f"  [{name}, weight {m['weight']:.2f}] blue {_q(q['m_blue'])} at {q['v_blue'][1] / 2.998e5:.3f}c; "
              f"purple {_q(q['m_purple'])} at {q['v_purple'][1] / 2.998e5:.3f}c; red {_q(q['m_red'])} at "
              f"{q['v_red'][1] / 2.998e5:.3f}c; total {_q(q['m_total'])}; purple+red {_q(q['m_purple_red'])}; "
              f"sigma {q['sigma'][1]:.3f}")
    print(f"  Villar+17: blue {v['blue']['m']} at {v['blue']['v']}c, purple {v['purple']['m']} at "
          f"{v['purple']['v']}c, red {v['red']['m']} at {v['red']['v']}c, sigma {v['sigma']}")
    print(f"  rms: best {kn['best_rms']:.3f} mag; Villar's parameters in our model {kn['villar_rms_our_model']:.3f} "
          f"(log L {kn['best_loglike']:.1f} vs {kn['villar_loglike_our_model']:.1f})")

    b = r["budget"]
    print("\n5. Europium budget")
    print(f"  Eu share of solar r-process matter (Kr-U): {100 * b['shares']['eu_share_kr_up']:.3f}%; "
          f"Eu's r-process share {b['shares']['eu_r_fraction']:.2f}")
    print(f"  Galaxy: {b['galaxy']['m_eu']:.1f} +- {b['galaxy']['m_eu_err']:.1f} Msun of r-process Eu "
          f"(disc giants' mean Eu level {b['mean_eu_h_linear_disc_giants']:.2f} x solar)")
    for reading, (lo, mid, hi) in b["per_event_eu"].items():
        bb = b["budgets"][reading]
        print(f"  {reading:11s}: {mid:.2e} Msun Eu per merger -> {bb['n_needed']:.2e} mergers needed, "
              f"{bb['rate_needed_per_myr']:.1f}/Myr; supply/demand at GWTC-4 low/mid/high rate: "
              + "/".join(f"{x:.2f}" for x in bb["supply_over_demand"]))
    print(f"  history factor (mean past rate / today): {b['history_factor']:.2f}")
    vc = b["villar_check"]
    print(f"  Villar+17's per-merger needs with the 2017 rate: light {vc['needed_light_2017']}, heavy "
          f"{vc['needed_heavy_2017']} Msun; with GWTC-4's: {vc['needed_light_now']:.3f}, {vc['needed_heavy_now']:.4f}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refit", action="store_true")
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args(argv)
    t0 = time.time()
    results, art = analysis.run(refit=args.refit)
    analysis.save_results(results)
    if not args.no_figures:
        import figures
        import interactive

        written = figures.make_all(results, art)
        written.append(interactive.write_explorer(results, art))
        print(f"\nwrote {len(written)} files to {FIGURES}")
    summary(results)
    print(f"\n(total {time.time() - t0:.0f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
