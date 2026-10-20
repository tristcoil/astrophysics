"""Run the whole project end to end.

    python run.py                 # fetch (if needed), every analysis, every figure
    python run.py --figures-only  # redraw figures from the last analysis run
    python run.py --refresh       # recompute the cached time-dilation fits too

The slow parts (the White et al. replication, the bootstrap and injection tests
of the global supernova fit) are cached in data/derived/; with a warm cache the
run takes about a minute and needs no network. From scratch it takes 30-60
minutes on 6 cores (set D7_WORKERS to change the number of processes).
"""

from __future__ import annotations

import argparse
import pickle
import sys
import warnings

warnings.filterwarnings("ignore")

import analysis
import figures
from config import DERIVED_DIR, FIGURE_DIR
from fetch_data import fetch_all

CACHE = DERIVED_DIR / "analysis_cache.pkl"


def _pm(v, e, spec=".3f"):
    return f"{v:{spec}} ± {e:{spec}}"


def print_summary(s):
    td, tc, to = s["time_dilation"], s["cmb_temperature"], s["tolman"]
    print("\n" + "=" * 78)
    print("IS THE UNIVERSE REALLY EXPANDING?  expansion vs tired light")
    print("=" * 78)

    g, w = td["global_fit"], td["white_replication"]
    pub = td["published"]["white24"]
    print(f"\n1. Supernova time dilation, D = (1 + z)^b   (expansion b = 1, tired light b = 0)")
    print(f"   {td['sample']['n_ia']} likely Type Ia of {td['sample']['n_des']} DES-SN5YR supernovae")
    print(f"   our global template fit:   b = {_pm(g['b'], g['b_err_bootstrap'])}  (bootstrap; "
          f"curvature {g['b_err_curvature']:.3f})")
    print(f"     tired light b = 0:       Δχ² = {g['delta_chi2_tired_light']:,.0f} worse "
          f"({g['sigma_from_b0']:.0f}σ from b)")
    print(f"   White et al. replication:  b = {_pm(w['b'], w['b_err_bootstrap'])}  "
          f"(published {pub['b'][0]} ± {pub['b'][1]} stat ± {pub['b_sys']} sys)")
    print(f"     same widths, Nelder-Mead from w = 1: b = {w['b_nelder_mead_from_1']:.3f}; "
          f"global minimum: b = {w['b_global_grid_minimum']:.3f}")
    print(f"     their own widths refitted: b = {_pm(w['white_widths_refit']['b'], w['white_widths_refit']['bootstrap_err'])}")
    print(f"     scatter minimum (method 1, riz): b = {w['dispersion_b_min_riz']:.2f} "
          f"(their curves: {w['dispersion_b_min_riz_white']:.2f})")
    print(f"     null test (reference not de-dilated): b = {_pm(w['null_test_b'], w['null_test_b_err'], '.2f')}")
    if "injections" in g:
        for b, v in g["injections"].items():
            print(f"   injection test, true b = {b:.1f}: recovered {v['mean_b']:.3f} ± {v['std_b']:.3f} "
                  f"({v['n']} simulations)")
    print("   robustness of the global fit:")
    for label, v in g["variants"].items():
        print(f"     {label:<28} b = {v['b']:.3f} ± {v['b_err_curvature']:.3f}  ({v['n_sne']} SNe)")

    b = tc["baseline"]
    m = tc["baseline_models"]
    print(f"\n2. CMB temperature, T = T0 (1 + z)^(1 - β)   (expansion β = 0, tired light β = 1)")
    print(f"   {tc['n_baseline']} measurements:  β = {b['beta']:+.4f} +{b['plus']:.4f} -{b['minus']:.4f}"
          f"   χ² = {b['chi2']:.1f} for {b['dof']} dof")
    print(f"   χ²: expansion {m['chi2_expansion']:.1f}, tired light {m['chi2_tired']:,.0f}")
    t0 = tc["t0_highz"]
    print(f"   T0 from distant measurements alone: {t0['T0']:.4f} ± {max(t0['plus'], t0['minus']):.4f} K "
          f"(FIRAS 2.72548)")
    for k, v in tc["noterdaeme11"].items():
        p = tc["published"]["noterdaeme11_table"][k]
        print(f"   Noterdaeme+11 data, {k:<6} β = {v['beta']:+.3f} ± {max(v['plus'], v['minus']):.3f}"
              f"   (published {p[0]:+.3f} ± {p[1]:.3f})")
    r = tc["riechers22_all_rows"]
    print(f"   Riechers+22 table: β = {r['beta'] * 1e3:.1f} ± {max(r['plus'], r['minus']) * 1e3:.1f} ×10⁻³ "
          f"(published 3.4 +8.1 -7.3 ×10⁻³)")

    print(f"\n3. Tolman surface-brightness test, dimming 2.5 log (1 + z)^n   "
          f"(expansion n = 4, tired light n = 1, both minus evolution)")
    for model, v in to["results"].items():
        print(f"   {v['label']:<48} n_R = {_pm(v['n_R'], v['n_R_err'], '.2f')}   "
              f"n_I = {_pm(v['n_I'], v['n_I_err'], '.2f')}")
    p = to["published"]
    print(f"   Lubin & Sandage 2001 (q0 = 1/2): n_R = {p['n_R'][0]} ± {p['n_R'][1]}, "
          f"n_I = {p['n_I'][0]} ± {p['n_I'][1]}; tired light {p['n_R_tired'][0]}, {p['n_I_tired'][0]}")
    print(f"\n   Runtime {s['runtime_s']:.0f} s.")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--figures-only", action="store_true")
    parser.add_argument("--no-figures", action="store_true")
    parser.add_argument("--refresh", action="store_true")
    args = parser.parse_args(argv)
    if args.figures_only:
        if not CACHE.exists():
            print("no cached analysis yet; run `python run.py` first")
            return 1
        summary, art = pickle.loads(CACHE.read_bytes())
    else:
        fetch_all()
        summary, art = analysis.run(refresh=args.refresh)
        analysis.save_results(summary)
        CACHE.write_bytes(pickle.dumps((summary, art)))
    if not args.no_figures:
        print("Figures")
        written = figures.make_all(*art)
        print(f"  wrote {len(written)} files to {FIGURE_DIR}")
    print_summary(summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
