"""Run the whole project end to end.

    python run.py                  # download what is missing, analyse, draw, summarise
    python run.py --figures-only   # redraw figures from the last analysis run

The first run downloads ~27 MB and measures ~450 SDSS-RM lags with our own
cross-correlation code (~30 minutes on 4 cores; cached afterwards). Once the
cache is warm everything runs offline in about 15 minutes (the LINMIX
regressions dominate).
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
import warnings

warnings.filterwarnings("ignore")

import analysis
import config as C
import explorer
import figures
from fetch_data import fetch_all

CACHE = C.CACHE_DIR / "analysis.pkl"


def _v(t, nd=2):
    if t is None:
        return "-"
    if isinstance(t, (int, float)):
        return f"{t:.{nd}f}"
    if len(t) == 3:
        return f"{t[0]:.{nd}f} (-{t[1]:.{nd}f}/+{t[2]:.{nd}f})"
    return f"{t[0]:.{nd}f} ± {t[1]:.{nd}f}"


def print_summary(s: dict) -> None:
    print("\n" + "=" * 78)
    print("RESULTS: reverberation mapping")
    print("=" * 78)
    print("\nNGC 5548 in 2014 (AGN STORM): rest-frame lags (days), ours vs Pei et al. 2017")
    for k, v in s["storm_lags"].items():
        print(f"  {k:<32} {_v(v['ours'])}   r_max {v['r_max']:.2f}   published {_v(v['published'])}")
    m = s["storm_mass"]
    print(f"  BLR radius {m['tau']:.2f} light-days = {m['R_au']:.0f} AU; virial product {m['vp'] / 1e7:.2f} ± "
          f"{m['vp_err'] / 1e7:.2f} x 1e7 Msun (Pei: 1.49 ± 0.49); M = {m['mbh_f447'] / 1e7:.2f} ± "
          f"{m['mbh_f447_err'] / 1e7:.2f} x 1e7 Msun with f = 4.47 (Pei: 6.66 ± 2.17)")
    print("\nThe BLR holiday: lag (obs. days), responsivity, flux lost in the anomaly; ours vs published")
    for k, v in s["holiday"].items():
        print(f"  {k:<5} lag {_v(v['tau_obs'])} ({_v(v['published_tau'])})   eta {_v(v['eta'], 3)} "
              f"({_v(v['published_eta'], 2)})   lost {100 * v['flost']:.1f}% ({100 * v['published_flost']:.0f}%)")
    print("\nNGC 5548, 1988-2008: seasonal H-beta lags (ours vs published)")
    for w in s["watch"]:
        print(f"  {w['year']:7.1f}  {_v(w['tau'], 1):>22}   published {_v(w['published'], 1):>22}   log L {w['logL']:.2f} "
              f"({w['published_logL']:.2f})")
    a = s["watch_agreement"]
    print(f"  {a['n_within_1sig']} of {a['n']} within 1 sigma; median |difference| {a['median_abs_diff_days']:.2f} d")
    b = s["breathing"]
    print(f"  NGC 5548's own R-L slope: {_v(b['bentz07_sample']['b'])} (Bentz+07: 0.73 ± 0.14); without Year 12 "
          f"{_v(b['bentz07_no_y12']['b'])} (0.66 ± 0.13); with the published 0.1-dex L errors "
          f"{_v(b['bentz07_sample_pub_errors']['b'])}; published lags and L with their errors "
          f"{_v(b['published_values']['b'])}; host flux ±10%: {b['host_systematic']['slope_more_host']:.2f}"
          f"-{b['host_systematic']['slope_less_host']:.2f}")
    sv = s["storm_vs_own_rl"]
    print(f"  2014 lag {sv['measured_tau']:.2f} d vs {sv['predicted_tau']:.1f} d expected from its own trend: "
          f"{sv['ratio']:.1f}x shorter")
    print("\nR-L relation (log tau = a + b log(L/1e44))")
    bv = s["bentz13_validation"]["linmix"]
    print(f"  Bentz+13 table through our LINMIX: a {_v(bv['a'], 3)}, b {_v(bv['b'], 3)}, scatter {bv['sig'][0]:.3f}"
          f" (published 1.527, 0.533, 0.19)")
    lr = s["local_rl"]
    print(f"  Database, {lr['n']} lags of {lr['n_obj']} AGN: a {_v(lr['all']['a'], 3)}, b {_v(lr['all']['b'], 3)}, "
          f"scatter {lr['all']['sig'][0]:.3f}")
    for k, v in s["sdssrm_rl"].items():
        p = v["published"]
        print(f"  SDSS-RM {k:<4} N={v['n']:3d}: a {_v(v['fit']['a'], 3)}, b {_v(v['fit']['b'], 3)}, scatter "
              f"{v['fit']['sig'][0]:.3f}   (Shen+24: {p['a'][0]}, {p['b'][0]}, {p['sig'][0]})")
    nh = s["sdssrm_rl_hb_no_host"]
    print(f"  SDSS-RM H-beta without host correction: a {_v(nh['a'], 3)}, b {_v(nh['b'], 3)}")
    print(f"  SDSS-RM H-beta lags vs local AGN at fixed slope 0.533: {s['sdssrm_offset_dex']:+.3f} dex")
    print("\nSDSS-RM with our own ICCF")
    for k, v in s["sdssrm_compare"].items():
        print(f"  {k:<4} measured {v['n_measured']:3d}; published detections {v['n_pub_det']:3d}, ours "
              f"{v['n_ours_det']:3d}, both {v['n_both_det']:3d}; median ratio to their ICCF lag "
              f"{v['median_ratio_iccf']:.3f}, {100 * v['frac_within_1sig_iccf']:.0f}% within 1 sigma")
    f = s["sdssrm_fpr"]
    print(f"  asymmetry (H-beta, r > 0.4): {f['pos_hi']} positive vs {f['neg_hi']} negative -> false-positive "
          f"rate ~{100 * f['fpr_hi']:.0f}%")
    for k, v in s["sdssrm_ours_rl"].items():
        if v:
            print(f"  R-L from our own {k} lags (N={v['n']}): b {_v(v['b'], 2)}, scatter {v['sig'][0]:.2f}")
    print("\nVirial factor and masses")
    fd = s["f_dynamical"]
    print(f"  <log f> from {fd['n']} dynamical models: {_v(fd['mean'], 3)}, scatter {_v(fd['scatter'], 3)} "
          f"(Shen+24: 0.62 ± 0.07, 0.31 ± 0.07)")
    fm = s["f_msigma"]
    print(f"  M-sigma joint fit ({fm['n_q']} galaxies + {fm['n_agn']} AGN): slope {_v(fm['beta'])}, log f "
          f"{_v(fm['logf'])}, scatter {fm['eps'][0]:.2f} (Woo+15 with KH13: 4.97, 0.65, 0.43)")
    v = s["virial"]
    print(f"  NGC 5548 virial test: slope {_v(v['fit_all']['b'])} over {v['n']} line/season pairs (gravity: -0.5); "
          f"log VP scatter {v['logvp_scatter']:.2f} dex")
    n = s["ngc5548_msigma"]
    print(f"  NGC 5548: log M = {n['ours_dyn']['logM']:.2f} (dynamical f), {n['ours_msigma']['logM']:.2f} (M-sigma f); "
          f"M-sigma predicts {n['msigma_prediction']:.2f} ({n['msigma_offset_dex']:+.2f} dex)")
    print()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figures-only", action="store_true", help="redraw figures from the cached analysis")
    ap.add_argument("--workers", type=int, default=4, help="processes for the SDSS-RM lags (first run only)")
    args = ap.parse_args()
    if args.figures_only:
        if not CACHE.exists():
            print("No cached analysis yet: run `python run.py` first.")
            return 1
        res = pickle.loads(CACHE.read_bytes())
    else:
        fetch_all()
        res = analysis.run(sdssrm_workers=args.workers)
        CACHE.write_bytes(pickle.dumps(res))
        C.RESULTS_PATH.write_text(json.dumps(res["summary"], indent=1, default=str))
    paths = figures.make_all(res)
    paths.append(explorer.write(res))
    print_summary(res["summary"])
    print("Figures:")
    for p in paths:
        if str(p).endswith((".png", ".html")):
            print("  ", p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
