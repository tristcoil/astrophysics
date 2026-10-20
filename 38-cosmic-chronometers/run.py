"""Run every fit, write data/results.json, and draw every figure.

    python run.py                 # analysis + figures (~3 minutes)
    python run.py --figures-only  # redraw from the cached analysis
    python run.py --no-figures    # numbers only

Needs the cache from ``python fetch_data.py``; after that it runs offline.
"""

from __future__ import annotations

import argparse
import json
import pickle
import time
import warnings

import numpy as np

import config as cfg

CACHE = cfg.DATA_DIR / "analysis_cache.pkl"


def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _jsonable(obj.tolist())
    if isinstance(obj, (np.floating, float)):
        return None if not np.isfinite(obj) else float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def _fmt(s, nd=2):
    return f"{s['median']:.{nd}f} +{s['plus']:.{nd}f} -{s['minus']:.{nd}f}"


def summary(r: dict) -> None:
    f = r["flat_lcdm"]["full"]
    print("\n=== Cosmic chronometers: 32 H(z) measurements, 0.07 < z < 1.965 ===")
    print(f"flat LCDM, full covariance:  H0 = {_fmt(f['H0'], 1)}   Om = {_fmt(f['Om'], 3)}"
          f"   chi2 = {f['best']['chi2']:.1f} / {f['best']['dof']}")
    print(f"   transition z_t = {_fmt(f['z_transition'])}   q0 = {_fmt(f['q0'])}"
          f"   age = {_fmt(f['age_gyr'], 1)} Gyr")
    for name in r["flat_lcdm"]:
        if name == "full":
            continue
        g = r["flat_lcdm"][name]
        print(f"   [{name}] H0 = {_fmt(g['H0'], 1)}  Om = {_fmt(g['Om'], 3)}")
    o = r["open_lcdm"]
    print(f"open LCDM: Om = {_fmt(o['Om'], 2)}, OL = {_fmt(o['OL'], 2)}, "
          f"P(accelerating today) = {o['p_accelerating_today']:.3f}")
    print("Gaussian process:")
    for name, g in r["gp"].items():
        print(f"   {name:45s} H0 = {_fmt(g['H0'], 1)}  q0 = {_fmt(g['q0'])}"
              f"  P(q0<0) = {g['p_accelerating_today']:.3f}  z_t = {_fmt(g['z_transition'])}"
              f" ({100 * g['frac_transition_in_data']:.0f}% of draws turn inside the data)")
    print("Two straight lines in adot(z):")
    for name, p in r["piecewise"].items():
        pr = p["profile"]
        print(f"   anchor {name:32s} z_t = {pr['best']:.2f} (68%: {pr['68'][0]:.2f}-{pr['68'][1]:.2f})"
              f"  line worse by dchi2 = {p['delta_chi2_line_minus_piecewise']:.1f}")
    m = r["moresco2016_replication"]
    print(f"Moresco 2016 replication: oLCDM z_t = {m['olcdm_profile']['best']:.2f} "
          f"(68% {m['olcdm_profile']['68']}), piecewise {m['piecewise_profile']['best']:.2f} "
          f"(68% {m['piecewise_profile']['68']}, 95% {m['piecewise_profile']['95']})")
    for label, sb in r["scaling_bias"].items():
        print(f"scaling bias ({label}, truth H0 = {sb['truth_H0']:.1f}): " + ", ".join(
            f"{k}: {v['bias']:+.2f} ± {v['std'] / np.sqrt(sb['n_mock']):.2f}"
            for k, v in sb.items() if isinstance(v, dict)))
    d = r["desi"]
    print(f"DESI DR2 H(z) with D3's r_d: chi2 vs CC GP = {d['chi2_vs_cc_gp']:.1f} / {d['n_vs_cc_gp']};"
          f" vs D3 LCDM {d['chi2_vs_D3']:.1f}, vs D1 LCDM {d['chi2_vs_D1']:.1f} (n = {d['n']})")
    for k, v in r["clock"]["moresco2012_reread"].items():
        print(f"Moresco 2012 re-read ({k}):", ", ".join(
            f"H({c['z']:.3f}) = {c['H']:.1f} ± {c['stat']:.1f} ± {c['syst']:.1f}" for c in v))
    fut = r["future"]
    for name, e in fut.items():
        if isinstance(e, dict) and "age_gyr" in e:
            extra = ""
            if e.get("big_rip_from_now_gyr"):
                extra = f"  Big Rip in {e['big_rip_from_now_gyr']:.1f} Gyr"
            if e.get("acceleration_end"):
                extra = f"  acceleration ends in {e['acceleration_end']['t_from_now_gyr']:.1f} Gyr"
            print(f"future: {name:62s} age {e['age_gyr']:.2f} Gyr{extra}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--figures-only", action="store_true")
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    t0 = time.time()
    if args.figures_only and CACHE.exists():
        results, art = pickle.loads(CACHE.read_bytes())
    else:
        import analysis

        results, art = analysis.run_all()
        cfg.RESULTS_PATH.write_text(json.dumps(_jsonable(results), indent=1))
        CACHE.write_bytes(pickle.dumps((results, art)))
        print(f"analysis: {time.time() - t0:.0f} s")
    summary(results)
    if not args.no_figures:
        import figures

        written = figures.make_all(results, art)
        print(f"\n{len(written)} figure files written to figures/")
    print(f"total: {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
