"""Run every fit, write data/results.json, and draw every figure.

    python run.py                 # analysis + figures (5-12 minutes; ~35 on a busy machine)
    python run.py --figures-only  # redraw from the cached analysis
    python run.py --no-figures    # numbers only
    python run.py --photometry    # re-measure the ZTF/PTF photometry first (~5 min)
    python run.py --animation     # also re-render the orbit movie (~2 min; made once otherwise)

Needs the cache from ``python fetch_data.py``; after that it runs offline. The
photometry tables are cached in data/derived/ and reused unless --photometry is given.
"""

from __future__ import annotations

import argparse
import json
import pickle
import time
import warnings

import numpy as np

import config as C

CACHE = C.DATA_DIR / "analysis_cache.pkl"


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


def summary(r: dict) -> None:
    print("\n=== ZTF J1539+5027 and LISA's white-dwarf pairs ===")
    p = r["ztf_photometry"]
    print(f"ZTF: {p['n_after_basic_cuts']} of {p['n_exposures_with_cutouts']} exposures used {p['per_band']}; "
          f"PTF: {r['ptf_photometry']['n_after_cuts']} images {r['ptf_photometry']['per_year']}")
    d = r["discovery"]
    print(f"discovery (2018 only): P = {d['P_best_s']:.4f} s, power {d['power']:.2f}, FAP {d['fap']:.1e}")
    print(f"season-frequency drift alone: P-dot = {r['frequency_drift_only']['Pdot']:.3e}")
    for key in ("ztf_only", "joint_ztf_ptf"):
        e = r[key]
        print(f"{key:14s}: P(T_REF) = {e['P_s']:.7f} +- {e['P_err_s']:.7f} s, P-dot = ({e['Pdot'] * 1e11:.4f} +- "
              f"{e['Pdot_err'] * 1e11:.4f})e-11, T_ecl = {e['T_ecl_bjd']:.7f} +- {e['T_ecl_err_s']:.2f} s, "
              f"chi2 {e['chi2']:.0f}/{e['dof']}")
    for k, v in r["ztf_variants"].items():
        if "Pdot" in v:
            print(f"  variant {k:16s}: P-dot = ({v['Pdot'] * 1e11:.4f} +- {v['Pdot_err'] * 1e11:.4f})e-11 "
                  f"(n = {v['n_points']})")
    cu = r["joint_cubic"]
    print(f"  cubic term (ZTF+PTF): f2 = {cu['f2_per_day3']:.2e} +- {cu['f2_err']:.2e} /d^3 (GR {cu['f2_gr_prediction']:.1e}); "
          f"braking index {cu['braking_index']:.0f} +- {cu['braking_index_err']:.0f} (GR: 11/3); delta chi2 {cu['delta_chi2']:.1f}")
    for q in r["ptf_vs_ztf_prediction"]:
        print(f"  PTF {q['year']}: phase vs ZTF-only prediction {q['dphi']:+.3f} +- {q['dphi_err']:.3f} cycles "
              f"(prediction uncertainty {q['prediction_err']:.3f}; one orbit = dP/dt off by {q['pdot_change_for_one_cycle']:.2e})")
    for name, c in r["literature_check_joint"].items():
        print(f"vs {name}: dT0 = {c['dT0_s']:+.2f} +- {c['dT0_err_s']:.2f} s ({c['dT0_sigma']:+.1f} sig); "
              f"P {c['P_ours_s']:.7f} vs {c['P_theirs_s']:.7f} ({c['dP_sigma']:+.1f} sig); "
              f"P-dot {c['dPdot_sigma']:+.1f} sig")
    cm = r["chirp_mass"]
    print(f"chirp mass (GR only) {cm['Mc_gr_only']:.4f} +- {cm['Mc_gr_only_err']:.4f}; from Burdge masses "
          f"{cm['Mc_from_burdge_masses']:.4f}; observed/GR(masses) = {cm['ratio_obs_over_gr']:.3f} +- "
          f"{cm['ratio_err']:.3f} (tides predicted {1 + cm['tide_prediction']:.3f})")
    ls = r["lisa_summary"]
    print(f"LISA: {ls['n_binaries']} binaries; SNR>7 in 3 months: {ls['n_snr3mo_gt7']} "
          f"({ls['n_verification_with_snr3mo_gt7']} of Kupfer's {ls['n_verification_kupfer']} verification); "
          f"SNR>7 in 4 yr: {ls['n_snr4yr_gt7']}")
    for label, p in r["population"].items():
        if "n_in_band" in p:
            print(f"population {label}: {p['n_wd_galaxy']:.3g} WDs, {p['n_in_band']:.3g} pairs in band, "
                  f"{p['n_resolved_snr7_4yr']:,} resolved (SNR>7, 4 yr), foreground edge "
                  f"{p['foreground_edge_Hz'] * 1e3:.2f} mHz, merger rate per WD {p['merger_rate_per_wd_yr']:.2e}/yr")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--figures-only", action="store_true")
    parser.add_argument("--no-figures", action="store_true")
    parser.add_argument("--photometry", action="store_true")
    parser.add_argument("--animation", action="store_true", help="re-render figures/j1539_orbit.mp4 (~2 min)")
    args = parser.parse_args()
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    t0 = time.time()
    if args.photometry:
        import ptf_phot
        import ztf_phot

        ztf_phot.run(force=True)
        ptf_phot.run(force=True)
        print(f"photometry: {time.time() - t0:.0f} s")
    import analysis

    if args.figures_only and CACHE.exists():
        results, art = pickle.loads(CACHE.read_bytes())
        analysis.post_process(results, art)
    else:
        results, art = analysis.run_all()
        print(f"analysis: {time.time() - t0:.0f} s")
    C.RESULTS_PATH.write_text(json.dumps(_jsonable(results), indent=1))
    CACHE.write_bytes(pickle.dumps((results, art)))
    summary(results)
    if not args.no_figures:
        import figures

        written = figures.make_all(results, art, animation_too=args.animation)
        print(f"\n{len(written)} figure files written to figures/")
    print(f"total: {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
