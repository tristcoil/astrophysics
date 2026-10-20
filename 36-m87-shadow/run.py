"""Run the whole project end to end.

    python fetch_data.py     # once: ~11 MB from the EHT's GitHub data releases
    python run.py            # every fit, image, number and figure, plus the explorer page
    python test_physics.py   # 18 self-checks

The first run performs about 40 nested-sampling fits, checks each for a missed
mode (redoing any fit that another data set's solution beats), and images four
nights twice (our imager and the EHT's own script): several hours on a few
cores, most of it the ring fits. Everything is cached in data/cache/, so a re-run
takes a few minutes and needs no network. ``--refit`` starts the fits over.
Every headline number is printed next to the published value it is graded against.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import warnings

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
warnings.filterwarnings("ignore")

import analysis  # noqa: E402
from config import M87_2017_DAYS, PUBLISHED, RESULTS_PATH  # noqa: E402


def _rule(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def pm(s: dict, d: int = 2) -> str:
    return f"{s['median']:.{d}f} (+{s['hi68'] - s['median']:.{d}f} / -{s['median'] - s['lo68']:.{d}f})"


def print_summary(res: dict) -> None:
    _rule("THE DATA (M87, April 2017, EHT-HOPS, scan-averaged)")
    for k, d in res["data_2017"].items():
        print(f"  {k}: {d['n_scan_avg']} visibilities in {d['n_scans']} scans from {d['n_10s']} 10-s records; "
              f"{d['n_cphase']} closure phases, {d['n_lcamp']} log closure amplitudes; longest baseline "
              f"{d['max_baseline_glambda']:.1f} Glambda; {d['trivial_cphase_n']} trivial closure phases: median "
              f"|psi| {d['trivial_cphase_median_abs_deg']:.1f} deg, chi2/N about zero {d['trivial_cphase_chi2_per_n']:.2f}")

    _rule("THE RING, NIGHT BY NIGHT (closure-only fits, nested sampling)")
    for k, e in res["m87_2017_fits"].items():
        f, x = e.get("fiducial"), e.get("xsring")
        p6 = e["paper6_xsring"]
        line = f"  {e['label']:26s}"
        if f:
            line += (f" ring d = {f['d']['median']:5.2f} uas, PA {f['pa_deg']['median']:5.1f}, "
                     f"chi2/N {f['reduced_chi2']['cphase']:.2f}/{f['reduced_chi2']['lcamp']:.2f}")
        if x:
            line += f" | crescent d = {x['d']['median']:5.2f} (Paper VI {p6['d'][0]:.1f}), PA {x['pa_deg']['median']:5.1f} ({p6['pa'][0]:.1f})"
        print(line)
    if "m87_2017_ring" in res:
        r = res["m87_2017_ring"]["d"]
        print(f"  blurred-ring (delta-ring) diameter, 8 data sets: {r['weighted_mean']:.2f} uas, night-to-night rms "
              f"{r['rms_scatter']:.2f}  | EHT 2025 (2017, delta-ring, all methods): 42.0-46.4")
    if "m87_2017_xsring" in res:
        r = res["m87_2017_xsring"]
        print(f"  crescent diameter: {r['d']['weighted_mean']:.2f} uas, rms {r['d']['rms_scatter']:.2f} | Paper VI "
              f"(same model): {PUBLISHED['m87_xsring_mean_uas']} (rms 0.69)")

    if res.get("model_comparison"):
        _rule("IS IT A RING? (April 11, low band; ln Z relative to the best model)")
        for k, v in sorted(res["model_comparison"].items(), key=lambda kv: -kv[1]["logz"]):
            red = v["reduced_chi2"]
            print(f"  {v['name']:34s} {v['n_params']:2d} params  Delta lnZ {v['delta_logz']:9.1f}  chi2/N "
                  f"{red.get('cphase', float('nan')):.2f}/{red.get('lcamp', float('nan')):.2f}")

    if "m87_mass" in res:
        _rule("FROM RING TO MASS (D = 16.8 Mpc)")
        m = res["m87_mass"]
        g = m["grmhd_calibrated"]
        print(f"  GRMHD-calibrated (alpha {g['alpha']}): theta_g = {pm(g['theta_g'])} uas -> M = {pm(g['mass_1e9'])} billion "
              f"suns | EHT 2019: 3.8 +/- 0.4 uas, 6.5 +/- 0.7")
        s = m["shadow_schwarzschild"]
        print(f"  if the ring were the shadow edge (alpha 10.39): theta_g = {pm(s['theta_g'])} -> M = {pm(s['mass_1e9'])}")
        for name, t in m["tension_sigma"].items():
            print(f"  ring vs {name}: {t:+.1f} sigma (distance-free, in theta_g)")

    if "sgra_test" in res:
        _rule("SGR A*: H1'S PREDICTION FROM THE STAR S2, TESTED")
        st = res["sgra_test"]
        print(f"  H1: theta_g = {st['theta_g_h1']:.3f} +/- {st['theta_g_h1_err']:.3f} uas -> shadow "
              f"{st['shadow_prediction_uas']:.2f} +/- {st['shadow_prediction_err']:.2f} uas; ring expected with the "
              f"EHT's calibration (alpha {st['alpha_eht_sgra']:.2f}): {st['ring_prediction_uas']:.1f} uas")
        for fam, label in (("amplitudes", "from the amplitudes beyond 2 Glambda (primary)"),
                           ("closures", "from closure quantities only (degenerate with a halo)")):
            f = st.get(fam)
            if not f or "d_peak" not in f:
                continue
            for b, v in f["best_by_band"].items():
                print(f"  {label}, {b} band ({v['model']}): peak diameter {pm(v['d_peak'])} uas, chi2/N "
                      + ", ".join(f"{k} {x:.2f}" for k, x in v["reduced_chi2"].items()))
            print(f"    combined {f['d_peak']['value']:.1f} +/- {f['d_peak']['err']:.1f} uas (all fits "
                  f"{f['d_peak_all_fits']['min']:.1f}-{f['d_peak_all_fits']['max']:.1f}); measured/predicted "
                  f"{f['measured_over_predicted']['value']:.3f} +/- {f['measured_over_predicted']['err']:.3f} "
                  f"| EHT 2022: 51.8 +/- 2.3")

    if "m87_2018_ring" in res:
        _rule("M87 A YEAR LATER (April 2018)")
        for k, e in res["m87_2018_fits"].items():
            f = e.get("fiducial")
            if f:
                print(f"  {e['label']:26s} d = {f['d']['median']:5.2f} uas, PA {f['pa_deg']['median']:5.1f}")
        r = res["m87_2018_ring"]["d"]
        print(f"  2018 ring {r['weighted_mean']:.2f} uas (rms {r['rms_scatter']:.2f}) | EHT 2024: 43.3 +1.5/-3.1")

    if res.get("image_features"):
        _rule("IMAGES (ring diameter and brightest direction from each image)")
        for day, e in res["image_features"].items():
            parts = [f"{k}: d {v['d']:.1f} uas, PA {v['pa_deg']:.0f}" for k, v in e.items()]
            print(f"  {M87_2017_DAYS[day]:9s} " + " | ".join(parts))
        print("  Paper VI image-domain diameters: eht-imaging 40.5, SMILI 41.5, DIFMAP 42.5 uas")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--refit", action="store_true", help="redo every nested-sampling fit")
    args = parser.parse_args(argv)
    analysis.ensure_dirs()
    for group in analysis.GROUPS:
        print(f"[fits: {group}]", flush=True)
        analysis.run_group(group, refit=args.refit)
    print("[mode check: did each fit find the best solution?]", flush=True)
    analysis.main(["verify"])
    print("[images]", flush=True)
    analysis.main(["images"])
    import validation

    print("[imaging validation on simulated data]", flush=True)
    with contextlib.redirect_stdout(io.StringIO()):
        validation.run("101", settings=analysis.IMAGER_SETTINGS, weights=analysis.IMAGER_WEIGHTS)
    print("[summary]", flush=True)
    res = analysis.summarize()
    import figures
    import viewer

    print("[figures]", flush=True)
    figures.main()
    print("  ", viewer.build(res))
    print_summary(json.loads(RESULTS_PATH.read_text()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
