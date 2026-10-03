"""Run the whole project end to end.

    python run.py                # everything the cached data allow
    python run.py --no-maps      # spectra only, even if the maps are present
    python run.py --no-fit       # skip the CAMB model curves and parameter fit
    python run.py --recompute    # ignore cached intermediate results

Needs ``python fetch_data.py`` first. With the cache warm, no network access
happens at all. The expensive steps -- the map spectra (~2 min), the simulated-
sky check (~1 min), the parameter fit (~2 min), the peak error bars (~1.5 min)
-- are cached in ``data/derived/``, so a second run takes seconds.
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from datetime import datetime, timezone

import numpy as np

warnings.filterwarnings("ignore")

import cmb_physics as cp  # noqa: E402
import spectra  # noqa: E402
from config import (DERIVED_DIR, FIGURE_DIR, H0_D1, H0_SH0ES, PRODUCTS,  # noqa: E402
                    PUBLISHED_TT_PEAKS, RESULTS_PATH, T_CMB_ERR_K, T_CMB_K)

PEAKS_CACHE = DERIVED_DIR / "peaks.json"
N_PEAKS = 7  # the 8th is marginal even for Planck, and unstable under our method


def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist() if obj.size < 5000 else f"<array of {obj.size}>"
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    return obj


def run_peaks(full, theory, *, force: bool) -> dict:
    if PEAKS_CACHE.exists() and not force:
        return json.loads(PEAKS_CACHE.read_text())
    print("  measuring peaks, with Monte Carlo error bars (cached afterwards)")
    data = spectra.measure_peaks(full, n_peaks=N_PEAKS)
    model = spectra.measure_theory_peaks(theory, n_peaks=N_PEAKS)
    harmonic = spectra.acoustic_scale(data)
    # The same series fitted to the model's peaks, with the data's error bars as
    # weights so both fits weigh the peaks identically.
    harmonic_model = spectra.acoustic_scale(
        [dict(m, ell_err=d["ell_err"]) for m, d in zip(model, data)])
    out = {"data": data, "model": model, "harmonic": harmonic, "harmonic_model": harmonic_model}
    PEAKS_CACHE.write_text(json.dumps(out, indent=2))
    return out


def curvature_for_h0(fit: dict, h0_values) -> list[dict]:
    """Curvature implied by our measured acoustic angle, for a range of assumed H0.

    The peak heights fix the physical densities (hence the sound horizon); the
    peak positions fix the angle it subtends. Distance to last scattering
    depends on both curvature and H0, so one angle fixes one combination --
    the CMB's "geometric degeneracy". This solves for the curvature that
    reproduces our angle under each assumed H0.
    """
    from scipy.optimize import brentq

    p, d = fit["parameters"], fit["derived"]
    ombh2, omch2 = p["ombh2"]["value"], p["omch2"]["value"]
    z_star = d["zstar"]["value"]
    target = d["thetastar100"]["value"] / 100.0

    def theta(omk, h):
        cosmo = cp.make_cosmology(h=h, ombh2=ombh2, omch2=omch2, omk=omk)
        return cp.angular_scale(cosmo, z_star, ombh2)["theta_star"]

    rows = []
    for h0 in h0_values:
        h = h0 / 100.0
        try:
            omk = brentq(lambda k: theta(k, h) - target, -0.4, 0.4, xtol=1e-6)
        except ValueError:
            omk = float("nan")
        rows.append({"H0": float(h0), "omk": float(omk)})
    return rows


def print_summary(res: dict) -> None:
    line = "=" * 80
    print(f"\n{line}\nRESULTS\n{line}")

    f = res["firas"]
    print("\nThe blackbody (COBE/FIRAS)")
    print(f"  T = {f['temperature_k']:.5f} K   (stat. ±{f['temperature_stat_err_k'] * 1e6:.0f} µK; "
          f"FIRAS calibration ±1 mK; Fixsen 2009: {T_CMB_K} ± {T_CMB_ERR_K} K)")
    print(f"  chi2 = {f['blackbody']['chi2']:.1f} for {f['blackbody']['dof']} dof; "
          f"Compton y = ({f['y'] * 1e6:+.2f} ± {f['y_err'] * 1e6:.2f})e-6; "
          f"peak {f['peak_ghz']:.1f} GHz")

    print(f"\nPlanck binning reproduced: max |ours / Planck's binned best fit - 1| = "
          f"{res['binning_check']['max_rel_diff']:.1e}")

    pk = res["peaks"]
    print("\nAcoustic peaks (Planck 2018 TT)")
    print(f"  {'':>3} {'ours: l':>16} {'D_l (uK^2)':>16}  {'model, same method':>20}  "
          f"{'Planck 2015 table':>22}")
    pub = {r[0]: r for r in PUBLISHED_TT_PEAKS}
    for d, m in zip(pk["data"], pk["model"]):
        r = pub[d["n"]]
        print(f"  {d['n']:>3} {d['ell']:>9.1f} ± {d['ell_err']:<4.1f} {d['dl']:>8.1f} ± "
              f"{d['dl_err']:<5.1f}  {m['ell']:>8.1f} {m['dl']:>9.1f}   "
              f"{r[1]:>8.1f} ± {r[2]:<4} {r[3]:>7.1f}")
    h, hm = pk["harmonic"], pk["harmonic_model"]
    print(f"  harmonic series l_n = l_A (n - phi): l_A = {h['ell_A']:.1f} ± {h['ell_A_err']:.1f}, "
          f"phi = {h['phi']:.3f} ± {h['phi_err']:.3f}  (model through same fit: "
          f"l_A = {hm['ell_A']:.1f}, phi = {hm['phi']:.3f})")
    print(f"  peak ratios {' : '.join(f'{x:.2f}' for x in h['ratios_to_first'])}"
          f"   (isocurvature would give 1 : 3 : 5 : 7)")

    if "maps" in res:
        m = res["maps"]
        print("\nOur spectrum from the sky maps (SMICA half-missions)")
        print(f"  sky fraction used {m['fsky_binary']:.3f}; map rms {m['rms_sum_uK']:.1f} µK "
              f"(noise {m['rms_diff_uK']:.1f} µK); predicted from the best-fit spectrum "
              f"{m['predicted_signal_rms_uK']:.1f} µK signal")
        print(f"  dT/T = {m['delta_t_over_t']:.1e}")
        for lo, hi in ((30, 500), (500, 1000), (1000, 1500), (1500, 2000), (2000, 2508)):
            c = m["comparison"]
            print(f"    l {lo:>4}-{hi:<4}: ours / Planck = {c[f'mean_ratio_{lo}_{hi}']:.4f}   "
                  f"chi2 vs best fit {c[f'chi2_{lo}_{hi}']:6.1f} / {c[f'n_{lo}_{hi}']}")
        c = m["comparison"]
        print(f"  high-l excess as point sources: D = A (l/3000)^2, A = {c['poisson_amp']:.1f} ± "
              f"{c['poisson_amp_err']:.1f} µK² (chi2 {c['poisson_chi2']:.0f} / {c['poisson_dof']})")
        v = m["validation"]
        print(f"  simulated sky: area-only correction biased {v['fsky_mean_bias_30_500'] * 100:+.1f}% "
              f"(l 30-500, worst bin {v['fsky_worst_bin'] * 100:+.1f}%); correlation method "
              f"{v['corr_mean_bias_30_2508'] * 100:+.2f}% mean, rms {v['corr_rms_30_2508'] * 100:.2f}%")

    if "fit" in res:
        fit = res["fit"]
        print(f"\nLCDM fit to the binned TT spectrum (CAMB; tau fixed at {fit['tau_fixed']})")
        print(f"  chi2 = {fit['chi2']:.1f} for {fit['dof']} dof   "
              f"(CAMB reproduces Planck's best-fit curve to "
              f"{res['camb_check']['rms_rel_diff'] * 100:.2f}% rms)")
        print(f"    {'':<14}{'ours':>22}{'Planck TT+lowE':>24}{'shift':>9}")
        for row in res["fit_comparison"]:
            print(f"    {row['name']:<14}{row['ours']:>13.5f} ± {row['ours_err']:<8.5f}"
                  f"{row['planck']:>13.5f} ± {row['planck_err']:<8.5f}"
                  f"{row['shift_sigma']:>+7.2f}σ")
        h0 = fit["derived"]["H0"]
        t = (H0_SH0ES[0] - h0["value"]) / np.hypot(H0_SH0ES[1], h0["error"])
        print(f"  H0 = {h0['value']:.2f} ± {h0['error']:.2f} km/s/Mpc vs SH0ES "
              f"{H0_SH0ES[0]} ± {H0_SH0ES[1]}: {t:.1f} sigma apart (the Hubble tension)")
        print(f"  ... and vs project D1's supernovae, {H0_D1[0]} ± {H0_D1[1]}: "
              f"{res['H0_tension_sigma_D1']:.1f} sigma apart")
        print("  curvature implied by our acoustic angle, for an assumed H0:")
        print("    " + ",  ".join(f"H0 {r['H0']:.1f} -> Omega_k {r['omk']:+.3f}"
                                  for r in res["curvature_vs_h0"]))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--no-maps", action="store_true", help="skip the sky-map analysis")
    parser.add_argument("--no-fit", action="store_true", help="skip CAMB (fit and model curves)")
    parser.add_argument("--no-figures", action="store_true", help="skip figure generation")
    parser.add_argument("--recompute", action="store_true",
                        help="recompute cached intermediate results")
    args = parser.parse_args(argv)

    missing = [p.filename for p in PRODUCTS if p.group != "maps" and not p.path.exists()]
    if missing:
        print(f"missing data: {missing}\nrun `python fetch_data.py` first", file=sys.stderr)
        return 1

    res: dict = {"_generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                 "_data_files": {p.filename: datetime.fromtimestamp(p.path.stat().st_mtime,
                                                                    timezone.utc).date().isoformat()
                                 for p in PRODUCTS if p.path.exists()}}

    print("[1] The blackbody")
    firas = spectra.load_firas()
    res["firas"] = firas_fit = spectra.fit_firas(firas)

    print("[2] Planck's power spectrum")
    full = spectra.load_tt_full()
    binned, bf_b = spectra.load_tt_binned()
    theory = spectra.load_theory()
    bestfit = spectra.load_bestfit_parameters()
    _, rebinned = cp.planck_bin(theory.ell, theory.tt)
    res["binning_check"] = {"max_rel_diff": float(np.max(np.abs(rebinned / bf_b - 1)))}
    res["planck_bestfit_chi2"] = spectra.binned_chi2(binned, bf_b)

    print("[3] Acoustic peaks")
    res["peaks"] = pk = run_peaks(full, theory, force=args.recompute)
    hi = full.select(30)
    envelope = spectra.damping_envelope(hi.ell.astype(float), hi.dl, hi.err)

    from figures import Context
    ctx = Context(firas=firas, firas_fit=firas_fit, full=full, binned=binned,
                  binned_bestfit=bf_b, theory=theory, bestfit=bestfit, peaks=pk["data"],
                  model_peaks=pk["model"], harmonic=pk["harmonic"], envelope=envelope)

    from fetch_data import have_maps
    if args.no_maps:
        print("[4] Sky maps: skipped (--no-maps)")
    elif not have_maps():
        print("[4] Sky maps: not downloaded -- run `python fetch_data.py` to include them")
    else:
        print("[4] Sky maps")
        import sky
        s = sky.compute_map_spectra(force=args.recompute)
        ctx.maps = s
        ctx.display = sky.compute_display_maps(force=args.recompute)
        ctx.validation = v = sky.validate_estimator(theory.ell, theory.tt, force=args.recompute)
        ctx.map_compare = comp = sky.compare_to_planck(s, binned, bf_b)
        lv, full_b = sky.bin_cl(v["ell"], v["cl_full"])
        _, simple_b = sky.bin_cl(v["ell"], v["cl_masked_fsky"])
        _, corr_b = sky.bin_cl(v["ell"], v["cl_masked"])
        low = lv < 500
        signal_rms = sky.predicted_rms(theory.ell, theory.tt, s["beam"])
        res["maps"] = {
            **{k: float(s[k]) for k in ("fsky_binary", "fsky_common", "fsky_eff", "w2",
                                         "rms_sum_uK", "rms_diff_uK", "rms_sum_fullsky_uK",
                                         "min_sum_uK", "max_sum_uK", "max_abs_fullsky_uK")},
            "dipole_removed_uK": np.linalg.norm(s["dipole_removed"], axis=1).tolist(),
            "predicted_signal_rms_uK": signal_rms,
            "measured_signal_rms_uK": float(np.sqrt(s["rms_sum_uK"] ** 2 - s["rms_diff_uK"] ** 2)),
            "delta_t_over_t": sky.delta_t_over_t(float(s["rms_sum_uK"])),
            "comparison": comp,
            "validation": {
                "fsky_mean_bias_30_500": float(np.mean(simple_b[low] / full_b[low]) - 1),
                "fsky_worst_bin": float(np.max(np.abs(simple_b / full_b - 1))),
                "corr_mean_bias_30_2508": float(np.mean(corr_b / full_b) - 1),
                "corr_rms_30_2508": float(np.sqrt(np.mean((corr_b / full_b - 1) ** 2))),
                "fsky_ratio": (simple_b / full_b).tolist(),
                "corr_ratio": (corr_b / full_b).tolist(),
            },
        }

    if args.no_fit:
        print("[5] CAMB fit: skipped (--no-fit)")
    else:
        try:
            import cosmology_fit as cf
        except ImportError:
            cf = None
            print("[5] CAMB not installed -- skipping the fit (pip install camb)")
        if cf is not None:
            print("[5] Model spectra and the parameter fit (CAMB)")
            res["camb_check"] = cf.reproduce_bestfit(bestfit, theory)
            ctx.responses = cf.parameter_responses(bestfit, force=args.recompute)
            ctx.fit = fit = cf.fit_tt(binned, force=args.recompute)
            res["fit"] = fit
            res["fit_comparison"] = cf.compare_to_planck(fit)
            res["curvature_vs_h0"] = curvature_for_h0(
                fit, [55.0, 60.0, 65.0, fit["derived"]["H0"]["value"], 70.0, H0_SH0ES[0], 80.0])
            h0 = fit["derived"]["H0"]
            res["H0_tension_sigma"] = float(
                (H0_SH0ES[0] - h0["value"]) / np.hypot(H0_SH0ES[1], h0["error"]))
            res["H0_tension_sigma_D1"] = float(
                (H0_D1[0] - h0["value"]) / np.hypot(H0_D1[1], h0["error"]))

    if not args.no_figures:
        print("[6] Figures")
        import figures
        written = figures.make_all(ctx)
        print(f"  wrote {len(written)} files to {FIGURE_DIR.name}/")

    RESULTS_PATH.write_text(json.dumps(_jsonable(res), indent=2))
    print_summary(res)
    print(f"\nFull results in {RESULTS_PATH.relative_to(RESULTS_PATH.parent.parent)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
