"""Run the whole project end to end.

    python run.py                 # everything; ~10 min once the profile fits are cached
    python run.py --rebuild       # redo every cached analysis step (profile fits ~25 min)
    python run.py --no-figures

Fetches whatever is not cached, fits every 21-cm line profile of the four
THINGS galaxies (cached after the first run), derives the rotation curves,
runs the cusp-core simulation for DDO 154, fits stellar and gas kinematic
axes for every MaNGA galaxy, writes ``data/results.json``, draws the figures
and prints a summary graded against the published values. Once the cache is
warm it runs offline.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings

# The fits are many tiny problems, vectorized by hand; multithreaded BLAS only
# adds contention (and MaNGA's worker pool multiplies it).
for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402

import analysis as A  # noqa: E402
import literature as L  # noqa: E402
import manga as MG  # noqa: E402
from config import GALAXIES, RESULTS_PATH  # noqa: E402


def _clean(obj):
    """Make analysis output JSON-safe (numpy scalars and arrays, tuples)."""
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [_clean(v) for v in obj.tolist()]
    if isinstance(obj, (np.floating, float)):
        return None if not np.isfinite(obj) else float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    return obj


def things_results(grs, cusp) -> dict:
    out = {}
    for key, gr in grs.items():
        tr = gr.primary
        r = gr.r_kpc()
        c = gr.comparisons
        outer = r >= 0.5 * np.nanmax(r)
        g = next(x for x in GALAXIES if x.key == key)
        out[key] = {
            "name": gr.name, "distance_mpc": gr.distance_mpc, "beam_arcsec": gr.beam[:2],
            "noise_mjy_per_beam": gr.noise_mjy, "channel_kms": gr.channel_kms, "pixel_arcsec": gr.pixel_arcsec,
            "spectra_fitted": gr.n_fitted, "pixels_in_velocity_field": gr.n_good,
            "rings": len(tr.radius), "ring_step_arcsec": gr.step_arcsec, "r_max_kpc": float(np.nanmax(r)),
            "centre_ra_dec": gr.centre_radec, "geometry": c["geometry"], "geometry_note": g.geometry_note,
            "v_outer_half_median": float(np.nanmedian(tr.vrot[outer])), "v_max": float(np.nanmax(tr.vrot)),
            "vs_deblok": c["deblok"], "vs_sparc": c["sparc"],
            "harmonics": gr.harmonic_summary, "hi": {k: v for k, v in gr.flux.items() if "spectrum" not in k},
            "notes": tr.notes,
            "rotation_curve": {"r_kpc": r, "v": tr.vrot, "err": tr.vrot_err, "v_approaching": tr.vrot_app,
                               "v_receding": tr.vrot_rec, "pa": tr.pa, "incl": tr.incl},
        }
    cusp_out = {k: v for k, v in cusp.items() if not k.startswith("rho_") and k != "models"}
    cusp_out["models"] = {name: {k: v for k, v in m.items()} for name, m in cusp["models"].items()}
    cusp_out["oh2011"] = {"alpha": L.OH2011_DDO154_ALPHA, "points": L.oh2011_ddo154_density()}
    return {"galaxies": out, "ddo154_cusp": cusp_out}


def print_things(res: dict) -> None:
    print("\nRotation curves from the THINGS cubes, against de Blok et al. (2008)")
    print(f"  {'galaxy':<10}{'rings':>6}{'R_max':>8}{'V_sys':>8}{'(dB08)':>8}{'<i>':>6}{'(dB08)':>7}{'<PA>':>7}"
          f"{'(dB08)':>8}{'centre':>8}{'dV med':>8}{'rms':>6}{'1 sigma':>8}")
    for key, g in res["galaxies"].items():
        geo, d = g["geometry"], g["vs_deblok"]["gh_vpeak"]
        print(f"  {g['name']:<10}{g['rings']:>6}{g['r_max_kpc']:>6.1f}kpc{geo['vsys']:>8.1f}{geo['vsys_deblok']:>8.1f}"
              f"{geo['incl_mean']:>6.1f}{geo['incl_deblok']:>7.1f}{geo['pa_mean']:>7.1f}{geo['pa_deblok']:>8.1f}"
              f"{geo['centre_offset_arcsec']:>7.1f}\"{d['median_diff']:>+8.2f}{d['rms_diff']:>6.2f}"
              f"{100 * d['within_1sigma']:>7.0f}%")
    print("\n  Which velocity recipe matches de Blok et al.? rms difference (km/s) / inner offset")
    kinds = list(A.VELOCITY_KINDS)
    print("  " + " " * 10 + "".join(f"{A.VELOCITY_KINDS[k][:22]:>24}" for k in kinds))
    for key, g in res["galaxies"].items():
        print(f"  {g['name']:<10}" + "".join(
            f"{g['vs_deblok'][k]['rms_diff']:>14.2f} / {g['vs_deblok'][k]['inner_median_diff']:>+6.1f}" for k in kinds))
    print("\n  Against SPARC's curves from older telescopes (project D4), radii at the THINGS distance")
    for key, g in res["galaxies"].items():
        s = g["vs_sparc"]
        print(f"  {g['name']:<10} {s['n']:>3} points: median {s['median_diff']:+.2f}, rms {s['rms_diff']:.2f} km/s "
              f"({s['refs']})")
    print("\n  Non-circular motions, median A_r (km/s): ours (RESWRI-style rings) vs Trachternach et al. 2008")
    for key, g in res["galaxies"].items():
        h, t = g["harmonics"]["reswri_rings"], g["harmonics"]["trachternach"]
        print(f"  {g['name']:<10} whole disc {h['A_r_median']:5.2f} vs {t['A_r']:5.2f};  inner kpc "
              f"{(h['A_r_1kpc_median'] or float('nan')):5.2f} vs {t['A_r_1kpc']:5.2f};  eps_pot "
              f"{h['eps_pot_mean']:+.3f} vs {t['eps_pot']:+.3f} +- {t['eps_pot_err']:.3f}")
    print("\n  HI: flux in the THINGS moment-0 map vs Walter et al.'s global spectrum")
    for key, g in res["galaxies"].items():
        hi = g["hi"]
        print(f"  {g['name']:<10} {hi['mom0_flux_jykms']:7.1f} vs {hi['walter_flux_jykms']:7.1f} Jy km/s "
              f"({100 * (hi['mom0_flux_jykms'] / hi['walter_flux_jykms'] - 1):+.1f}%);  M_HI = "
              f"{hi['m_hi_msun'] / 1e9:.2f}e9 Msun (SPARC, rescaled: {hi['m_hi_sparc_msun'] / 1e9:.2f}e9)")

    c = res["ddo154_cusp"]
    print("\nDDO 154: cusp or core?")
    for lab, k in (("our curve", "alpha_ours"), ("our curve, gas subtracted", "alpha_ours_dm"),
                   ("de Blok's curve, our code", "alpha_deblok")):
        print(f"  inner density slope, {lab:<26} alpha = {c[k]['alpha']:+.2f} +- {c[k]['alpha_err']:.2f}")
    print(f"  Oh et al. (2011), de Blok's curve:           alpha = {c['oh2011']['alpha'][0]:+.2f} +- "
          f"{c['oh2011']['alpha'][1]:.2f}")
    m = c["models"]["LCDM cusp (NFW)"]
    print(f"  a cold-dark-matter cusp (NFW, c = {c['lcdm_nfw']['c']:.1f}, V200 = {c['lcdm_nfw']['v200']:.1f} km/s): "
          f"alpha {m['alpha_input']['alpha']:+.2f} as built, {m['alpha_recovered']['alpha']:+.2f} after the beam "
          f"and our pipeline")
    ex = c["cusp_excess"]
    print(f"  inside 1 kpc the smeared cusp spins {ex['mean_excess_kms']:.1f} km/s ({ex['mean_ratio']:.2f}x) faster "
          f"than DDO 154: chi2 = {ex['chi2']:.0f} for {ex['n']} points")
    iso, nfw = c["halo"]["ours"]["ISO"], c["halo"]["deblok"]["NFW"]
    print(f"  cored halo: rho0 = {iso['params'][0] * 1e3:.1f} +- {iso['errors'][0] * 1e3:.1f} e-3 Msun/pc^3, R_c = "
          f"{iso['params'][1]:.2f} +- {iso['errors'][1]:.2f} kpc (Oh et al.: 34.5 +- 1.8, 1.30 +- 0.04)")
    print(f"  NFW on de Blok's curve: c = {nfw['params'][0]:.1f} +- {nfw['errors'][0]:.1f}, V200 = "
          f"{nfw['params'][1]:.1f} +- {nfw['errors'][1]:.1f} (Oh et al.: 5.3 +- 0.5, 58.3 +- 4.2)")


def print_manga(ms: dict) -> None:
    cen = ms["census"]
    print("\nMaNGA DR17: stars against gas")
    print(f"  {cen['n_streamed']:,} galaxies streamed; {cen['n_emission']:,} with emission lines; "
          f"{cen['n_measured']:,} with both axes measured; {cen['n_regular']:,} regular rotators")
    for sample in ("measured", "regular"):
        s = cen["by_sample"][sample]
        print(f"  {sample:<9} misaligned (>30 deg) {s['n_misaligned']:>4} = {100 * s['frac_misaligned']:.1f}%, "
              f"counter-rotating (>150 deg) {s['n_counter']:>4} = {100 * s['frac_counter']:.1f}%, "
              f"perpendicular {s['n_perpendicular']}")
        print("            by colour: " + ", ".join(
            f"{c} {100 * v['frac_misaligned']:.1f}% ({v['misaligned']}/{v['n']})" for c, v in s["by_colour"].items()))
    print("  Zhou et al. (2022): 456/6958 = 6.6% misaligned; star-forming 2.4%, green valley 7.3%, quiescent 15.6%")
    print("  Davis et al. (2011, ATLAS3D): 36 +- 5% of fast-rotating early types with ionised gas are misaligned")
    for g in ms["drift"]:
        rows = [r for r in g["rows"] if 0.8 <= r["r_re"] <= 1.5 and "v_star" in r and "v_gas" in r
                and "v_star_corrected" in r]
        if rows:
            lag = np.mean([r["v_gas"] - r["v_star"] for r in rows])
            pred = np.mean([r["v_star_corrected"] - r["v_star"] for r in rows])
            s_imp = np.mean([r["sigma_r_implied"] for r in rows])
            s_meas = np.mean([r["sigma_r"] for r in rows])
            print(f"  asymmetric drift, {g['plateifu']} (0.8-1.5 Re): stars lag gas by {lag:.0f} km/s; the Jeans "
                  f"equation with MaNGA's dispersions predicts {pred:.0f}; the lag implies sigma_R = {s_imp:.0f} "
                  f"km/s against {s_meas:.0f} measured")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--rebuild", action="store_true", help="redo every cached analysis step")
    parser.add_argument("--no-figures", action="store_true")
    parser.add_argument("--skip-manga", action="store_true", help="THINGS half only")
    args = parser.parse_args(argv)
    t0 = time.time()

    import fetch_data
    stages = [s for s in ("refs", "things", "manga") if fetch_data.missing(s)]
    if stages:
        print("Fetching missing data:", ", ".join(stages))
        fetch_data.main(sum((["--stage", s] for s in stages), []))

    print("THINGS: profile fits, tilted rings, harmonics, gas")
    grs = A.galaxy_results(rebuild=args.rebuild)
    print("DDO 154: density slopes, halo fits, beam-smearing simulation")
    cusp = A.cusp_test(grs["DDO154"], rebuild=args.rebuild)
    results = things_results(grs, cusp)
    print_things(results)

    ms = None
    if not args.skip_manga:
        print("\nMaNGA: kinematic position angles of stars and gas")
        ms = MG.run(rebuild=args.rebuild)
        results["manga"] = {"census": ms["census"], "drift": ms["drift"]}
        print_manga(ms)

    RESULTS_PATH.write_text(json.dumps(_clean(results), indent=1))
    if not args.no_figures:
        import figures as F
        import interactive
        print("\nFigures")
        written = [F.fig_spider(grs), F.fig_recipes(grs), F.fig_tilted_rings(grs), F.fig_rotation_curves(grs),
                   F.fig_noncircular(grs), F.fig_cusp_core(grs, cusp)]
        if ms:
            written += [F.fig_drift(ms), F.fig_counter_rotation(ms), F.fig_scorecard(results)]
        written.append([interactive.build(grs)])
        for w in written:
            print("  " + w[0].split("/")[-1].rsplit(".", 1)[0])
    print(f"\nDone in {(time.time() - t0) / 60:.1f} min. Numbers in {RESULTS_PATH.name}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
