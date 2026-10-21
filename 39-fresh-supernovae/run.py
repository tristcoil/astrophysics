"""Run D6 end to end: fits (cached), analysis, figures, explorer, summary.

    python run.py                 # everything (first run: ~15-25 min of SALT fits)
    python run.py --figures-only  # redraw figures and the explorer from results.json
    python run.py --workers 4     # fewer parallel SALT fits on a busy machine

Needs `python fetch_data.py` once (and `python fetch_live.py` for the live
points; without a snapshot the live part is skipped). No network after that.
"""

from __future__ import annotations

import argparse
import json
import sys
import time

import config


def summary(R: dict) -> str:
    s, land, sc = R["sample"], R["landing"], R["scatter"]
    cc = R["convention_check"]["fixed_x1_c_t0"]
    lines = [
        f"ZTF DR2: {s['n_released']} released, {s['n_basic_cuts']} pass the basic cuts with our "
        f"SALT3 fits (ZTF team: {s['published']['n_basic_cuts']}), {s['n_cosmo']} normal and new to D1, "
        f"{s['n_d1_test']} also inside D1's colour/stretch range.",
        f"D1's recipe: alpha {R['d1_recipe']['alpha']:.4f}, beta {R['d1_recipe']['beta']:.4f} "
        f"(rebuilds m_b_corr to {R['d1_recipe']['regression_rms']:.0e} mag); sncosmo vs SNANA "
        f"amplitude convention {cc['mean']:.3f} +- {cc['err']:.3f} mag (MAG_OFFSET 0.27).",
    ]
    for flavour in ("nominal", "full"):
        o = land[flavour]["offset_hubble_flow"]
        tz = land[flavour]["trend_vs_z_0.015-0.06"]
        lines.append(f"  D1 {flavour:7s}: offset from D1's curve {o['mean']:+.3f} +- {o['err']:.3f} mag "
                     f"(0.023-0.06, n={o['n']}); trend with z {tz['slope']:+.2f} +- {tz['slope_err']:.2f} mag/unit z")
    lines.append(f"  (naive SALT3 + 10.635 convention: {land['naive_salt3_offset']['mean']:+.3f} mag)")
    c = sc["dr2_clean_0.03-0.1"]
    lines.append(f"Scatter 0.03 < z < 0.1 (nMAD): unstandardised {c['raw_salt3']['nmad']:.3f}, ZTF-fitted "
                 f"{c['ztf_fit_salt3']['nmad']:.3f}, D1's recipe {c['d1_recipe_b21']['nmad']:.3f} mag "
                 f"(ZTF team: 0.33 -> 0.150)")
    zf = R["ztf_standardisation"]["with_mass_step"]
    lines.append("ZTF's own standardisation: " + ", ".join(
        f"{k} {zf['params'][k]:.3f} +- {zf['errors'][k]:.3f}" for k in ("alpha", "beta", "gamma"))
        + f", sigma_int {zf['sigma_int']:.3f}")
    b = R["scatter_budget"]["galaxy_redshifts"]
    lines.append(f"Random velocities: sigma_v = {b['sigma_v']:.0f} km/s ({b['sigma_v_68'][0]:.0f}-"
                 f"{b['sigma_v_68'][1]:.0f}), sigma_SN = {b['sigma_sn']:.3f} mag")
    fl = R.get("flows")
    if fl:
        d = fl["dipole"]
        lines.append(f"Bulk flow: {d['dipole_kms']:.0f} +- {d['dipole_err']:.0f} km/s towards (l, b) = "
                     f"({d['l']:.0f}, {d['b']:.0f}), dchi2 {d['delta_chi2_vs_offset']:.1f}; sigma_v after "
                     f"removing flows {fl['sigma_v_after_removing_flow']['sigma_v']:.0f} km/s")
    lv = R.get("live", {})
    if lv.get("available"):
        lines.append(f"Live snapshot {lv['meta'].get('snapshot')}: {lv['n_objects']} fitted, "
                     f"{lv['n_d1_test']} pass the cuts"
                     + (f", offset {lv['offset']['mean']:+.3f} +- {lv['offset']['err']:.3f}"
                        if "offset" in lv else ""))
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--figures-only", action="store_true")
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()
    t0 = time.time()
    if args.figures_only:
        R = json.loads(config.RESULTS_PATH.read_text())
    else:
        import analysis

        R = analysis.run(n_workers=args.workers)
    import figures
    import explorer

    written = figures.make_all(R)
    written.append(explorer.build(R))
    print(summary(R))
    print(f"\n{len(written)} files written in {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
