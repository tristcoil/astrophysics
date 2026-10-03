"""Run the whole project end to end.

    python run.py                  # the pinned epoch (2026-09-29 20:00 UTC)
    python run.py --epoch now      # the Sun as it is right now
    python run.py --epoch 2026-10-05T12:00

Downloads whatever is not already cached, measures everything, draws the figures,
the time-lapse and the interactive page, and prints every number the README quotes.
Once the cache is complete this runs offline.
"""

from __future__ import annotations

import argparse
import sys
import time
import warnings

warnings.filterwarnings("ignore")

import analysis  # noqa: E402
import config as C  # noqa: E402
import fetch_data  # noqa: E402


def _t(iso: str) -> str:
    return iso[:16].replace("T", " ")


def print_summary(r: dict, m: dict) -> None:
    line = "=" * 78
    total = sum(v["bytes"] for v in m["files"].values())
    print(f"\n{line}\nTHE SUN ON {_t(r['epoch'])} UTC   (cache: {len(m['files'])} files, "
          f"{total / 1e6:.0f} MB)\n{line}")
    s = r["snapshot"]
    print(f"  Sun-SDO distance {s['dsun_au']:.4f} AU; apparent radius {s['rsun_arcsec']:.1f}\"; "
          f"B0 = {s['b0_deg']:+.2f} deg; Carrington rotation {s['carrington_rotation']}")

    x = r["xrays"]
    print(f"\nX-RAYS (GOES-{x['satellite']}), the week to the epoch")
    print(f"  {x['n_events']} flares ({x['list_source']}): "
          + ", ".join(f"{k} {v}" for k, v in x["counts_by_class"].items() if v)
          + f"; largest {x['largest']['class']} at {_t(x['largest']['peak'])}")
    print(f"  background between flares: {x['background_class']}; daily minimum "
          f"{' '.join(x['daily_minimum_class'])}")
    cc = x["cross_check_swpc_vs_ncei"]
    for ch, v in cc.items():
        print(f"  real-time vs science flux ({ch}): median ratio {v['median_ratio']:.3f} "
              f"(5-95%: {v['p05_ratio']:.2f}-{v['p95_ratio']:.2f}, {v['n_minutes']} min)")

    f = r["flare"]
    if f.get("selection"):
        sel, en, loc = f["selection"], f["energetics"], f["location"]
        reg = loc.get("noaa") or {}
        print(f"\nTHE FLARE OF THE WEEK: {sel['class']} (SWPC real-time: {sel['swpc_class']}), "
              f"peak {_t(sel['peak'])}")
        print(f"  site: ({loc['hpc_x_arcsec']:.0f}\", {loc['hpc_y_arcsec']:.0f}\"), "
              f"{loc['radius_rsun']:.3f} R_sun from centre -> NOAA {reg.get('region')} "
              f"on the limb (longitude {reg.get('stonyhurst_lon_at_flare', float('nan')):.0f} deg)")
        print(f"  temperature peak {en['t_max_MK']:.1f} MK at {_t(en['time_t_max'])[11:]} "
              f"({f['t_peak_before_flux_peak_min']:+.0f} min); emission measure peak "
              f"{en['em_cm3']:.2e} cm^-3 at {_t(en['time_em_max'])[11:]} "
              f"({f['em_peak_after_flux_peak_min']:+.0f} min)")
        print(f"  hot-loop area {f['hot_area']['area_cm2']:.2e} cm^2 -> volume "
              f"{en['volume_cm3']:.2e} cm^3, n_e = {en['density_cm3']:.2e} cm^-3, thermal energy "
              f"{en['thermal_energy_erg']:.2e} erg, radiative cooling time "
              f"{en['radiative_cooling_time_min']:.0f} min")
        print(f"  cooling sequence (episode ends {_t(f['episode_end'])[11:]}, renewed heating):")
        for row in f["cooling_sequence"]:
            print(f"    {row['wavelength']:>5} A {row['ion']:<18} log T {row['log_t']:.1f}  "
                  f"{'>' if row['lower_limit'] else ' '}{row['delay_min']:6.1f} min  ({row['kind']})")

    ro = r["rotation"]
    a, b, c = ro["fits"]["ABC"]["coef"]
    ea, eb, ec = ro["fits"]["ABC"]["err"]
    print(f"\nDIFFERENTIAL ROTATION ({ro['n_frames']} magnetograms, {ro['n_pairs']} pairs "
          f"{ro['pair_separation_h']} h apart)")
    print(f"  omega = {a:.3f}({ea:.3f}) {b:+.2f}({eb:.2f}) sin^2 {c:+.2f}({ec:.2f}) sin^4 deg/day "
          f"sidereal; chi2 {ro['fits']['ABC']['chi2']:.1f}/{ro['fits']['ABC']['dof']}; "
          f"A systematic +-{ro['A_systematic']:.3f}")
    alts = ", ".join(f"{k} h: {v['A']:.3f}" for k, v in ro["alternative_separations"].items())
    print(f"  other pair spacings: {alts}")
    print(f"  equator: {ro['sidereal_period_equator_d']:.2f} d sidereal, "
          f"{ro['synodic_period_equator_d']:.2f} d as seen from Earth; 30 deg: "
          f"{ro['sidereal_period_30_d']:.2f} d; 45 deg: {ro['sidereal_period_45_d']:.2f} d")
    for k, v in ro["comparisons"].items():
        print(f"  vs {v['label']:<44} {v['mean_diff']:+.3f} +- {v['mean_diff_err']:.3f} deg/day, "
              f"chi2 {v['chi2']:.0f}/{v['n']}")

    rg = r["regions"]
    print(f"\nACTIVE REGIONS: {len(rg['noaa_regions'])} NOAA regions tracked")
    for s_ in rg["noaa_regions"]:
        om = (f"omega {s_['omega']:.2f}+-{s_['omega_err']:.2f} (profile {s_['omega_profile']:.2f})"
              if "omega" in s_ else "too short for a rate")
        print(f"  {s_['noaa']}  lat {s_['lat']:+5.1f}  tracked {s_['span_h']:4.0f} h  "
              f"max {s_['max_area_msh']:5.0f} MSH  leading {'+' if s_.get('leading_polarity', 0) > 0 else '-'}"
              f"  {om}")
    print(f"  Hale's law (Cycle 25: N leads +, S leads -): {rg['hale_obeyed']} of "
          f"{rg['hale_tested']} bipolar regions; Joy's-law tilt sign: {rg['joy_positive_tilts']} of "
          f"{rg['hale_tested']}")

    lb = r["limb"]
    if "continuum" in lb:
        lc = lb["continuum"]
        print(f"\nLIMB DARKENING (HMI continuum {_t(lc['date'])})")
        print(f"  fit u = {lc['u_fit']:.3f}, v = {lc['v_fit']:.3f}  (Allen: {lc['u_allen']:.3f}, "
              f"{lc['v_allen']:.3f});  I(limb)/I(centre) {lc['limb_ratio_fit']:.3f} (Allen "
              f"{lc['limb_ratio_allen']:.3f}, grey 0.400)")
        print(f"  data vs Allen: {100 * lc['rms_resid_allen']:.1f}% rms, max "
              f"{100 * lc['max_abs_resid_allen']:.1f}%; vs grey atmosphere "
              f"{100 * lc['rms_resid_eddington']:.1f}% rms; scattered light off the limb "
              f"{100 * lc['scattered_light_fraction']:.2f}%")
        print(f"  temperature: {lc['t_at_tau_0p1']:,.0f} K at tau 0.1 -> "
              f"{lc['t_at_tau_1']:,.0f} K at tau 1 (anchored at T_eff at tau 2/3)")
    co = lb["corona"]
    print(f"\nTHE CORONA ABOVE THE LIMB ({len(co['quiet_sectors'])} quiet sectors of "
          f"{co['sector_width_deg']:.0f} deg; fit {co['fit_range_mm'][0]:.0f}-"
          f"{co['fit_range_mm'][1]:.0f} Mm)")
    for wl, v in co["channels"].items():
        print(f"  {wl:>4} A {v['ion']:<16} formation {v['t_formation_MK']:.2f} MK -> hydrostatic "
              f"{v['t_quiet_median_MK']:.2f} MK (16-84%: {v['t_quiet_p16_MK']:.2f}-"
              f"{v['t_quiet_p84_MK']:.2f}); plain exponential {v['t_plain_exponential_MK']:.2f} MK")
    pr = lb["profiles"]
    print("  limb brightening (peak / disc centre): "
          + ", ".join(f"{k.split('_')[0]} {v:.2f}" for k, v in pr.items() if k.endswith("_limb_peak")))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--epoch", default=None, help="'now' or an ISO time (UTC)")
    parser.add_argument("--no-figures", action="store_true", help="numbers only")
    parser.add_argument("--no-movie", action="store_true", help="skip the time-lapse (slowest step)")
    args = parser.parse_args(argv)

    t0 = time.time()
    paths = fetch_data.fetch_all(C.parse_epoch(args.epoch), quiet=True)
    products, results = analysis.analyse(paths, quiet=True)
    analysis.save(paths, products, results)
    if not args.no_figures:
        import figures
        import interactive

        written = figures.make_all(products, results, paths.figures)
        written.append(interactive.build(products, results, paths.figures))
        if not args.no_movie:
            import movie

            sel = results["flare"].get("selection")
            written += movie.render(paths, products["xrs"], sel["peak"] if sel else None,
                                    paths.figures)
        print(f"wrote {len(written)} files to {paths.figures}")
    print_summary(results, fetch_data.load_manifest(paths))
    print(f"\n(done in {time.time() - t0:.0f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
