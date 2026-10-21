"""Build the whole project with one command.

    python run.py              # fetch what is missing, every analysis step, every figure, the
                               # interactive page, a summary
    python run.py --check      # also run the self-checks in test_physics.py
    python run.py --cassini    # also fetch Cassini's raw 2002 Doppler data (~1.3 GB) and refit it

The 1919 and 1922 measurements come from tables transcribed by hand from the
scanned papers (``transcribed/``, committed); nothing in them is downloaded.
Once ``data/`` is filled this runs offline in about 20 seconds; --cassini adds
about 4 minutes once its one-second Doppler model is cached.
"""

from __future__ import annotations

import argparse
import sys
import warnings


def summary(r: dict) -> None:
    print("\n" + "=" * 78)
    print("N2 — Bending starlight: from the 1919 eclipse to Cassini")
    print("=" * 78)
    fr = r["frame"]
    print(f"Table I vs Gaia DR3 (1919): rms {fr['rms_arcsec']:.2f}\", scale 1 {fr['scale'] - 1:+.1e} "
          f"± {fr['scale_err']:.1e}, rotation {fr['rotation_deg']:+.3f} ± {fr['rotation_err_deg']:.3f} deg")
    for site, s in fr["sites"].items():
        print(f"  {site:9s} 1919 predictions vs exact GR: max {s['pred_vs_table1_max_arcsec']:.3f}\"; "
              f"Sun {s['sun_offset_arcsec']:.0f}\" from where Table I implies")
    c = r["contacts"]
    print(f"Totality (JPL): Sobral {c['sobral']['second_contact_utc']}-{c['sobral']['third_contact_utc']} UT "
          f"(recorded 2nd contact {c['sobral']['recorded']['second_contact_utc']}); Principe "
          f"{c['principe']['second_contact_utc']}-{c['principe']['third_contact_utc']} (computed 1919: "
          f"{c['principe']['recorded']['second_contact_utc']}-{c['principe']['recorded']['third_contact_utc']})")
    p = r["plate1"]
    print(f"Plate 1: 7 stars identified to {p['rms_px'] / p['scale_px_per_arcmin'] * 60:.1f}\" rms; "
          f"exposure dated {p['dated_utc']} UT ± {p['sigma_t_s']:.0f} s"
          + (f" = 4-inch plate {p['exposure_match']['plate']}" if (p.get("exposure_match") or {}).get("plate") else ""))
    s4 = r["sobral_4inch"]
    rep = s4["replication"]
    print(f"Sobral 4-inch, 1920 reduction replicated: {rep['limb']:.2f}\" (RA {rep['limb_ra']:.2f}, Dec {rep['limb_dec']:.2f}); "
          f"published 1.98 ± 0.12 (p.e.)")
    for k, v in s4["joint_fit"].items():
        print(f"  joint least squares, {k:4s}: {v['limb']:.2f} ± {v['limb_err']:.2f}\" (1σ)")
    a = r["astrographic"]
    fa = a["fits"]["scale_common_free"]
    print(f"Sobral astrographic: plates' scale {a['replication']['limb_plates_scale']:.2f}\", refraction-only scale "
          f"{a['replication']['limb_refraction_scale']:.2f}\" (published 0.93 and 1.52); "
          f"common scale fitted {fa['limb']:.2f} ± {fa['limb_err']:.2f}\""
          + (f" (± {fa['limb_err_jackknife']:.2f} from a jackknife over plates)" if "limb_err_jackknife" in fa else ""))
    if "principe" in r:
        pr = r["principe"]
        print(f"Principe: replication {pr['replication']['limb']:.2f}\" (published 1.61 ± 0.30 p.e.); "
              f"our fit {pr['fit']['limb']:.2f} ± {pr['fit']['limb_err']:.2f}\"")
    cb = r.get("combined_1919") or {}
    if cb.get("sobral_two", {}).get("sigma_apart") is not None:
        al = cb["all"]
        print(f"1919 combined: the two Sobral instruments {cb['sobral_two']['sigma_apart']:.1f}σ apart; all three "
              f"{al['limb']:.2f} ± {al['limb_err_scaled']:.2f}\" after scaling (χ² = {al['chi2']:.1f} for {al['dof']} dof)")
    if "wallal" in r:
        w = r["wallal"]
        print(f"Wallal 1922 (five-foot cameras): our fit {w['fit']['limb']:.2f} ± {w['fit']['limb_err']:.2f}\" "
              f"(published 1.82 ± 0.15 p.e.)")
    cas = (r.get("cassini") or {}).get("summary") or {}
    if cas:
        print(f"Cassini 2002 raw Doppler: gamma - 1 = {cas['gamma_minus_1']:+.1e} ± {cas['gamma_err_stat']:.1e} (stat) "
              f"± {cas['gamma_err_syst']:.1e} (syst); total ± {cas['gamma_err_total']:.1e} "
              f"(Bertotti et al. 2003: (2.1 ± 2.3)e-5)")
    sh = r["shapiro"]
    print(f"Shapiro delay to Cassini: {sh['delay_2way_max_us']:.0f} µs round trip at {sh['b_min_rsun']:.2f} R_sun "
          f"({sh['closest_date']})")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true", help="also run test_physics.py")
    ap.add_argument("--cassini", action="store_true", help="fetch and refit Cassini's raw Doppler data")
    args = ap.parse_args(argv)
    warnings.filterwarnings("ignore")

    import analysis
    import fetch_data
    import figures

    fetch_data.fetch_all(cassini=args.cassini)
    if args.cassini:
        import cassini_doppler

        cassini_doppler.run_cassini()
    print("Analysis")
    results = analysis.run()
    print("Figures")
    figures.main([])
    import explorer

    print("Interactive page")
    explorer.main()
    summary(results)
    if args.check:
        import test_physics

        return test_physics.main()
    return 0


if __name__ == "__main__":
    sys.exit(main())
