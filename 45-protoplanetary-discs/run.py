"""Run the whole project end to end.

    python fetch_data.py     # once: ~2.5 GB (0.8 GB of continuum images, 1.7 GB CO cube)
    python run.py            # every measurement, number and figure, plus the explorer page
    python test_physics.py   # self-checks of the methods on synthetic inputs

The first run takes ~5-10 minutes on a few cores (the 19 disc pipelines, then the CO
fit's MCMC); results are cached in data/cache/, so later runs take about a minute.
``--refit`` starts everything over. Every headline number is printed next to the
published value it is graded against.
"""

from __future__ import annotations

import argparse
import json
import os
import warnings

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402

import analysis  # noqa: E402
from config import RESULTS_PATH  # noqa: E402


def _rule(title: str) -> None:
    print("\n" + "=" * 86)
    print(title)
    print("=" * 86)


def print_summary(r: dict) -> None:
    _rule("THE DATA: 20 DSHARP discs (ALMA 1.25 mm) and HL Tau (ALMA 1.0 mm)")
    im = r["images"]
    ratio = np.array([v["flux_mjy"] / v["flux_pub"] for v in im.values()])
    print(f"  total fluxes / DSHARP I Table 4: median {np.median(ratio):.3f} (range {ratio.min():.2f}-{ratio.max():.2f}); "
          f"noise {np.median([v['rms_ujy'] for v in im.values()]):.0f} uJy/beam (DSHARP I: "
          f"{np.median([v['rms_pub'] for v in im.values()]):.0f})")
    for k, v in im.items():
        if abs(v["flux_mjy"] / v["flux_pub"] - 1) > 0.04:
            why = ("the released image differs from the tabulated one" if abs(v["peak_mjy"] / v["peak_pub"] - 1) > 0.03
                   else "same image; the total depends on how much faint outer emission the aperture includes")
            print(f"    {k}: {v['flux_mjy']:.0f} vs {v['flux_pub']:.0f} mJy, peak {v['peak_mjy']:.2f} vs {v['peak_pub']:.2f} "
                  f"mJy/beam ({why})")

    _rule("DISTANCES (Gaia DR3, zero-point corrected; DSHARP used DR2)")
    for k, s in r["distances"].items():
        if s.get("distance_source", "").startswith("Gaia"):
            continue
        print(f"  {k}: {s['distance']:.1f} +/- {s['distance_err']:.1f} pc -- {s['distance_source']}")
    dd = [s["distance"] for k, s in r["distances"].items() if "parallax" in s]
    print(f"  {len(dd)} hosts from their own parallaxes: {min(dd):.0f}-{max(dd):.0f} pc")

    _rule("ORIENTATIONS: our own fits vs DSHARP II (and ALMA 2015 for HL Tau)")
    rows = r["geometry"]
    for k, g in rows.items():
        o, p = g["ours"], g["published"]
        print(f"  {g['name']:10s} {o['method']} i = {o['inc']:5.2f} +/- {o['inc_err']:4.2f} ({p['method']} {p['inc']:5.2f} "
              f"+/- {p['inc_err']}), PA = {o['pa']:6.2f} +/- {o['pa_err']:5.2f} ({p['pa']:6.2f}); "
              f"R_dust {g['r_dust_au_dsharp']:5.1f} au ({g['r_dust_pub']})")
    nm = sum(1 for g in rows.values() if (g["ours"]["inc_err"] or 99) < 10)
    n2m = sum(1 for g in rows.values() if (g["ours"]["inc_err"] or 99) < 10 and abs(g["d_inc_sigma"] or 0) <= 2)
    print(f"  -> {n2m} of the {nm} discs whose tilt we can measure (inclination error < 10 deg) agree within 2 sigma; "
          f"{len(rows) - nm} unconstrained")

    _rule("GAPS AND RINGS vs DSHARP II")
    fe = r["features"]
    for m, v in fe["by_method"].items():
        lab = {"E": "ellipse-fitted", "R": "profile extrema", "V": "found by eye (plateaus)", "-": "HL Tau (ALMA 2015)"}[m]
        md = f", median |offset| {v['median_abs_dr_mas']:.1f} mas = {v['median_abs_dr_beam']:.2f} beam" if v["found"] else ""
        print(f"  {lab:24s}: {v['found']:2d} of {v['n']:2d} recovered{md}")
    w = [x for x in r["widths"] if x.get("gap_width_pub") is not None and x["gap_width_au"] is not None
         and not x.get("gap_width_pub_upper")]
    dw = np.array([x["gap_width_au"] - x["gap_width_pub"] for x in w])
    dp = np.array([x["depth"] / x["depth_pub"] for x in r["widths"] if x.get("depth_pub")])
    print(f"  gap widths: {len(w)} compared, median |difference| {np.median(np.abs(dw)):.2f} au; "
          f"depths: median ratio {np.median(dp):.3f}")
    secure = [x for x in fe["extras"] if x["secure"]]
    print(f"  features we find that DSHARP II does not list: {len(fe['extras'])} ({len(secure)} above 5 sigma: "
          + ", ".join(f"{x['disc']}:{x['label']}" for x in secure) + ")")

    _rule("PLANETS THAT WOULD CARVE THE GAPS (DSHARP VII's relation, our widths)")
    pl = r["planets"]
    for row in pl["rows"]:
        if not row.get("found"):
            print(f"  {row['disc']:9s} {row['r_gap_pub']:4.0f} au: gap not found in our profile")
            continue
        m = row["masses"]["DSD1_0.001"]
        lo = min(v["m"] for kk, v in row["masses"].items() if kk.startswith("DSD1"))
        hi = max(v["m"] for kk, v in row["masses"].items() if kk.startswith("DSD1"))
        print(f"  {row['name']:9s} {row['r_gap_pub']:4.0f} au: Delta = {row['delta']:.3f} +/- {row['delta_err']:.3f} "
              f"(DSHARP VII {row['delta_pub']:.2f}) -> M_p = {m['m']:6.3f} M_J (DSHARP VII {m['pub']:5.2f}); "
              f"alpha 1e-4..1e-2: {lo:.3f}-{hi:.3f}")
    print(f"  -> masses: median ratio ours/DSHARP VII {pl['median_ratio']:.3f}, rms {pl['log_rms']:.3f} dex")

    _rule("SNOW LINES: do gaps sit at ice-line temperatures?")
    sn = r["snowlines"]
    for k, v in sn["fraction_in_range"].items():
        print(f"  {k:18s}: {v * 100:4.0f}% of {sn['n_gaps']} gaps (chance, uniform in log T: "
              f"{sn['uniform_log_t_expectation'][k] * 100:.0f}%)")

    _rule("SPIRAL ARMS: pitch angles vs DSHARP III")
    for k, s in r["spirals"].items():
        arms = ", ".join(f"{a}: {v['pitch']:.1f} +/- {v['pitch_err']:.1f}" for a, v in sorted(s["arms"].items()))
        pub = ", ".join(f"{a}: {v['pitch']:.1f}" for a, v in s["published"].items())
        print(f"  {k:8s} Fourier {s['fourier']['pitch']:.1f} deg | traced {arms} | DSHARP III {pub}")

    k = r["kinematics"]
    if "fit" in k:
        f, t = k["fit"], k["published"]
        _rule("HD 163296 WEIGHED BY ITS DISC (CO J=2-1) vs Teague et al. 2019")
        print(f"  near side PA {k['near_side_pa']:.0f} deg (NE), gas circulates {k['rotation_on_sky']} on the sky")
        print(f"  M* = {f['mstar']:.3f} +/- {f['errors']['mstar']:.3f} (stat) +/- {k['mstar_systematic']:.3f} (sys) Msun "
              f"at 101 pc  (Teague: {t['mstar'][0]} +/- {t['mstar'][1]});  {k['mstar_at_gaia']['mstar']:.3f} at Gaia's "
              f"{k['mstar_at_gaia']['distance']:.1f} pc  (DSHARP I tracks: 2.04)")
        print(f"  i = {f['inc']:.2f} ({t['inc_deg'][0]}), PA = {f['pa']:.2f} ({t['pa_deg'][0]}), v_sys = {f['vsys']:.4f} "
              f"({t['vlsr_kms'][0]}) km/s")
        print(f"  surface: z0 = {f['z0_mas']:.0f} +/- {f['errors']['z0']:.0f} mas ({t['z0_mas'][0]:.0f}), phi = "
              f"{f['phi']:.2f} ({t['phi'][0]}), z1 = {f['z1_mas']:.1f} ({t['z1_mas'][0]:.0f}), psi = {f['psi']:.2f} "
              f"({t['psi'][0]})")
        print(f"  residuals: median |dv| {k['resid_mad_ms']:.0f} m/s; jitter {f['jitter_ms']:.0f} m/s")
        print(f"  gas pressure minima from the rotation curve: {', '.join(f'{x:.0f}' for x in k['pressure_minima_au_101pc'])} au "
              "at 101 pc (Teague et al. 2018 planets: 83 and 137 au)")
        for n, kk in k["kinks"].items():
            print(f"  {n}: velocity residual at the kink {kk['resid_here'] * 1e3:+.0f} m/s vs {kk['annulus_mean'] * 1e3:+.0f} "
                  f"+/- {kk['annulus_std'] * 1e3:.0f} along its orbit (z = {kk['z_score']:+.1f})")

    _rule("DUST MASSES (optically thin, T = 20 K; Earth masses)")
    dm = r["dust_masses"]
    ms = sorted(((v["m_dust_earth"], k) for k, v in dm.items()), reverse=True)
    print("  " + ", ".join(f"{k} {m:.0f}" for m, k in ms))
    print("  (kappa = 2.3 cm^2/g at 230 GHz; with DSHARP's 0.43 cm^2/g opacity the masses are 5.4x higher)")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--refit", action="store_true", help="redo the disc pipelines and the CO fit")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args(argv)
    if args.refit:
        from config import CACHE_DIR

        (CACHE_DIR / "kinematics" / "hd163296.json").unlink(missing_ok=True)
    res = analysis.run(refit=args.refit)
    if not args.no_figures:
        import figures
        import interactive

        r = json.loads(RESULTS_PATH.read_text())
        for fn in (figures.gallery_same_scale, figures.gallery_annotated, figures.profiles_grid,
                   figures.geometry_comparison, figures.features_comparison, figures.planets_mass_separation,
                   figures.planet_mass_check, figures.snowlines, figures.hd163296_kinematics, figures.hd163296_channels,
                   figures.spirals_figure, figures.hltau_figure, figures.gap_to_planet):
            fn(r)
        interactive.build(r)
    print_summary(json.loads(RESULTS_PATH.read_text()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
