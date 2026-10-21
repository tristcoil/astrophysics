"""Run the whole project end to end.

    python run.py

Fetches whatever is not cached, registers every plate on Gaia, measures the
Crab's filament motions and rewinds them, fits Cas A's expansion in fifteen
Chandra epochs, tries the same on Cas A's optical shell, analyses Green's
catalogue and the ATNF pulsar motions, writes every figure, the animations and
the interactive page, saves data/results.json, and prints a summary that grades
our numbers against the published ones. Once the cache is warm it runs offline
(30-45 minutes, most of it the Crab's patch tracking and processing
variants; later runs reuse the cached vectors and take a few minutes).

    python run.py --refresh       # recompute every cached analysis product
    python run.py --no-figures    # numbers only
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import warnings

warnings.filterwarnings("ignore")
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")

import numpy as np

import config
import literature as L


def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def _have_data() -> bool:
    needed = [config.GAIA_DIR / "crab_gaia_dr3.tsv", config.GAIA_DIR / "casa_gaia_dr3.tsv",
              config.CATALOG_DIR / "green_snrs_2024oct.html", config.ATNF_DB]
    needed += [config.PLATE_DIR / f"{p.key}.fits" for p in config.PLATES if not p.survey.startswith("ps1-")]
    n_chandra = len(list(config.CHANDRA_DIR.glob("*_cntr_img2.fits.gz")))
    n_obs = sum(len(v) for v in config.CHANDRA_EPOCHS.values())
    return all(p.exists() for p in needed) and n_chandra >= n_obs


def summary(r: dict) -> None:
    c, x, o, k = r["crab"], r["casa_xray"], r["casa_optical"], r["kicks"]
    m = c["main"]
    p = c["pulsar"]
    print("\n================ L3: Rewind a supernova remnant ================")
    print("\nThe Crab (POSS-I red 1951 -> Pan-STARRS r 2012):")
    print(f"  filament vectors used ............ {c['pairs'][m['pair']]['n_used']}")
    print(f"  convergence date (fit) ........... {m['T']:.1f} +- {m['T_err_total']:.1f} CE "
          f"(stat {m['T_err_bootstrap']:.1f}, processing {m['T_err_processing']:.1f})")
    print(f"     Martin et al. 2025 ............ {L.MARTIN2025['T_fit']} +- {L.MARTIN2025['T_fit_err']}")
    print("     Nugent 1998 ................... 1130 +- 16;  Wyckoff & Murray 1977: 1120 +- 7")
    for name, T in m["T_red_pairs"].items():
        print(f"     {name:24s} ......... {T:.1f}")
    print(f"  closest approach ................. {c['pairs'][m['pair']]['T_closest_approach']:.0f} "
          f"(Martin et al. 2025: {L.MARTIN2025['T_closest']})")
    print(f"  expansion rate ................... {m['expansion_rate_pct_yr']:.4f} %/yr; "
          f"median filament {m['velocity_median_kms']:.0f} km/s at 2 kpc")
    a = c["acceleration"]
    print(f"  today's speed / mean since 1054 .. {a['speed_now_over_mean']:.3f} "
          f"(constant push: {100 * a['fraction_gained_constant_push']:.0f}% of today's speed gained)")
    print(f"  centre ........................... ({m['centre'][0]:.1f}, {m['centre'][1]:.1f})\" from the pulsar; "
          f"Martin et al. ({c['martin2025_centre'][0]:.1f}, {c['martin2025_centre'][1]:.1f})")
    print(f"  pulsar (Gaia DR3) back to 1054 ... {p['offset_1054_from_centre_arcsec']:.1f}\" from our centre; "
          f"{p['v_transverse_kms']:.0f} km/s across the sky")

    fs = x["summary"]["regions"]["forward_shock"]["all"]
    ej = x["summary"]["regions"]["ejecta_ring"]["all"]
    print("\nCas A in X-rays (Chandra, 15 epochs 2000-2019):")
    print(f"  forward shock .................... {fs['rate_pct_yr']:.3f} +- {fs['rate_err']:.3f} %/yr, "
          f"m = {fs['m']:.2f}, {fs['v_kms']:.0f} km/s")
    print(f"     Vink et al. 2022 .............. {L.VINK2022_FS['mean_rate']} %/yr (sector mean), m = "
          f"{L.VINK2022_FS['mean_m']}, {L.VINK2022_FS['mean_v_kms']:.0f} km/s")
    sec = x["fs_sector_comparison"]
    print(f"     sectors within 0.03 %/yr of Vink: {sec['n_within_0.03']} of {sec['n']}")
    print(f"  ejecta ring ...................... {ej['rate_pct_yr']:.3f} %/yr; sector pattern vs Vink's "
          f"reverse shock: r = {x['rs_sector_correlation']:.2f}")
    off = x["neutron_star"]["offset"]
    print(f"  neutron star offset .............. {off['arcsec']:.1f}\" at PA {off['pa_deg']:.0f} deg -> "
          f"{off['implied_v_kms']:.0f} km/s (Fesen et al. 2006: 7\" at 169 deg, ~350 km/s)")

    print("\nCas A in visible light (plates 1954, 1991; Pan-STARRS ~2012):")
    for name, pr in o["pairs"].items():
        print(f"  {name:20s} best stretch {pr['s']:.3f} (coasting shell: {pr['s_expected']:.3f})")
    cc = o["crab_control"]
    print(f"  Crab control: best stretch {cc['s_best']:.3f} (expected {cc['s_expected']:.3f})")

    st = k["stats"]
    print("\nPulsar kicks (ATNF v2.8.1):")
    print(f"  young pulsars (< 3 Myr): mean 2-D speed {st['young']['v_pec_mean']:.0f} +- "
          f"{st['young']['v_pec_mean_err']:.0f} km/s (n = {st['young']['n']}); "
          f"Hobbs et al. 2005: {L.HOBBS2005['v2d_young']:.0f} +- {L.HOBBS2005['v2d_young_err']:.0f}")
    print(f"  normal pulsars: {st['normal']['v_pec_mean']:.0f} +- {st['normal']['v_pec_mean_err']:.0f} km/s "
          f"(n = {st['normal']['n']}); Hobbs et al.: {L.HOBBS2005['v2d_normal']:.0f} +- {L.HOBBS2005['v2d_normal_err']:.0f}")
    print(f"  moving away from their remnant's centre (within 35 deg): {k['n_outward']} of {k['n_secure']} "
          f"secure associations (chance: p = {k['p_outward']:.1e})")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--refresh", action="store_true", help="recompute cached analysis products")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args(argv)

    if not _have_data():
        import fetch_data

        if fetch_data.main([]) != 0:
            return 1

    import casa_optical
    import casa_xray
    import crab
    import kicks
    import motions as M
    import remnants
    from scipy import stats

    print("The Crab ...")
    c = crab.analyse(refresh=args.refresh)
    print("Cas A, Chandra ...")
    xfits = casa_xray.run(refresh=args.refresh)
    xsum = casa_xray.summarise(xfits)
    ns = casa_xray.neutron_star(xfits)
    fs = xsum["regions"]["forward_shock"]
    diffs = [fs[str(pa)]["rate_pct_yr"] - v for pa, v in zip(L.VINK2022_FS["pa"], L.VINK2022_FS["rate"])]
    rs = xsum["regions"]["ejecta_ring"]
    rs_corr = float(np.corrcoef([rs[str(pa)]["rate_pct_yr"] for pa in L.VINK2022_RS["pa"]], L.VINK2022_RS["rate"])[0, 1])
    x = {"summary": xsum, "neutron_star": ns,
         "fs_sector_comparison": {"n": len(diffs), "n_within_0.03": int(np.sum(np.abs(diffs) < 0.03)),
                                  "median_diff": float(np.median(diffs)),
                                  "mean_ours": float(np.mean([fs[str(pa)]["rate_pct_yr"] for pa in L.VINK2022_FS["pa"]]))},
         "rs_sector_correlation": rs_corr}
    print("Cas A, optical ...")
    o = casa_optical.analyse(refresh=args.refresh)
    print("Green's catalogue and pulsar kicks ...")
    g = remnants.analyse()
    k = kicks.analyse()
    import pandas as pd

    homes = pd.DataFrame(k["homes"])
    secure = homes[~homes.doubtful & (homes.offset_arcmin > 2.0)]
    n_out = int((secure.angle < 35).sum())
    k["n_secure"], k["n_outward"] = int(len(secure)), n_out
    k["p_outward"] = float(stats.binom.sf(n_out - 1, len(secure), 35.0 / 180.0))

    results = {"crab": c, "casa_xray": x, "casa_optical": o, "remnants": g, "kicks": k}
    config.RESULTS_PATH.write_text(json.dumps(_jsonable(results), indent=1))

    if not args.no_figures:
        print("Figures ...")
        import figures as F
        import interactive

        cr = crab.Crab()
        v = cr.vectors(*crab.PAIRS[crab.MAIN])
        fit = M.fit_expansion(v)
        vi = v.subset(fit.inliers)
        frame, epochs = casa_xray.load()
        green = remnants.load_green()
        hist = remnants.historical(green)
        age_casa = xsum["age_at_ref"]
        m_points = [
            ("Cas A forward shock (X-ray)", age_casa, fs["all"]["m"], fs["all"]["m_err"], True),
            ("Cas A ejecta ring (X-ray)", age_casa, rs["all"]["m"], rs["all"]["m_err"], True),
            ("Vink et al. 2022 (forward shock)", age_casa, L.VINK2022_FS["mean_m"], 0.096, False),
            ("Crab filaments (optical)", c["acceleration"]["true_age_yr"], c["acceleration"]["speed_now_over_mean"],
             c["main"]["T_err_total"] / c["acceleration"]["apparent_age_yr"], True),
        ]
        F.fig_crab_epochs(cr, c)
        F.fig_crab_vectors(cr, c, vi)
        F.fig_crab_rewind(c, vi)
        F.fig_crab_stills(cr, c, vi)
        F.fig_casa_xray(frame, epochs, xsum)
        F.fig_casa_expansion(xsum)
        F.fig_casa_neutron_star(frame, epochs, ns)
        F.fig_casa_optical(o)
        F.fig_galaxy_remnants(green, hist, m_points)
        F.fig_kicks(k, kicks.load_pulsars(), homes)
        F.gif_crab(cr)
        F.gif_casa(frame, epochs)
        interactive.build(cr, c, vi, frame, epochs, xsum)
    summary(results)
    return 0


if __name__ == "__main__":
    sys.exit(main())
