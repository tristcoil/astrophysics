"""Run the whole project end to end.

    python run.py                # fetch what is missing, analyse, draw everything
    python run.py --no-movie     # skip the year movie (the slowest step)
    python run.py --offline      # use the cache only (no network)

The first run downloads ~2.7 GB (hours; the 55,000 thumbnails dominate) and
measures every image (~10 min on 6 cores). Later runs read the cached
photometry and take a few minutes; the movie adds ~15 min.
"""

from __future__ import annotations

import argparse
import time
import warnings

import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

import analysis
import epicdata
import photometry as ph
from config import FIGURE_DIR, RANDOM_SEED, RESULTS_PATH


def _rule(title: str) -> None:
    print(f"\n{title}\n" + "-" * len(title))


def print_summary(r: dict) -> None:
    print("=" * 78)
    print("THE EARTH AS AN EXOPLANET — DSCOVR/EPIC, 2015–2026")
    print("=" * 78)
    g = r["geometry"]
    _rule(f"Geometry: {g['n_images']:,} images on {g['n_days']:,} days")
    print(f"  phase angle Sun-Earth-DSCOVR {g['phase_deg']['min']:.1f}-{g['phase_deg']['max']:.1f} deg; "
          f"distance {g['dist_km']['min'] / 1e6:.2f}-{g['dist_km']['max'] / 1e6:.2f} million km; "
          f"latitude facing DSCOVR {g['subobs_lat']['min']:+.1f} to {g['subobs_lat']['max']:+.1f} deg")
    tr = r["transform"]
    _rule("Display images vs light (gamma)")
    print(f"  ocean-cloud mixing {tr['mixing_blue_vs_red']['gamma']:.2f} (blue), "
          f"{tr['mixing_green_vs_red']['gamma']:.2f} (green); Yang+2018 scatter "
          f"{tr['gamma_from_yang_std']}; adopted {tr['adopted']['gamma']} ({tr['adopted']['range']})")
    _rule("Length of the day (one year of 1080-pixel images)")
    for name, rot in r["rotation"].items():
        for ch in ("red", "green", "blue"):
            x = rot[ch]
            print(f"  {name} {ch:5s}: {x['period_h']:.4f} +- {x['period_err_h']:.4f} h "
                  f"({x['offset_s']:+.0f} +- {x['offset_err_s']:.0f} s); ACF peak {x['acf_peak_h']:.1f} h")
    if "long" in r:
        lg = r["long"]
        _rule(f"Eleven years of thumbnails ({lg['n_images']:,} images)")
        for ch in ("red", "blue"):
            x = lg[ch]
            print(f"  {ch}: day = 24 h {x['offset_s']:+.1f} +- {x['offset_err_s']:.1f} s "
                  f"(sidereal would be {(x['sidereal_period_h'] - 24) * 3600:+.0f} s); "
                  f"power at sidereal {x['power_at_sidereal']:.3f} vs peak {x['power_peak']:.3f}")
        p = r["phase"]["F_red"]
        print(f"  phase angle: {100 * p['slope_per_deg']:+.2f} +- {100 * p['slope_err']:.2f} %/deg (display); "
              f"Lambert {100 * r['phase']['lambert_slope_per_deg_at_median']:+.2f} %/deg")
    s = r["sampling"]
    _rule("How rarely can we look? (random week-long stretches, 680 nm)")
    import pandas as pd
    for key in ("2015-16 random", "2015-16 regular"):
        d = pd.DataFrame(s[key])
        print(f"  {key}: " + ", ".join(f"{row.interval_h:.0f} h: {row.p_detect:.2f}" for row in d.iloc[[0, 2, 4, 6, 9, 13]].itertuples()))
    _rule("Continents and clouds")
    for name, sur in r["surface"].items():
        pa = sur["pacific_vs_africa"]
        print(f"  {name}: Pacific/Africa red {pa['F_red']['pacific_over_africa']:.3f}, blue "
              f"{pa['F_blue']['pacific_over_africa']:.3f}, blue/red {pa['blue_red']['pacific_over_africa']:.3f}; "
              f"colour-land r = {sur['corr_with_land']['blue_red']:+.2f}; PC2-land r = {sur['pca']['corr_pc2_land']:+.2f}; "
              f"rotation share red {sur['rotation_vs_weather']['F_red']['rotation_share']:.0%}, "
              f"colour {sur['rotation_vs_weather']['blue_red']['rotation_share']:.0%}")
    for name, se in r["seasons"].items():
        print(f"  seasons {name}: brightest month {se['brightest_month']}, faintest {se['faintest_month']}, "
              f"range {se['seasonal_range_red']:.1%}")
    mo = r["moon"]
    _rule("The Moon")
    ev = mo["events"]
    print(f"  predicted 2015-2026: {ev['n_transits']} transits, {ev['n_occultations']} occultations, "
          f"{ev['n_near_misses']} near misses")
    im = mo.get("imaged", {})
    if im:
        print(f"  imaged by EPIC (login-gated L1A misc): {im['n_with_misc_granules']} of {im['n_outside_gaps']} "
              f"events outside archive gaps; random windows {im['n_controls_with_granules']} of "
              f"{im['n_control_windows']}")
    for key, t in mo["transits"].items():
        print(f"  {key}: display Moon/Earth red {t['display_red']['ratio']:.2f}, blue {t['display_blue']['ratio']:.2f}; "
              f"linear (gamma {analysis.GAMMA_NOMINAL}) red {t['moon_over_earth_red'][f'{analysis.GAMMA_NOMINAL:.1f}']:.2f}, "
              f"blue {t['moon_over_earth_blue'][f'{analysis.GAMMA_NOMINAL:.1f}']:.2f}; calibrated EPIC red "
              f"{mo['expected']['red']['ratio']:.2f}, blue {mo['expected']['blue']['ratio']:.2f}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--no-movie", action="store_true", help="skip the movie and animations")
    parser.add_argument("--no-figures", action="store_true", help="skip all figures")
    parser.add_argument("--offline", action="store_true", help="do not fetch anything")
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args(argv)

    t0 = time.time()
    if not args.offline:
        from fetch_data import ensure_data
        ensure_data(quiet=False)
    meta = epicdata.load_metadata()
    for kind in ("jpg", "thumbs"):
        imgs = [im for im in meta["image"] if epicdata.image_path(kind, im).exists()]
        ph.photometry_table(kind, imgs, workers=args.workers)
    rng = np.random.default_rng(RANDOM_SEED)
    results, art = analysis.run_all(rng)
    analysis.save_results(results)
    if not args.no_figures:
        import figures
        written = figures.make_all(results, art)
        if not args.no_movie:
            import movie
            written += movie.make_all(results, art)
        print(f"\nwrote {len(written)} figure files to {FIGURE_DIR}")
    print_summary(results)
    print(f"\nNumbers saved to {RESULTS_PATH}  ({time.time() - t0:.0f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
