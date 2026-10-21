"""Regenerate every result and figure of project N3 from the cached data.

    python fetch_data.py     # once: ~3.4 GB streamed, ~0.45 GB kept
    python run.py            # ~1-1.5 h the first time (null tests), minutes afterwards
    python run.py --force    # recompute the cached fits too

Writes data/results.json and figures/.
"""

from __future__ import annotations

import json
import sys
import time

import pandas as pd

import config as C


def main(force: bool = False) -> dict:
    t0 = time.time()
    import prepare
    if force or not all(prepare.cache_path(p).exists() for p in C.KEEP_SATS):
        print("== preparing satellite-day arrays")
        prepare.build()

    print("== GPS clock budget")
    import gps_budget
    path = C.CACHE_DIR / "gps_satellite_days.parquet"
    if force or not path.exists():
        gps_budget.satellite_days().to_parquet(path)
    gps = pd.read_parquet(path)
    gps_summary = gps_budget.summary(gps)

    print("== the Sun's redshift")
    import solar
    sun = solar.analyse(force=force)

    print("== Galileo redshift test")
    import redshift_test
    gal = redshift_test.run(force=force)
    for prn in ("E14", "E18"):
        redshift_test.once_per_rev(prn, force=force)
        redshift_test.folded(prn, force=force)
        redshift_test.clock_history(prn, force=force)

    print("== Sirius B and the ladder")
    import sirius
    import ladder
    wd = sirius.analyse()
    rungs = ladder.table(gal, sun, wd, gps_summary)

    print("== figures")
    import figures
    figures.gps_budget(gps, gps_summary)
    figures.solar(pd.read_parquet(C.CACHE_DIR / "solar_lines.parquet"),
                  pd.read_csv(C.CACHE_DIR / "solar_windows.csv"), sun)
    figures.galileo_clocks()
    figures.galileo_fold(gal)
    figures.galileo_years(gal)
    figures.galileo_sun_vs_perigee()
    figures.galileo_history()
    figures.galileo_budget_and_j2(gal)
    ladder.figure(rungs)
    print("== animation")
    import animation
    if force or not (C.FIGURE_DIR / "galileo_clock.mp4").exists():
        animation.make()

    results = {"gps": gps_summary, "sun": {k: v for k, v in sun.items() if k != "binned"},
               "galileo": {"combined": gal["combined"], "j2_tide": gal["j2_tide"],
                           **{prn: {k: v for k, v in gal["per_satellite"][prn].items() if k != "yearly"}
                              for prn in ("E14", "E18")}},
               "sirius_b": wd, "ladder": rungs.to_dict("records"),
               "runtime_s": time.time() - t0}
    C.RESULTS_PATH.write_text(json.dumps(results, indent=1, default=float))
    print(f"done in {time.time() - t0:.0f} s -> {C.RESULTS_PATH}")
    return results


if __name__ == "__main__":
    main(force="--force" in sys.argv)
