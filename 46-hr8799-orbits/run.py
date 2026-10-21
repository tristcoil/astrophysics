"""Run the whole project end to end.

    python run.py                  # everything (resumes from data/cache)
    python run.py --stage fits     # one stage: data, fits, orbitize, stability,
                                   #   betapic, prediction, hgca, figures, animation

Every expensive step is cached in data/cache and resumable, so a second run
only redraws figures and rewrites data/results.json (offline). From scratch
the whole thing takes roughly a day of CPU on 12 cores:

* 15 nested-sampling fits of HR 8799 plus one of beta Pic b: 5-25 min each on
  one core, run side by side (python fits.py KEY);
* N-body stability: ~0.3 CPU-s per 40-Myr draw (most systems fall apart in
  ~20,000 yr), ~570,000 draws; ~90 s for each stable system, which is then
  replayed three times (resonant angles, a nudged twin, 100 Myr forwards);
* orbitize!: OFTI for b ~1 h (1-mas GPI points make OFTI slow), MCMC for d, e
  and beta Pic b a few hours each.
"""

from __future__ import annotations

import argparse
import json
import time
import warnings

import numpy as np

warnings.filterwarnings("ignore")

import config as C

#: N-body campaigns: fit key -> (draws, mass prior, tag), in the order they run
STABILITY_RUNS = [
    ("W18/near1248_coplanar", 200_000, "wang", ""),
    ("W18/near1248_coplanar", 150_000, "hotstart", "_hotstart"),
    ("ALL/near1248_coplanar", 20_000, "wang", ""),
    ("W18/near1248", 30_000, "wang", ""),
    ("W18/coplanar", 30_000, "wang", ""),
    ("W18/unconstrained", 30_000, "wang", ""),
    ("ALL/coplanar", 10_000, "wang", ""),
]
#: the fit whose stable systems get the resonance / twin / future passes
MAIN_STABLE = ("W18/near1248_coplanar", "")
#: the orbits drawn in the animation and the sky figures (all data, 1998-2025)
FINAL_FIT = "ALL/near1248_coplanar"


def stage_data():
    import astrometry as A
    import fetch_data
    fetch_data.ensure_data()
    df = A.build()
    df.to_csv(C.DERIVED / "astrometry.csv", index=False)


def stage_fits(nproc_parallel: int = 6):
    """Launch every missing nested-sampling fit as its own process."""
    import subprocess
    import sys
    import fits
    todo = [k for k in list(fits.FITS) + ["W18/coplanar#2", "W18/near1248_coplanar#2"]
            if not fits.cache_path(k).exists()]
    running = []
    while todo or running:
        while todo and len(running) < nproc_parallel:
            k = todo.pop(0)
            log = open(C.CACHE / "logs" / f"{k.replace('/', '_').replace('#', '_')}.log", "w")
            running.append(subprocess.Popen([sys.executable, "-u", "fits.py", k], stdout=log, stderr=log,
                                            cwd=C.HERE))
        running = [p for p in running if p.poll() is None]
        time.sleep(5)


def stage_orbitize():
    """orbitize! cross-checks: OFTI for every planet, MCMC where OFTI fails (d, e)
    and for beta Pic b. Each runs in its own process."""
    import subprocess
    import sys
    jobs = [["ofti", p] for p in C.PLANETS] + [["mcmc", "d", "2"], ["mcmc", "e", "2"], ["betapic", "2"]]
    procs = []
    for j in jobs:
        log = open(C.CACHE / "logs" / f"orbitize_{'_'.join(j[:2])}.log", "a")
        procs.append(subprocess.Popen([sys.executable, "-u", "orbitize_fits.py", *j], stdout=log, stderr=log,
                                      cwd=C.HERE))
    for p in procs:
        p.wait()


def stage_stability(nproc: int = C.N_PROCESSES):
    """N-body campaigns, the most informative first (every step resumes from its cache)."""
    import stability as St
    key, tag = MAIN_STABLE

    def campaign(i):
        k, n, prior, t = STABILITY_RUNS[i]
        St.run(k, n, nproc=nproc, tag=t, mass_prior=prior, seed_offset=0 if prior == "wang" else 1)

    campaign(0)                                     # Wang's set-up
    campaign(1)                                     # the same with hot-start masses
    St.resonances(key, nproc=nproc, tag=tag)        # which angles librate
    campaign(2)                                     # the orbits fitted to all data
    St.twins(key, nproc=nproc, tag=tag)             # chaos: nudged copies
    St.future(key, years=100e6, nproc=nproc, tag=tag)
    for i in range(3, len(STABILITY_RUNS)):         # the other models
        campaign(i)
    St.resonances(key, nproc=nproc, tag="_hotstart")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all")
    a = ap.parse_args()
    stages = {"data": stage_data, "fits": stage_fits, "orbitize": stage_orbitize,
              "stability": stage_stability}
    t0 = time.time()
    if a.stage == "all":
        for name in ("data", "fits", "orbitize", "stability"):
            print(f"== {name}", flush=True)
            stages[name]()
        import report
        report.write_all()
    elif a.stage in stages:
        stages[a.stage]()
    else:
        import report
        getattr(report, f"stage_{a.stage}")()
    print(f"done in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
