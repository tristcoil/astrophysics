"""The whole M4 pipeline, in order. Every step caches its result, so a rerun
only does what is missing.

    python run.py              # everything
    python run.py --list       # the steps
    python run.py figures      # one or more steps by name

The nested-sampling fits run as `python switch_fit.py fit NAME [TAG]`
subprocesses (that is how their checkpoints were written; a killed fit
resumes from its checkpoint).
"""

from __future__ import annotations

import subprocess
import sys
import time

import numpy as np

import config as C

PY = sys.executable
SWITCH = ("SS Lac", "HS Hya", "V907 Sco")


def step_fetch():
    import fetch_data
    fetch_data.ensure_data()


def step_tests():
    subprocess.run([PY, "test_physics.py"], check=True)


def step_etv():
    import etv
    if not (C.DERIVED / "tic290061484" / "eclipse_times.csv").exists():
        etv.run()


def step_tic290_fit():
    import tic290_nbody as N
    for v in ("triple", "all"):
        r = N.run_fit(variant=v)
        print(f"TIC 290061484 timing fit ({v}): chi2 {r['chi2']:.0f} for {r['dof']} dof")


def step_profiles():
    import switch_fit as F
    for name in SWITCH:
        F.profiles(name)


#: the check runs of section 5: the other sense of motion, repeats with
#: another seed, and SS Lac outside the Kozai-Lidov window
CHECK_RUNS = {"V907 Sco": ("retro",), "HS Hya": ("rerun", "all", "all_retro", "all_pro", "all_rerun"),
              "SS Lac": ("quietpro", "quietretro")}


def _fit(name, tag="", nproc=1):
    import switch_fit as F
    if F.posterior_path(name, tag).exists():
        return
    subprocess.run([PY, "-W", "ignore", "switch_fit.py", "fit", name] + ([tag] if tag else [])
                   + ([f"--nproc={nproc}"] if nproc > 1 else []), check=True)


def step_switch_fits():
    """Plates + TESS for all three, ground-data profiles, then the check runs."""
    import switch_fit as F
    for name in SWITCH:
        _fit(name)
    for name in SWITCH:
        post = np.load(F.posterior_path(name), allow_pickle=True)["posterior"]
        F.ground_profiles(name, np.median(post, axis=0))
    for name, tags in CHECK_RUNS.items():
        for tag in tags:
            _fit(name, tag, nproc=4 if name == "SS Lac" else 1)


def step_nbody_checks():
    """The final solutions, HS Hya's retrograde alternative, and SS Lac's
    unrestricted fit (which fails: Kozai-Lidov)."""
    import switch_fit as F
    import switch_nbody as SN
    for name in SWITCH:
        SN.run_checks(name, tag=F.FINAL_TAG[name])
    SN.run_checks("HS Hya", tag="all_retro")
    unres = F.CACHE / "nbody_check_ss_lac_unrestricted.json"
    if not unres.exists():
        SN.run_checks("SS Lac", tag="")                 # writes nbody_check_ss_lac.json
        (F.CACHE / "nbody_check_ss_lac.json").rename(unres)


def step_kozai():
    import kozai as K
    for lab in ("newton", "gr", "tides", "tides_Bsoft", "tides_Aonly"):
        K.run(lab, gr=lab != "newton", tides=lab.startswith("tides"), k2=K.VARIANTS.get(lab))


def step_stability():
    import stability as St
    for kind in ("ratio_e", "ratio_i"):
        St.megno_map(kind)
    St.disruption_ensemble()


def step_results():
    import analysis
    analysis.write_all()


def step_figures():
    import figures
    for name, f in figures.FIGURES.items():
        print(name, f())


def step_animation():
    """Kostov et al.'s solution with the inner period calibrated: it fits our
    eclipse times to 3.3 min and, unlike our timing-only refit, keeps B's
    inclination as their light curves (B's own eclipses) fix it."""
    import animate
    import tic290_nbody as N
    print(animate.make(params=N.run_fit(variant="triple")["base"],
                       label="Kostov et al. (2024)'s solution, inner period calibrated to our eclipse times"))


def step_explorer():
    import interactive
    print(interactive.build())


STEPS = [("fetch", step_fetch), ("tests", step_tests), ("etv", step_etv), ("tic290_fit", step_tic290_fit),
         ("profiles", step_profiles), ("switch_fits", step_switch_fits), ("nbody_checks", step_nbody_checks),
         ("kozai", step_kozai), ("stability", step_stability), ("results", step_results),
         ("figures", step_figures), ("animation", step_animation), ("explorer", step_explorer)]

if __name__ == "__main__":
    args = sys.argv[1:]
    if args == ["--list"]:
        for name, f in STEPS:
            print(f"{name:14s} {((f.__doc__ or '').strip().splitlines() or [''])[0]}")
        sys.exit()
    names = args or [n for n, _ in STEPS]
    lookup = dict(STEPS)
    for n in names:
        if n not in lookup:
            sys.exit(f"unknown step {n!r}; see --list")
    for n in names:
        t0 = time.time()
        print(f"== {n}", flush=True)
        lookup[n]()
        print(f"   {n} done in {(time.time() - t0) / 60:.1f} min", flush=True)
