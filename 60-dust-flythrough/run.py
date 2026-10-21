"""Regenerate every result, figure, the viewer and the movie of project H6 from the cached data.

    python fetch_data.py     # once: ~8.1 GB (the two 3D maps are 7.4 of it)
    python run.py            # ~1.5 h the first time on a busy machine (wave fits), minutes afterwards
    python run.py --force    # recompute the cached products too
    python run.py --no-movie # skip rendering the fly-through movie (needs a GPU browser)

Writes data/results.json, figures/*.png, figures/dust_flythrough.html (offline viewer),
figures/dust_viewer/ (the same viewer with separate data files) and figures/dust_flythrough.mp4.
"""

from __future__ import annotations

import json
import sys
import time

import config as C
import literature as L


def main(force: bool = False, movie: bool = True) -> dict:
    t0 = time.time()
    print("== Cartesian grids of both maps", flush=True)
    import grids
    grids.build_main(force=force)
    grids.build_wide(force=force)

    print("== checks against the published paper and its samples", flush=True)
    import validation
    path = C.CACHE_DIR / "validation.json"
    if force or not path.exists():
        val = validation.run()
    else:
        val = json.loads(path.read_text())

    print("== the second map (Vergely et al. 2022)", flush=True)
    import compare
    path = C.CACHE_DIR / "compare.json"
    cmp_ = compare.run() if force or not path.exists() else json.loads(path.read_text())

    print("== the Local Bubble", flush=True)
    import bubble
    import history
    shell = bubble.shell_model(force=force)
    bub = {"oneill_comparison": bubble.compare_oneill(shell), "sphere_fit": bubble.sphere_fit(shell),
           "interior": bubble.interior_volume(shell)}
    hist = history.fit(force=force)
    bub["history"] = hist
    bub["wall_vs_model_sphere"] = history.shell_vs_model(hist["centre_pc"], hist["R_now_pc"][1])
    cl = bubble.clusters_on_wall(shell, rmax=250)
    young = cl[cl["age_myr"] < 20]
    older = cl[(cl["age_myr"] >= 60) & (cl["age_myr"] < 300)]
    bub["clusters_within_250pc"] = {
        "young_lt20myr_n": int(len(young)), "young_frac_within_30pc": float((young["offset"].abs() < 30).mean()),
        "older_60_300myr_n": int(len(older)), "older_frac_within_30pc": float((older["offset"].abs() < 30).mean())}

    print("== the Radcliffe Wave", flush=True)
    import wave
    wave.cloud_distances(force=force)
    fits = wave.run_fits(force=force)

    print("== the 3D viewer", flush=True)
    import numpy as np
    import viewer
    x, y, z, _ = wave.kon_curve(fits[wave.MAIN_FIT]["median"], 400)
    vinfo = viewer.build(np.c_[x, y, z])

    print("== figures", flush=True)
    import figures
    figures.projections(labelled=False)
    figures.projections(labelled=True)
    for f in figures.ALL[1:]:
        f()

    if movie:
        print("== the fly-through movie", flush=True)
        import movie as mv
        mv.make(force=force)

    results = {
        "validation": val, "second_map": {k: v for k, v in cmp_.items() if k != "elongation"},
        "elongation": {name: [{k: r[k] for k in ("r_lo", "r_hi", "corr_len_radial_pc", "corr_len_tangential_pc",
                                                 "elongation")} for r in cmp_["elongation"][name]]
                       for name in ("edenhofer", "vergely")},
        "local_bubble": bub,
        "radcliffe_wave": {tag: {k: v for k, v in f.items() if k not in ("members", "odds")} for tag, f in fits.items()},
        "viewer": vinfo,
        "literature": {"edenhofer": L.EDENHOFER, "oneill": L.ONEILL, "zucker2022": L.ZUCKER2022,
                       "alves2020": L.ALVES2020, "konietzka2024": L.KONIETZKA2024},
        "runtime_s": time.time() - t0,
    }
    C.RESULTS_PATH.write_text(json.dumps(results, indent=1, default=float))
    print(f"done in {time.time() - t0:.0f} s -> {C.RESULTS_PATH}")
    return results


if __name__ == "__main__":
    main(force="--force" in sys.argv, movie="--no-movie" not in sys.argv)
