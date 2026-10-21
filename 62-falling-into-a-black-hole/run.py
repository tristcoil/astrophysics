"""Regenerate every result and figure of project N5 from the cached data.

    python fetch_data.py     # once: ~0.3 GB (Gaia DR3 stars + summed-flux map, Hipparcos, colour tables)
    python run.py            # figures + results (~10 min); the movie is separate:
    python movie.py          # ~2,300 frames, ~1-2 h on 6 cores (resumable), then H.264 encoding
    python run.py --force    # also rebuild the sky texture and star list from the raw files

Writes data/results.json and figures/.
"""

from __future__ import annotations

import json
import sys
import time

import config as C


def main(force: bool = False) -> dict:
    t0 = time.time()
    import sky
    if force or not sky.POINTS_CACHE.exists():
        print("== star list (Gaia DR3 G < 11 + Hipparcos V < 5)", flush=True)
        sky.build_points()
    if force or not sky.TEXTURE_CACHE.exists():
        print("== sky texture from the Gaia summed-flux map", flush=True)
        sky.build_texture()

    print("== numbers for each black hole", flush=True)
    import analysis
    res = analysis.summary()
    bhs, uni = res["black_holes"], res["universal"]

    print("== symbolic checks (SymPy)", flush=True)
    import symbolic
    sym = symbolic.derive()

    print("== figures", flush=True)
    import figures
    import render as RD
    figures.clocks(bhs)
    figures.tides(bhs, uni)
    figures.friend(bhs)
    figures.penrose(bhs)
    figures.river()
    figures.shadow_size()
    figures.kerr_figure(bhs)
    scene = RD.Scene()
    figures.views_sequence(scene, bhs)
    figures.hover_vs_fall(scene)
    figures.near_singularity(scene, bhs)
    figures.poster_frame()

    meta_path = C.CACHE_DIR / "movie_meta.json"
    out = {"black_holes": bhs, "universal": uni, "symbolic": sym,
           "movie": json.loads(meta_path.read_text()) if meta_path.exists() else None,
           "runtime_s": time.time() - t0}
    C.RESULTS_PATH.write_text(json.dumps(out, indent=1, default=float))
    print(f"done in {time.time() - t0:.0f} s -> {C.RESULTS_PATH}")
    return out


if __name__ == "__main__":
    main(force="--force" in sys.argv)
