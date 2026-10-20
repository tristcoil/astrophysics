"""Regenerate every result and figure of project H4 from the cached data.

    python fetch_data.py     # once: ~0.4 GB from the Gaia archive (~15 min)
    python run.py            # ~30-40 min the first time (bootstrap, l_max scan)
    python run.py --force    # recompute every cached step
    python run.py --only baseline,lmax   # recompute just these steps

Writes data/results.json and figures/. Each analysis step is cached as JSON
in data/cache/, so a rerun that only redraws figures takes seconds.
"""

from __future__ import annotations

import argparse
import json
import time

import numpy as np

import config as C


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.bool_):
        return bool(o)
    raise TypeError(type(o))


def cached(name: str, func, force: bool):
    path = C.CACHE_DIR / f"{name}.json"
    if path.exists() and not force:
        return json.loads(path.read_text())
    t0 = time.time()
    print(f"== {name}", flush=True)
    out = func()
    path.write_text(json.dumps(out, default=_json_default) + "\n")
    print(f"   {name}: {time.time() - t0:.0f} s", flush=True)
    return json.loads(path.read_text())


def main(force: bool = False, only: set[str] | None = None, figures_too: bool = True) -> dict:
    t0 = time.time()
    import analysis
    import galaxy
    import prepare

    df = prepare.build(force=force and (only is None or "prepare" in only))
    print(f"{len(df):,} Gaia-CRF3 sources, {int(df['five_param'].sum()):,} five-parameter", flush=True)

    def f(name):
        return force and (only is None or name in only) or (only is not None and name in only)

    base = cached("baseline", lambda: analysis.baseline(df), f("baseline"))
    k_out = base["outlier_source_ids"]
    lmax = cached("lmax_scan", lambda: analysis.lmax_scan(df), f("lmax"))
    choices = cached("choices", lambda: analysis.choice_tests(df), f("choices"))
    subsets = cached("subsets", lambda: analysis.subset_fits(df), f("subsets"))
    dr16q = cached("dr16q", lambda: analysis.dr16q_redshift_fits(df, base["vsh_coefficients"]), f("dr16q"))
    errs = cached("error_model", lambda: analysis.error_model(df), f("error_model"))
    errfit = cached("error_model_refit",
                    lambda: analysis.error_model_refit(df, errs["k_model"]["coef"]), f("error_model"))
    contam = cached("contamination", lambda: analysis.contamination(df, k_out), f("contamination"))
    emerge = cached("emergence", lambda: analysis.emergence(df), f("emergence"))
    zp_path = C.CACHE_DIR / "zp_l21.npy"
    if not zp_path.exists():
        np.save(zp_path, analysis.lindegren_zero_point(df))
    zp = np.load(zp_path)
    cond = cached("conditioning", lambda: analysis.conditioning_tests(df, zp), f("conditioning"))
    redshift = cached("redshift_quaia", lambda: analysis.redshift_tests(analysis.attach_quaia(df), zp),
                      f("redshift"))
    redshift_null = cached("redshift_null", lambda: analysis.redshift_null_tests(analysis.attach_quaia(df), zp),
                           f("redshift"))
    extended = cached("extendedness", lambda: analysis.extendedness_tests(analysis.attach_quaia(df),
                                                                          base["vsh_coefficients"]),
                      f("extended"))
    frame = cached("frame_rotation", lambda: analysis.frame_rotation_tests(df), f("frame"))
    lookalikes = cached("lookalikes", lambda: analysis.lookalike_tests(df), f("lookalikes"))

    g_gal = np.array(base["bootstrap"]["glide_galactic"]["vec"])
    g_icrs = np.array(base["bootstrap"]["glide_icrs"]["vec"])
    w_icrs = np.array(base["bootstrap"]["rotation_icrs"]["vec"]) if "rotation_icrs" in base["bootstrap"] \
        else np.array(base["scaled"]["rotation_icrs"]["vec"])
    five = df["five_param"].to_numpy()
    outl = np.isin(df["source_id"].to_numpy(), np.asarray(k_out))

    clean = analysis.remove_nonglide(df, five & ~outl, base["vsh_coefficients"])

    def maps():
        import healpy as hp
        out = {}
        for variant, frame in (("raw", df), ("clean", clean)):
            for nside, radius in ((2, 45.0), (4, 30.0)):
                th, ph = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
                centres = np.column_stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
                out[f"{variant}_nside{nside}_r{radius:.0f}"] = {
                    "nside": nside, "radius_deg": radius,
                    "l": np.degrees(ph).tolist(), "b": (90 - np.degrees(th)).tolist(),
                    **analysis.cap_averages(frame, five, centres, radius, model_glide_icrs=g_icrs,
                                            exclude=outl),
                }
        out["nonglide_rms"] = analysis.nonglide_rms(df, five & ~outl, base["vsh_coefficients"])
        return out

    skymaps = cached("sky_maps", maps, f("maps"))
    profile = cached("drift_profile", lambda: {
        "raw_minus_rotation": analysis.drift_profile(df, five, g_icrs, w_icrs, exclude=outl),
        "clean": analysis.drift_profile(clean, five, g_icrs, None, exclude=outl)}, f("maps"))

    # The Galaxy: predictions, and what our acceleration implies.
    preds = cached("predictions", lambda: {
        "constants": galaxy.predictions(), "mcmillan17": galaxy.mcmillan17_acceleration(),
        "lmc": galaxy.lmc_pull()}, f("galaxy"))
    boot_gal = np.load(C.CACHE_DIR / "baseline_bootstrap_glide_gal.npy")
    rng = np.random.default_rng(C.SEED)
    gx = boot_gal[:, 0]
    # Draws of g_X with the bootstrap's spread (a Gaussian with its mean/sd).
    g_draws = rng.normal(g_gal[0], np.std(gx), 200_000)
    implied = cached("implications", lambda: galaxy.implications(g_draws), f("galaxy"))

    results = {
        "baseline": {k: base[k] for k in ("n_sample", "n_bootstrap", "bootstrap", "scaled", "sandwich",
                                          "formal", "bootstrap_over_scaled_err", "grade")},
        "lmax_scan": lmax, "choices": choices, "subsets": subsets, "dr16q": dr16q,
        "error_model": {k: v for k, v in errs.items() if not isinstance(v, dict) or "hist_counts" not in v}
        | {k: {kk: vv for kk, vv in v.items() if not kk.startswith("hist")}
           for k, v in errs.items() if isinstance(v, dict) and "hist_counts" in v},
        "error_model_refit": errfit, "contamination": contam,
        "conditioning": cond, "redshift_quaia": redshift, "redshift_null": redshift_null,
        "extendedness": extended,
        "frame_rotation": frame,
        "lookalikes": lookalikes,
        "drift_profile": profile, "predictions": preds, "implications": implied,
        "runtime_s": None,
    }

    if not figures_too:
        print(f"analysis done in {time.time() - t0:.0f} s (figures skipped)")
        return results
    print("== figures", flush=True)
    import figures
    figures.make_all(df, results, skymaps, emerge, boot_gal, errs)
    import interactive
    interactive.build(df, results, outl)

    results["runtime_s"] = round(time.time() - t0, 1)
    C.RESULTS_PATH.write_text(json.dumps(results, indent=1, default=_json_default) + "\n")
    print(f"done in {time.time() - t0:.0f} s")
    return results


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--force", action="store_true")
    p.add_argument("--only", default=None, help="comma-separated steps to recompute")
    p.add_argument("--no-figures", action="store_true", help="run the analysis steps only")
    a = p.parse_args()
    main(force=a.force, only=set(a.only.split(",")) if a.only else None, figures_too=not a.no_figures)
