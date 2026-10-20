"""Every measurement of the project, gathered for the figures and the summary.

Each function returns plain numbers (for ``results.json``) and, where a figure
needs them, arrays. Slow steps (the 4-D posterior grid, the synchrotron
spectra) are cached in ``data/derived/`` keyed by their settings.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict

import numpy as np

import catalogs
import config
import cosmology
import lengths
import radiomaps as rm
import spectra
from config import DERIVED_DIR, RESULTS_PATH


def _key(*parts) -> str:
    return hashlib.sha1(json.dumps(parts, sort_keys=True, default=str).encode()).hexdigest()[:12]


# --------------------------------------------------------------------------
# The census
# --------------------------------------------------------------------------

def census() -> dict:
    g = catalogs.all_giants()
    o = catalogs.oei2023()
    lp = g["lp_mpc"].to_numpy()
    z = g["z"].to_numpy()
    minang, zmin = cosmology.minimum_angle_arcmin(config.GIANT_THRESHOLD_MPC, config.PLANCK18)
    return {
        "n_giants": int(len(g)),
        "by_method": {k: int(v) for k, v in g["found_by"].value_counts().items()},
        "n_oei2023": int(len(o)),
        "oei2023_spec_z_fraction": float(o["z_spec"].mean()),
        "median_lp_mpc": float(np.median(lp)),
        "lp_p10_p90": [float(v) for v in np.percentile(lp, [10, 90])],
        "n_above": {str(t): int((lp >= t).sum()) for t in (1, 2, 3, 4, 5)},
        "median_z": float(np.nanmedian(z)),
        "z_p10_p90": [float(v) for v in np.nanpercentile(z, [10, 90])],
        "n_z_above_1": int((z > 1).sum()),
        "max_z": float(np.nanmax(z)),
        "largest_angle_arcmin": float(g["las_arcmin"].max()),
        "n_larger_than_moon": int((g["las_arcmin"] > config.MOON_DIAMETER_ARCMIN).sum()),
        "min_angle_0p7mpc_arcmin": minang, "min_angle_redshift": zmin,
        "cosmology": lengths.cosmology_effect(),
    }


# --------------------------------------------------------------------------
# The length distribution
# --------------------------------------------------------------------------

def posterior(z_max: float, *, refresh: bool = False):
    axes = lengths.default_axes()
    lp = lengths.oei_sample(z_max)
    path = DERIVED_DIR / f"posterior_z{z_max:g}_{_key(z_max, len(lp), float(lp.sum()), {k: v.tolist() for k, v in axes.items()}, config.LENGTH_MODEL)}.npz"
    if path.exists() and not refresh:
        d = np.load(path)
        return lengths.Posterior(d["xi"], d["half_width"], d["b_ref"], d["sigma_ref"], d["log_post"],
                                 int(d["n"]), float(d["z_max"]))
    print(f"    posterior grid for z < {z_max} ({lp.size} giants) ...", flush=True)
    post = lengths.grid_posterior(lp, z_max, axes=axes)
    np.savez_compressed(path, xi=post.xi, half_width=post.half_width, b_ref=post.b_ref,
                        sigma_ref=post.sigma_ref, log_post=_compact(post.log_post), n=post.n, z_max=post.z_max)
    return post


def _compact(log_post: np.ndarray) -> np.ndarray:
    """Relative log-posterior as float32, with everything below e^-40 of the
    peak set to -inf: identical results, and a cache of megabytes instead of
    half a gigabyte."""
    rel = (log_post - np.max(log_post)).astype(np.float32)
    rel[rel < -40] = -np.inf
    return rel


def length_distribution() -> tuple[dict, dict]:
    res, art = {}, {}
    for z_max in config.LENGTH_MODEL["z_max"]:
        lp = lengths.oei_sample(z_max)
        post = posterior(z_max)
        xi_naive, xi_naive_err = lengths.pareto_mle(lp)
        summ = {k: post.summary(k) for k in ("xi", "half_width", "b_ref", "sigma_ref")}
        res[str(z_max)] = {"n": int(lp.size), "naive_xi": xi_naive, "naive_xi_err": xi_naive_err,
                           "posterior": summ, "map": post.map_point(),
                           "published": config.LENGTH_MODEL["published"][z_max]}
        if z_max == 0.5:
            pred = lengths.posterior_predictive(post, n_draws=150)
            art["post"] = post
            art["pred"] = pred
            art["lp"] = lp
    # Porphyrion in 3D: deprojected length for the posterior tail index
    xi_mean = res["0.5"]["posterior"]["xi"]["mean"]
    lp_p = 6.43  # Oei et al. 2024
    res["porphyrion_deprojection"] = {
        "xi_used": xi_mean,
        **{k: float(v) for k, v in lengths.deprojection(lp_p, xi_mean).items()},
        "published_mean_mpc": 7.28, "published_median_mpc": 6.8,
        "by_xi": {f"{x:.1f}": float(lengths.deprojection(lp_p, x)["mean"]) for x in (-3.0, -3.5, -4.0, -4.5)},
    }
    # what the tail index says about giants of Porphyrion's size
    xi_s = res["0.5"]["posterior"]["xi"]
    res["tail"] = {"fraction_above_5mpc_of_giants": float(lengths.survival(5.0, xi_s["mean"])),
                   "fraction_above_porphyrion": float(lengths.survival(lp_p, xi_s["mean"])),
                   "fraction_above_1mpc": float(lengths.survival(1.0, xi_s["mean"]))}
    return res, art


# --------------------------------------------------------------------------
# Our own lengths from the images
# --------------------------------------------------------------------------

def _ext_summary(ext: rm.Extent, g: config.Giant) -> dict:
    lp = float(cosmology.projected_length_mpc(ext.las_arcmin, g.z, config.PLANCK18))
    p2p = (float(cosmology.projected_length_mpc(ext.peak_to_peak_arcmin, g.z, config.PLANCK18))
           if ext.peak_to_peak_arcmin else None)
    return {"beam_arcsec": ext.beam_arcsec, "sigma_ujy": ext.sigma_jy * 1e6, "las_arcmin": ext.las_arcmin,
            "lp_mpc": lp, "peak_to_peak_arcmin": ext.peak_to_peak_arcmin, "peak_to_peak_mpc": p2p,
            "pa_deg": ext.pa_deg, "flux_jy": ext.flux_jy, "ends": ext.las_endpoints,
            "complete": not any(s in ext.missing_seeds for s in g.seeds)}


def measure_all() -> tuple[list[dict], dict]:
    """Angular and projected lengths of every imaged giant, at 6", 20", 60"."""
    rows, art = [], {}
    for g in config.GIANTS:
        base20 = rm.load_lotss(g.key, 20)
        exts = {"20": rm.measure_extent(base20, g),
                "60": rm.measure_extent(rm.smooth_to(base20, 60.0), g)}
        try:
            base6 = rm.load_lotss(g.key, 6)
            exts["6"] = rm.measure_extent(base6, g)
            del base6
        except FileNotFoundError:
            pass  # mosaicked giants are fetched at 20" only
        if g.measure_beam not in (20.0, 60.0):
            exts[f"{g.measure_beam:g}"] = rm.measure_extent(rm.smooth_to(base20, g.measure_beam), g)
        adopted = exts[f"{g.measure_beam:g}"]
        row = {"key": g.key, "name": g.name, "z": g.z, "z_err": g.z_err, "z_kind": g.z_kind,
               "groups": list(g.groups), "measure_beam": g.measure_beam, "survey": g.survey,
               "at": {k: _ext_summary(e, g) for k, e in exts.items()},
               "adopted": _ext_summary(adopted, g),
               "claims": [{"source": c[0], "las_arcmin": c[1], "lp_mpc": c[2]} for c in g.claims]}
        # the length in each cosmology, and its uncertainty from the redshift
        row["lp_by_cosmology"] = {c.name: float(cosmology.projected_length_mpc(adopted.las_arcmin, g.z, c))
                                  for c in config.COSMOLOGIES}
        if g.alt_seeds:
            from dataclasses import replace

            gc = replace(g, seeds=g.alt_seeds, extra_seeds=())
            # the finest resolution that shows both ends: smoothing can merge
            # unrelated neighbours into the host's own emission
            ext_c = rm.measure_extent(base20, gc)
            if any(sd in ext_c.missing_seeds for sd in gc.seeds):
                ext_c = rm.measure_extent(rm.smooth_to(base20, 60.0), gc)
            row["conservative"] = _ext_summary(ext_c, gc)
            art.setdefault(g.key, {})["conservative"] = ext_c
        if g.z_err > 0:
            lo = cosmology.projected_length_mpc(adopted.las_arcmin, max(g.z - g.z_err, 1e-4), config.PLANCK18)
            hi = cosmology.projected_length_mpc(adopted.las_arcmin, g.z + g.z_err, config.PLANCK18)
            row["lp_z_err_mpc"] = float(abs(hi - lo) / 2)
        rows.append(row)
        art.setdefault(g.key, {}).update({"ext": exts, "adopted": adopted})
    return rows, art


def compare_releases(rows: list[dict]) -> dict:
    """The DR2 giants measured again, with the same seeds, on their DR3 20"
    images: does the newer release change their lengths or brightness?"""
    out = []
    for r in rows:
        g = config.giant_by_key(r["key"])
        if g.survey != "dr2":
            continue
        try:
            d3 = rm.load_lotss(g.key, 20, check=True)
        except FileNotFoundError:
            continue
        img = d3 if g.measure_beam < 60 else rm.smooth_to(d3, 60.0)
        e3 = rm.measure_extent(img, g)
        d2 = rm.load_lotss(g.key, 20)
        img2 = d2 if g.measure_beam < 60 else rm.smooth_to(d2, 60.0)
        e2 = rm.measure_extent(img2, g)
        ok = not any(sd in e3.missing_seeds for sd in g.seeds) and not any(sd in e2.missing_seeds for sd in g.seeds)
        # flux inside the same sky footprint: the DR2 mask, regridded onto DR3
        from reproject import reproject_interp

        m2, _ = reproject_interp((e2.mask.astype(float), img2.wcs), img.wcs, shape_out=img.data.shape, order="bilinear")
        m2 = np.nan_to_num(m2) > 0.5
        f3_same = float(np.nansum(np.where(m2, img.data, 0.0)) / img.beam_area_pix)
        out.append({"key": g.key, "name": g.name, "complete": bool(ok), "beam_arcsec": img.beam_arcsec,
                    "las_dr2": e2.las_arcmin, "las_dr3": e3.las_arcmin,
                    "sigma_dr2_ujy": e2.sigma_jy * 1e6, "sigma_dr3_ujy": e3.sigma_jy * 1e6,
                    "flux_dr2_jy": e2.flux_jy, "flux_dr3_same_footprint_jy": f3_same,
                    "flux_ratio_dr3_dr2": f3_same / e2.flux_jy if e2.flux_jy > 0 else np.nan})
    good = [o for o in out if o["complete"]]
    lr = np.array([o["las_dr3"] / o["las_dr2"] for o in good])
    fr = np.array([o["flux_ratio_dr3_dr2"] for o in good if np.isfinite(o["flux_ratio_dr3_dr2"])])
    nr = np.array([o["sigma_dr3_ujy"] / o["sigma_dr2_ujy"] for o in good])
    return {"giants": out, "n_complete": len(good),
            "las_ratio_median": float(np.median(lr)) if lr.size else np.nan,
            "las_ratio_range": [float(lr.min()), float(lr.max())] if lr.size else None,
            "flux_ratio_median": float(np.median(fr)) if fr.size else np.nan,
            "flux_ratio_range": [float(fr.min()), float(fr.max())] if fr.size else None,
            "noise_ratio_median": float(np.median(nr)) if nr.size else np.nan}


# --------------------------------------------------------------------------
# Lobe spectra and ages
# --------------------------------------------------------------------------

#: Lobes start this far from the host (arcmin, or 2 NVSS beams if larger), so
#: the flat-spectrum core and inner jets are not counted as lobe.
CORE_RADIUS_ARCMIN = 1.5


def ageing(g: config.Giant, ext: rm.Extent) -> tuple[dict, dict]:
    sm = spectra.spectral_maps(g)
    # the giant's footprint at 45": 3 sigma in the smoothed LoTSS map, using
    # the same seeds
    img45 = rm.RadioImage(sm.lotss, sm.wcs, sm.beam_arcsec, sm.pixel_arcsec, spectra.NU_LO, "LoTSS 45\"")
    ext45 = rm.measure_extent(img45, g, sigma=sm.sigma_lotss)
    mask = ext45.mask
    end_a, end_b = g.seeds
    along, across = spectra.axis_coordinates(sm, g, end_a)
    r_core = max(CORE_RADIUS_ARCMIN, 2 * sm.beam_arcsec / 60)
    kpc_per_arcmin = float(cosmology.kpc_per_arcmin(g.z, config.PLANCK18))
    out = {"sigma_lotss_mjy": sm.sigma_lotss * 1e3, "sigma_nvss_mjy": sm.sigma_nvss * 1e3,
           "core": spectra.region_flux(sm, mask & (np.hypot(along, across) < r_core)), "lobes": {}}
    profiles = {}
    b_ic = config.B_CMB0_UG * (1 + g.z) ** 2
    out["b_cmb_ug"] = b_ic
    for name, sign in (("A", 1), ("B", -1)):
        lobe = mask & (sign * along > r_core)
        if lobe.sum() < sm.beam_area_pix:
            continue
        reg = spectra.region_flux(sm, lobe)
        a_vals = sign * along[lobe]
        length_arcmin = float(a_vals.max() - a_vals.min()) + sm.beam_arcsec / 60
        width_arcmin = float(np.percentile(np.abs(across[lobe]), 90) * 2)
        width_arcmin = max(width_arcmin, sm.beam_arcsec / 60)
        alpha_b = reg["alpha"] if np.isfinite(reg["alpha"]) else -1.0
        b_eq = spectra.equipartition_field_ug(reg["flux_144_jy"], spectra.NU_LO, min(alpha_b, -0.55), g.z,
                                              length_arcmin * 60, width_arcmin * 60, width_arcmin * kpc_per_arcmin)
        reg.update({"length_arcmin": length_arcmin, "width_arcmin": width_arcmin, "b_eq_ug": b_eq,
                    "end_distance_arcmin": float(a_vals.max())})
        if np.isfinite(reg["alpha"]):
            nb_ci = spectra.break_from_index(reg["alpha"], config.ALPHA_INJ, "CI")
            reg["ci_break_mhz"] = nb_ci
            if not np.isfinite(nb_ci):
                reg["ci_note"] = ("flatter than injection: no ageing measurable" if reg["alpha"] > config.ALPHA_INJ
                                  else "steeper than continuous injection allows (NVSS likely misses flux)")
            reg["ci_age_myr"] = float(spectra.spectral_age_myr(nb_ci, b_eq, g.z)) if np.isfinite(nb_ci) else np.nan
            reg["ci_age_max_myr"] = (float(spectra.spectral_age_myr(nb_ci, spectra.max_age_field(g.z), g.z))
                                     if np.isfinite(nb_ci) else np.nan)
        out["lobes"][name] = reg
        # profile from the lobe's end toward the core, one beam per bin
        step = sm.beam_arcsec / 60
        edges = np.arange(r_core, a_vals.max() + step, step)
        prof = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = lobe & (sign * along >= lo) & (sign * along < hi)
            if m.sum() < 0.5 * sm.beam_area_pix:
                continue
            f = spectra.region_flux(sm, m)
            age = np.nan
            if np.isfinite(f["alpha"]):
                nb = spectra.break_from_index(f["alpha"], config.ALPHA_INJ, "JP")
                age = float(spectra.spectral_age_myr(nb, b_eq, g.z)) if np.isfinite(nb) else 0.0
            prof.append({"distance_from_host_arcmin": float((lo + hi) / 2),
                         "distance_from_end_arcmin": float(a_vals.max() - (lo + hi) / 2),
                         **{k: f[k] for k in ("flux_144_jy", "flux_1400_jy", "alpha", "alpha_err", "alpha_upper")},
                         "jp_age_myr": age})
        profiles[name] = prof
        out["lobes"][name]["profile"] = prof
        # advance speed: slope of age against distance from the end
        pts = [(p["distance_from_end_arcmin"] * kpc_per_arcmin, p["jp_age_myr"]) for p in prof
               if np.isfinite(p["jp_age_myr"]) and p["jp_age_myr"] > 0]
        if len(pts) >= 4:
            d, t = np.array(pts).T
            slope = np.polyfit(d, t, 1)[0]  # Myr per kpc
            if slope > 0:
                v_c = (1 / slope) * 3.0857e16 / 3.15576e13 / 2.998e5  # kpc/Myr -> c
                reg["advance_speed_c"] = float(v_c)
    return out, {"maps": sm, "mask45": mask, "along": along, "across": across}


def ageing_all(measured_art: dict) -> tuple[dict, dict]:
    res, art = {}, {}
    for g in config.giants_in("ageing"):
        r, a = ageing(g, measured_art[g.key]["adopted"])
        res[g.key] = r
        art[g.key] = a
    # 3C 236 lobes against Shulevski et al. (2019): NW is our end A (seed 1)
    if "3c236" in res:
        lobes = res["3c236"]["lobes"]
        res["3c236"]["shulevski2019"] = {
            "NW": {"flux_143_mjy": 5030, "flux_1400_mjy": 730, "ci_age_myr": 129.1, "b_ug": 1.28},
            "SE": {"flux_143_mjy": 3580, "flux_1400_mjy": 580, "ci_age_myr": 117.1, "b_ug": 1.28}}
        # their fluxes through our CI model, with their field: a check of the model itself
        chk = {}
        for lobe, f1, f2 in (("NW", 5.03, 0.73), ("SE", 3.58, 0.58)):
            a = np.log(f2 / f1) / np.log(1400 / 143)
            nb = spectra.break_from_index(a, -0.55, "CI")
            chk[lobe] = {"alpha": float(a), "ci_break_mhz": nb,
                         "ci_age_myr": float(spectra.spectral_age_myr(nb, 1.28, 0.1005))}
        res["3c236"]["model_check_with_their_fluxes"] = chk
        res["3c236"]["lobe_names"] = {"A": "NW", "B": "SE"}
    return res, art


# --------------------------------------------------------------------------

def _clean(obj):
    """Make results JSON-serialisable: arrays to lists, NaN to None."""
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items() if k != "values"}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _clean(obj.tolist())
    if isinstance(obj, (np.floating, float)):
        return None if not np.isfinite(obj) else float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def save_results(**sections) -> None:
    RESULTS_PATH.write_text(json.dumps(_clean(sections), indent=1))


# --------------------------------------------------------------------------
# Does the length-distribution machinery recover a known answer?
# --------------------------------------------------------------------------

COARSE_AXES = {"xi": np.arange(-6.0, -1.5 + 1e-9, 0.1), "half_width": np.arange(0, 4.01, 0.25),
               "b_ref": np.arange(50.0, 3501.0, 100.0), "sigma_ref": np.arange(0.2, 3.01, 0.2)}


def simulation_check(n_sims: int = 12, n_giants: int = 1473, z_max: float = 0.5, seed: int = 100) -> dict:
    """Simulated catalogues from the full forward model at Oei et al.'s
    published best-fit selection parameters and several true tail indices,
    each fitted exactly as the real data are (coarser grid, for speed)."""
    pub = config.LENGTH_MODEL["published"][z_max]
    truths = np.tile([-3.0, -3.5, -4.0], int(np.ceil(n_sims / 3)))[:n_sims]
    path = DERIVED_DIR / f"simcheck_{_key(n_sims, n_giants, z_max, seed, truths.tolist(), config.LENGTH_MODEL)}.json"
    if path.exists():
        return json.loads(path.read_text())
    out = []
    for i, xi_true in enumerate(truths):
        sim = lengths.simulate(n_giants, float(xi_true), pub["half_width_arcmin"][0], pub["b_ref"][0],
                               pub["sigma_ref"][0], z_max, seed=seed + i)
        post = lengths.grid_posterior(sim["lp"], z_max, axes=COARSE_AXES)
        s = post.summary("xi")
        out.append({"xi_true": float(xi_true), "xi_mean": s["mean"], "xi_sd": s["sd"], "xi_median": s["median"],
                    "xi_mode": s["mode"], "xi_p16": s["p16"], "xi_p84": s["p84"],
                    "naive_xi": lengths.pareto_mle(sim["lp"])[0]})
        print(f"    simulation {i + 1}/{n_sims}: true {xi_true:.1f}, recovered {s['mean']:.2f} +- {s['sd']:.2f}", flush=True)
    res = {"settings": {"n_giants": n_giants, "z_max": z_max, "half_width": pub["half_width_arcmin"][0],
                        "b_ref": pub["b_ref"][0], "sigma_ref": pub["sigma_ref"][0]}, "runs": out}
    for est in ("mean", "median", "mode"):
        pulls = np.array([(r[f"xi_{est}"] - r["xi_true"]) / r["xi_sd"] for r in out])
        res[f"{est}_pull"] = float(pulls.mean())
        res[f"{est}_rms_pull"] = float(np.sqrt(np.mean(pulls**2)))
        res[f"{est}_bias"] = float(np.mean([r[f"xi_{est}"] - r["xi_true"] for r in out]))
    res["coverage_68"] = float(np.mean([r["xi_p16"] <= r["xi_true"] <= r["xi_p84"] for r in out]))
    path.write_text(json.dumps(res, indent=1))
    return res


# --------------------------------------------------------------------------
# Catalogue redshifts against SDSS spectra of the same hosts
# --------------------------------------------------------------------------

def redshift_check() -> dict:
    """How often a giant's catalogue redshift disagrees with the SDSS
    spectrum taken within 3" of its host, and what that does to its length."""
    import hosts

    g = hosts.giant_hosts()
    s = g[g["has_spectrum"] & g["z_sdss"].notna()].copy()
    dz = np.abs(s["z"] - s["z_sdss"]) / (1 + s["z_sdss"])
    bad = s[dz > 0.1]
    lp_sdss = cosmology.projected_length_mpc(bad["las_arcmin"], bad["z_sdss"], config.PLANCK18)
    return {"n_with_spectrum": int(len(s)),
            "frac_within_0p01": float((dz <= 0.01).mean()),
            "n_off_by_10pct": int(len(bad)), "frac_off_by_10pct": float(len(bad) / len(s)),
            "n_off_by_30pct": int((dz > 0.3).sum()),
            "off_spectral_classes": {k: int(v) for k, v in bad["spec_class"].value_counts().items()},
            "off_by_method": {k: int(v) for k, v in bad["found_by"].value_counts().items()},
            "n_no_longer_giants": int((lp_sdss < config.GIANT_THRESHOLD_MPC).sum()),
            "n_catalogue_z_above_3": int((g["z"] > 3).sum()),
            "examples": [{"name": r["name"], "z_catalogue": float(r["z"]), "z_sdss": float(r["z_sdss"]),
                          "class": r["spec_class"]} for _, r in bad.sort_values("z").tail(4).iterrows()]}
