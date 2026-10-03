"""Every measurement in the project, from the cached data to ``results.json``.

``analyse(paths)`` returns two things: ``products``, the in-memory material the
figures are drawn from (maps, light curves, grids), and ``results``, every number
the README quotes, JSON-ready. The products are also pickled into
``data/<epoch>/derived/`` so figures can be redrawn without recomputing.
"""

from __future__ import annotations

import json
import pickle
import warnings
from datetime import datetime, timezone

import astropy.units as u
import numpy as np
import pandas as pd

import config as C
import fetch_data
import flare
import goes
import limb
import rotation
import sdo
import solar_physics as sp

warnings.filterwarnings("ignore")


def _log(msg: str, quiet: bool) -> None:
    if not quiet:
        print(msg, flush=True)


# --- The snapshot ---------------------------------------------------------------------------


def snapshot(paths: C.EpochPaths, manifest: dict) -> tuple[dict, dict]:
    """The Sun at the epoch: every channel as a DN/s map, and the magnetogram."""
    maps = {c.wavelength: sdo.load_aia(paths.data / manifest["snapshot"][str(c.wavelength)])
            for c in C.CHANNELS}
    hmi = sdo.load_map(paths.data / manifest["snapshot"]["hmi"])
    res = {
        "times": {str(wl): m.date.isot for wl, m in maps.items()} | {"hmi": hmi.date.isot},
        "dsun_au": float(maps[171].dsun.to_value(u.AU)),
        "rsun_arcsec": float(maps[171].rsun_obs.to_value(u.arcsec)),
        "b0_deg": float(maps[171].observer_coordinate.lat.to_value(u.deg)),
        "carrington_rotation": int(maps[171].meta.get("car_rot", 0)),
        "l0_deg": float(maps[171].meta.get("crln_obs", np.nan)),
    }
    return {"maps": maps, "hmi": hmi}, res


# --- Active regions and Hale's law ----------------------------------------------------------


def regions(paths: C.EpochPaths, manifest: dict, hmi_week: list, quiet: bool) -> tuple[list, dict]:
    """Detect, track and name the week's active regions; test Hale's law on them."""
    _log("  tracking active regions through the week", quiet)
    dets = [rotation.detect_regions(m) for m in hmi_week]
    tracks = rotation.link_tracks(dets)
    summaries = [rotation.track_summary(t, paths.week_start.replace(tzinfo=None)) for t in tracks]
    noaa = rotation.noaa_regions(paths.data / manifest["goes"]["regions"])
    rotation.name_tracks(summaries, noaa, paths.week_start, paths.epoch)
    # One entry per NOAA region: its longest track (fragments of the same region, a
    # few frames where it split in two, are dropped)
    best = {}
    for s in summaries:
        if "noaa" not in s:
            continue
        if s["noaa"] not in best or s["n"] > best[s["noaa"]]["n"]:
            best[s["noaa"]] = s
    named = sorted(best.values(), key=lambda s: s["noaa"])
    # Hale's law: leading polarity against hemisphere, for regions that are true
    # bipoles (the weaker polarity carries at least a fifth of the stronger's flux)
    hale = []
    for s in named:
        bal = s.get("flux_balance")
        if bal is None or not (0.2 <= bal <= 5.0):
            continue
        expected = C.CYCLE25_LEADING_POLARITY[s["hemisphere"]]
        hale.append({"noaa": s["noaa"], "hemisphere": s["hemisphere"], "lat": s["lat"],
                     "leading": s["leading_polarity"], "obeys": s["leading_polarity"] == expected,
                     "tilt_deg": s.get("tilt_deg"), "flux_balance": bal})
    res = {
        "n_tracks": len(summaries),
        "noaa_regions": [{k: v for k, v in s.items() if k != "track"} for s in named],
        "hale": hale,
        "hale_obeyed": int(sum(h["obeys"] for h in hale)),
        "hale_tested": len(hale),
        "joy_positive_tilts": int(sum((h["tilt_deg"] or 0) > 0 for h in hale)),
    }
    return summaries, res


# --- Rotation ----------------------------------------------------------------------------------


def differential_rotation(hmi_week: list, quiet: bool) -> tuple[list, dict]:
    _log("  remapping the week's magnetograms onto Carrington grids", quiet)
    from sunpy.coordinates import sun

    mid = hmi_week[len(hmi_week) // 2].date
    l0_ref = float(sun.L0(mid).to_value(u.deg))
    grids = [rotation.remap(m, l0_ref) for m in hmi_week]
    _log("  cross-correlating latitude bands", quiet)
    main = rotation.rotation_profile(grids)
    alts = {}
    base = C.ROT_PAIR_SEPARATION_H
    for sep in C.ROT_ALT_SEPARATIONS_H:
        C.ROT_PAIR_SEPARATION_H = sep
        try:
            r = rotation.rotation_profile(grids)
        finally:
            C.ROT_PAIR_SEPARATION_H = base
        if "ABC" in r["fits"]:
            lats = [b["lat"] for b in r["bands"]]
            alts[str(sep)] = {"A": r["fits"]["ABC"]["coef"][0], "A_err": r["fits"]["ABC"]["err"][0],
                              "B": r["fits"]["ABC"]["coef"][1], "C": r["fits"]["ABC"]["coef"][2],
                              "lat_range": [min(lats), max(lats)], "n_bands": len(lats),
                              "median_peak_corr": float(np.median([b["median_peak_corr"]
                                                                   for b in r["bands"]]))}
    a_all = [main["fits"]["ABC"]["coef"][0]] + [v["A"] for v in alts.values()]
    main["alternative_separations"] = alts
    main["A_systematic"] = float((max(a_all) - min(a_all)) / 2)
    main["pair_separation_h"] = base
    coef = main["fits"]["ABC"]["coef"]
    main["sidereal_period_equator_d"] = 360.0 / coef[0]
    w_earth = sp.earth_orbital_rate_deg_per_day(hmi_week[-1].date)
    main["earth_orbital_rate_deg_per_day"] = w_earth
    main["synodic_period_equator_d"] = sp.synodic_period_days(coef[0], w_earth)
    for lat in (30, 45):
        om = float(sp.rotation_law(lat, *coef))
        main[f"sidereal_period_{lat}_d"] = 360.0 / om
    return grids, main


# --- The corona and the limb ------------------------------------------------------------------


def limb_physics(paths: C.EpochPaths, manifest: dict, maps: dict, quiet: bool) -> dict:
    out = {}
    if manifest.get("continuum"):
        _log("  limb darkening on the HMI continuum image", quiet)
        out["continuum"] = limb.continuum_limb_darkening(paths.data / manifest["continuum"])
    _log("  radial profiles and coronal scale heights", quiet)
    out["profiles"] = limb.radial_profiles(maps)
    out["corona"] = limb.corona_scale_heights(maps)
    return out


# --- X-rays and the flare -------------------------------------------------------------------------


def xray_week(paths: C.EpochPaths) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    xrs = goes.load_xrs(paths)
    events = goes.flare_list(paths, xrs)
    valid = xrs["long"].dropna()
    daily_min = valid.resample("1D").min()
    counts = {}
    for letter, base in C.FLARE_CLASSES:
        counts[letter] = int(((events.peak_flux >= base) & (events.peak_flux < base * 10)).sum())
    res = {
        "satellite": xrs.attrs["satellite"],
        "list_source": events.attrs["source"],
        "n_events": int(len(events)),
        "counts_by_class": counts,
        "largest": {"class": events.loc[events.peak_flux.idxmax(), "class"],
                    "swpc_class": events.loc[events.peak_flux.idxmax(), "swpc_class"],
                    "peak": events.loc[events.peak_flux.idxmax(), "peak"].isoformat(),
                    "peak_flux": float(events.peak_flux.max())},
        "background_median": float(valid.median()),
        "background_class": goes.flare_class(float(valid.median())),
        "daily_minimum_class": [goes.flare_class(v) for v in daily_min.values],
        "coverage_fraction": float(valid.size / len(xrs)),
        "science_fraction": float((xrs["source"] == "ncei").mean()),
        "cross_check_swpc_vs_ncei": goes.cross_check(xrs),
        "events": [
            {"begin": r.begin.isoformat(), "peak": r.peak.isoformat(), "end": r.end.isoformat(),
             "class": r["class"], "swpc_class": r.swpc_class,
             "ncei_class": r.get("ncei_class"), "fluence_J_m2": r.fluence_J_m2}
            for _, r in events.iterrows()],
    }
    return xrs, events, res


def flare_of_the_week(paths: C.EpochPaths, manifest: dict, xrs: pd.DataFrame,
                      noaa: list, quiet: bool) -> tuple[dict, dict]:
    sel = manifest["flare"].get("selection")
    if not sel:
        return {}, {"selection": None}
    _log("  the flare of the week: thermal history, arcade light curves", quiet)
    table = paths.data / manifest["goes"]["response"]
    thermal = goes.flare_thermal(xrs, sel, table)
    mask = flare.arcade_mask(paths, manifest)
    curves = flare.light_curves(paths, manifest, mask)
    t_last = max(df.index[-1] for df in curves.values())
    t_end = flare.episode_end(xrs, sel, t_last)
    seq = flare.cooling_sequence(curves, sel, t_end)
    area = flare.hot_area_cm2(paths, manifest, pd.Timestamp(sel["peak"]))
    energy = flare.energetics(thermal, area, sel)
    loc = dict(manifest["flare"]["location"])
    loc["noaa"] = _flare_region(loc, sel, noaa)
    res = {"selection": sel, "location": loc, "episode_end": t_end.isoformat(),
           "cooling_sequence": seq, "hot_area": area, "energetics": energy,
           "mask_pixels": int(mask.sum()),
           "t_peak_before_flux_peak_min": float(
               (pd.Timestamp(energy["time_t_max"]) - pd.Timestamp(sel["peak"])).total_seconds() / 60),
           "em_peak_after_flux_peak_min": float(
               (pd.Timestamp(energy["time_em_max"]) - pd.Timestamp(sel["peak"])).total_seconds() / 60)}
    products = {"thermal": thermal, "curves": curves, "mask": mask, "episode_end": t_end}
    return products, res


def _flare_region(loc: dict, sel: dict, noaa: list) -> dict | None:
    """The NOAA region the flare belongs to.

    On the disc: the region whose reported position is nearest. Above the limb (the
    flare loops stand out against dark sky), the region then at that limb: its
    Stonyhurst longitude at the flare time within 20 degrees of 90 degrees, and its
    latitude within 10 degrees of the flare's projected latitude.
    """
    from astropy.time import Time
    from sunpy.coordinates import sun

    t = Time(sel["peak"].replace("+00:00", ""))
    l0 = sun.L0(t).to_value(u.deg)
    # the position angle of the flare, read as a latitude (exact for a source in the
    # plane of the sky, which a limb flare is)
    lat_proj = float(np.degrees(np.arctan2(loc["hpc_y_arcsec"], abs(loc["hpc_x_arcsec"]))))
    day = t.datetime.date().isoformat()
    best, best_score = None, 1e9
    for r in noaa:
        if r["observed_date"][:10] != day:
            continue
        stony = (r["carrington_longitude"] - l0 + 180) % 360 - 180
        if loc["radius_rsun"] >= 1.0:
            limb_side = 90.0 if loc["hpc_x_arcsec"] > 0 else -90.0
            if abs(stony - limb_side) > 20 or abs(r["latitude"] - lat_proj) > 10:
                continue
            score = abs(stony - limb_side) + abs(r["latitude"] - lat_proj)
        else:
            score = np.hypot(stony - loc.get("hgs_lon_deg", 0), r["latitude"] - loc.get("hgs_lat_deg", 0))
            if score > 10:
                continue
        if score < best_score:
            best, best_score = r, score
    if best is None:
        return None
    return {"region": int(best["region"]), "latitude": best["latitude"],
            "stonyhurst_lon_at_flare": float((best["carrington_longitude"] - l0 + 180) % 360 - 180)}


# --- Everything ---------------------------------------------------------------------------------


def analyse(paths: C.EpochPaths, *, quiet: bool = False) -> tuple[dict, dict]:
    manifest = fetch_data.load_manifest(paths)
    results = {"epoch": paths.epoch.isoformat(), "analysed_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    products = {}

    _log("snapshot", quiet)
    snap, results["snapshot"] = snapshot(paths, manifest)
    products.update(snap)

    _log("GOES X-rays", quiet)
    xrs, events, results["xrays"] = xray_week(paths)
    products.update(xrs=xrs, events=events)

    _log("the week's magnetograms", quiet)
    hmi_rels = [v for k, v in sorted(manifest["week"]["hmi"].items()) if v]
    hmi_week = [sdo.load_map(paths.data / r) for r in hmi_rels]
    grids, results["rotation"] = differential_rotation(hmi_week, quiet)
    tracks, results["regions"] = regions(paths, manifest, hmi_week, quiet)
    products.update(grids=grids, tracks=tracks, hmi_week_times=[m.date.isot for m in hmi_week],
                    synoptic=rotation.synoptic_map(grids), epoch_paths=paths)
    noaa = rotation.noaa_regions(paths.data / manifest["goes"]["regions"])
    # each tracked NOAA region's own rotation rate against the fitted profile
    coef = results["rotation"]["fits"]["ABC"]["coef"]
    for s in results["regions"]["noaa_regions"]:
        if "omega" in s:
            s["omega_profile"] = float(sp.rotation_law(s["lat"], *coef))

    _log("the limb and the corona", quiet)
    results["limb"] = limb_physics(paths, manifest, products["maps"], quiet)

    _log("the flare of the week", quiet)
    fl_products, results["flare"] = flare_of_the_week(paths, manifest, xrs, noaa, quiet)
    products["flare"] = fl_products
    return products, results


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, float):
        return None if not np.isfinite(o) else o
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (pd.Timestamp, datetime)):
        return o.isoformat()
    return o


def save(paths: C.EpochPaths, products: dict, results: dict) -> None:
    paths.results.write_text(json.dumps(_jsonable(results), indent=1) + "\n")
    paths.cache.mkdir(parents=True, exist_ok=True)
    keep = {k: v for k, v in products.items() if k not in ("grids",)}
    with open(paths.cache / "products.pkl", "wb") as fh:
        pickle.dump(keep, fh)


def load_cached(paths: C.EpochPaths) -> tuple[dict, dict] | None:
    pkl = paths.cache / "products.pkl"
    if not (pkl.exists() and paths.results.exists()):
        return None
    with open(pkl, "rb") as fh:
        products = pickle.load(fh)
    return products, json.loads(paths.results.read_text())
