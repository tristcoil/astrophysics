"""The measurements.

Each ``analyse_*`` function takes clean tables from ``catalog.py`` and returns
two things: a summary of plain numbers (saved to ``data/results.json`` and
quoted in the README) and the arrays the figures need. Nothing is plotted here.

Five analyses, in the order of the figures:

1. **The snapshot** -- what passed in the last month and what is due in the next two.
2. **Encounters** -- how often asteroids of each size pass inside the Moon's
   orbit, and how late they are found.
3. **Impact rates** -- a test of the d^2 law of encounter geometry, and its
   extrapolation from near misses to hits, compared with the fireball record
   and published estimates.
4. **The risk list** -- Sentry's objects on the Torino and Palermo scales,
   recomputed here from Sentry's own probabilities and energies.
5. **Apophis in 2029** -- the flyby's geometry, what it does to the orbit, and
   how well it is known.
"""

from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd

import neo_physics as nphys
from config import (
    AU_KM, EARTH_RADIUS_KM, FIREBALL_YEARS, GEO_RADIUS_KM, GPS_RADIUS_KM, HISTORY_YEARS,
    LARGE_REF_DIST_AU, LARGE_WINDOW, LD_KM, N_BOOTSTRAP, PUBLISHED_140M_RATES,
    RECENT_YEARS, RESULTS_PATH, SIZE_CLASS_LABELS, SMALL_REF_DIST_LD, SNAPSHOT_DAYS_AFTER,
    SNAPSHOT_DAYS_BEFORE, STUART_BINZEL_1KM_INTERVAL_YR, SURVEY_YEARS,
)

YEAR_S = 365.25 * 86400


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        return None if not np.isfinite(obj) else float(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, (pd.Timestamp,)):
        return obj.isoformat()
    return obj


def save_results(results: dict) -> None:
    tmp = RESULTS_PATH.with_name(RESULTS_PATH.name + ".part")
    tmp.write_text(json.dumps(_jsonable(results), indent=2, ensure_ascii=False))
    tmp.replace(RESULTS_PATH)


def _percentiles(samples) -> list[float]:
    return [float(v) for v in np.percentile(samples, [16, 50, 84])]


def describe_approach(row) -> dict:
    """The facts about one approach that a caption or table might quote."""
    return {
        "name": row["fullname"],
        "designation": row["des"],
        "time": row["time"].strftime("%Y-%m-%d %H:%M"),
        "dist_ld": float(row["dist_ld"]),
        "dist_km": float(row["dist_km"]),
        "altitude_km": float(row["dist_km"] - EARTH_RADIUS_KM),
        "diameter_m": float(row["d_best_m"]),
        "diameter_measured": bool(row["d_measured"]),
        "diameter_range_m": [float(row["d_est_lo_m"]), float(row["d_est_hi_m"])],
        "h": float(row["h"]),
        "v_rel_kms": float(row["v_rel"]),
        "first_observed": (row["first_obs"].strftime("%Y-%m-%d")
                           if pd.notna(row["first_obs"]) else None),
        "lag_days": None if pd.isna(row["lag_days"]) else int(row["lag_days"]),
    }


def _lag_counts(df: pd.DataFrame) -> dict:
    """First observed before, on the day of, or only after closest approach."""
    return {"before": int((df["lag_days"] < 0).sum()),
            "same_day": int((df["lag_days"] == 0).sum()),
            "after": int((df["lag_days"] > 0).sum())}


def _in_years(df: pd.DataFrame, years: tuple[int, int]) -> pd.DataFrame:
    return df[(df["year"] >= years[0]) & (df["year"] <= years[1])]


# ===========================================================================
# 1. The snapshot: last month and the next two
# ===========================================================================

def analyse_snapshot(ap: pd.DataFrame, now: pd.Timestamp) -> tuple[dict, dict]:
    lo = now - pd.Timedelta(days=SNAPSHOT_DAYS_BEFORE)
    hi = now + pd.Timedelta(days=SNAPSHOT_DAYS_AFTER)
    w = ap[(ap["time"] >= lo) & (ap["time"] <= hi)].copy()
    w["past"] = w["time"] < now
    past, future = w[w["past"]], w[~w["past"]]

    # Future approaches were all found in advance, by definition. How far?
    lead_years = -future["lag_days"] / 365.25
    summary = {
        "now_utc": now.isoformat(timespec="minutes"),
        "window": [lo.strftime("%Y-%m-%d"), hi.strftime("%Y-%m-%d")],
        "days_before": SNAPSHOT_DAYS_BEFORE, "days_after": SNAPSHOT_DAYS_AFTER,
        "past": {
            "count": len(past), "per_day": len(past) / SNAPSHOT_DAYS_BEFORE,
            "within_1ld": int((past["dist_ld"] < 1).sum()),
            "first_observed": _lag_counts(past),
            "first_observed_in_window": int((past["first_obs"] >= lo.normalize()).sum()),
            "median_diameter_m": float(past["d_best_m"].median()),
        },
        "future": {
            "count": len(future), "per_day": len(future) / SNAPSHOT_DAYS_AFTER,
            "within_1ld": int((future["dist_ld"] < 1).sum()),
            "found_over_a_year_ahead": int((lead_years > 1).sum()),
            "median_diameter_m": float(future["d_best_m"].median()),
        },
        "rate_ratio_past_to_future": (len(past) / SNAPSHOT_DAYS_BEFORE)
                                     / max(len(future) / SNAPSHOT_DAYS_AFTER, 1e-9),
        "largest": describe_approach(w.loc[w["d_best_m"].idxmax()]),
        "closest": describe_approach(w.loc[w["dist_ld"].idxmin()]),
        "closest_future": (describe_approach(future.loc[future["dist_ld"].idxmin()])
                           if len(future) else None),
        "largest_future": (describe_approach(future.loc[future["d_best_m"].idxmax()])
                           if len(future) else None),
    }
    return summary, {"frame": w, "now": now, "lo": lo, "hi": hi}


# ===========================================================================
# 2. Encounters: who passes inside the Moon's orbit, and when we notice
# ===========================================================================

def analyse_encounters(ap: pd.DataFrame, now: pd.Timestamp) -> tuple[dict, dict]:
    y0, y1 = SURVEY_YEARS
    n_years = y1 - y0 + 1
    survey = _in_years(ap, SURVEY_YEARS)
    inside = survey[survey["dist_ld"] < 1]

    per_class = {}
    for label in SIZE_CLASS_LABELS:
        n_in = int((inside["size_class"] == label).sum())
        n_all = int((survey["size_class"] == label).sum())
        per_class[label] = {"inside_1ld": n_in, "inside_1ld_per_year": n_in / n_years,
                            "within_005au_per_year": n_all / n_years}

    past = ap[ap["time"] < now]
    closest = past.nsmallest(10, "dist_km")

    h0, h1 = HISTORY_YEARS
    history = _in_years(ap, HISTORY_YEARS)
    history = history[history["dist_ld"] < 1]
    by_year = pd.crosstab(history["year"], np.sign(history["lag_days"])).reindex(
        index=range(h0, h1 + 1), columns=[-1, 0, 1], fill_value=0)
    by_year.columns = ["before", "same_day", "after"]

    recent = _in_years(ap, RECENT_YEARS)
    recent = recent[recent["dist_ld"] < 1]
    lag_q = recent["lag_days"].quantile([0.1, 0.25, 0.5, 0.75, 0.9])

    summary = {
        "survey_years": list(SURVEY_YEARS),
        "inside_1ld_total": len(inside),
        "inside_1ld_per_year": len(inside) / n_years,
        "inside_geo_total": int((survey["dist_km"] < GEO_RADIUS_KM).sum()),
        "within_005au_per_year": len(survey) / n_years,
        "per_class": per_class,
        "largest_inside_1ld": describe_approach(inside.loc[inside["d_best_m"].idxmax()]),
        "closest_known_passes": [describe_approach(r) for _, r in closest.iterrows()],
        "history_years": [h0, h1],
        "inside_1ld_by_year": {int(y): row.to_dict() for y, row in by_year.iterrows()},
        "recent_years": list(RECENT_YEARS),
        "recent_inside_1ld": len(recent),
        "recent_first_observed": _lag_counts(recent),
        "recent_fraction_after": float((recent["lag_days"] > 0).mean()),
        "recent_fraction_same_day": float((recent["lag_days"] == 0).mean()),
        "recent_lag_quantiles_days": {f"{q:.2f}": float(v) for q, v in lag_q.items()},
        "recent_within_3_days": float((recent["lag_days"].abs() <= 3).mean()),
        "recent_median_diameter_m": float(recent["d_best_m"].median()),
    }
    return summary, {"survey": survey, "by_year": by_year, "recent": recent}


# ===========================================================================
# 3. From near misses to impact rates
# ===========================================================================

def distance_index(dist, d_max) -> float:
    """Maximum-likelihood exponent k of N(<d) proportional to d^k, from every
    distance below ``d_max``.

    If N(<d) grows as d^k, then u = ln(d_max/d) is exponentially distributed
    with rate k, so k = (n - 1) / sum(u) (the n - 1 removes the small-sample
    bias). Uses every approach, with no binning.
    """
    d = np.asarray(dist, float)
    d = d[(d > 0) & (d <= d_max)]
    return (len(d) - 1) / np.log(d_max / d).sum()


def bootstrap_object_sums(des, columns: list, rng, n_boot: int = N_BOOTSTRAP,
                          chunk: int = 200) -> list[np.ndarray]:
    """Bootstrap sums of per-approach quantities, resampling *asteroids*.

    One asteroid on an Earth-like orbit can contribute a dozen approaches in
    three centuries, so approaches are not independent. Resampling whole
    asteroids with replacement keeps each one's approaches together. Every
    statistic used here is a ratio of such sums, so each asteroid's totals are
    computed once and the resampling is a matrix of picks.
    """
    codes, uniques = pd.factorize(pd.Series(des).to_numpy())
    n_obj = len(uniques)
    totals = [np.bincount(codes, weights=np.asarray(c, float), minlength=n_obj) for c in columns]
    out = [np.empty(n_boot) for _ in columns]
    for start in range(0, n_boot, chunk):
        stop = min(start + chunk, n_boot)
        pick = rng.integers(0, n_obj, size=(stop - start, n_obj))
        for o, t in zip(out, totals):
            o[start:stop] = t[pick].sum(axis=1)
    return out


def index_with_errors(df: pd.DataFrame, d_max: float, rng, n_boot: int = N_BOOTSTRAP) -> dict:
    u = np.log(d_max / df["dist_au"].to_numpy())
    n_boot_arr, u_boot = bootstrap_object_sums(df["des"], [np.ones_like(u), u], rng, n_boot)
    k_boot = (n_boot_arr - 1) / u_boot
    return {"k": float(distance_index(df["dist_au"], d_max)), "k_16_50_84": _percentiles(k_boot),
            "k_2p5_97p5": [float(v) for v in np.percentile(k_boot, [2.5, 97.5])],
            "approaches": len(df), "asteroids": int(df["des"].nunique())}


def weighted_rate_with_errors(df: pd.DataFrame, weights, years: float, rng,
                              n_boot: int = N_BOOTSTRAP) -> dict:
    """A rate that is a sum of per-approach weights over a time span, with a
    bootstrap-by-asteroid uncertainty."""
    weights = np.asarray(weights, float)
    rate = float(weights.sum() / years)
    (boot,) = bootstrap_object_sums(df["des"], [weights], rng, n_boot)
    lo, mid, hi = _percentiles(boot / years)
    return {"rate_per_year": rate, "rate_16_50_84": [lo, mid, hi],
            "interval_years": 1.0 / rate if rate > 0 else None,
            "interval_16_84": [1.0 / hi, 1.0 / lo] if lo > 0 else None,
            "approaches": len(df), "asteroids": int(df["des"].nunique())}


def extrapolation_weights(df: pd.DataFrame, ref_dist_km: float, target_km: float = EARTH_RADIUS_KM):
    """For each approach within ``ref_dist_km``: the fraction of the reference
    target (in the plane of aim) that a smaller target covers -- Earth itself
    by default, gravitational focusing included (Öpik's reasoning; see
    ``neo_physics.impact_weight``). Summed over approaches and divided by the
    time span, this turns counted near misses into an expected rate."""
    v = df["v_inf"].to_numpy()
    return (nphys.impact_parameter_km(target_km, v) / nphys.impact_parameter_km(ref_dist_km, v)) ** 2


def fireball_diameters_m(fb: pd.DataFrame, v_typical: float) -> np.ndarray:
    """Diameter implied by each fireball's impact energy, for a stony density,
    at its measured entry speed where the sensors recorded one."""
    v = pd.to_numeric(fb["vel"], errors="coerce").fillna(v_typical).to_numpy(float)
    e_kt = pd.to_numeric(fb["impact_kt"], errors="coerce").to_numpy(float)
    return nphys.diameter_for_energy_km(e_kt / 1e3, v) * 1e3


def analyse_impact_rates(ap: pd.DataFrame, fb: pd.DataFrame, rng,
                         apophis_dist_km: float) -> tuple[dict, dict]:
    """``apophis_dist_km``: Apophis's 2029 miss distance, the target for the
    'how often does something this big pass this close' estimate."""
    t_long = float(int(LARGE_WINDOW[1][:4]) - int(LARGE_WINDOW[0][:4]))
    ref_large_km = LARGE_REF_DIST_AU * AU_KM
    d_grid_ld = np.logspace(np.log10(0.01), np.log10(nphys.au_to_ld(LARGE_REF_DIST_AU)), 60)

    # --- (a) the d^2 test, one size class at a time --------------------------
    index, curves = {}, {}
    for label in SIZE_CLASS_LABELS:
        sub = ap[ap["size_class"] == label]
        index[label] = index_with_errors(sub, LARGE_REF_DIST_AU, rng)
        dist = np.sort(sub["dist_ld"].to_numpy())
        curves[label] = np.searchsorted(dist, d_grid_ld, side="right") / t_long

    # --- (b) large asteroids: extrapolate every approach within 0.05 au -----
    large = {}
    for name, d_min in (("1 km", 1000.0), ("340 m", 340.0), ("140 m", 140.0)):
        sub = ap[ap["d_best_m"] >= d_min]
        rate = weighted_rate_with_errors(sub, extrapolation_weights(sub, ref_large_km), t_long, rng)
        large[name] = {"min_diameter_m": d_min, **rate}
    sb_interval, sb_err = STUART_BINZEL_1KM_INTERVAL_YR
    large["1 km"]["published"] = {"source": "Stuart & Binzel 2004", "interval_years": sb_interval,
                                  "interval_err_years": sb_err}
    large["140 m"]["published"] = {"source": "compiled by Nugent et al. 2025",
                                   "rates_per_year": PUBLISHED_140M_RATES}

    # Apophis-sized asteroids passing as close as Apophis will: the same
    # extrapolation, with Apophis's 2029 distance as the target instead of Earth.
    apophis_q_km = float(apophis_dist_km)
    big = ap[ap["d_best_m"] >= 340.0]
    apophis_like = weighted_rate_with_errors(
        big, extrapolation_weights(big, ref_large_km, target_km=apophis_q_km), t_long, rng)
    apophis_like.update({"min_diameter_m": 340.0, "distance_km": apophis_q_km,
                         "published": "NASA: 'only every few thousand years on average'"})

    # --- (c) small asteroids: observed passes inside 1 LD, recent years -----
    recent = _in_years(ap, RECENT_YEARS)
    n_recent_years = RECENT_YEARS[1] - RECENT_YEARS[0] + 1
    ref_small_km = SMALL_REF_DIST_LD * LD_KM
    near = recent[recent["dist_km"] < ref_small_km]
    w_small = extrapolation_weights(near, ref_small_km)
    d_small = near["d_best_m"].to_numpy()
    d_grid_m = np.logspace(0, np.log10(300.0), 80)
    small_curve = np.array([w_small[d_small >= d].sum() for d in d_grid_m]) / n_recent_years
    small_count = np.array([(d_small >= d).sum() for d in d_grid_m])

    sensitivity = {}
    for ref_ld in (0.25, 0.5, 1.0, 2.0):
        sub = recent[recent["dist_ld"] < ref_ld]
        ww = extrapolation_weights(sub, ref_ld * LD_KM)
        dd = sub["d_best_m"].to_numpy()
        sensitivity[f"{ref_ld:g} LD"] = {f">= {d} m": float(ww[dd >= d].sum() / n_recent_years)
                                         for d in (5, 10)}

    # --- (d) the fireball record ---------------------------------------------
    fbw = fb[(fb["year"] >= FIREBALL_YEARS[0]) & (fb["year"] <= FIREBALL_YEARS[1])]
    n_fb_years = FIREBALL_YEARS[1] - FIREBALL_YEARS[0] + 1
    v_typical = float(fbw["vel"].median())
    d_fb = fireball_diameters_m(fbw, v_typical)
    fb_curve = np.array([(d_fb >= d).sum() for d in d_grid_m]) / n_fb_years
    e_grid_kt = np.logspace(-1.3, np.log10(fbw["impact_kt"].max()), 80)
    fb_energy_curve = np.array([(fbw["impact_kt"] >= e).sum() for e in e_grid_kt]) / n_fb_years
    blocks = {}
    for start in range(FIREBALL_YEARS[0], FIREBALL_YEARS[1] + 1, 5):
        stop = min(start + 4, FIREBALL_YEARS[1])
        sub = fbw[(fbw["year"] >= start) & (fbw["year"] <= stop)]
        blocks[f"{start}-{stop}"] = float((sub["impact_kt"] >= 1).sum() / (stop - start + 1))

    compare = {}
    for d in (3, 5, 10):
        ca = float(w_small[d_small >= d].sum() / n_recent_years)
        fr = float((d_fb >= d).sum() / n_fb_years)
        compare[f">= {d} m"] = {"close_approach_lower_limit": ca, "fireballs": fr,
                                "fraction_seen": ca / fr if fr > 0 else None}

    chelyabinsk = fbw.loc[fbw["impact_kt"].idxmax()]
    summary = {
        "distance_index": index,
        "large": large,
        "apophis_like_passes": apophis_like,
        "small_lower_limits": {
            "years": list(RECENT_YEARS), "ref_dist_ld": SMALL_REF_DIST_LD,
            "passes": len(near),
            "by_min_diameter": {f">= {d} m": float(w_small[d_small >= d].sum() / n_recent_years)
                                for d in (1, 3, 5, 10, 20)},
            "reference_distance_sensitivity": sensitivity,
        },
        "fireballs": {
            "years": list(FIREBALL_YEARS), "events": len(fbw),
            "per_year_at_least_kt": {f"{e:g}": float((fbw["impact_kt"] >= e).sum() / n_fb_years)
                                     for e in (0.1, 1, 10, 100)},
            "brown2002_per_year_at_least_kt": {f"{e:g}": float(nphys.brown2002_rate(e))
                                               for e in (0.1, 1, 10, 100)},
            "per_year_1kt_by_5yr_block": blocks,
            "velocity_measured_fraction": float(fbw["vel"].notna().mean()),
            "median_velocity_kms": v_typical,
            "largest": {"time": chelyabinsk["time"].strftime("%Y-%m-%d %H:%M"),
                        "impact_kt": float(chelyabinsk["impact_kt"]),
                        "lat": float(chelyabinsk["lat"]), "lon": float(chelyabinsk["lon"]),
                        "diameter_m": float(fireball_diameters_m(chelyabinsk.to_frame().T, v_typical)[0])},
        },
        "seen_vs_fireballs": compare,
        "focusing": {
            "mean_factor_large": float(np.mean(
                (nphys.capture_radius_km(ap.loc[ap["d_best_m"] >= 1000, "v_inf"]) / EARTH_RADIUS_KM) ** 2)),
            "escape_speed_kms": float(nphys.escape_speed()),
        },
    }
    art = {"d_grid_ld": d_grid_ld, "curves": curves, "t_long": t_long,
           "d_grid_m": d_grid_m, "small_curve": small_curve, "small_count": small_count,
           "fb_curve": fb_curve,
           "e_grid_kt": e_grid_kt, "fb_energy_curve": fb_energy_curve,
           "v_typical": v_typical}
    return summary, art


# ===========================================================================
# 4. The risk list
# ===========================================================================

def analyse_risk(objs: pd.DataFrame, vis: pd.DataFrame, removed: pd.DataFrame,
                 now_year: float) -> tuple[dict, dict]:
    vis = vis.copy()
    vis["years_ahead"] = vis["impact_year"] - now_year
    within = vis["years_ahead"] <= 100.0

    # Torino, recomputed from Sentry's own probability and energy
    vis["torino_ours"] = nphys.torino_scale(vis["ip"].to_numpy(), vis["energy"].to_numpy())
    rated = vis[within]
    beyond = vis[~within].groupby("fullname")["torino_ours"].max()

    # Palermo, recomputed. Sentry evidently counts at least a year to impact:
    # with that floor even an impact two days away reproduces its value.
    years = np.maximum(vis["years_ahead"], 1.0)
    vis["palermo_ours"] = nphys.palermo_scale(vis["ip"], vis["energy"], years)
    diff = (vis["palermo_ours"] - vis["ps"]).abs()
    far = vis["years_ahead"] > 2

    # One point per object: its most probable impact date
    top = vis.loc[vis.groupby("des")["ip"].idxmax()].set_index("des")
    per_obj = objs.set_index("des").join(top[["energy", "ip", "impact_year", "torino_ours"]],
                                         rsuffix="_top")
    per_obj = per_obj.rename(columns={"ip_top": "ip_max_date", "energy": "energy_mt"})
    per_obj["beyond_100yr"] = per_obj["impact_year"] - now_year > 100

    ranked = objs.sort_values("ps_cum", ascending=False)
    removed_by_year = removed["removed"].dt.year.value_counts().sort_index()
    famous = {}
    for des in ("99942", "2024 YR4"):
        hit = removed[removed["des"] == des]
        famous[des] = hit["removed"].iloc[0].strftime("%Y-%m-%d") if len(hit) else None

    now_ts = pd.Timestamp(year=int(now_year), month=1, day=1) + pd.Timedelta(days=(now_year % 1) * 365.25)
    stale = objs["last_obs"] < now_ts - pd.Timedelta(days=5 * 365.25)
    summary = {
        "objects": len(objs),
        "possible_impacts": len(vis),
        "torino_max_counts": objs["ts_max"].value_counts(dropna=False).rename(
            lambda v: "unrated" if pd.isna(v) else str(int(v))).to_dict(),
        "torino_check": {
            "rated_dates": int(len(rated)),
            "ours_all_zero": bool((rated["torino_ours"] == 0).all()),
            "agree_with_sentry": bool((rated["torino_ours"] == rated["ts"]).all()),
            "unrated_dates_beyond_100yr": int((~within).sum()),
            "unrated_match_sentry": bool(((~within) == vis["ts"].isna()).all()),
            "what_if_within_100yr": {k: int(v) for k, v in beyond.items() if v > 0},
        },
        "palermo_check": {
            "within_0p01": float((diff <= 0.01).mean()),
            "within_0p02": float((diff <= 0.02).mean()),
            "within_0p02_more_than_2yr_ahead": float((diff[far] <= 0.02).mean()),
            "max_diff_more_than_2yr_ahead": float(diff[far].max()),
            "max_diff": float(diff.max()),
        },
        "above_minus2": [
            {"name": r["fullname"], "ps_cum": float(r["ps_cum"]), "ip": float(r["ip"]),
             "diameter_m": float(r["diameter_m"]), "years": r["range"]}
            for _, r in ranked[ranked["ps_cum"] > -2].iterrows()],
        "top_palermo": [
            {"name": r["fullname"], "ps_cum": float(r["ps_cum"]), "ip": float(r["ip"]),
             "diameter_m": float(r["diameter_m"]), "years": r["range"]}
            for _, r in ranked.head(6).iterrows()],
        "most_probable": [
            {"name": r["fullname"], "ip": float(r["ip"]), "diameter_m": float(r["diameter_m"]),
             "years": r["range"], "last_observed": r["last_obs"].strftime("%Y-%m-%d")}
            for _, r in objs.nlargest(3, "ip").iterrows()],
        "median_diameter_m": float(objs["diameter_m"].median()),
        "at_least_140m": int((objs["diameter_m"] >= 140).sum()),
        "at_least_1km": int((objs["diameter_m"] >= 1000).sum()),
        "below_1mt_fraction": float((per_obj["energy_mt"] < 1).mean()),
        "not_seen_for_5_years_fraction": float(stale.mean()),
        "median_last_observed_year": float(objs["last_obs"].dt.year.median()),
        "removed_total": len(removed),
        "removed_first": removed["removed"].iloc[0].strftime("%Y-%m-%d"),
        "removed_by_year": {int(y): int(n) for y, n in removed_by_year.items()},
        "ever_listed": len(objs) + len(removed),
        "famous_removed": famous,
    }
    return summary, {"per_object": per_obj, "vis": vis, "removed": removed}


# ===========================================================================
# 5. Apophis, 13 April 2029
# ===========================================================================

def _closest_point(vec: pd.DataFrame) -> dict:
    """Minimum geocentric distance, refined between grid points by a parabola."""
    r = vec["r_km"].to_numpy()
    i = int(np.argmin(r))
    y0, y1, y2 = r[i - 1], r[i], r[i + 1]
    offset = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)  # in grid steps
    step_days = vec["jd_tdb"].iloc[i + 1] - vec["jd_tdb"].iloc[i]
    return {"index": i, "r_km": float(y1 - 0.25 * (y0 - y2) * offset),
            "jd_tdb": float(vec["jd_tdb"].iloc[i] + offset * step_days)}


def _angle_deg(a, b) -> float:
    return math.degrees(math.acos(float(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b))))


def analyse_apophis(cad: pd.DataFrame, record: dict, coarse: pd.DataFrame, fine: pd.DataFrame,
                    range_sigma: pd.DataFrame, elements: pd.DataFrame,
                    moon: pd.DataFrame) -> tuple[dict, dict]:
    from astropy.time import Time

    ca = _closest_point(fine)
    t_ca = Time(ca["jd_tdb"], format="jd", scale="tdb")
    i = ca["index"]
    pos = fine[["x", "y", "z"]].to_numpy()
    vel = fine[["vx", "vy", "vz"]].to_numpy()

    row2029 = cad.loc[(cad["time"].dt.year == 2029) & (cad["dist_au"] < 0.001)].iloc[0]
    q_cad, v_inf = float(row2029["dist_km"]), float(row2029["v_inf"])

    # How far the path is bent: velocity directions a day either side of closest approach
    near = coarse[(coarse["jd_tdb"] - ca["jd_tdb"]).abs() <= 1.0]
    bend = _angle_deg(near[["vx", "vy", "vz"]].iloc[0], near[["vx", "vy", "vz"]].iloc[-1])

    # The geostationary ring lies in the equatorial plane at 42,164 km.
    rho = np.hypot(pos[:, 0], pos[:, 1])
    to_ring = np.hypot(rho - GEO_RADIUS_KM, pos[:, 2])
    j = int(np.argmin(to_ring))
    crossing = np.where(np.diff(np.sign(pos[:, 2])) != 0)[0]

    # Formal uncertainty at closest approach (Horizons, UTC minutes)
    k = int(range_sigma["range_km"].idxmin())
    sigma3 = float(range_sigma["range_3sigma_km"].iloc[k])

    moon_i = int((moon["jd_tdb"] - ca["jd_tdb"]).abs().idxmin())
    moon_pos = moon[["x", "y", "z"]].iloc[moon_i].to_numpy()

    before = elements.iloc[0]
    after = elements.iloc[-1]

    def orbit(el):
        return {"date": el["time_tdb"].strftime("%Y-%m-%d"), "a_au": float(el["a_au"]),
                "e": float(el["e"]), "q_au": float(el["q_au"]), "Q_au": float(el["Q_au"]),
                "i_deg": float(el["i_deg"]), "period_days": float(el["period_days"]),
                "class": ("Aten" if el["a_au"] < 1 and el["Q_au"] > 0.983 else
                          "Apollo" if el["a_au"] > 1 and el["q_au"] < 1.017 else "other")}

    orb = record["orbit"]
    phys = {p["name"]: p for p in record.get("phys_par", [])}
    disc = record.get("discovery", {})
    nongrav = {p["name"]: p for p in orb.get("model_pars", [])}
    a2 = float(nongrav["A2"]["value"]) if "A2" in nongrav else None
    a2_sigma = float(nongrav["A2"]["sigma"]) if "A2" in nongrav else None
    a_osc = float(next(el["value"] for el in orb["elements"] if el["name"] == "a"))
    e_osc = float(next(el["value"] for el in orb["elements"] if el["name"] == "e"))

    widths = nphys.au_to_km(cad["dist_max_au"] - cad["dist_min_au"])  # 3-sigma span
    pre = cad[cad["time"].dt.year < 2029]
    post = cad[cad["time"].dt.year > 2029]
    last = cad.iloc[-1]

    summary = {
        "closest_approach": {
            "time_utc": t_ca.utc.strftime("%Y-%m-%d %H:%M:%S"),
            "time_tdb": t_ca.tdb.strftime("%Y-%m-%d %H:%M:%S"),
            "distance_km": ca["r_km"],
            "altitude_km": ca["r_km"] - EARTH_RADIUS_KM,
            "distance_ld": ca["r_km"] / LD_KM,
            "jpl_cad_distance_km": q_cad,
            "jpl_cad_3sigma_km": [float(nphys.au_to_km(row2029["dist_min_au"])),
                                  float(nphys.au_to_km(row2029["dist_max_au"]))],
            "horizons_3sigma_range_km": sigma3,
            "speed_kms": float(np.linalg.norm(vel[i])),
            "v_inf_kms": v_inf,
            "declination_deg": float(np.degrees(np.arcsin(pos[i, 2] / np.linalg.norm(pos[i])))),
        },
        "two_body": {
            "eccentricity": float(nphys.hyperbolic_eccentricity(q_cad, v_inf)),
            "turning_angle_deg": float(nphys.turning_angle_deg(q_cad, v_inf)),
            "perigee_speed_kms": float(nphys.speed_at_distance(q_cad, v_inf)),
            "measured_turn_deg_pm1day": bend,
        },
        "geostationary": {
            "radius_km": GEO_RADIUS_KM,
            "inside_geo_distance": bool(ca["r_km"] < GEO_RADIUS_KM),
            "min_distance_to_ring_km": float(to_ring[j]),
            "min_distance_time_tdb": fine["time_tdb"].iloc[j].strftime("%H:%M"),
            "equator_crossing_radius_km": [float(rho[c]) for c in crossing],
            "gps_radius_km": GPS_RADIUS_KM,
        },
        "moon": {"distance_km": float(np.linalg.norm(moon_pos)),
                 "distance_from_apophis_km": float(np.linalg.norm(moon_pos - pos[i]))},
        "orbit_before": orbit(before),
        "orbit_after": orbit(after),
        "record": {
            "optical_observations": int(orb["n_obs_used"]),
            "radar_delay": int(orb.get("n_del_obs_used") or 0),
            "radar_doppler": int(orb.get("n_dop_obs_used") or 0),
            "first_obs": orb["first_obs"], "last_obs": orb["last_obs"],
            "arc_days": int(float(orb["data_arc"])),
            "solution_date": orb["soln_date"],
            "discovery": disc.get("discovery"),
            "diameter_km": float(phys["diameter"]["value"]) if "diameter" in phys else None,
            "diameter_ref": phys.get("diameter", {}).get("ref"),
            "albedo": float(phys["albedo"]["value"]) if "albedo" in phys else None,
            "rotation_h": float(phys["rot_per"]["value"]) if "rot_per" in phys else None,
            "h": float(phys["H"]["value"]) if "H" in phys else None,
            "yarkovsky_a2_au_d2": a2,
            "yarkovsky_a2_sigma": a2_sigma,
            "yarkovsky_a2_significance": abs(a2 / a2_sigma) if a2 and a2_sigma else None,
            "yarkovsky_drift_m_per_year": (float(nphys.yarkovsky_drift_m_per_year(a2, a_osc, e_osc))
                                           if a2 is not None else None),
        },
        "approaches": {
            "listed": len(cad), "first": cad["time"].iloc[0].strftime("%Y-%m-%d"),
            "last": last["time"].strftime("%Y-%m-%d"),
            "max_3sigma_width_km_before_2029": float(widths[pre.index].max()),
            "max_3sigma_width_km_after_2029": float(widths[post.index].max()),
            "last_listed": {"date": last["time"].strftime("%Y-%m-%d"),
                            "dist_ld": float(last["dist_ld"]),
                            "dist_min_ld": float(last["dist_min_ld"]),
                            "dist_max_ld": float(last["dist_max_ld"]),
                            "time_sigma_days": float(last["t_sigma_days"])},
            "closest_after_2029": describe_approach_apophis(post.loc[post["dist_min_au"].idxmin()]),
        },
    }
    art = {"fine": fine, "coarse": coarse, "moon": moon, "ca": ca, "cad": cad,
           "elements": elements, "range_sigma": range_sigma}
    return summary, art


def describe_approach_apophis(row) -> dict:
    return {"date": row["time"].strftime("%Y-%m-%d"), "dist_ld": float(row["dist_ld"]),
            "dist_min_ld": float(row["dist_min_ld"]), "dist_max_ld": float(row["dist_max_ld"]),
            "dist_min_km": float(nphys.au_to_km(row["dist_min_au"])),
            "time_sigma_days": float(row["t_sigma_days"])}
