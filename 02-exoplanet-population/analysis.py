"""The measurements.

Four questions, one section each:

1. **Census** -- how the known population grew, by discovery method, and which
   single papers made the tallest spikes.
2. **Radius valley** -- where the gap between super-Earths and sub-Neptunes
   sits, how its visibility depends on how precisely radii are known, and how
   it tilts with orbital period.
3. **Composition** -- which small planets are too large to be bare rock, or
   even pure water.
4. **Habitable zone** -- which planets receive Earth-like starlight, and how
   many of those are small enough to be rocky.

Every ``analyse_*`` function returns ``(summary, artefacts)``: ``summary`` is a
plain dict of numbers (saved to ``data/results.json`` and printed by ``run.py``);
``artefacts`` holds arrays the figures need. Nothing here draws anything.

A word on what these numbers are. Every count below is a count of *discovered*
planets. Discovery is heavily biased -- towards short orbits, big planets,
bright stars, and whatever the surveys of the day could do -- so no count here
is an occurrence rate, and no ratio of counts says how common one kind of planet
is relative to another.
"""

from __future__ import annotations

import json
from datetime import datetime

import numpy as np
import pandas as pd

import planet_physics as pp
from catalog import Sample, symmetric_error
from config import (
    METHOD_ORDER, MR_EXCESS_SIGMA, MR_HIGHLIGHT, MR_ILLUSTRATIVE_CORE_MASS,
    N_BOOTSTRAP, RESULTS_PATH, ROCKY_RADIUS_MAX, HZ_HIGHLIGHT, HZ_TEFF_VALID_K,
    VALLEY_KDE_BANDWIDTH_DEX, VALLEY_PEAK_RANGE_REARTH, VALLEY_PIVOT_PERIOD_DAYS,
    VALLEY_SEARCH_RANGE_REARTH, VALLEY_SLOPE_GRID, VAN_EYLEN_INTERCEPT,
    VAN_EYLEN_SLOPE, ZENG_TEMPERATURE_K, ZENG_WATER_BY_TEMPERATURE,
)
from fetch_data import load_model, load_model_file


def _percentiles(values, q=(16, 50, 84)) -> list[float]:
    return [float(v) for v in np.percentile(np.asarray(values, dtype=float), q)]


# ===========================================================================
# 1. Census
# ===========================================================================

def analyse_census(sample: Sample, retrieved_utc: str) -> tuple[dict, dict]:
    """Discoveries per year by method, cumulative milestones, and the spikes."""
    df = sample.frame
    years = np.arange(int(df["disc_year"].min()), int(df["disc_year"].max()) + 1)
    counts = (pd.crosstab(df["disc_year"], df["method"])
              .reindex(index=years, columns=METHOD_ORDER, fill_value=0))
    totals = counts.sum(axis=1)
    cumulative = totals.cumsum()

    milestones = {}
    for n in (1000, 2000, 3000, 4000, 5000, 6000):
        if cumulative.iloc[-1] >= n:
            milestones[str(n)] = int(cumulative.index[np.argmax(cumulative.values >= n)])

    # The three tallest years, and the single paper behind most of each.
    spikes = []
    for year, n in totals.nlargest(3).items():
        refs = df.loc[df["disc_year"] == year, "disc_reference"].value_counts()
        spikes.append({
            "year": int(year), "planets": int(n),
            "reference": refs.index[0], "reference_planets": int(refs.iloc[0]),
        })
    spikes.sort(key=lambda s: s["year"])

    firsts = {}
    for method in METHOD_ORDER[:4]:
        rows = df[df["method"] == method].sort_values(["disc_year", "pl_name"])
        if len(rows):
            firsts[method] = {"planet": rows.iloc[0]["pl_name"], "year": int(rows.iloc[0]["disc_year"])}
    earliest = df.sort_values(["disc_year", "pl_name"]).iloc[0]

    other = df[df["method"] == "Other"]["discoverymethod"].value_counts()
    retrieved = datetime.fromisoformat(retrieved_utc)

    summary = {
        "planets": int(len(df)),
        "first_year": int(years[0]),
        "last_year": int(years[-1]),
        "partial_last_year": {"year": retrieved.year, "through": retrieved.date().isoformat()},
        "by_method": {m: int(counts[m].sum()) for m in METHOD_ORDER},
        "other_methods": {k: int(v) for k, v in other.items()},
        "cumulative_milestones": milestones,
        "spikes": spikes,
        "earliest": {"planet": earliest["pl_name"], "year": int(earliest["disc_year"]),
                     "method": earliest["discoverymethod"]},
        "first_by_method": firsts,
        "controversial_included": int(df["controversial"].sum()),
        "per_year": {int(y): {m: int(counts.loc[y, m]) for m in METHOD_ORDER} for y in years},
    }
    return summary, {"counts": counts, "totals": totals, "cumulative": cumulative}


def analyse_solar_analogs(sample: Sample, factor: float = 2.0) -> dict:
    """How many known planets resemble each Solar System planet?

    "Resemble" means within ``factor`` of its mass and of its orbital distance.
    Minimum masses count as masses here, which can only flatter the RV planets.
    The answer is a measure of detectability, not of how common such planets are.
    """
    df = sample.frame
    log_a = np.log10(df["pl_orbsmax"].to_numpy(float))
    log_m = np.log10(df["pl_bmasse"].to_numpy(float))
    tol = np.log10(factor)
    counts = {}
    for name, p in pp.SOLAR_SYSTEM.items():
        near = (np.abs(log_a - np.log10(p["a"])) < tol) & (np.abs(log_m - np.log10(p["mass"])) < tol)
        rows = df[near]
        counts[name] = {
            "planets": int(near.sum()),
            "by_method": {k: int(v) for k, v in rows["method"].value_counts().items()},
            # Only listed for the rare ones; Jupiter has dozens.
            "names": rows["pl_name"].tolist() if near.sum() <= 5 else [],
            "host_masses_msun": rows["st_mass"].tolist() if near.sum() <= 5 else [],
        }
    return {"factor": factor, "sample_size": int(len(df)), "analogs": counts}


# ===========================================================================
# 2. The radius valley
# ===========================================================================

#: Fine grid in log10(radius) on which densities are evaluated. It extends well
#: beyond the analysis window so that no planet near the window's edges is lost
#: from the density estimate.
LOG_DX = 0.002
LOG_GRID = np.arange(np.log10(0.5), np.log10(8.0) + LOG_DX / 2, LOG_DX)
_LOG_EDGES = np.concatenate([LOG_GRID - LOG_DX / 2, [LOG_GRID[-1] + LOG_DX / 2]])


def kde_log_radius(log_r, bandwidth: float = VALLEY_KDE_BANDWIDTH_DEX) -> np.ndarray:
    """Gaussian kernel density of log10(radius) on ``LOG_GRID``, per dex.

    A smoothed histogram: each planet contributes a small Gaussian bump instead
    of a hard-edged bin count, so the result does not depend on where bin edges
    happen to fall -- which matters when the feature being measured is only a
    few bins wide.

    It is computed by histogramming onto the fine grid (0.002 dex, 15 times
    finer than the default kernel) and convolving with the kernel. At that
    resolution the result is indistinguishable from summing a Gaussian per
    planet, and fast enough to repeat thousands of times for bootstrapping.
    """
    log_r = np.asarray(log_r, dtype=float)
    counts, _ = np.histogram(log_r, _LOG_EDGES)
    half = int(np.ceil(5 * bandwidth / LOG_DX))
    x = np.arange(-half, half + 1) * LOG_DX
    kernel = np.exp(-0.5 * (x / bandwidth) ** 2)
    kernel /= kernel.sum() * LOG_DX
    return np.convolve(counts, kernel, mode="same") / max(len(log_r), 1)


def find_valley(log_r, bandwidth: float = VALLEY_KDE_BANDWIDTH_DEX,
                window=VALLEY_SEARCH_RANGE_REARTH, peak_range=VALLEY_PEAK_RANGE_REARTH) -> dict:
    """Locate the valley and the two peaks either side of it.

    The valley is the lowest point of the density inside ``window``; each peak
    is the highest point between the valley and the ``peak_range`` limit on its
    side. ``contrast`` is the density at the valley divided by the mean of the
    two peak densities: 1 means no dip at all, lower means a deeper valley.

    A minimum that lands exactly on the edge of the search window is not a
    valley -- the density was still falling -- and is flagged ``on_edge``.
    """
    density = kde_log_radius(log_r, bandwidth)
    lo, hi = np.log10(window[0]), np.log10(window[1])
    inside = np.flatnonzero((LOG_GRID >= lo) & (LOG_GRID <= hi))
    i_valley = inside[np.argmin(density[inside])]

    below = np.flatnonzero((LOG_GRID >= np.log10(peak_range[0])) & (LOG_GRID < LOG_GRID[i_valley]))
    above = np.flatnonzero((LOG_GRID > LOG_GRID[i_valley]) & (LOG_GRID <= np.log10(peak_range[1])))
    i_low = below[np.argmax(density[below])]
    i_high = above[np.argmax(density[above])]

    return {
        "valley": float(10 ** LOG_GRID[i_valley]),
        "peak_low": float(10 ** LOG_GRID[i_low]),
        "peak_high": float(10 ** LOG_GRID[i_high]),
        "contrast": float(density[i_valley] / np.mean([density[i_low], density[i_high]])),
        "on_edge": bool(i_valley in (inside[0], inside[-1])),
    }


def fit_valley_tilt(log_r, log_p_rel, slopes) -> tuple[float, np.ndarray]:
    """The slope that makes the valley deepest.

    If the valley's radius changes with period as ``R ~ P**m``, then projecting
    every planet along that slope -- ``log R' = log R - m log(P / P_pivot)`` --
    lines the valley up across all periods, and the projected histogram shows
    the deepest dip. Projecting along the wrong slope smears it. So the slope is
    measured by scanning trial values and keeping the one with the lowest
    valley contrast.
    """
    contrasts = np.array([find_valley(log_r - m * log_p_rel)["contrast"] for m in slopes])
    return float(slopes[int(np.argmin(contrasts))]), contrasts


def _tier_label(lo: float, hi: float) -> str:
    return f"{lo * 100:.0f}–{hi * 100:.0f}%"


def analyse_valley(precise: Sample, all_precisions: Sample, rng) -> tuple[dict, dict]:
    df = precise.frame
    log_r = np.log10(df["pl_rade"].to_numpy())
    log_p = np.log10(df["pl_orbper"].to_numpy())
    log_p_rel = log_p - np.log10(VALLEY_PIVOT_PERIOD_DAYS)
    n = len(log_r)

    # --- location, and its bootstrap uncertainty ------------------------------
    main = find_valley(log_r)
    boot = np.array([find_valley(log_r[rng.integers(0, n, n)])["valley"]
                     for _ in range(N_BOOTSTRAP)])
    bandwidths = (0.02, 0.025, 0.03, 0.035, 0.04)
    sensitivity = {f"{bw:.3f}": find_valley(log_r, bw)["valley"] for bw in bandwidths}

    # --- tilt with orbital period -------------------------------------------------
    slopes = np.round(np.arange(VALLEY_SLOPE_GRID[0], VALLEY_SLOPE_GRID[1] + 1e-9,
                                VALLEY_SLOPE_GRID[2]), 3)
    best_slope, contrast_curve = fit_valley_tilt(log_r, log_p_rel, slopes)
    pivot_valley = find_valley(log_r - best_slope * log_p_rel)["valley"]
    # log R = m log P + a, with the valley at the pivot period fixing a.
    intercept = float(np.log10(pivot_valley) - best_slope * np.log10(VALLEY_PIVOT_PERIOD_DAYS))

    boot_slope, boot_intercept = [], []
    for _ in range(N_BOOTSTRAP):
        idx = rng.integers(0, n, n)
        m, _ = fit_valley_tilt(log_r[idx], log_p_rel[idx], slopes)
        rv = find_valley(log_r[idx] - m * log_p_rel[idx])["valley"]
        boot_slope.append(m)
        boot_intercept.append(np.log10(rv) - m * np.log10(VALLEY_PIVOT_PERIOD_DAYS))
    boot_slope = np.array(boot_slope)

    # --- does precision decide whether the valley is visible? ---------------------
    base = all_precisions.frame
    edges = np.concatenate([[0.0], base["radius_relerr"].quantile([1 / 3, 2 / 3]).to_numpy(), [np.inf]])
    tiers, tier_curves = [], []
    for k, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        sel = base[(base["radius_relerr"] > lo) & (base["radius_relerr"] <= hi)] if k else \
            base[base["radius_relerr"] <= hi]
        lr = np.log10(sel["pl_rade"].to_numpy())
        window = sel["pl_rade"].between(0.8, 4.5)
        v = find_valley(lr)
        density = kde_log_radius(lr)
        span = (LOG_GRID >= np.log10(0.8)) & (LOG_GRID <= np.log10(4.5))
        density = density / (density[span].sum() * LOG_DX)  # unit area over 0.8-4.5 R_e
        refs = sel.loc[window, "radius_reference"].value_counts(normalize=True).head(2)
        tiers.append({
            "label": _tier_label(sel["radius_relerr"].min(), sel["radius_relerr"].max())
            if k < 2 else f"> {lo * 100:.0f}%",
            "relerr_range": [float(sel["radius_relerr"].min()), float(sel["radius_relerr"].max())],
            "planets": int(len(sel)),
            "planets_0p8_to_4p5": int(window.sum()),
            "median_relerr_0p8_to_4p5": float(sel.loc[window, "radius_relerr"].median()),
            "median_radius_0p8_to_4p5": float(sel.loc[window, "pl_rade"].median()),
            "contrast": v["contrast"],
            "valley": v["valley"],
            "valley_found": not v["on_edge"],
            "kepler_fraction": float((sel.loc[window, "disc_facility"] == "Kepler").mean()),
            "top_radius_references": {r: float(f) for r, f in refs.items()},
        })
        tier_curves.append(density)

    # --- SIMULATION: blur the best tier to the worst tier's quoted errors ----------
    # Adds Gaussian noise so each precise radius carries an error drawn from the
    # least-precise tier's real error distribution. If quoted errors alone
    # explained the missing valley, the blurred contrast would match that tier's.
    best_sel = base[base["radius_relerr"] <= edges[1]]
    worst_sel = base[base["radius_relerr"] > edges[2]]
    worst_err = worst_sel.loc[worst_sel["pl_rade"].between(0.8, 4.5), "radius_relerr"].clip(upper=1.0).to_numpy()
    best_lr = np.log10(best_sel["pl_rade"].to_numpy())
    best_err = best_sel["radius_relerr"].to_numpy()
    blurred = []
    for _ in range(300):
        target = rng.choice(worst_err, len(best_lr))
        extra = np.sqrt(np.clip(target**2 - best_err**2, 0, None))
        # a fractional error e is an error of e / ln(10) in log10(R)
        blurred.append(find_valley(best_lr + rng.standard_normal(len(best_lr)) * extra / np.log(10))["contrast"])

    # --- the published line, for comparison -----------------------------------------
    median_p = float(np.median(df["pl_orbper"]))
    summary = {
        "sample_size": n,
        "median_period_days": median_p,
        "median_radius_relerr": float(df["radius_relerr"].median()),
        "bandwidth_dex": VALLEY_KDE_BANDWIDTH_DEX,
        "search_window_rearth": list(VALLEY_SEARCH_RANGE_REARTH),
        "valley_rearth": main["valley"],
        "valley_rearth_16_50_84": _percentiles(boot),
        "peak_super_earth_rearth": main["peak_low"],
        "peak_sub_neptune_rearth": main["peak_high"],
        "contrast": main["contrast"],
        "bandwidth_sensitivity": sensitivity,
        "tilt": {
            "slope": best_slope,
            "slope_16_50_84": _percentiles(boot_slope),
            "slope_2p5_97p5": _percentiles(boot_slope, (2.5, 97.5)),
            "fraction_negative": float(np.mean(boot_slope < 0)),
            "intercept": intercept,
            "intercept_16_50_84": _percentiles(boot_intercept),
            "valley_at_pivot_rearth": float(pivot_valley),
            "pivot_period_days": VALLEY_PIVOT_PERIOD_DAYS,
            "contrast_untilted": float(contrast_curve[np.argmin(np.abs(slopes))]),
            "contrast_best": float(contrast_curve.min()),
        },
        "van_eylen_2018": {
            "slope": VAN_EYLEN_SLOPE, "intercept": VAN_EYLEN_INTERCEPT,
            "valley_at_pivot_rearth": float(pp.valley_radius(VALLEY_PIVOT_PERIOD_DAYS, VAN_EYLEN_SLOPE, VAN_EYLEN_INTERCEPT)),
            "valley_at_median_period_rearth": float(pp.valley_radius(median_p, VAN_EYLEN_SLOPE, VAN_EYLEN_INTERCEPT)),
        },
        "precision_tiers": tiers,
        "blur_simulation": {
            "note": "SIMULATION on real radii: most precise third blurred to the least precise third's quoted errors",
            "contrast_mean": float(np.mean(blurred)),
            "contrast_std": float(np.std(blurred)),
            "realisations": len(blurred),
        },
    }
    artefacts = {
        "log_grid": LOG_GRID,
        "density": kde_log_radius(log_r),
        "tier_curves": tier_curves,
        "slopes": slopes,
        "contrast_curve": contrast_curve,
        "boot_slopes": boot_slope,
        "boot_intercepts": np.array(boot_intercept),
    }
    return summary, artefacts


# ===========================================================================
# 3. Composition: what the small planets could be made of
# ===========================================================================

class ModelCurve:
    """A tabulated mass-radius model, interpolated in log-log space."""

    def __init__(self, mass, radius):
        mass, radius = np.asarray(mass, float), np.asarray(radius, float)
        # Hot, light water worlds have extended steam envelopes that make the
        # tabulated radius *fall* with mass at first; interpolation needs the
        # monotonic branch above the smallest radius.
        start = int(np.argmin(radius))
        self.log_m = np.log10(mass[start:])
        self.log_r = np.log10(radius[start:])

    @property
    def mass_range(self):
        return 10 ** self.log_m[0], 10 ** self.log_m[-1]

    def radius(self, mass):
        return 10 ** np.interp(np.log10(mass), self.log_m, self.log_r, left=np.nan, right=np.nan)

    def slope(self, mass):
        """Local d(ln R)/d(ln M): how much a mass error moves the model radius."""
        lm = np.log10(mass)
        h = 0.02
        up = np.interp(lm + h, self.log_m, self.log_r, left=np.nan, right=np.nan)
        dn = np.interp(lm - h, self.log_m, self.log_r, left=np.nan, right=np.nan)
        return (up - dn) / (2 * h)


def _excess_sigma(df: pd.DataFrame, curve_radius: np.ndarray, curve_slope: np.ndarray) -> np.ndarray:
    """How many standard errors a planet sits above a model curve.

    The error combines the radius error with the mass error mapped through the
    curve: a 25% mass error moves a rocky model radius by ~7%, because radius
    grows only as mass^0.27 for rock.
    """
    r = df["pl_rade"].to_numpy(float)
    m = df["pl_bmasse"].to_numpy(float)
    sig_r = symmetric_error(r, df["radius_relerr"])
    sig_model = curve_radius * curve_slope * (symmetric_error(m, df["mass_relerr"]) / m)
    return (r - curve_radius) / np.hypot(sig_r, sig_model)


def analyse_composition(sample: Sample, valley_summary: dict) -> tuple[dict, dict]:
    df = sample.frame.copy()
    rock = ModelCurve(*load_model("rock"))
    df["rock_radius"] = rock.radius(df["pl_bmasse"])
    df["excess_rock"] = _excess_sigma(df, df["rock_radius"].to_numpy(), rock.slope(df["pl_bmasse"].to_numpy()))

    # Pure water at the tabulated temperature nearest each planet's own.
    temps = np.array(sorted(ZENG_WATER_BY_TEMPERATURE))
    teq = df["pl_eqt"].fillna(ZENG_TEMPERATURE_K).to_numpy(float)
    nearest = temps[np.argmin(np.abs(teq[:, None] - temps[None, :]), axis=1)]
    df["water_temperature"] = nearest
    df["water_radius"] = np.nan
    df["excess_water"] = np.nan
    for t in temps:
        sel = nearest == t
        if not sel.any():
            continue
        curve = ModelCurve(*load_model_file(ZENG_WATER_BY_TEMPERATURE[t]))
        m = df.loc[sel, "pl_bmasse"].to_numpy(float)
        df.loc[sel, "water_radius"] = curve.radius(m)
        df.loc[sel, "excess_water"] = _excess_sigma(df.loc[sel], curve.radius(m), curve.slope(m))

    # Which side of the radius valley? Use our own fitted, tilted valley line.
    tilt = valley_summary["tilt"]
    df["valley_radius"] = pp.valley_radius(df["pl_orbper"], tilt["slope"], tilt["intercept"])
    small = df[df["pl_rade"] < 4.0]
    below = small[small["pl_rade"] < small["valley_radius"]]
    above = small[small["pl_rade"] >= small["valley_radius"]]

    def frac(excess: pd.Series, test) -> float:
        """Fraction passing ``test``, among planets inside the model's mass range."""
        known = excess.dropna()
        return float(test(known).mean()) if len(known) else float("nan")

    core = MR_ILLUSTRATIVE_CORE_MASS
    bare = float(ModelCurve(*load_model("earthlike")).radius(core))
    h2_m, h2_r = load_model("h2_1pct")
    enveloped = float(10 ** np.interp(np.log10(core), np.log10(h2_m), np.log10(h2_r)))

    highlights = {}
    for name in MR_HIGHLIGHT:
        row = df[df["pl_name"] == name]
        if len(row):
            r = row.iloc[0]
            highlights[name] = {
                "mass": float(r["pl_bmasse"]), "radius": float(r["pl_rade"]),
                "density_gcc": float(pp.bulk_density_gcc(r["pl_bmasse"], r["pl_rade"])),
                "excess_rock_sigma": float(r["excess_rock"]),
            }

    summary = {
        "sample_size": int(len(df)),
        "small_planets": int(len(small)),
        "below_valley": int(len(below)),
        "above_valley": int(len(above)),
        "excess_threshold_sigma": MR_EXCESS_SIGMA,
        "below_consistent_with_rock": frac(below["excess_rock"], lambda z: z <= MR_EXCESS_SIGMA),
        "above_too_big_for_rock": frac(above["excess_rock"], lambda z: z > MR_EXCESS_SIGMA),
        "above_too_big_for_water": frac(above["excess_water"], lambda z: z > MR_EXCESS_SIGMA),
        "outside_model_mass_range": int(small["excess_rock"].isna().sum()),
        "median_teq_small_k": float(small["pl_eqt"].median()),
        "envelope_illustration": {
            "core_mass_mearth": core, "bare_rearth": bare,
            "with_1pct_h2_rearth": enveloped, "temperature_k": ZENG_TEMPERATURE_K,
        },
        "highlights": highlights,
    }
    return summary, {"frame": df}


# ===========================================================================
# 4. The habitable zone
# ===========================================================================

def analyse_habitable_zone(sample: Sample) -> tuple[dict, dict]:
    df = sample.frame.copy()
    df["hz"] = pp.hz_zone(df["pl_insol"], df["st_teff"])
    rocky = df["radius_measured"] & (df["pl_rade"] <= ROCKY_RADIUS_MAX)
    larger = df["radius_measured"] & (df["pl_rade"] > ROCKY_RADIUS_MAX)
    df["size_class"] = np.select([rocky, larger], ["rocky-size", "larger"], "no radius")

    in_hz = df[df["hz"] != ""]
    cons = df[df["hz"] == "conservative"]
    small = in_hz[in_hz["size_class"] == "rocky-size"].sort_values("pl_insol", ascending=False)

    highlights = {}
    for name in HZ_HIGHLIGHT:
        row = df[df["pl_name"] == name]
        if len(row):
            r = row.iloc[0]
            highlights[name] = {
                "insolation": float(r["pl_insol"]), "teff": float(r["st_teff"]),
                "radius": float(r["pl_rade"]) if r["radius_measured"] else None,
                "zone": r["hz"] or "outside", "size_class": r["size_class"],
            }

    summary = {
        "sample_size": int(len(df)),
        "in_optimistic_hz": int(len(in_hz)),
        "in_conservative_hz": int(len(cons)),
        "rocky_size_in_hz": int(len(small)),
        "rocky_size_in_conservative_hz": int((cons["size_class"] == "rocky-size").sum()),
        "no_radius_in_hz": int((in_hz["size_class"] == "no radius").sum()),
        "m_dwarf_fraction_hz": float(in_hz["m_dwarf_host"].mean()),
        "m_dwarf_fraction_all": float(df["m_dwarf_host"].mean()),
        "m_dwarf_fraction_rocky_hz": float(small["m_dwarf_host"].mean()) if len(small) else None,
        "rocky_size_threshold_rearth": ROCKY_RADIUS_MAX,
        "hosts_below_fit_range": sorted(in_hz.loc[in_hz["st_teff"] < HZ_TEFF_VALID_K[0], "hostname"].unique().tolist()),
        "rocky_size_planets": [
            {"planet": r["pl_name"], "insolation": float(r["pl_insol"]), "teff": float(r["st_teff"]),
             "radius": float(r["pl_rade"]), "zone": r["hz"]}
            for _, r in small.iterrows()
        ],
        "highlights": highlights,
    }
    return summary, {"frame": df}


# ===========================================================================
# Saving
# ===========================================================================

def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def save_results(results: dict) -> None:
    RESULTS_PATH.write_text(json.dumps(_jsonable(results), indent=2, ensure_ascii=False))
