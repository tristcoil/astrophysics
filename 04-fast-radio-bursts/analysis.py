"""The analyses: whether the delay really goes as 1/frequency^2, where FRBs are,
how far their DMs overshoot the Milky Way, how repeaters differ from one-offs,
and what DM against redshift says about the universe's baryons.

Each ``*_analysis`` function returns a JSON-able summary (the numbers quoted in
the README) and the arrays the figures need. ``run_all`` ties them together.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy import stats

import catalogs as cat
import frb_physics as fp
import macquart_fit as mf
from config import (
    CHIME_DEC_LIMIT_DEG,
    DM_HALO,
    DM_HALO_MAX,
    DM_HALO_RANGE,
    F_DIFFUSE_FIDUCIAL,
    FIT_GRID,
    KNOWN_HOST_OUTLIERS,
    MACQUART_2020_LOG10_F,
    MANIFEST_PATH,
    RESULTS_PATH,
)

ARTEFACTS_PATH = RESULTS_PATH.with_name("artefacts.pkl")

RNG_SEED = 20260929  # bootstrap resampling only; fixed so reruns are identical


# ---------------------------------------------------------------------------
# Small statistical helpers
# ---------------------------------------------------------------------------


def _bootstrap_median(x, n: int = 2000, seed: int = RNG_SEED) -> dict:
    """Median with a 68% bootstrap interval (resampling the real values)."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    rng = np.random.default_rng(seed)
    meds = np.median(rng.choice(x, size=(n, x.size), replace=True), axis=1)
    lo, hi = np.percentile(meds, [16, 84])
    return {"median": float(np.median(x)), "lo68": float(lo), "hi68": float(hi), "n": int(x.size)}


def _two_sample(a, b) -> dict:
    """KS and Anderson-Darling two-sample tests. scipy caps the AD p-value to
    [0.001, 0.25]; we report it as a bound when it hits the cap."""
    a = np.asarray(a, float)[np.isfinite(a)]
    b = np.asarray(b, float)[np.isfinite(b)]
    ks = stats.ks_2samp(a, b)
    ad = stats.anderson_ksamp([a, b])
    return {
        "ks_statistic": float(ks.statistic), "ks_pvalue": float(ks.pvalue),
        "ad_statistic": float(ad.statistic), "ad_pvalue_capped": float(ad.pvalue),
    }


def _rayleigh(angles_deg, m: int) -> dict:
    """Rayleigh test for the m-th harmonic of a circular distribution.

    Z = N R^2, with R the mean resultant length of the angles multiplied by m;
    under uniformity p ~ exp(-Z) (with the standard small-N correction).
    m = 1 catches a single hot spot in right ascension; m = 2 catches two -- the
    pattern a Galactic population makes, since the plane crosses a declination
    circle twice.
    """
    a = np.deg2rad(np.asarray(angles_deg, float)) * m
    n = a.size
    r = np.hypot(np.cos(a).sum(), np.sin(a).sum()) / n
    zstat = n * r * r
    p = np.exp(np.sqrt(1 + 4 * n + 4 * (n * n - (n * r) ** 2)) - (1 + 2 * n))
    return {"m": m, "n": int(n), "R": float(r), "Z": float(zstat), "pvalue": float(min(p, 1.0))}


# ---------------------------------------------------------------------------
# 0. Is the NE2001 map in Galactic coordinates?
# ---------------------------------------------------------------------------


def validate_ne2001_map(bursts: pd.DataFrame, ne: cat.Ne2001Map) -> dict:
    """Compare the map with CHIME's own per-burst NE2001 values, in both frames.

    CHIME computed NE2001 for every burst; ``dm_fitb - dm_exc_ne2001`` recovers
    it. If the map is Galactic, looking it up at (l, b) should reproduce those
    values; looking it up at (RA, Dec) instead should not.
    """
    b = bursts[bursts["dm_mw_ne2001"].notna() & bursts["gl"].notna()]
    out = {}
    for frame, lon, lat in (("galactic", b["gl"], b["gb"]), ("equatorial", b["ra"], b["dec"])):
        ratio = b["dm_mw_ne2001"].to_numpy() / ne.dm(lon.to_numpy(), lat.to_numpy())
        out[frame] = {
            "median_ratio_catalogue_over_map": float(np.median(ratio)),
            "p16": float(np.percentile(ratio, 16)), "p84": float(np.percentile(ratio, 84)),
        }
    out["n"] = int(len(b))
    out["frame_adopted"] = "galactic"
    hi = b["gb"].abs() > 60
    out["polar_median_ne2001_catalogue"] = float(b.loc[hi, "dm_mw_ne2001"].median())
    return out


# ---------------------------------------------------------------------------
# 0b. Is the delay really proportional to 1/frequency^2?
# ---------------------------------------------------------------------------


def dispersion_analysis() -> tuple[dict, dict]:
    """Measure the dispersion law directly on the two bright one-off bursts."""
    import waterfalls as wf
    from config import DISPERSION_BURSTS

    out, art = {}, {}
    for tns in DISPERSION_BURSTS:
        arr = wf.subband_arrival_times(tns)
        fit = wf.fit_dispersion_law(arr)
        art[tns] = {**arr, **fit}
        out[tns] = {k: v for k, v in fit.items() if k not in ("residuals_s",)}
    return out, art


# ---------------------------------------------------------------------------
# 1. Where on the sky?
# ---------------------------------------------------------------------------


def _galactic_latitude(ra_deg, dec_deg):
    """Galactic latitude b (signed, degrees) of equatorial coordinates, from the
    J2000 position of the north Galactic pole (RA 192.85948, Dec 27.12825).
    Numpy-only so it vectorises over millions of points."""
    ra_p, dec_p = np.deg2rad(192.85948), np.deg2rad(27.12825)
    ra, dec = np.deg2rad(ra_deg), np.deg2rad(dec_deg)
    sinb = np.sin(dec) * np.sin(dec_p) + np.cos(dec) * np.cos(dec_p) * np.cos(ra - ra_p)
    return np.rad2deg(np.arcsin(np.clip(sinb, -1, 1)))


def _scrambled_abs_b(dec_deg, n_ra: int = 720) -> np.ndarray:
    """|b| of every source's declination circle, sampled uniformly in RA.

    This is the "RA scrambling" null model: a transit telescope sweeps every
    right ascension equally over a year, so its exposure depends on declination
    only. Keeping each source's declination and spreading it uniformly in RA
    therefore gives exactly what an isotropic sky would look like through the
    same exposure -- with no model of the exposure needed at all.
    """
    ra = (np.arange(n_ra) + 0.5) * 360.0 / n_ra
    dec = np.asarray(dec_deg, float)
    return np.abs(_galactic_latitude(ra[None, :], dec[:, None])).ravel()


def sky_analysis(pop: pd.DataFrame, psr: pd.DataFrame) -> tuple[dict, dict]:
    frb = pop[pop["gl"].notna()]
    psr_vis = psr[psr["dec_deg"] >= CHIME_DEC_LIMIT_DEG]
    out, art = {}, {}

    edges = np.arange(0, 91, 5.0)
    for label, df, dec_col, gb_col in (("frb", frb, "dec", "gb"), ("pulsar_chime_sky", psr_vis, "dec_deg", "gb")):
        absb = df[gb_col].abs().to_numpy()
        null = _scrambled_abs_b(df[dec_col].to_numpy())
        cdf_null = lambda x, null_sorted=np.sort(null): np.searchsorted(null_sorted, x, side="right") / null_sorted.size
        ks = stats.kstest(absb, cdf_null)
        frac_obs = float(np.mean(absb < 10))
        frac_exp = float(np.mean(null < 10))
        out[label] = {
            "n": int(absb.size),
            "frac_abs_b_lt_10_observed": frac_obs,
            "frac_abs_b_lt_10_isotropic": frac_exp,
            "ks_vs_isotropic_statistic": float(ks.statistic),
            "ks_vs_isotropic_pvalue": float(ks.pvalue),
            "rayleigh_ra_m1": _rayleigh(df["ra" if dec_col == "dec" else "ra_deg"], 1),
            "rayleigh_ra_m2": _rayleigh(df["ra" if dec_col == "dec" else "ra_deg"], 2),
        }
        art[label] = {
            "hist_observed": np.histogram(absb, bins=edges)[0],
            "hist_isotropic": np.histogram(null, bins=edges)[0] * absb.size / null.size,
            "edges": edges,
        }
    # Pile-up toward the celestial pole is exposure, not astrophysics: CHIME sees
    # circumpolar sky on both upper and lower transits, and its beams are wider
    # in RA near the pole.
    polar = frb["dec"] > 70
    out["exposure"] = {
        "median_hours_dec_below_70": float(frb.loc[~polar, "exposure_h"].median()),
        "median_hours_dec_above_70": float(frb.loc[polar, "exposure_h"].median()),
        "fraction_sources_dec_above_70": float(polar.mean()),
        "fraction_sky_dec_above_70": float((1 - np.sin(np.deg2rad(70))) / (1 - np.sin(np.deg2rad(CHIME_DEC_LIMIT_DEG)))),
    }
    return out, art


# ---------------------------------------------------------------------------
# 2. Dispersion measures: more than the Milky Way can supply?
# ---------------------------------------------------------------------------


def dm_analysis(pop: pd.DataFrame, psr: pd.DataFrame, ne: cat.Ne2001Map) -> tuple[dict, dict]:
    frb = pop[pop["dm_fitb"].notna()].copy()
    frb["ratio"] = frb["dm_fitb"] / frb["dm_mw_ne2001"]
    frb["dm_mw_max"] = frb["dm_mw_ne2001"] + DM_HALO_MAX
    frb["dm_eg"] = frb["dm_fitb"] - frb["dm_mw_ne2001"] - DM_HALO
    frb["z_max"] = fp.z_from_mean_dm(frb["dm_eg"])

    p = psr[psr["dm"].notna() & psr["gl"].notna()].copy()
    p["dm_mw_ne2001"] = ne.dm(p["gl"].to_numpy(), p["gb"].to_numpy())
    p["ratio"] = p["dm"] / p["dm_mw_ne2001"]

    lowest = frb.nsmallest(5, "ratio")[["tns_name", "source", "dm_fitb", "dm_mw_ne2001", "ratio", "gb"]]
    # Bursts this test cannot place outside the Milky Way, and why: either a
    # sightline close to the plane, where NE2001 predicts a lot and is least
    # certain, or simply a low DM -- a nearby burst, not a Galactic one.
    below = frb[frb["dm_fitb"] <= frb["dm_mw_max"]]
    near_plane = below["gb"].abs() < 10
    below_summary = {
        "n": int(len(below)),
        "n_within_10deg_of_plane": int(near_plane.sum()),
        "n_high_latitude": int((~near_plane).sum()),
        "high_latitude_dm_range": [float(below.loc[~near_plane, "dm_fitb"].min()),
                                   float(below.loc[~near_plane, "dm_fitb"].max())]
        if (~near_plane).any() else None,
        "sources": below[["tns_name", "dm_fitb", "dm_mw_ne2001", "gb"]].to_dict(orient="records"),
    }

    def _psr_group(mask):
        g = p[mask]
        return {"n": int(len(g)), "frac_above_ne2001": float((g["ratio"] > 1).mean()),
                "frac_above_ne2001_plus_halo": float((g["dm"] > g["dm_mw_ne2001"] + DM_HALO_MAX).mean()),
                "median_ratio": float(g["ratio"].median())}

    field = ~p["magellanic"] & ~p["globular"]
    out = {
        "frb": {
            "n": int(len(frb)),
            "median_dm": float(frb["dm_fitb"].median()),
            "median_dm_mw_ne2001": float(frb["dm_mw_ne2001"].median()),
            "median_ratio": float(frb["ratio"].median()),
            "frac_ratio_gt_1": float((frb["ratio"] > 1).mean()),
            "frac_ratio_gt_2": float((frb["ratio"] > 2).mean()),
            "frac_ratio_gt_5": float((frb["ratio"] > 5).mean()),
            "frac_above_mw_plus_max_halo": float((frb["dm_fitb"] > frb["dm_mw_max"]).mean()),
            "n_below_mw_plus_max_halo": int((frb["dm_fitb"] <= frb["dm_mw_max"]).sum()),
            "below_mw_plus_max_halo": below_summary,
            "median_dm_eg": float(frb["dm_eg"].median()),
            "median_z_max": float(np.nanmedian(frb["z_max"])),
            "p90_z_max": float(np.nanpercentile(frb["z_max"], 90)),
            "max_dm": float(frb["dm_fitb"].max()),
            "max_dm_source": str(frb.loc[frb["dm_fitb"].idxmax(), "tns_name"]),
            "lowest_ratio": lowest.to_dict(orient="records"),
            "ymw16_median_ratio": float((frb["dm_fitb"] / frb["dm_mw_ymw16"]).median()),
        },
        "pulsar": {
            "milky_way_all": _psr_group(~p["magellanic"]),
            "galactic_field": _psr_group(field),
            "globular_cluster": _psr_group(p["globular"]),
            "magellanic_clouds": _psr_group(p["magellanic"]),
        },
    }
    return out, {"frb": frb, "psr": p}


# ---------------------------------------------------------------------------
# 3. Repeaters vs one-offs
# ---------------------------------------------------------------------------

MORPH_QUANTITIES = {
    "width_ms": "intrinsic width (Gaussian sigma, ms)",
    "bandwidth_mhz": "emission bandwidth (FWTM, MHz)",
    "peak_freq": "peak frequency (MHz)",
    "dm_exc_ne2001": "DM beyond the NE2001 disk (pc cm^-3)",
    "scat_ms_400": "scattering time at 400 MHz (ms)",
}


def repeater_analysis(bursts: pd.DataFrame) -> tuple[dict, dict]:
    samples = cat.morphology_samples(bursts)
    one, rep, rep_all = samples["one_off"], samples["repeater_first"], samples["repeater_all"]
    out = {"n_one_off": len(one), "n_repeater_sources": len(rep), "n_repeater_bursts": len(rep_all),
           "quantities": {}}
    for q, label in MORPH_QUANTITIES.items():
        out["quantities"][q] = {
            "label": label,
            "one_off": _bootstrap_median(one[q]),
            "repeater_first": _bootstrap_median(rep[q]),
            "repeater_all": _bootstrap_median(rep_all[q]),
            "test_one_off_vs_repeater_first": _two_sample(one[q], rep[q]),
        }
    out["band_filling_fraction"] = {
        "one_off": float(one["band_filling"].mean()),
        "repeater_first": float(rep["band_filling"].mean()),
        "repeater_all": float(rep_all["band_filling"].mean()),
    }
    out["scattering_detected_fraction"] = {
        "one_off": float((one["scat_time"] > 0).mean()),
        "repeater_first": float((rep["scat_time"] > 0).mean()),
    }
    out["multi_component_fraction"] = {
        "one_off": float((one["n_components"] > 1).mean()),
        "repeater_first": float((rep["n_components"] > 1).mean()),
    }
    # How prolific are the repeaters? (All bursts, no morphology cut.)
    counts = bursts[bursts["is_repeater"]].groupby("source").size().sort_values(ascending=False)
    out["most_active_repeaters"] = counts.head(5).to_dict()
    out["repeater_fraction_of_sources"] = float(bursts.groupby("source")["is_repeater"].first().mean())
    return out, samples


# ---------------------------------------------------------------------------
# 4. The Macquart relation
# ---------------------------------------------------------------------------


def macquart_analysis(ne: cat.Ne2001Map, *, quick: bool = False) -> tuple[dict, dict]:
    """Baseline fit plus the variants that probe its assumptions.

    Each variant changes exactly one thing relative to the baseline.
    """
    loc, flow = cat.load_localized_frbs(ne)
    names = loc["name"].to_numpy()
    z = loc["z"].to_numpy()
    dm_eg = loc["dm_eg"].to_numpy()

    base = mf.fit(names, z, dm_eg, label="baseline")
    fits = {"baseline": base}
    if not quick:
        mq_axes = mf.grid_axes({**FIT_GRID, "log10_F": MACQUART_2020_LOG10_F})
        fits["macquart2020_setup"] = mf.fit(names, z, dm_eg, label="F <= 0.32, no MW scatter (as Macquart+20)",
                                            axes=mq_axes, sigma_mw=0.0)
        fits["no_mw_scatter"] = mf.fit(names, z, dm_eg, label="no Milky Way scatter term", sigma_mw=0.0)
        for halo in DM_HALO_RANGE:
            d = dm_eg + (DM_HALO - halo)
            keep = d > 0
            fits[f"halo_{halo:.0f}"] = mf.fit(names[keep], z[keep], d[keep],
                                              label=f"Milky Way halo {halo:.0f}, not {DM_HALO:.0f}")
        keep = ~np.isin(names, KNOWN_HOST_OUTLIERS)
        fits["no_outlier"] = mf.fit(names[keep], z[keep], dm_eg[keep], label="without FRB 20190520B")
        keep = z >= 0.1
        fits["z_ge_0.1"] = mf.fit(names[keep], z[keep], dm_eg[keep], label="only z >= 0.1")

    z_grid = np.concatenate([np.linspace(0.01, 0.1, 10), np.linspace(0.12, 1.2, 28)])
    band = mf.predictive_band(base, z_grid)
    ob_planck = mf.omega_b_h70_planck()
    out = {
        "sample_cut_flow": flow,
        "n_by_survey": loc["survey"].value_counts().to_dict(),
        "z_disagreements": loc.loc[loc["z_disagree"], ["name", "z_base", "z_host"]].to_dict(orient="records"),
        "omega_b_h70_planck": ob_planck,
        "f_d_fiducial": F_DIFFUSE_FIDUCIAL,
        "fits": {k: v.to_json() for k, v in fits.items()},
    }
    return out, {"sample": loc, "fits": fits, "band": band}


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def run_all(*, quick: bool = False, verbose: bool = True) -> tuple[dict, dict]:
    def say(msg):
        if verbose:
            print(msg, flush=True)

    rows = cat.load_chime_rows()
    bursts = cat.chime_bursts(rows)
    sources = cat.chime_sources(bursts)
    pop = cat.population_sources(sources)
    psr = cat.load_pulsars()
    ne = cat.Ne2001Map()

    results = {
        "_generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "chime_cut_flow": cat.chime_cut_flow(rows, bursts, sources, pop),
        "ne2001_map_check": validate_ne2001_map(bursts, ne),
    }
    art = {"bursts": bursts, "sources": sources, "pop": pop, "psr": psr}

    say("  dispersion law on the example bursts ...")
    try:
        results["dispersion"], art["dispersion"] = dispersion_analysis()
    except FileNotFoundError:
        say("    (dynamic spectra not downloaded -- skipped)")
    say("  sky distribution ...")
    results["sky"], art["sky"] = sky_analysis(pop, psr)
    say("  dispersion measures ...")
    results["dm"], art["dm"] = dm_analysis(pop, psr, ne)
    say("  repeaters vs one-offs ...")
    results["repeaters"], art["repeaters"] = repeater_analysis(bursts)
    say("  Macquart relation fits (a few minutes) ...")
    results["macquart"], art["macquart"] = macquart_analysis(ne, quick=quick)
    if MANIFEST_PATH.exists():
        results["data_manifest"] = json.loads(MANIFEST_PATH.read_text())
    return results, art


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return o


def save_results(results: dict) -> None:
    RESULTS_PATH.write_text(json.dumps(_jsonable(results), indent=2))


def save_artefacts(art: dict) -> None:
    """Cache the figure inputs (tables, fit grids) so figures can be redrawn
    without rerunning the fits. A local convenience file, never shared."""
    import pickle

    # Only the baseline fit's posterior grid is drawn; the variants' grids
    # (~8 MB each) are summarised in results.json already.
    slim = dict(art)
    if "macquart" in art:
        slim["macquart"] = {**art["macquart"], "fits": {"baseline": art["macquart"]["fits"]["baseline"]}}
    with open(ARTEFACTS_PATH, "wb") as fh:
        pickle.dump(slim, fh, protocol=pickle.HIGHEST_PROTOCOL)


def load_cached() -> tuple[dict, dict]:
    import pickle

    if not (RESULTS_PATH.exists() and ARTEFACTS_PATH.exists()):
        raise FileNotFoundError("no cached analysis -- run `python run.py` once first")
    with open(ARTEFACTS_PATH, "rb") as fh:
        art = pickle.load(fh)
    return json.loads(RESULTS_PATH.read_text()), art
