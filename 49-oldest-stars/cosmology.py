"""The oldest stars against the age of the universe: the 1990s "age crisis", restaged.

A flat universe with matter density Omega_m and a cosmological constant has age

    t0 = 2 / (3 H0 sqrt(Omega_L)) * asinh( sqrt(Omega_L / Omega_m) ),

which reduces to 2 / (3 H0) with no dark energy (Omega_m = 1). In the early
1990s, H0 ~ 70-80 km/s/Mpc and Omega_m = 1 gave 8-9 Gyr, while globular
clusters looked 14-16 Gyr old. Dark energy (which makes the universe older for
the same H0) and Hipparcos distances (which made the clusters younger) ended it.

Here the clusters' ages from this project are set against the universe's age
from D1 (supernovae + Cepheids, H0 = 73.5), D3 (the CMB, H0 = 66.7) and D8
(cosmic chronometers). Stars must also have had time to form: the oldest
clusters formed at redshift z_f > 11 (as Valcin et al. assume), which adds
Delta t = t(z_f), 0.1-0.4 Gyr, to their ages. Requiring t0 >= t_GC + Delta t
turns the cluster ages into an upper limit on H0 for each Omega_m.

Radiation is ignored in t0 (it shortens the age by ~0.003 Gyr, far below the
cluster uncertainties).
"""

from __future__ import annotations

import json

import numpy as np

import config as C

HUBBLE_TIME_GYR = 977.792221  # 1/H0 in Gyr for H0 = 1 km/s/Mpc
METAL_POOR = ("M92", "M15", "NGC 6397")


def age_flat_lcdm(h0, om):
    """Age of a flat matter + Lambda universe, Gyr (vectorised)."""
    h0, om = np.asarray(h0, float), np.asarray(om, float)
    ol = 1.0 - om
    th = HUBBLE_TIME_GYR / h0
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = 2.0 / (3.0 * np.sqrt(ol)) * np.arcsinh(np.sqrt(ol / om))
    return np.where(ol > 1e-9, th * lam, th * 2.0 / 3.0)


def age_at_redshift(z, h0, om):
    """Cosmic time at redshift z in a flat matter + Lambda universe, Gyr."""
    ol = 1.0 - om
    th = HUBBLE_TIME_GYR / h0
    a = 1.0 / (1.0 + np.asarray(z, float))
    return th * 2.0 / (3.0 * np.sqrt(ol)) * np.arcsinh(np.sqrt(ol / om) * a**1.5)


def h0_max(t_needed, om):
    """Largest H0 for which the universe is at least t_needed old (inverts t0 ~ 1/H0)."""
    return age_flat_lcdm(1.0, om) / np.asarray(t_needed, float)


def load_universe_ages() -> dict:
    out = {}
    if C.D1_RESULTS.exists():
        d1 = json.loads(C.D1_RESULTS.read_text())
        cal = d1["calibrated"]
        out["D1"] = {"label": "D1: supernovae + Cepheids", "age": cal["age_gyr"]["median"],
                     "err": 0.5 * (cal["age_gyr"]["plus"] + cal["age_gyr"]["minus"]),
                     "H0": cal["H0"]["median"], "H0_err": cal["H0"]["std"]}
    if C.D3_FIT.exists():
        d3 = json.loads(C.D3_FIT.read_text())["derived"]
        out["D3"] = {"label": "D3: Planck CMB (LCDM)", "age": d3["age"]["value"], "err": d3["age"]["error"],
                     "H0": d3["H0"]["value"], "H0_err": d3["H0"]["error"], "Om": d3["omegam"]["value"],
                     "Om_err": d3["omegam"]["error"]}
    if C.D8_RESULTS.exists():
        d8 = json.loads(C.D8_RESULTS.read_text())["flat_lcdm"]["full"]
        out["D8"] = {"label": "D8: cosmic chronometers", "age": d8["age_gyr"]["median"],
                     "err": 0.5 * (d8["age_gyr"]["plus"] + d8["age_gyr"]["minus"]),
                     "H0": d8["H0"]["median"], "H0_err": 0.5 * (d8["H0"]["plus"] + d8["H0"]["minus"])}
    return out


def age_crisis(results: dict, n: int = 200_000, seed: int = C.RANDOM_SEED) -> tuple[dict, dict]:
    rng = np.random.default_rng(seed)
    clusters = results["clusters"]
    universe = load_universe_ages()

    # age of the oldest population: the three metal-poor clusters. Their errors are
    # mostly shared (distance scale, models), so the mean keeps the typical error.
    ages = np.array([clusters[c]["combined"]["age"] for c in METAL_POOR])
    errs = np.array([clusters[c]["combined"]["total_err"] for c in METAL_POOR])
    t_gc, t_gc_err = float(np.mean(ages)), float(np.mean(errs))
    oldest = max(METAL_POOR, key=lambda c: clusters[c]["combined"]["age"])

    # formation delay: z_f uniform in [11, 30] in D3's cosmology
    h0_ref = universe.get("D3", {}).get("H0", 67.4)
    om_ref = universe.get("D3", {}).get("Om", 0.315)
    zf = rng.uniform(11, 30, n)
    delay = age_at_redshift(zf, h0_ref, om_ref)
    t_need = t_gc + t_gc_err * rng.standard_normal(n) + delay

    out = {"t_gc": t_gc, "t_gc_err": t_gc_err, "metal_poor_clusters": list(METAL_POOR),
           "oldest_cluster": oldest,
           "delay_gyr": {"median": float(np.median(delay)), "p16": float(np.percentile(delay, 16)),
                         "p84": float(np.percentile(delay, 84))},
           "t_universe_min": {"median": float(np.median(t_need)), "std": float(np.std(t_need)),
                              "p5": float(np.percentile(t_need, 5))},
           "universe": universe}

    # H0 upper limits for the matter densities of D1 and D3 (and 0.3)
    limits = {}
    for label, om, om_err in (("Omega_m = 0.30", 0.30, 0.0),
                              ("D3's Omega_m", om_ref, universe.get("D3", {}).get("Om_err", 0.0)),
                              ("D1's Omega_m", 0.332, 0.018)):
        om_s = om + om_err * rng.standard_normal(n)
        hmax = h0_max(t_need, om_s)
        limits[label] = {"Om": om, "H0_max_median": float(np.median(hmax)),
                         "H0_max_95": float(np.percentile(hmax, 95)),
                         "H0_max_84": float(np.percentile(hmax, 84))}
    out["h0_limits"] = limits

    # tension between each universe age and the clusters (+ delay)
    tensions = {}
    for key, u in universe.items():
        diff = u["age"] - out["t_universe_min"]["median"]
        sig = np.hypot(u["err"], out["t_universe_min"]["std"])
        tensions[key] = {"universe_age": u["age"], "difference": float(diff), "sigma": float(diff / sig),
                         "p_universe_older": float(np.mean(u["age"] + u["err"] * rng.standard_normal(n)
                                                           > t_need))}
    out["tensions"] = tensions

    # the 1990s version: no dark energy
    out["eds"] = {f"H0={h}": float(age_flat_lcdm(h, 1.0)) for h in (50, 67.4, 73.5, 80)}
    out["eds_h0_max"] = float(np.median(h0_max(t_need, 1.0)))

    art = {"t_need": t_need, "delay": delay}
    return out, art
