"""The analysis: load the catalogue, sort pulsars into families, measure the populations.

Everything here works from the cached catalogue file and never touches the
network. It produces two things:

* a classified table, one row per pulsar, with our own derived ages, fields and
  luminosities alongside the catalogue's measurements -- the figures draw from it;
* a results dictionary holding every number quoted in the README, saved to
  ``data/results.json`` so the prose can be checked against the code.

Each claim the figures make ("binaries cluster with the millisecond pulsars",
"old pulsars drift out of the Galactic plane") is backed by a measurement here,
with a confidence interval where the number is a proportion.
"""

from __future__ import annotations

import json
import re
import warnings

import numpy as np
import pandas as pd
from scipy.stats import binomtest

import pulsar_physics as pp
from config import (
    AGE_BIN_EDGES_YR,
    BINARY_FRACTION_BIN_EDGES_S,
    CATALOGUE_DB,
    CRAB_BIRTH_YEAR,
    CRAB_JNAME,
    DISTANCE_CAP_KPC,
    MSP_MAX_PDOT,
    MSP_MAX_PERIOD_S,
    RESULTS_PATH,
)
from fetch_data import catalogue_version

#: Families, in the order they are drawn (most numerous first, so rarer classes
#: land on top).
CLASSES = ("normal", "rrat", "radio_quiet", "msp", "xins", "magnetar")

CLASS_NAMES = {
    "normal": "ordinary radio pulsars",
    "rrat": "rotating radio transients",
    "radio_quiet": "radio-quiet (gamma/X-ray only)",
    "msp": "millisecond pulsars",
    "xins": "X-ray isolated neutron stars",
    "magnetar": "magnetars",
}

#: Rotation-powered pulsars that are not millisecond pulsars -- the comparison
#: population for the binary fraction.
NON_RECYCLED = ("normal", "rrat", "radio_quiet")


# --- Loading and classifying --------------------------------------------------------


def load_catalogue(db_path=CATALOGUE_DB) -> pd.DataFrame:
    """Parse the catalogue with psrqpy, which also derives distances, galactic
    coordinates and the Shklovskii-corrected Pdot the same way ``psrcat`` does."""
    import psrqpy

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        query = psrqpy.QueryATNF(loadfromdb=str(db_path))
    return query.catalogue


def _type_tokens(field: str) -> set[str]:
    """``'AXP,HE[mgz+13]'`` -> ``{'AXP', 'HE'}``. Uncertain tags such as
    ``'RRAT(?)'`` are kept -- nine pulsars carry one, none in a way that matters."""
    tokens = set()
    for part in field.split(","):
        part = re.sub(r"\[.*?\]", "", part).replace("(?)", "").strip()
        if part:
            tokens.add(part)
    return tokens


def _assoc_kinds(field: str) -> set[str]:
    """``'GC:47Tuc(NGC104),XRS:[bgh+06]'`` -> ``{'GC', 'XRS'}``.

    psrqpy's own ASSOC column keeps only the *first* association, which would
    miss e.g. a gamma-ray pulsar that also sits in a supernova remnant, so we
    parse the full original string.
    """
    return {part.split(":")[0].strip() for part in field.split(",") if ":" in part}


def classify(cat: pd.DataFrame) -> pd.DataFrame:
    """One row per pulsar: measurements, family, flags, and our derived quantities."""
    types = cat["TYPE_ORIG"].fillna("").map(_type_tokens)
    assoc_str = cat["ASSOC_ORIG"].fillna("")
    assoc = assoc_str.map(_assoc_kinds)

    t = pd.DataFrame(index=cat.index)
    t["jname"] = cat["JNAME"]
    t["bname"] = cat["BNAME"]
    t["P"] = cat["P0"]
    t["Pdot"] = cat["P1"]
    #: Pdot with the Shklovskii term removed, where proper motion and distance exist.
    t["Pdot_int"] = cat["P1_I"]
    t["F0"], t["F1"], t["F2"] = cat["F0"], cat["F1"], cat["F2"]
    t["pepoch_mjd"] = cat["PEPOCH"]
    t["gl"], t["gb"] = cat["GL"], cat["GB"]
    t["dist_kpc"] = cat["DIST"]
    #: The same, but with the older NE2001 electron-density model standing in for
    #: YMW16 wherever there is no independent distance -- a robustness check.
    t["dist_ne2001_kpc"] = cat["DIST_A"].fillna(cat["DIST_DM1"])
    t["pm_mas_yr"] = cat["PMTOT"]
    t["pb_days"] = cat["PB"]
    t["companion"] = cat["BINCOMP"].fillna("")
    t["survey"] = cat["SURVEY"].fillna("")
    #: The catalogue lists every survey that detected a pulsar, discoverer first.
    t["discovered_by"] = t["survey"].str.split(",").str[0]

    is_axp = types.map(lambda s: "AXP" in s)
    is_xins = types.map(lambda s: "XINS" in s)
    is_nrad = types.map(lambda s: "NRAD" in s)
    is_rrat = types.map(lambda s: "RRAT" in s)
    t["high_energy"] = types.map(lambda s: "HE" in s)

    t["is_gc"] = assoc.map(lambda s: "GC" in s)
    t["cluster"] = assoc_str.str.extract(r"GC:([^,\[(]+)")[0]
    t["is_exgal"] = assoc.map(lambda s: "EXGAL" in s)
    t["host_galaxy"] = assoc_str.str.extract(r"EXGAL:([^,\[(]+)")[0]
    t["in_snr"] = assoc.map(lambda s: "SNR" in s)
    t["is_binary"] = (cat["PB"].notna() | cat["BINARY"].fillna("").ne("")
                      | cat["BINCOMP"].fillna("").ne(""))

    fast = t["P"] < MSP_MAX_PERIOD_S
    # True where Pdot is small, unmeasured, or measured inside a cluster (where it
    # says more about the cluster's gravity than about the pulsar's age).
    gentle = ~(t["Pdot"] > MSP_MAX_PDOT) | t["is_gc"]
    is_msp = fast & gentle & ~is_axp & ~is_xins
    t["klass"] = np.select(
        [is_axp, is_xins, is_msp, is_nrad, is_rrat],
        ["magnetar", "xins", "msp", "radio_quiet", "rrat"],
        default="normal",
    )

    # Our own derived quantities -- deliberately not the catalogue's, so that the
    # two can be compared (see `consistency_checks`).
    t["tau_yr"] = pp.characteristic_age_yr(t["P"], t["Pdot"])
    t["B_gauss"] = pp.surface_field_gauss(t["P"], t["Pdot"])
    t["edot"] = pp.spin_down_luminosity(t["P"], t["Pdot"])
    t["tau_int_yr"] = pp.characteristic_age_yr(t["P"], t["Pdot_int"])

    # Height above the Galactic plane. The Sun sits ~20 pc above the plane, which
    # is negligible next to the kiloparsec distances involved.
    t["dist_ok"] = (t["dist_kpc"].notna() & (t["dist_kpc"] < DISTANCE_CAP_KPC)
                    & ~t["is_exgal"])
    t["z_pc"] = np.abs(t["dist_kpc"] * np.sin(np.deg2rad(t["gb"]))) * 1e3
    t["z_ne2001_pc"] = np.abs(t["dist_ne2001_kpc"] * np.sin(np.deg2rad(t["gb"]))) * 1e3
    t["dist_ne2001_ok"] = t["dist_ne2001_kpc"].notna() & ~t["is_exgal"]

    #: The P-Pdot sample: a period and a positive Pdot, outside globular clusters,
    #: where the cluster's own gravity corrupts Pdot (see `globular_clusters`).
    t["in_ppdot"] = t["P"].notna() & (t["Pdot"] > 0) & ~t["is_gc"]
    return t


# --- Measurements ------------------------------------------------------------------


def proportion(k: int, n: int) -> dict:
    """k of n, with a 95% Wilson score interval (sensible even when k is 0 or n)."""
    if n == 0:
        return {"k": 0, "n": 0, "fraction": None, "ci95": [None, None]}
    ci = binomtest(int(k), int(n)).proportion_ci(confidence_level=0.95, method="wilson")
    return {"k": int(k), "n": int(n), "fraction": k / n, "ci95": [ci.low, ci.high]}


def _row(t: pd.DataFrame, idx) -> dict:
    r = t.loc[idx]
    out = {"jname": r["jname"], "P_s": r["P"], "klass": r["klass"]}
    for key in ("Pdot", "tau_yr", "B_gauss", "edot", "cluster"):
        v = r[key]
        if isinstance(v, str) or (v is not None and np.isfinite(v)):
            out[key] = v
    return out


def census(t: pd.DataFrame) -> dict:
    has_p = t["P"].notna()
    pd_ = t["Pdot"]
    return {
        "n_total": len(t),
        "n_with_period": int(has_p.sum()),
        "n_with_pdot": int(pd_.notna().sum()),
        "n_pdot_positive": int((pd_ > 0).sum()),
        "n_pdot_negative": int((pd_ < 0).sum()),
        "n_pdot_zero": int((pd_ == 0).sum()),
        "n_period_no_pdot": int((has_p & pd_.isna()).sum()),
        "n_globular_cluster": int(t["is_gc"].sum()),
        "n_clusters": int(t["cluster"].nunique()),
        "n_extragalactic": int(t["is_exgal"].sum()),
        "n_binary": int(t["is_binary"].sum()),
        "n_in_ppdot": int(t["in_ppdot"].sum()),
        "by_class_all": {c: int((t["klass"] == c).sum()) for c in CLASSES},
        "by_class_ppdot": {c: int(((t["klass"] == c) & t["in_ppdot"]).sum()) for c in CLASSES},
    }


def binary_fractions(t: pd.DataFrame) -> dict:
    """The recycling evidence: how often each population has a binary companion.

    Measured on the timed sample (a Pdot means a timing solution exists, so an
    orbit would have been noticed). Among pulsars found in the last few years
    without a solution yet, "no known companion" often just means "not solved".
    """
    timed = t["in_ppdot"]
    msp = timed & (t["klass"] == "msp")
    ordinary = timed & t["klass"].isin(NON_RECYCLED)
    out = {
        "msp": proportion((msp & t["is_binary"]).sum(), msp.sum()),
        "non_recycled": proportion((ordinary & t["is_binary"]).sum(), ordinary.sum()),
    }
    # The same, including pulsars without a timing solution, to show the effect.
    msp_all = t["P"].notna() & ~t["is_gc"] & (t["klass"] == "msp")
    out["msp_including_untimed"] = proportion((msp_all & t["is_binary"]).sum(), msp_all.sum())

    # As a function of spin period, over every rotation-powered timed pulsar.
    rot = timed & t["klass"].isin(NON_RECYCLED + ("msp",))
    edges = np.asarray(BINARY_FRACTION_BIN_EDGES_S)
    bins = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = rot & (t["P"] >= lo) & (t["P"] < hi)
        bins.append({"P_lo_s": lo, "P_hi_s": hi,
                     **proportion((m & t["is_binary"]).sum(), m.sum())})
    out["by_period"] = bins
    return out


#: Companion codes used by the catalogue, and what they mean.
COMPANIONS = {
    "UL": "ultra-light / planet (black widows)",
    "He": "helium white dwarf",
    "CO": "carbon-oxygen white dwarf",
    "NS": "neutron star",
    "MS": "main-sequence star",
}


def companions(t: pd.DataFrame) -> dict:
    """Spin period by companion type, for binaries outside globular clusters.

    The recycling picture predicts an ordering: the longer the mass transfer
    lasted, the faster the pulsar was spun up. Low-mass companions (which become
    helium white dwarfs) transfer matter for up to billions of years; massive ones
    (which explode and leave a neutron star) for a few thousand.
    """
    field = t["is_binary"] & ~t["is_gc"] & t["P"].notna()
    out = {}
    for code, name in COMPANIONS.items():
        m = field & (t["companion"] == code)
        p = t.loc[m, "P"].to_numpy()
        if len(p) == 0:
            continue
        q1, med, q3 = np.percentile(p, [25, 50, 75])
        out[code] = {"name": name, "n": int(m.sum()),
                     "median_P_ms": med * 1e3, "iqr_P_ms": [q1 * 1e3, q3 * 1e3]}
    return out


def death_valley(t: pd.DataFrame) -> dict:
    """Where radio pulsars sit relative to the two death lines.

    Only pulsars that shine in radio are counted: the death line is a statement
    about when radio emission stops, so magnetars and radio-quiet pulsars are
    irrelevant to it.
    """
    radio = t["in_ppdot"] & t["klass"].isin(("normal", "rrat"))
    upper = pp.death_line_pdot(t["P"], "dipole")
    lower = pp.death_line_pdot(t["P"], "multipole")
    above = radio & (t["Pdot"] >= upper)
    beyond = radio & (t["Pdot"] < lower)
    inside = radio & ~above & ~beyond
    beyond_names = t.loc[beyond].sort_values("P", ascending=False)["jname"].tolist()
    return {
        "n_radio": int(radio.sum()),
        "above_dipole_line": proportion(above.sum(), radio.sum()),
        "in_valley": proportion(inside.sum(), radio.sum()),
        "beyond_valley": proportion(beyond.sum(), radio.sum()),
        "beyond_valley_slowest": beyond_names[:8],
    }


def msp_ages(t: pd.DataFrame) -> dict:
    """Millisecond pulsars' characteristic ages, before and after the Shklovskii
    correction, against the age of the Universe."""
    m = t["in_ppdot"] & (t["klass"] == "msp")
    tau = t.loc[m, "tau_yr"]
    has_int = m & t["Pdot_int"].notna()
    pos_int = has_int & (t["Pdot_int"] > 0)
    tau_int = t.loc[pos_int, "tau_int_yr"]
    kin_frac = 1.0 - t.loc[has_int, "Pdot_int"] / t.loc[has_int, "Pdot"]
    return {
        "n": int(m.sum()),
        "median_tau_yr": float(np.median(tau)),
        "older_than_universe": proportion((tau > pp.AGE_OF_UNIVERSE_YR).sum(), len(tau)),
        "n_with_shklovskii": int(has_int.sum()),
        "median_kinematic_fraction_of_pdot": float(np.median(kin_frac)),
        "n_kinematic_exceeds_observed": int((has_int & (t["Pdot_int"] <= 0)).sum()),
        "older_than_universe_intrinsic": proportion(
            (tau_int > pp.AGE_OF_UNIVERSE_YR).sum(), len(tau_int)),
    }


def shklovskii_example(t: pd.DataFrame, jname: str = "J0437-4715") -> dict:
    """The nearest millisecond pulsar, where the effect is largest and best measured."""
    r = t.loc[t["jname"] == jname].iloc[0]
    shk = pp.shklovskii_pdot(r["P"], r["pm_mas_yr"], r["dist_kpc"])
    return {
        "jname": jname, "P_ms": r["P"] * 1e3, "pm_mas_yr": r["pm_mas_yr"],
        "dist_kpc": r["dist_kpc"], "pdot_observed": r["Pdot"], "pdot_shklovskii": shk,
        "pdot_intrinsic": r["Pdot"] - shk, "kinematic_fraction": shk / r["Pdot"],
        "tau_observed_yr": r["tau_yr"],
        "tau_intrinsic_yr": pp.characteristic_age_yr(r["P"], r["Pdot"] - shk),
    }


def globular_clusters(t: pd.DataFrame) -> dict:
    """Pulsars inside star clusters, whose Pdot measures the cluster's gravity.

    A pulsar accelerating towards us (falling through the cluster's potential on
    the far side) has its pulses compressed; that adds a term a P / c to Pdot,
    which in a dense cluster core outweighs the tiny intrinsic Pdot of a
    millisecond pulsar -- in either direction.
    """
    gc = t["is_gc"]
    has = gc & t["Pdot"].notna()
    counts = t.loc[gc, "cluster"].value_counts()
    return {
        "n": int(gc.sum()),
        "n_clusters": int(counts.size),
        "richest": {k: int(v) for k, v in counts.head(3).items()},
        "n_with_pdot": int(has.sum()),
        "n_positive_pdot": int((has & (t["Pdot"] > 0)).sum()),
        "n_zero_pdot": int((has & (t["Pdot"] == 0)).sum()),
        "negative_pdot": proportion((has & (t["Pdot"] < 0)).sum(), has.sum()),
        "n_negative_outside_clusters": int((~gc & (t["Pdot"] < 0)).sum()),
    }


def spin_distribution(t: pd.DataFrame) -> dict:
    """The two peaks of the spin-period distribution, over every pulsar with a
    period -- including those in clusters and those with no Pdot yet."""
    has = t["P"].notna()
    msp = has & (t["klass"] == "msp")
    slow = has & (t["klass"] != "msp")
    return {
        "n": int(has.sum()),
        "n_msp": int(msp.sum()),
        "median_msp_ms": float(t.loc[msp, "P"].median() * 1e3),
        "n_other": int(slow.sum()),
        "median_other_s": float(t.loc[slow, "P"].median()),
        "n_faster_than_2ms": int((has & (t["P"] < 0.002)).sum()),
        "n_slower_than_10s": int((has & (t["P"] > 10.0)).sum()),
    }


def heights(t: pd.DataFrame) -> dict:
    """Height above the Galactic plane against characteristic age.

    Neutron stars are born where massive stars live -- in the thin disc, within
    ~50 pc of the plane -- and receive a kick of hundreds of km/s in the supernova.
    If that picture is right, older pulsars should sit progressively farther from
    the plane. Distances are DM-based, so individual heights are uncertain by
    tens of percent; medians over hundreds of pulsars are robust to that.
    """
    edges = np.asarray(AGE_BIN_EDGES_YR)

    def by_age(ok, z_col, d_col):
        rows = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = ok & (t["tau_yr"] >= lo) & (t["tau_yr"] < hi)
            z = t.loc[m, z_col].to_numpy()
            if len(z) < 5:
                continue
            q1, med, q3 = np.percentile(z, [25, 50, 75])
            rows.append({"tau_lo_yr": lo, "tau_hi_yr": hi, "n": int(len(z)),
                         "median_z_pc": med, "iqr_z_pc": [q1, q3],
                         "median_dist_kpc": float(np.median(t.loc[m, d_col]))})
        return rows

    base = t["in_ppdot"] & t["klass"].isin(NON_RECYCLED)
    rows = by_age(base & t["dist_ok"], "z_pc", "dist_kpc")
    # The 25 kpc placeholders are mostly high-latitude, middle-aged pulsars, so
    # dropping them biases heights low. NE2001 gives every one of them a distance.
    rows_ne2001 = by_age(base & t["dist_ne2001_ok"], "z_ne2001_pc", "dist_ne2001_kpc")
    placeholder = base & t["dist_kpc"].ge(DISTANCE_CAP_KPC) & ~t["is_exgal"]

    msp = t["dist_ok"] & (t["klass"] == "msp") & ~t["is_gc"] & t["P"].notna()
    zm = t.loc[msp, "z_pc"].to_numpy()
    ordinary = t["dist_ok"] & t["klass"].isin(NON_RECYCLED) & t["P"].notna()

    def lat(m):
        b = np.abs(t.loc[m, "gb"])
        return {"n": int(m.sum()), "median_abs_b_deg": float(np.median(b)),
                "within_5deg": proportion((b < 5).sum(), len(b))}

    # One survey dominates the recent millisecond-pulsar haul: FAST's Galactic
    # Plane Pulsar Snapshot survey, which only looks within 10 deg of the plane.
    gpps = t["discovered_by"] == "fast_gpps"
    field_msp = (t["klass"] == "msp") & ~t["is_gc"] & t["P"].notna()
    return {
        "by_age": rows,
        "by_age_ne2001": rows_ne2001,
        "n_placeholder_distance": int(placeholder.sum()),
        "msp_found_by_fast_gpps": proportion((field_msp & gpps).sum(), field_msp.sum()),
        "msp": {"n": int(msp.sum()), "median_z_pc": float(np.median(zm)),
                "iqr_z_pc": [float(np.percentile(zm, 25)), float(np.percentile(zm, 75))],
                "median_dist_kpc": float(np.median(t.loc[msp, "dist_kpc"]))},
        "latitude": {"non_recycled": lat(ordinary), "msp": lat(msp),
                     "msp_not_fast_gpps": lat(msp & ~gpps)},
    }


def braking_indices(t: pd.DataFrame) -> dict:
    """Braking index n = f0 f2 / f1^2 wherever the catalogue has f2.

    Pure magnetic dipole braking predicts exactly 3. The honest finding is that
    for most pulsars the measured f2 is dominated by timing noise, and only the
    youngest give physically meaningful values. Magnetars are left out: their
    spin-down rate swings by factors of several between outbursts, so f2 there
    measures the outbursts, not a braking law.
    """
    rotation_powered = t["klass"].isin(NON_RECYCLED + ("msp",))
    has = (t["F2"].notna() & t["F1"].notna() & ~t["is_gc"] & (t["F1"] < 0)
           & rotation_powered)
    n = pd.Series(pp.braking_index(t.loc[has, "F0"], t.loc[has, "F1"], t.loc[has, "F2"]),
                  index=t.index[has])
    n_young = n[t.loc[has, "tau_yr"] < 1e4]
    rows = [{"jname": t.at[i, "jname"], "tau_yr": float(t.at[i, "tau_yr"]),
             "n": float(v)} for i, v in n_young.sort_values().items()]
    return {
        "n_measured": int(has.sum()),
        "abs_n_above_10": proportion((np.abs(n) > 10).sum(), len(n)),
        "physical_1_to_4": proportion(((n > 1) & (n < 4)).sum(), len(n)),
        "young": rows,
        "values": n.to_numpy().tolist(),
        "tau_yr": t.loc[has, "tau_yr"].to_numpy().tolist(),
    }


def crab_age_check(t: pd.DataFrame) -> dict:
    """The one pulsar whose true age is written in the historical record.

    Compare its characteristic age with the time since SN 1054 (at the epoch of
    the timing solution), then turn the discrepancy around: assuming the spin-down
    law held since birth, what birth period does the true age require?
    """
    from astropy.time import Time

    r = t.loc[t["jname"] == CRAB_JNAME].iloc[0]
    epoch_year = Time(r["pepoch_mjd"], format="mjd").decimalyear
    true_age = epoch_year - CRAB_BIRTH_YEAR
    n_meas = float(pp.braking_index(r["F0"], r["F1"], r["F2"]))
    return {
        "epoch_year": epoch_year,
        "true_age_yr": true_age,
        "characteristic_age_yr": float(r["tau_yr"]),
        "overestimate_percent": (r["tau_yr"] / true_age - 1.0) * 100.0,
        "braking_index": n_meas,
        "spin_down_age_measured_n_yr": float(pp.spin_down_age_yr(r["P"], r["Pdot"], n_meas)),
        "period_now_ms": r["P"] * 1e3,
        "birth_period_ms_n3": float(pp.birth_period_s(r["P"], r["Pdot"], true_age, 3.0)) * 1e3,
        "birth_period_ms_measured_n": float(
            pp.birth_period_s(r["P"], r["Pdot"], true_age, n_meas)) * 1e3,
    }


def spin_extremes(t: pd.DataFrame) -> dict:
    """Record holders, plus what the fastest spin proves about neutron stars."""
    has_p = t["P"].notna()
    samp = t["in_ppdot"]
    # A magnetar's characteristic age is not an age (its spin-down is driven by
    # outbursts), so "youngest" is judged among rotation-powered pulsars.
    rot = samp & t["klass"].isin(NON_RECYCLED + ("msp",))
    fastest = t.loc[has_p, "P"].idxmin()
    p = t.at[fastest, "P"]
    return {
        "fastest": _row(t, fastest),
        "fastest_spin_hz": 1.0 / p,
        "fastest_min_density_g_cm3": float(pp.minimum_mean_density(p)),
        "fastest_max_radius_km_1p4": float(pp.max_radius_km(1.0 / p, 1.4)),
        "fastest_equator_fraction_c_10km": float(pp.equatorial_speed_km_s(p, 10.0)
                                                 / (pp.C_CGS / 1e5)),
        "fastest_outside_clusters": _row(t, t.loc[has_p & ~t["is_gc"], "P"].idxmin()),
        "slowest": _row(t, t.loc[has_p, "P"].idxmax()),
        "strongest_field": _row(t, t.loc[samp, "B_gauss"].idxmax()),
        "weakest_field": _row(t, t.loc[samp, "B_gauss"].idxmin()),
        "youngest_rotation_powered": _row(t, t.loc[rot, "tau_yr"].idxmin()),
        "most_powerful": _row(t, t.loc[samp, "edot"].idxmax()),
    }


def consistency_checks(t: pd.DataFrame, cat: pd.DataFrame) -> dict:
    """Our derived quantities against the ones psrqpy computes for the catalogue.

    Two independent implementations of the same formulas should agree to
    rounding. Any disagreement means one of them has a units slip.
    """
    m = t["in_ppdot"] | (t["P"].notna() & (t["Pdot"] > 0))
    out = {}
    for ours, theirs in (("tau_yr", "AGE"), ("B_gauss", "BSURF"), ("edot", "EDOT")):
        a, b = t.loc[m, ours].to_numpy(), cat.loc[m, theirs].to_numpy()
        ok = np.isfinite(a) & np.isfinite(b)
        out[theirs] = {"n": int(ok.sum()),
                       "max_rel_diff": float(np.max(np.abs(a[ok] / b[ok] - 1.0)))}
    # Shklovskii: our term against the difference psrqpy applies to get P1_I.
    # (Zero proper motion is skipped: one pulsar, J1836-0436, has its motion
    # catalogued only in ecliptic coordinates, so its equatorial total reads 0.)
    s = (t["Pdot_int"].notna() & (t["pm_mas_yr"] > 0) & t["dist_kpc"].notna())
    ours = pp.shklovskii_pdot(t.loc[s, "P"], t.loc[s, "pm_mas_yr"], t.loc[s, "dist_kpc"])
    theirs = (t.loc[s, "Pdot"] - t.loc[s, "Pdot_int"]).to_numpy()
    out["SHKLOVSKII"] = {"n": int(s.sum()),
                         "max_rel_diff": float(np.max(np.abs(ours / theirs - 1.0)))}
    return out


# --- Driver ----------------------------------------------------------------------------


def analyse(db_path=CATALOGUE_DB):
    """Run every measurement. Returns (classified table, results dictionary)."""
    cat = load_catalogue(db_path)
    t = classify(cat)
    results = {
        "catalogue_version": catalogue_version(db_path),
        "census": census(t),
        "binary_fractions": binary_fractions(t),
        "companions": companions(t),
        "death_valley": death_valley(t),
        "msp_ages": msp_ages(t),
        "shklovskii_example": shklovskii_example(t),
        "globular_clusters": globular_clusters(t),
        "spin_distribution": spin_distribution(t),
        "heights": heights(t),
        "braking_indices": braking_indices(t),
        "crab": crab_age_check(t),
        "extremes": spin_extremes(t),
        "consistency": consistency_checks(t, cat),
    }
    return t, results


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
    return obj


def save_results(results: dict, path=RESULTS_PATH) -> None:
    slim = dict(results)
    # The per-pulsar braking indices are for the figure, not the record.
    slim["braking_indices"] = {k: v for k, v in results["braking_indices"].items()
                               if k not in ("values", "tau_yr")}
    path.write_text(json.dumps(_jsonable(slim), indent=2) + "\n")
