"""Every analysis of the project, from cached data to numbers.

Each ``analyse_*`` function returns (results, artefacts): ``results`` is
JSON-serialisable and lands in data/results.json; ``artefacts`` holds the
arrays the figures need.
"""

from __future__ import annotations

import json
import re
import zlib

import numpy as np
import pandas as pd

import blackcat
import gaia_orbits
import graveyard
import literature as lit
import masses
import outbursts
from config import (
    GAIA_BH_SOURCE_IDS, KPC_CM, L_EDD_PER_MSUN, N_MASS_SAMPLES, N_POP_SAMPLES,
    NS_MAX_CAUSAL, NS_MAX_KALOGERA_BAYM, POP_NLIVE, RANDOM_SEED, RESULTS_PATH,
)
from fetch_data import CATNOIR_PATH


# ---------------------------------------------------------------------------
# 1. Weighing: mass functions and masses of the X-ray binaries
# ---------------------------------------------------------------------------

def _inclination_draws(spec, q, rng):
    """Inclination draws for a BlackCAT spec; None if the mass is unbounded."""
    n = len(q)
    if spec is None or spec[0] == "max":
        return None
    if spec[0] == "I_ecl":  # isotropic from the lower limit up to the no-eclipse limit
        i_max = masses.max_inclination_no_eclipse(q)
        lo = np.cos(np.radians(i_max))
        hi = np.cos(np.radians(spec[1]))
        return np.degrees(np.arccos(rng.uniform(lo, hi)))
    return masses.draw(spec, n, rng, lo=0.0, hi=90.0)


def blackcat_mass_samples(a: dict, n: int, rng) -> dict:
    """f(M) from (P, K2), then M1 from (q, i), for one adopted BlackCAT system."""
    porb_d = masses.draw(("N", a["porb_h"][0] / 24, a["porb_h"][1] / 24), n, rng, lo=0.0)
    k2 = masses.draw(("N", *a["k2"]) if a["k2"][1] > 0 else ("fixed", a["k2"][0]), n, rng, lo=0.0)
    fm = masses.mass_function(porb_d, k2)
    q = masses.draw(a["q_spec"], n, rng, lo=0.0)
    incl = _inclination_draws(a["i_spec"], q, rng)
    out = {"fm": fm, "q": q}
    if incl is None:  # only a minimum mass: q = 0 at the highest allowed inclination
        i_lim = a["i_spec"][1] if a["i_spec"] is not None else 90.0
        out["m1_min"] = masses.m1_from_fm(fm, 0.0, i_lim)
    else:
        out["incl"] = incl
        out["m1"] = masses.m1_from_fm(fm, q, incl)
    return out


def analyse_weighing(rng) -> tuple[dict, dict]:
    systems = blackcat.load_online()
    adopted = blackcat.adopted_inputs(systems)
    vizier = blackcat.load_vizier_a4()
    rows, draws = [], {}
    for a in adopted:
        d = blackcat_mass_samples(a, N_MASS_SAMPLES, rng)
        draws[a["name"]] = d
        fm = masses.summarise(d["fm"])
        pub = a["fm_published"]
        pull = None
        if pub is not None and "plus" in pub:
            sig = np.hypot(0.5 * (fm["plus"] + fm["minus"]), pub["plus"])
            pull = (fm["median"] - pub["value"]) / sig
        row = {
            "id": a["id"], "name": a["name"], "common_name": a["common_name"],
            "sptype": a["sptype"],
            "porb_h": a["porb_h"], "k2_kms": a["k2"],
            "k2_ref": "; ".join(k for _, k in a["k2_ref"]),
            "q": masses.spec_summary(a["q_spec"]), "q_note": a["q_note"],
            "q_ref": "; ".join(k for _, k in a["q_ref"]),
            "i": ("unconstrained" if a["i_spec"] is None else
                  f"< {a['i_spec'][1]:g}" if a["i_spec"][0] == "max" else
                  f"I({a['i_spec'][1]:g}, i_ecl)" if a["i_spec"][0] == "I_ecl" else
                  masses.spec_summary(a["i_spec"])),
            "i_note": a["i_note"], "i_ref": "; ".join(k for _, k in a["i_ref"]),
            "fm_ours": fm, "fm_published": a["fm_published_raw"],
            "fm_ref": "; ".join(k for _, k in a["fm_ref"]),
            "fm_pull": pull,
            "m1_published": a["m1_published_raw"],
            "m1_ref": "; ".join(k for _, k in a["m1_ref"]),
            "notes": a["notes"],
            "in_population": "m1" in d,
        }
        if "m1" in d:
            row["m1_ours"] = masses.summarise(d["m1"])
        else:
            row["m1_min_ours"] = masses.summarise(d["m1_min"])
        rows.append(row)

    # Cyg X-1 recomputed from its orbit (Miller-Jones et al. 2021, Table 1)
    cyg = lit.PERSISTENT_BHS[0]
    n = N_MASS_SAMPLES
    fm_cyg = masses.mass_function(cyg["porb_d"], masses.draw(cyg["k_kms"], n, rng, lo=0), cyg["e"])
    m_cyg = masses.solve_m_bh(fm_cyg, masses.draw(cyg["m_companion"], n, rng, lo=0),
                              masses.draw(cyg["incl_deg"], n, rng, lo=0, hi=90))
    cyg_check = {"fm": masses.summarise(fm_cyg), "m_bh": masses.summarise(m_cyg),
                 "published": lit.CHECKS["cygx1_mass"][:2]}

    persistent = []
    for p in lit.PERSISTENT_BHS:
        s = m_cyg if p["name"] == "Cyg X-1" else masses.draw(p["m_bh"], n, rng, lo=0)
        draws[p["name"]] = {"m1": s}
        fm = masses.mass_function(p["porb_d"], p["k_kms"][1]) if p.get("k_kms") else None
        persistent.append({"name": p["name"], "host": p["host"], "type": p["type"],
                           "m_bh": masses.summarise(s), "m_published": masses.spec_summary(p["m_bh"]),
                           "fm_from_P_K": float(fm) if fm is not None else None, "ref": p["ref"]})
    gaia = []
    for g in lit.GAIA_BHS:
        s = masses.draw(g["m_bh"], n, rng, lo=0)
        draws[g["name"]] = {"m1": s}
        gaia.append({"name": g["name"], "m_bh": masses.summarise(s), "porb_d": g["porb_d"],
                     "feh": g["feh"], "ref": g["ref"]})

    v404 = next(r for r in rows if r["name"] == "GS 2023+338")
    n_dyn_2016 = len(vizier)
    names_2016 = set(vizier["Name"].str.upper())
    names_now = {r["name"].upper() for r in rows}
    res = {
        "n_objects_online": len(systems),
        "n_dynamical_online": len(rows),
        "n_dynamical_2016": n_dyn_2016,
        "new_since_2016": sorted(names_now - names_2016),
        "dropped_since_2016": sorted(names_2016 - names_now),
        "systems": rows,
        "persistent": persistent,
        "gaia": gaia,
        "cygx1_check": cyg_check,
        "v404_check": {"fm_ours": v404["fm_ours"], "published": lit.CHECKS["v404_fm"][:2]},
        "fm_above_2_9": [r["name"] for r in rows if r["fm_ours"]["median"] > NS_MAX_KALOGERA_BAYM],
        "fm_above_3_2": [r["name"] for r in rows if r["fm_ours"]["median"] > NS_MAX_CAUSAL],
        "fm_inconsistent": [{"name": r["name"], "pull": r["fm_pull"], "ours": r["fm_ours"]["median"],
                             "published": r["fm_published"]}
                            for r in rows if r["fm_pull"] is not None and abs(r["fm_pull"]) > 2],
        "vizier_2016": vizier.replace({np.nan: None}).to_dict(orient="records"),
    }
    return res, {"draws": draws, "adopted": adopted}


# ---------------------------------------------------------------------------
# 2. The mass gap: hierarchical population fits
# ---------------------------------------------------------------------------

def _catnoir_spec(field: dict):
    t, v, u = field.get("type"), field.get("value"), field.get("uncertainty") or {}
    if t is None:
        return None
    if t == "u":
        return ("U", u["down"], u["up"])
    if t == "an" and u.get("up") != u.get("down"):
        return ("SN", v, u["up"], u["down"])
    return ("N", v, u["up"])


def _norm_name(name: str) -> str:
    return re.sub(r"[^A-Z0-9+-]", "", name.upper().replace("GAIA", "GAIA "))


def catnoir_inputs(names: list[str]) -> dict:
    """CatNoir's (P, K, q, i) for the BlackCAT systems it also lists (by name).
    Its uniform inclinations are isotropic (Rocha et al. 2026, Sect. 3.1)."""
    cat = json.loads(CATNOIR_PATH.read_text())
    alias = {"1H 1659-487": "GX 339-4"}
    by = {_norm_name(c["name"]): c for c in cat}
    out = {}
    for name in names:
        c = by.get(_norm_name(alias.get(name, name)))
        if c is None:
            continue
        specs = {k: _catnoir_spec(c[k]) for k in ("P_orb", "K_cp", "q", "orb_angle")}
        if any(v is None for v in specs.values()):
            continue
        if specs["orb_angle"][0] == "U":
            specs["orb_angle"] = ("I",) + specs["orb_angle"][1:]
        out[name] = {"specs": specs, "m_bh": c["m_bh"], "type": c["Type"]}
    return out


def analyse_mass_gap(weigh_art: dict, rng, quiet=False) -> tuple[dict, dict]:
    n = N_POP_SAMPLES
    samples = {}
    samples["farr2011"] = {k: masses.sample_lmxb_mass(f, q, i, n, rng)
                           for k, (f, q, i) in lit.FARR2011_LOW.items()}
    samples["kreidberg2012"] = {k: masses.sample_lmxb_mass(f, q, i, n, rng)
                                for k, (f, q, i, _, _) in lit.KREIDBERG2012.items()}
    lmxb = {}
    for a in weigh_art["adopted"]:
        d = blackcat_mass_samples(a, n, rng)
        if "m1" in d:
            lmxb[a["name"]] = d["m1"]
    samples["blackcat_lmxb"] = lmxb
    combined = dict(lmxb)
    for name in [p["name"] for p in lit.PERSISTENT_BHS] + [g["name"] for g in lit.GAIA_BHS]:
        src = weigh_art["draws"][name]["m1"]
        combined[name] = rng.choice(src, n, replace=True)
    samples["blackcat_combined"] = combined
    cn = catnoir_inputs(list(lmxb))
    samples["catnoir_lmxb"] = {
        k: masses.sample_lmxb_mass(
            masses.mass_function(masses.draw(v["specs"]["P_orb"], n, rng, lo=0),
                                 masses.draw(v["specs"]["K_cp"], n, rng, lo=0)),
            v["specs"]["q"], v["specs"]["orb_angle"], n, rng)
        for k, v in cn.items()}

    fits, art = {}, {"samples": samples, "chains": {}}
    for label, smp in samples.items():
        fits[label] = {}
        for model in masses.MODELS:
            seed = RANDOM_SEED + zlib.crc32(f"{label}/{model}".encode()) % 1000
            r = masses.fit_population(smp, model, label=label, nlive=POP_NLIVE,
                                      seed=int(seed), quiet=quiet)
            art["chains"][(label, model)] = r.pop("_chain")
            art.setdefault("m1pct", {})[(label, model)] = r.pop("_m1pct")
            fits[label][model] = r
    # model comparison within each sample (Bayes factors relative to the best)
    for label in fits:
        best = max(fits[label].values(), key=lambda r: r["ln_evidence"])["ln_evidence"]
        for r in fits[label].values():
            r["ln_bayes_vs_best"] = r["ln_evidence"] - best

    per_system = {}
    for k in lmxb:
        row = {"blackcat": masses.summarise(lmxb[k])}
        if k in samples["catnoir_lmxb"]:
            row["catnoir_inputs"] = masses.summarise(samples["catnoir_lmxb"][k])
            row["catnoir_published"] = cn[k]["m_bh"]
        per_system[k] = row
    res = {
        "samples": {k: list(v) for k, v in samples.items()},
        "fits": fits,
        "per_system_blackcat_vs_catnoir": per_system,
        "checks": {"farr2011": lit.CHECKS["farr2011_m1pct"], "kreidberg2012": lit.CHECKS["kreidberg2012"],
                   "rocha2026": lit.CHECKS["rocha2026_m1pct_90"]},
    }
    return res, art


# ---------------------------------------------------------------------------
# 3. Inclination as the dominant systematic (Kreidberg et al. 2012)
# ---------------------------------------------------------------------------

def analyse_inclination(weigh_res: dict, weigh_art: dict, rng) -> dict:
    """A0620-00's history, and how much ignored disc light would put each
    black hole into the gap."""
    a0620 = next(a for a in weigh_art["adopted"] if a["name"] == "3A 0620-003")
    d = weigh_art["draws"]["3A 0620-003"]
    fm, q = np.median(d["fm"]), np.median(d["q"])
    history = {  # from Kreidberg et al. 2012, Sect. 1
        "Shahbaz et al. 1994": {"i": 37.0, "m_quoted": 16.0},
        "Cantrell et al. 2010": {"i": 51.0, "m_quoted": 6.6},
        "Froning & Robinson 2001": {"i": 74.0, "m_quoted": 4.0},
    }
    for v in history.values():
        v["m_ours"] = float(masses.m1_from_fm(fm, q, v["i"]))
    gap_edge = 5.0
    distance = []
    for r in weigh_res["systems"]:
        if not r["in_population"]:
            continue
        dd = weigh_art["draws"][r["name"]]
        fm_s, q_s, i_s = np.median(dd["fm"]), np.median(dd["q"]), np.median(dd["incl"])
        s5 = (fm_s * (1 + q_s) ** 2 / gap_edge) ** (1 / 3)
        if s5 >= 1:
            i5, disc = 90.0, None
        else:
            i5 = float(np.degrees(np.arcsin(s5)))
            # i_fit (the adopted value) would be what a model that ignores a disc
            # fraction d infers if the truth were i5: sin^2 i_fit = (1-d) sin^2 i5
            disc = float(1 - np.sin(np.radians(i_s)) ** 2 / s5**2) if i_s < i5 else 0.0
        distance.append({"name": r["name"], "common_name": r["common_name"],
                         "m1_median": r["m1_ours"]["median"], "i_median": float(i_s),
                         "i_for_5msun": i5, "delta_i": i5 - float(i_s), "disc_fraction_needed": disc})
    distance.sort(key=lambda x: x["delta_i"])
    return {"a0620": {"fm": float(fm), "q": float(q), "history": history,
                      "adopted_i": masses.spec_summary(a0620["i_spec"])},
            "veiling_examples": {f"{d:.1f}": float(masses.veiling_mass_bias(d))
                                 for d in (0.1, 0.2, 0.3, 0.4, 0.5)},
            "distance_to_gap": distance}


# ---------------------------------------------------------------------------
# 4. Light versus gravitational waves
# ---------------------------------------------------------------------------

def analyse_graveyard(weigh_art: dict, rng) -> tuple[dict, dict]:
    ns = graveyard.load_freire_ns()
    gw = graveyard.gw_selection_weights(graveyard.load_gw_events())
    bbh = gw[(gw["kind1"] == "BH") & (gw["kind2"] == "BH")].copy()
    for col in ("w_chirp", "w_vmax"):
        bbh[col] = bbh[col] / bbh[col].sum()
    comp_m = np.concatenate([bbh["m1"], bbh["m2"]])
    q = [0.1, 0.5, 0.9]
    gw_stats = {"observed": {"m1": list(np.quantile(bbh["m1"], q)), "all": list(np.quantile(comp_m, q))}}
    for col in ("w_chirp", "w_vmax"):
        w = bbh[col].values
        gw_stats[col] = {
            "m1": list(graveyard.weighted_quantile(bbh["m1"], w, q)),
            "all": list(graveyard.weighted_quantile(comp_m, np.concatenate([w, w]), q)),
            "frac_m1_above_20": float(w[bbh["m1"].values > 20].sum()),
            "n_eff": float(1.0 / np.sum(w**2)),  # effective number of events after weighting
        }
    gw_stats["observed"]["frac_m1_above_20"] = float(np.mean(bbh["m1"] > 20))
    # X-ray binaries: every system counts once (mixture of their mass distributions)
    xrb_names = [k for k, v in weigh_art["draws"].items()
                 if "m1" in v and not k.startswith("Gaia")]
    xrb_mix = np.concatenate([rng.choice(weigh_art["draws"][k]["m1"], 4000) for k in xrb_names])
    xrb_stats = {"names": xrb_names, "quantiles": list(np.quantile(xrb_mix, q)),
                 "frac_above_20": float(np.mean(xrb_mix > 20)),
                 "median_of_medians": float(np.median([np.median(weigh_art["draws"][k]["m1"])
                                                       for k in xrb_names]))}
    special = gw[gw["name"].str.match(r"GW(170817|190425|190814|200105|200115|230529|231123|190521_03|250114)")]
    res = {
        "n_ns_pulsar_timing": int(len(ns)),
        "ns_max": ns.loc[ns["mass"].idxmax()].to_dict(),
        "ns_min": ns.loc[ns["mass"].idxmin()].to_dict(),
        "n_ns_xrb": len(lit.NS_XRB),
        "n_gw_events": int(len(gw)),
        "gw_catalogs": gw["catalog"].value_counts().to_dict(),
        "n_bbh": int(len(bbh)),
        "gw_stats": gw_stats,
        "xrb_stats": xrb_stats,
        "special_events": special[["name", "m1", "m2", "mf", "dl_mpc", "snr"]].to_dict(orient="records"),
    }
    return res, {"ns": ns, "gw": gw, "bbh": bbh, "xrb_mix": xrb_mix}


# ---------------------------------------------------------------------------
# 5. Gaia's dormant black holes
# ---------------------------------------------------------------------------

def analyse_gaia(rng) -> tuple[dict, dict]:
    res, art = {}, {}
    for name, sid in GAIA_BH_SOURCE_IDS.items():
        g = next(x for x in lit.GAIA_BHS if x["name"] == name)
        r = gaia_orbits.analyse(name, sid, g["m_star"], rng)
        art[name] = {"draws": r.pop("_draws"), "row": r.pop("_row")}
        if name == "Gaia BH2":
            r1 = gaia_orbits.analyse(name, sid, ("fixed", 1.0), rng, n=100_000)
            r["m2_for_mstar_1"] = r1["m2"]
        r["check"] = lit.GAIA_DR3_ONLY_CHECKS[name]
        res[name] = r
    return res, art


# ---------------------------------------------------------------------------
# 6. Outbursts
# ---------------------------------------------------------------------------

def analyse_outbursts(weigh_art: dict, rng) -> tuple[dict, dict]:
    crab = outbursts.CrabYardstick()
    j1820, j1820_df = outbursts.j1820_analysis(crab, rng)
    gx, gx_df = outbursts.gx339_outbursts(crab)
    # MAXI J1820+070 calibrates the soft-to-hard luminosity (2-20 keV, VLBI
    # distance); applied to GX 339-4's transitions with its BlackCAT mass,
    # it gives a distance. A demonstration of 'transitions as standard candles'.
    s2h = j1820["transitions"]["soft_to_hard_mjd"]
    f_j = float(j1820_df.loc[np.abs(j1820_df["mjd"] - s2h) <= 1, "eflux_2-20"].mean())
    lam = outbursts.eddington_ratio_samples(f_j, rng)  # L(2-20)/L_Edd at the transition
    m_gx = rng.choice(weigh_art["draws"]["1H 1659-487"]["m1"], len(lam))
    f_gx = np.array([o["soft_to_hard_crab_2_20"] for o in gx if o["kind"] == "full"
                     and np.isfinite(o.get("soft_to_hard_crab_2_20", np.nan))])
    crab_e = sum(outbursts.crab_energy_flux(*b) for b in outbursts.MAXI_BANDS.values())
    f_gx_cgs = f_gx * crab_e  # Crab units -> erg/cm^2/s (Crab-like spectrum assumed)
    f_pick = rng.choice(f_gx_cgs, len(lam))  # transition-to-transition scatter included
    d_kpc = np.sqrt(lam * L_EDD_PER_MSUN * m_gx / (4 * np.pi * f_pick)) / KPC_CM
    lo, med, hi = np.quantile(d_kpc, [0.16, 0.5, 0.84])
    res = {
        "crab": crab.summary(),
        "j1820": j1820,
        "j1820_checks": {"peak_ledd": lit.CHECKS["j1820_peak_ledd"],
                         "soft_to_hard_ledd": lit.CHECKS["j1820_soft_to_hard_ledd"],
                         "soft_to_hard_flux_bol": lit.CHECKS["j1820_soft_to_hard_flux_cgs"],
                         "homan_transition": [58302, 58308]},
        "j1820_soft_to_hard_eflux_2_20": f_j,
        "gx339": gx,
        "gx339_distance_kpc": {"median": float(med), "plus": float(hi - med), "minus": float(med - lo),
                               "n_transitions": int(len(f_gx))},
    }
    return res, {"j1820": j1820_df, "gx339": gx_df, "crab": crab}


# ---------------------------------------------------------------------------

def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items() if not str(k).startswith("_")}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist())
    if isinstance(x, (pd.Timestamp,)):
        return str(x)
    return x


def save_results(results: dict) -> None:
    RESULTS_PATH.write_text(json.dumps(_jsonable(results), indent=1))
