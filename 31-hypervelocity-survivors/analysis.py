"""Every analysis step, from cached data to data/results.json.

Each step returns a JSON-able results dict and, where figures need them,
arrays ("artefacts"). Expensive steps cache their arrays in data/cache/.
"""

from __future__ import annotations

import json
import warnings

import numpy as np
from astropy.table import Table

import binary
import kinematics as K
import orbits
import population as P
import remnants as REM
import search as S
from config import (AUX_DIR, CACHE_DIR, GAIA_DIR, N_SAMPLES, N_SHUFFLES, N_TRACK_SAMPLES,
                    N_POP_STARS, REMNANT_LOOKBACK_YR, RESULTS_PATH, SN_IA_RATE_ERR,
                    SN_IA_RATE_PER_YR, V_CIRC_FLAT)
from literature import CENSUS, D62_REMNANT, REFERENCE, US708_DISTANCE_PRIOR_KPC

warnings.filterwarnings("ignore", category=RuntimeWarning)


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist())
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def save_results(results: dict) -> None:
    RESULTS_PATH.write_text(json.dumps(_jsonable(results), indent=1))


def known_rows() -> dict:
    t = Table.read(GAIA_DIR / "known.fits")
    return {int(r["source_id"]): r for r in t}


def dust() -> dict:
    return json.loads((AUX_DIR / "irsa_dust.json").read_text())


# --------------------------------------------------------------------------
# 1. The search
# --------------------------------------------------------------------------

def analyse_search() -> tuple[dict, dict]:
    search = S.load_search()
    pos = search[np.asarray(search["plx_snr"]) > 0]
    union = S.union_selection(search)
    eb = S.elbadry_check()
    rows = eb["rows"]
    observed = [r for r in rows if r["cls"] in ("runaway", "false positive")]
    res = dict(
        n_search_rows=len(search), n_positive=len(pos), n_negative=int(np.sum(np.asarray(search["plx_snr"]) < 0)),
        mirror=S.mirror_test(search),
        elbadry=dict(n_returned=eb["n_returned"], identical_to_table1=eb["identical_to_table1"],
                     missing=eb["missing"], extra=eb["extra"], n_observed=len(observed),
                     n_runaways=sum(r["cls"] == "runaway" for r in rows), rows=rows),
        union=dict(n=len(union), n_precise=int(np.sum(np.asarray(union["subset"]) == "plx/err >= 5")),
                   ids=[int(i) for i in union["source_id"]]),
        census=S.census_flags(search, union),
        shen=S.shen_dr2_vs_dr3(),
    )
    res["shift_test"], robust = S.shift_test(search)
    # how well does P(v_tan > 1000 | parallax) separate runaways from impostors?
    run = [r["p_fast1000"] for r in rows if r["cls"] == "runaway"]
    fp = [r["p_fast1000"] for r in rows if r["cls"] == "false positive"]
    res["elbadry"]["p_fast1000_runaways"] = run
    res["elbadry"]["p_fast1000_false"] = fp
    # the 5 highest v_tan_lower: all D6 stars (El-Badry et al.'s remark)
    top = sorted([r for r in rows if r["plx"] + r["plx_err"] > 0], key=lambda r: -r["v_tan_lower"])[:5]
    res["elbadry"]["top5_lower_all_runaways"] = all(r["cls"] == "runaway" for r in top)
    return res, dict(search=search, union=union, robust=robust)


# --------------------------------------------------------------------------
# 2. Kinematics of the known runaways
# --------------------------------------------------------------------------

def _prior_for(r):
    return ("gauss", US708_DISTANCE_PRIOR_KPC) if r.name == "US 708" else ("flat_MG", None)


def analyse_kinematics(rng) -> tuple[dict, dict]:
    rows, ebv = known_rows(), dust()
    res, samples = {}, {}
    for r in CENSUS:
        row = rows[r.source_id]
        e = ebv[str(r.source_id)]["refPixelValueSandF"]
        prior, gauss = _prior_for(r)
        fid = K.posterior(row, r.rv, r.rv_err, e, rng, N_SAMPLES, prior=prior, gauss=gauss,
                          vcirc=orbits.vcirc)
        ebconv = K.posterior(row, r.rv, r.rv_err, e, rng, N_SAMPLES, prior=prior, gauss=gauss,
                             vcirc=V_CIRC_FLAT, zero_point=False)
        vesc = orbits.escape_speed(fid["x"])
        out = dict(name=r.name, cls=r.cls, source_id=r.source_id, rv=r.rv, rv_err=r.rv_err,
                   G=float(row["phot_g_mean_mag"]), bp_rp=float(row["bp_rp"]),
                   plx=float(row["parallax"]), plx_err=float(row["parallax_error"]),
                   pm=float(np.hypot(row["pmra"], row["pmdec"])), ruwe=float(row["ruwe"]),
                   ebv=e, zero_point=fid["meta"]["zero_point"],
                   l=float(row["l"]), b=float(row["b"]), acceptance=fid["meta"]["acceptance"])
        for key in ("d", "v_tan", "v_tot", "v_ej", "z", "t_disk", "M_G0"):
            out[key] = K.summarize(fid[key])
            out[f"eb_conv_{key}"] = K.summarize(ebconv[key])
        out["bp_rp0"] = float(fid["bp_rp0"])
        out["v_esc"] = K.summarize(vesc)
        out["p_unbound"] = float(np.mean(fid["v_tot"] > vesc))
        out["eb23"] = r.eb23
        # prior sensitivity (fiducial conventions otherwise)
        if prior != "gauss":
            out["priors"] = {}
            for alt in ("flat_d", "edsd"):
                s_alt = K.posterior(row, r.rv, r.rv_err, e, rng, 2000, prior=alt, vcirc=orbits.vcirc)
                out["priors"][alt] = {k: K.summarize(s_alt[k]) for k in ("d", "v_tot", "v_ej")}
        res[r.name] = out
        samples[r.name] = fid
    return res, samples


# --------------------------------------------------------------------------
# 3. Orbits: launch speeds and birthplaces with the potential
# --------------------------------------------------------------------------

def analyse_orbits(samples, rng, n_orbit=1000) -> tuple[dict, dict]:
    res, art = {}, {}
    t_back = np.linspace(0, -15.0, 1501)
    t_plot = np.linspace(-3.0, 3.0, 301)
    for name, s in samples.items():
        k = rng.choice(s["d"].size, size=min(n_orbit, s["d"].size), replace=False)
        x, v = s["x"][:, k], s["v"][:, k]
        pos, vel = orbits.integrate(x, v, t_back)
        tc, xc, v_ej = orbits.launch_from_orbit(pos, vel, t_back)
        # speed lost climbing out: present speed vs speed at the plane crossing
        crossed = np.isfinite(tc)
        res[name] = dict(
            v_ej_potential=K.summarize(v_ej),
            v_ej_straight=K.summarize(s["v_ej"][k]),
            frac_crossed_15myr=float(crossed.mean()),
            t_cross=K.summarize(tc[crossed]) if crossed.sum() > 10 else None,
            R_birth=K.summarize(np.hypot(xc[0], xc[1])[crossed]) if crossed.sum() > 10 else None,
        )
        # trajectories for the figures: 120 samples, past and future
        kk = k[:120]
        fut, _ = orbits.integrate(s["x"][:, kk], s["v"][:, kk], np.linspace(0, 3.0, 151))
        past, _ = orbits.integrate(s["x"][:, kk], s["v"][:, kk], np.linspace(0, -3.0, 151))
        art[name] = dict(past=past.astype(np.float32), future=fut.astype(np.float32),
                         xc=xc[:, crossed][:, :300].astype(np.float32))
    return res, art


# --------------------------------------------------------------------------
# 4. Remnants
# --------------------------------------------------------------------------

def analyse_remnants(samples, rng) -> tuple[dict, dict]:
    green = REM.load_green()
    t_myr = np.linspace(0, -REMNANT_LOOKBACK_YR / 1e6, 301)
    t_yr = -t_myr * 1e6
    res, art = dict(n_remnants=len(green)), {}
    for name, s in samples.items():
        k = rng.choice(s["d"].size, size=min(N_TRACK_SAMPLES, s["d"].size), replace=False)
        pos, _ = orbits.integrate(s["x"][:, k], s["v"][:, k], t_myr)
        l, b, dist = REM.apparent_site_track(pos, t_myr, orbits.vcirc)
        matches = REM.match_samples(l, b, dist, t_yr, green, REMNANT_LOOKBACK_YR)
        best = [m for m in matches if m["frac_inside"] > 0.05][:5]
        # median track (circular mean in longitude)
        lm = np.degrees(np.angle(np.mean(np.exp(1j * np.radians(l)), axis=0))) % 360
        bm = np.median(b, axis=0)
        entry = dict(matches=best, closest=matches[:3])
        if best:
            top = best[0]
            j = int(np.flatnonzero(np.asarray(green["name"]) == top["name"])[0])
            sep, tt = REM.closest_approach(lm, bm, t_yr, [green["l"][j]], [green["b"][j]])
            ratio = float(sep[0] / green["radius_deg"][j])
            entry["best"] = dict(name=top["name"], sep_deg=float(sep[0]), t_yr=float(tt[0]),
                                 radius_deg=float(green["radius_deg"][j]), sep_over_radius=ratio,
                                 frac_inside=top["frac_inside"],
                                 d_star_kpc=top["d_star_median"])
        res[name] = entry
        art[name] = dict(l=l[:150].astype(np.float32), b=b[:150].astype(np.float32),
                         lm=lm, bm=bm, t_yr=t_yr)
    # Look-elsewhere: we searched every star's track, so the relevant chance
    # is that *any* of them passes some remnant as deep inside as the best
    # real match does.
    best_ratio = min((e["best"]["sep_over_radius"] for k, e in res.items()
                      if isinstance(e, dict) and "best" in e), default=None)
    if best_ratio is not None:
        per, full = {}, {}
        for name, a in art.items():
            full[name] = REM.chance_probability(a["lm"], a["bm"], green, rng, N_SHUFFLES, best_ratio)
            per[name] = full[name]["p_chance"]
            if "best" in res[name]:
                res[name]["chance"] = full[name]
        d6 = [n for n in per if next(r.cls for r in CENSUS if r.name == n) == "D6"]
        res["look_elsewhere"] = dict(
            ratio=best_ratio, p_per_star=per, n_shuffles=N_SHUFFLES,
            n_in_band={n: f["n_remnants_in_band"] for n, f in full.items()},
            p_any_d6=float(1 - np.prod([1 - per[n] for n in d6])),
            p_any_all=float(1 - np.prod([1 - p for p in per.values()])))
    # Shen et al. (2018)'s statement for D6-2, checked
    d62 = res.get("D6-2", {})
    j = int(np.flatnonzero(np.asarray(green["name"]) == D62_REMNANT["name"])[0])
    a = art["D6-2"]
    k90 = int(np.argmin(np.abs(a["t_yr"] - D62_REMNANT["lookback_yr"])))
    sep90 = REM._sep_deg(a["lm"][k90], a["bm"][k90], green["l"][j], green["b"][j])
    d62["shen_check"] = dict(sep_at_90kyr_deg=float(sep90), published_deg=D62_REMNANT["offset_deg"],
                             remnant_radius_deg=float(green["radius_deg"][j]),
                             remnant_distance_kpc=D62_REMNANT["distance_kpc"])
    # Is a remnant this big plausible at this age? At D6-2's distance its
    # radius is R pc; a Sedov-Taylor blast wave, R = 1.15 (E t^2 / rho)^(1/5),
    # reaches R at age t if the surrounding gas has density n. (By ~1e5 yr a
    # remnant is radiative and grows more slowly, so this n is a lower bound.)
    if "best" in d62:
        d_kpc = d62["best"]["d_star_kpc"]
        R_pc = np.radians(d62["best"]["radius_deg"]) * d_kpc * 1e3
        t_s = d62["best"]["t_yr"] * 3.156e7
        R_cm = R_pc * 3.0857e18
        rho = 1.15 ** 5 * 1e51 * t_s ** 2 / R_cm ** 5
        d62["sedov"] = dict(radius_pc=R_pc, age_yr=d62["best"]["t_yr"], energy_erg=1e51,
                            n_cm3=rho / (1.4 * 1.6726e-24))
    return res, dict(tracks=art, green=green)


# --------------------------------------------------------------------------
# 5. Donor masses from launch speeds
# --------------------------------------------------------------------------

def analyse_binary(kin) -> dict:
    res = dict(
        check=dict(v_06_on_10_bedard=float(binary.ejection_speed(0.6, 1.0, with_kick=False)),
                   v_06_on_10_eggleton=float(binary.ejection_speed(0.6, 1.0, "eggleton", with_kick=False)),
                   v_10_on_10_bedard=float(binary.ejection_speed(1.0, 1.0, with_kick=False)),
                   v_10_on_10_eggleton=float(binary.ejection_speed(1.0, 1.0, "eggleton", with_kick=False)),
                   shen18=[REFERENCE["shen18_vrun_0.6+1.0"], REFERENCE["shen18_vrun_1.0+1.0"]],
                   vmax_two_105=binary.max_speed(1.05, with_kick=False),
                   vmax_two_105_kick=binary.max_speed(1.05),
                   vmax_two_120=binary.max_speed(1.20, with_kick=False),
                   vmax_two_105_eggleton=binary.max_speed(1.05, "eggleton", with_kick=False),
                   eb23_vmax=REFERENCE["eb23_vmax_two_1.05"]),
        stars={})
    # El-Badry et al.'s Fig. 14 limits, recomputed at two donor temperatures
    res["eb23_fig14"] = {}
    for label, v in (("J0546+0836", 1298), ("J1332-3541", 1448), ("J1235-3752", 2126),
                     ("J0927-6335", 2372)):
        res["eb23_fig14"][label] = dict(
            v=v, teff10k=binary.min_donor_mass(v, 1.05, with_kick=False),
            teff20k=binary.min_donor_mass(v, 1.05, teff=2e4, with_kick=False),
            teff10k_kick=binary.min_donor_mass(v, 1.05))
    for name, k in kin.items():
        if k["cls"] != "D6":
            continue
        v_lo = k["v_ej"]["p2"]           # 2-sigma lower bound
        v_med = k["v_ej"]["median"]
        res["stars"][name] = dict(
            v_ej_2sigma_lower=v_lo, v_ej_median=v_med,
            mdon_min_105=binary.min_donor_mass(v_lo, 1.05),
            mdon_min_120=binary.min_donor_mass(v_lo, 1.20),
            mdon_min_105_cold=binary.min_donor_mass(v_lo, 1.05, "eggleton"),
            mdon_med_105=binary.min_donor_mass(v_med, 1.05),
            mdon_med_120=binary.min_donor_mass(v_med, 1.20))
    return res


# --------------------------------------------------------------------------
# 6. Population: survivors seen vs the Type Ia rate
# --------------------------------------------------------------------------

POP_TAUS_MYR = np.geomspace(0.01, 30.0, 25)
#: Absolute magnitude splitting the D6 stars into El-Badry et al.'s two
#: groups: bright and fast (M_G ~ 5-6.5) vs faint and nearby (M_G > 7).
FAINT_MG = 6.8


def analyse_population(kin, union_ids, rng, remnant_ages=None) -> tuple[dict, dict]:
    """1/V_max-style estimate: each survivor found by the union selection
    stands for a population of look-alikes (its own launch speed, absolute
    magnitude and colour) whose detection fraction the forward model gives.
    The implied birth rate of survivors is the sum over found stars of
    1 / (tau x detection fraction), for a visible lifetime tau.

    A star seen at flight time t has been visible for at least t, so tau >= t
    for its group, and the implied rate falls with tau: the fraction at the
    group's oldest flight time is an upper bound (within the model)."""
    from literature import REFERENCE
    perr = P.ParallaxErrors()
    found = {name: k for name, k in kin.items() if k["source_id"] in set(union_ids)}
    ages = {n: k["t_disk"]["median"] for n, k in found.items() if k["t_disk"]["median"] > 0}
    for n, a in (remnant_ages or {}).items():
        if n in found:
            ages[n] = a
    per_star = {}
    for name, k in found.items():
        v = k["v_ej"]["median"]
        lo, hi = max(v - k["v_ej"]["minus"], 300.0), v + k["v_ej"]["plus"]
        fr = []
        for tau in POP_TAUS_MYR:
            sim = P.simulate(v_range=(lo, hi), M_G0=k["M_G0"]["median"], color0=k["bp_rp0"],
                             tau_myr=tau, n=N_POP_STARS // 10, rng=rng, plx_errors=perr,
                             vcirc_fn=orbits.vcirc)
            fr.append(sim["frac"])
        per_star[name] = dict(cls=k["cls"], v_range=(lo, hi), M_G0=k["M_G0"]["median"],
                              color0=k["bp_rp0"], age_myr=ages.get(name), frac_found=np.array(fr))
        with np.errstate(divide="ignore"):
            per_star[name]["rate_per_yr"] = 1.0 / (POP_TAUS_MYR * 1e6 * np.maximum(np.array(fr), 1e-9))
    iax = REFERENCE["foley13_iax_per_ia"][0]
    groups = {
        "D6 (all)": ([n for n, p in per_star.items() if p["cls"] == "D6"], SN_IA_RATE_PER_YR),
        "D6 bright, fast": ([n for n, p in per_star.items() if p["cls"] == "D6" and p["M_G0"] < FAINT_MG],
                            SN_IA_RATE_PER_YR),
        "D6 faint": ([n for n, p in per_star.items() if p["cls"] == "D6" and p["M_G0"] >= FAINT_MG],
                     SN_IA_RATE_PER_YR),
        "LP40-365 (per Type Iax)": ([n for n, p in per_star.items() if p["cls"] == "LP40-365"],
                                    SN_IA_RATE_PER_YR * iax),
    }
    res = dict(taus_myr=POP_TAUS_MYR, rate=SN_IA_RATE_PER_YR, rate_err=SN_IA_RATE_ERR, iax_per_ia=iax,
               faint_mg=FAINT_MG, found=list(found), ages_myr=ages, groups={},
               per_star={n: dict(cls=p["cls"], M_G0=p["M_G0"], v_range=p["v_range"], age_myr=p["age_myr"],
                                 frac_found=p["frac_found"],
                                 fraction_of_parent_at_age=(float(np.exp(np.interp(
                                     np.log(p["age_myr"]), np.log(POP_TAUS_MYR),
                                     np.log(p["rate_per_yr"])))) /
                                     (SN_IA_RATE_PER_YR * (iax if p["cls"] == "LP40-365" else 1.0)))
                                 if p["age_myr"] else None)
                          for n, p in per_star.items()})
    for g, (names, parent_rate) in groups.items():
        if not names:
            continue
        total = np.sum([per_star[n]["rate_per_yr"] for n in names], axis=0)
        n = len(names)
        lo, hi = P.poisson_interval(n)
        frac = total / parent_rate
        known = [per_star[nm]["age_myr"] for nm in names if per_star[nm]["age_myr"]]
        oldest = max(known) if known else None
        entry = dict(names=names, n_found=n, parent_rate_per_yr=parent_rate, fraction=frac,
                     fraction_lo=frac * lo / n, fraction_hi=frac * hi / n, oldest_age_myr=oldest)
        if oldest is not None:
            entry["fraction_at_oldest"] = float(np.exp(np.interp(np.log(oldest), np.log(POP_TAUS_MYR),
                                                                 np.log(frac))))
            entry["fraction_at_oldest_range"] = (entry["fraction_at_oldest"] * lo / n,
                                                 entry["fraction_at_oldest"] * hi / n)
        res["groups"][g] = entry
    return res, dict(per_star=per_star)
