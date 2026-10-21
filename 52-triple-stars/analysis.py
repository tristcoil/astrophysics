"""Collect every number quoted in the README into data/results.json.

    python analysis.py

Reads the caches written by the expensive steps (fits, simulations) and
recomputes the cheap ones.
"""

from __future__ import annotations

import json

import numpy as np

import config as C


def _json(x):
    if isinstance(x, dict):
        return {k: _json(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_json(v) for v in x]
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, np.ndarray):
        return x.tolist()
    return x


def architecture_results():
    import architecture as A
    res = A.run_all()
    res["tokovinin2017"].pop("running", None)
    res["tokovinin2017"].pop("two_orbit_phis", None)
    for k in ("tokovinin2017_model",):
        res[k].pop("lnL", None)
    res["msc2026_alignment"]["model"].pop("lnL", None)
    res["tokovinin2006"].pop("periods", None)
    res["tokovinin2006"].pop("has_tertiary", None)
    return res


def proxima_results():
    import proxima as X
    out = {}
    hAB = None
    ab = C.SYSTEMS["alpha Cen"]["a21"]["AB_orbit"]
    i, O = np.radians(ab["i"]), np.radians(ab["Omega"])
    hAB = np.array([-np.sin(O) * np.sin(i), np.cos(O) * np.sin(i), np.cos(i)])
    for key, f in (("kervella2017_inputs", X.kervella2017), ("gaia_dr3_akeson2021", X.update_gaia_dr3)):
        r = f()
        s = X.summarize(r)
        mc = r["mc"]
        ii, OO = np.radians(mc["i"]), np.radians(mc["Omega"])
        hP = np.stack([-np.sin(OO) * np.sin(ii), np.cos(OO) * np.sin(ii), np.cos(ii)], -1)
        phi = np.degrees(np.arccos(np.clip(hP @ hAB, -1, 1)))
        q = np.percentile(phi, [16, 50, 84])
        s["mutual_inclination_vs_AB"] = dict(median=float(q[1]), lo=float(q[1] - q[0]), hi=float(q[2] - q[1]))
        s["plane_angle_vs_AB"] = float(min(q[1], 180 - q[1]))
        s["light_time_days"] = r["light_time_days"]
        out[key] = s
    out["published"] = C.SYSTEMS["alpha Cen"]["k17"]["out"]
    return out


def kozai_results():
    import kozai as K
    out = {"timescales": K.timescales(),
           "variant_bulge_timescales_yr": {k: K.timescales(v)["t_tide"] for k, v in K.VARIANTS.items()},
           "apsidal_constants": K.K2_APSIDAL}
    for lab in ("newton", "gr", "tides", "tides_Bsoft", "tides_Aonly"):
        path = K.OUT / f"algol_{lab}.json"
        if path.exists():
            r = json.loads(path.read_text())
            out[lab] = dict(touched_yr=r["touched_yr"], e_max=r["e_max"], years=r["t_yr"][-1],
                            J_range=[float(np.min(r["J"])), float(np.max(r["J"]))])
    return out


def stability_results():
    import stability as St
    out = {}
    q_out = C.SYSTEMS["TIC 290061484"]["masses"]["B"] / (C.SYSTEMS["TIC 290061484"]["masses"]["Aa"]
                                                       + C.SYSTEMS["TIC 290061484"]["masses"]["Ab"])
    for kind in ("ratio_e", "ratio_i"):
        path = St.OUT / f"megno_{kind}.npz"
        if not path.exists():
            continue
        d = np.load(path)
        m, b = d["megno"], d["broken"].astype(bool)
        rows = []
        for yi, y in enumerate(d["ys"]):
            reg = (m[yi] < 2.5) & ~b[yi]
            idx = next((i for i in range(len(d["ratios"])) if reg[i:].all()), None)
            crit = C.ma01_ratio(q_out, y, 9.17) if kind == "ratio_e" else C.ma01_ratio(q_out, 0.2011, y)
            rows.append(dict(y=float(y), regular_from_ratio=None if idx is None else float(d["ratios"][idx]),
                             ma01=float(crit)))
        out[kind] = dict(n_points=int(m.size), n_out=int(d["n_out"]), frac_regular=float(np.mean((m < 2.5) & ~b)),
                         frac_chaotic=float(np.mean((m >= 2.5) & ~b)), frac_broken=float(np.mean(b)),
                         boundary=rows)
    path = St.OUT / "disrupt_ratio2.50.json"
    if path.exists():
        r = json.loads(path.read_text())
        res = r["results"]
        ej = [x for x in res if x["ejected"]]
        n_ej = len(ej)
        frac = {k: sum(1 for x in ej if x["ejected"] == k) / n_ej for k in ("Aa", "Ab", "B")}
        out["disruption"] = dict(n=len(res), n_ejected=n_ej, survived=len(res) - n_ej, t_max_out=r["t_max_out"],
                                 fractions=frac, fraction_err={k: float(np.sqrt(v * (1 - v) / n_ej)) for k, v in frac.items()},
                                 theory=St.ejection_theory(),
                                 t_out_orbits_p10_p50_p90=[float(np.percentile([x["t_out_orbits"] for x in ej], q))
                                                           for q in (10, 50, 90)],
                                 binary_a_median=float(np.median([x["a_bin"] for x in ej])),
                                 binary_e_median=float(np.median([x["e_bin"] for x in ej])),
                                 v_inf_median={k: float(np.median([x["v_inf"] for x in ej if x["ejected"] == k]))
                                               for k in ("Aa", "Ab", "B")})
    sysT = C.SYSTEMS["TIC 290061484"]
    out["tic290061484_ma01_ratio"] = float(C.ma01_ratio(q_out, sysT["A_B"]["e"], sysT["pub"]["i_mut"]))
    out["tic290061484_actual_ratio"] = float(sysT["A_B"]["a_rsun"] / sysT["Aa_Ab"]["a_rsun"])
    return out


def _published_switching(name):
    pub = dict(C.SYSTEMS[name]["pub"])
    pub.pop("inclinations", None)
    return pub


def zasche_forecast():
    """V907 Sco: Zasche et al. (2023)'s precession model (their J, I, P_node,
    m3; quadrupole split of the angular momentum), anchored to their own
    2019-2022 inclinations, against the tilts TESS measured afterwards."""
    from scipy.optimize import minimize_scalar
    import switching as S
    import switch_fit as F
    s = C.SYSTEMS["V907 Sco"]
    pub = s["pub"]
    eps_in, _, nu, gam, _, _ = S.secular(s["m1"], s["m2"], pub["m3"], s["P_in"], s["P_out"], s["e_in"], s["e_out"],
                                         np.radians(pub["J"]))
    I = np.radians(pub["i_tot"])
    rate = 2 * np.pi / (pub["P_node"] * 365.25)
    tab = [r for r in pub["inclinations"] if r[0] > 58000]
    t = np.array([r[0] for r in tab]) + 2400000.0
    i_obs = np.array([r[1] for r in tab])
    i_err = np.array([r[2] for r in tab])

    def incl(h1, sgn, tt):
        h = h1 + sgn * rate * (np.asarray(tt) - t[1])
        return np.degrees(np.arccos(np.cos(I) * np.cos(eps_in) - np.sin(I) * np.sin(eps_in) * np.cos(h)))
    best = None
    for sgn in (-1, 1):
        for h0 in np.linspace(0, 2 * np.pi, 73):
            r = minimize_scalar(lambda h: np.sum(((incl(h, sgn, t) - i_obs) / i_err) ** 2),
                                bracket=(h0 - 0.05, h0 + 0.05))
            if best is None or r.fun < best[0]:
                best = (r.fun, r.x, sgn)
    chi2, h1, sgn = best
    P = json.loads((F.CACHE / "profiles_v907_sco.json").read_text())
    rs = (s["R1"] + s["R2"]) / s["a_in_rsun"]
    rows = []
    for pr in P:
        if not pr["detected"]:
            continue
        x = np.array(pr["x"])
        c = (np.array(pr["chi2"]) - pr["chi2_min"]) / pr["scale"]
        k = int(np.argmin(c))
        sl = slice(max(0, k - 2), k + 3)
        a2, a1, _ = np.polyfit(x[sl], c[sl], 2)
        xb = -a1 / (2 * a2)
        ib = float(np.degrees(np.arccos(xb * rs)))
        rows.append(dict(sector=pr["sector"], year=float(C.jd_to_year(pr["jd_mid"])), tilt_measured=ib,
                         tilt_forecast=float(incl(h1, sgn, pr["jd_mid"])),
                         after_paper=bool(pr["jd_mid"] > 2459830.0)))
    yrs = np.arange(2019, 2070, 0.1)
    ii = incl(h1, sgn, C.year_to_jd(yrs))
    return dict(anchor_chi2=float(chi2), n_anchor=int(len(t)), eps_in_deg=float(np.degrees(eps_in)),
                quadrupole_P_node_yr=float(2 * np.pi / abs(nu) / 365.25), sectors=rows,
                eclipses_end_model=float(yrs[ii > pub["i_limit"]].max()),
                note="tilts from our radii-fixed impact-parameter profiles; S39 sits 0.13 deg above Zasche's own value")


def tilt_history_stats(name, tag="", n=600, seed=C.RANDOM_SEED):
    """From posterior draws: when the orbit was edge-on between 1850 and 2000
    (unfolded tilt crossing 90 deg; fraction of draws that do), the mean tilt
    rate 1900-1950 (deg/yr, folded), the sum of fractional radii and the
    eclipse limit (deg)."""
    import switch_fit as F
    import switching as S
    post = np.load(F.posterior_path(name, tag), allow_pickle=True)["posterior"]
    rng = np.random.default_rng(seed)
    m = S.Model(name)
    yrs = np.arange(1850.0, 2000.0, 0.25)
    jd = C.year_to_jd(yrs)
    t90, rate, rsum = [], [], []
    for th in post[rng.choice(len(post), min(n, len(post)), replace=False)]:
        g = m.geometry(m.unpack(th))
        i = np.degrees(m.inclination(g, jd))
        k = np.where(np.diff(np.sign(i - 90)) != 0)[0]
        t90.append(yrs[k[0]] if len(k) else np.nan)
        f = 90 - np.abs(90 - np.degrees(m.inclination(g, C.year_to_jd(np.array([1900.0, 1950.0])))))
        rate.append((f[1] - f[0]) / 50.0)
        rsum.append(g["r1"] + g["r2"])
    q = lambda a: [float(v) for v in np.nanpercentile(np.asarray(a, float), [16, 50, 84])]
    return dict(t_edge_on=q(t90), frac_edge_on_found=float(np.mean(np.isfinite(t90))),
                tilt_rate_1900_1950=q(rate), r_sum=q(rsum), i_edge=q(np.degrees(np.arccos(rsum))))


def retrograde_comparison(name, tag, s_main):
    """The best solution with the other sense of orbital motion and the log
    odds of the final solution against it (switch_fit.alternative_solution)."""
    import switch_fit as F
    alt = F.alternative_solution(name)
    if alt is None:
        return None
    sr = F.summarize(name, tag=alt["tag"])
    keep_p = ("J", "I", "m3")
    keep_d = ("P_node_yr", "last_end", "next_start", "current_end", "i_2026", "season_len")
    return dict(tag=alt["tag"], logZ=sr["logZ"], logZ_err=sr["logZ_err"], lnB_final_over_alt=alt["lnB"],
                lnB_err=alt["lnB_err"], p_alt=float(1.0 / (1.0 + np.exp(alt["lnB"]))),
                params={k: sr["params"][k] for k in keep_p},
                derived={k: sr["derived"].get(k) for k in keep_d})


def ss_lac_quiet_comparison():
    """SS Lac fitted only outside the Kozai-Lidov window, prograde (J < 39.2)
    and retrograde (J > 140.8): both ranges carry the same isotropic prior
    mass, so the evidences compare directly."""
    import switch_fit as F
    out = {}
    for tag in ("quietpro", "quietretro"):
        if F.posterior_path("SS Lac", tag).exists():
            s = F.summarize("SS Lac", tag=tag)
            out[tag] = dict(logZ=s["logZ"], logZ_err=s["logZ_err"],
                            params={k: s["params"][k] for k in ("J", "I", "m3", "r_sum")},
                            derived={k: s["derived"].get(k) for k in ("P_node_yr", "eps_in", "last_end",
                                                                     "next_start", "i_1900", "i_2026",
                                                                     "season_len")})
    if len(out) == 2:
        lnB = out["quietpro"]["logZ"] - out["quietretro"]["logZ"]
        out["lnB_pro_over_retro"] = float(lnB)
        out["lnB_err"] = float(np.hypot(out["quietpro"]["logZ_err"], out["quietretro"]["logZ_err"]))
    return out


def switching_results():
    import switch_fit as F
    import switch_nbody as SN
    out = {}
    for name, tag in F.FINAL_TAG.items():
        if not F.posterior_path(name, tag).exists():
            continue
        s = F.summarize(name, tag=tag)
        r = dict(posterior=tag or "plates+TESS", logZ=s["logZ"], n_posterior=s["n_posterior"],
                 frac_prograde=s["frac_prograde"], params=s["params"], derived=s["derived"],
                 published=_published_switching(name))
        chk = F.CACHE / f"nbody_check_{F._slug(name)}{'_' + tag if tag else ''}.json"
        if chk.exists():
            res = [x for x in json.loads(chk.read_text()) if x.get("ok")]
            quiet = [x for x in res if SN.is_quiet(name, x)]
            r["nbody_check"] = dict(n=len(res), n_quiet=len(quiet),
                                    di_rms_median=float(np.median([x["di_rms"] for x in res])),
                                    di_rms_p95=float(np.percentile([x["di_rms"] for x in res], 95)),
                                    di_max_p95=float(np.percentile([x["di_max"] for x in res], 95)),
                                    e_max_max=float(np.max([x["e_max"] for x in res])))
        r["ground_check"] = F.ground_check(name, tag)
        r["tilt_history"] = tilt_history_stats(name, tag)
        r["retrograde"] = retrograde_comparison(name, tag, s)
        out[name] = r
    out["SS Lac"] = out.get("SS Lac", {})
    out["SS Lac"]["kozai_quiet_runs"] = ss_lac_quiet_comparison()
    unres = F.CACHE / "nbody_check_ss_lac_unrestricted.json"
    if F.posterior_path("SS Lac", "").exists():
        s0 = F.summarize("SS Lac", tag="")
        u = dict(logZ=s0["logZ"], params={k: s0["params"][k] for k in ("J", "I", "m3")},
                 derived={k: s0["derived"][k] for k in ("P_node_yr", "eps_in")})
        if unres.exists():
            res = [x for x in json.loads(unres.read_text()) if x.get("ok")]
            u["nbody_check"] = dict(n=len(res), n_quiet=sum(SN.is_quiet("SS Lac", x) for x in res),
                                    e_max_p5_p50_p95=[float(np.percentile([x["e_max"] for x in res], q))
                                                      for q in (5, 50, 95)],
                                    J_range_median=float(np.median([x["J_range"] for x in res])),
                                    di_rms_median=float(np.median([x["di_rms"] for x in res])))
        out["SS Lac"]["unrestricted"] = u
    if F.posterior_path("HS Hya", "").exists():
        s = F.summarize("HS Hya", tag="")
        out["HS Hya"]["without_ground_data"] = dict(params={k: s["params"][k] for k in ("J", "I", "m3")},
                                                     derived={k: s["derived"][k] for k in
                                                              ("P_node_yr", "last_end", "next_start", "i_2026")},
                                                     ground_check=F.ground_check("HS Hya", ""))
    out["V907 Sco"]["zasche_forecast"] = zasche_forecast()
    return out


def tic290_results():
    import tic290_nbody as N
    SYS = C.SYSTEMS["TIC 290061484"]
    obs = N.observed_times()
    eph = json.loads((C.DERIVED / "tic290061484" / "ephemeris.json").read_text())
    out = dict(n_eclipses=int(len(obs)), median_err_min=float(np.median(obs["err"]) * 1440),
               ephemeris=dict(P=eph["P"], T0=eph["T0"]), published=SYS["pub"])
    scan = N.etv_period_scan(obs)
    out["etv_period_scan"] = dict(P_best=scan["P_best"], P_lo=scan["P_lo"], P_hi=scan["P_hi"],
                                  chi2_scale=scan["chi2_scale"])
    fits = {}
    for v in ("triple", "all"):
        f = N.run_fit(variant=v)
        par = dict(f["base"], **f["best"])
        r = (np.array(f["t_obs"]) - np.array(f["t_model"])) * 1440
        i1, i2 = np.radians(SYS["Aa_Ab"]["i"]), np.radians(par["A_B.i"])
        dO = np.radians(par["A_B.Omega"] - SYS["Aa_Ab"]["Omega"])
        J = np.degrees(np.arccos(np.cos(i1) * np.cos(i2) + np.sin(i1) * np.sin(i2) * np.cos(dO)))
        fits[v] = dict(chi2=f["chi2"], dof=f["dof"], rms_min=float(np.sqrt(np.nanmean(r ** 2))),
                       chi2_start=f["chi2_start"], best=f["best"], err=f["err"], i_mut_deg=float(J),
                       mean_periods=N.outer_mean_periods(par))
    out["fits"] = fits
    base = N.run_fit(variant="triple")["base"]
    out["calibration"] = dict(P_osc_published=SYS["Aa_Ab"]["P"], P_osc_calibrated=base["Aa_Ab.P"],
                              tau_shift_d=base["Aa_Ab.tau"] - SYS["Aa_Ab"]["tau"],
                              chi2=N.run_fit(variant="triple")["chi2_start"])
    out["kostov_published_mean_periods"] = N.outer_mean_periods(None)
    return out


def write_all(sections=None):
    sections = sections or ["architecture", "proxima", "kozai", "stability", "switching", "tic290"]
    res = json.loads(C.RESULTS.read_text()) if C.RESULTS.exists() else {}
    funcs = {"architecture": architecture_results, "proxima": proxima_results, "kozai": kozai_results,
             "stability": stability_results, "switching": switching_results, "tic290": tic290_results}
    for s in sections:
        res[s] = _json(funcs[s]())
    C.RESULTS.write_text(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    import sys
    r = write_all(sys.argv[1:] or None)
    print(json.dumps({k: list(v.keys()) if isinstance(v, dict) else v for k, v in r.items()}, indent=1))
