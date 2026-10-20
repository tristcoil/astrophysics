"""Every number in the README, computed from the cached data.

Each ``*_results`` function returns a JSON-serialisable dict and caches the
heavier arrays the figures need in data/derived/. ``run_all`` writes
data/results.json.
"""

from __future__ import annotations

import json
import os
import warnings

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

import config as C
import literature as L


def _f(x, nd=6):
    """Plain numbers for JSON, kept to ``nd + 3`` significant figures (so tiny
    quantities such as a period derivative of 1e-9 survive)."""
    if isinstance(x, (np.floating, float)):
        x = float(x)
        return x if not np.isfinite(x) else float(f"{x:.{nd + 3}g}")
    if isinstance(x, (np.integer, int)):
        return int(x)
    if isinstance(x, dict):
        return {str(k): _f(v, nd) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [_f(v, nd) for v in x]
    return x


# ---------------------------------------------------------------------------
# Roche geometry
# ---------------------------------------------------------------------------


def roche_results() -> dict:
    import roche as R
    out = {"lagrange": {}, "lobes": []}
    for q in (1.0, 0.1):
        lp = R.lagrange_points(q)
        out["lagrange"][str(q)] = {k: v for k, v in lp.items() if k != "L4"}
    for q in (1.0, 0.5, 0.228, 0.1, 0.05):
        r1 = R.lobe_volume_radius(q, "primary", 64, 64)
        r2 = R.lobe_volume_radius(q, "secondary", 64, 64)
        out["lobes"].append({"q": q, "r1": r1, "r1_eggleton": R.eggleton_radius(1 / q),
                             "r2": r2, "r2_eggleton": R.eggleton_radius(q)})
    # validation: KIC 9832227 with Molnar et al.'s Binary Maker parameters vs Kepler's mean light curve
    import kic9832227 as K
    _, tpl = K.measure_all()  # builds (and caches) the Kepler template if needed
    m = L.MOLNAR2017
    cb = R.ContactBinary(m["q"], m["fillout"], t1=m["t1"], t2=m["t2"], lam_nm=640.0, ld=(0.40, 0.26))
    ph = np.linspace(0, 1, 201)
    model = cb.light_curve(m["incl"], ph, npix=220)
    obs = tpl(ph) / tpl(ph).max()
    np.savez(C.DERIVED_DIR / "roche_kic_model.npz", phase=ph, model=model, kepler=obs)
    out["kic_model"] = {"amp_model_mag": float(-2.5 * np.log10(model.min())),
                        "amp_kepler_mag": float(-2.5 * np.log10(obs.min())),
                        "rms_diff": float(np.std(model - obs))}
    return _f(out)


# ---------------------------------------------------------------------------
# V1309 Sco
# ---------------------------------------------------------------------------

#: Subsamples per season as in Tylenda et al.'s Fig. 2 (17 periods in all).
V1309_NSUB = {2002: 1, 2003: 1, 2004: 3, 2005: 2, 2006: 5, 2007: 5}


def v1309_results(n_boot: int = 200) -> dict:
    import v1309sco as V
    d = V.load_ogle()["I"]
    t, m, e = d["t"], d["mag"], d["err"]
    seas = [s for s in V.seasons(t) if 2002 <= V.jd_to_year(t[s].mean()) <= 2007.9]
    groups, labels = [], []
    for s in seas:
        yr = int(V.jd_to_year(t[s].mean()))
        for k, g in enumerate(V.subsamples(s, int(np.ceil(len(s) / V1309_NSUB[yr])))):
            groups.append(g)
            labels.append(f"{yr}{'abcde'[k] if V1309_NSUB[yr] > 1 else ''}")
    per = []
    for g, lab in zip(groups, labels):
        best, sig, _ = V.bootstrap_period(t[g], m[g], e[g], n_boot=n_boot)
        per.append({"label": lab, "t": best.t_mean, "t_min": best.t_min, "t_max": best.t_max, "n": best.n,
                    "period": best.period, "err_curv": best.period_err, "err_boot": sig,
                    "power": best.power, "half_power": best.half_period_power})
    tt = np.array([p["t"] for p in per])
    pp = np.array([p["period"] for p in per])
    ee = np.array([max(p["err_boot"], p["err_curv"]) for p in per])
    law, law_err, _, chi2 = V.fit_period_law(tt, pp, ee)
    # whole seasons too (the classic "one period per year" view)
    season_per = []
    for s in seas:
        r = V.best_period(t[s], m[s], e[s])
        season_per.append({"year": float(V.jd_to_year(r.t_mean)), "period": r.period, "err": r.period_err,
                           "power": r.power, "half_power": r.half_period_power, "n": r.n})
    # eclipse timing against the fitted law; refit the law to the phases
    t_ref = 2453500.0
    oc = V.oc_exponential(t, m, e, groups, tuple(law), t_ref)
    use = np.array([o["a2"] > 0.03 for o in oc])
    tph = np.array([o["t"] for o in oc])
    dph = np.array([o["dphase"] for o in oc])
    sph = np.array([o["dphase_err"] for o in oc])
    law_ph, law_ph_err, chi2_ph, dof_ph = V.fit_law_to_phases(tph[use], dph[use], sph[use], tuple(law), t_ref)
    # season light curves folded on the phase-fitted law (phase 0 = deeper minimum)
    p0, c, t0 = law_ph[:3]
    zero = law_ph[3]
    curves, shape = [], []
    for s, lab in zip(V.seasons(t), range(100)):
        yr = float(V.jd_to_year(t[s].mean()))
        if yr < 2002 or yr > 2008.9:
            continue
        if yr > 2008.0:
            s = s[t[s] < 2454550.0]  # 2008: before the eruption's rise
            if s.size < 20:
                continue
        ph = np.mod(V.cycles(t[s], p0, c, t0, t_ref) - zero, 1.0)
        # four-harmonic Fourier series plus a linear trend in time; the trend is removed for display
        tt = (t[s] - t[s].mean()) / max(1.0, np.ptp(t[s]))
        cols = [np.ones_like(ph), tt]
        for k in range(1, 5):
            cols += [np.cos(2 * np.pi * k * ph), np.sin(2 * np.pi * k * ph)]
        A = np.column_stack(cols)
        coef, *_ = np.linalg.lstsq(A / e[s][:, None], m[s] / e[s], rcond=None)
        detr = m[s] - coef[1] * tt
        zero_mag = np.median(detr)
        fph = np.linspace(0, 1, 201)
        fc = [np.ones_like(fph), np.zeros_like(fph)]
        for k in range(1, 5):
            fc += [np.cos(2 * np.pi * k * fph), np.sin(2 * np.pi * k * fph)]
        fit = np.column_stack(fc) @ coef - zero_mag
        curves.append({"year": yr, "phase": ph.tolist(), "dmag": (detr - zero_mag).tolist(),
                       "fit_phase": fph.tolist(), "fit_dmag": fit.tolist(), "n": int(s.size)})
    for o in oc:
        shape.append({"t": o["t"], "a1": o["a1"], "a2": o["a2"], "rms": o["rms"]})
    onset = V.eruption_onset(t, m)
    t_on10 = onset["t_on"] - onset["tau"] * np.log(10.0)  # excess flux reaches 10% of quiescence
    tyl = L.tylenda_periods()
    # our periods vs Tylenda's, point by point
    dp = pp - tyl["p"]
    tyl_err = 0.5 * (tyl["p_hi"] - tyl["p_lo"])
    comb = np.sqrt(ee**2 + tyl_err**2)
    np.savez(C.DERIVED_DIR / "v1309_curves.npz", allow_pickle=True, curves=np.array(curves, dtype=object))
    res = {
        "n_I": int(t.size), "n_pre_discovery": int(np.sum(t < L.TYLENDA2011["discovery_jd"])),
        "periods": per, "season_periods": season_per,
        "law_from_periods": {"p0": law[0], "c": law[1], "t0": law[2], "err": law_err.tolist(),
                             "t0_year": float(V.jd_to_year(law[2])), "t0_year_err": float(law_err[2] / C.YEAR_D),
                             "chi2": chi2, "dof": len(tt) - 3},
        "law_from_phases": {"p0": law_ph[0], "c": law_ph[1], "t0": law_ph[2], "zero": law_ph[3],
                            "err": law_ph_err.tolist(), "t0_year": float(V.jd_to_year(law_ph[2])),
                            "t0_year_err": float(law_ph_err[2] / C.YEAR_D), "chi2": chi2_ph, "dof": dof_ph},
        "oc": oc, "shape": shape,
        "onset": {"t_on": onset["t_on"], "tau_d": onset["tau"], "t_on10": float(t_on10),
                  "t_on10_year": float(V.jd_to_year(t_on10)), "t_on_year": float(V.jd_to_year(onset["t_on"])),
                  "err": onset["err"].tolist()},
        "period_drop_percent": float(100 * (pp[0] - pp[-1]) / pp[0]),
        "vs_tylenda": {"median_diff": float(np.median(dp)), "max_abs_diff": float(np.max(np.abs(dp))),
                       "chi2": float(np.sum((dp / comb) ** 2)), "n": int(dp.size),
                       "periods": tyl["p"].tolist(), "t": tyl["t"].tolist(), "err": tyl_err.tolist()},
    }
    return _f(res, 8)


# ---------------------------------------------------------------------------
# KIC 9832227
# ---------------------------------------------------------------------------


def kic_results() -> dict:
    import kic9832227 as K
    df, tpl = K.measure_all()
    out = {"n_timings": int(len(df)), "n_ours": int((df["kind"] == "ours").sum()),
           "years": [float(df["year"].min()), float(df["year"].max())]}
    nsvs = df[df["source"] == "NSVS"].iloc[0]
    nsvs_w = df[df["source"] == "NSVS (time as labelled)"].iloc[0]
    out["nsvs"] = {"oc_min": nsvs["oc_min"], "err_min": nsvs["err_min"], "oc_as_labelled_min": nsvs_w["oc_min"],
                   "shift_min": nsvs["oc_min"] - nsvs_w["oc_min"],
                   "expected_shift_min": (0.5 - K.P_REF) * 1440.0}
    # agreement with Kovacs et al. (2019) where both exist (same reference ephemeris)
    kov = {r["source"]: r for r in L.kovacs_table1()}
    pairs = []
    for _, r in df[df["kind"] == "ours"].iterrows():
        best = None
        for k, v in kov.items():
            if k.startswith(("Vulcan", "HAT", "WISE", "MLO", "ASAS-IV")):
                continue
            if abs(v["t_mid"] - r["t"]) < max(30.0, 0.5 * (r["t_max"] - r["t_min"])) and \
                    ("Kepler" in k) == (r["source"] == "Kepler") and ("NSVS" in k) == ("NSVS" in r["source"]):
                if "time as labelled" in r["source"]:
                    continue
                if best is None or abs(v["t_mid"] - r["t"]) < abs(best[1]["t_mid"] - r["t"]):
                    best = (k, v)
        if best:
            pairs.append({"ours": r["source"], "theirs": best[0], "year": r["year"], "d_min": r["oc_min"] - best[1]["oc"][3],
                          "err_min": float(np.hypot(r["err_min"], best[1]["eoc"][3]))})
    out["vs_kovacs"] = {"n": len(pairs), "median_diff_min": float(np.median([p["d_min"] for p in pairs])),
                        "rms_diff_min": float(np.std([p["d_min"] for p in pairs])), "pairs": pairs}

    def sel(upto=None, nsvs_label="NSVS", sources=None, drop_wasp2004=False, lit=True):
        d = df[df["source"] != ("NSVS (time as labelled)" if nsvs_label == "NSVS" else "NSVS")]
        if sources is not None:
            d = d[d["source"].isin([s for s in sources if s != "NSVS"] + [nsvs_label])]
        if drop_wasp2004:
            d = d[~((d["source"] == "SuperWASP") & (d["year"] < 2006))]
        if not lit:
            d = d[d["kind"] == "ours"]
        if upto:
            d = d[d["year"] <= upto]
        return d["t"].to_numpy(), d["oc_min"].to_numpy() / 1440, d["err_min"].to_numpy() / 1440

    grid = np.concatenate([np.arange(2016.8, 2030, 0.02), np.arange(2030, 2100, 0.25), np.arange(2100, 3000, 5.0)])
    profiles = {}
    for key, kw in (("molnar_like_labelled", dict(upto=2016.75, nsvs_label="NSVS (time as labelled)",
                                                    sources=K.MOLNAR_SOURCES, drop_wasp2004=True)),
                    ("molnar_like_correct", dict(upto=2016.75, sources=K.MOLNAR_SOURCES, drop_wasp2004=True)),
                    ("to2016_all_correct", dict(upto=2016.75)),
                    ("all_correct", dict()),
                    ("all_labelled", dict(nsvs_label="NSVS (time as labelled)"))):
        t, oc, e = sel(**kw)
        ty, nll, coef, jit = K.exponential_profile(t, oc, e, grid)
        j = int(np.nanargmin(nll))
        tau, p0 = K.exponential_tau(coef[j])
        A = np.column_stack([(t - K.T_REF), np.ones_like(t), (t - K.T_REF) ** 2])
        nq, cq, jq = K._wls_jitter(A, oc, e)
        ok = nll - nll[j] <= 0.5
        lo = float(ty[ok].min())
        hi = float(ty[ok].max())
        # one-sided 95% lower limit on t0 from the profile (Delta -lnL = 1.92)
        ok95 = nll - nll[j] <= 1.92
        try:
            lt = K.fit_oc(t, oc, e, "quad+ltte")
            ltte_jit = float(lt["jitter_d"] * 1440)
        except Exception:
            ltte_jit = float("nan")
        profiles[key] = {"n": int(t.size), "t0_best": float(ty[j]), "t0_lo": lo, "t0_hi": hi,
                         "jitter_grid": (jit * 1440).tolist(), "ltte_jitter_min": ltte_jit,
                         "t0_95_lo": float(ty[ok95].min()), "tau_d": float(tau), "p0": float(p0),
                         "jitter_min": float(jit[j] * 1440), "nll": float(nll[j]),
                         "quad_nll": float(nq), "quad_jitter_min": float(jq * 1440),
                         "delta_2lnL_exp_vs_quad": float(2 * (nq - nll[j])),
                         "grid": ty.tolist(), "profile": (nll - nll[j]).tolist()}
    out["exponential_profiles"] = profiles
    # full-data model comparison
    t, oc, e = sel()
    fits = {}
    for mname in ("quadratic", "cubic", "quad+ltte"):
        f = K.fit_oc(t, oc, e, mname)
        fits[mname] = {"params": f["params"].tolist(), "errors": f["errors"].tolist(), "jitter_min": f["jitter_d"] * 1440,
                       "nll": f["nll"], "bic": f["bic"], "k": f["k"]}
    fits["quadratic"]["pdot"] = K.pdot_from_quadratic(K.fit_oc(t, oc, e, "quadratic"))
    q3 = K.fit_oc(t, oc, e, "quad+ltte")
    fits["quad+ltte"]["pdot"] = float(2 * q3["params"][2] * K.P_REF / C.YEAR_D**2)
    fits["quad+ltte"]["amp_min"] = float(abs(q3["params"][3]) * 1440)
    fits["quad+ltte"]["period_yr"] = float(q3["params"][4])
    # light-travel-time interpretation: projected size of the pair's orbit and the mass function
    a12 = abs(q3["params"][3]) * C.DAY_S * 299792458.0 / C.AU_M
    p3 = q3["params"][4]
    fits["quad+ltte"]["a12sini_au"] = float(a12)
    fits["quad+ltte"]["mass_function_msun"] = float(a12**3 / p3**2)
    out["fits_all"] = fits
    # how well did the published models predict 2017-2025?
    late = df[(df["year"] > 2017.0) & (df["kind"] == "ours")]
    kov_pred = K.kovacs_model(late["t"].to_numpy())
    out["kovacs_prediction"] = {"max_miss_min": float(np.max(np.abs(late["oc_min"].to_numpy() - kov_pred))),
                                "miss_2024_min": float(np.interp(2024.6, late["year"], late["oc_min"] - kov_pred))}
    out["molnar"] = L.MOLNAR2017
    # period now vs Kepler era, from the quad+LTTE model's slope
    tt = np.array([K.jd_of(2011.0), K.jd_of(2025.0)])
    h = 5.0
    slope = (q3["f"](tt + h) - q3["f"](tt - h)) / (2 * h)
    out["period_change"] = {"p_2011": float(K.P_REF * (1 + slope[0])), "p_2025": float(K.P_REF * (1 + slope[1])),
                            "dp_s": float(K.P_REF * (slope[1] - slope[0]) * C.DAY_S)}
    return _f(out, 8)


def run_all(parts=("roche", "v1309", "kic")) -> dict:
    res = json.loads(C.RESULTS_PATH.read_text()) if C.RESULTS_PATH.exists() else {}
    funcs = {"roche": roche_results, "v1309": v1309_results, "kic": kic_results,
             "census": census_results, "pl": pl_results, "tcrb": tcrb_results}
    for p in parts:
        print(f"analysis: {p} ...", flush=True)
        res[p] = funcs[p]()
        C.RESULTS_PATH.write_text(json.dumps(res, indent=1))
    return res


def _gaia_sample():
    """ZTF binaries joined to their Gaia DR3 counterpart and AllWISE, with
    Galactic coordinates, dust-layer extinction and corrected parallaxes."""
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    import catalogs
    import plrelation as PL
    path = C.DERIVED_DIR / "ztf_gaia_sample.parquet"
    if path.exists():
        return pd.read_parquet(path)
    z = catalogs.read_ztf(types=C.ZTF_BINARY_TYPES)
    s = PL.build_sample(z, catalogs.read_xmatch(), catalogs.read_allwise())
    gal = SkyCoord(s["ra_gaia"].to_numpy() * u.deg, s["dec_gaia"].to_numpy() * u.deg).galactic
    s["l"], s["b"] = gal.l.deg, gal.b.deg
    s["av_layer"] = PL.dust_layer_av(s["dist_pc"].to_numpy(), s["b"].to_numpy())
    s.to_parquet(path)
    return s


def census_results() -> dict:
    import catalogs
    import census as CE
    z = catalogs.read_ztf()
    ew = CE.ew_sample(z)
    out = {"n_ztf": int(len(z)), "types": z["type"].value_counts().to_dict(), "n_ew": int(len(ew))}
    f = CE.fit_edge(ew["period"].to_numpy())
    mid, n = f["hist"]
    out["edge"] = {k: f[k] for k in ("p_edge", "logp_edge", "logp_edge_err", "width_dex", "width_err", "chi2", "dof")}
    out["quantiles"] = CE.quantiles(ew["period"].to_numpy())
    out["below"] = {"lt_0.22": int((ew["period"] < 0.22).sum()), "lt_0.20": int((ew["period"] < 0.20).sum()),
                    "lt_0.15": int((ew["period"] < 0.15).sum())}
    out["mag_slices"] = CE.mag_slices(ew)
    # the peak of the (magnitude-limited) log-P distribution
    h, e = np.histogram(ew["logp"], np.linspace(-0.9, 0.2, 221))
    out["peak_period_logp_hist"] = float(10 ** (0.5 * (e[np.argmax(h)] + e[np.argmax(h) + 1])))
    # Gaia distances: a volume-limited view, and what the sub-edge systems are
    s = _gaia_sample()
    sew = s[(s["type"] == "EW") & ~s["blended"] & (s["ruwe"] < 1.4)]
    near = sew[(sew["plx_snr"] > 5) & (sew["plx_corr"] > 1.0)]
    fn = CE.fit_edge(near["period"].to_numpy())
    out["edge_within_1kpc"] = {"n": int(len(near)), "p_edge": fn["p_edge"], "logp_edge_err": fn["logp_edge_err"],
                               "width_dex": fn["width_dex"]}
    hn, _ = np.histogram(near["logp"], e)
    out["peak_period_within_1kpc"] = float(10 ** (0.5 * (e[np.argmax(hn)] + e[np.argmax(hn) + 1])))
    good = sew[sew["plx_snr"] > 5].copy()
    good["M_G"] = good["g"] - 0.789 * good["av_layer"] - good["dm"]
    good["bp_rp0"] = good["bp_rp"] - 0.413 * good["av_layer"]
    sub = good[good["period"] < 0.21]
    norm = good[(good["period"] > 0.25) & (good["period"] < 0.30)]
    out["sub_edge"] = {"n_with_parallax": int(len(sub)), "median_M_G": float(np.median(sub["M_G"])),
                       "median_bp_rp0": float(np.median(sub["bp_rp0"])),
                       "frac_M_G_gt_5.5": float(np.mean(sub["M_G"] > 5.5)),
                       "frac_M_G_lt_3": float(np.mean(sub["M_G"] < 3.0)),
                       "ref_0.25_0.30_median_M_G": float(np.median(norm["M_G"])),
                       "ref_0.25_0.30_median_bp_rp0": float(np.median(norm["bp_rp0"]))}
    # Gaia DR3's own eclipsing-binary census (search floor at P = 0.2 d)
    v = catalogs.read_gaia_veb()
    ratio = v["DepthE2"] / v["DepthE1"]
    contact = (v["model_type"].str.contains("ELLIPSOIDAL") | ((v["sigPhaseGauss1"] > 0.08) &
               (v["sigPhaseGauss2"] > 0.08))) & (ratio > 0.7) & (v["period"] < 1.0)
    fg = CE.fit_edge(v.loc[contact, "period"].to_numpy(), logp_lo=-0.69)
    out["gaia"] = {"n_p_le_2d": int(len(v)), "p_min": float(v["period"].min()), "n_contact_like": int(contact.sum()),
                   "edge_contact_like": {"p_edge": fg["p_edge"], "logp_edge_err": fg["logp_edge_err"],
                                         "width_dex": fg["width_dex"]}}
    xm1 = s[["zid", "source_id", "period", "type"]]
    m = xm1[xm1["type"] == "EW"].merge(v[["source_id", "period"]], on="source_id", suffixes=("_ztf", "_gaia"))
    r = m["period_gaia"] / m["period_ztf"]
    out["ztf_vs_gaia"] = {"n_ew": int((s["type"] == "EW").sum()), "n_in_gaia_eb": int(len(m)),
                          "same": float(np.mean((r > 0.99) & (r < 1.01))),
                          "half": float(np.mean((r > 0.495) & (r < 0.505))),
                          "double": float(np.mean((r > 1.98) & (r < 2.02)))}
    np.savez(C.DERIVED_DIR / "census_hist.npz", edges=e, ztf=h, near=hn,
             gaia=np.histogram(np.log10(v.loc[contact, "period"]), e)[0],
             gaia_all=np.histogram(np.log10(v["period"]), e)[0])
    return _f(out)


def pl_results() -> dict:
    import plrelation as PL
    out = {"chen2018": PL.chen2018_recalibration()}
    s = _gaia_sample()
    lo, hi = L.CHEN2018_PL["logp_range"]
    base = (s["type"] == "EW") & ~s["blended"] & (s["ruwe"] < 1.4)
    w = s[base & s["w1_clean"] & (s["plx_snr"] > 10) & (s["dm_err"] < 0.11) & (s["logp"] > lo) & (s["logp"] < hi)].copy()
    w["M_W1"] = w["w1"] - PL.A_W1_OVER_AV * w["av_layer"] - w["dm"]
    w["M_G"] = w["g"] - PL.A_G_OVER_AV * w["av_layer"] - w["dm"]
    w["M_WG"] = w["w_g"] - w["dm"]
    e = np.hypot(w["dm_err"], w["w1_err"]).to_numpy()
    res = {}
    for band, col, err in (("W1", "M_W1", e), ("G", "M_G", w["dm_err"].to_numpy()), ("W_G", "M_WG", w["dm_err"].to_numpy())):
        for wt in (True, False):
            f = PL.fit_pl(w["logp"].to_numpy(), w[col].to_numpy(), err, logp0=0.0, weighted=wt)
            res[f"{band}_{'wls' if wt else 'ols'}"] = {k: f[k] for k in ("slope", "slope_err", "b0", "b0_err", "scatter", "n")}
    out["ztf_precise"] = res
    out["ztf_precise_median_dist_pc"] = float(w["dist_pc"].median())
    # the selection effect: the same fit in distance bins
    bins = []
    for d0, d1 in ((0, 500), (500, 800), (800, 1200), (1200, 2000), (2000, 4000)):
        k = (w["dist_pc"] >= d0) & (w["dist_pc"] < d1)
        if k.sum() < 50:
            continue
        f = PL.fit_pl(w["logp"][k].to_numpy(), w["M_W1"][k].to_numpy(), e[k.to_numpy()], logp0=0.0, weighted=False)
        bins.append({"d0": d0, "d1": d1, "n": int(k.sum()), "slope": f["slope"], "slope_err": f["slope_err"],
                     "b0": f["b0"], "scatter": f["scatter"]})
    out["distance_bins"] = bins
    # early-type (long-period) contact binaries
    et = s[base & s["w1_clean"] & (s["plx_snr"] > 10) & (s["dm_err"] < 0.11) & (s["logp"] > -0.25) & (s["logp"] < 0.1)]
    fe = PL.fit_pl(et["logp"].to_numpy(), (et["w1"] - PL.A_W1_OVER_AV * et["av_layer"] - et["dm"]).to_numpy(),
                   np.hypot(et["dm_err"], et["w1_err"]).to_numpy(), logp0=0.0, weighted=False)
    out["early_type"] = {k: fe[k] for k in ("slope", "slope_err", "b0", "scatter", "n")}
    # distance-indicator test: the Gaia-DR3 recalibration of Chen's bright stars,
    # applied to ZTF contact binaries, compared with their parallaxes (parallax space)
    cal = out["chen2018"]["W1_mean_dr3_ruwe"]
    tst = s[base & s["w1_clean"] & (s["plx_snr"] > 5) & (s["logp"] > lo) & (s["logp"] < hi)]
    a_w1 = PL.A_W1_OVER_AV * tst["av_layer"]
    pt = PL.parallax_space_test(tst["logp"].to_numpy(), tst["w1"].to_numpy(), tst["plx_corr"].to_numpy(),
                                tst["parallax_error"].to_numpy(), cal["slope"], cal["b0"], a_w1.to_numpy())
    out["parallax_test"] = {"n": int(len(tst)), "ratio": pt["ratio"], "ratio_err": pt["ratio_err"], "dmag": pt["dmag"]}
    # per-star distance precision from the calibrated relation
    out["distance_precision_percent"] = float(100 * 0.4605 * cal["scatter"])
    # Wesenheit distances for every ZTF contact binary with Gaia photometry, calibrated
    # on Chen's stars with DR3 parallaxes (no extinction needed)
    ch = pd.read_parquet(C.DERIVED_DIR / "chen2018_dr3.parquet")
    k = np.isfinite(ch["dm_dr3"]) & (ch["RUWE"] < 1.4) & np.isfinite(ch["BPmag"])
    wg = ch["g_dr3"] - PL.WESENHEIT_R * (ch["BPmag"] - ch["RPmag"]) - ch["dm_dr3"]
    fw = PL.fit_pl(ch["logp"][k].to_numpy(), wg[k].to_numpy(), np.zeros(int(k.sum())), logp0=0.0, weighted=False)
    out["wesenheit_calibration"] = {kk: fw[kk] for kk in ("slope", "slope_err", "b0", "b0_err", "scatter", "n")}
    allew = s[(s["type"] == "EW") & ~s["blended"] & np.isfinite(s["w_g"]) & (s["logp"] > lo) & (s["logp"] < hi)].copy()
    allew["dm_pl"] = allew["w_g"] - (fw["slope"] * allew["logp"] + fw["b0"])
    allew["d_pl_pc"] = 10 ** (0.2 * allew["dm_pl"] + 1)
    chk = allew[allew["plx_snr"] > 10]
    rr = chk["d_pl_pc"] / chk["dist_pc"]
    out["wesenheit_vs_parallax"] = {"n": int(len(chk)), "median_ratio": float(rr.median()),
                                    "scatter_percent": float(100 * 1.4826 * np.median(np.abs(rr - rr.median()))),
                                    "n_mapped": int(len(allew)), "median_d_pl": float(allew["d_pl_pc"].median()),
                                    "p90_d_pl": float(np.quantile(allew["d_pl_pc"], 0.9))}
    cols = ["zid", "ra", "dec", "l", "b", "logp", "period", "type", "w_g", "dm_pl", "d_pl_pc", "dist_pc", "plx_snr", "g",
            "bp_rp", "rmag", "gmag"]
    allew[cols].to_parquet(C.DERIVED_DIR / "pl_map.parquet")
    w[["logp", "M_W1", "M_G", "M_WG", "dist_pc", "bp_rp", "w1", "g", "av_layer", "dm_err"]].to_parquet(C.DERIVED_DIR / "pl_sample.parquet")
    return _f(out)


def tcrb_results() -> dict:
    import tcrb as TC
    df = TC.load_aavso(refresh=True)
    out = {"n_obs": int(len(df)), "n_bands": df["band"].value_counts().to_dict(),
           "first_jd": float(df["jd"].min()), "last_jd": float(df["jd"].max()),
           "first_year": float(df["year"].min()), "last_year": float(df["year"].max())}
    # cleaning statistics and the brightest rejected reports in the last year
    vis = TC.clean_daily(df, "Vis.")
    out["visual_rejected"] = int(vis["rejected"].sum())
    recent = vis[(vis["jd"] > df["jd"].max() - 365)]
    bad = recent[recent["rejected"]].sort_values("mag").head(5)
    out["recent_false_alarms"] = [{"jd": float(r.jd), "mag": float(r.mag), "obscode": r.obscode} for r in bad.itertuples()]
    # the last 30 days, by band: is it erupting?
    last = df[df["jd"] > df["jd"].max() - 30]
    out["last30"] = {b: {"median": float(g["mag"].median()), "n": int(len(g))} for b, g in last.groupby("band")
                     if len(g) >= 20}
    # eruptions: decline times after Schaefer's peaks (peak V = 2.0 assumed; AID misses the peaks)
    er = {}
    for yr in (1866, 1946):
        ec = TC.eruption_curve(df, yr, before=10, after=400)
        v = ec["Vis."]
        v = v[~v["rejected"]]
        d = TC.decline_times(v["days"].to_numpy()[v["days"] < 40], v["mag"].to_numpy()[v["days"] < 40],
                             peak_mag=L.SCHAEFER2023["peak_v"])
        sec = v[(v["days"] > 100) & (v["days"] < 260)]
        sb = TC.binned(sec.assign(rejected=False), 5.0, min_n=3)
        er[str(yr)] = {"n_vis": int(len(v)), "first_day": float(v["days"].min()), "t2": d["t_drop"][2.0],
                       "t3": d["t_drop"][3.0], "brightest_vis": float(v["mag"].min()),
                       "secondary_peak_day": float(sb["jd"].iloc[np.argmin(sb["mag"])] - L.SCHAEFER2023["peak_jd"][yr])
                       if len(sb) else None,
                       "secondary_peak_mag": float(sb["mag"].min()) if len(sb) else None,
                       "plateau_mag": float(np.median(v["mag"][(v["days"] > 40) & (v["days"] < 90)]))}
    out["eruptions"] = er
    # the decade before: visual brightness 1936-1946 vs 2015-2026 (common photometric system)
    vm = TC.seasonal_levels(df, "Vis.", 60.0)
    q1 = vm[(vm["year"] > 1925) & (vm["year"] < 1935)]["mag"].median()
    q2 = vm[(vm["year"] > 2005) & (vm["year"] < 2014)]["mag"].median()
    hs1 = vm[(vm["year"] > 1938) & (vm["year"] < 1944.5)]["mag"].median()
    hs2 = vm[(vm["year"] > 2016) & (vm["year"] < 2023)]["mag"].median()
    out["visual_levels"] = {"quiescence_1925_35": float(q1), "quiescence_2005_14": float(q2),
                            "high_1938_44": float(hs1), "high_2016_23": float(hs2)}
    for band in ("B", "V"):
        bm = TC.seasonal_levels(df, band, 30.0)
        out[f"{band}_levels"] = {
            "quiescence_2005_14": float(bm[(bm["year"] > 2005) & (bm["year"] < 2014)]["mag"].median()),
            "high_2016_23": float(bm[(bm["year"] > 2016) & (bm["year"] < 2023)]["mag"].median()),
            "dip_2023_24": float(bm[(bm["year"] > 2023.3) & (bm["year"] < 2024.3)]["mag"].max()),
            "recovery_2024_25": float(bm[(bm["year"] > 2024.4) & (bm["year"] < 2025.7)]["mag"].median()),
            "last_90d": float(bm[bm["jd"] > df["jd"].max() - 90]["mag"].median())}
    # ellipsoidal modulation in quiescence (V, 1990-2014) and in the high state
    for key, lo, hi in (("quiescence", 2447800.0, 2456700.0), ("high_state", 2457400.0, 2460000.0)):
        ph, resid, day = TC.ellipsoidal_fold(df, "V", lo, hi)
        if len(ph) < 100:
            continue
        coef, model = TC.fourier_fit(ph, resid, 2)
        a2 = float(np.hypot(coef[3], coef[4]))
        a1 = float(np.hypot(coef[1], coef[2]))

        def at(c, p):
            return c[0] + sum(c[2 * j - 1] * np.cos(2 * np.pi * j * p) + c[2 * j] * np.sin(2 * np.pi * j * p)
                              for j in (1, 2))
        # minimum at 0.75 (giant in front) minus minimum at 0.25 (inner, irradiated face in view);
        # bootstrap over days (optimistic: flickering correlates neighbouring days)
        rng = np.random.default_rng(C.SEED)
        boots = []
        for _ in range(300):
            kk = rng.integers(0, len(ph), len(ph))
            cb, _ = TC.fourier_fit(ph[kk], resid[kk], 2)
            boots.append(at(cb, 0.75) - at(cb, 0.25))
        out[f"ellipsoidal_{key}"] = {"n_days": int(len(ph)), "half_amp_2nd": a2, "half_amp_1st": a1,
                                     "full_amp_2nd": 2 * a2, "min_diff": float(at(coef, 0.75) - at(coef, 0.25)),
                                     "min_diff_err": float(np.std(boots))}
    # hazard from the recurrence record (as Pei et al. 2026)
    iv = L.PEI2026["intervals_yr"]
    last = 2000.0 + (L.SCHAEFER2023["peak_jd"][1946] - C.JD_J2000) / C.YEAR_D
    now = float(2000.0 + (df["jd"].max() - C.JD_J2000) / C.YEAR_D)
    hz = TC.recurrence_hazard(iv, last, now, 2027.0 - now)
    hz12 = TC.recurrence_hazard(iv, last, now, 1.0)
    # check: Pei et al.'s own numbers from 2026 July 11
    jul11 = 2000.0 + (2461232.5 - C.JD_J2000) / C.YEAR_D
    pei_check = TC.recurrence_hazard(iv, last, jul11, 2027.0 - jul11)["p"]
    pei_check12 = TC.recurrence_hazard(iv, last, jul11, 1.0)["p"]  # 2026 July 11 to 2027 July 11
    out["hazard"] = {"now_year": now, "mean_interval": hz["mean_interval"], "sd_interval": hz["sd_interval"],
                     "p_rest_of_2026": hz["p"], "p_next_12_months": hz12["p"],
                     "intervals": hz["intervals"].tolist(), "elapsed_yr": hz["elapsed"],
                     "check_rest_of_2026_from_jul11": pei_check, "check_12_months_from_jul11": pei_check12,
                     "pei2026_rest_of_2026_from_jul11": L.PEI2026["p_rest_of_2026"],
                     "pei2026_following_year": L.PEI2026["p_following_12_months"]}
    vm.to_parquet(C.DERIVED_DIR / "tcrb_visual_60d.parquet")
    return _f(out)
