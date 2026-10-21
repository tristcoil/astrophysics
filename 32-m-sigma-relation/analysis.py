"""Every number in the README, computed from the cached catalogues.

    python analysis.py      # runs the analysis and prints the summary (run.py also draws)

``run()`` returns one dictionary: the catalogues with derived columns (for the
figures), every fit, and a JSON-friendly ``summary``.
"""

from __future__ import annotations

import json
import warnings

import numpy as np
import pandas as pd
from scipy import optimize, stats

import catalogs as K
import config as C
import fitting as F
import literature as LIT
import physics as P

warnings.filterwarnings("ignore", category=FutureWarning)


LOG_PIVOT = np.log10(C.SIGMA_PIVOT)


def _fit_dict(fit: F.LineFit, pub: dict | None = None) -> dict:
    d = fit.summary()
    if pub is not None:
        d["published"] = {k: pub[k] for k in ("a", "ea", "b", "eb", "eps") if k in pub}
        d["db_sigma"] = (fit.b - pub["b"]) / pub["eb"]
        d["da_sigma"] = (fit.a - pub["a"]) / pub["ea"]
    return d


def residuals(fit: F.LineFit, df: pd.DataFrame, col: str = "logsig", pivot: float = LOG_PIVOT):
    return df["logM"].to_numpy(float) - fit.predict(df[col].to_numpy(float) - pivot)


def _scatter_units(fit: F.LineFit, logm, elogm, logsig, elogsig):
    """Residual of one galaxy in units of the relation's total scatter at that point."""
    res = logm - fit.predict(logsig - LOG_PIVOT)
    tot = np.sqrt(fit.eps**2 + elogm**2 + (fit.b * elogsig) ** 2)
    return res, res / tot


# ---------------------------------------------------------------------------
# Catalogue checks: did we read the papers right?
# ---------------------------------------------------------------------------

def catalog_checks(cat) -> dict:
    out = {}
    v, viz = cat["vdb16"], K.load_vdb16_vizier("table2")
    # The GitHub file lists a few galaxies twice (two methods); match on name and method.
    vv = v[v["selected"] >= 0].copy()
    viz["method"] = viz["Meth"].map({"star": "stars", "stars": "stars", "gas": "gas", "CO": "gas",
                                     "maser": "maser", "reverb": "reverb"})
    m = vv.merge(viz[["name", "method", "logBHMass", "logsigma"]], on=["name", "method"], how="inner")
    out["vdb16_vs_vizier"] = {
        "github_rows_with_data": int(len(vv)), "vizier_rows": int(len(viz)), "matched": int(len(m)),
        "max_dlogM": float(np.nanmax(np.abs(m["logM"] - m["logBHMass"]))),
        "max_dlogsig": float(np.nanmax(np.abs(m["logsig"] - m["logsigma"]))),
    }
    kh = cat["kh13"]
    fit = kh[kh["kh13_fit"]]
    out["kh13"] = {
        "rows": int(len(kh)), "ellipticals": int((kh["table"] == 2).sum()),
        "classical_rows": int((kh["bulge"] == "classical").sum()),
        "pseudobulges": int((kh["bulge"] == "pseudo").sum()),
        "fit_sample": int(len(fit)), "fit_sample_sigma_below_270": int((fit["sigma"] < 270).sum()),
        "mean_logM_error_fit_sample": float(fit["elogM"].mean()),
        "published": {"fit_sample_sigma_below_270": 29, "mean_logM_error": C.KH13_EQUAL_LOGM_ERR},
    }
    mm = cat["mm13"]
    u = mm.drop_duplicates("name")
    out["mm13"] = {"rows": int(len(mm)), "galaxies": int(u["name"].nunique()),
                   "early": int(u["early_type"].sum()), "late": int((~u["early_type"]).sum()),
                   "published": {"galaxies": 72, "early": 53, "late": 19}}
    # Our H1 Sgr A* mass, if the H1 project's results are on disk.
    h1 = C.MILKY_WAY["h1_results"]
    if h1.exists():
        r = json.loads(h1.read_text())["headline"]["vs_gravity2022"]["gm6"]
        out["h1_mass"] = {"results_file": r["ours"] * 1e6, "config": C.MILKY_WAY["M"],
                          "agree": bool(abs(r["ours"] * 1e6 - C.MILKY_WAY["M"]) < 1e3)}
    else:
        out["h1_mass"] = {"results_file": None, "config": C.MILKY_WAY["M"], "agree": None}
    return out


# ---------------------------------------------------------------------------
# The M-sigma relation: our fits and the three published ones
# ---------------------------------------------------------------------------

def msigma_fits(cat) -> dict:
    v, kh, mm = cat["vdb16"], cat["kh13"], cat["mm13"]
    vs = v[v["selected"] == 1]
    fits = {}
    print("  van den Bosch sample (230 galaxies, 49 upper limits)")
    fits["vdb16"] = F.fit_msigma(vs, model="structural", keep_chain=True)
    fits["vdb16_conditional"] = F.fit_msigma(vs, model="conditional")
    fits["vdb16_two_gaussians"] = F.fit_msigma(vs, model="structural", n_gauss=2)
    fits["vdb16_detections"] = F.fit_msigma(vs[~vs["upper"]], model="structural")
    sharp = vs.copy()  # upper limits as hard limits at the top of their range
    sharp.loc[sharp["upper"], "logM"] = sharp.loc[sharp["upper"], "logM_limit"]
    sharp.loc[sharp["upper"], "elogM"] = 0.0
    fits["vdb16_sharp_limits"] = F.fit_msigma(sharp, model="structural")
    fits["vdb16_fitexy_detections"] = F.fit_msigma(vs, method="fitexy")
    print("  Kormendy & Ho sample (45 ellipticals and classical bulges)")
    ks = kh[kh["kh13_fit"]]
    fits["kh13"] = F.fit_msigma(ks, model="structural", keep_chain=True)
    fits["kh13_eq3"] = F.fit_msigma(ks, method="fitexy", equal_logm_err=C.KH13_EQUAL_LOGM_ERR)
    fits["kh13_eq5"] = F.fit_msigma(ks, method="fitexy")
    fits["kh13_eq12"] = F.fit_msigma(ks[ks["sigma"] < 270], method="fitexy",
                                     equal_logm_err=C.KH13_EQUAL_LOGM_ERR)
    print("  McConnell & Ma sample (72 galaxies)")
    fits["mm13"] = F.fit_msigma(mm, method="fitexy")
    fits["mm13_early"] = F.fit_msigma(mm[mm["early_type"]], method="fitexy")
    fits["mm13_late"] = F.fit_msigma(mm[~mm["early_type"]], method="fitexy")
    fits["mm13_structural"] = F.fit_msigma(mm, model="structural")
    return fits


def reproduction_table(fits) -> list[dict]:
    rows = [
        ("Kormendy & Ho 2013, eq. 3", "kh13_eq3", LIT.KH13["msigma_eq3"], "least squares, equal log M errors"),
        ("Kormendy & Ho 2013, eq. 5", "kh13_eq5", LIT.KH13["msigma_eq5"], "least squares, individual errors"),
        ("Kormendy & Ho 2013, eq. 12", "kh13_eq12", LIT.KH13["msigma_eq12"], "sigma_e < 270 km/s only"),
        ("McConnell & Ma 2013, all", "mm13", LIT.MM13["all"], "MPFITEXY"),
        ("McConnell & Ma 2013, early types", "mm13_early", LIT.MM13["early"], "MPFITEXY"),
        ("McConnell & Ma 2013, late types", "mm13_late", LIT.MM13["late"], "MPFITEXY"),
        ("van den Bosch 2016", "vdb16", LIT.VDB16["msigma"], "Bayesian, censored limits"),
    ]
    out = []
    for label, key, pub, how in rows:
        f = fits[key]
        out.append({"label": label, "key": key, "how": how, "n": f.n, "n_upper": f.n_upper,
                    "a": f.a, "ea": f.ea, "b": f.b, "eb": f.eb, "eps": f.eps,
                    "pub_a": pub["a"], "pub_ea": pub["ea"], "pub_b": pub["b"], "pub_eb": pub["eb"],
                    "pub_eps": pub["eps"], "db_sigma": (f.b - pub["b"]) / pub["eb"],
                    "da_sigma": (f.a - pub["a"]) / pub["ea"]})
    return out


# ---------------------------------------------------------------------------
# How each way of measuring a black hole sits on the relation
# ---------------------------------------------------------------------------

def method_offsets(cat, fit: F.LineFit) -> dict:
    v = cat["vdb16"]
    vs = v[(v["selected"] == 1) & ~v["upper"]].copy()
    vs["res"] = residuals(fit, vs)
    out = {}
    for meth, g in vs.groupby("method"):
        out[meth] = {"n": int(len(g)), "mean": float(g["res"].mean()),
                     "sem": float(g["res"].std(ddof=1) / np.sqrt(len(g))),
                     "rms": float(np.sqrt(np.mean(g["res"] ** 2))),
                     "median_logsig": float(g["logsig"].median())}
    return out


# ---------------------------------------------------------------------------
# Classical bulges, ellipticals and pseudobulges (Kormendy & Ho's split)
# ---------------------------------------------------------------------------

def bulge_types(cat, fits) -> dict:
    kh = cat["kh13"].copy()
    ref = fits["kh13"]
    kh["res"] = residuals(ref, kh)
    out = {"reference": "our fit to Kormendy & Ho's 45 ellipticals and classical bulges"}
    fit_sample = kh[kh["kh13_fit"]]
    for kind in ("elliptical", "classical"):
        g = fit_sample[fit_sample["bulge"] == kind]
        out[kind] = {"n": int(len(g)), "mean": float(g["res"].mean()),
                     "sem": float(g["res"].std(ddof=1) / np.sqrt(len(g))),
                     "rms": float(np.sqrt(np.mean(g["res"] ** 2)))}
    pb = kh[(kh["bulge"] == "pseudo")]
    det = pb[~pb["upper"]]
    out["pseudo"] = {"n": int(len(pb)), "n_detected": int(len(det)), "mean": float(det["res"].mean()),
                     "sem": float(det["res"].std(ddof=1) / np.sqrt(len(det))),
                     "rms_about_own_mean": float(det["res"].std(ddof=1)),
                     "frac_below": float((det["res"] < 0).mean()),
                     "median_mass_ratio_to_classical": float(10 ** np.median(det["res"]))}
    # A two-sample test: are pseudobulge residuals drawn from the classical ones?
    ks_stat = stats.ks_2samp(det["res"], fit_sample["res"])
    out["pseudo"]["ks_pvalue"] = float(ks_stat.pvalue)
    # A free fit to the pseudobulges: is there any correlation at all?
    pfit = F.fit_msigma(det, model="structural")
    fits["pseudo"] = pfit
    out["pseudo"]["free_fit"] = pfit.summary()
    out["pseudo"]["pearson_r"] = float(stats.pearsonr(det["logsig"], det["logM"]).statistic)
    out["pseudo"]["pearson_p"] = float(stats.pearsonr(det["logsig"], det["logM"]).pvalue)
    # Mergers in progress and the "monsters" (Kormendy & Ho's sections 6.4-6.5).
    special = {}
    for name in ("NGC 1316", "NGC 5128", "NGC 2960", "IC 1481", "NGC 4382", "NGC 1277", "NGC 4486B",
                 "NGC 4342", "NGC 3842", "NGC 4889"):
        r = kh[kh["name"] == name].iloc[0]
        special[name] = {"res": float(r["res"]), "sigma": float(r["sigma"]), "logM": float(r["logM"]),
                         "upper": bool(r["upper"]),
                         "kind": "merger" if r["merger"] else ("monster" if name in ("NGC 1277", "NGC 4486B")
                                                               else "excluded")}
    out["special"] = special
    return out


def anchors(fits) -> dict:
    """Where our two model-independent black holes sit."""
    out = {}
    mw, m87 = C.MILKY_WAY, C.M87
    for key, fit in (("kh13", fits["kh13"]), ("vdb16", fits["vdb16"]), ("pseudo", fits["pseudo"])):
        res, z = _scatter_units(fit, np.log10(mw["M"]), mw["M_err"] / mw["M"] / np.log(10),
                                np.log10(mw["sigma"]), mw["sigma_err"] / mw["sigma"] / np.log(10))
        out.setdefault("milky_way", {})[key] = {
            "res": float(res), "z": float(z),
            "predicted_M": float(10 ** fit.predict(np.log10(mw["sigma"]) - LOG_PIVOT))}
        esig = 0.5 * (np.log10(m87["sigma"] + m87["sigma_err_hi"]) - np.log10(m87["sigma"] - m87["sigma_err_lo"]))
        res, z = _scatter_units(fit, np.log10(m87["M"]), m87["M_err"] / m87["M"] / np.log(10),
                                np.log10(m87["sigma"]), esig)
        out.setdefault("m87", {})[key] = {
            "res": float(res), "z": float(z),
            "predicted_M": float(10 ** fit.predict(np.log10(m87["sigma"]) - LOG_PIVOT))}
    # The Milky Way's sigma is uncertain by 20 km/s; what sigma would put it on the classical line?
    f = fits["kh13"]
    out["milky_way"]["sigma_on_classical_line"] = float(
        C.SIGMA_PIVOT * 10 ** ((np.log10(mw["M"]) - f.a) / f.b))
    out["milky_way"]["mass_used"] = mw["M"]
    out["m87"]["history"] = LIT.M87_HISTORY
    return out


# ---------------------------------------------------------------------------
# Why do the published slopes differ? From McConnell & Ma's 5.6 to Kormendy & Ho's 4.4
# ---------------------------------------------------------------------------

def slope_budget(cat) -> dict:
    kh, mm = cat["kh13"], cat["mm13"]
    khu = kh.drop_duplicates("name", keep="first").set_index("name")  # NGC 3998: stellar mass first
    pseudo = set(kh.loc[kh["bulge"] == "pseudo", "name"])
    kept = set(kh.loc[kh["kh13_fit"], "name"])
    excluded = set(kh.loc[kh["bulge"].isin(["elliptical", "classical"]) & ~kh["kh13_fit"], "name"]) - kept

    def fx(d, **kw):
        f = F.fit_msigma(d, method="fitexy", **kw)
        return {"n": int(d["name"].nunique()), "b": f.b, "eb": f.eb, "a": f.a, "eps": f.eps}

    steps = [("McConnell & Ma's 72 galaxies, their values", fx(mm))]
    d = mm[~mm["name"].isin(pseudo)]
    steps.append(("drop the 17 pseudobulges", fx(d)))
    dropped = sorted(excluded & set(d["name"]))
    d = d[~d["name"].isin(excluded)]
    steps.append((f"drop {len(dropped)} more that Kormendy & Ho set aside", fx(d)))
    d = d.copy()
    d["logsig"], d["elogsig"] = d["name"].map(khu["logsig"]), d["name"].map(khu["elogsig"])
    steps.append(("use their σ_e values", fx(d)))
    d = d.drop_duplicates("name").copy()
    d["weight"] = 1.0
    d["logM"], d["elogM"] = d["name"].map(khu["logM"]), d["name"].map(khu["elogM"])
    steps.append(("use their black hole masses", fx(d)))
    full = kh[kh["kh13_fit"]]
    added = sorted(set(full["name"]) - set(d["name"]))
    steps.append((f"add {len(added)} galaxies new in 2013", fx(full)))
    steps.append(("equal log M errors (their eq. 3)", fx(full, equal_logm_err=C.KH13_EQUAL_LOGM_ERR)))
    # One change at a time, each applied alone to McConnell & Ma's sample.
    single = {}
    base = fx(mm)["b"]
    single["drop pseudobulges"] = fx(mm[~mm["name"].isin(pseudo)])["b"] - base
    single["drop the galaxies Kormendy & Ho set aside"] = fx(mm[~mm["name"].isin(excluded)])["b"] - base
    d = mm.copy()
    d["logsig"], d["elogsig"] = d["name"].map(khu["logsig"]), d["name"].map(khu["elogsig"])
    single["use their σ_e values"] = fx(d)["b"] - base
    d = mm.drop_duplicates("name").copy()
    d["weight"] = 1.0
    d["logM"], d["elogM"] = d["name"].map(khu["logM"]), d["name"].map(khu["elogM"])
    single["use their black hole masses"] = fx(d)["b"] - base
    # Leverage: drop one galaxy at a time from Kormendy & Ho's sample.
    jack = []
    for name in full["name"]:
        jf = F.fit_msigma(full[full["name"] != name], method="fitexy", equal_logm_err=C.KH13_EQUAL_LOGM_ERR)
        jack.append({"name": name, "b": jf.b, "sigma": float(full.loc[full["name"] == name, "sigma"].iloc[0])})
    jack = sorted(jack, key=lambda r: r["b"])
    return {"steps": [{"label": lab, **v} for lab, v in steps], "dropped": dropped, "added": added,
            "one_at_a_time": single, "jackknife": jack}


# ---------------------------------------------------------------------------
# Black hole mass and bulge mass: the ~0.5% ratio
# ---------------------------------------------------------------------------

def bulge_mass(cat, fits) -> dict:
    kh = cat["kh13"]
    ks = kh[kh["kh13_fit"]].copy()
    x = ks["logMbulge"].to_numpy() - np.log10(C.MASS_PIVOT)
    ex = ks["elogMbulge"].to_numpy()
    y = ks["logM"].to_numpy()
    eq = F.fit_fitexy(x, ex, y, np.full(len(ks), C.KH13_EQUAL_LOGM_ERR))
    eq.pivot = C.MASS_PIVOT
    ours = F.fit_bayes(x, ex, y, ks["elogM"].to_numpy(), model="structural")
    ours.pivot = C.MASS_PIVOT
    fits["mbulge_eq10"], fits["mbulge"] = eq, ours
    ratio = 10 ** (y - ks["logMbulge"].to_numpy()) * 100.0
    mm = cat["mm13"]
    md = mm[mm["Mbulge"].notna()].drop_duplicates("name")
    xm = np.log10(md["Mbulge"].to_numpy()) - np.log10(C.MASS_PIVOT)
    mmf = F.fit_fitexy(xm, np.full(len(md), 0.24), md["logM"].to_numpy(), md["elogM"].to_numpy())
    mmf.pivot = C.MASS_PIVOT
    fits["mbulge_mm13"] = mmf
    pb = kh[(kh["bulge"] == "pseudo") & ~kh["upper"]]
    pratio = 10 ** (pb["logM"] - pb["logMbulge"]) * 100.0
    mon = kh[kh["name"].isin(["NGC 1277", "NGC 4486B", "NGC 4342"])]
    return {
        "eq10": _fit_dict(eq, LIT.KH13["mbulge_eq10"]),
        "ours": ours.summary(),
        "ratio_at_1e11_percent_eq10": float(10 ** (eq.a - 11) * 100),
        "ratio_at_1e11_percent_ours": float(10 ** (ours.a - 11) * 100),
        "ratio_at_1e11_percent_ours_err": [float(10 ** (ours.extra["a_lo"] - 11) * 100),
                                           float(10 ** (ours.extra["a_hi"] - 11) * 100)],
        "ratio_median_percent": float(np.median(ratio)),
        "ratio_p16_p84_percent": [float(np.percentile(ratio, 16)), float(np.percentile(ratio, 84))],
        "ratio_range_percent": [float(ratio.min()), float(ratio.max())],
        "pseudo_ratio_median_percent": float(np.median(pratio)),
        "monsters_percent": {r["name"]: float(10 ** (r["logM"] - r["logMbulge"]) * 100) for _, r in mon.iterrows()},
        "mm13_dynamical": {**mmf.summary(), "published": LIT.MM13["mbulge_dynamical"],
                           "ratio_at_1e11_percent": float(10 ** (mmf.a - 11) * 100), "n": int(len(md))},
    }


# ---------------------------------------------------------------------------
# The sphere of influence: how small, and can we see it?
# ---------------------------------------------------------------------------

def sphere_of_influence(cat) -> dict:
    v = cat["vdb16"]
    vs = v[(v["selected"] == 1)].copy()
    M, sig = 10 ** vs["logM"], 10 ** vs["logsig"]
    vs["r_infl_pc"] = P.r_influence_pc(M, sig)
    vs["theta_infl"] = P.angular_size_arcsec(vs["r_infl_pc"], vs["dist"])
    vs["re_pc"] = 1e3 * 10 ** vs["logRe"]
    vs["rs_pc"] = P.schwarzschild_radius_pc(M)
    cat["vdb16_scales"] = vs
    det = vs[~vs["upper"]]
    by_method = {}
    for meth, g in det.groupby("method"):
        by_method[meth] = {"n": int(len(g)), "theta_median": float(g["theta_infl"].median()),
                           "frac_below_resolution": float((g["theta_infl"] < C.RESOLUTION_ARCSEC).mean())}
    mw_r = float(P.r_influence_pc(C.MILKY_WAY["M"], C.MILKY_WAY["sigma"]))
    m87_r = float(P.r_influence_pc(C.M87["M"], C.M87["sigma"]))
    return {
        "r_infl_pc_median": float(det["r_infl_pc"].median()),
        "r_infl_pc_p05_p95": [float(np.percentile(det["r_infl_pc"], 5)), float(np.percentile(det["r_infl_pc"], 95))],
        "r_infl_over_re_median": float(np.nanmedian(det["r_infl_pc"] / det["re_pc"])),
        "r_infl_over_re_p05_p95": [float(np.nanpercentile(det["r_infl_pc"] / det["re_pc"], 5)),
                                   float(np.nanpercentile(det["r_infl_pc"] / det["re_pc"], 95))],
        "rs_over_r_infl_median": float(np.median(det["rs_pc"] / det["r_infl_pc"])),
        "rs_over_re_median": float(np.nanmedian(det["rs_pc"] / det["re_pc"])),
        "volume_fraction_median": float(np.nanmedian((det["r_infl_pc"] / det["re_pc"]) ** 3)),
        "by_method": by_method,
        "milky_way": {"r_infl_pc": mw_r, "theta_arcsec": float(P.angular_size_arcsec(mw_r, 0.00793)),
                      "rs_au": float(P.schwarzschild_radius_pc(C.MILKY_WAY["M"]) * 206264.806)},
        "m87": {"r_infl_pc": m87_r, "theta_arcsec": float(P.angular_size_arcsec(m87_r, C.M87["D_mpc"])),
                "rs_au": float(P.schwarzschild_radius_pc(C.M87["M"]) * 206264.806)},
        "min_mass_at_20mpc_sigma100": float(P.min_measurable_mass(100.0, 20.0)),
        "min_mass_at_20mpc_sigma200": float(P.min_measurable_mass(200.0, 20.0)),
    }


# ---------------------------------------------------------------------------
# Are black hole hosts ordinary galaxies? (ATLAS3D test of selection bias)
# ---------------------------------------------------------------------------

def host_bias(cat, fits) -> dict:
    """Do galaxies with measured black holes have higher sigma than their peers at fixed L_K?

    Shankar et al. (2016) found that, at fixed stellar mass, black hole hosts
    have higher sigma than SDSS galaxies. ATLAS3D lets us repeat the test with
    one homogeneous data set: every early-type galaxy within 42 Mpc, with
    integral-field sigma_e and 2MASS M_K, hosts and non-hosts measured alike.
    """
    t = cat["atlas3d"].copy()
    v, kh = cat["vdb16"], cat["kh13"]
    measured = set(v.loc[(v["selected"] >= 0) & (v["method"].isin(["stars", "gas", "maser"])), "name"]) | set(kh["name"])
    detected = set(v.loc[(v["selected"] >= 0) & ~v["upper"] & v["method"].isin(["stars", "gas", "maser"]), "name"]) \
        | set(kh.loc[~kh["upper"], "name"])
    t["host"] = t["name"].isin(measured)
    t["host_detected"] = t["name"].isin(detected)
    ok = t["logsig"].notna() & t["logLk"].notna()
    t = t[ok].copy()
    # The comparison relation: sigma_e against L_K for the galaxies *without* a
    # black hole measurement, a quadratic in log L_K (Faber-Jackson bends).
    ref = t[~t["host"]]
    coef = np.polyfit(ref["logLk"] - 11.0, ref["logsig"], 2)
    t["dlogsig"] = t["logsig"] - np.polyval(coef, t["logLk"] - 11.0)
    rms_ref = float(t.loc[~t["host"], "dlogsig"].std(ddof=1))
    cat["atlas3d_test"] = t
    out = {"n_atlas3d": int(len(t)), "n_hosts": int(t["host"].sum()), "n_hosts_detected": int(t["host_detected"].sum()),
           "fj_coef": [float(c) for c in coef], "rms_nonhost": rms_ref}
    rng = np.random.default_rng(C.RNG_SEED)
    for key, mask in (("hosts", t["host"]), ("hosts_detected", t["host_detected"])):
        d = t.loc[mask, "dlogsig"].to_numpy()
        boot = np.array([rng.choice(d, len(d)).mean() for _ in range(20000)])
        # Permutation test: how often does a random set of N ATLAS3D galaxies do as well?
        allres = t["dlogsig"].to_numpy()
        perm = np.array([rng.choice(allres, len(d), replace=False).mean() for _ in range(20000)])
        out[key] = {"n": int(len(d)), "mean_dlogsig": float(d.mean()), "sem": float(d.std(ddof=1) / np.sqrt(len(d))),
                    "boot_68": [float(np.percentile(boot, 16)), float(np.percentile(boot, 84))],
                    "perm_p": float((perm >= d.mean()).mean()),
                    "median_dlogsig": float(np.median(d))}
    # Matched comparison in luminosity bins (no functional form assumed).
    bins = np.arange(9.8, 12.01, 0.3)
    rows = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        b = t[(t["logLk"] >= lo) & (t["logLk"] < hi)]
        h, n = b[b["host"]], b[~b["host"]]
        if len(h) >= 2 and len(n) >= 2:
            rows.append({"lo": lo, "hi": hi, "n_host": int(len(h)), "n_other": int(len(n)),
                         "host_median": float(h["logsig"].median()), "other_median": float(n["logsig"].median())})
    out["bins"] = rows
    # Host fraction rises with sigma at fixed L_K?
    d = out["hosts"]["mean_dlogsig"]
    out["implied_mass_bias_dex_at_fixed_LK"] = float(fits["vdb16"].b * d)
    out["implied_mass_bias_factor"] = float(10 ** (fits["vdb16"].b * d))
    return out


def resolution_simulation(cat, fits, n_real: int = 200) -> dict:
    """LABELLED SIMULATION: can 'resolve the sphere of influence' bias the fit?

    Fixed real inputs: the sigma and distance of each van den Bosch galaxy
    whose mass came from stars or gas (the resolution-limited methods).
    Simulated: black hole masses drawn from an assumed true relation; a mass
    counts as measured if its sphere of influence spans at least an angle
    theta_lim, otherwise the galaxy yields an upper limit at that threshold.
    We fit the simulated measurements the way one fits real data (detections
    only, and with censored limits) and compare with the truth.

    Real measurements have no sharp limit (17% of stellar-dynamical and 37%
    of gas-dynamical detections are formally below 0.1"), so the main
    threshold is the *effective* one that reproduces the real fraction of
    detections among these galaxies; 0.1" is kept as a stricter case.
    """
    v = cat["vdb16"]
    s = v[(v["selected"] == 1) & v["method"].isin(["stars", "gas"])]
    x = s["logsig"].to_numpy() - LOG_PIVOT
    ex = np.nan_to_num(s["elogsig"].to_numpy())
    sig = 10 ** s["logsig"].to_numpy()
    dist = s["dist"].to_numpy()
    ey = np.full(len(s), 0.1)
    rng = np.random.default_rng(C.RNG_SEED + 1)
    frac_obs = float((~s["upper"]).mean())
    truth = fits["vdb16"]

    def frac_resolved(theta, a, b, eps):
        m_min = np.log10(P.min_measurable_mass(sig, dist, theta))
        return float(np.mean(stats.norm.sf((m_min - a - b * x) / eps)))

    theta_eff = optimize.brentq(lambda t: frac_resolved(t, truth.a, truth.b, truth.eps) - frac_obs, 1e-4, 1.0)

    def ml(xx, exx, yy, eyy, up, start):
        nll = lambda t: -float(F.log_posterior(np.array(t), xx, exx, yy, eyy, up, np.ones_like(xx),  # noqa: E731
                                               "conditional", 0)[0])
        return optimize.minimize(nll, start, method="Nelder-Mead",
                                 options={"xatol": 1e-5, "fatol": 1e-7, "maxiter": 3000}).x

    out = {"n_galaxies": int(len(s)), "observed_detected_fraction": frac_obs, "theta_effective": theta_eff,
           "realisations": n_real}
    truths = {"our_fit": (truth.a, truth.b, truth.eps), "kh13_like": (8.5, 4.4, 0.3)}
    for t_label, (a_t, b_t, e_t) in truths.items():
        for th_label, theta in (("effective", theta_eff), ("strict_0p1", C.RESOLUTION_ARCSEC)):
            m_min = np.log10(P.min_measurable_mass(sig, dist, theta))
            det_fits, cen_fits, frac = [], [], []
            for _ in range(n_real):
                y_true = a_t + b_t * x + e_t * rng.standard_normal(len(x))
                y_obs = y_true + 0.1 * rng.standard_normal(len(x))
                resolved = y_true >= m_min
                frac.append(resolved.mean())
                det_fits.append(ml(x[resolved], ex[resolved], y_obs[resolved], ey[resolved],
                                   np.zeros(resolved.sum(), bool), [a_t, b_t, e_t]))
                cen_fits.append(ml(x, ex, np.where(resolved, y_obs, m_min), ey, ~resolved, [a_t, b_t, e_t]))
            det_fits, cen_fits = np.array(det_fits), np.array(cen_fits)
            out[f"{t_label}__{th_label}"] = {
                "truth": {"a": a_t, "b": b_t, "eps": e_t}, "theta": theta, "frac_resolved": float(np.mean(frac)),
                "detections_only": {"a": float(det_fits[:, 0].mean()), "b": float(det_fits[:, 1].mean()),
                                    "eps": float(det_fits[:, 2].mean()), "b_std": float(det_fits[:, 1].std())},
                "censored": {"a": float(cen_fits[:, 0].mean()), "b": float(cen_fits[:, 1].mean()),
                             "eps": float(cen_fits[:, 2].mean()), "b_std": float(cen_fits[:, 1].std())},
            }
    return out


# ---------------------------------------------------------------------------
# Active galaxies: the virial factor <f> and why it is circular
# ---------------------------------------------------------------------------

def agn_virial_factor(cat, fits) -> dict:
    """Reverberation mapping measures the 'virial product' c tau dV^2 / G; M = f x that.

    Nothing in the echo fixes f. It is calibrated by assuming that active
    galaxies follow the same M-sigma relation as quiet ones. We redo that
    calibration against relations fitted *without* any reverberation masses.
    """
    v, vp = cat["vdb16"], cat["agn_vp"]
    rv = v[(v["method"] == "reverb") & v["logsig"].notna()].drop_duplicates("name")
    m = rv[["name", "logsig", "elogsig", "selected"]].merge(
        vp[["name", "logM", "elogM", "z"]].rename(columns={"logM": "logVP", "elogM": "elogVP"}), on="name")
    m = m[m["logVP"].notna()].copy()
    cat["agn_matched"] = m
    quiet = v[(v["selected"] == 1) & (v["method"] != "reverb")]
    print("  quiescent relation without reverberation masses")
    fits["vdb16_no_reverb"] = F.fit_msigma(quiet, model="structural")
    out = {"n_agn": int(len(m))}

    def calibrate(fit: F.LineFit, label: str):
        # log f_i = log M_pred(sigma_i) - log VP_i; fit the mean with an intrinsic
        # scatter, propagating the sigma errors through the relation's slope.
        d = fit.predict(m["logsig"] - LOG_PIVOT) - m["logVP"]
        err2 = m["elogVP"] ** 2 + (fit.b * m["elogsig"]) ** 2

        def nll(p):
            mu, s = p
            if s <= 0:
                return np.inf
            var = err2 + s**2
            return 0.5 * np.sum((d - mu) ** 2 / var + np.log(var))

        mu, s = optimize.minimize(nll, [0.6, 0.3], method="Nelder-Mead").x
        var = err2 + s**2
        emu = float(np.sqrt(1.0 / np.sum(1.0 / var)))
        # The relation's own zero-point uncertainty moves every AGN together.
        emu_tot = float(np.hypot(emu, fit.ea))
        return {"relation": label, "log_f": float(mu), "elog_f": emu_tot, "f": float(10**mu),
                "f_lo": float(10 ** (mu - emu_tot)), "f_hi": float(10 ** (mu + emu_tot)),
                "scatter_dex": float(s)}

    out["calibrations"] = [
        calibrate(fits["vdb16_no_reverb"], "all quiescent galaxies (van den Bosch, no reverberation masses)"),
        calibrate(fits["kh13"], "ellipticals + classical bulges (Kormendy & Ho sample)"),
        calibrate(fits["pseudo"], "pseudobulges (Kormendy & Ho sample)"),
    ]
    out["database_values"] = {k: {"f": f, "ef": ef} for k, (f, ef) in LIT.F_FACTOR.items()}
    return out


# ---------------------------------------------------------------------------
# What it implies: feedback and energy
# ---------------------------------------------------------------------------

def feedback(cat, fits) -> dict:
    m_king_200 = float(P.king_msigma(C.SIGMA_PIVOT))
    m_king_200_h = float(P.king_msigma(C.SIGMA_PIVOT, kappa=P.kappa_es(1.0), f_gas=0.16))
    out = {"king_M_at_200": m_king_200, "king_M_at_200_pure_H_fg016": m_king_200_h,
           "kappa": float(P.kappa_es()), "f_gas": C.F_GAS_COSMIC}
    for key in ("kh13", "vdb16", "mm13"):
        f = fits[key]
        out[f"ratio_observed_to_king_{key}"] = float(10**f.a / m_king_200)
        out[f"slope_minus_4_sigma_{key}"] = float((f.b - 4.0) / f.eb)
    kh = cat["kh13"]
    ks = kh[kh["kh13_fit"]]
    er = P.energy_ratio(10 ** ks["logM"], 10 ** ks["logMbulge"], ks["sigma"])
    out["energy_ratio_median"] = float(np.median(er))
    out["energy_ratio_p16_p84"] = [float(np.percentile(er, 16)), float(np.percentile(er, 84))]
    out["energy_ratio_min"] = float(er.min())
    cat["energy_ratio"] = er
    return out


# ---------------------------------------------------------------------------
# Same galaxy, different numbers: how much do compilations disagree?
# ---------------------------------------------------------------------------

def compilation_systematics(cat) -> dict:
    kh = cat["kh13"].drop_duplicates("name", keep="first")
    mm = cat["mm13"].groupby("name", as_index=False).agg(
        logM=("logM", "mean"), logsig=("logsig", "first"), dist=("dist", "first"))
    v = cat["vdb16"]
    vd = v[v["selected"] >= 0].sort_values("selected", ascending=False).drop_duplicates("name")
    m = kh[["name", "logM", "logsig", "dist", "bulge"]].merge(
        mm, on="name", suffixes=("_kh", "_mm")).merge(
        vd[["name", "logM", "logsig", "dist"]].rename(columns={"logM": "logM_v", "logsig": "logsig_v", "dist": "dist_v"}),
        on="name", how="left")
    m["dsig_mm"] = m["logsig_mm"] - m["logsig_kh"]
    m["dsig_v"] = m["logsig_v"] - m["logsig_kh"]
    m["dM_mm"] = m["logM_mm"] - m["logM_kh"]
    m["dM_v"] = m["logM_v"] - m["logM_kh"]
    m["dD_mm"] = np.log10(m["dist_mm"] / m["dist_kh"])
    m["dD_v"] = np.log10(m["dist_v"] / m["dist_kh"])
    cat["compilations_matched"] = m

    def robust(s):
        s = s.dropna()
        return {"n": int(len(s)), "median": float(s.median()), "rms": float(np.sqrt(np.mean(s**2))),
                "frac_differ_gt_0p05": float((s.abs() > 0.05).mean()),
                "max_abs": float(s.abs().max())}

    # Dynamical masses scale with assumed distance: M ~ D. Check on the pairs
    # whose distances differ, masses corrected for that alone.
    dd = m[(m["dD_mm"].abs() > 0.01)]
    slope = float(np.polyfit(dd["dD_mm"], dd["dM_mm"], 1)[0]) if len(dd) > 3 else np.nan
    big_sig = m.loc[m["dsig_mm"].abs() > 0.1, ["name", "logsig_kh", "logsig_mm", "bulge"]]
    return {"sigma_mm_vs_kh": robust(m["dsig_mm"]), "sigma_vdb_vs_kh": robust(m["dsig_v"]),
            "mass_mm_vs_kh": robust(m["dM_mm"]), "mass_vdb_vs_kh": robust(m["dM_v"]),
            "distance_mm_vs_kh": robust(m["dD_mm"]), "dlogM_dlogD_slope": slope, "n_distance_pairs": int(len(dd)),
            "largest_sigma_disagreements": [
                {"name": r["name"], "sigma_kh": float(10 ** r["logsig_kh"]), "sigma_mm": float(10 ** r["logsig_mm"]),
                 "bulge": r["bulge"]} for _, r in big_sig.iterrows()]}


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run() -> dict:
    print("Loading catalogues")
    cat = K.load_all()
    res = {"cat": cat}
    res["checks"] = catalog_checks(cat)
    print("Fitting the M-sigma relation")
    fits = msigma_fits(cat)
    res["fits"] = fits
    res["reproduction"] = reproduction_table(fits)
    res["methods"] = method_offsets(cat, fits["vdb16"])
    print("Bulges and pseudobulges")
    res["bulges"] = bulge_types(cat, fits)
    res["anchors"] = anchors(fits)
    print("Slope budget")
    res["slope_budget"] = slope_budget(cat)
    print("Bulge masses")
    res["bulge_mass"] = bulge_mass(cat, fits)
    res["scales"] = sphere_of_influence(cat)
    print("Selection: ATLAS3D test and resolution simulation")
    res["host_bias"] = host_bias(cat, fits)
    res["resolution_sim"] = resolution_simulation(cat, fits)
    print("Active galaxies")
    res["agn"] = agn_virial_factor(cat, fits)
    res["feedback"] = feedback(cat, fits)
    res["systematics"] = compilation_systematics(cat)
    res["summary"] = summarise(res)
    return res


def summarise(res) -> dict:
    fits = res["fits"]
    s = {k: v for k, v in res.items() if k not in ("cat", "fits", "summary")}
    s["fits"] = {k: f.summary() for k, f in fits.items()}
    return json.loads(json.dumps(s, default=lambda o: o.item() if hasattr(o, "item") else str(o)))


if __name__ == "__main__":
    r = run()
    print(json.dumps(r["summary"]["reproduction"], indent=1))
