"""Ensemble analysis: from 104 sets of global parameters to the science.

``run_analysis()`` assumes the per-star pipeline (``stars.analyse_all``), the
Sun (``sun.analyse_sun``) and 16 Cyg (``mainseq.analyse_16cyg``) have run, and
then:

1. collects our numax, Dnu, epsilon and Delta Pi_1 and joins them with Yu et
   al. (2018), APOKASC-3 and Vrard et al. (2016), star by star;
2. compares the global parameters with both catalogues and with each of
   APOKASC-3's ten pipelines;
3. fits Dnu = alpha numax^beta to our stars and to Yu et al.'s 16,094;
4. turns numax and Dnu into masses and radii on three scales:
   * **self-calibrated**: our own measurement of the Sun as the reference,
     APOGEE temperatures, no corrections;
   * **Yu scale**: Yu et al.'s solar values (3090, 135.1 uHz; 5777 K) and
     temperatures, no corrections -> compared with their "no correction"
     masses and radii;
   * **APOKASC-3 scale**: their solar values (3076, 135.1416 uHz; 5772 K),
     their per-star corrections f_numax and f_Dnu, APOGEE temperatures ->
     compared with their recommended masses and radii;
5. tests the radii against radii from Gaia parallaxes (APOKASC-3's
   R_Gaia, from 2MASS K, Gaia DR3 parallaxes and APOGEE temperatures);
6. classifies hydrogen-shell versus helium-core burning from Delta Pi_1 and
   compares with the catalogues' evolutionary states;
7. places the stars on B1's Gaia HR diagram.
"""

from __future__ import annotations

import json
import warnings

import numpy as np
import pandas as pd

import catalogs as K
import config as C
import literature as LIT
import period_spacing as P
import sample as S
import scaling as SC
import stars

DPI1_CACHE = C.CACHE_DIR / "dpi1.json"
RNG_SEED = 48


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def robust(x) -> dict:
    """Median, robust scatter (1.4826 MAD), mean, std, N of finite values."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return {"n": 0}
    med = float(np.median(x))
    return {"median": med, "scatter": float(1.4826 * np.median(np.abs(x - med))),
            "mean": float(np.mean(x)), "std": float(np.std(x, ddof=1)) if len(x) > 1 else 0.0,
            "n": int(len(x))}


def boot_median_err(x, n=2000, seed=0) -> float:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    rng = np.random.default_rng(seed)
    return float(np.std([np.median(rng.choice(x, len(x))) for _ in range(n)]))


# ---------------------------------------------------------------------------
# 1. collect
# ---------------------------------------------------------------------------


def period_spacings(kics, *, force: bool = False) -> dict:
    """Delta Pi_1 for every star from its cached spectrum (cached as JSON)."""
    if DPI1_CACHE.exists() and not force:
        d = json.loads(DPI1_CACHE.read_text())
        if all(str(k) in d for k in kics):
            return {int(k): v for k, v in d.items()}
    import multiprocessing as mp

    with mp.get_context("fork").Pool(6) as pool:
        res = dict(pool.map(_dpi1_worker, [int(k) for k in kics]))
    DPI1_CACHE.write_text(json.dumps({str(k): v for k, v in res.items()}, indent=1, default=float))
    return res


def _dpi1_worker(kic):
    from threadpoolctl import threadpool_limits

    with threadpool_limits(1):
        js, arr = stars.load_star(kic)
        return kic, P.measure_for_star(arr, js)


def collect() -> pd.DataFrame:
    smp = S.load_or_build_sample(quiet=True).set_index("kic")
    dpi = period_spacings(smp.index)
    rows = []
    for kic in smp.index:
        js, arr = stars.load_star(int(kic))
        if not js.get("ok"):
            continue
        d = dpi[int(kic)]
        rows.append({
            "kic": int(kic), "numax": js["numax"], "numax_err": js["numax_err"],
            "dnu": js["dnu"], "dnu_err": js["dnu_err"], "dnu_acf": js["dnu_acf"],
            "dnu_acf_err": js["dnu_acf_err"], "epsilon": js["epsilon"],
            "dnu_radial": js["radial"].get("dnu_radial", np.nan),
            "fwhm": js["envelope_fwhm"], "amp": js["amp_per_radial_mode"],
            "gran": js["gran_at_numax"], "delta_lnL": js["delta_lnL_envelope"],
            "span_days": js["span_days"], "duty_cycle": js["duty_cycle"],
            "dpi1_attempted": bool(d.get("attempted", False)),
            "dpi1_detected": bool(d.get("detected", False)),
            "dpi1": d.get("dpi1", np.nan), "dpi1_err": d.get("dpi1_err", np.nan),
            "dpi1_q": d.get("q", np.nan), "dpi1_height": d.get("height", np.nan),
            "dpi1_hypothesis": d.get("hypothesis", ""),
        })
    df = pd.DataFrame(rows).set_index("kic")
    df = df.join(smp[["selection", "numax_bin", "kepmag", "n_quarters", "state_apokasc", "state_yu"]])
    y = K.yu2018()
    a = K.apokasc3()
    v = K.vrard2016()
    df = df.join(y[["numax", "numax_err", "dnu", "dnu_err", "amp", "width", "gran", "teff",
                    "teff_err", "mass_nocorr", "mass_nocorr_err", "radius_nocorr",
                    "radius_nocorr_err", "mass_rgbcorr", "mass_clcorr", "radius_rgbcorr",
                    "radius_clcorr"]].add_suffix("_yu"))
    df = df.join(a[["numax", "numax_err", "dnu", "dnu_err", "f_dnu", "f_dnu_err", "f_numax",
                    "mass", "mass_err", "radius", "radius_err", "teff", "teff_err", "feh",
                    "radius_gaia", "radius_gaia_err", "logg_seis", "category"]].add_suffix("_ak"))
    df = df.join(v[["dpi1", "dpi1_err", "q", "state", "dnu"]].add_suffix("_vrard"))
    return df


# ---------------------------------------------------------------------------
# 2. global parameters vs the catalogues
# ---------------------------------------------------------------------------


def compare_global(df: pd.DataFrame) -> dict:
    out = {}
    for cat in ("yu", "ak"):
        for q in ("numax", "dnu"):
            r = df[q] / df[f"{q}_{cat}"] - 1.0
            z = (df[q] - df[f"{q}_{cat}"]) / np.hypot(df[f"{q}_err"], df[f"{q}_err_{cat}"])
            out[f"{q}_vs_{cat}"] = {**robust(r), "median_err": boot_median_err(r),
                                    "normalised_scatter": robust(z)["scatter"]}
            for st in ("RGB", "RC"):
                m = df["state_apokasc"] == st
                out[f"{q}_vs_{cat}_{st}"] = robust(r[m])
    out["dnu_acf_vs_ak"] = robust(df["dnu_acf"] / df["dnu_ak"] - 1)
    out["dnu_acf_vs_yu"] = robust(df["dnu_acf"] / df["dnu_yu"] - 1)
    for st in ("RGB", "RC"):
        m = df["state_apokasc"] == st
        out[f"dnu_acf_vs_ak_{st}"] = robust(df["dnu_acf"][m] / df["dnu_ak"][m] - 1)
    # our formal uncertainties, for the README
    out["numax_err_frac"] = robust(df["numax_err"] / df["numax"])
    out["dnu_err_frac"] = robust(df["dnu_err"] / df["dnu"])
    # oscillation amplitude, granulation and envelope width vs Yu (SYD convention)
    out["amp_vs_yu"] = robust(df["amp"] / df["amp_yu"] - 1)
    out["gran_vs_yu"] = robust(df["gran"] / df["gran_yu"] - 1)
    out["width_vs_yu"] = robust(df["fwhm"] / df["width_yu"] - 1)
    # where do we sit among APOKASC-3's ten pipelines?
    t9 = K.apokasc3_pipelines()
    common = df.index.intersection(t9.index)
    pipes = {}
    for p in K.APOKASC_PIPELINES:
        col = f"Nmax-{p}"
        if col in t9:
            r = t9.loc[common, col] / df.loc[common, "numax_ak"] - 1
            pipes[p] = {"numax": robust(r)}
        col = f"DNu-{p}"
        if col in t9:
            r = t9.loc[common, col] / df.loc[common, "dnu_ak"] - 1
            pipes.setdefault(p, {})["dnu"] = robust(r)
    pipes["this work"] = {"numax": robust(df.loc[common, "numax"] / df.loc[common, "numax_ak"] - 1),
                          "dnu": robust(df.loc[common, "dnu"] / df.loc[common, "dnu_ak"] - 1)}
    out["pipelines"] = pipes
    return out


# ---------------------------------------------------------------------------
# 3. Dnu-numax relation
# ---------------------------------------------------------------------------


def _powerlaw_fit(numax, dnu, *, n_boot=1000, seed=0) -> dict:
    x = np.log10(np.asarray(numax))
    y = np.log10(np.asarray(dnu))
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    b, a = np.polyfit(x, y, 1)
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(n_boot):
        i = rng.integers(0, len(x), len(x))
        bs.append(np.polyfit(x[i], y[i], 1))
    bs = np.array(bs)
    resid = y - (a + b * x)
    return {"alpha": float(10 ** a), "alpha_err": float(np.std(10 ** bs[:, 1])),
            "beta": float(b), "beta_err": float(np.std(bs[:, 0])),
            "scatter_dex": float(np.std(resid)), "n": int(len(x))}


def dnu_numax_relation(df: pd.DataFrame) -> dict:
    y = K.yu2018()
    out = {"this work": _powerlaw_fit(df["numax"], df["dnu"], seed=1),
           "Yu et al. 2018 catalogue (our fit)": _powerlaw_fit(y["numax"], y["dnu"], seed=2)}
    out["literature"] = LIT.DNU_NUMAX
    return out


# ---------------------------------------------------------------------------
# 4-5. masses and radii on three scales; the Gaia test
# ---------------------------------------------------------------------------


def _mc(df_row, numax_samples, *, numax_sun, dnu_sun, teff_sun, teff, teff_err,
        f_numax=1.0, f_dnu=1.0, f_dnu_err=0.0, seed=0):
    rng = np.random.default_rng(seed)
    n = 4000
    dn = df_row["dnu"] + df_row["dnu_err"] * rng.standard_normal(n)
    return SC.monte_carlo(rng.choice(numax_samples, n), dn, teff, teff_err,
                          numax_sun=numax_sun, dnu_sun=dnu_sun, teff_sun=teff_sun,
                          f_numax=f_numax, f_dnu=f_dnu, f_dnu_err=f_dnu_err, seed=seed + 1)


def masses_radii(df: pd.DataFrame, sun: dict) -> tuple[pd.DataFrame, dict]:
    scales = {
        "self": dict(numax_sun=sun["numax_sun"], dnu_sun=sun["dnu_sun"],
                     teff_sun=C.TEFF_SUN_NOMINAL),
        "yu": dict(numax_sun=C.NUMAX_SUN_SYD, dnu_sun=C.DNU_SUN_SYD, teff_sun=C.TEFF_SUN_YU),
        "ak": dict(numax_sun=C.NUMAX_SUN_APOKASC, dnu_sun=C.DNU_SUN_APOKASC,
                   teff_sun=C.TEFF_SUN_NOMINAL),
    }
    res = {s: [] for s in scales}
    for kic, row in df.iterrows():
        _, arr = stars.load_star(int(kic))
        ns = arr["numax_samples"].astype(float)
        seed = int(kic) % 100000
        res["self"].append(_mc(row, ns, teff=row["teff_ak"], teff_err=row["teff_err_ak"],
                               seed=seed, **scales["self"]))
        res["yu"].append(_mc(row, ns, teff=row["teff_yu"], teff_err=row["teff_err_yu"],
                             seed=seed, **scales["yu"]))
        fd = row["f_dnu_ak"] if np.isfinite(row["f_dnu_ak"]) else 1.0
        fde = row["f_dnu_err_ak"] if np.isfinite(row["f_dnu_err_ak"]) else 0.0
        fn = row["f_numax_ak"] if np.isfinite(row["f_numax_ak"]) else 1.0
        res["ak"].append(_mc(row, ns, teff=row["teff_ak"], teff_err=row["teff_err_ak"],
                             f_numax=fn, f_dnu=fd, f_dnu_err=fde, seed=seed, **scales["ak"]))
    short = {"mass": "M", "radius": "R", "logg": "logg", "luminosity": "L"}
    for s, lst in res.items():
        t = pd.DataFrame(lst, index=df.index)
        for c in t.columns:
            base, err = (c[:-4], "_err") if c.endswith("_err") else (c, "")
            df[f"{short[base]}_{s}{err}"] = t[c]
    out = {"scales": {k: {kk: float(vv) for kk, vv in v.items()} for k, v in scales.items()}}
    # star-by-star comparisons: each catalogue on its own scale
    out["yu_scale_vs_yu"] = {"mass": robust(df["M_yu"] / df["mass_nocorr_yu"] - 1),
                             "radius": robust(df["R_yu"] / df["radius_nocorr_yu"] - 1)}
    out["ak_scale_vs_apokasc"] = {"mass": robust(df["M_ak"] / df["mass_ak"] - 1),
                                  "radius": robust(df["R_ak"] / df["radius_ak"] - 1)}
    for st in ("RGB", "RC"):
        m = df["state_apokasc"] == st
        out["ak_scale_vs_apokasc"][st] = {"mass": robust(df["M_ak"][m] / df["mass_ak"][m] - 1),
                                          "radius": robust(df["R_ak"][m] / df["radius_ak"][m] - 1)}
    out["self_vs_apokasc"] = {"mass": robust(df["M_self"] / df["mass_ak"] - 1),
                              "radius": robust(df["R_self"] / df["radius_ak"] - 1)}
    # how much each ingredient moves the masses (median factor)
    out["uncertainty_frac"] = {s: {"mass": robust(df[f"M_{s}_err"] / df[f"M_{s}"]),
                                   "radius": robust(df[f"R_{s}_err"] / df[f"R_{s}"])}
                               for s in scales}
    out["correction_factors"] = {"f_dnu": robust(df["f_dnu_ak"]), "f_numax": robust(df["f_numax_ak"]),
                                 "f_dnu_RGB": robust(df["f_dnu_ak"][df["state_apokasc"] == "RGB"]),
                                 "f_dnu_RC": robust(df["f_dnu_ak"][df["state_apokasc"] == "RC"])}
    return df, out


def gaia_test(df: pd.DataFrame) -> dict:
    """Seismic radius / Gaia radius, by scale, evolutionary state and numax."""
    out = {}
    for s in ("self", "yu", "ak"):
        r = df[f"R_{s}"] / df["radius_gaia_ak"]
        out[s] = {"all": {**robust(r), "median_err": boot_median_err(r)}}
        for st in ("RGB", "RC"):
            m = df["state_apokasc"] == st
            out[s][st] = {**robust(r[m]), "median_err": boot_median_err(r[m])}
            for lo, hi in ((0, 20), (20, 50), (50, 300)):
                mm = m & (df["numax"] >= lo) & (df["numax"] < hi)
                out[s][f"{st}_{lo}-{hi}"] = robust(r[mm])
    # the solar numax that would put our uncorrected radii on Gaia's scale,
    # for the low-RGB and red clump stars where APOKASC-3 calibrate (R ~ numax/numax_sun)
    for st in ("RGB", "RC"):
        m = (df["state_apokasc"] == st) & (df["numax"] > (50 if st == "RGB" else 0))
        ratio = np.median(df.loc[m, "R_self"] / df.loc[m, "radius_gaia_ak"])
        out[f"effective_numax_sun_{st}"] = {"value": float(ratio * df.attrs.get("numax_sun", np.nan)),
                                            "n": int(m.sum())}
    return out


# ---------------------------------------------------------------------------
# 6. evolutionary states
# ---------------------------------------------------------------------------


def evolutionary_states(df: pd.DataFrame) -> dict:
    att = df[df["dpi1_attempted"]]
    det = att[att["dpi1_detected"]]
    ours = np.where(det["dpi1_hypothesis"] == "RGB", "RGB", "HeB")
    out = {"n_attempted": int(len(att)), "n_detected": int(len(det)),
           "n_not_attempted": int((~df["dpi1_attempted"]).sum())}
    for name, col, mapping in (("APOKASC-3", "state_apokasc", {"RGB": "RGB", "RC": "HeB"}),
                               ("Yu et al. 2018", "state_yu", {"RGB": "RGB", "RC": "HeB"}),
                               ("Vrard et al. 2016", "state_vrard", {"RGB": "RGB", "RC": "HeB", "2CL": "HeB"})):
        theirs = det[col].map(mapping)
        known = theirs.notna().to_numpy()
        agree = (ours[known] == theirs[known].to_numpy())
        conf = pd.crosstab(pd.Series(ours[known], name="ours"),
                           pd.Series(theirs[known].to_numpy(), name=name))
        out[name] = {"n_compared": int(known.sum()), "n_agree": int(agree.sum()),
                     "confusion": conf.to_dict(),
                     "disagree": [int(k) for k in det.index[known][~agree]]}
    v = det[det["dpi1_vrard"].notna()]
    r = v["dpi1"] / v["dpi1_vrard"] - 1
    out["dpi1_vs_vrard"] = {**robust(r), "n_within_2pc": int((np.abs(r) < 0.02).sum()),
                            "n_within_5pc": int((np.abs(r) < 0.05).sum()),
                            "outliers": {int(k): [float(a), float(b)] for k, a, b in
                                         zip(v.index[np.abs(r) >= 0.05], v["dpi1"][np.abs(r) >= 0.05],
                                             v["dpi1_vrard"][np.abs(r) >= 0.05])}}
    vq = v[np.abs(r) < 0.05]
    out["q_vs_vrard"] = robust(vq["dpi1_q"] - vq["q_vrard"])
    # stars Vrard measured but we did not detect, and the reverse
    inv = att[att["dpi1_vrard"].notna()]
    out["vrard_not_detected_by_us"] = [int(k) for k in inv.index[~inv["dpi1_detected"]]]
    out["detected_not_in_vrard"] = [int(k) for k in det.index[det["dpi1_vrard"].isna()]]
    return out


# ---------------------------------------------------------------------------
# 7. HR diagram (B1)
# ---------------------------------------------------------------------------


def gaia_photometry(df: pd.DataFrame) -> pd.DataFrame:
    """Absolute G magnitude and dereddened BP-RP of the sample.

    Parallaxes get the Lindegren et al. (2021) zero-point correction
    (``gaiadr3-zeropoint`` package); reddening is the Schlafly & Finkbeiner
    (2011) E(B-V) at each position, converted with Casagrande & VandenBerg
    (2018)'s mean Gaia extinction coefficients (their Table 2: R_G = 2.740,
    R_BP = 3.374, R_RP = 2.035).
    """
    from astropy.table import Table
    from zero_point import zpt

    from fetch_data import EXTINCTION_PATH, GAIA_SAMPLE_PATH

    zpt.load_tables()
    g = Table.read(GAIA_SAMPLE_PATH).to_pandas().set_index("source_id")
    e = pd.read_csv(EXTINCTION_PATH).set_index("source_id")
    smp = S.load_or_build_sample(quiet=True).set_index("kic")
    sid = smp.loc[df.index, "gaia_dr3"].astype("int64")
    gg = g.loc[sid.to_numpy()]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        zp = zpt.get_zpt(gg["phot_g_mean_mag"].to_numpy(), gg["nu_eff_used_in_astrometry"].to_numpy(),
                         gg["pseudocolour"].to_numpy(), gg["ecl_lat"].to_numpy(),
                         gg["astrometric_params_solved"].to_numpy(), _warnings=False)
    out = pd.DataFrame(index=df.index)
    out["parallax"] = gg["parallax"].to_numpy()
    out["parallax_err"] = gg["parallax_error"].to_numpy()
    out["parallax_zp"] = zp
    out["parallax_corr"] = out["parallax"] - out["parallax_zp"]
    out["ebv"] = e.loc[sid.to_numpy(), "ebv_sandf"].to_numpy()
    out["a_g"] = 2.740 * out["ebv"]
    out["bp_rp0"] = gg["bp_rp"].to_numpy() - (3.374 - 2.035) * out["ebv"]
    out["g_abs0"] = (gg["phot_g_mean_mag"].to_numpy() - out["a_g"]
                     + 5 * np.log10(out["parallax_corr"]) - 10)
    out["distance_pc"] = 1000.0 / out["parallax_corr"]
    out["ruwe"] = gg["ruwe"].to_numpy()
    return out


def b1_background():
    """B1's diagram: parallax > 10 mas with B1's main quality cuts."""
    from astropy.table import Table

    from fetch_data import link_b1_sample

    t = Table.read(link_b1_sample(quiet=True))
    d = {c: np.asarray(t[c], float) for c in ("parallax", "parallax_over_error", "ruwe",
                                               "phot_g_mean_mag", "bp_rp")}
    ok = (d["parallax_over_error"] > 10) & (d["ruwe"] < 1.4) & np.isfinite(d["bp_rp"])
    g_abs = d["phot_g_mean_mag"] + 5 * np.log10(d["parallax"]) - 10
    return d["bp_rp"][ok], g_abs[ok]


# ---------------------------------------------------------------------------


def sun_cycle(sun: dict) -> dict:
    """Does the Sun's numax follow the activity cycle? (Howe et al. 2020: ~25 uHz.)

    Yearly numax against the yearly mean sunspot number (WDC-SILSO), for all
    well-covered years and for cycle 23 alone (1996-2008), before the drift
    of the later years.
    """
    from fetch_data import fetch_sunspots

    ssn = pd.read_csv(fetch_sunspots(quiet=True), sep=";", header=None)
    ssn = dict(zip(np.floor(ssn[0]).astype(int), ssn[1].astype(float)))
    out = {}
    for label, years in (("all years", sun["years_used"]),
                         ("cycle 23 (1996-2008)", [y for y in sun["years_used"] if int(y) <= 2008])):
        yy = [int(y) for y in years]
        x = np.array([ssn[y] for y in yy])
        nm = np.array([sun["years"][str(y)]["numax"] if str(y) in sun["years"] else sun["years"][y]["numax"]
                       for y in yy])
        b, a = np.polyfit(x, nm, 1)
        r = float(np.corrcoef(x, nm)[0, 1])
        rng = np.random.default_rng(5)
        bs = []
        for _ in range(2000):
            i = rng.integers(0, len(x), len(x))
            bs.append(np.polyfit(x[i], nm[i], 1)[0])
        out[label] = {"slope_uhz_per_ssn": float(b), "slope_err": float(np.std(bs)), "pearson_r": r,
                      "n_years": len(yy), "shift_0_to_170": float(170 * b), "shift_err": float(170 * np.std(bs))}
    return out


def run_analysis(*, quiet: bool = False) -> dict:
    import mainseq
    import sun as SUN

    sun = SUN.analyse_sun(quiet=True)
    df = collect()
    df.attrs["numax_sun"] = sun["numax_sun"]
    res = {"n_stars": int(len(df))}
    res["global"] = compare_global(df)
    res["dnu_numax"] = dnu_numax_relation(df)
    df, mr = masses_radii(df, sun)
    res["masses_radii"] = mr
    res["sun"] = {k: v for k, v in sun.items() if k != "years"}
    res["sun"]["years"] = {y: {k: sun["years"][y][k] for k in ("numax", "numax_err", "dnu", "dnu_err",
                                                               "duty_cycle", "envelope_fwhm",
                                                               "amp_per_radial_mode")}
                           for y in sun["years"]}
    res["sun_cycle"] = sun_cycle(sun)
    res["gaia_radii"] = gaia_test(df)
    res["states"] = evolutionary_states(df)
    phot = gaia_photometry(df)
    df = df.join(phot)
    res["hr"] = {"distance_pc": robust(df["distance_pc"]), "ebv": robust(df["ebv"])}
    res["sixteen_cyg"] = sixteen_cyg_summary(mainseq.analyse_16cyg(quiet=True), sun)
    df.to_csv(C.DATA_DIR / "results_stars.csv")
    C.RESULTS_PATH.write_text(json.dumps(res, indent=1, default=float))
    return res


def sixteen_cyg_summary(ms: dict, sun: dict) -> dict:
    t1, _ = K.lund2017()
    out = {}
    for name, r in ms.items():
        row = t1[t1["kic"] == str(r["kic"])].iloc[0]
        lit = LIT.SIXTEEN_CYG[name]
        teff, teff_err = float(row["teff"]), float(row["teff_err"])
        rng = np.random.default_rng(r["kic"] % 10000)
        n = 4000
        nm = r["numax"] + r["numax_err"] * rng.standard_normal(n)
        dn = r["dnu"] + r["dnu_err"] * rng.standard_normal(n)
        sc = {}
        for label, kw in (("self", dict(numax_sun=sun["numax_sun"], dnu_sun=sun["dnu_sun"])),
                          ("syd", dict(numax_sun=C.NUMAX_SUN_SYD, dnu_sun=C.DNU_SUN_SYD))):
            sc[label] = SC.monte_carlo(nm, dn, teff, teff_err, teff_sun=C.TEFF_SUN_NOMINAL,
                                       seed=1, **kw)
        out[name] = {"numax": r["numax"], "numax_err": r["numax_err"], "dnu": r["dnu"],
                     "dnu_err": r["dnu_err"], "dnu_acf": r["dnu_acf"], "epsilon": r["epsilon"],
                     "numax_lund": float(row["numax"]), "dnu_lund": float(row["dnu"]),
                     "teff": teff, "teff_err": teff_err, "scaling": sc, "literature": lit,
                     "span_days": r["span_days"], "duty_cycle": r["duty_cycle"]}
    return out
