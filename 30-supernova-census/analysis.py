"""Every measurement of the project, from the cached data, in one place.

``run_all()`` returns a dictionary of results (numbers, tables, arrays for the
figures) and writes the slow intermediate products to data/derived/:

* census: what the Bright Transient Survey caught, and how complete it is;
* rates: luminosity functions, volumetric rates and type fractions (rates.py),
  checked against Perley et al. (2020), Li et al. (2011) and other surveys;
* hosts: where each kind explodes, on the SFR-mass plane (hosts.py);
* light curves: Type II plateaus -> energies and ejecta masses (lightcurves.py);
* progenitors: the red supergiant problem, three ways (progenitors.py).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import bts
import config as C
import hosts as H
import lightcurves as LC
import literature as lit
import progenitors as P
import rates as R

LC_CACHE = C.DERIVED_DIR / "lightcurve_fits.csv"
HOSTS_CACHE = C.DERIVED_DIR / "hosts.csv"


def _f(x):
    return float(x) if x is not None and np.isfinite(x) else None


# ---------------------------------------------------------------------------
# Census
# ---------------------------------------------------------------------------

def census(df: pd.DataFrame) -> dict:
    stat = bts.statistical_sample(df)
    full = df[stat]
    sn = full[full["family"].notna()]
    out = {
        "n_candidates_all": int(len(df)),
        "n_classified_all": int(df["classified"].sum()),
        "n_sn_like_classified_all": int(df["family"].notna().sum()),
        "peak_date_range": [df["peak_date"].min().date().isoformat(),
                            df["peak_date"].max().date().isoformat()],
        "n_stat_sample": int(stat.sum()),
        "n_stat_classified": int(full["classified"].sum()),
        "frac_stat_classified": float(full["classified"].mean()),
        "stat_by_family": sn["family"].value_counts().to_dict(),
        "stat_by_type": full.loc[full["classified"], "tns_type"].value_counts().to_dict(),
        "quality_bit_pass": {b: float(df["q_" + b].mean()) for b in C.QUALCODE_BITS},
        "frac_quality": float(df["quality"].mean()),
        "frac_purity_given_quality": float(df.loc[df["quality"], "purity"].mean()),
    }
    # the non-supernova contaminants the purity cut and classifications remove
    out["stat_non_sn"] = full.loc[full["classified"] & full["family"].isna(), "tns_type"] \
        .value_counts().to_dict()
    return out


def distance_summary(df: pd.DataFrame) -> dict:
    sn = df[df["family"].notna()]
    near = sn[sn["z_cmb"] * C.C_KMS < C.CF4_MAX_CZ_KMS]
    cf4 = sn[sn["mu_source"] == "Cosmicflows-4 group"]
    d = (cf4["M_peak"] - cf4["M_bts"]).to_numpy()
    ex = cf4.assign(dM=cf4["M_peak"] - cf4["M_bts"]).sort_values("dM").head(8)
    return {"n_below_cz": int(len(near)), "n_cf4": int(len(cf4)),
            "dM_median": float(np.median(d)), "dM_abs_median": float(np.median(np.abs(d))),
            "dM_extremes": ex[["iau_id", "tns_type", "redshift", "M_bts", "M_peak"]]
            .round(3).to_dict("records"),
            "sn2023ixf": sn.loc[sn["iau_id"] == "SN2023ixf", ["M_bts", "M_peak", "mu"]]
            .round(3).to_dict("records")}


# ---------------------------------------------------------------------------
# Rates
# ---------------------------------------------------------------------------

def imf_number_per_mass(m1: float, m2: float, imf: str = "kroupa") -> float:
    """Stars born with m1 < M < m2 per solar mass of stars formed (0.08-100 M_sun
    for Kroupa 2001; 0.1-100 M_sun for Salpeter 1955)."""
    m = np.geomspace(0.08 if imf == "kroupa" else 0.1, 100, 200001)
    if imf == "kroupa":
        xi = np.where(m < 0.5, m ** -1.3, 0.5 * m ** -2.3)
    else:
        xi = m ** -2.35
    mass = np.trapezoid(xi * m, m)
    sel = (m >= m1) & (m <= m2)
    return float(np.trapezoid(xi[sel], m[sel]) / mass)


def rate_analysis(df: pd.DataFrame, rng) -> dict:
    out = {"perley_replication": R.perley_replication(df, rng)}
    fe = R.f_ext_from_dust()
    out["f_ext"] = fe
    s_all, info = R.our_sample(df)
    out["seasons"] = info["seasons"]
    out["seasonal_check"] = R.seasonal_check(s_all, info, rng, C.PERLEY_F_SKY, fe)
    out["n_candidates"], out["n_used"] = info["n_candidates"], info["n_used"]

    def window_rates(first_season: str, label: str):
        seas = info["seasons"][info["seasons"]["season"] >= first_season]
        expo = float((seas["T_eff_yr"] * seas["f_rec"]).sum())
        s = R.apply_weights(s_all[s_all["season"] >= first_season], exposure_yr=expo,
                            f_sky=C.PERLEY_F_SKY, f_ext=fe)
        res = {"label": label, "first_season": first_season, "n": int(len(s)),
               "T_eff_yr": float(seas["T_eff_yr"].sum()), "exposure_yr": expo,
               "f_rec_mean": expo / float(seas["T_eff_yr"].sum()),
               "rates": R.rates_by_class(s, rng),
               "ratio_cc_ia": R.rate_ratio(s, s["is_cc"] & (s["M_peak"] < C.RATE_MLIM_CC),
                                           (s["cls"] == "Ia") & (s["M_peak"] < C.RATE_MLIM_IA), rng),
               "fractions": {f"M<{m}": R.volumetric_fractions(s, rng, m) for m in (-14.0, -15.0, -16.0)}}
        return s, res

    s_main, out["main"] = window_rates("2019-06", "June 2019 - October 2025 (stable survey)")
    # the core-collapse rate tied to independent SN Ia rates through our ratio
    rat = out["main"]["ratio_cc_ia"]
    anchored = {}
    for key, (ria, eia) in (("PTF (Frohmaier+2019)", lit.FROHMAIER_2019_IA),
                            ("ASAS-SN (Desai+2026)", lit.DESAI_2026_IA)):
        e_rat = 0.5 * (rat["hi68"] - rat["lo68"]) / rat["ratio"]
        v = rat["ratio"] * ria
        anchored[key] = {"rate": v, "err68": v * np.hypot(e_rat, eia / ria)}
    out["cc_anchored"] = anchored
    _, out["full"] = window_rates("2018-06", "June 2018 - October 2025 (all seasons)")
    out["sample_main"] = s_main

    # Luminosity functions (main window)
    lfs = {}
    for fam in ("Ia", "II", "Ibc", "IIn", "SLSN", "TDE"):
        lfs[fam] = R.luminosity_function(s_main, s_main["family"] == fam, rng=rng)
    lfs["CC"] = R.luminosity_function(s_main, s_main["is_cc"], rng=rng)
    out["lf"] = lfs
    out["mag_limited"] = R.magnitude_limited_fractions(df)
    out["malmquist"] = R.malmquist_table(s_all)

    # How much the nearby distances matter: the same rates with redshift
    # distances (BTS convention magnitudes, cosmological V_max)
    s_z = s_main.copy()
    s_z["M_peak"] = s_z["M_bts"]
    s_z["vmax"] = R.vmax_cosmological(s_z["M_peak"].to_numpy())
    expo = out["main"]["exposure_yr"]
    s_z = R.apply_weights(s_z, exposure_yr=expo, f_sky=C.PERLEY_F_SKY, f_ext=fe)
    out["main_redshift_distances"] = {k: R.rates_by_class(s_z, rng)[k] for k in ("Ia", "CC", "II")}

    # Core-collapse rate against star formation
    cc = out["main"]["rates"]["CC"]
    z_cc = float(np.median(s_main.loc[s_main["is_cc"], "z_cmb"]))
    sfrd_md14 = lit.madau_dickinson_sfrd(z_cc) * lit.KROUPA_OVER_SALPETER_SFR
    sfrd_j3 = lit.J3_SFRD_KROUPA[0] * ((1 + z_cc) / 1.1) ** 2.7
    k = {}
    for lo, hi in ((8, 50), (8, 18), (8, 25), (8, 100)):
        k[f"{lo}-{hi}"] = imf_number_per_mass(lo, hi, "kroupa")
    out["cc_vs_sfr"] = {
        "z_median_cc": z_cc, "sfrd_md14_kroupa": sfrd_md14, "sfrd_j3_kroupa": sfrd_j3,
        "k_measured_md14": cc["rate"] * 1e-9 / sfrd_md14,
        "k_measured_j3": cc["rate"] * 1e-9 / sfrd_j3,
        "k_measured_err68": [cc["err68"][0] * 1e-9 / sfrd_md14, cc["err68"][1] * 1e-9 / sfrd_md14],
        "k_kroupa": k, "k_salpeter_8_50": imf_number_per_mass(8, 50, "salpeter"),
    }
    return out


# ---------------------------------------------------------------------------
# Hosts
# ---------------------------------------------------------------------------

def host_analysis(df: pd.DataFrame, rng) -> dict:
    h = H.assign_hosts(df)
    h.to_csv(HOSTS_CACHE, index=False)
    m = df.merge(h, on="ztf_id", how="inner")
    t = m["peak_date"]
    base = (m["family"].notna() & m["quality"] & m["purity"] & (m["peakmag"] <= 19.0)
            & (t >= C.STAT_PEAK_START) & (t <= pd.Timestamp(C.STAT_PEAK_END) + pd.Timedelta(days=1)))
    sel = m[base & m["has_mpa"] & (m["z_sn"] >= C.HOST_Z_MIN) & (m["z_sn"] <= C.HOST_Z_MAX)].copy()
    gal = H.load_background()
    gm, gs, gz, gss = (10 ** gal["lgm_tot_p50"].to_numpy(), 10 ** gal["sfr_tot_p50"].to_numpy(),
                       gal["z"].to_numpy(), gal["ssfr"].to_numpy())
    out = {"n_in_sdss_fields": int(df["in_sdss_field"].sum()),
           "n_with_host": int(len(h)), "n_host_with_spectrum": int(h["has_spec"].sum()),
           "frac_spec_z_agrees": float((h["z_agrees"] & h["has_spec"]).sum() / h["has_spec"].sum()),
           "n_with_mpa": int(h["has_mpa"].sum()), "n_plane": int(len(sel)),
           "n_background": int(len(gal))}
    quench = {}
    for cut in (-11.0, -12.0):
        rows = {}
        for fam in ("Ia", "II", "Ibc", "IIn"):
            x = sel[sel["family"] == fam]
            if len(x) < 5:
                continue
            hz = x["z_sn"].to_numpy()
            k = int((x["ssfr"] < cut).sum())
            q = gss < cut
            rows[fam] = {"n": int(len(x)), "k": k, "obs": k / len(x),
                         "ci68": H.binomial_interval(k, len(x)),
                         "pred_mass": float(np.sum(H.zmatched_weights(gz, gm, hz) * q)),
                         "pred_sfr": float(np.sum(H.zmatched_weights(gz, gs, hz) * q)),
                         "pred_number": float(np.sum(H.zmatched_weights(gz, np.ones_like(gm), hz) * q))}
        quench[str(cut)] = rows
    out["quenched"] = quench
    # Type Ia per unit stellar mass: star-forming vs quenched galaxies
    ia = quench["-11.0"]["Ia"]
    f_obs_q, f_mass_q = ia["obs"], ia["pred_mass"]
    ratio = ((1 - f_obs_q) / (1 - f_mass_q)) / (f_obs_q / f_mass_q)
    nq = ia["k"]
    nsf = ia["n"] - nq
    out["ia_per_mass_sf_over_q"] = {"ratio": ratio,
                                    "err": ratio * np.sqrt(1 / max(nq, 1) + 1 / max(nsf, 1))}
    out["a_plus_b"] = H.fit_a_plus_b(sel[sel["family"] == "Ia"], gal)
    out["per_mass_Ia"] = H.per_mass_rates(sel[sel["family"] == "Ia"], gal,
                                          out["a_plus_b"]["log_B_over_A_yr"])
    out["per_mass_CC"] = H.per_mass_rates(sel[sel["is_cc"]], gal)
    # convert B/A into A and B with our own SN Ia rate and J3's densities
    out["plane_sample"] = sel
    out["background_sample"] = gal.sample(n=min(len(gal), 60000), random_state=1)
    cc = sel[sel["is_cc"]].sort_values("ssfr")
    out["lowest_ssfr_cc"] = cc.head(8)[["ztf_id", "iau_id", "tns_type", "z_sn", "host_ra", "host_dec",
                                        "lgm_tot_p50", "sfr_tot_p50", "ssfr", "bptclass", "d_dlr",
                                        "ra", "dec"]]
    out["perley_table2"] = perley_table2(df)
    return out


def perley_table2(df: pd.DataFrame) -> dict:
    """Perley et al. (2020) Table 2 with today's sample: host colour-magnitude
    classes from the BTS host photometry, SDSS fields, 0.015 < z < 0.1,
    quality cut, m < 18.5."""
    t = df["peak_date"]
    s = df[df["in_sdss_field"] & df["quality"] & (df["peakmag"] <= C.M_LIM)
           & (df["redshift"] > 0.015) & (df["redshift"] < 0.1) & df["family"].notna()
           & (t >= C.STAT_PEAK_START) & (t <= pd.Timestamp(C.STAT_PEAK_END) + pd.Timedelta(days=1))]
    Mi, gi = s["hostabs"], s["hostcol"]
    tilt = np.maximum(0.85 - 0.05 * (Mi + 18), 0.85)
    cls = np.select([Mi.isna(), Mi > -16, gi > tilt, gi > 0.85],
                    ["hostless/unmeasured", "subdwarf", "red sequence", "green valley"],
                    default="blue cloud")
    out = {}
    for grp, mask in (("Ia", s["family"] == "Ia"), ("CC", s["is_cc"])):
        c = pd.Series(cls[mask.to_numpy()]).value_counts()
        n = int(mask.sum())
        out[grp] = {k: {"n": int(v), "frac": v / n} for k, v in c.items()}
        out[grp]["n"] = n
    return out


# ---------------------------------------------------------------------------
# Light curves
# ---------------------------------------------------------------------------

def lightcurve_analysis(df: pd.DataFrame, s_main: pd.DataFrame, refresh: bool = False) -> dict:
    ii = df[df["tns_type"].isin(C.LC_TYPE_II_CLASSES) & df["mu"].notna()]
    if LC_CACHE.exists() and not refresh:
        res = pd.read_csv(LC_CACHE)
    else:
        res = LC.analyse_sample(ii)
        res.to_csv(LC_CACHE, index=False)
    res, nirel = LC.finish_explosions(res)
    res = res.merge(df[["ztf_id", "iau_id", "tns_type", "z_cmb", "peakmag", "M_peak", "quality"]],
                    on="ztf_id", how="left")
    # 1/V_max weights (by peak magnitude) to undo the survey's bias toward
    # luminous events, for volume-representative medians
    res["vmax"] = R.vmax_cosmological(res["M_peak"].to_numpy())
    ok = res[res["plateau_ok"]].copy()

    def wmed(x, w):
        o = np.argsort(x)
        c = np.cumsum(w[o]) / w.sum()
        return float(x[o][np.searchsorted(c, 0.5)])

    w = 1 / ok["vmax"].to_numpy()
    tail = ok["M_ni_source"] == "tail"
    stats = {
        "n_lc": int(len(res)), "n_explosion_epoch": int((res["phase_ref"] == "explosion").sum()),
        "n_fall_seen": int(res["fit_accepted"].fillna(False).astype(bool).sum()),
        "n_full": int(len(ok)), "n_ni_tail": int(tail.sum()),
        "t_pt_median": float(ok["fit_t_pt"].median()), "t_pt_mean": float(ok["fit_t_pt"].mean()),
        "t_pt_std": float(ok["fit_t_pt"].std()),
        "t_pt_vweighted_median": wmed(ok["fit_t_pt"].to_numpy(), w),
        "t_end_mean": float((ok["fit_t_pt"] - 2 * ok["fit_w0"]).mean()),
        "t_end_std": float((ok["fit_t_pt"] - 2 * ok["fit_w0"]).std()),
        "Mr50_median": float(ok["M_r50"].median()),
        "Mr50_vweighted_median": wmed(ok["M_r50"].to_numpy(), w),
        "L50_median": float(ok["L50"].median()),
        "ni_tail_median": float(ok.loc[tail, "M_ni"].median()),
        "ni_relation": nirel,
        "E51_median": float(ok["E51"].median()), "E51_q": ok["E51"].quantile([0.16, 0.84]).tolist(),
        "Mej_median": float(ok["M_ej"].median()), "Mej_q": ok["M_ej"].quantile([0.16, 0.84]).tolist(),
        "Mej_median_Rlo": float(ok["M_ej_Rlo"].median()), "Mej_median_Rhi": float(ok["M_ej_Rhi"].median()),
        "E51_median_Rlo": float(ok["E51_Rlo"].median()), "E51_median_Rhi": float(ok["E51_Rhi"].median()),
        "frac_Mej_above_15": float((ok["M_ej"] > 15).mean()),
    }
    # longer plateaus are fainter? (Perley et al. found the same for durations)
    from scipy.stats import spearmanr
    rho, p = spearmanr(ok["fit_t_pt"], ok["M_r50"])
    stats["spearman_tpt_Mr50"] = {"rho": float(rho), "p": float(p)}
    b, a = np.polyfit(np.log10(ok.loc[tail, "L50"]), np.log10(ok.loc[tail, "M_ni"]), 1) \
        if tail.sum() > 5 else (np.nan, np.nan)
    stats["ni_vs_L50_slope"] = float(b)
    return {"stats": stats, "table": res}


# ---------------------------------------------------------------------------
# Progenitors
# ---------------------------------------------------------------------------

def progenitor_analysis(df: pd.DataFrame, rng) -> dict:
    s15 = P.load_smartt2015()
    db = P.load_davies_beasor2018()
    det = s15[(s15["kind"] == "detection") & s15["type"].isin(["II-P", "II-L"])]
    lim = s15[s15["kind"] == "limit"]
    out = {"smartt2015_table": s15, "db18_table": db,
           "n_s15_det": int(len(det)), "n_s15_lim": int(len(lim)),
           "counting": P.counting_argument(n_below=int(len(det) + len(lim)
                                                          + (s15["type"] == "IIb").sum() + 1))}
    out["ml_s15_stars"] = P.fit_mass_range(det, lim, "M_SG", "M_SG_up", "M_SG_lo", "M_SG")
    detk = det.assign(M_K_up=det["M_SG_up"], M_K_lo=det["M_SG_lo"])
    out["ml_s15_kepler"] = P.fit_mass_range(detk, lim, "M_K", "M_K_up", "M_K_lo", "M_K")
    dd, dl = db[db["kind"] == "detection"], db[db["kind"] == "limit"]
    out["ml_db18"] = P.fit_mass_range(dd, dl, "M_init", "M_up", "M_lo", "M_init")
    g = P.salpeter_gamma_L()
    out["gamma_L"] = g
    out["lum_fit"] = P.fit_luminosity_range(db, g, 0.0)
    out["lum_fit_beasor"] = P.fit_luminosity_range(db, g, C.BEASOR24_EXTRA_ERR_DEX)
    out["brightest"] = P.brightest_of_n(db, g, rng=rng)
    out["brightest_beasor"] = P.brightest_of_n(db, g, extra_single_band=C.BEASOR24_EXTRA_ERR_DEX, rng=rng)
    out["nothing_brighter"] = P.nothing_brighter(db, g, rng=rng)
    out["nothing_brighter_beasor"] = P.nothing_brighter(db, g, extra_single_band=C.BEASOR24_EXTRA_ERR_DEX,
                                                        rng=rng)
    out["nothing_brighter_25msun"] = P.nothing_brighter(db, g, xhi=float(np.log10(
        10 ** ((np.log10(25.0) - C.DB18_MLR[0]) / C.DB18_MLR[1]))), rng=rng)
    # The "missing" stars as a probability, against the assumed upper mass:
    # with N progenitors drawn from Salpeter between 8 M_sun and M_max, the
    # chance that none was born above 18 M_sun
    n = int(len(dd) + len(dl))
    mmax = np.arange(18.0, 40.01, 0.5)
    f_above = np.array([P.salpeter_number_fraction(18.0, m, 8.0, m) if m > 18 else 0.0 for m in mmax])
    out["p_none_above18"] = {"N": n, "M_max": mmax, "expected": n * f_above,
                             "p_none": (1 - f_above) ** n}
    # SNe with progenitor studies that the BTS also caught
    names = list(P.RECENT_SED_PROGENITORS["sn"]) + ["SN2018aoq", "SN2020jfo", "SN2021yja",
                                                   "SN2022acko", "SN2024abfl", "SN2025pht",
                                                   "SN2019yvr"]
    out["in_bts"] = df[df["iau_id"].isin(names)][["iau_id", "ztf_id", "tns_type", "peakmag",
                                                    "M_peak", "mu_source"]]
    out["recent"] = P.RECENT_SED_PROGENITORS
    return out


# ---------------------------------------------------------------------------

def run_all(refresh_lc: bool = False) -> tuple[dict, pd.DataFrame]:
    rng = np.random.default_rng(C.RNG_SEED)
    df = bts.load()
    res = {"bts_check": bts.check_bts_convention(df), "census": census(df),
           "distances": distance_summary(df)}
    res["rates"] = rate_analysis(df, rng)
    res["hosts"] = host_analysis(df, rng)
    res["lightcurves"] = lightcurve_analysis(df, res["rates"]["sample_main"], refresh_lc)
    res["progenitors"] = progenitor_analysis(df, rng)
    return res, df
