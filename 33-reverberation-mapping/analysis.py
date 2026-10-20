"""Every measurement in the project, in one place.

``run()`` returns a dict with the raw products (light curves, lag
distributions, fits) for the figures and a JSON-friendly ``summary`` of every
headline number next to the published value it is checked against.
"""

from __future__ import annotations

import numpy as np

import config as C
import disk
import iccf
import lightcurves as LC
import literature as LIT
import masses as M
import ngc5548 as N
import physics as P
import regression as R
import rl
import sdssrm as S

Z1 = 1.0 + C.Z_NGC5548

#: The SDSS-RM lag shown in detail: RM305's Mg II, where the alias weighting
#: visibly rejects a half-year alias (32% of resamplings) and keeps the lag
#: Shen et al. also found (322 vs 311 days observed; their grade 4).
SHOWCASE = (305, "MG2")


def _round(x, nd=4):
    if isinstance(x, (float, np.floating)):
        return round(float(x), nd)
    if isinstance(x, (tuple, list)):
        return [_round(v, nd) for v in x]
    if isinstance(x, dict):
        return {k: _round(v, nd) for k, v in x.items()}
    if isinstance(x, (np.integer,)):
        return int(x)
    return x


def ours_rl(ours, line: str) -> dict | None:
    """R-L relation from the lags our own ICCF detected in SDSS-RM quasars."""
    t = rl.sdssrm_table().set_index("RMID")
    o = ours[(ours.line == line) & ours.ok & ours.detected].set_index("RMID")
    if len(o) < 8:
        return None
    cfg = rl.SDSSRM_LINE[line]
    L = t.loc[o.index, "L5100_agn"] if cfg["host"] else t.loc[o.index, cfg["lcol"]]
    z1 = 1 + t.loc[o.index, "z"]
    good = np.isfinite(np.log10(L.where(L > 0)))
    tau = (o.lag / z1)[good]
    x = np.log10(L[good]) - cfg["norm"]
    y = np.log10(tau)
    sy = R.log_err(tau, (o.lag_lo / z1)[good], (o.lag_hi / z1)[good])
    fit = R.linmix_fit(x, np.full(len(x), 1e-3), y, sy, seed=41)
    return dict(a=fit.a, b=fit.b, sig=fit.sig, n=fit.n, x=x.to_numpy(), y=y.to_numpy(), sy=np.asarray(sy))


def showcase(ours) -> dict:
    """Full ICCF products for one SDSS-RM lag (for the alias figure)."""
    rmid, line = SHOWCASE
    row = ours[(ours.RMID == rmid) & (ours.line == line)].iloc[0]
    cont = LC.sdssrm(rmid, "cont")
    ln = LC.sdssrm(rmid, S.SHORT[line])
    lags = iccf.lag_grid(np.floor(row.search_lo / 2) * 2, np.ceil(row.search_hi / 2) * 2, C.SDSSRM_LAG_STEP)
    res = iccf.fr_rss(cont.t, cont.y, cont.e, ln.t, ln.y, ln.e, lags, n_mc=1000, seed=rmid)
    w = iccf.weighted_lag(res, cont.t, cont.y, ln.t, C.ALIAS_K, C.ALIAS_GAP_DAYS, C.ALIAS_SMOOTH_DAYS)
    t = rl.sdssrm_table().set_index("RMID").loc[rmid]
    return dict(rmid=rmid, line=line, z=float(t.z), cont=cont, line_lc=ln, res=res, weighted=w,
                pub_pyroa=(t[f"{line}_lag"], t[f"{line}_lag_lo"], t[f"{line}_lag_hi"]),
                pub_iccf=(t[f"{line}_iccf"], t[f"{line}_iccf_lo"], t[f"{line}_iccf_hi"]))


def run(n_mc: int = C.N_MC, sdssrm_workers: int = 4) -> dict:
    out: dict = {}
    print("NGC 5548, AGN STORM 2014: lags ...")
    storm = N.storm_lags(n_mc)
    out["storm"] = storm
    hb_opt = storm["lags"][("Hb", "F5100", "full")]["cent"]
    out["storm_mass"] = N.storm_mass(hb_opt)
    out["storm_lum"] = N.storm_luminosity()
    print("NGC 5548: the BLR holiday ...")
    out["holiday"] = N.holiday(n_mc=500)
    print("NGC 5548: the accretion disk's own echo ...")
    dl = disk.disk_lags(n_mc=500)
    out["disk"] = dict(lags=dl, fits=disk.fits(dl))
    print("NGC 5548, 1988-2008: seasonal lags ...")
    watch = N.watch_seasons(n_mc)
    out["watch"] = watch
    print("R-L relations (LINMIX; a few minutes) ...")
    out["breathing"] = rl.ngc5548_breathing(watch["seasons"])
    out["bentz13"] = rl.bentz13_validation()
    out["local_rl"] = rl.local_rl()
    out["sdssrm_rl"] = rl.sdssrm_rl()
    out["combined_rl"] = rl.combined_hbeta(out["local_rl"]["sample"], out["sdssrm_rl"]["HB"]["sample"])
    print("SDSS-RM: our own ICCF lags ...")
    ours = S.run_all(workers=sdssrm_workers)
    out["sdssrm_ours"] = ours
    out["sdssrm_compare"] = S.compare(ours)
    out["sdssrm_fpr"] = S.false_positive_test(ours)
    out["sdssrm_ours_rl"] = {ln: ours_rl(ours, ln) for ln in ("HB", "MG2", "C4")}
    out["showcase"] = showcase(ours)
    print("Black hole masses ...")
    out["f_dyn"] = M.f_dynamical()
    out["f_msigma"] = M.f_msigma()
    out["virial"] = M.ngc5548_virial(hb_opt)
    sm = out["storm_mass"]
    out["ngc5548_msigma"] = M.ngc5548_on_msigma(sm["vp"], sm["vp_err"], out["f_msigma"], out["f_dyn"])
    out["local_masses"] = M.local_agn_masses(out["f_dyn"]["mean"][0])
    m_bh = out["ngc5548_msigma"]["ours_dyn"]["mbh"]
    out["disk"]["alpha_thin"] = {md: disk.thin_disk_alpha(m_bh, md) for md in (0.1, 1.0, 10.0)}
    out["disk"]["mbh"] = m_bh
    out["summary"] = summarize(out)
    return out


def summarize(o: dict) -> dict:
    s: dict = {}
    lags = o["storm"]["lags"]
    s["storm_lags"] = {f"{k[0]} vs {k[1]} ({k[2]})": dict(ours=v["cent"], r_max=v["r_max"],
                                                          published=LIT.PEI17_LAGS.get(
                                                              (k[0], k[1], "cent" if k[2] == "full" else "cent_" + k[2])))
                       for k, v in lags.items()}
    sm = o["storm_mass"]
    s["storm_mass"] = dict(tau=sm["tau"], tau_err=sm["tau_err"], R_au=sm["R_au"], vp=sm["vp"], vp_err=sm["vp_err"],
                           mbh_f447=sm["mbh"], mbh_f447_err=sm["mbh_err"], published_vp=LIT.PEI17_VP,
                           published_mbh=LIT.PEI17_MBH)
    s["storm_luminosity"] = o["storm_lum"]
    # The holiday
    h = o["holiday"]
    s["holiday"] = {k: dict(tau_obs=v["tau"], eta=v["eta"], flost=v["flost_red"], flost_err=v["flost_red_err"],
                            eta_anomaly=v.get("eta_red"),
                            published_tau=LIT.GOAD16_LAGS.get(k), published_eta=LIT.GOAD16_ETA.get(k, LIT.PEI17_HB_ETA),
                            published_flost=LIT.GOAD16_FLOST.get(k, LIT.PEI17_HB_FLOST))
                    for k, v in h.items() if not k.startswith("_")}
    # AGN Watch seasons
    seas = o["watch"]["seasons"]
    s["watch"] = [dict(year=round(x["year"], 1), jd=(x["jd_start"], x["jd_end"]), tau=x["tau"], published=x["pub_tau"],
                       tau_detrended=x["tau_detrended"], logL=x["logL"], published_logL=x["pub_logL"][0],
                       r_max=x["r_max"]) for x in seas]
    diffs = np.array([x["tau"][0] - x["pub_tau"][0] for x in seas])
    s["watch_agreement"] = dict(median_abs_diff_days=float(np.median(np.abs(diffs))),
                                n_within_1sig=int(np.sum([abs(x["tau"][0] - x["pub_tau"][0])
                                                          <= 0.5 * (x["pub_tau"][1] + x["pub_tau"][2]) for x in seas])),
                                n=len(seas))
    br = o["breathing"]
    s["breathing"] = dict(bentz07_sample=br["bentz07_sample"], bentz07_no_y12=br["bentz07_no_y12"],
                          all_watch=br["all_watch"], host_systematic=br["host_systematic"],
                          bentz07_sample_pub_errors=br["bentz07_sample_pub_errors"],
                          published_values=br["published_values"],
                          published=LIT.BENTZ07_NGC5548_SLOPE, published_no_y12=LIT.BENTZ07_NGC5548_SLOPE_NO_Y12)
    # How short was 2014 for NGC 5548's own R-L?
    f = br["bentz07_no_y12"]
    x14 = o["storm_lum"]["logL"] - 44.0
    pred14 = 10 ** (f["a"][0] + f["b"][0] * x14)
    s["storm_vs_own_rl"] = dict(logL_2014=o["storm_lum"]["logL"], predicted_tau=pred14, measured_tau=sm["tau"],
                                ratio=pred14 / sm["tau"])
    # R-L
    s["bentz13_validation"] = o["bentz13"]
    lr = o["local_rl"]
    s["local_rl"] = dict(all=lr["all"], per_object=lr["per_object"], fixed_0533=lr["fixed_0533"], n=lr["n"],
                         n_obj=lr["n_obj"], published=LIT.BENTZ13_RL)
    s["sdssrm_rl"] = {k: dict(fit=v["fit"], fixed=v["fixed"], n=v["n"], published=LIT.SHEN24_RL[rl.SDSSRM_LINE[k]["name"]])
                      for k, v in o["sdssrm_rl"].items() if k != "HB_no_host"}
    s["sdssrm_rl_hb_no_host"] = o["sdssrm_rl"]["HB_no_host"]
    s["sdssrm_offset_dex"] = o["sdssrm_rl"]["HB"]["fixed"]["a"][0] - lr["fixed_0533"]["a"][0]
    s["combined_rl"] = o["combined_rl"]
    s["sdssrm_compare"] = o["sdssrm_compare"]
    t78 = rl.sdssrm_table().set_index("RMID").loc[78]
    o78 = o["sdssrm_ours"][(o["sdssrm_ours"].RMID == 78) & (o["sdssrm_ours"].line == "HB")].iloc[0]
    s["rm078"] = dict(ours=(float(o78.lag), float(o78.lag_lo), float(o78.lag_hi)), pub_iccf=float(t78.HB_iccf),
                      pub_pyroa=float(t78.HB_lag))
    s["sdssrm_fpr"] = o["sdssrm_fpr"]
    s["sdssrm_ours_rl"] = {k: (None if v is None else dict(a=v["a"], b=v["b"], sig=v["sig"], n=v["n"]))
                           for k, v in o["sdssrm_ours_rl"].items()}
    sc = o["showcase"]
    s["showcase"] = dict(rmid=sc["rmid"], line=sc["line"], z=sc["z"], ours=(sc["weighted"]["lag"], *sc["weighted"]["err"]),
                         frac_rej=sc["weighted"]["frac_rej"], pub_iccf=sc["pub_iccf"], pub_pyroa=sc["pub_pyroa"])
    # Masses
    fd = o["f_dyn"]
    s["f_dynamical"] = dict(mean=fd["mean"], scatter=fd["scatter"], ml=fd["ml"], n=fd["n"], published=LIT.SHEN24_LOGF)
    fm = o["f_msigma"]
    s["f_msigma"] = dict(alpha=fm["alpha"], beta=fm["beta"], logf=fm["logf"], eps=fm["eps"], n_q=fm["n_q"],
                         n_agn=fm["n_agn"], quiescent_only=fm["quiescent_only"], published=LIT.WOO15_FIT)
    v = o["virial"]
    s["virial"] = dict(fit_all=v["fit_all"], fit_hb=v["fit_hb"], logvp_weighted=v["logvp_weighted"],
                       logvp_scatter=v["logvp_scatter"], n=len(v["sample"]), published_slope=LIT.BENTZ07_VIRIAL_SLOPE)
    s["ngc5548_msigma"] = o["ngc5548_msigma"]
    dk = o["disk"]
    s["disk"] = dict(lags={k: dict(lam=v["lam"], tau=v["tau"], published=v["published"]) for k, v in dk["lags"].items()},
                     fits=dk["fits"], alpha_thin=dk["alpha_thin"], mbh=dk["mbh"],
                     size_ratio=dk["fits"]["no_u_43"]["alpha"][0] / dk["alpha_thin"][0.1],
                     published_alpha43=disk.FAUSNAUGH16_ALPHA43, published_free=disk.FAUSNAUGH16_FREE)
    s["ngc5548_dynamical"] = dict(pancoast14=LIT.PANCOAST14_MBH)
    lm = o["local_masses"]
    s["local_masses"] = dict(n=len(lm), logM_min=float(lm.logM.min()), logM_max=float(lm.logM.max()))
    return _round(s)
