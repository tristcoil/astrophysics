"""NGC 5548: lags, mass and the 'BLR holiday' in 2014; 20 years of breathing.

Three analyses of the best-studied active galaxy:

1. **AGN STORM 2014** (Pei et al. 2017): H-beta and He II against the optical
   (5100 A) and far-UV (1158 and 1367 A, HST) continua, for the full campaign
   and its two halves, with Pei et al.'s ICCF settings.
2. **The BLR holiday** (Goad et al. 2016; Pei et al. 2017): for about two months
   the broad lines stopped following the continuum. We model each line's flux
   as a power of the lagged continuum flux, fit it where the two behave
   normally, and measure how much line flux went missing during the anomaly.
3. **AGN Watch, 1988-2008** (Peterson et al. 2002; 2013): the H-beta lag season
   by season from the 20-year light curves, against the AGN luminosity of that
   season: one black hole whose broad-line region grows and shrinks.
"""

from __future__ import annotations

import numpy as np

import agnmass as A
import config as C
import iccf
import lightcurves as LC
import literature as LIT
import physics as P

Z1 = 1.0 + C.Z_NGC5548


def _lag(cont: LC.LightCurve, line: LC.LightCurve, lags, n_mc, seed, detrend=True):
    yc = iccf.detrend(cont.t, cont.y) if detrend else cont.y
    yl = iccf.detrend(line.t, line.y) if detrend else line.y
    return iccf.fr_rss(cont.t, yc, cont.e, line.t, yl, line.e, lags, n_mc=n_mc, seed=seed)


def _rest(res: iccf.LagResult) -> dict:
    """Rest-frame summary of an FR/RSS run."""
    m, lo, hi = res.cent
    pm, plo, phi = res.peak
    return dict(cent=(m / Z1, lo / Z1, hi / Z1), peak=(pm / Z1, plo / Z1, phi / Z1),
                r_max=res.observed.r_max, success=res.success, n_mc=res.n_mc,
                cent_obs_curve=res.observed.tau_cent / Z1)


# ---------------------------------------------------------------------------
# 1. AGN STORM 2014
# ---------------------------------------------------------------------------


def storm_lightcurves() -> dict:
    s = LC.storm_spectroscopy("MDM", truncate=True)
    return {
        "F5100": s["F5100"], "Hb": s["HbSD"], "HeII": s["HeIISD"], "Hb_linear": s["Hb"],
        "F1158": LC.storm_hst_continuum(1157.5), "F1367": LC.storm_hst_continuum(1367.0),
    }


def storm_lags(n_mc: int = C.N_MC) -> dict:
    """Pei et al. 2017's lag table, recomputed: ICCF on -20..40 d (0.25 d steps),
    linear detrending, 1000 FR/RSS realizations, MDM spectra only (133 epochs)."""
    lcs = storm_lightcurves()
    lags = iccf.lag_grid(*C.STORM_LAG_RANGE, C.STORM_LAG_STEP)
    pairs = [("Hb", "F5100"), ("Hb", "F1158"), ("Hb", "F1367"), ("HeII", "F5100"), ("HeII", "F1158"),
             ("F5100", "F1158")]
    out, raw = {}, {}
    for k, (line, cont) in enumerate(pairs):
        for seg, (lo, hi) in {"full": (-np.inf, np.inf), "T1": (-np.inf, C.STORM_SPLIT_THJD),
                              "T2": (C.STORM_SPLIT_THJD, np.inf)}.items():
            c = lcs[cont].window(lo, hi)
            ln = lcs[line].window(lo, hi)
            res = _lag(c, ln, lags, n_mc, seed=100 + 7 * k + len(seg))
            out[(line, cont, seg)] = _rest(res)
            raw[(line, cont, seg)] = res
    # Without detrending, for the record (Pei et al. say detrending moved lags by ~0.01 d).
    res = _lag(lcs["F5100"], lcs["Hb"], lags, n_mc, seed=5, detrend=False)
    out[("Hb", "F5100", "full_nodetrend")] = _rest(res)
    # H-beta measured with a straight-line continuum instead of spectral decomposition.
    res = _lag(lcs["F5100"], lcs["Hb_linear"], lags, n_mc, seed=6)
    out[("Hb_linear", "F5100", "full")] = _rest(res)
    return {"lags": out, "raw": raw, "lcs": lcs}


def storm_mass(lag_rest: tuple, sigma=LIT.PEI17_SIGMA_RMS, log_f=C.LOG_F_WOO15) -> dict:
    """Virial product and mass from a rest-frame lag (value, -, +) and sigma_line(rms).

    The line width is the one ingredient taken from Pei et al.: it needs the
    nightly spectra, which are not public.
    """
    tau, lo, hi = lag_rest
    tau_err = 0.5 * (lo + hi)
    vp = float(P.virial_product(tau, sigma[0]))
    vp_err = float(P.virial_product_err(tau, tau_err, sigma[0], sigma[1]))
    f = 10**log_f
    return dict(tau=tau, tau_err=tau_err, R_lt_days=tau, R_au=tau * C.LIGHT_DAY_CM / C.AU_CM,
                sigma=sigma, vp=vp, vp_err=vp_err, log_f=log_f, mbh=f * vp, mbh_err=f * vp_err)


def storm_luminosity() -> dict:
    """NGC 5548's AGN luminosity during AGN STORM: mean 5100 A flux minus the
    starlight in the 5" x 15" aperture (Pei et al., from Bentz et al.'s HST image)."""
    s = LC.storm_spectroscopy("MDM", truncate=True)["F5100"]
    mean = float(np.mean(s.y))
    sem = float(np.std(s.y, ddof=1) / np.sqrt(len(s.y)))
    f_agn = mean - C.HOST_FLUX_STORM
    f_err = float(np.hypot(sem, C.HOST_FLUX_STORM_ERR))
    out = dict(F_total=mean, F_total_err=sem, F_agn=f_agn, F_agn_err=f_err, host_fraction=C.HOST_FLUX_STORM / mean)
    out["logL"] = float(P.log_lambda_l5100(f_agn))
    out["logL_err"] = float(0.4343 * f_err / f_agn)
    out["logL_total"] = float(P.log_lambda_l5100(mean))
    return out


# ---------------------------------------------------------------------------
# 2. The BLR holiday
# ---------------------------------------------------------------------------


def _segments_goad():
    seg = {k: (v[0] - 50000.0, v[1] - 50000.0) for k, v in LIT.GOAD16_SEGMENTS.items()}
    return seg


def _classify(t, segs):
    lab = np.full(len(t), "", dtype=object)
    for name, (lo, hi) in segs.items():
        lab[(t >= lo - 0.01) & (t <= hi + 0.01)] = name
    return lab


def _wls(x, y, e):
    """Weighted least squares y = A + eta x; returns (A, eta), (sA, seta)."""
    w = 1.0 / e**2
    X = np.column_stack([np.ones_like(x), x])
    cov = np.linalg.inv(X.T @ (X * w[:, None]))
    beta = cov @ (X.T @ (w * y))
    return beta, np.sqrt(np.diag(cov))


def holiday(n_mc: int = 500) -> dict:
    """Goad et al. 2016's anomaly analysis, recomputed for C IV, Si IV, He II, Ly-alpha
    and (as Pei et al. 2017 did) H-beta, all against the HST 1157.5 A continuum.

    For each line: (1) its lag over the first 75 days of the HST campaign;
    (2) shift the line back by that lag and remove the narrow-line flux;
    (3) fit log F_line = A + eta log F_cont on the epochs where line and
    continuum behave normally; (4) predict the line flux from the continuum and
    measure the fraction lost, F_lost = (F_pred - F_obs) / F_pred.
    """
    cont = LC.storm_hst_continuum(1157.5)
    uv = LC.storm_uv_lines()
    mdm = LC.storm_spectroscopy("MDM", truncate=False)
    lines = {k: uv[k] for k in ("Lya", "SiIV", "CIV", "HeII")}
    # H-beta: straight-line-continuum fluxes (in 1e-13 to match the UV units),
    # all 143 MDM epochs, as in Pei et al.'s figure 12.
    hb = mdm["Hb"]
    lines["Hb"] = LC.LightCurve(hb.t, hb.y / 100.0, hb.e / 100.0, "Hβ", "1e-13 erg/s/cm²")
    t_ref_end = C.HST_START_THJD + C.HOLIDAY_REFERENCE_DAYS
    lags = iccf.lag_grid(-10.0, 20.0, 0.25)
    out = {}
    for k, (name, ln) in enumerate(lines.items()):
        # (1) lag over the first 75 days
        c75 = cont.window(-np.inf, t_ref_end)
        l75 = ln.window(C.HST_START_THJD - 30.0, t_ref_end + 0.0) if name == "Hb" else ln.window(-np.inf, t_ref_end)
        res = iccf.fr_rss(c75.t, c75.y, c75.e, l75.t, l75.y, l75.e, lags, n_mc=n_mc, seed=300 + k)
        tau = res.cent[0]
        # (2) shift and remove narrow lines
        narrow = C.NARROW_FLUX[name]
        t_shift = ln.t - tau
        inside = (t_shift >= cont.t[0]) & (t_shift <= cont.t[-1])
        fl = ln.y[inside] - narrow
        el = ln.e[inside]
        ts = t_shift[inside]
        fc = np.interp(ts, cont.t, cont.y)
        ec = np.interp(ts, cont.t, cont.e)
        x, y = np.log10(fc), np.log10(fl)
        ey = 0.4343 * el / fl
        ex = 0.4343 * ec / fc
        # (3) segments: Goad's for the UV lines; Pei's (blue = before 6743) for H-beta
        if name == "Hb":
            b = C.HB_SEGMENTS_THJD
            segs = {"blue": (-np.inf, b[0]), "cyan": (b[0], b[1]), "red": (b[1], b[2]),
                    "magenta": (b[2], b[3]), "green": (b[3], np.inf)}
            fit_on = ("blue",)
            usable = ts <= C.STORM_TRUNCATE_THJD - tau  # last 10 MDM epochs badly calibrated
        else:
            segs = _segments_goad()
            fit_on = ("blue", "green")
            usable = np.ones(len(ts), bool)
        lab = _classify(ts, segs)
        sel = np.isin(lab, fit_on) & usable
        (a0, eta), (sa, seta) = _wls(x[sel], y[sel], np.hypot(ey[sel], ex[sel]))
        pred = 10 ** (a0 + eta * x)
        flost = (pred - fl) / pred
        red = (lab == "red") & usable
        ew = fl / fc  # equivalent width against the 1157 A continuum (units of 100 A)
        out[name] = dict(
            tau=res.cent, tau_r_max=res.observed.r_max, eta=(float(eta), float(seta)), A=float(a0),
            t_shift=ts, f_line=fl, f_cont=fc, f_pred=pred, f_lost=flost, label=lab, usable=usable,
            flost_red=float(np.mean(flost[red])) if red.any() else np.nan,
            flost_red_err=float(np.std(flost[red], ddof=1) / np.sqrt(red.sum())) if red.sum() > 1 else np.nan,
            n_red=int(red.sum()), ew=ew, fit_on=fit_on,
        )
        # The responsivity of the anomalous epochs alone (Pei et al.: 0.02 +/- 0.03 for H-beta)
        if red.sum() > 3:
            (_, eta_red), (_, seta_red) = _wls(x[red], y[red], np.hypot(ey[red], ex[red]))
            out[name]["eta_red"] = (float(eta_red), float(seta_red))
    out["_cont"] = cont
    return out


# ---------------------------------------------------------------------------
# 3. AGN Watch: 20 years of NGC 5548
# ---------------------------------------------------------------------------


def watch_campaigns() -> list[dict]:
    """The published H-beta campaigns of NGC 5548 inside the 1988-2008 light curves."""
    d = A.measurements()
    n = d[(d.object == "NGC5548") & (d.line == "Hb") & np.isfinite(d.tau_cent)].sort_values("jd_start")
    w = LC.watch()
    t0, t1 = w["Hb"].t[0], w["Hb"].t[-1]
    out = []
    for _, r in n.iterrows():
        if r.jd_start >= t0 - 1 and r.jd_end <= t1 + 1:
            out.append(dict(jd_start=r.jd_start, jd_end=r.jd_end, pub_tau=(r.tau_cent, r.tau_cent_lo, r.tau_cent_hi),
                            pub_sigma=(r.sigma_line, r.sigma_line_err), pub_logL=(r.logL5100, r.logL5100_err),
                            ref=r.ref_rm))
    return out


def watch_seasons(n_mc: int = C.N_MC) -> dict:
    """ICCF lag of each campaign from the 20-year light curves.

    As in the published yearly analyses (Peterson et al. 2002; 2004), the light
    curves are not detrended. The lag search runs from -T/4 to +T/3 for a
    campaign of length T (capped at -30..+80 d): at longer lags only a handful
    of points overlap and r can approach 1 by chance.
    """
    w = LC.watch()
    rows = []
    for k, camp in enumerate(watch_campaigns()):
        lo, hi = camp["jd_start"], camp["jd_end"]
        span = hi - lo
        lag_lo = max(C.WATCH_LAG_RANGE[0], -span / 4.0)
        lag_hi = min(C.WATCH_LAG_RANGE[1], span / 3.0)
        lags = iccf.lag_grid(np.floor(lag_lo), np.ceil(lag_hi), C.WATCH_LAG_STEP)
        c = w["F5100"].window(lo, hi)
        h = w["Hb"].window(lo, hi)
        res = iccf.fr_rss(c.t, c.y, c.e, h.t, h.y, h.e, lags, n_mc=n_mc, seed=500 + k)
        res_dt = _lag(c, h, lags, max(200, n_mc // 5), seed=900 + k, detrend=True)
        f_agn = float(np.mean(c.y))
        f_err = float(np.std(c.y, ddof=1) / np.sqrt(len(c.y)))
        logl = float(P.log_lambda_l5100(f_agn))
        rows.append(dict(**camp, n_cont=len(c), n_line=len(h), lag_range=(float(lags[0]), float(lags[-1])),
                         tau=_rest(res)["cent"], r_max=res.observed.r_max, success=res.success,
                         tau_detrended=_rest(res_dt)["cent"], f_agn=f_agn, f_agn_err=f_err,
                         logL=logl, logL_err=float(0.4343 * f_err / f_agn),
                         year=float(1858.87 + (lo + hi) / 2 / 365.25),
                         hb_mean=float(np.mean(h.y)), raw=res))
    return {"seasons": rows, "lcs": w}
