"""Black hole masses: the virial factor f, NGC 5548's mass, and the M-sigma relation.

A lag and a line width give the virial product VP = c tau sigma^2 / G. The
black hole mass is f x VP, where f (of order 1-10) absorbs the unknown shape
and orientation of the broad-line region. f is calibrated two ways:

* **Dynamical models** (Pancoast et al. 2014; Williams et al. 2018; Villafana
  et al. 2022; ...): fitting the full velocity-resolved echo of an individual
  AGN with a model of orbiting clouds gives its mass directly; mass / VP is
  that AGN's own f. Shen et al. (2024) compiled 30 such measurements.
* **The M-sigma relation**: assume reverberation-mapped AGN follow the same
  M-sigma relation as quiescent galaxies whose black holes were weighed by
  stellar or gas dynamics, and find the f that puts them there (Onken et al.
  2004; Woo et al. 2010, 2013, 2015).

We redo both: the first exactly as Shen et al. did, the second as Woo et al.
2015 did but with van den Bosch (2016)'s 160 quiescent galaxies instead of
Kormendy & Ho (2013)'s 84.
"""

from __future__ import annotations

import emcee
import numpy as np
import pandas as pd
from astropy.table import Table
from scipy import optimize

import agnmass as A
import config as C
import literature as LIT
import physics as P
import regression as R

VDB16 = C.DATA_DIR / "vdb16" / "BHcompilation.fits"
LOG_SIG0 = np.log10(200.0)


def quiescent_galaxies() -> pd.DataFrame:
    """van den Bosch (2016) black holes weighed by stars, gas, CO or masers in
    his fitted sample (SELECTED = 1), upper limits excluded."""
    t = Table.read(VDB16)
    typ = np.array([str(x).strip() for x in t["TYPE"]])
    m = (np.asarray(t["SELECTED"]) == 1) & np.isin(typ, ["star", "stars", "gas", "CO", "maser"]) \
        & (np.asarray(t["UPPERLIMIT"]) == 0)
    return pd.DataFrame(dict(name=[str(n).strip() for n in t["NAME"][m]], logM=np.asarray(t["MBH"][m], float),
                             logM_err=np.asarray(t["DMBH"][m], float), logsig=np.asarray(t["SIG"][m], float),
                             logsig_err=np.asarray(t["DSIG"][m], float), method=typ[m]))


def rm_agn_woo15() -> pd.DataFrame:
    """Woo et al. 2015 Table A1: log VP (sigma_line, rms) and log sigma_* of 29 AGN."""
    w = LIT.woo15_table().copy()
    w["logVP"] = np.log10(w.vp)
    w["logVP_err"] = 0.4343 * 0.5 * (w.vp_lo + w.vp_hi) / w.vp
    w["logsig"] = np.log10(w.sigma_star)
    w["logsig_err"] = 0.4343 * w.sigma_star_err / w.sigma_star
    return w


# ---------------------------------------------------------------------------
# f from dynamical models
# ---------------------------------------------------------------------------


def f_dynamical() -> dict:
    """Mean and intrinsic scatter of log f over Shen et al.'s 30 dynamical models."""
    t = LIT.shen24_f_table()
    mu, sig, ml = R.gaussian_mean_scatter(t.logf, t.logf_err)
    return dict(mean=mu, scatter=sig, ml=ml, n=len(t), table=t)


# ---------------------------------------------------------------------------
# f from the M-sigma relation (joint fit, Woo et al. 2013, 2015)
# ---------------------------------------------------------------------------


def _joint_lnlike(p, q, a):
    alpha, beta, logf, eps = p
    if eps <= 0 or not (0 < beta < 12) or not (-2 < logf < 3):
        return -np.inf
    sq = q.logM_err**2 + (beta * q.logsig_err) ** 2 + eps**2
    rq = q.logM - alpha - beta * (q.logsig - LOG_SIG0)
    sa = a.logVP_err**2 + (beta * a.logsig_err) ** 2 + eps**2
    ra = a.logVP + logf - alpha - beta * (a.logsig - LOG_SIG0)
    return -0.5 * (np.sum(rq**2 / sq + np.log(sq)) + np.sum(ra**2 / sa + np.log(sa)))


def f_msigma(quiescent: pd.DataFrame | None = None, agn: pd.DataFrame | None = None, nsteps: int = 4000) -> dict:
    """Joint fit of log M = alpha + beta log(sigma/200) to quiescent galaxies and
    to the AGN's log VP + log f, with one intrinsic scatter for both."""
    q = quiescent_galaxies() if quiescent is None else quiescent
    a = rm_agn_woo15() if agn is None else agn
    qa = {k: q[k].to_numpy() for k in ("logM", "logM_err", "logsig", "logsig_err")}
    aa = {k: a[k].to_numpy() for k in ("logVP", "logVP_err", "logsig", "logsig_err")}
    qd = pd.DataFrame(qa)
    ad = pd.DataFrame(aa)
    p0 = np.array([8.3, 5.0, 0.6, 0.4])
    ml = optimize.minimize(lambda p: -_joint_lnlike(p, qd, ad), p0, method="Nelder-Mead",
                           options=dict(maxiter=20000, xatol=1e-7, fatol=1e-9)).x
    rng = np.random.default_rng(17)
    nw = 32
    start = ml + 1e-3 * rng.standard_normal((nw, 4))
    s = emcee.EnsembleSampler(nw, 4, _joint_lnlike, args=(qd, ad))
    s.random_state = np.random.RandomState(17).get_state()
    s.run_mcmc(start, nsteps, progress=False)
    ch = s.get_chain(discard=nsteps // 3, thin=5, flat=True)

    def qn(v):
        lo, md, hi = np.percentile(v, [15.865525, 50, 84.134475])
        return float(md), float(md - lo), float(hi - md)

    # Quiescent galaxies alone, for the M-sigma line drawn in the figures
    qfit = R.fit_line(q.logsig - LOG_SIG0, q.logsig_err, q.logM, q.logM_err)
    return dict(alpha=qn(ch[:, 0]), beta=qn(ch[:, 1]), logf=qn(ch[:, 2]), eps=qn(ch[:, 3]), ml=tuple(ml),
                n_q=len(q), n_agn=len(a), quiescent_only=dict(a=qfit.a, b=qfit.b, sig=qfit.sig),
                quiescent=q, agn=a)


# ---------------------------------------------------------------------------
# NGC 5548
# ---------------------------------------------------------------------------


def ngc5548_virial(storm_lag_rest: tuple | None = None) -> dict:
    """Line width vs lag for every broad line and campaign of NGC 5548.

    If gravity dominates, every line in every season has the same tau sigma^2:
    log sigma = const - 0.5 log tau. Uses the database's centroid lags and rms
    widths; the 2014 H-beta lag can be replaced by ours.
    """
    d = A.measurements()
    n = d[(d.object == "NGC5548") & np.isfinite(d.tau_cent) & np.isfinite(d.sigma_line) & (d.tau_cent > 0)].copy()
    if storm_lag_rest is not None:
        m = (n.line == "Hb") & (n.jd_start > 56600)
        n.loc[m, ["tau_cent", "tau_cent_lo", "tau_cent_hi"]] = storm_lag_rest
    n["x"] = np.log10(n.tau_cent)
    n["sx"] = R.log_err(n.tau_cent, n.tau_cent_lo, n.tau_cent_hi)
    n["y"] = np.log10(n.sigma_line)
    n["sy"] = 0.4343 * n.sigma_line_err / n.sigma_line
    fit = R.linmix_fit(n.x, n.sx, n.y, n.sy, seed=31)
    fixed = R.fit_line(n.x, n.sx, n.y, n.sy, fix_b=-0.5)
    hb = n.line == "Hb"
    fit_hb = R.linmix_fit(n.x[hb], n.sx[hb], n.y[hb], n.sy[hb], seed=32)
    n["vp"] = P.virial_product(n.tau_cent, n.sigma_line)
    n["vp_err"] = P.virial_product_err(n.tau_cent, 0.5 * (n.tau_cent_lo + n.tau_cent_hi), n.sigma_line,
                                       n.sigma_line_err)
    w = 1.0 / (0.4343 * n.vp_err / n.vp) ** 2
    logvp_mean = float(np.sum(w * np.log10(n.vp)) / np.sum(w))
    return dict(sample=n, fit_all=dict(a=fit.a, b=fit.b, sig=fit.sig, n=fit.n),
                fit_hb=dict(a=fit_hb.a, b=fit_hb.b, sig=fit_hb.sig, n=fit_hb.n),
                fixed_half=dict(a=fixed.a, sig=fixed.sig), logvp_weighted=logvp_mean,
                logvp_scatter=float(np.std(np.log10(n.vp), ddof=1)))


def ngc5548_on_msigma(vp: float, vp_err: float, msig: dict, fdyn: dict) -> dict:
    """NGC 5548's mass for each virial factor, and where M-sigma would put it."""
    s, se = C.SIGMA_STAR_NGC5548, C.SIGMA_STAR_NGC5548_ERR
    out = dict(sigma_star=(s, se), vp=vp, vp_err=vp_err)
    for key, logf in (("woo15", (C.LOG_F_WOO15, 0.12)), ("shen24", LIT.SHEN24_LOGF["mean"]),
                      ("ours_dyn", (fdyn["mean"][0], 0.5 * (fdyn["mean"][1] + fdyn["mean"][2]))),
                      ("ours_msigma", (msig["logf"][0], 0.5 * (msig["logf"][1] + msig["logf"][2])))):
        m = 10 ** logf[0] * vp
        out[key] = dict(logf=logf, mbh=m, logM=float(np.log10(m)),
                        logM_err=float(np.hypot(0.4343 * vp_err / vp, logf[1])))
    a, b = msig["alpha"][0], msig["beta"][0]
    out["msigma_prediction"] = float(a + b * (np.log10(s) - LOG_SIG0))
    out["msigma_offset_dex"] = out["ours_msigma"]["logM"] - out["msigma_prediction"]
    out["msigma_offset_sigma"] = out["msigma_offset_dex"] / msig["eps"][0]
    return out


def local_agn_masses(log_f: float) -> pd.DataFrame:
    """One mass per local AGN: weighted mean log VP over its H-beta campaigns."""
    d = A.measurements()
    h = d[(d.line == "Hb") & np.isfinite(d.tau_cent) & np.isfinite(d.sigma_line) & (d.tau_cent > 0)].copy()
    h["logvp"] = np.log10(P.virial_product(h.tau_cent, h.sigma_line))
    h["logvp_err"] = 0.4343 * P.virial_product_err(h.tau_cent, 0.5 * (h.tau_cent_lo + h.tau_cent_hi), h.sigma_line,
                                                   h.sigma_line_err) / P.virial_product(h.tau_cent, h.sigma_line)
    h = h[np.isfinite(h.logvp_err) & (h.logvp_err > 0)]
    rows = []
    for name, g in h.groupby("object"):
        w = 1.0 / g.logvp_err**2
        lv = float(np.sum(w * g.logvp) / np.sum(w))
        rows.append(dict(object=name, n=len(g), logVP=lv, logVP_err=float(1 / np.sqrt(np.sum(w))), logM=lv + log_f,
                         z=float(g.z.iloc[0])))
    return pd.DataFrame(rows).sort_values("logM").reset_index(drop=True)
