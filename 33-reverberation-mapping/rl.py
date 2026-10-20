"""The radius-luminosity (R-L) relation, from three samples.

Bigger, brighter AGN have bigger broad-line regions: R ~ L^alpha with alpha
near 1/2, the value expected if every AGN's broad lines form where the
ionizing photon flux (and so the gas's ionization state) is about the same.
This relation is what turns a single spectrum into a black hole mass for
distant quasars (project I3), so its slope and scatter set their precision.

* **Local AGN** (the AGN Black Hole Mass Database): H-beta lags with
  host-galaxy starlight removed from L5100 using HST images. Fitted with
  LINMIX_ERR after checking that the same code reproduces Bentz et al. (2013)
  on their own table.
* **SDSS-RM** quasars (Shen et al. 2024): H-beta, Mg II and C IV lags of
  distant quasars from one survey, with their published lags, and with lags
  we measured ourselves (``sdssrm.py``).
* **NGC 5548 alone**: one black hole over 25 years. Does its own BLR follow
  the population's relation as its luminosity changes?
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from astropy.io import fits

import agnmass as A
import config as C
import literature as LIT
import regression as R

SDSSRM_SUMMARY = C.DATA_DIR / "sdssrm" / "summary.fits"

#: Luminosity normalization (log erg/s) and continuum column for each line.
SDSSRM_LINE = {
    "HB": dict(name="Hb", lcol="L5100_MEDIAN", norm=44.0, host=True, label="Hβ", lam=5100),
    "MG2": dict(name="MgII", lcol="L3000_MEDIAN", norm=45.0, host=False, label="Mg II", lam=3000),
    "C4": dict(name="CIV", lcol="L1350_MEDIAN", norm=45.0, host=False, label="C IV", lam=1350),
    "HA": dict(name="Ha", lcol="L5100_MEDIAN", norm=44.0, host=True, label="Hα", lam=5100),
}


def _fit_dict(f: R.LineFit) -> dict:
    return dict(a=f.a, b=f.b, sig=f.sig, n=f.n, x0=f.x0)


# ---------------------------------------------------------------------------
# Local AGN
# ---------------------------------------------------------------------------


def bentz13_validation() -> dict:
    """Our LINMIX run and our own likelihood fit on Bentz et al.'s 71 measurements."""
    b = LIT.bentz13_table()
    x = b.logL.values - 44.0
    y = np.log10(b.tau.values)
    sy = R.log_err(b.tau, b.tau_lo, b.tau_hi)
    sx = b.logL_err.values
    lm = R.linmix_fit(x, sx, y, sy)
    ml = R.fit_line(x, sx, y, sy)
    return dict(linmix=_fit_dict(lm), likelihood=_fit_dict(ml), n=len(b), n_obj=b.object.nunique())


def local_sample() -> pd.DataFrame:
    """Every H-beta campaign in the database with a lag and a host-corrected L5100."""
    d = A.hbeta_rl().copy()
    d["x"] = d.logL5100 - 44.0
    d["y"] = np.log10(d.tau_cent)
    d["sy"] = R.log_err(d.tau_cent, d.tau_cent_lo, d.tau_cent_hi)
    d["sx"] = d.logL5100_err
    return d


def local_rl() -> dict:
    d = local_sample()
    fit = R.linmix_fit(d.x, d.sx, d.y, d.sy)
    # One point per AGN (mean of its campaigns), the other way to count (Bentz et al. 2013).
    g = d.groupby("object").agg(x=("x", "mean"), y=("y", "mean"), sx=("sx", "mean"), sy=("sy", "mean"), n=("x", "size"))
    g["sy"] = g.sy / np.sqrt(g.n)
    fit1 = R.linmix_fit(g.x, g.sx, g.y, g.sy, seed=5)
    # Fixed slope 0.533 (for comparison with SDSS-RM's offset)
    fixed = R.fit_line(d.x, d.sx, d.y, d.sy, fix_b=0.533)
    return dict(all=_fit_dict(fit), per_object=_fit_dict(fit1), fixed_0533=_fit_dict(fixed),
                n=len(d), n_obj=d.object.nunique(), sample=d, per_object_sample=g)


# ---------------------------------------------------------------------------
# SDSS-RM
# ---------------------------------------------------------------------------


def sdssrm_table() -> pd.DataFrame:
    """The SDSS-RM summary as a DataFrame with the quantities we use."""
    d = fits.getdata(SDSSRM_SUMMARY, 1)
    out = dict(RMID=d["RMID"].astype(int), z=d["ZSYS"].astype(float), f_host=d["F_H_5100_YR1"].astype(float),
               lag_exp_obs=d["LAG_OBS_EXP"].astype(float))
    for key in ("L5100_MEDIAN", "L3000_MEDIAN", "L1350_MEDIAN"):
        out[key] = d[key].astype(float)
    out["search_lo"], out["search_hi"] = d["LAG_SEARCH"][:, 0].astype(float), d["LAG_SEARCH"][:, 1].astype(float)
    for ln in SDSSRM_LINE:
        out[f"{ln}_done"] = d[f"{ln}_LAG_DONE"].astype(int)
        out[f"{ln}_lag"] = d[f"{ln}_LAG_DET"].astype(float)
        out[f"{ln}_lag_lo"], out[f"{ln}_lag_hi"] = (d[f"{ln}_LAG_DET_ERR"][:, i].astype(float) for i in (0, 1))
        out[f"{ln}_iccf"] = d[f"{ln}_ICCF_LAG"].astype(float)
        out[f"{ln}_iccf_lo"], out[f"{ln}_iccf_hi"] = (d[f"{ln}_ICCF_LAG_ERR"][:, i].astype(float) for i in (0, 1))
        out[f"{ln}_iccf_rmax"] = d[f"{ln}_ICCF_LAG_RMAX"].astype(float)
        out[f"{ln}_iccf_rej"] = d[f"{ln}_ICCF_FRAC_REJ"].astype(float)
        out[f"{ln}_grade"] = d[f"{ln}_LAG_GRADE"].astype(float)
        out[f"{ln}_sig_rms"] = d[f"{ln}_W_RMS"][:, 1].astype(float)
        out[f"{ln}_sig_rms_err"] = d[f"{ln}_W_RMS_ERR"][:, 1].astype(float)
        out[f"{ln}_fwhm_mean"] = d[f"{ln}_W_MEAN"][:, 0].astype(float)
        out[f"{ln}_mrm"] = d[f"{ln}_MRM"].astype(float)
    df = pd.DataFrame(out)
    df["L5100_agn"] = df.L5100_MEDIAN * (1.0 - df.f_host.clip(lower=0.0))
    return df


def sdssrm_rl_sample(line: str, lag_col: str = "lag") -> pd.DataFrame:
    """Detected lags (Shen et al.'s fiducial PyROA lags by default) with luminosities."""
    t = sdssrm_table()
    cfg = SDSSRM_LINE[line]
    det = t[f"{line}_lag"] > 0
    s = t[det].copy()
    L = s.L5100_agn if cfg["host"] else s[cfg["lcol"]]
    s["logL"] = np.log10(L.where(L > 0))
    s["tau"] = s[f"{line}_{lag_col}"] / (1 + s.z)
    s["tau_lo"] = s[f"{line}_{lag_col}_lo"] / (1 + s.z)
    s["tau_hi"] = s[f"{line}_{lag_col}_hi"] / (1 + s.z)
    s = s[np.isfinite(s.logL) & (s.tau > 0)]
    s["x"] = s.logL - cfg["norm"]
    s["y"] = np.log10(s.tau)
    s["sy"] = R.log_err(s.tau, s.tau_lo, s.tau_hi)
    s["sx"] = 1e-3  # Shen et al. quote no luminosity errors; median-luminosity errors are tiny
    return s.reset_index(drop=True)


def sdssrm_rl() -> dict:
    out = {}
    for line in ("HB", "MG2", "C4"):
        s = sdssrm_rl_sample(line)
        fit = R.linmix_fit(s.x, s.sx, s.y, s.sy)
        fixed = R.fit_line(s.x, s.sx, s.y, s.sy, fix_b=LIT.SHEN24_RL[SDSSRM_LINE[line]["name"]]["b_fixed"])
        out[line] = dict(fit=_fit_dict(fit), fixed=_fit_dict(fixed), n=len(s), sample=s)
    # H-beta without removing host starlight: how much the starlight matters.
    s = sdssrm_rl_sample("HB")
    x_total = np.log10(s.L5100_MEDIAN) - 44.0
    out["HB_no_host"] = _fit_dict(R.linmix_fit(x_total, s.sx, s.y, s.sy, seed=13))
    return out


def combined_hbeta(local: pd.DataFrame, sdss: pd.DataFrame) -> dict:
    """Local AGN + SDSS-RM H-beta lags in one fit (Shen et al. warn that the
    survey's selection makes this fit illustrative only)."""
    x = np.concatenate([local.x, sdss.x])
    sx = np.concatenate([local.sx, sdss.sx])
    y = np.concatenate([local.y, sdss.y])
    sy = np.concatenate([local.sy, sdss.sy])
    return _fit_dict(R.linmix_fit(x, sx, y, sy, seed=11))


# ---------------------------------------------------------------------------
# NGC 5548 alone
# ---------------------------------------------------------------------------


def ngc5548_breathing(seasons: list[dict], storm: dict | None = None) -> dict:
    """NGC 5548's own R-L slope from our seasonal lags.

    Bentz et al. 2007 fitted AGN Watch years 1-13 plus 2005; we do the same
    with our lags and luminosities, then without Year 12 (2000), then add the
    later campaigns. The luminosity errors used are the statistical errors of
    each season's mean flux: the host-galaxy subtraction is common to all
    seasons, so its error is treated as a systematic (``host_systematic``).
    """
    rows = []
    for s in seasons:
        tau, lo, hi = s["tau"]
        rows.append(dict(year=s["year"], x=s["logL"] - 44.0, sx=max(s["logL_err"], 0.005), tau=tau, lo=lo, hi=hi,
                         jd=s["jd_start"], ref=s["ref"]))
    d = pd.DataFrame(rows)
    d["y"] = np.log10(d.tau.clip(lower=0.1))
    d["sy"] = R.log_err(d.tau.clip(lower=0.1), d.lo, d.hi)
    watch = d.jd < 53500  # AGN Watch years 1-13 + 2005 (Bentz et al. 2007's sample)
    y12 = (d.jd > 51500) & (d.jd < 51800)
    good = d.tau > 0
    out = dict(sample=d)
    out["bentz07_sample"] = _fit_dict(R.linmix_fit(d.x[watch & good], d.sx[watch & good], d.y[watch & good],
                                                   d.sy[watch & good], seed=21))
    m = watch & ~y12 & good
    out["bentz07_no_y12"] = _fit_dict(R.linmix_fit(d.x[m], d.sx[m], d.y[m], d.sy[m], seed=22))
    out["all_watch"] = _fit_dict(R.linmix_fit(d.x[good], d.sx[good], d.y[good], d.sy[good], seed=23))
    # The same fit with the published per-season luminosity errors (~0.1 dex,
    # mostly the host-starlight correction, which is common to all seasons):
    # LINMIX then corrects for "noise" in L that does not scatter season to
    # season, and the slope steepens to Bentz et al. 2007's value.
    sx_pub = np.array([s["pub_logL"][1] for s in seasons])
    out["bentz07_sample_pub_errors"] = _fit_dict(R.linmix_fit(d.x[watch & good], sx_pub[watch & good],
                                                              d.y[watch & good], d.sy[watch & good], seed=24))
    # And the published numbers alone (database lags, luminosities and errors):
    # with the same code this reproduces Bentz et al. 2007's slope.
    pub = np.array([s["pub_tau"] for s in seasons])
    xl_pub = np.array([s["pub_logL"][0] for s in seasons]) - 44.0
    m = watch.to_numpy()
    out["published_values"] = _fit_dict(R.linmix_fit(xl_pub[m], sx_pub[m], np.log10(pub[m, 0]),
                                                     R.log_err(pub[m, 0], pub[m, 1], pub[m, 2]), seed=25))
    # Host-subtraction systematic: the starlight in the AGN Watch 5" x 7.5"
    # aperture is 3.75 +/- 0.375 e-15 erg/s/cm2/A (Bentz et al. 2013, Table 13).
    # Subtracting 0.375 more (or less) from every season steepens (flattens)
    # the slope, because faint seasons lose a larger share of their flux.
    f_agn = np.array([s["f_agn"] for s in seasons])
    const = d.x.values - np.log10(f_agn)  # log(4 pi D^2 lambda) - 44
    slopes = {}
    for dh in (-0.375, +0.375):
        xx = np.log10(np.clip(f_agn - dh, 0.05, None)) + const
        slopes[dh] = R.fit_line(xx[m], d.sx[m], d.y[m], d.sy[m]).b[0]
    out["host_systematic"] = dict(host_flux=3.75, host_err=0.375, slope_less_host=slopes[-0.375],
                                  slope_more_host=slopes[0.375])
    return out
