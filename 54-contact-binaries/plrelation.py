"""Contact binaries as standard candles: period-luminosity relations with Gaia DR3.

A contact binary's two stars share one envelope, so their combined size and
temperature, hence luminosity, are tied to the orbital period (Kepler's law
fixes the separation; the stars fill it). Chen et al. (2018) calibrated the
relation on 183 nearby W UMa stars with TGAS parallaxes; Gaia DR3 gives us
thousands of ZTF contact binaries with parallaxes better than 5%.

Method notes:
* Parallaxes get the Lindegren et al. (2021) zero point (a function of
  magnitude, colour and ecliptic latitude; the ``gaiadr3-zeropoint`` code),
  and we test the sensitivity to it explicitly: a 0.02 mas error moves a
  1-kpc star by 0.04 mag.
* Fits are done in magnitude space only for precise parallaxes (sigma_M <
  0.11 mag), where the Lutz-Kelker/regression-dilution bias is small; we also
  fit in parallax space (predicted vs measured parallax), which has no such bias.
* Extinction: W1 is nearly free of it (A_W1 ~ 0.04 A_V); for G we use the
  extinction-free Wesenheit index W_G = G - 1.90 (BP - RP) and, as a check, a
  dust-layer model like Chen et al.'s.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import config as C

#: Lindegren et al. (2021) global parallax zero point for 5-parameter solutions
#: (QSO median), used as a constant here; its spread with colour and magnitude
#: is ~0.01-0.02 mas.
PARALLAX_ZP_MAS = -0.017

#: Wesenheit coefficient for Gaia: W_G = G - 1.90 (BP - RP) (Ripepi et al. 2019).
WESENHEIT_R = 1.90

#: Relative extinction (Wang & Chen 2019): A_W1/A_V = 0.039, A_G/A_V ~ 0.789.
A_W1_OVER_AV = 0.039
A_G_OVER_AV = 0.789


def build_sample(ztf: pd.DataFrame, xm: pd.DataFrame, allwise: pd.DataFrame | None = None) -> pd.DataFrame:
    """One Gaia counterpart per ZTF binary: the nearest within the radius,
    rejected if a second Gaia source within 1.5" is less than 2 mag fainter
    (blends would bias both the light curve and the luminosity)."""
    xm = xm.sort_values(["zid", "sep_arcsec"])
    g = xm.groupby("zid")
    first = g.head(1).set_index("zid")
    # brightness of the next neighbour
    second = g.nth(1).set_index("zid") if hasattr(g, "nth") else None
    blend = pd.Series(False, index=first.index)
    if second is not None and len(second):
        dm = second["g"] - first.loc[second.index, "g"]
        blend.loc[second.index] = (dm < 2.0).to_numpy()
    first["blended"] = blend
    s = ztf.set_index("zid").join(first, how="inner", rsuffix="_gaia")
    if allwise is not None:
        s = s.join(allwise.set_index("zid")[["w1", "w1_err", "w2", "w2_err", "allwise_sep", "w1_clean"]], how="left")
        s["w1_clean"] = s["w1_clean"].fillna(False).astype(bool)
    s["zp_l21"] = lindegren_zero_point(s)
    s["plx_corr"] = s["parallax"] - np.where(np.isfinite(s["zp_l21"]), s["zp_l21"], PARALLAX_ZP_MAS)
    s["plx_snr"] = s["plx_corr"] / s["parallax_error"]
    s["logp"] = np.log10(s["period"])
    s["w_g"] = s["g"] - WESENHEIT_R * s["bp_rp"]
    with np.errstate(invalid="ignore", divide="ignore"):
        s["dist_pc"] = 1000.0 / s["plx_corr"]
        s["dm"] = 5 * np.log10(s["dist_pc"]) - 5
        s["dm_err"] = 2.1715 * s["parallax_error"] / s["plx_corr"]
    return s.reset_index()


def lindegren_zero_point(s: pd.DataFrame) -> np.ndarray:
    """Lindegren et al. (2021) parallax zero point (mas), NaN outside its
    calibrated range (6 < G < 21, 1.1 < nu_eff < 1.9) or for 2-parameter solutions."""
    import warnings
    from astropy import units as u
    from astropy.coordinates import BarycentricMeanEcliptic, SkyCoord
    from zero_point import zpt
    zpt.load_tables()
    ecl = SkyCoord(s["ra_gaia"].to_numpy() * u.deg, s["dec_gaia"].to_numpy() * u.deg).transform_to(
        BarycentricMeanEcliptic())
    solved = s["solved"].fillna(0).astype(int).to_numpy()
    nu = s["nu_eff"].to_numpy(float)
    pc = s["pseudocolour"].to_numpy(float)
    g = s["g"].to_numpy(float)
    ok = np.isin(solved, (31, 95)) & np.isfinite(g) & (g > 6) & (g < 21)
    ok &= np.where(solved == 31, np.isfinite(nu) & (nu > 1.1) & (nu < 1.9), np.isfinite(pc) & (pc > 1.24) & (pc < 1.72))
    out = np.full(len(s), np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out[ok] = zpt.get_zpt(g[ok], np.where(np.isfinite(nu[ok]), nu[ok], 1.5), np.where(np.isfinite(pc[ok]), pc[ok], 1.5),
                              ecl.lat.deg[ok], solved[ok], _warnings=False)
    return out


def dust_layer_av(dist_pc, b_deg, rho0=0.54e-3, h=164.0):
    """Chen et al. (2018)'s dust-layer extinction: A_V = rho0 int_0^d exp(-r |sin b| / H) dr,
    rho0 in mag/pc (0.54 mag/kpc) and H = 164 pc (Nataf et al. 2013)."""
    sb = np.abs(np.sin(np.radians(b_deg)))
    sb = np.where(sb < 1e-4, 1e-4, sb)
    return rho0 * h / sb * (1 - np.exp(-dist_pc * sb / h))


def fit_pl(logp, mag, err, logp0=-0.4, clip=4.0, iters=5, weighted=True):
    """Straight-line PL fit M = a (log P - logp0) + b with an intrinsic scatter,
    iterated 4-sigma clipping; returns slope, intercept at log P = logp0, the
    intercept at log P = 0 (Chen et al.'s convention), scatter and errors.
    ``weighted=False`` gives ordinary least squares, which is what Chen et al.
    (2018) used (their Table 2 is reproduced exactly that way)."""
    keep = np.isfinite(logp) & np.isfinite(mag) & np.isfinite(err)
    if not weighted:
        err = np.zeros_like(np.asarray(err, float))
    for _ in range(iters):
        x, y, e = logp[keep] - logp0, mag[keep], err[keep]
        A = np.column_stack([x, np.ones_like(x)])
        sig_int = 0.1
        for _ in range(20):
            w = 1 / np.sqrt(e**2 + sig_int**2)
            coef, *_ = np.linalg.lstsq(A * w[:, None], y * w, rcond=None)
            r = y - A @ coef
            # update intrinsic scatter from the residual variance
            new = np.sqrt(max(np.mean(r**2) - np.mean(e**2), 1e-6))
            if abs(new - sig_int) < 1e-5:
                break
            sig_int = new
        r_all = mag - (coef[0] * (logp - logp0) + coef[1])
        mad = 1.4826 * np.median(np.abs(r_all[keep]))
        new_keep = keep & (np.abs(r_all) < clip * mad)
        if new_keep.sum() == keep.sum():
            break
        keep = new_keep
    w = 1 / np.sqrt(err[keep] ** 2 + sig_int**2)
    A = np.column_stack([logp[keep] - logp0, np.ones(keep.sum())])
    cov = np.linalg.inv((A * w[:, None]).T @ (A * w[:, None]))
    a, b = coef
    return {"slope": float(a), "slope_err": float(np.sqrt(cov[0, 0])), "zp": float(b), "zp_err": float(np.sqrt(cov[1, 1])),
            "b0": float(b - a * logp0), "b0_err": float(np.sqrt(cov[1, 1] + logp0**2 * cov[0, 0] - 2 * logp0 * cov[0, 1])),
            "scatter": float(np.std(r_all[keep])), "sigma_int": float(sig_int), "n": int(keep.sum()), "keep": keep,
            "logp0": logp0}


def parallax_space_test(logp, mag, plx, plx_err, slope, b0, a_lambda=0.0):
    """Predicted parallax (mas) from a PL relation vs the measured one.

    Returns the error-weighted mean ratio measured/predicted, which is free of
    the Lutz-Kelker bias of magnitude-space comparisons (Feast & Catchpole 1997).
    """
    m_abs = slope * logp + b0
    plx_pred = 10 ** (-0.2 * (mag - a_lambda - m_abs - 10.0)) * 1.0  # mas: d = 10^(0.2(m - M + 5)) pc
    ratio = plx / plx_pred
    w = (plx_pred / plx_err) ** 2
    mean = np.sum(w * ratio) / np.sum(w)
    err = 1 / np.sqrt(np.sum(w))
    return {"ratio": float(mean), "ratio_err": float(err), "dmag": float(-5 * np.log10(mean)),
            "plx_pred": plx_pred}


def chen2018_recalibration() -> dict:
    """Chen et al. (2018)'s 183 W UMa calibrators, first with their own TGAS
    numbers (to reproduce Table 2), then with Gaia DR3 parallaxes (Lindegren
    zero point) and AllWISE W1 for the same stars.

    Their Table 1 lists absolute magnitudes; the apparent W1 they used is
    M + DM_TGAS + A_W1 with A_W1 = 0.05 A_V (their Section 3.2).
    """
    import catalogs
    t1, _ = catalogs.read_chen2018()
    t1 = t1.reset_index(drop=True)
    t1["logp"] = np.log10(t1["Period"])
    out = {"n": int(len(t1))}
    # 1. their own numbers
    for col, key in (("W1Mag", "W1_mean"), ("Gmag", "G_mean"), ("W1MagMax", "W1_max"), ("VMagMax", "V_max"),
                     ("GMagMax", "G_max")):
        f = fit_pl(t1["logp"].to_numpy(), t1[col].to_numpy(), t1["e_DM"].to_numpy(), logp0=0.0, clip=99,
                   weighted=False)
        out[f"reproduce_{key}"] = {k: f[k] for k in ("slope", "slope_err", "b0", "b0_err", "scatter", "n")}
    # 2. Gaia DR3 parallaxes for the same stars
    gx = pd.read_csv(C.ZTF_DIR / "chen2018_table1_x_gaiadr3.csv").sort_values(["zid", "Gmag"])
    gx = gx.groupby("zid").head(1).set_index("zid")  # brightest Gaia source within 3"
    aw = pd.read_csv(C.ZTF_DIR / "chen2018_table1_x_allwise.csv", dtype={"ccf": str, "qph": str})
    aw = aw.sort_values(["zid", "angDist"]).groupby("zid").head(1).set_index("zid")
    s = t1.join(gx[["Source", "Plx", "e_Plx", "RUWE", "Gmag", "BPmag", "RPmag", "nueff", "pscol", "Solved", "RAdeg",
                    "DEdeg"]].rename(columns={"Gmag": "g_dr3"}), how="left")
    s = s.join(aw[["W1mag", "e_W1mag", "ccf", "qph"]], how="left")
    zp_in = pd.DataFrame({"ra_gaia": s["RAdeg"], "dec_gaia": s["DEdeg"], "solved": s["Solved"], "nu_eff": s["nueff"],
                          "pseudocolour": s["pscol"], "g": s["g_dr3"]})
    s["zp"] = lindegren_zero_point(zp_in)
    s["plx"] = s["Plx"] - np.where(np.isfinite(s["zp"]), s["zp"], PARALLAX_ZP_MAS)
    s["dm_dr3"] = 5 * np.log10(1000.0 / s["plx"]) - 5
    s["dm_dr3_err"] = 2.1715 * s["e_Plx"] / s["plx"]
    a_w1 = 0.05 * s["Av"]
    s["w1_app_chen"] = s["W1Mag"] + s["DM"] + a_w1
    out["dm_dr3_minus_tgas"] = {"median": float(np.nanmedian(s["dm_dr3"] - s["DM"])),
                                "scatter": float(1.4826 * np.nanmedian(np.abs(s["dm_dr3"] - s["DM"] -
                                                                              np.nanmedian(s["dm_dr3"] - s["DM"])))),
                                "median_err_tgas": float(np.median(s["e_DM"])),
                                "median_err_dr3": float(np.nanmedian(s["dm_dr3_err"]))}
    have = np.isfinite(s["dm_dr3"])
    good = have & (s["RUWE"] < 1.4)
    out["n_dr3"] = int(have.sum())
    out["n_dr3_ruwe"] = int(good.sum())
    for label, app, sel in (("W1_mean_dr3_all", s["w1_app_chen"], have),
                            ("W1_mean_dr3_ruwe", s["w1_app_chen"], good),
                            ("W1_mean_dr3_allwise_ruwe", s["W1mag"], good),
                            ("W1_mean_tgas_same_stars", s["w1_app_chen"], good)):
        dm = s["DM"] if "tgas" in label else s["dm_dr3"]
        M = app - a_w1 - dm
        k = sel & np.isfinite(M)
        f = fit_pl(s["logp"][k].to_numpy(), M[k].to_numpy(), np.zeros(int(k.sum())), logp0=0.0, clip=99,
                   weighted=False)
        out[label] = {kk: f[kk] for kk in ("slope", "slope_err", "b0", "b0_err", "scatter", "n")}
    s.to_parquet(C.DERIVED_DIR / "chen2018_dr3.parquet")
    return out
