"""The whole analysis, from cached data to a results dict for figures and README.

    import analysis; res = analysis.run()

Sections (each a function returning plain dicts/arrays):
 1. light-curve showcase and the Hertzsprung progression
 2. Leavitt's 1912 stars against OGLE's SMC
 3. Leavitt laws in V, I and W_I for the LMC and SMC
 4. why the Wesenheit index works (colour-magnitude diagram)
 5. SMC vs LMC relative distance and the metallicity term
 6. the tilted LMC disk and the deep SMC (geometry)
 7. Gaia zero point from Milky Way parallaxes -> absolute LMC distance
 8. double-mode Cepheids (Petersen diagram)
"""

from __future__ import annotations

import numpy as np

import config as C
import lightcurves as LC
import magellanic
import ogle


# ---------------------------------------------------------------------------
# 1. Light curves
# ---------------------------------------------------------------------------

#: Showcase periods (days) and modes; the star nearest each is picked by rule:
#: isolated in Gaia (no other source within 2"), no OGLE remark, most V epochs.
SHOWCASE = ((1.7, "1O"), (3.0, "F"), (6.0, "F"), (10.0, "F"), (20.0, "F"), (40.0, "F"))


def _curve(gal, sid, band, period):
    lc = ogle.read_lightcurve(gal, sid, band)
    if lc is None or len(lc[0]) < 15:
        return None
    phase, mag, grid, model, mean = LC.model_for_plot(*lc, period, max_order=8 if band == "I" else 5)
    return {"phase": phase, "mag": mag, "grid": grid, "model": model, "mean": mean}


def lc_examples(tab) -> list:
    out = []
    for p_target, mode in SHOWCASE:
        # the longest periods are rare: widen the window and allow one faint Gaia neighbour there
        width = 0.08 if p_target < 30 else 0.35
        s = (tab["mode"] == mode) & (np.abs(np.log(tab["P"] / p_target)) < width)
        s &= (tab["gaia_n_nb"] == 1) | ((p_target >= 30) & (tab["gaia_nb_dG"] > 3))
        s &= (tab["remark"] == "") & np.isfinite(tab["lc_V_n"]) & (tab["lc_I_n"] >= 300)
        if not s.any():
            continue
        k = np.where(s)[0][np.argmax(tab["lc_V_n"][s])]
        sid, P = tab["id"][k], tab["P"][k]
        out.append({"id": sid, "P": P, "mode": mode, "I": _curve("lmc", sid, "I", P), "V": _curve("lmc", sid, "V", P)})
    return out


#: Period bins (days) of the Hertzsprung ladder: median light-curve shape per bin.
LADDER_PERIODS = (5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.5, 14.0, 16.0, 18.5, 22.0)


def _bump_phase(grid, y, ph_max=0.0):
    """Phase of the strongest secondary brightness maximum or shoulder, else NaN.

    ``y`` is brightness (1 = maximum). A bump is a local maximum of y away from
    maximum light; failing that, a 'shoulder': a local minimum of |dy/dphase| on the
    declining or rising branch.
    """
    from scipy.signal import find_peaks
    yy = np.concatenate([y, y])  # wrap
    gg = np.concatenate([grid, grid + 1])
    peaks, prop = find_peaks(yy, prominence=0.01)
    cand = [(prop["prominences"][j], gg[k] % 1.0) for j, k in enumerate(peaks)
            if 0.12 < (gg[k] % 1.0) < 0.9]
    if cand:
        return max(cand)[1], "bump"
    dy = np.abs(np.gradient(yy, gg))
    troughs, prop = find_peaks(-dy, prominence=0.3)
    cand = [(prop["prominences"][j], gg[k] % 1.0) for j, k in enumerate(troughs)
            if 0.12 < (gg[k] % 1.0) < 0.9 and yy[k] < 0.95]
    if cand:
        return max(cand)[1], "shoulder"
    return np.nan, ""


def hertzsprung(tab, half_width=0.05) -> dict:
    """The Hertzsprung progression: median I-band light-curve shape at fixed periods.

    Every LMC fundamental-mode Cepheid within +-5% of each period (and with at
    least 200 epochs, no OGLE remark) gets its Fourier model, phased to maximum
    light and normalised to brightness 0 (minimum) .. 1 (maximum); the median per
    bin is drawn. The progression's centre, where the bump coincides with
    maximum light, is where the second Fourier amplitude ratio R21 bottoms out.
    """
    grid = np.linspace(0, 1, 400, endpoint=False)
    ladder = []
    for p0 in LADDER_PERIODS:
        s = (tab["mode"] == "F") & (np.abs(np.log(tab["P"] / p0)) < half_width) & (tab["lc_I_n"] >= 200)
        s &= tab["remark"] == ""
        curves = []
        for k in np.where(s)[0]:
            lc = ogle.read_lightcurve("lmc", tab["id"][k], "I")
            _, _, g, model, _ = LC.model_for_plot(*lc, tab["P"][k], max_order=8)
            curves.append(np.interp(grid, g, (model.max() - model) / (model.max() - model.min())))
        if len(curves) < 3:
            continue
        med = np.median(curves, axis=0)
        ph, kind = _bump_phase(grid, med)
        ladder.append({"P": p0, "n": len(curves), "curve": med, "bump_phase": ph, "bump_kind": kind})
    sF = (tab["mode"] == "F") & np.isfinite(tab["lc_I_R21"]) & (tab["P"] > 6) & (tab["P"] < 16)
    lp, r21 = np.log10(tab["P"][sF]), tab["lc_I_R21"][sF]
    centres = np.linspace(np.log10(6.5), np.log10(15), 60)
    med = np.array([np.median(r21[np.abs(lp - c) < 0.03]) for c in centres])
    phi = {}
    for mode in ("F", "1O"):
        m = (tab["mode"] == mode) & np.isfinite(tab["lc_I_phi21"]) & (tab["lc_I_n"] >= 100)
        phi[mode] = {"logP": tab["logP"][m], "phi21": tab["lc_I_phi21"][m], "R21": tab["lc_I_R21"][m]}
    return {"grid": grid, "ladder": ladder, "p_center_r21": float(10 ** centres[np.argmin(med)]),
            "r21_running": {"logP": centres, "R21": med}, "phi21": phi}


# ---------------------------------------------------------------------------
# 3. Leavitt laws
# ---------------------------------------------------------------------------

def leavitt_laws(tabs) -> dict:
    """Clipped fits per galaxy, mode and band, plus the checks against OGLE's own numbers."""
    import plfit
    out = {}
    for gal, t in tabs.items():
        out[gal] = magellanic.pl_table(t)
    lmc = tabs["lmc"]
    sF = lmc["mode"] == "F"
    # the same fit on OGLE's catalogue magnitudes (to show our photometry changes nothing)
    out["lmc_ogle_phot"] = plfit.fit_line(lmc["logP"], lmc["W_I_ogle"], keep=sF)
    # Jacyszyn-Dobrzeniecka et al. quote W_I = a log P + b: our intercept at log P = 0
    out["lmc_F_WI_logp0"] = plfit.fit_line(lmc["logP"], lmc["W_I"], keep=sF, x0=0.0)
    out["lmc_1O_WI_logp0"] = plfit.fit_line(lmc["logP"], lmc["W_I"], keep=lmc["mode"] == "1O", x0=0.0)
    # is the law straight? break at 10 days, LMC fundamental mode
    out["break_10d"] = {band: plfit.fit_broken(lmc["logP"], lmc[band], 1.0, keep=sF) for band in ("V", "I", "W_I")}
    # the SMC's short-period break at 2.5 days (log P = 0.4)
    smc = tabs["smc"]
    out["smc_break_2p5d"] = plfit.fit_broken(smc["logP"], smc["W_I"], 0.4, keep=smc["mode"] == "F")
    return out


# ---------------------------------------------------------------------------
# 4. Why the Wesenheit index works
# ---------------------------------------------------------------------------

def wesenheit_physics(lmc) -> dict:
    """Period-luminosity-colour fit: at fixed period, how much brighter is a bluer Cepheid?

    I = a + b (log P - 1) + c (V - I). If c equals the dust ratio R = A_I/E(V-I)
    = 1.55, then the instability strip's width and dust push stars along the
    same direction, and W_I = I - 1.55 (V - I) removes both.
    """
    import plfit
    sF = lmc["mode"] == "F"
    plc = plfit.fit_plane_pl(lmc["logP"], lmc["VI"], lmc["I"], keep=sF)
    # colour spread at fixed period (the strip's width in V - I)
    pc = plfit.fit_line(lmc["logP"], lmc["VI"], keep=sF)
    return {"plc": {k: v for k, v in plc.items() if k != "mask"}, "period_colour": {k: pc[k] for k in ("a", "b", "sigma", "n")}}


# ---------------------------------------------------------------------------
# 5-6. Geometry, and the SMC-LMC offset
# ---------------------------------------------------------------------------

#: Centres adopted for the disk fits: LMC (80.05, -69.30) as Pietrzynski et al.
#: 2019 and Jacyszyn-Dobrzeniecka et al. 2016; SMC (12.54, -73.11), the Cepheid
#: centroid of Ripepi et al. 2017 adopted by Graczyk et al. 2020.
LMC_CENTRE = (80.05, -69.30)
SMC_CENTRE = (12.54, -73.11)
#: Fundamental-mode SMC Cepheids shorter than ~2-2.5 d fall off the law (Bauer
#: et al. 1999; Sharpee et al. 2002); like OGLE (Udalski et al. 1999; Soszynski
#: et al. 2010), the SMC fits use log P > 0.4 for F.
SMC_LOGP_MIN_F = 0.4


def geometry(tabs, n_boot=200) -> dict:
    """Residuals from the W_I law as distances; the LMC's tilted plane; the SMC's depth."""
    import geometry as GEO
    import plfit
    out = {}
    for gal, centre in (("lmc", LMC_CENTRE), ("smc", SMC_CENTRE)):
        t = tabs[gal]
        dmu = np.full(len(t["id"]), np.nan)
        fits = {}
        for mode in ("F", "1O"):
            keep = t["mode"] == mode
            if gal == "smc" and mode == "F":
                keep &= t["logP"] > SMC_LOGP_MIN_F
            f = plfit.fit_line(t["logP"], t["W_I"], keep=keep)
            fits[mode] = f
            dmu[f["mask"]] = f["resid"][f["mask"]]
        ok = np.isfinite(dmu)
        disk = GEO.fit_disk(t["ra"][ok], t["dec"][ok], dmu[ok], *centre)
        boot = GEO.bootstrap_disk(t["ra"][ok], t["dec"][ok], dmu[ok], *centre, n_boot=n_boot, seed=C.SEED)
        _, _, _, rho, pa = GEO.sky_offsets(t["ra"], t["dec"], *centre)
        corr = 5 * np.log10(GEO.disk_distance_ratio(rho, pa, disk["incl"], disk["theta"]))
        out[gal] = {"dmu": dmu, "disk": {k: v for k, v in disk.items() if k not in ("mask", "resid")},
                    "boot": boot, "geo_corr": corr, "rho_deg": np.degrees(rho), "pa_deg": np.degrees(pa),
                    "pl": {m: {k: f[k] for k in ("a", "b", "sigma", "n")} for m, f in fits.items()}}
        # scatter of the law before and after removing the plane
        for mode in ("F", "1O"):
            keep = fits[mode]["mask"]
            r_after = fits[mode]["resid"] - corr
            out[gal][f"sigma_{mode}_before"] = float(np.std(fits[mode]["resid"][keep]))
            out[gal][f"sigma_{mode}_after"] = float(np.std(r_after[keep] - np.mean(r_after[keep])))
        # near/far: mean residual in the NE vs SW quadrants (PA 0-90 vs 180-270)
        pa = np.degrees(pa)
        ne = ok & (pa < 90) & (np.degrees(rho) > 0.8)
        sw = ok & (pa >= 180) & (pa < 270) & (np.degrees(rho) > 0.8)
        out[gal]["ne_minus_sw"] = float(np.nanmean(dmu[ne]) - np.nanmean(dmu[sw]))
    return out


def smc_lmc(tabs, geo) -> dict:
    """Relative distance SMC - LMC in W_I, and the metallicity term it implies.

    Both laws share the LMC's slope; each galaxy's Cepheids are first moved to
    its centre with its fitted plane, as the eclipsing-binary distances refer
    to the centres. A metal-poor Cepheid that is intrinsically fainter makes
    the SMC look farther: gamma = (dmu_Cepheids - dmu_binaries) / ([Fe/H]_SMC - [Fe/H]_LMC).
    """
    import literature as LIT
    import plfit
    out = {}
    for mode in ("F", "1O"):
        lmc, smc = tabs["lmc"], tabs["smc"]
        kl = lmc["mode"] == mode
        ks = smc["mode"] == mode
        if mode == "F":
            ks &= smc["logP"] > SMC_LOGP_MIN_F
            kl &= lmc["logP"] > SMC_LOGP_MIN_F  # same period range in both
        res = {}
        for label, cl, cs in (("raw", 0.0, 0.0), ("centred", geo["lmc"]["geo_corr"], geo["smc"]["geo_corr"])):
            fl = plfit.fit_line(lmc["logP"], lmc["W_I"] - cl, keep=kl)
            fs = plfit.fit_line(smc["logP"], smc["W_I"] - cs, keep=ks, slope=fl["b"])
            d = fs["a"] - fl["a"]
            err = float(np.hypot(fs["ea"], fl["ea"]))
            dz = LIT.FEH_SMC["feh"] - LIT.FEH_LMC["feh"]
            gamma = (d - LIT.SMC_DEB["dmu_smc_lmc"]) / dz
            g_err = float(np.hypot(err, LIT.SMC_DEB["dmu_err"]) / abs(dz))
            res[label] = {"dmu": float(d), "err": err, "sigma_smc": fs["sigma"], "sigma_lmc": fl["sigma"],
                          "n_smc": fs["n"], "n_lmc": fl["n"], "slope": fl["b"], "gamma": float(gamma),
                          "gamma_err": g_err}
        # central SMC only (within 1 deg of the centre, where most eclipsing binaries are)
        cen = ks & (geo["smc"]["rho_deg"] < 1.0)
        fl = plfit.fit_line(lmc["logP"], lmc["W_I"] - geo["lmc"]["geo_corr"], keep=kl)
        fs = plfit.fit_line(smc["logP"], smc["W_I"], keep=cen, slope=fl["b"])
        d_c = float(fs["a"] - fl["a"])
        e_c = float(np.hypot(fs["ea"], fl["ea"]))
        dz = LIT.FEH_SMC["feh"] - LIT.FEH_LMC["feh"]
        res["smc_central_1deg"] = {"dmu": d_c, "err": e_c, "n_smc": fs["n"], "sigma_smc": fs["sigma"],
                                   "gamma": (d_c - LIT.SMC_DEB["dmu_smc_lmc"]) / dz,
                                   "gamma_err": float(np.hypot(e_c, LIT.SMC_DEB["dmu_err"]) / abs(dz))}
        out[mode] = res
    # the same test in Gaia's Wesenheit (both Clouds are barely reddened, so the
    # dust coefficient hardly matters here), for R = 1.90 and R = 1.74
    for R in (C.R_WG, 1.74):
        wl = lmc["gaia_G_int"] - R * (lmc["gaia_BP_int"] - lmc["gaia_RP_int"]) - geo["lmc"]["geo_corr"]
        ws = smc["gaia_G_int"] - R * (smc["gaia_BP_int"] - smc["gaia_RP_int"])
        kl = (lmc["mode"] == "F") & (lmc["logP"] > SMC_LOGP_MIN_F) & (lmc["gaia_cstar_sig"] < 3)
        ks = (smc["mode"] == "F") & (smc["logP"] > SMC_LOGP_MIN_F) & (smc["gaia_cstar_sig"] < 3) & \
            (geo["smc"]["rho_deg"] < 1.0)
        fl = plfit.fit_line(lmc["logP"], wl, keep=kl)
        fs = plfit.fit_line(smc["logP"], ws, keep=ks, slope=fl["b"])
        d = fs["a"] - fl["a"]
        err = float(np.hypot(fs["ea"], fl["ea"]))
        dz = LIT.FEH_SMC["feh"] - LIT.FEH_LMC["feh"]
        out[f"gaia_WG_R{R:.2f}"] = {"dmu": float(d), "err": err, "n_smc": fs["n"], "n_lmc": fl["n"],
                                    "gamma": float((d - LIT.SMC_DEB["dmu_smc_lmc"]) / dz),
                                    "gamma_err": float(np.hypot(err, LIT.SMC_DEB["dmu_err"]) / abs(dz))}
    # the SMC's depth: extra scatter over the LMC's, as a distance spread
    s_smc = out["F"]["centred"]["sigma_smc"]
    s_lmc = out["F"]["centred"]["sigma_lmc"]
    extra = np.sqrt(max(s_smc ** 2 - s_lmc ** 2, 0.0))
    d_smc = 10 ** (0.2 * LIT.SMC_DEB["mu"] + 1) / 1e3
    out["depth_sigma_mag"] = float(extra)
    out["depth_sigma_kpc"] = float(d_smc * np.log(10) / 5 * extra)
    return out


# ---------------------------------------------------------------------------
# 7. The absolute distance: Gaia parallaxes of Milky Way Cepheids
# ---------------------------------------------------------------------------

def mw_metallicities(mw) -> tuple:
    """High-resolution [Fe/H] per Milky Way calibrator: C-MetaLL IX, else Trentin+24 flag 0."""
    import literature_data as LD
    met = LD.metallicities()
    feh = np.full(len(mw["name"]), np.nan)
    src = np.array([""] * len(feh), dtype=object)
    for j, sid in enumerate(mw["source_id"]):
        sid = int(sid)
        if sid in met["cmetall9"]:
            feh[j], src[j] = met["cmetall9"][sid][0], "cmetall9"
        elif sid in met["trentin24"] and met["trentin24"][sid][2] == "0":
            feh[j], src[j] = met["trentin24"][sid][0], "trentin24"
    return feh, src.astype(str)


def mw_vi_photometry(mw) -> tuple:
    """Literature mean V and I of the calibrators (C-MetaLL IX, else Trentin+24)."""
    import literature_data as LD
    V = np.full(len(mw["name"]), np.nan)
    I = np.full(len(mw["name"]), np.nan)
    for name in ("trentin2024_table1.tsv", "cmetall9_tablea1.tsv"):  # C-MetaLL IX last, so it wins
        t = LD.read_vizier_tsv(name)
        row = {int(s): k for k, s in enumerate(t["GaiaDR3"]) if s.isdigit()}
        for j, sid in enumerate(mw["source_id"]):
            k = row.get(int(sid))
            if k is None:
                continue
            try:
                v, i = float(t["Vmag"][k]), float(t["Imag"][k])
            except ValueError:
                continue
            if np.isfinite(v) and np.isfinite(i) and v > 0 and i > 0:
                V[j], I[j] = v, i
    return V, I


BASE_FREE = ("mu", "gamma", "delta", "R", "sigR", "k", "d1O")


def gaia_calibration(tabs, geo, mcmc=True) -> dict:
    import gaia_cal as G
    import literature as LIT
    lmc = tabs["lmc"]
    corr = geo["lmc"]["geo_corr"]
    clean = np.isfinite(lmc["gaia_W_G"]) & (lmc["gaia_cstar_sig"] < 3)

    def shapes_for(geo_corr=corr, cut=clean, bands=("gaia_G_int", "gaia_BP_int", "gaia_RP_int")):
        out = {}
        for m in ("F", "1O"):
            keep = (lmc["mode"] == m) & cut
            out[m] = G.LMCShape(lmc["logP"], lmc[bands[0]], lmc[bands[1]], lmc[bands[2]], geo_corr=geo_corr, keep=keep)
        return out

    shapes = shapes_for()
    mw = G.milky_way()
    feh, feh_src = mw_metallicities(mw)
    calib = (G.calibrator_mask(mw, "F") | G.calibrator_mask(mw, "1O")) & np.isfinite(feh)
    sig_lmc = float(np.sqrt((shapes["F"].sigma_w ** 2 * shapes["F"].n + shapes["1O"].sigma_w ** 2 * shapes["1O"].n)
                            / (shapes["F"].n + shapes["1O"].n)))
    base_fixed = {"lnsig": np.log(sig_lmc)}
    feh_lmc = LIT.FEH_LMC["feh"]

    def run_fit(sel=calib, free=BASE_FREE, fixed=None, shp=None, data=mw, sample=False):
        jf = G.JointFit(data, sel, feh, shp or shapes, feh_lmc, free=free, fixed={**base_fixed, **(fixed or {})})
        jf.fit()
        if sample:
            jf.sample()
        return jf

    base = run_fit(sample=mcmc)
    res = {"sigma_int_lmc": sig_lmc, "n_calib": int(calib.sum()), "feh_src": {s: int(np.sum(calib & (feh_src == s)))
                                                                           for s in ("cmetall9", "trentin24")},
           "lmc_shapes": {m: {"a_190": sh.a(C.R_WG), "b_190": sh.b(C.R_WG), "sigma": sh.sigma_w, "n": sh.n,
                              "e_a": sh.e_a, "a_g": sh.a_g, "b_g": sh.b_g, "a_c": sh.a_c, "b_c": sh.b_c}
                          for m, sh in shapes.items()},
           "baseline": base.summary()}
    # the LMC zero point's own error adds to mu in quadrature
    p = base.params(base.theta)
    res["baseline"]["e_mu_lmc_zero_point"] = float(shapes["F"].e_a)
    # quantities for the figures: each calibrator's absolute W from its parallax, and the model
    r, var, pred, pp = base.terms(base.theta)
    s = base.sel
    w = mw["int_average_g"][s] - pp["R"] * (mw["int_average_bp"][s] - mw["int_average_rp"][s])
    plx_c = mw["parallax"][s] - mw["zp_l21"][s] - pp["delta"]
    res["points"] = {"name": mw["name"][s], "mode": mw["mode"][s], "logP": mw["logP"][s], "W": w,
                     "plx": plx_c, "e_plx": mw["parallax_error"][s], "k": pp["k"], "feh": feh[s], "E": base.E,
                     "M_model": base.absolute_w(pp), "kept": base.mask, "G": mw["phot_g_mean_mag"][s],
                     "l": mw["l"][s], "b": mw["b"][s], "ra": mw["ra"][s], "dec": mw["dec"][s], "P": mw["P"][s]}
    if mcmc:
        res["chain"] = {n: base.chain[:, j] for j, n in enumerate(base.free)}

    # ---- variants: each changes one assumption ----
    V = {}
    nolnsig = tuple(f for f in BASE_FREE)
    V["intrinsic scatter free"] = run_fit(free=BASE_FREE + ("lnsig",))
    V["dust coefficient R fixed at 1.90"] = run_fit(free=tuple(f for f in BASE_FREE if f != "R"))
    V["parallax offset fixed: Riess+21 (14 uas)"] = run_fit(free=tuple(f for f in BASE_FREE if f != "delta"),
                                                             fixed={"delta": 0.014})
    V["no residual parallax offset (L21 only)"] = run_fit(free=tuple(f for f in BASE_FREE if f != "delta"),
                                                          fixed={"delta": 0.0})
    V["no metallicity term"] = run_fit(free=tuple(f for f in BASE_FREE if f != "gamma"), fixed={"gamma": 0.0})
    V["gamma fixed: Breuval+22 (-0.384)"] = run_fit(free=tuple(f for f in BASE_FREE if f != "gamma"),
                                                    fixed={"gamma": -0.384})
    V["gamma fixed: C-MetaLL IX (-0.527)"] = run_fit(free=tuple(f for f in BASE_FREE if f != "gamma"),
                                                     fixed={"gamma": -0.527})
    V["catalogue parallax errors (k = 1)"] = run_fit(free=tuple(f for f in BASE_FREE if f != "k"), fixed={"k": 1.0})
    V["fundamental mode only"] = run_fit(sel=calib & (mw["mode"] == "F"), free=tuple(f for f in BASE_FREE if f != "d1O"))
    V["first overtones only"] = run_fit(sel=calib & (mw["mode"] == "1O"), free=tuple(f for f in BASE_FREE if f != "d1O"),
                                        fixed={"d1O": 0.0})
    V["C-MetaLL IX metallicities only"] = run_fit(sel=calib & (feh_src == "cmetall9"))
    V["Gaia also classifies as DCEP"] = run_fit(sel=calib & (mw["gaia_type"] == "DCEP"))
    V["G > 8 (no very bright stars)"] = run_fit(sel=calib & (mw["phot_g_mean_mag"] > 8))
    V["BP-RP < 2 (less reddened)"] = run_fit(sel=calib & (mw["bp_rp"] < 2.0))
    V["LMC: no geometric correction"] = run_fit(shp=shapes_for(geo_corr=None))
    V["LMC: no crowding cut"] = run_fit(shp=shapes_for(cut=np.isfinite(lmc["gaia_W_G"])))
    # an independent photometric route: OGLE's W_I = I - 1.55 (V - I) with literature V, I
    Vm, Im = mw_vi_photometry(mw)
    mw_vi = dict(mw)
    mw_vi["int_average_g"], mw_vi["int_average_bp"], mw_vi["int_average_rp"] = Im, Vm, Im
    shapes_vi = shapes_for(cut=np.isfinite(lmc["W_I"]), bands=("I", "V", "I"))
    sel_vi = calib & np.isfinite(Vm) & np.isfinite(Im)
    V["OGLE W_I route (literature V, I; R fitted)"] = run_fit(sel=sel_vi, shp=shapes_vi, data=mw_vi,
                                                              fixed={"R": C.R_WI, "lnsig": np.log(shapes_vi["F"].sigma_w)})
    V["OGLE W_I route, R fixed at 1.55"] = run_fit(sel=sel_vi, shp=shapes_vi, data=mw_vi,
                                                   free=tuple(f for f in BASE_FREE if f != "R"),
                                                   fixed={"R": C.R_WI, "lnsig": np.log(shapes_vi["F"].sigma_w)})
    res["variants"] = {}
    for idx, (name, jf) in enumerate(V.items()):
        # every variant gets its own injection test, so each is bias-corrected for its own set-up
        inj = injection_test(jf, jf.shapes, feh_lmc, n_sims=40, seed=C.SEED + 100 + idx, tag=f"_v{idx:02d}")
        prm = {k: float(v) for k, v in jf.params(jf.theta).items()}
        res["variants"][name] = {"n": int(jf.mask.sum()), **prm, "mu_raw": prm["mu"], "mu_bias": inj["mu_bias"],
                                 "mu_bias_err": inj["mu_bias_err"], "mu": prm["mu"] - inj["mu_bias"]}
    res["injection"] = injection_test(base, shapes, feh_lmc)
    # the same at a different truth (literature-like gamma, R = 1.90, smaller offset): is the bias stable?
    res["injection_alt"] = injection_test(base, shapes, feh_lmc, n_sims=40, seed=C.SEED + 1, tag="_alt",
                                          truth_override={"gamma": -0.45, "R": 1.90, "delta": 0.015})
    res["n_vi"] = int(sel_vi.sum())
    return res


def injection_test(base, shapes, feh_lmc, n_sims=100, seed=C.SEED, truth_override=None, tag="") -> dict:
    """Is the calibration unbiased, and is its error bar honest?

    Each simulation keeps the real calibrators (periods, modes, metallicities,
    colours and reddenings, parallax errors and zero points) and takes a model
    as the truth: each star's true distance is its photometric distance under
    that model, and new Wesenheit magnitudes (scatter in magnitudes, as in
    nature) and parallaxes are drawn with the model's intrinsic scatter,
    dust-law scatter, error inflation and parallax offset. The fit is then
    rerun and compared with the truth.

    Finding: the fit returns mu_LMC ~0.03 mag too large. The cause is
    regression dilution (an errors-in-variables effect): in parallax space the
    model is plx = delta + s * 10^(-0.2 W), with the LMC distance in the scale s
    and the photometric scatter in the input W. Noise in an input flattens the
    fitted scale (s too small = mu too large) and the free offset delta absorbs
    the rest (+2-4 uas). With delta held fixed the bias vanishes; it does not
    depend on the parallax errors, and it shrinks with the photometric scatter.
    """
    import gaia_cal as G
    cache = C.CACHE_DIR / f"injection_test{tag}.npz"
    truth = base.params(base.theta)
    if truth_override:
        truth.update(truth_override)
    names = ("mu", "gamma", "delta", "R", "k")
    key = np.array([truth[k] for k in ("mu", "gamma", "delta", "R", "sigR", "k", "d1O")] + [n_sims, seed])
    if cache.exists():
        z = np.load(cache)
        if z["key"].shape == key.shape and np.allclose(z["key"], key):
            return _injection_summary(z["out"], truth, names)
    rng = np.random.default_rng(seed)
    M = base.absolute_w(truth)
    w_obs = base.g - truth["R"] * (base.bp - base.rp)
    mu_star = w_obs - M
    plx_true = 10 ** (-0.2 * (mu_star - 10.0))
    sig_mag = np.sqrt(np.exp(2 * truth["lnsig"]) + (truth["sigR"] * base.E) ** 2)
    fixed = {k: v for k, v in truth.items() if k not in base.free}
    out = []
    for _ in range(n_sims):
        g_sim = base.g + rng.normal(0, sig_mag)  # colours kept; the law's scatter goes into G
        plx_sim = plx_true + base.zp + truth["delta"] + rng.normal(0, truth["k"] * base.e_plx)
        mw = {"int_average_g": g_sim, "int_average_bp": base.bp, "int_average_rp": base.rp, "logP": base.logp,
              "parallax": plx_sim, "parallax_error": base.e_plx, "zp_l21": base.zp, "mode": base.mode}
        jf = G.JointFit(mw, np.ones(len(g_sim), bool), base.feh, shapes, feh_lmc, free=base.free, fixed=fixed)
        jf.fit()
        out.append([jf.params(jf.theta)[k] for k in names])
    out = np.array(out)
    np.savez(cache, key=key, out=out)
    return _injection_summary(out, truth, names)


def _injection_summary(out, truth, names) -> dict:
    res = {"n_sims": int(len(out)), "truth": {n: float(truth[n]) for n in names}}
    for j, n in enumerate(names):
        res[f"{n}_bias"] = float(np.mean(out[:, j]) - truth[n])
        res[f"{n}_bias_err"] = float(np.std(out[:, j], ddof=1) / np.sqrt(len(out)))
        res[f"{n}_spread"] = float(np.std(out[:, j], ddof=1))
    return res


# ---------------------------------------------------------------------------
# 8. Double-mode Cepheids (Petersen diagram)
# ---------------------------------------------------------------------------

def petersen(tabs) -> dict:
    out = {}
    for gal, t in tabs.items():
        f1 = np.isfinite(t["P_F"]) & np.isfinite(t["P_1O"])
        o12 = np.isfinite(t["P_1O"]) & np.isfinite(t["P_2O"])
        out[gal] = {"F1O": {"P": t["P_F"][f1], "ratio": t["P_1O"][f1] / t["P_F"][f1]},
                    "1O2O": {"P": t["P_1O"][o12], "ratio": t["P_2O"][o12] / t["P_1O"][o12]}}
    mw = ogle.galactic_list()
    f1 = np.isfinite(mw["P_F"]) & np.isfinite(mw["P_1O"])
    o12 = np.isfinite(mw["P_1O"]) & np.isfinite(mw["P_2O"])
    out["mw"] = {"F1O": {"P": mw["P_F"][f1], "ratio": mw["P_1O"][f1] / mw["P_F"][f1]},
                 "1O2O": {"P": mw["P_1O"][o12], "ratio": mw["P_2O"][o12] / mw["P_1O"][o12]}}
    for key in ("lmc", "smc", "mw"):
        for kind in ("F1O", "1O2O"):
            r = out[key][kind]["ratio"]
            out[key][kind]["median"] = float(np.median(r)) if len(r) else np.nan
            out[key][kind]["n"] = int(len(r))
    return out


# ---------------------------------------------------------------------------

def run(mcmc: bool = True) -> dict:
    tabs = magellanic.master()
    res = {"tabs": tabs}
    print("  light curves ...", flush=True)
    res["lc_check"] = {gal: magellanic.lightcurve_check(t) for gal, t in tabs.items()}
    res["lc_examples"] = lc_examples(tabs["lmc"])
    res["hertzsprung"] = hertzsprung(tabs["lmc"])
    print("  Leavitt laws, geometry ...", flush=True)
    res["laws"] = leavitt_laws(tabs)
    res["wesenheit"] = wesenheit_physics(tabs["lmc"])
    res["geometry"] = geometry(tabs)
    res["smc_lmc"] = smc_lmc(tabs, res["geometry"])
    print("  Gaia zero point ...", flush=True)
    res["gaia"] = gaia_calibration(tabs, res["geometry"], mcmc=mcmc)
    res["petersen"] = petersen(tabs)
    res["leavitt"] = leavitt1912(tabs["smc"])
    return res


def leavitt1912(smc) -> dict:
    import leavitt
    c = leavitt.compare(smc)
    c["aliases"] = leavitt.aliases(c["table"])
    t = c["table"]
    flagged = {a["H"] for a in c["aliases"]}
    c["n_agree_5digits"] = int(sum(1 for h in t["H"] if h not in flagged))
    c["n_ogle_smc"] = int(len(smc["id"]))
    return c


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

#: Variants that swap one modelling *assumption* for another defensible one;
#: the RMS of their shifts in mu_LMC is our systematic error. (Subsample splits,
#: such as F-only, are consistency checks: their shifts are mostly noise.)
SYS_VARIANTS = ("intrinsic scatter free", "dust coefficient R fixed at 1.90",
                "parallax offset fixed: Riess+21 (14 uas)", "catalogue parallax errors (k = 1)",
                "LMC: no geometric correction", "LMC: no crowding cut",
                "OGLE W_I route (literature V, I; R fitted)")


def _f(d, keys=("a", "b", "ea", "eb", "sigma", "n")):
    return {k: (float(d[k]) if k != "n" else int(d[k])) for k in keys}


def summarize(res) -> dict:
    import literature as LIT
    tabs = res["tabs"]
    s = {"counts": {gal: {m: int(np.sum(t["mode"] == m)) for m in np.unique(t["mode"])} for gal, t in tabs.items()}}
    s["counts"]["total_magellanic"] = int(sum(len(t["id"]) for t in tabs.values()))
    s["lc_check"] = res["lc_check"]
    s["hertzsprung_centre_r21_days"] = res["hertzsprung"]["p_center_r21"]
    laws = res["laws"]
    s["laws"] = {f"{gal}_{m}_{b}": _f(f) for gal in ("lmc", "smc") for (m, b), f in laws[gal].items()}
    # against Soszynski et al. 2015 Table 2 (same OGLE-IV data; intercept at log P = 0)
    s["s15_comparison"] = []
    for (gal, mode, band), (slope, zp, sig) in LIT.S15_TABLE2.items():
        f = laws[gal][(mode, band)]
        s["s15_comparison"].append({"gal": gal, "mode": mode, "band": band, "slope": float(f["b"]),
                                    "e_slope": float(f["eb"]), "zp": float(f["a"] - f["b"]), "sigma": float(f["sigma"]),
                                    "n": int(f["n"]), "s15_slope": slope, "s15_zp": zp, "s15_sigma": sig})
    s["lmc_F_WI_logp0"] = _f(laws["lmc_F_WI_logp0"])
    s["lmc_1O_WI_logp0"] = _f(laws["lmc_1O_WI_logp0"])
    s["lmc_F_WI_ogle_photometry"] = _f(laws["lmc_ogle_phot"])
    s["break_10d"] = {b: {"short": float(v["short"]["b"]), "e_short": float(v["short"]["eb"]),
                          "long": float(v["long"]["b"]), "e_long": float(v["long"]["eb"]), "n_long": v["n_long"],
                          "F": v["F"], "p": v["p"]} for b, v in laws["break_10d"].items()}
    v = laws["smc_break_2p5d"]
    s["smc_break_2p5d"] = {"short": float(v["short"]["b"]), "long": float(v["long"]["b"]), "p": v["p"]}
    s["wesenheit"] = res["wesenheit"]
    g = res["geometry"]
    s["geometry"] = {gal: {**{k: float(x) for k, x in g[gal]["disk"].items() if k != "n"}, "n": g[gal]["disk"]["n"],
                           "incl_boot": g[gal]["boot"]["incl_std"], "theta_boot": g[gal]["boot"]["theta_std"],
                           "ne_minus_sw": g[gal]["ne_minus_sw"],
                           **{k: float(g[gal][k]) for k in g[gal] if k.startswith("sigma_")}} for gal in ("lmc", "smc")}
    sl = res["smc_lmc"]
    s["smc_lmc"] = {k: v for k, v in sl.items()}
    # ---- Gaia: the absolute distance ----
    ga = res["gaia"]
    b = ga["baseline"]
    inj, inj_alt = ga["injection"], ga["injection_alt"]
    mu_raw = b["mu"]["median"]
    mu = mu_raw - inj["mu_bias"]  # corrected for the estimator's measured bias (regression dilution)
    stat = float(np.sqrt(b["mu"]["std"] ** 2 + ga["lmc_shapes"]["F"]["e_a"] ** 2 + inj["mu_bias_err"] ** 2))
    base_corr = b["map"]["mu"] - inj["mu_bias"]
    shifts = np.array([ga["variants"][n]["mu"] - base_corr for n in SYS_VARIANTS])
    bias_sys = abs(inj_alt["mu_bias"] - inj["mu_bias"])  # how much the bias depends on the assumed truth
    sys_ = float(np.sqrt(np.mean(shifts ** 2) + bias_sys ** 2))
    tot = float(np.hypot(stat, sys_))
    d = 10 ** (0.2 * mu + 1) / 1e3
    deb = LIT.LMC_DEB
    deb_tot = float(np.hypot(deb["stat"], deb["sys"]))
    s["lmc_distance"] = {
        "mu": mu, "mu_raw": mu_raw, "bias": inj["mu_bias"], "bias_err": inj["mu_bias_err"],
        "bias_alt_truth": inj_alt["mu_bias"], "bias_sys": bias_sys,
        "gamma_corrected": b["gamma"]["median"] - inj["gamma_bias"],
        "delta_corrected_uas": 1e3 * (b["delta"]["median"] - inj["delta_bias"]),
        "R_corrected": b["R"]["median"] - inj["R_bias"],
        "mu": mu, "stat": stat, "sys": sys_, "total": tot, "d_kpc": d, "d_err_kpc": d * np.log(10) / 5 * tot,
        "deb_mu": deb["mu"], "diff": mu - deb["mu"], "n_sigma": (mu - deb["mu"]) / float(np.hypot(tot, deb_tot)),
        "sys_shifts": dict(zip(SYS_VARIANTS, map(float, shifts))),
    }
    s["gaia"] = {"baseline": b, "variants": ga["variants"], "n_calib": ga["n_calib"], "feh_src": ga["feh_src"],
                 "sigma_int_lmc": ga["sigma_int_lmc"], "lmc_shapes": ga["lmc_shapes"], "n_vi": ga["n_vi"],
                 "injection": inj, "injection_alt": inj_alt}
    # ---- what it means for H0 (SH0ES, LMC as the only anchor) ----
    h = LIT.H0["riess22_lmc_only"]
    h0 = h["H0"] * 10 ** (-0.2 * (mu - h["mu_lmc"]))
    need = h["mu_lmc"] + 5 * np.log10(h["H0"] / LIT.H0["planck18"]["H0"])
    s["h0"] = {"H0_lmc_anchor_ours": float(h0), "H0_err_from_mu": float(h0 * np.log(10) / 5 * tot),
               "mu_needed_for_planck": float(need), "d_needed_kpc": float(10 ** (0.2 * need + 1) / 1e3),
               "sigma_from_needed": float((need - mu) / tot)}
    s["petersen"] = {k: {kind: {"median": v[kind]["median"], "n": v[kind]["n"]} for kind in ("F1O", "1O2O")}
                     for k, v in res["petersen"].items()}
    lv = res["leavitt"]
    t = lv["table"]
    s["leavitt"] = {"n_matched": lv["n_matched"], "modes": sorted(set(t["ogle_mode"])),
                    "slope_max": _f(lv["fit_max"]), "slope_min": _f(lv["fit_min"]), "slope_V_today": _f(lv["fit_v"]),
                    "slope_Vmax_today": _f(lv["fit_vmax"]), "her_stated_slope": -1 / 0.48,
                    "aliases": lv["aliases"], "n_agree_within_5_last_digits": lv["n_agree_5digits"],
                    "ratio_sample": lv["n_ogle_smc"] / 25.0 if "n_ogle_smc" in lv else None,
                    "stars": [{"HV": int(t["H"][j]), "ogle": str(t["ogle_id"][j]), "P1912": float(t["P"][j]),
                               "P_ogle": float(t["P_ogle"][j]), "dP_over_P": float(t["dP_over_P"][j])}
                              for j in range(len(t["H"]))]}
    return s


def print_summary(res) -> None:
    s = summarize(res)
    L = s["laws"]
    print("\n" + "=" * 78)
    print("RESULTS: Cepheids and the Leavitt law")
    print("=" * 78)
    c = s["counts"]
    print(f"OGLE-IV Magellanic Cepheids: {c['total_magellanic']} (LMC {sum(c['lmc'].values())}, SMC {sum(c['smc'].values())})")
    for gal in ("lmc", "smc"):
        for m in ("F", "1O"):
            lc = s["lc_check"][gal][m]
            print(f"  {gal.upper()} {m}: blind period = OGLE's for {100 * lc['blind_period_recovered']:.2f}%; "
                  f"<I> ours - OGLE {lc['dI_median']:+.4f} (MAD {lc['dI_mad']:.4f})")
    print("\nLMC Leavitt laws (our photometry), m = a + b (log P - 1):")
    for m in ("F", "1O"):
        for band in ("V", "I", "W_I"):
            f = L[f"lmc_{m}_{band}"]
            print(f"  {m:2s} {band:3s}: a = {f['a']:.3f}  b = {f['b']:.3f} ± {f['eb']:.3f}  sigma = {f['sigma']:.3f}  N = {f['n']}")
    f = s["lmc_F_WI_logp0"]
    print(f"  F W_I as a log P + b: {f['b']:.3f} log P + {f['a']:.3f} (Soszynski+15 Table 2: -3.314 log P + 15.888)")
    worst = max(s["s15_comparison"], key=lambda r: abs(r["slope"] - r["s15_slope"]))
    print(f"  all 12 Soszynski+15 fits: largest slope difference {worst['slope'] - worst['s15_slope']:+.3f} "
          f"({worst['gal']} {worst['mode']} {worst['band']}); largest zero-point difference "
          f"{max(abs(r['zp'] - r['s15_zp']) for r in s['s15_comparison']):.3f} mag")
    w = s["wesenheit"]["plc"]
    print(f"  Period-luminosity-colour: I = ... + {w['c']:.3f} (V - I)  (dust: 1.55); scatter {w['sigma']:.3f}")
    for gal in ("lmc", "smc"):
        g = s["geometry"][gal]
        print(f"{gal.upper()} plane: i = {g['incl']:.1f} ± {max(g['e_incl'], g['incl_boot']):.1f} deg, "
              f"line of nodes PA = {g['theta']:.1f} ± {max(g['e_theta'], g['theta_boot']):.1f}, near side PA "
              f"{g['near_side_pa']:.0f}; NE - SW = {g['ne_minus_sw']:+.3f} mag")
    sl = s["smc_lmc"]
    print(f"SMC - LMC (W_I, F, central SMC): {sl['F']['smc_central_1deg']['dmu']:.3f} ± "
          f"{sl['F']['smc_central_1deg']['err']:.3f} (binaries 0.500 ± 0.017); depth sigma {sl['depth_sigma_kpc']:.1f} kpc")
    print(f"SMC - LMC (Gaia W_G): {sl['gaia_WG_R1.90']['dmu']:.3f} -> gamma = {sl['gaia_WG_R1.90']['gamma']:+.2f} ± "
          f"{sl['gaia_WG_R1.90']['gamma_err']:.2f}")
    d = s["lmc_distance"]
    b = s["gaia"]["baseline"]
    print(f"\nLMC distance from {b['n']} Milky Way Cepheid parallaxes (Gaia DR3):")
    print(f"  mu = {d['mu']:.3f} ± {d['stat']:.3f} (stat) ± {d['sys']:.3f} (sys) = {d['d_kpc']:.2f} ± "
          f"{d['d_err_kpc']:.2f} kpc   [eclipsing binaries: 18.477 ± 0.026; {d['n_sigma']:+.1f} sigma]")
    inj = s["gaia"]["injection"]
    print(f"  raw fit {d['mu_raw']:.3f}; injection test ({inj['n_sims']} simulations): the fit returns mu "
          f"{d['bias']:+.3f} ± {d['bias_err']:.3f} too far (regression dilution; {d['bias_alt_truth']:+.3f} at a "
          f"literature-like truth); scatter of fits {inj['mu_spread']:.3f} vs MCMC error {b['mu']['std']:.3f}")
    print(f"  bias-corrected: gamma = {d['gamma_corrected']:+.3f} ± {b['gamma']['std']:.3f}, delta = "
          f"{d['delta_corrected_uas']:.1f} ± {1e3 * b['delta']['std']:.1f} uas, R = {d['R_corrected']:.3f} ± {b['R']['std']:.3f}")
    for n in ("k", "sigR", "d1O"):
        if n in b:
            print(f"  {n:6s} = {b[n]['median']:.3f} ± {b[n]['std']:.3f}")
    for name, v in s["gaia"]["variants"].items():
        print(f"    {name:46s} mu = {v['mu']:.3f} (raw {v['mu_raw']:.3f}, bias {v['mu_bias']:+.3f})  "
              f"(gamma {v['gamma']:+.2f}, delta {1e3 * v['delta']:+.0f} uas, R {v['R']:.2f})")
    h = s["h0"]
    print(f"\nH0 with the LMC as the only anchor (SH0ES 2022 rescaled): {h['H0_lmc_anchor_ours']:.1f} ± "
          f"{h['H0_err_from_mu']:.1f}; Planck's 67.4 would need mu_LMC = {h['mu_needed_for_planck']:.3f} "
          f"({h['d_needed_kpc']:.1f} kpc), {h['sigma_from_needed']:.1f} sigma from ours")
    p = s["petersen"]
    print(f"Double-mode period ratio P1O/PF: MW {p['mw']['F1O']['median']:.4f}, LMC {p['lmc']['F1O']['median']:.4f}, "
          f"SMC {p['smc']['F1O']['median']:.4f}")
    print()
