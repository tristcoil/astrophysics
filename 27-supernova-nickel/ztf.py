"""Arnett's rule for hundreds of Type Ia supernovae from the ZTF Bright Transient Survey.

For each supernova, from its ZTF g and r alert light curves (ALeRCE):

* the g-band peak (time and magnitude), the r magnitude at that moment, and
  the decline over the next 15 rest-frame days, Delta m15(g);
* the time of first light, from a "fireball" fit f ~ (t - t_fl)^2 to the
  earliest detections (non-detections forbid an earlier start);
* the peak bolometric luminosity, by scaling SN 2011fe: its bolometric peak
  (measured in this project) times the ratio of the two supernovae's
  combined g + r flux at g-band peak. The g + r flux of a Type Ia at peak is a
  steady fraction of its total light, so the ratio carries over; how steady
  is the main systematic (Watch out for, in the README). K-corrections come
  from redshifting SN 2011fe's own spectrum (Pereira et al. 2013);
* the nickel mass from Arnett's rule, with the rise time from first light.

Distances from redshift with H0 = 73 (the Cepheid scale, matching the M101
distance used for SN 2011fe).
"""

from __future__ import annotations

import csv
import io
import json

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits

import config as C
import extinction as E
import nickel as N
import photometry as P
from gp import LightCurveGP

COSMO = FlatLambdaCDM(H0=C.H0, Om0=C.OMEGA_M)
V_PEC_KMS = 300.0  # peculiar-velocity scatter, km/s

#: CMB dipole (Planck 2018): the Sun moves at 369.82 km/s toward (l, b) =
#: (264.021, 48.253) deg relative to the cosmic microwave background.
CMB_APEX = SkyCoord(l=264.021 * u.deg, b=48.253 * u.deg, frame="galactic")
V_SUN_CMB = 369.82

#: Quality cuts (observer frame, days relative to the BTS peak time).
MIN_PRE_PEAK_G = 3
MIN_POST_PEAK_G = 5
MAX_NONDET_GAP = 5.0
MIN_EARLY_POINTS = 3
#: Low-dust subsample: peak (g - r), Milky-Way-dereddened and K-corrected, at or
#: below this. SN 2011fe itself has -0.26. Redder supernovae are mostly dimmed
#: by dust in their host galaxies (M_g grows 2.5 mag per mag of g - r here).
LOW_DUST_GR = -0.15
#: Outliers: more than this many robust standard deviations from the
#: peak-magnitude vs decline-rate (Phillips) relation of the low-dust sample.
#: These are measurement failures (a wrong redshift, a bad peak), not physics:
#: the brightest implies 1.9 solar masses of nickel, more than a white dwarf has.
OUTLIER_SIGMA = 3.0
TYPES = {"SN Ia": "normal", "SN Ia-91T": "91T-like", "SN Ia-91bg": "91bg-like"}


# ---------------------------------------------------------------------------
# SN 2011fe as the calibrator
# ---------------------------------------------------------------------------


def pereira_spectra() -> list[dict]:
    """SN 2011fe spectrophotometry (Pereira et al. 2013), dereddened for its host too."""
    folder = C.DATA_DIR / "vizier" / "J_A+A_554_A27" / "sp"
    cfg = C.SUPERNOVAE["SN2011fe"]
    out = []
    for f in sorted(folder.glob("11fe*.fit")):
        with fits.open(f) as h:
            hd = h[0].header
            flux = h[0].data.astype(float)
            wave = hd["CRVAL1"] + hd["CDELT1"] * (np.arange(hd["NAXIS1"]) + 1 - hd["CRPIX1"])
            # Milky Way dust was already removed by the SNfactory; remove the host's.
            flux = flux * 10 ** (0.4 * E.a_lambda(wave, cfg["ebv_host"], cfg["rv_host"]))
            out.append({"mjd": hd["MJD"], "phase": hd["MJD"] - cfg["t0_mjd"], "wave": wave,
                        "flam": flux, "calib_err_mag": hd["FLXERROR"]})
    return sorted(out, key=lambda s: s["mjd"])


def _gr_proxy(mag_g, mag_r):
    """lambda f_lambda summed over ZTF g and r (erg/s/cm^2): an optical-flux proxy."""
    return (P.pivot("ztf_g") * P.mag_to_flam(mag_g, "ztf_g", "AB")
            + P.pivot("ztf_r") * P.mag_to_flam(mag_r, "ztf_r", "AB"))


def calibrate_2011fe(l_peak_bol: float, t_peak_bol: float) -> dict:
    """SN 2011fe at its g-band peak: g + r flux, K-corrections and the peak offset."""
    spec = pereira_spectra()
    ph = np.array([s["phase"] for s in spec])
    g = np.array([P.synthetic_mag(s["wave"], s["flam"], "ztf_g") for s in spec])
    r = np.array([P.synthetic_mag(s["wave"], s["flam"], "ztf_r") for s in spec])
    near = (ph > 8) & (ph < 30)
    c = np.polyfit(ph[near], g[near], 3)
    tt = np.linspace(10, 28, 1801)
    t_gpk = float(tt[np.argmin(np.polyval(c, tt))])
    g_pk = float(np.polyval(c, t_gpk))
    r_at = float(np.interp(t_gpk, ph, r))
    d_cm = 10 ** (0.2 * C.SUPERNOVAE["SN2011fe"]["mu"] + 1.0) * C.PC_CM
    proxy_lum = 4 * np.pi * d_cm**2 * _gr_proxy(g_pk, r_at)

    # K-corrections from the spectrum nearest g-band peak.
    s = spec[int(np.argmin(np.abs(ph - t_gpk)))]
    zs = np.linspace(0.0, 0.06, 13)
    kg, kr = [], []
    for z in zs:
        # f_obs(lambda) = L_lambda(lambda / (1+z)) / (4 pi d_L^2 (1+z)), so
        # K = m_obs - M - mu is the magnitude of the stretched, (1+z)-dimmed
        # spectrum minus that of the original (the d_L factors cancel).
        wave_obs = s["wave"] * (1 + z)
        flam_obs = s["flam"] / (1 + z)
        kg.append(P.synthetic_mag(wave_obs, flam_obs, "ztf_g") - P.synthetic_mag(s["wave"], s["flam"], "ztf_g"))
        kr.append(P.synthetic_mag(wave_obs, flam_obs, "ztf_r") - P.synthetic_mag(s["wave"], s["flam"], "ztf_r"))
    return {"t_gpeak": t_gpk, "g_peak": g_pk, "r_at_gpeak": r_at, "gr_color": g_pk - r_at,
            "proxy_lum": proxy_lum, "l_peak_bol": l_peak_bol, "t_peak_bol": t_peak_bol,
            "bol_minus_g_peak_days": t_peak_bol - t_gpk,
            "M_g": g_pk - C.SUPERNOVAE["SN2011fe"]["mu"],
            "k_z": zs, "k_g": np.array(kg), "k_r": np.array(kr),
            "spec_phase": s["phase"], "synthetic_lc": pd.DataFrame({"phase": ph, "g": g, "r": r})}


# ---------------------------------------------------------------------------
# One ZTF supernova
# ---------------------------------------------------------------------------


def load_bts() -> pd.DataFrame:
    rows = list(csv.DictReader(io.StringIO((C.DATA_DIR / "ztf" / "bts.csv").read_text().lstrip())))
    df = pd.DataFrame(rows)
    for col in ("peakt", "peakmag", "redshift", "A_V", "rise", "fade"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["peak_mjd"] = df["peakt"] + 2458000.0 - 2400000.5
    return df


def z_cmb(z_helio: float, ra_deg: float, dec_deg: float) -> float:
    c = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
    cos_angle = np.cos(c.separation(CMB_APEX).rad)
    v_proj = V_SUN_CMB * cos_angle / 299792.458
    return float((1 + z_helio) / (1 - v_proj) - 1)  # remove the solar motion


def _first_light(t, f, fe, band_idx, nondet, t_lo, t_hi, n_grid=1201):
    """Profile chi^2 over first-light time for f_b = A_b (t - t_fl)^2.

    For each trial t_fl the amplitudes A_b follow by linear least squares.
    Non-detections (5-sigma limits) later than t_fl add a penalty when the
    model would have been above the limit. Returns t_fl, its 1-sigma
    half-width (Delta chi^2 = 1) and the minimum chi^2.
    """
    nd_b, nd_t, nd_lim = nondet
    grid = np.linspace(t_lo, t_hi, n_grid)
    x = np.clip(t[None, :] - grid[:, None], 0, None) ** 2  # (n_grid, n_points)
    w = 1.0 / fe**2
    chi = np.zeros(n_grid)
    for b in np.unique(band_idx):
        m = band_idx == b
        amp = (x[:, m] * w[m] * f[m]).sum(1) / np.clip((x[:, m] ** 2 * w[m]).sum(1), 1e-30, None)
        chi += (w[m] * (f[m] - amp[:, None] * x[:, m]) ** 2).sum(1)
        mb = nd_b == b
        if mb.any():
            model = amp[:, None] * np.clip(nd_t[mb][None, :] - grid[:, None], 0, None) ** 2
            excess = np.clip(model - nd_lim[mb][None, :], 0, None) / (nd_lim[mb][None, :] / 5.0)
            chi += (excess**2).sum(1)
    i = int(np.argmin(chi))
    ok = chi <= chi[i] + 1.0
    return float(grid[i]), float(0.5 * (grid[ok].max() - grid[ok].min())), float(chi[i])


def measure(oid: str, bts_row: pd.Series, cal: dict) -> dict:
    """All measurements for one ZTF supernova; ``reason`` says why it failed if it did."""
    z = float(bts_row["redshift"])
    out = {"oid": oid, "iau": bts_row["IAUID"], "type": TYPES.get(bts_row["type"], bts_row["type"]),
           "z": z, "ok": False}
    path = C.DATA_DIR / "ztf" / "lightcurves" / f"{oid}.json"
    if not path.exists():
        return {**out, "reason": "no light curve"}
    lc = P.load_alerce(oid)
    raw = json.loads(path.read_text())
    tp = float(bts_row["peak_mjd"])
    # Only this supernova's season: ALeRCE also holds unrelated detections
    # (host variability, earlier transients) from other years.
    lc = lc[(lc["mjd"] > tp - 45.0) & (lc["mjd"] < tp + 120.0)]
    raw["non_detections"] = [x for x in raw["non_detections"] if tp - 75.0 < x["mjd"] < tp + 120.0]
    g = lc[(lc["band"] == "ztf_g") & (lc["err"] < 0.3)]
    r = lc[(lc["band"] == "ztf_r") & (lc["err"] < 0.3)]
    gp_rel = g["mjd"].values - tp
    if (np.sum((gp_rel > -20) & (gp_rel < -1)) < MIN_PRE_PEAK_G
            or np.sum((gp_rel > 1) & (gp_rel < 30)) < MIN_POST_PEAK_G or len(r) < 5):
        return {**out, "reason": "sparse light curve"}
    t_first = float(min(g["mjd"].min(), r["mjd"].min()))
    nd = [x for x in raw["non_detections"] if x["mjd"] < t_first and x.get("diffmaglim")]
    if not nd or t_first - max(x["mjd"] for x in nd) > MAX_NONDET_GAP:
        return {**out, "reason": "no recent non-detection"}

    # Peaks from GPs (the GP works in log time; an offset keeps phases positive).
    t_ref = t_first - 30.0
    window = lambda d: d[(d["mjd"] < tp + 60)]
    gw, rw = window(g), window(r)
    try:
        gpg = LightCurveGP(gw["mjd"] - t_ref, gw["mag"], gw["err"], ell_bounds=(0.01, 0.6))
        gpr = LightCurveGP(rw["mjd"] - t_ref, rw["mag"], rw["err"], ell_bounds=(0.01, 0.6))
    except Exception:  # noqa: BLE001 - degenerate light curves
        return {**out, "reason": "GP failed"}
    tt = np.arange(max(gw["mjd"].min(), tp - 15), min(gw["mjd"].max(), tp + 15), 0.05)
    if len(tt) < 10:
        return {**out, "reason": "peak not covered"}
    mg, sg = gpg.predict(tt - t_ref)
    i = int(np.argmin(mg))
    if i == 0 or i == len(tt) - 1:
        return {**out, "reason": "peak at edge of coverage"}
    t_gpk, g_pk, g_pk_err = float(tt[i]), float(mg[i]), float(sg[i])
    t15 = t_gpk + 15.0 * (1 + z)
    if np.min(np.abs(g["mjd"].values - t15)) > 4.0 or t15 > gw["mjd"].max():
        return {**out, "reason": "no data near +15 d"}
    dm15 = float(gpg.predict([t15 - t_ref])[0][0] - g_pk)
    if not gpr.covered([t_gpk - t_ref], frac=1.0, min_days=4.0)[0]:
        return {**out, "reason": "r not covered at g peak"}
    r_at, r_err = (float(v[0]) for v in gpr.predict([t_gpk - t_ref]))

    # Dust in the Milky Way (BTS gives A_V from SFD) and K-corrections.
    a_v = float(bts_row["A_V"]) if np.isfinite(bts_row["A_V"]) else 0.0
    a_g = a_v * float(E.ccm89(P.pivot("ztf_g"))[0])
    a_r = a_v * float(E.ccm89(P.pivot("ztf_r"))[0])
    k_g = float(np.interp(z, cal["k_z"], cal["k_g"]))
    k_r = float(np.interp(z, cal["k_z"], cal["k_r"]))
    gr0 = (g_pk - a_g - k_g) - (r_at - a_r - k_r)

    # First light from the early detections in both bands.
    early_t, early_f, early_e, early_b = [], [], [], []
    for b, d in ((0, g), (1, r)):
        fl = 10 ** (-0.4 * (d["mag"].values - 23.9))  # microJansky
        fe = fl * d["err"].values * np.log(10) / 2.5
        peak_f = 10 ** (-0.4 * ((g_pk if b == 0 else r_at) - 23.9))
        sel = (d["mjd"].values < t_gpk - 3.0) & (fl < 0.5 * peak_f)
        early_t += list(d["mjd"].values[sel]); early_f += list(fl[sel])
        early_e += list(np.hypot(fe[sel], 0.03 * fl[sel])); early_b += [b] * int(sel.sum())
    if len(early_t) < MIN_EARLY_POINTS:
        return {**out, "reason": "too few early points"}
    nd_b = np.array([x["fid"] - 1 for x in nd])
    nd_t = np.array([x["mjd"] for x in nd])
    nd_f = 10 ** (-0.4 * (np.array([x["diffmaglim"] for x in nd]) - 23.9))
    recent = nd_t > t_first - 30
    t_fl, t_fl_err, chi = _first_light(np.array(early_t), np.array(early_f), np.array(early_e),
                                       np.array(early_b), (nd_b[recent], nd_t[recent], nd_f[recent]),
                                       t_first - 25.0, t_first)
    t_rise_g = (t_gpk - t_fl) / (1 + z)

    # Luminosity: SN 2011fe's bolometric peak scaled by the g + r flux ratio.
    zc = z_cmb(z, *_radec(bts_row))
    d_l = COSMO.luminosity_distance(zc).to(u.cm).value
    proxy = 4 * np.pi * d_l**2 * _gr_proxy(g_pk - a_g - k_g, r_at - a_r - k_r)
    l_peak = cal["l_peak_bol"] * proxy / cal["proxy_lum"]
    mu_err = 5 / np.log(10) * V_PEC_KMS / (299792.458 * zc)
    l_err_frac = np.hypot(0.4 * np.log(10) * np.hypot(g_pk_err, mu_err), 0.03)
    t_rise_bol = t_rise_g + cal["bol_minus_g_peak_days"]
    m_ni = N.arnett_rule_mass(l_peak, t_rise_bol)
    q = N.power(t_rise_bol)["total"]
    dq = abs(N.power(t_rise_bol + 0.5)["total"] - N.power(t_rise_bol - 0.5)["total"])
    m_ni_err = m_ni * np.hypot(l_err_frac, dq / q * t_fl_err / (1 + z))
    return {**out, "ok": True, "reason": "", "z_cmb": zc, "t_gpeak": t_gpk, "g_peak": g_pk,
            "g_peak_err": g_pk_err, "r_at_gpeak": r_at, "gr_peak": gr0, "dm15_g": dm15,
            "M_g": g_pk - a_g - k_g - COSMO.distmod(zc).value, "t_first_light": t_fl,
            "t_first_light_err": t_fl_err, "t_rise_g": t_rise_g, "t_rise_bol": t_rise_bol,
            "L_peak": l_peak, "L_peak_err": l_peak * l_err_frac, "m_ni": m_ni, "m_ni_err": m_ni_err,
            "n_early": len(early_t), "a_v_mw": a_v, "bts_rise": float(bts_row["rise"]),
            "bts_fade": float(bts_row["fade"])}


def _radec(row) -> tuple[float, float]:
    c = SkyCoord(row["RA"], row["Dec"], unit=(u.hourangle, u.deg))
    return float(c.ra.deg), float(c.dec.deg)


def run(l_peak_bol_2011fe: float, t_peak_bol_2011fe: float) -> dict:
    """Measure every candidate; returns the per-object table and the calibration."""
    import fetch_data

    cal = calibrate_2011fe(l_peak_bol_2011fe, t_peak_bol_2011fe)
    bts = load_bts()
    cands = {r["ZTFID"] for r in fetch_data.bts_candidates(C.DATA_DIR / "ztf" / "bts.csv")}
    rows = []
    for _, row in bts[bts["ZTFID"].isin(cands)].iterrows():
        rows.append(measure(row["ZTFID"], row, cal))
    table = pd.DataFrame(rows)
    good = table[table["ok"]].copy()
    good["low_dust"] = good["gr_peak"] <= LOW_DUST_GR
    good["outlier"] = False
    ld = good["low_dust"].values
    keep = ld.copy()
    for _ in range(5):  # iterate: fit, clip at 3 robust sigma, refit
        p = np.polyfit(good["dm15_g"][keep], good["M_g"][keep], 1)
        res = good["M_g"] - np.polyval(p, good["dm15_g"])
        sig = 1.4826 * np.median(np.abs(res[keep] - np.median(res[keep])))
        new_keep = ld & (np.abs(res) <= OUTLIER_SIGMA * sig)
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep
    good.loc[ld & ~keep, "outlier"] = True
    good["phillips_resid"] = res
    return {"all": table, "good": good, "calibration": cal,
            "phillips": {"slope": float(p[0]), "intercept": float(p[1]), "sigma": float(sig)}}
