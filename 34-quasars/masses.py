"""Single-epoch black hole masses and Eddington ratios, recomputed from line widths.

The broad emission lines come from gas orbiting the black hole. Reverberation
mapping (project I2) found that the size of that region grows with the quasar's
luminosity, R ~ L^0.5, so one spectrum gives both a size (from L) and a speed
(from the line width), and the mass follows from M = f R v^2 / G. Packaged as

    log M = a + b log(lambda L_lambda / 1e44 erg/s) + 2 log(FWHM / km/s)

with (a, b) calibrated on reverberation-mapped AGN. Different lines are used at
different redshifts because only ~3,600-10,400 A is observed: H beta up to
z ~ 0.7-0.9, Mg II at 0.4 < z < 2.6, C IV beyond z ~ 1.4.

We recompute Wu & Shen's masses from their published widths and luminosities
(a check of the bookkeeping), then test the recipes against each other where
two lines are seen in the same spectrum.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import physics
from config import (BOLOMETRIC_CORRECTION, CIV_REST, COATMAN17, FWHM_RANGE, MASS_ERR_MAX,
                    VIRIAL)

C_KMS = 299792.458


def virial_mass(fwhm_kms, log_lambda_l, recipe: str):
    """log10 M/Msun from a line FWHM (km/s) and log10 lambda L_lambda (erg/s)."""
    a, b, _ = VIRIAL[recipe]
    fwhm = np.asarray(fwhm_kms, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return a + b * (np.asarray(log_lambda_l, float) - 44.0) + 2.0 * np.log10(fwhm)


def bolometric(log_lambda_l, wave: int):
    """log10 L_bol from a monochromatic luminosity with the Richards et al. 2006 corrections."""
    return np.asarray(log_lambda_l, float) + np.log10(BOLOMETRIC_CORRECTION[wave])


def eddington_ratio(log_lbol, log_mass, per_msun: float | None = None):
    per_msun = physics.eddington_per_msun() if per_msun is None else per_msun
    return np.asarray(log_lbol, float) - np.log10(per_msun) - np.asarray(log_mass, float)


def civ_blueshift_kms(civ_half_rest, z_fit, z_sys):
    """C IV blueshift (km/s, positive = towards the observer) relative to a systemic redshift.

    ``civ_half_rest`` is the wavelength bisecting the C IV flux, measured in the
    rest frame of ``z_fit``; it is moved to the frame of ``z_sys`` first.
    """
    lam = np.asarray(civ_half_rest, float) * (1 + np.asarray(z_fit, float)) / (1 + np.asarray(z_sys, float))
    return C_KMS * (CIV_REST - lam) / CIV_REST


def coatman_corrected_mass(fwhm_civ, blueshift_kms, log_l1350):
    """Coatman et al. 2017 eqs. 4 and 6: C IV mass with the wind-broadening removed."""
    denom = COATMAN17["slope"] * np.asarray(blueshift_kms, float) / 1e3 + COATMAN17["intercept"]
    fwhm_corr = np.asarray(fwhm_civ, float) / denom
    with np.errstate(divide="ignore", invalid="ignore"):
        return (COATMAN17["a"] + COATMAN17["b"] * (np.asarray(log_l1350, float) - 44.0)
                + 2.0 * np.log10(fwhm_corr))


def fit_civ_recalibration(blueshift_kms, fwhm_civ, fwhm_mgii, log_l1350, log_m_mgii, n_bins: int = 12):
    """Our own C IV recipe, put on the Mg II scale with quasars that show both lines.

    Like Coatman et al. (2017), but calibrated on DR16Q against Mg II instead of on
    luminous quasars against H alpha: (1) the C IV/Mg II width ratio is fitted as a
    straight line in blueshift, through the medians of blueshift bins; (2) the
    zero point a of log M = a + 0.53 log(L1350/1e44) + 2 log(FWHM_corr) is the
    median offset that makes the corrected C IV masses agree with Mg II.
    Returns dict(slope, intercept, a).
    """
    bs = np.asarray(blueshift_kms, float)
    ratio = np.asarray(fwhm_civ, float) / np.asarray(fwhm_mgii, float)
    good = np.isfinite(bs) & np.isfinite(ratio)
    bs, ratio = np.where(good, bs, np.nan), np.where(good, ratio, np.nan)
    edges = np.nanpercentile(bs, np.linspace(1, 99, n_bins + 1))
    idx = np.digitize(bs, edges) - 1
    xs, ys = [], []
    for k in range(n_bins):
        m = idx == k
        if m.sum() >= 30:
            xs.append(np.median(bs[m]) / 1e3)
            ys.append(np.median(ratio[m]))
    slope, intercept = np.polyfit(xs, ys, 1)
    denom = slope * bs / 1e3 + intercept
    with np.errstate(divide="ignore", invalid="ignore"):
        fwhm_corr = np.where(denom > 0.2, np.asarray(fwhm_civ, float) / denom, np.nan)
        base = 0.53 * (np.asarray(log_l1350, float) - 44.0) + 2.0 * np.log10(fwhm_corr)
    a = float(np.nanmedian(np.asarray(log_m_mgii, float) - base))
    return dict(slope=float(slope), intercept=float(intercept), a=a)


def recalibrated_civ_mass(blueshift_kms, fwhm_civ, log_l1350, cal: dict):
    denom = cal["slope"] * np.asarray(blueshift_kms, float) / 1e3 + cal["intercept"]
    with np.errstate(divide="ignore", invalid="ignore"):
        fwhm_corr = np.where(denom > 0.2, np.asarray(fwhm_civ, float) / denom, np.nan)
        return cal["a"] + 0.53 * (np.asarray(log_l1350, float) - 44.0) + 2.0 * np.log10(fwhm_corr)


def mass_error(fwhm, fwhm_err, logl_err, b):
    """Measurement error of log M from the width and luminosity errors (as Wu & Shen do:
    our estimate matches their LOGMBH_*_ERR wherever they give one)."""
    fwhm = np.asarray(fwhm, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.sqrt((2.0 * np.asarray(fwhm_err, float) / (fwhm * np.log(10))) ** 2
                       + (b * np.asarray(logl_err, float)) ** 2)


def _line_ok(fwhm, flux, logl, err=None):
    ok = (fwhm > FWHM_RANGE[0]) & (fwhm < FWHM_RANGE[1]) & (flux > 0) & (logl > 40)
    if err is not None:
        ok &= (err > 0) & (err < MASS_ERR_MAX)
    return ok


def recompute(ws: pd.DataFrame) -> pd.DataFrame:
    """Masses, bolometric luminosities and Eddington ratios for every quasar.

    Columns of the result (NaN where a line is unusable):
      m_hb, m_mgii, m_civ: our recomputation with Wu & Shen's recipes (VP06 H beta,
          Shen 2011 Mg II, VP06 C IV) from their FWHMs and continuum luminosities;
      m_mgii_vo09: Mg II with the Vestergaard & Osmer 2009 recipe;
      civ_blueshift: C IV blueshift relative to the Mg II systemic redshift;
      m_civ_c17: Coatman-corrected C IV mass where that blueshift is known;
      log_lbol: from L3000 where measured, else L5100 (z < 0.7) or L1350;
      m_fid, line_fid: the fiducial mass, by Wu & Shen's redshift ranges;
      log_edd: log L_bol / L_Edd with L_Edd = 1.257e38 M (pure hydrogen).
    """
    z = (ws["zbest"] if "zbest" in ws else ws["Z_SYS"]).to_numpy()
    out = pd.DataFrame(index=ws.index)
    hb_fwhm, hb_flux = ws["HBETA_BR_4"].to_numpy(), ws["HBETA_BR_2"].to_numpy()
    mg_fwhm, mg_flux = ws["MGII_BR_4"].to_numpy(), ws["MGII_BR_2"].to_numpy()
    c4_fwhm, c4_flux = ws["CIV_4"].to_numpy(), ws["CIV_2"].to_numpy()
    l5100, l3000, l1350 = (ws["LOGL5100"].to_numpy(), ws["LOGL3000"].to_numpy(),
                           ws["LOGL1350"].to_numpy())

    # Mass errors computed here rather than read from the catalogue: identical where
    # Wu & Shen give one, and available where they did not compute a mass (Mg II
    # at 2.0 < z < 2.4 in BOSS spectra, beyond their switch to C IV).
    e_hb = mass_error(hb_fwhm, ws["HBETA_BR_ERR_4"], ws["LOGL5100_ERR"], VIRIAL["HB_VP06"][1])
    e_mg = mass_error(mg_fwhm, ws["MGII_BR_ERR_4"], ws["LOGL3000_ERR"], VIRIAL["MGII_S11"][1])
    e_c4 = mass_error(c4_fwhm, ws["CIV_ERR_4"], ws["LOGL1350_ERR"], VIRIAL["CIV_VP06"][1])
    ok_hb = _line_ok(hb_fwhm, hb_flux, l5100, e_hb)
    ok_mg = _line_ok(mg_fwhm, mg_flux, l3000, e_mg)
    ok_c4 = _line_ok(c4_fwhm, c4_flux, l1350, e_c4)
    out["m_hb"] = np.where(ok_hb, virial_mass(hb_fwhm, l5100, "HB_VP06"), np.nan)
    out["m_mgii"] = np.where(ok_mg, virial_mass(mg_fwhm, l3000, "MGII_S11"), np.nan)
    out["m_mgii_vo09"] = np.where(ok_mg, virial_mass(mg_fwhm, l3000, "MGII_VO09"), np.nan)
    out["m_civ"] = np.where(ok_c4, virial_mass(c4_fwhm, l1350, "CIV_VP06"), np.nan)

    # C IV blueshift against the Mg II-based systemic redshift (Mg II sits within
    # ~100 km/s of systemic; Wu & Shen's ZSYS_LINES_4 includes that small offset).
    z_mg = ws["ZSYS_LINES_4"].to_numpy()
    has_zmg = ok_mg & (z_mg > 0) & (ws["ZSYS_LINES_ERR_4"].to_numpy() > 0)
    bs = civ_blueshift_kms(ws["CIV_1"].to_numpy(), ws["Z_FIT"].to_numpy(), z_mg)
    out["civ_blueshift"] = np.where(has_zmg & ok_c4, bs, np.nan)
    with np.errstate(invalid="ignore"):
        # Coatman's denominator reaches zero at -1,500 km/s: the recipe is not
        # meant for redshifted C IV, so it is applied only above -1,000 km/s.
        c17_ok = np.isfinite(out["civ_blueshift"]) & (out["civ_blueshift"] > -1000)
    out["m_civ_c17"] = np.where(c17_ok, coatman_corrected_mass(c4_fwhm, out["civ_blueshift"], l1350),
                                np.nan)
    out["fwhm_civ"] = np.where(ok_c4, c4_fwhm, np.nan)
    out["fwhm_mgii"] = np.where(ok_mg, mg_fwhm, np.nan)
    out["log_l1350"] = np.where(l1350 > 40, l1350, np.nan)

    # Bolometric luminosity: L3000 preferred (least host light, least reddening),
    # as Wu & Shen intended (their catalogue used L5100 at z < 0.7 by mistake).
    lb = np.full(len(ws), np.nan)
    for wave, arr in ((1350, l1350), (5100, l5100), (3000, l3000)):  # later wins
        good = arr > 40
        lb[good] = bolometric(arr[good], wave)
    out["log_lbol"] = lb

    boss = (ws["IF_BOSS_SDSS"] == "BOSS").to_numpy()
    z_civ = np.where(boss, 2.0, 1.9)
    fid = np.where(z < 0.7, out["m_hb"], np.where(z < z_civ, out["m_mgii"], out["m_civ"]))
    line = np.where(z < 0.7, "Hbeta", np.where(z < z_civ, "MgII", "CIV"))
    fid = np.where(z > 0, fid, np.nan)
    out["m_fid"] = fid
    out["line_fid"] = np.where(np.isfinite(fid), line, "")
    out["log_edd"] = eddington_ratio(out["log_lbol"], out["m_fid"])
    return out
