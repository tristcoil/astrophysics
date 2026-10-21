"""The Leavitt law in the Magellanic Clouds: one master table per galaxy and the fits.

The master table starts from OGLE's catalogue (modes, periods, positions) and
replaces OGLE's mean magnitudes with our own intensity means from the raw
light curves wherever our fit succeeded (the two agree to ~0.001 mag; see
``lightcurve_check``). Gaia DR3 photometry of the same stars is attached.
"""

from __future__ import annotations

import numpy as np

import config as C
import gaia_cal
import lightcurves as LC
import ogle
import plfit

#: Our intensity means are used if the light-curve fit kept at least this many points.
MIN_POINTS_I = 60
MIN_POINTS_V = 20


def master() -> dict:
    """{'lmc': table, 'smc': table}: OGLE catalogue + our light-curve fits + Gaia photometry."""
    fits = LC.load_fits()
    mc = gaia_cal.magellanic()
    out = {}
    for gal in C.GALAXIES:
        cat = dict(ogle.catalog(gal))
        n = len(cat["id"])
        ours = {k: np.full(n, np.nan) for k in ("I_mean", "V_mean", "I_P", "I_eP", "I_amp", "V_amp", "I_max",
                                                 "I_min", "V_max", "V_min", "I_n", "V_n", "I_R21", "I_phi21",
                                                 "I_R31", "I_phi31", "I_T_max", "V_T_max", "I_P_blind", "V_P")}
        for j, sid in enumerate(cat["id"]):
            r = fits.get(sid)
            if r is None:
                continue
            for key in ours:
                v = r.get(key, np.nan)
                ours[key][j] = np.nan if v is None else v
        # our period must agree with OGLE's (it fails for ~0.2%: daily aliases in
        # sparse light curves, and the eclipsing binary OGLE-LMC-CEP-0227)
        same_p = np.abs(ours["I_P"] - cat["P"]) / cat["P"] < 1e-3
        good_i = (ours["I_n"] >= MIN_POINTS_I) & same_p
        good_v = (ours["V_n"] >= MIN_POINTS_V) & same_p
        cat["I_ogle"], cat["V_ogle"] = cat["I"].copy(), cat["V"].copy()
        cat["I"] = np.where(good_i & np.isfinite(ours["I_mean"]), ours["I_mean"], cat["I_ogle"])
        cat["V"] = np.where(good_v & np.isfinite(ours["V_mean"]), ours["V_mean"], cat["V_ogle"])
        cat["phot_source"] = np.where(good_i & good_v, "ours", "ogle")
        cat["W_I"] = cat["I"] - C.R_WI * (cat["V"] - cat["I"])
        cat["W_I_ogle"] = cat["I_ogle"] - C.R_WI * (cat["V_ogle"] - cat["I_ogle"])
        cat.update({f"lc_{k}": v for k, v in ours.items()})
        cat.update({f"gaia_{k}": v for k, v in mc[gal].items()})
        cat["VI"] = cat["V"] - cat["I"]
        out[gal] = cat
    return out


def lightcurve_check(tab: dict) -> dict:
    """How our blind periods and intensity means compare with OGLE's catalogue (single-mode F and 1O)."""
    out = {}
    for mode in ("F", "1O"):
        s = (tab["mode"] == mode) & np.isfinite(tab["lc_I_P"]) & (tab["lc_I_n"] >= MIN_POINTS_I)
        rel = (tab["lc_I_P"][s] - tab["P"][s]) / tab["P"][s]
        blind_ok = np.abs(tab["lc_I_P_blind"][s] - tab["P"][s]) / tab["P"][s] < 1e-3
        dI = tab["lc_I_mean"][s] - tab["I_ogle"][s]
        sv = s & np.isfinite(tab["lc_V_mean"]) & np.isfinite(tab["V_ogle"]) & (tab["lc_V_n"] >= MIN_POINTS_V)
        dV = tab["lc_V_mean"][sv] - tab["V_ogle"][sv]
        out[mode] = {"n": int(s.sum()), "blind_period_recovered": float(np.mean(blind_ok)),
                     "n_blind_failed": int((~blind_ok).sum()),
                     "median_abs_dP_over_P": float(np.median(np.abs(rel))),
                     "frac_dP_over_P_below_1e-4": float(np.mean(np.abs(rel) < 1e-4)),
                     "dI_median": float(np.nanmedian(dI)), "dI_mad": float(1.4826 * np.nanmedian(np.abs(dI - np.nanmedian(dI)))),
                     "dV_median": float(np.nanmedian(dV)), "dV_mad": float(1.4826 * np.nanmedian(np.abs(dV - np.nanmedian(dV)))),
                     "n_V": int(sv.sum())}
    return out


def pl_table(tab: dict, modes=("F", "1O"), bands=("V", "I", "W_I"), keep=None) -> dict:
    """Clipped straight-line Leavitt laws per mode and band."""
    out = {}
    for mode in modes:
        s = tab["mode"] == mode
        if keep is not None:
            s &= keep
        for band in bands:
            f = plfit.fit_line(tab["logP"], tab[band], keep=s)
            out[(mode, band)] = f
    return out
