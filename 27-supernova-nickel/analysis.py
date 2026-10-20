"""Every measurement in the project: bolometric light curves and nickel masses.

    import analysis; r = analysis.run()

Steps, per supernova:
1. Bolometric light curve (``bolometric.py``): UV+optical integrated from the
   photometry; near-infrared measured where JHKs exist, otherwise added with
   the infrared fraction measured at other epochs (see ``nir_template_*``).
2. Type Ia: peak luminosity and Arnett's rule; a fit of Arnett's full model;
   the Katz integral; and the radioactive tail (nickel mass and the gamma-ray
   escape time T0).
3. Type II: the radioactive tail.
4. SN 2014J: nickel from INTEGRAL's gamma-ray lines (``gammarays.py``).
The ZTF sample (Arnett's rule for hundreds of Type Ia) is in ``ztf.py``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import bolometric as B
import config as C
import gammarays as G
import nickel as N
import literature as LIT
import photometry as P
import ztf as Z

N_MC = 200


def _grid(t_lo: float, t_hi: float, n: int) -> np.ndarray:
    return np.unique(np.round(np.geomspace(t_lo, t_hi, n), 2))


def distance_cm(sn: str) -> float:
    return 10 ** (0.2 * C.SUPERNOVAE[sn]["mu"] + 1.0) * C.PC_CM


def distance_scale_err(sn: str) -> float:
    """Fractional luminosity error from the distance alone (common to all epochs)."""
    return 0.4 * np.log(10) * C.SUPERNOVAE[sn]["mu_err"]


# ---------------------------------------------------------------------------
# Near-infrared templates
# ---------------------------------------------------------------------------

#: Assumed nebular-phase near-infrared colours (Vega) where only J was observed:
#: J - H = 0 +/- 0.3 and H - Ks = 0 +/- 0.3. SN 2011fe itself had J - H = 0.07
#: at day 224 (Weyant et al. 2017); late SNe Ia are dominated by [Fe II] lines
#: at 1.26 and 1.64 micron, which keeps J and H comparable.
NEBULAR_J_H, NEBULAR_H_K, NEBULAR_COLOUR_ERR = 0.0, 0.0, 0.3


def _anchor_from_j(bol: B.Bolometric, t: float, j_mag: float, j_err: float, n: int = 400,
                   seed: int = 3) -> tuple[float, float]:
    """Near-infrared fraction at an epoch with only a J measurement.

    Integrates from the reddest optical band through J and assumed H, Ks
    (colours above) plus a Rayleigh-Jeans tail; Monte Carlo over the J error
    and the colour assumptions.
    """
    rng = np.random.default_rng(seed)
    row = bol.luminosity([t], n_mc=50).iloc[0]
    red_band = row["bands"].split(",")
    red_band = [b for b in red_band if P.pivot(b) < B.NIR_MIN_WAVE][-1]
    m_red, s_red = bol.gps[red_band].predict(t)
    d2 = 4 * np.pi * bol.dist_cm**2
    fracs = []
    for _ in range(n):
        mj = j_mag + j_err * rng.standard_normal()
        mh = mj - (NEBULAR_J_H + NEBULAR_COLOUR_ERR * rng.standard_normal())
        mk = mh - (NEBULAR_H_K + NEBULAR_COLOUR_ERR * rng.standard_normal())
        mr = m_red[0] + s_red[0] * rng.standard_normal()
        lam = np.array([P.pivot(red_band), P.pivot("J"), P.pivot("H"), P.pivot("Ks")])
        f = np.array([P.mag_to_flam(mr - bol.extinction(lam[0]), red_band, bol.system[red_band]),
                      P.mag_to_flam(mj - bol.extinction(lam[1]), "J", "Vega"),
                      P.mag_to_flam(mh - bol.extinction(lam[2]), "H", "Vega"),
                      P.mag_to_flam(mk - bol.extinction(lam[3]), "Ks", "Vega")])
        l_nir = (np.trapezoid(f, lam) + f[-1] * lam[-1] / 3.0) * d2
        fracs.append(l_nir / (row["L_uvoi"] + l_nir))
    return float(np.mean(fracs)), float(np.std(fracs))


def nir_template_2011fe(bol: B.Bolometric, table: pd.DataFrame) -> dict:
    """Near-infrared fraction of SN 2011fe vs time: measured to day 70, LBT J at day 253."""
    meas = B.nir_fraction(table)
    meas = meas[meas["phase"] <= 75]
    shappee = P.load_shappee2017()
    j = shappee[(shappee["band"] == "J") & (shappee["phase"] < 300)].iloc[0]
    f253, e253 = _anchor_from_j(bol, float(j["phase"]), float(j["mag"]), float(j["err"]))
    # Smooth the measured fraction a little (GP interpolation wiggles at the few-% level).
    ph = meas["phase"].values
    fv = pd.Series(meas["f_nir"].values).rolling(5, center=True, min_periods=1).median().values
    phase = np.concatenate([ph, [float(j["phase"])]])
    value = np.concatenate([fv, [f253]])
    # Measured epochs: +/-0.02 (photometric systems, trapezoid between bands).
    # Between day 70 and 253 nothing was observed in the near-infrared; the
    # straight-line interpolation is given +/-0.07 at the middle of the gap.
    err = np.concatenate([np.full(len(ph), 0.02), [max(e253, 0.03)]])
    gap = np.array([0.5 * (ph[-1] + phase[-1])])
    phase = np.concatenate([phase, gap])
    value = np.concatenate([value, [np.interp(gap[0], [ph[-1], phase[-2]], [fv[-1], f253])]])
    err = np.concatenate([err, [0.07]])
    order = np.argsort(phase)
    return {"phase": phase[order], "value": value[order], "err": err[order],
            "measured": meas, "anchor": {"phase": float(j["phase"]), "f": f253, "err": e253}}


# ---------------------------------------------------------------------------
# Bolometric light curves
# ---------------------------------------------------------------------------


def bolometric_2011fe() -> dict:
    phot = P.load_supernova("SN2011fe")
    phot = phot[~((phot["band"] == "J") & phot["source"].str.contains("Shappee"))]
    bol = B.Bolometric("SN2011fe", phot)
    grid = np.concatenate([np.arange(1.5, 40.0, 0.5), _grid(40.0, 360.0, 70)])
    table = bol.luminosity(grid, n_mc=N_MC)
    tmpl = nir_template_2011fe(bol, table)
    table = B.apply_nir_fraction(table, tmpl["phase"], tmpl["value"], tmpl["err"])
    return {"table": table, "nir": tmpl, "bol": bol}


#: Alternative dust toward SN 2014J: A_V = 1.85, R_V = 1.6, the values Churazov
#: et al. (2014) used for their optical luminosity (E(B-V) = A_V / R_V).
SN2014J_ALT_DUST = (1.85 / 1.6, 1.6)


def bolometric_2014j(tmpl_2011fe: dict, host_dust=None) -> dict:
    """SN 2014J: optical + its own JHK to day 48; beyond that, SN 2011fe's infrared fraction.

    Swift UV is left out: through E(B-V) = 1.37 of R_V = 1.4 dust the UV
    correction is 7-10 magnitudes and meaningless.
    """
    phot = P.load_supernova("SN2014J")
    bands = ["U", "B", "V", "R", "I", "J", "H", "Ks"]
    bol = B.Bolometric("SN2014J", phot, bands=bands, host_dust=host_dust)
    grid = np.concatenate([np.arange(6.0, 40.0, 0.5), _grid(40.0, 128.0, 30)])
    table = bol.luminosity(grid, n_mc=N_MC)
    own = B.nir_fraction(table)
    last = own["phase"].max()
    # Own measurements up to day ~48; SN 2011fe's template (+0.03 extra
    # uncertainty for using a twin) afterwards.
    t_ph, t_v, t_e = tmpl_2011fe["phase"], tmpl_2011fe["value"], tmpl_2011fe["err"]
    after = t_ph > last
    phase = np.concatenate([own["phase"].values, t_ph[after]])
    value = np.concatenate([own["f_nir"].values, t_v[after]])
    err = np.concatenate([np.full(len(own), 0.02), np.hypot(t_e[after], 0.03)])
    table = B.apply_nir_fraction(table, phase, value, err)
    return {"table": table, "nir": {"phase": phase, "value": value, "err": err, "measured": own},
            "bol": bol}


def bolometric_2023ixf() -> dict:
    phot = P.load_supernova("SN2023ixf")
    # The first weeks of a Type II are hot (>12,000 K): blackbody UV tail then.
    bol = B.Bolometric("SN2023ixf", phot, hot_uv_days=20.0)
    grid = np.concatenate([np.arange(1.5, 30.0, 0.5), _grid(30.0, 362.0, 90)])
    table = bol.luminosity(grid, n_mc=N_MC)
    meas = B.nir_fraction(table)
    late = meas[meas["phase"] >= 120]
    f_late, f_late_err = float(late["f_nir"].mean()), 0.10
    # Measured to day ~145; afterwards held at the day 120-145 mean, +/-0.10
    # (dust formed in the ejecta around then, Singh et al. 2024, which could
    # push light further into the infrared).
    phase = np.concatenate([meas["phase"].values, [400.0]])
    value = np.concatenate([meas["f_nir"].values, [f_late]])
    err = np.concatenate([np.full(len(meas), 0.02), [f_late_err]])
    after = meas["phase"].max()
    err = np.where(phase > after, f_late_err, err)
    table = B.apply_nir_fraction(table, phase, value, err)
    return {"table": table, "nir": {"phase": phase, "value": value, "err": err, "measured": meas,
                                    "f_late": f_late}, "bol": bol}


def bolometric_1987a() -> dict:
    """SN 1987A: the published UVOIR curve (rescaled to our distance) plus our own UBVRI check."""
    scale = (10 ** (0.2 * C.SUPERNOVAE["SN1987A"]["mu"] + 1.0) / 1e3
             / C.SN1987A_PUBLISHED_DISTANCE_KPC) ** 2
    pub = P.load_sn1987a_uvoir()
    pub["L"] = 10 ** pub["log_l"] * scale
    # Up to day 600 the tabulated values are good to 0.03 dex absolute, 0.01 dex
    # relative (Suntzeff et al. 1991, as noted by Jeffery); 0.05 dex later.
    pub["L_err"] = pub["L"] * np.log(10) * np.where(pub["phase"] < 600, 0.01, 0.03)
    late = P.load_sn1987a_late_uvoir()
    late["L"] = 10 ** late["log_l"] * scale
    late["L_err"] = late["L"] * np.log(10) * 0.03
    phot = P.load_supernova("SN1987A")
    bol = B.Bolometric("SN1987A", phot, hot_uv_days=5.0)
    grid = _grid(1.0, 1500.0, 150)
    ours = bol.luminosity(grid, n_mc=N_MC)
    # Infrared fraction implied by the published UVOIR vs our UBVRI integral.
    both = pd.concat([pub[["phase", "L"]], late[late["phase"] > pub["phase"].max()][["phase", "L"]]])
    both = both.sort_values("phase")
    ours["L_pub_interp"] = np.interp(ours["phase"], both["phase"], both["L"], left=np.nan, right=np.nan)
    ours["f_outside_ubvri"] = 1.0 - ours["L_uvoi"] / ours["L_pub_interp"]
    return {"published": pub, "published_late": late, "ours": ours, "bol": bol,
            "distance_rescale": scale}


# ---------------------------------------------------------------------------
# Measurements
# ---------------------------------------------------------------------------


def peak(table: pd.DataFrame, t_lo: float, t_hi: float) -> dict:
    """Peak of the bolometric light curve: parabola through the points near the maximum."""
    d = table[(table["phase"] >= t_lo) & (table["phase"] <= t_hi)]
    i = int(np.argmax(d["L"].values))
    t_i = d["phase"].values[i]
    near = d[np.abs(d["phase"] - t_i) <= 3.0]
    c = np.polyfit(near["phase"], np.log10(near["L"]), 2)
    t_pk = -c[1] / (2 * c[0])
    l_pk = 10 ** np.polyval(c, t_pk)
    return {"t_peak": float(t_pk), "L_peak": float(l_pk),
            "L_peak_err_phot": float(np.interp(t_pk, d["phase"], d["L_err"]))}


def type_ia_measurements(sn: str, table: pd.DataFrame, fit_until: float) -> dict:
    """Arnett's rule, Arnett's model with leakage, the tail and the Katz integral."""
    cfg = C.SUPERNOVAE[sn]
    out = {"sn": sn, "distance_frac_err": distance_scale_err(sn), "t0_ref": cfg["t0_ref"]}
    pk = peak(table, 10.0, 30.0)
    out["peak"] = pk
    l_err = pk["L_peak"] * np.hypot(pk["L_peak_err_phot"] / pk["L_peak"], out["distance_frac_err"])
    m_rule = N.arnett_rule_mass(pk["L_peak"], pk["t_peak"])
    out["arnett_rule"] = {"m_ni": m_rule, "m_ni_err_phot_dist": m_rule * l_err / pk["L_peak"],
                          "alpha": 1.0, "t_rise": pk["t_peak"]}

    # Primary: Arnett's model with gamma-ray leakage, fitted from first light
    # through the early tail.
    af = N.fit_arnett(table["phase"], table["L"], table["L_err"], t_max_fit=fit_until)
    out["arnett_fit"] = af

    # Tail alone with T0 fixed at the Arnett-fit value: does the tail agree?
    tail = table[(table["phase"] >= 60.0) & (table["phase"] <= fit_until)]
    tf = N.fit_tail(tail["phase"], tail["L"], tail["L_err"], fit_t0=False,
                    full_trapping_t0=af["t0"])
    tf["window"] = [60.0, fit_until]
    out["tail"] = tf

    # Katz integral from first light to the end of the fit window.
    kz = table[table["phase"] <= fit_until]
    out["katz"] = {"m_ni": N.katz_ratio(kz["phase"].values, kz["L"].values, af["t0"]),
                   "t_max": float(kz["phase"].max())}

    # The late tail: nickel mass implied epoch by epoch for the fitted T0
    # with positrons fully trapped. A flat line means the model holds.
    late = table[table["phase"] >= 40.0]
    out["m_inst"] = {"phase": late["phase"].values,
                     "m": (late["L"] / N.deposition(late["phase"].values, af["t0"])).values,
                     "err": (late["L_err"] / N.deposition(late["phase"].values, af["t0"])).values}
    return out


def type_ii_measurements(sn: str, phase, lum, lum_err, windows: dict) -> dict:
    out = {"sn": sn, "distance_frac_err": distance_scale_err(sn)}
    for name, (lo, hi, fit_t0) in windows.items():
        sel = (phase >= lo) & (phase <= hi)
        tf = N.fit_tail(phase[sel], lum[sel], lum_err[sel], fit_t0=fit_t0, t0_guess=200.0)
        tf["window"] = [lo, hi]
        out[name] = tf
    return out


def gamma_2014j(t0_optical: float | None) -> dict:
    d = distance_cm("SN2014J")
    out = {"transparent": G.nickel_from_lines(d, None)}
    for prof in ("uniform", "exponential"):
        r = G.nickel_from_lines(d, prof)
        ej = r.pop("ejecta")
        r["t0_model_kappa025"] = ej.t0_days(0.025)
        r["t0_model_kappa030"] = ej.t0_days(0.03)
        r["v_scale_kms"] = ej.v_scale / 1e5
        r["ni_shell_kms"] = (ej.v_in / 1e5, ej.v_out / 1e5)
        out[prof] = r
    out["t0_optical"] = t0_optical
    out["distance_frac_err"] = distance_scale_err("SN2014J")
    out["diehl"] = G.diehl_points()
    return out


def run() -> dict:
    res = {}
    fe = bolometric_2011fe()
    res["SN2011fe"] = {"bolo": fe, "fit": type_ia_measurements("SN2011fe", fe["table"], 150.0)}
    j14 = bolometric_2014j(fe["nir"])
    res["SN2014J"] = {"bolo": j14, "fit": type_ia_measurements("SN2014J", j14["table"], 128.0)}
    # Robustness: (a) the other published dust law; (b) the peak only, with
    # T0 borrowed from SN 2011fe, since 2014J's late optical tail is the least
    # trustworthy part (bright M82 background, broad-band extinction effects).
    alt = bolometric_2014j(fe["nir"], host_dust=SN2014J_ALT_DUST)
    pk_alt = peak(alt["table"], 10.0, 30.0)
    t0_fe = res["SN2011fe"]["fit"]["arnett_fit"]["t0"]
    early = j14["table"][j14["table"]["phase"] <= 60.0]
    res["SN2014J"]["variants"] = {
        "alt_dust": {"host_dust": SN2014J_ALT_DUST, "peak": pk_alt,
                     "arnett_rule_m_ni": N.arnett_rule_mass(pk_alt["L_peak"], pk_alt["t_peak"]),
                     "table": alt["table"]},
        "peak_only_t0_from_2011fe": {
            **N.fit_arnett_fixed_t0(early["phase"], early["L"], early["L_err"], t0_fe),
            "t0_fixed": t0_fe},
    }

    ixf = bolometric_2023ixf()
    t = ixf["table"]
    res["SN2023ixf"] = {"bolo": ixf, "fit": type_ii_measurements(
        "SN2023ixf", t["phase"].values, t["L"].values, t["L_err"].values,
        {"tail_measured_nir": (95.0, 146.0, True),
         "tail_full_year": (95.0, 362.0, True),
         "tail_full_trapping_early": (95.0, 146.0, False)})}

    sn87 = bolometric_1987a()
    pub = pd.concat([sn87["published"], sn87["published_late"][sn87["published_late"]["phase"] > 935]])
    res["SN1987A"] = {"bolo": sn87, "fit": type_ii_measurements(
        "SN1987A", pub["phase"].values, pub["L"].values, pub["L_err"].values,
        {"tail_full_trapping": (130.0, 300.0, False),
         "tail_with_escape": (130.0, 600.0, True)})}

    res["gamma"] = gamma_2014j(res["SN2014J"]["fit"]["arnett_fit"]["t0"])
    pk = res["SN2011fe"]["fit"]["peak"]
    res["ztf"] = Z.run(pk["L_peak"], pk["t_peak"])
    res["summary"] = summary(res)
    return res


def _pm(v, e, nd=3):
    return {"value": round(float(v), nd), "err": None if e is None else round(float(e), nd)}


def summary(res: dict) -> dict:
    """The headline numbers, JSON-ready, with their statistical and distance errors."""
    out = {}
    for sn in ("SN2011fe", "SN2014J"):
        f = res[sn]["fit"]
        de = f["distance_frac_err"]
        af = f["arnett_fit"]
        out[sn] = {
            "L_peak_1e43": _pm(f["peak"]["L_peak"] / 1e43,
                               f["peak"]["L_peak"] / 1e43 * np.hypot(f["peak"]["L_peak_err_phot"] / f["peak"]["L_peak"], de)),
            "t_peak_days": round(f["peak"]["t_peak"], 2),
            "m_ni_arnett_rule": _pm(f["arnett_rule"]["m_ni"], f["arnett_rule"]["m_ni_err_phot_dist"]),
            "m_ni_arnett_fit": _pm(af["m_ni"], af["m_ni"] * np.hypot(af["m_ni_err"] / af["m_ni"], de)),
            "t0_days": _pm(af["t0"], af["t0_err"], 1),
            "tau_m_days": _pm(af["tau_m"], af["tau_m_err"], 2),
            "m_ni_tail_t0_fixed": _pm(f["tail"]["m_ni"], f["tail"]["m_ni"] * np.hypot(f["tail"]["m_ni_err"] / f["tail"]["m_ni"], de)),
            "m_ni_katz": round(f["katz"]["m_ni"], 3),
            "distance_frac_err_L": round(de, 3),
        }
    v = res["SN2014J"]["variants"]
    out["SN2014J"]["m_ni_arnett_rule_alt_dust"] = round(v["alt_dust"]["arnett_rule_m_ni"], 3)
    out["SN2014J"]["L_peak_alt_dust_1e43"] = round(v["alt_dust"]["peak"]["L_peak"] / 1e43, 3)
    out["SN2014J"]["m_ni_peak_only_t0_from_2011fe"] = _pm(v["peak_only_t0_from_2011fe"]["m_ni"],
                                                          v["peak_only_t0_from_2011fe"]["m_ni_err"])
    g = res["gamma"]
    de = g["distance_frac_err"]
    out["SN2014J_gamma"] = {
        k: {"m_ni": _pm(g[k]["combined"][0], g[k]["combined"][0] * np.hypot(g[k]["combined"][1] / g[k]["combined"][0], de)),
            "m_ni_847": round(g[k][847][0], 3), "m_ni_1238": round(g[k][1238][0], 3)}
        for k in ("transparent", "uniform", "exponential")}
    for k in ("uniform", "exponential"):
        out["SN2014J_gamma"][k]["escape_75d"] = {str(a): round(b, 3) for a, b in g[k]["escape_75d"].items()}
        out["SN2014J_gamma"][k]["t0_model_kappa025"] = round(g[k]["t0_model_kappa025"], 1)
    for sn, keys in (("SN2023ixf", ("tail_measured_nir", "tail_full_year", "tail_full_trapping_early")),
                     ("SN1987A", ("tail_full_trapping", "tail_with_escape"))):
        f = res[sn]["fit"]
        de = f["distance_frac_err"]
        out[sn] = {}
        for k in keys:
            m = f[k]["m_ni"]
            out[sn][k] = {"m_ni": _pm(m, m * np.hypot(f[k]["m_ni_err"] / m, de), 4),
                          "t0_days": None if f[k]["t0"] is None else _pm(f[k]["t0"], f[k]["t0_err"], 1),
                          "window": f[k]["window"], "rms_dex": round(f[k]["rms_dex"], 3)}
    z = res["ztf"]["good"]
    blue = z[z["low_dust"] & ~z["outlier"]]
    ph = res["ztf"]["phillips"]
    out["ztf"] = {
        "n_measured": int(len(z)), "n_candidates": int(len(res["ztf"]["all"])),
        "n_low_dust": int(len(blue)), "n_outliers": int(z["outlier"].sum()),
        "outliers": z[z["outlier"]][["iau", "M_g", "dm15_g", "m_ni"]].round(3).to_dict("records"),
        "blue_cut_gr": Z.LOW_DUST_GR, "phillips_scatter_mag": round(ph["sigma"], 3),
        "low_dust_types": blue["type"].value_counts().to_dict(),
        "m_ni_median": round(float(blue["m_ni"].median()), 3),
        "m_ni_p16": round(float(blue["m_ni"].quantile(0.16)), 3),
        "m_ni_p84": round(float(blue["m_ni"].quantile(0.84)), 3),
        "rise_bol_median": round(float(blue["t_rise_bol"].median()), 2),
        "rise_g_median_all": round(float(z["t_rise_g"].median()), 2),
        "phillips_slope_mag_per_mag": round(ph["slope"], 2),
        "corr_dm15_logmni": round(float(np.corrcoef(blue["dm15_g"], np.log10(blue["m_ni"]))[0, 1]), 2),
        "m_ni_median_by_cut": {str(c): round(float(z[(z["gr_peak"] <= c) & ~z["outlier"]]["m_ni"].median()), 3)
                               for c in (-0.10, -0.15, -0.20)},
        "calibration_2011fe": {k: round(float(res["ztf"]["calibration"][k]), 3)
                               for k in ("t_gpeak", "g_peak", "gr_color", "M_g", "bol_minus_g_peak_days")},
    }
    out["literature"] = LIT.NICKEL
    return out

