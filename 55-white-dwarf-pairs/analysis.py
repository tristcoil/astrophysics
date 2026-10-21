"""Every number in the README: ZTF J1539+5027's ephemeris, the LISA binaries, their
fates, and the Galactic foreground. ``run_all()`` returns (results, artefacts); the
results go to data/results.json, the artefacts feed figures.py.

Nothing published about ZTF J1539's *timing* enters our fits: the period is found
from scratch in the 2018 ZTF data, and Burdge et al.'s / Shaw et al.'s ephemerides
are only compared with at the end.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle

import config as C
import fates
import gw
import kupfer
import population as pop
import ptf_phot
import timing as T
import ztf_phot

BANDS = ("zg", "zr", "zi")
N_HARM = {"zg": 4, "zr": 4, "zi": 3, "ptf": 2}
#: Shaw et al. (2024, AJ), "Calibrating the clock of JWST": K. Burdge's updated ephemeris
#: from ~6 years of HiPERCAM/CHIMERA monitoring (first published there).
SHAW24 = {"T0_bjd_tdb": 59341.97656626 + 2400000.5, "T0_err_s": 0.022, "P_d": 0.0048008034992,
          "P_err_d": 1.3e-12, "Pdot": -2.36812e-11, "Pdot_err": 0.00044e-11}
SEASON_MIN_POINTS = 25


def _log(msg):
    print(msg, flush=True)


# ---------------------------------------------------------------------------
# Photometry
# ---------------------------------------------------------------------------


def load_ztf() -> pd.DataFrame:
    df = ztf_phot.run()
    df = df.copy()
    df["band"] = df["filtercode"]
    df["year"] = pd.to_datetime(df["shutopen_utc"]).dt.year
    df["flux"] = df["flux_ujy"]
    df["err"] = df["flux_err_ujy"]
    df["group"] = df["band"] + "_" + df["year"].astype(str)
    # cuts: a PSF footprint free of masked pixels, a measurable noise level, the
    # pipeline's own "usable" status, and no failed subtraction. The star's whole flux
    # is ~25-35 microjansky, so a difference flux beyond +-100 cannot be the star
    # (three exposures, -2,800 to -6,700: subtraction artefacts).
    ok = ((df["n_bad_in_psf"] == 0) & np.isfinite(df["flux"]) & np.isfinite(df["err"]) & (df["err"] > 0)
          & (df["status"] == 1) & (df["seeing_px"] < 4.0) & (df["flux"].abs() < 100.0)
          & (df["err"] < 4 * df.groupby("band")["err"].transform("median")))
    df["cut_basic"] = ~ok
    return df


def load_ptf() -> pd.DataFrame:
    df = ptf_phot.run().copy()
    df["band"] = "ptf"
    df["year"] = pd.to_datetime(df["utc_open"]).dt.year
    df["flux"] = df["flux_ujy"]
    df["err"] = df["flux_err_ujy"]
    df["group"] = "ptf_" + df["year"].astype(str)
    ok = np.isfinite(df["flux"]) & np.isfinite(df["err"]) & (df["n_zp"] >= 8) & (df["err"] < 4 * df["err"].median())
    df["cut_basic"] = ~ok
    return df


def bands_from(df, n_harm=N_HARM, sines=False, eclipses=True, groups="group"):
    out = []
    for b in [b for b in BANDS + ("ptf",) if b in set(df["band"])]:
        d = df[df["band"] == b]
        out.append(T.make_band(d, b, group_col=groups, n_harm=n_harm[b], sines=sines, eclipses=eclipses))
    return out


# ---------------------------------------------------------------------------
# Finding the period and its drift without outside help
# ---------------------------------------------------------------------------


def discover_period(df) -> dict:
    """Lomb-Scargle on the 2018 season (the ZTF high-cadence survey that found it)."""
    d = df[(df["year"] == 2018) & (df["band"].isin(["zg", "zr"]))]
    y = d["flux"] - d.groupby("band")["flux"].transform("median")
    freq = np.arange(20.0, 480.0, 0.0002)  # cycles/day: periods 3 min - 72 min
    ls = LombScargle(d["bjd_tdb"], y, d["err"])
    power = ls.power(freq, method="fast", assume_regular_frequency=True)
    i = int(np.argmax(power))
    # refine on a fine grid
    fine = np.arange(freq[i] - 0.002, freq[i] + 0.002, 1e-6)
    pf = ls.power(fine)
    j = int(np.argmax(pf))
    fap = float(ls.false_alarm_probability(power[i], minimum_frequency=freq[0], maximum_frequency=freq[-1]))
    order = np.argsort(power)[::-1]
    # strongest peak that is not within 0.02 c/d of the best: the 1-day alias
    alias = next(int(k) for k in order if abs(freq[k] - freq[i]) > 0.02)
    return {"f_best": float(fine[j]), "P_best_s": float(86400 / fine[j]), "power": float(pf[j]),
            "fap": fap, "n_points": int(len(d)), "t_mean": float(d["bjd_tdb"].mean()),
            "alias_f": float(freq[alias]), "alias_power": float(power[alias]),
            "freq": freq, "pgram": power}


def season_frequencies(df, f_guess) -> pd.DataFrame:
    """Constant-period fit per season: local orbital frequency and its error."""
    rows = []
    for yr, s in df.groupby("year"):
        if len(s) < 60:
            continue
        bands = bands_from(s, n_harm={"zg": 2, "zr": 2, "zi": 2, "ptf": 2}, sines=True, eclipses=False)
        t_mid = float(s["bjd_tdb"].mean())
        x_mid = t_mid - T.T_REF
        # coarse scan over frequency with a harmonic (linear) model: f fixed, phase absorbed
        grid = np.arange(f_guess - 0.004, f_guess + 0.004, 0.00001)
        chi = []
        for fq in grid:
            th = np.array([-fq * x_mid, fq, 0.0])
            chi.append(T.chi2(th, bands))
        chi = np.array(chi)
        k = int(np.argmin(chi))
        lo, hi = max(k - 5, 0), min(k + 6, grid.size)
        a, b, c = np.polyfit(grid[lo:hi] - grid[k], chi[lo:hi], 2)
        fbest = grid[k] - b / (2 * a)
        n = sum(bb.t.size for bb in bands)
        red = max(chi.min() / max(n - 10, 1), 1.0)
        rows.append({"year": int(yr), "t_mid": t_mid, "f_local": float(fbest), "f_err": float(np.sqrt(red / a)),
                     "n": int(n), "P_local_s": 86400 / float(fbest),
                     "P_err_s": 86400 / float(fbest) ** 2 * float(np.sqrt(red / a))})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Global ephemeris fits
# ---------------------------------------------------------------------------


def initial_theta(seasons: pd.DataFrame, bands) -> np.ndarray:
    c = np.polyfit(seasons["t_mid"] - T.T_REF, seasons["f_local"], 1, w=1 / seasons["f_err"])
    best = None
    for phi0 in np.linspace(0, 1, 50, endpoint=False):
        th = np.array([phi0, c[1], c[0]])
        v = T.chi2(th, bands)
        if best is None or v < best[0]:
            best = (v, th)
    return best[1]


def robust_fit(df, theta0, n_harm=N_HARM, sines=False, clip=5.0, max_rounds=4, extra=None):
    """Fit, reject >clip-sigma outliers (errors first rescaled per band), refit."""
    d = df.copy()
    th = np.asarray(theta0, float)
    for _ in range(max_rounds):
        bands = bands_from(d, n_harm=n_harm, sines=sines)
        th, cov, c2, dof = T.fit_ephemeris(bands, th if extra is None else np.r_[th[:3], extra])
        if extra is not None:
            extra = th[3:]
        # per-band error scale from the residuals, then clip
        keep_all = []
        for b in bands:
            _, _, r = T.solve_band(b, th)
            s = 1.4826 * np.median(np.abs(r))
            keep_all.append(np.abs(r) < clip * max(s, 1e-9))
        keep = np.concatenate(keep_all)
        idx = np.concatenate([d.index[d["band"] == b.name].to_numpy() for b in bands])
        if keep.all():
            break
        d = d.drop(index=idx[~keep])
    # phase zero must be the *deeper* eclipse (the hot star hidden). The two eclipse terms
    # have free amplitudes, so a fit can lock on half an orbit away, with the deep eclipse
    # absorbed by the shallow-eclipse term: if so, shift by half a cycle and refit.
    bands = bands_from(d, n_harm=n_harm, sines=sines)
    main = max(bands, key=lambda b: b.t.size)
    beta, _, _ = T.solve_band(main, th)
    names = T.shape_names(main.n_harm, main.sines, main.eclipses)
    coef = dict(zip(names, beta[main.G.shape[1]:]))
    flipped = coef["d1_primary"] < coef["d2_secondary"]
    if flipped:
        th = np.array(th, float)
        th[0] += 0.5
        th, cov, c2, dof = T.fit_ephemeris(bands, th)
    # final error rescaling to reduced chi^2 = 1 per band
    scales = {}
    for b in bands:
        _, _, r = T.solve_band(b, th)
        n_lin = b.G.shape[1] + len(T.shape_names(b.n_harm, b.sines, b.eclipses))
        scales[b.name] = float(np.sqrt(np.sum(r**2) / max(r.size - n_lin - 1, 1)))
    d = d.copy()
    d["err"] = d["err"] * d["band"].map(scales)
    bands = bands_from(d, n_harm=n_harm, sines=sines)
    th, cov, c2, dof = T.fit_ephemeris(bands, th)
    return {"theta": th, "cov": cov, "chi2": c2, "dof": dof, "data": d, "error_scales": scales,
            "n_rejected": int(len(df) - len(d)), "bands": bands, "half_cycle_flip_corrected": bool(flipped)}


def shape_coefficients(fit) -> dict:
    """Best-fit shape (non-offset) coefficients per band at the fitted theta."""
    out = {}
    for b in fit["bands"]:
        beta, _, _ = T.solve_band(b, fit["theta"])
        out[b.name] = {"coef": beta[b.G.shape[1]:], "offsets": dict(zip(b.group_names, beta[: b.G.shape[1]])),
                       "names": T.shape_names(b.n_harm, b.sines, b.eclipses)}
    return out


def pdot_landscape(df_fit, theta, pdots):
    """Profile chi^2 against P-dot: f0 and phi0 re-fitted at every trial P-dot."""
    bands = bands_from(df_fit)
    out = []
    th = np.array(theta, float)
    f0 = th[1]
    for pd_ in pdots:
        f1 = -pd_ * f0**2

        def res(p, f1=f1):
            return T.residuals(np.array([p[0], p[1], f1]), bands)

        best = None
        # two phase starts: a wrong P-dot can land half an orbit away
        for dphi in (0.0, 0.5):
            _, _, c2, _ = T.fit_ephemeris(None, np.array([th[0] + dphi, f0]), x_scale=np.array([1e-3, 1e-6]),
                                          residual_fn=res, n_lin=0)
            if best is None or c2 < best:
                best = c2
        out.append(best)
    return np.array(out)


def season_phase_offsets(df_fit, theta, shapes, ptf=None, shape_for_ptf="zr"):
    """Per season and band: phase shift of the fixed-shape template -> eclipse-time offset."""
    rows = []
    sets = [(df_fit, None)] + ([(ptf, shape_for_ptf)] if ptf is not None else [])
    for frame, borrowed in sets:
        for (band, yr), s in frame.groupby(["band", "year"]):
            if len(s) < SEASON_MIN_POINTS:
                continue
            name = borrowed or band
            sh = shapes[name]
            b = T.make_band(s, band, n_harm=N_HARM.get(name, 4))
            b.n_harm = N_HARM[name]
            best, err, shifts, vals = T.phase_offset_fit([b], theta, [sh["coef"]])
            t_mid = float(s["bjd_tdb"].mean())
            x = t_mid - T.T_REF
            f_here = theta[1] + theta[2] * x
            rows.append({"band": band, "year": int(yr), "t_mid": t_mid, "n": int(len(s)), "dphi": float(best),
                         "dphi_err": float(err), "dt_s": float(-best / f_here * 86400),
                         "dt_err_s": float(err / f_here * 86400)})
    return pd.DataFrame(rows)


def ephemeris_summary(theta, cov):
    e = T.theta_to_ephemeris(theta, cov)
    return {k: float(v) for k, v in e.items()}


def compare_literature(theta, cov) -> dict:
    """Our ephemeris at Burdge et al.'s and Shaw et al.'s reference epochs."""
    out = {}
    for name, T0, P_s, P_err_s, Pdot, Pdot_err, T0_err in (
            ("burdge2019", C.BURDGE19["T0_bjd_tdb"][0], C.BURDGE19["P_s"][0], C.BURDGE19["P_s"][1],
             C.BURDGE19["Pdot"][0], C.BURDGE19["Pdot"][1], C.BURDGE19["T0_bjd_tdb"][1] * 86400),
            ("shaw2024", SHAW24["T0_bjd_tdb"], SHAW24["P_d"] * 86400, SHAW24["P_err_d"] * 86400,
             SHAW24["Pdot"], SHAW24["Pdot_err"], SHAW24["T0_err_s"])):
        t_pred, n = T.predicted_eclipse_time(theta, T0)
        x = T0 - T.T_REF
        f_here = theta[1] + theta[2] * x
        # propagate our covariance to the eclipse time and period at that epoch
        g_t = np.array([-1 / f_here, -x / f_here, -0.5 * x**2 / f_here] + [0] * (len(theta) - 3)) * 86400
        g_P = np.array([0, -86400 / f_here**2, -86400 * x / f_here**2] + [0] * (len(theta) - 3))
        e = T.theta_to_ephemeris(theta, cov)
        out[name] = {"dT0_s": float((t_pred - T0) * 86400), "dT0_err_s": float(np.sqrt(g_t @ cov @ g_t)),
                     "their_T0_err_s": float(T0_err), "P_ours_s": float(86400 / f_here),
                     "P_ours_err_s": float(np.sqrt(g_P @ cov @ g_P)), "P_theirs_s": float(P_s),
                     "P_theirs_err_s": float(P_err_s), "Pdot_ours": e["Pdot"], "Pdot_ours_err": e["Pdot_err"],
                     "Pdot_theirs": float(Pdot), "Pdot_theirs_err": float(Pdot_err)}
        o = out[name]
        o["dT0_sigma"] = o["dT0_s"] / np.hypot(o["dT0_err_s"], o["their_T0_err_s"])
        o["dP_sigma"] = (o["P_ours_s"] - o["P_theirs_s"]) / np.hypot(o["P_ours_err_s"], o["P_theirs_err_s"])
        o["dPdot_sigma"] = (o["Pdot_ours"] - o["Pdot_theirs"]) / np.hypot(o["Pdot_ours_err"], o["Pdot_theirs_err"])
    return out


def derived_quantities(results, shapes) -> dict:
    """Numbers quoted in the README, derived from our ephemeris and Burdge et al.'s system."""
    e = results["joint_ztf_ptf"]
    P, pdot = e["P_s"], e["Pdot"]
    m1, m2 = C.BURDGE19["M1"][0], C.BURDGE19["M2"][0]
    a = C.BURDGE19["a_rsun"][0] * C.R_SUN
    out = {}
    # the total primary eclipse removes the hot star's light: its depth is the hot star's flux
    for band, sh in shapes.items():
        d1 = float(dict(zip(sh["names"], sh["coef"]))["d1_primary"])
        out[f"hot_star_flux_ujy_{band}"] = d1
        out[f"hot_star_mag_{band}"] = float(23.9 - 2.5 * np.log10(d1)) if d1 > 0 else None
    # separation shrinks as a ~ P^(2/3): da/dt = (2/3) a Pdot / P
    dadt = 2 / 3 * a * pdot / P
    out["da_dt_m_per_yr"] = float(dadt * C.YEAR)
    out["orbit_shrink_us_per_yr"] = float(pdot * C.YEAR * 1e6)
    # power radiated: -dE/dt with E = -G m1 m2 / (2a)  ->  G m1 m2 Pdot / (3 a P)
    gm1m2 = C.G * m1 * m2 * C.M_SUN**2
    l_gw = float(abs(gm1m2 * pdot / (3 * a * P)))
    out["L_gw_W"] = l_gw
    out["L_gw_Lsun"] = l_gw / 3.828e26
    # the hot white dwarf's own luminosity (Burdge et al.: R1, Teff1)
    r1 = C.BURDGE19["R1_rsun"][0]
    out["L_hot_star_Lsun"] = float(r1**2 * (48900.0 / 5772.0) ** 4)
    out["orbital_speed_km_s"] = float(2 * np.pi * a / P / 1e3)
    f = 2 / P
    out["A_strain_at_2p34kpc"] = float(gw.amplitude(f, gw.chirp_mass(m1, m2), C.BURDGE19["d_kpc"][0] * C.KPC))
    out["t_coalesce_yr_from_pdot"] = float(3 / 8 * P / abs(pdot) / C.YEAR)
    return out


def chirp_mass_analysis(pdot, pdot_err, P_s, rng) -> dict:
    """Chirp mass if the decay is pure GR; and P-dot predicted from Burdge's masses."""
    mc_gr = float(gw.chirp_mass_from_pdot(P_s, pdot))
    mc_err = float(0.6 * mc_gr * pdot_err / abs(pdot))
    tide = C.BURDGE19["tide_fraction"]
    mc_tide = float(gw.chirp_mass_from_pdot(P_s, pdot / (1 + tide)))
    # masses with asymmetric errors -> split normal draws
    def split_normal(v, ep, em, n):
        z = rng.standard_normal(n)
        return np.where(z > 0, v + z * ep, v + z * em)

    n = 200_000
    m1 = split_normal(*C.BURDGE19["M1"], n)
    m2 = split_normal(*C.BURDGE19["M2"], n)
    mc_masses = gw.chirp_mass(m1, m2)
    pdot_gr = gw.pdot_gw(P_s, mc_masses)
    ratio = pdot / pdot_gr
    return {"Mc_gr_only": mc_gr, "Mc_gr_only_err": mc_err, "Mc_if_6p7pct_tides": mc_tide,
            "Mc_from_burdge_masses": float(np.median(mc_masses)),
            "Mc_from_burdge_masses_err": float(np.std(mc_masses)),
            "Pdot_gr_from_masses": float(np.median(pdot_gr)), "Pdot_gr_from_masses_err": float(np.std(pdot_gr)),
            "ratio_obs_over_gr": float(np.median(ratio)), "ratio_err": float(np.std(ratio)),
            "tide_prediction": tide,
            "shklovskii_pdot": float(gw.shklovskii_pdot(P_s, np.hypot(C.TARGET_PMRA, C.TARGET_PMDEC),
                                                         C.BURDGE19["d_kpc"][0] * C.KPC))}


# ---------------------------------------------------------------------------
# The LISA binaries
# ---------------------------------------------------------------------------


def lisa_binaries(rng, j1539_mc=None, n_mc=4000) -> pd.DataFrame:
    df = kupfer.load()
    rows = []
    for r in df.to_dict("records"):
        f = 2.0 / r["P_s"]
        m1, m2 = r["m1"], r["m2"]
        mc = float(gw.chirp_mass(m1, m2))
        d = r["d_pc"] * C.PC
        incl = np.radians(r["incl"]) if r["incl_kind"] not in ("assumed",) else None
        A = float(gw.amplitude(f, mc, d))
        out = dict(r)
        out.update({"f_Hz": f, "Mc": mc, "A": A, "fdot_gw": float(gw.fdot_gw(f, mc)),
                    "t_coalesce_yr": float(gw.time_to_coalescence(r["P_s"], m1, m2) / C.YEAR)})
        for label, t_yr in (("3mo", 0.25), ("1yr", 1.0), ("4yr", 4.0)):
            out[f"snr_{label}"] = float(gw.snr(f, A, t_yr * C.YEAR, incl=incl, t_obs_yr=t_yr))
            out[f"snr_{label}_noconf"] = float(np.sqrt(4 * t_yr * C.YEAR * A**2 *
                                                      (0.8 if incl is None else gw.inclination_factor(incl))
                                                      / gw.sn_instrument(f)))
        # Monte Carlo over masses, distance and inclination uncertainties
        def draw(v, ep, em, kind, lo=None):
            if kind in ("range",):
                return rng.uniform(v - ep, v + ep, n_mc)
            if not np.isfinite(ep):
                ep = em = 0.0 if kind != "lower_limit" else 0.0
            z = rng.standard_normal(n_mc)
            x = np.where(z > 0, v + z * ep, v + z * em)
            return np.clip(x, lo if lo is not None else 1e-3, None)

        mm1 = draw(m1, r["m1_errp"], r["m1_errm"], r["m1_kind"], 0.05)
        mm2 = draw(m2, r["m2_errp"], r["m2_errm"], r["m2_kind"], 0.01)
        dd = draw(r["d_pc"], r["d_pc_errp"], r["d_pc_errm"], r["d_pc_kind"], 20.0) * C.PC
        if incl is None:
            ii = np.arccos(rng.uniform(-1, 1, n_mc))
        else:
            ii = np.radians(draw(r["incl"], r["incl_errp"] if np.isfinite(r["incl_errp"]) else 5.0,
                                 r["incl_errm"] if np.isfinite(r["incl_errm"]) else 5.0, r["incl_kind"], 0.0))
        AA = gw.amplitude(f, gw.chirp_mass(mm1, mm2), dd)
        s4 = gw.snr(f, AA, 4 * C.YEAR, incl=ii, t_obs_yr=4.0)
        s3 = gw.snr(f, AA, 0.25 * C.YEAR, incl=ii, t_obs_yr=0.25)
        out.update({"A_p16": float(np.percentile(AA, 16)), "A_p84": float(np.percentile(AA, 84)),
                    "snr_4yr_p16": float(np.percentile(s4, 16)), "snr_4yr_p84": float(np.percentile(s4, 84)),
                    "snr_3mo_p16": float(np.percentile(s3, 16)), "snr_3mo_p84": float(np.percentile(s3, 84)),
                    "p_snr3mo_gt7": float(np.mean(s3 > 7))})
        rows.append(out)
    res = pd.DataFrame(rows)
    if j1539_mc is not None:  # our own chirp mass from the orbital decay, for comparison
        i = res.index[res["name"] == "ZTFJ1539"][0]
        f = res.at[i, "f_Hz"]
        for dist_label, d_kpc in (("burdge_spectroscopic", C.BURDGE19["d_kpc"][0]),
                                  ("kupfer_parallax_prior", res.at[i, "d_pc"] / 1e3)):
            A = float(gw.amplitude(f, j1539_mc, d_kpc * C.KPC))
            res.at[i, f"ours_A_{dist_label}"] = A
            res.at[i, f"ours_snr4yr_{dist_label}"] = float(gw.snr(f, A, 4 * C.YEAR, incl=np.radians(84.15), t_obs_yr=4.0))
    return res


def fates_table(binaries: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for r in binaries.to_dict("records"):
        if r["type"] not in ("DWD", "sdB"):
            continue
        m_acc, m_don = max(r["m1"], r["m2"]), min(r["m1"], r["m2"])
        st = fates.stability(m_acc, m_don)
        tc = float(fates.time_to_contact(r["P_s"], m_acc, m_don) / C.YEAR)
        p_c = float(fates.contact_period(m_acc, m_don))
        mt = m_acc + m_don
        if st["label"] == "stable":
            fate = "AM CVn (stable mass transfer)"
        elif st["label"] == "unstable":
            fate = ("merger: possible Type Ia supernova" if mt > 1.2 else "merger into one star")
        else:
            fate = "AM CVn or merger (depends on tides)"
        if r["type"] == "sdB":
            fate = "hot subdwarf still burning helium: becomes a WD pair first"
        rows.append({"name": r["name"], "type": r["type"], "status": r["status"], "P_min": r["P_s"] / 60,
                     "m_acc": m_acc, "m_don": m_don, "q": st["q"], "M_total": mt, "stability": st["label"],
                     "mode": st["mode"], "q_stable_below": st["q_stable_below"],
                     "q_unstable_above": st["q_unstable_above"], "P_contact_min": p_c / 60,
                     "t_contact_yr": tc, "fate": fate, "masses_assumed": r["m1_kind"] == "assumed"})
    return pd.DataFrame(rows).sort_values("t_contact_yr").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


#: Kupfer et al. (2018) Tables 1-2: period (s), masses, inclination (deg; 60 = assumed),
#: distance (pc), and their amplitude (1e-23) and 4-year SNR from a full LISA simulation.
KUPFER18 = {
    "HM Cnc": (321.529, 0.55, 0.27, 38, 5000, 6.4, 211.1),
    "V407 Vul": (569.395, 0.8, 0.177, 60, 1786, 11.0, 169.7),
    "ES Cet": (620.21, 0.8, 0.161, 60, 1584, 10.7, 154.3),
    "AM CVn": (1028.73, 0.68, 0.125, 43, 299, 28.3, 101.2),
    "HP Lib": (1102.70, 0.645, 0.068, 30, 276, 17.5, 43.7),
    "SDSS J0651": (765.5, 0.247, 0.49, 86.9, 933, 16.2, 90.1),
    "SDSS J0935": (1188.0, 0.312, 0.75, 60, 645, 29.9, 44.9),
    "CD-30 11223": (4231.8, 0.54, 0.79, 82.9, 337, 41.5, 4.9),
}


def kupfer2018_check() -> dict:
    """Our amplitude and sky-averaged SNR with Kupfer et al. (2018)'s own inputs."""
    out = {}
    for name, (P, m1, m2, inc, d, a_pub, snr_pub) in KUPFER18.items():
        f = 2 / P
        A = float(gw.amplitude(f, gw.chirp_mass(m1, m2), d * C.PC))
        snr = float(gw.snr(f, A, 4 * C.YEAR, incl=np.radians(inc), t_obs_yr=4.0))
        out[name] = {"A_ours_1e-23": A / 1e-23, "A_kupfer18": a_pub, "snr4yr_ours": snr, "snr4yr_kupfer18": snr_pub,
                     "snr_ratio": snr / snr_pub}
    return out


def post_process(results, art):
    """Cheap comparisons recomputed from cached fits (so config edits show up)."""
    th, cov = art["joint_theta"], art["joint_cov"]
    results["literature_check_joint"] = compare_literature(th, cov)
    results["literature_check_ztf_only"] = compare_literature(art["ztf_theta"], art["ztf_cov"])
    results["kupfer2018_check"] = kupfer2018_check()
    if "shapes" in art and "joint_ztf_ptf" in results:
        results["derived"] = derived_quantities(results, art["shapes"])


def run_all():
    t0 = time.time()
    rng = np.random.default_rng(1539)
    results, art = {}, {}
    # --- photometry ---------------------------------------------------------
    ztf = load_ztf()
    ptf = load_ptf()
    results["ztf_photometry"] = {"n_exposures_with_cutouts": int(len(ztf)),
                                 "n_after_basic_cuts": int((~ztf["cut_basic"]).sum()),
                                 "per_band": ztf[~ztf["cut_basic"]]["band"].value_counts().to_dict(),
                                 "years": sorted(int(y) for y in ztf["year"].unique()),
                                 "median_err_ujy": ztf[~ztf["cut_basic"]].groupby("band")["err"].median().to_dict()}
    results["ptf_photometry"] = {"n_images": int(len(ptf)), "n_after_cuts": int((~ptf["cut_basic"]).sum()),
                                 "per_year": ptf[~ptf["cut_basic"]]["year"].value_counts().sort_index().to_dict(),
                                 "noise_factor": float(ptf["noise_factor"].iloc[0])}
    zgood = ztf[~ztf["cut_basic"]].copy()
    pgood = ptf[~ptf["cut_basic"]].copy()
    _log(f"photometry: ZTF {len(zgood)} points, PTF {len(pgood)} points ({time.time() - t0:.0f} s)")
    # --- discovery and season frequencies ----------------------------------
    disc = discover_period(zgood)
    art["periodogram"] = (disc.pop("freq"), disc.pop("pgram"))
    results["discovery"] = disc
    seasons = season_frequencies(zgood[zgood["band"].isin(["zg", "zr"])], disc["f_best"])
    results["season_frequencies"] = seasons.to_dict("records")
    art["seasons"] = seasons
    c = np.polyfit(seasons["t_mid"] - T.T_REF, seasons["f_local"], 1, w=1 / seasons["f_err"])
    results["frequency_drift_only"] = {"f1_per_day2": float(c[0]), "Pdot": float(-c[0] / c[1] ** 2)}
    _log(f"discovery P = {disc['P_best_s']:.4f} s; seasons {len(seasons)}; drift P-dot {-c[0] / c[1] ** 2:.3e}")
    # --- ZTF-only global fit -----------------------------------------------------
    bands0 = bands_from(zgood)
    th0 = initial_theta(seasons, bands0)
    fit = robust_fit(zgood, th0)
    art["ztf_fit"] = fit
    art["ztf_theta"], art["ztf_cov"] = fit["theta"], fit["cov"]
    results["ztf_only"] = {**ephemeris_summary(fit["theta"], fit["cov"]), "chi2": fit["chi2"], "dof": fit["dof"],
                           "n_points": int(len(fit["data"])), "n_rejected": fit["n_rejected"],
                           "error_scales": fit["error_scales"]}
    shapes = shape_coefficients(fit)
    art["shapes"] = shapes
    results["light_curve_shape"] = {b: dict(zip(v["names"], [float(x) for x in v["coef"]])) for b, v in shapes.items()}
    _log(f"ZTF-only: P-dot = {results['ztf_only']['Pdot']:.5e} +- {results['ztf_only']['Pdot_err']:.1e}")
    # --- robustness variants -------------------------------------------------------
    variants = {}
    for name, sub in (("g_only", fit["data"][fit["data"]["band"] == "zg"]),
                      ("r_only", fit["data"][fit["data"]["band"] == "zr"]),
                      ("infobits_zero", fit["data"][fit["data"]["infobits_meta"] == 0]),
                      ("no_2018", fit["data"][fit["data"]["year"] != 2018])):
        try:
            v = robust_fit(sub, fit["theta"])
            variants[name] = {**ephemeris_summary(v["theta"], v["cov"]), "n_points": int(len(v["data"]))}
        except Exception as exc:  # pragma: no cover - reported, not fatal
            variants[name] = {"error": str(exc)}
    v = robust_fit(fit["data"], fit["theta"], sines=True)
    variants["with_sine_terms"] = {**ephemeris_summary(v["theta"], v["cov"]), "n_points": int(len(v["data"]))}
    results["ztf_variants"] = variants
    _log(f"variants done ({time.time() - t0:.0f} s)")
    # --- P-dot landscape ------------------------------------------------------------
    pd_grid = np.linspace(-6e-11, 2e-11, 81)
    art["landscape_ztf"] = (pd_grid, pdot_landscape(fit["data"], fit["theta"], pd_grid))
    # --- PTF: per-season phases with the ZTF r template, then a joint fit ------------
    pt = pgood.copy()
    pt["exptime_s"] = pt["exptime_s"].astype(float)
    art["ptf_data"] = pt
    pred = ptf_vs_ztf_prediction(fit, pt, shapes)
    art["ptf_prediction"] = pred
    results["ptf_vs_ztf_prediction"] = pred.drop(columns="scan").to_dict("records")
    # joint fit: PTF fitted with the ZTF r-band shape (amplitude and offset free per season)
    joint = joint_fit(fit, pt, shapes)
    art["joint"] = joint
    art["joint_theta"], art["joint_cov"] = joint["theta"], joint["cov"]
    results["joint_ztf_ptf"] = {**ephemeris_summary(joint["theta"], joint["cov"]), "chi2": joint["chi2"],
                                "dof": joint["dof"], "ptf_error_scale": joint["ptf_error_scale"]}
    results["joint_cubic"] = joint["cubic"]
    _log(f"joint: P-dot = {results['joint_ztf_ptf']['Pdot']:.5e} +- {results['joint_ztf_ptf']['Pdot_err']:.1e}")
    # --- O-C points -----------------------------------------------------------------
    oc = season_phase_offsets(fit["data"], joint["theta"], shapes)
    oc_ptf = []
    for yr, s in pt.groupby("year"):
        b = T.make_band(s, "ptf", n_harm=N_HARM["zr"])
        best, err, _, _ = T.phase_offset_fit([b], joint["theta"], [shapes["zr"]["coef"]])
        t_mid = float(s["bjd_tdb"].mean())
        fh = joint["theta"][1] + joint["theta"][2] * (t_mid - T.T_REF)
        oc_ptf.append({"band": "ptf", "year": int(yr), "t_mid": t_mid, "n": int(len(s)), "dphi": float(best),
                       "dphi_err": float(err), "dt_s": float(-best / fh * 86400), "dt_err_s": float(err / fh * 86400)})
    oc = pd.concat([pd.DataFrame(oc_ptf), oc], ignore_index=True)
    art["oc"] = oc
    results["oc_points"] = oc.to_dict("records")
    # --- comparisons ------------------------------------------------------------------
    post_process(results, art)
    e = results["joint_ztf_ptf"]
    results["chirp_mass"] = chirp_mass_analysis(e["Pdot"], e["Pdot_err"], e["P_s"], rng)
    results["derived"] = derived_quantities(results, shapes)
    _log(f"timing done ({time.time() - t0:.0f} s)")
    # --- LISA -----------------------------------------------------------------------------
    binaries = lisa_binaries(rng, j1539_mc=results["chirp_mass"]["Mc_gr_only"])
    art["binaries"] = binaries
    results["lisa_binaries"] = binaries.to_dict("records")
    ver = binaries["status"] == "verification"
    # two bounds: Robson et al.'s Galactic foreground for that observing time (pessimistic
    # early on: their shortest fit is for 6 months), and instrument noise alone (optimistic)
    s3, s3n = binaries["snr_3mo"] > 7, binaries["snr_3mo_noconf"] > 7
    s4, s4n = binaries["snr_4yr"] > 7, binaries["snr_4yr_noconf"] > 7
    results["lisa_summary"] = {
        "n_binaries": int(len(binaries)), "n_verification_kupfer": int(ver.sum()),
        "n_snr3mo_gt7": int(s3.sum()), "n_verification_with_snr3mo_gt7": int((s3 & ver).sum()),
        "n_snr3mo_gt7_instrument_only": int(s3n.sum()),
        "n_verification_with_snr3mo_gt7_instrument_only": int((s3n & ver).sum()),
        "extra_snr3mo_gt7_instrument_only": binaries.loc[s3n & ~ver, "name"].tolist(),
        "n_snr4yr_gt7": int(s4.sum()), "n_snr4yr_gt7_instrument_only": int(s4n.sum()),
        "below_snr7_4yr_instrument_only": binaries.loc[~s4n, "name"].tolist(),
        "j1539_snr4yr": float(binaries.loc[binaries["name"] == "ZTFJ1539", "snr_4yr"].iloc[0]),
    }
    fat = fates_table(binaries)
    art["fates"] = fat
    results["fates"] = fat.to_dict("records")
    # --- the Galaxy -----------------------------------------------------------------------
    loc = pop.local_density()
    results["local_white_dwarfs"] = {k: v for k, v in loc.items() if k != "masses"}
    art["wd_masses_40pc"] = loc["masses"]
    pops = {}
    for label, kw in (("default", {"variant": "mcmillan17", "primary": "40pc"}),
                      ("korol22_check", {"variant": "korol22", "primary": "kepler15"})):
        f, pw, summ, sample = pop.simulate(**kw)
        res = pop.resolve(f, pw)
        summ.update({"n_resolved_snr7_4yr": res["n_resolved"], "foreground_edge_Hz": res["f_conf_edge"],
                     "iterations": res["iterations"]})
        # how many would SNR > 7 *without* any foreground (upper bound)
        summ["n_snr7_instrument_only"] = int(np.sum(4 * 4 * C.YEAR * pw / gw.sn_instrument(f) > 49))
        pops[label] = summ
        art[f"pop_{label}"] = {"f": res["f"], "s_conf": res["s_conf"], "s_inst": res["s_inst"],
                               "s_all_mean": res["s_all_sources_mean"], "counts": res["n_per_bin_1overT"],
                               "sample": sample, "resolved_f": f[res["resolved"]],
                               "resolved_pw": pw[res["resolved"]]}
        del f, pw, res
        _log(f"population {label}: {summ['n_in_band']:.3g} in band, {summ['n_resolved_snr7_4yr']} resolved "
             f"({time.time() - t0:.0f} s)")
    results["population"] = pops
    results["population"]["korol22_published"] = {k: v for k, v in pop.KOROL22.items() if k.endswith("paper")}
    results["runtime_s"] = time.time() - t0
    return results, art


def joint_fit(fit, pt, shapes):
    """ZTF + PTF, with PTF fitted with the ZTF r-band light-curve shape (scale and
    offset free per PTF season), its errors first rescaled to reduced chi^2 = 1."""
    ztf_bands = bands_from(fit["data"])
    ptf_band = T.make_band(pt, "ptf", n_harm=N_HARM["zr"])
    shape_r = shapes["zr"]["coef"]
    _, _, rp = T.solve_band(ptf_band, fit["theta"], shape_r)
    s = float(np.sqrt(np.sum(rp**2) / max(rp.size - ptf_band.G.shape[1] - 2, 1)))
    ptf_band.sig = ptf_band.sig * max(s, 1.0)

    def res(theta):
        r1 = T.residuals(theta, ztf_bands)
        _, _, r2 = T.solve_band(ptf_band, theta, shape_r)
        return np.concatenate([r1, r2])

    n_lin = sum(b.G.shape[1] + len(T.shape_names(b.n_harm, b.sines, b.eclipses)) for b in ztf_bands) \
        + ptf_band.G.shape[1] + 1
    th, cov, c2, dof = T.fit_ephemeris(None, fit["theta"], residual_fn=res, n_lin=n_lin)
    # is dP/dt itself changing? add a cubic phase term (f2 = d^2f/dt^2)
    th3, cov3, c23, dof3 = T.fit_ephemeris(None, np.r_[th, 0.0], residual_fn=res, n_lin=n_lin)
    f0, f1, f2 = th3[1], th3[2], th3[3]
    cubic = {"f2_per_day3": float(f2), "f2_err": float(np.sqrt(cov3[3, 3])),
             "f2_gr_prediction": float(11 / 3 * f1**2 / f0), "braking_index": float(f0 * f2 / f1**2),
             "braking_index_err": float(f0 * np.sqrt(cov3[3, 3]) / f1**2), "delta_chi2": float(c2 - c23)}
    return {"theta": th, "cov": cov, "chi2": c2, "dof": dof, "ptf_error_scale": s, "cubic": cubic}


def ptf_vs_ztf_prediction(fit, pt, shapes) -> pd.DataFrame:
    """Each PTF season's eclipse phase against what the ZTF-only ephemeris predicts.

    The prediction's own uncertainty comes from the ZTF covariance; a miscount by one
    orbit would show up as an offset of 1.0 cycles."""
    th, cov = fit["theta"], fit["cov"]
    rows = []
    for yr, s in pt.groupby("year"):
        b = T.make_band(s, "ptf", n_harm=N_HARM["zr"])
        best, err, shifts, vals = T.phase_offset_fit([b], th, [shapes["zr"]["coef"]], grid_half_width=0.5,
                                                     n_grid=1001)
        x = float(s["bjd_tdb"].mean()) - T.T_REF
        g = np.array([1.0, x, 0.5 * x**2])
        pred_err = float(np.sqrt(g @ cov[:3, :3] @ g))
        f0 = th[1]
        rows.append({"year": int(yr), "n": int(len(s)), "t_mid": x + T.T_REF, "dphi": float(best),
                     "dphi_err": float(err), "prediction_err": pred_err,
                     "pdot_change_for_one_cycle": float(2 / x**2 / f0**2),
                     "scan": (shifts, vals)})
    return pd.DataFrame(rows)
