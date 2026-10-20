"""Every number in the project, computed from the cached data; writes data/results.json.

    python analysis.py           # everything (first run ~25 min, mostly the Higgs-portal line)
    python analysis.py --force   # recompute the cached pieces too

Arrays for the figures go to data/derived/figdata.npz; headline numbers to data/results.json.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # works under python -I too

import numpy as np  # noqa: E402

import config as C  # noqa: E402
import freezeout as fo  # noqa: E402
import thermo  # noqa: E402

FIG = {}          # arrays for figures.py


def _log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def interp_loglog(x, xp, fp):
    return np.exp(np.interp(np.log(x), np.log(xp), np.log(fp)))


# --- 1. the plasma ---------------------------------------------------------------------------------


def section_eos() -> dict:
    names = ["ss2020", "ss2020_hi", "ss2020_lo", "ss2018", "lms2015", "p2", "p2_ideal", "ideal"]
    out = {}
    T = np.logspace(-5, 4, 1200)
    for n in names:
        if not thermo.available(n):
            continue
        e = thermo.load(n)
        out[n] = {"label": e.label, "g_s0": e.g_s0, "table": thermo.describe(e)}
        FIG[f"eos_{n}_grho"] = e.g_rho(T)
        FIG[f"eos_{n}_gs"] = e.g_s(T)
        FIG[f"eos_{n}_gstarhalf"] = e.gstar_half(T)
    FIG["eos_T"] = T
    Tt, t = thermo.load("ss2020").time_of(np.logspace(-6, 5, 3000))
    FIG["time_T"], FIG["time_t"] = Tt, t
    return out


# --- 2. freeze-out curves (the Kolb & Turner plot) ------------------------------------------------


def section_examples(sv_thermal_100: float) -> dict:
    m = 100.0
    svs = [1e-28, 1e-27, 1e-26, sv_thermal_100, 1e-25, 1e-24]
    x = np.geomspace(1.0, 1000.0, 1500)
    rows = []
    for k, sv in enumerate(svs):
        s = fo.solve(fo.Relic(m, a=sv, g=2), x_start=1.0, dense=True)
        lx = np.log(x)
        Y = np.exp(s.extra["dense"](np.clip(lx, np.log(1.0), np.log(fo.X_TAIL)))[0])
        # beyond x = 150 the dense solution is not defined; use the stored tail
        tail = x > fo.X_TAIL
        Y[tail] = np.exp(np.interp(np.log(x[tail]), np.log(s.x), np.log(s.Y)))
        FIG[f"yx_{k}"] = Y
        rows.append({"sigv": sv, "omega_h2": s.omega_h2, "x_f": s.x_f, "x_gamma_eq_H": s.x_gamma_eq_h,
                     "Y0": s.Y0})
    e = thermo.load("ss2020")
    FIG["yx_x"] = x
    FIG["yx_eq"] = np.exp(fo.ln_yeq(x, 2.0, e.g_s(m / x)))
    Y_target = C.OMEGA_C_H2 / fo.omega_h2_from_Y(m, 1.0, e.g_s0)
    return {"mass": m, "rows": rows, "Y_for_omega_0120": Y_target,
            "T_at_x_f_GeV": m / rows[3]["x_f"],
            "t_at_x_f_s": float(np.interp(math.log(m / rows[3]["x_f"]), np.log(FIG["time_T"][::-1]),
                                          FIG["time_t"][::-1]))}


# --- 3. the thermal cross section ---------------------------------------------------------------


def section_thermal() -> dict:
    out = {}
    masses = np.logspace(math.log10(C.MASS_MIN), math.log10(C.MASS_MAX), C.N_MASS)
    FIG["th_m"] = masses
    _log("thermal cross section: Majorana, Saikawa & Shirai EoS")
    main = fo.thermal_cross_section_grid(masses)
    FIG["th_main"] = main
    # Omega h^2 +- 1 sigma (Planck)
    FIG["th_om_lo"] = fo.thermal_cross_section_grid(masses, target=C.OMEGA_C_H2 + C.OMEGA_C_H2_ERR)
    FIG["th_om_hi"] = fo.thermal_cross_section_grid(masses, target=C.OMEGA_C_H2 - C.OMEGA_C_H2_ERR)
    for n in ["ss2020_hi", "ss2020_lo", "lms2015", "ss2018", "p2", "p2_ideal", "ideal",
              "constant_10.75", "constant_106.75"]:
        if not thermo.available(n):
            continue
        _log(f"  EoS {n}")
        # a constant g* stands for the plasma at freeze-out only: today's entropy stays the real one
        kw = {"g_s0": thermo.load("ss2020").g_s0} if n.startswith("constant") else {}
        FIG[f"th_eos_{n}"] = fo.thermal_cross_section_grid(masses, eos=n, **kw)
    for key, d in C.DM_TYPES.items():
        if key == "majorana":
            continue
        FIG[f"th_type_{key}"] = fo.thermal_cross_section_grid(masses, g=d["g"])
    FIG["th_pwave"] = fo.thermal_cross_section_grid(masses, wave="p")
    FIG["th_pwave_dirac"] = fo.thermal_cross_section_grid(masses, wave="p", g=4)
    FIG["th_pwave_real_scalar"] = fo.thermal_cross_section_grid(masses, wave="p", g=1)
    FIG["th_pwave_vector"] = fo.thermal_cross_section_grid(masses, wave="p", g=3)

    # x_f along the main curve
    xf = []
    for m, sv in zip(masses[::5], main[::5]):
        s = fo.solve(fo.Relic(float(m), a=float(sv)))
        xf.append((float(m), s.x_f, s.x_gamma_eq_h))
    FIG["th_xf"] = np.array(xf)

    def at(m, arr=main):
        return float(interp_loglog(m, masses, arr))

    sel = (masses >= 10) & (masses <= 1e4)
    lowsel = masses < 10
    imax = int(np.argmax(np.where(lowsel, main, 0)))
    out["headline"] = {
        "sigv_100GeV": at(100.0), "sigv_1GeV": at(1.0), "sigv_10GeV": at(10.0), "sigv_1TeV": at(1000.0),
        "sigv_10TeV": at(1e4), "range_10GeV_10TeV": [float(main[sel].min()), float(main[sel].max())],
        "low_mass_peak": {"sigv": float(main[imax]), "mass": float(masses[imax])},
        "dirac_physical_100GeV": 2 * at(100.0, FIG["th_type_dirac"]),
        "pwave_b_100GeV": at(100.0, FIG["th_pwave"]),
        "omega_band_100GeV": [at(100.0, FIG["th_om_lo"]), at(100.0, FIG["th_om_hi"])],
        "eos_band_100GeV": [at(100.0, FIG["th_eos_ss2020_hi"]), at(100.0, FIG["th_eos_ss2020_lo"])],
        "x_f_100GeV": float(np.interp(math.log(100.0), np.log(FIG["th_xf"][:, 0]), FIG["th_xf"][:, 1])),
    }
    # QCD feature: largest relative step, and the EoS sensitivity there
    rel = {n: (FIG[f"th_eos_{n}"] / main - 1).tolist() for n in ["lms2015", "ideal", "constant_10.75",
                                                                  "constant_106.75", "p2", "p2_ideal"]
           if f"th_eos_{n}" in FIG}
    out["eos_relative"] = {n: {"max_abs": float(np.max(np.abs(v))),
                               "at_mass": float(masses[int(np.argmax(np.abs(v)))]),
                               "at_100GeV": float(np.interp(math.log(100.0), np.log(masses), v))}
                           for n, v in rel.items()}

    # validation against Saikawa & Shirai 2020's own tables
    val = {}
    for key, d in C.DM_TYPES.items():
        for wave, prefix in (("s", "S"), ("p", "P")):
            tab = np.loadtxt(C.RAW_DIR / "relic" / f"ss2020_{prefix}_{d['ss2020'][2]}.dat")
            ok = (tab[:, 0] >= C.MASS_MIN) & (tab[:, 0] <= C.MASS_MAX)
            ms = tab[ok, 0][::6]
            ours = fo.thermal_cross_section_grid(ms, g=d["g"], wave=wave)
            r = ours / (tab[ok, 1][::6] * 1e-26)
            val[f"{key}_{wave}"] = {"n": int(len(ms)), "median_ratio": float(np.median(r)),
                                    "min_ratio": float(r.min()), "max_ratio": float(r.max())}
            if key == "majorana" and wave == "s":
                FIG["val_ss_m"], FIG["val_ss_ratio"] = ms, r
    out["validation_ss2020"] = val
    tab = np.loadtxt(C.RAW_DIR / "relic" / "ss2020_S_2.dat")
    FIG["ss2020_m"], FIG["ss2020_sv"] = tab[:, 0], tab[:, 1] * 1e-26
    FIG["ss2020_err_om"], FIG["ss2020_err_eos"] = tab[:, 2] * 1e-26, tab[:, 3] * 1e-26
    # EoS error: ours (half the hi-lo spread) against their column
    ok = (tab[:, 0] >= 1.0) & (tab[:, 0] <= 1e4)
    ours_err = 0.5 * np.abs(interp_loglog(tab[ok, 0], masses, FIG["th_eos_ss2020_lo"])
                            - interp_loglog(tab[ok, 0], masses, FIG["th_eos_ss2020_hi"]))
    out["eos_error_vs_ss2020"] = {"median_ratio": float(np.median(ours_err / (tab[ok, 3] * 1e-26)))}

    # Radau spot checks of the fast stepper
    spots = []
    for m in (0.03, 1.0, 5.0, 100.0, 3000.0, 5e4):
        sv = at(m)
        r = fo.solve(fo.Relic(m, a=sv)).omega_h2
        spots.append({"mass": m, "omega_radau": r, "rel_dev": r / C.OMEGA_C_H2 - 1})
    out["radau_spotchecks"] = spots

    # Steigman, Dasgupta & Beacom 2012 in their own set-up
    sdb = np.loadtxt(C.RAW_DIR / "relic" / "sdb12_m-sv-num.txt")
    kw = dict(eos="lms2015_gs_eq_grho", x_end=1000.0, nonrel_yeq=True, g_s0=3.909, target=C.OMEGA_H2_SDB12)
    ours = fo.thermal_cross_section_grid(sdb[:, 0], **kw)
    r = ours / (sdb[:, 1] * 1e-26)
    FIG["sdb12_m"], FIG["sdb12_sv"], FIG["sdb12_ratio"] = sdb[:, 0], sdb[:, 1] * 1e-26, r
    away = (sdb[:, 0] < 2.5) | (sdb[:, 0] > 6.5)
    out["validation_sdb12"] = {"median_ratio": float(np.median(r)), "max_abs_dev_away_from_qcd":
                               float(np.max(np.abs(r[away] - 1))), "max_abs_dev": float(np.max(np.abs(r - 1))),
                               "worst_mass": float(sdb[int(np.argmax(np.abs(r - 1))), 0])}

    # From SDB12's 2.2 to Saikawa & Shirai's 2.1: one change at a time, at 100 GeV
    steps = [("SDB12 set-up (WMAP-7 Omega h^2 = 0.11, Laine-Schroeder-type g with g*s = g*rho, "
              "stop at x = 1000, non-relativistic Y_eq)", kw),
             ("... Planck's Omega h^2 = 0.120", {**kw, "target": C.OMEGA_C_H2}),
             ("... integrate to T = 1 eV", {k: v for k, v in kw.items() if k != "x_end"} | {"target": C.OMEGA_C_H2}),
             ("... separate g*s and g*rho (Laine et al.)", dict(eos="lms2015", nonrel_yeq=True, g_s0=3.909)),
             ("... exact Maxwell-Boltzmann Y_eq", dict(eos="lms2015", g_s0=3.909)),
             ("... Saikawa & Shirai's equation of state (g*s0 = 3.931)", dict(eos="ss2020"))]
    dec = []
    for label, k in steps:
        sv = float(fo.thermal_cross_section_grid([100.0], **k)[0])
        dec.append({"step": label, "sigv_100GeV": sv})
    out["decomposition_100GeV"] = dec
    return out


# --- 4. analytic approximations ------------------------------------------------------------------


def section_analytic() -> dict:
    masses = np.array([0.3, 1.0, 3.0, 10.0, 100.0, 1e3, 1e4])
    svs = np.array([2e-27, 2.2e-26, 2e-25])
    rows = []
    for m in masses:
        for sv in svs:
            num = fo.solve(fo.Relic(float(m), a=float(sv)))
            kt = fo.kolb_turner(float(m), float(sv))
            im = fo.improved_analytic(float(m), float(sv))
            rows.append({"mass": float(m), "sigv": float(sv), "omega_num": num.omega_h2, "x_f_num": num.x_f,
                         "omega_kt": kt["omega_h2"], "x_f_kt": kt["x_f"], "omega_improved": im["omega_h2"],
                         "x_f_improved": im["x_f"]})
    FIG["an_rows"] = np.array([[r["mass"], r["sigv"], r["omega_num"], r["omega_kt"], r["omega_improved"],
                                r["x_f_num"], r["x_f_kt"], r["x_f_improved"]] for r in rows])
    # the thermal cross section from each method along the mass axis
    ms = np.logspace(-1, 4, 26)
    num = fo.thermal_cross_section_grid(ms)
    kt, im = [], []
    from scipy import optimize

    for m in ms:
        kt.append(math.exp(optimize.brentq(lambda l: math.log(fo.kolb_turner(m, math.exp(l))["omega_h2"] / C.OMEGA_C_H2),
                                           math.log(1e-27), math.log(1e-24))))
        im.append(math.exp(optimize.brentq(lambda l: math.log(fo.improved_analytic(m, math.exp(l))["omega_h2"] / C.OMEGA_C_H2),
                                           math.log(1e-27), math.log(1e-24))))
    FIG["an_m"], FIG["an_num"], FIG["an_kt"], FIG["an_im"] = ms, num, np.array(kt), np.array(im)
    r_kt, r_im = np.array(kt) / num, np.array(im) / num
    return {"grid": rows, "thermal_ratio_kt": [float(r_kt.min()), float(r_kt.max())],
            "thermal_ratio_improved": [float(r_im.min()), float(r_im.max())],
            "kt_at_100GeV": float(np.interp(2.0, np.log10(ms), r_kt)),
            "improved_at_100GeV": float(np.interp(2.0, np.log10(ms), r_im))}


# --- 5. neutrons ---------------------------------------------------------------------------------


def section_neutrons() -> dict:
    import neutrons as nu

    out = {}
    eta = nu.eta_from_omega_b()
    out["eta"] = eta
    runs = {}
    _log("neutrons: instant decoupling")
    r_inst = nu.run(eta=eta, background=nu.instant_background())
    _log("neutrons: Escudero background")
    r_esc = nu.run(eta=eta, background=nu.escudero_background())
    _log("neutrons: free neutrons only")
    r_free = nu.free_neutrons_only(background=nu.escudero_background())
    for name, r in (("instant", r_inst), ("escudero", r_esc)):
        runs[name] = {"Y_p": r.Y_p, "D_H": r.D_H, "He3_H": r.He3_H, "np_1s": r.np_ratio_1s,
                      "T_freeze_MeV": r.T_freeze * 1e3, "t_freeze_s": r.t_freeze, "t_nuc_s": r.t_nuc,
                      "T_nuc_keV": r.T_nuc * 1e6, "X_n_at_nuc": r.X_n_at_nuc}
    # sensitivity: neutron lifetime (bottle vs beam), extra radiation
    _log("neutrons: tau_n beam, N_eff variants")
    r_beam = nu.run(eta=eta, tau_n=C.TAU_N_BEAM, background=nu.escudero_background())
    runs["escudero_tau_beam"] = {"Y_p": r_beam.Y_p, "D_H": r_beam.D_H}
    for neff in (2.0, 4.0):
        r = nu.run(eta=eta, background=nu.instant_background(N_eff=neff))
        runs[f"instant_Neff{neff:g}"] = {"Y_p": r.Y_p, "D_H": r.D_H}
    out["runs"] = runs
    out["dYp_dNeff"] = 0.5 * (runs["instant_Neff4"]["Y_p"] - runs["instant_Neff2"]["Y_p"])
    out["dYp_beam_minus_bottle"] = runs["escudero_tau_beam"]["Y_p"] - runs["escudero"]["Y_p"]

    # curves for the figure
    FIG["n_t"], FIG["n_T"] = r_esc.t, r_esc.T_gamma
    for s in nu.SPECIES:
        FIG[f"n_Y_{s}"] = r_esc.Y[s]
    FIG["n_free_np"] = r_free.Y["n"] / r_free.Y["p"]
    FIG["n_eq_np"] = np.exp(-C.Q_NP / r_esc.T_gamma)
    bg = nu.escudero_background()
    lam = np.array([sum(nu.weak_rates(float(Tg), float(bg.Tn(t)))) for Tg, t in zip(r_esc.T_gamma[::10], r_esc.t[::10])])
    FIG["n_rate_t"], FIG["n_rate"] = r_esc.t[::10], lam
    FIG["n_H"] = bg.hubble(r_esc.t[::10])

    # PRyMordial: like-for-like and the corrections we leave out
    _log("PRyMordial runs (cached)")
    pr = {}
    cfgs = {"born_small_instant_a": dict(eta=eta, born=True, small=True, aTid=False),
            "full": dict(eta=eta),
            "corrections_small_instant_a": dict(eta=eta, born=False, small=True, aTid=False),
            "born_full_net": dict(eta=eta, born=True, small=False, aTid=True),
            "full_tau_beam": dict(eta=eta, tau_n=C.TAU_N_BEAM)}
    for k, cfg in cfgs.items():
        res = nu.prymordial(**cfg)
        if res is not None:
            pr[k] = {"Y_p": res["Yp_bbn"], "D_H": res["DoH"], "He3_H": res["He3oH"], "N_eff": res["Neff"]}
    out["prymordial"] = pr
    if pr:
        bgp = nu.prymordial_background(eta)
        r_p = nu.run(eta=eta, background=bgp)
        out["ours_on_prymordial_background"] = {"Y_p": r_p.Y_p, "D_H": r_p.D_H, "He3_H": r_p.He3_H}
    return out


# --- 6. recombination ------------------------------------------------------------------------------


def section_recombination() -> dict:
    import recombination as rc

    r = rc.run()
    FIG["rec_z"], FIG["rec_xe"], FIG["rec_saha"], FIG["rec_gam"] = r.z, r.x_e, r.x_saha, r.gamma_over_H
    FIG["rec_t"], FIG["rec_zc"], FIG["rec_xc"] = r.t_s, r.z_camb, r.x_e_camb
    xc = lambda z: float(np.interp(z, r.z_camb[::-1], r.x_e_camb[::-1]))  # noqa: E731
    xo = lambda z: float(np.interp(z, r.z[::-1], r.x_e[::-1]))  # noqa: E731
    return {"z_freeze": r.z_freeze, "x_e_200": r.x_e_200, "x_e_200_camb": xc(200.0),
            "x_e_1000": xo(1000.0), "x_e_1000_camb": xc(1000.0), "x_e_20": xo(20.0), "x_e_20_camb": xc(20.0),
            "saha_at_200": float(np.interp(200.0, r.z[::-1], r.x_saha[::-1]))}


# --- 7. dwarf galaxies -------------------------------------------------------------------------------


def section_dwarfs(force: bool = False) -> dict:
    import gamma_dwarfs as gd
    import pdf_curves as pc

    out = {"samples": {k: len(v) for k, v in gd.samples().items()}}
    th_m = np.logspace(math.log10(C.GAMMA_MASS_MIN), 4.0, 80)
    th = fo.thermal_cross_section_grid(th_m)
    FIG["dw_th_m"], FIG["dw_th"] = th_m, th

    def crossing(ms, ul):
        ms, ul = np.asarray(ms, float), np.asarray(ul, float)
        good = np.isfinite(ul)
        ms, ul = ms[good], ul[good]
        r = np.log(ul / interp_loglog(ms, th_m, th))
        k = np.flatnonzero((r[:-1] < 0) & (r[1:] >= 0))
        return [float(np.exp(np.interp(0, [r[i], r[i + 1]], [np.log(ms[i]), np.log(ms[i + 1])]))) for i in k]

    res = {}
    for ch in ("bb", "tau"):
        for src in ("pppc", "cosmixs"):
            _log(f"dwarfs: {ch} {src}")
            cube = gd.profile_cube(ch, src, force=force)
            for smp in ("benchmark", "measured"):
                lim = gd.stack_limits(cube, smp)
                key = f"{ch}_{src}_{smp}"
                FIG[f"dw_{key}_m"], FIG[f"dw_{key}_ul"] = np.array(lim["mass"]), np.array(lim["ul"])
                FIG[f"dw_{key}_ts"] = np.array(lim["ts_max"])
                j = int(np.nanargmax(lim["ts_max"]))
                res[key] = {"crossing_GeV": crossing(lim["mass"], lim["ul"]), "ts_max": lim["ts_max"][j],
                            "ts_max_mass": lim["mass"][j],
                            "ul_100GeV": float(interp_loglog(100.0, lim["mass"], lim["ul"]))}
        # McDaniel's own TS profiles, summed
        for smp in ("benchmark", "measured"):
            mcd = gd.mcd_stack_limits(smp, ch)
            key = f"{ch}_mcdsum_{smp}"
            m = np.array(mcd["mass"])
            ul = np.array(mcd["ul"])
            FIG[f"dw_{key}_m"], FIG[f"dw_{key}_ul"] = m, ul
            ok = m >= C.GAMMA_MASS_MIN
            ours = interp_loglog(m[ok], FIG[f"dw_{ch}_pppc_{smp}_m"], FIG[f"dw_{ch}_pppc_{smp}_ul"])
            ratio = ours / ul[ok]
            ratio = ratio[np.isfinite(ratio)]
            res[key] = {"crossing_GeV": crossing(m[ok], ul[ok]), "ours_over_theirs_median": float(np.median(ratio)),
                        "ours_over_theirs_range": [float(ratio.min()), float(ratio.max())]}
        # our reproduction with their conventions (bin centre)
        try:
            if not (C.DERIVED_DIR / f"dwarf_ts_{ch}_pppc_mcd.npz").exists():
                raise FileNotFoundError("bin-centre cube not computed (gamma_dwarfs.profile_cube(ch, 'pppc', 'mcd'))")
            cube = gd.profile_cube(ch, "pppc", "mcd", force=False)
            for smp in ("benchmark", "measured"):
                lim = gd.stack_limits(cube, smp)
                m = np.array(FIG[f"dw_{ch}_mcdsum_{smp}_m"])
                ok = m >= C.GAMMA_MASS_MIN
                theirs = FIG[f"dw_{ch}_mcdsum_{smp}_ul"][ok]
                mine = interp_loglog(m[ok], lim["mass"], lim["ul"])
                rr = mine / theirs
                rr = rr[np.isfinite(rr)]
                res[f"{ch}_mcdconv_{smp}"] = {"ours_over_theirs_median": float(np.median(rr)),
                                              "ours_over_theirs_range": [float(rr.min()), float(rr.max())]}
                FIG[f"dw_{ch}_mcdconv_{smp}_m"], FIG[f"dw_{ch}_mcdconv_{smp}_ul"] = np.array(lim["mass"]), np.array(lim["ul"])
        except Exception as e:  # pragma: no cover - optional check
            res[f"{ch}_mcdconv_error"] = str(e)
        # Ackermann et al. 2015 from their figure
        name = "bb" if ch == "bb" else "tautau"
        ack = pc.ackermann2015(C.LIT_DIR / "1503.02641" / f"p302_6yr_source_composite_limit_data_plot_{name}.pdf")
        a = ack["limit_6yr_pass8"]
        FIG[f"dw_{ch}_ack_m"], FIG[f"dw_{ch}_ack_ul"] = np.array(a["mass"]), np.array(a["sigmav"])
        FIG[f"dw_{ch}_ack_med"] = np.array(ack["median_expected"]["sigmav"])
        res[f"{ch}_ackermann2015"] = {"crossing_GeV": crossing(a["mass"], a["sigmav"]),
                                      "calibration_residual_dex": ack["calibration"]["max_tick_residual"]}
        t = ack["thermal_sdb12"]
        sdb = np.loadtxt(C.RAW_DIR / "relic" / "sdb12_m-sv-num.txt")
        tm, ts = np.array(t["mass"]), np.array(t["sigmav"])
        sel = (tm > sdb[0, 0]) & (tm < sdb[-1, 0])
        res[f"{ch}_ackermann2015"]["thermal_line_vs_sdb12_table"] = (
            ts[sel] / (np.exp(np.interp(np.log(tm[sel]), np.log(sdb[:, 0]), np.log(sdb[:, 1]))) * 1e-26)).tolist()
    out["limits"] = res
    # which dwarfs matter: TS of each at the benchmark limit, bb, 100 GeV
    cube = gd.profile_cube("bb", "pppc")
    tab = gd.dwarf_table()
    j = int(np.argmin(np.abs(cube["masses"] - 100.0)))
    lim = gd.stack_limits(cube, "benchmark")
    l_ul = math.log10(lim["ul"][j])
    k = int(np.argmin(np.abs(cube["lsv"] - l_ul)))
    contrib = []
    for i, n in enumerate(cube["names"]):
        contrib.append({"name": n, "label": tab[n]["label"], "log10J": tab[n]["log10J"], "sigmaJ": tab[n]["sigmaJ"],
                        "method": tab[n]["method"], "minus_ts_at_limit": float(-cube["ts"][i, j, k]),
                        "in_benchmark": n in gd.samples()["benchmark"], "in_measured": n in gd.samples()["measured"],
                        "ts_max_100GeV": float(cube["ts"][i, j].max())})
    out["dwarf_contributions_bb_100GeV"] = sorted(contrib, key=lambda r: -r["minus_ts_at_limit"])
    out["gce_calore2015"] = {"bb": {"mass": 48.7, "mass_err": [5.2, 6.4], "sigv": 1.75e-26, "sigv_err": [0.26e-26, 0.28e-26]},
                             "tau": {"mass": 9.96, "mass_err": [0.91, 1.05], "sigv": 0.337e-26,
                                     "sigv_err": [0.048e-26, 0.047e-26]},
                             "source": "Calore, Cholis, McCabe & Weniger 2015 (PRD 91, 063003), Table 2 (model F, "
                                       "rho_local = 0.4 GeV/cm^3, NFW gamma = 1.2)"}
    out["gce_daylan2016"] = {"bb_mass_range": [36, 51], "sigv_range": [1e-26, 3e-26],
                             "source": "Daylan et al. 2016 (Phys. Dark Univ. 12, 1), abstract; rho_local = 0.4"}
    return out


# --- 8. direct detection ----------------------------------------------------------------------------


def section_direct(force: bool = False) -> dict:
    import direct_detection as dd

    lz, xn = dd.lz2024(), dd.xenonnt2025()
    FIG["dd_lz_m"], FIG["dd_lz"] = lz["limit"]
    FIG["dd_lz_med"] = lz["median"][1]
    FIG["dd_xn_m"], FIG["dd_xn"] = xn["limit"]
    for k in dd.OHARE:
        FIG[f"dd_hist_{k}_m"], FIG[f"dd_hist_{k}"] = dd.ohare(k)
    FIG["dd_floor_m"], FIG["dd_floor"] = dd.neutrino_floor()
    for i, reg in enumerate(dd.dama_regions()):
        FIG[f"dd_dama{i}"] = reg
    _log("Higgs portal relic line (cached after the first run, ~10 min)")
    hp = dd.higgs_portal_line(force=force)
    m = np.array(hp["mass"])
    s = np.array(hp["sigma_si"])
    lam = np.array(hp["lambda"])
    FIG["hp_m"], FIG["hp_sigma"], FIG["hp_lambda"] = m, s, lam
    inv_ok = lam <= np.array(hp["lambda_inv_max"])
    win_lz = dd.allowed_windows(m, s, lz["limit"], extra_ok=inv_ok)
    win_xn = dd.allowed_windows(m, s, dd.xenonnt2025()["limit"], extra_ok=inv_ok)
    imin = int(np.argmin(lz["limit"][1]))
    big = m > 1000
    return {"lz2024_min": {"sigma": float(lz["limit"][1][imin]), "mass": float(lz["limit"][0][imin])},
            "lz2024_at_40GeV": float(dd.limit_at(40.0, lz["limit"])[0]),
            "xenonnt2025_min": {"sigma": float(np.min(xn["limit"][1])),
                                "mass": float(xn["limit"][0][int(np.argmin(xn["limit"][1]))])},
            "higgs_portal": {"allowed_by_lz2024": win_lz, "allowed_by_xenonnt2025": win_xn,
                             "sigma_si_high_mass_mean": float(np.mean(s[big])) if big.any() else None,
                             "lambda_1TeV": float(np.interp(math.log(1000.0), np.log(m), lam)),
                             "lambda_10TeV": float(np.interp(math.log(1e4), np.log(m), lam)),
                             "f_N": hp["f_N"], "gamma_h_MeV": 1e3 * hp["gamma_h_GeV"],
                             "lz_over_prediction_1TeV": float(dd.limit_at(1000.0, lz["limit"])[0]
                                                              / np.interp(math.log(1000.0), np.log(m), s))}}


# --- run ---------------------------------------------------------------------------------------------


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return o


def run_core(force: bool = False) -> dict:
    t0 = time.time()
    C.DERIVED_DIR.mkdir(parents=True, exist_ok=True)
    res = {"meta": {"run": dt.datetime.now().isoformat(timespec="seconds"), "omega_c_h2": C.OMEGA_C_H2,
                    "numpy": np.__version__}}
    _log("equation of state")
    res["eos"] = section_eos()
    res["thermal"] = section_thermal()
    _log("freeze-out examples")
    res["examples"] = section_examples(res["thermal"]["headline"]["sigv_100GeV"])
    _log("analytic approximations")
    res["analytic"] = section_analytic()
    res["neutrons"] = section_neutrons()
    _log("recombination")
    res["recombination"] = section_recombination()
    res["dwarfs"] = section_dwarfs(force=force)
    res["direct"] = section_direct(force=force)
    res["meta"]["seconds"] = time.time() - t0
    np.savez_compressed(C.DERIVED_DIR / "figdata.npz", **{k: np.asarray(v) for k, v in FIG.items()})
    C.RESULTS_PATH.write_text(json.dumps(_jsonable(res), indent=1))
    _log(f"wrote {C.RESULTS_PATH.relative_to(C.HERE)} ({time.time() - t0:.0f} s)")
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--force", action="store_true")
    run_core(force=ap.parse_args().force)
