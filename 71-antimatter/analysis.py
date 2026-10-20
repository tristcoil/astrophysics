"""Every result of the project, computed step by step and cached.

    python analysis.py           # use cached steps where present
    python analysis.py --force   # recompute everything (~25 min, mostly the fluid runs and BBN)

Each section writes its numbers into ``data/results.json`` and the curves the
figures need into ``data/derived/*.pkl``.
"""

from __future__ import annotations

import argparse
import json
import math
import pickle
import time

import numpy as np

import baryons
import config as cfg
import literature as lit

DOMAIN_DIR = cfg.DERIVED_DIR / "domains"
DOMAIN_DIR.mkdir(parents=True, exist_ok=True)


def _save_pickle(name: str, obj) -> None:
    with open(cfg.DERIVED_DIR / f"{name}.pkl", "wb") as f:
        pickle.dump(obj, f)


def load_pickle(name: str):
    with open(cfg.DERIVED_DIR / f"{name}.pkl", "rb") as f:
        return pickle.load(f)


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


# ---------------------------------------------------------------------------
# A. Photons and baryons; the excess back in time
# ---------------------------------------------------------------------------


def d3_ombh2() -> tuple[float, float, str]:
    try:
        r = json.loads(cfg.D3_RESULTS.read_text())
        p = r["fit"]["parameters"]["ombh2"]
        return float(p["value"]), float(p["error"]), "D3 results.json"
    except (OSError, KeyError, ValueError):
        return (*cfg.D3_OMBH2_FALLBACK, "D3 README (results file missing)")


def section_baryons() -> dict:
    import fetch_data as fd

    out = {"planck": baryons.numbers_today(cfg.OMBH2_PLANCK)}
    v, e, src = d3_ombh2()
    out["d3"] = {"ombh2": v, "ombh2_err": e, "source": src,
                 "eta": float(baryons.eta_from_ombh2(v)),
                 "eta_err": float(baryons.eta_from_ombh2(v) * e / v)}
    out["planck"]["eta_err"] = float(baryons.eta_from_ombh2(cfg.OMBH2_PLANCK_ERR))
    out["eta10_per_ombh2"] = baryons.eta10_coefficient()
    # the quark excess above the QCD crossover
    tt, grho, gs = fd.load_gstar()
    T = np.geomspace(0.16, 1e4, 120)
    g = np.interp(np.log(T), np.log(tt), gs)
    y_b = out["planck"]["Y_B"]
    excess = baryons.quark_excess(T, g, y_b)
    out["quark_excess"] = {
        "at_1GeV": float(np.interp(0.0, np.log(T), excess)),
        "at_200GeV": float(np.interp(math.log(200), np.log(T), excess)),
        "at_160MeV": float(excess[0]),
        "pairs_per_extra_quark_1GeV": float(0.5 / np.interp(0.0, np.log(T), excess)),
    }
    # the annihilation catastrophe and the real universe
    sym = baryons.nucleon_freezeout()
    real = baryons.nucleon_freezeout(eta=out["planck"]["eta"])
    s_over_ng = baryons.entropy_per_photon()
    i_f = int(np.argmin(np.abs(sym["Y_N"] / sym["Y_eq"] - 2)))
    variants = {sv: baryons.nucleon_freezeout(sigma_v_cm3s=sv)["Y_N"][-1] * s_over_ng
                for sv in (1.0e-15, 1.6e-15)}

    def ratio_at(t):
        return float(np.interp(-t, -real["T"], real["Y_Nbar"] / real["Y_N"]))

    out["catastrophe"] = {
        "sigma_v_cm3s": baryons.SIGMA_V_CM3S,
        "eta_symmetric": float(sym["Y_N"][-1] * s_over_ng),
        "Y_symmetric": float(sym["Y_N"][-1]),
        "T_freezeout_MeV": float(sym["T"][i_f] * 1e3),
        "x_freezeout": float(sym["x"][i_f]),
        "observed_over_symmetric": float(out["planck"]["eta"] / (sym["Y_N"][-1] * s_over_ng)),
        "eta_symmetric_range": [float(variants[1.6e-15]), float(variants[1.0e-15])],
        "real_antinucleons_per_nucleon": {f"{t * 1e3:.0f}MeV": ratio_at(t)
                                         for t in (0.05, 0.04, 0.03, 0.02)},
    }
    _save_pickle("baryons_curves", {"T_quark": T, "quark_excess": excess, "symmetric": sym,
                                    "real": real})
    return out


# ---------------------------------------------------------------------------
# B. Deuterium (BBN) against the CMB
# ---------------------------------------------------------------------------


def section_bbn(force: bool) -> dict:
    import bbn

    res, curves = bbn.analyse(force=force)
    cmb = []
    v, e, _ = d3_ombh2()
    cmb.append({"label": "D3 (this atlas): Planck TT", "ombh2": v, "err": e})
    for lab, val, err, src in lit.CMB_OMBH2:
        cmb.append({"label": lab, "ombh2": val, "err": err, "source": src})
    for c in cmb:
        c["eta"] = float(baryons.eta_from_ombh2(c["ombh2"]))
        c["eta_err"] = float(baryons.eta_from_ombh2(c["err"]))
    res["cmb"] = cmb
    p = res["rates"]["primat"]["pdg2026_12"]
    n = res["rates"]["nacreii"]["pdg2026_12"]
    for c in cmb:
        c["tension_primat_sigma"] = (c["ombh2"] - p["ombh2_median"]) / math.hypot(c["err"], p["ombh2_err"])
        c["tension_nacreii_sigma"] = (c["ombh2"] - n["ombh2_median"]) / math.hypot(c["err"], n["ombh2_err"])
    res["literature"] = [{"label": a, "ombh2": b, "err": c, "source": d}
                         for a, b, c, d in lit.BBN_OMBH2_LITERATURE]
    _save_pickle("bbn_curves", {k: {"eta": c.eta, "dh": c.dh, "yp": c.yp, "li7": c.li7,
                                    "he3": c.he3, "frac_theory": c.frac_theory,
                                    "mc_dh": c.mc_dh}
                                for k, c in curves.items()})
    return res


# ---------------------------------------------------------------------------
# C. Antimatter domains and the MeV sky
# ---------------------------------------------------------------------------


def _boundary_run(name: str, universe, force: bool, **kw):
    import domains as dm

    path = DOMAIN_DIR / f"sol_{name}.pkl"
    if path.exists() and not force:
        return pickle.load(open(path, "rb"))
    t0 = time.time()
    sol = dm.solve_boundary(universe, **kw)
    sol.settings["seconds"] = time.time() - t0
    pickle.dump(sol, open(path, "wb"))
    print(f"    boundary run {name}: {sol.settings['steps']} steps, {sol.settings['seconds']:.0f} s")
    return sol


def _band_average(E, F, lo, hi):
    out = []
    for a, b in zip(lo, hi):
        g = np.geomspace(a, b, 60)
        out.append(np.trapezoid(np.exp(np.interp(np.log(g), np.log(E), np.log(np.maximum(F, 1e-300)))), g)
                   / (b - a))
    return np.array(out)


def section_domains(force: bool) -> dict:
    import cdg98
    import domains as dm

    out = {}
    # -- photons from one annihilation ----------------------------------
    ph = dm.AnnihilationPhotons(lit.PBARP_CHANNELS, n_events=40000)
    n_gamma = lit.PBARP_MULTIPLICITY["n_gamma"][0]
    scale = n_gamma / ph.n_gamma

    class Photons:
        def __init__(self, norm):
            self.norm = norm

        def __call__(self, e):
            return self.norm * ph(e)

    bk = lit.BACKENSTOSS_DIGITISED
    c = np.cumsum(ph.dnde * np.gradient(ph.e_grid))
    c /= c[-1]
    out["photons"] = {"n_pi0_channels": ph.n_pi0, "n_pic_channels": ph.n_pic,
                      "n_gamma_phase_space": ph.n_gamma, "n_gamma_adopted": n_gamma,
                      "mean_E_MeV": ph.mean_e_gamma,
                      "median_E_MeV": float(np.interp(0.5, c, ph.e_grid)),
                      "frac_above_100MeV": float(1 - np.interp(100, ph.e_grid, c)),
                      "frac_above_200MeV": float(1 - np.interp(200, ph.e_grid, c)),
                      "backenstoss_ratio_30_500MeV": [
                          float(scale * ph(e) / np.interp(e, bk[:, 0], bk[:, 1]))
                          for e in (30, 70, 150, 300, 500)]}

    # -- reproducing Cohen, De Rujula & Glashow 1998 ----------------------
    fid = dm.cdg_fiducial()
    sol_cdg = _boundary_run("cdg_fiducial", fid, force, gbar_e=3.8, dedl_override_a=10.0)
    refJ = cdg98.annihilation_rate()["numerical"]
    refT = cdg98.temperatures()
    ys = np.array([1000, 500, 300, 200, 100, 70, 50, 40, 30, 25, 20.5])

    def at(sol_y, val, y):
        return np.exp(np.interp(np.log(y), np.log(sol_y[::-1]), np.log(val[::-1])))

    jr = np.exp(np.interp(np.log(ys), np.log(refJ["y"]), np.log(refJ["J_cm2_s"])))
    tr = np.interp(np.log(ys), np.log(refT["electron_heating"]["y"]), refT["electron_heating"]["T_eV"])
    out["cdg_check"] = {"y": ys, "J_ratio": at(sol_cdg.y, sol_cdg.J_cgs, ys) / jr,
                        "T_ours": at(sol_cdg.y, sol_cdg.T_far_ev, ys), "T_cdg": tr}
    g = cdg98.gamma_spectra()
    E_chk = np.array([0.3, 0.5, 1, 2, 3, 5, 10, 20, 30])
    ph38 = Photons(3.8 / ph.n_gamma)
    f_their_j = dm.diffuse_spectrum(fid, refJ["y"], refJ["J_cm2_s"] * 1e4, ph38, E_chk, d0_mpc=20.0)
    f_our_j = dm.diffuse_spectrum(fid, sol_cdg.y, sol_cdg.J, ph38, E_chk, d0_mpc=20.0)
    f_ref = np.exp(np.interp(np.log(E_chk), np.log(g["computed"]["E_MeV"]), np.log(g["computed"]["flux"])))
    out["cdg_check"].update({"E": E_chk, "spec_ratio_their_J": f_their_j / f_ref,
                             "spec_ratio_our_J": f_our_j / f_ref})

    # -- today's universe --------------------------------------------------
    u = dm.planck2018()
    runs = {"cdg": _boundary_run("planck_cdg_heating", u, force, heating="cdg"),
            "full": _boundary_run("planck_full_heating", u, force, heating="full"),
            "none": _boundary_run("planck_no_heating", u, force, heating=False)}
    photons = Photons(scale)
    Eg = np.geomspace(0.08, 120, 90)
    spectra = {k: dm.diffuse_spectrum(u, s.y, s.J, photons, Eg, d0_mpc=1000.0)
               for k, s in runs.items()}
    comptel = {w: lit.comptel(w) for w in ("weidenspointner", "kappadath")}
    limits = {}
    for heat, F in spectra.items():
        for w, cpt in comptel.items():
            pred = _band_average(Eg, F, cpt["E_lo"], cpt["E_hi"])
            ratio = pred / (cpt["I"] + 2 * cpt["err"])
            k = int(np.argmax(ratio))
            limits[f"{heat}/{w}"] = {"d0_min_gpc": float(ratio.max()),
                                     "band_MeV": [float(cpt["E_lo"][k]), float(cpt["E_hi"][k])],
                                     "ratio_per_band": ratio}
    # epochs of annihilation: CDG count 20 < y < 1100
    ywin = {}
    for ymin in (10.0, 30.0, 50.0):
        F = dm.diffuse_spectrum(u, runs["cdg"].y, runs["cdg"].J, photons, Eg, d0_mpc=1000.0,
                                y_min=ymin)
        cpt = comptel["weidenspointner"]
        ywin[f"y_min_{ymin:.0f}"] = float(np.max(_band_average(Eg, F, cpt["E_lo"], cpt["E_hi"])
                                                 / (cpt["I"] + 2 * cpt["err"])))
    # one flat boundary at distance r from us
    r20 = u.comoving_distance(20.0) / cfg.GPC
    r1100 = u.comoving_distance(1100.0) / cfg.GPC
    planes = {}
    for rp in (0.0, 8.0, r20, 11.5, 12.0, 12.5, 13.0, 13.5):
        F = dm.diffuse_spectrum(u, runs["cdg"].y, runs["cdg"].J, photons, Eg,
                                plane_distance_gpc=rp)
        cpt = comptel["weidenspointner"]
        planes[f"{rp:.2f}"] = float(np.max(_band_average(Eg, F, cpt["E_lo"], cpt["E_hi"])
                                           / (cpt["I"] + 2 * cpt["err"])))
    # how much matter would annihilate (per comoving volume, 20 < y < 1100)
    s = runs["cdg"]
    yy = np.geomspace(20, 1100, 2000)
    jj = np.exp(np.interp(np.log(yy), np.log(s.y[::-1]), np.log(s.J[::-1])))
    frac_1gpc = float(np.trapezoid(jj / (cfg.GPC * yy**3 * u.hubble(yy) * u.n_b(1.0)), yy))
    out["modern"] = {
        "eta": u.eta, "yp": u.yp, "h": u.h, "omega_m": u.omega_m,
        "limits": limits,
        "headline_d0_min_gpc": limits["cdg/weidenspointner"]["d0_min_gpc"],
        "y_window": ywin,
        "comoving_distance_y20_gpc": r20, "comoving_distance_y1100_gpc": r1100,
        "single_plane_excess": planes,
        "fraction_annihilated_d0_1gpc": frac_1gpc,
        "J_at": {k: {str(y): float(at(v.y, v.J_cgs, y)) for y in (1000, 200, 50, 20.5)}
                 for k, v in runs.items()},
        "T_at": {k: {str(y): float(at(v.y, v.T_far_ev, y)) for y in (1000, 200, 50, 20.5)}
                 for k, v in runs.items()},
        "run_seconds": {k: v.settings.get("seconds") for k, v in runs.items()},
    }
    # spectrum CDG would have predicted for 1 Gpc in their universe
    f_cdg_1gpc = g["computed"]["flux"] * 20.0 / 1000.0
    cpt = comptel["weidenspointner"]
    pred_cdg = _band_average(g["computed"]["E_MeV"], f_cdg_1gpc, cpt["E_lo"], cpt["E_hi"])
    out["cdg_1998_vs_final_comptel_d0_min_gpc"] = float(np.max(pred_cdg / (cpt["I"] + 2 * cpt["err"])))
    _save_pickle("domain_curves", {
        "Eg": Eg, "spectra_1gpc": spectra, "runs": {k: {"y": v.y, "J": v.J_cgs, "T": v.T_far_ev,
                                                        "D": v.depletion_m, "profiles": v.profiles}
                                                    for k, v in runs.items()},
        "cdg_run": {"y": sol_cdg.y, "J": sol_cdg.J_cgs, "T": sol_cdg.T_far_ev},
        "photons": {"E": ph.e_grid, "dnde": scale * ph.dnde, "backenstoss": bk},
        "cdg_ref": {"J": refJ, "T": refT, "spec": g},
        "cdg_spec_check": {"E": E_chk, "their_J": f_their_j, "our_J": f_our_j, "ref": f_ref}})
    return out


# ---------------------------------------------------------------------------
# D. Cosmic-ray antiprotons
# ---------------------------------------------------------------------------


def section_cosmic_rays(force: bool) -> dict:
    import cosmic_rays as cr
    import fetch_data as fd
    import pbar_production as pp

    lb = cr.LeakyBox()
    fit = lb.fit_bc(r_min=4.0)
    p = fit["params"]
    out = {"bc_fit": {"X0_10GV": float(p[0]), "X0_err": float(math.sqrt(fit["cov"][0, 0])),
                      "delta": float(p[1]), "delta_err": float(math.sqrt(fit["cov"][1, 1])),
                      "d_delta": float(p[2]), "d_delta_err": float(math.sqrt(fit["cov"][2, 2])),
                      "chi2": fit["chi2"], "ndof": fit["ndof"],
                      "X_at": {str(r): float(cr.grammage(r, *p)) for r in (10, 100, 1000)}}}
    R = np.geomspace(1.0, 600.0, 90)
    preds = {}
    for period, sub in (("2011-2015", "AMS02 (2011/05-2015/05)"), ("2011-2022", "AMS02 (2011/05-2022/06)")):
        d = cr.load_ratio("pbar_over_p", sub)
        phi = float(d.phi_mv.iloc[0]) / 1e3
        for name in ("II", "I"):
            pp.PARAMETRISATION["name"] = name
            path = cfg.DERIVED_DIR / "pbar_kernels.pkl"
            if path.exists():
                path.unlink()                     # kernels depend on the parametrisation
            model = lb.pbar_over_p(R, p, phi)
            at_data = lb.pbar_over_p(d.e_mean.values, p, phi)
            preds[f"{period}/{name}"] = {"R": R, "model": model, "phi": phi,
                                         "data_R": d.e_mean.values, "data": d.value.values,
                                         "data_err": d.err.values, "model_at_data": at_data}
        pp.PARAMETRISATION["name"] = "II"
    # model/data above 20 GV (where both parametrisations agree)
    summ = {}
    for k, v in preds.items():
        m = v["data_R"] > 20
        ratio = v["model_at_data"][m] / v["data"][m]
        w = 1 / (v["data_err"][m] / v["data"][m]) ** 2
        summ[k] = {"mean_model_over_data_above_20GV": float(np.sum(w * ratio) / np.sum(w)),
                   "rms": float(np.std(ratio)),
                   "chi2_above_20GV": float(np.sum(((v["model_at_data"][m] - v["data"][m])
                                                    / v["data_err"][m]) ** 2)),
                   "n_above_20GV": int(m.sum())}
    out["pbar_over_p"] = summ
    # high-energy values
    d16 = cr.load_ratio("pbar_over_p", "AMS02 (2011/05-2015/05)")
    hi = d16.e_mean > 60
    w = 1 / d16.err[hi] ** 2
    out["ams2016_pbar_p_above_60GV"] = {"mean": float(np.sum(w * d16.value[hi]) / np.sum(w)),
                                         "err": float(np.sum(w) ** -0.5)}
    # antihelium
    ah = fd.load_crdb("hebar_over_he")
    best = ah.loc[ah.value.idxmin()]
    out["antihelium_best_limit"] = {"value": float(best.value), "experiment": best.subexp,
                                    "R_range_GV": [float(best.e_low), float(best.e_high)],
                                    "ads": best.ads}
    _save_pickle("cr_curves", {"fit": fit, "preds": preds, "bc_model_R": np.geomspace(2, 3000, 80),
                               "bc_model": lb.bc_toa(np.geomspace(2, 3000, 80), p)})
    return out


# ---------------------------------------------------------------------------
# E. Galaxy clusters: Steigman's argument with Fermi-LAT
# ---------------------------------------------------------------------------


def section_clusters() -> dict:
    import domains as dm

    z = 0.296
    kT = 14.1                                    # keV (Markevitch 2006, 500 ks)
    # bremsstrahlung emissivity in the emitted 2-10 keV band, per n_e sum(n_i Z^2)
    # (Rybicki & Lightman eq. 5.14b, Born-approximation Gaunt factor)
    T_k = kT * 1e3 * cfg.EV / cfg.K_B
    E = np.linspace(2 * (1 + z), 10 * (1 + z), 400)                 # keV
    nu = E * 1e3 * cfg.EV / (2 * math.pi * cfg.HBAR)
    gaunt = math.sqrt(3) / math.pi * np.log(4 / 1.781 * kT / E)
    eps_nu = 6.8e-38 * T_k**-0.5 * np.exp(-E / kT) * np.maximum(gaunt, 0.1)  # erg s^-1 cm^3 Hz^-1
    lam_x = float(np.trapezoid(eps_nu, nu))                                  # erg s^-1 cm^3
    # thermal <sigma v> for p pbar with Coulomb focusing (Sommerfeld factor),
    # (sigma v)_0 = 44 mb x c from OBELIX
    sv0 = baryons.SIGMA_V_CM3S
    mu = 0.5 * cfg.M_P_EV * cfg.EV / cfg.C_LIGHT**2
    v = np.linspace(1e3, 3e7, 20000)                                         # m/s
    fmb = v**2 * np.exp(-mu * v**2 / (2 * kT * 1e3 * cfg.EV))
    xi = 2 * math.pi * cfg.ALPHA * cfg.C_LIGHT / v
    sommer = xi / (1 - np.exp(-xi))
    sv = sv0 * float(np.trapezoid(fmb * sommer, v) / np.trapezoid(fmb, v))
    # photons per annihilation in the emitted LAT band 0.2-100 GeV (observed)
    ph = dm.AnnihilationPhotons(lit.PBARP_CHANNELS, n_events=40000)
    norm = lit.PBARP_MULTIPLICITY["n_gamma"][0] / ph.n_gamma
    eg = np.geomspace(200 * (1 + z), 1e5 * (1 + z), 300)
    n_band = float(norm * np.trapezoid(ph(eg), eg))
    eg100 = np.geomspace(100 * (1 + z), 1e5 * (1 + z), 300)
    n_100 = float(norm * np.trapezoid(ph(eg100), eg100))
    comp = (1 + 2 * (0.25 / 4) / 0.75) * (1 + 4 * (0.25 / 4) / 0.75)       # n_e sum Z^2 n_i / n_p^2
    f_lat = 2.75e-9                 # Ackermann et al. 2010, Table 2: 0.2-100 GeV, 95% CL
    out = {"kT_keV": kT, "z": z, "lambda_x_2_10keV": lam_x, "sigma_v_thermal": sv,
           "sommerfeld_mean": sv / sv0, "n_gamma_0p2_100GeV": n_band, "n_gamma_above_100MeV": n_100,
           "fermi_limit_0p2_100GeV": f_lat}
    for label, fx in (("steigman_FX", 4.7e-12), ("FX_1e-11", 1.0e-11)):
        # F_N / F_E = (Ndot / L) (1 + z):  f = comp lam_x (F_N/F_E) / ((1+z) N_gamma <sv>)
        out[f"f_limit_{label}"] = comp * lam_x * (f_lat / fx) / ((1 + z) * n_band * sv)
    # Steigman 2008's own formula with the EGRET limit, as a check
    out["steigman2008_check"] = 2.6e-10 * (kT / 8.6) * (3.5e-8 / 4.7e-12)
    return out


# ---------------------------------------------------------------------------
# F. Sakharov's conditions
# ---------------------------------------------------------------------------


def section_sakharov() -> dict:
    from scipy import optimize

    import sakharov as sk

    c = lit.CKM["pdg2026"]
    J, ang = sk.jarlskog_from_wolfenstein(c["lambda"][0], c["A"][0], c["rhobar"][0], c["etabar"][0])
    lat = lit.LATTICE_EW

    def gamma_over_t3(T):
        b = lat["sphaleron_broken"]
        return np.exp(b["a"] * T - b["b"]) * T

    def hub(T):
        return sk.hubble_rate_gev(T, lat["freezeout_gstar"])

    t_star = optimize.brentq(lambda T: gamma_over_t3(T) - lat["freezeout_alpha"] * hub(T), 110, 159)
    meas = {}
    for lab, m in (("pdg2026", lit.QUARK_MASSES_PDG2026), ("running_mz", lit.QUARK_MASSES_MZ)):
        up, dn = [m["u"], m["c"], m["t"]], [m["d"], m["s"], m["b"]]
        meas[lab] = {f"T{T:g}": float(sk.jarlskog_measure(c["J"][0], up, dn, T))
                     for T in (100.0, round(t_star, 1), lat["T_c"][0])}
    mh = np.linspace(20, 130, 111)
    vt = np.array([sk.transition_strength(m)[1] for m in mh])
    # the simplest classic estimate v/T = 2E/lambda (no thermal logs)
    E = (2 * lit.EW_MASSES["W"] ** 3 + lit.EW_MASSES["Z"] ** 3) / (4 * math.pi * sk.V_EW**3)
    vt_simple = 2 * E / (mh**2 / (2 * sk.V_EW**2))
    m_washout_simple = float(np.interp(-1.0, -vt_simple, mh))
    out = {"J_from_wolfenstein": J, "J_pdg": c["J"][0], "ckm_angles": ang,
           "T_star_ours": float(t_star), "T_star_paper": lat["T_star"][0],
           "jarlskog_measure": meas, "Y_B_needed": baryons.numbers_today(cfg.OMBH2_PLANCK)["Y_B"],
           "v_over_T_at_125": float(np.interp(125.2, mh, vt)),
           "v_over_T_simple_at_125": float(np.interp(125.2, mh, vt_simple)),
           "m_h_washout_simple": m_washout_simple,
           "H_over_T_at_Tc": float(hub(lat["T_c"][0]) / lat["T_c"][0])}
    _save_pickle("sakharov_curves", {"mh": mh, "vt": vt, "vt_simple": vt_simple,
                                     "T": np.linspace(120, 175, 200)})
    return out


# ---------------------------------------------------------------------------


def run(force: bool = False, sections=None) -> dict:
    res = json.loads(cfg.RESULTS_PATH.read_text()) if cfg.RESULTS_PATH.exists() and not force else {}
    todo = sections or ["baryons", "bbn", "domains", "cosmic_rays", "clusters", "sakharov"]
    fns = {"baryons": lambda: section_baryons(), "bbn": lambda: section_bbn(force),
           "domains": lambda: section_domains(force), "cosmic_rays": lambda: section_cosmic_rays(force),
           "clusters": lambda: section_clusters(), "sakharov": lambda: section_sakharov()}
    for name in todo:
        if name in res and not force and sections is None:
            continue
        t0 = time.time()
        print(f"  analysis: {name} ...", flush=True)
        res[name] = fns[name]()
        print(f"  analysis: {name} done ({time.time() - t0:.0f} s)", flush=True)
        res["_generated"] = time.strftime("%Y-%m-%d %H:%M:%S")
        cfg.RESULTS_PATH.write_text(json.dumps(res, indent=1, default=_json_default))
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args()
    run(force=a.force, sections=a.only)
