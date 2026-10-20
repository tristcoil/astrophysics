"""Every number of the project: the equation of state, its checks, neutrino decoupling and the calendar.

    python analysis.py      # writes data/results.json and data/derived/curves.npz (about a minute)

Run order: equation of state (ideal and lattice builds, plus the variants that make our
systematic band) -> comparisons with Saikawa & Shirai 2018, Laine-Meyer-Schroeder 2015 and
Borsanyi et al. 2016 -> the QCD crossover in detail -> neutrino decoupling -> the calendar.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

import config as C  # noqa: E402
import cosmic_time as ct  # noqa: E402
import eos  # noqa: E402
import neutrinos  # noqa: E402
import particles  # noqa: E402
import qcd  # noqa: E402
import thermo  # noqa: E402

PI2 = np.pi ** 2


def _f(x, nd=6):
    """Round for JSON (floats only)."""
    if x is None:
        return None
    x = float(x)
    if not np.isfinite(x):
        return None
    return float(f"{x:.{nd}g}")


def _at(T, y, Tq):
    return float(np.interp(np.log(Tq), np.log(T), y))


# --- 1. The equation of state ----------------------------------------------------------------


QCD_VARIANTS = {
    "hrg_switch_110MeV": {"qcd_kwargs": {"t_hrg": 0.110}},
    "hrg_switch_130MeV": {"qcd_kwargs": {"t_hrg": 0.130}},
    "pqcd_blend_0.7GeV": {"qcd_kwargs": {"t_pqcd": 0.7}},
    "pqcd_blend_1.5GeV": {"qcd_kwargs": {"t_pqcd": 1.5}},
    # scale mu = (0.5 ... 2) x 2 pi T, each with q_c refitted so that pQCD meets the lattice at 1 GeV
    "pqcd_scale_x0.5_matched": {"qcd_kwargs": {"mu_factor": 0.5, "qc": None}},
    "pqcd_scale_x2_matched": {"qcd_kwargs": {"mu_factor": 2.0, "qc": None}},
    "pqcd_qc_matched": {"qcd_kwargs": {"qc": None}},
    "hadrons_with_f0_500": {"hadron_kwargs": {"keep_broad_scalars": True}},
    # charm and bottom thresholds of perturbative QCD at the pole masses (1.67, 4.78 GeV), q_c refitted
    "pole_quark_masses": {"qcd_kwargs": {"masses": {4: 1.67, 5: 4.78}, "qc": None}},
}


def equation_of_state(T):
    t0 = time.time()
    run = eos.decoupling_run()
    lat = eos.build(T, "lattice", run=run)
    ideal = eos.build(T, "ideal")
    pions = eos.build(T, "ideal_pions")
    variants = {}
    for name, kw in QCD_VARIANTS.items():
        v = eos.build(T, "lattice", run=run, **kw)
        variants[name] = {"g_rho": v["g_rho"], "g_s": v["g_s"]}
    lo = np.min([v["g_rho"] for v in variants.values()] + [lat["g_rho"]], axis=0)
    hi = np.max([v["g_rho"] for v in variants.values()] + [lat["g_rho"]], axis=0)
    lo_s = np.min([v["g_s"] for v in variants.values()] + [lat["g_s"]], axis=0)
    hi_s = np.max([v["g_s"] for v in variants.values()] + [lat["g_s"]], axis=0)
    print(f"  equation of state: {len(T)} temperatures, {len(QCD_VARIANTS)} variants ({time.time() - t0:.0f} s)")
    return lat, ideal, pions, variants, (lo, hi, lo_s, hi_s)


def half_temperatures(sm) -> dict:
    """Where each heavy species' share of g*s falls to half its massless value (the step labels)."""
    from scipy.optimize import brentq

    out = {}
    for sp in sm:
        if sp.mass <= 0:
            continue
        f = lambda lx: float(thermo.functions(np.exp(lx), sp.stat)["s"][0]
                             / (thermo.MASSLESS[sp.stat][0] + thermo.MASSLESS[sp.stat][1]) - 0.5)
        x_half = float(np.exp(brentq(f, np.log(1e-3), np.log(50))))
        out[sp.name] = {"mass_GeV": _f(sp.mass), "x_half": _f(x_half), "T_half_GeV": _f(sp.mass / x_half)}
    pion = thermo.Species("pi", "π", "meson", particles.pdg_mass(211), 3, +1, 211)
    f = lambda lx: float(thermo.functions(np.exp(lx), +1)["s"][0] / (4 * PI2 / 90) - 0.5)
    xh = float(np.exp(brentq(f, np.log(1e-3), np.log(50))))
    out["pions"] = {"mass_GeV": _f(pion.mass), "x_half": _f(xh), "T_half_GeV": _f(pion.mass / xh)}
    return out


def comparisons(T, lat, ideal, band) -> dict:
    ss, lms, b16 = eos.saikawa_shirai(), eos.laine_meyer_schroeder(), qcd.borsanyi2016_spline_table()
    S, Se = eos.on_grid(ss, "g_rho", T), eos.on_grid(ss, "g_rho_err", T)
    Ss, Sse = eos.on_grid(ss, "g_s", T), eos.on_grid(ss, "g_s_err", T)
    Lm, Lms = eos.on_grid(lms, "g_rho", T), eos.on_grid(lms, "g_s", T)
    out = {"ranges": {}}
    ranges = {"10 TeV - 200 GeV": (200, 1e4), "200 - 10 GeV": (10, 200), "10 - 1 GeV": (1, 10),
              "1 GeV - 200 MeV": (0.2, 1), "200 - 120 MeV": (0.12, 0.2), "120 - 10 MeV": (0.01, 0.12),
              "10 MeV - 10 keV": (1e-5, 0.01)}
    for name, (a, b) in ranges.items():
        m = (T >= a) & (T <= b)
        d_lat = (lat["g_rho"][m] - S[m]) / S[m]
        d_id = (ideal["g_rho"][m] - S[m]) / S[m]
        pull = (lat["g_rho"][m] - S[m]) / Se[m]
        d_lms = (lat["g_rho"][m] - Lm[m]) / Lm[m]
        out["ranges"][name] = {
            "max_abs_dev_lattice_vs_SS_pct": _f(100 * np.nanmax(np.abs(d_lat)), 3),
            "max_abs_pull_lattice_vs_SS": _f(np.nanmax(np.abs(pull)), 3),
            "max_abs_dev_ideal_vs_SS_pct": _f(100 * np.nanmax(np.abs(d_id)), 3),
            "max_abs_dev_lattice_vs_LMS_pct": _f(100 * np.nanmax(np.abs(d_lms)), 3),
            "SS_relative_error_pct_max": _f(100 * np.nanmax(Se[m] / S[m]), 3),
        }
    within = np.abs(lat["g_rho"] - S) <= Se
    out["fraction_of_grid_within_SS_error"] = _f(np.nanmean(within[np.isfinite(S)]), 4)
    within_s = np.abs(lat["g_s"] - Ss) <= Sse
    out["fraction_of_grid_within_SS_error_gs"] = _f(np.nanmean(within_s[np.isfinite(Ss)]), 4)
    # ideal vs S&S: the biggest misses (the QCD crossover)
    r = ideal["g_rho"] / S
    i = int(np.nanargmax(r))
    out["ideal_max_ratio_to_SS"] = {"T_GeV": _f(T[i]), "ratio": _f(r[i], 4)}
    j = int(np.nanargmin(r))
    out["ideal_min_ratio_to_SS"] = {"T_GeV": _f(T[j]), "ratio": _f(r[j], 4)}
    # Borsanyi 2016 knots
    knots = []
    for Tk, gk, gsk in zip(b16["T"], b16["g_rho"], b16["g_s"]):
        knots.append({"T_GeV": _f(Tk), "g_rho_B16": _f(gk, 5), "g_rho_ours": _f(_at(T, lat["g_rho"], Tk), 5),
                      "g_s_B16": _f(gsk, 5), "g_s_ours": _f(_at(T, lat["g_s"], Tk), 5),
                      "g_rho_SS": _f(_at(T, np.nan_to_num(S), Tk), 5)})
    out["borsanyi2016_knots"] = knots
    lo, hi, lo_s, hi_s = band
    out["our_band_halfwidth_pct_max"] = _f(100 * np.max((hi - lo) / 2 / lat["g_rho"]), 3)
    i = int(np.argmax((hi - lo) / lat["g_rho"]))
    out["our_band_halfwidth_pct_max_at_T"] = _f(T[i])
    return out


def checks(T, lat, ideal) -> dict:
    sm = particles.standard_model()
    g_all = sum(sp.g_eff_massless for sp in sm)
    run_inst = neutrinos.solve(collisions=False, qed_order=0, T_start=0.2)
    gi = neutrinos.low_T_gstar(run_inst)
    gl = {"g_rho": float(lat["g_rho"][0]), "g_s": float(lat["g_s"][0])}
    return {
        "g_high_count": {"ours": g_all, "spec": 106.75},
        "g_at_10TeV": {"ideal": _f(_at(T, ideal["g_rho"], 1e4), 6), "lattice": _f(_at(T, lat["g_rho"], 1e4), 6),
                       "SS2018": _f(_at(T, eos.on_grid(eos.saikawa_shirai(), "g_rho", T), 1e4), 6),
                       "SS2018_err": _f(_at(T, eos.on_grid(eos.saikawa_shirai(), "g_rho_err", T), 1e4), 3),
                       "LMS2015": _f(_at(T, eos.on_grid(eos.laine_meyer_schroeder(), "g_rho", T), 1e4), 6)},
        "g_few_MeV": {"ideal_at_5MeV": _f(_at(T, ideal["g_rho"], 5e-3), 5),
                      "ideal_at_3MeV": _f(_at(T, ideal["g_rho"], 3e-3), 5),
                      "massless_electron_limit": 2 + 7 / 8 * (4 + 6), "spec": 10.75},
        "today_instantaneous": {"g_rho": _f(gi["g_rho"], 6), "g_s": _f(gi["g_s"], 6),
                                "T_gamma_over_T_nu": _f(run_inst["ratio_e"], 7),
                                "exact_11_4": _f((11 / 4) ** (1 / 3), 7), "spec": [3.36, 3.91]},
        "today_full": {"g_rho": _f(gl["g_rho"], 6), "g_s": _f(gl["g_s"], 6),
                       "SS2018": [3.3830836, 3.9309363], "EJLS2025_no_osc": [3.3828, 3.9302],
                       "EJLS2025_osc_full": [3.3827, 3.9302]},
    }


# --- 2. The QCD crossover ---------------------------------------------------------------------


def qcd_crossover(T, lat, ideal) -> dict:
    out = {}
    hot, wb = qcd.hotqcd2014(), qcd.wb2014()
    hs = particles.hadrons()
    # where the entropy of QCD matter rises fastest (max of d(s/T^3)/dT): lattice tables and fit
    from scipy.interpolate import CubicSpline

    def steepest(Tt, sT3):
        sp = CubicSpline(Tt, sT3)
        Tf = np.linspace(Tt[0], Tt[-1], 4001)
        return float(Tf[np.argmax(sp(Tf, 1))])

    Tf = np.linspace(0.12, 0.30, 1801)
    d16 = qcd.borsanyi2016_delta(Tf, "2+1+1")
    # s/T^3 of the 2+1+1 fit needs p: integrate on the full strong-sector grid instead
    strong = qcd.strong_sector(T)
    sT3 = strong["s"]
    m = (T > 0.12) & (T < 0.3)
    out["steepest_entropy_rise_MeV"] = {
        "HotQCD2014_table": _f(1000 * steepest(hot["T"], hot["s"]), 4),
        "WB2014_table": _f(1000 * steepest(wb["T"], wb["s"]), 4),
        "B16_fit_our_integration": _f(1000 * steepest(T[m], sT3[m]), 4),
        "T_pc_HotQCD2019": 156.5, "T_c_Borsanyi2020": 158.0,
    }
    out["trace_anomaly_peak_MeV"] = {"HotQCD2014_table": _f(1000 * hot["T"][np.argmax(hot["delta"])], 4),
                                     "WB2014_table": _f(1000 * wb["T"][np.argmax(wb["delta"])], 4),
                                     "B16_fit": _f(1000 * Tf[np.argmax(d16)], 4)}
    # hadron gas against lattice (HotQCD 2014 s/T^3 with its asymmetric errors)
    h = qcd.hrg(hot["T"], hs)
    dev = (h["s"] - hot["s"]) / np.where(h["s"] > hot["s"], hot["s_hi"], hot["s_lo"])
    first = next((float(t) for t, d in zip(hot["T"], dev) if abs(d) > 1), None)
    out["hrg_vs_hotqcd"] = {"T_MeV": [_f(1000 * t, 4) for t in hot["T"][:12]],
                            "s_over_T3_hrg": [_f(x, 4) for x in h["s"][:12]],
                            "s_over_T3_hotqcd": [_f(x, 4) for x in hot["s"][:12]],
                            "pull": [_f(x, 3) for x in dev[:12]],
                            "first_T_beyond_1sigma_MeV": _f(1000 * first, 4) if first else None}
    wbh = qcd.hrg(wb["T"], hs)
    devw = (wbh["s"] - wb["s"]) / wb["s_err"]
    firstw = next((float(t) for t, d in zip(wb["T"], devw) if abs(d) > 1), None)
    out["hrg_vs_wb2014_first_T_beyond_1sigma_MeV"] = _f(1000 * firstw, 4) if firstw else None
    # how many resonances the hadron gas needs: g_rho of the strong sector at 150 MeV
    T150 = np.array([0.150])
    cuts = {}
    for mmax in (0.2, 0.6, 1.0, 1.5, 2.0, 3.0):
        hcut = [x for x in hs if x.mass <= mmax]
        cuts[f"{mmax:.1f}"] = {"states": float(sum(x.g for x in hcut)),
                               "g_rho_strong_150MeV": _f(qcd.hrg(T150, hcut)["rho"][0] * 30 / PI2, 4)}
    out["hrg_mass_cut_at_150MeV"] = cuts
    out["lattice_g_rho_strong_150MeV"] = _f(_at(T, strong["rho"] * 30 / PI2, 0.150), 4)
    out["free_quark_gluon_g_rho_150MeV"] = _f(16 + 7 / 8 * 12 * 2 + sum(
        thermo.gstar(sp, np.array([0.150]))["rho"][0]
        for sp in particles.standard_model() if sp.name in ("strange", "charm")), 4)
    out["hadron_list"] = particles.hadron_summary(hs)
    # the abrupt textbook switch at 156.5 MeV: g just above and below
    edge = eos.build(np.array([C.T_SWITCH_IDEAL * 1.0001, C.T_SWITCH_IDEAL * 0.9999]), "ideal")
    edge_pi = eos.build(np.array([C.T_SWITCH_IDEAL * 0.9999]), "ideal_pions")
    out["ideal_jump_at_switch"] = {"above": _f(edge["g_rho"][0], 4), "below": _f(edge["g_rho"][1], 4),
                                   "below_pions_only": _f(edge_pi["g_rho"][0], 4),
                                   "lattice": _f(_at(T, lat["g_rho"], C.T_SWITCH_IDEAL), 4)}
    return out


# --- 3. Neutrino decoupling ----------------------------------------------------------------------


def neutrino_decoupling() -> dict:
    t0 = time.time()
    rows = {}
    cases = {
        "instantaneous_textbook": dict(collisions=False, qed_order=0, T_start=0.2),
        "instantaneous_10MeV_qed": dict(collisions=False, qed_order=3),
        "fd_noqed_separate": dict(rates="fd", qed_order=0),
        "fd_qed2_separate": dict(rates="fd", qed_order=2),
        "fd_qed3_separate": dict(rates="fd", qed_order=3),
        "fd_qed3_common": dict(rates="fd", qed_order=3, separate=False),
        "fdme_qed3_separate": dict(rates="fd_me", qed_order=3),
        "fdme_qed3_common": dict(rates="fd_me", qed_order=3, separate=False),
        "fdme_qed2_common": dict(rates="fd_me", qed_order=2, separate=False),
        "fdme_qed3_separate_rates_x0.9": dict(rates="fd_me", rate_scale=0.9),
        "fdme_qed3_separate_rates_x1.1": dict(rates="fd_me", rate_scale=1.1),
    }
    runs = {}
    for name, kw in cases.items():
        r = neutrinos.solve(**kw)
        runs[name] = r
        g = neutrinos.low_T_gstar(r)
        rows[name] = {"T_gamma_over_T_nu_e": _f(r["ratio_e"], 6), "T_gamma_over_T_nu_mu": _f(r["ratio_mu"], 6),
                      "T_gamma_over_T_nu_common": _f(r["ratio_common"], 6), "N_eff": _f(r["N_eff"], 6),
                      "g_rho_today": _f(g["g_rho"], 6), "g_s_today": _f(g["g_s"], 6)}
    base = runs["fdme_qed3_separate"]
    td = neutrinos.decoupling_temperatures()
    ratio = base["ratio_common"]
    T_nu_inst = C.T0_K * (4 / 11) ** (1 / 3)
    n_flavour = lambda Tk: 3 * C.ZETA3 / (2 * PI2) * (C.K_B_GEV_PER_K * Tk / C.HBARC_GEV_M / 100) ** 3
    out = {
        "cases": rows,
        "reference": {"Escudero2020_fd_nlo_common": [1.3954, 3.049], "Escudero2020_fd_nlo_separate": [1.3941, 1.3962, 3.048],
                      "Escudero2020_fdme_nlo_common": [1.39578, 3.045],
                      "EJLS2025_no_osc": [3.0443, 3.3828, 3.9302], "EJLS2025_osc_mu0": [3.0453, 3.3832, 3.9307],
                      "Bennett2021": [3.0440, 0.0002], "Fortepiano_no_osc": 3.0435, "Fortepiano_osc": 3.0439},
        "decoupling_T_MeV": {k: _f(1e3 * v, 4) for k, v in td.items()},
        "T_nu_today_K": {"instantaneous": _f(T_nu_inst, 5),
                         "nu_e": _f(C.T0_K / base["ratio_e"], 5), "nu_mu_tau": _f(C.T0_K / base["ratio_mu"], 5)},
        "n_nu_today_cm3": {"per_flavour_instantaneous": _f(n_flavour(T_nu_inst), 5),
                           "all_instantaneous": _f(3 * n_flavour(T_nu_inst), 5),
                           "all_with_heating": _f(n_flavour(C.T0_K / base["ratio_e"]) + 2 * n_flavour(C.T0_K / base["ratio_mu"]), 5)},
        "n_gamma_today_cm3": _f(2 * C.ZETA3 / PI2 * (C.K_B_GEV_PER_K * C.T0_K / C.HBARC_GEV_M / 100) ** 3, 5),
        "N_eff_ours": _f(base["N_eff"], 5), "N_eff_ours_common": _f(runs["fdme_qed3_common"]["N_eff"], 5),
        "rate_sensitivity_dNeff_per_10pct": _f((runs["fdme_qed3_separate_rates_x0.9"]["N_eff"]
                                                - runs["fdme_qed3_separate_rates_x1.1"]["N_eff"]) / 2, 3),
        "T_gamma_over_T_nu_common_ours": _f(ratio, 6),
    }
    print(f"  neutrino decoupling: {len(cases)} runs ({time.time() - t0:.0f} s)")
    return out, runs


# --- 4. The calendar -----------------------------------------------------------------------------


def eta_baryon(cosmo=C.PLANCK18) -> float:
    """Baryon-to-photon ratio today from omega_b (mean mass per baryon for Y_p = 0.245)."""
    n_gamma = 2 * C.ZETA3 / PI2 * C.T0_GEV ** 3
    m_per_baryon = 0.938272 * (1 - 0.245) + 3.727379 / 4 * 0.245
    return cosmo["ombh2"] * ct.RHO_CRIT_H2 / m_per_baryon / n_gamma


def calendar(lat, cal, halves) -> dict:
    T, t = cal["T"], cal["t_s"]
    m_n_minus_p = particles.pdg_mass(2112) - particles.pdg_mass(2212)
    eta = eta_baryon()
    B_D = 2.224566e-3  # deuteron binding energy (GeV), CODATA/AME
    td = neutrinos.decoupling_temperatures()
    z_star = 1089.80  # Planck 2018 TT,TE,EE+lowE+lensing+BAO
    i_eq = cal["i_equality"]
    z_eq = T[i_eq] / C.T0_GEV * (cal["g_s"][i_eq] / cal["g_s"][0]) ** (1 / 3) - 1
    rho_over_n = lambda TT: float(np.interp(np.log(TT), np.log(lat["T"]), lat["rho"] / lat["T"]))  # GeV^3
    events = [
        ("ew_crossover", "Electroweak crossover: the Higgs field switches on", C.LIT["T_ew_donofrio2016"]["value"],
         "D'Onofrio & Rummukainen 2016"),
        ("top_annihilation", "Top quarks annihilate (half gone)", halves["top"]["T_half_GeV"], "this work"),
        ("W_annihilation", "W bosons gone (half)", halves["W"]["T_half_GeV"], "this work"),
        ("bottom_annihilation", "Bottom quarks annihilate (half gone)", halves["bottom"]["T_half_GeV"], "this work"),
        ("tau_annihilation", "Tau leptons annihilate (half gone)", halves["tau"]["T_half_GeV"], "this work"),
        ("charm_annihilation", "Charm quarks annihilate (half gone, free-quark estimate)", halves["charm"]["T_half_GeV"], "this work"),
        ("qcd_crossover", "Quark-hadron crossover: protons, neutrons and pions form",
         C.LIT["T_pc_hotqcd2019"]["value"], "HotQCD 2019"),
        ("muon_annihilation", "Muons annihilate (half gone)", halves["muon"]["T_half_GeV"], "this work"),
        ("pion_annihilation", "Pions decay and annihilate (half gone)", halves["pions"]["T_half_GeV"], "this work"),
        ("nu_mu_decoupling", "nu_mu and nu_tau stop exchanging energy (Gamma = H)", td["nu_mu"], "this work"),
        ("nu_e_decoupling", "nu_e stops exchanging energy (Gamma = H)", td["nu_e"], "this work"),
        ("np_freeze_out", "Neutron-proton ratio freezes near 1/6 (exp(-Q/T) = 1/6)", m_n_minus_p / np.log(6),
         "PDG BBN review (n/p = 1/6), Q from PDG masses"),
        ("electron_annihilation", "Electrons and positrons annihilate (half gone)", halves["electron"]["T_half_GeV"], "this work"),
        ("nuclei_form", "Deuterium survives: nuclei start to form (eta^-1 exp(-B_D/T) = 1)",
         B_D / np.log(1 / eta), "PDG BBN review criterion, our eta"),
        ("equality", "Matter-radiation equality", float(T[i_eq]), "this work (Planck 2018 densities)"),
        ("last_scattering", "Atoms form: the CMB is released (z* = 1089.80)", C.T0_GEV * (1 + z_star), "Planck 2018"),
        ("today", "Today", C.T0_GEV, "Fixsen 2009"),
    ]
    out = []
    for key, label, Tq, src in events:
        out.append({"key": key, "label": label, "T_GeV": _f(Tq), "T_K": _f(Tq / C.K_B_GEV_PER_K),
                    "t_s": _f(ct.time_at(cal, Tq)), "source": src,
                    "hubble_radius_m": _f(ct.value_at(cal, "hubble_radius_m", Tq)),
                    "particle_horizon_m": _f(ct.value_at(cal, "particle_horizon_m", Tq)),
                    "M_hubble_Msun": _f(ct.value_at(cal, "m_hubble_msun", Tq)),
                    "f_today_Hz": _f(ct.value_at(cal, "f_today_hz", Tq))})
    # uncertainty of the QCD crossover time from T_pc +- 1.5 MeV
    Tq, dT = C.LIT["T_pc_hotqcd2019"]["value"], C.LIT["T_pc_hotqcd2019"]["err"]
    t_qcd_err = (ct.time_at(cal, Tq - dT) - ct.time_at(cal, Tq + dT)) / 2
    # lab and nature
    mean_E_over_T = lambda TT: rho_over_n(TT) / float(np.interp(np.log(TT), np.log(lat["T"]), _number_density(lat)))
    E_beam = 6.8e3  # GeV per LHC proton (Run 3)
    from scipy.optimize import brentq

    T_lhc = brentq(lambda TT: mean_E_over_T(TT) * TT - E_beam, 100.0, 1e4)
    labs = [
        ("lhc_proton", "An LHC proton (6.8 TeV) = the mean particle energy", float(T_lhc), "CERN, Run 3 (13.6 TeV collisions)"),
        ("sun_surface", "The Sun's surface, 5,772 K", 5772.0 * C.K_B_GEV_PER_K, "IAU 2015 Resolution B3 (nominal T_eff)"),
        ("water_boils", "Water boils: 373 K", 373.15 * C.K_B_GEV_PER_K, "Loeb 2014's 'habitable epoch' (273-373 K)"),
        ("water_freezes", "Water freezes: 273 K", 273.15 * C.K_B_GEV_PER_K, "Loeb 2014's 'habitable epoch' (273-373 K)"),
        ("alice_qgp", "Quark-gluon plasma at the LHC: ALICE photons, 297 ± 43 MeV", 0.297, "ALICE 2016, PLB 754, 235 (0-20%)"),
        ("phenix_qgp", "Quark-gluon plasma at RHIC: PHENIX photons, 221 ± 27 MeV", 0.221, "PHENIX 2010, PRL 104, 132301 (0-20%)"),
        ("sun_core", "The centre of the Sun, 15.6 million K", 1.560e7 * C.K_B_GEV_PER_K, "B16-GS98 standard solar model (Vinyoles+2017)"),
    ]
    lab_rows = [{"key": k, "label": lab, "T_GeV": _f(Tq), "t_s": _f(ct.time_at(cal, Tq)), "source": s,
                 "mean_E_over_T": _f(mean_E_over_T(Tq)) if k == "lhc_proton" else None} for k, lab, Tq, s in labs]
    # matter - dark energy equality (rho_m = rho_Lambda) and the start of acceleration (q = 0)
    Lu = cal["late"]
    a_mL = (Lu["rho_m0"] / Lu["rho_L"]) ** (1 / 3)
    a_acc = (Lu["rho_m0"] / (2 * Lu["rho_L"])) ** (1 / 3)
    T_mL, T_acc = C.T0_GEV / a_mL, C.T0_GEV / a_acc
    late = {"matter_lambda_equality": {"z": _f(1 / a_mL - 1, 4), "t_Gyr": _f(ct.time_at(cal, T_mL) / C.YEAR_S / 1e9, 4)},
            "acceleration_begins": {"z": _f(1 / a_acc - 1, 4), "t_Gyr": _f(ct.time_at(cal, T_acc) / C.YEAR_S / 1e9, 4)}}
    return {"events": out, "labs": lab_rows, "late_epochs": late, "eta": _f(eta, 5), "z_eq": _f(z_eq, 5),
            "t_qcd_err_s": _f(t_qcd_err, 3), "age_today_Gyr": _f(cal["t_s"][0] / C.YEAR_S / 1e9, 7),
            "late_universe": {k: _f(v) for k, v in cal["late"].items()},
            "byrnes2018_MH_formula_check": {
                "T_MeV": [1.0, 156.5, 200.0],
                "ours_Msun": [_f(ct.value_at(cal, "m_hubble_msun", x)) for x in (1e-3, 0.1565, 0.2)],
                "formula_Msun": [_f(1.5e5 * (np.interp(np.log(x), np.log(cal["T"]), cal["g_rho"]) / 10.75) ** -0.5
                                    * (x * 1e3) ** -2) for x in (1e-3, 0.1565, 0.2)]}}


def _number_density(E) -> np.ndarray:
    """Number density of the bath (GeV^3) for <E> = rho/n: every species as an ideal gas."""
    T = E["T"]
    n = np.zeros_like(T)
    for sp in particles.standard_model():
        if sp.group in ("neutrino",):
            n += thermo.gas(sp, T, T * E["T_nu_mu_over_T"])["n"]
        else:
            n += thermo.gas(sp, T)["n"]
    return n


def eos_softening(T, lat) -> dict:
    """The dips of w and c_s^2 (the universe getting 'softer') at the QCD crossover and elsewhere."""
    ss, lms = eos.saikawa_shirai(), eos.laine_meyer_schroeder()
    out = {}
    for name, (a, b) in {"qcd": (0.03, 1.0), "ew": (30.0, 1e3), "electron": (1e-5, 3e-3)}.items():
        m = (T > a) & (T < b)
        i = np.argmin(lat["w"][m])
        j = np.argmin(lat["cs2"][m])
        ms = (ss["w_T"] > a) & (ss["w_T"] < b)
        ml = (lms["T"] > a) & (lms["T"] < b)
        row = {"w_min_ours": _f(lat["w"][m][i], 4), "T_w_min_ours_GeV": _f(T[m][i]),
               "cs2_min_ours": _f(lat["cs2"][m][j], 4), "T_cs2_min_ours_GeV": _f(T[m][j])}
        if ms.any():  # S&S tabulate w only above 10 MeV
            k = np.argmin(ss["w"][ms])
            row.update({"w_min_SS": _f(ss["w"][ms][k], 4), "T_w_min_SS_GeV": _f(ss["w_T"][ms][k])})
        if ml.any():
            l_ = np.argmin(lms["w"][ml])
            row.update({"w_min_LMS": _f(lms["w"][ml][l_], 4), "T_w_min_LMS_GeV": _f(lms["T"][ml][l_])})
        out[name] = row
    return out


# --- run everything ---------------------------------------------------------------------------------


def run_core() -> dict:
    t0 = time.time()
    C.DERIVED_DIR.mkdir(parents=True, exist_ok=True)
    T = eos.grid(C.T_MIN, C.T_MAX)
    lat, ideal, pions, variants, band = equation_of_state(T)
    sm = particles.standard_model()
    halves = half_temperatures(sm)
    res = {"meta": {"T_range_GeV": [C.T_MIN, C.T_MAX], "grid_points": len(T),
                    "pdg_table": "particle package table " + ", ".join(__import__("particle").Particle.table_names())
                    + f" (particle {__import__('particle').__version__})",
                    "masses_GeV": {sp.name: _f(sp.mass) for sp in sm}}}
    res["checks"] = checks(T, lat, ideal)
    res["comparisons"] = comparisons(T, lat, ideal, band)
    res["half_temperatures"] = halves
    res["qcd"] = qcd_crossover(T, lat, ideal)
    res["eos_softening"] = eos_softening(T, lat)
    nu, runs = neutrino_decoupling()
    res["neutrinos"] = nu
    cal = ct.calendar(lat)
    res["calendar"] = calendar(lat, cal, halves)
    for name, row in res["eos_softening"].items():
        row["M_hubble_at_w_min_Msun"] = _f(ct.value_at(cal, "m_hubble_msun", row["T_w_min_ours_GeV"]))
        row["t_at_w_min_s"] = _f(ct.time_at(cal, row["T_w_min_ours_GeV"]))
    # temperatures whose horizon-entry modes are in the pulsar-timing band today (2-100 nHz)
    f = cal["f_today_hz"]
    res["calendar"]["pta_band_T_GeV"] = [_f(np.exp(np.interp(np.log(x), np.log(f), np.log(cal["T"])))) for x in (2e-9, 1e-7)]
    print(f"  calendar: age today {res['calendar']['age_today_Gyr']} Gyr")

    # curves for the figures and the explorer
    ss, lms, b16, esc = eos.saikawa_shirai(), eos.laine_meyer_schroeder(), qcd.borsanyi2016_spline_table(), eos.escudero2020()
    strong = qcd.strong_sector(T)
    base = runs["fdme_qed3_separate"]
    inst = runs["instantaneous_textbook"]
    np.savez_compressed(
        C.DERIVED_DIR / "curves.npz",
        T=T, g_rho=lat["g_rho"], g_s=lat["g_s"], w=lat["w"], cs2=lat["cs2"],
        g_rho_ideal=ideal["g_rho"], g_s_ideal=ideal["g_s"], g_rho_pions=pions["g_rho"],
        band_lo=band[0], band_hi=band[1], band_lo_s=band[2], band_hi_s=band[3],
        strong_p=strong["p"], strong_rho=strong["rho"], strong_s=strong["s"], strong_delta=strong["delta"],
        strong_delta_hrg=strong["delta_hrg"], strong_p_hrg=strong["p_hrg"],
        parts_names=np.array(list(lat["parts_g_rho"].keys())),
        parts_g_rho=np.array(list(lat["parts_g_rho"].values())),
        ideal_parts_names=np.array(list(ideal["parts_g_rho"].keys())),
        ideal_parts_g_rho=np.array(list(ideal["parts_g_rho"].values())),
        ss_T=ss["T"], ss_g_rho=ss["g_rho"], ss_g_rho_err=ss["g_rho_err"], ss_g_s=ss["g_s"], ss_g_s_err=ss["g_s_err"],
        ss_w_T=ss["w_T"], ss_w=ss["w"], ss_w_err=ss["w_err"],
        lms_T=lms["T"], lms_g_rho=lms["g_rho"], lms_g_s=lms["g_s"], lms_w=lms["w"], lms_cs2=lms["cs2"],
        b16_T=b16["T"], b16_g_rho=b16["g_rho"], b16_g_s=b16["g_s"],
        esc_T=esc["T"], esc_ratio=esc["Tg_over_Tnu"], esc_t=esc["t"], esc_g_rho=esc["g_rho"], esc_g_s=esc["g_s"],
        nu_T=base["T_gamma"], nu_Te=base["T_nu_e"], nu_Tm=base["T_nu_mu"], nu_t=base["t_s"],
        nu_inst_T=inst["T_gamma"], nu_inst_Tn=inst["T_nu_e"],
        cal_T=cal["T"], cal_t=cal["t_s"], cal_H=cal["H_s"], cal_a=cal["a"], cal_RH=cal["hubble_radius_m"],
        cal_dp=cal["particle_horizon_m"], cal_MH=cal["m_hubble_msun"], cal_f=cal["f_today_hz"],
        cal_g_rho=cal["g_rho"], cal_rho_bath=cal["rho_bath"], cal_rho_m=cal["rho_m"],
        **{f"var_{k}": v["g_rho"] for k, v in variants.items()},
    )
    res["meta"]["runtime_s"] = round(time.time() - t0, 1)
    C.RESULTS_PATH.write_text(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    r = run_core()
    print(json.dumps(r["checks"], indent=1))
