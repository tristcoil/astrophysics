"""Every result of the project, computed step by step and cached.

Each step writes ``data/derived/<step>.json``; ``run(force=True)`` recomputes
everything. ``data/results.json`` collects the headline numbers.
"""

from __future__ import annotations

import json
import math
import time

import numpy as np

import config as cfg
import datasets as D
import figparse
import hadrons
import hrg
import lattice
import photons
import thermal_fit as F
import universe as U

NUCLEI = ("d", "dbar", "He3", "He3bar", "H3L", "H3Lbar", "He4", "He4bar")


def _cached(name: str, fn, force: bool):
    path = cfg.DERIVED_DIR / f"{name}.json"
    if path.exists() and not force:
        return json.loads(path.read_text())
    t = time.time()
    out = fn()
    path.write_text(json.dumps(out, indent=1, default=_json_default))
    print(f"  {name}: {time.time() - t:.0f} s")
    return json.loads(path.read_text())


def _json_default(o):
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def _brief(r: dict) -> dict:
    p = r["params"]
    out = {"T": p["T"]["value"], "T_err": p["T"]["error"], "muB": p["muB"]["value"],
           "muB_err": p["muB"]["error"], "V": p["V"]["value"], "V_err": p["V"]["error"],
           "V_EV": r["V_excluded_volume"], "chi2": r["chi2"], "ndf": r["ndf"]}
    if "gammaS" in p:
        out.update(gammaS=p["gammaS"]["value"], gammaS_err=p["gammaS"]["error"])
    if "muS" in p:
        out.update(muS=p["muS"]["value"])
    return out


# ---------------------------------------------------------------------------
# Steps
# ---------------------------------------------------------------------------


def step_compilation() -> dict:
    """Our 0-10% compilation against Thermal-FIST's file and Andronic's Fig. 1."""
    ds = D.alice("0-10", include_lambdabar=True)
    tf = D.thermal_fist_compilation()
    fig = figparse.andronic_fig1()
    rows = []
    for y in ds.yields:
        pdg = 310 if y.key == "K0S" else next(iter(y.targets))
        f = fig["entries"].get(y.key, {})
        rows.append({"key": y.key, "html": y.html, "value": y.value, "stat": y.stat,
                     "sys": y.sys, "err": y.err, "fitted": y.fitted, "source": y.source,
                     "note": y.note, "tf_value": tf.get(pdg, (None, None))[0],
                     "tf_err": tf.get(pdg, (None, None))[1], "fig_value": f.get("data"),
                     "fig_model": f.get("model_dot")})
    return {"rows": rows, "fig_axis_residual_decades": fig["axis_fit_residual_decades"]}


def step_main_fit() -> dict:
    ds = D.alice("0-10")
    r = F.fit(ds, F.FitSettings(), minos=True)
    gas = hrg.model()
    st = gas.solve_mu(r["params"]["T"]["value"] / 1e3, r["params"]["muB"]["value"] / 1e3,
                      ds.q_over_b)
    # composition of the gas at freeze-out (primordial)
    hl = gas.hl
    groups = {"pions": np.isin(np.abs(hl.pdg), [211, 111]),
              "kaons": np.isin(np.abs(hl.pdg), [321, 311]),
              "nucleons": np.isin(np.abs(hl.pdg), [2212, 2112]),
              "other ground states": np.isin(np.abs(hl.pdg), [221, 3122, 3212, 3222, 3112,
                                                               3322, 3312, 3334]),
              }
    res = ~np.any(np.array(list(groups.values())), axis=0) & ~hl.is_nucleus
    groups["resonances"] = res
    n_tot = st.n.sum()
    comp = {k: float(st.n[m].sum() / n_tot) for k, m in groups.items()}
    r["composition"] = comp
    r["n_species"] = int(len(hl))
    r["n_with_decays"] = int((~hl.stable).sum())
    r["charm_per_pair"] = hadrons.charm_feeddown_per_pair()
    return r


def step_profile() -> dict:
    return F.profile_T(D.alice("0-10"), F.FitSettings(), np.arange(140.0, 172.01, 0.5))


def step_variants() -> dict:
    ds = D.alice("0-10")
    dsl = D.alice("0-10", include_lambdabar=True)
    da = D.alice_andronic_compilation()
    no_p = tuple(y.key for y in ds.yields if y.fitted and y.key not in ("p", "pbar"))
    hadr = tuple(y.key for y in ds.yields if y.fitted and y.key not in NUCLEI)
    S = F.FitSettings
    runs = [
        ("default", "Default: PDG 2025 list, Breit-Wigner widths, quantum statistics, "
         "charm feed-down (10 pairs)", ds, S(), None, "model"),
        ("no_charm", "No charm feed-down", ds, S(charm_pairs=0.0), None, "model"),
        ("charm5", "Charm feed-down 5 pairs", ds, S(charm_pairs=5.0), None, "model"),
        ("charm15", "Charm feed-down 15 pairs", ds, S(charm_pairs=15.0), None, "model"),
        ("zero_width", "Zero-width resonances", ds, S(widths="zero"), None, "model"),
        ("pdg2014", "PDG 2014 hadron list", ds, S(list_name="PDG2014"), None, "model"),
        ("boltzmann", "Boltzmann statistics", ds, S(stats="boltzmann"), None, "model"),
        ("gammaS", "Strangeness saturation gamma_S free", ds, S(free_gammaS=True), None, "freedom"),
        ("lambdabar", "Lambda-bar added (as Andronic et al.)", dsl, S(), None, "data"),
        ("andronic_data", "Andronic et al.'s data compilation", da, S(), None, "data"),
        ("no_protons", "Without p and p-bar", ds, S(), no_p, "subset"),
        ("hadrons_only", "Hadrons only (no nuclei)", ds, S(), hadr, "subset"),
        ("nuclei_only", "Nuclei only (d, 3He, hypertriton, 4He and antinuclei)", ds, S(), NUCLEI,
         "subset"),
        ("d_he3_he4", "d, 3He, 4He and antinuclei, mu_B = 0 (as ALICE's 4He paper)", ds,
         S(fix_muB=0.0), ("d", "dbar", "He3", "He3bar", "He4", "He4bar"), "subset"),
        ("p_d_he3", "p, d, 3He and antiparticles, mu_B = 0 (as ALICE's nuclei paper)", ds,
         S(fix_muB=0.0), ("p", "pbar", "d", "dbar", "He3", "He3bar"), "subset"),
    ]
    out = {}
    for key, label, d, s, keys, kind in runs:
        r = F.fit(d, s, minos=False, keys=keys)
        out[key] = {"label": label, "kind": kind, **_brief(r),
                    "pulls": {y["key"]: y["pull"] for y in r["yields"] if y["fitted"]}}
    T0 = out["default"]["T"]
    model = [k for k, v in out.items() if v["kind"] == "model" and k != "default"]
    pairs = [("no_charm", None), ("charm5", "charm15"), ("zero_width", None),
             ("pdg2014", None), ("boltzmann", None)]
    contrib = {}
    for a, b in pairs:
        if b is None:
            contrib[a] = abs(out[a]["T"] - T0)
        else:
            contrib[a + "/" + b] = 0.5 * abs(out[a]["T"] - out[b]["T"])
    # the charm on/off difference is covered by the +-5-pair range: use that
    contrib.pop("no_charm")
    out["_model_systematic_T"] = {"components": contrib,
                                  "total": math.sqrt(sum(v * v for v in contrib.values())),
                                  "variants": model}
    return out


def step_andronic_comparison() -> dict:
    """Our model at Andronic et al.'s best fit vs their model, particle by particle."""
    fig = figparse.andronic_fig1()["entries"]
    ds = D.alice("0-10", include_lambdabar=True)
    A = cfg.ANDRONIC_2018
    out = {"their_params": {"T": A["T"], "muB": A["muB"], "V": A["V"]}, "variants": {}}
    keys = list(fig)
    theirs = np.array([fig[k]["model_dot"] for k in keys])
    nuc = [i for i, k in enumerate(keys) if k in NUCLEI]
    for label, s in (("default", F.FitSettings()), ("no charm", F.FitSettings(charm_pairs=0.0)),
                     ("zero width", F.FitSettings(widths="zero")),
                     ("PDG2014", F.FitSettings(list_name="PDG2014"))):
        m = F.YieldModel(ds, s)
        dens, st = m.densities(A["T"] / 1e3, A["muB"] / 1e3)
        idx = [[y.key for y in ds.yields].index(k) for k in keys]
        # volume from the nuclei (pure thermal: no decays, no charm) -> same V in both codes
        V = float(np.exp(np.mean(np.log(theirs[nuc] / dens[idx][nuc]))))
        ours = V * dens[idx] + m.c[idx]
        out["variants"][label] = {"V_matched": V,
                                  "ratio": {k: float(o / t) for k, o, t in zip(keys, ours, theirs)}}
    da = D.alice_andronic_compilation()
    out["fit_their_data"] = _brief(F.fit(da, F.FitSettings(), minos=False))
    out["fit_their_data_no_charm"] = _brief(F.fit(da, F.FitSettings(charm_pairs=0.0), minos=False))
    out["fig1"] = fig
    return out


def step_centrality() -> dict:
    npart = D.alice_npart()
    out = []
    for c in D.ALICE_CENTRALITIES:
        ds = D.alice(c)
        rec = {"centrality": c, "npart": npart[D._cent(c)], "charm_pairs": ds.charm_pairs}
        for lab, s in (("eq", F.FitSettings()), ("gammaS", F.FitSettings(free_gammaS=True))):
            r = F.fit(ds, s, minos=False)
            rec[lab] = _brief(r)
            rec[lab]["ratios"] = {y["key"]: y["ratio"] for y in r["yields"]}
            rec[lab]["errs"] = {y["key"]: y["err"] / y["model"] for y in r["yields"]}
        out.append(rec)
    return {"classes": out}


def step_star() -> dict:
    pub = D.star_published_freezeout()
    out = {"published": {str(k): v for k, v in pub.items()}, "replica": {}, "ours": {}}
    for e in (7.7, 39.0):
        ds = D.star_gce_dataset(e)
        _, model = D.star_gce_inputs(e)
        r = F.fit(ds, F.FitSettings(free_muS=True, free_gammaS=True), minos=False)
        names = {"pi+": "pi+", "pi-": "pi-", "K+": "ka+", "K-": "ka-", "p": "p", "pbar": "pbar",
                 "Lambda": "lambda", "Lambdabar": "lambdabar", "Xi-": "Xi", "Xibar+": "Xibar"}
        out["replica"][str(e)] = {**_brief(r), "model_ratio_to_THERMUS": {
            y["key"]: y["model"] / model[names[y["key"]]] for y in r["yields"]}}
    for e in D.STAR_ENERGIES:
        ds = D.star(e)
        rec = {}
        for lab, s in (("eq", F.FitSettings()), ("gammaS", F.FitSettings(free_gammaS=True))):
            r = F.fit(ds, s, minos=False)
            rec[lab] = _brief(r)
            rec[lab]["yields"] = [{k: y[k] for k in ("key", "html", "value", "err", "model",
                                                       "ratio", "pull")} for y in r["yields"]]
            rec[lab]["muS_fit"], rec[lab]["muQ_fit"] = r["muS"], r["muQ"]
        out["ours"][str(e)] = rec
    return out


def step_phase_diagram(main: dict, star: dict) -> dict:
    fig = figparse.andronic_fig4()
    pub = D.star_published_freezeout()
    groups = {"rgb(100%, 0%, 0%)|5": "andronic", "rgb(0%, 0%, 0%)|6": "star_fig",
              "rgb(0%, 0%, 100%)|3": "other_a", "rgb(100%, 0%, 100%)|4": "other_b",
              "rgb(0%, 0%, 0%)|13": "other_c"}
    pts = {}
    for p in fig["points"]:
        g = groups.get(p["group"])
        if g is None:            # a lone marker at (21 MeV, 100 MeV): a legend symbol
            continue
        T_err = max(p["T_hi"] - p["T"], p["T"] - p["T_lo"])
        pts.setdefault(g, []).append({"muB": p["muB"], "T": p["T"], "T_err": T_err,
                                      "muB_lo": p["muB_lo"], "muB_hi": p["muB_hi"]})
    # STAR's points in the figure vs STAR's own HEPData values
    match = []
    for q in pts.get("star_fig", []):
        e = min(pub, key=lambda k: abs(pub[k]["muB"] - q["muB"]))
        match.append({"energy": e, "dT": q["T"] - pub[e]["T"], "dmuB": q["muB"] - pub[e]["muB"]})
    mu = np.linspace(0, 300, 61)   # HotQCD: the expansion holds for mu_B <~ 2 T_pc
    lines = {s: {k: v.tolist() for k, v in lattice.tpc_line(mu, s).items()}
             for s in ("hotqcd", "borsanyi")}
    ours = [{"label": "ALICE 2.76 TeV", "sqrt_s": 2760.0, "T": main["params"]["T"]["value"],
             "T_err": main["params"]["T"]["error"], "muB": main["params"]["muB"]["value"],
             "muB_err": main["params"]["muB"]["error"]}]
    for e, rec in star["ours"].items():
        r = rec["eq"]
        ours.append({"label": f"STAR {float(e):g} GeV", "sqrt_s": float(e), "T": r["T"],
                     "T_err": r["T_err"], "muB": r["muB"], "muB_err": r["muB_err"]})
    T_path = np.concatenate([np.linspace(250, 150, 21), np.linspace(149, 60, 30),
                             np.linspace(58, 5, 54)])
    path = U.path(T_path)
    # chemical freeze-out at fixed energy per hadron, through our LHC point
    gas = hrg.model()
    st = gas.solve_mu(main["params"]["T"]["value"] / 1e3, main["params"]["muB"]["value"] / 1e3)
    tot = gas.totals(st)
    e_over_n = tot["e"] / tot["n"]
    from scipy.optimize import brentq
    en_line = []
    for m in np.linspace(0, 800, 33):
        def f(T):
            s2 = gas.solve_mu(T, m / 1e3, cfg.Q_OVER_B["Au"])
            t2 = gas.totals(s2)
            return t2["e"] / t2["n"] - e_over_n
        try:
            en_line.append((float(m), float(brentq(f, 0.03, 0.2, xtol=1e-7) * 1e3)))
        except ValueError:
            break
    return {"published": pts, "star_hepdata": {str(k): v for k, v in pub.items()},
            "star_figure_vs_hepdata": match, "lattice": lines, "ours": ours,
            "universe_path": {k: v.tolist() for k, v in path.items()},
            "e_over_n_GeV": e_over_n, "e_over_n_line": en_line,
            "nuclear_matter_muB": 923.0,
            "axis_residuals": {"x": fig["x_axis_residual"], "y_MeV": fig["y_axis_residual_MeV"]}}


def step_lattice() -> dict:
    gas = hrg.model()
    Ts = np.arange(100.0, 200.01, 2.5)
    hrg_rows = []
    for T in Ts:
        t = T / 1e3
        tot = gas.totals(gas.state(t))
        T4 = t ** 4 / cfg.HBARC ** 3
        chi = gas.susceptibilities(t)
        hrg_rows.append({"T": T, "p": tot["p"] / T4, "e": tot["e"] / T4,
                         "s": tot["s"] * cfg.HBARC ** 3 / t ** 3,
                         "I": (tot["e"] - 3 * tot["p"]) / T4, "chiBB": chi["BB"],
                         "chiQQ": chi["QQ"], "chiSS": chi["SS"]})
    hq = {k: v.tolist() for k, v in lattice.hotqcd_2014().items()}
    wb = {k: v.tolist() for k, v in lattice.wb_2014().items()}
    sus = {k: v.tolist() for k, v in lattice.hotqcd_2012_susceptibilities().items()}
    Tf = np.arange(100.0, 400.01, 1.0)
    ans = {k: v.tolist() for k, v in lattice.hotqcd_ansatz_thermo(Tf).items()}
    ans["T"] = Tf.tolist()
    # where does the hadron gas leave the lattice band (HotQCD 2014, pressure and s)?
    h_T = np.array([r["T"] for r in hrg_rows])
    first = {}
    for q in ("p", "e", "s"):
        hv = np.interp(hq["T"], h_T, [r[q] for r in hrg_rows])
        dev = (hv - np.array(hq[q])) / np.array(hq[q + "_err"])
        bad = [T for T, d in zip(hq["T"], dev) if T <= 200 and abs(d) > 1]
        first[q] = min(bad) if bad else None
    st = gas.state(cfg.HOTQCD_2019["Tpc"] / 1e3)
    tot = gas.totals(st)
    return {"hrg": hrg_rows, "hotqcd2014": hq, "wb2014": wb, "hotqcd2012_chi": sus,
            "hotqcd_ansatz": ans, "sb_p": lattice.SB_P, "hrg_leaves_hotqcd_band_at_MeV": first,
            "hrg_at_Tpc": {"e_GeV_fm3": tot["e"], "s_fm3": tot["s"], "n_fm3": tot["n"]},
            "ansatz_params": lattice.hotqcd_2014_ansatz_params()}


def step_universe(main: dict, phot: dict) -> dict:
    T_fit = main["params"]["T"]["value"]
    T_err = main["params"]["T"]["error"]
    cal = U.p2_calendar()
    temps = {"T_CF (our fit)": T_fit, "T_pc (HotQCD 2019)": cfg.HOTQCD_2019["Tpc"],
             "T_eff photons (ours)": phot["T_eff_MeV"]}
    ages = {}
    for k, T in temps.items():
        ages[k] = {"T_MeV": T, "age_s_P2": float(U.age_at(T)[0]) if cal else None,
                   "age_s_own": float(U.friedmann_age(np.array([T]))[0])}
    a_lo, a_hi = U.age_at(T_fit + T_err)[0], U.age_at(T_fit - T_err)[0]
    ages["T_CF (our fit)"]["age_err_s"] = float(0.5 * (a_hi - a_lo))
    mu_tpc = {l: U.chemical_potentials(cfg.HOTQCD_2019["Tpc"] / 1e3, l)
              for l in ("zero", "sphaleron")}
    # antimatter annihilation: where mu_B/T reaches 1 (antibaryons e^-2 below baryons)
    from scipy.optimize import brentq
    T_ann = brentq(lambda T: U.mu_b_nonlinear(T) / T - 1.0, 0.02, 0.06, xtol=1e-6)
    lhc_x = main["params"]["muB"]["value"] / T_fit
    return {"ages": ages, "P2_available": cal is not None, "eta": U.eta_from_omega_b(),
            "baryon_per_entropy": U.baryon_per_entropy(), "mu_at_Tpc": mu_tpc,
            "T_antibaryons_gone_MeV": T_ann * 1e3,
            "age_antibaryons_gone_s": float(U.age_at(T_ann * 1e3)[0]),
            "lhc_muB_over_T": lhc_x,
            "lhc_muB_over_T_err": main["params"]["muB"]["error"] / T_fit,
            "T_fit_kelvin": T_fit * cfg.MEV_TO_K}


def run(force: bool = False) -> dict:
    print("Analysis (cached steps reused)")
    comp = _cached("compilation", step_compilation, force)
    main = _cached("main_fit", step_main_fit, force)
    prof = _cached("profile", step_profile, force)
    var = _cached("variants", step_variants, force)
    andr = _cached("andronic_comparison", step_andronic_comparison, force)
    cent = _cached("centrality", step_centrality, force)
    star = _cached("star", step_star, force)
    phot = _cached("photons", photons.inverse_slope, force)
    lat = _cached("lattice", step_lattice, force)
    phase = _cached("phase_diagram", lambda: step_phase_diagram(main, star), force)
    univ = _cached("universe", lambda: step_universe(main, phot), force)
    p = main["params"]
    results = {
        "freeze_out": {"T_MeV": p["T"]["value"], "T_err": p["T"]["error"],
                       "T_minos": [p["T"].get("minos_lo"), p["T"].get("minos_hi")],
                       "T_model_sys": var["_model_systematic_T"]["total"],
                       "muB_MeV": p["muB"]["value"], "muB_err": p["muB"]["error"],
                       "V_fm3": p["V"]["value"], "V_err": p["V"]["error"],
                       "V_excluded_volume_fm3": main["V_excluded_volume"],
                       "V_excluded_volume_err": main["V_excluded_volume_err"],
                       "chi2": main["chi2"], "ndf": main["ndf"],
                       "proton_pull": [y["pull"] for y in main["yields"] if y["key"] == "p"][0]},
        "andronic_2018": cfg.ANDRONIC_2018,
        "fit_with_their_data": andr["fit_their_data"],
        "nuclei_only": var["nuclei_only"], "hotqcd_2019": cfg.HOTQCD_2019,
        "photons": {k: phot[k] for k in ("T_eff_MeV", "stat_MeV", "syst_MeV", "chi2", "ndf")},
        "universe": {"age_at_T_CF_us": univ["ages"]["T_CF (our fit)"]["age_s_P2"] * 1e6
                     if univ["P2_available"] else univ["ages"]["T_CF (our fit)"]["age_s_own"] * 1e6,
                     "muB_over_T_at_Tpc": univ["mu_at_Tpc"]["zero"]["muB_over_T"],
                     "T_antibaryons_gone_MeV": univ["T_antibaryons_gone_MeV"]},
        "star_replica": star["replica"],
        "n_species": main["n_species"],
    }
    cfg.RESULTS_PATH.write_text(json.dumps(results, indent=1, default=_json_default))
    return {"compilation": comp, "main": main, "profile": prof, "variants": var,
            "andronic": andr, "centrality": cent, "star": star, "photons": phot,
            "lattice": lat, "phase": phase, "universe": univ, "results": results}


if __name__ == "__main__":
    import sys
    run(force="--force" in sys.argv)
