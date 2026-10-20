"""The headline numbers: a JSON summary (data/results.json) and a printed report."""

from __future__ import annotations

import numpy as np

import literature as L


def summarize(res) -> dict:
    fit = res["sirius"]
    hs = fit["hier"]["summary"]
    rep = fit["replication"]
    ms = res["masses"]
    adopted = [k for k in ms if k.startswith("Bond")][0]
    g = res["gaia"]
    wd = res["wd"]
    ml = res["mlr"]
    s = ml["sample"]
    return {
        "sirius_data": {"rows": fit["n_rows"], "used": fit["n_ours_used"], "by_method": fit["n_by_method"],
                        "flagged_outliers": fit["screen"]["n_flagged"],
                        "flagged_also_bond_r1": fit["screen"]["n_flagged_bond_r1"], "bond_r1": fit["screen"]["n_bond_r1"]},
        "sirius_blind_search": dict(zip(("P", "T", "e", "a", "i", "Omega", "omega"), map(float, fit["ti"]["theta"]))),
        "sirius_replication": {n: (float(v), float(e)) for n, v, e in
                               zip(("P", "T", "e", "a", "i", "Omega", "omega"), rep["theta"], rep["err"])},
        "sirius_orbit": {k: (v["median"], v["std"]) for k, v in hs.items()},
        "sirius_hyper": {k: (v["median"], v["std"]) for k, v in fit["hier"]["hyper_summary"].items()},
        "sirius_alternatives": {k: {n: (v["summary"][n]["median"], v["summary"][n]["std"]) for n in ("P", "a", "e")}
                                for k, v in fit["alternatives"].items()},
        "sirius_masses": {k: {m: (v[m]["median"], v[m]["std"]) for m in ("Mtot", "MA", "MB")}
                          for k, v in ms.items() if k != "samples" and isinstance(v, dict) and "Mtot" in v},
        "sirius_a3P2_systematics_percent": {k: 100 * v for k, v in ms["a3P2_spread"].items()},
        "sirius_gaia": {"plx_dr3": g["obs"]["plx"], "zp": g["zp"], "d_obs": g["d_obs"], "s_pm": g["s_pm"],
                        "transits": g["n_transits"], "visibility_periods": g["n_vp"],
                        "bias_all_transits": g["full_fit"], "bias_given_pm": g["bias_given_pm"],
                        "plx_corrected": g["plx_corrected"], "frac_reproduce_all_three": g["frac_reproduce"],
                        "mc_plx_bias_percentiles": np.percentile(g["mc"]["plx_bias"], [16, 50, 84]).tolist()},
        "sirius_b": {"MB": (wd["MB"]["median"], wd["MB"]["std"]),
                     "R": {k: (v["R"]["median"], v["R"]["std"]) for k, v in wd["radii"].items()},
                     "L": {k: (v["L"]["median"], v["L"]["std"]) for k, v in wd["radii"].items()},
                     "vgr": (wd["vgr"]["median"], wd["vgr"]["std"]), "logg_dyn": (wd["logg_dyn"]["median"],
                                                                                wd["logg_dyn"]["std"]),
                     "M_from_R_montreal": wd["M_from_R"], "cooling_age_myr": wd["cooling_age_myr"],
                     "chandrasekhar_limit": float(wd["chandra"]["M"].max())},
        "mlr": {"n_grade12": int(len(s)), "by_parallax_source": s["plx_source"].replace("", "none").value_counts().to_dict(),
                "n_fit": ml["n_fit"], "n_fit_gaia": ml["n_fit_gaia"], "sigma_int_dex": ml["sigma_int"],
                "f_out": ml["f_out"], "n_outliers": ml["n_outliers"],
                "relation_knots": {f"MV={k:g}": float(10 ** v) for k, v in
                                   zip(ml["like"].knots, ml["like"].nodes(ml["best"]))},
                "sun_mass_at_MV_4.81": float(10 ** np.interp(4.81, ml["MV"], ml["band"][2]))},
        "gallery": [{k: g_[k] for k in ("name", "plx", "plx_src", "a_au", "Mtot", "Mtot_err", "MA", "MB", "mass_src")}
                    | {"P": g_["elements"][0]} for g_ in res["gallery"]],
    }


def print_summary(res) -> None:
    s = summarize(res)
    o = s["sirius_orbit"]
    print("\n=== Sirius: our orbit (hierarchical model, MCMC) vs Bond et al. 2017 ===")
    for k in ("P", "a", "e", "i", "Omega", "omega", "T"):
        b = L.BOND17_ELEMENTS.get(k if k != "T" else "T_besselian")
        print(f"  {k:6s} {o[k][0]:12.5f} ± {o[k][1]:.5f}   (Bond: {b[0]} ± {b[1]})")
    print("  blind Thiele-Innes start: P = %.3f, e = %.4f, a = %.4f" % (s["sirius_blind_search"]["P"],
                                                                      s["sirius_blind_search"]["e"],
                                                                      s["sirius_blind_search"]["a"]))
    print("  replication with Bond's weights: P = %.4f ± %.4f (Julian yr), a = %.5f ± %.5f" % (
        *s["sirius_replication"]["P"], *s["sirius_replication"]["a"]))
    for k, v in s["sirius_alternatives"].items():
        print(f"  {k:20s} P = {v['P'][0]:.3f} ± {v['P'][1]:.3f}   a = {v['a'][0]:.4f} ± {v['a'][1]:.4f}")
    print("\n=== Sirius masses (Msun) ===")
    for k, v in s["sirius_masses"].items():
        print(f"  {k:48s} Mtot {v['Mtot'][0]:.3f} ± {v['Mtot'][1]:.3f}   A {v['MA'][0]:.3f} ± {v['MA'][1]:.3f}"
              f"   B {v['MB'][0]:.3f} ± {v['MB'][1]:.3f}")
    print("  (Bond et al. 2017: Mtot 3.081 ± 0.034, A 2.063 ± 0.023, B 1.018 ± 0.011)")
    g = s["sirius_gaia"]
    print("\n=== Sirius B in Gaia DR3 ===")
    print(f"  DR3 parallax {g['plx_dr3']:.2f} mas (zp {g['zp']:+.3f}); minus Hipparcos: {g['d_obs']['plx']:+.2f} mas;"
          f" proper motion minus expectation ({g['d_obs']['pmra']:+.2f}, {g['d_obs']['pmdec']:+.2f}) mas/yr")
    print(f"  simulated 5-parameter parallax bias, 16/50/84%: {np.round(g['mc_plx_bias_percentiles'], 2)} mas;"
          f" given the proper motion: {g['bias_given_pm'][0]:+.2f} ± {g['bias_given_pm'][1]:.2f}")
    print(f"  Gaia DR3 corrected: {g['plx_corrected'][0]:.2f} ± {g['plx_corrected'][1]:.2f} mas "
          f"(Hipparcos 379.21 ± 1.58; adopted 378.9 ± 1.4)")
    b = s["sirius_b"]
    print("\n=== Sirius B, white dwarf ===")
    for k, v in b["R"].items():
        print(f"  R from {k:36s} {v[0]:.6f} ± {v[1]:.6f} Rsun   L = {b['L'][k][0]:.4f} ± {b['L'][k][1]:.4f} Lsun")
    print(f"  gravitational redshift {b['vgr'][0]:.1f} ± {b['vgr'][1]:.1f} km/s (Joyce et al. 2018: 80.65 ± 0.77)")
    print(f"  log g from M and R {b['logg_dyn'][0]:.3f} ± {b['logg_dyn'][1]:.3f}; Montreal mass for our radius "
          f"{b['M_from_R_montreal']:.3f}; cooling age {b['cooling_age_myr']:.0f} Myr; Chandrasekhar limit "
          f"{b['chandrasekhar_limit']:.3f}")
    m = s["mlr"]
    print("\n=== ORB6 mass-luminosity relation ===")
    print(f"  grade 1-2 orbits {m['n_grade12']}; parallax sources {m['by_parallax_source']}")
    print(f"  fit: {m['n_fit']} systems ({m['n_fit_gaia']} with Gaia parallaxes); intrinsic scatter "
          f"{m['sigma_int_dex']:.3f} dex; outlier fraction {m['f_out']:.3f} ({m['n_outliers']} systems)")
    print("  mass at knots: " + ", ".join(f"{k} {v:.2f}" for k, v in m["relation_knots"].items()))
    print(f"  a star like the Sun (M_V = 4.81): {m['sun_mass_at_MV_4.81']:.3f} Msun")
    print("\n=== Gallery ===")
    for gg in s["gallery"]:
        print(f"  {gg['name']:16s} P {gg['P']:7.2f} yr  a {gg['a_au']:6.1f} au  M {gg['Mtot']:.2f} ± {gg['Mtot_err']:.2f}"
              f"  ({gg['plx_src']})")
