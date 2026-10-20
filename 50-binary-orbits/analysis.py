"""Every number in the project, computed from the cached data (run.py calls ``run``)."""

from __future__ import annotations

import numpy as np
import pandas as pd

import config as C
import gaia_sirius as GS
import literature as L
import mlr
import orb6
import orbits
import sirius
import whitedwarf as WD

RNG_SEED = C.SEED


# ---------------------------------------------------------------------------
# Sirius: orbit, masses, Gaia, white dwarf
# ---------------------------------------------------------------------------


def sirius_masses(fit: dict) -> dict:
    """Masses from our orbit posterior with each parallax choice, and the error budget."""
    rng = np.random.default_rng(RNG_SEED)
    el = fit["hier"]["elements"]
    P, a = el[:, 0], el[:, 3]
    choices = {
        "Bond et al. 2017 adopted (Hipparcos + ground)": L.BOND17_PLX["adopted"],
        "Hipparcos (van Leeuwen 2007)": L.BOND17_PLX["HIP2"],
        "Ground-based (Gatewood & Gatewood 1978)": L.BOND17_PLX["GG78"],
    }
    out = {}
    for lab, plx in choices.items():
        m = sirius.masses_from_samples(P, a, rng, plx=plx, aA=L.BOND17_AA)
        out[lab] = {k: summary(v) for k, v in m.items()}
        out[lab]["plx"] = plx
        if lab.startswith("Bond"):
            out["samples"] = m
    # systematic spread of the orbit's a^3/P^2 across data choices
    k_all = np.median(a**3 / P**2)
    spread = {lab: float(np.median(v["draws"][:, 3] ** 3 / v["draws"][:, 0] ** 2) / k_all - 1)
              for lab, v in fit["alternatives"].items()}
    rep = fit["replication"]["theta"]
    spread["Bond et al. weights (replication)"] = float(rep[3] ** 3 / rep[0] ** 2 / k_all - 1)
    out["a3P2_spread"] = spread
    # linear error budget (as Bond et al.'s Table 7)
    cov = np.cov(np.vstack([P, a]))
    out["budget"] = sirius.error_budget((np.median(P), np.median(a)), cov, L.BOND17_PLX["adopted"], L.BOND17_AA)
    return out


def summary(x):
    x = np.asarray(x)
    q = np.percentile(x, [15.865, 50, 84.135])
    return {"median": float(q[1]), "lo": float(q[1] - q[0]), "hi": float(q[2] - q[1]), "std": float(np.std(x))}


def sirius_gaia(fit: dict, n_draw=20000) -> dict:
    """Sirius B in Gaia DR3 against what Hipparcos + our orbit predict; the curvature bias."""
    theta = np.median(fit["hier"]["elements"], axis=0)
    aA = L.BOND17_AA["value"]
    rec = pd.read_csv(C.GAIA_DIR / "sirius_b_dr3.csv").iloc[0]
    from zero_point import zpt
    zpt.load_tables()
    # (the package needs arrays under numpy 2)
    zp = float(zpt.get_zpt(np.array([rec["phot_g_mean_mag"]]), np.array([rec["nu_eff_used_in_astrometry"]]),
                           np.array([np.nan]), np.array([rec["ecl_lat"]]),
                           np.array([int(rec["astrometric_params_solved"])]))[0])
    plx_ref = L.BOND17_PLX["HIP2"]["value"] * 1e3
    exp = GS.expected_dr3(theta, aA, plx_ref)
    obs = {"plx": rec["parallax"] - zp, "pmra": rec["pmra"], "pmdec": rec["pmdec"],
           "e_plx": rec["parallax_error"], "e_pmra": rec["pmra_error"], "e_pmdec": rec["pmdec_error"]}
    leak = GS.hipparcos_leak_mc(theta, aA, n_draw=2000)
    d_obs = {"plx": obs["plx"] - plx_ref, "pmra": obs["pmra"] - exp["pmra"], "pmdec": obs["pmdec"] - exp["pmdec"]}
    # uncertainty of the *expected* proper motion: Hipparcos errors + leak scatter
    s_pmra = np.sqrt(L.HIP2_SIRIUS["pmra"][1] ** 2 + leak["pmra_leak"].std() ** 2 + obs["e_pmra"] ** 2)
    s_pmdec = np.sqrt(L.HIP2_SIRIUS["pmdec"][1] ** 2 + leak["pmdec_leak"].std() ** 2 + obs["e_pmdec"] ** 2)
    tr = pd.read_csv(C.CACHE_DIR / "sirius_b_transits_dr3.csv") if (C.CACHE_DIR / "sirius_b_transits_dr3.csv").exists() \
        else _transits()
    full = GS.simulate_fit(tr, rec["ra"], rec["dec"], theta, 1 - aA / theta[3], plx_ref)
    mc_path = C.CACHE_DIR / "sirius_b_subset_mc.csv"
    if mc_path.exists():
        mc = pd.read_csv(mc_path)
    else:
        mc = GS.subset_monte_carlo(tr, rec["ra"], rec["dec"], theta, 1 - aA / theta[3], plx_ref, n_draw=n_draw)
        mc.to_csv(mc_path, index=False)
    # condition on the proper-motion offsets only (they do not involve any parallax)
    w = np.exp(-0.5 * ((mc["pmra_bias"] - d_obs["pmra"]) / s_pmra) ** 2
               - 0.5 * ((mc["pmdec_bias"] - d_obs["pmdec"]) / s_pmdec) ** 2)
    w = w / w.sum()
    mean_b = float(np.sum(w * mc["plx_bias"]))
    sd_b = float(np.sqrt(np.sum(w * (mc["plx_bias"] - mean_b) ** 2)))
    n_eff = float(1 / np.sum(w**2))
    corrected = obs["plx"] - mean_b
    corrected_err = float(np.hypot(sd_b, obs["e_plx"]))
    return {"record": rec.to_dict(), "zp": zp, "obs": obs, "expected": exp, "d_obs": d_obs,
            "s_pm": (float(s_pmra), float(s_pmdec)), "leak": leak.describe().to_dict(),
            "transits": tr, "n_transits": len(tr), "n_vp": int(GS.visibility_periods(tr).max() + 1),
            "full_fit": {k: full[k] for k in ("plx_bias", "pmra_bias", "pmdec_bias", "rms_res", "n")},
            "mc": mc, "weights": w.to_numpy(), "bias_given_pm": (mean_b, sd_b), "n_eff": n_eff,
            "plx_corrected": (float(corrected), corrected_err),
            "frac_reproduce": float(np.mean(((mc["pmra_bias"] - d_obs["pmra"]) / s_pmra) ** 2
                                            + ((mc["pmdec_bias"] - d_obs["pmdec"]) / s_pmdec) ** 2
                                            + ((mc["plx_bias"] - d_obs["plx"]) / L.HIP2_SIRIUS["plx"][1]) ** 2 < 7.81))}


def _transits():
    rec = pd.read_csv(C.GAIA_DIR / "sirius_b_dr3.csv").iloc[0]
    tr = GS.transits(rec["ra"], rec["dec"])
    gaps = pd.read_csv(GS.GAPS, encoding="utf-8-sig")
    tr = GS.remove_gaps(tr, gaps)
    tr.to_csv(C.CACHE_DIR / "sirius_b_transits_dr3.csv", index=False)
    return tr


def sirius_b_white_dwarf(masses: dict) -> dict:
    """Radius from brightness + temperature, luminosity, redshift, mass-radius comparison."""
    MB = masses["samples"]["MB"]
    plx = L.BOND17_PLX["adopted"]
    teff, teff_err = L.BOND17_WD["Teff"]
    logg = L.BOND17_WD["logg"][0]
    radii = {}
    for lab, (mag, err, band) in {
        "V (HST/STIS, Barstow et al. 2005)": (*L.BARSTOW05["V_hst"], "V"),
        "V (ground, Holberg et al. 1984)": (*L.BARSTOW05["V_holberg84"], "V"),
        "Gaia G (DR3)": (8.524133, 0.03, "G3"),
    }.items():
        r = WD.photometric_radius(mag, err, plx["value"], plx["err"], teff, teff_err, logg, band=band)
        radii[lab] = {"R": summary(r["R"]), "L": summary(r["L"]), "samples": r["R"], "Lsamples": r["L"]}
    main = radii["V (HST/STIS, Barstow et al. 2005)"]
    n = min(len(MB), len(main["samples"]))
    vgr = WD.grav_redshift(MB[:n], main["samples"][:n])
    logg_dyn = np.log10(C.G_NEWTON * MB[:n] * C.M_SUN / (main["samples"][:n] * C.R_SUN) ** 2 * 100)
    mr = {T: WD.montreal_mr(T) for T in (10000.0, teff, 40000.0)}
    # the mass a Montreal C/O star of our radius would have at Sirius B's temperature
    m_t = mr[teff].sort_values("R")
    M_from_R = np.interp(main["R"]["median"], m_t["R"], m_t["M"])
    # cooling age at our mass
    age = np.interp(np.median(MB), mr[teff]["M"], mr[teff]["age"])
    return {"radii": radii, "vgr": summary(vgr), "logg_dyn": summary(logg_dyn), "MB": summary(MB),
            "montreal": mr, "M_from_R": float(M_from_R), "cooling_age_myr": float(age / 1e6),
            "chandra": WD.chandrasekhar_curve(2.0), "chandra_fe": WD.chandrasekhar_curve(56.0 / 26.0),
            "parsons": WD.parsons_table(), "visual_wds": WD.visual_binary_wds()}


# ---------------------------------------------------------------------------
# ORB6: gallery and the mass-luminosity relation
# ---------------------------------------------------------------------------


def gallery(sir_fit: dict, mlr_fit: dict) -> list:
    """The six famous orbits: elements, adopted distance, masses and their sources."""
    o = orb6.read_orb6()
    xm = pd.read_csv(C.GAIA_DIR / "orb6_xmatch_gaiadr3.csv")
    hip2 = mlr.read_hip2().set_index("hip")
    like, chain, best = mlr_fit["like"], mlr_fit["chain"], mlr_fit["best"]
    out = []
    for wds, disc, name in C.GALLERY:
        r = orb6.find(o, wds, disc)
        el = list(orb6.elements(r))
        entry = {"name": name, "wds": wds, "disc": disc, "grade": int(r["grade"]), "ref": r["ref"],
                 "v1": float(r["v1"]), "v2": float(r["v2"]), "hip": int(r["hip"])}
        if name == "Sirius":
            el = list(np.median(sir_fit["hier"]["elements"], axis=0))
            plx, eplx, src = L.BOND17_PLX["adopted"]["value"] * 1e3, L.BOND17_PLX["adopted"]["err"] * 1e3, \
                "Bond et al. 2017 (Hipparcos + ground)"
            fB = L.BOND17_AA["value"] / el[3]
            msrc = "our orbit + Sirius A's measured orbit (Bond et al. 2017)"
        elif name == "Procyon":
            plx, eplx, src = L.BOND15_PROCYON["plx"][0] * 1e3, L.BOND15_PROCYON["plx"][1] * 1e3, \
                "Bond et al. 2015 (weighted mean)"
            fB = L.BOND15_PROCYON["aA"][0] / el[3]
            msrc = "Procyon A's measured orbit (Girard et al. 2000)"
        elif name == "α Centauri":
            plx, eplx, src = L.AKESON21_ACEN["plx"][0] * 1e3, L.AKESON21_ACEN["plx"][1] * 1e3, \
                "Akeson et al. 2021 (orbital parallax)"
            fB = L.AKESON21_ACEN["MB"][0] / (L.AKESON21_ACEN["MA"][0] + L.AKESON21_ACEN["MB"][0])
            msrc = "Akeson et al. 2021 masses"
        else:
            gg = pd.read_csv(C.GAIA_DIR / "gallery_gaia.csv").drop_duplicates("source_id")
            if name == "61 Cygni":
                c = gg[gg["source_id"].isin([1872046609345556480, 1872046574983497216])]
                w = 1 / c["parallax_error"] ** 2
                plx, eplx = float(np.sum(w * c["parallax"]) / w.sum()), float(1 / np.sqrt(w.sum()))
                src = "Gaia DR3 (A and B)"
            elif name == "ξ Ursae Majoris":
                c = gg[gg["dr2_parallax"].notna()].iloc[0]
                plx, eplx = float(c["dr2_parallax"]), float(c["dr2_parallax_error"])
                src = "Gaia DR2 (ξ UMa B; DR3 has no parallax for either star)"
            else:
                plx, eplx = float(hip2.loc[r["hip"], "plx"]), float(hip2.loc[r["hip"], "e_plx"])
                src = "Hipparcos (van Leeuwen 2007)"
            mu = 5 * np.log10(plx / 1e3) + 5
            m1, m2 = 10 ** like.f(best, r["v1"] + mu), 10 ** like.f(best, r["v2"] + mu)
            fB = float(m2 / (m1 + m2))
            msrc = "mass ratio from our mass-luminosity relation"
        Mtot = float(orbits.total_mass(el[3], plx / 1e3, el[0]))
        entry.update({"elements": [float(x) for x in el], "plx": plx, "e_plx": eplx, "plx_src": src,
                      "a_au": float(el[3] / (plx / 1e3)), "Mtot": Mtot,
                      "Mtot_err": float(Mtot * orbits.total_mass_frac_err(
                          el[3], r["a_err_arcsec"] if np.isfinite(r["a_err_arcsec"]) else 0.01 * el[3],
                          plx, eplx, el[0], r["P_err_yr"] if np.isfinite(r["P_err_yr"]) else 0.0)),
                      "fB": float(fB), "MA": Mtot * (1 - fB), "MB": Mtot * fB, "mass_src": msrc})
        out.append(entry)
    return out


def run_mlr(cache=True) -> dict:
    s = mlr.build_sample()
    sel = mlr.fit_selection(s)
    like, chain, best = mlr.fit_mlr(s[sel], cache=C.CACHE_DIR / "mlr_chain.npz" if cache else None)
    selg = mlr.fit_selection(s, parallax="gaia")
    MV, band = mlr.relation_band(like, chain)
    q = like.outlier_prob(best)
    s["in_fit"] = sel
    s.loc[sel, "p_outlier"] = q
    resid, var = like.residuals(best)
    s.loc[sel, "resid"] = resid
    s.loc[sel, "resid_err"] = np.sqrt(var)
    n = len(like.knots)
    # Gaia vs Hipparcos parallaxes for the same systems (where both exist)
    both = s[(s["plx_p"].notna()) & (s["plx_hip2"].notna()) & s["solved_p"].isin([31, 95])].copy()
    both["dplx_sig"] = (both["plx_p"] - both["plx_hip2"]) / np.hypot(both["e_plx_p"], both["e_plx_hip2"])
    mm = mlr.mamajek_table()
    masses = np.geomspace(0.18, 8.0, 200)
    return {"sample": s, "like": like, "chain": chain, "best": best, "MV": MV, "band": band,
            "n_fit": int(sel.sum()), "n_fit_gaia": int(selg.sum()), "sigma_int": float(np.exp(best[n])),
            "f_out": float(1 / (1 + np.exp(-best[n + 1]))), "n_outliers": int((q > 0.5).sum()),
            "gaia_vs_hip": both, "mamajek": mm, "eker_curve": {"M": masses, **mlr.eker_curve(masses)},
            "eker_stars": mlr.eker_stars()}


def run() -> dict:
    print("Sirius: orbit fits ...", flush=True)
    fit = sirius.run()
    print("Sirius: masses ...", flush=True)
    masses = sirius_masses(fit)
    print("Sirius B in Gaia DR3 ...", flush=True)
    gaia = sirius_gaia(fit)
    print("Sirius B as a white dwarf ...", flush=True)
    wd = sirius_b_white_dwarf(masses)
    print("ORB6: mass-luminosity relation ...", flush=True)
    ml = run_mlr()
    print("ORB6: gallery ...", flush=True)
    gal = gallery(fit, ml)
    return {"sirius": fit, "masses": masses, "gaia": gaia, "wd": wd, "mlr": ml, "gallery": gal}
