"""Every result of the project, computed from the transcribed tables and the downloads.

    python analysis.py        # run everything, print a summary, write data/results.json

Steps (each a function returning plain numbers and small tables):

* ``frame``: Table I's standard coordinates against Gaia DR3 carried back to
  1919, and the 1919 predicted displacements against ours (JPL Sun, exact GR).
* ``contacts``: the Moon and the Sun from both stations (JPL Horizons): contact
  times against the ones the expeditions assumed (Sobral) or computed (Principe).
* ``plate1``: the paper's photograph: its stars identified and its exposure
  dated from the Moon's position.
* ``sobral_4inch``, ``astrographic``, ``principe``: the 1920 reductions
  replicated from the printed measures, and refitted.
* ``wallal``: the 1922 Lick five-foot cameras' radial displacements refitted.
* ``history``: a century of gamma; ``shapiro`` and ``cassini``: the 2002
  conjunction (geometry, and the fit of the raw Doppler data if it has run).
"""

from __future__ import annotations

import json
import warnings

import numpy as np
import pandas as pd

import config as C

warnings.filterwarnings("ignore")

#: probable error -> one sigma
PE_TO_SIGMA = 1.0 / 0.6745


# ---------------------------------------------------------------------------
# The star field and the Sun
# ---------------------------------------------------------------------------

def step_frame() -> dict:
    import stars as S

    f = S.fit_table1_frame()
    t1 = S.table1().set_index("no")
    out = {"tangent_ra_deg": f["ra0_deg"], "tangent_dec_deg": f["dec0_deg"], "rotation_deg": f["rotation_deg"],
           "rotation_err_deg": f["rotation_err_deg"], "scale": f["scale"], "scale_err": f["scale_err"],
           "rms_arcsec": f["rms_arcsec"],
           "resid_arcsec": {int(n): [float(a), float(b)] for n, a, b in zip(f["no"], f["resid_x_arcsec"], f["resid_y_arcsec"])},
           "sites": {}}
    for site in ("sobral", "principe"):
        pat = S.predicted_pattern(f, site, limb_value=C.DEFLECTION_LIMB_1919)
        imp = S.sun_from_table1(site)
        sun = S.sun_at(site, S.MID_TOTALITY[site])
        d = np.hypot(pat["dx"] - t1[f"{site}_dx"], pat["dy"] - t1[f"{site}_dy"])
        out["sites"][site] = {
            "sun_xy_horizons": [float(pat["sun_x"].iloc[0]), float(pat["sun_y"].iloc[0])],
            "sun_xy_implied_by_table1": [float(imp["sun_x"]), float(imp["sun_y"])],
            "sun_offset_arcsec": float(np.hypot(pat["sun_x"].iloc[0] - imp["sun_x"], pat["sun_y"].iloc[0] - imp["sun_y"])
                                       * C.UNIT_TABLE1_ARCMIN * 60),
            "limb_implied_by_table1": float(imp["limb_implied"]), "table1_pattern_rms": float(imp["rms_arcsec"]),
            "pred_vs_table1_max_arcsec": float(d.max()), "pred_vs_table1_rms_arcsec": float(np.sqrt(np.mean(d ** 2))),
            "sun_radius_arcmin": float(pat["sun_radius_arcmin"].iloc[0]), "sun_elevation_deg": float(pat["sun_el_deg"].iloc[0]),
            # GR at the limb as seen from the Earth that day (exact formula, Earth-Sun distance from JPL)
            "gr_limb_arcsec": float(S.deflection_arcsec(np.deg2rad(sun["radius_arcmin"] / 60), sun["delta_au"])),
            "earth_sun_au": sun["delta_au"],
        }
    hy = S.membership()
    out["hyades_members"] = [int(n) for n in hy.index[hy["hyades"]]]
    return out


def step_contacts() -> dict:
    """Second and third contact at each station from JPL's Moon and Sun."""
    import plate1 as P1

    tr = P1.moon_sun_track()
    out = {}
    for site, g in tr.groupby("site"):
        t = (g["t"] - g["t"].iloc[0]).dt.total_seconds().to_numpy()
        dx = (g["moon_ra"] - g["sun_ra"]).to_numpy() * np.cos(np.deg2rad(g["sun_dec"].to_numpy())) * 3600
        dy = (g["moon_dec"] - g["sun_dec"]).to_numpy() * 3600
        margin = (g["moon_diam_arcsec"] - g["sun_diam_arcsec"]).to_numpy() / 2 - np.hypot(dx, dy)
        tt = np.linspace(t[0], t[-1], 200001)
        mm = np.interp(tt, t, margin)
        cross = np.flatnonzero(np.diff(np.sign(mm)))
        c2, c3 = (g["t"].iloc[0] + pd.Timedelta(float(tt[i]), unit="s") for i in cross[:2])
        out[site] = {"second_contact_utc": c2.strftime("%H:%M:%S.%f")[:10], "third_contact_utc": c3.strftime("%H:%M:%S.%f")[:10],
                     "duration_s": (c3 - c2).total_seconds()}
    # what the expeditions assumed (p.300) and computed (p.314)
    out["sobral"]["recorded"] = {"second_contact_utc": C.SOBRAL_SECOND_CONTACT_UT[11:],
                                 "note": "assumed, p.300 (contacts not observed; no duration given)"}
    out["principe"]["recorded"] = {"second_contact_utc": "14:13:05", "third_contact_utc": "14:18:07",
                                   "note": "computed, p.314 (2h 13m 5s - 2h 18m 7s astronomical G.M.T.)"}
    return out


def step_plate1() -> dict:
    import plate1 as P1

    tr = P1.moon_sun_track()
    sob = tr[tr["site"] == "sobral"].reset_index(drop=True)
    centre = (np.deg2rad(sob["sun_ra"].iloc[len(sob) // 2]), np.deg2rad(sob["sun_dec"].iloc[len(sob) // 2]))
    im = P1.image()
    marks = P1.find_marks(im)
    ident = P1.identify(marks, centre)
    moon = P1.fit_moon(im)
    when = P1.date_exposure(ident, moon, centre, sob)
    out = {"rms_px": ident["rms_px"], "next_best_rms_px": ident["next_rms_px"],
           "scale_px_per_arcmin": ident["scale_px_per_arcmin"], "mirrored": ident["mirrored"],
           "stars": {int(k): {"x_px": float(v["x_px"]), "y_px": float(v["y_px"])} for k, v in ident["stars"].iterrows()},
           "moon": {"cx": float(moon["cx"]), "cy": float(moon["cy"]), "r_px": float(moon["r_px"]),
                    "r_arcmin": float(moon["r_px"] / ident["scale_px_per_arcmin"]),
                    "r_arcmin_horizons": float(when["moon_radius_pred_px"] / ident["scale_px_per_arcmin"]),
                    "edge_rms_px": moon["rms_px"]},
           "dated_utc": when["utc"].strftime("%H:%M:%S.%f")[:10], "sigma_t_s": when["sigma_t_s"],
           "miss_px": when["miss_px"], "speed_px_per_s": when["speed_px_per_s"],
           "affine": ident["affine"].tolist(), "centre_rad": list(centre)}
    out["exposure_match"] = _match_plate1_exposure(when["utc"], when["sigma_t_s"])
    return out


def _match_plate1_exposure(when: pd.Timestamp, sigma_t: float) -> dict | None:
    """Which 4-inch exposure (p.300) has its middle nearest the dated instant: with the 1919
    clock as recorded, and with second contact moved to JPL's time (8 s earlier)."""
    import plate1 as P1

    if not (C.TRANSCRIBED_DIR / "dyson1920_sobral_exposures.csv").exists():
        return None
    m = P1.match_exposure(when, sigma_t)
    m["jpl_clock"] = P1.match_exposure(when, sigma_t, clock_offset_s=-8.0)
    return m


# ---------------------------------------------------------------------------
# The 1919 reductions
# ---------------------------------------------------------------------------

def _limb(alpha_at_50: float, sun_radius_arcmin: float = C.SUN_RADIUS_TRUE_ARCMIN) -> float:
    """Deflection at 50' from the Sun's centre -> at the limb (1/r law). Our fits use the
    true apparent radius (15'.771); replications of the 1920 arithmetic use the paper's 15'.8."""
    return alpha_at_50 * C.UNIT_TABLE1_ARCMIN / sun_radius_arcmin


def step_sobral_4inch() -> dict:
    import eclipse1919 as E

    d = E.load_4inch()                       # with the two Table II values the 1920 reduction used
    rep = E.replicate_4inch(d)
    cmp = E.compare_printed_4inch(rep)
    sol = cmp[cmp["quantity"] == "sol_alpha"]
    rep_printed = E.replicate_4inch(E.load_4inch(as_printed=True))
    out = {"replication": {
        "alpha_ra_rev": rep["alpha_ra_rev"], "alpha_dec_rev": rep["alpha_dec_rev"], "alpha_rev": rep["alpha_rev"],
        "limb_ra": rep["limb_ra_1920"], "limb_dec": rep["limb_dec_1920"], "limb": rep["limb_1920"],
        "pe_fraction": rep["pe_fraction_1920"], "pe_fraction_by_coord": rep["pe_fraction_by_coord"],
        "limb_table2_as_printed": rep_printed["limb_1920"],
        "printed_vs_ours_n": int(len(cmp)), "printed_vs_ours_n_over_0.003": int((cmp["diff"].abs() > 0.003).sum()),
        "plate_alpha_max_abs_diff_rev": float(sol["diff"].abs().max()), "n_plate_alphas": int(len(sol)),
        "per_plate": rep["solutions"][["coord", "plate", "kind", "alpha", "scale", "orient"]].to_dict("records")},
        "table2_corrections": [{"plate": k[0], "star": k[1], "coord": k[2][:2], "used_1920": v}
                               for k, v in E.TABLE2_AS_REDUCED_1920.items()]}
    fits = {}
    for key, coords in (("both", ("RA", "Dec")), ("RA", ("RA",)), ("Dec", ("Dec",))):
        f = E.joint_fit_4inch(d, coords=coords)
        fits[key] = {"alpha_50": f["alpha"], "alpha_50_err": f["alpha_err"], "limb": _limb(f["alpha"]),
                     "limb_err": _limb(f["alpha_err"]), "chi2": f["chi2"], "dof": f["dof"],
                     "sigma_by_group": f["sigma_by_group"]}
    f = E.joint_fit_4inch(E.load_4inch(as_printed=True))
    fits["both_table2_as_printed"] = {"limb": _limb(f["alpha"]), "limb_err": _limb(f["alpha_err"])}
    out["joint_fit"] = fits
    jk = E.jackknife_4inch(d)
    n = len(jk) - 1
    vals = jk["alpha"].iloc[1:].to_numpy()
    def jk_err(prefix: str) -> float:
        v = np.array([_limb(a) for lo, a in zip(jk["left_out"], jk["alpha"]) if lo.startswith(prefix)])
        return float(np.sqrt((len(v) - 1) / len(v) * np.sum((v - v.mean()) ** 2)))

    out["jackknife"] = {"rows": [{"left_out": r.left_out, "limb": _limb(r.alpha), "limb_err": _limb(r.alpha_err)}
                                 for r in jk.itertuples()],
                        "spread_limb": float(_limb(np.ptp(vals))),
                        "err_plates": jk_err("plate"), "err_stars": jk_err("star")}
    return out


def step_astrographic() -> dict:
    import eclipse1919 as E

    d = E.load_astro()
    rep = E.replicate_astro(d)
    to_limb = C.REV_ASTRO * C.UNIT_TABLE1_ARCMIN / C.SUN_RADIUS_1919_ARCMIN
    out = {"replication": {"mean_coef_e": rep["mean_coef_e"], "mean_value_rev": rep["mean_value_rev"],
                           "limb_plates_scale": rep["limb_plates_scale"], "limb_refraction_scale": rep["limb_refraction_scale"],
                           "per_plate_lsq_limb_mean": rep["per_plate_lsq_limb_mean"], "table9_limb_mean": rep["table9_limb_mean"],
                           "per_plate": rep["per_plate"].to_dict("records")},
           "fits": {}}
    for label, kw in (("scale_common_free", {"scale": "common"}), ("scale_per_plate", {"scale": "per_plate"}),
                      ("scale_fixed_refraction", {"scale": "fixed", "e_fixed": 0.020}),
                      ("scale_fixed_plates", {"scale": "fixed", "e_fixed": 0.082})):
        f = E.fit_astro(d, **kw)
        row = {"limb": _limb(f["alpha"]), "limb_err": _limb(f["alpha_err"]), "sigma_arcsec": f["sigma"], "dof": f["dof"]}
        if "e_rev" in f:
            row.update({"e_rev": f["e_rev"], "e_err_rev": f["e_err_rev"], "corr_alpha_e": f["corr_alpha_e"]})
        out["fits"][label] = row
    # jackknife over the 16 plates, as for the 4-inch (one over the five stars means nothing here:
    # without star 5, the only one south of the Sun, the deflection is all but degenerate with
    # the scale, -1.3 +- 4.4")
    jk = np.array([_limb(E.fit_astro({**d, "Dy": d["Dy"].drop(index=p)}, scale="common")["alpha"]) for p in d["Dy"].index])
    out["fits"]["scale_common_free"]["limb_err_jackknife"] = float(np.sqrt((len(jk) - 1) / len(jk) * np.sum((jk - jk.mean()) ** 2)))
    e_grid = np.linspace(-0.02, 0.18, 201)
    out["alpha_vs_e"] = {"e_rev": e_grid.tolist(), "limb": [float(v) for v in rep["alpha_vs_e"](e_grid)]}
    out["to_limb_per_rev"] = to_limb
    return out


def step_principe() -> dict:
    import principe as PR

    d = PR.load()
    rep = PR.replicate(d)
    out = {"replication": {"limb": rep["limb"], "mean_four": rep["mean"], "mean_four_theta_ls": rep["mean_theta_ls"],
                           "per_comparison": rep["per_comparison"].to_dict("records")}}
    for key, kw in (("fit", {}), ("fit_1920_pattern", {"pattern": "1920"}), ("fit_free_scale", {"free_scale": True})):
        f = PR.fit(d, **kw)
        out[key] = {"limb": f["limb"], "limb_err": f["limb_err"], "limb_err_measures": f["limb_err_measures"],
                    "limb_err_scale": f["limb_err_scale"], "sigma_x_arcsec": f["sigma_x_units"] * PR.UNIT,
                    "sigma_y_arcsec": f["sigma_y_units"] * PR.UNIT, "dof": f["dof"],
                    "per_comparison": {c: {"limb": v["limb"], "limb_err": v["limb_err"], "theta": v["theta"]}
                                       for c, v in f["per_comparison"].items()},
                    "plates": f["plates"]}
    return out


def step_wallal() -> dict:
    import wallal as W

    return W.analyse()


def step_combined_1919() -> dict:
    """The three 1919 instruments combined (our fits). Honest errors: for the 4-inch and the
    astrographic, the larger of the formal error and a jackknife over plates (and, for the
    4-inch, over stars). The chi-square says whether they agree at all; for two instruments,
    sigma_apart = sqrt(chi2)."""
    r = json.loads(C.RESULTS_PATH.read_text())
    s4 = r["sobral_4inch"]
    v4 = s4["joint_fit"]["both"]["limb"]
    e4 = max(s4["joint_fit"]["both"]["limb_err"], s4["jackknife"]["err_plates"], s4["jackknife"]["err_stars"])
    vp, ep = r["principe"]["fit"]["limb"], r["principe"]["fit"]["limb_err"]
    fa = r["astrographic"]["fits"]["scale_common_free"]
    va, ea = fa["limb"], max(fa["limb_err"], fa.get("limb_err_jackknife", 0.0))

    def comb(vals):
        v = np.array([x[0] for x in vals])
        w = 1 / np.array([x[1] for x in vals]) ** 2
        m = float(np.sum(w * v) / np.sum(w))
        chi2 = float(np.sum(w * (v - m) ** 2))
        dof = len(v) - 1
        err = float(1 / np.sqrt(np.sum(w)))
        out = {"limb": m, "limb_err": err, "chi2": chi2, "dof": dof,
               "limb_err_scaled": err * max(1.0, np.sqrt(chi2 / dof)) if dof else err}
        if dof == 1:
            out["sigma_apart"] = float(np.sqrt(chi2))
        return out

    return {"inputs": {"sobral_4inch": [v4, e4], "principe": [vp, ep], "astrographic": [va, ea]},
            "selection_1920": comb([(v4, e4), (vp, ep)]),
            "all": comb([(v4, e4), (vp, ep), (va, ea)]),
            "sobral_two": comb([(v4, e4), (va, ea)])}


# ---------------------------------------------------------------------------
# History, Shapiro delay, Cassini
# ---------------------------------------------------------------------------

def step_history() -> dict:
    import history as H

    t = H.timeline()
    fr = H.frontier(t)
    return {"timeline": t.to_dict("records"), "frontier": fr.to_dict("records")}


def step_shapiro() -> dict:
    import shapiro

    s = shapiro.summary()
    return {k: (v if not isinstance(v, np.generic) else float(v)) for k, v in s.items()}


def step_cassini() -> dict | None:
    """The fit of Cassini's raw Doppler data (cassini_doppler.py), if it has been run."""
    path = C.CACHE_DIR / "cassini_fit.json"
    return json.loads(path.read_text()) if path.exists() else None


# ---------------------------------------------------------------------------

STEPS = {"frame": step_frame, "contacts": step_contacts, "plate1": step_plate1,
         "sobral_4inch": step_sobral_4inch, "astrographic": step_astrographic, "principe": step_principe,
         "combined_1919": step_combined_1919, "wallal": step_wallal,
         "history": step_history, "shapiro": step_shapiro, "cassini": step_cassini}


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (pd.Timestamp,)):
        return str(x)
    if isinstance(x, float) and not np.isfinite(x):
        return None
    return x


def run(only: list[str] | None = None) -> dict:
    """Run the steps (all, or ``only`` these) and merge them into data/results.json."""
    res = json.loads(C.RESULTS_PATH.read_text()) if (only and C.RESULTS_PATH.exists()) else {}
    for name, fn in STEPS.items():
        if only and name not in only:
            continue
        print(f"  {name}", flush=True)
        res[name] = _jsonable(fn())
        C.RESULTS_PATH.write_text(json.dumps(res, indent=1))     # later steps read earlier ones
    return res


if __name__ == "__main__":
    import sys

    r = run(sys.argv[1:] or None)
    print(json.dumps({k: (v if k not in ("history",) else "...") for k, v in r.items()}, indent=1)[:6000])
