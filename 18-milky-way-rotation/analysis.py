"""Every measurement in the project, from cached data to numbers.

``run()`` returns a dict of results (also written to ``data/results.json``)
and caches the arrays the figures need in ``data/products/``. Nothing here
draws; ``figures.py`` does that from these products.
"""

from __future__ import annotations

import json
from itertools import product

import numpy as np
from scipy.ndimage import gaussian_filter1d

import config as cfg
import galaxy
import outer
import reference
import survey
import terminal

PRODUCTS = cfg.PRODUCTS_DIR


def _stats(d: np.ndarray) -> dict:
    d = np.asarray(d, dtype=float)
    d = d[np.isfinite(d)]
    med = float(np.median(d))
    return {"n": int(d.size), "median": med, "mean": float(d.mean()), "rms": float(d.std()),
            "mad_sigma": float(1.4826 * np.median(np.abs(d - med))),
            "max_abs": float(np.abs(d).max())}


# --------------------------------------------------------------------------
# Terminal velocities
# --------------------------------------------------------------------------

def terminal_velocities() -> dict:
    lv = {"mean": survey.lv_diagram(-cfg.TERMINAL_B_HALFWIDTH_DEG, cfg.TERMINAL_B_HALFWIDTH_DEG),
          "median": survey.lv_diagram(-cfg.TERMINAL_B_HALFWIDTH_DEG, cfg.TERMINAL_B_HALFWIDTH_DEG,
                                      stat="median")}
    out = {}
    for quad, rng in (("Q1", cfg.Q1_RANGE), ("Q4", cfg.Q4_RANGE)):
        base = terminal.measure(lv["mean"], rng)
        variants = {}
        for stat, thr, w_in in product(("mean", "median"), cfg.TERMINAL_THRESHOLDS_K,
                                       cfg.TERMINAL_W_IN_VARIANTS):
            res = terminal.measure(lv[stat], rng, threshold=thr, w_in=w_in)
            variants[f"{stat}_{thr:g}K_{w_in:g}"] = res["v_terminal"]
        stack = np.vstack(list(variants.values()))
        sign = terminal.quadrant_sign(base["l"])
        half = np.array([terminal.half_power_edge(lv["mean"]["v"], lv["mean"]["T"][i], int(s))
                         for i, s in zip(np.searchsorted(lv["mean"]["l"], base["l"]), sign)])
        spread = np.nanstd(stack, axis=0)
        out[quad] = {**base,
                     "v_method_spread": spread,
                     "reliable": base["ok"] & (spread <= cfg.TERMINAL_MAX_METHOD_SPREAD),
                     "variant_names": np.array(list(variants)),
                     "variants": stack,
                     "v_half_power": half}
    np.savez_compressed(PRODUCTS / "terminal.npz",
                        **{f"{q}_{k}": v for q, d in out.items() for k, v in d.items()})
    return out


def compare_terminal(tv: dict) -> dict:
    q1, q4 = tv["Q1"], tv["Q4"]
    res = {}
    # --- First quadrant: point by point against the VLA survey ---
    m = reference.mgd16_q1()
    ours = np.interp(m["l"], q1["l"], q1["v_terminal"])
    step = np.median(np.diff(m["l"]))
    m_smooth = gaussian_filter1d(m["v"], cfg.BEAM_FWHM_DEG / 2.3548 / step, mode="nearest")
    res["q1_vs_mgd16"] = _stats(ours - m["v"])
    res["q1_vs_mgd16_smoothed"] = _stats(ours - m_smooth)
    res["q1_vs_mgd16_l_range"] = [float(m["l"].min()), float(m["l"].max())]

    # --- Fourth quadrant: the 10 values printed in MGD07's Table 1 ---
    l07 = np.array([339.695, 339.628, 339.562, 339.495, 339.428,
                    339.362, 339.295, 339.228, 339.162, 339.095])
    v07 = np.array([-135.22, -135.14, -135.44, -136.17, -135.59,
                    -137.87, -136.53, -135.63, -135.47, -134.11])
    res["q4_vs_mgd07_table"] = _stats(np.interp(l07, q4["l"], q4["v_terminal"]) - v07)

    # --- Linear fits |V_t| = v_circ (1 - |sin l|) + v' over MGD's ranges ---
    def linfit(q, lo, hi, s_lo=0.0, s_hi=1.0):
        s = np.abs(np.sin(np.radians(q["l"])))
        sel = q["reliable"] & (q["l"] >= lo) & (q["l"] <= hi) & (s >= s_lo) & (s <= s_hi)
        A = np.vstack([1 - s[sel], np.ones(sel.sum())]).T
        coef, *_ = np.linalg.lstsq(A, np.abs(q["v_terminal"][sel]), rcond=None)
        resid = np.abs(q["v_terminal"][sel]) - A @ coef
        cov = np.linalg.inv(A.T @ A) * resid.var(ddof=2)
        return coef, np.sqrt(np.diag(cov)), int(sel.sum())

    (a1, b1), _, _ = linfit(q1, 18.4, 67.0)
    (a4, b4), _, _ = linfit(q4, 285.0, 339.5)
    res["linear_fit"] = {"Q1": {"v_circ": a1, "v_prime": b1, "published": [172.5, 20.7]},
                         "Q4": {"v_circ": a4, "v_prime": b4, "published": [177.7, 18.4]}}
    # Oort A x R0 from the slope near the solar circle (MGD07 section 4.2)
    ar0 = {}
    for name, q, sub in (("Q4", q4, (285.0, 345.0)), ("Q1", q1, (15.0, 75.0))):
        for s_lo, s_hi in ((0.80, 0.95), (0.35, 0.95), (0.85, 0.95), (0.80, 0.90)):
            (a, _), (ea, _), n = linfit(q, *sub, s_lo, s_hi)
            ar0[f"{name}_{s_lo:.2f}-{s_hi:.2f}"] = {"AR0": a / 2, "err": ea / 2, "n": n}
    res["AR0"] = ar0
    res["AR0_published_MGD07_Q4"] = {"0.80-0.95": [110, 2], "0.35-0.95": [89, 4]}

    # --- Clemens (1985) CO terminal velocities (R0 8.5, Theta0 220) ---
    lq = q1["l"]
    v_co = galaxy.theta_clemens85(cfg.R0_IAU_KPC * np.sin(np.radians(lq))) \
        - cfg.THETA0_IAU_KMS * np.sin(np.radians(lq))
    sel = lq > 22
    res["q1_clemens_minus_ours_l_gt_22"] = _stats(v_co[sel] - q1["v_terminal"][sel])

    # --- North-south asymmetry ---
    L = np.arange(18.5, 67.0 + 1e-9, 1 / survey.PIX_PER_DEG)
    v1 = np.interp(L, q1["l"], q1["v_terminal"])
    v4 = -np.interp(360 - L, q4["l"], q4["v_terminal"])
    diff = v1 - v4
    diff_s = gaussian_filter1d(diff, 1.0 * survey.PIX_PER_DEG / 2.3548)
    res["ns_asymmetry"] = {
        "at_25": float(np.interp(25, L, diff_s)), "at_40": float(np.interp(40, L, diff_s)),
        "at_65": float(np.interp(65, L, diff_s)),
        "mean_55_67": float(diff[(L >= 55) & (L <= 67)].mean()),
        "range_55_67": [float(diff[(L >= 55) & (L <= 67)].min()),
                        float(diff[(L >= 55) & (L <= 67)].max())],
        "published": {"at_25": -15, "at_40": 10, "at_65": -10, "mean_55_67": -5.2},
    }
    # can an outward motion of the LSR (Delta v = 2 Pi0 cos l) explain it?
    c = np.cos(np.radians(L))
    pi0 = float(np.sum(diff * 2 * c) / np.sum((2 * c) ** 2))
    res["ns_asymmetry"]["lsr_outflow_fit_Pi0"] = pi0
    res["ns_asymmetry"]["rms_before"] = float(diff.std())
    res["ns_asymmetry"]["rms_after_Pi0"] = float((diff - 2 * pi0 * c).std())
    np.savez_compressed(PRODUCTS / "ns_asymmetry.npz", l=L, diff=diff, diff_smooth=diff_s)

    # --- Method sensitivity: how each variant shifts the curve ---
    sens = {}
    for quad in ("Q1", "Q4"):
        q = tv[quad]
        base = q["v_terminal"]
        sens[quad] = {name: float(np.nanmedian(np.abs(v) - np.abs(base)))
                      for name, v in zip(q["variant_names"], q["variants"])}
        sens[quad]["median_spread"] = float(np.nanmedian(q["v_method_spread"]))
        sens[quad]["half_power_minus_fit"] = _stats(np.abs(q["v_half_power"]) - np.abs(base))
    res["method_sensitivity"] = sens
    res["reliability"] = {
        quad: {"n": int(tv[quad]["l"].size), "reliable": int(tv[quad]["reliable"].sum()),
               "flagged_l": [float(x) for x in tv[quad]["l"][~tv[quad]["reliable"]]]}
        for quad in ("Q1", "Q4")}
    return res


# --------------------------------------------------------------------------
# Rotation curve
# --------------------------------------------------------------------------

def _binned(R, theta, edges):
    idx = np.digitize(R, edges) - 1
    med = np.full(edges.size - 1, np.nan)
    sd = np.full(edges.size - 1, np.nan)
    n = np.zeros(edges.size - 1, int)
    for i in range(edges.size - 1):
        sel = (idx == i) & np.isfinite(theta)
        n[i] = sel.sum()
        if n[i] >= 3:
            med[i] = np.median(theta[sel])
            sd[i] = theta[sel].std()
    return med, sd, n


def rotation_curve(tv: dict, *, n_mc: int = 400, seed: int = 2) -> dict:
    obs = galaxy.Observer.reid19()
    obs_iau = galaxy.Observer.textbook(cfg.R0_IAU_KPC, cfg.THETA0_IAU_KMS)
    out = {"quadrants": {}}
    edges = np.arange(2.0, 8.26, 0.25)
    centres = 0.5 * (edges[1:] + edges[:-1])
    rng = np.random.default_rng(seed)
    for quad in ("Q1", "Q4"):
        q = tv[quad]
        s = np.abs(np.sin(np.radians(q["l"])))
        good = q["reliable"] & (s <= 0.95)
        R = galaxy.tangent_radius(q["l"], obs.r0)
        th = galaxy.theta_from_terminal(q["l"], q["v_terminal"], obs)
        th_iau = galaxy.theta_from_terminal(q["l"], q["v_terminal"], obs_iau)
        med, sd, n = _binned(R[good], th[good], edges)
        # Galactic-constant uncertainty: redo with R0, the Sun's angular
        # velocity and U_sun drawn from their published uncertainties
        mc = np.full((n_mc, centres.size), np.nan)
        for k in range(n_mc):
            o = galaxy.Observer.reid19(
                r0=rng.normal(cfg.R0_KPC, cfg.R0_ERR_KPC),
                omega_sun=rng.normal(cfg.OMEGA_SUN, cfg.OMEGA_SUN_ERR),
                u_sun=rng.normal(cfg.U_SUN_KMS, cfg.U_SUN_ERR_KMS))
            Rk = galaxy.tangent_radius(q["l"], o.r0)
            tk = galaxy.theta_from_terminal(q["l"], q["v_terminal"], o)
            mc[k], _, _ = _binned(Rk[good], tk[good], edges)
        out["quadrants"][quad] = {
            "l": q["l"], "R": R, "theta": th, "theta_iau": th_iau, "good": good,
            "theta_err_method": q["v_method_spread"],
            "bin_centres": centres, "bin_median": med, "bin_scatter": sd, "bin_n": n,
            "bin_sys_lo": np.nanpercentile(mc, 16, axis=0),
            "bin_sys_hi": np.nanpercentile(mc, 84, axis=0),
        }
    # linear fits
    def fit(Rs, ths, lo, hi):
        sel = (Rs >= lo) & (Rs <= hi) & np.isfinite(ths)
        A = np.vstack([np.ones(sel.sum()), Rs[sel] - cfg.R0_KPC]).T
        coef, *_ = np.linalg.lstsq(A, ths[sel], rcond=None)
        resid = ths[sel] - A @ coef
        return coef, float(resid.std())

    fits = {}
    for name, quads in (("Q1", ["Q1"]), ("Q4", ["Q4"]), ("joint", ["Q1", "Q4"])):
        Rs = np.concatenate([out["quadrants"][qq]["R"][out["quadrants"][qq]["good"]] for qq in quads])
        ts = np.concatenate([out["quadrants"][qq]["theta"][out["quadrants"][qq]["good"]] for qq in quads])
        (t0, slope), rms = fit(Rs, ts, 4.0, 7.8)
        fits[name] = {"theta_at_R0": float(t0), "slope": float(slope), "rms": rms}
    out["linear_fits_4_7.8kpc"] = fits

    # MGD-style fit in the IAU convention: Theta/Theta0 = k0 R/R0 + k1
    Rs, ts = [], []
    for quad in ("Q1", "Q4"):
        q = out["quadrants"][quad]
        Rs.append(galaxy.tangent_radius(q["l"][q["good"]], cfg.R0_IAU_KPC))
        ts.append(q["theta_iau"][q["good"]])
    Rs, ts = np.concatenate(Rs), np.concatenate(ts)
    sel = (Rs >= 3.0) & (Rs <= 8.0)
    A = np.vstack([Rs[sel] / cfg.R0_IAU_KPC, np.ones(sel.sum())]).T
    (k0, k1), *_ = np.linalg.lstsq(A, ts[sel] / cfg.THETA0_IAU_KMS, rcond=None)
    out["mgd_style_fit"] = {"k0": float(k0), "k1": float(k1), "published_2016": [0.171, 0.889],
                            "published_2007": [0.186, 0.887]}
    return out


def masers() -> dict:
    m = reference.reid19_masers()
    r = reference.maser_rotation(m)
    sel = ((r["frac_plx_err"] < 0.2) & (r["R"] > 4)
           & ~np.isin(m.name, list(reference.REID19_OUTLIERS)))
    good = sel & (r["vrot_err"] < 20)
    urc = galaxy.theta_reid19(r["R"])
    d = r["vrot"] - urc
    # Perseus arm: kinematic distance (circular rotation) versus parallax
    obs = galaxy.Observer.reid19()
    per = (m.arm == "Per") & (r["frac_plx_err"] < 0.2)
    _, d_kin = galaxy.outer_distance(m.l[per], m.b[per], m.vlsr[per], galaxy.theta_reid19, obs)
    d_plx = 1 / m.plx[per]
    ratio = d_kin / d_plx
    np.savez_compressed(PRODUCTS / "masers.npz", name=m.name, arm=m.arm, l=m.l, b=m.b,
                        vlsr=m.vlsr, plx=m.plx, e_plx=m.e_plx, R=r["R"], vrot=r["vrot"],
                        vrad=r["vrad"], R_err=r["R_err"], vrot_err=r["vrot_err"],
                        selected=sel, good=good)
    return {
        "n_total": int(m.name.size), "n_selected": int(sel.sum()), "n_good": int(good.sum()),
        "vrot_minus_reid_urc": _stats(d[good]),
        "vrad_median_outward": float(np.median(r["vrad"][good])),
        "published_streaming": {"Us_toward_centre": [6.0, 1.4], "Vs": [-4.3, 5.6]},
        "perseus_kinematic_over_parallax": {"n": int(np.isfinite(ratio).sum()),
                                            "median": float(np.nanmedian(ratio)),
                                            "range": [float(np.nanmin(ratio)),
                                                      float(np.nanmax(ratio))]},
    }


def milky_way_in_context(rc: dict) -> dict:
    sp = reference.sparc()
    vflat = np.array([g["v_flat"] for g in sp.values()])
    vflat = vflat[np.isfinite(vflat)]
    eil = reference.eilers19()
    fits = rc["linear_fits_4_7.8kpc"]["joint"]
    masses = {}
    for R in (4.0, 6.0, 7.5):
        th = fits["theta_at_R0"] + fits["slope"] * (R - cfg.R0_KPC)
        masses[f"{R:g}"] = float(galaxy.enclosed_mass(R, th))
    masses["R0_reid_theta0"] = float(galaxy.enclosed_mass(cfg.R0_KPC, cfg.THETA0_KMS))
    i = np.argmax(eil["R"])
    masses["eilers_outermost"] = {"R": float(eil["R"][i]), "v": float(eil["vc"][i]),
                                  "M": float(galaxy.enclosed_mass(eil["R"][i], eil["vc"][i]))}
    return {
        "sparc_n_with_vflat": int(vflat.size),
        "sparc_fraction_faster_than_mw": float(np.mean(vflat > cfg.THETA0_KMS)),
        "sparc_fraction_faster_than_eilers": float(np.mean(vflat > 229.0)),
        "sparc_vflat_percentiles": [float(x) for x in np.percentile(vflat, [10, 50, 90])],
        "enclosed_mass_Msun": masses,
    }


# --------------------------------------------------------------------------
# The outer disc
# --------------------------------------------------------------------------

LBH06_ELLIPTIC = (8.67, -1.08)  # LBH06 Appendix, fitted to the LAB survey


def _lbh06_setup(elliptic=None) -> outer.OuterSetup:
    """Levine, Blitz & Heiles' own assumptions: R0 = 8.5 kpc, flat 220 km/s."""
    return outer.OuterSetup(theta_fn=lambda r: galaxy.theta_flat(r, cfg.THETA0_IAU_KMS),
                            obs=galaxy.Observer.textbook(cfg.R0_IAU_KPC, cfg.THETA0_IAU_KMS),
                            r_min_kpc=cfg.R0_IAU_KPC + 0.5, elliptic=elliptic)


def elliptic_fits(*, refresh: bool = False) -> dict:
    """The elliptical-orbit correction, re-derived from HI4PI (cached)."""
    path = PRODUCTS / "elliptic_fit.json"
    if path.exists() and not refresh:
        return json.loads(path.read_text())
    ours = outer.fit_elliptic(outer.OuterSetup())
    theirs = outer.fit_elliptic(_lbh06_setup())
    out = {"ours": {k: v for k, v in ours.items() if k != "trace"},
           "lbh06_constants": {k: v for k, v in theirs.items() if k != "trace"},
           "lbh06_published": list(LBH06_ELLIPTIC)}
    path.write_text(json.dumps(out, indent=1))
    return out


def outer_disc(*, refresh: bool = False) -> dict:
    ell = elliptic_fits(refresh=refresh)
    configs = {
        # our constants, with the elliptical-orbit correction fitted to our data
        "default": outer.OuterSetup(elliptic=(ell["ours"]["alpha"], ell["ours"]["beta"])),
        # pure circular orbits, to show what the correction does
        "circular": outer.OuterSetup(),
        # Levine, Blitz & Heiles' own assumptions and published correction,
        # for a like-for-like check of their warp measurements
        "lbh06": _lbh06_setup(LBH06_ELLIPTIC),
    }
    out = {}
    for name, st in configs.items():
        path = PRODUCTS / f"outer_{name}.npz"
        if path.exists() and not refresh:
            with np.load(path) as z:
                prod = {k: z[k] for k in z.files}
        else:
            prod = outer.build(st)
            np.savez_compressed(path, **prod)
        vs = outer.vertical_structure(prod)
        rings = []
        for lo in np.arange(11.0, 30.0, 1.0):
            w = outer.fit_warp(vs, lo, lo + 1.0)
            if w["ok"]:
                rings.append({"R": lo + 0.5, **{k: float(w[k]) for k in
                                                ("W0", "W1", "phi1", "W2", "phi2", "rms")}})
        rc = vs["r_centres"]
        valid = np.isfinite(vs["zbar"])
        # total surface density (all heights), averaged over the azimuths mapped
        sig_R = np.array([np.nanmean(np.where(valid[i], vs["sigma_all"][i], np.nan))
                          for i in range(rc.size)])
        disp_R = np.nanmedian(vs["disp"], axis=1)
        out[name] = {"rings": rings, "sigma_R": sig_R, "disp_R": disp_R, "r_centres": rc,
                     "total_mass": float(prod["rpz_mass"].sum())}
        np.savez_compressed(PRODUCTS / f"outer_{name}_vertical.npz", **vs)
    # exponential scale length of the surface density, 12-25 kpc
    for d in out.values():
        sel = (d["r_centres"] >= 12) & (d["r_centres"] <= 25) & (d["sigma_R"] > 0)
        slope, icpt = np.polyfit(d["r_centres"][sel], np.log(d["sigma_R"][sel]), 1)
        d["sigma_scale_length_kpc"] = float(-1 / slope)
    return out


def run(*, refresh_outer: bool = False) -> dict:
    tv = terminal_velocities()
    results = {"terminal": compare_terminal(tv)}
    rc = rotation_curve(tv)
    np.savez_compressed(PRODUCTS / "rotation.npz",
                        **{f"{q}_{k}": v for q, d in rc["quadrants"].items() for k, v in d.items()})
    results["rotation"] = {k: v for k, v in rc.items() if k != "quadrants"}
    results["rotation"]["binned"] = {
        q: {"R": d["bin_centres"].tolist(), "theta": d["bin_median"].tolist(),
            "scatter": d["bin_scatter"].tolist(), "sys_lo": d["bin_sys_lo"].tolist(),
            "sys_hi": d["bin_sys_hi"].tolist(), "n": d["bin_n"].tolist()}
        for q, d in rc["quadrants"].items()}
    results["masers"] = masers()
    results["context"] = milky_way_in_context(rc)
    od = outer_disc(refresh=refresh_outer)
    results["outer_elliptic_fit"] = elliptic_fits()
    results["outer"] = {name: {"rings": d["rings"], "total_mass": d["total_mass"],
                               "sigma_scale_length_kpc": d.get("sigma_scale_length_kpc"),
                               "sigma_R": dict(zip(np.round(d["r_centres"], 2).tolist(),
                                                   np.round(d["sigma_R"], 4).tolist())),
                               "disp_R": dict(zip(np.round(d["r_centres"], 2).tolist(),
                                                  np.round(d["disp_R"], 4).tolist()))}
                        for name, d in od.items()}
    cfg.RESULTS_PATH.write_text(json.dumps(_jsonable(results), indent=1))
    return results


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist())
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else round(float(x), 6)
    if isinstance(x, (np.integer,)):
        return int(x)
    return x
