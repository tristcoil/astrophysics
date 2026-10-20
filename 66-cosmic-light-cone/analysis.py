"""Every number in the README: horizons, checks, real objects on the light cone,
superluminal recession measured in data, and the special-relativity test.

    python analysis.py      # writes data/results.json and data/derived/cone_objects.npz

Sections:
  1. our expansion history against astropy and CAMB
  2. the headline numbers for Planck 2018, against the spec and Planck itself
  3. Davis & Lineweaver (2004): their Figure 1 curves and text, recomputed
  4. how well the horizons are known (Planck 2018 and D3 uncertainties)
  5. real objects on our past light cone
  6. DESI: recession speed at emission, without knowing the size of the ruler
  7. Pantheon+: today's recession speeds of the supernovae
  8. is redshift a Doppler shift? D&L's supernova test, then and now
  9. what we can still reach
 10. how the answers depend on the cosmology
"""

from __future__ import annotations

import json
import sys
import time
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # works under python -I too

import numpy as np  # noqa: E402
from scipy.linalg import cho_factor, cho_solve  # noqa: E402

import catalogs  # noqa: E402
import config as C  # noqa: E402
import dl2004  # noqa: E402
from cosmology import Universe  # noqa: E402

OBJECTS_PATH = C.DERIVED_DIR / "cone_objects.npz"


def _f(x):
    """JSON-friendly floats (inf/nan -> strings)."""
    if isinstance(x, dict):
        return {k: _f(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_f(v) for v in x]
    if isinstance(x, (np.floating, float)):
        x = float(x)
        return x if np.isfinite(x) else ("inf" if x > 0 else "-inf" if x < 0 else "nan")
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.ndarray):
        return _f(x.tolist())
    return x


# ---------------------------------------------------------------------------
# 1. Our expansion history against astropy and CAMB
# ---------------------------------------------------------------------------


def planck18_camb_params():
    import camb

    h = C.PLANCK18.H0 / 100.0
    return camb.set_params(H0=C.PLANCK18.H0, ombh2=C.PLANCK18.Ob * h * h,
                           omch2=(C.PLANCK18.Om - C.PLANCK18.Ob) * h * h, mnu=0.06, nnu=3.046,
                           num_massive_neutrinos=1, TCMB=2.7255, tau=0.0561, As=2.105e-9,
                           ns=0.9665, omk=0)


def cosmology_checks(u: Universe) -> dict:
    from astropy.cosmology import Planck18 as AP
    from scipy.integrate import quad

    z = np.array([0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 5.0, 10.0, 14.44, 100.0, 1089.8])
    out = {"z": z}
    out["astropy_chi_rel"] = u.chi_of_z(z) / C.GLY_PER_MPC / AP.comoving_distance(z).value - 1
    out["astropy_t_rel"] = u.t_of_z(z) / AP.age(z).value - 1
    out["astropy_H_rel"] = u.H_of_z(z) / C.KMSMPC_TO_PER_GYR / AP.H(z).value - 1
    eta0_ap = quad(lambda a: C.C_KM_S / (AP.H0.value * a * a * AP.efunc(1 / a - 1)), 0, 1,
                   epsabs=0, epsrel=1e-12, limit=1000, points=[1e-6, 1e-4, 1e-2])[0]
    out["astropy_particle_horizon_gly"] = eta0_ap * C.GLY_PER_MPC
    out["astropy_age_gyr"] = AP.age(0).value
    try:
        import camb

        res = camb.get_background(planck18_camb_params())
        d = res.get_derived_params()
        out["camb_version"] = camb.__version__
        out["camb_derived"] = d
        out["camb_tau0_gly"] = res.conformal_time(0) * C.GLY_PER_MPC
        out["camb_chi_rel"] = u.chi_of_z(z) / C.GLY_PER_MPC / res.comoving_radial_distance(z) - 1
        out["camb_t_rel"] = u.t_of_z(z) / res.physical_time(z) - 1
        out["camb_da_star_rel"] = float(u.da_of_z(d["zstar"]) / C.GLY_PER_MPC
                                        / res.angular_diameter_distance(d["zstar"]) - 1)
        out["camb_omnuh2"] = res.Params.omnuh2
        out["our_omnuh2"] = u.Onu0 * (C.PLANCK18.H0 / 100) ** 2
    except ImportError:
        out["camb_derived"] = None
    out["max_abs_astropy_chi_rel"] = float(np.max(np.abs(out["astropy_chi_rel"])))
    out["max_abs_astropy_t_rel"] = float(np.max(np.abs(out["astropy_t_rel"])))
    if out.get("camb_derived"):
        out["max_abs_camb_chi_rel"] = float(np.max(np.abs(out["camb_chi_rel"])))
        out["max_abs_camb_t_rel"] = float(np.max(np.abs(out["camb_t_rel"])))
    return out


# ---------------------------------------------------------------------------
# 2. Headline numbers
# ---------------------------------------------------------------------------


def object_row(u: Universe, z: float, **extra) -> dict:
    """Where and when a source at redshift z sent the light we see today."""
    z = float(z)
    chi = float(u.chi_of_z(z))
    t_last = float(u.t_last_seen(np.array([chi]))[0])
    return {**extra, "z": z,
            "t_emit_gyr": float(u.t_of_z(z)), "lookback_gyr": u.age - float(u.t_of_z(z)),
            "chi_gly": chi, "d_emit_gly": float(u.emission_distance(z)),
            "dl_gly": float(u.dl_of_z(z)), "light_travel_gly": float(u.light_travel_distance(z)),
            "v_now_c": float(u.v_now(z)), "v_then_c": float(u.v_then(z)),
            "reachable": bool(chi < u.eta_inf - u.eta0),
            "t_last_seen_gyr": t_last}


def headline(u: Universe, z_star: float) -> dict:
    s = u.summary()
    z_turn = s["z_light_cone_meets_hubble_sphere"]
    t_turn = float(u.t_of_z(z_turn))
    out = {"summary": s, "spec_checks": {}, "planck_checks": {}}
    for k, v in C.SPEC_CHECKS.items():
        ours = s[k]
        digits = 2 if k == "reachable_fraction" else (1 if "gly" in k else 2)
        out["spec_checks"][k] = {"spec": v, "ours": ours, "ours_rounded": round(ours, digits),
                                 "agrees": abs(round(ours, digits) - v) < 0.5 * 10 ** (-digits) + 1e-9
                                 or (k == "reachable_fraction" and abs(ours - v) < 0.01)}
    out["planck_checks"]["age_gyr"] = {"planck": C.PLANCK18_PUBLISHED["age_gyr"], "ours": u.age}
    out["planck_checks"]["z_eq"] = {"planck": C.PLANCK18_PUBLISHED["z_eq"], "ours": u.z_equality()}
    out["z_star"] = z_star
    out["cmb"] = object_row(u, z_star, name="CMB last scattering")
    out["light_cone_turn"] = {"z": z_turn, "t_gyr": t_turn, "width_gly": float(u.da_of_z(z_turn)),
                              "note": "every photon now arriving was receding from us until this "
                                      "moment, when the Hubble sphere overtook it"}
    out["z_acceleration_start"] = u.z_acceleration_start()
    out["t_acceleration_start_gyr"] = float(u.t_of_z(out["z_acceleration_start"]))
    out["hubble_sphere_now"] = {"z": s["z_hubble_sphere"], "chi_gly": s["hubble_distance_gly"],
                                "t_emit_gyr": float(u.t_of_z(s["z_hubble_sphere"]))}
    out["event_horizon_now"] = {"z": s["z_event_horizon"], "chi_gly": s["event_horizon_gly"],
                                "t_emit_gyr": float(u.t_of_z(s["z_event_horizon"]))}
    # The comoving Hubble radius 1/(aH): it shrinks while the expansion accelerates.
    tt = np.geomspace(1e-3, 300.0, 4000)
    rh = 1.0 / (u.a_of_t(tt) * u.H_of_t(tt))
    i = int(np.argmax(rh))
    out["comoving_hubble_radius_max"] = {"gly": float(rh[i]), "t_gyr": float(tt[i]),
                                         "z": float(1 / u.a_of_t(tt[i]) - 1)}
    out["future_visibility_limit_gly"] = u.eta_inf
    out["observable_now_of_ever_observable"] = (u.eta0 / u.eta_inf) ** 3
    return out


# ---------------------------------------------------------------------------
# 3. Davis & Lineweaver (2004)
# ---------------------------------------------------------------------------


def _model_x(u: Universe, panel: str, curve: str, y: np.ndarray) -> np.ndarray:
    """Model distance (Gly) of a named curve at vertical coordinate y (t, or conformal time)."""
    if panel == "c":
        eta = np.clip(y, 1e-9, u.eta_inf - 1e-9)
        a = u.a_of_eta(eta)
    else:
        t = np.clip(y, 1e-9, None)
        a, eta = u.a_of_t(t), u.eta_of_t(t)
    chi = {"light_cone": lambda: u.eta0 - eta, "event_horizon": lambda: u.eta_inf - eta,
           "particle_horizon": lambda: eta,
           "hubble_sphere": lambda: 1.0 / (a * u.H_of_a(a))}[curve]()
    return chi * a if panel == "a" else chi


def dl2004_check(universes: dict) -> dict:
    src = C.LIT_DIR / "dl2004_src"
    out = {"panels": {}, "text": {}}
    for fname, panel, title in [("fig1a", "a", "proper distance vs time"),
                                ("fig1b", "b", "comoving distance vs time"),
                                ("fig1c", "c", "comoving distance vs conformal time")]:
        cur = dl2004.figure1_curves(src / f"{fname}.eps")
        frame = cur["frame"]
        rec = {"title": title, "calibration_rms_page": [cur["calibration"].rms_x, cur["calibration"].rms_y],
               "curves": {}}
        for curve in ("light_cone", "event_horizon", "particle_horizon", "hubble_sphere"):
            rec["curves"][curve] = {}
            for key, u in universes.items():
                res = []
                for w in cur[curve]:
                    x, y = np.abs(w[:, 0]), w[:, 1]
                    keep = (y > 0.2) & (x < frame[2] - 0.5) & (y < frame[3] - 0.05)
                    if curve == "light_cone" and panel != "c":
                        keep &= y < u.age - 0.05
                    if curve == "light_cone" and panel == "c":
                        keep &= y < u.eta0 - 0.05
                    res.append(x[keep] - _model_x(u, panel, curve, y[keep]))
                r = np.concatenate(res)
                rec["curves"][curve][key] = {"n": int(r.size), "median_gly": float(np.median(r)),
                                             "max_abs_gly": float(np.max(np.abs(r))),
                                             "rms_gly": float(np.sqrt(np.mean(r**2)))}
        # The scale-factor axis on the right: where their a(t) puts each labelled a.
        rec["scale_factor_axis"] = []
        for y, a in cur["right_axis_y_to_value"]:
            if not np.isfinite(a):
                rec["scale_factor_axis"].append({"a": "inf", "their_y": y,
                                                 "ours_no_rad": universes["dl2004"].eta_inf})
                continue
            row = {"a": a, "their_y": y}
            for key, u in universes.items():
                row[key] = float(u.eta_of_a(a)) if panel == "c" else float(u.t_of_a(a))
            rec["scale_factor_axis"].append(row)
        out["panels"][fname] = rec
        if panel == "a":
            out["fig1a_worldlines"] = _worldline_fits(cur, universes["dl2004"])
        if panel == "c":
            u = universes["dl2004"]
            lc = np.concatenate(cur["light_cone"])
            eh = np.concatenate(cur["event_horizon"])
            # In conformal coordinates the light cone and event horizon are 45-degree
            # lines, x + y = eta0 and x + y = eta_inf: read their eta0, eta_inf directly.
            out["their_eta0_gly"] = float(np.median(np.abs(lc[:, 0]) + lc[:, 1]))
            out["their_eta_inf_gly"] = float(np.median(np.abs(eh[:, 0]) + eh[:, 1]))
            out["our_eta0_gly"] = u.eta0
            out["our_eta_inf_gly"] = u.eta_inf
    out["text"] = dl2004_text_check(universes["dl2004"])
    out["fig5"] = None  # filled by sr_test
    return out


def _worldline_fits(cur: dict, u: Universe) -> list:
    """Fit each dotted worldline of Fig. 1a as D(t) = a(t) chi; compare with chi(z) for its label."""
    labels = []
    for t in dl2004.parse_grace_eps(C.LIT_DIR / "dl2004_src" / "fig1a.eps").texts:
        try:
            v = float(t.s)
        except ValueError:
            continue
        if t.size < 0.021 and v in (1.0, 3.0, 10.0, 1000.0):
            labels.append(cur["calibration"].to_world(np.array([[t.x, t.y]]))[0].tolist() + [v])
    out = []
    for w in cur["worldlines"]:
        x, y = np.abs(w[:, 0]), w[:, 1]
        keep = (y > 1.0) & (x < cur["frame"][2] - 0.5)
        a = u.a_of_t(y[keep])
        chi = float(np.sum(x[keep] * a) / np.sum(a * a))
        resid = x[keep] - a * chi
        top_x = float(w[np.argmax(w[:, 1]), 0])
        # Each label sits beside the end of its line (the top, or the side where the
        # line leaves the frame): take the closest label within 4 Gly / 1.5 Gyr.
        lab, best = None, np.inf
        for L in labels:
            d = np.min(np.hypot((w[:, 0] - L[0]) / 4.0, (w[:, 1] - L[1]) / 1.5))
            if d < 1.0 and d < best:
                lab, best = L[2], d
        row = {"side": "right" if top_x > 0 else "left", "chi_fit_gly": chi,
               "rms_resid_gly": float(np.sqrt(np.mean(resid**2))), "label_z": lab}
        if lab is not None:
            row["chi_ours_gly"] = float(u.chi_of_z(lab))
        out.append(row)
    return sorted(out, key=lambda r: (r["chi_fit_gly"], r["side"]))


def dl2004_text_check(u: Universe) -> dict:
    """The numbers D&L state in their text, recomputed for their model."""
    z_turn = u.z_light_cone_meets_hubble_sphere()
    eta_half = 0.5 * u.eta_inf
    out = {
        "z_hubble_sphere": u.z_hubble_sphere_now(),
        "z_event_horizon": u.z_event_horizon_now(),
        "particle_horizon_gly": u.eta0,
        "age_gyr": u.age,
        "cmb_v_now_c": float(u.v_now(1100.0)),
        "cmb_v_then_c": float(u.v_then(1100.0)),
        "superluminal_emission_before_gyr": float(u.t_of_z(z_turn)),
        "symmetry_worldline_z": float(u.z_of_chi(eta_half)),
    }
    return {k: {"theirs": C.DL2004_TEXT[k][0], "where": C.DL2004_TEXT[k][1], "ours": v}
            for k, v in out.items()}


# ---------------------------------------------------------------------------
# 4. Uncertainties
# ---------------------------------------------------------------------------


def _ln_cov_from_published(s_lnA, s_lnB, s_target, combo):
    """Covariance of two log-parameters from their errors plus the error of one
    linear combination combo = (ca, cb) of them (e.g. ln Om = 3 ln(Om h^2) - 2 ln(Om h^3))."""
    ca, cb = combo
    cov_ab = (s_target**2 - ca**2 * s_lnA**2 - cb**2 * s_lnB**2) / (2 * ca * cb)
    return np.array([[s_lnA**2, cov_ab], [cov_ab, s_lnB**2]])


QUANTITIES = ("age_gyr", "particle_horizon_gly", "event_horizon_gly", "hubble_distance_gly",
              "z_hubble_sphere", "z_event_horizon", "z_min_angular_size", "reachable_fraction")


def propagate(base, Om_incl, h, cov_lnOm_lnh, label: str) -> dict:
    """Linear error propagation of the headline numbers in (ln Omega_m, ln h)."""
    omega_nu = Universe(base, n_cells=1500).Onu0 * (base.H0 / 100) ** 2

    def summary(dlnOm, dlnh):
        hh = h * np.exp(dlnh)
        p = base.but(H0=100 * hh, Om=Om_incl * np.exp(dlnOm) - omega_nu / hh**2)
        return Universe(p, n_cells=3000).summary()

    step = 2e-3
    s0 = summary(0, 0)
    sp = [summary(step, 0), summary(0, step)]
    sm = [summary(-step, 0), summary(0, -step)]
    out = {"label": label, "Om_incl_nu": Om_incl, "h": h, "cov_ln": cov_lnOm_lnh.tolist(),
           "corr": float(cov_lnOm_lnh[0, 1] / np.sqrt(cov_lnOm_lnh[0, 0] * cov_lnOm_lnh[1, 1]))}
    for q in QUANTITIES:
        J = np.array([(sp[i][q] - sm[i][q]) / (2 * step) for i in range(2)])
        out[q] = {"value": s0[q], "sigma": float(np.sqrt(J @ cov_lnOm_lnh @ J)),
                  "dlnOm": float(J[0]), "dlnh": float(J[1])}
    return out


def uncertainties() -> dict:
    P = C.PLANCK18_PUBLISHED
    # Planck 2018 (+BAO): errors of Om h^2 (0.14240 +- 0.00087), Om h^3 and Om fix the
    # 2x2 covariance; it reproduces the published sigma(H0) = 0.42 as a check.
    sa, sb = 0.00087 / 0.14240, P["Om_h3"][1] / P["Om_h3"][0]
    s_lnOm = P["Om_incl_nu"][1] / P["Om_incl_nu"][0]
    cov_ab = _ln_cov_from_published(sa, sb, s_lnOm, (3.0, -2.0))
    # (ln Om, ln h) = M (ln Om h^2, ln Om h^3) with M = [[3, -2], [-1, 1]]
    M = np.array([[3.0, -2.0], [-1.0, 1.0]])
    cov_planck = M @ cov_ab @ M.T
    out = {"planck18": propagate(C.PLANCK18, P["Om_incl_nu"][0], P["H0"][0] / 100, cov_planck,
                                 "Planck 2018 TT,TE,EE+lowE+lensing+BAO"),
           "planck_sigma_H0_check": {"published": P["H0"][1],
                                     "from_cov": float(np.sqrt(cov_planck[1, 1]) * P["H0"][0])}}
    # D3's own TT fit (09-planck-cmb): errors of H0, Om and Om h^2 fix its covariance.
    d3 = json.loads((C.ROOT / "09-planck-cmb/data/results.json").read_text())["fit"]["derived"] \
        if (C.ROOT / "09-planck-cmb/data/results.json").exists() else None
    if d3:
        H0, sH0 = d3["H0"]["value"], d3["H0"]["error"]
        Om, sOm = d3["omegam"]["value"], d3["omegam"]["error"]
        w, sw = d3["omegamh2"]["value"], d3["omegamh2"]["error"]
        s_lnh, s_lnw, s_lnOm = sH0 / H0, sw / w, sOm / Om
        # ln Om = ln w - 2 ln h  ->  Var(ln Om) = Var(ln w) + 4 Var(ln h) - 4 Cov(ln w, ln h)
        cov_wh = (s_lnw**2 + 4 * s_lnh**2 - s_lnOm**2) / 4
        cov_lnOm_lnh = cov_wh - 2 * s_lnh**2
        cov = np.array([[s_lnOm**2, cov_lnOm_lnh], [cov_lnOm_lnh, s_lnh**2]])
        base = C.OTHER_UNIVERSES[0]
        out["d3_tt"] = propagate(base, Om, H0 / 100, cov, "D3: our CAMB fit to Planck TT")
    return out


# ---------------------------------------------------------------------------
# 5. Real objects on the light cone
# ---------------------------------------------------------------------------


def objects_on_cone(u: Universe, pops: dict, marks: list, z_star: float, *,
                    n_keep: int = 6000, seed: int = 20261019) -> tuple[dict, dict]:
    """Per-population summaries, plus arrays (subsampled for the big catalogues) for figures."""
    rng = np.random.default_rng(seed)
    z_hs = u.z_hubble_sphere_now()
    z_turn = u.z_light_cone_meets_hubble_sphere()
    z_eh = u.z_event_horizon_now()
    t_edges = np.linspace(0.0, u.age, 139)  # 0.1-Gyr bins
    summ, arrays = {}, {"t_edges": t_edges}
    for key, p in pops.items():
        z = p.z
        t = u.t_of_z(z)
        chi = u.chi_of_z(z)
        rec = {"label": p.label, "source": p.source, "note": p.note, "n": p.n,
               "z_min": float(z.min()), "z_median": float(np.median(z)), "z_max": float(z.max()),
               "t_emit_min_gyr": float(t.min()), "t_emit_median_gyr": float(np.median(t)),
               "lookback_max_gyr": float(u.age - t.min()),
               "chi_max_gly": float(chi.max()), "d_emit_max_gly": float(np.max(chi / (1 + z))),
               "v_now_max_c": float(u.H0 * chi.max()),
               "v_then_max_c": float(np.max(u.v_then(z))),
               "n_superluminal_now": int(np.sum(z > z_hs)),
               "n_superluminal_at_emission": int(np.sum(z > z_turn)),
               "n_beyond_event_horizon": int(np.sum(z > z_eh)),
               "frac_reachable": float(np.mean(z < z_eh))}
        tl = u.t_last_seen(chi[z > z_eh]) if np.any(z > z_eh) else np.array([])
        if tl.size:
            rec["t_last_seen_median_gyr"] = float(np.median(tl))
            rec["t_last_seen_min_gyr"] = float(np.min(tl))
        summ[key] = rec
        arrays[f"{key}_t_hist"] = np.histogram(t, bins=t_edges)[0]
        idx = np.arange(p.n) if p.n <= n_keep else np.sort(rng.choice(p.n, n_keep, replace=False))
        arrays[f"{key}_z"] = z[idx]
        if p.z_lo is not None:
            arrays[f"{key}_z_lo"] = p.z_lo[idx]
            arrays[f"{key}_z_hi"] = p.z_hi[idx]
    rows = [object_row(u, z_star, name="CMB (last scattering)", kind="cmb",
                       ref="CAMB with Planck 2018 parameters; Planck: 1089.80 +- 0.21")]
    for L in marks:
        rows.append(object_row(u, L["z"], name=L["name"], kind=L["kind"], ref=L["ref"],
                               note=L.get("note", "")))
    # The most distant gravitational-wave merger and supernova of the catalogues, by name.
    gw = pops["gw_events"]
    k = int(np.argmax(gw.z))
    rows.append(object_row(u, gw.z[k], name=gw.names[k], kind="gw", ref=gw.source,
                           note=f"90% range z = {gw.z_lo[k]:.2f}-{gw.z_hi[k]:.2f}"))
    for name in ("GW150914", "GW170817"):
        if name in gw.names:
            j = gw.names.index(name)
            rows.append(object_row(u, gw.z[j], name=name, kind="gw", ref=gw.source,
                                   note="the first merger seen" if name == "GW150914"
                                   else "the neutron-star merger (C4)"))
    sn = pops["supernovae"]
    k = int(np.argmax(sn.z))
    rows.append(object_row(u, sn.z[k], name=f"SN {sn.names[k]}", kind="supernova", ref=sn.source,
                           note="most distant supernova in Pantheon+"))
    q = pops["highz_quasars"]
    k = int(np.argmax(q.z))
    rows.append(object_row(u, q.z[k], name=q.names[k], kind="quasar", ref=q.source))
    rows.sort(key=lambda r: -r["z"])
    return {"populations": summ, "landmarks": rows,
            "z_hubble_sphere": z_hs, "z_turn": z_turn, "z_event_horizon": z_eh}, arrays


# ---------------------------------------------------------------------------
# 6. DESI: superluminal at emission, measured without the ruler's size
# ---------------------------------------------------------------------------


def desi_measurements(u: Universe, r_drag_mpc: float, *, n_mc: int = 200_000,
                      seed: int = 7) -> dict:
    """v_then/c = H(z) D_A / c = D_M / ((1+z) D_H): DESI's Alcock-Paczynski ratio
    F_AP = D_M/D_H divided by (1+z). The sound horizon r_d cancels, and so does H0.
    It exceeds 1 exactly where D_A decreases with z (flat FRW): the same redshifts
    at which objects of fixed size look bigger the farther they are."""
    d = catalogs.desi_bao()
    rng = np.random.default_rng(seed)
    rows = []
    zs = sorted({z for z, k in zip(d["z"], d["kind"]) if k == "DM_over_rs"})
    for z in zs:
        i = [j for j in range(d["z"].size) if d["z"][j] == z and d["kind"][j] == "DM_over_rs"][0]
        k = [j for j in range(d["z"].size) if d["z"][j] == z and d["kind"][j] == "DH_over_rs"][0]
        dm, dh = d["value"][i], d["value"][k]
        cov = d["cov"][np.ix_([i, k], [i, k])]
        r = dm / ((1 + z) * dh)
        g = np.array([1 / ((1 + z) * dh), -dm / ((1 + z) * dh**2)])
        sig = float(np.sqrt(g @ cov @ g))
        draws = rng.multivariate_normal([dm, dh], cov, n_mc)
        rmc = draws[:, 0] / ((1 + z) * draws[:, 1])
        da_rd = dm / (1 + z)
        rows.append({
            "z": float(z), "tracer": d["tracer"][i], "DM_over_rd": dm, "DH_over_rd": dh,
            "corr_DM_DH": float(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])),
            "v_then_c": float(r), "v_then_sigma": sig, "v_then_mc_mean": float(rmc.mean()),
            "v_then_mc_std": float(rmc.std()), "p_superluminal_mc": float(np.mean(rmc > 1)),
            "sigma_from_c": float((r - 1) / sig),
            "lcdm_planck18": float(u.v_then(z)),
            "DA_over_rd": float(da_rd), "DA_over_rd_sigma": float(np.sqrt(cov[0, 0]) / (1 + z)),
            "DA_mpc": float(da_rd * r_drag_mpc),
            "theta_1mpc_arcmin": float(1.0 / (da_rd * r_drag_mpc) * 180 / np.pi * 60),
            "theta_1mpc_arcmin_sigma": float(np.sqrt(cov[0, 0]) / dm * 1.0 / (da_rd * r_drag_mpc)
                                             * 180 / np.pi * 60),
        })
    # Model-independent crossing: interpolate v_then linearly between the two DESI
    # redshifts that bracket 1 (Monte Carlo over both), and the bias of doing so,
    # measured by applying the same interpolation to the Planck curve.
    brack = [(a, b) for a, b in zip(rows[:-1], rows[1:]) if (a["v_then_c"] - 1) * (b["v_then_c"] - 1) <= 0]
    cross = None
    if brack:
        a, b = brack[0]
        ra = rng.normal(a["v_then_c"], a["v_then_sigma"], n_mc)
        rb = rng.normal(b["v_then_c"], b["v_then_sigma"], n_mc)
        good = (ra < 1) & (rb > 1)
        zc = a["z"] + (1 - ra[good]) * (b["z"] - a["z"]) / (rb[good] - ra[good])
        za_true, zb_true = u.v_then(a["z"]), u.v_then(b["z"])
        z_interp_planck = a["z"] + (1 - za_true) * (b["z"] - a["z"]) / (zb_true - za_true)
        cross = {"between": [a["z"], b["z"]], "z_cross_median": float(np.median(zc)),
                 "z_cross_lo68": float(np.percentile(zc, 16)), "z_cross_hi68": float(np.percentile(zc, 84)),
                 "frac_draws_bracketing": float(good.mean()),
                 "interpolation_bias_on_planck": float(z_interp_planck - u.z_light_cone_meets_hubble_sphere())}
    # Flat LCDM through the ratios alone: they depend only on Omega_m.
    zz = np.array([r["z"] for r in rows])
    rv = np.array([r["v_then_c"] for r in rows])
    rs = np.array([r["v_then_sigma"] for r in rows])
    grid = np.linspace(0.20, 0.42, 441)
    chi2 = []
    for om in grid:
        uu = Universe(C.PLANCK18.but(Om=om - 0.00141), n_cells=1200)
        chi2.append(np.sum(((rv - uu.v_then(zz)) / rs) ** 2))
    chi2 = np.array(chi2)
    j = int(np.argmin(chi2))
    inside = grid[chi2 < chi2[j] + 1]
    om_best = float(grid[j])
    u_best = Universe(C.PLANCK18.but(Om=om_best - 0.00141), n_cells=3000)
    u_lo = Universe(C.PLANCK18.but(Om=inside.min() - 0.00141), n_cells=3000)
    u_hi = Universe(C.PLANCK18.but(Om=inside.max() - 0.00141), n_cells=3000)
    fit = {"Om_best": om_best, "Om_lo68": float(inside.min()), "Om_hi68": float(inside.max()),
           "chi2_min": float(chi2[j]), "n_points": int(zz.size),
           "z_cross": u_best.z_light_cone_meets_hubble_sphere(),
           "z_cross_range68": sorted([u_lo.z_light_cone_meets_hubble_sphere(),
                                      u_hi.z_light_cone_meets_hubble_sphere()]),
           "desi_dr2_bao_only_Om": [0.2975, 0.0086]}
    return {"rows": rows, "crossing": cross, "lcdm_ratio_fit": fit, "r_drag_mpc": r_drag_mpc,
            "dv_point": {"z": float(d["z"][0]), "DV_over_rd": float(d["value"][0])}}


# ---------------------------------------------------------------------------
# 7. Pantheon+: today's recession speeds of the supernovae
# ---------------------------------------------------------------------------


def _cosmography_dl(z, q0=-0.55, j0=1.0):
    """H0 D_L / c to third order (the SH0ES Hubble-flow convention, Riess et al. 2016)."""
    return z * (1 + 0.5 * (1 - q0) * z - (1 - q0 - 3 * q0**2 + j0) * z**2 / 6.0)


def sn_sample():
    t = catalogs.pantheon_table()
    cov = catalogs.pantheon_cov()
    rows = np.flatnonzero(t["zHD"] > 0.01)  # the supernova-only cosmology sample (D1)
    return t, cov, rows


def sn_measurements(u: Universe, bins=None) -> dict:
    """v_now/c = H0 D_M / c = H0 D_L / (c (1+z_hel)) for the supernovae (flat FRW).

    Supernovae measure H0 D_L up to one offset A = M_B + 5 log10(c/H0) + 25. A is
    set by the Hubble flow (0.0233 < z < 0.15, SH0ES's 277 supernovae) with the
    standard third-order expansion, so no cosmological model enters beyond that.
    Bins are generalised-least-squares means of the magnitude residuals about a
    smooth reference curve, with the full covariance. The last bin starts where
    Planck's v_now passes c (z = 1.476), so it holds exactly the supernovae that
    should be receding faster than light today.
    """
    if bins is None:
        bins = (0.01, 0.1, 0.3, 0.6, 1.0, round(u.z_hubble_sphere_now(), 4), 2.3)
    t, cov, rows = sn_sample()
    z, zh, m = t["zHD"], t["zHEL"], t["m_b_corr"]
    hf = np.flatnonzero(t["USED_IN_SH0ES_HF"] == 1)
    Chf = cho_factor(cov[np.ix_(hf, hf)])
    y = m[hf] - 5 * np.log10((1 + zh[hf]) / (1 + z[hf]) * _cosmography_dl(z[hf]))
    one = np.ones(hf.size)
    A = float(one @ cho_solve(Chf, y) / (one @ cho_solve(Chf, one)))
    sA = float(1 / np.sqrt(one @ cho_solve(Chf, one)))
    # residuals about Planck's shape (any smooth reference works: it only interpolates in bins)
    ref = 5 * np.log10(u.H0 * u.dl_of_z(z[rows]) * (1 + zh[rows]) / (1 + z[rows]))
    resid = m[rows] - A - ref
    Cs = cov[np.ix_(rows, rows)]
    cf = cho_factor(Cs)
    zb = np.digitize(z[rows], bins) - 1
    nb = len(bins) - 1
    W = np.zeros((rows.size, nb))
    W[np.arange(rows.size), np.clip(zb, 0, nb - 1)] = 1.0
    F = W.T @ cho_solve(cf, W)
    Fi = np.linalg.inv(F)
    dm = Fi @ (W.T @ cho_solve(cf, resid))
    dm_cov = Fi + sA**2  # the offset's error moves all bins together
    out_bins = []
    var = np.diag(Cs)
    for b in range(nb):
        sel = zb == b
        # Where to quote the bin: its inverse-variance mean redshift. (The GLS weights
        # themselves reach into other bins through correlated systematics, so they do
        # not define a redshift.) The residual about the smooth reference is nearly
        # flat inside a bin, so the choice barely matters.
        zeff = float(np.sum(z[rows][sel] / var[sel]) / np.sum(1 / var[sel]))
        v_ref = float(u.v_now(zeff))
        v = v_ref * 10 ** (dm[b] / 5)
        sv = v * np.log(10) / 5 * np.sqrt(dm_cov[b, b])
        out_bins.append({"z_lo": bins[b], "z_hi": bins[b + 1], "n": int(sel.sum()), "z_eff": zeff,
                         "v_now_c": float(v), "v_now_sigma": float(sv), "planck18": v_ref,
                         "sigma_from_c": float((v - 1) / sv)})
    # The individual supernovae beyond the Hubble sphere (diagonal errors + offset).
    far = np.flatnonzero(z[rows] > bins[-2])
    singles = []
    for i in far:
        zi = z[rows][i]
        v = float(u.v_now(zi)) * 10 ** (resid[i] / 5)
        sv = v * np.log(10) / 5 * np.sqrt(var[i] + sA**2)
        singles.append({"cid": str(t["CID"][rows][i]), "z": float(zi), "v_now_c": v,
                        "v_now_sigma": float(sv), "planck18": float(u.v_now(zi))})
    return {"offset_A": A, "offset_sigma": sA, "n_hubble_flow": int(hf.size),
            "n_sample": int(rows.size), "bins": out_bins, "beyond_hubble_sphere": singles,
            "bin_corr": (dm_cov / np.sqrt(np.outer(np.diag(dm_cov), np.diag(dm_cov)))).tolist()}


# ---------------------------------------------------------------------------
# 8. Is redshift a Doppler shift? D&L's supernova test, then and now
# ---------------------------------------------------------------------------


def beta_sr(z):
    """Special-relativistic Doppler velocity for redshift z, in units of c."""
    s = (1 + np.asarray(z, float)) ** 2
    return (s - 1) / (s + 1)


def hubble_free_dl(model: str, z, u: Universe | None = None):
    """H0 D_L / c for the readings of redshift that D&L compare.

    "SR": D&L's special-relativistic reading: Hubble's law v = H D with the Doppler
          velocity, H the same at emission as now, and D_L = D (1+z).
    "linear": v = cz in Hubble's law, D_L = cz (1+z) / H.
    "milne": the empty universe, the self-consistent special-relativistic picture
          (galaxies coasting apart through flat spacetime): D_L = (c/H0) z (1 + z/2).
    otherwise: a Friedmann universe u.
    """
    z = np.asarray(z, float)
    if model == "SR":
        return beta_sr(z) * (1 + z)
    if model == "SR_as_plotted":  # what D&L's Fig. 5 actually draws: D_L = D
        return beta_sr(z)
    if model == "linear":
        return z * (1 + z)
    if model == "linear_as_plotted":  # D&L's Fig. 5 v = cz curve: D_L = cz / H
        return z
    if model == "milne":
        return z * (1 + z / 2)
    return u.H0 * u.dl_of_z(z)


def _chi2_profile(m, C_fac, shape):
    """chi^2 with the additive magnitude offset minimised analytically."""
    one = np.ones(m.size)
    Ci1 = cho_solve(C_fac, one)
    d = m - shape
    Cid = cho_solve(C_fac, d)
    off = (one @ Cid) / (one @ Ci1)
    return float(d @ Cid - (one @ Cid) ** 2 / (one @ Ci1)), float(off)


def sr_test(universes: dict) -> dict:
    out = {}
    # (a) D&L's own data: the 42 Perlmutter et al. (1999) supernovae of their Fig. 5.
    f5 = dl2004.figure5_points(C.LIT_DIR / "dl2004_src" / "fig5.eps")
    z, mobs = f5["z"], f5["m"]
    err = 0.5 * (f5["err_up"] + f5["err_dn"])
    curly_m = C.DL2004_TEXT["script_M"][0]
    u_dl = universes["dl2004"]
    shapes = {name: 5 * np.log10(C.C_KM_S * hubble_free_dl(model, z, uu)) for name, model, uu in
              [("lcdm_0.3_0.7", "frw", u_dl), ("SR", "SR", None), ("SR_as_plotted", "SR_as_plotted", None),
               ("linear", "linear", None), ("linear_as_plotted", "linear_as_plotted", None),
               ("milne", "milne", None), ("eds", "frw", Universe(C.EDS.but(H0=70.0)))]}
    # The plotted error bars give chi2 = 83 for their LCDM, not their 38.9: find the
    # error scaling that reproduces 38.9 (their errors must have been ~1.46x larger).
    from scipy.optimize import brentq

    chi_l = lambda f: float(np.sum(((mobs - shapes["lcdm_0.3_0.7"] - curly_m) / (f * err)) ** 2))
    scale = float(brentq(lambda f: chi_l(f) - C.DL2004_TEXT["chi2_lcdm"][0], 0.5, 5.0))
    rows = {}
    for name, shape in shapes.items():
        rec = {}
        for tag, e in (("plotted_err", err), ("scaled_err", scale * err)):
            w = 1 / e**2
            off = float(np.sum(w * (mobs - shape)) / np.sum(w))
            rec[f"chi2_fixedM_{tag}"] = float(np.sum(((mobs - shape - curly_m) / e) ** 2))
            rec[f"chi2_freeM_{tag}"] = float(np.sum(((mobs - shape - off) / e) ** 2))
            rec[f"best_M_{tag}"] = off
        rows[name] = rec
    for name in rows:
        for key in [k for k in rows[name] if k.startswith("chi2_")]:
            dchi = rows[name][key] - rows["lcdm_0.3_0.7"][key]
            rows[name]["d" + key] = dchi
            rows[name]["sigma_" + key[5:]] = float(np.sqrt(max(dchi, 0)))
    # Their plotted curves against our formulas (same M = -3.45): which recipe they drew.
    curve_check = {}
    for style, c in zip(f5["curve_styles"], f5["curves"]):
        c = c[(c[:, 0] > 0.05) & (c[:, 0] < 0.99)]
        best = None
        for name, model, uu in [("SR (text: D_L = D(1+z))", "SR", None),
                                ("SR (D_L = D)", "SR_as_plotted", None),
                                ("v = cz (text: D_L = cz(1+z)/H)", "linear", None),
                                ("v = cz (D_L = cz/H)", "linear_as_plotted", None),
                                ("LCDM (0.3, 0.7)", "frw", u_dl),
                                ("Einstein-de Sitter", "frw", Universe(C.EDS.but(H0=70.0))),
                                ("de Sitter", "frw", Universe(C.DE_SITTER.but(H0=70.0)))]:
            d = c[:, 1] - (5 * np.log10(C.C_KM_S * hubble_free_dl(model, c[:, 0], uu)) + curly_m)
            mx = float(np.max(np.abs(d)))
            if best is None or mx < best[1]:
                best = (name, mx)
        curve_check[f"{style}_{len(curve_check)}"] = {"style": style, "best_match": best[0],
                                                      "max_abs_mag": best[1]}
    out["dl2004_reproduction"] = {"n": int(z.size), "rows": rows, "curve_check": curve_check,
                                  "error_scale_for_their_chi2": scale,
                                  "points": {"z": z.tolist(), "m": mobs.tolist(), "err": err.tolist()},
                                  "theirs": {k: C.DL2004_TEXT[k][0] for k in
                                             ("chi2_lcdm", "sr_dchi2", "linear_dchi2", "sr_sigma",
                                              "linear_sigma")}}
    # (b) Pantheon+: 1,590 light curves with the full covariance, offset profiled.
    t, cov, rws = sn_sample()
    zhd, zhel, m = t["zHD"][rws], t["zHEL"][rws], t["m_b_corr"][rws]
    cf = cho_factor(cov[np.ix_(rws, rws)])
    corr = (1 + zhel) / (1 + zhd)  # D_L = (1+z_hel) D_M(z_HD), the Pantheon+ convention
    models = {}
    grid = np.linspace(0.20, 0.45, 251)
    chis = []
    for om in grid:
        uu = Universe(C.PLANCK18.but(Om=om, T_cmb=0.0, N_eff=0.0, m_nu=()), n_cells=1500)
        chis.append(_chi2_profile(m, cf, 5 * np.log10(uu.H0 * uu.dl_of_z(zhd) * corr))[0])
    chis = np.array(chis)
    k = int(np.argmin(chis))
    a, b, c_ = np.polyfit(grid[k - 3:k + 4], chis[k - 3:k + 4], 2)
    om_best = float(-b / (2 * a))
    u_best = Universe(C.PLANCK18.but(Om=om_best, T_cmb=0.0, N_eff=0.0, m_nu=()))
    ref = {"flat_lcdm_best": (u_best, "frw"), "planck18": (universes["planck18"], "frw"),
           "milne": (None, "milne"), "SR": (None, "SR"), "SR_as_plotted": (None, "SR_as_plotted"),
           "linear": (None, "linear"), "linear_as_plotted": (None, "linear_as_plotted"),
           "eds": (Universe(C.EDS), "frw"), "de_sitter": (Universe(C.DE_SITTER), "frw")}
    for name, (uu, model) in ref.items():
        shape = 5 * np.log10(hubble_free_dl(model, zhd, uu) * corr)
        chi2, off = _chi2_profile(m, cf, shape)
        models[name] = {"chi2": chi2, "offset": off}
    base = models["flat_lcdm_best"]["chi2"]
    for name in models:
        models[name]["dchi2"] = models[name]["chi2"] - base
        models[name]["sigma"] = float(np.sqrt(max(models[name]["dchi2"], 0)))
    # Binned residuals for the figure: GLS bin means relative to the Milne universe.
    bins = np.array([0.01, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.85, 1.0,
                     1.3, 2.3])
    milne = 5 * np.log10(hubble_free_dl("milne", zhd) * corr) + models["milne"]["offset"]
    resid = m - milne
    zb = np.clip(np.digitize(zhd, bins) - 1, 0, bins.size - 2)
    W = np.zeros((zhd.size, bins.size - 1))
    W[np.arange(zhd.size), zb] = 1
    F = W.T @ cho_solve(cf, W)
    Fi = np.linalg.inv(F)
    dm = Fi @ (W.T @ cho_solve(cf, resid))
    var = np.diag(cov[np.ix_(rws, rws)])  # quote each bin at its inverse-variance mean redshift
    zeff = np.array([np.sum(zhd[zb == j] / var[zb == j]) / np.sum(1 / var[zb == j])
                     for j in range(bins.size - 1)])
    out["pantheon"] = {"n": int(rws.size), "Om_best_flat": om_best, "models": models,
                       "binned_vs_milne": {"z_eff": zeff.tolist(), "dmag": dm.tolist(),
                                           "sigma": np.sqrt(np.diag(Fi)).tolist(),
                                           "n": np.bincount(zb, minlength=bins.size - 1).tolist()},
                       "d1_check": "05-supernova-cosmology: empty-universe chi2 1451.44, flat best 1402.92"}
    return out


# ---------------------------------------------------------------------------
# 9. What we can still reach
# ---------------------------------------------------------------------------


def reachability(u: Universe) -> dict:
    eta_half = 0.5 * u.eta_inf
    t_eq = float(u.t_of_eta(eta_half))
    tt = np.geomspace(1e-3, 300.0, 2000)
    frac = np.minimum(1.0, ((u.eta_inf - u.eta_of_t(tt)) / u.eta_of_t(tt)) ** 3)
    out = {"t_horizons_equal_gyr": t_eq, "z_horizons_equal": float(1 / u.a_of_eta(eta_half) - 1),
           "chi_horizons_equal_gly": eta_half,
           "fraction_now": float(((u.eta_inf - u.eta0) / u.eta0) ** 3)}
    for f in (0.5, 0.1, 0.01, 0.001):
        out[f"t_fraction_{f}_gyr"] = float(np.interp(-f, -frac, tt))
    # The event horizon shrinks in comoving terms: by how much since the light of
    # the Pantheon+ supernovae, say, set out? Everything inside 16.6 Gly now is
    # reachable; a galaxy at z = 1.852 crosses out today.
    out["eh_comoving_at"] = {f"{t:g}_gyr": float(u.eta_inf - u.eta_of_t(t)) for t in (1, 5, 10, 13.787, 20, 50, 100)}
    return out


# ---------------------------------------------------------------------------
# 10. Sensitivity to the cosmology
# ---------------------------------------------------------------------------


def sensitivity() -> list:
    rows = []
    for p in [C.PLANCK18, C.DL2004, C.DL2004_RAD, *C.OTHER_UNIVERSES]:
        uu = Universe(p)
        s = uu.summary()
        s["source"] = p.source
        rows.append(s)
    return rows


# ---------------------------------------------------------------------------


def run_core() -> dict:
    warnings.filterwarnings("ignore")
    t0 = time.time()
    u = Universe(C.PLANCK18)
    universes = {"planck18": u, "dl2004": Universe(C.DL2004), "dl2004_rad": Universe(C.DL2004_RAD)}
    print("1. expansion history vs astropy and CAMB", flush=True)
    checks = cosmology_checks(u)
    camb_d = checks.get("camb_derived") or {}
    z_star = camb_d.get("zstar", C.PLANCK18_PUBLISHED["z_star"][0])
    r_drag = camb_d.get("rdrag", C.PLANCK18_PUBLISHED["r_drag_mpc"][0])
    print("2. headline numbers", flush=True)
    head = headline(u, z_star)
    head["planck_checks"]["z_star"] = {"planck": C.PLANCK18_PUBLISHED["z_star"], "camb": z_star}
    head["planck_checks"]["r_drag_mpc"] = {"planck": C.PLANCK18_PUBLISHED["r_drag_mpc"], "camb": r_drag}
    print("3. Davis & Lineweaver 2004", flush=True)
    dl = dl2004_check({"dl2004": universes["dl2004"], "dl2004_rad": universes["dl2004_rad"]})
    print("4. uncertainties", flush=True)
    unc = uncertainties()
    print("5. real objects on the light cone", flush=True)
    pops = catalogs.all_populations(u)
    marks = catalogs.landmarks(pops["jades"])
    cone, arrays = objects_on_cone(u, pops, marks, z_star)
    gw = pops["gw_events"]
    cone["gw_redshift_vs_catalogue"] = {"median_ratio": float(np.median(gw.z / gw.z_catalogue)),
                                        "max_abs_diff": float(np.max(np.abs(gw.z - gw.z_catalogue)))}
    print("6. DESI", flush=True)
    desi = desi_measurements(u, r_drag)
    print("7. supernova recession speeds", flush=True)
    snv = sn_measurements(u)
    print("8. special relativity vs expansion: supernova test", flush=True)
    sr = sr_test(universes)
    print("9. reachability", flush=True)
    reach = reachability(u)
    print("10. sensitivity", flush=True)
    sens = sensitivity()
    results = {"_generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "cosmology_checks": checks, "headline": head, "dl2004": dl, "uncertainties": unc,
               "cone": cone, "desi": desi, "sn_velocities": snv, "sr_test": sr,
               "reachability": reach, "sensitivity": sens, "seconds": time.time() - t0}
    C.RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    C.RESULTS_PATH.write_text(json.dumps(_f(results), indent=1))
    OBJECTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OBJECTS_PATH, **arrays)
    print(f"done in {time.time() - t0:.0f} s -> {C.RESULTS_PATH.relative_to(C.HERE)}", flush=True)
    return json.loads(C.RESULTS_PATH.read_text())


if __name__ == "__main__":
    run_core()
