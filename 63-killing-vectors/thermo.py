"""Black-hole thermodynamics for every black hole this atlas has weighed.

A horizon is a **Killing horizon**: the surface where a Killing vector (for
Kerr, xi = d_t + Omega_H d_phi) turns from timelike to null. Its **surface
gravity** kappa says how fast that Killing vector's length shrinks towards the
horizon, i.e. the force per unit mass needed, from far away, to hold a body just
outside. Hawking (1975) showed that quantum fields make the horizon glow at

    T_H = hbar kappa / (2 pi c k_B)   (6.17e-8 K for one solar mass, no spin),

and Bekenstein's entropy is a quarter of the horizon area in Planck units,

    S = k_B A c^3 / (4 G hbar)        (1.05e77 k_B for one solar mass).

Hawking's area theorem (1971) says A never decreases in classical processes;
read as entropy, it is the second law. For Kerr (G = c = M = 1, spin chi):

    r+ = 1 + sqrt(1 - chi^2),  A = 8 pi r+,  kappa = sqrt(1 - chi^2) / (2 r+),
    Omega_H = chi / (2 r+),    M_irr = sqrt(A / 16 pi),

and M - M_irr is rotational energy that the Penrose process (or magnetic
fields) could in principle extract: up to 29% of the mass for chi = 1.

None of this radiation has ever been observed: every astrophysical black hole
is far colder than the 2.7 K microwave background, so today each one absorbs
more than it emits. The numbers are predictions.
"""

from __future__ import annotations

import json

import numpy as np
from scipy import integrate

import config as C
import literature as lit

# ---------------------------------------------------------------------------
# Kerr horizon quantities in SI units
# ---------------------------------------------------------------------------


def kerr_thermo(M_msun, chi=0.0) -> dict:
    """Horizon thermodynamics for mass M (solar masses) and spin chi."""
    M = np.asarray(M_msun, float) * C.M_SUN
    chi = np.asarray(chi, float)
    s = np.sqrt(1 - chi**2)
    rg = C.G * M / C.C**2                     # GM/c^2, m
    rp = rg * (1 + s)
    area = 8 * np.pi * rg**2 * (1 + s)
    kappa_hat = s / (2 * (1 + s))             # surface gravity in units of c^4/(GM)
    kappa = kappa_hat * C.C**4 / (C.G * M)    # m/s^2
    T = C.HBAR * kappa / (2 * np.pi * C.C * C.K_B)
    S = area * C.C**3 / (4 * C.G * C.HBAR)    # in units of k_B
    omega = chi / (2 * (1 + s)) * C.C**3 / (C.G * M)
    m_irr = np.sqrt((1 + s) / 2)              # M_irr / M
    return {"r_plus_km": rp / 1e3, "area_km2": area / 1e6, "kappa_ms2": kappa, "T_K": T, "S_kB": S,
            "Omega_H_rad_s": omega, "spin_rate_Hz": omega / (2 * np.pi), "M_irr_frac": m_irr,
            "E_rot_frac": 1 - m_irr}


def evaporation_time_yr(M_msun, alpha=lit.PAGE["alpha_photons_gravitons"]):
    """Lifetime G^2 M^3 / (3 alpha hbar c^4) of a non-spinning hole in empty space (Page 1976)."""
    M = np.asarray(M_msun, float) * C.M_SUN
    return C.G**2 * M**3 / (3 * alpha * C.HBAR * C.C**4) / C.YEAR


def blackbody_evaporation_time_yr(M_msun):
    """Textbook photon-only estimate 5120 pi G^2 M^3 / (hbar c^4)."""
    M = np.asarray(M_msun, float) * C.M_SUN
    return 5120 * np.pi * C.G**2 * M**3 / (C.HBAR * C.C**4) / C.YEAR


def mass_at_temperature_msun(T_K):
    """Non-spinning mass whose Hawking temperature is T."""
    return C.HBAR * C.C**3 / (8 * np.pi * C.G * C.K_B * np.asarray(T_K, float)) / C.M_SUN


# ---------------------------------------------------------------------------
# GW250114: the area theorem read as an entropy increase
# ---------------------------------------------------------------------------


def gw250114(path=None) -> dict:
    """Horizon areas, entropies and temperatures from the LVK posterior samples.

    Source-frame masses; component spins a_1, a_2 are the magnitudes (the area
    depends only on the magnitude). The samples come from waveform models that
    assume GR, so the area increase is a consequence of the fit, not a test of
    it: the test (inspiral and ringdown analysed separately) is project I7's.
    """
    import h5py

    path = path or C.RAW_DIR / C.GW250114["samples_file"]
    with h5py.File(path, "r") as f:
        ps = f["bilby-NRSur7dq4_prod-reweighted/posterior_samples"][:]
    m1, m2 = ps["mass_1_source"], ps["mass_2_source"]
    a1, a2 = ps["a_1"], ps["a_2"]
    mf, chif = ps["final_mass_source"], ps["final_spin"]
    t1, t2, tf = kerr_thermo(m1, a1), kerr_thermo(m2, a2), kerr_thermo(mf, chif)
    A_i = t1["area_km2"] + t2["area_km2"]
    A_f = tf["area_km2"]
    S_i = t1["S_kB"] + t2["S_kB"]
    S_f = tf["S_kB"]
    # The largest energy the merger could have radiated with the final area equal to the
    # initial one (Hawking 1971), for the measured final spin.
    mf_min = np.sqrt(A_i * 1e6 / (8 * np.pi * (1 + np.sqrt(1 - chif**2)))) * C.C**2 / C.G / C.M_SUN
    erad = m1 + m2 - mf

    def q(x):
        lo, med, hi = np.percentile(x, [5, 50, 95])
        return {"median": float(med), "plus90": float(hi - med), "minus90": float(med - lo)}

    return {
        "n_samples": int(len(m1)),
        "m1": q(m1), "m2": q(m2), "a1": q(a1), "a2": q(a2), "M_final": q(mf), "chi_final": q(chif),
        "area_initial_km2": q(A_i), "area_final_km2": q(A_f), "area_ratio": q(A_f / A_i),
        "area_increase_fraction_positive": float(np.mean(A_f > A_i)),
        "S_initial_kB": q(S_i), "S_final_kB": q(S_f), "delta_S_kB": q(S_f - S_i),
        "T_final_K": q(tf["T_K"]), "T_final_over_schwarzschild": q(tf["T_K"] / kerr_thermo(mf, 0.0)["T_K"]),
        "horizon_spin_rate_Hz": q(tf["spin_rate_Hz"]), "E_rot_final_frac": q(tf["E_rot_frac"]),
        "E_radiated_msun": q(erad), "E_radiated_frac": q(erad / (m1 + m2)),
        "max_radiated_by_area_theorem_frac": q((m1 + m2 - mf_min) / (m1 + m2)),
        "paper": lit.GW250114_PAPER,
    }


# ---------------------------------------------------------------------------
# The census: masses (and spins) from earlier projects
# ---------------------------------------------------------------------------


def _load_json(path):
    return json.loads(path.read_text()) if path.exists() else None


def census() -> dict:
    """Every black hole the atlas has weighed, with T_H, S, lifetime (and spin where known)."""
    spins = {s["i4_name"]: s for s in lit.XRB_SPINS if s["i4_name"]}
    rows, notes = [], []

    def add(name, group, M, M_lo=None, M_hi=None, chi=None, chi_note="", source=""):
        th = kerr_thermo(M, chi if chi is not None else 0.0)
        rows.append({"name": name, "group": group, "M_msun": float(M), "M_lo": M_lo, "M_hi": M_hi,
                     "chi": chi, "chi_note": chi_note, "T_K": float(th["T_K"]),
                     "T_K_nonspinning": float(kerr_thermo(M, 0.0)["T_K"]), "S_kB": float(th["S_kB"]),
                     "r_plus_km": float(th["r_plus_km"]), "E_rot_frac": float(th["E_rot_frac"]),
                     "t_evap_yr": float(evaporation_time_yr(M)), "source": source})

    # Gravitational-wave remnants
    gw = gw250114()
    add("GW250114 remnant", "merger remnant", gw["M_final"]["median"], chi=gw["chi_final"]["median"],
        chi_note="LVK posterior", source="LVK 2025 posterior samples (I7 not built)")
    c1 = _load_json(C.C1_RESULTS)
    if c1:
        cat = c1["published"]["catalogue"]
        add("GW150914 remnant", "merger remnant", cat["final_mass_source"]["value"],
            chi=cat["a_final"]["value"], chi_note="GWTC", source="C1 (02-gravitational-waves), GWTC values")
    # X-ray binaries and Gaia black holes (I4)
    i4 = _load_json(C.I4_RESULTS)
    if i4:
        w = i4["weighing"]
        for s in w["systems"]:
            if "m1_ours" not in s:
                continue
            m = s["m1_ours"]
            label = s["common_name"] or s["name"]
            sp_ = spins.get(s["name"])
            add(label, "X-ray binary", m["median"], m["minus"], m["plus"],
                chi=sp_["a"] if sp_ else None, chi_note=(sp_["method"] + ", " + sp_["kind"]) if sp_ else "",
                source="I4 (35-x-ray-binaries), our mass")
        for s in w["persistent"]:
            m = s["m_bh"]
            sp_ = spins.get(s["name"])
            add(s["name"], "X-ray binary", m["median"], m["minus"], m["plus"],
                chi=sp_["a"] if sp_ else None, chi_note=(sp_["method"] + ", " + sp_["kind"]) if sp_ else "",
                source="I4 (35-x-ray-binaries)")
        for s in w["gaia"]:
            m = s["m_bh"]
            add(s["name"], "dormant (Gaia)", m["median"], m["minus"], m["plus"], source="I4 (35-x-ray-binaries)")
    else:
        notes.append("I4 results missing: X-ray binaries left out")
    # Sgr A* (H1) and M87* (I5)
    h1 = _load_json(C.H1_RESULTS)
    if h1:
        g = h1["fits"]["joint_gr"]["parameters"]["gm6"]
        add("Sgr A*", "supermassive", g["median"] * 1e6, g["sd"] * 1e6, g["sd"] * 1e6,
            source="H1 (17-galactic-centre), joint GR fit")
    else:
        add("Sgr A*", "supermassive", C.FALLBACK["sgra_mass_1e6"][0] * 1e6, source="H1 fallback")
    i5 = _load_json(C.I5_RESULTS)
    if i5:
        m = i5["m87_mass"]["grmhd_calibrated"]["mass_1e9"]
        add("M87*", "supermassive", m["median"] * 1e9, (m["median"] - m["lo68"]) * 1e9,
            (m["hi68"] - m["median"]) * 1e9, source="I5 (36-m87-shadow), ring + GRMHD calibration")
    # Kormendy & Ho 2013 (I1): the heaviest and lightest dynamically measured
    kh = None
    if C.I1_KH13.exists():
        import pandas as pd
        kh = pd.read_csv(C.I1_KH13)
        kh = kh[(kh["M"] > 0) & (~kh["upper"].astype(bool))]
        for _, r_ in kh.sort_values("M").iloc[[0, -1]].iterrows():
            add(r_["name"], "supermassive", r_["M"], r_["M"] - r_["M_lo"], r_["M_hi"] - r_["M"],
                source="I1 (32-m-sigma-relation), Kormendy & Ho 2013")
    # Quasars (I3): Wu & Shen 2022 virial masses, as a distribution
    quasars = None
    if C.I3_WUSHEN.exists():
        import pyarrow.parquet as pq
        lm = pq.read_table(C.I3_WUSHEN, columns=["LOGMBH"]).column(0).to_numpy()
        lm = lm[lm > 0]
        quasars = {"n": int(lm.size), "logM_p05_p50_p95": [float(v_) for v_ in np.percentile(lm, [5, 50, 95])],
                   "logM_max": float(lm.max()), "hist_edges": np.linspace(6, 11, 51).tolist(),
                   "hist_counts": np.histogram(lm, bins=np.linspace(6, 11, 51))[0].tolist()}
    else:
        notes.append("I3 quasar masses missing")
    kh_all = None if kh is None else {"n": int(len(kh)), "M": kh["M"].tolist(), "names": kh["name"].tolist()}
    return {"rows": rows, "quasars": quasars, "kh13": kh_all, "notes": notes, "gw250114": gw}


# ---------------------------------------------------------------------------
# The far future: when does the microwave background get colder than a black hole?
# ---------------------------------------------------------------------------


def far_future(masses_msun: dict) -> dict:
    """Time from today until T_CMB = T_H, in Planck 2018 LCDM (radiation + matter + Lambda).

    T_CMB falls as 1/a; once Lambda dominates, a grows as e^(H_L t) with
    H_L = H0 sqrt(Omega_L). Also the de Sitter horizon's own (Gibbons-Hawking)
    temperature hbar H_L / (2 pi k_B): the floor below which nothing cools.
    """
    H0 = C.PLANCK18["H0"] * 1e3 / (1e6 * C.PC)  # 1/s
    Om, Or = C.PLANCK18["Omega_m"], C.PLANCK18["Omega_r"]
    OL = 1 - Om - Or
    HL = H0 * np.sqrt(OL)

    def years_until(a_target):
        f = lambda lna: 1 / np.sqrt(Or * np.exp(-4 * lna) + Om * np.exp(-3 * lna) + OL)
        val, _ = integrate.quad(f, 0, np.log(a_target), limit=500)
        return val / H0 / C.YEAR

    out = {}
    for name, M in masses_msun.items():
        T = float(kerr_thermo(M, 0.0)["T_K"])
        a_star = C.T_CMB / T
        out[name] = {"T_H_K": T, "a_when_cmb_equals_T_H": a_star, "years_from_now": years_until(a_star)}
    T_dS = C.HBAR * HL / (2 * np.pi * C.K_B)
    return {"black_holes": out, "H_Lambda_per_s": HL, "T_de_Sitter_K": T_dS,
            "mass_with_T_H_equal_T_dS_msun": float(mass_at_temperature_msun(T_dS)),
            "mass_with_T_H_equal_T_CMB_kg": float(mass_at_temperature_msun(C.T_CMB) * C.M_SUN)}
