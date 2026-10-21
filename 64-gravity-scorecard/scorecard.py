"""The scorecard: every test of gravity on one footing.

Each test becomes one row: the quantity whose general-relativity value is known
(a ratio measured/GR, or a violation parameter that GR sets to zero), the
measured *fractional deviation* from GR with a symmetric 1-sigma error (see
``literature.one_sigma`` for the conversions), and the strength of the gravity
it probes, GM/(r c^2) of the dominant field at the place the test happens
(Baker, Psaltis & Skordis 2015 use the same "potential" axis):

* laboratory and Earth-orbit tests: the Earth's field at the apparatus;
* solar-system tests: the Sun's field where the light passes or the planet orbits;
* binary pulsars: G(m1 + m2)/(a c^2) of the orbit; the triple pulsar: the outer
  white dwarf's field at the inner pair (the neutron star's own compactness,
  ~0.1-0.2, is what makes it a strong-field body; the axis shows the field it falls in);
* black holes: 1/3, the photon ring at r = 3GM/c^2; merging black holes: ~0.1-0.3;
* galaxies and cosmology: ~1e-6 (a galaxy within a few kpc) to ~1e-5 (clusters,
  the large-scale structure). Wide binaries: two suns 5,000 AU apart.

Published results come from the checked tables in ``literature/checks``; this
atlas's own measurements are read from the sibling projects' ``results.json``
(or, where a number lives only in a project's README, transcribed with its path).
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

import config as C
import literature as L
from gr_predictions import GM_EARTH, GM_SUN, M_SUN_M, R_EARTH, R_SUN

C2 = C.C_LIGHT**2


def phi_earth(r_m: float) -> float:
    return GM_EARTH / (r_m * C2)


def phi_sun(r_m: float) -> float:
    return GM_SUN / (r_m * C2)


def phi_mass(m_sun: float, r_m: float) -> float:
    return m_sun * M_SUN_M / r_m


def orbit_a(m_tot: float, p_days: float) -> float:
    """Relative semi-major axis (m) from Kepler's third law."""
    p = p_days * C.DAY
    return (GM_SUN * m_tot * p**2 / (4 * np.pi**2)) ** (1 / 3)


def binary_phi(m_tot: float, p_days: float) -> float:
    return phi_mass(m_tot, orbit_a(m_tot, p_days))


X = {
    "lab": phi_earth(R_EARTH),                                    # 6.96e-10
    "microscope": phi_earth(R_EARTH + 710e3),
    "gpb": phi_earth(R_EARTH + 642e3),
    "lares": phi_earth(12_270_020.0),
    "galileo": phi_earth(27_977.8e3),
    "gpa": phi_earth(R_EARTH + 10_000e3),
    "sun_limb_cassini": phi_sun(1.6 * R_SUN),
    "sun_vlbi": phi_sun(np.radians(5.0) * C.AU),                  # sources ~5 degrees from the Sun
    "sun_surface": phi_sun(R_SUN),
    "mercury": phi_sun(0.38709927 * C.AU),
    "earth_orbit": phi_sun(C.AU),
    "hulse_taylor": binary_phi(1.438 + 1.390, 0.322997448918),
    "double_pulsar": binary_phi(2.587052, 0.1022515592973),
    "j0348": binary_phi(2.01 + 0.172, 0.102424062722),
    "j1141": binary_phi(2.2892, 0.1976509593),
    "triple_outer": phi_mass(0.4101, orbit_a(1.4378 + 0.19751 + 0.4101, 327.257541)),
    "s2_peri": phi_mass(4.3e6, 120 * C.AU),
    "sirius_b": None,                                             # filled from N3's results
    "photon_ring": 1 / 3,
    "inspiral": 0.1,
    "merger": 0.2,
    "galaxy": phi_mass(1e11, 2e3 * 3.0856775814913673e16),          # Collett et al. 2018: ~2 kpc
    "cosmology": 1e-5,
    "wide_binary": phi_mass(1.5, 5_000 * C.AU),
}


def lit(name, rid, *, ratio_to=None, minus_one=False, scale=1.0, treat_as=None):
    """(deviation, sigma) from a checked row: ratio_to divides by GR's value,
    minus_one subtracts 1 from a ratio, scale multiplies (e.g. 1/2 for gamma -> (1+gamma)/2)."""
    v, e = L.one_sigma(name, rid, treat_as=treat_as)
    if ratio_to is not None:
        v, e = v / ratio_to, e / abs(ratio_to)
    if minus_one:
        v = v - 1.0
    return v * scale, e * abs(scale)


def literature_rows() -> list[dict]:
    R = []

    def add(family, label, year, ref, dev_sig, x, quantity, *, short=None, note="", source=""):
        dev, sig = dev_sig
        R.append({"family": family, "label": label, "short": short or label, "year": year, "ref": ref,
                  "deviation": dev, "sigma": sig, "x": x, "quantity": quantity, "atlas": "",
                  "note": note, "source": source})

    # --- equivalence principle -------------------------------------------------
    add("equivalence", "MICROSCOPE (Ti vs Pt in orbit)", 2022, "Touboul et al. 2022",
        L.one_sigma("wep", "microscope2022"), X["microscope"], "Eötvös ratio η", short="MICROSCOPE",
        source="wep.csv:microscope2022")
    add("equivalence", "Eöt-Wash torsion balance (Be vs Al)", 2012, "Wagner et al. 2012",
        L.one_sigma("wep", "wagner2012_BeAl_earth"), X["lab"], "Eötvös ratio η", short="Eöt-Wash",
        source="wep.csv:wagner2012_BeAl_earth")
    add("equivalence", "Atom interferometer (⁸⁵Rb vs ⁸⁷Rb)", 2020, "Asenbaum et al. 2020",
        L.one_sigma("wep", "asenbaum2020"), X["lab"], "Eötvös ratio η", short="atoms",
        source="wep.csv:asenbaum2020")
    v, e = L.one_sigma("wep", "alphag2023")
    add("equivalence", "Antihydrogen falls (ALPHA-g)", 2023, "Anderson et al. 2023", (v - 1, e), X["lab"],
        "a(antihydrogen)/g − 1", short="antimatter", source="wep.csv:alphag2023")
    add("equivalence", "Earth and Moon fall alike (lunar ranging)", 2021, "Biskupek et al. 2021",
        L.one_sigma("wep", "llr_biskupek2021"), X["earth_orbit"], "(m_G/m_I)_Earth − (m_G/m_I)_Moon",
        short="LLR Earth–Moon", note="total: composition + self-energy", source="wep.csv:llr_biskupek2021")
    v, e = L.one_sigma("pulsars", "tp_v25_Delta_PL3")
    add("equivalence", "Neutron star falls like a white dwarf (triple pulsar)", 2025, "Voisin et al. 2025",
        (v, e), X["triple_outer"], "Δ = (m_G/m_I)_NS − (m_G/m_I)_WD", short="triple pulsar",
        note="red-noise model; 'planet' model gives (0.51 ± 0.62)e-6", source="pulsars.csv:tp_v25_Delta_PL3")

    # --- weak field: the solar system and Earth orbit ---------------------------
    add("weak", "Shapiro delay, Cassini", 2003, "Bertotti, Iess & Tortora 2003",
        lit("beta_lunar", "PPN-B03", scale=0.5), X["sun_limb_cassini"], "(1+γ)/2 − 1", short="Cassini",
        source="beta_lunar.csv:PPN-B03")
    vlbi = gamma_best_vlbi()
    if vlbi is not None:
        add("weak", f"Light bending, VLBI ({vlbi['ref']})", vlbi["year"], vlbi["ref"], (vlbi["dev"], vlbi["sig"]),
            X["sun_vlbi"], "(1+γ)/2 − 1", short="VLBI deflection", source=vlbi["source"])
    p = L.ppn_constraints()
    v, e = p["perihelion_2g_minus_b"]
    add("weak", "Mercury's perihelion (MESSENGER)", 2017, "Park et al. 2017", (v / 3, e / 3), X["mercury"],
        "(2 + 2γ − β)/3 − 1", short="Mercury perihelion", note="β from Park et al. with Cassini's γ",
        source="beta_lunar.csv:BETA-P17,PPN-B03")
    add("weak", "Clock redshift, Galileo satellites", 2018, "Delva et al. 2018", (0.19e-5, 2.48e-5), X["galileo"],
        "redshift α", short="Galileo clocks", source="N3 ladder (58-gravitational-redshift/ladder.py)")
    add("weak", "Clock redshift, Gravity Probe A", 1980, "Vessot et al. 1980", (5e-6, 1.4e-4), X["gpa"],
        "redshift α", short="GP-A", source="N3 ladder")
    add("weak", "Redshift in a tower (Pound & Snider)", 1965, "Pound & Snider 1965", (-0.0010, 0.0076), X["lab"],
        "measured/GR − 1", short="Pound–Snider", source="N3 ladder")
    add("weak", "Geodetic precession, Gravity Probe B", 2011, "Everitt et al. 2011",
        lit("strongfield", "gpb2011_geodetic", ratio_to=-6606.1, minus_one=True), X["gpb"], "measured/GR − 1",
        short="GP-B geodetic", source="strongfield.csv:gpb2011_geodetic")
    add("weak", "Frame dragging, Gravity Probe B", 2011, "Everitt et al. 2011",
        lit("strongfield", "gpb2011_framedrag", ratio_to=-39.2, minus_one=True), X["gpb"], "measured/GR − 1",
        short="GP-B frame dragging", source="strongfield.csv:gpb2011_framedrag")
    add("weak", "Frame dragging, LAGEOS + LARES 2", 2026, "Ciufolini et al. 2026",
        lit("strongfield", "lares2_2026_nature", minus_one=True), X["lares"], "μ − 1", short="LARES 2",
        note="error budget disputed (Iorio 2025; a Nature referee)", source="strongfield.csv:lares2_2026_nature")
    add("weak", "Geodetic precession of the Moon's orbit (LLR)", 2004, "Williams, Turyshev & Boggs 2004",
        L.one_sigma("beta_lunar", "GEO-WTB04"), X["earth_orbit"], "K_gp", short="LLR geodetic",
        source="beta_lunar.csv:GEO-WTB04")
    add("weak", "Sun's gravitational redshift", 2020, "González Hernández et al. 2020", (638 / 633.1 - 1, 6 / 633.1),
        X["sun_surface"], "measured/GR − 1", short="Sun redshift", source="N3 ladder")

    # --- strong field and dynamics ---------------------------------------------
    add("strong", "Hulse–Taylor pulsar: orbital decay", 2016, "Weisberg & Huang 2016",
        lit("pulsars", "b1913_wh16_pbdot_ratio", minus_one=True), X["hulse_taylor"], "Ṗb/Ṗb(GR) − 1",
        short="Hulse–Taylor Ṗb", note="limited by the Galactic correction (distance)", source="pulsars.csv:b1913_wh16_pbdot_ratio")
    add("strong", "Double pulsar: orbital decay", 2021, "Kramer et al. 2021",
        lit("pulsars", "dp_k21_test_pbdot", minus_one=True), X["double_pulsar"], "Ṗb/Ṗb(GR) − 1",
        short="double pulsar Ṗb", source="pulsars.csv:dp_k21_test_pbdot")
    add("strong", "Double pulsar: Shapiro delay shape s", 2021, "Kramer et al. 2021",
        lit("pulsars", "dp_k21_test_s", minus_one=True), X["double_pulsar"], "s/s(GR) − 1",
        short="double pulsar s", source="pulsars.csv:dp_k21_test_s")
    add("strong", "Double pulsar: spin precession of B", 2008, "Breton et al. 2008",
        lit("pulsars", "dp_k21_test_OmB", minus_one=True), X["double_pulsar"], "Ω/Ω(GR) − 1",
        short="spin precession", source="pulsars.csv:dp_k21_test_OmB")
    add("strong", "Pulsar–white dwarf J0348+0432: orbital decay", 2013, "Antoniadis et al. 2013",
        lit("pulsars", "j0348_a13_ratio", minus_one=True), X["j0348"], "Ṗb/Ṗb(GR) − 1",
        short="J0348+0432", source="pulsars.csv:j0348_a13_ratio")
    add("strong", "S2's redshift at the Galactic Centre", 2019, "GRAVITY Collaboration 2019",
        lit("strongfield", "gravity2019_f_redshift", minus_one=True), X["s2_peri"], "f − 1", short="S2 redshift",
        source="strongfield.csv:gravity2019_f_redshift")
    add("strong", "S2's orbit precesses", 2024, "GRAVITY Collaboration 2024",
        lit("strongfield", "gravity2024_fsp", minus_one=True), X["s2_peri"], "f_SP − 1", short="S2 precession",
        source="strongfield.csv:gravity2024_fsp")
    add("strong", "M87* shadow size", 2019, "EHT Collaboration 2019 (Paper VI)",
        L.one_sigma("strongfield", "eht_m87_2019_delta_stellar"), X["photon_ring"], "δ (shadow deviation)",
        short="M87* shadow", note="needs the stellar-dynamics mass", source="strongfield.csv:eht_m87_2019_delta_stellar")
    add("strong", "Sgr A* shadow size", 2022, "EHT Collaboration 2022 (Paper VI)",
        L.one_sigma("strongfield", "eht_sgra_2022_delta_vlti"), X["photon_ring"], "δ (shadow deviation)",
        short="Sgr A* shadow", source="strongfield.csv:eht_sgra_2022_delta_vlti")
    add("strong", "Ringdown frequency, GW250114", 2026, "LVK 2026",
        L.one_sigma("strongfield", "gw250114_df220"), X["photon_ring"], "δf̂₂₂₀", short="GW250114 ringdown",
        source="strongfield.csv:gw250114_df220")
    add("strong", "Inspiral phase (0PN), GW250114", 2026, "LVK 2026",
        L.one_sigma("strongfield", "gw250114_phi0"), X["inspiral"], "δφ̂₀", short="GW250114 inspiral",
        source="strongfield.csv:gw250114_phi0")
    add("strong", "Inspiral–merger–ringdown consistency, GWTC-4", 2026, "LVK 2026 (GWTC-4.0 tests)",
        L.one_sigma("strongfield", "gwtc4_imr_dMf_hier"), X["merger"], "ΔM_f/M̄_f", short="IMR consistency",
        source="strongfield.csv:gwtc4_imr_dMf_hier")
    add("strong", "Ringdown damping time, GWTC-5 catalogue", 2026, "LVK 2026 (GWTC-5.0 tests)",
        L.one_sigma("strongfield", "gwtc5_ringdown_dtau220_hier"), X["photon_ring"], "δτ̂₂₂₀",
        short="ringdown damping", note="combined catalogue sits slightly off GR; LVK suspect systematics",
        source="strongfield.csv:gwtc5_ringdown_dtau220_hier")

    # --- galaxies and cosmology ------------------------------------------------
    add("cosmology", "Lensing vs dynamics of a galaxy (γ at 2 kpc)", 2018, "Collett et al. 2018",
        lit("strongfield", "collett2018_gamma", minus_one=True, scale=0.5), X["galaxy"], "(1+γ)/2 − 1",
        short="galaxy γ", source="strongfield.csv:collett2018_gamma")
    add("cosmology", "Lensing vs galaxy motions, E_G", 2010, "Reyes et al. 2010",
        lit("strongfield", "reyes2010_EG", ratio_to=0.408, minus_one=True), X["cosmology"], "E_G/E_G(GR) − 1",
        short="E_G", source="strongfield.csv:reyes2010_EG")
    add("cosmology", "Light deflection by large-scale structure, Σ0 (DESI)", 2025, "Ishak et al. 2025",
        L.one_sigma("strongfield", "desi2024_Sigma0"), X["cosmology"], "Σ0", short="Σ0 (DESI)",
        source="strongfield.csv:desi2024_Sigma0")
    add("cosmology", "Growth of structure, μ0 (DESI)", 2025, "Ishak et al. 2025",
        L.one_sigma("strongfield", "desi2024_mu0"), X["cosmology"], "μ0", short="μ0 (DESI)",
        source="strongfield.csv:desi2024_mu0")
    v, e = L.one_sigma("strongfield", "nguyen2023_growth_index")
    add("cosmology", "Growth index of structure", 2023, "Nguyen, Huterer & Wen 2023", ((v - 0.55) / 0.55, e / 0.55),
        X["cosmology"], "γ_growth/0.55 − 1", short="growth index",
        note="3.4σ; relies on Planck PR3 lensing (see strongfield_notes.md)", source="strongfield.csv:nguyen2023_growth_index")
    return R


def gamma_best_vlbi() -> dict | None:
    """Light bending from geodetic VLBI: Lambert & Le Poncin-Lafitte (2011), gamma = 0.99992 +- 0.00012
    (formal). Titov et al. (2018)'s global solution has a smaller error but sits 3 sigma off, which its
    authors attribute to systematics (gamma_notes.md), so it is shown only on the gamma-history figure."""
    if L.gamma_table() is None:
        return None
    v, e = L.one_sigma("gamma", "V2011-lambert")
    return {"dev": (v - 1) / 2, "sig": e / 2, "year": 2011, "ref": "Lambert & Le Poncin-Lafitte 2011",
            "source": "gamma.csv:V2011-lambert"}


# ---------------------------------------------------------------------------
# This atlas's own measurements
# ---------------------------------------------------------------------------

def _load(rel: str) -> dict | None:
    p = C.ATLAS_DIR / rel
    return json.loads(p.read_text()) if p.exists() else None


def atlas_rows(llr: dict | None = None) -> list[dict]:
    R = []

    def add(project, family, label, dev, sig, x, quantity, short, source, note=""):
        R.append({"family": family, "label": label, "short": short, "year": 2026, "ref": f"this atlas, {project}",
                  "deviation": dev, "sigma": sig, "x": x, "quantity": quantity, "atlas": project, "note": note,
                  "source": source})

    n3 = _load("58-gravitational-redshift/data/results.json")
    if n3:
        g = n3["galileo"]["combined"]
        add("N3", "weak", "Galileo clocks' redshift (atlas)", g["alpha"], g["total"], X["galileo"], "redshift α",
            "N3 Galileo", "58-gravitational-redshift/data/results.json:galileo.combined")
        s = n3["sun"]
        add("N3", "weak", "Sun's redshift in iron lines (atlas)", s["strong_mean"] / s["theory_earth"] - 1,
            s["strong_total_err"] / s["theory_earth"], X["sun_surface"], "measured/GR − 1", "N3 Sun",
            "58-gravitational-redshift/data/results.json:sun")
        sb = n3["sirius_b"]
        X["sirius_b"] = sb["potential_over_c2"]
        add("N3", "strong", "Sirius B's redshift (atlas)", sb["ratio"] - 1, sb["ratio_err"], sb["potential_over_c2"],
            "measured/GR − 1", "N3 Sirius B", "58-gravitational-redshift/data/results.json:sirius_b")
        add("N3", "weak", "J2 relativistic clock term (atlas)", 0.02, 0.04, X["galileo"], "measured/GR − 1",
            "N3 J2 term", "58-gravitational-redshift/README.md ('1.02 ± 0.04 of GR's size')")
    h1 = _load("17-galactic-centre/data/results.json")
    if h1:
        k = h1["relativity"]["keck_redshift"]
        add("H1", "strong", "S2's redshift, Keck positions (atlas)", k["f_median"] - 1, k["f_sd"], X["s2_peri"], "f − 1",
            "H1 S2 redshift", "17-galactic-centre/data/results.json:relativity.keck_redshift",
            note="0.70 ± 0.15 adding VLT positions")
    i5 = _load("36-m87-shadow/data/results.json")
    if i5:
        # Sgr A*: ring 53.9 +- 1.3 uas vs the EHT-calibrated ring for H1's mass (53.14 uas); add the
        # EHT's ~11% ring-to-shadow calibration spread (I5 README), which the 2.5% leaves out
        st = i5["sgra_test"]
        ring, ring_err = 53.9, 1.3
        dev = ring / st["ring_prediction_uas"] - 1
        sig = float(np.hypot(ring_err / st["ring_prediction_uas"], 0.11))
        add("I5", "strong", "Sgr A*'s ring vs S2's mass (atlas)", dev, sig, X["photon_ring"], "ring/prediction − 1",
            "I5 Sgr A*", "36-m87-shadow/README.md + results.json:sgra_test",
            note="±2.5% measurement ⊕ 11% EHT calibration")
        m = i5["m87_mass"]["grmhd_calibrated"]["theta_g"]
        th, th_err = m["median"], m["sd"]
        th_star, th_star_err = 3.62, 3.62 * 0.4 / 6.6          # Gebhardt et al. 2011 (I5 predictions)
        add("I5", "strong", "M87*'s ring vs stellar-dynamics mass (atlas)", th / th_star - 1,
            (th / th_star) * float(np.hypot(th_err / th, th_star_err / th_star)), X["photon_ring"],
            "θ_g(ring)/θ_g(stars) − 1", "I5 M87*", "36-m87-shadow/data/results.json:m87_mass")
    m3 = _load("53-wide-binaries/data/results.json")
    if m3:
        a = m3["variants"]["ecc_hwang"]["alpha"]
        add("M3", "cosmology", "Wide binaries: Newton vs MOND (atlas)", a["median"], 0.5 * (a["hi"] - a["lo"]),
            X["wide_binary"], "α (0 = Newton, 1 = MOND)", "M3 wide binaries",
            "53-wide-binaries/data/results.json:variants.ecc_hwang", note="in units of MOND's boost; 0.10 ± 0.03 with γ free")
    d10 = _load("40-isw-effect/data/results.json")
    if d10:
        add("D10", "cosmology", "ISW effect, NVSS × Planck (atlas)", d10["amplitudes"]["nvss"] - 1, 0.40, X["cosmology"],
            "A − 1", "D10 ISW", "40-isw-effect/README.md (A = 1.01 ± 0.40)")
    if llr:
        e = llr["EPM2021"]["params"]
        eta, s_eta = e["eta_nordtvedt"]["value"], e["eta_nordtvedt"]["boot"]
        add("N7", "equivalence", "Nordtvedt test in APOLLO residuals (this project)", -4.45e-10 * eta, 4.45e-10 * s_eta,
            X["earth_orbit"], "(m_G/m_I)_Earth − (m_G/m_I)_Moon", "N7 LLR (this work)", "data/results.json:llr.EPM2021",
            note=f"η = ({eta*1e4:.1f} ± {s_eta*1e4:.1f})e-4; residual test, ephemeris not refitted")
        g, s_g = e["gamma_minus_1"]["value"], e["gamma_minus_1"]["boot"]
        add("N7", "weak", "Sun's Shapiro delay in APOLLO ranges (this project)", g / 2, s_g / 2, X["earth_orbit"],
            "(1+γ)/2 − 1", "N7 LLR Shapiro", "data/results.json:llr.EPM2021")
    return R


def table(llr: dict | None = None) -> pd.DataFrame:
    df = pd.DataFrame(literature_rows() + atlas_rows(llr))
    df["pull"] = df["deviation"] / df["sigma"]
    df["precision"] = df["sigma"]
    return df


if __name__ == "__main__":
    t = table()
    with pd.option_context("display.width", 220, "display.max_rows", 200, "display.max_colwidth", 40):
        print(t[["family", "short", "year", "deviation", "sigma", "pull", "x", "atlas"]].to_string(
            float_format=lambda v: f"{v:.3g}"))
