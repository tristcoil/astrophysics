"""Published measurements typed in from their primary sources.

Small tables that exist only in a PhD thesis or a paper's TeX are typed here
with their exact source; larger machine-readable tables (Fermi-LAT on VizieR,
CRDB) are downloaded by ``fetch_data.py``. Each entry was checked against the
primary document (thesis page, table number) -- see the comments.

Units: gamma-ray intensities are *bin-averaged* differential intensities
I = dN/dE in photons cm^-2 s^-1 sr^-1 MeV^-1 unless stated otherwise.
"""

from __future__ import annotations

import math

import numpy as np

# ---------------------------------------------------------------------------
# The MeV gamma-ray background (where annihilation photons from z ~ 20-1100
# would land today)
# ---------------------------------------------------------------------------

#: COMPTEL, G. Weidenspointner 1999, PhD thesis, TU Munich, "The Origin of the
#: Cosmic Gamma-Ray Background in the COMPTEL Energy Range", Table 6.6
#: (thesis p. 249; mediatum.ub.tum.de/doc/602832). The analysis summarised in
#: Weidenspointner et al. 2000 (AIP Conf. Proc. 510, 467). All pointings at
#: |b| > 30 deg; Galactic diffuse emission and point sources NOT subtracted
#: (so these are upper bounds on any truly extragalactic signal).
#: Columns: E_lo, E_hi (MeV), I, sigma_stat, sigma_syst; the thesis's total
#: error is the *linear* sum. The 1.8-2.7 MeV systematic is printed 0.351e-4,
#: a typo for 0.351e-3 (its Table 6.5 has 3.51e-4, and sigma_tot = 0.559e-3).
COMPTEL_WEIDENSPOINTNER1999 = np.array([
    [0.8, 1.2, 6.183e-3, 0.719e-3, 2.639e-3],
    [1.2, 1.8, 1.384e-3, 0.557e-3, 0.966e-3],
    [1.8, 2.7, 1.014e-3, 0.208e-3, 0.351e-3],
    [2.7, 4.3, 1.909e-4, 0.925e-4, 0.723e-4],
    [4.3, 6.0, 9.710e-5, 1.820e-5, 4.372e-5],
    [6.0, 9.0, 4.217e-5, 0.778e-5, 1.678e-5],
    [9.0, 12.0, 2.076e-5, 0.195e-5, 0.458e-5],
    [12.0, 17.0, 1.288e-5, 0.114e-5, 0.275e-5],
    [17.0, 30.0, 5.662e-6, 0.721e-6, 1.437e-6],
])

#: COMPTEL, S. C. Kappadath 1998, PhD thesis, Univ. of New Hampshire,
#: "Measurement of the cosmic diffuse gamma-ray spectrum from 800 keV to 30 MeV",
#: Table VII.2 (thesis p. 189): an independent analysis (Virgo + South
#: Galactic Pole pointings), Galactic emission not subtracted. Columns: E_lo,
#: E_hi, I, stat, syst. The 12-17 MeV value is printed "11.13e-6" (= 1.113e-5,
#: confirmed by the thesis's 9-30 MeV integral of 1.5e-4).
COMPTEL_KAPPADATH1998 = np.array([
    [0.8, 1.2, 4.65e-3, 1.1e-3, 3.9e-3],
    [1.2, 1.8, 2.0e-3, 0.46e-3, 1.05e-3],
    [1.8, 2.7, 1.21e-3, 0.29e-3, 0.34e-3],
    [2.7, 4.2, 2.45e-4, 1.23e-4, 0.86e-4],
    [4.2, 6.0, 6.88e-5, 2.43e-5, 3.47e-5],
    [6.0, 9.0, 2.88e-5, 1.32e-5, 1.34e-5],
    [9.0, 12.0, 1.88e-5, 0.383e-5, 0.293e-5],
    [12.0, 17.0, 1.113e-5, 0.211e-5, 0.227e-5],
    [17.0, 30.0, 2.77e-6, 1.03e-6, 1.20e-6],
])

#: EGRET extragalactic background as re-derived by Strong, Moskalenko &
#: Reimer 2004 (ApJ 613, 956; arXiv:astro-ph/0405441, Table 3): bin-averaged
#: intensity with statistical and EGRET systematic errors in quadrature.
#: The last band (20-50 GeV) is a 1 sigma upper limit (error set to NaN).
EGRET_STRONG2004 = np.array([
    [30.0, 50.0, 8.40e-7, 2.54e-7],
    [50.0, 70.0, 5.30e-7, 0.80e-7],
    [70.0, 100.0, 2.22e-7, 0.22e-7],
    [100.0, 150.0, 8.96e-8, 0.91e-8],
    [150.0, 300.0, 2.61e-8, 0.26e-8],
    [300.0, 500.0, 6.00e-9, 0.63e-9],
    [500.0, 1000.0, 1.52e-9, 0.17e-9],
    [1000.0, 2000.0, 3.20e-10, 0.44e-10],
    [2000.0, 4000.0, 1.20e-10, 0.16e-10],
    [4000.0, 10000.0, 1.95e-11, 0.30e-11],
    [10000.0, 20000.0, 3.40e-12, 1.39e-12],
    [20000.0, 50000.0, 1.11e-12, np.nan],
])

#: Swift-BAT cosmic X-ray background, Ajello et al. 2008 (ApJ 689, 666;
#: arXiv:0808.3377, Table 1): E (keV), E^2 dN/dE (keV^2 cm^-2 s^-1 sr^-1 keV^-1),
#: error -- for the low-energy end of the plot only.
SWIFT_BAT_AJELLO2008 = np.array([
    [15.0, 43.92, 12.71], [17.0, 46.82, 10.35], [19.0, 48.58, 8.65], [21.0, 50.39, 7.42],
    [24.0, 50.96, 5.28], [28.1, 47.12, 2.90], [32.2, 45.26, 2.72], [36.2, 45.07, 2.52],
    [41.4, 44.00, 2.34], [47.8, 40.39, 2.07], [54.3, 37.33, 1.89], [62.0, 33.25, 1.63],
    [70.9, 32.99, 1.79], [81.2, 28.08, 2.17], [92.9, 26.25, 3.50], [106.1, 24.86, 1.61],
    [120.8, 23.18, 1.93], [137.2, 21.15, 2.32], [156.8, 18.95, 2.62], [181.1, 18.02, 2.96],
])

#: HEAO-1 A4 (Kinzer et al. 1997, ApJ 475, 361): no table, only their fit
#: dN/dE = 2.62e-3 (E / 100 keV)^-2.75 photons cm^-2 s^-1 keV^-1 sr^-1,
#: significant for ~80-350 keV.
HEAO1_A4_KINZER1997 = {"norm_per_kev": 2.62e-3, "e0_kev": 100.0, "index": 2.75,
                       "range_kev": (80.0, 350.0)}


def comptel(which: str = "weidenspointner"):
    """COMPTEL points: dict of arrays (E_lo, E_hi, E_mid [geometric], I, err).

    ``err`` for Weidenspointner is the thesis's linear stat + syst sum; for
    Kappadath the quadrature sum (as that thesis states for its total)."""
    if which == "weidenspointner":
        t = COMPTEL_WEIDENSPOINTNER1999
        err = t[:, 3] + t[:, 4]
    else:
        t = COMPTEL_KAPPADATH1998
        err = np.hypot(t[:, 3], t[:, 4])
    return {"E_lo": t[:, 0], "E_hi": t[:, 1], "E_mid": np.sqrt(t[:, 0] * t[:, 1]),
            "I": t[:, 2], "err": err, "stat": t[:, 3], "syst": t[:, 4]}


def egret():
    t = EGRET_STRONG2004
    return {"E_lo": t[:, 0], "E_hi": t[:, 1], "E_mid": np.sqrt(t[:, 0] * t[:, 1]),
            "I": t[:, 2], "err": t[:, 3], "is_ul": ~np.isfinite(t[:, 3])}


def swift_bat():
    t = SWIFT_BAT_AJELLO2008
    e_mev = t[:, 0] * 1e-3
    # keV^2 cm^-2 s^-1 sr^-1 keV^-1 -> E^2 I in MeV cm^-2 s^-1 sr^-1
    return {"E_mid": e_mev, "E2I": t[:, 1] * 1e-3, "E2I_err": t[:, 2] * 1e-3}


# ---------------------------------------------------------------------------
# Primordial deuterium: D/H in near-pristine gas clouds seen against quasars
# ---------------------------------------------------------------------------

#: The 12 systems averaged by the PDG's BBN review (2025 update = 2026
#: edition, Takahashi et al., IJMPA 41, 2630011; 25.08 +- 0.29 x 10^-6 after
#: scaling the error by S = 1.08). Columns: name, z_abs, metallicity [X/H]
#: (oxygen unless noted), D/H x 10^5, error minus, error plus, reference.
#: Values published as log10(D/H) +- s are converted exactly (asymmetric).
#: The first seven are Cooke, Pettini & Steidel's (2018, ApJ 855, 102) Table 3
#: "precision sample" (weighted mean 2.527 +- 0.030).
def _from_log(logv: float, s: float):
    v = 10 ** (logv + 5)
    return v, v - 10 ** (logv - s + 5), 10 ** (logv + s + 5) - v


DH_SYSTEMS = [
    # name, z, [X/H], (D/H, -err, +err) x 1e5, reference, in Cooke+18 precision sample
    ("HS 0105+1619", 2.53651, -1.771, _from_log(-4.589, 0.026), "Cooke+14 (O'Meara+01 data)", True),
    ("Q0913+072", 2.61829, -2.416, _from_log(-4.597, 0.018), "Cooke+14 (Pettini+08 data)", True),
    ("Q1243+307", 2.52564, -2.769, _from_log(-4.622, 0.015), "Cooke+18 (Kirkman+03 data)", True),
    ("J1358+0349", 2.85305, -2.804, _from_log(-4.582, 0.012), "Cooke+16, ApJ 830, 148", True),
    ("J1358+6522", 3.06726, -2.335, _from_log(-4.588, 0.012), "Cooke+14, ApJ 781, 31", True),
    ("J1419+0829", 3.04973, -1.922, _from_log(-4.601, 0.009), "Cooke+14 (Pettini & Cooke 12)", True),
    ("J1558-0031", 2.70242, -1.650, _from_log(-4.619, 0.026), "Cooke+14 (O'Meara+06 data)", True),
    ("PKS 1937-1009 (z=3.256)", 3.256, -1.87, _from_log(-4.610, 0.050),
     "Riemer-Sorensen+15, MNRAS 447, 2925", False),
    ("J1444+2919", 2.437, -2.042, (1.97, 0.28, 0.33),
     "Balashev+16, MNRAS 458, 2188 (stat+syst)", False),
    ("PKS 1937-1009 (z=3.572)", 3.572, -2.52, (2.608, 0.102, 0.102),
     "Guarneri+24, MNRAS 529, 839 (ESPRESSO+UVES)", False),
    ("Q1009+2956", 2.504, -2.5, _from_log(-4.606, 0.066),
     "Zavarygin+18, MNRAS 477, 5536 ([Si/H])", False),
    ("J1332+0052", 3.42, -1.71, (2.388, 0.078, 0.078), "Kislitsyn+24, MNRAS 528, 4068", False),
]

#: Published averages (for comparison).
DH_PUBLISHED = {
    "cooke2018": (2.527, 0.030, "Cooke, Pettini & Steidel 2018, eqs. 9-10 (7 systems)"),
    "pdg2024": (2.547, 0.029, "PDG 2024 BBN review, eq. 24.2 (11 systems, S = 1.137)"),
    "pdg2026": (2.508, 0.029, "PDG 2025 update / 2026 edition (12 systems, S = 1.08)"),
    "pettini_cooke2026": (2.510, 0.028, "Pettini & Cooke 2026, Ap&SS 371, 20 (8 systems)"),
}


def dh_weighted_mean(systems=None, *, symmetric: str = "average"):
    """Inverse-variance weighted mean of D/H (x 1e5), PDG style: each
    asymmetric error replaced by the mean of its two sides; the error is
    scaled by S = sqrt(chi^2 / (N - 1)) when S > 1."""
    rows = DH_SYSTEMS if systems is None else systems
    v = np.array([r[3][0] for r in rows])
    e = np.array([0.5 * (r[3][1] + r[3][2]) for r in rows])
    w = 1 / e**2
    mean = np.sum(w * v) / np.sum(w)
    err = np.sum(w) ** -0.5
    chi2 = float(np.sum(w * (v - mean) ** 2))
    s = max(1.0, math.sqrt(chi2 / (len(v) - 1)))
    return {"mean": float(mean), "err": float(err), "err_scaled": float(err * s), "S": s,
            "chi2": chi2, "n": len(v)}


def dh_log_mean(systems=None):
    """Cooke et al.'s way: weighted mean of log10(D/H) (symmetric log errors)."""
    rows = DH_SYSTEMS if systems is None else systems
    v = np.array([math.log10(r[3][0] * 1e-5) for r in rows])
    e = np.array([math.log10((r[3][0] + r[3][2]) / r[3][0]) for r in rows])
    w = 1 / e**2
    m = np.sum(w * v) / np.sum(w)
    s = np.sum(w) ** -0.5
    return {"log_mean": float(m), "log_err": float(s), "mean": 10 ** (m + 5),
            "err": 10 ** (m + 5) * math.log(10) * s}


# ---------------------------------------------------------------------------
# The baryon density from the CMB (Omega_b h^2), and from BBN in the literature
# ---------------------------------------------------------------------------

CMB_OMBH2 = [
    # label, value, error, source
    ("Planck 2018 TT+lowE", 0.02212, 0.00022, "Planck 2018 VI, A&A 641, A6, Table 2"),
    ("Planck 2018 TT,TE,EE+lowE+lensing", 0.02237, 0.00015, "Planck 2018 VI, Table 2"),
    ("ACT DR6 + Planck (P-ACT)", 0.02250, 0.00011, "Louis et al. 2025, JCAP 11, 062, Table 5"),
    ("SPT-3G D1 alone", 0.02221, 0.00020, "Camphuis et al. 2026, PRD 113, 083504, Table I"),
    ("SPT + ACT + Planck (CMB-SPA)", 0.022398, 0.000095, "Camphuis et al. 2026, Table I"),
]

BBN_OMBH2_LITERATURE = [
    ("Cooke+18, theory d(p,g) rate", 0.02166, math.hypot(0.00015, 0.00011), "ApJ 855, 102, eq. 12"),
    ("Pitrou+21 (PRIMAT, LUNA)", 0.02195, 0.00022, "MNRAS 502, 2474, eq. 9"),
    ("Pisanti+21 (PArthENoPE, LUNA)", 0.02233, 0.00036, "JCAP 04 (2021) 020, Table 3"),
    ("Yeh+26 (BBN + D)", 0.02227, 0.00043, "ApJ 1008, 241, Table 1"),
    ("Schoeneberg 2024 (recommended)", 0.02218, 0.00055, "JCAP 06 (2024) 006"),
]


# ---------------------------------------------------------------------------
# Sakharov's conditions in the Standard Model
# ---------------------------------------------------------------------------

#: CKM matrix, PDG review "CKM quark-mixing matrix" (Ceccucci, Ligeti, Sakai),
#: global fit (CKMfitter method). 2026 edition (Takahashi et al., IJMPA 41,
#: 2630011; review revised March 2026) and 2024 edition (PRD 110, 030001).
#: Values: (central, minus, plus). J = Jarlskog invariant.
CKM = {
    "pdg2026": {"lambda": (0.22517, 0.00068, 0.00068), "A": (0.826, 0.015, 0.017),
                "rhobar": (0.1576, 0.0091, 0.0092), "etabar": (0.3556, 0.0069, 0.0071),
                "J": (3.16e-5, 0.11e-5, 0.13e-5), "s12": 0.22517, "s13": 0.003763,
                "s23": 0.04189, "delta": 1.154},
    "pdg2024": {"lambda": (0.22501, 0.00068, 0.00068), "A": (0.826, 0.015, 0.016),
                "rhobar": (0.1591, 0.0094, 0.0094), "etabar": (0.3523, 0.0071, 0.0073),
                "J": (3.12e-5, 0.12e-5, 0.13e-5), "s12": 0.22501, "s13": 0.003732,
                "s23": 0.04183, "delta": 1.147},
}

#: Quark masses (GeV). PDG 2026 summary tables (u, d, s: MS-bar at 2 GeV;
#: c, b: m(m); t: direct measurements) -- errors there are 90% CL.
QUARK_MASSES_PDG2026 = {"u": 2.16e-3, "d": 4.70e-3, "s": 92.9e-3, "c": 1.2729, "b": 4.186,
                        "t": 172.60}
#: Running MS-bar masses at mu = M_Z in the full Standard Model, Huang & Zhou
#: 2021 (PRD 103, 016010, Table II) -- the masses quarks actually have at the
#: electroweak scale.
QUARK_MASSES_MZ = {"u": 1.23e-3, "d": 2.67e-3, "s": 53.16e-3, "c": 0.620, "b": 2.839,
                   "t": 168.26}

#: Electroweak bosons and the Higgs (PDG 2024/2026, GeV).
EW_MASSES = {"W": 80.369, "Z": 91.1876, "H": 125.20, "t": 172.57}

#: Lattice results on the electroweak transition with the real Higgs mass.
LATTICE_EW = {
    # smooth crossover temperature (D'Onofrio & Rummukainen 2016, PRD 93, 025003)
    "T_c": (159.5, 1.5),
    # sphaleron rate, D'Onofrio, Rummukainen & Tranberg 2014 (PRL 113, 141602):
    # symmetric phase Gamma/T^4 = (8.0 +- 1.3)e-7 ~ (18 +- 3) alpha_W^5 (eq. 8);
    # Higgs phase ln(Gamma/T^4) = (0.83 +- 0.01) T/GeV - (147.7 +- 1.9) for
    # 130 GeV < T < T_c (eq. 7, natural log); freeze-out where
    # Gamma(T*)/T*^3 = alpha H(T*), alpha = 0.1015, g* = 106.75 (eq. 9).
    "sphaleron_symmetric": (8.0e-7, 1.3e-7),
    "sphaleron_broken": {"a": 0.83, "a_err": 0.01, "b": 147.7, "b_err": 1.9,
                         "valid": (130.0, 159.5)},
    "freezeout_alpha": 0.1015,
    "freezeout_gstar": 106.75,
    "T_star": (131.7, 2.3),
    # end of the first-order line (Csikor, Fodor & Heitger 1999, PRL 82, 21:
    # 4D, SM; Laine & Rummukainen 1999, NPPS 73, 180: 3D, SM 72.3 +- 0.7)
    "m_h_endpoint": (72.4, 1.7),
    # order parameter v(T_c)/T_c in the full SM (m_t = 175 GeV) from
    # Kajantie et al. 1996, NPB 466, 189, Table 10
    "v_over_T": [(51.2, 0.689), (68.0, 0.575)],
}

#: Sphalerons must stop before they erase a baryon excess made at the
#: transition: v(T_c)/T_c >~ 1 (conventional; Morrissey & Ramsey-Musolf 2012,
#: NJP 14, 125003); a gauge-invariant benchmark 0.73 with a plausible range
#: 0.4-1.4 (Patel & Ramsey-Musolf 2011, JHEP 07, 029).
WASHOUT = {"conventional": 1.0, "benchmark": 0.73, "range": (0.4, 1.4)}

#: What is needed: the observed asymmetry Y_B = n_B/s ~ 8.7e-11, and the
#: literature's statement of the Standard Model's CP measure, ~1e-20
#: (e.g. Rubakov & Shaposhnikov 1996; Canetti, Drewes & Shaposhnikov 2012).
SM_CP_LITERATURE = 1e-20

#: Antihydrogen falls (ALPHA-g, Anderson et al. 2023, Nature 621, 716):
#: a = (0.75 +- 0.13 (stat + syst) +- 0.16 (simulation)) g; no gravity is
#: disfavoured at p = 2.9e-4, repulsive gravity at p < 1e-15.
ALPHA_G = {"a_over_g": 0.75, "stat_syst": 0.13, "simulation": 0.16, "p_no_gravity": 2.9e-4}


# ---------------------------------------------------------------------------
# Proton-antiproton annihilation at rest (the source spectrum of every
# annihilation signal)
# ---------------------------------------------------------------------------

#: Multi-pion channel fractions (x 1e-3) in liquid hydrogen, Crystal Barrel
#: column of Klempt, Batty & Richard 2005 (Phys. Rep. 413, 197), Table 7;
#: the six-prong channels from the bubble-chamber columns. Keys are
#: (number of charged pions, number of neutral pions). They sum to ~88%;
#: the rest is kaons and eta/omega final states.
PBARP_CHANNELS = {
    (2, 0): 3.14, (0, 2): 0.65,
    (2, 1): 67.0, (0, 3): 7.0,
    (4, 0): 56.0, (2, 2): 122.0, (0, 4): 3.1,
    (4, 1): 210.0, (2, 3): 133.0, (0, 5): 9.2,
    (6, 0): 20.0, (4, 2): 177.0, (2, 4): 36.0, (0, 6): 0.12,
    (6, 1): 17.0, (4, 3): 6.0, (2, 5): 13.0, (0, 7): 1.3,
}

#: Mean multiplicities at rest, KBR 2005 eq. 5.3, and the photon number
#: measured by Backenstoss et al. 1983 (NPB 228, 424) as quoted there.
PBARP_MULTIPLICITY = {"n_pi": (4.98, 0.13), "n_pic": (3.05, 0.04), "n_pi0": (1.93, 0.12),
                      "n_gamma": (3.93, 0.24)}

#: Low-energy annihilation cross-section times velocity: OBELIX measures
#: beta * sigma_ann = 40-48 mb for p_lab = 44-175 MeV/c (Bertin et al. 1996,
#: tabulated by Carbonell, Protasov & Zenoni 1997, PLB 397, 345), and the
#: partial-wave analysis of Zhou & Timmermans 2012 (PRC 86, 044003) gives
#: sigma_ann = 198.1 / 108.5 / 81.9 / 68.6 mb at 200 / 400 / 600 / 800 MeV/c,
#: i.e. beta * sigma = 41-46 mb: sigma_ann * v is nearly constant,
#: ~1.3e-15 cm^3/s, from 1 MeV to ~150 MeV of kinetic energy.
PBARP_ZT2012 = [(200.0, 198.1), (400.0, 108.5), (600.0, 81.9), (800.0, 68.6)]

#: The measured inclusive photon spectrum of p pbar annihilation at rest
#: (Backenstoss et al. 1983, NPB 228, 424, their analytic shape function), as
#: plotted by von Ballmoos 2014 (arXiv:1401.7258, Fig. 1) and DIGITISED from
#: that bitmap (a few per cent accuracy). Used only to validate our
#: phase-space spectrum. Columns: E (MeV), dN/dE (photons / MeV / annihilation).
BACKENSTOSS_DIGITISED = np.array([
    [5, 2.08e-3], [10, 4.93e-3], [20, 9.22e-3], [30, 12.0e-3], [50, 14.9e-3], [70, 15.6e-3],
    [100, 14.9e-3], [150, 12.2e-3], [200, 9.44e-3], [300, 5.30e-3], [400, 2.72e-3],
    [500, 1.24e-3], [600, 0.489e-3], [700, 0.158e-3], [800, 0.040e-3], [850, 0.016e-3],
])
