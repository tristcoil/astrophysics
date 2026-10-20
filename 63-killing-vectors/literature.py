"""Published values N6 compares with, each checked against the paper itself.

Verification codes in ``ver``: "T" read in the full text, "A" in the authors'
own abstract only. Checked on 2026-10-19.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Stellar black-hole spins from X-ray spectra (a* = cJ/GM^2)
# ---------------------------------------------------------------------------
# Two methods: CONTINUUM FITTING models the thermal glow of the inner disc,
# whose temperature and size depend on how close the disc's inner edge (the
# ISCO) sits; it needs the black hole's mass, distance and inclination.
# REFLECTION fits the shape of the iron line, smeared by the gas's speed and
# gravitational redshift near the ISCO. Both assume the disc ends at the ISCO
# and is aligned with the spin. Values are model-dependent and some have moved
# a lot between analyses (notes). One headline value per system: the
# continuum-fitting one where it exists.
#
# kind: "value" (lo/hi are the quoted errors), "lower" (a* > value), "upper".
# i4_name: the name under which project I4 (35-x-ray-binaries) lists the
# system's dynamical mass, if it does.

XRB_SPINS = [
    {"name": "GRS 1915+105", "i4_name": "GRS 1915+105", "kind": "lower", "a": 0.98, "cl": "conservative",
     "method": "continuum", "ref": "McClintock et al. 2006, ApJ 652, 518", "arxiv": "astro-ph/0606076", "ver": "T",
     "note": "Mills et al. 2021 (ApJ 914, 6) find 0.4-0.99 depending on the data set"},
    {"name": "Cyg X-1", "i4_name": "Cyg X-1", "kind": "lower", "a": 0.9985, "cl": "3 sigma",
     "method": "continuum", "ref": "Zhao et al. 2021, ApJ 908, 117", "arxiv": "2102.09093", "ver": "T",
     "note": "with a warm Comptonising layer Zdziarski et al. 2024 find much lower spins (0.82 +0.16 -0.74, or <~0.1)"},
    {"name": "LMC X-1", "i4_name": "LMC X-1", "kind": "value", "a": 0.92, "lo": 0.07, "hi": 0.05, "cl": "1 sigma",
     "method": "continuum", "ref": "Gou et al. 2009, ApJ 701, 1076", "arxiv": "0901.0920", "ver": "T", "note": ""},
    {"name": "M33 X-7", "i4_name": "M33 X-7", "kind": "value", "a": 0.84, "lo": 0.05, "hi": 0.05, "cl": "1 sigma",
     "method": "continuum", "ref": "Liu et al. 2008, ApJ 679, L37 (erratum 2010, ApJ 719, L109)",
     "arxiv": "0803.1834", "ver": "T",
     "note": "for the revised mass (11.4 Msun, i = 65 deg) Ramachandran et al. 2022 estimate a* ~ 0.6 (no refit)"},
    {"name": "4U 1543-47", "i4_name": "4U 1543-475", "kind": "value", "a": 0.80, "lo": 0.10, "hi": 0.10,
     "cl": "1 sigma (widened)", "method": "continuum",
     "ref": "Shafee et al. 2006, ApJ 636, L113 (as tabulated by McClintock et al. 2014)", "arxiv": "astro-ph/0508302",
     "ver": "T", "note": "later fits: 0.46 +- 0.12 (Chen et al. 2023), 0.65 +0.14 -0.24 (Yorgancioglu et al. 2023)"},
    {"name": "GRO J1655-40", "i4_name": "GRO J1655-40", "kind": "value", "a": 0.70, "lo": 0.10, "hi": 0.10,
     "cl": "1 sigma (widened)", "method": "continuum",
     "ref": "Shafee et al. 2006, ApJ 636, L113 (as tabulated by McClintock et al. 2014)", "arxiv": "astro-ph/0508302",
     "ver": "T", "note": "reflection gives 0.98 (Miller et al. 2009)"},
    {"name": "XTE J1550-564", "i4_name": "XTE J1550-564", "kind": "value", "a": 0.34, "lo": 0.28, "hi": 0.20,
     "cl": "1 sigma", "method": "continuum", "ref": "Steiner et al. 2011, MNRAS 416, 941", "arxiv": "1010.1013",
     "ver": "T", "note": "90%: -0.11 to 0.71; iron line 0.55 +0.15 -0.22"},
    {"name": "H1743-322", "i4_name": None, "kind": "value", "a": 0.2, "lo": 0.3, "hi": 0.3, "cl": "68%",
     "method": "continuum", "ref": "Steiner et al. 2012, ApJL 745, L7", "arxiv": "1111.2388", "ver": "T",
     "note": "mass assumed (~8 Msun)"},
    {"name": "LMC X-3", "i4_name": "LMC X-3", "kind": "value", "a": 0.21, "lo": 0.12, "hi": 0.12, "cl": "68%",
     "method": "continuum", "ref": "Steiner et al. 2014, ApJL 793, L29", "arxiv": "1402.0148", "ver": "T",
     "note": "90%: +0.18 -0.22"},
    {"name": "A0620-00", "i4_name": "3A 0620-003", "kind": "value", "a": 0.12, "lo": 0.19, "hi": 0.19,
     "cl": "1 sigma", "method": "continuum", "ref": "Gou et al. 2010, ApJL 718, L122", "arxiv": "1002.2211",
     "ver": "T", "note": ""},
    {"name": "Nova Mus 1991", "i4_name": "GRS 1124-684", "kind": "value", "a": 0.63, "lo": 0.19, "hi": 0.16,
     "cl": "1 sigma", "method": "continuum", "ref": "Chen et al. 2016, ApJ 825, 45", "arxiv": "1601.00615",
     "ver": "A", "note": ""},
    {"name": "MAXI J1820+070", "i4_name": "MAXI J1820+070", "kind": "value", "a": 0.14, "lo": 0.09, "hi": 0.09,
     "cl": "1 sigma", "method": "continuum", "ref": "Zhao et al. 2021, ApJ 916, 108", "arxiv": "2012.05544",
     "ver": "A", "note": "~0.3 for the jet's inclination, ~0.94 for the disc's (Fan et al. 2024)"},
    {"name": "MAXI J1305-704", "i4_name": "MAXI J1305-704", "kind": "value", "a": 0.87, "lo": 0.13, "hi": 0.07,
     "cl": "68.3%", "method": "continuum", "ref": "Feng et al. 2023, MNRAS 520, 5803", "arxiv": "2212.04653",
     "ver": "A", "note": ""},
    {"name": "4U 1630-47", "i4_name": None, "kind": "value", "a": 0.92, "lo": 0.04, "hi": 0.04, "cl": "99.7%",
     "method": "continuum", "ref": "Pahari et al. 2018, ApJ 867, 86", "arxiv": "1810.01275", "ver": "T",
     "note": "mass, distance and inclination assumed"},
    {"name": "GX 339-4", "i4_name": "1H 1659-487", "kind": "value", "a": 0.95, "lo": 0.08, "hi": 0.02,
     "cl": "1 sigma", "method": "reflection", "ref": "Parker et al. 2016, ApJL 821, L6", "arxiv": "1603.03777",
     "ver": "T", "note": "continuum fitting gives a* < 0.9 (Kolehmainen & Done 2010)"},
    {"name": "V404 Cyg", "i4_name": "GS 2023+338", "kind": "lower", "a": 0.92, "cl": "99% (statistical)",
     "method": "reflection", "ref": "Walton et al. 2017, ApJ 839, 110", "arxiv": "1609.01293", "ver": "A",
     "note": ""},
    {"name": "XTE J1650-500", "i4_name": "XTE J1650-500", "kind": "value", "a": 0.79, "lo": 0.01, "hi": 0.01,
     "cl": "1 sigma (statistical)", "method": "reflection", "ref": "Miller et al. 2009, ApJ 697, 900",
     "arxiv": "0902.2840", "ver": "T", "note": "Walton et al. 2012: 0.84-0.98 (90%)"},
    {"name": "Swift J1753.5-0127", "i4_name": None, "kind": "value", "a": 0.76, "lo": 0.15, "hi": 0.11,
     "cl": "90%", "method": "reflection", "ref": "Reis et al. 2009, MNRAS 395, 1257", "arxiv": "0902.1745",
     "ver": "T", "note": ""},
    {"name": "MAXI J1535-571", "i4_name": None, "kind": "value", "a": 0.994, "lo": 0.002, "hi": 0.002,
     "cl": "1 sigma", "method": "reflection", "ref": "Miller et al. 2018, ApJL 860, L28", "arxiv": "1806.04115",
     "ver": "A", "note": "other analyses: >0.84 (Xu et al. 2018), 0.67 (Sridhar et al. 2019)"},
]

# ---------------------------------------------------------------------------
# Hawking evaporation (Page 1976, Phys. Rev. D 13, 198; Page 2005, NJP 7, 203)
# ---------------------------------------------------------------------------
# dM/dt = -alpha hbar c^4 / (G^2 M^2), lifetime tau = G^2 M^3 / (3 alpha hbar c^4).
PAGE = {
    "alpha_photons_gravitons_neutrinos": 2.011e-4,   # Page 1976 eq. 29: 4 one-helicity neutrino kinds
    "alpha_photons_gravitons": 3.747e-5,             # Page 2005 eq. 19 (n_1 = n_2 = 1)
    "lifetime_1msun_yr_with_neutrinos": 2.16e66,     # Page 1976 eq. 27
    "photon_blackbody_coefficient": 5120.0,          # textbook tau = 5120 pi G^2 M^3/(hbar c^4)
    "ref": "Page 1976, Phys. Rev. D 13, 198; Page 2005, New J. Phys. 7, 203 (hep-th/0409024)",
}
#: Checks: T_H and S for one solar mass (CODATA 2018, GM_sun = 1.32712440018e20).
HAWKING_1MSUN = {"T_K": 6.1703e-8, "S_over_kB": 1.0489e77}

# ---------------------------------------------------------------------------
# GW250114 (LVK 2025, PRL 135, 111403; arXiv:2509.08054), 90% credible intervals
# ---------------------------------------------------------------------------
GW250114_PAPER = {
    "m1": (33.6, 1.2, 0.8), "m2": (32.2, 0.8, 1.3), "M_final": (62.7, 1.0, 1.1), "chi_final": (0.68, 0.01, 0.01),
    "area_law_sigma": 4.4, "waveform": "NRSur7dq4", "ver": "T",
    "area_note": ("the paper gives the horizon areas only in figures; the area theorem is confirmed at 4.4 sigma "
                  "using the inspiral before t = -40 t_M and the ringdown after the peak"),
}

# ---------------------------------------------------------------------------
# S2's Schwarzschild precession measured by GRAVITY (f_SP = 1 is GR)
# ---------------------------------------------------------------------------
GRAVITY_FSP = [
    {"label": "GRAVITY 2020 (S2)", "f": 1.10, "err": 0.19, "ref": "A&A 636, L5", "arxiv": "2004.07187",
     "prediction_arcmin": 12.10, "mass_1e6": 4.261, "r0_pc": 8246.7},
    {"label": "GRAVITY 2022 (S2 only)", "f": 0.85, "err": 0.16, "ref": "A&A 657, L12", "arxiv": "2112.07478"},
    {"label": "GRAVITY 2024 (S2 only)", "f": 0.918, "err": 0.128, "ref": "A&A 692, A242", "arxiv": "2409.12261"},
    {"label": "GRAVITY 2024 (four stars)", "f": 1.135, "err": 0.110, "ref": "A&A 692, A242", "arxiv": "2409.12261"},
]

# ---------------------------------------------------------------------------
# Shadows
# ---------------------------------------------------------------------------
SHADOW = {
    "kerr_diameter_range_M": (9.6, 10.4),   # EHT M87 Paper VI (ApJL 875, L6), over spin and inclination
    "m87_inclination_deg": 17.0,            # jet angle: Walker et al. 2018 (ApJ 855, 128); Mertens et al. 2016: 17.2 +- 3.3
    "sgra_inclination_deg_max": 30.0,       # EHT Sgr A* Paper V (ApJL 930, L16): favoured i <= 30 deg
    "m87_eht_theta_g_uas": (3.8, 0.4),      # EHT M87 Paper VI
    "sgra_eht_ring_uas": (51.8, 2.3),       # EHT Sgr A* Paper IV
}
