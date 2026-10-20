"""Published values used by the analysis, each with its source.

Everything here was read off the papers themselves (arXiv versions); the
README lists the exact tables. Measurement distributions use the notation of
Farr et al. (2011): ("N", mu, sigma) a Gaussian, ("SN", median, plus, minus) a
split Gaussian for asymmetric errors, ("U", a, b) uniform, and ("I", a, b) an
isotropic inclination, i.e. cos i uniform between angles a and b (degrees).
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Black holes outside BlackCAT (which lists only transients)
# ---------------------------------------------------------------------------

#: Persistent, wind- or overflow-fed black holes with massive companions.
#: Masses come from full light-curve + velocity models; we adopt them as
#: published (68% intervals), and recompute Cyg X-1 from its orbit as a check.
PERSISTENT_BHS = [
    {
        "name": "Cyg X-1", "host": "Milky Way", "type": "HMXB",
        "m_bh": ("SN", 21.2, 2.2, 2.3), "m_companion": ("SN", 40.6, 7.7, 7.1),
        "porb_d": 5.599829, "k_kms": ("SN", 75.21, 0.42, 0.41), "e": 0.0189,
        "incl_deg": ("SN", 27.51, 0.77, 0.57),
        "ref": "Miller-Jones et al. 2021, Science 371, 1046 (Table 1); "
               "P from Brocksopp et al. 1999, A&A 343, 861",
        "bibcode": "2021Sci...371.1046M",
    },
    {
        "name": "LMC X-1", "host": "LMC", "type": "HMXB",
        "m_bh": ("N", 10.91, 1.41), "m_companion": ("N", 31.79, 3.48),
        "porb_d": 3.90917, "k_kms": ("N", 71.61, 0.67), "incl_deg": ("N", 36.38, 1.92),
        "ref": "Orosz et al. 2009, ApJ 697, 573", "bibcode": "2009ApJ...697..573O",
    },
    {
        "name": "LMC X-3", "host": "LMC", "type": "persistent LMXB/IMXB",
        "m_bh": ("N", 6.98, 0.56), "m_companion": ("N", 3.63, 0.57),
        "porb_d": 1.7048089, "k_kms": ("N", 241.1, 6.2), "incl_deg": ("N", 69.24, 0.72),
        "ref": "Orosz et al. 2014, ApJ 794, 154", "bibcode": "2014ApJ...794..154O",
    },
    {
        "name": "M33 X-7", "host": "M33", "type": "HMXB",
        "m_bh": ("SN", 11.4, 3.3, 1.7), "m_companion": ("SN", 38.0, 22.0, 10.0),
        "porb_d": 3.453014, "k_kms": ("N", 108.9, 5.7), "incl_deg": None,
        "ref": "Ramachandran et al. 2022, A&A 667 (doi:10.1051/0004-6361/202243683); "
               "was 15.65 +- 1.45 in Orosz et al. 2007, Nature 449, 872",
        "bibcode": "doi:10.1051/0004-6361/202243683",
    },
]

#: Wolf-Rayet systems once counted among the heaviest stellar black holes.
#: Their masses rest on He II emission-line velocities that probably trace a
#: shadowed part of the wind rather than the star: IC 10 X-1's mass function
#: is then unknown (Laycock et al. 2015), and NGC 300 X-1's black hole may be
#: much lighter than published (Binder et al. 2021). Not used here.
REJECTED_BHS = {
    "IC 10 X-1": "Laycock, Maccarone & Christodoulou 2015, MNRAS 452, L31 (arXiv:1506.03882)",
    "NGC 300 X-1": "Binder et al. 2021, ApJ 910, 74 (arXiv:2102.07065)",
}

#: Dormant black holes found from Gaia astrometry (no X-rays).
GAIA_BHS = [
    {
        "name": "Gaia BH1", "m_bh": ("N", 9.62, 0.18), "m_star": ("N", 0.93, 0.05),
        "porb_d": 185.59, "distance_pc": 480, "feh": -0.2,
        "ref": "El-Badry et al. 2023a, MNRAS 518, 1057 (Table 1)",
        "bibcode": "2023MNRAS.518.1057E",
    },
    {
        "name": "Gaia BH2", "m_bh": ("N", 8.94, 0.34), "m_star": ("N", 1.07, 0.19),
        "porb_d": 1276.7, "distance_pc": 1160, "feh": -0.22,
        "ref": "El-Badry et al. 2023b, MNRAS 521, 4323 (Table 2)",
        "bibcode": "2023MNRAS.521.4323E",
    },
    {
        "name": "Gaia BH3", "m_bh": ("N", 32.70, 0.82), "m_star": ("N", 0.76, 0.05),
        "porb_d": 4253.1, "distance_pc": 590.6, "feh": -2.56,
        "ref": "Gaia Collaboration, Panuzzo et al. 2024, A&A 686, L2",
        "bibcode": "2024A&A...686L...2G",
    },
]

#: What El-Badry et al. derived from the Gaia DR3 solutions *alone*, the
#: check for our own Thiele-Innes conversion (gaia_orbits.py).
GAIA_DR3_ONLY_CHECKS = {
    "Gaia BH1": {
        "a0_mas": (3.00, 0.22), "incl_deg": (121.2, 2.8), "period_d": (185.77, 0.31),
        "e": (0.49, 0.07), "Omega_deg": (89.6, 3.7), "omega_deg": (-10.3, 11.6),
        "f_ast": (11.2, 1.9), "m2": (12.8, 2.0),
        "joint_m2": (9.62, 0.18), "rv_only_fm": (4.08, 0.08), "rv_only_k": (66.7, 0.6),
        "m_star": (0.93, 0.05),
        "ref": "El-Badry et al. 2023a, Table 1 ('astrometry only') and Sect. 4.3.1",
    },
    "Gaia BH2": {
        "period_d": (1352.0, 46.0), "f_ast": (6.72, 0.5), "m2_for_mstar_1": (8.4, 0.5),
        "joint_m2": (8.94, 0.34), "m_star": (1.07, 0.19),
        "ref": "El-Badry et al. 2023b, Sect. 1 (DR3 solution) and Table 2",
    },
}

# ---------------------------------------------------------------------------
# Neutron stars weighed in X-ray binaries (the same technique as the black
# holes): Ozel & Freire 2016, ARA&A 54, 401, Table 3. The burst-based
# masses of that table's LMXBs are not dynamical and are left out.
# ---------------------------------------------------------------------------

NS_XRB = [
    ("LMC X-4", 1.57, 0.11, 0.11), ("Cen X-3", 1.57, 0.16, 0.16),
    ("4U 1538-522", 1.02, 0.17, 0.17), ("SMC X-1", 1.21, 0.12, 0.12),
    ("SAX J1802.7-2017", 1.57, 0.25, 0.25), ("XTE J1855-026", 1.41, 0.24, 0.24),
    ("Vela X-1", 2.12, 0.16, 0.16), ("EXO 1722-363", 1.91, 0.45, 0.45),
    ("OAO 1657-415", 1.74, 0.30, 0.30), ("Her X-1", 1.07, 0.36, 0.36),
    ("Cyg X-2", 1.90, 0.22, 0.35),
]
NS_XRB_REF = "Ozel & Freire 2016, ARA&A 54, 401, Table 3 (Falanga et al. 2015; Orosz & Kuulkers 1999)"

# ---------------------------------------------------------------------------
# Input tables of the two classic mass-gap analyses, for replication
# ---------------------------------------------------------------------------

#: Farr et al. 2011, ApJ 741, 103, Table 1: 15 Roche-lobe-overflow systems
#: (low-mass donors) and 5 wind-fed, high-mass systems.
FARR2011_LOW = {
    "GRS 1915+105": (("N", 9.5, 3.0), ("N", 0.0857, 0.0284), ("N", 70, 2)),
    "XTE J1118+480": (("N", 6.44, 0.08), ("N", 0.0264, 0.004), ("N", 68, 2)),
    "XTE J1650-500": (("N", 2.73, 0.56), ("U", 0, 0.5), ("I", 50, 80)),
    "GRS 1009-45": (("N", 3.17, 0.12), ("N", 0.137, 0.015), ("I", 37, 80)),
    "A0620-00": (("N", 2.76, 0.036), ("N", 0.06, 0.004), ("N", 50.98, 0.87)),
    "GRO J0422+32": (("N", 1.13, 0.09), ("U", 0.076, 0.31), ("N", 45, 2)),
    "Nova Mus 1991": (("N", 3.01, 0.15), ("N", 0.128, 0.04), ("N", 54, 1.5)),
    "GRO J1655-40": (("N", 2.73, 0.09), ("N", 0.3663, 0.04025), ("N", 70.2, 1.9)),
    "4U 1543-47": (("N", 0.25, 0.01), ("U", 0.25, 0.31), ("N", 20.7, 1.5)),
    "XTE J1550-564": (("N", 7.73, 0.4), ("U", 0, 0.04), ("N", 74.7, 3.8)),
    "V4641 Sgr": (("N", 3.13, 0.13), ("U", 0.42, 0.45), ("N", 75, 2)),
    "GS 2023+338": (("N", 6.08, 0.06), ("U", 0.056, 0.063), ("I", 66, 70)),
    "GS 1354-64": (("N", 5.73, 0.29), ("N", 0.12, 0.04), ("I", 50, 80)),
    "Nova Oph 77": (("N", 4.86, 0.13), ("U", 0, 0.053), ("I", 60, 80)),
    "GS 2000+25": (("N", 5.01, 0.12), ("U", 0.035, 0.053), ("I", 43, 74)),
}
FARR2011_HIGH = {
    "Cyg X-1": (("N", 0.251, 0.007), ("N", 2.778, 0.386), ("I", 23, 38)),
    "M33 X-7": (("N", 0.46, 0.08), ("N", 4.47, 0.61), ("N", 74.6, 1)),
    "NGC 300 X-1": (("N", 2.6, 0.3), ("U", 1.05, 1.65), ("I", 60, 75)),
    "LMC X-1": (("N", 0.148, 0.004), ("N", 2.91, 0.49), ("N", 36.38, 2.02)),
    "IC 10 X-1": (("N", 7.64, 1.26), ("U", 0.7, 1.7), ("I", 75, 90)),
}

#: Kreidberg et al. 2012, ApJ 757, 36, Table 2: the same systems after a
#: closer look at how disc light biases inclinations. i_min assumes a donor
#: of normal mass for its spectral type; i_max is the no-eclipse limit.
#: (Their SAX J1819.3-2525 inclination, I(80, 90), is deliberately above its
#: no-eclipse limit: see their Sect. 5.12. We use every row as printed.)
KREIDBERG2012 = {
    "4U 1543-47": (("N", 0.25, 0.01), ("U", 0.25, 0.31), ("I", 24.0, 36.0), 19.8, 73.4),
    "A0620-00": (("N", 3.1, 0.04), ("N", 0.060, 0.004), ("N", 51.0, 0.9), 36.5, 79.8),
    "GRO J0422+32": (("N", 1.19, 0.02), ("N", 0.116, 0.08), ("N", 63.7, 5.2), 0.0, 90.0),
    "GRO J1655-40": (("N", 2.73, 0.09), ("N", 0.38, 0.05), ("N", 69.0, 3.0), 0.0, 90.0),
    "GRS 1009-45": (("N", 3.17, 0.12), ("N", 0.137, 0.015), ("N", 62.0, 5.1), 42.3, 76.9),
    "GRS 1124-683": (("N", 3.01, 0.15), ("N", 0.128, 0.04), ("I", 54.0, 65.0), 33.2, 80.1),
    "GRS 1915+105": (("N", 9.5, 3.0), ("N", 0.058, 0.033), ("N", 70.0, 2.0), 0.0, 90.0),
    "GS 1354-64": (("N", 5.73, 0.29), ("N", 0.12, 0.04), ("I", 27.2, 80.8), 27.2, 80.8),
    "GS 2000+25": (("N", 5.01, 0.12), ("N", 0.042, 0.012), ("I", 55.0, 65.0), 28.3, 86.7),
    "GS 2023+338": (("N", 6.08, 0.06), ("N", 0.060, 0.005), ("N", 80.1, 5.1), 35.4, 80.0),
    "H1705-250": (("N", 4.86, 0.13), ("U", 0.0, 0.053), ("I", 48.0, 90.0), 0.0, 90.0),
    "SAX J1819.3-2525": (("N", 2.74, 0.12), ("N", 0.67, 0.04), ("I", 80.0, 90.0), 44.8, 69.6),
    "XTE J1118+480": (("N", 6.27, 0.04), ("N", 0.024, 0.009), ("I", 68.0, 82.0), 21.8, 89.4),
    "XTE J1550-564": (("N", 7.65, 0.38), ("U", 0.031, 0.037), ("I", 57.7, 77.1), 26.5, 82.0),
    "XTE J1650-500": (("N", 2.73, 0.56), ("U", 0.0, 0.5), ("N", 75.2, 5.9), 0.0, 90.0),
    "XTE J1859+226": (("N", 4.5, 0.6), ("U", 0.0, 0.5), ("N", 60.0, 3.0), 0.0, 90.0),
}

#: BlackCAT name -> Kreidberg et al. 2012 name, to borrow i_min (donor of
#: normal mass) when BlackCAT gives only an upper limit on the inclination.
KREIDBERG_NAME = {
    "4U 1543-475": "4U 1543-47", "3A 0620-003": "A0620-00", "GRO J0422+32": "GRO J0422+32",
    "GRO J1655-40": "GRO J1655-40", "GRS 1009-45": "GRS 1009-45",
    "GRS 1124-684": "GRS 1124-683", "GRS 1915+105": "GRS 1915+105",
    "GS 1354-64": "GS 1354-64", "GS 2000+251": "GS 2000+25", "GS 2023+338": "GS 2023+338",
    "H 1705-250": "H1705-250", "SAX J1819.3-2525": "SAX J1819.3-2525",
    "XTE J1118+480": "XTE J1118+480", "XTE J1550-564": "XTE J1550-564",
    "XTE J1650-500": "XTE J1650-500", "XTE J1859+226": "XTE J1859+226",
}

#: Common names of BlackCAT entries (the catalogue uses discovery designations).
COMMON_NAME = {
    "GS 2023+338": "V404 Cyg", "3A 0620-003": "A0620-00", "1H 1659-487": "GX 339-4",
    "SAX J1819.3-2525": "V4641 Sgr", "GRS 1124-684": "GU Mus (Nova Mus 1991)",
    "H 1705-250": "Nova Oph 1977", "GS 2000+251": "QZ Vul",
}

# ---------------------------------------------------------------------------
# Published results we check our numbers against
# ---------------------------------------------------------------------------

CHECKS = {
    "farr2011_m1pct": {  # Table 8: 10/50/90% quantiles of M_1%, low-mass sample
        "power_law": (4.3, 6.1, 6.6), "gaussian": (2.9, 4.4, 5.5), "exponential": (5.3, 6.0, 6.5),
        "ref": "Farr et al. 2011, ApJ 741, 103, Table 8",
    },
    "kreidberg2012": "Kreidberg et al. 2012: with their Table 2, significant probability "
                     "for M_1% < 3 M_sun; without GRO J0422+32, M_1% > 4.3 at 90%",
    "rocha2026_m1pct_90": {"combined_power_law": 4.06, "lmxb_gaussian": 2.56,
                           "ref": "Rocha et al. 2026, arXiv:2609.14154"},
    "ozel2010": {"mean": 7.8, "sigma": 1.2, "ref": "Ozel et al. 2010, ApJ 725, 1918"},
    "v404_fm": (6.08, 0.06, "Casares & Charles 1994, MNRAS 271, L5"),
    "cygx1_mass": (21.2, 2.2, "Miller-Jones et al. 2021"),
    "j1820_peak_ledd": (0.15, 0.03, "Atri et al. 2020, MNRAS 493, L81 (for 10 M_sun)"),
    "j1820_soft_to_hard_ledd": (0.03, 0.01, "Atri et al. 2020, from a bolometric flux of "
                                             "(2.5 +- 0.4)e-8 erg/cm2/s"),
    "j1820_soft_to_hard_flux_cgs": (2.5e-8, 0.4e-8),
    "maccarone2003": "soft-to-hard transitions at 1-4% of L_Edd, mean 2% "
                     "(Maccarone 2003, A&A 409, 697)",
    "vahdat_motlagh2019": (0.0158, 0.0093, "Vahdat Motlagh et al. 2019, as quoted by Atri et al. 2020"),
    "atri2020_mass": (9.2, 1.3),
}
