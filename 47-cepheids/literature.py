"""Published numbers we grade ourselves against (answer keys), each with its source.

Every value was checked in the paper itself (table, equation or section given).
None of these is an input to a fit, except where a function says it *adopts*
one (the eclipsing-binary distances used to measure the metallicity term, and
the mean metallicities of the Clouds).
"""

# ---------------------------------------------------------------------------
# Geometric distances to the Clouds (detached eclipsing binaries)
# ---------------------------------------------------------------------------

#: Pietrzynski et al. 2019 (Nature 567, 200): 20 eclipsing binaries.
LMC_DEB = {"mu": 18.477, "stat": 0.004, "sys": 0.026, "d_kpc": 49.59, "d_stat": 0.09, "d_sys": 0.54,
           "centre": (80.05, -69.30), "plane_incl": 25.0, "plane_incl_err": 4.0, "plane_pa": 132.0,
           "plane_pa_err": 10.0, "ref": "Pietrzyński et al. 2019, Nature 567, 200"}

#: Graczyk et al. 2020 (ApJ 904, 13): 15 eclipsing binaries; Sect. 6 gives the
#: difference from the LMC. The binaries span ~10 kpc in depth; NE nearer.
SMC_DEB = {"mu": 18.977, "stat": 0.016, "sys": 0.028, "d_kpc": 62.44, "d_stat": 0.47, "d_sys": 0.81,
           "dmu_smc_lmc": 0.500, "dmu_err": 0.017, "centre": (12.54, -73.11),
           "ref": "Graczyk et al. 2020, ApJ 904, 13"}

# ---------------------------------------------------------------------------
# Mean Cepheid metallicities of the Clouds (adopted inputs)
# ---------------------------------------------------------------------------

#: Romaniello et al. 2022 (A&A 658, A29), 89 LMC Cepheids, HR spectra:
#: [Fe/H] = -0.409 +- 0.003, dispersion 0.076 (+0.1 dex systematic on the mean).
FEH_LMC = {"feh": -0.409, "err": 0.003, "disp": 0.076, "sys": 0.1, "ref": "Romaniello et al. 2022"}
#: SMC: -0.75 +- 0.05 as adopted by Breuval et al. 2021, 2022 (from Gieren et al.
#: 2018 / Romaniello et al. 2008). Breuval et al. 2022 warn a reanalysis on the
#: 2022 temperature scale could give ~ -0.90; Breuval et al. 2024 adopt -0.785.
FEH_SMC = {"feh": -0.75, "err": 0.05, "alt": (-0.785, -0.90), "ref": "Breuval et al. 2021/2022 (Gieren et al. 2018)"}

# ---------------------------------------------------------------------------
# The Leavitt law in the LMC (OGLE-IV), for our reproduction
# ---------------------------------------------------------------------------

#: Soszynski et al. 2015 (AcA 65, 297) Table 2: the OGLE-IV fits on the same
#: data we use, m = slope log P + zp (intercept at log P = 0), iterative 3-sigma
#: clipping, all stars of a mode, no period cut. (slope, zp, sigma)
S15_TABLE2 = {
    ("lmc", "F", "W_I"): (-3.314, 15.888, 0.077), ("lmc", "F", "I"): (-2.911, 16.822, 0.146),
    ("lmc", "F", "V"): (-2.690, 17.438, 0.208),
    ("lmc", "1O", "W_I"): (-3.431, 15.393, 0.081), ("lmc", "1O", "I"): (-3.260, 16.362, 0.162),
    ("lmc", "1O", "V"): (-3.142, 16.979, 0.227),
    ("smc", "F", "W_I"): (-3.460, 16.493, 0.155), ("smc", "F", "I"): (-3.115, 17.401, 0.215),
    ("smc", "F", "V"): (-2.898, 17.984, 0.266),
    ("smc", "1O", "W_I"): (-3.548, 15.961, 0.169), ("smc", "1O", "I"): (-3.299, 16.818, 0.222),
    ("smc", "1O", "V"): (-3.155, 17.368, 0.271),
}
S15_SLOPE_ERR = {("lmc", "F", "W_I"): 0.008, ("lmc", "1O", "W_I"): 0.007, ("smc", "F", "W_I"): 0.011,
                 ("smc", "1O", "W_I"): 0.017, ("lmc", "F", "I"): 0.014, ("lmc", "F", "V"): 0.018}
#: Jacyszyn-Dobrzeniecka et al. 2016 (AcA 66, 149) Table 2, W_I, F split at log P = 0.4.
JD16_WI = {("lmc", "F>0.4"): (-3.317, 15.890, 0.075, 2103), ("smc", "F>0.4"): (-3.315, 16.379, 0.144, 957),
           ("lmc", "1O"): (-3.414, 15.388, 0.079, 1931), ("smc", "1O"): (-3.540, 15.959, 0.170, 2010)}
#: Is the law straight? Ngeow et al. 2009 (ApJ 693, 691) Table 3, F-test of a
#: break at 10 d in OGLE-III LMC: p = 0.042 (V), 0.003 (I), 0.055 (W).
BREAK_10D_NGEOW09 = {"V": 0.042, "I": 0.003, "W": 0.055}

#: Ripepi et al. 2023 (A&A 674, A17) Table 2: LMC relations in Gaia bands
#: (apparent mags; intercept at log P = 0). W = G - 1.90 (BP - RP).
GAIA_LMC_PW = {"F": {"slope": -3.317, "e_slope": 0.007, "icpt": 15.998, "e_icpt": 0.005, "rms": 0.075, "n": 2447},
               "1O": {"slope": -3.476, "e_slope": 0.009, "icpt": 15.521, "e_icpt": 0.003, "rms": 0.081, "n": 1745},
               "ref": "Ripepi et al. 2023, Table 2"}
#: Cruz Reyes & Anderson 2023 (A&A 672, A85) Table 8, LMC W_G at P = 10 d
#: (their W_G uses 1.921, not 1.90).
CRA23_LMC_WG = {"at_logp1": 12.669, "err": 0.005, "slope": -3.347, "e_slope": 0.020, "n": 543}

# ---------------------------------------------------------------------------
# Geometry of the Clouds from Cepheids
# ---------------------------------------------------------------------------

#: Jacyszyn-Dobrzeniecka et al. 2016 Table 6: 4,190 OGLE Cepheids, plane fit.
LMC_GEOMETRY = {
    "jd16": {"incl": 24.2, "e_incl": 0.7, "pa": 151.4, "e_pa": 1.7, "centre": (80.05, -69.30),
             "ref": "Jacyszyn-Dobrzeniecka et al. 2016, Table 6"},
    "inno16": {"incl": 25.05, "e_incl": 0.55, "pa": 150.76, "e_pa": 0.20, "centre": (80.78, -69.30),
               "mu": 18.48, "e_mu": 0.10, "ref": "Inno et al. 2016, ApJ 832, 176, Table 5"},
    "near_side": "north-east",
}
SMC_GEOMETRY = {"near_side": "north-east", "centre": (12.54, -73.11),
                "notes": "line-of-sight extent 4-5x the transverse size (Jacyszyn-Dobrzeniecka et al. 2016); "
                         "east up to 20 kpc closer than west (Scowcroft et al. 2016); elongated 25-30 kpc "
                         "E/NE-SW (Ripepi et al. 2017)"}

# ---------------------------------------------------------------------------
# Gaia parallax offsets for Cepheids (our convention: true = L21-corrected - delta,
# so delta > 0 means the Lindegren-corrected parallaxes are still too large)
# ---------------------------------------------------------------------------

PARALLAX_OFFSET = {
    "riess21": {"delta_uas": 14.0, "err": 6.0, "ref": "Riess et al. 2021 (ApJL 908, L6) Table 2: zp = -14 +- 6 uas"},
    "riess22_clusters": {"delta_uas": 3.0, "err": 4.0, "ref": "Riess et al. 2022 (ApJ 938, 36): cluster stars, -3 +- 4 uas"},
    "cra23": {"delta_uas": 19.0, "err": 3.0, "ref": "Cruz Reyes & Anderson 2023 (A&A 672, A85) eq. 23-25: -19 +- 3 uas (W_G)"},
    "cmetall9": {"delta_uas": 16.0, "err": 4.0, "ref": "Ripepi et al. 2026 (C-MetaLL IX, A&A 708, A216): -16 +- 4 uas"},
}

# ---------------------------------------------------------------------------
# Absolute calibrations in the Gaia Wesenheit and the LMC distances they imply
# ---------------------------------------------------------------------------

GAIA_PWZ = {
    "ripepi22": {"alpha": -5.988, "e_alpha": 0.018, "beta": -3.176, "e_beta": 0.044, "gamma": -0.520,
                 "e_gamma": 0.090, "mu_lmc": 18.503, "e_mu_lmc": 0.046, "n": 435,
                 "mu_lmc_noZ": (18.69, 18.72),
                 "ref": "Ripepi et al. 2022 (A&A 659, A167), Table 2 case 8"},
    "breuval22": {"gamma": -0.384, "e_gamma": 0.051, "ref": "Breuval et al. 2022 (ApJ 939, 89) Table 5, W_G"},
    "trentin24": {"gamma": -0.573, "e_gamma": 0.066, "mu_lmc": 18.436, "e_mu_lmc": 0.028,
                  "ref": "Trentin et al. 2024 (C-MetaLL IV) Table 3 ID 35"},
    "cmetall9": {"gamma": -0.527, "e_gamma": 0.094, "delta_uas": 16.0, "mu_lmc": 18.524, "e_mu_lmc": 0.06,
                 "ref": "Ripepi et al. 2026 (C-MetaLL IX) Table 2"},
    "cra23": {"m_at_logp1_solar": -6.004, "e": 0.019, "slope": -3.242, "gamma_fixed": -0.384,
              "ref": "Cruz Reyes & Anderson 2023 eqs. 23-25 (W_G with 1.921)"},
}

#: Metallicity term in other Wesenheit indices (sign: M = ... + gamma [Fe/H]).
GAMMA_OTHER = {
    "breuval21_WVI": (-0.251, 0.057, "Breuval et al. 2021, W_VI with 1.526"),
    "breuval22_WVI": (-0.201, 0.071, "Breuval et al. 2022, W_VI with 1.387"),
    "riess22_WH": (-0.217, 0.046, "Riess et al. 2022 SH0ES, W_H, per [O/H]"),
    "breuval24_WH": (-0.234, 0.052, "Breuval et al. 2024, W_H"),
    "wielgorski17_WI": (-0.025, 0.067, "Wielgorski et al. 2017, W_I (LMC-SMC with eclipsing binaries)"),
}

# ---------------------------------------------------------------------------
# The Hubble constant and its anchors
# ---------------------------------------------------------------------------

H0 = {
    "riess22": {"H0": 73.04, "err": 1.04, "ref": "Riess et al. 2022 (ApJL 934, L7), three anchors"},
    "riess22_lmc_only": {"H0": 73.59, "err": 1.36, "mu_lmc": 18.477},
    "riess19_lmc_only": {"H0": 74.22, "err": 1.82, "ref": "Riess et al. 2019 (ApJ 876, 85)"},
    "planck18": {"H0": 67.4, "err": 0.5, "ref": "Planck 2018 VI (TT,TE,EE+lowE+lensing)"},
    "atlas_D1": {"H0": 73.55, "err": 1.02, "ref": "this atlas, D1 (05-supernova-cosmology)"},
    "atlas_D3": {"H0": 66.7, "err": 0.9, "ref": "this atlas, D3 (09-planck-cmb)"},
}


# ---------------------------------------------------------------------------
# Leavitt's stars revisited; light-curve physics
# ---------------------------------------------------------------------------

#: Breuval, Huang & Riess 2025 (PASP 137, 044001): Leavitt's 25 stars with
#: OGLE and Gaia. Slope of (max+min)/2 vs log P: -2.033 +- 0.098, scatter 0.224;
#: six periods off by more than 0.01 d: HV 1446, 1506, 1436, 1413, 1646, 821.
BREUVAL25 = {"slope": -2.033, "e_slope": 0.098, "scatter": 0.224,
             "discrepant": (1446, 1506, 1436, 1413, 1646, 821)}
#: HV 821 (BZ Tuc) period history (days): Leavitt 1912; Harvard 1966 (Payne-
#: Gaposchkin & Gaposchkin, via Pietrukowicz 2002); Berdnikov & Turner 1997;
#: ASAS (Pietrukowicz 2002); ASAS 1998-2009 (Breuval et al. 2025); Gaia DR3; OGLE-IV.
HV821_PERIODS = (("Leavitt 1912", 127.0), ("Harvard 1966", 127.78), ("Berdnikov & Turner 1997", 127.447),
                 ("ASAS (Pietrukowicz 2002)", 126.312), ("ASAS 1998-2009", 127.445), ("Gaia DR3", 128.134),
                 ("OGLE-IV 2010-2020", 128.197))
#: Centre of the Hertzsprung progression (bump at maximum light), LMC:
#: 10.5 +- 0.5 d (Beaulieu 1998), 11.2 +- 0.8 d (Welch et al. 1997, MACHO), models
#: 11.24 +- 0.46 d (Bono, Marconi & Stellingwerf 2000); Galaxy 10.0 +- 0.5 d.
HP_CENTRE = {"lmc_beaulieu98": (10.5, 0.5), "lmc_welch97": (11.2, 0.8), "models_bono00": (11.24, 0.46),
             "galaxy": (10.0, 0.5)}
#: Double-mode F/1O period ratio: SMC ~0.01 above the LMC, ~0.01 above the
#: Galaxy (Beaulieu 1998); Galactic median 0.708 (Lemasle et al. 2018, 19 stars).
PETERSEN = {"mw_lemasle18": 0.708}
