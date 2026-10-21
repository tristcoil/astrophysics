"""Published numbers: the answer keys we grade ourselves against, and the few we adopt.

Every value was read from the paper itself (TeX source, or the PDF's text for
Barstow et al. 2005); table or section given. A value is an *input* only where
a function says it adopts it (Sirius A's orbit size, the parallaxes, the
spectroscopic temperature of Sirius B).
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Sirius (Bond et al. 2017, ApJ 840, 70; arXiv:1703.10625)
# ---------------------------------------------------------------------------

#: Table 4: elements of the relative visual orbit (J2000). T0 is a Besselian year.
BOND17_ELEMENTS = {
    "P": (50.1284, 0.0043), "a": (7.4957, 0.0025), "i": (136.336, 0.040), "Omega": (45.400, 0.071),
    "T_besselian": (1994.5715, 0.0058), "e": (0.59142, 0.00037), "omega": (149.161, 0.075),
}
#: Table 5: absolute parallax (arcsec). Ground-based compilation (Gatewood &
#: Gatewood 1978), Hipparcos (van Leeuwen 2007), and their weighted mean.
BOND17_PLX = {
    "GG78": {"value": 0.3777, "err": 0.0031, "ref": "Gatewood & Gatewood 1978 (ground-based compilation)"},
    "HIP2": {"value": 0.37921, "err": 0.00158, "ref": "Hipparcos, van Leeuwen 2007"},
    "adopted": {"value": 0.3789, "err": 0.0014, "ref": "Bond et al. 2017 weighted mean"},
}
SIRIUS_PLX_ADOPTED = BOND17_PLX["adopted"]
#: Table 5: semi-major axis of Sirius A's orbit about the barycentre (arcsec).
#: GG78 measured it on Yerkes and Allegheny plates; Bond et al. refitted GG78's
#: photocentre positions (their Table 6) with the new relative orbit.
BOND17_AA = {"value": 2.4761, "err": 0.0045, "ref": "Bond et al. 2017 Table 5 (GG78 data refitted)"}
GG78_AA = {"value": 2.4904, "err": 0.0040, "ref": "Gatewood & Gatewood 1978"}
#: Table 6: dynamical masses (Msun).
BOND17_MASSES = {"Mtot": (3.081, 0.034), "MA": (2.063, 0.023), "MB": (1.018, 0.011)}
GG78_MASSES = {"Mtot": (3.196, 0.083), "MA": (2.143, 0.056), "MB": (1.053, 0.028)}
VDB60_MASSES = {"Mtot": (3.20, None), "MA": (2.15, None), "MB": (1.05, None)}
#: Sect. 6.2 footnote: orbit from the HST astrometry alone.
BOND17_HST_ONLY = {"MA": (2.068, 0.062), "MB": (1.020, 0.030)}
#: Sect. 5.2: reduced chi^2 of the HST measures in their final fit.
BOND17_HST_CHI2NU = 1.87
#: Table 7: error budget sigma(M_A), sigma(M_B) per input.
BOND17_BUDGET = {"parallax": (0.023, 0.011), "a": (0.0024, 0.0007), "a_A": (0.0018, 0.0018),
                 "P": (0.0004, 0.0002)}
#: Sect. 5.1: systemic radial velocity -7.70 km/s; -8.47 km/s once corrected
#: for Sirius A's gravitational redshift (the true space-motion component).
SIRIUS_RV_SYSTEM = {"observed": -7.70, "true": -8.47}
#: Sect. 7 (Barstow et al. 2017, conference paper): Sirius B from STIS spectra.
BOND17_WD = {"Teff": (25369.0, 46.0), "logg": (8.591, 0.016), "R": (0.008098, 0.000046),
             "L": (0.02448, 0.00033), "logg_dyn": (8.629, 0.007), "vgr_pred": (79.8, 1.0),
             "cooling_age_myr": 126.0}

#: Joyce et al. 2018 (MNRAS 481, 2361; arXiv:1809.01240): gravitational redshift
#: of Sirius B measured against Sirius A's H-alpha, and the mass it implies.
JOYCE18 = {"vgr": (80.65, 0.77), "M_gr": (1.017, 0.025), "Teff_G430L": 25922.0,
           "plx_dr2": (376.68e-3, 0.45e-3), "R_dr2": (0.00808, 0.00011), "M_gr_dr2": (1.0246, 0.0268)}

#: Barstow et al. 2005 (MNRAS 362, 1134; astro-ph/0506600) Table 2: V from the
#: HST STIS G430L spectrum vs the older ground-based estimate (Holberg et al.
#: 1984, used by Holberg et al. 1998); Teff and log g from the Balmer lines.
BARSTOW05 = {"V_hst": (8.528, 0.05), "V_holberg84": (8.44, 0.06), "Teff": (25193.0, 37.0),
             "logg": (8.556, 0.010), "MV": (11.427, 0.05)}

#: Hipparcos (1997) DMSA/O orbital solution for Sirius A's photocentre: elements
#: adopted from Gatewood & Gatewood 1978; the main-catalogue proper motion then
#: refers to the barycentre. P in days, T in JD - 2440000, a0 in mas.
HIP1_DMSA_O_SIRIUS = {"P_d": 18295.4, "T_jd2440000": -27123.7, "a0_mas": 2490.40, "e": 0.5923,
                      "omega": 327.27, "i": 136.53, "Omega": 44.86}
#: Hipparcos (van Leeuwen 2007 keeps the 1997 solution): position at J1991.25,
#: parallax and proper motion (mas, mas/yr).
HIP2_SIRIUS = {"ra": 101.28854105, "dec": -16.71314306, "plx": (379.21, 1.58),
               "pmra": (-546.01, 1.33), "pmdec": (-1223.07, 1.24), "epoch": 1991.25}

# ---------------------------------------------------------------------------
# Procyon (Bond et al. 2015, ApJ 813, 106; arXiv:1510.00485)
# ---------------------------------------------------------------------------

BOND15_PROCYON = {"plx": (0.2850, 0.0007), "aA": (1.232, 0.008),
                  "Mtot": (2.070, 0.016), "MA": (1.478, 0.012), "MB": (0.592, 0.006)}

# ---------------------------------------------------------------------------
# alpha Centauri (Akeson et al. 2021, AJ 162, 14; arXiv:2104.10086), table of orbital elements
# ---------------------------------------------------------------------------

AKESON21_ACEN = {"plx": (750.81e-3, 0.38e-3), "MA": (1.0788, 0.0029), "MB": (0.9092, 0.0025)}

# ---------------------------------------------------------------------------
# White dwarfs with measured masses and radii
# ---------------------------------------------------------------------------

#: Bond et al. 2017b (ApJ 848, 16; arXiv:1709.00478) Table 2: white dwarfs in
#: nearby visual binaries (Stein 2051 B: mass from astrometric microlensing).
#: (Teff, err), (R, err), (M, err), spectral type
VISUAL_BINARY_WDS = {
    "Procyon B": {"Teff": (7740, 50), "R": (0.01232, 0.00032), "M": (0.592, 0.006), "sp": "DQZ"},
    "Sirius B": {"Teff": (25369, 46), "R": (0.008098, 0.000046), "M": (1.018, 0.011), "sp": "DA2"},
    "Stein 2051 B": {"Teff": (7122, 181), "R": (0.0114, 0.0004), "M": (0.675, 0.051), "sp": "DC"},
    "40 Eri B": {"Teff": (17200, 110), "R": (0.01308, 0.00020), "M": (0.573, 0.018), "sp": "DA2.9"},
}

# ---------------------------------------------------------------------------
# Eker et al. 2018 (MNRAS 479, 5491; arXiv:1807.02568): eclipsing-binary relations
# ---------------------------------------------------------------------------

#: Table 4: six-piece mass-luminosity relation, log L = alpha log M + b
#: (mass range lower bound, upper bound, alpha, b, scatter sigma in log L).
EKER18_MLR = (
    (0.179, 0.45, 2.028, -0.976, 0.076),
    (0.45, 0.72, 4.572, -0.102, 0.109),
    (0.72, 1.05, 5.743, -0.007, 0.129),
    (1.05, 2.40, 4.329, 0.010, 0.140),
    (2.40, 7.0, 3.967, 0.093, 0.165),
    (7.0, 31.0, 2.865, 1.105, 0.152),
)
#: Table 5: R = 0.438 M^2 + 0.479 M + 0.075 (M <= 1.5 Msun);
#: log Teff = -0.170 (log M)^2 + 0.888 log M + 3.671 (M > 1.5 Msun).
EKER18_MRR = (0.438, 0.479, 0.075)
EKER18_MTR = (-0.170, 0.888, 3.671)
