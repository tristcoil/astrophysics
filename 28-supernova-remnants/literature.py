"""Published results this project compares against, each with its source.

Nothing here feeds a measurement: these numbers only grade ours. Values were
read from the papers themselves (arXiv versions), not from summaries.
"""

from __future__ import annotations

# --------------------------------------------------------------------------
# The Crab
# --------------------------------------------------------------------------

#: The "guest star" was first recorded on 1054 July 4 (Song dynasty records;
#: Stephenson & Green 2002, Historical Supernovae and their Remnants).
CRAB_SN_YEAR = 1054.5

#: Convergence dates of the Crab's filaments: (date, error, method, reference).
#: "fit" = proper motion against position (WM77, N98, Martin et al.);
#: "closest approach" = run positions back and find when they huddle closest.
CRAB_DATES = [
    ("Duncan 1939", 1172.0, None, "first estimate", "Duncan 1939, ApJ 89, 482"),
    ("Trimble 1968", 1140.0, 15.0, "closest approach", "Trimble 1968, AJ 73, 535"),
    ("Wyckoff & Murray 1977", 1120.0, 7.0, "fit (Trimble's data)", "Wyckoff & Murray 1977, MNRAS 180, 717"),
    ("Nugent 1998", 1130.0, 16.0, "fit", "Nugent 1998, PASP 110, 831"),
    ("Rudie et al. 2008", 1055.0, 24.0, "northern jet only", "Rudie, Fesen & Yamada 2008, MNRAS 384, 1200"),
    ("Bietenholz & Nugent 2015", 1091.0, 34.0, "Bayesian (Nugent's data)", "Bietenholz & Nugent 2015, MNRAS 454, 2416"),
    ("Martin et al. 2025", 1105.5, 0.5, "fit, 19,974 vectors", "Martin et al. 2025, MNRAS (arXiv:2502.19632)"),
]
#: Bietenholz & Nugent (2015) from VLA radio images of the synchrotron nebula.
CRAB_SYNCHROTRON_DATE = (1255.0, 27.0)

#: Martin et al. (2025), Table 5: their expansion centre relative to Duncan's
#: (1921) reference star (Gaia DR3 3403818176867563264), arcsec (east, north),
#: and their closest-approach date; uncertainty 1.2" on the centre.
MARTIN2025 = {
    "T_fit": 1105.5, "T_fit_err": 0.5,
    "T_closest": 1140.5, "T_closest_err": 1.4,
    "centre_from_duncan": (6.9, -9.3), "centre_from_duncan_err": (0.8, 0.9),
    "pulsar_1054_from_duncan_gaia": (7.8, -6.02),
    "n_vectors": 19974, "epochs": (2007.28, 2016.10, 2019.52),
}

#: Kaplan, Chatterjee, Gaensler & Anderson (2008, ApJ 677, 1201): HST proper
#: motion of the Crab pulsar, relative to the pulsar's surroundings.
KAPLAN2008_PM = (-11.7, 4.2)  # mas/yr, +-0.4 (stat) +-0.5 (sys) each

# --------------------------------------------------------------------------
# Cas A
# --------------------------------------------------------------------------

#: Thorstensen, Fesen & van den Bergh (2001, AJ 122, 297): proper motions of
#: optical knots on plates and CCD images from 1951 to 1999.
THORSTENSEN2001 = {
    "T_outer_knots": 1671.3, "T_outer_err": 0.9,  # 17 fast knots outside the blast wave
    "T_main_shell": 1662.3, "T_main_err": 1.7,  # 21 bright main-shell knots (decelerated)
    "centre_ra": "23 23 27.77", "centre_dec": "+58 48 49.4", "centre_err_arcsec": 0.4,
    "cco_offset_arcsec": 6.6, "cco_offset_err": 1.5,
    "centre_pa_from_cco_deg": 354.0,  # i.e. the neutron star lies at PA ~174 deg from the centre
    "cco_speed_kms": 330.0,  # at 3.4 kpc
}

#: Fesen et al. (2006, ApJ 645, 283): HST ACS proper motions of 1825 outer knots.
FESEN2006 = {
    "T_all": 1662.0, "T_all_err": 27.0,
    "T_bright72": 1672.0, "T_bright72_err": 18.0,
    "T_northwest": 1681.0, "T_northwest_err": 19.0,
    "cco_offset_arcsec": 7.0, "cco_pa_deg": 169.0, "cco_pa_err": 8.4, "cco_speed_kms": 350.0,
}

#: Holland-Ashford, Slane & Long (2024, ApJ; arXiv:2310.19879): the neutron
#: star's proper motion from 2000-2023 Chandra data. Method 1 uses the
#: archive's astrometry as is; method 2 re-registers 13 observations on Gaia
#: stars. Heliocentric, at 3.4 kpc.
HOLLAND_ASHFORD2024 = {
    "method1_kms": 275.0, "method1_err": 121.0, "method1_pa": 177.0, "method1_pa_err": 22.0,
    "method2_kms": 436.0, "method2_err": 89.0, "method2_pa": 158.0, "method2_pa_err": 12.0,
    "previous_10yr_kms": 570.0, "previous_10yr_err": 260.0,
}

#: Vink, Patnaude & Castro (2022, ApJ 929, 57), Table 3: forward-shock
#: expansion rates (%/yr) in 20-degree sectors, 4.2-6 keV Chandra images
#: 2000-2019, about Thorstensen's centre; m for an age of 333 yr.
VINK2022_FS = {
    "pa": [10, 30, 50, 70, 90, 110, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350],
    "rate": [0.1900, 0.2760, 0.1990, 0.2098, 0.2217, 0.2212, 0.2384, 0.2707, 0.2329, 0.1518,
             0.2076, 0.2331, 0.1934, 0.2251, 0.2376, 0.2274, 0.1963, 0.1958],
    "rate_err": [0.0012, 0.0009, 0.0012, 0.0011, 0.0016, 0.0016, 0.0008, 0.0020, 0.0024, 0.0014,
                 0.0030, 0.0017, 0.0013, 0.0009, 0.0010, 0.0008, 0.0016, 0.0020],
    "mean_rate": 0.218, "mean_rate_sd": 0.029, "mean_m": 0.727, "mean_v_kms": 5803.0,
    "radius_arcmin": 2.8, "age_yr": 333.0,
}
#: Vink et al. (2022), Table 4: the reverse-shock region (annulus 1.67'-2.18').
VINK2022_RS = {
    "pa": [10, 30, 50, 70, 90, 110, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350],
    "rate": [0.2275, 0.1364, 0.2140, 0.1814, 0.1663, 0.1748, 0.2053, 0.1704, 0.1479, 0.1188,
             0.1090, 0.0929, 0.0679, -0.1023, -0.0226, 0.1062, 0.1865, 0.2285],
    "mean_rate": 0.134, "mean_rate_sd": 0.084,
}
#: DeLaney & Rudnick (2003, ApJ 589, 818): the first Chandra measurement
#: (2000 vs 2002, 0.3-10 keV): forward-shock filaments 0.02-0.33 %/yr, median
#: 0.21; the bright ring's expansion from Einstein and ROSAT, 0.21 %/yr.
DELANEY2003 = {"fs_median": 0.21, "fs_min": 0.02, "fs_max": 0.33, "bright_ring": 0.21}

# --------------------------------------------------------------------------
# Neutron-star kicks
# --------------------------------------------------------------------------

#: Hobbs, Lorimer, Lyne & Kramer (2005, MNRAS 360, 974): 233 pulsar proper
#: motions. Mean 2-D speeds (km/s): normal (not recycled) pulsars 246 +- 22,
#: characteristic age < 3 Myr 307 +- 47, recycled 87 +- 13; the 3-D birth
#: speeds follow a Maxwellian with 1-D sigma = 265 km/s (mean 3-D 400 +- 40).
HOBBS2005 = {"v2d_normal": 246.0, "v2d_normal_err": 22.0, "v2d_young": 307.0, "v2d_young_err": 47.0,
             "v2d_recycled": 87.0, "v2d_recycled_err": 13.0, "sigma_1d": 265.0, "v3d_mean": 400.0}

# --------------------------------------------------------------------------
# Historical supernovae and their remnants
# --------------------------------------------------------------------------

#: (event, remnant, Green's name, year, distance kpc, distance range, distance
#: reference, note). Years are when the light reached us. The SN 393 and
#: SN 1181 identifications are disputed and drawn hollow.
HISTORICAL = [
    ("SN 185", "RCW 86", "G315.4-2.3", 185, 2.5, (2.3, 2.8),
     "Helder et al. 2013, ApJ 764, 11", "Type Ia; identification widely accepted"),
    ("SN 393?", "RX J1713.7-3946", "G347.3-0.5", 393, 1.0, (0.9, 1.3),
     "Fukui et al. 2003, PASJ 55, L61", "identification doubtful"),
    ("SN 1006", "SN 1006", "G327.6+14.6", 1006, 2.18, (2.10, 2.26),
     "Winkler, Gupta & Long 2003, ApJ 585, 324", "Type Ia"),
    ("SN 1054", "Crab Nebula", "G184.6-5.8", 1054, 2.0, (1.7, 2.3),
     "Trimble 1973; Gaia DR3 parallax", "the visible nebula is the pulsar wind nebula, not the blast wave"),
    ("SN 1181?", "3C 58", "G130.7+3.1", 1181, 2.0, (2.0, 3.2),
     "Kothes 2013, A&A 560, A18", "disputed: Pa 30 is now favoured (Ritter et al. 2021, ApJL 918, L33)"),
    ("SN 1572", "Tycho", "G120.1+1.4", 1572, 3.0, (2.3, 4.0),
     "Tian & Leahy 2011 (2.5 +- 0.5); Hayato et al. 2010 (4 +- 1)", "Type Ia (light echo: Krause et al. 2008)"),
    ("SN 1604", "Kepler", "G4.5+6.8", 1604, 5.1, (4.4, 5.9),
     "Sankrit et al. 2016, ApJ 817, 36", "Type Ia"),
    ("~1671", "Cas A", "G111.7-2.1", 1671, 3.4, (3.3, 3.7),
     "Reed et al. 1995, ApJ 440, 706", "Type IIb (light echo: Krause et al. 2008); no confirmed sighting"),
    ("~1900", "G1.9+0.3", "G1.9+0.3", 1900, 8.5, (8.0, 9.0),
     "Reynolds et al. 2008; Carlton et al. 2011, ApJL 737, L22", "hidden by dust near the Galactic centre"),
]
