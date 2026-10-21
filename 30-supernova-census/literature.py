"""Published numbers this project is graded against, each with its source.

Rates are in supernovae per Gpc^3 per year for H0 = 70 (h70^3), the
convention of all the surveys quoted. Values were checked against the papers
themselves (arXiv versions) when the project was built, October 2026.
"""

# --- Volumetric rates ------------------------------------------------------

#: Perley et al. 2020, ApJ 904, 35 (BTS II), Sect. 4.3. 95% intervals, stat only.
PERLEY_2020 = {
    "rate_Ia": (2.35e4, 0.24e4),               # M < -16.5
    "rate_CC": (10.1e4, 5.0e4, 3.5e4),         # M < -14; +, - (95%)
    "rate_SLSN": (5.6, 5.4, 2.8),              # M < -21
    "f_rec": (412, 690),
    "completeness_m18p5": 0.93,
    "red_sequence_hosts": {"Ia": (0.31, 0.02), "CC": (0.11, 0.02)},   # Table 2
    "n_Ia_passing_m18p5": 875, "n_CC_passing_m18p5": 313,             # Table 1
}
#: Frohmaier et al. 2019, MNRAS 486, 2308 (PTF, z <= 0.09): normal SNe Ia.
FROHMAIER_2019_IA = (2.43e4, 0.29e4)
#: Frohmaier et al. 2021, MNRAS 500, 5142 (PTF, <z> = 0.028).
FROHMAIER_2021 = {"rate_CC": (9.10e4, 1.56e4, 1.27e4), "rate_SESN": (2.41e4, 0.81e4, 0.64e4),
                  "rate_SLSN_I": (35.0, 25.0, 13.0)}
#: Desai et al. 2024, MNRAS 530, 5016 (ASAS-SN 2014-2017, M_V < -16.5, z = 0.024)
#: and Desai et al. 2026, arXiv:2602.00223 (ASAS-SN, a decade, z = 0.029).
DESAI_2024_IA = (2.28e4, 0.20e4)
DESAI_2026_IA = (2.55e4, 0.12e4)
#: Pessi et al. 2025, A&A 703, A34 (ASAS-SN II, 2014-2017, M_V < -14, z = 0.0149).
PESSI_2025 = {"rate_CC": (7.0e4, 1.0e4, 0.9e4), "frac_II": 0.632, "frac_SE": 0.323,
              "frac_interacting": 0.044, "rate_SLSN": (1.5, 4.4, 1.1)}
#: Das et al. 2025, arXiv:2502.19493 (ZTF Census of the Local Universe).
DAS_2025 = {"rate_IIP": (3.9e4, 0.4e4), "rate_LLIIP": (7.3e3, 0.6e3)}

# --- Fractions (Li et al. 2011, MNRAS 412, 1441, LOSS; Table 7, in %) ------

LI_2011_VOLUME = {"Ia": (24.1, 3.7, 3.5), "Ibc": (18.7, 3.5, 3.3), "II": (57.2, 4.3, 4.1)}
LI_2011_MAGLIM = {"Ia": (79.2, 4.2, 5.5), "Ibc": (4.1, 1.6, 1.3), "II": (16.6, 5.0, 3.9)}
#: within Type II, volume-limited: II-P 69.9, II-L 9.7, IIb 11.9, IIn 8.6
LI_2011_WITHIN_II = {"IIP": 69.9, "IIL": 9.7, "IIb": 11.9, "IIn": 8.6}
#: Table 6: mean observed peak absolute magnitudes (LOSS unfiltered ~ R band)
LI_2011_MEAN_M = {"Ia": (-18.49, 0.76), "Ibc": (-16.09, 1.24), "II": (-16.05, 1.37),
                  "IIP": (-15.66, 1.23), "IIL": (-17.44, 0.64), "IIb": (-16.65, 1.30),
                  "IIn": (-16.86, 1.61)}

# --- Cosmic star formation --------------------------------------------------

def madau_dickinson_sfrd(z):
    """Madau & Dickinson (2014) Eq. 15, Salpeter IMF (0.1-100 M_sun),
    M_sun / yr / Mpc^3."""
    return 0.015 * (1 + z) ** 2.7 / (1 + ((1 + z) / 2.9) ** 5.6)

#: Madau & Dickinson's conversion between IMFs: SFR(Kroupa) = 0.67 SFR(Salpeter).
KROUPA_OVER_SALPETER_SFR = 0.67
#: Project J3 (24-galaxy-life-cycle): SDSS SFR density at z = 0.1 (Kroupa) and
#: log stellar mass density (MPA-JHU masses), both h = 0.7.
J3_SFRD_KROUPA = (0.0145, 0.0002)
J3_LOG_RHO_STAR = 8.588

#: Core-collapse supernovae hidden from optical surveys by dust in their hosts:
#: Mattila et al. 2012, ApJ 756, 111 (~19% locally); Jencson et al. 2019,
#: ApJ 886, 40 (SPIRITS, 38.5%, 90% CI 16.6-64.5%).
MATTILA_2012_MISSED = 0.19
JENCSON_2019_MISSED = (0.385, 0.166, 0.645)

# --- Type Ia hosts ("A + B") -----------------------------------------------

#: Sullivan et al. 2006, ApJ 648, 868 (SNLS): A in SNe/yr/M_sun, B in SNe/yr/(M_sun/yr)
SULLIVAN_2006_AB = {"A": (5.3e-14, 1.1e-14), "B": (3.9e-4, 0.7e-4)}

# --- Type II light curves ---------------------------------------------------

#: Anderson et al. 2014, ApJ 786, 67: optically thick phase (explosion to end
#: of plateau), V band, 116 SNe II: mean 83.7 d, sigma 16.7 d.
ANDERSON_2014_OPTD = (83.7, 16.7)
#: Anderson 2019, A&A 628, A7: median 56Ni mass of 115 literature SNe II.
ANDERSON_2019_NI_MEDIAN = 0.032
#: Martinez et al. 2022, A&A 660, A41 (CSP-I, 53 SNe II, hydrodynamical fits):
#: ejecta 7.9-14.8 M_sun, energies 0.15-1.40 foe, 56Ni 0.006-0.069 M_sun;
#: initial-mass distribution of the gold sample cut off at 21.3 (+3.8/-0.4).
MARTINEZ_2022 = {"M_ej_range": (7.9, 14.8), "E51_range": (0.15, 1.40),
                 "Ni_range": (0.006, 0.069), "M_upper": (21.3, 3.8, 0.4)}

# --- Progenitors -------------------------------------------------------------

#: Smartt 2015, PASA 32, e016.
SMARTT_2015 = {"n_expected": 13, "p": 3e-5, "M_hi_STARS": (16.5, 2.5, 2.5),
               "M_lo_STARS": (9.5, 0.5, 2.0), "M_hi_KEPLER": (18.5, 3.0, 4.0),
               "M_lo_KEPLER": (10.0, 0.5, 1.5), "logL_limit": 5.1}
#: Davies & Beasor 2018, MNRAS 474, 2116: M_hi with their bolometric corrections.
DAVIES_BEASOR_2018 = {"M_hi": (19.0, 2.5, 1.3), "M_hi_95_upper": 27.0,
                      "M_hi_with_systematics": 25.0, "M_hi_with_systematics_95": 33.0}
#: Davies & Beasor 2020, MNRAS 493, 468: luminosity-space fit.
DAVIES_BEASOR_2020 = {"logL_hi_free": (5.20, 0.17, 0.11), "logL_hi_fixed_slope": (5.28, 0.12, 0.06),
                      "significance": "1-2 sigma", "M_hi": 19.2}
#: Kochanek 2020, MNRAS 493, 4945: Bayesian M_hi for the progenitor sample.
KOCHANEK_2020 = {"M_hi_eldridge": (19.01, 4.04, 2.04), "M_hi_sukhbold": (21.28, 4.52, 2.28)}
#: Beasor, Smith & Jencson 2024, arXiv:2410.14027.
BEASOR_2024 = {"offset_dex": 0.3, "extra_scatter_dex": 0.2,
               "verdict": "most luminous RSGs not ruled out at 3 sigma"}
#: Fang, Moriya & Maeda 2025, arXiv:2504.14502: nebular [O I] masses.
FANG_2025 = {"logL_cut": (5.21, 0.09, 0.07), "significance": "2-3 sigma",
             "M_upper": (20.63, 2.42, 1.64)}
