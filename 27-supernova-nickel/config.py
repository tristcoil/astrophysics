"""Project configuration: paths, constants, nuclear data and per-supernova inputs.

Everything we *assume* lives here, each value with its source: distances,
dust, explosion dates, nuclear decay data. The published nickel masses we
grade ourselves against are kept separately in ``literature.py`` and are
never used as inputs to a fit.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
FIGURE_DIR = PROJECT_DIR / "figures"
CACHE_DIR = DATA_DIR / "cache"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, FIGURE_DIR, CACHE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Physical constants (cgs)
# ---------------------------------------------------------------------------

C_CGS = 2.99792458e10  # speed of light, cm/s
MSUN_G = 1.98841e33  # solar mass, g (IAU 2015 nominal GM_sun / CODATA G)
LSUN = 3.828e33  # nominal solar luminosity, erg/s (IAU 2015 B3)
MBOL_SUN = 4.74  # absolute bolometric magnitude of the Sun (IAU 2015 B2)
MEV_ERG = 1.602176634e-6  # 1 MeV in erg
AMU_G = 1.66053906660e-24  # atomic mass unit, g
DAY_S = 86400.0
PC_CM = 3.0856775814913673e18
MPC_CM = 1e6 * PC_CM
SIGMA_T = 6.6524587321e-25  # Thomson cross-section, cm^2
ME_C2_KEV = 510.99895  # electron rest energy, keV
H_CGS = 6.62607015e-27  # Planck constant, erg s
K_CGS = 1.380649e-16  # Boltzmann constant, erg/K
SIGMA_SB = 5.670374419e-5  # Stefan-Boltzmann constant, erg/cm^2/s/K^4

# ---------------------------------------------------------------------------
# Nuclear data for the 56Ni -> 56Co -> 56Fe chain
# ---------------------------------------------------------------------------

#: Half-lives from the NNDC / ENSDF evaluation (Junde, Huang & Tuli 2011,
#: Nuclear Data Sheets 112, 1513), in days.
HALF_LIFE_NI56_D = 6.075
HALF_LIFE_CO56_D = 77.236

#: Mass of one 56Ni atom, in atomic mass units (AME2020).
MASS_NI56_AMU = 55.942128

#: Energy released per decay that can heat the ejecta, in MeV (Nadyozhin 1994,
#: ApJS 92, 527, Table 1; the standard values in the supernova literature).
#: Neutrinos escape freely and are left out.
#:   56Ni: electron capture only; 1.75 MeV in gamma rays per decay.
#:   56Co: 81% electron capture, 19% positron emission; 3.61 MeV in gamma rays
#:         per decay (the two 511 keV annihilation photons included) plus
#:         0.12 MeV kinetic energy of positrons, averaged over all decays.
E_GAMMA_NI_MEV = 1.75
E_GAMMA_CO_MEV = 3.61
E_POSITRON_CO_MEV = 0.12

#: The two strongest 56Co gamma-ray lines: energy (keV) and photons per decay
#: (NNDC / ENSDF, Junde et al. 2011). INTEGRAL detected both from SN 2014J.
CO56_LINES = {
    847: {"energy_kev": 846.770, "yield": 0.999399},
    1238: {"energy_kev": 1238.288, "yield": 0.6646},
}

# ---------------------------------------------------------------------------
# Cosmology for the ZTF sample (Arnett's rule for dozens of Type Ia)
# ---------------------------------------------------------------------------

#: Hubble constant on the Cepheid distance scale, so that ZTF distances match
#: the Cepheid distance to M101 used for SN 2011fe and SN 2023ixf (D1 found
#: 73.55 +/- 1.02 from Pantheon+; Riess et al. 2022 found 73.04 +/- 1.04).
H0 = 73.0
OMEGA_M = 0.3

# ---------------------------------------------------------------------------
# The four nearby supernovae
# ---------------------------------------------------------------------------

#: Distance moduli, dust and explosion epochs. ``t0_mjd`` is the explosion (or,
#: for the Type Ia, first-light) epoch that all phases are measured from.
SUPERNOVAE = {
    "SN2011fe": {
        "label": "SN 2011fe",
        "type": "Ia",
        "host": "M101",
        # Cepheid distance to M101 (Riess et al. 2016, ApJ 826, 56): 6.71 Mpc.
        # Shared with SN 2023ixf, so their luminosities scale together.
        "mu": 29.135,
        "mu_err": 0.045,
        "mu_ref": "Riess et al. 2016 (Cepheids)",
        # Milky Way: Schlafly & Finkbeiner (2011) recalibration of SFD.
        "ebv_mw": 0.0077,
        # Host: Na I D absorption is weak; Patat et al. (2013, A&A 549, A62)
        # put the host reddening at E(B-V) = 0.014 +/- 0.002.
        "ebv_host": 0.014,
        "rv_host": 3.1,
        # First light from the earliest PTF data (Nugent et al. 2011, Nature
        # 480, 344): 2011 Aug 23 16:29 UT +/- 20 min.
        "t0_mjd": 55796.687,
        "t0_ref": "Nugent et al. 2011 (first light)",
    },
    "SN2014J": {
        "label": "SN 2014J",
        "type": "Ia",
        "host": "M82",
        # TRGB distance used by Churazov et al. (2014): 3.53 +/- 0.26 Mpc
        # (Karachentsev & Kashibadze 2006; Dalcanton et al. 2009).
        "mu": 27.739,
        "mu_err": 0.16,
        "mu_ref": "3.53 +/- 0.26 Mpc (TRGB), as in Churazov et al. 2014",
        # SFD reddening toward M82 is contaminated by M82's own dust emission;
        # the foreground is E(B-V) = 0.05 (Dalcanton et al. 2009), as adopted
        # by Amanullah et al. (2014) and Foley et al. (2014).
        "ebv_mw": 0.05,
        # Host dust is heavy and unusual: E(B-V) = 1.37, R_V = 1.4 with a
        # Cardelli law fitted to HST UV-to-NIR colours (Amanullah et al. 2014,
        # ApJ 788, L21). This is the single largest uncertainty for 2014J.
        "ebv_host": 1.37,
        "rv_host": 1.4,
        # First light 2014 Jan 14.75 +/- 0.21 UT (Zheng et al. 2014, ApJ 783, L24).
        "t0_mjd": 56671.75,
        "t0_ref": "Zheng et al. 2014 (first light)",
    },
    "SN1987A": {
        "label": "SN 1987A",
        "type": "II-pec",
        "host": "LMC",
        # Detached eclipsing binaries (Pietrzynski et al. 2019, Nature 567,
        # 200): 49.59 kpc. The published UVOIR curves assumed 50 kpc.
        "mu": 18.477,
        "mu_err": 0.026,
        "mu_ref": "Pietrzynski et al. 2019 (LMC, eclipsing binaries)",
        # Total (Milky Way + LMC) reddening E(B-V) = 0.16 (Suntzeff 1997;
        # Walker & Suntzeff 1990 found 0.17 +/- 0.02). Not split by galaxy.
        "ebv_mw": 0.16,
        "ebv_host": 0.0,
        "rv_host": 3.1,
        # The Kamiokande-II / IMB neutrino burst: 1987 Feb 23 07:35:35 UT.
        "t0_mjd": 46849.316,
        "t0_ref": "neutrino burst (Hirata et al. 1987; Bionta et al. 1987)",
    },
    "SN2023ixf": {
        "label": "SN 2023ixf",
        "type": "II",
        "host": "M101",
        "mu": 29.135,
        "mu_err": 0.045,
        "mu_ref": "Riess et al. 2016 (Cepheids)",
        "ebv_mw": 0.0077,
        # Host: Na I D (Lundquist et al. 2023; Smith et al. 2023), as adopted
        # by Hsu et al. (2025) and Singh et al. (2024).
        "ebv_host": 0.031,
        "rv_host": 3.1,
        # Explosion epoch from the earliest amateur and survey photometry
        # (Li et al. 2024, Nature 627, 754), as adopted by Hsu et al. (2025).
        "t0_mjd": 60082.788,
        "t0_ref": "Li et al. 2024 (explosion)",
    },
}

#: Reddening and distance the published SN 1987A UVOIR curves were computed
#: with. We rescale them to our distance; their reddening is left as published.
SN1987A_PUBLISHED_DISTANCE_KPC = 50.0
