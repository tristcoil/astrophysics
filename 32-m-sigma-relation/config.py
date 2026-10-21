"""Project configuration: paths, constants, source URLs and our anchor points.

Everything we *assume* lives here, each value with its source. The published
fits we grade ourselves against (slopes, zero points, scatters) are kept
separately in ``literature.py`` and are never used as inputs to a fit.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
CACHE_DIR = DATA_DIR / "cache"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, RAW_DIR, CACHE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

#: Newton's G in the units of galaxy dynamics: pc (km/s)^2 / Msun.
#: From the IAU 2015 nominal solar GM (1.3271244e20 m^3 s^-2) and the parsec.
G_PC = 4.300917e-3
C_KMS = 299792.458  # speed of light, km/s
ARCSEC_PER_RAD = 206264.80625
PC_PER_MPC = 1e6

# cgs, for the feedback argument (King 2003) and Eddington luminosity
G_CGS = 6.67430e-8  # CODATA 2018
C_CGS = 2.99792458e10
MSUN_G = 1.98841e33  # IAU 2015 nominal GM_sun / CODATA G
SIGMA_T = 6.6524587321e-25  # Thomson cross-section, cm^2
M_P = 1.67262192369e-24  # proton mass, g
KM_CM = 1e5

#: Electron-scattering opacity of fully ionised gas, kappa = 0.2 (1 + X) cm^2/g.
#: X = 0.70 is the cosmic hydrogen mass fraction (kappa = 0.34); pure hydrogen
#: (X = 1, kappa = sigma_T / m_p = 0.40) is what King & Pounds (2015) use.
HYDROGEN_X = 0.70

#: Cosmic baryon fraction Omega_b / Omega_m from Planck 2018 (TT,TE,EE+lowE+
#: lensing: Omega_b h^2 = 0.02237, Omega_m h^2 = 0.1430) = 0.156. In King's
#: argument it is the gas fraction of the protogalaxy the black hole clears.
F_GAS_COSMIC = 0.02237 / 0.1430

#: Radiative efficiency of accretion: the fraction of accreted rest-mass energy
#: radiated. 0.1 is the standard value (Soltan 1982 argument; Yu & Tremaine 2002).
RADIATIVE_EFFICIENCY = 0.1

#: Absolute K-band magnitude of the Sun, as adopted by van den Bosch (2016) to
#: turn his 2MASS photometry into solar luminosities. We use the same value
#: for ATLAS3D's magnitudes so both are on one scale.
MK_SUN = 3.28

# ---------------------------------------------------------------------------
# Fitting conventions
# ---------------------------------------------------------------------------

#: All three compilations pivot the M-sigma relation at 200 km/s, where the
#: zero point and slope are least correlated: log M = a + b log(sigma/200).
SIGMA_PIVOT = 200.0
#: Bulge-mass and luminosity relations pivot at 1e11 (Msun or Lsun).
MASS_PIVOT = 1e11

#: Kormendy & Ho (2013) fit with every galaxy given the same error in log M,
#: 0.117 dex = the mean of their fitted sample (their Section 6.6, eqs. 2-3).
KH13_EQUAL_LOGM_ERR = 0.117

#: Angular resolution of the spectroscopy behind most dynamical masses:
#: HST/STIS slits are 0.1-0.2" wide and the diffraction limit at the 8500 A
#: calcium triplet is 0.09"; ground-based adaptive optics reaches 0.05-0.1".
RESOLUTION_ARCSEC = 0.1

#: Monte Carlo settings for the Bayesian fits (emcee).
MCMC_WALKERS = 32
MCMC_STEPS = 8000
MCMC_BURN = 2000
RNG_SEED = 20261004

# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------

#: van den Bosch (2016, ApJ 831, 134): 294 black hole masses with sigma_e,
#: 2MASS growth-curve photometry, method and references. The author's own
#: FITS file; the repository moved from remcovandenbosch/ to remco-space/
#: (GitHub redirects the old name).
VDB16_REPO = "https://raw.githubusercontent.com/remco-space/Black-Hole-Mass-compilation/master/"
VDB16_FILES = ["BHcompilation.fits", "README.md", "LICENSE"]
#: The same compilation as published in the journal (VizieR J/ApJ/831/134),
#: used only to check that the GitHub file is the published one.
VDB16_VIZIER = {
    "table2": "J/ApJ/831/134/table2",  # the 245 objects with complete data
    "table3": "J/ApJ/831/134/table3",  # the 33 without photometry
}

#: arXiv source bundles holding the tables of the two big reviews.
ARXIV = {
    # Kormendy & Ho (2013, ARA&A 51, 511). Table 2 (ellipticals) is TeX; Table 3
    # (spiral and S0 galaxies) exists only as a dvips PostScript file, whose
    # colours mark which galaxies the authors fitted.
    "1304.7762": {"label": "kh13", "keep": ["annrev4K-E.tex", "annrev4K-S.ps"]},
    # McConnell & Ma (2013, ApJ 764, 184): Table A1, 72 galaxies.
    "1211.2816": {"label": "mm13", "keep": ["arxiv_v3.tex"]},
}

#: The AGN Black Hole Mass Database (Bentz & Katz 2015, PASP 127, 67): black
#: hole masses from reverberation mapping. The site rescales every mass to a
#: chosen virial factor <f>; asking for f = 1 returns the raw "virial
#: products" c tau dV^2 / G, so we can calibrate f ourselves.
AGN_DB_URL = "http://www.astro.gsu.edu/AGNmass/"
AGN_DB_F = {"virial_products": "1", "grier13": "4.3"}

#: ATLAS3D (Cappellari et al. 2011a, MNRAS 413, 813; 2013a, MNRAS 432, 1709):
#: a volume-limited sample of all 260 early-type galaxies within 42 Mpc and
#: brighter than M_K = -21.5, with integral-field velocity dispersions.
ATLAS3D_URL = "https://groups.physics.ox.ac.uk/atlas3d/tables/"
ATLAS3D_FILES = {
    "sample": "Cappellari2011a_Atlas3D_Paper1_Table3.txt",  # distances, M_K
    "scaling": "Cappellari2013a_Atlas3D_Paper15_Table1.txt",  # sigma_e, R_e
}

# ---------------------------------------------------------------------------
# Anchor points: two black holes weighed without any galaxy-dynamics model
# ---------------------------------------------------------------------------

#: Sgr A*, from our own H1 project (17-galactic-centre): S2's orbit fitted to
#: both teams' positions and velocities. M = 3.957 +/- 0.042 (stat) +/- 0.02
#: (sys) million Msun at R0 = 7.93 kpc. If H1's results file is present we
#: read it and check this number against it.
MILKY_WAY = {
    "label": "Milky Way (Sgr A*, H1)",
    "M": 3.957e6,
    "M_err": (0.042e6**2 + 0.02e6**2) ** 0.5,
    "M_note": "S2 orbit, H1 (17-galactic-centre); GRAVITY 2022 gives 4.297e6 at 8.277 kpc",
    #: The bulge velocity dispersion is ill-defined for our own Galaxy, which
    #: we see from inside. 105 +/- 20 km/s is Kormendy & Ho's (2013) sigma_e;
    #: McConnell & Ma use 103 +/- 20 (Tremaine et al. 2002), van den Bosch 101.
    "sigma": 105.0,
    "sigma_err": 20.0,
    "h1_results": PROJECT_DIR.parent / "17-galactic-centre" / "data" / "results.json",
}

#: M87*, weighed from the size of its shadow by the Event Horizon Telescope
#: (EHT Collaboration 2019, ApJL 875, L6): M = 6.5 +/- 0.2 (stat) +/- 0.7 (sys)
#: billion Msun at D = 16.8 Mpc. sigma_e = 324 (+28 / -12) km/s from Kormendy &
#: Ho (2013). (Project I5 will reconstruct the shadow itself.)
M87 = {
    "label": "M87* (EHT shadow)",
    "M": 6.5e9,
    "M_err": (0.2e9**2 + 0.7e9**2) ** 0.5,
    "D_mpc": 16.8,
    "sigma": 324.0,
    "sigma_err_hi": 28.0,
    "sigma_err_lo": 12.0,
}
