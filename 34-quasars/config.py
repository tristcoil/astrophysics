"""Paths, URLs, constants, cosmologies and every sample cut, with the reason for each.

Everything that a reader might want to change, or check against a paper, lives
here rather than inside the analysis code.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
FIGURES = HERE / "figures"
DERIVED = DATA / "derived"          # compact tables built from the raw downloads
RESULTS_JSON = DERIVED / "results.json"
RESULTS_PICKLE = DERIVED / "results.pkl"

USER_AGENT = "astrophysics-atlas-I3/1.0 (hobby research; python-requests)"
# VizieR answers scripted clients, but CDS's front end has been known to treat
# unusual agents differently; a browser-like agent is what was tested.
BROWSER_AGENT = "Mozilla/5.0 (X11; Linux x86_64) astrophysics-atlas-I3/1.0"

# ---------------------------------------------------------------------------
# Remote sources
# ---------------------------------------------------------------------------

VIZIER_ASU = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"

#: Wu & Shen (2022, ApJS 263, 42) spectral properties of all 750,414 DR16Q
#: quasars. The May 2024 version adds a host-galaxy fraction; everything else is
#: unchanged from the published catalogue. 1.39 GB gzipped FITS.
WUSHEN_URL = "http://quasar.astro.illinois.edu/paper_data/DR16Q/dr16q_prop_May01_2024.fits.gz"
WUSHEN_GZ = DATA / "wushen22" / "dr16q_prop_May01_2024.fits.gz"
WUSHEN_SIZE = 1_391_052_577

#: Columns kept from the Wu & Shen catalogue. Line arrays hold six numbers:
#: [peak wavelength, 50%-flux centroid wavelength, flux, log L, FWHM, rest EW].
WUSHEN_COLUMNS_EXT1 = [
    "SDSS_NAME", "PLATE", "MJD", "FIBERID", "RA", "DEC", "IF_BOSS_SDSS",
    "Z_DR16Q", "SOURCE_Z_DR16Q", "Z_FIT", "Z_SYS", "Z_SYS_ERR", "EBV", "SN_MEDIAN_ALL",
    "CONTI_FIT", "CONTI_STAT",
    "FEII_UV_EW", "FEII_OPT_EW",
    "LOGL1350", "LOGL1350_ERR", "LOGL1700", "LOGL2500", "LOGL3000", "LOGL3000_ERR",
    "LOGL5100", "LOGL5100_ERR",
    "HBETA_BR", "HBETA_BR_ERR", "OIII5007", "MGII_BR", "MGII_BR_ERR", "MGII", "CIV", "CIV_ERR",
    "HB_COMP_STAT", "MGII_COMP_STAT", "CIV_COMP_STAT",
    "LOGLBOL", "LOGLBOL_ERR", "LOGMBH_HB", "LOGMBH_HB_ERR", "LOGMBH_MGII", "LOGMBH_MGII_ERR",
    "LOGMBH_CIV", "LOGMBH_CIV_ERR", "LOGMBH", "LOGMBH_ERR", "LOGLEDD_RATIO", "LOGLEDD_RATIO_ERR",
    "ZSYS_LINES", "ZSYS_LINES_ERR", "FHOST_5100",
]
WUSHEN_COLUMNS_EXT2 = [
    "BAL_PROB", "BI_CIV", "PSFMAG", "PSFMAGERR", "EXTINCTION", "M_I", "FIRST_MATCHED", "FIRST_FLUX",
]
WUSHEN_TABLE = DERIVED / "wushen22_columns.parquet"

#: DR16Q (Lyke et al. 2020) from VizieR VII/289: just the columns the Wu & Shen
#: table lacks (redshift source and the target-selection bits that say which
#: survey programme picked each quasar). Fetched in RA strips.
DR16Q_SOURCE = "VII/289/dr16q"
DR16Q_COLUMNS = ["recno", "SDSS", "RAJ2000", "DEJ2000", "Plate", "MJD", "Fiber", "z", "r_z",
                 "QSO", "Btarg1", "eBOSST0", "eBOSST1", "eBOSST2", "Atarg1", "Atarg2", "NspI"]
DR16Q_DIR = DATA / "dr16q_vizier"
DR16Q_TABLE = DERIVED / "dr16q_columns.parquet"
DR16Q_N_EXPECTED = 750_414

#: Shen et al. (2011) DR7 quasar properties (J/ApJS/194/45): the uniform flag,
#: M_i(z=2) and dereddened i for the luminosity function, plus DR7 masses.
SHEN11_SOURCE = "J/ApJS/194/45/catalog"
SHEN11_COLUMNS = ["recno", "SDSS", "RAJ2000", "DEJ2000", "z", "Plate", "Fiber", "MJD", "Uni",
                  "iMAG", "logLbol", "BAL", "FIRST", "logBH", "e_logBH", "logEdd", "zHW",
                  "imag0", "gmag0", "Delg-i"]
SHEN11_N_EXPECTED = 105_783
#: The DR7 quasar catalogue (Schneider et al. 2010, VII/260) for the morphology
#: flag (extended sources have their own selection function) and FIRST fluxes.
DR7Q_SOURCE = "VII/260/dr7qso"
DR7Q_COLUMNS = ["recno", "SDSS", "z", "imag", "Au", "F20cm", "Mph", "utf", "iTmag", "Plate",
                "Fiber", "MJD-sp", "RAJ2000", "DEJ2000"]
DR7_DIR = DATA / "dr7"
DR7_TABLE = DERIVED / "dr7_uniform_input.parquet"

#: Richards et al. (2006, AJ 131, 2766): the selection function (Table 1),
#: K-corrections (Table 4), their 15,343-quasar statistical sample (Table 5) and
#: binned luminosity function (Table 6).
RICHARDS06_TABLES = {
    "selfunc": "J/AJ/131/2766/table1",
    "kcorr": "J/AJ/131/2766/table4",
    "sample": "J/AJ/131/2766/table5",
    "qlf": "J/AJ/131/2766/table6",
}
#: Shen & Kelly (2012, ApJ 746, 169): binned DR7 virial mass function and
#: luminosity function from the same uniform sample.
SHENKELLY12_TABLES = {
    "bhmf": "J/ApJ/746/169/table2",
    "qlf": "J/ApJ/746/169/table3",
}
#: z > 5.7 quasars with Mg II widths and 3000 A luminosities: Farina et al.
#: (2022) Table 5 compiles them from the literature onto one footing; Table 1 is
#: their own X-shooter sample. Shen et al. (2019) and Yang et al. (2021) are two
#: of the main sources, kept for cross-checks. Banados et al. (2016) list every
#: z > 5.6 quasar known in March 2016.
HIGHZ_TABLES = {
    "farina22": "J/ApJ/941/106/table5",
    "farina22_xshooter": "J/ApJ/941/106/table1",
    "shen19": "J/ApJ/873/35/table3",
    "yang21": "J/ApJ/923/262/table2",
    "banados16": "J/ApJS/227/11/highzqso",
}
LIT_DIR = DATA / "literature"

#: Vanden Berk et al. (2001): the median composite of 2,204 SDSS quasars and the
#: line list, both in the paper's arXiv source.
VB01_ARXIV_SRC = "https://arxiv.org/src/astro-ph/0105231"
VB01_TAR = LIT_DIR / "vandenberk2001_src.tar.gz"

#: SDSS filter responses (Doi et al. 2010 era curves served with DR7): point
#: source through 1.3 airmasses, column 2.
SDSS_FILTER_URL = "https://classic.sdss.org/dr7/instruments/imager/filters/{band}.dat"
FILTER_DIR = DATA / "filters"
BANDS = "ugriz"

#: A patch of the stellar locus from SkyServer, to show where quasars of each
#: redshift fall among ordinary stars in colour space.
SKYSERVER_URL = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch"
STARS_CSV = DATA / "stars" / "stellar_locus_patch.csv"
STAR_QUERY = """
SELECT p.objID, p.ra, p.dec,
       p.psfMag_u - p.extinction_u AS u, p.psfMag_g - p.extinction_g AS g,
       p.psfMag_r - p.extinction_r AS r, p.psfMag_i - p.extinction_i AS i,
       p.psfMag_z - p.extinction_z AS z
FROM Star AS p
WHERE p.ra BETWEEN 180 AND 190 AND p.dec BETWEEN 5 AND 15
  AND p.clean = 1
  AND p.psfMag_i - p.extinction_i BETWEEN 15.0 AND 19.1
  AND p.psfMagErr_u < 0.1 AND p.psfMagErr_g < 0.05 AND p.psfMagErr_r < 0.05
  AND p.psfMagErr_i < 0.05 AND p.psfMagErr_z < 0.05
"""

#: SDSS spectra: BOSS/eBOSS reductions of DR16 (the ones DR16Q and Wu & Shen
#: used), "lite" files without per-exposure frames.
SPEC_URL = "https://data.sdss.org/sas/dr16/eboss/spectro/redux/v5_13_0/spectra/lite/{plate:04d}/spec-{plate:04d}-{mjd:05d}-{fiber:04d}.fits"
SPEC_DIR = DATA / "spectra"
SPEC_LIST = DERIVED / "composite_sample.parquet"

# ---------------------------------------------------------------------------
# Physical constants and conventions
# ---------------------------------------------------------------------------

#: Eddington luminosity per solar mass for ionised hydrogen, 4 pi G m_p c / sigma_T,
#: in erg/s. (Wu & Shen's catalogue uses 1.3e38; the 3% difference is noted.)
L_EDD_PER_MSUN = 1.2572e38
L_EDD_PER_MSUN_WUSHEN = 1.3e38
#: Salpeter (e-folding) time at the Eddington limit for a radiative efficiency
#: eps: t = eps / (1 - eps) * sigma_T c / (4 pi G m_p). The last factor is 450 Myr.
T_EDD_MYR = 450.5
RADIATIVE_EFFICIENCY = 0.1
L_SUN = 3.828e33  # erg/s (IAU 2015 nominal)

#: Bolometric corrections (Richards et al. 2006b mean SED), as in Shen 2011 and
#: Wu & Shen 2022: L_bol = BC x lambda L_lambda.
BOLOMETRIC_CORRECTION = {5100: 9.26, 3000: 5.15, 1350: 3.81}

#: Single-epoch virial calibrations: log M = a + b log(lambda L / 1e44 erg/s)
#: + 2 log(FWHM / km/s). Coefficients from Shen et al. 2011, eqs. 5-8.
VIRIAL = {
    "HB_VP06": (0.910, 0.50, 5100),    # Vestergaard & Peterson 2006
    "MGII_S11": (0.740, 0.62, 3000),   # Shen et al. 2011 (their fiducial Mg II)
    "MGII_VO09": (0.860, 0.50, 3000),  # Vestergaard & Osmer 2009
    "CIV_VP06": (0.660, 0.53, 1350),   # Vestergaard & Peterson 2006
}
#: Coatman et al. (2017, MNRAS 465, 2120): C IV FWHM corrected for its outflow-
#: driven blueshift (their eq. 4), FWHM_corr = FWHM / (0.41 blueshift/(1000 km/s)
#: + 0.62), used in their eq. 6: log M = 6.71 + 2 log(FWHM_corr/1000 km/s)
#: + 0.53 log(lambda L_1350 / 1e44). The blueshift is c (1549.48 - lambda_half)
#: / 1549.48, lambda_half bisecting the line flux, relative to systemic.
COATMAN17 = dict(slope=0.41, intercept=0.62, a=0.71, b=0.53)
CIV_REST = 1549.48

# ---------------------------------------------------------------------------
# Cosmologies
# ---------------------------------------------------------------------------

#: Our working cosmology: project D3's own fit to the Planck temperature
#: spectrum (09-planck-cmb: H0 = 66.66 +/- 0.90, Omega_m = 0.325 +/- 0.013,
#: flat; age 13.825 Gyr). Radiation and a 0.06 eV neutrino are included as in
#: Planck's base model, so ages at z ~ 20-30 are right.
D3_H0 = 66.66
D3_OM0 = 0.325
D3_AGE_GYR = 13.825
#: Richards et al. 2006, Shen et al. 2011 and Shen & Kelly 2012 all used this:
REF_H0, REF_OM0 = 70.0, 0.30

# ---------------------------------------------------------------------------
# The luminosity function
# ---------------------------------------------------------------------------

#: Richards et al. 2006: 15,343 quasars over 1622 deg2.
R06_AREA_DEG2 = 1622.0
#: Shen & Kelly 2012 (from Shen et al. 2007b's appendix): the DR7 uniform sample
#: (UNIFORM_TARGET = 1) covers 6248 deg2.
DR7_UNIFORM_AREA_DEG2 = 6248.0
#: Flux limits of the uniform sample (dereddened i; Shen et al. 2011 note 10).
#: Richards 2006 split the two branches at z = 3.0, the DR7 flag at z = 2.9.
I_LIMIT_LOWZ = 19.1
I_LIMIT_HIGHZ = 20.2
Z_SPLIT_R06 = 3.0
Z_SPLIT_DR7 = 2.9
I_BRIGHT_LIMIT = 15.0   # targets brighter than this were not selected
#: Floor on the selection function: Richards et al. 2006, section 3.3.
SELFUNC_FLOOR = 1.0 / 3.0
#: Binning of Richards et al. 2006, section 6.1.
QLF_Z_EDGES = np.array([0.30, 0.68, 1.06, 1.44, 1.82, 2.20, 2.60, 3.00, 3.50, 4.00, 4.50, 5.00])
QLF_M_BRIGHT, QLF_M_STEP = -30.0, 0.3   # bins run from -22.5 upward (brighter) in 0.3 mag
QLF_M_FAINT = -22.5
#: The "quasar era" curve: quasars this luminous are inside the flux limits at
#: 0.8 < z < 4.8 (Richards et al. 2006, section 3.4).
LUMINOUS_M_CUT = -27.6
#: Conversions of Richards et al. 2006: M_1450 = M_i(z=2) + 1.486 (eq. 3) and
#: M_i(z=0) = M_i(z=2) + 0.596 (eq. 1).
M1450_FROM_MI2 = 1.486

# ---------------------------------------------------------------------------
# Black hole masses
# ---------------------------------------------------------------------------

#: A mass is used only if its measurement error is below this (Shen & Kelly 2012
#: rejected errors > 0.5 dex) and the line fit is not obviously bad.
MASS_ERR_MAX = 0.5
#: Line widths outside this range are fits to noise or narrow-line contamination.
FWHM_RANGE = (1000.0, 20000.0)
#: The Mg II / C IV overlap used to test C IV masses.
CIV_MGII_Z = (1.7, 2.2)

# ---------------------------------------------------------------------------
# The composite spectrum
# ---------------------------------------------------------------------------

#: Redshift range and binning for choosing spectra: equal steps in log(1+z) so
#: every part of the rest-frame spectrum gets similar coverage.
COMPOSITE_Z_RANGE = (0.06, 4.6)
COMPOSITE_N_ZBINS = 46
COMPOSITE_PER_BIN = 100
COMPOSITE_SEED = 20261004
COMPOSITE_MIN_SN = 4.0          # median S/N per pixel of the whole spectrum
COMPOSITE_MAX_BAL_PROB = 0.5    # broad absorption-line quasars would dent the stack
#: Rest-frame grid of the stack: 1 A pixels, as Vanden Berk et al. 2001 used.
COMPOSITE_GRID = np.arange(650.0, 9600.0 + 0.5, 1.0)
#: Observed-frame wavelength range trusted in BOSS spectra (the ends are noisy).
BOSS_TRUSTED = (3650.0, 10200.0)

# ---------------------------------------------------------------------------
# High-redshift quasars and seeds
# ---------------------------------------------------------------------------

HIGHZ_MIN_Z = 5.7
#: Seed formation epochs bracketing the literature: Population III remnants at
#: z ~ 20-30, direct-collapse black holes at z ~ 15-20.
SEED_Z_START = (30.0, 20.0, 15.0)
SEED_CLASSES = {
    "Population III remnant": (10.0, 1.0e3),
    "Dense star cluster": (1.0e3, 1.0e4),
    "Direct collapse": (1.0e4, 1.0e6),
}
