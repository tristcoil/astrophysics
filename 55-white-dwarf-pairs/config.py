"""Project configuration: paths, data sources, analysis choices and reference values.

Everything tunable lives here, together with the published numbers we grade our
own results against. Published *timing* results (Burdge et al. 2019's P, P-dot and
T0) are only ever used for comparison, never as inputs to our own timing fit.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
ATLAS_DIR = PROJECT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
ZTF_DIR = RAW_DIR / "ztf"
ZTF_CUTOUT_DIR = ZTF_DIR / "cutouts"
PTF_DIR = RAW_DIR / "ptf"
PTF_CUTOUT_DIR = PTF_DIR / "cutouts"
LIT_DIR = RAW_DIR / "literature"
DERIVED_DIR = DATA_DIR / "derived"
FIGURE_DIR = PROJECT_DIR / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, RAW_DIR, ZTF_DIR, ZTF_CUTOUT_DIR, PTF_DIR, PTF_CUTOUT_DIR, LIT_DIR,
           DERIVED_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# The target: ZTF J153932.16+502738.8 (Burdge et al. 2019, Nature 571, 528)
# ---------------------------------------------------------------------------

TARGET_NAME = "ZTF J1539+5027"
#: Gaia DR3 1402814555998561280: position at epoch J2016.0 and proper motion.
#: (G = 20.40, BP-RP = -0.26; its parallax, -0.49 +- 0.57 mas, is useless.)
#: The nearest other Gaia source is 20" away (G = 17.5).
TARGET_GAIA_DR3 = 1402814555998561280
TARGET_RA_2016 = 234.88397130556  # deg, ICRS
TARGET_DEC_2016 = 50.46075607758  # deg
TARGET_PMRA = -2.337  # mas/yr (includes cos dec)
TARGET_PMDEC = -2.860  # mas/yr

# ---------------------------------------------------------------------------
# ZTF: IRSA image archive (IBE). No account needed; every exposure listed by the
# metadata search is past its proprietary period (checked 2026-10-09).
# ---------------------------------------------------------------------------

IRSA = "https://irsa.ipac.caltech.edu"
ZTF_SEARCH_URL = IRSA + "/ibe/search/ztf/products/sci"
ZTF_DATA_URL = IRSA + "/ibe/data/ztf/products/sci"
#: Cutout size in pixels (1.012"/px). The difference-image PSF stamp is 25x25 px;
#: the extra margin gives empty sky to measure the noise of a PSF-weighted flux.
ZTF_CUTOUT_PIX = 61
ZTF_DOWNLOAD_WORKERS = 12

# ---------------------------------------------------------------------------
# PTF: public Level 1 science images, 2009-2012 (iPTF data at this position, used
# by Burdge et al. for their 2016 epoch, are not in the public archive).
# ---------------------------------------------------------------------------

PTF_SEARCH_URL = IRSA + "/ibe/search/ptf/images/level1"
PTF_DATA_URL = IRSA + "/ibe/data/ptf/images/level1"
#: Large enough to hold field stars for an empirical PSF and a zero point.
PTF_CUTOUT_ARCMIN = 8.0

# ---------------------------------------------------------------------------
# Physical constants (CODATA 2018 / IAU 2015 nominal values, SI units)
# ---------------------------------------------------------------------------

G = 6.67430e-11
C_LIGHT = 2.99792458e8
GM_SUN = 1.3271244e20  # IAU 2015 nominal solar mass parameter
M_SUN = GM_SUN / G
R_SUN = 6.957e8
PC = 3.0856775814913673e16
KPC = 1e3 * PC
YEAR = 365.25 * 86400.0
DAY = 86400.0
T_SUN = GM_SUN / C_LIGHT**3  # 4.925490947 microseconds

# ---------------------------------------------------------------------------
# Reference values we check against (never fitted inputs)
# ---------------------------------------------------------------------------

#: Burdge et al. 2019, Table 1 (ZTF J1539+5027).
BURDGE19 = {
    "P_s": (414.7915404, 0.0000029),
    "Pdot": (-2.373e-11, 0.005e-11),
    "T0_bjd_tdb": (2458305.6827886, 0.0000012),
    "M1": (0.610, 0.017, 0.022),  # value, +err, -err (Msun)
    "M2": (0.210, 0.014, 0.015),
    "R1_rsun": (1.562e-2, 0.038e-2),
    "R2_rsun": (3.140e-2, 0.053e-2),
    "a_rsun": (11.218e-2, 0.081e-2),
    "incl_deg": (84.15, 0.64, 0.57),
    "d_kpc": (2.34, 0.14),
    "snr_4yr": (143, 14, 13),
    #: P-dot from the four PTF/iPTF epochs alone, and from CHIMERA+KPED alone.
    "Pdot_ptf_only": (-2.378e-11, 0.049e-11),
    "Pdot_chimera_kped_only": (-2.487e-11, 0.19e-11),
    #: Expected tidal contribution to the decay if both stars are synchronised.
    "tide_fraction": 0.067,
}
