"""Project configuration: paths, the data products, analysis settings, and the
published values our results are graded against.

Three kinds of data are used, all free and public:

1. **COBE/FIRAS** -- the spectrum of the CMB itself: how bright it is at each
   frequency. This is what says the background is a 2.7 K blackbody. (2 KB)
2. **Planck 2018 power spectra** -- the Planck team's own measurement of how
   much the temperature fluctuates on each angular scale, plus their best-fit
   model. (~0.4 MB)
3. **Planck 2018 sky maps** -- the SMICA CMB temperature map, as two independent
   half-mission maps, with the masks needed to analyse them. (~1.6 GB)

Group 3 is optional: every spectrum-only result runs without it.

Every number that comes from a paper rather than from data we download is kept
at the bottom of this file, with its exact source, so it is obvious which values
are ours and which are the literature's.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
PLANCK_DIR = DATA_DIR / "planck"
FIRAS_DIR = DATA_DIR / "firas"
#: Expensive intermediate products (apodised mask, map spectra, fit results).
#: Everything here is regenerable from the raw downloads.
DERIVED_DIR = DATA_DIR / "derived"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, PLANCK_DIR, FIRAS_DIR, DERIVED_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------
# Data products
# --------------------------------------------------------------------------

#: The ESA Planck Legacy Archive is the primary source. The NASA/IPAC IRSA
#: mirror serves identical files; it was ~10x slower from Europe in testing but,
#: unlike the PLA, supports HTTP range requests -- so it is the one that can
#: resume an interrupted download, and the only one that can serve a slice of a
#: file (which is how the beam is fetched without the 2 GB map around it).
PLA_COSMO = "https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID={}"
PLA_MAP = "https://pla.esac.esa.int/pla/aio/product-action?MAP.MAP_ID={}"
IRSA_COSMO = "https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/cosmoparams/{}"
IRSA_CMB = ("https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/"
            "component-maps/cmb/{}")
IRSA_MASKS = "https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/masks/{}"
LAMBDA_FIRAS = "https://lambda.gsfc.nasa.gov/data/cobe/firas/monopole_spec/{}"


@dataclass(frozen=True)
class Product:
    """One file to download and cache."""

    key: str
    filename: str
    #: "firas", "spectra" (small) or "maps" (large, optional).
    group: str
    #: Expected size in bytes; a download is only accepted if it matches.
    size_bytes: int
    #: Candidate URLs, tried in order.
    urls: tuple[str, ...]
    what: str

    @property
    def path(self) -> Path:
        return (FIRAS_DIR if self.group == "firas" else PLANCK_DIR) / self.filename


def _planck_spectrum(key, filename, size, what):
    return Product(key, filename, "spectra", size,
                   (PLA_COSMO.format(filename), IRSA_COSMO.format(filename)), what)


def _planck_map(key, filename, size, irsa, what):
    return Product(key, filename, "maps", size,
                   (PLA_MAP.format(filename), irsa.format(filename)), what)


PRODUCTS: list[Product] = [
    Product(
        "firas", "firas_monopole_spec_v1.txt", "firas", 2388,
        (LAMBDA_FIRAS.format("firas_monopole_spec_v1.txt"),),
        "COBE/FIRAS CMB monopole spectrum (Fixsen et al. 1996, Table 4)",
    ),
    _planck_spectrum(
        "tt_full", "COM_PowerSpect_CMB-TT-full_R3.01.txt", 170547,
        "Planck 2018 TT power spectrum, one point per multipole, l = 2-2508",
    ),
    _planck_spectrum(
        "tt_binned", "COM_PowerSpect_CMB-TT-binned_R3.01.txt", 7143,
        "Planck 2018 TT power spectrum in bins of 30 multipoles, l >= 30",
    ),
    _planck_spectrum(
        "theory", "COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt",
        205647, "Best-fit LCDM model spectra from the Planck 2018 baseline fit",
    ),
    _planck_spectrum(
        "bestfit", "COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum_R3.01.txt",
        8842, "Parameters of that best-fit model, including derived quantities",
    ),
    _planck_map(
        "hm1", "COM_CMB_IQU-smica_2048_R3.00_hm1.fits", 603990720, IRSA_CMB,
        "SMICA CMB map from the first half of the mission",
    ),
    _planck_map(
        "hm2", "COM_CMB_IQU-smica_2048_R3.00_hm2.fits", 603990720, IRSA_CMB,
        "SMICA CMB map from the second half of the mission",
    ),
    _planck_map(
        "mask_common", "COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits", 201335040, IRSA_MASKS,
        "Common temperature confidence mask (Galaxy + point sources), f_sky = 77.9%",
    ),
    _planck_map(
        "mask_hm", "COM_Mask_CMB-HM-Misspix-Mask-Int_2048_R3.00.fits", 201335040, IRSA_MASKS,
        "Pixels unobserved in either half-mission, f_sky = 96.0%",
    ),
]

#: The SMICA beam is not a separate file: it is the last HDU of the 2 GB
#: full-mission map. We fetch only that HDU with an HTTP range request to IRSA.
BEAM_SOURCE_URL = IRSA_CMB.format("COM_CMB_IQU-smica_2048_R3.00_full.fits")
BEAM_SOURCE_SIZE = 2013312960
BEAM_PATH = PLANCK_DIR / "smica_beam_R3.00.fits"


def product(key: str) -> Product:
    for p in PRODUCTS:
        if p.key == key:
            return p
    raise KeyError(f"unknown product {key!r}; known: {[p.key for p in PRODUCTS]}")


# --------------------------------------------------------------------------
# Analysis settings
# --------------------------------------------------------------------------

#: HEALPix resolution of the Planck maps: 12 * 2048**2 = 50,331,648 pixels of
#: 1.7 arcmin.
NSIDE = 2048
#: Highest multipole in our own harmonic transforms. Planck's spectrum stops at
#: 2508; going a little beyond lets us watch the noise take over.
LMAX = 3000

#: Mask apodisation. A mask with hard edges couples multipoles strongly (its own
#: harmonic transform falls off slowly), so we taper it smoothly to zero. The
#: taper is a Gaussian of this FWHM, applied after eroding the mask by about the
#: same amount so the taper lies inside the region we trust. 30 arcmin is the
#: compromise: long enough to tame the coupling, short enough that tapering
#: thousands of small point-source holes does not eat the sky.
APODIZATION_FWHM_DEG = 0.5

#: Resolution used for the whole-sky Mollweide images: 3.1 million pixels of
#: 6.9 arcmin, about the size of one pixel of the printed figure. Averaging the
#: 1.7-arcmin map down to it loses nothing visible but removes the noise speckle
#: that nearest-pixel sampling of the full-resolution map would produce.
DISPLAY_NSIDE = 512
#: Colour-scale limit for temperature maps, in microkelvin.
MAP_LIMIT_UK = 300.0

#: Rough peak positions, used only to centre the fitting windows -- each window
#: is wide enough that the fitted position does not depend on the guess.
PEAK_GUESSES = (220, 537, 811, 1121, 1444, 1776, 2081, 2395)
#: Half-width of each fitting window, in multipoles. A peak's position is pinned
#: by its steep flanks more than its flat top, so wider windows give smaller
#: errors -- until the Gaussian shape stops describing the peak. At +/-100 the
#: method is biased by at most ~1 multipole on Planck's noise-free model (tested
#: against the model's exact maxima), while its statistical error is ~2.
PEAK_HALF_WINDOWS = (100, 100, 100, 100, 80, 90, 110, 110)
#: From this peak on, the damping tail is divided out before fitting (see
#: ``spectra.damping_envelope``); the envelope's polynomial order is 5, with 4
#: and 6 bracketing the systematic error.
FIRST_ENVELOPE_PEAK = 5
ENVELOPE_ORDERS = (4, 5, 6)
#: Monte Carlo draws used to put error bars on peak positions and heights.
PEAK_MC_DRAWS = 300
RANDOM_SEED = 20180717  # Planck 2018 parameter paper arXiv date

#: Parameter fit: which parameters float, and where they start. Everything not
#: listed is held at the Planck 2018 baseline convention (flat, one massive
#: neutrino of 0.06 eV, N_eff = 3.046, helium from BBN).
FIT_PARAMETERS = ("ombh2", "omch2", "theta100", "logA", "ns")
FIT_START = {"ombh2": 0.0224, "omch2": 0.120, "theta100": 1.0410, "logA": 3.04, "ns": 0.965}
#: High-l temperature data constrain only the combination A_s exp(-2 tau), not
#: tau itself -- that needs large-scale polarisation, which we do not use. So
#: tau is held at Planck's value from their TT+lowE fit and we quote the
#: combination the data actually measure.
FIT_TAU = 0.0522
#: Multipole range fitted (Planck's binned high-l spectrum).
FIT_LMIN = 30
FIT_LMAX = 2508
#: Highest multipole CAMB computes. Lensing needs a margin above what we fit.
CAMB_LMAX = 2800


# --------------------------------------------------------------------------
# Published values we compare against -- NOT inputs to any measurement
# --------------------------------------------------------------------------

#: Present-day CMB temperature, Fixsen (2009, ApJ 707, 916): FIRAS recalibrated
#: and combined with other measurements.
T_CMB_K = 2.72548
T_CMB_ERR_K = 0.00057
#: FIRAS alone, Fixsen & Mather (2002, ApJ 581, 817), 1 sigma, calibration-dominated.
T_FIRAS_2002_K = 2.725
T_FIRAS_2002_ERR_K = 0.001
#: FIRAS 95% upper limits on spectral distortions, Fixsen et al. (1996, ApJ 473, 576).
FIRAS_Y_LIMIT_95 = 15e-6
FIRAS_MU_LIMIT_95 = 9e-5

#: Planck 2015 results I (A&A 594, A1, 2016), Table 7: TT peak positions and
#: amplitudes, from Gaussian fits to the 2015 CMB-only bandpowers.
#: (number, l, l_err, D_l [uK^2], D_l_err)
PUBLISHED_TT_PEAKS = (
    (1, 220.0, 0.5, 5717.0, 35.0),
    (2, 537.5, 0.7, 2582.0, 11.0),
    (3, 810.8, 0.7, 2523.0, 10.0),
    (4, 1120.9, 1.0, 1237.0, 4.0),
    (5, 1444.2, 1.1, 797.1, 3.1),
    (6, 1776.0, 5.0, 377.4, 2.9),
    (7, 2081.0, 25.0, 214.0, 4.0),
    (8, 2395.0, 24.0, 105.0, 4.0),
)

#: Planck 2018 parameter tables (PLA COM_CosmoParams_parameter-tables_R3.01,
#: baseline_params_table_2018_68pc.pdf), table 2.1 "base plikHM TT lowl lowE":
#: temperature plus low-l polarisation. The fair comparison for a fit to
#: temperature data. Mean and 68% interval.
PLANCK_TT_LOWE = {
    "ombh2": (0.02212, 0.00022),
    "omch2": (0.1206, 0.0021),
    "theta100": (1.04077, 0.00047),
    "tau": (0.0522, 0.0080),
    "logA": (3.040, 0.016),
    "ns": (0.9626, 0.0057),
    "H0": (66.88, 0.92),
    "omegal": (0.679, 0.013),
    "omegam": (0.321, 0.013),
    "omegamh2": (0.1434, 0.0020),
    "Ase2tau": (1.884, 0.014),  # 10^9 A_s exp(-2 tau)
    "sigma8": (0.8118, 0.0089),
    "age": (13.830, 0.037),
    "zstar": (1090.30, 0.41),
    "rstar": (144.46, 0.48),
    "thetastar100": (1.04097, 0.00046),
    "DMstar_gpc": (13.878, 0.044),
    "zeq": (3411.0, 48.0),
    "rdrag": (147.21, 0.48),
}

#: Same source, table 2.97 "base plikHM TTTEEE lowl lowE lensing": Planck's
#: headline result using everything. Context only.
PLANCK_BASELINE = {
    "ombh2": (0.02237, 0.00015),
    "omch2": (0.1200, 0.0012),
    "theta100": (1.04092, 0.00031),
    "tau": (0.0544, 0.0073),
    "logA": (3.044, 0.014),
    "ns": (0.9649, 0.0042),
    "H0": (67.36, 0.54),
    "omegam": (0.3153, 0.0073),
    "age": (13.797, 0.023),
}

#: The local distance-ladder value that disagrees with Planck: SH0ES,
#: Riess et al. (2022, ApJL 934, L7).
H0_SH0ES = (73.04, 1.04)
#: The same quantity measured in this repository by project D1
#: (05-supernova-cosmology/README.md, as of 2026-09-29): Pantheon+ supernovae,
#: Cepheid-calibrated, full covariance, flat LCDM. Copied rather than read from
#: D1's cache so this project runs on its own.
H0_D1 = (73.55, 1.02)
