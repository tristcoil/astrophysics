"""Project configuration: paths, data products, cosmology, analysis choices and the
published numbers our results are graded against.

The integrated Sachs-Wolfe (ISW) signal is a few microkelvin hidden under CMB
fluctuations of +-100 microkelvin, so the analysis has many knobs (filter sizes,
which structures to stack, multipole ranges). Turning knobs after looking at the
answer is the classic way to fool yourself. Every choice that touches the real
CMB map is therefore written down in ``PREREGISTRATION`` below, and frozen to
``data/preregistration.json`` (with a hash and a timestamp) before the first
real-sky stack of each catalogue is computed. Exploratory scans over other
choices are allowed afterwards, but they are labelled as exploratory and their
"look-elsewhere" cost is measured on simulations.

Every number that comes from a paper rather than from our own data is kept in
the ``LITERATURE`` block at the bottom, with its exact source.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
REPO_DIR = PROJECT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
PLANCK_DIR = RAW_DIR / "planck"
CATALOG_DIR = RAW_DIR / "catalogs"
BOSS_DIR = RAW_DIR / "boss"
#: Expensive intermediate products; all regenerable from the raw downloads.
DERIVED_DIR = DATA_DIR / "derived"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
PREREG_PATH = DATA_DIR / "preregistration.json"

#: D3 (09-planck-cmb) already downloaded the Planck maps. If its cache is
#: present we hard-link the files instead of downloading 1.6 GB again.
D3_PLANCK_DIR = REPO_DIR / "09-planck-cmb" / "data" / "planck"
D3_FIT_PATH = REPO_DIR / "09-planck-cmb" / "data" / "derived" / "tt_fit.json"
#: D2 (07-cosmic-web) cached the SDSS DR17 Main Galaxy Sample.
D2_MGS_PATH = REPO_DIR / "07-cosmic-web" / "data" / "mgs_dr17.npz"

for _d in (DATA_DIR, RAW_DIR, PLANCK_DIR, CATALOG_DIR, BOSS_DIR, DERIVED_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------
# Data products
# --------------------------------------------------------------------------

PLA_MAP = "https://pla.esac.esa.int/pla/aio/product-action?MAP.MAP_ID={}"
PLA_COSMO = "https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID={}"
IRSA_CMB = ("https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/"
            "component-maps/cmb/{}")
IRSA_MASKS = "https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/masks/{}"
IRSA_COSMO = "https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/cosmoparams/{}"
CDS_FTP = "https://cdsarc.cds.unistra.fr/ftp/{}"
SAS_LSS = "https://data.sdss.org/sas/dr12/boss/lss/{}"
WAYBACK_GRANETT = ("https://web.archive.org/web/20100613125329id_/"
                   "http://ifa.hawaii.edu/cosmowave/supervoids/{}")


@dataclass(frozen=True)
class Product:
    """One file to download and cache."""

    key: str
    filename: str
    #: Sub-directory of data/raw.
    subdir: str
    #: Expected size in bytes (None = unknown; accept any non-empty file).
    size_bytes: int | None
    urls: tuple[str, ...]
    what: str
    #: User-Agent to send. CDS refuses browser-like agents on its FTP mirror;
    #: the Wayback Machine and SAS are happy with anything.
    agent: str = "astrophysics-atlas-D10/1.0 (educational; python-requests)"

    @property
    def path(self) -> Path:
        return RAW_DIR / self.subdir / self.filename


def _cds(key, catalog, filename, size, what):
    return Product(key, f"{catalog.replace('/', '_')}_{filename}", "catalogs", size,
                   (CDS_FTP.format(f"{catalog}/{filename}"),), what, agent="curl/8.5")


PRODUCTS: list[Product] = [
    # --- Planck PR3 (identical to D3's files) ---------------------------------
    Product("hm1", "COM_CMB_IQU-smica_2048_R3.00_hm1.fits", "planck", 603990720,
            (PLA_MAP.format("COM_CMB_IQU-smica_2048_R3.00_hm1.fits"),
             IRSA_CMB.format("COM_CMB_IQU-smica_2048_R3.00_hm1.fits")),
            "SMICA CMB map, first half of the mission"),
    Product("hm2", "COM_CMB_IQU-smica_2048_R3.00_hm2.fits", "planck", 603990720,
            (PLA_MAP.format("COM_CMB_IQU-smica_2048_R3.00_hm2.fits"),
             IRSA_CMB.format("COM_CMB_IQU-smica_2048_R3.00_hm2.fits")),
            "SMICA CMB map, second half of the mission"),
    Product("mask_common", "COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits", "planck", 201335040,
            (PLA_MAP.format("COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits"),
             IRSA_MASKS.format("COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits")),
            "Common temperature confidence mask (Galaxy + point sources), f_sky = 77.9%"),
    Product("mask_hm", "COM_Mask_CMB-HM-Misspix-Mask-Int_2048_R3.00.fits", "planck", 201335040,
            (PLA_MAP.format("COM_Mask_CMB-HM-Misspix-Mask-Int_2048_R3.00.fits"),
             IRSA_MASKS.format("COM_Mask_CMB-HM-Misspix-Mask-Int_2048_R3.00.fits")),
            "Pixels unobserved in either half-mission"),
    Product("theory",
            "COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt",
            "planck", 205647,
            (PLA_COSMO.format("COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-"
                              "theory_R3.01.txt"),
             IRSA_COSMO.format("COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-"
                               "theory_R3.01.txt")),
            "Planck 2018 best-fit LCDM spectra: the input for simulated CMB skies"),
    # --- WMAP: an independent CMB map of the same sky ---------------------------
    Product("wmap_ilc", "wmap_ilc_9yr_v5.fits", "planck", 25174080,
            ("https://lambda.gsfc.nasa.gov/data/map/dr5/dfp/ilc/wmap_ilc_9yr_v5.fits",),
            "WMAP 9-year Internal Linear Combination map (1 deg FWHM, Nside 512)"),
    # --- Superstructure catalogues ------------------------------------------------
    Product("granett_voids", "granett2008_voidcat.txt", "catalogs", 3679,
            (WAYBACK_GRANETT.format("voidcat.txt"),),
            "Granett, Neyrinck & Szapudi (2008): 50 supervoids (Internet Archive copy)"),
    Product("granett_clusters", "granett2008_clustercat.txt", "catalogs", 3679,
            (WAYBACK_GRANETT.format("clustercat.txt"),),
            "Granett, Neyrinck & Szapudi (2008): 50 superclusters (Internet Archive copy)"),
    _cds("mao_readme", "J/ApJ/835/161", "ReadMe", None, "Mao et al. (2017) ReadMe"),
    _cds("mao_voids", "J/ApJ/835/161", "table1.dat", None,
         "Mao et al. (2017): 1,228 BOSS DR12 voids after quality cuts"),
    _cds("nh14_readme", "J/MNRAS/440/1248", "ReadMe", None, "Nadathur & Hotchkiss (2014) ReadMe"),
    *[_cds(f"nh14_{f}", "J/MNRAS/440/1248", f"{f}.dat", None,
           f"Nadathur & Hotchkiss (2014) DR7 structures: {f}")
      for f in ("lrgdt1v", "lrgdt1c", "lrgbt1v", "lrgbt1c", "lrgdt2v", "lrgbt2v",
                "dim1t1v", "dim2t1v", "bri1t1v", "bri2t1v",
                "dim1t1c", "dim2t1c", "bri1t1c", "bri2t1c")],
    # --- Galaxies and radio sources ---------------------------------------------
    Product("boss_north", "galaxy_DR12v5_CMASSLOWZTOT_North.fits.gz", "boss", 216513758,
            (SAS_LSS.format("galaxy_DR12v5_CMASSLOWZTOT_North.fits.gz"),),
            "BOSS DR12 large-scale-structure galaxies, northern cap"),
    Product("boss_south", "galaxy_DR12v5_CMASSLOWZTOT_South.fits.gz", "boss", 83199463,
            (SAS_LSS.format("galaxy_DR12v5_CMASSLOWZTOT_South.fits.gz"),),
            "BOSS DR12 large-scale-structure galaxies, southern cap"),
    Product("boss_mask_north", "mask_DR12v5_CMASSLOWZTOT_North.fits.gz", "boss", 2153013,
            (SAS_LSS.format("mask_DR12v5_CMASSLOWZTOT_North.fits.gz"),),
            "BOSS DR12 survey footprint (mangle polygons with completeness), north"),
    Product("boss_mask_south", "mask_DR12v5_CMASSLOWZTOT_South.fits.gz", "boss", 2102564,
            (SAS_LSS.format("mask_DR12v5_CMASSLOWZTOT_South.fits.gz"),),
            "BOSS DR12 survey footprint (mangle polygons with completeness), south"),
    _cds("nvss", "VIII/65", "nvss.dat.gz", 77687389,
         "NRAO VLA Sky Survey catalogue (Condon et al. 1998), 1,773,484 sources"),
]


def product(key: str) -> Product:
    for p in PRODUCTS:
        if p.key == key:
            return p
    raise KeyError(f"unknown product {key!r}")


# --------------------------------------------------------------------------
# Cosmology: D3's own fit to the Planck temperature spectrum
# --------------------------------------------------------------------------

#: Project D3 fitted a five-parameter LCDM model to Planck's TT band powers with
#: CAMB (09-planck-cmb/data/derived/tt_fit.json). The ISW prediction uses those
#: parameters, so "LCDM" here means "the universe D3 measured". tau is fixed to
#: Planck's TT+lowE value as in D3; one massive neutrino of 0.06 eV.
COSMO = {
    "H0": 66.65839714063198,
    "ombh2": 0.022210836834538488,
    "omch2": 0.12142427668531083,
    "ns": 0.9628705248868852,
    "logA": 3.0438697272087976,  # ln(10^10 A_s)
    "tau": 0.0522,
    "mnu": 0.06,
    "TCMB": 2.7255,
}
#: Derived in D3 (CAMB); checked against our own background code in tests.
COSMO_OMEGA_M = 0.32471058009602766

# --------------------------------------------------------------------------
# CMB maps
# --------------------------------------------------------------------------

#: Working resolution for everything on the CMB: 3.1 million pixels of 6.9'.
#: The ISW lives on degree scales; the smallest filter used is ~0.5 deg.
NSIDE = 512
#: Simulated CMB skies are drawn up to this multipole (3 Nside - 1).
LMAX_SIM = 1535
#: SMICA's effective beam (5' Gaussian) is kept in the sims, as is the pixel
#: window of the working resolution. Neither matters on degree scales.
SMICA_FWHM_ARCMIN = 5.0
#: A degraded pixel is kept only if all 16 of its full-resolution children are
#: unmasked; a partly masked pixel could carry Galactic or point-source residue.
MASK_KEEP_FRACTION = 1.0
#: Number of Gaussian CMB simulations used for every error bar and p-value.
N_SIMS = 1000
RANDOM_SEED = 19670101  # Sachs & Wolfe (1967), ApJ 147, 73: January 1967

# --------------------------------------------------------------------------
# Pre-registered analysis choices
# --------------------------------------------------------------------------

#: Every choice that determines a headline number. Frozen (hash + UTC time) to
#: data/preregistration.json before the corresponding real-sky stack is first
#: computed; the stacking code refuses to unblind a catalogue whose frozen entry
#: is missing or differs from this dictionary.
PREREGISTRATION = {
    "cmb": {
        "map": "Planck PR3 SMICA, mean of the two half-mission maps",
        "mask": "Planck common temperature mask x half-mission missing-pixel mask; a degraded "
                "pixel is kept only if all its children are unmasked",
        "nside": NSIDE,
        "monopole_dipole": "removed, fitted outside the mask",
        "errors": f"{N_SIMS} Gaussian simulations of the Planck 2018 best-fit LCDM spectrum, "
                  "stacked at the real positions; random positions in the survey footprint as a "
                  "cross-check",
        "filter": "compensated top hat: mean T inside radius R minus mean T in the ring "
                  "R..sqrt(2) R; masked pixels excluded; a structure is dropped if less than "
                  "half of either region survives the mask",
    },
    "granett": {
        "purpose": "replication of Granett, Neyrinck & Szapudi (2008) on Planck",
        "selection": "all 50 supervoids and 50 superclusters of the published lists",
        "radius_deg": 4.0,
        "statistic": "mean filtered T of voids, of clusters, and combined = "
                     "(mean_clusters - mean_voids) / 2",
        "note": "radius and N were chosen by Granett et al. after seeing WMAP; Planck sees the "
                "same CMB, so this is a replication, not an independent test",
    },
    "nh14": {
        "purpose": "replication of Hotchkiss et al. (2015) and an independent ISW test",
        "selection": "Type1 voids and Type1 superclusters of the lrgdim and lrgbright samples, "
                     "as published (no further cuts)",
        "filter_scale": "R = 0.6 x theta_eff (the optimum Hotchkiss et al. found in the Jubilee "
                        "simulation, also Cai et al. 2014); not tuned here",
        "statistic": "mean filtered T per sample and type; A_ISW from the matched-filter fit "
                     "of all four sets jointly to the BOSS-predicted LCDM signal",
        "matched_filter": "A = sum(p_i T_i / s_i^2) / sum(p_i^2 / s_i^2), p_i = filtered BOSS "
                          "prediction, s_i = per-structure scatter in the simulations; error = "
                          "scatter of the same estimator over the simulations",
    },
    "mao": {
        "purpose": "the largest void sample: BOSS DR12 (LOWZ + CMASS)",
        "selection": "all 1,228 voids of the quality-cut catalogue, as published",
        "filter_scale": "R = 0.6 x theta_v with theta_v = R_eff / comoving distance",
        "statistic": "mean filtered T, and A_ISW from the matched filter of the nh14 entry "
                     "with p_i from our BOSS reconstruction (a full covariance cannot be "
                     "estimated from 1,000 simulations for 1,228 correlated voids)",
    },
    "template": {
        "purpose": "the whole BOSS ISW map as a template (not just its voids)",
        "statistic": "A_ISW = <F t, d> / <F t, t> with t the predicted ISW map, d the CMB, "
                     "F the inverse of Planck's C_l^TT restricted to l = 10-100, Nside 64; "
                     "error from the simulations",
        "reconstruction": "Wiener filter of the BOSS DR12 CMASSLOWZTOT galaxies (0.15 < z < "
                          "0.75) with D3's linear P(k), Poisson noise per cell, LOWZ/CMASS "
                          "bias from their power spectra (k = 0.02-0.08 h/Mpc), Kaiser factor "
                          "b + f/3; potential smoothed by a 10 Mpc/h Gaussian",
        "systematics": "smoothing 5 and 20 Mpc/h, bias -/+10%, no Kaiser term, no Wiener filter",
    },
    "crosscorr": {
        "tracers": "NVSS (S >= 2.5 mJy, Dec >= -40 deg, |b| >= 10 deg, declination-band "
                   "density correction) and BOSS in two redshift slices",
        "nside": 64,
        "ell_range": [10, 100],
        "bin_width": 10,
        "statistic": "A_ISW: fit of C_l(galaxy x CMB) to CAMB's prediction with the tracer bias "
                     "fitted to its own auto-spectrum (l = 10-100); covariance from simulations, "
                     "rotated maps as a cross-check",
        "nvss_dndz": "Planck 2015 XXI: dN/dz ~ (z/0.33)^0.37 exp(-0.37 z/0.33), bias shape "
                     "0.90 [1 + 0.54 (1+z)^2] with its amplitude refitted",
        "boss_dndz": "the spectroscopic redshift histogram of each slice",
    },
    "coldspot": {
        "filter": "spherical Mexican-hat wavelet, R = 5 deg (Cruz et al. 2005)",
        "note": "the Cold Spot was found by searching; any p-value for it is a posteriori",
    },
}


def prereg_hash(entry: dict) -> str:
    """Stable SHA-256 of one pre-registration entry."""
    blob = json.dumps(entry, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()


# --------------------------------------------------------------------------
# Stacking
# --------------------------------------------------------------------------

#: Exploratory grid for the Granett catalogue (Granett et al.'s Table 1 spans
#: 3-5 deg; we scan wider). Fixed before looking, like everything else.
SCAN_RADII_DEG = tuple(round(1.0 + 0.5 * i, 1) for i in range(19))  # 1.0 ... 10.0
SCAN_N = (10, 20, 30, 40, 50)
#: Exploratory rescaling factors alpha in R = alpha * theta for the void
#: catalogues (0.6 is the pre-registered one).
SCAN_ALPHA = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5, 2.0)
#: A filter region must keep at least this fraction of its pixels.
MIN_UNMASKED = 0.5
#: Random positions drawn in the survey footprint for the cross-check errors.
N_RANDOM_POSITIONS = 20000
#: Side and pixel size of the stacked images (gnomonic projection).
STACK_IMAGE_DEG = 24.0
STACK_IMAGE_PIX = 97

# --------------------------------------------------------------------------
# BOSS reconstruction
# --------------------------------------------------------------------------

BOSS_ZMIN = 0.15
BOSS_ZMAX = 0.75
#: Split between LOWZ and CMASS for the bias measurement.
BOSS_ZSPLIT = 0.43
#: Resolution at which the mangle footprint is rasterised.
BOSS_MASK_NSIDE = 1024
#: Comoving grid: cells of GRID_CELL h^-1 Mpc, box padded so periodic images
#: of the survey sit far away.
GRID_CELL = 15.625  # 6000 / 384
GRID_SHAPE = (384, 384, 256)
#: Box origin offset along z (Mpc/h): the survey lies mostly north of the
#: equator, so the box is shifted up to centre it.
GRID_Z_CENTRE = 800.0
#: Sub-samples per cell axis when computing the expected galaxy count.
GRID_SUBSAMPLE = 2
#: Cells less than this fraction inside the survey are treated as unobserved.
GRID_MIN_COVERAGE = 0.5
#: Fiducial Gaussian smoothing of the density field and variations.
SMOOTHING_MPC_H = 10.0
SMOOTHING_VARIANTS = (5.0, 20.0)
#: k-range of the bias fit, h/Mpc.
BIAS_KMIN, BIAS_KMAX = 0.02, 0.08
#: Resolution of the predicted ISW map.
ISW_NSIDE = 256
#: Radial step of the line-of-sight integration, Mpc/h.
ISW_STEP_MPC_H = 5.0

# --------------------------------------------------------------------------
# Cross-correlations
# --------------------------------------------------------------------------

CC_NSIDE = 64
CC_LMIN, CC_LMAX, CC_DL = 10, 100, 10
NVSS_SMIN_MJY = 2.5
NVSS_DEC_MIN = -40.0
NVSS_BMIN = 10.0
#: Declination bands for the NVSS density correction.
NVSS_DEC_BAND_DEG = 1.0
BOSS_SLICES = ((0.15, 0.43), (0.43, 0.75))
#: Rotations of the CMB in Galactic longitude for the rotated-map errors.
N_ROTATIONS = 35  # every 10 deg, excluding 0

# --------------------------------------------------------------------------
# Cold Spot
# --------------------------------------------------------------------------

COLDSPOT_LB = (209.0, -57.0)  # Cruz et al. (2005); Vielva et al. (2004)
SMHW_SCALES_DEG = (3.0, 4.0, 5.0, 6.0, 8.0, 10.0)
SMHW_SCALE_DEG = 5.0

# --------------------------------------------------------------------------
# Published values we grade against
# --------------------------------------------------------------------------

LITERATURE = {
    "granett2008": {
        "ref": "Granett, Neyrinck & Szapudi 2008, ApJ 683, L99 (arXiv:0805.3695), WMAP5",
        "voids_uK": (-11.3, 3.1), "clusters_uK": (7.9, 3.1), "combined_uK": (9.6, 2.2),
        "combined_sigma": 4.4, "filter_rms_uK": 22.0,
        "table1": {(30, 4.0): 11.1, (50, 4.0): 9.6, (70, 4.0): 5.4, (50, 3.0): 8.4,
                   (50, 3.5): 9.3, (50, 4.5): 9.2, (50, 5.0): 7.8},
        "millennium_expectation_uK": 4.2,
    },
    "gkh2015": {
        "ref": "Granett, Kovacs & Hawken 2015, MNRAS 454, 2804 (arXiv:1507.03914), Table 4",
        "granett_voids_A": (6.30, 2.59), "granett_clusters_A": (4.51, 3.57),
        "cmass_template_A": (0.33, 0.17), "cmass_bias_sigma8_0.8": 2.18,
    },
    "hotchkiss2015": {
        "ref": "Hotchkiss et al. 2015, MNRAS 446, 1321 (arXiv:1405.3552), Planck 2013 SMICA, "
               "R = 0.6 theta_eff",
        "lrgdim_voids_uK": (0.14, 2.8, 70), "lrgdim_clusters_uK": (2.05, 1.9, 196),
        "lrgbright_voids_uK": (-2.6, 5.4, 13), "lrgbright_clusters_uK": (-4.0, 3.8, 39),
    },
    "nc2016": {
        "ref": "Nadathur & Crittenden 2016, ApJL 830, L19 (arXiv:1608.08638)",
        "A_ISW": (1.64, 0.53), "expected_snr_lcdm": 1.9,
        "n_voids": 2445, "n_superclusters": 29866,
    },
    "planck2015_xxi": {
        "ref": "Planck 2015 results XXI, A&A 594, A21 (arXiv:1502.01595), Table 2, SMICA",
        "NVSS_A": (0.95, 0.36), "SDSS_CMASS_LOWZ_A": (1.37, 0.56), "all_A": (1.00, 0.25),
        "nvss_dndz": {"z0": 0.33, "alpha": 0.37}, "nvss_nbar_sr": 1.584e5, "nvss_fsky": 0.73,
        "nvss_bias": "0.90 [1 + 0.54 (1+z)^2]",
    },
    "kovacs2022": {
        "ref": "Kovacs et al. (DES) 2022, MNRAS 510, 216 (arXiv:2112.07699)",
        "coldspot_centre_uK": -150.0, "eridanus_isw_centre_uK": -20.0,
        "fraction_explained": (0.10, 0.20), "supervoid_A_ISW": (5.2, 1.6),
        "eridanus_delta0": -0.2, "eridanus_R_mpc_h": 200.0, "eridanus_z": 0.15,
    },
    "szapudi2015": {
        "ref": "Szapudi et al. 2015, MNRAS 450, 288 (arXiv:1405.1566)",
        "z": (0.22, 0.03), "R_mpc_h": (220.0, 50.0), "delta_m": (-0.14, 0.04),
        "isw_estimate_uK": -20.0,
    },
    "finelli2016": {
        "ref": "Finelli et al. 2016, MNRAS 455, 1246 (arXiv:1405.1555): compensated profile "
               "delta(r) = -delta0 (1 - 2r^2/3r0^2) exp(-r^2/r0^2)",
        "r0_mpc_h": 195.0, "delta0": 0.25, "z0": 0.155,
    },
}
