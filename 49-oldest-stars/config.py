"""Project configuration: paths, data sources, analysis choices and reference values.

Everything tunable lives here, together with the published numbers we grade our
own results against. Published *ages* are only ever used for comparison. The
inputs that do enter the fits are the cluster distances (Baumgardt & Vasiliev
2021), reddenings and metallicities (Harris 2010 edition) and alpha-element
enhancements (spectroscopy), each with an uncertainty that the Monte Carlo
propagates. They are listed in CLUSTERS below with their sources.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
ATLAS_DIR = PROJECT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
HUGS_DIR = RAW_DIR / "hugs"
MIST_DIR = RAW_DIR / "mist"
BASTI_DIR = RAW_DIR / "basti"
DSED_DIR = RAW_DIR / "dsed"
LIT_DIR = RAW_DIR / "literature"
DERIVED_DIR = DATA_DIR / "derived"
FIGURE_DIR = PROJECT_DIR / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, RAW_DIR, HUGS_DIR, MIST_DIR, BASTI_DIR, DSED_DIR, LIT_DIR,
           DERIVED_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

#: HST UV Legacy Survey of Galactic Globular Clusters (HUGS; Piotto et al. 2015,
#: Nardiello et al. 2018), MAST High Level Science Product, DOI 10.17909/T9810F.
#: One catalogue per photometric method; method 1 is the best one for stars
#: bright enough to give a star-like profile in each exposure (Nardiello et al.
#: 2018, Sect. 2.2), which covers the turnoff, subgiant and giant branches here.
HUGS_URL = ("https://archive.stsci.edu/hlsps/hugs/{hugs}/"
            "hlsp_hugs_hst_wfc3-uvis-acs-wfc_{hugs}_multi_v1_catalog-meth{method}.txt")
HUGS_METHOD = 1

#: MIST v2.5 (Dotter et al. 2026, "MIST II", arXiv:2602.22012): MESA isochrones
#: with [alpha/Fe] = -0.2 ... +0.6, synthetic HST ACS/WFC photometry. The model-
#: grids page links "HST_ACSWF.txz", which is a 404; the file is HST_ACS_WFC.txz.
#: The tarball holds both rotation rates; we use the non-rotating models.
MIST_ISO_URL = "https://mist.science/data/tarballs_v2.5/isos/HST_ACS_WFC.txz"
#: MIST v2.5 bolometric-correction tables (Teff, log g, [Fe/H], [alpha/Fe], A_V at
#: R_V = 3.1). We use them for the extinction of every isochrone point, so the
#: reddening is applied to each star's own spectrum rather than a fixed ratio.
MIST_BC_URLS = {
    "HST_ACS_WFC": "https://mist.science/BC_tables/v2/HST_ACS_WFC.txz",
    "HST_WFC3": "https://mist.science/BC_tables/v2/HST_WFC3.txz",
}
MIST_ROTATION = "vvcrit0.0"

#: BaSTI-IAC (Hidalgo et al. 2018; Pietrinferni et al. 2021): pre-computed HST
#: ACS isochrones, one tarball per [Fe/H] node. "P04" = [alpha/Fe] = +0.4,
#: "P00" = scaled solar; O1D1E1 = core overshoot, atomic diffusion, Reimers
#: eta = 0.3; Y247 = the standard helium, Y275 a helium-rich variant.
BASTI_SITE = "http://basti-iac.oa-abruzzo.inaf.it"
#: The site's isochrone form reads its options from small text files; the one
#: for (grid, [Fe/H] node, photometric system) holds the tarball's link.
BASTI_OPTION = BASTI_SITE + "/textdataISOCS/{grid}{node}ACS.isc"
BASTI_GRIDS = {
    "P04": "P04O1D1E1Y247",  # [alpha/Fe] = +0.4, Y = 0.247 (baseline)
    "P00": "P00O1D1E1Y247",  # scaled solar, same physics (for [alpha/Fe] interpolation)
    "P04Y275": "P04O1D1E1Y275",  # helium-rich variant, for the helium test
}
#: [Fe/H] nodes (identical in all three grids) and their codes on the site.
BASTI_FEH_NODES = {
    -3.20: "FEHm320", -2.50: "FEHm250", -2.20: "FEHm220", -1.90: "FEHm190",
    -1.70: "FEHm170", -1.55: "FEHm155", -1.40: "FEHm140", -1.30: "FEHm130",
    -1.20: "FEHm120", -1.05: "FEHm105", -0.90: "FEHm090", -0.70: "FEHm070",
    -0.60: "FEHm060", -0.40: "FEHm040",
}

#: Dartmouth Stellar Evolution Database (DSED; Dotter et al. 2008), the 2012
#: isochrone grid in HST ACS/WFC Vega magnitudes: [alpha/Fe] = -0.2 ... +0.8.
DSED_URL = "https://rcweb.dartmouth.edu/stellar/isochrones/HST_ACSWF.tgz"

#: Harris (1996), 2010 edition: positions, [Fe/H], E(B-V), horizontal-branch
#: magnitudes for 157 Milky Way globular clusters.
HARRIS_URL = "https://physics.mcmaster.ca/~harris/mwgc.dat"

#: arXiv LaTeX sources we read tables from: {local name: (arXiv id, member)}.
ARXIV_EPRINT = "https://arxiv.org/e-print/{id}"
ARXIV_SOURCES = {
    # Baumgardt & Vasiliev (2021), MNRAS 505, 5957: final distances (with the
    # separate Gaia parallax, Gaia/HST kinematic, subdwarf, star-count and
    # literature distances) and the literature distance compilation.
    "bv21_findis_tab.tex": ("2105.09526", "findis_tab.tex"),
    "bv21_distab.tex": ("2105.09526", "distab.tex"),
    # VandenBerg et al. (2013), ApJ 775, 134: ages of 55 clusters, Table 2.
    "vandenberg2013.tex": ("1308.2257", "ms.tex"),
    # Milone et al. (2017), MNRAS 464, 3636: chromosome maps, Table 2 holds the
    # RGB widths and the first-generation fractions.
    "milone2017.tex": ("1610.00451", "paper9.tex"),
    # Massari, Koppelman & Helmi (2019), A&A 630, L4: which progenitor galaxy
    # each cluster came from (main disc / bulge, or accreted).
    "massari2019.tex": ("1906.08271", "AMR_ref_arxiv.tex"),
    # Valcin et al. (2025), JCAP: ages of 69 clusters from full-CMD fits (DSED).
    "valcin2025.tex": ("2503.19481", "main.tex"),
}

#: Earlier atlas projects this one reads (their own cached outputs).
D1_RESULTS = ATLAS_DIR / "05-supernova-cosmology" / "data" / "results.json"
D3_FIT = ATLAS_DIR / "09-planck-cmb" / "data" / "derived" / "tt_fit.json"
D8_RESULTS = ATLAS_DIR / "38-cosmic-chronometers" / "data" / "results.json"

# ---------------------------------------------------------------------------
# The clusters
# ---------------------------------------------------------------------------

#: The four clusters of the spec. Inputs that enter the fit:
#:   feh, ebv     -- Harris (2010 edition); Monte Carlo widths below.
#:   alpha        -- [alpha/Fe] = +0.4 +- 0.1, the value Ying et al. (2023) adopt
#:                   for M92 (Roederer & Sneden 2011) and Thompson et al. (2020)
#:                   for 47 Tuc; typical of halo stars and globular clusters.
#: Distances are read from Baumgardt & Vasiliev (2021) in catalog.py.
CORE_CLUSTERS = ["M92", "NGC 6397", "M15", "47 Tuc"]
CLUSTERS = {
    "M92": {"hugs": "ngc6341", "harris": "NGC 6341", "bv21": "NGC 6341"},
    "NGC 6397": {"hugs": "ngc6397", "harris": "NGC 6397", "bv21": "NGC 6397"},
    "M15": {"hugs": "ngc7078", "harris": "NGC 7078", "bv21": "NGC 7078"},
    "47 Tuc": {"hugs": "ngc0104", "harris": "NGC 104", "bv21": "NGC 104"},
}

ALPHA_PRIOR = (0.40, 0.10)
#: [Fe/H] uncertainty: the spread between spectroscopic metallicity scales
#: (Carretta et al. 2009 vs Kraft & Ivans 2003) is ~0.1 dex for these clusters.
FEH_SIGMA = 0.10
#: E(B-V) uncertainty: 10% of Harris's value, with a 0.01 mag floor.
EBV_REL_SIGMA = 0.10
EBV_MIN_SIGMA = 0.01
#: Systematic floor added in quadrature to the Baumgardt & Vasiliev distance
#: modulus. Their mean distances average 30+ measurements and quote ~0.9%
#: (0.018 mag), but the offsets they measured *between* methods are 0.021 mag
#: (RR Lyrae vs CMD fitting) and 0.025 mag (their subdwarfs vs CMD fitting):
#: errors shared by every cluster that averaging cannot shrink.
DM_SYS_FLOOR = 0.025
#: Ratio of total to selective extinction for the MIST bolometric corrections.
R_V = 3.1

# ---------------------------------------------------------------------------
# Analysis choices
# ---------------------------------------------------------------------------

#: Monte Carlo draws per cluster and model set.
N_MC = 1000
#: Bootstrap resamplings of each cluster's fiducial (ridge) line.
N_BOOT = 200
RANDOM_SEED = 20261006

#: Fine age grid on which each isochrone family is tabulated before fitting (Gyr).
AGE_GRID = (8.0, 16.0, 0.05)

#: Model sets: label -> human-readable name.
MODEL_SETS = {
    "MIST": "MIST v2.5 (MESA, 2026)",
    "BaSTI": "BaSTI-IAC (2021)",
    "DSED": "Dartmouth DSED (2008/2012)",
}

# ---------------------------------------------------------------------------
# Published values we grade our results against (never used as fit inputs)
# ---------------------------------------------------------------------------

#: Ying et al. (2023, AJ 166, 18): M92 from 20,000 Monte Carlo DSEP isochrone
#: sets fitted to the ACS survey photometry with Voronoi-binned Hess diagrams.
#: Their own main-sequence-fitting distance modulus is 14.72 +- 0.04, farther
#: than Baumgardt & Vasiliev's 14.65 +- 0.02: a 0.07 mag shorter distance makes
#: a cluster ~0.8 Gyr older.
YING2023_M92 = {"age": 13.80, "err": 0.75, "dm": 14.72, "dm_err": 0.04}

#: Other independent clocks for the same clusters.
OTHER_CLOCKS = {
    "NGC 6397": [
        # White-dwarf cooling sequence (Hansen et al. 2007, ApJ 671, 380); the
        # quoted +-0.47 Gyr is a 95% interval.
        {"label": "white-dwarf cooling (Hansen+07)", "age": 11.47, "err": 0.47 / 1.96},
    ],
    "47 Tuc": [
        # Two detached eclipsing binaries + CMD (Thompson et al. 2020, MNRAS 492, 4254).
        {"label": "eclipsing binaries (Thompson+20)", "age": 12.0, "err": 0.5},
        # White-dwarf cooling (Hansen et al. 2013, Nature 500, 51), 95% interval.
        {"label": "white-dwarf cooling (Hansen+13)", "age": 9.9, "err": 0.7 / 1.96},
    ],
}

#: Valcin et al. (2025, JCAP; arXiv:2503.19481): age of the oldest clusters from
#: 69 full-CMD fits, and the age of the universe it implies.
VALCIN2025 = {"t_gc": 13.39, "stat": 0.10, "sys": 0.23, "t_u": 13.57}
#: Valcin et al. (2026, arXiv:2603.04872): the same with two populations per cluster,
#: metal-poor ([Fe/H] < -1.5) subsample.
VALCIN2026 = {"t_gc": 13.61, "stat": 0.25, "sys": 0.23, "t_u": 13.81}
