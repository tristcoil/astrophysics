"""Project configuration: paths, data sources, priors, model settings.

Every choice that changes a number quoted in the README lives here, once:
which papers the astrometry comes from, which points are left out and why,
the priors of each orbit model, the N-body stability criteria, and the random
seeds. Change a setting here and every fit and figure follows it.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
RAW = DATA / "raw"              # downloads, never modified
DERIVED = DATA / "derived"      # parsed tables (astrometry.csv, ...)
CACHE = DATA / "cache"          # expensive results: nested-sampling runs, N-body runs
FIGURES = HERE / "figures"
MANIFEST = DATA / "manifest.json"
RESULTS = DATA / "results.json"

for _d in (DATA, RAW, DERIVED, CACHE, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = 8799

# ---------------------------------------------------------------------------
# Physical constants (IAU 2015 nominal values; the same ones astropy uses)
# ---------------------------------------------------------------------------

GM_SUN = 1.3271244e20            # m^3 s^-2
AU = 1.495978707e11              # m
DAY = 86400.0                    # s
YEAR = 365.25 * DAY              # Julian year, s
MJUP_MSUN = 1.2668653e17 / GM_SUN   # Jupiter mass in solar masses (9.546e-4)
#: G in au^3 Msun^-1 yr^-2 for the REBOUND simulations, from the constants
#: above (39.4769; not 4 pi^2, which belongs to the sidereal year).
G_AU_YR = GM_SUN * YEAR**2 / AU**3

# ---------------------------------------------------------------------------
# The HR 8799 system
# ---------------------------------------------------------------------------

PLANETS = ("b", "c", "d", "e")          # outermost first, as in the papers
INNER_TO_OUTER = ("e", "d", "c", "b")
HIP_ID = 114189
GAIA_DR3_SOURCE_ID = 2832463659640297472
#: Wang et al. (2018) quote tau, the epoch of periastron in units of the
#: period, after MJD 50000. We use the same reference so their tables compare
#: directly. (orbitize! defaults to MJD 58849.)
TAU_REF_MJD = 50000.0
#: Orbital phases are sampled as mean longitudes lambda = Omega + omega + M at
#: this epoch (2013 Jul 27, mid-way through the precise data): for these nearly
#: face-on orbits lambda is almost the measured position angle, so the posterior
#: is compact. A uniform prior on lambda is the same as Wang's uniform tau prior.
LAMBDA_REF_MJD = 56500.0
#: Epoch of the sampled sky states (positions and velocities) - see models.py.
STATE_REF_MJD = 56500.0
#: The date Wang et al. (2018) start their N-body integrations from.
NBODY_EPOCH_MJD = 56609.0

# ---------------------------------------------------------------------------
# Data sources. The astrometry is parsed from the LaTeX sources of the papers
# on arXiv (fetch_data.py), so every number traces to its publication.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ArxivSource:
    key: str          # short label used in the astrometry table
    arxiv: str
    citation: str
    files: tuple[str, ...]   # the files inside the e-print we read
    what: str

    @property
    def url(self) -> str:
        return f"https://arxiv.org/e-print/{self.arxiv}"

    @property
    def folder(self) -> Path:
        return RAW / "arxiv" / self.arxiv


ARXIV = {
    s.key: s
    for s in [
        ArxivSource("Konopacky+2016", "1604.08157", "Konopacky et al. 2016, AJ 152, 28",
                    ("tab2.tex", "ms.tex"),
                    "Keck/NIRC2 astrometry of all four planets, 2004-2014, one pipeline"),
        ArxivSource("Wang+2018", "1809.04107", "Wang et al. 2018, AJ 156, 192",
                    ("main_clean.tex",),
                    "GPI astrometry 2013-2016; orbit fits, stable orbits (our answer key)"),
        ArxivSource("DeRosa+2020", "1910.08659", "De Rosa et al. 2020, JATIS 6, 015006",
                    ("article.tex",),
                    "Revised GPI astrometric calibration and revised HR 8799 astrometry"),
        ArxivSource("Zurlo+2022", "2207.10684", "Zurlo et al. 2022, A&A 666, A133",
                    ("hr8799.tex",),
                    "VLT/SPHERE IRDIS + IFS astrometry 2014-2021, LBT/LUCI 2020"),
        ArxivSource("Soummer+2011", "1110.1382", "Soummer et al. 2011, ApJ 741, 55",
                    ("ms.tex",),
                    "HST/NICMOS 1998 astrometry of b, c and d"),
        ArxivSource("Lafreniere+2009", "0902.3247", "Lafreniere et al. 2009, ApJ 694, L148",
                    ("ms.tex",),
                    "HST/NICMOS 1998 detection of b (epoch of the 1998 images)"),
        ArxivSource("GRAVITY+2019", "1903.11903", "GRAVITY Collaboration (Lacour et al.) 2019, A&A 623, L11",
                    ("HR.tex",),
                    "VLTI/GRAVITY astrometry of e (2018), ~0.1 mas"),
        ArxivSource("Balmer+2025", "2503.13608", "Balmer et al. 2025, AJ 169, 209",
                    ("hr8799bar.tex",),
                    "JWST/NIRCam coronagraphy of all four planets, 2023 Nov 5"),
        ArxivSource("Nguyen+2026", "2609.10507", "Nguyen et al. 2026, arXiv:2609.10507",
                    ("main.tex",),
                    "JWST/NIRISS AMI positions (2023 Aug 3) and Keck/NIRC2 2025 Oct 9"),
        ArxivSource("Brandt+2021", "2105.12820", "Brandt et al. 2021, ApJL 915, L16",
                    ("main.tex",),
                    "First dynamical mass of HR 8799 e from the HGCA"),
        ArxivSource("Lagrange+2026", "2609.20996", "Lagrange et al. 2026, arXiv:2609.20996",
                    ("ms.tex",),
                    "A claimed fifth, inner companion from Gaia + HGCA"),
        ArxivSource("Wang+2016", "1607.05272", "Wang et al. 2016, AJ 152, 97",
                    ("*.tex",),
                    "beta Pic b orbit and its 2017 Hill-sphere transit prediction"),
        ArxivSource("Nielsen+2020", "1911.11273", "Nielsen et al. 2020, AJ, doi:10.3847/1538-3881/ab5b92",
                    ("*.tex",),
                    "beta Pic b astrometry (orbitize! example data) and dynamical mass"),
    ]
}

#: orbitize!'s bundled example data: beta Pic b relative astrometry (Nielsen
#: et al. 2020, Table 1) plus the CRIRES planet RV (Snellen et al. 2014).
ORBITIZE_BETAPIC = "betaPic.csv"
ORBITIZE_BETAPIC_URL = ("https://raw.githubusercontent.com/sblunt/orbitize/main/"
                        "orbitize/example_data/betaPic.csv")
#: orbitize! also bundles the Gaia EDR3 edition of the Hipparcos-Gaia Catalog
#: of Accelerations (Brandt 2021, ApJS 254, 42; CDS J/ApJS/254/42).
HGCA_URL = "https://cdsarc.cds.unistra.fr/ftp/J/ApJS/254/42/HGCA_vEDR3.fits"

#: whereistheplanet's current HR 8799 posterior (orbitize!, fitted to data that
#: include unpublished GRAVITY astrometry). Used only as an outside check of
#: predicted positions, never as data. Git LFS object (22 kB).
WITP_URL = ("https://media.githubusercontent.com/media/semaphoreP/whereistheplanet/"
            "master/data/post_hr8799.hdf5")

GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
VIZIER_TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"

# ---------------------------------------------------------------------------
# Astrometric data sets
# ---------------------------------------------------------------------------

#: Konopacky et al. (2016) left these out of their own orbit fits (table notes):
#: d on 2007 Aug 2 sits on the edge of the 1000 mas focal-plane mask, and the
#: 2009 Aug 1 epoch is two days after 2009 Jul 30. We follow them.
KONOPACKY_EXCLUDE = {("2007-08-02", "d"), ("2009-08-01", "b"), ("2009-08-01", "c"),
                     ("2009-08-01", "d"), ("2009-08-01", "e")}

#: Data set "W18": exactly what Wang et al. (2018) fitted, GPI with its
#: original calibration plus the Keck astrometry of Konopacky et al. (2016).
W18_SOURCES = ("NIRC2-K16", "GPI-W18")

#: Data set "ALL": every precise, self-calibrated set we found, 1998-2025.
#: GPI uses De Rosa et al.'s revised calibration. JWST/AMI positions (10-33 mas)
#: and other heterogeneous literature points are kept only as checks.
ALL_SOURCES = ("NICMOS-S11", "NIRC2-K16", "GPI-D20", "IRDIS-Z22", "IFS-Z22",
               "LUCI-Z22", "GRAVITY-L19", "NIRCam-B25", "NIRC2-N26")

#: Data set "PRE2017": the training set of the prediction test. Same as W18;
#: everything measured after 2016 Sep 19 is predicted, not fitted.
PREDICTION_CUTOFF_MJD = 57651.0   # 2016-09-20

#: Astrometric calibration of each instrument: prior widths of a north-angle
#: rotation (deg) and a plate-scale error (fraction) fitted as nuisance
#: parameters in the ALL fits ("fit offsets between instruments"). The widths
#: are the calibration uncertainties each team published. NICMOS errors
#: already include 0.5 deg of north uncertainty (Soummer et al. 2011), and
#: GRAVITY's baselines set its scale and orientation, so neither gets offsets.
CALIBRATION = {
    # source         sigma_north_deg  sigma_scale  reference
    "NIRC2-K16":  (0.009, 0.002 / 9.952, "Yelda et al. 2010"),
    "GPI-D20":    (0.14, 0.021 / 14.161, "De Rosa et al. 2020"),
    "IRDIS-Z22":  (0.08, 0.004, "Maire et al. 2016 (SPHERE)"),
    "IFS-Z22":    (0.08, 0.004, "Maire et al. 2016 (SPHERE)"),
    "LUCI-Z22":   (0.5, 0.01, "no published calibration; wide"),
    "NIRCam-B25": (0.1, 0.005, "JWST NIRCam; wide"),
    "NIRC2-N26":  (0.020, 0.005 / 9.971, "Service et al. 2016"),
}

# ---------------------------------------------------------------------------
# Priors (Wang et al. 2018, Section 3, unless noted)
# ---------------------------------------------------------------------------

W18_PARALLAX = (24.76, 0.64)        # mas, Gaia DR1 (Gaia Collaboration 2016)
W18_MSTAR = (1.52, 0.15)            # Msun, Baines et al. 2012 + 10%
SMA_RANGE = (1.0, 100.0)            # au, log-uniform
ECC_RANGE = (1e-6, 0.999)
#: Near 1:2:4:8 models: uniform priors on period ratios (Wang Sec. 3.3).
PERIOD_RATIO_RANGES = {"b/c": (1.8, 2.4), "c/d": (1.8, 2.4), "d/e": (1.8, 2.2)}
LOW_E_MAX = 0.2                     # "Near 1:2:4:8 Coplanar Low-e" model

#: ALL fits: Gaia DR3 parallax 24.4620 +- 0.0455 mas. HR 8799 (G ~ 5.8) is
#: brighter than Lindegren et al.'s (2021) zero-point recipe covers (G > 6),
#: so instead of correcting we widen the error to 0.10 mas, which covers a
#: typical bright-star zero point (a few hundredths of a mas) and the known
#: underestimate of DR3 errors for bright stars. At 0.4% it is negligible
#: next to the orbits' own uncertainties.
DR3_PARALLAX = (24.462, 0.10)
#: ALL fits: a broad stellar-mass prior so the planets weigh their star.
ALL_MSTAR = (1.5, 0.3)

# ---------------------------------------------------------------------------
# Orbit models
# ---------------------------------------------------------------------------

MODELS = ("unconstrained", "coplanar", "near1248", "near1248_coplanar")
MODEL_LABEL = {
    "unconstrained": "Unconstrained",
    "coplanar": "Coplanar",
    "near1248": "Near 1:2:4:8",
    "near1248_coplanar": "Near 1:2:4:8 coplanar, low e",
}

#: Nested sampling (dynesty) settings.
NLIVE = 1000
DLOGZ = 0.05
N_PROCESSES = 12

# ---------------------------------------------------------------------------
# N-body stability test (Wang et al. 2018, Section 4.1)
# ---------------------------------------------------------------------------

STABILITY_YEARS = 40e6              # integrate the present backwards 40 Myr
DT_FRACTION_OF_PE = 0.01            # fixed step: 1% of e's initial period
EJECTION_AU = 500.0
#: Planet masses: c, d, e equal, uniform 4-11 MJup; b uniform 3 MJup .. M_cde.
MASS_CDE_RANGE = (4.0, 11.0)
MASS_B_MIN = 3.0
#: Luminosity-based (hot-start) masses for re-weighting the stable orbits,
#: Baraffe et al. 2003 tracks at 42 +- 5 Myr as computed by Wang et al. (2018).
HOTSTART_MASS = {"b": (6.0, 0.7), "cde": (8.7, 1.0)}
#: Snapshots for resonant angles: every 200 yr, as Wang et al. did.
SNAPSHOT_YEARS = 200.0
