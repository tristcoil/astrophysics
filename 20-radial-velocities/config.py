"""Project configuration: paths, data sources, targets, instrument epochs, cuts.

Every selection criterion used anywhere in the project is defined here, once,
so that each number quoted in the README or printed by ``run.py`` traces back to
a single line. Change a cut here and every figure and summary follows it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
CACHE_DIR = DATA_DIR / "cache"  # expensive intermediate results (MCMC chains)
FIGURE_DIR = PROJECT_DIR / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"  # provenance of every download
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, RAW_DIR, CACHE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

CDS_FTP = "https://cdsarc.cds.unistra.fr/ftp/"


@dataclass(frozen=True)
class CdsCatalog:
    """A VizieR catalogue whose files are fetched whole from the CDS FTP mirror."""

    key: str
    catalog: str  # VizieR identifier, e.g. "J/AJ/153/208"
    files: tuple[str, ...]
    citation: str
    description: str

    def url(self, filename: str) -> str:
        return f"{CDS_FTP}{self.catalog}/{filename}"

    def local(self, filename: str) -> Path:
        return RAW_DIR / f"{self.key}_{filename}"


#: Keck/HIRES velocities of the Lick-Carnegie Exoplanet Survey (1996-2014).
KECK = CdsCatalog(
    key="keck_butler2017",
    catalog="J/AJ/153/208",
    files=("ReadMe", "table1.dat.gz", "table2.dat", "table3.dat"),
    citation="Butler et al. 2017, AJ 153, 208",
    description="60,949 Keck/HIRES velocities of 1,624 stars, with S and H-alpha "
                "activity indices; plus the team's own list of periodic signals",
)

#: The same Keck velocities with small nightly zero-point errors removed. Used
#: only as a cross-check: our fits carry their own instrument offsets.
KECK_TALOR = CdsCatalog(
    key="keck_talor2019",
    catalog="J/MNRAS/484/L8",
    files=("ReadMe", "table1.dat.gz"),
    citation="Tal-Or et al. 2019, MNRAS 484, L8",
    description="Keck/HIRES velocities corrected for nightly zero-point drifts "
                "and the August 2004 jump",
)

#: The 25-year Lick planet search with the Hamilton spectrograph (1987-2011).
LICK = CdsCatalog(
    key="lick_fischer2014",
    catalog="J/ApJS/210/5",
    files=("ReadMe", "table1.dat", "table2.dat"),
    citation="Fischer, Marcy & Spronck 2014, ApJS 210, 5",
    description="14,000+ Lick/Hamilton velocities of 386 stars, tagged with the "
                "CCD dewar in use",
)

#: HARPS velocities re-derived with SERVAL and corrected for nightly zero
#: points (Trifonov et al. 2020); Ver. 02 adds Ca II H&K (Perdelwitz et al.
#: 2024). One CSV per star. The old MPIA address redirects here.
RVBANK_URL = ("https://exo-restart.com/wp-content/uploads/Stellar_parameters/"
              "Ver_02/{name}_RVs/{name}_HARPS_RVBank_ver02.csv")
RVBANK_CITATION = ("Trifonov et al. 2020, A&A 636, A74; Perdelwitz et al. 2024, "
                   "A&A 683, A125 (HARPS RVBank Ver. 02)")

#: The ELODIE archive at Observatoire de Haute-Provence: the spectrograph
#: 51 Peg b was discovered with. The cross-correlation table gives one velocity
#: per spectrum; the spectrum's FITS header gives its barycentric time.
ELODIE_CCF_URL = "http://atlas.obs-hp.fr/elodie/fE.cgi?n=e501&c=o&o={name}&a=csv"
ELODIE_HEADER_URL = "http://atlas.obs-hp.fr/elodie/fE.cgi?n=e500&c=i&z=fd&o=elodie:{night}/{imanum}"
ELODIE_CITATION = "Moultaka et al. 2004, PASP 116, 693 (the ELODIE archive)"
#: Seconds between header requests: the archive is a small public service.
ELODIE_POLITE_DELAY_S = 0.3

#: NASA Exoplanet Archive TAP service. ``pscomppars`` for the population,
#: ``ps`` for every published solution of our own planets.
NASA_TAP_URL = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
NASA_POPULATION_COLUMNS = [
    "pl_name", "hostname", "discoverymethod", "disc_year", "pl_controv_flag",
    "pl_orbper", "pl_bmasse", "pl_bmasseerr1", "pl_bmasseerr2", "pl_bmasselim",
    "pl_bmassprov", "pl_bmasse_reflink", "st_mass",
]
NASA_PS_COLUMNS = [
    "pl_name", "hostname", "default_flag", "pl_refname", "pl_orbper",
    "pl_orbpererr1", "pl_orbpererr2", "pl_rvamp", "pl_rvamperr1", "pl_rvamperr2",
    "pl_orbeccen", "pl_orbeccenerr1", "pl_msinij", "pl_msinie", "pl_msinieerr1",
    "pl_msinieerr2", "pl_bmasse", "pl_bmassprov", "pl_rade", "pl_radeerr1",
    "pl_radeerr2", "pl_orbincl", "pl_tranmid", "pl_tranmiderr1", "st_mass",
    "st_masserr1", "st_rad", "st_raderr1",
]

#: TESS light curves (SPOC 2-minute cadence) for the one transiting planet.
TESS_AUTHOR = "SPOC"
TESS_EXPTIME_S = 120
TESS_DIR = RAW_DIR / "tess"

#: Mass-radius curves of Zeng et al. (2019, PNAS 116, 9723): what radius a
#: planet of given mass and composition would have. Models, not data.
ZENG_BASE_URL = "https://lweb.cfa.harvard.edu/~lzeng/tables/"
ZENG_MODELS = {
    "iron": ("massradiusFe.txt", "pure iron"),
    "earthlike": ("massradiusEarthlikeRocky.txt", "Earth-like rock (32% iron)"),
    "rock": ("massradiusmgsio3.txt", "pure rock (no iron)"),
}

#: Limb-darkening coefficients for the TESS band (quadratic law, ATLAS model
#: atmospheres, least-squares method), used for 55 Cnc e's transit.
CLARET = CdsCatalog(
    key="claret2017",
    catalog="J/A+A/600/A30",
    files=("ReadMe", "table25.dat.gz"),
    citation="Claret 2017, A&A 600, A30",
    description="TESS limb-darkening coefficients from ATLAS model atmospheres",
)


# ---------------------------------------------------------------------------
# Instrument epochs: where one instrument's zero point stops being another's
# ---------------------------------------------------------------------------

#: HIRES got a new CCD mosaic in August 2004. The Keck table has a 61-day gap
#: with no observations of any star from 2004-08-23 to 2004-10-23; the split
#: sits inside it. (Tal-Or et al. 2019 measured the jump at 1.5 +- 0.1 m/s.)
HIRES_UPGRADE_BJD = 2453270.0

#: HARPS fibres were exchanged in May 2015. RVBank splits its SERVAL templates
#: between 2015-05-16 and 2015-05-21 (read off stars observed on both sides);
#: none of our targets was observed within months of either date.
HARPS_FIBRE_BJD = 2457163.0

#: Lick's Doppler precision improved from ~10 to ~3 m/s by 1995 (Fischer et
#: al. 2014). Before 1995 four dewars were in use at once; afterwards one dewar
#: at a time (13: 1995-97, 6: 1998-2001, 8: 2002-11), each its own zero point.
LICK_EARLY_END_BJD = 2449718.5  # 1995-01-01


# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Star:
    key: str
    name: str  # display name
    archive_host: str  # NASA Exoplanet Archive hostname
    ra_deg: float  # J2000, used to convert Lick JD(UTC) to BJD(TDB)
    dec_deg: float
    keck: str | None = None  # name in the Keck/HIRES table
    lick: str | None = None  # name in the Lick table
    rvbank: str | None = None  # name in HARPS RVBank
    elodie: str | None = None  # name in the ELODIE archive
    tess_target: str | None = None


STARS: dict[str, Star] = {
    "51peg": Star("51peg", "51 Pegasi", "51 Peg", 344.366583, 20.768833,
                  keck="HD217014", lick="217014", rvbank="GJ882", elodie="HD217014"),
    "55cnc": Star("55cnc", "55 Cancri", "55 Cnc", 133.149212, 28.330821,
                  keck="HD75732", lick="75732", rvbank="HD75732", tess_target="55 Cnc"),
    "proxima": Star("proxima", "Proxima Centauri", "Proxima Cen", 217.428953, -62.679484,
                    rvbank="GJ551"),
    "gj581": Star("gj581", "Gliese 581", "GJ 581", 229.861767, -7.722267,
                  rvbank="GJ581"),
}


# ---------------------------------------------------------------------------
# Quality cuts
# ---------------------------------------------------------------------------

#: HARPS RVBank: keep only spectra SERVAL flags as normal (flag 0). 51 Peg has
#: 50 spectra flagged 128 ("too high S/N", i.e. near saturation).
RVBANK_GOOD_FLAG = 0
#: Drop any velocity whose quoted error is this many times the star's median
#: error for that instrument (clouds, twilight, wrong exposure time).
MAX_ERROR_OVER_MEDIAN = 5.0
#: ELODIE: keep spectra taken with simultaneous thorium ("OBTH"), whose
#: instrumental drift is measured, cross-correlated with the K0 mask used for
#: most of the 51 Peg spectra. Four spectra used the F0 mask, whose zero point
#: differs, and 19 were taken without the thorium reference.
ELODIE_MASK = "R37000K0"
ELODIE_IMATYP_PREFIX = "OBTH"
#: HARPS line-shape cut: drop spectra whose SERVAL differential line width
#: (dLW) or chromatic index (CRX) is this many robust standard deviations from
#: the star's median (per instrument epoch). Distorted line profiles --
#: contamination, very low signal -- give velocities that are not the star's
#: motion: Proxima's worst, on 2017-09-24/25, have dLW ~ 350 and CRX ~ -700 to
#: -1060 (typical: +-10 and 20-40) and velocities 37-39 m/s off.
LINE_SHAPE_NSIGMA = 6.0
#: Proxima flares: drop spectra whose H-alpha index is this many robust
#: standard deviations above the median. Flares brighten the chromosphere and
#: distort the line profiles the velocities come from.
FLARE_HALPHA_NSIGMA = 4.0


# ---------------------------------------------------------------------------
# Periodograms and fitting
# ---------------------------------------------------------------------------

#: Shortest period searched. 55 Cnc e orbits in 0.74 d, so the search must go
#: below one day -- the Keck team's own search stopped at 1 d and found the
#: 2.8 d alias instead (Butler et al. 2017, table 2, note 6).
PERIOD_MIN_D = 0.5
#: Frequency grid oversampling: points per 1/baseline.
OVERSAMPLE = 10
#: A signal is "significant" below this analytic false-alarm probability.
FAP_THRESHOLD = 1e-3

#: Stellar-jitter prior: uniform from 0 to this, in m/s.
JITTER_MAX = 100.0
#: Eccentricity is limited to e < E_MAX (closer orbits would cross the star).
E_MAX = 0.95

#: MCMC (emcee) settings. Chains are cached in data/cache/.
MCMC_WALKERS_PER_DIM = 3
MCMC_STEPS = 12000
MCMC_BURN = 4000
MCMC_THIN = 10
RANDOM_SEED = 19951006  # 51 Peg b announced in Florence, 6 October 1995

#: Adopted host-star masses (M_sun), errors and sources, all rows of the NASA
#: archive's stellarhosts table (checked against the cached table at run time).
#: m sin i scales as M_star^(2/3), so a 10% mass error is a 7% planet-mass error.
STELLAR_MASS = {
    "51peg": (1.07, 0.05, "Cont et al. 2026"),
    # Interferometric radius + isochrones; also used by Bourrier et al. 2018,
    # whose 55 Cnc e results we compare with. Published masses range from
    # 0.873 (Dai et al. 2019) to 1.015 (Moutou et al. 2026) -- a systematic
    # of up to +7.6% on every 55 Cnc planet mass, reported alongside.
    "55cnc": (0.905, 0.015, "von Braun et al. 2011"),
    "proxima": (0.1221, 0.0022, "Faria et al. 2022"),
    "gj581": (0.295, 0.010, "von Stauffenberg et al. 2024"),
}
STELLAR_MASS_ALTERNATIVE_55CNC = (1.015, 0.051, "Moutou et al. 2026")
#: 55 Cnc's radius (R_sun) from interferometry, for 55 Cnc e's size.
STELLAR_RADIUS_55CNC = (0.943, 0.010, "von Braun et al. 2011")
#: 55 Cnc's atmosphere, for limb darkening.
TEFF_LOGG_FEH_55CNC = (5196.0, 4.45, 0.31, "von Braun et al. 2011")


# ---------------------------------------------------------------------------
# Published results we compare against (not in the NASA archive's ps table)
# ---------------------------------------------------------------------------

#: 51 Peg b's own orbital velocity, from its molecular absorption lines seen at
#: high resolution (CRIRES/VLT). With our stellar K it gives the true mass
#: ratio, with no sin i: m_p / M_star = K_star / K_p.
KP_51PEG = {
    "Brogi et al. 2013 (CO + H2O, ApJ 767, 27)": (134.1, 1.8, 1.8),  # km/s, +err, -err
    "Birkby et al. 2017 (H2O, AJ 153, 138)": (133.0, 4.3, 3.5),
}

#: The discovery paper's own solution (Mayor & Queloz 1995, Nature 378, 355).
#: Also in the archive's ps table; kept here so the figure caption can quote it.
MQ95 = {"P": 4.2293, "P_err": 0.0011, "K": 59.0, "K_err": 3.0, "msini_mjup": 0.47, "N": 35}
#: ...and the date it was announced (Florence, 6 October 1995). ELODIE spectra
#: taken before this are "the discovery data" in the figures.
DISCOVERY_ANNOUNCED_BJD = 2449996.5

#: The retracted GJ 581 d. First reported at 83.4 d (Udry et al. 2007, A&A 469,
#: L43), moved to 66.8 d once that was recognised as a one-year alias (Mayor et
#: al. 2009, A&A 507, 487). Hatzes (2013, AN 334, 616, table 4) fitted HARPS
#: + HIRES with an eccentric orbit for d and got the values below (a circular
#: fit gives K = 1.74 +- 0.17 m/s); he still counted four planets. Robertson et
#: al. (2014, Science 345, 440) found a 130 +- 2 d stellar rotation in H-alpha
#: and showed d is an artefact of activity.
GJ581D_PUBLISHED = {"P": 66.64, "P_err": 0.08, "K": 2.20, "K_err": 0.19,
                    "P_udry2007": 83.4, "P_mayor2009": 66.8}
GJ581_ROTATION_D = 130.0

#: Proxima Centauri's rotation period from ASAS photometry: 82.5 d (Kiraga &
#: Stepien 2007), 83.2 d (Suarez Mascareno et al. 2016).
PROXIMA_ROTATION_D = 83.0
