"""Project configuration: paths, data sources, epochs, and every adopted number.

Every number this project *chooses* rather than measures lives here, next to
the reason for it: the distances that turn arcseconds per year into km/s, the
image epochs, the sizes of the patches we track, the rotation curve that turns
a pulsar's proper motion into a kick. Published results we only compare
against live in ``literature.py``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
PLATE_DIR = DATA_DIR / "plates"  # DSS photographic plates and Pan-STARRS images
GAIA_DIR = DATA_DIR / "gaia"
CHANDRA_DIR = DATA_DIR / "chandra"
CATALOG_DIR = DATA_DIR / "catalogs"
DERIVED_DIR = DATA_DIR / "derived"  # registered images, motion vectors: regenerable caches
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, PLATE_DIR, GAIA_DIR, CHANDRA_DIR, CATALOG_DIR, DERIVED_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

USER_AGENT = "astrophysics-project-atlas/L3 (hobby research; python-urllib)"

#: Julian year in days and the km/s carried by 1 mas/yr at 1 kpc.
JULIAN_YEAR_DAYS = 365.25
KMS_PER_MASYR_KPC = 4.740470


# --------------------------------------------------------------------------
# The two remnants
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Target:
    """A supernova remnant and the box of sky we fetch around it."""

    key: str
    name: str
    ra: float  # deg (ICRS): centre of every cutout
    dec: float
    box_arcmin: float  # side of the square cutouts
    distance_kpc: float
    distance_ref: str
    distance_range_kpc: tuple[float, float]


#: The Crab: cutouts centred on the pulsar (Gaia DR3 3403818172572314624 at
#: epoch 2016.0). The nebula is about 7' x 5', so a 12' box leaves a ring of
#: field stars all around it to register the plates on.
#: Distance: 2.0 kpc, the conventional value from Trimble (1973)'s expansion
#: parallax. Gaia DR3 measures the pulsar's parallax as 0.511 +- 0.079 mas,
#: 1/parallax = 1.96 kpc, i.e. 1.7-2.3 kpc at 1 sigma: consistent.
CRAB = Target("crab", "Crab Nebula (M1, SN 1054)", 83.63305926, 22.01449737, 12.0,
              2.0, "Trimble 1973, PASP 85, 579; Gaia DR3 parallax 0.511 +- 0.079 mas", (1.7, 2.3))

#: Cas A: cutouts centred on the optical expansion centre of Thorstensen, Fesen
#: & van den Bergh (2001, AJ 122, 297): 23h23m27.77s +58d48m49.4s (J2000).
#: Distance 3.4 (+0.3 -0.1) kpc from the radial velocities and proper motions
#: of its knots (Reed et al. 1995, ApJ 440, 706), the value every comparison
#: paper here uses.
CASA = Target("casa", "Cassiopeia A", 350.8657083, 58.8137222, 9.0,
              3.4, "Reed et al. 1995, ApJ 440, 706", (3.3, 3.7))

TARGETS = {t.key: t for t in (CRAB, CASA)}

#: Gaia DR3 source_ids that matter by name.
CRAB_PULSAR_GAIA = 3403818172572314624
#: The 16th-magnitude star 4" NE of the pulsar that Duncan (1921) used as the
#: reference for the first proper motions of the Crab; Martin et al. (2025)
#: quote their expansion centre relative to it.
DUNCAN_STAR_GAIA = 3403818176867563264


# --------------------------------------------------------------------------
# The images: photographic plates and modern CCDs
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Plate:
    """One image of one remnant: a scanned plate or a modern stack."""

    target: str
    key: str  # file stem under data/plates/
    survey: str  # DSS survey code, or "ps1-<band>"
    label: str  # human description for figures
    band: str  # "red", "blue", "ir": which pairs are compared
    lines: str  # the emission lines inside the bandpass


DSS_URL = "https://archive.stsci.edu/cgi-bin/dss_search"
PS1_FILENAMES_URL = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
PS1_FITSCUT_URL = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"

#: The Crab's plates, oldest first. Dates come from each file's own header
#: (DATE-OBS); the comments give them for reading convenience. POSS-I's red
#: plate was scanned with 25 um pixels (1.7"), its blue plate with 15 um (1.0").
CRAB_PLATES = [
    # 1951-11-07, Palomar 48-inch Schmidt, 103aE emulsion + red plexiglass: 6200-6700 A
    Plate("crab", "crab_poss1_red", "poss1_red", "POSS-I red (1951)", "red", "H-alpha, [N II], [S II]"),
    # 1951-11-07, same night, 103aO emulsion, no filter: 3300-5000 A
    Plate("crab", "crab_poss1_blue", "poss1_blue", "POSS-I blue (1951)", "blue", "[O II], H-beta, [O III] (edge)"),
    # 1983-01-07, Palomar Quick-V survey, IIaD + W12: 4900-6000 A
    Plate("crab", "crab_quickv", "quickv", "Quick-V (1983)", "visual", "[O III] 5007, [O I] 6300 (edge)"),
    # 1991-10-17, POSS-II, IIIaF + RG610: 6100-6900 A
    Plate("crab", "crab_poss2_red", "poss2ukstu_red", "POSS-II red (1991)", "red", "H-alpha, [N II], [S II]"),
    # 1992-12-16, POSS-II, IIIaJ + GG395: 3950-5400 A
    Plate("crab", "crab_poss2_blue", "poss2ukstu_blue", "POSS-II blue (1992)", "blue", "H-beta, [O III]"),
    # 1998-01-27, POSS-II, IV-N + RG9: 7000-9000 A
    Plate("crab", "crab_poss2_ir", "poss2ukstu_ir", "POSS-II near-IR (1998)", "ir", "[O II] 7320, [S III] 9069 (edge)"),
    # Pan-STARRS1 3pi stacks, 2009-2014 (effective epoch from the stack header)
    Plate("crab", "crab_ps1_g", "ps1-g", "Pan-STARRS g (~2012)", "blue", "H-beta, [O III]"),
    Plate("crab", "crab_ps1_r", "ps1-r", "Pan-STARRS r (~2012)", "red", "[O I], H-alpha, [N II], [S II]"),
    Plate("crab", "crab_ps1_i", "ps1-i", "Pan-STARRS i (~2012)", "ir", "[O II] 7320, [Ar III] 7136"),
]

#: Cas A's plates. Its optical knots glow in [S II], [O II], [O III] and
#: [S III]; heavy dust extinction (A_V ~ 5-6 mag) leaves nothing in the blue.
CASA_PLATES = [
    # 1954-08-07, POSS-I E
    Plate("casa", "casa_poss1_red", "poss1_red", "POSS-I red (1954)", "red", "[S II] 6716/31 (edge), H-alpha, [N II]"),
    # 1991-08-07, POSS-II F
    Plate("casa", "casa_poss2_red", "poss2ukstu_red", "POSS-II red (1991)", "red", "[S II], H-alpha, [N II]"),
    # 1991-08-05, POSS-II N
    Plate("casa", "casa_poss2_ir", "poss2ukstu_ir", "POSS-II near-IR (1991)", "ir", "[O II] 7320, [S III] 9069"),
    Plate("casa", "casa_ps1_r", "ps1-r", "Pan-STARRS r (~2012)", "red", "[O I], H-alpha, [N II], [S II]"),
    Plate("casa", "casa_ps1_i", "ps1-i", "Pan-STARRS i (~2012)", "ir", "[O II] 7320, [Ar III] 7136"),
    Plate("casa", "casa_ps1_z", "ps1-z", "Pan-STARRS z (~2012)", "ir", "[S III] 9069, 9531"),
]

PLATES = CRAB_PLATES + CASA_PLATES


# --------------------------------------------------------------------------
# Gaia DR3: the reference frame for every plate
# --------------------------------------------------------------------------

#: VizieR's ASU interface to Gaia DR3 (I/355/gaiadr3). Its TAP service and the
#: ESA archive both failed on 2026-10-04 (503 and a 90-s timeout); this plain
#: cone search answered in under a second.
VIZIER_ASU_URL = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
GAIA_COLUMNS = "Source,RA_ICRS,DE_ICRS,Plx,e_Plx,pmRA,e_pmRA,pmDE,e_pmDE,Gmag,BP-RP,RUWE"
GAIA_EPOCH = 2016.0
#: Cone radii: the corner of each square cutout, plus a margin.
GAIA_RADIUS_ARCMIN = {"crab": 9.0, "casa": 7.0}
GAIA_GMAX = 20.0


# --------------------------------------------------------------------------
# Chandra: Cas A in X-rays
# --------------------------------------------------------------------------

CHANDRA_FTP = "https://cxc.cfa.harvard.edu/cdaftp/byobsid"

#: Every public ACIS-S imaging observation of the whole remnant (no grating,
#: full frame, > 15 ks), grouped into epochs: observations a few days apart
#: are summed. From the HEASARC ``chanmaster`` table, queried 2026-10-04. The
#: pipeline's "centre image" product is a 1024 x 1024 map of every event at
#: native 0.492" pixels (all energies) -- 0.7 MB, against 523 MB for the event
#: list -- and it covers all of Cas A except the tip of the north-east jet.
#: Observations aimed at the neutron star itself (6690, 13783, 16946, 17639,
#: 22426, 23248, 26655, 27099, 27100) read only a 100- or 256-row strip of the
#: chip and are not used. 2025-26 observations are still proprietary.
#: The keys are labels; the epoch used is the exposure-weighted mean date read
#: from the image headers.
CHANDRA_EPOCHS: dict[str, tuple[int, ...]] = {
    "2000.08": (114,),
    "2002.10": (1952,),
    "2004.11": (5196,),
    "2004.29": (4634, 4635, 4636, 4637, 4638, 4639, 5319, 5320),  # the 1 Ms "Very Large Project"
    "2007.93": (9117, 9773),
    "2009.84": (10935, 12020),
    "2010.83": (10936, 13177),
    "2012.37": (14229,),
    "2013.38": (14480,),
    "2014.36": (14481,),
    "2015.33": (14482,),
    "2016.81": (19903, 18344),
    "2017.38": (19604,),
    "2018.37": (19605,),
    "2019.37": (19606,),
}
#: The deep 2004 image every other epoch is compared with, as in Vink et al. (2022).
CHANDRA_REFERENCE_EPOCH = "2004.29"
CHANDRA_PIXEL_ARCSEC = 0.492


# --------------------------------------------------------------------------
# Catalogues
# --------------------------------------------------------------------------

#: Green's catalogue of Galactic supernova remnants, 2024 October version
#: (Green 2025, JApA 46, 14): 310 remnants.
GREEN_URL = "https://www.mrao.cam.ac.uk/surveys/snrs/snrs.data.html"
GREEN_REFERENCE = "Green 2025, J. Astrophys. Astron. 46, 14 (2024 October version)"

#: ATNF Pulsar Catalogue (Manchester et al. 2005, AJ 129, 1993), release 2.8.1
#: (2026-06-18), the same pinned release as project C2.
ATNF_VERSION = "2.8.1"
ATNF_URL = f"https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.v{ATNF_VERSION}.tar.gz"
ATNF_DB = CATALOG_DIR / f"psrcat_v{ATNF_VERSION}.db"


# --------------------------------------------------------------------------
# Measuring motions in images
# --------------------------------------------------------------------------

#: Common grid every Crab image is resampled onto: a tangent plane centred on
#: the pulsar, north up, east left.
CRAB_GRID_ARCSEC = 0.5
CRAB_GRID_SIZE_ARCMIN = 11.0

#: Plate solutions: a cubic polynomial in each coordinate (10 terms), fitted to
#: Gaia stars wound back to the plate's date. Cubic is what Martin et al. (2025)
#: used for their oldest images; Schmidt plates need at least that.
PLATE_POLY_ORDER = 3
#: Stars used to register: faint enough not to be saturated blobs on the
#: plates, bright enough to centroid. Brighter stars (G < 14) sit 0.1-0.3"
#: off on every plate, the photographic "magnitude equation" (guiding and
#: coma act differently on saturated images); the filaments are faint, so we
#: register on the faint stars.
REG_GMIN, REG_GMAX = 14.0, 18.5

#: Filament tracking: square patches of this size, on a grid with this step,
#: searched over shifts up to this size. Crab filaments move up to ~0.25"/yr,
#: i.e. < 16" in the 60 years from 1951 to 2012.
PATCH_ARCSEC = 20.0
PATCH_STEP_ARCSEC = 6.0
SEARCH_ARCSEC = 22.0
#: Every image is blurred to this resolution (FWHM) before comparing: a bit
#: coarser than the worst plate (star images on POSS-I measure 4.3-4.5").
COMMON_FWHM_ARCSEC = 5.0
#: Structure larger than this is removed before correlating (the smooth
#: synchrotron glow differs between bands; filaments do not).
HIGHPASS_ARCSEC = 12.0


# --------------------------------------------------------------------------
# The Galaxy: turning pulsar proper motions into kicks
# --------------------------------------------------------------------------

#: Hobbs et al. (2005) removed Galactic rotation with a flat rotation curve,
#: R0 = 8.5 kpc, V0 = 225 km/s, and the Sun's peculiar motion of 15.6 km/s
#: towards (l, b) = (48.8, 26.3) deg. We use the same model to compare like
#: with like.
HOBBS_R0_KPC = 8.5
HOBBS_V0_KMS = 225.0
HOBBS_SUN_PECULIAR = (15.6, 48.8, 26.3)  # km/s, l, b
#: Their "young pulsars" have characteristic ages below 3 Myr.
YOUNG_PULSAR_AGE_YR = 3e6


@dataclass(frozen=True)
class Source:
    """One file we download, and why."""

    filename: str  # path under data/
    url: str
    what: str
    reference: str
    kind: str = "text"  # payload check: "fits", "gz", "text", "html", "tsv"
    pause: float = 0.0
    timeout: float = 300.0
    params: dict = field(default_factory=dict)
