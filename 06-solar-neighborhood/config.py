"""Project configuration: paths, sample cuts, and the stars worth naming.

Every tunable number in the project lives here, each with the reason it has the
value it has. Nothing in the analysis should need a magic constant of its own.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
FIGURE_DIR = PROJECT_DIR / "figures"

#: Raw archive downloads. The Gaia and GCNS queries are split into sky chunks
#: (see fetch_data.py), each cached on its own so an interrupted download resumes.
GAIA_CHUNK_DIR = DATA_DIR / "gaia_chunks"
GCNS_CHUNK_DIR = DATA_DIR / "gcns_chunks"
GAIA_PATH = DATA_DIR / "gaia_dr3_100pc.fits"
GCNS_PATH = DATA_DIR / "gcns_100pc.fits"
HIPPARCOS_PATH = DATA_DIR / "hipparcos_100pc.fits"
NAKED_EYE_PATH = DATA_DIR / "hipparcos_naked_eye.fits"
CENSUS_10PC_PATH = DATA_DIR / "census_10pc.fits"
DWARF_SEQUENCE_PATH = DATA_DIR / "mamajek_dwarf_sequence.txt"
FAMOUS_PATH = DATA_DIR / "simbad_famous.json"
NAMES_PATH = DATA_DIR / "simbad_names.fits"
EXOPLANET_HOSTS_PATH = DATA_DIR / "exoplanet_hosts.fits"
CLUSTER_REFERENCE_PATH = DATA_DIR / "simbad_clusters.json"
PROVENANCE_PATH = DATA_DIR / "provenance.json"

#: Derived products written by the analysis.
CATALOGUE_PATH = DATA_DIR / "catalogue.fits"
RESULTS_PATH = DATA_DIR / "results.json"
MAP_HTML_PATH = FIGURE_DIR / "solar_neighborhood_3d.html"

for _d in (DATA_DIR, FIGURE_DIR, GAIA_CHUNK_DIR, GCNS_CHUNK_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# The sample
# --------------------------------------------------------------------------

#: The edge of the map. 100 pc is 326 light years: far enough to contain a
#: third of a million stars and the nearest open cluster, close enough that
#: interstellar dust is negligible and Gaia sees almost everything.
MAX_DISTANCE_PC = 100.0
MIN_PARALLAX_MAS = 1000.0 / MAX_DISTANCE_PC

#: Parallax must be measured to better than 10%. At that precision distance =
#: 1/parallax is an honest estimator (its bias is under 1%), so no Bayesian
#: distance prior is needed and every distance can be traced straight to data.
MIN_PARALLAX_OVER_ERROR = 10.0

#: Rybizki et al. (2022) trained a classifier to recognise spurious Gaia
#: astrometric solutions -- mostly faint sources in crowded fields whose
#: transits were mismatched, which can masquerade as large parallaxes. They
#: recommend fidelity_v2 > 0.5 as "good". This is used instead of the common
#: RUWE < 1.4 cut, which would also throw away 16% of the real nearby stars
#: (measured on this sample), most of them genuine binaries: fine for a clean
#: HR diagram, wrong for a census.
MIN_FIDELITY = 0.5

#: Hipparcos stars are only used where Gaia has no usable measurement, and only
#: with the same 10% precision requirement.
MIN_HIPPARCOS_PARALLAX_OVER_ERROR = 10.0

#: Components of a Hipparcos double (same CCDM number) whose parallaxes agree
#: within this many sigma are treated as a physical pair, at one distance.
PAIR_CONSISTENCY_SIGMA = 3.0

#: Match radius when looking for a Gaia counterpart of a Hipparcos star that the
#: official cross-match left unpaired, after propagating the Hipparcos position
#: from epoch 1991.25 to Gaia's 2016.0.
HIPPARCOS_MATCH_RADIUS_ARCSEC = 3.0

#: Inner zone where every star gets looked up by name in SIMBAD, and shown with
#: full detail in the 3D map. 20 pc is the radius of the most complete modern
#: census of the solar neighbourhood (Kirkpatrick et al. 2024).
NAMED_RADIUS_PC = 20.0

#: Dark-sky naked-eye limit, in Johnson V. 6.0 is the conventional figure; a
#: truly pristine site reaches ~6.5. Used only to ask which stars you could
#: actually see.
NAKED_EYE_V = 6.0

# --------------------------------------------------------------------------
# Classification
# --------------------------------------------------------------------------

#: A star is counted as evolved (subgiant or giant) when it is this many
#: magnitudes brighter than a main-sequence dwarf of the same colour. An
#: unresolved pair of identical stars sits only 0.75 mag above the sequence, so
#: 1.5 mag is comfortably beyond what a binary can do...
GIANT_OFFSET_MAG = 1.5

#: ...and it must also be intrinsically luminous. Every subgiant and giant is
#: brighter than M_G ~ 4.5, whereas young or triple M dwarfs can also sit 1.5 mag
#: above the main sequence but are always far fainter than that.
GIANT_MAX_MG = 4.5

#: White dwarfs sit in their own sequence far below the main sequence. This line
#: runs through the empty gap between the two, M_G = slope * (BP-RP) + intercept;
#: everything fainter is a white dwarf. See the HR-diagram figure for where it falls.
WD_LINE_SLOPE = 3.0
WD_LINE_INTERCEPT = 10.0

#: White dwarfs are never very red: 99.9% of the Gaia Catalogue of Nearby Stars
#: white dwarfs have BP-RP < 1.88. Without this ceiling, brown dwarfs -- whose
#: faint BP flux gives them a spuriously blue BP-RP near 3 -- fall below the
#: line above and would be counted as white dwarfs.
WD_MAX_BP_RP = 2.0

#: Below this absolute magnitude a Gaia source is treated as an ultracool (L/T)
#: dwarf regardless of its colour. L0 dwarfs sit at M_G ~ 16.6 (Mamajek table),
#: and Gaia's BP flux for such faint, red objects is unreliable enough that
#: colour alone would misplace them.
ULTRACOOL_MG = 16.5

#: A source with no Gaia colour but fainter than this is still classifiable:
#: nothing but a late-M or cooler dwarf is that intrinsically faint and that
#: common (M8V is M_G ~ 15.2). Its temperature is then read from its absolute
#: magnitude. Brighter colourless sources stay unclassified.
LATE_M_MG = 15.0

#: A faint star whose BP/RP photometry is contaminated by a neighbour can land
#: in an impossible place in the HR diagram, far below the main sequence. We
#: drop its colour -- and read its temperature from its absolute magnitude
#: instead -- only when both are true: the corrected excess factor C* exceeds
#: this many sigma (Riello et al. 2021) AND the star sits more than
#: CONTAMINATED_BELOW_MS_MAG below the main sequence, where no normal star can.
#: 95% of the stars in that impossible region carry the flag.
CONTAMINATION_SIGMA = 3.0
CONTAMINATED_BELOW_MS_MAG = 2.0

#: Gaia colours of the very brightest stars suffer from saturation: against
#: Hipparcos B-V they scatter by 0.10 mag at G = 3-4 and 0.27 mag below G = 3,
#: versus 0.03 mag for G > 5. For Gaia stars brighter than this we take the
#: temperature from Hipparcos B-V instead, when Hipparcos observed the star.
BRIGHT_COLOUR_G = 4.0

# --------------------------------------------------------------------------
# Structure finding
# --------------------------------------------------------------------------

#: Neighbour-count radius for finding overdensities in 3D. The Hyades core
#: radius is ~3 pc, so this is matched to the structure we expect to find.
CLUSTER_SEARCH_RADIUS_PC = 3.0

#: The Hyades' tidal radius is ~9 pc (Roeser et al. 2011). Stars within this
#: distance of a detected cluster centre are counted as "in the cluster" for the
#: headline count, and are masked out of the vertical density profile.
CLUSTER_MASK_RADIUS_PC = 10.0


@dataclass(frozen=True)
class FamousStar:
    """A star worth labelling on the map."""

    #: Label shown on the map.
    label: str
    #: Identifier SIMBAD resolves.
    simbad: str
    #: One line of context: why a reader would care about this particular star.
    note: str


#: Stars inside the map. Chosen because a reader will have heard of them, and
#: each makes a point: the nearest stars, the brightest stars, famous planetary
#: systems, and the objects that Gaia itself cannot see.
FAMOUS_STARS: list[FamousStar] = [
    FamousStar("Proxima Centauri", "Proxima Centauri", "The nearest star to the Sun; a red dwarf with a planet in its habitable zone"),
    FamousStar("Alpha Centauri A", "alpha Centauri A", "Brighter star of the pair Proxima orbits; a near-twin of the Sun"),
    FamousStar("Alpha Centauri B", "alpha Centauri B", "Second star of the Alpha Centauri pair"),
    FamousStar("Barnard's Star", "Barnard's star", "Fastest-moving star across our sky"),
    FamousStar("Luhman 16", "Luhman 16", "The nearest brown dwarfs: a pair of failed stars"),
    FamousStar("WISE 0855-0714", "WISE J085510.83-071442.5", "Coldest known brown dwarf, near the freezing point of water: invisible to Gaia"),
    FamousStar("Wolf 359", "Wolf 359", "The third-nearest star system: a small, flaring red dwarf"),
    FamousStar("Sirius", "Sirius", "The brightest star in the night sky"),
    FamousStar("Epsilon Eridani", "eps Eri", "A young orange dwarf with a debris disc and a giant planet"),
    FamousStar("Tau Ceti", "tau Cet", "The nearest single Sun-like star"),
    FamousStar("61 Cygni", "61 Cyg A", "The first star to have its distance measured (Bessel, 1838)"),
    FamousStar("Procyon", "Procyon", "Bright F star with a white dwarf companion"),
    FamousStar("Altair", "Altair", "A star spinning so fast it bulges at the equator"),
    FamousStar("Vega", "Vega", "The historical zero point of the magnitude scale"),
    FamousStar("Fomalhaut", "Fomalhaut", "Bright star ringed by a dusty debris belt"),
    FamousStar("TRAPPIST-1", "TRAPPIST-1", "Ultracool dwarf with seven Earth-sized planets"),
    FamousStar("Arcturus", "Arcturus", "The brightest star in the northern sky, an orange giant"),
    FamousStar("51 Pegasi", "51 Peg", "First planet found around a Sun-like star (1995)"),
    FamousStar("Pollux", "Pollux", "The nearest giant star"),
    FamousStar("Capella", "Capella", "A pair of giant stars orbiting each other"),
    FamousStar("Aldebaran", "Aldebaran", "Red giant that lies in front of the Hyades, not in it"),
    FamousStar("Regulus", "Regulus", "Hot blue B star spinning close to break-up; the brightest in Leo"),
    FamousStar("Achernar", "Achernar", "A B star spinning so fast it is visibly squashed"),
    FamousStar("Spica", "Spica", "A hot, massive binary; the brightest star in Virgo"),
    FamousStar("Canopus", "Canopus", "Second-brightest star in the sky, right at this map's edge"),
]

#: Famous stars beyond the map's edge. Shown as signposts on the boundary: their
#: absence from the map is itself the lesson, since most of the brightest stars
#: in our sky are rare, luminous giants that are visible from far away.
BEYOND_STARS: list[FamousStar] = [
    FamousStar("Polaris", "Polaris", "The North Star, a Cepheid pulsating supergiant"),
    FamousStar("Betelgeuse", "Betelgeuse", "Red supergiant; its distance is notoriously uncertain"),
    FamousStar("Antares", "Antares", "Red supergiant at the heart of Scorpius"),
    FamousStar("Rigel", "Rigel", "Blue supergiant in Orion"),
    FamousStar("Deneb", "Deneb", "One of the most luminous stars visible to the naked eye"),
    FamousStar("Pleiades", "M45", "The Seven Sisters, a young open cluster"),
    FamousStar("Kepler-10 (project A1)", "Kepler-10", "The rocky planet host from project A1"),
    FamousStar("HAT-P-7 (project A1)", "HAT-P-7", "The hot-Jupiter host from project A1"),
]

#: Star clusters whose literature distances we compare against, by SIMBAD name.
REFERENCE_CLUSTERS = {"Hyades": "Hyades", "Coma Berenices": "Melotte 111"}
