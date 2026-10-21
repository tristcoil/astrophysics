"""Project configuration: paths, data sources, assumptions, and the showcase galaxies.

Every number the analysis *assumes* rather than measures lives in this file, so
the assumptions can be read in one place and changed in one place. There are
fewer of them than you might expect: the whole visible-matter prediction for a
galaxy has one free parameter, the stellar mass-to-light ratio, and we fix it
from stellar-population physics rather than fitting it to the rotation curve.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------
# Data sources
# --------------------------------------------------------------------------

#: The SPARC database at Case Western Reserve University. HTTPS only -- the
#: server refuses plain-HTTP connections on port 80.
SPARC_URL = "https://astroweb.cwru.edu/SPARC/"

SPARC_PAPER = "Lelli, McGaugh & Schombert 2016, AJ 152, 157"


@dataclass(frozen=True)
class Source:
    """One file we download, and why."""

    #: File name under ``data/``.
    filename: str
    url: str
    #: What it holds and what we use it for.
    what: str
    #: The paper to cite for it.
    reference: str
    #: ``True`` if the analysis is built on it; ``False`` if it is a published
    #: result used only to grade our own at the end.
    input: bool = True


SOURCES: list[Source] = [
    Source(
        "SPARC_Lelli2016c.mrt", SPARC_URL + "SPARC_Lelli2016c.mrt",
        "Table 1: 175 galaxies -- distance, inclination, 3.6um luminosity, HI mass, "
        "flat rotation speed, quality flag",
        SPARC_PAPER,
    ),
    Source(
        "MassModels_Lelli2016c.mrt", SPARC_URL + "MassModels_Lelli2016c.mrt",
        "Table 2: 3391 rotation-curve points, each with the gas, stellar-disc and bulge "
        "contributions computed from the observed light and 21-cm gas",
        SPARC_PAPER,
    ),
    Source(
        "Bulges.mrt", SPARC_URL + "Bulges.mrt",
        "Bulge luminosities at 3.6um, needed to add up each galaxy's total visible mass",
        SPARC_PAPER,
    ),
    Source(
        "RAR.mrt", SPARC_URL + "RAR.mrt",
        "The published radial acceleration relation: 2693 (g_bar, g_obs) points",
        "Lelli et al. 2017, ApJ 836, 152", input=False,
    ),
    Source(
        "RARbins.mrt", SPARC_URL + "RARbins.mrt",
        "The published radial acceleration relation, binned",
        "Lelli et al. 2017, ApJ 836, 152", input=False,
    ),
    Source(
        "BTFR_Lelli2019.mrt", SPARC_URL + "BTFR_Lelli2019.mrt",
        "Published baryonic masses and flat velocities (baryonic Tully-Fisher relation)",
        "Lelli et al. 2019, MNRAS 484, 3267", input=False,
    ),
    Source(
        "MaximumDiskData.mrt", SPARC_URL + "MaximumDiskData.mrt",
        "Published maximum-disc mass-to-light ratios",
        "Starkman et al. 2018, MNRAS 480, 2292", input=False,
    ),
    Source(
        "planetary_factsheet.html", "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
        "NASA Planetary Fact Sheet: planet distances, eccentricities and mean orbital speeds",
        "D. R. Williams, NASA Space Science Data Coordinated Archive",
    ),
]


# --------------------------------------------------------------------------
# Assumptions
# --------------------------------------------------------------------------

#: Stellar mass-to-light ratios at 3.6 microns, in Msun/Lsun. These are the
#: SPARC defaults, from stellar population synthesis models (McGaugh &
#: Schombert 2014; Schombert & McGaugh 2014). The near-infrared is used
#: precisely because the mass-to-light ratio there barely depends on a galaxy's
#: age or star formation history -- in optical light it varies several-fold.
#: Bulges are older and redder, so get a higher value.
UPSILON_DISK = 0.5
UPSILON_BULGE = 0.7

#: Galaxy-to-galaxy scatter in the true mass-to-light ratio, in dex (McGaugh,
#: Lelli & Schombert 2016). Used only in error bars.
UPSILON_SCATTER_DEX = 0.11

#: The mass-to-light ratios stellar population models consider plausible at
#: 3.6um. Anything far outside this band is not a population of stars anyone
#: knows how to make. Used to judge the maximum-disc test.
UPSILON_PLAUSIBLE = (0.2, 0.8)

#: Hubble constant for halo virial quantities (R200, M200), km/s/Mpc. Matches
#: the value SPARC assumed for its Hubble-flow distances, so the two are
#: consistent. Only the halo *sizes and masses* depend on it, not the fits.
H0 = 73.0

#: Fractional uncertainty on 21-cm HI masses from flux calibration.
HI_MASS_ERROR = 0.10


# --------------------------------------------------------------------------
# Sample selection for the population-level results
# --------------------------------------------------------------------------

#: The quality cuts of McGaugh, Lelli & Schombert (2016), which the radial
#: acceleration relation was published with. Using the same cuts is what makes
#: our result directly comparable to theirs.
#:
#: * Quality flag 3 marks 12 galaxies whose rotation curves are too asymmetric
#:   or disturbed to trace an equilibrium orbit.
#: * Below 30 degrees inclination a disc is so close to face-on that the
#:   correction from line-of-sight speed to rotation speed (divide by sin i)
#:   blows up small inclination errors.
#: * Points whose velocity is uncertain by more than 10% are dropped.
MAX_QUALITY_FLAG = 2
MIN_INCLINATION_DEG = 30.0
MAX_FRACTIONAL_VELOCITY_ERROR = 0.10

#: For the maximum-disc test, how many standard errors of slack the stars get:
#: they may overshoot the observed inner curve by this much, and need only reach
#: the observed outer speed minus this much. Two sigma is deliberately generous
#: to the no-dark-matter hypothesis -- the innermost points of 21-cm rotation
#: curves carry the largest systematic errors (the radio beam smears the steep
#: inner rise), and one of them should not be what rules a heavy disc out. At
#: this setting our maximum-disc mass-to-light ratios sit ~5% *above* the
#: published ones (Starkman et al. 2018); at 1 sigma they sit ~6% below.
MAX_DISK_TOLERANCE_SIGMA = 2.0


# --------------------------------------------------------------------------
# The showcase galaxies
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Showcase:
    """A galaxy featured individually in the figures."""

    #: SPARC identifier, e.g. "NGC3198".
    name: str
    #: What kind of galaxy it is, in a few words.
    kind: str
    #: Why it earns its place.
    why: str


#: Ordered from the most massive to the least. Chosen to span what SPARC spans:
#: giant to dwarf, bulge-dominated to bulgeless, star-dominated to
#: gas-dominated, and high to low surface brightness. All have quality flag 1 or
#: 2 and inclinations well above 30 degrees.
SHOWCASE: list[Showcase] = [
    Showcase(
        "UGC02885", "giant spiral",
        "One of the largest spiral galaxies known. Vera Rubin measured its rotation "
        "in 1980; Benne Holwerda, who imaged it with Hubble in 2020, nicknamed it "
        "Rubin's Galaxy in her honour.",
    ),
    Showcase(
        "NGC2841", "bulge-dominated spiral",
        "A massive spiral with a large central bulge -- the gallery's case where the "
        "bulge, not the disc, dominates the inner rotation.",
    ),
    Showcase(
        "NGC3198", "the classic",
        "Subject of van Albada, Bahcall, Begeman & Sancisi (1985), whose split of its "
        "21-cm rotation curve into disc and dark halo became the textbook picture.",
    ),
    Showcase(
        "NGC2403", "compact spiral",
        "Rotates at the same speed as UGC 128, but packs its stars 16 times more "
        "densely -- half of a telling pair.",
    ),
    Showcase(
        "UGC00128", "diffuse spiral",
        "NGC 2403's twin in rotation speed, with its stars spread 16 times more "
        "thinly. The other half of the pair: dark matter dominates it everywhere.",
    ),
    Showcase(
        "DDO154", "gas-rich dwarf",
        "A dwarf irregular with over ten times more mass in gas than in stars. "
        "Carignan & Freeman (1988) titled their paper on it \"DDO 154 -- a 'dark' "
        "galaxy?\".",
    ),
]

#: The galaxy the flagship figure is built around.
FLAGSHIP = "NGC3198"


def showcase_by_name(name: str) -> Showcase:
    for s in SHOWCASE:
        if s.name.lower() == name.lower():
            return s
    raise KeyError(f"{name!r} is not a showcase galaxy; known: {[s.name for s in SHOWCASE]}")
