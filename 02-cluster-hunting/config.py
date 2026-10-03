"""Project configuration: paths, the clusters we hunt, and every tunable number.

The clusters are a difficulty ramp, as in A1. The first is the textbook case:
nearby, rich, and moving across the sky so differently from the background that
it cannot be missed. The last is six times further away, and moves almost like
the field stars around it. Running one identical pipeline over all three is the
point: it shows the method doing the work, not the choice of target.

The three also span a factor of ~35 in age, which is what makes the final
comparison figure work: the same kind of star cluster, photographed at three
different moments of its life.

What the configuration deliberately does NOT contain is any cluster's proper
motion or distance. We tell the pipeline where on the sky to look -- all three
are visible to the naked eye or binoculars, so that is not a secret -- and it has
to find the cluster's motion and distance in the data itself. Published values
are fetched separately and used only at the very end, to grade the result.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
GAIA_DIR = DATA_DIR / "gaia"
LITERATURE_DIR = DATA_DIR / "literature"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, GAIA_DIR, LITERATURE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)


@dataclass(frozen=True)
class ReferenceAge:
    """A published age, used only to annotate figures and grade the story.

    Cluster ages are model-dependent and the literature disagrees at the 10-30%
    level, so each carries its spread across methods rather than one number.
    """

    myr: float
    low_myr: float
    high_myr: float
    method: str
    citation: str


@dataclass(frozen=True)
class Cluster:
    """One star cluster to hunt for."""

    #: Common name, used in titles.
    name: str
    #: Short lowercase key, used for file names.
    key: str
    #: Name in the Hunt & Reffert (2023) Gaia DR3 cluster catalogue.
    catalogue_name: str
    #: Name in the Gaia Collaboration (2018) DR2 cluster tables.
    gc2018_name: str
    #: Centre of the cone we download. This is where we point the telescope,
    #: not a kinematic input: it tells the pipeline nothing about motion.
    ra_deg: float
    dec_deg: float
    #: Cone radius. Chosen to enclose the cluster's tidal radius -- the distance
    #: beyond which the Galaxy's pull beats the cluster's own -- with some margin.
    radius_deg: float
    #: One-line reason this cluster earns its place in the ramp.
    why: str
    difficulty: str
    #: Faint limit. Below G ~ 18 Gaia's colours and astrometry degrade quickly.
    g_max: float = 18.0
    #: Number of Gaussians used to describe the field-star population in the
    #: membership model. The field is lumpy, so it needs several.
    n_field_components: int = 6
    #: Named stars to label, resolved to Gaia DR3 source IDs through SIMBAD.
    famous_stars: tuple[str, ...] = ()
    reference_age: ReferenceAge | None = None


#: The clusters, in the order the pipeline runs them.
CLUSTERS: list[Cluster] = [
    Cluster(
        name="Pleiades",
        key="pleiades",
        catalogue_name="Melotte_22",
        gc2018_name="Pleiades",
        ra_deg=56.75,
        dec_deg=24.12,
        radius_deg=5.0,
        why=(
            "The Seven Sisters. Close (~135 pc), rich (>1000 members), and moving "
            "at ~50 mas/yr where the field stars around it move at ~5, so its "
            "members stand apart from everything else in proper-motion space."
        ),
        difficulty="easy",
        famous_stars=(
            "Alcyone", "Atlas", "Electra", "Maia", "Merope",
            "Taygeta", "Pleione", "Celaeno", "Asterope",
        ),
        reference_age=ReferenceAge(
            myr=112, low_myr=107, high_myr=117,
            method="lithium depletion boundary",
            citation="Dahm (2015), ApJ 813, 108",
        ),
    ),
    Cluster(
        name="Praesepe",
        key="praesepe",
        catalogue_name="NGC_2632",
        gc2018_name="Praesepe",
        ra_deg=130.10,
        dec_deg=19.67,
        radius_deg=4.0,
        why=(
            "The Beehive. At a similar distance to the Pleiades but receding at "
            "~34 km/s, so perspective makes its members' motions converge "
            "noticeably across the field -- a real test of the geometry."
        ),
        difficulty="medium",
        reference_age=ReferenceAge(
            myr=670, low_myr=590, high_myr=800,
            method="gyrochronology; isochrone fits span 590-800 Myr",
            citation=(
                "Douglas et al. (2019), ApJ 879, 100; "
                "Gossage et al. (2018); Brandt & Huang (2015)"
            ),
        ),
    ),
    Cluster(
        name="M67",
        key="m67",
        catalogue_name="NGC_2682",
        gc2018_name="NGC2682",
        ra_deg=132.85,
        dec_deg=11.82,
        radius_deg=1.0,
        why=(
            "One of the oldest open clusters known, about the Sun's age and "
            "composition. Six times further than the Pleiades, so its motion "
            "is small and its parallax sits where the field stars are densest."
        ),
        difficulty="hard",
        reference_age=ReferenceAge(
            myr=4000, low_myr=3460, high_myr=4200,
            method="asteroseismology (3.46 Gyr) and gyrochronology (4.2 Gyr)",
            citation="Stello et al. (2016), ApJ 832, 133; Barnes et al. (2016), ApJ 823, 16",
        ),
    ),
]


# --------------------------------------------------------------------------
# Data quality
# --------------------------------------------------------------------------

#: Renormalised unit weight error: how well a single-star model fits the
#: astrometry. 1.0 is a perfect fit; above ~1.4 the source is probably an
#: unresolved binary or otherwise problematic (Lindegren et al. 2021).
RUWE_GOOD = 1.4

#: Photometric consistency filter: keep sources whose corrected BP/RP flux
#: excess C* lies within this many sigma of zero (Riello et al. 2021). Catches
#: blends and nebular contamination -- which the Pleiades, sitting in a
#: reflection nebula, genuinely has.
CSTAR_NSIGMA = 3.0

#: Riello et al. (2021) advise not applying the C* filter brighter than this,
#: where detector saturation dominates the excess.
CSTAR_G_MIN = 4.0

#: Field stars shown in the comparison HR diagram need distances good to 10%,
#: otherwise the diagram is smeared by parallax noise rather than astrophysics.
FIELD_PARALLAX_OVER_ERROR = 10.0

#: The *measured* width of the field main sequence uses a stricter 5% cut, so
#: that parallax noise (at most 0.11 mag here) cannot pass for astrophysics.
WIDTH_PARALLAX_OVER_ERROR = 20.0

#: Colour range for measuring main-sequence width: G, K and early-M dwarfs,
#: where every cluster here has a well-populated sequence.
WIDTH_BP_RP_RANGE = (0.8, 2.2)

#: Stars within this many magnitudes of the cluster's ridge line count towards
#: the width, which keeps giants (2+ mag brighter) out of a dwarf statistic.
WIDTH_WINDOW_MAG = 1.5

#: Gaia's parallaxes share systematic errors of order 0.01 mas across patches
#: of sky a few degrees wide (Lindegren et al. 2021; Maiz Apellaniz et al.
#: 2021). No number of stars averages this away, so it floors the precision of
#: any cluster distance.
PARALLAX_SYSTEMATIC_MAS = 0.010

# --------------------------------------------------------------------------
# Finding the knot in proper-motion space
# --------------------------------------------------------------------------

#: Histogram bin for the proper-motion density map, mas/yr.
KNOT_BIN_MASYR = 0.25
#: The knot is found like a star in an image: count stars inside a small
#: circle ("aperture") at every point of proper-motion space, and compare with
#: the density in a ring around it. Radii in mas/yr.
KNOT_APERTURE_MASYR = 1.5
KNOT_ANNULUS_MASYR = (3.0, 6.0)
#: A candidate needs at least this many stars in its aperture. Stops a chance
#: handful of stars in an empty corner from registering as a huge contrast.
KNOT_MIN_STARS = 25
#: How far from the field's median motion to search, mas/yr.
KNOT_SEARCH_HALFWIDTH_MASYR = 100.0
#: How many of the highest-contrast places to examine.
KNOT_CANDIDATES = 5
#: A candidate is a cluster only if its stars share a distance: the robust
#: spread of their parallaxes must be below this fraction of the median. The
#: Pleiades knot comes in at ~3%; the field's own core at ~100%.
KNOT_MAX_PARALLAX_SPREAD = 0.25

# --------------------------------------------------------------------------
# Membership model
# --------------------------------------------------------------------------

#: A star is called a member when the model gives it at least this probability.
MEMBERSHIP_THRESHOLD = 0.5

#: Stars within this many mas/yr of the knot seed the cluster component.
SEED_RADIUS_MASYR = 3.0

#: The model is fitted to stars within this distance of the knot in
#: proper-motion space. Membership is a local question: a star moving 30 mas/yr
#: differently from a cluster whose members agree to ~1 mas/yr is not a
#: candidate, so the field only needs describing nearby. This is 15x the
#: Pleiades' internal spread; stars outside it get membership probability 0.
FIT_REGION_MASYR = 15.0

#: A fitted Gaussian counts as part of the cluster's halo only if its parallax
#: spread is below this fraction of the cluster's parallax (+/-15% in
#: distance). Field components in the same patch of motion space span 25-50%.
HALO_MAX_DEPTH_FRACTION = 0.15

#: Extreme-deconvolution iteration controls.
XD_MAX_ITER = 400
XD_TOL = 1e-8
#: Small variance floor added to each Gaussian so none can collapse onto a
#: single star, in (mas/yr)^2 or mas^2.
XD_REGULARISATION = 1e-4

#: Rounds of "fit, re-derive the cluster's space velocity, re-correct the
#: proper motions for perspective". The geometry converges in two.
PERSPECTIVE_ROUNDS = 2

# --------------------------------------------------------------------------
# Reference points
# --------------------------------------------------------------------------

#: The Sun in Gaia's bands (Casagrande & VandenBerg 2018, MNRAS 479, L102).
SUN_ABS_G = 4.67
SUN_BP_RP = 0.82


def cluster_by_key(key: str) -> Cluster:
    for c in CLUSTERS:
        if key.lower() in (c.key, c.name.lower()):
            return c
    raise KeyError(f"unknown cluster {key!r}; known: {[c.key for c in CLUSTERS]}")
