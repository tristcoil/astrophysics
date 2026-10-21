"""Project configuration: paths and the target list.

Targets are chosen as a difficulty ramp. The first is unmistakable -- you can see
its transits by eye in the raw light curve. The last is a genuinely hard
detection that no amount of squinting will reveal, and which only emerges once
many transits are folded together. Running the identical pipeline over both is
the point: it shows that the method, not the eye, is doing the work.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
CACHE_DIR = DATA_DIR / "mast"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
STELLAR_PARAMS_PATH = DATA_DIR / "stellar_params.json"

for _d in (DATA_DIR, CACHE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)


@dataclass(frozen=True)
class Target:
    """One star to search for transits around."""

    #: Name resolvable by MAST and by the NASA Exoplanet Archive.
    name: str
    #: Kepler Input Catalog identifier, for provenance.
    kic: str
    #: Name of the planet whose published parameters we compare against.
    planet: str
    #: Short label used in figure titles.
    label: str
    #: One-line reason this target earns its place in the ramp.
    why: str
    #: Difficulty of the detection: how obvious the transit is.
    difficulty: str
    #: Quarters to download. ``None`` means every available quarter.
    quarters: tuple[int, ...] | None = None
    #: Savitzky-Golay detrending window in days. Must comfortably exceed the
    #: transit duration or the filter will absorb the transit itself.
    flatten_window_days: float = 0.75
    #: Period search bounds in days for the BLS periodogram.
    period_min_days: float = 0.5
    period_max_days: float = 15.0
    #: Trial transit durations in days handed to BLS.
    durations_days: tuple[float, ...] = (0.05, 0.075, 0.1, 0.15, 0.2)
    #: Half-width of the phase window plotted around mid-transit, in days.
    zoom_halfwidth_days: float = 0.35
    notes: dict = field(default_factory=dict)


#: The targets, in the order the pipeline runs them.
TARGETS: list[Target] = [
    Target(
        name="HAT-P-7",
        kic="KIC 10666592",
        planet="HAT-P-7 b",
        label="HAT-P-7 b",
        why=(
            "A hot Jupiter with a ~0.6% transit -- deep enough to see by eye in the "
            "raw photometry, which makes it the ideal case for checking the pipeline "
            "against something unambiguous."
        ),
        difficulty="easy",
        flatten_window_days=0.75,
        period_min_days=0.4,
        period_max_days=6.0,
        # Trial durations bracket the ~3.9 h transit without assuming it.
        durations_days=(0.08, 0.12, 0.16, 0.20, 0.26),
        zoom_halfwidth_days=0.35,
    ),
    Target(
        name="Kepler-10",
        kic="KIC 11904151",
        planet="Kepler-10 b",
        label="Kepler-10 b",
        why=(
            "The first unambiguously rocky exoplanet. Its transit is only ~190 ppm "
            "-- 0.019%, invisible in any single transit -- so it is recoverable "
            "only by folding well over a thousand orbits together."
        ),
        difficulty="hard",
        # Measured depth converges above ~1 day; shorter windows let the filter
        # nibble at the transit and cost a few percent of depth.
        flatten_window_days=1.0,
        period_min_days=0.4,
        period_max_days=6.0,
        # Trial durations bracket the ~1.8 h transit without assuming it.
        durations_days=(0.03, 0.05, 0.07, 0.09),
        zoom_halfwidth_days=0.15,
    ),
]

#: Bond albedo assumed when quoting equilibrium temperatures. Stated explicitly
#: because T_eq is meaningless without it.
ASSUMED_ALBEDO = 0.3

#: Quadratic limb-darkening coefficients used only for the first-order depth
#: correction, typical of a Sun-like star in the broad Kepler bandpass.
LIMB_DARKENING_U1 = 0.4
LIMB_DARKENING_U2 = 0.26


def target_by_name(name: str) -> Target:
    for t in TARGETS:
        if t.name.lower() == name.lower():
            return t
    raise KeyError(f"unknown target {name!r}; known: {[t.name for t in TARGETS]}")
