"""Shared matplotlib styling for the project's figures.

The palette is the atlas-wide one (D1, D8): a colourblind-safe categorical set,
checked with a validator rather than by eye. Each colour has one job across
every figure:

* the three isochrone sets: MIST blue, BaSTI orange, DSED aqua (slots 1-3,
  which pass the colour-blindness checks for all pairs at once);
* the age of the universe from earlier atlas projects: D3 (CMB) orange as in
  D8's figures, D1 (supernovae) violet, D8 (chronometers) aqua -- these never
  share a panel with the model-set colours;
* first- and second-generation stars in the chromosome maps: aqua and magenta,
  Milone et al.'s own convention. That pair is in the 6-8 dE colour-blind floor
  band, so the two groups also differ in marker shape and are labelled directly.

Stars in colour-magnitude diagrams are ink, as a density; ages that are ordered
(10, 12, 14 Gyr ...) take one hue stepped light to dark.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

# --- Categorical slots, assigned in fixed order (never cycled past slot 8) ---
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
YELLOW = "#eda100"
MAGENTA = "#e87ba4"
GREEN = "#008300"
VIOLET = "#4a3aa7"
RED = "#e34948"

SERIES = [BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED]

# --- Chrome and ink ---
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

# Semantic roles used across the figures.
C_STARS = INK
C_OURS = BLUE  # our combined cluster ages
MODEL_COLOURS = {"MIST": BLUE, "BaSTI": ORANGE, "DSED": AQUA}
C_SN = VIOLET  # D1: supernovae on the Cepheid distance ladder
C_CMB = ORANGE  # D3: the cosmic microwave background (Planck), assuming LCDM
C_CC = AQUA  # D8: cosmic chronometers
C_PUBLISHED = MUTED  # published numbers we check against
C_1G = AQUA
C_2G = MAGENTA
#: Literature markers: shape, not colour, tells the sources apart.
LIT_MARKERS = {"Ying et al. 2023": "D", "VandenBerg et al. 2013": "s", "Valcin et al. 2025": "^",
               "white-dwarf cooling (Hansen+07)": "v", "white-dwarf cooling (Hansen+13)": "v",
               "eclipsing binaries (Thompson+20)": "P"}

# Sequential blue ramp (ordinal steps for ordered ages; lightest still >= 2:1).
BLUE_150 = "#b7d3f6"
BLUE_250 = "#86b6ef"
BLUE_300 = "#6da7ec"
BLUE_450 = "#2a78d6"
BLUE_550 = "#1c5cab"
BLUE_650 = "#104281"
AGE_RAMP = [BLUE_250, BLUE_300, BLUE_450, BLUE_550, BLUE_650]
CONTOUR_FILLS = [BLUE_150, BLUE_300, BLUE]
# Diverging pair for signed quantities (differential reddening): blue <-> red, gray middle.
DIVERGING = ["#1c5cab", "#f0efec", "#b83232"]  # blue pole, neutral midpoint, red pole of similar lightness

# Neutral wash for regions of a plot that are excluded on physical grounds.
EXCLUDED = "#f0efec"


def use_style() -> None:
    """Apply the project-wide rcParams. Call once at import time in figure code."""
    mpl.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "savefig.bbox": "tight",
            "savefig.dpi": 200,
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans"],
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.titleweight": "bold",
            "axes.titlelocation": "left",
            "axes.titlepad": 10,
            "axes.labelsize": 9,
            "axes.labelcolor": INK_SECONDARY,
            "axes.edgecolor": AXIS,
            "axes.linewidth": 0.8,
            "axes.grid": True,
            "axes.axisbelow": True,
            # Only the top/right spines go; a chart still needs its two rules.
            "axes.spines.top": False,
            "axes.spines.right": False,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "grid.linestyle": "-",  # solid hairlines; dashed grids read as thresholds
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "legend.frameon": False,
            "legend.fontsize": 8,
            "lines.linewidth": 1.6,
            "lines.markersize": 4,
            "text.color": INK,
            "figure.titlesize": 12,
            "figure.titleweight": "bold",
        }
    )


def annotate(ax, text: str, xy, xytext, *, color: str = INK_SECONDARY, **kwargs):
    """A callout with a thin leader line, in the project's annotation style."""
    return ax.annotate(
        text,
        xy=xy,
        xytext=xytext,
        color=color,
        fontsize=8,
        arrowprops={
            "arrowstyle": "-",
            "color": color,
            "linewidth": 0.8,
            "shrinkA": 2,
            "shrinkB": 4,
        },
        **kwargs,
    )


def caption(fig, text: str, y: float = -0.02) -> None:
    """Footnote under a figure: data provenance, caveats, the honest small print."""
    fig.text(
        0.5,
        y,
        text,
        ha="center",
        va="top",
        fontsize=7.5,
        color=MUTED,
        wrap=True,
    )


def save(fig, path_stem, *, formats=("png", "svg")) -> list:
    """Save a figure to every requested format, returning the paths written."""
    written = []
    for fmt in formats:
        out = f"{path_stem}.{fmt}"
        fig.savefig(out, format=fmt)
        written.append(out)
    plt.close(fig)
    return written
