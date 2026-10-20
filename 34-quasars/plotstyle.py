"""Shared matplotlib styling for the project's figures.

Copied from project J3 (24-galaxy-life-cycle), itself from D2 and A1. Categorical colours are a
colourblind-safe set, checked with a validator rather than by eye: the pairs
this project puts side by side (blue/red galaxies; our measurement/published)
keep a colour-vision-deficiency difference of dE >= 8 (OKLab x100) in every
pairing, so the figures stay readable for deuteranopia/protanopia.
The roles below are this project's.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

# --- Categorical slots ---
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
YELLOW = "#eda100"
MAGENTA = "#e87ba4"
GREEN = "#008300"
VIOLET = "#4a3aa7"
RED = "#e34948"

# --- Sequential ramp (one hue, light -> dark) for densities and maps ---
BLUE_RAMP = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
             "#2a78d6", "#1c5cab", "#104281", "#0d366b"]

# --- Chrome and ink ---
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

# Semantic roles (validated with the dataviz palette script, light surface):
# {blue, orange, aqua} pass all-pairs (CVD dE >= 9.2), {blue, violet, orange}
# pass all-pairs (>= 13.0). Aqua is below 3:1 contrast, so aqua marks are always
# direct-labelled.
C_OURS = BLUE             # our own measurement
C_PUBLISHED = VIOLET      # the literature value we compare against
C_MODEL = ORANGE          # fitted model curves
C_ANNOT = INK_SECONDARY

#: one colour per emission line, the same in every figure
LINE_COLOURS = {"Hbeta": AQUA, "MgII": BLUE, "CIV": ORANGE}
LINE_LABELS = {"Hbeta": "H\u03b2", "MgII": "Mg II", "CIV": "C IV"}

#: survey programmes in DR16Q, fixed order; the heterogeneous rest is grey
PROGRAMME_COLOURS = {"eBOSS CORE": BLUE, "BOSS CORE": ORANGE, "SDSS-I/II": AQUA,
                     "Other BOSS/eBOSS": "#b4b2aa"}

#: ordinal ramp (validated --ordinal) for ordered classes such as luminosity bins
ORDINAL_BLUE = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]


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
    kwargs.setdefault("fontsize", 8)
    return ax.annotate(
        text,
        xy=xy,
        xytext=xytext,
        color=color,
        arrowprops={
            "arrowstyle": "-",
            "color": color,
            "linewidth": 0.8,
            "shrinkA": 2,
            "shrinkB": 3,
        },
        **kwargs,
    )


def caption(fig, text: str, y: float = -0.02) -> None:
    """Footnote under a figure: data provenance, caveats, the honest small print.

    Wrapped explicitly to the figure's width (~15 characters per inch at this
    size) rather than by matplotlib, whose wrapping can overrun the edges
    once a tight bounding box is applied.
    """
    import textwrap

    width = int(fig.get_figwidth() * 15.5)
    fig.text(0.5, y, textwrap.fill(text, width), ha="center", va="top", fontsize=7.5,
             color=MUTED, linespacing=1.35)


def save(fig, path_stem, *, formats=("png", "svg"), dpi: int | None = None) -> list:
    """Save a figure to every requested format, returning the paths written.

    Artists created with ``rasterized=True`` (the galaxy clouds) are embedded
    in the SVG as images at ``dpi``; a vector file with 50,000 individual dots
    would run to tens of megabytes and render slowly.
    """
    written = []
    for fmt in formats:
        out = f"{path_stem}.{fmt}"
        fig.savefig(out, format=fmt, dpi=dpi or mpl.rcParams["savefig.dpi"])
        written.append(out)
    plt.close(fig)
    return written
