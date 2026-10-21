"""Shared matplotlib styling for the project's figures.

The categorical colours are the colour-blind-safe set used across the atlas:
adjacent slots keep a CVD colour difference of dE >= 8 (OKLab x100), so the
figures stay readable for deuteranopia/protanopia and in greyscale print. The
four "environment" colours below were validated in this order (orange, blue,
yellow, aqua) in both light and dark mode; yellow and aqua sit below 3:1
contrast on the light surface, so every chart that uses them labels its series
directly.
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

# --- Sequential ramp (blue), light -> dark, for densities ---
BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
             "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]

# --- Chrome and ink ---
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

# --- Map furniture ---
LAND = "#ebeae4"  # one step off the surface: land recedes, data sits on top
COAST = "#c3c2b7"
OCEAN = SURFACE

# Semantic roles. How a meteorite reached the catalogue is the organising idea
# of the project, so each route keeps one colour in every figure.
C_FALL = ORANGE  # seen to fall
C_ANTARCTIC = BLUE  # collected on Antarctic ice
C_DESERT = YELLOW  # collected in a hot-desert collection area
C_CHANCE = AQUA  # found on its own, anywhere else
ENV_COLORS = {"fall": C_FALL, "antarctic": C_ANTARCTIC, "desert": C_DESERT, "chance": C_CHANCE}

C_PUBLISHED = VIOLET  # a literature value we compare against
C_FIT = INK_SECONDARY  # a fitted model line

#: Meteorite families, in slot order (a stacked bar puts neighbours side by side,
#: and slot order is the order that was validated for adjacency).
FAMILY_COLORS = {
    "ordinary": BLUE, "carbonaceous": ORANGE, "other_chondrite": AQUA,
    "achondrite": YELLOW, "stony_iron": MAGENTA, "iron": GREEN,
}


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


def annotate(ax, text: str, xy, xytext, *, color: str = INK_SECONDARY, fontsize: float = 8,
             relpos=(0.5, 0.5), **kwargs):
    """A callout with a thin leader line, in the project's annotation style.

    ``relpos`` is where on the label the leader line starts, in fractions of the
    label's box -- (0, 0) is its lower-left corner.
    """
    return ax.annotate(
        text,
        xy=xy,
        xytext=xytext,
        color=color,
        fontsize=fontsize,
        arrowprops={
            "arrowstyle": "-",
            "color": color,
            "linewidth": 0.8,
            "shrinkA": 2,
            "shrinkB": 4,
            "relpos": relpos,
        },
        **kwargs,
    )


def caption(fig, text: str) -> None:
    """Footnote under a figure: data provenance, caveats, the honest small print."""
    fig.text(
        0.5,
        -0.02,
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
