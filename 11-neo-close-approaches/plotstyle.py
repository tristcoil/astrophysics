"""Shared matplotlib styling for the project's figures.

Identical in look to projects A1 and A2, so the portfolio reads as one.

The categorical colours are a colourblind-safe set: adjacent slots keep a CVD
colour difference of dE >= 8 (OKLab x100), so the figures stay readable for
deuteranopia/protanopia and in greyscale print. That guarantee covers
*adjacent* slots -- stacked bars, lines. Scatter plots, where every pair of
colours can meet, use at most the first three slots, which were checked on all
pairs (``validate_palette.js --pairs all``: worst CVD dE 9.2).
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

# Semantic roles. Naming them once keeps the same quantity the same colour in
# every figure, and matches A1/A2: blue for the data a figure is about, orange
# for our own estimate, violet for a published result we compare against.
C_CONTEXT = "#d6d5ce"  # the whole population, drawn behind a highlighted subset
C_DATA = BLUE  # the approaches / objects a figure is about
C_FIT = ORANGE  # our own estimate
C_PUBLISHED = VIOLET  # a published result we compare against
C_EARTH = "#86b6ef"  # Earth drawn to scale (blue ramp step 250)

#: When an asteroid was first observed relative to its closest approach. The
#: first three categorical slots, validated on all pairs for scatter use; aqua
#: sits below 3:1 contrast, so every use carries a legend and direct labels.
LAG_COLORS = {"before": BLUE, "same_day": AQUA, "after": ORANGE}
LAG_LABELS = {"before": "first seen before closest approach",
              "same_day": "first seen on the day of closest approach",
              "after": "first seen only afterwards"}

#: Three-step ordinal blue ramp (steps 250 -> 450 -> 650), darker = larger.
SIZE_RAMP_3 = ["#86b6ef", "#2a78d6", "#104281"]

#: Torino scale colour code (white / green / yellow / orange / red), drawn as
#: washes of the reserved status colours; every region also carries its number.
TORINO_WASH = {"white": SURFACE, "green": "#0ca30c", "yellow": "#fab219",
               "orange": "#ec835a", "red": "#d03b3b"}


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


def annotate(ax, text: str, xy, xytext, *, color: str = INK_SECONDARY, relpos=(0.5, 0.5),
             fontsize: float = 8, **kwargs):
    """A callout with a thin leader line, in the project's annotation style.

    ``relpos`` is where on the text's bounding box the leader starts, in box
    fractions: (0, 0) is bottom-left, (1, 1) top-right. The default centre
    makes leaders from left-aligned text run diagonally, so set it to the
    corner or edge nearest the point being labelled.
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
