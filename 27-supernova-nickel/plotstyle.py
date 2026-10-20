"""Shared matplotlib styling for the project's figures.

The categorical colours are a colourblind-safe set: adjacent slots keep a CVD
colour difference of dE >= 8 (OKLab x100), so the figures stay readable for
deuteranopia/protanopia and in greyscale print. The four supernova colours
(blue, violet, orange, aqua) appear together as points, so they were checked
for *every* pair with the dataviz validator, not eyeballed: worst colour-blind
pair dE 9.2, worst normal-vision pair 16.3.

(Same palette and chrome as 01-exoplanet-transit and 05-supernova-cosmology.)
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

# Semantic roles. Each supernova keeps its colour in every figure; the four
# together were validated all-pairs (worst colour-blind dE 9.2, normal-vision
# 16.3). Aqua is below 3:1 contrast on the surface, so every supernova is also
# labelled directly (never colour alone).
SN_COLORS = {
    "SN2011fe": BLUE,
    "SN2014J": VIOLET,
    "SN1987A": ORANGE,
    "SN2023ixf": AQUA,
}
SN_MARKERS = {"SN2011fe": "o", "SN2014J": "s", "SN1987A": "D", "SN2023ixf": "^"}
C_MODEL = INK_SECONDARY  # reference physics curves (full trapping, decay power)
C_LIT = MUTED  # published values
C_REST = "#bdbcb5"  # de-emphasised points (e.g. dusty ZTF supernovae)

# Sequential blue ramp (for ordered quantities).
BLUE_150 = "#b7d3f6"
BLUE_300 = "#6da7ec"
BLUE_550 = "#1c5cab"

# Neutral wash for regions of a plot that are excluded or unobserved.
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
