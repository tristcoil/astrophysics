"""Shared matplotlib styling for the project's figures.

Identical in look to projects A1 and D4, so the projects read as one
portfolio, plus the colour maps this project's images need.

The categorical colours are a colourblind-safe set: adjacent slots keep a CVD
colour difference of dE >= 8 (OKLab x100), so the figures stay readable for
deuteranopia/protanopia and in greyscale print. That guarantee covers lines;
scatter plots, where any two colours can meet, use at most three hues, and only
triples that pass the stricter all-pairs check (blue/orange/aqua,
blue/orange/violet). Aqua, yellow and magenta sit below 3:1 contrast on the
page, so every line drawn in them also carries a direct text label.
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
# every panel of every figure: blue for the northern (first-quadrant) half of
# the Galaxy, orange for the southern (fourth-quadrant) half, violet for
# published results we compare against, aqua for the masers.
C_OBSERVED = INK
C_OBSERVED_ERR = MUTED
C_NORTH = BLUE      # first quadrant, 0 < l < 90
C_SOUTH = ORANGE    # fourth quadrant, 270 < l < 360
C_MASER = AQUA      # Reid et al. 2019 parallax masers (our computation)
C_PUBLISHED = VIOLET
C_MODEL = MUTED     # a reference model curve, recessive
C_CONTEXT = "#d6d5ce"  # the whole population, drawn behind a highlighted subset
C_ARM = ORANGE      # spiral-arm tracks drawn over images

# --- Colour maps for images (sequential = one hue, light to dark; diverging =
# two opposite hues through a neutral grey) ---
_BLUE_RAMP = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#2a78d6",
              "#1c5cab", "#104281", "#0d366b", "#07203f"]
GREY_RAMP = ["#fcfcfb", "#e1e0d9", "#c3c2b7", "#a3a29a", "#898781", "#6d6c67",
             "#52514e", "#353432", "#1d1d1b", "#0b0b0b"]


def _cmap(name, colors):
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(name, colors, N=256)


CMAP_BLUE = _cmap("h2_blue", _BLUE_RAMP)
CMAP_GREY = _cmap("h2_grey", GREY_RAMP)
CMAP_DIVERGING = _cmap("h2_div", ["#0d366b", "#2a78d6", "#9ec5f4", "#f0efec",
                                  "#f3a9a8", "#e34948", "#8a1c1b"])


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
