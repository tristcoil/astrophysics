"""Shared matplotlib styling for the project's figures.

Identical in look to projects A1 and A2, so the portfolio reads as one.

The categorical colours are a colourblind-safe set: adjacent slots keep a CVD
colour difference of dE >= 8 (OKLab x100). Scatter plots are stricter, since
every pair of colours can meet there; this project never puts more than four
hues in one scatter, and the four it uses together -- blue, orange, aqua,
violet -- pass the all-pairs check (worst CVD dE 9.2, normal-vision dE 16.3,
run with the dataviz validator). Aqua is below 3:1 contrast on the light
surface, so wherever it appears it is also named by a direct label.

Densities -- a million asteroids do not fit in a scatter plot -- use one
sequential hue, blue, light to dark, with the empty surface meaning zero.
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
# every figure, and matches A1 and A2: orange for our own measurement, violet
# for the value it is compared against.
C_CONTEXT = "#d6d5ce"  # the whole population, drawn behind a highlighted subset
C_DATA = BLUE  # the asteroids a figure is about
C_MEASURED = ORANGE  # our own measurement (a gap centre, a fitted ratio)
C_PREDICTED = VIOLET  # the prediction or published value it is compared to
C_PLANET = INK  # the planets, for scale

#: Populations, for the maps. First three slots only: they validate all-pairs.
C_BELT = BLUE
C_HILDA = ORANGE
C_TROJAN = AQUA

#: Sequential blue ramp, steps 100 -> 700 of the reference palette. Density
#: maps prepend the surface colour so that "no asteroids" is simply blank.
BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
             "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281",
             "#0d366b"]


def density_cmap():
    """Surface -> light blue -> dark blue, for log-scaled count maps."""
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("belt_density", BLUE_RAMP, N=256)
    cmap.set_under(SURFACE)
    cmap.set_bad(SURFACE)
    return cmap


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
