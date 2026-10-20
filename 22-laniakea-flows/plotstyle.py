"""Shared matplotlib styling for the project's figures (copied from G1 via I18).

Identical in look to the other projects', so the atlas reads as one portfolio.

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
# every panel of every figure, and matches the other projects where the roles
# overlap: orange for our own measurements, violet for the published values we
# compare against.
C_OBSERVED = INK  # what was measured: the loudest thing on every plot
C_OBSERVED_ERR = MUTED
C_OURS = ORANGE  # our own measurement
C_PUBLISHED = VIOLET  # a published result we compare against
C_CONTEXT = "#d6d5ce"  # background population, drawn behind a highlighted subset

# Flows: receding red, approaching blue, no motion a neutral grey. Both arms
# are single-hue ramps with monotone lightness (validated with the dataviz
# skill's ordinal check; the pale ends deliberately recede into the grey
# midpoint, as a continuous diverging scale should).
DIVERGING = ["#0d366b", "#1c5cab", "#3987e5", "#86b6ef", "#cde2fb", "#f0efec",
             "#f9d4d2", "#f0928d", "#e34948", "#a92b2b", "#6b1717"]
MID_GREY = "#f0efec"
C_RECEDE = "#e34948"
C_APPROACH = "#1c5cab"
# Mass: one neutral ramp from the page to ink, so it never competes with the
# flow colours drawn over it.
DENSITY = ["#fcfcfb", "#e8e7e1", "#c9c8c0", "#9d9c94", "#6b6a64", "#3a3936", "#0b0b0b"]
C_BASIN = ORANGE  # our own basin of attraction (our measurement)
C_STREAM = "#3a3936"  # streamlines
C_LCDM = "#9d9c94"  # Lambda-CDM expectation bands
C_PUBLISHED_TABLE = MUTED  # results from CF4's table as published (SDSS sign reversed)
C_CORRECTED = BLUE  # results from the corrected catalogue

# Image panels sit on black, as the pictures themselves do.
IMAGE_BG = "#000000"


def velocity_cmap():
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("flow", DIVERGING)


def density_cmap():
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("mass", DENSITY)


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
        fig.savefig(out, format=fmt, facecolor=fig.get_facecolor())
        written.append(out)
    plt.close(fig)
    return written
