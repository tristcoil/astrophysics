"""Shared matplotlib styling for the project's figures.

Identical in look to the other atlas projects, so they read as one portfolio. The
categorical colours are a colourblind-safe set: adjacent slots keep a CVD colour
difference of dE >= 8 (OKLab x100), so line charts stay readable for
deuteranopia/protanopia and in greyscale print. Scatter plots, where any two
colours can meet, use at most three hues (blue/orange/aqua or blue/orange/violet,
the triples that pass the stricter all-pairs check). Aqua, yellow and magenta sit
below 3:1 contrast on the page, so anything drawn in them also carries a text label.

Magnitudes -- the sunspot area in the butterfly diagram -- get a single-hue
sequential ramp, orange stepped from light to dark at constant OKLCH hue 40.7 deg
(lightness 0.96 -> 0.36, strictly monotone), never a rainbow.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

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
WASH = "#f0efec"  # the neutral wash for shaded episodes (the diverging midpoint grey)

# --- Sequential: sunspot area (OKLCH hue 40.7, L 0.96 -> 0.36) ---
SUN_RAMP = ["#ffeee8", "#ffccba", "#fea98a", "#fe8050", "#e76530",
            "#ca4f19", "#aa3e0a", "#88320b", "#642912"]
#: The butterfly colormap starts one step in, so a single small spot group still
#: shows against the page; empty bins stay the page colour.
SUN_CMAP = LinearSegmentedColormap.from_list("sun", SUN_RAMP[1:])

# Semantic roles used across the solar-cycle figures. Naming them once keeps the
# same quantity the same colour in every panel.
C_MONTHLY = MUTED      # monthly means: recessive, the noisy raw record
C_SMOOTHED = BLUE      # the 13-month smoothed number: the signal we read
C_FOCUS = ORANGE       # the one thing a panel is about (Cycle 25, a fitted law)
C_LITERATURE = VIOLET  # a published value or law we compare against
C_NORTH = BLUE
C_SOUTH = ORANGE


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


def caption(fig, text: str, y: float = -0.02) -> None:
    """Footnote under a figure: data provenance, caveats, the honest small print."""
    fig.text(0.5, y, text, ha="center", va="top", fontsize=7.5, color=MUTED, wrap=True)


def panel_label(ax, text: str) -> None:
    """Bold panel title, left-aligned, in the axes title slot."""
    ax.set_title(text, loc="left")


def save(fig, path_stem, *, formats=("png", "svg")) -> list:
    """Save a figure to every requested format, returning the paths written."""
    written = []
    for fmt in formats:
        out = f"{path_stem}.{fmt}"
        fig.savefig(out, format=fmt)
        written.append(out)
    plt.close(fig)
    return written
