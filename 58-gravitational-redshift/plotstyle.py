"""Shared matplotlib styling for the project's figures.

The categorical colours are the atlas's colourblind-safe set (the same palette
and chrome as 50-binary-orbits and 27-supernova-nickel). Scatter plots let any
two points sit side by side, so they use only the first three slots: blue,
orange and aqua were run through the dataviz validator *all pairs* in both
themes (worst colour-blind dE 9.2, normal vision 24.0). Aqua is below 3:1
contrast on the light surface, so every category also gets its own marker
shape and a legend entry: colour is never the only cue.

Roles used across the figures:

* satellites: E14 (GSAT-0202) = blue circles, E18 (GSAT-0201) = orange squares,
  nominal Galileo controls = aqua diamonds (or grey context);
* clock budget: gravity (height) = blue, speed (time dilation) = orange,
  their sum = ink;
* general relativity's prediction = ink line; published results = grey;
  our measurements = blue/orange.
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
CONTEXT = "#b9b8b1"  # de-emphasised points (context, not the story)
WASH = "#ecebe6"  # scatter band about a fit

# Sequential blue ramp (ordered quantities) and the diverging neutral.
BLUE_150 = "#b7d3f6"
BLUE_300 = "#6da7ec"
BLUE_550 = "#1c5cab"
NEUTRAL = "#f0efec"

SAT_STYLE = {
    "E14": {"color": BLUE, "marker": "o", "label": "E14 (GSAT-0202)"},
    "E18": {"color": ORANGE, "marker": "s", "label": "E18 (GSAT-0201)"},
    "control": {"color": AQUA, "marker": "D", "label": "nominal Galileo (control)"},
}
GRAVITY = BLUE
SPEED = ORANGE
ANCHOR_STYLE = {"color": INK, "marker": "*", "ms": 15, "mec": SURFACE, "mew": 1.0}


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


def annotate(ax, text: str, xy, xytext, *, color: str = INK_SECONDARY, fontsize: float = 8, **kwargs):
    """A callout with a thin leader line, in the project's annotation style."""
    return ax.annotate(
        text, xy=xy, xytext=xytext, color=color, fontsize=fontsize,
        arrowprops={"arrowstyle": "-", "color": MUTED, "linewidth": 0.7, "shrinkA": 1, "shrinkB": 3},
        **kwargs,
    )


def caption(fig, text: str, y: float = -0.02) -> None:
    """Footnote under a figure: data provenance, caveats, the honest small print."""
    fig.text(0.5, y, text, ha="center", va="top", fontsize=7.5, color=MUTED, wrap=True)


def save(fig, path_stem, *, formats=("png", "svg")) -> list:
    """Save a figure to every requested format, returning the paths written."""
    written = []
    for fmt in formats:
        out = f"{path_stem}.{fmt}"
        fig.savefig(out, format=fmt)
        written.append(out)
    plt.close(fig)
    return written
