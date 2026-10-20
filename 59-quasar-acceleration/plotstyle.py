"""Shared matplotlib styling for the project's figures.

The categorical colours are the atlas's colourblind-safe set (the same palette
and chrome as 58-gravitational-redshift and 50-binary-orbits). Scatter-like
plots use only the first three slots, which the dataviz validator passes for
*all* pairs in both themes; aqua is below 3:1 contrast on the light surface,
so every category also gets its own marker shape and a legend entry: colour is
never the only cue.

Roles used across the figures:

* our measurement = blue; Klioner et al. (2021)'s published Gaia EDR3 result =
  ink (open markers); predictions from Galactic constants = orange; VLBI
  measurements = aqua;
* on sky maps: averaged quasar proper motions = ink arrows, the fitted glide
  (the acceleration's pattern) = orange arrows; Galactic Centre = a star.
* diverging maps: blue <-> red about a neutral grey (dataviz reference pair).
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# --- Categorical slots, assigned in fixed order ---
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
YELLOW = "#eda100"
MAGENTA = "#e87ba4"
GREEN = "#008300"
VIOLET = "#6250d6"
RED = "#e34948"
SERIES = [BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED]

# --- Chrome and ink ---
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
CONTEXT = "#b9b8b1"
WASH = "#ecebe6"
NEUTRAL = "#f0efec"

# Sequential blue ramp and the diverging blue <-> red pair (dataviz reference).
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
RED_RAMP = ["#fbd5d4", "#f5a9a8", "#ee7d7c", "#e34948", "#c13332", "#9a2525", "#701a1a"]
SEQ_BLUE = LinearSegmentedColormap.from_list("seq_blue", [SURFACE] + BLUE_RAMP)
DIVERGING = LinearSegmentedColormap.from_list(
    "blue_red", list(reversed(BLUE_RAMP))[:-1] + [NEUTRAL] + RED_RAMP[1:])

OURS = BLUE
PUBLISHED = INK
PREDICTED = ORANGE
VLBI = AQUA
DATA_ARROW = INK
MODEL_ARROW = ORANGE


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
            "grid.linestyle": "-",
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
