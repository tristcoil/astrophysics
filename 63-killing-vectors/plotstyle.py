"""Shared matplotlib styling for N6's figures.

The categorical colours are the atlas's colour-blind-safe set (the same palette
and chrome as 58-gravitational-redshift and 50-binary-orbits). Scatter plots,
where any two series may sit side by side, use only the first three slots
(blue, orange, aqua), which pass the dataviz validator on all pairs in both
themes. Aqua is below 3:1 contrast on the light surface, so every series also
gets its own marker or line style and a legend entry: colour is never the only cue.

Roles used across the figures:

* general relativity / Killing quantities = blue; Newtonian expressions = orange;
* prograde orbits = blue, retrograde = orange; the horizon = ink;
* Kerr = blue, Zipoy-Voorhees (no hidden symmetry) = orange;
* measurements = ink points with error bars; model or posterior bands = wash;
* continuum-fitting spins = blue circles, reflection spins = orange squares,
  gravitational-wave remnants = aqua diamonds.
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

# Sequential blue ramp and the diverging pair (blue <-> red, grey midpoint).
BLUE_100 = "#cde2fb"
BLUE_150 = "#b7d3f6"
BLUE_300 = "#6da7ec"
BLUE_550 = "#1c5cab"
BLUE_700 = "#0d366b"
NEUTRAL = "#f0efec"
SEQ_BLUE = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
            "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]

GR = BLUE
NEWTON = ORANGE
PRO = BLUE
RETRO = ORANGE


def seq_cmap():
    return mpl.colors.LinearSegmentedColormap.from_list("seq_blue", SEQ_BLUE)


def diverging_cmap():
    return mpl.colors.LinearSegmentedColormap.from_list("div", ["#104281", "#3987e5", NEUTRAL, "#e66767", "#a1292a"])


def use_style() -> None:
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
            "mathtext.default": "regular",
        }
    )


def annotate(ax, text: str, xy, xytext, *, color: str = INK_SECONDARY, fontsize: float = 8, **kwargs):
    return ax.annotate(
        text, xy=xy, xytext=xytext, color=color, fontsize=fontsize,
        arrowprops={"arrowstyle": "-", "color": MUTED, "linewidth": 0.7, "shrinkA": 1, "shrinkB": 3},
        **kwargs,
    )


def caption(fig, text: str, y: float = -0.02) -> None:
    """Footnote under a figure: data provenance, caveats, the honest small print."""
    fig.text(0.5, y, text, ha="center", va="top", fontsize=7.5, color=MUTED, wrap=True)


def save(fig, path_stem, *, formats=("png", "svg")) -> list:
    written = []
    for fmt in formats:
        out = f"{path_stem}.{fmt}"
        fig.savefig(out, format=fmt)
        written.append(out)
    plt.close(fig)
    return written
