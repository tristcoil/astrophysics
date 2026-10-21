"""Shared matplotlib styling for the project's figures.

The categorical colours are the atlas's colourblind-safe set (same palette and
chrome as 27-supernova-nickel and 05-supernova-cosmology). Scatter plots let
any two points sit side by side, so they may use only the first three slots:
blue, orange and aqua were run through the dataviz validator *all pairs* in
both themes (worst colour-blind dE 9.2, normal vision 24.0). Aqua is below 3:1
contrast on the light surface, so every category also gets its own marker
shape and a legend entry: colour is never the only cue.

Roles used across the figures:

* measurement method: stars = blue circles, gas = orange squares,
  masers = aqua diamonds; reverberation masses are calibrated *on* M-sigma, so
  they are context, drawn grey (triangles); upper limits are grey carets.
* bulge type (Kormendy & Ho): ellipticals = blue, classical bulges = orange,
  pseudobulges = aqua; masses Kormendy & Ho set aside = light grey.
* our two model-independent black holes (Milky Way, M87*): ink stars.
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

METHOD_STYLE = {
    "stars": {"color": BLUE, "marker": "o", "label": "stellar dynamics"},
    "gas": {"color": ORANGE, "marker": "s", "label": "gas dynamics (ionised gas, CO)"},
    "maser": {"color": AQUA, "marker": "D", "label": "water masers"},
    "reverb": {"color": CONTEXT, "marker": "^", "label": "reverberation (calibrated on M-σ)"},
}
BULGE_STYLE = {
    "elliptical": {"color": BLUE, "marker": "o", "label": "ellipticals"},
    "classical": {"color": ORANGE, "marker": "s", "label": "classical bulges"},
    "pseudo": {"color": AQUA, "marker": "D", "label": "pseudobulges"},
}
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
