"""Shared matplotlib styling for the project's figures.

The categorical colours are the atlas's colourblind-safe set (same palette and
chrome as 32-m-sigma-relation, 27-supernova-nickel and 05-supernova-cosmology).
Only the first three slots appear side by side in any one plot; they passed the
dataviz validator for all pairs in both themes (worst colour-blind dE 9.2).
Aqua is below 3:1 contrast on the light surface, so every category also gets
its own marker shape or a direct label: colour is never the only cue.

Roles used across the figures:

* **this project's measurements**: blue (the eRASS1 counts as we fit them)
* **the published eRASS1 result** (Ghirardini et al. 2024): orange
* **the cosmic microwave background** (Planck 2018, and project D3): aqua
* the data themselves (cluster counts): ink; context: light grey
* ordered quantities (sigma_8 steps, redshift): one blue ramp, light to dark
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

# --- Chrome and ink ---
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
CONTEXT = "#b9b8b1"
WASH = "#ecebe6"

# Sequential blue ramp (ordered quantities) and the diverging pair.
BLUE_100 = "#d6e6fb"
BLUE_150 = "#b7d3f6"
BLUE_300 = "#6da7ec"
BLUE_550 = "#1c5cab"
BLUE_750 = "#0e3566"
NEUTRAL = "#f0efec"
ORANGE_550 = "#c24a1a"

ROLE = {"ours": BLUE, "published": ORANGE, "cmb": AQUA, "data": INK, "context": CONTEXT}

#: Sequential map for redshift / exposure (one hue, light to dark).
SEQ_BLUE = LinearSegmentedColormap.from_list("seq_blue", [BLUE_100, BLUE_300, BLUE_550, BLUE_750])
#: Neutral sequential map for backgrounds (exposure map).
SEQ_GREY = LinearSegmentedColormap.from_list("seq_grey", ["#f4f3ef", "#d9d8d1", "#a9a8a1", "#6f6e69"])
#: Diverging residual map: blue (fewer than predicted) - neutral - orange (more).
DIVERGING = LinearSegmentedColormap.from_list("div", [BLUE_550, BLUE_300, NEUTRAL, "#f3a582", ORANGE_550])


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
            "axes.titlesize": 10.5,
            "axes.titleweight": "bold",
            "axes.titlelocation": "left",
            "axes.titlepad": 9,
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
    return ax.annotate(
        text, xy=xy, xytext=xytext, color=color, fontsize=fontsize,
        arrowprops={"arrowstyle": "-", "color": MUTED, "linewidth": 0.7, "shrinkA": 1, "shrinkB": 3},
        **kwargs,
    )


def caption(fig, text: str, y: float = -0.02) -> None:
    fig.text(0.5, y, text, ha="center", va="top", fontsize=7.5, color=MUTED, wrap=True)


def save(fig, path_stem, *, formats=("png", "svg")) -> list:
    written = []
    for fmt in formats:
        out = f"{path_stem}.{fmt}"
        fig.savefig(out, format=fmt)
        written.append(out)
    plt.close(fig)
    return written
