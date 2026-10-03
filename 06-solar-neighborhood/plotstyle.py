"""Shared matplotlib styling for this project's figures: a dark, night-sky theme.

Why dark, when project A1 is light: most of these figures show stars in their
true colours, and a hot star's true colour is a pale blue-white that vanishes
on a white page. Every value here is the dark-mode step of the same documented
palette A1 uses; the categorical slots were validated for the dark surface
(adjacent CVD dE >= 8.4, all-pairs for the first three slots, >= 3:1 contrast).
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# --- Categorical slots (dark-surface steps), assigned in fixed order ---
BLUE = "#3987e5"
ORANGE = "#d95926"
AQUA = "#199e70"
YELLOW = "#c98500"
MAGENTA = "#d55181"
GREEN = "#008300"
VIOLET = "#9085e9"
RED = "#e66767"
SERIES = [BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED]

# --- Chrome and ink ---
SURFACE = "#1a1a19"
PAGE = "#0d0d0d"
INK = "#ffffff"
INK_SECONDARY = "#c3c2b7"
MUTED = "#898781"
GRID = "#2c2c2a"
AXIS = "#383835"

# --- Sequential ramp (blue), dark-mode orientation: near-zero recedes into the
# surface, high density is light. Steps 700 -> 100 of the documented ramp.
BLUE_RAMP = ["#0d366b", "#104281", "#184f95", "#1c5cab", "#256abf", "#2a78d6", "#3987e5",
             "#5598e7", "#6da7ec", "#86b6ef", "#9ec5f4", "#b7d3f6", "#cde2fb"]
DENSITY_CMAP = LinearSegmentedColormap.from_list("density_dark", BLUE_RAMP)

# Semantic roles, fixed across every figure so a quantity keeps its colour.
C_GAIA = BLUE
C_HIPPARCOS = ORANGE
C_EMPHASIS = BLUE
C_CONTEXT = MUTED       # de-emphasised marks: the documented muted grey
C_PREDICTION = MUTED


def use_style() -> None:
    """Apply the project-wide rcParams. Call once at import time in figure code."""
    mpl.rcParams.update({
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
        "axes.titlecolor": INK,
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
        "xtick.labelcolor": MUTED,
        "ytick.labelcolor": MUTED,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.frameon": False,
        "legend.fontsize": 8,
        "legend.labelcolor": INK_SECONDARY,
        "lines.linewidth": 2.0,
        "lines.solid_capstyle": "round",
        "lines.markersize": 4,
        "text.color": INK,
        "figure.titlesize": 12,
        "figure.titleweight": "bold",
    })


def annotate(ax, text: str, xy, xytext, *, color: str = INK_SECONDARY, **kwargs):
    """A callout with a thin leader line."""
    return ax.annotate(text, xy=xy, xytext=xytext, color=color, fontsize=8,
                       arrowprops={"arrowstyle": "-", "color": color, "linewidth": 0.8,
                                   "shrinkA": 2, "shrinkB": 4}, **kwargs)


def caption(fig, text: str, y: float = -0.02) -> None:
    """Footnote under a figure: provenance, caveats, the honest small print."""
    fig.text(0.5, y, text, ha="center", va="top", fontsize=7.5, color=MUTED, wrap=True)


def subtitle(ax, text: str) -> None:
    """A one-line takeaway under the title, in secondary ink."""
    ax.text(0.0, 1.0, text, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=8.5, color=INK_SECONDARY)


def temperature_key(fig, rect, colours, temps=(2500, 3000, 4000, 5000, 5800, 7500, 10000)):
    """A horizontal strip showing which colour means which surface temperature.

    Star colour is the one place a multi-hue scale is allowed here: it is not a
    palette choice but the physical colour of a blackbody at each temperature.
    That makes a key obligatory, so nobody has to guess.
    """
    ax = fig.add_axes(rect)
    t = np.geomspace(2300, 12000, 400)
    ax.imshow(colours(t)[None, :, :], aspect="auto", extent=(np.log10(t[0]), np.log10(t[-1]), 0, 1))
    ax.set_yticks([])
    ax.set_xticks(np.log10(temps))
    ax.set_xticklabels([("Sun\n" if tt == 5800 else "") + f"{tt:,} K" for tt in temps], fontsize=7)
    ax.tick_params(length=2, colors=MUTED)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("star colour = surface temperature", fontsize=7.5, color=INK_SECONDARY,
                 fontweight="normal", loc="left", pad=3)
    return ax


def save(fig, path_stem, *, formats=("png", "svg")) -> list:
    """Save a figure to every requested format, returning the paths written."""
    written = []
    for fmt in formats:
        out = f"{path_stem}.{fmt}"
        fig.savefig(out, format=fmt)
        written.append(out)
    plt.close(fig)
    return written
