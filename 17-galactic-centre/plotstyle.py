"""Shared matplotlib styling for the project's figures (the atlas house style, as in A1-E2).

The categorical colours are a colourblind-safe set: adjacent slots keep a CVD
colour difference of dE >= 8 (OKLab x100), so the figures stay readable for
deuteranopia/protanopia and in greyscale print.
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

# Semantic roles. Each team keeps one colour in every figure; our own joint
# result is ink, published values are violet (as in the other projects), and
# general relativity's prediction is aqua. Time is shown with a sequential
# colormap, never with the categorical slots.
C_KECK = BLUE  # UCLA / Keck data and fits
C_VLT = ORANGE  # MPE / VLT data and fits
C_JOINT = INK  # both teams together (our headline)
C_PUBLISHED = VIOLET
C_GR = AQUA  # Einstein's prediction
C_NEWTON = MUTED  # Newton's prediction, drawn dashed
C_CONTEXT = MUTED
TIME_CMAP = "viridis"


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
