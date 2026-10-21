"""Shared matplotlib styling for the project's figures.

Identical in look to the other projects' (A1, D4), so the atlas reads as one
portfolio. The categorical colours are the atlas's validated colourblind-safe
set: adjacent slots keep a CVD colour difference of dE >= 8 (OKLab x100), and
scatter plots, where any two colours can meet, use at most the three hues that
pass the stricter all-pairs check (blue/orange/aqua). Aqua, yellow and magenta
sit below 3:1 contrast on the page, so anything drawn in them also carries a
direct text label.
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

# Semantic roles: one meaning per colour in every panel of every figure.
C_OBSERVED = INK  # what was measured: lensing masses, dispersions, our Einstein radii
C_OBSERVED_ERR = MUTED  # error bars, recessive
C_STARS = BLUE  # the stars (Chabrier IMF where an IMF is needed)
C_STARS_HEAVY = AQUA  # the stars with a Salpeter IMF (heavier, more low-mass stars)
C_DARK = ORANGE  # dark matter, and our own models and fits
C_PUBLISHED = VIOLET  # a published result we compare against
C_CONTEXT = "#d6d5ce"  # the whole population, drawn behind a highlighted subset

# Images: the lens galaxy in grey, the lensed light in a warm single-hue ramp.
CMAP_IMAGE = "gray"
CMAP_RING = "magma"
CMAP_RESIDUAL = "RdBu_r"  # diverging, neutral at zero


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
            "image.origin": "lower",
        }
    )


def annotate(ax, text: str, xy, xytext, *, color: str = INK_SECONDARY, **kwargs):
    """A callout with a thin leader line, in the project's annotation style."""
    return ax.annotate(text, xy=xy, xytext=xytext, color=color, fontsize=8,
                       arrowprops={"arrowstyle": "-", "color": color, "linewidth": 0.8,
                                   "shrinkA": 2, "shrinkB": 4}, **kwargs)


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
