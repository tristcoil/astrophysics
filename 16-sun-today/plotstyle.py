"""Shared matplotlib styling for the project's figures.

The charts share the atlas-wide look of the earlier projects (A1, C2, D4): the same
colourblind-safe categorical order, surfaces and ink. Solar images get a dark
counterpart, because the Sun sits on black sky and a pale frame around it would
glare; its ink steps are the same palette's dark-mode values.
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

# --- Chrome and ink (light charts) ---
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

# --- Dark surface for solar images ---
DARK_PAGE = "#0d0d0d"
DARK_INK = "#ffffff"
DARK_INK_SECONDARY = "#c3c2b7"
DARK_MUTED = "#898781"
DARK_BLUE = "#3987e5"
DARK_ORANGE = "#d95926"
DARK_AQUA = "#199e70"

# Semantic roles, kept the same across figures: our measurement is ink/blue, the
# model we fit is orange, a published comparison is violet.
C_DATA = INK
C_OURS = BLUE
C_FIT = ORANGE
C_PUBLISHED = VIOLET
C_CONTEXT = "#d6d5ce"


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


def dark_figure(*args, **kwargs):
    """A figure on the dark page, for solar images."""
    fig = plt.figure(*args, facecolor=DARK_PAGE, **kwargs)
    return fig


def annotate(ax, text: str, xy, xytext, *, color: str = INK_SECONDARY, fontsize: float = 8,
             **kwargs):
    """A callout with a thin leader line, in the project's annotation style."""
    return ax.annotate(
        text, xy=xy, xytext=xytext, color=color, fontsize=fontsize,
        arrowprops={"arrowstyle": "-", "color": color, "linewidth": 0.8,
                    "shrinkA": 2, "shrinkB": 4},
        **kwargs,
    )


def caption(fig, text: str, *, color: str = MUTED, y: float = -0.02) -> None:
    """Footnote under a figure: data provenance, caveats, the honest small print."""
    fig.text(0.5, y, text, ha="center", va="top", fontsize=7.5, color=color, wrap=True)


def save(fig, path_stem, *, formats=("png", "svg"), dpi=None) -> list:
    """Save a figure to every requested format, returning the paths written."""
    written = []
    for fmt in formats:
        out = f"{path_stem}.{fmt}"
        kw = {"format": fmt, "facecolor": fig.get_facecolor()}
        if dpi:
            kw["dpi"] = dpi
        fig.savefig(out, **kw)
        written.append(out)
    plt.close(fig)
    return written
