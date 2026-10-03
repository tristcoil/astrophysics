"""Shared matplotlib styling for the project's figures.

The categorical colours are a colourblind-safe set: adjacent slots keep a CVD
colour difference of dE >= 8 (OKLab x100), so the figures stay readable for
deuteranopia/protanopia and in greyscale print.
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

# Semantic roles used across the cluster figures. Naming them once keeps the
# same thing the same colour in every panel. Most panels are an "emphasis"
# form: the cluster in the accent hue, the field as recessive grey context.
C_MEMBER = BLUE      # cluster members
C_FIELD = MUTED      # field stars drawn as points
C_REJECTED = ORANGE  # picked by motion alone, rejected by parallax
C_MODEL = ORANGE     # a prediction drawn against measured data
C_PUBLISHED = VIOLET

#: Clusters on a shared plot, in the fixed categorical order. Scatter plots
#: put every pair side by side, and only the first three slots pass the
#: all-pairs colour-vision check -- which is exactly the number needed.
CLUSTER_COLOURS = [BLUE, ORANGE, AQUA]

# Sequential ramps (one hue, light -> dark) for density maps.
BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
             "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
CMAP_BLUE = LinearSegmentedColormap.from_list("ramp_blue", BLUE_RAMP)
#: Neutral ramp for context densities -- the field, drawn so it recedes.
CMAP_CONTEXT = LinearSegmentedColormap.from_list(
    "ramp_context", ["#eceae4", "#c3c2b7", "#a3a198", "#77756f", "#52514e", "#2b2a27"])


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
    """Footnote under a figure: data provenance, caveats, the honest small print.

    Placed just below the lowest thing already drawn (axis labels included),
    so the gap between plot and caption is the same in every figure.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    lowest = min(ax.get_tightbbox(renderer).y0 for ax in fig.axes)
    fig.text(
        0.5,
        lowest / fig.bbox.height - 0.015,
        text,
        ha="center",
        va="top",
        fontsize=7.5,
        color=MUTED,
        wrap=True,
    )


#: Resolution of the rasterised layers (star clouds, density maps) inside SVGs,
#: which exist for their vector text and lines. Thousands of individually drawn
#: vector markers made the star-field figures ~3 MB each; as 150 dpi rasters
#: they cost a fraction of that and look the same at any normal size.
SVG_RASTER_DPI = 150


def save(fig, path_stem, *, formats=("png", "svg")) -> list:
    """Save a figure to every requested format, returning the paths written."""
    written = []
    for fmt in formats:
        out = f"{path_stem}.{fmt}"
        fig.savefig(out, format=fmt, **({"dpi": SVG_RASTER_DPI} if fmt == "svg" else {}))
        written.append(out)
    plt.close(fig)
    return written
