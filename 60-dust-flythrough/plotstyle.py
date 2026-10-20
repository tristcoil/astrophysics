"""Shared matplotlib styling for the project's figures.

The categorical colours are the atlas's colourblind-safe set (same palette and
chrome as the other projects). Overlays on the dust maps use only the first
three slots (blue, orange, aqua), which pass the validator *all pairs* in both
themes; every overlay also has its own marker shape and a legend entry, so
colour is never the only cue.

Roles used across the figures:

* dust: one-hue grey ramp, light = little dust, dark = more (like the
  figures of Edenhofer et al. 2024, so ours can be compared side by side);
  the stretch (linear or log, and where it saturates) is stated on each panel;
* young star clusters = blue dots; the Radcliffe Wave model = orange line;
  molecular-cloud distances (Zucker et al. 2020) = aqua diamonds;
* the Local Bubble's shell = blue (its own map, or an outline);
* map differences = diverging blue <-> red around a grey zero;
* the Sun = ink star.
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
NEUTRAL = "#f0efec"

# Sequential blue ramp (ordered quantities other than dust).
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

#: Dust: a single grey hue from the surface to near-black.
DUST = LinearSegmentedColormap.from_list("dust", [SURFACE, "#c9c8c2", "#86857f", "#3d3c39", INK])
DUST.set_bad(SURFACE)
#: Sequential blue for magnitudes (shell distance, etc.).
SEQ_BLUE = LinearSegmentedColormap.from_list("seq_blue", BLUE_RAMP)
SEQ_BLUE.set_bad(NEUTRAL)
#: Diverging blue <-> grey <-> red for differences.
DIVERGING = LinearSegmentedColormap.from_list(
    "diverging", ["#184f95", "#3987e5", "#b7d3f6", NEUTRAL, "#f4b8b7", "#e34948", "#a32d2c"])
DIVERGING.set_bad(SURFACE)

SUN_STYLE = {"color": INK, "marker": "*", "ms": 11, "mec": SURFACE, "mew": 0.8, "ls": "none"}
CLUSTER_STYLE = {"color": BLUE, "marker": "o", "s": 9, "linewidths": 0.4, "edgecolors": SURFACE}
CLOUD_STYLE = {"color": AQUA, "marker": "D", "s": 12, "linewidths": 0.4, "edgecolors": SURFACE}
WAVE_STYLE = {"color": ORANGE, "lw": 1.8}


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
            "image.origin": "lower",
        }
    )


def map_axes(ax, *, grid: bool = False) -> None:
    """Axes that hold an image: no grid over the data, all four spines as a thin frame."""
    ax.grid(grid)
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color(AXIS)


def annotate(ax, text: str, xy, xytext, *, color: str = INK_SECONDARY, fontsize: float = 8, **kwargs):
    """A callout with a thin leader line, in the project's annotation style."""
    return ax.annotate(
        text, xy=xy, xytext=xytext, color=color, fontsize=fontsize,
        arrowprops={"arrowstyle": "-", "color": MUTED, "linewidth": 0.7, "shrinkA": 1, "shrinkB": 3},
        **kwargs,
    )


def caption(fig, text: str, y: float = -0.02) -> None:
    """Footnote under a figure: data provenance, caveats, the honest small print (wrapped to the figure width)."""
    import textwrap
    width = max(60, int(fig.get_figwidth() * 17))
    fig.text(0.5, y, "\n".join(textwrap.wrap(text, width)), ha="center", va="top", fontsize=7.5, color=MUTED)


def save(fig, path_stem, *, formats=("png",)) -> list:
    """Save a figure to every requested format, returning the paths written."""
    written = []
    for fmt in formats:
        out = f"{path_stem}.{fmt}"
        fig.savefig(out, format=fmt)
        written.append(out)
    plt.close(fig)
    return written
