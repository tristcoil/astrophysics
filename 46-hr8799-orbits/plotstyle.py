"""Shared matplotlib styling (the same look as the other atlas projects).

Colour roles in this project:

* **Planets** are the categorical series, in fixed order from the inside out:
  e = slot 1 (blue), d = slot 2 (orange), c = slot 3 (aqua), b = slot 4
  (yellow). On the sky the four orbits are nested rings, so only neighbours
  ever meet (e-d, d-c, c-b): the adjacent-pair checks apply, and they pass in
  both modes (dataviz validator: worst CVD dE 9.1 light / 8.4 dark, normal
  vision 22.9 / 19.8). Aqua and yellow are below 3:1 contrast on the light
  surface, so every planet is also direct-labelled with its letter. Panels in
  which all four would overlap are split into small multiples instead (the
  all-pairs check fails for yellow vs orange).
* **Measured positions** are drawn in ink with a surface ring; the planet
  colour belongs to orbits and posteriors.
* **Published values** (Wang et al. 2018 and others) are violet, ours are ink,
  where a figure compares the two without planet colours.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

# --- Categorical slots (validated reference palette) --------------------------
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
YELLOW = "#eda100"
MAGENTA = "#e87ba4"
GREEN = "#008300"
VIOLET = "#4a3aa7"
RED = "#e34948"

# --- Chrome and ink ---------------------------------------------------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
CONTEXT = "#d6d5ce"

#: planet -> colour, fixed in every figure (inner to outer: slots 1-4)
PLANET_COLOR = {"e": BLUE, "d": ORANGE, "c": AQUA, "b": YELLOW}
#: the same slots stepped for the dark surface (interactive page)
PLANET_COLOR_DARK = {"e": "#3987e5", "d": "#d95926", "c": "#199e70", "b": "#c98500"}
C_OURS = INK
C_PUBLISHED = VIOLET


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
            "svg.fonttype": "none",
        }
    )


def annotate(ax, text: str, xy, xytext, *, color: str = INK_SECONDARY, relpos=(0.5, 0.5),
             fontsize: float = 8, **kwargs):
    """A callout with a thin leader line, in the project's annotation style."""
    return ax.annotate(
        text, xy=xy, xytext=xytext, color=color, fontsize=fontsize,
        arrowprops={"arrowstyle": "-", "color": color, "linewidth": 0.8,
                    "shrinkA": 2, "shrinkB": 4, "relpos": relpos},
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
