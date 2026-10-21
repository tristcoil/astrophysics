"""Shared matplotlib styling, identical in look to projects A1-A3.

Colour roles in this project (fixed in every figure):

* **Our measurement vs a published one**: orange vs violet, as in A1-A3.
* **Classes of compact object** (graveyard, mass distributions): black holes
  weighed by light (X-ray binaries) aqua, black holes heard in gravitational
  waves blue, Gaia's dormant black holes green, neutron stars weighed by light
  yellow, neutron stars in gravitational waves magenta. Blue, orange and
  aqua, the first three slots of the validated reference palette, pass the
  all-pairs colour-vision checks; where more classes share a panel they sit
  in separate columns or rows and carry their own marker shape and label,
  so identity never rests on colour alone. Aqua, yellow and magenta fall
  below 3:1 contrast on the surface, so their marks are always labelled.
* **X-ray bands** (light curves): a sequential blue ramp from soft to hard.
* **Model curves and annotations** are ink.
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

# --- Sequential blue ramp (reference palette steps) ----------------------------
BLUE_RAMP = ["#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

# --- Chrome and ink ---------------------------------------------------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
CONTEXT = "#d6d5ce"  # background population / de-emphasised data
GAP_FILL = "#efeee9"  # the shaded mass gap

# --- Roles --------------------------------------------------------------------------
C_OURS = ORANGE
C_PUBLISHED = VIOLET
C_MODEL = INK
C_BH_LIGHT = AQUA
C_BH_GW = BLUE
C_BH_GAIA = GREEN
C_NS_LIGHT = YELLOW
C_NS_GW = MAGENTA

#: X-ray band -> colour (soft to hard on one blue ramp)
BAND_COLOURS = {"2-4": BLUE_RAMP[0], "4-10": BLUE_RAMP[2], "10-20": BLUE_RAMP[3], "bat": BLUE_RAMP[5]}
BAND_LABELS = {"2-4": "MAXI 2–4 keV", "4-10": "MAXI 4–10 keV", "10-20": "MAXI 10–20 keV",
               "bat": "Swift/BAT 15–50 keV"}


def use_style() -> None:
    """Apply the project-wide rcParams."""
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
