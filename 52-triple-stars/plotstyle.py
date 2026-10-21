"""Shared matplotlib styling (the same look as the other atlas projects).

Colour roles in this project:

* **The three eclipse-switching systems** are categorical slots 1-3: SS Lac
  blue, HS Hya orange, V907 Sco aqua. The first three slots pass the
  all-pairs colour-vision checks in both modes, so they can share a plot.
* **Our measurements** (plate and TESS points, our fits) are ink; **published
  values and published models** are violet, as in the other projects.
* **Magnitudes** (MEGNO chaos indicator, survival time) use the one-hue blue
  ramp, light to dark.
* The Kozai experiments use three fixed slots (red: point masses; yellow:
  + relativity; blue: + tides), always direct-labelled.
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

#: the three eclipse-switching systems (fixed slots 1-3: these three validate
#: all-pairs in both modes, so they may share a scatter plot)
SYSTEM_COLOR = {"SS Lac": BLUE, "HS Hya": ORANGE, "V907 Sco": AQUA}
SYSTEM_COLOR_DARK = {"SS Lac": "#3987e5", "HS Hya": "#d95926", "V907 Sco": "#199e70"}
#: the three Kozai experiments (no new physics -> more physics)
RUN_COLOR = {"newton": RED, "gr": YELLOW, "tides": BLUE}
#: stars of TIC 290061484
STAR_COLOR = {"Aa": BLUE, "Ab": ORANGE, "B": AQUA, "C": YELLOW}
C_OURS = INK
C_PUBLISHED = VIOLET
C_DASCH = INK_SECONDARY
#: sequential ramp (blue) for magnitudes such as MEGNO
BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
             "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]


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
