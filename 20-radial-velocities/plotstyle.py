"""Shared matplotlib styling, identical in look to projects A1 and A2.

Colour roles in this project:

* **Spectrographs** are the categorical series. Velocity plots are scatter
  plots, where every pair of colours can meet, so the four families use a
  subset of the palette that passes the all-pairs checks (validated with the
  dataviz validator, ``--pairs all``, light surface: worst CVD dE 9.2, worst
  normal-vision dE 16.3). Aqua sits below 3:1 contrast on the surface, so each
  family also has its own marker shape and a legend entry: identity never rests
  on colour alone.
* **Model curves** are ink, so no data colour does double duty.
* In comparison figures without instrument colours, orange is our own
  measurement and violet a published one, as in A1 and A2.
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
CONTEXT = "#d6d5ce"  # background population / de-emphasised data

# --- Roles --------------------------------------------------------------------------
C_MODEL = INK
C_OURS = ORANGE
C_PUBLISHED = VIOLET
C_SIGNAL = BLUE  # a single-series periodogram or light curve

#: Spectrograph family -> (colour, marker). Fixed in every figure.
FAMILY_STYLE = {
    "ELODIE": (VIOLET, "s"),
    "Lick": (ORANGE, "D"),
    "HIRES": (AQUA, "^"),
    "HARPS": (BLUE, "o"),
}
FAMILY_LABEL = {
    "ELODIE": "ELODIE (Haute-Provence)",
    "Lick": "Lick / Hamilton",
    "HIRES": "Keck / HIRES",
    "HARPS": "HARPS (La Silla)",
}


def family(instrument: str) -> str:
    """'Lick dewar 13' -> 'Lick', 'HARPS pre-2015' -> 'HARPS', ..."""
    return instrument.split()[0]


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
