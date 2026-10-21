"""Shared matplotlib styling: the atlas's validated palette and chrome.

Colour follows the *entity* in every figure:

* this project's own measurements -> ink (stars or filled circles, so the
  colour of a marker never has to be decoded to find "ours");
* the model (hadron resonance gas) -> blue;
* Andronic et al. 2018's published fits -> orange; STAR's published fits ->
  aqua; other analyses -> neutral grey with their own marker shapes;
* lattice QCD -> violet (bands stepped light for 68%);
* the early universe -> blue (it never shares a figure with the model curve).

Blue, orange, aqua and violet pass the dataviz validator for all pairs on
the light surface (worst colour-blind dE 9.2, normal 16+); aqua is below 3:1
contrast, so everything drawn in it is also labelled in words.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
VIOLET = "#4a3aa7"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
NEUTRAL = "#8b8a84"

C_MODEL = BLUE
C_ANDRONIC = ORANGE
C_STAR = AQUA
C_LATTICE = VIOLET
C_OTHER = NEUTRAL
#: the early universe (never drawn together with the model, which is also blue)
C_UNIVERSE = BLUE

#: light steps for bands (68% washes)
WASH = {BLUE: "#cde2fb", ORANGE: "#fbe1d6", AQUA: "#d5f0e5", VIOLET: "#dcd8f4",
        NEUTRAL: "#e6e5e0"}
#: sequential blue ramp (light = low, dark = high)
SEQ_BLUE = ["#eef4fc", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95",
            "#0d366b"]


def use_style() -> None:
    mpl.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "savefig.bbox": "tight", "savefig.dpi": 200,
        "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 9,
        "axes.titlesize": 10.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
        "axes.titlepad": 8, "axes.labelsize": 9, "axes.labelcolor": INK_SECONDARY,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "axes.grid": True,
        "axes.axisbelow": True, "axes.spines.top": False, "axes.spines.right": False,
        "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-",
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "xtick.direction": "out", "ytick.direction": "out",
        "legend.frameon": False, "legend.fontsize": 8, "lines.linewidth": 1.6,
        "lines.markersize": 4, "text.color": INK, "figure.titlesize": 12,
        "figure.titleweight": "bold", "mathtext.default": "regular",
    })


def lighten(color: str, amount: float) -> tuple:
    c, s = to_rgb(color), to_rgb(SURFACE)
    return tuple((1 - amount) * a + amount * b for a, b in zip(c, s))


def caption(fig, text: str, y: float = -0.02) -> None:
    fig.text(0.5, y, text, ha="center", va="top", fontsize=7.2, color=MUTED, wrap=True)


def save(fig, path_stem, *, formats=("png", "svg")) -> list:
    written = []
    for fmt in formats:
        out = f"{path_stem}.{fmt}"
        fig.savefig(out, format=fmt)
        written.append(out)
    plt.close(fig)
    return written
