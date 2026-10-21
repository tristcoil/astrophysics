"""Shared matplotlib styling: the atlas's validated palette and chrome.

Colour follows the *entity* in every figure, and only three hues ever overlap:

* the three data combinations the models are graded against --
  Planck 2018 + BK18 (our combination) -> blue, ACT DR6 P-ACT-LB + BK18 ->
  orange, SPT-3G CMB-SPA + DESI (x our BK18 r) -> aqua. These are the first
  three slots of the validated categorical order, which pass the colour-blind
  and normal-vision checks for *all* pairs (worst CVD dE 9.2, normal 24.0 on the
  light surface), so overlapping contours stay distinguishable. Aqua is below
  3:1 contrast, so every region is also labelled in words. Credible regions are
  one hue stepped light (95%) to darker (68%), never a rainbow.
* everything else -- Planck alone, the inflation models, BK18's own curves --
  is drawn in neutral ink and greys, told apart by line style, markers and
  direct labels.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
YELLOW = "#eda100"
MAGENTA = "#e87ba4"
GREEN = "#008300"
VIOLET = "#6250d6"
RED = "#e34948"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
NEUTRAL = "#8b8a84"
NEUTRAL_LIGHT = "#d8d7d1"

C_PLANCK = BLUE
C_ACT = ORANGE
C_SPT = AQUA

#: Sequential blue ramp (light = low).
SEQ_BLUE = ["#eef4fc", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95",
            "#0d366b"]


def lighten(color: str, amount: float) -> tuple:
    """Mix a colour toward the surface (amount 0 = unchanged, 1 = surface)."""
    c, s = to_rgb(color), to_rgb(SURFACE)
    return tuple((1 - amount) * a + amount * b for a, b in zip(c, s))


#: Light/dark steps of each data hue for 95% / 68% filled regions.
RAMPS = {
    BLUE: ("#b7d3f6", "#6da7ec"),
    ORANGE: ("#f8c9b4", "#f19a74"),
    AQUA: (lighten(AQUA, 0.70), lighten(AQUA, 0.35)),
    NEUTRAL: ("#e6e5e0", "#c9c8c2"),
}


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


def annotate(ax, text, xy, xytext, *, color=INK_SECONDARY, fontsize=7.5, **kw):
    return ax.annotate(text, xy=xy, xytext=xytext, color=color, fontsize=fontsize,
                       arrowprops={"arrowstyle": "-", "color": color, "linewidth": 0.7,
                                   "shrinkA": 2, "shrinkB": 3}, **kw)


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
