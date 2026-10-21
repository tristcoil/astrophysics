"""Shared matplotlib styling: the atlas's validated palette and chrome.

Colour follows the *entity* in every figure and never changes between them:

* the two clocks of the early universe: Big Bang nucleosynthesis (deuterium,
  3 minutes) -> orange, the CMB (380,000 years) -> violet;
* predictions from our own calculations -> blue (bands one hue, light to
  dark), the annihilation of matter with antimatter in particular;
* a fourth family where needed (other experiments, alternative inputs) ->
  aqua, always direct-labelled (it is under 3:1 contrast);
* measurements we compare against -> ink markers, other experiments -> grey.

Blue, orange, aqua and violet (#4a3aa7) pass the dataviz validator for all
pairs on the light surface (worst colour-blind dE 9.2, normal vision 16.3).
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
EXCLUDED = "#e9e8e3"

C_BBN = ORANGE
C_CMB = VIOLET
C_MODEL = BLUE
C_ALT = AQUA
C_DATA = INK
C_OTHER = "#8b8a84"

#: Light / dark steps of each hue for 95% / 68% filled regions.
RAMPS = {
    BLUE: ("#b7d3f6", "#6da7ec"),
    ORANGE: ("#f8c9b4", "#f19a74"),
    VIOLET: ("#cfc9f0", "#8f84e0"),
    AQUA: ("#bfe9d8", "#6fd0ab"),
    INK: ("#c9c8c2", "#6f6e69"),
}

#: Sequential blue ramp (light = small, dark = large).
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
    """Mix a colour toward the surface (amount 0 = unchanged, 1 = surface)."""
    c, s = to_rgb(color), to_rgb(SURFACE)
    return tuple((1 - amount) * a + amount * b for a, b in zip(c, s))


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
