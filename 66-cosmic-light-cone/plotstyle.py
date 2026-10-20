"""Shared matplotlib styling: the atlas's validated palette, applied by entity.

Colour follows the *thing* in every figure:

* our past light cone          -> ink (the one curve every figure is about)
* the Hubble sphere            -> blue   (and "recession speed when the light left")
* the event horizon            -> orange (and "recession speed today")
* the particle horizon         -> aqua
* worldlines, comparison models -> muted grey

Blue, orange and aqua were checked with the dataviz validator for all pairs, in
light and dark mode: worst CVD separation dE 9.2, normal-vision 24.0. Aqua sits
below 3:1 contrast on the light surface, so every curve carries a direct label.
DESI points share the blue of the quantity they measure (v at emission), the
supernovae the orange of theirs (v today).
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

C_CONE = INK
C_HUBBLE = BLUE
C_EVENT = ORANGE
C_PARTICLE = AQUA
C_WORLD = MUTED
C_THEN = BLUE
C_NOW = ORANGE

#: Light washes for regions (10%-ish of the hue).
WASH_BLUE = "#e6f0fb"
WASH_ORANGE = "#fdeee7"


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
        "lines.solid_capstyle": "round", "lines.solid_joinstyle": "round",
    })


def lighten(color: str, amount: float) -> tuple:
    """Mix a colour toward the surface (amount 0 = unchanged, 1 = surface)."""
    c, s = to_rgb(color), to_rgb(SURFACE)
    return tuple((1 - amount) * a + amount * b for a, b in zip(c, s))


def annotate(ax, text, xy, xytext, *, color=INK_SECONDARY, **kw):
    kw.setdefault("fontsize", 7.5)
    return ax.annotate(text, xy=xy, xytext=xytext, color=color,
                       arrowprops={"arrowstyle": "-", "color": color, "linewidth": 0.7,
                                   "shrinkA": 2, "shrinkB": 3}, **kw)


def caption(fig, text: str, y: float = -0.02) -> None:
    fig.text(0.5, y, text, ha="center", va="top", fontsize=7.2, color=MUTED, wrap=True)


def save(fig, path_stem, *, formats=("png",)) -> list:
    written = []
    for fmt in formats:
        out = f"{path_stem}.{fmt}"
        fig.savefig(out, format=fmt)
        written.append(out)
    plt.close(fig)
    return written
