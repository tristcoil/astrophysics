"""Shared matplotlib styling for the project's figures.

The palette is the atlas-wide one (D1, B4, B6): a colourblind-safe categorical set,
checked with the dataviz validator rather than by eye. Slots 1-3 (blue, orange,
aqua) pass the colour-blindness checks for all pairs at once, so they carry every
scatter plot; each category also has its own marker shape, so colour is never the
only cue. Roles:

* photometric bands: ZTF g = blue circles, ZTF r = orange squares, ZTF i = aqua
  triangles; PTF R (2009-2012, a different telescope) = ink diamonds;
* LISA binaries by type: AM CVn = blue circles, detached double white dwarfs =
  orange squares, hot subdwarf + white dwarf = aqua triangles, the ultracompact
  X-ray binary = grey pentagon, labelled; filled = Kupfer et al.'s verification binaries
  (detectable within 3 months), open = detectable within 4 years;
* ZTF J1539+5027 itself: an ink star wherever it appears;
* published ephemerides and curves we check against: muted grey lines;
* LISA's instrument noise: ink; the Galactic foreground: violet (our model) and
  muted (Robson et al.'s population-synthesis fit).
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

# --- Categorical slots, assigned in fixed order (never cycled past slot 8) ---
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
YELLOW = "#eda100"
MAGENTA = "#e87ba4"
GREEN = "#008300"
VIOLET = "#4a3aa7"
RED = "#e34948"

# --- Chrome and ink ---
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
CONTEXT = "#b9b8b1"  # de-emphasised points (context, not the story)
WASH = "#ecebe6"  # uncertainty bands

BAND_STYLE = {
    "zg": {"color": BLUE, "marker": "o", "label": "ZTF g"},
    "zr": {"color": ORANGE, "marker": "s", "label": "ZTF r"},
    "zi": {"color": AQUA, "marker": "^", "label": "ZTF i"},
    "ptf": {"color": INK, "marker": "D", "label": "PTF R"},
}
TYPE_STYLE = {
    "AMCVn": {"color": BLUE, "marker": "o", "label": "AM CVn (mass-transferring)"},
    "DWD": {"color": ORANGE, "marker": "s", "label": "detached white-dwarf pair"},
    "sdB": {"color": AQUA, "marker": "^", "label": "hot subdwarf + white dwarf"},
    # one object only: secondary ink + its own shape + a direct label (a fourth hue,
    # yellow, would sit next to orange, a pair the palette validator fails for scatters)
    "UCXB": {"color": INK_SECONDARY, "marker": "p", "label": "ultracompact X-ray binary"},
}
C_J1539 = INK
C_INSTRUMENT = INK
C_FOREGROUND_OURS = VIOLET
C_FOREGROUND_ROBSON = MUTED
C_PUBLISHED = MUTED

# Sequential blue ramp (ordered quantities: time to contact, density).
BLUE_150 = "#b7d3f6"
BLUE_300 = "#6da7ec"
BLUE_450 = "#2a78d6"
BLUE_550 = "#1c5cab"
BLUE_650 = "#104281"
SEQ_BLUES = [BLUE_150, BLUE_300, BLUE_450, BLUE_550, BLUE_650]
NEUTRAL = "#f0efec"


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
        }
    )


def annotate(ax, text: str, xy, xytext, *, color: str = INK_SECONDARY, **kwargs):
    """A callout with a thin leader line, in the project's annotation style."""
    return ax.annotate(text, xy=xy, xytext=xytext, color=color, fontsize=8,
                       arrowprops={"arrowstyle": "-", "color": color, "linewidth": 0.8, "shrinkA": 2,
                                   "shrinkB": 4}, **kwargs)


def caption(fig, text: str, y: float = -0.02) -> None:
    """Footnote under a figure: data provenance, caveats, the honest small print."""
    fig.text(0.5, y, text, ha="center", va="top", fontsize=7.5, color=MUTED, wrap=True)


def save(fig, path_stem, *, formats=("png", "svg")) -> list:
    """Save a figure to every requested format, returning the paths written."""
    written = []
    for fmt in formats:
        out = f"{path_stem}.{fmt}"
        fig.savefig(out, format=fmt)
        written.append(out)
    plt.close(fig)
    return written
