"""Shared matplotlib styling for the project's figures.

The categorical colours are the atlas's colourblind-safe set (same palette and
chrome as 50-binary-orbits and 27-supernova-nickel). Scatter plots let any two
points sit side by side, so they may use only the first three slots: blue,
orange and aqua were run through the dataviz validator *all pairs* in both
themes (worst colour-blind dE 9.2, normal vision 24.0). Aqua is below 3:1
contrast on the light surface, so every category also gets its own marker
shape and a legend entry: colour is never the only cue.

Roles used across the figures:

* 1919 instruments: Sobral 4-inch = blue circles, Principe = orange squares,
  Sobral astrographic = aqua diamonds. The 1922 Wallal stars get their own
  panel (blue circles) rather than a fourth colour.
* theory: general relativity = ink solid line; the Newtonian half deflection
  = ink dashed line; no deflection = the zero grid line.
* the century of gamma: starlight (optical) = blue, radio waves = orange,
  time delay (radar and spacecraft) = aqua; marker shape separates eclipses
  from space astrometry, and radio interferometry from VLBI.
"""

from __future__ import annotations

import textwrap

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
SERIES = [BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED]

# --- Chrome and ink ---
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
CONTEXT = "#b9b8b1"  # de-emphasised points (context, not the story)
WASH = "#ecebe6"  # scatter band about a fit

# Sequential blue ramp (ordered quantities) and the diverging neutral.
BLUE_150 = "#b7d3f6"
BLUE_300 = "#6da7ec"
BLUE_550 = "#1c5cab"
NEUTRAL = "#f0efec"

INSTRUMENT_STYLE = {
    "4inch": {"color": BLUE, "marker": "o", "label": "Sobral 4-inch (1919)"},
    "principe": {"color": ORANGE, "marker": "s", "label": "Pr\u00edncipe (1919)"},
    "astro": {"color": AQUA, "marker": "D", "label": "Sobral astrographic (1919)"},
}
GR_LINE = {"color": INK, "linewidth": 1.6, "linestyle": "-"}
NEWTON_LINE = {"color": INK_SECONDARY, "linewidth": 1.3, "linestyle": (0, (5, 3))}
ANCHOR_STYLE = {"color": INK, "marker": "*", "ms": 15, "mec": SURFACE, "mew": 1.0}


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
        }
    )


def annotate(ax, text: str, xy, xytext, *, color: str = INK_SECONDARY, fontsize: float = 8, **kwargs):
    """A callout with a thin leader line, in the project's annotation style."""
    return ax.annotate(
        text, xy=xy, xytext=xytext, color=color, fontsize=fontsize,
        arrowprops={"arrowstyle": "-", "color": MUTED, "linewidth": 0.7, "shrinkA": 1, "shrinkB": 3},
        **kwargs,
    )


#: caption characters per inch of width at 7.5 pt, which leaves a margin on either side
CAPTION_CHARS_PER_INCH = 16
#: space between the lowest artist and the caption's first line (inches)
CAPTION_GAP_INCH = 0.15


def caption(fig, text: str, y: float | None = None) -> None:
    """Footnote under a figure: data provenance, caveats, the honest small print. Wrapped here,
    not with wrap=True (which re-wraps to the cropped width on a tight save and leaves the PNG
    caption no margin); save() re-wraps it to the width of the rest of the figure and centres it
    under that, its top at ``y`` (figure fraction) or, by default, just below everything else.
    No-break spaces keep a number with its error ("0.96 +- 0.19") or its power of ten on one line."""
    nb, pm, times = "\u00a0", "\u00b1", "\u00d7"
    text = text.replace(f" {pm} ", f"{nb}{pm}{nb}").replace(f" {times} ", f"{nb}{times}{nb}").replace(
        f"{times} ", f"{times}{nb}")
    t = fig.text(0.5, -0.02 if y is None else y, textwrap.fill(text, int(fig.get_figwidth() * CAPTION_CHARS_PER_INCH)),
                 ha="center", va="top", fontsize=7.5, color=MUTED)
    fig._plotstyle_caption = (t, text, y is None)


def _fit_caption(fig) -> None:
    """Wrap the caption to the width of everything else as a tight save will crop it (tick labels
    can reach past the figure's edges), centre it under that and, unless caption() was given a
    height, put it CAPTION_GAP_INCH below the lowest artist."""
    cap = getattr(fig, "_plotstyle_caption", None)
    if cap is None:
        return
    t, text, auto_y = cap
    t.set_visible(False)
    dpi, save_dpi = fig.dpi, mpl.rcParams["savefig.dpi"]
    fig.set_dpi(dpi if save_dpi == "figure" else save_dpi)     # measure text as the PNG draws it
    bb = fig.get_tightbbox(fig.canvas.get_renderer())          # inches, without the caption
    fig.set_dpi(dpi)
    t.set_visible(True)
    t.set_text(textwrap.fill(text, int(bb.width * CAPTION_CHARS_PER_INCH)))
    t.set_x((bb.x0 + bb.x1) / 2 / fig.get_figwidth())
    if auto_y:
        t.set_y((bb.y0 - CAPTION_GAP_INCH) / fig.get_figheight())


def save(fig, path_stem, *, formats=("png", "svg")) -> list:
    """Save a figure to every requested format, returning the paths written."""
    _fit_caption(fig)
    written = []
    for fmt in formats:
        out = f"{path_stem}.{fmt}"
        fig.savefig(out, format=fmt)
        written.append(out)
    plt.close(fig)
    return written
