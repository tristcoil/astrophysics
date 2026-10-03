"""Shared matplotlib styling for the project's figures.

The categorical colours are the repository's validated colourblind-safe set:
adjacent slots keep a CVD colour difference of dE >= 8 (OKLab x100), and the
first three slots pass the check for every pair, so scatter plots never carry
more than three hues. Density maps use a single-hue ramp (light -> dark); maps
of a signed quantity (a contrast, a velocity offset) use the blue <-> red
diverging scale through a neutral grey, built as in project D3.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

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
C_MASK = "#e1e0d9"

# Semantic roles, named once so the same thing is the same colour everywhere.
C_OURS = BLUE          # our measurement / our members
C_REFERENCE = ORANGE   # a published value we compare against
C_THIRD = AQUA         # a third series when needed (another reference)
C_CONTEXT = MUTED      # field stars, background, context
C_MODEL = ORANGE       # a model curve drawn against data

BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
             "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
DIVERGING_MID = "#f0efec"
#: Sequential blue for densities of the thing we are measuring.
CMAP_BLUE = LinearSegmentedColormap.from_list("ramp_blue", [SURFACE] + BLUE_RAMP)
#: Neutral ramp for context densities (the field), drawn so it recedes.
CMAP_CONTEXT = LinearSegmentedColormap.from_list(
    "ramp_context", [SURFACE, "#eceae4", "#c3c2b7", "#a3a198", "#77756f", "#52514e", "#2b2a27"])


def use_style() -> None:
    """Apply the project-wide rcParams. Call once at import time in figure code."""
    mpl.rcParams.update({
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
        "image.cmap": "ramp_blue",
    })
    for cm in (CMAP_BLUE, CMAP_CONTEXT):
        if cm.name not in mpl.colormaps:
            mpl.colormaps.register(cm)


def annotate(ax, text: str, xy, xytext, *, color: str = INK_SECONDARY, **kwargs):
    """A callout with a thin leader line, in the project's annotation style."""
    return ax.annotate(text, xy=xy, xytext=xytext, color=color, fontsize=8,
                       arrowprops={"arrowstyle": "-", "color": color, "linewidth": 0.8,
                                   "shrinkA": 2, "shrinkB": 4}, **kwargs)


def caption(fig, text: str) -> None:
    """Footnote under a figure: data provenance, caveats, the honest small print."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    lowest = min(ax.get_tightbbox(renderer).y0 for ax in fig.axes)
    fig.text(0.5, lowest / fig.bbox.height - 0.012, text, ha="center", va="top", fontsize=7.5,
             color=MUTED, wrap=True)


#: Rasterised layers (star clouds, density maps) inside SVGs are stored at this
#: resolution: the SVGs exist for their vector text and lines.
SVG_RASTER_DPI = 150


def save(fig, path_stem, *, formats=("png", "svg")) -> list:
    written = []
    for fmt in formats:
        out = f"{path_stem}.{fmt}"
        fig.savefig(out, format=fmt, **({"dpi": SVG_RASTER_DPI} if fmt == "svg" else {}))
        written.append(out)
    plt.close(fig)
    return written


# --------------------------------------------------------------------------
# Diverging scale (from project D3): blue arm = the palette's blue ramp, red
# arm matched to it step for step in OKLab lightness.
# --------------------------------------------------------------------------


def _srgb_to_linear(c):
    import numpy as np
    c = np.asarray(c, float)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _linear_to_srgb(c):
    import numpy as np
    c = np.clip(np.asarray(c, float), 0.0, 1.0)
    return np.where(c <= 0.0031308, 12.92 * c, 1.055 * c ** (1 / 2.4) - 0.055)


def hex_to_oklab(hex_color: str):
    import numpy as np
    rgb = _srgb_to_linear([int(hex_color[i:i + 2], 16) / 255 for i in (1, 3, 5)])
    m1 = np.array([[0.4122214708, 0.5363325363, 0.0514459929],
                   [0.2119034982, 0.6806995451, 0.1073969566],
                   [0.0883024619, 0.2817188376, 0.6299787005]])
    m2 = np.array([[0.2104542553, 0.7936177850, -0.0040720468],
                   [1.9779984951, -2.4285922050, 0.4505937099],
                   [0.0259040371, 0.7827717662, -0.8086757660]])
    return m2 @ np.cbrt(m1 @ rgb)


def oklab_to_rgb(lab):
    import numpy as np
    m2i = np.array([[1.0, 0.3963377774, 0.2158037573],
                    [1.0, -0.1055613458, -0.0638541728],
                    [1.0, -0.0894841775, -1.2914855480]])
    m1i = np.array([[4.0767416621, -3.3077115913, 0.2309699292],
                    [-1.2684380046, 2.6097574011, -0.3413193965],
                    [-0.0041960863, -0.7034186147, 1.7076147010]])
    lin = m1i @ (m2i @ np.asarray(lab, float)) ** 3
    return _linear_to_srgb(lin), bool(np.all((lin >= -1e-4) & (lin <= 1 + 1e-4)))


def _mirror_ramp(ramp, target_hex):
    import numpy as np
    t = hex_to_oklab(target_hex)
    hue = np.arctan2(t[2], t[1])
    out = []
    for h in ramp:
        lab = hex_to_oklab(h)
        chroma = float(np.hypot(lab[1], lab[2]))
        for _ in range(60):
            rgb, ok = oklab_to_rgb([lab[0], chroma * np.cos(hue), chroma * np.sin(hue)])
            if ok:
                break
            chroma *= 0.97
        out.append(tuple(float(v) for v in rgb))
    return out


def diverging_cmap(name: str = "h3_diverging"):
    """Blue (negative) <-> grey (zero) <-> red (positive)."""
    from matplotlib.colors import to_rgb

    blue = [to_rgb(h) for h in BLUE_RAMP]
    red = _mirror_ramp(BLUE_RAMP, RED)
    colors = list(reversed(blue)) + [to_rgb(DIVERGING_MID)] + red
    cmap = LinearSegmentedColormap.from_list(name, colors, N=512)
    cmap.set_bad(C_MASK)
    return cmap
