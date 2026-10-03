"""Shared matplotlib styling for the project's figures.

The categorical colours are a colourblind-safe set: adjacent slots keep a CVD
colour difference of dE >= 8 (OKLab x100), so the figures stay readable for
deuteranopia/protanopia and in greyscale print.
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

SERIES = [BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED]

# --- Chrome and ink ---
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

# Semantic roles used across the CMB figures. Naming them once keeps the same
# quantity the same colour in every panel.
C_RAW = MUTED  # unbinned points, naive spectra: recessive context
C_DATA = BLUE  # the data being read (Planck's binned spectrum, FIRAS points)
C_OURS = BLUE  # a measurement this project made
C_MODEL = ORANGE  # a model curve (LCDM best fit, blackbody)
C_PUBLISHED = VIOLET  # a literature value we compare against
C_NOISE = YELLOW  # instrument noise (always direct-labelled: sub-3:1 contrast)
C_MASK = MUTED  # masked sky
C_ALT = AQUA  # an alternative hypothesis shown for contrast (direct-labelled)

# Neutral midpoint of the diverging temperature scale.
DIVERGING_MID = "#f0efec"

# The palette's blue ramp, steps 100 -> 700 (light -> dark).
BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
             "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]


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
            # Only the top/right spines go; a chart still needs its two rules.
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


def annotate(ax, text: str, xy, xytext, *, color: str = INK_SECONDARY, **kwargs):
    """A callout with a thin leader line, in the project's annotation style."""
    return ax.annotate(
        text,
        xy=xy,
        xytext=xytext,
        color=color,
        fontsize=8,
        arrowprops={
            "arrowstyle": "-",
            "color": color,
            "linewidth": 0.8,
            "shrinkA": 2,
            "shrinkB": 4,
        },
        **kwargs,
    )


def caption(fig, text: str) -> None:
    """Footnote under a figure: data provenance, caveats, the honest small print."""
    fig.text(
        0.5,
        -0.02,
        text,
        ha="center",
        va="top",
        fontsize=7.5,
        color=MUTED,
        wrap=True,
    )


def save(fig, path_stem, *, formats=("png", "svg")) -> list:
    """Save a figure to every requested format, returning the paths written."""
    written = []
    for fmt in formats:
        out = f"{path_stem}.{fmt}"
        fig.savefig(out, format=fmt)
        written.append(out)
    plt.close(fig)
    return written


# --------------------------------------------------------------------------
# Colour maths (for building the diverging temperature scale)
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
    """sRGB hex -> OKLab (Ottosson 2020), the space the palette validator uses."""
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
    """OKLab -> linear-light-clipped sRGB triple in [0, 1], plus an in-gamut flag."""
    import numpy as np
    m2i = np.array([[1.0, 0.3963377774, 0.2158037573],
                    [1.0, -0.1055613458, -0.0638541728],
                    [1.0, -0.0894841775, -1.2914855480]])
    m1i = np.array([[4.0767416621, -3.3077115913, 0.2309699292],
                    [-1.2684380046, 2.6097574011, -0.3413193965],
                    [-0.0041960863, -0.7034186147, 1.7076147010]])
    lin = m1i @ (m2i @ np.asarray(lab, float)) ** 3
    return _linear_to_srgb(lin), bool(np.all((lin >= -1e-4) & (lin <= 1 + 1e-4)))


def _mirror_ramp(ramp: list[str], target_hex: str) -> list[tuple]:
    """A ramp in the hue of ``target_hex`` with the same OKLab lightness steps as ``ramp``.

    Each step keeps the source step's lightness and chroma (reduced only as
    far as needed to stay inside the sRGB gamut), so the two arms of a
    diverging scale are equally loud at equal distances from the midpoint.
    """
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


def diverging_cmap(name: str = "cmb_diverging"):
    """Blue <-> red temperature scale through a neutral grey midpoint.

    Cold is blue and hot is red, matching the palette's documented diverging
    pair. The blue arm is the palette's blue ramp; the red arm is computed to
    match it step for step in OKLab lightness, so neither sign looks louder.
    Lightness falls monotonically from the midpoint along both arms, which
    ``check_diverging`` verifies.
    """
    from matplotlib.colors import LinearSegmentedColormap, to_rgb

    blue = [to_rgb(h) for h in BLUE_RAMP]  # light -> dark
    red = _mirror_ramp(BLUE_RAMP, RED)  # light -> dark
    colors = list(reversed(blue)) + [to_rgb(DIVERGING_MID)] + red
    cmap = LinearSegmentedColormap.from_list(name, colors, N=512)
    cmap.set_bad(C_MASK)
    return cmap


def check_diverging(cmap, n: int = 257) -> dict:
    """Lightness must rise monotonically toward the midpoint from both ends."""
    import numpy as np
    from matplotlib.colors import to_hex

    ls = np.array([hex_to_oklab(to_hex(cmap(x)))[0] for x in np.linspace(0, 1, n)])
    mid = n // 2
    left_ok = bool(np.all(np.diff(ls[: mid + 1]) > -1e-3))
    right_ok = bool(np.all(np.diff(ls[mid:]) < 1e-3))
    return {"monotone": left_ok and right_ok, "L_ends": (float(ls[0]), float(ls[-1])),
            "L_mid": float(ls[mid])}


# --------------------------------------------------------------------------
# Axes for power spectra
# --------------------------------------------------------------------------


def _hybrid_functions(lmin: float, lmax: float, split: float, log_fraction: float):
    """Forward/inverse maps for an axis that is logarithmic below ``split``, linear above."""
    import numpy as np

    a = log_fraction / np.log10(split / lmin)
    b = (1.0 - log_fraction) / (lmax - split)

    def forward(x):
        x = np.maximum(np.asarray(x, float), lmin / 10)
        return np.where(x <= split, a * np.log10(x / lmin), log_fraction + b * (x - split))

    def inverse(y):
        y = np.asarray(y, float)
        return np.where(y <= log_fraction, lmin * 10 ** (y / a), split + (y - log_fraction) / b)

    return forward, inverse


def set_ell_axis(ax, *, lmin: float = 2.0, lmax: float = 2510.0, split: float = 30.0,
                 log_fraction: float = 0.18) -> None:
    """Planck's hybrid multipole axis: logarithmic below l = 30, linear above.

    The largest scales (a handful of multipoles) and the acoustic peaks
    (thousands) cannot share a linear axis without the former vanishing, nor a
    logarithmic one without the peaks crowding into the right edge. The Planck
    team's convention squeezes the first decade into ~a fifth of the width.
    ``log_fraction=0`` gives a plain linear axis from ``lmin``.
    """
    from matplotlib.ticker import FixedLocator, NullFormatter

    if log_fraction <= 0:
        ax.set_xlim(lmin, lmax)
        ax.xaxis.set_major_locator(FixedLocator([500, 1000, 1500, 2000, 2500]))
        ax._ell_scale = None
        return
    fns = _hybrid_functions(lmin, lmax, split, log_fraction)
    ax.set_xscale("function", functions=fns)
    ax._ell_scale = fns
    ax.set_xlim(lmin, lmax)
    ax.xaxis.set_major_locator(FixedLocator([2, 10, 30, 500, 1000, 1500, 2000, 2500]))
    ax.xaxis.set_minor_locator(FixedLocator([3, 4, 5, 6, 7, 8, 9, 20, 250, 750, 1250, 1750,
                                             2250]))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.axvline(split, color=AXIS, lw=0.8, zorder=0)


def add_degree_axis(ax, ticks=(90, 10, 1, 0.5, 0.2, 0.1), label="Angular scale") -> None:
    """Top axis in degrees: multipole l corresponds to about 180 deg / l.

    Built as a twin axis carrying the parent's exact scale, with ticks placed at
    l = 180/theta. (matplotlib's secondary_xaxis does not compose correctly with
    a custom FuncScale parent -- it put "1 deg" near l = 5 instead of l = 180.)
    """
    from matplotlib.ticker import FixedLocator, NullLocator

    top = ax.twiny()
    fns = getattr(ax, "_ell_scale", None)
    if fns is not None:
        top.set_xscale("function", functions=fns)
    lo, hi = ax.get_xlim()
    top.set_xlim(lo, hi)
    keep = [t for t in ticks if lo <= 180.0 / t <= hi]
    top.xaxis.set_major_locator(FixedLocator([180.0 / t for t in keep]))
    top.xaxis.set_minor_locator(NullLocator())
    top.set_xticklabels([f"{t:g}°" for t in keep])
    top.tick_params(axis="x", colors=MUTED, labelsize=8)
    top.grid(False)
    for spine in top.spines.values():
        spine.set_visible(False)
    top.set_xlabel(label, color=INK_SECONDARY, fontsize=8.5, labelpad=6)
    ax.callbacks.connect("xlim_changed", lambda a: top.set_xlim(a.get_xlim()))
