"""Shared matplotlib styling: the atlas's validated palette, assigned by entity.

Colour follows the *thing* in every figure:

* this project's own calculation          -> ink (the curve every figure is about)
* the textbook ideal-gas picture          -> blue
* the main published reference (Saikawa & Shirai 2018; Escudero 2020; HotQCD) -> orange
* a second published reference (Laine-Meyer-Schroeder 2015; WB 2014)        -> aqua
* anything else                           -> neutral grey with its own marker shape

Blue, orange and aqua pass the dataviz validator for all pairs in light mode (worst CVD
separation 9.2, normal vision 24.0); aqua is under 3:1 contrast on the light surface, so every
curve carries a direct label. Stacked areas use the eight-slot order on adjacent pairs.
"""

from __future__ import annotations

import matplotlib as mpl
import numpy as np
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
SLOTS = [BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED]

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
WASH = "#f0efec"

C_OURS = INK
C_IDEAL = BLUE
C_REF1 = ORANGE
C_REF2 = AQUA
C_OTHER = MUTED


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


def human_time(t_s: float) -> str:
    """'19.5 µs', '1.4 s', '1.8 min', '51,300 yr' ..."""
    year = 365.25 * 86400
    if t_s < 1e-12:
        return f"{t_s * 1e15:.3g} fs"
    if t_s < 1e-9:
        return f"{t_s * 1e12:.3g} ps"
    if t_s < 1e-6:
        return f"{t_s * 1e9:.3g} ns"
    if t_s < 1e-3:
        return f"{t_s * 1e6:.3g} µs"
    if t_s < 1:
        return f"{t_s * 1e3:.3g} ms"
    if t_s < 120:
        return f"{t_s:.3g} s"
    if t_s < 7200:
        return f"{t_s / 60:.3g} min"
    if t_s < 2 * 86400:
        return f"{t_s / 3600:.3g} h"
    if t_s < 2 * year:
        return f"{t_s / 86400:.3g} days"
    if t_s < 1e6 * year:
        return f"{t_s / year:,.0f} yr"
    if t_s < 1e9 * year:
        return f"{t_s / year / 1e6:.3g} Myr"
    return f"{t_s / year / 1e9:.4g} Gyr"


def human_energy(T_gev: float) -> str:
    if T_gev >= 1e3:
        return f"{T_gev / 1e3:.3g} TeV"
    if T_gev >= 1:
        return f"{T_gev:.3g} GeV"
    if T_gev >= 1e-3:
        return f"{T_gev * 1e3:.3g} MeV"
    if T_gev >= 1e-6:
        return f"{T_gev * 1e6:.3g} keV"
    if T_gev >= 1e-11:
        return f"{T_gev * 1e9:.3g} eV"
    return f"{T_gev * 1e12:.3g} meV"


SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def sci(x: float, digits: int = 2) -> str:
    """1.9e-06 -> '1.9 × 10⁻⁶'."""
    if x == 0:
        return "0"
    e = int(np.floor(np.log10(abs(x))))
    m = x / 10 ** e
    if -2 <= e <= 3:
        return f"{x:.{digits}g}"
    return f"{m:.{digits}g} × 10{str(e).translate(SUPERSCRIPT)}"
