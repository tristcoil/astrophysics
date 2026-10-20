"""A periodic table of origins, coloured by what GALAH's stars show.

In the style of Johnson (2019, Science 363, 474), but every coloured cell is
a measurement from this project: the share of the element that arrived on
the prompt (core-collapse) timescale versus the delayed (Type Ia / AGB)
timescale in solar-like disc stars, from the two-process model. Hydrogen and
helium are labelled with the Big Bang's share (project K1's subject), which
GALAH cannot measure. Everything else is grey: not measured here.
"""

from __future__ import annotations

import numpy as np

SYMBOLS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
           "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce "
           "Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "
           "Fr Ra Ac Th Pa U").split()
#: No long-lived isotope: made in stars but decayed away.
UNSTABLE = {"Tc", "Pm", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Pa"}


def position(z: int):
    """(row, column) in an 18-column table; lanthanides (57-71) and actinides
    (89-103) on two rows below, offset by half a row."""
    if z == 1:
        return 0, 0
    if z == 2:
        return 0, 17
    if z <= 10:
        return 1, (z - 3) if z <= 4 else (z - 5 + 12)
    if z <= 18:
        return 2, (z - 11) if z <= 12 else (z - 13 + 12)
    if z <= 36:
        return 3, z - 19
    if z <= 54:
        return 4, z - 37
    if 57 <= z <= 71:
        return 7.4, z - 57 + 2
    if 89 <= z <= 103:
        return 8.4, z - 89 + 2
    if z <= 56:
        return 5, z - 55
    if z <= 86:
        return 5, z - 86 + 17
    return 6, z - 87 if z <= 88 else z - 104 + 3


def draw(ax, values: dict, *, compare: dict | None = None, uncertain: set | None = None,
         prompt_color="#eb6834", delayed_color="#2a78d6", bigbang_color="#4a3aa7",
         grey="#e1e0d9", ink="#0b0b0b", muted="#898781", bigbang: dict | None = None):
    """Draw the table on ``ax``.

    values: element -> prompt fraction (0-1; clipped) for measured elements.
    compare: element -> published prompt (or r-process) fraction, drawn as a tick.
    uncertain: elements whose value differs a lot between stellar classes (hatched).
    bigbang: element -> Big Bang share, for H and He.
    """
    from matplotlib.patches import Rectangle

    compare = compare or {}
    uncertain = uncertain or set()
    bigbang = bigbang or {}
    w, h = 0.94, 0.94
    for z, sym in enumerate(SYMBOLS, start=1):
        r, c = position(z)
        x0, y0 = c, -r
        if sym in bigbang:
            f = bigbang[sym]
            ax.add_patch(Rectangle((x0, y0 - h), w, h * f, color=bigbang_color, lw=0))
            ax.add_patch(Rectangle((x0, y0 - h + h * f), w, h * (1 - f), color=grey, lw=0))
            txt = "white"
        elif sym in values and np.isfinite(values[sym]):
            f = float(np.clip(values[sym], 0, 1))
            ax.add_patch(Rectangle((x0, y0 - h), w, h * f, color=prompt_color, lw=0))
            ax.add_patch(Rectangle((x0, y0 - h + h * f), w, h * (1 - f), color=delayed_color, lw=0))
            if sym in uncertain:
                ax.add_patch(Rectangle((x0, y0 - h), w, h, fill=False, hatch="////",
                                       edgecolor="white", lw=0))
            if sym in compare and np.isfinite(compare[sym]):
                yc = y0 - h + h * float(np.clip(compare[sym], 0, 1))
                ax.plot([x0 + 0.06, x0 + w - 0.06], [yc, yc], color=ink, lw=1.6,
                        solid_capstyle="butt")
            txt = "white"
        else:
            ax.add_patch(Rectangle((x0, y0 - h), w, h, color=grey, lw=0))
            txt = muted if sym not in UNSTABLE else "#b5b3ab"
        ax.text(x0 + 0.08, y0 - 0.1, str(z), ha="left", va="top", fontsize=5.2, color=txt)
        ax.text(x0 + w / 2, y0 - 0.52, sym, ha="center", va="center", fontsize=10.5,
                fontweight="bold", color=txt if txt != "white" else "white")
        if sym in values and np.isfinite(values.get(sym, np.nan)):
            v = values[sym]
            ax.text(x0 + w / 2, y0 - 0.84, ("≥100%" if v > 1.0 + 1e-9 else f"{100 * v:.0f}%"),
                    ha="center", va="center", fontsize=5.6, color="white")
    ax.text(2.47, -5.48, "57–71", ha="center", va="center", fontsize=7, color=muted)
    ax.text(2.47, -6.48, "89–103", ha="center", va="center", fontsize=7, color=muted)
    ax.set_xlim(-0.2, 18.1)
    ax.set_ylim(-9.5, 0.2)
    ax.set_aspect("equal")
    ax.axis("off")
