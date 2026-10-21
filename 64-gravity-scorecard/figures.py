"""Every figure of the project, drawn from cached results (run.py calls them in turn)."""

from __future__ import annotations

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.lines import Line2D                    # noqa: E402
from matplotlib.patches import FancyArrowPatch, Ellipse  # noqa: E402

import config as C                                     # noqa: E402
import gr_predictions as G                             # noqa: E402
import literature as L                                 # noqa: E402
import plotstyle as S                                  # noqa: E402

S.use_style()
OUT = C.FIGURE_DIR
DARK = {S.ORANGE: "#a3441f", S.AQUA: "#137a55", S.BLUE: S.BLUE_550, S.VIOLET: S.VIOLET}   # text versions


def _family_handles(extra_atlas=True):
    hs = [Line2D([], [], ls="", marker=st["marker"], color=st["color"], mec=S.SURFACE, mew=0.6, ms=7,
                 label=st["label"]) for st in S.FAMILY_STYLE.values()]
    if extra_atlas:
        hs.append(Line2D([], [], ls="", marker="*", color=S.INK, mec=S.SURFACE, mew=0.6, ms=11,
                         label="measured in this atlas"))
    return hs


def _point(ax, x, y, row, **kw):
    if row["atlas"]:
        ax.plot(x, y, **S.ATLAS_STYLE, zorder=5)
    else:
        st = S.FAMILY_STYLE[row["family"]]
        ax.plot(x, y, ls="", marker=st["marker"], color=st["color"], ms=7, mec=S.SURFACE, mew=0.7, zorder=4)


# ---------------------------------------------------------------------------
# 1. The scorecard: precision against the strength of gravity
# ---------------------------------------------------------------------------

#: label placement (dx points, dy points, ha); default is to the right
LABEL_POS = {
    # top panel (1e-3 .. 1)
    "antimatter": (-7, 0, "right"), "GP-B frame dragging": (7, 0, "left"), "N3 J2 term": (7, 0, "left"),
    "Pound–Snider": (-7, 0, "right"), "GP-B geodetic": (7, 0, "left"), "LARES 2": (-7, 0, "right"),
    "M3 wide binaries": (7, 0, "left"), "LLR geodetic": (7, 0, "left"), "N7 LLR Shapiro": (7, 0, "left"),
    "Sun redshift": (-7, 0, "right"), "N3 Sun": (-7, -6, "right"), "Hulse–Taylor Ṗb": (-7, 0, "right"),
    "galaxy γ": (-7, 6, "right"), "J0348+0432": (-7, 0, "right"), "spin precession": (-7, 0, "right"),
    "E_G": (7, -2, "left"), "growth index": (7, -7, "left"), "Σ0 (DESI)": (7, 6, "left"),
    "μ0 (DESI)": (7, 2, "left"), "D10 ISW": (7, 0, "left"), "N3 Sirius B": (7, 0, "left"),
    "S2 redshift": (7, 0, "left"), "S2 precession": (-7, 0, "right"), "H1 S2 redshift": (7, 0, "left"),
    "GW250114 inspiral": (-7, 0, "right"), "IMR consistency": (-7, 0, "right"),
    "GW250114 ringdown": (-7, 0, "right"), "ringdown damping": (-7, 0, "right"),
    "M87* shadow": (-7, 3, "right"), "Sgr A* shadow": (-7, -2, "right"), "I5 Sgr A*": (7, -3, "left"),
    "I5 M87*": (7, 3, "left"),
    # bottom panel (1e-16 .. 1e-3)
    "MICROSCOPE": (7, 0, "left"), "Eöt-Wash": (-7, 0, "right"), "atoms": (-7, 0, "right"),
    "LLR Earth–Moon": (7, -5, "left"), "N7 LLR (this work)": (7, 5, "left"), "triple pulsar": (7, 0, "left"),
    "Cassini": (7, 0, "left"), "VLBI deflection": (7, 0, "left"), "Mercury perihelion": (7, 0, "left"),
    "Galileo clocks": (-7, -3, "right"), "N3 Galileo": (7, 2, "left"), "GP-A": (-7, 0, "right"),
    "double pulsar Ṗb": (7, 0, "left"), "double pulsar s": (7, 0, "left"),
}
#: display-only sideways spread for tests that share a field strength
SPREAD = {"Eöt-Wash": 0.80, "atoms": 0.80, "antimatter": 0.80, "Pound–Snider": 0.80, "GW250114 ringdown": 0.80,
          "ringdown damping": 0.80, "I5 Sgr A*": 1.25, "I5 M87*": 1.25, "E_G": 0.90, "Σ0 (DESI)": 1.12,
          "μ0 (DESI)": 1.12, "growth index": 0.90, "D10 ISW": 1.25, "N7 LLR (this work)": 1.25,
          "N7 LLR Shapiro": 1.25, "N3 Galileo": 1.25, "H1 S2 redshift": 1.15, "S2 redshift": 0.87, "S2 precession": 0.87}
SPLIT = 1e-3


def fig_scorecard(table: pd.DataFrame):
    fig = plt.figure(figsize=(12.0, 9.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.15, 1.0], hspace=0.05)
    ax_t = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1], sharex=ax_t)
    for _, r in table.iterrows():
        x = r["x"] * SPREAD.get(r["short"], 1.0)
        y = r["sigma"]
        ax = ax_t if y >= SPLIT else ax_b
        _point(ax, x, y, r)
        dx, dy, ha = LABEL_POS.get(r["short"], (7, 0, "left"))
        ax.annotate(r["short"], (x, y), xytext=(dx, dy), textcoords="offset points", ha=ha, va="center",
                    fontsize=7.4, color=S.INK if r["atlas"] else S.INK_SECONDARY,
                    fontweight="bold" if r["atlas"] else "normal")
    for ax in (ax_t, ax_b):
        ax.set_xscale("log")
        ax.set_yscale("log")
    ax_t.set_xlim(1e-12, 3)
    ax_t.set_ylim(SPLIT, 1.0)
    ax_b.set_ylim(1e-16, SPLIT)
    ax_t.spines["bottom"].set_visible(False)
    ax_t.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    # break marks
    for ax, y in ((ax_t, 0.0), (ax_b, 1.0)):
        ax.plot([-0.012, 0.012], [y - 0.012, y + 0.012], transform=ax.transAxes, color=S.AXIS, lw=0.9, clip_on=False)
    regimes = [(1.4e-12, "two suns\n5,000 AU apart"), (1.3e-10, "Earth and\nits satellites"),
               (4e-9, "the Moon,\nthe planets"), (1.1e-6, "the Sun's edge,\nbinary pulsars"),
               (1.6e-5, "galaxies,\ncosmology"), (1.9e-4, "S2,\nSirius B"), (0.07, "black\nholes")]
    for xr, txt in regimes:
        ax_b.text(xr, 1.4e-16, txt, fontsize=7.4, color=S.MUTED, ha="left", va="bottom")
    ax_b.set_xlabel("Strength of gravity probed, GM / (r c²)   (weak  →  strong)")
    fig.text(0.065, 0.5, "How precisely each test agrees with general relativity (1σ, as a fraction of GR's effect,\n"
             "or the violation GR forbids)   —   lower is more precise", rotation=90, ha="center", va="center",
             fontsize=9, color=S.INK_SECONDARY)
    ax_t.set_title("Is Einstein still right? A scorecard of tests of gravity")
    ax_t.text(1.0, 1.01, "upper panel: 10⁻³ to 1 · lower panel: 10⁻¹⁶ to 10⁻³ (note the break)", transform=ax_t.transAxes,
              ha="right", va="bottom", fontsize=7.4, color=S.MUTED)
    ax_b.legend(handles=_family_handles(), loc="lower right", bbox_to_anchor=(1.0, 0.10), fontsize=8)
    S.caption(fig, "Each point is one test: the 1σ uncertainty of its fractional deviation from GR (or of a violation parameter "
              "GR sets to zero), placed at GM/(rc²) of the field it happens in. Published results were checked against the "
              "primary papers (literature/checks/*.csv); stars are measurements made in this atlas (N3, H1, I5, M3, D10) and "
              "in this project (N7). Tests that share a field strength are spread slightly sideways. All agree with GR within "
              "2σ except the growth of cosmic structure (3.4σ); see the report card.", y=0.035)
    return S.save(fig, OUT / "01_scorecard")


# ---------------------------------------------------------------------------
# 2. The report card: deviation in units of its own error
# ---------------------------------------------------------------------------

def fig_report_card(table: pd.DataFrame):
    order = {"equivalence": 0, "weak": 1, "strong": 2, "cosmology": 3}
    t = table.assign(o=table["family"].map(order)).sort_values(["o", "x"]).reset_index(drop=True)
    n = len(t)
    fig = plt.figure(figsize=(10.0, 0.24 * n + 3.4))
    ax = fig.add_axes([0.27, 0.19, 0.42, 0.74])
    ax.axvspan(-2, 2, color=S.WASH, lw=0)
    ax.axvspan(-1, 1, color=S.GRID, lw=0)
    ax.axvline(0, **S.GR_LINE)
    for i, r in t.iterrows():
        y = n - 1 - i
        p = float(np.clip(r["pull"], -4.4, 4.4))
        _point(ax, p, y, r)
        ax.text(-4.75, y, r["short"], ha="right", va="center", fontsize=7.6,
                color=S.INK if r["atlas"] else S.INK_SECONDARY, fontweight="bold" if r["atlas"] else "normal")
        ax.text(4.75, y, f"{r['pull']:+.1f}σ", ha="left", va="center", fontsize=7.2, color=S.MUTED)
    ax.set_xlim(-4.5, 4.5)
    ax.set_ylim(-0.8, n + 0.6)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("(measured − GR) / 1σ uncertainty")
    ax.text(0.04, n - 0.1, "GR", ha="left", va="bottom", fontsize=8, color=S.INK)
    fig.suptitle("The report card: how far each test lands from Einstein", x=0.05, ha="left", y=0.985,
                 fontsize=11, fontweight="bold")
    ax2 = fig.add_axes([0.27, 0.05, 0.42, 0.085])
    pulls = t["pull"].to_numpy()
    bins = np.arange(-4.5, 4.51, 0.5)
    ax2.hist(np.clip(pulls, -4.4, 4.4), bins=bins, color=S.BLUE_300, edgecolor=S.SURFACE, lw=1.2)
    xx = np.linspace(-4.5, 4.5, 300)
    ax2.plot(xx, len(pulls) * 0.5 * np.exp(-0.5 * xx**2) / np.sqrt(2 * np.pi), color=S.INK, lw=1.0)
    ax2.set_xlim(-4.5, 4.5)
    ax2.set_yticks([])
    ax2.spines["left"].set_visible(False)
    within1 = int(np.sum(np.abs(pulls) <= 1))
    within2 = int(np.sum(np.abs(pulls) <= 2))
    ax2.set_title(f"{within1} of {len(pulls)} within 1σ (if GR is right, expect {0.683 * len(pulls):.0f}); "
                  f"{within2} within 2σ (expect {0.954 * len(pulls):.0f}); curve: normal distribution",
                  fontsize=8, fontweight="normal", loc="left")
    fig.legend(handles=_family_handles(), loc="upper left", bbox_to_anchor=(0.77, 0.93), fontsize=7.8)
    fig.text(0.77, 0.70, "Darker band: ±1σ\nlighter band: ±2σ\n(values clipped at ±4.4σ)", fontsize=7.6,
             color=S.INK_SECONDARY, va="top")
    return S.save(fig, OUT / "02_report_card")


# ---------------------------------------------------------------------------
# 3. The weak equivalence principle, 1687-2022
# ---------------------------------------------------------------------------

def fig_equivalence_principle():
    t = L.wep_timeline()
    styles = {"lab": ("o", S.BLUE, "torsion balances, pendulums, free fall, MICROSCOPE"), "space": ("o", S.BLUE, None),
              "lab_doubtful": ("o", S.BLUE, None), "atom": ("D", S.AQUA, "atom interferometers"),
              "llr": ("s", S.ORANGE, "lunar laser ranging (the Earth and Moon as wholes)")}
    fig, ax = plt.subplots(figsize=(11.0, 6.4))
    for series, (mk, col, lab) in styles.items():
        d = t[t["series"] == series]
        for _, r in d.iterrows():
            hollow = series == "lab_doubtful"
            ms = 9 if series == "space" else 6.5
            ax.plot(r["year"], r["precision"], ls="", marker="v" if r["kind"] == "limit" else mk, ms=ms,
                    color=col, mfc=S.SURFACE if hollow else col, mec=col if hollow else S.SURFACE,
                    mew=1.0 if hollow else 0.6, zorder=4)
    lab = t[t["series"].isin(["lab", "space"])].sort_values("year")
    rec, best = [], np.inf
    for _, r in lab.iterrows():
        if r["precision"] < best:
            best = r["precision"]
            rec.append((r["year"], best))
    ry, rv = zip(*rec)
    ax.step(list(ry) + [2026], list(rv) + [rv[-1]], where="post", color=S.BLUE_300, lw=1.0, zorder=1)
    notes = {"newton1687": ("Newton, pendulums", (8, 0), "left"), "bessel1832": ("Bessel", (8, 0), "left"),
             "eotvos1890": ("Eötvös, torsion balance", (8, 0), "left"),
             "epf1922": ("Eötvös, Pekár & Fekete", (8, 4), "left"),
             "renner1935": ("Renner (errors\ntoo small; open)", (-8, -2), "right"),
             "rkd1964": ("Roll, Krotkov & Dicke", (-8, 0), "right"), "bp1972": ("Braginsky & Panov", (-8, 0), "right"),
             "niebauer1987": ("free fall", (-8, 0), "right"),
             "llr_shapiro1976": ("first lunar test", (8, -9), "left"),
             "wagner2012_BeAl_earth": ("Eöt-Wash 2012", (-14, -9), "right"),
             "asenbaum2020": ("atoms, Stanford 2020", (8, 0), "left"),
             "llr_biskupek2021": ("lunar ranging 2021", (8, 2), "left"),
             "microscope2017": ("MICROSCOPE 2017", (8, 0), "left"),
             "microscope2022": ("MICROSCOPE 2022 (satellite)", (8, 0), "left")}
    for rid, (txt, off, ha) in notes.items():
        r = t[t["id"] == rid].iloc[0]
        ax.annotate(txt, (r["year"], r["precision"]), xytext=off, textcoords="offset points", fontsize=7.6,
                    ha=ha, va="center", color=S.INK_SECONDARY)
    am = L.antimatter()
    ax.text(0.015, 0.035, f"Antimatter falls too: antihydrogen accelerates at {am['a_over_g']:.2f} ± {am['sigma']:.2f} g "
            "(ALPHA-g 2023)\n— down, not up, but only to ~20%: a test of gravity's sign, not of its precision.",
            transform=ax.transAxes, fontsize=7.8, color=S.INK_SECONDARY, va="bottom")
    ax.set_yscale("log")
    ax.set_xlim(1670, 2090)
    ax.set_xticks([1700, 1750, 1800, 1850, 1900, 1950, 2000])
    ax.set_ylim(5e-16, 5e-2)
    ax.set_xlabel("year of publication")
    ax.set_ylabel("precision of η = 2|a₁ − a₂| / |a₁ + a₂|   (1σ, or a stated limit ▼)")
    ax.set_title("Do all bodies fall alike? Three centuries of the equivalence principle")
    hs = [Line2D([], [], ls="", marker=m, color=c, mec=S.SURFACE, ms=7, label=lab) for m, c, lab in
          [(v[0], v[1], v[2]) for v in styles.values() if v[2]]]
    hs.append(Line2D([], [], ls="", marker="v", color=S.MUTED, ms=7, label="early limit (a sensitivity statement)"))
    hs.append(Line2D([], [], color=S.BLUE_300, lw=1.0, label="record of the material tests"))
    ax.legend(handles=hs, loc="upper right", fontsize=8)
    S.caption(fig, "Values from literature/checks/wep.csv, each checked against the original paper; 1σ uncertainties "
              "(statistical and systematic combined), 95% intervals divided by 1.96. Lunar ranging compares the Earth and the "
              "Moon as whole bodies (composition plus self-gravity), so it is shown as its own series.", y=-0.01)
    return S.save(fig, OUT / "03_equivalence_principle")


# ---------------------------------------------------------------------------
# 4. A century of light bending and time delay: (1+gamma)/2
# ---------------------------------------------------------------------------

GAMMA_METHODS = {
    "optical_eclipse": ("o", S.BLUE, "starlight at total eclipses"),
    "radio_interferometry": ("s", S.ORANGE, "radio interferometers (quasars passing the Sun)"),
    "vlbi": ("D", S.AQUA, "VLBI (continent-wide radio networks)"),
    "space_astrometry": ("^", S.VIOLET, "Hipparcos (optical, from space)"),
    "time_delay": ("P", S.INK, "radar and spacecraft time delay"),
}


def fig_gamma_history():
    t = L.gamma_timeline()
    t["m"] = t["method"].where(~t["method"].str.startswith("time_delay"), "time_delay")
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(10.8, 8.8), sharex=True, gridspec_kw={"height_ratios": [1.25, 1.0]})
    for key, (mk, col, lab) in GAMMA_METHODS.items():
        d = t[t["m"] == key]
        for _, r in d.iterrows():
            hollow = r["flag"] in ("secondary", "rereduction", "flagged")
            kw = dict(ls="", marker=mk, ms=6.5, color=col, mec=col if hollow else S.SURFACE,
                      mfc=S.SURFACE if hollow else col, mew=1.0 if hollow else 0.6)
            ax.errorbar(r["year"], r["value"], yerr=r["sigma"], elinewidth=0.9, capsize=0, **kw)
            bx.plot(r["year"], r["sigma"], **kw)
    ax.axhline(1.0, **S.GR_LINE)
    ax.axhline(0.5, color=S.MUTED, lw=0.8)
    ax.text(2024, 1.015, "Einstein (GR)", fontsize=7.8, ha="right", va="bottom", color=S.INK)
    ax.text(2024, 0.515, "half: Soldner's 1801 Newtonian value, Einstein's own 1911 value", fontsize=7.6,
            ha="right", va="bottom", color=S.MUTED)
    notes = [("E1919-sobral-4inch", "Eddington's 1919 eclipse\n(Sobral and Principe)", (-4, 30), "left"),
             ("E1919-sobral-combined-harvey1979", "the 1919 plates\nre-measured in 1979", (1923.0, 0.66), "data"),
             ("E1929-takengon", "1929", (8, 0), "left"), ("E1936-ussr-mikhailov-1959", "1936", (8, 0), "left"),
             ("E1947-bocaiuva", "1947", (8, 0), "left"),
             ("E2017-wyoming-bruns", "2017: a small telescope\nand a camera", (-6, 26), "right")]
    for rid, txt, off, ha in notes:
        r = t[t["id"] == rid].iloc[0]
        if ha == "data":
            ax.annotate(txt, (r["year"], r["value"]), xytext=off, textcoords="data", fontsize=7.4, ha="left",
                        color=S.INK_SECONDARY, arrowprops={"arrowstyle": "-", "color": S.MUTED, "lw": 0.6})
            continue
        ax.annotate(txt, (r["year"], r["value"]), xytext=off, textcoords="offset points", fontsize=7.4,
                    ha=ha, color=S.INK_SECONDARY,
                    arrowprops={"arrowstyle": "-", "color": S.MUTED, "lw": 0.6} if abs(off[1]) > 10 else None)
    ax.set_ylim(0.4, 1.95)
    ax.set_ylabel("(1 + γ)/2   (GR = 1)")
    ax.set_title("A century of bending light: how much does the Sun curve space?")
    bx.set_yscale("log")
    bx.set_ylim(5e-6, 0.6)
    bx.set_ylabel("1σ uncertainty of (1 + γ)/2")
    bx.set_xlabel("year of the observations")
    bnotes = [("T1976-viking", "Viking landers (Mars)", (8, -4), "left"), ("A1997-hipparcos", "Hipparcos", (8, 3), "left"),
              ("T2002-cassini", "Cassini", (8, 0), "left"), ("V2011-lambert", "VLBI 2011", (-8, 0), "right"),
              ("V2018-titov-global", "VLBI 2018: 3σ off GR; its authors\nsuspect systematics (open)", (2024.5, 9e-3), "data"),
              ("R1975-nrao-fomalont", "radio interferometers", (-8, -4), "right")]
    for rid, txt, off, ha in bnotes:
        r = t[t["id"] == rid].iloc[0]
        if ha == "data":
            bx.annotate(txt, (r["year"], r["sigma"]), xytext=off, textcoords="data", fontsize=7.4, ha="right",
                        va="center", color=S.INK_SECONDARY, arrowprops={"arrowstyle": "-", "color": S.MUTED, "lw": 0.6})
            continue
        bx.annotate(txt, (r["year"], r["sigma"]), xytext=off, textcoords="offset points", fontsize=7.4,
                    ha=ha, va="center", color=S.INK_SECONDARY,
                    arrowprops={"arrowstyle": "-", "color": S.MUTED, "lw": 0.6} if abs(off[1]) > 10 else None)
    hs = [Line2D([], [], ls="", marker=mk, color=col, mec=S.SURFACE, ms=7, label=lab) for mk, col, lab in GAMMA_METHODS.values()]
    hs.append(Line2D([], [], ls="", marker="o", color=S.MUTED, mfc=S.SURFACE, ms=7,
                     label="open: re-reduction, secondary source, or flagged"))
    ax.legend(handles=hs, loc="upper right", fontsize=7.6)
    bx.set_xlim(1912, 2026)
    S.caption(fig, "Values from literature/checks/gamma.csv (every row checked against the original paper where it could be "
              "reached). Errors converted to 1σ: 'probable errors' (the custom until about 1950) × 1.483; γ and γ − 1 halved to "
              "(1+γ)/2. Will (2014, Fig. 5) draws the early probable errors unconverted and quotes Lebach et al. (1995) at twice "
              "its published error; neither is repeated here.", y=-0.01)
    return S.save(fig, OUT / "04_gamma_history")


# ---------------------------------------------------------------------------
# 5. gamma and beta: which theories survive
# ---------------------------------------------------------------------------

def _band_lines(ax, x, lo, hi, color):
    ax.plot(x, lo, color=color, lw=0.7, alpha=0.7)
    ax.plot(x, hi, color=color, lw=0.7, alpha=0.7)


def fig_gamma_beta(llr_gamma: tuple[float, float] | None = None):
    p = L.ppn_constraints()
    u = 1e-5
    fig, ax = plt.subplots(figsize=(7.8, 7.2))
    xg = np.linspace(-12, 12, 400)
    yy = np.array([-15, 15])
    g, sg = (v / u for v in p["cassini_gamma_minus_1"])
    ax.axvspan(g - sg, g + sg, color=S.BLUE, alpha=0.28, lw=0)
    for k in (-2, 2):
        ax.plot([g + k * sg] * 2, yy, color=S.BLUE, lw=0.7, alpha=0.7)
    c, sc = (v / u for v in p["perihelion_2g_minus_b"])
    ax.fill_between(xg, 2 * xg - c - sc, 2 * xg - c + sc, color=S.ORANGE, alpha=0.26, lw=0)
    _band_lines(ax, xg, 2 * xg - c - 2 * sc, 2 * xg - c + 2 * sc, S.ORANGE)
    e, se = (v / u for v in p["nordtvedt_eta"])
    ax.fill_between(xg, (e + xg - se) / 4, (e + xg + se) / 4, color=S.AQUA, alpha=0.30, lw=0)
    _band_lines(ax, xg, (e + xg - 2 * se) / 4, (e + xg + 2 * se) / 4, S.AQUA)
    (pg, psg), (pb, psb) = p["epm2011"]
    ax.errorbar(pg / u, pb / u, xerr=psg / u, yerr=psb / u, fmt="D", color=S.VIOLET, ms=6, mec=S.SURFACE,
                elinewidth=1.0, capsize=0, zorder=5)
    ax.plot(0, 0, marker="*", color=S.INK, ms=16, mec=S.SURFACE, zorder=6)
    ax.annotate("general relativity\nγ = β = 1", (0, 0), xytext=(-6.5, 8.0), fontsize=8.5, color=S.INK,
                arrowprops={"arrowstyle": "-", "color": S.MUTED, "lw": 0.7})
    ax.text(g + 0.4, -13.2, "Cassini time delay\nγ − 1 = (2.1 ± 2.3)×10⁻⁵", fontsize=7.8, color=DARK[S.BLUE])
    ax.text(6.4, 9.5, "Mercury's perihelion\n(2 + 2γ − β)/3\n(MESSENGER)", fontsize=7.8, color=DARK[S.ORANGE])
    ax.text(-11.6, -4.9, "Nordtvedt effect in lunar ranging\nη = 4β − γ − 3 = (−0.2 ± 1.1)×10⁻⁴", fontsize=7.8,
            color=DARK[S.AQUA])
    ax.annotate("planetary ephemeris EPM2011\n(β and γ fitted together)", (pg / u, pb / u), xytext=(3.6, -10.0),
                fontsize=7.8, color=S.VIOLET, arrowprops={"arrowstyle": "-", "color": S.MUTED, "lw": 0.7})
    for om in (4e4, 1e5):
        ax.plot(-1 / (2 + om) / u, 0, marker="|", color=S.INK_SECONDARY, ms=10, zorder=5)
    ax.annotate("Brans–Dicke theory\nwith ω = 40,000 and 100,000", (-1 / (2 + 4e4) / u, 0), xytext=(-11.6, 2.4),
                fontsize=7.6, color=S.INK_SECONDARY, arrowprops={"arrowstyle": "-", "color": S.MUTED, "lw": 0.7})
    llr_note = ""
    if llr_gamma is not None:
        llr_note = (f" This project's lunar ranging alone gives γ − 1 = ({llr_gamma[0] / 1e-3:+.1f} ± "
                    f"{llr_gamma[1] / 1e-3:.1f})×10⁻³, far off this scale.")
    ax.set_xlim(-12, 12)
    ax.set_ylim(-15, 15)
    ax.set_xlabel("γ − 1   (×10⁻⁵): how much space a unit of mass curves")
    ax.set_ylabel("β − 1   (×10⁻⁵): how nonlinear gravity is")
    ax.set_title("Pinning down γ and β with three solar-system tests")
    S.caption(fig, "Shaded: 1σ; thin lines: 2σ. The perihelion band uses Park et al. (2017)'s β, fitted with Cassini's γ; "
              "the lunar band assumes the weak equivalence principle holds exactly (Hofmann & Müller 2018). Values from "
              "literature/checks/beta_lunar.csv and gamma.csv." + llr_note, y=-0.01)
    return S.save(fig, OUT / "05_gamma_beta")


# ---------------------------------------------------------------------------
# 6. Binary pulsars: each measured effect, through GR, is a line of possible masses
# ---------------------------------------------------------------------------

def _omdot_2pn(m1, m2, pb_days, e, omdot_lt_deg_yr=0.0):
    """Periastron advance with the 2PN correction of Damour & Schaefer (1988), plus a Lense-Thirring term."""
    pk = G.pk_parameters(m1, m2, pb_days, e)
    M = m1 + m2
    xa, xb = m1 / M, m2 / M
    beta2 = (G.T_SUN * M * 2 * np.pi / (pb_days * C.DAY)) ** (2 / 3)
    f_o = (39 / 4 * xa**2 + 27 / 4 * xb**2 + 15 * xa * xb) / (1 - e**2) - (13 / 4 * xa**2 + xb**2 / 4 + 13 / 3 * xa * xb)
    return pk["omega_dot_deg_yr"] * (1 + f_o * beta2) + omdot_lt_deg_yr


def _band(ax, M1, M2, Z, val, sig, color, lw=1.1, alpha=0.25):
    ax.contour(M1, M2, Z, levels=[val], colors=[color], linewidths=lw, linestyles="solid")
    if sig is not None and sig > 0:
        lo, hi = sorted([val - sig, val + sig])
        ax.contourf(M1, M2, Z, levels=[lo, hi], colors=[color], alpha=alpha)


def _y_on_curve(fun, x, ylo, yhi, target):
    """y in [ylo, yhi] where fun(x, y) = target (bisection), for placing labels on curves."""
    a, b = ylo, yhi
    fa = fun(x, a) - target
    for _ in range(80):
        mid = 0.5 * (a + b)
        fm = fun(x, mid) - target
        if np.sign(fm) == np.sign(fa):
            a, fa = mid, fm
        else:
            b = mid
    return 0.5 * (a + b)


def fig_binary_pulsars():
    pv = L.pulsar_value
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 6.2))
    # --- Hulse-Taylor (Weisberg & Huang 2016) ---
    ax = axes[0]
    pb, e = pv("b1913_wh16_pb")[0], pv("b1913_wh16_e")[0]
    x = pv("b1913_wh16_x")[0]
    M1, M2 = np.meshgrid(np.linspace(1.0, 1.9, 900), np.linspace(1.0, 1.8, 800))
    pk = G.pk_parameters(M1, M2, pb, e, x)
    om, s_om = pv("b1913_wh16_omdot")
    ga, s_ga = pv("b1913_wh16_gamma")
    pdi, s_pdi = pv("b1913_wh16_pbdot_intr")
    _band(ax, M1, M2, pk["omega_dot_deg_yr"], om, s_om, S.ORANGE)
    _band(ax, M1, M2, pk["gamma_s"], ga, s_ga, S.BLUE)
    _band(ax, M1, M2, pk["pbdot"], pdi, s_pdi, S.AQUA, alpha=0.35)
    _band(ax, M1, M2, pk["s"], 0.68, 0.08, S.VIOLET, alpha=0.10)
    mp, mc = pv("b1913_wh16_m1")[0], pv("b1913_wh16_m2")[0]
    ax.plot(mp, mc, marker="o", ms=5, color=S.INK, mfc="none", mew=1.0)
    ax.set_xlim(1.0, 1.9)
    ax.set_ylim(1.0, 1.8)
    ax.set_aspect("equal")
    ax.set_xlabel("pulsar mass (M☉)")
    ax.set_ylabel("companion mass (M☉)")
    ax.set_title("Hulse–Taylor pulsar B1913+16 (Weisberg & Huang 2016)", loc="left")
    f_om = lambda a, b: G.pk_parameters(a, b, pb, e)["omega_dot_deg_yr"]     # noqa: E731
    f_ga = lambda a, b: G.pk_parameters(a, b, pb, e)["gamma_s"]              # noqa: E731
    f_pd = lambda a, b: G.pk_parameters(a, b, pb, e)["pbdot"]                # noqa: E731
    f_s = lambda a, b: G.pk_parameters(a, b, pb, e, x)["s"]                  # noqa: E731
    for txt, col, xl, fun, target, xyt in [("periastron advance ω̇\n= 4.226585°/yr", S.ORANGE, 1.72, f_om, om, (1.55, 1.04)),
                                           ("Einstein delay γ = 4.307 ms", S.BLUE, 1.62, f_ga, ga, (1.52, 1.56)),
                                           ("orbital decay Ṗb = −2.398×10⁻¹²\n(gravitational waves)", S.AQUA, 1.13, f_pd, pdi,
                                            (1.15, 1.73)),
                                           ("Shapiro shape s = 0.68\n(band ±0.08)", S.VIOLET, 1.06, f_s, 0.68, (1.02, 1.42))]:
        yl = _y_on_curve(fun, xl, 1.0, 1.8, target)
        ax.annotate(txt, (xl, yl), xytext=xyt, textcoords="data", fontsize=7.4,
                    color=DARK[col], arrowprops={"arrowstyle": "-", "color": S.MUTED, "lw": 0.6})
    ax.annotate(f"all cross at {mp:.3f} + {mc:.3f} M☉", (mp, mc), xytext=(1.50, 1.20), fontsize=7.8, color=S.INK,
                arrowprops={"arrowstyle": "-", "color": S.MUTED, "lw": 0.7})
    # --- double pulsar (Kramer et al. 2021) ---
    ax = axes[1]
    pb, e, x = pv("dp_k21_pb")[0], pv("dp_k21_eT")[0], pv("dp_k21_x")[0]
    M1, M2 = np.meshgrid(np.linspace(1.20, 1.46, 1300), np.linspace(1.12, 1.38, 1300))
    pk = G.pk_parameters(M1, M2, pb, e, x)
    lt = pv("dp_k21_omdot_LT")[0]
    om, s_om = pv("dp_k21_omdot")
    ga, s_ga = pv("dp_k21_gammaE")
    pd_gw, s_pd = pv("dp_k21_pbdot_gw")
    pd_35 = pv("dp_k21_pbdot_gr35")[0]
    r_obs, s_r = pv("dp_k21_r")
    s_obs = pv("dp_k21_sini")[0]
    R, s_R = pv("dp_k21_R")
    _band(ax, M1, M2, _omdot_2pn(M1, M2, pb, e, lt), om, s_om, S.ORANGE)
    _band(ax, M1, M2, pk["gamma_s"] * 1e3, ga, s_ga, S.BLUE)
    _band(ax, M1, M2, pk["pbdot"] + pd_35, pd_gw, s_pd, S.AQUA, alpha=0.35)
    _band(ax, M1, M2, pk["r_s"] * 1e6, r_obs, s_r, S.MUTED, alpha=0.18)
    _band(ax, M1, M2, pk["s"], s_obs, 1e-5, S.VIOLET)
    _band(ax, M1, M2, M1 / M2, R, s_R, S.INK_SECONDARY, lw=0.9, alpha=0.12)
    ma, mb = pv("dp_k21_mA")[0], pv("dp_k21_mB")[0]
    ax.set_xlim(1.20, 1.46)
    ax.set_ylim(1.12, 1.38)
    ax.set_aspect("equal")
    ax.set_xlabel("mass of pulsar A (M☉)")
    ax.set_ylabel("mass of pulsar B (M☉)")
    ax.set_title("Double pulsar J0737−3039A/B (Kramer et al. 2021)", loc="left")
    f_om = lambda a, b: _omdot_2pn(a, b, pb, e, lt)                         # noqa: E731
    f_ga = lambda a, b: G.pk_parameters(a, b, pb, e)["gamma_s"] * 1e3       # noqa: E731
    f_pd = lambda a, b: G.pk_parameters(a, b, pb, e)["pbdot"] + pd_35       # noqa: E731
    f_s = lambda a, b: G.pk_parameters(a, b, pb, e, x)["s"]                 # noqa: E731
    for txt, col, xl, fun, target, off in [
            ("ω̇ (with 2PN + frame dragging)", S.ORANGE, 1.235, f_om, om, (1.248, 1.356)),
            ("Einstein delay γ_E", S.BLUE, 1.225, f_ga, ga, (1.205, 1.233)),
            ("orbital decay Ṗb", S.AQUA, 1.40, f_pd, pd_gw, (1.405, 1.215)),
            ("Shapiro shape s", S.VIOLET, 1.245, f_s, s_obs, (1.255, 1.168))]:
        yl = _y_on_curve(fun, xl, 1.12, 1.38, target)
        ax.annotate(txt, (xl, yl), xytext=off, textcoords="data", fontsize=7.4, color=DARK[col],
                    arrowprops={"arrowstyle": "-", "color": S.MUTED, "lw": 0.6})
    ax.text(1.205, 1.2535, "Shapiro range r", fontsize=7.4, color=S.INK_SECONDARY, va="bottom")
    ax.text(1.245, 1.125, "mass ratio R = x_B/x_A (theory-independent)", fontsize=7.4, color=S.INK_SECONDARY)
    ins = ax.inset_axes([0.60, 0.60, 0.37, 0.36])
    z1 = np.linspace(ma - 0.0012, ma + 0.0012, 700)
    Z1, Z2 = np.meshgrid(z1, np.linspace(mb - 0.0012, mb + 0.0012, 700))
    pz = G.pk_parameters(Z1, Z2, pb, e, x)
    _band(ins, Z1, Z2, _omdot_2pn(Z1, Z2, pb, e, lt), om, s_om, S.ORANGE, lw=0.9)
    _band(ins, Z1, Z2, pz["gamma_s"] * 1e3, ga, s_ga, S.BLUE, lw=0.9)
    _band(ins, Z1, Z2, pz["pbdot"] + pd_35, pd_gw, s_pd, S.AQUA, lw=0.9, alpha=0.35)
    _band(ins, Z1, Z2, pz["s"], s_obs, 1e-5, S.VIOLET, lw=0.9)
    ins.plot(ma, mb, marker="o", ms=4, color=S.INK, mfc="none")
    ins.set_xlim(z1[0], z1[-1])
    ins.set_ylim(mb - 0.0012, mb + 0.0012)
    ins.set_xticks([ma - 0.001, ma, ma + 0.001])
    ins.set_yticks([mb - 0.001, mb, mb + 0.001])
    ins.tick_params(labelsize=6)
    ins.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.3f"))
    ins.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.3f"))
    ins.set_title(f"zoomed ×100: {ma:.4f} + {mb:.4f} M☉", fontsize=7, fontweight="normal")
    S.caption(fig, "Each measured relativistic effect, read through general relativity's formula, allows a line of "
              "(pulsar, companion) masses; bands are ±1σ (most are thinner than the line). GR passes if all the lines cross "
              "at one point. Lines computed here from the published timing parameters (literature/checks/pulsars.csv); "
              "the Hulse–Taylor decay is after the Galactic correction, the double pulsar's after its kinematic corrections.",
              y=-0.01)
    return S.save(fig, OUT / "06_binary_pulsars")


# ---------------------------------------------------------------------------
# 7. The triple pulsar: what a strong-equivalence-principle violation would do
# ---------------------------------------------------------------------------

def fig_triple_pulsar(sig: dict, wn: dict):
    import triple_pulsar as T
    sim = T.simulate()
    t = sim["t"]
    odd = 0.5 * (sim["plus"] - sim["minus"])
    A = T.nuisance_basis(t)
    x, *_ = np.linalg.lstsq(A, odd, rcond=None)
    differential = T.DELTA_SIM * (1 + T.R14["m_p"] / T.R14["m_i"])
    w = 2 * 2 * np.pi / T.R14["P_I"] - 2 * np.pi / T.R14["P_O"]
    sig_curve = (x[-2] * np.cos(w * t) + x[-1] * np.sin(w * t)) / differential   # seconds per unit Delta
    resid = (odd - A @ x) / differential
    mjd, err = T.load_toas()
    fig = plt.figure(figsize=(15.2, 4.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.95, 1.45, 0.95], wspace=0.62)
    # (a) schematic
    ax = fig.add_subplot(gs[0])
    ax.add_patch(Ellipse((0.0, 0.25), 1.3, 1.0, fill=False, color=S.MUTED, lw=1.0))
    ax.plot(-0.65, 0.25, marker="o", ms=12, color=S.INK_SECONDARY, mec=S.SURFACE)
    ax.plot(0.65, 0.25, marker="o", ms=7, color=S.BLUE, mec=S.SURFACE)
    ax.text(-0.65, 0.05, "neutron star", ha="center", va="top", fontsize=7.4, color=S.INK_SECONDARY)
    ax.text(0.65, 0.07, "white\ndwarf", ha="center", va="top", fontsize=7.4, color=S.BLUE_550)
    ax.plot(2.35, 0.25, marker="o", ms=9, color=S.VIOLET, mec=S.SURFACE)
    ax.text(2.35, 0.05, "outer\nwhite dwarf", ha="center", va="top", fontsize=7.4, color=S.VIOLET)
    for x0, lab, ln in ((-0.65, "g (1 + Δ)", 0.75), (0.65, "g", 0.55)):
        ax.add_patch(FancyArrowPatch((x0, 0.55), (x0 + ln, 0.55), arrowstyle="-|>", mutation_scale=11,
                                     color=S.INK, lw=1.1))
        ax.text(x0 + ln / 2, 0.62, lab, ha="center", va="bottom", fontsize=7.8, color=S.INK)
    ax.text(-1.1, -0.42, "Both fall toward the outer star (g ≈ 1.7 mm/s²,\n1.2 AU away, orbit 327 days). A neutron "
            "star's\nbinding energy is ~10% of its mass; if that\nmade it fall faster, the 1.6-day inner orbit "
            "would\nbe stretched toward the outer star, and the\nstretch would turn with it.",
            fontsize=7.4, color=S.INK_SECONDARY, va="top")
    ax.set_xlim(-1.15, 2.75)
    ax.set_ylim(-1.45, 1.05)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("A neutron star and a white dwarf\nfalling together (not to scale)", loc="left")
    # (b) the signature at nanosecond scale
    ax = fig.add_subplot(gs[1])
    t0 = t[0] + 200.0
    sel = (t >= t0) & (t <= t0 + 5.0)
    dlt = 1.46e-6
    ax.axhspan(-wn["amp_sigma_ns"], wn["amp_sigma_ns"], color=S.WASH, lw=0,
               label=f"all {len(mjd):,} arrival times together: ±{wn['amp_sigma_ns']:.0f} ns (1σ, if the noise were white)")
    ax.plot(t[sel] - t0, sig_curve[sel] * dlt * 1e9, color=S.ORANGE, lw=1.4,
            label=f"signal if Δ = 1.46×10⁻⁶ (today's best limit): ±{sig['amplitude_ns_per_1e-6'] * 1.46:.0f} ns")
    ax.axhline(0, color=S.AXIS, lw=0.8)
    med = np.median(err) * 1e6
    ax.text(0.05, 128, f"one arrival time is uncertain by ±{med:.1f} µs, {med * 1e3 / 140:.0f}× the height of this plot",
            fontsize=7.4, color=S.INK_SECONDARY, va="top")
    ax.set_xlim(0, 5)
    ax.set_ylim(-140, 140)
    ax.set_xlabel("days")
    ax.set_ylabel("change in pulse arrival times (ns)")
    ax.legend(loc="lower left", fontsize=7.2, framealpha=0.95, frameon=True, facecolor=S.SURFACE, edgecolor=S.SURFACE)
    ax.set_title(f"The signal: {sig['amplitude_ns_per_1e-6']:.0f} ns per Δ = 10⁻⁶,\nrepeating every "
                 f"{2 * np.pi / w:.3f} d (2 f_inner − f_outer)", loc="left")
    # (c) published limits vs our white-noise expectation
    ax = fig.add_subplot(gs[2])
    rows = [("Archibald et al. 2018", 2.6e-6, S.MUTED), ("Voisin et al. 2020", 2.05e-6, S.MUTED),
            ("Voisin et al. 2025, red noise", 2.29e-6, S.MUTED), ("Voisin et al. 2025, planet", 1.46e-6, S.MUTED),
            ("white noise only, their\n12,474 times (this work)", 1.96 * wn["delta_sigma_white"], S.BLUE)]
    for i, (lab, v, col) in enumerate(rows):
        ax.barh(i, v * 1e6, color=col, height=0.55)
        ax.text(v * 1e6 + 0.05, i, f"{v * 1e6:.2f}", va="center", fontsize=7.6, color=S.INK_SECONDARY)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=7.4)
    ax.invert_yaxis()
    ax.set_xlim(0, 3.2)
    ax.set_xlabel("95% limit on |Δ| (×10⁻⁶)")
    ax.set_title("How closely the neutron star\nfalls like a white dwarf", loc="left")
    S.caption(fig, "Labelled simulation with real initial conditions: Ransom et al. (2014)'s three-body solution integrated with "
              "REBOUND, with and without Δ (only the difference in how the two inner bodies fall is kept; a common change is "
              "absorbed by the masses). After fitting out everything a timing model absorbs, the signal is "
              f"{sig['amplitude_ns_per_1e-6']:.1f} ns per Δ = 10⁻⁶ (Archibald et al. 2018: 30 ns). Arrival times: Voisin et al. "
              "(2025), Zenodo 10.5281/zenodo.13899771; limits from literature/checks/pulsars.csv.", y=-0.03)
    return S.save(fig, OUT / "07_triple_pulsar")


# ---------------------------------------------------------------------------
# 8. Lunar laser ranging with APOLLO
# ---------------------------------------------------------------------------

def fig_lunar_ranging(tables: dict, res: dict, c_nord: float):
    m = tables["EPM2021"]
    fit = res["EPM2021"]
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.6))
    ax = axes[0, 0]
    refl_style = {"Apollo 15": (S.BLUE, "o"), "Apollo 11": (S.ORANGE, "s"), "Apollo 14": (S.AQUA, "D"),
                  "Lunokhod 1": (S.VIOLET, "^"), "Lunokhod 2": (S.INK_SECONDARY, "v")}
    for r, (col, mk) in refl_style.items():
        d = m[m["reflector"] == r]
        ax.plot(d["years"], d["rt_s"], ls="", marker=mk, ms=2.2, color=col, label=r, alpha=0.8)
    ax.set_xlabel("year")
    ax.set_ylabel("round trip, telescope → Moon → telescope (s)")
    ax.set_title(f"{len(m):,} laser round trips to the Moon (APOLLO): 2.33-2.68 s", loc="left")
    ax.legend(loc="lower left", fontsize=7, markerscale=2.5, ncols=3, framealpha=0.9, frameon=True,
              facecolor=S.SURFACE, edgecolor=S.SURFACE, title="reflector", title_fontsize=7)
    ax.set_ylim(2.25, 2.70)
    ax = axes[0, 1]
    mm = m["sigma_s"].to_numpy() * C.C_LIGHT / 2 * 1e3
    ax.hist(mm, bins=np.logspace(np.log10(0.5), np.log10(40), 45), color=S.BLUE_300, edgecolor=S.SURFACE, lw=0.8)
    ax.set_xscale("log")
    med = np.median(mm)
    ax.axvline(med, color=S.INK, lw=1.0)
    top = ax.get_ylim()[1]
    ax.text(med * 1.08, top * 0.95, f"median {med:.1f} mm\n= {np.median(m['sigma_s']) * 1e12:.0f} ps of a\n2.5-s trip: "
            "1 part in 10¹¹", fontsize=7.6, va="top", color=S.INK)
    rms = fit["rms_m"] * 1e3
    ax.axvline(rms, color=S.ORANGE, lw=1.0)
    ax.text(rms * 0.93, top * 0.95, f"our model\nfits to {rms / 10:.1f} cm", fontsize=7.6, va="top", ha="right",
            color=DARK[S.ORANGE])
    ax.set_xlabel("quoted uncertainty of one normal point (mm, one way)")
    ax.set_ylabel("normal points")
    ax.set_title("Measured to millimetres", loc="left")
    ax = axes[1, 0]
    u = m[m["used"]]
    ax.plot(u["years"], u["post_resid_m"] * 100, ls="", marker="o", ms=1.8, color=S.BLUE, alpha=0.5)
    night = (u["year"] * 10000 + u["month"] * 100 + u["day"]).to_numpy()
    g = pd.DataFrame({"t": u["years"].to_numpy(), "r": u["post_resid_m"].to_numpy() * 100, "n": night}).groupby("n").mean()
    ax.plot(g["t"], g["r"], ls="", marker="o", ms=3.0, color=S.INK, mfc=S.SURFACE, mew=0.7, label="nightly mean")
    ax.axhline(0, **S.GR_LINE)
    ax.set_ylim(-12, 12)
    ax.set_xlabel("year")
    ax.set_ylabel("measured − predicted range (cm)")
    ax.set_title(f"Measured minus predicted: {fit['rms_m'] * 100:.1f} cm rms in 385,000 km", loc="left")
    ax.legend(loc="lower left", fontsize=7.4)
    ax = axes[1, 1]
    p = fit["params"]
    eta = p["eta_nordtvedt"]["value"]
    D = u["D_rad"].to_numpy()
    r_plus = (u["post_resid_m"].to_numpy() + c_nord * eta * np.cos(D)) * 1e3
    bins = np.linspace(0, 2 * np.pi, 19)
    idx = np.digitize(D, bins) - 1
    nights = pd.Series(night)
    cen, mean, err = [], [], []
    for k in range(len(bins) - 1):
        s = idx == k
        if s.sum() < 10:
            continue
        nm = pd.Series(r_plus[s]).groupby(nights[s].to_numpy()).mean()
        cen.append(0.5 * (bins[k] + bins[k + 1]))
        mean.append(nm.mean())
        err.append(nm.std(ddof=1) / np.sqrt(len(nm)))
    ax.errorbar(np.degrees(cen), mean, yerr=err, fmt="o", ms=4, color=S.INK, mfc=S.SURFACE, elinewidth=0.8, capsize=0,
                label="APOLLO residuals by lunar phase (night means ± error)")
    dd = np.linspace(0, 2 * np.pi, 300)
    ax.plot(np.degrees(dd), c_nord * 1e-3 * np.cos(dd) * 1e3, color=S.ORANGE, lw=1.3,
            label=f"if η were 0.001: {c_nord:.0f} mm × cos D (our simulation)")
    ax.plot(np.degrees(dd), c_nord * eta * np.cos(dd) * 1e3, color=S.BLUE, lw=1.3,
            label=f"best fit: η = ({eta * 1e4:.1f} ± {p['eta_nordtvedt']['boot'] * 1e4:.1f})×10⁻⁴")
    ax.axhline(0, **S.GR_LINE)
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_xticklabels(["new Moon", "first quarter", "full Moon", "last quarter", "new Moon"], fontsize=7.4)
    ax.set_ylabel("range residual (mm)")
    ax.set_ylim(-24, 22)
    ax.set_title("The Nordtvedt test: does the Moon's orbit lean toward the Sun?", loc="left")
    ax.legend(loc="lower left", fontsize=7.2)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    S.caption(fig, "Our own range model (JPL's 2020 recipe: Earth rotation, tides, atmosphere, the Sun's 7.5-m Shapiro delay, "
              "relativistic time and length transformations) with the EPM2021 lunar ephemeris. Fitted: a range bias per "
              "APOLLO hardware group, the telescope's position and drift, and the reflector coordinates. The ephemeris assumes "
              "GR and was fitted to lunar ranges, so the η fit is a consistency test, not an independent measurement.", y=0.03)
    return S.save(fig, OUT / "08_lunar_ranging")


# ---------------------------------------------------------------------------
# 9. Where is the Moon? Four ephemerides and the APOLLO data
# ---------------------------------------------------------------------------

def fig_moon_along_track(curves: dict, res: dict):
    fig, ax = plt.subplots(figsize=(10.8, 5.8))
    styles = {"DE421 (JPL 2008)": S.BLUE, "DE430 (JPL 2013)": S.AQUA,
              "DE441 (JPL 2020, via Horizons)": S.VIOLET, "EPM2021 (IAA RAS)": S.ORANGE}
    for name, (t, tau) in curves.items():
        ax.plot(t, tau, color=styles[name], lw=1.4, label=name)
    yrs = res["DE440"]["along_track_by_year"]
    ty = np.array(sorted(yrs)) + 0.5
    off = np.array([yrs[k]["offset_ms"] for k in sorted(yrs)])
    er = np.array([yrs[k]["err_ms"] for k in sorted(yrs)])
    ax.errorbar(ty, off, yerr=er, fmt="o", ms=5, color=S.INK, mfc=S.SURFACE, elinewidth=0.9, capsize=0,
                label="what the APOLLO ranges prefer (our model, year by year)", zorder=5)
    ax.axhline(0, **S.GR_LINE)
    ax.text(2000.3, 0.10, "DE440 (JPL 2020): the reference", fontsize=7.8, color=S.INK)
    ax.set_xlim(2000, 2022.5)
    ax.set_ylim(-3.2, 2.3)
    ax.set_xlabel("year")
    ax.set_ylabel("Moon's lead along its orbit, relative to DE440 (ms)")
    ax.set_title("Where exactly is the Moon? Ephemerides disagree by a metre, and APOLLO can tell", loc="left")
    ax.legend(loc="lower left", fontsize=7.8)
    ax.text(0.99, 0.60, "1 ms of orbital motion ≈ 1.0 m along the orbit;\nin range it shows as 5.6 cm × sin(mean anomaly) "
            "per ms", transform=ax.transAxes, ha="right", va="center", fontsize=7.6, color=S.INK_SECONDARY)
    S.caption(fig, "Each line: the time shift that best maps DE440's geocentric Moon onto the other ephemeris (15-day steps). "
              "DE441 differs from DE440 by a pure shift (1 cm left over); the others also differ by up to 0.6-1.1 m in other "
              "directions, mostly a constant radial offset that fitted reflector and station coordinates absorb. Points: the "
              "shift the APOLLO ranges prefer when our model uses DE440 (fitted jointly with the station, biases and "
              "reflectors; bootstrap errors over nights).", y=-0.01)
    return S.save(fig, OUT / "09_moon_along_track")


# ---------------------------------------------------------------------------
# 10. The Nordtvedt signal, simulated
# ---------------------------------------------------------------------------

def fig_nordtvedt(amp: dict):
    import nordtvedt as N
    t, y, nonlin = N.template()
    D, l, F = N.mean_elongation(t)
    fig = plt.figure(figsize=(13.0, 4.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.65, 1.0], wspace=0.55)
    ax = fig.add_subplot(gs[0])
    sel = (t > 800) & (t < 1000)
    ax.plot(t[sel] - 800, y[sel], color=S.ORANGE, lw=1.4, label="Earth–Moon distance change per unit η (simulation)")
    ax.plot(t[sel] - 800, amp["cosD"] * np.cos(D[sel]), color=S.INK, lw=0.9, label=f"{amp['cosD']:.2f} m × cos D")
    ax.axhline(0, color=S.AXIS, lw=0.8)
    ax.set_ylim(-19, 19)
    ax.set_xlabel("days (from March 2006)")
    ax.set_ylabel("δr per unit η (m)")
    ax.set_title("If the Earth fell toward the Sun a little faster than the Moon", loc="left")
    ax.legend(loc="lower center", fontsize=7.8, ncols=2)
    ax = fig.add_subplot(gs[1])
    lit = [("Nordtvedt 1968 ('about 12 m')", 12.0), ("pre-1995 standard estimate", 9.3),
           ("Damour & Vokrouhlický 1996", 13.10), ("Will 2014 (Nordtvedt 1995)", 13.1),
           ("Williams et al. 2009 (numerical)", 13.3), ("Williams et al. 2012 (numerical)", 2.992e13 * 4.45e-10 / 1e3),
           ("this work (REBOUND, DE440 start)", amp["cosD"])]
    for i, (lab, v) in enumerate(lit):
        ours = lab.startswith("this work")
        ax.plot(v, i, marker="*" if ours else "o", ms=11 if ours else 6, color=S.INK if ours else S.BLUE, mec=S.SURFACE)
        ax.text(v + 0.15, i, f"{v:.2f}", va="center", fontsize=7.4, color=S.INK_SECONDARY)
    ax.set_yticks(range(len(lit)))
    ax.set_yticklabels([x[0] for x in lit], fontsize=7.4)
    ax.invert_yaxis()
    ax.set_xlim(8.5, 14.5)
    ax.set_xlabel("signal amplitude per unit η (m)")
    ax.set_title("How big is the signal?", loc="left")
    S.caption(fig, "Labelled simulation with real initial conditions: the Sun, planets, Earth and Moon from DE440 integrated with "
              "REBOUND (IAS15), with an extra force making the Earth's and the Moon's fall toward the Sun differ by "
              "η × (4.64 − 0.19)×10⁻¹⁰ (their gravitational self-energies), switched on over two years. Literature values from "
              "literature/checks/beta_lunar.csv (coefficients converted with the same 4.45×10⁻¹⁰).", y=-0.03)
    return S.save(fig, OUT / "10_nordtvedt_signal")
