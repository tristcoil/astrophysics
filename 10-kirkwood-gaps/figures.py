"""Every figure. Each function takes the analysis output and writes PNG + SVG.

Numbering follows the story: the gaps (01), measuring them (02), are they
real (03), the belt's architecture (04), families (05), where everything is
today (06), how resonances protect the Hildas and Trojans (07), and the
labelled N-body simulation (08).
"""

from __future__ import annotations

import json
from fractions import Fraction

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

import orbit_physics as op
import plotstyle as ps
from config import (
    EARTH_APHELION_AU,
    FIGURE_DIR,
    JUPITER_PERIHELION_AU,
    MARS_APHELION_AU,
    PROVENANCE_PATH,
)

ps.use_style()


def _retrieved() -> str:
    prov = json.loads(PROVENANCE_PATH.read_text())
    return prov["MPCORB.DAT.gz"]["retrieved_utc"][:10]


def _mpc_credit(n_shown: int | None = None) -> str:
    shown = f"; {n_shown:,} orbits shown" if n_shown else ""
    return (f"Data: Minor Planet Center, MPCORB.DAT retrieved {_retrieved()}{shown}. "
            "Osculating elements at epoch 2026 June 9 TT.")


def _hist(a_sorted: np.ndarray, lo: float, hi: float, width: float):
    edges = np.arange(lo, hi + width / 2, width)
    counts = np.diff(np.searchsorted(a_sorted, edges))
    return edges, counts


def _step_fill(ax, edges, counts, color=ps.C_DATA, alpha=0.16, lw=1.0, label=None,
               zorder=2):
    """A histogram drawn as a thin outline over a light wash of the same hue."""
    ax.stairs(counts, edges, fill=True, color=color, alpha=alpha, zorder=zorder,
              linewidth=0)
    ax.stairs(counts, edges, color=color, linewidth=lw, zorder=zorder + 0.1, label=label)


def _ratio_axis(ax, ticks: list[Fraction], *, label: str = "orbits per Jupiter orbit"):
    """A top axis in period ratio -- the same quantity as the bottom axis,
    converted by Kepler's third law, so the resonances land on round numbers."""
    def fwd(a):
        return op.period_ratio(np.clip(np.asarray(a, dtype=float), 0.05, None))

    def inv(r):
        return op.a_from_period_ratio(np.clip(np.asarray(r, dtype=float), 1e-3, None))

    sec = ax.secondary_xaxis("top", functions=(fwd, inv))
    sec.set_xticks([float(t) for t in ticks])
    sec.set_xticklabels([f"{t.numerator}" if t.denominator == 1
                         else f"{t.numerator}/{t.denominator}" for t in ticks])
    sec.set_xlabel(label, color=ps.INK_SECONDARY)
    sec.tick_params(axis="x", colors=ps.MUTED, labelsize=8, length=3)
    sec.spines["top"].set_visible(True)
    sec.spines["top"].set_color(ps.AXIS)
    return sec


def _res_ticks(ax, resonances, *, y0=0.965, y1=1.0, color=ps.C_PREDICTED, lw=1.2):
    """Short ticks hanging from the top of the plot at resonance locations."""
    tr = ax.get_xaxis_transform()
    for p, q in resonances:
        a = op.resonance_a(p, q)
        ax.plot([a, a], [y0, y1], transform=tr, color=color, lw=lw, solid_capstyle="butt",
                zorder=5, clip_on=False)


# ---------------------------------------------------------------------------
# 01 -- the Kirkwood gaps
# ---------------------------------------------------------------------------

def fig_kirkwood_gaps(res: dict, arr: dict) -> list[str]:
    a = arr["a_osc_sorted"]
    fig = plt.figure(figsize=(10.5, 8.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.45, 1.0], hspace=0.42)

    # (a) The main belt, linear.
    ax = fig.add_subplot(gs[0])
    edges, counts = _hist(a, 2.0, 3.4, 0.004)
    _step_fill(ax, edges, counts)
    ax.set_xlim(2.0, 3.4)
    ax.set_ylim(0, counts.max() * 1.5)
    ax.set_xlabel("semi-major axis (AU) — the average distance from the Sun")
    ax.set_ylabel("asteroids per 0.004 AU")
    ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter("{x:,.0f}"))
    ticks = [Fraction(4), Fraction(7, 2), Fraction(3), Fraction(8, 3), Fraction(5, 2),
             Fraction(7, 3), Fraction(9, 4), Fraction(2)]
    _ratio_axis(ax, ticks)
    _res_ticks(ax, [(t.numerator, t.denominator) for t in ticks])

    named = {r["label"]: r["osculating"] for r in res["gaps"]["named"]}
    ymax = ax.get_ylim()[1]
    labels = {
        "3:1": "3:1", "5:2": "5:2", "7:3": "7:3", "2:1": "2:1",
        "4:1": "4:1", "8:3": "8:3", "9:4": "9:4",
    }
    for lab, text in labels.items():
        p, q = map(int, lab.split(":"))
        x = op.resonance_a(p, q)
        k = np.searchsorted(edges, x) - 1
        local = counts[max(k - 3, 0):k + 4].max() if 0 <= k < len(counts) else 0
        strong = lab in ("3:1", "5:2", "7:3", "2:1")
        y = min(local + ymax * 0.04, ymax * 0.7) if not strong else ymax * 0.7
        ax.text(x, y, text, ha="center", va="bottom", fontsize=9 if strong else 7.5,
                color=ps.INK if strong else ps.INK_SECONDARY,
                fontweight="bold" if strong else "normal")
    ax.set_title("(a) The Kirkwood gaps: the main belt is empty wherever an orbit "
                 "would keep step with Jupiter's")
    ax.text(2.012, counts.max() * 1.44,
            "Each gap sits where an asteroid's orbit would fit a simple\n"
            "ratio with Jupiter's. In the 3:1 gap it would circle the Sun\n"
            "3 times per Jupiter orbit, so the same alignments recur every\n"
            "11.9 years and Jupiter's tugs add up instead of averaging out.",
            ha="left", va="top", fontsize=8, color=ps.INK_SECONDARY, linespacing=1.3)

    # (b) The whole span out to Jupiter, log scale: gaps and resonant groups.
    ax = fig.add_subplot(gs[1])
    edges, counts = _hist(a, 1.6, 5.6, 0.01)
    _step_fill(ax, edges, np.maximum(counts, 0.8), alpha=0.16)
    ax.set_yscale("log")
    ax.set_xlim(1.6, 5.6)
    ax.set_ylim(0.8, counts.max() * 12)
    ax.set_xlabel("semi-major axis (AU)")
    ax.set_ylabel("asteroids per 0.01 AU")
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
        lambda v, _: f"{v:,.0f}" if v >= 1 else ""))
    ticks_b = [Fraction(4), Fraction(3), Fraction(5, 2), Fraction(2), Fraction(3, 2),
               Fraction(4, 3), Fraction(1)]
    _ratio_axis(ax, ticks_b)
    _res_ticks(ax, [(t.numerator, t.denominator) for t in ticks_b])
    groups = [
        ("Hungarias", 1.94, None),
        ("main belt", 2.75, None),
        ("Cybeles", 3.47, None),
        ("Hildas\nlocked in the 3:2", op.resonance_a(3, 2), "3:2"),
        ("the tiny\nThule group (4:3)", op.resonance_a(4, 3), "4:3"),
        ("Trojans\nshare Jupiter's orbit", op.resonance_a(1, 1), "1:1"),
    ]
    for text, x, _ in groups:
        k = np.searchsorted(edges, x) - 1
        local = counts[max(k - 2, 0):k + 3].max()
        ax.text(x, max(local, 1) * 2.2, text, ha="center", va="bottom", fontsize=8,
                color=ps.INK_SECONDARY, linespacing=1.1)
    ax.set_title("(b) Out to Jupiter: some resonances empty the belt, others hold "
                 "whole populations in place")
    ps.caption(fig, _mpc_credit(len(a)) + " Top axes: how many orbits an asteroid "
               "at that distance completes per orbit of Jupiter (Kepler's third law, "
               f"Jupiter's period {op.JUPITER_PERIOD_DAYS:.2f} d); violet ticks mark "
               "the resonances.")
    return ps.save(fig, FIGURE_DIR / "01_kirkwood_gaps")


# ---------------------------------------------------------------------------
# 02 -- measuring the gaps
# ---------------------------------------------------------------------------

ZOOM_GAPS = ["3:1", "8:3", "5:2", "7:3", "9:4", "11:5"]


def fig_gap_measurements(res: dict, arr: dict) -> list[str]:
    a = arr["a_osc_sorted"]
    det = {g["label"]: g for g in res["gaps"]["detections"] if g["kind"] == "gap"
           and g["credible"]}
    named = {r["label"]: r for r in res["gaps"]["named"]}

    fig = plt.figure(figsize=(11, 8.2))
    gs = fig.add_gridspec(2, 6, height_ratios=[1.0, 1.05], hspace=0.55, wspace=0.55)
    for i, lab in enumerate(ZOOM_GAPS):
        ax = fig.add_subplot(gs[0, i])
        g = det[lab]
        p, q = g["p"], g["q"]
        a_res = op.resonance_a(p, q)
        half = max(1.6 * g["fwhm"], 0.012)
        edges, counts = _hist(a, a_res - half, a_res + half, 0.001)
        _step_fill(ax, edges, counts, lw=0.9)
        ax.axvline(a_res, color=ps.C_PREDICTED, lw=1.4, zorder=4)
        ax.axvspan(g["center"] - g["center_err"], g["center"] + g["center_err"],
                   color=ps.C_MEASURED, alpha=0.25, lw=0, zorder=3)
        ax.axvline(g["center"], color=ps.C_MEASURED, lw=1.4, zorder=4)
        ax.set_xlim(a_res - half, a_res + half)
        ax.set_ylim(0, counts.max() * 1.12)
        ax.set_xticks([round(a_res - half / 2, 3), round(a_res + half / 2, 3)])
        ax.tick_params(axis="x", labelsize=7)
        ax.tick_params(axis="y", labelsize=7)
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
        ax.set_title(f"{lab}", fontsize=10, loc="center")
        d = 1000 * (g["center"] - a_res)
        e = 1000 * g["center_err"]
        ax.text(0.5, -0.17, f"{d:+.1f} ± {e:.1f} mAU", transform=ax.transAxes,
                ha="center", va="top", fontsize=8, color=ps.INK)
        if i == 0:
            ax.set_ylabel("asteroids per 0.001 AU", fontsize=8)
    fig.text(0.07, 0.955, "(a) Each gap found blind, then compared with where Kepler's "
             "third law puts the resonance", fontsize=11, fontweight="bold", ha="left")
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], color=ps.C_MEASURED, lw=2, label="measured centre"),
                        Line2D([], [], color=ps.C_PREDICTED, lw=2,
                               label="predicted resonance")],
               loc="upper left", bbox_to_anchor=(0.065, 0.94), ncol=2, fontsize=8)
    fig.text(0.93, 0.925, "under each panel: measured − predicted, in mAU "
             "(1 mAU = 0.001 AU ≈ 150,000 km)", fontsize=8, color=ps.INK_SECONDARY,
             ha="right")

    # (b) Offsets for every gap, osculating and proper, against gap widths.
    ax = fig.add_subplot(gs[1, :])
    rows = ["4:1", "7:2", "3:1", "8:3", "5:2", "7:3", "9:4", "11:5", "2:1"]
    y = np.arange(len(rows))[::-1]
    for yi, lab in zip(y, rows):
        p, q = map(int, lab.split(":"))
        a_res = op.resonance_a(p, q)
        g_det = det.get(lab)
        osc = named[lab]["osculating"] if lab in named else g_det
        prop = named[lab].get("proper") if lab in named else None
        if osc is None:
            continue
        w = osc["fwhm"] / 2 * 1000
        edge = osc.get("walls", {}).get("is_edge", False)
        if not edge:
            ax.barh(yi, 2 * w, left=-w, height=0.56, color=ps.C_CONTEXT, alpha=0.6, lw=0,
                    zorder=1, label="each gap's half-width, for scale" if lab == "3:1"
                    else None)
        ax.errorbar(1000 * (osc["center"] - a_res), yi + 0.13, xerr=1000 * osc["center_err"],
                    fmt="o", color=ps.C_MEASURED, ms=6, mec=ps.SURFACE, mew=1.2,
                    elinewidth=1.3, capsize=0, zorder=4,
                    label="osculating elements (MPC)" if lab == "3:1" else None)
        if prop is not None and not edge:
            ax.errorbar(1000 * (prop["center"] - a_res), yi - 0.13,
                        xerr=1000 * prop["center_err"], fmt="s", color=ps.BLUE, ms=5.5,
                        mec=ps.SURFACE, mew=1.2, elinewidth=1.3, capsize=0, zorder=4,
                        label="proper elements (AstDyS)" if lab == "3:1" else None)
        note = ("  edge of the belt — no symmetric centre" if edge else
                "  weak; missed by the blind search" if g_det is None else "")
        ax.text(-62, yi, f"{lab}{note}", va="center", ha="left", fontsize=8.5,
                color=ps.INK if not (edge or g_det is None) else ps.MUTED)
    ax.axvline(0, color=ps.C_PREDICTED, lw=1.4, zorder=2, label="predicted resonance")
    ax.set_xlim(-62, 62)
    ax.set_yticks([])
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xlabel("measured centre − predicted resonance (mAU)")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4, fontsize=8)
    chance = res["gaps"]["chance_alignment"]
    if chance["hits"] == 0:
        odds = f"under 1 in {chance['trials'] / 3 / 1e6:.0f} million"
    else:
        inv = 1 / chance["probability"]
        odds = (f"about 1 in {inv / 1e6:.1f} million" if inv > 1e6 else
                f"about 1 in {inv:,.0f}")
        odds += f" ({chance['hits']} hits in {chance['trials'] / 1e6:.0f} million trials)"
    ax.set_title("(b) How far each gap sits from its resonance, in two independent "
                 "element sets", pad=30)
    ps.caption(fig, _mpc_credit(len(a)) + " Proper elements: AstDyS synthetic catalogue "
               "(Knežević & Milani), June 2024. Gaps found by a blind search for dips with "
               "a wall on both sides; centre = mean of a model fit and the smoothed "
               "histogram's deepest point, uncertainty = their disagreement or the fit's "
               "sensitivity to window and background, whichever is larger. If gaps fell at "
               f"random, the chance that the four strongest would lie this close to "
               f"resonances with denominators ≤ 3 is {odds}.")
    return ps.save(fig, FIGURE_DIR / "02_gap_measurements")


# ---------------------------------------------------------------------------
# 03 -- are the gaps real?
# ---------------------------------------------------------------------------

def _ridge(ax, rows, *, lo=2.0, hi=3.4, gap_lines=True):
    """Small multiples stacked vertically: each row normalised to its own
    peak, so shapes compare regardless of sample size."""
    n = len(rows)
    for i, (label, a_sorted, width, rug) in enumerate(rows):
        base = n - 1 - i
        if rug:
            sel = a_sorted[(a_sorted > lo) & (a_sorted < hi)]
            ax.vlines(sel, base + 0.05, base + 0.55, color=ps.C_DATA, lw=0.9)
        else:
            edges, counts = _hist(a_sorted, lo, hi, width)
            h = counts / counts.max() * 0.9
            ax.stairs(h + base, edges, baseline=base, fill=True, color=ps.C_DATA,
                      alpha=0.16, lw=0)
            ax.stairs(h + base, edges, baseline=None, color=ps.C_DATA, lw=0.9)
        ax.axhline(base, color=ps.AXIS, lw=0.6, zorder=0)
        ax.text(lo - 0.03, base + 0.3, label, ha="right", va="center", fontsize=8,
                color=ps.INK_SECONDARY, linespacing=1.15)
    if gap_lines:
        for p, q in ((3, 1), (5, 2), (7, 3), (2, 1)):
            ax.axvline(op.resonance_a(p, q), color=ps.C_PREDICTED, lw=0.8, alpha=0.7,
                       zorder=0)
    ax.set_xlim(lo, hi)
    ax.set_ylim(-0.1, n)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.grid(False)


def fig_is_it_real(res: dict, arr: dict, d: dict) -> list[str]:
    q = arr["quality"]
    fig = plt.figure(figsize=(12, 8.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1], height_ratios=[1.15, 1],
                          hspace=0.42, wspace=0.42)

    # (a) By size.
    ax = fig.add_subplot(gs[0, 0])
    rows = []
    for row in res["reality"]["by_H"]:
        lo, hi = row["H_range"]
        sel = q & (d["H"] >= lo) & (d["H"] < hi)
        dmin, dmax = row["diameter_km_at_albedo_0.14"]
        if lo < 0:
            lab = f"H < {hi:g}\nwider than ~{dmin:.0f} km\n{row['n']:,} asteroids"
        elif hi >= 30:
            lab = f"H > {lo:g}\nunder ~{dmax:.1f} km\n{row['n']:,} asteroids"
        else:
            lab = f"H {lo:g}–{hi:g}\n~{dmin:.1f}–{dmax:.1f} km\n{row['n']:,} asteroids"
        rows.append((lab, np.sort(d["a"][sel]), 0.008 if row["n"] < 20000 else 0.004, False))
    _ridge(ax, rows)
    ax.set_xlabel("semi-major axis (AU)")
    ax.set_title("(a) Big or small, the gaps are in the same places")

    # (b) By era: the catalogue in numbering order.
    ax = fig.add_subplot(gs[0, 1])
    numbered = d["number"] > 0
    order = np.argsort(d["number"][numbered])
    a_num = d["a"][numbered][order]
    eras = res["reality"]["by_era"]
    rows = []
    for era in eras:
        n = era["first_n"]
        s = np.sort(a_num[:n])
        if n == 91:
            lab = "the 91 known\nby the end of 1866"
            rows.append((lab, s, None, True))
        else:
            lab = f"first {n:,}\nnumbered" if n < 800_000 else f"all {n:,}\nnumbered"
            width = {1000: 0.02, 10_000: 0.01}.get(n, 0.004)
            rows.append((lab, s, width, False))
    _ridge(ax, rows)
    ax.set_ylim(-0.1, len(rows) + 0.1)
    k66 = eras[0]
    ax.set_xlabel("semi-major axis (AU)")
    ax.set_title("(b) How many asteroids does it take to see the gaps?")
    ax.text(3.38, len(rows) - 0.2,
            f"none in the core of the 3:1 gap, where {k66['3:1_expected']:.1f} were expected",
            ha="right", va="center", fontsize=7.5, color=ps.INK_SECONDARY)

    # (c) Discovery bias.
    ax = fig.add_subplot(gs[1, :])
    good = q & np.isfinite(d["H"]) & (d["H"] < 30)
    hh, xe, ye = np.histogram2d(d["a"][good], d["H"][good], bins=[520, 170],
                                range=[[1.6, 5.5], [4, 21]])
    im = ax.pcolormesh(xe, ye, hh.T, cmap=ps.density_cmap(), norm=LogNorm(vmin=1),
                       rasterized=True, zorder=1)
    faint = np.array(res["reality"]["faint_limit_95pct"])
    ax.plot(faint[:, 0], faint[:, 1], color=ps.C_MEASURED, lw=2, zorder=3,
            label="95% of known asteroids are brighter than this\n(orbits that never "
                  "cross Mars')")
    ax.set_ylim(21, 4)
    ax.set_xlim(1.6, 5.5)
    ax.set_xlabel("semi-major axis (AU)")
    ax.set_ylabel("absolute magnitude H\n(bigger number = smaller asteroid)")
    ax.legend(loc="lower right", fontsize=8)
    cb = fig.colorbar(im, ax=ax, pad=0.01, fraction=0.03)
    cb.set_label("asteroids per cell", color=ps.INK_SECONDARY)
    cb.outline.set_visible(False)
    ax.set_title("(c) What the catalogue misses: the farther out, the bigger an "
                 "asteroid must be to have been found")
    for p, qq in ((3, 1), (5, 2), (7, 3), (2, 1), (3, 2)):
        ax.axvline(op.resonance_a(p, qq), color=ps.C_PREDICTED, lw=0.8, alpha=0.6, zorder=2)
    ps.caption(fig, _mpc_credit() + " (a) orbits with U ≤ 2, split by absolute magnitude; "
               "diameters assume an albedo of 0.14 and are uncertain by a factor of ~2. "
               "(b) numbered asteroids in numbering order; asteroids (1)–(91) were "
               "discovered by 1866 Nov 4, (92) in 1867 (JPL Small-Body Database). Violet "
               "lines: the 3:1, 5:2, 7:3 and 2:1 resonances.")
    return ps.save(fig, FIGURE_DIR / "03_is_it_real")


# ---------------------------------------------------------------------------
# 04 -- the architecture of the belt: a-e and a-i
# ---------------------------------------------------------------------------

def fig_architecture(res: dict, arr: dict, d: dict) -> list[str]:
    q = arr["quality"]
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 9.4), sharex=True,
                             gridspec_kw={"hspace": 0.16})
    a = d["a"][q]
    cmap = ps.density_cmap()
    lo, hi = 1.6, 5.5
    for ax, yv, ylim, ylabel in (
            (axes[0], d["e"][q], (0, 0.45), "eccentricity (0 = circle)"),
            (axes[1], d["incl"][q], (0, 36), "inclination to the ecliptic (deg)")):
        hh, xe, ye = np.histogram2d(a, yv, bins=[780, 300], range=[[lo, hi], list(ylim)])
        im = ax.pcolormesh(xe, ye, hh.T, cmap=cmap, norm=LogNorm(vmin=1, vmax=hh.max()),
                           rasterized=True, zorder=1)
        ax.set_xlim(lo, hi)
        ax.set_ylim(*ylim)
        ax.set_ylabel(ylabel)
        ax.grid(False)
        for p_, q_ in ((4, 1), (3, 1), (5, 2), (7, 3), (2, 1), (3, 2), (4, 3), (1, 1)):
            ax.axvline(op.resonance_a(p_, q_), color=ps.C_PREDICTED, lw=0.7, alpha=0.55,
                       zorder=2)
    cb = fig.colorbar(im, ax=axes, pad=0.012, fraction=0.025, aspect=40)
    cb.set_label("asteroids per cell", color=ps.INK_SECONDARY)
    cb.outline.set_visible(False)

    # Planet-crossing boundaries on the a-e plane.
    ax = axes[0]
    xs = np.linspace(lo, hi, 400)
    e_mars = op.e_mars_crossing(xs, MARS_APHELION_AU)
    e_jup = op.e_aphelion_crossing(xs, JUPITER_PERIHELION_AU)
    ax.plot(xs, e_mars, color=ps.INK_SECONDARY, lw=1.1, zorder=3)
    ok = e_jup < 0.45
    ax.plot(xs[ok], e_jup[ok], color=ps.INK_SECONDARY, lw=1.1, zorder=3)
    ps.annotate(ax, "above this line, the orbit\ncrosses Mars' orbit",
                xy=(2.3, float(op.e_mars_crossing(2.3, MARS_APHELION_AU))),
                xytext=(1.7, 0.42), relpos=(0.1, 0))
    ps.annotate(ax, "above this line, it\nreaches Jupiter's orbit",
                xy=(4.35, float(op.e_aphelion_crossing(4.35, JUPITER_PERIHELION_AU))),
                xytext=(4.45, 0.34), relpos=(0, 0.5))
    ax.set_title("(a) Eccentricity: the belt fills the space between the planets' orbits "
                 "— and no further", pad=20)
    labels_e = [("Hungarias", 1.93, 0.035), ("Cybeles", 3.48, 0.03),
                ("Hildas", 3.97, 0.02), ("Trojans", 5.2, 0.015)]
    for text, x, yv in labels_e:
        ax.text(x, yv, text, ha="center", va="bottom", fontsize=8, color=ps.INK)
    for p_, q_ in ((4, 1), (3, 1), (5, 2), (7, 3), (2, 1), (3, 2), (4, 3), (1, 1)):
        x = op.resonance_a(p_, q_)
        # 5:2 and 7:3 sit close together; nudge their labels apart.
        dx = {(5, 2): -0.025, (7, 3): 0.025}.get((p_, q_), 0.0)
        ax.text(x + dx, 0.453, f"{p_}:{q_}", ha="center", va="bottom",
                fontsize=7.5, color=ps.INK_SECONDARY)

    ax = axes[1]
    ax.set_xlabel("semi-major axis (AU)")
    ax.set_title("(b) Inclination: families form clumps; secular resonances carve the "
                 "inner edge")
    labels_i = [("Hungarias", 1.92, 25.0), ("Phocaeas", 2.37, 26.5), ("Trojans", 5.2, 31),
                ("Hildas", 3.97, 13.5)]
    for text, x, yv in labels_i:
        ax.text(x, yv, text, ha="center", va="bottom", fontsize=8, color=ps.INK)
    ps.annotate(ax, "inner edge: the ν6\nsecular resonance", xy=(2.29, 9.3),
                xytext=(2.02, 12.8), relpos=(1, 0.5), fontsize=8)
    ps.caption(fig, _mpc_credit(int(q.sum())) + " Violet lines: mean-motion resonances with "
               "Jupiter. The planet-crossing curves use Mars' aphelion (1.666 AU) and "
               "Jupiter's perihelion (4.951 AU). The clumps in (b) are collisional families, "
               "blurred here by osculating elements; figure 5 shows them in proper elements. "
               "ν6: where an asteroid's perihelion would precess at Saturn's rate.")
    return ps.save(fig, FIGURE_DIR / "04_belt_architecture")


# ---------------------------------------------------------------------------
# 05 -- families: osculating vs proper
# ---------------------------------------------------------------------------

def fig_families(res: dict, arr: dict, d: dict, proper: dict) -> list[str]:
    num = proper["number"]
    ok = (num > 0) & (num <= d["number"].max())
    idx = np.full(int(d["number"].max()) + 1, -1)
    idx[d["number"][d["number"] > 0]] = np.flatnonzero(d["number"] > 0)
    j = idx[num[ok]]
    have = j >= 0
    a_osc = d["a"][j[have]]
    s_osc = np.sin(np.radians(d["incl"][j[have]]))
    a_pr = proper["a"][ok][have]
    s_pr = proper["sini"][ok][have]
    n = int(have.sum())

    fig = plt.figure(figsize=(13, 6.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.52], wspace=0.08)
    rng = [[2.1, 3.3], [0.0, 0.34]]
    bins = [600, 340]
    h1, xe, ye = np.histogram2d(a_osc, s_osc, bins=bins, range=rng)
    h2, _, _ = np.histogram2d(a_pr, s_pr, bins=bins, range=rng)
    vmax = max(h1.max(), h2.max())
    cmap = ps.density_cmap()
    axes = []
    for k, (h, title) in enumerate((
            (h1, "(a) Osculating: each orbit as it is today"),
            (h2, "(b) Proper: each orbit averaged over millions of years"))):
        ax = fig.add_subplot(gs[0, k], sharey=axes[0] if axes else None)
        im = ax.pcolormesh(xe, ye, h.T, cmap=cmap, norm=LogNorm(vmin=1, vmax=vmax),
                           rasterized=True)
        ax.set_xlim(*rng[0])
        ax.set_ylim(*rng[1])
        ax.grid(False)
        ax.set_xlabel("semi-major axis (AU)")
        ax.set_title(title)
        for p_, q_ in ((3, 1), (5, 2), (7, 3), (2, 1)):
            ax.axvline(op.resonance_a(p_, q_), color=ps.C_PREDICTED, lw=0.7, alpha=0.5)
        axes.append(ax)
    axes[0].set_ylabel("sine of inclination")
    plt.setp(axes[1].get_yticklabels(), visible=False)

    # Family names on the proper panel, at each family's median position.
    fams = res["families"]["top"]
    offsets = {
        "Eos": (0.10, 0.045), "Hertha": (-0.17, 0.055), "Eunomia": (0.13, 0.05),
        "Vesta": (-0.13, 0.06), "Massalia": (-0.10, -0.018), "Koronis": (0.02, 0.05),
        "Astraea": (0.02, -0.035), "Themis": (0.03, 0.05), "Agnia": (0.04, -0.038),
        "Maria": (-0.12, 0.05), "Veritas": (0.08, 0.03), "Klumpkea": (0.06, 0.035),
    }
    for f in fams:
        dx, dy = offsets.get(f["name"], (0.05, 0.04))
        ps.annotate(axes[1], f["name"], xy=(f["a_median"], f["sini_median"]),
                    xytext=(f["a_median"] + dx, f["sini_median"] + dy),
                    relpos=(0.5, 0) if dy > 0 else (0.5, 1), fontsize=8,
                    color=ps.INK, ha="center")

    # (c) How much tighter each family is in proper elements.
    ax = fig.add_subplot(gs[0, 2])
    ordered = sorted([f for f in fams if not f["secular_resonant"]],
                     key=lambda f: f["sini_sharpening"])
    yy = np.arange(len(ordered))
    ax.hlines(yy, 1, [f["sini_sharpening"] for f in ordered], color=ps.C_CONTEXT, lw=2)
    ax.plot([f["sini_sharpening"] for f in ordered], yy, "o", color=ps.C_MEASURED, ms=7,
            mec=ps.SURFACE, mew=1.5)
    for yi, f in zip(yy, ordered):
        ax.text(f["sini_sharpening"] * 1.12, yi, f"{f['sini_sharpening']:.1f}×",
                va="center", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_yticks(yy, [f["name"] for f in ordered], fontsize=8)
    ax.set_xscale("log")
    ax.set_xlim(1, 30)
    ax.set_xticks([1, 2, 5, 10, 20], ["1", "2", "5", "10", "20"])
    ax.set_xlabel("spread in inclination,\nosculating ÷ proper")
    ax.set_title("(c) How much tighter", fontsize=10)
    ax.grid(axis="y", visible=False)
    ax.yaxis.tick_right()
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)

    ps.caption(fig, f"The same {n:,} numbered asteroids in both panels. (a) {_mpc_credit()} "
               "(b) Proper elements: AstDyS synthetic catalogue (Knežević & Milani), June "
               "2024; family names from the AstDyS classification (Milani et al. 2014, "
               "updated). (c) Interquartile range of sin i among each family's numbered "
               "members. Violet lines: the 3:1, 5:2, 7:3 and 2:1 resonances.")
    return ps.save(fig, FIGURE_DIR / "05_families")


# ---------------------------------------------------------------------------
# 06 -- where everything is today
# ---------------------------------------------------------------------------

PLANET_ORDER = ("Mercury", "Venus", "Earth", "Mars", "Jupiter")


def _planets_overlay(ax, planets, *, rotate_deg: float = 0.0, labels: bool = True):
    rot = np.radians(rotate_deg)
    c, s = np.cos(rot), np.sin(rot)
    for name in PLANET_ORDER:
        b = planets["bodies"][name]
        x, y = op.ellipse_xy(b["a_au"], b["e"], b["node_deg"] + b["peri_deg"])
        ax.plot(c * x - s * y, s * x + c * y, color=ps.C_PLANET, lw=0.7, alpha=0.55,
                zorder=4)
        px, py = b["xyz_au"][0], b["xyz_au"][1]
        rx, ry = c * px - s * py, s * px + c * py
        ax.plot(rx, ry, "o", color=ps.C_PLANET, ms=5 if name != "Jupiter" else 8,
                mec=ps.SURFACE, mew=1.5, zorder=6)
        if labels and name in ("Earth", "Mars", "Jupiter"):
            ax.text(rx + 0.18, ry + 0.18, name, fontsize=8, color=ps.INK, zorder=7)
    ax.plot(0, 0, marker="*", color=ps.INK, ms=11, mec=ps.SURFACE, mew=1, zorder=6)


def fig_snapshot(res: dict, arr: dict, d: dict) -> list[str]:
    snap = arr["snapshot"]
    planets = arr["planets"]
    x, y = snap["x"].astype(float), snap["y"].astype(float)
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 6.9), gridspec_kw={"wspace": 0.12})

    ax = axes[0]
    lim = 6.0
    hh, xe, ye = np.histogram2d(x, y, bins=700, range=[[-lim, lim], [-lim, lim]])
    ax.pcolormesh(xe, ye, hh.T, cmap=ps.density_cmap(), norm=LogNorm(vmin=1),
                  rasterized=True, zorder=1)
    _planets_overlay(ax, planets)
    J = res["resonant"]["jupiter"]
    r = np.hypot(*planets["bodies"]["Jupiter"]["xyz_au"][:2])
    for sign, name in ((+1, "L4"), (-1, "L5")):
        th = np.radians(J["theta"] + sign * 60)
        x0, y0 = r * np.cos(th), r * np.sin(th)
        # Each label goes just beyond its cloud, pushed clear of the frame.
        tx = float(np.clip(x0 * 1.13, -5.2, 4.9))
        ty = float(np.clip(y0 * 1.13 + (-2.6 if abs(x0) > 4 else 0), -5.8, 5.7))
        ps.annotate(ax, f"{name} Trojans", xy=(x0 * 1.02, y0 * 1.02), xytext=(tx, ty),
                    fontsize=8, relpos=(0.5, 1) if ty < y0 else (0.5, 0.5),
                    ha="center", va="center")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("x (AU), toward the vernal equinox")
    ax.set_ylabel("y (AU)")
    ax.grid(False)
    date = "29 September 2026"
    ax.set_title(f"(a) Every known minor planet, {date}")

    # (b) Hildas and Trojans, rotated so Jupiter is straight up.
    ax = axes[1]
    rot = 90.0 - J["theta"]
    c, s = np.cos(np.radians(rot)), np.sin(np.radians(rot))
    xr, yr = c * x - s * y, s * x + c * y
    q = arr["quality"]
    belt = q & (d["a"] > 2.0) & (d["a"] < 3.3)
    hb, xe, ye = np.histogram2d(xr[belt], yr[belt], bins=500, range=[[-lim, lim], [-lim, lim]])
    ax.pcolormesh(xe, ye, np.where(hb.T > 0, 1, np.nan), cmap=matplotlib.colors.ListedColormap(
        [ps.C_CONTEXT]), rasterized=True, zorder=1, alpha=0.55)
    a_, e_, i_ = d["a"], d["e"], d["incl"]
    from config import HILDA_A, HILDA_E_MAX, HILDA_I_MAX_DEG, TROJAN_A
    hil = q & (a_ > HILDA_A[0]) & (a_ < HILDA_A[1]) & (e_ < HILDA_E_MAX) & (i_ < HILDA_I_MAX_DEG)
    tro = q & (a_ > TROJAN_A[0]) & (a_ < TROJAN_A[1])
    ax.scatter(xr[hil], yr[hil], s=1.2, color=ps.C_HILDA, lw=0, alpha=0.7, rasterized=True,
               zorder=3, label=f"Hildas, in the 3:2 resonance ({hil.sum():,})")
    ax.scatter(xr[tro], yr[tro], s=1.2, color=ps.C_TROJAN, lw=0, alpha=0.7, rasterized=True,
               zorder=3, label=f"Trojans, sharing Jupiter's orbit ({tro.sum():,})")
    _planets_overlay(ax, planets, rotate_deg=rot, labels=False)
    rj = np.hypot(*planets["bodies"]["Jupiter"]["xyz_au"][:2])
    ax.text(0.35, rj + 0.12, "Jupiter", ha="left", va="bottom", fontsize=9, color=ps.INK)
    tr = res["resonant"]["trojans"]["L4_L5_by_H"][0]
    for xt, name, n in ((-5.3, "L4 Trojans, leading", tr["L4"]),
                        (5.3, "L5 Trojans, trailing", tr["L5"])):
        ax.text(xt, 5.35, f"{name}\n{n:,}", ha="center", va="center", fontsize=8,
                color=ps.INK_SECONDARY)
    ax.annotate("", xy=(-2.9, 5.35), xytext=(-1.1, 5.35),
                arrowprops={"arrowstyle": "->", "color": ps.INK_SECONDARY, "lw": 0.9})
    ax.text(-2.0, 5.5, "Jupiter's motion", ha="center", va="bottom", fontsize=7.5,
            color=ps.INK_SECONDARY)
    ps.annotate(ax, "the Hildas crowd into a\ntriangle, its corners\nopposite and 60° either\n"
                "side of Jupiter", xy=(0.4, -4.35), xytext=(1.9, -5.55), relpos=(0, 0.5),
                fontsize=7.5)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.grid(False)
    ax.legend(loc="lower left", fontsize=7.5, markerscale=6, handletextpad=0.3)
    ax.set_title("(b) Hildas and Trojans, turned so Jupiter is at the top")
    ps.caption(fig, _mpc_credit() + " Each orbit carried by two-body motion from its epoch "
               f"to {date} 0h TT (a median of {snap['dt_days_median']:.0f} days); planet "
               "positions and orbits from JPL Horizons. Viewed from ecliptic north; the "
               "planets and asteroids circle anticlockwise. Grey in (b): the main belt.")
    return ps.save(fig, FIGURE_DIR / "06_snapshot")


# ---------------------------------------------------------------------------
# 07 -- resonance lock
# ---------------------------------------------------------------------------

def fig_resonance_lock(res: dict, arr: dict) -> list[str]:
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.8), gridspec_kw={"wspace": 0.38})
    edges = np.linspace(-180, 180, 37)
    hil = res["resonant"]["hildas"]

    ax = axes[0]
    h1, _ = np.histogram(arr["hilda_sigma"], edges)
    h0, _ = np.histogram(arr["control_sigma"], edges)
    ax.stairs(h0 / h0.sum() * 100, edges, color=ps.MUTED, lw=1.2,
              label=f"Cybele zone, not in resonance ({hil['control_n']:,})")
    ax.stairs(h1 / h1.sum() * 100, edges, fill=True, color=ps.C_HILDA, alpha=0.2, lw=0)
    ax.stairs(h1 / h1.sum() * 100, edges, color=ps.C_HILDA, lw=1.6,
              label=f"Hildas ({hil['n']:,})")
    ax.set_xlim(-180, 180)
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_xlabel(r"resonant angle $\sigma = 3\lambda_J - 2\lambda - \varpi$ (deg)")
    ax.set_ylabel("% of asteroids per 10°")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_title("(a) The Hildas' 3:2 lock, measured")
    ax.text(0.97, 0.62, f"{hil['frac_sigma_within_90']:.1%} of Hildas\nwithin ±90° of 0,\n"
            f"vs {hil['control_frac_sigma_within_90']:.1%} for\nthe control group",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
            color=ps.INK_SECONDARY)
    ax.set_ylim(0, (h1 / h1.sum() * 100).max() * 1.35)

    ax = axes[1]
    tro = res["resonant"]["trojans"]
    ht, _ = np.histogram(arr["trojan_dlam"], edges)
    ax.stairs(ht, edges, fill=True, color=ps.C_TROJAN, alpha=0.2, lw=0)
    ax.stairs(ht, edges, color=ps.C_TROJAN, lw=1.6)
    ax.set_xlim(-180, 180)
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
    ax.set_xlabel(r"longitude measured from Jupiter's, $\lambda - \lambda_J$ (deg)")
    ax.set_ylabel("Trojans per 10°")
    ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter("{x:,.0f}"))
    all_ = tro["L4_L5_by_H"][0]
    for x0, name, n in ((60, "L4 (ahead)", all_["L4"]), (-60, "L5 (behind)", all_["L5"])):
        k = np.searchsorted(edges, x0) - 1
        ax.text(x0, ht[max(k - 1, 0):k + 2].max() * 1.04, f"{name}\n{n:,}", ha="center",
                va="bottom", fontsize=8, color=ps.INK_SECONDARY)
    ax.set_ylim(0, ht.max() * 1.35)
    ax.set_title("(b) Trojans: 60° ahead and 60° behind")

    ax = axes[2]
    # Left to right: everything, then ever brighter (bigger, more complete).
    rows = [tro["L4_L5_by_H"][0]] + [r for r in tro["L4_L5_by_H"] if r["H_max"] < 50]
    xs = np.arange(len(rows))
    lit, lit_err = tro["literature_ratio"], tro["literature_ratio_err"]
    ax.axhspan(lit - lit_err, lit + lit_err, color=ps.C_PREDICTED, alpha=0.15, lw=0)
    ax.axhline(lit, color=ps.C_PREDICTED, lw=1.4,
               label="Vokrouhlický et al. 2024,\nbias-corrected, H < 15")
    ax.errorbar(xs, [r["ratio"] for r in rows], yerr=[r["ratio_err"] for r in rows], fmt="o",
                color=ps.C_MEASURED, ms=7, mec=ps.SURFACE, mew=1.5, elinewidth=1.4,
                capsize=0, zorder=4, label="this project (raw counts)")
    ax.set_xticks(xs, [f"H < {r['H_max']:g}" if r["H_max"] < 50 else "all" for r in rows],
                  fontsize=8)
    ax.set_xlim(-0.5, len(rows) - 0.5)
    ax.set_ylim(1.0, 2.0)
    ax.set_ylabel("L4 / L5 number ratio")
    ax.set_xlabel("which Trojans are counted (only bigger ones to the right)")
    ax.legend(loc="upper right", fontsize=7.5)
    ax.set_title("(c) More Trojans lead than trail")
    ax.text(0.5, rows[0]["ratio"] - 0.07, "faint Trojans: raw ratio\ninflated by discovery bias",
            ha="center", va="top", fontsize=7.5, color=ps.INK_SECONDARY)
    ps.caption(fig, _mpc_credit() + " Mean longitudes λ and longitudes of perihelion ϖ at "
               "2026 Sep 29.0 TT; Jupiter's from JPL Horizons. Hildas: 3.7 < a < 4.2 AU, "
               "e < 0.3, i < 20°; control: 3.3 < a < 3.6 AU; Trojans: 4.95 < a < 5.45 AU. "
               "Error bars: binomial, 1σ.")
    return ps.save(fig, FIGURE_DIR / "07_resonance_lock")


# ---------------------------------------------------------------------------
# 08 -- the simulation
# ---------------------------------------------------------------------------

def fig_simulation(res: dict, sim: dict, arr: dict) -> list[str]:
    a0 = sim["ic_a"]
    t = sim["t_years"]
    removed = np.isfinite(sim["t_removed"])
    fig = plt.figure(figsize=(13, 8.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.45, wspace=0.22)

    # (a) initial vs surviving, against the real belt.
    ax = fig.add_subplot(gs[0, :])
    w = 0.01
    edges = np.arange(2.0, 3.4 + w / 2, w)
    n0, _ = np.histogram(a0, edges)
    n1, _ = np.histogram(a0[~removed], edges)
    ax.stairs(n0, edges, color=ps.MUTED, lw=1.1, label="at the start: spread evenly")
    ax.stairs(n1, edges, fill=True, color=ps.C_DATA, alpha=0.2, lw=0)
    ax.stairs(n1, edges, color=ps.C_DATA, lw=1.5,
              label=f"still there after {t[-1] / 1e6:g} million years")
    ax.set_xlim(2.0, 3.4)
    ax.set_ylim(0, n0.max() * 1.5)
    ax.set_xlabel("starting semi-major axis (AU)")
    ax.set_ylabel(f"test particles per {w} AU")
    _ratio_axis(ax, [Fraction(4), Fraction(3), Fraction(5, 2), Fraction(7, 3), Fraction(2)])
    _res_ticks(ax, [(4, 1), (3, 1), (5, 2), (7, 3), (2, 1)])
    ax.legend(loc="upper right", fontsize=8, ncol=2)
    ax.set_title(f"(a) SIMULATION: {len(a0):,} test particles, Jupiter and Saturn — "
                 "the gaps open by themselves")

    # (b) eccentricity histories.
    ax = fig.add_subplot(gs[1, 0])
    picks = _example_particles(sim)
    colors = [ps.ORANGE, ps.BLUE, ps.AQUA]
    for (k, text), col in zip(picks, colors):
        e = sim["e"][:, k].astype(float)
        good = np.isfinite(e) & (t > 0)
        ax.plot(t[good] / 1e3, e[good], color=col, lw=1.3, label=text)
        if removed[k]:
            last = np.flatnonzero(good)[-1]
            ax.plot(t[last] / 1e3, e[last], "x", color=col, ms=8, mew=2)
    from matplotlib.lines import Line2D
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Line2D([], [], color=ps.INK_SECONDARY, marker="x", ls="none", ms=7, mew=2))
    labels.append("removed: perihelion inside Mars' orbit")
    ax.set_xscale("log")
    ax.set_xlabel("time (thousands of years, log scale)")
    ax.set_ylabel("eccentricity")
    ax.set_xlim(1, t[-1] / 1e3)
    ax.set_xticks([1, 10, 100, 1000], ["1", "10", "100", "1,000"])
    ax.set_ylim(0, 0.75)
    ax.legend(handles, labels, loc="upper left", fontsize=7.5)
    ax.set_title("(b) Inside a resonance, Jupiter pumps the orbit's eccentricity")

    # (c) survival fraction.
    ax = fig.add_subplot(gs[1, 1])
    w = 0.02
    edges = np.arange(2.0, 3.4 + w / 2, w)
    n0, _ = np.histogram(a0, edges)
    n1, _ = np.histogram(a0[~removed], edges)
    frac = np.where(n0 > 0, n1 / np.maximum(n0, 1), np.nan)
    ax.stairs(frac * 100, edges, color=ps.C_DATA, lw=1.5)
    ax.stairs(frac * 100, edges, fill=True, color=ps.C_DATA, alpha=0.15, lw=0)
    for p_, q_ in ((4, 1), (3, 1), (5, 2), (7, 3), (2, 1)):
        x = op.resonance_a(p_, q_)
        ax.axvline(x, color=ps.C_PREDICTED, lw=0.9, alpha=0.7)
        ax.text(x, 103, f"{p_}:{q_}", ha="center", va="bottom", fontsize=8,
                color=ps.INK_SECONDARY)
    ax.set_xlim(2.0, 3.4)
    ax.set_ylim(0, 100)
    ax.set_xlabel("starting semi-major axis (AU)")
    ax.set_ylabel("% surviving")
    ps.annotate(ax, "2:1: not yet cleared —\nescape needs e > 0.5",
                xy=(op.resonance_a(2, 1), 92), xytext=(3.02, 30), relpos=(0.5, 1),
                fontsize=7.5, ha="center")
    ps.annotate(ax, "inner edge: Mars-crossing\nand the ν6 resonance", xy=(2.12, 10),
                xytext=(2.36, 12), relpos=(0, 0.5), fontsize=7.5, ha="left")
    ax.set_title(f"(c) Survival after {t[-1] / 1e6:g} Myr", pad=18)
    s = res.get("simulation") or {}
    ps.caption(fig, "A SIMULATION, not data. Sun, Jupiter and Saturn start from their JPL "
               "Horizons positions on 2026 Sep 29; each test particle takes the "
               "eccentricity, inclination and orientation of a real asteroid (MPCORB) but a "
               "semi-major axis drawn evenly from 2.0–3.4 AU. REBOUND WHFast integrator, "
               f"0.1-year steps; planets' energy conserved to {s.get('energy_drift_max', 0):.0e}. "
               "A particle is removed when its perihelion falls inside Mars' aphelion (Mars "
               "itself is not simulated).")
    return ps.save(fig, FIGURE_DIR / "08_simulation")


def _example_particles(sim: dict) -> list[tuple[int, str]]:
    """Three illustrative histories: removed from the 3:1, removed from the
    5:2, and a survivor far from any strong resonance."""
    a0 = sim["ic_a"]
    t_rem = sim["t_removed"]
    out = []
    for (p, q), text in (((3, 1), "started in the 3:1"), ((5, 2), "started in the 5:2")):
        c = op.resonance_a(p, q)
        cand = np.flatnonzero((np.abs(a0 - c) < 0.006) & np.isfinite(t_rem))
        if len(cand):
            # The one removed closest to the median removal time: typical, not extreme.
            med = np.median(t_rem[cand])
            k = cand[np.argmin(np.abs(t_rem[cand] - med))]
            out.append((int(k), f"{text} (a = {a0[k]:.3f} AU)"))
    safe = np.flatnonzero((np.abs(a0 - 2.65) < 0.02) & ~np.isfinite(t_rem))
    if len(safe):
        k = safe[0]
        out.append((int(k), f"started between gaps (a = {a0[k]:.3f} AU)"))
    return out


def make_all(res: dict, arr: dict, d: dict, proper: dict | None, sim: dict | None) -> list[str]:
    written = []
    written += fig_kirkwood_gaps(res, arr)
    written += fig_gap_measurements(res, arr)
    written += fig_is_it_real(res, arr, d)
    written += fig_architecture(res, arr, d)
    if proper is not None and "families" in res:
        written += fig_families(res, arr, d, proper)
    written += fig_snapshot(res, arr, d)
    written += fig_resonance_lock(res, arr)
    if sim is not None:
        written += fig_simulation(res, sim, arr)
    return written
