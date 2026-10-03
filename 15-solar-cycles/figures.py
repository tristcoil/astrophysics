"""Every figure, drawn from the analysis tables and the results dictionary.

Figures are numbered in reading order:

    01_sunspot_record     275 years of sunspots, 25 cycles -- the keeper
    02_butterfly          where the spots are: 152 years of sunspot latitudes
    03_zone_drift         the equatorward drift of the sunspot zones (Sporer's law)
    04_cycle_rules        cycle lengths, the Waldmeier effect, Gnevyshev-Ohl
    05_cycle25            the current cycle against its official forecast
    06_rotation           the Sun's differential rotation, from tracked spot groups
    07_four_centuries     the Maunder minimum and the group-number reconstructions
    08_caveats            where the record is thin, and where its seams are
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import gaussian_kde

import plotstyle as ps
import solar_physics as sp
from config import (
    DALTON_MINIMUM,
    FIGURE_DIR,
    LIT_BALTHASAR_1986,
    LIT_HATHAWAY_DRIFT,
    LIT_PANEL_2019,
    MANIFEST_PATH,
    MAUNDER_MINIMUM,
    RGO_LAST_YEAR,
    USAF_AREA_FACTOR_HATHAWAY,
)

ps.use_style()


def _manifest() -> dict:
    try:
        return json.loads(MANIFEST_PATH.read_text())
    except (OSError, ValueError):
        return {}


def _credit_sn() -> str:
    m = _manifest()
    when = f", retrieved {m['retrieved_utc'][:10]}" if m else ""
    return f"Sunspot number: WDC-SILSO, Royal Observatory of Belgium, Brussels (version 2.0{when})."


def _credit_groups() -> str:
    return ("Sunspot groups: Royal Greenwich Observatory 1874–1976, USAF/NOAA 1977–, "
            "compiled by D. H. Hathaway (solarcyclescience.com).")


def _t(ym: str) -> float:
    y, m = map(int, ym.split("-"))
    return y + (m - 0.5) / 12.0


_MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def _month_name(ym: str) -> str:
    y, m = map(int, ym.split("-"))
    return f"{_MONTHS[m - 1]} {y}"


def _lat_ticks(ax, top: int = 40, step: int = 20):
    ticks = np.arange(-top, top + 1, step)
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{abs(v)}°{'N' if v > 0 else 'S' if v < 0 else ''}" for v in ticks])


def _overlap(a, b) -> bool:
    return not (a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1])


#: Offsets (points) and alignments tried, in order, when labelling a point.
_SPOTS = [(5, 3, "left", "bottom"), (5, -3, "left", "top"), (-5, 3, "right", "bottom"),
          (-5, -3, "right", "top"), (0, 7, "center", "bottom"), (0, -7, "center", "top"),
          (8, 0, "left", "center"), (-8, 0, "right", "center"), (12, 10, "left", "bottom"),
          (-12, 10, "right", "bottom"), (12, -10, "left", "top"), (-12, -10, "right", "top")]


def _label_points(ax, xs, ys, texts, *, marker_pts: float = 3.5, fontsize: float = 7.5,
                  avoid=()) -> None:
    """Label scatter points without letting labels collide with each other, with
    any point, or with the extra boxes in ``avoid`` (display coordinates)."""
    fig = ax.figure
    renderer = fig.canvas.get_renderer()
    pts = ax.transData.transform(np.column_stack([xs, ys]))
    r = marker_pts * fig.dpi / 72.0
    taken = [(x - r, y - r, x + r, y + r) for x, y in pts] + list(avoid)
    frame = ax.get_window_extent(renderer)
    for x, y, text in zip(xs, ys, texts):
        chosen = None
        for dx, dy, ha, va in _SPOTS:
            t = ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points", ha=ha, va=va,
                            fontsize=fontsize, color=ps.INK_SECONDARY, zorder=6)
            bb = t.get_window_extent(renderer)
            box = (bb.x0, bb.y0, bb.x1, bb.y1)
            inside = bb.x0 >= frame.x0 and bb.x1 <= frame.x1 and bb.y0 >= frame.y0 and bb.y1 <= frame.y1
            if inside and not any(_overlap(box, b) for b in taken):
                chosen = box
                break
            t.remove()
        if chosen is None:  # nowhere clean: take the first spot anyway
            dx, dy, ha, va = _SPOTS[0]
            t = ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points", ha=ha, va=va,
                            fontsize=fontsize, color=ps.INK_SECONDARY, zorder=6)
            bb = t.get_window_extent(renderer)
            chosen = (bb.x0, bb.y0, bb.x1, bb.y1)
        taken.append(chosen)


def _text_box(ax, artist):
    bb = artist.get_window_extent(ax.figure.canvas.get_renderer())
    return (bb.x0, bb.y0, bb.x1, bb.y1)


# --- 01 The record ---------------------------------------------------------------------------------


def figure_record(tables, results) -> list:
    m = tables["monthly"]
    yearly = tables["yearly"]
    cyc = tables["cycles"]
    fig, ax = plt.subplots(figsize=(13.6, 5.6))
    top = 430

    ax.axvspan(*DALTON_MINIMUM, color=ps.WASH, lw=0, zorder=0)
    ax.text(np.mean(DALTON_MINIMUM), 176, "Dalton minimum:\nthe two weakest\ncycles on record",
            ha="center", va="bottom", fontsize=8, color=ps.INK_SECONDARY)

    early = yearly[yearly["t"] < 1749.5]
    ax.step(early["t"], early["sn"], where="mid", color=ps.INK_SECONDARY, lw=0.9, alpha=0.7, zorder=2)
    ax.plot(m["t"], m["sn"], color=ps.C_MONTHLY, lw=0.45, alpha=0.6, zorder=2)
    ax.plot(m["t"], m["smoothed"], color=ps.C_SMOOTHED, lw=1.8, zorder=3)

    for c in cyc.itertuples():
        ax.text(c.t_max, c.max_sn + 14 if c.cycle != 3 else 280, str(c.cycle),
                ha="center", va="bottom", fontsize=7.5, color=ps.INK_SECONDARY, zorder=4)

    big = cyc.loc[cyc["max_sn"].idxmax()]
    ps.annotate(ax, f"Cycle 19, the strongest:\n{big['max_sn']:.0f} ({_month_name(big['max_ym'])})",
                (big["t_max"], big["max_sn"]), (1916, 372), relpos=(1, 0.5), ha="right", va="center")
    c25 = cyc.iloc[-1]
    ax.scatter([c25["t_max"]], [c25["max_sn"]], s=34, color=ps.C_FOCUS, edgecolor=ps.SURFACE,
               linewidth=1.5, zorder=5)
    ps.annotate(ax, f"Cycle 25 peak:\n{c25['max_sn']:.1f} ({_month_name(c25['max_ym'])})",
                (c25["t_max"], c25["max_sn"]), (2003, 318), relpos=(0.5, 0), ha="center", va="bottom")

    ax.set_xlim(1698, 2029)
    ax.set_ylim(0, top)
    ax.set_xticks(np.arange(1700, 2026, 25))
    ax.set_ylabel("sunspot number")
    ax.set_title("275 years of sunspots: the Sun's activity cycle, counted since 1749")
    handles = [Line2D([], [], color=ps.C_MONTHLY, lw=1.0, alpha=0.8, label="monthly mean"),
               Line2D([], [], color=ps.C_SMOOTHED, lw=1.8, label="13-month smoothed (cycles numbered)"),
               Line2D([], [], color=ps.INK_SECONDARY, lw=0.9, alpha=0.7, drawstyle="steps-mid",
                      label="yearly means before 1749 (Wolf's reconstruction)")]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 1.0))
    r = results["record"]
    ps.caption(fig, f"{_credit_sn()} {r['n_cycles']} cycles; mean length "
               f"{r['lengths']['mean_yr']:.1f} ± {r['lengths']['sd_yr']:.1f} years. Values before 1818 rest on "
               "sparse observations (no daily record exists before then); the last five months are provisional.")
    return ps.save(fig, FIGURE_DIR / "01_sunspot_record")


# --- 02 The butterfly diagram --------------------------------------------------------------------------


def figure_butterfly(tables, results) -> list:
    grid = tables["butterfly_grid"]
    te = tables["butterfly_t_edges"]
    le = tables["butterfly_lat_edges"]
    cyc = tables["cycles"]
    mon = tables["hemispheric_area_monthly"]

    fig = plt.figure(figsize=(13.6, 8.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[3.3, 1.3], width_ratios=[1, 0.011],
                          hspace=0.07, wspace=0.012)
    ax = fig.add_subplot(gs[0, 0])
    axa = fig.add_subplot(gs[1, 0], sharex=ax)
    cax = fig.add_subplot(gs[0, 1])

    z = np.ma.masked_invalid(np.where(grid > 0, grid, np.nan)).T
    norm = LogNorm(vmin=0.1, vmax=1000.0, clip=True)
    mesh = ax.pcolormesh(te, le, z, cmap=ps.SUN_CMAP, norm=norm, shading="flat", rasterized=True)
    ax.axhline(0, color=ps.AXIS, lw=0.6, zorder=1)
    ax.grid(False)
    ax.set_ylim(-50, 50)
    _lat_ticks(ax, 40, 20)
    ax.set_ylabel("latitude on the Sun")
    for c in cyc[cyc["cycle"] >= 12].itertuples():
        ax.text(c.t_max, 51.5, str(c.cycle), ha="center", va="bottom", fontsize=7.5,
                color=ps.INK_SECONDARY)
    ax.text(1874.2, 51.5, "cycle", ha="left", va="bottom", fontsize=7.5, color=ps.MUTED)
    ax.set_title("The butterfly diagram: every sunspot group since May 1874, by latitude", pad=22)
    plt.setp(ax.get_xticklabels(), visible=False)
    cb = fig.colorbar(mesh, cax=cax)
    cb.set_label("spot area per 1° latitude strip, mean per day over one solar rotation\n"
                 "(millionths of the Sun's visible hemisphere)", fontsize=8, color=ps.INK_SECONDARY)
    cb.ax.tick_params(labelsize=7.5, colors=ps.MUTED)
    cb.outline.set_edgecolor(ps.AXIS)

    t = mon["t"].to_numpy()
    axa.fill_between(t, 0, mon["north"].to_numpy() / 1e4, color=ps.SUN_RAMP[2], lw=0, step="mid")
    axa.fill_between(t, 0, -mon["south"].to_numpy() / 1e4, color=ps.SUN_RAMP[2], lw=0, step="mid")
    axa.plot(t, mon["north_s"] / 1e4, color=ps.SUN_RAMP[6], lw=1.4)
    axa.plot(t, -mon["south_s"] / 1e4, color=ps.SUN_RAMP[6], lw=1.4)
    axa.axhline(0, color=ps.AXIS, lw=0.6)
    lim = 0.45
    axa.set_ylim(-lim, lim)
    ticks = np.array([-0.4, -0.2, 0, 0.2, 0.4])
    axa.set_yticks(ticks)
    axa.set_yticklabels([f"{abs(v):.1f}%" for v in ticks])
    axa.text(1875, lim * 0.92, "north", va="top", fontsize=8, color=ps.INK_SECONDARY)
    axa.text(1875, -lim * 0.92, "south", va="bottom", fontsize=8, color=ps.INK_SECONDARY)
    axa.text(2027, lim * 0.92, "monthly mean (light), 13-month smoothed (dark)", ha="right", va="top",
             fontsize=7.5, color=ps.MUTED)
    axa.set_ylabel("spot area\n(% of hemisphere)", fontsize=8.5)
    axa.set_xlim(1873.5, 2027.5)
    axa.set_xticks(np.arange(1880, 2021, 10))

    seam = RGO_LAST_YEAR + 1
    for a in (ax, axa):
        a.axvline(seam, color=ps.INK_SECONDARY, lw=0.6, alpha=0.5, zorder=3)
    axa.text(seam - 0.6, -lim * 0.92, "Royal Greenwich Observatory", ha="right", va="bottom",
             fontsize=7, color=ps.INK_SECONDARY)
    axa.text(seam + 0.6, -lim * 0.92, "US Air Force / NOAA", ha="left", va="bottom",
             fontsize=7, color=ps.INK_SECONDARY)

    cen = results["butterfly"]["census"]
    ps.caption(fig, f"{_credit_groups()} {cen['records']:,} daily records of {cen['groups']:,} groups on "
               f"{cen['days_observed']:,} observed days. After 1976 the areas are multiplied by "
               f"{USAF_AREA_FACTOR_HATHAWAY} to match RGO's scale, and latitudes are whole degrees. "
               f"{cen['records_beyond_50deg']} records beyond 50° are not shown.", y=0.02)
    return ps.save(fig, FIGURE_DIR / "02_butterfly")


# --- 03 The drift of the sunspot zones ------------------------------------------------------------------


def figure_zone_drift(tables, results) -> list:
    cen = tables["centroids"]
    d = results["butterfly"]["drift_law"]
    fit = d["fit"]
    fig, axes = plt.subplots(1, 3, figsize=(14.4, 5.0),
                             gridspec_kw={"width_ratios": [1.55, 1.0, 1.0], "wspace": 0.3})

    ax = axes[0]
    yrs = cen["dt_months"] / 12.0
    ax.scatter(yrs, cen["lat"], s=5, color=ps.MUTED, alpha=0.25, lw=0, zorder=1, rasterized=True)
    bins = np.arange(-2.0, 14.01, 0.5)
    g = cen.groupby(pd.cut(yrs, bins), observed=True)["lat"].agg(["mean", "std", "size"])
    g = g[g["size"] >= 10]
    mids = np.array([iv.mid for iv in g.index])
    ax.errorbar(mids, g["mean"], yerr=g["std"] / np.sqrt(g["size"]), fmt="o", ms=4.5,
                color=ps.C_SMOOTHED, mec=ps.SURFACE, mew=1.0, elinewidth=1.0, capsize=0, zorder=4)
    tt = np.linspace(-2, 14, 300)
    ax.plot(tt, sp.drift_law(tt * 12, fit["lambda0_deg"], fit["tau_months"]), color=ps.C_FOCUS,
            lw=2.0, zorder=3)
    ax.plot(tt, sp.drift_law(tt * 12, LIT_HATHAWAY_DRIFT["lambda0_deg"], LIT_HATHAWAY_DRIFT["tau_months"]),
            color=ps.C_LITERATURE, lw=1.2, zorder=3)
    handles = [
        Line2D([], [], color=ps.C_FOCUS, lw=2.0,
               label=f"our fit: {fit['lambda0_deg']:.1f}° × exp(−t / {fit['tau_months']:.0f} months)"),
        Line2D([], [], color=ps.C_LITERATURE, lw=1.2, label="Hathaway (2011): 28° × exp(−t / 90 months)"),
        Line2D([], [], marker="o", ls="", ms=4.5, color=ps.C_SMOOTHED, mec=ps.SURFACE,
               label="mean of each half-year, ± its error"),
        Line2D([], [], marker="o", ls="", ms=3, color=ps.MUTED, alpha=0.6,
               label=f"one hemisphere, one solar rotation ({fit['n']:,})"),
    ]
    ax.legend(handles=handles, loc="upper right")
    ax.set_xlim(-2, 14)
    ax.set_ylim(0, 40)
    ax.set_xlabel("years since the cycle began (t₀, from fitting its sunspot numbers)")
    ax.set_ylabel("latitude of the sunspot zone (area-weighted centroid)")
    ax.set_title("Each cycle's spots march toward the equator")

    ax = axes[1]
    st = d["strength_dependence"]["first_3_years"]
    per = pd.DataFrame(st["per_cycle"])
    f = st["fit"]
    xx = np.linspace(95, 300, 10)
    ax.plot(xx, f["a"] + f["b"] * xx, color=ps.C_FOCUS, lw=1.6, zorder=2)
    ax.scatter(per["amplitude"], per["mean_lat_first_3yr"], s=38, color=ps.C_SMOOTHED,
               edgecolor=ps.SURFACE, linewidth=1.5, zorder=3)
    ax.set_xlim(95, 300)
    ax.set_ylim(17.8, 24.2)
    c = st["correlation"]
    note = ax.text(0.97, 0.04, f"r = {c['r']:.2f}  (p = {c['p']:.3f}, {c['n']} cycles)\n"
                   f"slope {f['b'] * 100:.1f}° per 100 in amplitude", transform=ax.transAxes, ha="right",
                   va="bottom", fontsize=8, color=ps.INK_SECONDARY)
    fig.canvas.draw()
    _label_points(ax, per["amplitude"], per["mean_lat_first_3yr"], [str(v) for v in per["cycle"]],
                  avoid=[_text_box(ax, note)])
    ax.set_xlabel("cycle amplitude\n(peak of the smoothed sunspot number)")
    ax.set_ylabel("mean spot latitude in the cycle's first 3 years")
    ax.set_title("Stronger cycles start higher")

    ax = axes[2]
    s = tables["spots"]
    cyc = tables["cycles"].set_index("cycle")
    first = s.groupby("track").agg(cycle=("cycle", "first"), t_first=("t", "min"), lat=("alat", "mean"))
    new_lat, old_lat = [], []
    for n in range(13, int(cyc.index.max()) + 1):
        tmin = cyc.at[n, "t_min"]
        near = first[(first["t_first"] > tmin - 1) & (first["t_first"] < tmin + 1)]
        new_lat.append(near.loc[near["cycle"] == n, "lat"])
        old_lat.append(near.loc[near["cycle"] == n - 1, "lat"])
    new_lat, old_lat = pd.concat(new_lat), pd.concat(old_lat)
    bins = np.arange(0, 42, 2)
    for vals, color in ((old_lat, ps.C_SOUTH), (new_lat, ps.C_NORTH)):
        ax.hist(vals, bins=bins, color=color, alpha=0.18, histtype="stepfilled", lw=0)
        ax.hist(vals, bins=bins, color=color, histtype="step", lw=1.6)
    ov = results["butterfly"]["overlap"]
    ymax = ax.get_ylim()[1] * 1.25
    ax.set_ylim(0, ymax)
    ax.text(12.5, ymax * 0.80, f"old cycle\n{len(old_lat)} groups,\nmedian {np.median(old_lat):.0f}°",
            fontsize=8, color=ps.INK_SECONDARY, va="top")
    ax.text(31.5, ymax * 0.62, f"new cycle\n{len(new_lat)} groups,\nmedian {np.median(new_lat):.0f}°",
            fontsize=8, color=ps.INK_SECONDARY, va="top")
    ax.text(0.98, 0.98, f"first new-cycle group: a median {ov['median_first_new_months_before_min']:.0f} "
            f"months\nbefore minimum; last old-cycle group: {ov['median_last_old_months_after_min']:.0f} after",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=ps.MUTED)
    ax.set_xlim(0, 40)
    ax.set_xlabel("latitude of groups appearing\nwithin a year of a minimum")
    ax.set_ylabel("sunspot groups (13 minima since 1890)")
    ax.set_title("At minimum, two cycles share the Sun")

    ps.caption(fig, f"{_credit_groups()} Groups are assigned to cycles by latitude and date "
               f"(where cycles overlap, by the drift law itself). Fit uncertainty from resampling whole "
               f"cycles: λ₀ {fit['lambda0_ci68'][0]:.1f}–{fit['lambda0_ci68'][1]:.1f}°, τ "
               f"{fit['tau_ci68'][0]:.0f}–{fit['tau_ci68'][1]:.0f} months (68%).", y=-0.04)
    return ps.save(fig, FIGURE_DIR / "03_zone_drift")


# --- 04 The rules of the cycle -----------------------------------------------------------------------------


def figure_rules(tables, results) -> list:
    cyc = tables["cycles"]
    ru = results["rules"]
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 9.4), gridspec_kw={"hspace": 0.42, "wspace": 0.22})
    fig.canvas.draw()

    # (a) lengths, split into rise and decline
    ax = axes[0, 0]
    x = cyc["cycle"].to_numpy()
    rise = cyc["rise_yr"].to_numpy()
    decline = cyc["decline_yr"].to_numpy()
    ax.bar(x, rise, width=0.62, color=ps.BLUE, lw=0, label="rise: minimum to maximum")
    ax.bar(x, np.nan_to_num(decline), bottom=rise, width=0.62, color=ps.ORANGE, lw=0,
           label="decline: maximum to next minimum")
    ax.hlines(rise[:-1], x[:-1] - 0.31, x[:-1] + 0.31, color=ps.SURFACE, lw=1.4)  # 2px gap
    L = results["record"]["lengths"]
    ax.axhline(L["mean_yr"], color=ps.INK_SECONDARY, lw=0.8)
    ax.text(25.9, L["mean_yr"], f"mean\n{L['mean_yr']:.1f} yr", fontsize=8, color=ps.INK_SECONDARY,
            ha="left", va="center", clip_on=False)
    ax.text(25, rise[-1] + 0.25, "so\nfar", fontsize=7, color=ps.MUTED, ha="center", va="bottom")
    ax.set_xticks(x)
    ax.tick_params(axis="x", labelsize=7)
    ax.set_xlim(0.3, 25.7)
    ax.set_ylim(0, 16.5)
    ax.set_xlabel("cycle")
    ax.set_ylabel("years")
    ax.legend(loc="upper right", ncol=2)
    ax.set_title(f"Cycles last {L['shortest']['yr']:.1f} to {L['longest']['yr']:.1f} years, and rise "
                 "faster than they fall")

    # (b) the Waldmeier effect
    ax = axes[0, 1]
    w = ru["waldmeier"]
    f = ru["waldmeier_fit"]
    xx = np.linspace(2.6, 7.2, 10)
    ax.plot(xx, f["a"] + f["b"] * xx, color=ps.C_FOCUS, lw=1.6, zorder=2)
    ax.scatter(cyc["rise_yr"], cyc["max_sn"], s=36, color=ps.C_SMOOTHED, edgecolor=ps.SURFACE,
               linewidth=1.5, zorder=3)
    ax.set_xlim(2.5, 7.4)
    ax.set_ylim(60, 300)
    note = ax.text(0.97, 0.96, f"r = {w['r']:.2f} (95%: {w['r_ci95'][0]:.2f} to {w['r_ci95'][1]:.2f}), "
                   f"all {w['n']} cycles\nfrom Cycle 10 on: r = {ru['waldmeier_from_cycle10']['r']:.2f}",
                   transform=ax.transAxes, ha="right", va="top", fontsize=8, color=ps.INK_SECONDARY)
    _label_points(ax, cyc["rise_yr"], cyc["max_sn"], [str(v) for v in cyc["cycle"]],
                  avoid=[_text_box(ax, note)])
    ax.set_xlabel("rise time, minimum to maximum (years)")
    ax.set_ylabel("amplitude (smoothed sunspot number at maximum)")
    ax.set_title("Waldmeier effect: strong cycles rise fast")

    # (c) a long cycle is followed by a weak one
    ax = axes[1, 0]
    c = ru["length_vs_next_amplitude"]
    f = ru["length_vs_next_fit"]
    lengths = cyc["length_yr"].to_numpy()[:-1]
    nxt = cyc["max_sn"].to_numpy()[1:]
    xx = np.linspace(8.7, 13.9, 10)
    ax.plot(xx, f["a"] + f["b"] * xx, color=ps.C_FOCUS, lw=1.6, zorder=2)
    ax.scatter(lengths, nxt, s=36, color=ps.C_SMOOTHED, edgecolor=ps.SURFACE, linewidth=1.5, zorder=3)
    ax.set_xlim(8.6, 14.0)
    ax.set_ylim(60, 300)
    note = ax.text(0.97, 0.96, f"r = {c['r']:.2f} (95%: {c['r_ci95'][0]:.2f} to {c['r_ci95'][1]:.2f}), "
                   f"{c['n']} pairs", transform=ax.transAxes, ha="right", va="top", fontsize=8,
                   color=ps.INK_SECONDARY)
    _label_points(ax, lengths, nxt, [f"{n}→{n + 1}" for n in range(1, len(lengths) + 1)],
                  avoid=[_text_box(ax, note)], fontsize=7)
    ax.set_xlabel("length of a cycle (years)")
    ax.set_ylabel("amplitude of the next cycle")
    ax.set_title("A long cycle is followed by a weak one")

    # (d) Gnevyshev-Ohl
    ax = axes[1, 1]
    go = ru["gnevyshev_ohl"]
    for k, p in enumerate(go["pairs"]):
        color = ps.MUTED if p["holds"] else ps.C_FOCUS
        ax.plot([k - 0.17, k + 0.17], [p["even_total"] / 1e3, p["odd_total"] / 1e3], color=color, lw=1.4,
                zorder=2)
        ax.scatter([k - 0.17], [p["even_total"] / 1e3], s=30, color=ps.MUTED, edgecolor=ps.SURFACE,
                   linewidth=1.2, zorder=3)
        ax.scatter([k + 0.17], [p["odd_total"] / 1e3], s=34,
                   color=ps.C_SMOOTHED if p["odd_complete"] else ps.SURFACE,
                   edgecolor=ps.C_SMOOTHED, linewidth=1.3, zorder=4)
        if not p["holds"] and p["odd_complete"]:
            ax.text(k, max(p["even_total"], p["odd_total"]) / 1e3 + 0.6, "broken", ha="center",
                    fontsize=7.5, color=ps.INK_SECONDARY)
        if not p["odd_complete"]:
            ax.text(k + 0.17, p["odd_total"] / 1e3 + 0.7, "so far", ha="center", fontsize=7, color=ps.MUTED)
    ax.set_xticks(range(len(go["pairs"])))
    ax.set_xticklabels([f"{p['even']}→{p['odd']}" for p in go["pairs"]], fontsize=7.5)
    ax.set_xlim(-0.6, len(go["pairs"]) - 0.4)
    ax.set_ylim(0, 19)
    ax.set_xlabel("cycle pair (even → odd)")
    ax.set_ylabel("sunspot number summed over the cycle (thousands)")
    handles = [Line2D([], [], marker="o", ls="", color=ps.MUTED, label="even cycle"),
               Line2D([], [], marker="o", ls="", color=ps.C_SMOOTHED, label="the odd cycle after it"),
               Line2D([], [], color=ps.C_FOCUS, lw=1.4, label="rule broken")]
    ax.legend(handles=handles, loc="lower right")
    ax.set_title(f"Gnevyshev–Ohl rule: odd beats even in {go['n_hold']} of {go['n_complete']} pairs")

    ps.caption(fig, f"{_credit_sn()} Cycles run minimum to minimum of the 13-month smoothed series; "
               "Cycle 25 is not over. With so few cycles, correlations are fragile: their 95% intervals come "
               "from resampling cycles.", y=0.035)
    return ps.save(fig, FIGURE_DIR / "04_cycle_rules")


# --- 05 Cycle 25 --------------------------------------------------------------------------------------------


def figure_cycle25(tables, results) -> list:
    m = tables["monthly"]
    cyc = tables["cycles"]
    c25 = cyc.iloc[-1]
    r = results["cycle25"]
    panel = tables["panel"]
    cur = tables["current_forecast"]
    avg = tables["average_cycle"]

    fig = plt.figure(figsize=(13.6, 9.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], hspace=0.38, wspace=0.2)
    ax = fig.add_subplot(gs[0, :])
    axh = fig.add_subplot(gs[1, 0])
    axs = fig.add_subplot(gs[1, 1])

    sel = m["t"] >= 2018
    t_avg = c25["t_min"] + avg["months"].to_numpy() / 12.0
    ax.fill_between(t_avg, avg["mean"] - avg["sd"], avg["mean"] + avg["sd"], color=ps.WASH, lw=0, zorder=0)
    ax.plot(t_avg, avg["mean"], color=ps.AXIS, lw=1.4, zorder=1)
    tp = sp.decimal_year(panel["date"].to_numpy())
    ok = panel["envelope_low"].notna()
    ax.fill_between(tp[ok], panel.loc[ok, "envelope_low"], panel.loc[ok, "envelope_high"],
                    color=ps.C_LITERATURE, alpha=0.12, lw=0, zorder=1)
    ax.plot(tp, panel["ssn"], color=ps.C_LITERATURE, lw=1.6, zorder=2)
    tc = sp.decimal_year(cur["date"].to_numpy())
    ax.plot(tc, cur["predicted_ssn"], color=ps.C_LITERATURE, lw=1.1, alpha=0.9, zorder=2, ls=(0, (1, 1.5)))
    ax.plot(m.loc[sel, "t"], m.loc[sel, "sn"], color=ps.C_MONTHLY, lw=0.8, alpha=0.85, zorder=2)
    ax.plot(m.loc[sel, "t"], m.loc[sel, "smoothed"], color=ps.C_SMOOTHED, lw=2.2, zorder=3)

    ax.scatter([c25["t_max"]], [c25["max_sn"]], s=40, color=ps.C_SMOOTHED, edgecolor=ps.SURFACE,
               linewidth=1.5, zorder=5)
    ps.annotate(ax, f"observed peak: {c25['max_sn']:.1f}, {_month_name(c25['max_ym'])}",
                (c25["t_max"], c25["max_sn"]), (2025.6, 232), relpos=(0, 0.5), ha="left", va="center")
    tpk = _t(LIT_PANEL_2019["max_date"])
    ax.scatter([tpk], [LIT_PANEL_2019["max_ssn"]], s=36, color=ps.C_LITERATURE, edgecolor=ps.SURFACE,
               linewidth=1.5, zorder=5)
    ps.annotate(ax, f"forecast peak: {LIT_PANEL_2019['max_ssn']:.0f}, {_month_name(LIT_PANEL_2019['max_date'])}",
                (tpk, LIT_PANEL_2019["max_ssn"]), (2027.0, 172), relpos=(0, 0.5), ha="left", va="center")
    handles = [
        Line2D([], [], color=ps.C_MONTHLY, lw=1.0, label="monthly mean"),
        Line2D([], [], color=ps.C_SMOOTHED, lw=2.2, label="13-month smoothed"),
        (Patch(color=ps.C_LITERATURE, alpha=0.18, lw=0), Line2D([], [], color=ps.C_LITERATURE, lw=1.6)),
        Line2D([], [], color=ps.C_LITERATURE, lw=1.1, ls=(0, (1, 1.5))),
        (Patch(color=ps.WASH, lw=0), Line2D([], [], color=ps.AXIS, lw=1.4)),
    ]
    labels = ["monthly mean", "13-month smoothed",
              "official forecast, NOAA/NASA panel (Dec 2019), with its ±10%, ±6-month range",
              "NOAA SWPC's current forecast for the decline",
              "average of Cycles 10–24 ± 1σ, started at the same minimum"]
    ax.legend(handles, labels, loc="upper left", handlelength=2.6)
    ax.set_xlim(2018.8, 2031.0)
    ax.set_ylim(0, 290)
    ax.set_xticks(np.arange(2019, 2031, 1))
    ax.set_ylabel("sunspot number")
    ax.set_title(f"Cycle 25 peaked {r['max_ratio_to_panel'] - 1:.0%} above its forecast and "
                 f"{r['max_months_early']} months early — yet it was an ordinary cycle")

    # (b) what the early rise said
    hind = pd.DataFrame(r["hindcast"])
    x = hind["t_available"].to_numpy()
    axh.axhline(c25["max_sn"], color=ps.C_SMOOTHED, lw=1.2)
    axh.axhline(LIT_PANEL_2019["max_ssn"], color=ps.C_LITERATURE, lw=1.2)
    axh.vlines(x, hind["lo95"], hind["hi95"], color=ps.INK_SECONDARY, lw=0.9)
    axh.vlines(x, hind["lo68"], hind["hi68"], color=ps.INK_SECONDARY, lw=3.2)
    axh.scatter(x, hind["prediction"], s=34, color=ps.SURFACE, edgecolor=ps.INK, linewidth=1.4, zorder=4)
    axh.set_xticks(x)
    axh.set_xticklabels([f"{_month_name(a)}\n+{lead} months" for a, lead in
                         zip(hind["available"], hind["lead_months"])], fontsize=7.5)
    axh.text(2024.35, c25["max_sn"] + 3, "observed peak", fontsize=8, color=ps.INK_SECONDARY, va="bottom",
             ha="right")
    axh.text(2024.35, LIT_PANEL_2019["max_ssn"] - 3, "panel forecast", fontsize=8, color=ps.INK_SECONDARY,
             va="top", ha="right")
    axh.set_xlim(2021.1, 2024.4)
    axh.set_ylim(40, 300)
    axh.set_xlabel("when the prediction could have been made (and how long after minimum)")
    axh.set_ylabel("predicted cycle maximum")
    fe = r["panel_first_excluded_68"]
    axh.set_title("What the early rise said about the peak")
    axh.text(0.985, 0.03, f"calibrated on Cycles 10–24; bars: 68% and 95% ranges.\nBy {_month_name(fe['available'])} "
             f"the forecast of 115 was outside\nthe 68% range — but the peak came higher still.",
             transform=axh.transAxes, ha="right", va="bottom", fontsize=7.5, color=ps.MUTED)

    # (c) the two hemispheres
    h = tables["hemispheric"]
    ok = h["north_s"].notna()
    axs.plot(h.loc[ok, "t"], h.loc[ok, "north_s"], color=ps.C_NORTH, lw=1.8)
    axs.plot(h.loc[ok, "t"], h.loc[ok, "south_s"], color=ps.C_SOUTH, lw=1.8)
    for row in results["hemispheres"]["silso"]:
        for key, color in (("north", ps.C_NORTH), ("south", ps.C_SOUTH)):
            axs.scatter([_t(row[f"{key}_peak"])], [row[f"{key}_sn"]], s=26, color=color,
                        edgecolor=ps.SURFACE, linewidth=1.2, zorder=4)
    for c in cyc[cyc["t_min"] > 1995].itertuples():
        axs.axvline(c.t_min, color=ps.AXIS, lw=0.7, zorder=0)
        axs.text(c.t_min + 0.3, 3, f"Cycle {c.cycle}", fontsize=7.5, color=ps.MUTED, va="bottom")
    last = h[ok].iloc[-1]
    axs.text(last["t"] + 0.3, last["north_s"] + 2, "north", color=ps.INK_SECONDARY, fontsize=8, va="bottom")
    axs.text(last["t"] + 0.3, last["south_s"] - 2, "south", color=ps.INK_SECONDARY, fontsize=8, va="top")
    rows = {row["cycle"]: row for row in results["hemispheres"]["silso"]}
    axs.set_xlim(1992, 2029.5)
    axs.set_ylim(0, 125)
    axs.set_ylabel("hemispheric sunspot number, smoothed")
    axs.set_title("The two hemispheres peak at different times")
    axs.text(2004.2, 121, f"south peaked {rows[24]['south_minus_north_yr']:.1f} yr after\nnorth in "
             f"Cycle 24, {rows[25]['south_minus_north_yr']:.1f} yr after in\nCycle 25 — the whole-Sun "
             "peak came\nwith the south's", va="top", fontsize=7.5, color=ps.MUTED)

    ps.caption(fig, f"{_credit_sn()} Hemispheric numbers exist since 1992. Forecasts: NOAA/NASA Solar Cycle 25 "
               f"Prediction Panel (2019) and NOAA SWPC. Cycle 25's peak ranks {r['rank']}th of 25; the mean "
               f"peak of Cycles 10–24 is {r['mean_max_calibration_cycles']:.0f}.", y=0.03)
    return ps.save(fig, FIGURE_DIR / "05_cycle25")


# --- 06 Differential rotation -------------------------------------------------------------------------------


def figure_rotation(tables, results) -> list:
    p = tables["rotation_pairs"]
    b = tables["rotation_binned"]
    ro = results["rotation"]
    fit = ro["fit"]
    lit = LIT_BALTHASAR_1986
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.3), gridspec_kw={"width_ratios": [1.45, 1.0], "wspace": 0.22})

    ax = axes[0]
    lat = np.linspace(0, 40, 200)
    ax.plot(lat, sp.rotation_law(lat, lit["A"], lit["B"]), color=ps.C_LITERATURE, lw=1.3, zorder=3)
    ax.plot(lat, sp.rotation_law(lat, fit["A"], fit["B"]), color=ps.C_FOCUS, lw=2.0, zorder=3)
    ax.errorbar(b["lat_mid"], b["mean"], yerr=b["sem"], fmt="o", ms=5, color=ps.C_SMOOTHED,
                mec=ps.SURFACE, mew=1.2, elinewidth=1.0, zorder=4)
    ax.set_xlim(0, 40)
    ax.set_ylim(13.2, 14.85)
    ax.set_xlabel("latitude (either hemisphere)")
    ax.set_ylabel("rotation rate, degrees per day (sidereal)")
    ax.set_title("The Sun is not solid: its equator laps the higher latitudes")
    handles = [
        Line2D([], [], marker="o", ls="", ms=5, color=ps.C_SMOOTHED, mec=ps.SURFACE,
               label=f"mean of {ro['n_pairs']:,} day-to-day moves of {ro['n_groups']:,} groups, per 2.5°"),
        Line2D([], [], color=ps.C_FOCUS, lw=2.0,
               label=f"our fit: Ω = {fit['A']:.3f} − {abs(fit['B']):.2f} sin²(latitude)  "
                     f"(±{fit['A_err_boot']:.3f}, ±{fit['B_err_boot']:.2f})"),
        Line2D([], [], color=ps.C_LITERATURE, lw=1.3,
               label=f"Balthasar et al. (1986), same plates: Ω = {lit['A']} − {abs(lit['B'])} sin²(latitude)"),
    ]
    ax.legend(handles=handles, loc="lower left")
    per = ro["periods"]
    for key, lat0 in (("0deg", 0), ("15deg", 15), ("30deg", 30), ("40deg", 40)):
        om = per[key]["omega"]
        dx, dy, ha = {0: (10, 16, "left"), 15: (0, 16, "center"), 30: (0, 16, "center"),
                      40: (-14, 30, "right")}[lat0]
        ax.annotate(f"{per[key]['sidereal_days']:.1f} days\nper turn", (lat0, om), xytext=(dx, dy),
                    textcoords="offset points", ha=ha, fontsize=7.5, color=ps.INK_SECONDARY)

    ax = axes[1]
    bands = ((0, 5, ps.BLUE), (15, 20, ps.ORANGE), (28, 35, ps.AQUA))
    grid = np.linspace(11.0, 17.5, 400)
    for lo, hi, color in bands:
        v = p.loc[p["lat"].abs().between(lo, hi), "omega"].to_numpy()
        kde = gaussian_kde(v, bw_method=0.12)
        ax.plot(grid, kde(grid), color=color, lw=1.8)
        med = float(np.median(v))
        ax.vlines(med, 0, kde([med])[0], color=color, lw=0.9)
    handles = [Line2D([], [], color=c, lw=1.8,
                      label=f"{lo}–{hi}°:  {p.loc[p['lat'].abs().between(lo, hi), 'omega'].median():.2f}")
               for lo, hi, c in bands]
    ax.legend(handles=handles, loc="upper right", title="latitude: median °/day", title_fontsize=8)
    ax.set_xlim(11.0, 17.5)
    ax.set_ylim(0, None)
    ax.set_xlabel("rotation rate from a single day's move (degrees per day)")
    ax.set_ylabel("share of moves (density)")
    ax.set_title("One day's move is noisy; the shift is not")

    lap = ro["equator_to_40deg_lap_days"]
    ps.caption(fig, f"{_credit_groups()} Royal Greenwich Observatory positions 1874–1976 (to 0.1°), "
               f"groups within {ro['max_cmd_deg']:.0f}° of the central meridian, consecutive days. Errors from "
               f"resampling groups. The equator gains a full lap on 40° latitude every {lap:.0f} days.", y=-0.02)
    return ps.save(fig, FIGURE_DIR / "06_rotation")


# --- 07 Four centuries ------------------------------------------------------------------------------------------


def figure_four_centuries(tables, results) -> list:
    y = tables["yearly"]
    gns = tables["group_numbers"]
    lr = results["long_record"]["series"]
    fig = plt.figure(figsize=(13.6, 8.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1.0], width_ratios=[1.3, 1.0], hspace=0.42, wspace=0.34)
    ax = fig.add_subplot(gs[0, :])
    ax_m = fig.add_subplot(gs[1, 0])
    ax_r = fig.add_subplot(gs[1, 1])

    style = {"HS98": ps.ORANGE, "SS16": ps.AQUA, "U16": ps.YELLOW}
    short = {"HS98": "Hoyt & Schatten", "SS16": "Svalgaard & Schatten\n(backbone)",
             "U16": "Usoskin et al.\n(active days)"}
    for (lo, hi), name in ((MAUNDER_MINIMUM, "Maunder minimum"), (DALTON_MINIMUM, "Dalton minimum")):
        ax.axvspan(lo, hi, color=ps.WASH, lw=0, zorder=0)
        ax.text((lo + hi) / 2, 292, name, ha="center", va="top", fontsize=8, color=ps.INK_SECONDARY)
    for key, color in style.items():
        g = gns[key]
        ax.plot(g["t"], g["sn_equiv"], color=color, lw=1.3, zorder=3)
    ax.plot(y["t"], y["sn"], color=ps.C_SMOOTHED, lw=1.8, zorder=4)
    handles = [Line2D([], [], color=ps.C_SMOOTHED, lw=1.8, label="sunspot number, yearly (SILSO v2.0)")]
    for key, color in style.items():
        s = lr[key]
        handles.append(Line2D([], [], color=color, lw=1.5,
                              label=f"{s['label']}: group number × {s['scale_to_sn_v2']:.1f}"))
    ax.legend(handles=handles, loc="lower left", ncol=2, bbox_to_anchor=(0.0, 1.0))
    ax.set_xlim(1608, 2028)
    ax.set_ylim(0, 300)
    ax.set_xticks(np.arange(1625, 2026, 25))
    ax.set_ylabel("sunspot number (or its equivalent)")
    ax.set_title("Four centuries of sunspots: the telescope record back to 1610", pad=40)

    g = gns["HS98"].set_index("t")
    s2 = gns["SS16"].set_index("t")
    sel = (g.index >= 1640) & (g.index <= 1720)
    yrs = g.index[sel]
    ax_m.axvspan(MAUNDER_MINIMUM[0], MAUNDER_MINIMUM[1] + 1, color=ps.WASH, lw=0, zorder=0)
    ax_m.bar(yrs - 0.2, s2.reindex(yrs)["sn_equiv"], width=0.4, color=ps.AQUA, lw=0,
             label="Svalgaard & Schatten (backbone)")
    ax_m.bar(yrs + 0.2, g.loc[sel, "sn_equiv"], width=0.4, color=ps.ORANGE, lw=0, label="Hoyt & Schatten")
    ax_m.set_xlim(1639, 1721)
    ax_m.set_ylim(0, 150)
    ax_m.set_ylabel("sunspot number equivalent")
    ax_m.legend(loc="upper left")
    hs, ss = lr["HS98"], lr["SS16"]
    ax_m.set_title("Seventy years of near silence, 1645–1715")
    ax_m.text(1646.5, 118, f"Hoyt & Schatten: no spot seen in {hs['maunder_years_with_no_spot_seen']}\n"
              f"of {hs['maunder_years']} years, though those years have a\nmedian of "
              f"{hs['maunder_median_obs_in_blank_years']:.0f} recorded observations each.\n"
              f"Mean level {hs['maunder_mean_sn_equiv']:.0f} (Hoyt & Schatten) or "
              f"{ss['maunder_mean_sn_equiv']:.0f} (backbone);\n"
              f"since 1749 the average is {results['long_record']['modern_mean_sn_1749_2025']:.0f}.",
              va="top", fontsize=7.5, color=ps.MUTED)

    keys = [k for k in ("HS98", "U16", "SS16") if "early_ratio_1749_1849" in lr[k]]
    vals = [lr[k]["early_ratio_1749_1849"] for k in keys]
    ypos = np.arange(len(keys))
    ax_r.barh(ypos, vals, height=0.5, color=[style[k] for k in keys], lw=0)
    ax_r.axvline(1.0, color=ps.C_SMOOTHED, lw=1.4)
    for yv, v in zip(ypos, vals):
        ax_r.text(v + 0.02, yv, f"{v:.2f}", va="center", fontsize=8, color=ps.INK)
    ax_r.set_yticks(ypos)
    ax_r.set_yticklabels([short[k] for k in keys], fontsize=8)
    ax_r.set_xlim(0, 1.35)
    ax_r.set_xlabel("mean activity 1749–1849, as a multiple of the sunspot number's")
    ax_r.grid(axis="y", visible=False)
    ax_r.set_title("How strong were the early cycles?\nIt depends on the reconstruction")

    ps.caption(fig, f"{_credit_sn()} Group numbers from the SILSO archive: Hoyt & Schatten (1998), "
               "Svalgaard & Schatten (2016), Usoskin et al. (2016; their own method 1749–1899). Each is scaled "
               "to the sunspot number over 1874–1976. Yearly sunspot numbers before 1749 are Wolf's reconstruction.",
               y=0.035)
    return ps.save(fig, FIGURE_DIR / "07_four_centuries")


# --- 08 Caveats ----------------------------------------------------------------------------------------------------


def figure_caveats(tables, results) -> list:
    d = tables["daily"]
    per = tables["area_per_decade"]
    seam = results["butterfly"]["seam"]
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.0), gridspec_kw={"wspace": 0.22})

    ax = axes[0]
    cov = d.groupby("year")["sn"].apply(lambda v: v.notna().mean())
    cov = cov[cov.index < 1870]
    ax.bar(cov.index, cov.values * 100, width=0.7, color=ps.BLUE, lw=0)
    ax.set_xlim(1816.5, 1869.5)
    ax.set_ylim(0, 125)
    ax.set_yticks(np.arange(0, 101, 20))
    ax.set_ylabel("days with an observation (%)")
    ax.set_xlabel("year")
    ax.set_title("Before 1849, many days went unobserved")
    ax.text(0.01, 0.98, "No daily record at all before 1818; the monthly values back to 1749\n"
            "rest on scattered observers and interpolation.", transform=ax.transAxes, va="top",
            fontsize=7.5, color=ps.MUTED)

    ax = axes[1]
    color = {"RGO": ps.BLUE, "USAF": ps.ORANGE, "NOAA": ps.ORANGE}
    mid = (per["first_year"] + per["last_year"] + 1) / 2.0
    width = (per["last_year"] - per["first_year"] + 1) * 0.8
    ax.bar(mid, per["relative_area_per_sn"], width=width, color=[color[s] for s in per["source"]], lw=0)
    ax.axhline(1.0, color=ps.INK_SECONDARY, lw=0.8)
    ax.axhline(1 / seam["factor"], color=ps.ORANGE, lw=0.8)
    ax.text(2031, 1.0 + 0.015, "RGO\naverage", fontsize=7.5, color=ps.INK_SECONDARY, va="bottom", ha="right")
    usaf = per[per["source"] == "USAF"].iloc[0]
    ax.text((usaf["first_year"] + usaf["last_year"] + 1) / 2, usaf["relative_area_per_sn"] / 2, "1977–81",
            rotation=90, ha="center", va="center", fontsize=7, color="white")
    ax.set_xlim(1868, 2032)
    ax.set_ylim(0, 1.6)
    ax.set_yticks(np.arange(0, 1.21, 0.2))
    ax.set_ylabel("spot area per sunspot number\n(relative to RGO, at equal sunspot number)")
    ax.set_xlabel("year")
    ax.set_title("A seam in 1977: the same Sun, smaller spot areas")
    handles = [Patch(color=ps.BLUE, label="Royal Greenwich Observatory"),
               Patch(color=ps.ORANGE, label="US Air Force / NOAA"),
               Line2D([], [], color=ps.ORANGE, lw=0.8, label=f"USAF/NOAA level (÷{seam['factor']:.2f})")]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 0.84), ncol=2)
    lit = seam["literature"]
    ax.text(0.01, 0.98, f"Across {seam['window'][0]}–{seam['window'][1]}, USAF/NOAA areas are "
            f"{seam['factor']:.2f}× smaller\n(68%: {seam['factor_ci68'][0]:.2f}–{seam['factor_ci68'][1]:.2f}); "
            f"{seam['factor_usaf_1977_1981']:.2f}× for 1977–81 alone. Published:\n{lit['hathaway']} (Hathaway), "
            f"{lit['balmaceda_corrected']} (Balmaceda et al. 2009).", transform=ax.transAxes, va="top",
            fontsize=7.5, color=ps.MUTED)

    ps.caption(fig, f"{_credit_sn()} {_credit_groups()} Area per sunspot number is fitted to monthly means "
               "(months with sunspot number above 20) with one common slope; each bar is one period's level "
               "relative to the RGO mean.", y=-0.03)
    return ps.save(fig, FIGURE_DIR / "08_caveats")


def make_all(tables, results) -> list:
    written = []
    for f in (figure_record, figure_butterfly, figure_zone_drift, figure_rules, figure_cycle25,
              figure_rotation, figure_four_centuries, figure_caveats):
        written += f(tables, results)
    return written
