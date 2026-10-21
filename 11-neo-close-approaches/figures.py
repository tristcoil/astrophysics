"""Figure generation.

Every figure is drawn from the analysis outputs; nothing here re-derives a
result, so what is plotted is exactly what was measured.

Five figures, in reading order:

    01  this month's close approaches -- what passed, what is coming, and why they differ
    02  who comes close              -- size vs miss distance, passes inside the Moon's
                                        orbit, and when each asteroid was noticed
    03  from near misses to impacts  -- the d^2 law, and how often Earth is hit
    04  the risk list                -- NASA's list on the Torino scale; how it churns
    05  Apophis, 13 April 2029       -- the flyby to scale, the new orbit, why we know

Design notes, as in A1/A2: log axes carry plain-language tick labels rather than
powers of ten; text is always ink, never the colour of a data series; every
figure carries a caption with its provenance and the caveat a reader needs.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, NullLocator

import neo_physics as nphys
import plotstyle as ps
from config import (
    DENSITY_KG_M3, EARTH_RADIUS_KM, FIGURE_DIR, GEO_RADIUS_KM, GPS_RADIUS_KM, LD_KM,
    SIZE_CLASS_LABELS,
)

ps.use_style()

GEO_LD = GEO_RADIUS_KM / LD_KM
EARTH_LD = EARTH_RADIUS_KM / LD_KM
#: Text over data gets a surface-coloured backing so it stays legible.
BACKING = {"facecolor": ps.SURFACE, "edgecolor": "none", "pad": 1.0, "alpha": 0.88}


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _provenance(results) -> str:
    return ("Data: NASA/JPL Solar System Dynamics and CNEOS — close-approach, small-body, Sentry, "
            f"fireball and Horizons services, retrieved {results['snapshot']['now_utc'][:10]}.")


def _num(x, _pos=None) -> str:
    """Plain-number tick labels for log axes: 0.01, 0.1, 1.5, 10, 1,000."""
    if x >= 1 and float(x).is_integer():
        return f"{x:,.0f}"
    return f"{x:.10f}".rstrip("0").rstrip(".")


def _plain_log(ax, which="both"):
    for axis in ([ax.xaxis, ax.yaxis] if which == "both" else [getattr(ax, f"{which}axis")]):
        axis.set_major_formatter(FuncFormatter(_num))
        axis.set_minor_formatter(NullFormatter())


def _metres(x, _pos=None) -> str:
    return f"{x / 1000:g} km" if x >= 1000 else f"{x:g} m"


def _sig(x: float, digits: int = 2) -> str:
    """Round to significant figures for display: 33,136 -> '33,000'."""
    if x <= 0:
        return "0"
    rounded = round(x, -int(np.floor(np.log10(abs(x)))) + digits - 1)
    return f"{rounded:,.0f}" if rounded >= 1 else f"{rounded:g}"


def _save(fig, stem: str) -> list:
    return ps.save(fig, FIGURE_DIR / stem)


def _size_points(d_m):
    """Marker area (pt^2) growing with log diameter, 1 m -> 6, 1 km -> 256.
    Always shown with a size key, never as the only way to read a size."""
    x = np.clip(np.log10(np.asarray(d_m, float)) / 3.0, 0.0, 1.1)
    return 6.0 + 250.0 * x**2


def _reference_line(ax, y, label, *, x_text, ha="right", va="bottom"):
    ax.axhline(y, color=ps.AXIS, linewidth=0.8, zorder=1)
    ax.text(x_text, y, label, ha=ha, va=va, fontsize=7.5, color=ps.INK_SECONDARY, zorder=6,
            bbox=BACKING)


def _lag_category(lag_days) -> np.ndarray:
    lag = np.asarray(lag_days, float)
    return np.select([lag < 0, lag == 0], ["before", "same_day"], "after")


def _fmt_date(ts, year=True) -> str:
    ts = pd.Timestamp(ts)
    return f"{ts.day} {ts.strftime('%b')}" + (f" {ts.year}" if year else "")


def _short_name(fullname: str) -> str:
    """'101955 Bennu (1999 RQ36)' -> 'Bennu'; '29075 (1950 DA)' -> '1950 DA';
    '(2010 RF12)' -> '2010 RF12'."""
    name = fullname.strip()
    if "(" in name:
        before, inside = name.split("(", 1)
        words = before.split()
        proper = [w for w in words if not w.isdigit()]
        return " ".join(proper) if proper else inside.rstrip(")")
    return name


def _text_angle(ax, x0, y0, x1, y1) -> float:
    """On-screen angle of the segment (x0, y0)-(x1, y1), for text along a line."""
    (a, b), (c, d) = ax.transData.transform([(x0, y0), (x1, y1)])
    return float(np.degrees(np.arctan2(d - b, c - a)))


# ===========================================================================
# 01  This month's close approaches
# ===========================================================================

def figure_snapshot(results: dict, art: dict) -> list:
    """Timeline of every approach within 0.05 au from a month before the data
    were retrieved to two months after, with daily counts above."""
    s = results["snapshot"]
    w = art["frame"].copy()
    now, lo, hi = art["now"], art["lo"], art["hi"]
    w["cat"] = _lag_category(w["lag_days"])
    past, fut = s["past"], s["future"]

    fig = plt.figure(figsize=(11.5, 7.8))
    gs = GridSpec(2, 1, height_ratios=[1, 3.3], hspace=0.07, top=0.88, bottom=0.08, figure=fig)
    axc = fig.add_subplot(gs[0])
    ax = fig.add_subplot(gs[1], sharex=axc)

    # --- daily counts ------------------------------------------------------
    days = pd.date_range(lo.normalize(), hi.normalize(), freq="D")
    bottom = np.zeros(len(days))
    for cat in ("before", "same_day", "after"):
        daily = (w[w["cat"] == cat].groupby(w["time"].dt.normalize()).size()
                 .reindex(days, fill_value=0).to_numpy())
        axc.bar(days + pd.Timedelta(hours=12), daily, width=0.78, bottom=bottom,
                color=ps.LAG_COLORS[cat], edgecolor=ps.SURFACE, linewidth=0.5, zorder=3)
        bottom += daily
    top = bottom.max() * 1.55
    axc.set_ylim(0, top)
    axc.set_ylabel("Approaches\nper day")
    axc.tick_params(labelbottom=False)
    axc.yaxis.set_major_locator(FixedLocator([0, 5, 10]))
    axc.text(lo + (now - lo) / 2, top * 0.97, f"{past['per_day']:.1f} a day", ha="center",
             va="top", fontsize=8, color=ps.INK)
    axc.text(now + (hi - now) / 2, top * 0.97, f"{fut['per_day']:.1f} a day predicted",
             ha="center", va="top", fontsize=8, color=ps.INK)
    axc.text(now + pd.Timedelta(hours=14), top * 0.6, f"data retrieved\n{_fmt_date(now)}",
             ha="left", va="center", fontsize=7.5, color=ps.INK_SECONDARY)

    # --- the timeline ------------------------------------------------------
    for cat in ("before", "same_day", "after"):
        sub = w[w["cat"] == cat]
        ax.scatter(sub["time"], sub["dist_ld"], s=_size_points(sub["d_best_m"]),
                   color=ps.LAG_COLORS[cat], edgecolor=ps.SURFACE, linewidth=0.9,
                   alpha=0.92, zorder=4)
    ax.set_yscale("log")
    ax.set_ylim(0.012, 40)
    ax.yaxis.set_major_locator(FixedLocator([0.01, 0.1, 1, 10]))
    _plain_log(ax, "y")
    ax.set_ylabel("Miss distance from Earth's centre (lunar distances, LD)")
    x_text = hi - pd.Timedelta(days=0.5)
    _reference_line(ax, 1.0, "the Moon's average distance (384,400 km)", x_text=x_text)
    _reference_line(ax, GEO_LD, "geostationary satellites (42,164 km)", x_text=x_text)
    _reference_line(ax, EARTH_LD, "Earth's surface", x_text=x_text)
    ax.set_xlim(lo, hi)
    ax.xaxis.set_major_locator(mdates.DayLocator(bymonthday=[1, 15]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: _fmt_date(mdates.num2date(v), False)))

    for a in (axc, ax):
        a.axvline(now, color=ps.INK_SECONDARY, linewidth=1.0, zorder=5)
    ax.text(now - pd.Timedelta(hours=14), 32, "← already passed", ha="right", va="top",
            fontsize=8, color=ps.INK_SECONDARY)
    ax.text(now + pd.Timedelta(hours=14), 32, "predicted →", ha="left", va="top",
            fontsize=8, color=ps.INK_SECONDARY)

    # --- named objects -----------------------------------------------------
    c = s["closest"]
    t = pd.Timestamp(c["time"])
    ps.annotate(ax, f"{_short_name(c['name'])}, ~{c['diameter_m']:.0f} m: passed {c['dist_km']:,.0f} km\n"
                    "from Earth's centre, half the distance\nof geostationary satellites",
                (t, c["dist_ld"]), (lo + pd.Timedelta(days=1), 0.03), relpos=(0.35, 1),
                ha="left", va="top", fontsize=7.5, color=ps.INK)
    b = s["largest"]
    t = pd.Timestamp(b["time"])
    ps.annotate(ax, f"{_short_name(b['name'])}, ~{b['diameter_m']:.0f} m, the largest:\n"
                    f"{b['dist_ld']:.1f} LD on {_fmt_date(t, False)}; known since "
                    f"{b['first_observed'][:4]}",
                (t, b["dist_ld"]), (t, 4.4), relpos=(0.5, 1),
                ha="center", va="top", fontsize=7.5, color=ps.INK)
    cf = s["closest_future"]
    if cf:
        t = pd.Timestamp(cf["time"])
        ps.annotate(ax, f"{_short_name(cf['name'])}, ~{cf['diameter_m']:.0f} m: the closest predicted\n"
                        f"pass, {cf['dist_ld']:.2f} LD on {_fmt_date(t, False)}; found "
                        f"{cf['first_observed'][:4]}",
                    (t, cf["dist_ld"]), (t + pd.Timedelta(days=5), 0.36), relpos=(0, 0.5),
                    ha="left", va="center", fontsize=7.5, color=ps.INK)

    # --- legends: colour (identity) and size ----------------------------------
    handles = [Line2D([], [], marker="o", linestyle="", markersize=7, color=ps.LAG_COLORS[k],
                      markeredgecolor=ps.SURFACE, label=ps.LAG_LABELS[k])
               for k in ("before", "same_day", "after")]
    leg = ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(1.0, 0.075), fontsize=7.5,
                    borderpad=0.6)
    ax.add_artist(leg)
    sizes = [3, 30, 300]
    size_handles = [Line2D([], [], marker="o", linestyle="", color=ps.MUTED,
                           markeredgecolor=ps.SURFACE, markersize=np.sqrt(_size_points(d)),
                           label=_metres(d)) for d in sizes]
    ax.legend(handles=size_handles, loc="lower right", bbox_to_anchor=(0.66, 0.075), ncol=3,
              title="estimated size", title_fontsize=7.5, fontsize=7.5, columnspacing=1.2,
              handletextpad=0.5, borderpad=0.6)

    fo = past["first_observed"]
    axc.set_title(f"{past['count']} known asteroids passed Earth in the last {s['days_before']} days; "
                  f"{fut['count']} are predicted for the next {s['days_after']}", pad=26)
    axc.text(0, 1.07, f"The future only looks quiet: {past['first_observed_in_window']} of the "
                      f"{past['count']} were unknown a month earlier, and {fo['after']} were first "
                      f"seen only after they had passed.",
             transform=axc.transAxes, fontsize=9, color=ps.INK_SECONDARY, va="bottom")
    ps.caption(fig, f"{_provenance(results)} Every approach of a known near-Earth asteroid within "
                    f"0.05 au (19.5 LD). Sizes come from brightness assuming a typical albedo of 0.14, "
                    f"so true sizes may be 0.75–1.7× those shown, unless measured. 'First seen' is the "
                    f"earliest observation in JPL's orbit fit, which can be a pre-discovery image found "
                    f"later, so 'only afterwards' is if anything an undercount.")
    return _save(fig, "01_this_month")


# ===========================================================================
# 02  Who comes close, and when we notice
# ===========================================================================

def _rate_label(rate: float) -> str:
    if rate >= 1.5:
        return f"{rate:.0f} a year"
    n = round(rate * 10)
    return "none in ten years" if n == 0 else f"{n} in ten years"


def figure_encounters(results: dict, art: dict) -> list:
    e = results["encounters"]
    survey, by_year = art["survey"], art["by_year"]
    inside = survey[survey["dist_ld"] < 1]
    y0, y1 = e["survey_years"]

    fig = plt.figure(figsize=(12, 10))
    gs = GridSpec(2, 2, width_ratios=[1.4, 1], height_ratios=[1.3, 1], hspace=0.38, wspace=0.3,
                  top=0.9, bottom=0.07, figure=fig)
    ax = fig.add_subplot(gs[0, 0])
    axr = fig.add_subplot(gs[0, 1])
    axh = fig.add_subplot(gs[1, :])

    # --- (a) size vs miss distance ------------------------------------------
    ax.scatter(survey["d_best_m"], survey["dist_ld"], s=3, color=ps.C_DATA, alpha=0.3,
               linewidth=0, zorder=3, rasterized=True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.3, 8000)
    ax.set_ylim(0.012, 40)
    ax.xaxis.set_major_locator(FixedLocator([1, 10, 100, 1000]))
    ax.xaxis.set_major_formatter(FuncFormatter(_metres))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_locator(FixedLocator([0.01, 0.1, 1, 10]))
    _plain_log(ax, "y")
    ax.set_xlabel("Estimated diameter")
    ax.set_ylabel("Miss distance (lunar distances)")
    _reference_line(ax, 1.0, "the Moon", x_text=7000)
    _reference_line(ax, GEO_LD, "geostationary satellites", x_text=7000)
    _reference_line(ax, EARTH_LD, "Earth's surface", x_text=7000)
    ax.text(0.4, 11, "Too faint to see:\nsmall asteroids are found\nonly when they pass close",
            fontsize=8, color=ps.INK_SECONDARY, va="top", bbox=BACKING)
    ax.text(5000, 0.035, "Large asteroids\nrarely come this close", fontsize=8,
            color=ps.INK_SECONDARY, ha="right")

    big = e["largest_inside_1ld"]
    ps.annotate(ax, f"{_short_name(big['name'])}, ~{big['diameter_m']:.0f} m ({big['time'][:4]}):\n"
                    "the largest inside the Moon's\ndistance this decade",
                (big["diameter_m"], big["dist_ld"]), (330, 0.33), relpos=(0, 0.5), fontsize=7.5,
                color=ps.INK, ha="left", va="center")
    florence = survey.loc[survey["d_best_m"].idxmax()]
    ps.annotate(ax, f"{_short_name(florence['fullname'])}, {florence['d_best_m'] / 1000:.1f} km "
                    f"({florence['year']})",
                (florence["d_best_m"], florence["dist_ld"]), (1500, 30), relpos=(0.5, 0),
                fontsize=7.5, color=ps.INK, ha="center")
    c0 = e["closest_known_passes"][0]
    ps.annotate(ax, f"{_short_name(c0['name'])}, ~{c0['diameter_m']:.1f} m ({c0['time'][:4]}): "
                    f"{c0['altitude_km']:.0f} km\nabove the surface, the closest known miss",
                (c0["diameter_m"], c0["dist_ld"]), (2.2, 0.021), relpos=(0, 0.5), fontsize=7.5,
                color=ps.INK, ha="left", va="center")
    ax.set_title(f"(a) Every known approach within 0.05 au, {y0}–{y1}  ({len(survey):,})",
                 loc="left")

    # --- (b) passes inside the Moon's distance, by size ---------------------
    labels = list(reversed(SIZE_CLASS_LABELS))  # small at the bottom
    rows = np.arange(len(labels))
    for yv, label in zip(rows, labels):
        rate = e["per_class"][label]["inside_1ld_per_year"]
        members = inside[inside["size_class"] == label]
        text = _rate_label(rate)
        if len(members) == 1:
            text += f" ({_short_name(members['fullname'].iloc[0])})"
        if rate > 0:
            axr.plot([0.03, rate], [yv, yv], color=ps.GRID, linewidth=1.0, zorder=2)
            axr.scatter([rate], [yv], s=64, color=ps.C_DATA, edgecolor=ps.SURFACE,
                        linewidth=1.5, zorder=3)
            axr.text(rate * 1.35, yv, text, va="center", fontsize=8, color=ps.INK)
        else:
            axr.text(0.035, yv, text, va="center", fontsize=8, color=ps.INK_SECONDARY)
    axr.set_xscale("log")
    axr.set_xlim(0.03, 900)
    axr.xaxis.set_major_locator(FixedLocator([0.1, 1, 10, 100]))
    _plain_log(axr, "x")
    axr.set_yticks(rows, labels)
    axr.set_ylim(-0.6, len(labels) - 0.4)
    axr.grid(axis="y", visible=False)
    axr.set_xlabel("Passes per year")
    axr.set_title(f"(b) Passes inside the Moon's distance,\n{y0}–{y1}, by size", loc="left")

    # --- (c) discovery timing, year by year ---------------------------------
    years = by_year.index.to_numpy()
    bottom = np.zeros(len(years))
    for cat in ("before", "same_day", "after"):
        vals = by_year[cat].to_numpy()
        axh.bar(years, vals, bottom=bottom, width=0.72, color=ps.LAG_COLORS[cat],
                edgecolor=ps.SURFACE, linewidth=0.7, label=ps.LAG_LABELS[cat], zorder=3)
        bottom += vals
    axh.set_xlim(years[0] - 0.7, years[-1] + 0.7)
    axh.set_ylim(0, bottom.max() * 1.12)
    axh.set_ylabel("Known passes inside\nthe Moon's distance")
    axh.legend(loc="upper left", fontsize=7.5)
    ry0, ry1 = e["recent_years"]
    axh.text(years[-1], bottom[-1] + bottom.max() * 0.015, f"{int(bottom[-1])}", ha="center",
             va="bottom", fontsize=8, color=ps.INK)
    axh.text(0.34, 0.95,
             f"{ry0}–{ry1}: {e['recent_fraction_after']:.0%} of these asteroids were first seen only "
             f"after their closest\napproach, {e['recent_fraction_same_day']:.0%} on the day itself. "
             f"The rise is our telescopes improving,\nnot the sky getting busier.",
             transform=axh.transAxes, ha="left", va="top", fontsize=8.5, color=ps.INK)
    axh.set_title("(c) Year by year: when was each asteroid that passed inside the Moon's distance "
                  "first seen?", loc="left")

    fig.suptitle("Small asteroids pass inside the Moon's orbit every few days; large ones almost never",
                 x=0.07, ha="left", y=0.975)
    ps.caption(fig, f"{_provenance(results)} Sizes from brightness assuming albedo 0.14 unless measured "
                    f"(true sizes may be 0.75–1.7× those shown). Distances are from Earth's centre; "
                    f"1 LD = 384,400 km. The close-approach table leaves out impacts: the few tiny "
                    f"asteroids discovered hours before hitting Earth, such as 2008 TC3, are not shown.")
    return _save(fig, "02_who_comes_close")


# ===========================================================================
# 03  From near misses to impacts
# ===========================================================================

def _interval_label(rate, _pos=None) -> str:
    if rate >= 1:
        return f"{rate:,.0f} per year"
    years = 1 / rate
    names = {1e2: "1 per century", 1e6: "1 per million years", 1e8: "1 per 100 million years"}
    for value, name in names.items():
        if abs(np.log10(years) - np.log10(value)) < 0.01:
            return name
    return f"1 per {_sig(years)} years"


def figure_impacts(results: dict, art: dict) -> list:
    ir = results["impact_rates"]
    fig = plt.figure(figsize=(12, 11))
    gs = GridSpec(2, 2, height_ratios=[1, 1.45], hspace=0.36, wspace=0.28, top=0.9,
                  bottom=0.07, figure=fig)
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])
    ax = fig.add_subplot(gs[1, :])

    # --- (a) share of approaches within d, three size classes ---------------
    d = art["d_grid_ld"]
    shown = [("< 10 m", ps.SIZE_RAMP_3[0]), ("20 – 50 m", ps.SIZE_RAMP_3[1]),
             ("≥ 1 km", ps.SIZE_RAMP_3[2])]
    axa.set_xscale("log")
    axa.set_yscale("log")
    axa.set_xlim(0.004, 25)
    axa.set_ylim(1e-6, 1.6)
    axa.plot(d, (d / d[-1]) ** 2, color=ps.INK_SECONDARY, linewidth=1.0, zorder=2)
    angle = _text_angle(axa, 0.05, (0.05 / d[-1]) ** 2, 0.5, (0.5 / d[-1]) ** 2)
    axa.text(0.06, (0.06 / d[-1]) ** 2 * 0.45, "d²: geometry alone", rotation=angle,
             rotation_mode="anchor", fontsize=7.5, color=ps.INK_SECONDARY, ha="left", va="top")
    label_at = {"< 10 m": None, "20 – 50 m": None, "≥ 1 km": (2.2, 1.1e-3)}
    for label, color in shown:
        curve = art["curves"][label]
        frac = curve / curve[-1]
        ok = frac > 0
        axa.plot(d[ok], frac[ok], color=color, linewidth=2.0, zorder=3)
        k = ir["distance_index"][label]["k"]
        i0 = np.argmax(ok)
        if label_at[label] is None:
            axa.text(d[i0] * 0.85, frac[i0], f"{label}  (k = {k:.2f})", ha="right", va="center",
                     fontsize=7.5, color=ps.INK)
        else:
            axa.text(*label_at[label], f"{label}  (k = {k:.2f})", ha="left", va="center",
                     fontsize=7.5, color=ps.INK)
    axa.xaxis.set_major_locator(FixedLocator([0.01, 0.1, 1, 10]))
    axa.yaxis.set_major_locator(FixedLocator([1e-6, 1e-4, 1e-2, 1]))
    axa.yaxis.set_major_formatter(FuncFormatter(lambda v, _: {1e-6: "1 in a million",
                                                              1e-4: "1 in 10,000",
                                                              1e-2: "1 in 100", 1: "all"}.get(v, "")))
    axa.xaxis.set_major_formatter(FuncFormatter(_num))
    axa.xaxis.set_minor_formatter(NullFormatter())
    axa.yaxis.set_minor_formatter(NullFormatter())
    axa.set_xlabel("Miss distance d (lunar distances)")
    axa.set_ylabel("Share of approaches within 0.05 au\nthat pass closer than d")
    axa.set_title("(a) Big asteroids follow the d² law;\nsmall ones don't", loc="left")

    # --- (b) the distance index, every size class ----------------------------
    labels = SIZE_CLASS_LABELS
    x = np.arange(len(labels))
    ks = np.array([ir["distance_index"][lab]["k"] for lab in labels])
    lo = np.array([ir["distance_index"][lab]["k_16_50_84"][0] for lab in labels])
    hi = np.array([ir["distance_index"][lab]["k_16_50_84"][2] for lab in labels])
    axb.axhline(2.0, color=ps.INK_SECONDARY, linewidth=1.0, zorder=2)
    axb.text(len(labels) - 0.5, 2.05, "k = 2: geometry alone (Öpik)", ha="right", va="bottom",
             fontsize=7.5, color=ps.INK_SECONDARY)
    axb.errorbar(x, ks, yerr=[ks - lo, hi - ks], fmt="none", ecolor=ps.C_DATA, elinewidth=1.4,
                 capsize=0, zorder=3)
    axb.scatter(x, ks, s=56, color=ps.C_DATA, edgecolor=ps.SURFACE, linewidth=1.5, zorder=4)
    for xi, k in zip(x, ks):
        axb.text(xi + 0.14, k, f"{k:.2f}", fontsize=7.5, color=ps.INK, va="center")
    axb.set_xticks(x, [lab.replace(" – ", "–") for lab in labels], fontsize=7.5)
    axb.set_xlim(-0.5, len(labels) - 0.4)
    axb.set_ylim(0, 2.6)
    axb.grid(axis="x", visible=False)
    axb.set_ylabel("Distance index k,  N(< d) ∝ dᵏ")
    axb.set_xlabel("Asteroid size")
    axb.text(3.6, 0.35, "k below 2: distant passes\nof small asteroids go unseen",
             fontsize=7.5, color=ps.INK_SECONDARY, ha="center")
    axb.set_title("(b) The index falls as asteroids get smaller", loc="left")

    # --- (c) how often Earth is hit -----------------------------------------
    v_typ = art["v_typical"]
    d_m = art["d_grid_m"]

    def energy_to_d(e_mt):
        return nphys.diameter_for_energy_km(e_mt, v_typ, DENSITY_KG_M3) * 1e3

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.6, 20000)
    ax.set_ylim(1e-9, 200)

    # published curves (violet)
    e_pal = np.logspace(-4, 8.3, 200)
    ax.plot(energy_to_d(e_pal), nphys.background_impact_rate(e_pal), color=ps.C_PUBLISHED,
            linewidth=1.2, zorder=2)
    ax.text(19000, 1.2e-6, "background rate behind\nthe Palermo scale\n(Chesley et al. 2002)",
            fontsize=7.5, color=ps.INK, ha="right", va="bottom", bbox=BACKING, zorder=7)
    e_brown = np.logspace(-1, 1, 50) / 1e3  # 0.1-10 kt, the range Brown et al. fitted
    ax.plot(energy_to_d(e_brown), nphys.brown2002_rate(e_brown * 1e3), color=ps.C_PUBLISHED,
            linewidth=1.8, zorder=5, solid_capstyle="round")
    ax.text(2.4, 11, "fit to satellite-detected\nbolides, Brown et al. 2002", fontsize=7.5,
            color=ps.INK, ha="left", va="bottom")

    # the fireball record (blue)
    fb = art["fb_curve"]
    ok = fb > 0
    ax.step(d_m[ok], fb[ok], where="post", color=ps.C_DATA, linewidth=2.0, zorder=4)
    f = ir["fireballs"]
    n_fb_years = f["years"][1] - f["years"][0] + 1
    ax.scatter([f["largest"]["diameter_m"]], [1 / n_fb_years], s=50, color=ps.C_DATA,
               edgecolor=ps.SURFACE, linewidth=1.5, zorder=5)
    ps.annotate(ax, f"Chelyabinsk, 2013: ~{f['largest']['diameter_m']:.0f} m, "
                    f"{f['largest']['impact_kt']:.0f} kt\n(the largest fireball on record)",
                (f["largest"]["diameter_m"], 1 / n_fb_years), (60, 0.25), relpos=(0, 0.5),
                fontsize=7.5, color=ps.INK, ha="left", va="center")
    ax.text(0.68, 60, f"fireballs actually seen,\n{f['years'][0]}–{f['years'][1]} "
                      f"({f['events']:,} events)", fontsize=7.5, color=ps.INK, ha="left",
            va="bottom")

    # our estimates (orange): small asteroids, a lower limit ...
    sm, count = art["small_curve"], art["small_count"]
    ok = (sm > 0) & (count >= 5)  # stop where fewer than five passes remain
    ax.step(d_m[ok], sm[ok], where="post", color=ps.C_FIT, linewidth=2.0, zorder=4)
    sl = ir["small_lower_limits"]
    ax.text(0.68, 0.045, f"known asteroids passing inside the\nMoon's distance, "
                         f"{sl['years'][0]}–{sl['years'][1]}, scaled\ndown to Earth's size "
                         f"(a lower limit)", fontsize=7.5, color=ps.INK, ha="left", va="top")
    sv = ir["seen_vs_fireballs"]
    ax.text(0.68, 3.2, f"the gap: small asteroids\nthat pass unseen. Surveys\ncatch "
                       f"~{sv['>= 10 m']['fraction_seen']:.0%} of 10 m rocks,\n"
                       f"~{sv['>= 3 m']['fraction_seen']:.0%} of 3 m ones",
            fontsize=7.5, color=ps.INK_SECONDARY, ha="left", va="top")

    # ... and large asteroids, from centuries of computed approaches
    large = ir["large"]
    keys = ("140 m", "340 m", "1 km")
    xs = np.array([large[k]["min_diameter_m"] for k in keys])
    ys = np.array([large[k]["rate_per_year"] for k in keys])
    lo = np.array([large[k]["rate_16_50_84"][0] for k in keys])
    hi = np.array([large[k]["rate_16_50_84"][2] for k in keys])
    ax.errorbar(xs, ys, yerr=[ys - lo, hi - ys], fmt="none", ecolor=ps.C_FIT, elinewidth=1.5,
                zorder=5)
    ax.scatter(xs, ys, s=60, color=ps.C_FIT, edgecolor=ps.SURFACE, linewidth=1.5, zorder=6)
    placements = {"140 m": ((0.84, 1.0), "right", "center"),
                  "340 m": ((1.0, 0.45), "center", "top"),
                  "1 km": ((0.84, 1.0), "right", "center")}
    for xv, yv, key in zip(xs, ys, keys):
        (fx, fy), ha, va = placements[key]
        ax.text(xv * fx, yv * fy, "1 per\n" + _interval_label(yv)[6:], fontsize=7.5, color=ps.INK,
                ha=ha, va=va)
    ax.text(125, 3.2e-4, "every approach within 0.05 au by a known\nasteroid, 1900–2200, "
                         "scaled down to Earth's size", fontsize=7.5, color=ps.INK, ha="left",
            va="bottom", bbox=BACKING, zorder=7)

    # published points at 140 m and 1 km (violet), offset slightly to the right
    sb = large["1 km"]["published"]
    rate_sb = 1 / sb["interval_years"]
    rate_sb_lo = 1 / (sb["interval_years"] + sb["interval_err_years"])  # longer interval
    rate_sb_hi = 1 / (sb["interval_years"] - sb["interval_err_years"])
    ax.errorbar([1000 * 1.18], [rate_sb], yerr=[[rate_sb - rate_sb_lo], [rate_sb_hi - rate_sb]],
                fmt="D", color=ps.C_PUBLISHED, markersize=5.5, markeredgecolor=ps.SURFACE,
                elinewidth=1.4, zorder=6)
    ax.text(1300, 6e-7, f"Stuart & Binzel 2004:\n1 per {_sig(sb['interval_years'], 1)} years",
            fontsize=7.5, color=ps.INK, va="top", ha="left")
    pub = large["140 m"]["published"]["rates_per_year"]
    vals = np.array(list(pub.values()))
    ax.plot([140 * 1.18] * 2, [vals.min(), vals.max()], color=ps.C_PUBLISHED, linewidth=1.4,
            zorder=5)
    ax.scatter([140 * 1.18] * len(vals), vals, marker="D", s=26, color=ps.C_PUBLISHED,
               edgecolor=ps.SURFACE, linewidth=1.0, zorder=6)
    ax.text(140 * 1.36, 5.5e-5, "published estimates,\n2018–2025", fontsize=7.5,
            color=ps.INK, va="center")

    ax.xaxis.set_major_locator(FixedLocator([1, 10, 100, 1000, 10000]))
    ax.xaxis.set_major_formatter(FuncFormatter(_metres))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_locator(FixedLocator([100, 1, 1e-2, 1e-4, 1e-6, 1e-8]))
    ax.yaxis.set_major_formatter(FuncFormatter(_interval_label))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("Asteroid diameter (at least)")
    ax.set_ylabel("How often Earth is hit")
    ax.set_title("(c) How often Earth is hit by an asteroid at least this large", loc="left")

    fig.suptitle("From near misses to hits: counting close approaches measures the impact rate",
                 x=0.07, ha="left", y=0.975)
    ps.caption(fig, f"{_provenance(results)} Scaling assumes asteroids cross Earth's neighbourhood "
                    f"evenly spread across the plane of aim (so the number within distance d grows as "
                    f"d²), plus gravitational focusing, which enlarges Earth's target by "
                    f"{ir['focusing']['mean_factor_large']:.1f}× on average for km-sized asteroids. "
                    f"Fireball energies (JPL) and the published energy curves are turned into sizes for "
                    f"a stony density of 2,600 kg/m³ at {v_typ:.1f} km/s, the median fireball speed. "
                    f"Asteroid sizes from brightness, albedo 0.14 unless measured. Orange bars: 68% "
                    f"bootstrap ranges (resampling asteroids).")
    return _save(fig, "03_near_misses_to_impacts")


# ===========================================================================
# 04  The risk list
# ===========================================================================

_TORINO_BAND = np.array([0, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4])  # level -> colour band


def figure_risk(results: dict, art: dict) -> list:
    k = results["risk"]
    po = art["per_object"]

    fig = plt.figure(figsize=(12.5, 7.0))
    gs = GridSpec(1, 2, width_ratios=[1.35, 1], wspace=0.26, top=0.86, bottom=0.1, figure=fig)
    ax = fig.add_subplot(gs[0])
    axr = fig.add_subplot(gs[1])

    # --- (a) the Torino chart, computed from its definition -----------------
    lp = np.linspace(-11, 0, 700)
    le = np.linspace(-5.5, 8.5, 700)
    P, E = np.meshgrid(lp, le)
    level = nphys.torino_scale(10.0**P, 10.0**E)
    band = _TORINO_BAND[level]
    colours = [ps.TORINO_WASH[c] for c in ("white", "green", "yellow", "orange", "red")]
    ax.pcolormesh(10.0**lp, 10.0**le, band, cmap=ListedColormap(colours), vmin=-0.5, vmax=4.5,
                  alpha=0.2, shading="auto", zorder=0, rasterized=True)
    ax.contour(10.0**lp, 10.0**le, level, levels=np.arange(0.5, 10.5, 1.0), colors=ps.AXIS,
               linewidths=0.7, zorder=1)
    for text, (p, e) in {"0": (4e-11, 8e7), "1": (2e-6, 1e6), "2": (2e-4, 1e6),
                         "3": (0.1, 35), "4": (0.022, 5e3), "5": (0.3, 3e4), "6": (2e-3, 1e7),
                         "7": (0.12, 1e7)}.items():
        ax.text(p, e, text, fontsize=11, fontweight="bold", color=ps.INK_SECONDARY, ha="center",
                va="center", zorder=2)
    ax.text(0.1, 5, "2024 YR4\n(Jan–Feb 2025)", fontsize=6.5, color=ps.INK_SECONDARY,
            ha="center", va="center")
    ax.text(0.03, 600, "Apophis\n(Dec 2004)", fontsize=6.5, color=ps.INK_SECONDARY,
            ha="center", va="center")
    ax.annotate("8–10, certain collision (≥ 99%):\na sliver at the right edge", xy=(0.99, 1e8),
                xytext=(0.28, 1.6e8), fontsize=7, color=ps.INK_SECONDARY, ha="right",
                va="center", bbox=BACKING,
                arrowprops={"arrowstyle": "-", "color": ps.INK_SECONDARY, "linewidth": 0.8})
    ax.text(0.6, 0.6, "below 1 Mt: level 0\nat any probability", fontsize=7.5,
            color=ps.INK_SECONDARY, ha="right", va="top", bbox=BACKING)

    rated = po[~po["beyond_100yr"]]
    beyond = po[po["beyond_100yr"]]
    ax.scatter(rated["ip_max_date"], rated["energy_mt"], s=12, color=ps.C_DATA, alpha=0.55,
               linewidth=0, zorder=3, rasterized=True)
    ax.scatter(beyond["ip_max_date"], beyond["energy_mt"], s=44, facecolor=ps.SURFACE,
               edgecolor=ps.C_DATA, linewidth=1.6, zorder=4)
    for _, row in beyond.sort_values("energy_mt", ascending=False).iterrows():
        year = int(row["impact_year"])
        y_text = 1.5e6 if year > 2500 else 3e3
        ax.annotate(f"{_short_name(row['fullname'])}, {row['diameter_m']:,.0f} m:\n1 in "
                    f"{_sig(1 / row['ip_max_date'])} for {year},\nbeyond the scale's "
                    f"100 years\n(it would rate {int(row['torino_ours'])})",
                    xy=(row["ip_max_date"], row["energy_mt"]), xytext=(1.5e-11, y_text),
                    fontsize=7.5, color=ps.INK, va="center", ha="left", bbox=BACKING,
                    arrowprops={"arrowstyle": "-", "color": ps.INK_SECONDARY, "linewidth": 0.8,
                                "shrinkB": 4})
    top = po.loc[po["ip"].idxmax()]
    ax.annotate(f"{_short_name(top['fullname'])}: the most likely impact, "
                f"{top['ip_max_date']:.0%} in {int(top['impact_year'])},\nbut only "
                f"~{top['diameter_m']:.0f} m across: it would burst high in the air",
                xy=(top["ip_max_date"], top["energy_mt"]), xytext=(1e-7, 2e-5),
                fontsize=7.5, color=ps.INK, va="center", bbox=BACKING,
                arrowprops={"arrowstyle": "-", "color": ps.INK_SECONDARY, "linewidth": 0.8,
                            "shrinkB": 4})

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-11, 1)
    ax.set_ylim(10**-5.5, 10**8.5)
    p_ticks = {1e-9: "1 in a billion", 1e-6: "1 in a million", 1e-3: "1 in 1,000", 1: "certain"}
    ax.xaxis.set_major_locator(FixedLocator(list(p_ticks)))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: p_ticks.get(v, "")))
    e_ticks = {1e-4: "0.1 kt", 1e-2: "10 kt", 1: "1 Mt", 1e2: "100 Mt", 1e4: "10,000 Mt",
               1e6: "1 million Mt", 1e8: "100 million Mt"}
    ax.yaxis.set_major_locator(FixedLocator(list(e_ticks)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: e_ticks.get(v, "")))
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_minor_locator(NullLocator())
    ax.grid(False)
    ax.set_xlabel("Impact probability (each object's most likely date)")
    ax.set_ylabel("Impact energy (TNT equivalent)")
    ax.set_title(f"(a) All {k['objects']:,} objects on the list rate 0 on the Torino scale",
                 loc="left")

    # --- (b) objects removed each year -----------------------------------------
    by_year = pd.Series(k["removed_by_year"])
    years = by_year.index.to_numpy()
    last = years.max()
    full = years < last
    axr.bar(years[full], by_year.to_numpy()[full], width=0.72, color=ps.C_DATA,
            edgecolor=ps.SURFACE, linewidth=0.7, zorder=3)
    axr.bar([last], [by_year[last]], width=0.72, color=ps.C_CONTEXT, edgecolor=ps.SURFACE,
            linewidth=0.7, zorder=3)
    axr.text(last, by_year[last] + 6, "to\ndate", ha="center", va="bottom", fontsize=7,
             color=ps.INK_SECONDARY)
    removed = k["famous_removed"]
    peak_year = int(by_year.idxmax())
    placements = {"99942": (peak_year + 1.0, by_year.max() - 5, "left", "top"),
                  "2024 YR4": (last + 0.7, 215, "right", "bottom")}
    for des, name in (("99942", "Apophis"), ("2024 YR4", "2024 YR4")):
        date = removed.get(des)
        if not date:
            continue
        year = int(date[:4])
        x_t, y_t, ha, va = placements[des]
        axr.annotate(f"{name} removed\n{_fmt_date(date)}", xy=(year, by_year.get(year, 0)),
                     xytext=(x_t, y_t), fontsize=7.5, color=ps.INK, ha=ha, va=va,
                     arrowprops={"arrowstyle": "-", "color": ps.INK_SECONDARY, "linewidth": 0.8,
                                 "shrinkB": 3, "relpos": (1, 0) if ha == "right" else (0, 1)})
    axr.set_xlim(years.min() - 0.8, years.max() + 0.8)
    axr.set_ylim(0, by_year.max() * 1.08)
    axr.set_ylabel("Objects removed from the list")
    axr.text(0.02, 0.97,
             f"{k['removed_total']:,} removed since {k['removed_first'][:4]};\n"
             f"{k['objects']:,} listed now. New discoveries\njoin with uncertain orbits; "
             f"further\nobservations usually remove them.\n"
             f"{k['not_seen_for_5_years_fraction']:.0%} of today's list have not been\nseen for "
             f"five years or more.",
             transform=axr.transAxes, va="top", fontsize=8, color=ps.INK)
    axr.set_title("(b) Objects removed from the list, per year", loc="left")

    fig.suptitle("NASA's list of possible impacts: small, unlikely, and usually temporary",
                 x=0.07, ha="left", y=0.975)
    ps.caption(fig, f"{_provenance(results)} Sentry impact-monitoring table: "
                    f"{k['possible_impacts']:,} possible impact dates for {k['objects']:,} objects. "
                    f"Torino regions drawn from the scale's definition (Binzel 2000, revised 2005) "
                    f"and applied here to Sentry's own probabilities and energies; our ratings match "
                    f"Sentry's for all {k['torino_check']['rated_dates']:,} dates within 100 years. "
                    f"Past ratings of Apophis (level 4) and 2024 YR4 (level 3) are labels only.")
    return _save(fig, "04_risk_list")


# ===========================================================================
# 05  Apophis, 13 April 2029
# ===========================================================================

def _utc_from_tdb(jd_tdb) -> pd.DatetimeIndex:
    from astropy.time import Time
    return pd.DatetimeIndex(Time(np.asarray(jd_tdb), format="jd", scale="tdb").utc.datetime64)


def figure_apophis(results: dict, art: dict) -> list:
    a = results["apophis"]
    fine, cad, el = art["fine"], art["cad"], art["elements"]
    ca = art["ca"]
    c = a["closest_approach"]

    fig = plt.figure(figsize=(12.5, 11.5))
    gs = GridSpec(2, 2, height_ratios=[1.15, 1], hspace=0.26, wspace=0.2, top=0.92,
                  bottom=0.07, figure=fig)
    ax = fig.add_subplot(gs[0, 0])
    axo = fig.add_subplot(gs[0, 1])
    axh = fig.add_subplot(gs[1, :])

    # --- (a) the flyby, in its own plane, to scale ---------------------------
    r = fine[["x", "y", "z"]].to_numpy()
    v = fine[["vx", "vy", "vz"]].to_numpy()
    i = ca["index"]
    e1 = r[i] / np.linalg.norm(r[i])
    e2 = v[i] - (v[i] @ e1) * e1
    e2 /= np.linalg.norm(e2)
    px, py = r @ e1 / 1e3, r @ e2 / 1e3  # thousand km
    lim = 62
    ax.add_patch(Circle((0, 0), EARTH_RADIUS_KM / 1e3, color=ps.C_EARTH, zorder=3))
    ax.text(0, 0, "Earth", ha="center", va="center", fontsize=8, color=ps.INK, zorder=4)
    for radius, label, angle in ((GEO_RADIUS_KM, "geostationary distance, 42,164 km", 222),
                                 (GPS_RADIUS_KM, "GPS satellites, 26,560 km", 238)):
        ax.add_patch(Circle((0, 0), radius / 1e3, fill=False, edgecolor=ps.AXIS, linewidth=0.9,
                            zorder=2))
        t = np.radians(angle)
        ax.text(radius / 1e3 * np.cos(t), radius / 1e3 * np.sin(t), label, fontsize=7,
                color=ps.INK_SECONDARY, ha="center", va="center", rotation=angle + 90,
                bbox={"facecolor": ps.SURFACE, "edgecolor": "none", "pad": 0.5})
    ax.plot(px, py, color=ps.C_DATA, linewidth=2.2, zorder=4, solid_capstyle="round")
    # Hour marks: the grid point nearest each whole UTC hour (the grid is in TDB,
    # which runs 69 s ahead, so no grid point falls exactly on the hour).
    utc = _utc_from_tdb(fine["jd_tdb"])
    hours = pd.date_range(utc[0].ceil("h"), utc[-1].floor("h"), freq="h")
    for hour in hours:
        j = int(np.argmin(np.abs((utc - hour).total_seconds())))
        if abs(px[j]) < lim - 4 and abs(py[j]) < lim - 4:
            ax.scatter([px[j]], [py[j]], s=22, color=ps.C_DATA, edgecolor=ps.SURFACE,
                       linewidth=1.2, zorder=5)
            ax.text(px[j] + 2.4, py[j], hour.strftime("%H:%M UTC"), fontsize=7, color=ps.INK,
                    va="center")
    ax.scatter([px[i]], [py[i]], s=64, color=ps.C_DATA, edgecolor=ps.SURFACE, linewidth=1.5,
               zorder=6)
    ps.annotate(ax, f"closest, {c['time_utc'][11:16]} UTC:\n{c['distance_km']:,.0f} km from "
                    f"Earth's centre,\n{c['altitude_km']:,.0f} km above the surface.\n"
                    f"Uncertainty ±{c['horizons_3sigma_range_km']:.1f} km (3σ),\nfar thinner "
                    f"than this line",
                (px[i], py[i]), (8, -33), relpos=(0.5, 1), fontsize=7.5, color=ps.INK,
                ha="center", va="top", bbox=BACKING)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda val, _: f"{val:,.0f}"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda val, _: f"{val:,.0f}"))
    ax.set_xlabel("thousand km")
    g = a["geostationary"]
    ax.text(-lim + 2, lim - 2,
            f"Earth's gravity bends the path by {a['two_body']['measured_turn_deg_pm1day']:.0f}°.\n"
            f"Apophis passes {c['declination_deg']:.0f}° north of the equator,\nso it never "
            f"crosses the ring of geostationary\nsatellites: the nearest point of the ring\nis "
            f"{g['min_distance_to_ring_km']:,.0f} km away. The Moon is "
            f"{a['moon']['distance_km'] / 1e3:,.0f},000 km\naway that night, far off this chart.",
            fontsize=7.5, color=ps.INK, va="top", bbox=BACKING)
    ax.set_title("(a) The flyby, to scale, in the plane of Apophis's path", loc="left")

    # --- (b) the orbit before and after --------------------------------------
    nu = np.linspace(0, 360, 721)
    t = np.linspace(0, 2 * np.pi, 400)
    axo.plot(np.cos(t), np.sin(t), color=ps.INK_SECONDARY, linewidth=1.0, zorder=2)
    axo.text(np.cos(np.radians(-58)) * 1.03, np.sin(np.radians(-58)) * 1.03, "Earth's orbit",
             fontsize=7.5, color=ps.INK_SECONDARY, ha="left", va="top")
    for row, colour in ((el.iloc[0], ps.C_DATA), (el.iloc[-1], ps.C_FIT)):
        x, y, _ = nphys.orbit_xyz(row["a_au"], row["e"], row["i_deg"], row["node_deg"],
                                  row["peri_deg"], nu)
        axo.plot(x, y, color=colour, linewidth=2.0, zorder=3)
    axo.scatter([0], [0], s=90, marker="*", color=ps.YELLOW, edgecolor=ps.INK_SECONDARY,
                linewidth=0.6, zorder=4)
    axo.text(0.06, -0.06, "Sun", fontsize=7.5, color=ps.INK, va="top")
    near = el.iloc[(el["jd_tdb"] - ca["jd_tdb"]).abs().argmin()]
    fx, fy, _ = nphys.orbit_xyz(near["a_au"], near["e"], near["i_deg"], near["node_deg"],
                                near["peri_deg"], near["true_anomaly_deg"])
    axo.scatter([fx], [fy], s=40, color=ps.INK, edgecolor=ps.SURFACE, linewidth=1.2, zorder=5)
    ps.annotate(axo, "the flyby,\n13 Apr 2029", (fx, fy), (-1.47, -0.62), relpos=(1, 0.5),
                fontsize=7.5, color=ps.INK, ha="left", va="center")
    b, af = a["orbit_before"], a["orbit_after"]
    handles = [
        Line2D([], [], color=ps.C_DATA, linewidth=2,
               label=f"before: a {b['period_days']:.0f}-day year, a = {b['a_au']:.3f} au "
                     f"({b['class']}: mostly inside Earth's orbit)"),
        Line2D([], [], color=ps.C_FIT, linewidth=2,
               label=f"after: a {af['period_days']:.0f}-day year, a = {af['a_au']:.3f} au "
                     f"({af['class']}: mostly outside)"),
    ]
    axo.legend(handles=handles, loc="lower center", fontsize=7.5, labelspacing=0.8)
    axo.set_xlim(-1.5, 1.5)
    axo.set_ylim(-1.75, 1.45)
    axo.set_aspect("equal")
    axo.set_xlabel("au (ecliptic plane, seen from the north)")
    axo.set_title("(b) Earth's gravity moves Apophis to a new orbit", loc="left")

    # --- (c) every listed approach, with its 3-sigma range -------------------
    yrs = cad["time"].dt.year + (cad["time"].dt.dayofyear - 0.5) / 365.25
    dist, dmin, dmax = cad["dist_ld"], cad["dist_min_ld"], cad["dist_max_ld"]
    is2029 = ((cad["time"].dt.year == 2029) & (cad["dist_ld"] < 1)).to_numpy()
    axh.vlines(yrs, dmin, dmax, color=ps.C_DATA, linewidth=2.0, alpha=0.55, zorder=3)
    axh.scatter(yrs[~is2029], dist[~is2029], s=18, color=ps.C_DATA, edgecolor=ps.SURFACE,
                linewidth=0.8, zorder=4)
    axh.scatter(yrs[is2029], dist[is2029], s=70, color=ps.C_DATA, edgecolor=ps.SURFACE,
                linewidth=1.5, zorder=5)
    axh.set_yscale("log")
    axh.set_ylim(0.009, 400)
    axh.set_xlim(1900, 2125)
    axh.yaxis.set_major_locator(FixedLocator([0.01, 0.1, 1, 10, 100]))
    _plain_log(axh, "y")
    axh.set_ylabel("Miss distance (lunar distances)")
    _reference_line(axh, 1.0, "the Moon", x_text=1902, ha="left")
    _reference_line(axh, EARTH_LD, "Earth's surface", x_text=1902, ha="left")
    today = pd.Timestamp(results["snapshot"]["now_utc"])
    today_yr = today.year + (today.dayofyear - 0.5) / 365.25
    axh.axvline(today_yr, color=ps.INK_SECONDARY, linewidth=0.9, zorder=2)
    axh.text(today_yr - 0.8, 300, "today", fontsize=7.5, color=ps.INK_SECONDARY, va="top",
             ha="right")
    ps.annotate(axh, f"13 Apr 2029: {c['distance_ld']:.3f} LD, a tenth of the way to the Moon",
                (2029.28, c["distance_ld"]), (1990, 0.03), relpos=(1, 0.5), fontsize=7.5,
                color=ps.INK, ha="right", va="center")
    ca_list = a["approaches"]
    cl = ca_list["closest_after_2029"]
    ps.annotate(axh, f"{cl['date'][:4]}: anywhere from {cl['dist_min_ld']:.1f} to "
                     f"{cl['dist_max_ld']:.0f} LD (3σ),\ntiming uncertain by ±{cl['time_sigma_days']:.0f} "
                     f"days",
                (int(cl["date"][:4]) + 0.28, cl["dist_min_ld"]), (2111, 2.6), relpos=(1, 0.5),
                fontsize=7.5, color=ps.INK, ha="right", va="center", bbox=BACKING)
    rec = a["record"]
    axh.text(1903, 0.8,
             f"Why the 2029 distance is known to a few km: the orbit is fitted to "
             f"{rec['optical_observations']:,} telescope\npositions, {rec['radar_delay']} radar "
             f"ranges and {rec['radar_doppler']} radar speed measurements ({rec['first_obs'][:4]}–"
             f"{rec['last_obs'][:4]}). Even the push\nof sunlight is measured (the Yarkovsky effect, "
             f"which shrinks the orbit by {abs(rec['yarkovsky_drift_m_per_year']):.0f} m a year).\n"
             f"After 2029 the flyby amplifies every small uncertainty; JPL's list ends in "
             f"{ca_list['last'][:4]}.\nApophis left NASA's risk list on "
             f"{_fmt_date(results['risk']['famous_removed']['99942'])}.",
             ha="left", va="top", fontsize=7.5, color=ps.INK,
             bbox={"facecolor": ps.SURFACE, "edgecolor": ps.GRID, "pad": 5})
    axh.set_title("(c) Every Earth approach of Apophis within 0.5 au that JPL lists, "
                  "with its 3σ range", loc="left")

    fig.suptitle("Apophis, 13 April 2029: closer than geostationary satellites orbit, "
                 "and no danger", x=0.07, ha="left", y=0.975)
    ps.caption(fig, f"{_provenance(results)} Path from JPL Horizons (solution JPL#220), geometric "
                    f"positions every minute; times converted from TDB to UTC. Distances in (a) are to "
                    f"scale; the circles mark distances from Earth's centre, while geostationary "
                    f"satellites orbit over the equator, not in this plane. (b) osculating orbits on "
                    f"{b['date']} and {af['date']}. Bars in (c): JPL's 3σ minimum and maximum "
                    f"distances.")
    return _save(fig, "05_apophis_2029")


def make_all(results: dict, artefacts: dict) -> list:
    written = []
    written += figure_snapshot(results, artefacts["snapshot"])
    written += figure_encounters(results, artefacts["encounters"])
    written += figure_impacts(results, artefacts["impact_rates"])
    written += figure_risk(results, artefacts["risk"])
    written += figure_apophis(results, artefacts["apophis"])
    return written
