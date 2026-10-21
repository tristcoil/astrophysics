"""Every static figure, drawn from analysis results and artefacts (never recomputed).

Brightness is the dot's "apparent albedo" F in display-linear units (see
``lightcurve.py``), usually shown relative to its mean. Colour channels keep
their literal colours, each with its own line style and wavelength label.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from PIL import Image

import epicdata
import plotstyle as ps
from config import FIGURE_DIR, GALLERY_DIR
from transform import GAMMA_NOMINAL, GAMMA_RANGE, YANG2018

ps.use_style()

CH = ("red", "green", "blue")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def dot_colour(fr, fg, fb, gamma: float = GAMMA_NOMINAL, gain: float = 1.0):
    """The colour EPIC's camera would give Earth if it were one uniform pixel:
    the disc's mean linear light per channel, re-encoded for display."""
    lin = np.clip(gain * np.array([fr, fg, fb], float), 0, 1)
    return tuple(lin ** (1.0 / gamma))


def earth_image(image: str, size: int = 360, kind: str = "jpg") -> np.ndarray:
    with Image.open(epicdata.image_path(kind, image)) as im:
        return np.asarray(im.convert("RGB").resize((size, size), Image.LANCZOS))


def _image_axis(ax, img, title=None, fontsize=8):
    ax.imshow(img)
    ax.set_axis_off()
    if title:
        ax.set_title(title, fontsize=fontsize, fontweight="normal", loc="center", pad=3)


def _rel(x):
    x = np.asarray(x, float)
    return 100 * (x / np.nanmean(x) - 1)


def _channel_lines(ax, x, frame, cols, *, rel=True, lw=1.4, ms=0, alpha=1.0, labels=True):
    for ch in CH:
        col, ls, mk, lab = ps.CHANNEL_STYLE[ch]
        y = frame[cols[ch]]
        y = _rel(y) if rel else y
        ax.plot(x, y, color=col, ls=ls, lw=lw, marker=mk if ms else None, ms=ms,
                alpha=alpha, label=lab if labels else None)


# ---------------------------------------------------------------------------
# 1. one day of Earth, and the same day as a dot
# ---------------------------------------------------------------------------

def fig_one_day(results, art) -> list:
    jpg = art["jpg"]
    w = jpg[jpg["ok"] & (jpg["t"].dt.year == 2024)]
    per_day = w.groupby(w["t"].dt.date).size()
    day = per_day[per_day >= 21].index[len(per_day[per_day >= 21]) // 3]  # a full summer day
    d = w[w["t"].dt.date == day].sort_values("t")
    pick = np.unique(np.linspace(0, len(d) - 1, 8).round().astype(int))
    fig = plt.figure(figsize=(11.5, 6.4))
    gs = GridSpec(3, len(pick), height_ratios=[1.0, 0.36, 1.3], hspace=0.08, top=0.88, bottom=0.14,
                  figure=fig)
    for i, k in enumerate(pick):
        r = d.iloc[k]
        _image_axis(fig.add_subplot(gs[0, i]), earth_image(r["image"], 300),
                    f"{r['t']:%H:%M} UTC\nfacing {r['lon']:+.0f}°", fontsize=7.5)
        ax = fig.add_subplot(gs[1, i])
        ax.set_facecolor("black")
        ax.add_patch(plt.Circle((0.5, 0.5), 0.36, color=dot_colour(r["F_red"], r["F_green"], r["F_blue"])))
        ax.set_xlim(0, 1), ax.set_ylim(0, 1), ax.set_aspect("equal"), ax.set_xticks([]), ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        if i == 0:
            ax.set_ylabel("as a dot", fontsize=8, rotation=0, ha="right", va="center", labelpad=6)
    ax = fig.add_subplot(gs[2, :])
    hours = (d["t"] - pd.Timestamp(day)).dt.total_seconds() / 3600
    _channel_lines(ax, hours, d, {c: f"F_{c}" for c in CH}, ms=3)
    for k in pick:
        ax.axvline(hours.iloc[k], color=ps.GRID, lw=0.8, zorder=0)
    ax.set_xlim(0, 24)
    ax.set_xlabel(f"hours after 00:00 UTC on {day:%d %B %Y}")
    ax.set_ylabel("brightness of the dot (% from the day's mean)")
    ax.legend(loc="upper left", ncol=3)
    fig.suptitle("One day of Earth from 1.5 million km, and the same day as a single point of light",
                 x=0.125, ha="left", y=0.985)
    fig.text(0.125, 0.935, "Middle row: the colour EPIC's camera would record if the whole planet were one pixel",
             fontsize=8, color=ps.INK_SECONDARY)
    ps.caption(fig, "DSCOVR/EPIC natural-colour images (NASA EPIC Team). Brightness = the disc's total light "
               f"in display-linear units (gamma {GAMMA_NOMINAL}), divided by the area of the full disc. "
               "As Earth turns, Africa and Asia (bright land) replace the dark Pacific.", y=0.045)
    return ps.save(fig, FIGURE_DIR / "fig1_one_day_as_a_dot")


# ---------------------------------------------------------------------------
# 2. the length of the day
# ---------------------------------------------------------------------------

def fig_rotation(results, art) -> list:
    jpg = art["jpg"]
    w = jpg[jpg["ok"] & (jpg["t"].dt.year == 2024)]
    rot = results["rotation"]
    fig = plt.figure(figsize=(12, 8.2))
    gs = GridSpec(2, 2, height_ratios=[1, 1], hspace=0.42, wspace=0.22, figure=fig)
    ax = fig.add_subplot(gs[0, :])
    daily = w.groupby(w["t"].dt.floor("D"))[["F_red", "F_green", "F_blue"]].mean()
    ax.plot(w["t"], _rel(w["F_red"]), ls="none", marker="o", ms=1.2, color=ps.CONTEXT, zorder=1,
            label="every image, 680 nm")
    _channel_lines(ax, daily.index, daily, {c: f"F_{c}" for c in CH}, lw=1.3)
    ax.set_ylabel("brightness (% from the year's mean)")
    n = rot["2024"]["n_images"]
    ax.set_title(f"(a) 2024: {n:,} images of the dot. The daily means (lines) wander with the seasons; "
                 "single images swing around them as Earth turns", fontsize=10)
    ax.legend(loc="upper center", ncol=4)
    # (b) periodogram
    ax = fig.add_subplot(gs[1, 0])
    a = art["rotation"][("2024", "red")]
    ax.plot(a["f"], a["p"], color=ps.RED, lw=1.0, label="680 nm, 2024")
    for k, lab in ((1, "24 h"), (2, "12 h"), (3, "8 h"), (4, "6 h")):
        p = rot["2024"]["red"]["harmonics"][lab.replace(" ", "")]["power"]
        xy_text = (k + 0.12, 0.85 * p) if k == 1 else (k + 0.15, p + 0.022)
        ps.annotate(ax, lab, (k, p), xy_text, relpos=(0, 0.5))
    ax.set_xlim(0, 5)
    ax.set_xlabel("frequency (cycles per day)")
    ax.set_ylabel("Lomb–Scargle power")
    r24 = rot["2024"]["red"]
    r16 = rot["2015-16"]["red"]
    ax.set_title(f"(b) Periodogram: one peak at 24 h, and its harmonics", fontsize=10)
    ins = ax.inset_axes([0.5, 0.5, 0.47, 0.38])
    for name, ls in (("2024", "-"), ("2015-16", "--")):
        b = art["rotation"][(name, "red")]
        sel = (b["f"] > 0.985) & (b["f"] < 1.015)
        ins.plot(24 / b["f"][sel], b["p"][sel], color=ps.RED, ls=ls, lw=1.0, label=name)
    ins.axvline(24.0, color=ps.INK_SECONDARY, lw=0.8)
    ins.set_xlabel("period (h)", fontsize=7)
    ins.tick_params(labelsize=7)
    ins.legend(fontsize=6.5, loc="upper right")
    ins.set_title(f"P = {r24['period_h']:.3f} ± {r24['period_err_h']:.3f} h (2024)\n"
                  f"P = {r16['period_h']:.3f} ± {r16['period_err_h']:.3f} h (2015–16)",
                  fontsize=7, fontweight="normal", loc="left", pad=2)
    # (c) autocorrelation
    ax = fig.add_subplot(gs[1, 1])
    a = art["rotation"][("2024", "red")]
    lag_d = a["acf_lag_h"] / 24
    ax.plot(lag_d, a["acf"], color=ps.RED, lw=1.2)
    for k in range(1, 12):
        ax.axvline(k, color=ps.GRID, lw=0.8, zorder=0)
    ax.axhline(0, color=ps.AXIS, lw=0.8)
    ax.set_xlim(0, 12)
    ax.set_xlabel("time lag (days)")
    ax.set_ylabel("autocorrelation, 680 nm")
    ax.set_title(f"(c) Autocorrelation: it repeats every {r24['acf_peak_h']:.1f} h (r = {r24['acf_r']:.2f})",
                 fontsize=10)
    ps.caption(fig, "Lomb–Scargle periodograms of the 680 nm brightness; errors from leaving out one month at "
               "a time. The harmonics at 12, 8 and 6 h come from the uneven layout of continents and oceans, "
               "as Jiang et al. (2018) found.", y=0.01)
    return ps.save(fig, FIGURE_DIR / "fig2_length_of_the_day")


# ---------------------------------------------------------------------------
# 3. how rarely could we look?
# ---------------------------------------------------------------------------

def fig_sampling(results, art) -> list:
    s = results["sampling"]
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    jx = [24 / float(k) for k in s["jiang2018"]]
    ax.plot(jx, list(s["jiang2018"].values()), color=ps.C_PUBLISHED, marker="D", ms=5, lw=1.2,
            label="Jiang et al. 2018 (L1B 780 nm, 2015–2017)")
    styles = {"2015-16 random": (ps.C_OURS, "-", "o"), "2024 random": (ps.C_OURS, "--", "s"),
              "2015-16 regular": (ps.INK_SECONDARY, "-", "v"), "2024 regular": (ps.INK_SECONDARY, "--", "^")}
    for key, (col, ls, mk) in styles.items():
        rows = pd.DataFrame(s[key])
        what = "random times" if "random" in key else "fixed schedule"
        ax.plot(rows["per_day"], rows["p_detect"], color=col, ls=ls, marker=mk, ms=4, lw=1.3,
                label=f"this work, {key.split()[0]}, {what}")
    ax.axvline(2, color=ps.GRID, lw=1.0, zorder=0)
    ax.text(2.05, 0.03, "2 per day:\nNyquist limit\nfor a 24-h signal", fontsize=7.5, color=ps.MUTED)
    ax.set_xscale("log")
    ax.set_xticks([0.67, 1, 2, 3, 5, 8, 13])
    ax.set_xticklabels(["0.67", "1", "2", "3", "5", "8", "13"])
    ax.set_xlabel("images per day kept (one week of data)")
    ax.set_ylabel("chance the strongest peak is the 24-h rotation")
    ax.set_ylim(0, 1.03)
    ax.legend(loc="lower right", fontsize=7.5)
    ax.set_title("How rarely could an alien astronomer look and still find Earth's day?", fontsize=11)
    ps.caption(fig, "Each point: 300 random week-long stretches of the 680 nm light curve, thinned to the "
               "given rate. 'Random times' keeps a random subset (Jiang et al.'s procedure); 'fixed schedule' "
               "keeps one image per interval, which aliases below ~4 images a day. Detected = strongest "
               "peak within half a frequency-resolution element (1/14 per day) of one cycle per day.", y=-0.04)
    return ps.save(fig, FIGURE_DIR / "fig3_sampling")


# ---------------------------------------------------------------------------
# 4. continents and oceans in the dot's colour
# ---------------------------------------------------------------------------

def fig_continents(results, art) -> list:
    sur = results["surface"]
    fig = plt.figure(figsize=(12.5, 8.6))
    gs = GridSpec(3, 3, width_ratios=[1.35, 1.35, 1.0], height_ratios=[0.8, 1, 1],
                  hspace=0.45, wspace=0.3, figure=fig)
    p24 = art["surface"][("2024", "profile")]
    p16 = art["surface"][("2015-16", "profile")]
    x = p24.index.to_numpy()
    ax0 = fig.add_subplot(gs[0, :2])
    ax0.fill_between(x, 100 * p24[("land_flux", "mean")], color=ps.CONTEXT, step="mid")
    ax0.plot(x, 100 * p24[("land_flux", "mean")], color=ps.INK_SECONDARY, lw=1.0, drawstyle="steps-mid")
    ax0.set_ylabel("land's share\nof the light (%)")
    ax0.set_title("(a) Land's share of the light, by the longitude facing DSCOVR", fontsize=10)
    ax0.set_xlim(-180, 180)
    ax1 = fig.add_subplot(gs[1, :2], sharex=ax0)
    for p, ls, mk, lab in ((p24, "-", "o", "2024"), (p16, "--", "s", "2015–16")):
        ax1.plot(p.index, _rel(p[("blue_red", "mean")]), color=ps.C_OURS, ls=ls, marker=mk, ms=3.5,
                 label=lab)
    ax1.set_ylabel("blue / red light\n(% from mean)")
    pa = sur["2024"]["pacific_vs_africa"]["blue_red"]["pacific_over_africa"]
    ax1.set_title(f"(b) The dot's colour: the Pacific-facing Earth is {100 * (pa - 1):.0f}% bluer "
                  "than the Africa-facing one", fontsize=10)
    ax1.legend(loc="upper center", ncol=2)
    ax2 = fig.add_subplot(gs[2, :2], sharex=ax0)
    for ch in CH:
        col, ls, mk, lab = ps.CHANNEL_STYLE[ch]
        v = p24[(f"F_{ch}", "mean")].to_numpy()
        v30 = (np.roll(v, 1) + v + np.roll(v, -1)) / 3.0  # 30-degree running mean, wrapping
        ax2.plot(x, _rel(v30), color=col, ls=ls, marker=mk, ms=3.5, label=lab)
    ax2.set_ylabel("brightness, 30° mean\n(% from mean)")
    ax2.set_xlabel("longitude facing DSCOVR (degrees east)")
    pr = sur["2024"]["pacific_vs_africa"]["F_red"]["pacific_over_africa"]
    pb = sur["2024"]["pacific_vs_africa"]["F_blue"]["pacific_over_africa"]
    ax2.set_title(f"(c) Brightness, 2024: Pacific side {100 * (1 - pr):.0f}% fainter at 680 nm, "
                  f"{100 * (pb - 1):+.0f}% at 443 nm", fontsize=10)
    ax2.legend(loc="lower center", ncol=3)
    for a in (ax0, ax1, ax2):
        for lo in (-160, 20):
            a.axvspan(lo - 20, lo + 20, color=ps.GRID, alpha=0.5, zorder=0, lw=0)
    ax2.set_xticks(range(-180, 181, 60))
    # (d) colour vs land
    ax = fig.add_subplot(gs[1:, 2])
    w = art["jpg"]
    w = w[w["ok"] & (w["t"].dt.year == 2024)]
    ax.plot(100 * w["land_flux"], w["blue_red"], ls="none", marker="o", ms=1.5, color=ps.C_OURS, alpha=0.25)
    fit = sur["2024"]["blue_red_vs_land_fit"]
    xx = np.linspace(0, 0.6, 10)
    ax.plot(100 * xx, fit["intercept"] + fit["slope"] * xx, color=ps.INK, lw=1.2)
    r = sur["2024"]["corr_with_land"]["blue_red"]
    ax.set_xlabel("land's share of the light (%)")
    ax.set_ylabel(f"blue / red light (gamma {GAMMA_NOMINAL})")
    ax.set_title(f"(d) Colour tracks land: r = {r:.2f}\n(every 2024 image)", fontsize=10)
    # hemispheres, with their dots
    sub_gs = gs[0, 2].subgridspec(2, 2, height_ratios=[1, 0.35], hspace=0.05, wspace=0.1)
    for j, (rng, name) in enumerate((((-180, -140), "Pacific side"), ((0, 40), "Africa side"))):
        sub = w[(w["lon"] >= rng[0]) & (w["lon"] < rng[1]) & (w["t"].dt.month == 7)]
        r_ = sub.iloc[len(sub) // 2]
        _image_axis(fig.add_subplot(sub_gs[0, j]), earth_image(r_["image"], 220), name, fontsize=8)
        axd = fig.add_subplot(sub_gs[1, j])
        axd.set_facecolor("black")
        axd.add_patch(plt.Circle((0.5, 0.5), 0.4, color=dot_colour(r_["F_red"], r_["F_green"], r_["F_blue"])))
        axd.set_xlim(0, 1), axd.set_ylim(0, 1), axd.set_aspect("equal")
        axd.set_xticks([]), axd.set_yticks([])
        for sp in axd.spines.values():
            sp.set_visible(False)
    ps.caption(fig, "Means in 10° bins of the longitude at the centre of the disc. Shaded bands: the "
               "Pacific (−180° to −140°) and Africa (0° to 40°) views compared in the titles. Land share "
               "from NOAA's 0.25° land-sea mask, weighted by sunlight and viewing angle (cos θ₀ cos θ). (c) is "
               "smoothed over 30°: frame-to-frame tone adjustments (README) make single 10° bins jump by ~3%. "
               f"Linear light at gamma {GAMMA_NOMINAL}; the colour contrast is {100 * (sur['2024']['pacific_over_africa_g1.0']['blue_red'] - 1):.0f}% "
               f"in raw display units and {100 * (sur['2024'][f'pacific_over_africa_g{GAMMA_RANGE[1]:.1f}']['blue_red'] - 1):.0f}% at gamma {GAMMA_RANGE[1]}.",
               y=0.02)
    return ps.save(fig, FIGURE_DIR / "fig4_continents_in_the_colour")


# ---------------------------------------------------------------------------
# 5. clouds: rotation vs weather, and the principal components
# ---------------------------------------------------------------------------

def fig_weather(results, art) -> list:
    sur = results["surface"]
    fig = plt.figure(figsize=(12, 4.8))
    gs = GridSpec(1, 3, width_ratios=[1.1, 1, 1], wspace=0.32, figure=fig)
    ax = fig.add_subplot(gs[0])
    labels = ["680 nm", "551 nm", "443 nm", "blue/red"]
    keys = ["F_red", "F_green", "F_blue", "blue_red"]
    y = np.arange(len(keys))
    for off, name, mk in ((-0.17, "2015-16", "s"), (0.17, "2024", "o")):
        rot = [100 * sur[name]["rotation_vs_weather"][k]["rotation_share"] for k in keys]
        ax.barh(y + off, rot, height=0.3, color=ps.C_OURS if name == "2024" else ps.CONTEXT,
                label=f"{name}")
        for yi, v in zip(y + off, rot):
            ax.text(v + 1.5, yi, f"{v:.0f}%", va="center", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 112)
    ax.set_xlabel("share of the hour-to-hour flicker that repeats\nwith longitude (the rest is weather)")
    ax.set_title("(a) Rotation vs weather", fontsize=10)
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 1.0))
    pcs = art["surface"][("2024", "pcs")]
    pca = sur["2024"]["pca"]
    ax = fig.add_subplot(gs[1])
    ax.plot(100 * pcs["land"], pcs["pc2"], ls="none", marker="o", ms=1.5, color=ps.C_OURS, alpha=0.25)
    ax.set_xlabel("land's share of the light (%)")
    ax.set_ylabel("2nd principal component")
    ax.set_title(f"(b) PC2 ({100 * pca['explained'][1]:.0f}% of variance: red vs blue)\n"
                 f"is the surface: r = {pca['corr_pc2_land']:.2f} with land", fontsize=10)
    ax = fig.add_subplot(gs[2])
    ax.plot(100 * pcs["land"], pcs["pc1"], ls="none", marker="o", ms=1.5, color=ps.CONTEXT, alpha=0.4)
    ax.set_xlabel("land's share of the light (%)")
    ax.set_ylabel("1st principal component")
    ax.set_title(f"(c) PC1 ({100 * pca['explained'][0]:.0f}%: all colours together)\n"
                 f"is clouds and ice: r = {pca['corr_pc1_land']:+.2f} with land", fontsize=10)
    ps.caption(fig, "(a) After dividing out a 15-day running mean, the part of the variance explained by "
                    "the mean longitude profile (per quarter) vs the remainder; the remainder also holds the "
                    "~3% (linear) frame-to-frame tone noise of the images. (b, c) Principal components "
                    "of the three standardised log-brightnesses, 2024 (Fan et al. 2019 found the same split "
                    "with ten bands).", y=-0.05)
    return ps.save(fig, FIGURE_DIR / "fig5_clouds_and_weather")


# ---------------------------------------------------------------------------
# 6. seasons
# ---------------------------------------------------------------------------

def fig_seasons(results, art) -> list:
    track = art["track"]
    fig = plt.figure(figsize=(12.5, 8.4))
    gs = GridSpec(3, 2, width_ratios=[2.3, 1], height_ratios=[1, 1, 0.8], hspace=0.55, wspace=0.12,
                  figure=fig)
    months = ([1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335],
              ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    ax = fig.add_subplot(gs[0, 0])
    daily = track.groupby(track["t"].dt.floor("D"))["lat"].mean()
    ax.plot(daily.index, daily, ls="none", marker="o", ms=1.2, color=ps.C_OURS)
    ax.axhline(23.44, color=ps.GRID, lw=0.8), ax.axhline(-23.44, color=ps.GRID, lw=0.8)
    ax.set_ylabel("latitude at disc centre (°)")
    ax.set_title("(a) The latitude facing DSCOVR, 2015–2026: north tilts in every June,\n"
                 "south every December", fontsize=10)
    ax = fig.add_subplot(gs[1, 0])
    for name, ls in (("2024", "-"), ("2015-16", "--")):
        d = art["seasons"][name]
        doy = np.asarray(d.index.dayofyear)
        order = np.argsort(doy)
        for ch in ("red", "blue"):
            col, _, mk, lab = ps.CHANNEL_STYLE[ch]
            y = d[f"F_{ch}"].rolling(7, center=True, min_periods=3).mean()
            ax.plot(doy[order], _rel(y.to_numpy())[order], color=col, ls=ls, lw=1.2, label=f"{lab}, {name}")
    ax.set_xlim(1, 366)
    ax.set_xticks(*months)
    ax.set_ylabel("daily-mean brightness\n(% from mean, 7-day mean)")
    ax.set_title("(b) Brightest in December–January, a second peak in May–June\n"
                 "(as Yang et al. 2018 found for 2015–16)", fontsize=10)
    ax.legend(loc="upper center", ncol=2, fontsize=7)
    ax = fig.add_subplot(gs[2, 0])
    for name, ls in (("2024", "-"), ("2015-16", "--")):
        d = art["seasons"][name]
        doy = np.asarray(d.index.dayofyear)
        order = np.argsort(doy)
        ax.plot(doy[order], 100 * d["polar"].to_numpy()[order], color=ps.INK_SECONDARY, ls=ls, lw=1.1,
                label=name)
    ax.set_xlim(1, 366)
    ax.set_xticks(*months)
    ax.set_ylabel("polar share\nof the light (%)")
    ax.set_title("(c) Light from poleward of 60°: the ice caps swing into view at the solstices", fontsize=10)
    ax.legend(loc="upper center", ncol=2)
    jpg = art["jpg"]
    for i, (month, day_, label) in enumerate(((6, 21, "21 June 2024: the Arctic in view"),
                                              (12, 21, "21 December 2024: Antarctica in view"))):
        sub = jpg[(jpg["t"].dt.year == 2024) & (jpg["t"].dt.month == month) & (jpg["t"].dt.day == day_)]
        r = sub.iloc[np.argmin(np.abs(sub["lon"] - 20))]
        ax = fig.add_subplot(gs[i, 1])
        _image_axis(ax, earth_image(r["image"], 320))
        ax.text(0.5, -0.04, f"{label}\n(disc centre at {r['lat']:+.0f}°)", transform=ax.transAxes,
                ha="center", va="top", fontsize=8, color=ps.INK_SECONDARY)
    ps.caption(fig, "The tilt seen by DSCOVR combines Earth's 23.4° axial tilt with DSCOVR's own loop "
                    "around L1. Polar share: as the land share in Fig. 4, for latitudes beyond ±60°.", y=0.02)
    return ps.save(fig, FIGURE_DIR / "fig6_seasons")


# ---------------------------------------------------------------------------
# 7. phase angle and DSCOVR's orbit
# ---------------------------------------------------------------------------

def fig_phase(results, art) -> list:
    track = art["track"]
    ph_ = results["phase"]
    fig = plt.figure(figsize=(12.5, 4.9))
    gs = GridSpec(1, 3, width_ratios=[1, 1.25, 1], wspace=0.42, figure=fig)
    ax = fig.add_subplot(gs[0])
    tt = track.iloc[::5]
    yr = tt["t"].dt.year + tt["t"].dt.dayofyear / 366
    sc = ax.scatter(tt["y_km"] / 1e3, tt["z_km"] / 1e3, c=yr, s=1.5, cmap="Blues", vmin=2013, vmax=2027)
    ax.plot(0, 0, marker="+", color=ps.INK, ms=10)
    ax.text(5, 8, "Sun–Earth line", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_aspect("equal")
    ax.set_xlabel("in the ecliptic (thousand km)")
    ax.set_ylabel("north–south (thousand km)")
    ax.set_title("(a) DSCOVR's loop around L1,\nseen from Earth", fontsize=10)
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("year")
    ax = fig.add_subplot(gs[1])
    daily = track.groupby(track["t"].dt.floor("D"))["phase_deg"].mean()
    daily = daily.asfreq("D")  # gaps stay gaps
    ax.plot(daily.index, daily, lw=0.7, color=ps.C_OURS)
    g = results["geometry"]["phase_deg"]
    ax.set_ylabel("phase angle (°)")
    ax.set_title(f"(b) The Sun–Earth–DSCOVR angle: {g['min']:.1f}°–{g['max']:.1f}°,\n"
                 f"cycling every ~{results['geometry']['phase_period_d']:.0f} days", fontsize=10)
    ax = fig.add_subplot(gs[2])
    d = art["phase"]
    yv = 100 * (np.exp(d["F_red_season_free"]) / np.exp(d["F_red_season_free"]).mean() - 1)
    ax.plot(d["phase"], yv, ls="none", marker="o", ms=1.5, color=ps.C_OURS, alpha=0.3)
    bins = np.arange(1.5, 13, 1.0)
    idx = np.digitize(d["phase"], bins)
    med = [np.median(yv[idx == i]) if (idx == i).sum() > 20 else np.nan for i in range(1, len(bins))]
    ax.plot(0.5 * (bins[1:] + bins[:-1]), med, color=ps.INK, marker="o", ms=4, lw=0, label="median per degree")
    s = ph_["F_red"]
    xx = np.linspace(2, 12.4, 20)
    ax.plot(xx, 100 * s["slope_per_deg"] * (xx - d["phase"].mean()), color=ps.INK, lw=1.2,
            label=f"fit: {100 * s['slope_per_deg']:+.2f} ± {100 * s['slope_err']:.2f} %/°")
    from geometry import lambert_phase
    lam = lambert_phase(xx)
    ax.plot(xx, 100 * (lam / lambert_phase(np.array([d['phase'].mean()]))[0] - 1), color=ps.C_PUBLISHED,
            lw=1.2, ls="--", label=f"Lambert sphere: {100 * ph_['lambert_slope_per_deg_at_median']:+.2f} %/°")
    ax.set_xlabel("phase angle (°)")
    ax.set_ylabel("daily brightness, seasons removed (%)")
    ax.set_ylim(-12, 12)
    ax.set_title("(c) Brighter at small phase angles\n(680 nm, 2015–2026, display units)", fontsize=10)
    ax.legend(loc="lower left", fontsize=7)
    ps.caption(fig, "(a) DSCOVR's position across the Sun–Earth line from the EPIC API's J2000 vectors. "
                    "(c) Daily means of 120-pixel thumbnails in display units, after a regression on the "
                    "phase angle, three seasonal harmonics, a linear trend and the processing version.", y=-0.04)
    return ps.save(fig, FIGURE_DIR / "fig7_phase_angle")


# ---------------------------------------------------------------------------
# 8. solar or sidereal day?
# ---------------------------------------------------------------------------

def fig_day_length(results, art) -> list:
    lg = results["long"]
    fig = plt.figure(figsize=(12, 4.6))
    gs = GridSpec(1, 2, width_ratios=[1.1, 1], wspace=0.25, figure=fig)
    ax = fig.add_subplot(gs[0])
    a = art["long"]["red"]
    period_s = 86400 / a["f"]
    ax.plot(period_s - 86400, a["p"], color=ps.RED, lw=1.0)
    sid = (lg["red"]["sidereal_period_h"] - 24) * 3600
    ax.axvline(0, color=ps.INK, lw=0.9)
    ax.axvline(sid, color=ps.INK_SECONDARY, lw=0.9, ls="--")
    ax.text(5, 0.9 * a["p"].max(), "solar day\n24 h 00 m 00 s", color=ps.INK_SECONDARY, fontsize=8)
    ax.text(sid + 5, 0.55 * a["p"].max(), "sidereal day\n23 h 56 m 04 s", color=ps.INK_SECONDARY, fontsize=8)
    ax.set_xlim(-300, 120)
    ax.set_xlabel("period minus 24 h (seconds)")
    ax.set_ylabel("Lomb–Scargle power, 680 nm")
    o, e = lg["red"]["offset_s"], lg["red"]["offset_err_s"]
    ax.set_title(f"(a) {lg['span_days'] / 365.25:.1f} years: one turn every 24 h {o:+.1f} ± {e:.1f} s",
                 fontsize=10)
    ax = fig.add_subplot(gs[1])
    b = art["long"][("red", "broad")]
    ax.plot(b["f"], b["p"], color=ps.RED, lw=0.6)
    ax.set_xscale("log")
    ax.set_xlim(1 / 600, 4.5)
    for f_, lab in ((1 / 365.25, "1 year"), (2 / 365.25, "½ year"),
                    (1 / results["geometry"]["phase_period_d"], f"{results['geometry']['phase_period_d']:.0f} d: phase angle"),
                    (1, "1 day"), (2, "12 h"), (3, "8 h")):
        ax.axvline(f_, color=ps.GRID, lw=0.8, zorder=0)
        ax.text(f_ * 1.06, 0.92 * b["p"].max(), lab, fontsize=7.5, color=ps.INK_SECONDARY, rotation=90, va="top")
    ax.set_xlabel("frequency (cycles per day)")
    ax.set_ylabel("Lomb–Scargle power")
    ax.set_title("(b) The whole spectrum: orbit, seasons and rotation", fontsize=10)
    n_res = (24 - lg["red"]["sidereal_period_h"]) * 3600 / lg["red"]["resolution_s"]
    ps.caption(fig, "DSCOVR keeps station on the Sun–Earth line, so it follows Earth around the Sun and "
                    "sees one turn per solar day. An observer fixed among the stars would see the sidereal "
                    f"day, 236 s shorter: {n_res:.0f} frequency-resolution elements away at this baseline. "
                    "Errors from leaving out one calendar year at a time. 120-pixel thumbnails, display units.",
               y=-0.03)
    return ps.save(fig, FIGURE_DIR / "fig8_solar_or_sidereal_day")


# ---------------------------------------------------------------------------
# 9. the Moon
# ---------------------------------------------------------------------------

def _crop_transit(path, cx, cy, half):
    with Image.open(path) as im:
        im = im.convert("RGB")
        box = (int(cx - half), int(cy - half), int(cx + half), int(cy + half))
        return np.asarray(im.crop(box).resize((400, 400), Image.LANCZOS))


def fig_moon(results, art) -> list:
    from moon import TRANSITS
    mo = results["moon"]
    tab = art["moon"]["transits"]
    fig = plt.figure(figsize=(12.5, 8.6))
    gs = GridSpec(2, 3, height_ratios=[1, 1], hspace=0.42, wspace=0.3, figure=fig)
    for i, key in enumerate(("2015", "2016", "2021")):
        t = tab[(tab["transit"] == key) & (tab["n_moon_on_space"] == 0)]
        r = t.iloc[len(t) // 2]
        s = mo["transits"][key]
        img = _crop_transit(GALLERY_DIR / TRANSITS[key] / r["file"], r["earth_x_px"], r["earth_y_px"],
                            1.08 * r["earth_r_px"])
        _image_axis(fig.add_subplot(gs[0, i]), img,
                    f"{r['t']:%d %b %Y %H:%M} UTC · Moon at phase {s['moon_phase_deg']:.1f}°", fontsize=8.5)
    # (b) display values
    ax = fig.add_subplot(gs[1, 0])
    y = np.arange(3)
    s21 = mo["transits"]["2021"]
    for off, what, col in ((-0.17, "earth", ps.CONTEXT), (0.17, "moon", ps.C_OURS)):
        vals = [s21[f"display_{ch}"][what] for ch in CH]
        ax.barh(y + off, vals, height=0.3, color=col, label="Earth (disc mean)" if what == "earth" else "Moon (far side)")
        for yi, v in zip(y + off, vals):
            ax.text(v + 2, yi, f"{v:.0f}", va="center", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_yticks(y, ["680 nm", "551 nm", "443 nm"])
    ax.invert_yaxis()
    ax.set_xlim(0, 125)
    ax.set_xlabel("display value (0–255), 11 Feb 2021")
    ax.set_title("(a) In the pictures: the Moon is darker\nand redder than Earth", fontsize=10)
    ax.legend(loc="lower right", fontsize=7.5)
    # (c) ratio vs gamma
    ax = fig.add_subplot(gs[1, 1])
    for ch in ("red", "blue"):
        col, ls, mk, lab = ps.CHANNEL_STYLE[ch]
        rr = s21[f"moon_over_earth_{ch}"]
        g = np.array([float(k) for k in rr])
        ax.plot(g, list(rr.values()), color=col, ls=ls, marker=mk, ms=4, label=f"measured, {lab}")
        exp = mo["expected"][ch]["ratio"]
        ax.axhline(exp, color=col, lw=0.9, ls=ls, alpha=0.6)
        ax.text(2.62, exp, f"calibrated\nEPIC {exp:.2f}", fontsize=7, color=ps.INK_SECONDARY, va="center")
    ax.axvspan(*GAMMA_RANGE, color=ps.GRID, alpha=0.6, lw=0, zorder=0)
    ax.set_xlim(0.95, 3.0)
    ax.set_xlabel("display exponent gamma assumed")
    ax.set_ylabel("Moon / Earth, linear light")
    ax.set_title("(b) Back in linear light, the ratio depends\non the unknown display curve", fontsize=10)
    ax.legend(loc="upper right", fontsize=7)
    # (d) event timeline
    ax = fig.add_subplot(gs[1, 2])
    ev = art["moon"]["events"]
    kinds = {"transit": (ps.C_OURS, "o", 1.0), "occultation": (ps.INK_SECONDARY, "s", 0.0),
             "near miss": (ps.CONTEXT, ".", -1.0)}
    for k, (col, mk, yv) in kinds.items():
        e = ev[ev["kind"] == k]
        ax.plot(e["t_min_sep"], np.full(len(e), yv), ls="none", marker=mk, ms=7, color=col,
                mec=ps.SURFACE, mew=0.6, label=f"{k} ({len(e)})")
    public = pd.to_datetime(["2015-07-16", "2016-07-05", "2021-02-11", "2015-09-27", "2020-10-02", "2025-05-11"])
    for p in public:
        e = ev.iloc[np.argmin(np.abs(ev["t_min_sep"] - p))]
        yv = kinds[e["kind"]][2]
        ax.plot(e["t_min_sep"], yv, marker="o", ms=12, mfc="none", mec=ps.INK, mew=1.0)
    n_k = ev["kind"].value_counts()
    ax.set_yticks([1, 0, -1], [f"in front\n(transit, {n_k.get('transit', 0)})",
                               f"behind\n(occultation, {n_k.get('occultation', 0)})",
                               f"in the frame\nonly ({n_k.get('near miss', 0)})"])
    ax.set_ylim(-1.6, 1.6)
    ax.set_title("(c) Every time the Moon entered EPIC's view,\npredicted from JPL Horizons", fontsize=10)
    ax.text(0.01, 0.97, "ringed: in EPIC's public galleries", transform=ax.transAxes, fontsize=7,
            color=ps.INK_SECONDARY, va="top")
    sl = mo["stray_light"].get("red", {})
    ps.caption(fig, "Frames from the EPIC team's galleries (natural colour). DSCOVR is sunward of Earth, so in "
                    "a transit it sees the Moon's sunlit far side. (b) Calibrated values: EPIC's far side "
                    "(Gorkavyi et al. 2026) over Earth's July-2016 mean (Yang et al. 2018); the band is our "
                    f"adopted gamma range. Stray light inside the camera brightens Moon terrain in front of "
                    f"Earth by {100 * sl.get('excess_display', np.nan):.0f}% (display) over the same terrain "
                    "in front of space.", y=0.01)
    return ps.save(fig, FIGURE_DIR / "fig9_moon")


# ---------------------------------------------------------------------------
# 10. the display transform
# ---------------------------------------------------------------------------

def fig_transform(results, art) -> list:
    tr = results["transform"]
    fig = plt.figure(figsize=(12.5, 4.8))
    gs = GridSpec(1, 3, width_ratios=[1, 1, 1.15], wspace=0.33, figure=fig)
    ax = fig.add_subplot(gs[0])
    vr, vb = art["transform"]["mixing_vr"], art["transform"]["mixing_vb"]
    H, xe, ye = np.histogram2d(vr, vb, bins=[np.arange(0, 257, 2)] * 2)
    ax.imshow(np.log10(H.T + 1), origin="lower", extent=[0, 256, 0, 256], cmap="Greys", aspect="auto")
    ridge = np.array(tr["mixing_blue_vs_red"]["ridge"])
    ax.plot(ridge[:, 0], ridge[:, 1], color=ps.C_OURS, marker="o", ms=2.5, lw=1.0, label="median ridge")
    o, c = ridge[1], ridge[-1]
    ax.plot([o[0], c[0]], [o[1], c[1]], color=ps.INK_SECONDARY, lw=0.9, ls="--", label="straight (gamma = 1)")
    ax.set_xlabel("red display value")
    ax.set_ylabel("blue display value")
    ax.set_xlim(0, 256), ax.set_ylim(0, 256)
    g = tr["mixing_blue_vs_red"]["gamma"]
    ax.set_title(f"(a) Ocean–cloud mixtures bend:\ngamma = {g:.2f} fits", fontsize=10)
    ax.legend(loc="upper left", fontsize=7)
    ax = fig.add_subplot(gs[1])
    gs_ = np.array([float(k) for k in tr["yang_comparison"]])
    for ch in ("red", "blue"):
        col, ls, mk, lab = ps.CHANNEL_STYLE[ch]
        v = [100 * tr["yang_comparison"][k][ch]["rel_std"] for k in tr["yang_comparison"]]
        ax.plot(gs_, v, color=col, ls=ls, marker=mk, ms=4, label=f"this work, {lab}")
        ax.axhline(100 * YANG2018["daily_rel_std"][ch], color=ps.C_PUBLISHED, ls=ls, lw=1.0)
    ax.text(1.02, 8.2, "Yang et al. 2018, 680 nm", fontsize=7, color=ps.C_PUBLISHED)
    ax.text(1.02, 4.2, "Yang et al. 2018, 443 nm", fontsize=7, color=ps.C_PUBLISHED)
    ax.set_xlabel("display exponent gamma assumed")
    ax.set_ylabel("day-to-day scatter of the dot (%)")
    ax.set_title("(b) Calibrated variability,\n20–31 July 2016", fontsize=10)
    ax.legend(loc="lower right", fontsize=7)
    ax = fig.add_subplot(gs[2])
    rows = []
    rows.append(("ocean–cloud mixing (blue)", tr["mixing_blue_vs_red"]["gamma"], tr["mixing_blue_vs_red"]["gamma"]))
    rows.append(("ocean–cloud mixing (green)", tr["mixing_green_vs_red"]["gamma"], tr["mixing_green_vs_red"]["gamma"]))
    oc = tr["ocean_cloud_colour"]
    rows.append(("ocean vs cloud colour", *oc["gamma_range"]))
    ys, yp = tr["gamma_from_yang_std"], tr["gamma_from_yang_p2p"]
    rows.append(("Yang 2018 scatter", min(v for v in ys.values() if v), max(v for v in ys.values() if v)))
    yv = [v for v in yp.values() if v]
    if yv:
        rows.append(("Yang 2018 peak-to-peak", min(yv), max(yv)))
    do = tr["dark_ocean"]
    rows.append(("darkest ocean (Rayleigh)", min(d["gamma_min"] for d in do.values()),
                 max(d["gamma_max"] for d in do.values())))
    gm = [v for k in results["moon"]["gamma_matching_moon"].values() for v in k.values() if v]
    if gm:
        rows.append(("the Moon (EPIC far side)", min(gm), max(gm)))
    else:
        rows.append(("the Moon (EPIC far side)", 0.95, 1.15))
    for i, (lab, lo, hi) in enumerate(rows):
        col = ps.C_PUBLISHED if "Moon" in lab else ps.C_OURS
        ax.plot([lo, hi], [i, i], color=col, lw=4, solid_capstyle="round")
        ax.plot([(lo + hi) / 2], [i], marker="o", color=col, ms=5)
        ax.text(0.85, i - 0.28, lab, fontsize=7.5, color=ps.INK_SECONDARY, va="bottom")
    ax.axvspan(*GAMMA_RANGE, color=ps.GRID, alpha=0.6, lw=0, zorder=0)
    ax.axvline(GAMMA_NOMINAL, color=ps.INK, lw=1.0)
    ax.set_yticks([])
    ax.set_ylim(len(rows) - 0.5, -0.9)
    ax.set_xlim(0.8, 3.2)
    ax.set_xlabel("display exponent gamma")
    ax.set_title(f"(c) Adopted: gamma = {GAMMA_NOMINAL} (band {GAMMA_RANGE[0]}–{GAMMA_RANGE[1]});\n"
                 "the Moon disagrees", fontsize=10)
    ps.caption(fig, "The natural-colour images are 'color and brightness adjusted' (NASA EPIC) by an "
                    "unpublished curve. We model it as display = (linear light)^(1/gamma) and test gamma "
                    "against physics and against calibrated published numbers (README, 'Display images are "
                    "not photometry').", y=-0.04)
    return ps.save(fig, FIGURE_DIR / "fig10_display_transform")


def make_all(results, art) -> list:
    written = []
    for f in (fig_one_day, fig_rotation, fig_sampling, fig_continents, fig_weather, fig_seasons,
              fig_transform, fig_moon):
        written += f(results, art)
    if "phase" in results and "long" in results:
        written += fig_phase(results, art)
        written += fig_day_length(results, art)
    return written
