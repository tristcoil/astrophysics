"""Figure generation.

Every figure is drawn from the samples and the analysis outputs; nothing here
re-derives a result, so what is plotted is exactly what was measured.

Six figures, in reading order:

    01  discoveries per year, by method   -- the history of the field
    02  radius vs orbital period          -- the shape of the known population
    03  mass vs orbital distance, by method -- why each method finds what it finds
    04  the radius valley                 -- where it is, what hides it, how it tilts
    05  mass vs radius                    -- what small planets are made of
    06  the habitable zone                -- who receives Earth-like starlight

Design notes: log axes carry plain-number tick labels (0.1, 1, 10, 1,000) rather
than powers of ten; the whole population is drawn as recessive grey wherever a
subset is the subject; text is always ink, never the colour of a data series;
and every figure carries a caption with its provenance and the caveat a reader
needs in order not to over-read it.
"""

from __future__ import annotations

from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

import planet_physics as pp
import plotstyle as ps
from config import (
    CONTROVERSIAL_FAMOUS, FIGURE_DIR, HZ_HIGHLIGHT, MR_HIGHLIGHT, METHOD_ORDER,
    VAN_EYLEN_INTERCEPT, VAN_EYLEN_SLOPE, ZENG_MODELS, ZENG_TEMPERATURE_K,
)
from fetch_data import load_model

ps.use_style()

#: What the single paper behind each tall bar did. Keyed by the archive's own
#: reference name, so a future archive update that changes the spikes simply
#: falls back to the bare reference.
SPIKE_NOTES = {
    "Rowe et al. 2014": "Kepler planets in multi-planet\nsystems, validated in one paper",
    "Morton et al. 2016": "Kepler candidates validated\nstatistically, in one paper",
    "Valizadegan et al. 2022": "Kepler candidates validated by\nthe ExoMiner neural network",
}


def _provenance(meta) -> str:
    return (f"Data: NASA Exoplanet Archive, Planetary Systems Composite Parameters "
            f"(doi:10.26133/NEA13), retrieved {meta['retrieved_utc'][:10]}.")


def _num(x, _pos=None) -> str:
    """Plain-number tick labels for log axes: 0.01, 0.1, 1.5, 10, 1,000."""
    if x >= 1 and float(x).is_integer():
        return f"{x:,.0f}"
    return f"{x:.10f}".rstrip("0").rstrip(".")


def _plain_log(ax, which="both"):
    for axis in ([ax.xaxis, ax.yaxis] if which == "both" else [getattr(ax, f"{which}axis")]):
        axis.set_major_formatter(FuncFormatter(_num))
        axis.set_minor_formatter(NullFormatter())


def _decades(lo_exp: int, hi_exp: int):
    return FixedLocator([10.0**k for k in range(lo_exp, hi_exp + 1)])


def _solar_system(ax, xkey: str, ykey: str, *, labels: bool = True, names=None,
                  size: float = 26, offsets: dict | None = None, fontsize: float = 7.5):
    """Solar System planets as ink diamonds -- the familiar yardstick."""
    offsets = offsets or {}
    for name, p in pp.SOLAR_SYSTEM.items():
        if names is not None and name not in names:
            continue
        x, y = p[xkey], p[ykey]
        ax.scatter([x], [y], marker="D", s=size, color=ps.C_SOLAR,
                   edgecolor=ps.SURFACE, linewidth=0.9, zorder=7)
        if labels:
            dx, dy = offsets.get(name, (5, 3))
            ax.annotate(name, (x, y), xytext=(dx, dy), textcoords="offset points",
                        fontsize=fontsize, color=ps.INK, zorder=8,
                        ha="left" if dx >= 0 else "right")


def _save(fig, stem: str) -> list:
    return ps.save(fig, FIGURE_DIR / stem)


# ===========================================================================
# 01  Discoveries per year
# ===========================================================================

def figure_discoveries(summary: dict, art: dict, meta: dict) -> list:
    """Stacked bars: planets discovered each year, by method.

    The three tallest bars are each dominated by a single paper that validated
    hundreds of Kepler candidates at once -- the history of the field is partly
    a history of statistics.
    """
    counts = art["counts"]
    years = counts.index.to_numpy()
    fig, ax = plt.subplots(figsize=(10, 5.6))

    bottom = np.zeros(len(years))
    for method in METHOD_ORDER:
        values = counts[method].to_numpy()
        ax.bar(years, values, bottom=bottom, width=0.76, color=ps.METHOD_COLORS[method],
               edgecolor=ps.SURFACE, linewidth=0.7,
               label=f"{method}  ({summary['by_method'][method]:,})")
        bottom += values
    totals = bottom

    ax.set_xlim(years[0] - 0.8, years[-1] + 0.8)
    ax.set_ylim(0, totals.max() * 1.16)
    ax.set_ylabel("Planets discovered that year")
    ax.grid(axis="x", visible=False)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.set_title(f"{summary['planets']:,} confirmed exoplanets, discovered "
                 f"{summary['first_year']}–{summary['last_year']}")

    # --- the three spikes, each mostly one paper -------------------------------
    # (text x offset in years, text top as a fraction of the tallest bar,
    #  alignment, leader start on the text box) for the spikes in year order.
    peak = totals.max()
    layout = [(-0.9, 0.62, "right", (1, 1)), (-1.3, 0.93, "right", (1, 1)),
              (-3.4, 0.80, "left", (0.55, 0))]
    for spike, (dx, fy, ha, rel) in zip(summary["spikes"], layout):
        y, n = spike["year"], spike["planets"]
        note = SPIKE_NOTES.get(spike["reference"], "")
        text = (f"{y}: {spike['reference']}\n{spike['reference_planets']:,} of the year's "
                f"{n:,}" + (f" —\n{note}" if note else ""))
        ps.annotate(ax, text, xy=(y, n), xytext=(y + dx, peak * fy), relpos=rel,
                    ha=ha, va="top", linespacing=1.35)

    # --- the early history, too small to see as bars ---------------------------
    early = summary["earliest"]
    first_rv = summary["first_by_method"].get("Radial velocity")
    ps.annotate(ax, f"{early['year']}: the first exoplanets —\norbiting a pulsar "
                f"({early['planet'].rsplit(' ', 1)[0]})",
                xy=(early["year"], 6), xytext=(early["year"] - 0.15, peak * 0.34),
                relpos=(0, 0), ha="left", va="bottom", linespacing=1.35)
    if first_rv:
        ps.annotate(ax, f"{first_rv['year']}: {first_rv['planet']} — the first\n"
                    f"around a Sun-like star",
                    xy=(first_rv["year"], 6), xytext=(first_rv["year"] - 0.15, peak * 0.2),
                    relpos=(0, 0), ha="left", va="bottom", linespacing=1.35)

    # --- the missions behind the bars ------------------------------------------
    for year, label in ((2009, "Kepler\nlaunches"), (2018, "TESS\nlaunches")):
        top = totals[years == year][0]
        ps.annotate(ax, label, xy=(year, top + 8), xytext=(year, top + peak * 0.13),
                    relpos=(0.5, 0), ha="center", va="bottom", linespacing=1.2)

    partial = summary["partial_last_year"]
    through = datetime.fromisoformat(partial["through"]).strftime("%-d %b")
    ax.text(years[-1], totals[-1] + peak * 0.012, f"to\n{through}",
            ha="center", va="bottom", fontsize=7, color=ps.MUTED, linespacing=1.1)

    milestones = summary["cumulative_milestones"]
    running = " · ".join(f"{int(k):,} by {v}" for k, v in milestones.items()
                         if k in ("1000", "3000", "5000", "6000"))
    ax.legend(loc="upper left", title="Discovery method  (total)", title_fontsize=8,
              alignment="left", borderaxespad=0.2)
    ax.text(0.012, 0.66, f"Running total:\n{running}", transform=ax.transAxes,
            fontsize=8, color=ps.INK_SECONDARY, va="top", linespacing=1.4)

    ps.caption(fig, f"{_provenance(meta)}  Includes the {summary['controversial_included']} "
               f"planets the archive flags as controversial. {partial['year']} is partial "
               f"(to {partial['through']}). A discovery year is the year of the paper; "
               f"these are counts of detections, which says as much about the surveys "
               f"as about the planets.")
    return _save(fig, "01_discoveries_by_year")


# ===========================================================================
# 02  Radius vs period
# ===========================================================================

def figure_radius_period(sample, meta: dict) -> list:
    """Every planet with a measured radius, against its orbital period."""
    df = sample.frame
    xlim, ylim = (0.1, 2e5), (0.3, 35)
    shown = df[df["pl_orbper"].between(*xlim) & df["pl_rade"].between(*ylim)]
    n_off = len(df) - len(shown)

    fig, ax = plt.subplots(figsize=(10, 6.4))
    ax.scatter(shown["pl_orbper"], shown["pl_rade"], s=6, color=ps.C_DATA, alpha=0.45,
               linewidth=0, rasterized=True, zorder=3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.xaxis.set_major_locator(_decades(-1, 5))
    ax.yaxis.set_major_locator(FixedLocator([0.5, 1, 2, 4, 10, 20]))
    _plain_log(ax)
    ax.set_xlabel("Orbital period (days)")
    ax.set_ylabel("Planet radius (Earth radii)")
    frac_transit = (df["method"] == "Transit").mean()
    ax.set_title(f"{len(df):,} planets with a measured radius — "
                 f"{frac_transit:.0%} of them found by transits")

    _solar_system(ax, "period", "radius", offsets={
        "Mercury": (6, -9), "Venus": (-6, -10), "Earth": (6, 3), "Mars": (-6, -10),
        "Jupiter": (6, 2), "Saturn": (6, -9), "Uranus": (6, 2), "Neptune": (6, -9),
    })

    # Region labels sit in empty space, in ink; the data carry the colour.
    ps.annotate(ax, "Hot Jupiters:\ngiants orbiting\nin days — the first\nsurprise, in 1995",
                xy=(3.2, 14), xytext=(0.115, 17.5), relpos=(1, 0.5), ha="left", va="center",
                color=ps.INK, linespacing=1.3)
    ps.annotate(ax, "The 'Neptune desert':\nNeptune-sized planets\nare rare this close in",
                xy=(1.6, 5.2), xytext=(0.115, 5.6), relpos=(1, 0.5), ha="left", va="center",
                linespacing=1.3)
    ps.annotate(ax, "Super-Earths and sub-Neptunes:\nthe planets found most often",
                xy=(60, 2.2), xytext=(700, 2.2), relpos=(0, 0.5), ha="left", va="center",
                color=ps.INK, linespacing=1.3)
    ps.annotate(ax, "Earth-like orbits: a transit is rare\n(0.5% chance of alignment) and shallow\n"
                "(84 parts per million), so few are found.\nEmpty here means undetected, not absent.",
                xy=(365, 1.0), xytext=(1300, 0.78), relpos=(0, 1), ha="left", va="top",
                linespacing=1.35)

    ps.caption(fig, f"{_provenance(meta)}  Radii calculated by the archive from masses "
               f"(via an empirical mass–radius relation) are excluded, as are limits and "
               f"controversial planets; {n_off} planets fall outside the axes. Transit "
               f"probability and depth both fall steeply with orbital period, so this is a map "
               f"of what is findable, not of what exists.")
    return _save(fig, "02_radius_period")


# ===========================================================================
# 03  Mass vs distance, one panel per method
# ===========================================================================

#: Lens stars for the microlensing band: red dwarfs of 0.1-0.8 solar masses
#: halfway to a source star 8 kpc away in the Galactic bulge.
LENS_MASS_RANGE_MSUN = (0.1, 0.8)
LENS_DISTANCE_PC, SOURCE_DISTANCE_PC = 4000.0, 8000.0

#: Kept to ~26 characters a line so each box stays in the empty lower-right
#: corner of its panel, clear of the Solar System diamonds.
METHOD_NOTES = {
    "Transit": ("Sees the planet's shadow.\n"
                "Needs a nearly edge-on\n"
                "orbit (10% chance at\n"
                "0.05 AU, 0.5% at 1 AU)\n"
                "and several dips within\n"
                "a mission's few years."),
    "Radial velocity": ("Measures the star's\n"
                        "wobble, bigger for heavier\n"
                        "planets. Line: a 1 m/s\n"
                        "wobble of a Sun-like star;\n"
                        "Earth moves the Sun by\n"
                        "0.09 m/s. Masses are\n"
                        "mostly minimums (M sin i)."),
    "Microlensing": ("A star's gravity briefly\n"
                     "magnifies a background\n"
                     "star; a planet near its\n"
                     "Einstein ring (shaded,\n"
                     "{lo:.1f}–{hi:.1f} AU for red-dwarf\n"
                     "lenses) adds a blip.\n"
                     "Needs mass, not light."),
    "Direct imaging": ("Photographs the planet.\n"
                       "Needs it far from its\n"
                       "star's glare and still\n"
                       "hot from formation:\n"
                       "young, massive, wide."),
    "Other": ("Pulsar timing (the first\n"
              "planets, 1992), transit-\n"
              "timing variations,\n"
              "astrometry, eclipse\n"
              "timing and others."),
}

_NOTE_BOX = {"boxstyle": "square,pad=0.35", "facecolor": ps.SURFACE, "edgecolor": "none",
             "alpha": 0.9}
_NOTE_SIZE = 6.6


def figure_methods(sample, census: dict, analogs: dict, meta: dict) -> list:
    """Small multiples: where in mass and distance each method finds planets."""
    df = sample.frame
    xlim, ylim = (3e-3, 3e4), (0.01, 2e4)
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 8.2), sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.36, wspace=0.08)

    for ax, method in zip(axes.flat, METHOD_ORDER + ["Solar System"]):
        ax.scatter(df["pl_orbsmax"], df["pl_bmasse"], s=4, color=ps.C_CONTEXT, linewidth=0,
                   rasterized=True, zorder=2)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.xaxis.set_major_locator(_decades(-2, 4))
        ax.yaxis.set_major_locator(_decades(-2, 4))
        _plain_log(ax)
        ax.tick_params(labelsize=7)

        if method == "Solar System":
            _solar_system(ax, "a", "mass", size=34, offsets={
                "Mercury": (-6, 4), "Venus": (-6, 5), "Earth": (5, -10), "Mars": (5, -8),
                "Uranus": (-6, -10), "Neptune": (5, 4), "Saturn": (5, -10), "Jupiter": (-6, 5),
            })
            ax.set_title("The Solar System, on the same axes", fontsize=10, pad=17)
            ax.text(0.0, 1.02, f"grey in every panel: all {len(df):,} planets with a measured mass",
                    transform=ax.transAxes, fontsize=7.5, color=ps.MUTED, va="bottom")
            a = analogs["analogs"]
            n = {k: v["planets"] for k, v in a.items()}
            text = (f"Known planets within ×{analogs['factor']:.0f}\n"
                    f"of each one's mass\nand orbit size:\n"
                    f"Jupiter {n['Jupiter']} · Saturn {n['Saturn']}\n"
                    f"Earth {n['Earth']} · Venus {n['Venus']}\n"
                    f"Mercury {n['Mercury']} · Mars {n['Mars']}\n"
                    f"Uranus {n['Uranus']} · Neptune {n['Neptune']}")
            ax.text(0.97, 0.03, text, transform=ax.transAxes, fontsize=_NOTE_SIZE,
                    color=ps.INK_SECONDARY, ha="right", va="bottom", linespacing=1.3,
                    bbox=_NOTE_BOX, zorder=8)
            continue

        sel = df["method"] == method
        ax.scatter(df.loc[sel, "pl_orbsmax"], df.loc[sel, "pl_bmasse"], s=9,
                   color=ps.METHOD_COLORS[method], alpha=0.85, linewidth=0, zorder=4)
        _solar_system(ax, "a", "mass", labels=False, size=18)
        ax.set_title(method, fontsize=10, pad=17)
        ax.text(0.0, 1.02, f"{sel.sum():,} with a measured mass, of "
                f"{census['by_method'][method]:,} found",
                transform=ax.transAxes, fontsize=7.5, color=ps.MUTED, va="bottom")

        note = METHOD_NOTES[method]
        if method == "Microlensing":
            lo, hi = (float(pp.einstein_radius_au(m, LENS_DISTANCE_PC, SOURCE_DISTANCE_PC))
                      for m in LENS_MASS_RANGE_MSUN)
            ax.axvspan(lo, hi, color=ps.METHOD_COLORS[method], alpha=0.14, linewidth=0, zorder=1)
            note = note.format(lo=lo, hi=hi)
        if method == "Radial velocity":
            a = np.logspace(np.log10(xlim[0]), np.log10(xlim[1]), 200)
            m_line = pp.rv_minimum_mass_mearth(1.0, a, 1.0)
            ax.plot(a, m_line, color=ps.INK_SECONDARY, linewidth=1.0, zorder=5)
            ax.text(2000, float(pp.rv_minimum_mass_mearth(1.0, 2000, 1.0)) * 1.6, "1 m/s",
                    fontsize=7, color=ps.INK_SECONDARY, ha="center", va="bottom")
        if method == "Other":
            psr = df[sel & df["pl_name"].str.startswith("PSR B1257+12")]
            if len(psr):
                smallest = psr.loc[psr["pl_bmasse"].idxmin()]
                ps.annotate(ax, f"{smallest['pl_name']}\n({smallest['pl_bmasse']:.2f} M⊕)",
                            xy=(smallest["pl_orbsmax"], smallest["pl_bmasse"]),
                            xytext=(0.0045, 0.018), relpos=(1, 0.5), ha="left", va="center",
                            linespacing=1.25, fontsize=6.8)
        ax.text(0.97, 0.03, note, transform=ax.transAxes, fontsize=_NOTE_SIZE,
                color=ps.INK_SECONDARY, ha="right", va="bottom", linespacing=1.3,
                bbox=_NOTE_BOX, zorder=8)

    for ax in axes[1, :]:
        ax.set_xlabel("Orbital distance (semi-major axis, AU)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Planet mass (Earth masses)")
    fig.suptitle("Why each method finds what it finds", x=0.125, ha="left", y=0.975)

    # The Earth matches get a sentence only when the data support it.
    earth = analogs["analogs"]["Earth"]
    hosts = [m for m in earth["host_masses_msun"] if m == m]  # drop NaN
    earth_note = ""
    if (earth["planets"] and set(earth["by_method"]) == {"Microlensing"}
            and hosts and max(hosts) < 0.2):
        earth_note = (f" The {earth['planets']} Earth matches were found by microlensing, "
                      f"around stars of {min(hosts):.2f}–{max(hosts):.2f} solar masses: "
                      f"Earth-like in mass and orbit size, but far colder.")
    lo_m, hi_m = LENS_MASS_RANGE_MSUN
    ps.caption(fig, f"{_provenance(meta)}  Masses the archive derived from radii via an "
               f"empirical relation are excluded, as are mass limits and controversial "
               f"planets. Radial-velocity masses are mostly minimum masses (M sin i); "
               f"microlensing separations are projected distances. Einstein radius for a "
               f"{lo_m}–{hi_m} solar-mass lens {LENS_DISTANCE_PC / 1000:.0f} kpc away, in front of "
               f"a source {SOURCE_DISTANCE_PC / 1000:.0f} kpc away. A match is within a factor of "
               f"{analogs['factor']:.0f} in both mass and semi-major axis.{earth_note}")
    return _save(fig, "03_methods_mass_distance")


# ===========================================================================
# 04  The radius valley
# ===========================================================================

def figure_radius_valley(sample, valley: dict, art: dict, meta: dict) -> list:
    df = sample.frame
    fig = plt.figure(figsize=(12, 9.8))
    gs = GridSpec(2, 5, figure=fig, height_ratios=[1, 1.05], hspace=0.62, wspace=0.9)
    ax_a = fig.add_subplot(gs[0, :3])
    ax_b = fig.add_subplot(gs[0, 3:])
    ax_c = fig.add_subplot(gs[1, :])
    radius_ticks = FixedLocator([0.8, 1, 1.5, 2, 3, 4])

    # --- (a) the histogram ---------------------------------------------------------
    width_dex = 0.035
    edges = 10 ** np.arange(np.log10(0.7), np.log10(5.0) + 1e-9, width_dex)
    counts, _ = np.histogram(df["pl_rade"], edges)
    centres = np.sqrt(edges[1:] * edges[:-1])
    ax_a.bar(edges[:-1], counts, width=np.diff(edges), align="edge", color=ps.C_DATA,
             alpha=0.85, edgecolor=ps.SURFACE, linewidth=0.7, zorder=3,
             label="planets per bin")
    ax_a.errorbar(centres, counts, yerr=np.sqrt(counts), fmt="none", ecolor=ps.INK_SECONDARY,
                  elinewidth=0.8, capsize=0, zorder=4, label="±√N counting error")
    scale = len(df) * width_dex
    ax_a.plot(10 ** art["log_grid"], art["density"] * scale, color=ps.C_FIT, linewidth=1.8,
              zorder=5, label="smoothed (kernel density)")
    lo, _mid, hi = valley["valley_rearth_16_50_84"]
    top = counts.max()
    y_top = top * 1.55
    line_top = y_top * 0.70  # the valley marker stops below its own label
    ax_a.fill_between([lo, hi], 0, line_top, color=ps.C_FIT, alpha=0.18, linewidth=0, zorder=2)
    ax_a.plot([valley["valley_rearth"]] * 2, [0, line_top], color=ps.C_FIT, linewidth=1.0, zorder=2)

    ax_a.set_xscale("log")
    ax_a.set_xlim(0.7, 5.0)
    ax_a.xaxis.set_major_locator(radius_ticks)
    _plain_log(ax_a, "x")
    ax_a.set_ylim(0, y_top)
    ax_a.set_xlabel("Planet radius (Earth radii)")
    ax_a.set_ylabel(f"Planets per bin (bins {10**width_dex - 1:.0%} wide in radius)")
    ax_a.set_title(f"(a) A gap in {len(df):,} precisely measured planets")
    ax_a.grid(axis="x", visible=False)
    # Right-aligned against the marker, clear of the sub-Neptune label.
    ax_a.text(valley["valley_rearth"] * 1.012, line_top * 0.97,
              f"valley {valley['valley_rearth']:.2f} R⊕\n68% range {lo:.2f}–{hi:.2f}",
              ha="right", va="bottom", fontsize=8, color=ps.INK, linespacing=1.3)
    def bar_top(r: float) -> float:
        """Tallest bar within ~10% of radius r, so a label clears the bars under it."""
        return float(counts[np.abs(np.log10(centres / r)) < 0.045].max())

    for key, name in (("peak_super_earth_rearth", "super-Earths"),
                      ("peak_sub_neptune_rearth", "sub-Neptunes")):
        ax_a.text(valley[key], bar_top(valley[key]) + 0.06 * top,
                  f"{name}\npeak {valley[key]:.2f} R⊕",
                  ha="center", va="bottom", fontsize=8, color=ps.INK_SECONDARY, linespacing=1.25)
    ax_a.plot([1.5, 2.0], [y_top * 0.93] * 2, color=ps.C_PUBLISHED, linewidth=2.2,
              solid_capstyle="butt", zorder=5)
    ax_a.text(2.03, y_top * 0.93, "Fulton et al. 2017: gap at 1.5–2.0 R⊕", ha="left",
              va="center", fontsize=7.5, color=ps.INK_SECONDARY)
    ax_a.legend(loc="upper left")

    # --- (b) precision tiers ----------------------------------------------------------
    tiers = valley["precision_tiers"]
    ramp = [ps.COMPOSITION_RAMP[0], ps.COMPOSITION_RAMP[2], ps.COMPOSITION_RAMP[4]]
    x = 10 ** art["log_grid"]
    for tier, curve, color in zip(tiers, art["tier_curves"], ramp):
        verdict = (f"valley floor at {tier['contrast']:.0%} of the peaks"
                   if tier["valley_found"] else "no valley")
        ax_b.plot(x, curve, color=color, linewidth=2.0,
                  label=f"radii to {tier['label']} ({tier['planets_0p8_to_4p5']} planets)\n{verdict}")
    ax_b.set_xscale("log")
    ax_b.set_xlim(0.8, 4.5)
    ax_b.xaxis.set_major_locator(radius_ticks)
    _plain_log(ax_b, "x")
    ax_b.set_ylim(0, max(c.max() for c in art["tier_curves"]) * 1.6)
    ax_b.set_xlabel("Planet radius (Earth radii)")
    ax_b.set_ylabel("Relative number of planets")
    ax_b.set_title("(b) Blurry radii hide it")
    ax_b.grid(axis="x", visible=False)
    ax_b.legend(loc="upper left", fontsize=7.2, handlelength=1.6, labelspacing=0.6,
                title="The same cuts, split into thirds by\nhow precisely each radius is known",
                title_fontsize=7.2, alignment="left")

    # --- (c) the tilt -------------------------------------------------------------------
    tilt = valley["tilt"]
    ax_c.scatter(df["pl_orbper"], df["pl_rade"], s=9, color=ps.C_DATA, alpha=0.5, linewidth=0,
                 rasterized=True, zorder=3)
    p = np.logspace(np.log10(0.5), np.log10(100), 200)
    lines = np.array([pp.valley_radius(p, m, a) for m, a in zip(art["boot_slopes"], art["boot_intercepts"])])
    band_lo, band_hi = np.percentile(lines, [16, 84], axis=0)
    ax_c.fill_between(p, band_lo, band_hi, color=ps.C_FIT, alpha=0.18, linewidth=0, zorder=2)
    s16, _s50, s84 = tilt["slope_16_50_84"]
    ax_c.plot(p, pp.valley_radius(p, tilt["slope"], tilt["intercept"]), color=ps.C_FIT,
              linewidth=2.0, zorder=5,
              label=f"this project: slope {tilt['slope']:+.2f} (68%: {s16:+.2f} to {s84:+.2f})")
    ve = valley["van_eylen_2018"]
    ax_c.plot(p, pp.valley_radius(p, VAN_EYLEN_SLOPE, VAN_EYLEN_INTERCEPT), color=ps.C_PUBLISHED,
              linewidth=2.0, zorder=5,
              label=f"Van Eylen et al. 2018: slope {ve['slope']:+.2f} (117 planets)")
    ax_c.set_xscale("log")
    ax_c.set_yscale("log")
    ax_c.set_xlim(0.5, 100)
    ax_c.set_ylim(0.8, 4.5)
    ax_c.yaxis.set_major_locator(radius_ticks)
    ax_c.xaxis.set_major_locator(FixedLocator([0.5, 1, 2, 5, 10, 20, 50, 100]))
    _plain_log(ax_c)
    ax_c.set_xlabel("Orbital period (days)")
    ax_c.set_ylabel("Planet radius (Earth radii)")
    # Title, explanation and legend all sit above the axes, so nothing covers data.
    ax_c.set_title("(c) The valley tilts: smaller radii on longer orbits", pad=50)
    ax_c.text(0.0, 1.02,
              f"The slope is negative in {tilt['fraction_negative']:.0%} of bootstrap resamples — "
              f"what atmospheric\nescape driven by the star's heat predicts. Rocky planets "
              f"forming late,\nwithout gas, would tilt it the other way (Van Eylen et al. 2018).",
              transform=ax_c.transAxes, fontsize=7.8, color=ps.INK_SECONDARY, va="bottom",
              linespacing=1.4)
    ax_c.legend(loc="lower right", bbox_to_anchor=(1.0, 1.005), fontsize=7.8,
                title="valley line: log R = slope × log P + const", title_fontsize=7.8,
                alignment="right", borderaxespad=0.0)

    ps.caption(fig, f"{_provenance(meta)}  Transiting planets with P < 100 d around FGK "
               f"main-sequence stars (Teff 4700–6500 K, log g > 3.9), radii measured by the "
               f"archive's source papers to better than 10% for (a) and (c); (b) uses the same "
               f"cuts at every precision. Counts of detected planets, not occurrence rates: "
               f"small planets are harder to detect, so the true super-Earth peak is taller "
               f"relative to the sub-Neptune one. Shaded bands are bootstrap 68% ranges.")
    return _save(fig, "04_radius_valley")


# ===========================================================================
# 05  Mass vs radius
# ===========================================================================

def figure_mass_radius(sample_frame, comp: dict, meta: dict) -> list:
    df = sample_frame
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 6.2),
                                     gridspec_kw={"width_ratios": [1, 1.12], "wspace": 0.2})
    models = {k: load_model(k) for k in ZENG_MODELS}
    labels = {k: v[1] for k, v in ZENG_MODELS.items()}

    # --- (a) the whole range ---------------------------------------------------------------
    xa, ya = (0.1, 2e4), (0.3, 30)
    ax_a.scatter(df["pl_bmasse"], df["pl_rade"], s=7, color=ps.INK_SECONDARY, alpha=0.55,
                 linewidth=0, rasterized=True, zorder=3)
    for key, color in (("earthlike", ps.COMPOSITION_RAMP[1]), ("hydrogen", ps.C_HYDROGEN)):
        m, r = models[key]
        ax_a.plot(m, r, color=color, linewidth=1.8, zorder=4, label=labels[key])
    ax_a.axvline(pp.DEUTERIUM_LIMIT_MEARTH, color=ps.AXIS, linewidth=0.9, zorder=1)
    ax_a.text(pp.DEUTERIUM_LIMIT_MEARTH * 1.12, 0.36, "13 Jupiter masses:\nbrown dwarfs beyond",
              fontsize=7, color=ps.MUTED, va="bottom", linespacing=1.2)
    _solar_system(ax_a, "mass", "radius", offsets={
        "Mercury": (5, -8), "Venus": (-5, 5), "Earth": (5, -9), "Mars": (5, -3),
        "Jupiter": (6, -10), "Saturn": (-6, 5), "Uranus": (-6, 3), "Neptune": (6, -8),
    })
    ax_a.set_xscale("log")
    ax_a.set_yscale("log")
    ax_a.set_xlim(*xa)
    ax_a.set_ylim(*ya)
    ax_a.xaxis.set_major_locator(_decades(-1, 4))
    ax_a.yaxis.set_major_locator(FixedLocator([0.5, 1, 2, 4, 10, 20]))
    _plain_log(ax_a)
    ax_a.set_xlabel("Planet mass (Earth masses)")
    ax_a.set_ylabel("Planet radius (Earth radii)")
    ax_a.set_title(f"(a) {len(df):,} planets weighed and measured")
    ax_a.add_patch(Rectangle((0.4, 0.5), 30 - 0.4, 4.2 - 0.5, fill=False, edgecolor=ps.MUTED,
                             linewidth=0.8, zorder=2))
    ax_a.text(0.42, 4.35, "(b)", fontsize=7.5, color=ps.MUTED)
    ax_a.text(900, 23.5, "Inflated hot Jupiters: bigger than\neven a ball of cold hydrogen",
              fontsize=7.5, color=ps.INK_SECONDARY, ha="center", va="bottom", linespacing=1.3)
    ax_a.text(1700, 7.3, "Giants stop growing: extra mass\ncompresses them instead",
              fontsize=7.5, color=ps.INK_SECONDARY, ha="center", va="top", linespacing=1.3)
    ax_a.legend(loc="upper left", title="Models (Zeng et al. 2019)", title_fontsize=8,
                alignment="left")

    # --- (b) small planets and what they could be made of ---------------------------
    xb, yb = (0.4, 30), (0.5, 4.2)
    small = df[df["pl_bmasse"].between(*xb) & df["pl_rade"].between(*yb)]
    ax_b.errorbar(small["pl_bmasse"], small["pl_rade"],
                  xerr=small["pl_bmasse"] * small["mass_relerr"],
                  yerr=small["pl_rade"] * small["radius_relerr"],
                  fmt="none", ecolor=ps.C_CONTEXT, elinewidth=0.7, zorder=2)
    ax_b.scatter(small["pl_bmasse"], small["pl_rade"], s=12, color=ps.INK_SECONDARY,
                 edgecolor=ps.SURFACE, linewidth=0.5, zorder=4)

    curves = [("iron", ps.COMPOSITION_RAMP[0]), ("earthlike", ps.COMPOSITION_RAMP[1]),
              ("rock", ps.COMPOSITION_RAMP[2]), ("water50", ps.COMPOSITION_RAMP[3]),
              ("water100", ps.COMPOSITION_RAMP[4]), ("h2_1pct", ps.C_HYDROGEN)]
    for key, color in curves:
        m, r = models[key]
        keep = (m >= xb[0] * 0.9) & (m <= xb[1] * 1.1)
        ax_b.plot(m[keep], r[keep], color=color, linewidth=1.8, zorder=3, label=labels[key])

    _solar_system(ax_b, "mass", "radius", names=("Venus", "Earth", "Uranus", "Neptune"),
                  offsets={"Venus": (-6, 5), "Earth": (-6, 5), "Uranus": (-6, -10), "Neptune": (6, -8)})
    for name in MR_HIGHLIGHT:
        row = df[df["pl_name"] == name]
        if not len(row):
            continue
        m, r = float(row["pl_bmasse"].iloc[0]), float(row["pl_rade"].iloc[0])
        ax_b.scatter([m], [r], s=30, color=ps.INK, edgecolor=ps.SURFACE, linewidth=1.0, zorder=6)
    offsets = {"Kepler-10 b": (-60, 32), "55 Cnc e": (40, -34), "K2-18 b": (-60, 22),
               "GJ 1214 b": (-50, 22), "TRAPPIST-1 e": (-12, -34), "LHS 1140 b": (22, -44)}
    for name in MR_HIGHLIGHT:
        row = df[df["pl_name"] == name]
        if not len(row):
            continue
        m, r = float(row["pl_bmasse"].iloc[0]), float(row["pl_rade"].iloc[0])
        suffix = " — found in A1" if name == "Kepler-10 b" else ""
        dx, dy = offsets.get(name, (10, 10))
        ps.annotate(ax_b, name + suffix, xy=(m, r), xytext=(dx, dy), textcoords="offset points",
                    color=ps.INK, ha="left" if dx >= 0 else "right")

    ax_b.set_xscale("log")
    ax_b.set_xlim(*xb)
    ax_b.set_ylim(*yb)
    ax_b.xaxis.set_major_locator(FixedLocator([0.5, 1, 2, 5, 10, 20]))
    _plain_log(ax_b, "x")
    ax_b.set_xlabel("Planet mass (Earth masses)")
    ax_b.set_ylabel("Planet radius (Earth radii)")
    ax_b.set_title("(b) Small planets against composition models")
    # Framed so the hydrogen curve passes behind the legend, not through its text.
    ax_b.legend(loc="upper left", title=f"Models (water and hydrogen at {ZENG_TEMPERATURE_K} K)",
                title_fontsize=7.5, fontsize=7.5, alignment="left", frameon=True,
                facecolor=ps.SURFACE, edgecolor="none", framealpha=0.95)

    # What 1% hydrogen does to a rocky core: an arrow from bare to enveloped.
    env = comp["envelope_illustration"]
    core, r0, r1 = env["core_mass_mearth"], env["bare_rearth"], env["with_1pct_h2_rearth"]
    ax_b.annotate("", xy=(core, r1), xytext=(core, r0), zorder=7,
                  arrowprops={"arrowstyle": "-|>", "color": ps.INK, "linewidth": 1.4,
                              "shrinkA": 4, "shrinkB": 4, "mutation_scale": 10})
    for r in (r0, r1):
        ax_b.scatter([core], [r], s=40, facecolor=ps.SURFACE, edgecolor=ps.INK,
                     linewidth=1.3, zorder=8)
    ps.annotate(ax_b, f"Add 1% hydrogen by mass to a {core:.0f} M⊕ rocky core\n"
                f"and its radius grows from {r0:.2f} to {r1:.2f} R⊕ —\n"
                f"the two peaks of the radius valley",
                xy=(core, (r0 + r1) / 2), xytext=(29, 0.62), relpos=(0.3, 1),
                ha="right", va="bottom", color=ps.INK, fontsize=7.5, linespacing=1.35)

    ax_b.text(0.015, 0.70,
              f"Below the radius valley:\n"
              f"{comp['below_consistent_with_rock']:.0%} of {comp['below_valley']} fit bare rock.\n"
              f"Above it, of {comp['above_valley']}:\n"
              f"{comp['above_too_big_for_rock']:.0%} are too big for rock,\n"
              f"{comp['above_too_big_for_water']:.0%} too big even for\n"
              f"pure water — they need\nhydrogen and helium.",
              transform=ax_b.transAxes, ha="left", va="top", fontsize=7.3,
              color=ps.INK_SECONDARY, linespacing=1.4, zorder=8,
              bbox={"boxstyle": "square,pad=0.4", "facecolor": ps.SURFACE, "edgecolor": ps.GRID})

    ps.caption(fig, f"{_provenance(meta)}  Transiting planets with a measured true mass "
               f"(uncertainty < 25%) and radius (< 10%); archive-calculated values, limits and "
               f"controversial planets excluded. Model curves: Zeng et al. 2019 (PNAS 116, 9723); "
               f"the hydrogen curve turns up at low mass because a light core's weak gravity lets "
               f"a hot envelope swell. 'Too big' means more than "
               f"{comp['excess_threshold_sigma']:.0f}σ above a curve, using the pure-water model "
               f"nearest each planet's own equilibrium temperature. Mass and radius alone cannot "
               f"tell a water world from a rocky core with a thin hydrogen envelope.")
    return _save(fig, "05_mass_radius")


# ===========================================================================
# 06  The habitable zone
# ===========================================================================

def figure_habitable_zone(frame, hz: dict, meta: dict) -> list:
    df = frame
    xlim, ylim = (5.0, 0.1), (2150, 7400)
    t_lo, t_hi = pp.KOPPARAPU_TEFF_RANGE
    # Legend and explanation live in a right-hand margin: this plot has no
    # empty corner to spare.
    fig, ax = plt.subplots(figsize=(10.2, 7.0))

    # The fits stop at 2,600 K. Planets of cooler stars are classified with the
    # 2,600 K limits, so the bands continue flat below it -- otherwise TRAPPIST-1
    # would appear outside a zone it is counted in. The seam is marked.
    teff = np.linspace(ylim[0], t_hi, 400)
    s = {k: pp.hz_flux_limit(k, teff) for k in pp.KOPPARAPU_2014}
    ax.fill_betweenx(teff, s["early_mars"], s["recent_venus"], color=ps.C_HZ, alpha=0.10,
                     linewidth=0, zorder=1)
    ax.fill_betweenx(teff, s["maximum_greenhouse"], s["runaway_greenhouse"], color=ps.C_HZ,
                     alpha=0.16, linewidth=0, zorder=1)
    for key in pp.KOPPARAPU_2014:
        ax.plot(s[key], teff, color=ps.C_HZ, linewidth=0.8, alpha=0.7, zorder=2)
    ax.plot([s["recent_venus"][0], s["early_mars"][0]], [t_lo, t_lo], color=ps.C_HZ,
            linewidth=0.8, alpha=0.7, zorder=2)
    ax.text(s["early_mars"][0] * 0.93, t_lo, f"fits end at {t_lo:,.0f} K;\nheld constant below",
            fontsize=6.8, color=ps.INK_SECONDARY, ha="left", va="center", linespacing=1.2)

    # Edge names at the top, pushed outwards so neighbouring pairs never touch.
    top = -1
    for key, name, ha, nudge in (("recent_venus", "recent\nVenus", "right", 1.03),
                                 ("runaway_greenhouse", "runaway\ngreenhouse", "left", 0.97),
                                 ("maximum_greenhouse", "maximum\ngreenhouse", "right", 1.03),
                                 ("early_mars", "early\nMars", "left", 0.97)):
        ax.text(s[key][top] * nudge, t_hi - 20, name, fontsize=6.8, color=ps.INK_SECONDARY,
                ha=ha, va="top", linespacing=1.1)
    mid = np.sqrt(s["maximum_greenhouse"][top] * s["runaway_greenhouse"][top])
    ax.text(mid, t_hi + 90, "conservative habitable zone", fontsize=8, color=ps.INK,
            ha="center", va="bottom", fontweight="bold")

    shown = df[df["pl_insol"].between(xlim[1], xlim[0]) & df["st_teff"].between(*ylim)]
    out = shown[shown["hz"] == ""]
    ax.scatter(out["pl_insol"], out["st_teff"], s=7, color=ps.C_CONTEXT, linewidth=0, zorder=3)
    inz = shown[shown["hz"] != ""]
    larger = inz[inz["size_class"] == "larger"]
    noradius = inz[inz["size_class"] == "no radius"]
    rocky = inz[inz["size_class"] == "rocky-size"]
    ax.scatter(larger["pl_insol"], larger["st_teff"], s=20, color=ps.MUTED, edgecolor=ps.SURFACE,
               linewidth=0.8, zorder=4)
    ax.scatter(noradius["pl_insol"], noradius["st_teff"], s=20, facecolor="none",
               edgecolor=ps.INK_SECONDARY, linewidth=0.8, zorder=4)
    ax.scatter(rocky["pl_insol"], rocky["st_teff"], s=46, color=ps.C_DATA, edgecolor=ps.SURFACE,
               linewidth=1.2, zorder=6)

    # Solar System anchors: the Sun's temperature and each planet's 1/a^2.
    for name in ("Venus", "Earth", "Mars"):
        s_pl = pp.SOLAR_SYSTEM[name]["insolation"]
        ax.scatter([s_pl], [pp.T_SUN], marker="D", s=34, color=ps.C_SOLAR, edgecolor=ps.SURFACE,
                   linewidth=0.9, zorder=7)
        ax.annotate(name, (s_pl, pp.T_SUN), xytext=(0, 9), textcoords="offset points",
                    fontsize=8, color=ps.INK, ha="center", zorder=8)

    # TRAPPIST-1's planets share one host temperature: letter them, name them once.
    trappist = inz[inz["hostname"] == "TRAPPIST-1"].sort_values("pl_insol", ascending=False)
    for _, r in trappist.iterrows():
        ax.annotate(r["pl_name"].split()[-1], (r["pl_insol"], r["st_teff"]), xytext=(0, -7),
                    textcoords="offset points", fontsize=8, color=ps.INK, ha="center",
                    va="top", zorder=8)
    if len(trappist):
        ax.text(np.sqrt(trappist["pl_insol"].max() * trappist["pl_insol"].min()),
                trappist["st_teff"].iloc[0] - 240, "TRAPPIST-1's planets (2,566 K)", fontsize=8,
                color=ps.INK, ha="center", va="top")

    offsets = {
        "Proxima Cen b": (40, -14), "TOI-700 d": (34, 16), "LHS 1140 b": (28, 12),
        "K2-18 b": (0, 52), "Teegarden's Star b": (-28, -18), "Kepler-442 b": (32, 14),
        "Kepler-62 f": (30, 22),
    }
    for name in HZ_HIGHLIGHT:
        row = df[df["pl_name"] == name]
        if not len(row) or name.startswith("TRAPPIST-1"):
            continue
        r = row.iloc[0]
        if not (xlim[1] <= r["pl_insol"] <= xlim[0]):
            continue
        dx, dy = offsets.get(name, (14, 14))
        if dx == 0:  # straight above: centred text, leader from its bottom
            ha, rel = "center", (0.5, 0)
        else:
            ha, rel = ("left", (0, 0.5)) if dx > 0 else ("right", (1, 0.5))
        ps.annotate(ax, name, xy=(r["pl_insol"], r["st_teff"]), xytext=(dx, dy),
                    textcoords="offset points", color=ps.INK, relpos=rel, ha=ha,
                    va="bottom" if dx == 0 else "center")

    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.xaxis.set_major_locator(FixedLocator([0.1, 0.2, 0.5, 1, 2, 5]))
    _plain_log(ax, "x")
    ax.set_xlabel("Starlight received (Earth = 1)   ← hotter · colder →")
    ax.set_ylabel("Host star temperature (K)")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.set_title(f"{hz['in_optimistic_hz']} planets in the habitable zone; "
                 f"{hz['rocky_size_in_hz']} small enough to plausibly be rocky", pad=24)

    handles = [
        Line2D([], [], marker="o", linestyle="", markersize=7, color=ps.C_DATA,
               markeredgecolor=ps.SURFACE,
               label=f"in the zone, radius ≤ {hz['rocky_size_threshold_rearth']} R⊕"),
        Line2D([], [], marker="o", linestyle="", markersize=5.5, color=ps.MUTED,
               markeredgecolor=ps.SURFACE, label="in the zone, larger"),
        Line2D([], [], marker="o", linestyle="", markersize=5.5, markerfacecolor="none",
               markeredgecolor=ps.INK_SECONDARY, label="in the zone, radius unknown (mass only)"),
        Line2D([], [], marker="o", linestyle="", markersize=4, color=ps.C_CONTEXT,
               label="outside the zone"),
        Line2D([], [], marker="D", linestyle="", markersize=5, color=ps.C_SOLAR,
               markeredgecolor=ps.SURFACE, label="Solar System (around a 5,772 K Sun)"),
        Patch(facecolor=ps.C_HZ, alpha=0.26, label="conservative zone"),
        Patch(facecolor=ps.C_HZ, alpha=0.10, label="optimistic extension"),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=7.5,
              borderaxespad=0.0)

    ax.text(1.02, 0.0,
            f"{hz['m_dwarf_fraction_rocky_hz']:.0%} of the small planets in\n"
            f"the zone orbit red dwarfs\n(cooler than 3,900 K).\n\n"
            f"A red dwarf's habitable zone is\n"
            f"so close that an orbit takes days\n"
            f"to weeks, so transits recur often\n"
            f"and are deep. Around a Sun-like\n"
            f"star the orbit takes a year, and\n"
            f"an Earth's transit dims it by only\n"
            f"84 parts per million.\n\n"
            f"Selection, not nature: this says\n"
            f"where planets are easy to find.",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=7.5,
            color=ps.INK_SECONDARY, linespacing=1.4)

    hidden = " and ".join(CONTROVERSIAL_FAMOUS)
    ps.caption(fig, f"{_provenance(meta)}  Habitable-zone limits: Kopparapu et al. 2013, "
                    f"2014 (Earth-mass planet), fitted for {t_lo:,.0f}–{t_hi:,.0f} K. Insolation as "
                    f"listed by the archive, or computed by it as L/a². In the zone means starlight "
                    f"compatible with liquid surface water given a suitable atmosphere — not that a "
                    f"planet is habitable. {hidden} are absent because the archive flags them as "
                    f"controversial.")
    return _save(fig, "06_habitable_zone")


def make_all(samples: dict, results: dict, artefacts: dict, meta: dict) -> list:
    written = []
    written += figure_discoveries(results["census"], artefacts["census"], meta)
    written += figure_radius_period(samples["radius_period"], meta)
    written += figure_methods(samples["mass_separation"], results["census"],
                              results["solar_analogs"], meta)
    written += figure_radius_valley(samples["valley"], results["valley"], artefacts["valley"], meta)
    written += figure_mass_radius(artefacts["composition"]["frame"], results["composition"], meta)
    written += figure_habitable_zone(artefacts["habitable_zone"]["frame"], results["habitable_zone"], meta)
    return written
