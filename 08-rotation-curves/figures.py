"""Figure generation.

Every figure is drawn from the analysis results and artefacts, never by
re-running the analysis, so what is plotted is exactly what was measured.

One colour rule runs through all of them, so a colour means the same thing in
every panel: **ink** is what was observed; **blue** is what the visible matter
predicts; **orange** is the unseen mass and our fits to it; **violet** is a
published result we compare against; aqua, yellow and magenta are the visible
components -- gas, stellar disc, bulge. The palette is colourblind-safe
(validated with a checker, not by eye), and the lighter hues always carry text
labels because they fall below 3:1 contrast on the page.

Figure text writes the stellar mass-to-light ratio as "M/L". The usual symbol,
the Greek capital upsilon, is indistinguishable from a Latin Y in this font.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator, MultipleLocator, NullFormatter

import plotstyle as ps
import rotation_physics as rp
from config import FIGURE_DIR, UPSILON_BULGE, UPSILON_DISK, UPSILON_PLAUSIBLE

ps.use_style()

PROVENANCE = ("Data: SPARC (Lelli, McGaugh & Schombert 2016, AJ 152, 157) — Spitzer 3.6 µm "
              "photometry of the stars, 21-cm radio maps of the gas, and HI/Hα rotation curves "
              "compiled from the literature.")

AU_KM = 1.495978707e8
LW_THIN = 1.3
LW_MAIN = 2.2
DASHED = (0, (5, 2.5))
DOTTED = (0, (1, 1.8))


# --------------------------------------------------------------------------
# Building blocks
# --------------------------------------------------------------------------

def _title(ax, title: str, subtitle: str | None = None) -> None:
    """Bold title with an optional regular-weight subtitle under it.

    The title's padding grows with the number of subtitle lines so the two
    never overlap.
    """
    lines = subtitle.count("\n") + 1 if subtitle else 0
    ax.set_title(title, pad=8 + 10.8 * lines)
    if subtitle:
        ax.annotate(subtitle, xy=(0, 1), xycoords="axes fraction", xytext=(0, 5),
                    textcoords="offset points", ha="left", va="bottom", fontsize=8.2,
                    color=ps.INK_SECONDARY, linespacing=1.3)


def _label(ax, text, xy, xytext, *, ha="left", va="center", leader=True, fontsize=8,
           textcoords="data"):
    """A direct label in secondary ink, tied to its point by a hairline if moved."""
    ax.annotate(text, xy=xy, xytext=xytext, textcoords=textcoords, ha=ha, va=va,
                fontsize=fontsize, color=ps.INK_SECONDARY,
                arrowprops=({"arrowstyle": "-", "color": ps.MUTED, "linewidth": 0.6,
                             "shrinkA": 1, "shrinkB": 3} if leader else None))


def _observed(ax, r, v, e, *, label="observed rotation", ms=4.6, zorder=6):
    """Measured speeds: ink points with a surface ring, recessive error bars."""
    return ax.errorbar(r, v, yerr=e, fmt="o", color=ps.C_OBSERVED, ecolor=ps.C_OBSERVED_ERR,
                       elinewidth=0.9, capsize=0, markersize=ms, markeredgecolor=ps.SURFACE,
                       markeredgewidth=0.9, zorder=zorder, label=label)


def _end_labels(ax, items, x_text, *, min_gap):
    """Direct labels at line ends, nudged apart with leader lines if they collide.

    ``items`` are ``(x_end, y_end, text)``. Labels are sorted by height and
    pushed apart to at least ``min_gap`` (data units); each is tied back to its
    line by a hairline.
    """
    items = sorted(items, key=lambda t: t[1])
    ys = [t[1] for t in items]
    for i in range(1, len(ys)):
        ys[i] = max(ys[i], ys[i - 1] + min_gap)
    for (x0, y0, text), y in zip(items, ys):
        _label(ax, text, (x0, y0), (x_text, y))


def _gap(ax, r, v_bar, v_obs, *, label="the gap: unseen mass"):
    """Wash between what the visible matter predicts and what is observed."""
    return ax.fill_between(r, v_bar, v_obs, where=v_obs >= v_bar, interpolate=True,
                           color=ps.C_DARK, alpha=0.12, linewidth=0, zorder=1, label=label)


def _rotation_panel(ax, art, *, components=True, gap=True):
    """The standard picture: observed speeds against the visible matter's prediction."""
    g, comp = art["galaxy"], art["components"]
    c = g.curve
    handles = []
    if gap:
        handles.append(_gap(ax, c.r, comp["baryons"], c.v_obs))
    if components:
        handles += ax.plot(c.r, comp["gas"], color=ps.C_GAS, lw=LW_THIN, label="gas", zorder=3)
        handles += ax.plot(c.r, comp["disk"], color=ps.C_DISK, lw=LW_THIN,
                           label=f"stellar disc (M/L = {UPSILON_DISK:g})", zorder=3)
        if g.has_bulge:
            handles += ax.plot(c.r, comp["bulge"], color=ps.C_BULGE, lw=LW_THIN,
                               label=f"bulge (M/L = {UPSILON_BULGE:g})", zorder=3)
    handles += ax.plot(c.r, comp["baryons"], color=ps.C_VISIBLE, lw=LW_MAIN,
                       label="all visible matter", zorder=4, solid_capstyle="round")
    handles.append(_observed(ax, c.r, c.v_obs, c.e_v_obs))
    ax.set_ylim(0, float(np.max(c.v_obs + c.e_v_obs)) * 1.18)
    ax.set_xlim(0, c.r[-1] * 1.04)
    return handles


def _sci(value: float, digits: int = 1) -> str:
    """1.23e10 -> '1.2 × 10¹⁰', for annotations."""
    exponent = int(np.floor(np.log10(abs(value))))
    sup = str(exponent).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{value / 10**exponent:.{digits}f} × 10{sup}"


def _log_ticks(axis):
    axis.set_major_locator(LogLocator(base=10))
    axis.set_minor_formatter(NullFormatter())


# --------------------------------------------------------------------------
# 1. The expectation, and the surprise
# --------------------------------------------------------------------------

def figure_kepler_vs_galaxy(result, art, planets):
    """The Solar System falls off as 1/sqrt(r); a spiral galaxy does not.

    Same encoding in both panels -- ink for what is observed, blue for what the
    visible mass predicts -- so the only thing that changes between them is the
    answer.
    """
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.2, 4.8), gridspec_kw={"wspace": 0.24})

    # --- the Solar System ---------------------------------------------------
    a = np.array([p.a_km for p in planets]) / AU_KM
    v = np.array([p.v_mean_kms for p in planets])
    rr = np.linspace(0.22, 42.0, 500)
    a1.plot(rr, np.sqrt(rp.GM_SUN_KM3_S2 / (rr * AU_KM)), color=ps.C_VISIBLE, lw=LW_MAIN,
            zorder=3)
    for p, x, y in zip(planets, a, v):
        face = ps.SURFACE if p.dwarf else ps.C_OBSERVED
        a1.plot(x, y, "o", color=face, markeredgecolor=ps.C_OBSERVED if p.dwarf else ps.SURFACE,
                markeredgewidth=1.0 if p.dwarf else 0.9, markersize=6.2, zorder=5)
        if p.dwarf:
            offset, ha, va, name = (9, 7), "right", "bottom", f"{p.name}\n(dwarf planet)"
        elif x < 2:
            offset, ha, va, name = (8, 0), "left", "center", p.name
        elif p.name == "Neptune":
            # Below the curve, leaving the space above it for Pluto's label.
            offset, ha, va, name = (5, -7), "left", "top", p.name
        else:
            offset, ha, va, name = (5, 6), "left", "bottom", p.name
        a1.annotate(name, (x, y), xytext=offset, textcoords="offset points", ha=ha, va=va,
                    fontsize=8, color=ps.INK_SECONDARY)
    a1.set_xlim(0, 42)
    a1.set_ylim(0, 54)
    a1.set_xlabel("Distance from the Sun (AU)")
    a1.set_ylabel("Mean orbital speed (km/s)")
    _title(a1, "The Solar System: speeds fall as 1/√r",
           "99.9% of its mass is in the Sun, so the pull fades with distance")
    a1.legend(handles=[
        Line2D([], [], color=ps.C_OBSERVED, marker="o", ls="", markersize=6,
               markeredgecolor=ps.SURFACE, label="planets (NASA fact sheet)"),
        Line2D([], [], color=ps.C_VISIBLE, lw=LW_MAIN, label="Kepler's law, v = √(GM☉ / r)"),
    ], loc="upper right")

    # --- a spiral galaxy ------------------------------------------------------
    c = art["galaxy"].curve
    comp = art["components"]
    _gap(a2, c.r, comp["baryons"], c.v_obs, label="the gap")
    a2.plot(c.r, comp["baryons"], color=ps.C_VISIBLE, lw=LW_MAIN, zorder=4,
            label="predicted from the visible matter")
    # Beyond the last point nearly all the visible mass is enclosed, so its
    # prediction keeps falling the way the planets do.
    r_tail = np.linspace(c.r[-1], c.r[-1] * 1.3, 60)
    a2.plot(r_tail, result.v_bar_last * np.sqrt(c.r[-1] / r_tail), color=ps.C_VISIBLE,
            lw=1.4, ls=(0, (3, 2.5)), zorder=4, label="…and beyond, falling as 1/√r")
    _observed(a2, c.r, c.v_obs, c.e_v_obs, label="observed (Doppler shifts of the gas)")
    a2.set_xlim(0, c.r[-1] * 1.3)
    a2.set_ylim(0, 185)
    a2.set_xlabel("Distance from the galaxy's centre (kpc)")
    a2.set_ylabel("Orbital speed (km/s)")
    _title(a2, f"A spiral galaxy, {result.label}: they don't",
           f"at {result.r_last_kpc:.0f} kpc it orbits at {result.v_obs_last:.0f} km/s, "
           f"{result.v_obs_last / result.v_bar_last:.1f}× what its stars and gas can hold")
    handles, labels = a2.get_legend_handles_labels()
    order = [labels.index(k) for k in ("observed (Doppler shifts of the gas)",
                                       "predicted from the visible matter",
                                       "…and beyond, falling as 1/√r", "the gap")]
    a2.legend([handles[i] for i in order], [labels[i] for i in order], loc="lower right")

    ps.caption(fig, "Left: planet distances and mean orbital speeds from the NASA Planetary Fact "
                    "Sheet; the curve uses only the Sun's mass. Right: "
                    f"{PROVENANCE} The blue curve is computed from where the stars and gas are, "
                    f"with a stellar mass-to-light ratio of {UPSILON_DISK:g} taken from "
                    "stellar-population models; nothing in it is fitted to the rotation.")
    return ps.save(fig, FIGURE_DIR / "01_kepler_vs_galaxy")


# --------------------------------------------------------------------------
# 2. Anatomy of one galaxy
# --------------------------------------------------------------------------

def figure_anatomy(result, art):
    """The flagship galaxy taken apart: components, halo, and the mass budget."""
    c = art["galaxy"].curve
    comp = art["components"]
    model = result.better_halo
    halo = art["halo_curves"][model]
    fit = result.halo_fits[model]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.6, 5.0),
                                 gridspec_kw={"wspace": 0.30, "width_ratios": [1.25, 1.0]})

    # --- rotation curve, taken apart -------------------------------------------
    a1.plot(c.r, comp["gas"], color=ps.C_GAS, lw=LW_THIN, zorder=3)
    a1.plot(c.r, comp["disk"], color=ps.C_DISK, lw=LW_THIN, zorder=3)
    a1.plot(c.r, comp["baryons"], color=ps.C_VISIBLE, lw=LW_MAIN, zorder=4)
    a1.plot(art["r_fine"], halo["halo"], color=ps.C_DARK, lw=LW_MAIN, zorder=4)
    a1.plot(art["r_fine"], halo["total"], color=ps.INK_SECONDARY, lw=1.0, zorder=5)
    _observed(a1, c.r, c.v_obs, c.e_v_obs)
    x_end = c.r[-1]
    _end_labels(a1, [
        (x_end, float(comp["gas"][-1]), "gas"),
        (x_end, float(comp["disk"][-1]), f"stars (M/L = {UPSILON_DISK:g})"),
        (x_end, float(comp["baryons"][-1]), "all visible matter"),
        (x_end, float(halo["halo"][-1]), "dark halo (fitted)"),
        (x_end, float(halo["total"][-1]) + 9, "visible + dark"),
    ], x_text=x_end * 1.05, min_gap=11)
    a1.set_xlim(0, x_end * 1.36)
    a1.set_ylim(0, 178)
    a1.set_xlabel("Distance from the galaxy's centre (kpc)")
    a1.set_ylabel("Orbital speed (km/s)")
    shape = "cored (pseudo-isothermal)" if model == "ISO" else "NFW"
    _title(a1, f"{result.label}: what holds the orbits up",
           f"a {shape} dark halo added to the visible matter matches the observed\nspeeds to "
           f"{fit['rms_residual_kms']:.1f} km/s rms; components add in quadrature")
    a1.legend(handles=[
        Line2D([], [], color=ps.C_OBSERVED, marker="o", ls="", markersize=5,
               markeredgecolor=ps.SURFACE, label="observed"),
        Line2D([], [], color=ps.C_VISIBLE, lw=LW_MAIN, label="visible matter"),
        Line2D([], [], color=ps.C_DARK, lw=LW_MAIN, label="dark halo"),
        Line2D([], [], color=ps.C_GAS, lw=LW_THIN, label="gas"),
        Line2D([], [], color=ps.C_DISK, lw=LW_THIN, label="stars"),
    ], loc="lower right", ncol=2)

    # --- the mass budget -----------------------------------------------------
    m_dyn = art["mass_dynamical"] / 1e10
    m_all = result.baryonic_mass / 1e10
    a2.axhline(m_all, color=ps.C_VISIBLE, lw=LW_MAIN, zorder=3)
    a2.plot(c.r, m_dyn, "o", color=ps.C_OBSERVED, markersize=4.6, markeredgecolor=ps.SURFACE,
            markeredgewidth=0.9, zorder=6)
    a2.annotate(f"every star and gas cloud in the galaxy:\n{_sci(result.baryonic_mass)} M☉",
                (x_end * 0.98, m_all), xytext=(0, 5), textcoords="offset points", ha="right",
                va="bottom", fontsize=8, color=ps.INK_SECONDARY)
    k = int(np.searchsorted(c.r, 0.45 * x_end))
    _label(a2, "mass needed inside each radius\n(from the observed speeds)", (c.r[k], m_dyn[k]),
           (c.r[k] - 1.5, m_dyn[k] + 2.5), ha="right", va="bottom")
    r_x = result.needs_more_than_all_visible_from_kpc
    if r_x:
        a2.axvline(r_x, color=ps.MUTED, lw=0.8, zorder=0)
        a2.annotate(f"beyond {r_x:.0f} kpc the orbits need\nmore than all of it", (r_x, m_dyn.max()),
                    xytext=(5, 0), textcoords="offset points", va="top", fontsize=8,
                    color=ps.INK_SECONDARY)
    ratio = result.dynamical_over_all_visible_last
    x_arrow = x_end + 1.6
    a2.annotate("", xy=(x_arrow, m_dyn[-1]), xytext=(x_arrow, m_all),
                arrowprops={"arrowstyle": "<->", "color": ps.INK_SECONDARY, "lw": 0.9,
                            "shrinkA": 0, "shrinkB": 0})
    a2.text(x_arrow + 1.0, 0.5 * (m_dyn[-1] + m_all), f"{ratio:.1f}×", fontsize=10.5,
            fontweight="bold", color=ps.INK, va="center")
    a2.set_xlim(0, x_end * 1.18)
    a2.set_ylim(0, m_dyn.max() * 1.1)
    a2.set_xlabel("Distance from the galaxy's centre (kpc)")
    a2.set_ylabel("Mass (10¹⁰ M☉)")
    _title(a2, "The mass budget",
           f"inside {result.r_last_kpc:.0f} kpc, gravity needs {_sci(result.dynamical_mass_last)} "
           f"M☉:\n{ratio:.1f}× all the visible mass in the galaxy")

    ps.caption(fig, f"{PROVENANCE} Right: mass needed = V²r/G, exact for a round mass "
                    "distribution and a good approximation where the (roughly round) dark halo "
                    f"dominates. Visible mass = {UPSILON_DISK:g} × stellar light + 1.33 × atomic "
                    "hydrogen, summed over the whole galaxy. The halo is fitted with the stars "
                    f"held at M/L = {UPSILON_DISK:g}.")
    return ps.save(fig, FIGURE_DIR / "02_ngc3198_anatomy")


# --------------------------------------------------------------------------
# 3. Six galaxies
# --------------------------------------------------------------------------

def figure_gallery(results, arts):
    """Six galaxies spanning giant to dwarf: the gap opens in every one."""
    fig, axes = plt.subplots(2, 3, figsize=(12.6, 8.6),
                             gridspec_kw={"hspace": 0.66, "wspace": 0.26})
    fig.subplots_adjust(top=0.83)
    handles = None
    for ax, r in zip(axes.flat, results):
        art = arts[r.name]
        h = _rotation_panel(ax, art)
        if handles is None or len(h) > len(handles):
            handles = h
        where = ("dark matter dominates everywhere" if r.dark_dominates_everywhere
                 else f"dark matter dominates beyond {r.dark_dominates_from_kpc:.0f} kpc")
        _title(ax, f"{r.label} · {r.kind}",
               f"{where}\nat the edge: {r.mass_ratio_last:.1f}× the visible matter's pull")
        ax.set_xlabel("Radius (kpc)")
        ax.set_ylabel("Speed (km/s)")
    # One shared legend names every component (figure 2 labels them directly),
    # led by the data and the visible-matter prediction it is compared with.
    by_label = {h.get_label(): h for h in handles}
    order = ["observed rotation", "all visible matter", "gas",
             f"stellar disc (M/L = {UPSILON_DISK:g})", f"bulge (M/L = {UPSILON_BULGE:g})",
             "the gap: unseen mass"]
    handles = [by_label[k] for k in order if k in by_label]
    fig.legend(handles=handles, loc="upper center", ncol=len(handles), bbox_to_anchor=(0.5, 0.95),
               frameon=False, fontsize=8.5)
    fig.suptitle("Six galaxies, from a giant spiral to a gas-rich dwarf", y=0.99, x=0.075,
                 ha="left")
    ps.caption(fig, f"{PROVENANCE} Stars at M/L = {UPSILON_DISK:g} (disc) and {UPSILON_BULGE:g} "
                    "(bulge) M☉/L☉ at 3.6 µm. Components add in quadrature, not linearly, so "
                    "the blue curve is not the sum of the thin ones. Each panel has its own "
                    "scale: the galaxies span a factor of 6 in speed and 700 in visible mass.")
    return ps.save(fig, FIGURE_DIR / "03_gallery")


# --------------------------------------------------------------------------
# 4. Could the stars just be heavier?
# --------------------------------------------------------------------------

def figure_heavier_stars(results, arts, population):
    """The maximum-disc test, for three galaxies and then for all of them."""
    picks = [r for r in results if r.name in ("NGC3198", "UGC00128", "DDO154")]
    fig = plt.figure(figsize=(12.6, 11.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.55], hspace=0.42, wspace=0.28,
                          top=0.87)

    for i, r in enumerate(picks):
        ax = fig.add_subplot(gs[0, i])
        art = arts[r.name]
        c = art["galaxy"].curve
        md = r.max_disk
        curves = art["max_disk_curves"]
        v_fid = art["components"]["baryons"]
        ax.plot(c.r, v_fid, color=ps.C_VISIBLE, lw=LW_MAIN, zorder=4)
        ax.plot(c.r, curves["max"], color=ps.C_VISIBLE, lw=1.5, ls=DASHED, zorder=4)
        ax.plot(c.r, curves["needed"], color=ps.C_VISIBLE, lw=1.5, ls=DOTTED, zorder=4)
        _observed(ax, c.r, c.v_obs, c.e_v_obs, ms=4.2)
        top = max(float(np.max(curves["needed"])), float(np.max(c.v_obs + c.e_v_obs)))
        ax.set_ylim(0, top * 1.14)
        end = c.r[-1]
        ax.set_xlim(0, end * 1.42)
        # The margin holds labels, not data: no ticks out there.
        ax.set_xticks([t for t in ax.get_xticks() if 0 <= t <= end * 1.02])
        _title(ax, r.label,
               f"the edge needs M/L = {md['upsilon_disk_needed']:.1f}; "
               f"the centre allows {md['upsilon_disk_max']:.2f}")
        ax.set_xlabel("Radius (kpc)")
        if i == 0:
            ax.set_ylabel("Speed (km/s)")
        kn = int(np.argmax(curves["needed"]))
        ax.annotate(f"M/L = {md['upsilon_disk_needed']:.1f}: spins the centre\n"
                    f"{md['inner_excess']:.1f}× faster than observed",
                    (c.r[kn], curves["needed"][kn]), xytext=(8, 0), textcoords="offset points",
                    va="center", fontsize=7.8, color=ps.INK_SECONDARY)
        # Name the three curves where they end, clear of the data.
        _end_labels(ax, [
            (end, float(curves["needed"][-1]), f"M/L {md['upsilon_disk_needed']:.1f}"),
            (end, float(curves["max"][-1]), f"M/L {md['upsilon_disk_max']:.2f}"),
            (end, float(v_fid[-1]), f"M/L {UPSILON_DISK:g}"),
        ], x_text=end * 1.07, min_gap=0.075 * top * 1.14)
    fig.legend(handles=[
        Line2D([], [], color=ps.C_OBSERVED, marker="o", ls="", markersize=5,
               markeredgecolor=ps.SURFACE, label="observed"),
        Line2D([], [], color=ps.C_VISIBLE, lw=LW_MAIN,
               label=f"visible matter, stars at M/L = {UPSILON_DISK:g} (stellar populations)"),
        Line2D([], [], color=ps.C_VISIBLE, lw=1.5, ls=DASHED,
               label="stars as heavy as the inner curve allows (maximum disc)"),
        Line2D([], [], color=ps.C_VISIBLE, lw=1.5, ls=DOTTED,
               label="stars heavy enough to hold the outermost point"),
    ], loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.955), frameon=False, fontsize=8.5)
    fig.suptitle("Could the stars just be heavier than we think?", y=0.99, x=0.075, ha="left")

    # --- every galaxy -----------------------------------------------------------
    ax = fig.add_subplot(gs[1, :2])
    md = population["max_disk"]
    rows = md["rows"]
    x = np.array([row["upsilon_disk_max"] for row in rows])
    y = np.array([row["upsilon_disk_needed"] for row in rows])
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    names = [row["name"] for row, k in zip(rows, ok) if k]
    x, y = x[ok], y[ok]
    lo_ax, hi_ax = 0.05, 300
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.fill_between([lo_ax, hi_ax], [lo_ax, hi_ax], [hi_ax, hi_ax], color=ps.C_DARK, alpha=0.07,
                    linewidth=0, zorder=0)
    ax.plot([lo_ax, hi_ax], [lo_ax, hi_ax], color=ps.INK_SECONDARY, lw=0.9, zorder=2)
    p_lo, p_hi = UPSILON_PLAUSIBLE
    ax.add_patch(plt.Rectangle((p_lo, p_lo), p_hi - p_lo, p_hi - p_lo, facecolor=ps.MUTED,
                               alpha=0.20, edgecolor="none", zorder=1))
    ax.plot(x, y, "o", color=ps.C_OBSERVED, markersize=5.0, markeredgecolor=ps.SURFACE,
            markeredgewidth=0.8, zorder=5)
    label_of = {r.name: r.label for r in picks}
    for n, xi, yi in zip(names, x, y):
        if n in label_of:
            _label(ax, label_of[n], (xi, yi), (xi * 0.28, yi * 2.2), ha="right")
    ax.text(30, 30 * 1.35, "needed = allowed", rotation=45, rotation_mode="anchor",
            fontsize=7.8, color=ps.INK_SECONDARY, transform_rotates_text=True)
    ax.set_xlim(lo_ax, hi_ax)
    ax.set_ylim(lo_ax, hi_ax)
    fmt = FuncFormatter(lambda val, _: f"{val:g}")
    ax.xaxis.set_major_formatter(fmt)
    ax.yaxis.set_major_formatter(fmt)
    ax.set_xlabel("Heaviest stellar M/L the inner rotation curve allows (M☉/L☉)")
    ax.set_ylabel("Stellar M/L needed to hold the outermost point")
    ax.set_aspect("equal", adjustable="box")
    ax.set_anchor("W")
    _title(ax, "The same test for every galaxy in the sample",
           f"{md['n_galaxies']} galaxies; the visible matter gets {md['tolerance_sigma']:.0f}σ of "
           "slack at both ends")

    n_exc = len(md["exceptions"])
    plaus = md["n_explained_with_plausible_upsilon"]
    short = [e for e in md["exceptions"] if UPSILON_PLAUSIBLE[0] <= e["upsilon_disk_needed"]
             <= UPSILON_PLAUSIBLE[1]]
    notes = fig.add_subplot(gs[1, 2])
    notes.axis("off")
    blocks = [
        ("Above the diagonal (shaded)",
         f"{md['n_no_single_upsilon_works']} of {md['n_galaxies']} galaxies. No mass-to-light "
         "ratio, however extreme, lets stars and gas explain both the centre and the edge: "
         "stars heavy enough for the edge overshoot the centre."),
        ("Below the diagonal",
         f"{n_exc} galaxies where a single ratio would do. {n_exc - plaus} of them need stars "
         f"heavier than any known stellar population (M/L above {UPSILON_PLAUSIBLE[1]:g})."),
        ("Inside the grey square",
         f"{plaus} galaxies where a plausible M/L works. All have short rotation curves "
         f"(at most {max(e['n_points'] for e in short) if short else 0} points, ending within "
         f"{max(e['r_last_over_r_disk'] for e in short) if short else 0:.0f} disc scale "
         "lengths): measured too close in to test."),
        ("Why this is robust",
         "An error in a galaxy's distance or tilt rescales both axes by the same factor, so "
         "it slides a point along the diagonal and cannot carry it across."),
    ]
    y0 = 0.98
    for head, body in blocks:
        notes.text(0.0, y0, head, fontsize=9, fontweight="bold", color=ps.INK, va="top",
                   transform=notes.transAxes)
        wrapped = _wrap(body, 44)
        notes.text(0.0, y0 - 0.055, wrapped, fontsize=8.5, color=ps.INK_SECONDARY, va="top",
                   linespacing=1.35, transform=notes.transAxes)
        y0 -= 0.08 + 0.047 * (wrapped.count("\n") + 1)

    ps.caption(fig, f"{PROVENANCE} M/L is the stellar mass-to-light ratio at 3.6 µm, scaled "
                    "together for disc and bulge; the gas mass is measured directly and is not "
                    "scaled. Maximum disc: stars made as heavy as possible without the visible "
                    "matter alone exceeding the observed speed (by more than 2σ) anywhere.")
    return ps.save(fig, FIGURE_DIR / "04_heavier_stars")


def _wrap(text: str, width: int) -> str:
    import textwrap

    return "\n".join(textwrap.wrap(text, width))


# --------------------------------------------------------------------------
# 5. The dark haloes
# --------------------------------------------------------------------------

def figure_halos(results, arts):
    """What the unseen mass must supply, and two halo shapes fitted to it."""
    fig, axes = plt.subplots(2, 3, figsize=(12.6, 8.6),
                             gridspec_kw={"hspace": 0.66, "wspace": 0.26})
    fig.subplots_adjust(top=0.83)
    for ax, r in zip(axes.flat, results):
        art = arts[r.name]
        c = art["galaxy"].curve
        dm = art["dark_required"]
        ax.errorbar(c.r, dm["v"], yerr=[dm["v"] - dm["v_lo"], dm["v_hi"] - dm["v"]], fmt="o",
                    color=ps.C_OBSERVED, ecolor=ps.C_OBSERVED_ERR, elinewidth=0.9, capsize=0,
                    markersize=4.2, markeredgecolor=ps.SURFACE, markeredgewidth=0.8, zorder=6)
        for model, style, lw in (("ISO", "-", LW_MAIN), ("NFW", DASHED, 1.6)):
            ax.plot(art["r_fine"], art["halo_curves"][model]["halo"], color=ps.C_DARK, lw=lw,
                    ls=style, zorder=4)
        iso, nfw = r.halo_fits["ISO"], r.halo_fits["NFW"]
        d = nfw["derived"]
        _title(ax, r.label,
               f"misfit (rms): cored {iso['rms_residual_kms']:.1f}, NFW "
               f"{nfw['rms_residual_kms']:.1f} km/s\nNFW c = {nfw['params']['concentration']:.1f}"
               f"; ΛCDM expects {d['concentration_lcdm']:.1f} "
               f"({d['concentration_offset_sigma']:+.1f}σ)")
        ax.set_ylim(min(0.0, float(np.nanmin(dm["v_lo"])) * 1.1),
                    float(np.nanmax(dm["v_hi"])) * 1.15)
        ax.set_xlim(0, c.r[-1] * 1.04)
        ax.axhline(0, color=ps.AXIS, lw=0.8, zorder=1)
        ax.set_xlabel("Radius (kpc)")
        ax.set_ylabel("Speed (km/s)")
    fig.legend(handles=[
        Line2D([], [], color=ps.C_OBSERVED, marker="o", ls="", markersize=5,
               markeredgecolor=ps.SURFACE,
               label="what the unseen mass must supply: √(V_observed² − V_visible²)"),
        Line2D([], [], color=ps.C_DARK, lw=LW_MAIN, label="cored halo (pseudo-isothermal)"),
        Line2D([], [], color=ps.C_DARK, lw=1.6, ls=DASHED,
               label="cuspy halo (NFW, the shape cold-dark-matter simulations produce)"),
    ], loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.955), frameon=False, fontsize=8.5)
    fig.suptitle("The dark haloes: how much, and what shape", y=0.99, x=0.075, ha="left")
    ps.caption(fig, f"{PROVENANCE} Stars held at M/L = {UPSILON_DISK:g}; each halo shape has two "
                    "free parameters. Points below zero mark radii where the visible matter alone "
                    "already slightly over-predicts. ΛCDM concentrations: Dutton & Macciò (2014), "
                    "median for the fitted halo mass, scatter 0.11 dex. Error bars are random "
                    "errors only.")
    return ps.save(fig, FIGURE_DIR / "05_dark_halos")


# --------------------------------------------------------------------------
# 6. The radial acceleration relation
# --------------------------------------------------------------------------

def figure_rar(population, pop_art, comparison, results):
    """Every point of every galaxy: observed acceleration against the visible matter's."""
    from sparc import load_published_rar

    data = pop_art["rar"]
    rar = population["rar"]
    bins = population["rar_bins"]
    _, pub_bins = load_published_rar()
    gd = rar["g_dagger"]

    fig, (ax, axr) = plt.subplots(2, 1, figsize=(8.4, 9.8), sharex=True,
                                  gridspec_kw={"height_ratios": [3.2, 1.0], "hspace": 0.07})
    fig.subplots_adjust(bottom=0.07)
    xg = 10.0 ** data.log_gbar
    yg = 10.0 ** data.log_gobs
    ax.plot(xg, yg, "o", color=ps.C_CONTEXT, markersize=2.6, markeredgewidth=0, rasterized=True,
            zorder=1, label=f"{len(data):,} points from {rar['n_galaxies']} galaxies")
    xx = np.logspace(-12.6, -8.3, 300)
    ax.plot(xx, xx, color=ps.C_VISIBLE, lw=1.8, zorder=3,
            label="if the visible matter were everything")
    ax.plot(xx, rp.rar_mcgaugh(xx, gd), color=ps.C_DARK, lw=LW_MAIN, zorder=4,
            label=f"fitted relation, g† = {gd / 1e-10:.2f} × 10⁻¹⁰ m/s²")

    # Two galaxies with the same rotation speed but very different densities:
    # the dense one lives at high accelerations, the diffuse one at low.
    labels = {r.name: r.label for r in results}
    for name, marker, face, text in (
        ("NGC2403", "o", ps.C_OBSERVED, "bright, compact spiral"),
        ("UGC00128", "s", ps.SURFACE, "diffuse spiral, same speed"),
    ):
        m = data.galaxy == name
        ax.plot(xg[m], yg[m], marker, color=face, markeredgecolor=ps.C_OBSERVED,
                markeredgewidth=0.9, markersize=4.4, zorder=6,
                label=f"{labels.get(name, name)}: {text}")
    ax.plot(10.0**pub_bins.log_gbar, 10.0**pub_bins.log_gobs, "D", color=ps.C_PUBLISHED,
            markersize=7.0, markerfacecolor="none", markeredgewidth=1.3, zorder=7,
            label="published bin means (Lelli et al. 2017)")
    bx = np.array([10.0**v if v is not None else np.nan for v in bins["log_gbar"]])
    by = np.array([10.0**v if v is not None else np.nan for v in bins["log_gobs"]])
    ax.plot(bx, by, "D", color=ps.C_OBSERVED, markersize=2.8, zorder=8, label="our bin means")

    # The gap between the two lines is the mass discrepancy.
    x0 = 1.6e-12
    ax.annotate("", xy=(x0, rp.rar_mcgaugh(x0, gd)), xytext=(x0, x0),
                arrowprops={"arrowstyle": "<->", "color": ps.INK_SECONDARY, "lw": 0.9,
                            "shrinkA": 0, "shrinkB": 0})
    ax.text(x0 * 1.18, np.sqrt(x0 * rp.rar_mcgaugh(x0, gd)) * 0.72,
            "the missing pull:\nunseen mass", fontsize=8, color=ps.INK_SECONDARY, va="center")
    ax.axvline(gd, color=ps.MUTED, lw=0.8, zorder=0)
    ax.text(gd * 1.12, 10**-11.8, "g†", fontsize=9, color=ps.INK_SECONDARY, va="bottom")
    ax.set_xscale("log")
    ax.set_yscale("log")
    _log_ticks(ax.xaxis)
    _log_ticks(ax.yaxis)
    ax.set_xlim(10**-12.4, 10**-8.4)
    ax.set_ylim(10**-11.9, 10**-8.4)
    ax.set_ylabel("Observed acceleration, V²/r (m/s²)")
    _title(ax, "One relation for every galaxy",
           "at every radius of every galaxy: the pull the orbit needs, against the pull "
           "the visible matter supplies")
    ax.legend(loc="upper left", fontsize=7.8)

    resid = data.log_gobs - np.log10(rp.rar_mcgaugh(xg, gd))
    s = rar["residual_rms_dex"]
    axr.axhspan(-s, s, color=ps.MUTED, alpha=0.14, lw=0, zorder=0)
    axr.plot(xg, resid, "o", color=ps.C_CONTEXT, markersize=2.4, markeredgewidth=0,
             rasterized=True, zorder=1)
    axr.axhline(0, color=ps.C_DARK, lw=1.6, zorder=3)
    axr.text(10**-12.3, -0.72, f"band: ±{s:.2f} dex rms, a factor of {10**s:.2f}", fontsize=8,
             color=ps.INK_SECONDARY, va="bottom")
    axr.set_xscale("log")
    _log_ticks(axr.xaxis)
    axr.set_ylim(-0.8, 0.8)
    axr.yaxis.set_major_locator(MultipleLocator(0.4))
    axr.set_ylabel("Residual (dex)")
    axr.set_xlabel("Acceleration the visible matter alone supplies (m/s²)")
    c = comparison["rar"]
    ps.caption(fig, f"{PROVENANCE} Quality cuts of McGaugh, Lelli & Schombert (2016): quality flag "
                    "1–2, inclination ≥ 30°, velocity known to 10%. Our points reproduce "
                    f"{c['points_matched_within_0p015']:,} of the {c['n_points_published']:,} "
                    "published ones to within 0.015 dex. For scale: Earth's surface gravity is "
                    "9.8 m/s², about a hundred billion times g†.")
    return ps.save(fig, FIGURE_DIR / "06_radial_acceleration")


# --------------------------------------------------------------------------
# 7. The baryonic Tully-Fisher relation
# --------------------------------------------------------------------------

#: Where each showcase galaxy's label goes (speed km/s, mass Msun): in the empty
#: space above or below the relation, tied to its point by a leader line.
BTFR_LABELS = {
    "UGC02885": (175, 9e11),
    "NGC2841": (345, 1.5e10),
    "NGC3198": (330, 3.5e9),
    "NGC2403": (230, 8e8),
    "UGC00128": (62, 8e10),
    "DDO154": (26, 3e9),
}


def figure_btfr(population, pop_art, results):
    """Total visible mass against flat rotation speed -- and why the gas matters."""
    d = pop_art["btfr"]
    fit = population["btfr"]
    st = fit["stellar_only"]
    v = 10.0**d.log_vflat
    m_bar = 10.0**d.log_mbar
    m_star = 10.0**d.log_mstar

    fig, ax = plt.subplots(figsize=(8.6, 6.6))
    for vi, a, b in zip(v, m_star, m_bar):
        ax.plot([vi, vi], [a, b], color=ps.C_CONTEXT, lw=0.8, zorder=1)
    h_star, = ax.plot(v, m_star, "o", color=ps.SURFACE, markeredgecolor=ps.MUTED,
                      markeredgewidth=0.9, markersize=4.6, zorder=2)
    h_bar = ax.errorbar(v, m_bar, xerr=[v - 10 ** (d.log_vflat - d.e_log_vflat),
                                        10 ** (d.log_vflat + d.e_log_vflat) - v],
                        fmt="o", color=ps.C_OBSERVED, ecolor=ps.C_CONTEXT, elinewidth=0.7,
                        capsize=0, markersize=4.8, markeredgecolor=ps.SURFACE,
                        markeredgewidth=0.8, zorder=4)
    vv = np.logspace(np.log10(15), np.log10(400), 100)
    h_fit, = ax.plot(vv, 10 ** (fit["slope"] * (np.log10(vv) - fit["pivot"])
                                + fit["zero_point_at_pivot"]),
                     color=ps.C_DARK, lw=LW_MAIN, zorder=3)
    lbl = {r.name: r.label for r in results}
    for n, vi, mb in zip(d.name, v, m_bar):
        if n in BTFR_LABELS:
            tx, ty = BTFR_LABELS[n]
            _label(ax, lbl.get(n, n), (vi, mb), (tx, ty), ha="center",
                   va="bottom" if ty > mb else "top")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(18, 420)
    ax.set_ylim(1e6, 2e12)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.xaxis.set_minor_formatter(FuncFormatter(
        lambda x, _: f"{x:g}" if x in (20, 30, 50, 200, 300) else ""))
    ax.set_xlabel("Speed along the flat outer part of the rotation curve (km/s)")
    ax.set_ylabel("Mass (M☉)")
    _title(ax, "The baryonic Tully–Fisher relation",
           f"{fit['n_galaxies']} galaxies: the flat speed predicts the visible mass, once the "
           "gas is counted\n(grey lines add each galaxy's gas to its stars)")
    ax.legend([h_bar, h_star, h_fit],
              ["stars + gas (all visible mass)", "stars only",
               f"fit to stars + gas: $M \\propto V^{{{fit['slope']:.2f}}}$"], loc="upper left")
    ax.text(0.98, 0.04, f"scatter in mass about the line: {fit['rms_vertical_dex']:.2f} dex with "
                        f"gas, {st['rms_vertical_dex']:.2f} dex stars only\nintrinsic scatter "
                        f"perpendicular to the line: {fit['sigma_int_dex']:.3f} dex",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color=ps.INK_SECONDARY)
    ps.caption(fig, f"{PROVENANCE} Visible mass = {UPSILON_DISK:g} × disc light + "
                    f"{UPSILON_BULGE:g} × bulge light + 1.33 × atomic hydrogen (the 1.33 adds "
                    "helium). The fit allows intrinsic scatter perpendicular to the line, as "
                    "Lelli et al. (2019) did. Mass error bars (mostly from distances) are omitted "
                    "for clarity.")
    return ps.save(fig, FIGURE_DIR / "07_tully_fisher")


# --------------------------------------------------------------------------
# 0. Grading ourselves
# --------------------------------------------------------------------------

def figure_comparison(comparison, population):
    """Our numbers against the published ones, as percent differences."""
    c_rar, c_btfr, c_md = comparison["rar"], comparison["btfr"], comparison["max_disk"]
    rows = [
        # label, ours, published, published fractional uncertainty (or None)
        ("RAR acceleration scale g†", c_rar["g_dagger_ours"], c_rar["g_dagger_published"],
         0.24 / 1.20),
        ("↳ published points, our fitting code", c_rar["g_dagger_published_data_refit"],
         c_rar["g_dagger_published"], 0.24 / 1.20),
        ("RAR scatter (rms)", c_rar["rms_ours"], c_rar["rms_published"], None),
        ("RAR scatter (core width)", c_rar["gaussian_width_ours"],
         c_rar["gaussian_width_published"], None),
        ("RAR points passing the cuts", c_rar["n_points_ours"], c_rar["n_points_published"], None),
        ("Tully–Fisher slope", c_btfr["slope_ours"], c_btfr["slope_published"],
         c_btfr["slope_err_published"] / c_btfr["slope_published"]),
        ("↳ published table, our fitting code", c_btfr["slope_published_data_refit"],
         c_btfr["slope_published"], c_btfr["slope_err_published"] / c_btfr["slope_published"]),
        ("Galaxy visible masses (median)", 10 ** c_btfr["log_mbar_median_diff_dex"], 1.0, None),
        ("Maximum-disc M/L (median)", c_md["median_ratio_ours_over_published"], 1.0, None),
    ]
    labels = [r[0] for r in rows]
    diffs = np.array([100.0 * (r[1] / r[2] - 1.0) for r in rows])
    unc = [100.0 * r[3] if r[3] else None for r in rows]

    fig, ax = plt.subplots(figsize=(9.4, 0.46 * len(rows) + 2.2))
    ypos = np.arange(len(rows))
    for y, u in zip(ypos, unc):
        if u:
            ax.add_patch(plt.Rectangle((-u, y - 0.32), 2 * u, 0.64, facecolor=ps.C_PUBLISHED,
                                       alpha=0.12, edgecolor="none", zorder=1))
    ax.barh(ypos, diffs, height=0.36, color=ps.C_DARK, zorder=3)
    ax.axvline(0.0, color=ps.INK_SECONDARY, linewidth=1.0, zorder=4)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=8.5, color=ps.INK_SECONDARY)
    ax.tick_params(axis="y", length=0)
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.set_xlim(-26, 26)
    for y, v in zip(ypos, diffs):
        off = 5 if v >= 0 else -5
        ax.annotate(f"{v:+.1f}%".replace("-", "−"), xy=(v, y), xytext=(off, 0),
                    textcoords="offset points",
                    va="center", ha="left" if v >= 0 else "right", fontsize=8,
                    color=ps.INK_SECONDARY)
    ax.set_xlabel("Difference from the published value (%)")
    _title(ax, "Our numbers against the literature",
           "violet band: the published uncertainty, where one is quoted "
           "(for g†, its ±20% systematic)")
    ps.caption(fig, "Published: McGaugh, Lelli & Schombert 2016 (RAR fit); Lelli et al. 2017 (RAR "
                    "points); Lelli et al. 2019 (Tully–Fisher); Starkman et al. 2018 (maximum "
                    "disc). The rows marked ↳ run the published data through our own fitting "
                    "code. For the Tully–Fisher slope it returns the published value exactly; "
                    "for g† it lands where our own data do, so that difference lies in fitting "
                    "details we could not reproduce, not in the data.")
    return ps.save(fig, FIGURE_DIR / "00_comparison")


# --------------------------------------------------------------------------
# Any other galaxy, on request
# --------------------------------------------------------------------------

def figure_single_galaxy(result, art):
    """Rotation curve and required halo for one galaxy chosen with --galaxy."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.0, 4.6), gridspec_kw={"wspace": 0.26})
    handles = _rotation_panel(a1, art)
    a1.legend(handles=handles, loc="lower right", fontsize=7.6)
    a1.set_xlabel("Distance from the galaxy's centre (kpc)")
    a1.set_ylabel("Orbital speed (km/s)")
    _title(a1, f"{result.label} ({result.hubble_type})",
           f"at {result.r_last_kpc:.0f} kpc: {result.mass_ratio_last:.1f}× the visible matter's "
           "pull")
    c = art["galaxy"].curve
    dm = art["dark_required"]
    a2.errorbar(c.r, dm["v"], yerr=[dm["v"] - dm["v_lo"], dm["v_hi"] - dm["v"]], fmt="o",
                color=ps.C_OBSERVED, ecolor=ps.C_OBSERVED_ERR, elinewidth=0.9, capsize=0,
                markersize=4.2, markeredgecolor=ps.SURFACE, markeredgewidth=0.8, zorder=6,
                label="required from the data")
    for model, style, lw in (("ISO", "-", LW_MAIN), ("NFW", DASHED, 1.6)):
        f = result.halo_fits[model]
        a2.plot(art["r_fine"], art["halo_curves"][model]["halo"], color=ps.C_DARK, lw=lw,
                ls=style, zorder=4,
                label=f"{'cored' if model == 'ISO' else 'NFW'} halo, "
                      f"misfit {f['rms_residual_kms']:.1f} km/s rms")
    a2.axhline(0, color=ps.AXIS, lw=0.8)
    a2.set_xlim(0, c.r[-1] * 1.04)
    a2.legend(loc="lower right", fontsize=7.6)
    a2.set_xlabel("Distance from the galaxy's centre (kpc)")
    a2.set_ylabel("Speed (km/s)")
    _title(a2, "What the unseen mass must supply")
    ps.caption(fig, f"{PROVENANCE} Stars at M/L = {UPSILON_DISK:g} (disc), "
                    f"{UPSILON_BULGE:g} (bulge).")
    return ps.save(fig, FIGURE_DIR / f"galaxy_{result.name.lower()}")


def make_all(results, arts, population, pop_art, comparison, flagship):
    """Every standard figure; returns the files written."""
    from sparc import load_planets

    flag = next(r for r in results if r.name == flagship)
    written = []
    written += figure_comparison(comparison, population)
    written += figure_kepler_vs_galaxy(flag, arts[flagship], load_planets())
    written += figure_anatomy(flag, arts[flagship])
    written += figure_gallery(results, arts)
    written += figure_heavier_stars(results, arts, population)
    written += figure_halos(results, arts)
    written += figure_rar(population, pop_art, comparison, results)
    written += figure_btfr(population, pop_art, results)
    return written
