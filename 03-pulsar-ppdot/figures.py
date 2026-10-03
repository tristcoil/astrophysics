"""Every figure, drawn from the classified table and the results dictionary.

Figures are numbered in reading order:

    01_ppdot            the P-Pdot diagram -- the keeper
    02_spin_periods     how fast neutron stars spin, and how fast they could
    03_galaxy           where pulsars are, and how they drift out of the disc
    04_recycling        the evidence that millisecond pulsars were spun up
    05_caveats          where the diagram's numbers stop meaning what they say
    ppdot_interactive   the P-Pdot diagram with every pulsar hoverable (HTML)
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import plotstyle as ps
import pulsar_physics as pp
from analysis import NON_RECYCLED
from config import FAMOUS, FIGURE_DIR, MANIFEST_PATH, MSP_MAX_PERIOD_S

ps.use_style()

#: How each family is drawn. Only three hues are used in any one scatter -- the
#: most a colour-blind-safe palette supports when every pair of classes can end
#: up side by side -- so the remaining distinctions ride on marker shape.
STYLE = {
    "ordinary": dict(color=ps.MUTED, marker="o", s=5, alpha=0.45, z=2,
                     label="ordinary radio pulsars"),
    "msp": dict(color=ps.BLUE, marker="o", s=9, alpha=0.9, z=4,
                label="millisecond pulsars"),
    "radio_quiet": dict(color=ps.AQUA, marker="o", s=15, alpha=0.95, z=5,
                        label="radio-quiet pulsars (γ-ray/X-ray only)"),
    "xins": dict(color=ps.AQUA, marker="s", s=20, alpha=0.95, z=6,
                 label="X-ray isolated neutron stars"),
    "magnetar": dict(color=ps.ORANGE, marker="D", s=20, alpha=0.95, z=6,
                     label="magnetars"),
}
#: Rotating radio transients are radio pulsars found through single bright
#: pulses; physically they belong with the ordinary pulsars and are drawn so.
FAMILY_OF = {"normal": "ordinary", "rrat": "ordinary", "msp": "msp",
             "radio_quiet": "radio_quiet", "xins": "xins", "magnetar": "magnetar"}

C_GUIDE = ps.AXIS  # constant-age / constant-field lines
C_VALLEY = ps.GRID  # the death-valley wash


def _provenance() -> str:
    try:
        m = json.loads(MANIFEST_PATH.read_text())
        return (f"Data: ATNF Pulsar Catalogue v{m['version']} (Manchester et al. 2005, "
                f"AJ 129, 1993), retrieved {m['retrieved_utc'][:10]}.")
    except (OSError, KeyError, ValueError):
        return "Data: ATNF Pulsar Catalogue (Manchester et al. 2005, AJ 129, 1993)."


def _family(t):
    return t["klass"].map(FAMILY_OF)


_SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def _sci(x: float, digits: int = 1) -> str:
    """7.2e13 -> '7.2×10¹³', for captions."""
    exp = int(np.floor(np.log10(abs(x))))
    return f"{x / 10**exp:.{digits}f}×10{str(exp).translate(_SUPERSCRIPT)}"


def _screen_angle(ax, x0, y0, x1, y1) -> float:
    """Angle on screen, in degrees, of the segment between two data points."""
    (sx0, sy0), (sx1, sy1) = ax.transData.transform([(x0, y0), (x1, y1)])
    return float(np.degrees(np.arctan2(sy1 - sy0, sx1 - sx0)))


def _line_label(ax, func, x, text, *, dy_pts=3.0, color=ps.MUTED, fontsize=7, **kw):
    """Label a guide line at ``x``, rotated to lie along it, nudged off the line."""
    y = func(x)
    ang = _screen_angle(ax, x / 1.2, func(x / 1.2), x * 1.2, func(x * 1.2))
    ax.annotate(text, (x, y), xytext=(0, dy_pts), textcoords="offset points",
                rotation=ang, rotation_mode="anchor", ha="center", va="bottom",
                fontsize=fontsize, color=color, **kw)


def _scatter(ax, t, mask, style, **overrides):
    st = {**style, **overrides}
    return ax.scatter(t.loc[mask, "P"], t.loc[mask, "Pdot"], s=st["s"], c=st["color"],
                      marker=st["marker"], alpha=st["alpha"], linewidths=0.3,
                      edgecolors=ps.SURFACE, zorder=st["z"], rasterized=False)


def _ppdot_axes(ax, xlim, ylim):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("spin period P (seconds)")
    ax.set_ylabel("period derivative Ṗ (seconds per second)")


# --- 01 The P-Pdot diagram ----------------------------------------------------------------


def figure_ppdot(t, results) -> list:
    fig, ax = plt.subplots(figsize=(9.6, 7.8))
    xlim, ylim = (7e-4, 3e2), (1e-22, 1e-8)
    _ppdot_axes(ax, xlim, ylim)
    ax.grid(False)  # the guide lines are the grid on this plot
    P = np.logspace(np.log10(xlim[0]), np.log10(xlim[1]), 500)

    # Death valley: between the pure-dipole and curved-field death lines.
    upper = pp.death_line_pdot(P, "dipole")
    lower = pp.death_line_pdot(P, "multipole")
    ax.fill_between(P, lower, upper, color=C_VALLEY, alpha=0.7, lw=0, zorder=1)
    for y in (upper, lower):
        ax.plot(P, y, color=ps.MUTED, lw=0.8, ls=(0, (4, 2.5)), zorder=1)

    # Constant characteristic age (slope +1) and constant surface field (slope -1),
    # each labelled where its line runs through empty space.
    # (The 1 kyr line meets the death line at the right edge, so it is labelled
    # at the left, where it runs past the Crab.)
    for age, text, x in ((1e3, "τ = 1 kyr", 0.012), (1e5, "100 kyr", 175.0),
                         (1e7, "10 Myr", 175.0), (1e9, "1 Gyr", 175.0)):
        ax.plot(P, pp.pdot_for_age(P, age), color=C_GUIDE, lw=0.7, zorder=1)
        _line_label(ax, lambda xx, a=age: pp.pdot_for_age(xx, a), x, text)
    for field, text, x in ((1e8, "B = 10⁸ G", 0.022), (1e10, "10¹⁰ G", 40.0),
                           (1e12, "10¹² G", 150.0), (1e14, "10¹⁴ G", 0.3)):
        ax.plot(P, pp.pdot_for_field(P, field), color=C_GUIDE, lw=0.7, zorder=1)
        _line_label(ax, lambda xx, b=field: pp.pdot_for_field(xx, b), x, text)

    fam = _family(t)
    samp = t["in_ppdot"]
    handles = []
    for key in ("ordinary", "msp", "radio_quiet", "xins", "magnetar"):
        m = samp & (fam == key)
        st = STYLE[key]
        _scatter(ax, t, m, st)
        handles.append(Line2D([], [], ls="", marker=st["marker"], color=st["color"],
                              markersize=np.sqrt(st["s"]) + 1.5,
                              label=f"{st['label']} ({m.sum():,})"))
    binaries = samp & t["is_binary"]
    ax.scatter(t.loc[binaries, "P"], t.loc[binaries, "Pdot"], s=34, facecolors="none",
               edgecolors=ps.INK_SECONDARY, linewidths=0.55, zorder=7)
    handles.append(Line2D([], [], ls="", marker="o", markerfacecolor="none",
                          markeredgecolor=ps.INK_SECONDARY, markersize=6.5,
                          label=f"ringed: has a binary companion ({binaries.sum():,})"))

    # Region labels, in ink -- the legend carries colour identity.
    region = dict(fontsize=10, color=ps.INK, fontweight="bold", zorder=8)
    ax.text(1.1e-3, 1.2e-18, "millisecond\npulsars", **region, ha="left", va="center")
    ax.text(0.03, 3e-15, "ordinary\npulsars", **region, ha="center", va="center")
    ax.text(1.5, 3e-9, "magnetars", **region, ha="center")
    ax.text(40, 3e-19, "pulsar\ngraveyard", **{**region, "color": ps.INK_SECONDARY},
            ha="center", va="center")
    _line_label(ax, lambda x: np.sqrt(pp.death_line_pdot(x, "dipole")
                                      * pp.death_line_pdot(x, "multipole")),
                28.0, "death valley", dy_pts=-4, color=ps.INK_SECONDARY, fontsize=8.5,
                fontstyle="italic")

    # The famous ones: ringed in ink so the eye can find them, then labelled.
    for jname, spec in FAMOUS.items():
        r = t.loc[t["jname"] == jname]
        if r.empty:
            continue
        x, y = float(r["P"].iloc[0]), float(r["Pdot"].iloc[0])
        st = STYLE[FAMILY_OF[r["klass"].iloc[0]]]
        ax.scatter([x], [y], s=30, c=st["color"], marker=st["marker"], edgecolors=ps.INK,
                   linewidths=0.9, zorder=9)
        text = spec["label"] + (f"\n{spec['note']}" if spec["note"] else "")
        ps.annotate(ax, text, (x, y), spec["offset"], textcoords="offset points",
                    color=ps.INK_SECONDARY, ha=spec["ha"], va="center", zorder=9)
    # The Double Pulsar's two members are one system: tie them together.
    a = t.loc[t["jname"] == "J0737-3039A"].iloc[0]
    b = t.loc[t["jname"] == "J0737-3039B"].iloc[0]
    ax.plot([a["P"], b["P"]], [a["Pdot"], b["Pdot"]], color=ps.INK_SECONDARY, lw=0.6,
            ls=(0, (1, 2)), zorder=3)

    leg = ax.legend(handles=handles, loc="upper left", fontsize=7.5, handletextpad=0.3,
                    borderaxespad=0.4, labelspacing=0.45, frameon=True,
                    facecolor=ps.SURFACE, edgecolor="none", framealpha=1.0)
    leg.set_zorder(10)
    n = results["census"]
    gc = results["globular_clusters"]
    ax.set_title(f"The P–Ṗ diagram: {n['n_in_ppdot']:,} neutron stars sorted by spin "
                 f"and spin-down")
    ps.caption(fig, (
        f"{_provenance()} Every object with a measured positive spin-down is shown, except "
        f"{gc['n_positive_pdot']} globular-cluster pulsars whose Ṗ is dominated by the "
        f"cluster's gravity ({gc['negative_pdot']['k']} more have negative Ṗ and could not "
        f"be plotted at all). Rotating radio transients are included with the ordinary "
        f"pulsars. Guide lines assume magnetic-dipole braking: τ = P/2Ṗ, "
        f"B = 3.2×10¹⁹ √(PṖ) G (I = 10⁴⁵ g cm², R = 10 km). Death valley after Chen & "
        f"Ruderman (1993), lines from Zhang, Harding & Muslimov (2000)."))
    return ps.save(fig, FIGURE_DIR / "01_ppdot")


# --- 02 Spin periods ------------------------------------------------------------------------


def _reciprocal(x):
    """1/x for the frequency axis. matplotlib probes the axis ends, including 0."""
    x = np.array(x, dtype=float)
    out = np.full_like(x, np.inf)
    nonzero = x != 0
    out[nonzero] = 1.0 / x[nonzero]
    return out


HALO = [pe.withStroke(linewidth=2.5, foreground=ps.SURFACE)]


def figure_spin_periods(t, results) -> list:
    fig, ax = plt.subplots(figsize=(9.6, 5.0))
    has = t["P"].notna()
    fam = _family(t)
    bins = np.logspace(-3, 2, 76)  # 15 bins per decade
    groups = (("msp", "millisecond pulsars"), ("ordinary", "ordinary radio pulsars"),
              ("radio_quiet", "radio-quiet and X-ray isolated"), ("magnetar", "magnetars"))
    data, colors, labels = [], [], []
    for key, label in groups:
        m = has & ((fam == key) | ((key == "radio_quiet") & (fam == "xins")))
        data.append(t.loc[m, "P"].to_numpy())
        colors.append(STYLE[key]["color"])
        labels.append(f"{label} ({m.sum():,})")
    counts, _, _ = ax.hist(data, bins=bins, stacked=True, color=colors, label=labels,
                           edgecolor=ps.SURFACE, linewidth=0.6, zorder=3)
    ymax = float(np.max(counts[-1])) * 1.25
    ax.set_ylim(0, ymax)
    ax.set_xscale("log")
    ax.set_xlim(7e-4, 1.2e2)
    ax.set_xlabel("spin period P (seconds)")
    ax.set_ylabel("pulsars per bin (15 bins per decade)")
    ax.grid(axis="x", visible=False)
    top = ax.secondary_xaxis("top", functions=(_reciprocal, _reciprocal))
    top.set_xlabel("spins per second (Hz)", color=ps.INK_SECONDARY)
    top.tick_params(colors=ps.MUTED, labelsize=8)

    # How fast could a neutron star spin before flinging off its own equator?
    lo = 1.0 / pp.max_spin_frequency_hz(1.4, 10.0)
    hi = 1.0 / pp.max_spin_frequency_hz(1.4, 14.0)
    ax.axvspan(lo, hi, color=ps.GRID, alpha=0.9, lw=0, zorder=1)
    ax.text(np.sqrt(lo * hi), ymax * 0.97, "break-\nup\nspin", ha="center", va="top",
            fontsize=7.5, color=ps.INK_SECONDARY)

    ex = results["extremes"]
    sd = results["spin_distribution"]
    ps.annotate(ax, f"fastest: {ex['fastest']['jname']}\n"
                    f"{ex['fastest_spin_hz']:.0f} spins per second",
                (ex["fastest"]["P_s"], 2.0), (1.25e-3, 0.45 * ymax), textcoords="data",
                ha="left", va="bottom")
    ps.annotate(ax, f"slowest: {ex['slowest']['jname']}\none turn every "
                    f"{ex['slowest']['P_s']:.0f} s",
                (ex["slowest"]["P_s"], 2.0), (ex["slowest"]["P_s"], 0.3 * ymax),
                textcoords="data", ha="center", va="bottom")
    ax.text(sd["median_msp_ms"] / 1e3, 0.66 * ymax,
            f"millisecond pulsars\nmedian {sd['median_msp_ms']:.1f} ms",
            ha="center", fontsize=8.5, color=ps.INK)
    ax.text(sd["median_other_s"], 0.97 * ymax,
            f"ordinary pulsars\nmedian {sd['median_other_s']:.2f} s",
            ha="center", va="top", fontsize=8.5, color=ps.INK)
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.86), fontsize=7.5)
    ax.set_title("Two populations: pulsars spin either hundreds of times a second "
                 "or about once", pad=28)
    ps.caption(fig, (
        f"{_provenance()} All {sd['n']:,} catalogued pulsars with a known period, including "
        f"those in globular clusters and those without a measured Ṗ. Break-up band: the "
        f"mass-shedding limit ν_max = 1045 (M/M☉)^½ (10 km/R)^³ᐟ² Hz (Lattimer & Prakash "
        f"2004) for M = 1.4 M☉ and radii 10–14 km. The {ex['fastest_spin_hz']:.0f} Hz pulsar "
        f"requires a mean density above {_sci(ex['fastest_min_density_g_cm3'])} g/cm³ to hold "
        f"together, and a radius below {ex['fastest_max_radius_km_1p4']:.1f} km for "
        f"1.4 M☉ (as found by Hessels et al. 2006)."))
    return ps.save(fig, FIGURE_DIR / "02_spin_periods")


# --- 03 The Galaxy ------------------------------------------------------------------------------


def _sky_xy(gl_deg, gb_deg):
    """Galactic coordinates to Mollweide axes: l = 0 in the centre, increasing to
    the left, as astronomers draw the sky."""
    l_wrapped = (np.asarray(gl_deg, dtype=float) + 180.0) % 360.0 - 180.0
    return -np.radians(l_wrapped), np.radians(np.asarray(gb_deg, dtype=float))


def _sky_box(ax, l0, l1, b0, b1, **kw):
    """Outline a Galactic longitude-latitude box (must not straddle l = 180)."""
    ls = np.linspace(l0, l1, 80)
    bs = np.linspace(b0, b1, 20)
    edge_l = np.concatenate([ls, np.full_like(bs, l1), ls[::-1], np.full_like(bs, l0)])
    edge_b = np.concatenate([np.full_like(ls, b0), bs, np.full_like(ls, b1), bs[::-1]])
    x, y = _sky_xy(edge_l, edge_b)
    ax.plot(x, y, **kw)


def figure_galaxy(t, results) -> list:
    fig = plt.figure(figsize=(10, 10.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.05, 1], hspace=0.34, bottom=0.06,
                          top=0.96)
    ax = fig.add_subplot(gs[0], projection="mollweide")
    ax2 = fig.add_subplot(gs[1])
    h = results["heights"]
    lat = h["latitude"]

    # (a) the sky
    fam = _family(t)
    has = t["gl"].notna()
    young = has & fam.isin(("ordinary", "radio_quiet", "xins")) & ~t["is_gc"]
    msp = has & (fam == "msp") & ~t["is_gc"]
    gc = has & t["is_gc"]
    mag = has & (fam == "magnetar")
    for m, kw in ((young, dict(s=3, c=ps.MUTED, alpha=0.5, lw=0, zorder=2,
                                 label=f"ordinary and radio-quiet pulsars ({young.sum():,})")),
                  (msp, dict(s=7, c=ps.BLUE, alpha=0.85, lw=0, zorder=3,
                             label=f"millisecond pulsars ({msp.sum():,})")),
                  (gc, dict(s=16, facecolors="none", edgecolors=ps.BLUE, marker="^",
                            lw=0.7, zorder=4,
                            label=f"pulsars in globular clusters ({gc.sum():,})")),
                  (mag, dict(s=18, c=ps.ORANGE, marker="D", lw=0.3, edgecolors=ps.SURFACE,
                             zorder=5, label=f"magnetars ({mag.sum():,})"))):
        x, y = _sky_xy(t.loc[m, "gl"], t.loc[m, "gb"])
        ax.scatter(x, y, **kw)
    # Where FAST's Galactic Plane Pulsar Snapshot survey looks (Han et al. 2021).
    box = dict(color=ps.INK_SECONDARY, lw=0.8, ls=(0, (3, 2)), zorder=6)
    _sky_box(ax, 30, 90, -10, 10, **box)
    _sky_box(ax, 145, 179.9, -10, 10, **box)
    _sky_box(ax, 180.1, 215, -10, 10, **box)
    ticks = np.arange(-150, 151, 30)
    ax.set_xticks(np.radians(ticks))
    ax.set_xticklabels([f"{(-d) % 360:.0f}°" for d in ticks], fontsize=7, color=ps.INK_SECONDARY)
    for label in ax.get_xticklabels():
        label.set_path_effects(HALO)
        label.set_zorder(8)
    ax.tick_params(axis="y", labelsize=7, colors=ps.MUTED)
    ax.grid(True, color=ps.GRID, lw=0.6)

    def sky_label(text, gl, gb, dx, dy, ha="left"):
        x, y = _sky_xy(gl, gb)
        ps.annotate(ax, text, (float(x), float(y)), (dx, dy), textcoords="offset points",
                    ha=ha, va="center", zorder=9)

    fast = h["msp_found_by_fast_gpps"]
    sky_label(f"FAST plane survey: found {fast['k']} of the\n{fast['n']} millisecond pulsars "
              f"outside clusters", 60, 10, 0, 74, ha="center")
    sky_label("Terzan 5\n(49 pulsars)", 3.84, 1.69, 40, 56)
    sky_label("47 Tuc (27 pulsars), with\nthe SMC behind it", 305.9, -44.9, 34, -24)
    sky_label("Large Magellanic\nCloud", 280.5, -32.9, 26, 34)
    sky_label("M15", 65.0, -27.3, -24, -16, ha="right")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.04), ncol=4, fontsize=7.5,
              handletextpad=0.2, columnspacing=1.2)
    ax.set_title("(a) The pulsar sky in Galactic coordinates: the young hug the Milky Way's "
                 "disc; millisecond pulsars spread wider", fontsize=10)

    # (b) height above the plane against characteristic age
    ok = t["in_ppdot"] & t["dist_ok"] & t["klass"].isin(NON_RECYCLED)
    okm = t["in_ppdot"] & t["dist_ok"] & (t["klass"] == "msp")
    z_floor = 1.0
    ax2.scatter(t.loc[ok, "tau_yr"], t.loc[ok, "z_pc"].clip(lower=z_floor), s=4,
                c=ps.MUTED, alpha=0.4, lw=0, zorder=2,
                label=f"ordinary and radio-quiet pulsars ({ok.sum():,})")
    ax2.scatter(t.loc[okm, "tau_yr"], t.loc[okm, "z_pc"].clip(lower=z_floor), s=6,
                c=ps.BLUE, alpha=0.6, lw=0, zorder=2,
                label=f"millisecond pulsars ({okm.sum():,})")
    rows = h["by_age"]
    centres = np.array([np.sqrt(r["tau_lo_yr"] * r["tau_hi_yr"]) for r in rows])
    med = np.array([r["median_z_pc"] for r in rows])
    q1 = np.array([r["iqr_z_pc"][0] for r in rows])
    q3 = np.array([r["iqr_z_pc"][1] for r in rows])
    ax2.errorbar(centres, med, yerr=[med - q1, q3 - med], color=ps.INK, lw=1.6,
                 marker="o", ms=6, mfc=ps.INK, mec=ps.SURFACE, mew=1.5, capsize=0,
                 elinewidth=1.0, zorder=5,
                 label="median height of ordinary pulsars per age bin (bar: middle 50%)")
    for c, m_, r in zip(centres, med, rows):
        ax2.annotate(f"{m_:.0f} pc\nn = {r['n']}", (c, m_), xytext=(9, -2),
                     textcoords="offset points", fontsize=7, color=ps.INK_SECONDARY,
                     va="top", zorder=6, path_effects=HALO)
    # The same medians with a different model of the Galaxy's electrons.
    alt = h["by_age_ne2001"]
    ax2.plot([np.sqrt(r["tau_lo_yr"] * r["tau_hi_yr"]) for r in alt],
             [r["median_z_pc"] for r in alt], color=ps.INK_SECONDARY, lw=1.0,
             ls=(0, (3, 2)), marker="o", ms=5, mfc=ps.SURFACE, mec=ps.INK_SECONDARY,
             mew=1.0, zorder=4, label="the same medians with NE2001 distances")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(3e2, 3e11)
    ax2.set_ylim(z_floor * 0.8, 3e4)
    ax2.set_xlabel("characteristic age τ = P/2Ṗ (years)")
    ax2.set_ylabel("height above or below the Galactic plane |z| (pc)")
    ax2.axvline(pp.AGE_OF_UNIVERSE_YR, color=ps.AXIS, lw=0.8, zorder=1)
    ax2.text(pp.AGE_OF_UNIVERSE_YR * 1.15, 1.5e4, "age of the\nUniverse", fontsize=7,
             color=ps.MUTED, va="top")
    ax2.legend(loc="upper left", fontsize=7.5)
    ax2.set_title("(b) Born in the disc, flung out of it: heights grow with age for "
                  "tens of millions of years, then level off", fontsize=10)
    n_floor = int(((t.loc[ok | okm, "z_pc"]) < z_floor).sum())
    ps.caption(fig, (
        f"{_provenance()} (a) Mollweide projection, Galactic longitude increasing to the "
        f"left; dashed outlines mark the sky covered by FAST's Galactic Plane Pulsar "
        f"Snapshot survey (|b| < 10°). The map shows where surveys looked as much as where "
        f"pulsars are: leaving out that one survey's discoveries, millisecond pulsars' "
        f"median latitude rises from {lat['msp']['median_abs_b_deg']:.1f}° to "
        f"{lat['msp_not_fast_gpps']['median_abs_b_deg']:.1f}°, against "
        f"{lat['non_recycled']['median_abs_b_deg']:.1f}° for ordinary pulsars. (b) "
        f"|z| = d·|sin b| from the catalogue's best distance (mostly dispersion measure "
        f"through the YMW16 electron-density model, uncertain by tens of percent per "
        f"pulsar). YMW16 gives {h['n_placeholder_distance']} of these pulsars a "
        f"placeholder 25 kpc and they are left out; the dashed line repeats the medians "
        f"with the older NE2001 model, which places them all. Magellanic Cloud pulsars "
        f"excluded; {n_floor} pulsars within {z_floor:.0f} pc of the plane are drawn on "
        f"the bottom edge. Low-latitude surveys cannot see far above the plane at small "
        f"distances, which holds down the oldest, faintest bin. Characteristic age is "
        f"only a rough clock (figure 05)."))
    return ps.save(fig, FIGURE_DIR / "03_galaxy")


# --- 04 Recycling -----------------------------------------------------------------------------


def figure_recycling(t, results) -> list:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.4),
                                   gridspec_kw=dict(width_ratios=[1, 1.2], wspace=0.34))
    bf = results["binary_fractions"]

    # (a) binary fraction against spin period
    rows = [r for r in bf["by_period"] if r["n"] > 0]
    centres = np.array([np.sqrt(r["P_lo_s"] * r["P_hi_s"]) for r in rows])
    frac = np.array([r["fraction"] for r in rows]) * 100
    lo = np.array([r["ci95"][0] for r in rows]) * 100
    hi = np.array([r["ci95"][1] for r in rows]) * 100
    ax1.axvspan(1e-3, MSP_MAX_PERIOD_S, color=ps.GRID, alpha=0.6, lw=0, zorder=0)
    ax1.text(np.sqrt(1e-3 * MSP_MAX_PERIOD_S), 104, "millisecond pulsars", ha="center",
             fontsize=8, color=ps.INK_SECONDARY)
    ax1.errorbar(centres, frac, yerr=[frac - lo, hi - frac], color=ps.BLUE, lw=0,
                 elinewidth=1.6, marker="o", ms=7, mec=ps.SURFACE, mew=1.5, zorder=4)
    for c, h_, r in zip(centres, hi, rows):
        ax1.annotate(f"{r['k']}/{r['n']}", (c, h_), xytext=(0, 4), textcoords="offset points",
                     ha="center", fontsize=7, color=ps.INK_SECONDARY)
    ax1.text(0.97, 0.62,
             f"millisecond pulsars: {bf['msp']['fraction']:.0%} in binaries\n"
             f"(95% interval {bf['msp']['ci95'][0]:.0%}–{bf['msp']['ci95'][1]:.0%})\n\n"
             f"ordinary pulsars: {bf['non_recycled']['fraction']:.1%}\n"
             f"({bf['non_recycled']['ci95'][0]:.1%}–{bf['non_recycled']['ci95'][1]:.1%})",
             transform=ax1.transAxes, ha="right", va="top", fontsize=8, color=ps.INK)
    ax1.set_xscale("log")
    ax1.set_xlim(1e-3, 100)
    ax1.set_ylim(0, 110)
    ax1.set_yticks(range(0, 101, 20))
    ax1.set_xlabel("spin period P (seconds)")
    ax1.set_ylabel("pulsars with a binary companion (%)")
    ax1.set_title("(a) How often a pulsar has a companion", fontsize=10)

    # (b) spin period by companion type
    comp = results["companions"]
    order = ["UL", "He", "MS", "CO", "NS"]
    names = {"UL": "ultra-light star\n(black widows)", "He": "helium\nwhite dwarf",
             "MS": "ordinary star\n(main sequence)", "CO": "carbon–oxygen\nwhite dwarf",
             "NS": "neutron star"}
    rng = np.random.default_rng(1)
    field = t["is_binary"] & ~t["is_gc"] & t["P"].notna()
    for i, code in enumerate(order):
        m = field & (t["companion"] == code)
        p = t.loc[m, "P"].to_numpy()
        y = i + rng.uniform(-0.2, 0.2, size=len(p))
        ax2.scatter(p, y, s=12, c=ps.BLUE, alpha=0.55, lw=0, zorder=3)
        medv = comp[code]["median_P_ms"] / 1e3
        ax2.plot([medv, medv], [i - 0.3, i + 0.3], color=ps.INK, lw=2, zorder=4,
                 solid_capstyle="round")
        ax2.text(medv * 1.12, i - 0.3, f"median {comp[code]['median_P_ms']:.1f} ms",
                 ha="left", va="bottom", fontsize=7.5, color=ps.INK, path_effects=HALO,
                 zorder=5)
    ax2.set_yticks(range(len(order)))
    ax2.set_yticklabels([f"{names[c]}\n({comp[c]['n']})" for c in order], fontsize=8,
                        color=ps.INK_SECONDARY)
    ax2.set_ylim(len(order) - 0.45, -0.6)
    ax2.set_xscale("log")
    ax2.set_xlim(1e-3, 10)
    ax2.set_xlabel("spin period P (seconds)")
    ax2.grid(axis="y", visible=False)
    # Two individual systems, one at each end of the story, labelled in the gaps
    # between rows.
    ms = order.index("MS")
    for jname, text, xy_text, relpos in (
            ("J1023+0038", "J1023+0038 flips between radio pulsar and accreting\n"
                           "X-ray binary: recycling caught in the act",
             (1.1e-3, ms - 0.62), (0.0, 0.0)),
            ("J1302-6350", "B1259−63: a young pulsar orbiting a massive star\n"
                           "that has not yet transferred any mass",
             (0.1, ms + 0.52), (0.5, 0.5))):
        r = t.loc[t["jname"] == jname].iloc[0]
        ax2.scatter([r["P"]], [ms], s=40, facecolors="none", edgecolors=ps.INK, lw=1.0,
                    zorder=6)
        ps.annotate(ax2, text, (r["P"], ms), xy_text, textcoords="data", ha="left",
                    va="center", fontsize=7, zorder=6, relpos=relpos)
    ax2.set_title("(b) Spin period by type of companion", fontsize=10)
    ps.caption(fig, (
        f"{_provenance()} (a) Rotation-powered pulsars with a timing solution (a measured "
        f"Ṗ), outside globular clusters; bars are 95% Wilson intervals. Counting pulsars "
        f"without a timing solution lowers the millisecond-pulsar figure to "
        f"{bf['msp_including_untimed']['fraction']:.0%}, because unsolved orbits count as "
        f"'no companion'. (b) Binaries outside globular clusters, by the catalogue's "
        f"companion type; bars mark medians. The ordering is the recycling story: a helium "
        f"white dwarf is what remains of a low-mass star that fed the pulsar for up to "
        f"billions of years; a carbon–oxygen white dwarf, of a heavier star with a shorter "
        f"feeding phase; a neutron-star companion, of a massive star that transferred mass "
        f"only briefly before exploding."))
    return ps.save(fig, FIGURE_DIR / "04_recycling")


# --- 05 Caveats ---------------------------------------------------------------------------------


def figure_caveats(t, results) -> list:
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 10),
                             gridspec_kw=dict(hspace=0.45, wspace=0.3))
    (ax1, ax2), (ax3, ax4) = axes
    ages = results["msp_ages"]

    # (a) the Shklovskii effect on millisecond pulsars, period in milliseconds
    m = t["in_ppdot"] & (t["klass"] == "msp") & t["Pdot_int"].notna()
    ylim = (3e-23, 3e-18)
    pos = m & (t["Pdot_int"] > 0)
    neg = m & (t["Pdot_int"] <= 0)
    p_ms = t["P"] * 1e3
    for i in t.index[pos]:
        ax1.plot([p_ms[i], p_ms[i]], [t.at[i, "Pdot"], t.at[i, "Pdot_int"]],
                 color=ps.AXIS, lw=0.7, zorder=2)
    for i in t.index[neg]:
        ax1.annotate("", xy=(p_ms[i], ylim[0] * 1.2), xytext=(p_ms[i], t.at[i, "Pdot"]),
                     arrowprops=dict(arrowstyle="-|>", color=ps.MUTED, lw=0.6,
                                     mutation_scale=6), zorder=2)
    ax1.scatter(p_ms[m], t.loc[m, "Pdot"], s=14, facecolors="none", edgecolors=ps.MUTED,
                lw=0.7, zorder=3, label="observed Ṗ")
    ax1.scatter(p_ms[pos], t.loc[pos, "Pdot_int"], s=12, c=ps.BLUE, lw=0, zorder=4,
                label="after removing the Shklovskii term")
    ex = results["shklovskii_example"]
    ps.annotate(ax1, f"{ex['jname']}, the nearest millisecond pulsar:\n"
                     f"{ex['kinematic_fraction']:.0%} of its apparent spin-down is its "
                     f"motion", (ex["P_ms"], ex["pdot_observed"]), (1.45, 8e-19),
                textcoords="data", ha="left", va="center", fontsize=7.5)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(1.3, 32)
    ax1.set_ylim(*ylim)
    ax1.set_xticks([2, 3, 5, 10, 20, 30])
    ax1.set_xticklabels(["2", "3", "5", "10", "20", "30"])
    ax1.minorticks_off()
    ax1.set_xlabel("spin period P (milliseconds)")
    ax1.set_ylabel("period derivative Ṗ")
    ax1.legend(loc="lower right", fontsize=7.5, frameon=True, facecolor=ps.SURFACE,
               edgecolor="none", framealpha=1.0)
    ax1.set_title("(a) Motion masquerades as spin-down", fontsize=10)

    # (b) globular clusters: signed Pdot
    # Every measured Pdot, of either sign -- the P-Pdot sample would be circular
    # here, since it keeps only the positive ones.
    field = ~t["is_gc"] & (t["klass"] == "msp") & t["Pdot"].notna()
    gcm = t["is_gc"] & t["Pdot"].notna() & (t["P"] < MSP_MAX_PERIOD_S)
    rng = np.random.default_rng(2)
    rows = ((field, ps.BLUE, "in the Galactic field"),
            (gcm, ps.ORANGE, "in globular clusters"))
    for i, (mask, color, _) in enumerate(rows):
        v = t.loc[mask, "Pdot"].to_numpy()
        ax2.scatter(v, i + rng.uniform(-0.25, 0.25, len(v)), s=12, c=color, alpha=0.7,
                    lw=0, zorder=3)
    ax2.set_xscale("symlog", linthresh=1e-22, linscale=0.6)
    ax2.set_xlim(-3e-18, 3e-18)
    ax2.set_xticks([-1e-18, -1e-20, -1e-22, 0, 1e-22, 1e-20, 1e-18])
    ax2.axvline(0, color=ps.INK_SECONDARY, lw=0.8, zorder=2)
    ax2.set_yticks([0, 1])
    labels = []
    for mask, _, where in rows:
        v = t.loc[mask, "Pdot"]
        labels.append(f"millisecond pulsars\n{where}\n{(v < 0).sum()} of {len(v)} "
                      f"negative")
    ax2.set_yticklabels(labels, fontsize=8, color=ps.INK_SECONDARY)
    ax2.set_ylim(1.6, -0.6)
    ax2.grid(axis="y", visible=False)
    ax2.set_xlabel("observed Ṗ (symmetric log scale)")
    ax2.text(-2e-18, -0.5, "← apparently spinning up", fontsize=7.5, color=ps.INK_SECONDARY,
             ha="left", va="top")
    ax2.text(2e-18, -0.5, "slowing down →", fontsize=7.5, color=ps.INK_SECONDARY,
             ha="right", va="top")
    ax2.set_title("(b) In a star cluster, gravity can fake a spin-up", fontsize=10)

    # (c) braking indices
    bi = results["braking_indices"]
    n = np.asarray(bi["values"])
    tau = np.asarray(bi["tau_yr"])
    phys = (n > 0) & (n < 3) & (tau < 1e4)
    ax3.scatter(tau[~phys], n[~phys], s=10, c=ps.MUTED, alpha=0.6, lw=0, zorder=3,
                label="all other measurements")
    ax3.scatter(tau[phys], n[phys], s=24, c=ps.BLUE, lw=0.4, edgecolors=ps.SURFACE,
                zorder=4, label="young pulsars (τ < 10 kyr) with 0 < n < 3")
    ax3.axhline(3, color=ps.INK_SECONDARY, lw=0.9, zorder=2)
    ax3.text(2e8, 3.8, "n = 3: pure magnetic dipole", fontsize=7.5, color=ps.INK_SECONDARY,
             va="bottom", ha="center", path_effects=HALO, zorder=5)
    ax3.set_xscale("log")
    ax3.set_yscale("symlog", linthresh=1.0, linscale=1.0)
    ax3.set_xlim(1.5e2, 3e10)
    ax3.set_ylim(-3e6, 3e6)
    ax3.set_xlabel("characteristic age τ (years)")
    ax3.set_ylabel("braking index n = f·f̈ / ḟ²  (symmetric log scale)")
    young = [r for r in bi["young"] if 0 < r["n"] < 3]
    crab = results["crab"]
    ax3.legend(loc="lower left", fontsize=7.5, frameon=True, facecolor=ps.SURFACE,
               edgecolor="none", framealpha=1.0)
    ax3.text(0.02, 0.98,
             f"{bi['abs_n_above_10']['fraction']:.0%} of {bi['n_measured']} measurements "
             f"have |n| > 10: timing noise\nand glitches, not braking. The {len(young)} "
             f"young ones span {min(r['n'] for r in young):.1f}–"
             f"{max(r['n'] for r in young):.2f}\n(the Crab: {crab['braking_index']:.3f}), "
             f"all below the dipole's 3.",
             transform=ax3.transAxes, fontsize=7.5, color=ps.INK, va="top", zorder=6,
             bbox=dict(facecolor=ps.SURFACE, edgecolor="none", alpha=0.9, pad=2))
    ax3.set_title("(c) Braking index: only the youngest pulsars give a clean value",
                  fontsize=10)

    # (d) the Crab's clock
    epoch = crab["epoch_year"]
    marks = (
        (epoch - crab["spin_down_age_measured_n_yr"],
         f"birth date implied by the\nmeasured n = {crab['braking_index']:.2f}", ps.MUTED, 1),
        (epoch - crab["characteristic_age_yr"], "birth date implied\nby τ = P/2Ṗ", ps.MUTED, -1),
        (1054.5, "actual birth: the supernova\nrecorded by Chinese astronomers",
         ps.BLUE, 1),
        (epoch, "timing epoch", ps.INK, -1),
    )
    ax4.axhline(0, color=ps.AXIS, lw=1.0, zorder=1)
    ax4.plot([1054.5, epoch], [0, 0], color=ps.BLUE, lw=3, solid_capstyle="round", zorder=2)
    ax4.text((1054.5 + epoch) / 2, 0.08, f"true age {crab['true_age_yr']:.0f} yr",
             ha="center", va="bottom", fontsize=7.5, color=ps.INK)
    for year, text, color, side in marks:
        ax4.scatter([year], [0], s=60, c=color, edgecolors=ps.SURFACE, lw=1.5, zorder=4)
        ax4.annotate(f"{text}\n{year:.0f} AD" if year < 1900 else f"{text}\n{year:.1f}",
                     (year, 0), xytext=(0, 30 * side), textcoords="offset points",
                     ha="center", va="bottom" if side > 0 else "top", fontsize=7.5,
                     color=ps.INK_SECONDARY,
                     arrowprops=dict(arrowstyle="-", color=ps.AXIS, lw=0.8))
    ax4.set_xlim(-150, 2200)
    ax4.set_ylim(-1, 1)
    ax4.set_yticks([])
    ax4.spines["left"].set_visible(False)
    ax4.grid(False)
    ax4.set_xlabel("calendar year (AD)")
    ax4.text(0.5, 0.02,
             f"To be {crab['true_age_yr']:.0f} years old at the epoch, the Crab must have been "
             f"born spinning at\nP₀ ≈ {crab['birth_period_ms_n3']:.0f} ms (if n = 3) or "
             f"{crab['birth_period_ms_measured_n']:.0f} ms (if n = {crab['braking_index']:.2f})"
             f" -- not much faster than its {crab['period_now_ms']:.1f} ms then.",
             transform=ax4.transAxes, ha="center", va="bottom", fontsize=7.5, color=ps.INK)
    ax4.set_title(f"(d) The Crab's clock runs fast: τ = {crab['characteristic_age_yr']:,.0f} yr, "
                  f"true age {crab['true_age_yr']:.0f} yr", fontsize=10)
    ps.caption(fig, (
        f"{_provenance()} (a) Shklovskii term P·μ²·d/c from the catalogue's proper motions "
        f"and distances, for the {ages['n_with_shklovskii']} millisecond pulsars that have "
        f"them; Galactic-acceleration terms are not removed. For "
        f"{ages['n_kinematic_exceeds_observed']} (arrows) the motion term exceeds the whole "
        f"observed Ṗ, a sign that their distance or an unmodelled acceleration is off. "
        f"(b) Cluster pulsars' observed Ṗ includes a·P/c, where a is the line-of-sight "
        f"acceleration in the cluster's gravity. (c) n from the catalogue's frequency "
        f"derivatives, rotation-powered pulsars outside clusters. (d) Ages at the epoch of "
        f"the catalogue's timing solution (Lyne et al. 2015); birth periods assume the "
        f"spin-down law held unchanged since birth."))
    return ps.save(fig, FIGURE_DIR / "05_caveats")


def make_all(t, results) -> list:
    written = []
    written += figure_ppdot(t, results)
    written += figure_spin_periods(t, results)
    written += figure_galaxy(t, results)
    written += figure_recycling(t, results)
    written += figure_caveats(t, results)
    written.append(interactive_ppdot(t, results))
    return written


# --- The interactive P-Pdot diagram ------------------------------------------------------------

#: Light and dark steps of the same hues (the dark ones are re-stepped for a dark
#: surface, not inverted -- see the palette notes in plotstyle.py).
THEMES = {
    "light": dict(surface=ps.SURFACE, ink=ps.INK, ink2=ps.INK_SECONDARY, muted=ps.MUTED,
                  grid=ps.GRID, axis=ps.AXIS, valley="rgba(225,224,217,0.7)",
                  ordinary=ps.MUTED, msp=ps.BLUE, quiet=ps.AQUA, magnetar=ps.ORANGE,
                  ring=ps.INK_SECONDARY),
    "dark": dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#898781",
                 grid="#2c2c2a", axis="#4a4a46", valley="rgba(56,56,53,0.7)",
                 ordinary="#898781", msp="#3987e5", quiet="#199e70", magnetar="#d95926",
                 ring="#c3c2b7"),
}

_COMPANION_NAMES = {"He": "helium white dwarf", "CO": "carbon-oxygen white dwarf",
                    "NS": "neutron star", "MS": "main-sequence star",
                    "UL": "ultra-light star or planet", "HeT": "helium white dwarf",
                    "ULT": "ultra-light star or planet"}


def _age_text(yr: float) -> str:
    for scale, unit in ((1e9, "Gyr"), (1e6, "Myr"), (1e3, "kyr")):
        if yr >= scale:
            return f"{yr / scale:.3g} {unit}"
    return f"{yr:.3g} yr"


def _hover(r, family_label: str) -> str:
    name = f"<b>PSR {r['jname']}</b>"
    if isinstance(r["bname"], str) and r["bname"]:
        name += f" ({r['bname']})"
    period = (f"{r['P'] * 1e3:.4g} ms" if r["P"] < 1 else f"{r['P']:.4g} s")
    lines = [name, family_label,
             f"P = {period}   Ṗ = {_sci(r['Pdot'], 2)}",
             f"age τ = {_age_text(r['tau_yr'])}   field B = {_sci(r['B_gauss'])} G",
             f"spin-down power Ė = {_sci(r['edot'])} erg/s"]
    if r["is_binary"]:
        orbit = f"orbit {r['pb_days']:.3g} d" if np.isfinite(r["pb_days"]) else "binary"
        comp = _COMPANION_NAMES.get(r["companion"], "")
        lines.append(orbit + (f", companion: {comp}" if comp else ""))
    if r["dist_ok"]:
        lines.append(f"distance ≈ {r['dist_kpc']:.2g} kpc")
    if r["jname"] in FAMOUS:
        spec = FAMOUS[r["jname"]]
        lines.append(f"<i>{spec['label']}{': ' + spec['note'] if spec['note'] else ''}</i>")
    return "<br>".join(lines)


def interactive_ppdot(t, results) -> str:
    """A self-contained HTML page: the P-Pdot diagram with every pulsar hoverable."""
    import plotly.graph_objects as go

    L = THEMES["light"]
    samp = t["in_ppdot"]
    fam = _family(t)
    P = np.logspace(np.log10(7e-4), np.log10(3e2), 300)
    fig = go.Figure()
    styles = []  # per trace: {theme: restyle dict}, so dark mode can swap colours

    def add(trace, light: dict, dark: dict):
        fig.add_trace(trace)
        styles.append({"light": light, "dark": dark})

    D = THEMES["dark"]
    dash = dict(width=1, dash="dash")
    add(go.Scatter(x=P, y=pp.death_line_pdot(P, "multipole"), mode="lines",
                   line=dict(color=L["muted"], **dash), hoverinfo="skip", showlegend=False),
        {"line.color": L["muted"]}, {"line.color": D["muted"]})
    add(go.Scatter(x=P, y=pp.death_line_pdot(P, "dipole"), mode="lines", fill="tonexty",
                   fillcolor=L["valley"], line=dict(color=L["muted"], **dash),
                   name="death valley", hoverinfo="skip"),
        {"line.color": L["muted"], "fillcolor": L["valley"]},
        {"line.color": D["muted"], "fillcolor": D["valley"]})
    for age in (1e3, 1e5, 1e7, 1e9):
        add(go.Scatter(x=P, y=pp.pdot_for_age(P, age), mode="lines",
                       line=dict(color=L["axis"], width=1), hoverinfo="skip",
                       showlegend=False),
            {"line.color": L["axis"]}, {"line.color": D["axis"]})
    for field in (1e8, 1e10, 1e12, 1e14):
        add(go.Scatter(x=P, y=pp.pdot_for_field(P, field), mode="lines",
                       line=dict(color=L["axis"], width=1), hoverinfo="skip",
                       showlegend=False),
            {"line.color": L["axis"]}, {"line.color": D["axis"]})

    classes = (("ordinary", "ordinary radio pulsars", "ordinary radio pulsar", "circle", 4,
                "ordinary"),
               ("msp", "millisecond pulsars", "millisecond pulsar", "circle", 6, "msp"),
               ("radio_quiet", "radio-quiet (γ-ray/X-ray only)",
                "radio-quiet pulsar, seen only in γ-rays or X-rays", "circle", 7, "quiet"),
               ("xins", "X-ray isolated neutron stars", "X-ray isolated neutron star",
                "square", 8, "quiet"),
               ("magnetar", "magnetars", "magnetar", "diamond", 9, "magnetar"))
    for key, legend, label, symbol, size, colour in classes:
        d = t.loc[samp & (fam == key)]
        add(go.Scatter(x=d["P"], y=d["Pdot"], mode="markers", name=f"{legend} ({len(d):,})",
                       marker=dict(color=L[colour], size=size, symbol=symbol,
                                   opacity=0.75 if key == "ordinary" else 0.95,
                                   line=dict(width=0)),
                       text=[_hover(r, label) for _, r in d.iterrows()],
                       hovertemplate="%{text}<extra></extra>"),
            {"marker.color": L[colour]}, {"marker.color": D[colour]})
    b = t.loc[samp & t["is_binary"]]
    add(go.Scatter(x=b["P"], y=b["Pdot"], mode="markers",
                   name=f"ringed: has a binary companion ({len(b):,})",
                   marker=dict(symbol="circle-open", size=11, color=L["ring"],
                               line=dict(width=1)), hoverinfo="skip"),
        {"marker.color": L["ring"]}, {"marker.color": D["ring"]})

    def axis(title, lo, hi):
        return dict(type="log", range=[np.log10(lo), np.log10(hi)], title=title,
                    exponentformat="power", showgrid=False, zeroline=False,
                    linecolor=L["axis"], tickcolor=L["axis"], ticks="outside",
                    tickfont=dict(color=L["muted"]), title_font=dict(color=L["ink2"]))

    fig.update_layout(
        height=760, margin=dict(l=70, r=20, t=20, b=60),
        paper_bgcolor=L["surface"], plot_bgcolor=L["surface"],
        font=dict(family="system-ui, -apple-system, 'Segoe UI', sans-serif", size=12,
                  color=L["ink"]),
        xaxis=axis("spin period P (seconds)", 7e-4, 3e2),
        yaxis=axis("period derivative Ṗ (seconds per second)", 1e-22, 1e-8),
        # (A filled trace would otherwise make plotly reverse the legend.)
        legend=dict(x=0.01, y=0.99, bgcolor=L["surface"], font=dict(size=11),
                    traceorder="normal"),
        hoverlabel=dict(bgcolor=L["surface"], bordercolor=L["axis"],
                        font=dict(color=L["ink"], size=12)),
        hovermode="closest",
    )
    theme_layout = {
        name: {"paper_bgcolor": th["surface"], "plot_bgcolor": th["surface"],
               "font.color": th["ink"], "legend.bgcolor": th["surface"],
               "xaxis.linecolor": th["axis"], "yaxis.linecolor": th["axis"],
               "xaxis.tickcolor": th["axis"], "yaxis.tickcolor": th["axis"],
               "xaxis.tickfont.color": th["muted"], "yaxis.tickfont.color": th["muted"],
               "xaxis.title.font.color": th["ink2"], "yaxis.title.font.color": th["ink2"],
               "hoverlabel.bgcolor": th["surface"], "hoverlabel.bordercolor": th["axis"],
               "hoverlabel.font.color": th["ink"]}
        for name, th in THEMES.items()}
    themes_js = json.dumps({name: {"layout": theme_layout[name],
                                   "traces": [s[name] for s in styles]}
                            for name in THEMES})
    plot_html = fig.to_html(full_html=False, include_plotlyjs=True, div_id="ppdot",
                            config={"displaylogo": False, "responsive": True})

    n = results["census"]
    bf = results["binary_fractions"]
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pulsar P–Ṗ diagram</title>
<style>
  :root {{ color-scheme: light; --surface: {L['surface']}; --ink: {L['ink']};
          --ink2: {L['ink2']}; --muted: {L['muted']}; }}
  @media (prefers-color-scheme: dark) {{
    :root {{ color-scheme: dark; --surface: {D['surface']}; --ink: {D['ink']};
            --ink2: {D['ink2']}; --muted: {D['muted']}; }}
  }}
  body {{ margin: 0; background: var(--surface); color: var(--ink);
         font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }}
  main {{ max-width: 1180px; margin: 0 auto; padding: 24px 16px 40px; }}
  h1 {{ font-size: 1.35rem; margin: 0 0 6px; }}
  p {{ margin: 0 0 10px; color: var(--ink2); max-width: 80ch; }}
  .note {{ font-size: 0.85rem; color: var(--muted); }}
</style>
</head>
<body>
<main>
<h1>The pulsar P–Ṗ diagram</h1>
<p>{n['n_in_ppdot']:,} neutron stars, placed by how fast they spin (P, across) and how fast
they slow down (Ṗ, up). Hover over any point for its name and the age, magnetic field and
spin-down power that follow from those two numbers. Drag to zoom; double-click to reset.</p>
<p>Guide lines: rising diagonals are constant characteristic age (1 kyr, 100 kyr, 10 Myr,
1 Gyr from top); falling diagonals are constant surface field (10⁸, 10¹⁰, 10¹², 10¹⁴ G from
bottom). The shaded band is the death valley, where radio emission is expected to stop.
{bf['msp']['fraction']:.0%} of millisecond pulsars are in binaries, against
{bf['non_recycled']['fraction']:.1%} of ordinary pulsars.</p>
{plot_html}
<p class="note">{_provenance()} Globular-cluster pulsars are left out: their Ṗ measures
the cluster's gravity as much as their spin-down. Ages, fields and powers assume
magnetic-dipole braking with I = 10⁴⁵ g cm² and R = 10 km. Made by run.py in
03-pulsar-ppdot.</p>
</main>
<script>
(function () {{
  var gd = document.getElementById("ppdot");
  var THEMES = {themes_js};
  function apply(dark) {{
    var th = THEMES[dark ? "dark" : "light"];
    th.traces.forEach(function (s, i) {{ Plotly.restyle(gd, s, [i]); }});
    Plotly.relayout(gd, th.layout);
  }}
  var mq = window.matchMedia("(prefers-color-scheme: dark)");
  if (mq.matches) apply(true);
  mq.addEventListener("change", function (e) {{ apply(e.matches); }});
}})();
</script>
</body>
</html>
"""
    out = FIGURE_DIR / "ppdot_interactive.html"
    out.write_text(page, encoding="utf-8")
    return str(out)
