"""Figure generation.

Every figure is drawn from what ``analysis.py`` returned, never by re-running
any of it, so what is plotted is exactly what was measured.

Design notes. Star density is drawn as a 2D histogram on a logarithmic,
single-hue scale: a scatter plot of a quarter of a million points is a solid
blob, and a log scale is what lets a lone giant and a bin holding 400 red dwarfs
share one picture. Empty bins are left transparent. Categorical colours are the
project's validated colourblind-safe slots, never more than three per panel.
Both "extra" axes on the HR diagram -- temperature on top, brightness relative
to the Sun on the right -- are unit conversions of the axis they sit opposite,
not second data scales.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import plotstyle as ps
from config import FIGURE_DIR, SUN_ABS_G, SUN_BP_RP

ps.use_style()

PROVENANCE = (
    "Data: ESA Gaia DR3 (gaiadr3.gaia_source), retrieved from the Gaia@AIP mirror of the "
    "release. Temperatures: Pecaut & Mamajek (2013), 2022 update."
)

#: Extent of the HR diagram panels: colour (x) and absolute magnitude (y).
HRD_XLIM = (-0.7, 5.4)
HRD_YLIM = (19.5, -3.8)   # inverted: bright at the top, as astronomers draw it


# --------------------------------------------------------------------------
# Shared pieces
# --------------------------------------------------------------------------

def density_image(ax, x, y, *, dx=0.02, dy=0.05, xlim=HRD_XLIM, ylim=HRD_YLIM,
                  vmax=None, cmap=None, zorder=1):
    """Star counts per (dx x dy) cell, log-scaled, empty cells transparent."""
    xe = np.arange(xlim[0], xlim[1] + dx, dx)
    ye = np.arange(min(ylim), max(ylim) + dy, dy)
    h, _, _ = np.histogram2d(x, y, bins=[xe, ye])
    h = np.ma.masked_less(h.T, 1)
    mesh = ax.pcolormesh(xe, ye, h, cmap=cmap or ps.density_cmap(),
                         norm=LogNorm(vmin=1, vmax=vmax or max(h.max(), 2)),
                         rasterized=True, zorder=zorder, linewidth=0)
    return mesh


def temperature_ticks(dwarf):
    """Colour -> temperature for main-sequence stars, from Mamajek's table.

    Only the part of the table where colour rises monotonically with spectral
    type is used (through M8V); beyond it the Gaia colour of ultracool dwarfs
    stops being a clean thermometer.
    """
    colour, teff = dwarf["Bp-Rp"], dwarf["Teff"]
    ok = np.isfinite(colour) & np.isfinite(teff)
    colour, teff = colour[ok], teff[ok]
    keep = np.concatenate([[True], np.diff(np.maximum.accumulate(colour)) > 0])
    colour, teff = colour[keep], teff[keep]
    wanted = [10000, 7000, 5000, 4000, 3500, 3000, 2700]
    ticks = [float(np.interp(t, teff[::-1], colour[::-1])) for t in wanted
             if teff.min() <= t <= teff.max()]
    labels = [f"{t:,}" for t in wanted if teff.min() <= t <= teff.max()]
    return ticks, labels


def spectral_class_edges(dwarf):
    """BP-RP at which each spectral class begins, for dwarfs (A0V, F0V, ...)."""
    edges = {}
    for cls in "AFGKM":
        i = np.where(dwarf["SpT"] == f"{cls}0V")[0]
        if len(i) and np.isfinite(dwarf["Bp-Rp"][i[0]]):
            edges[cls] = float(dwarf["Bp-Rp"][i[0]])
    return edges


def hrd_frame(ax, dwarf, *, xlim=HRD_XLIM, ylim=HRD_YLIM, top=True, right=True,
              classes=True):
    """Axes, limits and the two conversion axes of an HR diagram."""
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("Colour  BP − RP  (mag)   ← hotter · cooler →")
    ax.set_ylabel("Absolute magnitude  M$_G$   ← fainter · brighter →")

    if top and dwarf is not None:
        ticks, labels = temperature_ticks(dwarf)
        top_ax = ax.secondary_xaxis("top")
        top_ax.set_xticks(ticks, labels)
        top_ax.set_xlabel("Surface temperature of a main-sequence star (K)",
                          color=ps.INK_SECONDARY)
        top_ax.tick_params(colors=ps.MUTED)
        top_ax.spines["top"].set_visible(False)

    if right:
        lum = [1e4, 1e2, 1.0, 1e-2, 1e-4, 1e-5]
        names = ["10,000×", "100×", "1× Sun", "1/100", "1/10,000", "1/100,000"]
        mags = [SUN_ABS_G - 2.5 * np.log10(v) for v in lum]
        inside = [(m, n) for m, n in zip(mags, names) if min(ylim) <= m <= max(ylim)]
        right_ax = ax.secondary_yaxis("right")
        right_ax.set_yticks([m for m, _ in inside], [n for _, n in inside])
        right_ax.set_ylabel("Brightness relative to the Sun (G band)",
                            color=ps.INK_SECONDARY)
        right_ax.tick_params(colors=ps.MUTED)
        right_ax.spines["right"].set_visible(False)

    if classes and dwarf is not None:
        edges = spectral_class_edges(dwarf)
        # Letters sit just inside the top of the plot, centred on each class.
        bounds = [xlim[0]] + list(edges.values()) + [4.9]
        letters = ["B"] + list(edges.keys())
        y = ylim[1] + 0.02 * (ylim[0] - ylim[1])
        for letter, a, b in zip(letters, bounds[:-1], bounds[1:]):
            ax.text(0.5 * (a + b), y, letter, ha="center", va="top", fontsize=8.5,
                    color=ps.INK_SECONDARY, fontweight="bold")
        for x in edges.values():
            ax.plot([x, x], [ylim[1], ylim[1] + 0.045 * (ylim[0] - ylim[1])],
                    color=ps.AXIS, linewidth=0.8, zorder=2)


def mark_sun(ax, *, label=True, size=9):
    """The Sun, at its Gaia colour and absolute magnitude."""
    ax.plot([SUN_BP_RP], [SUN_ABS_G], "o", color=ps.C_SUN, markersize=size,
            markeredgecolor=ps.SURFACE, markeredgewidth=2, zorder=8)
    ax.plot([SUN_BP_RP], [SUN_ABS_G], ".", color=ps.SURFACE, markersize=2.5, zorder=9)
    if label:
        ps.annotate(ax, "The Sun", xy=(SUN_BP_RP, SUN_ABS_G),
                    xytext=(SUN_BP_RP - 0.95, SUN_ABS_G + 1.2), color=ps.INK,
                    ha="right", va="center", fontsize=8.5, fontweight="bold")


def colorbar(fig, mesh, ax, label):
    cb = fig.colorbar(mesh, ax=ax, pad=0.09, fraction=0.035, aspect=35)
    cb.set_label(label, color=ps.INK_SECONDARY)
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=ps.MUTED, labelsize=7.5)
    return cb


# --------------------------------------------------------------------------
# The sky: star counts over the whole catalogue
# --------------------------------------------------------------------------

#: Galactic (l, b) in degrees, resolved by name through CDS Sesame.
SKY_LABELS = [
    # (text, l, b, text offset in plot degrees (dx, dy), alignment)
    ("Galactic centre, behind\nthe bulge's dust", 0.0, 0.0, (-38, -34), "center"),
    ("Large Magellanic Cloud\n(a satellite galaxy, 50 kpc)", 280.47, -32.89, (38, -8), "left"),
    ("Small Magellanic Cloud", 302.81, -44.33, (34, -20), "left"),
    ("ω Centauri, a globular\ncluster of ~10 million stars", 309.10, 14.97, (42, 20), "left"),
    ("Andromeda Galaxy (M31)", 121.17, -21.57, (-18, -34), "center"),
    ("Dust clouds in Ophiuchus,\nhiding the stars behind them", 353.69, 17.69, (-40, 40), "center"),
]


def _plot_x(l_deg):
    """Galactic longitude -> Mollweide x in radians, l increasing leftwards."""
    return np.radians(-(((l_deg + 180.0) % 360.0) - 180.0))


def figure_sky(res, art):
    """Where Gaia's 1.8 billion stars are: a map from star counts alone."""
    sky = art["sky"]
    stats = res["sky"]
    fig = plt.figure(figsize=(11.5, 6.6))
    ax = fig.add_subplot(projection="mollweide")
    dens = sky["density"]
    norm = LogNorm(vmin=np.percentile(dens, 0.5), vmax=dens.max(), clip=True)
    cmap = ps.density_cmap("sky")
    mesh = ax.pcolormesh(np.radians(sky["x_edges"]), np.radians(sky["y_edges"]), dens,
                         cmap=cmap, norm=norm, rasterized=True, linewidth=0)

    # Graticule labels, relabelled for longitude increasing to the left.
    ax.set_xticks(np.radians([-120, -60, 0, 60, 120]))
    ax.set_xticklabels([])
    for l_deg in (120, 60, 0, 300, 240):
        ax.text(_plot_x(l_deg), np.radians(2.0), f"{l_deg}°", ha="center", va="bottom",
                fontsize=7.5, color=ps.SURFACE, fontweight="bold", zorder=6)
    ax.set_yticks(np.radians([-60, -30, 0, 30, 60]))
    ax.set_yticklabels(["−60°", "−30°", "b = 0°", "+30°", "+60°"], color=ps.MUTED)
    ax.grid(color=ps.SURFACE, linewidth=0.5, alpha=0.55)

    for text, l_deg, b_deg, (dx, dy), ha in SKY_LABELS:
        x, y = _plot_x(l_deg), np.radians(b_deg)
        ps.annotate(ax, text, xy=(x, y), xytext=(x + np.radians(dx), y + np.radians(dy)),
                    color=ps.INK, ha=ha, va="center", fontsize=8, linespacing=1.3)

    cb = fig.colorbar(mesh, ax=ax, orientation="horizontal", fraction=0.04, pad=0.06,
                      aspect=45, shrink=0.6)
    cb.set_label("Gaia sources per square degree", color=ps.INK_SECONDARY)
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=ps.MUTED, labelsize=7.5)

    ax.set_title(
        f"The Milky Way from star counts alone — {stats['n_catalogue'] / 1e9:.2f} billion "
        f"Gaia sources, Galactic coordinates", pad=14)
    ps.caption(
        fig,
        f"Data: ESA Gaia DR3 (gaiadr3.gaia_source_lite). Counted by the archive for a random "
        f"{stats['sample_fraction']:.1%} of the catalogue (random_index < "
        f"{stats['n_sampled']:,}) and scaled up; the fine speckle at high latitude is the "
        f"sampling noise of that subsample, not structure. Half of all Gaia sources lie within "
        f"{stats['half_of_stars_within_abs_b_deg']:.1f}° of the Galactic plane, a band covering "
        f"{stats['sky_fraction_within_that_band']:.0%} of the sky. Pale lanes along the plane "
        f"are dust blocking the light of the stars behind it. Galactic longitude increases to "
        f"the left, the convention for a sky seen from inside.",
    )
    return ps.save(fig, FIGURE_DIR / "05_sky_map")


# --------------------------------------------------------------------------
# 1. The HR diagram -- the keeper
# --------------------------------------------------------------------------

#: Named stars to label on the main diagram: text position in (colour, mag)
#: relative to the star, and a short note. Kept to a few that each teach
#: something; the rest are in the README table.
HERO_STARS = {
    "61 Cygni A": (0.75, -1.45, "first stellar distance, 1838"),
    "Barnard's Star": (-1.2, 0.55, None),
    "TRAPPIST-1": (-0.55, -2.9, "seven planets"),
    "40 Eridani B": (0.75, 0.5, "the first white dwarf known"),
    "van Maanen's Star": (1.1, -0.1, None),
}


def figure_hr_diagram(res, art):
    """The annotated HR diagram of every well-measured star within 100 pc."""
    import analysis as an

    d, kept, dwarf = art["d"], art["kept"], art["dwarf"]
    grid, ridge, _ = art["ridge"]
    x, y = d["bp_rp"][kept], d["abs_g"][kept]

    fig, ax = plt.subplots(figsize=(10.5, 11.0))
    mesh = density_image(ax, x, y, dx=0.015, dy=0.04)
    hrd_frame(ax, dwarf)
    mark_sun(ax)

    lbl = dict(color=ps.INK, fontsize=8.5, linespacing=1.35, ha="left", va="center")
    f = res["features"]

    def on_ridge(c, offset=0.0):
        return (c, float(an.ridge_at(np.array([c]), grid, ridge)[0]) + offset)

    # Upper left: the hot end of the main sequence.
    ps.annotate(ax, "Hot, massive, short-lived\n(A and F stars)", xy=on_ridge(0.2),
                xytext=(-0.62, -1.4), **lbl)
    # Upper right: the evolved stars, stacked in the empty corner.
    ps.annotate(ax, "Red giant branch", xy=(1.38, -1.0), xytext=(2.3, -2.6), **lbl)
    rc = f.get("red_clump", {})
    if rc.get("ok"):
        ps.annotate(ax, "Red clump: helium burning in the core\n"
                    f"M$_G$ = {rc['abs_g']:.2f} measured here",
                    xy=(rc["bp_rp"], rc["abs_g"]), xytext=(2.3, -0.5), **lbl)
    ps.annotate(ax, "Subgiants: hydrogen used up in the core,\nswelling towards giants",
                xy=(1.0, 3.2), xytext=(2.3, 1.7), **lbl)
    # Right: the main sequence and what lies along it.
    ps.annotate(ax, "MAIN SEQUENCE\nstars fusing hydrogen in their cores,\n"
                "where they spend ~90% of their lives", xy=on_ridge(1.15),
                xytext=(2.45, 4.2), **lbl)
    bs = f.get("binary_sequence", {})
    if bs.get("ok"):
        ps.annotate(ax, "Unresolved binaries: two stars in one dot,\n"
                    f"up to 0.75 mag brighter (measured {abs(bs['offset']):.2f})",
                    xy=on_ridge(1.75, bs["offset"]), xytext=(3.05, 7.0), **lbl)
    ps.annotate(ax, "Red dwarfs (M type):\nthe commonest stars in the Galaxy",
                xy=on_ridge(2.8), xytext=(3.5, 9.6), **lbl)
    ps.annotate(ax, "The faintest stars, close to the\nhydrogen-burning limit: below it,\n"
                "brown dwarfs", xy=(5.05, 15.95), xytext=(3.45, 17.6), **lbl)
    # Lower left: the white dwarfs.
    ps.annotate(ax, "WHITE DWARFS\nexposed cores of dead stars:\nno fusion, cooling for "
                "billions of years", xy=(0.3, 13.1), xytext=(-0.62, 17.2), **lbl)

    # Named stars: an ink ring so they read on top of any density.
    for name, (dx, dy, note) in HERO_STARS.items():
        rec = res["named_stars"].get(name, {})
        if not rec.get("kept"):
            continue
        sx, sy = rec["bp_rp"], rec["abs_g"]
        ax.plot([sx], [sy], "o", markersize=5.5, markerfacecolor="none",
                markeredgecolor=ps.INK, markeredgewidth=1.1, zorder=7)
        text = name if not note else f"{name}\n{note}"
        ps.annotate(ax, text, xy=(sx, sy), xytext=(sx + dx, sy + dy),
                    color=ps.INK_SECONDARY, ha="left" if dx > 0 else "right",
                    va="center", fontsize=7.5)

    # The colour bar sits in the empty space below the main sequence.
    cax = ax.inset_axes([0.36, 0.085, 0.24, 0.014])
    cb = fig.colorbar(mesh, cax=cax, orientation="horizontal")
    cb.set_label("stars per cell (0.015 × 0.04 mag), log scale", fontsize=7.5,
                 color=ps.INK_SECONDARY)
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=ps.MUTED, labelsize=7)

    ax.set_title(f"The Hertzsprung–Russell diagram of the solar neighbourhood — "
                 f"{res['n_final']:,} stars within 100 pc", pad=34)
    ps.caption(
        fig,
        f"{PROVENANCE}  Every star with Gaia parallax > 10 mas measured to better than 10%, "
        f"RUWE < {res['cuts']['max_ruwe']}, and clean BP/RP photometry (see figure 02). "
        f"No extinction correction: within 100 pc the Sun sits inside the Local Bubble, "
        f"nearly free of dust. The Sun, which Gaia cannot observe, is placed from synthetic "
        f"photometry (Casagrande & VandenBerg 2018). Gaia saturates on the very brightest "
        f"stars, so Sirius, Vega, Arcturus and a few others are absent; the faintest red and "
        f"white dwarfs are incomplete.",
    )
    return ps.save(fig, FIGURE_DIR / "01_hr_diagram")


# --------------------------------------------------------------------------
# 2. What the quality cuts remove
# --------------------------------------------------------------------------

#: Rejected stars, grouped by the first cut they failed (indices into
#: ``analysis.quality_cuts``). Three groups: the most a scatter panel can
#: carry with colours that stay distinct for colour-blind readers.
REJECT_GROUPS = [
    ("parallax measured to worse than 10%", (0,), ps.C_CUT_PARALLAX),
    ("RUWE ≥ 1.4: not a clean single-star fit", (1,), ps.C_CUT_RUWE),
    ("colour missing, contaminated, or too faint in BP", (2, 3, 4), ps.C_CUT_PHOT),
]


def context_cmap():
    """Grey ramp for 'the kept stars, for reference' underneath other data."""
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("context", [ps.GRID, ps.MUTED, ps.INK_SECONDARY])
    cmap.set_bad(alpha=0.0)
    cmap.set_under(alpha=0.0)
    return cmap


def sky_axes(ax, title):
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.grid(color=ps.GRID, linewidth=0.5)
    ax.set_title(title, fontsize=9.5, loc="left")


def sky_density_grid(l_deg, b_deg, *, level=4, width=720):
    """Stars per square degree from individual positions, on the sky-map grid.

    Counts go into equal-area HEALPix pixels in Galactic coordinates (level 4:
    3,072 pixels of 13.4 square degrees), then are sampled onto the same
    regular longitude-latitude grid as the all-sky map.
    """
    import astropy.units as u
    from astropy.coordinates import Galactic, SkyCoord
    from astropy_healpix import HEALPix

    hpx = HEALPix(nside=2 ** level, order="nested", frame=Galactic())
    pix = hpx.lonlat_to_healpix(np.asarray(l_deg) * u.deg, np.asarray(b_deg) * u.deg)
    counts = np.bincount(pix, minlength=hpx.npix).astype(float)
    density = counts / hpx.pixel_area.to_value(u.deg ** 2)
    x_edges = np.linspace(-180.0, 180.0, width + 1)
    y_edges = np.linspace(-90.0, 90.0, width // 2 + 1)
    X, Y = np.meshgrid(0.5 * (x_edges[:-1] + x_edges[1:]), 0.5 * (y_edges[:-1] + y_edges[1:]))
    grid = density[hpx.lonlat_to_healpix(np.mod(-X, 360.0) * u.deg, Y * u.deg)]
    return x_edges, y_edges, grid


def figure_quality_cuts(res, art):
    """Where the rejected stars would have gone, and where they are on the sky."""
    d, kept, ff = art["d"], art["kept"], art["first_fail"]
    fig = plt.figure(figsize=(14.5, 11.2))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.2, 0.8], hspace=0.34)
    gs = outer[0].subgridspec(1, 3, wspace=0.3)
    bottom = outer[1].subgridspec(1, 3, width_ratios=[0.85, 1.2, 1.2], wspace=0.22)
    xk, yk = d["bp_rp"][kept], d["abs_g"][kept]

    for j, (label, idx, colour) in enumerate(REJECT_GROUPS):
        ax = fig.add_subplot(gs[0, j])
        density_image(ax, xk, yk, dx=0.03, dy=0.08, cmap=context_cmap())
        rej = np.isin(ff, idx)
        placed = rej & np.isfinite(d["bp_rp"]) & np.isfinite(d["abs_g"])
        ax.plot(d["bp_rp"][placed], d["abs_g"][placed], ".", color=colour,
                markersize=1.3, alpha=0.3, rasterized=True, zorder=3, linestyle="none")
        hrd_frame(ax, None, top=False, right=False, classes=False)
        ax.set_xlabel("BP − RP")
        ax.set_ylabel("M$_G$" if j else "Absolute magnitude  M$_G$")
        n_rej, n_placed = int(rej.sum()), int(placed.sum())
        note = (f"; {n_rej - n_placed:,} without a colour are not drawn"
                if n_placed < n_rej else "")
        ax.set_title(f"Removed: {label}\n{n_rej:,} stars{note}", fontsize=9, loc="left")
        ax.plot([], [], "s", color=colour, markersize=6, label="removed by this cut")
        ax.plot([], [], "s", color=ps.MUTED, markersize=6, label="kept (all panels)")
        ax.legend(loc="lower left", fontsize=7.5, handletextpad=0.3)

    # --- the cut flow, as a waterfall ----------------------------------------
    ax = fig.add_subplot(bottom[0, 0])
    flow = res["cut_flow"]
    colour_of = {}
    for _, idx, colour in REJECT_GROUPS:
        for i in idx:
            colour_of[i] = colour
    y = np.arange(len(flow))
    ax.barh(y[0], flow[0]["remaining"], height=0.5, color=ps.MUTED, zorder=3)
    for i, step in enumerate(flow[1:], start=1):
        ax.barh(y[i], step["remaining"], height=0.5, color=ps.MUTED, zorder=3)
        ax.barh(y[i], step["removed"], left=step["remaining"], height=0.5,
                color=colour_of[i - 1], zorder=3)
        ax.annotate(f"−{step['removed']:,}", xy=(flow[i - 1]["remaining"], y[i]),
                    xytext=(4, 0), textcoords="offset points", va="center",
                    fontsize=7.5, color=ps.INK_SECONDARY)
    ax.annotate(f"{flow[-1]['remaining']:,} kept", xy=(flow[-1]["remaining"], y[-1]),
                xytext=(0, -15), textcoords="offset points", ha="left", va="center",
                fontsize=8, color=ps.INK, fontweight="bold")
    labels = ["archive query:\nparallax > 10 mas"] + [
        s["cut"].replace(" (|C*|", "\n(|C*|").replace(": not", ":\nnot") for s in flow[1:]]
    ax.set_yticks(y, labels, fontsize=7.5)
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.set_xlim(0, flow[0]["remaining"] * 1.28)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(
        lambda v, _: f"{v / 1e3:,.0f}k" if flow[0]["remaining"] > 20_000 else f"{v:,.0f}"))
    ax.set_xlabel("Stars remaining after each cut")
    ax.set_title("The cut flow", fontsize=9.5, loc="left")

    # --- the sky: kept vs rejected-for-parallax, on one density scale --------
    rej = ff == 0
    maps = [(kept, context_cmap(), f"Kept: {int(kept.sum()):,} stars, spread evenly\n"
                                   "over the sky, as nearby stars should be"),
            (rej, ps.density_cmap("rejected"), f"Removed for parallax: {int(rej.sum()):,},\n"
                                              "piled towards the Galactic centre")]
    grids = [sky_density_grid(d["l"][m], d["b"][m]) for m, _, _ in maps]
    positive = np.concatenate([g[2][g[2] > 0].ravel() for g in grids])
    norm = LogNorm(vmin=max(positive.min(), 0.05), vmax=positive.max(), clip=True)
    for k, ((m, cmap, title), (xe, ye, grid)) in enumerate(zip(maps, grids)):
        ax = fig.add_subplot(bottom[0, 1 + k], projection="mollweide")
        mesh = ax.pcolormesh(np.radians(xe), np.radians(ye), np.ma.masked_equal(grid, 0),
                             cmap=cmap, norm=norm, rasterized=True, linewidth=0)
        ax.set_anchor("N")
        sky_axes(ax, title)
        if k == 1:
            x0, y0 = _plot_x(280.47), np.radians(-32.89)   # the Large Magellanic Cloud
            ps.annotate(ax, "Magellanic Clouds, 50,000 pc away:\nspurious parallaxes in a crowd",
                        xy=(x0, y0), xytext=(x0 - np.radians(95), y0 - np.radians(38)),
                        color=ps.INK, fontsize=7.5, ha="center", va="center")
        cb = fig.colorbar(mesh, ax=ax, orientation="horizontal", fraction=0.05, pad=0.04,
                          aspect=30, shrink=0.8)
        cb.set_label("stars per square degree (same scale in both maps)", fontsize=7.5,
                     color=ps.INK_SECONDARY)
        cb.outline.set_visible(False)
        cb.ax.tick_params(colors=ps.MUTED, labelsize=7)

    fig.suptitle("What the quality cuts remove, and why it matters", x=0.125, ha="left",
                 y=0.955)
    ps.caption(
        fig,
        f"{PROVENANCE}  Every candidate is attributed to the first cut it fails, so the counts "
        f"add up. Upper row: grey is the final diagram; coloured points are where the removed "
        f"stars would have landed. Lower right: Galactic coordinates, centre at the middle, "
        f"longitude increasing to the left. A star really within 100 pc has no reason to lie "
        f"towards the Galactic centre; a spurious parallax in a crowded field does.",
    )
    return ps.save(fig, FIGURE_DIR / "02_quality_cuts")


# --------------------------------------------------------------------------
# 3. Kinematics: how fast the stars move
# --------------------------------------------------------------------------

def figure_kinematics(res, art):
    """The HR diagram coloured by speed, split by speed, and Parenago's curve."""
    from scipy.stats import binned_statistic_2d

    d, kept, classes, dwarf = art["d"], art["kept"], art["classes"], art["dwarf"]
    grid, ridge, _ = art["ridge"]
    kin = res["kinematics"]

    fig = plt.figure(figsize=(14.5, 12.0))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.3, 1.0], hspace=0.34, wspace=0.3)

    # --- (a) median tangential velocity per cell ------------------------------
    ax = fig.add_subplot(gs[0, :2])
    x, y, v = d["bp_rp"][kept], d["abs_g"][kept], d["v_tan"][kept]
    xe = np.arange(HRD_XLIM[0], HRD_XLIM[1] + 0.05, 0.05)
    ye = np.arange(HRD_YLIM[1], HRD_YLIM[0] + 0.15, 0.15)
    med = binned_statistic_2d(x, y, v, "median", bins=[xe, ye]).statistic
    cnt = binned_statistic_2d(x, y, v, "count", bins=[xe, ye]).statistic
    med = np.ma.masked_where(~(cnt >= 3), med).T
    cmap = ps.density_cmap("speed")
    cmap.set_over(ps.BLUE_RAMP[-1])
    vmax = 150
    mesh = ax.pcolormesh(xe, ye, med, cmap=cmap, vmin=0, vmax=vmax, rasterized=True,
                         linewidth=0)
    hrd_frame(ax, dwarf, right=False)
    mark_sun(ax, label=False)
    ps.annotate(ax, "Halo subdwarfs:\nold, metal-poor stars\npassing through the disk,\n"
                "fainter than disk stars\nof the same colour", xy=(1.75, 9.4),
                xytext=(-0.6, 7.9), color=ps.INK, fontsize=8, ha="left", va="center",
                linespacing=1.3)
    for name, (tx, ty) in {"Groombridge 1830": (1.9, 4.9), "Kapteyn's Star": (1.1, 17.6)}.items():
        rec = res["named_stars"].get(name, {})
        if not rec.get("kept"):
            continue
        ax.plot([rec["bp_rp"]], [rec["abs_g"]], "o", markersize=5.5, markerfacecolor="none",
                markeredgecolor=ps.INK, markeredgewidth=1.1, zorder=7)
        text = f"{name}: {rec['v_tan_km_s']:.0f} km/s across the sky"
        if rec.get("radial_velocity_km_s") is not None:
            text += (f",\n{abs(rec['radial_velocity_km_s']):.0f} km/s along the line of sight "
                     f"(Gaia RVS): {rec['total_speed_km_s']:.0f} km/s in all")
        ps.annotate(ax, text, xy=(rec["bp_rp"], rec["abs_g"]), xytext=(tx, ty),
                    color=ps.INK_SECONDARY, fontsize=7.5, ha="left", va="center")
    cb = fig.colorbar(mesh, ax=ax, pad=0.02, fraction=0.03, aspect=35, extend="max")
    cb.set_label("Median speed across the sky, v$_{tan}$ (km/s)", color=ps.INK_SECONDARY)
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=ps.MUTED, labelsize=7.5)
    ax.set_title("(a) The same diagram, coloured by how fast its stars move", pad=30)

    # --- (b) Parenago's discontinuity -----------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    rows = kin["parenago_rows"]
    if rows:
        c = np.array([r["colour"] for r in rows])
        s = np.array([r["S"] for r in rows])
        e = np.array([r["S_err"] for r in rows])
        ax.errorbar(c, s, yerr=e, fmt="o", color=ps.BLUE, markersize=5, elinewidth=1,
                    capsize=0, markeredgecolor=ps.SURFACE, markeredgewidth=1.2, zorder=4,
                    label="main-sequence stars, per colour bin")
        p = kin.get("parenago")
        if p:
            xb = p["break_colour"]
            ax.axvline(xb, color=ps.ORANGE, linewidth=1.6, zorder=3,
                       label=f"slope changes at BP−RP = {xb:.2f}")
            ax.axvline(p["literature_bp_rp"], color=ps.INK_SECONDARY, linewidth=1.0,
                       linestyle=(0, (1, 2)), zorder=2,
                       label=f"B−V = 0.6 (Binney+ 2000) = {p['literature_bp_rp']:.2f}")
        ax.legend(loc="lower right", fontsize=7.5)
    ax.set_xlabel("Colour  BP − RP  (mag)")
    ax.set_ylabel("Velocity dispersion S (km/s)")
    ax.set_title("(b) Parenago's discontinuity", pad=30)
    ax.text(0.30, 0.50, "Blue of the break every star is young:\nhot stars burn out before "
            "they age.\nRed of it, stars of every age mix,\nand the old ones have been "
            "stirred up.", transform=ax.transAxes, ha="left", va="center", fontsize=7.5,
            color=ps.INK_SECONDARY, linespacing=1.4)

    # --- (c-e) small multiples by velocity class ------------------------------
    for j, (name, m) in enumerate(classes.items()):
        ax = fig.add_subplot(gs[1, j])
        sel = kept & m
        density_image(ax, d["bp_rp"][sel], d["abs_g"][sel], dx=0.03, dy=0.08)
        ok = np.isfinite(ridge)
        ax.plot(grid[ok], ridge[ok], color=ps.INK_SECONDARY, linewidth=0.9, zorder=4,
                label="main-sequence ridge, all stars")
        hrd_frame(ax, None, top=False, right=False, classes=False)
        ax.set_xlabel("BP − RP")
        ax.set_ylabel("M$_G$")
        ax.set_title(f"({'cde'[j]}) {name}\n{int(sel.sum()):,} stars", fontsize=9.5, loc="left")
        if j == 0:
            ax.legend(loc="lower left", fontsize=7.5)

    fig.suptitle("Motion and age: the fast-moving stars are the old ones", x=0.125,
                 ha="left", y=0.96)
    ps.caption(
        fig,
        f"{PROVENANCE}  v$_{{tan}}$ = 4.74 μ/ϖ is the speed across the sky relative to the Sun; "
        f"the third component, along the line of sight, is not included, so a star moving "
        f"mostly towards or away from us looks slow. Velocity classes follow Gaia "
        f"Collaboration, Babusiaux et al. (2018). Each small panel is scaled to its own "
        f"density. S is the proper-motion velocity dispersion of Dehnen & Binney (1998), "
        f"computed for main-sequence stars within 0.5 mag of the ridge.",
    )
    return ps.save(fig, FIGURE_DIR / "03_kinematics")


# --------------------------------------------------------------------------
# 4. Features up close, each measured against a published value
# --------------------------------------------------------------------------

def _hist_panel(ax, edges, counts, colour):
    ax.stairs(counts, edges, fill=True, color=colour, alpha=0.12, zorder=2)
    ax.stairs(counts, edges, color=colour, linewidth=1.4, zorder=3)


def figure_features(res, art):
    """Four features of the diagram, each measured and checked."""
    import analysis as an

    d, kept = art["d"], art["kept"]
    f = res["features"]
    fig = plt.figure(figsize=(13.5, 12.0))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.26)

    # --- (a) the twin-binary sequence -------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    bs = f["binary_sequence"]
    if not bs.get("ok"):
        raise RuntimeError("binary-sequence measurement failed; see results.json")
    edges = np.array(bs["hist_edges"])
    kept_frac = np.array(bs["hist"]) / np.sum(bs["hist"])
    rej_frac = np.array(bs["rejected_hist"]) / np.sum(bs["rejected_hist"])
    ax.stairs(kept_frac, edges, color=ps.BLUE, linewidth=1.5, zorder=4,
              label=f"kept: {bs['n_stars']:,} stars")
    ax.stairs(rej_frac, edges, color=ps.C_CUT_RUWE, linewidth=1.5, zorder=3,
              label=f"removed by the RUWE cut: {bs['n_rejected']:,}")
    ax.set_yscale("log")
    ax.set_ylim(5e-5, 0.2)
    ax.axvline(bs["predicted"], color=ps.INK_SECONDARY, linewidth=1.0,
               linestyle=(0, (1, 2)), zorder=2)
    peak_y = kept_frac[np.searchsorted(edges, bs["offset"]) - 1]
    ax.plot([bs["offset"]], [peak_y * 1.9], "v", color=ps.INK, markersize=7, zorder=5)
    ps.annotate(ax, f"twin peak: {bs['offset']:+.2f} ± {bs['offset_err']:.2f} mag\n"
                f"two identical stars: {bs['predicted']:+.3f}",
                xy=(bs["offset"], peak_y * 2.3), xytext=(-1.55, 0.035), color=ps.INK,
                fontsize=7.5, ha="left", va="center", linespacing=1.4)
    ax.text(0.03, 0.97, "single stars (the ridge) at 0", transform=ax.transAxes,
            ha="left", va="top", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_xlim(-1.6, 0.8)
    ax.set_xlabel("M$_G$ − main-sequence ridge  (mag;  negative = brighter)")
    ax.set_ylabel("Fraction of stars per 0.02 mag (log scale)")
    ax.legend(loc="upper right", fontsize=7.5)
    ax.set_title(f"(a) Two stars in one dot — K and early-M dwarfs, BP−RP "
                 f"{bs['colour_range'][0]}–{bs['colour_range'][1]}", loc="left", fontsize=9.5)

    # --- (b) the red clump --------------------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    rc = f["red_clump"]
    if rc.get("ok"):
        sel = (kept & (d["bp_rp"] > 0.95) & (d["bp_rp"] < 1.55)
               & (d["abs_g"] > -1.5) & (d["abs_g"] < 2.5))
        e = np.arange(-1.5, 2.5 + 0.1, 0.1)
        h, _ = np.histogram(d["abs_g"][sel], bins=e)
        _hist_panel(ax, e, h, ps.BLUE)
        xs = np.linspace(-1.5, 2.5, 400)
        ax.plot(xs, an._gauss_line(xs, *rc["fit"]), color=ps.INK, linewidth=1.6, zorder=5,
                label="Gaussian + sloping background fit")
        # Comparison strip: our value and the published ones, with 1-sigma bars.
        top = h.max()
        rows = [("this work", f"{rc['abs_g']:.2f} ± {rc['abs_g_err']:.2f}", rc["abs_g"],
                 rc["abs_g_err"], ps.INK)] + [
            (r["ref"], f"{r['abs_g']:.3f} ± {r['err']:.3f}", r["abs_g"], r["err"],
             ps.INK_SECONDARY) for r in rc["literature"]]
        for j, (lab, txt, val, err, col) in enumerate(rows):
            yv = top * (1.55 - 0.13 * j)
            ax.errorbar([val], [yv], xerr=[err], fmt="o", color=col, markersize=4.5,
                        elinewidth=1.4, capsize=0, zorder=6)
            ax.text(val - 0.12, yv, f"{lab}: {txt}", ha="left", va="center",
                    fontsize=7.5, color=col, transform=ax.transData)
        ax.set_ylim(0, top * 1.75)
        ax.legend(loc="upper left", fontsize=7.5)
        ax.set_title(f"(b) The red clump, a standard candle: {rc['n_giants']} giants, "
                     f"BP−RP 0.95–1.55", loc="left", fontsize=9.5)
    ax.set_xlabel("Absolute magnitude  M$_G$")
    ax.set_ylabel("Giants per 0.1 mag")
    ax.set_xlim(2.5, -1.5)

    # --- (c) the Jao gap ----------------------------------------------------------
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from scipy import ndimage

    jg = f["jao_gap"]
    sub = gs[1, 0].subgridspec(2, 1, height_ratios=[1.25, 1], hspace=0.42)
    ax = fig.add_subplot(sub[0])
    xe = np.arange(2.20, 2.70 + 0.005, 0.01)
    ye = np.arange(9.55, 10.75 + 0.005, 0.01)
    h2, _, _ = np.histogram2d(d["bp_rp"][kept], d["abs_g"][kept], bins=[xe, ye])
    # Light smoothing, a little more along colour (the gap runs nearly across
    # the sequence) than along M_G, so noise averages out along the gap.
    smooth = ndimage.gaussian_filter(h2, (2.0, 1.2))
    # Expected density: the same map smoothed heavily along M_G, which erases
    # any feature as narrow as the gap but keeps the main sequence itself.
    expect = ndimage.gaussian_filter(h2, (2.0, 12))
    ratio = np.ma.masked_where(expect < 0.25 * expect.max(), smooth / np.maximum(expect, 1e-9))
    diverging = LinearSegmentedColormap.from_list(
        "fewer-more", [ps.RED, "#f0efec", ps.BLUE])
    diverging.set_bad(alpha=0.0)
    mesh = ax.pcolormesh(xe, ye, ratio.T, cmap=diverging, norm=Normalize(0.6, 1.4),
                         rasterized=True, linewidth=0)
    xs = np.array([2.22, 2.68])
    ax.plot(xs, an.jao_dr2_gap(xs) + 0.08, color=ps.INK, linewidth=0.9,
            linestyle=(0, (1, 2)))
    ax.plot(xs, an.jao_dr2_gap(xs) - 0.08, color=ps.INK, linewidth=0.9,
            linestyle=(0, (1, 2)))
    ax.text(2.69, 9.60, "between the dotted lines: the gap's\npublished (DR2) position, "
            "±0.08 mag", fontsize=7.5, color=ps.INK, ha="right", va="top")
    ax.set_xlim(xe[0], xe[-1])
    ax.set_ylim(ye[-1], ye[0])
    ax.set_ylabel("M$_G$")
    ax.tick_params(labelbottom=True)
    ax.set_xlabel("BP − RP   (red: fewer stars than the smooth expectation)", fontsize=8)
    ax.set_title("(c) The Jao gap: where red dwarfs become fully convective", loc="left",
                 fontsize=9.5)
    cax = ax.inset_axes([0.03, 0.16, 0.28, 0.035])
    cb = fig.colorbar(mesh, cax=cax, orientation="horizontal", ticks=[0.6, 1.0, 1.4])
    cb.ax.set_xticklabels(["−40%", "as expected", "+40%"], fontsize=6.5, color=ps.MUTED)
    cb.outline.set_visible(False)

    ax = fig.add_subplot(sub[1])
    c = np.array(jg["centres"])
    h = np.array(jg["hist"])
    base = np.array(jg["baseline"])
    bw = c[1] - c[0]
    ax.stairs(h, np.append(c - bw / 2, c[-1] + bw / 2), color=ps.BLUE, linewidth=1.0,
              zorder=3, label="stars per 0.01 mag, along the gap")
    ax.plot(c, base, color=ps.INK, linewidth=1.4, zorder=4, label="smooth expectation")
    ax.axvspan(jg["centre"] - 0.025, jg["centre"] + 0.025, color=ps.RED, alpha=0.15,
               zorder=1, linewidth=0)
    ax.text(jg["centre"] + 0.05, base.max() * 1.05,
            f"{jg['deficit']:.0%} fewer stars, {jg['significance']:.1f}σ\n"
            f"(controls: ≤ {jg['max_control_significance']:.1f}σ)",
            fontsize=7.5, color=ps.INK, ha="left", va="top")
    ax.set_xlim(c[0], c[-1])
    ax.set_ylim(0, base.max() * 1.35)
    ax.set_xlabel(f"M$_G$ − {jg['slope']} × (BP−RP − {jg['pivot']}):  the diagram "
                  f"straightened along the gap", fontsize=8)
    ax.set_ylabel("Stars")
    ax.legend(loc="upper left", fontsize=7)

    # --- (d) the white dwarf sequence ---------------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    wdm = art["regions"]["white_dwarfs"] & kept
    density_image(ax, d["bp_rp"][wdm], d["abs_g"][wdm], dx=0.01, dy=0.025,
                  xlim=(-0.6, 1.4), ylim=(16.0, 10.0))
    wb = f["white_dwarf_branches"]
    row = min(wb["slices"], key=lambda r: abs(r["abs_g_lo"] - 12.8))
    ymid = 0.5 * (row["abs_g_lo"] + row["abs_g_hi"])
    ps.annotate(ax, "A: hydrogen atmospheres,\nthe 0.6 M$_\\odot$ track", xy=(row["a_colour"], ymid),
                xytext=(row["a_colour"] + 0.35, ymid - 1.3), color=ps.INK, fontsize=7.5,
                ha="left", va="center")
    ps.annotate(ax, "B: linked to helium\natmospheres", xy=(row["b_colour"], ymid),
                xytext=(row["b_colour"] - 0.45, ymid + 1.3), color=ps.INK, fontsize=7.5,
                ha="center", va="center")
    ax.set_xlim(-0.6, 1.4)
    ax.set_ylim(16.0, 10.0)
    ax.set_xlabel("Colour  BP − RP  (mag)")
    ax.set_ylabel("Absolute magnitude  M$_G$")
    ax.set_title(f"(d) {int(wdm.sum()):,} white dwarfs: two cooling sequences, "
                 f"{wb['median_separation']:.2f} mag apart", loc="left", fontsize=9.5)

    fig.suptitle("Four features of the diagram, each measured against a published value",
                 x=0.125, ha="left", y=0.955)
    ps.caption(
        fig,
        f"{PROVENANCE}  (a) Offsets from the single-star ridge: the RUWE cut removes most "
        f"unresolved binaries, but near-equal twins, whose light centre barely wobbles, "
        f"survive it. (b) Published calibrations differ by 0.055 mag; the clump's brightness "
        f"depends on the age and metallicity mix. (c) Top: stars counted in 0.01 mag cells and "
        f"divided by the same counts smoothed along M$_G$. Bottom: the published DR2 slope of "
        f"the gap is adopted and only its position fitted; the same search in five control "
        f"colour ranges calibrates what noise alone produces. (d) Branch labels follow Gaia "
        f"Collaboration, Babusiaux et al. (2018).",
    )
    return ps.save(fig, FIGURE_DIR / "04_features")


def ndimage_smooth(y, sigma=1.0):
    from scipy import ndimage

    return ndimage.gaussian_filter1d(np.asarray(y, dtype=float), sigma)


def make_all(res, art):
    written = []
    written += figure_hr_diagram(res, art)
    written += figure_quality_cuts(res, art)
    written += figure_kinematics(res, art)
    written += figure_features(res, art)
    if "sky" in art:
        written += figure_sky(res, art)
    return written
