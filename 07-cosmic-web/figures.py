"""Figure generation.

Every figure is drawn from what ``analysis.analyse`` returned, never by
re-running any of it, so what is plotted is exactly what was measured.

Design notes: galaxies are drawn in the colour they are named for -- red-
sequence galaxies red, blue-cloud galaxies blue, a pair that stays distinct
under every colour-vision deficiency (checked with a validator, not by eye).
When colour is not the point, galaxies are neutral ink, so that two panels
being compared differ only in what they show. Point clouds are rasterised
inside the SVGs; everything else stays vector. Each figure carries a caption
with its data provenance and the caveat a reader needs.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, ConnectionPatch, Polygon

import cosmo_physics as cp
import literature as lit
import plotstyle as ps
from config import (
    DENSITY_DV_KMS,
    EVOLUTION_Q,
    FIGURE_DIR,
    MGS_R_LIMIT,
    SLICE_DEC_HALF_WIDTH,
    SLICE_RA_RANGE,
    VL_MAG_LIMIT_H,
    VL_Z_MAX,
    VL_Z_MIN,
    WEDGE_Z_MAX,
)

ps.use_style()

PROVENANCE = ("Data: Sloan Digital Sky Survey DR17, Main Galaxy Sample (r < 17.77), via SkyServer. "
              "Distances: comoving, Planck 2018 cosmology (H0 = 67.7 km/s/Mpc).")

#: Callout text over a point cloud gets a surface-coloured backing so it stays
#: legible without hiding what is underneath entirely.
BACKING = {"boxstyle": "round,pad=0.25", "fc": ps.SURFACE, "ec": "none", "alpha": 0.85}

#: Half-sizes of the finger-of-God zoom window, Mpc: across and along the line of sight.
ZOOM_HALF_WIDTH = 14.0
ZOOM_HALF_DEPTH = 36.0


# ---------------------------------------------------------------------------
# The wedge: a slice of the universe seen from above
# ---------------------------------------------------------------------------


def _plane(ra, d, ra_mid):
    """Slice coordinates: the observer at the origin, the slice's central RA
    pointing up, RA increasing to the left as on the sky seen from inside."""
    ang = np.radians(np.asarray(ra, dtype=float) - ra_mid)
    d = np.asarray(d, dtype=float)
    return -d * np.sin(ang), d * np.cos(ang)


def _arc(ax, d, ra_mid, ra_lo=SLICE_RA_RANGE[0], ra_hi=SLICE_RA_RANGE[1], **kw):
    th = np.linspace(ra_lo, ra_hi, 400)
    x, y = _plane(th, np.full_like(th, d), ra_mid)
    ax.plot(x, y, **kw)


def wedge_frame(ax, ra_mid: float, d_max: float, rings_z, *, ring_labels: bool = True,
                ra_labels: bool = True, scalebar: bool = True, observer: bool = True,
                label_size: float = 7.5, right_margin: float = 1.16) -> None:
    """Rings of constant redshift, spokes of constant right ascension, a scale bar."""
    lo, hi = SLICE_RA_RANGE
    for z in rings_z:
        d = float(cp.comoving_distance(z))
        _arc(ax, d, ra_mid, color=ps.GRID, lw=0.7, zorder=0)
        if ring_labels:
            x, y = _plane(lo - 1.2, d, ra_mid)
            ax.text(x, y, f"z = {z:g}\n{d:,.0f} Mpc", fontsize=label_size, color=ps.MUTED,
                    ha="left", va="center", linespacing=1.15)
    _arc(ax, d_max, ra_mid, color=ps.AXIS, lw=0.8, zorder=0)
    for ra in (lo, hi):
        x, y = _plane([ra, ra], [0.0, d_max], ra_mid)
        ax.plot(x, y, color=ps.AXIS, lw=0.8, zorder=0)
    if ra_labels:
        for hour in range(int(np.ceil(lo / 15.0)), int(hi // 15) + 1):
            ra = hour * 15.0
            x0, y0 = _plane(ra, d_max, ra_mid)
            x1, y1 = _plane(ra, d_max * 1.03, ra_mid)
            ax.plot([x0, x1], [y0, y1], color=ps.AXIS, lw=0.8)
            xt, yt = _plane(ra, d_max * 1.07, ra_mid)
            ax.text(xt, yt, f"{hour}h", fontsize=label_size, color=ps.MUTED,
                    ha="center", va="center")
        xt, yt = _plane(ra_mid, d_max * 1.14, ra_mid)
        ax.text(xt, yt, "right ascension", fontsize=label_size, color=ps.MUTED,
                ha="center", va="center")
    if observer:
        ax.plot(0, 0, marker="o", ms=4, color=ps.INK, zorder=5)
        ax.text(0, -0.03 * d_max, "Milky Way\n(we are here)", fontsize=label_size,
                color=ps.INK_SECONDARY, ha="center", va="top", linespacing=1.1)
    if scalebar:
        x0, y0 = -0.93 * d_max, 0.03 * d_max
        ax.plot([x0, x0 + 100], [y0, y0], color=ps.INK_SECONDARY, lw=1.6,
                solid_capstyle="butt")
        ax.text(x0 + 50, y0 + 0.018 * d_max, "100 Mpc", fontsize=label_size,
                color=ps.INK_SECONDARY, ha="center", va="bottom")
        ax.text(x0 + 50, y0 - 0.018 * d_max,
                f"= {100 * cp.MPC_IN_LY / 1e6:.0f} million\nlight-years",
                fontsize=label_size - 0.5, color=ps.MUTED, ha="center", va="top",
                linespacing=1.1)
    half = np.radians(0.5 * (hi - lo))
    ax.set_xlim(-d_max * np.sin(half) * 1.08, d_max * np.sin(half) * right_margin)
    ax.set_ylim(-0.1 * d_max, d_max * 1.18)
    ax.set_aspect("equal")
    ax.axis("off")


def _galaxy_colors(red):
    rgba = np.empty((len(red), 4))
    rgba[:] = plt.matplotlib.colors.to_rgba(ps.C_BLUE_GALAXY)
    rgba[red] = plt.matplotlib.colors.to_rgba(ps.C_RED_GALAXY)
    return rgba


def _shuffled(n, seed=7):
    """Draw order for point clouds: random, so neither colour hides the other."""
    return np.random.default_rng(seed).permutation(n)


def _blue_cmap():
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("blue_ramp", ps.BLUE_RAMP[1:])
    cmap.set_bad(ps.SURFACE)
    return cmap


def _lognorm(vmin=1, vmax=None):
    from matplotlib.colors import LogNorm

    return LogNorm(vmin=vmin, vmax=vmax)


def _colorbar(fig, mappable, ax, label, **kw):
    cb = fig.colorbar(mappable, ax=ax, pad=0.015, fraction=0.04, **kw)
    cb.set_label(label, color=ps.INK_SECONDARY)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=7, colors=ps.MUTED)
    return cb


# ---------------------------------------------------------------------------
# 0. Where SDSS looked
# ---------------------------------------------------------------------------


def figure_footprint(res, art):
    """The survey on the sky, with the slice picked out and the Milky Way's
    disc -- the zone the survey could not see through -- marked."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    g = art["g"]
    rng = np.random.default_rng(3)
    pick = rng.choice(len(g["ra"]), size=min(120_000, len(g["ra"])), replace=False)

    def lon(ra):  # RA 12h at the centre, increasing to the left, as on the sky
        return np.radians(180.0 - np.asarray(ra, dtype=float))

    fig = plt.figure(figsize=(10, 5.4))
    ax = fig.add_subplot(projection="mollweide")
    ax.scatter(lon(g["ra"][pick]), np.radians(g["dec"][pick]), s=0.15, color=ps.MUTED,
               linewidths=0, alpha=0.6, rasterized=True)
    sl = art["slice"]
    ax.scatter(lon(g["ra"][sl]), np.radians(g["dec"][sl]), s=0.12, color=ps.BLUE,
               linewidths=0, rasterized=True)

    # The Galactic plane, split where it wraps around the map's edge.
    lgal = np.linspace(0, 360, 1441)
    c = SkyCoord(l=lgal * u.deg, b=np.zeros_like(lgal) * u.deg, frame="galactic").icrs
    x, y = lon(c.ra.deg), np.radians(c.dec.deg)
    x = np.where(x <= -np.pi, x + 2 * np.pi, x)
    breaks = np.flatnonzero(np.abs(np.diff(x)) > np.pi) + 1
    for seg in np.split(np.arange(len(x)), breaks):
        ax.plot(x[seg], y[seg], color=ps.ORANGE, lw=1.2)
    ax.text(lon(88), np.radians(-22), "the Milky Way's disc:\nits dust hides\nwhat lies behind",
            color=ps.ORANGE, fontsize=8, ha="center", va="center")
    ps.annotate(ax, f"the slice mapped in Fig. 01:\n{2 * SLICE_DEC_HALF_WIDTH:g}° thick, "
                "RA 8.5h–16.8h", xy=(lon(222), np.radians(-1.5)),
                xytext=(lon(262), np.radians(-35)), color=ps.BLUE, fontsize=8, ha="center")
    ax.set_xticks(lon(np.arange(30, 360, 30)))
    ax.set_xticklabels([])
    for ra in np.arange(60, 360, 60):
        ax.text(lon(ra), np.radians(-52), f"{ra // 15}h", color=ps.MUTED, fontsize=7,
                ha="center", va="center")
    ax.tick_params(axis="y", labelsize=7, colors=ps.MUTED)
    ax.grid(True, color=ps.GRID, lw=0.6)
    area = res["mask"]["surveyed_deg2"]
    ax.set_title(f"Where SDSS looked: {res['n_galaxies']:,} galaxies with measured redshifts "
                 f"over {area:,.0f} square degrees", fontsize=11, pad=14)
    ps.caption(fig, "Equatorial coordinates: right ascension across (increasing to the left, as "
               "on the sky), declination up. Grey: a random 120,000 of the main-sample galaxies; "
               "blue: the equatorial slice. The survey looks out of the Milky Way's disc, above "
               "its north pole (the big patch) and through three stripes in the south. "
               + PROVENANCE, y=0.07)
    return ps.save(fig, FIGURE_DIR / "00_footprint")


# ---------------------------------------------------------------------------
# 1. The cosmic web
# ---------------------------------------------------------------------------


def _fog_cluster(art):
    """The slice cluster to zoom in on: the richest named one that is clearly
    stretched along the line of sight."""
    cands = [c for c in art["groups_in_slice"] if c["elongation"] > 3.5 and c["z"] < 0.1]
    named = [c for c in cands if c["abell"]]
    return (named or cands or art["groups_in_slice"])[0]


def _finger_panel(ax, art, cluster, half_width=ZOOM_HALF_WIDTH, half_depth=ZOOM_HALF_DEPTH,
                  slab=4.0):
    """Galaxies around a cluster: across vs along the line of sight, in Mpc.

    Uses every galaxy of the catalogue (not only the slice) within +/-``slab``
    Mpc of the cluster's centre in declination, so the cluster is not cut in
    half by the slice's edge. The line of sight is vertical, the Milky Way
    straight down, and right ascension increases to the left as in the map.
    """
    g = art["g"]
    d_c = cluster["distance_mpc"]
    dra = ((g["ra"] - cluster["ra"] + 180.0) % 360.0) - 180.0
    across = -d_c * np.radians(dra) * np.cos(np.radians(cluster["dec"]))
    updown = d_c * np.radians(g["dec"] - cluster["dec"])
    along = g["d_mpc"] - d_c
    m = (np.abs(across) < half_width) & (np.abs(along) < half_depth) & (np.abs(updown) < slab)
    order = _shuffled(int(m.sum()))
    ax.scatter(across[m][order], along[m][order], s=7, c=_galaxy_colors(art["red"][m])[order],
               linewidths=0, alpha=0.9, rasterized=True)
    ax.set_xlim(-half_width, half_width)
    ax.set_ylim(-half_depth, half_depth)
    ax.set_aspect("equal")
    ax.set_xlabel("across the line of sight (Mpc)")
    ax.set_ylabel("along the line of sight: distance from its redshift (Mpc)")
    ax.annotate("", xy=(0.5, 0.015), xytext=(0.5, 0.075), xycoords="axes fraction",
                arrowprops={"arrowstyle": "-|>", "color": ps.INK_SECONDARY, "lw": 1.0})
    ax.text(0.53, 0.045, "toward\nthe Milky Way", transform=ax.transAxes, fontsize=7,
            color=ps.INK_SECONDARY, va="center")
    return int(m.sum())


def figure_wedge(res, art):
    """The flagship: every galaxy in a 4-degree slice, coloured by type, with a
    zoom onto one finger of God."""
    g, sl, red = art["g"], art["slice"], art["red"]
    ra_mid = art["ra_mid"]
    d_max = float(cp.comoving_distance(WEDGE_Z_MAX))
    x, y = _plane(g["ra"][sl], g["d_mpc"][sl], ra_mid)
    order = _shuffled(len(x))
    cols = _galaxy_colors(red[sl])

    fig = plt.figure(figsize=(15.5, 8.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[4.2, 1], wspace=0.08)
    ax = fig.add_subplot(gs[0])
    ax.scatter(x[order], y[order], s=1.0, c=cols[order], linewidths=0, alpha=0.85,
               rasterized=True, zorder=2)
    wedge_frame(ax, ra_mid, d_max, (0.05, 0.10, 0.15, 0.20))

    # The Sloan Great Wall: outline the structure we measured.
    wall = art["wall"]
    walls = res["great_wall"]["structures"]
    if walls and wall["map"]["structure"] is not None:
        w = walls[0]
        xe, ye = wall["map"]["x_edges"], wall["map"]["y_edges"]
        xc = 0.5 * (xe[:-1] + xe[1:]) / cp.H_LITTLE
        yc = 0.5 * (ye[:-1] + ye[1:]) / cp.H_LITTLE
        from scipy.ndimage import gaussian_filter

        shape = gaussian_filter(wall["map"]["structure"].astype(float), 1.0)
        ax.contour(xc, yc, shape.T, levels=[0.5], colors=ps.INK, linewidths=1.0, zorder=3)
        ra_c = 0.5 * (w["ra_range"][0] + w["ra_range"][1])
        tip = _plane(ra_c - 8, w["distance_range_h"][1] / cp.H_LITTLE, ra_mid)
        ps.annotate(ax, "the Sloan Great Wall (outlined): its densest\nconnected stretch is "
                    f"{w['length_mpc']:,.0f} Mpc ({w['length_mpc'] * cp.MPC_IN_LY / 1e9:.1f} "
                    "billion light-years) end to end",
                    xy=tip, xytext=_plane(ra_c - 30, 0.66 * d_max, ra_mid), fontsize=8.5,
                    ha="center", color=ps.INK, bbox=BACKING, zorder=6)

    v = res["largest_void"]
    vx, vy = _plane(v["ra_center"], v["distance_h"] / cp.H_LITTLE, ra_mid)
    ax.add_patch(Circle((vx, vy), v["radius_mpc"], fill=False, ec=ps.INK, lw=1.0, zorder=3))
    ps.annotate(ax, f"a void: no bright galaxy in\na circle {2 * v['radius_mpc']:.0f} Mpc across",
                xy=(vx + v["radius_mpc"] * 0.75, vy - v["radius_mpc"] * 0.65),
                xytext=_plane(v["ra_center"] - 25, 0.13 * d_max, ra_mid),
                fontsize=8, ha="center", color=ps.INK, bbox=BACKING, zorder=6)

    xt, yt = -0.97 * d_max, 1.12 * d_max
    ps.annotate(ax, "Out here only the most luminous galaxies\nare bright enough to see: the map "
                "thins\nbecause the survey does, not the universe\n(Fig. 02)",
                xy=_plane(233, 0.9 * d_max, ra_mid), xytext=(xt, yt), fontsize=8,
                ha="left", va="top", color=ps.INK_SECONDARY)

    fog = _fog_cluster(art)
    fx, fy = _plane(fog["ra"], fog["distance_mpc"], ra_mid)
    # The zoom window, drawn where it sits: its long side along our line of sight.
    los = np.array([fx, fy]) / np.hypot(fx, fy)
    across = np.array([los[1], -los[0]])  # pointing right on the map, like the zoom's x
    corners = [np.array([fx, fy]) + a * across + b * los
               for a, b in ((-ZOOM_HALF_WIDTH, -ZOOM_HALF_DEPTH), (ZOOM_HALF_WIDTH, -ZOOM_HALF_DEPTH),
                            (ZOOM_HALF_WIDTH, ZOOM_HALF_DEPTH), (-ZOOM_HALF_WIDTH, ZOOM_HALF_DEPTH))]
    ax.add_patch(Polygon(corners, closed=True, fill=False, ec=ps.INK, lw=1.0, zorder=4))

    n_red, n_blue = int(red[sl].sum()), int((~red[sl]).sum())
    handles = [
        Line2D([], [], ls="", marker="o", ms=5, color=ps.C_RED_GALAXY,
               label=f"red-sequence galaxies ({n_red:,})\nold stars, no longer forming new ones"),
        Line2D([], [], ls="", marker="o", ms=5, color=ps.C_BLUE_GALAXY,
               label=f"blue-cloud galaxies ({n_blue:,})\nstill forming stars"),
    ]
    ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(1.0, 0.0), fontsize=8.5,
              labelspacing=0.9)
    ax.set_title(f"The cosmic web: {sl.sum():,} galaxies in a {2 * SLICE_DEC_HALF_WIDTH:g}°-thick "
                 "slice of the universe", fontsize=14, loc="left", pad=4)

    axz = fig.add_subplot(gs[1])
    _finger_panel(axz, art, fog)
    name = fog["abell"] or "a rich cluster"
    axz.set_title(f"Zoom: {name}", fontsize=10, loc="left")
    axz.text(0.03, 0.985, f"a \"finger of God\": {fog['n']} bright\nmembers orbiting at "
             f"~{round(fog['sigma_v'], -1):,.0f} km/s\nlook {fog['elongation']:.0f}x deeper than wide,\n"
             "all pointing at us.\nThe cluster is really\n~round, a few Mpc across.",
             transform=axz.transAxes, fontsize=7.5, va="top", color=ps.INK, bbox=BACKING)
    fig.add_artist(ConnectionPatch(xyA=tuple(corners[1] + 0.5 * (corners[2] - corners[1])),
                                   coordsA=ax.transData,
                                   xyB=(0.0, 0.5), coordsB=axz.transAxes,
                                   color=ps.INK_SECONDARY, lw=0.8))
    ps.caption(fig, "Each dot is a galaxy, placed by its direction on the sky (angle) and the "
               "distance implied by its redshift (radius). Galaxies within 2° of the celestial "
               "equator, right ascension 128°-252°; the slice's thickness grows from ~7 Mpc to "
               "~60 Mpc with distance, so structures are seen in projection. Colour class from "
               "rest-frame g-r (Fig. 06). The zoom uses every catalogued galaxy within ±4 Mpc of "
               "the cluster in declination. " + PROVENANCE, y=0.06)
    return ps.save(fig, FIGURE_DIR / "01_cosmic_web", dpi=300)


# ---------------------------------------------------------------------------
# 2. Why the map thins out
# ---------------------------------------------------------------------------


def figure_selection(res, art):
    """The flux limit, shown two ways: counts vs a uniform universe, and the
    empty corner of the luminosity-redshift plane."""
    g = art["g"]
    z_edges = art["z_edges"]
    zc = 0.5 * (z_edges[:-1] + z_edges[1:])
    obs, exp = art["nz_observed"], art["nz_expected"]

    fig = plt.figure(figsize=(13, 5.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.1, 1], width_ratios=[1, 1.05],
                          hspace=0.08, wspace=0.25)
    ax = fig.add_subplot(gs[0, 0])
    axr = fig.add_subplot(gs[1, 0], sharex=ax)
    ax.stairs(obs, z_edges, color=ps.C_OURS, lw=1.2, label="observed galaxies")
    ax.plot(zc, exp, color=ps.C_MODEL, lw=1.8,
            label="expected if galaxies were spread evenly\n(same flux limit, no structure)")
    ax.set_ylabel("galaxies per Δz = 0.005")
    ax.set_xlim(0, 0.25)
    ax.set_ylim(0, 1.3 * obs.max())
    ax.legend(loc="upper right")
    ax.tick_params(labelbottom=False)
    ax.set_title("A: the drop-off is predicted by the flux limit alone")
    ps.annotate(ax, "fewer and fewer galaxies\nare bright enough\nto be seen",
                xy=(0.17, float(np.interp(0.17, zc, exp))), xytext=(0.185, 0.42 * obs.max()),
                fontsize=8)

    def ratio(o, e, min_expected=50):
        return np.where(e > min_expected, o / np.maximum(e, 1), np.nan)

    r_all = ratio(obs, exp)
    r_sl = ratio(art["nz_slice_observed"], art["nz_slice_expected"], 30)
    axr.axhline(1.0, color=ps.AXIS, lw=0.8)
    axr.stairs(r_sl, z_edges, color=ps.VIOLET, lw=1.2)
    axr.stairs(r_all, z_edges, color=ps.C_OURS, lw=1.2)
    axr.set_ylabel("observed ÷ expected")
    axr.set_xlabel("redshift z")
    axr.set_ylim(0, 2.6)
    k = int(np.nanargmax(np.where((zc > 0.06) & (zc < 0.1), r_sl, np.nan)))
    ps.annotate(axr, "the Sloan Great Wall", xy=(zc[k], r_sl[k]), xytext=(zc[k] + 0.035, 2.25),
                fontsize=8, color=ps.INK)
    axr.text(0.251, 1.45, "equatorial slice", color=ps.VIOLET, fontsize=7.5, ha="right")
    axr.text(0.251, 0.55, "whole survey", color=ps.C_OURS, fontsize=7.5, ha="right")

    axm = fig.add_subplot(gs[:, 1])
    sel = (g["z"] < 0.25) & (g["M_lit"] < -16.5) & (g["M_lit"] > -23.5)
    h = axm.hist2d(g["z"][sel], g["M_lit"][sel], bins=[250, 175],
                   range=[[0, 0.25], [-23.5, -16.5]], cmap=_blue_cmap(), norm=_lognorm(),
                   rasterized=True)
    zz = np.linspace(0.005, 0.25, 200)
    med_col = _running_median(g["z"], g["gr_obs"], zz)
    m_lim = cp.to_literature_magnitude(
        MGS_R_LIMIT - cp.distance_modulus(zz) - cp.k_correction("r", zz, med_col),
        zz, g["band_offset"], EVOLUTION_Q)
    axm.plot(zz, m_lim, color=ps.C_MODEL, lw=1.6)
    ps.annotate(axm, "flux limit r = 17.77:\nfainter than this, invisible",
                xy=(0.2, float(np.interp(0.2, zz, m_lim))), xytext=(0.155, -17.6),
                fontsize=8, color=ps.INK, bbox=BACKING)
    axm.add_patch(plt.Rectangle((VL_Z_MIN, -23.5), VL_Z_MAX - VL_Z_MIN, VL_MAG_LIMIT_H + 23.5,
                                fill=False, ec=ps.INK, lw=1.2))
    axm.text(VL_Z_MIN + 0.003, -23.35, "volume-limited sample:\nevery galaxy brighter than "
             f"{VL_MAG_LIMIT_H:g}\nis visible throughout this box", fontsize=8, va="top",
             color=ps.INK, bbox=BACKING)
    axm.set_xlim(0, 0.25)
    axm.set_ylim(-16.5, -23.5)
    axm.set_xlabel("redshift z")
    axm.set_ylabel("luminosity: absolute magnitude  M_r − 5 log h\n(more negative = brighter)")
    axm.set_title("B: luminous galaxies are seen far away, faint ones only nearby")
    axm.grid(False)
    _colorbar(fig, h[3], axm, "galaxies per cell")

    ps.caption(fig, "A: redshift distribution of all 656,795 galaxies, against the expectation "
               "for a universe without structure, built from the galaxies' own luminosities by "
               "the 1/Vmax method (Schmidt 1968). Their ratio isolates structure: walls push it "
               "up, voids pull it down; in the thin equatorial slice a single wall dominates. "
               "B: absolute magnitudes K-corrected and evolved to z = 0.1 as in Zehavi et al. "
               "(2011); the flux-limit curve is for the median galaxy colour. " + PROVENANCE,
               y=0.0)
    return ps.save(fig, FIGURE_DIR / "02_selection")


def _running_median(x, y, at, width=0.01):
    out = np.empty(len(at))
    for i, a in enumerate(at):
        m = np.abs(x - a) < width / 2
        out[i] = np.median(y[m]) if m.sum() > 20 else np.nan
    good = np.isfinite(out)
    return np.interp(at, at[good], out[good])


# ---------------------------------------------------------------------------
# 3. Flux-limited vs volume-limited
# ---------------------------------------------------------------------------


def figure_volume_limited(res, art):
    """The same slice, twice: everything the survey saw, and a fair census --
    then the density of each against distance, over the whole survey."""
    g, sl, vl, mask = art["g"], art["slice"], art["vl"], art["mask"]
    ra_mid = art["ra_mid"]
    d_max = float(cp.comoving_distance(0.12))
    in_view = sl & (g["z"] < 0.12)

    fig = plt.figure(figsize=(16, 5.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.3, 1.3, 0.9], wspace=0.22)
    panels = (
        (in_view, f"A: flux-limited — all {in_view.sum():,} galaxies the survey saw"),
        (in_view & vl, f"B: volume-limited — the {(in_view & vl).sum():,} brighter than "
                       f"M = {VL_MAG_LIMIT_H:g}"),
    )
    for k, (m, title) in enumerate(panels):
        ax = fig.add_subplot(gs[k])
        x, y = _plane(g["ra"][m], g["d_mpc"][m], ra_mid)
        ax.scatter(x, y, s=0.5 if k == 0 else 0.9, color=ps.C_GALAXY, linewidths=0,
                   alpha=0.8, rasterized=True)
        wedge_frame(ax, ra_mid, d_max, (0.05, 0.10), ring_labels=k == 0, scalebar=k == 0,
                    label_size=7)
        ax.set_title(title, fontsize=10, loc="left")
        if k == 1:
            d_edge = float(cp.comoving_distance(VL_Z_MAX))
            _arc(ax, d_edge, ra_mid, color=ps.ORANGE, lw=1.2, zorder=3)
            xt, yt = _plane(ra_mid, d_edge + 0.05 * d_max, ra_mid)
            ax.text(xt, yt, f"sample edge, z = {VL_Z_MAX:.3f}", color=ps.ORANGE, fontsize=7,
                    ha="center", va="bottom", bbox=BACKING)

    # Density against distance over the whole survey (interior footprint).
    axd = fig.add_subplot(gs[2])
    inside = mask.contains(g["ra"], g["dec"], "interior")
    omega = mask.solid_angle_sr()
    d_lo, d_hi = float(cp.comoving_distance(VL_Z_MIN)), float(cp.comoving_distance(VL_Z_MAX))
    curves = {}
    for key, m, edges, col, mfc, lab in (
            ("flux", inside, np.arange(40.0, 521.0, 20.0), ps.C_GALAXY, ps.SURFACE,
             "flux-limited: every galaxy"),
            ("volume", inside & vl, np.linspace(d_lo, d_hi, 19), ps.INK, ps.INK,
             "volume-limited")):
        n, _ = np.histogram(g["d_mpc"][m], bins=edges)
        dens = n / cp.shell_volume(omega, edges[:-1], edges[1:]) * 1e6
        centres = 0.5 * (edges[:-1] + edges[1:])
        axd.plot(centres, dens, "o-", color=col, mfc=mfc, ms=4, lw=1.3, label=lab)
        curves[key] = (centres, dens)
    axd.set_yscale("log")
    axd.set_xlabel("distance (Mpc)")
    axd.set_ylabel("galaxies per million Mpc³")
    axd.axvline(d_hi, color=ps.ORANGE, lw=1.0)
    axd.legend(loc="upper right", fontsize=7.5)
    axd.set_title("C: density against distance, whole survey", fontsize=10, loc="left")
    axd.set_box_aspect(0.95)

    fc, fd = curves["flux"]
    near = float(np.mean(fd[(fc > 90) & (fc < 150)]))
    far = float(np.mean(fd[(fc > d_hi - 40) & (fc < d_hi + 20)]))
    vc, vd = curves["volume"]
    spread = float(np.std(np.log(vd)))
    ps.caption(fig, "A and B: the equatorial slice of Fig. 01, out to z = 0.12. The "
               "volume-limited sample keeps only galaxies luminous enough to be seen anywhere "
               f"out to z = {VL_Z_MAX:.3f}, so its density means density, and the walls and voids "
               "are all still there. C: over the whole survey, the flux-limited density falls "
               f"{near / far:.0f}-fold between ~120 and ~{d_hi - 10:.0f} Mpc, while the "
               f"volume-limited density stays level (point-to-point scatter ±{100 * spread:.0f}%: "
               "real structures). Every measurement in this project uses the volume-limited "
               "sample. " + PROVENANCE, y=0.04)
    return ps.save(fig, FIGURE_DIR / "03_flux_vs_volume_limited", dpi=250)


# ---------------------------------------------------------------------------
# 4. Is it random?
# ---------------------------------------------------------------------------


def figure_random(res, art):
    """Real galaxies next to the same number of points placed at random."""
    g, sl, vl = art["g"], art["slice"], art["vl"]
    ra_mid = art["ra_mid"]
    m = sl & vl
    rnd = art["slice_random"]
    d_max = float(cp.comoving_distance(0.115))

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.9), gridspec_kw={"wspace": 0.05})
    for ax, (ra, d, title) in zip(axes, (
            (g["ra"][m], g["d_mpc"][m], f"A: the real universe — {m.sum():,} galaxies"),
            (rnd["ra"], rnd["d_h"] / cp.H_LITTLE,
             "B: the same number of points, placed at random"))):
        x, y = _plane(ra, d, ra_mid)
        ax.scatter(x, y, s=1.0, color=ps.C_GALAXY, linewidths=0, alpha=0.85, rasterized=True)
        wedge_frame(ax, ra_mid, d_max, (0.05, 0.10), ring_labels=ax is axes[0],
                    scalebar=ax is axes[0], label_size=7)
        ax.set_title(title, fontsize=10.5, loc="left")
    cl = res["clustering"]
    f = cl["fits"]["diagonal"]
    ps.caption(fig, "Both panels share the volume-limited sample's selection: same sky, same "
               f"distance range (z = {VL_Z_MIN}-{VL_Z_MAX:.3f}), same number of points, spread "
               "uniformly through the volume on the right. Only the arrangement differs. Measured "
               f"(Fig. 05): a galaxy is ~{1 + _xi_at(f, 1.0 * cp.H_LITTLE):.0f}x more likely than "
               "chance to have a neighbour 1 Mpc away, and still twice as likely at "
               f"r0 = {f['r0'] / cp.H_LITTLE:.1f} Mpc. " + PROVENANCE, y=0.05)
    return ps.save(fig, FIGURE_DIR / "04_real_vs_random", dpi=250)


def _xi_at(fit, r_h):
    """Real-space xi(r) of a power-law fit to w_p, at r in h^-1 Mpc."""
    return (r_h / fit["r0"]) ** (-fit["gamma"])


# ---------------------------------------------------------------------------
# 5. Clustering, measured
# ---------------------------------------------------------------------------


def figure_clustering(res, art):
    """w_p(r_p) against Zehavi et al. (2011), and the 2D xi(r_p, pi) map."""
    clus = art["clustering"]
    r = clus["result"]
    fits = res["clustering"]["fits"]
    fd, fc = fits["diagonal"], fits["full_cov"]
    z = lit.ZEHAVI_MR20

    fig = plt.figure(figsize=(13.5, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.1], width_ratios=[1.15, 1],
                          hspace=0.07, wspace=0.2)
    ax = fig.add_subplot(gs[0, 0])
    axr = fig.add_subplot(gs[1, 0], sharex=ax)

    lo, hi = fd["rp_range"]
    rr = np.logspace(np.log10(lo), np.log10(hi), 100)
    ax.plot(rr, cp.wp_power_law(rr, fd["r0"], fd["gamma"]), color=ps.C_MODEL, lw=1.4, zorder=2,
            label=f"power-law fit: r0 = {fd['r0']:.2f} ± {fd['r0_err']:.2f} h⁻¹ Mpc, "
                  f"γ = {fd['gamma']:.2f}")
    ax.errorbar(lit.ZEHAVI_WP_RP, lit.ZEHAVI_WP_MR20, yerr=lit.ZEHAVI_WP_MR20_ERR, fmt="s",
                ms=6.5, mfc="none", mec=ps.C_PUBLISHED, ecolor=ps.C_PUBLISHED, elinewidth=1,
                capsize=0, zorder=3, label="Zehavi et al. (2011), SDSS DR7, same sample cuts")
    ax.errorbar(r.rp, r.wp, yerr=r.wp_err, fmt="o", ms=4.5, color=ps.C_OURS,
                ecolor=ps.C_OURS, elinewidth=1, capsize=0, zorder=4,
                label=f"this project: {r.n_data:,} galaxies, {r.n_random:,} random points")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("projected correlation  w_p(r_p)  (h⁻¹ Mpc)")
    ax.tick_params(labelbottom=False)
    ax.legend(loc="lower left", fontsize=7.5)
    ax.set_title("A: how strongly galaxies cluster, by separation")
    ax.axvspan(0.1, 0.2, color=ps.GRID, lw=0, zorder=0)
    ax.set_xlim(0.1, 55)
    ax.text(0.97, 0.97,
            "Power-law fits, r < 20 h⁻¹ Mpc (ours; Zehavi et al.)\n"
            f"independent errors:  r0 = {fd['r0']:.2f} ± {fd['r0_err']:.2f};  "
            f"{z['r0_diag']:.2f} ± {z['r0_diag_err']:.2f}\n"
            f"correlated errors:    r0 = {fc['r0']:.2f} ± {fc['r0_err']:.2f};  "
            f"{z['r0']:.2f} ± {z['r0_err']:.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=ps.INK,
            linespacing=1.5)

    ratio = r.wp / np.interp(np.log(r.rp), np.log(lit.ZEHAVI_WP_RP), lit.ZEHAVI_WP_MR20)
    axr.axhline(1.0, color=ps.AXIS, lw=0.8)
    zrel = lit.ZEHAVI_WP_MR20_ERR / lit.ZEHAVI_WP_MR20
    axr.fill_between(lit.ZEHAVI_WP_RP, 1 - zrel, 1 + zrel, color=ps.C_PUBLISHED, alpha=0.12,
                     lw=0)
    axr.errorbar(r.rp, ratio,
                 yerr=r.wp_err / np.interp(np.log(r.rp), np.log(lit.ZEHAVI_WP_RP),
                                           lit.ZEHAVI_WP_MR20),
                 fmt="o", ms=4, color=ps.C_OURS, elinewidth=1, capsize=0)
    axr.set_ylim(0.6, 1.3)
    axr.set_xlabel("separation across the line of sight  r_p  (h⁻¹ Mpc)")
    axr.set_ylabel("ours ÷ theirs")
    axr.text(0.105, 0.64, "fibre\ncollisions", fontsize=7, color=ps.MUTED, va="bottom")

    # The 2D map, mirrored into four quadrants as is conventional.
    from matplotlib.colors import BoundaryNorm
    from scipy.ndimage import gaussian_filter

    axm = fig.add_subplot(gs[:, 1])
    xi = gaussian_filter(np.nan_to_num(r.xi_map), 0.8)
    full = np.block([[xi[::-1, ::-1], xi[::-1, :]], [xi[:, ::-1], xi]]).T
    ext = clus["binning"].map_max
    grid = np.linspace(-ext, ext, full.shape[0] + 1)
    centres = 0.5 * (grid[:-1] + grid[1:])
    levels = [0.25, 0.5, 1, 2, 4, 8, 16, 32]
    cmap = _blue_cmap()
    cf = axm.contourf(centres, centres, full, levels=levels, cmap=cmap, extend="max",
                      norm=BoundaryNorm(levels, cmap.N, extend="max"))
    axm.contour(centres, centres, full, levels=levels, colors=ps.INK, linewidths=0.6)
    th = np.linspace(0, 2 * np.pi, 300)
    for lev in levels:
        rad = fd["r0"] * lev ** (-1.0 / fd["gamma"])
        if rad < ext:
            axm.plot(rad * np.cos(th), rad * np.sin(th), color=ps.ORANGE, lw=0.9)
    axm.set_aspect("equal")
    axm.set_xlim(-ext, ext)
    axm.set_ylim(-ext, ext)
    axm.set_xlabel("across the line of sight  r_p  (h⁻¹ Mpc)")
    axm.set_ylabel("along the line of sight  π  (h⁻¹ Mpc)")
    axm.set_title("B: the same pairs, split by direction")
    axm.grid(False)
    ps.annotate(axm, "fingers of God:\nstretched along the line\nof sight by the orbital\n"
                "speeds of cluster galaxies", xy=(0.6, 13.0), xytext=(6.5, 21.0), fontsize=7.5,
                color=ps.INK, bbox=BACKING)
    ps.annotate(axm, "large scales: squashed,\nbecause galaxies are falling\ntogether "
                "(Kaiser effect)", xy=(-17.5, -6.5), xytext=(-29.0, -26.5), fontsize=7.5,
                color=ps.INK, bbox=BACKING)
    handles = [Line2D([], [], color=ps.INK, lw=0.8, label="measured ξ contours"),
               Line2D([], [], color=ps.ORANGE, lw=0.9,
                      label="if galaxies did not move: circles")]
    axm.legend(handles=handles, loc="upper left", fontsize=7, framealpha=1,
               facecolor=ps.SURFACE)
    cb = _colorbar(fig, cf, axm, "ξ: excess of pairs over random", ticks=levels,
                   spacing="uniform")
    cb.ax.set_yticklabels([f"{lv:g}" for lv in levels])

    ps.caption(fig, "Landy-Szalay estimator on the volume-limited sample (Mr − 5 log h < −20, "
               "0.02 < z < 0.106), random catalogue 3x larger with the same angular mask and "
               "uniform in volume; errors from 50 jackknife regions. w_p integrates ξ(r_p, π) to "
               "π = 60 h⁻¹ Mpc, which cancels the velocity distortions seen in B. The line in A is "
               "the fit with independent errors; the correlated-errors fit is statistically "
               "preferred but, as for Zehavi et al., sits below the points (see README). Shaded: "
               "fibre-collision scales, uncorrected here. B lightly smoothed (0.8 h⁻¹ Mpc); circles "
               "use the fit in A. " + PROVENANCE, y=0.0)
    return ps.save(fig, FIGURE_DIR / "05_clustering")


# ---------------------------------------------------------------------------
# 6. Colours
# ---------------------------------------------------------------------------


def figure_colors(res, art):
    """The colour-magnitude diagram: two populations, and the line between."""
    g = art["g"]
    d = art["divider"]
    sel = (g["z"] < 0.2) & (g["gr0"] > 0.1) & (g["gr0"] < 1.1) & (g["M"] < -17) & (g["M"] > -24)

    fig, (ax, axh) = plt.subplots(1, 2, figsize=(13, 5.2),
                                  gridspec_kw={"width_ratios": [1.5, 1], "wspace": 0.25})
    h = ax.hist2d(g["M"][sel], g["gr0"][sel], bins=[200, 160], range=[[-24, -17], [0.1, 1.1]],
                  cmap=_blue_cmap(), norm=_lognorm(), rasterized=True)
    mm = np.linspace(-23.5, -17.5, 50)
    ax.plot(mm, d["slope"] * mm + d["intercept"], color=ps.INK, lw=1.4)
    for b in d["bins"]:
        ax.plot(b["M_center"], b["divider"], "o", ms=4.5, color=ps.INK,
                mfc=ps.INK if b["bimodal"] else ps.SURFACE)
    ax.text(-22.4, 0.95, "red sequence", color=ps.C_RED_GALAXY, fontsize=10.5,
            fontweight="bold", ha="center")
    ax.text(-18.3, 0.2, "blue cloud", color=ps.C_BLUE_GALAXY, fontsize=10.5,
            fontweight="bold", ha="center", bbox=BACKING)
    sign = "−" if d["intercept"] < 0 else "+"
    ps.annotate(ax, f"divider:\ng − r = {d['slope']:.3f} × M_r {sign} {abs(d['intercept']):.2f}",
                xy=(-22.0, d["slope"] * -22.0 + d["intercept"]), xytext=(-22.6, 0.3),
                fontsize=8, color=ps.INK, bbox=BACKING, ha="center")
    ax.set_xlim(-17, -24)
    ax.set_xlabel("absolute magnitude M_r  (more negative = more luminous)")
    ax.set_ylabel("rest-frame colour  g − r  (larger = redder)")
    ax.set_title(f"A: {sel.sum():,} galaxies, z < 0.2")
    ax.grid(False)
    _colorbar(fig, h[3], ax, "galaxies per cell")

    used = [b for b in d["bins"] if b["bimodal"]]
    b = max(used, key=lambda b: b["n"])
    m = sel & (np.abs(g["M"] - b["M_center"]) < 0.25)
    binw = 0.01
    axh.hist(g["gr0"][m], bins=np.arange(0.1, 1.1 + binw / 2, binw), color=ps.MUTED,
             histtype="stepfilled", alpha=0.45, lw=0)
    xx = np.linspace(0.1, 1.1, 400)
    for mu, sd, w, col in ((b["blue_mean"], b["blue_sd"], 1 - b["red_weight"], ps.C_BLUE_GALAXY),
                           (b["red_mean"], b["red_sd"], b["red_weight"], ps.C_RED_GALAXY)):
        axh.plot(xx, m.sum() * binw * w * np.exp(-0.5 * ((xx - mu) / sd) ** 2)
                 / (sd * np.sqrt(2 * np.pi)), color=col, lw=1.6)
    axh.axvline(b["divider"], color=ps.INK, lw=1.2)
    axh.text(b["divider"] + 0.012, 0.97, "equal odds:\nred or blue", fontsize=8, va="top",
             transform=axh.get_xaxis_transform())
    axh.set_xlabel("rest-frame colour  g − r")
    axh.set_ylabel("galaxies per 0.01 mag")
    axh.set_title(f"B: one column of A, M_r = {b['M_center'] + 0.25:.2f} to "
                  f"{b['M_center'] - 0.25:.2f}")
    ps.caption(fig, "Colours are model-magnitude g − r, corrected for Milky Way dust and "
               "K-corrected to the rest frame (Chilingarian et al. 2010, 2012). In each 0.5-mag "
               "column a two-Gaussian model is fitted and the colour of equal odds found (points); "
               "the divider is a line through those columns where the two populations are "
               "separated by a real dip (filled points). Brighter than M_r ≈ −21 the blue "
               "population is only a shoulder on the red one, so the line is extrapolated there. "
               + PROVENANCE, y=0.0)
    return ps.save(fig, FIGURE_DIR / "06_colors")


# ---------------------------------------------------------------------------
# 7. Environment
# ---------------------------------------------------------------------------


def figure_environment(res, art):
    """Red fraction against local density, and where each type lives."""
    env = res["environment"]
    g, sl, vl, red = art["g"], art["slice"], art["vl"], art["red"]
    ra_mid = art["ra_mid"]
    fig = plt.figure(figsize=(15.5, 5.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.85, 1.25, 1.25], wspace=0.06)
    ax = fig.add_subplot(gs[0])
    for key, col, mfc, lab in (("bright", ps.INK, ps.INK, "luminous (M < −21)"),
                               ("faint", ps.MUTED, ps.SURFACE, "fainter (−21 < M < −20)")):
        e = env[key]
        ax.errorbar(e["sigma_median"], e["red_fraction"], yerr=e["red_fraction_err"], fmt="o-",
                    color=col, mfc=mfc, ms=5, lw=1.3, elinewidth=1, capsize=0, label=lab)
    ax.set_xscale("log")
    ax.set_xlabel("local density Σ₅  (galaxies per h⁻² Mpc²)")
    ax.set_ylabel("fraction of galaxies that are red")
    ax.set_ylim(0, 1)
    ax.legend(loc="lower right", fontsize=7.5, title="magnitude M_r − 5 log h",
              title_fontsize=7.5)
    ax.set_title("A: crowded places are red places")
    ax.text(0.03, 0.97, "sparsest ← → densest\n(each point: 1/8 of the galaxies)",
            transform=ax.transAxes, fontsize=7.5, va="top", color=ps.MUTED)

    m = sl & vl
    d_max = float(cp.comoving_distance(0.115))
    for k, (mask, col, title) in enumerate((
            (m & red, ps.C_RED_GALAXY, f"B: red-sequence galaxies ({(m & red).sum():,})"),
            (m & ~red, ps.C_BLUE_GALAXY, f"C: blue-cloud galaxies ({(m & ~red).sum():,})"))):
        axw = fig.add_subplot(gs[k + 1])
        x, y = _plane(g["ra"][mask], g["d_mpc"][mask], ra_mid)
        axw.scatter(x, y, s=1.2, color=col, linewidths=0, alpha=0.85, rasterized=True)
        wedge_frame(axw, ra_mid, d_max, (0.05, 0.10), ring_labels=k == 0, scalebar=k == 0,
                    label_size=7)
        axw.set_title(title, fontsize=10, loc="left")
    ps.caption(fig, "A: Σ₅ = 5 / (π d₅²), with d₅ the projected distance to the 5th-nearest "
               f"neighbour within ±{DENSITY_DV_KMS:.0f} km/s, for volume-limited galaxies away "
               "from the survey's edges; equal-count bins, binomial errors (smaller than the "
               "points). The rise holds at fixed luminosity, so it is not merely that crowded "
               "regions host brighter galaxies. B, C: the slice's volume-limited galaxies by "
               "colour class. Both trace the same web — it is the web of all matter — but the red "
               "ones are concentrated into its densest knots, the blue ones spread more evenly "
               "along filaments and into sparse regions. " + PROVENANCE, y=0.03)
    return ps.save(fig, FIGURE_DIR / "07_environment", dpi=250)


# ---------------------------------------------------------------------------
# Everything
# ---------------------------------------------------------------------------


def make_all(res, art) -> list:
    written = []
    for fn in (figure_footprint, figure_wedge, figure_selection, figure_volume_limited,
               figure_random, figure_clustering, figure_colors, figure_environment):
        written += fn(res, art)
        print(f"  {fn.__name__}")
    import interactive

    written.append(interactive.cosmic_web_3d(res, art))
    print("  interactive 3D map")
    return written
