"""Figure generation.

Every figure is drawn from the analysis artefacts, never by re-running the
analysis, so what is plotted is exactly what was measured.

Design notes. Most panels are an *emphasis* form: cluster members in the accent
blue, the field as recessive grey context, because the story is always "this
group, against everything else". Densities use one-hue sequential ramps. The
three-cluster figures use the first three categorical slots, the most an
all-pairs scatter can carry and stay colour-vision safe (checked with a
validator, not by eye). Text is always in ink tones, never in a data colour.
Each figure carries a caption with its provenance and the caveat a reader
needs in order not to over-read it.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Ellipse, Rectangle

import gaia_physics as gp
import plotstyle as ps
from config import FIGURE_DIR, MEMBERSHIP_THRESHOLD, SUN_ABS_G, SUN_BP_RP, WIDTH_BP_RP_RANGE

ps.use_style()

PROVENANCE = "Data: ESA Gaia DR3, read from the CDS VizieR mirror (I/355/gaiadr3)."

#: A very light neutral ramp: star-count shading that must stay behind the dots.
CMAP_SHADE = LinearSegmentedColormap.from_list("shade", ["#f6f5f2", "#dedcd5", "#c3c2b7"])

#: How names read in running text ("the Pleiades", but "M67").
TITLE_NAME = {"pleiades": "the Pleiades"}

BOX = {"boxstyle": "round,pad=0.3", "facecolor": ps.SURFACE, "edgecolor": "none", "alpha": 0.85}


def _title_name(cluster):
    return TITLE_NAME.get(cluster.key, cluster.name)


def _size_by_mag(g, g_faint, g_bright, s_max=70.0, s_min=0.5):
    """Marker area from magnitude, so bright stars look bright -- like a star chart."""
    g = np.asarray(g, dtype=float)
    frac = np.clip((g_faint - g) / max(g_faint - g_bright, 1e-6), 0, 1)
    return s_min + (s_max - s_min) * frac ** 2.2


def _sky_axes(ax, radius, labels=True):
    """A flat map of the sky, east to the left as it appears from the ground."""
    ax.add_patch(Circle((0, 0), radius, fill=False, edgecolor=ps.AXIS, linewidth=0.8))
    ax.set_xlim(radius * 1.04, -radius * 1.04)
    ax.set_ylim(-radius * 1.04, radius * 1.04)
    ax.set_aspect("equal")
    ax.grid(False)
    if labels:
        ax.set_xlabel("Degrees east of centre (east is left, as on the sky)")
        ax.set_ylabel("Degrees north of centre")


def _place_labels(ax, points, fontsize=7.0, color=None, avoid=()):
    """Label points without overlaps: try positions around each point, keep the first clear one.

    ``points`` is a list of (x, y, text) in data coordinates; ``avoid`` holds
    display-space boxes that must stay clear (legends, notes).
    """
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    placed = list(avoid)
    # Each point's own marker is an obstacle too.
    for x, y, _ in points:
        px, py = ax.transData.transform((x, y))
        placed.append(_box(px - 4, py - 4, px + 4, py + 4))
    options = [((7, 0), "left", "center"), ((-7, 0), "right", "center"),
               ((0, 7), "center", "bottom"), ((0, -7), "center", "top"),
               ((7, 7), "left", "bottom"), ((7, -7), "left", "top"),
               ((-7, 7), "right", "bottom"), ((-7, -7), "right", "top"),
               ((16, 12), "left", "bottom"), ((-16, -12), "right", "top"),
               ((16, -12), "left", "top"), ((-16, 12), "right", "bottom")]
    for x, y, text in points:
        chosen = None
        for (dx, dy), ha, va in options:
            ann = ax.annotate(text, xy=(x, y), xytext=(dx, dy), textcoords="offset points",
                              ha=ha, va=va, fontsize=fontsize, color=color or ps.INK_SECONDARY,
                              zorder=7, arrowprops={"arrowstyle": "-", "color": ps.MUTED,
                                                    "linewidth": 0.5, "shrinkA": 0,
                                                    "shrinkB": 3} if abs(dx) > 10 else None)
            bb = ann.get_window_extent(renderer)
            if not any(bb.overlaps(p) for p in placed):
                chosen = ann
                placed.append(bb)
                break
            ann.remove()
        if chosen is None:  # nothing clear: accept the first option rather than drop it
            ann = ax.annotate(text, xy=(x, y), xytext=(7, 0), textcoords="offset points",
                              ha="left", va="center", fontsize=fontsize,
                              color=color or ps.INK_SECONDARY, zorder=7)
            placed.append(ann.get_window_extent(renderer))


def _box(x0, y0, x1, y1):
    from matplotlib.transforms import Bbox
    return Bbox([[x0, y0], [x1, y1]])


def _named(art, result, cluster):
    """Named stars with sky and HR-diagram positions, and whether we selected them.

    Sky positions are relative to the cone centre, like every other sky panel.
    Stars our test does not select are placed in the HR diagram at the
    cluster's distance, for reference, and drawn hollow.
    """
    named = art.get("named")
    if named is None or not len(named):
        return []
    df = art["df"].set_index("source_id")
    dm = result.distance["distance_modulus"]
    out = []
    for _, row in named.iterrows():
        sid = int(row["source_id"])
        xi, eta = gp.gnomonic(row["ra"], row["dec"], cluster.ra_deg, cluster.dec_deg)
        member = sid in df.index and bool(df.loc[sid, "member"])
        if member:
            mg, status = float(df.loc[sid, "abs_g"]), "member"
        else:
            mg = float(row["phot_g_mean_mag"]) - dm
            status = "no Gaia motion" if not np.isfinite(row["parallax"]) else "not selected"
        out.append({"name": str(row["name"]), "xi": float(xi), "eta": float(eta),
                    "g": float(row["phot_g_mean_mag"]), "bp_rp": float(row["bp_rp"]),
                    "abs_g": mg, "status": status})
    return out


def _hrd_axes(ax, ylim=(15.5, -3.5), xlim=(-0.7, 4.0)):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("Colour  BP − RP  (mag; blue, hot stars left)")
    ax.set_ylabel("Absolute magnitude  $M_G$  (bright up)")


def _hrd_limits(df, members):
    """Axis limits that fit this cluster's members, with room for labels."""
    mg = df.loc[members, "abs_g"].to_numpy()
    mg = mg[np.isfinite(mg)]
    top = max(float(np.min(mg)) - 1.2, -3.8)
    bottom = min(float(np.percentile(mg, 99.8)) + 1.6, 15.5)
    return (bottom, top), (-0.7, 4.0)


def _field_hexbin(ax, df, mask, xlim, ylim, gridsize=140):
    hb = ax.hexbin(df.loc[mask, "bp_rp"], df.loc[mask, "abs_g"], gridsize=gridsize,
                   extent=(xlim[0], xlim[1], min(ylim), max(ylim)), bins="log", mincnt=1,
                   cmap=ps.CMAP_CONTEXT, linewidths=0)
    hb.set_rasterized(True)
    return hb


def _white_dwarfs(df, members):
    s = members & df["bp_rp"].between(-0.8, 0.8).to_numpy() & (df["abs_g"] > 9.5).to_numpy()
    return df.loc[s, ["bp_rp", "abs_g"]].to_numpy()


def _corner(ax, x, y):
    """Axes-fraction anchor in the corner nearest a data point."""
    fx, fy = ax.transAxes.inverted().transform(ax.transData.transform((x, y)))
    right, top = fx > 0.5, fy > 0.5
    return ((0.97 if right else 0.03), (0.97 if top else 0.03),
            "right" if right else "left", "top" if top else "bottom")


# --------------------------------------------------------------------------
# 1. The story
# --------------------------------------------------------------------------

def figure_story(result, art, cluster):
    """Four panels: the field, the motions, the members on the sky, their HR diagram.

    The whole project in one figure. Panel 1 is roughly what a telescope
    shows: a crowded field where only the brightest members stand out. Panel 2
    is what Gaia's motions show: one tight knot. Panels 3 and 4 follow that
    knot back onto the sky and into the HR diagram -- neither of which was
    used to pick the members -- where it turns out to be concentrated and to
    form one thin sequence.
    """
    df, knot = art["df"], art["knot"]
    members = df["member"].to_numpy()
    r = cluster.radius_deg
    xi, eta = df["xi"].to_numpy(), df["eta"].to_numpy()
    g = df["phot_g_mean_mag"].to_numpy()
    fig = plt.figure(figsize=(12.5, 12.0))
    gs = fig.add_gridspec(2, 2, hspace=0.28, wspace=0.26, top=0.91, bottom=0.07)

    # --- 1. the field ---------------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    hb = ax.hexbin(xi, eta, gridsize=48, extent=(-r, r, -r, r), mincnt=1, cmap=CMAP_SHADE,
                   linewidths=0)
    hb.set_rasterized(True)
    n_dots = min(len(df), 1500)
    bright = np.argsort(g)[:n_dots]
    g_cut = float(g[bright].max())
    order = bright[np.argsort(-g[bright])]
    ax.scatter(xi[order], eta[order], s=_size_by_mag(g[order], g_cut + 0.3, g[order].min()),
               color=ps.INK, linewidths=0, zorder=3, rasterized=True)
    _sky_axes(ax, r)
    ax.set_title(f"1 · The field: {len(df):,} stars within {r:g}°")
    ax.text(0.02, 0.02, f"shading: all {len(df):,} stars to G = {cluster.g_max:g}\n"
            f"dots: the {n_dots:,} brightest (G < {g_cut:.1f})", transform=ax.transAxes,
            fontsize=8, color=ps.INK_SECONDARY, va="bottom")

    # --- 2. the motions -------------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    fx, fy = knot["field_median_pm"]
    kx, ky = result.motion["pmra_members_mean"], result.motion["pmdec_members_mean"]
    cx, cy = 0.5 * (fx + kx), 0.5 * (fy + ky)
    half = 0.5 * max(abs(kx - fx), abs(ky - fy)) + 14
    bins = np.arange(-half, half + 0.4, 0.4)
    h, ex, ey = np.histogram2d(df["pmra"], df["pmdec"], bins=[bins + cx, bins + cy])
    h = np.ma.masked_less(h, 1)
    mesh = ax.pcolormesh(ex, ey, h.T, cmap=ps.CMAP_CONTEXT, norm=LogNorm(vmin=1, vmax=h.max()),
                         shading="flat", rasterized=True)
    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy - half, cy + half)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.add_patch(Circle((kx, ky), 3.2, fill=False, edgecolor=ps.C_MEMBER, linewidth=1.6))
    for text, (px, py), colour in (
            (f"{cluster.name}: {result.membership['n']:,} stars\nmoving together", (kx, ky), ps.INK),
            ("the field: mostly distant stars,\nwhose motions look small\n"
             "because they are far away", (fx, fy), ps.INK_SECONDARY)):
        tx, ty, ha, va = _corner(ax, px, py)
        ax.annotate(text, xy=(px, py), xytext=(tx, ty), textcoords="axes fraction", ha=ha,
                    va=va, fontsize=8, color=colour, bbox=BOX,
                    arrowprops={"arrowstyle": "-", "color": colour, "linewidth": 0.8,
                                "shrinkA": 2, "shrinkB": 18 if colour == ps.INK else 4})
    ax.set_xlabel("Proper motion east, $\\mu_{\\alpha*}$ (mas/yr)")
    ax.set_ylabel("Proper motion north, $\\mu_\\delta$ (mas/yr)")
    ax.set_title("2 · Their motions: one knot stands out")
    cb = fig.colorbar(mesh, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("stars per 0.4 × 0.4 mas/yr cell", color=ps.INK_SECONDARY)
    cb.outline.set_visible(False)

    # --- 3. members on the sky -------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    rng = np.random.default_rng(1)
    field = np.where(~members)[0]
    show = rng.choice(field, size=min(len(field), 40000), replace=False)
    ax.scatter(xi[show], eta[show], s=0.3, color=ps.GRID, linewidths=0, rasterized=True)
    mo = np.where(members)[0]
    mo = mo[np.argsort(-g[mo])]
    msize = np.maximum(_size_by_mag(g[mo], cluster.g_max, g[mo].min(), s_max=90), 4)
    ax.scatter(xi[mo], eta[mo], s=msize, color=ps.C_MEMBER, edgecolors=ps.SURFACE,
               linewidths=0.5, zorder=3, rasterized=True)
    _sky_axes(ax, r)
    ax.set_title(f"3 · The knot's stars on the sky: {result.membership['n']:,} members")
    ax.text(0.02, 0.02, "chosen by motion and parallax only —\nposition was never used, "
            "yet they crowd together", transform=ax.transAxes, fontsize=8,
            color=ps.INK_SECONDARY, va="bottom", bbox=BOX, zorder=6)
    stars = _named(art, result, cluster)
    if stars:
        _sisters_inset(ax, stars, df, members, g, cluster)

    # --- 4. HR diagram --------------------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    ylim, xlim = _hrd_limits(df, members)
    _field_hexbin(ax, df, art["field_hrd"], xlim, ylim, gridsize=120)
    good = art["clean_members"]
    ax.scatter(df.loc[good, "bp_rp"], df.loc[good, "abs_g"], s=7, color=ps.C_MEMBER,
               edgecolors=ps.SURFACE, linewidths=0.3, zorder=3, rasterized=True)
    _hrd_axes(ax, ylim=ylim, xlim=xlim)
    hrd = result.hr_diagram
    ax.set_title("4 · Their HR diagram: one thin sequence")
    ax.text(0.97, 0.97,
            f"main-sequence width\n{cluster.name}: {hrd['width_members_mag']:.3f} mag\n"
            f"field stars: {hrd['width_field_mag']:.2f} mag\n"
            f"→ {hrd['width_ratio']:.1f}× thinner", transform=ax.transAxes, ha="right",
            va="top", fontsize=8.5, color=ps.INK, linespacing=1.5)
    wd = _white_dwarfs(df, members)
    if len(wd):
        label = ("EGGR 25, the Pleiades'\nonly white dwarf" if cluster.key == "pleiades"
                 else f"{len(wd)} white dwarfs")
        wx, wy = np.median(wd, axis=0)
        ps.annotate(ax, label, xy=(wx, wy), xytext=(xlim[0] + 0.08, wy - 2.6), color=ps.INK)
    handles = [Line2D([], [], marker="o", ls="", color=ps.C_MEMBER, markeredgecolor=ps.SURFACE,
                      label=f"{cluster.name} members"),
               Line2D([], [], marker="h", ls="", color="#a3a198", markersize=8,
                      label="field stars in the same patch of sky")]
    ax.legend(handles=handles, loc="lower left")

    fig.suptitle(f"Finding {_title_name(cluster)} by how its stars move", x=0.08, ha="left",
                 y=0.975)
    ps.caption(fig, f"{PROVENANCE} Members come from a mixture model of proper motion and "
               f"parallax; colour, brightness and position on the sky were never used to "
               f"select them, so the concentration in panel 3 and the thin sequence in panel 4 "
               f"are results, not assumptions. Absolute magnitudes use each member's "
               f"parallax, refined by the cluster's own distance; field stars use their own "
               f"parallax (better than 10%). No extinction correction (A_V ≈ 0.1 here).")
    return ps.save(fig, FIGURE_DIR / f"{cluster.key}_01_story")


def _sisters_inset(ax, stars, df, members, g, cluster):
    """A zoom on the named stars, which crowd into the centre of the full map."""
    xs = np.array([s["xi"] for s in stars])
    ys = np.array([s["eta"] for s in stars])
    pad = 0.28
    x0, x1 = xs.min() - pad, xs.max() + pad
    y0, y1 = ys.min() - pad, ys.max() + pad
    side = max(x1 - x0, y1 - y0)
    xc, yc = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
    x0, x1, y0, y1 = xc - side / 2, xc + side / 2, yc - side / 2, yc + side / 2
    ins = ax.inset_axes([0.585, 0.585, 0.41, 0.41])
    xi, eta = df["xi"].to_numpy(), df["eta"].to_numpy()
    inside = members & (xi > x0) & (xi < x1) & (eta > y0) & (eta < y1)
    idx = np.where(inside)[0]
    idx = idx[np.argsort(-g[idx])]
    ins.scatter(xi[idx], eta[idx], s=np.maximum(_size_by_mag(g[idx], cluster.g_max, 2.5, s_max=110), 5),
                color=ps.C_MEMBER, edgecolors=ps.SURFACE, linewidths=0.5, zorder=3,
                rasterized=True)
    for s in stars:
        if s["status"] != "member":
            ins.scatter([s["xi"]], [s["eta"]], s=55, facecolors="none", edgecolors=ps.INK,
                        linewidths=1.0, zorder=4)
    ins.set_xlim(x1, x0)
    ins.set_ylim(y0, y1)
    ins.set_aspect("equal")
    ins.set_xticks([])
    ins.set_yticks([])
    ins.grid(False)
    for spine in ins.spines.values():
        spine.set_visible(True)
        spine.set_color(ps.AXIS)
    ins.set_facecolor(ps.SURFACE)
    ins.text(0.03, 0.97, "zoom: the named stars\nhollow = not confirmed", transform=ins.transAxes,
             ha="left", va="top", fontsize=7, color=ps.INK_SECONDARY, bbox=BOX, zorder=8)
    ax.add_patch(Rectangle((x0, y0), side, side, fill=False, edgecolor=ps.INK_SECONDARY,
                           linewidth=0.8, zorder=6))
    _place_labels(ins, [(s["xi"], s["eta"], s["name"]) for s in stars], fontsize=7)




# --------------------------------------------------------------------------
# 2. Membership evidence
# --------------------------------------------------------------------------

def figure_membership(result, art, cluster):
    """How the members were chosen, and the evidence that the choice is right."""
    df, mix3 = art["df"], art["mix3"]
    members = df["member"].to_numpy()
    kin = df["kinematic_candidate"].to_numpy()
    rejected = kin & ~members
    region = df["fit_region"].to_numpy()
    plx0 = result.distance["parallax_mas"]
    fig, axs = plt.subplots(2, 2, figsize=(12, 10.2))
    fig.subplots_adjust(hspace=0.36, wspace=0.26, top=0.9, bottom=0.08)

    # --- motion space, in the cluster's frame ---------------------------
    ax = axs[0, 0]
    core = mix3.core
    mx, my = mix3.means[core][:2]
    sd = np.sqrt(np.max(np.diag(mix3.covs[mix3.cluster[-1]][:2, :2])))
    w = float(np.clip(4.5 * sd, 1.5, 8.0))
    near = (region & (np.abs(df["pmra_c"] - mx) < w).to_numpy()
            & (np.abs(df["pmdec_c"] - my) < w).to_numpy())
    fld = near & ~members & ~rejected
    ax.scatter(df.loc[fld, "pmra_c"], df.loc[fld, "pmdec_c"], s=4, color=ps.C_FIELD, alpha=0.6,
               linewidths=0, rasterized=True, label="field")
    ax.scatter(df.loc[near & members, "pmra_c"], df.loc[near & members, "pmdec_c"], s=5,
               color=ps.C_MEMBER, linewidths=0, rasterized=True, label="members")
    ax.scatter(df.loc[near & rejected, "pmra_c"], df.loc[near & rejected, "pmdec_c"], s=16,
               color=ps.C_REJECTED, edgecolors=ps.SURFACE, linewidths=0.5, zorder=4,
               label="same motion, wrong distance")
    for j in mix3.cluster:
        vals, vecs = np.linalg.eigh(mix3.covs[j][:2, :2])
        angle = np.degrees(np.arctan2(vecs[1, 1], vecs[0, 1]))
        for k in (1, 2):
            ax.add_patch(Ellipse(mix3.means[j][:2], 2 * k * np.sqrt(vals[1]), 2 * k * np.sqrt(vals[0]),
                                 angle=angle, fill=False, linewidth=1.1 if j == core else 0.7,
                                 edgecolor=ps.INK if j == core else ps.INK_SECONDARY, zorder=5))
    ax.set_xlim(mx - w, mx + w)
    ax.set_ylim(my - w, my + w)
    ax.set_aspect("equal")
    ax.set_xlabel("$\\mu_{\\alpha*}$ in the cluster's frame (mas/yr)")
    ax.set_ylabel("$\\mu_\\delta$ in the cluster's frame (mas/yr)")
    ax.set_title("The model in motion space")
    ax.legend(loc="lower left", markerscale=1.5, framealpha=0.85, frameon=True,
              facecolor=ps.SURFACE, edgecolor="none")
    halo = "; thin: its halo" if len(mix3.cluster) > 1 else ""
    ax.text(0.98, 0.98, f"ellipses: 1σ and 2σ of the cluster's\nown spread (bold: core{halo}),\n"
            "with measurement blur removed", transform=ax.transAxes, ha="right", va="top",
            fontsize=7.5, color=ps.INK_SECONDARY, bbox=BOX, zorder=6)

    # --- parallax: the independent check ---------------------------------
    ax = axs[0, 1]
    lo, hi = 0.0, max(2.2 * plx0, plx0 + 3)
    bins = np.linspace(lo, hi, 120)
    ax.hist(df.loc[kin & members, "parallax_corr"], bins=bins, color=ps.C_MEMBER,
            label=f"confirmed ({result.kinematic_only['n_confirmed_by_parallax']:,})")
    ax.hist(df.loc[rejected, "parallax_corr"], bins=bins, color=ps.C_REJECTED,
            label=f"rejected ({result.kinematic_only['n_rejected_by_parallax']:,})")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_xlabel("Parallax (mas; larger = nearer)")
    ax.set_ylabel("Stars per bin (log scale)")
    ax.set_title(f"Motion picked {result.kinematic_only['n']:,}; parallax checks them")
    ax.legend(loc="upper right")
    ax.annotate(f"{plx0:.3f} mas = {result.distance['distance_pc']:.1f} pc",
                xy=(plx0, 1.0), xycoords=("data", "axes fraction"), xytext=(-8, -6),
                textcoords="offset points", ha="right", va="top", fontsize=8,
                color=ps.INK_SECONDARY, bbox=BOX)
    ax.axvline(plx0, color=ps.INK_SECONDARY, linewidth=0.8, zorder=0)

    # --- how decisive the model is ---------------------------------------
    ax = axs[1, 0]
    p = df.loc[region, "p_member"].to_numpy()
    ax.hist(p[p > 0.01], bins=np.linspace(0.01, 1.0, 50), color=ps.C_MEMBER)
    ax.axvline(MEMBERSHIP_THRESHOLD, color=ps.INK_SECONDARY, linewidth=0.8)
    ax.set_yscale("log")
    ax.set_xlabel("Membership probability")
    ax.set_ylabel("Stars per bin (log scale)")
    m = result.membership
    ax.set_title("Most stars are clearly in or clearly out")
    ax.text(0.03, 0.97, f"{m['n_ambiguous']} stars between 10% and 90%\n"
            f"expected field stars among the members: {m['expected_contamination']:.0f}\n"
            f"expected members missed: {m['expected_missed']:.0f}",
            transform=ax.transAxes, va="top", fontsize=8, color=ps.INK_SECONDARY)
    ax.annotate("threshold: 50%", xy=(MEMBERSHIP_THRESHOLD, 0.55),
                xycoords=("data", "axes fraction"), xytext=(4, 0), textcoords="offset points",
                fontsize=8, color=ps.INK_SECONDARY)

    # --- concentration on the sky ------------------------------------------
    ax = axs[1, 1]
    prof = art["profile"]
    ok = prof["density"] > 0
    ax.errorbar(prof["r_mid"][ok], prof["density"][ok], yerr=prof["density_err"][ok], fmt="o",
                color=ps.C_MEMBER, markersize=5, markeredgecolor=ps.SURFACE, elinewidth=1.0,
                capsize=0, label="members")
    ax.axhline(prof["contamination_density"], color=ps.C_REJECTED, linewidth=1.6,
               label="expected field contamination")
    ax.set_yscale("log")
    ax.set_xlim(0, cluster.radius_deg)
    ax.set_xlabel("Distance from the cluster centre on the sky (deg)")
    ax.set_ylabel("Members per square degree")
    ax.set_title("Members crowd towards the centre")
    r50 = result.structure["r50_deg"]
    ax.axvline(r50, color=ps.INK_SECONDARY, linewidth=0.8)
    ax.annotate(f"half the members within\n{r50:.2f}° = {result.structure['r50_pc']:.1f} pc",
                xy=(r50, 0.5), xycoords=("data", "axes fraction"), xytext=(6, 0),
                textcoords="offset points", fontsize=8, color=ps.INK_SECONDARY, va="center")
    ax.legend(loc="upper right")

    fig.suptitle(f"{cluster.name} — who belongs, and how we know", x=0.08, ha="left", y=0.975)
    ps.caption(fig, f"{PROVENANCE} Top left: proper motions corrected for perspective (see "
               f"the 3D figure). Top right: stars chosen by motion alone, before parallax is "
               f"used; the ones parallax rejects are field stars that happen to share the "
               f"motion at a different distance. Bottom: probabilities from the final model "
               f"(motion + parallax); the contamination level is the model's own estimate, "
               f"spread evenly over the field.")
    return ps.save(fig, FIGURE_DIR / f"{cluster.key}_02_membership")


# --------------------------------------------------------------------------
# 3. The HR diagram
# --------------------------------------------------------------------------

def _box_median(df, sel, lo_c, hi_c, lo_m, hi_m, min_stars=3):
    s = (sel & df["bp_rp"].between(lo_c, hi_c).to_numpy()
         & df["abs_g"].between(lo_m, hi_m).to_numpy())
    if s.sum() < min_stars:
        return None
    return float(df.loc[s, "bp_rp"].median()), float(df.loc[s, "abs_g"].median()), int(s.sum())


def _annotate_features(ax, df, members, art, cluster):
    """Label only features that are really there, located from the data."""
    feats = []
    wd = _box_median(df, members, -0.8, 0.8, 9.5, 16, min_stars=1)
    if wd:
        text = ("EGGR 25 — the Pleiades' only\nknown white dwarf" if cluster.key == "pleiades"
                else f"{wd[2]} white dwarfs")
        # Above the white dwarfs, where the diagram is empty; below is the legend.
        feats.append((text, wd[:2], (ax.get_xlim()[0] + 0.1, wd[1] - 2.6)))
    giants = _box_median(df, members, 0.95, 2.5, -2.5, 2.3)
    if giants:
        feats.append((f"red giants ({giants[2]})", giants[:2], (giants[0] + 0.7, giants[1] - 0.8)))
    if cluster.key == "m67":
        for text, box, off in (("blue stragglers", (-0.3, 0.62, -1, 2.8), (-0.25, -1.2)),
                               ("red clump", (1.05, 1.35, 0.0, 1.2), (0.7, 0.8)),
                               ("subgiants", (0.8, 1.02, 2.3, 3.4), (0.75, 0.35)),
                               ("turnoff", (0.6, 0.8, 2.6, 3.6), (-0.8, 0.9))):
            f = _box_median(df, members, *box)
            if f:
                feats.append((text, f[:2], (f[0] + off[0], f[1] + off[1])))
    rc, rm = art["ridge"]
    k = np.argmin(np.abs(rc - 1.7))
    room_right = ax.get_xlim()[1] - rc[k] > 1.6
    feats.append(("main sequence", (rc[k], rm[k]),
                  (rc[k] + 0.8, rm[k] - 1.3) if room_right else (rc[k] - 1.45, rm[k] - 0.9)))
    for text, xy, xytext in feats:
        ps.annotate(ax, text, xy=xy, xytext=xytext, color=ps.INK, ha="left")


def _display_ridge(df, good):
    """The members' sequence over its whole length, for drawing over the field."""
    c = df.loc[good, "bp_rp"].to_numpy()
    m = df.loc[good, "abs_g"].to_numpy()
    lo = float(np.percentile(c, 1))
    hi = float(np.percentile(c, 99))
    import analysis
    return analysis.ridge_line(c, m, lo, hi, step=0.1, min_stars=4)


def figure_hrd(result, art, cluster):
    """The members' HR diagram against the field's, and the width of each."""
    df = art["df"]
    members = df["member"].to_numpy()
    good = art["clean_members"]
    rc, rm = art["ridge"]
    hrd = result.hr_diagram
    fig = plt.figure(figsize=(15.5, 6.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.15, 0.8], wspace=0.24, bottom=0.12)
    ylim, xlim = _hrd_limits(df, members)
    if cluster.key == "m67":
        xlim = (-0.3, 2.8)

    # --- members --------------------------------------------------------
    ax = fig.add_subplot(gs[0])
    flagged = members & ~good
    ax.scatter(df.loc[flagged, "bp_rp"], df.loc[flagged, "abs_g"], s=11, facecolors="none",
               edgecolors=ps.C_MEMBER, linewidths=0.6, label="members, photometry flagged")
    ax.scatter(df.loc[good, "bp_rp"], df.loc[good, "abs_g"], s=8, color=ps.C_MEMBER,
               edgecolors=ps.SURFACE, linewidths=0.3, label="members", zorder=3)
    grid = np.linspace(rc.min(), rc.max(), 200)
    ax.plot(grid, np.interp(grid, rc, rm) - 0.753, color=ps.INK_SECONDARY, linewidth=0.9,
            zorder=2, label="where equal-mass binaries sit")
    stars = _named(art, result, cluster)
    for s in stars:
        if s["status"] != "member":
            ax.scatter([s["bp_rp"]], [s["abs_g"]], s=34, marker="D", facecolors="none",
                       edgecolors=ps.INK, linewidths=0.9, zorder=4)
    _hrd_axes(ax, ylim=ylim, xlim=xlim)
    _annotate_features(ax, df, members, art, cluster)
    ax.set_title(f"{cluster.name}: {result.membership['n']:,} members")
    leg = ax.legend(loc="lower left")
    bright = [s for s in stars if s["g"] < 4.5]
    if bright:
        fig.canvas.draw()
        _place_labels(ax, [(s["bp_rp"], s["abs_g"], s["name"] + ("" if s["status"] == "member"
                                                              else f" ({s['status']})"))
                           for s in sorted(bright, key=lambda s: s["abs_g"])],
                      avoid=[leg.get_window_extent()])

    # --- the field ---------------------------------------------------------
    ax = fig.add_subplot(gs[1])
    _field_hexbin(ax, df, art["field_hrd"], xlim, ylim, gridsize=110)
    dc, dm = _display_ridge(df, good)
    ax.plot(dc, dm, color=ps.C_MEMBER, linewidth=2.0, label=f"{cluster.name}'s sequence (median)")
    _hrd_axes(ax, ylim=ylim, xlim=xlim)
    ax.set_title(f"Field stars in the same patch of sky ({art['field_hrd'].sum():,})")
    ax.legend(loc="lower left")
    ax.text(0.97, 0.97, "all ages, compositions and\ndistances mixed together",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color=ps.INK_SECONDARY)

    # --- widths -------------------------------------------------------------
    ax = fig.add_subplot(gs[2])
    bins = np.linspace(-1.5, 1.5, 61)
    wf, wm = art["width_field"]["residuals"], art["width_members"]["residuals"]
    ax.hist(wf - np.median(wf), bins=bins, density=True, color="#c3c2b7",
            label=f"field: {hrd['width_field_mag']:.2f} mag")
    ax.hist(wm - np.median(wm), bins=bins, density=True, histtype="step", linewidth=2.0,
            color=ps.C_MEMBER, label=f"{cluster.name}: {hrd['width_members_mag']:.3f} mag")
    ax.set_xlabel("Distance from the sequence's spine (mag)\n← fainter · brighter →")
    ax.set_ylabel("Fraction of stars (normalised)")
    ax.set_title(f"{hrd['width_ratio']:.1f}× thinner than the field")
    ax.invert_xaxis()
    ax.legend(loc="upper left", title="robust width (1σ)", title_fontsize=8)
    ax.text(0.03, 0.62, f"dwarfs with {WIDTH_BP_RP_RANGE[0]} < BP−RP < {WIDTH_BP_RP_RANGE[1]};\n"
            f"field stars with parallax\nbetter than 5%",
            transform=ax.transAxes, fontsize=7.5, color=ps.INK_SECONDARY, va="top")

    ps.caption(fig, f"{PROVENANCE} A cluster's stars share an age, a composition and a "
               f"distance, so they trace one line; the field mixes all three. The cluster's "
               f"residual width comes from unresolved binaries, stellar spots and rotation, "
               f"and patchy dust. Hollow circles: members with RUWE > 1.4 or inconsistent "
               f"BP/RP photometry, left out of the width. Diamonds: named stars our test "
               f"cannot confirm, placed at the cluster's distance. No extinction correction.")
    return ps.save(fig, FIGURE_DIR / f"{cluster.key}_03_hr_diagram")


# --------------------------------------------------------------------------
# 4. The cluster in three dimensions
# --------------------------------------------------------------------------

def figure_3d(result, art, cluster):
    """Perspective across the sky, and nearer members moving faster.

    Only drawn where the cluster's depth can be seen in the parallaxes (the
    Pleiades and Praesepe, not M67).
    """
    if not result.structure["depth_resolvable"]:
        return []
    df = art["df"]
    m, per = result.motion, result.perspective
    fig, axs = plt.subplots(1, 2, figsize=(13, 6.4))
    fig.subplots_adjust(wspace=0.28, top=0.84, bottom=0.12)

    # --- perspective: measured vs predicted, in cells on the sky ------------
    ax = axs[0]
    r = cluster.radius_deg
    cells = art["perspective_cells"]
    scale = per["step_deg"] / 1.3  # degrees of arrow per mas/yr
    for c in cells:
        x, y = c["xi"], c["eta"]
        for vec, colour, lw in ((c["predicted"], ps.C_MODEL, 2.4), (c["observed"], ps.C_MEMBER, 1.3)):
            # xi already increases eastwards, and the inverted axis puts east on
            # the left, so the east component is added as it is.
            ax.annotate("", xy=(x + vec[0] * scale, y + vec[1] * scale), xytext=(x, y),
                        arrowprops={"arrowstyle": "-|>", "color": colour, "linewidth": lw,
                                    "shrinkA": 0, "shrinkB": 0, "mutation_scale": 9})
        ax.plot([x], [y], "o", color=ps.INK, markersize=2.5, zorder=5)
    _sky_axes(ax, r)
    ax.plot([-r * 0.95, -r * 0.95 + scale], [-r * 0.93, -r * 0.93], color=ps.INK_SECONDARY,
            linewidth=1.4, solid_capstyle="butt")
    ax.text(-r * 0.95 + 0.5 * scale, -r * 0.88, "1 mas/yr", fontsize=8, ha="center",
            color=ps.INK_SECONDARY)
    ax.legend(handles=[Line2D([], [], color=ps.C_MEMBER, lw=1.6, label="measured (median of members)"),
                       Line2D([], [], color=ps.C_MODEL, lw=2.4, label="predicted from one shared velocity")],
              loc="upper right", frameon=True, facecolor=ps.SURFACE, edgecolor="none",
              framealpha=0.85)
    ax.set_title("Perspective: each part of the cluster's motion,\nminus the centre's")
    shrink, fan = gp.perspective_terms(cluster.dec_deg, m["pmra"], m["pmdec"],
                                       result.distance["parallax_mas"], m["rv_kms"])
    if shrink > fan:
        note = (f"{'receding' if m['rv_kms'] > 0 else 'approaching'} fast, so it appears to "
                f"{'shrink' if m['rv_kms'] > 0 else 'grow'}:\nmotions tilt "
                f"{'towards' if m['rv_kms'] > 0 else 'away from'} the centre")
    else:
        cp = m["convergent_point_deg"]
        sep = float(gp.angular_separation(cluster.ra_deg, cluster.dec_deg, *cp))
        pa = gp.position_angle(cluster.ra_deg, cluster.dec_deg, *cp)
        compass = ["north", "north-east", "east", "south-east", "south", "south-west", "west",
                   "north-west"][int(((pa + 22.5) % 360) // 45)]
        note = (f"parallel motions in perspective: all aim at one\nvanishing point, "
                f"{sep:.0f}° away to the {compass}")
    ax.text(0.02, 0.98, note, transform=ax.transAxes, fontsize=8, color=ps.INK_SECONDARY,
            va="top")
    ax.text(0.02, 0.02,
            f"pattern: {per['rms_offset_masyr']:.2f} mas/yr rms\n"
            f"after subtracting the prediction: {per['rms_after_model_masyr']:.2f}\n"
            f"measurement noise alone: {per['rms_cell_noise_masyr']:.2f}",
            transform=ax.transAxes, fontsize=8, color=ps.INK, va="bottom", linespacing=1.5)

    # --- nearer members move faster ------------------------------------------
    ax = axs[1]
    ds = m["depth_slope"]["along_motion"]
    sel = (df["member"] & df["ruwe_ok"] & (df["parallax_corr"] / df["parallax_error"] > 50)).to_numpy()
    mu = np.array([m["pmra"], m["pmdec"]])
    along = df.loc[sel, ["pmra_c", "pmdec_c"]].to_numpy(float) @ (mu / np.hypot(*mu))
    plx = df.loc[sel, "parallax_corr"].to_numpy(float)
    ax.scatter(plx, along, s=7, color=ps.C_MEMBER, alpha=0.7, linewidths=0, rasterized=True,
               label=f"members with parallax better than 2% ({sel.sum():,})")
    plx0 = result.distance["parallax_mas"]
    xs = np.linspace(np.percentile(plx, 0.5), np.percentile(plx, 99.5), 50)
    base = float(np.hypot(*mu))
    ax.plot(xs, base + ds["predicted"] * (xs - plx0), color=ps.C_MODEL, linewidth=2.6,
            label=f"prediction: slope {ds['predicted']:.2f}")
    ax.plot(xs, base + ds["slope"] * (xs - plx0), color=ps.INK, linewidth=1.1,
            label=f"measured: slope {ds['slope']:.2f} ± {ds['err']:.2f}")
    ax.set_xlabel("Parallax (mas; right = nearer)")
    ax.set_ylabel("Proper motion along the direction of travel (mas/yr)")
    ax.set_title("Nearer members move faster across the sky")
    ax.legend(loc="upper left")
    lo, hi = np.percentile(along, [0.5, 99.5])
    pad = 0.3 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)
    s = result.structure
    ax.text(0.98, 0.03, f"depth along our line of sight: σ = {s['depth_sigma_pc']:.1f} pc\n"
            f"width across the sky: σ = {s['width_sigma_pc']:.1f} pc\n→ roughly round",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color=ps.INK,
            linespacing=1.5)

    direction = "receding" if m["rv_kms"] > 0 else "approaching"
    fig.suptitle(f"{cluster.name} in three dimensions: {direction} at {abs(m['rv_kms']):.1f} km/s, "
                 f"crossing our view at {m['tangential_speed_kms']:.1f} km/s",
                 x=0.08, ha="left", y=0.99)
    ps.caption(fig, f"{PROVENANCE} Left: stars sharing one velocity show slightly different "
               f"proper motions across a wide field, because we see them from different "
               f"angles. Arrows show each cell's difference from the centre's motion: measured "
               f"(median of its members, RUWE < 1.4) and predicted from the 3D velocity (mean "
               f"motion, parallax, and Gaia radial velocities). Right: proper motion is velocity "
               f"divided by distance, so one shared velocity makes nearer members move "
               f"proportionally faster; the measured slope is corrected for parallax scatter.")
    return ps.save(fig, FIGURE_DIR / f"{cluster.key}_04_three_d")


# --------------------------------------------------------------------------
# 00. Across clusters
# --------------------------------------------------------------------------

def figure_age_ladder(results, arts, clusters, literature):
    """Three clusters of different ages on one HR diagram."""
    fig, ax = plt.subplots(figsize=(9.5, 9.8))
    fig.subplots_adjust(top=0.93, bottom=0.08)
    label_at = {"pleiades": (-0.55, -1.2), "praesepe": (0.62, -1.0), "m67": (1.28, 2.55)}
    wds = []
    for colour, r, art, c in zip(ps.CLUSTER_COLOURS, results, arts, clusters):
        df = art["df"]
        good = art["clean_members"] | (df["member"] & (df["phot_g_mean_mag"] <= 6)).to_numpy()
        av = (literature["clusters"][c.key]["hr23"] or {}).get("AV50") or 0.0
        a_g, e_bprp = gp.extinction_from_av(av)
        col, mag = df.loc[good, "bp_rp"] - e_bprp, df.loc[good, "abs_g"] - a_g
        ax.scatter(col, mag, s=9, color=colour, edgecolors=ps.SURFACE, linewidths=0.3,
                   zorder=3, label=f"{c.name} ({good.sum():,} stars)", rasterized=True)
        wd = (col < 0.8) & (mag > 9.5)
        wds += list(zip(col[wd], mag[wd]))
        age = c.reference_age
        age_txt = (f"{age.myr / 1000:.1f} billion years" if age.myr >= 1000
                   else f"{age.myr:.0f} million years")
        x, y = label_at.get(c.key, (0, 0))
        ax.text(x, y, f"{c.name}\n~{age_txt}", fontsize=10, color=ps.INK,
                fontweight="bold", va="center")
        ax.plot([x - 0.06], [y], "o", color=colour, markersize=9,
                markeredgecolor=ps.SURFACE, markeredgewidth=1.5)
    ax.plot([SUN_BP_RP], [SUN_ABS_G], "o", markersize=10, markerfacecolor="none",
            markeredgecolor=ps.INK, markeredgewidth=1.2, zorder=6)
    ax.plot([SUN_BP_RP], [SUN_ABS_G], "o", markersize=2.2, color=ps.INK, zorder=6)
    ps.annotate(ax, "the Sun (4.6 billion years)\nsits on M67's sequence",
                xy=(SUN_BP_RP + 0.06, SUN_ABS_G + 0.08), xytext=(1.35, 5.9), color=ps.INK,
                ha="left")
    ps.annotate(ax, "turnoff: the top of the main sequence\nburns away as a cluster ages",
                xy=(0.72, 2.9), xytext=(-0.55, 6.5), color=ps.INK_SECONDARY, ha="left")
    if wds:
        wx, wy = np.median(np.array(wds), axis=0)
        ps.annotate(ax, "white dwarfs: stars that\nhave already died", xy=(wx, wy),
                    xytext=(wx + 0.5, wy - 1.8), color=ps.INK_SECONDARY, ha="left")
    _hrd_axes(ax, ylim=(13.0, -3.2), xlim=(-0.7, 3.8))
    ax.set_xlabel("Colour  (BP − RP)$_0$  (mag, corrected for dust)")
    ax.set_ylabel("Absolute magnitude  $M_{G,0}$  (bright up)")
    ax.set_title("Three clusters, three ages: stars leave the main sequence from the top down")
    ax.legend(loc="lower left")
    ages = "; ".join(f"{c.name}: {c.reference_age.citation} ({c.reference_age.method})"
                     for c in clusters)
    ps.caption(fig, f"{PROVENANCE} Members with clean photometry, found by this project's "
               f"pipeline; dust corrections use Hunt & Reffert (2023) A_V with mean Gaia "
               f"coefficients. M67 is six times further away, so its faintest stars are "
               f"beyond Gaia's limit here. Ages are published values, not measured here: {ages}.")
    return ps.save(fig, FIGURE_DIR / "00_age_ladder")


def figure_comparison(results):
    """Our measurements against the literature, in units of the combined uncertainty."""
    rows = [("pmra", "Proper motion, RA", "HR23"), ("pmdec", "Proper motion, Dec", "HR23"),
            ("parallax_uncorrected", "Parallax", "HR23"), ("distance_pc", "Distance", "HR23"),
            ("rv_kms", "Radial velocity", "HR23"), ("rv_kms_gc2018", "Radial velocity", "GC2018")]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.2),
                                  gridspec_kw={"width_ratios": [1.5, 1.0], "wspace": 0.55})
    fig.subplots_adjust(top=0.82, bottom=0.14)
    ypos = np.arange(len(rows))
    ax.axvspan(-2, 2, color="#f0efec", zorder=0)
    ax.axvspan(-1, 1, color="#e1e0d9", zorder=0)
    ax.axvline(0, color=ps.INK_SECONDARY, linewidth=0.9, zorder=1)
    lim = 4.0
    for k, (colour, r) in enumerate(zip(ps.CLUSTER_COLOURS, results)):
        for i, (key, _, _) in enumerate(rows):
            c = r.comparison.get(key)
            if not c or c.get("ours_err") is None or c.get("published_err") is None:
                continue
            pull = c["difference"] / np.hypot(c["ours_err"], c["published_err"])
            y = i + (k - 1) * 0.22
            off_scale = abs(pull) > lim - 0.2
            shown = float(np.clip(pull, -lim + 0.2, lim - 0.2))
            ax.plot([shown], [y], "<" if off_scale and pull < 0 else (">" if off_scale else "o"),
                    color=colour, markersize=8, markeredgecolor=ps.SURFACE,
                    markeredgewidth=1.5, zorder=4, label=r.cluster if i == 0 else None)
            if off_scale:
                ax.annotate(f"{pull:+.1f}σ (off scale)", xy=(shown, y),
                            xytext=(9 if pull < 0 else -9, 0), textcoords="offset points",
                            ha="left" if pull < 0 else "right", va="center", fontsize=8,
                            color=ps.INK)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{label}  ({src})" for _, label, src in rows])
    ax.invert_yaxis()
    ax.set_xlim(-lim, lim)
    ax.set_xlabel("Our value minus published, in combined standard errors (σ)\n"
                  "shaded: within 1σ and 2σ")
    ax.set_title("Our measurements vs the literature")
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.07), ncol=3)

    labels = ["their members we recover", "our members they list"]
    for k, (colour, r) in enumerate(zip(ps.CLUSTER_COLOURS, results)):
        mv = r.membership_vs_published
        for i, v in enumerate((100 * mv["recovered_fraction"], 100 * mv["ours_confirmed_fraction"])):
            y = i + (k - 1) * 0.22
            ax2.plot([v], [y], "o", color=colour, markersize=8, markeredgecolor=ps.SURFACE,
                     markeredgewidth=1.5, zorder=4)
            ax2.annotate(f"{v:.0f}%", xy=(v, y), xytext=(8, 0), textcoords="offset points",
                         va="center", fontsize=8, color=ps.INK_SECONDARY)
    ax2.set_yticks([0, 1])
    ax2.set_yticklabels(labels)
    ax2.set_ylim(1.5, -0.5)
    ax2.set_xlim(85, 102)
    ax2.set_xticks([85, 90, 95, 100])
    ax2.set_xlabel("Percent")
    ax2.set_title("Star by star, vs Hunt & Reffert (2023)")
    ax2.grid(axis="y", visible=False)

    ps.caption(fig, "HR23 = Hunt & Reffert (2023, Gaia DR3); GC2018 = Gaia Collaboration, "
               "Babusiaux et al. (2018, Gaia DR2). Motions compare plain averages of members' "
               "measured proper motions, as HR23 compute theirs; the parallax row compares "
               "means before zero-point correction, since HR23 correct only their distances. "
               "GC2018's radial velocities come from Gaia DR2, whose calibration differs from "
               "DR3's. Membership is compared within our cone and magnitude limit.")
    return ps.save(fig, FIGURE_DIR / "00_comparison")


def make_all(results, arts, clusters, literature):
    written = []
    for r, art, c in zip(results, arts, clusters):
        written += figure_story(r, art, c)
        written += figure_membership(r, art, c)
        written += figure_hrd(r, art, c)
        written += figure_3d(r, art, c)
    if len(results) == 3:
        written += figure_age_ladder(results, arts, clusters, literature)
        written += figure_comparison(results)
    return written
