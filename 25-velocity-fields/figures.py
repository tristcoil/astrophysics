"""Every static figure of the project (PNG + SVG in ``figures/``).

Colour roles (``plotstyle.py``): ink is our measurement, violet a published
THINGS result, aqua SPARC's curves from older telescopes, orange models we
fit or simulate; blue/red are the approaching/receding halves of a disc.
Velocity maps use a diverging blue-grey-red scale centred on the systemic
velocity, intensity maps a single-hue grey ramp.
"""

from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Ellipse  # noqa: E402

import plotstyle as P  # noqa: E402
from config import FIGURE_DIR, GALAXIES, DERIVED_DIR, galaxy_by_name  # noqa: E402

P.use_style()

ORDER = [g.key for g in GALAXIES]
NAMES = {g.key: g.name for g in GALAXIES}
KIND_LABELS = {"gh_vpeak": "Gauss-Hermite, peak", "gh_v0": "Gauss-Hermite, centre",
               "gauss_v": "single Gaussian", "peak_v": "brightest channel", "mom1_things": "first moment"}


def _maps(key):
    with np.load(DERIVED_DIR / f"{key}_maps.npz") as f:
        return {k: f[k] for k in f.files}


def _crop(mask, margin=8):
    rows, cols = np.nonzero(mask)
    return (max(rows.min() - margin, 0), rows.max() + margin + 1, max(cols.min() - margin, 0), cols.max() + margin + 1)


def _extent(gr, frame, r0, r1, c0, c1, unit=60.0):
    """imshow extent in arcmin east/north of the dynamical centre (east to the left)."""
    x0, y0 = gr.primary.centre
    xs = ((np.array([c0, c1]) - 0.5 - frame.cref) * frame.pixel - x0) / unit
    ys = ((np.array([r0, r1]) - 0.5 - frame.rref) * frame.pixel - y0) / unit
    # x (west) to the right means east to the left: label as offsets east, decreasing
    return [-xs[0], -xs[1], ys[0], ys[1]]


def _frame_for(gr):
    import analysis as A
    import things
    return A.frame_of(things.load_cube(galaxy_by_name(gr.name)))


def _beam(ax, gr, x, y, unit=60.0):
    bmaj, bmin, bpa = gr.beam
    ax.add_patch(Ellipse((x, y), bmin / unit, bmaj / unit, angle=bpa, facecolor=P.INK_SECONDARY,
                         edgecolor="none", alpha=0.8))


def _rings(ax, gr, every=8, unit=60.0, color=P.INK, lw=0.6):
    tr = gr.primary
    t = np.linspace(0, 2 * np.pi, 200)
    for k in range(every - 1, len(tr.radius), every):
        pa, inc, r = np.radians(tr.pa[k]), np.radians(tr.incl[k]), tr.radius[k]
        along, across = r * np.cos(t), r * np.sin(t) * np.cos(inc)
        x = -along * np.sin(pa) - across * np.cos(pa)  # invert disc_coordinates' rotation
        y = along * np.cos(pa) - across * np.sin(pa)
        ax.plot(-x / unit, y / unit, color=color, lw=lw, alpha=0.55)


def _vel_levels(vsys, vmax, step):
    n = int(np.ceil(vmax / step))
    return vsys + step * np.arange(-n, n + 1)


def _colourbar(fig, im, ax, label):
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label(label, color=P.INK_SECONDARY)
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=P.MUTED, labelsize=7)
    return cb


def _map_axes(ax, title):
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("arcmin east of centre")
    ax.set_ylabel("arcmin north")
    ax.grid(False)
    ax.set_aspect("equal")


# --------------------------------------------------------------------------
# 01 The spider diagram
# --------------------------------------------------------------------------

def fig_spider(grs, key="NGC3198"):
    gr = grs[key]
    m = _maps(key)
    fr = _frame_for(gr)
    vf, model, resid = m["velocity"], m["model"], m["residual"]
    r0, r1, c0, c1 = _crop(np.isfinite(vf), 12)
    ext = _extent(gr, fr, r0, r1, c0, c1)
    vsys = gr.primary.vsys
    span = float(np.nanpercentile(np.abs(vf - vsys), 99))
    cmap = P.velocity_cmap()

    fig, axes = plt.subplots(1, 4, figsize=(17, 5.3), constrained_layout=True)
    mom0 = m["mom0"][r0:r1, c0:c1]
    from massmodel import hi_surface_density
    sig = hi_surface_density(np.where(mom0 > 0, mom0, np.nan), gr.beam)
    im = axes[0].imshow(sig, origin="lower", extent=ext, cmap=P.intensity_cmap(),
                        vmin=0, vmax=np.nanpercentile(sig, 99.5))
    _colourbar(fig, im, axes[0], "HI surface density seen on the sky (Msun/pc², no helium)")
    _map_axes(axes[0], "a  The hydrogen")
    for ax, data, title in ((axes[1], vf, "b  Measured velocity field"),
                            (axes[2], model, "c  Tilted-ring model")):
        im = ax.imshow(data[r0:r1, c0:c1], origin="lower", extent=ext, cmap=cmap,
                       vmin=vsys - span, vmax=vsys + span)
        levels = _vel_levels(vsys, span, 20)
        xx = np.linspace(ext[0], ext[1], c1 - c0)
        yy = np.linspace(ext[2], ext[3], r1 - r0)
        ax.contour(xx, yy, np.ma.masked_invalid(data[r0:r1, c0:c1]), levels=levels, colors=P.INK,
                   linewidths=0.5, alpha=0.8)
        ax.contour(xx, yy, np.ma.masked_invalid(data[r0:r1, c0:c1]), levels=[vsys], colors=P.INK,
                   linewidths=1.3)
        _map_axes(ax, title)
    _colourbar(fig, im, axes[2], "line-of-sight velocity (km/s)")
    _rings(axes[0], gr, color=P.C_MODEL, lw=0.9)
    for ax in axes[1:]:
        ax.set_ylabel("")
    im = axes[3].imshow(resid[r0:r1, c0:c1], origin="lower", extent=ext, cmap=cmap, vmin=-15, vmax=15)
    _colourbar(fig, im, axes[3], "measured minus model (km/s)")
    _map_axes(axes[3], "d  What circular rotation leaves out")
    for ax in axes:
        _beam(ax, gr, ext[0] - 0.06 * (ext[0] - ext[1]), ext[2] + 0.05 * (ext[3] - ext[2]))
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[2], ext[3])
    tr = gr.primary
    fig.suptitle(f"{NAMES[key]}: the 'spider diagram' of a rotating disc, from {gr.n_good['gh_vpeak']:,} fitted "
                 f"21-cm line profiles", x=0.01, ha="left")
    P.caption(fig, f"THINGS natural-weighted VLA cube (Walter et al. 2008). b: velocity where each pixel's Gauss-Hermite "
              f"fit peaks; contours every 20 km/s, the thick one at the systemic velocity {vsys:.1f} km/s. a: orange "
              f"ellipses are every 8th of our {len(tr.radius)} tilted rings. c: our rings' circular-rotation model. d: b minus c. "
              f"Ellipse in the corner: the {gr.beam[0]:.0f}\" x {gr.beam[1]:.0f}\" beam. East is left, north up.")
    return P.save(fig, FIGURE_DIR / "01_spider_diagram")


# --------------------------------------------------------------------------
# 02 Five ways to read a line profile
# --------------------------------------------------------------------------

def _example_profiles(key, n=2):
    """High-S/N pixels whose Gauss-Hermite fits are strongly skewed."""
    import things
    g = galaxy_by_name(NAMES[key])
    maps = things.velocity_fields(g)
    cube = things.load_cube(g)
    ok = things.quality_mask(maps, cube)
    snr = np.where(ok, maps["gh_amp"] / cube.noise, 0)
    h3 = np.where(ok, np.abs(maps["gh_h3"]), 0)
    score = np.where((snr > 15) & (h3 > 0.15), snr * h3, 0)
    picks = []
    flat = np.argsort(score.ravel())[::-1]
    for idx in flat:
        r, c = np.unravel_index(idx, score.shape)
        if score[r, c] <= 0:
            break
        if all(np.hypot(r - rr, c - cc) > 60 for rr, cc in picks):
            picks.append((r, c))
        if len(picks) == n:
            break
    out = []
    for r, c in picks:
        spec = np.asarray(cube.data[:, r, c], dtype=float)
        out.append({"row": r, "col": c, "v": cube.velocity, "spec": spec, "noise": cube.noise,
                    **{k: float(maps[k][r, c]) for k in ("gh_amp", "gh_v0", "gh_sigma", "gh_h3", "gh_vpeak",
                                                         "gauss_amp", "gauss_v", "gauss_sigma", "peak_v",
                                                         "mom1_things")}})
    return out


def fig_recipes(grs, profile_key="NGC2903"):
    import things
    profiles = _example_profiles(profile_key, 2)
    fig = plt.figure(figsize=(14, 6.8), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.25, 1, 1])
    marker_colours = {"gh_vpeak": P.INK, "gh_v0": P.ORANGE, "gauss_v": P.BLUE, "peak_v": P.AQUA,
                      "mom1_things": P.VIOLET}
    for i, pr in enumerate(profiles):
        ax = fig.add_subplot(gs[i, 0])
        v = pr["v"]
        ax.step(v, pr["spec"] * 1e3, where="mid", color=P.MUTED, lw=1.0, label="data")
        vv = np.linspace(v.min(), v.max(), 800)
        ax.plot(vv, things.gauss_hermite(vv, pr["gh_amp"], pr["gh_v0"], pr["gh_sigma"], pr["gh_h3"]) * 1e3,
                color=P.INK, lw=1.6, label=f"Gauss-Hermite fit, h3 = {pr['gh_h3']:+.2f}")
        ax.plot(vv, things.gauss_hermite(vv, pr["gauss_amp"], pr["gauss_v"], pr["gauss_sigma"], 0) * 1e3,
                color=P.BLUE, lw=1.0, alpha=0.8, label="single-Gaussian fit")
        for k, lab in KIND_LABELS.items():
            ax.axvline(pr[k], color=marker_colours[k], lw=1.1, ymin=0.0, ymax=0.82,
                       label=f"{lab}: {pr[k]:.1f}")
        top = np.nanmax(pr["spec"]) * 1e3
        ax.set_xlim(pr["gh_vpeak"] - 6 * pr["gh_sigma"], pr["gh_vpeak"] + 6 * pr["gh_sigma"])
        ax.set_ylim(-3 * pr["noise"] * 1e3, top * 1.25)
        ax.set_xlabel("heliocentric velocity (km/s)")
        ax.set_ylabel("flux density (mJy/beam)")
        ax.set_title(f"{'ab'[i]}  A lopsided 21-cm profile in {NAMES[profile_key]}", fontsize=10)
        ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=7, title="velocities (km/s)",
                  title_fontsize=7)

    kinds = list(KIND_LABELS)
    y = np.arange(len(kinds))[::-1]
    keys = [k for k in ORDER if k in grs]
    offs = {key: (j - (len(keys) - 1) / 2) * 0.16 for j, key in enumerate(keys)}
    for col, (metric, xlabel, title) in enumerate((
            ("rms_diff", "rms difference over the whole curve (km/s)", "c  How closely each recipe\nreproduces de Blok et al."),
            ("inner_median_diff", "median offset in the inner 15% of the curve (km/s)", "d  ...and the bias where\nthe curve rises"))):
        ax = fig.add_subplot(gs[:, col + 1])
        for key in keys:
            comp = grs[key].comparisons["deblok"]
            ax.scatter([comp[k][metric] for k in kinds], y + offs[key], s=46, color=P.GALAXY_COLOURS[key],
                       label=NAMES[key], zorder=3, edgecolor=P.SURFACE, linewidth=1.5)
        if metric == "inner_median_diff":
            ax.axvline(0, color=P.AXIS, lw=0.9)
        ax.set_yticks(y, [KIND_LABELS[k] for k in kinds] if col == 0 else [""] * len(kinds))
        ax.set_xlabel(xlabel)
        ax.set_title(title, fontsize=10)
        ax.grid(True, axis="x")
        ax.grid(False, axis="y")
        ax.set_ylim(-0.6, len(kinds) - 0.4)
        if col == 0:
            ax.set_xlim(0, None)
            ax.legend(loc="upper right", fontsize=8)
    P.caption(fig, "a, b: two of the most skewed high-S/N profiles in the barred spiral NGC 2903, where gas at "
                   "different speeds shares a line of sight; each recipe reports a different 'velocity'. c, d: our "
                   "full pipeline run on each recipe's velocity field and compared point by point with de Blok et al. "
                   "(2008). They fitted Gauss-Hermite profiles without saying which velocity they read off: the peak "
                   "reproduces their curves (as does the brightest channel), the centre parameter sits 0.5-1.5 km/s "
                   "low, and the first moment drags the inner curve down by up to 12 km/s.")
    return P.save(fig, FIGURE_DIR / "02_velocity_recipes")


# --------------------------------------------------------------------------
# 03 Tilted rings: warps and the two halves of each disc
# --------------------------------------------------------------------------

def fig_tilted_rings(grs):
    keys = [k for k in ORDER if k in grs]
    fig, axes = plt.subplots(3, len(keys), figsize=(4.2 * len(keys), 9.5), sharex="col",
                             constrained_layout=True, gridspec_kw={"height_ratios": [1, 1, 1.6]})
    axes = np.atleast_2d(axes).reshape(3, len(keys))
    for j, key in enumerate(keys):
        gr = grs[key]
        tr = gr.primary
        r = gr.r_kpc()
        ax = axes[0, j]
        ax.errorbar(r, tr.free_pa, yerr=tr.free_pa_err, fmt="o", ms=2.5, color=P.MUTED, elinewidth=0.6,
                    label="each ring free")
        ax.plot(r, tr.pa, color=P.INK, lw=1.6, label="adopted (smoothed)")
        ax.set_title(f"{NAMES[key]}", fontsize=11)
        if j == 0:
            ax.set_ylabel("position angle (deg)")
            ax.legend(fontsize=7)
        lo, hi = np.nanpercentile(np.concatenate([tr.free_pa, tr.pa]), [2, 98])
        ax.set_ylim(lo - 4, hi + 4)
        ax = axes[1, j]
        ax.errorbar(r, tr.free_incl, yerr=tr.free_incl_err, fmt="o", ms=2.5, color=P.MUTED, elinewidth=0.6)
        ax.plot(r, tr.incl, color=P.INK, lw=1.6)
        lo, hi = np.nanpercentile(np.concatenate([tr.free_incl, tr.incl]), [2, 98])
        ax.set_ylim(lo - 4, hi + 4)
        if j == 0:
            ax.set_ylabel("inclination (deg)")
        ax = axes[2, j]
        ax.plot(r, tr.vrot_app, color=P.C_APPROACHING, lw=1.2, label="approaching half")
        ax.plot(r, tr.vrot_rec, color=P.C_RECEDING, lw=1.2, label="receding half")
        ax.errorbar(r, tr.vrot, yerr=tr.vrot_err, fmt="o", ms=3, color=P.INK, ecolor=P.MUTED, elinewidth=0.7,
                    label="both halves", zorder=3)
        ax.set_xlabel(f"radius (kpc, at {gr.distance_mpc:g} Mpc)")
        ax.set_ylim(0, np.nanmax(tr.vrot + tr.vrot_err) * 1.12)
        if j == 0:
            ax.set_ylabel("rotation speed (km/s)")
            ax.legend(fontsize=7, loc="lower right")
        fixed = galaxy_by_name(NAMES[key]).geometry_fixed_inside_arcsec
        if fixed:
            for a in axes[:, j]:
                a.axvspan(0, fixed * gr.kpc_per_arcsec, color=P.C_CONTEXT, alpha=0.35, lw=0)
            axes[0, j].text(fixed * gr.kpc_per_arcsec / 2, axes[0, j].get_ylim()[1], "bar:\ngeometry\nheld",
                            fontsize=7, color=P.INK_SECONDARY, ha="center", va="top")
    fig.suptitle("Tilted rings: each ring's own tilt and orientation, and each half of the disc separately",
                 x=0.01, ha="left")
    P.caption(fig, "Grey points: position angle and inclination fitted freely ring by ring (centre and systemic "
                   "velocity fixed); black: the smoothed profiles adopted for the rotation curve. Bottom: rotation "
                   "with the adopted rings, error bars following de Blok et al. (2008): the dispersion around each "
                   "ring plus a quarter of the difference between the halves.")
    return P.save(fig, FIGURE_DIR / "03_tilted_rings")


# --------------------------------------------------------------------------
# 04 Our rotation curves against the published ones
# --------------------------------------------------------------------------

def fig_rotation_curves(grs):
    keys = [k for k in ORDER if k in grs]
    fig = plt.figure(figsize=(13, 9.5), constrained_layout=True)
    outer = fig.add_gridspec(2, 2)
    for i, key in enumerate(keys):
        gr = grs[key]
        tr = gr.primary
        sub = outer[i // 2, i % 2].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
        ax, axr = fig.add_subplot(sub[0]), fig.add_subplot(sub[1])
        r = gr.r_kpc()
        db, sp = gr.deblok, gr.sparc
        ax.errorbar(sp.r_kpc, sp.v, yerr=sp.err, fmt="s", ms=3.5, color=P.C_OLDER, mfc="none", elinewidth=0.7,
                    label=f"SPARC (D4): {gr.comparisons['sparc']['refs'].replace(',', ', ')}")
        ax.errorbar(db.r_kpc, db.v, yerr=db.err, fmt="o", ms=3.5, color=P.C_PUBLISHED, mfc="none",
                    elinewidth=0.6, label="de Blok et al. 2008 (THINGS)")
        ax.errorbar(r, tr.vrot, yerr=tr.vrot_err, fmt="o", ms=2.6, color=P.INK, ecolor=P.MUTED, elinewidth=0.6,
                    label="ours (THINGS cube)", zorder=4)
        ax.set_ylim(0, np.nanmax(np.concatenate([tr.vrot, db.v])) * 1.18)
        ax.set_xlim(0, max(np.nanmax(r), db.r_kpc.max(), sp.r_kpc.max()) * 1.03)
        ax.set_ylabel("rotation speed (km/s)")
        ax.set_title(f"{'abcd'[i]}  {NAMES[key]}", fontsize=10)
        ax.tick_params(labelbottom=False)
        ax.legend(fontsize=7, loc="lower right")
        ok = np.isfinite(tr.vrot)
        vi = np.interp(db.r_kpc, r[ok], tr.vrot[ok], left=np.nan, right=np.nan)
        axr.axhline(0, color=P.AXIS, lw=0.8)
        axr.fill_between(db.r_kpc, -db.err, db.err, color=P.C_PUBLISHED, alpha=0.12, lw=0)
        axr.plot(db.r_kpc, vi - db.v, "o", ms=2.6, color=P.INK)
        c = gr.comparisons["deblok"]["gh_vpeak"]
        axr.text(0.99, 0.92, f"median {c['median_diff']:+.2f}, rms {c['rms_diff']:.2f} km/s",
                 transform=axr.transAxes, ha="right", va="top", fontsize=7, color=P.INK_SECONDARY)
        lim = max(6, np.nanpercentile(np.abs(vi - db.v), 98) * 1.3)
        axr.set_ylim(-lim, lim)
        axr.set_xlim(ax.get_xlim())
        axr.set_ylabel("ours − de Blok")
        axr.set_xlabel(f"radius (kpc; all at {gr.distance_mpc:g} Mpc)")
    fig.suptitle("Rotation curves measured from the raw cubes, over the published THINGS curves and older data",
                 x=0.01, ha="left")
    P.caption(fig, "Black: our tilted-ring curves from the Gauss-Hermite velocity fields. Violet: de Blok et al. "
                   "(2008), read point by point from the vector figures of their paper (shaded band below: their "
                   "error bars). Aqua: the curves SPARC compiled from older radio and optical data, radii rescaled to "
                   "the THINGS distance. NGC 2903's inner rise differs between telescopes because of its bar.")
    return P.save(fig, FIGURE_DIR / "04_rotation_curves")


# --------------------------------------------------------------------------
# 05 What circular rotation leaves out: residual fields and harmonics
# --------------------------------------------------------------------------

def fig_noncircular(grs):
    keys = [k for k in ORDER if k in grs]
    fig = plt.figure(figsize=(4.3 * len(keys), 8.6), constrained_layout=True)
    gs = fig.add_gridspec(2, len(keys), height_ratios=[1.35, 1])
    cmap = P.velocity_cmap()
    for j, key in enumerate(keys):
        gr = grs[key]
        m = _maps(key)
        fr = _frame_for(gr)
        r0, r1, c0, c1 = _crop(np.isfinite(m["residual"]), 10)
        ext = _extent(gr, fr, r0, r1, c0, c1)
        ax = fig.add_subplot(gs[0, j])
        im = ax.imshow(m["residual"][r0:r1, c0:c1], origin="lower", extent=ext, cmap=cmap, vmin=-20, vmax=20)
        _map_axes(ax, f"{'abcd'[j]}  {NAMES[key]}")
        if j:
            ax.set_ylabel("")
        _beam(ax, gr, ext[0] - 0.06 * (ext[0] - ext[1]), ext[2] + 0.05 * (ext[3] - ext[2]))
        fixed = galaxy_by_name(NAMES[key]).geometry_fixed_inside_arcsec
        if fixed:
            t = np.linspace(0, 2 * np.pi, 100)
            ax.plot(fixed / 60 * np.cos(t), fixed / 60 * np.sin(t), color=P.INK, lw=0.8)
            ax.text(0, fixed / 60 * 1.08, "bar", ha="center", va="bottom", fontsize=7, color=P.INK)
        if j == len(keys) - 1:
            _colourbar(fig, im, ax, "measured minus circular model (km/s)")

        ax = fig.add_subplot(gs[1, j])
        h = gr.extra.get("harmonics_reswri", gr.harmonics)
        r = h.radius * gr.kpc_per_arcsec
        ax.plot(r, h.a_r, "o", ms=2.6, color=P.INK, label="ours, ring by ring")
        hs = gr.harmonic_summary
        tr = hs["trachternach"]
        ours = hs["reswri_rings"]
        if tr:
            ax.axhspan(tr["A_r_q"][0], tr["A_r_q"][1], color=P.C_PUBLISHED, alpha=0.13, lw=0)
            ax.axhline(tr["A_r"], color=P.C_PUBLISHED, lw=1.4, label=f"Trachternach et al. 2008: median {tr['A_r']:.1f}")
        ax.axhline(ours["A_r_median"], color=P.INK, lw=1.4, label=f"ours: median {ours['A_r_median']:.1f}")
        ax.axvspan(0, 1.0, color=P.C_CONTEXT, alpha=0.35, lw=0)
        ax.text(0.5, 0.97, "inner\nkpc", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=7,
                color=P.INK_SECONDARY)
        ax.set_yscale("log")
        ax.set_ylim(0.3, 100)
        ax.set_xlabel("radius (kpc)")
        if j == 0:
            ax.set_ylabel("non-circular speed A_r (km/s, in the disc)")
        ax.legend(fontsize=7, loc="upper right")
        ax.set_title(f"median inside 1 kpc: ours {ours['A_r_1kpc_median']:.1f}"
                     + (f", theirs {tr['A_r_1kpc']:.1f}" if tr and np.isfinite(tr['A_r_1kpc']) else ""),
                     fontsize=8, color=P.INK_SECONDARY, loc="left")
    fig.suptitle("What circular rotation leaves out: streaming along spiral arms and along a bar", x=0.01, ha="left")
    P.caption(fig, "Top: each galaxy's velocity field minus its tilted-ring model (same rings as the rotation curves). "
                   "Spiral arms show as arcs of +-5-10 km/s; NGC 2903's bar as a strong pattern inside the circle. "
                   "Bottom: the quadratic sum of the non-circular harmonic terms (s1, c2, s2, c3, s3) around each ring, "
                   "rings refitted without weights as GIPSY's RESWRI does, against Trachternach et al.'s median and "
                   "quartiles (violet).")
    return P.save(fig, FIGURE_DIR / "05_noncircular_motions")


# --------------------------------------------------------------------------
# 06 The cusp-core test in DDO 154
# --------------------------------------------------------------------------

def fig_cusp_core(grs, cusp):
    import literature as L
    import massmodel as MM
    gr = grs["DDO154"]
    tr = gr.primary
    r, v, e = gr.r_kpc(), tr.vrot, tr.vrot_err
    fig = plt.figure(figsize=(15, 5.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.15, 0.85])

    ax = fig.add_subplot(gs[0])
    cuspm, obsm = cusp["models"]["LCDM cusp (NFW)"], cusp["models"]["observed curve"]
    lc = cusp["lcdm_nfw"]
    rr = np.linspace(0.02, r.max(), 300)
    ax.plot(rr, MM.v_nfw(rr, lc["c"], lc["v200"]), color=P.C_MODEL, lw=1.0, alpha=0.6,
            label=f"cold-dark-matter cusp: NFW, c = {lc['c']:.1f}, matched beyond 3 kpc")
    ax.plot(cuspm["r_kpc"], cuspm["v_recovered"], "s-", ms=3.5, color=P.C_MODEL, lw=1.4,
            label="the same cusp seen through the VLA beam and our pipeline")
    iso = cusp["halo"]["ours"]["ISO"]
    ax.plot(rr, MM.v_iso(rr, *iso["params"]), color=P.C_PUBLISHED, lw=1.2,
            label=f"cored halo fit: R_c = {iso['params'][1]:.2f} kpc")
    ax.plot(obsm["r_kpc"], obsm["v_recovered"], color=P.MUTED, lw=1.0,
            label="our curve seen through the beam and pipeline")
    ax.errorbar(r, v, yerr=e, fmt="o", ms=3.2, color=P.INK, ecolor=P.MUTED, elinewidth=0.7, zorder=4,
                label="DDO 154, measured")
    ax.set_xlim(0, 3.0)
    ax.set_ylim(0, 50)
    ax.set_xlabel("radius (kpc)")
    ax.set_ylabel("rotation speed (km/s)")
    ex = cusp["cusp_excess"]
    ax.set_title(f"a  Inside 1 kpc a cusp would spin {ex['mean_ratio']:.1f}x faster", fontsize=10)
    ax.legend(fontsize=7, loc="lower right")

    ax = fig.add_subplot(gs[1])
    rho_r, rho = cusp["rho_ours"]
    ax.plot(rho_r, rho, "o", ms=3.4, color=P.INK, label="DDO 154, our curve", zorder=4)
    db_r, db_rho = cusp["rho_deblok"]
    ax.plot(db_r, db_rho, "o", ms=3.4, mfc="none", color=P.C_PUBLISHED, label="de Blok et al.'s curve, same code")
    oh = L.oh2011_ddo154_density()
    ax.plot(oh[:, 0], oh[:, 1], "D", ms=4.2, color=P.C_PUBLISHED, label="Oh et al. 2011's fitted points")
    rc = cuspm["r_kpc"]
    ax.plot(rc, MM.density_smoothed(rc, cuspm["v_recovered"], cusp["window"]), "s-", ms=3, lw=1.1,
            color=P.C_MODEL, label="the cusp through beam and pipeline")
    x = np.array([0.12, 0.45])
    for slope, lab in ((0, "core: rho ~ R⁰"), (-1, "cusp: rho ~ R⁻¹")):
        ax.plot(x, 0.004 * (x / 0.12) ** slope, color=P.INK_SECONDARY, lw=1.0)
        ax.text(x[-1] * 1.08, 0.004 * (x[-1] / 0.12) ** slope, lab, fontsize=7, color=P.INK_SECONDARY,
                va="center")
    ax.axvline(cusp["break_kpc"], color=P.AXIS, lw=0.8)
    ax.text(cusp["break_kpc"] * 1.05, 0.25, "break radius", fontsize=7, color=P.INK_SECONDARY, va="top")
    ax.set_xscale("log")
    ax.set_yscale("log")
    from matplotlib.ticker import FixedLocator, NullFormatter, FuncFormatter
    ax.xaxis.set_major_locator(FixedLocator([0.1, 0.2, 0.5, 1, 2, 5]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda val, _: f"{val:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(0.1, 8)
    ax.set_ylim(5e-4, 0.3)
    ax.set_xlabel("radius (kpc)")
    ax.set_ylabel("dark matter density, minimum disc (Msun/pc³)")
    ax.set_title("b  Density from the rotation curve", fontsize=10)
    ax.legend(fontsize=7, loc="lower left")

    ax = fig.add_subplot(gs[2])
    rows = [("ours", cusp["alpha_ours"], P.INK),
            ("ours, gas subtracted", cusp["alpha_ours_dm"], P.INK),
            ("de Blok's curve,\nour code", cusp["alpha_deblok"], P.C_PUBLISHED),
            ("Oh et al. 2011", {"alpha": L.OH2011_DDO154_ALPHA[0], "alpha_err": L.OH2011_DDO154_ALPHA[1]},
             P.C_PUBLISHED),
            ("cusp, as built", cuspm["alpha_input"], P.C_MODEL),
            ("cusp, through beam\nand pipeline", cuspm["alpha_recovered"], P.C_MODEL)]
    y = np.arange(len(rows))[::-1]
    for yi, (lab, a, col) in zip(y, rows):
        ax.errorbar(a["alpha"], yi, xerr=a["alpha_err"], fmt="o", ms=6, color=col, elinewidth=1.2, capsize=0)
        ax.text(a["alpha"], yi + 0.22, f"{a['alpha']:+.2f}", ha="center", fontsize=7, color=P.INK_SECONDARY)
    ax.set_yticks(y, [r_[0] for r_ in rows])
    ax.axvline(0, color=P.AXIS, lw=0.8)
    ax.axvline(-1, color=P.AXIS, lw=0.8)
    ax.set_xlim(-1.7, 0.6)
    ax.set_xlabel("inner density slope alpha (inside 0.68 kpc)")
    ax.set_title("c  Core (0) or cusp (-1)?", fontsize=10)
    ax.grid(False, axis="y")
    fig.suptitle("DDO 154: the cusp-core problem, measured directly and checked against beam smearing", x=0.01,
                 ha="left")
    P.caption(fig, "a: our rotation curve against a cold-dark-matter (NFW) halo with the concentration simulations "
                   f"predict for its mass (Dutton & Maccio 2014), matched to the curve beyond 3 kpc. We built model "
                   f"data cubes of both curves (labelled simulation), smoothed them with the real beam, added noise and "
                   f"ran the identical pipeline: smearing pulls the cusp down by 1-5 km/s, yet it stays "
                   f"{ex['mean_excess_kms']:.1f} km/s above the data on average (chi2 = {ex['chi2']:.0f} for {ex['n']} "
                   f"points). b, c: densities from rho = (2 V/R dV/dR + V²/R²) / 4 pi G (all mass dark: the cusp's "
                   f"best case), slopes fitted inside Oh et al.'s break radius.")
    return P.save(fig, FIGURE_DIR / "06_cusp_core")


# --------------------------------------------------------------------------
# MaNGA helpers
# --------------------------------------------------------------------------

def _bin_map(ax, b, which, vlim=None, title=None):
    """A binned MaNGA velocity field, interpolated across the bin centres."""
    import manga as MG
    sel, vel, _ = MG._component(b, which) if which in ("star", "gas") else (None, None, None)
    sel = sel & np.isfinite(vel)
    x, y = b["x_arcsec"][sel], b["y_arcsec"][sel]
    v = vel[sel] - np.median(vel[sel])
    vlim = vlim or float(np.nanpercentile(np.abs(v), 95))
    tpc = ax.tripcolor(x, y, v, shading="gouraud", cmap=P.velocity_cmap(), vmin=-vlim, vmax=vlim)
    ax.set_aspect("equal")
    lim = 1.05 * max(np.max(np.abs(x)), np.max(np.abs(y)))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    if title:
        ax.set_title(title, fontsize=8, loc="center")
    return tpc, vlim


# --------------------------------------------------------------------------
# 07 Stars lag the gas: asymmetric drift
# --------------------------------------------------------------------------

def fig_drift(manga_summary):
    import manga_io as MIO
    gals = manga_summary["drift"]
    fig = plt.figure(figsize=(13, 3.6 * len(gals)), constrained_layout=True)
    gs = fig.add_gridspec(len(gals), 3, width_ratios=[1, 1, 1.9])
    for i, g in enumerate(gals):
        b = MIO.load_bins(g["plateifu"])
        rows = g["rows"]
        _, vl = _bin_map(fig.add_subplot(gs[i, 1]), b, "gas", title="ionized gas (H-alpha)")
        tpc, _ = _bin_map(fig.add_subplot(gs[i, 0]), b, "star", vlim=vl, title=f"stars\n{g['plateifu']}")
        ax = fig.add_subplot(gs[i, 2])
        r = np.array([row["r_re"] for row in rows])

        def col(key):
            return np.array([row.get(key, np.nan) for row in rows])

        ax.plot(r, col("v_gas"), "o-", ms=4, color=P.C_GAS, lw=1.5, label="gas")
        ax.plot(r, col("v_star"), "o-", ms=4, color=P.C_STARS, lw=1.5, label="stars, as measured")
        ax.plot(r, col("v_star_corrected"), "o-", ms=4, mfc=P.SURFACE, color=P.C_STARS, lw=1.0,
                label="stars + asymmetric drift")
        ax.plot(r, col("sigma_r"), "s-", ms=3, color=P.MUTED, lw=1.0, label="stars' radial dispersion, from MaNGA")
        ax.plot(r, col("sigma_r_implied"), "s-", ms=3, mfc=P.SURFACE, color=P.MUTED, lw=1.0,
                label="radial dispersion the lag implies")
        ax.set_xlim(0, 2.0)
        ax.set_ylim(0, None)
        ax.set_ylabel("km/s")
        ax.set_title(f"{g['plateifu']}: inclination {g['incl']:.0f} deg, log M* = {g['log_mstar']:.1f}",
                     fontsize=9)
        if i == len(gals) - 1:
            ax.set_xlabel("radius (effective radii)")
        if i == 0:
            ax.legend(fontsize=7, loc="lower right")
    fig.suptitle("Stars rotate slower than gas, but by less than their measured random motions predict", x=0.01,
                 ha="left")
    P.caption(fig, "MaNGA DR17 Voronoi-binned maps (systemic velocity removed; same colour scale for stars and gas "
                   "in each galaxy). Right: rotation in rings with the gas's kinematic position angle and "
                   "an inclination from the photometric axis ratio. Open orange: the stars corrected with the Jeans "
                   "equation for an exponential disc (Binney & Tremaine eq. 4.228), assuming sigma_z = 0.6 sigma_R and "
                   "sigma_phi² = sigma_R²/2; it overshoots the gas. Run backwards, the measured lag implies radial "
                   "dispersions about half of MaNGA's, whose disc values (45-70 km/s along the line of sight) sit at "
                   "the instrument's ~70 km/s resolution.")
    return P.save(fig, FIGURE_DIR / "07_asymmetric_drift")


# --------------------------------------------------------------------------
# 08 Gas that rotates against its stars
# --------------------------------------------------------------------------

#: Zhou et al. (2022), MaNGA MPL-10: misaligned (> 30 deg) fractions among
#: emission-line galaxies, after visual inspection, by SED-based class.
ZHOU2022 = {"blue": ("star-forming", 0.024), "green": ("green valley", 0.073), "red": ("quiescent", 0.156),
          "all": ("all", 456 / 6958)}


def _binomial_ci(k, n, z=1.0):
    """Wilson interval for a proportion."""
    if n == 0:
        return np.nan, np.nan
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return centre - half, centre + half


def fig_counter_rotation(manga_summary, n_gallery=8):
    import manga_io as MIO
    cen = manga_summary["census"]
    fig = plt.figure(figsize=(15, 11), constrained_layout=True)
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 0.06, 1.15])

    ax = fig.add_subplot(gs[0, 0])
    st = cen["by_sample"]["regular"]
    edges = np.arange(0, 181, 10)
    ax.bar(edges[:-1] + 5, st["hist"], width=9.2, color=P.C_OBSERVED)
    ax.set_yscale("log")
    for t, lab in ((30, "misaligned →"), (150, "counter-\nrotating →")):
        ax.axvline(t, color=P.AXIS, lw=0.9)
        ax.text(t + 2, max(st["hist"]) * 0.9, lab, fontsize=7, color=P.INK_SECONDARY, va="top")
    ax.set_xlim(0, 180)
    ax.set_xticks([0, 30, 60, 90, 120, 150, 180])
    ax.set_xlabel("angle between the stars' and the gas's rotation (deg)")
    ax.set_ylabel("galaxies")
    ax.set_title(f"a  {st['n']:,} galaxies with both rotations measured", fontsize=10)

    ax = fig.add_subplot(gs[0, 1])
    classes = ["blue", "green", "red"]
    xs = np.arange(len(classes))
    for off, sample, colr, lab in ((-0.12, "measured", P.MUTED, "ours, automatic cuts (as Zhou et al.)"),
                                   (0.12, "regular", P.INK, "ours, regular rotators only")):
        bc = cen["by_sample"][sample]["by_colour"]
        f = [bc[c]["frac_misaligned"] for c in classes]
        ci = [_binomial_ci(bc[c]["misaligned"], bc[c]["n"]) for c in classes]
        ax.errorbar(xs + off, f, yerr=[[fi - lo for fi, (lo, hi) in zip(f, ci)], [hi - fi for fi, (lo, hi) in zip(f, ci)]],
                    fmt="o", ms=6, color=colr, label=lab, capsize=0)
    ax.plot(xs + 0.3, [ZHOU2022[c][1] for c in classes], "D", ms=6, color=P.C_PUBLISHED,
            label="Zhou et al. 2022, MaNGA (by star formation rate, vetted by eye)")
    ax.errorbar(xs[2] + 0.42, 0.36, yerr=0.05, fmt="s", ms=6, mfc="none", color=P.C_PUBLISHED, capsize=0,
                label="ATLAS3D: fast-rotating early types with gas (Davis et al. 2011)")
    ax.set_xticks(xs, ["blue cloud\n(star-forming)", "green valley", "red sequence\n(quiescent)"])
    ax.set_ylabel("fraction with gas misaligned by > 30 deg")
    ax.set_ylim(0, None)
    ax.legend(fontsize=7, loc="upper left")
    ax.set_title("b  Captured gas is commonest in red galaxies", fontsize=10)

    ax = fig.add_subplot(gs[0, 2])
    for sample, colr, lab in (("measured", P.MUTED, "automatic cuts"), ("regular", P.INK, "regular rotators")):
        sm = cen["by_sample"][sample]
        mb = np.array(sm["mass_bins"])
        mid = 0.5 * (mb[1:] + mb[:-1])
        f = np.array([np.nan if v is None else v for v in sm["frac_misaligned_by_mass"]], dtype=float)
        ax.plot(mid, f, "o-", ms=5, color=colr, lw=1.3, label=lab)
    ax.set_xlabel("stellar mass, log M* (Msun)")
    ax.set_ylabel("fraction misaligned by > 30 deg")
    ax.set_ylim(0, None)
    ax.legend(fontsize=7)
    ax.set_title("c  ...and in more massive galaxies", fontsize=10)

    sub = gs[2, :].subgridspec(2, n_gallery, hspace=0.05, wspace=0.05)
    gallery = cen["counter_rotators"][:n_gallery]
    for j, g in enumerate(gallery):
        b = MIO.load_bins(g["plateifu"])
        a1, a2 = fig.add_subplot(sub[0, j]), fig.add_subplot(sub[1, j])
        _, vl = _bin_map(a1, b, "star", title=f"{g['plateifu']}\nstars")
        _bin_map(a2, b, "gas", vlim=vl, title=f"gas  ({g['dpa']:.0f} deg)")
    head = fig.add_subplot(gs[1, :])
    head.axis("off")
    head.text(0.0, 0.0, f"d  {len(gallery)} of the {len(cen['counter_rotators'])} counter-rotators among the regular "
                        "rotators (the cleanest): blue approaching, red receding, one scale for stars (top) and gas "
                        "(bottom) in each galaxy", fontsize=10, fontweight="bold", transform=head.transAxes)
    fig.suptitle("Gas that rotates against its stars in MaNGA DR17: the signature of gas captured from outside",
                 x=0.01, ha="left")
    P.caption(fig, "Kinematic position angles fitted to stars and H-alpha gas within 1.5 effective radii with "
                   "Krajnovic et al.'s method (our faster implementation of pafit, identical results), with the sense "
                   "of rotation added. Colour classes from NUV - r (< 4 blue, 4-5 green, > 5 red; NASA-Sloan Atlas). "
                   "Error bars: 68% binomial intervals. Zhou et al. classified by SED star formation rate and inspected "
                   "their candidates by eye, removing mergers and irregulars.")
    return P.save(fig, FIGURE_DIR / "08_counter_rotation")


# --------------------------------------------------------------------------
# 00 Graded against the literature
# --------------------------------------------------------------------------

def fig_scorecard(results):
    """Signed percentage differences from published values (D4's layout)."""
    gal = results["galaxies"]
    rows = []  # (label, ours_minus_published_pct, band_lo_pct, band_hi_pct, value text)
    for key in ORDER:
        if key not in gal:
            continue
        g = gal[key]
        d = g["vs_deblok"]["gh_vpeak"]
        vmax = g["v_max"]
        rows.append((f"{g['name']}: rotation curve vs de Blok et al. 2008", 100 * d["median_diff"] / vmax,
                     None, None, f"median {d['median_diff']:+.2f} km/s, rms {d['rms_diff']:.2f}"))
    for key in ORDER:
        if key not in gal:
            continue
        g = gal[key]
        hi = g["hi"]
        rows.append((f"{g['name']}: HI flux vs THINGS global spectrum", 100 * (hi["mom0_flux_jykms"] / hi["walter_flux_jykms"] - 1),
                     None, None, f"{hi['mom0_flux_jykms']:.0f} vs {hi['walter_flux_jykms']:.0f} Jy km/s"))
    for key in ORDER:
        if key not in gal:
            continue
        g = gal[key]
        h, t = g["harmonics"]["reswri_rings"], g["harmonics"]["trachternach"]
        rows.append((f"{g['name']}: non-circular A_r vs Trachternach et al.", 100 * (h["A_r_median"] / t["A_r"] - 1),
                     100 * (t["A_r_q"][0] / t["A_r"] - 1), 100 * (t["A_r_q"][1] / t["A_r"] - 1),
                     f"{h['A_r_median']:.2f} vs {t['A_r']:.2f} km/s"))
    c = results["ddo154_cusp"]["halo"]
    oh = {"ISO": ((34.5e-3, 1.8e-3), (1.30, 0.04)), "NFW": ((5.3, 0.5), (58.3, 4.2))}
    for curve, label in (("ours", "our curve"), ("deblok", "de Blok's curve, our code")):
        for kind, names in (("ISO", ("core density", "core radius")), ("NFW", ("concentration", "V200"))):
            p = c[curve][kind]["params"]
            for j in range(2):
                ref, err = oh[kind][j]
                rows.append((f"DDO 154 {kind} {names[j]} ({label}) vs Oh et al. 2011", 100 * (p[j] / ref - 1),
                             -100 * err / ref, 100 * err / ref, f"{p[j]:.3g} vs {ref:.3g}"))
    if "manga" in results:
        s = results["manga"]["census"]["by_sample"]
        for sample, lab in (("regular", "regular rotators"),):
            f = s[sample]["frac_misaligned"]
            ref = ZHOU2022["all"][1]
            err = np.sqrt(456) / 6958
            rows.append((f"MaNGA misaligned fraction ({lab}) vs Zhou et al. 2022", 100 * (f / ref - 1),
                         -100 * err / ref, 100 * err / ref, f"{100 * f:.1f}% vs {100 * ref:.1f}%"))

    fig, ax = plt.subplots(figsize=(13, 0.42 * len(rows) + 1.8), constrained_layout=True)
    y = np.arange(len(rows))[::-1]
    for yi, (lab, val, lo, hi, txt) in zip(y, rows):
        if lo is not None:
            ax.add_patch(plt.Rectangle((lo, yi - 0.36), hi - lo, 0.72, color=P.C_PUBLISHED, alpha=0.13, lw=0))
        ax.barh(yi, val, height=0.5, color=P.C_MODEL)
        ax.text(val + (1.2 if val >= 0 else -1.2), yi, f"{val:+.1f}%  ({txt})", va="center",
                ha="left" if val >= 0 else "right", fontsize=7.5, color=P.INK_SECONDARY)
    ax.axvline(0, color=P.INK, lw=1.0)
    ax.set_yticks(y, [r[0] for r in rows])
    ax.set_xlabel("difference from the published value (%)")
    lim = max(60, max(abs(r[1]) for r in rows) * 1.5)
    ax.set_xlim(-lim, lim)
    ax.grid(True, axis="x")
    ax.grid(False, axis="y")
    ax.set_title("Our numbers against the literature", loc="left")
    P.caption(fig, "Violet bands: the published uncertainty or quartile range, where one is quoted. "
                   "Rotation-curve rows: median of (ours − published) as a percentage of the galaxy's peak speed. "
                   "HI flux: the THINGS moment-0 map summed, against Walter et al.'s global spectrum. Non-circular "
                   "motions: median over rings of the harmonic amplitude, rings refitted unweighted as RESWRI does. "
                   "DDO 154 halo fits use minimum disc, as Oh et al. did. MaNGA: our objective 'regular rotator' cut "
                   "stands in for Zhou et al.'s visual inspection" + (
                       f" (without any vetting the fraction is "
                       f"{100 * results['manga']['census']['by_sample']['measured']['frac_misaligned']:.1f}%)."
                       if "manga" in results else "."))
    return P.save(fig, FIGURE_DIR / "00_comparison")
