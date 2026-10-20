"""Static figures. Each function reads cached products and writes figures/<name>.png."""

from __future__ import annotations

import json

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patheffects as pe
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm
from matplotlib.patches import Circle

import catalogs
import config as C
import grids
import plotstyle as ps
from dustmap import to_av

ps.use_style()
HALO = [pe.withStroke(linewidth=2.4, foreground=ps.SURFACE)]

#: Curated labels for the face-on map: name -> (x, y, text offset in points).
FACE_LABELS = {          # Corona Australis and the Coalsack overlap Ophiuchus/Lupus and Chamaeleon from above
    "Taurus": (-10, -14), "Perseus": (4, 7), "California": (-60, 2), "Orion": (8, -6),
    "Mon R2": (6, -10), "Mon OB1": (6, 2), "Vela": (6, -2), "Chamaeleon": (6, -8), "Lupus": (8, -10),
    "Ophiuchus": (8, 5), "Aquila Rift": (8, 2), "Serpens": (6, 4),
    "Cepheus": (6, 4), "North America": (6, -2), "Cygnus X": (6, -12),
}


def _extent(xs, ys):
    hx, hy = (xs[1] - xs[0]) / 2, (ys[1] - ys[0]) / 2
    return (xs[0] - hx, xs[-1] + hx, ys[0] - hy, ys[-1] + hy)


def _sun(ax, **kw):
    ax.plot([0], [0], **{**ps.SUN_STYLE, **kw})


def _colorbar(fig, im, ax, label, **kw):
    kw.setdefault("pad", 0.02)
    cb = fig.colorbar(im, ax=ax, **kw)
    cb.set_label(label, color=ps.INK_SECONDARY)
    cb.outline.set_edgecolor(ps.AXIS)
    cb.ax.tick_params(colors=ps.MUTED, labelsize=7.5)
    return cb


def mollweide(ax, hmap, *, nest=False, cmap=None, norm=None, xsize=1400):
    """Draw a HEALPix map in Mollweide projection on an ordinary axes (l = 0 in the centre,
    longitude increasing to the left), with faint graticule lines every 60/30 degrees."""
    proj = hp.projector.MollweideProj(xsize=xsize)
    nside = hp.npix2nside(len(hmap))
    img = proj.projmap(np.asarray(hmap, dtype=float), lambda x, y, z: hp.vec2pix(nside, x, y, z, nest=nest))
    img = np.where(img == hp.UNSEEN, np.nan, img)
    cmap = cmap.with_extremes(bad=(0, 0, 0, 0))
    ax.imshow(img, origin="lower", cmap=cmap, norm=norm, extent=(-2, 2, -1, 1), interpolation="none")
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(2 * np.cos(t), np.sin(t), color=ps.AXIS, lw=0.8)
    for lat in (-60, -30, 0, 30, 60):
        theta = np.radians(lat)
        aux = theta
        for _ in range(30):                       # Mollweide auxiliary angle
            aux -= (2 * aux + np.sin(2 * aux) - np.pi * np.sin(theta)) / (2 + 2 * np.cos(2 * aux) + 1e-12)
        ax.plot([-2 * np.cos(aux), 2 * np.cos(aux)], [np.sin(aux)] * 2, color=ps.AXIS, lw=0.35, alpha=0.6)
    yy = np.linspace(-np.pi / 2, np.pi / 2, 200)
    aux = yy.copy()
    for _ in range(30):
        aux -= (2 * aux + np.sin(2 * aux) - np.pi * np.sin(yy)) / (2 + 2 * np.cos(2 * aux) + 1e-12)
    for lon in (-120, -60, 0, 60, 120):
        ax.plot(-2 * np.radians(lon) / np.pi * np.cos(aux), np.sin(aux), color=ps.AXIS, lw=0.35, alpha=0.6)
    for lon, txt in ((0, "l = 0°"), (90, "90°"), (-90, "270°")):
        ax.annotate(txt, (-lon / 90.0, 0), xytext=(0, -9), textcoords="offset points", ha="center",
                    fontsize=6.5, color=ps.MUTED)
    ax.set_xlim(-2.02, 2.02)
    ax.set_ylim(-1.02, 1.02)
    ax.set_aspect("equal")
    ax.axis("off")
    return ax.images[-1]


def _labels(ax, rmax=1230, which=None):
    lab = catalogs.cloud_labels().set_index("complex")
    for name, off in FACE_LABELS.items():
        if (which and name not in which) or name not in lab.index:
            continue
        x, y = lab.loc[name, "x"], lab.loc[name, "y"]
        if np.hypot(x, y) > rmax:
            continue
        ax.plot(x, y, "o", ms=2.4, color=ps.ORANGE, mec="none")
        ax.annotate(name, (x, y), xytext=off, textcoords="offset points", fontsize=6.8,
                    color=ps.INK_SECONDARY, path_effects=HALO)


def projections(labelled: bool = True) -> str:
    """Face-on and edge-on views of the main map (the layout of Edenhofer et al. 2024, Fig. 4)."""
    g = grids.build_main()
    xs, zs = g["xs"], g["zs"]
    zsel = np.abs(zs) <= 400
    pxy = to_av(g["proj_xy"])
    pxz = to_av(g["proj_xz"])[:, zsel]
    pyz = to_av(g["proj_yz"])[:, zsel]
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    vmax_xy = np.quantile(pxy[np.hypot(X, Y) < C.R_OUTER_MAIN], 0.999)
    vmax_ed = np.quantile(np.r_[pxz.ravel(), pyz.ravel()], 0.999)

    fig = plt.figure(figsize=(8.2, 12.4))
    gs = fig.add_gridspec(3, 1, height_ratios=[2.5, 0.8, 0.8], hspace=0.3)
    ax = fig.add_subplot(gs[0])
    im = ax.imshow(pxy.T, extent=_extent(xs, xs), cmap=ps.DUST, norm=Normalize(0, vmax_xy),
                   interpolation="none")
    ax.add_patch(Circle((0, 0), C.R_OUTER_MAIN, fill=False, ec=ps.AXIS, lw=0.6))
    ps.map_axes(ax)
    ax.set_aspect("equal")
    ax.set_xlabel("x [pc]  (toward the Galactic Centre →)")
    ax.set_ylabel("y [pc]  (direction of Galactic rotation →)")
    ax.set_title("Dust within 1.25 kpc of the Sun, seen from above the Galactic plane")
    _colorbar(fig, im, ax, f"A_V through |z| ≤ 400 pc [mag]  (linear, saturates at {vmax_xy:.2f})",
              shrink=0.85, extend="max")
    _sun(ax)
    if labelled:
        _labels(ax)
    for k, (data, label, xl) in enumerate(((pxz, "side view, integrated over y", "x [pc]"),
                                           (pyz, "side view, integrated over x", "y [pc]"))):
        ax2 = fig.add_subplot(gs[1 + k])
        im2 = ax2.imshow(data.T, extent=_extent(xs, zs[zsel]), cmap=ps.DUST,
                         norm=Normalize(0, vmax_ed), interpolation="none", aspect="equal")
        ps.map_axes(ax2)
        ax2.set_xlabel(xl)
        ax2.set_ylabel("z [pc]")
        ax2.set_title(label, fontsize=9.5)
        _sun(ax2, ms=8)
        _colorbar(fig, im2, ax2, f"A_V [mag] (sat. {vmax_ed:.1f})", extend="max")
    ps.caption(fig, "Edenhofer et al. (2024) posterior mean, resampled onto 2.5 pc cells (log-density "
               "interpolation, as the authors' interp2box.py) and integrated; A_V = 2.8 × ZGR23 extinction. "
               "The map starts at 69 pc (the hole around the Sun) and ends at 1.25 kpc."
               + (" Orange dots: median positions of the Zucker et al. (2020) cloud sightlines." if labelled else ""),
               y=0.065)
    name = "map_projections_labelled" if labelled else "map_projections"
    ps.save(fig, C.FIGURE_DIR / name)
    return name


def fig4_check() -> str:
    """Our face-on map against the published Fig. 4, read back pixel by pixel."""
    d = np.load(C.CACHE_DIR / "fig4_pixels.npz")
    val = json.loads((C.CACHE_DIR / "validation.json").read_text())
    theirs, ours, px, py = d["theirs"], d["ours"], d["px"], d["py"]
    PX, PY = np.meshgrid(px, py)
    R = np.hypot(PX, PY)
    ext = (px[0], px[-1], py[-1], py[0])
    vmax = val["fig4"]["bar_vmax"]
    fig = plt.figure(figsize=(12.5, 8.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.62], hspace=0.38, wspace=0.28)
    pair = []
    for k, (img, title) in enumerate(((theirs, "Published Fig. 4, read back to A_V"),
                                      (ours, "Ours, same box and stretch"))):
        ax = fig.add_subplot(gs[0, k])
        im = ax.imshow(np.where(R < 1245, img, np.nan), extent=ext, origin="upper", cmap=ps.DUST,
                       norm=Normalize(0, vmax), interpolation="none")
        ps.map_axes(ax)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("x [pc]")
        if k == 0:
            ax.set_ylabel("y [pc]")
        pair.append(ax)
    _colorbar(fig, im, pair, "A_V [mag]", shrink=0.9)
    ax = fig.add_subplot(gs[0, 2])
    ratio = np.where((R < 1200) & (theirs > 0.15) & (theirs < 0.92 * vmax) & (ours < 0.92 * vmax), ours / theirs, np.nan)
    im = ax.imshow(ratio, extent=ext, origin="upper", cmap=ps.DIVERGING, norm=TwoSlopeNorm(1.0, 0.85, 1.15),
                   interpolation="none")
    ps.map_axes(ax)
    ax.set_title("Ours ÷ published", fontsize=10)
    ax.set_xlabel("x [pc]")
    _colorbar(fig, im, ax, "ratio", shrink=0.9)
    # the ratio against distance, and the samples test
    ax = fig.add_subplot(gs[1, :2])
    rb = val["fig4"]["ratio_by_distance"]
    rc = [0.5 * (r["r_lo"] + r["r_hi"]) for r in rb]
    ax.plot(rc, [r["ratio"] for r in rb], "o-", color=ps.BLUE, label="our map ÷ published Fig. 4 (median per 200-pc ring)")
    hold = val["holder"]["pairs"]
    ax.plot([h["r_pc"] for h in hold], [h["ratio_of_sums"] for h in hold], "s", color=ps.ORANGE, ms=7,
            label="interpolated mean ÷ mean of interpolated samples (12 posterior samples)")
    ax.errorbar([h["r_pc"] for h in hold], [h["ratio_median"] for h in hold],
                yerr=[[h["ratio_median"] - h["ratio_p16"] for h in hold], [h["ratio_p84"] - h["ratio_median"] for h in hold]],
                fmt="none", ecolor=ps.ORANGE, elinewidth=1, capsize=0, alpha=0.6)
    ax.axhline(1, color=ps.AXIS, lw=0.8)
    ax.set_xlabel("distance from the Sun [pc]")
    ax.set_ylabel("ratio")
    ax.set_ylim(0.98, 1.13)
    ax.legend(loc="upper left")
    ax.set_title("The gap is the log-interpolation of the mean vs. of the samples", fontsize=10)
    ax = fig.add_subplot(gs[1, 2])
    sel = (R < 1200) & (theirs < 0.92 * vmax) & (ours < 0.92 * vmax)
    hb = ax.hexbin(theirs[sel], ours[sel], gridsize=70, bins="log", cmap=ps.SEQ_BLUE, mincnt=1, linewidths=0)
    ax.plot([0, vmax], [0, vmax], color=ps.INK, lw=0.8)
    ax.set_xlabel("published A_V [mag]")
    ax.set_ylabel("ours [mag]")
    ax.set_title(f"r = {val['fig4']['pearson_r']:.3f}, slope {val['fig4']['slope']:.3f}", fontsize=10)
    ps.caption(fig, "The published figure (Edenhofer et al. 2024, arXiv source) was converted back to A_V with its own "
               f"colour bar (linear 0 to {vmax:.3f} mag; tick fit rms {1e3 * val['fig4']['tick_rms']:.1f} mmag). "
               "Its 12 posterior samples were fetched by byte range for four pairs of distance shells; at the shells' "
               "voxel centres their mean equals the released mean map to 3e-7.", y=0.02)
    ps.save(fig, C.FIGURE_DIR / "check_published_fig4")
    return "check_published_fig4"


def sky_extinction() -> str:
    """All-sky A_V to 250, 500, 750 and 1250 pc (cf. Edenhofer et al. 2024, Fig. 2)."""
    from deredden import Extinction
    ext = Extinction()
    nside = 128
    l, b = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.0), gridspec_kw={"hspace": 0.18, "wspace": 0.08})
    for ax, d in zip(axes.ravel(), (250, 500, 750, 1248)):
        av = C.ZGR23_TO_AV * ext.zgr23(l, b, np.full_like(l, float(d)))
        vmax = float(np.quantile(av, 0.999))
        im = mollweide(ax, av, cmap=ps.DUST, norm=Normalize(0, vmax))
        ax.set_title(f"A_V from the Sun to {d if d < 1248 else 1250} pc", fontsize=10)
        _colorbar(fig, im, ax, "A_V [mag]", orientation="horizontal", shrink=0.55, pad=0.03, extend="max")
    ps.caption(fig, "Galactic coordinates, l = 0 in the centre, longitude increasing to the left; each panel saturates at "
               "its 99.9% quantile. Integrated extinction from the map's voxels (each constant between its radial "
               "boundaries, plus the map's own integral inside 69 pc). Cf. Edenhofer et al. (2024), Fig. 2.", y=0.05)
    ps.save(fig, C.FIGURE_DIR / "sky_extinction")
    return "sky_extinction"


def bubble_wall() -> str:
    """The Local Bubble's wall: distance on the sky, a cut through the plane, and the O'Neill check."""
    import bubble
    m = bubble.shell_model()
    st = bubble.compare_oneill(m)
    ref = catalogs.oneill_shell()
    fig = plt.figure(figsize=(12.5, 9.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.05, 1.0], hspace=0.32, wspace=0.3)
    # (a) Mollweide map of the wall distance (from our map, NESTED -> RING for healpy plotting)
    ax = fig.add_subplot(gs[0, :2])
    im = mollweide(ax, m["d_peak"].to_numpy(), nest=True, cmap=ps.SEQ_BLUE, norm=Normalize(70, 400))
    ax.set_title("Distance from the Sun to the bubble's wall (first dust peak)", fontsize=10)
    _colorbar(fig, im, ax, "wall distance [pc] (dark: the open 'chimney' toward the north)", orientation="horizontal",
              shrink=0.5, pad=0.03, extend="max")
    # (b) mid-plane slice with the wall where it crosses |z| < 10 pc and young clusters
    g = grids.build_main()
    xs = g["xs"]
    ax = fig.add_subplot(gs[0, 2])
    sl = C.ZGR23_TO_AV * g["slice_z0"]
    im = ax.imshow(np.where(g["nan_z0"], np.nan, sl).T, extent=_extent(xs, xs), cmap=ps.DUST,
                   norm=LogNorm(1e-4, 3e-2), interpolation="none")
    near = (np.abs(m["z"]) < 10) & (np.hypot(m["x"], m["y"]) < 450)
    ax.plot(m.loc[near, "x"], m.loc[near, "y"], ".", ms=1.6, color=ps.BLUE, alpha=0.7, label="wall peak (|z| < 10 pc)")
    c = bubble.clusters_on_wall(m, rmax=450)
    yc = c[(c["age_myr"] < 20) & (np.abs(c["z"]) < 60)]
    ax.scatter(yc["x"], yc["y"], s=14, color=ps.ORANGE, marker="D", edgecolors=ps.SURFACE, linewidths=0.4,
               label="clusters < 20 Myr, |z| < 60 pc", zorder=3)
    ax.add_patch(Circle((0, 0), C.R_INNER, fill=False, ec=ps.MUTED, lw=0.6, ls="-"))
    _sun(ax, ms=8)
    ax.set_xlim(-450, 450)
    ax.set_ylim(-450, 450)
    ps.map_axes(ax)
    ax.set_aspect("equal")
    ax.set_xlabel("x [pc]")
    ax.set_ylabel("y [pc]")
    ax.set_title("Mid-plane density and the wall", fontsize=10)
    ax.legend(loc="lower left", fontsize=6.5)
    _colorbar(fig, im, ax, "A_V per pc [mag/pc] (log)", shrink=0.8)
    # (c-e) distributions: ours vs O'Neill et al. 2024
    for k, (col_o, col_r, label, rng, logx) in enumerate((("d_peak", "d", "wall distance [pc]", (60, 650), False),
                                                         ("thick", "thick", "wall thickness [pc]", (2, 400), True),
                                                         ("n_peak", "npeak", "peak density n_H [cm⁻³]", (0.01, 1000), True))):
        ax = fig.add_subplot(gs[1, k])
        bins = np.geomspace(*rng, 60) if logx else np.arange(rng[0] - 0.5, rng[1], 10.0)
        ax.hist(ref[col_r].dropna(), bins=bins, histtype="stepfilled", color=ps.CONTEXT, alpha=0.7,
                label="O'Neill et al. 2024 (released)")
        ax.hist(m[col_o].dropna(), bins=bins, histtype="step", color=ps.BLUE, lw=1.4, label="ours, same recipe")
        if logx:
            ax.set_xscale("log")
        ax.set_xlabel(label)
        ax.set_ylabel("sightlines")
        if k == 0:
            ax.legend(loc="upper right", fontsize=7)
            ax.set_title(f"{100 * st['frac_same_peak_within_1pc']:.1f}% of sightlines agree within 1 pc", fontsize=9.5)
        if k == 1:
            ax.set_title(f"median {np.nanmedian(m['thick']):.0f} pc (theirs {np.nanmedian(ref['thick']):.0f})", fontsize=9.5)
        if k == 2:
            ax.set_title(f"shell mass {st['mass_total_ours'] / 1e5:.2f}e5 M☉ (theirs 6.0e5)", fontsize=9.5)
    ps.caption(fig, "Wall = first peak along each of 786,432 sightlines (O'Neill et al. 2024's recipe: profile sampled every "
               "1 pc from 69 to 650 pc, Gaussian smoothing 7 pc, prominence > 2.04e-6 ZGR23/pc); 29 sightlines have no "
               "wall within 650 pc. n_H = 1653 cm⁻³ per ZGR23 mag/pc. Wall peak distances: 83% identical to O'Neill et "
               "al.'s, 15% differ by ±1 pc, symmetrically.", y=0.03)
    ps.save(fig, C.FIGURE_DIR / "bubble_wall")
    return "bubble_wall"


def bubble_history() -> str:
    """Zucker et al. (2022)'s expanding shell, refitted, against the dust wall today."""
    import history
    res = history.fit()
    s = np.load(C.CACHE_DIR / "bubble_history_samples.npy")
    m = history.ExpansionModel()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axes[0]
    t = np.linspace(-15.5, 0, 300)
    rng = np.random.default_rng(1)
    for th in s[rng.choice(len(s), 120, replace=False)]:
        tt = t - th[0]
        ax.plot(t, history.radius(tt, th[1], th[2]), color=ps.BLUE, alpha=0.05, lw=1)
    med = [np.median(s[:, i]) for i in range(4)]
    ax.plot(t, history.radius(t - med[0], med[1], med[2]), color=ps.BLUE, lw=2, label="shell radius (our fit)")
    r_birth = [np.linalg.norm(b["birth"] - m.centre(med[0])) for b in m.born]
    ax.scatter(m.t_birth + rng.uniform(-0.15, 0.15, len(m.t_birth)), r_birth, s=22, color=ps.ORANGE, marker="D",
               edgecolors=ps.SURFACE, linewidths=0.5, zorder=3, label="clusters at birth (distance from the centre)")
    ax.set_xlabel("time [Myr]  (0 = today)")
    ax.set_ylabel("radius [pc]")
    ax.set_xlim(-15.6, 0.3)
    ax.legend(loc="upper left")
    ax.set_title("Fitting the shell's growth to where young clusters were born", fontsize=10)
    txt = (f"first supernova {-res['t_exp_myr'][1]:.1f} Myr ago (Zucker+22: 14.4)\n"
           f"today R = {res['R_now_pc'][1]:.0f} ± {0.5 * (res['R_now_pc'][2] - res['R_now_pc'][0]):.0f} pc (165 ± 6), "
           f"v = {res['v_now_kms'][1]:.1f} km/s (6.7)")
    ax.text(0.98, 0.04, txt, transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5, color=ps.INK_SECONDARY)
    ax = axes[1]
    import bubble
    sh = bubble.shell_model()
    ok = sh["d_peak"].notna() & (sh["d_peak"] < 400)
    rr = np.linalg.norm(sh.loc[ok, ["x", "y", "z"]].to_numpy() - np.asarray(res["centre_pc"]), axis=1)
    ax.hist(rr, bins=np.linspace(0, 450, 91), color=ps.CONTEXT, histtype="stepfilled", density=True,
            label="dust wall today (each sightline's first peak)")
    ax.axvspan(res["R_now_pc"][0], res["R_now_pc"][2], color=ps.BLUE, alpha=0.25, lw=0)
    ax.axvline(res["R_now_pc"][1], color=ps.BLUE, lw=1.6, label="model shell today")
    ax.set_xlabel("distance from the explosion centre [pc]")
    ax.set_ylabel("fraction of sightlines per pc")
    ax.legend(loc="upper right", fontsize=7.5)
    ax.set_title("A sphere is only a sketch of the real wall", fontsize=10)
    ps.caption(fig, "Model and priors of Zucker et al. (2022): R(t) ∝ (Δt_SNe n0)^(-1/5) t^(3/5) (El-Badry et al. 2019), "
               "centre midway between UCL and LCC at the first explosion, 22 cluster birth positions from their galpy "
               "tracebacks (Harvard Dataverse). Centre today (heliocentric): "
               f"({res['centre_pc'][0]:.0f}, {res['centre_pc'][1]:.0f}, {res['centre_pc'][2]:.0f}) pc.", y=-0.04)
    ps.save(fig, C.FIGURE_DIR / "bubble_history")
    return "bubble_history"


def cloud_distances() -> str:
    """Cloud distances measured in the map vs Zucker et al. (2020)."""
    import wave
    cd = wave.cloud_distances()
    z = cd[(cd["set"] == "zucker2020") & (cd["d50"] < 1950)]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    ax = axes[0]
    same = np.abs(np.log(z["d_map"] / z["d50"])) < np.log(1.3)
    ax.errorbar(z["d50"][same], z["d_map"][same], xerr=[z["d50"][same] - z["d16"][same], z["d84"][same] - z["d50"][same]],
                fmt="o", ms=3.5, color=ps.BLUE, ecolor=ps.BLUE_RAMP[1], elinewidth=0.8, label="same cloud (within 30%)")
    ax.plot(z["d50"][~same], z["d_map"][~same], "D", ms=4, color=ps.ORANGE, mec=ps.SURFACE, mew=0.4,
            label="map's dominant cloud is another one")
    lim = [60, 2100]
    ax.plot(lim, lim, color=ps.INK, lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("Zucker et al. (2020) cloud distance [pc]")
    ax.set_ylabel("distance in the Edenhofer map [pc]")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_title(f"{len(z)} of the 326 cloud sightlines (those within 1.95 kpc)", fontsize=10)
    ax = axes[1]
    r = z["d_map"][same] / z["d50"][same]
    ax.scatter(z["d50"][same], r, s=10, color=ps.BLUE, alpha=0.7)
    edges = [60, 200, 400, 700, 1100, 2000]
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (z["d50"][same] >= lo) & (z["d50"][same] < hi)
        ax.plot([lo, hi], [np.median(r[s])] * 2, color=ps.INK, lw=2)
    ax.axhline(1, color=ps.AXIS, lw=0.8)
    ax.axhspan(0.95, 1.05, color=ps.NEUTRAL, lw=0, zorder=0)
    ax.set_xscale("log")
    ax.set_ylim(0.7, 1.35)
    ax.set_xlabel("Zucker et al. (2020) distance [pc]")
    ax.set_ylabel("map ÷ Zucker")
    ax.set_title("Map clouds sit 2-4% farther within 700 pc, more beyond", fontsize=10)
    ps.caption(fig, "Map distance = where half of the dominant cloud's extinction has accumulated along the sightline "
               "(pixels within 15′; main map to 1240 pc, 2 kpc map beyond). Grey band: Zucker et al.'s ~5% systematic. "
               "Black bars: medians per distance bin.", y=-0.03)
    ps.save(fig, C.FIGURE_DIR / "cloud_distances")
    return "cloud_distances"


def extinction_checks() -> str:
    """The dereddening tool against cluster photometry, SFD and the dustmaps package."""
    import pandas as pd
    val = json.loads((C.CACHE_DIR / "validation.json").read_text())["extinction"]
    c = pd.read_parquet(C.CACHE_DIR / "extinction_clusters.parquet")
    sfd = np.load(C.CACHE_DIR / "extinction_sfd.npz")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    ax = axes[0]
    ccm = 0.900                      # A_550(ZGR23 curve) / A_V(CCM R_V = 3.1) for the same reddening
    ax.errorbar(ccm * c["AV50"], c["A550_map"], xerr=[ccm * (c["AV50"] - c["AV16"]), ccm * (c["AV84"] - c["AV50"])],
                fmt="o", ms=2.5, color=ps.BLUE, ecolor=ps.BLUE_RAMP[0], elinewidth=0.6, alpha=0.8)
    ax.plot([0, 5], [0, 5], color=ps.INK, lw=0.8)
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 4)
    ax.set_xlabel("Hunt & Reffert (2023) cluster A_V × 0.900  [mag]")
    ax.set_ylabel("map: A at 550 nm to the cluster [mag]")
    big = c["AV50"] > 0.3
    ax.set_title(f"{len(c)} clusters, 100-1200 pc: median ratio {np.median(c['A550_map'][big] / (ccm * c['AV50'][big])):.3f}",
                 fontsize=10)
    ax = axes[1]
    ax.hexbin(sfd["sf11"], sfd["ebv_map"], gridsize=80, bins="log", cmap=ps.SEQ_BLUE, mincnt=1,
              extent=(0, 0.15, 0, 0.15), linewidths=0)
    ax.plot([0, 0.15], [0, 0.15], color=ps.INK, lw=0.8)
    fitl = val["sfd_high_latitude"]["fit_map_vs_sf11"]
    xx = np.linspace(0, 0.15, 10)
    ax.plot(xx, np.polyval(fitl, xx), color=ps.ORANGE, lw=1.4, label=f"fit: {fitl[0]:.3f} x + {fitl[1]:.4f}")
    ax.set_xlim(0, 0.15)
    ax.set_ylim(0, 0.15)
    ax.set_xlabel("SFD (1998) E(B-V) × 0.86 (Schlafly & Finkbeiner 2011)")
    ax.set_ylabel("map: E(B-V) to 1.25 kpc")
    ax.legend(loc="upper left")
    ax.set_title(f"|b| > 30°: {val['sfd_high_latitude']['n_pixels']} sightlines, r = {val['sfd_high_latitude']['pearson']:.2f}",
                 fontsize=10)
    ps.caption(fig, "Left: cluster A_V from Hunt & Reffert's isochrone network assume Cardelli et al. (1989) with R_V = 3.1; "
               "the ZGR23 curve behind the map is closest to R_V = 2.76, which assigns 0.900× the A_V to the same reddening. "
               "Right: E(B-V) of the map = (A440 - A550) of ZGR23's curve = 0.966 E; the +0.008 offset is the map's floor "
               "where there is almost no dust. dustmaps' integrated mode agrees with ours half a voxel further out "
               f"(99th-percentile ratio {val['dustmaps_integrated']['p99_ratio_shifted']:.6f}).", y=-0.03)
    ps.save(fig, C.FIGURE_DIR / "extinction_checks")
    return "extinction_checks"


def second_map() -> str:
    """Edenhofer et al. (2024) against Vergely, Lallement & Cox (2022)."""
    cmp_ = json.loads((C.CACHE_DIR / "compare.json").read_text())
    pr = np.load(C.CACHE_DIR / "compare_projections.npz")
    sl = np.load(C.CACHE_DIR / "compare_sightlines.npz")
    xs = pr["xs"]
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    mask = np.hypot(X, Y) > 1240
    fig = plt.figure(figsize=(12.5, 9.0))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 0.85], hspace=0.32, wspace=0.3)
    vmax = np.quantile(C.ZGR23_TO_AV * pr["ours"][~mask], 0.999)
    pair = []
    for k, (img, title) in enumerate(((pr["ours"], "Edenhofer et al. 2024"), (pr["theirs"], "Vergely et al. 2022"))):
        ax = fig.add_subplot(gs[0, k])
        pair.append(ax)
        im = ax.imshow(np.where(mask, np.nan, C.ZGR23_TO_AV * img).T, extent=_extent(xs, xs), cmap=ps.DUST,
                       norm=Normalize(0, vmax), interpolation="none")
        ps.map_axes(ax)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("x [pc]")
        if k == 0:
            ax.set_ylabel("y [pc]")
        _sun(ax, ms=7)
    _colorbar(fig, im, pair, "A_V through |z| ≤ 400 pc [mag]", orientation="horizontal", shrink=0.6, pad=0.16)
    ax = fig.add_subplot(gs[0, 2])
    ok = sl["inside"] & (sl["ours"] > 0.04)
    ax.hexbin(sl["ours"][ok], sl["theirs"][ok], gridsize=60, bins="log", cmap=ps.SEQ_BLUE, mincnt=1, xscale="log",
              yscale="log", extent=(np.log10(0.04), 1.3, np.log10(0.04), 1.3), linewidths=0)
    ax.plot([0.04, 20], [0.04, 20], color=ps.INK, lw=0.8)
    ax.set_xlabel("Edenhofer: A_V to (l, b, d) [mag]")
    ax.set_ylabel("Vergely [mag]")
    ax.set_title(f"integrated extinction: median ratio {cmp_['sightlines']['median_ratio_vergely_over_ours']:.3f}",
                 fontsize=10)
    ax = fig.add_subplot(gs[1, :])
    el = cmp_["elongation"]
    for name, color, marker in (("edenhofer", ps.BLUE, "o"), ("vergely", ps.ORANGE, "s")):
        rows = el[name]
        rc = [0.5 * (r["r_lo"] + r["r_hi"]) for r in rows]
        ratio = [r["elongation"] for r in rows]
        ax.plot(rc, ratio, marker=marker, color=color, label=f"{'Edenhofer et al. 2024' if name == 'edenhofer' else 'Vergely et al. 2022'}")
        for x_, rr, row in zip(rc, ratio, rows):
            if not np.isfinite(rr):
                ax.annotate("radial length\n> lag range", (x_, 3.05), ha="center", fontsize=6.5, color=color)
    ax.axhline(1, color=ps.AXIS, lw=0.8)
    ax.set_ylim(0.5, 3.4)
    ax.set_xlabel("distance from the Sun [pc]")
    ax.set_ylabel("correlation length along ÷ across the sightline")
    ax.legend(loc="upper left")
    ax.set_title("Fingers of god: how much each map stretches clouds along the line of sight", fontsize=10)
    ps.caption(fig, "Vergely et al.'s 10 pc-resolution cube (A0 at 550 nm; ÷ 2.658 for ZGR23 units). Their grid aligns with "
               f"ours to ({cmp_['offset']['shift_x_pc']:.1f}, {cmp_['offset']['shift_y_pc']:.1f}) pc. Elongation: "
               "correlation of ln(density) after removing 30-pc-smoothed trends, both maps on 5 pc cells, points with "
               "|z| < 150 pc.", y=0.03)
    ps.save(fig, C.FIGURE_DIR / "second_map")
    return "second_map"


def uncertainty() -> str:
    """Relative posterior uncertainty of single voxels in the mid-plane."""
    g = grids.build_main()
    xs = g["xs"]
    rel = np.where(g["nan_z0"] | (g["slice_z0"] <= 0), np.nan, g["std_z0"] / np.maximum(g["slice_z0"], 1e-12))
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    im = ax.imshow(rel.T, extent=_extent(xs, xs), cmap=ps.SEQ_BLUE, norm=Normalize(0.2, 1.2), interpolation="none")
    ps.map_axes(ax)
    ax.set_aspect("equal")
    _sun(ax, ms=8)
    ax.set_xlabel("x [pc]")
    ax.set_ylabel("y [pc]")
    ax.set_title("Single-voxel uncertainty in the mid-plane (std ÷ mean)", fontsize=10)
    _colorbar(fig, im, ax, "relative posterior standard deviation", shrink=0.85)
    med = np.nanmedian(rel)
    ps.caption(fig, f"Median {med:.2f}. Sums over many voxels are much better determined (the paper quotes ~10% for its "
               "projections), although neighbouring voxels are positively correlated in the posterior samples "
               "(+0.3 to +0.5 between adjacent shells or pixels, measured on the released samples).", y=0.01)
    ps.save(fig, C.FIGURE_DIR / "uncertainty_midplane")
    return "uncertainty_midplane"


def _wave_fit(tag):
    path = C.CACHE_DIR / f"wave_fit_{tag}.json"
    return json.loads(path.read_text()) if path.exists() else None


def wave_faceon() -> str:
    """The Radcliffe Wave from above: the 2 kpc map, the clouds, the fitted model, young clusters."""
    import wave
    w = grids.build_wide()
    xs = w["xs"]
    img = to_av(w["proj_xy"])
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    vmax = np.quantile(img[np.hypot(X, Y) < 1980], 0.995)
    fit_map = _wave_fit(wave.MAIN_FIT)
    cd = wave.cloud_distances()
    z = cd[cd["set"] == "zucker2020"]
    sets = wave.cloud_sets()
    fig, ax = plt.subplots(figsize=(9.0, 9.0))
    ax.imshow(np.where(np.hypot(X, Y) < 1985, img, np.nan).T, extent=_extent(xs, xs), cmap=ps.DUST,
              norm=Normalize(0, vmax), interpolation="none")
    ps.map_axes(ax)
    ax.set_aspect("equal")
    for R in (500, 1000, 1500, 2000):
        ax.add_patch(Circle((0, 0), R, fill=False, ec=ps.AXIS, lw=0.5))
    hr = catalogs.hunt_reffert()
    yc = hr[hr["reliable"] & (hr["age_myr"] < 30) & (hr["dist50"] < 2000)]
    ax.scatter(yc["x"], yc["y"], s=8, color=ps.BLUE, edgecolors="none", alpha=0.7, label="clusters < 30 Myr (Hunt & Reffert 2023)")
    if fit_map is not None:
        mem = np.asarray(fit_map["members"], dtype=bool)
        idx = sets["map_index"]
        P = z.loc[idx, ["xm", "ym"]].to_numpy()
        ax.scatter(P[~mem, 0], P[~mem, 1], s=10, marker="D", color=ps.CONTEXT, edgecolors="none",
                   label="cloud sightlines, distance from the map: off the wave")
        ax.scatter(P[mem, 0], P[mem, 1], s=16, marker="D", color=ps.AQUA, edgecolors=ps.SURFACE, linewidths=0.4,
                   label="... on the wave (mixture-model members)")
        x, y, zz, _ = wave.kon_curve(fit_map["median"], 600)
        ax.plot(x, y, color=ps.ORANGE, lw=2.2, label="our fit (cloud distances from the map)")
    x, y, zz, _ = wave.kon_curve(wave.KON_PUBLISHED["clouds"], 600)
    ax.plot(x, y, color=ps.INK, lw=1.0, label="Konietzka et al. 2024 (clouds fit)")
    _sun(ax)
    ax.set_xlim(-2000, 2000)
    ax.set_ylim(-2000, 2000)
    ax.set_xlabel("x [pc]  (toward the Galactic Centre →)")
    ax.set_ylabel("y [pc]")
    ax.set_title("The Radcliffe Wave seen from above (2 kpc map, |z| < 500 pc)", fontsize=10.5)
    ax.legend(loc="lower right", fontsize=7, framealpha=0.85, frameon=True, edgecolor=ps.AXIS)
    for name in ("CMa OB1", "Orion", "Perseus", "Taurus", "Cepheus", "North America", "Cygnus X"):
        lab = catalogs.cloud_labels().set_index("complex")
        if name in lab.index:
            ax.annotate(name, (lab.loc[name, "x"], lab.loc[name, "y"]), xytext=(-8, 8), textcoords="offset points",
                        ha="right", fontsize=7, color=ps.INK_SECONDARY, path_effects=HALO)
    ps.caption(fig, "Background: the 2 kpc validation map of Edenhofer et al. (2024), A_V integrated over |z| < 500 pc "
               "(linear grey, saturating at the 99.5% quantile; the authors advise caution with its fine detail beyond "
               "1.25 kpc). Rings every 500 pc.", y=0.04)
    ps.save(fig, C.FIGURE_DIR / "wave_faceon")
    return "wave_faceon"


def wave_curtain() -> str:
    """The dust along the wave, unrolled: it rises and falls with no fitting involved."""
    import wave
    fit_map = _wave_fit(wave.MAIN_FIT)
    th = fit_map["median"] if fit_map is not None else wave.KON_PUBLISHED["clouds"]
    path = C.CACHE_DIR / "wave_curtain.npz"
    key = json.dumps(th, sort_keys=True)
    if path.exists() and str(np.load(path)["key"]) == key:
        cu = dict(np.load(path))
    else:
        cu = wave.curtain(th)
        np.savez_compressed(path, key=key, **{k: v for k, v in cu.items() if isinstance(v, np.ndarray)})
    col = to_av(cu["column"])
    fig, axes = plt.subplots(2, 1, figsize=(13, 7.6), gridspec_kw={"hspace": 0.45, "height_ratios": [1, 1]})
    for k, ax in enumerate(axes):
        vmax = np.quantile(col, 0.995)
        ax.imshow(col.T, origin="lower", extent=_extent(cu["s"], cu["z"]), cmap=ps.DUST,
                  norm=Normalize(0, vmax), interpolation="none", aspect="equal")
        ps.map_axes(ax)
        x, y, zz, s = wave.kon_curve(th, 800)
        ax.plot(s, zz, color=ps.ORANGE, lw=2.0, label="our wave model")
        xk, yk, zk, sk = wave.kon_curve(wave.KON_PUBLISHED["clouds"], 800)
        sp = wave.project_to_baseline(th, np.c_[xk, yk, zk])
        ax.plot(sp[:, 0], sp[:, 2], color=ps.INK, lw=0.9, label="Konietzka et al. 2024")
        ax.axhline(0, color=ps.AXIS, lw=0.6)
        ax.set_xlim(cu["s"][0], cu["s"][-1])
        ax.set_ylabel("z [pc]")
        ax.set_xlabel("distance along the wave's baseline [pc]  (Canis Major → Cygnus)")
        if k == 0:
            cd = wave.cloud_distances()
            z = cd[cd["set"] == "zucker2020"]
            P = wave.project_to_baseline(th, z[["xm", "ym", "zm"]].dropna().to_numpy())
            near = np.abs(P[:, 1]) < 150
            ax.scatter(P[near, 0], P[near, 2], s=14, marker="D", color=ps.AQUA, edgecolors=ps.SURFACE, linewidths=0.4,
                       label="cloud sightlines within 150 pc (map distances)", zorder=3)
            ax.set_title("Dust within 150 pc of the wave's baseline, unrolled and seen from the side (true aspect ratio)",
                         fontsize=10)
        else:
            hr = catalogs.hunt_reffert()
            yc = hr[hr["reliable"] & (hr["age_myr"] < 30)]
            P = wave.project_to_baseline(th, yc[["x", "y", "z"]].to_numpy())
            near = (np.abs(P[:, 1]) < 150) & (P[:, 0] > cu["s"][0]) & (P[:, 0] < cu["s"][-1])
            ax.scatter(P[near, 0], P[near, 2], s=12, color=ps.BLUE, edgecolors=ps.SURFACE, linewidths=0.3,
                       label="star clusters < 30 Myr within 150 pc", zorder=3)
            ax.set_title("The same with young star clusters", fontsize=10)
        ax.legend(loc="upper left", fontsize=7, ncols=3, framealpha=0.85, frameon=True, edgecolor=ps.AXIS)
    ps.caption(fig, "Each column: A_V through a 300-pc-wide band centred on the wave's in-plane baseline (2 kpc map of "
               "Edenhofer et al. 2024, sampled every 5 pc). No fit enters the image; the curves are models.", y=0.03)
    ps.save(fig, C.FIGURE_DIR / "wave_curtain")
    return "wave_curtain"


def render_modes() -> str:
    """One view of the viewer in its three render modes and two stretches: rendering is a choice."""
    import subprocess
    from PIL import Image
    out = C.CACHE_DIR / "render_modes"
    page = C.FIGURE_DIR / "dust_flythrough.html"
    jobs = [(2, 0.0, "Solid (the movie): opacity and fade are choices"),
            (0, 0.0, "Glow: brightness = dust column (log stretch)"),
            (1, 0.0, "Starlight: white light dimmed by 10^(-0.4 A)"),
            (0, 1.5, "Glow, stretched 4.5x: faint dust looks like structure")]
    frames = []
    for mode, expo, title in jobs:
        sub = out / f"m{mode}_e{expo:+.1f}"
        subprocess.run(["node", str(C.PROJECT_DIR / "tools" / "render_frames.js"), str(page), str(sub), "--mode", str(mode),
                        "--times", "52", "--size", "960x540", "--exposure", str(expo)], check=True,
                       capture_output=True)
        frames.append((Image.open(sub / f"t52_m{mode}.png").convert("RGB"), title))
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 7.6), gridspec_kw={"hspace": 0.12, "wspace": 0.03})
    for ax, (img, title) in zip(axes.ravel(), frames):
        ax.imshow(np.asarray(img), origin="upper")
        ax.set_title(title, fontsize=9.5)
        ax.axis("off")
    ps.caption(fig, "The same camera, data and overlays (looking down on the Local Bubble from 900 pc, 50 s into the "
               "tour); only the mapping from dust density to light changes. Glow and Starlight are fixed physical "
               "mappings; the stretch and Solid's opacity are choices, stated in the viewer.", y=0.06)
    ps.save(fig, C.FIGURE_DIR / "render_modes")
    return "render_modes"


ALL = [projections, fig4_check, sky_extinction, bubble_wall, bubble_history, cloud_distances, extinction_checks,
       second_map, uncertainty, wave_faceon, wave_curtain, render_modes]


if __name__ == "__main__":
    import sys
    names = sys.argv[1:]
    for f in ALL:
        if not names or f.__name__ in names:
            print("figure", f.__name__, flush=True)
            if f is projections:
                projections(labelled=False)
                projections(labelled=True)
            else:
                f()
