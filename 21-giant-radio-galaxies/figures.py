"""Every figure of the project. Each function writes PNG + SVG into figures/."""

from __future__ import annotations

import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS
from matplotlib.patches import Ellipse, Rectangle
from PIL import Image

import catalogs
import config
import cosmology
import plotstyle as ps
import radiomaps as rm
from config import FIGURE_DIR, IMAGE_DIR

warnings.filterwarnings("ignore", category=RuntimeWarning)
ps.use_style()

DARK = "#0b0b0b"


# --------------------------------------------------------------------------
# Image helpers
# --------------------------------------------------------------------------

def legacy_image(g: config.Giant, *, zoom: bool = False) -> tuple[np.ndarray, WCS]:
    """The optical JPEG (Legacy Surveys, or Pan-STARRS1 near the pole) as a
    float RGB array (origin lower) and its WCS.

    Both services cut a TAN projection centred on the requested position,
    north up, east left, spanning the requested field; we rebuild that WCS
    from the request itself (config.legacy_url / config.ps1_url)."""
    path = IMAGE_DIR / f"{g.key}_{'host' if zoom else 'legacy'}.jpg"
    rgb = np.asarray(Image.open(path).convert("RGB"), float)[::-1] / 255.0
    size_arcsec = config.HOST_ZOOM_ARCSEC if zoom else g.size_arcmin * 60
    n = rgb.shape[0]
    pixscale = size_arcsec / n
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [g.ra, g.dec]
    w.wcs.crpix = [(rgb.shape[1] + 1) / 2, (n + 1) / 2]
    w.wcs.cdelt = [-pixscale / 3600, pixscale / 3600]
    return rgb, w


def radio_layer(data: np.ndarray, sigma: float, *, floor: float = 2.0, top: float | None = None) -> np.ndarray:
    """Radio intensity mapped to 0..1 with an arcsinh stretch: faint lobes
    visible, bright hotspots not saturated."""
    top = top if top is not None else np.nanpercentile(data[data > 5 * sigma], 99.5) if np.any(data > 5 * sigma) else 50 * sigma
    s = np.arcsinh(np.clip(data - floor * sigma, 0, None) / (2 * sigma)) / np.arcsinh(max(top, 10 * sigma) / (2 * sigma))
    return np.nan_to_num(np.clip(s, 0, 1))


def overlay(optical: np.ndarray, radio01: np.ndarray, color=ps.C_RADIO, dim: float = 0.8) -> np.ndarray:
    """Screen-blend a single-hue radio glow over a (dimmed) optical image."""
    c = np.array(matplotlib.colors.to_rgb(color))
    glow = radio01[..., None] * c[None, None, :]
    base = optical * dim
    return 1 - (1 - base) * (1 - glow)


def _crop_box(w: WCS, shape, ra, dec, half_arcmin: float):
    x, y = w.world_to_pixel_values(ra, dec)
    pix = abs(w.pixel_scale_matrix[1, 1]) * 60
    h = half_arcmin / pix
    x0, x1 = int(max(0, x - h)), int(min(shape[1], x + h))
    y0, y1 = int(max(0, y - h)), int(min(shape[0], y + h))
    return x0, x1, y0, y1


def _scale_bar(ax, x0, y0, length_px, label, color="white"):
    ax.plot([x0, x0 + length_px], [y0, y0], color=color, lw=1.6, solid_capstyle="butt")
    ax.text(x0 + length_px / 2, y0 + 0.012 * ax.get_ylim()[1], label, color=color, ha="center", va="bottom",
            fontsize=7)


def _midpoint(ends):
    (r1, d1), (r2, d2) = ends
    v1 = np.array([np.cos(np.radians(d1)) * np.cos(np.radians(r1)), np.cos(np.radians(d1)) * np.sin(np.radians(r1)), np.sin(np.radians(d1))])
    v2 = np.array([np.cos(np.radians(d2)) * np.cos(np.radians(r2)), np.cos(np.radians(d2)) * np.sin(np.radians(r2)), np.sin(np.radians(d2))])
    v = (v1 + v2) / np.linalg.norm(v1 + v2)
    return float(np.degrees(np.arctan2(v[1], v[0])) % 360), float(np.degrees(np.arcsin(v[2])))


def display_map(g: config.Giant) -> tuple[rm.RadioImage, float]:
    """The radio image we show for a giant: 20" for compact lobes; for faint,
    diffuse ones the 20" map smoothed to 40" (60" for the faintest, measured at
    90-120"), a compromise between detail and showing the lobes at all."""
    base = rm.load_lotss(g.key, 20)
    if g.measure_beam >= 90:
        img = rm.smooth_to(base, 60.0)
    elif g.measure_beam >= 60:
        img = rm.smooth_to(base, 40.0)
    else:
        img = base
    return img, rm.noise(img)


def display_extent(g: config.Giant, art: dict) -> rm.Extent:
    """The outline that decides what is shown bright: the inclusive 60" one,
    or the giant's own measured outline where that is smoother still."""
    a = art[g.key]
    return a["adopted"] if g.measure_beam > 60 else a["ext"]["60"]


# --------------------------------------------------------------------------
# 00: the gallery (the keeper)
# --------------------------------------------------------------------------

def giant_mask_on(g: config.Giant, ext: rm.Extent, wcs_out: WCS, shape) -> np.ndarray:
    """The giant's measured 3-sigma footprint, regridded and grown by one beam
    so its faint edges are not clipped."""
    from scipy import ndimage
    from reproject import reproject_interp

    base = rm.load_lotss(g.key, 20)
    grow = int(round(ext.beam_arcsec / base.pixel_arcsec))
    m = ndimage.binary_dilation(ext.mask, iterations=max(grow, 1)).astype(float)
    out, _ = reproject_interp((m, base.wcs), wcs_out, shape_out=shape, order="bilinear")
    return np.nan_to_num(out) > 0.5


def figure_gallery(measured: list[dict], art: dict, group: str = "gallery") -> list[str]:
    rows = {r["key"]: r for r in measured}
    giants = config.giants_in(group)
    giants.sort(key=lambda g: -rows[g.key]["adopted"]["lp_mpc"])
    ncol = 5
    nrow = int(np.ceil(len(giants) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(17, 3.65 * nrow + 0.9), facecolor=DARK)
    for ax in np.ravel(axes):
        ax.set_facecolor(DARK)
        ax.axis("off")
    for ax, g in zip(np.ravel(axes), giants):
        row = rows[g.key]
        opt, wopt = legacy_image(g)
        img, sig = display_map(g)
        rad = rm.reproject_to(img, wopt, opt.shape[:2])
        layer = radio_layer(rad, sig)
        mask = giant_mask_on(g, display_extent(g, art), wopt, opt.shape[:2])
        # the giant at full strength; every other radio source faded
        rgb = overlay(opt, layer * np.where(mask, 1.0, 0.35))
        ends = row["adopted"]["ends"]
        cra, cdec = _midpoint(ends)
        half = 0.5 * row["adopted"]["las_arcmin"] * 1.22 + 0.8
        x0, x1, y0, y1 = _crop_box(wopt, opt.shape, cra, cdec, half)
        ax.imshow(rgb[y0:y1, x0:x1], origin="lower", interpolation="lanczos")
        hx, hy = wopt.world_to_pixel_values(g.ra, g.dec)
        ax.scatter([hx - x0], [hy - y0], s=60, facecolors="none", edgecolors="white", linewidths=0.8)
        # scale bar of 1 Mpc and the Milky Way to scale
        pix_arcmin = abs(wopt.wcs.cdelt[1]) * 60
        mpc_px = cosmology.angle_arcmin(1.0, g.z, config.PLANCK18) / pix_arcmin
        w_px, h_px = x1 - x0, y1 - y0
        bar = mpc_px if mpc_px < 0.6 * w_px else mpc_px / 2
        bar_label = "1 Mpc" if mpc_px < 0.6 * w_px else "0.5 Mpc"
        bx, by = 0.06 * w_px, 0.07 * h_px
        ax.plot([bx, bx + bar], [by, by], color="white", lw=1.5, solid_capstyle="butt")
        ax.text(bx + bar / 2, by + 0.025 * h_px, bar_label, color="white", ha="center", va="bottom", fontsize=7.5)
        mw_px = config.MILKY_WAY_DIAMETER_KPC / 1e3 * mpc_px
        mx = bx + bar + 0.08 * w_px
        dot = max(mw_px, 0.006 * w_px)  # never smaller than a visible dot
        ax.add_patch(Ellipse((mx, by), dot, max(mw_px * 0.3, 0.006 * w_px), color="white", lw=0))
        ax.text(mx + dot / 2 + 0.02 * w_px, by, "Milky Way", color="white", fontsize=6.5, va="center")
        ax.set_xlim(0, w_px)
        ax.set_ylim(0, h_px)
        zt = f"z = {g.z:.3f}" if g.z < 0.1 else f"z = {g.z:.2f}"
        ax.set_title(f"{g.name}", color="white", fontsize=10.5, loc="left", pad=4)
        ax.text(0.0, 1.0, "", transform=ax.transAxes)
        ax.text(0.02, 0.97, f"{row['adopted']['lp_mpc']:.1f} Mpc   ({row['adopted']['las_arcmin']:.0f}')\n{zt}",
                transform=ax.transAxes, color="white", fontsize=7.5, va="top", ha="left", linespacing=1.3)
    if group == "dr3":
        title = "Beyond DR2: giants that only LoTSS DR3 shows (144 MHz, orange) over the optical sky"
        source = ("Radio: LoTSS DR3 (Shimwell et al. 2026), 20\" resolution, smoothed to 40\" for the faint, diffuse "
                  "giants and 60\" for the faintest (NGC 2300, J1526+2620); NGC 2300 is a mosaic of four cutouts. Optical: "
                  "DESI Legacy Imaging Surveys DR10; Pan-STARRS1 for NGC 2300, beyond the Legacy Surveys' reach near the "
                  "celestial pole.")
        out = FIGURE_DIR / "10_dr3_gallery"
    else:
        title = "Giant radio galaxies: LOFAR at 144 MHz (orange) over Legacy Surveys optical images"
        source = ("Radio: LoTSS DR2 (Shimwell et al. 2022), 20\" resolution, smoothed to 40\" for the faint, diffuse "
                  "giants. Optical: DESI Legacy Imaging Surveys DR10.")
        out = FIGURE_DIR / "00_gallery"
    fig.suptitle(title, color="white", x=0.01, ha="left", fontsize=13, y=0.995)
    fig.text(0.01, 0.005, "Each panel shows one giant at its own scale; the white bar is 1 (or 0.5) Mpc and the white mark "
             "beside it the Milky Way's 30-kpc disc at the same distance (a dot where it would be smaller). Circles mark "
             "host galaxies. The giant's own emission is bright; unrelated radio sources are faded. Lengths are our own "
             "3-sigma measurements (Planck 2018 cosmology). " + source,
             color="#bdbcb5", fontsize=7.5, ha="left", va="bottom", wrap=True)
    fig.subplots_adjust(left=0.005, right=0.995, top=0.94, bottom=0.05, wspace=0.03, hspace=0.12)
    return ps.save(fig, out)


# --------------------------------------------------------------------------
# 01: the size lineup
# --------------------------------------------------------------------------

LINEUP_KPC_PER_PX = 6.0


def lineup_strip(g: config.Giant, row: dict, ext: rm.Extent, max_half_height_mpc: float = 0.9
                 ) -> tuple[np.ndarray, float, float]:
    """The giant's own radio emission (its measured footprint only) resampled
    to a fixed physical scale, turned so that the axis through its two lobe
    seeds is horizontal. Returns the image, its length and its height (Mpc),
    the height trimmed to the emission."""
    img, sig = display_map(g)
    # centred between the two lobe seeds, so both lobes sit on the centre line
    cra, cdec = _midpoint(g.seeds)
    pa = rm.position_angle(*g.seeds[0], *g.seeds[1])
    kpc_per_arcsec = cosmology.kpc_per_arcmin(g.z, config.PLANCK18) / 60
    pix_arcsec = LINEUP_KPC_PER_PX / kpc_per_arcsec
    length_mpc = row["adopted"]["lp_mpc"] + 0.8
    nx = int(np.ceil(length_mpc * 1e3 / LINEUP_KPC_PER_PX))
    ny = int(np.ceil(2 * max_half_height_mpc * 1e3 / LINEUP_KPC_PER_PX))
    w = rm.tangent_wcs(cra, cdec, pix_arcsec, max(nx, ny), rotation_deg=pa)
    w.wcs.crpix = [(nx + 1) / 2, (ny + 1) / 2]
    data = rm.reproject_to(img, w, (ny, nx))
    mask = giant_mask_on(g, ext, w, (ny, nx))
    layer = radio_layer(data, sig) * mask
    rows_on = np.where(layer.max(axis=1) > 0.02)[0]
    if rows_on.size:
        c = (ny - 1) / 2
        half = max(abs(rows_on.min() - c), abs(rows_on.max() - c)) + 2
        lo, hi = int(max(0, c - half)), int(min(ny, c + half + 1))
        layer = layer[lo:hi]
    return layer, length_mpc, layer.shape[0] * LINEUP_KPC_PER_PX / 1e3


def figure_lineup(measured: list[dict], art: dict) -> list[str]:
    rows = {r["key"]: r for r in measured}
    keys = [g.key for g in config.GIANTS if "gallery" in g.groups or "dr3" in g.groups]
    keys.sort(key=lambda k: -rows[k]["adopted"]["lp_mpc"])
    strips = [(k, *lineup_strip(config.giant_by_key(k), rows[k], display_extent(config.giant_by_key(k), art)))
              for k in keys]
    max_len = max(s[2] for s in strips)
    gap = 0.12
    total = sum(s[3] for s in strips) + gap * len(strips)
    fig = plt.figure(figsize=(14, 1.6 + 0.92 * (total + 1.0) * 14 * 0.66 / (max_len + 0.4)), facecolor=DARK)
    ax = fig.add_axes([0.2, 0.08, 0.66, 0.86])
    ax.set_facecolor(DARK)
    c = np.array(matplotlib.colors.to_rgb(ps.C_RADIO))
    y_top = total
    for k, s, length, height in strips:
        g = config.giant_by_key(k)
        y0 = y_top - height
        x_extent = s.shape[1] * LINEUP_KPC_PER_PX / 1e3
        rgba = np.zeros(s.shape + (4,))
        rgba[..., :3] = c
        rgba[..., 3] = s
        ax.imshow(rgba, origin="lower", extent=[-x_extent / 2, x_extent / 2, y0, y0 + height],
                  interpolation="bilinear")
        lp = rows[k]["adopted"]["lp_mpc"]
        label = g.name + ("  (DR3)" if g.survey == "dr3" else "")
        ax.text(-max_len / 2 - 0.15, y0 + height / 2, label, color="white", ha="right", va="center", fontsize=8.5)
        ax.text(max_len / 2 + 0.15, y0 + height / 2, f"{lp:.2f} Mpc   z = {g.z:.2f}", color="#bdbcb5",
                ha="left", va="center", fontsize=7.5)
        y_top = y0 - gap
    # the Local Group to scale
    yl = -0.55
    sep = config.ANDROMEDA_DISTANCE_KPC / 1e3
    ax.add_patch(Ellipse((-sep / 2, yl), config.MILKY_WAY_DIAMETER_KPC / 1e3, 0.008, color="white"))
    ax.add_patch(Ellipse((sep / 2, yl), config.ANDROMEDA_DIAMETER_KPC / 1e3, 0.012, color="white"))
    ax.plot([-sep / 2, sep / 2], [yl - 0.12, yl - 0.12], color="#bdbcb5", lw=0.8)
    ax.text(0, yl - 0.16, f"{sep * 1e3:.0f} kpc", color="#bdbcb5", ha="center", va="top", fontsize=7)
    ax.text(-sep / 2, yl + 0.06, "Milky Way", color="white", ha="center", va="bottom", fontsize=7.5)
    ax.text(sep / 2, yl + 0.06, "Andromeda", color="white", ha="center", va="bottom", fontsize=7.5)
    ax.text(-max_len / 2 - 0.15, yl, "The Local Group, to scale", color="white", ha="right", va="center", fontsize=8.5)
    ax.set_xlim(-max_len / 2 - 0.1, max_len / 2 + 0.1)
    ax.set_ylim(yl - 0.4, total)
    ax.set_aspect("equal")
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_color("#52514e")
    ax.tick_params(axis="x", colors="#bdbcb5")
    ax.set_xlabel("projected length (Mpc)", color="#bdbcb5")
    ax.grid(False)
    fig.suptitle("Giants lined up by size: each radio galaxy drawn at the same physical scale, axis horizontal",
                 color="white", x=0.01, ha="left", fontsize=12.5)
    fig.text(0.01, 0.01, f"LoTSS DR2 (DR3 where marked) radio emission resampled to {LINEUP_KPC_PER_PX:.0f} kpc per pixel along each "
             "giant's axis (Planck 2018 cosmology); only each giant's own emission is shown. Projected lengths are lower "
             "limits to the true, 3D lengths. Claimed giants whose ends may be unrelated sources are in figure 03.",
             color="#bdbcb5", fontsize=7.5, ha="left", va="bottom")
    return ps.save(fig, FIGURE_DIR / "01_lineup")


# --------------------------------------------------------------------------
# 02: Porphyrion in full
# --------------------------------------------------------------------------

QUASAR_NEAR_PORPHYRION = (232.387625, 60.264583, 0.79836)  # SDSS J152933.03+601552.5 (BOSS redshift)


def _radio_panel(ax, img, sig, *, cmap="gray_r", vmax_sigma=60, contour_mask=None):
    shown = np.arcsinh(np.clip(img.data, -3 * sig, None) / (1.5 * sig))
    ax.imshow(shown, origin="lower", cmap=cmap, vmin=-1, vmax=np.arcsinh(vmax_sigma / 1.5), interpolation="nearest")
    if contour_mask is not None:
        ax.contour(contour_mask.astype(float), levels=[0.5], colors=[ps.ORANGE], linewidths=0.9)


def figure_porphyrion(measured: list[dict], art: dict, lres: dict) -> list[str]:
    g = config.giant_by_key("porphyrion")
    row = {r["key"]: r for r in measured}["porphyrion"]
    ext = art["porphyrion"]["adopted"]
    img = rm.load_lotss("porphyrion", 20)
    sig = rm.noise(img)
    fig = plt.figure(figsize=(16, 6.2))
    ax1 = fig.add_axes([0.03, 0.1, 0.33, 0.8])
    _radio_panel(ax1, img, sig, contour_mask=ext.mask)
    (r1, d1), (r2, d2) = row["adopted"]["ends"]
    x, y = img.world_to_pixel(np.array([r1, r2]), np.array([d1, d2]))
    ax1.plot(x, y, color=ps.BLUE, lw=1.2)
    hx, hy = img.world_to_pixel(g.ra, g.dec)
    ax1.scatter(hx, hy, s=80, facecolors="none", edgecolors=ps.BLUE, linewidths=1.2)
    ps.annotate(ax1, "host galaxy", (hx[0], hy[0]), (hx[0] + 60, hy[0] - 25), color=ps.BLUE)
    nx_, ny_ = img.world_to_pixel(*g.seeds[0])
    sx_, sy_ = img.world_to_pixel(*g.seeds[1])
    ps.annotate(ax1, "northern lobe", (nx_[0], ny_[0]), (nx_[0] + 40, ny_[0] + 30))
    ps.annotate(ax1, "southern hotspot", (sx_[0], sy_[0]), (sx_[0] - 120, sy_[0] - 30))
    bx, by = img.world_to_pixel(232.52, 60.405)
    ps.annotate(ax1, "unrelated bright source", (bx[0], by[0]), (bx[0] - 40, by[0] + 15))
    mpc_px = cosmology.angle_arcmin(1.0, g.z, config.PLANCK18) * 60 / img.pixel_arcsec
    ax1.plot([20, 20 + mpc_px], [20, 20], color=ps.INK, lw=1.5)
    ax1.text(20 + mpc_px / 2, 27, "1 Mpc", ha="center", fontsize=8)
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_title("LoTSS DR2, 144 MHz, 20\" resolution", fontsize=10)
    for sp in ax1.spines.values():
        sp.set_visible(True)
    # host close-up
    ax2 = fig.add_axes([0.375, 0.1, 0.24, 0.8])
    opt, wopt = legacy_image(g, zoom=True)
    ax2.imshow(np.clip(opt * 1.15, 0, 1), origin="lower", interpolation="nearest")
    im6 = rm.load_lotss("porphyrion", 6)
    s6 = rm.noise(im6)
    r6 = rm.reproject_to(im6, wopt, opt.shape[:2])
    ax2.contour(r6, levels=s6 * np.array([4, 8, 16]), colors=[ps.ORANGE], linewidths=0.8)
    for (ra, dec, lab, col, dxy) in ((g.ra, g.dec, f"host, z = {g.z:.3f}\n(Keck spectrum)", ps.BLUE, (20, -55)),
                                     (QUASAR_NEAR_PORPHYRION[0], QUASAR_NEAR_PORPHYRION[1],
                                      f"quasar, z = {QUASAR_NEAR_PORPHYRION[2]:.3f}\n(foreground, unrelated)", "white", (15, 35))):
        px, py = wopt.world_to_pixel_values(ra, dec)
        ax2.scatter(px, py, s=260, facecolors="none", edgecolors=col, linewidths=1.2)
        ax2.annotate(lab, (px, py), (px + dxy[0], py + dxy[1]), color=col, fontsize=8,
                     arrowprops={"arrowstyle": "-", "color": col, "lw": 0.8})
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_title("The host: Legacy Surveys DR10, 80\" square;\nLoTSS 6\" contours (4, 8, 16 sigma)", fontsize=10)
    # how long is it?
    ax3 = fig.add_axes([0.66, 0.14, 0.32, 0.72])
    dep = lres["porphyrion_deprojection"]
    las = row["adopted"]["las_arcmin"]
    lp_qso = float(cosmology.projected_length_mpc(las, QUASAR_NEAR_PORPHYRION[2], config.PLANCK18))
    entries = [
        ("Mostert et al. 2024 (citizen-science\nassociation, photometric z = 0.916)", 8.39, None, ps.C_PUBLISHED),
        ("Hardcastle et al. 2023 (host and z uncertain)", 6.9, None, ps.C_PUBLISHED),
        ("Oei et al. 2024, Nature", 6.43, 0.05, ps.C_PUBLISHED),
        ("Andernach et al. 2025 (13.2')", 6.16, None, ps.C_PUBLISHED),
        (f"This work: {las:.1f}' on LoTSS 20\"", row["adopted"]["lp_mpc"], row.get("lp_z_err_mpc"), ps.C_OURS),
        ("This work, if the quasar were the host", lp_qso, None, ps.C_OURS),
    ]
    for i, (lab, v, e, col) in enumerate(entries):
        yy = len(entries) - i
        ax3.errorbar(v, yy, xerr=e, fmt="o", color=col, ms=6, capsize=2)
        ax3.text(2.05, yy, lab, va="center", ha="left", fontsize=7.5, color=ps.INK_SECONDARY)
    yy = 0
    ax3.errorbar(dep["median"], yy, xerr=[[dep["median"] - dep["p16"]], [dep["p84"] - dep["median"]]],
                 fmt="D", color=ps.C_OURS, ms=6, capsize=2)
    ax3.text(2.05, yy, f"true 3D length (median, 68%) for\nrandom orientation, xi = {dep['xi_used']:.2f}",
             va="center", ha="left", fontsize=7.5, color=ps.INK_SECONDARY)
    ax3.set_xlim(2.0, 9.0)
    ax3.set_xticks([5, 6, 7, 8, 9])
    ax3.set_ylim(-0.8, len(entries) + 0.8)
    ax3.set_yticks([])
    ax3.set_xlabel("length (Mpc)")
    ax3.set_title("How long is Porphyrion? Projected lengths", fontsize=10)
    ax3.spines["left"].set_visible(False)
    for v in (6.43,):
        ax3.axvline(v, color=ps.GRID, lw=0.8, zorder=0)
    fig.suptitle("Porphyrion: the longest black hole jets known, 7.5 billion years ago", x=0.03, ha="left",
                 fontsize=13)
    ps.caption(fig, "Lengths differ because catalogues associate different emission with the giant and use different "
               "redshifts. Our length is the 3-sigma extent of the emission joining the two lobe seeds "
               "(Planck 2018 cosmology, z = 0.896). Oei et al. 2024: projected 6.43 Mpc, 3D estimate 6.8 (median) / 7.28 (mean) Mpc.")
    return ps.save(fig, FIGURE_DIR / "02_porphyrion")


# --------------------------------------------------------------------------
# 03: the record claims
# --------------------------------------------------------------------------

def _short_name(name: str) -> str:
    """ILTJ003621.03+253500.9 -> J0036+2535, the usual short form."""
    if name.startswith("ILTJ") and len(name) > 15:
        sign = "+" if "+" in name else "-"
        ra, dec = name[4:].split(sign)
        return f"J{ra[:4]}{sign}{dec[:4]}"
    return name


def figure_records(measured: list[dict], art: dict) -> list[str]:
    rows = {r["key"]: r for r in measured}
    keys = [g.key for g in config.giants_in("records") if "gallery" not in g.groups]
    n = len(keys)
    ncol = 4
    nrow = int(np.ceil((n + 1) / ncol))
    fig = plt.figure(figsize=(16, 4.1 * nrow + 0.6))
    for i, k in enumerate(keys):
        g = config.giant_by_key(k)
        row = rows[k]
        ax = fig.add_subplot(nrow, ncol, i + 1)
        img = rm.smooth_to(rm.load_lotss(k, 20), 60.0)
        sig = rm.noise(img)
        _radio_panel(ax, img, sig, contour_mask=art[k]["adopted"].mask, vmax_sigma=40)
        (r1, d1), (r2, d2) = row["adopted"]["ends"]
        x, y = img.world_to_pixel(np.array([r1, r2]), np.array([d1, d2]))
        ax.plot(x, y, color=ps.BLUE, lw=1.0)
        hx, hy = img.world_to_pixel(g.ra, g.dec)
        ax.scatter(hx, hy, s=60, facecolors="none", edgecolors=ps.BLUE, linewidths=1.1)
        for s in g.seeds:
            sx, sy = img.world_to_pixel(*s)
            ax.scatter(sx, sy, marker="x", s=40, color=ps.VIOLET, linewidths=1.2)
        c = g.claims[0]
        ax.set_title(f"{g.name}" + (f"  ({_short_name(g.name)})" if g.name.startswith("ILTJ") else ""), fontsize=9)
        txt = (f"claimed: {c[2]:.2f} Mpc ({c[1]}'), {c[0]}\nsame association, ours: "
               f"{row['adopted']['las_arcmin']:.1f}' -> {row['adopted']['lp_mpc']:.2f} Mpc (z = {g.z:.3f})")
        if "conservative" in row:
            cons = row["conservative"]
            ax.contour(art[k]["conservative"].mask.astype(float), levels=[0.5], colors=[ps.AQUA], linewidths=1.0)
            txt += f"\nhost's own emission only (aqua): {cons['las_arcmin']:.1f}' -> {cons['lp_mpc']:.2f} Mpc"
        ax.text(0.02, 0.98, txt, transform=ax.transAxes, va="top", fontsize=6.8, color=ps.INK,
                bbox={"facecolor": "white", "alpha": 0.85, "lw": 0})
        ax.set_xticks([])
        ax.set_yticks([])
    # summary: claims against our measurement of the same association
    ax = fig.add_subplot(nrow, ncol, n + 1)
    allk = [g.key for g in config.giants_in("records")]
    for j, k in enumerate(sorted(allk, key=lambda kk: rows[kk]["adopted"]["lp_mpc"])):
        r = rows[k]
        for c in r["claims"]:
            ax.scatter(c["lp_mpc"], j, color=ps.C_PUBLISHED, s=18, zorder=3)
        ax.scatter(r["adopted"]["lp_mpc"], j, color=ps.C_OURS, s=28, marker="D", zorder=4)
        if "conservative" in r:
            ax.scatter(r["conservative"]["lp_mpc"], j, facecolors="none", edgecolors=ps.AQUA, s=34, marker="D",
                       linewidths=1.4, zorder=4)
            ax.plot([r["conservative"]["lp_mpc"], r["adopted"]["lp_mpc"]], [j, j], color=ps.GRID, lw=1, zorder=1)
        ax.text(0.15, j, _short_name(config.giant_by_key(k).name), fontsize=7, va="center", ha="left",
                color=ps.INK_SECONDARY)
    ax.set_xlim(0.0, 9.3)
    ax.set_xticks(range(2, 10))
    ax.set_yticks([])
    ax.set_xlabel("projected length (Mpc)")
    ax.set_title("claims (violet), ours for the same association (orange),\nthe host's own emission only (aqua)", fontsize=9)
    ax.spines["left"].set_visible(False)
    fig.suptitle("Giants claimed to exceed 5 Mpc: what LoTSS shows", x=0.01, ha="left", fontsize=13)
    ps.caption(fig, "LoTSS DR2 smoothed to 60\" for display. Orange outline: emission joined to the claimed ends (violet "
               "crosses, placed by us on the components each catalogue's length implies); blue: host and our largest extent; "
               "aqua: the emission clearly tied to the host. Every conservative length still exceeds 2 Mpc: these are real "
               "giants, but whether they are record breakers rests on detached sources that may be unrelated.")
    fig.subplots_adjust(hspace=0.25, wspace=0.08, top=0.93)
    return ps.save(fig, FIGURE_DIR / "03_record_claims")


# --------------------------------------------------------------------------
# 04: the census
# --------------------------------------------------------------------------

METHOD_COLORS = {"LOFAR Galaxy Zoo (citizen science)": ps.BLUE, "machine learning (Mostert et al.)": ps.ORANGE,
                 "by eye (Oei et al. 2023)": ps.AQUA, "earlier literature": ps.YELLOW}


def _plain_log_ticks(ax, axis: str, ticks):
    """Log axis labelled with plain numbers, no minor-tick labels."""
    from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator

    a = ax.xaxis if axis == "x" else ax.yaxis
    a.set_major_locator(FixedLocator(ticks))
    a.set_minor_formatter(NullFormatter())
    a.set_minor_locator(NullLocator())
    (ax.set_xticklabels if axis == "x" else ax.set_yticklabels)([f"{t:g}" for t in ticks])


def figure_census(census: dict) -> list[str]:
    g = catalogs.all_giants()
    fig = plt.figure(figsize=(16, 10.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], height_ratios=[1, 1], left=0.05, right=0.97,
                          top=0.9, bottom=0.07, wspace=0.18, hspace=0.32)
    # sky map, centred on RA 180 deg so the main LoTSS region is not split
    ax = fig.add_subplot(gs[0, 0], projection="mollweide")
    for method, col in METHOD_COLORS.items():
        m = g["found_by"] == method
        x = -np.radians(g.loc[m, "ra"].to_numpy() - 180.0)
        ax.scatter(x, np.radians(g.loc[m, "dec"]), s=1.3, color=col, alpha=0.6, lw=0)
    ax.set_xticks(np.radians(np.arange(-150, 151, 30)))
    ax.set_xticklabels([f"{(180 - d) % 360:.0f}°" for d in np.arange(-150, 151, 30)], fontsize=7)
    ax.grid(True, color=ps.GRID, lw=0.5)
    ax.set_title(f"All {len(g):,} known giants on the sky (equatorial; right ascension increases to the left)",
                 fontsize=10)
    # length vs redshift
    ax = fig.add_subplot(gs[0, 1])
    for method, col in METHOD_COLORS.items():
        m = g["found_by"] == method
        ax.scatter(g.loc[m, "z"], g.loc[m, "lp_mpc"], s=2, color=col, alpha=0.45, lw=0)
    zz = np.geomspace(0.003, 3, 400)
    for phi, lab, zl in ((5.0, "5': the by-eye search threshold", 1.6), (config.MOON_DIAMETER_ARCMIN, "31': the full Moon", 0.03)):
        lpc = cosmology.projected_length_mpc(phi, zz, config.PLANCK18)
        ax.plot(zz, lpc, color=ps.INK_SECONDARY, lw=1)
        yl = float(cosmology.projected_length_mpc(phi, zl, config.PLANCK18))
        ax.text(zl * (0.75 if phi < 10 else 1.12), yl * (1.08 if phi < 10 else 0.95), lab, fontsize=7.5,
                color=ps.INK_SECONDARY, ha="center" if phi < 10 else "left")
    ax.axhline(config.GIANT_THRESHOLD_MPC, color=ps.INK, lw=0.8)
    ax.text(0.0035, 0.73, "0.7 Mpc: a giant", fontsize=7.5)
    for k, lp in (("porphyrion", 6.43), ("alcyoneus", 4.99), ("3c236", 4.88), ("ngc315", 1.19)):
        gg = config.giant_by_key(k)
        ax.scatter(gg.z, lp, s=45, facecolors="none", edgecolors=ps.INK, lw=1)
        ax.annotate(gg.name, (gg.z, lp), (gg.z * 1.15, lp * 1.1), fontsize=7.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.003, 3)
    ax.set_ylim(0.6, 10)
    _plain_log_ticks(ax, "x", [0.003, 0.01, 0.03, 0.1, 0.3, 1, 3])
    _plain_log_ticks(ax, "y", [0.7, 1, 2, 3, 5, 10])
    ax.set_xlabel("redshift z")
    ax.set_ylabel("projected length (Mpc)")
    ax.set_title("Length against redshift: the survey sees a wedge", fontsize=10)
    # redshift distributions, with the shared legend
    ax = fig.add_subplot(gs[1, 0])
    bins = np.linspace(0, 2.5, 51)
    for method, col in METHOD_COLORS.items():
        m = g["found_by"] == method
        ax.hist(g.loc[m, "z"].dropna(), bins=bins, histtype="step", color=col, lw=1.6,
                label=f"{method}: {m.sum():,}")
    ax.set_xlabel("redshift z")
    ax.set_ylabel("giants per bin of 0.05")
    ax.legend(loc="upper right", fontsize=8, title="first found by", title_fontsize=8)
    ax.set_title(f"Redshifts: median {census['median_z']:.2f}; by eye finds the nearby ones, the crowd and the machine "
                 "the distant ones", fontsize=10)
    # the same angle in three universes
    ax = fig.add_subplot(gs[1, 1])
    zz = np.linspace(0.01, 2, 300)
    base = cosmology.angular_diameter_distance(zz, config.PLANCK18)
    cc = census["cosmology"]
    for c, col, dy in zip(config.COSMOLOGIES, (ps.INK, ps.C_OURS, ps.C_PUBLISHED), (-0.006, 0.006, 0.006)):
        ratio = cosmology.angular_diameter_distance(zz, c) / base
        ax.plot(zz, ratio, color=col, lw=1.6)
        ax.text(1.98, float(ratio[-1]) + dy, f"{c.name} (H0 = {c.H0:.1f}): {cc[c.name]['n_giants']:,} giants",
                fontsize=7.5, va="center", ha="right", color=ps.INK_SECONDARY)
    ax.set_xlim(0, 2)
    ax.set_xlabel("redshift z")
    ax.set_ylabel("length relative to Planck 2018")
    ax.set_title("Every length depends on the cosmology: how many giants stay above 0.7 Mpc", fontsize=10)
    fig.suptitle("The census of giant radio galaxies (Mostert et al. 2024)", x=0.01, ha="left", fontsize=13)
    ps.caption(fig, "Colours mark how each giant was first found (legend, bottom left). Lengths: Planck 2018 cosmology, "
               "as published. Curves (top right): the length a fixed angle corresponds to at each redshift. A 0.7-Mpc "
               f"giant never looks smaller than {census['min_angle_0p7mpc_arcmin']:.2f}' (at z = "
               f"{census['min_angle_redshift']:.2f}). Bottom right: lengths in this atlas's own cosmologies (D3, D1); the "
               "catalogue lists only giants above 0.7 Mpc for Planck, so counts can only fall.")
    return ps.save(fig, FIGURE_DIR / "04_census")


# --------------------------------------------------------------------------
# 05: the length distribution
# --------------------------------------------------------------------------

def figure_lengths(lres: dict, lart: dict) -> list[str]:
    import lengths

    post, pred, lp = lart["post"], lart["pred"], lart["lp"]
    r = lres["0.5"]
    fig = plt.figure(figsize=(16, 9.4))
    # survival function
    ax = fig.add_axes([0.06, 0.55, 0.42, 0.38])
    x, esf = lengths.empirical_survival(lp)
    ax.step(x, esf, where="post", color=ps.C_OBSERVED, lw=1.6, label=f"observed: {lp.size:,} giants at z < 0.5")
    band = np.percentile(pred["survival"], [16, 50, 84], axis=0)
    ax.fill_between(pred["lp"], band[0], band[2], color=ps.C_OURS, alpha=0.25, lw=0)
    ax.plot(pred["lp"], band[1], color=ps.C_OURS, lw=1.6, label="our fit with selection effects")
    xi_n = r["naive_xi"]
    ax.plot(pred["lp"], lengths.survival(pred["lp"], xi_n), color=ps.INK_SECONDARY, lw=1.2,
            label=f"pure power law, no selection: xi = {xi_n:.2f}")
    xi_i = r["posterior"]["xi"]["mean"]
    ax.plot(pred["lp"], lengths.survival(pred["lp"], xi_i), color=ps.C_PUBLISHED, lw=1.2,
            label=f"intrinsic distribution, xi = {xi_i:.2f}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.7, 8)
    ax.set_ylim(3e-4, 1.2)
    ax.set_xlabel("projected length l_p (Mpc)")
    ax.set_ylabel("fraction of giants longer than l_p")
    _plain_log_ticks(ax, "x", [0.7, 1, 2, 3, 5, 7])
    ax.legend(loc="lower left", fontsize=8)
    ax.set_title("How rare are long giants? Survival function of Oei et al.'s giants", fontsize=10)
    # posterior of xi
    ax = fig.add_axes([0.56, 0.55, 0.4, 0.38])
    xg, m = post.marginal("xi")
    ax.plot(xg, m, color=ps.C_OURS, lw=1.8, label=f"ours: {r['posterior']['xi']['mean']:.2f} +- {r['posterior']['xi']['sd']:.2f}")
    pub = r["published"]["xi"]
    yp = m.max() * 0.55
    ax.errorbar(pub[0], yp, xerr=pub[1], fmt="o", color=ps.C_PUBLISHED, capsize=3,
                label=f"Oei et al. 2023: {pub[0]} +- {pub[1]}")
    ax.axvline(xi_n, color=ps.INK_SECONDARY, lw=1)
    ax.text(xi_n + 0.05, m.max() * 0.9, "no selection\neffects", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_xlim(-5.5, -2.0)
    ax.set_xlabel("tail index xi  (P(L > l) ~ l^(xi + 1))")
    ax.set_ylabel("posterior density")
    ax.legend(loc="upper left", fontsize=8)
    ax.set_title(f"The tail index, after modelling selection (z < 0.5, {r['n']:,} giants)", fontsize=10)
    # completeness
    ax = fig.add_axes([0.06, 0.08, 0.42, 0.36])
    comp = pred["completeness"] / pred["completeness"].max(axis=1, keepdims=True)  # each draw's shape
    cb = np.percentile(comp, [16, 50, 84], axis=0)
    ax.fill_between(pred["lp"], cb[0], cb[2], color=ps.C_OURS, alpha=0.25, lw=0)
    ax.plot(pred["lp"], cb[1], color=ps.C_OURS, lw=1.6)
    ax.set_xscale("log")
    ax.set_xlim(0.7, 20)
    ax.set_ylim(0, 1.05)
    _plain_log_ticks(ax, "x", [0.7, 1, 2, 5, 10, 20])
    ax.set_xlabel("projected length (Mpc)")
    ax.set_ylabel("completeness (relative to the best-found length)")
    ax.set_title("Which giants LoTSS can find: small ones are missed for their angle,\nlong ones for their faint lobes",
                 fontsize=10)
    # degeneracy
    ax = fig.add_axes([0.56, 0.08, 0.4, 0.36])
    j = post.joint("xi", "b_ref")
    ax.contour(post.b_ref, post.xi, j, levels=_hpd_levels(j, (0.95, 0.68)), colors=[ps.C_OURS], linewidths=[1, 1.8])
    ax.set_ylim(-5.0, -1.8)
    ax.set_xlabel("median lobe surface brightness b_ref (Jy/deg^2, at 0.7 Mpc, z = 0)")
    ax.set_ylabel("tail index xi")
    ax.set_title("Why the uncertainty is large: brighter lobes (fewer missed giants)\nallow a steeper tail (68%, 95% regions)",
                 fontsize=10)
    fig.suptitle("The length distribution of giants", x=0.01, ha="left", fontsize=13)
    ps.caption(fig, "Model and sample as Oei et al. (2023, A&A 672, A163): Pareto-distributed lengths, random orientations, a "
               "fuzzy 5' angular threshold and lognormal lobe surface brightness falling as l^-2 (1+z)^-4, with their priors; "
               "posterior evaluated on a 4-D grid. Bands: 68% of posterior draws.")
    return ps.save(fig, FIGURE_DIR / "05_length_distribution")


def _hpd_levels(p: np.ndarray, masses) -> list[float]:
    flat = np.sort(p.ravel())[::-1]
    cum = np.cumsum(flat) / flat.sum()
    return sorted(float(flat[np.searchsorted(cum, mm)]) for mm in masses)


# --------------------------------------------------------------------------
# 06: lobe spectra and ages (3C 236)
# --------------------------------------------------------------------------

def figure_ageing(ares: dict, aart: dict, key: str = "3c236") -> list[str]:
    g = config.giant_by_key(key)
    res, art = ares[key], aart[key]
    sm = art["maps"]
    fig = plt.figure(figsize=(16, 6.4))
    ax = fig.add_axes([0.02, 0.1, 0.27, 0.76])
    shown = np.arcsinh(np.clip(sm.lotss, 0, None) / (2 * sm.sigma_lotss))
    ax.imshow(shown, origin="lower", cmap="gray_r", vmin=0, vmax=np.nanpercentile(shown, 99.8), interpolation="nearest")
    a = np.where(art["mask45"], sm.alpha, np.nan)
    im = ax.imshow(a, origin="lower", cmap="magma", vmin=-1.5, vmax=-0.4, interpolation="nearest")
    ax.contour(art["mask45"].astype(float), levels=[0.5], colors=[ps.INK_SECONDARY], linewidths=0.6)
    cax = fig.add_axes([0.295, 0.18, 0.008, 0.6])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("spectral index, 144 MHz to 1.4 GHz", fontsize=8)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"{g.name}: spectral index where both surveys\ndetect the lobes (darker = steeper = older)", fontsize=10)
    names = res.get("lobe_names", {"A": "A", "B": "B"})
    kpa = float(cosmology.kpc_per_arcmin(g.z, config.PLANCK18))
    ax2 = fig.add_axes([0.4, 0.1, 0.26, 0.76])
    ax3 = fig.add_axes([0.73, 0.1, 0.25, 0.76])
    for lobe, col in (("A", ps.BLUE), ("B", ps.ORANGE)):
        if lobe not in res["lobes"]:
            continue
        prof = res["lobes"][lobe]["profile"]
        d = np.array([p["distance_from_host_arcmin"] for p in prof]) * kpa / 1e3
        al = np.array([np.nan if p["alpha"] is None else p["alpha"] for p in prof], float)
        ae = np.array([np.nan if p["alpha_err"] is None else p["alpha_err"] for p in prof], float)
        up = np.array([np.nan if p.get("alpha_upper") is None else p["alpha_upper"] for p in prof], float)
        ax2.errorbar(d, al, yerr=ae, fmt="o", color=col, ms=4, capsize=1.5, label=f"{names[lobe]} lobe")
        ok = np.isfinite(up) & (up < -0.45)
        ax2.scatter(d[ok], up[ok], marker="v", color=col, s=22)
        de = np.array([p["distance_from_end_arcmin"] for p in prof]) * kpa
        age = np.array([np.nan if p["jp_age_myr"] is None else p["jp_age_myr"] for p in prof], float)
        lobe_res = res["lobes"][lobe]
        ax3.plot(de, age, "o", color=col, ms=4, label=f"{names[lobe]} lobe, B = {lobe_res['b_eq_ug']:.2f} uG")
        ci = lobe_res.get("ci_age_myr")
        if ci is not None and np.isfinite(ci):
            ax3.axhline(ci, color=col, lw=1.0)
            ax3.text(2950, ci + 1.5, f"{names[lobe]} lobe as a whole (CI): {ci:.0f} Myr", color=col, fontsize=7.5,
                     ha="right", va="bottom")
    ax2.axhline(config.ALPHA_INJ, color=ps.INK_SECONDARY, lw=0.8)
    ax2.text(0.05, config.ALPHA_INJ + 0.02, "assumed injection index", fontsize=7.5, color=ps.INK_SECONDARY)
    ax2.set_ylim(-1.65, -0.45)
    ax2.set_xlabel("distance from the host (Mpc, projected)")
    ax2.set_ylabel("spectral index alpha (S ~ nu^alpha)")
    ax2.legend(loc="lower right", fontsize=8)
    ax2.set_title("Spectra steepen toward the host: older electrons\nleft behind as the lobes advance (v: NVSS limits)",
                  fontsize=10)
    ax3.set_xlim(0, 3000)
    ax3.set_xlabel("distance from the lobe's outer end (kpc)")
    ax3.set_ylabel("radiative age (Myr)")
    ax3.legend(loc="upper left", fontsize=8)
    ax3.set_title("Synchrotron ages along the lobes (JP model)", fontsize=10)
    fig.suptitle(f"How old are {g.name}'s lobes? LoTSS (144 MHz) against NVSS (1.4 GHz) at 45\"", x=0.01, ha="left",
                 fontsize=13)
    ps.caption(fig, f"JP: Jaffe-Perola model, one burst of injection, for each strip one beam wide; CI: continuous "
               f"injection, for the whole lobe. Injection index {config.ALPHA_INJ}, minimum-energy field, inverse-Compton "
               f"losses on the CMB (B_CMB = {res['b_cmb_ug']:.2f} uG at z = {g.z:.3f}). Model ages: halving the field "
               "changes them by ~20%. NVSS misses smooth structure larger than ~15', so the faintest, oldest plasma may be "
               "steeper than shown.")
    return ps.save(fig, FIGURE_DIR / "06_lobe_ages")


# --------------------------------------------------------------------------
# 07: hosts and black holes
# --------------------------------------------------------------------------

def figure_hosts(hres: dict) -> list[str]:
    import hosts

    g = hosts.giant_hosts()
    c = hosts.control_hosts()
    fig = plt.figure(figsize=(16, 5.8))
    ax = fig.add_axes([0.05, 0.13, 0.27, 0.72])
    for i, b in enumerate(hres["bins"]):
        if "giants_median" not in b:
            continue
        xc = i
        ax.errorbar(xc - 0.1, b["control_median"], fmt="o", color=ps.C_CONTROL, ms=6)
        ax.errorbar(xc + 0.1, b["giants_median"], yerr=b["difference_err_dex"], fmt="o", color=ps.C_GIANT, ms=6,
                    capsize=2)
        ax.text(xc, 8.72, f"{b['n_giants']} / {b['n_control']}", ha="center", fontsize=7, color=ps.MUTED)
    ax.set_xticks(range(len(hres["bins"])))
    ax.set_xticklabels([f"{b['z'][0]:.2f}-{b['z'][1]:.2f}" for b in hres["bins"]])
    ax.set_xlabel("redshift range (giants / ordinary radio galaxies)")
    ax.set_ylabel("median log10 black hole mass (Msun)")
    ax.set_ylim(8.7, 9.15)
    ax.plot([], [], "o", color=ps.C_GIANT, label="giant hosts")
    ax.plot([], [], "o", color=ps.C_CONTROL, label="ordinary radio galaxies (L144 > 1e24 W/Hz)")
    ax.legend(loc="upper left", fontsize=8)
    m = hres["matched"]
    ax.set_title(f"Same black holes: matched in z and radio power,\ngiants differ by {m['median_difference_dex']:+.3f} +- "
                 f"{m['err_dex']:.3f} dex ({m['n_giants']:,} giants)", fontsize=10)
    ax = fig.add_axes([0.39, 0.13, 0.27, 0.72])
    s = g[g["sigma_ok"] & (g["z"] < 0.5)]
    ax.scatter(s["lp_mpc"], s["log_mbh"], s=4, color=ps.C_GIANT, alpha=0.4, lw=0)
    ax.set_xscale("log")
    _plain_log_ticks(ax, "x", [0.7, 1, 2, 3, 5])
    ax.set_xlabel("projected length (Mpc)")
    ax.set_ylabel("log10 black hole mass (Msun)")
    mv = hres["mbh_vs_length"]
    ax.set_title(f"A heavier black hole does not make a longer giant\n(Spearman rho = {mv['spearman_rho']:.3f}, "
                 f"p = {mv['p']:.2f}, {mv['n']:,} giants at z < 0.5)", fontsize=10)
    ax = fig.add_axes([0.72, 0.13, 0.26, 0.72])
    bins = np.linspace(7.8, 10.2, 49)
    ax.hist(c.loc[c["sigma_ok"], "log_mbh"], bins=bins, density=True, histtype="step", color=ps.C_CONTROL, lw=1.5,
            label="ordinary radio galaxies")
    ax.hist(g.loc[g["sigma_ok"], "log_mbh"], bins=bins, density=True, histtype="step", color=ps.C_GIANT, lw=1.8,
            label="giants")
    ga = hres["giants_all"]
    ax.set_xlabel("log10 black hole mass (Msun)")
    ax.set_ylabel("density")
    ax.legend(loc="upper left", fontsize=8)
    ax.set_title(f"{ga['n_with_mbh']:,} giant hosts with SDSS velocity dispersions:\nmedian "
                 f"{ga['median_mbh_1e9']:.2f} x 10^9 Msun", fontsize=10)
    fig.suptitle("The black holes behind the giants (M-sigma relation, Kormendy & Ho 2013)", x=0.01, ha="left", fontsize=13)
    ps.caption(fig, "Black hole masses from SDSS DR17 stellar velocity dispersions (70-420 km/s, error < 20%) with "
               "M = 3.09e8 Msun (sigma/200 km/s)^4.38; each mass carries ~0.3 dex intrinsic scatter. Ordinary radio galaxies: "
               "resolved LoTSS DR2 sources smaller than 0.7 Mpc with L144 > 1e24 W/Hz (Hardcastle et al. 2023).")
    return ps.save(fig, FIGURE_DIR / "07_hosts")


# --------------------------------------------------------------------------
# 08: environments
# --------------------------------------------------------------------------

def figure_environment(bres: dict, lumchk: dict, eres: dict) -> list[str]:
    import environment as env

    t = catalogs.oei2024()
    fig = plt.figure(figsize=(16, 5.8))
    ax = fig.add_axes([0.05, 0.13, 0.28, 0.72])
    for key, col, lab in (("giants", ps.C_GIANT, "luminous giants"), ("radiogalaxies", ps.C_CONTROL, "general radio galaxies")):
        v = np.sort(np.exp(t[key]["mu_fixed"] + t[key]["sigma2_fixed"] / 2))
        ax.step(v, np.arange(1, v.size + 1) / v.size, where="post", color=col, lw=1.6,
                label=f"{lab} (n = {v.size:,}; mean {v.mean():.2f})")
    ax.set_xscale("log")
    ax.set_xlabel("matter density around the host, 1 + delta (BORG SDSS)")
    ax.set_ylabel("cumulative fraction")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title(f"Oei et al.'s reconstruction: giants in denser places\n(KS p = {bres['fixed']['ks_p']:.1e})", fontsize=10)
    ax = fig.add_axes([0.42, 0.13, 0.3, 0.72])
    items = []
    lf = lumchk["fixed"]
    items.append(("BORG: giants / all radio galaxies", lf["giants_mean"] / lf["all_rg_mean"], None, ps.C_PUBLISHED))
    items.append((f"BORG: giants / radio galaxies > 1e24 W/Hz (n = {lf['n_bright']})",
                  lf["giants_mean"] / lf["bright_rg_mean"], None, ps.C_PUBLISHED))
    for k, lab in (("vs_radiogalaxies_z", "SDSS: giants / radio galaxies, same z"),
                   ("vs_radiogalaxies_z_L", "SDSS: giants / radio galaxies, same z and power"),
                   ("vs_lookalikes", "SDSS: giants / radio-quiet look-alikes")):
        v = eres.get(k)
        if v:
            items.append((f"{lab} ({v['n_giants']})", v["ratio"], v["err"], ps.C_OURS))
    for i, (lab, v, e, col) in enumerate(items):
        y = len(items) - i
        ax.errorbar(v, y, xerr=e, fmt="o", color=col, ms=6, capsize=2)
        ax.text(0.52, y + 0.28, lab, fontsize=7.5, color=ps.INK_SECONDARY)
    ax.axvline(1.0, color=ps.INK, lw=0.8)
    ax.set_xlim(0.5, 1.7)
    ax.set_ylim(0.3, len(items) + 0.8)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("density ratio (1 = same surroundings)")
    ax.set_title("The excess comes from radio power, not size", fontsize=10)
    ax = fig.add_axes([0.79, 0.13, 0.19, 0.72])
    classes = ("cluster", "filament", "sheet", "void")
    wg, wr = bres["web_giants"], bres["web_radiogalaxies"]
    y = np.arange(len(classes))
    ax.barh(y + 0.18, [wg[c] for c in classes], height=0.34, color=ps.C_GIANT, label="giants")
    ax.barh(y - 0.18, [wr[c] for c in classes], height=0.34, color=ps.C_CONTROL, label="radio galaxies")
    ax.set_yticks(y)
    ax.set_yticklabels(classes)
    ax.invert_yaxis()
    ax.set_xlabel("mean probability")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title("Cosmic web class (BORG tidal field)", fontsize=10)
    fig.suptitle("Where giants live", x=0.01, ha="left", fontsize=13)
    ps.caption(fig, "BORG: per-object densities of Oei et al. (2024, A&A 686, A137), smoothed on 2.9 Mpc/h, z < 0.16; their "
               "control radio galaxies given 144-MHz powers by matching to LoTSS DR2. SDSS: our counts of Main Galaxy Sample "
               f"galaxies within {config.ENV_RADIUS_MPC:g} Mpc and +-{config.ENV_DV_KMS:g} km/s, compared only at equal redshift.")
    return ps.save(fig, FIGURE_DIR / "08_environment")


# --------------------------------------------------------------------------
# 09: validation
# --------------------------------------------------------------------------

def figure_validation(census: dict, measured: list[dict], sim: dict, ares: dict, hres_oei: dict, bres: dict) -> list[str]:
    import hosts

    fig = plt.figure(figsize=(16, 9.6))
    # (a) catalogue lengths reproduced
    ax = fig.add_axes([0.05, 0.57, 0.26, 0.33])
    for name, df, col in (("Mostert et al. 2024", catalogs.all_giants(), ps.BLUE), ("Oei et al. 2023", catalogs.oei2023(), ps.ORANGE)):
        ok = np.isfinite(df["z"]) & (df["z"] > 0)
        r = cosmology.projected_length_mpc(df.loc[ok, "las_arcmin"], df.loc[ok, "z"], config.PLANCK18) / df.loc[ok, "lp_mpc"]
        ax.hist((r - 1) * 100, bins=np.linspace(-0.15, 0.15, 61), histtype="step", lw=1.5, color=col,
                label=f"{name}: median {np.median(r - 1) * 100:+.3f}%")
    ax.set_xlabel("our length / catalogue length - 1 (%)")
    ax.set_ylabel("giants")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_title("Catalogue lengths reproduced from angle and z\n(Planck 2018; spread = rounding of the angles)", fontsize=9.5)
    # (b) our angles vs catalogue angles
    ax = fig.add_axes([0.38, 0.57, 0.26, 0.33])
    for r in measured:
        for c in r["claims"]:
            if c["las_arcmin"] is None:
                continue
            for res_key, mk, col in (("20", "o", ps.BLUE), ("6", "s", ps.AQUA)):
                at = r["at"].get(res_key)
                if at and at["complete"]:
                    ax.scatter(c["las_arcmin"], at["las_arcmin"], marker=mk, s=22, color=col, zorder=3)
            ax.scatter(c["las_arcmin"], r["adopted"]["las_arcmin"], marker="D", s=14, facecolors="none",
                       edgecolors=ps.C_OURS, zorder=4)
    lim = [5, 75]
    ax.plot(lim, lim, color=ps.INK_SECONDARY, lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    _plain_log_ticks(ax, "x", [5, 10, 20, 40, 70])
    _plain_log_ticks(ax, "y", [5, 10, 20, 40, 70])
    ax.set_xlabel("catalogue angular length (arcmin)")
    ax.set_ylabel("our 3-sigma angular length (arcmin)")
    ax.scatter([], [], marker="o", color=ps.BLUE, label="LoTSS 20\"")
    ax.scatter([], [], marker="s", color=ps.AQUA, label="LoTSS 6\"")
    ax.scatter([], [], marker="D", facecolors="none", edgecolors=ps.C_OURS, label="adopted (20\" or 60\")")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_title("Our lengths against the catalogues' (same giants)", fontsize=9.5)
    # (c) simulations
    ax = fig.add_axes([0.71, 0.57, 0.26, 0.33])
    runs = sim["runs"]
    for i, r in enumerate(runs):
        x = r["xi_true"] + (i % 4 - 1.5) * 0.04
        ax.errorbar(x, r["xi_mean"], yerr=r["xi_sd"], fmt="o", color=ps.C_OURS, ms=4, capsize=1.5)
        ax.scatter(x, r["xi_mode"], marker="_", s=60, color=ps.INK)
        ax.scatter(x, r["naive_xi"], marker="x", s=20, color=ps.INK_SECONDARY)
    ax.plot([-4.3, -2.7], [-4.3, -2.7], color=ps.INK_SECONDARY, lw=0.8)
    ax.set_xlabel("true tail index (simulated catalogue)")
    ax.set_ylabel("recovered")
    ax.errorbar([], [], yerr=[], fmt="o", color=ps.C_OURS, label="posterior mean +- sd")
    ax.scatter([], [], marker="_", color=ps.INK, label="posterior mode")
    ax.scatter([], [], marker="x", color=ps.INK_SECONDARY, label="no selection modelled")
    ax.legend(loc="lower right", fontsize=7.5)
    ax.set_title(f"The fit recovers known answers ({len(runs)} simulated surveys;\nmode bias "
                 f"{sim['mode_bias']:+.2f}, 68% intervals cover {sim['coverage_68'] * 100:.0f}%)", fontsize=9.5)
    # (d) 3C 236 against Shulevski et al.
    ax = fig.add_axes([0.05, 0.08, 0.26, 0.33])
    r3 = ares.get("3c236")
    if r3:
        lit = r3["shulevski2019"]
        names = r3["lobe_names"]
        for lobe in ("A", "B"):
            nm = names[lobe]
            ours = r3["lobes"][lobe]
            for f_ours, e_ours, f_lit, mk in ((ours["flux_144_jy"], ours["flux_144_err"], lit[nm]["flux_143_mjy"] / 1e3, "o"),
                                              (ours["flux_1400_jy"], ours["flux_1400_err"], lit[nm]["flux_1400_mjy"] / 1e3, "s")):
                ax.errorbar(f_lit, f_ours, yerr=e_ours, fmt=mk, color=ps.BLUE if lobe == "A" else ps.ORANGE, capsize=2)
        ax.plot([0.3, 10], [0.3, 10], color=ps.INK_SECONDARY, lw=0.8)
        ax.set_xscale("log")
        ax.set_yscale("log")
        _plain_log_ticks(ax, "x", [0.5, 1, 2, 5, 10])
        _plain_log_ticks(ax, "y", [0.5, 1, 2, 5, 10])
        ax.set_xlabel("Shulevski et al. 2019 lobe flux density (Jy)")
        ax.set_ylabel("ours (Jy)")
        ax.text(0.35, 6, "circles 144 MHz, squares 1.4 GHz\nblue NW lobe, orange SE lobe", fontsize=7.5)
        mc = r3["model_check_with_their_fluxes"]
        ax.set_title(f"3C 236's lobes. Our CI model on their fluxes: NW {mc['NW']['ci_age_myr']:.0f} Myr\n"
                     f"(theirs 129), SE {mc['SE']['ci_age_myr']:.0f} Myr (theirs 117)", fontsize=9.5)
    # (e) black hole masses vs Oei et al.
    ax = fig.add_axes([0.38, 0.08, 0.26, 0.33])
    o = catalogs.oei2023()
    g = hosts.giant_hosts()
    from scipy.spatial import cKDTree

    def xyz(ra, dec):
        ra, dec = np.radians(ra), np.radians(dec)
        return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])

    d, i = cKDTree(xyz(g["ra"], g["dec"])).query(xyz(o["ra"], o["dec"]))
    ok = (d < np.radians(2 / 3600)) & np.isfinite(g["mbh"].to_numpy()[i]) & np.isfinite(o["mbh_1e9"].to_numpy())
    ax.scatter(np.log10(o["mbh_1e9"].to_numpy()[ok] * 1e9), np.log10(g["mbh"].to_numpy()[i][ok]), s=5,
               color=ps.BLUE, alpha=0.6, lw=0)
    ax.plot([7.8, 10.6], [7.8, 10.6], color=ps.INK_SECONDARY, lw=0.8)
    ax.set_xlabel("log10 M_BH, Oei et al. 2023 (SDSS DR12)")
    ax.set_ylabel("log10 M_BH, ours (SDSS DR17)")
    ax.set_title(f"Black hole masses: {hres_oei['n_both']} hosts, median offset {hres_oei['median_log_ratio']:+.3f} dex,\n"
                 f"scatter {hres_oei['scatter_dex']:.3f} dex", fontsize=9.5)
    # (f) BORG gamma fits
    ax = fig.add_axes([0.71, 0.08, 0.26, 0.33])
    rows = []
    for method in ("fixed", "flex"):
        for key in ("giants", "radiogalaxies"):
            pk, pt, pm = bres["published"][method][key]
            rows.append((f"{key[:-1] if key == 'giants' else 'radio galaxies'} ({method})", bres[method][key]["gamma_k"], pk,
                         bres[method][key]["gamma_mean"], pm))
    for j, (lab, k_ours, k_pub, m_ours, m_pub) in enumerate(rows):
        y = len(rows) - j
        ax.scatter(m_ours, y, color=ps.C_OURS, s=30, zorder=3)
        ax.scatter(m_pub, y, color=ps.C_PUBLISHED, s=30, marker="s", zorder=3)
        ax.text(2.0, y + 0.25, f"{lab}: k = {k_ours:.2f} (published {k_pub})", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_xlim(2, 9)
    ax.set_ylim(0.4, len(rows) + 0.8)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("mean density 1 + delta of the gamma fit")
    ax.scatter([], [], color=ps.C_OURS, label="ours, from their table")
    ax.scatter([], [], color=ps.C_PUBLISHED, marker="s", label="Oei et al. 2024")
    ax.legend(loc="lower right", fontsize=7.5)
    ax.set_title("Cosmic Web statistics recomputed", fontsize=9.5)
    fig.suptitle("Checks: every published number we lean on, recomputed", x=0.01, ha="left", fontsize=13)
    return ps.save(fig, FIGURE_DIR / "09_validation")


# --------------------------------------------------------------------------
# 11: DR2 against DR3
# --------------------------------------------------------------------------

def figure_releases(rel: dict) -> list[str]:
    rows = [o for o in rel["giants"] if o["complete"]]
    rows.sort(key=lambda o: o["las_dr2"])
    fig = plt.figure(figsize=(16, 6.0))
    names = [o["name"] for o in rows]
    y = np.arange(len(rows))
    for i, (key, xlabel, title, ref) in enumerate((
            ("las", "DR3 / DR2 length (same seeds, same resolution)", "Lengths", 1.0),
            ("flux", "DR3 / DR2 flux density (same sky footprint)", "Brightness", 1.0),
            ("noise", "DR3 / DR2 image noise", "Depth", 1.0))):
        ax = fig.add_axes([0.13 + i * 0.29, 0.12, 0.24, 0.74])
        if key == "las":
            v = np.array([o["las_dr3"] / o["las_dr2"] for o in rows])
        elif key == "flux":
            v = np.array([o["flux_ratio_dr3_dr2"] for o in rows])
        else:
            v = np.array([o["sigma_dr3_ujy"] / o["sigma_dr2_ujy"] for o in rows])
        ax.scatter(v, y, color=ps.C_OURS, s=30, zorder=3)
        ax.axvline(ref, color=ps.INK, lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(names if i == 0 else [], fontsize=8)
        ax.set_xlabel(xlabel)
        med = np.median(v)
        ax.set_title(f"{title}: median ratio {med:.3f}", fontsize=10)
        if key == "flux":
            ax.axvspan(1 - config.LOTSS_FLUX_SCALE_ERR, 1 + config.LOTSS_FLUX_SCALE_ERR, color=ps.GRID, alpha=0.6, lw=0,
                       zorder=0)
    fig.suptitle("The same giants in LoTSS DR2 and DR3", x=0.01, ha="left", fontsize=13)
    ps.caption(fig, "Each DR2 giant measured again on its DR3 20\" image with the same lobe seeds (60\" for the diffuse "
               "ones). Flux densities are summed over the same patch of sky (our DR2 3-sigma outline); grey band: DR2's "
               "10% flux-scale uncertainty. DR3 is reprocessed from the same and newer LOFAR data (Shimwell et al. 2026).")
    return ps.save(fig, FIGURE_DIR / "11_dr2_vs_dr3")
