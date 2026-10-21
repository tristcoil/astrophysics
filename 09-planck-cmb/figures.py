"""Every figure in the project.

Colour roles are fixed across all figures (see ``plotstyle``): blue is a
measurement -- the data being read or a number we measured -- orange is a
model curve, violet is a published value we compare against, grey is raw or
recessive context. Temperature maps use a diverging blue-red scale through
neutral grey, and masked sky is a darker flat grey.

Line charts are written as PNG and SVG. Sky maps are written as PNG only: an
SVG of a raster image is the same image wrapped in XML, several MB larger.
"""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize, to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle

import cmb_physics as cp
import plotstyle as ps
from config import (FIGURE_DIR, FIRAS_Y_LIMIT_95, H0_D1, H0_SH0ES, MAP_LIMIT_UK, PLANCK_BASELINE,
                    PLANCK_TT_LOWE, PUBLISHED_TT_PEAKS, T_CMB_ERR_K, T_CMB_K)

ps.use_style()

PROV_PLANCK = "Data: Planck 2018 (PR3), ESA Planck Legacy Archive."
PROV_FIRAS = ("Data: COBE/FIRAS monopole spectrum (Fixsen et al. 1996, Fixsen & Mather 2002), "
              "NASA LAMBDA.")

DL = r"$D_\ell$"
DL_UNITS = r"$D_\ell$  ($\mu$K$^2$)"
LABELS = {
    "ombh2": r"$\Omega_b h^2$  (baryons)",
    "omch2": r"$\Omega_c h^2$  (dark matter)",
    "theta100": r"$100\,\theta_{\rm MC}$  (acoustic angle)",
    "ns": r"$n_s$  (primordial tilt)",
    "H0": r"$H_0$  (km/s/Mpc)",
    "omegam": r"$\Omega_m$  (matter fraction)",
    "omegal": r"$\Omega_\Lambda$  (dark energy fraction)",
    "age": "Age  (Gyr)",
    "rstar": r"$r_*$  (sound horizon, Mpc)",
    "zstar": r"$z_*$  (last scattering)",
    "Ase2tau": r"$10^9 A_s e^{-2\tau}$  (amplitude)",
}


@dataclass
class Context:
    """Everything the figures need, assembled by ``run.py``."""

    firas: object
    firas_fit: dict
    full: object
    binned: object
    binned_bestfit: np.ndarray
    theory: object
    bestfit: dict
    peaks: list
    model_peaks: list
    harmonic: dict
    envelope: object = None
    maps: dict | None = None  # spectra from the sky maps
    display: dict | None = None
    validation: dict | None = None
    map_compare: dict | None = None
    responses: dict | None = None
    fit: dict | None = None


def _save(fig, name: str, *, svg: bool = True):
    return ps.save(fig, FIGURE_DIR / name, formats=("png", "svg") if svg else ("png",))


def _ring_scatter(ax, x, y, *, color, s=26, zorder=5, label=None, marker="o"):
    """Markers with a surface-coloured ring, per the mark spec."""
    return ax.scatter(x, y, s=s, color=color, edgecolor=ps.SURFACE, linewidth=1.2,
                      zorder=zorder, label=label, marker=marker)


def _label_box(ax, x, y, text, **kw):
    """Text that must stay legible on top of a busy map."""
    kw.setdefault("fontsize", 7.5)
    kw.setdefault("color", ps.INK)
    return ax.text(x, y, text, bbox={"boxstyle": "round,pad=0.25", "facecolor": ps.SURFACE,
                                     "edgecolor": "none", "alpha": 0.85}, **kw)


# --------------------------------------------------------------------------
# 01 -- the blackbody
# --------------------------------------------------------------------------


def figure_blackbody(ctx: Context):
    f, fit = ctx.firas, ctx.firas_fit
    t = fit["temperature_k"]
    sky = np.asarray(fit["sky_mjy_sr"])
    resid = np.asarray(fit["residual_kjy_sr"])

    fig, (ax, axr) = plt.subplots(2, 1, figsize=(7.4, 6.2), sharex=True,
                                  gridspec_kw={"height_ratios": [3.0, 1.35], "hspace": 0.08})
    nu = np.linspace(20, 720, 800)
    ax.plot(nu, cp.planck_intensity(nu, t), color=ps.C_MODEL, lw=2.0, zorder=3,
            label=f"Blackbody at T = {t:.4f} K (our fit)")
    _ring_scatter(ax, f.ghz, sky, color=ps.C_DATA, s=24, label="FIRAS measurements")
    ax.set_ylabel("Brightness  (MJy / sr)")
    ax.set_ylim(0, 500)
    ax.set_title("The CMB is a near-perfect blackbody at 2.725 K")
    peak = fit["peak_ghz"]
    ps.annotate(ax, f"peak at {peak:.1f} GHz  (wavelength {299.792458 / peak:.2f} mm)",
                xy=(peak, fit["peak_mjy_sr"]), xytext=(185, 455))
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.84), frameon=False)
    ax.text(0.97, 0.55,
            "43 frequencies. The 1σ errors (5–150 kJy/sr)\n"
            "are 0.001–0.04% of the peak: far smaller\n"
            "than the dots. The lower panel magnifies\n"
            "the differences 2,000-fold.",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color=ps.INK_SECONDARY)
    sec = ax.secondary_xaxis("top", functions=(lambda x: 299.792458 / np.maximum(x, 1e-3),
                                               lambda x: 299.792458 / np.maximum(x, 1e-3)))
    sec.set_xticks([5, 2, 1, 0.7, 0.5])
    sec.set_xticklabels(["5 mm", "2 mm", "1 mm", "0.7 mm", "0.5 mm"])
    sec.tick_params(colors=ps.MUTED, labelsize=8)
    sec.spines["top"].set_visible(False)

    # Residuals, with the shape a distortion at FIRAS's published limit would have.
    axr.axhline(0, color=ps.AXIS, lw=1.0, zorder=1)
    axr.plot(nu, FIRAS_Y_LIMIT_95 * cp.y_distortion(nu, t) * 1e3, color=ps.C_PUBLISHED,
             lw=1.4, zorder=2)
    axr.text(80, -80, "violet: a Compton-y distortion at FIRAS's 95% limit, |y| = 1.5×10⁻⁵",
             fontsize=7.5, color=ps.INK_SECONDARY)
    axr.errorbar(f.ghz, resid, yerr=f.sigma, fmt="none", ecolor=ps.C_DATA, elinewidth=0.9,
                 alpha=0.7, zorder=3)
    _ring_scatter(axr, f.ghz, resid, color=ps.C_DATA, s=14, zorder=4)
    axr.set_ylim(-120, 120)
    axr.set_ylabel("Data − blackbody\n(kJy / sr)")
    axr.set_xlabel("Frequency  (GHz)")
    axr.set_xlim(40, 700)
    y, ye = fit["y"], fit["y_err"]
    axr.text(0.02, 0.86, f"χ² = {fit['blackbody']['chi2']:.0f} for {fit['blackbody']['dof']} "
             f"degrees of freedom · best-fit y = ({y * 1e6:+.1f} ± {ye * 1e6:.1f})×10⁻⁶",
             transform=axr.transAxes, fontsize=7.5, color=ps.INK_SECONDARY)

    ps.caption(fig, f"{PROV_FIRAS}  The temperature's statistical error is "
               f"±{fit['temperature_stat_err_k'] * 1e6:.0f} µK, but FIRAS's absolute scale rests "
               f"on its calibrator's thermometry (±1 mK). Modern combined value: "
               f"{T_CMB_K} ± {T_CMB_ERR_K} K (Fixsen 2009). The y error is statistical only.")
    return _save(fig, "01_blackbody")


# --------------------------------------------------------------------------
# 02, 03 -- the sky
# --------------------------------------------------------------------------


def _mollweide_rgba(hmap, *, limit, mask=None, xsize=2400, cmap=None):
    """Project a HEALPix map to a Mollweide RGBA image (Galactic, l increasing left)."""
    import healpy as hp

    cmap = cmap or ps.diverging_cmap()
    proj = hp.projector.MollweideProj(xsize=xsize)
    nside = hp.npix2nside(len(hmap))
    to_pix = lambda x, y, z: hp.vec2pix(nside, x, y, z)  # noqa: E731
    img = proj.projmap(np.asarray(hmap, float), to_pix)
    outside = ~np.isfinite(img)
    rgba = cmap(Normalize(-limit, limit)(np.where(outside, 0.0, img)))
    if mask is not None:
        mimg = proj.projmap(np.asarray(mask, float), to_pix)
        rgba[(mimg < 0.5) & ~outside] = to_rgba(ps.C_MASK)
    rgba[outside] = (0, 0, 0, 0)
    return rgba, proj


def _draw_mollweide(ax, rgba, proj):
    ext = proj.get_extent()
    ax.imshow(rgba, origin="lower", extent=ext, interpolation="nearest")
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(2 * np.cos(th), np.sin(th), color=ps.AXIS, lw=0.8)
    ax.set_xlim(ext[0] - 0.02, ext[1] + 0.02)
    ax.set_ylim(ext[2] - 0.02, ext[3] + 0.02)
    ax.set_aspect("equal")
    ax.axis("off")


def _zoom_outline(proj):
    """The zoom patch's boundary, projected onto a Mollweide projector."""
    import healpy as hp

    from sky import ZOOM_CENTER, ZOOM_PIXELS, ZOOM_SIZE_DEG
    g = hp.projector.GnomonicProj(rot=(ZOOM_CENTER[0], ZOOM_CENTER[1], 0.0), xsize=ZOOM_PIXELS,
                                  ysize=ZOOM_PIXELS, reso=ZOOM_SIZE_DEG * 60 / ZOOM_PIXELS)
    x0, x1, y0, y1 = g.get_extent()
    s = np.linspace(0, 1, 60)
    xs = np.r_[x0 + (x1 - x0) * s, np.full(60, x1), x1 - (x1 - x0) * s, np.full(60, x0)]
    ys = np.r_[np.full(60, y0), y0 + (y1 - y0) * s, np.full(60, y1), y1 - (y1 - y0) * s]
    lon, lat = g.xy2ang(xs, ys, lonlat=True)
    return proj.ang2xy(lon, lat, lonlat=True)


def _colorbar(fig, rect, *, limit, label, ticks):
    cax = fig.add_axes(rect)
    sm = plt.cm.ScalarMappable(norm=Normalize(-limit, limit), cmap=ps.diverging_cmap())
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=ticks)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=8, colors=ps.INK_SECONDARY, length=0)
    cb.set_label(label, fontsize=8.5, color=ps.INK_SECONDARY)
    return cb


def figure_sky(ctx: Context):
    d, m = ctx.display, ctx.maps
    rgba, proj = _mollweide_rgba(d["sum"], limit=MAP_LIMIT_UK, mask=d["mask"])
    fig = plt.figure(figsize=(10.0, 6.4))
    ax = fig.add_axes([0.02, 0.17, 0.96, 0.72])
    _draw_mollweide(ax, rgba, proj)

    # The zoom patch, labelled from outside the ellipse.
    x, y = _zoom_outline(proj)
    ax.plot(x, y, color=ps.INK, lw=1.0)
    ps.annotate(ax, "zoomed in:\nnext figure", xy=(float(np.max(x)), float(np.mean(y))),
                xytext=(1.55, -0.92), color=ps.INK)
    ax.text(0.0, 0.03, "Galactic plane: masked", ha="center", va="center", fontsize=8,
            color=ps.SURFACE, fontweight="bold")
    # Planck's scanning left some sky unobserved in one half of the mission.
    xw, yw = proj.ang2xy(352.0, 47.0, lonlat=True)
    ps.annotate(ax, "not observed in both\nhalves of the mission",
                xy=(float(xw), float(yw)), xytext=(1.45, 0.86), color=ps.INK)

    fig.text(0.02, 0.955, "The cosmic microwave background: the universe at 380,000 years old",
             fontsize=12, fontweight="bold", color=ps.INK)
    fig.text(0.02, 0.925,
             f"Planck 2018 SMICA map, Galactic coordinates. Hot spots red, cold spots blue. "
             f"Typical deviation {m['rms_sum_uK']:.0f} µK on a {T_CMB_K:.4f} K background: "
             f"{m['rms_sum_uK'] / (T_CMB_K * 1e6) * 1e5:.0f} parts in 100,000.",
             fontsize=8.5, color=ps.INK_SECONDARY)
    _colorbar(fig, [0.30, 0.115, 0.40, 0.022], limit=MAP_LIMIT_UK,
              label="Temperature minus 2.7255 K  (µK)",
              ticks=[-300, -200, -100, 0, 100, 200, 300])
    fig.patches.append(Rectangle((0.73, 0.115), 0.018, 0.022, transform=fig.transFigure,
                                 facecolor=ps.C_MASK, edgecolor="none"))
    fig.text(0.753, 0.126, "masked", fontsize=8, color=ps.INK_SECONDARY, va="center")
    ps.caption(fig, f"{PROV_PLANCK} Average of the two half-mission SMICA maps, averaged to "
               f"6.9-arcmin pixels for display, monopole and residual dipole removed. Grey: "
               f"Planck's confidence mask (Galactic plane, bright sources: what remains there "
               f"after foreground removal is not CMB), plus sky that one of the two half-missions "
               f"did not observe -- the wedge and streaks. Unmasked, the map reaches "
               f"±{m['max_abs_fullsky_uK']:.0f} µK in the Galactic centre.")
    return _save(fig, "02_sky", svg=False)


def figure_zoom_noise(ctx: Context):
    d, m = ctx.display, ctx.maps
    cmap = ps.diverging_cmap()
    half = d["zoom_size_deg"] / 2
    ext = (half, -half, -half, half)  # l increases to the left, as on the sky

    fig = plt.figure(figsize=(10.0, 10.6))
    fig.text(0.05, 0.972, "Zoom in: the spots are real sky, not noise", fontsize=12,
             fontweight="bold", color=ps.INK)
    zl, zb = d["zoom_center"]
    fig.text(0.05, 0.935,
             f"The same {d['zoom_size_deg']:.0f}°×{d['zoom_size_deg']:.0f}° patch (centred on "
             f"Galactic l = {zl:.0f}°, b = {zb:.0f}°) at Planck's full resolution, both on the "
             f"same ±{MAP_LIMIT_UK:.0f} µK scale.\nThe typical spot is about a degree across: "
             f"that is the first acoustic peak, seen directly.",
             fontsize=8.5, color=ps.INK_SECONDARY, va="center")
    axes = [fig.add_axes([0.05, 0.50, 0.42, 0.397]), fig.add_axes([0.53, 0.50, 0.42, 0.397])]
    titles = ("The sky:  (hm1 + hm2) / 2", "Noise only:  (hm1 − hm2) / 2")
    for ax, key, title in zip(axes, ("zoom_sum", "zoom_diff"), titles):
        ax.imshow(d[key], origin="lower", extent=ext, cmap=cmap, vmin=-MAP_LIMIT_UK,
                  vmax=MAP_LIMIT_UK, interpolation="nearest")
        ax.set_title(title, fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(False)
    # Scale references: one degree, and the full Moon (0.52 deg across). The
    # x axis runs right-to-left, like the sky seen from inside.
    ax = axes[0]
    ax.add_patch(Rectangle((-5.4, -9.75), -4.3, 2.35, facecolor=ps.SURFACE, alpha=0.88,
                           edgecolor="none", zorder=3))
    ax.plot([-6.0, -7.0], [-8.2, -8.2], color=ps.INK, lw=2.5, solid_capstyle="butt", zorder=4)
    ax.text(-6.5, -9.2, "1°", ha="center", va="center", fontsize=8.5, color=ps.INK, zorder=4)
    ax.add_patch(Circle((-8.6, -8.2), 0.26, facecolor="none", edgecolor=ps.INK, lw=1.3,
                        zorder=4))
    ax.text(-8.6, -9.2, "Moon", ha="center", va="center", fontsize=8, color=ps.INK, zorder=4)
    _label_box(axes[1], 0.0, 0.0, f"rms {np.std(d['zoom_diff']):.0f} µK\n"
               f"vs {np.std(d['zoom_sum']):.0f} µK of sky", fontsize=9, ha="center",
               va="center", color=ps.INK_SECONDARY)
    _colorbar(fig, [0.35, 0.462, 0.30, 0.013], limit=MAP_LIMIT_UK, label="µK  (both patches)",
              ticks=[-300, -150, 0, 150, 300])

    # The full-sky noise map on a scale 20x tighter, to show its structure.
    limit = 15.0
    fig.text(0.05, 0.395, f"The whole noise map, on a ±{limit:.0f} µK scale (20× tighter)",
             fontsize=10, fontweight="bold", color=ps.INK)
    rgba, proj = _mollweide_rgba(d["diff"], limit=limit, mask=d["mask"], xsize=2000)
    axn = fig.add_axes([0.05, 0.10, 0.90, 0.28])
    _draw_mollweide(axn, rgba, proj)
    for lon, lat, lab, dx in ((96.38, 29.81, "north ecliptic pole", 0.10),
                              (276.38, -29.81, "south ecliptic pole", 0.10)):
        x, y = proj.ang2xy(lon, lat, lonlat=True)
        axn.plot(x, y, marker="+", ms=10, mew=1.6, color=ps.INK)
        _label_box(axn, x + dx, y + 0.06, lab)
    _colorbar(fig, [0.35, 0.075, 0.30, 0.013], limit=limit, label="µK",
              ticks=[-15, -10, -5, 0, 5, 10, 15])
    ps.caption(fig, f"{PROV_PLANCK} The two half-mission maps saw the same sky with independent "
               f"noise, so their half-difference contains no sky at all. Planck scanned in great "
               f"circles through the ecliptic poles, so the poles were crossed most often and are "
               f"least noisy. Noise rms in 1.7-arcmin pixels: {m['rms_diff_uK']:.1f} µK, against "
               f"{m['rms_sum_uK']:.0f} µK of sky.")
    return _save(fig, "03_zoom_and_noise", svg=False)


# --------------------------------------------------------------------------
# 04 -- the power spectrum
# --------------------------------------------------------------------------


def figure_power_spectrum(ctx: Context):
    full, binned, th = ctx.full, ctx.binned, ctx.theory
    fig, ax = plt.subplots(figsize=(10.0, 6.0))
    ps.set_ell_axis(ax)

    # Cosmic variance: the scatter a perfect experiment would still see.
    cv = cp.cosmic_variance_fraction(th.ell)
    ax.fill_between(th.ell, th.tt * (1 - cv), th.tt * (1 + cv), color=ps.C_MODEL, alpha=0.13,
                    lw=0, zorder=1)
    hi = full.ell >= 30
    ax.scatter(full.ell[hi], full.dl[hi], s=1.2, color=ps.C_RAW, alpha=0.45, lw=0, zorder=2)
    ax.plot(th.ell, th.tt, color=ps.C_MODEL, lw=1.8, zorder=3)
    lo = ~hi
    ax.errorbar(full.ell[lo], full.dl[lo], yerr=[full.err_lo[lo], full.err_hi[lo]], fmt="none",
                ecolor=ps.C_DATA, elinewidth=1.0, zorder=4)
    _ring_scatter(ax, full.ell[lo], full.dl[lo], color=ps.C_DATA, s=18)
    ax.errorbar(binned.ell, binned.dl, yerr=binned.err, fmt="none", ecolor=ps.C_DATA,
                elinewidth=1.0, zorder=4)
    _ring_scatter(ax, binned.ell, binned.dl, color=ps.C_DATA, s=18)

    ax.set_ylim(0, 6600)
    ax.set_xlabel(r"Multipole  $\ell$   (logarithmic below 30, linear above)")
    ax.set_ylabel(r"$D_\ell = \ell(\ell+1)\,C_\ell\,/\,2\pi$   ($\mu$K$^2$)")
    ps.add_degree_axis(ax)
    ax.set_title("The CMB power spectrum: how much the temperature varies on each angular scale",
                 pad=34)

    p = {q["n"]: q for q in ctx.peaks}
    ps.annotate(ax, "Largest scales, bigger\nthan the sound horizon:\na nearly flat plateau\n"
                "(Sachs–Wolfe effect)", xy=(12, 1050), xytext=(2.3, 4100))
    ps.annotate(ax, f"1st peak, ℓ ≈ {p[1]['ell']:.0f} (≈ {180 / p[1]['ell']:.1f}°):\n"
                "the sound horizon seen across\nthe universe: its angle says\nspace is flat",
                xy=(p[1]["ell"], p[1]["dl"]), xytext=(420, 5250))
    ps.annotate(ax, "2nd peak lower than 1st and 3rd:\nbaryons add inertia, favouring\n"
                "compressions over rarefactions", xy=(p[2]["ell"], p[2]["dl"]),
                xytext=(700, 3900))
    ps.annotate(ax, "3rd peak ≈ 2nd: set by the\ndark-matter density",
                xy=(p[3]["ell"], p[3]["dl"]), xytext=(1150, 2900))
    ps.annotate(ax, "Damping tail: photons diffuse,\nsmoothing out the finest ripples",
                xy=(1900, 330), xytext=(1780, 1500))

    handles = [
        Line2D([], [], color=ps.C_DATA, marker="o", ls="none", ms=5,
               label="Planck 2018 (bins of 30 in ℓ; single ℓ below 30)"),
        Line2D([], [], color=ps.C_RAW, marker="o", ls="none", ms=3, alpha=0.6,
               label="Planck 2018, every multipole"),
        Line2D([], [], color=ps.C_MODEL, lw=2, label="ΛCDM best fit (Planck 2018)"),
        Rectangle((0, 0), 1, 1, color=ps.C_MODEL, alpha=0.13,
                  label="Cosmic variance: the sky offers only 2ℓ+1 modes per ℓ"),
    ]
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=7.5)
    ps.caption(fig, f"{PROV_PLANCK} Planck's CMB-only TT spectrum (foregrounds subtracted by "
               f"Planck's likelihood model); ℓ < 30 from the Commander pipeline, with asymmetric "
               f"errors. Model: base ΛCDM best fit to TT,TE,EE + low-ℓ + lensing. A multipole ℓ "
               f"corresponds to an angle of roughly 180°/ℓ (top axis).")
    return _save(fig, "04_power_spectrum")


# --------------------------------------------------------------------------
# 05 -- the peaks
# --------------------------------------------------------------------------


def figure_peaks(ctx: Context):
    full, th = ctx.full, ctx.theory
    fig = plt.figure(figsize=(11.0, 5.2))
    ax = fig.add_axes([0.06, 0.13, 0.58, 0.72])
    axh = fig.add_axes([0.73, 0.13, 0.25, 0.72])

    m = (full.ell >= 100) & (full.ell <= 2450)
    ax.scatter(full.ell[m], full.dl[m], s=1.5, color=ps.C_RAW, alpha=0.5, lw=0, zorder=1)
    mt = (th.ell >= 100) & (th.ell <= 2450)
    ax.plot(th.ell[mt], th.tt[mt], color=ps.C_MODEL, lw=1.4, zorder=2)
    if ctx.envelope is not None:
        le = np.arange(1300, 2451)
        ax.plot(le, ctx.envelope(le), color=ps.INK_SECONDARY, lw=0.9, zorder=2)
        ps.annotate(ax, "smooth damping envelope,\ndivided out to find peaks 5–7",
                    xy=(1620, float(ctx.envelope(1620))), xytext=(1560, 1750))
    pub = {row[0]: row for row in PUBLISHED_TT_PEAKS}
    for q in ctx.peaks:
        n = q["n"]
        ax.errorbar(q["ell"], q["dl"], xerr=q["ell_err"], yerr=q["dl_err"], fmt="none",
                    ecolor=ps.C_OURS, elinewidth=1.3, zorder=4)
        _ring_scatter(ax, q["ell"], q["dl"], color=ps.C_OURS, s=40, zorder=5)
        r = pub[n]
        ax.scatter(r[1], r[3], s=46, facecolor="none", edgecolor=ps.C_PUBLISHED, lw=1.4,
                   zorder=6)
        ax.text(q["ell"], q["dl"] + max(260, 0.07 * q["dl"]), str(n), ha="center",
                fontsize=9, color=ps.INK, fontweight="bold")
    ax.set_xlim(100, 2450)
    ax.set_ylim(0, 6300)
    ax.set_xlabel(r"Multipole  $\ell$")
    ax.set_ylabel(DL_UNITS)
    ax.set_title("Seven acoustic peaks, located in Planck's spectrum")
    handles = [
        Line2D([], [], color=ps.C_OURS, marker="o", ls="none", ms=6, label="Peaks we measured"),
        Line2D([], [], color=ps.C_PUBLISHED, marker="o", ls="none", ms=6, mfc="none",
               label="Planck 2015 published peaks"),
        Line2D([], [], color=ps.C_MODEL, lw=1.6, label="ΛCDM best fit"),
        Line2D([], [], color=ps.C_RAW, marker="o", ls="none", ms=3, label="Planck data, every ℓ"),
    ]
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=7.5)

    # The harmonic series.
    h = ctx.harmonic
    n = np.array([q["n"] for q in ctx.peaks], float)
    ell = np.array([q["ell"] for q in ctx.peaks])
    err = np.array([q["ell_err"] for q in ctx.peaks])
    nn = np.linspace(0.5, 7.5, 50)
    axh.plot(nn, h["ell_A"] * (nn - h["phi"]), color=ps.C_MODEL, lw=1.6, zorder=2)
    iso = ell[0] * (2 * np.arange(1, 5) - 1)
    axh.plot(np.arange(1, 5), iso, color=ps.C_ALT, lw=1.4, marker="s", ms=4, zorder=2)
    axh.text(3.25, iso[2] + 150, "if the oscillations\nhad started from\nzero (isocurvature):\n"
             "peaks at 1 : 3 : 5 : 7", fontsize=7.5, color=ps.INK_SECONDARY, ha="right",
             va="bottom")
    axh.errorbar(n, ell, yerr=err, fmt="none", ecolor=ps.C_OURS, elinewidth=1.2, zorder=3)
    _ring_scatter(axh, n, ell, color=ps.C_OURS, s=34)
    axh.set_xlabel(r"Peak number  $n$")
    axh.set_ylabel(r"Peak position  $\ell_n$")
    axh.set_xlim(0.5, 7.5)
    axh.set_ylim(0, 3000)
    axh.set_title("A harmonic series")
    axh.text(0.97, 0.03,
             f"$\\ell_n = \\ell_A\\,(n - \\varphi)$\n$\\ell_A$ = {h['ell_A']:.0f} ± "
             f"{h['ell_A_err']:.0f}\n$\\varphi$ = {h['phi']:.2f} ± {h['phi_err']:.2f}\n"
             f"spacing ⇒ θ = {h['theta_deg']:.3f}°",
             transform=axh.transAxes, ha="right", va="bottom", fontsize=8, color=ps.INK_SECONDARY)
    ps.caption(fig, f"{PROV_PLANCK} Peaks 1–4: position and height of the maximum of a Gaussian-"
               f"plus-line fit to every-ℓ data within ±100 of each peak. Peaks 5–7: the same, "
               f"after dividing by a smooth fit to the damping tail, which otherwise hides them "
               f"(the 7th is not even a local maximum of D_ℓ). Errors from 300 noisy re-runs of "
               f"the whole measurement; for peaks 5–7 they include the spread over envelope "
               f"choices. Planck 2015 values: Planck Collaboration I (2016), Table 7.")
    return _save(fig, "05_peaks")


# --------------------------------------------------------------------------
# 06 -- our own spectrum from the maps
# --------------------------------------------------------------------------


def figure_map_spectrum(ctx: Context):
    s, c, v = ctx.maps, ctx.map_compare, ctx.validation
    bf_b, th = ctx.binned_bestfit, ctx.theory
    ell = s["ell"]

    fig = plt.figure(figsize=(10.0, 10.4))
    ax = fig.add_axes([0.10, 0.50, 0.86, 0.40])
    axd = fig.add_axes([0.10, 0.27, 0.86, 0.18])
    axv = fig.add_axes([0.10, 0.07, 0.86, 0.14])
    for a in (ax, axd, axv):
        ps.set_ell_axis(a)

    # The naive spectrum, binned like everything else above l = 30.
    naive_dl = cp.dl_from_cl(ell, s["naive_raw"])
    lnb, naive_b = cp.planck_bin(ell[:2509], naive_dl[:2509])
    ax.plot(np.r_[ell[2:30], lnb], np.r_[naive_dl[2:30], naive_b], color=ps.C_RAW, lw=1.3,
            zorder=2)
    ax.plot(c["l_eff"], c["noise"], color=ps.C_NOISE, lw=1.8, zorder=2)
    ax.plot(th.ell, th.tt, color=ps.C_MODEL, lw=1.4, zorder=3)
    low = (ell >= 2) & (ell < 30)
    ax.scatter(ell[low], cp.dl_from_cl(ell[low], s["cross"][low]), s=16, color=ps.C_OURS,
               edgecolor=ps.SURFACE, lw=1.0, zorder=5)
    ax.errorbar(c["l_eff"], c["cross"], yerr=c["err"], fmt="none", ecolor=ps.C_OURS,
                elinewidth=1.0, zorder=4)
    _ring_scatter(ax, c["l_eff"], c["cross"], color=ps.C_OURS, s=18)
    ax.set_ylim(0, 6600)
    ax.set_ylabel(DL_UNITS)
    ax.set_title("Our own power spectrum, measured from the raw sky maps", pad=34)
    ps.add_degree_axis(ax)
    ax.tick_params(labelbottom=False)
    ps.annotate(ax, "noise level of the map\n(measured from hm1 − hm2)",
                xy=(2400, float(np.interp(2400, c["l_eff"], c["noise"]))), xytext=(1900, 2300))
    ps.annotate(ax, "naive: whole sky,\nno corrections at all",
                xy=(760, float(np.interp(760, lnb, naive_b))), xytext=(560, 900))
    handles = [
        Line2D([], [], color=ps.C_OURS, marker="o", ls="none", ms=5,
               label="Ours: hm1 × hm2 cross-spectrum, masked and corrected"),
        Line2D([], [], color=ps.C_MODEL, lw=2, label="ΛCDM best fit (Planck 2018)"),
        Line2D([], [], color=ps.C_RAW, lw=1.3, label="Naive spectrum of the whole-sky map"),
        Line2D([], [], color=ps.C_NOISE, lw=2, label="Noise level of the sum map"),
    ]
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=7.5)

    # Difference from Planck's best fit, in uK^2: the same units at every scale.
    diff = c["cross"] - bf_b
    axd.axhline(0, color=ps.AXIS, lw=1.0)
    lp = np.linspace(30, 2508, 400)
    axd.plot(lp, c["poisson_amp"] * (lp / 3000) ** 2, color=ps.C_ALT, lw=1.6, zorder=2)
    ps.annotate(axd, f"faint point sources: $D_\\ell \\propto \\ell^2$,\n"
                f"{c['poisson_amp']:.0f} ± {c['poisson_amp_err']:.0f} µK² at ℓ = 3000 (fitted)",
                xy=(2100, c["poisson_amp"] * (2100 / 3000) ** 2), xytext=(1450, 150))
    axd.errorbar(c["l_eff"], diff, yerr=c["err"], fmt="none", ecolor=ps.C_OURS, elinewidth=1.0)
    _ring_scatter(axd, c["l_eff"], diff, color=ps.C_OURS, s=14)
    axd.set_ylim(-260, 260)
    axd.set_ylabel("Ours − Planck\nbest fit  (µK²)")
    axd.tick_params(labelbottom=False)

    # The estimator on a simulated sky.
    lv, full_b = cp.planck_bin(v["ell"][:2509], cp.dl_from_cl(v["ell"], v["cl_full"])[:2509])
    _, simple_b = cp.planck_bin(v["ell"][:2509],
                                cp.dl_from_cl(v["ell"], v["cl_masked_fsky"])[:2509])
    _, corr_b = cp.planck_bin(v["ell"][:2509], cp.dl_from_cl(v["ell"], v["cl_masked"])[:2509])
    axv.axhline(0, color=ps.AXIS, lw=1.0)
    axv.plot(lv, 100 * (simple_b / full_b - 1), color=ps.C_RAW, lw=1.2, marker="o", ms=2.5)
    axv.plot(lv, 100 * (corr_b / full_b - 1), color=ps.C_OURS, lw=1.4, marker="o", ms=2.5)
    i = int(np.argmax(simple_b / full_b))
    axv.text(lv[i] + 40, 100 * (simple_b[i] / full_b[i] - 1) - 1.5,
             "dividing by sky area alone fills in\nthe troughs (+11% at the first)",
             fontsize=7.5, color=ps.INK_SECONDARY, va="top")
    axv.text(2480, -4.8, "correlation-function correction (used above)", fontsize=7.5,
             color=ps.INK_SECONDARY, ha="right")
    axv.set_ylim(-6, 14)
    axv.set_ylabel("SIMULATED sky:\nmethod error (%)")
    axv.set_xlabel(r"Multipole  $\ell$")

    ps.caption(fig, f"{PROV_PLANCK} SMICA half-mission maps; common confidence mask with "
               f"missing pixels removed (f_sky = {s['fsky_binary']:.3f}), tapered over 30′; mask "
               f"coupling removed via correlation functions; beam and pixel window divided out; "
               f"Knox error bars. Above ℓ ≈ 1300 our spectrum exceeds Planck's CMB-only spectrum "
               f"by an amount matching unresolved point sources, which Planck's likelihood models "
               f"and removes; we do not. Bottom: a Gaussian sky drawn from the best-fit spectrum, "
               f"processed identically and compared with its own full-sky spectrum.")
    return _save(fig, "06_map_spectrum")


# --------------------------------------------------------------------------
# 07 -- what the peaks measure
# --------------------------------------------------------------------------


def figure_physics(ctx: Context):
    r, binned = ctx.responses, ctx.binned
    ell = r["ell"]
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6), sharey=True,
                             gridspec_kw={"hspace": 0.30, "wspace": 0.06})
    light, dark = ps.BLUE_RAMP[4], ps.BLUE_RAMP[11]
    panels = [
        ("ombh2", r"Baryons, $\Omega_b h^2$",
         "More baryons: the odd peaks (compressions)\nrise against the even ones.",
         lambda v: rf"$\Omega_b h^2$ = {v:.4f}"),
        ("omch2", r"Dark matter, $\Omega_c h^2$",
         "Less dark matter: gravity decays while the\nwaves oscillate, driving all peaks higher.",
         lambda v: rf"$\Omega_c h^2$ = {v:.3f}"),
        ("omk", r"Curvature, $\Omega_k$",
         "Curvature bends the light's path: a closed\nuniverse (Ω_k < 0) moves peaks left.",
         lambda v: rf"$\Omega_k$ = {v:+.2f}"),
        ("ns", r"Primordial tilt, $n_s$",
         "Tilt: more (n_s > 1) or less (n_s < 1) power\non small scales than on large ones.",
         lambda v: rf"$n_s$ = {v:.2f}"),
    ]
    for ax, (key, title, note, fmt) in zip(axes.flat, panels):
        ps.set_ell_axis(ax)
        ax.plot(ell[2:], r[f"{key}_0"][2:], color=light, lw=1.5, zorder=3)
        ax.plot(ell[2:], r[f"{key}_1"][2:], color=dark, lw=1.5, zorder=3)
        ax.plot(ell[2:], r["base"][2:], color=ps.C_MODEL, lw=1.7, zorder=4)
        ax.errorbar(binned.ell, binned.dl, yerr=binned.err, fmt="none", ecolor=ps.INK_SECONDARY,
                    elinewidth=0.8, zorder=5)
        ax.scatter(binned.ell, binned.dl, s=6, color=ps.INK_SECONDARY, zorder=6, lw=0)
        ax.set_title(title, fontsize=10)
        ax.text(0.03, 0.97, note, transform=ax.transAxes, fontsize=7.5, color=ps.INK_SECONDARY,
                va="top")
        v0, v1 = float(r[f"{key}_0_value"]), float(r[f"{key}_1_value"])
        handles = [Line2D([], [], color=light, lw=2, label=fmt(v0)),
                   Line2D([], [], color=ps.C_MODEL, lw=2, label="Planck best fit"),
                   Line2D([], [], color=dark, lw=2, label=fmt(v1))]
        ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=7.2)
        ax.set_ylim(0, 8200)
    for ax in axes[1]:
        ax.set_xlabel(r"Multipole  $\ell$")
    for ax in axes[:, 0]:
        ax.set_ylabel(DL_UNITS)
    fig.suptitle("What the peaks measure: change one ingredient, watch the spectrum respond",
                 x=0.06, ha="left", y=0.965)
    ps.caption(fig, f"{PROV_PLANCK} Dark grey points: Planck's binned TT spectrum. Model spectra "
               f"computed with CAMB, varying one parameter from Planck's best fit. Baryons, dark "
               f"matter and tilt are varied with the acoustic angle held fixed, so peaks keep their "
               f"places and only heights respond; curvature is varied at fixed H₀, since it acts by "
               f"changing the angle. Light and dark blue: the lower and higher value.")
    return _save(fig, "07_what_the_peaks_measure")


# --------------------------------------------------------------------------
# 08 -- the fit
# --------------------------------------------------------------------------


def _decimals(v: float, err: float) -> int:
    """Enough decimals to show two significant figures of the error."""
    return int(max(0, min(6, 1 - np.floor(np.log10(abs(err))))))


def figure_fit(ctx: Context):
    fit, binned = ctx.fit, ctx.binned
    model_b = np.asarray(fit["model_binned"])
    fig = plt.figure(figsize=(11.4, 7.8))
    axr = fig.add_axes([0.07, 0.66, 0.90, 0.24])
    axp = fig.add_axes([0.21, 0.08, 0.29, 0.46])
    axh = fig.add_axes([0.68, 0.08, 0.29, 0.46])

    # Residuals.
    ps.set_ell_axis(axr, lmin=30, lmax=2510, log_fraction=0.0)
    axr.axhline(0, color=ps.AXIS, lw=1.0)
    axr.errorbar(binned.ell, binned.dl - model_b, yerr=binned.err, fmt="none",
                 ecolor=ps.C_DATA, elinewidth=1.0)
    _ring_scatter(axr, binned.ell, binned.dl - model_b, color=ps.C_DATA, s=16)
    axr.set_ylabel("Planck data − our\nbest-fit model  (µK²)")
    axr.set_xlabel(r"Multipole  $\ell$")
    axr.set_title(f"Our five-parameter ΛCDM fit to Planck's TT spectrum: χ² = {fit['chi2']:.1f} "
                  f"for {fit['dof']} degrees of freedom")
    axr.set_ylim(-230, 230)

    # Parameters, as shifts from Planck TT+lowE in units of Planck's error.
    keys = ["ombh2", "omch2", "theta100", "ns", "Ase2tau", "H0", "omegam", "age", "rstar"]
    vals = {**fit["parameters"], **fit["derived"]}
    ys = np.arange(len(keys))[::-1]
    axp.axvspan(-1, 1, color=ps.C_PUBLISHED, alpha=0.10, lw=0)
    axp.axvline(0, color=ps.C_PUBLISHED, lw=1.2)
    for y, k in zip(ys, keys):
        ref, ref_err = PLANCK_TT_LOWE[k]
        shift = (vals[k]["value"] - ref) / ref_err
        e = vals[k]["error"] / ref_err
        axp.errorbar(shift, y, xerr=e, fmt="none", ecolor=ps.C_OURS, elinewidth=1.4)
        _ring_scatter(axp, shift, y, color=ps.C_OURS, s=32)
        dec = _decimals(vals[k]["value"], vals[k]["error"])
        axp.text(3.05, y, f"{vals[k]['value']:.{dec}f} ± {vals[k]['error']:.{dec}f}",
                 va="center", ha="left", fontsize=7.5, color=ps.INK_SECONDARY)
    axp.set_yticks(ys)
    axp.set_yticklabels([LABELS[k] for k in keys], fontsize=8)
    axp.set_xlim(-2.8, 2.8)
    axp.set_ylim(-0.7, len(keys) + 0.2)
    axp.grid(axis="y", visible=False)
    axp.set_xlabel("Shift from Planck TT+lowE  (in Planck's σ)")
    axp.set_title("Ours vs Planck's own fit", fontsize=10)
    axp.text(0, len(keys) - 0.35, "Planck ±1σ", ha="center", fontsize=7.5,
             color=ps.INK_SECONDARY)

    # H0: the tension, measured on both sides within this repository. Blue marks
    # our own measurements (this project, and project D1); violet, published ones.
    h = fit["derived"]["H0"]
    rows = [("D3, this project: CMB (our fit to Planck TT)", h["value"], h["error"], ps.C_OURS),
            ("Planck 2018, TT+lowE", *PLANCK_TT_LOWE["H0"], ps.C_PUBLISHED),
            ("Planck 2018, TT,TE,EE+lowE+lensing", *PLANCK_BASELINE["H0"], ps.C_PUBLISHED),
            ("D1, this repository: Pantheon+ supernovae", *H0_D1, ps.C_OURS),
            ("SH0ES 2022: Cepheids + supernovae", *H0_SH0ES, ps.C_PUBLISHED)]
    ys_h = [5.4, 4.4, 3.4, 1.4, 0.4]  # a gap separates the early and late universe
    for y, (lab, val, err, col) in zip(ys_h, rows):
        axh.errorbar(val, y, xerr=err, fmt="none", ecolor=col, elinewidth=1.6)
        _ring_scatter(axh, val, y, color=col, s=38)
        axh.text(63.2, y + 0.33, lab, fontsize=8, color=ps.INK_SECONDARY, va="center")
        axh.text(val + err + 0.25, y, f"{val:.2f} ± {err:.2f}", fontsize=7.5,
                 color=ps.INK_SECONDARY, va="center")
    for y, lab in ((6.15, "EARLY UNIVERSE: the CMB"), (2.15, "TODAY: the distance ladder")):
        axh.text(63.2, y, lab, fontsize=7.5, color=ps.INK, fontweight="bold", va="center")
    axh.axhline(2.55, color=ps.AXIS, lw=0.8)
    axh.set_yticks([])
    axh.set_xlim(63, 77.5)
    axh.set_ylim(-0.2, 6.5)
    axh.grid(axis="y", visible=False)
    axh.set_xlabel(r"$H_0$  (km/s/Mpc)")
    t_d1 = (H0_D1[0] - h["value"]) / np.hypot(H0_D1[1], h["error"])
    axh.set_title(f"The Hubble tension: {t_d1:.1f}σ between D3 and D1", fontsize=10)
    ps.caption(fig, f"{PROV_PLANCK} Fit: CAMB model against Planck's binned CMB-only TT spectrum "
               f"(83 bins, ℓ = 30–2508), independent Gaussian errors; τ fixed at "
               f"{fit['tau_fixed']} (Planck's TT+lowE value) because temperature data above "
               f"ℓ = 30 constrain only A_s e^(−2τ). Derived parameters assume a flat universe with "
               f"one 0.06 eV neutrino. Planck values: 2018 parameter tables, 68% limits. D1: "
               f"project D1 in this repository (05-supernova-cosmology). SH0ES: Riess et al. "
               f"(2022).")
    return _save(fig, "08_fit")


# --------------------------------------------------------------------------
# 00 -- scorecard
# --------------------------------------------------------------------------


def figure_summary(ctx: Context):
    rows = []  # (group, label, ours, ours_err, ref, ref_err, fmt)
    f = ctx.firas_fit
    rows.append(("vs Fixsen (2009)", r"CMB temperature $T_0$ (K)", f["temperature_k"], 0.001,
                 T_CMB_K, T_CMB_ERR_K, ".4f"))
    pub = {r[0]: r for r in PUBLISHED_TT_PEAKS}
    for q in ctx.peaks[:4]:
        r = pub[q["n"]]
        rows.append(("vs Planck 2015 peak table", f"Peak {q['n']} position ℓ", q["ell"],
                     q["ell_err"], r[1], r[2], ".1f"))
    for q in ctx.peaks[:4]:
        r = pub[q["n"]]
        rows.append(("vs Planck 2015 peak table", f"Peak {q['n']} height (µK²)", q["dl"],
                     q["dl_err"], r[3], r[4], ".0f"))
    if ctx.fit is not None:
        vals = {**ctx.fit["parameters"], **ctx.fit["derived"]}
        for k, spec in (("ombh2", ".5f"), ("omch2", ".4f"), ("theta100", ".5f"), ("ns", ".4f"),
                        ("H0", ".2f"), ("omegam", ".3f"), ("age", ".3f")):
            ref, ref_err = PLANCK_TT_LOWE[k]
            rows.append(("vs Planck 2018 TT+lowE", LABELS[k], vals[k]["value"], vals[k]["error"],
                         ref, ref_err, spec))

    fig, ax = plt.subplots(figsize=(9.6, 0.36 * len(rows) + 2.3))
    ys = np.arange(len(rows))[::-1]
    ax.axvspan(-1, 1, color=ps.GRID, alpha=0.6, lw=0)
    ax.axvline(0, color=ps.INK_SECONDARY, lw=1.0)
    for y, (_, lab, v, e, ref, re, spec) in zip(ys, rows):
        shift = (v - ref) / np.hypot(e, re)
        ax.plot([0, shift], [y, y], color=ps.C_OURS, lw=2.0, solid_capstyle="round")
        _ring_scatter(ax, shift, y, color=ps.C_OURS, s=30)
        ax.text(3.12, y, f"{v:{spec}}  vs  {ref:{spec}}", va="center", fontsize=7.6,
                color=ps.INK_SECONDARY)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[1] for r in rows], fontsize=8)
    ax.grid(axis="y", visible=False)
    ax.set_xlim(-3, 3)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xlabel("Difference from the published value, in units of the combined uncertainty")
    ax.set_title("Our measurements vs the literature")
    ax.text(3.12, len(rows) - 0.55, "ours  vs  published", fontsize=7.6, color=ps.INK,
            fontweight="bold")
    groups = [r[0] for r in rows]
    for i in range(len(rows)):
        if i == 0 or groups[i] != groups[i - 1]:
            if i:
                ax.axhline(ys[i] + 0.5, color=ps.AXIS, lw=0.8)
            ax.text(-2.95, ys[i] + 0.36, groups[i], fontsize=7.5, color=ps.INK_SECONDARY,
                    style="italic", va="center")
    ps.caption(fig, "Grey band: agreement within 1σ. T₀ reference: Fixsen (2009); our uncertainty "
               "is FIRAS's calibration (±1 mK), not the fit's statistical error. Peaks: Planck 2015 "
               "results I, Table 7 -- a different estimator applied to the 2015 data; through our "
               "own method, Planck's model and the data agree within errors (see 05_peaks). "
               "Parameters: Planck 2018 TT+lowE, 68% limits.")
    return _save(fig, "00_summary")


def make_all(ctx: Context) -> list:
    written = []
    written += figure_blackbody(ctx)
    written += figure_power_spectrum(ctx)
    written += figure_peaks(ctx)
    if ctx.display is not None and ctx.maps is not None:
        written += figure_sky(ctx)
        written += figure_zoom_noise(ctx)
    if ctx.maps is not None and ctx.map_compare is not None and ctx.validation is not None:
        written += figure_map_spectrum(ctx)
    if ctx.responses is not None:
        written += figure_physics(ctx)
    if ctx.fit is not None:
        written += figure_fit(ctx)
    written += figure_summary(ctx)
    return written
