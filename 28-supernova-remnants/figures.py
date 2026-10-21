"""Every figure of the project, from the cached analysis products.

Each ``fig_*`` function writes one PNG + SVG pair into ``figures/`` and
returns the paths. Images use the same grey stretch throughout (a percentile
asinh stretch); two epochs are overlaid as magenta (earlier) and green (later),
white where they agree.
"""

from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse
from scipy import ndimage

import config
import literature as L
import plotstyle as S

S.use_style()
FIG = config.FIGURE_DIR


# --------------------------------------------------------------------------
# Image helpers
# --------------------------------------------------------------------------

def stretch(img: np.ndarray, lo: float = 1.0, hi: float = 99.7, soft: float = 0.15, region=None) -> np.ndarray:
    """Percentile asinh stretch to 0-1 (computed over ``region`` if given)."""
    ref = img[region] if region is not None else img
    ref = ref[np.isfinite(ref)]
    a, b = np.percentile(ref, [lo, hi])
    x = np.clip((img - a) / (b - a), 0, 1)
    return np.arcsinh(x / soft) / np.arcsinh(1 / soft)


def blur_to(img: np.ndarray, fwhm_now: float, fwhm: float, step: float) -> np.ndarray:
    """Blur to ``fwhm`` (arcsec); holes (masked star cores) are filled from
    their surroundings first, not with the sky level."""
    import motions as M

    bad = ~np.isfinite(img)
    filled = M.inpaint(img, bad, 2.0 / step) if bad.any() else img
    sig = np.sqrt(max(fwhm**2 - fwhm_now**2, 0.0)) / 2.3548 / step
    return ndimage.gaussian_filter(filled, sig) if sig > 0 else filled


def two_epoch_rgb(early: np.ndarray, late: np.ndarray) -> np.ndarray:
    """Magenta = early, green = late, white where both are bright."""
    return np.dstack([early, late, early])


def crop(img: np.ndarray, half_px: int) -> np.ndarray:
    c = img.shape[0] // 2
    return img[c - half_px : c + half_px, c - half_px : c + half_px]


def scale_bar(ax, length_arcsec: float, label: str, *, x0=0.06, y0=0.06, color="white", extent=None):
    """Horizontal bar of ``length_arcsec`` in an axes whose x runs in arcsec
    (east to the left, so the bar is drawn from the right end)."""
    xl = ax.get_xlim()
    yl = ax.get_ylim()
    xs = xl[0] + x0 * (xl[1] - xl[0])
    ys = yl[0] + y0 * (yl[1] - yl[0])
    sign = np.sign(xl[1] - xl[0])
    ax.plot([xs, xs + sign * length_arcsec], [ys, ys], color=color, lw=2, solid_capstyle="butt")
    ax.text(xs + sign * length_arcsec / 2, ys + 0.025 * (yl[1] - yl[0]), label, color=color,
            ha="center", va="bottom", fontsize=8)


def compass(ax, *, x0=0.9, y0=0.08, size=0.08, color="white"):
    """North-up, east-left arrows in axes coordinates."""
    ax.annotate("", xy=(x0, y0 + size), xytext=(x0, y0), xycoords="axes fraction",
                arrowprops={"arrowstyle": "-|>", "color": color, "lw": 1})
    ax.annotate("", xy=(x0 - size, y0), xytext=(x0, y0), xycoords="axes fraction",
                arrowprops={"arrowstyle": "-|>", "color": color, "lw": 1})
    ax.text(x0, y0 + size + 0.015, "N", color=color, ha="center", va="bottom", fontsize=7,
            transform=ax.transAxes)
    ax.text(x0 - size - 0.015, y0, "E", color=color, ha="right", va="center", fontsize=7,
            transform=ax.transAxes)


def image_axes(ax):
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


# --------------------------------------------------------------------------
# The Galaxy's remnants
# --------------------------------------------------------------------------

def fig_galaxy_remnants(green, hist, m_points: list) -> list[str]:
    """Green's 310 remnants on the sky, size against age for the historical
    ones, and the expansion parameters we measured.

    ``m_points``: (label, age_yr, m, m_err, is_ours) tuples."""
    import remnants as R

    fig = plt.figure(figsize=(11.0, 8.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.25], hspace=0.38, wspace=0.28)
    ax = fig.add_subplot(gs[0, :])
    lw = ((green.l + 180.0) % 360.0) - 180.0  # -180..180, centre at the Galactic centre
    types = [("S", "shell", S.BLUE), ("C", "composite", S.ORANGE), ("F", "filled-centre (pulsar wind nebula)", S.AQUA),
             ("", "type unknown", S.MUTED)]
    for code, label, col in types:
        sel = green["type"] == code
        size = np.clip(green.size_arcmin[sel], 3, 300)
        ax.scatter(lw[sel], green.b[sel], s=4 + 2.2 * np.sqrt(size), color=col, alpha=0.85, lw=0.4,
                   edgecolor=S.SURFACE, label=f"{label} ({sel.sum()})", zorder=3)
    # Label offsets (deg) chosen by hand so the crowded inner Galaxy stays legible.
    offsets = {"RCW 86": (8, 9), "RX J1713.7-3946": (14, 16), "SN 1006": (-4, 7), "Crab Nebula": (-6, -14),
               "3C 58": (24, 12), "Tycho": (-10, 17), "Kepler": (24, 15), "Cas A": (2, -15), "G1.9+0.3": (-14, -16)}
    for r in hist.itertuples():
        g = green[green.name == r.green].iloc[0]
        lx = ((g.l + 180) % 360) - 180
        ox, oy = offsets.get(r.remnant, (6, 6))
        ax.annotate(f"{r.remnant}{'?' if r.disputed else ''} ({r.event.rstrip('?').replace('SN ', '')})",
                    xy=(lx, g.b), xytext=(lx + ox, g.b + oy), fontsize=7.5, color=S.INK_SECONDARY, ha="center",
                    va="center", arrowprops={"arrowstyle": "-", "color": S.MUTED, "lw": 0.6})
    ax.set_xlim(180, -180)
    ax.set_ylim(-28, 28)
    ax.set_aspect(3.0)
    ax.set_xticks(np.arange(-180, 181, 30))
    ax.set_xticklabels([f"{int(v % 360)}°" for v in np.arange(-180, 181, 30)])
    ax.set_xlabel("Galactic longitude (the Galactic centre is at 0°)")
    ax.set_ylabel("Galactic latitude\n(stretched ×3)")
    ax.set_title(f"1. The Galaxy's {len(green)} known supernova remnants (Green's catalogue, 2024 October)", pad=24)
    ax.legend(loc="lower left", ncol=4, fontsize=7, markerscale=1.3, bbox_to_anchor=(0.0, 1.0))

    # Size against age.
    ax2 = fig.add_subplot(gs[1, 0])
    t = np.logspace(1.5, 4.0, 200)
    for v, ls in ((10000, ":"), (5000, "--")):
        ax2.plot(t, R.free_radius_pc(t, v), color=S.MUTED, ls=ls, lw=1.1,
                 label=f"coasting at {v:,} km/s".replace(",", " "))
    for n, lw_, al in ((0.1, 0.9, 0.6), (1.0, 1.4, 1.0)):
        ax2.plot(t, R.sedov_radius_pc(t, 1e51, n), color=S.C_OURS, lw=lw_, alpha=al,
                 label=f"Sedov-Taylor, gas density {n:g} atom/cm$^3$")
    ax2.legend(loc="upper left", fontsize=7)
    label_pos = {"RCW 86": (1.1, 1.0, "left"), "RX J1713.7-3946": (1.1, 0.88, "left"),
                 "SN 1006": (0.9, 1.12, "right"), "Crab Nebula": (1.12, 0.8, "left"), "3C 58": (1.1, 1.14, "left"),
                 "Tycho": (0.9, 1.15, "right"), "Kepler": (1.1, 0.84, "left"), "Cas A": (0.87, 1.0, "right"),
                 "G1.9+0.3": (1.0, 0.78, "center")}
    for r in hist.itertuples():
        pwn = r.remnant in ("Crab Nebula", "3C 58")
        fc = S.SURFACE if (r.disputed or pwn) else S.C_OBSERVED
        ax2.errorbar(r.age, r.radius_pc, yerr=[[r.radius_pc - r.radius_lo], [r.radius_hi - r.radius_pc]],
                     fmt="s" if pwn else "o", ms=6, mfc=fc, mec=S.C_OBSERVED, ecolor=S.C_OBSERVED_ERR,
                     lw=1, capsize=2, zorder=4)
        fx, fy, ha = label_pos.get(r.remnant, (1.07, 1.1, "left"))
        name = r.remnant.replace("RX J1713.7-3946", "RX J1713")
        ax2.text(r.age * fx, r.radius_pc * fy, f"{name} ({r.event.replace('SN ', '')})", fontsize=7,
                 color=S.INK_SECONDARY, ha=ha, va="center", zorder=5,
                 bbox={"boxstyle": "square,pad=0.1", "fc": S.SURFACE, "ec": "none", "alpha": 0.8})
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(60, 1.2e4)
    ax2.set_ylim(0.5, 60)
    ax2.set_xticks([100, 300, 1000, 3000, 10000])
    ax2.set_xticklabels(["100", "300", "1000", "3000", "10 000"])
    ax2.set_yticks([1, 3, 10, 30])
    ax2.set_yticklabels(["1", "3", "10", "30"])
    ax2.set_xlabel("age in 2025 (years)")
    ax2.set_ylabel("radius (pc)")
    ax2.set_title("2. Remnants of supernovae seen in history")

    ax3 = fig.add_subplot(gs[1, 1])
    ax3.axhline(1.0, color=S.MUTED, lw=0.8, ls="--")
    ax3.axhline(0.4, color=S.C_OURS, lw=1.2)
    ax3.text(160, 1.02, "free expansion (coasting), m = 1", ha="left", va="bottom", fontsize=7,
             color=S.INK_SECONDARY)
    ax3.text(160, 0.42, "Sedov-Taylor blast wave, m = 0.4", ha="left", va="bottom", fontsize=7,
             color=S.INK_SECONDARY)
    for i, (label, age, m, m_err, ours) in enumerate(m_points):
        col = S.C_OURS if ours else S.C_PUBLISHED
        ax3.errorbar(age * (1.0 + 0.04 * (i % 3 - 1)), m, yerr=m_err, fmt="o" if ours else "D", ms=6, color=col,
                     mec=col, ecolor=col, lw=1, capsize=2, zorder=4)
        ax3.annotate(label, (age * (1.0 + 0.04 * (i % 3 - 1)), m), xytext=(8, -3), textcoords="offset points",
                     fontsize=7, color=S.INK_SECONDARY)
    ax3.set_xscale("log")
    ax3.set_xlim(150, 3000)
    ax3.set_xticks([200, 300, 500, 1000, 2000])
    ax3.set_xticklabels(["200", "300", "500", "1000", "2000"])
    ax3.minorticks_off()
    ax3.set_ylim(0, 1.3)
    ax3.set_xlabel("age (years)")
    ax3.set_ylabel("expansion parameter m = (speed now) / (mean speed)")
    ax3.set_title("3. Speedometer readings: how much each has slowed")
    S.caption(fig, "Sizes from Green's catalogue (geometric mean of the axes); distances and their ranges from the "
              "literature (see literature.py). Hollow: disputed identification, or (squares) a pulsar wind nebula "
              "whose size is not that of a blast wave. Sedov-Taylor curves for 10^51 erg. Panel 3: orange = this "
              "project, violet = published.", y=0.02)
    return S.save(fig, FIG / "09_galaxy_remnants")


def fig_kicks(k: dict, pulsars, homes) -> list[str]:
    """Pulsar speeds, and pulsars flying away from their remnants."""
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.3), gridspec_kw={"wspace": 0.32})
    ax = axes[0]
    galactic = pulsars[~pulsars.cluster & ~pulsars.extragalactic]
    young = galactic[galactic.young]
    rec = galactic[galactic.recycled]
    bins = np.arange(0, 1601, 50)
    ax.hist(young.v_pec, bins=bins, color=S.C_OBSERVED, alpha=0.85, label=f"young pulsars, < 3 Myr ({len(young)})",
            histtype="stepfilled", edgecolor=S.SURFACE, lw=0.5)
    ax.hist(rec.v_pec, bins=bins, color=S.BLUE, alpha=0.9, histtype="step", lw=1.4,
            label=f"recycled (millisecond) pulsars ({len(rec)})")
    v = np.linspace(0, 1600, 400)
    sig = L.HOBBS2005["sigma_1d"]
    ray = v / sig**2 * np.exp(-0.5 * v**2 / sig**2) * len(young) * 50
    ax.plot(v, ray, color=S.C_PUBLISHED, lw=1.6, label="Maxwellian kicks, σ = 265 km/s\n(Hobbs et al. 2005), seen in 2-D")
    st = k["stats"]["young"]
    ax.axvline(st["v_pec_mean"], color=S.C_OURS, lw=1.2)
    ax.set_xlabel("speed across the sky, relative to local Galactic rotation (km/s)")
    ax.set_ylabel("number of pulsars")
    ax.set_title("1. Neutron stars are born moving")
    ax.legend(fontsize=7)
    ymax = ax.get_ylim()[1]
    ax.text(st["v_pec_mean"] + 25, ymax * 0.55, f"young: mean\n{st['v_pec_mean']:.0f} ± {st['v_pec_mean_err']:.0f} km/s",
            fontsize=7.5, color=S.INK_SECONDARY)

    ax = axes[1]
    good = homes[~homes.doubtful & (homes.offset_arcmin > 2.0)]
    rest = homes[~homes.index.isin(good.index)]
    ax.scatter(good.angle, good.offset_over_radius, color=S.C_OBSERVED, s=36, zorder=3,
               label="secure, > 2′ off centre")
    ax.scatter(rest.angle, rest.offset_over_radius, facecolor=S.SURFACE, edgecolor=S.C_OBSERVED, s=30, zorder=3,
               label="doubtful, or too central")
    for r in homes.itertuples():
        ax.annotate(r.snr.split("(")[0].replace("_", " "), (r.angle, r.offset_over_radius), xytext=(4, 3),
                    textcoords="offset points", fontsize=6.5, color=S.INK_SECONDARY)
    ax.axvspan(0, 35, color=S.C_OURS, alpha=0.08, lw=0)
    ax.set_xlim(-5, 185)
    ax.set_xticks([0, 45, 90, 135, 180])
    ax.set_xlabel("angle between the pulsar's motion and the\nline from the remnant's centre to the pulsar (deg)")
    ax.set_ylabel("distance from the centre / remnant radius")
    ax.set_title("2. ...flying away from their birthplace")
    ax.legend(fontsize=7, loc="upper right")

    ax = axes[2]
    ok = homes[np.isfinite(homes.char_age_yr)]
    sec = ~ok.doubtful & (ok.offset_arcmin > 2.0)
    ax.scatter(ok.char_age_yr[sec], ok.travel_yr[sec], color=S.C_OBSERVED, s=30, zorder=3)
    ax.scatter(ok.char_age_yr[~sec], ok.travel_yr[~sec], facecolor=S.SURFACE, edgecolor=S.C_OBSERVED, s=26, zorder=3)
    show = {"S147": (6, 2), "CTB80": (6, -8), "Crab": (6, -8), "Vela": (-24, 6), "PUPPIS_A": (-30, 7),
            "G5.4-1.2": (6, 3), "3C58": (-28, -2)}
    for r in ok.itertuples():
        key = r.snr.split("(")[0]
        if key in show:
            ax.annotate(key.replace("_", " "), (r.char_age_yr, r.travel_yr), xytext=show[key],
                        textcoords="offset points", fontsize=7, color=S.INK_SECONDARY)
    lim = [5e2, 1e9]
    ax.plot(lim, lim, color=S.MUTED, ls="--", lw=1)
    ax.text(1.2e5, 2.6e5, "travel time =\nspin-down age", fontsize=7, color=S.INK_SECONDARY, ha="right")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*lim)
    ax.set_ylim(5e2, 1e6)
    ax.set_xlabel("characteristic (spin-down) age P / 2Ṗ (yr)")
    ax.set_ylabel("time to fly from the remnant's centre (yr)")
    ax.set_title("3. Two clocks for a pulsar's age")
    S.caption(fig, "ATNF Pulsar Catalogue v2.8.1: 702 pulsars with proper motions and distances. Galactic rotation "
              "and the Sun's motion removed as in Hobbs et al. (2005). Remnant centres from Green's catalogue "
              "(good to ~1'). Globular-cluster and Magellanic Cloud pulsars excluded.", y=-0.06)
    out = S.save(fig, FIG / "10_pulsar_kicks")
    return out


# --------------------------------------------------------------------------
# The Crab
# --------------------------------------------------------------------------

CRAB_HALF = 270.0  # arcsec shown either side of the pulsar


def _crab_display(crab, keys, fwhm: float = 4.5, soft: float = 300.0) -> dict:
    """Registered images blurred to one resolution and histogram-matched to
    the Pan-STARRS r image (asinh-stretched, softening ``soft`` counts), so
    plates and CCD look alike. Pan-STARRS pixels far below the sky (masked
    cores in the stack) are filled from their surroundings first."""
    from skimage.exposure import match_histograms

    import astrometry as A

    g = crab.grid
    half = int(CRAB_HALF / g.step)

    def ps1_clean(img):
        sky, sig = np.nanmedian(img[: int(60 / g.step)]), A.robust_sigma(img[: int(60 / g.step)])
        return np.where(img < sky - 8 * sig, np.nan, img)

    ref_reg = crab.regs["crab_ps1_r"]
    ref = crop(np.arcsinh(blur_to(ps1_clean(ref_reg.image), ref_reg.fwhm_arcsec, fwhm, g.step) / soft), half)
    region = crop(crab.region, half)
    out = {}
    for k in keys:
        r = crab.regs[k]
        img = ps1_clean(r.image) if "ps1" in k else r.image
        im = crop(blur_to(img, r.fwhm_arcsec, fwhm, g.step), half)
        if "ps1" in k:
            im = np.arcsinh(im / soft)
        im = match_histograms(im, ref)
        a, b = np.percentile(im[~region], 5), np.percentile(im[region], 99.95)
        out[k] = np.clip((im - a) / (b - a), 0, 1)
    return out


def fig_crab_epochs(crab, res) -> list[str]:
    """The Crab on Palomar plates in 1951 and 1991 and in Pan-STARRS (~2012)."""
    keys = ["crab_poss1_red", "crab_poss2_red", "crab_ps1_r"]
    disp = _crab_display(crab, keys)
    ext = [CRAB_HALF, -CRAB_HALF, -CRAB_HALF, CRAB_HALF]
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 10.4), gridspec_kw={"wspace": 0.04, "hspace": 0.12})
    for ax, k in zip(axes.flat[:3], keys):
        ax.imshow(disp[k], cmap=S.CMAP_IMAGE, extent=ext, origin="lower")
        image_axes(ax)
        reg = crab.regs[k]
        ax.set_title(f"{reg.plate.label.split(' (')[0]}, {reg.epoch:.1f}", fontsize=10)
        scale_bar(ax, 60.0, "1′")
        compass(ax)
    ax = axes.flat[3]
    ax.imshow(two_epoch_rgb(disp["crab_poss1_red"], disp["crab_ps1_r"]), extent=ext, origin="lower")
    image_axes(ax)
    ax.set_title("1951 (magenta) and 2012 (green) together", fontsize=10)
    scale_bar(ax, 60.0, "1′")
    fig.suptitle("Sixty years of the Crab Nebula", y=0.94)
    S.caption(fig, "Each image is registered on Gaia DR3 stars moved back to its own date (residuals 0.14\" for "
              "1951, 0.09\" for 1991, 0.02\" for Pan-STARRS), blurred to 4.5\" and matched in brightness "
              "distribution. In the overlay the stars line up (the deeper 2012 image draws them bigger, hence green "
              "halos) while the nebula's green rim is 60 years of expansion: 7% in size, 6-12\" at the edges. The "
              "1951 plate is saturated in the bright centre. Bands: H-alpha, [N II], [S II].", y=0.075)
    return S.save(fig, FIG / "01_crab_sixty_years")


def fig_crab_vectors(crab, res, v) -> list[str]:
    """Every tracked filament's motion, the expansion centre and the pulsar."""
    g = crab.grid
    half = int(CRAB_HALF / g.step)
    bg = crop(np.arcsinh(blur_to(crab.regs["crab_ps1_r"].image, 1.24, 2.0, g.step) / 30.0), half)
    bg = stretch(bg, lo=1, hi=99.8, soft=0.6)
    ext = [CRAB_HALF, -CRAB_HALF, -CRAB_HALF, CRAB_HALF]
    m = res["main"]
    fit_c = m["centre"]
    fig = plt.figure(figsize=(12.0, 7.6))
    ax = fig.add_axes([0.03, 0.08, 0.6, 0.84])
    ax.imshow(bg, cmap="Greys", extent=ext, origin="lower", vmin=-0.1, vmax=2.2)
    sub = slice(None, None, 2)
    yrs = 100.0
    ax.quiver(v.xi[sub], v.eta[sub], v.mu_xi[sub] * yrs, v.mu_eta[sub] * yrs, angles="xy", scale_units="xy",
              scale=1.0, color=S.INK, width=0.0016, headwidth=3.5, headlength=4)
    image_axes(ax)
    ax.set_xlim(CRAB_HALF, -CRAB_HALF)
    ax.set_ylim(-CRAB_HALF, CRAB_HALF)
    ax.set_title(f"1. {len(v.xi)} filament motions, 1951 → 2012, drawn as the distance each covers in 100 years",
                 fontsize=10)
    scale_bar(ax, 60.0, "1′", color=S.INK)
    compass(ax, color=S.INK)
    ax.plot(*fit_c, "+", color=S.C_OURS, ms=14, mew=2.2)
    # Inset: the centre.
    ins = fig.add_axes([0.66, 0.30, 0.32, 0.55])
    p = res["pulsar"]
    t_track = np.linspace(1054.5, 2016.0, 100)
    xs = p["pmra_mas_yr"] * (t_track - 2016.0) / 1000.0
    ys = p["pmdec_mas_yr"] * (t_track - 2016.0) / 1000.0
    ins.plot(xs, ys, color=S.C_PULSAR, lw=2)
    ins.plot(0, 0, "o", color=S.C_PULSAR, ms=7, mec=S.INK, mew=0.6)
    ins.plot(*p["position_1054"], "o", color=S.SURFACE, ms=7, mec=S.C_PULSAR, mew=1.6)
    ins.annotate("Crab pulsar today", (0, 0), xytext=(-4, 8), fontsize=8, color=S.INK,
                 arrowprops={"arrowstyle": "-", "color": S.MUTED, "lw": 0.7})
    ins.annotate("where Gaia's proper motion\nputs it in 1054", p["position_1054"], xytext=(24, 13), fontsize=8,
                 color=S.INK, arrowprops={"arrowstyle": "-", "color": S.MUTED, "lw": 0.7})
    err = np.hypot(m["centre_err_bootstrap"], m["centre_err_processing"])
    ins.add_patch(Ellipse(fit_c, 2 * err[0], 2 * err[1], fill=False, ec=S.C_OURS, lw=1.6))
    ins.plot(*fit_c, "+", color=S.C_OURS, ms=12, mew=2)
    ins.annotate(f"expansion centre,\nthis project ({m['T']:.0f} CE)", fit_c, xytext=(27, -12), fontsize=8,
                 color=S.INK, arrowprops={"arrowstyle": "-", "color": S.MUTED, "lw": 0.7})
    mc = res["martin2025_centre"]
    ins.add_patch(Circle(mc, 1.2, fill=False, ec=S.C_PUBLISHED, lw=1.4))
    ins.plot(*mc, "x", color=S.C_PUBLISHED, ms=9, mew=2)
    ins.annotate("Martin et al. 2025\n(1105.5 CE)", mc, xytext=(25, -24), fontsize=8, color=S.INK,
                 arrowprops={"arrowstyle": "-", "color": S.MUTED, "lw": 0.7})
    ins.set_xlim(30, -30)
    ins.set_ylim(-30, 30)
    ins.set_aspect("equal")
    ins.set_xlabel("arcsec east of the pulsar (2016)")
    ins.set_ylabel("arcsec north")
    ins.set_title("2. Where it all started", fontsize=10)
    S.caption(fig, f"Vectors: POSS-I red plate (1951.85) to Pan-STARRS r (2011.86), every second one of "
              f"{len(v.xi)} shown. Ellipse: 1σ (block bootstrap and processing choices combined). The pulsar's "
              f"Gaia DR3 motion, {p['mu_mas_yr']:.1f} mas/yr towards PA {p['pa_deg']:.0f}°, run back to 1054 "
              f"lands {p['offset_1054_from_centre_arcsec']:.1f}\" from our centre.", y=0.02)
    return S.save(fig, FIG / "02_crab_motions")


def fig_crab_rewind(res, v) -> list[str]:
    """Proper motion grows with distance: the slope is the age."""
    import motions as M

    m = res["main"]
    pr = res["pairs"]
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.4), gridspec_kw={"hspace": 0.38, "wspace": 0.62})
    ax = axes[0, 0]
    T, r, pa = M.per_vector_dates(v, *m["centre"])
    dx, dy = v.xi - m["centre"][0], v.eta - m["centre"][1]
    mu_r = (v.mu_xi * dx + v.mu_eta * dy) / r * 1000.0
    ax.scatter(r, mu_r, s=4, color=S.C_OBSERVED, alpha=0.35, lw=0)
    rr = np.linspace(0, 260, 50)
    tmid = pr[m["pair"]]["t_mid"]
    ax.plot(rr, rr / (tmid - m["T"]) * 1000, color=S.C_OURS, lw=2, label=f"best fit: all leave the centre in {m['T']:.0f}")
    ax.plot(rr, rr / (tmid - L.CRAB_SN_YEAR) * 1000, color=S.C_HISTORY, lw=1.6, ls="--",
            label="if they had coasted since 1054")
    ax.set_xlim(0, 250)
    ax.set_ylim(-20, 260)
    ax.set_xlabel("distance from the expansion centre (arcsec)")
    ax.set_ylabel("outward proper motion (mas/yr)")
    ax.set_title("1. The farther out, the faster")
    ax.legend(fontsize=7.5, loc="upper left")
    d = config.CRAB.distance_kpc
    sec = ax.secondary_yaxis("right", functions=(lambda x: x * config.KMS_PER_MASYR_KPC * d,
                                                 lambda y: y / config.KMS_PER_MASYR_KPC / d))
    sec.set_ylabel(f"speed at {d:g} kpc (km/s)")

    ax = axes[0, 1]
    rows = []
    for name, date, err, method, ref in L.CRAB_DATES:
        rows.append((name, date, err or 0.0, S.C_PUBLISHED, "D", method))
    for name in ("blue 1992-2012", "blue 1951-1992", "blue 1951-2012", "red 1991-2012", "red 1951-1991"):
        p = pr[name]
        rows.append((f"ours, {name}", p["T_fit"], p["T_fit_err"], S.MUTED if "blue" in name else S.C_OBSERVED,
                     "o", f"{p['n_used']} vectors"))
    rows.append((f"ours, {m['pair']} (main)", m["T"], m["T_err_total"], S.C_OURS, "o",
                 f"{pr[m['pair']]['n_used']} vectors"))
    for i, (name, date, err, col, mk, note) in enumerate(rows):
        ax.errorbar(date, i, xerr=err, fmt=mk, color=col, ms=6, capsize=2, lw=1.2)
        ax.text(1268, i, note, fontsize=6.5, va="center", ha="right", color=S.MUTED)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=7.5)
    ax.axvline(L.CRAB_SN_YEAR, color=S.C_HISTORY, lw=1.4)
    ax.text(L.CRAB_SN_YEAR + 3, len(rows) - 0.4, "guest star\n1054", fontsize=7.5, color=S.INK_SECONDARY, va="top")
    ax.set_xlim(1030, 1270)
    ax.set_xlabel("year the filaments converge (CE)")
    ax.set_title("2. Rewound to 1106, not 1054")
    ax.grid(axis="y", visible=False)

    ax = axes[1, 0]
    w = res["wedges"]
    ax.errorbar(w["pa"], w["T"], yerr=w["T_err"], fmt="o", color=S.C_OBSERVED, ms=5, capsize=2, lw=1)
    ax.axhline(m["T"], color=S.C_OURS, lw=1.5)
    ax.axhline(L.CRAB_SN_YEAR, color=S.C_HISTORY, lw=1.4)
    ax.text(362, L.CRAB_SN_YEAR + 4, "1054", fontsize=7.5, color=S.INK_SECONDARY, ha="right")
    ax.text(362, m["T"] + 4, f"all directions: {m['T']:.0f}", fontsize=7.5, color=S.INK_SECONDARY, ha="right")
    for pa, lab in ((0, "N"), (90, "E"), (180, "S"), (270, "W")):
        ax.text(pa, 1192, lab, ha="center", fontsize=9, color=S.INK_SECONDARY, fontweight="bold")
    ax.set_xlim(-5, 365)
    ax.set_ylim(1040, 1200)
    ax.set_xticks(range(0, 361, 45))
    ax.set_xlabel("position angle from the centre (deg, east of north)")
    ax.set_ylabel("convergence date in that direction (CE)")
    ax.set_title("3. Direction by direction (20 wedges)")

    ax = axes[1, 1]
    rad = res["radial"]
    for lo, hi, Tm, n in zip(rad["r_in"], rad["r_out"], rad["T_median"], rad["n"]):
        ax.plot([lo, hi], [Tm, Tm], color=S.C_OBSERVED, lw=3, solid_capstyle="butt")
        ax.text(0.5 * (lo + hi), Tm + 5, f"{n}", ha="center", fontsize=7, color=S.MUTED)
    ax.axhline(L.CRAB_SN_YEAR, color=S.C_HISTORY, lw=1.4)
    ax.axhline(m["T"], color=S.C_OURS, lw=1.5)
    ax.set_xlim(0, 260)
    ax.set_ylim(1040, 1200)
    ax.set_xlabel("distance from the centre (arcsec)")
    ax.set_ylabel("median convergence date of the vectors (CE)")
    ax.set_title("4. Inner structure looks younger")
    S.caption(fig, "Panel 2: violet = published; numbers at right = vectors used. Red pairs track the filaments' "
              "H-alpha/[N II]/[S II]; blue plates also record more of the smooth synchrotron nebula and give later "
              "dates. Panel 4: per-vector dates (distance / outward speed); numbers = vectors per ring.", y=0.02)
    return S.save(fig, FIG / "03_crab_rewind")


def fig_crab_stills(crab, res, v) -> list[str]:
    """The filaments run backwards in time, and how tightly they huddle."""
    import motions as M

    m = res["main"]
    p = res["pulsar"]
    tmid = v.tmid
    years = [2012, 1800, 1500, 1300, round(m["T"])]
    fig = plt.figure(figsize=(14.0, 3.3))
    g = crab.grid
    half = int(CRAB_HALF / g.step)
    reg = crop(crab.region, half)
    for i, yr in enumerate(years):
        ax = fig.add_axes([0.005 + i * 0.148, 0.08, 0.14, 0.72])
        x = v.xi - v.mu_xi * (tmid - yr)
        y = v.eta - v.mu_eta * (tmid - yr)
        ax.contour(np.linspace(CRAB_HALF, -CRAB_HALF, reg.shape[1]), np.linspace(-CRAB_HALF, CRAB_HALF, reg.shape[0]),
                   reg.astype(float), levels=[0.5], colors=[S.C_CONTEXT], linewidths=1)
        ax.scatter(x, y, s=1.2, color=S.C_OBSERVED, lw=0)
        px = p["pmra_mas_yr"] * (yr - 2016.0) / 1000.0
        py = p["pmdec_mas_yr"] * (yr - 2016.0) / 1000.0
        ax.plot(px, py, "o", color=S.C_PULSAR, ms=4, mec=S.INK, mew=0.4)
        ax.plot(*m["centre"], "+", color=S.C_OURS, ms=8, mew=1.6)
        ax.set_xlim(CRAB_HALF, -CRAB_HALF)
        ax.set_ylim(-CRAB_HALF, CRAB_HALF)
        ax.set_aspect("equal")
        image_axes(ax)
        ax.set_title(f"{yr}" + (": they meet" if yr == round(m["T"]) else ""), fontsize=10)
    ax = fig.add_axes([0.77, 0.2, 0.21, 0.6])
    t_ca, grid_T, spread = M.closest_approach(v, T_grid=np.arange(950.0, 1400.0, 2.0))
    ax.plot(grid_T, spread, color=S.C_OBSERVED, lw=1.6)
    ax.axvline(m["T"], color=S.C_OURS, lw=1.4)
    ax.axvline(t_ca, color=S.C_OBSERVED, lw=1, ls="--")
    ax.axvline(L.CRAB_SN_YEAR, color=S.C_HISTORY, lw=1.4)
    ymax = spread.max()
    ax.text(m["T"] - 4, ymax * 0.92, f"fit {m['T']:.0f}", ha="right", fontsize=7.5, color=S.INK_SECONDARY)
    ax.text(t_ca + 4, ymax * 0.92, f"tightest\n{t_ca:.0f}", ha="left", fontsize=7.5, color=S.INK_SECONDARY)
    ax.text(L.CRAB_SN_YEAR + 4, ymax * 0.45, "1054", fontsize=7.5, color=S.INK_SECONDARY)
    ax.set_xlabel("year (CE)")
    ax.set_ylabel("rms spread (arcsec)")
    ax.set_title("How tightly they huddle", fontsize=10)
    fig.suptitle("The Crab's filaments run backwards at their present speeds", y=1.0)
    S.caption(fig, f"Dots: the {len(v.xi)} tracked filaments, each moved back along its own measured motion. Grey "
              "outline: the nebula today. Aqua: the pulsar (Gaia DR3). Orange cross: the fitted centre. Right: "
              "errors in the motions (~13 mas/yr, x900 yr) keep the cloud ~25\" wide and push its tightest moment "
              "later than the fit, the bias of the 'closest approach' method.", y=-0.04)
    return S.save(fig, FIG / "04_crab_run_backwards")


def gif_crab(crab, path=None) -> str:
    """Blink: 1951 -> 1991 -> 2012, same stretch, year stamped."""
    import imageio.v3 as iio
    from PIL import Image, ImageDraw, ImageFont

    keys = ["crab_poss1_red", "crab_poss2_red", "crab_ps1_r"]
    disp = _crab_display(crab, keys)
    frames = []
    for k in keys:
        # Arrays have east at column 0 and north at the last row: flip rows only.
        im = (np.clip(disp[k][::-1], 0, 1) * 255).astype(np.uint8)
        im = Image.fromarray(im).resize((540, 540), Image.LANCZOS).convert("RGB")
        d = ImageDraw.Draw(im)
        try:
            font = ImageFont.truetype("DejaVuSans-Bold.ttf", 26)
        except OSError:
            font = ImageFont.load_default()
        d.text((14, 10), f"{crab.regs[k].epoch:.0f}", fill=(255, 255, 255), font=font)
        frames.append(np.asarray(im))
    path = path or (FIG / "crab_blink.gif")
    iio.imwrite(path, frames, duration=[900, 900, 900], loop=0)
    return str(path)


# --------------------------------------------------------------------------
# Cas A in X-rays
# --------------------------------------------------------------------------

def _xray_norm(counts: np.ndarray, sigma: float = 0.8) -> np.ndarray:
    """Counts per unit total, lightly smoothed and log-stretched (0-1)."""
    im = ndimage.gaussian_filter(counts, sigma) / counts.sum() * 1e6
    return np.clip(np.log10(np.clip(im, 1e-2, None)) / np.log10(np.percentile(im, 99.9)), 0, 1)


def fig_casa_xray(frame, epochs, summ) -> list[str]:
    import casa_xray as X

    pix = config.CHANDRA_PIXEL_ARCSEC
    cx, cy = frame.radec_to_box(X.CENTRE_RA, X.CENTRE_DEC)
    cx, cy = float(cx), float(cy)
    half = int(200 / pix)
    sl = (slice(int(cy) - half, int(cy) + half), slice(int(cx) - half, int(cx) + half))
    ext = [(cx - (int(cx) - half)) * pix, -(int(cx) + half - cx) * pix,
           -(cy - (int(cy) - half)) * pix, (int(cy) + half - cy) * pix]
    ref = next(e for e in epochs if e.label == config.CHANDRA_REFERENCE_EPOCH)
    first, last = epochs[0], epochs[-1]
    fig = plt.figure(figsize=(13.0, 5.2))
    ax = fig.add_axes([0.02, 0.1, 0.3, 0.8])
    ax.imshow(_xray_norm(ref.counts)[sl], cmap=S.CMAP_XRAY, extent=ext, origin="lower")
    for r0, ls, lab in ((X.FORWARD_SHOCK, "-", "forward-shock annulus"), (X.EJECTA_RING, ":", "ejecta-ring annulus")):
        for rr in r0:
            ax.add_patch(Circle((0, 0), rr * 60, fill=False, ec="white", lw=0.8, ls=ls, alpha=0.8))
    for pa in np.arange(0, 360, X.SECTOR_DEG):
        t = np.radians(pa)
        ax.plot([2.4 * 60 * np.sin(t), 3.0 * 60 * np.sin(t)], [2.4 * 60 * np.cos(t), 3.0 * 60 * np.cos(t)],
                color="white", lw=0.5, alpha=0.6)
    ax.plot(0, 0, "+", color="white", ms=10, mew=1.5)
    image_axes(ax)
    ax.set_xlim(200, -200)
    ax.set_ylim(-200, 200)
    scale_bar(ax, 60.0, "1′")
    compass(ax)
    ax.set_title(f"1. Chandra, 2004 ({ref.exposure / 1e6:.2f} Ms, all energies)", fontsize=10)
    ax.text(0.03, 0.95, "forward shock: 2.4-3.0′\nejecta ring: 1.67-2.18′", transform=ax.transAxes, color="white",
            fontsize=7.5, va="top")

    ax = fig.add_axes([0.34, 0.1, 0.3, 0.8])
    a = _xray_norm(first.counts)[sl]
    b = _xray_norm(last.counts)[sl]
    ax.imshow(two_epoch_rgb(a, b), extent=ext, origin="lower")
    image_axes(ax)
    ax.set_xlim(200, -200)
    ax.set_ylim(-200, 200)
    scale_bar(ax, 60.0, "1′")
    ax.set_title(f"2. {first.year:.0f} (magenta) and {last.year:.0f} (green)", fontsize=10)

    ax = fig.add_axes([0.71, 0.14, 0.27, 0.72])
    yy, xx = np.mgrid[: frame.shape[0], : frame.shape[1]]
    dx, dy = xx - cx, yy - cy
    r = np.hypot(dx, dy) * pix
    pa = np.degrees(np.arctan2(-dx, dy)) % 360.0
    sel = (pa >= 135) & (pa < 160)
    edges = np.arange(120, 200, 1.0)
    mid = 0.5 * (edges[:-1] + edges[1:])
    half_pts = []
    for ep, col, lab in ((first, S.MAGENTA, f"{first.year:.0f}"), (last, S.GREEN, f"{last.year:.0f}")):
        prof = np.array([ep.counts[sel & (r >= lo) & (r < hi)].mean() for lo, hi in zip(edges[:-1], edges[1:])])
        prof = ndimage.uniform_filter1d(prof, 3)
        rim = (mid > 140) & (mid < 175)
        bg = np.median(prof[mid > 185])
        peak = prof[rim].max()
        prof = (prof - bg) / (peak - bg)
        ax.plot(mid, prof, color=col, lw=1.8, label=lab)
        i_pk = np.flatnonzero(rim)[np.argmax(prof[rim])]
        beyond = np.flatnonzero((np.arange(len(prof)) > i_pk) & (prof < 0.5))
        half_pts.append(np.interp(0.5, prof[beyond[0]:beyond[0] - 2:-1], mid[beyond[0]:beyond[0] - 2:-1]))
    ax.annotate("", xy=(half_pts[1], 0.5), xytext=(half_pts[0], 0.5),
                arrowprops={"arrowstyle": "->", "color": S.INK, "lw": 1.3})
    ax.text(half_pts[1] + 1.5, 0.53, f"{half_pts[1] - half_pts[0]:.1f}\" in {last.year - first.year:.0f} yr",
            fontsize=7.5, color=S.INK)
    ax.set_ylim(-0.1, 1.15)
    ax.set_xlabel("distance from the expansion centre (arcsec)")
    ax.set_ylabel("brightness, scaled to the rim's peak")
    ax.set_title("3. The south-east rim moves out (PA 135-160°)", fontsize=10)
    ax.legend(fontsize=7.5)
    S.caption(fig, "Chandra ACIS-S 'centre images' (all energies, 0.492\" pixels) from the Chandra Data Archive. "
              "The cross marks the optical expansion centre of Thorstensen et al. (2001). In the overlay, green "
              "fringes on the outside of the rim are 19 years of expansion.", y=0.02)
    return S.save(fig, FIG / "05_casa_xray")


def fig_casa_expansion(summ) -> list[str]:
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6), gridspec_kw={"wspace": 0.55})
    ax = axes[0]
    reg = summ["regions"]
    t_ref = summ["t_ref"]
    for key, sub, col, lab in (("forward_shock", "all", S.C_OURS, "forward shock (2.4-3.0′)"),
                               ("ejecta_ring", "all", S.C_OBSERVED, "ejecta ring (1.67-2.18′)")):
        rec = reg[key][sub]
        t = np.array(rec["t"])
        s = (np.array(rec["s"]) - 1) * 100
        ax.plot(t, s, "o", color=col, ms=5)
        tt = np.linspace(1999, 2020, 10)
        ax.plot(tt, rec["rate_pct_yr"] * (tt - t_ref), color=col, lw=1.4,
                label=f"{lab}: {rec['rate_pct_yr']:.3f} %/yr")
    ax.axhline(0, color=S.MUTED, lw=0.6)
    ax.axvline(t_ref, color=S.MUTED, lw=0.6, ls=":")
    ax.text(t_ref + 0.3, -1.0, "reference:\n2004 deep image", fontsize=7, color=S.INK_SECONDARY)
    ax.set_xlabel("year")
    ax.set_ylabel("size relative to 2004 (%)")
    ax.set_title("1. Fifteen epochs, one straight line")
    ax.legend(fontsize=7.5, loc="upper left")

    age = summ["age_at_ref"]
    for ax, key, pub, title in ((axes[1], "forward_shock", L.VINK2022_FS, "2. Forward shock, sector by sector"),
                                (axes[2], "ejecta_ring", L.VINK2022_RS, "3. Ejecta ring, sector by sector")):
        pas = pub["pa"]
        ours = [reg[key][str(p)]["rate_pct_yr"] for p in pas]
        errs = [reg[key][str(p)]["rate_err"] for p in pas]
        ax.errorbar(pas, ours, yerr=errs, fmt="o-", color=S.C_OURS, ms=5, lw=1.2, capsize=2,
                    label="this project (all energies)")
        lab = "Vink et al. 2022 (4.2-6 keV)" + ("" if key == "forward_shock" else ", reverse-shock region")
        ax.plot(pas, pub["rate"], "D--", color=S.C_PUBLISHED, ms=4.5, lw=1.0, label=lab)
        ax.axhline(0, color=S.MUTED, lw=0.6)
        ax.set_xticks(range(0, 361, 45))
        ax.set_xlim(-5, 365)
        ax.set_xlabel("position angle (deg, east of north)")
        ax.set_ylabel("expansion rate (% per year)")
        ax.set_title(title)
        ax.legend(fontsize=7, loc="lower left")
        sec = ax.secondary_yaxis("right", functions=(lambda x: x / 100 * age, lambda m: m * 100 / age))
        sec.set_ylabel(f"expansion parameter m (age {age:.0f} yr)")
        for pa, lab_ in ((0, "N"), (90, "E"), (180, "S"), (270, "W")):
            ax.text(pa, 0.315, lab_, ha="center", fontsize=8.5, color=S.INK_SECONDARY, fontweight="bold")
        ax.set_ylim(-0.13, 0.33)
    axes[1].axvspan(200, 280, color=S.C_CONTEXT, alpha=0.4, lw=0)
    axes[1].text(240, 0.035, "ejecta inside\nthe annulus", ha="center", fontsize=7, color=S.INK_SECONDARY)
    S.caption(fig, "Expansion rate = fractional growth per year about the optical expansion centre; m = rate × age "
              "(free expansion m = 1, Sedov-Taylor m = 0.4). Errors are statistical, scaled by the scatter about "
              "each line; pointing differences between epochs are fitted per epoch.", y=-0.04)
    return S.save(fig, FIG / "06_casa_expansion")


def fig_casa_neutron_star(frame, epochs, ns) -> list[str]:
    import casa_xray as X

    pix = config.CHANDRA_PIXEL_ARCSEC
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5), gridspec_kw={"wspace": 0.33})
    ax = axes[0]
    ref = next(e for e in epochs if e.label == config.CHANDRA_REFERENCE_EPOCH)
    cx, cy = frame.radec_to_box(X.CENTRE_RA, X.CENTRE_DEC)
    cx, cy = float(cx), float(cy)
    half = int(14 / pix)
    sl = (slice(int(cy) - half, int(cy) + half), slice(int(cx) - half, int(cx) + half))
    ext = [(cx - (int(cx) - half) + 0.5) * pix, -(int(cx) + half - cx - 0.5) * pix,
           -(cy - (int(cy) - half) + 0.5) * pix, (int(cy) + half - cy - 0.5) * pix]
    ax.imshow(np.sqrt(ndimage.gaussian_filter(ref.counts, 0.6)[sl]), cmap=S.CMAP_XRAY, extent=ext,
              origin="lower", interpolation="nearest")
    off = ns["offset"]
    ax.add_patch(Circle((0, 0), 0.4, fill=False, ec="white", lw=1.2))
    ax.plot(0, 0, "+", color="white", ms=10, mew=1.4)
    ax.annotate("", xy=(off["east_arcsec"], off["north_arcsec"]), xytext=(0, 0),
                arrowprops={"arrowstyle": "->", "color": "white", "lw": 1.4})
    ax.text(13, 5.5, "expansion centre\n(fastest optical knots,\nThorstensen et al. 2001)", color="white",
            fontsize=7.5, ha="left")
    ax.text(13, -9.0, f"neutron star: {off['arcsec']:.1f}\" away\nat PA {off['pa_deg']:.0f}° → {off['implied_v_kms']:.0f} km/s"
            f"\nif it left the centre\n{off['age_yr']:.0f} years ago", color="white", fontsize=7.5, va="top", ha="left")
    ax.set_xlim(14, -14)
    ax.set_ylim(-14, 14)
    ax.set_xlabel("arcsec east of the centre")
    ax.set_ylabel("arcsec north")
    ax.grid(False)
    ax.set_title("1. A neutron star off-centre", fontsize=10)

    ax = axes[1]
    t = np.array([e["year"] for e in ns["epochs"]])
    for tag, mk, col, lab in (("archive", "o", S.MUTED, "archive astrometry"),
                              ("corrected", "s", S.C_OURS, "aligned on the remnant")):
        ky = "y" if tag == "archive" else "y_corr"
        Y_ = np.array([e[ky] for e in ns["epochs"]])
        n_ = (Y_ - Y_.mean()) * pix * 1000
        ax.plot(t, n_, mk, color=col, ms=4.5, label=f"north, {lab}")
        rec = ns[tag]
        ax.plot(t, rec["mu_north_mas_yr"] * (t - t.mean()), color=col, lw=1.2)
    ax.set_xlabel("year")
    ax.set_ylabel("northward position (mas, about the mean)")
    ax.set_title("2. Its position, 2000-2019", fontsize=10)
    ax.legend(fontsize=7)

    ax = axes[2]
    est = [
        ("offset from centre\n(this project)", off["implied_v_kms"], off["pa_deg"], 0, 0, S.C_OURS),
        ("aligned on remnant\n(this project; biased)", ns["corrected"]["v_kms"], ns["corrected"]["pa_deg"],
         ns["corrected"]["v_err_kms"], ns["corrected"]["pa_err_deg"], S.C_OURS),
        ("archive astrometry\n(this project)", ns["archive"]["v_kms"], ns["archive"]["pa_deg"],
         ns["archive"]["v_err_kms"], ns["archive"]["pa_err_deg"], S.MUTED),
        ("Holland-Ashford+ 2024,\nmethod 1", L.HOLLAND_ASHFORD2024["method1_kms"], L.HOLLAND_ASHFORD2024["method1_pa"],
         L.HOLLAND_ASHFORD2024["method1_err"], L.HOLLAND_ASHFORD2024["method1_pa_err"], S.C_PUBLISHED),
        ("Holland-Ashford+ 2024,\nmethod 2 (Gaia-registered)", L.HOLLAND_ASHFORD2024["method2_kms"],
         L.HOLLAND_ASHFORD2024["method2_pa"], L.HOLLAND_ASHFORD2024["method2_err"],
         L.HOLLAND_ASHFORD2024["method2_pa_err"], S.C_PUBLISHED),
        ("Fesen et al. 2006\n(offset)", L.FESEN2006["cco_speed_kms"], L.FESEN2006["cco_pa_deg"], 0,
         L.FESEN2006["cco_pa_err"], S.C_PUBLISHED),
    ]
    for i, (lab, v, pa, ve, pae, col) in enumerate(est):
        y = len(est) - 1 - i
        ax.errorbar(v, y, xerr=ve or None, fmt="o", color=col, ms=6, capsize=2)
        ax.text(1180, y, f"PA {pa:.0f}°" + (f" ± {pae:.0f}" if pae else ""), fontsize=7, va="center", ha="right",
                color=S.INK_SECONDARY)
    ax.set_yticks(range(len(est)))
    ax.set_yticklabels([e[0] for e in est[::-1]], fontsize=7)
    ax.set_xlim(0, 1200)
    ax.set_xlabel("speed across the sky at 3.4 kpc (km/s)")
    ax.set_title("3. How fast? Six estimates", fontsize=10)
    ax.grid(axis="y", visible=False)
    S.caption(fig, "The 'aligned' track removes Chandra's pointing scatter by fitting each epoch's shift to the "
              "remnant, but Cas A expands faster in the east than in the west, and that lopsidedness leaks into the "
              "shift: the neutron star then appears to drift west. A clean measurement needs X-ray sources with Gaia "
              "counterparts, which the 8′ centre images lack.", y=-0.06)
    return S.save(fig, FIG / "07_casa_neutron_star")


def gif_casa(frame, epochs, path=None) -> str:
    """All fifteen Chandra epochs of Cas A, 2000-2019, at a fixed stretch."""
    import imageio.v3 as iio
    from PIL import Image, ImageDraw, ImageFont
    import casa_xray as X

    cx, cy = frame.radec_to_box(X.CENTRE_RA, X.CENTRE_DEC)
    half = int(200 / config.CHANDRA_PIXEL_ARCSEC)
    sl = (slice(int(cy) - half, int(cy) + half), slice(int(cx) - half, int(cx) + half))
    cmap = plt.get_cmap(S.CMAP_XRAY)
    frames = []
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", 24)
    except OSError:
        font = ImageFont.load_default()
    for ep in epochs:
        im = _xray_norm(ep.counts)[sl]
        rgb = (cmap(im)[..., :3] * 255).astype(np.uint8)[::-1]  # rows flipped: north up; east already left
        rgb = np.where(ep.footprint[sl][::-1, :, None], rgb, 0).astype(np.uint8)
        pil = Image.fromarray(rgb).resize((560, 560), Image.LANCZOS)
        d = ImageDraw.Draw(pil)
        d.text((14, 10), f"{ep.year:.1f}", fill=(255, 255, 255), font=font)
        frames.append(np.asarray(pil))
    path = path or (FIG / "casa_xray_movie.gif")
    iio.imwrite(path, frames + frames[-1:] * 2, duration=350, loop=0)
    return str(path)


# --------------------------------------------------------------------------
# Cas A in visible light: an honest failure, with the Crab as control
# --------------------------------------------------------------------------

def fig_casa_optical(opt: dict) -> list[str]:
    import plates

    regs = plates.register_all("casa")
    grid = plates.grid_for(config.CASA)
    half = int(190 / grid.step)
    keys = ["casa_poss1_red", "casa_poss2_red", "casa_ps1_r"]
    fig = plt.figure(figsize=(13.5, 7.6))
    for i, k in enumerate(keys):
        ax = fig.add_axes([0.01 + i * 0.205, 0.5, 0.2, 0.42])
        r = regs[k]
        im = crop(blur_to(r.image, r.fwhm_arcsec, 4.5, grid.step), half)
        if "ps1" in k:
            im = np.arcsinh(im / 30.0)
        xi, eta = grid.coords
        rr = crop(np.hypot(xi, eta), half)
        shell = (rr > 55) & (rr < 175)
        a, b = np.percentile(im, 20), np.percentile(im[shell], 99.7)
        ax.imshow(np.clip((im - a) / (b - a), 0, 1) ** 0.7, cmap=S.CMAP_IMAGE, origin="lower",
                  extent=[190, -190, -190, 190])
        ax.set_xlim(190, -190)
        ax.set_ylim(-190, 190)
        ax.add_patch(Circle((0, 0), 55, fill=False, ec=S.ORANGE, lw=0.6, ls="--"))
        ax.add_patch(Circle((0, 0), 175, fill=False, ec=S.ORANGE, lw=0.6, ls="--"))
        image_axes(ax)
        ax.set_title(f"{r.plate.label.split(' (')[0]}, {r.epoch:.1f}", fontsize=10)
        scale_bar(ax, 60.0, "1′")
    ax = fig.add_axes([0.66, 0.53, 0.32, 0.36])
    cc = opt["crab_control"]
    ax.plot(cc["s"], cc["ncc"], color=S.C_OBSERVED, lw=1.8)
    ax.axvline(cc["s_expected"], color=S.C_OURS, lw=1.4)
    ax.text(cc["s_expected"] + 0.003, min(cc["ncc"]) + 0.05, f"predicted by the\nfilament vectors: {cc['s_expected']:.3f}",
            fontsize=7.5, color=S.INK_SECONDARY)
    ax.set_xlabel("stretch applied to the 1951 image")
    ax.set_ylabel("match with 2012 (correlation)")
    ax.set_title("1. Control: the Crab, 1951 → 2012", fontsize=10)

    ax = fig.add_axes([0.07, 0.08, 0.86, 0.32])
    styles = {"red 1954-1991": (S.C_OBSERVED, "-"), "red 1991-2012": (S.BLUE, "-"),
              "red 1954-2012": (S.MUTED, "-"), "near-IR 1991-2012": (S.BLUE, "--")}
    for name, p in opt["pairs"].items():
        col, ls = styles.get(name, (S.MUTED, ":"))
        ax.plot(p["scan_s"], p["scan_ncc"], color=col, ls=ls, lw=1.6, label=f"Cas A, {name}")
        ax.axvline(p["s_expected"], color=col, ls=":", lw=1.2)
    ax.set_xlabel("stretch applied to the earlier image (dotted: what a shell coasting since 1671 would show)")
    ax.set_ylabel("match (correlation)")
    ax.set_title("2. Cas A: every pair matches best with no stretch at all", fontsize=10)
    ax.legend(fontsize=7.5, ncol=2, loc="upper right")
    S.caption(fig, "Each image registered on Gaia stars at its own date; stars masked out of every comparison. Dashed "
              "circles: the shell annulus compared (55-175\" from the optical expansion centre). On survey plates "
              "Cas A's knots are faint and short-lived, so what the epochs share is structure that does not move "
              "(slow circumstellar flocculi, foreground gas). The same code recovers the Crab's expansion.", y=0.0)
    return S.save(fig, FIG / "08_casa_optical")
