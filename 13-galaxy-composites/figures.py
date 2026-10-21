"""Figure generation.

Every figure is drawn from the analysis results and images, never by re-running
the measurements, so what is plotted is exactly what was measured.

Colour in the *charts* follows one rule, so it means the same thing
everywhere: **ink** is what was measured, **orange** is our own measurement or
choice, **violet** a published value we compare against, and the SDSS filters
are drawn in the display channel each one feeds (g blue, r green, i red), always
with a direct label. The colour-index maps use a diverging scale centred on the
colour of the Sun: blue is bluer than sunlight, red is redder, grey is
sunlight-coloured. Colour in the *pictures* is the subject of the project, and
each picture's caption says exactly how it was made.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm, TwoSlopeNorm
from matplotlib.patches import Circle, Rectangle

import color_physics as cp
import composite as cm
import plotstyle as ps
from config import (
    BANDS,
    FIGURE_DIR,
    GALLERY_DIR,
    LUPTON_Q,
    LUPTON_SATURATION,
    MESOPIC_LOWER_CD_M2,
    MESOPIC_UPPER_CD_M2,
)

ps.use_style()

PROVENANCE = ("Images: SDSS DR8 g, r, i (Aihara et al. 2011), via the Legacy Surveys viewer "
              "(Dey et al. 2019). Colour science: CIE 1931 colour-matching functions; the Sun's "
              "spectrum from HST CALSPEC.")

DASHED = (0, (5, 2.5))
DOTTED = (0, (1, 1.8))

#: The logarithmic stretch's black level, nanomaggies/arcsec^2 (mu_r ~ 26.7):
#: below the pixel noise, as a log stretch must be to show the faint outskirts.
LOG_FLOOR = 0.02

#: The diverging scale for colour indices: blue arm, sunlight-grey midpoint, red arm.
CMAP_COLOR_INDEX = LinearSegmentedColormap.from_list(
    "bluer_redder", ["#0d366b", "#2a78d6", "#86b6ef", "#f0efec", "#f0a19a", "#e34948", "#8f2321"])


# --------------------------------------------------------------------------
# Building blocks
# --------------------------------------------------------------------------

def _title(ax, title: str, subtitle: str | None = None, *, fontsize=None) -> None:
    lines = subtitle.count("\n") + 1 if subtitle else 0
    ax.set_title(title, pad=8 + 10.8 * lines, **({"fontsize": fontsize} if fontsize else {}))
    if subtitle:
        ax.annotate(subtitle, xy=(0, 1), xycoords="axes fraction", xytext=(0, 5),
                    textcoords="offset points", ha="left", va="bottom", fontsize=8.2,
                    color=ps.INK_SECONDARY, linespacing=1.3)


def _label(ax, text, xy, xytext, *, ha="left", va="center", leader=True, fontsize=8,
           textcoords="data", color=None):
    ax.annotate(text, xy=xy, xytext=xytext, textcoords=textcoords, ha=ha, va=va,
                fontsize=fontsize, color=color or ps.INK_SECONDARY,
                arrowprops=({"arrowstyle": "-", "color": ps.MUTED, "linewidth": 0.6,
                             "shrinkA": 1, "shrinkB": 3} if leader else None))


def _downsample(rgb: np.ndarray, max_px: int) -> np.ndarray:
    """Block-average an image so figure files stay small (the full-resolution
    pictures are saved separately in figures/gallery/)."""
    n = rgb.shape[0]
    k = int(np.ceil(n / max_px))
    if k <= 1:
        return rgb
    h, w = (rgb.shape[0] // k) * k, (rgb.shape[1] // k) * k
    x = rgb[:h, :w]
    return x.reshape(h // k, k, w // k, k, *x.shape[2:]).mean(axis=(1, 3))


def _scalebar_arcmin(fov_arcmin: float) -> float:
    for s in (0.5, 1, 2, 5, 10):
        if s >= fov_arcmin / 8:
            return s
    return 10


def _image(ax, rgb, *, fov_arcmin: float, max_px: int = 700, scalebar: bool = True,
           cmap=None, norm=None, label: str | None = None):
    """Show an image (FITS row order) with north up, a scale bar, and no chart chrome."""
    img = _downsample(rgb, max_px)
    n = img.shape[0]
    ax.imshow(img, origin="lower", cmap=cmap, norm=norm, interpolation="lanczos",
              extent=(0, n, 0, n))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_facecolor(ps.IMAGE_BG)
    if scalebar:
        sb = _scalebar_arcmin(fov_arcmin)
        length = n * sb / fov_arcmin
        x0, y0 = 0.05 * n, 0.05 * n
        ax.plot([x0, x0 + length], [y0, y0], color="white", lw=1.6, solid_capstyle="butt")
        ax.text(x0 + length / 2, y0 + 0.02 * n, f"{sb:g}′", color="white", fontsize=7,
                ha="center", va="bottom")
    if label:
        ax.text(0.03, 0.97, label, transform=ax.transAxes, color="white", fontsize=8,
                ha="left", va="top", fontweight="bold")
    return n


def _caption(fig, text, y=-0.01):
    fig.text(0.5, y, text, ha="center", va="top", fontsize=7.5, color=ps.MUTED, wrap=True)


def _save(fig, name: str, *, svg_dpi: int | None = None) -> list[str]:
    """PNG at full resolution; the SVG keeps text and lines as vectors, and for
    picture-heavy figures embeds the pictures at a lower resolution
    (``svg_dpi``) so the file stays a few MB."""
    import matplotlib.pyplot as plt_
    png = FIGURE_DIR / f"{name}.png"
    fig.savefig(png, format="png")
    svg = FIGURE_DIR / f"{name}.svg"
    fig.savefig(svg, format="svg", **({"dpi": svg_dpi} if svg_dpi else {}))
    plt_.close(fig)
    return [str(png), str(svg)]


def _band_label(b: str) -> str:
    return f"{b}  ({cp.pivot_wavelength(f'sdss_{b}') / 10:.0f} nm)"


# --------------------------------------------------------------------------
# 1. How a colour picture is made
# --------------------------------------------------------------------------

def figure_how_made(art: dict, ref: dict) -> list[str]:
    fr = art["frame"]
    g = fr.galaxy
    w = ref["white_balance"]["ours"]["weights"]
    soft = cm.softening_for(LUPTON_Q, LUPTON_SATURATION)
    fig, axes = plt.subplots(1, 4, figsize=(13.2, 3.55), gridspec_kw={"wspace": 0.18})
    channel_name = {"g": "blue", "r": "green", "i": "red"}
    for ax, b in zip(axes[:3], BANDS):
        grey = cm.asinh_curve(np.clip(fr.sb[b] * w[b], 0, None), soft, LUPTON_SATURATION)
        _image(ax, np.clip(grey, 0, 1), fov_arcmin=g.fov_arcmin, cmap="gray",
               norm=plt.Normalize(0, 1), scalebar=(b == "g"))
        _title(ax, f"{b} filter, {cp.pivot_wavelength(f'sdss_{b}') / 10:.0f} nm",
               f"one greyscale image → the {channel_name[b]} channel", fontsize=10)
    _image(axes[3], art["composite"]["standard"], fov_arcmin=g.fov_arcmin, scalebar=False)
    _title(axes[3], "Combined", "i → red, r → green, g → blue", fontsize=10)
    for x, sym in ((0.318, "+"), (0.521, "+"), (0.724, "=")):
        fig.text(x, 0.48, sym, fontsize=18, color=ps.INK_SECONDARY, ha="center", va="center")
    fig.suptitle(f"A telescope sees in greyscale: {g.name}, the {g.common}, "
                 "through three filters", x=0.125, ha="left", y=1.03)
    _caption(fig, "Each panel is a calibrated SDSS image (brightness in nanomaggies per square "
             "arcsecond), stretched with the same arcsinh curve. The galaxy's old stars are "
             "brighter in the red i filter, its young ones in blue g. Colour appears only when the "
             "three are assigned to a screen's red, green and blue. " + PROVENANCE, y=0.07)
    return _save(fig, "01_how_a_picture_is_made", svg_dpi=90)


# --------------------------------------------------------------------------
# 2. The stretch
# --------------------------------------------------------------------------

def figure_stretch(art: dict, ref: dict, profile_display: dict) -> list[str]:
    fr = art["frame"]
    g = fr.galaxy
    w = ref["white_balance"]["ours"]["weights"]
    R, G, B = cm.channels(fr.sb, w)
    soft = cm.softening_for(LUPTON_Q, LUPTON_SATURATION)
    panels = [
        ("Linear", "brightness proportional to light:\nonly the cores show",
         cm.linear_rgb(R, G, B, saturation=LUPTON_SATURATION)),
        ("Logarithmic", "every factor of 10 gets equal space:\nthe sky noise turns grey, colours wash out",
         cm.log_rgb(R, G, B, floor=LOG_FLOOR, saturation=LUPTON_SATURATION)),
        ("Arcsinh, channel by channel", "bright parts burn out towards white:\ncolour lost where the galaxy is brightest",
         cm.per_channel_asinh(R, G, B, softening=soft, saturation=LUPTON_SATURATION)),
        ("Arcsinh on total light (Lupton et al. 2004)", "every channel scaled by the same factor:\ncolour kept everywhere",
         art["composite"]["standard"]),
    ]
    fig = plt.figure(figsize=(13.0, 8.9))
    gs = fig.add_gridspec(2, 3, wspace=0.16, hspace=0.42)
    spots = [gs[0, 0], gs[0, 1], gs[0, 2], gs[1, 0]]
    for spot, (title, sub, rgb) in zip(spots, panels):
        ax = fig.add_subplot(spot)
        _image(ax, rgb, fov_arcmin=g.fov_arcmin, max_px=520)
        _title(ax, title, sub, fontsize=10)

    # transfer curves
    ax = fig.add_subplot(gs[1, 1])
    mu = np.linspace(26.5, 14.5, 400)
    inten = cp.mag_to_nmgy(mu)
    logc = lambda i: np.clip((np.log10(np.clip(i, LOG_FLOOR, None)) - np.log10(LOG_FLOOR))
                             / (np.log10(LUPTON_SATURATION) - np.log10(LOG_FLOOR)), 0, 1)
    asin = lambda i: np.clip(cm.asinh_curve(i, soft, LUPTON_SATURATION), 0, 1)
    curves = [("linear", ps.BLUE, np.clip(inten / LUPTON_SATURATION, 0, 1)),
              ("arcsinh", ps.ORANGE, asin(inten)),
              ("log", ps.AQUA, logc(inten))]
    for name, color, y in curves:
        ax.plot(mu, y, color=color, lw=2.0, label=name)
    ax.set_xlim(26.5, 14.5)
    ax.set_ylim(-0.02, 1.05)
    sky_sigma = float(art["sky"]["r"]["std"])
    mu_noise = float(cp.mag_arcsec2(sky_sigma))
    ax.axvline(mu_noise, color=ps.MUTED, lw=0.8, ls=DOTTED)
    ax.text(mu_noise - 0.15, 0.45, "pixel noise\n(1σ)", fontsize=7.5, color=ps.INK_SECONDARY, va="center")
    ax.axvline(art["peak_mu_r"], color=ps.MUTED, lw=0.8, ls=DOTTED)
    ax.text(art["peak_mu_r"] + 0.15, 0.06, "core of\n" + g.name, fontsize=7.5, color=ps.INK_SECONDARY,
            ha="right")
    _label(ax, "linear", (17.0, float(cp.mag_to_nmgy(17.0) / LUPTON_SATURATION)), (18.6, 0.2))
    _label(ax, "arcsinh", (21.5, float(asin(cp.mag_to_nmgy(21.5)))), (19.8, 0.34))
    _label(ax, "log", (20.0, float(logc(cp.mag_to_nmgy(20.0)))), (18.6, 0.93))
    ax.set_xlabel("surface brightness (r mag per arcsec², brighter →)")
    ax.set_ylabel("screen brightness (0 = black, 1 = white)")
    ax.legend(loc="upper left", fontsize=7.5, bbox_to_anchor=(0.0, 0.86))
    _title(ax, "The three curves", "how sky brightness becomes\nscreen brightness", fontsize=10)

    # colour fidelity
    ax = fig.add_subplot(gs[1, 2])
    r = profile_display["r_over_r25"]
    ax.plot(r, profile_display["true"], color=ps.INK, lw=2.4, label="measured (the truth)")
    ax.plot(r, profile_display["lupton"], color=ps.ORANGE, lw=1.5, ls=DASHED,
            label="arcsinh on total light")
    ax.plot(r, profile_display["per_channel"], color=ps.BLUE, lw=1.6, label="arcsinh per channel")
    sun_gr = ref["sun"]["colors"]["g-r"]
    ax.axhline(sun_gr, color=ps.MUTED, lw=0.8, ls=DOTTED)
    ax.text(0.01, sun_gr + 0.006, "the Sun's colour (shown as white)", fontsize=7.5,
            color=ps.INK_SECONDARY, ha="left", va="bottom")
    ax.set_xlabel(f"distance from the centre of {g.name} (fraction of its radius)")
    ax.set_ylabel("g − r colour read back from the picture (mag)")
    ax.legend(loc="upper right", fontsize=7.5)
    _title(ax, "Which picture tells the truth?", "the colour each picture displays,\nring by ring",
           fontsize=10)

    fig.suptitle(f"The stretch problem: {g.name} spans a factor of ~10,000 in brightness; a screen "
                 "shows 256 levels", x=0.125, ha="left", y=0.995)
    _caption(fig, "All four pictures use the same three calibrated images and the same white balance; "
             "only the stretch differs. The arcsinh curve is linear for faint light and logarithmic for "
             "bright light. Applied to each channel separately it squeezes the bright bulge's channels "
             "together, draining its colour; applied to the total light with one common factor per pixel "
             "(Lupton et al. 2004, the algorithm behind SDSS's own colour images) it keeps every "
             "pixel's colour -- the dashed orange line lies on the black one. " + PROVENANCE, y=0.035)
    return _save(fig, "02_the_stretch", svg_dpi=90)


# --------------------------------------------------------------------------
# 3. Filters versus the eye
# --------------------------------------------------------------------------

def figure_filters_vs_eye(ref: dict, eyeband: dict, m51_spectrum) -> list[str]:
    fig, (ax_eye, ax_f) = plt.subplots(2, 1, figsize=(11.5, 8.2), sharex=True,
                                       gridspec_kw={"hspace": 0.34, "height_ratios": [1, 1.25]})
    cl, cmf = cp.load_cmf()
    names = [("x̄  (mostly the red-sensitive cones)", "-"), ("ȳ  (brightness; peaks at 555 nm)", DASHED),
             ("z̄  (the blue-sensitive cones)", DOTTED)]
    for k, (name, ls) in enumerate(names):
        ax_eye.plot(cl / 10, cmf[:, k] / cmf[:, 1].max(), color=ps.C_EYE, lw=1.6, ls=ls, label=name)
    sl, sv = cp.load_scotopic()
    ax_eye.plot(sl / 10, sv, color=ps.VIOLET, lw=1.4, label="rods: night vision (no colour)")
    # the visible spectrum as a colour strip
    vis = np.linspace(380, 720, 341)
    def spectral(v):
        # a pure wavelength lies outside any screen's gamut; clip, the usual way
        rgb = np.clip(cp.xyz_to_display_rgb(cmf[np.argmin(np.abs(cl / 10 - v))]), 0, None)
        return cp.srgb_encode(rgb / rgb.max())
    strip = np.array([spectral(v) for v in vis])
    ax_eye.imshow(strip[None, :, :], extent=(380, 720, 1.83, 1.93), aspect="auto", zorder=3)
    ax_eye.text(550, 1.97, "the colours of single wavelengths (approximate, as a screen can show them)",
                fontsize=7.5, color=ps.INK_SECONDARY, ha="center", va="bottom")
    ax_eye.set_ylim(0, 2.1)
    ax_eye.set_ylabel("sensitivity (ȳ peak = 1)")
    ax_eye.legend(loc="upper right", fontsize=7.5)
    _title(ax_eye, "The human eye",
           "the CIE 1931 colour-matching functions: any two lights that give the same three numbers look identical")

    for b in "ugriz":
        fl, fr = cp.load_filter(f"sdss_{b}")
        lw, alpha = (2.2, 1.0) if b in BANDS else (1.0, 0.55)
        ax_f.plot(fl / 10, fr / fr.max(), color=ps.C_BAND[b], lw=lw, alpha=alpha)
        peak = fl[np.argmax(fr)] / 10
        what = {"g": "→ blue", "r": "→ green", "i": "→ red", "u": "(not used)", "z": "(not used)"}[b]
        ax_f.text(peak, 1.03, f"{b} {what}", ha="center", va="bottom", fontsize=8,
                  color=ps.INK_SECONDARY, fontweight="bold" if b in BANDS else "normal")
    sun_l, sun_f = cp.load_sun()
    sel = (sun_l > 3000) & (sun_l < 11000)
    ax_f.fill_between(sun_l[sel] / 10, 0, sun_f[sel] / sun_f[sel].max() * 0.9, color=ps.C_CONTEXT,
                      alpha=0.8, lw=0, zorder=0, label="the Sun's spectrum")
    lam, flam = m51_spectrum
    sel = (lam > 3000) & (lam < 11000)
    ax_f.plot(lam[sel] / 10, flam[sel] / np.percentile(flam[sel], 99.5) * 0.85, color=ps.INK, lw=0.8,
              label="M51's integrated spectrum (Brown et al. 2014)")
    ax_f.axvline(656.3, color=ps.MUTED, lw=0.8, ls=DOTTED)
    _label(ax_f, "Hα, 656 nm: hydrogen's deep-red glow falls in r, which is shown as green",
           (656.3, 1.13), (672, 1.2), va="center")
    ax_f.set_ylim(0, 1.28)
    ax_f.set_xlim(300, 1100)
    ax_f.set_xlabel("wavelength (nm)")
    ax_f.set_ylabel("filter response / spectrum (scaled)")
    ax_f.legend(loc="center right", fontsize=7.5, bbox_to_anchor=(1.0, 0.62))
    rel = eyeband
    _title(ax_f, "The SDSS camera",
           f"what the eye makes of each filter's light (ȳ averaged over the filter): g {rel['g']:.2f}, "
           f"r {rel['r']:.2f}, i {rel['i']:.4f} -- the i filter is essentially invisible")
    for ax in (ax_eye, ax_f):
        ax.axvspan(380, 700, color=ps.GRID, alpha=0.25, lw=0, zorder=0)
    fig.suptitle("Why survey colours are not the colours you would see", x=0.125, ha="left", y=0.995)
    _caption(fig, "Filter curves: SDSS system response (Doi et al. 2010, via the SVO Filter Profile Service). "
             "Eye: CIE 1931 2° colour-matching functions and the CIE 1951 scotopic (rod) function. "
             "Shaded band: 380-700 nm. The i filter samples light the eye barely registers, and r "
             "(displayed as green) contains the deep-red glow of hydrogen.", y=0.045)
    return _save(fig, "03_filters_vs_the_eye")


# --------------------------------------------------------------------------
# 4. Validation
# --------------------------------------------------------------------------

def figure_validation(val: dict, ref: dict) -> list[str]:
    fig, axes = plt.subplots(2, 3, figsize=(13.6, 8.4), gridspec_kw={"wspace": 0.32, "hspace": 0.62})

    # (a) the survey's own JPEG
    ax = axes[0, 0]
    bins = np.arange(-60.5, 61.5, 1)
    raw = np.concatenate([v["diff_raw"] for v in val["jpeg_pixels"]])
    comp = np.concatenate([v["diff_compressed"] for v in val["jpeg_pixels"]])
    ax.hist(raw, bins=bins, color=ps.C_CONTEXT, label="our rendering vs theirs")
    ax.hist(comp, bins=bins, color=ps.C_OURS, label="…after the same JPEG compression")
    ax.set_yscale("log")
    ax.set_xlabel("difference in 8-bit value (ours − survey's)")
    ax.set_ylabel("pixel values")
    ax.legend(fontsize=7.5, loc="upper right")
    frac = np.mean(np.abs(comp) <= 1)
    _title(ax, "We read the pixels as the survey does",
           f"{frac * 100:.3f}% of {comp.size / 1e6:.1f} M values within ±1\n"
           f"({len(val['jpeg'])} galaxies; the survey's JPEG quality {val['jpeg'][0]['jpeg_quality']})", fontsize=10)

    # (b, c) stars
    stars = val["stars"]
    t = stars["table"]
    import pandas as pd
    df = pd.DataFrame(t)
    for ax, cat, name in ((axes[0, 1], "sdss", "SDSS's PSF photometry"), (axes[0, 2], "gaia", "Gaia's spectra")):
        d = df[df["catalog"] == cat]
        for b in BANDS:
            ax.scatter(d[f"cat_{b}"], d[f"ours_{b}"] - d[f"cat_{b}"], s=9, color=ps.C_BAND[b],
                       alpha=0.55, lw=0, label=f"{b}")
        ax.axhline(0, color=ps.MUTED, lw=0.8)
        ax.set_ylim(-0.3, 0.3)
        ax.set_xlabel("catalogue magnitude")
        ax.set_ylabel("ours − catalogue (mag)")
        ax.legend(fontsize=7.5, loc="lower left", ncol=3, markerscale=1.6)
        s = stars[cat]
        sub = ", ".join(f"{b} {s[b]['median']:+.3f}" for b in BANDS)
        _title(ax, f"Stars: our images vs {name}",
               f"{s['r']['n']} isolated stars; median offsets {sub}\n"
               f"colours: g−r {s['g-r']['median']:+.3f}, r−i {s['r-i']['median']:+.3f}", fontsize=10)

    # (d) Brown et al. apertures
    ax = axes[1, 0]
    rows = val["brown"]
    names = [r["name"].replace("_", " ").replace("NGC 0", "NGC ") for r in rows]
    y = np.arange(len(rows))
    ax.axvspan(-0.05, 0.05, color=ps.GRID, alpha=0.6, lw=0)
    ax.set_xlim(-0.16, 0.16)
    for k, b in enumerate(BANDS):
        ax.scatter([r["mag_diff"][b] for r in rows], y + (k - 1) * 0.22, s=22, color=ps.C_BAND[b],
                   label=b, zorder=3, edgecolor=ps.SURFACE, lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=7.5)
    ax.invert_yaxis()
    ax.axvline(0, color=ps.MUTED, lw=0.8)
    ax.set_xlabel("ours − Brown et al. 2014 (mag)")
    ax.legend(fontsize=7.5, loc="lower right", ncol=1)
    diffs = np.array([[r["mag_diff"][b] for b in BANDS] for r in rows])
    cdiff = np.array([[r["our_colors"][c] - r["published_colors"][c] for c in ("g-r", "r-i")] for r in rows])
    _title(ax, "Whole galaxies: Brown et al.'s apertures",
           f"median {np.median(diffs):+.3f} mag; colours agree to {np.median(np.abs(cdiff)):.3f} mag\n"
           f"(median); grey: their ±0.05 mag error floor", fontsize=10)

    # (e) stitching
    ax = axes[1, 1]
    x = np.arange(3)
    for k, st in enumerate(val["stitch"]):
        ax.bar(x + (k - 0.5) * 0.36, [st[b]["typical_diff_over_noise"] for b in BANDS], width=0.34,
               color=[ps.C_OURS, ps.C_PUBLISHED][k], label=st["galaxy"].upper())
    ax.axhline(1.0, color=ps.MUTED, lw=0.8, ls=DASHED)
    ax.text(2.45, 1.02, "the pixel noise itself", fontsize=7.5, color=ps.INK_SECONDARY, ha="right", va="bottom")
    ax.set_xticks(x)
    ax.set_xticklabels(list(BANDS))
    ax.set_ylim(0, 1.25)
    ax.set_xlabel("SDSS filter")
    ax.set_ylabel("typical |stitched − single| / pixel noise")
    ax.legend(fontsize=7.5, loc="upper left")
    shifts = [max(abs(v) for v in st[b]["shift_px"]) for st in val["stitch"] for b in BANDS]
    worst_flux = max(abs(st[b]["total_flux_ratio"] - 1) for st in val["stitch"] for b in BANDS)
    peaks = max(st[b]["bright_pixels_median_abs_frac"] for st in val["stitch"] for b in BANDS)
    _title(ax, "Stitching tiles changes nothing that matters",
           f"offset ≤ {max(shifts):.3f} px; light conserved to {worst_flux * 100:.3f}%;\n"
           f"bright pixels change by {peaks * 100:.1f}% (median)", fontsize=10)

    # (f) the Sun
    ax = axes[1, 2]
    sun = ref["sun"]
    bands5 = list(sun["abs_mag_ab_published"])
    ax.axhspan(-0.02, 0.02, color=ps.GRID, alpha=0.6, lw=0)
    ax.scatter(range(len(bands5)), [sun["abs_mag_diff"][b] for b in bands5], s=40, color=ps.C_OURS,
               zorder=3, edgecolor=ps.SURFACE)
    ax.axhline(0, color=ps.MUTED, lw=0.8)
    ax.set_xticks(range(len(bands5)))
    ax.set_xticklabels(bands5)
    ax.set_ylim(-0.05, 0.05)
    ax.set_xlabel("SDSS filter")
    ax.set_ylabel("ours − Willmer 2018 (mag)")
    _title(ax, "The Sun through the SDSS filters",
           "synthetic photometry of its measured spectrum;\nthese colours define 'white' in our pictures",
           fontsize=10)

    fig.suptitle("Checking the measurements before trusting the pictures", x=0.125, ha="left", y=1.0)
    _caption(fig, "Top left: the viewer's SDSS JPEG recipe applied to our copy of the pixels, compared with "
             "the JPEG it served. Top middle/right: 7.4\" aperture photometry of isolated stars against SDSS "
             "DR17 PSF magnitudes and Gaia DR3 synthetic SDSS magnitudes (GSPC). Bottom left: our "
             "photometry inside the rectangular apertures of Brown et al. (2014, Table 2). Bottom middle: "
             "tiles stitched with reproject vs the same frame fetched in one piece. Bottom right: Willmer "
             "(2018, ApJS 236, 47).", y=0.035)
    return _save(fig, "04_validation")


# --------------------------------------------------------------------------
# 5. The gallery
# --------------------------------------------------------------------------

def figure_gallery(results: list, arts: dict) -> list[str]:
    n = len(results)
    ncol = 4
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(13.2, 3.75 * nrow),
                             gridspec_kw={"wspace": 0.04, "hspace": 0.3})
    for ax, res in zip(axes.flat, results):
        art = arts[res.key]
        _image(ax, art["composite"]["standard"], fov_arcmin=art["frame"].galaxy.fov_arcmin, max_px=480)
        ax.set_title(f"{res.name} · {res.common}", fontsize=9.5, pad=4)
        ax.text(0.5, -0.03, res.kind, transform=ax.transAxes, ha="center", va="top", fontsize=7.8,
                color=ps.INK_SECONDARY)
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    fig.suptitle("Twelve galaxies, one camera, one recipe", x=0.125, ha="left", y=0.995)
    fig.text(0.125, 0.965, "Same filters (SDSS g, r, i → blue, green, red), same white point (the Sun "
             "is white), same arcsinh stretch in physical units: equal surface brightness is equal "
             "screen brightness in every panel.", fontsize=8.5, color=ps.INK_SECONDARY, ha="left")
    _caption(fig, f"Stretch: Lupton et al. (2004) arcsinh, Q = {LUPTON_Q:g}, saturating at "
             f"{LUPTON_SATURATION:g} nanomaggies/arcsec² ({float(cp.mag_arcsec2(LUPTON_SATURATION)):.1f} "
             "mag/arcsec²). Scale bars in arcminutes. " + PROVENANCE, y=0.02)
    return _save(fig, "05_gallery", svg_dpi=90)


def save_gallery_pngs(arts: dict) -> list[str]:
    """Full-resolution pictures, one JPEG each, north up, plus a small
    thumbnail for the README's table."""
    from PIL import Image

    out = []
    for key, art in arts.items():
        out.append(cm.save_jpeg(art["composite"]["standard"], GALLERY_DIR / f"{key}.jpg"))
        out.append(cm.save_jpeg(art["composite"]["true"], GALLERY_DIR / f"{key}_eye.jpg"))
        thumb = Image.fromarray(cm.to_uint8(art["composite"]["standard"])).resize((150, 150), Image.LANCZOS)
        path = GALLERY_DIR / f"{key}_thumb.jpg"
        thumb.save(path, quality=88)
        out.append(str(path))
    return out


def save_swatches(results: list) -> list[str]:
    """A small square of each galaxy's eye colour (sunlight = white), for the README."""
    from PIL import Image

    out = []
    for res in results:
        hexs = res.eye_from_image["hex"].lstrip("#")
        rgb = tuple(int(hexs[k:k + 2], 16) for k in (0, 2, 4))
        img = Image.new("RGB", (56, 56), (26, 26, 25))
        img.paste(Image.new("RGB", (44, 44), rgb), (6, 6))
        path = GALLERY_DIR / f"{res.key}_swatch.png"
        img.save(path)
        out.append(str(path))
    return out


# --------------------------------------------------------------------------
# 6. What "white" means
# --------------------------------------------------------------------------

def figure_white_point(art: dict, ref: dict) -> list[str]:
    fr = art["frame"]
    g = fr.galaxy
    wb = ref["white_balance"]
    v, f = wb["viewer"], wb["flat"]
    viewer_title = (f"White = a {v['neutral_star']} star" if v["neutral_star_distance_mag"] < 0.05
                    else "White = the survey viewer's default")
    viewer_sub = (f"grey = g−r {v['neutral_color']['g-r']:.2f}, r−i {v['neutral_color']['r-i']:.2f}: "
                  "no real star has\nthis colour, so most stars take on a cyan tint")
    cases = [("ours", "White = the Sun (G2V)", "our choice: grey means 'the colour of\nsunlight'"),
             ("viewer", viewer_title, viewer_sub),
             ("flat", "White = equal light in every filter",
              f"grey = g−r 0, r−i 0: close to an {f['neutral_star']} star,\nhotter than the Sun")]
    star_types = ["B0V", "A0V", "F0V", "G2V", "K0V", "M0V"]
    fig = plt.figure(figsize=(13.2, 5.9))
    gs = fig.add_gridspec(2, 3, height_ratios=[5, 0.9], hspace=0.06, wspace=0.08)
    for k, (key, title, sub) in enumerate(cases):
        w = wb[key]["weights"]
        rgb = cm.composite(fr.sb, w, Q=LUPTON_Q, saturation=LUPTON_SATURATION)
        ax = fig.add_subplot(gs[0, k])
        _image(ax, rgb, fov_arcmin=g.fov_arcmin, max_px=560, scalebar=(k == 0))
        _title(ax, title, sub, fontsize=10)
        axs = fig.add_subplot(gs[1, k])
        axs.set_xlim(-0.5, len(star_types) - 0.5)
        axs.set_ylim(-0.9, 0.7)
        axs.axis("off")
        axs.set_facecolor(ps.IMAGE_BG)
        axs.add_patch(Rectangle((-0.5, -0.45), len(star_types), 1.1, color="black", zorder=0))
        for j, sp in enumerate(star_types):
            col = ref["stars"][sp]
            # a star's display colour: its flux ratios through the same weights
            f_r = 1.0
            f_g = 10 ** (-0.4 * col["g-r"])
            f_i = 10 ** (0.4 * col["r-i"])
            rgbv = np.array([w["i"] * f_i, w["r"] * f_r, w["g"] * f_g])
            rgbv = rgbv / rgbv.max()
            axs.add_patch(Circle((j, 0.1), 0.3, color=rgbv, zorder=2))
            axs.text(j, -0.62, sp, ha="center", va="center", fontsize=7.5, color=ps.INK_SECONDARY)
    fig.suptitle(f"There is no neutral choice: {g.name} with three white balances", x=0.125, ha="left",
                 y=0.975)
    _caption(fig, "Same images, same stretch; only the relative weights of the g, r and i images differ. "
             "Beneath each picture: how stars of each spectral type would appear with those weights. "
             "What looks white is whatever the weights say is white.", y=0.06)
    return _save(fig, "06_what_white_means", svg_dpi=90)


# --------------------------------------------------------------------------
# 7. True colour
# --------------------------------------------------------------------------

def figure_true_color(art: dict, results: list, ref: dict, eye_summary: dict) -> list[str]:
    from matplotlib.lines import Line2D

    fr = art["frame"]
    g = fr.galaxy
    fig = plt.figure(figsize=(13.4, 8.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.3, 0.55], hspace=0.34, wspace=0.2)
    ax = fig.add_subplot(gs[0, 0])
    _image(ax, art["composite"]["standard"], fov_arcmin=g.fov_arcmin, max_px=560)
    _title(ax, "As survey pictures show it", "i, r, g → red, green, blue\n(sunlight = white)", fontsize=10)
    ax = fig.add_subplot(gs[0, 1])
    _image(ax, art["composite"]["true"], fov_arcmin=g.fov_arcmin, max_px=560)
    _title(ax, "As an eye would see it, if it could", "g, r, i → CIE XYZ → your screen\n(sunlight = white)",
           fontsize=10)

    # u'v' diagram near the white point
    ax = fig.add_subplot(gs[0, 2])
    temps, luv = cp.planckian_locus()
    upvp = np.column_stack([luv[:, 0], luv[:, 1] * 1.5])  # CIE 1960 v -> 1976 v'
    sel = (temps > 2500) & (temps < 30000)
    ax.plot(upvp[sel, 0], upvp[sel, 1], color=ps.MUTED, lw=1.2, zorder=1)
    for T in (3000, 4000, 5000, 6500, 10000):
        k = np.argmin(np.abs(temps - T))
        ax.plot(*upvp[k], "o", color=ps.MUTED, ms=2.5)
        ax.text(upvp[k, 0] + 0.003, upvp[k, 1] + 0.001, f"{T:,} K", fontsize=6.8, color=ps.INK_SECONDARY,
                ha="left", va="top")
    for sp, col in ref["stars"].items():
        if "Vega" in sp:
            continue
        p = cp.xyz_to_upvp(cp.xy_to_xyz(col["x"], col["y"]))
        ax.plot(*p, "o", ms=4, color=ps.C_CONTEXT, zorder=2)
    for res in results:
        e = res.eye_from_image
        p = cp.xyz_to_upvp(cp.xy_to_xyz(e["x"], e["y"]))
        if res.eye_from_spectrum:
            s_ = res.eye_from_spectrum
            q = cp.xyz_to_upvp(cp.xy_to_xyz(s_["x"], s_["y"]))
            ax.plot([p[0], q[0]], [p[1], q[1]], color=ps.MUTED, lw=0.6, zorder=2)
            ax.plot(*q, "D", ms=4.5, color=ps.C_PUBLISHED, mec=ps.SURFACE, mew=0.6, zorder=3)
        ax.plot(*p, "o", ms=6.5, color=ps.C_OURS, mec=ps.SURFACE, mew=0.8, zorder=4)
    sun = cp.xyz_to_upvp(cp.sun_white_xyz())
    ax.plot(*sun, "*", ms=11, color=ps.INK, zorder=5)
    _label(ax, "the Sun", tuple(sun), (0.182, 0.492))
    d65 = cp.xyz_to_upvp(cp.d65_white_xyz())
    ax.plot(*d65, "+", ms=9, color=ps.INK_SECONDARY, zorder=5, mew=1.2)
    _label(ax, "D65: a screen's white", tuple(d65), (0.178, 0.452), ha="left")
    handles = [Line2D([], [], marker="o", ls="", color=ps.C_OURS, ms=6, label="galaxy, from our images"),
               Line2D([], [], marker="D", ls="", color=ps.C_PUBLISHED, ms=4.5,
                      label="same galaxy, from its spectrum"),
               Line2D([], [], marker="o", ls="", color=ps.C_CONTEXT, ms=4, label="stars, O to M (Pickles 1998)"),
               Line2D([], [], color=ps.MUTED, lw=1.2, label="glowing blackbodies")]
    ax.legend(handles=handles, fontsize=7, loc="lower right")
    ax.set_xlim(0.175, 0.265)
    ax.set_ylim(0.425, 0.53)
    ax.set_aspect("equal")
    ax.set_xlabel("u′  (CIE 1976)")
    ax.set_ylabel("v′")
    _title(ax, "Every galaxy is a shade of white", "colour with brightness removed; one just-\n"
           "noticeable difference ≈ 0.004 here", fontsize=10)

    # swatches of every galaxy's total light
    ax = fig.add_subplot(gs[1, :])
    n = len(results)
    ax.set_xlim(-0.55, n - 0.45)
    ax.set_ylim(-0.95, 1.2)
    ax.axis("off")
    ax.add_patch(Rectangle((-0.55, -0.95), n + 0.1, 2.15, color="#1a1a19", zorder=0))
    order = sorted(results, key=lambda r: -r.eye_from_image["cct_K"])
    for j, res in enumerate(order):
        e = res.eye_from_image
        ax.add_patch(Rectangle((j - 0.36, 0.1), 0.72, 0.72, color=e["hex"], zorder=2))
        if res.eye_from_spectrum:
            ax.add_patch(Rectangle((j - 0.36, -0.22), 0.72, 0.26, color=res.eye_from_spectrum["hex"], zorder=2))
        ax.text(j, 0.9, res.name, ha="center", va="bottom", fontsize=8, color="white")
        ax.text(j, -0.34, f"{e['cct_K']:,.0f} K", ha="center", va="top", fontsize=7.5, color="#c3c2b7")
        ax.text(j, -0.62, e["hex"], ha="center", va="top", fontsize=6.8, color="#898781")
    _title(ax, "The colour of each galaxy's total light, hottest to coolest",
           "large squares: from our images (inside the D25 ellipse); strips beneath: from the galaxy's "
           "integrated spectrum, where one exists; numbers: the blackbody temperature with the closest colour")

    fig.suptitle(f"True colour: what {g.name} would look like to a human eye", x=0.125, ha="left", y=0.965)
    _caption(fig, "Eye colours come from a 3×3 matrix mapping g, r, i fluxes to CIE XYZ, fitted to the integrated "
             "spectra of 129 galaxies (Brown et al. 2014); on spectra left out of the fit it errs by "
             f"Δu′v′ = {eye_summary['loo_duv_median']:.4f} (median), about one just-noticeable difference. "
             "Chromatic adaptation (Bradford) makes sunlight the screen's white. Brightness uses the same "
             "arcsinh curve as the survey-style picture, so the two pictures differ only in colour. "
             "Swatches sit on a dark surround, as galaxies do.", y=0.045)
    return _save(fig, "07_true_color", svg_dpi=90)


# --------------------------------------------------------------------------
# 8. Too faint for colour
# --------------------------------------------------------------------------

#: Short names for the everyday scenes in figure 8.
SCENE_SHORT = {
    "white paper in full sunlight": "paper in sunlight",
    "white paper on an office desk": "paper on a desk",
    "white paper by the light of a full Moon": "paper by moonlight",
    "the natural night sky, no Moon": "the night sky itself",
}


def figure_too_faint(art: dict, results: list, context: dict, dark_sky: float) -> list[str]:
    from matplotlib.lines import Line2D
    from scipy.ndimage import gaussian_filter

    fr = art["frame"]
    g = fr.galaxy
    fig = plt.figure(figsize=(13.2, 6.3))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.3], wspace=0.22)
    ax = fig.add_subplot(gs[0, 0])
    lum = art["luminance"]
    n = _image(ax, np.clip(lum, 1e-6, None), fov_arcmin=g.fov_arcmin, cmap="gray",
               norm=LogNorm(vmin=2e-5, vmax=0.05), max_px=700)
    # contours of the galaxy's light only: foreground stars are filled in first
    clean = lum.copy()
    clean[art["star_mask"]] = np.nan
    sm = gaussian_filter(np.nan_to_num(clean, nan=0.0), 2.0) / np.maximum(
        gaussian_filter((~np.isnan(clean)).astype(float), 2.0), 1e-3)
    sm = _downsample(sm, 700)
    ext = (0, n, 0, n)
    ax.contour(sm, levels=[dark_sky], colors=[ps.BLUE], linewidths=1.0, origin="lower", extent=ext)
    ax.contour(sm, levels=[MESOPIC_LOWER_CD_M2], colors=[ps.ORANGE], linewidths=1.5, origin="lower", extent=ext)
    ax.legend(handles=[Line2D([], [], color=ps.ORANGE, lw=1.5, label="0.005 cd/m²: cones begin to respond"),
                       Line2D([], [], color=ps.BLUE, lw=1.0, label="as bright as the natural night sky")],
              loc="lower right", fontsize=7, facecolor="#1a1a19", labelcolor="white", framealpha=0.85)
    _title(ax, f"{g.name}'s light, as an eye would receive it",
           "luminance (log scale, foreground stars left out of the contours):\ncolour vision needs more light "
           "than exists outside the tiny orange rings", fontsize=10)

    ax = fig.add_subplot(gs[0, 1])
    ax.set_xscale("log")
    xmin, xmax = 1e-6, 1e5
    ax.set_xlim(xmin, xmax)
    ax.axvspan(xmin, MESOPIC_LOWER_CD_M2, color=ps.GRID, alpha=0.9, lw=0)
    ax.axvspan(MESOPIC_LOWER_CD_M2, MESOPIC_UPPER_CD_M2, color=ps.GRID, alpha=0.4, lw=0)
    order = sorted(results, key=lambda r: r.luminance["peak_cd_m2"])
    top = len(order)
    for j, res in enumerate(order):
        L = res.luminance
        ax.plot([L["median_in_d25_cd_m2"], L["peak_cd_m2"]], [j, j], color=ps.C_OURS, lw=1.8,
                solid_capstyle="round")
        ax.plot(L["peak_cd_m2"], j, "o", color=ps.C_OURS, ms=6.5, mec=ps.SURFACE)
        ax.plot(L["median_in_d25_cd_m2"], j, "|", color=ps.C_OURS, ms=9, mew=1.5)
    yscene = top + 0.6
    for k, (name, value) in enumerate(sorted(context.items(), key=lambda kv: kv[1])):
        ax.plot(value, yscene, "s", color=ps.INK, ms=6)
        ax.text(value, yscene + (0.45 if k % 2 == 0 else -0.45), SCENE_SHORT.get(name, name), fontsize=7.5,
                ha="center", va="bottom" if k % 2 == 0 else "top", color=ps.INK_SECONDARY)
    ax.axhline(top - 0.2, color=ps.AXIS, lw=0.6)
    for x, text in ((3e-5, "rods only:\nno colour"), (0.17, "mesopic:\nfaint colour"), (700, "photopic:\nfull colour")):
        ax.text(x, top + 1.9, text, ha="center", va="bottom", fontsize=7.8, color=ps.INK_SECONDARY)
    ax.set_yticks(list(range(top)) + [yscene])
    ax.set_yticklabels([r.name for r in order] + ["for scale"], fontsize=8)
    ax.set_ylim(-0.8, top + 2.9)
    ax.set_xlabel("luminance (candela per square metre)")
    ax.grid(axis="y", visible=False)
    _title(ax, "Twelve galaxies against everyday light",
           "each bar runs from the median brightness inside the galaxy's D25 ellipse (tick) to its\n"
           "brightest core (dot), blurred by ~1.5\" of seeing as any eyepiece view would be", fontsize=10)
    fig.suptitle("No telescope can show you a galaxy's colours", x=0.125, ha="left", y=1.0)
    _caption(fig, "A telescope collects more light but spreads it over a bigger image, so the brightness per area "
             "you see can never exceed what the naked eye would get: these luminances are the most any eyepiece "
             "can deliver. Vision thresholds: CIE 191:2010 (mesopic range 0.005-5 cd/m²). Night sky: 22.0 V "
             "mag/arcsec² (Falchi et al. 2016). Paper: 80% reflective, lit by 100,000 lux (sun), 500 lux (desk), "
             "0.25 lux (full Moon).", y=0.035)
    return _save(fig, "08_too_faint_for_colour")


# --------------------------------------------------------------------------
# 9. What the colours mean
# --------------------------------------------------------------------------

def figure_colour_meaning(art: dict, results: list, ref: dict, color_map: np.ndarray) -> list[str]:
    fr = art["frame"]
    g = fr.galaxy
    fig = plt.figure(figsize=(13.2, 6.3))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.2], wspace=0.42)
    ax = fig.add_subplot(gs[0, 0])
    sun_gr = ref["sun"]["colors"]["g-r"]
    norm = TwoSlopeNorm(vmin=sun_gr - 0.6, vcenter=sun_gr, vmax=sun_gr + 0.6)
    img = _downsample(color_map, 700)
    n = img.shape[0]
    im = ax.imshow(np.ma.masked_invalid(img), origin="lower", cmap=CMAP_COLOR_INDEX, norm=norm,
                   interpolation="nearest", extent=(0, n, 0, n))
    ax.set_facecolor("#1a1a19")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    ticks = {sp: ref["stars"][sp]["g-r"] for sp in ("A0V", "F0V", "G2V", "K0V", "K0III", "M0III")}
    ticks = {k: v for k, v in ticks.items() if norm.vmin <= v <= norm.vmax}
    cb.set_ticks(list(ticks.values()))
    cb.set_ticklabels([f"{k}  {v:.2f}" for k, v in ticks.items()], fontsize=7)
    cb.set_label("g − r (mag); ticks: single stars (G2V = the Sun)", fontsize=8)
    _title(ax, f"{g.name}: the colour of every patch", "blue = bluer than sunlight, red = redder;\n"
           "shown only where the signal is strong enough to measure", fontsize=10)

    ax = fig.add_subplot(gs[0, 1])
    show = ["m87", "m104", "m51", "m74", "ngc4631", "ngc5866"]
    colors = [ps.BLUE, ps.ORANGE, ps.AQUA, ps.YELLOW, ps.MAGENTA, ps.GREEN]
    by_key = {r.key: r for r in results}
    for key, color in zip(show, colors):
        if key not in by_key:
            continue
        p = by_key[key].profile
        x = np.array(p["r_over_r25"])
        y = np.array(p["g-r"])
        e = np.array(p["g-r_err"])
        ok = np.isfinite(y) & (e < 0.05)
        ax.plot(x[ok], y[ok], color=color, lw=1.9, label=f"{by_key[key].name} ({by_key[key].kind})")
    for sp in ("F5V", "G2V", "K0V", "K0III"):
        v = ref["stars"][sp]["g-r"]
        ax.axhline(v, color=ps.MUTED, lw=0.7, ls=DOTTED)
        ax.text(1.215, v, sp + (" (the Sun)" if sp == "G2V" else ""), fontsize=7.3,
                color=ps.INK_SECONDARY, va="center", ha="left")
    ax.set_xlim(0, 1.2)
    ax.set_xlabel("distance from the centre (fraction of the galaxy's D25 radius)")
    ax.set_ylabel("g − r (mag), Milky Way dust removed")
    ax.legend(fontsize=7.3, loc="lower left")
    _title(ax, "Red centres, bluer discs", "the colour of the light in elliptical rings; dotted lines: single stars "
           "of each type\n-- a galaxy's light is the sum of billions of them", fontsize=10)
    fig.suptitle("What the colours mean: old stars are red, young stars are blue, dust reddens",
                 x=0.125, ha="left", y=1.0)
    _caption(fig, "Colours from the calibrated images, corrected for reddening by Milky Way dust (E(B−V) from "
             "Schlegel et al. 1998 with the Schlafly & Finkbeiner 2011 coefficients). Rings follow each "
             "galaxy's RC3 ellipse, after subtracting the residual sky level; rings with bootstrap errors above "
             "0.05 mag, fainter than 24 r mag/arcsec², or less than ten times brighter than that residual are "
             "dropped. Star colours: "
             "synthetic photometry of the Pickles (1998) spectral library.", y=0.03)
    return _save(fig, "09_what_colours_mean")


# --------------------------------------------------------------------------
# 12. When a pipeline erases a galaxy
# --------------------------------------------------------------------------

def figure_pipeline_caveat(caveat: dict, sdss_rgb: np.ndarray, viewer_jpeg: np.ndarray, fov_arcmin: float) -> list[str]:
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.7), gridspec_kw={"wspace": 0.3, "width_ratios": [1, 1, 1.0]})
    _image(axes[0], sdss_rgb, fov_arcmin=fov_arcmin, max_px=520)
    _title(axes[0], "M51 in SDSS", "the images this project uses", fontsize=10)
    _image(axes[1], viewer_jpeg[::-1] / 255.0, fov_arcmin=fov_arcmin, max_px=520, scalebar=False)
    _title(axes[1], "M51 in Legacy Surveys DR9-north", "the survey viewer's own picture: cyan, with\n"
           "rectangular seams where image tiles meet", fontsize=10)
    ax = axes[2]
    ret = caveat["retained"]
    bands = ["g", "r", "z"]
    ax.bar(range(3), [100 * ret[b] for b in bands], color=[ps.C_BAND["g"], ps.C_BAND["r"], ps.C_BAND["z"]],
           width=0.6)
    for k, b in enumerate(bands):
        ax.text(k, 100 * ret[b] + 2, f"{100 * ret[b]:.0f}%", ha="center", va="bottom", fontsize=9,
                color=ps.INK)
    ax.axhline(100, color=ps.MUTED, lw=0.8, ls=DASHED)
    ax.set_xticks(range(3))
    ax.set_xticklabels([f"{b}" for b in bands])
    ax.set_ylim(0, 125)
    ax.set_ylabel("light inside M51's D25 ellipse (% of expected)")
    ax.set_xlabel("filter")
    ax.grid(axis="x", visible=False)
    _title(ax, "Most of the red light is gone", "expected: SDSS for g and r; for z, SDSS i and\n"
           f"M51's i − z colour ({caveat['i_minus_z_brown']:.2f}, Brown et al. 2014)", fontsize=10)
    fig.suptitle("Why the gallery uses SDSS: a survey pipeline can erase a big galaxy", x=0.125, ha="left",
                 y=1.04)
    _caption(fig, "Sky subtraction models the background with smooth functions over each image or tile. A galaxy "
             "bigger than that scale looks like background, and gets subtracted with it. In DR9-north's z band "
             "this removed about two-thirds of M51's light, which turns the galaxy cyan (z drives the red "
             "channel). The Legacy Surveys team built custom mosaics for such galaxies (the Siena Galaxy Atlas); "
             "SDSS's DR8 sky model (Blanton et al. 2011), designed with large galaxies in mind, keeps them. Data: Legacy Surveys DR9 "
             "(Dey et al. 2019).", y=0.02)
    return _save(fig, "12_when_a_pipeline_erases_a_galaxy", svg_dpi=90)
