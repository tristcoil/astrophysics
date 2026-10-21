"""M51 across the spectrum, and through Hubble's eyes.

Every image here comes from a different telescope, each with its own pixel
grid, orientation and resolution: GALEX's ultraviolet tiles, a ground-based
H-alpha image, Spitzer's infrared mosaics, a VLA radio map, Hubble's ACS
mosaic. ``reproject`` puts each one onto M51's optical grid, pixel for pixel,
which is what lets them be compared or combined at all. How well they line up
afterwards is measured, not assumed (``alignment``): foreground stars with
Gaia positions are centroided in every map that shows them.

The colours in this module's composites are *assigned*, not seen:
ultraviolet, infrared and radio light have no colour. The convention is
chromatic order -- shorter wavelengths bluer -- and every panel names the light
it shows.
"""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

import color_physics as cp
import composite as cm
import imaging
import plotstyle as ps
from config import BANDS, DATA_DIR, LUPTON_Q, LUPTON_SATURATION, galaxy_by_key, target_wcs

M51 = galaxy_by_key("m51")
MW = DATA_DIR / "multiwavelength"


@dataclass(frozen=True)
class Layer:
    key: str
    filename: str
    #: which plane of a data cube to use (the GALEX cutout holds NUV then FUV)
    plane: int | None
    title: str
    band: str
    traces: str
    #: approximate resolution (FWHM, arcsec), for the panel labels
    resolution_arcsec: float
    #: whether foreground stars are visible, so alignment can be measured
    has_stars: bool = True


LAYERS = [
    Layer("fuv", "m51_galex.fits", 1, "Far ultraviolet", "GALEX FUV · 0.15 µm",
          "hot young stars (< ~100 Myr old)", 4.2, has_stars=False),
    Layer("nuv", "m51_galex.fits", 0, "Near ultraviolet", "GALEX NUV · 0.23 µm",
          "young stars", 5.3),
    Layer("halpha", "ngc5194_halpha_sub.fits", None, "Hydrogen-alpha", "KPNO · 0.656 µm, starlight subtracted",
          "gas lit by the youngest stars (< ~10 Myr)", 1.5, has_stars=False),
    Layer("irac1", "ngc5194_irac1.fits", None, "Near infrared", "Spitzer IRAC · 3.6 µm",
          "old stars: most of the stellar mass", 1.7),
    Layer("irac4", "ngc5194_irac4.fits", None, "Mid infrared", "Spitzer IRAC · 8.0 µm",
          "dust and PAH molecules warmed by starlight", 2.0),
    Layer("mips24", "ngc5194_mips24.fits", None, "Mid infrared", "Spitzer MIPS · 24 µm",
          "warm dust in star-forming regions", 6.0, has_stars=False),
    Layer("hi", "ngc5194_things_mom0.fits", None, "Radio", "VLA · 21 cm",
          "cold hydrogen gas, the raw material", 11.0, has_stars=False),
]
LAYER = {l.key: l for l in LAYERS}


def _strip_stray_sip(hdr: fits.Header) -> fits.Header:
    """Drop SIP distortion keywords that the header does not actually invoke.

    By the FITS convention SIP polynomials apply only when CTYPE ends in
    "-SIP". The SINGS Spitzer mosaics are already distortion-corrected but
    carry leftover SIP keywords from their input frames without the suffix;
    astropy would apply them anyway, and polynomials meant for a 256-pixel
    detector, evaluated across a 2,400-pixel mosaic, move positions by
    arcseconds. So they are removed unless CTYPE asks for them.
    """
    hdr = hdr.copy()
    if not str(hdr.get("CTYPE1", "")).endswith("-SIP"):
        for key in list(hdr.keys()):
            if key in ("A_ORDER", "B_ORDER", "AP_ORDER", "BP_ORDER", "A_DMAX", "B_DMAX") or \
                    (key[:2] in ("A_", "B_") and key[2:].replace("_", "").isdigit()) or \
                    (key[:3] in ("AP_", "BP_") and key[3:].replace("_", "").isdigit()):
                del hdr[key]
    return hdr


def _read(filename: str, plane: int | None) -> tuple[np.ndarray, WCS, fits.Header]:
    path = MW / filename
    if not path.exists():
        raise FileNotFoundError(path.name)
    with fits.open(path) as hdul:
        hdu = next(h for h in hdul if h.data is not None and np.asarray(h.data).ndim >= 2)
        data = np.asarray(hdu.data, dtype=np.float64)
        hdr = _strip_stray_sip(hdu.header)
    w = WCS(hdr).celestial
    data = np.squeeze(data)
    if data.ndim == 3:
        data = data[plane if plane is not None else 0]
    return data, w, hdr


def reprojected(key: str, target: WCS | None = None, shape=None) -> np.ndarray:
    """One map, reprojected onto M51's optical grid (or any other)."""
    from reproject import reproject_interp

    lay = LAYER[key]
    data, w, _ = _read(lay.filename, lay.plane)
    target = target if target is not None else target_wcs(M51)
    shape = shape or (M51.size_px, M51.size_px)
    out, _ = reproject_interp((data, w), target, shape_out=shape, order="bilinear")
    return out


# --------------------------------------------------------------------------
# Alignment
# --------------------------------------------------------------------------

def _centroid(img, x, y, half):
    """First-moment centroid of a star in a box, with the background and noise
    taken from the box's border so a star filling the box is not mistaken for
    noise. Two passes: the second is centred on the first's result."""
    ny, nx = img.shape
    for _ in range(2):
        x0, x1 = int(round(x)) - half, int(round(x)) + half + 1
        y0, y1 = int(round(y)) - half, int(round(y)) + half + 1
        if x0 < 0 or y0 < 0 or x1 > nx or y1 > ny:
            return None
        cut = img[y0:y1, x0:x1]
        if not np.all(np.isfinite(cut)):
            return None
        border = np.concatenate([cut[0], cut[-1], cut[1:-1, 0], cut[1:-1, -1]])
        bkg = np.median(border)
        noise = 1.4826 * np.median(np.abs(border - bkg))
        w = np.clip(cut - bkg, 0, None)
        if noise <= 0 or (cut - bkg).max() < 10 * noise or w.sum() <= 0:
            return None
        yy, xx = np.mgrid[y0:y1, x0:x1]
        x, y = float((w * xx).sum() / w.sum()), float((w * yy).sum() / w.sum())
    return x, y


MAX_PM_MAS_YR = 10.0


def alignment(images: dict[str, np.ndarray], frame: imaging.Frame, *, gmax: float = 17.5) -> dict:
    """Where foreground stars land in each reprojected map, against Gaia.

    Only maps in which stars are visible are tested; the ultraviolet,
    H-alpha (stars subtracted), 24 micron and radio maps are too sparse in
    stars and are placed by their own astrometry.
    """
    stars = imaging.foreground_stars(imaging.load_gaia(M51))
    # slow stars only: these images were taken ~2003-2005, Gaia's positions are
    # for 2016.0, and a star moving 10 mas/yr shifts 0.12" in the meantime
    stars = stars[(stars["phot_g_mean_mag"] < gmax) & (stars["pm"] < MAX_PM_MAS_YR)]
    xs, ys = frame.world_to_pixel(stars["ra"].to_numpy(), stars["dec"].to_numpy())
    out = {}
    for key, img in images.items():
        fwhm = LAYER[key].resolution_arcsec if key in LAYER else 1.5
        half = int(np.ceil(max(4.0, 2.0 * fwhm) / frame.pixscale))
        dx, dy = [], []
        for x, y in zip(xs, ys):
            c = _centroid(img, x, y, half)
            if c is not None:
                dx.append(c[0] - x)
                dy.append(c[1] - y)
        if len(dx) >= 3:
            dx, dy = np.array(dx) * frame.pixscale, np.array(dy) * frame.pixscale
            out[key] = {"n_stars": len(dx), "median_east_arcsec": float(-np.median(dx)),
                        "median_north_arcsec": float(np.median(dy)),
                        "median_total_arcsec": float(np.hypot(np.median(dx), np.median(dy))),
                        "scatter_arcsec": float(1.4826 * np.median(np.hypot(dx - np.median(dx),
                                                                            dy - np.median(dy))))}
        else:
            out[key] = {"n_stars": len(dx)}
    return out


# --------------------------------------------------------------------------
# Hubble
# --------------------------------------------------------------------------

HST_FILTERS = {"F435W": ("b", "acs_f435w"), "F555W": ("v", "acs_f555w"), "F814W": ("i", "acs_f814w")}

#: The mosaic was drizzled to 0.05"/pixel; the files used here are 4x4 block
#: *averages* of it (MAST readme), so each value is still per original 0.05"
#: pixel -- just as the SDSS cutouts are per native 0.396" pixel.
HST_NATIVE_PIXSCALE = 0.05
#: SDSS filter each Hubble filter is compared with.
HST_VS_SDSS = {"F435W": "g", "F555W": "r", "F814W": "i"}


def _hst(filter_name: str, offset_arcsec: tuple[float, float] = (0.0, 0.0)) -> tuple[np.ndarray, WCS, float]:
    """An ACS mosaic in nanomaggies per square arcsecond, with its WCS and pixel scale.

    The values are count rates (electrons per second; the header's BUNIT says
    ELECTRONS, but the photometry check in ``hubble_photometry_check`` settles
    that they are rates) per original 0.05" pixel. PHOTFLAM turns a count rate
    into f_lambda, and PHOTPLAM, the filter's pivot wavelength, gives f_nu.
    """
    letter, _ = HST_FILTERS[filter_name]
    path = MW / f"h_m51_{letter}_s20_drz_sci.fits"
    if not path.exists():
        raise FileNotFoundError(path.name)
    with fits.open(path) as hdul:
        hdu = next(h for h in hdul if h.data is not None and np.asarray(h.data).ndim == 2)
        data = np.asarray(hdu.data, dtype=np.float64)
        hdr = hdul[0].header.copy()
        hdr.update(hdu.header)
    w = WCS(_strip_stray_sip(hdu.header)).celestial
    if any(offset_arcsec):
        # register: move the whole image by (east, north) arcsec on the sky
        east, north = offset_arcsec
        w.wcs.crval = [w.wcs.crval[0] - east / 3600.0 / np.cos(np.radians(w.wcs.crval[1])),
                       w.wcs.crval[1] - north / 3600.0]
    photflam = float(hdr["PHOTFLAM"])
    photplam = float(hdr["PHOTPLAM"])
    pixscale = float(np.sqrt(np.abs(np.linalg.det(w.pixel_scale_matrix)))) * 3600.0
    fnu = data * photflam * photplam**2 / cp.C_ANGSTROM_S  # erg/s/cm2/Hz per pixel
    return fnu / cp.nanomaggy_in_cgs() / HST_NATIVE_PIXSCALE**2, w, pixscale


def hubble_registration(frame: imaging.Frame) -> dict:
    """How far Hubble's mosaic is from SDSS on the sky, by cross-correlation.

    The F814W mosaic is reprojected onto M51's frame, blurred to SDSS's
    ~1.4" seeing, and phase-correlated with the SDSS i image (the closest SDSS
    filter) over the region both cover. Thousands of features -- arms, knots,
    dust lanes -- enter, so the answer is far steadier than a handful of
    stars. Returns the offset of Hubble relative to SDSS, (east, north) in
    arcsec: the correction to apply.
    """
    from reproject import reproject_interp
    from scipy.ndimage import gaussian_filter

    from analysis import _subpixel_shift

    img, w, _ = _hst("F814W")
    on_frame, _ = reproject_interp((img, w), frame.wcs, shape_out=frame.shape, order="bilinear")
    ok = np.isfinite(on_frame)
    ys, xs = np.nonzero(ok)
    y0, y1, x0, x1 = ys.min() + 20, ys.max() - 20, xs.min() + 20, xs.max() - 20
    n = min(y1 - y0, x1 - x0) // 2 * 2
    cy, cx = (y0 + y1) // 2, (x0 + x1) // 2
    cut = (slice(cy - n // 2, cy + n // 2), slice(cx - n // 2, cx + n // 2))
    sigma_px = (1.4 / 2.355) / frame.pixscale
    a = gaussian_filter(np.nan_to_num(on_frame), sigma_px)[cut]
    b = frame.sb["i"][cut]
    a = np.arcsinh(a / np.percentile(a, 90))
    b = np.arcsinh(b / np.percentile(b, 90))
    dy, dx = _subpixel_shift(a - a.mean(), b - b.mean())
    return {"east_arcsec": float(-dx * frame.pixscale), "north_arcsec": float(dy * frame.pixscale),
            "region_px": int(n)}


def hubble_photometry_check(frame: imaging.Frame, offset) -> dict:
    """Hubble's surface brightness against SDSS's, pixel by pixel, over the
    bright parts of M51 both cover -- compared with the ratio M51's own
    integrated spectrum predicts for each pair of filters (they differ in
    wavelength, so the ratio is not 1)."""
    from reproject import reproject_interp

    lam, flam, _ = cp.load_brown().spectra["NGC_5194"]
    out = {}
    for filt, band in HST_VS_SDSS.items():
        img, w, _ = _hst(filt, offset)
        on, _ = reproject_interp((img, w), frame.wcs, shape_out=frame.shape, order="bilinear")
        ok = np.isfinite(on) & (frame.sb[band] > 3.0)
        measured = float(np.median(on[ok] / frame.sb[band][ok]))
        predicted = float(cp.band_flux(lam, flam, HST_FILTERS[filt][1]) / cp.band_flux(lam, flam, f"sdss_{band}"))
        out[f"{filt}/{band}"] = {"measured": measured, "predicted_from_spectrum": predicted,
                                 "measured_over_predicted": measured / predicted, "n_pixels": int(ok.sum())}
    return out


def hubble_white_weights() -> dict[str, float]:
    """Weights that make the Sun neutral through ACS F435W, F555W, F814W."""
    lam, flam = cp.load_sun()
    f = {k: cp.band_flux(lam, flam, v[1]) for k, v in HST_FILTERS.items()}
    return {k: f["F555W"] / f[k] for k in f}


def make_figures(ref: dict, eye) -> list[str]:
    for lay in LAYERS:
        if not (MW / lay.filename).exists():
            raise FileNotFoundError(lay.filename)
    written = figure_multiwavelength(ref)
    try:
        written += figure_hubble(ref)
    except FileNotFoundError as err:
        print(f"  (Hubble figure skipped: {err})")
    return written


def _stretch_mono(img: np.ndarray, *, noise: float, top: float) -> np.ndarray:
    """Arcsinh stretch for a single map: linear below ~2x the noise, logarithmic above,
    saturating at ``top``."""
    soft = 2.0 * noise
    return np.clip(np.arcsinh(np.clip(img, 0, None) / soft) / np.arcsinh(top / soft), 0, 1)


def _noise_and_top(img: np.ndarray, galaxy_mask: np.ndarray) -> tuple[float, float, float]:
    """(background, noise, top) of a map, from pixels away from the galaxy.

    Some maps are blanked (exactly zero or NaN) wherever nothing was detected
    -- the 21-cm moment map, for instance -- so blanked pixels are ignored; if
    no honest noise estimate is left, the stretch falls back to 1/300 of the
    top of the map.
    """
    finite = np.isfinite(img) & (img != 0)
    sky = img[finite & ~galaxy_mask]
    med = float(np.median(sky)) if sky.size > 100 else 0.0
    noise = float(1.4826 * np.median(np.abs(sky - med))) if sky.size > 100 else 0.0
    top = float(np.percentile(img[finite & galaxy_mask] - med, 99.7))
    if not noise > 0:
        noise = top / 300.0
    return med, noise, top


def figure_multiwavelength(ref: dict) -> list[str]:
    import figures as fg

    frame = imaging.load_frame(M51)
    main = imaging.rc3_entry(M51)
    gal = imaging.d25_mask(frame, main, scale=1.1) | imaging.d25_mask(
        frame, imaging.rc3_entry(M51, "NGC 5195"), scale=1.1)
    maps = {lay.key: reprojected(lay.key) for lay in LAYERS}
    align = alignment({"sdss_i": frame.sb["i"], **{k: maps[k] for k in ("nuv", "irac1", "irac4")}}, frame)

    std = cm.composite(frame.sb, ref["white_balance"]["ours"]["weights"], Q=LUPTON_Q,
                       saturation=LUPTON_SATURATION)
    show = ["fuv", "halpha", "optical", "irac1", "irac4", "mips24", "hi", "panchromatic"]
    fig, axes = plt.subplots(2, 4, figsize=(13.4, 7.9), gridspec_kw={"wspace": 0.06, "hspace": 0.42})
    stretched = {}
    for key in ("fuv", "halpha", "irac1", "irac4", "mips24", "hi"):
        img = maps[key]
        med, noise, top = _noise_and_top(img, gal)
        stretched[key] = _stretch_mono(np.nan_to_num(img - med), noise=noise, top=top)
    # a composite from three kinds of light, in chromatic order
    pan = np.stack([stretched["irac4"], stretched["irac1"], stretched["fuv"]], axis=-1)
    for ax, key in zip(axes.flat, show):
        if key == "optical":
            fg._image(ax, std, fov_arcmin=M51.fov_arcmin, max_px=520, scalebar=False)
            ax.set_title("Visible light", fontsize=9.5, pad=26, loc="left")
            ax.text(0, 1.015, "SDSS g, r, i · 0.47-0.75 µm\n~1.4\" resolution", transform=ax.transAxes,
                    fontsize=7.5, color=ps.INK_SECONDARY, va="bottom", linespacing=1.25)
            ax.text(0.5, -0.03, "starlight: all ages, reddened by dust", transform=ax.transAxes, ha="center",
                    va="top", fontsize=7.5, color=ps.INK_SECONDARY)
            continue
        if key == "panchromatic":
            fg._image(ax, pan, fov_arcmin=M51.fov_arcmin, max_px=520, scalebar=False)
            ax.set_title("Three lights at once", fontsize=9.5, pad=26, loc="left")
            ax.text(0, 1.015, "blue: far ultraviolet · green: 3.6 µm\nred: 8 µm (colours assigned)",
                    transform=ax.transAxes, fontsize=7.5, color=ps.INK_SECONDARY, va="bottom", linespacing=1.25)
            ax.text(0.5, -0.03, "young stars, old stars, dust",
                    transform=ax.transAxes, ha="center", va="top", fontsize=7.5, color=ps.INK_SECONDARY)
            continue
        lay = LAYER[key]
        fg._image(ax, stretched[key], fov_arcmin=M51.fov_arcmin, max_px=520, cmap="gray",
                  norm=plt.Normalize(0, 1), scalebar=(key == "fuv"))
        ax.set_title(lay.title, fontsize=9.5, pad=26, loc="left")
        ax.text(0, 1.015, f"{lay.band}\n~{lay.resolution_arcsec:g}\" resolution", transform=ax.transAxes,
                fontsize=7.5, color=ps.INK_SECONDARY, va="bottom", linespacing=1.25)
        ax.text(0.5, -0.03, lay.traces, transform=ax.transAxes, ha="center", va="top", fontsize=7.5,
                color=ps.INK_SECONDARY, wrap=False)
    fig.suptitle("M51 across the spectrum: seven telescopes, one grid", x=0.125, ha="left", y=0.995)
    names = {"sdss_i": "SDSS i", "nuv": "GALEX NUV", "irac1": "3.6 µm", "irac4": "8 µm"}
    tested = [f"{names.get(k, k)} {v['median_total_arcsec']:.2f}\" ({v['n_stars']} stars)"
              for k, v in align.items() if "median_total_arcsec" in v]
    fg._caption(fig, "Each map was reprojected from its own archive grid onto the optical frame with reproject. "
                "Alignment, from slow-moving foreground stars against their Gaia positions: " + "; ".join(tested)
                + ". "
                "Greyscale maps: arcsinh stretch, linear up to twice each map's noise. Data: GALEX (via the "
                "Legacy Surveys viewer); SINGS Hα and Spitzer (Kennicutt et al. 2003; Dale et al. 2007); "
                "THINGS (Walter et al. 2008); SDSS DR8.", y=0.04)
    written = fg._save(fig, "10_m51_across_the_spectrum", svg_dpi=90)
    figure_multiwavelength.alignment = align
    return written


def figure_hubble(ref: dict) -> list[str]:
    """SDSS against Hubble on the same patch of sky, same filters' worth of colour."""
    import figures as fg
    from reproject import reproject_interp

    frame = imaging.load_frame(M51)
    weights = hubble_white_weights()
    # a 2.6' x 2.6' box on M51's northern arms and the bridge to NGC 5195, at 0.2"/pixel
    box = 2.6
    zoom = galaxy_by_key("m51")
    zg = type(zoom)(zoom.key, zoom.name, zoom.common, zoom.ngc, 202.4740, 47.2140, box, 0.2,
                    zoom.kind, zoom.why)
    zw = target_wcs(zg)
    shape = (zg.size_px, zg.size_px)
    reg = hubble_registration(frame)
    offset = (reg["east_arcsec"], reg["north_arcsec"])
    hst = {}
    for filt in HST_FILTERS:
        img, w, _ = _hst(filt, offset)
        hst[filt], _ = reproject_interp((img, w), zw, shape_out=shape, order="bilinear")
    sdss = {}
    for b in BANDS:
        sdss[b], _ = reproject_interp((frame.sb[b], frame.wcs), zw, shape_out=shape, order="bicubic")
    sdss_rgb = cm.composite(sdss, ref["white_balance"]["ours"]["weights"], Q=LUPTON_Q,
                            saturation=LUPTON_SATURATION)
    # alignment against Gaia, before and after registering to SDSS
    checks = {}
    for label, off in (("hst_f814w_as_published", (0.0, 0.0)), ("hst_f814w_registered", offset)):
        img814, w814, _ = _hst("F814W", off)
        on_frame, _ = reproject_interp((img814, w814), frame.wcs, shape_out=frame.shape, order="bilinear")
        checks[label] = on_frame
    figure_hubble.alignment = {"registration_to_sdss": reg, **alignment(checks, frame),
                               "photometry_vs_sdss": hubble_photometry_check(frame, offset)}
    R, G, B = hst["F814W"] * weights["F814W"], hst["F555W"], hst["F435W"] * weights["F435W"]
    hst_rgb = cm.lupton_rgb(np.nan_to_num(R), np.nan_to_num(G), np.nan_to_num(B),
                            softening=cm.softening_for(LUPTON_Q, LUPTON_SATURATION), saturation=LUPTON_SATURATION)

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 6.6), gridspec_kw={"wspace": 0.04})
    fg._image(axes[0], sdss_rgb, fov_arcmin=box, max_px=800)
    fg._title(axes[0], "From the ground: SDSS", "2.5 m telescope, 54 s per filter, ~1.4\" seeing; g, r, i",
              fontsize=10)
    fg._image(axes[1], hst_rgb, fov_arcmin=box, max_px=800)
    fg._title(axes[1], "From orbit: Hubble ACS", "2.4 m telescope above the air, 23-45 min per filter; "
              "F435W, F555W, F814W\n(here block-averaged to 0.2\"/pixel: ~0.4\" detail)", fontsize=10)
    fig.suptitle("The same 2.6′ of M51 and the same recipe, from the ground and from orbit", x=0.125, ha="left",
                 y=1.0)
    phot = figure_hubble.alignment["photometry_vs_sdss"]
    agree = ", ".join(f"{k.split('/')[0]} {100 * (v['measured_over_predicted'] - 1):+.0f}%" for k, v in phot.items())
    fg._caption(fig, "Both pictures: white = the Sun (each through its own filters) and the same arcsinh stretch in "
                "physical units (nanomaggies per arcsec²), on one 0.2\"/pixel grid. Hubble's filters span a wider "
                "range (435/555/814 nm vs SDSS's 470/618/749 nm), which stretches the colour differences: its "
                "redder look is the filters, not the galaxy. No H-alpha image is mixed in. The Hubble mosaic was "
                f"registered to SDSS by cross-correlating the galaxy ({reg['east_arcsec']:.2f}\" east, "
                f"{reg['north_arcsec']:.2f}\" north; Gaia stars agree). Its brightness matches SDSS's to within "
                f"what M51's spectrum predicts for each filter pair: {agree}. The two diagonal streaks in the SDSS "
                "panel are satellite trails in the SDSS exposures. Hubble: Mutchler et al. 2005 (Hubble Heritage, "
                "MAST HLSP). SDSS DR8.", y=0.04)
    return fg._save(fig, "11_hubble_vs_ground", svg_dpi=90)
