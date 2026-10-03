"""The pipeline: from calibrated images to honest pictures, and to what their colours mean.

In order:

1. **Reference points** (``reference``): the Sun's colours through the SDSS
   filters (checked against Willmer 2018), which fix what "white" means in our
   pictures; the colours of every kind of star, to read galaxy colours against.
2. **The eye's colour matrix** (``eye_matrix``): a linear map from g, r, i
   fluxes to what the human eye sees (CIE XYZ), learned from 129 real galaxy
   spectra and tested on spectra it never saw.
3. **Validation**: do we read the pixels the way the survey does
   (``validate_viewer_jpeg``)? Does stitching tiles change them
   (``validate_stitching``)? Are the images photometrically right
   (``validate_star_photometry``, against SDSS's catalogue and against Gaia)?
4. **Each galaxy** (``analyse_galaxy``): the composite, integrated colours
   (checked against Brown et al. 2014), colour profiles, the colour an eye
   would see, and whether the eye could see colour at all.

Published values enter only in comparisons, after our own numbers exist.
"""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass, field, is_dataclass
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from PIL import Image

import color_physics as cp
import composite as cm
import imaging
from config import (
    BANDS,
    BROWN_APERTURES,
    CONTEXT_SCENES,
    DATA_DIR,
    EXTINCTION_PER_EBV_SFD,
    LUPTON_Q,
    LUPTON_SATURATION,
    MESOPIC_LOWER_CD_M2,
    MESOPIC_UPPER_CD_M2,
    PUBLISHED,
    RESULTS_PATH,
    VIEWER_SDSS_RGB,
    Galaxy,
    central_tile,
    galaxy_by_key,
)

RNG_SEED = 20040201  # Lupton et al. 2004, PASP 116, 133 -- February 2004


# --------------------------------------------------------------------------
# 1. Reference points: the Sun and the stars
# --------------------------------------------------------------------------

def _colors_of(lam, flam) -> dict[str, float]:
    f = {b: cp.band_flux(lam, flam, f"sdss_{b}") for b in BANDS}
    return {"g-r": -2.5 * np.log10(f["g"] / f["r"]), "r-i": -2.5 * np.log10(f["r"] / f["i"])}


def _eye_color(xyz, sun_white) -> dict:
    """Everything we report about how a light looks: chromaticity, colour
    temperature, and a displayable swatch (sunlight = white)."""
    x, y = cp.xyz_to_xy(xyz)
    cct, duv = cp.correlated_color_temperature(xyz)
    _, hexs, in_gamut = cp.swatch(xyz, adapt_from=sun_white)
    _, hex_d65, _ = cp.swatch(xyz)
    return {"x": float(x), "y": float(y), "cct_K": cct, "duv": duv, "hex": hexs,
            "hex_d65_white": hex_d65, "in_srgb_gamut": in_gamut}


def reference() -> dict:
    """The Sun, the white balance, the survey viewer's white balance, the stars."""
    sun_lam, sun_flam = cp.load_sun()
    sun_white = cp.sun_white_xyz()
    mags = cp.sun_absolute_magnitudes()
    pub = PUBLISHED["sun_abs_mag_ab"]["values"]
    weights = cp.solar_white_weights()
    viewer_scales = {b: VIEWER_SDSS_RGB["scales"][b][1] for b in BANDS}
    viewer_weights = {b: viewer_scales[b] / viewer_scales["r"] for b in BANDS}

    stars = {}
    for sp, (lam, flam) in cp.pickles_stars().items():
        stars[sp] = {**_colors_of(lam, flam), **_eye_color(cp.spectrum_xyz(lam, flam), sun_white)}
    vega_lam, vega_flam = cp.load_vega()
    stars["Vega (A0V, CALSPEC)"] = {**_colors_of(vega_lam, vega_flam),
                                    **_eye_color(cp.spectrum_xyz(vega_lam, vega_flam), sun_white)}

    def nearest_star(col):
        return min((k for k in stars if "Vega" not in k),
                   key=lambda k: np.hypot(stars[k]["g-r"] - col["g-r"], stars[k]["r-i"] - col["r-i"]))

    def star_distance(col):
        k = nearest_star(col)
        return float(np.hypot(stars[k]["g-r"] - col["g-r"], stars[k]["r-i"] - col["r-i"]))

    viewer_neutral = cm.neutral_color_of(viewer_weights)
    flat_neutral = cm.neutral_color_of({b: 1.0 for b in BANDS})
    sun_xyz = cp.spectrum_xyz(sun_lam, sun_flam)
    # Illuminance of sunlight at 1 AU and its luminance on a white page
    e_sun = cp.illuminance_lux(sun_lam, sun_flam)
    return {
        "sun": {
            "abs_mag_ab": mags,
            "abs_mag_ab_published": pub,
            "published_reference": PUBLISHED["sun_abs_mag_ab"]["reference"],
            "abs_mag_diff": {b: mags[b] - pub[b] for b in pub},
            "colors": _colors_of(sun_lam, sun_flam),
            "eye": _eye_color(sun_xyz, sun_white),
            "eye_d65_xy": list(map(float, cp.xyz_to_xy(sun_xyz))),
            "illuminance_1au_lux": e_sun,
        },
        "white_balance": {
            "ours": {"weights": weights, "neutral_color": cm.neutral_color_of(weights),
                     "neutral_star": "G2V (the Sun)"},
            "viewer": {"weights": viewer_weights, "neutral_color": viewer_neutral,
                       "neutral_star": nearest_star(viewer_neutral),
                       "neutral_star_distance_mag": star_distance(viewer_neutral)},
            "flat": {"weights": {b: 1.0 for b in BANDS}, "neutral_color": flat_neutral,
                     "neutral_star": nearest_star(flat_neutral),
                     "neutral_star_distance_mag": star_distance(flat_neutral)},
        },
        "stars": stars,
    }


# --------------------------------------------------------------------------
# 2. From three filters to the eye
# --------------------------------------------------------------------------

@dataclass
class EyeModel:
    color: cp.ColorMatrix
    #: fluxes (maggies) @ scotopic -> the rod-weighted integral (V' instead of ybar)
    scotopic: np.ndarray
    scotopic_loo_frac: np.ndarray
    summary: dict


def eye_matrix() -> EyeModel:
    """Learn fluxes -> XYZ from the Brown et al. (2014) galaxy spectra.

    For each of the 129 spectra: synthetic g, r, i fluxes through the SDSS
    curves, and XYZ through the CIE functions. The fit is judged by
    leave-one-out prediction, and then on a harder, independent test: 14 kinds
    of star (Pickles 1998), whose spectra are unlike any galaxy's.
    """
    brown = cp.load_brown()
    names, F, XYZ, VP, src_frac = [], [], [], [], []
    for name, (lam, flam, src) in sorted(brown.spectra.items()):
        names.append(name)
        F.append([cp.band_flux(lam, flam, f"sdss_{b}") for b in BANDS])
        XYZ.append(cp.spectrum_xyz(lam, flam))
        VP.append(cp.spectrum_scotopic(lam, flam))
        vis = (lam > 3800) & (lam < 7800)
        src_frac.append(float(np.mean(src[vis] == 1)))
    F, XYZ, VP = np.array(F), np.array(XYZ), np.array(VP)
    model = cp.fit_color_matrix(F, XYZ, names, bands=BANDS)

    # rods: V' integral from the fluxes, same normalisation (per unit Y)
    Fn, Vn = F / XYZ[:, 1:2], VP / XYZ[:, 1]
    scot, *_ = np.linalg.lstsq(Fn, Vn, rcond=None)
    loo = np.empty(len(Fn))
    for k in range(len(Fn)):
        keep = np.arange(len(Fn)) != k
        s_k, *_ = np.linalg.lstsq(Fn[keep], Vn[keep], rcond=None)
        loo[k] = Fn[k] @ s_k / Vn[k] - 1.0

    # independent test: stars
    star_err = {}
    for sp, (lam, flam) in cp.pickles_stars().items():
        f = np.array([cp.band_flux(lam, flam, f"sdss_{b}") for b in BANDS])
        star_err[sp] = float(cp.delta_upvp(model.xyz(f), cp.spectrum_xyz(lam, flam)))

    worst = int(np.argmax(model.loo_duv))
    summary = {
        "n_training_spectra": len(names),
        "median_optical_fraction_measured": float(np.median(src_frac)),
        "loo_duv_median": float(np.median(model.loo_duv)),
        "loo_duv_p90": float(np.percentile(model.loo_duv, 90)),
        "loo_duv_max": float(model.loo_duv.max()),
        "loo_duv_worst_galaxy": names[worst],
        "loo_duv_worst_morphology": brown.morphology.get(names[worst], ""),
        "stars_duv": star_err,
        "stars_duv_median": float(np.median(list(star_err.values()))),
        "scotopic_loo_rms_frac": float(np.sqrt(np.mean(loo**2))),
        "luminance_loo_rms_frac": float(np.sqrt(np.mean(model.loo_dY**2))),
        "luminance_loo_max_frac": float(np.max(np.abs(model.loo_dY))),
        "matrix": model.matrix.tolist(),
        "jnd_duv": 0.004,
    }
    return EyeModel(model, scot, loo, summary)


def reference_xyz_per_intensity(model: EyeModel) -> float:
    """Y of a source with the Sun's colours and 1 nanomaggy/arcsec^2 in r.

    This puts the eye's luminance on the same scale as the composite's
    intensity, so the true-colour picture can share the composite's stretch."""
    lam, flam = cp.load_sun()
    f = np.array([cp.band_flux(lam, flam, f"sdss_{b}") for b in BANDS])
    f = f / f[BANDS.index("r")] * 1e-9  # scale so r = 1 nanomaggy
    return float(model.color.xyz(f)[1])


# --------------------------------------------------------------------------
# 3. Validation
# --------------------------------------------------------------------------

def estimate_jpeg_quality(path) -> int:
    """The libjpeg quality setting that reproduces a file's luminance
    quantisation table (the standard IJG scaling of the Annex K table)."""
    std = np.array([16, 11, 10, 16, 24, 40, 51, 61, 12, 12, 14, 19, 26, 58, 60, 55,
                    14, 13, 16, 24, 40, 57, 69, 56, 14, 17, 22, 29, 51, 87, 80, 62,
                    18, 22, 37, 56, 68, 109, 103, 77, 24, 35, 55, 64, 81, 104, 113, 92,
                    49, 64, 78, 87, 103, 121, 120, 101, 72, 92, 95, 98, 112, 100, 103, 99])
    table = np.array(sorted(Image.open(path).quantization.items())[0][1])
    best, best_err = 75, np.inf
    for q in range(1, 101):
        scale = 5000 / q if q < 50 else 200 - 2 * q
        t = np.clip(np.floor((std * scale + 50) / 100), 1, 255)
        err = np.abs(np.sort(t) - np.sort(table)).sum()  # zigzag order differs: compare sorted
        if err < best_err:
            best, best_err = q, err
    return best


def validate_viewer_jpeg(g: Galaxy) -> dict:
    """Render the central tile with the viewer's own recipe and compare with the
    JPEG the viewer served for exactly the same pixels.

    Two comparisons: our 8-bit rendering straight against their JPEG, and our
    rendering after passing it through JPEG compression at the viewer's quality
    setting -- if the recipe is identical, that second difference should be
    only the tiny variation between JPEG encoders.
    """
    imgs, _, jpg = imaging.load_tile(g)
    t = central_tile(g)
    q = estimate_jpeg_quality(DATA_DIR / f"cutouts/{g.key}/tile_{t.iy}_{t.ix}_viewer.jpg")
    ours = cm.to_uint8(cm.viewer_sdss_rgb({b: imgs[b] for b in BANDS}, **VIEWER_SDSS_RGB),
                       flip=False, mode="floor")
    # compress in display orientation (north up), as the viewer did, then back
    ours_jpeg = cm.jpeg_roundtrip(ours[::-1], q)[::-1]
    theirs = jpg.astype(float)

    def stats(a):
        d = a.astype(float) - theirs
        return {"mean_diff": float(d.mean()), "median_abs_diff": float(np.median(np.abs(d))),
                "p99_abs_diff": float(np.percentile(np.abs(d), 99)),
                "frac_within_2": float(np.mean(np.abs(d) <= 2)),
                "corr": [float(np.corrcoef(a[..., k].ravel(), theirs[..., k].ravel())[0, 1]) for k in range(3)]}

    # what compression alone does to our image, for scale
    self_d = ours_jpeg.astype(float) - ours.astype(float)
    summary = {"galaxy": g.key, "tile": t.filename, "jpeg_quality": q,
               "raw": stats(ours), "compressed": stats(ours_jpeg),
               "compression_alone_median_abs": float(np.median(np.abs(self_d))),
               "compression_alone_p99_abs": float(np.percentile(np.abs(self_d), 99))}
    pixels = {"diff_raw": (ours.astype(np.int16) - jpg.astype(np.int16)).ravel(),
              "diff_compressed": (ours_jpeg.astype(np.int16) - jpg.astype(np.int16)).ravel()}
    return summary, pixels


def _subpixel_shift(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """(dy, dx) that best aligns a with b: phase correlation plus a parabola
    through the peak for the sub-pixel part."""
    F = np.fft.fft2(a) * np.conj(np.fft.fft2(b))
    cc = np.fft.ifft2(F / np.maximum(np.abs(F), 1e-30)).real
    iy, ix = np.unravel_index(np.argmax(cc), cc.shape)
    ny, nx = cc.shape

    def refine(m1, c0, p1):
        den = m1 - 2 * c0 + p1
        return 0.5 * (m1 - p1) / den if den != 0 else 0.0

    dy = refine(cc[(iy - 1) % ny, ix], cc[iy, ix], cc[(iy + 1) % ny, ix])
    dx = refine(cc[iy, (ix - 1) % nx], cc[iy, ix], cc[iy, (ix + 1) % nx])
    wrap = lambda i, n: i - n if i > n // 2 else i
    return float(wrap(iy, ny) + dy), float(wrap(ix, nx) + dx)


def validate_stitching(key: str) -> dict:
    """Our stitched mosaic against the same frame fetched as one piece.

    Three questions: do the two line up (sub-pixel shift by phase
    correlation)? How much does a typical pixel change, relative to the
    noise? And is the light conserved?  The extra interpolation smooths the
    noise a little and rounds off the sharpest star peaks by a few per cent,
    so the typical change and the bright-pixel change are reported separately.
    """
    g = galaxy_by_key(key)
    frame = imaging.load_frame(g)
    single = imaging.load_single_cutout(g)
    out = {"galaxy": key}
    n = frame.shape[0]
    c = slice(n // 5, 4 * n // 5)
    for b in BANDS:
        a, s = frame.native[b], single[b]
        d = a - s
        noise = 1.4826 * np.median(np.abs(s - np.median(s)))
        bright = s > 20 * noise
        clipped = [np.clip(x[c, c] - np.median(x[c, c]), -5 * noise, 50 * noise) for x in (a, s)]
        dy, dx = _subpixel_shift(*clipped)
        out[b] = {
            "shift_px": [dy, dx],
            "typical_diff_over_noise": float(1.4826 * np.median(np.abs(d - np.median(d))) / noise),
            "frac_pixels_beyond_5_noise": float(np.mean(np.abs(d) > 5 * noise)),
            "bright_pixels_median_abs_frac": float(np.median(np.abs(d[bright] / s[bright]))),
            "bright_pixels_p99_abs_frac": float(np.percentile(np.abs(d[bright] / s[bright]), 99)),
            "total_flux_ratio": float(a.sum() / s.sum()),
        }
    return out


def _robust(x) -> dict:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    med = float(np.median(x))
    return {"n": int(x.size), "median": med, "robust_sigma": float(1.4826 * np.median(np.abs(x - med)))}


def star_photometry_rows(fr: imaging.Frame) -> list[dict]:
    """Aperture photometry of the field's stars on our images, next to two
    catalogues' magnitudes for the same stars.

    SDSS's PSF magnitudes come from the same raw data processed by the survey;
    Gaia's synthetic SDSS magnitudes come from a different telescope entirely
    (GSPC; only stars inside Gaia's validated colour-magnitude range).
    """
    g = fr.galaxy
    sdss = imaging.load_sdss_stars(g)
    gaia = imaging.load_gaia(g)
    ok = np.ones(len(gaia), bool)
    for b in BANDS:
        ok &= gaia[f"{b}_sdss_flag"].fillna(0).to_numpy() == 1
    ok &= gaia["r_sdss_mag"].between(14.5, 17.6).to_numpy()
    gx, gy = fr.world_to_pixel(gaia["ra"].to_numpy(), gaia["dec"].to_numpy())
    rows = []
    for cat, df, magcols in (("sdss", sdss, {b: f"psfMag_{b}" for b in BANDS}),
                             ("gaia", gaia[ok], {b: f"{b}_sdss_mag" for b in BANDS})):
        if df.empty:
            continue
        ph = imaging.star_photometry(fr, df["ra"].to_numpy(), df["dec"].to_numpy())
        for k, (_, star) in enumerate(df.iterrows()):
            p = ph.iloc[k]
            if not p["inside"]:
                continue
            dist = np.hypot(gx - p["x"], gy - p["y"]) * fr.pixscale
            dist = dist[dist > 0.5]  # the star itself
            rec = {"galaxy": g.key, "catalog": cat, "neighbour_arcsec": float(dist.min()) if dist.size else np.inf,
                   "bkg_r": float(p["bkg_r"])}
            for b in BANDS:
                rec[f"cat_{b}"] = float(star[magcols[b]])
                rec[f"ours_{b}"] = float(p[f"mag_{b}"])
            rows.append(rec)
    return rows


def summarize_star_photometry(rows: list[dict], sky_noise_r: float) -> dict:
    """Median offsets and robust scatter, from clean stars only.

    Clean: no other Gaia source within 15" (confusion), and a local background
    within 2 pixel-noise units of zero -- a star sitting on bright galaxy light
    has an annulus that cannot see the true sky.
    """
    df = pd.DataFrame(rows)
    clean = (df["neighbour_arcsec"] > 15) & (df["bkg_r"].abs() < 2 * sky_noise_r)
    for b in BANDS:
        clean &= np.isfinite(df[f"ours_{b}"])
    out = {"n_measured": int(len(df)), "n_clean": int(clean.sum()),
           "clean_background_limit_nmgy_arcsec2": 2 * sky_noise_r,
           "table": df[clean].to_dict("list")}
    for cat in ("sdss", "gaia"):
        sel = clean & (df["catalog"] == cat)
        out[cat] = {b: _robust(df.loc[sel, f"ours_{b}"] - df.loc[sel, f"cat_{b}"]) for b in BANDS}
        for c1, c2 in (("g", "r"), ("r", "i")):
            ours = df.loc[sel, f"ours_{c1}"] - df.loc[sel, f"ours_{c2}"]
            cat_c = df.loc[sel, f"cat_{c1}"] - df.loc[sel, f"cat_{c2}"]
            out[cat][f"{c1}-{c2}"] = _robust(ours - cat_c)
    return out


# --------------------------------------------------------------------------
# 4. One galaxy
# --------------------------------------------------------------------------

def dust_ebv_sfd(g: Galaxy) -> float:
    """E(B-V) from the original SFD map at the galaxy's position (IRSA)."""
    text = (DATA_DIR / f"catalogs/{g.key}_dust.xml").read_text()
    m = re.search(r"<refPixelValueSFD>\s*([0-9.]+)", text)
    return float(m.group(1))


@dataclass
class GalaxyResult:
    key: str
    name: str
    common: str
    kind: str
    rc3_name: str
    d25_arcmin: float
    axis_ratio: float
    pa_deg: float
    ebv_sfd: float
    n_foreground_stars: int
    masked_fraction_in_d25: float
    sky: dict
    #: total AB magnitudes inside the D25 ellipse (stars masked, not corrected)
    mags_d25: dict
    colors_d25: dict
    colors_d25_dereddened: dict
    #: same for companions listed in config
    companions: dict
    #: Brown et al. (2014): their aperture photometry, and colours of their spectrum
    brown: dict
    #: what an eye would see: from our image and (if available) from the spectrum
    eye_from_image: dict
    eye_from_spectrum: dict
    #: colour profile along the disc
    profile: dict
    #: surface brightness and luminance
    luminance: dict
    composite: dict = field(default_factory=dict)


def _galaxy_masks(fr: imaging.Frame):
    """(stars mask, main-galaxy D25 mask, {companion: D25 mask}, anything-galaxy mask)."""
    g = fr.galaxy
    stars = imaging.foreground_stars(imaging.load_gaia(g))
    smask = imaging.star_mask(fr, stars)
    main = imaging.rc3_entry(g)
    d25 = imaging.d25_mask(fr, main)
    comp = {}
    for c in g.companions:
        comp[c] = imaging.d25_mask(fr, imaging.rc3_entry(g, c))
    anything = np.zeros_like(d25)
    for _, row in imaging.load_rc3(g.key).iterrows():
        if np.isfinite(row["d25_arcmin"]):
            anything |= imaging.d25_mask(fr, row, scale=1.3)
    return stars, smask, main, d25, comp, anything


def _mags(fr, mask, images: dict | None = None) -> dict:
    images = images if images is not None else fr.sb
    return {b: float(cp.nmgy_to_mag(imaging.aperture_flux(images[b], mask, fr.pixscale))) for b in BANDS}


def _colors(m: dict) -> dict:
    return {"g-r": m["g"] - m["r"], "r-i": m["r"] - m["i"]}


#: A ring's colour is reported only where the galaxy is at least this many
#: times brighter (in r) than the residual sky level measured around it: there
#: an error in the sky as large as the residual itself would move the ring's
#: light by no more than 10%.
PROFILE_MIN_SKY_RATIO = 10.0
#: ...and in any case brighter than this in r (nanomaggies/arcsec^2; 24.0
#: mag/arcsec^2). Where the residual happens to measure near zero, the sky
#: level is still uncertain at the ~0.02 level seen in other frames.
PROFILE_MIN_SB_R = 0.25


def color_profile(fr: imaging.Frame, entry, keep: np.ndarray, *, n_bins: int = 24,
                  rmax_frac: float = 1.2, ebv: float = 0.0, sky: dict | None = None) -> dict:
    """g-r and r-i in elliptical annuli, out to 1.2 x the D25 radius.

    Fluxes are summed in each annulus first and the colour taken from the sums,
    so faint noisy pixels average out instead of producing wild pixel colours.
    The residual sky level measured away from the galaxy is subtracted first,
    and rings too faint for that subtraction to be trusted are left out (see
    ``PROFILE_MIN_SKY_RATIO``): the residual is itself red, and would otherwise
    redden the faint outskirts. Errors come from resampling the pixels in each
    annulus.
    """
    sky = sky or {b: {"median": 0.0} for b in BANDS}
    offset = {b: sky[b]["median"] for b in BANDS}
    floor_r = max(PROFILE_MIN_SKY_RATIO * abs(offset["r"]), PROFILE_MIN_SB_R)
    q = entry["axis_ratio"] if np.isfinite(entry["axis_ratio"]) else 1.0
    r_ell = imaging.elliptical_radius(fr, entry["ra"], entry["dec"], q, entry["pa_deg"])
    r25 = 0.5 * entry["d25_arcmin"] * 60.0
    edges = np.linspace(0.0, rmax_frac, n_bins + 1)
    rng = np.random.default_rng(RNG_SEED)
    red = {b: EXTINCTION_PER_EBV_SFD[b] * ebv for b in BANDS}
    out = {"r_over_r25": [], "g-r": [], "r-i": [], "g-r_err": [], "r-i_err": [], "mu_r": []}
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = keep & (r_ell >= lo * r25) & (r_ell < hi * r25)
        if sel.sum() < 30:
            continue
        vals = {b: fr.sb[b][sel] - offset[b] for b in BANDS}
        tot = {b: vals[b].sum() for b in BANDS}
        if min(tot.values()) <= 0 or tot["r"] / sel.sum() < floor_r:
            continue
        col = {"g-r": -2.5 * np.log10(tot["g"] / tot["r"]) - (red["g"] - red["r"]),
               "r-i": -2.5 * np.log10(tot["r"] / tot["i"]) - (red["r"] - red["i"])}
        boots = {"g-r": [], "r-i": []}
        n = sel.sum()
        for _ in range(60):
            idx = rng.integers(0, n, n)
            s = {b: vals[b][idx].sum() for b in BANDS}
            if min(s.values()) <= 0:
                continue
            boots["g-r"].append(-2.5 * np.log10(s["g"] / s["r"]))
            boots["r-i"].append(-2.5 * np.log10(s["r"] / s["i"]))
        out["r_over_r25"].append(0.5 * (lo + hi))
        out["mu_r"].append(float(cp.mag_arcsec2(np.median(vals["r"]))))
        for c in ("g-r", "r-i"):
            out[c].append(float(col[c]))
            out[f"{c}_err"].append(float(np.std(boots[c])) if len(boots[c]) > 5 else np.nan)
    return out


def analyse_galaxy(g: Galaxy, eye: EyeModel, ref: dict) -> tuple[GalaxyResult, dict]:
    """Everything measured on one galaxy. Returns the numbers and the images
    (as arrays) the figures need."""
    fr = imaging.load_frame(g)
    stars, smask, main, d25, comp, anything = _galaxy_masks(fr)
    ebv = dust_ebv_sfd(g)
    red = {b: EXTINCTION_PER_EBV_SFD[b] * ebv for b in BANDS}

    # sky: away from every catalogued galaxy and every star
    sky = imaging.sky_statistics(fr, anything | imaging.star_mask(fr, stars, scale=1.5))
    # aperture photometry on images whose foreground stars are replaced by the
    # surrounding galaxy light (dropping them would drop the galaxy under them too)
    filled = imaging.fill_star_masks(fr, stars, smask)
    keep_main = d25.copy()
    for m in comp.values():
        keep_main &= ~m
    mags = _mags(fr, keep_main, filled)
    colors = _colors(mags)
    dered = {"g-r": colors["g-r"] - (red["g"] - red["r"]), "r-i": colors["r-i"] - (red["r"] - red["i"])}
    companions = {}
    for c, m in comp.items():
        cm_ = _mags(fr, m, filled)
        companions[c] = {"mags": cm_, "colors": _colors(cm_)}

    # Brown et al.: their photometry, and synthetic colours of their spectrum
    brown = {}
    if g.brown:
        atlas = cp.load_brown()
        for bname in g.brown:
            lam, flam, _ = atlas.spectra[bname]
            ph = atlas.photometry[bname]
            length, width, pa = BROWN_APERTURES[bname]
            ra_b, dec_b = atlas.position[bname]
            rect = imaging.rectangle_mask(fr, ra_b, dec_b, length, width, pa)
            ours = _mags(fr, rect, filled)
            dropped = _mags(fr, rect & ~smask)
            brown[bname] = {
                "aperture": {"length_arcsec": length, "width_arcsec": width, "pa_deg": pa,
                             "ra": ra_b, "dec": dec_b,
                             "inside_frame": bool(rect.sum() * fr.pixscale**2 > 0.98 * length * width)},
                "published_mags": {b: ph[b] for b in BANDS},
                "published_colors": {"g-r": ph["g"] - ph["r"], "r-i": ph["r"] - ph["i"]},
                "our_mags": ours,
                "our_colors": _colors(ours),
                "mag_diff": {b: ours[b] - ph[b] for b in BANDS},
                "mag_diff_if_star_pixels_dropped": {b: dropped[b] - ph[b] for b in BANDS},
                "star_masked_fraction": float((rect & smask).sum() / max(rect.sum(), 1)),
                "spectrum_colors": _colors_of(lam, flam),
                "eye": _eye_color(cp.spectrum_xyz(lam, flam), ref["sun_white_xyz"]),
            }

    # the eye: integrated light inside D25
    flux = np.array([imaging.aperture_flux(filled[b], keep_main, fr.pixscale) for b in BANDS]) * 1e-9
    eye_img = _eye_color(eye.color.xyz(flux), ref["sun_white_xyz"])
    eye_spec = {}
    if g.brown:
        eye_spec = brown[g.brown[0]]["eye"]
        eye_img["duv_vs_spectrum"] = float(cp.delta_upvp(
            eye.color.xyz(flux), cp.xy_to_xyz(eye_spec["x"], eye_spec["y"])))

    # luminance of every pixel
    fluxes = np.stack([fr.sb[b] for b in BANDS], axis=-1) * 1e-9  # maggies / arcsec^2
    Y = eye.color.xyz(fluxes)[..., 1]
    lum = cp.luminance_cd_m2(Y)
    # the rods' view: the same fluxes through the scotopic map fitted in eye_matrix
    scot = cp.scotopic_luminance_cd_m2(fluxes @ eye.scotopic)
    # the seeing-limited peak: brightest 3x3 pixel average near the centre
    x0, y0 = fr.world_to_pixel(main["ra"], main["dec"])
    ix, iy = int(round(float(x0))), int(round(float(y0)))
    half = int(round(10.0 / fr.pixscale))
    box = lum[max(0, iy - half):iy + half + 1, max(0, ix - half):ix + half + 1]
    from scipy.ndimage import uniform_filter
    peak = float(uniform_filter(box, size=3).max())
    in_d25 = lum[keep_main & ~smask]
    luminance = {
        "peak_cd_m2": peak,
        "peak_mu_r": float(cp.mag_arcsec2(uniform_filter(fr.sb["r"], 3)[max(0, iy - half):iy + half + 1,
                                                                      max(0, ix - half):ix + half + 1].max())),
        "median_in_d25_cd_m2": float(np.median(in_d25)),
        "p90_in_d25_cd_m2": float(np.percentile(in_d25, 90)),
        "frac_d25_above_mesopic_lower": float(np.mean(in_d25 > MESOPIC_LOWER_CD_M2)),
        "frac_d25_above_photopic": float(np.mean(in_d25 > MESOPIC_UPPER_CD_M2)),
        "frac_d25_brighter_than_dark_sky": float(np.mean(in_d25 > _dark_sky_luminance())),
        "scotopic_peak_cd_m2": float(uniform_filter(scot[max(0, iy - half):iy + half + 1,
                                                        max(0, ix - half):ix + half + 1], 3).max()),
    }

    keep_profile = ~smask
    for m in comp.values():
        keep_profile &= ~m
    profile = color_profile(fr, main, keep_profile, ebv=ebv, sky=sky)

    res = GalaxyResult(
        key=g.key, name=g.name, common=g.common, kind=g.kind, rc3_name=str(main["name"]),
        d25_arcmin=float(main["d25_arcmin"]), axis_ratio=float(main["axis_ratio"]),
        pa_deg=float(main["pa_deg"]) if np.isfinite(main["pa_deg"]) else float("nan"),
        ebv_sfd=ebv, n_foreground_stars=int(len(stars)),
        masked_fraction_in_d25=float(np.mean(smask[d25])),
        sky=sky, mags_d25=mags, colors_d25=colors, colors_d25_dereddened=dered,
        companions=companions, brown=brown, eye_from_image=eye_img, eye_from_spectrum=eye_spec,
        profile=profile, luminance=luminance,
    )
    art = {"frame": fr, "star_mask": smask, "d25": d25, "companions": comp, "luminance": lum,
           "scotopic_luminance": scot, "rc3": main}
    return res, art


def _dark_sky_luminance() -> float:
    """The natural moonless sky, 22.0 V mag/arcsec^2 (Falchi et al. 2016), which
    they give as 174 micro-cd/m^2."""
    return 1.74e-4


def context_luminances() -> dict[str, float]:
    out = {}
    for name, spec in CONTEXT_SCENES.items():
        if "illuminance_lux" in spec:
            out[name] = cp.paper_luminance(spec["illuminance_lux"])
        else:
            out[name] = _dark_sky_luminance()
    return out


# --------------------------------------------------------------------------
# Supporting measurements for the figures
# --------------------------------------------------------------------------

def eye_visibility_per_band() -> dict[str, float]:
    """How visible each filter's light is: the eye's brightness response ybar,
    averaged over the filter's photon-weighted response (1 = peak eye response)."""
    cl, cmf = cp.load_cmf()
    ybar = cmf[:, 1] / cmf[:, 1].max()
    out = {}
    for b in BANDS:
        fl, fr = cp.load_filter(f"sdss_{b}")
        grid = np.linspace(fl.min(), fl.max(), 4000)
        resp = np.interp(grid, fl, fr) * grid
        out[b] = float(np.trapezoid(resp * np.interp(grid, cl, ybar, left=0, right=0), grid)
                       / np.trapezoid(resp, grid))
    return out


def display_color_profile(fr: imaging.Frame, entry, keep: np.ndarray, weights: dict, noise_r: float,
                          n_bins: int = 20, rmax_frac: float = 1.0) -> dict:
    """The g-r colour each stretch *displays*, ring by ring, against the truth.

    For every pixel we invert the displayed red/green/blue values back into
    the colour they imply (``composite.displayed_color``), then take the
    ring's median. A colour-preserving stretch returns the measured colour; a
    per-channel one drifts towards whatever colour the stretch treats as grey.
    """
    R, G, B = cm.channels(fr.sb, weights)
    soft = cm.softening_for(LUPTON_Q, LUPTON_SATURATION)
    shown = {
        "lupton": cm.lupton_rgb(R, G, B, softening=soft, saturation=LUPTON_SATURATION),
        "per_channel": cm.per_channel_asinh(R, G, B, softening=soft, saturation=LUPTON_SATURATION),
    }
    q = entry["axis_ratio"] if np.isfinite(entry["axis_ratio"]) else 1.0
    r_ell = imaging.elliptical_radius(fr, entry["ra"], entry["dec"], q, entry["pa_deg"])
    r25 = 0.5 * entry["d25_arcmin"] * 60.0
    edges = np.linspace(0, rmax_frac, n_bins + 1)
    out = {"r_over_r25": [], "true": [], "lupton": [], "per_channel": []}
    disp = {k: cm.displayed_color(v, weights)["g-r"] for k, v in shown.items()}
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = keep & (r_ell >= lo * r25) & (r_ell < hi * r25)
        if sel.sum() < 30:
            continue
        tot_g, tot_r = fr.sb["g"][sel].sum(), fr.sb["r"][sel].sum()
        if tot_g <= 0 or tot_r <= 0:
            continue
        # the displayed colour of the ring: median over pixels bright enough to have one
        bright = sel & (fr.sb["r"] > 5 * noise_r) & (fr.sb["g"] > 0)
        out["r_over_r25"].append(0.5 * (lo + hi))
        out["true"].append(float(np.median(-2.5 * np.log10(fr.sb["g"][bright] / fr.sb["r"][bright])))
                           if bright.sum() > 10 else np.nan)
        for k in ("lupton", "per_channel"):
            v = disp[k][bright]
            v = v[np.isfinite(v)]
            out[k].append(float(np.median(v)) if v.size > 10 else np.nan)
    return out


def color_index_map(fr: imaging.Frame, smask: np.ndarray, ebv: float, *, min_snr: float = 8.0,
                    sky: dict | None = None) -> np.ndarray:
    """g-r of every pixel, lightly smoothed, dust-corrected, residual sky
    removed, NaN where too noisy."""
    from scipy.ndimage import gaussian_filter

    sky = sky or {b: {"median": 0.0} for b in BANDS}
    sig = 1.0  # pixels: about the seeing
    g_s = gaussian_filter(fr.sb["g"] - sky["g"]["median"], sig)
    r_s = gaussian_filter(fr.sb["r"] - sky["r"]["median"], sig)
    # noise of the smoothed images, from the frame's own sky statistics
    noise_g = np.std(fr.sb["g"][~smask & (fr.sb["r"] < np.percentile(fr.sb["r"], 40))]) / (2 * np.sqrt(np.pi) * sig)
    noise_r = np.std(fr.sb["r"][~smask & (fr.sb["r"] < np.percentile(fr.sb["r"], 40))]) / (2 * np.sqrt(np.pi) * sig)
    ok = (g_s > min_snr * noise_g) & (r_s > min_snr * noise_r) & ~smask
    with np.errstate(divide="ignore", invalid="ignore"):
        gr = -2.5 * np.log10(g_s / r_s) - (EXTINCTION_PER_EBV_SFD["g"] - EXTINCTION_PER_EBV_SFD["r"]) * ebv
    return np.where(ok, gr, np.nan)


def pipeline_caveat(ref: dict) -> dict:
    """How much of M51's light survives in the Legacy Surveys DR9-north images.

    The same galaxy, the same sky, a different survey pipeline. We reproject
    our SDSS g, r, i images onto the DR9 cutout's grid and compare the light
    inside M51's D25 ellipse. SDSS has no z image in the viewer, so the z-band
    light we *expect* comes from SDSS i and M51's i - z colour in Brown et al.'s
    photometry. DR9-north's g, r and z come from different cameras and filters
    than SDSS's (BASS g, r; MzLS z), which differ by a few per cent -- far less
    than the effect being measured.
    """
    from reproject import reproject_interp

    g = galaxy_by_key("m51")
    path = DATA_DIR / "multiwavelength" / "m51_lsdr9north.fits"
    imgs, w = imaging.read_cutout(path)
    shape = imgs["g"].shape
    fr = imaging.load_frame(g)
    sdss = {b: reproject_interp((fr.sb[b], fr.wcs), w, shape_out=shape, order="bicubic")[0] for b in BANDS}
    pixscale = float(np.sqrt(np.abs(np.linalg.det(w.pixel_scale_matrix)))) * 3600.0
    from dataclasses import replace
    g_dr9 = replace(g, pixscale=pixscale, fov_arcmin=shape[0] * pixscale / 60.0)
    stub = imaging.Frame(g_dr9, w, {b: imgs.get(b, sdss["r"]) for b in BANDS}, sdss)
    main = imaging.rc3_entry(g)
    ap = imaging.d25_mask(stub, main) & np.all([np.isfinite(sdss[b]) for b in BANDS], axis=0)
    native_dr9 = 0.262  # arcsec: DR9 pixel values are nanomaggies per native pixel
    dr9_flux = {b: float(np.nansum(imgs[b][ap]) * (pixscale / native_dr9) ** 2) for b in ("g", "r", "z")}
    sdss_flux = {b: float(np.nansum(sdss[b][ap]) * pixscale**2) for b in BANDS}
    ph = cp.load_brown().photometry["NGC_5194"]
    i_minus_z = ph["i"] - ph["z"]
    expected_z = sdss_flux["i"] * 10 ** (0.4 * i_minus_z)
    return {
        "aperture": "RC3 D25 ellipse of NGC 5194",
        "dr9_flux_nmgy": dr9_flux, "sdss_flux_nmgy": sdss_flux,
        "i_minus_z_brown": i_minus_z,
        "retained": {"g": dr9_flux["g"] / sdss_flux["g"], "r": dr9_flux["r"] / sdss_flux["r"],
                     "z": dr9_flux["z"] / expected_z},
        "images": {"dr9": imgs, "wcs": w},
    }


# --------------------------------------------------------------------------
# Composites for the gallery
# --------------------------------------------------------------------------

def make_composites(fr: imaging.Frame, ref: dict, eye: EyeModel) -> dict[str, np.ndarray]:
    weights = ref["white_balance"]["ours"]["weights"]
    std = cm.composite(fr.sb, weights, Q=LUPTON_Q, saturation=LUPTON_SATURATION)
    true, diag = cm.true_color_rgb(fr.sb, eye.color, white_xyz=ref["sun_white_xyz"],
                                   reference_xyz_per_intensity=ref["y_per_intensity"],
                                   Q=LUPTON_Q, saturation=LUPTON_SATURATION)
    return {"standard": std, "true": true, "true_diag": diag}


# --------------------------------------------------------------------------
# Saving
# --------------------------------------------------------------------------

def _jsonable(obj):
    if is_dataclass(obj):
        return _jsonable(asdict(obj))
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        return None if not np.isfinite(obj) else float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return _jsonable(obj.tolist())
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def save_results(**sections) -> None:
    payload = {"generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    payload.update({k: _jsonable(v) for k, v in sections.items()})
    RESULTS_PATH.write_text(json.dumps(payload, indent=1))
