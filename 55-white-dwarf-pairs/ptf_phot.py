"""Forced photometry of ZTF J1539+5027 on Palomar Transient Factory images (2009-2012).

PTF imaged the field 124 times in R (60-s exposures). No public difference images
exist, so we measure the star directly on the calibrated science images, the way
Burdge et al. used these epochs (they had difference images; we do not):

* aperture photometry at the star's Gaia position (whether or not it is detected),
  with a small aperture (0.8 x the image FWHM, which maximises the signal-to-noise of
  a faint point source) and a local sky annulus; sub-pixel aperture weights;
* every Gaia star of G = 14.5-18.5 that is isolated and not flagged variable is
  measured the same way; their instrumental-minus-catalogue magnitudes give the
  image's zero point in the same aperture, so seeing and aperture losses cancel;
* flux errors from the sky noise and the star's photon noise, then rescaled by how
  much the comparison stars of G = 19-20 scatter from image to image;
* time: shutter open (UTC-OBS) + 30 s, converted to BJD_TDB as for ZTF.

The fluxes are on the Gaia G scale (microjansky of a G-band AB-like zero point): the
colour difference between the PTF R filter and G only rescales them by a constant,
which the timing model absorbs.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning

import config as C
import fetch_data as F
from ztf_phot import PALOMAR, bjd_tdb, target_coord

OUT_PATH = C.DERIVED_DIR / "ptf_photometry.csv"
SKY_ANNULUS = (9.0, 16.0)  # px
APERTURE_FWHM = 0.8
SUBSAMPLE = 5


def read_field_stars() -> pd.DataFrame:
    # VizieR TSV: header line, units line, dashes line, then the rows
    df = pd.read_csv(F.PTF_GAIA_FIELD, sep="\t", comment="#")
    df = df.iloc[2:].reset_index(drop=True)
    for c in ("RA_ICRS", "DE_ICRS", "Gmag", "BP-RP", "RUWE", "pmRA", "pmDE"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["Source"] = df["Source"].astype(str).str.strip()
    df["VarFlag"] = df["VarFlag"].astype(str).str.strip()
    return df


def aperture_weights(shape, x0, y0, radius):
    """Fractional pixel coverage of a circle (sub-sampled SUBSAMPLE x SUBSAMPLE)."""
    ny, nx = shape
    r_box = int(np.ceil(radius)) + 2
    xi0, yi0 = int(round(x0)), int(round(y0))
    xs = np.arange(max(xi0 - r_box, 0), min(xi0 + r_box + 1, nx))
    ys = np.arange(max(yi0 - r_box, 0), min(yi0 + r_box + 1, ny))
    off = (np.arange(SUBSAMPLE) + 0.5) / SUBSAMPLE - 0.5
    gx = xs[:, None] + off[None, :]
    gy = ys[:, None] + off[None, :]
    dx2 = (gx - x0) ** 2  # (nxs, sub)
    dy2 = (gy - y0) ** 2
    inside = (dy2[:, None, :, None] + dx2[None, :, None, :]) <= radius**2  # (nys, nxs, sub, sub)
    w = inside.mean(axis=(2, 3))
    return ys, xs, w


def aperture_sum(img, x0, y0, radius, sky_in, sky_out):
    ny, nx = img.shape
    ys, xs, w = aperture_weights(img.shape, x0, y0, radius)
    cut = img[np.ix_(ys, xs)]
    good = np.isfinite(cut)
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(xx - x0, yy - y0)
    ann = (r >= sky_in) & (r <= sky_out) & np.isfinite(img)
    if ann.sum() < 30 or not good.all():
        return np.nan, np.nan, np.nan, np.nan
    sky_vals = img[ann]
    med = np.median(sky_vals)
    # one pass of 3-sigma clipping against stars in the annulus
    rms = 1.4826 * np.median(np.abs(sky_vals - med))
    keep = np.abs(sky_vals - med) < 3 * rms
    sky = np.median(sky_vals[keep])
    rms = np.std(sky_vals[keep])
    npix = w.sum()
    flux = float(np.sum(w * (cut - sky)))
    return flux, float(sky), float(rms), float(npix)


def measure_image(path, stars: pd.DataFrame) -> dict | None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        with fits.open(path) as hdul:
            img = hdul[0].data.astype(float)
            hdr = hdul[0].header
            wcs = WCS(hdr)
    gain = float(hdr.get("GAIN", 1.8))
    fwhm_px = float(hdr.get("MEDFWHM", np.nan)) / float(hdr.get("PIXSCALE", 1.01))
    if not np.isfinite(fwhm_px) or fwhm_px <= 0:
        fwhm_px = float(hdr.get("SEEING", 2.0))
    radius = float(np.clip(APERTURE_FWHM * fwhm_px, 1.5, 4.0))
    t_open = Time(hdr["UTC-OBS"], format="isot", scale="utc")
    exptime = float(hdr.get("AEXPTIME", hdr["EXPTIME"]))
    tmid = t_open + exptime / 2 * u.s
    epoch = tmid.jyear
    # comparison stars (positions propagated with proper motion)
    dt = epoch - 2016.0
    ra = stars["RA_ICRS"] + stars["pmRA"].fillna(0) * dt / 3.6e6 / np.cos(np.radians(stars["DE_ICRS"]))
    dec = stars["DE_ICRS"] + stars["pmDE"].fillna(0) * dt / 3.6e6
    sx, sy = wcs.world_to_pixel(SkyCoord(ra.to_numpy() * u.deg, dec.to_numpy() * u.deg))
    ny, nx = img.shape
    rows = []
    for i in range(len(stars)):
        if not (20 <= sx[i] <= nx - 21 and 20 <= sy[i] <= ny - 21):
            continue
        f, sky, rms, npix = aperture_sum(img, sx[i], sy[i], radius, *SKY_ANNULUS)
        rows.append((stars["Source"].iloc[i], stars["Gmag"].iloc[i], f, sky, rms, npix))
    comp = pd.DataFrame(rows, columns=["source", "G", "flux", "sky", "rms", "npix"])
    tx, ty = wcs.world_to_pixel(target_coord(epoch))
    ft, skyt, rmst, npixt = aperture_sum(img, float(tx), float(ty), radius, *SKY_ANNULUS)
    return {"image": path.name, "utc_open": hdr["UTC-OBS"], "exptime_s": exptime, "tmid_mjd_utc": tmid.mjd,
            "fwhm_px": fwhm_px, "aperture_px": radius, "gain": gain, "target_counts": ft,
            "target_sky": skyt, "target_sky_rms": rmst, "target_npix": npixt, "x": float(tx), "y": float(ty),
            "imagezpt": hdr.get("IMAGEZPT", np.nan), "phtcalfl": hdr.get("PHTCALFL", -1),
            "airmass": hdr.get("AIRMASS", np.nan), "comp": comp}


def run(force: bool = False) -> pd.DataFrame:
    if OUT_PATH.exists() and not force:
        return pd.read_csv(OUT_PATH)
    stars = read_field_stars()
    target = SkyCoord(C.TARGET_RA_2016 * u.deg, C.TARGET_DEC_2016 * u.deg)
    pos = SkyCoord(stars["RA_ICRS"].to_numpy() * u.deg, stars["DE_ICRS"].to_numpy() * u.deg)
    stars = stars[pos.separation(target).arcsec > 3].reset_index(drop=True)  # not the target itself
    # isolation: no other Gaia source within 8"
    pos = SkyCoord(stars["RA_ICRS"].to_numpy() * u.deg, stars["DE_ICRS"].to_numpy() * u.deg)
    idx, sep, _ = pos.match_to_catalog_sky(pos, nthneighbor=2)
    stars["isolated"] = sep.arcsec > 8
    exposures = pd.read_csv(F.PTF_EXPOSURES)
    results, comps = [], []
    for row in exposures.to_dict("records"):
        path = F.ptf_cutout_path(row)
        if not path.exists():
            continue
        m = measure_image(path, stars)
        if m is None:
            continue
        comp = m.pop("comp")
        comp["image"] = m["image"]
        comps.append(comp)
        results.append(m)
    df = pd.DataFrame(results)
    comp = pd.concat(comps, ignore_index=True).merge(
        stars[["Source", "isolated", "VarFlag", "RUWE", "BP-RP"]], left_on="source", right_on="Source")
    # zero point per image from bright, isolated, non-variable comparison stars
    zp_rows = []
    for img, g in comp.groupby("image"):
        ok = (g["isolated"] & (g["VarFlag"] != "VARIABLE") & (g["G"].between(14.5, 18.5))
              & (g["flux"] > 0) & (g["RUWE"] < 1.4))
        zps = g.loc[ok, "G"] + 2.5 * np.log10(g.loc[ok, "flux"])
        if len(zps) < 5:
            zp_rows.append((img, np.nan, np.nan, int(len(zps))))
            continue
        med = np.median(zps)
        mad = 1.4826 * np.median(np.abs(zps - med))
        keep = np.abs(zps - med) < 3 * max(mad, 0.01)
        zp_rows.append((img, float(np.median(zps[keep])), float(np.std(zps[keep]) / np.sqrt(keep.sum())),
                        int(keep.sum())))
    zp = pd.DataFrame(zp_rows, columns=["image", "zp_g", "zp_err", "n_zp"])
    df = df.merge(zp, on="image")
    # flux on the G AB-like scale: m = zp - 2.5 log10(counts) -> microjansky
    to_ujy = 10 ** ((23.9 - df["zp_g"]) / 2.5)
    noise_counts = np.sqrt(df["target_npix"] * df["target_sky_rms"] ** 2
                           + np.clip(df["target_counts"], 0, None) / df["gain"])
    df["flux_ujy"] = df["target_counts"] * to_ujy
    df["flux_err_white_ujy"] = noise_counts * to_ujy
    # empirical noise: scatter of faint constant comparison stars (G 19-20) between images
    comp = comp.merge(zp[["image", "zp_g"]], on="image")
    comp["flux_ujy"] = comp["flux"] * 10 ** ((23.9 - comp["zp_g"]) / 2.5)
    faint = comp[comp["isolated"] & comp["G"].between(19.0, 20.3) & (comp["VarFlag"] != "VARIABLE")]
    ratios = []
    for src, g in faint.groupby("source"):
        if len(g) < 30:
            continue
        g = g.merge(df[["image", "flux_err_white_ujy", "target_sky_rms", "target_npix", "zp_g"]], on="image",
                    suffixes=("", "_t"))
        pred = np.sqrt(g["target_npix"] * g["target_sky_rms"] ** 2) * 10 ** ((23.9 - g["zp_g"]) / 2.5)
        resid = g["flux_ujy"] - np.median(g["flux_ujy"])
        ratios.append(1.4826 * np.median(np.abs(resid)) / np.median(pred))
    factor = float(np.median(ratios)) if ratios else 1.0
    df["noise_factor"] = factor
    df["flux_err_ujy"] = df["flux_err_white_ujy"] * factor
    times = Time(df["tmid_mjd_utc"].to_numpy(), format="mjd", scale="utc")
    df["bjd_tdb"] = bjd_tdb(times, target_coord(2011.0), PALOMAR)
    df = df.sort_values("bjd_tdb").reset_index(drop=True)
    df.to_csv(OUT_PATH, index=False)
    comp.to_csv(C.DERIVED_DIR / "ptf_comparison_stars.csv", index=False)
    return df


if __name__ == "__main__":
    out = run(force=True)
    print(out[["utc_open", "fwhm_px", "aperture_px", "n_zp", "zp_g", "flux_ujy", "flux_err_ujy"]].to_string())
    print("noise factor", out["noise_factor"].iloc[0])
