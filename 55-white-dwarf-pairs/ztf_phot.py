"""Forced PSF photometry of ZTF J1539+5027 on ZTF difference images.

For every ZTF exposure that covers the star we have a 61x61-pixel cutout of the
difference image (the exposure minus a deep reference image, PSF-matched by IPAC's
pipeline) and the PSF of that difference image. The star's flux *relative to the
reference* is measured by fitting the PSF at the star's known position, whether or
not the star is detected: the faint phases, including the eclipse, count as much
as the bright ones. This is what IPAC's forced-photometry service does; we do it
ourselves on public files.

Per exposure:

* time: mid-point of the shutter-open interval (header SHUTOPEN/SHUTCLSD, UTC),
  converted to Barycentric Dynamical Time at the Solar System barycentre (BJD_TDB)
  for the star's direction and Palomar's position. (OBSJD in the headers is ~0.7 s
  *before* the shutter opens, and the header "HJD" is not heliocentric.)
* flux: sum(P D) / sum(P^2), P the PSF interpolated to the star's sub-pixel position,
  D the difference image minus its local background (median of a 10-25 px annulus);
  calibrated to microjansky with the exposure's zero point MAGZP (AB mag).
* error: the pixel noise (robust RMS of the annulus) propagated through the PSF
  weights, then scaled up by the factor by which fluxes measured the same way at
  ~20 empty positions in the cutout scatter more than that (correlated noise from
  the image resampling and PSF matching).
"""

from __future__ import annotations

import json
import time
import warnings

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning
from scipy.interpolate import RectBivariateSpline

import config as C
import fetch_data as F

PALOMAR = EarthLocation.from_geodetic(lon=-116.8598 * u.deg, lat=33.3573 * u.deg, height=1668 * u.m)
OUT_PATH = C.DERIVED_DIR / "ztf_forced_photometry.csv"
ANNULUS = (10.0, 25.0)  # px, background and noise
N_EMPTY = 24  # empty positions per cutout for the correlated-noise factor
EMPTY_MIN_SEP = 9.0  # px from the star


def target_coord(epoch_jyear: float) -> SkyCoord:
    """The star's ICRS position at ``epoch_jyear`` (Gaia DR3 J2016.0 + proper motion)."""
    dt = epoch_jyear - 2016.0
    dec = C.TARGET_DEC_2016 + C.TARGET_PMDEC * dt / 3.6e6
    ra = C.TARGET_RA_2016 + C.TARGET_PMRA * dt / 3.6e6 / np.cos(np.radians(C.TARGET_DEC_2016))
    return SkyCoord(ra * u.deg, dec * u.deg)


def bjd_tdb(times_utc: Time, coord: SkyCoord, location: EarthLocation) -> np.ndarray:
    """Barycentric Julian Date (TDB) for UTC times observed from ``location``."""
    t = Time(times_utc, location=location)
    ltt = t.light_travel_time(coord, kind="barycentric")
    return (t.tdb + ltt).jd


def psf_interpolator(psf: np.ndarray) -> RectBivariateSpline:
    """Bicubic interpolator of a PSF stamp, coordinates relative to its centroid."""
    ny, nx = psf.shape
    y = np.arange(ny) - (ny - 1) / 2
    x = np.arange(nx) - (nx - 1) / 2
    return RectBivariateSpline(y, x, psf, kx=3, ky=3)


def psf_on_grid(interp, half: float, x0: float, y0: float, shape) -> np.ndarray:
    """PSF centred at (x0, y0) evaluated on the pixel grid of ``shape`` (zero outside)."""
    ny, nx = shape
    yy = np.arange(ny) - y0
    xx = np.arange(nx) - x0
    model = interp(yy, xx, grid=True)
    model[np.abs(yy) > half, :] = 0.0
    model[:, np.abs(xx) > half] = 0.0
    return np.clip(model, 0.0, None)


def forced_flux(diff: np.ndarray, psf_model: np.ndarray) -> float:
    w = psf_model > 0
    good = w & np.isfinite(diff)
    p = psf_model[good]
    return float(np.sum(p * diff[good]) / np.sum(p * p)), float(np.sum(p * p)), int(np.sum(w & ~np.isfinite(diff)))


def measure_one(diff_path, psf_path, rng: np.random.Generator) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        with fits.open(diff_path) as hdul:
            hdu = next(h for h in hdul if h.data is not None)
            data = hdu.data.astype(float)
            hdr = hdu.header
            wcs = WCS(hdr)
    psf = fits.getdata(psf_path).astype(float)
    psf /= psf.sum()
    t0 = Time(hdr["SHUTOPEN"], format="isot", scale="utc")
    t1 = Time(hdr["SHUTCLSD"], format="isot", scale="utc")
    tmid = t0 + (t1 - t0) / 2
    coord = target_coord(tmid.jyear)
    x0, y0 = wcs.world_to_pixel(coord)
    ny, nx = data.shape
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(xx - x0, yy - y0)
    ann = (r >= ANNULUS[0]) & (r <= ANNULUS[1]) & np.isfinite(data)
    bkg = float(np.median(data[ann])) if ann.sum() > 50 else np.nan
    mad = float(1.4826 * np.median(np.abs(data[ann] - bkg))) if ann.sum() > 50 else np.nan
    interp = psf_interpolator(psf)
    half = (psf.shape[0] - 1) / 2
    model = psf_on_grid(interp, half, x0, y0, data.shape)
    sub = data - bkg
    flux, sum_p2, n_bad = forced_flux(sub, model)
    err_white = mad / np.sqrt(sum_p2)
    # the same estimator at empty positions measures the real (correlated) noise
    empties = []
    tries = 0
    while len(empties) < N_EMPTY and tries < 400:
        tries += 1
        xe, ye = rng.uniform(half, nx - 1 - half), rng.uniform(half, ny - 1 - half)
        if np.hypot(xe - x0, ye - y0) < EMPTY_MIN_SEP:
            continue
        if any(np.hypot(xe - a, ye - b) < 6.0 for a, b in empties):
            continue
        empties.append((xe, ye))
    empty_flux = []
    for xe, ye in empties:
        m = psf_on_grid(interp, half, xe, ye, data.shape)
        fe, _, nbe = forced_flux(sub, m)
        if nbe == 0:
            empty_flux.append(fe)
    empty_flux = np.asarray(empty_flux)
    if empty_flux.size >= 8:
        med = np.median(empty_flux)
        rms = 1.4826 * np.median(np.abs(empty_flux - med))
    else:
        rms = np.nan
    zp = hdr["MAGZP"]
    scale = 10 ** ((23.9 - zp) / 2.5)  # DN -> microjansky
    return {
        "shutopen_utc": hdr["SHUTOPEN"], "exptime_s": float((t1 - t0).sec),
        "tmid_mjd_utc": tmid.mjd,
        "x": float(x0), "y": float(y0), "nx": nx, "ny": ny,
        "flux_dn": flux, "err_white_dn": err_white, "empty_rms_dn": rms,
        "bkg_dn": bkg, "pix_rms_dn": mad, "n_bad_in_psf": n_bad,
        "magzp": zp, "magzprms": hdr.get("MAGZPRMS", np.nan), "dn_to_ujy": scale,
        "seeing_px": hdr.get("SEEING", np.nan), "maglim": hdr.get("MAGLIM", np.nan),
        "airmass": hdr.get("AIRMASS", np.nan), "infobits": int(hdr.get("INFOBITS", -1)),
        "status": int(hdr.get("STATUS", -1)), "programid": int(hdr.get("PROGRMID", -1)),
        "qcomment": str(hdr.get("QCOMMENT", "")).strip(),
    }


def run(force: bool = False) -> pd.DataFrame:
    if OUT_PATH.exists() and not force:
        return pd.read_csv(OUT_PATH)
    exposures = pd.read_csv(F.ZTF_EXPOSURES)
    rng = np.random.default_rng(1539)
    rows = []
    t_start = time.time()
    for i, row in enumerate(exposures.to_dict("records")):
        diff_path, psf_path = F.ztf_cutout_paths(row)
        if not (diff_path.exists() and psf_path.exists()):
            continue
        try:
            m = measure_one(diff_path, psf_path, rng)
        except Exception as exc:  # a damaged cutout must not stop the run
            print(f"  skipped {diff_path.name}: {exc}")
            continue
        m.update({"filefracday": row["filefracday"], "field": row["field"], "ccdid": row["ccdid"],
                  "qid": row["qid"], "filtercode": row["filtercode"], "ipac_gid": row["ipac_gid"],
                  "infobits_meta": row["infobits"]})
        rows.append(m)
        if (i + 1) % 500 == 0:
            print(f"  ZTF photometry: {i + 1}/{len(exposures)} ({time.time() - t_start:.0f} s)", flush=True)
    df = pd.DataFrame(rows)
    times = Time(df["tmid_mjd_utc"].to_numpy(), format="mjd", scale="utc")
    df["bjd_tdb"] = bjd_tdb(times, target_coord(2022.0), PALOMAR)
    # correlated-noise factor: empty-position scatter / white-noise prediction, per band
    df["noise_factor"] = df["empty_rms_dn"] / df["err_white_dn"]
    factors = {}
    for band, g in df.groupby("filtercode"):
        factors[band] = float(np.nanmedian(g["noise_factor"]))
    df["flux_err_dn"] = df["err_white_dn"] * df["filtercode"].map(factors)
    df["flux_ujy"] = df["flux_dn"] * df["dn_to_ujy"]
    df["flux_err_ujy"] = df["flux_err_dn"] * df["dn_to_ujy"]
    df = df.sort_values("bjd_tdb").reset_index(drop=True)
    df.to_csv(OUT_PATH, index=False)
    (C.DERIVED_DIR / "ztf_noise_factors.json").write_text(json.dumps(factors, indent=1))
    return df


if __name__ == "__main__":
    out = run(force=True)
    print(out.groupby("filtercode")[["flux_ujy", "flux_err_ujy"]].describe().T)
