"""TESS photometry: our own aperture photometry on full-frame-image cutouts,
and SPOC 2-minute light curves.

Cutouts (TESScut) hold SPOC-calibrated, background-subtracted pixels. We sum
a circular aperture around the target, subtract any residual background
(sigma-clipped median of an annulus, per cadence) and keep cadences with
QUALITY == 0. In crowded fields the aperture also collects neighbours' light,
which dilutes eclipses; ``dilution`` estimates the target's share of the
aperture flux from Gaia DR3 with a Gaussian approximation to TESS's
point-spread function.
"""

from __future__ import annotations

import glob
from dataclasses import dataclass

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

import config as C

BTJD0 = 2457000.0


def _slug(name):
    return name.lower().replace(" ", "_").replace("-", "_")


@dataclass
class Cutout:
    sector: int
    time: np.ndarray          # BJD_TDB - 2457000
    flux: np.ndarray          # (n, ny, nx) e-/s
    flux_err: np.ndarray
    quality: np.ndarray
    wcs: WCS
    exptime_days: float


def read_cutout(path) -> Cutout:
    with fits.open(path) as h:
        d = h[1].data
        sector = int(h[0].header["SECTOR"])
        wcs = WCS(h[2].header)
        dt = float(h[1].header.get("TIMEDEL", np.nan))
        return Cutout(sector, np.array(d["TIME"], float), np.array(d["FLUX"], float),
                      np.array(d["FLUX_ERR"], float), np.array(d["QUALITY"], int), wcs, dt)


def cutouts(name: str) -> list[Cutout]:
    paths = sorted(glob.glob(str(C.RAW / "tess" / _slug(name) / "tess-s*_astrocut.fits")))
    return [read_cutout(p) for p in paths]


def aperture_photometry(cut: Cutout, ra: float, dec: float, r_ap: float = 2.0,
                        annulus=(4.5, 7.0)) -> tuple[pd.DataFrame, dict]:
    x0, y0 = cut.wcs.world_to_pixel_values(ra, dec)
    ny, nx = cut.flux.shape[1:]
    yy, xx = np.mgrid[0:ny, 0:nx]
    rr = np.hypot(xx - x0, yy - y0)
    ap = rr <= r_ap
    ann = (rr >= annulus[0]) & (rr <= annulus[1])
    good = (cut.quality == 0) & np.isfinite(cut.time)
    good &= np.all(np.isfinite(cut.flux[:, ap]), axis=1)
    f = cut.flux[good]
    e = cut.flux_err[good]
    # residual background: clipped median of the annulus pixels
    a = f[:, ann]
    med = np.nanmedian(a, axis=1)
    mad = 1.4826 * np.nanmedian(np.abs(a - med[:, None]), axis=1)
    clip = np.abs(a - med[:, None]) < 3 * np.maximum(mad[:, None], 1e-9)
    bkg = np.array([np.median(row[c]) if c.any() else m for row, c, m in zip(a, clip, med)])
    flux = f[:, ap].sum(axis=1) - ap.sum() * bkg
    err = np.sqrt((e[:, ap] ** 2).sum(axis=1))
    df = pd.DataFrame({"btjd": cut.time[good], "flux": flux, "err": err, "bkg": bkg})
    df = df[np.isfinite(df["flux"]) & (df["flux"] > 0)]
    meta = dict(sector=cut.sector, x0=float(x0), y0=float(y0), n_ap=int(ap.sum()), r_ap=r_ap,
                exptime_days=cut.exptime_days, n=int(len(df)), aperture=ap)
    return df.reset_index(drop=True), meta


def tess_mag_from_gaia(g, bp_rp):
    """Stassun et al. (2019, AJ 158, 138) TESS magnitude from Gaia G and BP-RP."""
    c = np.nan_to_num(np.asarray(bp_rp, float), nan=1.0)
    return g - 0.00522555 * c**3 + 0.0891337 * c**2 - 0.633923 * c + 0.0324473


def dilution(name: str, cut: Cutout, meta: dict, sigma_px: float = 0.75,
             target_radius_arcsec: float = 2.0) -> dict:
    """Share of the aperture flux that comes from the target, from Gaia DR3
    neighbours (G < 19 within 4') and a Gaussian PSF of width sigma_px."""
    s = C.SYSTEMS[name]
    g = pd.read_csv(C.RAW / "tess" / _slug(name) / "gaia_dr3_field.csv")
    g.columns = [c.lower() for c in g.columns]
    tmag = tess_mag_from_gaia(g["phot_g_mean_mag"], g["bp_rp"])
    flux = 10 ** (-0.4 * (tmag - 10.0))
    x, y = cut.wcs.world_to_pixel_values(g["ra"].to_numpy(), g["dec"].to_numpy())
    ap = meta["aperture"]
    ny, nx = ap.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    # fraction of each star's light in the aperture: integrate the Gaussian
    # over each pixel (erf differences), sum over aperture pixels
    from scipy.special import erf
    frac = np.zeros(len(g))
    s2 = np.sqrt(2) * sigma_px
    for k in range(len(g)):
        if not (-3 < x[k] < nx + 2 and -3 < y[k] < ny + 2):
            continue
        fx = 0.5 * (erf((xx + 0.5 - x[k]) / s2) - erf((xx - 0.5 - x[k]) / s2))
        fy = 0.5 * (erf((yy + 0.5 - y[k]) / s2) - erf((yy - 0.5 - y[k]) / s2))
        frac[k] = np.sum((fx * fy)[ap])
    sep = 3600 * np.hypot((g["ra"] - s["ra"]) * np.cos(np.radians(s["dec"])), g["dec"] - s["dec"])
    tgt = sep < target_radius_arcsec
    contrib = flux * frac
    share = float(contrib[tgt].sum() / contrib.sum())
    top = np.argsort(-contrib.to_numpy() if hasattr(contrib, "to_numpy") else -contrib)[:5]
    return dict(target_share=share, n_target_sources=int(tgt.sum()),
                brightest_contaminants=[dict(sep_arcsec=float(sep[k]), tmag=float(tmag[k]),
                                             share=float(contrib[k] / contrib.sum())) for k in top])


def spoc_lightcurves(name: str) -> list[tuple[int, pd.DataFrame, dict]]:
    out = []
    for p in sorted(glob.glob(str(C.RAW / "tess" / _slug(name) / "tess*_lc.fits"))):
        with fits.open(p) as h:
            d = h[1].data
            sec = int(h[0].header["SECTOR"])
            hdr = dict(CROWDSAP=h[1].header.get("CROWDSAP"), FLFRCSAP=h[1].header.get("FLFRCSAP"))
            good = (d["QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"]) & np.isfinite(d["TIME"])
            df = pd.DataFrame({"btjd": np.array(d["TIME"][good], float),
                               "flux": np.array(d["PDCSAP_FLUX"][good], float),
                               "err": np.array(d["PDCSAP_FLUX_ERR"][good], float),
                               "sap": np.array(d["SAP_FLUX"][good], float)})
            out.append((sec, df, hdr))
    return out
