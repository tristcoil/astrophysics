"""The real sky behind the black hole, from Gaia DR3 and Hipparcos.

Two layers, so that no star is counted twice:

* **Point stars**: the 1.25 million Gaia DR3 sources with G < 11, plus
  Hipparcos for the naked-eye stars brighter than V = 5, whose Gaia photometry
  is saturated or missing (their Gaia counterparts, found through the DR3
  Hipparcos-2 cross-match, are dropped). Each star is lensed individually.
* **The glow**: every fainter Gaia source (G >= 11, 1.8 billion of them)
  summed per 3.4-arcminute HEALPix pixel on the archive's server. Converted to
  surface brightness and colour, resampled into an equirectangular texture in
  Galactic coordinates, with a mipmap pyramid so that a screen pixel that sees
  a large patch of sky (near the shadow's edge, where rays wind round the
  black hole) gets that patch's average, not one random sample.

Brightness unit: V-band flux of a V = 0 star = 1. Surface brightness: that per
steradian. This is the sky as seen from the Sun (see the README for why we use
it at Sgr A*).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import healpy as hp
import numpy as np
from astropy.coordinates import Galactic, ICRS, SkyCoord
from astropy.io import fits
import astropy.units as u

import colour
import config as C

STARS_PATH = C.GAIA_DIR / "gaia_dr3_g_lt_11.fits"
FLUX_PATH = C.GAIA_DIR / f"gaia_dr3_flux_map_L{C.FLUX_MAP_LEVEL}.npz"
HIP_PATH = C.GAIA_DIR / "hipparcos_main.fits"
HIPXM_PATH = C.GAIA_DIR / "gaiadr3_hipparcos2_best_neighbour.fits"
POINTS_CACHE = C.CACHE_DIR / "sky_points.npz"
TEXTURE_CACHE = C.CACHE_DIR / f"sky_texture_{C.TEXTURE_WIDTH}.npz"

#: A few naked-eye stars to label in figures (Hipparcos numbers).
NAMED_STARS = {32349: "Sirius", 30438: "Canopus", 69673: "Arcturus", 91262: "Vega",
               24608: "Capella", 24436: "Rigel", 37279: "Procyon", 27989: "Betelgeuse",
               80763: "Antares", 21421: "Aldebaran", 65474: "Spica", 49669: "Regulus",
               102098: "Deneb", 97649: "Altair", 113368: "Fomalhaut", 71683: "Alpha Centauri",
               11767: "Polaris", 90185: "Kaus Australis", 92855: "Nunki", 86032: "Rasalhague",
               17702: "Pleiades (Alcyone)"}


def _rot_icrs_to_gal() -> np.ndarray:
    """3x3 rotation taking ICRS unit vectors to Galactic unit vectors (from astropy)."""
    pts = SkyCoord(ra=[0.0, 90.0, 0.0] * u.deg, dec=[0.0, 0.0, 90.0] * u.deg, frame=ICRS)
    gal = pts.transform_to(Galactic())
    cols = gal.cartesian.xyz.value          # columns = images of x, y, z
    return cols


R_ICRS_TO_GAL = _rot_icrs_to_gal()


def radec_to_galvec(ra_deg, dec_deg) -> np.ndarray:
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    v = np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=-1)
    return v @ R_ICRS_TO_GAL.T


def galvec_to_lb(v):
    v = np.asarray(v)
    lon = np.arctan2(v[..., 1], v[..., 0]) % (2 * np.pi)
    lat = np.arcsin(np.clip(v[..., 2], -1.0, 1.0))
    return lon, lat


# ---------------------------------------------------------------------------
# Point stars
# ---------------------------------------------------------------------------

@dataclass
class Points:
    vec: np.ndarray        # (N, 3) Galactic unit vectors (float64)
    vmag: np.ndarray       # V magnitude
    logT: np.ndarray       # log10 colour temperature (K)
    source: np.ndarray     # 0 = Gaia, 1 = Hipparcos
    hip: np.ndarray        # Hipparcos number or 0


def build_points() -> Points:
    with fits.open(STARS_PATH) as h:
        d = h[1].data
        sid = np.asarray(d["source_id"], np.int64)
        ra, dec = np.asarray(d["ra"], float), np.asarray(d["dec"], float)
        gmag = np.asarray(d["phot_g_mean_mag"], float)
        bp, rp = np.asarray(d["phot_bp_mean_mag"], float), np.asarray(d["phot_rp_mean_mag"], float)
    with fits.open(HIP_PATH) as h:
        d = h[1].data
        hip = np.asarray(d["hip"], np.int64)
        hra, hde = np.asarray(d["ra"], float), np.asarray(d["de"], float)
        vmag = np.asarray(d["vmag"], float)
        bv = np.asarray(d["b_v"], float)
    with fits.open(HIPXM_PATH) as h:
        xm_sid = np.asarray(h[1].data["source_id"], np.int64)
        xm_hip = np.asarray(h[1].data["original_ext_source_id"], np.int64)

    bright = np.isfinite(vmag) & (vmag < C.HIP_V_LIMIT) & np.isfinite(hra) & np.isfinite(hde)
    bright_hip = set(hip[bright].tolist())
    drop_sid = set(xm_sid[np.isin(xm_hip, list(bright_hip))].tolist())
    keep = ~np.isin(sid, list(drop_sid))

    bprp = bp - rp
    bprp = np.where(np.isfinite(bprp), bprp, 0.82)          # (no colour: assume Sun-like; ~0.1% of stars)
    v_gaia = colour.v_from_gaia(gmag, bprp)
    lt_gaia = colour.logT_from_bprp(bprp)

    bv_h = np.where(np.isfinite(bv), bv, 0.65)
    vec = np.concatenate([radec_to_galvec(ra[keep], dec[keep]), radec_to_galvec(hra[bright], hde[bright])])
    pts = Points(vec=vec,
                 vmag=np.concatenate([v_gaia[keep], vmag[bright]]),
                 logT=np.concatenate([lt_gaia[keep], colour.logT_from_bv(bv_h[bright])]),
                 source=np.concatenate([np.zeros(keep.sum(), np.int8), np.ones(bright.sum(), np.int8)]),
                 hip=np.concatenate([np.zeros(keep.sum(), np.int64), hip[bright]]))
    np.savez(POINTS_CACHE, vec=pts.vec.astype(np.float64), vmag=pts.vmag.astype(np.float32),
             logT=pts.logT.astype(np.float32), source=pts.source, hip=pts.hip,
             n_gaia_dropped=len(drop_sid), n_hip=int(bright.sum()))
    return pts


@lru_cache(maxsize=1)
def points() -> Points:
    if not POINTS_CACHE.exists():
        return build_points()
    z = np.load(POINTS_CACHE)
    return Points(vec=z["vec"], vmag=z["vmag"].astype(float), logT=z["logT"].astype(float),
                  source=z["source"], hip=z["hip"])


# ---------------------------------------------------------------------------
# The glow of the fainter stars: HEALPix -> equirectangular texture + mipmaps
# ---------------------------------------------------------------------------

def flux_map_physical():
    """Per HEALPix pixel: V-band surface brightness (V=0 flux per sr) and log10 colour temperature."""
    z = np.load(FLUX_PATH)
    sg, sbp, srp = (z[k].astype(np.float64) for k in ("sum_g", "sum_bp", "sum_rp"))
    nside = 2**int(z["level"])
    omega = hp.nside2pixarea(nside)
    with np.errstate(divide="ignore", invalid="ignore"):
        gmag = -2.5 * np.log10(sg) + C.ZP_G
        bprp = -2.5 * np.log10(sbp / srp) + (C.ZP_BP - C.ZP_RP)
    ok_col = np.isfinite(bprp)
    bprp = np.where(ok_col, bprp, np.nanmedian(bprp[ok_col]))
    v = colour.v_from_gaia(gmag, bprp)
    lum = np.where(sg > 0, 10.0 ** (-0.4 * v), 0.0) / omega
    lt = colour.logT_from_bprp(bprp)
    return nside, lum, lt


def build_texture(width: int = C.TEXTURE_WIDTH):
    nside, lum, lt = flux_map_physical()
    height = width // 2
    lon = (np.arange(width) + 0.5) / width * 2 * np.pi
    lat = np.pi / 2 - (np.arange(height) + 0.5) / height * np.pi     # row 0 = north
    out_lum = np.empty((height, width), np.float32)
    out_lt = np.empty((height, width), np.float32)
    rot = R_ICRS_TO_GAL.T            # Galactic -> ICRS
    lum_lt = lum * lt
    step = 256
    for i0 in range(0, height, step):
        la = lat[i0:i0 + step]
        LON, LAT = np.meshgrid(lon, la)
        vg = np.stack([np.cos(LAT) * np.cos(LON), np.cos(LAT) * np.sin(LON), np.sin(LAT)], axis=-1)
        vi = vg @ rot.T
        theta = np.arccos(np.clip(vi[..., 2], -1, 1))
        phi = np.arctan2(vi[..., 1], vi[..., 0]) % (2 * np.pi)
        a = hp.get_interp_val(lum, theta.ravel(), phi.ravel(), nest=True)
        b = hp.get_interp_val(lum_lt, theta.ravel(), phi.ravel(), nest=True)
        out_lum[i0:i0 + step] = a.reshape(LON.shape)
        out_lt[i0:i0 + step] = (b / np.maximum(a, 1e-300)).reshape(LON.shape)
    np.savez(TEXTURE_CACHE, lum=out_lum, lt=out_lt)
    return out_lum, out_lt


class SkyTexture:
    """Mipmapped equirectangular sky in Galactic coordinates (row 0 = b = +90)."""

    def __init__(self, lum: np.ndarray, lt: np.ndarray):
        self.levels = [(lum.astype(np.float32), (lum * lt).astype(np.float32))]
        while self.levels[-1][0].shape[1] > 64:
            a, b = self.levels[-1]
            h, w = a.shape
            a2 = a.reshape(h // 2, 2, w // 2, 2).mean(axis=(1, 3))
            b2 = b.reshape(h // 2, 2, w // 2, 2).mean(axis=(1, 3))
            self.levels.append((a2, b2))
        self.width = lum.shape[1]
        self.texel0 = 2 * np.pi / self.width

    @classmethod
    def load(cls):
        if not TEXTURE_CACHE.exists():
            lum, lt = build_texture()
        else:
            z = np.load(TEXTURE_CACHE)
            lum, lt = z["lum"], z["lt"]
        return cls(lum, lt)

    def _bilinear(self, level: int, lon, lat):
        a, b = self.levels[level]
        h, w = a.shape
        x = lon / (2 * np.pi) * w - 0.5
        y = (np.pi / 2 - lat) / np.pi * h - 0.5
        x0 = np.floor(x).astype(np.int64)
        y0 = np.floor(y).astype(np.int64)
        fx, fy = x - x0, y - y0
        x0w, x1w = x0 % w, (x0 + 1) % w
        y0c, y1c = np.clip(y0, 0, h - 1), np.clip(y0 + 1, 0, h - 1)
        out = []
        for arr in (a, b):
            v = ((1 - fx) * (1 - fy) * arr[y0c, x0w] + fx * (1 - fy) * arr[y0c, x1w]
                 + (1 - fx) * fy * arr[y1c, x0w] + fx * fy * arr[y1c, x1w])
            out.append(v)
        return out

    def sample(self, vec, footprint):
        """Luminance per sr and log10 T at Galactic unit vectors ``vec`` (..., 3).

        ``footprint`` is the angular size (radians) of the patch of sky each sample
        stands for; the mip level is chosen to match it (trilinear blend).
        """
        lon, lat = galvec_to_lb(vec)
        lev = np.log2(np.maximum(footprint, 1e-12) / self.texel0)
        lev = np.clip(lev, 0, len(self.levels) - 1)
        l0 = np.floor(lev).astype(int)
        frac = lev - l0
        lum = np.zeros(lon.shape, np.float64)
        llt = np.zeros(lon.shape, np.float64)
        for L in np.unique(l0):
            m = l0 == L
            a0, b0 = self._bilinear(L, lon[m], lat[m])
            if L + 1 < len(self.levels):
                a1, b1 = self._bilinear(L + 1, lon[m], lat[m])
                f = frac[m]
                a0, b0 = (1 - f) * a0 + f * a1, (1 - f) * b0 + f * b1
            lum[m], llt[m] = a0, b0
        lt = np.where(lum > 0, llt / np.maximum(lum, 1e-300), 3.7)
        return lum, lt
