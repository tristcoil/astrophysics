"""The survey footprint, and how deep the survey is across it.

A cluster count is only as good as its denominator: how much sky was searched,
and how sensitive the search was at each point. This module rebuilds that from
public pieces, on HEALPix pixels of 0.11 deg (order 9, nested):

* **Exposure**: the eRASS1 0.2-2.3 keV vignetted exposure map. eROSITA scans
  great circles through the ecliptic poles, so a typical patch of sky got
  ~150 s while the south ecliptic pole got tens of kiloseconds.
* **Where the clusters could be found**: the eROSITA-DE half of the sky, inside
  the Legacy Surveys DR10 g, r and z coverage (clusters are confirmed and
  given redshifts from optical images), south of Dec 32.375 (DR10-South).
* **What Bulbul et al. masked**: the Magellanic Clouds, the Virgo cluster,
  bright supernova remnants and a few smaller areas. No mask file is public,
  but their Fig. 2 draws the cleaned footprint. We read it off the image.
  The map there is drawn with round markers that spill ~0.5 deg past every
  edge, so we shrink the read-off area by the marker radius, calibrated so
  that the whole map has the published 13,116 deg^2. The cosmology footprint
  we then build (12,7xx deg^2) is an independent check against the published 12,791.
* **Galactic absorption** (HI4PI N_H) and the **optical depth limit**
  ZVLIM (the highest redshift at which the Legacy Surveys images are deep
  enough to count a cluster's galaxies; from the nearest catalogued source).
"""

from __future__ import annotations

import numpy as np
import healpy as hp
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from scipy import ndimage
from scipy.spatial import cKDTree

import config as C

ORDER = 9
NSIDE = 2**ORDER
NPIX = hp.nside2npix(NSIDE)
PIX_AREA = hp.nside2pixarea(NSIDE, degrees=True)
CACHE = C.CACHE_DIR / "footprint.npz"


# ---------------------------------------------------------------------------
# Maps from public files
# ---------------------------------------------------------------------------


def exposure_map() -> np.ndarray:
    """eRASS1 exposure (s) on nested HEALPix order 9, from the HiPS Allsky file.

    The Allsky image holds the 768 order-3 tiles, 27 per row, each 64x64
    pixels. With the image flipped top-to-bottom, tile pixel (row i, column j)
    is face pixel (x0 + i, y0 + j) of HEALPix order 9. That orientation was
    found by trying all eight and matching the catalogue's own per-cluster
    exposures (EXP): it reproduces them to a median 0.5% (test_physics.py).
    """
    img = fits.getdata(C.RAW_DIR / "erass1_exposure_024_Norder3_Allsky.fits").astype(float)[::-1]
    w, ncol, nside_t = 64, 27, 8
    out = np.full(NPIX, np.nan)
    ii, jj = np.meshgrid(np.arange(w), np.arange(w), indexing="ij")
    for t in range(12 * nside_t**2):
        r, c = divmod(t, ncol)
        xt, yt, face = hp.pix2xyf(nside_t, t, nest=True)
        pix = hp.xyf2pix(NSIDE, xt * w + ii, yt * w + jj, face, nest=True)
        out[pix] = img[r * w:(r + 1) * w, c * w:(c + 1) * w]
    out[~(out > 0)] = np.nan
    return out


def moc_mask(path) -> np.ndarray:
    """A MOC (NUNIQ cells) as a boolean nested map at order 9."""
    uniq = fits.getdata(path, 1)["UNIQ"].astype(np.int64)
    order = (np.floor(np.log2(uniq / 4) / 2)).astype(np.int64)
    ipix = uniq - 4 * 4**order
    m = np.zeros(NPIX, bool)
    for o in np.unique(order):
        sel = ipix[order == o]
        if o <= ORDER:
            f = 4 ** (ORDER - o)
            for p in sel:
                m[p * f:(p + 1) * f] = True
        else:
            m[sel // 4 ** (o - ORDER)] = True
    return m


def pixel_coords():
    ra, dec = hp.pix2ang(NSIDE, np.arange(NPIX), nest=True, lonlat=True)
    return ra, dec


def nh_map(ra, dec) -> np.ndarray:
    """HI4PI column density (cm^-2), sampled at the given positions."""
    with fits.open(C.RAW_DIR / "hi4pi_nhi_car_0p25deg.fits") as h:
        data, wcs = h[0].data.astype(float), WCS(h[0].header)
    x, y = wcs.world_to_pixel_values(ra, dec)
    return ndimage.map_coordinates(data, [y, x], order=1, mode="wrap")


# ---------------------------------------------------------------------------
# Reading Bulbul et al. (2024) Fig. 2
# ---------------------------------------------------------------------------


#: The Mollweide ellipse of Fig. 2 (left), fitted to its border (rms 2.6 px);
#: the graticule then falls within 1-2 px (0.1 deg) of where it should.
FIG2A_ELLIPSE = (2561.79, 1614.2, 1937.7, 968.69)


def _fig2a_classes() -> np.ndarray:
    """Which pixels of the footprint figure are covered sky (boolean image)."""
    import matplotlib.image as mpimg

    im = mpimg.imread(C.RAW_DIR / "bulbul2024_fig2a_footprint.png")[..., :3]
    H, W, _ = im.shape
    mx, mn = im.max(axis=2), im.min(axis=2)
    del im
    white = mn > 0.93
    grey = (mn > 0.5) & (mx < 0.85) & ((mx - mn) < 0.04)  # graticule
    dark = mx < 0.3  # border, labels
    del mx, mn
    known = ~(grey | dark)
    covered = known & ~white
    del white, grey, dark
    # Graticule lines and labels hide what is under them: give each hidden
    # pixel the majority class of the visible pixels within 10 px.
    frac_cov = ndimage.uniform_filter(covered.astype(np.float32), 21)
    frac_known = ndimage.uniform_filter(known.astype(np.float32), 21)
    filled = np.where(known, covered, frac_cov > 0.5 * frac_known)
    del frac_cov, frac_known, covered, known
    cx, cy, a, b = FIG2A_ELLIPSE
    yy, xx = np.ogrid[0:H, 0:W]
    inside = ((xx - cx) / a) ** 2 + ((yy - cy) / b) ** 2 < 1
    filled &= inside
    # Rendering noise: drop holes and islands smaller than 30 px^2 (0.3 deg^2).
    for target in (False, True):
        lab, n = ndimage.label((filled == target) & inside)
        sizes = np.bincount(lab.ravel())
        tiny_labels = np.where(sizes < 30)[0]
        tiny_labels = tiny_labels[tiny_labels > 0]
        filled[np.isin(lab, tiny_labels)] = not target
        del lab
    return filled


def _mollweide_pixels(ra, dec, ellipse):
    """Image pixel of each sky position (RA increases to the left, RA 0 central)."""
    cx, cy, a, b = ellipse
    phi = np.radians(dec)
    th = phi.copy()
    for _ in range(60):
        th -= (2 * th + np.sin(2 * th) - np.pi * np.sin(phi)) / (2 + 2 * np.cos(2 * th) + 1e-12)
    lam = -np.radians(np.where(ra > 180, ra - 360, ra))
    return cx + lam * np.cos(th) / np.pi * a, cy - np.sin(th) * b


def fig2a_masks(ra, dec, radii_px) -> dict:
    """Bulbul et al.'s cleaned footprint, shrunk by each trial marker radius."""
    filled = _fig2a_classes()
    H, W = filled.shape
    px, py = _mollweide_pixels(ra, dec, FIG2A_ELLIPSE)
    iy = np.clip(np.round(py).astype(np.int32), 0, H - 1)
    ix = np.clip(np.round(px).astype(np.int32), 0, W - 1)
    del px, py
    out = {}
    for radius in radii_px:
        if radius > 0:
            r = np.arange(-radius, radius + 1)
            disk = np.add.outer(r**2, r**2) <= radius**2
            eroded = ndimage.binary_erosion(filled, structure=disk)
        else:
            eroded = filled
        out[radius] = eroded[iy, ix]
        del eroded
    return out


def digitize_fig2b():
    """The 'eRASS1, clean' and 'total' curves of Bulbul et al. Fig. 2 (right).

    Both are drawn in the same blue; 'clean' is the lower curve. Axes from the
    major ticks: 235.5 px per decade of exposure (0.1 ks at x = 209) and
    110.6 px per decade of area (10^4 deg^2 at y = 71). Returns exposure (s)
    and cumulative area (deg^2) for each curve.
    """
    import matplotlib.image as mpimg

    rgb = (mpimg.imread(C.RAW_DIR / "bulbul2024_fig2b_area_exposure.png")[..., :3] * 255).round()
    blue = (np.abs(rgb[..., 0] - 76) < 30) & (np.abs(rgb[..., 1] - 76) < 30) & (rgb[..., 2] > 225)
    blue[:175, 470:] = False  # the legend
    clean, total = [], []
    for x in range(160, 820):
        rows = np.where(blue[20:733, x])[0] + 20
        if rows.size == 0:
            continue
        runs = np.split(rows, np.where(np.diff(rows) > 2)[0] + 1)
        centres = [r.mean() for r in runs]
        t = 10 ** (-1 + (x - 209.0) / 235.5) * 1e3
        if len(centres) >= 2:
            clean.append((t, 10 ** (4 - (max(centres) - 71.0) / 110.6)))
            total.append((t, 10 ** (4 - (min(centres) - 71.0) / 110.6)))
    return np.array(clean), np.array(total)


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def _zvlim_map(ra, dec, sel):
    """Optical depth limit at each selected pixel: ZVLIM of the nearest source.

    ZVLIM (Kluge et al. 2024) depends only on how deep the Legacy Surveys
    images are at a position, so every catalogued source samples it. With
    ~11,700 sources the typical spacing is ~1 deg.
    """
    K = Table.read(C.RAW_DIR / "eRASS1_clusters_optical.fits")
    run = np.asarray(K["RUN"]).astype(str)
    zv = np.asarray(K["ZVLIM_02"], float)
    ok = np.isfinite(zv) & np.char.startswith(run, "legacy_dr10")
    vec_src = np.array(hp.ang2vec(np.asarray(K["RA"], float)[ok], np.asarray(K["DEC"], float)[ok], lonlat=True))
    vec_pix = np.array(hp.ang2vec(ra[sel], dec[sel], lonlat=True))
    _, j = cKDTree(vec_src).query(vec_pix)
    out = np.full(ra.size, np.nan)
    out[sel] = zv[ok][j]
    return out


def build(radius_px: int | None = None, verbose: bool = True) -> dict:
    """Assemble the footprint and per-pixel survey conditions (cached)."""
    if CACHE.exists() and radius_px is None:
        return dict(np.load(CACHE))
    ra, dec = pixel_coords()
    expo = exposure_map()
    de = moc_mask(C.RAW_DIR / "erass1_exposure_024_Moc.fits") & np.isfinite(expo)
    ls = np.ones(NPIX, bool)
    for band in "grz":
        ls &= moc_mask(C.RAW_DIR / f"ls_dr10_{band}_moc8.fits")
    south = dec < C.DEC_MAX_DR10_SOUTH

    # Calibrate the marker radius on the published primary-footprint area.
    radii = range(0, 10) if radius_px is None else [radius_px]
    masks = fig2a_masks(ra, dec, radii)
    trials = {r: float((m & de).sum() * PIX_AREA) for r, m in masks.items()}
    if radius_px is None:
        radius_px = min(trials, key=lambda r: abs(trials[r] - C.AREA_PRIMARY_DEG2))
        if verbose:
            print("  Fig. 2 read-off area by marker radius (px -> deg^2):",
                  {k: round(v) for k, v in trials.items()}, "->", radius_px)
    clean = masks[radius_px]
    del masks
    primary = clean & de
    cosmo = primary & ls & south
    nh = np.full(NPIX, np.nan)
    nh[primary] = nh_map(ra[primary], dec[primary])
    zvlim = _zvlim_map(ra, dec, cosmo)
    out = {
        "nside": NSIDE, "radius_px": radius_px,
        "exposure": expo.astype(np.float32), "nh": nh.astype(np.float32), "zvlim": zvlim.astype(np.float32),
        "mask_de": de, "mask_ls": ls, "mask_clean": clean, "mask_primary": primary, "mask_cosmo": cosmo,
        "area_primary": primary.sum() * PIX_AREA, "area_cosmo": cosmo.sum() * PIX_AREA,
        "area_base": (de & ls & south).sum() * PIX_AREA,
        "trial_radii": np.array(list(trials)), "trial_areas": np.array(list(trials.values())),
    }
    np.savez_compressed(CACHE, **out)
    return out


def conditions_at(fp: dict, ra, dec):
    """Exposure, N_H and ZVLIM of the footprint pixels at given positions."""
    pix = hp.ang2pix(NSIDE, np.asarray(ra, float), np.asarray(dec, float), lonlat=True, nest=True)
    return fp["exposure"][pix], fp["nh"][pix], fp["zvlim"][pix], fp["mask_cosmo"][pix]


def survey_classes(fp: dict, t_edges=None, nh_edges=None, zv_edges=None) -> dict:
    """Group the cosmology footprint into classes of (exposure, N_H, ZVLIM).

    The expected counts need the selection averaged over the sky. Pixels with
    the same exposure, absorption and optical depth are interchangeable, so
    the sky collapses to a few hundred weighted classes.
    """
    m = fp["mask_cosmo"] & np.isfinite(fp["exposure"]) & np.isfinite(fp["nh"]) & np.isfinite(fp["zvlim"])
    t, nh, zv = (np.asarray(fp[k], float)[m] for k in ("exposure", "nh", "zvlim"))
    t_edges = np.geomspace(40, 30000, 46) if t_edges is None else t_edges
    nh_edges = np.geomspace(5e19, 2e21, 11) if nh_edges is None else nh_edges
    zv_edges = np.linspace(0.3, 1.3, 21) if zv_edges is None else zv_edges
    it = np.clip(np.digitize(t, t_edges) - 1, 0, len(t_edges) - 2)
    inh = np.clip(np.digitize(nh, nh_edges) - 1, 0, len(nh_edges) - 2)
    izv = np.clip(np.digitize(zv, zv_edges) - 1, 0, len(zv_edges) - 2)
    key = (it * (len(nh_edges) - 1) + inh) * (len(zv_edges) - 1) + izv
    uniq, inv, counts = np.unique(key, return_inverse=True, return_counts=True)
    w = counts / counts.sum()
    mean = lambda v: np.bincount(inv, weights=v) / counts  # noqa: E731
    return {"t": mean(t), "nh": mean(nh), "zvlim": mean(zv), "weight": w,
            "area_deg2": m.sum() * PIX_AREA, "n_pixels": int(m.sum())}
