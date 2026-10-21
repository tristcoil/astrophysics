"""Measure how filaments move between two registered images, and rewind them.

The recipe, for two images of the same remnant taken years apart:

1. **Remove the stars.** Stars sit still (after registration) while the
   nebula moves; left in, they pull every correlation towards zero shift.
   Each Gaia star is masked out to a radius that grows with its brightness,
   and the hole is filled from the surrounding nebula.
2. **Make the images comparable.** Blur both to the same resolution, then
   replace each pixel by its rank among the nebula's pixels. A photographic
   plate records density, roughly the logarithm of the light, while a CCD
   records the light itself; any such monotonic difference disappears in
   ranks. Then subtract a smooth version, keeping the filaments and dropping
   the broad glow of synchrotron light, whose share differs between bands.
3. **Track patches.** Cut a small square from the first image, slide it over
   the second, and find where it matches best (normalised cross-correlation),
   to a fraction of a pixel. Then slide back from the match and keep the
   vector only if we return to where we started.
4. **Rewind.** An undecelerated, homologous expansion puts every filament at
   p(t) = c + v (t - T): all of them were at the centre c at time T. Fit c and
   T to all the vectors at once.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage, optimize
from skimage.feature import match_template

import astrometry as A
import config


# --------------------------------------------------------------------------
# 1-2. Preparing images
# --------------------------------------------------------------------------

def star_radii(reg, grid: A.Grid, gaia: A.GaiaStars, *, threshold: float = 1.0,
               bins=np.arange(9.0, 21.0, 1.0)) -> tuple[np.ndarray, np.ndarray]:
    """How far each star's light reaches in this image, by magnitude.

    For stars in one-magnitude bins, take the median radial profile (after
    subtracting each star's local background, 25-35" out) and find where it
    drops below ``threshold`` times the image's pixel noise. Deep CCD stacks
    show a G = 15 star out to ~8"; a Schmidt plate's saturated G = 11 star
    spreads over ~20". Returns (bin centres, radius in arcsec); bins without
    enough stars are filled by extrapolating the trend."""
    img = reg.image.astype(float)
    noise = A.robust_sigma(img - A.background(img, int(20 / grid.step)))
    ra, dec = gaia.at_epoch(reg.epoch)
    ok = np.isfinite(ra)
    xi, eta = A.to_tangent(ra[ok], dec[ok], grid.ra0, grid.dec0)
    col, row = grid.to_index(xi, eta)
    g = gaia.g[ok]
    rmax = int(35 / grid.step)
    yy, xx = np.mgrid[-rmax : rmax + 1, -rmax : rmax + 1]
    rad = np.hypot(xx, yy) * grid.step
    edges = np.arange(0, 35.5, max(grid.step, 0.5))
    centres, radii = [], []
    for lo in bins:
        sel = (g >= lo) & (g < lo + 1)
        profiles = []
        for c, r in zip(col[sel], row[sel]):
            ic, ir = int(round(c)), int(round(r))
            if ir - rmax < 0 or ic - rmax < 0 or ir + rmax >= grid.n or ic + rmax >= grid.n:
                continue
            st = img[ir - rmax : ir + rmax + 1, ic - rmax : ic + rmax + 1]
            if not np.all(np.isfinite(st)):
                continue
            st = st - np.median(st[(rad > 25) & (rad < 35)])
            prof = [np.median(st[(rad >= a) & (rad < b)]) for a, b in zip(edges[:-1], edges[1:])]
            profiles.append(prof)
        if len(profiles) < 5:
            continue
        med = np.median(np.array(profiles), axis=0)
        below = np.flatnonzero(med < threshold * noise)
        centres.append(lo + 0.5)
        radii.append(edges[below[0]] if len(below) else edges[-1])
    centres, radii = np.array(centres), np.array(radii)
    # Bins with too few stars (the brightest) come from a log-linear fit of
    # radius against G; then force the radius never to grow with G.
    good = radii > 0
    slope, icept = np.polyfit(centres[good], np.log(np.maximum(radii[good], grid.step)), 1)
    out = np.exp(icept + slope * (bins + 0.5))
    for c, r in zip(centres, radii):
        out[np.isclose(bins + 0.5, c)] = r
    out = np.maximum.accumulate(out[::-1])[::-1]
    return bins + 0.5, out


def star_mask(grid: A.Grid, gaia: A.GaiaStars, epoch: float, fwhm_arcsec: float, *,
              gmax: float = 20.0, radius_model: tuple | None = None, pad: float = 0.0) -> np.ndarray:
    """True where a Gaia star's light falls, at the image's own epoch.

    With ``radius_model`` (from ``star_radii``) each star's radius is read off
    the image's own profile; without it, a rough rule (one FWHM at G = 19,
    growing 35% of a FWHM per magnitude) is used. ``pad`` (arcsec) is added,
    e.g. half the blur applied afterwards. Capped at 40"."""
    ra, dec = gaia.at_epoch(epoch)
    ra = np.where(np.isfinite(ra), ra, gaia.ra)  # stars without proper motions stay at their 2016 place
    dec = np.where(np.isfinite(dec), dec, gaia.dec)
    sel = gaia.g < gmax
    xi, eta = A.to_tangent(ra[sel], dec[sel], grid.ra0, grid.dec0)
    col, row = grid.to_index(xi, eta)
    if radius_model is None:
        radius = np.clip(fwhm_arcsec * (1.0 + 0.35 * (19.0 - gaia.g[sel])), 1.2 * fwhm_arcsec, 30.0)
    else:
        radius = np.interp(gaia.g[sel], *radius_model)
        radius = np.maximum(radius, fwhm_arcsec)
    radius = np.minimum(radius + pad, 40.0)
    mask = np.zeros((grid.n, grid.n), bool)
    for c, r, rad in zip(col, row, radius / grid.step):
        r0, r1 = int(max(0, np.floor(r - rad))), int(min(grid.n, np.ceil(r + rad) + 1))
        c0, c1 = int(max(0, np.floor(c - rad))), int(min(grid.n, np.ceil(c + rad) + 1))
        if r0 >= r1 or c0 >= c1:
            continue
        yy, xx = np.mgrid[r0:r1, c0:c1]
        mask[r0:r1, c0:c1] |= (xx - c) ** 2 + (yy - r) ** 2 <= rad**2
    return mask


def inpaint(img: np.ndarray, mask: np.ndarray, sigma_px: float) -> np.ndarray:
    """Fill masked pixels from their surroundings by normalised convolution,
    repeated with a growing kernel until even the largest holes are filled."""
    out = img.copy()
    good = (~mask) & np.isfinite(img)
    filled = good.copy()
    sigma = sigma_px
    while not filled.all() and sigma < 64 * sigma_px:
        num = ndimage.gaussian_filter(np.where(good, img, 0.0), sigma)
        den = ndimage.gaussian_filter(good.astype(float), sigma)
        fill_here = (~filled) & (den > 0.05)
        out[fill_here] = num[fill_here] / den[fill_here]
        filled |= fill_here
        sigma *= 1.6
    out[~filled] = np.nanmedian(img[good])
    return out


def rank_transform(img: np.ndarray, region: np.ndarray) -> np.ndarray:
    """Each pixel's rank (0-1) among the pixels of ``region``; pixels outside
    the region are ranked against it too (clipped to 0-1)."""
    ref = np.sort(img[region & np.isfinite(img)])
    out = np.searchsorted(ref, img) / len(ref)
    return np.clip(out, 0.0, 1.0)


@dataclass
class Prepared:
    """An image ready for correlation, and what was done to it."""

    key: str
    epoch_map: np.ndarray
    image: np.ndarray  # rank-transformed, high-passed
    smooth: np.ndarray  # rank-transformed, before high-pass (for figures)
    stars: np.ndarray  # the star mask used (already widened by half the blur)
    valid: np.ndarray  # pixels the correlation may use: no star light, inside the image


def prepare(reg, grid: A.Grid, gaia: A.GaiaStars, region: np.ndarray, *,
            fwhm: float = config.COMMON_FWHM_ARCSEC, highpass: float = config.HIGHPASS_ARCSEC,
            asinh_scale: float | None = None) -> Prepared:
    """Steps 1-2 of the module docstring for one registered image.

    ``asinh_scale`` compresses a linear (CCD) image before blurring, so a few
    very bright knots do not dominate the blurred light the way they never did
    on a saturating plate; the rank transform makes the exact choice
    unimportant."""
    img = reg.image.astype(float)
    radius_model = star_radii(reg, grid, gaia)
    if asinh_scale:
        img = np.arcsinh(img / asinh_scale)
    stars = star_mask(grid, gaia, reg.epoch, max(reg.fwhm_arcsec, 1.5), radius_model=radius_model,
                      pad=0.5 * fwhm)
    filled = inpaint(img, stars | ~np.isfinite(img), sigma_px=2.0 / grid.step)
    blur = np.sqrt(max(fwhm**2 - reg.fwhm_arcsec**2, 0.0)) / 2.3548 / grid.step
    filled = ndimage.gaussian_filter(filled, blur) if blur > 0 else filled
    ranked = rank_transform(filled, region)
    hp = ranked - ndimage.gaussian_filter(ranked, highpass / 2.3548 / grid.step)
    valid = ~stars & np.isfinite(reg.image)
    return Prepared(reg.plate.key, reg.epoch_map, hp, ranked, stars, valid)


# --------------------------------------------------------------------------
# 3. Tracking patches
# --------------------------------------------------------------------------

def _subpixel_peak(cc: np.ndarray, iy: int, ix: int) -> tuple[float, float, float]:
    """Peak position to a fraction of a pixel: a 2-D quadratic through the
    3 x 3 neighbourhood of the best integer shift. Returns (dy, dx, curvature),
    where curvature is the mean second derivative (negative at a maximum)."""
    if iy < 1 or ix < 1 or iy >= cc.shape[0] - 1 or ix >= cc.shape[1] - 1:
        return np.nan, np.nan, np.nan
    z = cc[iy - 1 : iy + 2, ix - 1 : ix + 2]
    dxx = z[1, 2] - 2 * z[1, 1] + z[1, 0]
    dyy = z[2, 1] - 2 * z[1, 1] + z[0, 1]
    dxy = (z[2, 2] - z[2, 0] - z[0, 2] + z[0, 0]) / 4.0
    gx = (z[1, 2] - z[1, 0]) / 2.0
    gy = (z[2, 1] - z[0, 1]) / 2.0
    H = np.array([[dyy, dxy], [dxy, dxx]])
    det = np.linalg.det(H)
    if det <= 0 or dxx >= 0 or dyy >= 0:
        return np.nan, np.nan, np.nan
    off = -np.linalg.solve(H, [gy, gx])
    if np.any(np.abs(off) > 1.0):
        return np.nan, np.nan, np.nan
    return float(off[0]), float(off[1]), float(0.5 * (dxx + dyy))


def masked_ncc(tpl: np.ndarray, tmask: np.ndarray, win: np.ndarray, wmask: np.ndarray, *,
               min_overlap: float = 0.6) -> np.ndarray:
    """Normalised cross-correlation of a template inside a window, ignoring
    masked pixels in either (Padfield 2012, IEEE Trans. Image Proc. 21, 2706).

    For every placement s of the template in the window, the correlation
    coefficient is computed over only the pixels valid in both; each needed
    sum is itself a correlation, done with FFTs. Placements where fewer than
    ``min_overlap`` of the template's pixels overlap valid window pixels get
    -1. Returns an array of shape window - template + 1, like
    ``skimage.feature.match_template``."""
    from scipy.signal import fftconvolve

    mt = tmask.astype(float)
    mw = wmask.astype(float)
    t = np.where(tmask, tpl, 0.0)
    w = np.where(wmask, win, 0.0)

    def corr(a, b):  # sum_x a(x) b(x + s) for every valid placement s
        return fftconvolve(b, a[::-1, ::-1], mode="valid")

    n = corr(mt, mw)
    st = corr(t, mw)
    sw = corr(mt, w)
    stt = corr(t * t, mw)
    sww = corr(mt, w * w)
    stw = corr(t, w)
    with np.errstate(invalid="ignore", divide="ignore"):
        num = stw - st * sw / n
        den = np.sqrt(np.clip(stt - st**2 / n, 0, None) * np.clip(sww - sw**2 / n, 0, None))
        cc = num / den
    cc[(n < min_overlap * tmask.size) | ~np.isfinite(cc)] = -1.0
    return np.clip(cc, -1.0, 1.0)


def _match(template: np.ndarray, search: np.ndarray, tmask=None, smask=None):
    """Best match of ``template`` inside ``search``: (dy, dx) of the match's
    corner relative to the search window's centre placement, peak value,
    curvature, and the best correlation outside the main peak."""
    if tmask is None:
        cc = match_template(search, template)
    else:
        cc = masked_ncc(template, tmask, search, smask)
    iy, ix = np.unravel_index(np.argmax(cc), cc.shape)
    peak = float(cc[iy, ix])
    fy, fx, curv = _subpixel_peak(cc, iy, ix)
    cy, cx = (cc.shape[0] - 1) / 2.0, (cc.shape[1] - 1) / 2.0
    # Second-best peak: mask a disc around the best one (2 FWHM-ish radius).
    yy, xx = np.ogrid[: cc.shape[0], : cc.shape[1]]
    away = (yy - iy) ** 2 + (xx - ix) ** 2 > (template.shape[0] / 4.0) ** 2
    second = float(cc[away].max()) if away.any() else -1.0
    return iy + fy - cy, ix + fx - cx, peak, curv, second


@dataclass
class Vectors:
    """Motion vectors between two images, in standard coordinates (arcsec)."""

    pair: tuple[str, str]
    xi: np.ndarray  # mid-epoch position
    eta: np.ndarray
    dxi: np.ndarray  # displacement from image 1 to image 2 (arcsec, east positive)
    deta: np.ndarray  # (north positive)
    t1: np.ndarray  # epoch of image 1 at the patch (Julian year)
    t2: np.ndarray
    peak: np.ndarray  # normalised cross-correlation at the match
    second: np.ndarray  # best correlation away from the match
    back_err: np.ndarray  # arcsec: where the reverse match lands relative to the start

    @property
    def dt(self) -> np.ndarray:
        return self.t2 - self.t1

    @property
    def tmid(self) -> np.ndarray:
        return 0.5 * (self.t1 + self.t2)

    @property
    def mu_xi(self) -> np.ndarray:  # arcsec / yr
        return self.dxi / self.dt

    @property
    def mu_eta(self) -> np.ndarray:
        return self.deta / self.dt

    def subset(self, keep: np.ndarray) -> "Vectors":
        return Vectors(self.pair, *(getattr(self, f)[keep] for f in
                                    ("xi", "eta", "dxi", "deta", "t1", "t2", "peak", "second", "back_err")))

    def to_dict(self) -> dict:
        return {f: getattr(self, f) for f in ("xi", "eta", "dxi", "deta", "t1", "t2", "peak", "second", "back_err")}


def track(im1: Prepared, im2: Prepared, grid: A.Grid, points_mask: np.ndarray, *,
          patch: float = config.PATCH_ARCSEC, step: float = config.PATCH_STEP_ARCSEC,
          search: float = config.SEARCH_ARCSEC, min_std: float = 0.0, masked: bool = True,
          min_valid: float = 0.6) -> Vectors:
    """Track square patches of image 1 into image 2 (step 3 of the docstring).

    Patch centres lie on a regular grid inside ``points_mask``. Each patch is
    matched over shifts up to ``search`` arcsec in each direction; the match
    is then reversed (image 2 -> image 1) from where it landed. With
    ``masked`` (the default) star pixels are left out of every correlation
    rather than filled in, so stars that sit still cannot pull the match
    towards zero shift; patches with less than ``min_valid`` of their pixels
    free of stars are skipped."""
    half = int(round(patch / grid.step / 2))
    srch = int(round(search / grid.step))
    stride = max(1, int(round(step / grid.step)))
    a, b = im1.image, im2.image
    n = grid.n
    rows, cols = np.mgrid[half + srch : n - half - srch : stride, half + srch : n - half - srch : stride]
    rows, cols = rows.ravel(), cols.ravel()
    keep = points_mask[rows, cols]
    rows, cols = rows[keep], cols[keep]
    out = {k: [] for k in ("r", "c", "dy", "dx", "peak", "second", "back")}
    va, vb = im1.valid, im2.valid
    cut = lambda arr, r, c, h: arr[r - h : r + h, c - h : c + h]
    for r, c in zip(rows, cols):
        tpl = cut(a, r, c, half)
        if tpl.std() <= min_std:
            continue
        if masked:
            tm = cut(va, r, c, half)
            if tm.mean() < min_valid:
                continue
            dy, dx, peak, curv, second = _match(tpl, cut(b, r, c, half + srch), tm, cut(vb, r, c, half + srch))
        else:
            dy, dx, peak, curv, second = _match(tpl, cut(b, r, c, half + srch))
        if not np.isfinite(dy):
            continue
        # Reverse: the patch of image 2 at the match, searched for in image 1.
        r2, c2 = int(round(r + dy)), int(round(c + dx))
        if r2 - half - srch < 0 or c2 - half - srch < 0 or r2 + half + srch > n or c2 + half + srch > n:
            continue
        if masked:
            tm2 = cut(vb, r2, c2, half)
            if tm2.mean() < min_valid:
                continue
            bdy, bdx, *_ = _match(cut(b, r2, c2, half), cut(a, r2, c2, half + srch), tm2, cut(va, r2, c2, half + srch))
        else:
            bdy, bdx, *_ = _match(cut(b, r2, c2, half), cut(a, r2, c2, half + srch))
        back = np.hypot((r2 + bdy) - r, (c2 + bdx) - c) * grid.step
        out["r"].append(r), out["c"].append(c), out["dy"].append(dy), out["dx"].append(dx)
        out["peak"].append(peak), out["second"].append(second), out["back"].append(back)
    o = {k: np.array(v, float) for k, v in out.items()}
    # Patch centre (row r, column c) is pixel r - 0.5 + 0.5 ... the template
    # spans [r - half, r + half), whose centre is at r - 0.5.
    r_mid = o["r"] - 0.5 + 0.5 * o["dy"]
    c_mid = o["c"] - 0.5 + 0.5 * o["dx"]
    cgrid = (n - 1) / 2.0
    xi = -(c_mid - cgrid) * grid.step
    eta = (r_mid - cgrid) * grid.step
    ri, ci = np.clip(np.round(o["r"]).astype(int), 0, n - 1), np.clip(np.round(o["c"]).astype(int), 0, n - 1)
    r2i = np.clip(np.round(o["r"] + o["dy"]).astype(int), 0, n - 1)
    c2i = np.clip(np.round(o["c"] + o["dx"]).astype(int), 0, n - 1)
    t1 = im1.epoch_map[ri, ci]
    t2 = im2.epoch_map[r2i, c2i]
    return Vectors((im1.key, im2.key), xi, eta, -o["dx"] * grid.step, o["dy"] * grid.step, t1, t2,
                   o["peak"], o["second"], o["back"])


def quality_cut(v: Vectors, *, min_peak: float = 0.5, min_margin: float = 0.05,
                max_back: float = 1.0) -> np.ndarray:
    """Vectors we trust: a strong, unique match that reverses consistently."""
    return (v.peak >= min_peak) & (v.peak - v.second >= min_margin) & (v.back_err <= max_back)


# --------------------------------------------------------------------------
# 4. Rewinding: the expansion centre and the convergence date
# --------------------------------------------------------------------------

@dataclass
class ExpansionFit:
    xi_c: float  # arcsec
    eta_c: float
    T: float  # convergence date (Julian year)
    n: int
    rms_mu: float  # arcsec/yr: scatter of the vectors about the model
    inliers: np.ndarray
    cov: np.ndarray  # covariance of (xi_c, eta_c, T) from the fit

    @property
    def errors(self) -> np.ndarray:
        return np.sqrt(np.diag(self.cov))


def _residuals(params, xi, eta, mux, mue, tmid):
    xc, ec, T = params
    age = tmid - T
    return np.concatenate([mux - (xi - xc) / age, mue - (eta - ec) / age])


def fit_expansion(v: Vectors, *, clip: float = 4.0, T0: float = 1100.0) -> ExpansionFit:
    """Least-squares fit of p_mid - c = mu (t_mid - T) to all vectors.

    This is the "fit" method of Wyckoff & Murray (1977), Nugent (1998) and
    Martin et al. (2025): proper motion against position is a straight line
    whose slope is 1/age and whose zero point is the centre. Outliers beyond
    ``clip`` times the robust scatter are rejected and the fit repeated."""
    keep = np.ones(len(v.xi), bool)
    p = np.array([0.0, 0.0, T0])
    for _ in range(10):
        res = optimize.least_squares(_residuals, p, args=(v.xi[keep], v.eta[keep], v.mu_xi[keep],
                                                          v.mu_eta[keep], v.tmid[keep]))
        p = res.x
        r = _residuals(p, v.xi, v.eta, v.mu_xi, v.mu_eta, v.tmid)
        r2 = np.hypot(r[: len(v.xi)], r[len(v.xi):])
        sigma = 1.4826 * np.median(r2[keep]) / np.sqrt(2 * np.log(2))
        new = r2 < clip * sigma * np.sqrt(2)
        if np.array_equal(new, keep):
            break
        keep = new
    J = res.jac
    dof = max(1, 2 * keep.sum() - 3)
    s2 = np.sum(res.fun**2) / dof
    cov = np.linalg.inv(J.T @ J) * s2
    return ExpansionFit(float(p[0]), float(p[1]), float(p[2]), int(keep.sum()), float(np.sqrt(s2)), keep, cov)


def closest_approach(v: Vectors, *, T_grid=None, centre=None) -> tuple[float, np.ndarray, np.ndarray]:
    """The other classical method (Trimble 1968): run every filament
    backwards along its own motion and find when they huddle closest.

    Returns the date of minimum spread and the rms spread (arcsec) on
    ``T_grid``. Noise in the motions inflates the spread at early times, which
    biases this date late -- the reason Martin et al. (2025) preferred the fit."""
    T_grid = np.arange(900.0, 1400.0, 1.0) if T_grid is None else T_grid
    spread = []
    for T in T_grid:
        x = v.xi - v.mu_xi * (v.tmid - T)
        y = v.eta - v.mu_eta * (v.tmid - T)
        if centre is None:
            spread.append(np.sqrt(np.mean((x - np.median(x)) ** 2 + (y - np.median(y)) ** 2)))
        else:
            spread.append(np.sqrt(np.mean((x - centre[0]) ** 2 + (y - centre[1]) ** 2)))
    spread = np.array(spread)
    return float(T_grid[np.argmin(spread)]), T_grid, spread


def per_vector_dates(v: Vectors, xi_c: float, eta_c: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Each vector's own convergence date from its radial motion, its
    distance from the centre (arcsec) and its position angle (deg, east of north)."""
    dx, dy = v.xi - xi_c, v.eta - eta_c
    r = np.hypot(dx, dy)
    mu_r = (v.mu_xi * dx + v.mu_eta * dy) / np.where(r > 0, r, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        T = v.tmid - r / mu_r
    pa = np.degrees(np.arctan2(dx, dy)) % 360.0
    return T, r, pa


def block_bootstrap(v: Vectors, fitter, *, block_arcsec: float = 30.0, n_boot: int = 300,
                    seed: int = 1) -> np.ndarray:
    """Bootstrap over independent blocks of sky: overlapping patches give
    correlated vectors, so resample whole blocks, not single vectors."""
    rng = np.random.default_rng(seed)
    bx = np.floor(v.xi / block_arcsec).astype(int)
    by = np.floor(v.eta / block_arcsec).astype(int)
    labels = bx * 10007 + by
    blocks = np.unique(labels)
    members = {b: np.flatnonzero(labels == b) for b in blocks}
    out = []
    for _ in range(n_boot):
        pick = rng.choice(blocks, size=len(blocks), replace=True)
        idx = np.concatenate([members[b] for b in pick])
        out.append(fitter(v.subset(idx)))
    return np.array(out)
