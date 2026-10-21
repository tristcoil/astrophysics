"""Put old photographic plates and modern images on one Gaia-defined sky grid.

A Schmidt plate's geometry is not a perfect tangent-plane projection: the
plate was bent over a curved focal surface, the emulsion shifts, and the
scanner adds its own distortions. The header WCS that comes with a DSS scan is
good to about an arcsecond -- far too coarse when a filament moves only 6" in
40 years. So every image gets its own plate solution here:

1. find the stars in it and centroid them;
2. take Gaia DR3 stars and move each one along its own proper motion back to
   the date the image was taken (a star moving 20 mas/yr has moved 1.3"
   between 1951 and Gaia's 2016 reference epoch);
3. fit a cubic polynomial from pixel coordinates to the tangent plane, with
   outlier rejection, and record the residuals.

Every image is then resampled onto the same north-up grid of standard
coordinates (xi, eta) around a fixed tangent point, so that after
registration a pixel means the same piece of sky in 1951 and in 2012.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from datetime import datetime, timedelta

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage

import config

ARCSEC = 1.0 / 3600.0


# --------------------------------------------------------------------------
# Epochs
# --------------------------------------------------------------------------

def julian_year(when: datetime) -> float:
    """Julian epoch of a calendar date-time: J2000.0 + days / 365.25.

    Plate dates are UT; the ~1-minute difference between UT and the
    Terrestrial Time that defines Julian epochs is irrelevant here. (astropy's
    Time objects warn that UTC is undefined before 1960, hence plain
    arithmetic.)"""
    return 2000.0 + (when - datetime(2000, 1, 1, 12, 0, 0)).total_seconds() / 86400.0 / 365.25


def dss_epoch(header: fits.Header) -> float:
    """Date of a DSS plate, as a Julian year.

    DSS headers sometimes write minutes as '60' ('1951-11-07T09:60:00',
    the POSS-I red Crab plate), which no date parser accepts, so the fields
    are parsed by hand and added up. The exposure (minutes) is added to the
    start time's half-way point.
    """
    raw = header["DATE-OBS"].strip()
    m = re.match(r"(\d{4})-(\d{2})-(\d{2})(?:T(\d{1,2}):(\d{1,2})(?::(\d{1,2}(?:\.\d*)?))?)?", raw)
    if not m:
        raise ValueError(f"cannot parse DATE-OBS {raw!r}")
    y, mo, d = (int(m.group(i)) for i in (1, 2, 3))
    hh = int(m.group(4) or 0)
    mm = int(m.group(5) or 0)
    ss = float(m.group(6) or 0.0)
    start = datetime(y, mo, d) + timedelta(seconds=hh * 3600 + mm * 60 + ss)
    mid = start + timedelta(minutes=0.5 * float(header.get("EXPOSURE", 0.0)))
    return julian_year(mid)


def ps1_epoch(header: fits.Header) -> tuple[float, float, float]:
    """Exposure-weighted mean date of a Pan-STARRS stack, with the first and last.

    A stack's header lists its input warps (``INP_nnnn``), whose names begin
    with the exposure's night, ``o5180g...`` meaning MJD 55180, and their
    exposure times (``EXP_nnnn``)."""
    mjds, weights = [], []
    for k in range(int(header["NINPUTS"])):
        name = header[f"INP_{k:04d}"]
        m = re.match(r"o(\d{4})g", name)
        if not m:
            raise ValueError(f"unexpected warp name {name!r}")
        mjds.append(50000.0 + int(m.group(1)) + 0.5)  # mid-night: the hour is not in the name
        weights.append(float(header[f"EXP_{k:04d}"]))
    mjds, weights = np.array(mjds), np.array(weights)
    mean = np.sum(mjds * weights) / np.sum(weights)
    to_year = lambda mjd: 2000.0 + (mjd - 51544.5) / 365.25
    return to_year(mean), to_year(mjds.min()), to_year(mjds.max())


# --------------------------------------------------------------------------
# Gaia stars at any epoch
# --------------------------------------------------------------------------

@dataclass
class GaiaStars:
    source_id: np.ndarray
    ra: np.ndarray  # deg at epoch 2016.0
    dec: np.ndarray
    pmra: np.ndarray  # mas/yr (pmra includes cos dec); NaN for 2-parameter solutions
    pmdec: np.ndarray
    parallax: np.ndarray
    g: np.ndarray
    bp_rp: np.ndarray
    ruwe: np.ndarray

    def __len__(self) -> int:
        return len(self.ra)

    def at_epoch(self, epoch: float) -> tuple[np.ndarray, np.ndarray]:
        """Positions wound along each star's proper motion to ``epoch``.

        Linear motion on the sphere is fine here: over 65 years the largest
        proper motions in these fields (~100 mas/yr) move a star by 6.5", and
        the curvature terms are ~1e-9 of that. Stars without a proper motion
        (2-parameter Gaia solutions) are returned as NaN."""
        dt = epoch - config.GAIA_EPOCH
        dec = self.dec + self.pmdec * dt / 3.6e6
        ra = self.ra + self.pmra * dt / 3.6e6 / np.cos(np.radians(self.dec))
        return ra, dec


def read_gaia(target: str) -> GaiaStars:
    """Read the VizieR TSV written by ``fetch_data.py``."""
    path = config.GAIA_DIR / f"{target}_gaia_dr3.tsv"
    lines = [ln for ln in path.read_text().splitlines() if ln and not ln.startswith("#")]
    header = lines[0].split("\t")
    rows = [ln.split("\t") for ln in lines[3:]]  # skip the units and dashes lines

    def col(name, dtype=float):
        i = header.index(name)
        out = []
        for r in rows:
            v = r[i].strip() if i < len(r) else ""
            out.append(dtype(v) if v else np.nan)
        return np.array(out)

    return GaiaStars(col("Source", int), col("RA_ICRS"), col("DE_ICRS"), col("pmRA"), col("pmDE"),
                     col("Plx"), col("Gmag"), col("BP-RP"), col("RUWE"))


# --------------------------------------------------------------------------
# Tangent plane
# --------------------------------------------------------------------------

def to_tangent(ra, dec, ra0: float, dec0: float):
    """Gnomonic (TAN) projection: standard coordinates xi (east), eta (north) in arcsec."""
    ra, dec = np.radians(ra), np.radians(dec)
    a0, d0 = np.radians(ra0), np.radians(dec0)
    cosc = np.sin(d0) * np.sin(dec) + np.cos(d0) * np.cos(dec) * np.cos(ra - a0)
    xi = np.cos(dec) * np.sin(ra - a0) / cosc
    eta = (np.cos(d0) * np.sin(dec) - np.sin(d0) * np.cos(dec) * np.cos(ra - a0)) / cosc
    return np.degrees(xi) * 3600.0, np.degrees(eta) * 3600.0


def from_tangent(xi, eta, ra0: float, dec0: float):
    """Inverse gnomonic projection: (xi, eta) in arcsec back to RA, Dec in degrees."""
    x, y = np.radians(np.asarray(xi) / 3600.0), np.radians(np.asarray(eta) / 3600.0)
    a0, d0 = np.radians(ra0), np.radians(dec0)
    rho = np.hypot(x, y)
    c = np.arctan(rho)
    with np.errstate(invalid="ignore", divide="ignore"):
        dec = np.arcsin(np.cos(c) * np.sin(d0) + np.where(rho > 0, y * np.sin(c) * np.cos(d0) / rho, 0.0))
    ra = a0 + np.arctan2(x * np.sin(c), rho * np.cos(d0) * np.cos(c) - y * np.sin(d0) * np.sin(c))
    return np.degrees(ra) % 360.0, np.degrees(dec)


# --------------------------------------------------------------------------
# Finding and centroiding stars
# --------------------------------------------------------------------------

def background(img: np.ndarray, size_px: int) -> np.ndarray:
    """A smooth background: median filter on a coarse grid, then bilinear up.

    Median-filtering the full image is slow; every ``step`` pixels is enough
    for a background that only has to follow the nebula's broad glow."""
    step = max(1, size_px // 4)
    filled = np.where(np.isfinite(img), img, np.nanmedian(img))
    small = filled[::step, ::step]
    small = ndimage.median_filter(small, size=max(3, size_px // step), mode="reflect")
    zoom = (img.shape[0] / small.shape[0], img.shape[1] / small.shape[1])
    return ndimage.zoom(small, zoom, order=1)[: img.shape[0], : img.shape[1]]


def robust_sigma(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def refine_centroids(img: np.ndarray, x0: np.ndarray, y0: np.ndarray, radius: float,
                     max_iterations: int = 40, tol: float = 1e-3) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Iterative, Gaussian-weighted centroids (a cheap stand-in for a PSF fit).

    Weighting by a Gaussian of the star's own size keeps neighbours and the
    nebula from pulling the centroid. For a Gaussian star of width sigma and
    a weight of width w, each iteration shrinks the error by sigma^2 / (sigma^2
    + w^2) -- a half when w = sigma -- so we iterate until the step is below
    ``tol`` pixels (stopping early would leave a fraction of the starting
    offset in every centroid). Returns x, y and a quality flag (False if the
    centroid ran away, did not converge, or the window left the image)."""
    ny, nx = img.shape
    r = int(np.ceil(2.5 * radius))
    yy, xx = np.mgrid[-r : r + 1, -r : r + 1]
    x, y = x0.astype(float).copy(), y0.astype(float).copy()
    ok = np.ones(len(x), bool)
    for i in range(len(x)):
        converged = False
        for _ in range(max_iterations):
            ix, iy = int(round(x[i])), int(round(y[i]))
            if ix - r < 0 or iy - r < 0 or ix + r >= nx or iy + r >= ny:
                break
            stamp = img[iy - r : iy + r + 1, ix - r : ix + r + 1]
            if not np.all(np.isfinite(stamp)):
                break
            stamp = stamp - np.median(np.concatenate([stamp[0], stamp[-1], stamp[:, 0], stamp[:, -1]]))
            dx, dy = xx - (x[i] - ix), yy - (y[i] - iy)
            w = np.exp(-0.5 * (dx * dx + dy * dy) / radius**2) * np.clip(stamp, 0, None)
            s = w.sum()
            if s <= 0:
                break
            nxp = ix + (w * xx).sum() / s
            nyp = iy + (w * yy).sum() / s
            step = np.hypot(nxp - x[i], nyp - y[i])
            x[i], y[i] = nxp, nyp
            if step < tol:
                converged = True
                break
        ok[i] = converged
    moved = np.hypot(x - x0, y - y0)
    ok &= moved < 2.0 * radius
    return x, y, ok


# --------------------------------------------------------------------------
# Plate solutions
# --------------------------------------------------------------------------

def _poly_terms(u: np.ndarray, v: np.ndarray, order: int) -> np.ndarray:
    """Design matrix of all monomials u^i v^j with i + j <= order."""
    cols = [u**i * v ** (k - i) for k in range(order + 1) for i in range(k, -1, -1)]
    return np.stack(cols, axis=-1)


@dataclass
class PlateSolution:
    """Pixel <-> standard coordinates for one image, from Gaia stars at its epoch."""

    key: str
    epoch: float
    ra0: float
    dec0: float
    order: int
    pix_center: tuple[float, float]
    pix_scale: float  # normalisation of pixel coordinates (not a plate scale)
    coef_xi: np.ndarray
    coef_eta: np.ndarray
    inv_order: int
    sky_scale: float
    coef_x: np.ndarray
    coef_y: np.ndarray
    n_stars: int
    rms_arcsec: float  # per-coordinate rms of the star residuals
    residuals: dict = field(default_factory=dict)  # star table: x, y, G, dxi, deta

    def to_sky(self, x, y):
        u = (np.asarray(x) - self.pix_center[0]) / self.pix_scale
        v = (np.asarray(y) - self.pix_center[1]) / self.pix_scale
        A = _poly_terms(u, v, self.order)
        return A @ self.coef_xi, A @ self.coef_eta

    def to_pixel(self, xi, eta):
        p = np.asarray(xi) / self.sky_scale
        q = np.asarray(eta) / self.sky_scale
        A = _poly_terms(p, q, self.inv_order)
        return A @ self.coef_x, A @ self.coef_y

    def pixel_scale_arcsec(self) -> float:
        """Mean linear pixel size at the image centre, from the Jacobian."""
        cx, cy = self.pix_center
        xi0, eta0 = self.to_sky(cx, cy)
        xi1, eta1 = self.to_sky(cx + 1, cy)
        xi2, eta2 = self.to_sky(cx, cy + 1)
        det = abs((xi1 - xi0) * (eta2 - eta0) - (xi2 - xi0) * (eta1 - eta0))
        return float(np.sqrt(det))


def _fit_clipped(A: np.ndarray, b1: np.ndarray, b2: np.ndarray, clip: float = 3.0, max_iter: int = 10):
    """Least squares for two right-hand sides sharing one design matrix,
    rejecting points whose 2-D residual exceeds ``clip`` times the robust rms."""
    keep = np.ones(len(b1), bool)
    for _ in range(max_iter):
        c1, *_ = np.linalg.lstsq(A[keep], b1[keep], rcond=None)
        c2, *_ = np.linalg.lstsq(A[keep], b2[keep], rcond=None)
        r = np.hypot(A @ c1 - b1, A @ c2 - b2)
        sigma = 1.4826 * np.median(r[keep]) / np.sqrt(2 * np.log(2))  # Rayleigh median -> sigma
        new = r < clip * max(sigma, 1e-6) * np.sqrt(2)
        if np.array_equal(new, keep):
            break
        keep = new
    return c1, c2, keep


def solve_plate(key: str, img: np.ndarray, wcs: WCS, epoch: float, gaia: GaiaStars,
                ra0: float, dec0: float, *, star_radius_px: float, order: int = config.PLATE_POLY_ORDER,
                gmin: float = config.REG_GMIN, gmax: float = config.REG_GMAX,
                match_arcsec: float = 4.0, exclude=None) -> PlateSolution:
    """Fit a polynomial plate solution to Gaia stars wound back to ``epoch``.

    The header ``wcs`` only has to be good to a few arcseconds: it predicts
    where each Gaia star should be, the star is centroided there, and the
    measured centroids are what the polynomial is fitted to.
    ``exclude(xi, eta) -> bool`` drops stars on top of the nebula, whose
    centroids the bright filaments could bias.
    """
    ra, dec = gaia.at_epoch(epoch)
    good = np.isfinite(ra) & (gaia.g > gmin) & (gaia.g < gmax)
    # Isolated stars only: no other Gaia star within 6" brighter than G + 3.
    xi_all, eta_all = to_tangent(gaia.ra, gaia.dec, ra0, dec0)
    from scipy.spatial import cKDTree

    tree = cKDTree(np.c_[xi_all, eta_all])
    crowded = np.zeros(len(gaia), bool)
    for i in np.flatnonzero(good):
        for j in tree.query_ball_point([xi_all[i], eta_all[i]], r=6.0):
            if j != i and gaia.g[j] < gaia.g[i] + 3.0:
                crowded[i] = True
                break
    good &= ~crowded
    if exclude is not None:
        good &= ~exclude(xi_all, eta_all)
    idx = np.flatnonzero(good)
    xp, yp = wcs.world_to_pixel_values(ra[idx], dec[idx])
    xp, yp = np.asarray(xp), np.asarray(yp)
    inside = (xp > 5) & (yp > 5) & (xp < img.shape[1] - 6) & (yp < img.shape[0] - 6)
    idx, xp, yp = idx[inside], xp[inside], yp[inside]

    bkg_free = img - background(img, int(8 * star_radius_px))
    # First pass: centroid with a generous window, then find the global offset.
    xc, yc, ok = refine_centroids(bkg_free, xp, yp, radius=star_radius_px * 1.5)
    off_x = np.median((xc - xp)[ok])
    off_y = np.median((yc - yp)[ok])
    xc, yc, ok = refine_centroids(bkg_free, xp + off_x, yp + off_y, radius=star_radius_px)
    pixel_arcsec = np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix))) * 3600.0
    ok &= np.hypot(xc - xp - off_x, yc - yp - off_y) * pixel_arcsec < match_arcsec
    idx, xc, yc = idx[ok], xc[ok], yc[ok]

    xi, eta = to_tangent(ra[idx], dec[idx], ra0, dec0)
    cx, cy = (img.shape[1] - 1) / 2.0, (img.shape[0] - 1) / 2.0
    s = max(img.shape) / 2.0
    A = _poly_terms((xc - cx) / s, (yc - cy) / s, order)
    c_xi, c_eta, keep = _fit_clipped(A, xi, eta)
    dxi = xi - A @ c_xi
    deta = eta - A @ c_eta
    rms = float(np.sqrt(0.5 * np.mean(dxi[keep] ** 2 + deta[keep] ** 2)))

    # Inverse map from a dense grid through the forward map (better conditioned
    # than fitting it to the stars alone); one order higher than the forward.
    gy, gx = np.mgrid[0 : img.shape[0] : 16, 0 : img.shape[1] : 16]
    gx, gy = gx.ravel().astype(float), gy.ravel().astype(float)
    Ag = _poly_terms((gx - cx) / s, (gy - cy) / s, order)
    gxi, geta = Ag @ c_xi, Ag @ c_eta
    sky_scale = float(np.max(np.abs(np.r_[gxi, geta])))
    inv_order = order + 2
    Ai = _poly_terms(gxi / sky_scale, geta / sky_scale, inv_order)
    c_x, *_ = np.linalg.lstsq(Ai, gx, rcond=None)
    c_y, *_ = np.linalg.lstsq(Ai, gy, rcond=None)
    inv_err = float(np.max(np.hypot(Ai @ c_x - gx, Ai @ c_y - gy)))
    if inv_err > 0.02:
        raise RuntimeError(f"{key}: inverse plate solution off by {inv_err:.3f} px")

    return PlateSolution(
        key=key, epoch=epoch, ra0=ra0, dec0=dec0, order=order, pix_center=(cx, cy), pix_scale=s,
        coef_xi=c_xi, coef_eta=c_eta, inv_order=inv_order, sky_scale=sky_scale, coef_x=c_x, coef_y=c_y,
        n_stars=int(keep.sum()), rms_arcsec=rms,
        residuals={"source_id": gaia.source_id[idx], "x": xc, "y": yc, "g": gaia.g[idx],
                   "bp_rp": gaia.bp_rp[idx], "dxi": dxi, "deta": deta, "used": keep})


# --------------------------------------------------------------------------
# The common grid
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Grid:
    """A square north-up grid of standard coordinates around (ra0, dec0).

    Pixel (i, j) of an array on this grid is at xi = -(j - c) * step (east is
    to the left, as on the sky) and eta = (i - c) * step, c = (n - 1) / 2."""

    ra0: float
    dec0: float
    step: float  # arcsec per pixel
    n: int

    @property
    def coords(self) -> tuple[np.ndarray, np.ndarray]:
        c = (self.n - 1) / 2.0
        j = np.arange(self.n)
        xi_1d = -(j - c) * self.step
        eta_1d = (j - c) * self.step
        return np.meshgrid(xi_1d, eta_1d)

    def to_index(self, xi, eta):
        """Fractional (column, row) on the grid of standard coordinates."""
        c = (self.n - 1) / 2.0
        return c - np.asarray(xi) / self.step, c + np.asarray(eta) / self.step

    def extent(self) -> list[float]:
        """imshow extent in arcsec, east to the left: [left, right, bottom, top]."""
        h = (self.n / 2.0) * self.step
        return [h, -h, -h, h]

    def wcs(self) -> WCS:
        w = WCS(naxis=2)
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.crval = [self.ra0, self.dec0]
        w.wcs.crpix = [(self.n + 1) / 2.0, (self.n + 1) / 2.0]
        w.wcs.cdelt = [-self.step / 3600.0, self.step / 3600.0]
        return w


def resample(img: np.ndarray, sol: PlateSolution, grid: Grid, order: int = 3) -> np.ndarray:
    """The image on the common grid (NaN where it has no data).

    Cubic-spline interpolation of a well-sampled image; NaNs are filled with
    the median before the spline and the footprint is resampled separately,
    so an edge does not ring into the data."""
    xi, eta = grid.coords
    x, y = sol.to_pixel(xi.ravel(), eta.ravel())
    finite = np.isfinite(img)
    filled = np.where(finite, img, np.nanmedian(img))
    out = ndimage.map_coordinates(filled.astype(float), [y, x], order=order, mode="constant", cval=np.nan)
    foot = ndimage.map_coordinates(finite.astype(float), [y, x], order=1, mode="constant", cval=0.0)
    out[foot < 0.999] = np.nan
    return out.reshape(grid.n, grid.n)


def stamp_fwhm(img: np.ndarray, x: np.ndarray, y: np.ndarray, half: int = 10) -> float:
    """Median FWHM (pixels) of point sources at (x, y), from second moments
    of background-subtracted stamps, clipped at zero."""
    out = []
    yy, xx = np.mgrid[-half : half + 1, -half : half + 1]
    for xc, yc in zip(x, y):
        ix, iy = int(round(xc)), int(round(yc))
        if ix - half < 0 or iy - half < 0 or ix + half >= img.shape[1] or iy + half >= img.shape[0]:
            continue
        st = img[iy - half : iy + half + 1, ix - half : ix + half + 1]
        if not np.all(np.isfinite(st)):
            continue
        st = st - np.median(np.concatenate([st[0], st[-1], st[:, 0], st[:, -1]]))
        st = np.clip(st, 0, None)
        # Restrict to the core: pixels above 10% of the peak, to ignore the wings.
        st = np.where(st > 0.1 * st.max(), st, 0.0)
        s = st.sum()
        if s <= 0:
            continue
        mx, my = (st * xx).sum() / s, (st * yy).sum() / s
        var = (st * ((xx - mx) ** 2 + (yy - my) ** 2)).sum() / s / 2.0
        out.append(2.3548 * np.sqrt(var))
    # The 10%-of-peak cut removes the Gaussian's wings and shrinks the
    # measured second moment. For a 2-D Gaussian cut at r^2 = 2 sigma^2 ln 10
    # the kept variance is sigma^2 (1 - 0.1 (1 + ln 10)) / 0.9 = 0.744 sigma^2.
    corr = 1.0 / np.sqrt((1.0 - 0.1 * (1.0 + np.log(10.0))) / 0.9)
    return float(np.median(out) * corr) if out else float("nan")
