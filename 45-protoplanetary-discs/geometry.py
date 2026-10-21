"""Disc orientation: from sky offsets to radius and azimuth in the disc, and back.

A thin circular disc tilted by inclination ``i`` (0 = face-on) appears as an ellipse whose
long axis points along the position angle ``PA`` (east of north) and whose short axis is
shortened by cos i. "Deprojecting" undoes the tilt: stretch the short axis by 1/cos i
and every ring becomes a circle.

Azimuth convention (DSHARP II, Appendix A): theta = 90 deg points along the PA direction
of the major axis, theta = 0 along PA + 90 deg, and theta increases clockwise on the sky
(east to the left). For PA = 0, theta = 0 is east and theta = 90 is north.

Two ways of measuring the orientation, as in DSHARP II:

* **E, ellipse fits** -- for each ring (or gap), find the radius of peak (or minimum)
  brightness along spokes about one beam apart, map those points back to the sky, and fit
  an ellipse (centre, semi-major axis, cos i, PA) with the scatter about it as a free
  parameter. The disc's geometry is the weighted mean over its rings.
* **G, a 2-D Gaussian** -- for discs without well-resolved rings: fit an elliptical
  Gaussian convolved with the beam; the deconvolved axis ratio gives cos i.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import emcee
import numpy as np
from scipy import ndimage, optimize

from alma_images import Image


# ---------------------------------------------------------------------------
# Coordinate transforms
# ---------------------------------------------------------------------------
def disc_polar(x, y, x0: float, y0: float, inc_deg: float, pa_deg: float):
    """Sky offsets (arcsec; x east, y north) -> deprojected radius (arcsec), azimuth (deg)."""
    pa, cosi = np.radians(pa_deg), np.cos(np.radians(inc_deg))
    dx, dy = np.asarray(x) - x0, np.asarray(y) - y0
    y_maj = dx * np.sin(pa) + dy * np.cos(pa)  # along the PA direction
    x_min = (dx * np.cos(pa) - dy * np.sin(pa)) / cosi  # along PA + 90, deprojected
    return np.hypot(x_min, y_maj), np.degrees(np.arctan2(y_maj, x_min))


def sky_from_polar(r, theta_deg, x0: float, y0: float, inc_deg: float, pa_deg: float):
    """Inverse of :func:`disc_polar`."""
    pa, cosi, th = np.radians(pa_deg), np.cos(np.radians(inc_deg)), np.radians(theta_deg)
    x_min = np.asarray(r) * np.cos(th) * cosi
    y_maj = np.asarray(r) * np.sin(th)
    return x0 + x_min * np.cos(pa) + y_maj * np.sin(pa), y0 - x_min * np.sin(pa) + y_maj * np.cos(pa)


def sky_pa_of_azimuth(theta_deg, pa_deg: float, inc_deg: float):
    """Position angle (deg E of N) on the sky of the direction at disc azimuth theta."""
    x, y = sky_from_polar(1.0, theta_deg, 0.0, 0.0, inc_deg, pa_deg)
    return np.degrees(np.arctan2(x, y)) % 360.0


def angle_diff_180(a, b):
    """Smallest difference between two axis orientations (period 180 deg)."""
    return (np.asarray(a) - np.asarray(b) + 90.0) % 180.0 - 90.0


_SPLINES: dict[int, tuple[np.ndarray, np.ndarray]] = {}


def sample_image(img: Image, x, y, order: int = 3):
    """Interpolate the image at sky offsets (arcsec) with a cubic spline. The spline
    coefficients of each image are computed once and cached (the expensive step)."""
    col = (np.asarray(x) - img.x[0]) / (img.x[1] - img.x[0])
    row = (np.asarray(y) - img.y[0]) / (img.y[1] - img.y[0])
    if order != 3:
        return ndimage.map_coordinates(img.data, [row, col], order=order, mode="constant", cval=np.nan)
    cached = _SPLINES.get(id(img))
    if cached is None or cached[0] is not img.data:
        cached = (img.data, ndimage.spline_filter(img.data, order=3))
        _SPLINES[id(img)] = cached
    return ndimage.map_coordinates(cached[1], [row, col], order=3, mode="constant", cval=np.nan, prefilter=False)


@dataclass
class DiscGeometry:
    x0: float  # arcsec, east offset of the disc centre from the phase centre
    y0: float
    inc: float  # deg
    pa: float  # deg
    method: str = ""
    x0_err: float = np.nan
    y0_err: float = np.nan
    inc_err: float = np.nan
    pa_err: float = np.nan
    n_features: int = 0

    def polar(self, x, y):
        return disc_polar(x, y, self.x0, self.y0, self.inc, self.pa)

    def sky(self, r, theta):
        return sky_from_polar(r, theta, self.x0, self.y0, self.inc, self.pa)

    def as_dict(self) -> dict:
        return {k: (float(v) if isinstance(v, (float, np.floating)) else v) for k, v in asdict(self).items()}


# ---------------------------------------------------------------------------
# G: a 2-D Gaussian convolved with the beam
# ---------------------------------------------------------------------------
def _cov(sig_maj: float, sig_min: float, pa_deg: float) -> np.ndarray:
    pa = np.radians(pa_deg)
    e1 = np.array([np.sin(pa), np.cos(pa)])
    e2 = np.array([np.cos(pa), -np.sin(pa)])
    return sig_maj**2 * np.outer(e1, e1) + sig_min**2 * np.outer(e2, e2)


FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))


@dataclass
class GaussianFit:
    x0: float
    y0: float
    peak: float
    fwhm_maj: float  # deconvolved, arcsec
    fwhm_min: float
    pa: float
    inc: float
    errors: dict = field(default_factory=dict)


def fit_gaussian(img: Image, x0: float, y0: float, radius: float) -> GaussianFit:
    """Fit an elliptical Gaussian, convolved with the synthesized beam, to the pixels within
    ``radius`` of (x0, y0), sampled one beam apart so that the residuals are independent."""
    beam_cov = _cov(img.bmaj / FWHM, img.bmin / FWHM, img.bpa)
    step = max(1, int(round(img.beam_mean / img.pix)))
    xx, yy = img.grids()
    sel = np.zeros_like(img.data, dtype=bool)
    sel[::step, ::step] = True
    sel &= np.hypot(xx - x0, yy - y0) < radius
    xs, ys, zs = xx[sel], yy[sel], img.data[sel]

    def model(p):
        cx, cy, amp, ls1, ls2, pa = p
        c = _cov(np.exp(ls1), np.exp(ls2), pa) + beam_cov
        ci = np.linalg.inv(c)
        dx, dy = xs - cx, ys - cy
        return amp * np.exp(-0.5 * (ci[0, 0] * dx * dx + 2 * ci[0, 1] * dx * dy + ci[1, 1] * dy * dy))

    # starting guess from intensity-weighted moments of the bright pixels
    w = np.clip(zs, 0, None) * (zs > 0.1 * zs.max())
    mx, my = np.average(xs, weights=w), np.average(ys, weights=w)
    cxx = np.average((xs - mx) ** 2, weights=w)
    cyy = np.average((ys - my) ** 2, weights=w)
    cxy = np.average((xs - mx) * (ys - my), weights=w)
    evals, evecs = np.linalg.eigh(np.array([[cxx, cxy], [cxy, cyy]]))
    pa0 = np.degrees(np.arctan2(evecs[0, 1], evecs[1, 1])) % 180.0
    s1, s2 = np.sqrt(max(evals[1], 1e-6)), np.sqrt(max(evals[0], 1e-6))
    p0 = [mx, my, zs.max(), np.log(s1), np.log(s2), pa0]
    res = optimize.least_squares(lambda p: model(p) - zs, p0, x_scale="jac")
    cx, cy, amp, ls1, ls2, pa = res.x
    s_maj, s_min = np.exp(ls1), np.exp(ls2)
    if s_min > s_maj:
        s_maj, s_min, pa = s_min, s_maj, pa + 90.0
    # formal errors from the Jacobian, scaled by the residual variance
    dof = max(1, zs.size - len(p0))
    s2 = np.sum(res.fun**2) / dof
    try:
        cov = np.linalg.inv(res.jac.T @ res.jac) * s2
        err = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        err = np.full(len(p0), np.nan)
    ratio = s_min / s_maj
    inc = np.degrees(np.arccos(ratio))
    # error on the inclination from those on the two widths
    d_ratio = ratio * np.hypot(err[3], err[4])
    inc_err = np.degrees(d_ratio / max(np.sqrt(1 - ratio**2), 1e-3))
    return GaussianFit(cx, cy, amp, s_maj * FWHM, s_min * FWHM, pa % 180.0, inc,
                       {"x0": err[0], "y0": err[1], "pa": err[5], "inc": inc_err})


# ---------------------------------------------------------------------------
# E: ellipses fitted to the ridge (ring) or trough (gap) of an annular feature
# ---------------------------------------------------------------------------
def ellipse_distance(params, xs, ys):
    """First-order (Sampson) orthogonal distance of points from an ellipse.

    params = (x0, y0, a, cos i, PA). The deprojected radius r_d of each point is compared
    with the semi-major axis a; dividing by |grad r_d| turns the radial mismatch into a
    distance on the sky (along the minor axis a radial offset shrinks by cos i). Exact to
    first order in (r_d - a)/a, which is < 1% for ring points.
    """
    p = np.asarray(params, dtype=float)
    if p.ndim == 2:  # one row of parameters per walker: broadcast over points
        x0, y0, a, cosi, pa = (p[:, k, None] for k in range(5))
    else:
        x0, y0, a, cosi, pa = p[:5]
    pa = np.radians(pa)
    dx, dy = xs - x0, ys - y0
    y_maj = dx * np.sin(pa) + dy * np.cos(pa)
    x_min = (dx * np.cos(pa) - dy * np.sin(pa)) / cosi
    rd = np.hypot(x_min, y_maj)
    grad = np.hypot(x_min / (rd * cosi), y_maj / rd)
    return (rd - a) / grad


@dataclass
class RidgePoints:
    theta: np.ndarray  # azimuths sampled (deg)
    r: np.ndarray  # radius of the extremum along each spoke (arcsec, in the trial geometry)
    x: np.ndarray  # sky position of the extremum
    y: np.ndarray
    n_tried: int


def extract_ridge(img: Image, geom: DiscGeometry, kind: str, r_lo: float, r_hi: float, dtheta_deg: float,
                  noise: float, exclude: tuple = (), min_contrast: float = 3.0,
                  step_frac: float = 0.1) -> RidgePoints:
    """Locate a ring's brightness maximum (kind 'B') or a gap's minimum ('D') along spokes.

    Each spoke runs from r_lo to r_hi (arcsec, deprojected in the trial geometry ``geom``);
    the window should stop short of the neighbouring features. Spokes whose extremum sits
    at the window's edge, or stands out by less than ``min_contrast`` times the noise, are
    dropped (DSHARP II: "if an extremum cannot be identified ... that azimuthal angle is
    excluded").
    """
    thetas = np.arange(-180.0, 180.0, dtheta_deg)
    rs = np.arange(max(r_lo, 0.3 * img.beam_mean), r_hi + 1e-9, step_frac * img.pix)
    out_t, out_r, n_tried = [], [], 0
    for th in thetas:
        if any(lo < th < hi for lo, hi in exclude):
            continue
        n_tried += 1
        xs, ys = geom.sky(rs, th)
        v = sample_image(img, xs, ys)
        if not np.all(np.isfinite(v)):
            continue
        k = int(np.argmax(v) if kind == "B" else np.argmin(v))
        if k == 0 or k == rs.size - 1:
            continue  # no extremum inside the window
        ends = min(v[0], v[-1]) if kind == "B" else max(v[0], v[-1])
        contrast = (v[k] - ends) if kind == "B" else (ends - v[k])
        if contrast < min_contrast * noise:
            continue
        y1, y2, y3 = v[k - 1], v[k], v[k + 1]
        denom = y1 - 2 * y2 + y3
        shift = 0.5 * (y1 - y3) / denom if denom != 0 else 0.0
        out_t.append(th)
        out_r.append(rs[k] + np.clip(shift, -1, 1) * (rs[1] - rs[0]))
    t, r = np.array(out_t), np.array(out_r)
    x, y = geom.sky(r, t)
    return RidgePoints(t, r, x, y, n_tried)


@dataclass
class EllipseFit:
    x0: float
    y0: float
    a: float  # semi-major axis, arcsec
    inc: float
    pa: float
    sigma: float  # scatter of the points about the ellipse, arcsec
    errors: dict
    n_points: int
    n_tried: int
    samples: np.ndarray | None = None


def fit_ellipse(pts: RidgePoints, guess: DiscGeometry, a_guess: float, pa_range: tuple[float, float],
                rng: np.random.Generator, nwalkers: int = 32, nsteps: int = 2000, burn: int = 600,
                inc_prior: tuple[float, float] | None = None, pa_prior: tuple[float, float] | None = None
                ) -> EllipseFit:
    """MCMC fit of (x0, y0, a, cos i, PA, ln sigma^2) as in DSHARP II eq. (1): orthogonal
    distances d_i, likelihood prod N(d_i; 0, sigma). Uniform priors; PA restricted to
    ``pa_range`` (the axis is only defined modulo 180 deg)."""
    xs, ys = pts.x, pts.y
    lo_pa, hi_pa = pa_range

    def log_prob(p):
        x0, y0, a, cosi, pa, lns2 = p.T
        ok = (cosi > 0.0) & (cosi < 1.0) & (pa >= lo_pa) & (pa < hi_pa) & (a > 0) & (lns2 > -25) & (lns2 < -2)
        out = np.full(p.shape[0], -np.inf)
        if not np.any(ok):
            return out
        q = p[ok]
        d = ellipse_distance(q[:, :5], xs, ys)
        s2 = np.exp(q[:, 5])[:, None]
        lp = -0.5 * np.sum(d**2 / s2 + np.log(2 * np.pi * s2), axis=1)
        if inc_prior is not None:
            inc = np.degrees(np.arccos(q[:, 3]))
            lp += -0.5 * ((inc - inc_prior[0]) / inc_prior[1]) ** 2
        if pa_prior is not None:
            lp += -0.5 * (angle_diff_180(q[:, 4], pa_prior[0]) / pa_prior[1]) ** 2
        out[ok] = lp
        return out

    pa0 = lo_pa + ((guess.pa - lo_pa) % 180.0)
    start = np.array([guess.x0, guess.y0, a_guess, np.cos(np.radians(guess.inc)), pa0])
    res = optimize.least_squares(lambda p: ellipse_distance(p, xs, ys), start,
                                 bounds=([-np.inf, -np.inf, 0, 1e-3, lo_pa - 90], [np.inf, np.inf, np.inf, 0.9999,
                                                                                    hi_pa + 90]))
    best = res.x.copy()
    best[4] = lo_pa + ((best[4] - lo_pa) % 180.0)
    s2 = max(np.mean(res.fun**2), 1e-10)
    ndim = 6
    p0 = np.concatenate([best, [np.log(s2)]])
    scale = np.array([1e-4, 1e-4, 1e-4, 1e-3, 0.2, 0.05])
    walkers = p0 + scale * rng.standard_normal((nwalkers, ndim))
    walkers[:, 3] = np.clip(walkers[:, 3], 1e-3, 0.999)
    walkers[:, 4] = np.clip(walkers[:, 4], lo_pa + 1e-6, hi_pa - 1e-6)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, vectorize=True)
    sampler.random_state = np.random.RandomState(int(rng.integers(2**31)))  # emcee keeps its own RNG
    sampler.run_mcmc(walkers, nsteps, progress=False)
    ch = sampler.get_chain(discard=burn, flat=True)
    inc = np.degrees(np.arccos(ch[:, 3]))
    pct = lambda v: np.percentile(v, [16, 50, 84])  # noqa: E731
    stats = {name: pct(v) for name, v in (("x0", ch[:, 0]), ("y0", ch[:, 1]), ("a", ch[:, 2]), ("inc", inc),
                                         ("pa", ch[:, 4]), ("sigma", np.exp(0.5 * ch[:, 5])))}
    errors = {k: 0.5 * (v[2] - v[0]) for k, v in stats.items()}
    return EllipseFit(stats["x0"][1], stats["y0"][1], stats["a"][1], stats["inc"][1], stats["pa"][1],
                      stats["sigma"][1], errors, len(xs), pts.n_tried, ch[:: max(1, len(ch) // 2000)])


def weighted_geometry(fits_: list[EllipseFit], pa_range: tuple[float, float]) -> DiscGeometry:
    """Inverse-variance weighted mean of several rings' ellipses (DSHARP II's E geometry)."""

    def wmean(vals, errs):
        w = 1.0 / np.asarray(errs) ** 2
        return float(np.sum(w * vals) / np.sum(w)), float(1.0 / np.sqrt(np.sum(w)))

    lo = pa_range[0]
    pas = np.array([lo + ((f.pa - lo) % 180.0) for f in fits_])
    x0 = wmean([f.x0 for f in fits_], [f.errors["x0"] for f in fits_])
    y0 = wmean([f.y0 for f in fits_], [f.errors["y0"] for f in fits_])
    inc = wmean([f.inc for f in fits_], [f.errors["inc"] for f in fits_])
    pa = wmean(pas, [f.errors["pa"] for f in fits_])
    return DiscGeometry(x0[0], y0[0], inc[0], pa[0], "E", x0[1], y0[1], inc[1], pa[1], len(fits_))
