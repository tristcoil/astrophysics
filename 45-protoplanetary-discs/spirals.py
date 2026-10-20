"""Spiral arms in the dust: how tightly are they wound?

Three DSHARP discs (Elias 27, IM Lup, WaOph 6) show two-armed spirals in the millimetre
continuum. A spiral's *pitch angle* mu is the angle between the arm and a circle:
tan(mu) = |dr / (r dtheta)|. For a logarithmic spiral, r = R0 exp(b theta), the pitch
angle is the same everywhere: mu = arctan|b|. Spirals launched by a planet are tightly
wound far from it; spirals from gravitational instability in a massive disc have pitch
angles of ~10-15 deg. (Galaxy arms are logarithmic spirals too, with mu ~ 5-30 deg.)

Two independent measurements:

1. **Tracing** (DSHARP III's method): deproject, subtract the azimuthal median at each
   radius (leaving the "non-axisymmetric residual"), and along rays of fixed azimuth find
   the residual's local maxima; link them into two arms and fit each with a log spiral.
2. **Log-spiral Fourier transform** (as used for galaxies, e.g. Puerari & Dottori 1992):
   with u = ln r, a log spiral is a straight line in (u, theta), so the 2-D Fourier
   amplitude A(m, p) = |sum I(u, theta) exp(-i (m theta + p u))| peaks at p = -m/b. For a
   two-armed (m = 2) pattern this gives the mean pitch angle without tracing anything.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage, optimize

from alma_images import Image
from geometry import DiscGeometry, sample_image


@dataclass
class PolarImage:
    r: np.ndarray  # arcsec
    theta: np.ndarray  # deg (DSHARP convention)
    image: np.ndarray  # (n_r, n_theta)
    residual: np.ndarray  # image minus the azimuthal median at each radius


def polar_image(img: Image, geom: DiscGeometry, r_max: float, dr: float, dtheta: float = 1.0) -> PolarImage:
    r = np.arange(dr / 2, r_max, dr)
    th = np.arange(-180.0, 180.0, dtheta)
    rr, tt = np.meshgrid(r, th, indexing="ij")
    x, y = geom.sky(rr, tt)
    im = sample_image(img, x, y)
    med = np.nanmedian(im, axis=1, keepdims=True)
    return PolarImage(r, th, im, im - med)


# ---------------------------------------------------------------------------
# 1. Tracing
# ---------------------------------------------------------------------------
@dataclass
class ArmPoints:
    theta: np.ndarray  # unwrapped polar angle, deg
    r: np.ndarray  # arcsec
    r_err: np.ndarray


def ray_maxima(pol: PolarImage, r_lo: float, r_hi: float, dtheta: float, min_prominence: float) -> list[tuple]:
    """Local maxima of the residual along rays spaced dtheta apart: (theta, r, height)."""
    out = []
    sel = (pol.r >= r_lo) & (pol.r <= r_hi)
    r = pol.r[sel]
    for th in np.arange(-180.0, 180.0, dtheta):
        j = int(np.argmin(np.abs(pol.theta - th)))
        prof = pol.residual[sel, j]
        if not np.all(np.isfinite(prof)):
            continue
        from scipy.signal import find_peaks

        pk, props = find_peaks(prof, prominence=min_prominence)
        for p in pk:
            if 0 < p < r.size - 1:
                y1, y2, y3 = prof[p - 1], prof[p], prof[p + 1]
                den = y1 - 2 * y2 + y3
                sh = 0.5 * (y1 - y3) / den if den != 0 else 0.0
                out.append((float(pol.theta[j]), float(r[p] + np.clip(sh, -1, 1) * (r[1] - r[0])), float(prof[p])))
    return out


def link_arm(points: list[tuple], b: float, theta0: float, r0: float, tol_lnr: float) -> ArmPoints:
    """Follow one arm: at each ray keep the maximum closest (in ln r) to the log spiral
    ln r = ln r0 + b (theta - theta0), trying the ray's angle in every winding."""
    by_theta: dict[float, list] = {}
    for th, r, h in points:
        by_theta.setdefault(th, []).append((r, h))
    t_out, r_out = [], []
    for th, cands in by_theta.items():
        best = None
        for wrap in (-720.0, -360.0, 0.0, 360.0, 720.0):
            tt = th + wrap
            pred = np.log(r0) + b * np.radians(tt - theta0)
            for r, h in cands:
                d = abs(np.log(r) - pred)
                if d < tol_lnr and (best is None or d < best[0]):
                    best = (d, tt, r)
        if best is not None:
            t_out.append(best[1])
            r_out.append(best[2])
    o = np.argsort(t_out)
    return ArmPoints(np.array(t_out)[o], np.array(r_out)[o], np.full(len(o), np.nan))


@dataclass
class SpiralFit:
    r0: float  # arcsec at theta = 0
    b: float
    pitch: float  # deg
    r0_err: float
    b_err: float
    pitch_err: float
    n: int
    theta_range: tuple[float, float]
    extra: dict = field(default_factory=dict)


def fit_log_spiral(arm: ArmPoints, sigma_r: np.ndarray) -> SpiralFit:
    """Weighted least squares of r = R0 exp(b theta), theta in radians (DSHARP III eq. A1)."""
    th = np.radians(arm.theta)

    def resid(p):
        return (arm.r - p[0] * np.exp(p[1] * th)) / sigma_r

    b0, lnr0 = np.polyfit(th, np.log(arm.r), 1)
    res = optimize.least_squares(resid, [np.exp(lnr0), b0])
    dof = max(1, arm.r.size - 2)
    s2 = max(1.0, np.sum(res.fun**2) / dof)  # inflate errors if the scatter exceeds the beam-based errors
    cov = np.linalg.inv(res.jac.T @ res.jac) * s2
    r0, b = res.x
    e_r0, e_b = np.sqrt(np.diag(cov))
    return SpiralFit(r0, b, np.degrees(np.arctan(abs(b))), e_r0, e_b, np.degrees(e_b / (1 + b * b)), int(arm.r.size),
                     (float(arm.theta.min()), float(arm.theta.max())), {"chi2_red": float(np.sum(res.fun**2) / dof)})


def radial_sigma(arm: ArmPoints, geom: DiscGeometry, beam: float) -> np.ndarray:
    """1-sigma radial uncertainty of a traced point: the beam's sigma, stretched by the
    deprojection at that azimuth (along the minor axis a sky offset s is s / cos i in the
    disc) -- as in DSHARP III."""
    sig = beam / 2.3548
    th = np.radians(arm.theta)
    cosi = np.cos(np.radians(geom.inc))
    # radial unit vector (x_min, y_maj) = (cos th, sin th); on the sky the x_min part shrinks by cos i
    stretch = 1.0 / np.hypot(np.cos(th) * cosi, np.sin(th))
    return sig * stretch


# ---------------------------------------------------------------------------
# 2. Logarithmic-spiral Fourier transform
# ---------------------------------------------------------------------------
def log_spiral_spectrum(pol: PolarImage, r_lo: float, r_hi: float, m: int = 2, p_grid=None):
    """A(m, p) of the non-axisymmetric residual between r_lo and r_hi.

    The residual is resampled uniformly in u = ln r and weighted by r^2 (the area element
    r dr dtheta = r^2 du dtheta), normalised by the total weight."""
    if p_grid is None:
        p_grid = np.linspace(-40.0, 40.0, 1601)
    u = np.linspace(np.log(r_lo), np.log(r_hi), 400)
    th = np.radians(pol.theta)
    res = np.nan_to_num(pol.residual)
    resid_u = np.array([np.interp(np.exp(u), pol.r, res[:, j]) for j in range(th.size)]).T  # (n_u, n_theta)
    w = np.exp(2 * u)[:, None] * np.clip(resid_u, 0, None)  # positive residual = arm emission
    tot = w.sum()
    phase_theta = np.exp(-1j * m * th)[None, :]
    a_theta = (w * phase_theta).sum(axis=1)  # (n_u,)
    amp = np.array([np.abs(np.sum(a_theta * np.exp(-1j * p * u))) / tot for p in p_grid])
    return p_grid, amp


def pitch_from_spectrum(p_grid, amp, m: int = 2) -> tuple[float, float, float]:
    """(p_max, b = -m/p_max, pitch angle in deg) at the strongest peak (|p| > 1)."""
    ok = np.abs(p_grid) > 1.0
    k = np.argmax(np.where(ok, amp, -np.inf))
    # refine with a parabola
    if 0 < k < amp.size - 1:
        y1, y2, y3 = amp[k - 1], amp[k], amp[k + 1]
        den = y1 - 2 * y2 + y3
        sh = 0.5 * (y1 - y3) / den if den != 0 else 0.0
    else:
        sh = 0.0
    p = p_grid[k] + sh * (p_grid[1] - p_grid[0])
    b = -m / p
    return float(p), float(b), float(np.degrees(np.arctan(abs(b))))
