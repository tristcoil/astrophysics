"""Weighing HD 163296 with its own disc: the CO velocity field.

Carbon monoxide gas emits at 230.538 GHz; the Doppler shift of that line maps the
line-of-sight velocity of the gas, pixel by pixel. A disc in Keplerian rotation has
v_phi = sqrt(G M* / r), so the velocity map weighs the star. Two complications:

* 12CO is optically thick, so we see it from an elevated "emission surface" z(r) -- the
  layer about 2-3 pressure scale heights up where the line becomes opaque -- not from
  the midplane. On a flared disc this bends the classic butterfly pattern. The surface is
  fitted together with the mass, as z(r) = z0 (r/1")^phi - z1 (r/1")^psi (Teague et al.
  2019), and the gas there orbits a little slower: v^2 = G M r^2 / (r^2 + z^2)^(3/2).
* The DSHARP cube has 0.32 km/s channels (0.64 km/s resolution after Hanning smoothing),
  coarse compared with the planet-induced deviations we would like to see (~0.1-0.3 km/s).
  The line centre in each pixel is therefore found by fitting a parabola through the
  brightest channel and its two neighbours (Teague & Foreman-Mackey 2018), which recovers
  centres to a fraction of a channel.

Geometry conventions. PA is the position angle (east of north) of the red-shifted half of
the major axis. In the disc frame, (r, phi) are cylindrical coordinates with phi = 0 along
that half; a point at (r, phi, z) appears on the sky at
    x_maj = r cos(phi),  y_min = r sin(phi) cos(i) + z sin(i),
along and across the major axis, and its line-of-sight velocity is
    v = v_sys + v_phi(r, z) cos(phi) sin(i).
with 0 < i < 90 deg and z measured towards the observer's side of the disc (the bright
"front" surface is at z = +h(r)). ``near`` = +1 puts the near half of the minor axis at
PA - 90 deg, -1 at PA + 90 deg (y_min is mirrored). Together PA (0-360) and ``near`` cover
all four ways a disc can sit on the sky -- which half of the major axis recedes, and which
half of the minor axis is closer to us -- and so fix its sense of rotation. Only the
elevated emission surface tells the two near sides apart: the front surface appears
shifted towards the far side of the disc.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import emcee
import numpy as np
from scipy import optimize

from alma_images import Cube
from config import AU_CGS, G_CGS, M_SUN_CGS


# ---------------------------------------------------------------------------
# Line centres
# ---------------------------------------------------------------------------
@dataclass
class LineMap:
    v0: np.ndarray  # km/s, NaN where masked
    dv0: np.ndarray  # 1-sigma, km/s
    peak: np.ndarray  # Jy/beam
    x: np.ndarray
    y: np.ndarray
    rms: float


def channel_rms(cube: Cube, n_edge: int = 8) -> float:
    """Noise from the line-free channels at both ends of the band."""
    edge = np.concatenate([cube.data[:n_edge].ravel(), cube.data[-n_edge:].ravel()])
    edge = edge[np.isfinite(edge)]
    return float(1.4826 * np.median(np.abs(edge - np.median(edge))))


def line_centres(cube: Cube, rms: float, v_lo: float, v_hi: float, clip: float = 5.0) -> LineMap:
    """Quadratic interpolation of each spectrum's brightest channel (bettermoments' method),
    restricted to channels v_lo..v_hi (the disc's velocity range)."""
    sel = np.where((cube.velocity >= v_lo) & (cube.velocity <= v_hi))[0]
    data = cube.data[sel[0]:sel[-1] + 1]
    vel = cube.velocity[sel[0]:sel[-1] + 1]
    dv = float(np.median(np.diff(vel)))
    k = np.nanargmax(np.nan_to_num(data, nan=-np.inf), axis=0)
    k = np.clip(k, 1, data.shape[0] - 2)
    rows, cols = np.indices(k.shape)
    ym, y0, yp = data[k - 1, rows, cols], data[k, rows, cols], data[k + 1, rows, cols]
    den = ym - 2.0 * y0 + yp
    num = ym - yp
    with np.errstate(invalid="ignore", divide="ignore"):
        x = np.clip(0.5 * num / den, -1.0, 1.0)
        peak = y0 - 0.25 * num * x
        var_x = rms**2 * ((den - num) ** 2 + (den + num) ** 2 + 4.0 * num**2) / (4.0 * den**4)
    v0 = vel[k] + x * dv
    dv0 = np.sqrt(var_x) * dv
    good = (peak > clip * rms) & (den < 0) & np.isfinite(v0) & (dv0 < 0.5 * dv)
    v0 = np.where(good, v0, np.nan)
    dv0 = np.where(good, dv0, np.nan)
    return LineMap(v0, dv0, np.where(good, peak, np.nan), cube.x, cube.y, rms)


# ---------------------------------------------------------------------------
# The flared Keplerian disc
# ---------------------------------------------------------------------------
PARAMS = ("x0", "y0", "inc", "pa", "mstar", "vsys", "z0", "phi", "z1", "psi")


@dataclass
class DiscModel:
    x0: float  # arcsec
    y0: float
    inc: float  # deg, 0-90 (the near side is chosen by `near`)
    pa: float  # deg, red-shifted major axis
    mstar: float  # M_sun, for distance d
    vsys: float  # km/s
    z0: float  # arcsec at r = 1"
    phi: float
    z1: float = 0.0
    psi: float = 4.0
    distance: float = 101.0
    near: int = 1  # +1: near side of the minor axis at PA - 90 deg; -1: at PA + 90 deg

    @classmethod
    def from_vector(cls, p, distance: float, near: int = 1) -> "DiscModel":
        return cls(*[float(v) for v in p], distance=distance, near=near)

    @property
    def near_side_pa(self) -> float:
        """Position angle (deg) of the near side of the minor axis."""
        return (self.pa - 90.0 * self.near) % 360.0

    @property
    def rotation_on_sky(self) -> str:
        """Sense of rotation on an image with north up and east left. For near = +1 the frame
        (major axis, minor axis, towards us) is right-handed and the gas circulates
        counter-clockwise about the normal facing us, i.e. counter-clockwise as we see it;
        mirroring the minor axis (near = -1) reverses that."""
        return "counter-clockwise" if self.near == 1 else "clockwise"

    def vector(self) -> np.ndarray:
        return np.array([getattr(self, k) for k in PARAMS])

    def height(self, r):
        r = np.maximum(np.asarray(r, dtype=float), 1e-6)
        return np.clip(self.z0 * r**self.phi - self.z1 * r**self.psi, 0.0, None)

    def disc_coords(self, x, y, surface: bool = True, n_iter: int = 10):
        """Sky offsets -> (r, phi, z) on the front emission surface (or the midplane)."""
        pa, inc = np.radians(self.pa), np.radians(self.inc)
        dx, dy = np.asarray(x) - self.x0, np.asarray(y) - self.y0
        x_maj = dx * np.sin(pa) + dy * np.cos(pa)
        y_min = self.near * (dx * np.cos(pa) - dy * np.sin(pa))  # mirrored when the near side is at PA + 90
        r = np.hypot(x_maj, y_min / np.cos(inc))
        z = np.zeros_like(r)
        if surface:
            for _ in range(n_iter):
                z = self.height(r)
                r = np.hypot(x_maj, (y_min - z * np.sin(inc)) / np.cos(inc))
            z = self.height(r)
        yy = (y_min - z * np.sin(inc)) / np.cos(inc)
        return r, np.arctan2(yy, x_maj), z

    def sky(self, r, phi, z=0.0):
        pa, inc = np.radians(self.pa), np.radians(self.inc)
        x_maj = r * np.cos(phi)
        y_min = self.near * (r * np.sin(phi) * np.cos(inc) + z * np.sin(inc))
        return self.x0 + x_maj * np.sin(pa) + y_min * np.cos(pa), self.y0 + x_maj * np.cos(pa) - y_min * np.sin(pa)

    def v_kep(self, r, z=0.0):
        """Keplerian speed (km/s) at cylindrical radius r and height z (arcsec)."""
        rc = np.asarray(r) * self.distance * AU_CGS
        zc = np.asarray(z) * self.distance * AU_CGS
        return np.sqrt(G_CGS * self.mstar * M_SUN_CGS * rc**2 / (rc**2 + zc**2) ** 1.5) / 1e5

    def velocity(self, x, y, surface: bool = True):
        r, phi, z = self.disc_coords(x, y, surface)
        return self.vsys + self.v_kep(r, z) * np.cos(phi) * np.sin(np.radians(self.inc)), r


@dataclass
class KinematicFit:
    model: DiscModel
    errors: dict
    chi2_red: float
    jitter: float
    n_points: int
    r_range: tuple[float, float]
    samples: np.ndarray | None = None
    extra: dict = field(default_factory=dict)


def _sample_points(lm: LineMap, guess: DiscModel, r_range: tuple[float, float], spacing: float):
    """Pixels one beam apart (so residuals are nearly independent) within r_range."""
    step = max(1, int(round(spacing / abs(lm.x[1] - lm.x[0]))))
    xx, yy = np.meshgrid(lm.x, lm.y)
    sub = (slice(step // 2, None, step), slice(step // 2, None, step))
    xs, ys, v, dv = xx[sub].ravel(), yy[sub].ravel(), lm.v0[sub].ravel(), lm.dv0[sub].ravel()
    r, _, _ = guess.disc_coords(xs, ys, surface=False)
    ok = np.isfinite(v) & np.isfinite(dv) & (r > r_range[0]) & (r < r_range[1])
    return xs[ok], ys[ok], v[ok], dv[ok]


def fit_rotation(lm: LineMap, guess: DiscModel, r_range: tuple[float, float], spacing: float,
                 rng: np.random.Generator, fix_taper: bool = False, mcmc: bool = True,
                 inc_prior: tuple[float, float] | None = None, centre_box: float = 0.05,
                 nwalkers: int = 48, nsteps: int = 3000, burn: int = 1500, clip_sigma: float = 5.0) -> KinematicFit:
    """Least squares, a 5-sigma clip of outliers (cloud absorption, the far side), then an
    MCMC with a jitter term added in quadrature to each pixel's error.

    Velocities constrain roughly M sin^2 i, so the inclination trades off against the mass;
    ``inc_prior`` = (mean, sigma) adds a Gaussian prior on it (Teague et al. 2019 used the
    continuum's 46.7 +/- 0.1 deg). The surface's two power laws are bounded so that they
    cannot cancel each other (0 < z0 < 0.6", 0.5 < phi < 2, 0 < z1 < 0.5", 2 < psi < 8), and the
    centre is kept within ``centre_box`` arcsec of the starting guess (the continuum centre):
    the star sits there to within a few mas, while the gas velocities, with their meridional
    flows and spirals, can otherwise pull a fitted centre far off. The least-squares stage
    uses a soft-L1 loss so that a few hundred structured pixels cannot steer it."""
    xs, ys, v, dv = _sample_points(lm, guess, r_range, spacing)
    names = list(PARAMS[:8]) if fix_taper else list(PARAMS)
    base = guess.vector()

    def vec_to_model(p):
        full = base.copy()
        full[: len(names)] = p
        if fix_taper:
            full[8] = 0.0
        return DiscModel.from_vector(full, guess.distance, guess.near)

    def resid(p, xs, ys, v, sig):
        mv, _ = vec_to_model(p).velocity(xs, ys)
        out = (v - mv) / sig
        if inc_prior is not None:
            out = np.append(out, (p[2] - inc_prior[0]) / inc_prior[1])
        return out

    p0 = base[: len(names)].copy()
    lo = np.array([guess.x0 - centre_box, guess.y0 - centre_box, 1.0, -720.0, 0.05, -50.0, 0.0, 0.5, 0.0, 2.0][: len(names)])
    hi = np.array([guess.x0 + centre_box, guess.y0 + centre_box, 89.0, 720.0, 10.0, 50.0, 0.6, 2.0, 0.5, 8.0][: len(names)])
    p0 = np.clip(p0, lo + 1e-6, hi - 1e-6)
    sig = np.hypot(dv, 0.05)
    keep = np.ones(v.size, bool)
    for _ in range(3):
        res = optimize.least_squares(resid, p0, args=(xs[keep], ys[keep], v[keep], sig[keep]), bounds=(lo, hi),
                                     x_scale="jac", loss="soft_l1", f_scale=2.0)
        p0 = res.x
        r_all = resid(p0, xs, ys, v, sig)[: v.size]
        scale = 1.4826 * np.median(np.abs(r_all[keep]))
        keep = np.abs(r_all) < clip_sigma * scale
    r1 = resid(p0, xs[keep], ys[keep], v[keep], 1.0)[: int(keep.sum())]
    jitter = float(np.sqrt(max(np.median(r1**2) - np.median(dv[keep] ** 2), 1e-6)))
    xs, ys, v, dv = xs[keep], ys[keep], v[keep], dv[keep]
    best = p0
    samples = None
    if mcmc:
        ndim = len(names) + 1

        def log_prob(theta):
            p, lnj = theta[:-1], theta[-1]
            if np.any(p < lo) or np.any(p > hi) or not (-8 < lnj < 1):
                return -np.inf
            s2 = dv**2 + np.exp(2 * lnj)
            mv, _ = vec_to_model(p).velocity(xs, ys)
            lp = -0.5 * np.sum((v - mv) ** 2 / s2 + np.log(s2))
            if inc_prior is not None:
                lp += -0.5 * ((p[2] - inc_prior[0]) / inc_prior[1]) ** 2
            return lp

        start = np.concatenate([best, [np.log(jitter)]])
        scale = np.concatenate([np.maximum(np.abs(best) * 1e-4, 1e-4), [0.01]])
        walkers = start + scale * rng.standard_normal((nwalkers, ndim))
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob)
        sampler.random_state = np.random.RandomState(int(rng.integers(2**31)))  # emcee keeps its own RNG
        sampler.run_mcmc(walkers, nsteps, progress=False)
        samples = sampler.get_chain(discard=burn, flat=True)
        best = np.median(samples[:, :-1], axis=0)
        jitter = float(np.exp(np.median(samples[:, -1])))
        errors = {n: float(np.std(samples[:, i])) for i, n in enumerate(names)}
        errors["jitter"] = float(np.std(np.exp(samples[:, -1])))
    else:
        j = res.jac
        try:
            cov = np.linalg.pinv(j.T @ j) * np.mean(res.fun**2)
        except np.linalg.LinAlgError:
            cov = np.full((len(names), len(names)), np.nan)
        errors = {n: float(np.sqrt(cov[i, i])) for i, n in enumerate(names)}
    model = vec_to_model(best)
    mv, _ = model.velocity(xs, ys)
    chi2 = float(np.sum((v - mv) ** 2 / (dv**2 + jitter**2)) / max(1, v.size - len(names) - 1))
    return KinematicFit(model, errors, chi2, jitter, int(v.size), r_range, samples,
                        dict(n_clipped=int((~keep).sum()), names=names))


# ---------------------------------------------------------------------------
# Rotation curve: v_phi and v_r in annuli
# ---------------------------------------------------------------------------
def rotation_curve(lm: LineMap, model: DiscModel, edges: np.ndarray, beam: float, rng: np.random.Generator,
                   n_boot: int = 200) -> dict:
    """In each annulus of the emission surface fit v = A + B cos(phi) + C sin(phi):
    v_phi = B / sin i, v_r = C / sin i. Errors by bootstrapping beam-sized blocks."""
    xx, yy = np.meshgrid(lm.x, lm.y)
    r, phi, z = model.disc_coords(xx, yy)
    ok = np.isfinite(lm.v0)
    sin_i = np.sin(np.radians(model.inc))
    out = {k: [] for k in ("r", "z", "v_phi", "v_phi_err", "v_r", "v_r_err", "v_sys", "v_kep", "n_beams")}
    block = max(1, int(round(beam / abs(lm.x[1] - lm.x[0]))))
    bid = (np.arange(xx.shape[0])[:, None] // block) * 10000 + (np.arange(xx.shape[1])[None, :] // block)
    for a, b in zip(edges[:-1], edges[1:]):
        sel = ok & (r >= a) & (r < b)
        if sel.sum() < 30:
            continue
        ph, v, w = phi[sel], lm.v0[sel], bid[sel]
        design = np.column_stack([np.ones_like(ph), np.cos(ph), np.sin(ph)])
        coef, *_ = np.linalg.lstsq(design, v, rcond=None)
        blocks = np.unique(w)
        boot = []
        for _ in range(n_boot):
            pick = rng.choice(blocks, blocks.size)
            idx = np.concatenate([np.where(w == q)[0] for q in pick[: min(pick.size, 400)]])
            c, *_ = np.linalg.lstsq(design[idx], v[idx], rcond=None)
            boot.append(c)
        boot = np.array(boot)
        rm = 0.5 * (a + b)
        zm = float(np.median(z[sel]))
        out["r"].append(rm)
        out["z"].append(zm)
        out["v_phi"].append(coef[1] / sin_i)
        out["v_phi_err"].append(np.std(boot[:, 1]) / sin_i)
        out["v_r"].append(coef[2] / sin_i)
        out["v_r_err"].append(np.std(boot[:, 2]) / sin_i)
        out["v_sys"].append(coef[0])
        out["v_kep"].append(float(model.v_kep(rm, zm)))
        out["n_beams"].append(int(blocks.size))
    return {k: np.array(v) for k, v in out.items()}
