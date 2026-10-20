"""The Radcliffe Wave: a 2.7-kpc chain of dust clouds that undulates through the Galactic plane.

Two published models of its shape, both a curved baseline in the plane with
a damped vertical sine wave along it:

* **Alves et al. (2020, Nature 578, 237)**: baseline quadratic in x for y and
  z through three 3D anchor points; d = straight-line distance from the first
  anchor; dz = A exp[-delta (d/kpc)^2] sin[(2 pi d / P)(1 + (d/d_max)/gamma) + phi].
  Their Extended Data Table 2 values (delta = 6.68) put the largest
  displacement at only ~40 pc, not the 160 pc they quote, so the printed
  equation and table cannot both be right; we fit the printed form anyway.
* **Konietzka et al. (2024, Nature 628, 62)**: baseline y(x) quadratic through
  three anchors in the plane; s = arc length along it from the first anchor;
  z(s) = zeta(s) sin(2 pi Lambda(s) + phi) with Lambda(s) = s / (p - gamma s)
  (a wavelength that shrinks along the wave, lambda = (p - gamma s)^2 / p) and
  a Lorentzian envelope zeta(s) = -A / (1 + ((s - s0)/delta)^2).

The data are cloud positions: the 326 sightlines of Zucker et al. (2020) with
their published distances (what both papers fitted), and the same sightlines
with distances measured in the Edenhofer et al. (2024) map. The likelihood is
a mixture: a cloud either belongs to a tube of Gaussian radius sigma around
the curve (density spread evenly along its length) or to a uniform background
filling the clouds' bounding box. Fits use dynesty (as both papers did).
"""

from __future__ import annotations

import json
import multiprocessing as mp

import healpy as hp
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

import catalogs
import config as C
from dustmap import DustSphere, lbd_to_xyz

N_CURVE = 900
FIT_CURVE = 300        # points along the curve inside the likelihood (spacing ~8 pc, vs a ~50 pc tube)

# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------

KON_NAMES = ["x0", "y0", "x1", "y1", "x2", "y2", "A", "p", "phi", "s0", "delta", "gamma"]
#: Konietzka et al. (2024) Extended Data Table 2: column 3 (molecular clouds only)
#: and column 6 (clouds + clusters).
KON_PUBLISHED = {
    "clouds": dict(x0=-853.13, y0=-797.19, x1=-276.15, y1=18.27, x2=289.05, y2=1386.24, p=4773.99,
                   A=219.88, delta=797.78, s0=482.71, gamma=1.49, phi=-0.12),
    "clouds+clusters": dict(x0=-852.99, y0=-807.94, x1=-275.13, y1=30.48, x2=293.84, y2=1387.37,
                            p=4779.9, A=218.64, delta=740.01, s0=544.99, gamma=1.46, phi=-0.15),
}
ALV_NAMES = ["x0", "y0", "z0", "x1", "y1", "z1", "x2", "y2", "z2", "P", "A", "phi", "gamma", "delta"]
#: Alves et al. (2020) Extended Data Table 2 (medians).
ALV_PUBLISHED = dict(x0=-910, y0=-860, z0=-30, x1=-270, y1=-20, z1=-10, x2=290, y2=1400, z2=30,
                     P=3560, A=160, phi=2.89, gamma=0.50, delta=6.68)


def _quad(xa, ya, x):
    """The quadratic through three points (Lagrange form) and its slope, evaluated at x."""
    x0, x1, x2 = xa
    y0, y1, y2 = ya
    d0, d1, d2 = (x0 - x1) * (x0 - x2), (x1 - x0) * (x1 - x2), (x2 - x0) * (x2 - x1)
    y = y0 * (x - x1) * (x - x2) / d0 + y1 * (x - x0) * (x - x2) / d1 + y2 * (x - x0) * (x - x1) / d2
    dy = y0 * (2 * x - x1 - x2) / d0 + y1 * (2 * x - x0 - x2) / d1 + y2 * (2 * x - x0 - x1) / d2
    return y, dy


def _arclength(x, slope):
    f = np.sqrt(1 + slope**2)
    s = np.empty_like(x)
    s[0] = 0.0
    np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(x), out=s[1:])
    return s


def kon_curve(th: dict, n: int = N_CURVE):
    """Konietzka model -> (x, y, z, s) along the baseline (s = in-plane arc length)."""
    x = np.linspace(th["x0"], th["x2"], n)
    y, slope = _quad((th["x0"], th["x1"], th["x2"]), (th["y0"], th["y1"], th["y2"]), x)
    s = _arclength(x, slope)
    lam = s / (th["p"] - th["gamma"] * s)
    zeta = -th["A"] / (1 + ((s - th["s0"]) / th["delta"]) ** 2)
    z = zeta * np.sin(2 * np.pi * lam + th["phi"])
    return x, y, z, s


def kon_wavelength(th: dict, s):
    return (th["p"] - th["gamma"] * np.asarray(s)) ** 2 / th["p"]


def alves_curve(th: dict, n: int = N_CURVE):
    """Alves model as printed -> (x, y, z, d) with d the distance from the first anchor."""
    xa = (th["x0"], th["x1"], th["x2"])
    x = np.linspace(th["x0"], th["x2"], n)
    y, _ = _quad(xa, (th["y0"], th["y1"], th["y2"]), x)
    zc, _ = _quad(xa, (th["z0"], th["z1"], th["z2"]), x)
    d = np.sqrt((x - th["x0"]) ** 2 + (y - th["y0"]) ** 2 + (zc - th["z0"]) ** 2)
    dz = th["A"] * np.exp(-th["delta"] * (d / 1000) ** 2) * np.sin(
        (2 * np.pi * d / th["P"]) * (1 + (d / d[-1]) / th["gamma"]) + th["phi"])
    return x, y, zc + dz, d


CURVES = {"konietzka": (kon_curve, KON_NAMES), "alves": (alves_curve, ALV_NAMES)}


def describe(model: str, th: dict, pts: np.ndarray | None = None, member: np.ndarray | None = None) -> dict:
    """Amplitude, wavelengths and length of a model curve.

    Length: 3D arc length between the projections of the first and last member
    clouds onto the curve (Alves et al.'s definition), or of the whole curve.
    """
    curve, _ = CURVES[model]
    x, y, z, s = curve(th, 4000)
    xyz = np.c_[x, y, z]
    seg = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(xyz, axis=0), axis=1))]
    i0, i1 = 0, len(x) - 1
    if pts is not None:
        sel = pts if member is None else pts[member]
        idx = np.argmin(((sel[:, None, :] - xyz[None]) ** 2).sum(-1), axis=1)
        i0, i1 = idx.min(), idx.max()
    out = {"amplitude_pc": float(np.max(np.abs(z[i0:i1 + 1]))), "length_pc": float(seg[i1] - seg[i0]),
           "trough_xyz": [float(v) for v in xyz[i0 + np.argmin(z[i0:i1 + 1])]],
           "crest_xyz": [float(v) for v in xyz[i0 + np.argmax(z[i0:i1 + 1])]],
           "in_plane_length_pc": float(s[i1] - s[i0]) if model == "konietzka" else None}
    if model == "konietzka":
        lam = kon_wavelength(th, s[i0:i1 + 1])
        out.update({"wavelength_start_pc": float(lam[0]), "wavelength_end_pc": float(lam[-1]),
                    "wavelength_mean_pc": float(np.mean(lam))})
    return out


# ---------------------------------------------------------------------------
# Likelihood and priors
# ---------------------------------------------------------------------------

try:
    from numba import njit

    @njit(cache=True, fastmath=True)
    def _min_d2(pts, curve):
        """Squared distance from each point to the nearest curve vertex (compiled loop)."""
        out = np.empty(pts.shape[0])
        for i in range(pts.shape[0]):
            best = 1e300
            for j in range(curve.shape[0]):
                dx = pts[i, 0] - curve[j, 0]
                dy = pts[i, 1] - curve[j, 1]
                dz = pts[i, 2] - curve[j, 2]
                d = dx * dx + dy * dy + dz * dz
                if d < best:
                    best = d
            out[i] = best
        return out
except ImportError:                                   # same result, ~10x slower
    def _min_d2(pts, curve):
        return ((pts[:, None, :] - curve[None]) ** 2).sum(-1).min(1)


#: Background clouds are spread evenly over this cross-section (pc^2) perpendicular to the wave:
#: 2 kpc across the plane x 0.7 kpc in height, the slab the cloud catalogue samples. A modelling
#: choice, like Alves et al.'s constant 1e-7 background term.
BACKGROUND_AREA = 2000.0 * 700.0


class TubeMixture:
    """Each cloud is on the wave (Gaussian scatter sigma across the curve) or in the background.

    Only a cloud's offset from the curve is modelled, not where along the curve it sits: the cloud
    catalogue samples the wave very unevenly (24 sightlines in Orion, 4 in CMa OB1), so spreading
    the wave's clouds evenly along its length would reward short curves through crowded regions.
    On the wave: a 2D Gaussian across the tube; background: uniform over BACKGROUND_AREA.
    """

    def __init__(self, model: str, pts: np.ndarray, priors: dict):
        self.model = model
        self.curve, self.names = CURVES[model]
        self.pts = np.asarray(pts, dtype=float)
        self.log_aout = float(np.log(BACKGROUND_AREA))
        self.priors = priors                       # name -> ("n", mu, sd) | ("u", lo, hi) | ("ln", mu, sd) | ("lu", lo, hi)
        self.all_names = self.names + ["sigma", "f"]
        self.ndim = len(self.all_names)

    def prior_transform(self, u):
        """Unit cube -> parameters; vectorised (scipy.stats' ppf per parameter cost ~3 ms a call)."""
        if not hasattr(self, "_pt"):
            from scipy.special import ndtri
            kinds = [self.priors[n][0] for n in self.all_names]
            a = np.array([self.priors[n][1] for n in self.all_names], dtype=float)
            b = np.array([self.priors[n][2] for n in self.all_names], dtype=float)
            masks = {k: np.array([kk == k for kk in kinds]) for k in ("n", "u", "ln", "lu")}
            self._pt = (ndtri, a, b, masks)
        ndtri, a, b, mk = self._pt
        th = np.empty_like(u)
        g = a + b * ndtri(np.clip(u, 1e-15, 1 - 1e-15))          # normal and log-normal
        lin = a + (b - a) * u                                     # uniform and log-uniform
        th[mk["n"]] = g[mk["n"]]
        th[mk["ln"]] = np.exp(g[mk["ln"]])
        th[mk["u"]] = lin[mk["u"]]
        th[mk["lu"]] = np.exp(lin[mk["lu"]])
        return th

    def unpack(self, theta) -> dict:
        return dict(zip(self.all_names, theta))

    def member_terms(self, theta):
        th = self.unpack(theta)
        if not (th["x0"] < th["x1"] < th["x2"]):
            return None
        x, y, z, s = self.curve(th, FIT_CURVE)
        if not np.all(np.isfinite(z)):
            return None
        if self.model == "konietzka" and th["p"] - th["gamma"] * s[-1] <= 50.0:
            return None                              # the period would reach zero before the wave ends
        xyz = np.empty((x.size, 3))
        xyz[:, 0], xyz[:, 1], xyz[:, 2] = x, y, z
        d2 = _min_d2(self.pts, xyz)
        sig, f = th["sigma"], th["f"]
        log_wave = np.log1p(-f) - d2 / (2 * sig**2) - np.log(2 * np.pi * sig**2)
        log_out = np.log(f) - self.log_aout
        return log_wave, log_out

    def loglike(self, theta) -> float:
        terms = self.member_terms(theta)
        if terms is None:
            return -1e300
        log_wave, log_out = terms
        return float(np.logaddexp(log_wave, log_out).sum())


def priors_konietzka(kind: str = "published") -> dict:
    """Konietzka et al.'s Extended Data Table 1, or broad priors that let the data speak."""
    if kind == "published":
        p = {"x0": ("n", -850, 100), "y0": ("n", -790, 100), "x1": ("n", -255, 100), "y1": ("n", 23, 100),
             "x2": ("n", 313, 100), "y2": ("n", 1368, 100), "p": ("n", 4650, 300), "A": ("n", 215, 20),
             "delta": ("ln", 6.5, 0.1), "s0": ("n", 680, 200), "gamma": ("ln", 0.35, 0.07),
             "phi": ("n", -0.09, 0.5)}
    else:
        p = {"x0": ("n", -850, 200), "y0": ("n", -790, 200), "x1": ("n", -255, 200), "y1": ("n", 23, 200),
             "x2": ("n", 313, 200), "y2": ("n", 1368, 200), "p": ("u", 2000, 9000), "A": ("u", 20, 500),
             "delta": ("lu", np.log(100), np.log(5000)), "s0": ("u", -500, 2500),
             "gamma": ("u", 0.2, 2.2), "phi": ("u", -np.pi, np.pi)}
    p.update({"sigma": ("lu", 3.0, 5.5), "f": ("u", 0.05, 0.95)})
    return p


def priors_alves(kind: str = "published") -> dict:
    """Alves et al.'s Extended Data Table 1 (their f prior U(0.15, 0.4) widened: our
    background model differs), or broad priors on the damping and period terms."""
    p = {"x0": ("n", -900, 100), "y0": ("n", -900, 100), "z0": ("n", 0, 50), "x1": ("n", -300, 100),
         "y1": ("n", 0, 100), "z1": ("n", 0, 50), "x2": ("n", 300, 100), "y2": ("n", 1400, 100),
         "z2": ("n", 0, 50), "P": ("n", 3500, 300), "A": ("n", 170, 20), "phi": ("n", 2.9, 0.5),
         "gamma": ("ln", -0.5, 0.5), "delta": ("ln", -0.5, 0.5)}
    if kind == "broad":
        p.update({"P": ("u", 1000, 8000), "A": ("u", 20, 500), "phi": ("u", -np.pi, np.pi),
                  "gamma": ("lu", np.log(0.05), np.log(50)), "delta": ("lu", np.log(0.01), np.log(30))})
    p.update({"sigma": ("lu", 3.0, 5.5), "f": ("u", 0.05, 0.95)})
    return p


# ---------------------------------------------------------------------------
# Fitting
# ---------------------------------------------------------------------------

_MODEL: TubeMixture | None = None


def _ll(theta):
    return _MODEL.loglike(theta)


def _pt(u):
    return _MODEL.prior_transform(u)


def fit(model: str, pts: np.ndarray, priors: dict, tag: str, *, nlive: int = 500,
        workers: int = 8, force: bool = False, seed: int = 20261019, maxiter: int | None = None) -> dict:
    """Nested-sampling fit, cached as JSON (posterior summary) + npz (equal-weight samples)."""
    path = C.CACHE_DIR / f"wave_fit_{tag}.json"
    if path.exists() and not force:
        return json.loads(path.read_text())
    import dynesty
    from dynesty.utils import resample_equal
    global _MODEL
    _MODEL = TubeMixture(model, pts, priors)
    rng = np.random.default_rng(seed)
    ctx = mp.get_context("fork")
    with ctx.Pool(workers, initializer=_one_thread) as pool:
        sampler = dynesty.NestedSampler(_ll, _pt, _MODEL.ndim, nlive=nlive, sample="rslice",
                                        pool=pool, queue_size=workers, rstate=rng)
        sampler.run_nested(dlogz=0.5, print_progress=False, maxiter=maxiter)
    res = sampler.results
    np.savez(C.CACHE_DIR / f"wave_dynesty_{tag}.npz", samples=res.samples, logwt=res.logwt, logl=res.logl,
             logz=res.logz, logzerr=res.logzerr, ncall=res.ncall)
    samples = resample_equal(res.samples, res.importance_weights(), rstate=rng)
    np.save(C.CACHE_DIR / f"wave_samples_{tag}.npy", samples)
    names = _MODEL.all_names
    summ = {n: [float(v) for v in np.percentile(samples[:, i], [16, 50, 84])] for i, n in enumerate(names)}
    best = samples[np.argmax([_MODEL.loglike(t) for t in samples[:400]])]
    # membership: posterior-mean odds of tube vs background (Alves et al.'s <R_i>)
    sub = samples[rng.choice(len(samples), min(300, len(samples)), replace=False)]
    odds = np.zeros(len(pts))
    for t in sub:
        lw, lo = _MODEL.member_terms(t)
        odds += np.exp(np.clip(lw - lo, -50, 50))
    odds /= len(sub)
    med = {n: summ[n][1] for n in names}
    derived = [describe(model, _MODEL.unpack(t), pts, odds > 1) for t in sub[:150]]
    der = {k: [float(v) for v in np.percentile([d[k] for d in derived], [16, 50, 84])]
           for k in derived[0] if derived[0][k] is not None and not isinstance(derived[0][k], list)}
    out = {"model": model, "tag": tag, "n_clouds": int(len(pts)), "logz": float(res.logz[-1]),
           "logz_err": float(res.logzerr[-1]), "ncall": int(np.sum(res.ncall)),
           "posterior": summ, "median": med, "best": dict(zip(names, map(float, best))),
           "derived": der, "members": (odds > 1).tolist(), "odds": odds.tolist()}
    path.write_text(json.dumps(out, indent=1))
    np.save(C.CACHE_DIR / f"wave_samples_{tag}.npy", samples)
    _MODEL = None
    return out


def _one_thread():
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except ImportError:
        pass


def cloud_sets() -> dict:
    """Point sets the models are fitted to (pc, heliocentric Galactic Cartesian)."""
    cd = cloud_distances()
    z = cd[cd["set"] == "zucker2020"]
    mapped = z[z["d_map"].notna() & (z["d50"] < 1950)]
    hr = catalogs.hunt_reffert()
    young = hr[hr["reliable"] & (hr["age_myr"] < 30) & (hr["dist50"] < 2000)]
    return {"zucker": z[["x", "y", "z"]].to_numpy(float),
            "map": mapped[["xm", "ym", "zm"]].to_numpy(float),
            "map_index": mapped.index.to_numpy(),
            "clusters": young[["x", "y", "z"]].to_numpy(float)}


#: The fits we run: tag -> (model, priors, point set).
FITS = {
    "kon_pub_zucker": ("konietzka", "published", "zucker"),
    "kon_pub_map": ("konietzka", "published", "map"),
    "alv_pub_zucker": ("alves", "published", "zucker"),
    "alv_pub_map": ("alves", "published", "map"),
    "kon_broad_zucker": ("konietzka", "broad", "zucker"),
    "kon_broad_map": ("konietzka", "broad", "map"),
    "alv_broad_zucker": ("alves", "broad", "zucker"),
}
# Not fitted: young clusters alone ("clusters" set). Within 2 kpc they are dominated by groups off the
# wave (Sco-Cen, the rest of the Local Arm); a trial fit wandered 4-5 sigma from its priors onto a
# curve through the solar neighbourhood. Konietzka et al. selected wave clusters first and used
# their velocities; we only overlay the clusters on the dust (figures.wave_curtain).
#: The wave model drawn in the figures and the viewer: Konietzka et al.'s form and priors,
#: fitted to cloud distances measured in the Edenhofer map.
MAIN_FIT = "kon_pub_map"


def run_fits(tags=None, force: bool = False) -> dict:
    sets = cloud_sets()
    out = {}
    for tag in tags or FITS:
        model, prior_kind, data = FITS[tag]
        priors = (priors_konietzka if model == "konietzka" else priors_alves)(prior_kind)
        print(f"== wave fit {tag}", flush=True)
        out[tag] = fit(model, sets[data], priors, tag, force=force)
        post = out[tag]["posterior"]
        print("   " + ", ".join(f"{k} {v[1]:.3g}" for k, v in post.items()), flush=True)
    return out


def curtain(th: dict, *, which: str = "2kpc", half_width: float = 150.0, ds: float = 5.0,
            dz: float = 5.0, zmax: float = 400.0, s_pad: float = 300.0) -> dict:
    """The dust along the wave's baseline, unrolled: column density in (s, z).

    For the Konietzka-model baseline ``th`` (y quadratic in x, arc length s),
    the map is sampled on the vertical surface that follows the baseline and on
    parallel surfaces up to ``half_width`` pc on either side, then integrated
    across: the result is extinction (ZGR23 mag) per (s, z) cell as seen
    perpendicular to the wave. No fitting is involved.
    """
    x, y, _, s = kon_curve(th, 4000)
    s_grid = np.arange(-s_pad, s[-1] + s_pad, ds)
    # extend the baseline straight beyond its ends
    xs = np.interp(s_grid, s, x, left=np.nan, right=np.nan)
    ys = np.interp(s_grid, s, y, left=np.nan, right=np.nan)
    t0 = np.array([x[1] - x[0], y[1] - y[0]]); t0 /= np.linalg.norm(t0)
    t1 = np.array([x[-1] - x[-2], y[-1] - y[-2]]); t1 /= np.linalg.norm(t1)
    lo, hi = s_grid < 0, s_grid > s[-1]
    xs[lo], ys[lo] = x[0] + s_grid[lo] * t0[0], y[0] + s_grid[lo] * t0[1]
    xs[hi], ys[hi] = x[-1] + (s_grid[hi] - s[-1]) * t1[0], y[-1] + (s_grid[hi] - s[-1]) * t1[1]
    tx, ty = np.gradient(xs), np.gradient(ys)
    norm = np.hypot(tx, ty)
    nx, ny = ty / norm, -tx / norm                     # in-plane normal
    us = np.arange(-half_width, half_width + 1e-9, ds)
    zs = np.arange(-zmax, zmax + 1e-9, dz)
    sphere = DustSphere(which)
    img = np.zeros((s_grid.size, zs.size))
    for u in us:
        X = (xs + u * nx)[:, None] * np.ones_like(zs)[None, :]
        Y = (ys + u * ny)[:, None] * np.ones_like(zs)[None, :]
        Z = np.ones_like(xs)[:, None] * zs[None, :]
        rho = sphere.density_xyz(X.ravel(), Y.ravel(), Z.ravel()).reshape(X.shape)
        img += np.nan_to_num(rho) * ds
    return {"s": s_grid, "z": zs, "column": img, "half_width": half_width, "which": which,
            "baseline_x": xs, "baseline_y": ys, "normal": np.c_[nx, ny]}


def project_to_baseline(th: dict, pts: np.ndarray, pad: float = 1500.0) -> np.ndarray:
    """(s, u, z) of points: arc length along the baseline, in-plane offset, height.

    The baseline is continued straight beyond its ends (as in ``curtain``), so
    points past the ends get s < 0 or s > the baseline's length.
    """
    x, y, _, s = kon_curve(th, 4000)
    t0 = np.array([x[1] - x[0], y[1] - y[0]]); t0 /= np.linalg.norm(t0)
    t1 = np.array([x[-1] - x[-2], y[-1] - y[-2]]); t1 /= np.linalg.norm(t1)
    ext = np.linspace(0, pad, 600)[1:]
    x = np.r_[x[0] - ext[::-1] * t0[0], x, x[-1] + ext * t1[0]]
    y = np.r_[y[0] - ext[::-1] * t0[1], y, y[-1] + ext * t1[1]]
    s = np.r_[-ext[::-1], s, s[-1] + ext]
    d2 = (pts[:, 0, None] - x[None]) ** 2 + (pts[:, 1, None] - y[None]) ** 2
    i = np.argmin(d2, axis=1)
    tx, ty = np.gradient(x), np.gradient(y)
    norm = np.hypot(tx, ty)
    u = (pts[:, 0] - x[i]) * (ty[i] / norm[i]) + (pts[:, 1] - y[i]) * (-tx[i] / norm[i])
    return np.c_[s[i], u, pts[:, 2]]


# ---------------------------------------------------------------------------
# Cloud distances measured in the dust map
# ---------------------------------------------------------------------------

def cloud_distances(force: bool = False) -> pd.DataFrame:
    """For each Zucker et al. (2020) sightline, the distance of its dominant dust cloud.

    Profiles average the pixels within 15' of the sightline; the main map is
    used to 1240 pc and the 2 kpc map beyond. The cloud is the peak holding the
    most extinction (between its two neighbouring valleys) and its distance is
    where half of that extinction has accumulated - the analogue of the single
    dust screen that Zucker et al. fitted to each sightline.
    """
    path = C.CACHE_DIR / "cloud_distances.parquet"
    if path.exists() and not force:
        return pd.read_parquet(path)
    df = pd.concat([catalogs.zucker2020().assign(set="zucker2020"),
                    catalogs.bridging().assign(set="bridging")], ignore_index=True)
    lnd = np.arange(np.log(70.0), np.log(1985.0), 0.002)   # the 2 kpc map ends at 1989 pc
    d = np.exp(lnd)
    main, wide = DustSphere("main"), DustSphere("2kpc")
    use_main = d <= 1240.0
    out = []
    for _, row in df.iterrows():
        vec = hp.ang2vec(row["l"], row["b"], lonlat=True)
        pix = hp.query_disc(256, vec, np.radians(15 / 60), nest=True, inclusive=True)
        prof = np.empty_like(d)
        prof[use_main] = np.nanmean(main.column_profiles(pix, d[use_main]), axis=1)
        prof[~use_main] = np.nanmean(wide.column_profiles(pix, d[~use_main]), axis=1)
        prof = np.nan_to_num(prof)
        sm = gaussian_filter1d(prof, 5)                 # 1% in distance
        dd = np.gradient(d)
        peaks, props = find_peaks(sm, prominence=0.05 * sm.max())
        if peaks.size == 0:
            out.append((np.nan, np.nan, np.nan, 0))
            continue
        valleys = [p + np.argmin(sm[p:q + 1]) for p, q in zip(peaks[:-1], peaks[1:])]
        bounds = np.r_[0, valleys, len(d) - 1].astype(int)
        areas = [np.sum((prof * dd)[bounds[k]:bounds[k + 1] + 1]) for k in range(peaks.size)]
        k = int(np.argmax(areas))
        seg = (prof * dd)[bounds[k]:bounds[k + 1] + 1]
        cum = np.cumsum(seg)
        dist = np.interp(0.5 * cum[-1], cum, d[bounds[k]:bounds[k + 1] + 1])
        out.append((dist, d[peaks[k]], C.ZGR23_TO_AV * areas[k], peaks.size))
    res = pd.DataFrame(out, columns=["d_map", "d_map_peak", "av_cloud", "n_peaks"])
    df = pd.concat([df.reset_index(drop=True), res], axis=1)
    df["xm"], df["ym"], df["zm"] = lbd_to_xyz(df["l"], df["b"], df["d_map"])
    df.to_parquet(path)
    return df
