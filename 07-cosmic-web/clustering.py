"""Measuring clustering: is the galaxy distribution random or not?

The eye is very good at finding patterns, including in pure noise, so "it
looks like a web" is not evidence. The standard quantitative test is the
**two-point correlation function** xi(r): the excess probability, over a
random distribution, of finding two galaxies a distance r apart. For randomly
scattered points xi = 0 at every r. For galaxies, xi(1 Mpc) is ~100 or more.

Measuring it needs three ingredients, all here:

* **A survey mask** -- which parts of the sky were actually surveyed.
* **A random catalogue** -- points scattered with no clustering at all, but
  through exactly the same mask and over exactly the same volume, so that
  anything the survey geometry does to the pair counts, it does to both.
* **Pair counts** -- galaxy-galaxy (DD), galaxy-random (DR) and random-random
  (RR) pairs, binned by separation. The Landy & Szalay (1993) estimator
  combines them: xi = (DD - 2 DR + RR) / RR, each suitably normalised.

Separations are split across (r_p) and along (pi) the line of sight, because
galaxies' own orbital motions distort only the line-of-sight direction (see
:func:`cosmo_physics.rp_pi`). Integrating over pi gives the *projected*
correlation function w_p(r_p), which those distortions cannot touch.
"""

from __future__ import annotations

import multiprocessing as mp
import os
from dataclasses import dataclass

import numpy as np
from astropy_healpix import HEALPix
from scipy.spatial import cKDTree

# ---------------------------------------------------------------------------
# The angular mask
# ---------------------------------------------------------------------------


@dataclass
class Mask:
    """Which HEALPix pixels the survey covered.

    ``surveyed`` marks pixels holding a normal density of galaxies; ``interior``
    additionally requires all eight neighbours to be surveyed, trimming the
    ragged, partially covered edge of the footprint.
    """

    nside: int
    counts: np.ndarray
    surveyed: np.ndarray
    interior: np.ndarray

    @property
    def hp(self) -> HEALPix:
        return HEALPix(nside=self.nside, order="nested")

    @property
    def pixel_area_deg2(self) -> float:
        return float(self.hp.pixel_area.to_value("deg2"))

    def area_deg2(self, which: str = "interior") -> float:
        return float(getattr(self, which).sum()) * self.pixel_area_deg2

    def solid_angle_sr(self, which: str = "interior") -> float:
        return self.area_deg2(which) * (np.pi / 180.0) ** 2

    def pixel(self, ra, dec) -> np.ndarray:
        import astropy.units as u

        return self.hp.lonlat_to_healpix(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg)

    def contains(self, ra, dec, which: str = "interior") -> np.ndarray:
        return getattr(self, which)[self.pixel(ra, dec)]


def build_mask(ra, dec, nside: int, min_fraction: float) -> Mask:
    """Derive the survey footprint from where the galaxies are.

    A pixel is "surveyed" if it holds at least ``min_fraction`` of the median
    count of occupied pixels. With ~65 galaxies expected per pixel, a genuinely
    surveyed pixel falling below 30% of that by chance -- even looking through
    a void -- essentially never happens; a pixel on the survey's edge that is
    only partly covered does.

    Deriving the mask from the data is a simplification. The SDSS team's own
    masks come from the tiling geometry and record completeness (~90-95%)
    sector by sector; ours treats every surveyed pixel as equally complete.
    """
    hp = HEALPix(nside=nside, order="nested")
    import astropy.units as u

    pix = hp.lonlat_to_healpix(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg)
    counts = np.bincount(pix, minlength=hp.npix)
    median = np.median(counts[counts > 0])
    surveyed = counts >= min_fraction * median

    neighbours = hp.neighbours(np.arange(hp.npix))  # (8, npix); -1 where none
    ok = np.where(neighbours >= 0, surveyed[np.clip(neighbours, 0, None)], True)
    interior = surveyed & ok.all(axis=0)
    return Mask(nside=nside, counts=counts, surveyed=surveyed, interior=interior)


def random_sky(n: int, mask: Mask, rng: np.random.Generator,
               which: str = "interior") -> tuple[np.ndarray, np.ndarray]:
    """``n`` points uniformly distributed over the masked part of the sphere.

    Uniform on a sphere means uniform in RA and in sin(dec) -- not in dec,
    which would pile points up at the poles. Points are drawn over the whole
    sphere and those outside the mask rejected.
    """
    keep_ra, keep_dec, have = [], [], 0
    frac = getattr(mask, which).mean()
    while have < n:
        m = int((n - have) / frac * 1.1) + 1000
        ra = rng.uniform(0.0, 360.0, m)
        dec = np.degrees(np.arcsin(rng.uniform(-1.0, 1.0, m)))
        ok = mask.contains(ra, dec, which)
        keep_ra.append(ra[ok])
        keep_dec.append(dec[ok])
        have += int(ok.sum())
    return np.concatenate(keep_ra)[:n], np.concatenate(keep_dec)[:n]


def random_volume_distances(n: int, d_min: float, d_max: float,
                            rng: np.random.Generator) -> np.ndarray:
    """Distances for points spread uniformly through a spherical shell.

    Volume grows as d^2, so uniform-in-volume means d^3 uniform between the
    limits -- far more points at large distance than at small.
    """
    u = rng.uniform(0.0, 1.0, n)
    return (d_min**3 + u * (d_max**3 - d_min**3)) ** (1.0 / 3.0)


def jackknife_regions(mask: Mask, n_regions: int, seed: int,
                      which: str = "interior") -> np.ndarray:
    """Cut the footprint into ``n_regions`` compact, near-equal-area patches.

    Returns a region label per HEALPix pixel (-1 outside the mask). Done by
    k-means on the pixels' unit vectors, which gives compact patches of similar
    size. Repeating a measurement with each patch left out in turn estimates its
    uncertainty *including* cosmic variance -- the fact that a finite survey
    only samples a finite piece of the universe, which pure counting statistics
    would ignore.
    """
    from sklearn.cluster import KMeans

    good = np.flatnonzero(getattr(mask, which))
    lon, lat = mask.hp.healpix_to_lonlat(good)
    lon, lat = lon.to_value("rad"), lat.to_value("rad")
    xyz = np.column_stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    labels = KMeans(n_clusters=n_regions, n_init=4, random_state=seed).fit_predict(xyz)
    out = np.full(mask.hp.npix, -1, dtype=np.int64)
    out[good] = labels
    return out


# ---------------------------------------------------------------------------
# Pair counting
# ---------------------------------------------------------------------------


@dataclass
class Binning:
    """Separation bins shared by every pair count."""

    rp_edges: np.ndarray  # log-spaced, h^-1 Mpc
    pi_max: float
    pi_bin: float
    map_bin: float  # linear bins of the 2D xi(r_p, pi) map
    map_max: float

    @property
    def n_rp(self) -> int:
        return len(self.rp_edges) - 1

    @property
    def n_pi(self) -> int:
        return int(round(self.pi_max / self.pi_bin))

    @property
    def n_map(self) -> int:
        return int(round(self.map_max / self.map_bin))

    @property
    def r_max(self) -> float:
        rp_max = max(self.rp_edges[-1], self.map_max)
        return float(np.hypot(rp_max, max(self.pi_max, self.map_max)))


@dataclass
class PairCounts:
    """Pair counts, with enough bookkeeping to redo them leaving out any region."""

    wp: np.ndarray        # (n_rp, n_pi) all pairs
    wp_region: np.ndarray  # (n_regions, n_rp, n_pi) pairs touching each region
    map2d: np.ndarray     # (n_map, n_map) linear r_p x pi grid
    rp_sum: np.ndarray    # (n_rp,) sum of r_p over pairs with pi < pi_max
    n_pairs_examined: int


# Arrays shared with worker processes. They are set in the parent before the
# pool forks, so children inherit them without copying or pickling.
_SHARED: dict = {}


def _count_chunk(bounds: tuple[int, int]):
    lo, hi = bounds
    s = _SHARED
    p1, d1sq, r1 = s["p1"][lo:hi], s["d1sq"][lo:hi], s["reg1"][lo:hi]
    bn: Binning = s["binning"]

    sdm = cKDTree(p1, balanced_tree=False, compact_nodes=False).sparse_distance_matrix(
        s["tree2"], bn.r_max, output_type="ndarray")
    i, j, v = sdm["i"], sdm["j"], sdm["v"]
    del sdm
    if s["auto"]:
        # Each pair turns up twice across the chunks, once from each end; keep
        # the copy whose second member comes later in the list. That counts
        # every pair exactly once and never pairs a point with itself.
        keep = j > i + lo
        i, j, v = i[keep], j[keep], v[keep]

    # Line-of-sight split from distances alone (exact): for a pair at radii
    # D1, D2 with separation s, the midpoint vector l = p1 + p2 obeys
    # s.l = D2^2 - D1^2 and |l|^2 = 2 D1^2 + 2 D2^2 - s^2.
    a, b, v2 = d1sq[i], s["d2sq"][j], v * v
    pi = np.abs(b - a) / np.sqrt(2.0 * (a + b) - v2)
    rp2 = np.clip(v2 - pi * pi, 0.0, None)
    del a, b, v2

    n_rp, n_pi, n_reg = bn.n_rp, bn.n_pi, s["n_regions"]
    edges2 = bn.rp_edges**2
    m = (pi < bn.pi_max) & (rp2 >= edges2[0]) & (rp2 < edges2[-1])
    k_rp = np.searchsorted(edges2, rp2[m], side="right") - 1
    k = k_rp * n_pi + (pi[m] / bn.pi_bin).astype(np.int64)
    wp = np.bincount(k, minlength=n_rp * n_pi)
    rp_sum = np.bincount(k_rp, weights=np.sqrt(rp2[m]), minlength=n_rp)

    # Region bookkeeping for the jackknife: a pair touches the region of each
    # of its members (once if both are in the same one).
    ra_, rb_ = r1[i[m]], s["reg2"][j[m]]
    wp_reg = (np.bincount(ra_ * (n_rp * n_pi) + k, minlength=n_reg * n_rp * n_pi)
              + np.bincount(rb_ * (n_rp * n_pi) + k, minlength=n_reg * n_rp * n_pi))
    same = ra_ == rb_
    wp_reg -= np.bincount(ra_[same] * (n_rp * n_pi) + k[same], minlength=n_reg * n_rp * n_pi)

    mm = (pi < bn.map_max) & (rp2 < bn.map_max**2)
    km = ((np.sqrt(rp2[mm]) / bn.map_bin).astype(np.int64) * bn.n_map
          + (pi[mm] / bn.map_bin).astype(np.int64))
    map2d = np.bincount(km, minlength=bn.n_map * bn.n_map)
    return wp, wp_reg, map2d, rp_sum, len(v)


def count_pairs(p1, reg1, binning: Binning, n_regions: int, p2=None, reg2=None,
                *, workers: int | None = None, target_pairs_per_task: float = 3e6,
                progress: str | None = None) -> PairCounts:
    """Count pairs in (r_p, pi) bins; auto-pairs if ``p2`` is None.

    ``p1``/``p2`` are comoving positions in h^-1 Mpc with the observer at the
    origin; ``reg1``/``reg2`` are each point's jackknife region.

    Work is split into chunks of the first catalogue and spread over worker
    processes. Each chunk's pairs are enumerated by scipy's k-d tree out to the
    largest separation any bin needs, then binned with numpy. Memory per worker
    stays around ``target_pairs_per_task`` x 100 bytes.
    """
    auto = p2 is None
    if auto:
        p2, reg2 = p1, reg1
    p1 = np.ascontiguousarray(p1, dtype=np.float64)
    p2 = np.ascontiguousarray(p2, dtype=np.float64)
    tree2 = cKDTree(p2, balanced_tree=False, compact_nodes=False)

    # How many points each task takes: pairs per point, measured exactly on a
    # sample of points, sets the chunk size that keeps memory bounded.
    sample = np.random.default_rng(0).choice(len(p1), size=min(len(p1), 300), replace=False)
    per_point = float(np.mean(tree2.query_ball_point(p1[sample], binning.r_max,
                                                     return_length=True)))
    chunk = int(np.clip(target_pairs_per_task / max(per_point, 1.0), 16, 20000))
    tasks = [(lo, min(lo + chunk, len(p1))) for lo in range(0, len(p1), chunk)]

    _SHARED.clear()
    _SHARED.update(
        p1=p1, p2=p2, d1sq=(p1**2).sum(1), d2sq=(p2**2).sum(1),
        reg1=np.asarray(reg1, dtype=np.int64), reg2=np.asarray(reg2, dtype=np.int64),
        binning=binning, n_regions=n_regions, auto=auto, tree2=tree2,
    )
    if workers is None:
        workers = max(1, min(8, (os.cpu_count() or 2) // 2))

    n_rp, n_pi, n_map = binning.n_rp, binning.n_pi, binning.n_map
    wp = np.zeros(n_rp * n_pi)
    wp_reg = np.zeros(n_regions * n_rp * n_pi)
    map2d = np.zeros(n_map * n_map)
    rp_sum = np.zeros(n_rp)
    examined = 0

    def _accumulate(result):
        nonlocal examined
        a, b, c, d, e = result
        np.add(wp, a, out=wp)
        np.add(wp_reg, b, out=wp_reg)
        np.add(map2d, c, out=map2d)
        np.add(rp_sum, d, out=rp_sum)
        examined += e

    if workers > 1 and len(tasks) > 1:
        ctx = mp.get_context("fork")
        with ctx.Pool(workers) as pool:
            for n_done, result in enumerate(pool.imap_unordered(_count_chunk, tasks), 1):
                _accumulate(result)
                if progress and n_done % max(1, len(tasks) // 10) == 0:
                    print(f"      {progress}: {n_done / len(tasks):4.0%}", flush=True)
    else:
        for result in map(_count_chunk, tasks):
            _accumulate(result)
    _SHARED.clear()

    return PairCounts(
        wp=wp.reshape(n_rp, n_pi),
        wp_region=wp_reg.reshape(n_regions, n_rp, n_pi),
        map2d=map2d.reshape(n_map, n_map),
        rp_sum=rp_sum,
        n_pairs_examined=int(examined),
    )


# ---------------------------------------------------------------------------
# Estimators
# ---------------------------------------------------------------------------


def landy_szalay(dd, dr, rr, n_d: float, n_r: float):
    """xi = (DD' - 2 DR' + RR') / RR', primes denoting counts per possible pair.

    Normalising by the number of possible pairs -- n(n-1)/2 for auto-pairs,
    n_d n_r for cross-pairs -- makes catalogues of different sizes comparable.
    Bins with no random pairs have no defined xi and are returned as NaN.
    """
    dd_n = dd / (n_d * (n_d - 1) / 2.0)
    dr_n = dr / (n_d * n_r)
    rr_n = rr / (n_r * (n_r - 1) / 2.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(rr > 0, (dd_n - 2.0 * dr_n + rr_n) / rr_n, np.nan)


def projected_correlation(xi_rp_pi, pi_bin: float):
    """w_p(r_p) = 2 * integral of xi(r_p, pi) d(pi) from 0 to pi_max."""
    return 2.0 * np.nansum(xi_rp_pi, axis=-1) * pi_bin


@dataclass
class ClusteringResult:
    rp: np.ndarray            # pair-weighted mean r_p of each bin
    wp: np.ndarray
    wp_err: np.ndarray        # jackknife, diagonal
    wp_cov: np.ndarray        # jackknife covariance
    wp_jackknife: np.ndarray  # (n_regions, n_rp)
    xi_map: np.ndarray        # (n_map, n_map) xi on a linear r_p x pi grid
    n_data: int
    n_random: int


def measure(dd: PairCounts, dr: PairCounts, rr: PairCounts,
            n_d_region: np.ndarray, n_r_region: np.ndarray,
            binning: Binning) -> ClusteringResult:
    """Landy-Szalay w_p(r_p) with jackknife errors, plus the 2D xi map.

    For the jackknife, removing region k removes every pair touching it, and
    the catalogue sizes shrink by the region's members.
    """
    n_d, n_r = float(n_d_region.sum()), float(n_r_region.sum())
    wp_full = projected_correlation(landy_szalay(dd.wp, dr.wp, rr.wp, n_d, n_r), binning.pi_bin)

    n_reg = len(n_d_region)
    wp_jk = np.empty((n_reg, binning.n_rp))
    for k in range(n_reg):
        xi_k = landy_szalay(dd.wp - dd.wp_region[k], dr.wp - dr.wp_region[k],
                            rr.wp - rr.wp_region[k],
                            n_d - n_d_region[k], n_r - n_r_region[k])
        wp_jk[k] = projected_correlation(xi_k, binning.pi_bin)
    mean = wp_jk.mean(axis=0)
    diff = wp_jk - mean
    cov = (n_reg - 1) / n_reg * diff.T @ diff

    dd_rp = dd.wp.sum(axis=1)
    rp_mean = np.where(dd_rp > 0, dd.rp_sum / np.maximum(dd_rp, 1), np.nan)
    geometric = np.sqrt(binning.rp_edges[:-1] * binning.rp_edges[1:])
    rp_mean = np.where(np.isfinite(rp_mean), rp_mean, geometric)

    xi_map = landy_szalay(dd.map2d, dr.map2d, rr.map2d, n_d, n_r)
    return ClusteringResult(rp=rp_mean, wp=wp_full, wp_err=np.sqrt(np.diag(cov)),
                            wp_cov=cov, wp_jackknife=wp_jk, xi_map=xi_map,
                            n_data=int(n_d), n_random=int(n_r))


def fit_power_law(rp, wp, cov, rp_range, n_jackknife: int, *,
                  diagonal: bool = False) -> dict:
    """Fit w_p with the projection of xi(r) = (r / r0)^-gamma.

    With the full jackknife covariance the fit accounts for the strong
    correlation between neighbouring bins (the same structures feed them all).
    The inverse of a covariance estimated from ``n_jackknife`` samples is
    biased high; the Hartlap et al. (2007) factor corrects it.
    ``diagonal=True`` ignores the correlations; Zehavi et al. (2011) quote
    both kinds of fit.
    """
    from scipy.optimize import minimize

    from cosmo_physics import wp_power_law

    sel = (rp >= rp_range[0]) & (rp <= rp_range[1]) & np.isfinite(wp)
    x, y = rp[sel], wp[sel]
    c = cov[np.ix_(sel, sel)]
    if diagonal:
        icov = np.diag(1.0 / np.diag(c))
    else:
        hartlap = (n_jackknife - int(sel.sum()) - 2) / (n_jackknife - 1)
        icov = hartlap * np.linalg.inv(c)

    def chi2(theta):
        r0, g = theta
        if r0 <= 0 or g <= 1.05:
            return 1e30
        d = y - wp_power_law(x, r0, g)
        return float(d @ icov @ d)

    best = minimize(chi2, x0=[5.0, 1.8], method="Nelder-Mead",
                    options={"xatol": 1e-5, "fatol": 1e-6, "maxiter": 4000})
    r0, g = best.x
    # Errors from the curvature of chi^2 (a numerical Hessian).
    h = np.array([1e-3 * r0, 1e-3 * g])
    hess = np.zeros((2, 2))
    for a in range(2):
        for b in range(2):
            ea, eb = np.eye(2)[a] * h[a], np.eye(2)[b] * h[b]
            hess[a, b] = (chi2(best.x + ea + eb) - chi2(best.x + ea - eb)
                          - chi2(best.x - ea + eb) + chi2(best.x - ea - eb)) / (4 * h[a] * h[b])
    try:
        pcov = 2.0 * np.linalg.inv(hess)
        errs = np.sqrt(np.clip(np.diag(pcov), 0, None))
    except np.linalg.LinAlgError:
        errs = np.array([np.nan, np.nan])
    dof = int(sel.sum()) - 2
    return {"r0": float(r0), "r0_err": float(errs[0]), "gamma": float(g),
            "gamma_err": float(errs[1]), "chi2": float(best.fun), "dof": dof,
            "rp_range": list(rp_range), "diagonal": diagonal, "n_points": int(sel.sum())}
