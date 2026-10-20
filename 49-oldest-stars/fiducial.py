"""The ridge line of a cluster's colour-magnitude diagram, in "natural" coordinates.

A cluster's main sequence, turnoff, subgiant branch and lower giant branch form
one curved sequence: steep on the main sequence, flat on the subgiant branch,
steep again on the giant branch. To trace its centre (the ridge or fiducial
line) everywhere with one recipe, we measure each star's position *along* a
guide curve (s) and its signed offset *across* it (t), in a plane where colour
is stretched by a factor K so that both axes carry similar information. In bins
of s, the median t is the ridge; its scatter over sqrt(N) is the statistical
uncertainty. Starting from a model isochrone slid onto the data, two passes
replace the guide by the data's own (smoothed) ridge, so the result does not
depend on the starting model.

The age fit (agefit.py) compares a model isochrone with this ridge along the
same normal directions.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import gaussian_filter1d

#: Colour stretch: one magnitude of colour counts as K magnitudes of brightness.
K = 5.0


@dataclass
class Ridge:
    s: np.ndarray  # arc length along the guide (scaled units)
    c: np.ndarray  # colour of each ridge point (F606W - F814W)
    m: np.ndarray  # magnitude (F606W)
    sig: np.ndarray  # statistical uncertainty across the sequence (scaled units)
    nx: np.ndarray  # unit normal of the guide at the point (scaled plane)
    ny: np.ndarray
    n: np.ndarray  # stars in the bin
    i_to: int  # index of the turnoff (bluest ridge point)
    width: np.ndarray  # robust spread of stars across the sequence (scaled)

    @property
    def c_to(self) -> float:
        return float(self.c[self.i_to])

    @property
    def m_to(self) -> float:
        return float(self.m[self.i_to])

    def subset(self, keep: np.ndarray) -> "Ridge":
        idx = np.flatnonzero(keep)
        i_to = int(np.argmin(abs(idx - self.i_to)))
        return Ridge(self.s[idx], self.c[idx], self.m[idx], self.sig[idx], self.nx[idx], self.ny[idx],
                     self.n[idx], i_to, self.width[idx])


def _polyline_frame(gx: np.ndarray, gy: np.ndarray):
    """Arc length, unit tangents and unit normals of a polyline."""
    dx, dy = np.gradient(gx), np.gradient(gy)
    norm = np.hypot(dx, dy)
    tx, ty = dx / norm, dy / norm
    s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(gx), np.diff(gy)))])
    return s, tx, ty, -ty, tx  # normal = tangent rotated by +90 degrees


def project(gx, gy, px, py):
    """(s, t) of points (px, py) relative to the polyline (gx, gy), nearest segment."""
    s_nodes = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(gx), np.diff(gy)))])
    ax, ay = gx[:-1], gy[:-1]
    bx, by = np.diff(gx), np.diff(gy)
    seg_len2 = bx**2 + by**2
    s_out = np.zeros(len(px))
    t_out = np.zeros(len(px))
    for start in range(0, len(px), 20000):  # chunks bound the memory
        sl = slice(start, start + 20000)
        qx = px[sl, None] - ax[None, :]
        qy = py[sl, None] - ay[None, :]
        u = np.clip((qx * bx + qy * by) / seg_len2, 0.0, 1.0)
        rx, ry = qx - u * bx, qy - u * by
        d2 = rx**2 + ry**2
        j = np.argmin(d2, axis=1)
        rows = np.arange(d2.shape[0])
        uu = u[rows, j]
        seg = np.sqrt(seg_len2[j])
        # signed offset: positive on the left of the direction of travel
        cross = (bx[j] * ry[rows, j] - by[j] * rx[rows, j]) / seg
        s_out[sl] = s_nodes[j] + uu * seg
        t_out[sl] = cross
    return s_out, t_out


def slide_guide(m, c, iso_m, iso_c, *, dm_range=(-0.6, 0.6), dc_range=(-0.08, 0.08)):
    """Shift a model isochrone in colour and magnitude onto the densest part of the data."""
    best, best_n = (0.0, 0.0), -1
    gx0, gy0 = K * iso_c, iso_m
    sample = np.random.default_rng(0).choice(len(m), size=min(len(m), 6000), replace=False)
    px, py = K * c[sample], m[sample]
    for dm in np.arange(dm_range[0], dm_range[1] + 1e-9, 0.05):
        for dc in np.arange(dc_range[0], dc_range[1] + 1e-9, 0.02):
            _, t = project(gx0 + K * dc, gy0 + dm, px, py)
            n = int(np.sum(abs(t) < 0.06))
            if n > best_n:
                best_n, best = n, (dm, dc)
    # refine on a finer grid around the best coarse shift
    dm0, dc0 = best
    for dm in np.arange(dm0 - 0.05, dm0 + 0.051, 0.01):
        for dc in np.arange(dc0 - 0.02, dc0 + 0.0201, 0.004):
            _, t = project(gx0 + K * dc, gy0 + dm, px, py)
            n = int(np.sum(abs(t) < 0.04))
            if n > best_n:
                best_n, best = n, (dm, dc)
    return iso_m + best[0], iso_c + best[1]


def ridge_line(m, c, guide_m, guide_c, *, ds: float = 0.06, tmax: float = 0.15,
               min_stars: int = 25, passes: int = 3) -> Ridge:
    """Ridge line of the stars (m, c), refined from an initial guide curve."""
    gx, gy = K * np.asarray(guide_c, float), np.asarray(guide_m, float)
    px, py = K * np.asarray(c, float), np.asarray(m, float)
    for p in range(passes):
        s_star, t_star = project(gx, gy, px, py)
        s_nodes, tx, ty, nx, ny = _polyline_frame(gx, gy)
        # stars beyond either end of the guide project onto its end point: drop them
        inside = (s_star > s_nodes[0] + 1e-6) & (s_star < s_nodes[-1] - 1e-6)
        near = inside & (abs(t_star) < tmax)
        rows = []
        # adaptive bins: at least ds long and at least min_stars stars
        order = np.argsort(s_star[near])
        s_sorted = s_star[near][order]
        t_sorted = t_star[near][order]
        a = 0
        while a < len(s_sorted):
            b = int(np.searchsorted(s_sorted, s_sorted[a] + ds))
            b = max(b, a + min_stars)
            if b > len(s_sorted):
                break
            t = t_sorted[a:b]
            keep = np.ones(len(t), bool)
            for _ in range(5):  # iterative 2.5-sigma clipping about the median
                med = np.median(t[keep])
                mad = 1.4826 * np.median(abs(t[keep] - med))
                new = abs(t - med) < 2.5 * max(mad, 2e-3)
                if new.sum() == keep.sum():
                    break
                keep = new
            med = np.median(t[keep])
            mad = 1.4826 * np.median(abs(t[keep] - med))
            sc = np.median(s_sorted[a:b][keep])
            rows.append((sc, med, 1.2533 * mad / np.sqrt(keep.sum()), keep.sum(), mad))
            a = b
        rows = np.array(rows)
        sc, tm = rows[:, 0], rows[:, 1]
        # ridge point = guide point at s + t along the guide's normal there
        gxs, gys = np.interp(sc, s_nodes, gx), np.interp(sc, s_nodes, gy)
        nxs, nys = np.interp(sc, s_nodes, nx), np.interp(sc, s_nodes, ny)
        rx, ry = gxs + tm * nxs, gys + tm * nys
        if p < passes - 1:
            # the smoothed ridge becomes the next guide
            gx = gaussian_filter1d(rx, 1.0, mode="nearest")
            gy = gaussian_filter1d(ry, 1.0, mode="nearest")
    norm = np.hypot(nxs, nys)
    c_r, m_r = rx / K, ry
    i_to = int(np.argmin(np.where(m_r > np.median(m_r) - 3, c_r, np.inf)))
    return Ridge(s=sc, c=c_r, m=m_r, sig=rows[:, 2], nx=nxs / norm, ny=nys / norm,
                 n=rows[:, 3].astype(int), i_to=i_to, width=rows[:, 4])


def turnoff_and_sgb_point(c: np.ndarray, m: np.ndarray, offset: float = 0.05):
    """Turnoff (bluest point) and the subgiant point `offset` mag redder than it.

    The classic age indicator of Chaboyer et al. (1996) and Ying et al. (2023):
    the turnoff itself is vertical in the CMD, so its magnitude is ill-defined,
    but the point on the subgiant branch 0.05 mag redder is sharp.
    Points must be ordered from the lower main sequence towards the giant branch.
    """
    ok = np.isfinite(c) & np.isfinite(m)
    c, m = c[ok], m[ok]
    i_to = int(np.argmin(c))
    # parabola through the bluest points refines the turnoff colour
    lo, hi = max(i_to - 3, 0), min(i_to + 4, len(c))
    if hi - lo >= 3:
        p = np.polyfit(m[lo:hi], c[lo:hi], 2)
        if p[0] > 0:
            m_to = -p[1] / (2 * p[0])
            c_to = np.polyval(p, m_to)
        else:
            m_to, c_to = m[i_to], c[i_to]
    else:
        m_to, c_to = m[i_to], c[i_to]
    after = np.arange(i_to, len(c))
    target = c_to + offset
    j = after[np.argmax(c[after] >= target)] if np.any(c[after] >= target) else None
    if j is None or j == 0:
        return float(c_to), float(m_to), np.nan
    f = (target - c[j - 1]) / (c[j] - c[j - 1])
    m_sgb = m[j - 1] + f * (m[j] - m[j - 1])
    return float(c_to), float(m_to), float(m_sgb)
