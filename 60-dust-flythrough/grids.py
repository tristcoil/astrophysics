"""Cartesian products of the dust maps: projections, slices and viewer volumes.

The main map is sampled once on 2.5 pc cells (x, y within 1250 pc, |z| < 600
pc), plane by plane, and everything is accumulated from those planes:

* ``proj_xy``: extinction integrated over |z| <= 400 pc (the box of
  Edenhofer et al. 2024, Fig. 4), ``proj_xz`` / ``proj_yz``: integrated over y
  / x through the whole volume. Units: ZGR23 mag (x 2.8 for A_V).
* ``slice_z0``: density in the Galactic mid-plane (mean of the two planes at
  z = +-1.25 pc); ``std_z0`` its posterior standard deviation.
* ``vol5``: 5 pc cells, the 2.5 pc cells block-averaged 2x2x2 (averaging keeps
  column densities right; point-sampling a coarse grid would alias thin
  filaments).
* ``inner``: the 2.5 pc cells within |x|, |y| < 450 pc, |z| < 300 pc, the
  volume the fly-through crosses (Taurus, Perseus, Orion).

Cells inside 69 pc and outside the map are NaN in the sampler; they are
stored as 0 here (``nan_*`` masks record where). The 2 kpc validation map
gets the same treatment at 5 pc (x, y within 2000 pc, |z| < 500 pc).
"""

from __future__ import annotations

import time

import numpy as np

import config as C
from dustmap import DustSphere

MAIN = {"h": 2.5, "xy": 1250.0, "z": 600.0, "proj_z": 400.0, "inner_xy": 450.0, "inner_z": 300.0}
WIDE = {"h": 5.0, "xy": 2000.0, "z": 500.0}


def centres(half: float, h: float) -> np.ndarray:
    """Cell centres of [-half, half) in steps of h."""
    return np.arange(-half + h / 2, half, h)


def _block2(a: np.ndarray) -> np.ndarray:
    n0, n1 = a.shape[0] // 2, a.shape[1] // 2
    return a[:2 * n0, :2 * n1].reshape(n0, 2, n1, 2).mean(axis=(1, 3))


def build_main(force: bool = False) -> dict:
    path = C.CACHE_DIR / "grid_main.npz"
    if path.exists() and not force:
        return dict(np.load(path))
    t0 = time.time()
    g = MAIN
    h = g["h"]
    xs, zs = centres(g["xy"], h), centres(g["z"], h)
    nx, nz = xs.size, zs.size
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    X, Y = X.ravel(), Y.ravel()
    sphere = DustSphere("main")
    proj_xy = np.zeros((nx, nx))
    proj_xz = np.zeros((nx, nz))
    proj_yz = np.zeros((nx, nz))
    vol5 = np.zeros((nx // 2, nx // 2, nz // 2), dtype=np.float32)
    ii = np.flatnonzero(np.abs(xs) < g["inner_xy"])
    kk = np.flatnonzero(np.abs(zs) < g["inner_z"])
    inner = np.zeros((ii.size, ii.size, kk.size), dtype=np.float32)
    nan_xy = np.zeros((nx, nx), dtype=bool)
    k0 = np.argmin(np.abs(zs + h / 2))              # plane at z = -h/2
    slice_z0 = np.zeros((nx, nx), dtype=np.float32)
    for k, z in enumerate(zs):
        plane = sphere.density_xyz(X, Y, np.full_like(X, z)).reshape(nx, nx)
        nan = np.isnan(plane)
        plane = np.where(nan, 0.0, plane).astype(np.float32)
        if abs(z) <= g["proj_z"]:
            proj_xy += plane * h
        proj_xz[:, k] = plane.sum(axis=1) * h
        proj_yz[:, k] = plane.sum(axis=0) * h
        vol5[:, :, k // 2] += _block2(plane) / 2
        if k in (k0, k0 + 1):
            slice_z0 += plane / 2
            nan_xy |= nan
        if kk[0] <= k <= kk[-1]:
            inner[:, :, k - kk[0]] = plane[np.ix_(ii, ii)]
        if k % 40 == 0:
            print(f"   main grid: plane {k + 1}/{nz} ({time.time() - t0:.0f} s)", flush=True)
    sphere.free("log_mean")
    zc = np.full(nx * nx, -h / 2)
    std_z0 = 0.5 * (sphere.density_xyz(X, Y, zc, what="std") + sphere.density_xyz(X, Y, -zc, what="std"))
    out = {"xs": xs, "zs": zs, "h": np.float64(h), "proj_xy": proj_xy, "proj_xz": proj_xz,
           "proj_yz": proj_yz, "slice_z0": slice_z0, "std_z0": std_z0.reshape(nx, nx).astype(np.float32),
           "nan_z0": nan_xy, "vol5": vol5, "vol5_xs": centres(g["xy"], 2 * h), "vol5_zs": centres(g["z"], 2 * h),
           "inner": inner, "inner_xs": xs[ii], "inner_zs": zs[kk]}
    np.savez(path, **out)
    print(f"   main grid done in {time.time() - t0:.0f} s -> {path.name}", flush=True)
    return out


def build_wide(force: bool = False) -> dict:
    """The 2 kpc validation map on 5 pc cells, for the far end of the Radcliffe Wave."""
    path = C.CACHE_DIR / "grid_2kpc.npz"
    if path.exists() and not force:
        return dict(np.load(path))
    t0 = time.time()
    g = WIDE
    h = g["h"]
    xs, zs = centres(g["xy"], h), centres(g["z"], h)
    nx, nz = xs.size, zs.size
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    X, Y = X.ravel(), Y.ravel()
    sphere = DustSphere("2kpc")
    proj_xy = np.zeros((nx, nx))
    proj_xz = np.zeros((nx, nz))
    proj_yz = np.zeros((nx, nz))
    vol10 = np.zeros((nx // 2, nx // 2, nz // 2), dtype=np.float32)
    for k, z in enumerate(zs):
        plane = sphere.density_xyz(X, Y, np.full_like(X, z)).reshape(nx, nx)
        plane = np.where(np.isnan(plane), 0.0, plane).astype(np.float32)
        proj_xy += plane * h
        proj_xz[:, k] = plane.sum(axis=1) * h
        proj_yz[:, k] = plane.sum(axis=0) * h
        vol10[:, :, k // 2] += _block2(plane) / 2
        if k % 40 == 0:
            print(f"   2 kpc grid: plane {k + 1}/{nz} ({time.time() - t0:.0f} s)", flush=True)
    out = {"xs": xs, "zs": zs, "h": np.float64(h), "proj_xy": proj_xy, "proj_xz": proj_xz,
           "proj_yz": proj_yz, "vol10": vol10, "vol10_xs": centres(g["xy"], 2 * h),
           "vol10_zs": centres(g["z"], 2 * h)}
    np.savez(path, **out)
    print(f"   2 kpc grid done in {time.time() - t0:.0f} s -> {path.name}", flush=True)
    return out


if __name__ == "__main__":
    import sys
    build_main(force="--force" in sys.argv)
