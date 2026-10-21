"""Reconstruction: undoing the bulk flows that blur the BAO peak.

Over billions of years galaxies have moved ~5-10 h^-1 Mpc from where they
formed, pulled by the same density field they trace. Those motions smear the
100 h^-1 Mpc acoustic peak. Because the motions are caused by the density
field we can see, they can be estimated and undone (Eisenstein et al. 2007).

In linear theory (Zel'dovich) a galaxy's displacement Psi obeys
del . Psi = -delta_matter. Galaxies trace matter with a bias b, and in
redshift space their line-of-sight motions add a term f (Psi . r^) r^, so

    del . Psi + (f/b) del . [(Psi . r^) r^] = -delta_galaxies / b

(Nusser & Davis 1994; Padmanabhan et al. 2012). We solve it on a 3D mesh with
FFTs, treating the line-of-sight term iteratively (Burden et al. 2015), which
handles the survey's wide angles exactly instead of assuming a fixed line of
sight. Then every galaxy and every random point is moved back by
Psi + f (Psi . r^) r^ ("RecSym", the convention DESI used), and the correlation
function is measured again as xi = (DD - 2DS + SS) / RR, where S are the
shifted randoms and R the original ones (Padmanabhan et al. 2012).

Settings follow DESI DR1 for LRGs: b = 2.0, f = 0.8, 15 h^-1 Mpc Gaussian
smoothing. The mesh has 8 h^-1 Mpc cells; fields are float32 to save memory.
"""

from __future__ import annotations

import time

import numpy as np
import scipy.fft as sfft

CELL = 8.0
SMOOTH = 15.0
BIAS = 2.0
GROWTH = 0.8
NITER = 8   # the line-of-sight term converges as (f/b)^n: 0.07% after 8
PAD = 120.0
RAN_MIN = 0.01  # cells with fewer randoms than this x the mean count as empty


class Mesh:
    """A cubic-cell mesh in comoving h^-1 Mpc enclosing a set of points."""

    def __init__(self, pos: np.ndarray, cell: float = CELL, pad: float = PAD):
        lo, hi = pos.min(0) - pad, pos.max(0) + pad
        self.cell = cell
        self.origin = lo
        self.shape = tuple(sfft.next_fast_len(int(np.ceil((h - l) / cell)), real=True)
                           for l, h in zip(lo, hi))
        self.size = int(np.prod(self.shape))

    def _corners(self, pos):
        u = (pos - self.origin) / self.cell - 0.5  # cell centres at (i + 1/2) cell
        i0 = np.floor(u).astype(np.int64)
        f = (u - i0).astype(np.float32)
        n0, n1, n2 = self.shape
        for dx in (0, 1):
            for dy in (0, 1):
                for dz in (0, 1):
                    ix = np.clip(i0[:, 0] + dx, 0, n0 - 1)
                    iy = np.clip(i0[:, 1] + dy, 0, n1 - 1)
                    iz = np.clip(i0[:, 2] + dz, 0, n2 - 1)
                    w = ((f[:, 0] if dx else 1 - f[:, 0]) * (f[:, 1] if dy else 1 - f[:, 1])
                         * (f[:, 2] if dz else 1 - f[:, 2]))
                    yield (ix * n1 + iy) * n2 + iz, w

    def paint(self, pos, weight, chunk: int = 4_000_000) -> np.ndarray:
        """Cloud-in-cell assignment of weighted points."""
        out = np.zeros(self.size, np.float32)
        for start in range(0, len(pos), chunk):
            p, wt = pos[start:start + chunk], weight[start:start + chunk]
            for idx, w in self._corners(p):
                out += np.bincount(idx, weights=w * wt, minlength=self.size).astype(np.float32)
        return out.reshape(self.shape)

    def read(self, field: np.ndarray, pos, chunk: int = 4_000_000) -> np.ndarray:
        """Trilinear (cloud-in-cell) interpolation of a mesh field at points."""
        flat = field.ravel()
        out = np.empty(len(pos), np.float32)
        for start in range(0, len(pos), chunk):
            acc = np.zeros(min(chunk, len(pos) - start), np.float32)
            for idx, w in self._corners(pos[start:start + chunk]):
                acc += flat[idx] * w
            out[start:start + chunk] = acc
        return out

    def kvectors(self):
        n0, n1, n2 = self.shape
        kx = (2 * np.pi * sfft.fftfreq(n0, d=self.cell)).astype(np.float32)[:, None, None]
        ky = (2 * np.pi * sfft.fftfreq(n1, d=self.cell)).astype(np.float32)[None, :, None]
        kz = (2 * np.pi * sfft.rfftfreq(n2, d=self.cell)).astype(np.float32)[None, None, :]
        return kx, ky, kz

    def unit_los(self):
        """r^ at every cell centre (the observer sits at the coordinate origin)."""
        axes = [self.origin[a] + (np.arange(self.shape[a], dtype=np.float32) + 0.5) * self.cell
                for a in range(3)]
        x, y, z = axes[0][:, None, None], axes[1][None, :, None], axes[2][None, None, :]
        r = np.sqrt(x**2 + y**2 + z**2)
        return x / r, y / r, z / r


def displacement(mesh: Mesh, delta_k: np.ndarray, kvec, workers: int) -> list[np.ndarray]:
    """Psi with del . Psi = -delta: Psi_k = i k delta_k / k^2 (real-space components)."""
    kx, ky, kz = kvec
    k2 = kx**2 + ky**2 + kz**2
    k2[0, 0, 0] = 1.0
    out = []
    for kc in (kx, ky, kz):
        comp = 1j * kc * delta_k / k2
        comp[0, 0, 0] = 0.0
        out.append(sfft.irfftn(comp, s=mesh.shape, workers=workers).astype(np.float32))
    return out


def reconstruct(data_pos, data_w, rand_pos, rand_w, *, f: float = GROWTH, b: float = BIAS,
                smooth: float = SMOOTH, niter: int = NITER, workers: int = 8,
                verbose: bool = True) -> tuple[np.ndarray, np.ndarray, dict]:
    """Displacements (Psi + f (Psi . r^) r^) at the data and random positions."""
    t0 = time.monotonic()
    mesh = Mesh(rand_pos)
    kvec = mesh.kvectors()
    kx, ky, kz = kvec
    gauss = np.exp(-0.5 * (kx**2 + ky**2 + kz**2) * smooth**2).astype(np.float32)

    def smoothed(field):
        return sfft.irfftn(sfft.rfftn(field, workers=workers) * gauss, s=mesh.shape,
                           workers=workers).astype(np.float32)

    D = smoothed(mesh.paint(data_pos, data_w))
    R = smoothed(mesh.paint(rand_pos, rand_w))
    alpha = data_w.sum() / rand_w.sum()
    thr = RAN_MIN * R[R > 0].mean()
    delta = np.where(R > thr, D / (alpha * np.where(R > thr, R, 1.0)) - 1.0, 0.0).astype(np.float32)
    del D, R
    delta_k = sfft.rfftn(delta, workers=workers)
    del delta
    rhat = mesh.unit_los()
    psi = displacement(mesh, delta_k / b, kvec, workers)
    for _ in range(niter - 1):
        # del . [(Psi . r^) r^] by FFT, then re-solve with the redshift-space term moved over
        proj = psi[0] * rhat[0] + psi[1] * rhat[1] + psi[2] * rhat[2]
        div_k = sum(1j * kc * sfft.rfftn(proj * rc, workers=workers)
                    for kc, rc in zip((kx, ky, kz), rhat))
        del proj
        psi = displacement(mesh, (delta_k + f * div_k) / b, kvec, workers)
        del div_k

    def shift(pos):
        comps = np.column_stack([mesh.read(c, pos) for c in psi])
        r = np.linalg.norm(pos, axis=1, keepdims=True)
        los = pos / r
        return comps + f * np.sum(comps * los, axis=1, keepdims=True) * los

    data_shift, rand_shift = shift(data_pos), shift(rand_pos)
    info = {"mesh": mesh.shape, "cell": mesh.cell, "alpha": float(alpha),
            "rms_shift_data": float(np.sqrt(np.mean(np.sum(data_shift**2, 1)))),
            "seconds": time.monotonic() - t0}
    if verbose:
        print(f"      mesh {mesh.shape}, rms displacement {info['rms_shift_data']:.1f} h^-1 Mpc, "
              f"{info['seconds']:.0f} s", flush=True)
    return data_shift, rand_shift, info
