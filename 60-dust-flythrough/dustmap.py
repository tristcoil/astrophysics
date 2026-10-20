"""Read, sample and integrate the Edenhofer et al. (2024) 3D dust map.

The map is a stack of HEALPix spheres (Nside 256, NESTED ordering, 14' pixels)
at logarithmically spaced distances: 516 spheres from 69 to 1250 pc in the
main map. Each voxel holds the posterior mean (and standard deviation) of the
*differential* extinction rho, in ZGR23 extinction per parsec (multiply by
2.8 for A_V per parsec). Inside 69 pc the file gives only the extinction
integrated from the Sun to 69 pc, one number per pixel.

Two conventions, both taken from the authors' own tools:

* **Density at a point** (``density_lbd``, ``density_xyz``, ``to_grid``): the
  *logarithm* of the density is interpolated bilinearly on the sphere
  (healpy's four-pixel weights) and linearly in distance between sphere
  centres, then exponentiated. This is what the release's ``interp2box.py``
  and the ``dustmaps`` package do. Outside the first and last sphere centre
  the density is undefined (NaN).
* **Integrated extinction** (``extinction``): each voxel is a constant density
  between its radial boundaries (the file's "radial pixel boundaries"), added
  to the inner 69-pc term; the four-pixel weights interpolate the integrated
  extinction itself. ``dustmaps``' integrated mode instead places the
  cumulative sum to each voxel's *outer* boundary at the voxel's *centre*,
  i.e. half a voxel (0.2-3.5 pc) too close; ``extinction_compare`` measures
  the size of that offset.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property

import healpy as hp
import numpy as np
from astropy.io import fits

import config as C


@dataclass
class DustSphere:
    """One reconstruction ("main" or "2kpc") as HEALPix times distance."""

    which: str = "main"
    path: object = None                        # override the file (tests use a small synthetic map)
    nside: int = field(init=False)
    nest: bool = field(init=False)
    radii: np.ndarray = field(init=False)      # sphere centres (pc), increasing
    bounds: np.ndarray = field(init=False)     # radial voxel boundaries (pc), len(radii) + 1
    inner: np.ndarray = field(init=False)      # extinction integrated within bounds[0]
    inner_std: np.ndarray | None = field(init=False)
    r_inner: float = field(init=False)         # radius of that inner sphere (pc)
    units: str = field(init=False)

    def __post_init__(self):
        if self.path is None:
            self.path = C.EDENHOFER_DIR / C.EDENHOFER_FILES[self.which][0]
        self._hdul = fits.open(self.path, memmap=True, lazy_load_hdus=False)
        inner = inner_std = None
        for hdu in self._hdul:
            name = hdu.name.lower()
            if name == "mean":
                self.nside = int(hdu.header["NSIDE"])
                self.nest = hdu.header["ORDERING"].lower().startswith("nest")
                self.units = hdu.header.get("CUNIT", "")
                self._mean_hdu = hdu
            elif name == "std.":
                self._std_hdu = hdu
            elif isinstance(hdu, fits.BinTableHDU):
                col = hdu.data.names[0]
                if col == "radial pixel centers":
                    self.radii = np.asarray(hdu.data[col], dtype=np.float64)
                elif col == "radial pixel boundaries":
                    self.bounds = np.asarray(hdu.data[col], dtype=np.float64)
            elif name.startswith("mean of integrated inner") or name.startswith("integrated inner"):
                ctype = (hdu.header.get("CTYPE") or hdu.header.get("CTYPE1")).lower()
                self.r_inner = float(ctype.split("within")[1].replace("pc", ""))
                inner = np.asarray(hdu.data, dtype=np.float32)
            elif name.startswith("std. of integrated inner"):
                inner_std = np.asarray(hdu.data, dtype=np.float32)
        self.inner = inner
        self.inner_std = inner_std
        assert self.bounds.size == self.radii.size + 1
        assert np.all(np.diff(self.radii) > 0)

    # ------------------------------------------------------------------
    # Raw arrays
    # ------------------------------------------------------------------

    @property
    def npix(self) -> int:
        return 12 * self.nside**2

    @property
    def nshell(self) -> int:
        return self.radii.size

    @cached_property
    def mean(self) -> np.ndarray:
        """Posterior-mean density, (nshell, npix) float32 in native byte order (1.6 GB)."""
        return np.asarray(self._mean_hdu.data, dtype=np.float32)

    @cached_property
    def std(self) -> np.ndarray:
        """Posterior standard deviation, same shape as ``mean``."""
        return np.asarray(self._std_hdu.data, dtype=np.float32)

    @cached_property
    def log_mean(self) -> np.ndarray:
        """ln(mean), flattened, for log-space interpolation (read shell by shell)."""
        out = np.empty((self.nshell, self.npix), dtype=np.float32)
        src = self._mean_hdu.data
        for i in range(self.nshell):
            np.log(src[i], out=out[i])
        return out.ravel()

    def free(self, *names: str) -> None:
        """Drop cached arrays to save memory."""
        for name in names or ("mean", "std", "log_mean", "cumulative"):
            self.__dict__.pop(name, None)

    @cached_property
    def cumulative(self) -> np.ndarray:
        """Extinction from the Sun to each radial boundary, (nshell + 1, npix) float32."""
        cum = np.empty((self.nshell + 1, self.npix), dtype=np.float32)
        acc = self.inner.astype(np.float64)
        cum[0] = acc
        widths = np.diff(self.bounds)
        for i in range(self.nshell):
            acc += widths[i] * self.mean[i]
            cum[i + 1] = acc
        return cum

    # ------------------------------------------------------------------
    # Interpolation
    # ------------------------------------------------------------------

    def _weights(self, l_deg: np.ndarray, b_deg: np.ndarray):
        pix, w = hp.get_interp_weights(self.nside, l_deg, b_deg, nest=self.nest, lonlat=True)
        return pix, w

    def density_lbd(self, l_deg, b_deg, d_pc, *, chunk: int = 2_000_000,
                    what: str = "mean") -> np.ndarray:
        """Density (ZGR23 mag/pc) at Galactic (l, b, distance), NaN outside the map.

        ``what="std"`` interpolates the variance linearly (as interp2box.py)
        and returns its square root.
        """
        l_deg, b_deg, d_pc = (np.asarray(a, dtype=np.float64) for a in np.broadcast_arrays(l_deg, b_deg, d_pc))
        shape = d_pc.shape
        l_deg, b_deg, d_pc = l_deg.ravel(), b_deg.ravel(), d_pc.ravel()
        out = np.full(d_pc.size, np.nan, dtype=np.float32)
        if what == "mean":
            table = self.log_mean
        else:
            table = (self.std.astype(np.float32) ** 2).ravel()
        for s in range(0, d_pc.size, chunk):
            sl = slice(s, s + chunk)
            d = d_pc[sl]
            ok = (d >= self.radii[0]) & (d <= self.radii[-1])
            if not ok.any():
                continue
            ir = np.searchsorted(self.radii, d[ok], side="right") - 1
            ir = np.clip(ir, 0, self.nshell - 2)
            t = (d[ok] - self.radii[ir]) / (self.radii[ir + 1] - self.radii[ir])
            pix, w = self._weights(l_deg[sl][ok], b_deg[sl][ok])
            lo = (table[ir * self.npix + pix] * w).sum(axis=0)
            hi = (table[(ir + 1) * self.npix + pix] * w).sum(axis=0)
            val = (1.0 - t) * lo + t * hi
            res = np.exp(val) if what == "mean" else np.sqrt(np.maximum(val, 0.0))
            tmp = out[sl]
            tmp[ok] = res
            out[sl] = tmp
        return out.reshape(shape)

    def density_xyz(self, x, y, z, **kw) -> np.ndarray:
        l, b, d = xyz_to_lbd(x, y, z)
        return self.density_lbd(l, b, d, **kw)

    def extinction(self, l_deg, b_deg, d_pc) -> np.ndarray:
        """Extinction (ZGR23 mag) from the Sun to distance ``d_pc`` along (l, b).

        Inside ``r_inner`` (69 pc) the inner term is scaled linearly with
        distance (the file gives only its total). Beyond the last boundary the
        value saturates at the map's total.
        """
        l_deg, b_deg, d_pc = (np.asarray(a, dtype=np.float64) for a in np.broadcast_arrays(l_deg, b_deg, d_pc))
        shape = d_pc.shape
        l_deg, b_deg, d_pc = l_deg.ravel(), b_deg.ravel(), d_pc.ravel()
        pix, w = self._weights(l_deg, b_deg)
        k = np.searchsorted(self.bounds, d_pc, side="right") - 1        # voxel containing d
        inside = k < 0
        beyond = k >= self.nshell
        kk = np.clip(k, 0, self.nshell - 1)
        base = self.cumulative[kk[None, :], pix]
        dens = self.mean[kk[None, :], pix]
        part = np.clip(d_pc - self.bounds[kk], 0.0, None)
        a = base + dens * part
        a = np.where(beyond[None, :], self.cumulative[-1][pix], a)
        a = np.where(inside[None, :], self.inner[pix] * (np.clip(d_pc, 0, None) / self.bounds[0]), a)
        return (a * w).sum(axis=0).reshape(shape)

    def column_profiles(self, pix: np.ndarray, d_grid: np.ndarray) -> np.ndarray:
        """Density at pixel *centres* along their lines of sight, (len(d_grid), len(pix)).

        At a pixel centre the four-pixel weights reduce to the pixel itself,
        so this is the log-linear radial interpolation of one voxel column,
        exactly what dustmaps returns there.
        """
        d_grid = np.asarray(d_grid, dtype=np.float64)
        ir = np.clip(np.searchsorted(self.radii, d_grid, side="right") - 1, 0, self.nshell - 2)
        t = ((d_grid - self.radii[ir]) / (self.radii[ir + 1] - self.radii[ir]))[:, None]
        lm = self.log_mean
        lo = lm[ir[:, None] * self.npix + pix[None, :]]
        hi = lm[(ir[:, None] + 1) * self.npix + pix[None, :]]
        out = np.exp((1.0 - t) * lo + t * hi).astype(np.float32)
        bad = (d_grid < self.radii[0]) | (d_grid > self.radii[-1])
        out[bad] = np.nan
        return out

    # ------------------------------------------------------------------
    # Cartesian grids
    # ------------------------------------------------------------------

    def to_grid(self, extent, shape, *, what: str = "mean", verbose: bool = False) -> np.ndarray:
        """Resample onto a regular heliocentric Galactic (x, y, z) grid.

        ``extent`` = ((x0, x1), (y0, y1), (z0, z1)) are the centres of the first
        and last voxels (endpoints included, as in interp2box.py); ``shape`` =
        (nx, ny, nz). Returns float32 (nx, ny, nz); NaN outside the map.
        """
        axes = [np.linspace(lo, hi, n) for (lo, hi), n in zip(extent, shape)]
        out = np.empty(tuple(shape), dtype=np.float32)
        X, Y = np.meshgrid(axes[0], axes[1], indexing="ij")
        X, Y = X.ravel(), Y.ravel()
        for iz, z in enumerate(axes[2]):
            out[:, :, iz] = self.density_xyz(X, Y, np.full_like(X, z), what=what).reshape(shape[0], shape[1])
            if verbose and iz % 20 == 0:
                print(f"   grid plane {iz + 1}/{shape[2]}", flush=True)
        return out


def xyz_to_lbd(x, y, z):
    """Heliocentric Galactic Cartesian (pc) -> (l deg, b deg, distance pc)."""
    x, y, z = (np.asarray(a, dtype=np.float64) for a in (x, y, z))
    rho = np.hypot(x, y)
    return (np.degrees(np.arctan2(y, x)) % 360.0, np.degrees(np.arctan2(z, rho)),
            np.sqrt(rho**2 + z**2))


def lbd_to_xyz(l_deg, b_deg, d_pc):
    l, b = np.radians(l_deg), np.radians(b_deg)
    d = np.asarray(d_pc, dtype=np.float64)
    return d * np.cos(b) * np.cos(l), d * np.cos(b) * np.sin(l), d * np.sin(b)


def to_av(e_zgr23):
    """ZGR23 extinction (or per-pc density) -> A_V at 540 nm (Edenhofer's convention)."""
    return C.ZGR23_TO_AV * np.asarray(e_zgr23)


def to_nh(rho_zgr23_per_pc):
    """ZGR23 per pc -> hydrogen nuclei per cm^3 (O'Neill et al. 2024, eq. 4)."""
    return C.ZGR23_TO_NH * np.asarray(rho_zgr23_per_pc)
