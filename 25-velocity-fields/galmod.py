"""Model galaxies, "observed" the way THINGS observed the real ones (labelled simulation).

A radio telescope sees a galaxy through a beam: each pixel's spectrum mixes
gas from a patch ~10-15 arcsec across. Where the rotation speed climbs
steeply, the mix holds gas at many velocities, and the fitted velocity is
pulled toward the slower gas: *beam smearing* flattens inner rotation curves,
which is exactly what can make a dense cusp of dark matter look like a
shallow core. The test (after de Blok et al. 2008, Sect. 3.3, who used the
GIPSY task GALMOD) is to build a model data cube with a rotation curve we
choose, smooth it with the real beam, add noise like the real cube's, and run
it through the identical pipeline: Gauss-Hermite fits, the same quality cuts,
the same rings. Whatever the pipeline does to a known curve, it does to the
real one.

The model disc has the observed gas distribution (ring averages of the THINGS
moment-0 map), the measured geometry, a sech^2 vertical layer and a Gaussian
velocity dispersion. Each sky pixel is subsampled 3x3 and integrated along
the line of sight through the thick layer; each channel collects the emission
falling within it (erf differences, not point samples).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import fftconvolve
from scipy.special import erf

import things
import tiltedring as TR


@dataclass
class ModelSpec:
    centre: tuple[float, float]  # arcsec in the map frame
    vsys: float
    pa: float
    incl: float
    r_arcsec: np.ndarray
    vrot: np.ndarray  # km/s at r_arcsec (V = 0 at R = 0 is implied)
    s0: np.ndarray  # face-on integrated intensity, Jy/beam km/s, at r_arcsec
    dispersion: float  # km/s
    z0_arcsec: float  # sech^2 scale height


def beam_kernel(bmaj, bmin, bpa_deg, pixel) -> np.ndarray:
    """Elliptical Gaussian, unit sum, on the pixel grid (FWHMs in arcsec)."""
    sx, sy = bmin / pixel / 2.3548, bmaj / pixel / 2.3548  # before rotation: major along y
    half = int(np.ceil(4 * max(sx, sy)))
    yy, xx = np.mgrid[-half:half + 1, -half:half + 1].astype(float)
    t = np.radians(bpa_deg)
    # position angle runs from north (y) toward east (-x)
    u = xx * np.cos(t) + yy * np.sin(t)
    w = -xx * np.sin(t) + yy * np.cos(t)
    k = np.exp(-0.5 * ((u / sx) ** 2 + (w / sy) ** 2))
    return k / k.sum()


def make_cube(spec: ModelSpec, shape, x_of_col, y_of_row, velocity, *, sub=3, n_z=12,
              block=24) -> np.ndarray:
    """Noise-free model cube (channel, y, x) in Jy/beam before beam smoothing.

    ``x_of_col``/``y_of_row`` give the arcsec coordinate of each pixel
    column/row. Built in blocks of rows to keep memory modest.
    """
    ny, nx = shape
    pix = abs(x_of_col[1] - x_of_col[0])
    offs = (np.arange(sub) - (sub - 1) / 2) / sub * pix
    inc, pa = np.radians(spec.incl), np.radians(spec.pa)
    zq, wq = np.polynomial.legendre.leggauss(n_z)
    z = zq * 4 * spec.z0_arcsec
    wz = wq / np.cosh(z / spec.z0_arcsec) ** 2 if spec.z0_arcsec > 0 else wq
    wz = wz / wz.sum()
    edges = np.concatenate([[velocity[0] - 0.5 * (velocity[1] - velocity[0])],
                            0.5 * (velocity[1:] + velocity[:-1]),
                            [velocity[-1] + 0.5 * (velocity[-1] - velocity[-2])]])
    width = np.abs(np.diff(edges))
    r_tab = np.concatenate([[0.0], spec.r_arcsec])
    v_tab = np.concatenate([[0.0], spec.vrot])
    r_max = spec.r_arcsec.max()
    out = np.zeros((len(velocity), ny, nx))
    for y0 in range(0, ny, block):
        rows = y_of_row[y0:y0 + block]
        xs = (x_of_col[None, :, None, None] + offs[None, None, None, :])
        ys = (rows[:, None, None, None] + offs[None, None, :, None])
        xs = np.broadcast_to(xs, (len(rows), nx, sub, sub)).reshape(len(rows), nx, -1)
        ys = np.broadcast_to(ys, (len(rows), nx, sub, sub)).reshape(len(rows), nx, -1)
        dx, dy = xs - spec.centre[0], ys - spec.centre[1]
        along = -dx * np.sin(pa) + dy * np.cos(pa)  # sky offset along the receding major axis
        across = -dx * np.cos(pa) - dy * np.sin(pa)  # sky offset along the minor axis
        acc = np.zeros((len(velocity), len(rows), nx))
        for zz, ww in zip(z, wz):
            # where this line of sight crosses height zz, in the disc plane
            w_plane = (across - zz * np.sin(inc)) / np.cos(inc)
            r = np.hypot(along, w_plane)
            cos_t = np.where(r > 0, along / np.where(r > 0, r, 1), 0.0)
            s0 = np.interp(r, spec.r_arcsec, spec.s0, left=spec.s0[0], right=0.0)
            s0[r > r_max] = 0.0
            if not s0.any():
                continue
            mom0 = ww * s0 / np.cos(inc)  # projected intensity of this slice, Jy/beam km/s
            vlos = spec.vsys + np.interp(r, r_tab, v_tab) * np.sin(inc) * cos_t
            cdf = 0.5 * (1 + erf((edges[:, None, None, None] - vlos[None]) / (np.sqrt(2) * spec.dispersion)))
            frac = np.abs(np.diff(cdf, axis=0)) / width[:, None, None, None]  # per km/s
            acc += np.mean(frac * mom0[None], axis=-1)
        out[:, y0:y0 + block, :] = acc
    return out


def observe(model: np.ndarray, beam, pixel, noise, rng) -> np.ndarray:
    """Smooth each channel with the beam and add beam-correlated noise."""
    kern = beam_kernel(*beam, pixel)
    smooth = np.stack([fftconvolve(ch, kern, mode="same") for ch in model])
    if noise > 0:
        white = rng.standard_normal(model.shape)
        corr = np.stack([fftconvolve(ch, kern, mode="same") for ch in white])
        corr *= noise / corr.std()
        smooth = smooth + corr
    return smooth


@dataclass
class Recovery:
    radius: np.ndarray
    v_input: np.ndarray
    v_recovered: np.ndarray
    v_recovered_err: np.ndarray
    n_pixels: int


class _NoiseInfo:
    def __init__(self, noise, dv):
        self.noise, self.channel_width = noise, dv


def recover(cube_obs: np.ndarray, velocity, noise, x, y, spec: ModelSpec, ring_step, *,
            mom0_noiseless: np.ndarray) -> Recovery:
    """Run the real pipeline on a model cube; return input vs recovered V(R).

    Geometry is held at the model's own (known) values, so differences come
    from the line profiles and the beam alone, as in de Blok et al.'s test.
    """
    candidates = mom0_noiseless > 3 * noise * abs(velocity[1] - velocity[0])
    maps = things.fit_cube(cube_obs, velocity, noise, candidates)
    maps["mom0_things"][:] = mom0_noiseless
    ok = things.quality_mask(maps, _NoiseInfo(noise, abs(velocity[1] - velocity[0])))
    vf = np.where(ok, maps["gh_vpeak"], np.nan)
    rows, cols = np.nonzero(np.isfinite(vf))
    xs, ys, vs = x[cols], y[rows], vf[rows, cols]
    radius = ring_step * np.arange(1, int(spec.r_arcsec.max() / ring_step) + 1)
    vin = np.interp(radius, np.concatenate([[0.0], spec.r_arcsec]), np.concatenate([[0.0], spec.vrot]))
    vrec, verr = np.full(len(radius), np.nan), np.full(len(radius), np.nan)
    geo = {"x0": spec.centre[0], "y0": spec.centre[1], "vsys": spec.vsys, "pa": spec.pa, "incl": spec.incl}
    for k, r in enumerate(radius):
        vrec[k], _, disp, n = TR.fit_vrot_fixed(xs, ys, vs, geo, r, ring_step)
        verr[k] = disp
    return Recovery(radius, vin, vrec, verr, int(len(vs)))
