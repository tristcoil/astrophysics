"""Render what an observer near a Schwarzschild black hole sees of the real sky.

For every pixel: viewing direction -> angle chi from the direction to the
black hole -> (raymap) the angle phi the ray swept around the hole, and the
frequency shift g -> the sky direction it came from -> the Gaia sky there,
shifted by g (colour.py). Sky that a pixel sees compressed (rays that wound
round the photon sphere) is averaged over the right patch via the mipmap.

Point stars go the other way: for each star, which viewing directions see it?
Spherical symmetry gives at most one image per winding: the primary
(phi = Theta), the secondary on the opposite side (phi = 2 pi - Theta), and
the faint tertiary (2 pi + Theta). Each image gets its magnification (the
ratio of solid angles, d Omega_seen / d Omega_sky = sin chi / |sin phi dphi/dchi|),
and its brightness and colour change with g.

Geometry: z_hat points from the black hole to the observer. The sky is the
celestial sphere at infinity in Galactic coordinates; the black hole is placed
in front of the direction of Sgr A* as seen from the Sun.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import ndimage

import colour
import config as C
import raymap
import sky as S

ARCSEC2_PER_SR = (180.0 / math.pi * 3600.0) ** 2
#: Luminance per steradian shown at display value 1 (before tone mapping).
I_REF = 10.0 ** (-0.4 * C.EXPOSURE_MAG_ARCSEC2) * ARCSEC2_PER_SR
#: Lensing magnification cap for point images (a star exactly behind the hole
#: would become an Einstein ring of formally infinite magnification).
MU_MAX = 1.0e3


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def sgra_direction() -> np.ndarray:
    return S.radec_to_galvec(C.SGRA_RA_DEG, C.SGRA_DEC_DEG)


@dataclass
class Camera:
    """A stereographic (or equidistant full-sky) camera in the observer's frame."""
    width: int
    height: int
    forward: np.ndarray
    up: np.ndarray
    fov_deg: float = 150.0          # horizontal field of view (stereographic)
    projection: str = "stereographic"

    def basis(self):
        f = unit(self.forward)
        r = unit(np.cross(f, self.up))
        u = np.cross(r, f)
        return f, r, u

    def directions(self):
        """Unit viewing directions (H, W, 3), pixel angular size (H, W), and validity mask."""
        f, r, u = self.basis()
        x = (np.arange(self.width) + 0.5) - self.width / 2
        y = self.height / 2 - (np.arange(self.height) + 0.5)
        X, Y = np.meshgrid(x, y)
        rho = np.hypot(X, Y)
        az = np.arctan2(Y, X)
        valid = np.ones(rho.shape, bool)
        if self.projection == "stereographic":
            fpx = (self.width / 2) / (2.0 * math.tan(math.radians(self.fov_deg) / 4.0))
            theta = 2.0 * np.arctan(rho / (2.0 * fpx))
            dtheta = 1.0 / (fpx * (1.0 + (rho / (2.0 * fpx)) ** 2))
        elif self.projection == "equidistant":      # whole sky in a disc: rim = straight behind
            rmax = min(self.width, self.height) / 2
            theta = rho / rmax * math.pi
            dtheta = np.full(rho.shape, math.pi / rmax)
            valid = theta <= math.pi
            theta = np.minimum(theta, math.pi)
        else:
            raise ValueError(self.projection)
        st = np.sin(theta)
        d = (np.cos(theta)[..., None] * f + st[..., None] * (np.cos(az)[..., None] * r + np.sin(az)[..., None] * u))
        tang = np.where(rho > 0, st / np.maximum(rho, 1e-9), dtheta)
        pix = np.maximum(dtheta, tang)
        return d, pix, valid

    def project(self, d):
        """Pixel coordinates (x, y) of unit directions d (N, 3); NaN if outside the image."""
        f, r, u = self.basis()
        cz = d @ f
        theta = np.arccos(np.clip(cz, -1.0, 1.0))
        az = np.arctan2(d @ u, d @ r)
        if self.projection == "stereographic":
            fpx = (self.width / 2) / (2.0 * math.tan(math.radians(self.fov_deg) / 4.0))
            rho = 2.0 * fpx * np.tan(theta / 2.0)
            rho = np.where(theta < math.pi * 0.999, rho, np.nan)
        else:
            rho = theta / math.pi * (min(self.width, self.height) / 2)
        x = self.width / 2 + rho * np.cos(az) - 0.5
        y = self.height / 2 - rho * np.sin(az) - 0.5
        inside = (x > -1) & (x < self.width) & (y > -1) & (y < self.height)
        return np.where(inside, x, np.nan), np.where(inside, y, np.nan), theta

    def pixel_solid_angle(self, theta):
        """Solid angle of one pixel at angle theta from the centre."""
        if self.projection == "stereographic":
            fpx = (self.width / 2) / (2.0 * math.tan(math.radians(self.fov_deg) / 4.0))
            # conformal: d theta / d rho = 1 / (fpx (1 + (rho/2f)^2)) in both directions, rho/2f = tan(theta/2)
            dth = 1.0 / (fpx * (1.0 + np.tan(theta / 2.0) ** 2))
            return dth * dth
        rmax = min(self.width, self.height) / 2
        dth = math.pi / rmax
        return dth * dth * np.where(theta > 1e-6, np.sin(theta) / np.maximum(theta, 1e-6), 1.0)


class Scene:
    """The sky (points + glow), loaded once and shared by all frames."""

    def __init__(self):
        self.texture = S.SkyTexture.load()
        self.points = S.points()
        self.flux = 10.0 ** (-0.4 * self.points.vmag)
        self._axis = None

    def star_geometry(self, zh):
        """Each star's angle Theta from the fall axis and its azimuth unit vector (cached per axis)."""
        key = tuple(np.round(zh, 12))
        if self._axis != key:
            vec = self.points.vec
            cos_t = vec @ zh
            perp = vec - cos_t[:, None] * zh
            nperp = np.linalg.norm(perp, axis=1)
            good = nperp > 1e-12
            e_star = np.zeros_like(vec)
            e_star[good] = perp[good] / nperp[good, None]
            self._geom = (np.arccos(np.clip(cos_t, -1.0, 1.0)), e_star, good)
            self._axis = key
        return self._geom


def _shadow_coverage(chi, chi_sh, pix):
    """Fraction of a pixel outside the shadow (analytic anti-aliasing of the edge)."""
    return np.clip((chi - chi_sh) / np.maximum(pix, 1e-9) + 0.5, 0.0, 1.0)


def render_view(scene: Scene, cam: Camera, table: raymap.RayTable, zhat: np.ndarray,
                stars: bool = True, glow: bool = True, star_mag_limit: float = 99.0):
    """Linear-RGB image (H, W, 3) in units of I_REF (display value before tone mapping)."""
    d, pix, valid = cam.directions()
    H, W = pix.shape
    img = np.zeros((H, W, 3), np.float64)
    zh = unit(zhat)
    cos_chi = -(d @ zh)
    chi = np.arccos(np.clip(cos_chi, -1.0, 1.0))
    sin_chi = np.sin(chi)
    if glow:
        phi, g, dphi = table.lookup(chi.ravel())
        phi, g, dphi = phi.reshape(chi.shape), g.reshape(chi.shape), dphi.reshape(chi.shape)
        ok = np.isfinite(phi) & valid
        e_alpha = (d + cos_chi[..., None] * zh) / np.maximum(sin_chi, 1e-12)[..., None]
        # straight ahead / behind: any perpendicular works (phi = 0 or pi there)
        s = unit(np.cos(phi)[..., None] * zh + np.sin(phi)[..., None] * e_alpha)
        stretch = np.maximum(np.abs(dphi), np.abs(np.sin(phi)) / np.maximum(sin_chi, 1e-9))
        footprint = pix * stretch
        lum = np.zeros(chi.shape)
        lt = np.full(chi.shape, 3.7)
        lum[ok], lt[ok] = scene.texture.sample(s[ok], footprint[ok])
        rgb = colour.shifted_rgb(lum[ok], lt[ok], g[ok]) / I_REF
        img[ok] = rgb
        cov = _shadow_coverage(chi, table.chi_sh, pix)
        img *= cov[..., None]
        img[~valid] = 0.0
    if stars:
        img += render_stars(scene, cam, table, zh, star_mag_limit)
    return img


def star_images(scene: Scene, table: raymap.RayTable, zh: np.ndarray, max_order: int = 2,
                min_flux: float = 0.0):
    """Images of every point star: yields (order, d_hat, mu, g, idx) for images whose
    lensed, shifted luminance flux * mu * Y(gT)/Y(T) exceeds ``min_flux``."""
    theta, e_star, good = scene.star_geometry(zh)
    for order in range(max_order + 1):
        if order % 2 == 0:
            target = theta + math.pi * order
            sign = 1.0
        else:
            target = math.pi * (order + 1) - theta
            sign = -1.0
        chi = table.invert(target)
        ok = np.isfinite(chi) & good
        idx = np.flatnonzero(ok)
        if idx.size == 0:
            continue
        c = chi[idx]
        _, g, dphi = table.lookup(c)
        sin_phi = np.abs(np.sin(target[idx]))
        mu = np.minimum(np.sin(c) / np.maximum(sin_phi * np.abs(dphi), 1e-300), MU_MAX)
        bright = scene.flux[idx] * mu * colour.luminance_ratio(scene.points.logT[idx], g) > min_flux
        idx, c, g, mu = idx[bright], c[bright], g[bright], mu[bright]
        d = -np.cos(c)[:, None] * zh + np.sin(c)[:, None] * (sign * e_star[idx])
        yield order, d, mu, g, idx


def render_stars(scene: Scene, cam: Camera, table: raymap.RayTable, zh: np.ndarray,
                 mag_limit: float = 99.0, max_order: int = 2, threshold: float = 2e-3):
    """Point-star layer in display units. Images fainter than ``threshold`` display units
    in the camera's smallest pixels are skipped (they would round to black)."""
    H, W = cam.height, cam.width
    layer = np.zeros((3, H * W), np.float64)
    omega_min = float(cam.pixel_solid_angle(np.array([math.radians(cam.fov_deg) / 2]))[0]) \
        if cam.projection == "stereographic" else float(cam.pixel_solid_angle(np.array([math.pi / 2]))[0])
    min_flux = threshold * omega_min * I_REF
    lt_all = scene.points.logT
    for order, d, mu, g, idx in star_images(scene, table, zh, max_order, min_flux):
        if mag_limit < 99:
            keep = scene.points.vmag[idx] < mag_limit
            d, mu, g, idx = d[keep], mu[keep], g[keep], idx[keep]
        x, y, theta = cam.project(d)
        ok = np.isfinite(x)
        if not ok.any():
            continue
        x, y, theta = x[ok], y[ok], theta[ok]
        k = idx[ok]
        rgb = colour.shifted_rgb(scene.flux[k] * mu[ok], lt_all[k], g[ok])
        rgb /= (cam.pixel_solid_angle(theta) * I_REF)[:, None]
        x0 = np.floor(x).astype(int)
        y0 = np.floor(y).astype(int)
        fx, fy = x - x0, y - y0
        for dx, dy, w in ((0, 0, (1 - fx) * (1 - fy)), (1, 0, fx * (1 - fy)),
                          (0, 1, (1 - fx) * fy), (1, 1, fx * fy)):
            xx, yy = x0 + dx, y0 + dy
            m = (xx >= 0) & (xx < W) & (yy >= 0) & (yy < H)
            flat = yy[m] * W + xx[m]
            for ch in range(3):
                layer[ch] += np.bincount(flat, weights=rgb[m, ch] * w[m], minlength=H * W)
    layer = layer.reshape(3, H, W)
    # A camera-like point-spread function: a sharp core plus faint wide wings (glare),
    # so that very bright images bloom as they would on film or in an eye.
    out = np.empty((H, W, 3))
    for ch in range(3):
        out[..., ch] = 0.97 * ndimage.gaussian_filter(layer[ch], 0.7) + 0.03 * ndimage.gaussian_filter(layer[ch], 5.0)
    return out


#: Tone curve: display brightness = asinh(L / TONE_K) / asinh(TONE_LMAX / TONE_K), applied to
#: the luminance L (display units) with every pixel's colour kept (channels scaled together,
#: and pulled back below 1 by the largest channel), then sRGB-encoded. A logarithmic-ish
#: stretch is needed because a falling observer's sky spans ~6 decades of brightness.
TONE_K = 0.1
TONE_LMAX = 100.0
_LUMA = np.array([0.2126, 0.7152, 0.0722])


def _srgb_encode(y):
    y = np.clip(y, 0.0, 1.0)
    return np.where(y <= 0.0031308, 12.92 * y, 1.055 * np.power(y, 1 / 2.4) - 0.055)


def _srgb_decode(e):
    e = np.clip(e, 0.0, 1.0)
    return np.where(e <= 0.04045, e / 12.92, ((e + 0.055) / 1.055) ** 2.4)


def tone_map(img, exposure: float = 1.0):
    """Display sRGB (uint8) from linear values in units of I_REF."""
    x = np.clip(np.asarray(img, float) * exposure, 0.0, None)
    L = x @ _LUMA
    target = _srgb_decode(np.minimum(np.arcsinh(L / TONE_K) / np.arcsinh(TONE_LMAX / TONE_K), 1.0))
    out = x * (target / np.maximum(L, 1e-30))[..., None]
    m = out.max(axis=-1, keepdims=True)
    out = np.where(m > 1.0, out / np.maximum(m, 1e-30), out)
    return (_srgb_encode(out) * 255 + 0.5).astype(np.uint8)


def default_frame():
    """The fall's geometry: the faller comes from the Sun's side and looks at Sgr A*.

    Returns zhat (black hole -> observer) and an 'up' vector (Galactic north,
    made perpendicular), both Galactic unit vectors.
    """
    ahead = sgra_direction()
    zhat = -ahead
    north = np.array([0.0, 0.0, 1.0])
    up = unit(north - (north @ ahead) * ahead)
    return zhat, up
