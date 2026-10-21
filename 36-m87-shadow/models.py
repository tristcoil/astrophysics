"""Geometric source models with exact Fourier transforms.

Every model is defined in the image (what the source looks like) and evaluated
in the Fourier domain (what the interferometer measures), using

    V(u, v) = integral I(x, y) exp(+2 pi i (u x + v y)) dx dy,

x east, y north, position angles east of north: the convention of the EHT data
and of ehtim (``test_physics.py`` checks every model against a brute-force
transform of its image computed with ehtim's own Fourier matrix).

The workhorse is the **mG-ring** (Johnson et al. 2020; Sgr A* Paper IV eqs.
10-17): an infinitely thin ring of diameter d whose brightness varies around it
as sum_k beta_k exp(i k phi), blurred by a Gaussian of FWHM w, plus a Gaussian
"floor" of emission at its centre. Its transform is a sum of Bessel functions:

    V_ring(rho, phi_u) = sum_k beta_k J_k(pi d rho) exp(i k (phi_u + pi/2))

so a thin ring's amplitude falls to its first null where pi d rho = 2.405,
i.e. at rho = 0.766/d: 3.6 billion wavelengths for a 43-microarcsecond ring.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import j0, j1, jv

from config import RAD_PER_UAS

FWHM_TO_SIGMA = 1.0 / np.sqrt(8.0 * np.log(2.0))
ARCSEC_IN_UAS = 1e6


# ---------------------------------------------------------------------------
# Parameters and priors
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Param:
    """A model parameter with a uniform, log-uniform or periodic prior."""

    name: str
    lo: float
    hi: float
    kind: str = "uniform"  # "uniform" | "log" | "angle" (radians, periodic)
    label: str = ""

    def from_unit(self, x):
        if self.kind == "log":
            return np.exp(np.log(self.lo) + x * (np.log(self.hi) - np.log(self.lo)))
        return self.lo + x * (self.hi - self.lo)


def _bessel_stack(z: np.ndarray, kmax: int) -> np.ndarray:
    """J_0..J_kmax at z >= 0 (rows = order).

    Fast j0/j1, then the recurrence J_{k+1} = (2k/z) J_k - J_{k-1} (stable for
    k < z). For small z, where the recurrence loses digits, the leading series
    term (z/2)^k / k! is used instead (relative error < z^2/(4(k+1)) < 1e-4).
    20x faster than scipy's general jv; tests check it against jv.
    """
    z = np.asarray(z, dtype=float)
    out = np.empty((kmax + 1,) + z.shape)
    out[0] = j0(z)
    if kmax >= 1:
        out[1] = j1(z)
    if kmax >= 2:
        safe = np.where(z > 0, z, 1.0)
        small = z < 0.02
        fact = 1.0
        for k in range(1, kmax):
            out[k + 1] = (2.0 * k / safe) * out[k] - out[k - 1]
            fact *= k + 1
            out[k + 1] = np.where(small, (0.5 * z) ** (k + 1) / fact, out[k + 1])
    return out


# ---------------------------------------------------------------------------
# Building blocks (visibility of unit-flux shapes)
# ---------------------------------------------------------------------------
def gaussian_vis(rho, fwhm_uas):
    """Circular Gaussian of unit flux and FWHM `fwhm_uas`."""
    w = fwhm_uas * RAD_PER_UAS
    return np.exp(-((np.pi * w * rho) ** 2) / (4.0 * np.log(2.0)))


def elliptical_gaussian_vis(u, v, fwhm_maj_uas, fwhm_min_uas, pa_rad):
    """Elliptical Gaussian of unit flux; major axis at position angle pa (E of N)."""
    s, c = np.sin(pa_rad), np.cos(pa_rad)
    u_maj = u * s + v * c  # spatial frequency along the major axis
    u_min = u * c - v * s
    a = fwhm_maj_uas * RAD_PER_UAS
    b = fwhm_min_uas * RAD_PER_UAS
    return np.exp(-(np.pi**2) / (4.0 * np.log(2.0)) * ((a * u_maj) ** 2 + (b * u_min) ** 2))


def shift_phase(u, v, x0_uas, y0_uas):
    """Phase factor of a component displaced (x0 east, y0 north)."""
    return np.exp(2j * np.pi * (u * x0_uas + v * y0_uas) * RAD_PER_UAS)


def disk_vis(rho, d_uas):
    """Uniform disk of unit flux and diameter d: 2 J1(z)/z with z = pi d rho."""
    z = np.pi * d_uas * RAD_PER_UAS * rho
    z = np.where(z == 0, 1e-12, z)
    return 2.0 * j1(z) / z


def mring_vis(rho, phi_u, d_uas, betas):
    """Thin ring of unit flux, I(phi) proportional to 1 + 2 Re sum_k beta_k e^{ik phi}.

    `betas` are complex beta_1..beta_m (beta_0 = 1). phi_u = atan2(u, v).
    """
    z = np.pi * d_uas * RAD_PER_UAS * rho
    m = len(betas)
    J = _bessel_stack(z, m)
    out = J[0].astype(complex)
    for k, b in enumerate(betas, start=1):
        e = np.exp(1j * k * (phi_u + np.pi / 2.0))
        # beta_{-k} = conj(beta_k) and J_{-k} = (-1)^k J_k
        out = out + b * J[k] * e + np.conj(b) * ((-1) ** k) * J[k] * np.conj(e)
    return out


# ---------------------------------------------------------------------------
# Image-domain renderers (for figures and for the brute-force tests)
# ---------------------------------------------------------------------------
def render_mring(x, y, d_uas, betas, w_uas, n_az: int = 720):
    """Image of a blurred m-ring on a grid (x east, y north, uas), unit total flux."""
    phi = (np.arange(n_az) + 0.5) * 2 * np.pi / n_az
    weight = np.ones(n_az)
    for k, b in enumerate(betas, start=1):
        weight = weight + 2.0 * np.real(b * np.exp(1j * k * phi))
    weight = weight / n_az
    xr, yr = 0.5 * d_uas * np.sin(phi), 0.5 * d_uas * np.cos(phi)
    sig = max(w_uas, 1e-3) * FWHM_TO_SIGMA
    img = np.zeros_like(x, dtype=float)
    for wk, xk, yk in zip(weight, xr, yr):
        img += wk * np.exp(-((x - xk) ** 2 + (y - yk) ** 2) / (2 * sig**2))
    return img / (2 * np.pi * sig**2)


def render_gaussian(x, y, fwhm_uas, x0=0.0, y0=0.0):
    sig = fwhm_uas * FWHM_TO_SIGMA
    return np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / (2 * sig**2)) / (2 * np.pi * sig**2)


def render_elliptical_gaussian(x, y, fwhm_maj, fwhm_min, pa_rad, x0=0.0, y0=0.0):
    s, c = np.sin(pa_rad), np.cos(pa_rad)
    dx, dy = x - x0, y - y0
    along = dx * s + dy * c
    across = dx * c - dy * s
    sa, sb = fwhm_maj * FWHM_TO_SIGMA, fwhm_min * FWHM_TO_SIGMA
    return np.exp(-0.5 * ((along / sa) ** 2 + (across / sb) ** 2)) / (2 * np.pi * sa * sb)


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------
class Model:
    """Base class: a parameter list, a prior transform, visibilities and an image."""

    name = "model"
    description = ""

    def __init__(self, params: list[Param]):
        self.params = params
        self.names = [p.name for p in params]
        self.ndim = len(params)
        self._index = {n: i for i, n in enumerate(self.names)}

    def prior_transform(self, cube):
        return np.array([p.from_unit(c) for p, c in zip(self.params, cube)])

    def unpack(self, theta) -> dict:
        return {n: theta[i] for i, n in enumerate(self.names)}

    def vis(self, theta, u, v):  # pragma: no cover - abstract
        raise NotImplementedError

    def image(self, theta, x, y):  # pragma: no cover - abstract
        raise NotImplementedError

    # extra parameters every closure-only fit needs to explain the intra-site baselines
    @staticmethod
    def large_scale_params() -> list[Param]:
        return [
            Param("f_large", 0.0, 5.0, label="large-scale flux / compact flux"),
            Param("fwhm_large", 1e4, 1e7, "log", label="large-scale FWHM (uas)"),
        ]


def _ordered_gaussians_transform(model, theta, cube):
    """Nuisance Gaussians ordered in x: sorted uniforms (order statistics) for their x positions,
    which removes the label-switching degeneracy without changing the prior's shape."""
    n = getattr(model, "n_gauss", 0)
    if n > 1:
        idx = [model._index[f"g{j}_x"] for j in range(n)]
        uu = np.asarray(cube)[idx]
        xs = np.empty(n)
        xs[-1] = uu[-1] ** (1.0 / n)
        for k in range(n - 2, -1, -1):
            xs[k] = xs[k + 1] * uu[k] ** (1.0 / (k + 1))
        theta[idx] = -model.gauss_box + 2 * model.gauss_box * xs
    return theta


class MGRing(Model):
    """mG-ring of order m, optionally with the large-scale nuisance component.

    Parameters: d (diameter, uas); w (FWHM of the blurring kernel, uas); for each
    mode k, |beta_k| and arg beta_k; f_floor and w_floor (central Gaussian: flux
    relative to the ring, FWHM). The ring has unit flux: closure quantities do not
    know the absolute flux scale. The brightness position angle is eta = -arg beta_1.
    """

    def __init__(self, m: int = 1, *, floor: bool = True, large: bool = True, blur: bool = True,
                 d_range=(10.0, 100.0), w_max: float = 50.0, beta_max: float = 0.5,
                 total_flux: bool = False, n_gauss: int = 0, gauss_box: float = 100.0,
                 floor_flux_max: float = 1.0, floor_fwhm_range=(1.0, 100.0)):
        self.m, self.floor, self.large, self.blur, self.total_flux = m, floor, large, blur, total_flux
        self.n_gauss, self.gauss_box = n_gauss, gauss_box
        ps = [Param("d", *d_range, label="ring diameter d (uas)")]
        if blur:
            ps.append(Param("w", 0.0, w_max, label="ring FWHM w (uas)"))
        for k in range(1, m + 1):
            ps.append(Param(f"a{k}", 0.0, beta_max, label=f"|beta_{k}|"))
            ps.append(Param(f"p{k}", -np.pi, np.pi, "angle", label=f"arg beta_{k}"))
        if floor:
            ps.append(Param("f_floor", 0.0, floor_flux_max, label="floor flux / ring flux"))
            ps.append(Param("w_floor", *floor_fwhm_range, label="floor FWHM (uas)"))
        # "nuisance" Gaussians for emission the ring does not capture (Paper VI app. B.2);
        # ordered in x to remove the label-switching degeneracy
        for j in range(n_gauss):
            ps += [Param(f"g{j}_f", 0.0, 2.0, label=f"Gaussian {j + 1} flux / ring flux"),
                   Param(f"g{j}_x", -gauss_box, gauss_box, label=f"Gaussian {j + 1} x (uas, E)"),
                   Param(f"g{j}_y", -gauss_box, gauss_box, label=f"Gaussian {j + 1} y (uas, N)"),
                   Param(f"g{j}_w", 2.0, 100.0, label=f"Gaussian {j + 1} FWHM (uas)"),
                   Param(f"g{j}_q", 0.05, 1.0, label=f"Gaussian {j + 1} axis ratio"),
                   Param(f"g{j}_pa", 0.0, np.pi, label=f"Gaussian {j + 1} PA")]
        if large:
            ps.extend(self.large_scale_params())
        if total_flux:
            ps.append(Param("flux", 0.1, 2.0, label="compact flux (Jy)"))
        super().__init__(ps)
        self.name = (f"mG-ring m={m}" + ("" if floor else " (no floor)") + ("" if blur else " (thin)")
                     + (f" + {n_gauss} Gaussian" + ("s" if n_gauss > 1 else "") if n_gauss else ""))

    def prior_transform(self, cube):
        return _ordered_gaussians_transform(self, Model.prior_transform(self, cube), cube)

    def gaussians_vis(self, p, u, v):
        out = 0.0
        for j in range(self.n_gauss):
            w = p[f"g{j}_w"]
            out = out + p[f"g{j}_f"] * elliptical_gaussian_vis(u, v, w, w * p[f"g{j}_q"], p[f"g{j}_pa"]) * shift_phase(
                u, v, p[f"g{j}_x"], p[f"g{j}_y"])
        return out

    def betas(self, p: dict):
        return [p[f"a{k}"] * np.exp(1j * p[f"p{k}"]) for k in range(1, self.m + 1)]

    def compact_vis(self, theta, u, v):
        p = self.unpack(theta)
        rho = np.hypot(u, v)
        phi_u = np.arctan2(u, v)
        out = mring_vis(rho, phi_u, p["d"], self.betas(p))
        if self.blur:
            out = out * gaussian_vis(rho, p["w"])
        if self.floor:
            out = out + p["f_floor"] * gaussian_vis(rho, p["w_floor"])
        if self.n_gauss:
            out = out + self.gaussians_vis(p, u, v)
        if self.total_flux:
            out = out * p["flux"] / self.compact_flux_units(p)
        return out

    def compact_flux_units(self, p) -> float:
        """Total compact flux in units of the ring's flux."""
        return 1.0 + (p["f_floor"] if self.floor else 0.0) + sum(p[f"g{j}_f"] for j in range(self.n_gauss))

    def vis(self, theta, u, v):
        out = self.compact_vis(theta, u, v)
        if self.large:
            p = self.unpack(theta)
            out = out + p["f_large"] * gaussian_vis(np.hypot(u, v), p["fwhm_large"])
        return out

    def image(self, theta, x, y, *, min_blur: float = 2.0):
        p = self.unpack(theta)
        w = p["w"] if self.blur else 0.0
        img = render_mring(x, y, p["d"], self.betas(p), max(w, min_blur))
        if self.floor:
            img = img + p["f_floor"] * render_gaussian(x, y, p["w_floor"])
        for j in range(self.n_gauss):
            img = img + p[f"g{j}_f"] * render_elliptical_gaussian(
                x, y, p[f"g{j}_w"], p[f"g{j}_w"] * p[f"g{j}_q"], p[f"g{j}_pa"], p[f"g{j}_x"], p[f"g{j}_y"])
        return img

    def derived(self, theta) -> dict:
        p = self.unpack(theta)
        out = {"d": p["d"]}
        if self.blur:
            out["w"] = p["w"]
            out["fractional_width"] = p["w"] / p["d"]
        if self.m >= 1:
            out["asymmetry"] = p["a1"]  # |beta_1|: 0 = uniform ring, 0.5 = all on one side
            out["pa_deg"] = np.degrees(-p["p1"]) % 360.0  # brightest direction, E of N
        total = self.compact_flux_units(p)
        out["ring_fraction"] = 1.0 / total  # share of the compact flux in the ring itself
        if self.floor:
            out["floor_fraction"] = p["f_floor"] / total
        if self.n_gauss:
            out["gauss_fraction"] = sum(p[f"g{j}_f"] for j in range(self.n_gauss)) / total
        return out


class SymmetricShape(Model):
    """A point-symmetric shape: 'gaussian', 'disk' or 'ring' (thin), all centred.

    Point-symmetric images have closure phases of exactly 0 or 180 degrees, so the
    M87 closure phases alone can reject them.
    """

    def __init__(self, shape: str, *, large: bool = True, size_range=(5.0, 150.0)):
        self.shape, self.large = shape, large
        ps = [Param("size", *size_range, label={"gaussian": "FWHM (uas)", "disk": "diameter (uas)",
                                                 "ring": "diameter (uas)"}[shape])]
        if large:
            ps.extend(self.large_scale_params())
        super().__init__(ps)
        self.name = {"gaussian": "Gaussian", "disk": "uniform disk", "ring": "thin ring"}[shape]

    def vis(self, theta, u, v):
        p = self.unpack(theta)
        rho = np.hypot(u, v)
        if self.shape == "gaussian":
            out = gaussian_vis(rho, p["size"]).astype(complex)
        elif self.shape == "disk":
            out = disk_vis(rho, p["size"]).astype(complex)
        else:
            out = j0(np.pi * p["size"] * RAD_PER_UAS * rho).astype(complex)
        if self.large:
            out = out + p["f_large"] * gaussian_vis(rho, p["fwhm_large"])
        return out

    def image(self, theta, x, y):
        p = self.unpack(theta)
        r = np.hypot(x, y)
        if self.shape == "gaussian":
            return render_gaussian(x, y, p["size"])
        if self.shape == "disk":
            img = (r <= p["size"] / 2).astype(float)
            return img / max(img.sum(), 1)
        return render_mring(x, y, p["size"], [], 2.0)


class Crescent(Model):
    """Top-hat crescent (Kamruddin & Dexter 2013): a uniform disk minus a smaller,
    offset uniform disk, optionally blurred.

    Parameters: r_out (uas), psi = 1 - r_in/r_out (fractional thickness), tau
    (0 = inner disk touching the outer edge, 1 = concentric), pa (direction of the
    thickest side, E of N), sigma (Gaussian blur, standard deviation, uas).
    The 'diameter' is d = r_out + r_in, twice the mean radius, as in Paper VI.
    """

    def __init__(self, *, large: bool = True, blur: bool = True):
        self.large, self.blur = large, blur
        ps = [
            Param("r_out", 5.0, 60.0, label="outer radius (uas)"),
            Param("psi", 0.0, 1.0, label="fractional thickness"),
            Param("tau", 0.0, 1.0, label="symmetry"),
            Param("pa", -np.pi, np.pi, "angle", label="PA of thick side"),
        ]
        if blur:
            ps.append(Param("sigma", 0.0, 20.0, label="blur sigma (uas)"))
        if large:
            ps.extend(self.large_scale_params())
        super().__init__(ps)
        self.name = "crescent"

    def geometry(self, p):
        r_out = p["r_out"]
        r_in = r_out * (1.0 - p["psi"])
        offset = p["psi"] * r_out * (1.0 - p["tau"])  # inner disk moved away from the thick side
        return r_out, r_in, offset

    def vis(self, theta, u, v):
        p = self.unpack(theta)
        r_out, r_in, offset = self.geometry(p)
        rho = np.hypot(u, v)
        area_out, area_in = np.pi * r_out**2, np.pi * r_in**2
        # inner disk centred at distance `offset` opposite to the thick side
        x0, y0 = -offset * np.sin(p["pa"]), -offset * np.cos(p["pa"])
        out = area_out * disk_vis(rho, 2 * r_out) - area_in * disk_vis(rho, 2 * r_in) * shift_phase(u, v, x0, y0)
        out = out / (area_out - area_in)
        if self.blur:
            out = out * gaussian_vis(rho, p["sigma"] / FWHM_TO_SIGMA)
        if self.large:
            out = out + p["f_large"] * gaussian_vis(rho, p["fwhm_large"])
        return out

    def image(self, theta, x, y):
        p = self.unpack(theta)
        r_out, r_in, offset = self.geometry(p)
        x0, y0 = -offset * np.sin(p["pa"]), -offset * np.cos(p["pa"])
        img = ((np.hypot(x, y) <= r_out) & (np.hypot(x - x0, y - y0) > r_in)).astype(float)
        if self.blur and p["sigma"] > 0:
            from scipy.ndimage import gaussian_filter
            pix = abs(x[0, 1] - x[0, 0])
            img = gaussian_filter(img, p["sigma"] / pix)
        return img / max(img.sum(), 1e-30)

    def derived(self, theta) -> dict:
        p = self.unpack(theta)
        r_out, r_in, offset = self.geometry(p)
        return {"d": r_out + r_in, "fractional_width": (r_out - r_in) / (0.5 * (r_out + r_in)),
                "pa_deg": np.degrees(p["pa"]) % 360.0}


class GaussianSum(Model):
    """Two elliptical Gaussians (a comparably flexible model with no ring)."""

    def __init__(self, *, large: bool = True):
        self.large = large
        ps = [
            Param("w1", 2.0, 120.0, label="FWHM major 1"), Param("q1", 0.05, 1.0, label="axis ratio 1"),
            Param("pa1", 0.0, np.pi, label="PA 1"),
            Param("f2", 0.0, 1.0, label="flux 2 / flux 1"), Param("x2", -80.0, 80.0, label="x2"),
            Param("y2", -80.0, 80.0, label="y2"), Param("w2", 2.0, 120.0, label="FWHM major 2"),
            Param("q2", 0.05, 1.0, label="axis ratio 2"), Param("pa2", 0.0, np.pi, label="PA 2"),
        ]
        if large:
            ps.extend(self.large_scale_params())
        super().__init__(ps)
        self.name = "two Gaussians"

    def vis(self, theta, u, v):
        p = self.unpack(theta)
        out = elliptical_gaussian_vis(u, v, p["w1"], p["w1"] * p["q1"], p["pa1"]).astype(complex)
        out = out + p["f2"] * elliptical_gaussian_vis(u, v, p["w2"], p["w2"] * p["q2"], p["pa2"]) * shift_phase(
            u, v, p["x2"], p["y2"])
        if self.large:
            out = out + p["f_large"] * gaussian_vis(np.hypot(u, v), p["fwhm_large"])
        return out

    def image(self, theta, x, y):
        p = self.unpack(theta)
        return (render_elliptical_gaussian(x, y, p["w1"], p["w1"] * p["q1"], p["pa1"])
                + p["f2"] * render_elliptical_gaussian(x, y, p["w2"], p["w2"] * p["q2"], p["pa2"], p["x2"], p["y2"]))


def _disk_ft_and_slope(rho, k_xi, R_uas):
    """Unit-brightness disk of radius R: its transform D(rho) and F[xi * disk].

    D(rho) = R J1(2 pi rho R) / rho, and multiplying the image by xi (the
    coordinate along the crescent's axis) takes a derivative in the Fourier
    domain: F[xi f] = (1 / 2 pi i) dF/dk_xi (this transform's + sign convention).
    Lengths in radians, rho in wavelengths.
    """
    R = R_uas * RAD_PER_UAS
    x = 2 * np.pi * rho * R
    J = _bessel_stack(x, 2)
    D = R * J[1] / rho
    dD_drho = -(R / rho**2) * J[1] + (np.pi * R**2 / rho) * (J[0] - J[2])
    xiD = (k_xi / rho) * dD_drho / (2j * np.pi)
    return D, xiD


class XSRing(Model):
    """The 'xs-ring' generalised crescent of M87 Paper VI (appendix B.1.1), re-derived here.

    A uniform disk of radius r_out minus a smaller disk of radius
    r_in = r_out (1 - psi) whose centre is moved r0 = psi r_out (1 - tau) towards the
    thin side; a linear brightness gradient ('slash') along the axis, from 1 at the
    thick edge to beta at the thin edge; a uniform 'floor' of brightness gamma * beta
    filling the hole; the whole blurred by a Gaussian of standard deviation sigma. The
    thick side points to position angle `pa`. Paper VI's crescent diameter is
    d = r_out + r_in (twice the mean radius), and its GRMHD calibration
    alpha = d / theta_g = 11.56 applies to exactly this model and definition.

    Fitted, as in Paper VI, with nuisance Gaussians and a large-scale component.
    """

    def __init__(self, *, n_gauss: int = 3, large: bool = True, gauss_box: float = 100.0):
        self.n_gauss, self.large, self.gauss_box = n_gauss, large, gauss_box
        ps = [
            Param("r_out", 5.0, 60.0, label="outer radius (uas)"),
            Param("psi", 0.0, 1.0, label="fractional thickness psi"),
            Param("tau", 0.0, 1.0, label="symmetry tau"),
            Param("pa", -np.pi, np.pi, "angle", label="PA of the thick side"),
            Param("sigma", 0.0, 30.0, label="blur sigma (uas)"),
            Param("gamma", 0.0, 1.0, label="floor / thin-edge brightness"),
            Param("beta", np.exp(-5.0), np.exp(5.0), "log", label="thin / thick edge brightness"),
        ]
        for j in range(n_gauss):
            ps += [Param(f"g{j}_f", 0.0, 2.0), Param(f"g{j}_x", -gauss_box, gauss_box),
                   Param(f"g{j}_y", -gauss_box, gauss_box), Param(f"g{j}_w", 2.0, 100.0),
                   Param(f"g{j}_q", 0.05, 1.0), Param(f"g{j}_pa", 0.0, np.pi)]
        if large:
            ps.extend(self.large_scale_params())
        super().__init__(ps)
        self.name = "xs-ring crescent" + (f" + {n_gauss} Gaussians" if n_gauss else "")

    def prior_transform(self, cube):
        return _ordered_gaussians_transform(self, Model.prior_transform(self, cube), cube)

    gaussians_vis = MGRing.gaussians_vis  # a plain function of (self, p, u, v): safe to share

    def geometry(self, p):
        r_out = p["r_out"]
        r_in = r_out * (1.0 - p["psi"])
        r0 = p["psi"] * r_out * (1.0 - p["tau"])
        return r_out, r_in, r0

    def crescent_vis(self, p, u, v):
        """Unit-flux slashed crescent with floor, blurred (no nuisance components)."""
        r_out, r_in, r0 = self.geometry(p)
        beta, gamma = p["beta"], p["gamma"]
        rho = np.hypot(u, v)
        rho = np.where(rho == 0, 1e-3, rho)
        k_xi = u * np.sin(p["pa"]) + v * np.cos(p["pa"])  # spatial frequency along the axis
        slope = (1.0 - beta) / (2.0 * r_out * RAD_PER_UAS)  # brightness gradient per radian
        xi_c = -r0 * RAD_PER_UAS  # inner disk centre, towards the thin side
        s_c = 0.5 * (1.0 + beta) + slope * xi_c  # slash brightness at the inner disk centre
        D_out, xiD_out = _disk_ft_and_slope(rho, k_xi, r_out)
        D_in, xiD_in = _disk_ft_and_slope(rho, k_xi, max(r_in, 1e-6))
        shift = np.exp(2j * np.pi * k_xi * xi_c)
        K = gamma * beta
        V = 0.5 * (1.0 + beta) * D_out + slope * xiD_out - shift * ((s_c - K) * D_in + slope * xiD_in)
        R_out, R_in = r_out * RAD_PER_UAS, r_in * RAD_PER_UAS
        flux = np.pi * (0.5 * (1.0 + beta) * R_out**2 - (s_c - K) * R_in**2)
        V = V / flux
        sig = p["sigma"] * RAD_PER_UAS
        return V * np.exp(-2.0 * np.pi**2 * sig**2 * rho**2)

    def compact_vis(self, theta, u, v):
        p = self.unpack(theta)
        out = self.crescent_vis(p, u, v)
        if self.n_gauss:
            out = out + self.gaussians_vis(p, u, v)
        return out

    def vis(self, theta, u, v):
        out = self.compact_vis(theta, u, v)
        if self.large:
            p = self.unpack(theta)
            out = out + p["f_large"] * gaussian_vis(np.hypot(u, v), p["fwhm_large"])
        return out

    def crescent_image(self, p, x, y):
        r_out, r_in, r0 = self.geometry(p)
        s_, c_ = np.sin(p["pa"]), np.cos(p["pa"])
        xi = x * s_ + y * c_  # coordinate along the axis (towards the thick side)
        zeta = x * c_ - y * s_
        slash = 0.5 * (1 + p["beta"]) + (1 - p["beta"]) / (2 * r_out) * xi
        inner = np.hypot(xi + r0, zeta) <= r_in
        img = np.where((np.hypot(xi, zeta) <= r_out) & ~inner, slash, 0.0)
        img = img + np.where(inner, p["gamma"] * p["beta"], 0.0)
        if p["sigma"] > 0:
            from scipy.ndimage import gaussian_filter
            img = gaussian_filter(img, p["sigma"] / abs(x[0, 1] - x[0, 0]))
        return img / img.sum()

    def image(self, theta, x, y):
        p = self.unpack(theta)
        img = self.crescent_image(p, x, y)
        for j in range(self.n_gauss):
            img = img + p[f"g{j}_f"] * render_elliptical_gaussian(
                x, y, p[f"g{j}_w"], p[f"g{j}_w"] * p[f"g{j}_q"], p[f"g{j}_pa"], p[f"g{j}_x"], p[f"g{j}_y"]) * (
                abs(x[0, 1] - x[0, 0]) ** 2)
        return img

    def brightness_offset(self, p) -> float:
        """Centre of light of the crescent along its axis (uas; + towards the thick side).

        First moment of the slashed outer disk minus the slashed inner disk, plus the
        floor: the slash makes one side brighter, so the bright side need not be the
        thick side.
        """
        r_out, r_in, r0 = self.geometry(p)
        beta, K = p["beta"], p["gamma"] * p["beta"]
        slope = (1.0 - beta) / (2.0 * r_out)
        xi_c = -r0
        s_c = 0.5 * (1.0 + beta) + slope * xi_c
        flux = np.pi * (0.5 * (1.0 + beta) * r_out**2 - (s_c - K) * r_in**2)
        m1 = slope * np.pi * r_out**4 / 4.0 - ((s_c - K) * xi_c * np.pi * r_in**2 + slope * np.pi * r_in**4 / 4.0)
        return m1 / flux

    def derived(self, theta) -> dict:
        p = self.unpack(theta)
        r_out, r_in, r0 = self.geometry(p)
        total = 1.0 + sum(p[f"g{j}_f"] for j in range(self.n_gauss))
        bright = np.degrees(p["pa"]) + (0.0 if self.brightness_offset(p) >= 0 else 180.0)
        return {"d": r_out + r_in, "pa_thick_deg": np.degrees(p["pa"]) % 360.0, "pa_deg": bright % 360.0,
                "psi": p["psi"], "sharpness": p["sigma"] / (r_out + r_in), "ring_fraction": 1.0 / total}
