"""Roche geometry of close binaries, and light curves computed from it.

Units: the separation is 1, the more massive star ("primary", mass M1) sits at
the origin and the companion (mass M2 = q M1, q <= 1) at (1, 0, 0). The frame
co-rotates with the circular, synchronous orbit; z is the rotation axis.

The dimensionless (Kopal) potential, larger deeper inside the stars, is

    Omega(x, y, z) = 1/r1 + q (1/r2 - x) + (1 + q)/2 (x^2 + y^2),

with r1, r2 the distances to the two stars. It differs from the true potential
Phi = -G M1/a [Omega - q^2/(2(1 + q))] only by a scale and an offset, so its
equipotentials are the surfaces of constant Phi and |grad Omega| is the local
effective gravity in units of G M1 / a^2.

* Lagrange points: grad Omega = 0 on the x axis; L1 between the stars, L2
  behind the companion, L3 behind the primary.
* A contact binary's common surface is the equipotential
  Omega_f = Omega(L1) - f [Omega(L1) - Omega(L2)], f the usual fill-out
  factor: f = 0 just touching at L1, f = 1 at the outer
  critical surface, beyond which gas leaves through L2.
* ``render`` ray-traces the surface as seen from any direction, with von
  Zeipel/Lucy gravity darkening, a limb-darkening law and Planck intensities,
  and sums the image into a flux. ``light_curve`` does that around the orbit.
  This is a deliberately small model (no reflection effect, no spots unless
  asked for, blackbody intensities): enough to show *why* contact binaries
  look the way they do, and to fit their mean light curves at the 1% level.
"""

from __future__ import annotations

import math

import numba as nb
import numpy as np
from scipy.optimize import brentq

# ---------------------------------------------------------------------------
# Potential and Lagrange points
# ---------------------------------------------------------------------------


def omega(x, y, z, q):
    """Kopal's dimensionless Roche potential (array-friendly)."""
    r1 = np.sqrt(x * x + y * y + z * z)
    r2 = np.sqrt((x - 1.0) ** 2 + y * y + z * z)
    return 1.0 / r1 + q * (1.0 / r2 - x) + 0.5 * (1.0 + q) * (x * x + y * y)


def grad_omega(x, y, z, q):
    r1 = np.sqrt(x * x + y * y + z * z)
    r2 = np.sqrt((x - 1.0) ** 2 + y * y + z * z)
    gx = -x / r1**3 - q * (x - 1.0) / r2**3 - q + (1.0 + q) * x
    gy = -y / r1**3 - q * y / r2**3 + (1.0 + q) * y
    gz = -z / r1**3 - q * z / r2**3
    return gx, gy, gz


def _domega_dx_axis(x, q):
    return -np.sign(x) / x**2 - q * np.sign(x - 1.0) / (x - 1.0) ** 2 - q + (1.0 + q) * x


def lagrange_points(q: float) -> dict:
    """x positions and potentials of L1, L2 (behind the companion) and L3."""
    eps = 1e-9
    x1 = brentq(_domega_dx_axis, eps, 1.0 - eps, args=(q,), xtol=1e-14)
    x2 = brentq(_domega_dx_axis, 1.0 + eps, 3.0, args=(q,), xtol=1e-14)
    x3 = brentq(_domega_dx_axis, -3.0, -eps, args=(q,), xtol=1e-14)
    out = {}
    for name, x in (("L1", x1), ("L2", x2), ("L3", x3)):
        out[name] = x
        out[f"omega_{name}"] = float(omega(x, 0.0, 0.0, q))
    # L4/L5 are at the apexes of equilateral triangles in the orbital plane.
    out["L4"] = (0.5, math.sqrt(3) / 2)
    return out


def fillout_potential(q: float, f: float) -> float:
    """The potential of a common envelope with fill-out factor ``f``.

    f < 0 describes detached stars inside their lobes (then the formula is
    just a convenient label); 0 <= f <= 1 is a contact binary.
    """
    lp = lagrange_points(q)
    return lp["omega_L1"] - f * (lp["omega_L1"] - lp["omega_L2"])


def fillout_factor(q: float, omega_surf: float) -> float:
    lp = lagrange_points(q)
    return (lp["omega_L1"] - omega_surf) / (lp["omega_L1"] - lp["omega_L2"])


def eggleton_radius(q: float) -> float:
    """Eggleton (1983) volume-equivalent Roche-lobe radius of a star of mass ratio
    q = M_star / M_other, in units of the separation (accurate to ~1%)."""
    q23 = q ** (2.0 / 3.0)
    return 0.49 * q23 / (0.6 * q23 + math.log(1.0 + q ** (1.0 / 3.0)))


# ---------------------------------------------------------------------------
# Surfaces: radii along rays, lobe volumes, meshes
# ---------------------------------------------------------------------------


def _surface_radius(q, omega_surf, cx, ux, uy, uz, rmax, dr=0.004):
    """Distance from the star centre (cx, 0, 0) along the unit vector u to the
    first crossing of the equipotential ``omega_surf``, or ``rmax`` if the ray
    reaches it first (e.g. the plane through L1).

    The search steps outwards and stops at the first exit: far out in the
    orbital plane Omega rises again (centrifugal term), so a bisection over
    the whole ray could lock onto that outer region instead.
    """
    lo = 1e-6
    r = dr
    while r < rmax:
        if omega(cx + r * ux, r * uy, r * uz, q) < omega_surf:
            break
        lo = r
        r += dr
    else:
        return rmax
    hi = min(r, rmax)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if omega(cx + mid * ux, mid * uy, mid * uz, q) - omega_surf > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def lobe_volume_radius(q: float, which: str = "primary", n_mu: int = 96, n_phi: int = 96) -> float:
    """Volume-equivalent radius of the critical (L1) lobe of one star, by direct
    integration of r^3/3 over the sphere of directions around its centre.

    The lobe is cut by the plane x = x_L1 through the inner Lagrange point; for
    the critical equipotential that plane is where the two lobes meet.
    """
    lp = lagrange_points(q)
    om = lp["omega_L1"]
    xl1 = lp["L1"]
    cx = 0.0 if which == "primary" else 1.0
    mu, w_mu = np.polynomial.legendre.leggauss(n_mu)
    phis = (np.arange(n_phi) + 0.5) * 2 * np.pi / n_phi
    vol = 0.0
    for m, wm in zip(mu, w_mu):
        s = math.sqrt(max(0.0, 1.0 - m * m))
        for ph in phis:
            ux, uy, uz = m, s * math.cos(ph), s * math.sin(ph)  # polar axis along x
            # distance to the cutting plane along this ray (only rays heading towards L1)
            if which == "primary":
                rplane = (xl1 - cx) / ux if ux > 1e-12 else 10.0
            else:
                rplane = (xl1 - cx) / ux if ux < -1e-12 else 10.0
            rmax = min(rplane, 2.0)
            r = _surface_radius(q, om, cx, ux, uy, uz, rmax)
            vol += wm * (2 * np.pi / n_phi) * r**3 / 3.0
    return (3.0 * vol / (4.0 * np.pi)) ** (1.0 / 3.0)


def _stellar_region(q: float, omega_surf: float, n: int = 120, pad: float = 0.06):
    """Grid of Omega and a mask of the voxels inside the stars' surface(s).

    Omega >= omega_surf also holds far out in the orbital plane, where the
    centrifugal term dominates; only the connected regions that contain a
    star's centre are kept.
    """
    from scipy import ndimage
    lp = lagrange_points(q)
    x0, x1 = lp["L3"] - pad, lp["L2"] + pad
    ymax = 0.85
    xs = np.linspace(x0, x1, n)
    dx = xs[1] - xs[0]
    ys = np.arange(-ymax, ymax + 0.5 * dx, dx)
    zs = ys.copy()
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    with np.errstate(divide="ignore", invalid="ignore"):
        F = omega(X, Y, Z, q)
    F = np.nan_to_num(F, nan=1e6, posinf=1e6)
    labels, _ = ndimage.label(F >= omega_surf)
    keep = set()
    for cxs in (0.0, 1.0):
        idx = (np.abs(xs - cxs).argmin(), np.abs(ys).argmin(), np.abs(zs).argmin())
        if labels[idx] > 0:
            keep.add(labels[idx])
    stars = np.isin(labels, list(keep))
    return xs, ys, zs, F, stars


def surface_mesh(q: float, omega_surf: float, n: int = 120):
    """Triangle mesh of the stars' equipotential surface (marching cubes).

    Returns vertices (N, 3) and faces (M, 3).
    """
    from skimage.measure import marching_cubes
    xs, ys, zs, F, stars = _stellar_region(q, omega_surf, n=n)
    # Voxels outside the stars' own regions are pushed below the level, so the
    # outer (centrifugal) region produces no surface.
    G = np.where(stars, np.minimum(F, 1e3), np.minimum(F, omega_surf - 1.0))
    G[~stars & (F >= omega_surf)] = omega_surf - 1.0
    verts, faces, _, _ = marching_cubes(G, level=omega_surf,
                                        spacing=(xs[1] - xs[0], ys[1] - ys[0], zs[1] - zs[0]))
    verts += np.array([xs[0], ys[0], zs[0]])
    return verts, faces


def bounding_sphere(q: float, omega_surf: float, n: int = 100):
    """Centre (on the x axis) and radius of a sphere that encloses the stars'
    surface and no other region with Omega >= omega_surf."""
    xs, ys, zs, F, stars = _stellar_region(q, omega_surf, n=n)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    ix = np.nonzero(stars.any(axis=(1, 2)))[0]
    cx = 0.5 * (xs[ix[0]] + xs[ix[-1]])
    d = np.sqrt((X - cx) ** 2 + Y**2 + Z**2)
    rb = d[stars].max() + 2.0 * (xs[1] - xs[0])
    intruders = (~stars) & (F >= omega_surf) & (d < rb)
    if intruders.any():
        raise ValueError(f"q={q}: the outer centrifugal region reaches inside the bounding sphere "
                         f"(r = {rb:.3f}); surface potential too low (fill-out > 1?)")
    return float(cx), float(rb)


def orbital_plane_contours(q: float, levels, extent=(-1.6, 2.6, -1.6, 1.6), n: int = 800):
    """Omega on a grid in the orbital plane, for drawing equipotentials."""
    xs = np.linspace(extent[0], extent[1], n)
    ys = np.linspace(extent[2], extent[3], int(n * (extent[3] - extent[2]) / (extent[1] - extent[0])))
    X, Y = np.meshgrid(xs, ys)
    with np.errstate(divide="ignore"):
        F = omega(X, Y, 0.0 * X, q)
    return X, Y, F


# ---------------------------------------------------------------------------
# Ray-traced images and light curves
# ---------------------------------------------------------------------------

#: Lucy (1967) gravity-darkening exponent for convective envelopes: T ~ g^(beta/4).
BETA_CONVECTIVE = 0.32
#: Planck constant times c over k_B, in m K.
_HC_K = 1.438776877e-2


@nb.njit(cache=True, fastmath=True)
def _omega_nb(x, y, z, q):
    r1 = math.sqrt(x * x + y * y + z * z)
    r2 = math.sqrt((x - 1.0) * (x - 1.0) + y * y + z * z)
    return 1.0 / r1 + q * (1.0 / r2 - x) + 0.5 * (1.0 + q) * (x * x + y * y)


@nb.njit(cache=True, fastmath=True)
def _grad_nb(x, y, z, q):
    r1 = math.sqrt(x * x + y * y + z * z)
    r2 = math.sqrt((x - 1.0) * (x - 1.0) + y * y + z * z)
    r13 = r1 * r1 * r1
    r23 = r2 * r2 * r2
    gx = -x / r13 - q * (x - 1.0) / r23 - q + (1.0 + q) * x
    gy = -y / r13 - q * y / r23 + (1.0 + q) * y
    gz = -z / r13 - q * z / r23
    return gx, gy, gz


@nb.njit(cache=True, parallel=True, fastmath=True)
def _trace(q, om_surf, sx, sy, sz, e1x, e1y, e1z, e2x, e2y, e2z, cx, us, vs, rb, step,
           x_split, t1, t2, g_ref1, g_ref2, beta, lam_m, ld_a, ld_b,
           spot_lon, spot_lat, spot_rad, spot_fac, spot_star):
    """Intensity image: one ray per (u, v) pixel, cast from the observer side.

    Each ray is clipped to the sphere of radius ``rb`` around (cx, 0, 0), which
    encloses the stars' surface and nothing else with Omega >= om_surf (far out
    in the orbital plane the centrifugal term makes Omega large again).
    Returns the image (relative surface brightness per pixel) and a hit mask.
    Temperatures: T = T_star (g / g_ref_star)^(beta/4) on each side of the
    plane x = x_split (the neck), so the two stars may differ in temperature.
    Limb darkening: I(mu)/I(1) = 1 - a (1 - mu) - b (1 - mu)^2.
    One circular spot (angular radius spot_rad, temperature factor spot_fac)
    may be placed on a star (1 or 2) at (spot_lon, spot_lat): longitude 0 faces
    the companion, +90 deg points along +y (the leading side for phase 0.25).
    """
    nu, nv = us.size, vs.size
    img = np.zeros((nv, nu))
    hit = np.zeros((nv, nu), dtype=np.bool_)
    c2 = _HC_K / lam_m
    for j in nb.prange(nv):
        for i in range(nu):
            d2 = us[i] * us[i] + vs[j] * vs[j]
            if d2 >= rb * rb:
                continue
            h = math.sqrt(rb * rb - d2)
            # entry point on the observer side of the sphere; march along -s for 2h
            ox = cx + us[i] * e1x + vs[j] * e2x + h * sx
            oy = us[i] * e1y + vs[j] * e2y + h * sy
            oz = us[i] * e1z + vs[j] * e2z + h * sz
            nsteps = int(2.0 * h / step) + 2
            t_prev = 0.0
            found = False
            t = 0.0
            for k in range(1, nsteps):
                t = min(k * step, 2.0 * h)
                if _omega_nb(ox - t * sx, oy - t * sy, oz - t * sz, q) >= om_surf:
                    lo, hi = t_prev, t
                    for _ in range(30):
                        mid = 0.5 * (lo + hi)
                        if _omega_nb(ox - mid * sx, oy - mid * sy, oz - mid * sz, q) >= om_surf:
                            hi = mid
                        else:
                            lo = mid
                    t = 0.5 * (lo + hi)
                    found = True
                    break
                t_prev = t
            if not found:
                continue
            px = ox - t * sx
            py = oy - t * sy
            pz = oz - t * sz
            gx, gy, gz = _grad_nb(px, py, pz, q)
            g = math.sqrt(gx * gx + gy * gy + gz * gz)
            # outward normal is -grad Omega (Omega decreases outwards)
            mu = -(gx * sx + gy * sy + gz * sz) / g
            if mu <= 0.0:
                mu = 1e-6
            if px < x_split:
                temp = t1 * (g / g_ref1) ** (0.25 * beta)
                star = 1
                cxs = 0.0
            else:
                temp = t2 * (g / g_ref2) ** (0.25 * beta)
                star = 2
                cxs = 1.0
            if spot_star == star and spot_rad > 0.0:
                dx = px - cxs
                rr = math.sqrt(dx * dx + py * py + pz * pz)
                ux, uy, uz = dx / rr, py / rr, pz / rr
                sgn = 1.0 if star == 1 else -1.0
                cl = math.cos(spot_lat)
                vx = sgn * cl * math.cos(spot_lon)
                vy = cl * math.sin(spot_lon)
                vz = math.sin(spot_lat)
                if ux * vx + uy * vy + uz * vz >= math.cos(spot_rad):
                    temp *= spot_fac
            planck = 1.0 / (math.exp(c2 / temp) - 1.0)
            ld = 1.0 - ld_a * (1.0 - mu) - ld_b * (1.0 - mu) * (1.0 - mu)
            img[j, i] = planck * ld
            hit[j, i] = True
    return img, hit


def view_vectors(incl_deg: float, phase: float):
    """Unit vector towards the observer and the two sky-plane axes, in the
    co-rotating frame. Phase 0: the companion is in front of the primary
    (the primary is eclipsed); the system turns anticlockwise for the observer."""
    i = math.radians(incl_deg)
    ph = 2 * math.pi * phase
    s = np.array([math.sin(i) * math.cos(ph), -math.sin(i) * math.sin(ph), math.cos(i)])
    zaxis = np.array([0.0, 0.0, 1.0])
    e2 = zaxis - s * (zaxis @ s)
    n2 = np.linalg.norm(e2)
    e2 = e2 / n2 if n2 > 1e-12 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(e2, s)
    return s, e1, e2


class ContactBinary:
    """A Roche-geometry binary that can be imaged and turned into a light curve.

    Parameters: mass ratio q (<= 1), fill-out factor f (contact for 0..1; a
    negative f gives two detached stars whose surfaces are the corresponding
    equipotential, which is only an approximation for the smaller star),
    temperatures of the two components (K), the effective wavelength of the
    bandpass (nm), gravity-darkening exponent and quadratic limb darkening.
    """

    def __init__(self, q, f, t1=5800.0, t2=None, lam_nm=640.0, beta=BETA_CONVECTIVE,
                 ld=(0.40, 0.25), spot=None):
        self.q = float(q)
        self.f = float(f)
        self.lp = lagrange_points(self.q)
        self.om = fillout_potential(self.q, self.f)
        self.t1 = float(t1)
        self.t2 = float(t2 if t2 is not None else t1)
        self.lam_m = lam_nm * 1e-9
        self.beta = float(beta)
        self.ld = ld
        # spot: (star 1|2, longitude deg, latitude deg, radius deg, temperature factor)
        self.spot = spot
        # reference gravities: at each star's pole, where von Zeipel's law is anchored
        self.g_ref1 = self._pole_gravity(0.0)
        self.g_ref2 = self._pole_gravity(1.0)
        # a sphere around the stars' surface, inside which every ray is marched
        self.cx, self.rb = bounding_sphere(self.q, self.om)

    def _pole_gravity(self, cx):
        r = _surface_radius(self.q, self.om, cx, 0.0, 0.0, 1.0, 1.5)
        gx, gy, gz = grad_omega(cx, 0.0, r, self.q)
        return float(np.sqrt(gx * gx + gy * gy + gz * gz))

    def render(self, incl_deg, phase, npix=160, step=0.004):
        """Image of the binary on the sky (rows = sky 'north' along the
        projected rotation axis) plus the pixel area, for summing into a flux."""
        s, e1, e2 = view_vectors(incl_deg, phase)
        half = self.rb
        us = np.linspace(-half, half, npix)
        vs = np.linspace(-half, half, npix)
        if self.spot is None:
            sp = (0, 0.0, 0.0, 0.0, 1.0)
        else:
            star, lon, lat, rad, fac = self.spot
            sp = (int(star), math.radians(lon), math.radians(lat), math.radians(rad), float(fac))
        img, hit = _trace(self.q, self.om, s[0], s[1], s[2], e1[0], e1[1], e1[2], e2[0], e2[1], e2[2],
                          self.cx, us, vs, self.rb, step, self.lp["L1"], self.t1, self.t2,
                          self.g_ref1, self.g_ref2, self.beta, self.lam_m, self.ld[0], self.ld[1],
                          sp[1], sp[2], sp[3], sp[4], sp[0])
        pix_area = (us[1] - us[0]) * (vs[1] - vs[0])
        return img, hit, pix_area

    def flux(self, incl_deg, phase, npix=160):
        img, _, area = self.render(incl_deg, phase, npix=npix)
        return float(img.sum() * area)

    def light_curve(self, incl_deg, phases, npix=160, normalise=True):
        fl = np.array([self.flux(incl_deg, p, npix=npix) for p in np.atleast_1d(phases)])
        if normalise:
            fl = fl / fl.max()
        return fl
