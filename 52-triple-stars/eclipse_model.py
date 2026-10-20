"""Light curves of an eclipsing binary whose orbit precesses.

Geometry
--------
The third star makes the inner orbit's angular momentum precess on a cone
around the system's total angular momentum J. We describe the inner orbit in
the *invariable frame* (z' along J): inclination ``eps`` to the invariable
plane, node ``h`` (regressing at a constant rate), argument of periastron
``g`` (advancing: apsidal motion) and mean longitude ``lam = h + g + M``.
The sky frame follows by a rotation of angle ``I`` (the inclination of J to
the line of sight) about the invariable plane's line of nodes; z is toward
the observer. Rotating about the line of sight does not change a light
curve, so J's position angle on the sky is irrelevant here.

The inclination of the binary to the sky then follows Torres (2001, eq. 2):

    cos i = cos I cos eps - sin I sin eps cos(h)   (h measured from that node)

but the full 3-D positions also give the eclipse *times*: as the orbital
plane turns, the moment of conjunction slides along the orbit. This is what
makes old eclipses fall at "wrong" phases of a fixed ephemeris.

Photometry
----------
Two spheres with quadratic limb darkening; the fraction of a star's light
hidden by the other (an opaque disc) is integrated numerically over the
eclipsed star's radius: closed form for the rings fully covered, Gauss-Legendre
for the rings the occulter's edge crosses (tested in test_physics.py against
the analytic uniform-disc overlap and a brute-force pixel integration). Reflection and
ellipsoidal light are not modelled: they are below the plate noise, and for
TESS we fit them out around each eclipse.
"""

from __future__ import annotations

import numba as nb
import numpy as np

@nb.njit(cache=True, fastmath=False)
def _alpha(r, z, p):
    """Half-angle of the arc of the circle of radius r (centred on the eclipsed
    star) that lies inside the occulting disc of radius p at distance z."""
    if z <= 0.0:
        return np.pi if r < p else 0.0
    if r + z <= p:
        return np.pi
    if r >= z + p or z >= r + p:
        return 0.0
    c = (r * r + z * z - p * p) / (2.0 * r * z)
    if c >= 1.0:
        return 0.0
    if c <= -1.0:
        return np.pi
    return np.arccos(c)


_NSEG = 96
_xs, _ws = np.polynomial.legendre.leggauss(_NSEG)


@nb.njit(cache=True)
def _cum_flux(r, u1, u2):
    """Flux of a quadratically limb-darkened disc inside radius r (r <= 1),
    in units where the whole disc gives pi (1 - u1/3 - u2/6)."""
    a, b, c = 1.0 - u1 - u2, u1 + 2.0 * u2, -u2
    m = np.sqrt(max(0.0, 1.0 - r * r))
    F1 = a / 2.0 + b / 3.0 + c / 4.0
    Fm = a * m * m / 2.0 + b * m ** 3 / 3.0 + c * m ** 4 / 4.0
    return 2.0 * np.pi * (F1 - Fm)


@nb.njit(cache=True)
def occulted_fraction(z, p, u1, u2):
    """Fraction of the flux of a limb-darkened star (radius 1) hidden by an
    opaque disc of radius p whose centre is at projected distance z.

    Rings fully inside the occulter are integrated in closed form; the rings
    the occulter's edge crosses (|z - p| < r < z + p) by Gauss-Legendre in
    s = 1 - mu, which resolves the limb."""
    if z >= 1.0 + p:
        return 0.0
    tot = np.pi * (1.0 - u1 / 3.0 - u2 / 6.0)
    full = 0.0
    if z < p:
        ra = min(1.0, p - z)
        full = _cum_flux(ra, u1, u2)
        if ra >= 1.0:
            return 1.0
    r_lo = abs(z - p)
    r_hi = min(1.0, z + p)
    if r_hi <= r_lo:
        return full / tot
    s_lo = 1.0 - np.sqrt(max(0.0, 1.0 - r_lo * r_lo))
    s_hi = 1.0 - np.sqrt(max(0.0, 1.0 - r_hi * r_hi))
    part = 0.0
    half = 0.5 * (s_hi - s_lo)
    mid = 0.5 * (s_hi + s_lo)
    for k in range(_NSEG):
        sk = mid + half * _xs[k]
        mu = 1.0 - sk
        r = np.sqrt(max(0.0, 1.0 - mu * mu))
        I = 1.0 - u1 * sk - u2 * sk * sk
        # r dr = mu ds; the ring's covered arc is 2 r alpha -> weight 2 alpha r dr
        part += _ws[k] * half * I * mu * 2.0 * _alpha(r, z, p)
    return (full + part) / tot


@nb.njit(cache=True)
def _kepler_E(M, e):
    E = M if e < 0.8 else np.pi
    for _ in range(50):
        f = E - e * np.sin(E) - M
        d = f / (1.0 - e * np.cos(E))
        E -= d
        if abs(d) < 1e-13:
            break
    return E


@nb.njit(cache=True)
def sky_positions(t, t_ref, lam0, n_lam, h0, hdot, g0, gdot, e, eps, I):
    """Position (X, Y, Z) of the secondary relative to the primary in units of
    the semimajor axis; Z > 0 means the secondary is closer to us.

    Angles in radians, rates in radians per day, times in days."""
    nt = t.shape[0]
    out = np.empty((nt, 3))
    ce, se = np.cos(eps), np.sin(eps)
    cI, sI = np.cos(I), np.sin(I)
    for j in range(nt):
        dt = t[j] - t_ref
        h = h0 + hdot * dt
        g = g0 + gdot * dt
        lam = lam0 + n_lam * dt
        M = (lam - h - g) % (2.0 * np.pi)
        if e > 0.0:
            E = _kepler_E(M, e)
            cosE, sinE = np.cos(E), np.sin(E)
            xo = cosE - e
            yo = np.sqrt(1.0 - e * e) * sinE
        else:
            xo = np.cos(M)
            yo = np.sin(M)
        # rotate by g (in the orbital plane), eps (about the node), h (about z')
        cg, sg = np.cos(g), np.sin(g)
        x1 = cg * xo - sg * yo
        y1 = sg * xo + cg * yo
        x2 = x1
        y2 = ce * y1
        z2 = se * y1
        ch, sh = np.cos(h), np.sin(h)
        x3 = ch * x2 - sh * y2
        y3 = sh * x2 + ch * y2
        z3 = z2
        # invariable frame -> sky (rotation by I about the x' axis)
        out[j, 0] = x3
        out[j, 1] = cI * y3 - sI * z3
        out[j, 2] = sI * y3 + cI * z3
    return out


@nb.njit(cache=True)
def binary_flux(xyz, r1, r2, L1, L2, L3, u1a, u2a, u1b, u2b):
    """Normalised flux for sky positions xyz (units of a). r1, r2 fractional
    radii; L1, L2 surface-integrated luminosities in the band; L3 third light
    (same units)."""
    nt = xyz.shape[0]
    f = np.empty(nt)
    tot = L1 + L2 + L3
    for j in range(nt):
        rho = np.sqrt(xyz[j, 0] ** 2 + xyz[j, 1] ** 2)
        lost = 0.0
        if rho < r1 + r2:
            if xyz[j, 2] > 0.0:      # secondary in front: primary eclipse
                lost = L1 * occulted_fraction(rho / r1, r2 / r1, u1a, u2a)
            else:                    # primary in front
                lost = L2 * occulted_fraction(rho / r2, r1 / r2, u1b, u2b)
        f[j] = (tot - lost) / tot
    return f


def inclination_to_sky(t, t_ref, h0, hdot, eps, I):
    """Inclination of the inner orbit to the plane of the sky (radians)."""
    h = h0 + hdot * (np.asarray(t, float) - t_ref)
    # orbit normal n = R_x(I) (sin eps sin h, -sin eps cos h, cos eps); i = angle to z
    nz = np.cos(I) * np.cos(eps) - np.sin(I) * np.sin(eps) * np.cos(h)
    return np.arccos(np.clip(nz, -1.0, 1.0))


def planck_ratio(t1: float, t2: float, wavelength_nm: float) -> float:
    """Surface-brightness ratio B_lambda(T2)/B_lambda(T1) at one wavelength."""
    c2 = 1.438776877e7      # nm K
    x1 = c2 / (wavelength_nm * t1)
    x2 = c2 / (wavelength_nm * t2)
    return float(np.expm1(x1) / np.expm1(x2))


#: effective wavelengths: photographic blue plates (APASS-calibrated, B) and TESS
LAMBDA_NM = {"B": 440.0, "V": 550.0, "TESS": 800.0}
