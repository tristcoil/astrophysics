"""Trace the runaways through a model of the Milky Way's gravity.

We use McMillan (2017)'s Milky Way potential (bulge, thin and thick stellar
discs, two gas discs and a dark halo fitted to rotation-curve and other data),
as implemented in galpy, and integrate each posterior sample's orbit backwards
and forwards in time with galpy's C integrator.

For a star moving at 2,000 km/s the Galaxy hardly matters over a million
years: the pull near the Sun, v_c^2 / R ~ (233 km/s)^2 / 8 kpc, changes its
speed by only ~7 km/s per Myr. The integration matters for the slower,
older runaways (LP 40-365 stars, US 708) and for the escape speed, which
says whether a star will ever come back.

Coordinates: positions and velocities enter and leave in astropy's
Galactocentric frame (Sun at x = -R0, Galactic rotation toward +y at the Sun).
galpy measures azimuth in the direction of rotation, which here means
mirroring x; the McMillan potential is axisymmetric, so the mirror changes
nothing physical and is undone on the way out.
"""

from __future__ import annotations

from functools import lru_cache

import astropy.units as u
import numpy as np

#: McMillan (2017)'s own scales, which galpy's implementation is defined in.
RO_KPC, VO_KMS = 8.21, 233.1


@lru_cache(maxsize=1)
def potential():
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from galpy.potential.mwpotentials import McMillan17
    return McMillan17


def phi(R, z):
    """Potential (km/s)^2 at cylindrical (R, z) kpc; zero at infinity.
    (galpy's McMillan17 has physical output switched on, so natural units
    are requested explicitly and scaled here, once.)"""
    from galpy.potential import evaluatePotentials
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.asarray(evaluatePotentials(potential(), np.asarray(R) / RO_KPC,
                                             np.asarray(z) / RO_KPC,
                                             use_physical=False)) * VO_KMS ** 2


@lru_cache(maxsize=1)
def _phi_grid():
    """Phi tabulated on a log-R x log-|z| grid (galpy evaluates McMillan17's
    expansion components point by point in Python, ~2 ms each, so a few
    thousand direct calls per star would take minutes). Cached on disk."""
    from scipy.interpolate import RegularGridInterpolator
    from config import CACHE_DIR
    path = CACHE_DIR / "mcmillan17_phi_grid.npz"
    if path.exists():
        g = np.load(path)
        lr, lz, tab = g["lr"], g["lz"], g["phi"]
    else:
        lr = np.linspace(np.log(0.05), np.log(100.0), 90)
        lz = np.linspace(np.log(0.002), np.log(100.0), 70)
        R, Z = np.meshgrid(np.exp(lr), np.exp(lz), indexing="ij")
        tab = np.vectorize(lambda r, z: float(phi(r, z)))(R, Z)
        np.savez(path, lr=lr, lz=lz, phi=tab)
    return RegularGridInterpolator((lr, lz), tab, bounds_error=False, fill_value=None)


def phi_fast(R, z):
    """Interpolated Phi (km/s)^2; within 0.1% of phi() (test_physics.py)."""
    R = np.clip(np.asarray(R, float), 0.05, 100.0)
    az = np.clip(np.abs(np.asarray(z, float)), 0.002, 100.0)
    return _phi_grid()(np.column_stack([np.log(R).ravel(), np.log(az).ravel()])).reshape(R.shape)


def escape_speed(x) -> np.ndarray:
    """sqrt(-2 Phi) at Galactocentric positions x (3, n) kpc, km/s."""
    R = np.hypot(x[0], x[1])
    return np.sqrt(-2 * phi_fast(R, x[2]))


@lru_cache(maxsize=1)
def _vcirc_table():
    from galpy.potential import vcirc
    R = np.geomspace(0.05, 60.0, 400)
    v = np.array([vcirc(potential(), r / RO_KPC, use_physical=False) for r in R]) * VO_KMS
    return R, v


def vcirc(R) -> np.ndarray:
    """Midplane circular speed (km/s) at radius R (kpc), interpolated."""
    Rt, vt = _vcirc_table()
    return np.interp(np.asarray(R, float), Rt, vt)


def _to_galpy(x, v):
    xg, yg, z = -x[0], x[1], x[2]
    vxg, vyg, vz = -v[0], v[1], v[2]
    R = np.hypot(xg, yg)
    phi_ = np.arctan2(yg, xg)
    vR = (xg * vxg + yg * vyg) / R
    vT = (xg * vyg - yg * vxg) / R
    return np.array([R / RO_KPC, vR / VO_KMS, vT / VO_KMS, z / RO_KPC, vz / VO_KMS, phi_]).T


def integrate(x, v, t_myr) -> tuple[np.ndarray, np.ndarray]:
    """Integrate orbits from present positions x and velocities v (3, n) over
    times t_myr (1D, starting at 0; negative = past). Returns positions and
    velocities, shape (3, n, nt), in the astropy Galactocentric frame."""
    import warnings
    from galpy.orbit import Orbit
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        o = Orbit(_to_galpy(x, v), ro=RO_KPC, vo=VO_KMS)
        o.integrate(np.asarray(t_myr) * u.Myr, potential(), method="dop853_c", progressbar=False)
        xs, ys, zs = (o.x(np.asarray(t_myr) * u.Myr), o.y(np.asarray(t_myr) * u.Myr),
                      o.z(np.asarray(t_myr) * u.Myr))
        vxs, vys, vzs = (o.vx(np.asarray(t_myr) * u.Myr), o.vy(np.asarray(t_myr) * u.Myr),
                         o.vz(np.asarray(t_myr) * u.Myr))
    pos = np.array([-np.atleast_2d(xs), np.atleast_2d(ys), np.atleast_2d(zs)])
    vel = np.array([-np.atleast_2d(vxs), np.atleast_2d(vys), np.atleast_2d(vzs)])
    return pos, vel


def last_disc_crossing(pos, vel, t_myr):
    """Most recent midplane crossing in the past for each orbit.

    pos, vel: (3, n, nt) over t_myr going from 0 into the past. Returns the
    time (Myr, positive, nan if none within the window), and position and
    velocity at the crossing (linear interpolation between steps)."""
    z = pos[2]
    s = np.sign(z)
    flip = (s[:, 1:] * s[:, :-1]) <= 0
    has = flip.any(axis=1)
    k = np.where(has, flip.argmax(axis=1), 0)
    n = np.arange(z.shape[0])
    z0, z1 = z[n, k], z[n, k + 1]
    f = np.where(z1 != z0, z0 / (z0 - z1), 0.0)
    t0, t1 = t_myr[k], t_myr[k + 1]
    tc = -(t0 + f * (t1 - t0))
    xc = pos[:, n, k] + f * (pos[:, n, k + 1] - pos[:, n, k])
    vc = vel[:, n, k] + f * (vel[:, n, k + 1] - vel[:, n, k])
    tc = np.where(has, tc, np.nan)
    return tc, xc, vc


def launch_from_orbit(pos, vel, t_myr):
    """Launch speed with the potential: velocity at the last midplane crossing
    minus the circular velocity there. Stars still heading for the plane use
    their present position and velocity (as in the straight-line method)."""
    from kinematics import rotation_velocity
    tc, xc, vc = last_disc_crossing(pos, vel, t_myr)
    approaching = ~np.isfinite(tc)
    xc = np.where(approaching, pos[:, :, 0], xc)
    vc = np.where(approaching, vel[:, :, 0], vc)
    v_ej = np.linalg.norm(vc - rotation_velocity(xc, vcirc), axis=0)
    return tc, xc, v_ej
