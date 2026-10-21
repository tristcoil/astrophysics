"""Libraries of wide-binary orbits in Newtonian gravity, MOND, and in between.

A family of force laws runs from Newton (alpha = 0) to QUMOND with the
Galactic external field (alpha = 1):

    g_alpha = g_N + alpha (g_QUMOND - g_N),

so alpha is a "dial" like Banik et al.'s (2024) alpha_grav. (Theirs
interpolates the predicted distributions rather than the force; both are 0
for Newton and 1 for MOND. Near alpha = 0 both scale with the asymptotic
velocity boost, sqrt(1 + 0.42 alpha) - 1.)

As in Banik & Zhao (2018), each binary is a test particle orbiting a point
mass M (the pair's total mass) in a field that is axisymmetric about the
Galaxy's pull (+z here). In units of the MOND radius r_M = sqrt(GM/a0) and
GM = 1 the orbit depends only on a/r_M, e and the tilt of the orbital pole
to the external field. Orbits start at r = a perpendicular to the field with
speed sqrt(-r.g) and radial velocity fraction e (their definition of a and e
for a non-Keplerian force), run for 20 revolutions with RK4 and time step
0.01 (r^3/GM)^(1/2), and are sampled uniformly in time.
"""

from __future__ import annotations

import numpy as np
from numba import njit, prange

from config import DERIVED_DIR
import qumond

#: Library grid (see module docstring)
LOG_A_GRID = np.round(np.arange(-2.0, 1.2001, 0.1), 3)       # log10(a / r_M)
E_GRID = np.round(np.arange(0.01, 1.0, 0.02), 3)              # 50 values, 0.01 ... 0.99
COSI_GRID = (np.arange(6) + 0.5) / 6.0                         # pole tilt to the field
N_STORE = 400                                                  # states kept per orbit
N_REV = 20
DT_FAC = 0.01
ALPHA_GRID = np.array([-0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5])


@njit(cache=True, fastmath=True)
def _accel(x, y, z, alpha, lnr0, dlnr, nlnr, dth, nth, R, U):
    r2 = x * x + y * y + z * z
    r = np.sqrt(r2)
    lnr = np.log(r)
    mu = z / r
    if mu > 1.0:
        mu = 1.0
    elif mu < -1.0:
        mu = -1.0
    th = np.arccos(mu)
    fi = (lnr - lnr0) / dlnr
    if fi < 0.0:
        fi = 0.0
    elif fi > nlnr - 1.000001:
        fi = nlnr - 1.000001
    i = int(fi)
    ti = fi - i
    fj = th / dth
    if fj > nth - 1.000001:
        fj = nth - 1.000001
    j = int(fj)
    tj = fj - j
    Rv = ((1 - ti) * ((1 - tj) * R[i, j] + tj * R[i, j + 1])
          + ti * ((1 - tj) * R[i + 1, j] + tj * R[i + 1, j + 1]))
    Uv = ((1 - ti) * ((1 - tj) * U[i, j] + tj * U[i, j + 1])
          + ti * ((1 - tj) * U[i + 1, j] + tj * U[i + 1, j + 1]))
    Ra = 1.0 + alpha * (Rv - 1.0)
    Ua = alpha * Uv
    r3 = r2 * r
    r4 = r2 * r2
    ax = -Ra * x / r3 + Ua * x * z / r4
    ay = -Ra * y / r3 + Ua * y * z / r4
    az = -Ra * z / r3 - Ua * (x * x + y * y) / r4
    return ax, ay, az


@njit(cache=True)
def _orbit(a, e, cosi, alpha, lnr0, dlnr, nlnr, dth, nth, R, U, n_rev, dt_fac, n_store,
           seed, out):
    """Integrate one orbit; fill out[n_store, 6] with states sampled uniformly
    in time; return (r_min, r_max, revolutions completed, steps)."""
    np.random.seed(seed)
    # start on the x axis (perpendicular to the external field along z)
    x, y, z = a, 0.0, 0.0
    ax, ay, az = _accel(x, y, z, alpha, lnr0, dlnr, nlnr, dth, nth, R, U)
    v = np.sqrt(max(-(x * ax + y * ay + z * az), 1e-300))
    sini = np.sqrt(max(1.0 - cosi * cosi, 0.0))
    # tangential direction in the y-z plane so that the pole (x^ x v_t) makes
    # angle i with z: v_t = cos(i) y^ - ... gives pole = cos(i) z^ + sin(i) y^
    vt = v * np.sqrt(max(1.0 - e * e, 0.0))
    vx = v * e
    vy = vt * cosi
    vz = -vt * sini
    # first pass: total time for n_rev revolutions (accumulated angle)
    max_steps = 50_000_000
    ang = 0.0
    t = 0.0
    steps = 0
    rmin = a
    rmax = a
    # keep the trajectory in a growing buffer of (t, state)
    cap = 200_000
    buf = np.empty((cap, 7))
    nbuf = 0
    while ang < 2.0 * np.pi * n_rev and steps < max_steps:
        r = np.sqrt(x * x + y * y + z * z)
        dt = dt_fac * r * np.sqrt(r)
        # RK4
        k1x, k1y, k1z = vx, vy, vz
        a1x, a1y, a1z = _accel(x, y, z, alpha, lnr0, dlnr, nlnr, dth, nth, R, U)
        k2x, k2y, k2z = vx + 0.5 * dt * a1x, vy + 0.5 * dt * a1y, vz + 0.5 * dt * a1z
        a2x, a2y, a2z = _accel(x + 0.5 * dt * k1x, y + 0.5 * dt * k1y, z + 0.5 * dt * k1z,
                               alpha, lnr0, dlnr, nlnr, dth, nth, R, U)
        k3x, k3y, k3z = vx + 0.5 * dt * a2x, vy + 0.5 * dt * a2y, vz + 0.5 * dt * a2z
        a3x, a3y, a3z = _accel(x + 0.5 * dt * k2x, y + 0.5 * dt * k2y, z + 0.5 * dt * k2z,
                               alpha, lnr0, dlnr, nlnr, dth, nth, R, U)
        k4x, k4y, k4z = vx + dt * a3x, vy + dt * a3y, vz + dt * a3z
        a4x, a4y, a4z = _accel(x + dt * k3x, y + dt * k3y, z + dt * k3z,
                               alpha, lnr0, dlnr, nlnr, dth, nth, R, U)
        nx = x + dt / 6.0 * (k1x + 2 * k2x + 2 * k3x + k4x)
        ny = y + dt / 6.0 * (k1y + 2 * k2y + 2 * k3y + k4y)
        nz = z + dt / 6.0 * (k1z + 2 * k2z + 2 * k3z + k4z)
        vx = vx + dt / 6.0 * (a1x + 2 * a2x + 2 * a3x + a4x)
        vy = vy + dt / 6.0 * (a1y + 2 * a2y + 2 * a3y + a4y)
        vz = vz + dt / 6.0 * (a1z + 2 * a2z + 2 * a3z + a4z)
        r_old = r
        rn = np.sqrt(nx * nx + ny * ny + nz * nz)
        c = (x * nx + y * ny + z * nz) / (r_old * rn)
        if c > 1.0:
            c = 1.0
        ang += np.arccos(c)
        x, y, z = nx, ny, nz
        t += dt
        steps += 1
        if rn < rmin:
            rmin = rn
        if rn > rmax:
            rmax = rn
        if nbuf == cap:
            nb = np.empty((cap * 2, 7))
            nb[:cap] = buf
            buf = nb
            cap *= 2
        buf[nbuf, 0] = t
        buf[nbuf, 1] = x
        buf[nbuf, 2] = y
        buf[nbuf, 3] = z
        buf[nbuf, 4] = vx
        buf[nbuf, 5] = vy
        buf[nbuf, 6] = vz
        nbuf += 1
        if rn > 1e4 or rn < 1e-6:
            break
    # stratified uniform-in-time samples
    T = t
    k = 0
    for s in range(n_store):
        ts = (s + np.random.random()) / n_store * T
        while k < nbuf - 1 and buf[k, 0] < ts:
            k += 1
        for c6 in range(6):
            out[s, c6] = buf[k, 1 + c6]
    return rmin, rmax, ang / (2.0 * np.pi), steps


@njit(parallel=True, cache=True)
def _library(log_a, e_grid, cosi_grid, alpha, lnr0, dlnr, nlnr, dth, nth, R, U, n_rev, dt_fac,
             n_store, seed0, states, meta):
    na, ne, ni = len(log_a), len(e_grid), len(cosi_grid)
    n = na * ne * ni
    for k in prange(n):
        ia = k // (ne * ni)
        ie = (k // ni) % ne
        ii = k % ni
        out = np.empty((n_store, 6))
        rmin, rmax, revs, steps = _orbit(10.0 ** log_a[ia], e_grid[ie], cosi_grid[ii], alpha,
                                         lnr0, dlnr, nlnr, dth, nth, R, U, n_rev, dt_fac,
                                         n_store, seed0 + k, out)
        for s in range(n_store):
            for c in range(6):
                states[ia, ie, ii, s, c] = out[s, c]
        meta[ia, ie, ii, 0] = rmin
        meta[ia, ie, ii, 1] = rmax
        meta[ia, ie, ii, 2] = revs
        meta[ia, ie, ii, 3] = steps


def field_tables(yNe: float | None = None):
    """(lnr0, dlnr, nlnr, dth, nth, R, U) for the integrator, from qumond."""
    if yNe is None:
        yNe = qumond.galactic_field()["gNe_a0"]
    f = qumond.cached_field(yNe)
    lnr = f.lnr
    theta = np.arccos(np.clip(f.mu, -1, 1))
    assert np.allclose(np.diff(lnr), lnr[1] - lnr[0]) and np.allclose(np.diff(theta), theta[1] - theta[0])
    return (float(lnr[0]), float(lnr[1] - lnr[0]), len(lnr), float(theta[1] - theta[0]),
            len(theta), np.ascontiguousarray(f.ratio_r), np.ascontiguousarray(f.u_theta))


def build_library(alpha: float, *, yNe: float | None = None, quiet: bool = False):
    """Orbit library for one alpha, cached as data/derived/orbits_alpha<..>.npz."""
    import time
    tag = "" if yNe is None else f"_yNe{yNe:.4f}"
    path = DERIVED_DIR / f"orbits_alpha{alpha:+.2f}{tag}.npz"
    if path.exists():
        z = np.load(path)
        return dict(states=z["states"], meta=z["meta"], alpha=float(z["alpha"]))
    tabs = field_tables(yNe)
    na, ne, ni = len(LOG_A_GRID), len(E_GRID), len(COSI_GRID)
    states = np.zeros((na, ne, ni, N_STORE, 6), np.float64)
    meta = np.zeros((na, ne, ni, 4), np.float64)
    t0 = time.time()
    _library(LOG_A_GRID, E_GRID, COSI_GRID, float(alpha), *tabs, N_REV, DT_FAC, N_STORE,
             int(1000 * (alpha + 10)), states, meta)
    if not quiet:
        print(f"    orbit library alpha = {alpha:+.2f}: {na * ne * ni} orbits in "
              f"{time.time() - t0:.0f} s (median {np.median(meta[..., 3]):.0f} steps)", flush=True)
    np.savez(path, states=states.astype(np.float32), meta=meta, alpha=alpha,
             log_a=LOG_A_GRID, e=E_GRID, cosi=COSI_GRID)
    return dict(states=states.astype(np.float32), meta=meta, alpha=alpha)


if __name__ == "__main__":
    import time
    tabs = field_tables()
    # quick accuracy test: Newtonian (alpha = 0) orbit energy and period
    out = np.empty((N_STORE, 6))
    for e in (0.0, 0.5, 0.9, 0.99):
        t0 = time.time()
        rmin, rmax, revs, steps = _orbit(0.5, e, 0.5, 0.0, *tabs, N_REV, DT_FAC, N_STORE, 1, out)
        r = np.linalg.norm(out[:, :3], axis=1)
        v2 = np.sum(out[:, 3:] ** 2, axis=1)
        E = 0.5 * v2 - 1.0 / r
        print(f"e={e}: r_min {rmin:.4f} (exact {0.5 * (1 - e):.4f}) r_max {rmax:.4f} "
              f"(exact {0.5 * (1 + e):.4f}); energy {E.mean():.6f} +- {E.std():.1e} "
              f"(exact {-1 / (2 * 0.5):.6f}); mean r/a {r.mean() / 0.5:.4f} (exact {1 + e * e / 2:.4f}); "
              f"{steps} steps {time.time() - t0:.2f} s")
