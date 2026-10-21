"""Numerical orbit model for LAGEOS-type satellites (compiled with numba).

Accelerations, all in the GCRS, following the IERS Conventions (2010):

* Earth's gravity field to degree and order 30 (GOCO06s, time-variable terms
  evaluated per arc), plus solid-Earth, ocean and pole tides as time-varying
  coefficients interpolated from a grid (``tides.py``). The field is evaluated in
  the Earth-fixed frame with the recursion of Montenbruck & Gill (2000, sec. 3.2).
* The Sun, Moon, Venus, Mars, Jupiter and Saturn as point masses (DE440).
* Solar radiation pressure on a sphere (cannonball, scale factor C_R) with a
  conical Earth shadow.
* General relativity (IERS 2010 eq. 10.12, PPN beta = gamma = 1):
  Schwarzschild (Earth's mass), de Sitter (geodetic precession in the Sun's field,
  ~19.2 mas/yr about the ecliptic pole) and Lense-Thirring (frame dragging by the
  Earth's spin). The last two carry scale factors mu_dS and mu_LT, so frame
  dragging can be switched off (mu_LT = 0) and then measured from the data.
* Empirical along-track accelerations (constant + once per revolution) and a
  constant radial one, fitted per arc: they absorb the small forces we do not
  model (thermal re-radiation, Earth albedo, drag by charged particles). There is
  deliberately no cross-track empirical term: a once-per-revolution cross-track
  force would turn the orbital plane and soak up the very drift we measure.

The integrator is a fixed-step 8th-order Runge-Kutta (Dormand-Prince DOP853
coefficients from SciPy). Partial derivatives for the least-squares fits come
from variational equations with the gravity gradient of the central field plus J2
(a two-body Jacobian alone is 12% off after one day and 50% after two).
"""

from __future__ import annotations

import numpy as np
import scipy.integrate._ivp.dop853_coefficients as _dop
from numba import njit

RK_A = np.ascontiguousarray(_dop.A[:12, :12])
RK_B = np.ascontiguousarray(_dop.B)
RK_C = np.ascontiguousarray(_dop.C[:12])

#: Parameters with partials: C_R, along-track (constant, cos u, sin u), radial constant.
PARAM_NAMES = ("cr", "t0", "tc", "ts", "r0")
NPAR = len(PARAM_NAMES)
NSTATE = 6 + 36 + 6 * NPAR


@njit(cache=True)
def _interp_lin(grid, t0, dt, t):
    x = (t - t0) / dt
    k = int(np.floor(x))
    if k < 0:
        k = 0
    if k > grid.shape[0] - 2:
        k = grid.shape[0] - 2
    w = x - k
    return k, w


@njit(cache=True)
def _lagrange4(grid, t0, dt, t):
    """Cubic Lagrange interpolation of rows of a uniformly sampled (N, d) array."""
    x = (t - t0) / dt
    k = int(np.floor(x)) - 1
    if k < 0:
        k = 0
    if k > grid.shape[0] - 4:
        k = grid.shape[0] - 4
    u = x - k
    w0 = -(u - 1.0) * (u - 2.0) * (u - 3.0) / 6.0
    w1 = u * (u - 2.0) * (u - 3.0) / 2.0
    w2 = -u * (u - 1.0) * (u - 3.0) / 2.0
    w3 = u * (u - 1.0) * (u - 2.0) / 6.0
    out = np.empty(grid.shape[1:])
    out[...] = w0 * grid[k] + w1 * grid[k + 1] + w2 * grid[k + 2] + w3 * grid[k + 3]
    return out


@njit(cache=True)
def geopotential(x, y, z, Cnm, Snm, nmax, gm, re):
    """Acceleration (Earth-fixed) from unnormalised coefficients, central term included."""
    N1 = nmax + 2
    r2 = x * x + y * y + z * z
    rho = re * re / r2
    x0 = re * x / r2
    y0 = re * y / r2
    z0 = re * z / r2
    V = np.zeros((N1 + 1, N1 + 1))
    W = np.zeros((N1 + 1, N1 + 1))
    V[0, 0] = re / np.sqrt(r2)
    V[1, 0] = z0 * V[0, 0]
    for n in range(2, N1 + 1):
        V[n, 0] = ((2 * n - 1) * z0 * V[n - 1, 0] - (n - 1) * rho * V[n - 2, 0]) / n
    for m in range(1, N1 + 1):
        V[m, m] = (2 * m - 1) * (x0 * V[m - 1, m - 1] - y0 * W[m - 1, m - 1])
        W[m, m] = (2 * m - 1) * (x0 * W[m - 1, m - 1] + y0 * V[m - 1, m - 1])
        if m + 1 <= N1:
            V[m + 1, m] = (2 * m + 1) * z0 * V[m, m]
            W[m + 1, m] = (2 * m + 1) * z0 * W[m, m]
        for n in range(m + 2, N1 + 1):
            V[n, m] = ((2 * n - 1) * z0 * V[n - 1, m] - (n + m - 1) * rho * V[n - 2, m]) / (n - m)
            W[n, m] = ((2 * n - 1) * z0 * W[n - 1, m] - (n + m - 1) * rho * W[n - 2, m]) / (n - m)
    ax = 0.0
    ay = 0.0
    az = 0.0
    for m in range(nmax + 1):
        for n in range(m, nmax + 1):
            C = Cnm[n, m]
            if m == 0:
                ax -= C * V[n + 1, 1]
                ay -= C * W[n + 1, 1]
                az -= (n + 1) * C * V[n + 1, 0]
            else:
                S = Snm[n, m]
                fac = 0.5 * (n - m + 1) * (n - m + 2)
                ax += 0.5 * (-C * V[n + 1, m + 1] - S * W[n + 1, m + 1]) \
                    + fac * (C * V[n + 1, m - 1] + S * W[n + 1, m - 1])
                ay += 0.5 * (-C * W[n + 1, m + 1] + S * V[n + 1, m + 1]) \
                    + fac * (-C * W[n + 1, m - 1] + S * V[n + 1, m - 1])
                az += (n - m + 1) * (-C * V[n + 1, m] - S * W[n + 1, m])
    k = gm / (re * re)
    return k * ax, k * ay, k * az


@njit(cache=True)
def shadow_factor(r, sun, re, rsun):
    """Fraction of the solar disk visible from r (conical shadow, Montenbruck & Gill 3.4.2)."""
    d = sun - r
    dn = np.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2)
    rn = np.sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2)
    a = np.arcsin(min(1.0, rsun / dn))
    b = np.arcsin(min(1.0, re / rn))
    cosc = -(r[0] * d[0] + r[1] * d[1] + r[2] * d[2]) / (rn * dn)
    c = np.arccos(max(-1.0, min(1.0, cosc)))
    if c >= a + b:
        return 1.0
    if c < b - a:
        return 0.0
    if c < a - b:                      # Sun larger than Earth's disk: annular
        return 1.0 - b * b / (a * a)
    xx = (c * c + a * a - b * b) / (2.0 * c)
    yy = np.sqrt(max(0.0, a * a - xx * xx))
    area = a * a * np.arccos(max(-1.0, min(1.0, xx / a))) \
        + b * b * np.arccos(max(-1.0, min(1.0, (c - xx) / b))) - c * yy
    return 1.0 - area / (np.pi * a * a)


@njit(cache=True)
def acceleration(t, r, v, M_f, M_i, Qg, Wg, tug, Cs, Ss, dCg, dSg, bodies, gmb, dsg, prm, out_dp):
    """Total GCRS acceleration at TT seconds t (from the arc epoch) and dA/dparams.

    M_f: float settings [gm, re, c, j_per_mass, am, p1au, au, rsun, mu_lt, mu_ds,
         t_rot0, dt_rot, t_tide0, dt_tide, t_body0, dt_body]
    M_i: int settings [nmax, ntide, nbodies]
    """
    gm = M_f[0]
    re = M_f[1]
    c = M_f[2]
    jpm = M_f[3]
    am = M_f[4]
    p1au = M_f[5]
    au = M_f[6]
    rsun = M_f[7]
    mu_lt = M_f[8]
    mu_ds = M_f[9]
    nmax = M_i[0]
    ntide = M_i[1]
    nb = M_i[2]

    # --- Earth orientation: GCRS -> ITRS -------------------------------------
    k, w = _interp_lin(Qg, M_f[10], M_f[11], t)
    Q = (1.0 - w) * Qg[k] + w * Qg[k + 1]
    Wm = (1.0 - w) * Wg[k] + w * Wg[k + 1]
    tu = (1.0 - w) * tug[k] + w * tug[k + 1]
    era = 2.0 * np.pi * ((tu % 1.0) + 0.7790572732640 + 0.00273781191135448 * tu)
    ce = np.cos(era)
    se = np.sin(era)
    # r_cirs = Q^T r
    rc = np.empty(3)
    for i in range(3):
        rc[i] = Q[0, i] * r[0] + Q[1, i] * r[1] + Q[2, i] * r[2]
    rt = np.empty(3)
    rt[0] = ce * rc[0] + se * rc[1]
    rt[1] = -se * rc[0] + ce * rc[1]
    rt[2] = rc[2]
    ri = np.empty(3)
    for i in range(3):
        ri[i] = Wm[0, i] * rt[0] + Wm[1, i] * rt[1] + Wm[2, i] * rt[2]

    # --- gravity field with tides --------------------------------------------
    kt, wt = _interp_lin(dCg, M_f[12], M_f[13], t)
    Cn = Cs.copy()
    Sn = Ss.copy()
    for n in range(ntide + 1):
        for m in range(n + 1):
            Cn[n, m] += (1.0 - wt) * dCg[kt, n, m] + wt * dCg[kt + 1, n, m]
            Sn[n, m] += (1.0 - wt) * dSg[kt, n, m] + wt * dSg[kt + 1, n, m]
    gx, gy, gz = geopotential(ri[0], ri[1], ri[2], Cn, Sn, nmax, gm, re)
    # back: a_tirs = Wm a_itrs ; a_cirs = R3(-era) a_tirs ; a_gcrs = Q a_cirs
    at = np.empty(3)
    for i in range(3):
        at[i] = Wm[i, 0] * gx + Wm[i, 1] * gy + Wm[i, 2] * gz
    ac = np.empty(3)
    ac[0] = ce * at[0] - se * at[1]
    ac[1] = se * at[0] + ce * at[1]
    ac[2] = at[2]
    a = np.empty(3)
    for i in range(3):
        a[i] = Q[i, 0] * ac[0] + Q[i, 1] * ac[1] + Q[i, 2] * ac[2]

    # --- third bodies -----------------------------------------------------------
    B = _lagrange4(bodies, M_f[14], M_f[15], t)      # (nb, 3): Sun first, Moon second
    for j in range(nb):
        s = B[j]
        d0 = s[0] - r[0]
        d1 = s[1] - r[1]
        d2 = s[2] - r[2]
        dn3 = (d0 * d0 + d1 * d1 + d2 * d2) ** 1.5
        sn3 = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]) ** 1.5
        a[0] += gmb[j] * (d0 / dn3 - s[0] / sn3)
        a[1] += gmb[j] * (d1 / dn3 - s[1] / sn3)
        a[2] += gmb[j] * (d2 / dn3 - s[2] / sn3)

    # --- solar radiation pressure ------------------------------------------------
    sun = B[0]
    d = r - sun
    dn = np.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2)
    nu = shadow_factor(r, sun, re, rsun)
    fs = am * p1au * (au / dn) ** 2 * nu / dn
    for i in range(3):
        out_dp[i, 0] = fs * d[i]                 # d a / d C_R
        a[i] += prm[0] * fs * d[i]

    # --- relativity -----------------------------------------------------------------
    rn = np.sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2)
    v2 = v[0] ** 2 + v[1] ** 2 + v[2] ** 2
    rv = r[0] * v[0] + r[1] * v[1] + r[2] * v[2]
    f = gm / (c * c * rn ** 3)
    for i in range(3):
        a[i] += f * ((4.0 * gm / rn - v2) * r[i] + 4.0 * rv * v[i])
    # Earth's spin axis = CIP = third column of Q (GCRS <- CIRS)
    J0 = jpm * Q[0, 2]
    J1 = jpm * Q[1, 2]
    J2_ = jpm * Q[2, 2]
    rJ = r[0] * J0 + r[1] * J1 + r[2] * J2_
    hx = r[1] * v[2] - r[2] * v[1]
    hy = r[2] * v[0] - r[0] * v[2]
    hz = r[0] * v[1] - r[1] * v[0]
    vJx = v[1] * J2_ - v[2] * J1
    vJy = v[2] * J0 - v[0] * J2_
    vJz = v[0] * J1 - v[1] * J0
    flt = mu_lt * 2.0 * f
    a[0] += flt * (3.0 / (rn * rn) * hx * rJ + vJx)
    a[1] += flt * (3.0 / (rn * rn) * hy * rJ + vJy)
    a[2] += flt * (3.0 / (rn * rn) * hz * rJ + vJz)
    # de Sitter: 2 Omega_dS x v
    om = _lagrange4(dsg, M_f[14], M_f[15], t)
    a[0] += mu_ds * 2.0 * (om[1] * v[2] - om[2] * v[1])
    a[1] += mu_ds * 2.0 * (om[2] * v[0] - om[0] * v[2])
    a[2] += mu_ds * 2.0 * (om[0] * v[1] - om[1] * v[0])

    # --- empirical accelerations (RTN) -------------------------------------------
    hn = np.sqrt(hx * hx + hy * hy + hz * hz)
    Rx, Ry, Rz = r[0] / rn, r[1] / rn, r[2] / rn
    Nx, Ny, Nz = hx / hn, hy / hn, hz / hn
    Tx = Ny * Rz - Nz * Ry
    Ty = Nz * Rx - Nx * Rz
    Tz = Nx * Ry - Ny * Rx
    # ascending-node direction z x N, normalised
    ex, ey = -Ny, Nx
    en = np.sqrt(ex * ex + ey * ey)
    ex /= en
    ey /= en
    cosu = Rx * ex + Ry * ey
    # (N x e_node) . R
    bx = -Nz * ey
    by = Nz * ex
    bz = Nx * ey - Ny * ex
    sinu = Rx * bx + Ry * by + Rz * bz
    out_dp[0, 1] = Tx
    out_dp[1, 1] = Ty
    out_dp[2, 1] = Tz
    out_dp[0, 2] = cosu * Tx
    out_dp[1, 2] = cosu * Ty
    out_dp[2, 2] = cosu * Tz
    out_dp[0, 3] = sinu * Tx
    out_dp[1, 3] = sinu * Ty
    out_dp[2, 3] = sinu * Tz
    out_dp[0, 4] = Rx
    out_dp[1, 4] = Ry
    out_dp[2, 4] = Rz
    for i in range(3):
        a[i] += prm[1] * out_dp[i, 1] + prm[2] * out_dp[i, 2] + prm[3] * out_dp[i, 3] \
            + prm[4] * out_dp[i, 4]
    return a


@njit(cache=True)
def _j2_accel(x, y, z, gm, re, j2):
    r2 = x * x + y * y + z * z
    r = np.sqrt(r2)
    k = -1.5 * j2 * gm * re * re / (r2 * r2 * r)
    f = 5.0 * z * z / r2
    return k * x * (1.0 - f), k * y * (1.0 - f), k * z * (3.0 - f)


@njit(cache=True)
def _deriv(t, Y, M_f, M_i, Qg, Wg, tug, Cs, Ss, dCg, dSg, bodies, gmb, dsg, prm, with_var):
    dY = np.zeros(Y.shape[0])
    r = Y[0:3]
    v = Y[3:6]
    dp = np.zeros((3, NPAR))
    a = acceleration(t, r, v, M_f, M_i, Qg, Wg, tug, Cs, Ss, dCg, dSg, bodies, gmb, dsg, prm, dp)
    dY[0:3] = v
    dY[3:6] = a
    if with_var:
        # gravity gradient of the central field plus J2 (J2 axis taken along the
        # GCRS z-axis, 0.1 deg from the true pole: plenty for Gauss-Newton partials)
        rn = np.sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2)
        g = -M_f[0] / rn ** 3
        G = np.empty((3, 3))
        for i in range(3):
            for j in range(3):
                G[i, j] = g * ((1.0 if i == j else 0.0) - 3.0 * r[i] * r[j] / (rn * rn))
        j2 = -Cs[2, 0]
        dlt = 1.0
        for j in range(3):
            rp = r.copy()
            rm = r.copy()
            rp[j] += dlt
            rm[j] -= dlt
            ap = _j2_accel(rp[0], rp[1], rp[2], M_f[0], M_f[1], j2)
            am_ = _j2_accel(rm[0], rm[1], rm[2], M_f[0], M_f[1], j2)
            G[0, j] += (ap[0] - am_[0]) / (2 * dlt)
            G[1, j] += (ap[1] - am_[1]) / (2 * dlt)
            G[2, j] += (ap[2] - am_[2]) / (2 * dlt)
        ncol = 6 + NPAR
        for col in range(ncol):
            base = 6 + col * 6               # column-major blocks of 6 rows
            for i in range(3):
                dY[base + i] = Y[base + 3 + i]
            for i in range(3):
                s = 0.0
                for j in range(3):
                    s += G[i, j] * Y[base + j]
                if col >= 6:
                    s += dp[i, col - 6]
                dY[base + 3 + i] = s
    return dY


@njit(cache=True)
def propagate(y0, t_start, h, nsteps, stride, M_f, M_i, Qg, Wg, tug, Cs, Ss, dCg, dSg,
              bodies, gmb, dsg, prm, with_var):
    """Fixed-step DOP853 from t_start; returns times and states every ``stride`` steps."""
    n = y0.shape[0]
    nout = nsteps // stride + 1
    T = np.empty(nout)
    Yout = np.empty((nout, n))
    y = y0.copy()
    t = t_start
    T[0] = t
    Yout[0] = y
    K = np.empty((12, n))
    j = 1
    for step in range(1, nsteps + 1):
        K[0] = _deriv(t, y, M_f, M_i, Qg, Wg, tug, Cs, Ss, dCg, dSg, bodies, gmb, dsg, prm, with_var)
        for s in range(1, 12):
            ys = y.copy()
            for q in range(s):
                if RK_A[s, q] != 0.0:
                    ys += h * RK_A[s, q] * K[q]
            K[s] = _deriv(t + RK_C[s] * h, ys, M_f, M_i, Qg, Wg, tug, Cs, Ss, dCg, dSg, bodies,
                          gmb, dsg, prm, with_var)
        for q in range(12):
            y += h * RK_B[q] * K[q]
        t = t_start + step * h
        if step % stride == 0:
            T[j] = t
            Yout[j] = y
            j += 1
    return T[:j], Yout[:j]
