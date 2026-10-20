"""Compiled (numba) likelihood core: the same mathematics as kepler.py and
models.OrbitModel.loglike_phys, about 50x faster. test_physics.py checks that
both give the same log-likelihood to 1e-9.

Element array layout: el[k] = (a, e, inc, aop, pan, tau) for planets in
config.INNER_TO_OUTER order (e, d, c, b) - or any number of planets, inner
first; masses in Msun in the same order.
"""

from __future__ import annotations

import math

import numpy as np
from numba import njit

import kepler as K

TWO_PI = 2.0 * math.pi
K_GAUSS = float(K.K_GAUSS)
LOG_2PI = math.log(2.0 * math.pi)


@njit(cache=True)
def _kepler_E(M, e):
    E = M + 0.85 * e * (1.0 if math.sin(M) >= 0.0 else -1.0)
    for _ in range(60):
        s = math.sin(E)
        c = math.cos(E)
        f = E - e * s - M
        fp = 1.0 - e * c
        dE = -f / (fp - 0.5 * f * e * s / fp)
        E += dE
        if abs(dE) < 1e-12:
            break
    return E


@njit(cache=True)
def positions(epochs, el, plx, mstar, masses, tau_ref):
    """(ra, dec) in mas, arrays (n_planets, n_epochs), with the Jacobi star offset."""
    n = epochs.shape[0]
    npl = el.shape[0]
    ra = np.zeros((npl, n))
    dec = np.zeros((npl, n))
    m_inside = mstar
    for k in range(npl):
        a, e, inc, aop, pan, tau = el[k, 0], el[k, 1], el[k, 2], el[k, 3], el[k, 4], el[k, 5]
        mk = m_inside + masses[k]
        P = TWO_PI * math.sqrt(a * a * a / mk) / K_GAUSS
        co, so = math.cos(aop), math.sin(aop)
        cO, sO = math.cos(pan), math.sin(pan)
        ci = math.cos(inc)
        Ax, Ay = co * sO + so * ci * cO, co * cO - so * ci * sO
        Bx, By = -so * sO + co * ci * cO, -so * cO - co * ci * sO
        sq = math.sqrt(1.0 - e * e)
        for j in range(n):
            frac = (epochs[j] - tau_ref) / P
            frac -= math.floor(frac)
            M = TWO_PI * (frac - tau)
            M -= TWO_PI * math.floor(M / TWO_PI)
            E = _kepler_E(M, e)
            X = a * (math.cos(E) - e)
            Y = a * sq * math.sin(E)
            x = (X * Ax + Y * Bx) * plx
            y = (X * Ay + Y * By) * plx
            for q in range(k):                      # + sum_j (m_j / M_{k-1}) X_j
                w = masses[q] / m_inside
                x += w * ra[q, j]
                y += w * dec[q, j]
            ra[k, j] = x
            dec[k, j] = y
        m_inside = mk
    return ra, dec


@njit(cache=True)
def loglike(el, plx, mstar, masses, tau_ref, epochs, pidx, eidx, is_seppa, obs, err, corr,
            jit_row, off_row, prior_prec, logdet_prior):
    """Gaussian log-likelihood with jitter and analytically marginalised
    per-source (scale, rotation) offsets. ``off_row[i]`` is the offset-source
    index of row i, or -1; ``jit_row[i]`` the jitter (mas) of row i."""
    ra, dec = positions(epochs, el, plx, mstar, masses, tau_ref)
    n_off = prior_prec.shape[0]
    F = np.zeros((n_off, 2, 2))
    b = np.zeros((n_off, 2))
    lnl = 0.0
    for i in range(pidx.shape[0]):
        mra = ra[pidx[i], eidx[i]]
        mdec = dec[pidx[i], eidx[i]]
        if is_seppa[i]:
            msep = math.sqrt(mra * mra + mdec * mdec)
            mpa = math.atan2(mra, mdec)
            r0 = obs[i, 0] - msep
            dpa = math.radians(obs[i, 1]) - mpa + math.pi
            dpa = dpa - TWO_PI * math.floor(dpa / TWO_PI) - math.pi
            r1 = dpa * obs[i, 0]
            a00, a01, a10, a11 = msep, 0.0, 0.0, obs[i, 0]
        else:
            r0 = obs[i, 0] - mra
            r1 = obs[i, 1] - mdec
            a00, a01, a10, a11 = mra, mdec, mdec, -mra
        s0 = math.sqrt(err[i, 0] ** 2 + jit_row[i] ** 2)
        s1 = math.sqrt(err[i, 1] ** 2 + jit_row[i] ** 2)
        c01 = corr[i] * s0 * s1
        det = s0 * s0 * s1 * s1 - c01 * c01
        i00, i11, i01 = s1 * s1 / det, s0 * s0 / det, -c01 / det
        lnl += -0.5 * (i00 * r0 * r0 + i11 * r1 * r1 + 2.0 * i01 * r0 * r1 + math.log(det) + 2.0 * LOG_2PI)
        q = off_row[i]
        if q >= 0:
            # N^-1 A (2x2) and accumulate F = A^T N^-1 A, b = A^T N^-1 r
            n00 = i00 * a00 + i01 * a10
            n01 = i00 * a01 + i01 * a11
            n10 = i01 * a00 + i11 * a10
            n11 = i01 * a01 + i11 * a11
            F[q, 0, 0] += a00 * n00 + a10 * n10
            F[q, 0, 1] += a00 * n01 + a10 * n11
            F[q, 1, 0] += a01 * n00 + a11 * n10
            F[q, 1, 1] += a01 * n01 + a11 * n11
            b[q, 0] += n00 * r0 + n10 * r1
            b[q, 1] += n01 * r0 + n11 * r1
    for q in range(n_off):
        h00 = F[q, 0, 0] + prior_prec[q, 0]
        h11 = F[q, 1, 1] + prior_prec[q, 1]
        h01 = F[q, 0, 1]
        dh = h00 * h11 - h01 * h01
        # b^T H^-1 b
        quad = (h11 * b[q, 0] ** 2 - 2.0 * h01 * b[q, 0] * b[q, 1] + h00 * b[q, 1] ** 2) / dh
        lnl += 0.5 * quad - 0.5 * (logdet_prior[q] + math.log(dh))
    return lnl


@njit(cache=True)
def bound_region(x, y, vx, vy, plx, mu, acap):
    """For a sky state (mas, mas/yr): the largest line-of-sight distance z_max
    (au) of an orbit with semi-major axis below acap, and the projected radius
    rho (au) and sky speed squared (au/yr)^2. z_max < 0 means no such orbit."""
    rho = math.sqrt(x * x + y * y) / plx
    vs2 = (vx * vx + vy * vy) / (plx * plx)
    r_max = 2.0 / (vs2 / mu + 1.0 / acap)
    if r_max <= rho:
        return -1.0, rho, vs2
    return math.sqrt(r_max * r_max - rho * rho), rho, vs2


@njit(cache=True)
def vz_max(z, rho, vs2, mu, acap):
    """Largest |vz| (au/yr) of a bound orbit with a < acap at depth z (au)."""
    r = math.sqrt(rho * rho + z * z)
    w = mu * (2.0 / r - 1.0 / acap) - vs2
    return math.sqrt(w) if w > 0.0 else 0.0


@njit(cache=True)
def fill_depths(theta, u, idx, fit_masses, mjup, G, acap):
    """Conditional prior transform of the non-coplanar models: z ~ U(0, z_max),
    vz ~ U(-vz_max(z), +vz_max(z)), so every draw is a bound orbit with a < acap.
    Fills theta[z], theta[vz] (mas, mas/yr) in place from the unit cube u."""
    plx = theta[idx[35]]
    m_inside = theta[idx[36]]
    for k in range(4):
        base = 8 * k
        mk = m_inside + (theta[idx[base + 6]] * mjup if fit_masses else 0.0)
        mu = G * mk
        zmax, rho, vs2 = bound_region(theta[idx[base]], theta[idx[base + 1]], theta[idx[base + 2]],
                                      theta[idx[base + 3]], plx, mu, acap)
        if zmax <= 0.0:
            theta[idx[base + 4]] = 0.0
            theta[idx[base + 5]] = 0.0
        else:
            z = u[idx[base + 4]] * zmax
            vm = vz_max(z, rho, vs2, mu, acap)
            theta[idx[base + 4]] = z * plx
            theta[idx[base + 5]] = (2.0 * u[idx[base + 5]] - 1.0) * vm * plx
        m_inside = mk


@njit(cache=True)
def loglike_state(theta, idx, kind, fit_masses, mjup, G, tau_ref, t_state, emax, emin, amin, amax,
                  acap, ratio_bounds, log_vbox, jit_idx, row_jit, epochs, pidx, eidx, is_seppa, obs, err,
                  corr, off_row, prior_prec, logdet_prior):
    """ln L + ln(target prior / sampling prior) for the state-vector parametrisation
    (models.OrbitModel; same numbers as its numpy loglike_reference).

    ``idx`` (int array, -1 = absent): per planet k (inner to outer) at 8k + j the
    positions of x, y, vx, vy, z, vz (j = 0..5) and mass (j = 6); [32] cos_i,
    [33] Omega (coplanar models); [35] plx; [36] mstar. ``kind``: 0 unconstrained,
    1 coplanar, 2 near1248, 3 near1248_coplanar.

    States outside the prior's support (unbound, a or e out of range, period
    ratios outside their windows) return a large finite penalty that shrinks
    towards the support, -1e12 (1 + violation), instead of a flat -inf: nested
    sampling then compresses smoothly into the support. These regions have
    zero prior mass, so they add nothing to the evidence.
    """
    coplanar = kind == 1 or kind == 3
    ratios = kind == 2 or kind == 3
    plx = theta[idx[35]]
    mstar = theta[idx[36]]
    if plx <= 0.0 or mstar <= 0.0:
        return -1e300
    masses = np.zeros(4)
    if fit_masses:
        for k in range(4):
            masses[k] = theta[idx[8 * k + 6]] * mjup
    ci = 1.0
    ti = 0.0
    cO = 1.0
    sO = 0.0
    if coplanar:
        ci = theta[idx[32]]
        if ci <= 0.0 or ci > 1.0:
            return -1e300
        ti = math.sqrt(1.0 - ci * ci) / ci
        cO = math.cos(theta[idx[33]])
        sO = math.sin(theta[idx[33]])
    el = np.zeros((4, 6))
    lp = 0.0
    viol = 0.0
    m_inside = mstar
    prevP = 0.0
    log_plx = math.log(plx)
    for k in range(4):
        base = 8 * k
        x = theta[idx[base]]
        y = theta[idx[base + 1]]
        vx = theta[idx[base + 2]]
        vy = theta[idx[base + 3]]
        if coplanar:
            z = ti * (cO * x - sO * y)
            vz = ti * (cO * vx - sO * vy)
        else:
            z = theta[idx[base + 4]]
            vz = theta[idx[base + 5]]
        mk = m_inside + masses[k]
        mu = G * mk
        if not coplanar:
            # (z, vz) were drawn uniformly in the bound region: divide out that density
            zmax, rho, vs2 = bound_region(x, y, vx, vy, plx, mu, acap)
            vm = vz_max(z / plx, rho, vs2, mu, acap) if zmax > 0.0 else 0.0
            if zmax <= 0.0 or vm <= 0.0 or z > zmax * plx * (1 + 1e-9) or abs(vz) > vm * plx * (1 + 1e-9):
                # graded: how far the sky speed exceeds the bound limit at rho
                vlim2 = mu * (2.0 / rho - 1.0 / acap)
                viol += 1.0 + max(0.0, vs2 / vlim2 - 1.0 if vlim2 > 0.0 else 1.0)
                m_inside = mk
                prevP = 0.0
                continue
            lp += math.log(zmax * plx) + math.log(2.0 * vm * plx)
        rx, ry, rz = x / plx, y / plx, z / plx
        ux, uy, uz = vx / plx, vy / plx, vz / plx
        rn = math.sqrt(rx * rx + ry * ry + rz * rz)
        v2 = ux * ux + uy * uy + uz * uz
        inv_a = 2.0 / rn - v2 / mu
        if inv_a <= 0.0:
            # unbound: measure how far, skip this planet's elements
            viol += 1.0 - inv_a * rn
            m_inside = mk
            prevP = 0.0
            continue
        a = 1.0 / inv_a
        hx = ry * uz - rz * uy
        hy = rz * ux - rx * uz
        hz = rx * uy - ry * ux
        hn = math.sqrt(hx * hx + hy * hy + hz * hz)
        ex = (uy * hz - uz * hy) / mu - rx / rn
        ey = (uz * hx - ux * hz) / mu - ry / rn
        ez = (ux * hy - uy * hx) / mu - rz / rn
        e = math.sqrt(ex * ex + ey * ey + ez * ez)
        if e > emax or e < emin:
            viol += abs(e - emax) if e > emax else (emin - e) / emin
            m_inside = mk
            prevP = 0.0
            continue
        inc = math.acos(min(1.0, max(-1.0, -hz / hn)))
        pan = math.atan2(-hy, hx) % TWO_PI
        Nx, Ny = -hy, hx
        Nn = math.sqrt(Nx * Nx + Ny * Ny)
        cosw = (Nx * ex + Ny * ey) / (Nn * e)
        aop = math.acos(min(1.0, max(-1.0, cosw)))
        if ez < 0.0:
            aop = TWO_PI - aop
        cosnu = (ex * rx + ey * ry + ez * rz) / (e * rn)
        nu = math.acos(min(1.0, max(-1.0, cosnu)))
        if rx * ux + ry * uy + rz * uz < 0.0:
            nu = TWO_PI - nu
        E = 2.0 * math.atan(math.sqrt((1.0 - e) / (1.0 + e)) * math.tan(0.5 * nu))
        M = (E - e * math.sin(E)) % TWO_PI
        P = TWO_PI * math.sqrt(a * a * a / mk) / K_GAUSS
        frac = (t_state - tau_ref) / P
        frac -= math.floor(frac)
        tau = frac - M / TWO_PI
        tau -= math.floor(tau)
        el[k, 0] = a
        el[k, 1] = e
        el[k, 2] = inc
        el[k, 3] = aop
        el[k, 4] = pan
        el[k, 5] = tau
        # Wang et al.'s prior in the sampled coordinates
        if ratios and k > 0:
            if prevP <= 0.0:                      # inner planet invalid: no ratio to test
                viol += 1.0
                prevP = P
                m_inside = mk
                continue
            rr = P / prevP
            lo, hi = ratio_bounds[k, 0], ratio_bounds[k, 1]
            if rr < lo or rr > hi:
                viol += (lo - rr if rr < lo else rr - hi) / (hi - lo)
                prevP = P
                m_inside = mk
                continue
            lp += math.log(1.5 * rr / (a * (hi - lo)))
        else:
            if a < amin or a > amax:
                viol += (amin - a) / amin if a < amin else (a - amax) / amax
                prevP = P
                m_inside = mk
                continue
            lp -= math.log(a * math.log(amax / amin))
        lp -= math.log(emax - emin)
        if coplanar:
            lp += (-2.0 * math.log(TWO_PI) - 4.0 * log_plx - 2.0 * math.log(ci) - math.log(mu)
                   - math.log(e) + math.log(2.0) + 0.5 * math.log(1.0 - e * e))
        else:
            lp += (math.log(2.0) - 3.0 * math.log(TWO_PI) - 6.0 * log_plx - 1.5 * math.log(mu)
                   - 0.5 * math.log(a) - math.log(e))
        prevP = P
        m_inside = mk
    if viol > 0.0:
        return -1e12 * (1.0 + viol)
    if coplanar:
        lp -= math.log(TWO_PI)
    n = pidx.shape[0]
    jit_row = np.zeros(n)
    for i in range(n):
        if row_jit[i] >= 0:
            jit_row[i] = 10.0 ** theta[jit_idx[row_jit[i]]]
    val = loglike(el, plx, mstar, masses, tau_ref, epochs, pidx, eidx, is_seppa, obs, err, corr,
                  jit_row, off_row, prior_prec, logdet_prior)
    if not math.isfinite(val):
        return -1e300
    return val + lp + log_vbox
