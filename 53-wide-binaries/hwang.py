"""Eccentricities of individual wide binaries from the v-r angle
(Hwang, Ting & Zakamska 2022), reimplemented.

The v-r angle gamma is the angle between a pair's sky-projected separation
vector and its relative proper motion. Circular orbits are seen near 90 deg,
very eccentric ones near 0 or 180 deg (radial motion), with projection
blurring the relation. With random orientations and phases, p(gamma | e) is
computed by Monte Carlo on a grid e = 0, 0.01, ..., 0.99; a pair's posterior
is

    p(e | gamma_obs) ~ p(e) * integral p(gamma | e) N(gamma_obs; gamma, sigma_gamma) dgamma,

with their separation-dependent prior p(e) = (1 + alpha) e^alpha,
alpha(w) = 1.25 tanh((w - 1.87)/0.88) + 0.12 (w = log10 s/AU; alpha = 0 below
w = 1.78). Their catalogue lists the posterior mode and the 68% highest
posterior density interval; so does ours. We need this to give *mock*
binaries the same kind of eccentricity estimates that Chae's method uses on
real ones, and we check it against their catalogue for real pairs.
"""

from __future__ import annotations

import numpy as np

from config import DERIVED_DIR

E_GRID = np.round(np.arange(0.0, 0.995, 0.01), 3)        # 100 values, 0 ... 0.99
G_EDGES = np.linspace(0.0, 180.0, 181)
G_MID = 0.5 * (G_EDGES[1:] + G_EDGES[:-1])


def alpha_fit(s_au):
    w = np.log10(s_au)
    a = 1.25 * np.tanh((w - 1.87) / 0.88) + 0.12
    return np.where(w < 1.78, 0.0, a)


def vr_angle_grid(n_per_e: int = 400_000, seed: int = 2022):
    """p(gamma | e) per degree, shape (len(E_GRID), 180); cached."""
    path = DERIVED_DIR / "hwang_p_gamma_e.npy"
    if path.exists():
        return np.load(path)
    rng = np.random.default_rng(seed)
    P = np.zeros((len(E_GRID), len(G_MID)))
    for k, e in enumerate(E_GRID):
        M = 2 * np.pi * rng.random(n_per_e)
        E = M + e * np.sin(M)
        for _ in range(60):
            E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        x = np.stack([np.cos(E) - e, np.sqrt(1 - e * e) * np.sin(E), np.zeros(n_per_e)], 1)
        v = np.stack([-np.sin(E), np.sqrt(1 - e * e) * np.cos(E), np.zeros(n_per_e)], 1)
        # isotropic orientation: random rotation, observe along z
        A = rng.standard_normal((n_per_e, 3, 3))
        Q, R = np.linalg.qr(A)
        Q = Q * np.sign(np.diagonal(R, axis1=1, axis2=2))[:, None, :]
        xs = np.einsum("nij,nj->ni", Q, x)[:, :2]
        vs = np.einsum("nij,nj->ni", Q, v)[:, :2]
        c = np.sum(xs * vs, 1) / (np.linalg.norm(xs, axis=1) * np.linalg.norm(vs, axis=1))
        g = np.degrees(np.arccos(np.clip(c, -1, 1)))
        P[k], _ = np.histogram(g, bins=G_EDGES, density=True)
    np.save(path, P)
    return P


def vr_angle(ra1, dec1, ra2, dec2, pmra1, pmdec1, pmra2, pmdec2, epmra1, epmdec1, epmra2, epmdec2):
    """Observed v-r angle (deg), its error (deg) and Delta mu / sigma
    (small-angle tangent plane at star 1; errors without correlations)."""
    dra = (ra2 - ra1) * np.cos(np.radians(dec1))
    ddec = dec2 - dec1
    dmx, dmy = pmra2 - pmra1, pmdec2 - pmdec1
    sx2, sy2 = epmra1**2 + epmra2**2, epmdec1**2 + epmdec2**2
    rr = np.hypot(dra, ddec)
    mm = np.hypot(dmx, dmy)
    c = (dra * dmx + ddec * dmy) / (rr * mm)
    g = np.degrees(np.arccos(np.clip(c, -1, 1)))
    # error on the angle: proper-motion error perpendicular to Delta mu
    ux, uy = -dmy / mm, dmx / mm
    s_perp = np.sqrt(ux**2 * sx2 + uy**2 * sy2)
    s_par = np.sqrt((dmx / mm) ** 2 * sx2 + (dmy / mm) ** 2 * sy2)
    return g, np.degrees(s_perp / mm), mm / s_par


def posterior_summary(g_obs, g_err, s_au, P=None, chunk: int = 20000):
    """Mode and 68% HPD interval of p(e | gamma_obs) for each pair."""
    P = vr_angle_grid() if P is None else P
    n = len(g_obs)
    mode, lo, hi = np.empty(n), np.empty(n), np.empty(n)
    alpha = alpha_fit(s_au)
    de = 0.01
    for a in range(0, n, chunk):
        sl = slice(a, min(a + chunk, n))
        go, ge = g_obs[sl, None], np.maximum(g_err[sl, None], 0.5)
        K = np.exp(-0.5 * ((G_MID[None, :] - go) / ge) ** 2)
        # reflect the Gaussian at 0 and 180 deg so no probability leaks out
        K += np.exp(-0.5 * ((-G_MID[None, :] - go) / ge) ** 2)
        K += np.exp(-0.5 * ((360 - G_MID[None, :] - go) / ge) ** 2)
        L = K @ P.T                                         # (m, n_e)
        e = np.maximum(E_GRID, 1e-4)
        prior = (1 + alpha[sl, None]) * e[None, :] ** alpha[sl, None]
        post = L * prior
        post /= post.sum(1, keepdims=True)
        mode[sl] = E_GRID[np.argmax(post, 1)]
        # HPD: include grid points from the highest density down to 68%
        order = np.argsort(-post, axis=1)
        ps = np.take_along_axis(post, order, 1)
        cum = np.cumsum(ps, 1)
        k68 = np.argmax(cum >= 0.68, axis=1)
        for i in range(post.shape[0]):
            sel = E_GRID[order[i, :k68[i] + 1]]
            lo[a + i], hi[a + i] = sel.min(), min(sel.max() + de, 1.0) if sel.max() >= 0.99 else sel.max()
    return mode, lo, hi
