"""The three-dimensional shape of the Magellanic Clouds from Cepheid distances.

A Cepheid's offset from the Leavitt law (in the reddening-free Wesenheit
index) is mostly a distance offset: 0.1 mag brighter = 4.7% nearer. Across
the face of the LMC these offsets reveal its tilted disk.

Geometry (our own derivation; same form as van der Marel & Cioni 2001, eq. 8).
Put the galaxy centre at distance D0 in the direction (alpha0, delta0). A star
at angular distance rho from the centre and position angle PA (north through
east) lies in a plane disk, inclined by i to the sky with its line of nodes at
position angle Theta, if its distance is

    D / D0 = cos i / (cos i cos rho - sin i sin rho sin(PA - Theta)).

With this sign the half of the disk around PA = Theta + 90 deg is the far side
and PA = Theta - 90 deg the near side.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares


def sky_offsets(ra, dec, ra0, dec0):
    """Orthographic offsets (east, north) in radians, the cosine of rho, and PA (rad, N through E)."""
    a, d = np.radians(ra), np.radians(dec)
    a0, d0 = np.radians(ra0), np.radians(dec0)
    x = np.cos(d) * np.sin(a - a0)
    y = np.sin(d) * np.cos(d0) - np.cos(d) * np.sin(d0) * np.cos(a - a0)
    z = np.cos(d) * np.cos(d0) * np.cos(a - a0) + np.sin(d) * np.sin(d0)
    rho = np.arccos(np.clip(z, -1, 1))
    pa = np.arctan2(x, y) % (2 * np.pi)
    return x, y, z, rho, pa


def disk_distance_ratio(rho, pa, incl_deg, theta_deg):
    i, th = np.radians(incl_deg), np.radians(theta_deg)
    return np.cos(i) / (np.cos(i) * np.cos(rho) - np.sin(i) * np.sin(rho) * np.sin(pa - th))


def fit_disk(ra, dec, dmu, ra0, dec0, *, nsigma=3.0, start=(25.0, 150.0), iters=20):
    """Fit an inclined plane to distance-modulus offsets ``dmu`` (relative to the centre).

    The model is dmu = mu_off + 5 log10(D/D0)(i, Theta); mu_off absorbs any
    zero-point offset. Iteratively clipped at ``nsigma``. Returns the fit,
    its covariance-based errors and the kept mask.
    """
    _, _, _, rho, pa = sky_offsets(ra, dec, ra0, dec0)
    ok = np.isfinite(dmu)
    mask = ok.copy()
    p = np.array([0.0, *start])

    def resid(par, m):
        off, inc, th = par
        return dmu[m] - (off + 5 * np.log10(disk_distance_ratio(rho[m], pa[m], inc, th)))

    for _ in range(iters):
        res = least_squares(resid, p, args=(mask,), bounds=([-1, 0, -360], [1, 80, 720]))
        p = res.x
        r = dmu - (p[0] + 5 * np.log10(disk_distance_ratio(rho, pa, p[1], p[2])))
        s = np.std(r[mask], ddof=3)
        new = ok & (np.abs(r) < nsigma * s)
        if np.array_equal(new, mask):
            break
        mask = new
    J = res.jac
    cov = np.linalg.inv(J.T @ J) * s ** 2
    err = np.sqrt(np.diag(cov))
    return {"offset": p[0], "incl": p[1], "theta": p[2] % 360, "e_offset": err[0], "e_incl": err[1],
            "e_theta": err[2], "sigma": s, "n": int(mask.sum()), "mask": mask, "resid": r,
            "near_side_pa": (p[2] - 90) % 360}


def bootstrap_disk(ra, dec, dmu, ra0, dec0, n_boot=200, seed=0, **kw):
    """Bootstrap spread of (incl, theta): stars resampled with replacement."""
    rng = np.random.default_rng(seed)
    out = []
    n = len(ra)
    for _ in range(n_boot):
        k = rng.integers(0, n, n)
        f = fit_disk(ra[k], dec[k], dmu[k], ra0, dec0, **kw)
        out.append((f["incl"], f["theta"]))
    out = np.array(out)
    return {"incl_std": float(np.std(out[:, 0])), "theta_std": float(np.std(out[:, 1]))}


def cartesian(ra, dec, dist, ra0, dec0, d0):
    """Galaxy-centred Cartesian coordinates (kpc): x east, y north, z away from us."""
    x, y, z, _, _ = sky_offsets(ra, dec, ra0, dec0)
    return dist * x, dist * y, dist * z - d0


def distance_from_mu(mu):
    """Distance in kpc from a distance modulus."""
    return 10 ** (0.2 * np.asarray(mu) + 1) / 1e3


def mu_from_distance(d_kpc):
    return 5 * np.log10(np.asarray(d_kpc) * 1e3) - 5
