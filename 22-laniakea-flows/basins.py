"""Basins of attraction: which way does each patch of the universe fall?

In linear theory the peculiar velocity is the downhill gradient of a velocity
potential psi, with del^2 psi = beta* H0 delta*_g (the continuity equation):
v = -grad psi. Every point, followed downhill, ends at a minimum of psi (an
attractor), and the set of points that end at the same minimum is a *basin of
attraction*. Tully et al. (2014) named the one containing the Milky Way
Laniakea. Basins are the catchment areas of a landscape, so they are found
here exactly like river catchments on a terrain map: each grid cell points to
its lowest neighbour, and following the pointers leads every cell to a pit.

Two choices matter and are varied: the smoothing of the potential (small
smoothing splits basins at every noise ripple; large smoothing merges them), and
whether the uniform flow from mass beyond the survey is included. Tully et al.
used the "local" field, with the tidal flow from distant structures removed,
which is what the 2MRS internal field is.
"""

from __future__ import annotations

import numpy as np
from scipy import fft

import catalogs
import config
import cosmology as cg
import density as dn


def potential(delta: np.ndarray, beta: float, cell: float, *, extra_smooth: float = 0.0,
              vext_gal: np.ndarray | None = None) -> np.ndarray:
    """Velocity potential psi (km/s x h^-1 Mpc) with v = -grad psi, no periodic images.

    psi(x) = -(beta H0 / 4 pi) int delta(x') / |x - x'| d^3x'; a uniform external
    flow adds -V_ext . x.
    """
    if extra_smooth > 0:
        delta = dn.gaussian_smooth(delta, extra_smooth, cell)
    n = delta.shape[0]
    m = 2 * n
    idx = np.arange(m)
    idx = np.where(idx < n, idx, idx - m) * cell
    X, Y, Z = np.meshgrid(idx, idx, idx, indexing="ij")
    rr = np.sqrt(X**2 + Y**2 + Z**2)
    with np.errstate(divide="ignore"):
        kern = np.where(rr > 0, 1.0 / rr, 0.0)
    # The self-cell: the mean of 1/r over a cube of side cell is ~2.38/cell.
    kern[0, 0, 0] = 2.38 / cell
    del X, Y, Z, rr
    big = np.zeros((m, m, m), dtype=np.float32)
    big[:n, :n, :n] = delta
    psi = fft.irfftn(fft.rfftn(big, workers=4) * fft.rfftn(kern.astype(np.float32), workers=4), s=big.shape,
                     workers=4)[:n, :n, :n]
    psi = -beta * 100.0 / (4 * np.pi) * cell**3 * psi.astype(np.float64)
    if vext_gal is not None:
        c, _ = dn.grid_axes(n)
        Xg, Yg, Zg = np.meshgrid(c, c, c, indexing="ij")
        psi = psi - (vext_gal[0] * Xg + vext_gal[1] * Yg + vext_gal[2] * Zg)
    return psi


def watershed(psi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Steepest-descent catchment basins on a 3D grid (26 neighbours).

    Returns (labels, minima): labels[i,j,k] is the flat index of the local
    minimum that cell drains to; minima lists the flat indices of all pits.
    """
    n = psi.shape
    flat = psi.ravel()
    best_val = flat.copy()
    best_idx = np.arange(flat.size)
    ii, jj, kk = np.meshgrid(np.arange(n[0]), np.arange(n[1]), np.arange(n[2]), indexing="ij")
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                if dx == dy == dz == 0:
                    continue
                i2, j2, k2 = ii + dx, jj + dy, kk + dz
                ok = (i2 >= 0) & (i2 < n[0]) & (j2 >= 0) & (j2 < n[1]) & (k2 >= 0) & (k2 < n[2])
                nb = np.where(ok, np.ravel_multi_index((np.clip(i2, 0, n[0] - 1), np.clip(j2, 0, n[1] - 1),
                                                         np.clip(k2, 0, n[2] - 1)), n), 0).ravel()
                # Downhill steps are compared by slope, so diagonal neighbours are not favoured.
                dist = np.sqrt(dx * dx + dy * dy + dz * dz)
                val = np.where(ok.ravel(), flat[nb], np.inf)
                slope = (val - flat) / dist
                cur = (best_val - flat) / np.where(best_idx == np.arange(flat.size), 1.0,
                                                   _dist_of(best_idx, n, flat.size))
                better = slope < cur
                best_val = np.where(better, val, best_val)
                best_idx = np.where(better, nb, best_idx)
    ptr = best_idx
    for _ in range(64):
        nxt = ptr[ptr]
        if np.array_equal(nxt, ptr):
            break
        ptr = nxt
    minima = np.unique(ptr)
    return ptr.reshape(n), minima


def _dist_of(best_idx, shape, size):
    i0, j0, k0 = np.unravel_index(np.arange(size), shape)
    i1, j1, k1 = np.unravel_index(best_idx, shape)
    d = np.sqrt((i1 - i0) ** 2.0 + (j1 - j0) ** 2.0 + (k1 - k0) ** 2.0)
    return np.where(d > 0, d, 1.0)


def basin_report(labels: np.ndarray, psi: np.ndarray, *, beta: float, rmax: float = config.RECON_RMAX,
                 h_for_mass: float = 0.75) -> dict:
    """Properties of the basin containing the origin, and of the other large basins."""
    n = labels.shape[0]
    c, cell = dn.grid_axes(n)
    X, Y, Z = np.meshgrid(c, c, c, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2)
    origin = labels[n // 2, n // 2, n // 2]
    # Average the eight cells around the origin (the origin sits on a cell corner).
    corner = labels[n // 2 - 1:n // 2 + 1, n // 2 - 1:n // 2 + 1, n // 2 - 1:n // 2 + 1].ravel()
    vals, counts = np.unique(corner, return_counts=True)
    ours = int(vals[np.argmax(counts)])
    inside = R < rmax
    lab_in = labels[inside]
    uniq, cnt = np.unique(lab_in, return_counts=True)
    order = np.argsort(cnt)[::-1]
    big = []
    for u, k in zip(uniq[order][:12], cnt[order][:12]):
        pos = np.array([X.ravel()[u], Y.ravel()[u], Z.ravel()[u]])
        sg = catalogs.gal_to_sg(pos)
        l, b = catalogs.lonlat(pos)
        touches = bool(np.any((labels == u) & (R > rmax - 10)))
        big.append({"label": int(u), "volume": float(k * cell**3), "attractor_gal": pos.tolist(),
                    "attractor_sg": sg.tolist(), "attractor_r": float(np.linalg.norm(pos)),
                    "attractor_lb": [float(l), float(b)], "touches_edge": touches, "is_ours": int(u) == ours})
    mask = (labels == ours) & inside
    vol = float(mask.sum() * cell**3)
    pos = np.array([X.ravel()[ours], Y.ravel()[ours], Z.ravel()[ours]])
    l, b = catalogs.lonlat(pos)
    f = cg.growth_rate(config.PLANCK18.Om0)
    rho_crit_h = 2.775e11  # h^2 Msun / Mpc^3, i.e. Msun per (h^-1 Mpc)^3 times h^-1
    mass = config.PLANCK18.Om0 * rho_crit_h * vol / h_for_mass
    # Extent along each supergalactic axis.
    pts = np.stack([X[mask], Y[mask], Z[mask]], axis=1)
    sgp = catalogs.gal_to_sg(pts)
    return {
        "n_basins_inside": int(len(uniq)), "largest": big, "ours_label": ours,
        "ours_attractor_gal": pos.tolist(), "ours_attractor_sg": catalogs.gal_to_sg(pos).tolist(),
        "ours_attractor_lb": [float(l), float(b)], "ours_attractor_r": float(np.linalg.norm(pos)),
        "ours_volume": vol, "ours_equiv_diameter": float(2 * (3 * vol / (4 * np.pi)) ** (1 / 3)),
        "ours_mass_msun": mass, "ours_mass_h": h_for_mass,
        "ours_touches_edge": bool(np.any(mask & (R > rmax - 10))),
        "ours_sg_extent": [[float(sgp[:, i].min()), float(sgp[:, i].max())] for i in range(3)],
        "growth_rate_f": f, "bias_bstar": f / beta,
    }


def landmark_membership(labels: np.ndarray, ours: int) -> dict:
    """Which named clusters lie in our basin."""
    import flows

    lm = flows.landmarks()
    out = {}
    for _, row in lm.iterrows():
        pos = catalogs.unit_vector(row["l"], row["b"]) * row["r"]
        lab = dn.sample_field(labels.astype(float), pos[None, :])  # nearest-ish; labels are integers
        c, cell = dn.grid_axes(labels.shape[0])
        idx = np.clip(np.round((pos + config.RECON_HALF) / cell - 0.5).astype(int), 0, labels.shape[0] - 1)
        out[row["name"]] = bool(labels[idx[0], idx[1], idx[2]] == ours)
    return out


def all_basins(rec: dict, beta: float, vext_gal: np.ndarray) -> dict:
    """Basins at each smoothing scale, with and without the external flow."""
    delta = dn.field_for(rec, beta)
    _, cell = dn.grid_axes()
    out = {}
    for sm in config.BASIN_SMOOTH:
        extra = np.sqrt(max(sm**2 - config.RECON_SMOOTH**2, 0.0))
        for with_ext in (False, True):
            psi = potential(delta, beta, cell, extra_smooth=extra, vext_gal=vext_gal if with_ext else None)
            labels, minima = watershed(psi)
            rep = basin_report(labels, psi, beta=beta)
            rep["members"] = landmark_membership(labels, rep["ours_label"])
            key = f"{sm:g}{'_ext' if with_ext else ''}"
            out[key] = {"report": rep, "labels": labels, "psi": psi}
    return out
