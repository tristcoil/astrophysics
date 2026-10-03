"""The outer Galaxy seen from above: hydrogen mapped by kinematic distance.

Outside the solar circle a line of sight meets each orbit only once, so if the
gas moves on circular orbits its velocity tells us its distance. Every
brightness measurement T_b(l, b, v) then becomes gas at a known point in
space, and the whole survey can be redrawn as a map of the disc seen from
above. Inside the solar circle this fails, because a line of sight crosses
each orbit twice (near and far), so the map covers only R > R0.

Method (after Levine, Blitz & Heiles 2006, "LBH06", who did this with the
older LAB survey):

1. Walk out along every sightline in steps of 0.1 kpc. At each step, the
   rotation curve predicts the velocity the gas there should have.
2. Read the brightness at that velocity from the spectrum. Brightness per unit
   velocity is column density per unit velocity; multiplying by how fast the
   velocity changes with distance (|dv/ds|) gives the gas density there.
   Where the gas is bright enough to be partly opaque we correct with an
   assumed spin temperature of 155 K (LBH06's choice).
3. Add each step's hydrogen mass to a map cell (x, y) and to a
   (R, azimuth, z) grid, from which the disc's mean height (the warp) and
   thickness (the flare) follow.

Directions within 15 deg of the Galactic centre and anticentre are skipped:
there, orbits cross the line of sight almost perpendicularly, velocities are
tiny, and random motions of a few km/s would move gas by kiloparsecs.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

import config as cfg
import galaxy
import survey

M_H_G = 1.6735575e-24          # hydrogen atom mass, g
MSUN_PER_PC3_IN_G_CM3 = 6.7702e-23
KPC_CM = 3.0857e21
LB_EXCLUDE_DEG = 15.0          # LBH06 exclusion around l = 0 and 180
SPIN_TEMPERATURE_K = 155.0     # LBH06


@dataclass
class OuterSetup:
    theta_fn: object = galaxy.theta_reid19
    obs: galaxy.Observer = field(default_factory=galaxy.Observer.reid19)
    ds_kpc: float = 0.1
    s_max_kpc: float = 45.0
    r_min_kpc: float = cfg.R0_KPC + 0.5
    spin_temperature: float | None = SPIN_TEMPERATURE_K
    #: LBH06's empirical radial velocity v_R = vpi(R) cos(phi) that removes the
    #: jump in density across l = 180; (alpha, beta) of their Appendix, or None
    #: for pure circular orbits.
    elliptic: tuple[float, float] | None = None
    # face-on Cartesian map
    xy_half: float = 32.0
    xy_cell: float = 0.25
    # polar grid for the vertical structure
    r_edges: np.ndarray = field(default_factory=lambda: np.arange(8.0, 36.01, 0.5))
    phi_edges: np.ndarray = field(default_factory=lambda: np.arange(0.0, 360.01, 2.0))
    z_edges: np.ndarray = field(default_factory=lambda: np.arange(-8.0, 8.001, 0.1))
    batch: int = 16


def _velocity_along(l_deg, b_deg, d, st: OuterSetup):
    """Model LSRK velocity along a sightline (in-plane distance d), with the
    optional LBH06 elliptical-orbit term."""
    v = galaxy.v_lsr(l_deg, b_deg, d, st.theta_fn, st.obs)
    if st.elliptic is not None:
        a, bq = st.elliptic
        x, y = galaxy.xy_from_ld(l_deg, d, st.obs.r0)
        R = np.hypot(x, y)
        q = (R - st.obs.r0) / st.obs.r0
        v_pi = a * q + bq * q ** 2
        phi = np.arctan2(x, -y)  # LBH06 azimuth: 0 beyond the centre, 90 on the l = 90 side
        lr = np.radians(l_deg)
        proj = (x * np.sin(lr) - y * np.cos(lr)) / R  # radial unit vector . sightline
        v = v + v_pi * np.cos(phi) * proj * np.cos(np.radians(b_deg))
    return v


def _in_ranges(lon, l_ranges):
    if l_ranges is None:
        return np.ones(lon.shape, dtype=bool)
    return np.any([(lon >= lo) & (lon <= hi) for lo, hi in l_ranges], axis=0)


def build(st: OuterSetup | None = None, *, progress: bool = False,
          l_ranges: list[tuple[float, float]] | None = None, names: list[str] | None = None,
          cache: dict | None = None) -> dict:
    """Map the outer disc. ``l_ranges`` restricts to longitude ranges (deg),
    ``names`` to particular survey cubes; ``cache`` (a dict) keeps the cubes'
    selected columns in memory between calls (for fitting loops)."""
    st = st or OuterSetup()
    s = np.arange(st.ds_kpc / 2, st.s_max_kpc, st.ds_kpc)
    n_xy = int(round(2 * st.xy_half / st.xy_cell))
    nR, nP, nZ = st.r_edges.size - 1, st.phi_edges.size - 1, st.z_edges.size - 1
    xy_m = np.zeros(n_xy * n_xy)
    xy_mz = np.zeros(n_xy * n_xy)
    rpz = np.zeros(nR * nP * nZ)
    pix_sr = np.radians(1 / survey.PIX_PER_DEG) ** 2
    to_msun_pc3 = M_H_G / MSUN_PER_PC3_IN_G_CM3

    for tile in survey.tiles(names or (cfg.HI4PI_ROW_E + cfg.HI4PI_ROW_F_WARP)):
        cols, kglob = survey.own_columns(tile)
        rows = survey.own_rows(tile)
        lon = (kglob / survey.PIX_PER_DEG) % 360.0
        keep = ((np.abs((lon + 180) % 360 - 180) > LB_EXCLUDE_DEG)
                & (np.abs(lon - 180) > LB_EXCLUDE_DEG) & _in_ranges(lon, l_ranges))
        if not keep.any():
            continue
        if cache is not None and tile.name in cache:
            cube = cache[tile.name]
        else:
            cube = np.asarray(tile.data[:, rows[0]:rows[-1] + 1, :][:, :, cols[keep]],
                              dtype=np.float32)                      # (v, b, l)
            cube = np.nan_to_num(cube, nan=0.0)
            if cache is not None:
                cache[tile.name] = cube
        v_axis = tile.vel
        v0, dv, nv = v_axis[0], v_axis[1] - v_axis[0], v_axis.size
        lat = tile.lat[rows[0]:rows[-1] + 1]
        sinb = np.sin(np.radians(lat))[None, :, None]
        cosb = np.cos(np.radians(lat))[None, :, None]
        ls = lon[keep]
        for j0 in range(0, ls.size, st.batch):
            lb = ls[j0:j0 + st.batch][:, None, None]              # (c, 1, 1)
            jj = np.arange(j0, min(j0 + st.batch, ls.size))[:, None, None]
            d = s[None, None, :] * cosb                              # (1, b, s)
            vmod = _velocity_along(lb, lat[None, :, None], d, st)    # (c, b, s)
            dvds = np.abs(np.gradient(vmod, st.ds_kpc, axis=2))
            x, y = galaxy.xy_from_ld(lb, d, st.obs.r0)
            R = np.hypot(x, y)
            z = np.broadcast_to(s[None, None, :] * sinb + 0.0055, R.shape)  # Sun 5.5 pc up
            fi = (vmod - v0) / dv
            inside = (fi >= 0) & (fi <= nv - 1) & (R >= st.r_min_kpc)
            i0 = np.clip(np.floor(fi).astype(int), 0, nv - 2)
            f = np.clip(fi - i0, 0.0, 1.0)
            bidx = np.arange(lat.size)[None, :, None]
            T = cube[i0, bidx, jj] * (1 - f) + cube[i0 + 1, bidx, jj] * f
            # Negative brightness (noise) is kept, not clipped: clipping would
            # add a positive floor that at R = 30 kpc is a third of the signal.
            if st.spin_temperature is not None:
                ts = st.spin_temperature
                Tp = np.clip(T, 0.0, 0.95 * ts)
                T = np.where(T > 0, -ts * np.log1p(-Tp / ts), T)
            n = cfg.C_HI * T * dvds / KPC_CM                         # atoms / cm^3
            vol_pc3 = (s[None, None, :] * 1e3) ** 2 * pix_sr * cosb * (st.ds_kpc * 1e3)
            mass = np.where(inside, n * to_msun_pc3 * vol_pc3, 0.0)  # Msun, can be < 0
            ix = np.floor((x + st.xy_half) / st.xy_cell).astype(int)
            iy = np.floor((y + st.xy_half) / st.xy_cell).astype(int)
            ok = inside & (ix >= 0) & (ix < n_xy) & (iy >= 0) & (iy < n_xy)
            flat = iy[ok] * n_xy + ix[ok]
            xy_m += np.bincount(flat, weights=mass[ok], minlength=xy_m.size)
            xy_mz += np.bincount(flat, weights=(mass * z)[ok], minlength=xy_m.size)
            phi = np.degrees(np.arctan2(x, -y)) % 360.0
            ir = np.searchsorted(st.r_edges, R) - 1
            ip = np.searchsorted(st.phi_edges, phi) - 1
            iz = np.searchsorted(st.z_edges, z) - 1
            ok = (inside & (ir >= 0) & (ir < nR) & (ip >= 0) & (ip < nP)
                  & (iz >= 0) & (iz < nZ))
            flat = (ir[ok] * nP + ip[ok]) * nZ + iz[ok]
            rpz += np.bincount(flat, weights=mass[ok], minlength=rpz.size)
        if progress:
            print(f"    {tile.name} done", flush=True)
    cell_pc2 = (st.xy_cell * 1e3) ** 2
    xy_m = xy_m.reshape(n_xy, n_xy)
    xy_mz = xy_mz.reshape(n_xy, n_xy)
    with np.errstate(invalid="ignore", divide="ignore"):
        zmean = np.where(xy_m > 0, xy_mz / xy_m, np.nan)
    return {
        "xy_sigma": xy_m / cell_pc2,  # Msun / pc^2 (HI only, no helium)
        "xy_zmean": zmean,
        "xy_edges": np.linspace(-st.xy_half, st.xy_half, n_xy + 1),
        "rpz_mass": rpz.reshape(nR, nP, nZ),
        "r_edges": st.r_edges, "phi_edges": st.phi_edges, "z_edges": st.z_edges,
    }


#: The two sectors either side of the anticentre exclusion whose surface
#: densities LBH06 matched to fix the elliptical-orbit correction.
ANTICENTRE_SECTORS = ((155.0, 165.0), (195.0, 205.0))


def anticentre_profiles(prod: dict, r0: float, r_bins: np.ndarray) -> tuple:
    """Mean surface density in radius bins for the two anticentre sectors."""
    e = prod["xy_edges"]
    x = 0.5 * (e[1:] + e[:-1])
    X, Y = np.meshgrid(x, x)
    R = np.hypot(X, Y)
    lcell, _ = galaxy.ld_from_xy(X, Y, r0)
    out = []
    for lo, hi in ANTICENTRE_SECTORS:
        sec = (lcell >= lo) & (lcell <= hi)
        prof = []
        for a, b in zip(r_bins[:-1], r_bins[1:]):
            sel = sec & (R >= a) & (R < b)
            prof.append(prod["xy_sigma"][sel].mean() if sel.sum() >= 2 else np.nan)
        out.append(np.array(prof))
    return tuple(out)


def anticentre_mismatch(sig_a: np.ndarray, sig_b: np.ndarray) -> float:
    """LBH06's modified chi^2: sum over radii of the fractional difference,
    (a - b) / min(a, b), squared."""
    ok = np.isfinite(sig_a) & np.isfinite(sig_b) & (sig_a > 0) & (sig_b > 0)
    frac = (sig_a[ok] - sig_b[ok]) / np.minimum(sig_a[ok], sig_b[ok])
    return float(np.sum(frac ** 2))


def fit_elliptic(st: OuterSetup, *, r_bins: np.ndarray | None = None,
                 start: tuple[float, float] = (8.67, -1.08)) -> dict:
    """Re-derive LBH06's elliptical-orbit correction v_R = v_pi(R) cos(phi),
    v_pi = alpha (R - R0)/R0 + beta ((R - R0)/R0)^2, for our Galactic
    constants, by making the surface density match across the anticentre."""
    from dataclasses import replace

    from scipy.optimize import minimize

    r_bins = np.arange(10.0, 30.01, 1.0) if r_bins is None else r_bins
    names = ["CAR_E08.fits", "CAR_E09.fits", "CAR_E10.fits", "CAR_E11.fits"]
    cache: dict = {}
    trace = []

    def cost(p):
        trial = replace(st, elliptic=(float(p[0]), float(p[1])),
                        r_edges=np.arange(8.0, 12.01, 4.0))  # polar grid not needed here
        prod = build(trial, l_ranges=list(ANTICENTRE_SECTORS), names=names, cache=cache)
        a, b = anticentre_profiles(prod, st.obs.r0, r_bins)
        c = anticentre_mismatch(a, b)
        trace.append((float(p[0]), float(p[1]), c))
        return c

    none = cost((0.0, 0.0))
    res = minimize(cost, np.array(start), method="Nelder-Mead",
                   options={"xatol": 0.1, "fatol": 1e-3, "maxfev": 60})
    return {"alpha": float(res.x[0]), "beta": float(res.x[1]), "chi2": float(res.fun),
            "chi2_circular": none, "n_evals": len(trace), "trace": trace}


def vertical_structure(prod: dict, *, clip_sigma: float = 2.0, iterations: int = 3) -> dict:
    """Surface density, mean height and thickness on the (R, phi) grid.

    The disc is isolated with LBH06's dispersion filter: estimate the mean
    height and dispersion, keep only gas within ``clip_sigma`` dispersions of
    it (dropping clouds and spurs off the disc), and repeat. The window is
    found on a lightly smoothed, non-negative copy of each vertical profile
    (noise can be negative); the sums inside the window use the unclipped
    data, so noise averages out instead of adding a floor.
    """
    m = prod["rpz_mass"]
    ze = prod["z_edges"]
    zc = 0.5 * (ze[1:] + ze[:-1])
    re, pe = prod["r_edges"], prod["phi_edges"]
    area = 0.5 * (re[1:, None] ** 2 - re[:-1, None] ** 2) * np.radians(np.diff(pe))[None, :] * 1e6
    kernel = np.ones(3) / 3
    ms = np.clip(np.apply_along_axis(lambda p: np.convolve(p, kernel, mode="same"), 2, m), 0, None)
    keep = np.ones(m.shape, dtype=bool)
    with np.errstate(invalid="ignore", divide="ignore"):
        for _ in range(iterations):
            w = np.where(keep, ms, 0.0)
            tot = w.sum(axis=2)
            zbar = (w * zc).sum(axis=2) / tot
            disp = np.sqrt((w * (zc - zbar[..., None]) ** 2).sum(axis=2) / tot)
            keep = np.abs(zc - zbar[..., None]) <= clip_sigma * disp[..., None]
        mk = np.where(keep, m, 0.0)
        tot_k = mk.sum(axis=2)
        zbar_k = (mk * zc).sum(axis=2) / tot_k
        var_k = (mk * (zc - zbar_k[..., None]) ** 2).sum(axis=2) / tot_k
        disp_k = np.sqrt(np.where(var_k > 0, var_k, np.nan))
    return {"sigma_all": m.sum(axis=2) / area, "sigma": tot_k / area,
            "zbar": zbar_k, "disp": disp_k, "window": clip_sigma * disp,
            "r_centres": 0.5 * (re[1:] + re[:-1]), "phi_centres": 0.5 * (pe[1:] + pe[:-1])}


def arm_contrast(vs: dict, *, smooth_deg: float = 5.0, reference_deg: float = 40.0) -> np.ndarray:
    """Surface density relative to its surroundings along the ring, on the
    (R, phi) grid: a local arm-contrast map.

    Kinematic distances smear every sightline's noise into a radial streak, but
    spiral arms (pitch angles ~10 deg) run mostly along azimuth. Smoothing along
    azimuth (FWHM ``smooth_deg``) therefore suppresses the streaks and keeps the
    arms; dividing by a much wider azimuthal average (``reference_deg``)
    removes the radial decline and the large-scale lopsidedness, leaving
    arm-scale over- and under-densities.
    """
    from scipy.ndimage import gaussian_filter

    sig = vs["sigma_all"]
    good = np.isfinite(sig) & (sig != 0)
    dphi = np.median(np.diff(vs["phi_centres"]))
    s0, w = np.where(good, sig, 0.0), good.astype(float)

    def smooth(fwhm_deg):
        sig_bins = (0.6, fwhm_deg / 2.3548 / dphi)
        mode = ("nearest", "wrap")
        return gaussian_filter(s0, sig_bins, mode=mode) / np.maximum(
            gaussian_filter(w, sig_bins, mode=mode), 1e-9)

    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(good, smooth(smooth_deg) / smooth(reference_deg), np.nan)


def fit_warp(vs: dict, r_lo: float, r_hi: float, *, min_frac: float = 0.3,
             phi_exclude: float = LB_EXCLUDE_DEG) -> dict:
    """Fit LBH06's warp model to one ring: z(phi) = W0 + W1 sin(phi - phi1)
    + W2 sin(2 phi - phi2), as a linear least-squares problem in the sine and
    cosine terms. The excluded wedges near the Sun-centre line are skipped."""
    rc, pc = vs["r_centres"], vs["phi_centres"]
    sel_r = (rc >= r_lo) & (rc < r_hi)
    z = np.nanmean(vs["zbar"][sel_r], axis=0) * 1e3  # pc
    phi = np.radians(pc)
    # the excluded wedges sit near phi = 0 (beyond the centre) and 180 (anticentre),
    # but in azimuth they are narrower than 15 deg far out; drop what is empty
    good = np.isfinite(z)
    near_line = (np.minimum(np.abs(pc), np.abs(360 - pc)) < 3) | (np.abs(pc - 180) < 3)
    good &= ~near_line
    if good.mean() < min_frac:
        return {"ok": False}
    A = np.vstack([np.ones(good.sum()), np.sin(phi[good]), np.cos(phi[good]),
                   np.sin(2 * phi[good]), np.cos(2 * phi[good])]).T
    coef, *_ = np.linalg.lstsq(A, z[good], rcond=None)
    w0, s1, c1, s2, c2 = coef
    # a sin(phi) + c cos(phi) = W sin(phi - phi1) with W cos(phi1) = a, -W sin(phi1) = c
    w1, phi1 = np.hypot(s1, c1), np.degrees(np.arctan2(-c1, s1)) % 360
    w2, phi2 = np.hypot(s2, c2), np.degrees(np.arctan2(-c2, s2)) % 360
    resid = z[good] - A @ coef
    return {"ok": True, "W0": w0, "W1": w1, "phi1": phi1, "W2": w2, "phi2": phi2,
            "rms": float(np.sqrt(np.mean(resid ** 2))), "phi": pc, "z": z, "good": good}


def lbh06_warp_amplitudes(r_lbh):
    """LBH06 Table 1: W_n = k0 + k1 (R - Rk) + k2 (R - Rk)^2 (pc), for R >= Rk,
    in their scale (R0 = 8.5 kpc, flat 220 km/s rotation)."""
    r = np.asarray(r_lbh, dtype=float)
    out = {}
    for name, (rk, k0, k1, k2) in {"W0": (15, -66, 150, -0.47), "W1": (10, 9, 197, -3.1),
                                   "W2": (15, -70, 171, -5.3)}.items():
        out[name] = np.where(r >= rk, k0 + k1 * (r - rk) + k2 * (r - rk) ** 2, np.nan)
    return out
