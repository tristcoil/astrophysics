"""Where the mass is: the 2MRS galaxy density field and the gravity it exerts.

The 2MASS Redshift Survey counts every galaxy brighter than Ks = 11.75 over
91% of the sky. In the near infrared a galaxy's light traces its old stars and
so, roughly, its mass. Turning the catalogue into a density field takes four
corrections, following Carrick et al. (2015):

1. **Selection.** A flux-limited survey sees only the brightest galaxies far
   away. Each galaxy's luminosity is weighted up by the light expected from
   the fainter galaxies missed at its distance, using a luminosity function
   fitted to the survey itself (Schechter form, STY maximum likelihood).
2. **Bias.** Brighter galaxies cluster more strongly (b/b* = 0.73 + 0.24 L/L*,
   Westover 2007), and far away only bright ones remain, so the field is
   rescaled at each distance to the clustering of a typical (L*) galaxy.
3. **The Milky Way's shadow.** Within 5 deg of the Galactic plane (8 deg toward
   the bulge) 2MRS has no galaxies; the strips just above and below are cloned
   into it (Yahil et al. 1991; Lavaux & Hudson 2011).
4. **Redshift is not distance.** Galaxies fall toward mass, so their redshifts
   misplace them. Fingers of God are collapsed by putting group members at
   their group's mean redshift (Tully 2015 groups), and the rest is iterated:
   compute the velocity field, correct each distance, recompute -- while
   switching gravity on slowly (Yahil et al. 1991).

In linear theory the peculiar velocity field is fixed by the density field:
    v(x) = (beta* H0 / 4 pi) int delta*_g(x') (x' - x) / |x' - x|^3 d^3x',
with beta* = f / b*, the growth rate divided by the bias of L* galaxies. The
integral is the same one that gives Newtonian gravity: galaxies move in the
direction they are pulled. It is evaluated by FFT convolution on a grid padded
to twice its size, so the box has no periodic images.
"""

from __future__ import annotations

import functools
import time

import numpy as np
import pandas as pd
from scipy import fft, ndimage, optimize, special

import catalogs
import config
import cosmology as cg
from config import C_KMS

OM = config.PLANCK18.Om0
H100 = 100.0
#: Absolute magnitude of the Sun in Ks (Vega), Willmer 2018 (ApJS 236, 47): 3.27.
M_SUN_K = 3.27
#: 2MRS magnitudes are isophotal (20 mag/arcsec^2); total magnitudes are
#: brighter by 0.20 +- 0.04 on average (Kochanek et al. 2001). Only used to
#: compare our luminosity density with published totals: every result that
#: divides luminosities by the mean luminosity density is independent of it.
K20_TO_TOTAL = -0.20
#: Galaxies fainter than this are not used to fit the luminosity function
#: (the Schechter form fails among dwarfs, and nearby distances are poor).
LF_M_FAINT = -19.0
LF_CZ_MIN = 1000.0


# --------------------------------------------------------------------------
# The galaxy sample
# --------------------------------------------------------------------------

@functools.lru_cache(maxsize=1)
def galaxies() -> pd.DataFrame:
    """2MRS galaxies with redshifts, group-averaged LG-frame redshifts and completeness weights."""
    g = catalogs.twomrs_with_groups()
    g = g[g["has_cz"] & (g["Kcmag"] <= config.TWOMRS_KLIM)].copy()
    # Fingers of God: members of a Tully (2015) group share the group's mean
    # LG-frame redshift (mean over the 2MRS members we matched).
    grp = g[g["nest"] > 0].groupby("nest")["cz_lg"].agg(["mean", "size"])
    g["cz_group"] = g["cz_lg"]
    has = g["nest"] > 0
    g.loc[has, "cz_group"] = g.loc[has, "nest"].map(grp["mean"]).to_numpy()
    g["group_size"] = 1
    g.loc[has, "group_size"] = g.loc[has, "nest"].map(grp["size"]).to_numpy()
    # The Local Group itself, and anything so close that its redshift says
    # nothing about its distance, is left out.
    g = g[g["cz_group"] > config.TWOMRS_CZ_MIN].copy()
    g["completeness_w"] = completeness_weights(g)
    return g.reset_index(drop=True)


def completeness_weights(g: pd.DataFrame, nside: int = 8) -> np.ndarray:
    """1 / (fraction of 2MRS targets with a redshift) in each 7-degree sky pixel.

    2.4% of the 2MRS targets lack a redshift, mostly near the Galactic plane.
    """
    import healpy as hp

    allg = catalogs.twomrs()
    allg = allg[allg["Kcmag"] <= config.TWOMRS_KLIM]
    pix_all = hp.vec2pix(nside, *allg[["gx", "gy", "gz"]].to_numpy().T)
    n_all = np.bincount(pix_all, minlength=hp.nside2npix(nside))
    n_z = np.bincount(pix_all[allg["has_cz"].to_numpy()], minlength=hp.nside2npix(nside))
    with np.errstate(divide="ignore", invalid="ignore"):
        comp = np.where(n_all > 0, n_z / np.maximum(n_all, 1), 1.0)
    comp = np.clip(comp, 0.5, 1.0)
    pix = hp.vec2pix(nside, *g[["gx", "gy", "gz"]].to_numpy().T)
    return 1.0 / comp[pix]


# --------------------------------------------------------------------------
# Luminosity function, weights, bias
# --------------------------------------------------------------------------

def absolute_magnitude(K, r, z_cos):
    """M - 5 log10 h for apparent Ks magnitude K at comoving distance r (h^-1 Mpc)."""
    dl = (1.0 + np.asarray(z_cos)) * np.asarray(r)
    return np.asarray(K) - (5.0 * np.log10(np.maximum(dl, 1e-3)) + 25.0) - config.KCORR_COEFF * np.log10(1.0 + np.asarray(z_cos))


def limiting_magnitude(r, z_cos, mlim: float = config.TWOMRS_KLIM):
    return absolute_magnitude(mlim, r, z_cos)


def fit_luminosity_function(K, r, z_cos) -> dict:
    """Schechter (alpha, M*) by STY maximum likelihood (Sandage, Tammann & Yahil 1979).

    The probability of galaxy i's absolute magnitude given its distance is
    Phi(M_i) / int_{-inf}^{min(M_lim,i, M_faint)} Phi(M) dM; maximised over
    (alpha, M*). Independent of the (unknown) clustering along each sightline.
    """
    M = absolute_magnitude(K, r, z_cos)
    Mlim = np.minimum(limiting_magnitude(r, z_cos), LF_M_FAINT)
    use = M < LF_M_FAINT
    M, Mlim = M[use], Mlim[use]
    grid = np.linspace(-28.0, LF_M_FAINT, 4001)

    def nll(p):
        alpha, Ms = p
        x = 10 ** (0.4 * (Ms - grid))
        phi = x ** (alpha + 1) * np.exp(-x)
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (phi[1:] + phi[:-1]) * np.diff(grid))])
        norm = np.interp(Mlim, grid, cum)
        xi = 10 ** (0.4 * (Ms - M))
        lphi = (alpha + 1) * np.log(xi) - xi
        return -np.sum(lphi - np.log(norm))

    sol = optimize.minimize(nll, x0=[-1.0, -23.5], method="Nelder-Mead",
                            options={"xatol": 1e-5, "fatol": 1e-6, "maxiter": 2000})
    # Errors from the curvature of the likelihood.
    hess = _numerical_hessian(nll, sol.x, [0.01, 0.01])
    cov = np.linalg.inv(hess)
    return {"alpha": float(sol.x[0]), "Mstar": float(sol.x[1]), "alpha_err": float(np.sqrt(cov[0, 0])),
            "Mstar_err": float(np.sqrt(cov[1, 1])), "n": int(use.sum())}


def _numerical_hessian(f, x, h):
    x = np.asarray(x, float)
    n = len(x)
    H = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            ei, ej = np.eye(n)[i] * h[i], np.eye(n)[j] * h[j]
            H[i, j] = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) / (4 * h[i] * h[j])
    return H


def luminosity_weight(Mlim, lf: dict):
    """Total light per unit observed light at a distance where galaxies fainter than Mlim are missed:
    Gamma(alpha + 2) / Gamma(alpha + 2, L_lim / L*)."""
    x = 10 ** (0.4 * (lf["Mstar"] - np.asarray(Mlim)))
    return 1.0 / special.gammaincc(lf["alpha"] + 2.0, x)


def bias_factor(Mlim, lf: dict):
    """Effective bias of the luminosity-weighted field relative to L* galaxies:
    0.73 + 0.24 Gamma(a+3, x) / Gamma(a+2, x) for b/b* = 0.73 + 0.24 L/L*."""
    a = lf["alpha"]
    x = 10 ** (0.4 * (lf["Mstar"] - np.asarray(Mlim)))
    ratio = (special.gamma(a + 3) * special.gammaincc(a + 3, x)) / (special.gamma(a + 2) * special.gammaincc(a + 2, x))
    return config.BIAS_A + config.BIAS_B * ratio


# --------------------------------------------------------------------------
# The zone of avoidance
# --------------------------------------------------------------------------

def zoa_clones(g: pd.DataFrame, rng_seed: int = 7) -> pd.DataFrame:
    """Fill the masked Galactic plane with copies of the strips beside it.

    For |l| >= 30 deg, galaxies at 5 <= |b| < 10 deg are copied 5 deg toward the
    plane; toward the bulge (|l| < 30 deg), those at 8 <= |b| < 16 deg are copied
    8 deg. Each copy keeps its longitude, redshift and magnitude, and carries
    weight (mask area)/(source area), which is 1.008 and 1.020.
    """
    l = (g["GLON"].to_numpy() + 180.0) % 360.0 - 180.0
    b = g["GLAT"].to_numpy()
    bulge = np.abs(l) < config.ZOA_L_BULGE
    width = np.where(bulge, config.ZOA_B_BULGE, config.ZOA_B)
    src = (np.abs(b) >= width) & (np.abs(b) < 2 * width)
    c = g[src].copy()
    w = width[src]
    bnew = c["GLAT"].to_numpy() - np.sign(c["GLAT"].to_numpy()) * w
    ratio = np.sin(np.radians(w)) / (np.sin(np.radians(2 * w)) - np.sin(np.radians(w)))
    c["GLAT"] = bnew
    n = catalogs.unit_vector(c["GLON"].to_numpy(), bnew)
    c[["gx", "gy", "gz"]] = n
    c[["sgx", "sgy", "sgz"]] = catalogs.gal_to_sg(n)
    c["completeness_w"] = c["completeness_w"] * ratio
    c["clone"] = True
    c["source_index"] = np.flatnonzero(src)
    return c


# --------------------------------------------------------------------------
# Grids, smoothing and linear-theory velocities
# --------------------------------------------------------------------------

@functools.lru_cache(maxsize=2)
def grid_axes(n: int = config.RECON_NGRID, half: float = config.RECON_HALF):
    cell = 2 * half / n
    centres = -half + cell * (np.arange(n) + 0.5)
    return centres, cell


def cic_assign(pos: np.ndarray, weights: np.ndarray, n: int = config.RECON_NGRID, half: float = config.RECON_HALF):
    """Cloud-in-cell assignment of weighted points (N, 3) to an n^3 grid spanning [-half, half]."""
    _, cell = grid_axes(n, half)
    x = (pos + half) / cell - 0.5
    i0 = np.floor(x).astype(int)
    f = x - i0
    out = np.zeros((n, n, n))
    for dx in (0, 1):
        wx = f[:, 0] if dx else 1 - f[:, 0]
        for dy in (0, 1):
            wy = f[:, 1] if dy else 1 - f[:, 1]
            for dz in (0, 1):
                wz = f[:, 2] if dz else 1 - f[:, 2]
                ii, jj, kk = i0[:, 0] + dx, i0[:, 1] + dy, i0[:, 2] + dz
                ok = (ii >= 0) & (ii < n) & (jj >= 0) & (jj < n) & (kk >= 0) & (kk < n)
                np.add.at(out, (ii[ok], jj[ok], kk[ok]), (weights * wx * wy * wz)[ok])
    return out


def gaussian_smooth(field: np.ndarray, sigma: float, cell: float) -> np.ndarray:
    """Gaussian smoothing in Fourier space (zero-padded, so nothing wraps around)."""
    n = field.shape[0]
    big = np.zeros((2 * n,) * 3, dtype=np.float32)
    big[:n, :n, :n] = field
    k = 2 * np.pi * fft.fftfreq(2 * n, d=cell)
    kz = 2 * np.pi * fft.rfftfreq(2 * n, d=cell)
    F = fft.rfftn(big, workers=4)
    F *= np.exp(-0.5 * sigma**2 * k[:, None, None] ** 2)
    F *= np.exp(-0.5 * sigma**2 * k[None, :, None] ** 2)
    F *= np.exp(-0.5 * sigma**2 * kz[None, None, :] ** 2)
    return fft.irfftn(F, s=big.shape, workers=4)[:n, :n, :n].astype(np.float64)


@functools.lru_cache(maxsize=2)
def _kernel_ffts(n: int, cell: float):
    """FFTs of the three components of x/|x|^3 on the doubled grid (Hockney's method)."""
    m = 2 * n
    idx = np.arange(m)
    idx = np.where(idx < n, idx, idx - m) * cell
    X, Y, Z = np.meshgrid(idx, idx, idx, indexing="ij")
    r3 = (X**2 + Y**2 + Z**2) ** 1.5
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = np.where(r3 > 0, 1.0 / r3, 0.0)
    out = []
    for comp in (X, Y, Z):
        out.append(fft.rfftn((comp * inv).astype(np.float32), workers=4))
    del X, Y, Z, r3, inv
    return tuple(out)


def velocity_field(delta: np.ndarray, beta: float, cell: float) -> np.ndarray:
    """Linear-theory velocity (3, n, n, n) in km/s from a density contrast grid.

    v(x) = beta H0 / (4 pi) sum_cells delta(x') (x' - x) / |x' - x|^3 dV.
    The kernel is evaluated at cell centres (the self-cell contributes zero),
    which is accurate for a field smoothed over more than a cell.
    """
    n = delta.shape[0]
    kx, ky, kz = _kernel_ffts(n, cell)
    big = np.zeros((2 * n,) * 3, dtype=np.float32)
    big[:n, :n, :n] = delta
    D = fft.rfftn(big, workers=4)
    out = np.empty((3, n, n, n))
    # (delta * K)(x) = sum delta(x') K(x - x'); with K(y) = y/|y|^3 this is the
    # pull *toward* x from x', i.e. the negative of what we need: v points from
    # x toward the mass at x'.
    pref = -beta * H100 / (4 * np.pi) * cell**3
    for a, K in enumerate((kx, ky, kz)):
        out[a] = pref * fft.irfftn(D * K, s=big.shape, workers=4)[:n, :n, :n]
    return out


def sample_field(field: np.ndarray, pos: np.ndarray, half: float = config.RECON_HALF) -> np.ndarray:
    """Trilinear interpolation of grid field(s) at positions (N, 3) in h^-1 Mpc."""
    n = field.shape[-1]
    _, cell = grid_axes(n, half)
    coords = ((np.asarray(pos) + half) / cell - 0.5).T
    if field.ndim == 3:
        return ndimage.map_coordinates(field, coords, order=1, mode="nearest")
    return np.stack([ndimage.map_coordinates(f, coords, order=1, mode="nearest") for f in field], axis=-1)


# --------------------------------------------------------------------------
# The reconstruction
# --------------------------------------------------------------------------

def _radius_grid(n=config.RECON_NGRID, half=config.RECON_HALF):
    c, _ = grid_axes(n, half)
    X, Y, Z = np.meshgrid(c, c, c, indexing="ij")
    return np.sqrt(X**2 + Y**2 + Z**2)


def density_contrast(pos, lum, comp_w, r, lf: dict, *, rmax: float = config.RECON_RMAX,
                     smooth: float = config.RECON_SMOOTH, return_parts: bool = False):
    """Bias-normalised, smoothed luminosity-density contrast delta*_g on the grid.

    Each galaxy contributes L x (luminosity weight at its distance) x
    (completeness weight). The mean is taken over the sphere r < rmax; outside
    it the contrast is set to zero (the mean density).
    """
    z = cg.z_from_comoving(r, OM)
    Mlim = limiting_magnitude(r, z)
    wl = luminosity_weight(Mlim, lf)
    keep = r < rmax
    rho = cic_assign(pos[keep], (lum * wl * comp_w)[keep])
    _, cell = grid_axes()
    R = _radius_grid()
    inside = R < rmax
    rho_mean = rho[inside].sum() / (4.0 / 3.0 * np.pi * rmax**3)  # light per unit volume
    delta = np.where(inside, rho / (rho_mean * cell**3) - 1.0, 0.0)
    zc = cg.z_from_comoving(R, OM)
    bfac = bias_factor(limiting_magnitude(np.maximum(R, 0.5), zc), lf)
    delta_star = np.where(inside, delta / bfac, 0.0)
    sm = gaussian_smooth(delta_star, smooth, cell)
    if return_parts:
        return sm, {"rho_mean": rho_mean, "delta_raw": delta, "bias_grid_example": float(np.median(bfac[inside]))}
    return sm


def galaxy_luminosity(K):
    """Ks luminosity (h^-2 Lsun, isophotal) at unit distance-modulus bookkeeping:
    L = 10^(-0.4 (M - M_sun)); with M computed per galaxy in ``reconstruct``."""
    return np.asarray(K)


def reconstruct(*, betas=config.RECON_BETAS, verbose: bool = True) -> dict:
    """Iterative real-space reconstruction of 2MRS for a ramp of beta* values.

    Returns the luminosity function, and for each beta* the averaged galaxy
    distances (2MRS galaxies, ZoA clones inherit their source's distance) and
    the LG velocity it predicts. Fields are recomputed on demand from the
    distances (``field_for``). Cached on disk.
    """
    path = config.DERIVED_DIR / "reconstruction.npz"
    g = galaxies()
    if path.exists():
        d = np.load(path, allow_pickle=True)
        if len(d["r0"]) == len(g) and np.allclose(d["betas"], betas):
            return {k: d[k] for k in d.files} | {"lf": d["lf"].item()}
    n_hat = g[["gx", "gy", "gz"]].to_numpy()
    z_obs = g["cz_group"].to_numpy() / C_KMS
    r0 = cg.comoving_distance(z_obs, OM)
    lf = fit_luminosity_function(g["Kcmag"].to_numpy()[g["cz_group"] > LF_CZ_MIN],
                                 r0[g["cz_group"] > LF_CZ_MIN], z_obs[g["cz_group"] > LF_CZ_MIN])
    if verbose:
        print(f"    luminosity function: alpha {lf['alpha']:.3f} +- {lf['alpha_err']:.3f}, "
              f"M* - 5 log h = {lf['Mstar']:.3f} +- {lf['Mstar_err']:.3f} ({lf['n']:,} galaxies)")
    clones = zoa_clones(g)
    src = clones["source_index"].to_numpy()
    n_all = np.vstack([n_hat, clones[["gx", "gy", "gz"]].to_numpy()])
    K_all = np.concatenate([g["Kcmag"].to_numpy(), clones["Kcmag"].to_numpy()])
    cw_all = np.concatenate([g["completeness_w"].to_numpy(), clones["completeness_w"].to_numpy()])
    z_all = np.concatenate([z_obs, z_obs[src]])
    _, cell = grid_axes()
    r = r0.copy()
    history = []
    out_r, out_vlg = [], []
    t0 = time.time()
    for it, beta in enumerate(betas):
        r_all = np.concatenate([r, r[src]])
        pos = n_all * r_all[:, None]
        zc = cg.z_from_comoving(r_all, OM)
        lum = 10 ** (-0.4 * (absolute_magnitude(K_all, r_all, zc) - M_SUN_K))
        delta = density_contrast(pos, lum, cw_all, r_all, lf)
        if beta > 0:
            v = velocity_field(delta, beta, cell)
            v0 = sample_field(v, np.zeros((1, 3)))[0]
            vg = sample_field(v, n_hat * r[:, None])
            u = np.sum((vg - v0) * n_hat, axis=1)
            z_cos = (1.0 + z_obs) / (1.0 + u / C_KMS) - 1.0
            r_new = cg.comoving_distance(np.maximum(z_cos, 1e-5), OM)
        else:
            v0 = np.zeros(3)
            r_new = r0.copy()
        history.append(r_new)
        history = history[-config.RECON_AVERAGE_LAST:]
        r = np.mean(history, axis=0)
        out_r.append(r.astype(np.float32))
        out_vlg.append(v0)
        if verbose and (it % 5 == 0 or it == len(betas) - 1):
            print(f"    beta* = {beta:.2f}: v(LG) = {np.linalg.norm(v0):6.0f} km/s, "
                  f"median |dr| {np.median(np.abs(r_new - r0)):.2f} h^-1 Mpc ({time.time() - t0:.0f} s)", flush=True)
    res = {"betas": np.array(betas), "r": np.array(out_r), "v_lg": np.array(out_vlg), "r0": r0,
           "lf": lf, "n_clones": len(clones)}
    np.savez(path, **{k: v for k, v in res.items() if k != "lf"}, lf=np.array(lf, dtype=object))
    return res


def field_for(rec: dict, beta: float, *, smooth: float = config.RECON_SMOOTH, return_parts: bool = False):
    """The smoothed delta*_g grid for the reconstruction step closest to ``beta``."""
    g = galaxies()
    i = int(np.argmin(np.abs(rec["betas"] - beta)))
    r = rec["r"][i].astype(float)
    clones = zoa_clones(g)
    src = clones["source_index"].to_numpy()
    n_all = np.vstack([g[["gx", "gy", "gz"]].to_numpy(), clones[["gx", "gy", "gz"]].to_numpy()])
    K_all = np.concatenate([g["Kcmag"].to_numpy(), clones["Kcmag"].to_numpy()])
    cw_all = np.concatenate([g["completeness_w"].to_numpy(), clones["completeness_w"].to_numpy()])
    r_all = np.concatenate([r, r[src]])
    zc = cg.z_from_comoving(r_all, OM)
    lum = 10 ** (-0.4 * (absolute_magnitude(K_all, r_all, zc) - M_SUN_K))
    return density_contrast(n_all * r_all[:, None], lum, cw_all, r_all, rec["lf"], smooth=smooth,
                            return_parts=return_parts)
