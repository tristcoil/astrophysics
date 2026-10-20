"""The ISW sky that LCDM predicts from 1.2 million BOSS galaxies.

Galaxies trace matter, matter sets the gravitational potential, and the
potential's decay sets the ISW signal. So from a 3D map of galaxies we can
*predict* the ISW temperature pattern on the sky -- no CMB involved. Stacking
that prediction at the void positions gives the signal LCDM expects from each
catalogue, measured with exactly the same filters as the real CMB.

The steps
---------
1. **Footprint.** BOSS describes where it looked as 33,072 spherical polygons
   (the mangle format: each polygon is an intersection of caps on the sphere)
   with a completeness weight -- the fraction of targets that got a redshift.
   We rasterise them onto a HEALPix grid.
2. **Expected counts.** n(z), the mean number density at each redshift, comes
   from the galaxies themselves; times the completeness, it says how many
   galaxies each cell of a 3D grid *should* hold if the universe were smooth.
3. **Overdensity.** delta_g = (galaxies in cell) / (expected) - 1, with BOSS's
   systematics and close-pair weights. Unobserved cells are set to the mean
   (delta = 0): no information.
4. **Bias.** Galaxies cluster more strongly than matter: delta_g = b delta_m.
   We measure b from the galaxies' own power spectrum against CAMB's linear
   matter spectrum for D3's cosmology, separately for LOWZ (z < 0.43) and CMASS.
   Redshift-space distortions (galaxies falling into structures look squashed
   along the line of sight) add f/3 to the effective bias of the monopole.
5. **Potential.** delta_0 = delta_g / [D (b + f/3)] is the matter overdensity
   extrapolated to today; Poisson's equation in Fourier space gives
   Phi_0 = -(3/2) Omega_m (H0/c)^2 delta_0 / k^2.
6. **ISW sky.** For every direction, integrate K(chi) Phi_0 along the line of
   sight (cosmology.isw_kernel). In a matter-only universe K = 0 and the
   predicted map is identically zero.
"""

from __future__ import annotations

import time

import numpy as np

import cmb
import cosmology
from config import (BIAS_KMAX, BIAS_KMIN, BOSS_MASK_NSIDE, BOSS_ZMAX, BOSS_ZMIN, BOSS_ZSPLIT,
                    DERIVED_DIR, GRID_CELL, GRID_MIN_COVERAGE, GRID_SHAPE, GRID_SUBSAMPLE,
                    GRID_Z_CENTRE, ISW_NSIDE, ISW_STEP_MPC_H, SMOOTHING_MPC_H, product)

GAL_CACHE = DERIVED_DIR / "boss_galaxies.npz"
FOOT_CACHE = DERIVED_DIR / "boss_footprint.npz"
GRID_CACHE = DERIVED_DIR / "boss_grid.npz"
BIAS_CACHE = DERIVED_DIR / "boss_bias.npz"
ISW_CACHE = DERIVED_DIR / "boss_isw_{tag}.npy"

CAPS = ("north", "south")
DZ_NZ = 0.005  # redshift bin of n(z)
#: Regions with their own n(z). BOSS's first observing chunks used a faulty
#: LOWZ target selection; the replacement samples (LOWZE2 in chunk 2, LOWZE3 in
#: chunks 3-6, Reid et al. 2016) are sparser at z < 0.43 -- BOSS's own NZ column
#: gives them separate n(z) too. Left uncorrected they would look like
#: underdensities hundreds of Mpc across.
GROUPS = ("north", "north_lowze2", "north_lowze3", "south")


def _group(cap_index, chunk):
    chunk = np.asarray(chunk)
    g = np.where(cap_index == 0, 0, 3)
    g = np.where((cap_index == 0) & (chunk == 2), 1, g)
    g = np.where((cap_index == 0) & (chunk >= 3) & (chunk <= 6), 2, g)
    return g.astype(np.int8)


def _hp():
    import healpy as hp
    return hp


def _log(msg, t0=None):
    extra = f"  ({time.monotonic() - t0:.1f} s)" if t0 is not None else ""
    print(f"    {msg}{extra}", flush=True)


# --------------------------------------------------------------------------
# 1. Galaxies and footprint
# --------------------------------------------------------------------------


def galaxies() -> dict:
    """RA, Dec, z, total weight and squared weight of every galaxy, both caps."""
    if GAL_CACHE.exists():
        d = np.load(GAL_CACHE)
        return {k: d[k] for k in d.files}
    from astropy.io import fits

    out = {k: [] for k in ("ra", "dec", "z", "w", "wfkp", "nz_boss", "cap", "group")}
    for ic, cap in enumerate(CAPS):
        with fits.open(product(f"boss_{cap}").path) as h:
            d = h[1].data
            z = np.asarray(d["Z"], float)
            sel = (z >= BOSS_ZMIN) & (z < BOSS_ZMAX)
            # BOSS's recommended total weight (Reid et al. 2016, Anderson et al.
            # 2014): imaging systematics x (close pairs + redshift failures - 1).
            w = d["WEIGHT_SYSTOT"] * (d["WEIGHT_CP"] + d["WEIGHT_NOZ"] - 1.0)
            out["ra"].append(np.asarray(d["RA"], float)[sel])
            out["dec"].append(np.asarray(d["DEC"], float)[sel])
            out["z"].append(z[sel])
            out["w"].append(np.asarray(w, float)[sel])
            out["wfkp"].append(np.asarray(d["WEIGHT_FKP"], float)[sel])
            out["nz_boss"].append(np.asarray(d["NZ"], float)[sel])
            out["cap"].append(np.full(sel.sum(), ic, np.int8))
            out["group"].append(_group(ic, np.asarray(d["ICHUNK"])[sel]))
    out = {k: np.concatenate(v) for k, v in out.items()}
    np.savez(GAL_CACHE, **out)
    return out


def _mangle_simple_pixel(ra, dec, res: int):
    """Pixel number in mangle's 'simple' scheme (equal-area rows from the north pole).

    Checked against the PIXEL column of all 33,072 BOSS polygons: 99.99% of
    polygon midpoints land in their stated pixel (the rest sit on a boundary).
    """
    n = 2**res
    z = np.sin(np.radians(dec))
    i = np.clip(np.floor((1 - z) * n / 2).astype(np.int64), 0, n - 1)
    j = np.clip(np.floor((np.asarray(ra) % 360) / 360 * n).astype(np.int64), 0, n - 1)
    return (4**res - 1) // 3 + n * i + j


def footprint() -> dict:
    """Completeness of each HEALPix pixel (equatorial, Nside 1024) in each cap.

    A point is inside a mangle polygon if it is inside all of the polygon's
    caps; cap (x, cm) contains point p if 1 - p.x < cm (cm > 0), or the
    complement if cm < 0. The North and South files hold the same 33,072
    polygons with complementary weights, so the geometry is rasterised once.
    """
    if FOOT_CACHE.exists():
        d = np.load(FOOT_CACHE)
        return {k: d[k] for k in d.files}
    from astropy.io import fits

    hp = _hp()
    t0 = time.monotonic()
    tabs = [fits.open(product(f"boss_mask_{cap}").path)[1].data for cap in CAPS]
    poly = tabs[0]
    weights = np.stack([np.asarray(t["WEIGHT"], float) for t in tabs])
    assert np.allclose(tabs[0]["CMCAPS"], tabs[1]["CMCAPS"]), "mask files differ in geometry"
    used = weights.max(axis=0) > 0
    npix = hp.nside2npix(BOSS_MASK_NSIDE)
    # At Nside 1024 whole rings and columns of HEALPix pixel centres fall
    # *exactly* on the boundaries of mangle's cells (and thus on polygon edges),
    # where rounding makes both neighbouring polygons reject the point: one-pixel
    # cracks along every boundary, 2% of the area. A random nudge of ~0.02"
    # puts each such point on one side or the other.
    vec = np.asarray(hp.pix2vec(BOSS_MASK_NSIDE, np.arange(npix))).T
    vec += np.random.default_rng(1).normal(0.0, 1e-7, vec.shape)
    vec /= np.linalg.norm(vec, axis=1, keepdims=True)
    dec = np.degrees(np.arcsin(np.clip(vec[:, 2], -1, 1)))
    ra = np.degrees(np.arctan2(vec[:, 1], vec[:, 0])) % 360
    res = 6  # "pixelization 6s" in the file header
    mpix = _mangle_simple_pixel(ra, dec, res)
    order = np.argsort(mpix, kind="stable")
    sorted_m = mpix[order]
    poly_id = np.full(npix, -1, np.int32)
    ppix = np.asarray(poly["PIXEL"])
    xcaps = np.asarray(poly["XCAPS"], float)
    cmcaps = np.asarray(poly["CMCAPS"], float)
    ncaps = np.asarray(poly["NCAPS"])
    usecaps = np.asarray(poly["USE_CAPS"]).astype(np.int64)
    for mp in np.unique(ppix[used]):
        lo, hi = np.searchsorted(sorted_m, [mp, mp + 1])
        cand = order[lo:hi]
        if len(cand) == 0:
            continue
        v = vec[cand]
        for k in np.flatnonzero((ppix == mp) & used):
            inside = np.ones(len(cand), bool)
            for c in range(ncaps[k]):
                if not (usecaps[k] >> c) & 1:
                    continue
                cd = 1.0 - v @ xcaps[k, c]
                cm = cmcaps[k, c]
                inside &= (cd < cm) if cm >= 0 else (cd > -cm)
            poly_id[cand[inside & (poly_id[cand] < 0)]] = k
    comp = np.zeros((2, npix), np.float32)
    ok = poly_id >= 0
    comp[:, ok] = weights[:, poly_id[ok]]
    chunk = np.zeros(npix, np.int32)
    chunk[ok] = np.asarray(poly["CHUNK"])[poly_id[ok]]
    cap_of = np.where(comp[0] > 0, 0, 1)
    group = np.where(ok & (comp.max(axis=0) > 0), _group(cap_of, chunk), -1).astype(np.int8)
    pix_area = hp.nside2pixarea(BOSS_MASK_NSIDE, degrees=True)
    str_area = {cap: float((poly["STR"] * (weights[i] > 0)).sum() * (180 / np.pi) ** 2)
                for i, cap in enumerate(CAPS)}
    for i, cap in enumerate(CAPS):
        _log(f"{cap}: rasterised area {np.sum(comp[i] > 0) * pix_area:.0f} deg^2 "
             f"(polygon areas sum to {str_area[cap]:.0f}); effective "
             f"{comp[i].sum() * pix_area:.0f} deg^2", t0)
    out = {"comp": comp, "group": group, "pix_area_deg2": np.array(pix_area)}
    np.savez(FOOT_CACHE, **out)
    return out


def nz_table(gal: dict, foot: dict) -> dict:
    """Mean comoving number density n(z) per region, in (h/Mpc)^3, from the galaxies."""
    hp = _hp()
    edges = np.arange(BOSS_ZMIN, BOSS_ZMAX + DZ_NZ / 2, DZ_NZ)
    chi = cosmology.LCDM.chi(edges)
    shell = 4 / 3 * np.pi * np.diff(chi**3)
    sr = hp.nside2pixarea(BOSS_MASK_NSIDE)
    comp = foot["comp"].max(axis=0)
    out = {"edges": edges}

    def density(groups):
        sel_pix = np.isin(foot["group"], groups)
        frac_sky = comp[sel_pix].sum() * sr / (4 * np.pi)  # completeness-weighted
        sel = np.isin(gal["group"], groups)
        counts, _ = np.histogram(gal["z"][sel], edges, weights=gal["w"][sel])
        return counts / (shell * frac_sky), frac_sky * 4 * np.pi * (180 / np.pi) ** 2

    # CMASS (z >= 0.43) was targeted uniformly over the whole north: one n(z)
    # for the cap, as in BOSS's NZ column. Only LOWZ gets region-specific n(z).
    north_all, _ = density([0, 1, 2])
    zc = 0.5 * (edges[1:] + edges[:-1])
    for ig, grp in enumerate(GROUPS):
        own, area = density([ig])
        if grp != "south":
            own = np.where(zc < BOSS_ZSPLIT, own, north_all)
        out[grp] = own
        out[f"area_{grp}"] = area
    return out


# --------------------------------------------------------------------------
# 2-3. The 3D grid
# --------------------------------------------------------------------------


def grid_axes():
    nx, ny, nz = GRID_SHAPE
    x0 = np.array([-nx / 2 * GRID_CELL, -ny / 2 * GRID_CELL, GRID_Z_CENTRE - nz / 2 * GRID_CELL])
    return x0


def radec_to_xyz(ra, dec, chi):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.stack([chi * np.cos(dec) * np.cos(ra), chi * np.cos(dec) * np.sin(ra),
                     chi * np.sin(dec)], axis=-1)


def build_grid() -> dict:
    """Galaxy counts and expected counts on the grid, split at z = 0.43 (LOWZ/CMASS)."""
    if GRID_CACHE.exists():
        d = np.load(GRID_CACHE)
        return {k: d[k] for k in d.files}
    hp = _hp()
    t0 = time.monotonic()
    gal, foot = galaxies(), footprint()
    nzt = nz_table(gal, foot)
    L = cosmology.LCDM
    x0 = grid_axes()
    shape = GRID_SHAPE
    ncell = int(np.prod(shape))

    # galaxies -> cells (nearest grid point), split by redshift
    xyz = radec_to_xyz(gal["ra"], gal["dec"], L.chi(gal["z"]))
    idx = np.floor((xyz - x0) / GRID_CELL).astype(np.int64)
    inbox = np.all((idx >= 0) & (idx < np.array(shape)), axis=1)
    assert inbox.all(), "galaxies outside the box"
    flat = np.ravel_multi_index(idx.T, shape)
    hi = gal["z"] >= BOSS_ZSPLIT
    out = {}
    for name, sel in (("lo", ~hi), ("hi", hi)):
        out[f"ng_{name}"] = np.bincount(flat[sel], gal["w"][sel], ncell).astype(np.float32)
        out[f"w2_{name}"] = np.bincount(flat[sel], gal["w"][sel] ** 2, ncell).astype(np.float32)
    _log(f"{len(gal['z']):,} galaxies placed on a {shape} grid of {GRID_CELL} Mpc/h cells", t0)

    # expected counts by sub-sampling each candidate cell
    nexp = {"lo": np.zeros(ncell, np.float32), "hi": np.zeros(ncell, np.float32)}
    cover = np.zeros(ncell, np.float32)
    chi_lo, chi_hi = L.chi(BOSS_ZMIN), L.chi(BOSS_ZMAX)
    zc = 0.5 * (nzt["edges"][1:] + nzt["edges"][:-1])
    s = GRID_SUBSAMPLE
    offs = (np.stack(np.meshgrid(*[(np.arange(s) + 0.5) / s] * 3, indexing="ij"), -1)
            .reshape(-1, 3) * GRID_CELL)
    vcell = GRID_CELL**3
    ix, iy = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]), indexing="ij")
    comp_all = foot["comp"].max(axis=0)
    for kz in range(shape[2]):
        corner = np.stack([x0[0] + ix * GRID_CELL, x0[1] + iy * GRID_CELL,
                           np.full(ix.shape, x0[2] + kz * GRID_CELL)], -1).reshape(-1, 3)
        r = np.linalg.norm(corner + GRID_CELL / 2, axis=1)
        cand = np.flatnonzero((r > chi_lo - GRID_CELL) & (r < chi_hi + GRID_CELL))
        if len(cand) == 0:
            continue
        pts = (corner[cand, None, :] + offs[None]).reshape(-1, 3)
        rr = np.linalg.norm(pts, axis=1)
        zz = L.z_of_chi(rr)
        dec = np.degrees(np.arcsin(pts[:, 2] / rr))
        ra = np.degrees(np.arctan2(pts[:, 1], pts[:, 0])) % 360
        pix = hp.ang2pix(BOSS_MASK_NSIDE, ra, dec, lonlat=True)
        inz = (zz >= BOSS_ZMIN) & (zz < BOSS_ZMAX)
        jz = np.clip(np.searchsorted(nzt["edges"], zz) - 1, 0, len(zc) - 1)
        c = comp_all[pix]
        grp = foot["group"][pix]
        dens = np.zeros(len(pts))
        for ig, name in enumerate(GROUPS):
            sel_g = grp == ig
            if sel_g.any():
                dens[sel_g] = c[sel_g] * np.interp(zz[sel_g], zc, nzt[name])
        dens *= inz
        seen = (c > 0) & inz
        npt = len(offs)
        flatc = np.ravel_multi_index((ix.ravel()[cand], iy.ravel()[cand],
                                      np.full(len(cand), kz)), shape)
        lo_mask = zz < BOSS_ZSPLIT
        for name, m in (("lo", lo_mask), ("hi", ~lo_mask)):
            nexp[name][flatc] = (dens * m).reshape(-1, npt).mean(axis=1) * vcell
        cover[flatc] = seen.reshape(-1, npt).mean(axis=1)
    out["nexp_lo"], out["nexp_hi"], out["cover"] = nexp["lo"], nexp["hi"], cover
    tot_exp = nexp["lo"].sum() + nexp["hi"].sum()
    tot_gal = out["ng_lo"].sum() + out["ng_hi"].sum()
    _log(f"expected {tot_exp:,.0f} vs weighted galaxies {tot_gal:,.0f}; "
         f"{np.sum(cover >= GRID_MIN_COVERAGE):,} cells observed", t0)
    np.savez(GRID_CACHE, **out)
    return out


def overdensity(grid: dict, which: str = "all") -> np.ndarray:
    """delta_g on the grid (0 where unobserved). which: 'all', 'lo' or 'hi'."""
    if which == "all":
        ng = grid["ng_lo"] + grid["ng_hi"]
        ne = grid["nexp_lo"] + grid["nexp_hi"]
    else:
        ng, ne = grid[f"ng_{which}"], grid[f"nexp_{which}"]
    ok = (grid["cover"] >= GRID_MIN_COVERAGE) & (ne > 0)
    delta = np.zeros_like(ng)
    delta[ok] = ng[ok] / ne[ok] - 1.0
    return delta.reshape(GRID_SHAPE)


# --------------------------------------------------------------------------
# 4. Bias from the power spectrum
# --------------------------------------------------------------------------


def _kgrid():
    nx, ny, nz = GRID_SHAPE
    kx = 2 * np.pi * np.fft.fftfreq(nx, d=GRID_CELL)
    ky = 2 * np.pi * np.fft.fftfreq(ny, d=GRID_CELL)
    kz = 2 * np.pi * np.fft.rfftfreq(nz, d=GRID_CELL)
    return kx, ky, kz


def power_spectrum(grid: dict, which: str, p0: float = 1e4, dk: float = 0.005):
    """FKP power spectrum monopole of one redshift sample (Feldman, Kaiser & Peacock 1994).

    F = w (n_g - n_exp) / sqrt(I2) with weights w = 1/(1 + n P0); P(k) = |F(k)|^2
    minus the shot noise, divided by the nearest-grid-point assignment window.
    """
    import scipy.fft as sfft

    vcell = GRID_CELL**3
    ng, ne, w2 = grid[f"ng_{which}"], grid[f"nexp_{which}"], grid[f"w2_{which}"]
    ok = (grid["cover"] >= GRID_MIN_COVERAGE) & (ne > 0)
    nbar = np.where(ok, ne / vcell, 0.0)
    wf = np.where(ok, 1.0 / (1.0 + nbar * p0), 0.0)
    i2 = np.sum(wf**2 * nbar**2) * vcell
    shot = np.sum(wf**2 * w2 * ok) / i2
    field = (wf * np.where(ok, ng - ne, 0.0) / vcell / np.sqrt(i2)).reshape(GRID_SHAPE)
    fk = sfft.rfftn(field.astype(np.float32), workers=4) * vcell
    kx, ky, kz = _kgrid()
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij", sparse=True)
    kk = np.sqrt(KX**2 + KY**2 + KZ**2)
    win = (np.sinc(KX * GRID_CELL / (2 * np.pi)) * np.sinc(KY * GRID_CELL / (2 * np.pi))
           * np.sinc(KZ * GRID_CELL / (2 * np.pi))) ** 2
    pk = (np.abs(fk) ** 2 - shot) / win
    # rfft stores half the modes; weight the kz>0 planes twice
    mult = np.full(kz.shape, 2.0)
    mult[0] = 1.0
    if GRID_SHAPE[2] % 2 == 0:
        mult[-1] = 1.0
    mult = np.broadcast_to(mult, pk.shape)
    edges = np.arange(dk / 2, 0.3, dk)
    ib = np.digitize(kk.ravel(), edges)
    nb = len(edges) + 1
    psum = np.bincount(ib, (pk * mult).ravel(), nb)
    nmod = np.bincount(ib, mult.ravel(), nb)
    ksum = np.bincount(ib, (kk * mult).ravel(), nb)
    sl = slice(1, len(edges))
    zeff_w = wf**2 * nbar**2
    return {"k": ksum[sl] / nmod[sl], "pk": psum[sl] / nmod[sl], "nmodes": nmod[sl],
            "shot": shot, "i2": i2, "zeff_weights_sum": float(zeff_w.sum())}


def _zeff(grid, which):
    """Effective redshift of a sample: mean z weighted by (w n)^2 over cells."""
    vcell = GRID_CELL**3
    ne = grid[f"nexp_{which}"]
    ok = (grid["cover"] >= GRID_MIN_COVERAGE) & (ne > 0)
    nbar = ne / vcell
    wt = np.where(ok, (nbar / (1 + nbar * 1e4)) ** 2, 0.0)
    x0 = grid_axes()
    idx = np.flatnonzero(wt > 0)
    i, j, k = np.unravel_index(idx, GRID_SHAPE)
    pos = x0 + (np.stack([i, j, k], -1) + 0.5) * GRID_CELL
    z = cosmology.LCDM.z_of_chi(np.linalg.norm(pos, axis=1))
    return float(np.sum(wt[idx] * z) / np.sum(wt[idx]))


def measure_bias(*, force: bool = False) -> dict:
    """Linear bias of LOWZ and CMASS from their power spectra (Kaiser monopole)."""
    if BIAS_CACHE.exists() and not force:
        d = np.load(BIAS_CACHE, allow_pickle=True)
        return d["result"].item()
    grid = build_grid()
    L = cosmology.LCDM
    result = {}
    for which in ("lo", "hi"):
        ps = power_spectrum(grid, which)
        zeff = _zeff(grid, which)
        k, z, plin, sigma8 = cosmology.linear_power([zeff])
        iz = int(np.argmin(np.abs(z - zeff)))
        pl = np.interp(ps["k"], k, plin[iz])
        f = float(L.f(zeff))
        sel = (ps["k"] >= BIAS_KMIN) & (ps["k"] <= BIAS_KMAX)
        # Gaussian errors: var(P) = 2 (P + shot)^2 / N_modes, inflated by 1/f_vol
        # (the survey fills only part of the box, so modes are correlated).
        ratio = ps["pk"][sel] / pl[sel]
        err = np.sqrt(2.0 / ps["nmodes"][sel]) * (ps["pk"][sel] + ps["shot"]) / pl[sel]
        fill = np.mean(grid["cover"] >= GRID_MIN_COVERAGE)
        err /= np.sqrt(fill)
        wts = 1 / err**2
        amp = float(np.sum(wts * ratio) / np.sum(wts))  # b^2 + 2bf/3 + f^2/5
        amp_err = float(1 / np.sqrt(np.sum(wts)))
        b = (-2 * f / 3 + np.sqrt(4 * f**2 / 9 - 4 * (f**2 / 5 - amp))) / 2
        db = amp_err / (2 * b + 2 * f / 3)
        chi2 = float(np.sum(wts * (ratio - amp) ** 2))
        result[which] = {"zeff": zeff, "b": float(b), "b_err": float(db), "f": f,
                         "sigma8_0": sigma8, "sigma8_z": sigma8 * float(L.D(zeff)),
                         "b_sigma8": float(b) * sigma8 * float(L.D(zeff)),
                         "kaiser_amp": amp, "chi2": chi2, "nbins": int(sel.sum()),
                         "k": ps["k"], "pk": ps["pk"], "plin": pl, "shot": ps["shot"]}
        _log(f"{'LOWZ ' if which == 'lo' else 'CMASS'} z_eff = {zeff:.3f}: b = {b:.3f} +- {db:.3f} "
             f"(b sigma8(z) = {result[which]['b_sigma8']:.3f}; chi2 {chi2:.1f}/{sel.sum() - 1})")
    np.savez(BIAS_CACHE, result=np.array(result, dtype=object))
    return result


# --------------------------------------------------------------------------
# 5-6. Potential and the ISW sky
# --------------------------------------------------------------------------


def matter_field(grid, bias: dict, *, kaiser: bool = True, bias_scale: float = 1.0):
    """delta_0: the matter overdensity extrapolated to today, on the grid.

    Each sample's bias is assumed to evolve as b(z) D(z) = const (bias that
    keeps b sigma8(z) fixed, as measured for BOSS); the Kaiser term adds f(z)/3.
    """
    L = cosmology.LCDM
    delta = overdensity(grid).ravel()
    x0 = grid_axes()
    idx = np.flatnonzero(delta != 0)
    i, j, k = np.unravel_index(idx, GRID_SHAPE)
    pos = x0 + (np.stack([i, j, k], -1) + 0.5) * GRID_CELL
    z = L.z_of_chi(np.linalg.norm(pos, axis=1))
    D = L.D(z)
    bD = np.where(z < BOSS_ZSPLIT,
                  bias["lo"]["b"] * L.D(bias["lo"]["zeff"]),
                  bias["hi"]["b"] * L.D(bias["hi"]["zeff"])) * bias_scale
    denom = bD + (L.f(z) * D / 3 if kaiser else 0.0)
    out = np.zeros_like(delta)
    out[idx] = delta[idx] / denom
    return out.reshape(GRID_SHAPE)


def response_and_noise(grid, bias: dict, *, kaiser: bool = True, bias_scale: float = 1.0):
    """Per-cell response B (delta_g = B delta_0 + noise) and noise variance, flattened.

    B = D(z) [b(z) + f(z)/3] with b(z) D(z) held at each sample's measured value;
    the noise of a cell is Poisson in its *expected* count: var(delta_g) =
    kappa / N_exp, where kappa = sum(w^2) / sum(w) accounts for the weights.
    Unobserved cells get B = 0 (no data).
    """
    L = cosmology.LCDM
    ne = grid["nexp_lo"] + grid["nexp_hi"]
    ok = (grid["cover"] >= GRID_MIN_COVERAGE) & (ne > 0)
    idx = np.flatnonzero(ok)
    x0 = grid_axes()
    i, j, k = np.unravel_index(idx, GRID_SHAPE)
    pos = x0 + (np.stack([i, j, k], -1) + 0.5) * GRID_CELL
    z = L.z_of_chi(np.linalg.norm(pos, axis=1))
    bD = np.where(z < BOSS_ZSPLIT,
                  bias["lo"]["b"] * L.D(bias["lo"]["zeff"]),
                  bias["hi"]["b"] * L.D(bias["hi"]["zeff"])) * bias_scale
    resp = np.zeros(ne.shape, np.float64)
    resp[idx] = bD + (L.f(z) * L.D(z) / 3 if kaiser else 0.0)
    gal = galaxies()
    kappa = float(np.sum(gal["w"] ** 2) / np.sum(gal["w"]))
    noise = np.full(ne.shape, np.inf)
    noise[idx] = kappa / ne[idx]
    return resp, noise, kappa


def wiener_filter(grid, bias: dict, *, kaiser: bool = True, bias_scale: float = 1.0,
                  tol: float = 1e-4, maxiter: int = 300):
    """Wiener-filtered matter field delta_0 (today) given the galaxy overdensities.

    The galaxy field is d = B delta_0 + n: a response B per cell (bias, growth,
    redshift-space distortions; zero outside the survey) and Poisson noise n of
    variance N per cell. With the signal covariance S = P_lin(k) from CAMB, the
    minimum-variance estimate of delta_0 -- its conditional mean given d -- is
        delta_WF = (S^-1 + B N^-1 B)^-1 B N^-1 d.
    Writing delta = S^(1/2) x turns this into the well-conditioned system
        (1 + S^(1/2) W S^(1/2)) x = S^(1/2) B d / N,    W = B^2 / N,
    solved by conjugate gradients with FFTs (S is diagonal in Fourier space,
    W in configuration space). The filter suppresses modes the galaxies cannot
    measure -- above all the largest scales, where Poisson noise would
    otherwise dominate the potential (Phi ~ delta / k^2) -- and fills the
    unobserved sky with the field's conditional expectation, which fades to
    zero away from the survey.
    """
    import scipy.fft as sfft

    t0 = time.monotonic()
    vcell = GRID_CELL**3
    resp, noise, kappa = response_and_noise(grid, bias, kaiser=kaiser, bias_scale=bias_scale)
    d = overdensity(grid).ravel().astype(np.float64)
    w = np.where(np.isfinite(noise), resp**2 / noise, 0.0).reshape(GRID_SHAPE).astype(np.float32)
    rhs_real = np.where(np.isfinite(noise), resp * d / noise, 0.0).reshape(GRID_SHAPE)
    k, _, plin, _ = cosmology.linear_power([0.0])
    kx, ky, kz = _kgrid()
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij", sparse=True)
    kk = np.sqrt(KX**2 + KY**2 + KZ**2)
    pk = np.interp(np.log(np.maximum(kk, 1e-6)), np.log(k), plin[0], left=0.0)
    sqrt_s = np.sqrt(pk / vcell).astype(np.float32)  # per-cell normalisation
    sqrt_s[0, 0, 0] = 0.0

    def ssqrt(x):
        return sfft.irfftn(sfft.rfftn(x, workers=4) * sqrt_s, s=GRID_SHAPE, workers=4)

    def apply_a(x):
        return x + ssqrt(w * ssqrt(x))

    b = ssqrt(rhs_real.astype(np.float32)).astype(np.float32)
    x = np.zeros_like(b)
    r = b.copy()
    p = r.copy()
    rr = float(np.vdot(r, r))
    bnorm = np.sqrt(float(np.vdot(b, b)))
    for it in range(maxiter):
        ap = apply_a(p).astype(np.float32)
        alpha = rr / float(np.vdot(p, ap))
        x += alpha * p
        r -= alpha * ap
        rr_new = float(np.vdot(r, r))
        if np.sqrt(rr_new) < tol * bnorm:
            break
        p = r + (rr_new / rr) * p
        rr = rr_new
    _log(f"Wiener filter: CG converged in {it + 1} iterations "
         f"(residual {np.sqrt(rr_new) / bnorm:.1e}; noise factor kappa = {kappa:.3f})", t0)
    return ssqrt(x).astype(np.float32)


def potential(delta0, smoothing: float = SMOOTHING_MPC_H, cosmo=cosmology.LCDM):
    """Phi_0 = -(3/2) Omega_m (H0/c)^2 delta_0 / k^2, Gaussian-smoothed, via FFT."""
    import scipy.fft as sfft

    dk = sfft.rfftn(delta0.astype(np.float32), workers=4)
    kx, ky, kz = _kgrid()
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij", sparse=True)
    k2 = KX**2 + KY**2 + KZ**2
    with np.errstate(divide="ignore", invalid="ignore"):
        fac = -cosmo.poisson_factor / k2 * np.exp(-0.5 * k2 * smoothing**2)
    fac[0, 0, 0] = 0.0
    dk *= fac.astype(np.complex64)
    return sfft.irfftn(dk, s=GRID_SHAPE, workers=4).astype(np.float32)


def isw_map(phi0, nside: int = ISW_NSIDE, cosmo=cosmology.LCDM, chunk: int = 4000) -> np.ndarray:
    """Predicted ISW temperature (uK) in Galactic HEALPix pixels, by ray tracing Phi_0."""
    from scipy.ndimage import map_coordinates

    hp = _hp()
    x0 = grid_axes()
    chi_end = float(cosmology.LCDM.chi(BOSS_ZMAX)) + 400.0
    chi = np.arange(ISW_STEP_MPC_H / 2, chi_end, ISW_STEP_MPC_H)
    kern = cosmo.isw_kernel(cosmo.z_of_chi(chi)) * ISW_STEP_MPC_H * cmb.T_CMB_UK
    npix = hp.nside2npix(nside)
    vg = np.asarray(hp.pix2vec(nside, np.arange(npix)))
    rot = hp.Rotator(coord=["G", "C"])
    ve = np.asarray(rot(vg)).T  # equatorial unit vectors
    out = np.zeros(npix, np.float32)
    for s in range(0, npix, chunk):
        v = ve[s:s + chunk]
        pts = v[:, None, :] * chi[None, :, None]
        coords = ((pts - x0) / GRID_CELL - 0.5).reshape(-1, 3).T
        vals = map_coordinates(phi0, coords, order=1, mode="grid-wrap").reshape(len(v), len(chi))
        out[s:s + chunk] = vals @ kern
    return out


def predicted_isw(tag: str = "fiducial", *, force: bool = False, **kw) -> np.ndarray:
    """Cached predicted ISW map for one set of choices.

    tags: 'fiducial'; 'smooth5', 'smooth20' (smoothing); 'nokaiser'; 'raw20' (no Wiener filter);
    'bias_lo', 'bias_hi' (bias -/+ 10%).
    """
    path = DERIVED_DIR / f"boss_isw_{tag}.npy"
    if path.exists() and not force:
        return np.load(path)
    opts = {"fiducial": {}, "smooth5": {"smoothing": 5.0}, "smooth20": {"smoothing": 20.0},
            "nokaiser": {"kaiser": False}, "bias_lo": {"bias_scale": 0.9},
            "bias_hi": {"bias_scale": 1.1},
            # no Wiener filter: raw delta_g / B, smoothed -- noisy on large scales
            "raw20": {"raw": True, "smoothing": 20.0}}[tag]
    t0 = time.monotonic()
    grid = build_grid()
    bias = measure_bias()
    if opts.get("raw"):
        d0 = matter_field(grid, bias)
    else:
        d0 = wiener_filter(grid, bias, kaiser=opts.get("kaiser", True),
                           bias_scale=opts.get("bias_scale", 1.0))
    del grid
    phi = potential(d0, smoothing=opts.get("smoothing", SMOOTHING_MPC_H))
    del d0
    m = isw_map(phi)
    _log(f"ISW map '{tag}': rms {m.std():.2f} uK, range {m.min():.1f} .. {m.max():.1f} uK", t0)
    np.save(path, m)
    return m
