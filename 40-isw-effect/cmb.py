"""The CMB maps we search for the ISW signal, and the simulated skies that tell us
what pure chance looks like.

The data
--------
* **Planck SMICA** (2018): the CMB with the Galaxy and other foregrounds removed
  by combining Planck's nine frequency channels. We average the two
  half-mission maps (the same sky observed twice), which is the full-mission
  sky; noise is irrelevant here (12.8 uK per 1.7' pixel becomes < 0.1 uK over a
  degree-sized filter).
* **WMAP 9-year ILC**: a different satellite, different instrument, different
  foreground method -- but the same sky. Agreement between the two tests the
  foreground cleaning, not the statistics: the CMB itself is the noise in this
  problem, and both satellites see the same noise.

The simulations
---------------
A simulated sky is a random Gaussian field with Planck's best-fit power spectrum
-- statistically identical to the real CMB, but knowing nothing about where the
supervoids are. Stacking simulated skies at the real void positions shows how
big a "signal" chance alone produces. Simulation number i is always drawn from
seed RANDOM_SEED + i, so every module that needs sky i regenerates the same
sky; that lets us measure how correlated two different ISW measurements are.
"""

from __future__ import annotations

import time

import numpy as np

from config import (COSMO, DERIVED_DIR, LMAX_SIM, NSIDE, RANDOM_SEED, SMICA_FWHM_ARCMIN,
                    product)

K_TO_UK = 1e6
CMB_CACHE = DERIVED_DIR / "cmb_maps.npz"
#: At the cross-correlation resolution (Nside 64, 55' pixels) a pixel is kept if
#: at least this fraction of its Nside-512 children survive the mask; requiring
#: all of them would delete every pixel holding one of the mask's thousands of
#: small point-source holes.
LOWRES_KEEP = 0.8


def _hp():
    import healpy as hp
    return hp


def _log(msg, t0=None):
    extra = f"  ({time.monotonic() - t0:.1f} s)" if t0 is not None else ""
    print(f"    {msg}{extra}", flush=True)


def degrade(m, mask, nside_out, keep_fraction):
    """Average a masked map to lower resolution.

    Each output pixel is the mean of its unmasked children; it is kept only if
    at least ``keep_fraction`` of its children are unmasked. Returns
    (map, boolean mask).
    """
    hp = _hp()
    w = mask.astype(np.float64)
    num = hp.ud_grade(np.where(mask, m, 0.0).astype(np.float64), nside_out)
    frac = hp.ud_grade(w, nside_out)
    keep = frac >= keep_fraction - 1e-9
    out = np.zeros_like(num)
    out[keep] = num[keep] / frac[keep]
    return out.astype(np.float32), keep


def _remove_dipole(m, mask):
    """Subtract the monopole and dipole fitted by least squares to unmasked pixels.

    (healpy's remove_dipole recognises masked pixels by comparing with its UNSEEN
    sentinel; after a float32 round trip the comparison fails and the "dipole"
    gets fitted to values of 1e30. A three-line fit is safer.)
    """
    hp = _hp()
    nside = hp.get_nside(m)
    idx = np.flatnonzero(mask)
    vec = np.asarray(hp.pix2vec(nside, idx)).T
    design = np.column_stack([np.ones(len(idx)), vec])
    coef, *_ = np.linalg.lstsq(design, m[idx].astype(np.float64), rcond=None)
    out = np.zeros_like(m, dtype=np.float32)
    out[idx] = (m[idx] - design @ coef).astype(np.float32)
    return out


def prepare_maps(*, force: bool = False) -> dict:
    """Build (or load) the analysis maps: SMICA and WMAP at Nside 512 and 64."""
    if CMB_CACHE.exists() and not force:
        d = np.load(CMB_CACHE)
        return {k: d[k] for k in d.files}
    hp = _hp()
    t0 = time.monotonic()
    _log("reading the two SMICA half-mission maps (Nside 2048)")
    m = np.zeros(hp.nside2npix(2048), np.float32)
    for key in ("hm1", "hm2"):
        m += 0.5 * hp.read_map(product(key).path, field=0, dtype=np.float32)
    m *= np.float32(K_TO_UK)
    mask = (hp.read_map(product("mask_common").path, dtype=np.float32) > 0.5)
    mask &= (hp.read_map(product("mask_hm").path, dtype=np.float32) > 0.5)
    _log(f"full-resolution mask keeps {mask.mean():.1%} of the sky", t0)
    smica, mask512 = degrade(m, mask, NSIDE, 1.0)
    del m
    smica = _remove_dipole(smica, mask512)
    _log(f"Nside {NSIDE}: mask keeps {mask512.mean():.1%}; rms {smica[mask512].std():.1f} uK", t0)

    wmap = hp.read_map(product("wmap_ilc").path, dtype=np.float32) * np.float32(1e3)  # mK -> uK
    if hp.get_nside(wmap) != NSIDE:
        wmap = hp.ud_grade(wmap, NSIDE)
    wmap = _remove_dipole(wmap, mask512)
    _log(f"WMAP ILC rms {wmap[mask512].std():.1f} uK (1 deg beam)", t0)

    smica64, mask64 = degrade(smica, mask512, 64, LOWRES_KEEP)
    wmap64, _ = degrade(wmap, mask512, 64, LOWRES_KEEP)
    out = {"smica": smica, "wmap": wmap, "mask": mask512,
           "smica64": smica64, "wmap64": wmap64, "mask64": mask64}
    np.savez(CMB_CACHE, **out)
    return out


# --------------------------------------------------------------------------
# Theory spectrum and simulations
# --------------------------------------------------------------------------


def theory_cl(lmax: int = LMAX_SIM) -> np.ndarray:
    """Planck 2018 best-fit LCDM C_l^TT in uK^2, l = 0..lmax (lensed)."""
    a = np.loadtxt(product("theory").path)
    ell, dl = a[:, 0].astype(int), a[:, 1]
    cl = np.zeros(lmax + 1)
    sel = ell <= lmax
    cl[ell[sel]] = dl[sel] * 2 * np.pi / (ell[sel] * (ell[sel] + 1))
    return cl


def beam(lmax: int, nside: int, fwhm_arcmin: float = SMICA_FWHM_ARCMIN) -> np.ndarray:
    hp = _hp()
    b = hp.gauss_beam(np.radians(fwhm_arcmin / 60), lmax=lmax)
    return b * hp.pixwin(nside, lmax=lmax)[: lmax + 1]


_ELL_CACHE: dict[int, np.ndarray] = {}


def _alm_ell(lmax: int):
    if lmax not in _ELL_CACHE:
        _ELL_CACHE[lmax] = _hp().Alm.getlm(lmax)
    return _ELL_CACHE[lmax]


def sim_alm(i: int, cl: np.ndarray | None = None, lmax: int = LMAX_SIM) -> np.ndarray:
    """Harmonic coefficients of simulated sky number ``i`` (no beam, uK).

    Drawn from our own seeded generator rather than healpy's synalm, which
    uses numpy's global random state: sky i must be reproducible anywhere.
    """
    if cl is None:
        cl = theory_cl(lmax)
    ell, m = _alm_ell(lmax)
    rng = np.random.default_rng(RANDOM_SEED + i)
    g = rng.standard_normal((2, len(ell)))
    amp = np.sqrt(cl[ell])
    alm = amp * (g[0] + 1j * g[1]) / np.sqrt(2.0)
    zero = m == 0
    alm[zero] = amp[zero] * g[0, zero]
    return alm


def alm_to_map(alm, nside: int, lmax: int = LMAX_SIM, fl=None) -> np.ndarray:
    """Beam- and pixel-window-smoothed map at ``nside`` from full-lmax alms.

    Only multipoles up to 3 nside - 1 are used; the rest are cut, so maps at
    every resolution come from the same sky.
    """
    hp = _hp()
    lcut = min(lmax, 3 * nside - 1)
    a = hp.resize_alm(alm, lmax, lmax, lcut, lcut) if lcut < lmax else alm
    a = hp.almxfl(a, beam(lcut, nside) if fl is None else fl[: lcut + 1])
    return hp.alm2map(a, nside, lmax=lcut).astype(np.float32)


def sim_map(i: int, nside: int = NSIDE, cl=None) -> np.ndarray:
    return alm_to_map(sim_alm(i, cl), nside)


def gal_to_vec(l_deg, b_deg):
    hp = _hp()
    return np.asarray(hp.ang2vec(np.asarray(l_deg), np.asarray(b_deg), lonlat=True))


def equatorial_to_galactic(ra, dec):
    """(l, b) in degrees using healpy's rotator (no astropy needed in hot loops)."""
    hp = _hp()
    r = hp.Rotator(coord=["C", "G"])
    th, ph = r(np.radians(90 - np.asarray(dec)), np.radians(np.asarray(ra)))
    return np.degrees(ph) % 360, 90 - np.degrees(th)


T_CMB_UK = COSMO["TCMB"] * 1e6
