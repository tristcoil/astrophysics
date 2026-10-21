"""The CMB Cold Spot and the Eridanus supervoid.

In 2004 an unusually cold, ~10-degree patch was found in WMAP's southern sky at
Galactic (l, b) = (209, -57) deg (Vielva et al. 2004; Cruz et al. 2005). It
stood out when the map was filtered with a spherical Mexican-hat wavelet
(SMHW) of about 5 degrees -- a filter that responds to spots of that size
surrounded by a warmer ring. A decade later Szapudi et al. (2015) found a
supervoid in front of it, the Eridanus supervoid, about 200 h^-1 Mpc in
radius at z ~ 0.15-0.2. Could the void's ISW shadow *be* the Cold Spot?

Two questions, kept separate:

1. **How unusual is the Cold Spot?** We filter the real map and 1,000 simulated
   skies identically and ask how often chance makes a spot as cold. The answer
   is honest only up to a point: the spot was noticed first and the filter
   chosen afterwards (a posteriori), and with a Gaussian or top-hat filter it
   looks much less special (Zhang & Huterer 2010) -- so we show both.
2. **How much could the supervoid explain?** We compute the linear ISW profile
   of the published supervoid models with our own code and D3's cosmology,
   filter it exactly like the data, and compare.
"""

from __future__ import annotations

import time

import numpy as np

import cmb
import cosmology
from config import (COLDSPOT_LB, DERIVED_DIR, LITERATURE, N_SIMS, SMHW_SCALE_DEG,
                    SMHW_SCALES_DEG)

NSIDE_CS = 128
LMAX_CS = 3 * NSIDE_CS - 1
CS_SIMS = DERIVED_DIR / "sims_coldspot.npz"
#: Gaussian filter compared with the wavelet (FWHM, deg).
GAUSS_FWHM_DEG = 10.0
#: Positions are scored only if this far from the Galactic plane and if the
#: mask removes < 10% of the sky within 2R.
BMIN_DEG = 20.0


def _hp():
    import healpy as hp
    return hp


def smhw_profile(theta_rad, r_deg):
    """Spherical Mexican-hat wavelet (Martinez-Gonzalez et al. 2002), unnormalised shape."""
    R = np.radians(r_deg)
    y = 2 * np.tan(theta_rad / 2)
    norm = R * np.sqrt(1 + R**2 / 2 + R**4 / 4)
    return (1 + (y / 2) ** 2) ** 2 * (2 - (y / R) ** 2) * np.exp(-(y**2) / (2 * R**2)) / norm


def kernel_bl(kind: str, scale_deg: float, lmax: int = LMAX_CS):
    """Harmonic transfer function of an isotropic filter (SMHW of scale R, or Gaussian FWHM)."""
    hp = _hp()
    if kind == "gauss":
        return hp.gauss_beam(np.radians(scale_deg), lmax=lmax)
    th = np.linspace(0, np.radians(min(8 * scale_deg, 179.0)), 20000)
    return hp.beam2bl(smhw_profile(th, scale_deg), th, lmax)


class Analyzer:
    """Inpaint, filter and score maps at Nside 128, identically for data and sims."""

    def __init__(self):
        hp = _hp()
        d = cmb.prepare_maps()
        self.map128, mask = cmb.degrade(d["smica"], d["mask"], NSIDE_CS, 0.8)
        self.mask = mask
        npix = hp.nside2npix(NSIDE_CS)
        self.nbrs = hp.get_all_neighbours(NSIDE_CS, np.arange(npix))  # (8, npix)
        self.holes = np.flatnonzero(~mask)
        _, b = hp.pix2ang(NSIDE_CS, np.arange(npix), lonlat=True)
        self.filters = {("smhw", r): kernel_bl("smhw", r) for r in SMHW_SCALES_DEG}
        self.filters[("gauss", GAUSS_FWHM_DEG)] = kernel_bl("gauss", GAUSS_FWHM_DEG)
        # where a filter's footprint is mostly unmasked
        self.region = {}
        for key in self.filters:
            r = key[1] if key[0] == "smhw" else GAUSS_FWHM_DEG / 2.355
            th = np.linspace(0, np.radians(2 * r), 2000)
            tophat = hp.beam2bl(np.ones_like(th), th, LMAX_CS)
            frac = hp.alm2map(hp.almxfl(hp.map2alm(mask.astype(float), lmax=LMAX_CS),
                                        tophat / tophat[0]), NSIDE_CS, lmax=LMAX_CS)
            self.region[key] = (frac > 0.9) & (np.abs(b) >= BMIN_DEG)

    def inpaint(self, m, n_iter: int = 60):
        m = np.where(self.mask, m, 0.0).astype(np.float64)
        nb = self.nbrs[:, self.holes]
        valid = nb >= 0
        for _ in range(n_iter):
            vals = np.where(valid, m[np.maximum(nb, 0)], 0.0)
            m[self.holes] = vals.sum(0) / np.maximum(valid.sum(0), 1)
        return m

    def filtered(self, m):
        hp = _hp()
        alm = hp.map2alm(self.inpaint(m), lmax=LMAX_CS)
        return {k: hp.alm2map(hp.almxfl(alm, bl), NSIDE_CS, lmax=LMAX_CS)
                for k, bl in self.filters.items()}

    def score(self, m):
        """Per filter: minimum coefficient in units of the rms, and where it is."""
        hp = _hp()
        out = {}
        for k, f in self.filtered(m).items():
            reg = self.region[k]
            v = f[reg]
            sd = v.std()
            i = np.argmin(v)
            pix = np.flatnonzero(reg)[i]
            l, b = hp.pix2ang(NSIDE_CS, pix, lonlat=True)
            kurt = float(np.mean(((v - v.mean()) / sd) ** 4) - 3)
            out[k] = {"min_sigma": float(v[i] / sd), "l": float(l), "b": float(b), "rms": float(sd),
                      "kurtosis": kurt}
        return out


def coldspot_value(an: Analyzer, m, key):
    """Filtered value (in rms units) at the literature Cold Spot position."""
    hp = _hp()
    f = an.filtered(m)[key]
    pix = hp.ang2pix(NSIDE_CS, *COLDSPOT_LB, lonlat=True)
    return float(f[pix] / f[an.region[key]].std())


def run_sims(n_sims: int = N_SIMS, *, force: bool = False):
    if CS_SIMS.exists() and not force:
        return np.load(CS_SIMS, allow_pickle=True)["sims"].item()
    an = Analyzer()
    cl = cmb.theory_cl()
    acc = {k: {"min_sigma": [], "kurtosis": []} for k in an.filters}
    t0 = time.monotonic()
    for i in range(n_sims):
        m = cmb.alm_to_map(cmb.sim_alm(i, cl), NSIDE_CS)
        for k, v in an.score(m).items():
            acc[k]["min_sigma"].append(v["min_sigma"])
            acc[k]["kurtosis"].append(v["kurtosis"])
        if (i + 1) % 100 == 0:
            print(f"    cold-spot sim {i + 1}/{n_sims}  ({time.monotonic() - t0:.0f} s)", flush=True)
    out = {k: {kk: np.asarray(vv) for kk, vv in v.items()} for k, v in acc.items()}
    np.savez(CS_SIMS, sims=np.array(out, dtype=object))
    return out


# --------------------------------------------------------------------------
# Profiles and the supervoid
# --------------------------------------------------------------------------


def radial_profile(m, mask, nside, centre_lb, edges_deg):
    """Mean temperature in rings around a centre (masked pixels skipped)."""
    hp = _hp()
    vec = hp.ang2vec(*centre_lb, lonlat=True)
    out = []
    prev = np.array([], int)
    for r in edges_deg[1:]:
        disc = hp.query_disc(nside, vec, np.radians(r))
        ring = np.setdiff1d(disc, prev, assume_unique=True)
        ring = ring[mask[ring]]
        out.append(m[ring].mean() if len(ring) else np.nan)
        prev = disc
    return np.array(out)


SUPERVOID_MODELS = {
    # name: (depth > 0, r0 in h^-1 Mpc, z, depth quoted "today" or "at z"), all with the
    # compensated profile of Finelli et al. (2016). Finelli et al. quote the depth
    # extrapolated to today; DES and Szapudi et al. the underdensity at the void.
    "Finelli+16": (LITERATURE["finelli2016"]["delta0"], LITERATURE["finelli2016"]["r0_mpc_h"],
                   LITERATURE["finelli2016"]["z0"], "today"),
    "DES 2022": (-LITERATURE["kovacs2022"]["eridanus_delta0"],
                 LITERATURE["kovacs2022"]["eridanus_R_mpc_h"],
                 LITERATURE["kovacs2022"]["eridanus_z"], "at z"),
    "Szapudi+15": (-LITERATURE["szapudi2015"]["delta_m"][0], LITERATURE["szapudi2015"]["R_mpc_h"][0],
                   LITERATURE["szapudi2015"]["z"][0], "at z"),
}


def supervoid_profile(delta0, r0, z, theta_deg, quoted="at z", cosmo=cosmology.LCDM):
    """Linear ISW profile (uK) of a compensated supervoid of depth delta0 and scale r0 at z.

    The potential calculation wants the depth extrapolated to today
    (delta = D(z) delta_0); a depth quoted at the void's redshift is divided by D(z).
    """
    r = np.linspace(0, 8 * r0, 40001)
    d = cosmology.finelli_profile(r, delta0, r0)
    if quoted == "at z":
        d = d / cosmo.D(z)
    return cosmology.spherical_isw_profile(theta_deg, z, r, d, cosmo)


def profile_to_map(theta_deg, profile, centre_lb, nside=NSIDE_CS):
    hp = _hp()
    vec = np.asarray(hp.pix2vec(nside, np.arange(hp.nside2npix(nside))))
    c = hp.ang2vec(*centre_lb, lonlat=True)
    ang = np.degrees(np.arccos(np.clip(c @ vec, -1, 1)))
    return np.interp(ang, theta_deg, profile, right=0.0)
