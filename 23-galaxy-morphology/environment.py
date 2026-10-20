"""Environment: how crowded is the neighbourhood of each galaxy?

The measure is project D2's: Sigma_5, the surface density of galaxies out to
the 5th-nearest neighbour. Neighbours are counted within +/-1000 km/s along
the line of sight -- a cylinder rather than a sphere, because inside a cluster
the galaxies' own orbital motions smear the cluster out along the line of
sight (the "fingers of God"), and a sphere would lose its members.

Only *tracer* galaxies count as neighbours: D2's volume-limited sample, every
galaxy brighter than a fixed luminosity out to z = 0.106. With a fixed set of
tracers, Sigma_5 means the same thing at every distance. Counting every
catalogued galaxy instead would make nearby regions look crowded merely
because the survey sees fainter galaxies there.

Typical values: ~0.05 galaxies per (h^-1 Mpc)^2 in voids, ~0.4 at the median,
several per (h^-1 Mpc)^2 in the cores of rich clusters.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from astropy import units as u
from astropy_healpix import HEALPix
from scipy.spatial import cKDTree

import physics as ph
from config import (
    DENSITY_DV_KMS, DENSITY_NTH, EVOLUTION_Q, TRACER_M_LIT, TRACER_R_LIMIT, TRACER_Z_RANGE,
)

#: D2's footprint mask: HEALPix pixels of 0.84 deg^2, "surveyed" if they hold
#: at least 30% of the median pixel's galaxy count.
MASK_NSIDE = 64
MASK_MIN_FRACTION = 0.3


@dataclass
class Mask:
    """Which HEALPix pixels the survey covered (as in D2's ``clustering.Mask``)."""

    nside: int
    counts: np.ndarray
    surveyed: np.ndarray
    interior: np.ndarray

    @property
    def hp(self) -> HEALPix:
        return HEALPix(nside=self.nside, order="nested")

    def pixel(self, ra, dec) -> np.ndarray:
        return self.hp.lonlat_to_healpix(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg)

    def contains(self, ra, dec, which: str = "interior") -> np.ndarray:
        return getattr(self, which)[self.pixel(ra, dec)]

    def edge_distance_deg(self, ra, dec) -> np.ndarray:
        """Angular distance to the nearest *unsurveyed* pixel's edge, degrees.

        Measured to the nearest unsurveyed pixel centre, less the pixel's
        circumscribed radius, so it errs short: a galaxy is never called safely
        inside when it is not.
        """
        hp = self.hp
        lon, lat = hp.healpix_to_lonlat(np.flatnonzero(~self.surveyed))
        tree = cKDTree(_unit(lon.to_value(u.deg), lat.to_value(u.deg)))
        chord, _ = tree.query(_unit(ra, dec))
        sep = np.degrees(2.0 * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0)))
        # Pixel "radius": half the diagonal of a pixel of this area, generously.
        radius = np.sqrt(hp.pixel_area.to_value(u.deg**2)) * 0.75
        return sep - radius


def _unit(ra, dec) -> np.ndarray:
    ra, dec = np.radians(np.asarray(ra, float)), np.radians(np.asarray(dec, float))
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def build_mask(ra, dec, nside: int = MASK_NSIDE, min_fraction: float = MASK_MIN_FRACTION) -> Mask:
    """Derive the survey footprint from where the galaxies are (D2's method)."""
    hp = HEALPix(nside=nside, order="nested")
    pix = hp.lonlat_to_healpix(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg)
    counts = np.bincount(pix, minlength=hp.npix)
    median = np.median(counts[counts > 0])
    surveyed = counts >= min_fraction * median
    neighbours = hp.neighbours(np.arange(hp.npix))  # (8, npix); -1 where none
    ok = np.where(neighbours >= 0, surveyed[np.clip(neighbours, 0, None)], True)
    return Mask(nside=nside, counts=counts, surveyed=surveyed,
                interior=surveyed & ok.all(axis=0))


# ---------------------------------------------------------------------------
# Tracers
# ---------------------------------------------------------------------------

def derive_mgs(mgs: dict) -> dict:
    """D2's derived quantities for the Main Galaxy Sample."""
    z = mgs["z"].astype(float)
    r = (mgs["petroMag_r"] - mgs["extinction_r"]).astype(float)
    gr = (mgs["dered_g"] - mgs["dered_r"]).astype(float)
    kr = ph.k_correction("r", z, gr, "gr")
    out = {
        "specObjID": mgs["specObjID"], "ra": mgs["ra"].astype(float),
        "dec": mgs["dec"].astype(float), "z": z, "r": r, "gr_obs": gr, "kr": kr,
        "gr0": ph.rest_frame_colour("g", "r", gr, z),
        "M": ph.absolute_magnitude(r, z, kr),
    }
    out["band_offset"] = ph.bandshift_offset(kr, z)
    out["M_lit"] = ph.to_literature_magnitude(out["M"], z, out["band_offset"], EVOLUTION_Q)
    out["d_h"] = ph.comoving_distance_h(z)
    out["xyz_h"] = ph.cartesian(out["ra"], out["dec"], out["d_h"])
    return out


def select_tracers(m: dict) -> np.ndarray:
    """D2's volume-limited sample: M_r - 5 log h < -20, r < 17.6, 0.02 <= z < 0.106."""
    return ((m["r"] < TRACER_R_LIMIT) & (m["z"] >= TRACER_Z_RANGE[0])
            & (m["z"] < TRACER_Z_RANGE[1]) & (m["M_lit"] < TRACER_M_LIT))


# ---------------------------------------------------------------------------
# Sigma_5
# ---------------------------------------------------------------------------

def _los_split(d1sq, d2sq, s2):
    """Pair separation across (r_p) and along (pi) the line of sight.

    The line of sight is the direction to the pair's midpoint; with |p1|^2,
    |p2|^2 and |p1 - p2|^2 known, pi = | |p2|^2 - |p1|^2 | / |p1 + p2| and
    |p1 + p2|^2 = 2 (|p1|^2 + |p2|^2) - |p1 - p2|^2.
    """
    pi = np.abs(d2sq - d1sq) / np.sqrt(np.clip(2.0 * (d1sq + d2sq) - s2, 1e-12, None))
    rp = np.sqrt(np.clip(s2 - pi * pi, 0.0, None))
    return rp, pi


def local_density(target_xyz, target_z, target_id, tracer_xyz, tracer_id, *,
                  n_th: int = DENSITY_NTH, dv_kms: float = DENSITY_DV_KMS, k: int = 128,
                  chunk: int = 20000) -> dict:
    """Sigma_N around each target, counting tracers only (D2's estimator).

    A target that is itself a tracer is not its own neighbour (matched by the
    SDSS spectrum ID). Candidates are the ``k`` nearest tracers in 3D; if fewer
    than N of them fall inside the cylinder the galaxy sits in a very sparse
    region, d_N is only bounded from below, and it is flagged.

    Returns d_N (h^-1 Mpc, projected), Sigma_N (h^2 Mpc^-2) and the flag.
    """
    tree = cKDTree(tracer_xyz)
    t_d2 = (tracer_xyz**2).sum(1)
    window = ph.los_window(dv_kms, target_z)
    n = len(target_xyz)
    out_d = np.empty(n)
    lower = np.zeros(n, dtype=bool)
    for lo in range(0, n, chunk):
        hi = min(lo + chunk, n)
        dist, idx = tree.query(target_xyz[lo:hi], k=k + 1)
        own = tracer_id[idx] == np.asarray(target_id[lo:hi])[:, None]
        d1sq = (target_xyz[lo:hi] ** 2).sum(1)
        rp, pi = _los_split(d1sq[:, None], t_d2[idx], dist**2)
        rp = np.where((pi < window[lo:hi, None]) & ~own, rp, np.inf)
        rp.sort(axis=1)
        dn = rp[:, n_th - 1]
        bad = ~np.isfinite(dn)
        bound = np.sqrt(np.clip(dist[:, -1] ** 2 - window[lo:hi] ** 2, 0.0, None))
        dn[bad] = bound[bad]
        out_d[lo:hi] = dn
        lower[lo:hi] = bad
    return {"d_n": out_d, "sigma": n_th / (np.pi * out_d**2), "lower_bound": lower}


def usable_density(sigma_result: dict, z, ra, dec, mask: Mask) -> np.ndarray:
    """Galaxies whose Sigma_5 is not compromised by an edge of the survey.

    Three edges: the tracer volume's near and far redshift limits (the
    +/-1000 km/s cylinder must fit inside), and the footprint on the sky (the
    circle of radius d_5 must not reach unsurveyed sky, where missing
    neighbours would make a crowded region look sparse).
    """
    z = np.asarray(z, float)
    d = ph.comoving_distance_h(z)
    window = ph.los_window(DENSITY_DV_KMS, z)
    d_lo = ph.comoving_distance_h(TRACER_Z_RANGE[0])
    d_hi = ph.comoving_distance_h(TRACER_Z_RANGE[1])
    radial_ok = (d - window > d_lo) & (d + window < d_hi)
    d5_deg = np.degrees(sigma_result["d_n"] / d)
    sky_ok = mask.contains(ra, dec, "interior") & (mask.edge_distance_deg(ra, dec) > d5_deg)
    return radial_ok & sky_ok & ~sigma_result["lower_bound"]
