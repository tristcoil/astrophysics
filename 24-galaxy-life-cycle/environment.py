"""Where a galaxy lives: survey footprint, fibre collisions and local density.

Three jobs:

1. **The footprint.** Which part of the sky the Main Galaxy Sample covers,
   mapped on a HEALPix grid from the *targets* -- every galaxy the survey
   selected, observed or not -- so the map does not inherit the holes that
   fibre collisions leave in the spectroscopic sample.

2. **Fibre collisions.** Two fibres cannot be placed closer than 55" on one
   plate, so some galaxies in close pairs and groups were never observed.
   Each observed galaxy gets a weight 1 / TSR, where the target sampling rate
   TSR is the fraction of targets within 55" that got a spectrum (Peng et al.
   2010). Groups and clusters, where collisions bite, are then not
   under-counted.

3. **Local density.** The projected density out to the 5th-nearest
   neighbour, counting only neighbours within +/-1000 km/s along the line of
   sight (as in D2 and Peng et al. 2010), expressed as an overdensity 1 + delta
   relative to the mean density of the tracer galaxies. Galaxies whose
   5th-neighbour circle reaches past the survey edge would look falsely
   isolated; they are flagged and left out.
"""

from __future__ import annotations

import healpy as hp
import numpy as np
from scipy.spatial import cKDTree

import cosmo_physics as cp
from config import (
    DENSITY_DV_KMS,
    DENSITY_NTH,
    EDGE_MIN_FRACTION,
    FIBRE_COLLISION_ARCSEC,
    MASK_NSIDE,
)

# ---------------------------------------------------------------------------
# Targets: deduplicate, observed or not
# ---------------------------------------------------------------------------


def unique_targets(t: dict) -> dict:
    """One row per target, flagged observed or not.

    In DR17, ``sdssTiledTargetAll`` links a target to the modern photometry
    (``photoObjID``) only through its spectrum: every target with a link was
    observed, and the targets that never got a fibre have photoObjID = 0 and
    no DR17 magnitude. No two entries lie within 1.5" of each other, so each
    unlinked entry is a distinct unobserved target. A few thousand linked
    entries share a photoObjID (pieces of one galaxy that the older target
    photometry split apart); those are merged.
    """
    pid = t["photoObjID"]
    linked = pid > 0
    order = np.argsort(pid[linked], kind="stable")
    pid_s = pid[linked][order]
    first = np.concatenate([[True], pid_s[1:] != pid_s[:-1]])
    grp = np.cumsum(first) - 1
    obs_linked = np.zeros(first.sum(), bool)
    np.logical_or.at(obs_linked, grp, t["specObjID"][linked][order] > 0)
    flag_linked = np.zeros(first.sum(), np.int16)
    np.maximum.at(flag_linked, grp, t["spec_flag"][linked][order])
    il = np.flatnonzero(linked)[order][first]
    iu = np.flatnonzero(~linked)
    idx = np.concatenate([il, iu])
    r = (t["petroMag_r"][idx] - t["extinction_r"][idx]).astype(float)
    r[len(il):] = np.nan  # unobserved targets: no DR17 photometry linked
    observed = np.concatenate([obs_linked, t["specObjID"][iu] > 0])
    spec_flag = np.concatenate([flag_linked, t["spec_flag"][iu]])
    return {"photoObjID": pid[idx], "ra": t["ra"][idx], "dec": t["dec"][idx], "r": r,
            "observed": observed, "spec_flag": spec_flag, "n_entries": int(len(pid)),
            "n_merged": int(linked.sum() - first.sum()),
            "n_unlinked": int((~linked).sum()),
            "n_unlinked_observed": int(np.sum(t["specObjID"][iu] > 0))}


# ---------------------------------------------------------------------------
# Footprint
# ---------------------------------------------------------------------------


def footprint(ra, dec, nside_coarse: int = 128, nside: int = MASK_NSIDE,
              min_fraction: float = 0.3) -> dict:
    """HEALPix footprint from target positions.

    At nside 128 (0.21 deg^2 pixels, ~20 targets each) a pixel is in if it
    holds >= 30% of the median occupied pixel's count -- generous enough to keep
    real voids, strict enough to drop the ragged edge. The mask is then
    upgraded to ``nside`` (0.052 deg^2) for edge-distance calculations.
    """
    npix = hp.nside2npix(nside_coarse)
    pix = hp.ang2pix(nside_coarse, ra, dec, lonlat=True)
    counts = np.bincount(pix, minlength=npix)
    med = np.median(counts[counts > 0])
    coarse = counts >= min_fraction * med
    fine = hp.ud_grade(coarse.astype(float), nside) > 0.5
    pix_area = hp.nside2pixarea(nside_coarse, degrees=True)
    return {"mask": fine, "nside": nside, "coarse": coarse, "nside_coarse": nside_coarse,
            "counts": counts, "area_deg2": float(coarse.sum() * pix_area),
            "median_count": float(med)}


def in_footprint(fp: dict, ra, dec) -> np.ndarray:
    return fp["mask"][hp.ang2pix(fp["nside"], ra, dec, lonlat=True)]


def _disc_pattern(n_rings: int = 5) -> np.ndarray:
    """Unit-disc sample points on concentric rings, equal area per point."""
    pts = [(0.0, 0.0)]
    for k in range(1, n_rings + 1):
        rad = k / n_rings
        n = 6 * k
        ang = (np.arange(n) + 0.5 * (k % 2)) * 2 * np.pi / n
        pts += list(zip(rad * np.cos(ang), rad * np.sin(ang)))
    return np.array(pts)


def disc_coverage(fp: dict, ra, dec, radius_deg) -> np.ndarray:
    """Fraction of a disc of given angular radius that lies in the footprint.

    Evaluated on 91 points spread uniformly over each disc (a fixed ring
    pattern; accurate to ~1%, ample for a 90% threshold).
    """
    pat = _disc_pattern()
    ra = np.asarray(ra, float)
    dec = np.asarray(dec, float)
    rad = np.asarray(radius_deg, float)
    cosd = np.cos(np.radians(dec))
    cov = np.zeros(len(ra))
    for dx, dy in pat:
        d2 = np.clip(dec + dy * rad, -90, 90)
        r2 = ra + dx * rad / np.clip(cosd, 1e-3, None)
        cov += fp["mask"][hp.ang2pix(fp["nside"], r2 % 360.0, d2, lonlat=True)]
    return cov / len(pat)


# ---------------------------------------------------------------------------
# Fibre collisions
# ---------------------------------------------------------------------------


def target_sampling_rate(tgt: dict, radius_arcsec: float = FIBRE_COLLISION_ARCSEC) -> np.ndarray:
    """For each target: fraction of targets within 55" (itself included) observed."""
    xyz = cp.unit_vectors(tgt["ra"], tgt["dec"])
    tree = cKDTree(xyz)
    chord = 2 * np.sin(np.radians(radius_arcsec / 3600.0) / 2)
    pairs = tree.query_pairs(chord, output_type="ndarray")
    n = len(xyz)
    obs = tgt["observed"].astype(float)
    n_all = np.ones(n)
    n_obs = obs.copy()
    for a, b in (pairs.T, pairs.T[::-1]):
        np.add.at(n_all, a, 1.0)
        np.add.at(n_obs, a, obs[b])
    return n_obs / n_all


# ---------------------------------------------------------------------------
# Local density
# ---------------------------------------------------------------------------


def _los_split(d1sq, d2sq, v2):
    """(r_p, pi) from squared radii and squared separation (as in D2)."""
    pi = np.abs(d2sq - d1sq) / np.sqrt(np.clip(2.0 * (d1sq + d2sq) - v2, 1e-12, None))
    rp = np.sqrt(np.clip(v2 - pi * pi, 0.0, None))
    return rp, pi


def nth_neighbour_density(xyz_gal, z_gal, xyz_trc, w_trc, *, n_th: int = DENSITY_NTH,
                          dv_kms: float = DENSITY_DV_KMS, k: int = 200,
                          self_match_tol: float = 1e-6, chunk: int = 20000) -> dict:
    """Weighted Sigma_N around each galaxy, counting tracer galaxies.

    For each galaxy: tracers within +/-dv along the line of sight, sorted by
    projected distance; d_N is where their cumulative weight (1 / TSR) first
    reaches N. The galaxy itself, if it is a tracer, is not counted.
    Sigma_N = N / (pi d_N^2), in Mpc^-2 (comoving).

    Candidates are the k nearest tracers in 3D; if their weight inside the
    cylinder never reaches N, d_N is only a lower bound and is flagged.
    """
    tree = cKDTree(xyz_trc)
    d2t = (xyz_trc**2).sum(1)
    window = cp.los_window_mpc(dv_kms, z_gal)
    n = len(xyz_gal)
    d_n = np.empty(n)
    bound = np.zeros(n, bool)
    for lo in range(0, n, chunk):
        hi = min(lo + chunk, n)
        dist, idx = tree.query(xyz_gal[lo:hi], k=k)
        d2g = (xyz_gal[lo:hi] ** 2).sum(1)
        rp, pi = _los_split(d2g[:, None], d2t[idx], dist**2)
        inside = (pi < window[lo:hi, None]) & (dist > self_match_tol)
        rp = np.where(inside, rp, np.inf)
        order = np.argsort(rp, axis=1)
        rp_s = np.take_along_axis(rp, order, 1)
        w_s = np.where(np.isfinite(rp_s), w_trc[np.take_along_axis(idx, order, 1)], 0.0)
        cw = np.cumsum(w_s, axis=1)
        reach = cw >= n_th
        has = reach.any(1)
        j = np.argmax(reach, axis=1)
        dn = rp_s[np.arange(hi - lo), j]
        # lower bound for galaxies whose k candidates never reach N
        lb = np.sqrt(np.clip(dist[:, -1] ** 2 - window[lo:hi] ** 2, 0.0, None))
        d_n[lo:hi] = np.where(has, dn, lb)
        bound[lo:hi] = ~has
    return {"d_n": d_n, "sigma": n_th / (np.pi * d_n**2), "lower_bound": bound}


def mean_tracer_surface_density(n_tracer_weighted: float, volume_mpc3: float, z) -> np.ndarray:
    """Expected Sigma for a uniform field: n_bar times the cylinder length 2*window."""
    nbar = n_tracer_weighted / volume_mpc3
    return nbar * 2.0 * cp.los_window_mpc(DENSITY_DV_KMS, z)


def edge_ok(fp: dict, ra, dec, d_n_mpc, z, min_fraction: float = EDGE_MIN_FRACTION):
    """True where >= 90% of the 5th-neighbour circle lies inside the footprint."""
    radius_deg = np.degrees(np.asarray(d_n_mpc) / cp.comoving_distance(z))
    cov = disc_coverage(fp, ra, dec, radius_deg)
    return cov >= min_fraction, cov
