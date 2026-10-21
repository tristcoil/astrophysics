"""The analysis: from a list of redshifts to a measured cosmic web.

Each step is a function below, run in this order by :func:`analyse`:

1. :func:`derive` -- distances, K-corrections, absolute magnitudes, colours.
2. Samples -- the flux-limited slice (for the map) and the volume-limited
   sample (for anything that is measured).
3. :func:`color_divider` -- where the red sequence ends and the blue cloud
   begins, found from the colour distribution itself.
4. :func:`local_density` -- how crowded each galaxy's neighbourhood is.
5. :func:`friends_of_friends` -- clusters and groups, and the fingers of God.
6. :func:`homogeneous_nz` -- the redshift distribution a perfectly uniform
   universe would give this survey, which is what shows the thinning of the
   map to be a selection effect.
7. :func:`clustering_measurement` -- the two-point correlation function,
   graded against Zehavi et al. (2011).
8. :func:`great_wall` and :func:`largest_void` -- the two structures the
   slice is famous for, measured.

Figures are drawn from what this module returns, never by re-running it.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

import clustering as cl
import cosmo_physics as cp
import literature as lit
from config import (
    ABELL_MATCH_RADIUS_H,
    ABELL_MATCH_DZ,
    DENSITY_DV_KMS,
    DENSITY_NTH,
    EVOLUTION_Q,
    FETCH_Z_MAX,
    FETCH_Z_MIN,
    FIT_RP_RANGE,
    FOF_B_PAR,
    FOF_B_PERP,
    MASK_MIN_FRACTION,
    MASK_NSIDE,
    MGS_R_LIMIT,
    N_JACKKNIFE,
    PAIRCOUNT_DIR,
    PI_MAX,
    RANDOM_FACTOR,
    RESULTS_PATH,
    RP_LOG_EDGES,
    SEED,
    SLICE_DEC_HALF_WIDTH,
    SLICE_RA_RANGE,
    VL_MAG_LIMIT_H,
    VL_R_LIMIT,
    VL_Z_MAX,
    VL_Z_MIN,
    WEDGE_Z_MAX,
    XI2D_BIN,
    XI2D_MAX,
)

# ---------------------------------------------------------------------------
# 1. Derived quantities
# ---------------------------------------------------------------------------


def derive(arrays: dict) -> dict:
    """Everything that follows from redshift, magnitudes and position.

    Returns a dict of per-galaxy arrays. Magnitudes are corrected for Milky Way
    dust throughout; ``M`` is the rest-frame r-band absolute magnitude for
    Planck 2018's H0, and ``M_lit`` the same luminosity in the convention of
    Zehavi et al. (2011) -- see :func:`cosmo_physics.to_literature_magnitude`.
    """
    z = arrays["z"].astype(float)
    r = (arrays["petroMag_r"] - arrays["extinction_r"]).astype(float)
    gr = (arrays["dered_g"] - arrays["dered_r"]).astype(float)
    kr = cp.k_correction("r", z, gr)
    kg = cp.k_correction("g", z, gr)
    d = cp.comoving_distance(z)
    g = {
        "ra": arrays["ra"].astype(float),
        "dec": arrays["dec"].astype(float),
        "z": z,
        "zerr_kms": arrays["zErr"].astype(float) * cp.C_KMS,
        "r": r,
        "gr_obs": gr,
        "gr_err": np.hypot(arrays["modelMagErr_g"], arrays["modelMagErr_r"]).astype(float),
        "kr": kr,
        "gr0": gr - kg + kr,
        "M": cp.absolute_magnitude(r, z, kr),
        "d_mpc": d,
        "d_h": d * cp.H_LITTLE,
    }
    g["xyz_h"] = cp.cartesian(g["ra"], g["dec"], g["d_h"])
    g["band_offset"] = band_shift_offset(g)
    g["M_lit"] = cp.to_literature_magnitude(g["M"], z, g["band_offset"], EVOLUTION_Q)
    return g


def band_shift_offset(g: dict) -> float:
    """Mean difference between magnitudes in the z = 0.1-shifted r band and ours.

    For a galaxy observed at exactly z = 0.1 the observed r filter *is* the
    shifted band, whose K-correction is -2.5 log10(1.1) for any spectrum. So the
    offset between the two conventions is our K_r(0.1) + 2.5 log10(1.1),
    evaluated over the galaxies that sit at z ~ 0.1.
    """
    near = (g["z"] > 0.095) & (g["z"] < 0.105)
    return float(np.median(g["kr"][near]) + 2.5 * np.log10(1.1))


# ---------------------------------------------------------------------------
# 2. Samples
# ---------------------------------------------------------------------------


def select_slice(g: dict, z_max: float = WEDGE_Z_MAX) -> np.ndarray:
    """The equatorial wedge: a thin slice of sky, every galaxy in it."""
    return ((np.abs(g["dec"]) < SLICE_DEC_HALF_WIDTH)
            & (g["ra"] >= SLICE_RA_RANGE[0]) & (g["ra"] < SLICE_RA_RANGE[1])
            & (g["z"] < z_max))


def select_volume_limited(g: dict) -> np.ndarray:
    """Every galaxy more luminous than the threshold, out to where the
    faintest of them still makes the flux limit."""
    return ((g["r"] < VL_R_LIMIT) & (g["z"] >= VL_Z_MIN) & (g["z"] < VL_Z_MAX)
            & (g["M_lit"] < VL_MAG_LIMIT_H))


def completeness_check(g: dict) -> dict:
    """Is the volume-limited sample really complete out to its far edge?

    Take a galaxy exactly at the luminosity threshold, move it to the sample's
    outer redshift, and ask how bright it would look. If it is still brighter
    than the flux limit, every galaxy above the threshold is visible throughout
    the volume. Redder galaxies have larger K-corrections, so they are the
    stringent case; we test across the observed colour range.
    """
    offset = g["band_offset"]
    # Invert to_literature_magnitude at z_max for M_lit = threshold.
    m_thresh = (VL_MAG_LIMIT_H - offset - EVOLUTION_Q * (VL_Z_MAX - 0.1)
                + 5.0 * np.log10(cp.H_LITTLE))
    near = (g["z"] > VL_Z_MAX - 0.005) & (g["z"] < VL_Z_MAX + 0.005)
    colours = np.percentile(g["gr_obs"][near], [5, 50, 95, 99.5])
    r_at_zmax = [float(m_thresh + cp.distance_modulus(VL_Z_MAX)
                       + cp.k_correction("r", VL_Z_MAX, c)) for c in colours]
    return {
        "threshold_M_planck18": float(m_thresh),
        "observed_gr_percentiles_5_50_95_99.5": [float(c) for c in colours],
        "r_at_zmax": r_at_zmax,
        "r_limit": VL_R_LIMIT,
        "complete": bool(max(r_at_zmax) < VL_R_LIMIT),
        "margin_mag": float(VL_R_LIMIT - max(r_at_zmax)),
    }


# ---------------------------------------------------------------------------
# 3. Colours: the red sequence and the blue cloud
# ---------------------------------------------------------------------------


def color_divider(M, gr0, sel, bin_width: float = 0.5) -> dict:
    """Find the colour that best separates red and blue galaxies, vs luminosity.

    Galaxy colours are bimodal: a tight *red sequence* of galaxies that have
    stopped forming stars (their light is old, red stars) and a broad *blue
    cloud* of galaxies still forming them (young, hot, blue stars dominate).
    In each 0.5-magnitude bin we fit two Gaussians to the rest-frame g - r
    distribution and take the colour where a galaxy is equally likely to
    belong to either. A straight line through those points is the divider.

    Bins where the two Gaussians do not produce a real dip in the combined
    distribution -- which happens at the bright end, where almost nothing is
    blue and the fit splits the red sequence instead -- are not used.
    """
    from sklearn.mixture import GaussianMixture

    M, gr0 = np.asarray(M)[sel], np.asarray(gr0)[sel]
    ok = (gr0 > 0.0) & (gr0 < 1.2)
    M, gr0 = M[ok], gr0[ok]
    edges = np.arange(np.floor(np.percentile(M, 0.5) * 2) / 2,
                      np.ceil(np.percentile(M, 99.5) * 2) / 2 + 1e-9, bin_width)
    bins = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (M >= lo) & (M < hi)
        if m.sum() < 2000:
            continue
        gm = GaussianMixture(2, random_state=SEED).fit(gr0[m].reshape(-1, 1))
        mu = gm.means_.ravel()
        blue, red = np.argsort(mu)
        grid = np.linspace(mu[blue], mu[red], 2001).reshape(-1, 1)
        post_red = gm.predict_proba(grid)[:, red]
        valley = float(grid[np.argmin(np.abs(post_red - 0.5)), 0])
        dens = np.exp(gm.score_samples(np.array([[mu[blue]], [valley], [mu[red]]])))
        sd = np.sqrt(gm.covariances_.ravel())
        bins.append({
            "M_center": float((lo + hi) / 2), "n": int(m.sum()),
            "blue_mean": float(mu[blue]), "red_mean": float(mu[red]),
            "blue_sd": float(sd[blue]), "red_sd": float(sd[red]),
            "red_weight": float(gm.weights_[red]), "divider": valley,
            "bimodal": bool(dens[1] < 0.9 * min(dens[0], dens[2])),
        })
    use = [b for b in bins if b["bimodal"]]
    slope, intercept = np.polyfit([b["M_center"] for b in use], [b["divider"] for b in use], 1)
    return {"slope": float(slope), "intercept": float(intercept), "bins": bins,
            "n_bins_used": len(use)}


def is_red(M, gr0, divider: dict) -> np.ndarray:
    return np.asarray(gr0) > divider["slope"] * np.asarray(M) + divider["intercept"]


# ---------------------------------------------------------------------------
# 4. Environment
# ---------------------------------------------------------------------------


def _los_split(d1sq, d2sq, v2):
    """(r_p, pi) from squared radii and squared separation; see clustering."""
    pi = np.abs(d2sq - d1sq) / np.sqrt(np.clip(2.0 * (d1sq + d2sq) - v2, 1e-12, None))
    rp = np.sqrt(np.clip(v2 - pi * pi, 0.0, None))
    return rp, pi


def local_density(xyz, z, n_th: int = DENSITY_NTH, dv_kms: float = DENSITY_DV_KMS,
                  k: int = 96, chunk: int = 20000) -> dict:
    """Sigma_N: surface density out to the N-th nearest neighbour.

    For each galaxy, consider the others within +/-``dv_kms`` along the line of
    sight -- a cylinder rather than a sphere, so that a cluster galaxy's
    neighbours are not lost to the fingers-of-God smearing -- and find the
    projected distance d_N to the N-th nearest. Then Sigma_N = N / (pi d_N^2).

    Candidates are the ``k`` nearest neighbours in 3D. If fewer than N of those
    fall inside the cylinder, the galaxy is in a very sparse region and d_N is
    only known to be larger than a bound; such values are flagged.
    """
    tree = cKDTree(xyz)
    d2 = (xyz**2).sum(1)
    window = cp.los_window(dv_kms, z)
    out_d = np.empty(len(xyz))
    lower = np.zeros(len(xyz), dtype=bool)
    for lo in range(0, len(xyz), chunk):
        hi = min(lo + chunk, len(xyz))
        dist, idx = tree.query(xyz[lo:hi], k=k + 1)
        dist, idx = dist[:, 1:], idx[:, 1:]  # drop the galaxy itself
        rp, pi = _los_split(d2[lo:hi, None], d2[idx], dist**2)
        rp = np.where(pi < window[lo:hi, None], rp, np.inf)
        rp.sort(axis=1)
        dn = rp[:, n_th - 1]
        bad = ~np.isfinite(dn)
        # Lower bound: the N-th cylinder neighbour lies beyond the k-th 3D one.
        bound = np.sqrt(np.clip(dist[:, -1] ** 2 - window[lo:hi] ** 2, 0.0, None))
        dn[bad] = bound[bad]
        out_d[lo:hi] = dn
        lower[lo:hi] = bad
    return {"d_n": out_d, "sigma": n_th / (np.pi * out_d**2), "lower_bound": lower}


def red_fraction_by_density(sigma, red, sel, n_bins: int = 8) -> dict:
    """Fraction of red galaxies in bins of local density (equal-count bins)."""
    s, r = np.asarray(sigma)[sel], np.asarray(red)[sel]
    edges = np.quantile(s, np.linspace(0, 1, n_bins + 1))
    idx = np.clip(np.searchsorted(edges, s, side="right") - 1, 0, n_bins - 1)
    frac, err, centre, n = [], [], [], []
    for b in range(n_bins):
        m = idx == b
        f = r[m].mean()
        frac.append(float(f))
        err.append(float(np.sqrt(f * (1 - f) / m.sum())))
        centre.append(float(np.median(s[m])))
        n.append(int(m.sum()))
    return {"sigma_median": centre, "red_fraction": frac, "red_fraction_err": err,
            "n": n, "edges": [float(e) for e in edges]}


# ---------------------------------------------------------------------------
# 5. Groups and clusters
# ---------------------------------------------------------------------------


def friends_of_friends(xyz, mean_separation: float, b_perp: float = FOF_B_PERP,
                       b_par: float = FOF_B_PAR) -> np.ndarray:
    """Group finder: link any two galaxies closer than the linking lengths.

    Friends of friends are joined too, so a group is a chain of links. The
    linking lengths are fractions of the mean inter-galaxy separation, and the
    one along the line of sight is ~5x longer, because a cluster's members are
    spread out in redshift by their orbital speeds (hundreds of km/s, i.e.
    several Mpc of apparent depth) while their true separations are ~1 Mpc.
    Returns a group label per galaxy.
    """
    l_perp, l_par = b_perp * mean_separation, b_par * mean_separation
    tree = cKDTree(xyz)
    pairs = tree.query_pairs(np.hypot(l_perp, l_par), output_type="ndarray")
    d2 = (xyz**2).sum(1)
    v2 = ((xyz[pairs[:, 0]] - xyz[pairs[:, 1]]) ** 2).sum(1)
    rp, pi = _los_split(d2[pairs[:, 0]], d2[pairs[:, 1]], v2)
    linked = pairs[(rp < l_perp) & (pi < l_par)]
    n = len(xyz)
    graph = coo_matrix((np.ones(len(linked)), (linked[:, 0], linked[:, 1])), shape=(n, n))
    _, labels = connected_components(graph, directed=False)
    return labels


def _mean_direction(ra, dec):
    v = cp.cartesian(ra, dec, np.ones(len(ra))).mean(axis=0)
    ra_c = np.degrees(np.arctan2(v[1], v[0])) % 360.0
    dec_c = np.degrees(np.arcsin(v[2] / np.linalg.norm(v)))
    return float(ra_c), float(dec_c)


def _angular_sep_deg(ra1, dec1, ra2, dec2):
    ra1, dec1, ra2, dec2 = map(np.radians, (ra1, dec1, ra2, dec2))
    c = (np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2))
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def cluster_dispersion(g: dict, sky_tree: cKDTree, ra: float, dec: float, z: float,
                       r_aperture_h: float = ABELL_MATCH_RADIUS_H, window_kms: float = 3000.0,
                       nsigma: float = 3.0) -> dict:
    """Velocity dispersion of a cluster, from every catalogued galaxy around it.

    The friends-of-friends membership that finds a cluster is not a good one
    for measuring its dispersion: chains of links break in the velocity tails,
    trimming exactly the fastest members and biasing the dispersion low by
    20-40%. So, as is standard (e.g. Colless & Dunn 1996), we take all galaxies
    within one Abell radius on the sky and +/-3000 km/s, and iteratively discard
    those more than 3 sigma from the median velocity. The standard deviation of
    a Gaussian truncated at +/-3 sigma is 0.9866 of the true one; we divide by
    that.
    """
    d = float(cp.comoving_distance_h(z))
    theta = r_aperture_h / d
    centre = cp.cartesian([ra], [dec], [1.0])[0]
    cand = np.asarray(sky_tree.query_ball_point(centre, 2.0 * np.sin(theta / 2.0)), dtype=int)
    v = cp.C_KMS * (g["z"][cand] - z) / (1.0 + z)
    m = np.abs(v) < window_kms
    for _ in range(50):
        centre_v = np.median(v[m])
        s = np.std(v[m], ddof=1)
        new = (np.abs(v - centre_v) < nsigma * s) & (np.abs(v) < window_kms)
        if (new == m).all():
            break
        m = new
    z_c = z + (1.0 + z) * float(np.median(v[m])) / cp.C_KMS
    sigma, err = cp.velocity_dispersion(g["z"][cand][m])
    return {"sigma_v": sigma / 0.9866, "sigma_v_err": err / 0.9866, "n": int(m.sum()),
            "z": float(z_c)}


def group_catalog(g: dict, members: np.ndarray, labels: np.ndarray,
                  abell: pd.DataFrame | None, min_members: int = 10) -> list[dict]:
    """Properties of every group with at least ``min_members`` galaxies.

    ``members`` indexes the galaxies the finder ran on; ``labels`` their groups.
    The "finger" is how much deeper a group looks along the line of sight than
    across it -- for a real, roughly round cluster, that ratio is caused
    entirely by the members' orbital motions. The velocity dispersion comes
    from :func:`cluster_dispersion`, not from the finder's members.
    """
    sky_tree = cKDTree(cp.cartesian(g["ra"], g["dec"], np.ones(len(g["ra"]))))
    counts = np.bincount(labels)
    rich = np.flatnonzero(counts >= min_members)
    out = []
    for lab in rich[np.argsort(-counts[rich])]:
        idx = members[labels == lab]
        ra_c, dec_c = _mean_direction(g["ra"][idx], g["dec"][idx])
        z_c = float(np.median(g["z"][idx]))
        disp = cluster_dispersion(g, sky_tree, ra_c, dec_c, z_c)
        z_c, sigma, sigma_err = disp["z"], disp["sigma_v"], disp["sigma_v_err"]
        # Transverse size: rms projected distance from the centre.
        d_c = float(cp.comoving_distance_h(z_c))
        theta = np.radians(_angular_sep_deg(g["ra"][idx], g["dec"][idx], ra_c, dec_c))
        r_perp = float(np.sqrt(np.mean((d_c * theta) ** 2)))
        depth = float(np.std(g["d_h"][idx], ddof=1))
        grp = {
            "n": int(len(idx)), "ra": ra_c, "dec": dec_c, "z": z_c,
            "distance_h": d_c, "distance_mpc": d_c / cp.H_LITTLE,
            "sigma_v": sigma, "sigma_v_err": sigma_err, "n_dispersion": disp["n"],
            "r_perp_rms_h": r_perp, "los_depth_rms_h": depth,
            "elongation": depth / r_perp if r_perp > 0 else float("nan"),
            "abell": None, "abell_sep_arcmin": None, "abell_z": None,
        }
        if abell is not None:
            grp.update(match_abell(ra_c, dec_c, z_c, abell))
        out.append(grp)
    return out


def match_abell(ra, dec, z, abell: pd.DataFrame) -> dict:
    """Name a group after the nearest Abell cluster on the sky, if any.

    A match needs the Abell centre within one Abell radius of the group's
    centre (``ABELL_MATCH_RADIUS_H``, projected at the group's distance) and,
    where the Abell cluster has a measured redshift, agreement within
    ``ABELL_MATCH_DZ``. Most Abell clusters have no redshift in the
    catalogue; those match on position alone, and the output records it.
    """
    sep = _angular_sep_deg(abell["ra"].to_numpy(), abell["dec"].to_numpy(), ra, dec)
    max_sep = np.degrees(ABELL_MATCH_RADIUS_H / float(cp.comoving_distance_h(z)))
    zab = abell["z"].to_numpy()
    ok = (sep < max_sep) & (~np.isfinite(zab) | (np.abs(zab - z) < ABELL_MATCH_DZ))
    if not ok.any():
        return {}
    best = np.flatnonzero(ok)[np.argmin(sep[ok])]
    return {"abell": f"Abell {int(abell['aco'].iloc[best])}",
            "abell_sep_arcmin": float(sep[best] * 60.0),
            "abell_z": None if not np.isfinite(zab[best]) else float(zab[best])}


# ---------------------------------------------------------------------------
# 6. Why the map thins out
# ---------------------------------------------------------------------------


def zmax_visible(M, gr_obs, r_limit: float, z_obs, z_ceiling: float = FETCH_Z_MAX):
    """Largest redshift at which each galaxy would still pass the flux limit.

    Solves r_limit = M + (m - M)(z) + K_r(z) by bisection, holding the observed
    colour fixed when evaluating the K-correction -- a small approximation over
    the range that matters.
    """
    lo = np.asarray(z_obs, dtype=float).copy()
    hi = np.full_like(lo, z_ceiling)
    M, gr_obs = np.asarray(M, dtype=float), np.asarray(gr_obs, dtype=float)
    visible_at_ceiling = (M + cp.distance_modulus(hi) + cp.k_correction("r", hi, gr_obs)) < r_limit
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        bright = (M + cp.distance_modulus(mid) + cp.k_correction("r", mid, gr_obs)) < r_limit
        lo = np.where(bright, mid, lo)
        hi = np.where(bright, hi, mid)
    return np.where(visible_at_ceiling, z_ceiling, 0.5 * (lo + hi))


def homogeneous_nz(z, zmax, z_edges, z_min: float = FETCH_Z_MIN):
    """Redshift histogram a *uniform* universe would give this survey.

    Each observed galaxy stands for a population with space density 1 / V_max,
    where V_max is the volume in which it could have been seen (Schmidt 1968).
    Spreading that density evenly through its own V_max and summing over all
    galaxies gives the expected counts per redshift bin if galaxies were spread
    uniformly -- everything the flux limit does, and nothing the cosmic web
    does. By construction it has the same total as the data.
    """
    dz = cp.comoving_distance(np.asarray(z_edges))
    d_min = float(cp.comoving_distance(z_min))
    dmax = cp.comoving_distance(np.asarray(zmax))
    vmax = dmax**3 - d_min**3  # the solid angle cancels
    expected = np.empty(len(z_edges) - 1)
    for k in range(len(z_edges) - 1):
        lo = np.clip(dz[k] ** 3, d_min**3, None)
        hi = np.minimum(dz[k + 1] ** 3, dmax**3)
        expected[k] = np.sum(np.clip(hi - lo, 0.0, None) / vmax)
    return expected


# ---------------------------------------------------------------------------
# 7. Clustering
# ---------------------------------------------------------------------------


def _binning() -> cl.Binning:
    lo, hi, step = RP_LOG_EDGES
    return cl.Binning(rp_edges=10 ** np.arange(lo, hi + step / 2, step), pi_max=PI_MAX,
                      pi_bin=2.0, map_bin=XI2D_BIN, map_max=XI2D_MAX)


def _cache_key(xyz, n_random: int) -> str:
    h = hashlib.sha1()
    h.update(np.ascontiguousarray(xyz).tobytes())
    h.update(json.dumps([n_random, SEED, RP_LOG_EDGES, PI_MAX, XI2D_BIN, XI2D_MAX,
                         MASK_NSIDE, MASK_MIN_FRACTION, N_JACKKNIFE]).encode())
    return h.hexdigest()[:16]


def clustering_measurement(g: dict, sample: np.ndarray, mask: cl.Mask, *,
                           workers: int | None = None, quiet: bool = False) -> dict:
    """w_p(r_p) and xi(r_p, pi) of the volume-limited sample.

    Randoms are scattered uniformly over the same sky (the interior mask) and
    uniformly in volume over the same distance range, so they share the
    sample's selection exactly and differ from it only in having no structure.
    Pair counts are cached: they are the only slow step.
    """
    rng = np.random.default_rng(SEED)
    xyz = g["xyz_h"][sample]
    n_d = len(xyz)
    n_r = RANDOM_FACTOR * n_d
    ra_r, dec_r = cl.random_sky(n_r, mask, rng)
    d_min, d_max = float(cp.comoving_distance_h(VL_Z_MIN)), float(cp.comoving_distance_h(VL_Z_MAX))
    xyz_r = cp.cartesian(ra_r, dec_r, cl.random_volume_distances(n_r, d_min, d_max, rng))

    region_of_pixel = cl.jackknife_regions(mask, N_JACKKNIFE, SEED)
    reg_d = region_of_pixel[mask.pixel(g["ra"][sample], g["dec"][sample])]
    reg_r = region_of_pixel[mask.pixel(ra_r, dec_r)]
    assert (reg_d >= 0).all() and (reg_r >= 0).all(), "sample must lie inside the mask"

    b = _binning()
    path = PAIRCOUNT_DIR / f"pairs_{_cache_key(xyz, n_r)}.npz"
    fields = ("wp", "wp_region", "map2d", "rp_sum")
    if path.exists():
        with np.load(path) as f:
            counts = {k: cl.PairCounts(**{fld: f[f"{k}_{fld}"] for fld in fields},
                                       n_pairs_examined=int(f[f"{k}_n"]))
                      for k in ("dd", "dr", "rr")}
        if not quiet:
            print("    pair counts: cached")
    else:
        import time

        counts = {}
        for key, args in (("dd", (xyz, reg_d, b, N_JACKKNIFE)),
                          ("dr", (xyz, reg_d, b, N_JACKKNIFE, xyz_r, reg_r)),
                          ("rr", (xyz_r, reg_r, b, N_JACKKNIFE))):
            t0 = time.time()
            counts[key] = cl.count_pairs(*args, workers=workers,
                                         progress=None if quiet else key.upper())
            if not quiet:
                print(f"    {key.upper()}: {counts[key].n_pairs_examined:.3g} pairs examined "
                      f"in {time.time() - t0:.0f}s")
        np.savez_compressed(path, **{f"{k}_{fld}": getattr(c, fld)
                                     for k, c in counts.items() for fld in fields},
                            **{f"{k}_n": c.n_pairs_examined for k, c in counts.items()})

    res = cl.measure(counts["dd"], counts["dr"], counts["rr"],
                     np.bincount(reg_d, minlength=N_JACKKNIFE).astype(float),
                     np.bincount(reg_r, minlength=N_JACKKNIFE).astype(float), b)

    fits = {
        "full_cov": cl.fit_power_law(res.rp, res.wp, res.wp_cov, FIT_RP_RANGE, N_JACKKNIFE),
        "diagonal": cl.fit_power_law(res.rp, res.wp, res.wp_cov, FIT_RP_RANGE, N_JACKKNIFE,
                                     diagonal=True),
        "diagonal_zehavi_range": cl.fit_power_law(res.rp, res.wp, res.wp_cov,
                                                  lit.ZEHAVI_MR20["fit_range"], N_JACKKNIFE,
                                                  diagonal=True),
    }
    # Point-by-point comparison with Zehavi et al.'s published w_p. Our bins
    # are theirs, so the pair-weighted separations nearly coincide; their
    # points beyond our largest bin are not compared (NaN).
    inside = lit.ZEHAVI_WP_RP <= b.rp_edges[-1]
    ours_at_theirs = np.interp(np.log(lit.ZEHAVI_WP_RP), np.log(res.rp), res.wp)
    ratio = np.where(inside, ours_at_theirs / lit.ZEHAVI_WP_MR20, np.nan)
    return {
        "result": res, "binning": b, "fits": fits,
        "ratio_to_zehavi": ratio, "xyz_random": xyz_r, "ra_random": ra_r, "dec_random": dec_r,
        "n_data": n_d, "n_random": n_r,
        "area_deg2": mask.area_deg2(), "volume_h3": float(cl_volume(mask, d_min, d_max)),
    }


def cl_volume(mask: cl.Mask, d_min: float, d_max: float) -> float:
    return float(cp.shell_volume(mask.solid_angle_sr(), d_min, d_max))


# ---------------------------------------------------------------------------
# 8. The Sloan Great Wall and the biggest void in the slice
# ---------------------------------------------------------------------------


def slice_plane(ra, d, ra_mid: float):
    """Project slice galaxies onto the plane of the equator, with the slice's
    central right ascension pointing "up" (the +y axis)."""
    ang = np.radians(np.asarray(ra) - ra_mid)
    return np.column_stack([-np.asarray(d) * np.sin(ang), np.asarray(d) * np.cos(ang)])


def density_map(xy, extent, cell: float, smooth: float):
    """Gaussian-smoothed surface density of points on a regular grid."""
    from scipy.ndimage import gaussian_filter

    (x0, x1), (y0, y1) = extent
    nx, ny = int(np.ceil((x1 - x0) / cell)), int(np.ceil((y1 - y0) / cell))
    h, xe, ye = np.histogram2d(xy[:, 0], xy[:, 1], bins=[nx, ny], range=[[x0, x1], [y0, y1]])
    return gaussian_filter(h, smooth / cell), xe, ye


def great_wall(g: dict, sel: np.ndarray, randoms_xy: np.ndarray, *, cell: float = 2.0,
               smooth: float = 5.0, thresholds=(1.0, 2.0)) -> dict:
    """Measure the biggest connected overdensity in the volume-limited slice.

    Galaxies and random points (same selection, no structure) are each turned
    into a smoothed density map; their ratio is the density relative to the
    mean, 1 + delta. Grid cells above a threshold are joined into connected
    structures, and the largest one reported -- with its end-to-end length,
    the largest distance between any two of its cells.

    The answer depends on the threshold (a lower one joins more structure into
    the wall), which is why two thresholds are reported, and on the smoothing.
    "The length of a wall" has no unique definition; Gott et al. measured theirs
    by eye along its curve.
    """
    from scipy.ndimage import label

    ra_mid = 0.5 * (SLICE_RA_RANGE[0] + SLICE_RA_RANGE[1])
    xy = slice_plane(g["ra"][sel], g["d_h"][sel], ra_mid)
    extent = ((-330.0, 330.0), (-10.0, 330.0))
    dg, xe, ye = density_map(xy, extent, cell, smooth)
    dr, _, _ = density_map(randoms_xy, extent, cell, smooth)
    norm = len(randoms_xy) / len(xy)
    valid = dr > 0.2 * np.median(dr[dr > 0])  # well inside the surveyed wedge
    with np.errstate(divide="ignore", invalid="ignore"):
        one_plus_delta = np.where(valid, dg * norm / dr, np.nan)

    xc = 0.5 * (xe[:-1] + xe[1:])
    yc = 0.5 * (ye[:-1] + ye[1:])
    X, Y = np.meshgrid(xc, yc, indexing="ij")
    out = {"smoothing_h": smooth, "cell_h": cell, "structures": []}
    outline = None
    for thr in thresholds:
        lab, n = label(np.nan_to_num(one_plus_delta) > 1.0 + thr)
        if n == 0:
            continue
        sizes = np.bincount(lab.ravel())[1:]
        big = int(np.argmax(sizes)) + 1
        if outline is None:
            outline = lab == big  # the lowest threshold's structure, for figures
        pts = np.column_stack([X[lab == big], Y[lab == big]])
        # End-to-end length: the largest distance between two cells of it.
        from scipy.spatial import ConvexHull

        hull = pts[ConvexHull(pts).vertices]
        dmat = np.sqrt(((hull[:, None, :] - hull[None, :, :]) ** 2).sum(-1))
        i, j = np.unravel_index(np.argmax(dmat), dmat.shape)
        ends = hull[[i, j]]
        ra_ends = (ra_mid - np.degrees(np.arctan2(ends[:, 0], ends[:, 1])))
        d_all = np.hypot(pts[:, 0], pts[:, 1])
        ra_all = ra_mid - np.degrees(np.arctan2(pts[:, 0], pts[:, 1]))
        out["structures"].append({
            "threshold_delta": thr,
            "area_h2": float(len(pts) * cell**2),
            "length_h": float(dmat[i, j]),
            "length_mpc": float(dmat[i, j] / cp.H_LITTLE),
            "ra_range": [float(ra_all.min()), float(ra_all.max())],
            "ra_ends": [float(r) for r in ra_ends],
            "distance_range_h": [float(d_all.min()), float(d_all.max())],
            "median_distance_h": float(np.median(d_all)),
            "median_z": float(cp.redshift_at_distance(np.median(d_all) / cp.H_LITTLE)),
            "mean_overdensity": float(np.nanmean(one_plus_delta[lab == big]) - 1.0),
            "ends_xy_h": ends.tolist(),
        })
    out["map"] = {"one_plus_delta": one_plus_delta, "x_edges": xe, "y_edges": ye,
                  "structure": outline}
    return out


def slice_coverage(mask: cl.Mask, ra_mid: float, d_lo: float, d_hi: float,
                   cell: float = 1.0):
    """Which cells of the slice plane lie inside surveyed space.

    A cell is covered if it is within the distance range and, at its right
    ascension, the survey mask covers at least 90% of the slice's declination
    band. Returns the grid and the boolean coverage map.
    """
    x = np.arange(-330.0, 330.0, cell) + cell / 2
    y = np.arange(0.0, 330.0, cell) + cell / 2
    X, Y = np.meshgrid(x, y, indexing="ij")
    D = np.hypot(X, Y)
    RA = ra_mid - np.degrees(np.arctan2(X, Y))
    ra_grid = np.arange(SLICE_RA_RANGE[0], SLICE_RA_RANGE[1] + 0.05, 0.1)
    dec_grid = np.linspace(-SLICE_DEC_HALF_WIDTH, SLICE_DEC_HALF_WIDTH, 41)
    RR, DD = np.meshgrid(ra_grid, dec_grid, indexing="ij")
    frac = mask.contains(RR.ravel(), DD.ravel(), "surveyed").reshape(RR.shape).mean(axis=1)
    cover_ra = np.interp(RA, ra_grid, frac, left=0.0, right=0.0)
    covered = (D > d_lo) & (D < d_hi) & (cover_ra >= 0.9)
    return X, Y, covered


def largest_void(xy, mask: cl.Mask, ra_mid: float, d_lo: float, d_hi: float,
                 *, cell: float = 1.0) -> dict:
    """The largest circle in the slice containing no galaxy (in projection).

    Every point of a fine grid gets its distance to the nearest galaxy; the
    point with the largest such distance is the centre of the largest empty
    circle. The circle must lie entirely inside surveyed space -- otherwise
    the edge of the survey, or a hole in it, would masquerade as a void.
    It is measured in the 2D slice, so a void has to be empty through the
    slice's full thickness to count; voids in 3D are larger.
    """
    from scipy.ndimage import distance_transform_edt

    X, Y, covered = slice_coverage(mask, ra_mid, d_lo, d_hi, cell)
    room = distance_transform_edt(covered) * cell  # distance to uncovered space
    grid = np.column_stack([X[covered], Y[covered]])
    d_gal, _ = cKDTree(xy).query(grid)
    d_gal = np.where(d_gal <= room[covered], d_gal, 0.0)
    k = int(np.argmax(d_gal))
    c = grid[k]
    return {"center_xy_h": c.tolist(), "radius_h": float(d_gal[k]),
            "radius_mpc": float(d_gal[k] / cp.H_LITTLE),
            "distance_h": float(np.hypot(*c)),
            "ra_center": float(ra_mid - np.degrees(np.arctan2(c[0], c[1])))}


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def analyse(arrays: dict, meta: dict, abell: pd.DataFrame | None, *,
            workers: int | None = None, quiet: bool = False):
    """Run every step. Returns (results for the record, artefacts for figures)."""
    say = (lambda *a: None) if quiet else print
    rng = np.random.default_rng(SEED + 1)

    say("  deriving distances, K-corrections, luminosities")
    g = derive(arrays)
    mask = cl.build_mask(g["ra"], g["dec"], MASK_NSIDE, MASK_MIN_FRACTION)

    sl = select_slice(g)
    vl = select_volume_limited(g)
    vl_in = vl & mask.contains(g["ra"], g["dec"], "interior")
    completeness = completeness_check(g)

    say("  colours: fitting the red sequence / blue cloud divider")
    divider = color_divider(g["M"], g["gr0"], g["z"] < 0.2)
    red = is_red(g["M"], g["gr0"], divider)

    say("  environment: local density around each galaxy")
    vl_idx = np.flatnonzero(vl)
    dens = local_density(g["xyz_h"][vl_idx], g["z"][vl_idx])
    sigma = np.full(len(g["z"]), np.nan)
    sigma[vl_idx] = dens["sigma"]
    # Keep galaxies away from the survey's edges, where missing neighbours
    # would make a crowded region look sparse.
    window = cp.los_window(DENSITY_DV_KMS, g["z"])
    d_lo = cp.comoving_distance_h(VL_Z_MIN)
    d_hi = cp.comoving_distance_h(VL_Z_MAX)
    env_sel = vl_in & (g["d_h"] - window > d_lo) & (g["d_h"] + window < d_hi)
    env_sel &= np.isfinite(sigma)
    lum_split = -21.0
    env = {
        "all": red_fraction_by_density(sigma, red, env_sel),
        "bright": red_fraction_by_density(sigma, red, env_sel & (g["M_lit"] < lum_split)),
        "faint": red_fraction_by_density(sigma, red, env_sel & (g["M_lit"] >= lum_split)),
        "luminosity_split_M_lit": lum_split,
        "n": int(env_sel.sum()),
        "n_lower_bound": int((dens["lower_bound"] & env_sel[vl_idx]).sum()),
    }

    say("  groups: friends-of-friends")
    n_mean = vl_in.sum() / cl_volume(mask, d_lo, d_hi)
    mean_sep = float(n_mean ** (-1.0 / 3.0))
    labels = friends_of_friends(g["xyz_h"][vl_idx], mean_sep)
    groups = group_catalog(g, vl_idx, labels, abell)
    in_slice = [grp for grp in groups
                if abs(grp["dec"]) < SLICE_DEC_HALF_WIDTH
                and SLICE_RA_RANGE[0] <= grp["ra"] < SLICE_RA_RANGE[1]]

    say("  selection: what a uniform universe would look like")
    z_edges = np.arange(0.0, 0.3001, 0.005)
    zmax_all = zmax_visible(g["M"], g["gr_obs"], MGS_R_LIMIT, g["z"])
    nz_expected = homogeneous_nz(g["z"], zmax_all, z_edges)
    nz_observed = np.histogram(g["z"], bins=z_edges)[0]
    sl_all = select_slice(g, z_max=FETCH_Z_MAX)
    nz_slice_expected = homogeneous_nz(g["z"][sl_all], zmax_all[sl_all], z_edges)
    nz_slice_observed = np.histogram(g["z"][sl_all], bins=z_edges)[0]

    say("  clustering: pair counts (the slow step; cached afterwards)")
    clus = clustering_measurement(g, vl_in, mask, workers=workers, quiet=quiet)

    say("  the Great Wall and the largest void in the slice")
    ra_mid = 0.5 * (SLICE_RA_RANGE[0] + SLICE_RA_RANGE[1])
    vl_slice = vl & sl
    # Randoms for the slice: same selection as the volume-limited sample, but
    # over the whole surveyed sky, then cut to the slice.
    n_rs = 40 * int(vl_slice.sum())
    ra_s = rng.uniform(*SLICE_RA_RANGE, 4 * n_rs)
    dec_s = np.degrees(np.arcsin(rng.uniform(np.sin(np.radians(-SLICE_DEC_HALF_WIDTH)),
                                             np.sin(np.radians(SLICE_DEC_HALF_WIDTH)), 4 * n_rs)))
    keep = mask.contains(ra_s, dec_s, "surveyed")
    ra_s, dec_s = ra_s[keep][:n_rs], dec_s[keep][:n_rs]
    d_s = cl.random_volume_distances(len(ra_s), float(d_lo), float(d_hi), rng)
    rand_xy = slice_plane(ra_s, d_s, ra_mid)
    wall = great_wall(g, vl_slice, rand_xy)
    void = largest_void(slice_plane(g["ra"][vl_slice], g["d_h"][vl_slice], ra_mid),
                        mask, ra_mid, float(d_lo), float(d_hi))

    # A random catalogue with exactly as many points as the volume-limited
    # slice, for the side-by-side "is it random?" figure.
    pick = rng.choice(len(ra_s), size=int(vl_slice.sum()), replace=False)
    slice_random = {"ra": ra_s[pick], "d_h": d_s[pick]}

    results = {
        "_generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "catalog": meta,
        "n_galaxies": int(len(g["z"])),
        "median_z": float(np.median(g["z"])),
        "median_zerr_kms": float(np.median(g["zerr_kms"])),
        "mask": {"nside": MASK_NSIDE, "surveyed_deg2": mask.area_deg2("surveyed"),
                 "interior_deg2": mask.area_deg2("interior")},
        "slice": {"n": int(sl.sum()), "dec_half_width": SLICE_DEC_HALF_WIDTH,
                  "ra_range": list(SLICE_RA_RANGE), "z_max": WEDGE_Z_MAX,
                  "n_volume_limited": int(vl_slice.sum())},
        "band_shift_offset_mag": g["band_offset"],
        "volume_limited": {
            "n": int(vl.sum()), "n_interior": int(vl_in.sum()),
            "z_range": [VL_Z_MIN, VL_Z_MAX], "mag_limit_lit": VL_MAG_LIMIT_H,
            "number_density_h3": float(n_mean),
            "number_density_zehavi_h3": lit.ZEHAVI_MR20["number_density"],
            "mean_separation_h": mean_sep,
            "completeness": completeness,
            "red_fraction": float(red[vl_in].mean()),
            "red_fraction_zehavi": lit.ZEHAVI_MR20["n_red"] / lit.ZEHAVI_MR20["n_gal"],
        },
        "color_divider": divider,
        "environment": env,
        "groups": {"n_groups_10plus": len(groups), "richest": groups[:15],
                   "coma": next((c for c in groups if c["abell"] == "Abell 1656"), None),
                   "coma_colless_dunn": lit.COLLESS_COMA,
                   "richest_in_slice": in_slice[:8],
                   "linking_lengths_h": [FOF_B_PERP * mean_sep, FOF_B_PAR * mean_sep]},
        "clustering": {
            "n_data": clus["n_data"], "n_random": clus["n_random"],
            "area_deg2": clus["area_deg2"], "volume_h3": clus["volume_h3"],
            "rp": clus["result"].rp.tolist(), "wp": clus["result"].wp.tolist(),
            "wp_err": clus["result"].wp_err.tolist(),
            "fits": clus["fits"],
            "ratio_to_zehavi": clus["ratio_to_zehavi"].tolist(),
            "zehavi": {k: v for k, v in lit.ZEHAVI_MR20.items()},
        },
        "great_wall": {k: v for k, v in wall.items() if k != "map"},
        "largest_void": void,
    }

    art = {
        "g": g, "mask": mask, "slice": sl, "vl": vl, "vl_in": vl_in, "red": red,
        "sigma": sigma, "env_sel": env_sel, "labels": labels, "vl_idx": vl_idx,
        "groups": groups, "groups_in_slice": in_slice,
        "z_edges": z_edges, "nz_observed": nz_observed, "nz_expected": nz_expected,
        "nz_slice_observed": nz_slice_observed, "nz_slice_expected": nz_slice_expected,
        "clustering": clus, "wall": wall, "void": void, "slice_random": slice_random,
        "divider": divider, "env": env, "ra_mid": ra_mid,
    }
    return results, art


def _jsonable(x):
    if isinstance(x, dict):
        return {k: _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def save_results(results: dict) -> None:
    RESULTS_PATH.write_text(json.dumps(_jsonable(results), indent=2))
