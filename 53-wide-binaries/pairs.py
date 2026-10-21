"""Find wide binaries in Gaia DR3 the way Pittordis & Sutherland (2023) did.

Banik et al. (2024) analysed Pittordis & Sutherland's (2023, hereafter PS23)
EDR3 wide-binary sample, which is not public, so we rebuild it from the
Gaia DR3 query in gaia_queries.py (DR3 astrometry is identical to EDR3's):

1. Parent sample: parallax > 10/3 mas, G < 17, |b| > 15 deg (2,101,920 stars
   in PS23).
2. Pairs with projected separation <= 50 kAU at the mean distance, distances
   agreeing within min(8 pc, 4 sigma), and projected velocity difference
   <= 3 km/s with both stars put at the mean distance (92,677 in PS23).
3. Drop every pair that shares a star with another pair (83,884).
4. Drop four open-cluster boxes (PS23 table 1; 83,100).
5. Drop pairs with a third, possibly faint (G < 20) comoving star closer than
   2/3 of the pair's separation (82,402).
6. Arenou et al. (2018) astrometric chi^2 cut on both stars (73,159), applied
   later because it needs columns of the full Gaia rows.

The same search run against a copy of the parent sample shifted by 0.5 deg
on the sky finds only chance alignments, with the same statistics as the
real ones (El-Badry, Rix & Heintz 2021, Section 3.1): that copy calibrates
the line-of-sight contamination of our samples.
"""

from __future__ import annotations

import numpy as np
from astropy.table import Table
from scipy.spatial import cKDTree

from config import DERIVED_DIR, GAIA_DIR, K_PM

# PS23 table 1: open-cluster exclusion boxes (RA, Dec limits in degrees)
CLUSTER_BOXES = {
    "Blanco 1": ((0.0, 2.0), (-31.0, -29.0)),
    "Pleiades": ((54.0, 59.0), (22.0, 27.0)),
    "Praesepe": ((129.0, 131.0), (18.0, 22.0)),
    "Upper Sco": ((238.0, 248.0), (-29.0, -19.0)),
}

PS23 = dict(plx_min=10.0 / 3.0, g_max=17.0, b_min=15.0, rsky_max_au=50_000.0, dd_max_pc=8.0,
            dd_nsigma=4.0, dv_max_kms=3.0, faint_g_max=20.0, faint_plx_nsigma=4.0,
            faint_frac=2.0 / 3.0, faint_min_arcsec=0.5, faint_dv_max_kms=5.0)
SHIFT_DEG = 0.5


def unit_vectors(ra_deg, dec_deg):
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def shifted_positions(ra_deg, dec_deg, shift_deg: float = SHIFT_DEG):
    """Move every star by ``shift_deg`` along a great circle: north in
    declination, wrapping over the pole (ra + 180) for the very few stars
    within ``shift_deg`` of it."""
    dec = dec_deg + shift_deg
    ra = ra_deg.copy()
    over = dec > 90.0
    dec[over] = 180.0 - dec[over]
    ra[over] = (ra[over] + 180.0) % 360.0
    return ra, dec


def load_lean() -> Table:
    t = Table.read(GAIA_DIR / "nearby_lean.fits")
    return t


def _pair_search(xyz_a, xyz_b, plx_a, theta_max_rad, same: bool):
    """Index pairs (i in a, j in b) closer than theta_max_rad[i] on the sky.
    With ``same`` (a is b) each unordered pair is returned once."""
    tree = cKDTree(xyz_b)
    chord = 2.0 * np.sin(np.minimum(theta_max_rad, np.pi) / 2.0)
    lists = tree.query_ball_point(xyz_a, chord, workers=8)
    ii = np.repeat(np.arange(len(lists)), [len(x) for x in lists])
    jj = np.fromiter((j for x in lists for j in x), dtype=np.int64, count=len(ii))
    if same:
        keep = jj > ii
        ii, jj = ii[keep], jj[keep]
    return ii, jj


def pair_quantities(a: dict, b: dict) -> dict:
    """PS23's pair quantities (both stars at the mean distance)."""
    cosang = np.clip(np.sum(a["xyz"] * b["xyz"], axis=1), -1, 1)
    theta = np.arccos(cosang)                                  # rad
    theta_arcsec = np.degrees(theta) * 3600.0
    d_a, d_b = 1000.0 / a["parallax"], 1000.0 / b["parallax"]
    sd_a = d_a * a["parallax_error"] / a["parallax"]
    sd_b = d_b * b["parallax_error"] / b["parallax"]
    dbar = 0.5 * (d_a + d_b)
    dmu = np.hypot(a["pmra"] - b["pmra"], a["pmdec"] - b["pmdec"])
    return dict(theta_arcsec=theta_arcsec, d_a=d_a, d_b=d_b, sd=np.hypot(sd_a, sd_b), dbar=dbar,
                rsky_au=theta_arcsec * dbar, dv_kms=K_PM * dmu * dbar / 1000.0)


def _arrays(t: Table, idx, ra=None, dec=None) -> dict:
    out = {k: np.asarray(t[k][idx]) for k in ("source_id", "parallax", "parallax_error", "pmra",
                                              "pmdec", "phot_g_mean_mag", "b")}
    r = np.asarray(t["ra"][idx]) if ra is None else ra[idx]
    d = np.asarray(t["dec"][idx]) if dec is None else dec[idx]
    out["xyz"] = unit_vectors(r, d)
    out["ra"], out["dec"] = r, d
    return out


def in_cluster_boxes(ra, dec):
    m = np.zeros(len(ra), bool)
    for (r0, r1), (d0, d1) in CLUSTER_BOXES.values():
        m |= (ra >= r0) & (ra <= r1) & (dec >= d0) & (dec <= d1)
    return m


def find_ps23_pairs(t: Table, *, shifted: bool = False, params: dict = PS23, quiet: bool = False):
    """Steps 1-5 above. Returns a dict of pair arrays (star 'a' is the one
    with the smaller index into the lean table) and a cut-flow record."""
    p = params
    plx = np.asarray(t["parallax"])
    G = np.asarray(t["phot_g_mean_mag"])
    b = np.asarray(t["b"])
    parent = (plx > p["plx_min"]) & (G < p["g_max"]) & (np.abs(b) > p["b_min"])
    faint = (plx > p["plx_min"]) & (G < p["faint_g_max"]) & (np.abs(b) > p["b_min"] - 2.0)
    idx_parent = np.flatnonzero(parent)
    flow = {"parent stars": int(parent.sum())}
    ra0, dec0 = np.asarray(t["ra"]), np.asarray(t["dec"])
    if shifted:
        ra_s, dec_s = shifted_positions(ra0, dec0)
    A = _arrays(t, idx_parent)
    Bsrc = _arrays(t, idx_parent, ra_s, dec_s) if shifted else A
    # search radius: 50 kAU at the smallest mean distance a partner allowed by
    # the 8 pc distance cut could give (d - 4 pc); the exact cut on r_sky at
    # the mean distance follows
    d_a = 1000.0 / A["parallax"]
    d_min = np.maximum(d_a - 0.5 * p["dd_max_pc"], 1.0)
    theta_max = np.minimum(np.radians(p["rsky_max_au"] / d_min / 3600.0), np.radians(10.0))
    ii, jj = _pair_search(A["xyz"], Bsrc["xyz"], A["parallax"], theta_max, same=not shifted)
    if shifted:
        # copying and shifting doubles the possible pairs (N(N+1) instead of
        # N(N+1)/2), so keep each with probability 1/2 (El-Badry et al. 2021)
        ii, jj = ii[ii != jj], jj[ii != jj]
        rng = np.random.default_rng(4435257)
        keep = rng.random(len(ii)) < 0.5
        ii, jj = ii[keep], jj[keep]
    a = {k: v[ii] for k, v in A.items()}
    bb = {k: v[jj] for k, v in Bsrc.items()}
    q = pair_quantities(a, bb)
    ok = ((q["rsky_au"] <= p["rsky_max_au"]) & (np.abs(q["d_a"] - q["d_b"]) <= p["dd_max_pc"])
          & (np.abs(q["d_a"] - q["d_b"]) <= p["dd_nsigma"] * q["sd"]) & (q["dv_kms"] <= p["dv_max_kms"])
          & (q["theta_arcsec"] > 0.0))
    ii, jj = ii[ok], jj[ok]
    a = {k: v[ok] for k, v in a.items()}
    bb = {k: v[ok] for k, v in bb.items()}
    q = {k: v[ok] for k, v in q.items()}
    flow["candidate pairs (step 2)"] = int(ok.sum())
    # step 3: no star in more than one pair (a shifted copy is a different
    # 'star', so the shared-star test runs over the original indices)
    allidx = np.concatenate([idx_parent[ii], idx_parent[jj] + (len(t) if shifted else 0)])
    u, c = np.unique(allidx, return_counts=True)
    multi = set(u[c > 1].tolist())
    lone = np.array([(x not in multi) and (y not in multi) for x, y in
                     zip(idx_parent[ii], idx_parent[jj] + (len(t) if shifted else 0))], bool)
    flow["no shared star (step 3)"] = int(lone.sum())
    # step 4: open-cluster boxes (either star)
    box = in_cluster_boxes(a["ra"], a["dec"]) | in_cluster_boxes(bb["ra"], bb["dec"])
    keep = lone & ~box
    flow["outside cluster boxes (step 4)"] = int(keep.sum())
    # step 5: faint comoving third stars
    sel = np.flatnonzero(keep)
    F = _arrays(t, np.flatnonzero(faint))
    third = _faint_companions(a, bb, q, sel, F, p, shifted=shifted, idx_parent=idx_parent,
                              ii=ii, jj=jj, faint_idx=np.flatnonzero(faint))
    keep[sel[third]] = False
    flow["no faint third star (step 5)"] = int(keep.sum())
    if not quiet:
        print("   ", "shifted" if shifted else "real", flow, flush=True)
    out = dict(source_id_a=a["source_id"], source_id_b=bb["source_id"], idx_a=idx_parent[ii],
               idx_b=idx_parent[jj], **{k: v for k, v in q.items()}, keep=keep,
               lone=lone, in_box=box)
    return out, flow


def _faint_companions(a, bb, q, sel, F, p, *, shifted, idx_parent, ii, jj, faint_idx):
    """For pairs ``sel``, is there a faint-sample star (other than the two
    members) within 2/3 of the pair separation of either member, beyond 0.5",
    with parallax within 4 sigma and projected velocity difference <= 5 km/s?
    For the shifted catalogue the companion search runs around each member in
    its own sky frame (member b lives in the shifted frame)."""
    treeF = cKDTree(F["xyz"])
    hits = np.zeros(len(sel), bool)
    frac = p["faint_frac"]
    for star, gidx in (("a", idx_parent[ii]), ("b", idx_parent[jj])):
        S = a if star == "a" else bb
        xyz = S["xyz"][sel]
        if shifted and star == "b":
            # look around the member's true position (faint stars are unshifted)
            xyz = unit_vectors(*_unshift(S["ra"][sel], S["dec"][sel]))
        rad = np.radians(frac * q["theta_arcsec"][sel] / 3600.0)
        lists = treeF.query_ball_point(xyz, 2.0 * np.sin(rad / 2.0), workers=8)
        for k, lst in enumerate(lists):
            if not lst:
                continue
            j = np.asarray(lst)
            fid = faint_idx[j]
            others = (fid != idx_parent[ii][sel[k]]) & (fid != idx_parent[jj][sel[k]])
            j = j[others]
            if not len(j):
                continue
            ang = np.degrees(np.arccos(np.clip(F["xyz"][j] @ xyz[k], -1, 1))) * 3600.0
            plx_s, eplx_s = S["parallax"][sel[k]], None
            eplx_s = _plx_err(S, sel[k])
            dplx = np.abs(F["parallax"][j] - plx_s)
            ok_plx = dplx <= p["faint_plx_nsigma"] * np.hypot(F["parallax_error"][j], eplx_s)
            dmu = np.hypot(F["pmra"][j] - S["pmra"][sel[k]], F["pmdec"][j] - S["pmdec"][sel[k]])
            dv = K_PM * dmu / plx_s
            if np.any((ang > p["faint_min_arcsec"]) & ok_plx & (dv <= p["faint_dv_max_kms"])):
                hits[k] = True
    return hits


def _plx_err(S, k):
    return S["parallax_error"][k]


def _unshift(ra, dec, shift_deg: float = SHIFT_DEG):
    dec2 = dec - shift_deg
    ra2 = ra.copy()
    under = dec2 < -90.0
    dec2[under] = -180.0 - dec2[under]
    ra2[under] = (ra2[under] + 180.0) % 360.0
    return ra2, dec2


def cached_pairs(shifted: bool = False, quiet: bool = False):
    """find_ps23_pairs() on the cached lean table, saved as
    data/derived/ps23_pairs[_shifted].npz with its cut flow."""
    import json
    name = "ps23_pairs_shifted" if shifted else "ps23_pairs"
    path = DERIVED_DIR / f"{name}.npz"
    if path.exists():
        z = np.load(path)
        return {k: z[k] for k in z.files}, json.loads((DERIVED_DIR / f"{name}_flow.json").read_text())
    t = load_lean()
    out, flow = find_ps23_pairs(t, shifted=shifted, quiet=quiet)
    np.savez_compressed(path, **out)
    (DERIVED_DIR / f"{name}_flow.json").write_text(json.dumps(flow, indent=1))
    return out, flow


if __name__ == "__main__":
    import time
    t0 = time.time()
    real, flow = cached_pairs(False)
    print(flow, f"{time.time() - t0:.0f} s")
    t0 = time.time()
    sh, flow_s = cached_pairs(True)
    print(flow_s, f"{time.time() - t0:.0f} s")
