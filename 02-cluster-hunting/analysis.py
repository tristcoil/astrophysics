"""The cluster-hunting pipeline.

Given a cached Gaia DR3 cone around a cluster, we:

1. **Prepare** the sample: correct Gaia's parallax zero point, flag
   contaminated photometry, project positions onto a flat map of the sky.
2. **Discover** the cluster as the most concentrated knot in proper-motion
   space. Its motion is something we find, never something we look up.
3. **Correct for perspective.** Stars sharing one space velocity show
   slightly different proper motions across a field this wide; we move every
   star into the cluster's own frame.
4. **Model membership** with extreme deconvolution, twice: once from proper
   motions alone ("kinematics alone"), once adding parallax.
5. **Confirm** with evidence the model never saw. Position on the sky and
   colour/brightness are never used to select members, so a concentrated sky
   distribution and a thin HR-diagram sequence are genuine results rather than
   consequences of the selection.
6. **Measure** the cluster: distance, 3D space velocity, size, depth, internal
   velocity spread, and the width of its main sequence against the field's.
7. **Grade** everything against Hunt & Reffert (2023) and the Gaia
   Collaboration (2018), consulted only at this last step.

The output is a JSON results file, plus in-memory tables the figure code
reuses, so plotting never repeats the analysis.
"""

from __future__ import annotations

import json
import warnings
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.signal import fftconvolve

import gaia_physics as gp
import membership as mb
from config import (
    CSTAR_G_MIN, CSTAR_NSIGMA, FIELD_PARALLAX_OVER_ERROR, FIT_REGION_MASYR,
    HALO_MAX_DEPTH_FRACTION, KNOT_ANNULUS_MASYR,
    KNOT_APERTURE_MASYR, KNOT_BIN_MASYR, KNOT_CANDIDATES, KNOT_MAX_PARALLAX_SPREAD,
    KNOT_MIN_STARS, KNOT_SEARCH_HALFWIDTH_MASYR,
    MEMBERSHIP_THRESHOLD, PARALLAX_SYSTEMATIC_MAS, PERSPECTIVE_ROUNDS, RESULTS_PATH,
    RUWE_GOOD, SEED_RADIUS_MASYR, WIDTH_BP_RP_RANGE, WIDTH_PARALLAX_OVER_ERROR,
    WIDTH_WINDOW_MAG, XD_MAX_ITER, XD_REGULARISATION, XD_TOL, Cluster,
)

warnings.filterwarnings("ignore", category=RuntimeWarning)


def _robust_sigma(x) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if len(x) else float("nan")


# --------------------------------------------------------------------------
# Step 1 -- prepare
# --------------------------------------------------------------------------

def prepare(df: pd.DataFrame, cluster: Cluster) -> pd.DataFrame:
    """Add the derived columns every later step uses."""
    df = df.copy()
    zp, extrapolated = gp.parallax_zero_point(df)
    df["zero_point"] = zp
    df["zp_extrapolated"] = extrapolated
    df["parallax_corr"] = df["parallax"] - zp
    df["sep_deg"] = gp.angular_separation(df["ra"], df["dec"], cluster.ra_deg, cluster.dec_deg)
    df["xi"], df["eta"] = gp.gnomonic(df["ra"], df["dec"], cluster.ra_deg, cluster.dec_deg)
    df["clean_phot"] = gp.clean_photometry(
        df["phot_bp_rp_excess_factor"], df["bp_rp"], df["phot_g_mean_mag"],
        nsigma=CSTAR_NSIGMA, g_min=CSTAR_G_MIN,
    )
    df["ruwe_ok"] = df["ruwe"] < RUWE_GOOD
    return df


# --------------------------------------------------------------------------
# Step 2 -- find the knot
# --------------------------------------------------------------------------

def _disk_and_ring(bin_size: float):
    r_ap = KNOT_APERTURE_MASYR / bin_size
    r_in, r_out = (r / bin_size for r in KNOT_ANNULUS_MASYR)
    n = int(np.ceil(r_out))
    yy, xx = np.mgrid[-n:n + 1, -n:n + 1]
    rr = np.hypot(xx, yy)
    return (rr <= r_ap).astype(float), ((rr >= r_in) & (rr <= r_out)).astype(float)


def discover_knot(pmra, pmdec, parallax) -> dict:
    """Find the cluster as a clump of shared motion. Nothing is assumed.

    The method is the one used to find a star in an image. At every point of
    proper-motion space, count the stars inside a small circle (the
    "aperture"), and estimate from a ring around it how many a smooth field
    would have put there. The ratio is the *contrast*, and the few highest-
    contrast places are the candidates.

    Motion alone finds two kinds of clump. One is a cluster: stars sharing a
    velocity. The other is the field's own core -- the pile-up of distant
    stars whose proper motions are all small simply because they are far away.
    In the Pleiades field that pile-up reaches a contrast of ~9. What tells
    them apart is distance: a cluster's stars share a parallax, the pile-up's
    span a factor of several. So each candidate's parallaxes are checked, and
    the highest-contrast candidate whose stars share a distance is the knot.
    Motion proposes; distance confirms.
    """
    pmra, pmdec, parallax = (np.asarray(a, dtype=float) for a in (pmra, pmdec, parallax))
    b, h = KNOT_BIN_MASYR, KNOT_SEARCH_HALFWIDTH_MASYR
    cx, cy = float(np.median(pmra)), float(np.median(pmdec))
    ex = np.arange(cx - h, cx + h + b / 2, b)
    ey = np.arange(cy - h, cy + h + b / 2, b)
    counts, _, _ = np.histogram2d(pmra, pmdec, bins=[ex, ey])

    disk, ring = _disk_and_ring(b)
    on = np.clip(np.rint(fftconvolve(counts, disk, mode="same")), 0, None)
    off = np.clip(fftconvolve(counts, ring, mode="same"), 0, None)
    expected = np.maximum(off * disk.sum() / ring.sum(), 1.0)
    contrast = on / expected
    score = np.where(on >= KNOT_MIN_STARS, contrast, 0.0)
    xc, yc = 0.5 * (ex[1:] + ex[:-1]), 0.5 * (ey[1:] + ey[:-1])

    candidates = []
    remaining = score.copy()
    for _ in range(KNOT_CANDIDATES):
        i, j = np.unravel_index(np.argmax(remaining), remaining.shape)
        if remaining[i, j] <= 0:
            break
        cand = _describe(i, j, xc, yc, on, expected, contrast)
        # Re-centre on the stars themselves: the median motion in the aperture.
        pa, pd_ = cand["pmra"], cand["pmdec"]
        for _ in range(5):
            near = np.hypot(pmra - pa, pmdec - pd_) < KNOT_APERTURE_MASYR
            pa, pd_ = float(np.median(pmra[near])), float(np.median(pmdec[near]))
        plx_med = float(np.median(parallax[near]))
        spread = _robust_sigma(parallax[near]) / abs(plx_med) if plx_med else np.inf
        cand.update(pmra=pa, pmdec=pd_, n_aperture=int(near.sum()),
                    parallax_median=plx_med, parallax_spread_fraction=float(spread),
                    shares_distance=bool(spread < KNOT_MAX_PARALLAX_SPREAD))
        candidates.append(cand)
        remaining[np.hypot(xc[:, None] - xc[i], yc[None, :] - yc[j])
                  < 2 * KNOT_ANNULUS_MASYR[1]] = 0.0

    coherent = [c for c in candidates if c["shares_distance"]]
    if not coherent:
        raise RuntimeError("no clump in proper-motion space has stars sharing a distance")
    best = coherent[0]
    others = [c for c in candidates if c is not best]
    return {
        **best,
        "rank_by_contrast": candidates.index(best) + 1,
        "runner_up": others[0] if others else None,
        "candidates": candidates,
        "field_median_pm": (cx, cy),
        "map_counts": counts,
        "map_contrast": np.where(on >= KNOT_MIN_STARS, contrast, np.nan),
        "map_extent": (ex[0], ex[-1], ey[0], ey[-1]),
    }


def _describe(i, j, xc, yc, on, expected, contrast) -> dict:
    n, lam = float(on[i, j]), float(expected[i, j])
    # Poisson likelihood-ratio significance of n counts where lam were expected
    # (the "Cash" statistic) -- valid far into the tail, unlike sqrt(n).
    sig = float(np.sqrt(max(2.0 * (n * np.log(n / lam) - (n - lam)), 0.0))) if n > lam else 0.0
    return {"pmra": float(xc[i]), "pmdec": float(yc[j]), "count": int(n),
            "expected": lam, "contrast": float(contrast[i, j]), "significance": sig}


# --------------------------------------------------------------------------
# Step 3 -- perspective
# --------------------------------------------------------------------------

def robust_rv(df: pd.DataFrame, mask) -> dict:
    """Mean radial velocity of the selected stars, resistant to binaries.

    Gaia measured radial velocities for the brighter stars (G <~ 15). A
    member in a close binary shows its orbital velocity on top of the
    cluster's, so a plain mean is dragged around by a few extreme values.
    Iterative 3-sigma clipping around the median removes them.

    Two quality steps first. Stars with hot spectral templates (>= 8500 K)
    are left out: their velocities are less precise and follow a different
    magnitude trend (Blomme et al. 2023). The rest get the Katz et al. (2023)
    correction for fainter stars.
    """
    sub = df.loc[mask]
    rv = (sub["radial_velocity"] - gp.rv_magnitude_correction(sub["grvs_mag"])).to_numpy(float)
    err = sub["radial_velocity_error"].to_numpy(float)
    cool = (sub["rv_template_teff"] < 8500).to_numpy()
    ok = np.isfinite(rv) & np.isfinite(err) & (err < 5.0) & cool
    rv = rv[ok]
    if len(rv) < 5:
        return {"rv": float("nan"), "err": float("nan"), "n": int(len(rv)), "n_used": 0}
    keep = np.ones(len(rv), bool)
    for _ in range(10):
        med, sig = np.median(rv[keep]), _robust_sigma(rv[keep])
        new = np.abs(rv - med) < 3 * sig
        if (new == keep).all():
            break
        keep = new
    sel = rv[keep]
    return {"rv": float(np.median(sel)),
            "err": float(1.2533 * np.std(sel) / np.sqrt(len(sel))),
            "scatter": float(np.std(sel)), "n": int(len(rv)), "n_used": int(keep.sum())}


def apply_frame(df, cluster: Cluster, pmra0, pmdec0, plx0, rv0):
    """Move every star's proper motion into the cluster's frame."""
    velocity = gp.space_velocity(cluster.ra_deg, cluster.dec_deg, pmra0, pmdec0, plx0, rv0)
    dpa, dpd = gp.perspective_offset(df["ra"].to_numpy(), df["dec"].to_numpy(),
                                     cluster.ra_deg, cluster.dec_deg, velocity, plx0)
    df["dpmra_persp"], df["dpmdec_persp"] = dpa, dpd
    df["pmra_c"] = df["pmra"] - dpa
    df["pmdec_c"] = df["pmdec"] - dpd
    return velocity


# --------------------------------------------------------------------------
# Step 4 -- membership
# --------------------------------------------------------------------------

def _design(df, with_parallax: bool):
    cols = ["pmra_c", "pmdec_c"] + (["parallax_corr"] if with_parallax else [])
    return df[cols].to_numpy(float), gp.astrometric_covariance(df, with_parallax=with_parallax)


def fit_membership(df, region, seed, n_field: int, *, with_parallax: bool, init=None):
    """Fit the mixture to the stars in ``region``; everyone else gets p = 0."""
    x, s = _design(df[region], with_parallax)
    if init is None:
        # A halo is only identifiable with parallax, so the motion-only model
        # describes the cluster with a single Gaussian.
        init = mb.initialise(x, seed[region], n_field, halo=with_parallax)
    mix = mb.fit(x, s, init, max_iter=XD_MAX_ITER, tol=XD_TOL,
                 regularisation=XD_REGULARISATION)
    mix.cluster = mb.cluster_components(mix, max_depth_fraction=HALO_MAX_DEPTH_FRACTION)
    p = np.zeros(len(df))
    p[region] = mb.membership_probability(x, s, mix)
    return mix, p


def _check_cluster_component(mix, knot, label):
    """The cluster Gaussian must stay compact and on the knot, or the fit is suspect."""
    c = mix.core
    widths = [np.linalg.det(cv[:2, :2]) for cv in mix.covs]
    drift = float(np.hypot(mix.means[c][0] - knot["pmra"], mix.means[c][1] - knot["pmdec"]))
    if np.argmin(widths) != c or drift > 3.0:
        print(f"    WARNING [{label}]: cluster component drifted {drift:.1f} mas/yr "
              f"or is not the most compact -- inspect the fit")


# --------------------------------------------------------------------------
# Step 6 -- measurements
# --------------------------------------------------------------------------

def ridge_line(color, mag, lo, hi, step=0.1, min_stars=5):
    """Median absolute magnitude in colour bins: the spine of a sequence."""
    edges = np.arange(lo, hi + step / 2, step)
    centres, medians = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (color >= a) & (color < b)
        if sel.sum() >= min_stars:
            centres.append(0.5 * (a + b))
            medians.append(float(np.median(mag[sel])))
    return np.array(centres), np.array(medians)


def sequence_width(color, mag, ridge_c, ridge_m):
    """Robust scatter of a sample about a ridge line, in magnitudes.

    Stars more than WIDTH_WINDOW_MAG from the ridge are left out, which keeps
    giants out of a main-sequence statistic. The spread is measured around the
    sample's own median offset, so a field that sits systematically off the
    cluster's ridge is not penalised for the offset -- only for its spread.
    """
    inside = (color >= ridge_c.min()) & (color <= ridge_c.max())
    resid = mag[inside] - np.interp(color[inside], ridge_c, ridge_m)
    resid = resid[np.abs(resid) < WIDTH_WINDOW_MAG]
    return {"sigma": _robust_sigma(resid), "offset": float(np.median(resid)),
            "n": int(len(resid)), "residuals": resid}


def depth_and_motion(df, members, plx0) -> dict:
    """The cluster in three dimensions: its depth, and whether nearer members move faster.

    **Depth.** A cluster's extent along our line of sight shows up as a spread
    of parallaxes. Whether that spread is visible at all depends on distance:
    if the cluster is roughly round, its depth should match its width on the
    sky, which fixes the parallax spread to expect. That is compared with the
    noise -- the members' parallax errors plus Gaia's ~0.01 mas small-scale
    systematics. For the Pleiades the expected spread is ~7x the noise; for
    M67, six times further away, it is ~6x *smaller* than the noise, and any
    spread measured there is noise, not depth.

    **Motion.** A star's proper motion is its transverse velocity divided by
    its distance: mu = v_t * parallax / K. If all members share one velocity,
    a member 3% nearer shows 3% more proper motion, so the slope
    d(mu)/d(parallax) must equal mu/parallax. Fitted by least squares on the
    members with the best parallaxes, and divided by the "reliability ratio"
    (real spread / observed spread) -- the standard correction for a fit whose
    x-values are themselves noisy, which otherwise flattens the slope.

    Width and depth use the same stars and the same robust statistic, so the
    comparison between them is fair.
    """
    sel = (members & (df["parallax_corr"] / df["parallax_error"] > 50) & df["ruwe_ok"]).to_numpy()
    if sel.sum() < 30:
        return {"resolvable": False, "n": int(sel.sum()), "reason": "too few precise members"}
    d_pc = 1000.0 / plx0
    w = df.loc[sel, "parallax_corr"].to_numpy(float)
    e = df.loc[sel, "parallax_error"].to_numpy(float)
    keep = np.abs(w - np.median(w)) < 4 * _robust_sigma(w)
    width_deg = np.sqrt(0.5 * (_robust_sigma(df.loc[sel, "xi"]) ** 2
                               + _robust_sigma(df.loc[sel, "eta"]) ** 2))
    width_pc = float(np.radians(width_deg) * d_pc)
    expected = plx0 * width_pc / d_pc
    noise = float(np.sqrt(np.median(e[keep] ** 2) + PARALLAX_SYSTEMATIC_MAS ** 2))
    out = {"n": int(keep.sum()), "width_sigma_pc": width_pc,
           "expected_parallax_spread_mas": float(expected), "parallax_noise_mas": noise,
           "resolvable": bool(expected > 2 * noise)}
    if not out["resolvable"]:
        out["reason"] = (f"depth hidden by noise: a round cluster would spread parallaxes by "
                         f"{expected:.4f} mas, against {noise:.3f} mas of noise")
        return out
    real = max(_robust_sigma(w[keep]) ** 2 - np.median(e[keep] ** 2), 0.0)
    out["depth_sigma_pc"] = float(np.sqrt(real) / plx0 * d_pc)
    reliability = float(1.0 - np.mean(e[keep] ** 2) / np.var(w[keep]))
    out["reliability"] = reliability
    for comp in ("pmra_c", "pmdec_c"):
        y = df.loc[sel, comp].to_numpy(float)[keep]
        coef, cov = np.polyfit(w[keep] - plx0, y, 1, cov=True)
        out[comp] = {"slope": float(coef[0] / reliability),
                     "err": float(np.sqrt(cov[0, 0]) / reliability)}
    # The same test as one number: proper motion along the cluster's direction
    # of travel, whose predicted slope is simply |mu| / parallax.
    mu = np.array([df.loc[sel, "pmra_c"].median(), df.loc[sel, "pmdec_c"].median()])
    along = (df.loc[sel, ["pmra_c", "pmdec_c"]].to_numpy(float) @ (mu / np.hypot(*mu)))[keep]
    coef, cov = np.polyfit(w[keep] - plx0, along, 1, cov=True)
    out["along_motion"] = {"slope": float(coef[0] / reliability),
                           "err": float(np.sqrt(cov[0, 0]) / reliability),
                           "predicted": float(np.hypot(*mu) / plx0)}
    return out


def perspective_cells(df, members, radius_deg, pm0, min_stars=12) -> dict:
    """Test the perspective geometry on the members' measured motions.

    The field is cut into cells; in each, the members' median offset from the
    centre's motion is *measured*, and the offset the shared 3D velocity
    *predicts* is averaged. If the geometry is right, subtracting the
    prediction should leave only cell-to-cell noise. Reported as RMS offsets
    before and after, beside the noise expected from each cell's own scatter.
    Uses members with RUWE < 1.4, whose motions are trustworthy.
    """
    sel0 = (members & df["ruwe_ok"]).to_numpy()
    xi, eta = df["xi"].to_numpy(), df["eta"].to_numpy()
    dx = df["pmra"].to_numpy() - pm0[0]
    dy = df["pmdec"].to_numpy() - pm0[1]
    step = radius_deg / 3.5
    edges = np.arange(-radius_deg, radius_deg + step, step)
    cells = []
    for x0 in edges[:-1]:
        for y0 in edges[:-1]:
            s = sel0 & (xi >= x0) & (xi < x0 + step) & (eta >= y0) & (eta < y0 + step)
            n = int(s.sum())
            if n < min_stars:
                continue
            obs = np.array([np.median(dx[s]), np.median(dy[s])])
            err = 1.2533 * np.array([_robust_sigma(dx[s]), _robust_sigma(dy[s])]) / np.sqrt(n)
            pred = np.array([df["dpmra_persp"].to_numpy()[s].mean(),
                             df["dpmdec_persp"].to_numpy()[s].mean()])
            cells.append({"xi": float(xi[s].mean()), "eta": float(eta[s].mean()), "n": n,
                          "observed": obs.tolist(), "predicted": pred.tolist(),
                          "error": err.tolist()})
    if len(cells) < 3:
        return {"cells": cells, "n_cells": len(cells)}
    obs = np.array([c["observed"] for c in cells])
    pred = np.array([c["predicted"] for c in cells])
    err = np.array([c["error"] for c in cells])
    return {
        "cells": cells, "n_cells": len(cells), "step_deg": float(step),
        "rms_offset_masyr": float(np.sqrt(np.mean(np.sum(obs ** 2, axis=1)))),
        "rms_after_model_masyr": float(np.sqrt(np.mean(np.sum((obs - pred) ** 2, axis=1)))),
        "rms_cell_noise_masyr": float(np.sqrt(np.mean(np.sum(err ** 2, axis=1)))),
        "chi2_per_dof_after_model": float(np.mean(((obs - pred) / err) ** 2)),
    }


def radial_profile(df, members, p, centre, radius_deg, n_bins=12):
    """Surface density of members in rings, and the field-contamination floor."""
    sep = gp.angular_separation(df["ra"].to_numpy(), df["dec"].to_numpy(), *centre)
    edges = np.linspace(0, radius_deg, n_bins + 1)
    area = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    n_mem, _ = np.histogram(sep[members], bins=edges)
    n_all, _ = np.histogram(sep, bins=edges)
    contam = float(np.sum(1 - p[members])) / (np.pi * radius_deg ** 2)
    return {"r_mid": 0.5 * (edges[1:] + edges[:-1]), "edges": edges,
            "density": n_mem / area, "density_err": np.sqrt(n_mem) / area,
            "field_density": n_all / area, "contamination_density": contam,
            "r50_deg": float(np.median(sep[members]))}


# --------------------------------------------------------------------------
# Results
# --------------------------------------------------------------------------

@dataclass
class ClusterResult:
    cluster: str
    key: str
    difficulty: str
    n_stars: int
    knot: dict
    perspective: dict
    kinematic_only: dict
    membership: dict
    motion: dict
    distance: dict
    structure: dict
    hr_diagram: dict
    famous_stars: dict
    comparison: dict = field(default_factory=dict)
    membership_vs_published: dict = field(default_factory=dict)


def analyse(cluster: Cluster, df_raw: pd.DataFrame, literature: dict | None = None,
            published_members: pd.DataFrame | None = None,
            named: pd.DataFrame | None = None, *, quiet: bool = False):
    """Run the whole pipeline on one cluster. Returns (result, artefacts)."""
    say = (lambda *a: None) if quiet else print
    say(f"\n=== {cluster.name} ({cluster.difficulty} case) ===")

    # 1. prepare ----------------------------------------------------------
    df = prepare(df_raw, cluster)
    say(f"  {len(df):,} Gaia DR3 stars with G < {cluster.g_max} within "
        f"{cluster.radius_deg} deg")

    # 2. discover ---------------------------------------------------------
    knot = discover_knot(df["pmra"], df["pmdec"], df["parallax_corr"])
    say(f"  knot found at ({knot['pmra']:+.2f}, {knot['pmdec']:+.2f}) mas/yr: "
        f"{knot['count']} stars where {knot['expected']:.1f} expected "
        f"(contrast {knot['contrast']:.0f}x), parallaxes agree to "
        f"{knot['parallax_spread_fraction']:.0%}")
    for cand in knot["candidates"]:
        if cand["pmra"] != knot["pmra"]:
            kind = "shares a distance" if cand["shares_distance"] else "does not share a distance"
            say(f"    other candidate ({cand['pmra']:+.1f}, {cand['pmdec']:+.1f}): contrast "
                f"{cand['contrast']:.1f}x, parallax spread "
                f"{cand['parallax_spread_fraction']:.0%} -- {kind}")

    # Seed: stars near the knot in motion whose parallaxes agree with each
    # other. Used only to start the model and set the perspective geometry.
    near = np.hypot(df["pmra"] - knot["pmra"], df["pmdec"] - knot["pmdec"]) < SEED_RADIUS_MASYR
    plx_near = df.loc[near, "parallax_corr"]
    seed = (near & (np.abs(df["parallax_corr"] - plx_near.median())
                    < 3 * _robust_sigma(plx_near))).to_numpy()
    pm0 = np.array([knot["pmra"], knot["pmdec"]])
    plx0 = float(np.median(df.loc[seed, "parallax_corr"]))
    rv = robust_rv(df, seed)
    region = (np.hypot(df["pmra"] - knot["pmra"], df["pmdec"] - knot["pmdec"])
              < FIT_REGION_MASYR).to_numpy()
    df["fit_region"] = region
    say(f"  modelling the {region.sum():,} stars within {FIT_REGION_MASYR:.0f} mas/yr "
        f"of the knot")

    # 3 + 4. perspective and membership, iterated ---------------------------
    mix3 = None
    for rnd in range(PERSPECTIVE_ROUNDS):
        rv0 = rv["rv"] if np.isfinite(rv["rv"]) else 0.0
        velocity = apply_frame(df, cluster, pm0[0], pm0[1], plx0, rv0)
        mix3, p3 = fit_membership(df, region, seed, cluster.n_field_components,
                                  with_parallax=True, init=mix3)
        _check_cluster_component(mix3, knot, "motion+parallax")
        members = p3 > MEMBERSHIP_THRESHOLD
        c = mix3.core
        pm0, plx0 = mix3.means[c][:2].copy(), float(mix3.means[c][2])
        rv = robust_rv(df, members)
        say(f"  round {rnd + 1}: {members.sum()} members, {mix3.n_iter} iterations; "
            f"pm = ({pm0[0]:+.3f}, {pm0[1]:+.3f}) mas/yr, parallax {plx0:.4f} mas, "
            f"RV {rv['rv']:+.2f} km/s from {rv['n_used']} stars")

    # Final frame from the final parameters, then the "motion only" model in it.
    velocity = apply_frame(df, cluster, pm0[0], pm0[1], plx0, rv["rv"])
    mix3, p3 = fit_membership(df, region, seed, cluster.n_field_components,
                              with_parallax=True, init=mix3)
    mix2, p2 = fit_membership(df, region, seed, cluster.n_field_components,
                              with_parallax=False)
    _check_cluster_component(mix2, knot, "motion only")
    df["p_motion"] = p2
    df["p_member"] = p3
    members = (p3 > MEMBERSHIP_THRESHOLD)
    df["member"] = members
    kin = p2 > MEMBERSHIP_THRESHOLD
    df["kinematic_candidate"] = kin
    # Report the final model's own numbers. (The frame it was fitted in came
    # from the previous round; once converged they differ by ~0.001 mas/yr.)
    c = mix3.core
    pm0, plx0 = mix3.means[c][:2].copy(), float(mix3.means[c][2])
    rv = robust_rv(df, members)
    velocity = gp.space_velocity(cluster.ra_deg, cluster.dec_deg, pm0[0], pm0[1], plx0, rv["rv"])

    x3, s3 = _design(df[region], with_parallax=True)
    df["parallax_post"] = np.nan
    df.loc[region, "parallax_post"] = mb.posterior_mean(x3, s3, mix3)[:, 2]
    say(f"  motion only: {kin.sum()} candidates; motion + parallax: {members.sum()} members")

    persp_members = np.hypot(df.loc[members, "dpmra_persp"], df.loc[members, "dpmdec_persp"])
    perspective = {
        "rms_members_masyr": float(np.sqrt(np.mean(persp_members ** 2))),
        "max_members_masyr": float(persp_members.max()),
        "max_field_masyr": float(np.hypot(df["dpmra_persp"], df["dpmdec_persp"]).max()),
    }
    cells = perspective_cells(df, members, cluster.radius_deg, pm0)
    perspective.update({k: v for k, v in cells.items() if k != "cells"})

    # 5. confirmations ----------------------------------------------------
    # (a) Stars picked by motion alone: does adding parallax confirm them?
    # Those it rejects are field stars that happen to share the cluster's
    # motion while sitting at a different distance.
    cov3 = mix3.covs[c]
    kinematic_only = {
        "n": int(kin.sum()),
        "expected_contamination": float(np.sum(1 - p2[kin])),
        "n_confirmed_by_parallax": int((kin & members).sum()),
        "n_rejected_by_parallax": int((kin & ~members).sum()),
        "n_members_missed_by_motion_alone": int((members & ~kin).sum()),
        "pm_dispersion_masyr": [float(np.sqrt(mix2.covs[mix2.core][i, i])) for i in (0, 1)],
    }

    membership = {
        "n": int(members.sum()),
        "expected_contamination": float(np.sum(1 - p3[members])),
        "expected_missed": float(np.sum(p3[~members])),
        "n_clear": int((p3 > 0.9).sum()),
        "n_ambiguous": int(((p3 > 0.1) & (p3 < 0.9)).sum()),
        "converged": bool(mix3.converged),
        "iterations": int(mix3.n_iter),
        "has_halo": len(mix3.cluster) > 1,
    }
    n_region = int(region.sum())
    membership["n_core"] = float(mix3.weights[c] * n_region)
    if len(mix3.cluster) > 1:
        halo = list(mix3.cluster[1:])
        membership["n_halo"] = float(mix3.weights[halo].sum() * n_region)
        # Weighted mean parallax spread of the halo part(s), as a depth.
        halo_plx_sd = float(np.average([np.sqrt(mix3.covs[j][2, 2]) for j in halo],
                                       weights=mix3.weights[halo]))
        membership["halo_parallax_sd_mas"] = halo_plx_sd
        membership["halo_depth_sigma_pc"] = float(1000.0 / plx0 * halo_plx_sd / plx0)

    # 6. measurements -----------------------------------------------------
    mean_cov = mb.mean_uncertainty(s3, mix3, p3[region])
    err = np.sqrt(np.diag(mean_cov))
    plx_raw = float(np.average(df.loc[members, "parallax"], weights=p3[members]))
    zp_mean = float(np.average(df.loc[members, "zero_point"], weights=p3[members]))
    plx_err_total = float(np.hypot(err[2], PARALLAX_SYSTEMATIC_MAS))
    d_pc = 1000.0 / plx0
    distance = {
        "parallax_mas": plx0,
        "parallax_err_stat": float(err[2]),
        "parallax_err_sys": PARALLAX_SYSTEMATIC_MAS,
        "parallax_uncorrected_mas": plx_raw,
        "zero_point_mean_mas": zp_mean,
        "distance_pc": d_pc,
        "distance_err_pc": d_pc * plx_err_total / plx0,
        "distance_modulus": float(gp.distance_modulus(d_pc)),
    }

    # Intrinsic spread. The 3D Gaussian knows nearer members move faster across
    # the sky; conditioning on parallax removes that depth effect and leaves
    # the internal velocity spread.
    s_pm = cov3[:2, :2] - np.outer(cov3[:2, 2], cov3[2, :2]) / cov3[2, 2]
    sigma_pm = np.sqrt(np.diag(s_pm))
    sigma_v = gp.K_KMS * float(np.sqrt(np.mean(np.diag(s_pm)))) / plx0
    uvw = gp.galactic_velocity(velocity)
    cp = gp.convergent_point(velocity)
    motion = {
        # The model's motion refers to the cone centre, in the cluster frame.
        "pmra": float(pm0[0]), "pmdec": float(pm0[1]),
        "pmra_err": float(err[0]), "pmdec_err": float(err[1]),
        # Catalogues usually quote the plain average of members' measured
        # motions instead, which refers to wherever the members are centred.
        # Across the Pleiades that differs by ~0.05 mas/yr, so both are kept.
        "pmra_members_mean": float(np.average(df.loc[members, "pmra"], weights=p3[members])),
        "pmdec_members_mean": float(np.average(df.loc[members, "pmdec"], weights=p3[members])),
        "rv_kms": rv["rv"], "rv_err_kms": rv["err"], "rv_n": rv["n_used"],
        "rv_scatter_kms": rv.get("scatter", float("nan")),
        "velocity_icrs_kms": [float(v) for v in velocity],
        "speed_kms": float(np.linalg.norm(velocity)),
        "uvw_kms": [float(v) for v in uvw],
        "convergent_point_deg": list(cp),
        "tangential_speed_kms": float(gp.tangential_velocity(np.hypot(*pm0), plx0)),
        "pm_dispersion_masyr": [float(v) for v in sigma_pm],
        "velocity_dispersion_kms": sigma_v,
        "predicted_slope": [float(pm0[0] / plx0), float(pm0[1] / plx0)],
    }
    three_d = depth_and_motion(df, members, plx0)
    motion["depth_slope"] = {k: three_d[k] for k in ("pmra_c", "pmdec_c", "along_motion",
                                                     "reliability") if k in three_d}

    centre = (float(np.median(df.loc[members, "ra"])), float(np.median(df.loc[members, "dec"])))
    profile = radial_profile(df, members, p3, centre, cluster.radius_deg)
    structure = {
        "centre_deg": list(centre),
        "r50_deg": profile["r50_deg"],
        "r50_pc": float(np.radians(profile["r50_deg"]) * d_pc),
        "core_parallax_dispersion_mas": float(np.sqrt(cov3[2, 2])),
        "depth_resolvable": three_d["resolvable"],
        "width_sigma_pc": three_d.get("width_sigma_pc"),
        "depth_sigma_pc": three_d.get("depth_sigma_pc"),
        "expected_parallax_spread_mas": three_d.get("expected_parallax_spread_mas"),
        "parallax_noise_mas": three_d.get("parallax_noise_mas"),
        "depth_note": three_d.get("reason"),
        "centre_density_per_deg2": float(profile["density"][0]),
        "contamination_density_per_deg2": profile["contamination_density"],
    }

    # HR diagram: members at their best-estimate distance; field stars at their
    # own. Photometry is used here for the first time.
    df["abs_g"] = np.where(members, gp.absolute_magnitude(df["phot_g_mean_mag"], df["parallax_post"]),
                           gp.absolute_magnitude(df["phot_g_mean_mag"], df["parallax_corr"]))
    good = df["clean_phot"] & df["ruwe_ok"] & (df["phot_g_mean_mag"] > 6)
    width_sel = good & df["bp_rp"].between(*WIDTH_BP_RP_RANGE)
    cmem = members & width_sel.to_numpy()
    ridge_c, ridge_m = ridge_line(df.loc[cmem, "bp_rp"].to_numpy(), df.loc[cmem, "abs_g"].to_numpy(),
                                  *WIDTH_BP_RP_RANGE)
    w_cluster = sequence_width(df.loc[cmem, "bp_rp"].to_numpy(), df.loc[cmem, "abs_g"].to_numpy(),
                               ridge_c, ridge_m)
    field_sel = (width_sel & (p3 < 0.1)
                 & (df["parallax_corr"] / df["parallax_error"] > WIDTH_PARALLAX_OVER_ERROR)).to_numpy()
    w_field = sequence_width(df.loc[field_sel, "bp_rp"].to_numpy(), df.loc[field_sel, "abs_g"].to_numpy(),
                             ridge_c, ridge_m)
    # How much of each width is just parallax noise?
    mag_err_field = 2.1715 * (df.loc[field_sel, "parallax_error"] / df.loc[field_sel, "parallax_corr"])
    hr_diagram = {
        "n_members_clean": int((members & good.to_numpy()).sum()),
        "width_members_mag": w_cluster["sigma"],
        "width_field_mag": w_field["sigma"],
        "width_ratio": float(w_field["sigma"] / w_cluster["sigma"]),
        "n_width_members": w_cluster["n"],
        "n_width_field": w_field["n"],
        "field_offset_mag": w_field["offset"],
        "field_parallax_noise_mag": float(np.median(mag_err_field)),
        "colour_range": list(WIDTH_BP_RP_RANGE),
    }
    say(f"  main sequence width: members {w_cluster['sigma']:.3f} mag vs field "
        f"{w_field['sigma']:.3f} mag ({hr_diagram['width_ratio']:.1f}x thinner)")

    famous = _famous_status(cluster, df, named, literature)

    result = ClusterResult(
        cluster=cluster.name, key=cluster.key, difficulty=cluster.difficulty,
        n_stars=len(df),
        knot={k: v for k, v in knot.items() if not k.startswith("map_")},
        perspective=perspective, kinematic_only=kinematic_only, membership=membership,
        motion=motion, distance=distance, structure=structure, hr_diagram=hr_diagram,
        famous_stars=famous,
    )
    lit = (literature or {}).get("clusters", {}).get(cluster.key, {}) or {}
    if published_members is not None:
        rt = (lit.get("hr23") or {}).get("rt")
        result.membership_vs_published = compare_members(df, published_members, cluster,
                                                         centre, rt)
    if lit:
        result.comparison = compare(result, lit)

    artefacts = {
        "df": df, "knot": knot, "mix2": mix2, "mix3": mix3, "profile": profile,
        "perspective_cells": cells["cells"],
        "ridge": (ridge_c, ridge_m), "width_members": w_cluster, "width_field": w_field,
        "named": named, "published_members": published_members,
        "field_hrd": (good & (p3 < 0.1)
                      & (df["parallax_corr"] / df["parallax_error"] > FIELD_PARALLAX_OVER_ERROR)
                      ).to_numpy(),
        "clean_members": members & good.to_numpy(),
    }
    return result, artefacts


def _famous_status(cluster, df, named, literature) -> dict:
    """What happened to each named star, and why."""
    if named is None or not len(named):
        return {}
    out = {}
    by_id = df.set_index("source_id")
    for _, row in named.iterrows():
        sid = int(row["source_id"])
        entry = {"source_id": sid, "g": float(row["phot_g_mean_mag"]),
                 "bp_rp": float(row["bp_rp"]) if np.isfinite(row["bp_rp"]) else None}
        if sid in by_id.index:
            s = by_id.loc[sid]
            entry.update(in_sample=True, p_member=float(s["p_member"]),
                         p_motion=float(s["p_motion"]), ruwe=float(s["ruwe"]),
                         member=bool(s["member"]))
        else:
            entry.update(in_sample=False, member=False,
                         reason="Gaia DR3 has a position but no parallax or proper motion "
                                "(2-parameter solution): too bright to measure")
        out[str(row["name"])] = entry
    return out


# --------------------------------------------------------------------------
# Step 7 -- grading
# --------------------------------------------------------------------------

def _pair(ours, ours_err, published, source, published_err=None):
    diff = ours - published
    d = {"ours": ours, "ours_err": ours_err, "published": published,
         "published_err": published_err, "source": source, "difference": diff}
    if published:
        d["percent_diff"] = 100.0 * diff / abs(published)
    return d


def compare(result: ClusterResult, lit: dict) -> dict:
    """Our measurements beside the literature, labelled with where each came from."""
    out = {}
    hr, gc = lit.get("hr23") or {}, lit.get("gc2018") or {}
    m, dist, st = result.motion, result.distance, result.structure
    if hr:
        src = "Hunt & Reffert 2023 (DR3)"
        # Like for like: their means are plain averages of members' motions.
        out["pmra"] = _pair(m["pmra_members_mean"], m["pmra_err"], hr["pmRA"], src,
                            hr.get("e_pmRA"))
        out["pmdec"] = _pair(m["pmdec_members_mean"], m["pmdec_err"], hr["pmDE"], src,
                             hr.get("e_pmDE"))
        # HR23's Plx is a plain mean of members' catalogue parallaxes, so it is
        # compared with our *uncorrected* mean; their distance includes the
        # zero-point correction, as ours does.
        out["parallax_uncorrected"] = _pair(dist["parallax_uncorrected_mas"],
                                            dist["parallax_err_stat"], hr["Plx"], src,
                                            hr.get("e_Plx"))
        out["distance_pc"] = _pair(dist["distance_pc"], dist["distance_err_pc"],
                                   hr["dist50"], src,
                                   0.5 * (hr["dist84"] - hr["dist16"]))
        if hr.get("RV") is not None:
            out["rv_kms"] = _pair(m["rv_kms"], m["rv_err_kms"], hr["RV"], src, hr.get("e_RV"))
        # Half-member radius, measured identically on both member lists: the
        # median sky distance from the same centre, for everyone projected
        # inside the published tidal radius. (Their catalogue's own r50 uses
        # only members inside the tidal radius in 3D, which excludes the halo
        # in front of and behind the cluster -- a different quantity.)
        same = result.membership_vs_published.get("r50_same_method")
        if same:
            out["r50_deg"] = _pair(same["ours"], None, same["published"],
                                   src + " member list, measured our way")
    if gc:
        src = "Gaia Collaboration 2018 (DR2)"
        if gc.get("RV") is not None:
            out["rv_kms_gc2018"] = _pair(m["rv_kms"], m["rv_err_kms"], gc["RV"], src,
                                         gc.get("e_RV"))
        out["parallax_gc2018"] = _pair(dist["parallax_mas"], dist["parallax_err_stat"],
                                       gc["plx"], src + ", before DR3 zero-point work",
                                       gc.get("e_plx"))
        if gc.get("U") is not None:
            v_gc = np.array([gc["U"], gc["V"], gc["W"]], dtype=float)
            v_ours = np.array(m["velocity_icrs_kms"])
            out["speed_kms"] = _pair(m["speed_kms"], None, float(np.linalg.norm(v_gc)), src)
            out["velocity_difference_kms"] = {
                "ours": float(np.linalg.norm(v_ours - v_gc)), "source": src,
                "note": "length of the difference between the two 3D velocity vectors"}
            cp_ours, cp_gc = m["convergent_point_deg"], (gc["RAconv"], gc["DEconv"])
            out["convergent_point_offset_deg"] = {
                "ours": float(gp.angular_separation(*cp_ours, *cp_gc)), "source": src,
                "published": list(cp_gc), "ours_point": cp_ours}
    return out


def compare_members(df, published: pd.DataFrame, cluster: Cluster, centre,
                    rt_deg: float | None = None) -> dict:
    """Overlap between our members and Hunt & Reffert's, star by star.

    Compared only where both could have found a star: inside our cone and
    brighter than our magnitude limit.
    """
    pub = published.copy()
    sep = gp.angular_separation(pub["RA_ICRS"], pub["DE_ICRS"], cluster.ra_deg, cluster.dec_deg)
    pub = pub[(sep < cluster.radius_deg) & (pub["Gmag"] < cluster.g_max)]
    ours = set(df.loc[df["member"], "source_id"].astype(np.int64))
    in_sample = set(df["source_id"].astype(np.int64))
    theirs = set(pub["source_id"].astype(np.int64))
    both = ours & theirs
    missed = theirs - ours
    extra = ours - theirs
    p_missed = df.set_index("source_id").reindex(sorted(missed & in_sample))["p_member"]
    out = {
        "published_in_footprint": len(theirs),
        "ours": len(ours),
        "in_both": len(both),
        "recovered_fraction": len(both) / len(theirs) if theirs else float("nan"),
        "ours_confirmed_fraction": len(both) / len(ours) if ours else float("nan"),
        "missed": len(missed),
        "missed_not_in_our_sample": len(missed - in_sample),
        "missed_median_p": float(p_missed.median()) if len(p_missed) else None,
        "extra": len(extra),
    }
    ext = df[df["source_id"].isin(extra)]
    if len(ext):
        out["extra_median_g"] = float(ext["phot_g_mean_mag"].median())
        out["extra_median_sep_deg"] = float(ext["sep_deg"].median())
        out["extra_median_p"] = float(ext["p_member"].median())
    both_df = df[df["source_id"].isin(both)]
    if len(both_df):
        out["in_both_median_sep_deg"] = float(both_df["sep_deg"].median())
    if rt_deg:
        sep_pub = gp.angular_separation(pub["RA_ICRS"].to_numpy(), pub["DE_ICRS"].to_numpy(),
                                        *centre)
        mem = df[df["member"]]
        sep_ours = gp.angular_separation(mem["ra"].to_numpy(), mem["dec"].to_numpy(), *centre)
        out["r50_same_method"] = {"ours": float(np.median(sep_ours[sep_ours < rt_deg])),
                                  "published": float(np.median(sep_pub[sep_pub < rt_deg])),
                                  "rt_deg": float(rt_deg)}
    return out


def save_results(results: list[ClusterResult]) -> None:
    payload = {
        "_generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "results": [asdict(r) for r in results],
    }
    RESULTS_PATH.write_text(json.dumps(payload, indent=2, default=_json_default))


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)
