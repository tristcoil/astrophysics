"""The analysis: from catalogue rows to the galaxy life cycle.

Run in this order by :func:`analyse`:

1. :func:`survey_geometry` -- footprint, area, spectroscopic completeness and
   fibre-collision weights, all from the list of Main Galaxy Sample targets.
2. :func:`local_sample` -- 0.02 < z < 0.085 galaxies with 1/V_max weights.
3. :func:`main_sequence` -- Renzini & Peng's ridge line, the width of the
   star-forming peak, the green valley, and who is quenched.
4. :func:`environment_analysis` -- local overdensity and Peng et al.'s
   mass-quenching / environment-quenching separability test.
5. :func:`agn_analysis` -- BPT classes across the SFR-mass plane.
6. :func:`gswlc_comparison` -- the same galaxies with an independent method.
7. :func:`census` -- the star-formation-rate density and stellar mass density
   of today's universe, our anchor for the cosmic history.
8. :func:`cosmic_history` -- Madau & Dickinson's compilation, refitted and
   integrated.

Figures are drawn from what this module returns, never by re-running it.
"""

from __future__ import annotations

import json

import healpy as hp
import numpy as np
from scipy.optimize import least_squares

import cosmo_physics as cp
import environment as env
import galaxies as gx
import literature as lit
import mangle
import sfh
from config import (
    FETCH_Z_MAX,
    IMF_TO_SALPETER_MASS,
    IMF_TO_SALPETER_SFR,
    KROUPA_TO_CHABRIER_DEX,
    LOCAL_Z_RANGE,
    MASK_NSIDE,
    MGS_R_LIMIT,
    QUENCH_LOGM_EDGES,
    RESULTS_PATH,
    SEED,
    TILING_DIR,
    VMAX_R_BRIGHT,
    VMAX_R_FAINT,
)
from fetch_data import SPEC_FLAG, load_galaxies, load_gswlc, load_targets, parse_md14_tables

#: Luminosity evolution in the r band, mag per unit redshift (Blanton et al.
#: 2003), used to keep the density tracers' luminosity threshold fixed in
#: evolution-corrected terms (as D2 did following Zehavi et al. 2011).
EVOLUTION_Q = 1.62


def say(msg: str) -> None:
    print(msg, flush=True)


def match_ids(a, b) -> np.ndarray:
    """Index into b for each element of a (-1 where absent). b must be unique."""
    order = np.argsort(b)
    pos = np.searchsorted(b, a, sorter=order)
    pos = np.clip(pos, 0, len(b) - 1)
    hit = b[order][pos] == a
    return np.where(hit, order[pos], -1)


# ---------------------------------------------------------------------------
# 1. Survey geometry
# ---------------------------------------------------------------------------


def survey_geometry(g: dict, tgt_raw: dict) -> dict:
    """Footprint, area, completeness, and per-galaxy collision weights.

    The footprint is the survey's own: the tiled sectors of the legacy tiling
    polygons (mangle.py), whose areas sum to the official spectroscopic area.
    A HEALPix map of it (nside 512, 0.013 deg^2 pixels) serves for edge
    corrections. An area inferred from target counts is kept as a cross-check.
    """
    tu = env.unique_targets(tgt_raw)
    poly = mangle.load_polygons(TILING_DIR / "sdss_tiling_polygons.fits")
    sec = mangle.load_sectors(TILING_DIR / "sdss_sectors.fits")
    tiled = np.isin(poly["sector"], sec["sector"][sec["ntiles"] >= 1])
    area = float(poly["str"][tiled].sum() * (180.0 / np.pi) ** 2)

    def inside(ra, dec):
        k = mangle.which_polygon(poly, ra, dec)
        return (k >= 0) & tiled[np.clip(k, 0, None)]

    tu["in_fp"] = inside(tu["ra"], tu["dec"])
    tu["tsr"] = env.target_sampling_rate(tu)
    fine = mangle.tiled_mask({**poly, "use_caps": np.where(tiled, poly["use_caps"], 0),
                              "ncaps": np.where(tiled, poly["ncaps"], 0)}, MASK_NSIDE)
    fp = {"mask": fine, "nside": MASK_NSIDE}

    # Link galaxies to their target entries
    idx = match_ids(g["objID"], tu["photoObjID"])
    g["target_idx"] = idx
    g["tsr"] = np.where(idx >= 0, tu["tsr"][np.clip(idx, 0, None)], np.nan)
    g["in_fp"] = inside(g["ra"], g["dec"])
    in_sample = np.zeros(len(tu["ra"]), bool)
    in_sample[idx[idx >= 0]] = True
    tu["in_sample"] = in_sample

    def window_stats(r_lo, r_hi):
        """Targets in a DR17 magnitude window, counting the unobserved ones.

        Unobserved targets have no DR17 magnitude, so they are shared out in
        proportion: the same fraction of them is assumed to fall in the window
        as of the observed, linked targets.
        """
        fpm = tu["in_fp"]
        has_r = np.isfinite(tu["r"])
        win = fpm & has_r & (tu["r"] >= r_lo) & (tu["r"] < r_hi)
        n_obs_win = int((win & tu["observed"]).sum())
        f_win = n_obs_win / int((fpm & has_r & tu["observed"]).sum())
        n_unobs = int((fpm & ~tu["observed"]).sum())
        n_tot = n_obs_win + n_unobs * f_win
        ow = win & tu["observed"]
        ins = tu["in_sample"][ow]
        flag = tu["spec_flag"][ow]
        # Not galaxies of our sample by nature: stars, quasars, redshifts
        # outside 0.002-0.3. Failed redshifts and the rest are incompleteness.
        legit = np.isin(flag, [SPEC_FLAG["star"], SPEC_FLAG["qso"], SPEC_FLAG["z_low"],
                               SPEC_FLAG["z_high"]]) & ~ins
        badz = (flag == SPEC_FLAG["bad_z"]) & ~ins
        return {"n_targets": float(n_tot), "n_observed": n_obs_win,
                "n_unobserved_share": float(n_unobs * f_win),
                "fibre_completeness": float(n_obs_win / n_tot),
                "frac_observed_not_in_sample": float(1 - ins.mean()),
                "frac_not_sample_galaxies": float(legit.mean()),
                "frac_stars": float(np.mean((flag == SPEC_FLAG["star"]) & ~ins)),
                "frac_qso": float(np.mean((flag == SPEC_FLAG["qso"]) & ~ins)),
                "frac_bad_z": float(badz.mean()),
                "frac_missing_other": float(np.mean(~ins & ~legit & ~badz)),
                "frac_window": float(f_win)}

    # Cross-check: an area from target counts / density in a well-filled
    # interior. It comes out ~6% low, because choosing "well-filled" pixels
    # selects overdense sky (see README).
    tfp = env.footprint(tu["ra"], tu["dec"])
    nside_c = tfp["nside_coarse"]
    pix = hp.ang2pix(nside_c, tu["ra"], tu["dec"], lonlat=True)
    counts = np.bincount(pix[tu["in_fp"]], minlength=hp.nside2npix(nside_c))
    interior = tfp["coarse"].copy()
    for _ in range(2):  # strip two rings of pixels off the edge
        nb = hp.get_all_neighbours(nside_c, np.flatnonzero(interior))
        bad = ~np.all(np.where(nb >= 0, interior[np.clip(nb, 0, None)], False), axis=0)
        interior[np.flatnonzero(interior)[bad]] = False
    pix_area = hp.nside2pixarea(nside_c, degrees=True)
    n_interior = counts[interior].sum()
    density = n_interior / (interior.sum() * pix_area)
    area_from_targets = int(tu["in_fp"].sum()) / density

    geom = {
        "n_target_entries": tu["n_entries"], "n_targets_unique": int(len(tu["ra"])),
        "n_targets_merged": tu["n_merged"], "n_targets_unlinked": tu["n_unlinked"],
        "n_targets_unlinked_observed": tu["n_unlinked_observed"],
        "n_targets_unobserved": int((~tu["observed"]).sum()),
        "fibre_completeness_all": float(tu["observed"][tu["in_fp"]].mean()),
        "area_tiled_polygons_deg2": area,
        "area_mask_deg2": float(fine.sum() * hp.nside2pixarea(MASK_NSIDE, degrees=True)),
        "area_healpix_deg2": tfp["area_deg2"],
        "area_from_target_density_deg2": float(area_from_targets),
        "target_density_mean_deg2": float(tu["in_fp"].sum() / area),
        "frac_targets_outside_footprint": float(1 - tu["in_fp"].mean()),
        "interior_area_deg2": float(interior.sum() * pix_area),
        "target_density_interior_deg2": float(density),
        "window_vmax": window_stats(VMAX_R_BRIGHT, VMAX_R_FAINT),
        "window_full": window_stats(0.0, MGS_R_LIMIT),
        "n_galaxies_matched_to_targets": int((idx >= 0).sum()),
        "n_galaxies": int(len(idx)),
        "frac_galaxies_in_footprint": float(g["in_fp"].mean()),
    }
    geom["area_deg2"] = area
    geom["solid_angle_sr"] = geom["area_deg2"] * (np.pi / 180.0) ** 2
    return {"geom": geom, "targets": tu, "footprint": fp}


def completeness_weights(g: dict, sel: np.ndarray, geom: dict) -> np.ndarray:
    """Per-galaxy weight for unobserved targets: (1/TSR), renormalised.

    1/TSR restores galaxies lost to fibre collisions where they were lost.
    Collisions are not the only loss (some tiles were never observed), so the
    weights are scaled so that their sum over the V_max window equals the
    number of targets there, minus the observed targets that turned out not to
    belong in the sample (stars, quasars, redshifts outside 0.002-0.3), assumed
    to be the same fraction among the unobserved. Failed redshifts count as
    missing galaxies.
    """
    w = np.where(np.isfinite(g["tsr"]) & (g["tsr"] > 0), 1.0 / g["tsr"], 1.0)
    win = (g["r"] >= VMAX_R_BRIGHT) & (g["r"] < VMAX_R_FAINT) & g["in_fp"]
    st = geom["window_vmax"]
    expected = st["n_targets"] * (1.0 - st["frac_not_sample_galaxies"])
    scale = expected / w[win].sum()
    out = np.zeros(len(w))
    out[sel] = w[sel] * scale
    geom["completeness_weight_scale"] = float(scale)
    return out


# ---------------------------------------------------------------------------
# 2. 1/V_max samples
# ---------------------------------------------------------------------------


def vmax_sample(g: dict, z_range, omega_sr: float) -> tuple[np.ndarray, np.ndarray]:
    """Galaxies inside the z range and flux window, with their V_max (Mpc^3)."""
    sel = (g["z"] >= z_range[0]) & (g["z"] < z_range[1]) & g["in_fp"] \
        & (g["r"] >= VMAX_R_BRIGHT) & (g["r"] < VMAX_R_FAINT)
    vmax = np.full(len(g["z"]), np.nan)
    zmin, zmax = cp.vmax_limits(g["M_r"][sel], g["gr_obs"][sel], z_range,
                                VMAX_R_FAINT, VMAX_R_BRIGHT)
    vmax[sel] = cp.vmax(zmin, zmax, omega_sr)
    return sel, vmax


# ---------------------------------------------------------------------------
# 3. Main sequence, green valley, quenched fraction
# ---------------------------------------------------------------------------


def _sfr_floor(logm):
    """Below this the star-forming peak is never sought (see ridge_line)."""
    return 0.76 * logm - 8.9


def main_sequence(g: dict, sel: np.ndarray, w: np.ndarray, *, logm=None, logsfr=None,
                  agn=None, label: str = "MPA-JHU", valley_override=None) -> dict:
    """Ridge line (AGN excluded and not), SF peak width, valley, quenched flags.

    ``valley_override``: use this Delta_MS for the quenched / green-valley
    division instead of the valley found here (for GSWLC, whose distribution
    has no second peak to find a valley between; see README).
    """
    logm = g["logm"] if logm is None else logm
    logsfr = g["logsfr"] if logsfr is None else logsfr
    agn = gx.is_agn(g["bpt"]) if agn is None else agn
    ok = sel & np.isfinite(logm) & np.isfinite(logsfr) & (w > 0)
    res = {"label": label, "n": int(ok.sum())}

    noagn = ok & ~agn
    ms = gx.ridge_line(logm[noagn], logsfr[noagn], w[noagn], sfr_floor=_sfr_floor)
    ms_all = gx.ridge_line(logm[ok], logsfr[ok], w[ok], sfr_floor=_sfr_floor)
    sfr_w = w * 10 ** np.clip(logsfr, -5, 4)
    ms_sfrw = gx.ridge_line(logm[noagn], logsfr[noagn], sfr_w[noagn], sfr_floor=_sfr_floor)
    res["ridge"] = ms
    # The ridge is clean only where the star-forming peak dominates its mass
    # column; above ~10^10.5 it becomes a shoulder on the quenched population.
    good = np.isfinite(ms["peak"]) & (ms["logm"] < 10.5)
    c_lo, cov_lo = np.polyfit(ms["logm"][good], ms["peak"][good], 1, cov=True)
    res["ridge_lowmass"] = {"slope": float(c_lo[0]), "intercept": float(c_lo[1]),
                            "slope_err": float(np.sqrt(cov_lo[0, 0])), "range": (9.0, 10.5)}
    res["ridge_with_agn"] = {k: ms_all[k] for k in ("slope", "intercept", "slope_err",
                                                     "intercept_err")}
    res["ridge_sfr_weighted"] = {k: ms_sfrw[k] for k in ("slope", "intercept")}
    # offset of the number x SFR ridge at log M = 10.5 (R&P: ~0.2 dex)
    res["sfr_weighted_offset_dex"] = float(gx.ms_sfr(10.5, ms_sfrw) - gx.ms_sfr(10.5, ms))

    dms = logsfr - gx.ms_sfr(logm, ms)
    res["dms"] = dms
    # Width of the SF peak at log M = 10.4-10.6 (Renzini & Peng's Fig. 5 slice)
    sl = ok & (logm > 10.4) & (logm < 10.6)
    res["sf_peak_10p5"] = gx.sf_peak_width(dms[sl], w[sl])
    # Width over the whole fitted mass range
    rng = ok & (logm > 9.0) & (logm < 11.0)
    res["sf_peak_all"] = gx.sf_peak_width(dms[rng], w[rng])

    # The valley: overall and per mass bin
    vall = gx.valley(dms[ok & (logm > 9.5) & (logm < 11.5)],
                     w[ok & (logm > 9.5) & (logm < 11.5)])
    res["valley"] = {"found": vall["found"], "position": vall["position"],
                     "depth": vall["depth_vs_sf_peak"], "q_peak": vall["q_peak"],
                     "centres": vall["centres"], "smoothed": vall["smoothed"]}
    lo, hi, step = QUENCH_LOGM_EDGES
    edges = np.arange(lo, hi + step / 2, 0.4)
    per = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = ok & (logm >= a) & (logm < b)
        if m.sum() < 500:
            continue
        v = gx.valley(dms[m], w[m])
        per.append({"logm": 0.5 * (a + b), "found": v["found"], "position": v["position"],
                    "depth": v["depth_vs_sf_peak"], "q_peak": v["q_peak"],
                    "n": int(m.sum())})
    res["valley_by_mass"] = per
    vpos = vall["position"] if valley_override is None else valley_override
    res["valley_used"] = float(vpos)
    res["quenched"] = ok & (dms < vpos)
    # Green valley: a band 0.5 dex wide centred on the valley floor
    res["green_valley"] = ok & (dms >= vpos - 0.25) & (dms < vpos + 0.25)
    res["star_forming"] = ok & (dms >= vpos + 0.25)
    res["ok"] = ok

    # Quenched fraction vs mass
    centres, fq, err = [], [], []
    for a in np.arange(lo, hi - step / 2, step):
        m = ok & (logm >= a) & (logm < a + step)
        if m.sum() < 50:
            continue
        f, e = gx.weighted_fraction(res["quenched"][m], w[m])
        centres.append(a + step / 2)
        fq.append(f)
        err.append(e)
    res["fq_vs_mass"] = {"logm": centres, "fq": fq, "err": err}
    # mass where half are quenched
    fq_a = np.array(fq)
    cross = np.flatnonzero((fq_a[:-1] < 0.5) & (fq_a[1:] >= 0.5))
    res["logm_half_quenched"] = (float(np.interp(0.5, fq_a[cross[0]:cross[0] + 2],
                                                 centres[cross[0]:cross[0] + 2]))
                                 if len(cross) else float("nan"))
    return res


def color_divider(g: dict, sel: np.ndarray, w: np.ndarray) -> dict:
    """Red/blue split in rest-frame g - r vs stellar mass, found from the data.

    For Peng et al.'s colour-based red fraction: in mass bins, the valley of the
    weighted colour distribution; a straight line through the valley floors.
    """
    pts = []
    for a in np.arange(9.2, 11.2, 0.2):
        m = sel & (g["logm"] >= a) & (g["logm"] < a + 0.2)
        if m.sum() < 1000:
            continue
        edges = np.arange(0.2, 1.1, 0.01)
        h, _ = np.histogram(g["gr0"][m], bins=edges, weights=w[m])
        from scipy.ndimage import gaussian_filter1d
        sm = gaussian_filter1d(h.astype(float), 3)
        c = 0.5 * (edges[1:] + edges[:-1])
        # the two highest local maxima; the minimum between them
        peaks = [j for j in range(1, len(sm) - 1) if sm[j] >= sm[j - 1] and sm[j] >= sm[j + 1]]
        if len(peaks) < 2:
            continue
        p2 = sorted(sorted(peaks, key=lambda j: sm[j])[-2:])
        j = p2[0] + np.argmin(sm[p2[0]:p2[1] + 1])
        pts.append((a + 0.1, c[j]))
    pts = np.array(pts)
    coef = np.polyfit(pts[:, 0], pts[:, 1], 1)
    return {"logm": pts[:, 0].tolist(), "valley_gr": pts[:, 1].tolist(),
            "slope": float(coef[0]), "intercept": float(coef[1])}


# ---------------------------------------------------------------------------
# 4. Environment
# ---------------------------------------------------------------------------


def tracer_threshold(g: dict, z_far: float, r_lim: float = VMAX_R_FAINT) -> dict:
    """Faintest evolution-corrected M_r complete at z_far for 99.5% of colours.

    The tracer galaxies define density, so they must be a volume-limited
    sample: every galaxy above a luminosity threshold, everywhere. Following
    D2, the threshold is set so that even a red galaxy at the threshold
    (the 99.5th percentile of observed g - r near z_far) would still be
    brighter than the flux limit at the far edge.
    """
    near = (g["z"] > z_far - 0.005) & (g["z"] < z_far + 0.005)
    gr_red = float(np.percentile(g["gr_obs"][near], 99.5))
    k_red = float(cp.k_correction("r", z_far, gr_red))
    m_lim_at_far = r_lim - float(cp.distance_modulus(z_far)) - k_red
    # evolution-corrected to z = 0: M0 = M + Q z
    return {"M0_limit": m_lim_at_far + EVOLUTION_Q * z_far, "gr_red": gr_red,
            "k_red": k_red, "z_far": z_far}


def environment_analysis(g: dict, geo: dict, ms: dict, w_local: np.ndarray,
                         local: np.ndarray, red_flag: np.ndarray) -> dict:
    """Overdensity around every local galaxy, and the quenching separability test."""
    fp = geo["footprint"]
    geom = geo["geom"]
    z_lo, z_hi = LOCAL_Z_RANGE
    pad = (1 + z_hi) * env.DENSITY_DV_KMS / cp.C_KMS + 0.001
    tz = (z_lo - pad, z_hi + pad)
    thr = tracer_threshold(g, tz[1])
    M0 = g["M_r"] + EVOLUTION_Q * g["z"]
    trc = (g["z"] >= tz[0]) & (g["z"] < tz[1]) & g["in_fp"] & (M0 < thr["M0_limit"]) \
        & (g["r"] < VMAX_R_FAINT)
    w_trc_all = np.where(np.isfinite(g["tsr"]) & (g["tsr"] > 0), 1.0 / g["tsr"], 1.0)
    xyz = cp.cartesian(g["ra"], g["dec"], g["d_c"])
    # mean density of tracers, and its flatness with redshift
    vol = cp.vmax(tz[0], tz[1], geom["solid_angle_sr"])
    n_trc_w = float(w_trc_all[trc].sum())
    zb = np.linspace(tz[0], tz[1], 8)
    nz = []
    for a, b in zip(zb[:-1], zb[1:]):
        m = trc & (g["z"] >= a) & (g["z"] < b)
        nz.append(float(w_trc_all[m].sum() / cp.vmax(a, b, geom["solid_angle_sr"])))
    gal = local & np.isfinite(g["logm"])
    dens = env.nth_neighbour_density(xyz[gal], g["z"][gal], xyz[trc], w_trc_all[trc])
    sigma_mean = env.mean_tracer_surface_density(n_trc_w, vol, g["z"][gal])
    one_plus_delta = np.full(len(g["z"]), np.nan)
    one_plus_delta[gal] = dens["sigma"] / sigma_mean
    edge_ok, cov = env.edge_ok(fp, g["ra"][gal], g["dec"][gal], dens["d_n"], g["z"][gal])
    usable = np.zeros(len(g["z"]), bool)
    usable[gal] = edge_ok & ~dens["lower_bound"]
    res = {
        "tracer": {**thr, "n": int(trc.sum()), "n_weighted": n_trc_w,
                   "nbar_mpc3": n_trc_w / vol, "z_range": tz,
                   "n_of_z_mpc3": nz, "n_of_z_edges": zb.tolist()},
        "one_plus_delta": one_plus_delta, "usable": usable,
        "n_local": int(gal.sum()), "n_usable": int(usable.sum()),
        "frac_edge_dropped": float(1 - edge_ok.mean()),
        "median_d5_mpc": float(np.median(dens["d_n"])),
    }
    res["quench"] = separability(g["logm"], one_plus_delta, ms["quenched"], w_local,
                                 usable & ms["ok"])
    res["red"] = separability(g["logm"], one_plus_delta, red_flag, w_local,
                              usable & ms["ok"])
    return res


def separability(logm, opd, flag, w, sel, *, mass_edges=None, dens_edges=None) -> dict:
    """Peng et al.'s test: is f(m, rho) = eps_m + eps_rho - eps_m eps_rho?

    1. Weighted quenched (or red) fraction on a grid of mass x log(1+delta).
    2. eps_rho(rho; m) relative to the sparsest density bin, and eps_m(m; rho)
       relative to the lowest mass bin (their eqs. 3-4). Separability means
       eps_rho does not depend on m and eps_m does not depend on rho.
    3. Fit the separable model f = 1 - exp(-(rho/p1)^p2 - (m/p3)^p4) to the
       grid (Peng et al.'s eq. 7) and compare the fit with a free grid.
    """
    mass_edges = np.arange(9.0, 11.41, 0.3) if mass_edges is None else mass_edges
    ld = np.log10(opd)
    dens_edges = np.arange(-0.6, 1.81, 0.3) if dens_edges is None else dens_edges
    nm, nd = len(mass_edges) - 1, len(dens_edges) - 1
    f = np.full((nm, nd), np.nan)
    e = np.full((nm, nd), np.nan)
    n = np.zeros((nm, nd), int)
    for i in range(nm):
        for j in range(nd):
            m = sel & (logm >= mass_edges[i]) & (logm < mass_edges[i + 1]) \
                & (ld >= dens_edges[j]) & (ld < dens_edges[j + 1])
            n[i, j] = int(m.sum())
            if n[i, j] >= 30:
                f[i, j], e[i, j] = gx.weighted_fraction(flag[m], w[m])
    mc = 0.5 * (mass_edges[1:] + mass_edges[:-1])
    dc = 0.5 * (dens_edges[1:] + dens_edges[:-1])

    # relative efficiencies
    eps_rho = (f - f[:, :1]) / (1 - f[:, :1])
    eps_m = (f - f[:1, :]) / (1 - f[:1, :])

    # separable fit
    good = np.isfinite(f) & (e > 0)
    M, D = np.meshgrid(mc, dc, indexing="ij")

    def model(p, logm_, logd_):
        lp1, p2, lp3, p4 = p
        return 1 - np.exp(-(10 ** (logd_ - lp1)) ** p2 - (10 ** (logm_ - lp3)) ** p4)

    def resid(p):
        return ((f - model(p, M, D)) / e)[good]

    fit = least_squares(resid, [1.84, 0.6, 10.56, 0.8],
                        bounds=([-1, 0.05, 9.0, 0.05], [4, 5, 12.5, 5]))
    chi2 = float(np.sum(fit.fun**2))
    dof = int(good.sum() - 4)
    jac = fit.jac
    try:
        cov = np.linalg.inv(jac.T @ jac) * max(chi2 / dof, 1.0)
        perr = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        perr = np.full(4, np.nan)
    # residuals of the separable model, as a fraction
    resid_f = np.where(good, f - model(fit.x, M, D), np.nan)
    # Non-separable alternative: the same model but with p1 a free linear
    # function of mass (environment effect allowed to depend on mass)
    def model2(p, logm_, logd_):
        lp1, s1, p2, lp3, p4 = p
        lp1m = lp1 + s1 * (logm_ - 10.0)
        return 1 - np.exp(-(10 ** (logd_ - lp1m)) ** p2 - (10 ** (logm_ - lp3)) ** p4)

    fit2 = least_squares(lambda p: ((f - model2(p, M, D)) / e)[good],
                         [fit.x[0], 0.0, fit.x[1], fit.x[2], fit.x[3]],
                         bounds=([-1, -3, 0.05, 9.0, 0.05], [4, 3, 5, 12.5, 5]))
    cov2 = np.linalg.inv(fit2.jac.T @ fit2.jac) * max(float(np.sum(fit2.fun**2)) / (dof - 1), 1)
    return {
        "mass_centres": mc.tolist(), "dens_centres": dc.tolist(),
        "mass_edges": mass_edges.tolist(), "dens_edges": dens_edges.tolist(),
        "f": f, "err": e, "n": n, "eps_rho": eps_rho, "eps_m": eps_m,
        "fit": {"log_p1": float(fit.x[0]), "p2": float(fit.x[1]),
                "log_p3": float(fit.x[2]), "p4": float(fit.x[3]),
                "err": perr.tolist(), "chi2": chi2, "dof": dof},
        "resid": resid_f,
        "rms_resid": float(np.sqrt(np.nanmean(resid_f**2))),
        "max_abs_resid": float(np.nanmax(np.abs(resid_f))),
        "mass_dependent_env": {"slope_logp1_per_dex": float(fit2.x[1]),
                               "slope_err": float(np.sqrt(cov2[1, 1])),
                               "chi2": float(np.sum(fit2.fun**2))},
    }


# ---------------------------------------------------------------------------
# 5. AGN across the plane
# ---------------------------------------------------------------------------


def agn_analysis(g: dict, ms: dict, w: np.ndarray, *, dms=None, label="MPA-JHU") -> dict:
    """Fraction of galaxies with AGN-powered lines, by position on the plane."""
    ok = ms["ok"]
    dms = ms["dms"] if dms is None else dms
    agn = gx.is_agn(g["bpt"])
    out = {"label": label}
    # by Delta_MS in mass bins
    rows = []
    for a, b in ((9.5, 10.0), (10.0, 10.5), (10.5, 11.0), (11.0, 11.5)):
        m = ok & (g["logm"] >= a) & (g["logm"] < b) & np.isfinite(dms)
        edges = np.arange(-2.6, 1.01, 0.2)
        cen, fr, er = [], [], []
        for lo_, hi_ in zip(edges[:-1], edges[1:]):
            mm = m & (dms >= lo_) & (dms < hi_)
            if mm.sum() < 40:
                continue
            f, e = gx.weighted_fraction(agn[mm], w[mm])
            cen.append(0.5 * (lo_ + hi_))
            fr.append(f)
            er.append(e)
        if fr:
            rows.append({"logm": (a, b), "dms": cen, "f_agn": fr, "err": er,
                         "dms_at_peak": cen[int(np.argmax(fr))], "peak": max(fr)})
    out["by_dms"] = rows
    # in the three zones
    zones = {}
    for name, flag in (("star_forming", ms["star_forming"]), ("green_valley", ms["green_valley"]),
                       ("quenched", ms["quenched"])):
        m = flag & (g["logm"] > 10.0) & (g["logm"] < 11.0)
        zones[name] = gx.weighted_fraction(agn[m], w[m])
    out["zones_logm10_11"] = zones
    # class mix overall
    mix = {}
    for c, name in gx.BPT_NAMES.items():
        mix[name] = float(np.sum(w[ok] * (g["bpt"][ok] == c)) / np.sum(w[ok]))
    out["class_mix_weighted"] = mix
    return out


def bpt_agreement(g: dict, sel: np.ndarray) -> dict:
    """Our BPT classes against MPA-JHU's own bptclass column."""
    ours = g["bpt"][sel]
    theirs = g["bptclass_mpa"][sel]
    # MPA-JHU: 1 SF, 2 low-S/N SF, 3 composite, 4 AGN, 5 low-S/N LINER, -1 none
    table = {}
    for c, name in gx.BPT_NAMES.items():
        m = ours == c
        if m.sum() == 0:
            continue
        vals, cnt = np.unique(theirs[m], return_counts=True)
        table[name] = {int(v): int(k) for v, k in zip(vals, cnt)}
    full = (ours == gx.BPT_SF) | (ours == gx.BPT_COMPOSITE) | (ours == gx.BPT_SEYFERT) \
        | (ours == gx.BPT_LINER)
    map_ours = {gx.BPT_SF: 1, gx.BPT_COMPOSITE: 3}
    agree_sf = float(np.mean(theirs[ours == gx.BPT_SF] == 1))
    agree_comp = float(np.mean(theirs[ours == gx.BPT_COMPOSITE] == 3))
    agree_agn = float(np.mean(theirs[(ours == gx.BPT_SEYFERT) | (ours == gx.BPT_LINER)] == 4))
    del map_ours
    return {"table": table, "n_full": int(full.sum()), "agree_sf": agree_sf,
            "agree_composite": agree_comp, "agree_agn": agree_agn}


# ---------------------------------------------------------------------------
# 6. GSWLC
# ---------------------------------------------------------------------------


def attach_gswlc(g: dict, name: str = "GSWLC-X2.dat.gz") -> dict:
    """GSWLC-2 masses and SFRs matched by plate-MJD-fibre."""
    gs = load_gswlc(name)
    key_g = (g["plate"].astype(np.int64) * 100000 + (g["mjd"].astype(np.int64) - 50000)) \
        * 1000 + g["fiberID"].astype(np.int64)
    key_s = (gs["plate"].to_numpy(np.int64) * 100000 + (gs["mjd"].to_numpy(np.int64) - 50000)) \
        * 1000 + gs["fiberid"].to_numpy(np.int64)
    uniq, first = np.unique(key_s, return_index=True)
    idx = match_ids(key_g, uniq)
    has = idx >= 0
    j = first[np.clip(idx, 0, None)]
    def col(c):
        v = gs[c].to_numpy(float)[j]
        return np.where(has, v, np.nan)
    flag_sed = np.where(has, gs["flag_sed"].to_numpy()[j], -1)
    logm = col("logmstar")
    logsfr = col("logsfr")
    good = has & (flag_sed == 0) & (logm > 0) & (logsfr > -90)
    return {"logm": np.where(good, logm, np.nan), "logsfr": np.where(good, logsfr, np.nan),
            "has": has, "good": good, "uv_survey": col("uv_survey"),
            "flag_midir": col("flag_midir"), "flag_mgs": col("flag_mgs"),
            "file": name, "n_rows": int(len(gs)), "n_matched": int(has.sum()),
            "n_good": int(good.sum())}


def gswlc_comparison(g: dict, gs: dict, sel: np.ndarray) -> dict:
    """MPA-JHU vs GSWLC: masses, and sSFRs normalised by the same mass.

    As in Salim et al. (2016): MPA-JHU values moved from Kroupa to Chabrier
    (-0.025 dex) and both sSFRs divided by the GSWLC mass, so the comparison
    isolates the SFRs.
    """
    m = sel & gs["good"] & g["sfr_ok"]
    dm = (g["logm"] + KROUPA_TO_CHABRIER_DEX - gs["logm"])[m]
    ssfr_g = (gs["logsfr"] - gs["logm"])[m]
    ssfr_mpa = (g["logsfr"] + KROUPA_TO_CHABRIER_DEX - gs["logm"])[m]
    d = ssfr_mpa - ssfr_g
    bpt = g["bpt"][m]
    out = {"n": int(m.sum()), "dlogm_median": float(np.median(dm)),
           "dlogm_std": float(np.std(dm)),
           "dlogm_robust_sigma": float(1.4826 * np.median(np.abs(dm - np.median(dm))))}
    groups = {"star_forming": bpt == gx.BPT_SF, "agn": gx.is_agn(bpt),
              "weak_or_none": (bpt == gx.BPT_WEAK) | (bpt == gx.BPT_NONE)}
    for name, gm in groups.items():
        dd = d[gm]
        out[f"dssfr_{name}"] = {"n": int(gm.sum()), "median": float(np.median(dd)),
                                "std": float(np.std(dd)),
                                "robust_sigma": float(1.4826 * np.median(np.abs(dd - np.median(dd))))}
    # mass offset vs GSWLC sSFR (Salim: 0.03-0.13 dex, sSFR dependent)
    edges = np.arange(-13.0, -8.4, 0.5)
    cen, med = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        mm = (ssfr_g >= a) & (ssfr_g < b)
        if mm.sum() > 200:
            cen.append(0.5 * (a + b))
            med.append(float(np.median(dm[mm])))
    out["dlogm_vs_ssfr"] = {"ssfr": cen, "median": med}
    out["arrays"] = {"ssfr_gswlc": ssfr_g, "ssfr_mpa": ssfr_mpa, "bpt": bpt, "dlogm": dm}
    return out


# ---------------------------------------------------------------------------
# 7. Census of today's universe
# ---------------------------------------------------------------------------


def jackknife_regions(ra, dec, nside: int = 4) -> np.ndarray:
    return hp.ang2pix(nside, ra, dec, lonlat=True)


def census(g: dict, geo: dict, *, z_range=(0.02, 0.2), logm=None, logsfr=None,
           valid=None, imf: str = "kroupa", label: str = "MPA-JHU", sfr_linear=None,
           zone_flags=None) -> dict:
    """Star-formation-rate density and stellar mass density by 1/V_max.

    rho = sum_i w_i X_i / V_max,i, with w_i the completeness weight. Galaxies
    without a valid value are represented by those with one: the sum is
    divided by the weighted fraction of valid galaxies. The SFR density is
    quoted at the sample's SFR-weighted mean redshift, and also evolved to
    z = 0.1 with SFR ∝ (1+z)^beta for beta = 0 and 3, as Brinchmann et al.
    (2004) did.
    """
    geom = geo["geom"]
    omega = geom["solid_angle_sr"]
    sel, vmax = vmax_sample(g, z_range, omega)
    w = completeness_weights(g, sel, geom)
    logm = g["logm"] if logm is None else logm
    logsfr = g["logsfr"] if logsfr is None else logsfr
    valid = (g["sfr_ok"] if valid is None else valid) & sel
    frac_valid = float(np.sum(w[valid] / vmax[valid]) / np.sum(w[sel] / vmax[sel]))
    wv = np.where(valid, w / vmax, 0.0)
    sfr = np.where(valid, 10 ** np.clip(logsfr, -9, 5), 0.0)
    mass = np.where(valid, 10 ** np.clip(logm, 0, 13), 0.0)
    if sfr_linear is not None:  # e.g. the likelihood mean instead of the median
        sfr = np.where(valid, sfr_linear, 0.0)
    out = {"label": label, "imf": imf, "z_range": z_range, "n": int(valid.sum()),
           "frac_valid": frac_valid}
    rho_sfr = np.sum(wv * sfr) / frac_valid
    rho_m = np.sum(wv * mass) / frac_valid
    z_eff = float(np.sum(wv * sfr * g["z"]) / np.sum(wv * sfr))
    out["rho_sfr"] = float(rho_sfr)
    out["rho_star"] = float(rho_m)
    out["z_eff_sfr"] = z_eff
    out["z_eff_mass"] = float(np.sum(wv * mass * g["z"]) / np.sum(wv * mass))
    for beta in (0, 3):
        ev = ((1.1) / (1 + g["z"])) ** beta
        out[f"rho_sfr_z0.1_beta{beta}"] = float(np.sum(wv * sfr * ev) / frac_valid)
    out["rho_sfr_z0.1"] = 0.5 * (out["rho_sfr_z0.1_beta0"] + out["rho_sfr_z0.1_beta3"])
    out["rho_sfr_salpeter"] = out["rho_sfr_z0.1"] * IMF_TO_SALPETER_SFR[imf]
    out["rho_star_salpeter"] = out["rho_star"] * IMF_TO_SALPETER_MASS[imf]

    # jackknife over sky regions
    reg = jackknife_regions(g["ra"], g["dec"])
    tu = geo["targets"]
    treg = jackknife_regions(tu["ra"], tu["dec"])
    tsel = tu["in_fp"]
    n_t = np.bincount(treg[tsel], minlength=hp.nside2npix(4))
    regions = np.flatnonzero(n_t > 0.02 * n_t.max())
    js, jm = [], []
    for k in regions:
        keep = reg != k
        scale = n_t.sum() / (n_t.sum() - n_t[k])  # area lost -> V_max shrinks
        fv = np.sum(np.where(keep & valid, w / vmax, 0)) / np.sum(np.where(keep & sel, w / vmax, 0))
        ev = ((1.1) / (1 + g["z"]))
        s0 = np.sum(np.where(keep, wv * sfr, 0))
        s3 = np.sum(np.where(keep, wv * sfr * ev**3, 0))
        js.append(0.5 * (s0 + s3) * scale / fv)
        jm.append(np.sum(np.where(keep, wv * mass, 0)) * scale / fv)
    js, jm = np.array(js), np.array(jm)
    nreg = len(regions)
    out["rho_sfr_err"] = float(np.sqrt((nreg - 1) / nreg * np.sum((js - js.mean()) ** 2)))
    out["rho_star_err"] = float(np.sqrt((nreg - 1) / nreg * np.sum((jm - jm.mean()) ** 2)))
    out["n_jackknife"] = int(nreg)

    # who forms the stars: by mass
    edges = np.arange(8.0, 12.21, 0.2)
    cen, sfr_share, mass_share = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        mm = valid & (logm >= a) & (logm < b)
        cen.append(0.5 * (a + b))
        sfr_share.append(float(np.sum(wv[mm] * sfr[mm]) / frac_valid))
        mass_share.append(float(np.sum(wv[mm] * mass[mm]) / frac_valid))
    out["by_mass"] = {"logm": cen, "rho_sfr": sfr_share, "rho_star": mass_share}
    tot = np.sum(sfr_share)
    cum = np.cumsum(sfr_share) / tot
    out["logm_half_sfr"] = float(np.interp(0.5, cum, np.array(edges[1:])))
    cumm = np.cumsum(mass_share) / np.sum(mass_share)
    out["logm_half_mass"] = float(np.interp(0.5, cumm, np.array(edges[1:])))
    if zone_flags is not None:
        zs = {}
        for name, flag in zone_flags.items():
            mm = valid & flag
            zs[name] = {"sfr_share": float(np.sum(wv[mm] * sfr[mm]) / np.sum(wv * sfr)),
                        "mass_share": float(np.sum(wv[mm] * mass[mm]) / np.sum(wv * mass))}
        out["zones"] = zs
    return out


# ---------------------------------------------------------------------------
# 8. Cosmic history
# ---------------------------------------------------------------------------


def cosmic_history(n_boot: int = 1000) -> dict:
    sfrd, smd = parse_md14_tables()
    pts = sfh.table1_points(sfrd)
    pub = dict(lit.MD14_FIT)
    ours = sfh.fit_md(pts)
    resid_pub = sfh._residuals([np.log10(pub["A"]), pub["B"], pub["C"], pub["D"]], pts["z"],
                               pts["logpsi"], pts["err_plus"], pts["err_minus"])
    boot = sfh.bootstrap_fit(pts, n=n_boot)
    boot_summ = [sfh.summarise(dict(zip("ABCD", b))) for b in boot]
    keys = ["z_peak", "t_peak_gyr", "rho_star_z0", "z_half_mass", "frac_before_peak",
            "z_last_quarter", "frac_before_reionization_end", "mean_stellar_age_gyr",
            "efold_curve_gyr", "psi_z0"]
    boot_stats = {k: {"p16": float(np.percentile([b[k] for b in boot_summ], 16)),
                      "p50": float(np.percentile([b[k] for b in boot_summ], 50)),
                      "p84": float(np.percentile([b[k] for b in boot_summ], 84))}
                  for k in keys}
    sig = 0.5 * (pts["err_plus"] + pts["err_minus"])
    # stellar mass density vs Table 2
    h_pub = sfh.history(pub)
    h_ours = sfh.history(ours)
    smd_z = 0.5 * (smd["z_lo"].to_numpy() + smd["z_hi"].to_numpy())
    pred_pub = np.interp(smd_z, h_pub["z"], h_pub["rho_star"])
    smd_resid = smd["log_value"].to_numpy() - np.log10(pred_pub)
    zbins = [(0, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 3.0), (3.0, 8.5)]
    smd_by_z = []
    for a, b in zbins:
        m = (smd_z >= a) & (smd_z < b)
        smd_by_z.append({"z": (a, b), "n": int(m.sum()),
                         "median_resid_dex": float(np.median(smd_resid[m]))})
    omega_b_h2_d3 = 0.02221  # project D3's Planck TT fit
    summ_pub = sfh.summarise(pub)
    return {
        "points": pts, "smd": smd, "smd_z": smd_z, "smd_resid": smd_resid,
        "smd_by_z": smd_by_z,
        "published": {**pub, "chi2": float(np.sum(resid_pub**2)), "dof": len(pts["z"]) - 4,
                      "summary": summ_pub,
                      "summary_planck": sfh.summarise(pub, sfh.COSMO_PLANCK)},
        "ours": {**ours, "summary": sfh.summarise(ours), "boot": boot_stats,
                 "boot_params": boot},
        "efold_data_weighted": sfh.efold_time_gyr(pts["z"], pts["logpsi"], sig),
        "efold_data_unweighted": sfh.efold_time_gyr(pts["z"], pts["logpsi"]),
        "high_z_slope": sfh.high_z_slope(pts["z"], pts["logpsi"]),
        "baryon_fraction": {
            "md14_omega_b": summ_pub["rho_star_z0"] / sfh.baryon_density(0.045 * 0.49),
            "d3_planck": summ_pub["rho_star_z0"] / sfh.baryon_density(omega_b_h2_d3),
            "d3_planck_chabrier": summ_pub["rho_star_z0"] * 0.61
            / sfh.baryon_density(omega_b_h2_d3),
            "omega_b_h2_d3": omega_b_h2_d3,
        },
        "history_pub": h_pub, "history_ours": h_ours,
    }


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def analyse(*, quiet: bool = False) -> dict:
    log = (lambda m: None) if quiet else say
    log("Loading catalogues")
    gal, gmeta = load_galaxies()
    tgt, tmeta = load_targets()
    g = gx.derive(gal)
    res = {"meta": {"galaxies": gmeta, "targets": tmeta}, "n_galaxies": int(len(g["z"]))}
    res["mpa_coverage"] = {"has_mpa": float(g["has_mpa"].mean()),
                           "mass_ok": float(g["mass_ok"].mean()),
                           "sfr_ok": float(g["sfr_ok"].mean()),
                           "reliable": float(np.mean(g["reliable"] == 1)),
                           "n_ssfr_above_minus8": int(np.sum(
                               g["mass_ok"] & (g["logsfr"] > -9) & (g["logsfr"] < 4)
                               & (g["logsfr"] - g["logm"] >= -8.0)))}

    log("Survey geometry: footprint, completeness, fibre collisions")
    geo = survey_geometry(g, tgt)
    res["geometry"] = geo["geom"]

    log("Local sample, 0.02 < z < 0.085, with 1/V_max weights")
    local, vmax = vmax_sample(g, LOCAL_Z_RANGE, geo["geom"]["solid_angle_sr"])
    local &= g["sfr_ok"] & (g["reliable"] == 1)
    wc = completeness_weights(g, local, geo["geom"])
    w_local = np.where(local, wc / vmax, 0.0)
    res["local"] = {"n": int(local.sum()), "z_range": LOCAL_Z_RANGE}

    log("Main sequence (MPA-JHU)")
    ms = main_sequence(g, local, w_local)
    res["ms"] = ms
    log("Colour divider")
    cd = color_divider(g, local, w_local)
    red = g["gr0"] > cd["slope"] * g["logm"] + cd["intercept"]
    res["color_divider"] = cd
    res["red_fraction_local"] = gx.weighted_fraction(red[local], w_local[local])
    res["quenched_fraction_local"] = gx.weighted_fraction(ms["quenched"][local],
                                                          w_local[local])
    # how often colour and sSFR disagree
    m = local & ms["ok"]
    res["red_but_starforming"] = gx.weighted_fraction((red & ~ms["quenched"])[m]
                                                      , w_local[m])
    res["blue_but_quenched"] = gx.weighted_fraction((~red & ms["quenched"])[m], w_local[m])

    log("BPT classes")
    res["bpt_agreement"] = bpt_agreement(g, local)
    res["agn"] = agn_analysis(g, ms, w_local)

    log("Environment: tracers, 5th-neighbour densities, separability")
    res["env"] = environment_analysis(g, geo, ms, w_local, local, red)

    log("GSWLC-2")
    gs = attach_gswlc(g)
    res["gswlc_match"] = {k: gs[k] for k in ("file", "n_rows", "n_matched", "n_good")}
    res["gswlc_match"]["frac_local_matched"] = float(np.mean(gs["good"][local]))
    res["gswlc_cmp"] = gswlc_comparison(g, gs, local)
    gloc = local & gs["good"]
    ms_g = main_sequence(g, gloc, w_local, logm=gs["logm"], logsfr=gs["logsfr"],
                         label="GSWLC-2", valley_override=ms["valley"]["position"])
    res["ms_gswlc"] = ms_g
    dms_g = np.where(gloc, gs["logsfr"] - gx.ms_sfr(gs["logm"], ms_g["ridge"]), np.nan)
    agn_g = agn_analysis({**g, "logm": gs["logm"]}, ms_g, w_local, dms=dms_g,
                         label="GSWLC-2")
    res["agn_gswlc"] = agn_g
    # GSWLC-M2 (deeper UV) for the quenched fraction cross-check
    gsm = attach_gswlc(g, "GSWLC-M2.dat.gz")
    mloc = local & gsm["good"]
    ms_m = main_sequence(g, mloc, w_local, logm=gsm["logm"], logsfr=gsm["logsfr"],
                         label="GSWLC-M2", valley_override=ms["valley"]["position"])
    res["ms_gswlc_m2"] = ms_m

    log("Census at z ~ 0.1")
    zones = {"star_forming": ms["star_forming"], "green_valley": ms["green_valley"],
             "quenched": ms["quenched"]}
    res["census"] = census(g, geo)
    # likelihood mean instead of median (from the five percentiles)
    g["sfr_mean"] = gx.likelihood_mean_sfr(g["logsfr_p2p5"], g["logsfr_lo"], g["logsfr"],
                                           g["logsfr_hi"], g["logsfr_p97p5"])
    res["census_mean"] = census(g, geo, sfr_linear=g["sfr_mean"],
                                label="MPA-JHU (likelihood mean)")
    res["census_mean_local"] = census(g, geo, z_range=LOCAL_Z_RANGE, sfr_linear=g["sfr_mean"],
                                      zone_flags=zones, label="MPA-JHU mean 0.02-0.085")
    res["census_gswlc"] = census(g, geo, logm=gs["logm"], logsfr=gs["logsfr"],
                                 valid=gs["good"], imf="chabrier", label="GSWLC-X2")
    # GSWLC covers 90% of SDSS: rescale by the matched fraction of the window
    sel_win, _ = vmax_sample(g, (0.02, 0.2), geo["geom"]["solid_angle_sr"])
    res["census_gswlc"]["gswlc_coverage"] = float(np.mean(gs["has"][sel_win]))
    res["census_local"] = census(g, geo, z_range=LOCAL_Z_RANGE, zone_flags=zones,
                                 label="MPA-JHU 0.02-0.085")

    log("Cosmic star-formation history")
    res["sfh"] = cosmic_history()

    # keep arrays needed by the figures
    res["_g"] = g
    res["_w_local"] = w_local
    res["_local"] = local
    res["_red"] = red
    res["_gswlc"] = gs
    res["_dms_gswlc"] = dms_g
    res["_geo"] = geo
    return res


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items() if not str(k).startswith("_")}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        if x.size > 400:
            return None
        return _jsonable(x.tolist())
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if hasattr(x, "to_dict"):
        return None
    return x


def save_results(res: dict) -> None:
    RESULTS_PATH.write_text(json.dumps(_jsonable(res), indent=1))
