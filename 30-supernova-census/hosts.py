"""Which galaxies host which supernovae?

Core-collapse supernovae come from stars that live less than ~50 Myr, so they
should appear only where stars are forming. Type Ia supernovae come from white
dwarfs, which can wait billions of years, so they should appear in old,
"quenched" galaxies too. We test this on the star-formation-rate vs stellar-
mass plane of project J3, with the MPA-JHU measurements for SDSS galaxies.

**Host assignment.** For every supernova we take all SDSS galaxies within
30 projected kpc (at most 90") and pick the one with the smallest
*directional light radius* distance: the separation divided by the galaxy's
elliptical radius in the direction of the supernova (Sullivan et al. 2006;
Gupta et al. 2016). A galaxy farther than 4 such radii hosts nothing. If the
host has a spectrum, its redshift must agree with the supernova's.

**What a host distribution should look like.** If supernovae of some kind
occurred at a rate proportional to stellar mass, their hosts would be drawn
from the galaxy population weighted by mass; if proportional to star
formation, weighted by SFR. The "A + B" model (Mannucci et al. 2005;
Scannapieco & Bildsten 2005; Sullivan et al. 2006) gives Type Ia supernovae a
rate A M + B SFR: a delayed component that follows old stars and a prompt one
that follows young stars. Because the survey finds supernovae out to
different distances, each weighted galaxy population is re-weighted to the
hosts' own redshift distribution before comparing.
"""

from __future__ import annotations

import glob

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

import config as C


# ---------------------------------------------------------------------------
# Host assignment
# ---------------------------------------------------------------------------

def load_candidates() -> pd.DataFrame:
    pos = pd.read_csv(C.SDSS_HOST_DIR / "positions.csv")
    parts = [pd.read_csv(f, dtype={"objID": str, "specObjID": str})
             for f in sorted(glob.glob(str(C.SDSS_HOST_DIR / "batch_*.csv")))]
    parts = [p for p in parts if len(p)]
    cand = pd.concat(parts, ignore_index=True)
    cand["ztf_id"] = pos["ZTFID"].to_numpy()[cand["sn_index"].to_numpy()]
    cand["sn_ra"] = pos["ra"].to_numpy()[cand["sn_index"].to_numpy()]
    cand["sn_dec"] = pos["dec"].to_numpy()[cand["sn_index"].to_numpy()]
    return cand


def dlr_distance(cand: pd.DataFrame) -> np.ndarray:
    """Separation / (elliptical half-light radius toward the supernova).

    The shape comes from SDSS's r-band model fit that describes the galaxy
    better (de Vaucouleurs if fracDeV > 0.5, else exponential): effective
    radius, axis ratio b/a and position angle (degrees east of north)."""
    dev = cand["fracDeV_r"].to_numpy() > 0.5
    r_e = np.where(dev, cand["deVRad_r"], cand["expRad_r"]).astype(float)
    ab = np.where(dev, cand["deVAB_r"], cand["expAB_r"]).astype(float)
    phi = np.radians(np.where(dev, cand["deVPhi_r"], cand["expPhi_r"]).astype(float))
    # direction from the galaxy centre to the supernova, east of north
    dra = (cand["sn_ra"] - cand["ra"]).to_numpy() * np.cos(np.radians(cand["dec"].to_numpy()))
    ddec = (cand["sn_dec"] - cand["dec"]).to_numpy()
    pa_sn = np.arctan2(dra, ddec)
    theta = pa_sn - phi
    a = np.clip(r_e, 0.2, None)          # semi-major axis (arcsec)
    b = a * np.clip(ab, 0.05, 1.0)
    r_dir = a * b / np.sqrt((b * np.cos(theta)) ** 2 + (a * np.sin(theta)) ** 2)
    return cand["sep_arcsec"].to_numpy() / r_dir


def assign_hosts(df: pd.DataFrame) -> pd.DataFrame:
    """One row per supernova in SDSS imaging: the host (if any) and its
    MPA-JHU mass and SFR (if it has a usable spectrum)."""
    cand = load_candidates()
    cand["d_dlr"] = dlr_distance(cand)
    cand = cand[np.isfinite(cand["d_dlr"])]
    best = cand.sort_values("d_dlr").groupby("ztf_id", as_index=False).first()
    best = best[best["d_dlr"] < C.HOST_MAX_DLR]
    z_sn = df.set_index("ztf_id")["redshift"]
    best["z_sn"] = z_sn.reindex(best["ztf_id"]).to_numpy()
    spec_z = pd.to_numeric(best["spec_z"], errors="coerce")
    best["has_spec"] = spec_z.notna()
    best["z_agrees"] = (spec_z - best["z_sn"]).abs() < C.HOST_DZ_MAX
    for col in ("lgm_tot_p16", "lgm_tot_p50", "lgm_tot_p84", "sfr_tot_p16", "sfr_tot_p50",
                "sfr_tot_p84", "specsfr_tot_p50"):
        best[col] = pd.to_numeric(best[col], errors="coerce")
    # MPA-JHU flags failed fits with -9999 or -99
    best["has_mpa"] = (best["has_spec"] & best["z_agrees"] & (best["lgm_tot_p50"] > 5)
                       & (best["sfr_tot_p50"] > -10))
    best["ssfr"] = specific_sfr(best)
    # the SN's offset in units of the host's light radius, and in kpc
    best["host_spec_z"] = spec_z
    keep = ["ztf_id", "objID", "ra", "dec", "sep_arcsec", "d_dlr", "petroMag_r", "extinction_r",
            "modelMag_g", "modelMag_r", "modelMag_i", "fracDeV_r", "petroR50_r", "specObjID",
            "host_spec_z", "z_sn", "has_spec", "z_agrees", "has_mpa", "lgm_tot_p16", "lgm_tot_p50",
            "lgm_tot_p84", "sfr_tot_p16", "sfr_tot_p50", "sfr_tot_p84", "specsfr_tot_p50", "ssfr",
            "bptclass"]
    return best[keep].rename(columns={"ra": "host_ra", "dec": "host_dec"})


# ---------------------------------------------------------------------------
# The galaxy population
# ---------------------------------------------------------------------------

def specific_sfr(t: pd.DataFrame) -> np.ndarray:
    """MPA-JHU's median specific SFR. Its specsfr_tot_p50 column is -9999 for
    some galaxies whose mass and SFR are both valid (a trap: those would all
    look "quenched"); there we use log SFR - log M instead."""
    s = t["specsfr_tot_p50"].to_numpy(float)
    alt = t["sfr_tot_p50"].to_numpy(float) - t["lgm_tot_p50"].to_numpy(float)
    return np.where(s > -20, s, alt)


def load_background() -> pd.DataFrame:
    parts = [pd.read_csv(f) for f in sorted(glob.glob(str(C.SDSS_BACKGROUND_DIR / "plates_*.csv")))]
    g = pd.concat([p for p in parts if len(p)], ignore_index=True)
    g = g[(g["lgm_tot_p50"] > 5) & (g["sfr_tot_p50"] > -10)]
    # Repeat spectra of the same galaxy: keep one (rounded position)
    key = (g["ra"].round(4).astype(str) + "_" + g["dec"].round(4).astype(str))
    g = g.loc[~key.duplicated()].copy()
    # Main Galaxy Sample targets only (primtarget bits 0x40 | 0x80 | 0x100):
    # the flux-limited sample the hosts' spectra also come from.
    g = g[(g["primtarget"].astype(np.int64) & 0x1C0) != 0]
    g = g[(g["z"] >= C.HOST_Z_MIN) & (g["z"] <= C.HOST_Z_MAX)].copy()
    g["ssfr"] = specific_sfr(g)
    return g.reset_index(drop=True)


def zmatched_weights(gal_z: np.ndarray, base_w: np.ndarray, host_z: np.ndarray,
                     zbins=np.arange(C.HOST_Z_MIN, C.HOST_Z_MAX + 1e-9, 0.01)) -> np.ndarray:
    """Re-weight galaxies so their weighted redshift distribution matches the
    hosts': within each redshift slice the weights are scaled to the number of
    hosts found there."""
    gi = np.digitize(gal_z, zbins) - 1
    hi = np.digitize(host_z, zbins) - 1
    w = np.zeros_like(base_w, dtype=float)
    for k in range(len(zbins) - 1):
        nh = np.sum(hi == k)
        tot = base_w[gi == k].sum()
        if nh and tot > 0:
            w[gi == k] = base_w[gi == k] / tot * nh
    return w / w.sum()


def quenched_fraction(ssfr: np.ndarray, w: np.ndarray | None = None) -> float:
    q = ssfr < C.LOG_SSFR_QUENCHED
    return float(np.average(q, weights=w)) if w is not None else float(q.mean())


def binomial_interval(k: int, n: int, cl: float = 0.68) -> tuple[float, float]:
    """Clopper-Pearson interval for k successes in n trials."""
    from scipy.stats import beta
    a = (1 - cl) / 2
    lo = beta.ppf(a, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(1 - a, k + 1, n - k) if k < n else 1.0
    return float(lo), float(hi)


def fit_a_plus_b(hosts_ia: pd.DataFrame, gal: pd.DataFrame,
                 ssfr_bins=np.arange(-13.0, -8.49, 0.25)) -> dict:
    """Maximum-likelihood ratio B/A of the Type Ia "A + B" model from the
    sSFR distribution of Type Ia hosts.

    Model: the probability that a Type Ia host falls in sSFR bin j is
    proportional to sum over galaxies in that bin of (M + (B/A) SFR), with
    galaxies re-weighted to the hosts' redshift distribution. Only the ratio
    is constrained by where the hosts are; the absolute A and B need a rate."""
    gm = 10 ** gal["lgm_tot_p50"].to_numpy()
    gs = 10 ** gal["sfr_tot_p50"].to_numpy()
    gz = gal["z"].to_numpy()
    gss = gal["ssfr"].to_numpy()
    hz = hosts_ia["z_sn"].to_numpy()
    hss = hosts_ia["ssfr"].to_numpy()
    gidx = np.digitize(gss, ssfr_bins) - 1
    hidx = np.digitize(hss, ssfr_bins) - 1
    nb = len(ssfr_bins) - 1
    counts = np.bincount(hidx.clip(0, nb - 1), minlength=nb)[:nb]

    def model_p(log_ratio):
        w = zmatched_weights(gz, gm + 10 ** log_ratio * gs, hz)
        p = np.bincount(gidx.clip(0, nb - 1), weights=w, minlength=nb)[:nb]
        return p / p.sum()

    def nll(log_ratio):
        p = np.clip(model_p(log_ratio), 1e-12, None)
        return -np.sum(counts * np.log(p))

    res = minimize_scalar(nll, bounds=(7.0, 12.0), method="bounded")
    best = res.x
    grid = np.linspace(best - 1.0, best + 1.0, 201)
    prof = np.array([nll(g) for g in grid])
    inside = grid[prof - prof.min() < 0.5]
    return {"log_B_over_A_yr": float(best), "lo68": float(inside.min()), "hi68": float(inside.max()),
            "n_hosts": int(len(hosts_ia)), "ssfr_bins": ssfr_bins.tolist(), "counts": counts.tolist(),
            "model_counts": (model_p(best) * counts.sum()).tolist(),
            "mass_only_counts": (model_p(-5.0) * counts.sum()).tolist(),
            "sfr_only_counts": (model_p(15.0) * counts.sum()).tolist()}


def per_mass_rates(hosts_f: pd.DataFrame, gal: pd.DataFrame, log_b_over_a: float | None = None,
                   bins=(-13.0, -12.0, -11.5, -11.0, -10.5, -10.0, -9.5, -8.5)) -> pd.DataFrame:
    """Supernovae per unit stellar mass in bins of host specific SFR, relative
    to the average galaxy (mass-weighted, re-weighted to the hosts' redshifts):
    a value of 1 means "as many as the stars there would suggest"; a rate that
    follows stellar mass is flat, one that follows star formation rises in
    proportion to sSFR. With log_b_over_a, also the A + B model's prediction
    (rate per mass ~ 1 + (B/A) sSFR); without it, the SFR-only prediction."""
    bins = np.asarray(bins)
    w = zmatched_weights(gal["z"].to_numpy(), 10 ** gal["lgm_tot_p50"].to_numpy(),
                         hosts_f["z_sn"].to_numpy())
    gs = gal["ssfr"].to_numpy()
    s_lin = 10 ** gs
    model_w = w * (1 + 10 ** log_b_over_a * s_lin) if log_b_over_a is not None else w * s_lin
    model_w = model_w / model_w.sum()
    ig = np.digitize(gs, bins) - 1
    ih = np.digitize(hosts_f["ssfr"].to_numpy(), bins) - 1
    n_tot = len(hosts_f)
    rows = []
    for k in range(len(bins) - 1):
        n = int((ih == k).sum())
        fm = float(w[ig == k].sum())
        mid = float(np.average(gs[ig == k], weights=w[ig == k])) if fm > 0 else np.nan
        rel = (n / n_tot) / fm if fm > 0 else np.nan
        lo, hi = binomial_interval(n, n_tot)
        rows.append({"ssfr_lo": bins[k], "ssfr_hi": bins[k + 1], "ssfr_mid": mid, "n": n,
                     "mass_share": fm, "rel": rel, "rel_lo": lo / fm if fm else np.nan,
                     "rel_hi": hi / fm if fm else np.nan,
                     "model": float(model_w[ig == k].sum() / fm) if fm else np.nan})
    return pd.DataFrame(rows)
