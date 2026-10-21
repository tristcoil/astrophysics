"""All measurements, from the cached catalogues to ``data/results.json``.

1. :func:`derive`: luminosities, rest-frame colours, sizes, masses.
2. :func:`colour_divider`: where the red sequence ends and the blue cloud
   begins (project D2's fit, redone on the same data as a consistency check).
3. Samples: :func:`volume_limit` and :func:`mass_completeness` make fair
   censuses out of a flux-limited survey.
4. :func:`redshift_bias`: why the votes need debiasing.
5. :func:`wrong_colours`: red spirals and blue ellipticals, with the
   definitions of Masters et al. (2010) and Schawinski et al. (2009).
6. :func:`bars`: bar fractions, Masters et al. (2011)'s way and ours.
7. :func:`spiral_arms`: arm numbers, Hart et al. (2016)'s way.
8. :func:`environment_trends`: morphology-density and colour-density at
   fixed stellar mass, and which of the two drives the other.
9. :func:`experts`: the volunteers against Nair & Abraham (2010).

Published values (``literature.py``) are only consulted for grading.
"""

from __future__ import annotations

import json
import time

import numpy as np

import environment as env
import morphology as mo
import physics as ph
from config import (
    COMPLETENESS_R, MASS_BINS, MASS_COMPLETENESS_LOSS, N_DENSITY_BINS, RESULTS_PATH, SEED,
    VL_Z_RANGE,
)

H70 = 0.70


# ---------------------------------------------------------------------------
# 1. Derived quantities
# ---------------------------------------------------------------------------

def derive(c: dict) -> dict:
    """Everything that follows from redshift, magnitudes and image measurements.

    Magnitudes are corrected for Milky Way dust. ``M`` is the rest-frame
    r-band absolute magnitude for Planck 2018's H0, exactly as in project D2;
    ``M70`` the same for H0 = 70, the convention of most papers compared with.
    Stellar masses are MPA-JHU's, which assume H0 = 70, as published.
    """
    z = c["z"].astype(float)
    good = (c["in_sdss"] & (c["zWarning"] == 0) & np.isfinite(z) & (z > 0.0005) & (z < 0.25)
            & np.isfinite(c["petroMag_r"]) & np.isfinite(c["dered_g"]) & np.isfinite(c["dered_r"]))
    zz = np.where(good, z, 0.05)  # placeholder keeps the formulas finite; masked by `good`
    r = c["petroMag_r"] - c["extinction_r"]
    gr = c["dered_g"] - c["dered_r"]
    ur = c["dered_u"] - c["dered_r"]
    gi = c["dered_g"] - c["dered_i"]
    kr = ph.k_correction("r", zz, gr, "gr")
    d = {
        "good": good,
        "z": z,
        "r": r,
        "gr_obs": gr,
        "kr": kr,
        "M": ph.absolute_magnitude(r, zz, kr),
        "gr0": ph.rest_frame_colour("g", "r", gr, zz),
        "ur0": ph.rest_frame_colour("u", "r", ur, zz),
        "gi0": ph.rest_frame_colour("g", "i", gi, zz),
        "R50_kpc": c["petroR50_r"] * ph.kpc_per_arcsec(zz),
        "R90_arcsec": c["petroR90_r"].astype(float),
        "ba": c["expAB_r"].astype(float),
        "fracdev": c["fracDeV_r"].astype(float),
        "mass": c["lgm_tot_p50"].astype(float),
        "ssfr": c["specsfr_tot_p50"].astype(float),
        "u_err": c["modelMagErr_u"].astype(float),
    }
    d["M70"] = d["M"] + ph.h_magnitude_shift(H70)
    d["log_ab"] = -np.log10(np.clip(d["ba"], 1e-3, None))
    for k in ("M", "M70", "gr0", "ur0", "gi0", "R50_kpc"):
        d[k] = np.where(good, d[k], np.nan)
    return d


# ---------------------------------------------------------------------------
# 2. The colour divider (project D2's method, on D2's sample)
# ---------------------------------------------------------------------------

def colour_divider(M, gr0, sel, bin_width: float = 0.5) -> dict:
    """The colour that best separates red and blue galaxies, against luminosity.

    Identical to D2's ``color_divider``: in each 0.5-mag bin two Gaussians are
    fitted to the rest-frame g - r distribution, the divider is where a galaxy
    is equally likely to belong to either, and a straight line is fitted
    through the bins where the two populations leave a real dip between them.
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
        gm = GaussianMixture(2, random_state=20260929).fit(gr0[m].reshape(-1, 1))
        mu = gm.means_.ravel()
        blue, red = np.argsort(mu)
        grid = np.linspace(mu[blue], mu[red], 2001).reshape(-1, 1)
        post_red = gm.predict_proba(grid)[:, red]
        valley = float(grid[np.argmin(np.abs(post_red - 0.5)), 0])
        dens = np.exp(gm.score_samples(np.array([[mu[blue]], [valley], [mu[red]]])))
        bins.append({"M_center": float((lo + hi) / 2), "n": int(m.sum()), "divider": valley,
                     "bimodal": bool(dens[1] < 0.9 * min(dens[0], dens[2]))})
    use = [b for b in bins if b["bimodal"]]
    slope, intercept = np.polyfit([b["M_center"] for b in use], [b["divider"] for b in use], 1)
    return {"slope": float(slope), "intercept": float(intercept), "bins": bins}


def is_red(M, gr0, divider: dict) -> np.ndarray:
    return np.asarray(gr0) > divider["slope"] * np.asarray(M) + divider["intercept"]


# ---------------------------------------------------------------------------
# 3. Fair censuses
# ---------------------------------------------------------------------------

def red_kcorrection(d: dict, z_edges=np.arange(0.0, 0.2001, 0.005), pct: float = 99.5):
    """K_r(z) for the reddest galaxies: the faint-end limit of any volume-limited cut.

    At each redshift the reddest galaxies have the largest K-corrections, so
    they are the first to drop below the flux limit. We take the 99.5th
    percentile of observed g - r in each redshift slice.
    """
    g = d["good"]
    z, gr = d["z"][g], d["gr_obs"][g]
    centres, k = [], []
    for lo, hi in zip(z_edges[:-1], z_edges[1:]):
        m = (z >= lo) & (z < hi)
        if m.sum() < 200:
            continue
        zc = 0.5 * (lo + hi)
        centres.append(zc)
        k.append(float(ph.k_correction("r", zc, np.percentile(gr[m], pct), "gr")))
    return np.array(centres), np.array(k)


def faint_limit(z, kz) -> np.ndarray:
    """Faintest M_r (Planck) seen at redshift z by a sample complete to COMPLETENESS_R."""
    zc, kc = kz
    return COMPLETENESS_R - ph.distance_modulus(z) - np.interp(z, zc, kc)


def volume_limit(d: dict, kz, z_range=VL_Z_RANGE) -> dict:
    """The main volume-limited sample: complete for M_r < M_lim over z_range."""
    m_lim = float(faint_limit(z_range[1], kz))
    sel = d["good"] & (d["z"] >= z_range[0]) & (d["z"] < z_range[1]) & (d["M"] < m_lim)
    return {"z_range": list(z_range), "M_lim": m_lim, "M70_lim": m_lim + ph.h_magnitude_shift(H70),
            "n": int(sel.sum()), "sel": sel}


def mass_completeness(d: dict, kz, bins=MASS_BINS, loss: float = MASS_COMPLETENESS_LOSS,
                      z_floor: float = 0.02, z_cap: float = VL_Z_RANGE[1]) -> list[dict]:
    """For each stellar-mass bin, the redshift range over which it is complete.

    At fixed mass the faintest galaxies are the reddest (old stars give little
    light per unit mass). For a bin starting at mass m, find the absolute
    magnitude that all but the faintest ``loss`` of its galaxies are brighter
    than, measured at low redshift where the census is complete; the bin is
    complete out to the redshift where that magnitude meets the flux limit.
    Iterated, because each new limit adds galaxies to the reference sample.
    """
    out = []
    g = d["good"] & np.isfinite(d["mass"])
    for lo, hi in zip(bins[:-1], bins[1:]):
        near_lo = g & (d["mass"] >= lo) & (d["mass"] < lo + 0.1)
        z_max = 0.03
        for _ in range(5):
            ref = near_lo & (d["z"] < z_max) & (d["z"] >= 0.005)
            if ref.sum() < 50:
                break
            m_faint = float(np.quantile(d["M"][ref], 1.0 - loss))
            zs = np.linspace(0.005, 0.25, 4901)
            ok = faint_limit(zs, kz) > m_faint
            new = float(zs[ok][-1]) if ok.any() else 0.005
            if abs(new - z_max) < 1e-3:
                z_max = new
                break
            z_max = new
        z_hi = min(z_max, z_cap)
        sel = g & (d["mass"] >= lo) & (d["mass"] < hi) & (d["z"] >= z_floor) & (d["z"] < z_hi)
        out.append({"mass_lo": lo, "mass_hi": hi, "M_faint": m_faint, "z_complete": z_max,
                    "z_range": [z_floor, z_hi], "usable": bool(z_hi > z_floor + 0.005),
                    "n": int(sel.sum()), "sel": sel})
    return out


# ---------------------------------------------------------------------------
# 4. Classification bias
# ---------------------------------------------------------------------------

def redshift_bias(c: dict, d: dict, kz, vl: dict, *, m_range=(-21.8, -21.4),
                  r50_range=(3.0, 6.0)) -> dict:
    """Distant galaxies look smoother: the bias, and how well each correction removes it.

    (a) One kind of galaxy at different distances: galaxies of one luminosity
        and one physical size should have the same mix of shapes at every
        redshift here (the range spans ~1 Gyr of cosmic time), so any trend of
        their mean vote fractions with redshift is a bias of the classification.
        The redshift range stops where the luminosity bin stops being complete.
    (b) The flatness test: in the volume-limited sample, the fraction of
        galaxies in each class should not depend on redshift either. Each class
        is read three ways: the raw (consistency-weighted) votes, Willett et
        al. (2013)'s debiasing and Hart et al. (2016)'s.
    """
    zs = np.linspace(0.01, 0.2, 1901)
    z_hi = float(zs[faint_limit(zs, kz) > m_range[1]][-1])
    sel = (d["good"] & (d["M"] >= m_range[0]) & (d["M"] < m_range[1])
           & (d["R50_kpc"] >= r50_range[0]) & (d["R50_kpc"] < r50_range[1]))
    edges = np.arange(0.02, z_hi + 1e-9, 0.01)
    if z_hi - edges[-1] > 0.004:
        edges = np.append(edges, z_hi)
    else:
        edges[-1] = z_hi
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = sel & (d["z"] >= lo) & (d["z"] < hi)
        row = {"z": float(np.median(d["z"][m])), "z_lo": float(lo), "z_hi": float(hi),
               "n": int(m.sum())}
        for name, key in (("raw", "wf_features"), ("hart16", "p_features"),
                          ("willett13", "w13_features")):
            v = c[key][m].astype(float)
            row[name] = float(np.nanmean(v))
            row[name + "_err"] = float(np.nanstd(v) / np.sqrt(np.isfinite(v).sum()))
        rows.append(row)
    one_kind = {"M_range": list(m_range), "R50_kpc_range": list(r50_range), "z_max": z_hi,
                "bins": rows}

    # (b) Class fractions against redshift in the volume-limited sample.
    readers = {"raw": "wf_", "willett13": "w13_", "hart16": "p_"}
    nq = c["n_bar"] + c["n_nobar"]
    v = vl["sel"]
    z_edges = np.array([0.03, 0.04, 0.05, 0.06, 0.07, 0.085])
    flat = {}
    for rname, pre in readers.items():
        g = {k: np.nan_to_num(c[pre + k].astype(float)) for k in
             ("smooth", "features", "notedgeon", "spiral", "bar")}
        disk = g["features"] * g["notedgeon"]
        classes = {
            "smooth (p > 0.8)": (v, g["smooth"] > 0.8),
            "features (p > 0.8)": (v, g["features"] > 0.8),
            "spiral (product > 0.5)": (v, disk * g["spiral"] > 0.5),
            "barred, of disks": (v & (disk > 0.5) & (nq >= 5) & (c["expAB_r"] > 0.5),
                                 g["bar"] >= 0.5),
            "bar majority, of oblique disks": (v & mo.classes(c)["oblique"], g["bar"] >= 0.5),
        }
        out = {}
        for cname, (base, flag) in classes.items():
            fr = []
            for lo, hi in zip(z_edges[:-1], z_edges[1:]):
                m = base & (d["z"] >= lo) & (d["z"] < hi)
                fr.append(frac(flag[m]))
            f = np.array([x["f"] for x in fr])
            out[cname] = {"bins": fr, "drift": float(f[-1] - f[0]),
                          "rel_spread": float((f.max() - f.min()) / f.mean())}
        flat[rname] = out
    return {"one_kind": one_kind, "flatness": {"z_edges": z_edges, "classes": flat}}


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def frac(flags) -> dict:
    p, e, n = mo.binomial_fraction(flags)
    return {"f": p, "err": e, "n": n}


def _clean(x):
    """Make numpy objects JSON-friendly."""
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items() if not isinstance(v, np.ndarray)
                or v.ndim == 0 or v.size < 2000}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return _clean(x.tolist())
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, float) and not np.isfinite(x):
        return None
    return x


def save_results(results: dict) -> None:
    out = {k: v for k, v in results.items() if not k.startswith("_")}
    out["_generated_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    RESULTS_PATH.write_text(json.dumps(_clean(out), indent=1))


# ---------------------------------------------------------------------------
# 5. Galaxies of the "wrong" colour
# ---------------------------------------------------------------------------

def _gaussian_red_sequence(colour, M, sel, bin_width: float = 0.33, clip: float = 2.5):
    """Mean and width of the red sequence in bins of magnitude (iteratively clipped).

    Schawinski et al. (2009) fitted a Gaussian to the u - r colours of early
    types in 0.33-mag bins. Blue outliers would drag a plain fit bluer, so the
    fit is repeated keeping galaxies within 2.5 sigma, which converges on the
    sequence itself.
    """
    x, y = np.asarray(M)[sel], np.asarray(colour)[sel]
    edges = np.arange(np.floor(x.min() / bin_width) * bin_width, x.max() + bin_width, bin_width)
    centres, mu, sd, n = [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        v = y[(x >= lo) & (x < hi)]
        if v.size < 30:
            continue
        m, s = np.median(v), 1.4826 * np.median(np.abs(v - np.median(v)))
        for _ in range(20):
            keep = np.abs(v - m) < clip * s
            m_new, s_new = v[keep].mean(), v[keep].std(ddof=1)
            if abs(m_new - m) < 1e-5 and abs(s_new - s) < 1e-5:
                break
            m, s = m_new, s_new
        centres.append(0.5 * (lo + hi)); mu.append(m); sd.append(s); n.append(v.size)
    return np.array(centres), np.array(mu), np.array(sd), np.array(n)


def wrong_colours(c: dict, d: dict, cls: dict, divider: dict, vl: dict, kz) -> dict:
    """How many spirals are red, and how many ellipticals blue?"""
    red = is_red(d["M"], d["gr0"], divider)
    v = vl["sel"]
    pp = mo.path_probabilities(c)
    out = {"divider": {"slope": divider["slope"], "intercept": divider["intercept"]}}

    # Our census, in the main volume-limited sample, with D2's red/blue line.
    out["vl"] = {
        "smooth_blue": frac(~red[v & cls["smooth"]]),
        "features_red": frac(red[v & cls["features"]]),
        "spiral_red": frac(red[v & cls["spiral"]]),
        "faceon_spiral_red": frac(red[v & cls["spiral"] & (d["log_ab"] < 0.2)]),
        "edgeon_red": frac(red[v & cls["edgeon"]]),
        "red_smooth_fraction": frac(cls["smooth"][v & red]),
        "red_features_fraction": frac(cls["features"][v & red]),
        "blue_features_fraction": frac(cls["features"][v & ~red]),
        # Likelihood-weighted: every galaxy counts, in proportion to its votes.
        "red_fraction_of_smooth_weighted": mo.weighted_fraction(red[v], c["p_smooth"][v]),
        "red_fraction_of_features_weighted": mo.weighted_fraction(red[v], c["p_features"][v]),
    }

    # Masters et al. (2010): face-on, disk-dominated spirals with visible arms.
    def masters10(m_cut: float, z_range) -> dict:
        s = (d["good"] & (d["z"] >= z_range[0]) & (d["z"] < z_range[1]) & (d["M70"] < m_cut)
             & (pp["spiral"] >= 0.8) & (d["log_ab"] < 0.2) & (d["fracdev"] <= 0.5))
        line = 0.63 - 0.02 * (d["M70"] + 20.0)
        is_red_m10 = d["gr0"] > line
        bar_ok = s & cls["disk"]
        return {"z_range": list(z_range), "M70_cut": m_cut, "n": int(s.sum()),
                "red": frac(is_red_m10[s]),
                "bar_fraction_red": frac(cls["barred"][bar_ok & is_red_m10]),
                "bar_fraction_blue": frac(cls["barred"][bar_ok & ~is_red_m10]),
                "ssfr_median_red": float(np.nanmedian(d["ssfr"][s & is_red_m10])),
                "ssfr_median_blue": float(np.nanmedian(d["ssfr"][s & ~is_red_m10])),
                "_sel": s, "_red": is_red_m10}

    m_lim70 = vl["M70_lim"]
    # (a) their luminosity cut, over the redshifts where GZ2 is complete for it;
    zs = np.linspace(0.01, 0.1, 901)
    z_ok = zs[faint_limit(zs, kz) + ph.h_magnitude_shift(H70) > -20.17]
    out["masters10_their_cut"] = masters10(-20.17, (0.02, float(z_ok[-1])))
    # (b) their redshift range, with the brighter luminosity cut GZ2 needs there.
    out["masters10_their_z"] = masters10(m_lim70, (0.03, 0.085))

    # Schawinski et al. (2009): blue early types, 3 sigma below the red sequence in u - r.
    s = (d["good"] & (d["z"] >= 0.02) & (d["z"] < 0.05) & (d["M70"] < -20.7) & cls["smooth"]
         & np.isfinite(d["ur0"]))
    centres, mu, sd, n = _gaussian_red_sequence(d["ur0"], d["M70"], s)
    mu_i = np.interp(d["M70"], centres, mu)
    sd_i = np.interp(d["M70"], centres, sd)
    blue_etg = s & (d["ur0"] < mu_i - 3.0 * sd_i)
    out["schawinski09"] = {"n_early": int(s.sum()), "n_blue": int(blue_etg.sum()),
                           "blue": frac(blue_etg[s]),
                           "red_sequence": {"M70": centres, "mu": mu, "sd": sd, "n": n},
                           "_sel": s, "_blue": blue_etg}
    # The same question with the volunteers' raw votes, as GZ1 used them.
    s_raw = (d["good"] & (d["z"] >= 0.02) & (d["z"] < 0.05) & (d["M70"] < -20.7)
             & (np.nan_to_num(c["wf_smooth"]) > 0.8) & np.isfinite(d["ur0"]))
    out["schawinski09"]["blue_raw_votes"] = frac((d["ur0"] < mu_i - 3.0 * sd_i)[s_raw])
    return out


# ---------------------------------------------------------------------------
# 6. Bars
# ---------------------------------------------------------------------------

def bars(c: dict, d: dict, cls: dict, divider: dict, vl: dict, mass_info: list) -> dict:
    out = {}
    # Masters et al. (2011), with their exact cuts on the final GZ2 catalogue.
    n_q = c["n_bar"] + c["n_nobar"]
    s = (d["good"] & (d["z"] > 0.01) & (d["z"] < 0.06) & (d["M70"] < -19.38)
         & (d["log_ab"] < 0.3) & (n_q >= 10))
    barred_m11 = c["n_bar"] >= c["n_nobar"]
    out["masters11"] = {"n": int(s.sum()), "barred": frac(barred_m11[s])}
    colour_edges = np.arange(0.2, 1.01, 0.05)
    rows = []
    for lo, hi in zip(colour_edges[:-1], colour_edges[1:]):
        m = s & (d["gr0"] >= lo) & (d["gr0"] < hi)
        if m.sum() >= 30:
            rows.append({"gr": 0.5 * (lo + hi), **frac(barred_m11[m])})
    out["masters11"]["vs_colour"] = rows
    out["masters11"]["blue_disks"] = frac(barred_m11[s & (d["gr0"] < 0.6)])
    out["masters11"]["red_disks"] = frac(barred_m11[s & (d["gr0"] > 0.7)])
    out["masters11"]["red_bulgy_disks"] = frac(barred_m11[s & (d["gr0"] > 0.7)
                                                          & (d["fracdev"] > 0.5)])

    # Ours: debiased votes, bar-eligible disks of the main volume-limited sample.
    v = vl["sel"] & cls["disk"]
    out["vl"] = {"n": int(v.sum()), "barred": frac(cls["barred"][v]),
                 "barred_p03": frac((np.nan_to_num(c["p_bar"]) >= 0.3)[v]),
                 "barred_raw": frac((np.nan_to_num(c["wf_bar"]) >= 0.5)[v])}
    red = is_red(d["M"], d["gr0"], divider)
    out["vl"]["red"] = frac(cls["barred"][v & red])
    out["vl"]["blue"] = frac(cls["barred"][v & ~red])
    rows = []
    for lo, hi in zip(colour_edges[:-1], colour_edges[1:]):
        m = v & (d["gr0"] >= lo) & (d["gr0"] < hi)
        if m.sum() >= 30:
            rows.append({"gr": 0.5 * (lo + hi), **frac(cls["barred"][m])})
    out["vl"]["vs_colour"] = rows
    rows = []
    for b in mass_info:
        if not b["usable"]:
            continue
        m = b["sel"] & cls["disk"]
        rows.append({"mass": 0.5 * (b["mass_lo"] + b["mass_hi"]),
                     "mass_range": [b["mass_lo"], b["mass_hi"]],
                     "all": frac(cls["barred"][m]), "red": frac(cls["barred"][m & red]),
                     "blue": frac(cls["barred"][m & ~red])})
    out["vs_mass"] = rows
    return out


# ---------------------------------------------------------------------------
# 7. Spiral arms
# ---------------------------------------------------------------------------

def spiral_arms(c: dict, d: dict, cls: dict, sigma: np.ndarray | None = None) -> dict:
    """Hart et al. (2016)'s samples and arm-number demographics."""
    arm, _ = mo.arm_number(c)
    lum = d["good"] & (d["z"] > 0.03) & (d["z"] < 0.085) & (d["M70"] <= -21.0)
    mass_lim = lum & (d["mass"] >= 10.6)
    out = {"n_lum_limited": int(lum.sum()), "n_mass_limited": int(mass_lim.sum())}
    for name, base in (("lum", lum), ("mass", mass_lim)):
        sp = base & cls["spiral"] & (arm >= 0)
        counts = np.bincount(arm[sp], minlength=5)
        out[f"n_spirals_{name}"] = int(sp.sum())
        out[f"arm_fractions_{name}"] = {lab: float(counts[i] / sp.sum())
                                        for i, lab in enumerate(mo.ARM_LABELS)}
        out[f"m2_{name}"] = frac(arm[sp] == 1)
    sp = mass_lim & cls["spiral"] & (arm >= 0)
    out["_spirals_mass"] = sp
    out["_arm"] = arm
    # Against stellar mass: 0.1-dex bins (Hart et al.: 20 bins).
    edges = np.quantile(d["mass"][sp], np.linspace(0, 1, 21))
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = sp & (d["mass"] >= lo) & (d["mass"] <= hi)
        rows.append({"mass": float(np.median(d["mass"][m])), "n": int(m.sum()),
                     **{f"f{lab}": frac(arm[m] == i)["f"] for i, lab in enumerate(mo.ARM_LABELS)},
                     **{f"e{lab}": frac(arm[m] == i)["err"] for i, lab in enumerate(mo.ARM_LABELS)}})
    out["vs_mass"] = rows
    # Against g - i colour: 20 equal-count bins, as Hart et al.
    gi = d["gi0"]
    ok = sp & np.isfinite(gi)
    edges = np.quantile(gi[ok], np.linspace(0, 1, 21))
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = ok & (gi >= lo) & (gi <= hi)
        rows.append({"gi": float(np.median(gi[m])), "gi_sd": float(np.std(gi[m])), "n": int(m.sum()),
                     **{f"f{lab}": frac(arm[m] == i)["f"] for i, lab in enumerate(mo.ARM_LABELS)},
                     **{f"e{lab}": frac(arm[m] == i)["err"] for i, lab in enumerate(mo.ARM_LABELS)}})
    out["vs_colour"] = rows
    out["mean_gi_by_arms"] = {lab: float(np.nanmean(gi[ok & (arm == i)]))
                              for i, lab in enumerate(mo.ARM_LABELS)}
    out["sd_gi_by_arms"] = {lab: float(np.nanstd(gi[ok & (arm == i)]))
                            for i, lab in enumerate(mo.ARM_LABELS)}
    # Mass-matched: are many-armed spirals bluer at the same stellar mass?
    rng = np.random.default_rng(SEED)
    two = np.flatnonzero(ok & (arm == 1))
    many = np.flatnonzero(ok & (arm >= 2))
    mbins = np.arange(10.6, 11.8, 0.1)
    matched = []
    for lo, hi in zip(mbins[:-1], mbins[1:]):
        a = many[(d["mass"][many] >= lo) & (d["mass"][many] < hi)]
        b = two[(d["mass"][two] >= lo) & (d["mass"][two] < hi)]
        if len(a) and len(b):
            matched.append(rng.choice(b, size=len(a), replace=True))
    matched = np.concatenate(matched)
    out["mass_matched_gi"] = {"many_armed": float(np.nanmean(gi[many])),
                              "two_armed_matched": float(np.nanmean(gi[matched])),
                              "difference": float(np.nanmean(gi[many]) - np.nanmean(gi[matched])),
                              "pooled_sd": float(np.nanstd(gi[ok & (arm >= 1)]))}
    if sigma is not None:
        okd = sp & np.isfinite(sigma)
        edges = np.quantile(np.log10(sigma[okd]), np.linspace(0, 1, 11))
        rows = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = okd & (np.log10(sigma) >= lo) & (np.log10(sigma) <= hi)
            rows.append({"log_sigma": float(np.median(np.log10(sigma[m]))), **frac(arm[m] == 1)})
        out["m2_vs_density"] = rows
    return out


# ---------------------------------------------------------------------------
# 8. Environment: morphology-density versus colour-density
# ---------------------------------------------------------------------------

def densities(c: dict, d: dict, mgs: dict) -> dict:
    """Sigma_5 for every GZ2 galaxy inside the tracer volume (D2's measure)."""
    from config import TRACER_Z_RANGE

    m = env.derive_mgs(mgs)
    tracers = env.select_tracers(m)
    mask = env.build_mask(m["ra"], m["dec"])
    tgt = d["good"] & (d["z"] > TRACER_Z_RANGE[0]) & (d["z"] < TRACER_Z_RANGE[1])
    ra, dec = c["ra_sdss"][tgt], c["dec_sdss"][tgt]
    xyz = ph.cartesian(ra, dec, ph.comoving_distance_h(d["z"][tgt]))
    res = env.local_density(xyz, d["z"][tgt], c["specobjid"][tgt],
                            m["xyz_h"][tracers], m["specObjID"][tracers])
    usable = env.usable_density(res, d["z"][tgt], ra, dec, mask)
    sigma = np.full(len(d["z"]), np.nan)
    idx = np.flatnonzero(tgt)
    sigma[idx[usable]] = res["sigma"][usable]
    # The GZ2 galaxies that are themselves tracers: a check that our Sigma_5 is D2's.
    is_tracer = np.isin(c["specobjid"], m["specObjID"][tracers])
    return {"sigma": sigma, "n_tracers": int(tracers.sum()), "band_offset": m["band_offset"],
            "n_target": int(tgt.sum()), "n_usable": int(usable.sum()),
            "n_lower_bound": int(res["lower_bound"].sum()), "is_tracer": is_tracer,
            "mask_area_deg2": float(mask.interior.sum() * mask.hp.pixel_area.to_value("deg2")),
            "_mgs": m, "_tracers": tracers, "_mask": mask}


def _cell_stats(red, ws, smooth_clean, feat_clean, w=None) -> dict:
    """Colour and morphology fractions in one (mass, density) cell.

    ``ws`` is each galaxy's early-type likelihood p_smooth / (p_smooth +
    p_features): galaxies count fractionally, as in Bamford et al. (2009).
    """
    red = red.astype(float)
    wd = 1.0 - ws
    out = {"n": int(red.size), "f_red": float(red.mean()), "w_smooth": float(ws.mean()),
           "f_smooth_clean": float(smooth_clean.mean()), "f_features_clean": float(feat_clean.mean()),
           "red_given_smooth": float(np.sum(ws * red) / np.sum(ws)),
           "red_given_disk": float(np.sum(wd * red) / np.sum(wd)),
           "smooth_given_red": float(np.sum(ws * red) / max(red.sum(), 1e-9)),
           "smooth_given_blue": float(np.sum(ws * (1 - red)) / max((1 - red).sum(), 1e-9))}
    if smooth_clean.sum():
        out["red_given_smooth_clean"] = float(red[smooth_clean].mean())
    if feat_clean.sum():
        out["red_given_features_clean"] = float(red[feat_clean].mean())
    return out


def _decompose(lo: dict, hi: dict) -> dict:
    """Exact two-way split of a change between the sparsest and densest bins.

    f_red = w f_red|E + (1 - w) f_red|D, with w the early-type fraction. The
    change from low to high density splits exactly (evaluating each factor at
    the midpoint) into a *morphology-mix* term -- more early types, which are
    redder -- and a *colour-at-fixed-morphology* term -- each type itself
    redder. The converse split does the same for the early-type fraction.
    """
    w0, w1 = lo["w_smooth"], hi["w_smooth"]
    rE0, rE1, rD0, rD1 = lo["red_given_smooth"], hi["red_given_smooth"], lo["red_given_disk"], hi["red_given_disk"]
    wbar = 0.5 * (w0 + w1)
    morph = (w1 - w0) * (0.5 * (rE0 + rE1) - 0.5 * (rD0 + rD1))
    colour = wbar * (rE1 - rE0) + (1 - wbar) * (rD1 - rD0)
    f0, f1 = lo["f_red"], hi["f_red"]
    sR0, sR1, sB0, sB1 = lo["smooth_given_red"], hi["smooth_given_red"], lo["smooth_given_blue"], hi["smooth_given_blue"]
    fbar = 0.5 * (f0 + f1)
    colour_mix = (f1 - f0) * (0.5 * (sR0 + sR1) - 0.5 * (sB0 + sB1))
    morph_fixed = fbar * (sR1 - sR0) + (1 - fbar) * (sB1 - sB0)
    return {"delta_red": f1 - f0, "red_from_morphology_mix": morph,
            "red_at_fixed_morphology": colour,
            "delta_smooth": w1 - w0, "smooth_from_colour_mix": colour_mix,
            "smooth_at_fixed_colour": morph_fixed}


def environment_trends(c: dict, d: dict, cls: dict, divider: dict, mass_info: list,
                       sigma: np.ndarray, n_boot: int = 300) -> dict:
    red = is_red(d["M"], d["gr0"], divider)
    ps, pf = np.nan_to_num(c["p_smooth"]), np.nan_to_num(c["p_features"])
    ws = np.where(ps + pf > 0, ps / np.maximum(ps + pf, 1e-9), 0.5)
    usable = np.isfinite(sigma)
    in_any = np.zeros(len(sigma), bool)
    for b in mass_info:
        if b["usable"]:
            in_any |= b["sel"]
    logs = np.log10(np.where(usable, sigma, np.nan))
    edges = np.quantile(logs[in_any & usable], np.linspace(0, 1, N_DENSITY_BINS + 1))
    rng = np.random.default_rng(SEED)
    out = {"log_sigma_edges": edges, "bins": []}
    for b in mass_info:
        if not b["usable"]:
            continue
        sel = b["sel"] & usable
        k = np.clip(np.searchsorted(edges, logs, side="right") - 1, 0, N_DENSITY_BINS - 1)
        cells = []
        for j in range(N_DENSITY_BINS):
            m = sel & (k == j)
            st = _cell_stats(red[m], ws[m], cls["smooth"][m], cls["features"][m])
            st["log_sigma_median"] = float(np.median(logs[m]))
            cells.append(st)
        dec = _decompose(cells[0], cells[-1])
        # Bootstrap the decomposition: resample galaxies within the two end cells.
        lo_i = np.flatnonzero(sel & (k == 0))
        hi_i = np.flatnonzero(sel & (k == N_DENSITY_BINS - 1))
        boots = []
        for _ in range(n_boot):
            a = rng.choice(lo_i, lo_i.size)
            bb = rng.choice(hi_i, hi_i.size)
            boots.append(_decompose(
                _cell_stats(red[a], ws[a], cls["smooth"][a], cls["features"][a]),
                _cell_stats(red[bb], ws[bb], cls["smooth"][bb], cls["features"][bb])))
        dec_err = {key: float(np.std([bt[key] for bt in boots])) for key in dec}
        out["bins"].append({"mass_range": [b["mass_lo"], b["mass_hi"]], "z_range": b["z_range"],
                            "n": int(sel.sum()), "cells": cells, "decomposition": dec,
                            "decomposition_err": dec_err})
    return out


# ---------------------------------------------------------------------------
# 9. Volunteers versus experts
# ---------------------------------------------------------------------------

def experts(c: dict, d: dict, na10, max_sep_arcsec: float = 2.0) -> dict:
    """GZ2 against Nair & Abraham (2010)'s single-expert classifications."""
    from scipy.spatial import cKDTree

    def unit(r, de):
        r, de = np.radians(np.asarray(r, float)), np.radians(np.asarray(de, float))
        return np.column_stack([np.cos(de) * np.cos(r), np.cos(de) * np.sin(r), np.sin(de)])

    ra = np.where(np.isfinite(c["ra_sdss"]), c["ra_sdss"], c["ra"])
    dec = np.where(np.isfinite(c["dec_sdss"]), c["dec_sdss"], c["dec"])
    tree = cKDTree(unit(ra, dec))
    chord = 2 * np.sin(np.radians(max_sep_arcsec / 3600) / 2)
    dist, idx = tree.query(unit(na10["_RA"], na10["_DE"]), distance_upper_bound=chord)
    hit = np.isfinite(dist)
    gi = idx[hit]
    tt = na10["TT"].to_numpy()[hit]
    bar = na10["Bar"].to_numpy()[hit].astype(int)
    out = {"n_na10": int(len(na10)), "n_matched": int(hit.sum())}
    E, S0, late = tt == -5, (tt >= -3) & (tt <= 0), (tt >= 1) & (tt < 99)
    for name, key in (("hart16", "p_"), ("willett13", "w13_"), ("raw", "wf_")):
        sm = np.asarray(c[key + "smooth"][gi], float)
        fe = np.asarray(c[key + "features"][gi], float)
        clean = sm > 0.8
        out[name] = {
            "E_median_psmooth": float(np.nanmedian(sm[E])),
            "E_psmooth_sd": float(np.nanstd(sm[E])),
            "E_max_pfeatures": float(np.nanmax(fe[E])),
            "late_median_pfeatures": float(np.nanmedian(fe[late])),
            "late_pfeatures_sd": float(np.nanstd(fe[late])),
            "smooth08_n": int(clean.sum()),
            "smooth08_late": float(late[clean].mean()),
            "smooth08_late_or_s0": float((late | S0)[clean].mean()),
        }
    # Bars, among oblique disks (Willett et al. 2013, Table 3 thresholds).
    oblique = ((np.nan_to_num(c["wf_features"][gi]) > 0.430)
               & (np.nan_to_num(c["wf_notedgeon"][gi]) >= 0.715)
               & ((c["n_bar"][gi] + c["n_nobar"][gi]) >= 20))
    barred = bar > 0
    strong, inter, weak = (bar & 2) > 0, (bar & 4) > 0, (bar & 8) > 0
    out["n_na10_barred"] = int((na10["Bar"] > 0).sum())
    out["n_oblique"] = int(oblique.sum())
    for name, key in (("hart16", "p_bar"), ("willett13", "w13_bar"), ("raw", "wf_bar")):
        pb = np.nan_to_num(np.asarray(c[key][gi], float))
        o = oblique
        out[name]["bars"] = {
            "barred_below_05": float((pb[o & barred] < 0.5).mean()),
            "unbarred_above_05": float((pb[o & ~barred] > 0.5).mean()),
            "strong_above_08": float((pb[o & strong] > 0.8).mean()),
            "intermediate_above_08": float((pb[o & inter] > 0.8).mean()),
            "strong_above_05": float((pb[o & strong] >= 0.5).mean()),
            "intermediate_above_05": float((pb[o & inter] >= 0.5).mean()),
            "weak_above_08": float((pb[o & weak] > 0.8).mean()),
            "weak_above_05": float((pb[o & weak] >= 0.5).mean()),
            "unbarred_above_03": float((pb[o & ~barred] >= 0.3).mean()),
            "strong_or_inter_above_03": float((pb[o & (strong | inter)] >= 0.3).mean()),
            "weak_above_03": float((pb[o & weak] >= 0.3).mean()),
        }
    out["_match"] = {"gz2_index": gi, "tt": tt, "bar": bar, "oblique": oblique}
    return out
