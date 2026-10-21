"""Peculiar velocities of the Cosmicflows-4 groups, and maps of them.

A group's peculiar velocity is what is left of its redshift after the
expansion is taken out: positive means it recedes faster than its distance
predicts (it is moving away from us), negative that it lags (moving toward us).
On a map, galaxies falling toward a great mass show up as positive on our side
of it and negative behind it; a void shows the opposite, a hole that galaxies
stream out of.

We use the Watkins & Feldman (2015) estimator with the natural logarithm, with
CF4's own H0 = 75 and curvature factor. The groups are *placed* at their
redshift distance, not their measured distance: a 20% distance error moves a
group by 20 h^-1 Mpc, but its redshift places it within a few h^-1 Mpc (its
peculiar velocity), so maps drawn this way do not smear structures along the
line of sight (Nusser & Tully 2026).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import catalogs
import config
import cosmology as cg
from config import C_KMS

#: Velocity noise that is not distance error: random motions of galaxies
#: within their groups and nonlinear flows (km/s). Added in quadrature when
#: averaging velocities, so that a few groups with tiny errors cannot dominate.
SIGMA_STAR = 250.0
#: Outlier rejection for maps: genuine peculiar velocities allowed for, and the
#: number of standard deviations beyond which a distance is a catastrophic error.
MAP_SIGMA_V = 400.0
MAP_CLIP = 4.0


def group_velocities(version: str = "corrected") -> pd.DataFrame:
    """Every CF4 group with our peculiar velocities and redshift-space position.

    ``version`` "corrected" uses group distances rebuilt with the SDSS sign fixed
    (``cf4fix``); "published" uses CF4's Table 4 as it stands.
    """
    if version == "corrected":
        import cf4fix

        g = cf4fix.corrected_groups().copy()
    else:
        g = catalogs.cf4_groups().copy()
    H0 = config.CF4_H0_VPEC
    d = 10 ** ((g["DMzp"] - 25.0) / 5.0)  # luminosity distance, Mpc
    cz = g["V3k"].to_numpy(float)
    f = cg.cf4_curvature_factor(cz / C_KMS)
    g["d_mpc"] = d
    g["v_direct"] = cg.vpec_direct(cz, d, H0)
    g["v_wf"] = cg.vpec_wf(cz, d, H0)
    g["v_wf_log10"] = cg.vpec_wf(cz, d, H0, log=np.log10)
    g["v_ramp"] = cg.vpec_cf4_ramp(cz, g["Vls"], d, H0)
    g["v_ramp_log10"] = cg.vpec_cf4_ramp(cz, g["Vls"], d, H0, log=np.log10)
    # First-order error of the WF estimator: dv = cz_mod * (ln 10 / 5) * e_DM.
    g["e_v"] = f * cz * cg.A_LN10_5 * g["e_DMzp"]
    g["r_z"] = cg.comoving_distance(g["z_cmb"].clip(lower=0), config.PLANCK18.Om0)
    for ax in "xyz":
        g[f"SG{ax.upper()}_z"] = g["r_z"] * g[f"sg{ax}"]
    # Local Sheet frame, for nearby maps: the same estimator with Vls instead of V3k.
    g["v_ls"] = cg.vpec_wf(g["Vls"].to_numpy(float), d, H0)
    g["ok"] = (cz > 0) & np.isfinite(g["v_wf"]) & (g["e_DMzp"] > 0)
    # Catastrophic distance errors: a distance modulus more than 4 sigma from
    # the redshift's, where sigma combines the distance error and 400 km/s of
    # genuine motion. Kept out of the maps only; the bulk-flow fits follow
    # Nusser & Tully and use everything.
    with np.errstate(divide="ignore", invalid="ignore"):
        dmu = g["DMzp"] - (5 * np.log10(f * cz / H0) + 25.0)
        sig = np.sqrt(g["e_DMzp"] ** 2 + (5 / np.log(10) * MAP_SIGMA_V / np.maximum(cz, 1.0)) ** 2)
    g["outlier"] = (np.abs(dmu) > MAP_CLIP * sig) | ~np.isfinite(dmu)
    g["map_ok"] = g["ok"] & ~g["outlier"]
    return g


def estimator_check(g: pd.DataFrame | None = None) -> dict:
    """Reproduce CF4's three published velocity columns and show the ln 10 issue.

    Table 4's Vpds column is Eq. 9 with H0 = 75. Its Vpwf column equals
    Eq. 10 evaluated with a base-10 logarithm, whereas Watkins & Feldman's
    estimator uses the natural logarithm (their derivation, -cz ln(1 - v/cz) ~ v,
    holds only for ln). The published values are therefore smaller by ln 10 =
    2.303 beyond the 3,000 km/s ramp. Two independent checks decide which is
    right: the direct estimate, which needs no logarithm at all, agrees with
    the natural-log version for small velocities; and the scatter of
    well-measured groups matches the direct estimate only for ln.
    """
    g = group_velocities("published") if g is None else g
    ok = g["ok"] & (g["V3k"] > 1000)
    res = {}
    for col, ours in [("Vpds", "v_direct"), ("Vpwf", "v_wf_log10"), ("Vpec", "v_ramp_log10")]:
        dv = (g.loc[ok, ours] - g.loc[ok, col]).abs()
        res[f"{col}_reproduced"] = {"median_abs_diff": float(dv.median()), "p99_abs_diff": float(dv.quantile(0.99)),
                                    "max_abs_diff": float(dv.max()), "ours": ours}
    dv = (g.loc[ok, "v_wf"] - g.loc[ok, "Vpwf"]).abs()
    res["Vpwf_vs_natural_log"] = {"median_abs_diff": float(dv.median())}
    # Ratio to the direct estimate where the velocity is small against cz (first order).
    small = ok & (g["v_direct"].abs() < 0.05 * g["V3k"]) & (g["v_direct"].abs() > 100)
    res["ratio_to_direct_small_v"] = {
        "n": int(small.sum()),
        "natural_log": float(np.median(g.loc[small, "v_wf"] / g.loc[small, "v_direct"])),
        "cf4_column": float(np.median(g.loc[small, "Vpwf"] / g.loc[small, "v_direct"])),
    }
    good = g["ok"] & (g["e_DMzp"] <= config.HUBBLE_GOOD_EDM) & (g["V3k"] > config.CF4_RAMP_KMS)

    def spread(v):
        v = np.asarray(v)
        return {"std": float(np.std(v)), "robust": float(1.4826 * np.median(np.abs(v - np.median(v)))),
                "median": float(np.median(v))}

    res["well_measured_beyond_ramp"] = {
        "n": int(good.sum()),
        "cf4_Vpec_column": spread(g.loc[good, "Vpec"]),
        "direct": spread(g.loc[good, "v_direct"]),
        "wf_natural_log": spread(g.loc[good, "v_wf"]),
        "expected_error_median": float(np.median(g.loc[good, "e_v"])),
        "published_sigma": config.PUBLISHED["cf4_vpec_sigma"]["value"],
    }
    return res


# --------------------------------------------------------------------------
# Sky maps
# --------------------------------------------------------------------------

def _healpix_sg_vectors(nside: int) -> np.ndarray:
    import healpy as hp

    x, y, z = hp.pix2vec(nside, np.arange(hp.nside2npix(nside)))
    return np.stack([x, y, z], axis=1)  # interpreted as supergalactic vectors


def smoothed_sky(vec_sg: np.ndarray, values: np.ndarray, errors: np.ndarray, *, nside: int, radius_deg: float,
                 min_count: int = 8) -> dict:
    """Error-weighted mean of ``values`` within ``radius_deg`` of each HEALPix pixel.

    Weights are 1 / (error^2 + SIGMA_STAR^2). Pixels with fewer than
    ``min_count`` objects in their disc are left empty (NaN).
    """
    from scipy.spatial import cKDTree

    pix = _healpix_sg_vectors(nside)
    tree = cKDTree(vec_sg)
    chord = 2 * np.sin(np.radians(radius_deg) / 2)
    w = 1.0 / (errors**2 + SIGMA_STAR**2)
    mean = np.full(len(pix), np.nan)
    err = np.full(len(pix), np.nan)
    count = np.zeros(len(pix), int)
    for i, nb in enumerate(tree.query_ball_point(pix, chord)):
        count[i] = len(nb)
        if len(nb) >= min_count:
            ww = w[nb]
            mean[i] = np.sum(ww * values[nb]) / np.sum(ww)
            err[i] = 1.0 / np.sqrt(np.sum(ww))
    return {"mean": mean, "err": err, "count": count}


def sky_maps(g: pd.DataFrame | None = None) -> list[dict]:
    """Smoothed peculiar-velocity sky maps, one per redshift-distance shell."""
    g = group_velocities() if g is None else g
    out = []
    for lo, hi in config.SKY_SHELLS:
        s = g[g["map_ok"] & (g["r_z"] >= lo) & (g["r_z"] < hi)]
        vec = s[["sgx", "sgy", "sgz"]].to_numpy()
        m = smoothed_sky(vec, s["v_wf"].to_numpy(), s["e_v"].to_numpy(), nside=config.SKY_NSIDE,
                         radius_deg=config.SKY_SMOOTH_DEG)
        m.update({"shell": (lo, hi), "n": len(s), "frame": "CMB"})
        out.append(m)
    # The nearest shell again in the Local Sheet frame, where our own motion
    # (shared by our neighbours) no longer dominates and local flows show.
    lo, hi = config.SKY_SHELLS[0]
    s = g[g["map_ok"] & (g["Vls"] > 300) & (g["r_z"] >= 3.0) & (g["r_z"] < hi)]
    m = smoothed_sky(s[["sgx", "sgy", "sgz"]].to_numpy(), s["v_ls"].to_numpy(), s["e_v"].to_numpy(),
                     nside=config.SKY_NSIDE, radius_deg=config.SKY_SMOOTH_DEG)
    m.update({"shell": (3.0, hi), "n": len(s), "frame": "Local Sheet"})
    out.append(m)
    return out


def sky_dipole_check(g: pd.DataFrame | None = None) -> list[dict]:
    """Direction of the largest positive smoothed velocity in each shell: where
    the groups of that shell are, on average, moving away from us fastest."""
    maps = sky_maps(g)
    import healpy as hp

    out = []
    for m in maps:
        i = int(np.nanargmax(m["mean"]))
        j = int(np.nanargmin(m["mean"]))
        vmax = hp.pix2vec(config.SKY_NSIDE, i)
        vmin = hp.pix2vec(config.SKY_NSIDE, j)
        out.append({"shell": m["shell"], "frame": m["frame"],
                    "max": float(m["mean"][i]), "max_sg": [float(a) for a in catalogs.lonlat(np.array(vmax))],
                    "max_gal": [float(a) for a in catalogs.lonlat(catalogs.sg_to_gal(np.array(vmax)))],
                    "min": float(m["mean"][j]), "min_sg": [float(a) for a in catalogs.lonlat(np.array(vmin))],
                    "min_gal": [float(a) for a in catalogs.lonlat(catalogs.sg_to_gal(np.array(vmin)))]})
    return out


def landmarks() -> pd.DataFrame:
    """Named clusters with Galactic and supergalactic positions at their CMB-frame redshift distance."""
    rows = []
    for lm in config.LANDMARKS:
        l, b = catalogs.equatorial_to_galactic(lm.ra, lm.dec)
        n = catalogs.unit_vector(l, b)
        z_cmb = catalogs.helio_to_frame(lm.cz_helio / C_KMS, n, config.SUN_CMB)
        r = float(cg.comoving_distance(z_cmb, config.PLANCK18.Om0)[0])
        sg = catalogs.gal_to_sg(n)[0]
        sgl, sgb = catalogs.lonlat(sg)
        rows.append({"name": lm.name, "l": float(l[0]), "b": float(b[0]), "SGL": float(sgl), "SGB": float(sgb),
                     "cz_cmb": float(z_cmb[0] * C_KMS), "r": r, "SGX": r * sg[0], "SGY": r * sg[1], "SGZ": r * sg[2]})
    return pd.DataFrame(rows)
