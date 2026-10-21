"""Where giants live in the cosmic web.

Two independent looks:

1. **Oei et al.'s (2024) reconstruction.** BORG SDSS infers the full matter
   density field (dark matter included) from the SDSS galaxies by running
   cosmological initial conditions forward with gravity. Oei et al. read off
   the density around 260 luminous giants and 1,443 ordinary LoTSS DR1 radio
   galaxies (z < 0.16). We recompute their statistics from their published
   per-object values.

2. **Our own count of neighbours.** Around each host we count SDSS Main Galaxy
   Sample galaxies within 2 Mpc (projected) and +-1000 km/s: a direct, model-
   free tracer of how crowded a galaxy's surroundings are. The MGS is flux
   limited, so the expected count falls with redshift; we therefore only ever
   compare objects at the same redshift (and, where known, the same radio
   luminosity, because more powerful radio galaxies are known to live in
   denser places). We also compare with radio-quiet galaxies of the same
   redshift, brightness and colour.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

import catalogs
import config
import cosmology
from config import DATA_DIR

C_KMS = 299_792.458


# --------------------------------------------------------------------------
# 1. The BORG SDSS densities of Oei et al. (2024)
# --------------------------------------------------------------------------

def borg_densities() -> dict:
    """Per-object mean relative densities and the population statistics."""
    t = catalogs.oei2024()
    out = {}
    for method in ("fixed", "flex"):
        res = {}
        for key, df in t.items():
            mean = np.exp(df[f"mu_{method}"] + df[f"sigma2_{method}"] / 2)  # E[1 + Delta] of a lognormal
            k, _, theta = stats.gamma.fit(mean, floc=0)
            res[key] = {"n": int(len(df)), "mean": float(mean.mean()), "median": float(mean.median()),
                        "gamma_k": float(k), "gamma_theta": float(theta), "gamma_mean": float(k * theta),
                        "values": mean.to_numpy()}
        ks = stats.ks_2samp(res["giants"]["values"], res["radiogalaxies"]["values"])
        res["ks_p"] = float(ks.pvalue)
        out[method] = res
    # dynamical environment: mean probability of each class
    for key, df in t.items():
        out[f"web_{key}"] = {c: float(df[f"p_{c}"].mean()) for c in ("cluster", "filament", "sheet", "void")}
        out[f"bcg_fraction_{key}"] = float(df["is_bcg"].mean())
    out["published"] = {"fixed": {"giants": (2.5, 1.9, 4.8), "radiogalaxies": (2.1, 1.8, 3.8)},
                        "flex": {"giants": (4.7, 1.6, 7.5), "radiogalaxies": (3.6, 1.8, 6.5)},
                        "ks_p": "<~1e-6"}
    return out


# --------------------------------------------------------------------------
# 2. Neighbour counts in the SDSS Main Galaxy Sample
# --------------------------------------------------------------------------

MGS_CACHE = DATA_DIR / "derived" / "mgs.pkl"


@lru_cache(maxsize=None)
def mgs() -> pd.DataFrame:
    if MGS_CACHE.exists():
        return pd.read_pickle(MGS_CACHE)
    files = sorted((DATA_DIR / "sdss" / "mgs").glob("mgs_ra*.csv"))
    if not files:
        raise FileNotFoundError("SDSS Main Galaxy Sample missing; run fetch_data.py --only sdss")
    df = pd.concat([pd.read_csv(f, skiprows=1) for f in files], ignore_index=True)
    df = df.drop_duplicates("specObjID").reset_index(drop=True)
    df["r"] = df["petroMag_r"] - df["extinction_r"]
    df["g_r"] = df["dered_g"] - df["dered_r"]
    df.to_pickle(MGS_CACHE)
    return df


def _xyz(ra, dec) -> np.ndarray:
    ra, dec = np.radians(np.asarray(ra, float)), np.radians(np.asarray(dec, float))
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


@lru_cache(maxsize=1)
def _mgs_tree():
    m = mgs()
    return cKDTree(_xyz(m["ra"], m["dec"])), m["z"].to_numpy(), m["ra"].to_numpy(), m["dec"].to_numpy()


def neighbour_counts(ra, dec, z, *, radius_mpc: float = config.ENV_RADIUS_MPC,
                     dv_kms: float = config.ENV_DV_KMS) -> np.ndarray:
    """Number of MGS galaxies within ``radius_mpc`` (proper, projected) and
    +-dv of each position, excluding anything within 3" (the object itself)."""
    tree, mz, _, _ = _mgs_tree()
    ra, dec, z = (np.atleast_1d(np.asarray(v, float)) for v in (ra, dec, z))
    da = cosmology.angular_diameter_distance(z, config.COSMOLOGY)
    theta = radius_mpc / da  # radians
    chord = 2 * np.sin(theta / 2)
    self_chord = 2 * np.sin(np.radians(3 / 3600) / 2)
    pts = _xyz(ra, dec)
    out = np.zeros(len(z), int)
    for i in range(len(z)):
        idx = tree.query_ball_point(pts[i], chord[i])
        if not idx:
            continue
        idx = np.asarray(idx)
        dz = np.abs(mz[idx] - z[i]) * C_KMS / (1 + z[i])
        sep = np.linalg.norm(tree.data[idx] - pts[i], axis=1)
        out[i] = int(((dz < dv_kms) & (sep > self_chord)).sum())
    return out


def in_footprint(ra, dec, z, *, frac: float = 0.8) -> np.ndarray:
    """Whether the 2-Mpc circle is mostly inside the MGS footprint: we test a
    ring of 16 points at the full radius and require >= ``frac`` to have an
    MGS galaxy within 0.5 deg (the survey's coverage, not the structure)."""
    tree, _, _, _ = _mgs_tree()
    ra, dec, z = (np.atleast_1d(np.asarray(v, float)) for v in (ra, dec, z))
    theta = np.degrees(config.ENV_RADIUS_MPC / cosmology.angular_diameter_distance(z, config.COSMOLOGY))
    ok = np.zeros(len(z), bool)
    probe = 2 * np.sin(np.radians(0.5) / 2)
    for i in range(len(z)):
        ang = np.linspace(0, 2 * np.pi, 16, endpoint=False)
        pr = ra[i] + theta[i] * np.sin(ang) / np.cos(np.radians(dec[i]))
        pd_ = dec[i] + theta[i] * np.cos(ang)
        d, _ = tree.query(_xyz(pr, pd_), k=1, distance_upper_bound=probe)
        ok[i] = np.mean(np.isfinite(d)) >= frac
    return ok


def environment_samples() -> dict:
    """Hosts with SDSS spectroscopic redshifts in the MGS volume, with counts."""
    import hosts

    zlo, zhi = config.ENV_Z_RANGE
    g = hosts.giant_hosts()
    g = g[g["z_sdss"].between(zlo, zhi) & (np.abs(g["z_sdss"] - g["z"]) < 0.01)].copy()
    g["z_env"] = g["z_sdss"]
    rg = catalogs.radio_galaxies()
    c = rg[rg["has_sdss_z"] & (rg["L144"] >= config.RADIO_LOUD_W_HZ) & ~rg["is_giant_candidate"]
           & rg["z"].between(zlo, zhi)].copy()
    c["z_env"] = c["z"]
    out = {}
    for key, df, ra, dec in (("giants", g, "ra", "dec"), ("radiogalaxies", c, "host_ra", "host_dec")):
        df = df[in_footprint(df[ra], df[dec], df["z_env"])].copy()
        df["n_nb"] = neighbour_counts(df[ra], df[dec], df["z_env"])
        out[key] = df
    # radio-quiet look-alikes: MGS galaxies matched to each giant host that is
    # itself an MGS galaxy, in redshift, brightness and colour
    m = mgs()
    tree, mz, mra, mdec = _mgs_tree()
    gh = out["giants"]
    d, j = tree.query(_xyz(gh["ra"], gh["dec"]))
    in_mgs = d < 2 * np.sin(np.radians(2 / 3600) / 2)
    gh = gh.assign(host_r=np.where(in_mgs, m["r"].to_numpy()[j], np.nan),
                   host_gr=np.where(in_mgs, m["g_r"].to_numpy()[j], np.nan))
    out["giants"] = gh
    return out


def compare_environments(samples: dict | None = None, seed: int = 11) -> dict:
    s = samples or environment_samples()
    g, c = s["giants"], s["radiogalaxies"]
    rng = np.random.default_rng(seed)
    res = {"n_giants": int(len(g)), "n_radiogalaxies": int(len(c))}

    def matched_ratio(sel_giants, match_fn):
        """Mean neighbour count of giants over that of their matched ordinary
        radio galaxies (each giant weighted equally), with a bootstrap error."""
        num, den = [], []
        for _, row in sel_giants.iterrows():
            mask = match_fn(row)
            if mask.sum() >= 5:
                num.append(row["n_nb"])
                den.append(c.loc[mask, "n_nb"].mean())
        if not num:
            return None
        num, den = np.array(num, float), np.array(den, float)
        boot = [num[k].mean() / den[k].mean() for k in (rng.integers(0, num.size, num.size) for _ in range(2000))]
        return {"n_giants": int(num.size), "ratio": float(num.mean() / den.mean()), "err": float(np.std(boot)),
                "giant_mean": float(num.mean()), "control_mean": float(den.mean())}

    zc = c["z_env"].to_numpy()
    lc = np.log10(c["L144"].to_numpy())
    res["vs_radiogalaxies_z"] = matched_ratio(g, lambda r: np.abs(zc - r["z_env"]) < 0.005)
    gl = g[np.isfinite(g["L144"])]
    res["vs_radiogalaxies_z_L"] = matched_ratio(
        gl, lambda r: (np.abs(zc - r["z_env"]) < 0.01) & (np.abs(lc - np.log10(r["L144"])) < 0.2))
    # look-alikes
    m = mgs()
    mz, mr, mgr = m["z"].to_numpy(), m["r"].to_numpy(), m["g_r"].to_numpy()
    gq = g[np.isfinite(g["host_r"])]
    num, den = [], []
    for _, row in gq.iterrows():
        sel = np.where((np.abs(mz - row["z_env"]) < 0.005) & (np.abs(mr - row["host_r"]) < 0.25)
                       & (np.abs(mgr - row["host_gr"]) < 0.05))[0]
        if sel.size >= 5:
            pick = rng.choice(sel, size=min(sel.size, 30), replace=False)
            cnt = neighbour_counts(m["ra"].to_numpy()[pick], m["dec"].to_numpy()[pick], mz[pick])
            num.append(row["n_nb"])
            den.append(cnt.mean())
    if num:
        num, den = np.array(num, float), np.array(den, float)
        boot = [num[k].mean() / den[k].mean() for k in (rng.integers(0, len(num), len(num)) for _ in range(2000))]
        res["vs_lookalikes"] = {"n_giants": int(len(num)), "ratio": float(num.mean() / den.mean()),
                                "err": float(np.std(boot)), "giant_mean": float(num.mean()),
                                "control_mean": float(den.mean())}
    # does the environment change with length? (z-matched residual counts)
    resid = []
    for _, row in g.iterrows():
        mask = np.abs(zc - row["z_env"]) < 0.005
        if mask.sum() >= 5:
            resid.append((row["lp_mpc"], row["n_nb"] / max(c.loc[mask, "n_nb"].mean(), 1e-3)))
    if resid:
        lp, rr = np.array(resid).T
        rho, p = stats.spearmanr(lp, rr)
        res["length_vs_relative_density"] = {"n": int(lp.size), "spearman_rho": float(rho), "p": float(p),
                                             "ratio_lp_gt_2mpc": float(np.mean(rr[lp > 2])) if (lp > 2).any() else np.nan,
                                             "ratio_lp_lt_1mpc": float(np.mean(rr[lp < 1])) if (lp < 1).any() else np.nan,
                                             "n_gt_2mpc": int((lp > 2).sum())}
    return res


def borg_luminosity_check() -> dict:
    """Oei et al.'s own densities, with the control radio galaxies split by
    radio luminosity (from the LoTSS DR2 value-added catalogue, matched within
    3" of the host): are giants denser than radio galaxies *as luminous*?"""
    t = catalogs.oei2024()
    rg = catalogs.radio_galaxies()
    tree = cKDTree(_xyz(rg["host_ra"], rg["host_dec"]))
    lum = {}
    for key, df in t.items():
        d, i = tree.query(_xyz(df["ra"], df["dec"]))
        ok = d < 2 * np.sin(np.radians(3 / 3600) / 2)
        lum[key] = np.where(ok, rg["L144"].to_numpy()[i], np.nan)
    out = {"n_matched": {k: int(np.isfinite(v).sum()) for k, v in lum.items()}}
    for method in ("fixed", "flex"):
        g = np.exp(t["giants"][f"mu_{method}"] + t["giants"][f"sigma2_{method}"] / 2).to_numpy()
        r = np.exp(t["radiogalaxies"][f"mu_{method}"] + t["radiogalaxies"][f"sigma2_{method}"] / 2).to_numpy()
        bright = lum["radiogalaxies"] >= config.RADIO_LOUD_W_HZ
        faint = lum["radiogalaxies"] < config.RADIO_LOUD_W_HZ
        out[method] = {"giants_mean": float(g.mean()), "all_rg_mean": float(r.mean()),
                       "bright_rg_mean": float(r[bright].mean()), "n_bright": int(bright.sum()),
                       "faint_rg_mean": float(r[faint].mean()), "n_faint": int(faint.sum()),
                       "ks_p_giants_vs_bright": float(stats.ks_2samp(g, r[bright]).pvalue),
                       "ks_p_giants_vs_all": float(stats.ks_2samp(g, r).pvalue)}
    return out
