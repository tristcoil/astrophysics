"""The sign of the SDSS Fundamental Plane distances in Cosmicflows-4.

CF4's largest component, 34,045 SDSS early-type galaxies (61% of the
catalogue), comes from Howlett et al. (2022), who publish for each galaxy a
log-distance ratio eta = log10(d_z / d_H): d_z is the distance its (group)
redshift implies, d_H the distance the Fundamental Plane implies. A galaxy
with eta > 0 is *closer* than its redshift says, so its distance modulus is

    DM = DM_z - 5 eta.

In the CF4 table as published (journal electronic table, VizieR J/ApJ/944/94),
the SDSS distance moduli instead follow DM = DM_z(own redshift) + 4.98 eta +
0.012, with r = 0.9975 over all 34,045 galaxies: the deviation enters with the
opposite sign. This module documents the evidence and builds a corrected
catalogue:

1. The meaning of eta is fixed by geometry, not by any column description.
   Howlett's effective radii are computed at the group's redshift distance
   (we reproduce them to 0.0002 dex), and eta equals the Fundamental Plane
   residual in size (slope 1.02, r = 0.91): a galaxy that looks too big at its
   redshift distance has eta > 0, and must be closer than that distance.
2. Nusser & Tully (2026), with a later CF4 file from the Extragalactic Distance
   Database, find SDSS bulk-flow components opposite in sign to ours from the
   published table. With the sign corrected, our fits reproduce all 24 of their
   standardised All-Individual components to 0.24 sigma rms (``bulkflow``).

The corrected SDSS moduli are DM_z(z_group) - 4.98 eta + 0.012 + dzp, where
the zero-point shift dzp re-levels SDSS against the Tully-Fisher, SN Ia, SBF
and SN II distances of the same groups (CF4 set the SDSS level while the sign
was reversed). Galaxies with another distance method as well are recombined
by inverse variance, and groups are rebuilt from their members the same way
(which reproduces CF4's own group moduli to 0.004 mag).
"""

from __future__ import annotations

import functools

import numpy as np
import pandas as pd
from scipy import stats

import catalogs
import config
import cosmology as cg

OTHER_METHODS = ["DMtf", "DMsnIa", "DMsnII", "DMsbf", "DMtrgb", "DMceph", "DMmas"]
H0_REF = 75.0


@functools.lru_cache(maxsize=1)
def howlett() -> pd.DataFrame:
    path = config.DATA_DIR / "howlett2022" / "SDSS_PV_public.dat"
    with open(path) as fh:
        names = fh.readline().lstrip("#").split()
    return pd.read_csv(path, sep=r"\s+", comment="#", names=names)


@functools.lru_cache(maxsize=1)
def _matched() -> pd.DataFrame:
    """CF4 galaxies joined to Howlett's catalogue by PGC number."""
    h = howlett()[["PGC", "zcmb", "zhelio", "zcmb_group", "logdist", "logdist_err", "logdist_corr",
                   "logdist_corr_err", "IDgroupT17", "r", "s", "i", "deVRad_r", "deVAB_r"]]
    g = catalogs.cf4_galaxies().merge(h, on="PGC", how="left", suffixes=("", "_h"))
    g["in_sdss"] = g["logdist_corr"].notna() & g["DMfp"].notna()
    return g


def evidence() -> dict:
    """The numbers behind the correction (all recomputed from the files)."""
    g = _matched()
    h = howlett()
    s = g[g["in_sdss"]]
    dmz_own = cg.dm_from_redshift(s["zcmb"], s["zhelio"], H0_REF)
    slope, icpt, r, _, se = stats.linregress(s["logdist_corr"], s["DMfp"] - dmz_own)
    res = {"n_matched": int(len(s)), "cf4_vs_eta": {"slope": float(slope), "slope_err": float(se),
                                                   "intercept": float(icpt), "r": float(r),
                                                   "rms": float(np.std(s["DMfp"] - dmz_own - icpt - slope * s["logdist_corr"]))}}
    # Howlett's sizes: r = log10(theta sqrt(b/a) D_A(z_group) x 1000) in h^-1 kpc.
    theta = np.radians(h["deVRad_r"] / 3600.0) * np.sqrt(h["deVAB_r"])
    for zc in ("zcmb", "zcmb_group"):
        da = cg.comoving_distance(h[zc]) / (1 + h[zc])
        d = h["r"] - np.log10(theta * da * 1000.0)
        res[f"size_at_{zc}"] = {"median": float(np.median(d)), "mad": float(np.median(np.abs(d - np.median(d))))}
    # eta against the Fundamental Plane residual.
    A = np.column_stack([h["s"], h["i"], np.ones(len(h))])
    coef, *_ = np.linalg.lstsq(A, h["r"], rcond=None)
    resid = h["r"] - A @ coef
    sl, _, rr, _, _ = stats.linregress(resid, h["logdist_corr"])
    res["fp_plane"] = {"a": float(coef[0]), "b": float(coef[1]), "c": float(coef[2]), "rms": float(resid.std())}
    res["eta_vs_fp_residual"] = {"slope": float(sl), "r": float(rr)}
    # Within rich groups, eta does not follow each member's own redshift.
    big = h[h["IDgroupT17"] > 0].groupby("IDgroupT17").filter(lambda x: len(x) >= 10)
    x = pd.Series(np.log10(cg.comoving_distance(big["zcmb"])), index=big.index)
    x = x - x.groupby(big["IDgroupT17"]).transform("mean")
    y = big["logdist_corr"] - big.groupby("IDgroupT17")["logdist_corr"].transform("mean")
    sl2, _, _, _, se2 = stats.linregress(x, y)
    res["eta_within_group_slope_vs_own_logdz"] = {"slope": float(sl2), "err": float(se2),
                                                   "n_gal": int(len(big)), "n_groups": int(big["IDgroupT17"].nunique())}
    res["eta_mean_by_z"] = h.groupby(pd.cut(h["zcmb"], np.arange(0, 0.105, 0.01)), observed=True)["logdist_corr"].mean().round(4).tolist()
    return res


def _sdss_pure(g: pd.DataFrame) -> pd.Series:
    """SDSS galaxies whose CF4 modulus rests on the Fundamental Plane alone."""
    return g["in_sdss"] & g[OTHER_METHODS].isna().all(axis=1)


def _uncalibrated_sdss_dm(g: pd.DataFrame) -> pd.Series:
    ev = evidence()["cf4_vs_eta"]
    dmz_grp = cg.dm_from_redshift(g["zcmb_group"], g["zhelio"], H0_REF)
    return dmz_grp + ev["intercept"] - ev["slope"] * g["logdist_corr"]


def _combine_methods(row_dm: dict) -> tuple[float, float]:
    vals = np.array([v for v, e in row_dm.values()])
    w = np.array([1 / e**2 for v, e in row_dm.values()])
    return float(np.sum(w * vals) / w.sum()), float(1 / np.sqrt(w.sum()))


def shared_group_offset(dm_sdss: pd.Series, g: pd.DataFrame | None = None) -> dict:
    """Weighted mean of (SDSS-FP group modulus - other-method group modulus) over shared groups."""
    g = _matched() if g is None else g
    pure = _sdss_pure(g)
    others = ["DMtf", "DMsnIa", "DMsnII", "DMsbf"]
    has_other = g[others].notna().any(axis=1)
    rows = []
    for gid, idx in g.groupby("1PGC").groups.items():
        sub = g.loc[idx]
        f = sub[pure.loc[idx]]
        o = sub[has_other.loc[idx]]
        if len(f) == 0 or len(o) == 0:
            continue
        wf = 1 / f["e_DM"] ** 2
        dmf = float(np.sum(wf * dm_sdss.loc[f.index]) / wf.sum())
        ef = float(1 / np.sqrt(wf.sum()))
        vals, ws = [], []
        for _, row in o.iterrows():
            v, e = _combine_methods({m: (row[m], row["e_" + m]) for m in others if pd.notna(row[m])})
            vals.append(v)
            ws.append(1 / e**2)
        vals, ws = np.array(vals), np.array(ws)
        rows.append((dmf - np.sum(ws * vals) / ws.sum(), np.hypot(ef, 1 / np.sqrt(ws.sum())), sub["Vcmb"].median()))
    t = pd.DataFrame(rows, columns=["d", "e", "cz"])
    w = 1 / t["e"] ** 2
    mean = float(np.sum(w * t["d"]) / w.sum())
    return {"n_groups": int(len(t)), "offset": mean, "err": float(1 / np.sqrt(w.sum())),
            "chi2_dof": float(np.sum(w * (t["d"] - mean) ** 2) / (len(t) - 1))}


@functools.lru_cache(maxsize=1)
def zero_point() -> dict:
    g = _matched()
    pure = _sdss_pure(g)
    published = shared_group_offset(g["DM"].where(pure), g)
    raw = _uncalibrated_sdss_dm(g)
    corrected = shared_group_offset(raw.where(pure), g)
    return {"published": published, "corrected_uncalibrated": corrected, "dzp": -corrected["offset"],
            "dzp_err": corrected["err"]}


@functools.lru_cache(maxsize=1)
def corrected_galaxies() -> pd.DataFrame:
    """CF4 Table 2 with sign-corrected, re-levelled SDSS moduli (columns DM, e_DM replaced;
    the published values kept as DM_published, DMfp_published)."""
    g = _matched().copy()
    zp = zero_point()["dzp"]
    new_fp = _uncalibrated_sdss_dm(g) + zp
    g["DM_published"], g["DMfp_published"] = g["DM"], g["DMfp"]
    sd = g["in_sdss"]
    g.loc[sd, "DMfp"] = new_fp[sd]
    pure = _sdss_pure(g)
    g.loc[pure, "DM"] = g.loc[pure, "DMfp"]
    mixed = sd & ~pure
    for i in np.flatnonzero(mixed.to_numpy()):
        row = g.iloc[i]
        parts = {m: (row[m], row["e_" + m]) for m in ["DMfp"] + OTHER_METHODS if pd.notna(row[m]) and row["e_" + m] > 0}
        g.iloc[i, g.columns.get_loc("DM")] = _combine_methods(parts)[0]
    g["sdss_corrected"] = sd
    return g


@functools.lru_cache(maxsize=1)
def corrected_groups() -> pd.DataFrame:
    """CF4 groups with moduli rebuilt from the corrected galaxies (inverse-variance means).

    Rebuilt from the published galaxies the same way, group moduli match CF4's
    Table 4 to 0.004 mag; that small constant offset is added back so the
    distance scale is exactly CF4's.
    """
    grp = catalogs.cf4_groups().copy()
    gal = corrected_galaxies()
    w = 1 / gal["e_DM"] ** 2
    agg = pd.DataFrame({"sw": w, "swdm": w * gal["DM"], "swdm0": w * gal["DM_published"],
                        "sdss": gal["sdss_corrected"].astype(int), "grp": gal["1PGC"]}).groupby("grp").sum()
    rebuilt0 = agg["swdm0"] / agg["sw"]
    rebuilt = agg["swdm"] / agg["sw"]
    m = grp.set_index("1PGC")
    offset = float(np.nanmedian(m["DMzp"] - rebuilt0.reindex(m.index)))
    grp["DMzp_published"] = grp["DMzp"]
    new = (rebuilt + offset).reindex(grp["1PGC"]).to_numpy()
    has_sdss = agg["sdss"].reindex(grp["1PGC"]).fillna(0).to_numpy() > 0
    grp.loc[has_sdss, "DMzp"] = new[has_sdss]
    grp["sdss_corrected"] = has_sdss
    grp["Dist"] = 10 ** ((grp["DMzp"] - 25.0) / 5.0)
    grp.attrs["rebuild_offset"] = offset
    return grp
