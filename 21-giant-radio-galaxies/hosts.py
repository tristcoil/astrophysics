"""The galaxies and black holes that make giants.

A galaxy's central black hole cannot be weighed directly beyond the nearest
few dozen galaxies, but its mass tracks how fast the stars in the galaxy's
bulge move: the M-sigma relation. SDSS spectra measure that stellar velocity
dispersion sigma for most of the brighter galaxies they observed. With
Kormendy & Ho's (2013) calibration,

    M_BH = 3.09e8 Msun (sigma / 200 km/s)^4.38    (scatter 0.29 dex),

a galaxy with sigma = 250 km/s holds a black hole of ~8e8 Msun.

Question: do giants need bigger black holes than ordinary radio galaxies? We
compare the hosts of giants with the hosts of radio-loud LoTSS galaxies that
are *not* giants, matched in redshift (SDSS observes different galaxies at
different distances, so an unmatched comparison would compare survey
selection, not physics) and, where we know it, in radio luminosity.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd
from scipy import stats

import catalogs
import config
import fetch_data
from config import DATA_DIR, MSIGMA


def mbh_from_sigma(sigma_kms) -> np.ndarray:
    """Black hole mass (Msun) from stellar velocity dispersion, KH13 Eq. 7."""
    return MSIGMA["norm_1e9"] * 1e9 * (np.asarray(sigma_kms, float) / 200.0) ** MSIGMA["slope"]


def _read_batches(tag: str) -> pd.DataFrame:
    files = sorted((DATA_DIR / "sdss" / f"spec_{tag}").glob("batch*.csv"))
    if not files:
        raise FileNotFoundError(f"no SDSS spectrum batches for {tag}; run fetch_data.py --only sdss")
    frames = [pd.read_csv(f, skiprows=1) for f in files]
    frames = [f for f in frames if len(f)]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


@lru_cache(maxsize=None)
def spectra_for(tag: str) -> pd.DataFrame:
    """One row per queried position: the nearest good SDSS spectrum, if any.

    Batches record each position's index in the query list; the list is
    rebuilt here exactly as the fetch built it, and every match is re-checked
    by position (within 3.5"), so a changed list cannot silently mismatch.
    """
    pos = fetch_data.host_positions()[tag]
    ra = np.array([p[0] for p in pos])
    dec = np.array([p[1] for p in pos])
    sp = _read_batches(tag)
    sp = sp[(sp["zWarning"] == 0)].copy()
    sp = sp.sort_values(["id", "distance"]).drop_duplicates("id")
    ids = sp["id"].to_numpy(int)
    ok = ids < len(pos)
    sp = sp[ok]
    ids = ids[ok]
    dra = (sp["ra"].to_numpy() - ra[ids]) * np.cos(np.radians(dec[ids])) * 3600
    ddec = (sp["dec"].to_numpy() - dec[ids]) * 3600
    sp = sp[np.hypot(dra, ddec) < 3.5]
    out = pd.DataFrame({"ra": ra, "dec": dec})
    sp = sp.set_index("id")
    for c in ("z", "class", "subClass", "velDisp", "velDispErr", "snMedian", "instrument", "specObjID"):
        out[c if c != "class" else "spec_class"] = sp[c].reindex(out.index)
    out["has_spectrum"] = out["z"].notna()
    good = ((out["spec_class"] == "GALAXY") & out["velDisp"].between(*config.SIGMA_RANGE_KMS)
            & (out["velDispErr"] > 0) & (out["velDispErr"] < config.SIGMA_MAX_REL_ERR * out["velDisp"]))
    out["sigma_ok"] = good
    out["mbh"] = np.where(good, mbh_from_sigma(out["velDisp"]), np.nan)
    out["log_mbh"] = np.log10(out["mbh"])
    return out


def giant_hosts() -> pd.DataFrame:
    """All known giants with their hosts' SDSS spectra and black hole masses."""
    g = catalogs.all_giants().reset_index(drop=True)
    sp = spectra_for("giants")
    g = g.join(sp[["z", "spec_class", "velDisp", "velDispErr", "sigma_ok", "mbh", "log_mbh", "has_spectrum"]]
               .rename(columns={"z": "z_sdss"}))
    # luminosity and stellar mass from the value-added catalogue, where the
    # giant is a catalogue source (those found by LOFAR Galaxy Zoo)
    rg = catalogs.radio_galaxies().set_index("source")
    names = g["name"].fillna("").str.strip()
    g["L144"] = rg["L144"].reindex(names).to_numpy()
    g["logmass_vac"] = rg["logmass"].reindex(names).to_numpy()
    return g


def control_hosts() -> pd.DataFrame:
    """Radio-loud, resolved LoTSS galaxies smaller than giants, with SDSS spectra."""
    rg = catalogs.radio_galaxies()
    sel = (rg["has_sdss_z"] & (rg["L144"] >= config.RADIO_LOUD_W_HZ) & ~rg["is_giant_candidate"]
           & (rg["z"] < 0.7))
    rg = rg[sel].reset_index(drop=True)
    sp = spectra_for("radiogalaxies")
    return rg.join(sp[["spec_class", "velDisp", "velDispErr", "sigma_ok", "mbh", "log_mbh"]])


def compare_with_oei2023() -> dict:
    """Our DR17 black hole masses against Oei et al.'s DR12-based ones."""
    o = catalogs.oei2023()
    g = giant_hosts()
    from scipy.spatial import cKDTree

    def xyz(ra, dec):
        ra, dec = np.radians(ra), np.radians(dec)
        return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])

    t = cKDTree(xyz(g["ra"], g["dec"]))
    d, i = t.query(xyz(o["ra"], o["dec"]))
    close = d < np.radians(2 / 3600)
    ours = g["mbh"].to_numpy()[i]
    theirs = o["mbh_1e9"].to_numpy() * 1e9
    both = close & np.isfinite(ours) & np.isfinite(theirs)
    r = np.log10(ours[both] / theirs[both])
    return {"n_theirs": int(np.isfinite(theirs).sum()), "n_both": int(both.sum()),
            "median_log_ratio": float(np.median(r)), "scatter_dex": float(1.4826 * np.median(np.abs(r - np.median(r)))),
            "their_median_1e9": float(np.nanmedian(theirs) / 1e9),
            "their_p10_p90_1e9": [float(v) for v in np.nanpercentile(theirs, [10, 90]) / 1e9]}


Z_BINS = ((0.02, 0.15), (0.15, 0.3), (0.3, 0.45), (0.45, 0.7))


def compare_black_holes() -> dict:
    """Giant hosts against ordinary radio galaxy hosts, in redshift bins, and
    for the subset with luminosities, matched in luminosity too."""
    g = giant_hosts()
    c = control_hosts()
    gz = g["z"].to_numpy()
    out = {"bins": []}
    for lo, hi in Z_BINS:
        a = g.loc[(gz >= lo) & (gz < hi) & g["sigma_ok"], "log_mbh"].dropna().to_numpy()
        b = c.loc[(c["z"] >= lo) & (c["z"] < hi) & c["sigma_ok"], "log_mbh"].dropna().to_numpy()
        row = {"z": [lo, hi], "n_giants": int(a.size), "n_control": int(b.size)}
        if a.size >= 10 and b.size >= 10:
            row.update({"giants_median": float(np.median(a)), "control_median": float(np.median(b)),
                        "difference_dex": float(np.median(a) - np.median(b)),
                        "difference_err_dex": float(_boot_median_diff(a, b)),
                        "ks_p": float(stats.ks_2samp(a, b).pvalue)})
        out["bins"].append(row)
    # luminosity- and redshift-matched (giants known to the value-added catalogue)
    gm = g[g["sigma_ok"] & np.isfinite(g["L144"])]
    cm = c[c["sigma_ok"]]
    diffs, pairs = [], []
    rng = np.random.default_rng(3)
    lc = np.log10(cm["L144"].to_numpy())
    zc = cm["z"].to_numpy()
    mc = cm["log_mbh"].to_numpy()
    for _, row in gm.iterrows():
        m = (np.abs(zc - row["z"]) < 0.02) & (np.abs(lc - np.log10(row["L144"])) < 0.15)
        if m.sum() >= 3:
            diffs.append(row["log_mbh"] - np.median(mc[m]))
            pairs.append(int(m.sum()))
    diffs = np.array(diffs)
    if diffs.size:
        boot = [np.median(rng.choice(diffs, diffs.size)) for _ in range(2000)]
        out["matched"] = {"n_giants": int(diffs.size), "median_difference_dex": float(np.median(diffs)),
                          "err_dex": float(np.std(boot)), "median_controls_per_giant": float(np.median(pairs))}
    # does a heavier black hole make a longer giant?
    s = g[g["sigma_ok"] & (g["z"] < 0.5)]
    rho, p = stats.spearmanr(s["log_mbh"], np.log10(s["lp_mpc"]))
    out["mbh_vs_length"] = {"n": int(len(s)), "spearman_rho": float(rho), "p": float(p)}
    out["giants_all"] = {"n_with_mbh": int(g["sigma_ok"].sum()),
                         "median_mbh_1e9": float(np.nanmedian(g.loc[g["sigma_ok"], "mbh"]) / 1e9),
                         "p10_p90_1e9": [float(v) for v in np.nanpercentile(g.loc[g["sigma_ok"], "mbh"], [10, 90]) / 1e9]}
    return out


def _boot_median_diff(a, b, n=2000, seed=7) -> float:
    rng = np.random.default_rng(seed)
    d = [np.median(rng.choice(a, a.size)) - np.median(rng.choice(b, b.size)) for _ in range(n)]
    return float(np.std(d))


def compare_stellar_masses() -> dict:
    """Stellar masses (value-added catalogue SED fits) of giants and controls,
    matched in redshift and luminosity."""
    g = giant_hosts()
    rg = catalogs.radio_galaxies()
    c = rg[(rg["L144"] >= config.RADIO_LOUD_W_HZ) & ~rg["is_giant_candidate"] & rg["logmass"].notna()]
    gm = g[np.isfinite(g["logmass_vac"]) & np.isfinite(g["L144"])]
    lc, zc, mc = np.log10(c["L144"].to_numpy()), c["z"].to_numpy(), c["logmass"].to_numpy()
    diffs = []
    for _, row in gm.iterrows():
        m = (np.abs(zc - row["z"]) < 0.02) & (np.abs(lc - np.log10(row["L144"])) < 0.15)
        if m.sum() >= 3:
            diffs.append(row["logmass_vac"] - np.median(mc[m]))
    diffs = np.array(diffs)
    rng = np.random.default_rng(5)
    boot = [np.median(rng.choice(diffs, diffs.size)) for _ in range(2000)] if diffs.size else [np.nan]
    return {"n_giants": int(diffs.size), "median_difference_dex": float(np.median(diffs)) if diffs.size else np.nan,
            "err_dex": float(np.std(boot)),
            "giants_median_logmass": float(np.median(gm["logmass_vac"])) if len(gm) else np.nan}
