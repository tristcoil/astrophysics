"""The Bright Transient Survey catalogue: parse, decode, and put distances on it.

The Sample Explorer's CSV is text: limits are written as ">21.0", missing
values as "-", and the survey's quality and purity cuts are packed into digit
strings. ``load()`` turns it into typed columns, decodes every flag, maps the
TNS classifications to families, and recomputes absolute magnitudes two ways:

* **the BTS convention** (``M_bts``): heliocentric redshift, H0 = 70, Om = 0.3,
  K = 2.5 log10(1+z). It reproduces the explorer's own M_abs column (checked
  in ``check_bts_convention``), so the replication of Perley et al. uses it.
* **our distances** (``M_peak``): the same, except that redshifts are moved to
  the CMB frame and hosts closer than cz = 4,000 km/s get a redshift-
  independent Cosmicflows-4 distance (distances.py), because there a peculiar
  velocity of a few hundred km/s is a large fraction of the redshift.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

import config as C
from distances import cosmology, distance_modulus, cmb_frame_redshift, cf4_distances
from fetch_data import BTS_ALL, BTS_SDSS, read_bts_csv

NUMERIC = ("peakt", "peakmag", "peakabs", "duration", "rise", "fade", "redshift", "hostabs",
           "hostcol", "b", "A_V")


def _num(s: pd.Series) -> tuple[pd.Series, pd.Series]:
    """'>21.03' -> (21.03, True); '-' -> (NaN, False)."""
    s = s.astype(str).str.strip()
    is_limit = s.str.startswith(">") | s.str.startswith("<")
    val = pd.to_numeric(s.str.lstrip("<>"), errors="coerce")
    return val, is_limit


def load(with_distances: bool = True) -> pd.DataFrame:
    raw = read_bts_csv(BTS_ALL)
    df = pd.DataFrame({"ztf_id": raw["ZTFID"], "iau_id": raw["IAUID"],
                       "tns_type": raw["type"].replace("", "-")})
    for col in NUMERIC:
        val, lim = _num(raw[col])
        df[col] = val
        if col in ("duration", "rise", "fade"):
            df[col + "_is_limit"] = lim
    df.loc[df["redshift"] <= 0, "redshift"] = np.nan
    df["peak_filter"] = raw["peakfilt"]
    c = SkyCoord(raw["RA"].values, raw["Dec"].values, unit=(u.hourangle, u.deg))
    df["ra"], df["dec"] = c.ra.deg, c.dec.deg
    df["peak_jd"] = df["peakt"] + 2458000.0
    df["peak_date"] = pd.to_datetime(Time(df["peak_jd"].values, format="jd").isot)

    # --- the survey's own cuts, decoded bit by bit (config.QUALCODE_BITS) ---
    for code, names in (("qualcode", C.QUALCODE_BITS), ("purecode", C.PURECODE_BITS)):
        s = raw[code].astype(str).str.strip('"')
        for i, name in enumerate(names):
            df["q_" + name] = s.str[i] == "1"
    df["quality"] = raw["qualcode"].astype(str).str.strip('"') == "1" * len(C.QUALCODE_BITS)
    df["purity"] = df["q_sn_like_lightcurve"] | df["q_galaxy_match"]
    df["in_sdss_field"] = df["ztf_id"].isin(set(read_bts_csv(BTS_SDSS)["ZTFID"]))

    # --- classes ---
    df["cls"] = df["tns_type"].map(C.CLASS_MAP)
    df["family"] = df["cls"].map(C.FAMILY)
    df["classified"] = df["tns_type"] != "-"
    df["is_cc"] = df["cls"].isin(C.CORE_COLLAPSE_CLASSES)

    # --- absolute magnitudes in the BTS convention ---
    a_band = np.where(df["peak_filter"] == "g", C.A_G_PER_AV, C.A_R_PER_AV) * df["A_V"]
    df["A_peak"] = a_band
    z = df["redshift"].to_numpy()
    k = 2.5 * np.log10(1 + z)
    df["mu_bts"] = distance_modulus(z)
    df["M_bts"] = df["peakmag"] - df["mu_bts"] - df["A_peak"] + k

    if with_distances:
        add_distances(df)
    return df


def add_distances(df: pd.DataFrame) -> None:
    """Our adopted distance moduli (h70 scale) and absolute magnitudes."""
    z_helio = df["redshift"].to_numpy()
    z_cmb = cmb_frame_redshift(z_helio, df["ra"].to_numpy(), df["dec"].to_numpy())
    df["z_cmb"] = z_cmb
    mu = distance_modulus(z_cmb)
    df["mu_source"] = np.where(np.isfinite(mu), "redshift (CMB frame)", "none")

    cf4 = cf4_distances(df)
    have = cf4["mu_h70"].notna()
    df["cf4_pgc"] = cf4["pgc"]
    df["cf4_sep_kpc"] = cf4["sep_kpc"]
    df["mu_cf4"] = cf4["mu_h70"]
    mu = np.where(have, cf4["mu_h70"], mu)
    df.loc[have, "mu_source"] = "Cosmicflows-4 group"
    df["mu"] = mu
    k = 2.5 * np.log10(1 + df["z_cmb"])
    df["M_peak"] = df["peakmag"] - df["mu"] - df["A_peak"] + k


def check_bts_convention(df: pd.DataFrame) -> dict:
    """How well our BTS-convention magnitudes reproduce the explorer's column."""
    ok = df["peakabs"].notna() & df["M_bts"].notna() & (df["redshift"] > 0.002)
    d = (df.loc[ok, "M_bts"] - df.loc[ok, "peakabs"]).to_numpy()
    return {"n": int(ok.sum()), "median": float(np.median(d)),
            "p99_abs": float(np.percentile(np.abs(d), 99)),
            "frac_within_0p01": float(np.mean(np.abs(d) < 0.01))}


def statistical_sample(df: pd.DataFrame, start=C.STAT_PEAK_START, end=C.STAT_PEAK_END,
                       m_lim=C.M_LIM) -> pd.Series:
    """Perley et al.'s statistical sample: quality + purity, peak in the window,
    brighter than the classification limit (boolean mask over all candidates)."""
    t = df["peak_date"]
    return (df["quality"] & df["purity"] & (t >= pd.Timestamp(start))
            & (t <= pd.Timestamp(end) + pd.Timedelta(days=1)) & (df["peakmag"] <= m_lim))


def classification_completeness(df: pd.DataFrame, mask: pd.Series,
                                bins=np.arange(12.0, 19.01, 0.25)) -> pd.DataFrame:
    """Fraction of statistical-sample candidates that have a classification,
    per bin of peak magnitude. (Unclassified candidates are assumed to have the
    same type mix as classified ones at the same magnitude.)"""
    sub = df[mask]
    idx = np.digitize(sub["peakmag"], bins) - 1
    rows = []
    for i in range(len(bins) - 1):
        s = sub[idx == i]
        if len(s):
            rows.append({"m_lo": bins[i], "m_hi": bins[i + 1], "n": len(s),
                         "n_classified": int(s["classified"].sum()),
                         "frac": float(s["classified"].mean())})
    return pd.DataFrame(rows)


def fcl_model(m: np.ndarray, table: pd.DataFrame) -> np.ndarray:
    """Classification efficiency at peak magnitude m: the measured fraction,
    smoothed by a monotonic (non-increasing) fit, interpolated in m."""
    mid = 0.5 * (table["m_lo"] + table["m_hi"]).to_numpy()
    f = table["frac"].to_numpy()
    w = table["n"].to_numpy()
    # pool-adjacent-violators for a non-increasing sequence
    vals, wts, cnt = list(f), list(w.astype(float)), [1] * len(f)
    i = 0
    while i < len(vals) - 1:
        if vals[i] < vals[i + 1]:
            tot = wts[i] + wts[i + 1]
            vals[i] = (vals[i] * wts[i] + vals[i + 1] * wts[i + 1]) / tot
            wts[i] = tot
            cnt[i] += cnt[i + 1]
            del vals[i + 1], wts[i + 1], cnt[i + 1]
            i = max(i - 1, 0)
        else:
            i += 1
    smooth = np.repeat(vals, cnt)
    return np.interp(m, mid, smooth, left=smooth[0], right=smooth[-1])


def fcl_perley(m: np.ndarray) -> np.ndarray:
    """Perley et al.'s assumed efficiency: 1 below 17.2, linearly to 0.9 at 18.5."""
    m = np.asarray(m, float)
    slope = (1.0 - C.PERLEY_FCL_AT_LIMIT) / (C.M_LIM - C.PERLEY_FCL_BRIGHT_MAG)
    return np.clip(1.0 - slope * (m - C.PERLEY_FCL_BRIGHT_MAG), C.PERLEY_FCL_AT_LIMIT, 1.0)
