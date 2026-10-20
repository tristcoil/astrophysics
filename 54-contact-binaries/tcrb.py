"""T CrB: a white dwarf feeding on a red giant, waiting for its third recorded nova.

Everything here comes from the AAVSO International Database (1866 to today),
fetched through the public VSX API:

* the two eruptions, 1866 and 1946, overlaid: peak, decline times t2 and t3,
  and the unique secondary maximum months later;
* the decade before each eruption: the "high state" (more accretion, a bluer
  and brighter system) and the pre-eruption dip of 1945, against 2015-2026;
* the red giant's tidal distortion: an ellipsoidal light curve with two
  minima per 227.6-day orbit;
* how often the next eruption has been forecast, and what the recurrence
  record alone says about when it may come.

Raw visual estimates include occasional false "eruption" reports (a visual
magnitude of 4.9 was entered on 2026 October 9, the same hour CCDs measured
V = 10.1), so every daily bin is cleaned against its own median.
"""

from __future__ import annotations

import gzip
import json

import numpy as np
import pandas as pd

import config as C
import literature as L

#: Bands kept from the AID and what they stand for.
BANDS = {"Vis.": "visual", "V": "Johnson V", "B": "Johnson B", "TG": "DSLR green (V-like)",
         "CV": "unfiltered, V zero point", "R": "Cousins R", "I": "Cousins I", "TB": "DSLR blue"}


#: Individual AID records we exclude, with the reason. Kept short and explicit.
SUSPECT_RECORDS = {
    49969718: ("J. Meek, visual 5.1 at JD 2431856.0 (1946 Feb 5): five days before an eruption that rose within "
               "hours (Peltier saw 9.9 on Feb 2); the brightness matches the decline ten days later, so most "
               "likely a date-entry error"),
}


def load_aavso(cache: bool = True, refresh: bool = False) -> pd.DataFrame:
    """All AAVSO observations of T CrB, de-duplicated by observation Id.

    The parsed table is cached in data/derived/; ``refresh`` rebuilds it from
    the raw JSON files (the analysis does, so the cache always matches).

    Columns: jd, mag, err, band, obscode, valflag, mtype. Observations flagged
    discrepant by AAVSO (ValFlag 'T') are dropped here; everything else is kept
    for the cleaning step.
    """
    path = C.DERIVED_DIR / "tcrb_aavso.parquet"
    if cache and path.exists() and not refresh:
        return pd.read_parquet(path)
    rows = []
    for p in sorted(C.AAVSO_DIR.glob("tcrb_*.json.gz")):
        with gzip.open(p) as fh:
            for o in json.load(fh):
                rows.append((o.get("Id"), o.get("JD"), o.get("Mag"), o.get("Uncertainty"), o.get("Band"),
                             o.get("ObsCode"), o.get("ValFlag"), o.get("MType"), o.get("ObsType")))
    df = pd.DataFrame(rows, columns=["id", "jd", "mag", "err", "band", "obscode", "valflag", "mtype", "obstype"])
    df = df.drop_duplicates("id")
    df = df[~pd.to_numeric(df["id"], errors="coerce").isin(list(SUSPECT_RECORDS))]
    df["jd"] = pd.to_numeric(df["jd"], errors="coerce")
    df["mag"] = pd.to_numeric(df["mag"], errors="coerce")
    df["err"] = pd.to_numeric(df["err"], errors="coerce")
    df = df[np.isfinite(df["jd"]) & np.isfinite(df["mag"]) & (df["valflag"] != "T")]
    df = df[df["band"].isin(BANDS)].sort_values("jd").reset_index(drop=True)
    df["year"] = 2000.0 + (df["jd"] - C.JD_J2000) / C.YEAR_D
    if cache:
        df.to_parquet(path)
    return df


def clean_daily(df: pd.DataFrame, band: str, max_dev: float = 1.0, min_observers: int = 3) -> pd.DataFrame:
    """Observations of one band with gross outliers flagged.

    Each estimate is compared with that calendar day's consensus: the median
    over observers of each observer's own median, so that one observer who
    submits several identical bad estimates (as on 2026 October 9) counts
    once. Days with fewer than ``min_observers`` observers are not tested.
    A same-day consensus follows a real eruption (T CrB rose 7 mag in a day
    in 1946), which a running median over several days would not.
    The rejected rows are kept, with the boolean column 'rejected', for plots.
    """
    d = df[df["band"] == band].copy()
    d["day"] = np.floor(d["jd"] - 0.5).astype(np.int64)
    per_obs = d.groupby(["day", "obscode"])["mag"].median().reset_index()
    cons = per_obs.groupby("day")["mag"].agg(["median", "size"])
    cons = cons[cons["size"] >= min_observers]["median"]
    d["ref"] = d["day"].map(cons)
    d["rejected"] = d["ref"].notna() & (np.abs(d["mag"] - d["ref"]) > max_dev)
    return d


def binned(d: pd.DataFrame, width_days: float = 1.0, min_n: int = 1) -> pd.DataFrame:
    """Median magnitude in bins of ``width_days`` (rejected points excluded)."""
    d = d[~d["rejected"]] if "rejected" in d else d
    b = np.floor(d["jd"] / width_days).astype(np.int64)
    g = d.groupby(b)
    out = pd.DataFrame({"jd": g["jd"].median(), "mag": g["mag"].median(), "n": g["mag"].size(),
                        "mad": g["mag"].agg(lambda x: 1.4826 * np.median(np.abs(x - np.median(x))))})
    out = out[out["n"] >= min_n].reset_index(drop=True)
    out["year"] = 2000.0 + (out["jd"] - C.JD_J2000) / C.YEAR_D
    return out


# ---------------------------------------------------------------------------
# Eruptions
# ---------------------------------------------------------------------------


def eruption_curve(df: pd.DataFrame, year: int, before: float = 30.0, after: float = 400.0) -> dict:
    """Light curve around one eruption, in days from Schaefer's (2023) time of peak."""
    t_pk = L.SCHAEFER2023["peak_jd"][year]
    out = {}
    for band in ("Vis.", "V", "B"):
        d = clean_daily(df[(df["jd"] > t_pk - before) & (df["jd"] < t_pk + after)], band, max_dev=1.5)
        if len(d) == 0:
            continue
        d = d.assign(days=d["jd"] - t_pk)
        out[band] = d
    return out


def decline_times(days, mag, peak_mag=None, drops=(2.0, 3.0)):
    """Times (days after peak) when a smoothed decline first falls ``drop`` mag
    below peak. ``peak_mag`` defaults to the brightest binned point."""
    o = np.argsort(days)
    days, mag = np.asarray(days)[o], np.asarray(mag)[o]
    # half-day medians smooth the visual scatter without shifting the decline
    b = np.floor(days * 2) / 2
    ub = np.unique(b)
    md = np.array([np.median(mag[b == x]) for x in ub])
    if peak_mag is None:
        peak_mag = md.min()
        t_peak = ub[np.argmin(md)]
    else:
        t_peak = 0.0  # the assumed peak (Schaefer's time of peak light)
    res = {}
    for drop in drops:
        after = ub > t_peak
        k = np.nonzero(after & (md >= peak_mag + drop))[0]
        if k.size == 0:
            res[drop] = np.nan
            continue
        j = k[0]
        # the crossing must be bracketed by data: a bin brighter than the threshold just before it
        if j == 0 or md[j - 1] >= peak_mag + drop or ub[j - 1] < t_peak:
            res[drop] = np.nan
            continue
        x0, x1, y0, y1 = ub[j - 1], ub[j], md[j - 1], md[j]
        res[drop] = float(x0 + (peak_mag + drop - y0) * (x1 - x0) / (y1 - y0)) if y1 != y0 else float(x1)
    return {"peak_mag": float(peak_mag), "t_peak": float(t_peak), "t_drop": res}


# ---------------------------------------------------------------------------
# The decade before: high state and dip
# ---------------------------------------------------------------------------


def seasonal_levels(df: pd.DataFrame, band: str, width_days: float = 30.0) -> pd.DataFrame:
    """Monthly medians of cleaned observations in one band."""
    return binned(clean_daily(df, band), width_days=width_days, min_n=5)


# ---------------------------------------------------------------------------
# Ellipsoidal modulation
# ---------------------------------------------------------------------------


def ellipsoidal_fold(df: pd.DataFrame, band: str = "V", jd_lo: float = 2450000.0, jd_hi: float = 2457000.0,
                     p_orb: float | None = None, t_ref: float | None = None, smooth_days: float = 400.0):
    """Fold detrended photometry on the orbit (Fekel et al. 2000 ephemeris).

    Each point has a running median over ``smooth_days`` (longer than the
    orbit) subtracted, which removes the slow accretion changes but keeps
    the orbital signal. Phase 0 is the maximum of the giant's radial velocity
    (Fekel et al.), so the giant is in front of the white dwarf at phase 0.75.
    """
    p_orb = p_orb or L.SCHAEFER2023["p_orb_fekel"]
    t_ref = t_ref or L.SCHAEFER2023["t_vmax_fekel"]
    d = clean_daily(df[(df["jd"] > jd_lo) & (df["jd"] < jd_hi)], band)
    d = d[~d["rejected"]]
    day = binned(d, 1.0)
    trend = np.array([np.median(day["mag"][np.abs(day["jd"] - t) < 0.5 * smooth_days]) for t in day["jd"]])
    resid = day["mag"].to_numpy() - trend
    phase = np.mod((day["jd"].to_numpy() - t_ref) / p_orb, 1.0)
    return phase, resid, day


def fourier_fit(phase, y, order=2):
    A = [np.ones_like(phase)]
    for k in range(1, order + 1):
        A += [np.cos(2 * np.pi * k * phase), np.sin(2 * np.pi * k * phase)]
    A = np.column_stack(A)
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return coef, A @ coef


# ---------------------------------------------------------------------------
# When? The recurrence record as an empirical hazard
# ---------------------------------------------------------------------------


def recurrence_hazard(intervals, last_eruption_year, now_year, horizon_years):
    """Probability of an eruption between now and now + horizon, given none so
    far, if the next interval is normal with the mean and sample standard
    deviation of the past (effective) intervals, as in Pei et al. (2026).

    An empirical indicator from three intervals, not physics. Exact normal
    integrals, no sampling.
    """
    from scipy.stats import norm
    iv = np.asarray(intervals, float)
    mu, sd = iv.mean(), iv.std(ddof=1)
    a = (now_year - last_eruption_year - mu) / sd
    b = (now_year + horizon_years - last_eruption_year - mu) / sd
    p = (norm.cdf(b) - norm.cdf(a)) / norm.sf(a)
    return {"mean_interval": float(mu), "sd_interval": float(sd), "p": float(p), "intervals": iv,
            "elapsed": float(now_year - last_eruption_year)}
