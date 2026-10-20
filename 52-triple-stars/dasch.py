"""DASCH light curves: Harvard plate photometry, 1885-1992.

Each light curve is one row per plate exposure on which the star was measured,
calibrated against APASS (``magcal_magdep``: B band for the blue-sensitive
emulsions that dominate the collection; DASCH's recommended magnitude). The
time stamp ``date_jd`` is the heliocentric JD at mid-exposure, with an
estimated uncertainty ``time_accuracy_days``.

Cleaning (documented in the README):
* daschlab's standard rejections (AFLAGS bits listed in config.py),
* plus bad plate quality, uncertain dates and multiple-exposure blends,
* DASCH's own ``reject_flag`` and missing calibrated magnitudes,
* time stamps worse than 10% of the eclipse half-width.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

import config as C


def _slug(name):
    return name.lower().replace(" ", "_").replace("-", "_")


def read_raw(name: str, label: str = "target") -> pd.DataFrame:
    path = C.RAW / "dasch" / _slug(name) / f"lc_{label}.json"
    blob = json.loads(path.read_text())
    rows = blob["rows"]
    hdr = rows[0].split(",")
    data = [line.split(",") for line in rows[1:]]
    df = pd.DataFrame(data, columns=hdr[:len(data[0])] if data else hdr)
    num = ["date_jd", "limiting_mag_local", "magcal_magdep", "magcal_local_rms", "magcal_magdep_rms",
           "time_accuracy_days", "aflags", "bflags", "reject_flag", "plate_center_dist_deg", "fwhm_pix",
           "blended_mag", "extinction", "magcal_local_error"]
    for c in num:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df["source"] = blob["source"]["ref_text"]
    return df


def clean(df: pd.DataFrame, eclipse_hw_days: float | None = None) -> tuple[pd.DataFrame, dict]:
    """Apply the cuts; return the kept rows and a bookkeeping dict."""
    n0 = len(df)
    log = {"rows": n0}
    m = df["magcal_magdep"].notna() & (df["magcal_magdep"] < 50)
    log["no_magnitude"] = int((~m).sum())
    af = df["aflags"].fillna(0).astype(np.int64)
    bad_std = (af & C.DASCH_BAD_AFLAGS) != 0
    bad_extra = (af & C.DASCH_EXTRA_BAD_AFLAGS) != 0
    log["standard_aflags"] = int((m & bad_std).sum())
    log["extra_aflags"] = int((m & ~bad_std & bad_extra).sum())
    keep = m & ~bad_std & ~bad_extra
    rej = df["reject_flag"].fillna(0) != 0
    log["reject_flag"] = int((keep & rej).sum())
    keep &= ~rej
    if eclipse_hw_days is not None:
        tmax = C.DASCH_TIME_ACCURACY_FRACTION * eclipse_hw_days
        badt = df["time_accuracy_days"].fillna(9) > tmax
        log["time_accuracy"] = int((keep & badt).sum())
        keep &= ~badt
    out = df[keep].copy()
    out["year"] = C.jd_to_year(out["date_jd"])
    out["err"] = out["magcal_local_rms"].where(out["magcal_local_rms"] > 0, np.nan)
    log["kept"] = int(len(out))
    return out.sort_values("date_jd").reset_index(drop=True), log


def load(name: str, label: str = "target") -> tuple[pd.DataFrame, dict]:
    hw = C.SYSTEMS[name].get("eclipse_hw_days")
    return clean(read_raw(name, label), hw)


def comparison_scatter(name: str) -> dict:
    """Robust scatter of the comparison stars, per plate series: the noise the
    eclipse fit should expect for a constant star of similar brightness."""
    out = {}
    folder = C.RAW / "dasch" / _slug(name)
    for path in sorted(folder.glob("lc_comp*.json")):
        label = path.stem[3:]
        d, _ = load(name, label)
        med = np.median(d["magcal_magdep"])
        mad = 1.4826 * np.median(np.abs(d["magcal_magdep"] - med))
        out[label] = dict(n=len(d), median=float(med), robust_sigma=float(mad),
                          median_err=float(np.nanmedian(d["err"])))
    return out
