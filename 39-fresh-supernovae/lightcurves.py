"""Light curves: read ZTF SN Ia DR2 files and ALeRCE forced photometry into one
common format that the SALT fitter takes.

Common format (a pandas DataFrame, one row per observation):
    mjd     observation time [MJD]
    band    sncosmo passband name: ztfg, ztfr, ztfi
    flux    flux in units where `zp` is the AB magnitude of a flux of 1
    fluxerr its 1-sigma error
    zp      AB zero point of `flux` (30 for DR2: "blinded, near 30")
    good    True if no bad-quality flag is set

Why fluxes and not magnitudes: forced photometry measures the flux at the
supernova's position on every image, including before the explosion and long
after it has faded. Those measurements scatter around zero and can be negative.
A magnitude cannot represent them, but they still carry information (they pin
down when the supernova was *not* there), so the fit uses fluxes throughout.
"""

from __future__ import annotations

import gzip
import io
import re
from pathlib import Path

import numpy as np
import pandas as pd

import config

#: DR2 photometry flags (Rigault et al. 2025, Table 3; ztfcosmo). Bits that
#: mark a point as bad, per the release: 1 flux_err == 0, 2 chi2/dof > 3 in the
#: PSF fit, 4 cloudy, 8 bad pixels/processing (infobits), 16 shallow image
#: (5-sigma limit brighter than 19.3 mag). Other bits are informational.
DR2_BAD_FLAGS = 1 | 2 | 4 | 8 | 16

#: ZTF filter id (ALeRCE `fid`) -> sncosmo passband.
FID_TO_BAND = {1: "ztfg", 2: "ztfr", 3: "ztfi"}


def dr2_lightcurve_path(name: str) -> Path:
    return config.ZTF_LC_DIR / f"{name}_lc.csv.gz"


def read_dr2_lightcurve(name: str) -> tuple[pd.DataFrame, dict]:
    """One DR2 light curve in the common format, plus its header values."""
    with gzip.open(dr2_lightcurve_path(name), "rt") as f:
        text = f.read()
    header = _parse_dr2_header(text)
    raw = pd.read_csv(io.StringIO(text), sep=r"\s+", comment="#")
    lc = pd.DataFrame({
        "mjd": raw["mjd"].to_numpy(float),
        "band": raw["filter"].to_numpy(str),
        "flux": raw["flux"].to_numpy(float),
        "fluxerr": raw["flux_err"].to_numpy(float),
        "zp": raw["ZP"].to_numpy(float),
        "good": (raw["flag"].to_numpy(int) & DR2_BAD_FLAGS) == 0,
        "flag": raw["flag"].to_numpy(int),
    })
    # An error of zero is unphysical (flag 1 catches it, but be safe).
    lc.loc[~(lc["fluxerr"] > 0), "good"] = False
    return lc, header


def _parse_dr2_header(text: str) -> dict:
    """The few header values we use: the initial peak-time guess and the
    photometry version. Everything else in the header is kept as raw text."""
    header = {"raw": "\n".join(l for l in text.splitlines() if l.startswith("#"))}
    m = re.search(r"#\s*mjd_max_init:\s*([0-9.]+)", text)
    header["mjd_max_init"] = float(m.group(1)) if m else np.nan
    m = re.search(r"#\s*phot_source:\s*(.+)", text)
    header["phot_source"] = m.group(1).strip() if m else ""
    return header


def alerce_forced_to_lightcurve(rows: list[dict], zp: float = 30.0) -> pd.DataFrame:
    """ALeRCE's forced photometry (from ZTF alert packets) in the common format.

    ALeRCE stores each forced measurement as a magnitude of |flux| with a sign
    flag `isdiffpos` (+1/-1) and the image's own zero point `magzpsci`; we turn
    that back into a signed flux on a common zero point. Points whose
    processing status (`procstatus`) is not 0 are marked bad, as are those
    without a zero point.
    """
    out = []
    for r in rows:
        extra = r.get("extra_fields") or {}
        mag, emag = r.get("mag"), r.get("e_mag")
        band = FID_TO_BAND.get(r.get("fid"))
        if band is None or mag is None or emag is None:
            continue
        sign = 1.0 if (r.get("isdiffpos") or 1) > 0 else -1.0
        flux = sign * 10 ** (-0.4 * (mag - zp))
        fluxerr = abs(flux) * emag * np.log(10) / 2.5
        procstatus = str(extra.get("procstatus", "0")).strip()
        good = (procstatus in ("0", "")) and extra.get("magzpsci") is not None \
            and np.isfinite(flux) and fluxerr > 0
        out.append({"mjd": float(r["mjd"]), "band": band, "flux": flux,
                    "fluxerr": fluxerr, "zp": zp, "good": bool(good), "flag": 0})
    lc = pd.DataFrame(out, columns=["mjd", "band", "flux", "fluxerr", "zp", "good", "flag"])
    return lc.sort_values("mjd").drop_duplicates(["mjd", "band"]).reset_index(drop=True)


def subtract_baseline(lc: pd.DataFrame, t_peak: float, z: float,
                      before: float = -25.0, min_points: int = 3) -> tuple[pd.DataFrame, dict]:
    """Remove a constant flux offset per band, measured before the explosion.

    Difference images can leave a small constant residual at the supernova's
    position (an imperfect reference image). DR2 light curves are already
    corrected using years of baseline; alert forced photometry is not, so we
    take the weighted mean of points more than 25 rest-frame days before peak,
    where a Type Ia supernova has not yet exploded. Bands without enough such
    points are left alone.
    """
    lc = lc.copy()
    offsets = {}
    phase = (lc["mjd"] - t_peak) / (1 + z)
    for band in lc["band"].unique():
        pre = lc["good"] & (lc["band"] == band) & (phase < before)
        if pre.sum() >= min_points:
            w = 1 / lc.loc[pre, "fluxerr"] ** 2
            off = float(np.sum(w * lc.loc[pre, "flux"]) / np.sum(w))
            err = float(1 / np.sqrt(np.sum(w)))
            lc.loc[lc["band"] == band, "flux"] -= off
            offsets[band] = (off, err, int(pre.sum()))
    return lc, offsets


def detected_phases(lc: pd.DataFrame, t0: float, z: float,
                    window=(-10.0, 40.0), snr_min: float = 5.0) -> dict:
    """Light-curve coverage in the DR2 sense (Rigault et al. 2025, Sect. 3).

    Counts "phases": distinct (night, band) pairs with a detection at
    S/N >= 5 inside the rest-frame window. Several detections in the same band
    on the same night count once.
    """
    phase = (lc["mjd"] - t0) / (1 + z)
    det = lc["good"] & (lc["flux"] / lc["fluxerr"] >= snr_min) & \
        (phase >= window[0]) & (phase <= window[1])
    d = lc.loc[det].assign(night=np.floor(lc.loc[det, "mjd"]).astype(int),
                           phase=phase[det])
    unique = d.drop_duplicates(["night", "band"])
    return {
        "n_phases": int(len(unique)),
        "n_before": int((unique["phase"] < 0).sum()),
        "n_after": int((unique["phase"] > 0).sum()),
        "n_bands": int(unique["band"].nunique()),
        "bands": "".join(sorted(b[-1] for b in unique["band"].unique())),
    }


def good_sampling(cov: dict) -> bool:
    """DR2's "good light curve sampling" cut: at least 7 detected phases,
    2 before and 2 after peak, in at least 2 bands."""
    return cov["n_phases"] >= 7 and cov["n_before"] >= 2 and cov["n_after"] >= 2 \
        and cov["n_bands"] >= 2
