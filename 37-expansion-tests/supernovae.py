"""Assemble the DES supernova sample for the time-dilation test.

Two views of the same 1,635 supernovae:

* ``salt3_view()`` -- exactly the data White et al. (2024) used: the points of
  DES's SALT3 light-curve fits (the ``LCPLOT`` file), normalised by the fitted
  model's peak flux in each band. Those points were selected by SALT3's fit
  window of -15 to +45 *rest-frame* days, so this view already contains one
  (1 + z) assumption. It is used only to replicate White et al.

* ``photometry_view()`` -- every measurement from the raw DES photometry within
  a fixed *observer-frame* window around peak, so the selection cannot prefer
  any amount of time dilation. Our own global fit uses this view.

Both keep only likely Type Ia supernovae (probability > 0.5, as White et al.).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from config import BANDS, DERIVED_DIR, DES_FILTERS, FLUXERR_MAX, PROB_IA_MIN
import fetch_data

#: SNANA's photometry rejection mask for DES: bits 8-512 mark measurements
#: that the DES analyses never fit (the LCPLOT file contains none of them).
PHOTFLAG_REJECT = 8 | 16 | 32 | 64 | 128 | 256 | 512  # = 1016

#: DES fields: C3 and X3 are the deep fields (longer exposures, so fainter
#: supernovae are found and classified); the other eight are shallow.
DEEP_FIELDS = ("C3", "X3")


def supernova_table() -> pd.DataFrame:
    """One row per DES-SN5YR cosmology supernova, with selection columns.

    ``z`` is the heliocentric redshift: time dilation is set by the redshift we
    observe, including our own motion and the host's.
    """
    meta = fetch_data.load_metadata()
    des = meta[meta["IDSURVEY"] == 10].copy()
    des["CID"] = des["CID"].astype(np.int64)
    des["prob_ia"] = 1.0 - des["PROBCC_BEAMS"]
    head = fetch_data.load_des_head().set_index("SNID")
    des = des.join(head[["PTROBS_MIN", "PTROBS_MAX", "HOSTGAL_LOGMASS", "MWEBV"]],
                   on="CID", rsuffix="_head")
    sn = pd.DataFrame({
        "cid": des["CID"].to_numpy(),
        "z": des["zHEL"].to_numpy(float),
        "pkmjd": des["PKMJD"].to_numpy(float),
        "x1": des["x1"].to_numpy(float),
        "c": des["c"].to_numpy(float),
        "prob_ia": des["prob_ia"].to_numpy(float),
        "host_logmass": des["HOST_LOGMASS"].to_numpy(float),
        "ptr_min": des["PTROBS_MIN"].to_numpy(np.int64),
        "ptr_max": des["PTROBS_MAX"].to_numpy(np.int64),
    })
    sn["is_ia"] = sn["prob_ia"] > PROB_IA_MIN
    return sn.sort_values("cid").reset_index(drop=True)


def salt3_view(sn: pd.DataFrame | None = None):
    """White et al.'s data: LCPLOT points used by SALT3, and model peak fluxes.

    Returns ``(points, peaks)``. ``points`` has one row per measurement with
    ``t = MJD - PKMJD`` (observer days), ``flux`` and ``err`` divided by the
    SALT3 model's maximum in that band. ``peaks`` gives that maximum per
    supernova and band. Points with a raw flux error above 20 are dropped.
    """
    if sn is None:
        sn = supernova_table()
    lc = fetch_data.load_lcplot()
    lc = lc[lc["SNID"].isin(sn.loc[sn.is_ia, "cid"])]
    model = lc[lc["DATA_MODEL"] == 0]
    peaks = model.groupby(["SNID", "BAND"])["FLUXCAL"].max().rename("peak").reset_index()
    data = lc[(lc["DATA_MODEL"] == 1) & (lc["FLUXCAL_ERR"] <= FLUXERR_MAX)]
    data = data.merge(peaks, on=["SNID", "BAND"], how="inner")
    data = data.merge(sn[["cid", "z", "pkmjd"]], left_on="SNID", right_on="cid")
    points = pd.DataFrame({
        "cid": data["SNID"].to_numpy(np.int64),
        "band": data["BAND"].to_numpy(str),
        "z": data["z"].to_numpy(float),
        "t": (data["MJD"] - data["pkmjd"]).to_numpy(float),
        "flux": (data["FLUXCAL"] / data["peak"]).to_numpy(float),
        "err": (data["FLUXCAL_ERR"] / data["peak"]).to_numpy(float),
    })
    return points, peaks.rename(columns={"SNID": "cid", "BAND": "band"})


def photometry_view(sn: pd.DataFrame | None = None, *, window=(-40.0, 140.0),
                    max_err: float | None = None) -> pd.DataFrame:
    """Raw DES photometry of the likely Type Ia supernovae near peak.

    ``window`` is in observer-frame days relative to the SALT3 peak date, the
    same for every supernova. Fluxes stay in FLUXCAL units (zero point 27.5);
    the global fit gives each supernova and band its own amplitude. Cached.
    """
    if sn is None:
        sn = supernova_table()
    tag = f"phot_{window[0]:+.0f}_{window[1]:+.0f}_{max_err}"
    cache = DERIVED_DIR / f"{tag}.parquet"
    if cache.exists():
        return pd.read_parquet(cache)

    phot = fetch_data.load_des_phot()
    ia = sn[sn.is_ia]
    blocks = []
    for row in ia.itertuples():
        block = phot.iloc[row.ptr_min - 1: row.ptr_max]
        block = block[block["MJD"] > 0]  # drop SNANA's end-of-light-curve marker
        t = block["MJD"].to_numpy() - row.pkmjd
        keep = ((t >= window[0]) & (t <= window[1])
                & ((block["PHOTFLAG"].to_numpy() & PHOTFLAG_REJECT) == 0)
                & np.isin(block["BAND"].to_numpy(), BANDS)
                & (block["FLUXCALERR"].to_numpy() > 0))
        if max_err is not None:
            keep &= block["FLUXCALERR"].to_numpy() <= max_err
        b = block[keep]
        blocks.append(pd.DataFrame({
            "cid": row.cid, "band": b["BAND"].to_numpy(str), "field": b["FIELD"].to_numpy(str),
            "z": row.z, "t": t[keep], "flux": b["FLUXCAL"].to_numpy(float),
            "err": b["FLUXCALERR"].to_numpy(float),
        }))
    points = pd.concat(blocks, ignore_index=True)
    points.to_parquet(cache)
    return points


def rest_wavelength(band, z):
    """Central wavelength of a DES band in the supernova's frame (nm).

    This uses only the fact that wavelengths arrive stretched by (1 + z), which
    tired light also accepts; it carries no assumption about time dilation.
    """
    lam = np.vectorize(lambda b: DES_FILTERS[b][0])(np.asarray(band))
    return lam / (1.0 + np.asarray(z))
