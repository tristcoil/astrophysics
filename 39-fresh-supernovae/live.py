"""The live element: supernovae classified in the last few weeks, fitted
exactly like the DR2 sample and placed on the same Hubble diagram.

Differences from DR2 that matter, and how they are handled:
  * Photometry: ZTF's forced photometry shipped inside alert packets (via
    ALeRCE), on the survey's real calibration (each image's own zero point).
    DR2's fluxes are on a deliberately blinded zero point "near 30", so the
    two samples need not share an absolute magnitude scale.
  * Baseline: alert forced photometry has no multi-year baseline correction,
    so a constant per-band offset is measured on pre-explosion epochs and
    removed (lightcurves.subtract_baseline).
  * Redshifts: the host galaxy's when TNS lists one; otherwise the redshift
    the classifier reported, usually from the supernova's own spectral
    features and rounded to 2-3 decimals. Those get sigma_z = 0.005.
  * Light curves are young: many are near or before peak, when stretch is
    poorly constrained. The same quality cuts apply, so those simply do not
    pass yet; re-running fetch_live.py later lets them mature.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

import config
from lightcurves import alerce_forced_to_lightcurve, subtract_baseline
from saltfit import _guess_t0, fit_lightcurve

SIGMA_Z_TNS = 5e-3


def latest_snapshot() -> Path | None:
    snaps = sorted(p for p in config.LIVE_DIR.iterdir()
                   if p.is_dir() and (p / "objects.json").exists())
    return snaps[-1] if snaps else None


def load_snapshot(snap: Path) -> tuple[dict, list[dict]]:
    meta = json.loads((snap / "snapshot.json").read_text()) if (snap / "snapshot.json").exists() else {}
    objects = json.loads((snap / "objects.json").read_text())
    out = []
    for name, o in objects.items():
        if not o.get("oid"):
            continue
        path = snap / "alerce" / f"{o['oid']}_forced.json"
        rows = json.loads(path.read_text()) if path.exists() else None
        if not rows:
            continue
        o = dict(o)
        o["lc_raw"] = alerce_forced_to_lightcurve(rows)
        out.append(o)
    return meta, out


def fit_snapshot(snap: Path | None = None) -> tuple[dict, pd.DataFrame, dict]:
    """Fit every supernova in the snapshot. Returns (meta, table, light curves)."""
    from sample import mw_ebv
    from cosmology import z_hel_to_cmb

    snap = snap or latest_snapshot()
    if snap is None:
        return {}, pd.DataFrame(), {}
    meta, objs = load_snapshot(snap)
    meta["snapshot"] = snap.name
    rows, lcs = [], {}
    for o in objs:
        z_host, z_tns = o.get("host_redshift"), o.get("redshift")
        if z_host:
            z, z_err, z_group = z_host, 1e-4, "galaxy"
        elif z_tns:
            z, z_err, z_group = z_tns, SIGMA_Z_TNS, "classification spectrum"
        else:
            continue
        lc = o["lc_raw"]
        if lc["good"].sum() < 5:
            continue
        t_peak = _guess_t0(lc, z)
        lc, offsets = subtract_baseline(lc, t_peak, z)
        ebv = float(mw_ebv([o["ra"]], [o["dec"]])[0])
        fit = fit_lightcurve(lc, z, ebv, t_peak, name=o["tns_name"])
        fit.update({
            "tns_name": o["tns_name"], "oid": o["oid"], "ra": o["ra"], "dec": o["dec"],
            "redshift": z, "redshift_err": z_err, "z_group": z_group,
            "z_cmb": float(z_hel_to_cmb([z], [o["ra"]], [o["dec"]])[0]),
            "mwebv": ebv, "type": o.get("type"), "discovery_date": o.get("discovery_date"),
            "classifying_groups": o.get("classifying_groups"),
            "baseline_bands": "".join(sorted(b[-1] for b in offsets)),
            "n_epochs": int(lc["good"].sum()),
            "last_mjd": float(lc["mjd"].max()),
        })
        rows.append(fit)
        lcs[o["tns_name"]] = lc
    table = pd.DataFrame(rows)
    if not table.empty:
        table = table.set_index("tns_name", drop=False)
    return meta, table, lcs
