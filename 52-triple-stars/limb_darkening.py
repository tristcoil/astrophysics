"""Quadratic limb-darkening coefficients (u1, u2) from Claret's ATLAS tables:
TESS band from Claret (2017, A&A 600, A30, Table 25) and Johnson B (the
photographic plates, calibrated to APASS B) from Claret (2000, A&A 363, 1081).
Solar metallicity, microturbulence 2 km/s; bilinear interpolation in
(log g, Teff)."""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator

import config as C


@lru_cache(maxsize=None)
def _tess():
    path = C.RAW / "vizier" / "J_AA_600_A30" / "table25.dat"
    df = pd.read_fwf(path, colspecs=[(0, 5), (6, 12), (13, 17), (19, 22), (23, 31), (33, 41)],
                     names=["logg", "teff", "z", "xi", "a", "b"], header=None)
    return df[(df["z"] == 0.0) & (df["xi"] == 2.0)]


@lru_cache(maxsize=None)
def _ubvri(band="B"):
    path = C.RAW / "vizier" / "J_AA_363_1081" / "atlasco.dat"
    col = {"B": (68, 75), "V": (77, 84)}[band]
    cols = [(0, 1), (3, 4), (4, 9), (10, 16), (17, 21), col]
    df = pd.read_fwf(path, colspecs=cols, names=["coef", "vt", "logg", "teff", "z", "m"], header=None)
    df = df[(df["z"] == 0.0) & (df["vt"] == 2)]
    a = df[df["coef"] == "a"][["logg", "teff", "m"]].rename(columns={"m": "a"})
    b = df[df["coef"] == "b"][["logg", "teff", "m"]].rename(columns={"m": "b"})
    return a.merge(b, on=["logg", "teff"])


def quadratic(teff: float, logg: float, band: str) -> tuple[float, float]:
    df = _tess() if band == "TESS" else _ubvri(band)
    pts = np.column_stack([df["logg"], np.log10(df["teff"])])
    out = []
    for c in ("a", "b"):
        f = LinearNDInterpolator(pts, df[c].to_numpy())
        out.append(float(f(logg, np.log10(teff))))
    return out[0], out[1]
