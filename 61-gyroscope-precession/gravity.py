"""The Earth's static gravity field, read from ICGEM ``.gfc`` files.

Coefficients are fully normalised Stokes coefficients C_nm, S_nm. Models with
time-variable terms (GOCO06s: ``gfct`` + ``trnd`` per year + annual ``acos``,
``asin`` about a reference epoch; EIGEN-6S4 v2: piecewise ``gfct``/``trnd`` with
validity intervals) are evaluated at the requested epoch.
"""

from __future__ import annotations

import datetime as dt
import functools
from dataclasses import dataclass

import numpy as np

import config as C


@dataclass
class Field:
    name: str
    gm: float
    radius: float
    nmax: int
    tide_system: str
    cbar: np.ndarray      # (nmax+1, nmax+1), fully normalised
    sbar: np.ndarray


def _float(s: str) -> float:
    return float(s.replace("D", "E").replace("d", "e"))


def _yyyymmdd_to_mjd(s: str) -> float:
    s = s.strip()
    day = dt.datetime.strptime(s[:8], "%Y%m%d")
    frac = float(s[8:]) if len(s) > 8 else 0.0
    return (day - dt.datetime(1858, 11, 17)).days + frac


@functools.lru_cache(maxsize=8)
def _parse(name: str, nmax: int):
    path = C.GRAVITY_DIR / f"{name.replace(' ', '_')}.gfc"
    head = {}
    static = {}
    tv = []           # (kind, n, m, c, s, t0, t1 or period)
    in_head = True
    with open(path, errors="replace") as fh:
        for line in fh:
            if in_head:
                parts = line.split()
                if line.startswith("end_of_head"):
                    in_head = False
                elif len(parts) >= 2:
                    head.setdefault(parts[0], parts[1])
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            kind = parts[0]
            n, m = int(parts[1]), int(parts[2])
            if n > nmax or m > nmax:
                continue
            c, s = _float(parts[3]), _float(parts[4])
            if kind == "gfc":
                static[(n, m)] = (c, s)
            elif kind in ("gfct", "trnd", "acos", "asin"):
                extra = parts[7:]
                tv.append((kind, n, m, c, s, extra))
    return head, static, tv


def load(name: str = C.GRAVITY_MODEL, nmax: int = C.GRAVITY_DEGREE, mjd: float | None = None) -> Field:
    """The field truncated at degree nmax, with time-variable terms evaluated at ``mjd``."""
    head, static, tv = _parse(name, nmax)
    cbar = np.zeros((nmax + 1, nmax + 1))
    sbar = np.zeros((nmax + 1, nmax + 1))
    for (n, m), (c, s) in static.items():
        cbar[n, m], sbar[n, m] = c, s
    if tv:
        if mjd is None:
            raise ValueError(f"{name} has time-variable terms: give an epoch")
        # GOCO06s style: gfct with one reference epoch, trnd per year, acos/asin with a period
        # EIGEN-6S4 style: gfct/trnd with validity intervals [t0, t1)
        for kind, n, m, c, s, extra in tv:
            if kind == "gfct":
                if len(extra) >= 2:          # interval form
                    t0, t1 = _yyyymmdd_to_mjd(extra[0]), _yyyymmdd_to_mjd(extra[1])
                    if not (t0 <= mjd < t1):
                        continue
                cbar[n, m] += c
                sbar[n, m] += s
        ref_goco = 55_197.0             # 2010-01-01, GOCO06s reference epoch
        for kind, n, m, c, s, extra in tv:
            if kind == "trnd":
                if len(extra) >= 2:
                    t0, t1 = _yyyymmdd_to_mjd(extra[0]), _yyyymmdd_to_mjd(extra[1])
                    if not (t0 <= mjd < t1):
                        continue
                    dt_yr = (mjd - t0) / 365.25
                else:
                    dt_yr = (mjd - ref_goco) / 365.25
                cbar[n, m] += c * dt_yr
                sbar[n, m] += s * dt_yr
            elif kind in ("acos", "asin"):
                if len(extra) >= 3:          # EIGEN-6S4: t0 t1 period
                    t0, t1 = _yyyymmdd_to_mjd(extra[0]), _yyyymmdd_to_mjd(extra[1])
                    if not (t0 <= mjd < t1):
                        continue
                    period = float(extra[2])
                    phase = 2 * np.pi * (mjd - t0) / (period * 365.25)
                else:
                    period = float(extra[0]) if extra else 1.0
                    phase = 2 * np.pi * (mjd - ref_goco) / (period * 365.25)
                f = np.cos(phase) if kind == "acos" else np.sin(phase)
                cbar[n, m] += c * f
                sbar[n, m] += s * f
    return Field(name=name, gm=_float(head["earth_gravity_constant"]),
                 radius=_float(head["radius"]), nmax=nmax,
                 tide_system=head.get("tide_system", "unknown"), cbar=cbar, sbar=sbar)


def normalisation(nmax: int) -> np.ndarray:
    """N_nm with C_nm(unnormalised) = N_nm * Cbar_nm."""
    from math import factorial
    N = np.zeros((nmax + 2, nmax + 2))
    for n in range(nmax + 2):
        for m in range(n + 1):
            N[n, m] = np.sqrt((2 - (m == 0)) * (2 * n + 1) * factorial(n - m) / factorial(n + m))
    return N


def legendre_bar(nmax: int, sinphi: float) -> np.ndarray:
    """Fully normalised associated Legendre functions Pbar_nm(sin phi), n <= nmax."""
    cosphi = np.sqrt(max(0.0, 1.0 - sinphi * sinphi))
    P = np.zeros((nmax + 1, nmax + 1))
    P[0, 0] = 1.0
    if nmax >= 1:
        P[1, 1] = np.sqrt(3.0) * cosphi
    for m in range(2, nmax + 1):
        P[m, m] = np.sqrt((2 * m + 1) / (2 * m)) * cosphi * P[m - 1, m - 1]
    for m in range(0, nmax):
        P[m + 1, m] = np.sqrt(2 * m + 3) * sinphi * P[m, m]
    for m in range(0, nmax + 1):
        for n in range(m + 2, nmax + 1):
            a = np.sqrt((4 * n * n - 1) / (n * n - m * m))
            b = np.sqrt(((n - 1) ** 2 - m * m) / (4 * (n - 1) ** 2 - 1))
            P[n, m] = a * (sinphi * P[n - 1, m] - b * P[n - 2, m])
    return P
