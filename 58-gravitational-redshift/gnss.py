"""GNSS orbits and clocks: read SP3 files, put the orbits in an inertial frame, and
compute what general relativity says the satellite clocks should do.

The physics in one line: a clock at distance r from the Earth's centre, moving at
speed v in the geocentric inertial frame (GCRS), ticks at

    d(tau)/dt = 1 - [U(r) + v^2/2] / c^2,

relative to coordinate time t (TCG), where U is the (positive) Newtonian
potential: the Earth's monopole, its flattening (J2) and the tides of the Moon
and Sun. Ground clocks on the geoid tick at 1 - W0/c^2, so a satellite clock
gains (W0 - U - v^2/2)/c^2 on them per second.

GNSS products follow the IGS convention: the satellite clocks are published
*without* the once-per-orbit term a Keplerian orbit adds,

    dt_rel = -2 (r . v) / c^2       (= -(2/c^2) sqrt(GM a) e sin E),

which every receiver adds back. ``r . v`` is the same in the rotating (ITRF)
and inertial frames because omega x r is perpendicular to r.
"""

from __future__ import annotations

import datetime as dt
import functools
import gzip
from dataclasses import dataclass, field
from pathlib import Path

import erfa
import numpy as np

import config as C

# ---------------------------------------------------------------------------
# Reading SP3
# ---------------------------------------------------------------------------

BAD_CLOCK = 999_999.0   # SP3 "no value" for clocks (us)


@dataclass
class SP3Day:
    """One day of SP3 records. Epochs in GPS seconds since 1980-01-06 00:00."""
    path: Path
    series: str
    gps_s: np.ndarray                       # (N,)
    pos: dict = field(default_factory=dict)  # prn -> (N, 3) metres, ITRF; NaN if missing
    clk: dict = field(default_factory=dict)  # prn -> (N,) seconds; NaN if missing
    flags: dict = field(default_factory=dict)  # prn -> list of str (SP3c columns 75-80)

    @property
    def date(self) -> dt.date:
        return C.GPS_EPOCH + dt.timedelta(seconds=float(self.gps_s[0]))


def _epoch_to_gps_s(line: str) -> float:
    y, mo, d, h, mi = (int(line[3:7]), int(line[8:10]), int(line[11:13]),
                       int(line[14:16]), int(line[17:19]))
    sec = float(line[20:31])
    days = (dt.date(y, mo, d) - C.GPS_EPOCH).days
    return days * 86_400.0 + h * 3600.0 + mi * 60.0 + sec


def read_sp3(path: Path) -> SP3Day:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        lines = fh.read().splitlines()
    head = lines[0]
    if head[2] not in "PV" or not head.startswith("#"):
        raise ValueError(f"{path}: not an SP3 file")
    time_sys = next((l[9:12] for l in lines if l.startswith("%c")), "GPS")
    if time_sys.strip() not in ("GPS", "ccc", ""):
        raise ValueError(f"{path}: time system {time_sys!r}, expected GPS")
    epochs: list[float] = []
    recs: dict[str, dict[int, tuple]] = {}
    for line in lines:
        if line.startswith("*"):
            epochs.append(_epoch_to_gps_s(line))
        elif line.startswith("P") and epochs:
            prn = line[1:4]
            x, y, z = float(line[4:18]), float(line[18:32]), float(line[32:46])
            c = float(line[46:60]) if len(line) >= 60 and line[46:60].strip() else BAD_CLOCK
            recs.setdefault(prn, {})[len(epochs) - 1] = (x, y, z, c, line[74:80] if len(line) > 74 else "")
    n = len(epochs)
    day = SP3Day(path=Path(path), series=Path(path).name.split("_")[0], gps_s=np.array(epochs))
    for prn, by_epoch in recs.items():
        pos = np.full((n, 3), np.nan)
        clk = np.full(n, np.nan)
        flags = [""] * n
        for i, (x, y, z, c, fl) in by_epoch.items():
            if not (x == 0.0 and y == 0.0 and z == 0.0):
                pos[i] = (x * 1e3, y * 1e3, z * 1e3)
            if abs(c) < BAD_CLOCK:
                clk[i] = c * 1e-6
            flags[i] = fl
        day.pos[prn], day.clk[prn], day.flags[prn] = pos, clk, flags
    return day


def sp3_files(series_dir: Path = C.SP3_DIR) -> list[Path]:
    """Daily filtered SP3 files, one per day (MGNFIN preferred over OPSFIN)."""
    by_day: dict[str, Path] = {}
    for p in sorted(series_dir.glob("*/*.sp3.gz")):
        series, ydoy = p.name.split(".")[0].split("_")
        if ydoy not in by_day or C.ESA_SERIES.index(series) < C.ESA_SERIES.index(
                by_day[ydoy].name.split("_")[0]):
            by_day[ydoy] = p
    return [by_day[k] for k in sorted(by_day)]


# ---------------------------------------------------------------------------
# Time scales and Earth orientation
# ---------------------------------------------------------------------------

JD_GPS_EPOCH = 2_444_244.5      # 1980-01-06 00:00


@functools.lru_cache(maxsize=1)
def _eop_table():
    from astropy.utils import iers
    t = iers.IERS_A.open(str(C.EOP_DIR / "finals2000A.all"))
    mjd = np.asarray(t["MJD"].value, float)
    xp = np.asarray(t["PM_x"].value, float)        # arcsec
    yp = np.asarray(t["PM_y"].value, float)
    dut1 = np.asarray(t["UT1_UTC"].value, float)   # s
    good = np.isfinite(xp) & np.isfinite(yp) & np.isfinite(dut1)
    return mjd[good], xp[good], yp[good], dut1[good]


def _tai_utc(mjd_utc: np.ndarray) -> np.ndarray:
    """TAI-UTC (s) from the leap-second table in config (GPS-UTC + 19 s)."""
    out = np.zeros_like(np.asarray(mjd_utc, float))
    for start, gps_utc in C.LEAP_SECONDS_GPS:
        mjd0 = (start - dt.date(1858, 11, 17)).days
        out = np.where(mjd_utc >= mjd0, gps_utc + 19.0, out)
    if np.any(out == 0):
        raise ValueError("date before the leap-second table")
    return out


def eop(mjd_tai: np.ndarray):
    """Polar motion (rad) and UT1-TAI (s) interpolated from IERS finals2000A.

    The table gives UT1-UTC, which jumps by 1 s at leap seconds; UT1-TAI is
    smooth. TAI = GPS + 19 s exactly, so no leap-second bookkeeping is needed
    for the satellite epochs themselves.
    """
    mjd, xp, yp, dut1 = _eop_table()
    sel = mjd >= 57_204          # 2015-07-01, where our leap-second table starts
    mjd, xp, yp, dut1 = mjd[sel], xp[sel], yp[sel], dut1[sel]
    ut1_tai = dut1 - _tai_utc(mjd)
    as2rad = np.pi / 180 / 3600
    return (np.interp(mjd_tai, mjd, xp) * as2rad,
            np.interp(mjd_tai, mjd, yp) * as2rad,
            np.interp(mjd_tai, mjd, ut1_tai))


def time_args(gps_s: np.ndarray):
    """Two-part Julian dates (TT, UT1) and polar motion for an array of GPS seconds."""
    gps_s = np.asarray(gps_s, float)
    days = np.floor(gps_s / 86_400.0)
    frac = (gps_s - days * 86_400.0) / 86_400.0
    tt1 = JD_GPS_EPOCH + days
    tt2 = frac + (19.0 + 32.184) / 86_400.0          # TT = TAI + 32.184 s = GPS + 51.184 s
    mjd_tai = (tt1 - 2_400_000.5) + frac + 19.0 / 86_400.0
    xp, yp, ut1_tai = eop(mjd_tai)
    ut2 = frac + (19.0 + ut1_tai) / 86_400.0        # UT1 = TAI + (UT1 - TAI)
    return tt1, tt2, tt1, ut2, xp, yp


def terrestrial_to_celestial(gps_s: np.ndarray) -> np.ndarray:
    """Rotation matrices R (N,3,3) with r_GCRS = R @ r_ITRS (IAU 2006/2000A, IERS EOP)."""
    tt1, tt2, ut1, ut2, xp, yp = time_args(gps_s)
    rc2t = erfa.c2t06a(tt1, tt2, ut1, ut2, xp, yp)
    return np.transpose(rc2t, (0, 2, 1))


def sun_moon_gcrs(gps_s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Geocentric positions (m) of the Sun and the Moon; arcminute accuracy is plenty for tides."""
    tt1, tt2 = time_args(gps_s)[:2]
    pvh, _ = erfa.epv00(tt1, tt2)              # Earth, heliocentric (au), ICRS axes
    sun = -pvh["p"] * C.AU
    moon = erfa.moon98(tt1, tt2)["p"] * C.AU   # GCRS (au)
    return sun, moon


# ---------------------------------------------------------------------------
# Derivatives of uniformly sampled data (Fornberg finite-difference weights)
# ---------------------------------------------------------------------------

def fornberg_weights(x0: float, x: np.ndarray, m: int) -> np.ndarray:
    """Weights for the m-th derivative at x0 from samples at x (Fornberg 1988)."""
    n = len(x)
    c = np.zeros((n, m + 1))
    c1, c4 = 1.0, x[0] - x0
    c[0, 0] = 1.0
    for i in range(1, n):
        mn = min(i, m)
        c2, c5, c4 = 1.0, c4, x[i] - x0
        for j in range(i):
            c3 = x[i] - x[j]
            c2 *= c3
            if j == i - 1:
                for k in range(mn, 0, -1):
                    c[i, k] = c1 * (k * c[i - 1, k - 1] - c5 * c[i - 1, k]) / c2
                c[i, 0] = -c1 * c5 * c[i - 1, 0] / c2
            for k in range(mn, 0, -1):
                c[j, k] = (c4 * c[j, k] - k * c[j, k - 1]) / c3
            c[j, 0] = c4 * c[j, 0] / c3
        c1 = c2
    return c[:, m]


@functools.lru_cache(maxsize=None)
def _stencil(npts: int, offset: int) -> np.ndarray:
    """First-derivative weights on unit-spaced nodes 0..npts-1, evaluated at node ``offset``."""
    return fornberg_weights(float(offset), np.arange(npts, dtype=float), 1)


@functools.lru_cache(maxsize=8)
def _diff_matrix(n: int, npts: int) -> np.ndarray:
    """(n, n) matrix of first-derivative weights (unit spacing): centred stencils
    in the interior, shifted (one-sided) ones near the ends."""
    D = np.zeros((n, n))
    half = npts // 2
    for i in range(n):
        start = min(max(i - half, 0), n - npts)
        D[i, start:start + npts] = _stencil(npts, i - start)
    return D


def derivative(values: np.ndarray, h: float, npts: int = C.LAGRANGE_POINTS) -> np.ndarray:
    """d/dt of uniformly sampled values (axis 0) using npts-point Lagrange stencils.

    Samples containing NaN poison only the stencils that touch them.
    """
    D = _diff_matrix(len(values), npts)
    vals = np.asarray(values, float)
    out = np.tensordot(D, np.nan_to_num(vals), axes=(1, 0)) / h
    if np.isnan(vals).any():
        bad = np.isnan(vals).reshape(len(vals), -1).any(axis=1)
        touched = (np.abs(D[:, bad]) > 0).any(axis=1)
        out[touched] = np.nan
    return out


@functools.lru_cache(maxsize=None)
def _interval_weights(npts: int, offset: int) -> np.ndarray:
    """Weights w with  integral_{offset}^{offset+1} p(x) dx = sum_j w_j f(j), where p
    interpolates f at the unit-spaced nodes 0..npts-1 (exact for degree < npts)."""
    x = np.arange(npts, dtype=float) - (npts - 1) / 2       # centred nodes: better conditioning
    a, b = offset - (npts - 1) / 2, offset + 1 - (npts - 1) / 2
    V = np.vander(x, npts, increasing=True)                  # V[j, k] = x_j**k
    k = np.arange(npts)
    moments = (b ** (k + 1) - a ** (k + 1)) / (k + 1)        # integral of x**k over [a, b]
    return np.linalg.solve(V.T, moments)


@functools.lru_cache(maxsize=8)
def _cumint_matrix(n: int, npts: int) -> np.ndarray:
    """(n, n) matrix giving the cumulative integral (unit spacing, zero at the first sample)."""
    W = np.zeros((n - 1, n))
    half = npts // 2
    for i in range(n - 1):
        start = min(max(i - half + 1, 0), n - npts)
        W[i, start:start + npts] = _interval_weights(npts, i - start)
    return np.vstack([np.zeros(n), np.cumsum(W, axis=0)])


def cumulative_integral(f: np.ndarray, h: float, npts: int = 10) -> np.ndarray:
    """Cumulative integral of uniformly sampled f (axis 0), zero at the first sample.

    Each 5-min interval is integrated exactly over the degree-9 Lagrange polynomial
    through the 10 nearest samples (error far below 0.01 ps for GNSS orbits; cumulative
    Simpson leaves a 0.1 ps odd/even ripple).
    """
    M = _cumint_matrix(len(f), npts)
    return np.tensordot(M, np.asarray(f, float), axes=(1, 0)) * h


# ---------------------------------------------------------------------------
# Relativity along the orbit
# ---------------------------------------------------------------------------

@dataclass
class OrbitDay:
    """One satellite over one day, in the GCRS, with the clock and the GR model terms."""
    prn: str
    gps_s: np.ndarray        # (N,)
    r_itrs: np.ndarray       # (N,3) m
    r: np.ndarray            # (N,3) m, GCRS
    v: np.ndarray            # (N,3) m/s, GCRS
    clk: np.ndarray          # (N,) s, as published (IGS convention)
    terms: dict              # integrands (dimensionless) and their integrals (s)


def potentials(r_itrs: np.ndarray, r_gcrs: np.ndarray, sun: np.ndarray, moon: np.ndarray) -> dict:
    """Positive Newtonian potentials (m^2/s^2) at the satellite."""
    rr = np.linalg.norm(r_itrs, axis=1)
    sinphi = r_itrs[:, 2] / rr
    u_mono = C.GM_EARTH / rr
    u_j2 = -C.GM_EARTH * C.J2_EARTH * C.R_EARTH_EQ**2 / rr**3 * 0.5 * (3 * sinphi**2 - 1)
    u_tide = np.zeros_like(rr)
    for body, gm in ((sun, C.GM_SUN), (moon, C.GM_MOON)):
        dist_b = np.linalg.norm(body, axis=1)
        dist_sb = np.linalg.norm(body - r_gcrs, axis=1)
        # tidal part: remove the constant and the uniform-field (free-fall) terms
        u_tide += gm * (1 / dist_sb - 1 / dist_b - np.einsum("ij,ij->i", r_gcrs, body) / dist_b**3)
    return {"mono": u_mono, "j2": u_j2, "tide": u_tide}


def relativity_terms(h: float, r_itrs: np.ndarray, r: np.ndarray, sun: np.ndarray,
                     moon: np.ndarray) -> tuple[np.ndarray, dict]:
    """GCRS velocity and the proper-time terms for positions sampled every h seconds.

    Returns v (N,3) and a dict with, for each part of d(tau)/dt -- "mono" (GM/r),
    "j2" (flattening), "tide" (Moon + Sun) and "kin" (v^2/2) -- the rate
    ("rate_<part>", dimensionless) and the accumulated proper-time deviation
    ("<part>", seconds, mean rate removed: constant rates are absorbed by each
    day's clock fit), plus the IGS convention term -2 r.v/c^2.
    """
    v = derivative(r, h)
    pot = potentials(r_itrs, r, sun, moon)
    c2 = C.C_LIGHT**2
    integrands = {
        "mono": pot["mono"] / c2,
        "j2": pot["j2"] / c2,
        "tide": pot["tide"] / c2,
        "kin": 0.5 * np.einsum("ij,ij->i", v, v) / c2,
    }
    terms = {}
    for key, f in integrands.items():
        terms[key] = -cumulative_integral(f - f.mean(), h)
        terms["rate_" + key] = f
    terms["kepler_convention"] = -2.0 * np.einsum("ij,ij->i", r, v) / c2
    return v, terms


def orbit_day(day: SP3Day, prn: str, rot: np.ndarray | None = None,
              sun_moon: tuple | None = None) -> OrbitDay | None:
    """Orbit, clock and GR terms of one satellite over one SP3 day (None if the orbit has gaps)."""
    if prn not in day.pos:
        return None
    r_itrs = day.pos[prn]
    if np.isnan(r_itrs).any():
        # orbit gaps are rare; interpolating over them would hide problems, so skip the day
        return None
    gps_s = day.gps_s
    h = float(np.median(np.diff(gps_s)))
    if not np.allclose(np.diff(gps_s), h):
        raise ValueError(f"{day.path}: non-uniform epochs")
    if rot is None:
        rot = terrestrial_to_celestial(gps_s)
    r = np.einsum("nij,nj->ni", rot, r_itrs)
    sun, moon = sun_moon if sun_moon is not None else sun_moon_gcrs(gps_s)
    v, terms = relativity_terms(h, r_itrs, r, sun, moon)
    return OrbitDay(prn=prn, gps_s=gps_s, r_itrs=r_itrs, r=r, v=v,
                    clk=day.clk.get(prn, np.full(len(gps_s), np.nan)), terms=terms)


# ---------------------------------------------------------------------------
# Orbital elements and geometry
# ---------------------------------------------------------------------------

def kepler_elements(r: np.ndarray, v: np.ndarray, gm: float = C.GM_EARTH) -> dict:
    """Osculating elements from GCRS position/velocity (arrays (N,3))."""
    rn = np.linalg.norm(r, axis=1)
    v2 = np.einsum("ij,ij->i", v, v)
    rv = np.einsum("ij,ij->i", r, v)
    a = 1.0 / (2.0 / rn - v2 / gm)
    hvec = np.cross(r, v)
    hn = np.linalg.norm(hvec, axis=1)
    evec = np.cross(v, hvec) / gm - r / rn[:, None]
    e = np.linalg.norm(evec, axis=1)
    inc = np.arccos(hvec[:, 2] / hn)
    node = np.cross(np.array([0.0, 0.0, 1.0]), hvec)
    nn = np.linalg.norm(node, axis=1)
    raan = np.arctan2(node[:, 1], node[:, 0])
    argp = np.arccos(np.clip(np.einsum("ij,ij->i", node, evec) / (nn * e), -1, 1))
    argp = np.where(evec[:, 2] < 0, 2 * np.pi - argp, argp)
    # eccentric anomaly from e cos E = 1 - r/a, e sin E = r.v / sqrt(GM a)
    ecosE = 1.0 - rn / a
    esinE = rv / np.sqrt(gm * a)
    E = np.arctan2(esinE, ecosE)
    M = E - esinE
    # argument of latitude
    u = np.arctan2(r[:, 2] / np.sin(inc),
                   (r[:, 0] * np.cos(raan) + r[:, 1] * np.sin(raan)))
    return {"a": a, "e": e, "i": inc, "raan": raan, "argp": argp, "E": E, "M": M, "u": u,
            "r": rn, "h": hvec / hn[:, None]}


def sun_geometry(r: np.ndarray, v: np.ndarray, sun: np.ndarray) -> dict:
    """Sun elevation above the orbit plane (beta) and whether the satellite is in Earth's shadow."""
    hvec = np.cross(r, v)
    hhat = hvec / np.linalg.norm(hvec, axis=1)[:, None]
    shat = sun / np.linalg.norm(sun, axis=1)[:, None]
    beta = np.arcsin(np.clip(np.einsum("ij,ij->i", hhat, shat), -1, 1))
    along = np.einsum("ij,ij->i", r, shat)
    perp = np.linalg.norm(r - along[:, None] * shat, axis=1)
    shadow = (along < 0) & (perp < C.R_EARTH_EQ)      # cylindrical shadow model
    # angle of the satellite from the Sun direction projected in the orbit plane
    s_in = shat - np.einsum("ij,ij->i", shat, hhat)[:, None] * hhat
    s_in /= np.linalg.norm(s_in, axis=1)[:, None]
    rhat = r / np.linalg.norm(r, axis=1)[:, None]
    cosmu = np.einsum("ij,ij->i", rhat, s_in)
    sinmu = np.einsum("ij,ij->i", np.cross(s_in, rhat), hhat)
    mu = np.arctan2(sinmu, cosmu)       # orbit angle measured from the Sun's projection
    return {"beta": beta, "shadow": shadow, "mu": mu}
