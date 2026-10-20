"""Time scales, Earth orientation, and the ILRS SP3 orbit files.

The weekly orbits come in the Earth-fixed frame (SLRF2020, a realisation of the
ITRF). Frame dragging and the bulge-driven drift are motions of the orbital
plane *in inertial space*, so every orbit is rotated into the Geocentric
Celestial Reference System (GCRS) with the IAU 2006/2000A precession-nutation
model (via ERFA) and the IERS Earth orientation parameters:

    r_GCRS = Q(t) R3(-ERA) W(t) r_ITRS

Q: precession-nutation (CIO-based, X and Y corrected by the C04 offsets dX, dY);
ERA: Earth rotation angle from UT1; W: polar motion. UT1 matters most here: the
node of an orbit is an angle around the Earth's axis, and an error of 1 ms in
UT1 would rotate it by 15 milliarcseconds. The C04 series is built from VLBI of
quasars, independent of satellite laser ranging.

Times: SP3 epochs are UTC. We work in TT (= UTC + 37 s + 32.184 s over the whole
2022-2026 span), as seconds since a reference epoch.
"""

from __future__ import annotations

import datetime as dt
import functools
import gzip
from dataclasses import dataclass
from pathlib import Path

import erfa
import numpy as np

import config as C

AS2RAD = np.pi / 180.0 / 3600.0
MJD0 = 2_400_000.5


def tai_minus_utc(mjd_utc) -> np.ndarray:
    mjd_utc = np.asarray(mjd_utc, float)
    out = np.full(mjd_utc.shape, np.nan)
    for start, value in C.TAI_MINUS_UTC:
        out = np.where(mjd_utc >= (start - dt.date(1858, 11, 17)).days, value, out)
    if np.any(np.isnan(out)):
        raise ValueError("date before the leap-second table in config")
    return out


def utc_to_tt_mjd(mjd_utc) -> np.ndarray:
    return np.asarray(mjd_utc, float) + (tai_minus_utc(mjd_utc) + C.TT_MINUS_TAI) / C.SECONDS_PER_DAY


def tt_mjd_to_utc_mjd(mjd_tt) -> np.ndarray:
    first = np.asarray(mjd_tt, float) - (37.0 + C.TT_MINUS_TAI) / C.SECONDS_PER_DAY
    return np.asarray(mjd_tt, float) - (tai_minus_utc(first) + C.TT_MINUS_TAI) / C.SECONDS_PER_DAY


def date_to_mjd(d: dt.date | dt.datetime) -> float:
    if isinstance(d, dt.datetime):
        return (d - dt.datetime(1858, 11, 17)).total_seconds() / C.SECONDS_PER_DAY
    return float((d - dt.date(1858, 11, 17)).days)


def mjd_to_year(mjd) -> np.ndarray:
    """Decimal (Julian) year, for plotting."""
    return 2000.0 + (np.asarray(mjd, float) - 51_544.5) / 365.25


# ---------------------------------------------------------------------------
# Earth orientation parameters
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=1)
def eop_table() -> dict:
    """Daily EOP: C04 (to its last day), then the rapid series (finals2000A).

    Returns MJD (UTC), xp, yp (rad), UT1-TAI (s, smooth across leap seconds),
    dX, dY (rad) and a flag saying which series each day comes from.
    """
    c04 = np.loadtxt(C.EOP_DIR / "eopc04.1962-now", comments="#", usecols=(4, 5, 6, 7, 8, 9))
    mjd, xp, yp, dut1, dx, dy = c04.T
    keep = mjd >= 59_000          # 2020-05-31: plenty of margin before our span
    mjd, xp, yp, dut1, dx, dy = (a[keep] for a in (mjd, xp, yp, dut1, dx, dy))
    source = np.zeros(mjd.size, int)

    # finals2000A.all (fixed columns, IERS Bulletin A/B): MJD 8-15, x 19-27, y 38-46,
    # UT1-UTC 59-68, dX 98-106, dY 117-125 (milliarcsec for dX, dY).
    rows = []
    for line in (C.EOP_DIR / "finals2000A.all").read_text().splitlines():
        if len(line) < 125 or not line[58:68].strip() or not line[97:106].strip():
            continue
        m = float(line[7:15])
        if m <= mjd[-1]:
            continue
        rows.append((m, float(line[18:27]), float(line[37:46]), float(line[58:68]),
                     float(line[97:106]) / 1e3, float(line[116:125]) / 1e3))
    if rows:
        f = np.array(rows)
        mjd = np.concatenate([mjd, f[:, 0]])
        xp, yp, dut1 = (np.concatenate([a, f[:, i]]) for a, i in ((xp, 1), (yp, 2), (dut1, 3)))
        dx, dy = (np.concatenate([a, f[:, i]]) for a, i in ((dx, 4), (dy, 5)))
        source = np.concatenate([source, np.ones(len(rows), int)])
    return {"mjd": mjd, "xp": xp * AS2RAD, "yp": yp * AS2RAD,
            "ut1_tai": dut1 - tai_minus_utc(mjd), "dx": dx * AS2RAD, "dy": dy * AS2RAD,
            "source": source, "c04_end": float(mjd[source == 0][-1])}


def eop_at(mjd_utc) -> dict:
    """EOP interpolated to the given UTC epochs.

    Polar motion and UT1 vary smoothly between the daily samples; cubic (Lagrange,
    4-point) interpolation is used as the IERS recommends for UT1 and the pole.
    """
    t = eop_table()
    x = np.asarray(mjd_utc, float)
    i = np.clip(np.searchsorted(t["mjd"], x) - 2, 0, t["mjd"].size - 4)
    out = {}
    for key in ("xp", "yp", "ut1_tai", "dx", "dy"):
        y = t[key]
        acc = np.zeros_like(x)
        for j in range(4):
            w = np.ones_like(x)
            for k in range(4):
                if k != j:
                    w *= (x - t["mjd"][i + k]) / (t["mjd"][i + j] - t["mjd"][i + k])
            acc += w * y[i + j]
        out[key] = acc
    return out


def itrs_to_gcrs_matrices(mjd_utc) -> np.ndarray:
    """Rotation matrices R (N,3,3) with r_GCRS = R @ r_ITRS (IAU 2006/2000A + C04)."""
    mjd_utc = np.atleast_1d(np.asarray(mjd_utc, float))
    e = eop_at(mjd_utc)
    tai = mjd_utc + tai_minus_utc(mjd_utc) / C.SECONDS_PER_DAY
    tt2 = tai + C.TT_MINUS_TAI / C.SECONDS_PER_DAY
    ut12 = tai + e["ut1_tai"] / C.SECONDS_PER_DAY
    x, y = erfa.xy06(MJD0, tt2)
    x = x + e["dx"]
    y = y + e["dy"]
    s = erfa.s06(MJD0, tt2, x, y)
    rc2i = erfa.c2ixys(x, y, s)
    era = erfa.era00(MJD0, ut12)
    rpom = erfa.pom00(e["xp"], e["yp"], erfa.sp00(MJD0, tt2))
    rc2t = erfa.c2tcio(rc2i, era, rpom)
    return np.transpose(rc2t, (0, 2, 1))


def itrs_to_gcrs(mjd_utc, r_itrs, v_itrs=None):
    """Positions (and velocities) from the Earth-fixed frame to the GCRS.

    v_GCRS = Q R3(-ERA) [W v + omega z x (W r)]: the Earth's rotation adds
    omega x r. The rotation rate is the IERS nominal one corrected by the
    length of day implied by UT1 (relative change ~1e-8, irrelevant here).
    """
    R = itrs_to_gcrs_matrices(mjd_utc)
    r = np.einsum("nij,nj->ni", R, r_itrs)
    if v_itrs is None:
        return r
    e = eop_at(mjd_utc)
    tt2 = (np.asarray(mjd_utc, float) + (tai_minus_utc(mjd_utc) + C.TT_MINUS_TAI) / C.SECONDS_PER_DAY)
    W = np.transpose(erfa.pom00(e["xp"], e["yp"], erfa.sp00(MJD0, tt2)), (0, 2, 1))  # ITRS->TIRS
    r_tirs = np.einsum("nij,nj->ni", W, r_itrs)
    v_tirs = np.einsum("nij,nj->ni", W, v_itrs)
    v_tirs = v_tirs + C.OMEGA_EARTH * np.stack([-r_tirs[:, 1], r_tirs[:, 0], np.zeros(len(r_tirs))], axis=1)
    # R = Q R3(-ERA) W^T... compose: GCRS<-TIRS = R @ W^T (W is orthogonal)
    R_tirs = np.einsum("nij,nkj->nik", R, W)
    v = np.einsum("nij,nj->ni", R_tirs, v_tirs)
    return r, v


def gcrs_to_itrs(mjd_utc, r_gcrs, v_gcrs):
    """Inverse of itrs_to_gcrs: GCRS states to Earth-fixed positions and velocities."""
    R = itrs_to_gcrs_matrices(mjd_utc)
    r_itrs = np.einsum("nji,nj->ni", R, r_gcrs)
    e = eop_at(mjd_utc)
    tt2 = (np.asarray(mjd_utc, float) + (tai_minus_utc(mjd_utc) + C.TT_MINUS_TAI) / C.SECONDS_PER_DAY)
    W = np.transpose(erfa.pom00(e["xp"], e["yp"], erfa.sp00(MJD0, tt2)), (0, 2, 1))  # ITRS->TIRS
    R_tirs = np.einsum("nij,nkj->nik", R, W)                                           # GCRS<-TIRS
    r_tirs = np.einsum("nij,nj->ni", W, r_itrs)
    v_tirs = np.einsum("nji,nj->ni", R_tirs, v_gcrs) - C.OMEGA_EARTH * np.stack(
        [-r_tirs[:, 1], r_tirs[:, 0], np.zeros(len(r_tirs))], axis=1)
    return r_itrs, np.einsum("nji,nj->ni", W, v_tirs)


# ---------------------------------------------------------------------------
# SP3 files
# ---------------------------------------------------------------------------

@dataclass
class Orbit:
    """One weekly orbit: UTC epochs (MJD) and Earth-fixed states (m, m/s)."""

    sat: str
    ac: str
    week: str
    path: Path
    mjd_utc: np.ndarray
    r_itrs: np.ndarray
    v_itrs: np.ndarray
    header: str

    @property
    def mjd_tt(self) -> np.ndarray:
        return utc_to_tt_mjd(self.mjd_utc)


def read_sp3(path: Path) -> Orbit:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", errors="replace") as fh:
        lines = fh.read().splitlines()
    if not lines or not lines[0].startswith("#c"):
        raise ValueError(f"{path.name}: not an SP3-c file")
    time_system = lines[12][9:12] if len(lines) > 12 else "UTC"
    if time_system.strip() not in ("UTC", "ccc", ""):
        raise ValueError(f"{path.name}: time system {time_system!r} not UTC")
    epochs, pos, vel = [], [], []
    header = []
    for line in lines:
        if line.startswith("*  "):
            y, mo, d, h, mi = (int(v) for v in line[2:31].split()[:5])
            sec = float(line[2:31].split()[5])
            epochs.append(date_to_mjd(dt.datetime(y, mo, d, h, mi)) + sec / C.SECONDS_PER_DAY)
        elif line.startswith("P") and epochs:
            pos.append([float(line[4 + 14 * k:18 + 14 * k]) for k in range(3)])
        elif line.startswith("V") and epochs:
            vel.append([float(line[4 + 14 * k:18 + 14 * k]) for k in range(3)])
        elif not epochs:
            header.append(line)
    pos = np.array(pos) * 1e3            # km -> m
    vel = np.array(vel) * 1e-1           # dm/s -> m/s
    n = min(len(epochs), len(pos), len(vel))
    parts = path.name.split(".")
    return Orbit(sat=parts[2], ac=parts[0], week=parts[3], path=path,
                 mjd_utc=np.array(epochs[:n]), r_itrs=pos[:n], v_itrs=vel[:n],
                 header="\n".join(header))


def orbit_files(sat: str, ac: str) -> list[Path]:
    """The cached weekly SP3 files of one satellite and AC, oldest first (latest version each)."""
    best = {}
    for p in sorted((C.SP3_DIR / sat / ac).glob("*.sp3.gz")):
        week = p.name.split(".")[3]
        ver = int(p.name.split(".")[4][1:])
        if week not in best or ver > best[week][0]:
            best[week] = (ver, p)
    return [best[w][1] for w in sorted(best)]


# ---------------------------------------------------------------------------
# Orbital elements
# ---------------------------------------------------------------------------

def kepler_elements(r: np.ndarray, v: np.ndarray, gm: float = C.GM_EARTH) -> dict:
    """Osculating Keplerian elements (angles in rad) of GCRS states (N,3)."""
    r = np.atleast_2d(r)
    v = np.atleast_2d(v)
    h = np.cross(r, v)
    hn = np.linalg.norm(h, axis=1)
    rn = np.linalg.norm(r, axis=1)
    inc = np.arccos(np.clip(h[:, 2] / hn, -1, 1))
    node = np.arctan2(h[:, 0], -h[:, 1])
    energy = 0.5 * np.sum(v * v, axis=1) - gm / rn
    a = -gm / (2 * energy)
    evec = np.cross(v, h) / gm - r / rn[:, None]
    e = np.linalg.norm(evec, axis=1)
    nvec = np.stack([np.cos(node), np.sin(node), np.zeros_like(node)], axis=1)
    w = np.arctan2(np.einsum("ij,ij->i", np.cross(nvec, evec), h / hn[:, None]),
                   np.einsum("ij,ij->i", nvec, evec))
    u = np.arctan2(np.einsum("ij,ij->i", np.cross(nvec, r), h / hn[:, None]),
                   np.einsum("ij,ij->i", nvec, r))
    return {"a": a, "e": e, "i": inc, "node": node, "argp": w, "arglat": u,
            "h": h / hn[:, None]}
