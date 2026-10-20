"""How good are two-line elements? CelesTrak's current TLEs against the precise orbits.

SGP4 propagates each satellite's latest element set back to its most recent ILRS
weekly orbit; the position differences (radial, along-track, cross-track) show
why TLEs, fine for watching the 125 deg/yr bulge-driven drift, cannot see frame
dragging, which moves the orbital plane by about 1.8 m a year.
"""

from __future__ import annotations

import numpy as np

import config as C
import frames as F
import nodes


def compare(sat: str, ac: str = C.PRIMARY_AC) -> dict:
    from skyfield.api import EarthSatellite, load
    tle = sorted(C.TLE_DIR.glob(f"{sat}_*.tle"))[-1].read_text().splitlines()
    ts = load.timescale(builtin=True)
    es = EarthSatellite(tle[1], tle[2], tle[0].strip(), ts)
    path = nodes.arc_files(sat, ac)[-1]
    o = F.read_sp3(path)
    sel = slice(None, None, 15)                       # every 30 min
    r, v = F.itrs_to_gcrs(o.mjd_utc[sel], o.r_itrs[sel], o.v_itrs[sel])
    t = ts.tt_jd(o.mjd_tt[sel] + 2_400_000.5)
    p = es.at(t).position.m.T                         # GCRS (skyfield applies TEME -> GCRS)
    d = p - r
    R, T, N = (np.einsum("ij,ij->i", d, B) for B in _rtn(r, v))
    epoch_days = float(o.mjd_utc[-1] - (es.epoch.tt - 2_400_000.5))
    return {"sat": sat, "sp3": path.name, "days_from_tle_epoch": epoch_days,
            "rms_m": {"radial": float(np.sqrt(np.mean(R ** 2))), "along": float(np.sqrt(np.mean(T ** 2))),
                      "cross": float(np.sqrt(np.mean(N ** 2)))},
            "max_m": float(np.max(np.linalg.norm(d, axis=1)))}


def _rtn(r, v):
    Rh = r / np.linalg.norm(r, axis=1)[:, None]
    Nh = np.cross(r, v)
    Nh /= np.linalg.norm(Nh, axis=1)[:, None]
    return Rh, np.cross(Nh, Rh), Nh


if __name__ == "__main__":
    for s in ("lageos1", "lares2"):
        print(compare(s))
