"""Cross-checks against JPL Horizons, which serves DE441 (DE440's long-span twin).

1. **Our elements = Horizons' elements.** Horizons' own osculating elements for
   Mercury at 109 epochs (1560-2640) against ours computed from DE440 state
   vectors with GM = GM_sun + GM_mercury and the J2000 ecliptic -- the two
   conventions Horizons states in its headers.
2. **DE440 = DE441 here.** Mercury's position from both ephemerides at the same
   epochs: JPL fitted both to the same observations; DE441 drops a lunar
   core-mantle damping term so it can be integrated over 30,000 years.
3. **Our integration vs DE441 over six millennia** (1000 BC - AD 5000): the
   direction of Mercury's perihelion and its position, every 25 years.
"""

from __future__ import annotations

import numpy as np

import ephemeris as eph
import orbit_physics as op
from config import ARCSEC_PER_RAD, CENTURY_DAYS, J2000_JD
from fetch_data import horizons_files


def _rows(path) -> list[list[str]]:
    text = path.read_text()
    if "$$SOE" not in text:
        raise ValueError(f"{path.name}: no data block")
    body = text.split("$$SOE")[1].split("$$EOE")[0]
    return [[c.strip() for c in ln.split(",")] for ln in body.strip().splitlines()]


def load_vectors(label: str) -> dict:
    """JD (TDB), heliocentric ecliptic-J2000 position (au) and velocity (au/day)."""
    rows = [r for f in horizons_files(label) for r in _rows(f)]
    jd = np.array([float(r[0]) for r in rows])
    xyz = np.array([[float(x) for x in r[2:8]] for r in rows])
    order = np.argsort(jd)
    return {"jd": jd[order], "r": xyz[order, :3].T, "v": xyz[order, 3:].T}


def load_elements(label: str) -> dict:
    rows = [r for f in horizons_files(label) for r in _rows(f)]
    cols = ["EC", "QR", "IN", "OM", "W", "Tp", "N", "MA", "TA", "A", "AD", "PR"]
    data = np.array([[float(x) for x in r[2:2 + len(cols)]] for r in rows])
    jd = np.array([float(r[0]) for r in rows])
    order = np.argsort(jd)
    out = {c: data[order, i] for i, c in enumerate(cols)}
    out["jd"] = jd[order]
    return out


def _wrap_deg(x):
    return (x + 180.0) % 360.0 - 180.0


def elements_check() -> dict:
    """Horizons' elements vs ours from DE440 at the same epochs."""
    hz = load_elements("mercury_elements_span")
    r, v = eph.heliocentric_state("mercury", hz["jd"])
    rot = eph.rot_icrf_to_ecliptic()
    el = op.elements(rot @ r, rot @ v, eph.gm("sun") + eph.gm("mercury"))
    d = {"a_au": el["a"] - hz["A"], "e": el["e"] - hz["EC"],
         "inc_arcsec": (np.degrees(el["inc"]) - hz["IN"]) * 3600,
         "Omega_arcsec": _wrap_deg(np.degrees(el["Omega"]) - hz["OM"]) * 3600,
         "omega_arcsec": _wrap_deg(np.degrees(el["omega"]) - hz["W"]) * 3600}
    return {"n_epochs": int(len(hz["jd"])),
            "max_abs_difference": {k: float(np.abs(x).max()) for k, x in d.items()},
            "years": [float(2000 + (hz["jd"][0] - J2000_JD) / 365.25),
                      float(2000 + (hz["jd"][-1] - J2000_JD) / 365.25)]}


def de440_vs_de441() -> dict:
    hz = load_vectors("mercury_vectors_span")
    r, _ = eph.heliocentric_state("mercury", hz["jd"])
    r = eph.rot_icrf_to_ecliptic() @ r
    dr = np.linalg.norm(r - hz["r"], axis=0) * eph.au_km()
    return {"n_epochs": int(len(hz["jd"])), "max_position_difference_km": float(dr.max()),
            "rms_position_difference_km": float(np.sqrt((dr**2).mean())),
            "difference_km_by_year": {f"{2000 + (j - J2000_JD) / 365.25:.0f}": float(x)
                                      for j, x in zip(hz["jd"][::10], dr[::10])}}


def millennia_vs_de441(run_name: str) -> dict:
    """Epoch-by-epoch comparison of a six-millennium run with DE441."""
    import nbody

    hz = load_vectors("mercury_vectors_long")
    rec = nbody.load(run_name)
    t_hz = hz["jd"] - J2000_JD
    idx = np.searchsorted(rec["t"], t_hz)
    idx = np.clip(idx, 0, len(rec["t"]) - 1)
    if np.abs(rec["t"][idx] - t_hz).max() > 1e-6:
        raise ValueError("run is not sampled at the Horizons epochs")
    rot = eph.rot_icrf_to_ecliptic()
    s = rec["mercury"][idx]
    mu = eph.gm("sun") + eph.gm("mercury")
    ours = op.elements(rot @ s[:, :3].T, rot @ s[:, 3:].T, mu)
    theirs = op.elements(hz["r"], hz["v"], mu)
    cross = np.cross(theirs["ehat"].T, ours["ehat"].T).T
    dpsi = np.arctan2((cross * theirs["hhat"]).sum(0), (ours["ehat"] * theirs["ehat"]).sum(0))
    dpos = np.linalg.norm(rot @ s[:, :3].T - hz["r"], axis=0) * eph.au_km()
    years = 2000 + t_hz / 365.25
    return {"years": years.tolist(), "perihelion_offset_arcsec": (dpsi * ARCSEC_PER_RAD).tolist(),
            "position_offset_km": dpos.tolist(),
            "offset_rate_arcsec_per_cy": float(np.polyfit(t_hz / CENTURY_DAYS, dpsi * ARCSEC_PER_RAD, 1)[0])}
