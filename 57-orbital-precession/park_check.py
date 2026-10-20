"""Park et al. (2017)'s own recipe for "the" precession rate, applied to our run.

Park et al. (AJ 153, 121, Sect. 4) extract Mercury's rate from JPL's
MESSENGER-era ephemeris by

1. fixing Mercury's *mean* orbital plane over AD 1000-3000 (pole at ICRF
   RA 280.9876 deg, Dec 61.4481 deg),
2. measuring the perihelion angle in that plane over those 2,000 years,
3. fitting B + rate t + Q t^2 + sum_i (S_i sin nu_i t + C_i cos nu_i t), epoch 2000,
   with 14 frequencies nu_i (combinations of Mercury's, Venus's, the Earth's,
   Jupiter's and Saturn's mean motions; their Table 2) - giving 575.3100"/cy.

Our six-millennium full-physics run covers that window, so the same recipe can
be applied to it. The mean motions are measured from the run itself (a slow
beat like 2 n_Mercury - 5 n_Venus, period 5.7 yr, needs them to ~1e-7).
"""

from __future__ import annotations

import numpy as np

import budget
import ephemeris as eph
import nbody
import orbit_physics as op
from config import ARCSEC_PER_RAD, CENTURY_DAYS, YEAR_DAYS

PARK_POLE_RA_DEG, PARK_POLE_DEC_DEG = 280.9876, 61.4481
PARK_WINDOW_YR = (-1000.0, 1000.0)  # AD 1000 - 3000, relative to 2000
#: Park et al. (2017) Table 2: (k_Mercury, k_Venus, k_Earth, k_Jupiter, k_Saturn)
PARK_FREQUENCIES = [(0, 2, 0, 0, 0), (0, 1, 0, 0, 0), (1, -2, 0, 0, 0), (2, -3, 0, 0, 0),
                    (1, -3, 0, 0, 0), (2, -4, 0, 0, 0), (2, -5, 0, 0, 0), (1, 0, -2, 0, 0),
                    (1, 0, -4, 0, 0), (0, 0, 0, 3, 0), (0, 0, 0, 2, 0), (0, 0, 0, 1, 0),
                    (1, 0, 0, -2, 0), (0, 0, 0, 0, 2)]
PARK_RATE = 575.3100
PARK_Q = -0.04478132  # their t^2 coefficient, "/cy^2


def mean_frame() -> np.ndarray:
    """Rows: x, y, z unit vectors of Park et al.'s mean Mercury orbit frame (ICRF)."""
    z = eph.unit_vector_radec(PARK_POLE_RA_DEG, PARK_POLE_DEC_DEG)
    y = np.cross(z, [1.0, 0.0, 0.0])
    y /= np.linalg.norm(y)
    x = np.cross(y, z)
    return np.vstack([x, y, z])


def mean_motions(rec: dict) -> dict:
    """Mean motions (rad/day) of Mercury, Venus, Earth-Moon, Jupiter, Saturn from a run."""
    out = {}
    keys = rec["meta"]["bodies"]
    for k in ("mercury", "venus", "emb", "jupiter", "saturn"):
        i = keys.index(k)
        s = rec["outer"][:, i, :]
        r, v = s[:, :3].T, s[:, 3:].T
        mu = eph.gm("sun") + eph.gm(k)
        # mean longitude = longitude of the e-vector + mean anomaly, frame-independent enough
        el = op.elements(r, v, mu)
        a, e = el["a"], el["e"]
        rn = np.linalg.norm(r, axis=0)
        cosE = (1 - rn / a) / e
        sinE = (r * v).sum(0) / (e * np.sqrt(mu * a))
        M = np.arctan2(sinE, cosE) - e * sinE
        lam = np.unwrap(el["varpi"] + M)
        t = rec["t_outer"]
        out[k] = float(np.polyfit(t, lam, 1)[0])
    return out


def recipe(run_name: str = "millennia_full", periodic: bool = True) -> dict:
    rec = nbody.load(run_name)
    t = rec["t"]
    m = (t >= PARK_WINDOW_YR[0] * YEAR_DAYS) & (t <= PARK_WINDOW_YR[1] * YEAR_DAYS)
    pp = _perihelia(t[m], rec["mercury"][m])
    F = mean_frame()
    eh = pp["ehat"] @ F.T  # components in the mean frame
    ang = np.unwrap(np.arctan2(eh[:, 1], eh[:, 0]))
    tc = pp["t_days"] / CENTURY_DAYS
    cols = [np.ones_like(tc), tc, tc**2]
    nm = mean_motions(rec)
    order = ("mercury", "venus", "emb", "jupiter", "saturn")
    if periodic:
        for combo in PARK_FREQUENCIES:
            nu = sum(k * nm[b] for k, b in zip(combo, order))  # rad/day
            cols += [np.sin(nu * pp["t_days"]), np.cos(nu * pp["t_days"])]
    A = np.vstack(cols).T
    coef, *_ = np.linalg.lstsq(A, ang, rcond=None)
    resid = ang - A @ coef
    return {"rate": float(coef[1] * ARCSEC_PER_RAD), "Q": float(coef[2] * ARCSEC_PER_RAD),
            "rms_resid_arcsec": float(resid.std() * ARCSEC_PER_RAD), "n_passages": int(len(tc)),
            "periodic_terms": periodic,
            "amplitudes_arcsec": ([float(np.hypot(coef[3 + 2 * i], coef[4 + 2 * i]) * ARCSEC_PER_RAD)
                                   for i in range(len(PARK_FREQUENCIES))] if periodic else None),
            "mean_motions_deg_per_day": {k: float(np.degrees(v)) for k, v in nm.items()}}


def _perihelia(t_days, states):
    """Perihelion passages: time of minimum distance and the osculating e-vector there."""
    from analysis import perihelion_passages
    return perihelion_passages(t_days, states)


def analyse() -> dict:
    with_p = recipe(periodic=True)
    without = recipe(periodic=False)
    ast = budget.asteroid_belt_estimate()
    extra = ast["main_belt_rate"] + ast["kuiper_belt_rate_at_44au"]
    return {"with_periodic_terms": with_p, "quadratic_only": without,
            "plus_asteroid_estimate": with_p["rate"] + extra,
            "park2017": {"rate": PARK_RATE, "Q": PARK_Q},
            "difference": with_p["rate"] + extra - PARK_RATE}
