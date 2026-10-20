"""Spin precession in binary pulsars, and the measured-versus-predicted ladder.

A pulsar is a gyroscope: its spin axis is fixed in inertial space unless
something turns it. In a tight binary, general relativity makes the axis precess
about the orbit's angular momentum (geodetic, or "spin-orbit", precession) at
(Barker & O'Connell 1975)

    Omega_p = (2 pi / Pb)^(5/3) T_sun^(2/3) m_c (4 m_p + 3 m_c) / [2 (m_p + m_c)^(4/3)] / (1 - e^2)

with masses in solar masses, Pb the orbital period and T_sun = G M_sun / c^3.
The same physics, at 3,000 times lower speed, turns GP-B's gyroscopes by 6.6"/yr.
As the pulsar's axis precesses, its beam sweeps across our line of sight
differently: pulse profiles change shape, polarisation swings, eclipses change.
"""

from __future__ import annotations

import numpy as np

import config as C
import literature as L

T_SUN = 4.925490947e-6          # s, G M_sun / c^3 (IAU nominal)


def geodetic_rate_deg_yr(pb_d: float, e: float, mp: float, mc: float) -> float:
    """Barker-O'Connell spin precession rate of the star with mass mp (deg/yr)."""
    n = 2 * np.pi / (pb_d * C.SECONDS_PER_DAY)
    om = n ** (5 / 3) * T_SUN ** (2 / 3) * mc * (4 * mp + 3 * mc) / (2 * (mp + mc) ** (4 / 3)) / (1 - e * e)
    return float(np.degrees(om) * C.JULIAN_YEAR_S)


def predictions() -> dict:
    return {name: geodetic_rate_deg_yr(s["pb_d"], s["e"], s["mp"], s["mc"])
            for name, s in L.PULSAR_SYSTEMS.items()}


def ladder(gpb_pred: dict | None = None, slr_result: dict | None = None) -> list[dict]:
    """Measured / predicted for every precession measurement on the ladder.

    rate = GR's predicted precession rate (deg/yr); kind = 'geodetic' or 'frame'.
    """
    pred = predictions()
    rows = []
    mas = 1 / 3.6e6
    # GP-B: geodetic and frame dragging (Everitt et al. 2011)
    ns, ns_e, we, we_e = L.GPB_RESULTS["Joint"]
    rows.append({"label": "GP-B geodetic", "kind": "geodetic", "rate": 6606.1 * mas,
                 "ratio": ns / L.GPB_GR_PREDICTION["ns"], "err_lo": ns_e / 6606.1,
                 "err_hi": ns_e / 6606.1, "year": 2011, "ref": "Everitt et al. 2011"})
    rows.append({"label": "GP-B frame dragging", "kind": "frame", "rate": 39.2 * mas,
                 "ratio": we / L.GPB_GR_PREDICTION["we"], "err_lo": we_e / 39.2,
                 "err_hi": we_e / 39.2, "year": 2011, "ref": "Everitt et al. 2011"})
    # the Moon's orbit round the Sun (lunar laser ranging), best result
    best = L.LLR_DE_SITTER[-1]
    rows.append({"label": "Moon's orbit (LLR)", "kind": "geodetic", "rate": 19.2 * mas,
                 "ratio": 1 + best["dev"], "err_lo": best["err"], "err_hi": best["err"],
                 "year": best["year"], "ref": best["ref"]})
    # pulsars
    for m in L.PULSAR_PRECESSION:
        gr = pred[m["system"]]
        rows.append({"label": f"PSR {m['system']}", "kind": "geodetic", "rate": gr,
                     "ratio": m["value"] / gr, "err_lo": m["minus"] / gr, "err_hi": m["plus"] / gr,
                     "year": int(m["ref"][-4:]), "ref": m["ref"], "gr_published": m["gr"]})
    # laser-ranged satellites: frame dragging of the orbital planes (node rate ~31 mas/yr)
    for m in L.FRAME_DRAGGING_SLR:
        err = m.get("total", float(np.hypot(m["stat"], m["syst"])))
        rows.append({"label": f"{m['sats']} ({m['year']})", "kind": "frame",
                     "rate": 2 * L.LT_NODE_RATE["lageos1"] * mas if "LARES 2" in m["sats"]
                     else L.LT_NODE_RATE["lageos1"] * mas,
                     "ratio": m["mu"], "err_lo": err, "err_hi": err, "year": m["year"],
                     "ref": m["ref"]})
    if slr_result is not None:
        rows.append({"label": "this project: LAGEOS + LARES 2", "kind": "frame",
                     "rate": slr_result["rate_mas_yr"] * mas, "ratio": slr_result["mu"],
                     "err_lo": slr_result["err"], "err_hi": slr_result["err"],
                     "year": 2026, "ref": "61-gyroscope-precession", "ours": True})
    return rows
