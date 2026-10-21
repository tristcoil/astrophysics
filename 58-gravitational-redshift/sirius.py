"""Sirius B: a white dwarf's gravitational redshift, predicted from our own mass and radius.

Light leaving a white dwarf climbs out of a potential ~400 times deeper than the
Sun's surface: GR predicts a redshift v = GM/(R c), tens of km/s. M1
(50-binary-orbits) weighed Sirius B from its 50-year orbit (Bond et al.'s 2,367
measures) and sized it from its brightness and temperature; its prediction is
compared here with the redshift Joyce et al. (2018) measured in Sirius B's
H-alpha line with HST/STIS, after removing its orbital velocity.

The exact Schwarzschild redshift, 1/sqrt(1 - 2GM/(Rc^2)) - 1, differs from
GM/(Rc^2) by only 0.03 km/s here.
"""

from __future__ import annotations

import json

import numpy as np

import config as C

KMS = C.C_LIGHT / 1e3
R_SUN_M = C.R_SUN

#: Joyce et al. 2018 (MNRAS 481, 2361), HST/STIS H-alpha: Sirius B's gravitational redshift.
JOYCE2018 = (80.65, 0.77)

#: Bond et al. 2017 (ApJ 840, 70) adopted mass and radius, for comparison.
BOND2017 = {"M": (1.018, 0.011), "R": (0.008098, 0.000046)}

#: Used only if M1's results file is missing (its published values).
M1_FALLBACK = {"MB": (1.0150, 0.0114), "R": (0.007895, 0.000184), "vgr": (81.81, 1.04)}


def v_gr(M_sun: float, R_sun: float) -> float:
    """GM/(Rc) in km/s."""
    return C.GM_SUN * M_sun / (R_sun * R_SUN_M * C.C_LIGHT) / 1e3


def v_exact(M_sun: float, R_sun: float) -> float:
    x = C.GM_SUN * M_sun / (R_sun * R_SUN_M * C.C_LIGHT**2)
    return (1 / np.sqrt(1 - 2 * x) - 1) * KMS


def analyse() -> dict:
    if C.M1_RESULTS.exists():
        m1 = json.loads(C.M1_RESULTS.read_text())["sirius_b"]
        MB, vgr = tuple(m1["MB"]), tuple(m1["vgr"])
        R = tuple(m1["R"]["V (HST/STIS, Barstow et al. 2005)"])
        R_gaia = tuple(m1["R"]["Gaia G (DR3)"])
        source = "50-binary-orbits/data/results.json"
    else:
        MB, R, vgr = M1_FALLBACK["MB"], M1_FALLBACK["R"], M1_FALLBACK["vgr"]
        R_gaia = (np.nan, np.nan)
        source = "M1 fallback values (results file missing)"
    meas, emeas = JOYCE2018
    ratio = meas / vgr[0]
    eratio = ratio * np.hypot(emeas / meas, vgr[1] / vgr[0])
    return {
        "source": source, "M_B": MB, "R_B_hst": R, "R_B_gaia": R_gaia,
        "predicted_kms": vgr, "predicted_gaia_radius_kms": v_gr(MB[0], R_gaia[0]) if np.isfinite(R_gaia[0]) else None,
        "predicted_bond2017_kms": v_gr(BOND2017["M"][0], BOND2017["R"][0]),
        "exact_minus_weak_field_kms": v_exact(MB[0], R[0]) - v_gr(MB[0], R[0]),
        "measured_kms": JOYCE2018, "ratio": ratio, "ratio_err": eratio,
        "z_measured": meas / KMS, "potential_over_c2": vgr[0] / KMS,
    }
