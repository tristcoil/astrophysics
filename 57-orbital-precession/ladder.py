"""The precession ladder: one law from an Earth satellite to a black hole.

Every rung is an orbit whose periapsis is seen to turn, placed by how deep in
gravity it lies, x = GM / (c^2 a (1 - e^2)) (M = total mass, a = semi-major
axis of the relative orbit). To first post-Newtonian order, general relativity
predicts a turn of 6 pi x radians per orbit for *any* two bodies, whatever
their mass ratio. Each rung's "measured" value is the relativistic part of its
precession as published (or as derived in this project), with everything
classical (planetary tugs, tides, stellar spins, the Earth's bulge) removed.

Beyond the last rung the formula's weak-field form fails: for a circular orbit
the exact Schwarzschild result is 2 pi [(1 - 6x)^(-1/2) - 1] per orbit, which
diverges at the innermost stable circular orbit, x = 1/6.
"""

from __future__ import annotations

import numpy as np

import binaries
import ephemeris as eph
import pulsars
from config import AU_M, C_SI, GM_EARTH_SI, GM_SUN_SI, J2000_JD, PUBLISHED
from orbit_physics import elements, ppn_precession_factor


def _x_planet(key: str) -> tuple[float, float, float]:
    r, v = eph.heliocentric_state(key, J2000_JD)
    mu = eph.gm("sun") + eph.gm(key)
    el = elements(r, v, mu)
    a, e = float(el["a"][0]), float(el["e"][0])
    return mu / (eph.c_au_per_day() ** 2 * a * (1 - e**2)), a, e


def rungs(binary_res: dict | None = None, pulsar_res: dict | None = None,
          s2_res: dict | None = None) -> list[dict]:
    """All measured rungs with x, the GR prediction per orbit, and measured/predicted."""
    binary_res = binary_res or binaries.analyse()
    pulsar_res = pulsar_res or pulsars.analyse()
    out = []

    def add(name, group, x, ratio, err, source, note="", label_side="right"):
        out.append({"name": name, "group": group, "x": float(x), "gr_per_orbit": 6 * np.pi * x,
                    "ratio": float(ratio), "ratio_err": float(err), "source": source,
                    "note": note, "label_side": label_side})

    # --- the Solar System -------------------------------------------------
    pk, cg = PUBLISHED["park2017"], PUBLISHED["cassini_gamma"]
    ratio = ppn_precession_factor(1 + cg["gamma_minus_1"], 1 + pk["beta_minus_1"])
    err = np.hypot(2 * cg["err"], pk["beta_minus_1_err"]) / 3
    x, _, _ = _x_planet("mercury")
    add("Mercury", "Solar System", x, ratio, err,
        "Park et al. 2017 (beta) with Cassini's gamma: (2 + 2 gamma - beta)/3",
        "MESSENGER ranging")
    ic = PUBLISHED["icarus"]
    x = GM_SUN_SI / (C_SI**2 * ic["a_au"] * AU_M * (1 - ic["e"] ** 2))
    add("Icarus (asteroid)", "Solar System", x, ic["ratio"], ic["err"], ic["source"],
        "radar + optical, 1949-1969; error as the authors' realistic estimate")
    lg = PUBLISHED["lageos2"]
    x = GM_EARTH_SI / (C_SI**2 * lg["a_km"] * 1e3 * (1 - lg["e"] ** 2))
    add("LAGEOS II (satellite)", "Solar System", x, lg["eps"], np.hypot(lg["eps_stat"], lg["eps_sys"]),
        lg["source"], "perigee of a laser-ranged satellite, Earth's field")

    # --- binary stars -----------------------------------------------------
    for name, s in binary_res["systems"].items():
        if not s["used_in_test"]:
            continue
        add(name, "Binary stars", s["compactness"], s["ratio_ours"], s["ratio_err"],
            "Baroch et al. 2021 (Table 9): observed minus classical apsidal motion")
    ew = PUBLISHED["ew_ori"]
    # EW Ori: masses and period from Claret et al. 2021 (paper II), Table 1; e = 0.0758
    # (Clausen et al. 2010), as their Appendix A adopts.
    x_ew, _ = binaries.compactness(1.173 + 1.123, 6.9368442, 0.0758)
    gr_ew = np.degrees(6 * np.pi * x_ew)
    add("EW Ori", "Binary stars", x_ew, ew["wdot_rel_mea"] / gr_ew, ew["err"] / gr_ew, ew["source"])
    dh = binaries.di_her()
    add("DI Her", "Binary stars", dh["compactness"], dh["ratio"], dh["ratio_err"],
        "Claret et al. 2021 observed rate minus their tidal + misaligned-spin terms",
        "spins tipped ~75-80 deg; with aligned spins the ratio would be "
        f"{dh['ratio_if_aligned']:.2f}")

    # --- binary pulsars ---------------------------------------------------
    ht = pulsar_res["b1913"]
    add("Hulse-Taylor pulsar", "Binary pulsars", ht["compactness"], ht["ratio"], ht["ratio_err"],
        "Weisberg & Huang 2016; masses from gamma and the orbital decay (this project)")
    dp = pulsar_res["j0737"]
    add("Double pulsar", "Binary pulsars", dp["compactness"], dp["s_pbdot"]["ratio"],
        dp["s_pbdot"]["ratio_err"],
        "Kramer et al. 2021; masses from the Shapiro shape and orbital decay (this project)")

    # --- the Galactic Centre ---------------------------------------------
    g = PUBLISHED["gravity_s2"][0]
    o = PUBLISHED["s2_orbit_gravity2020"]
    gm = GM_SUN_SI * o["mass_1e6"] * 1e6
    a = o["a_mas"] / 1000 * o["r0_pc"] * AU_M
    x = gm / (C_SI**2 * a * (1 - o["e"] ** 2))
    add("S2 at Sgr A*", "Black hole", x, g["f_sp"], g["err"], g["source"],
        "VLTI/GRAVITY astrometry 2016-2019 + 27 years of adaptive optics")
    return out


def predictions_only() -> list[dict]:
    """Other planets' relativistic precession: predicted, measured only as bounds."""
    out = []
    for key, name in (("venus", "Venus"), ("emb", "Earth"), ("mars", "Mars"), ("jupiter", "Jupiter")):
        x, a, e = _x_planet(key)
        out.append({"name": name, "x": float(x), "gr_per_orbit": float(6 * np.pi * x)})
    return out


def exact_circular(x: np.ndarray) -> np.ndarray:
    """Exact Schwarzschild periapsis advance per orbit for a near-circular orbit, rad."""
    x = np.asarray(x, float)
    out = np.full_like(x, np.nan)
    ok = x < 1 / 6
    out[ok] = 2 * np.pi * (1 / np.sqrt(1 - 6 * x[ok]) - 1)
    return out
