"""What the Sun's acceleration should be, and what the measured one tells us.

The Sun moves on a nearly circular orbit around the Galactic Centre, so the
Galaxy's gravity pulls it towards the centre with

    a = V_c^2 / R0

where R0 is the distance to the centre and V_c the speed of a *circular* orbit
at the Sun's radius. Careful: V_c is not the Sun's own orbital speed. The Sun
moves ~12 km/s faster than a circular orbit (its peculiar motion); gravity
depends only on where the Sun is, not on how fast it goes, so the force is set
by V_c (Klioner et al. 2021, Sect. 3.1).

Two routes to V_c from sibling projects:

* **H1 + Sgr A*.** Sgr A* sits at rest at the centre, so its apparent proper
  motion, -6.411 mas/yr along the Galactic plane (Reid & Brunthaler 2020), is
  pure reflex of the Sun's orbit: the Sun's angular speed is 4.74 x 6.411 =
  30.39 km/s/kpc. Times H1's R0 = 7.93 kpc that is 241 km/s; minus the Sun's
  12.24 km/s peculiar motion, V_c = 229 km/s.
* **H2's adopted constants** (Reid et al. 2019): R0 = 8.15 kpc, Theta0 = 236 km/s.

Beyond the centripetal pull: the bar and spiral arms add ~0.04-0.1 km/s/Myr
sideways, the Galactic disc pulls the Sun (25 pc above the plane) down by
~0.15 km/s/Myr (McMillan 2017 potential), and the Large Magellanic Cloud pulls
by 0.18-0.45 km/s/Myr (Klioner et al. 2021, Sects. 3.2-3.5).

Run backwards, the measured acceleration is a speedometer and a ruler:
V_c = sqrt(a R0) with H1's R0; R0 itself from a and Sgr A*'s proper motion
alone; the length of the Galactic year 2 pi sqrt(R0 / a); and the mass
inside the Sun's orbit, V_c^2 R0 / G (for a spherical Galaxy).
"""

from __future__ import annotations

import math

import numpy as np

import config as C

#: (km/s)^2 / kpc -> uas/yr of drift.
KMS2_KPC_TO_UAS_YR = C.KMS2_PER_KPC_TO_M_S2 / C.UAS_YR_TO_M_S2


def omega_sun(mu_l_mas_yr: float = C.MU_SGRA_L_MAS_YR) -> float:
    """The Sun's angular speed around the Galactic Centre, km/s/kpc."""
    return -C.K_KMS * mu_l_mas_yr


def accel_uas_yr(vc_kms, r0_kpc):
    return np.asarray(vc_kms) ** 2 / np.asarray(r0_kpc) * KMS2_KPC_TO_UAS_YR


def _summary(samples: np.ndarray) -> dict:
    lo, med, hi = np.percentile(samples, [15.865, 50, 84.135])
    return {"median": float(med), "mean": float(np.mean(samples)), "sd": float(np.std(samples)),
            "lo68": float(lo), "hi68": float(hi)}


def predictions(n: int = 200_000, seed: int = C.SEED) -> dict:
    """Centripetal drift predicted by several sets of Galactic constants
    (Monte Carlo over their Gaussian uncertainties)."""
    rng = np.random.default_rng(seed)
    norm = lambda m, s: rng.normal(m, s, n)  # noqa: E731
    out = {}

    om = -C.K_KMS * norm(C.MU_SGRA_L_MAS_YR, C.MU_SGRA_L_ERR)
    vpec = norm(C.V_SUN_PEC_KMS, C.V_SUN_PEC_ERR_KMS)

    r0 = norm(C.R0_H1_KPC, C.R0_H1_ERR_KPC)
    vc = om * r0 - vpec
    out["h1_sgra"] = {
        "label": "H1's R0 + Sgr A*'s proper motion",
        "inputs": f"R0 = {C.R0_H1_KPC} +- {C.R0_H1_ERR_KPC:.3f} kpc (H1); mu_l(Sgr A*) = "
                  f"{C.MU_SGRA_L_MAS_YR} mas/yr; V_sun = {C.V_SUN_PEC_KMS} km/s",
        "vc_kms": _summary(vc), "g_uas_yr": _summary(accel_uas_yr(vc, r0)),
    }

    r0 = norm(C.R0_H1_KPC, C.R0_H1_ERR_KPC)
    vc = norm(C.OMEGA_SUN_H2, C.OMEGA_SUN_H2_ERR) * r0 - norm(10.7, 1.8)
    out["h1_h2"] = {
        "label": "H1's R0 + H2's solar angular speed",
        "inputs": f"R0 = {C.R0_H1_KPC} kpc (H1); Omega_sun = {C.OMEGA_SUN_H2} +- "
                  f"{C.OMEGA_SUN_H2_ERR} km/s/kpc and V_sun = 10.7 km/s (H2, from Reid et al. 2019)",
        "vc_kms": _summary(vc), "g_uas_yr": _summary(accel_uas_yr(vc, r0)),
    }

    r0 = norm(C.R0_H2_KPC, C.R0_H2_ERR_KPC)
    vc = norm(C.THETA0_H2_KMS, C.THETA0_H2_ERR_KMS)
    out["h2_reid"] = {
        "label": "H2's constants (Reid et al. 2019)",
        "inputs": f"R0 = {C.R0_H2_KPC} +- {C.R0_H2_ERR_KPC} kpc, Theta0 = {C.THETA0_H2_KMS} +- "
                  f"{C.THETA0_H2_ERR_KMS} km/s",
        "vc_kms": _summary(vc), "g_uas_yr": _summary(accel_uas_yr(vc, r0)),
    }

    for key, ref in (("gravity2019", "gravity2019_r0"), ("gravity2022", "gravity2022_r0")):
        m, s = C.REFERENCE[ref]["r0_kpc"]
        r0 = norm(m, s)
        vc = om * r0 - vpec
        out[key] = {
            "label": f"GRAVITY {key[-4:]} R0 + Sgr A*'s proper motion",
            "inputs": f"R0 = {m} kpc ({C.REFERENCE[ref]['ref']})",
            "vc_kms": _summary(vc), "g_uas_yr": _summary(accel_uas_yr(vc, r0)),
        }

    e = C.REFERENCE["eilers2019"]
    vc = norm(e["vc_kms"][0], math.hypot(e["vc_kms"][1], e["vc_sys_frac"] * e["vc_kms"][0]))
    out["eilers2019"] = {
        "label": "Eilers et al. 2019 (red giants)",
        "inputs": f"V_c = {e['vc_kms'][0]} km/s at R0 = {e['r0_kpc']} kpc (+-{e['vc_sys_frac']:.0%} sys.)",
        "vc_kms": _summary(vc), "g_uas_yr": _summary(accel_uas_yr(vc, e["r0_kpc"])),
    }
    return out


def mcmillan17_acceleration() -> dict:
    """Gravitational acceleration at the Sun in McMillan's (2017) Milky Way
    model (its own R0 = 8.21 kpc), radial and vertical, via galpy."""
    from galpy.potential import evaluateRforces, evaluatezforces
    from galpy.potential.mwpotentials import McMillan17
    ro, vo = 8.21, 233.1
    scale = vo**2 / ro  # natural force units -> (km/s)^2/kpc
    out = {}
    for z in (0.0, C.Z_SUN_KPC):
        fr = float(evaluateRforces(McMillan17, 1.0, z / ro, use_physical=False)) * scale
        fz = float(evaluatezforces(McMillan17, 1.0, z / ro, use_physical=False)) * scale
        out[f"z_{z * 1000:.0f}pc"] = {
            "a_R_uas_yr": -fr * KMS2_KPC_TO_UAS_YR,  # towards the centre (positive)
            "a_z_uas_yr": fz * KMS2_KPC_TO_UAS_YR,   # negative = towards the plane
            "a_R_kms_myr": -fr * KMS2_KPC_TO_UAS_YR * C.UAS_YR_TO_KMS_MYR,
            "a_z_kms_myr": fz * KMS2_KPC_TO_UAS_YR * C.UAS_YR_TO_KMS_MYR,
        }
    out["vc_kms"] = vo
    out["r0_kpc"] = ro
    return out


def lmc_pull() -> dict:
    """Direct pull of the LMC on the Solar System, in galactic (X, Y, Z)."""
    l, b = math.radians(C.LMC_L_DEG), math.radians(C.LMC_B_DEG)
    unit = np.array([math.cos(b) * math.cos(l), math.cos(b) * math.sin(l), math.sin(b)])
    out = {}
    for m in C.LMC_MASS_RANGE_MSUN:
        a = C.G_KPC_KMS2_MSUN * m / C.LMC_D_KPC**2 * KMS2_KPC_TO_UAS_YR
        out[f"{m:.1e}"] = {"uas_yr": (a * unit).tolist(),
                           "kms_myr": (a * C.UAS_YR_TO_KMS_MYR * unit).tolist()}
    return out


def r0_from_acceleration(a_uas_yr, om_kms_kpc, vpec_kms):
    """Solve (Omega R0 - V_pec)^2 / R0 = a for R0 (the physical, larger root)."""
    a = np.asarray(a_uas_yr) / KMS2_KPC_TO_UAS_YR  # (km/s)^2/kpc
    om, v = np.asarray(om_kms_kpc), np.asarray(vpec_kms)
    bq = 2 * om * v + a
    return (bq + np.sqrt(np.maximum(bq**2 - 4 * om**2 * v**2, 0))) / (2 * om**2)


def implications(g_samples: np.ndarray, seed: int = C.SEED) -> dict:
    """From samples of the acceleration towards the Galactic Centre (uas/yr):
    the circular speed, R0, the Galactic year and the enclosed mass."""
    rng = np.random.default_rng(seed + 1)
    n = g_samples.size
    a_kpc = g_samples / KMS2_KPC_TO_UAS_YR  # (km/s)^2/kpc
    r0 = rng.normal(C.R0_H1_KPC, C.R0_H1_ERR_KPC, n)
    vc = np.sqrt(np.maximum(a_kpc, 0) * r0)
    period_myr = 2 * np.pi * r0 * C.KPC_M / (vc * 1e3) / (C.JULIAN_YEAR_S * 1e6)
    mass = vc**2 * r0 / C.G_KPC_KMS2_MSUN
    om = -C.K_KMS * rng.normal(C.MU_SGRA_L_MAS_YR, C.MU_SGRA_L_ERR, n)
    vpec = rng.normal(C.V_SUN_PEC_KMS, C.V_SUN_PEC_ERR_KMS, n)
    r0_from_a = r0_from_acceleration(g_samples, om, vpec)
    a_si = g_samples * C.UAS_YR_TO_M_S2
    return {
        "vc_kms_with_h1_r0": _summary(vc),
        "galactic_year_myr": _summary(period_myr),
        "mass_inside_sun_msun": _summary(mass),
        "r0_kpc_from_sgra_motion": _summary(r0_from_a),
        "velocity_change_mm_s_per_year": _summary(a_si * C.JULIAN_YEAR_S * 1e3),
    }
