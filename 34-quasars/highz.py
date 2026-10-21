"""The growth problem: how did black holes reach a billion suns within a billion years?

For every quasar beyond z = 5.7 with a measured Mg II line we recompute its mass
with one recipe (Shen et al. 2011, the one Wu & Shen use for DR16Q), then ask:
starting at redshift z_start, growing without pause at the Eddington limit with
10% radiative efficiency, how big must the seed have been?

    M(t) = M_seed exp(t / t_e),  t_e = eps / (1 - eps) * 450 Myr = 50 Myr.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import fetch_data
import literature as lit
import masses
import physics
from config import HIGHZ_MIN_Z, HIGHZ_TABLES, RADIATIVE_EFFICIENCY, SEED_Z_START


def load() -> pd.DataFrame:
    """Farina et al. 2022 Table 5: 114 quasars at z > 5.7 with Mg II widths and L3000."""
    t = fetch_data.load_vizier_table(HIGHZ_TABLES["farina22"])
    t = t[t["zSys"] >= HIGHZ_MIN_Z].copy()
    out = pd.DataFrame({
        "name": t["ID"].str.strip(),
        "z": t["zSys"].astype(float),
        "M1450": t["M1450"].astype(float),
        "fwhm_mgii": t["FWHMMgII"].astype(float),
        "log_l3000": np.log10(t["L3000"].astype(float)) + 46.0,   # table unit: 1e39 W
        "log_lbol_pub": np.log10(t["Lbol"].astype(float)) + 46.0,
        "log_mbh_pub": np.log10(t["MBH-S11"].astype(float)) + 9.0,
        "edd_pub": t["LEdd-S11"].astype(float),
        "ref": t["Ref"].str.strip(),
    }).reset_index(drop=True)
    out["log_mbh"] = masses.virial_mass(out["fwhm_mgii"], out["log_l3000"], "MGII_S11")
    out["log_mbh_vo09"] = masses.virial_mass(out["fwhm_mgii"], out["log_l3000"], "MGII_VO09")
    out["log_lbol"] = masses.bolometric(out["log_l3000"], 3000)
    out["log_edd"] = masses.eddington_ratio(out["log_lbol"], out["log_mbh"])
    out["age_myr"] = physics.age_myr(out["z"])
    return out


def seed_table(hz: pd.DataFrame, eps: float = RADIATIVE_EFFICIENCY,
               z_starts=SEED_Z_START) -> pd.DataFrame:
    """Minimum seed for each quasar: Eddington-limited, uninterrupted growth from z_start."""
    out = hz[["name", "z", "log_mbh", "log_edd", "age_myr"]].copy()
    m = 10.0 ** out["log_mbh"]
    for zs in z_starts:
        dt = out["age_myr"] - physics.age_myr(zs)
        out[f"dt_{int(zs)}"] = dt
        out[f"log_seed_{int(zs)}"] = np.log10(physics.seed_mass(m, dt, eps, 1.0))
        # growing at the Eddington ratio we observe now, all the time
        lam = 10.0 ** out["log_edd"]
        out[f"log_seed_obs_{int(zs)}"] = np.log10(physics.seed_mass(m, dt, eps, lam))
    out["efolds_from_100"] = np.log(m / 100.0)
    out["efolds_available_30"] = out["dt_30"] / physics.efold_time_myr(eps)
    return out


def growth_track(log_m_obs: float, z_obs: float, z_start: float = 30.0, eps=RADIATIVE_EFFICIENCY,
                 n: int = 300):
    """Mass vs cosmic time going back from the observed quasar to z_start (Eddington-limited)."""
    t_obs = float(physics.age_myr(z_obs))
    t = np.linspace(float(physics.age_myr(z_start)), t_obs, n)
    logm = log_m_obs - (t_obs - t) / physics.efold_time_myr(eps) / np.log(10.0)
    return t, logm


def j0313_checks() -> dict:
    """J0313-1806 three ways: Wang et al.'s numbers, and their width with our recipes."""
    w = lit.WANG21_J0313
    log_l3000 = np.log10(w["lambda_l3000"])
    vo09 = float(masses.virial_mass(w["fwhm_mgii"], log_l3000, "MGII_VO09"))
    s11 = float(masses.virial_mass(w["fwhm_mgii"], log_l3000, "MGII_S11"))
    age = float(physics.age_myr(w["z"]))
    seeds = {}
    for zs in (30.0, 20.0, 15.0):
        dt = age - float(physics.age_myr(zs))
        seeds[int(zs)] = dict(dt_myr=dt, seed_vo09=float(physics.seed_mass(10**vo09, dt)),
                              seed_published=float(physics.seed_mass(w["mbh"], dt)))
    return dict(age_myr=age, log_mbh_vo09_from_wang_width=vo09, log_mbh_s11_from_wang_width=s11,
                log_mbh_published=float(np.log10(w["mbh"])), seeds=seeds)


def jwst_seeds(eps=RADIATIVE_EFFICIENCY, z_start: float = 30.0) -> pd.DataFrame:
    rows = []
    for bh in lit.JWST_BLACK_HOLES:
        age = float(physics.age_myr(bh["z"]))
        dt = age - float(physics.age_myr(z_start))
        rows.append(dict(name=bh["name"], z=bh["z"], log_mbh=bh["log_mbh"], age_myr=age,
                         log_seed_edd=float(np.log10(physics.seed_mass(10 ** bh["log_mbh"], dt, eps))),
                         log_seed_obs=float(np.log10(physics.seed_mass(10 ** bh["log_mbh"], dt, eps,
                                                                       bh["edd_ratio"]))),
                         edd_ratio=bh["edd_ratio"], ref=bh["ref"], method=bh["method"]))
    return pd.DataFrame(rows)


def efficiency_sensitivity(log_m: float, z: float, z_start: float = 30.0,
                           eps_grid=(0.057, 0.1, 0.2, 0.3, 0.42)) -> list[dict]:
    """Seed mass vs radiative efficiency (0.057 non-spinning, 0.42 maximally spinning BH)."""
    dt = float(physics.age_myr(z) - physics.age_myr(z_start))
    return [dict(eps=e, efold_myr=physics.efold_time_myr(e),
                 log_seed=float(np.log10(physics.seed_mass(10 ** log_m, dt, e))))
            for e in eps_grid]
