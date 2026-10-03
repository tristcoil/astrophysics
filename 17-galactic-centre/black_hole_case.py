"""Why the mass S2 orbits must be a black hole.

S2's orbit weighs about four million suns and shows they sit inside its
pericentre, ~120 AU from the centre. Could that be a cluster of dark objects
(brown dwarfs, white dwarfs, neutron stars, stellar black holes) instead of
one black hole? Maoz (1998, ApJ 494, L181; astro-ph/9710309) turned the
question into arithmetic, which we redo here with S2's numbers:

* A cluster of N objects **evaporates**: weak encounters keep kicking members
  above escape speed. For equal masses this takes ~300 median relaxation
  times (Spitzer & Thuan 1972), and the relaxation time *shrinks* as the
  objects get heavier (fewer, bigger bodies scatter each other harder).
* A cluster also destroys itself by **physical collisions**, which win for
  light, puffy objects (planets, brown dwarfs), boosted by gravitational
  focusing.

Each kind of object must lose to one process or the other. Maoz made every
assumption as generous to the cluster as possible: equal masses (no mass
segregation), the least concentrated profile (Plummer), zero-temperature
(smallest possible) objects, and he ignored core collapse, which would shorten
lifetimes ~20x further. The longest lifetime over all object types is the
cluster's best case; if even that is far below the Galaxy's ~10 Gyr age, a
cluster would be caught in a fleeting moment, which is implausible.

Following Maoz's Table 1, the cluster's half-mass is half the measured mass,
and its half-mass radius is at most the radius the data probe -- here S2's
pericentre. Two exotic options survive any such argument (Maoz sect. 3): a
cluster of very light *primordial* black holes (evaporation slows as the
constituents get lighter, and they cannot collide), and a boson star; both
need physics nobody has seen.
"""

from __future__ import annotations

import numpy as np

import orbit_physics as op

G_PC3_MSUN_YR2 = 4.498_502e-15  # G in pc^3 / (M_sun yr^2)
R_SUN_CM = 6.957e10
CM_PER_PC = op.PC_M * 100.0
KM_PER_PC = op.PC_M / 1e3
SEC_PER_YR = op.YEAR_S
GALAXY_AGE_YR = 1.0e10

# Nauenberg (1972) mass-radius relation for zero-temperature white dwarfs (mu_e = 2).
MU_E = 2.0
M3_CHANDRA = 5.816 / MU_E**2  # Chandrasekhar mass, M_sun

#: Mass ranges of each class of dark object, M_sun (Maoz 1998, section 2.1).
CLASSES = {
    "planets": (1e-5, 3e-3),
    "brown dwarfs": (3e-3, 0.09),
    "white dwarfs": (3e-3, 1.4),
    "neutron stars": (1.4, 3.0),
    "stellar black holes": (3.0, 100.0),
}
NEUTRON_STAR_RADIUS_KM = 12.0


def radius_cm(cls: str, m):
    """Radius of a cold object of mass m (M_sun), cm."""
    m = np.asarray(m, dtype=float)
    if cls in ("planets", "brown dwarfs"):
        # Zapolsky & Salpeter (1969) for cold bodies of cosmic composition.
        return 2.2e9 * m ** (-1.0 / 3.0) * (1.0 + (m / 0.0032) ** -0.5) ** (-4.0 / 3.0)
    if cls == "white dwarfs":
        x = m / M3_CHANDRA
        return 1.57e9 * np.sqrt(np.clip(1.0 - x ** (4.0 / 3.0), 0.0, None)) / (MU_E * x ** (1.0 / 3.0))
    if cls == "neutron stars":
        return np.full_like(m, NEUTRON_STAR_RADIUS_KM * 1e5)
    if cls in ("stellar black holes", "primordial black holes"):
        return 2.0 * m * op.GM_SUN_KM3S2 / op.C_KMS**2 * 1e5
    raise KeyError(cls)


def evaporation_time_yr(m_half, rho_half, m_star):
    """Maoz eq. (4): t_evap = 300 t_rh for a Plummer cluster of equal masses."""
    n_half = np.asarray(m_half) / np.asarray(m_star)
    return 4.3e4 * n_half / np.log(0.8 * n_half) * (np.asarray(rho_half) / 1e8) ** -0.5


def collision_time_yr(m_half, rho_half, m_star, r_star_cm):
    """Maoz eq. (5): t_coll = [16 sqrt(pi) n sigma r*^2 (1 + G m*/(2 sigma^2 r*))]^-1.

    With the Plummer central dispersion sigma_0 = 2^(-1/4) G^(1/2) Mh^(1/3) rho_h^(1/6)
    and n = rho_h / m*.
    """
    sigma = 2.0**-0.25 * np.sqrt(G_PC3_MSUN_YR2) * np.asarray(m_half) ** (1 / 3) * np.asarray(rho_half) ** (1 / 6)
    n = np.asarray(rho_half) / np.asarray(m_star)
    r = np.asarray(r_star_cm) / CM_PER_PC
    focus = 1.0 + G_PC3_MSUN_YR2 * np.asarray(m_star) / (2.0 * sigma**2 * r)
    return 1.0 / (16.0 * np.sqrt(np.pi) * n * sigma * r**2 * focus)


def lifetime_curves(m_half, rho_half, *, n_grid: int = 400) -> dict:
    """Lifetime min(t_evap, t_coll) against constituent mass, for each class."""
    out = {}
    for cls, (lo, hi) in CLASSES.items():
        m = np.logspace(np.log10(lo), np.log10(hi), n_grid)
        t_ev = evaporation_time_yr(m_half, rho_half, m)
        if cls == "stellar black holes":
            t_co = np.full_like(m, np.inf)  # black holes cannot collide physically
        else:
            t_co = collision_time_yr(m_half, rho_half, m, radius_cm(cls, m))
        out[cls] = {"m": m, "t_evap": t_ev, "t_coll": t_co, "lifetime": np.minimum(t_ev, t_co)}
    # The exotic loophole: black holes too light to form from stars.
    m = np.logspace(-8, np.log10(3.0), n_grid)
    out["primordial black holes"] = {"m": m, "t_evap": evaporation_time_yr(m_half, rho_half, m),
                                     "t_coll": np.full_like(m, np.inf)}
    out["primordial black holes"]["lifetime"] = out["primordial black holes"]["t_evap"]
    return out


def max_lifetime(m_half, rho_half) -> dict:
    """The cluster's best case over every astrophysical kind of dark object."""
    curves = lifetime_curves(m_half, rho_half)
    best = {"lifetime_yr": 0.0}
    for cls, c in curves.items():
        if cls == "primordial black holes":
            continue
        k = int(np.argmax(c["lifetime"]))
        if c["lifetime"][k] > best["lifetime_yr"]:
            best = {"lifetime_yr": float(c["lifetime"][k]), "class": cls, "m_star": float(c["m"][k])}
    return best


def primordial_bh_mass_limit(m_half, rho_half, age_yr=GALAXY_AGE_YR) -> float:
    """Heaviest primordial black hole whose cluster would survive ``age_yr``."""
    m = np.logspace(-12, 1, 4000)
    t = evaporation_time_yr(m_half, rho_half, m)
    ok = t > age_yr
    return float(m[ok].max()) if np.any(ok) else float("nan")


def density_bound_msun_pc3(gm6, r_au):
    """Mean density if the whole mass sits inside radius r: M / (4/3 pi r^3)."""
    r_pc = np.asarray(r_au) / op.PC_AU
    return np.asarray(gm6) * 1e6 / (4.0 / 3.0 * np.pi * r_pc**3)


def density_needed_for_lifetime(m_half, age_yr=GALAXY_AGE_YR) -> float:
    """The half-mass density below which some dark cluster could last ``age_yr``."""
    lo, hi = 1e3, 1e22
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        if max_lifetime(m_half, mid)["lifetime_yr"] > age_yr:
            lo = mid
        else:
            hi = mid
    return float(np.sqrt(lo * hi))


#: Maoz (1998) Table 1: the Galactic Centre as known in 1997 (Eckart & Genzel
#: 1997), before any star's orbit had been traced.
MAOZ_1997_GALAXY = {"m_half": 1.3e6, "r_half_pc": 0.008, "rho_half": 6e11}


def case_for_black_hole(joint_fit: dict, names: list[str]) -> dict:
    """The dark-cluster test with our posterior for M and S2's pericentre."""
    chain = joint_fit["chain"]
    p = {n: chain[:, k] for k, n in enumerate(names)}
    a = op.semi_major_axis_au(p["gm6"], p["period"])
    rp = a * (1.0 - p["e"])
    mass = p["gm6"] * 1e6
    # Maoz's convention: half the mass inside a half-mass radius no larger than
    # the region probed -- here S2's pericentre.
    m_half = 0.5 * mass
    r_half_pc = rp / op.PC_AU
    rho_half = 3.0 * m_half / (4.0 * np.pi * r_half_pc**3)
    mh, rh = float(np.median(m_half)), float(np.median(rho_half))
    now = max_lifetime(mh, rh)
    then = max_lifetime(MAOZ_1997_GALAXY["m_half"], MAOZ_1997_GALAXY["rho_half"])
    # Robustness: lower 2.5% of density (the most generous to the cluster).
    rho_lo = float(np.percentile(rho_half, 2.5))
    now_lo = max_lifetime(mh, rho_lo)
    rs_au = op.schwarzschild_radius_au(np.median(p["gm6"]))
    return {
        "m_half_msun": mh,
        "r_half_pc": float(np.median(r_half_pc)),
        "r_half_au": float(np.median(rp)),
        "rho_half_msun_pc3": rh,
        "rho_half_lo95": rho_lo,
        "rho_mean_inside_pericentre": float(np.median(density_bound_msun_pc3(p["gm6"], rp))),
        "max_lifetime_now": now,
        "max_lifetime_now_generous": now_lo,
        "max_lifetime_1997": then,
        "maoz_1997_input": MAOZ_1997_GALAXY,
        "density_gain_since_1997": rh / MAOZ_1997_GALAXY["rho_half"],
        "density_needed_for_10_gyr": density_needed_for_lifetime(mh),
        "primordial_bh_mass_limit_msun": primordial_bh_mass_limit(mh, rh),
        "primordial_bh_mass_limit_1997_msun": primordial_bh_mass_limit(
            MAOZ_1997_GALAXY["m_half"], MAOZ_1997_GALAXY["rho_half"]),
        "schwarzschild_radius_au": float(rs_au),
        "schwarzschild_radius_km": float(rs_au * op.AU_KM),
        "pericentre_over_rs": float(np.median(rp) / rs_au),
    }
