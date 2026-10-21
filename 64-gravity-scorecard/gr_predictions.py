"""General relativity's prediction for every classic test, computed here from first principles.

Each function returns the prediction with the inputs it used, so the scorecard
can show "GR says X, the experiment measured Y". Where an experiment's own GR
prediction depends on details we do not model (Gravity Probe B's exact orbit and
guide star, a binary pulsar's full timing model), the published prediction is
used for the ratio and ours is shown beside it as a check.

Constants: IAU 2015 nominal solar values, IERS 2010 Earth values, DE440 GMs.
"""

from __future__ import annotations

import numpy as np

import config as C

G_NEWTON = 6.67430e-11                 # CODATA 2018
GM_SUN = C.GM_SUN * 1e9                # m^3/s^2
GM_EARTH = 3.986004418e14              # m^3/s^2 (IERS 2010)
R_SUN = 6.957e8                        # m, IAU 2015 nominal radius
R_EARTH = 6_378_136.6                  # m
OMEGA_EARTH = 7.292115e-5              # rad/s
EARTH_C_MR2 = 0.3307                   # Earth's polar moment of inertia, C/(M R^2)
M_SUN_M = GM_SUN / C.C_LIGHT**2        # GM_sun/c^2 = 1476.6 m
ARCSEC = np.pi / 180 / 3600
MAS = ARCSEC / 1e3
YEAR = C.JULIAN_YEAR
CENTURY = 100 * YEAR


def mercury_perihelion() -> dict:
    """6 pi GM / (c^2 a (1 - e^2)) per orbit, in arcsec per Julian century."""
    a = 0.38709927 * C.AU          # J2000 mean elements (Standish & Williams, Expl. Suppl. 3rd ed.)
    e = 0.20563593
    period = 87.9691 * C.DAY
    per_orbit = 6 * np.pi * M_SUN_M / (a * (1 - e**2))
    rate = per_orbit / period * CENTURY / ARCSEC
    return {"value": rate, "unit": "arcsec/century", "formula": "6πGM/(c²a(1−e²)) per orbit",
            "ppn_factor": "(2 + 2γ − β)/3", "inputs": {"a_AU": 0.38709927, "e": e, "P_days": 87.9691}}


def planet_perihelia() -> dict:
    """The same for Venus, Earth, Mars (arcsec/century)."""
    el = {"Venus": (0.72333566, 0.00677672, 224.701), "Earth": (1.00000261, 0.01671123, 365.256),
          "Mars": (1.52371034, 0.09339410, 686.980)}
    out = {}
    for k, (a, e, p) in el.items():
        per = 6 * np.pi * M_SUN_M / (a * C.AU * (1 - e**2))
        out[k] = per / (p * C.DAY) * CENTURY / ARCSEC
    return out


def light_deflection_limb() -> dict:
    """(1+gamma)/2 * 4GM/(c^2 R): the deflection of a ray grazing the Sun."""
    val = 4 * M_SUN_M / R_SUN / ARCSEC
    return {"value": val, "unit": "arcsec", "formula": "4GM/(c²R☉)", "ppn_factor": "(1 + γ)/2",
            "newton_half": val / 2}


def shapiro_mars_conjunction() -> dict:
    """Round-trip radar delay for a signal grazing the Sun to Mars at superior conjunction."""
    r_e, r_m, b = 1.0 * C.AU, 1.52 * C.AU, R_SUN
    dt = 2 * (2 * M_SUN_M / C.C_LIGHT) * np.log(4 * r_e * r_m / b**2)
    return {"value": dt * 1e6, "unit": "microseconds (round trip)", "formula": "2(1+γ)GM/c³ ln(4 r_E r_P / b²)",
            "ppn_factor": "(1 + γ)/2"}


def cassini_doppler_peak() -> dict:
    """Peak fractional frequency shift of a radio link passing the Sun (Bertotti et al. 2003):
    y = 2(1 + gamma) (GM/c^3) (db/dt)/b, at b = 1.6 R_sun with db/dt ~ Earth's orbital speed."""
    b = 1.6 * R_SUN
    bdot = 29.8e3
    y = 2 * 2 * (M_SUN_M / C.C_LIGHT) * bdot / b
    return {"value": y, "unit": "fractional frequency", "formula": "2(1+γ)(GM/c³)(db/dt)/b", "ppn_factor": "(1 + γ)/2"}


def redshift_tower(height_m: float = 22.6) -> dict:
    g = 9.80
    return {"value": g * height_m / C.C_LIGHT**2, "unit": "fractional", "formula": "gh/c²", "inputs": {"h_m": height_m}}


def redshift_sun() -> dict:
    z = M_SUN_M / R_SUN
    earth = M_SUN_M / C.AU
    return {"value": z * C.C_LIGHT / 1e3, "seen_from_earth_kms": (z - earth) * C.C_LIGHT / 1e3,
            "unit": "km/s", "formula": "GM/(Rc²) × c"}


def gpb() -> dict:
    """Gravity Probe B: orbit-averaged precessions for an ideal circular polar orbit at 642 km.

    Geodetic: (3/2) (GM/(c^2 r)) n.  Frame dragging for a gyro pointing along the
    orbit plane's normal-perpendicular direction averages to G I w / (2 c^2 r^3).
    (The published GR values, 6606.1 and 39.2 mas/yr, use the real orbit, the
    Earth's oblateness and IM Pegasi's direction.)
    """
    r = R_EARTH + 642e3
    n = np.sqrt(GM_EARTH / r**3)
    geo = 1.5 * GM_EARTH / (C.C_LIGHT**2 * r) * n
    lt = EARTH_C_MR2 * GM_EARTH * R_EARTH**2 * OMEGA_EARTH / (2 * C.C_LIGHT**2 * r**3)
    return {"geodetic_mas_yr": geo * YEAR / MAS, "frame_dragging_mas_yr": lt * YEAR / MAS,
            "published_GR": {"geodetic": 6606.1, "frame_dragging": 39.2}, "ppn_factors": {"geodetic": "(1 + 2γ)/3",
                                                                                         "frame_dragging": "(1 + γ + α1/4)/2"}}


def lageos_node_lense_thirring(a_m: float = 12_270_020.0, e: float = 0.0) -> dict:
    """Lense-Thirring drift of an orbit's node: 2 G J / (c^2 a^3 (1 - e^2)^1.5) (rad/s)."""
    j = EARTH_C_MR2 * GM_EARTH * R_EARTH**2 * OMEGA_EARTH
    rate = 2 * j / (C.C_LIGHT**2 * a_m**3 * (1 - e**2) ** 1.5)
    return {"value": rate * YEAR / MAS, "unit": "mas/yr", "formula": "2GJ/(c²a³(1−e²)^{3/2})"}


def lunar_geodetic_precession() -> dict:
    """de Sitter precession of the Earth-Moon system around the Sun: (3/2)(GM/(c^2 a)) n."""
    n = 2 * np.pi / (365.25636 * C.DAY)
    rate = 1.5 * M_SUN_M / C.AU * n
    return {"value": rate * YEAR / MAS, "unit": "mas/yr", "formula": "(3/2)(GM☉/(c²a))n"}


def nordtvedt_amplitude(sim_coefficient_m: float) -> dict:
    return {"value": sim_coefficient_m, "unit": "m per unit eta (cos D)", "GR": 0.0,
            "literature": {"Nordtvedt 1995 / Will 2014": 13.1}}


def schwarzschild_shadow(mass_msun: float, dist_m: float) -> dict:
    """Shadow (photon-capture) diameter 2*sqrt(27) GM/(c^2 D) in micro-arcseconds."""
    theta_g = mass_msun * M_SUN_M / dist_m
    return {"theta_g_uas": theta_g / ARCSEC * 1e6, "shadow_uas": 2 * np.sqrt(27) * theta_g / ARCSEC * 1e6}


def s2_precession(mass_msun: float = 4.30e6, a_au: float = 1031.0, e: float = 0.884) -> dict:
    """Schwarzschild precession of S2 per orbit (arcmin)."""
    per = 6 * np.pi * mass_msun * M_SUN_M / (a_au * C.AU * (1 - e**2))
    return {"value": per / ARCSEC / 60, "unit": "arcmin per orbit", "inputs": {"M": mass_msun, "a_AU": a_au, "e": e}}


# ---------------------------------------------------------------------------
# Binary pulsars: post-Keplerian parameters from the two masses (Damour & Deruelle;
# Taylor & Weisberg 1989; Lorimer & Kramer 2005, eqs. 8.48-8.52)
# ---------------------------------------------------------------------------

T_SUN = GM_SUN / C.C_LIGHT**3          # 4.925490947 us


def pk_parameters(mp, mc, pb_days, e, x_lts=None):
    """GR's post-Keplerian parameters for pulsar mass mp and companion mc (solar masses).

    Returns omega_dot (deg/yr), gamma (s), pbdot (dimensionless), r (s), s (sin i).
    """
    pb = pb_days * C.DAY
    M = mp + mc
    nb = 2 * np.pi / pb
    omdot = 3 * (pb / (2 * np.pi)) ** (-5 / 3) * (T_SUN * M) ** (2 / 3) / (1 - e**2)            # rad/s
    gam = (pb / (2 * np.pi)) ** (1 / 3) * e * T_SUN ** (2 / 3) * M ** (-4 / 3) * mc * (mp + 2 * mc)
    fe = (1 + 73 / 24 * e**2 + 37 / 96 * e**4) / (1 - e**2) ** 3.5
    pbdot = -192 * np.pi / 5 * (pb / (2 * np.pi)) ** (-5 / 3) * fe * T_SUN ** (5 / 3) * mp * mc * M ** (-1 / 3)
    r = T_SUN * mc
    out = {"omega_dot_deg_yr": np.degrees(omdot) * YEAR, "gamma_s": gam, "pbdot": pbdot, "r_s": r}
    if x_lts is not None:
        out["s"] = x_lts * (pb / (2 * np.pi)) ** (-2 / 3) * T_SUN ** (-1 / 3) * M ** (2 / 3) / mc
    return out


def table() -> list[dict]:
    """Every classic prediction with a short label (used by the README and scorecard)."""
    rows = []
    m = mercury_perihelion()
    rows.append({"test": "Mercury's perihelion advance", "gr": m["value"], "unit": m["unit"], "formula": m["formula"]})
    ld = light_deflection_limb()
    rows.append({"test": "Starlight bent at the Sun's limb", "gr": ld["value"], "unit": ld["unit"], "formula": ld["formula"]})
    rows.append({"test": "  (Newton/Soldner, 'half' value)", "gr": ld["newton_half"], "unit": "arcsec", "formula": "2GM/(c²R☉)"})
    sd = shapiro_mars_conjunction()
    rows.append({"test": "Radar echo delay, Mars behind the Sun", "gr": sd["value"], "unit": sd["unit"], "formula": sd["formula"]})
    cd = cassini_doppler_peak()
    rows.append({"test": "Cassini radio Doppler, peak", "gr": cd["value"], "unit": cd["unit"], "formula": cd["formula"]})
    rt = redshift_tower()
    rows.append({"test": "Pound-Rebka redshift (22.6 m)", "gr": rt["value"], "unit": rt["unit"], "formula": rt["formula"]})
    rs = redshift_sun()
    rows.append({"test": "Sun's gravitational redshift (from Earth)", "gr": rs["seen_from_earth_kms"] * 1e3, "unit": "m/s",
                 "formula": "(GM/R − GM/AU)/c"})
    g = gpb()
    rows.append({"test": "Gravity Probe B geodetic (ideal orbit)", "gr": g["geodetic_mas_yr"], "unit": "mas/yr",
                 "formula": "(3/2)(GM/c²r)n"})
    rows.append({"test": "Gravity Probe B frame dragging (ideal orbit)", "gr": g["frame_dragging_mas_yr"], "unit": "mas/yr",
                 "formula": "GIω/(2c²r³)"})
    ln = lageos_node_lense_thirring()
    rows.append({"test": "LAGEOS node, Lense-Thirring", "gr": ln["value"], "unit": "mas/yr", "formula": ln["formula"]})
    lg = lunar_geodetic_precession()
    rows.append({"test": "Lunar orbit geodetic precession", "gr": lg["value"], "unit": "mas/yr", "formula": lg["formula"]})
    s2 = s2_precession()
    rows.append({"test": "S2 pericentre precession (Sgr A*)", "gr": s2["value"], "unit": s2["unit"], "formula": "6πGM/(c²a(1−e²))"})
    m87 = schwarzschild_shadow(6.5e9, 16.8e6 * 3.0856775814913673e16)
    rows.append({"test": "M87* shadow, 6.5e9 Msun at 16.8 Mpc", "gr": m87["shadow_uas"], "unit": "µas",
                 "formula": "2√27 GM/(c²D)"})
    sgr = schwarzschild_shadow(4.0e6, 8.0e3 * 3.0856775814913673e16)
    rows.append({"test": "Sgr A* shadow, 4.0e6 Msun at 8.0 kpc", "gr": sgr["shadow_uas"], "unit": "µas",
                 "formula": "2√27 GM/(c²D)"})
    return rows


if __name__ == "__main__":
    for r in table():
        print(f"{r['test']:48s} {r['gr']:12.6g} {r['unit']:28s} {r['formula']}")
    print("Venus/Earth/Mars perihelia:", {k: round(v, 3) for k, v in planet_perihelia().items()})
    # Hulse-Taylor with Weisberg & Huang (2016) masses as a check of the formulae
    print("B1913+16 check:", pk_parameters(1.438, 1.390, 0.322997448918, 0.6171340))
