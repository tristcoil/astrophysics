"""The numbers: what a fall into each of the atlas's black holes would be like.

For every black hole (masses measured in H1, I4, I5; GW250114's remnant and
Gaia BH1 from their papers) this computes, in seconds, kilometres and g:

* sizes: Schwarzschild radius, GM/c^3;
* clocks: time from 10 r_s to the horizon and from the horizon to the
  singularity for a body falling from rest far away (4/3 GM/c^3 inside), the
  longest anyone can last inside (pi GM/c^3, released from rest just above the
  horizon), and how fast a far-away friend sees the faller's light redden
  (e-folding time 4 GM/c^3) and fade (GM/c^3);
* tides: stretch across a 2 m body at the horizon, where it reaches 10 g
  (spaghettification), and how long before the singularity that is;
* hovering: the rocket acceleration needed to hover at the photon sphere;
* spin (Kerr), where measured: horizons, innermost stable orbit, the energy
  released by gas reaching it, and the plunge from there to the horizon.

The view (render.py) depends only on r in units of GM/c^2 and is the same for
all of them; only the clocks and the tides differ.
"""

from __future__ import annotations

import math

import numpy as np

import config as C
import infall
import kerr
import raymap


def final_plunge(a: float, dr: float = 0.1) -> dict:
    """Gas that has drifted inside the innermost stable orbit: from 90% of r_isco to the horizon.

    How long gas lingers near r_isco itself depends on how hard it was nudged (the
    orbit is marginally stable); the last part of the fall does not (checked for
    nudges of 0.001-0.1 M: proper time and turns agree to 1%).
    """
    p = kerr.isco_plunge(a, dr=dr)
    r = p["r"]
    i0 = int(np.flatnonzero(r <= 0.9 * p["r_isco"])[0])
    ih = int(np.flatnonzero(r <= kerr.horizons(a)[0])[0])
    return {"tau": float(p["tau"][ih] - p["tau"][i0]), "v": float(p["v"][ih] - p["v"][i0]),
            "turns": float((p["phi"][ih] - p["phi"][i0]) / (2 * math.pi))}


def per_black_hole() -> dict:
    out = {}
    for key in C.BH_ORDER:
        bh = C.black_holes()[key]
        m, mi, pl = bh["mass"]
        gm = C.GM_SUN * m
        tu = C.gm_over_c3_s(m)                       # GM/c^3
        lu = C.gm_over_c2_m(m)                       # GM/c^2
        rs_m = 2 * lu
        a_10g = C.SURVIVAL_G * C.G0
        r10 = infall.radius_for_stretch(m, C.BODY_LENGTH_M, a_10g, C.GM_SUN)
        r10_M = r10 / lu
        # by the raindrop's clock: time from feeling 10 g to the singularity
        tau_10g = float(infall.rain_tau_left(r10_M)) * tu
        row = {
            "name": bh["name"], "long_name": bh["long_name"], "mass_msun": m, "mass_minus": mi, "mass_plus": pl,
            "mass_source": bh["source"], "note": bh["note"], "spin": bh["spin"],
            "rs_km": rs_m / 1e3, "gm_c3_s": tu,
            "fall_10rs_to_horizon_s": float(infall.rain_tau_left(20.0) - infall.rain_tau_left(2.0)) * tu,
            "horizon_to_singularity_rain_s": 4.0 / 3.0 * tu,
            "horizon_to_singularity_max_s": math.pi * tu,
            "redshift_efold_s": 4.0 * tu,
            "fade_efold_s": tu,
            "stretch_at_horizon_g": float(infall.tidal_stretch_si(rs_m, m, C.BODY_LENGTH_M, C.GM_SUN)) / C.G0,
            "squeeze_at_horizon_g": float(infall.tidal_squeeze_si(rs_m, m, 0.4, C.GM_SUN)) / C.G0,
            "r_10g_km": r10 / 1e3, "r_10g_rs": r10 / rs_m,
            "spaghettified_outside_horizon": bool(r10 > rs_m),
            "time_from_10g_to_singularity_s": tau_10g,
            "hover_accel_at_photon_sphere_g": float(infall.hover_acceleration(3.0)) * C.C_LIGHT**4 / gm / C.G0,
            "hover_accel_at_2rs_g": float(infall.hover_acceleration(4.0)) * C.C_LIGHT**4 / gm / C.G0,
        }
        if bh["spin"] is not None:
            a = float(bh["spin"])
            rp, rm = kerr.horizons(a)
            ri = kerr.isco(a)
            fp = final_plunge(a)
            row["kerr"] = {"a": a, "r_plus_M": rp, "r_minus_M": rm, "r_isco_M": ri,
                           "r_isco_km": ri * lu / 1e3, "efficiency": kerr.efficiency(a),
                           "final_plunge_tau_M": fp["tau"], "final_plunge_tau_s": fp["tau"] * tu,
                           "final_plunge_v_s": fp["v"] * tu, "final_plunge_turns": fp["turns"],
                           "omega_horizon_hz": a / (rp**2 + a * a) / tu / (2 * math.pi),
                           "ergosphere_equator_km": 2.0 * lu / 1e3}
        out[key] = row
    return out


def universal() -> dict:
    """Facts that do not depend on the mass: angles, factors, and two critical masses."""
    rs_ratio = {}
    edge_h = raymap.shadow_edge(2.0)
    _, _, _, g = raymap.photon_constants(2.0, np.array([-1.0, math.cos(edge_h)]))
    m_crit = infall.critical_mass_for_horizon_stretch(C.BODY_LENGTH_M, C.SURVIVAL_G * C.G0, C.GM_SUN, C.C_LIGHT)
    # Hills mass for the Sun: tidal radius R_sun (M/M_sun)^(1/3) equal to the horizon 2GM/c^2
    m_hills = (C.R_SUN * C.C_LIGHT**2 / (2 * C.GM_SUN)) ** 1.5
    for r in (100.0, 20.0, 6.0, 3.0, 2.0, 1.0, 0.1):
        rs_ratio[f"{r / 2:g}"] = {
            "falling_shadow_deg": math.degrees(raymap.shadow_edge(r)),
            "hovering_shadow_deg": (math.degrees(raymap.RayTable(r, raymap.velocity_rel_rain(r, "static")).chi_sh)
                                    if r > 2 else None)}
    return {
        "shadow_at_horizon_deg": math.degrees(edge_h),
        "g_from_straight_behind_at_horizon": float(g[0]),
        "g_at_shadow_edge_at_horizon": float(g[1]),
        "shadow_by_radius": rs_ratio,
        "proper_time_inside_rain_M": 4.0 / 3.0,
        "proper_time_inside_max_M": math.pi,
        "critical_mass_10g_at_horizon_msun": m_crit,
        "hills_mass_sun_msun": m_hills,
    }


def summary() -> dict:
    return {"black_holes": per_black_hole(), "universal": universal()}


if __name__ == "__main__":
    import json
    print(json.dumps(summary(), indent=1, default=float))
