"""Mercury's perihelion budget: what turns the orbit, and by how much.

The method is subtraction. Each N-body run switches on a different set of
forces; Mercury's in-plane perihelion angle ``psi`` is measured in every run
and fitted with a slowly changing rate (a quadratic in time over 1550-2650,
whose linear coefficient is the rate at J2000). Then

* a planet's share   = the rate with only that planet tugging (Sun + Mercury + it);
* the planets' total = the rate with all of them (and how far that is from the
  sum of the shares measures how non-additive the tugs are);
* GR's share         = (all planets + GR) - (all planets), fitted on the
  difference of the two series, in which the planets' wiggles cancel;
* the Sun's J2 and spin likewise, from Sun + Mercury runs.

The answer key is DE440 itself, analysed in exactly the same way.
"""

from __future__ import annotations

import numpy as np

import ephemeris as eph
import nbody
import orbit_physics as op
from config import (ARCSEC_PER_RAD, BODY, CENTURY_DAYS, J2000_JD, PERTURBER_KEYS, RUN, RUNS,
                    SUN_POLE_DEC_DEG, SUN_POLE_RA_DEG, WHFAST_RUNS)

MU_MERCURY = None  # set lazily (needs the header)


def mu(key: str) -> float:
    return eph.gm("sun") + eph.gm(key)


def mercury_series(t_days: np.ndarray, states: np.ndarray) -> dict:
    """Angle series for Mercury from (N, 6) heliocentric states, three ways."""
    r, v = states[:, :3].T, states[:, 3:].T
    ecl = op.angle_series(r, v, mu("mercury"), eph.rot_icrf_to_ecliptic())
    eq = op.elements(r, v, mu("mercury"))
    return {"t_cy": t_days / CENTURY_DAYS, "psi": ecl["psi"], "varpi_ecl": ecl["varpi"],
            "Omega_ecl": ecl["Omega"], "inc_ecl": ecl["inc"], "varpi_icrf": eq["varpi"],
            "Omega_icrf": eq["Omega"], "e": ecl["e"], "a": ecl["a"], "ehat": eq["ehat"],
            "hhat": eq["hhat"]}


def rates(series: dict) -> dict:
    t = series["t_cy"]
    out = {k: op.fit_secular(t, series[k]) for k in ("psi", "varpi_ecl", "varpi_icrf", "Omega_ecl")}
    out["psi_stability"] = op.rate_stability(t, series["psi"])
    return out


def de440_series(t_days: np.ndarray) -> dict:
    jd = J2000_JD + t_days
    r, v = eph.heliocentric_state("mercury", jd)
    return mercury_series(t_days, np.vstack([r, v]).T)


def run_series(name: str) -> dict:
    rec = nbody.load(name)
    return mercury_series(rec["t"], rec["mercury"])


def difference_rate(a: dict, b: dict) -> dict:
    """Rate of psi_a - psi_b (same sample times): one force's share, cleanly."""
    if not np.allclose(a["t_cy"], b["t_cy"]):
        raise ValueError("runs sampled at different times")
    fit = op.fit_secular(a["t_cy"], a["psi"] - b["psi"])
    return {"rate": fit["rate"], "curvature": fit["curvature"],
            "rms_resid_arcsec": fit["rms_resid_arcsec"]}


def direction_offset(a: dict, b: dict) -> np.ndarray:
    """Epoch-by-epoch angle (rad) from b's perihelion direction to a's, about b's normal."""
    cross = np.cross(b["ehat"].T, a["ehat"].T).T
    return np.arctan2((cross * b["hhat"]).sum(0), (a["ehat"] * b["ehat"]).sum(0))


# --------------------------------------------------------------------------
# Analytic expectations
# --------------------------------------------------------------------------


def j2000_elements(key: str) -> dict:
    r, v = eph.heliocentric_state(key, J2000_JD)
    el = op.elements(r, v, mu(key))
    ecl = op.elements(eph.rot_icrf_to_ecliptic() @ r, eph.rot_icrf_to_ecliptic() @ v, mu(key))
    pole = eph.unit_vector_radec(SUN_POLE_RA_DEG, SUN_POLE_DEC_DEG)
    inc_sun = float(np.arccos(np.clip(el["hhat"][:, 0] @ pole, -1, 1)))
    return {"a": float(el["a"][0]), "e": float(el["e"][0]), "varpi_ecl": float(ecl["varpi"][0]),
            "inc_ecl": float(ecl["inc"][0]), "inc_sun_equator": inc_sun}


def analytic_expectations() -> dict:
    c = eph.c_au_per_day()
    m = j2000_elements("mercury")
    gm_s = eph.gm("sun")
    phys = eph.sun_physical()
    g_spin = gm_s * phys["moment_of_inertia_factor"] * phys["radius_au"] ** 2 * phys["omega_rad_per_day"]
    out = {
        "mercury_j2000": {**m, "inc_ecl_deg": np.degrees(m["inc_ecl"]),
                          "inc_sun_equator_deg": np.degrees(m["inc_sun_equator"])},
        "schwarzschild": op.schwarzschild_rate(mu("mercury"), m["a"], m["e"], c),
        "schwarzschild_per_orbit_arcsec": op.schwarzschild_per_orbit(mu("mercury"), m["a"], m["e"], c)
        * ARCSEC_PER_RAD,
        "j2": op.j2_in_plane_rate(mu("mercury"), m["a"], m["e"], phys["j2"], phys["radius_au"],
                                  m["inc_sun_equator"]),
        "lense_thirring": op.lense_thirring_in_plane_rate(g_spin, m["a"], m["e"], c,
                                                          m["inc_sun_equator"]),
        "sun_angular_momentum_kg_m2_s": g_spin / 6.67430e-11 * (eph.au_km() * 1e3) ** 5 / 86400.0**3,
        "laplace_lagrange": {},
    }
    for k in PERTURBER_KEYS:
        p = j2000_elements(k)
        out["laplace_lagrange"][k] = op.laplace_lagrange_rate(
            gm_s, eph.gm("mercury"), m["a"], m["e"], m["varpi_ecl"], eph.gm(k), p["a"], p["e"],
            p["varpi_ecl"])
    out["laplace_lagrange_total"] = sum(v["total"] for v in out["laplace_lagrange"].values())
    # Gauss's method: exact in e and i, first order in the masses.
    r_m, v_m = eph.heliocentric_state("mercury", J2000_JD)
    out["gauss"] = {}
    for k in PERTURBER_KEYS:
        r_p, v_p = eph.heliocentric_state(k, J2000_JD)
        out["gauss"][k] = op.gauss_secular_rate(mu("mercury"), r_m, v_m, eph.gm(k), r_p, v_p)["rate"]
    out["gauss_total"] = sum(out["gauss"].values())
    # GR for every planet, for the N-body comparison.
    out["schwarzschild_by_planet"] = {}
    for k in ["mercury", *PERTURBER_KEYS]:
        p = j2000_elements(k)
        out["schwarzschild_by_planet"][k] = op.schwarzschild_rate(mu(k), p["a"], p["e"], c)
    return out


def asteroid_belt_estimate() -> dict:
    """The 343 main-belt asteroids DE440 includes (and our runs leave out).

    Each acts on Mercury, to first order, like a ring along its orbit: Mercury's
    perihelion turns at n/4 (m/M) alpha^2 b_{3/2}^{(1)}(alpha) per asteroid
    (Laplace-Lagrange's free term; alpha = a_mercury / a_asteroid). Masses from
    DE440's header (MA0001-MA1467), orbits from JPL's small-body database.
    """
    import json
    from fetch_data import SBDB_ASTEROIDS

    c = eph.constants()
    rows = json.loads(SBDB_ASTEROIDS.read_text())["data"]
    a_by_number = {int(r[1]): float(r[2]) for r in rows}
    m = j2000_elements("mercury")
    n = op.mean_motion(mu("mercury"), m["a"])
    total, biggest = 0.0, []
    for key, gm in c.items():
        if not key.startswith("MA") or int(key[2:]) >= 8000:
            continue
        num = int(key[2:])
        alpha = m["a"] / a_by_number[num]
        rate = op.per_century(0.25 * n * gm / mu("mercury") * alpha**2 *
                              op.laplace_coefficient(1.5, 1, alpha))
        total += rate
        biggest.append((rate, num))
    biggest.sort(reverse=True)
    kbo_gm = sum(v for k, v in c.items() if k.startswith("MA") and int(k[2:]) >= 8000)
    kbo = op.per_century(0.75 * n * kbo_gm / mu("mercury") * (m["a"] / 44.0) ** 3)
    return {"main_belt_rate": total, "n_asteroids": len(biggest),
            "top5": [{"number": num, "rate": r} for r, num in biggest[:5]],
            "kuiper_belt_rate_at_44au": kbo,
            "main_belt_mass_msun": sum(v for k, v in c.items()
                                       if k.startswith("MA") and int(k[2:]) < 8000) / c["GMS"]}


# --------------------------------------------------------------------------
# Other planets: GR's share from the coarse samples
# --------------------------------------------------------------------------


def planet_psi(name: str, key: str) -> tuple[np.ndarray, np.ndarray]:
    rec = nbody.load(name)
    i = rec["meta"]["bodies"].index(key)
    s = rec["outer"][:, i, :]
    el = op.angle_series(s[:, :3].T, s[:, 3:].T, mu(key))
    return rec["t_outer"] / CENTURY_DAYS, el["psi"]


def gr_by_planet() -> dict:
    out = {}
    for k in ["mercury", *PERTURBER_KEYS]:
        t1, p1 = planet_psi("all_gr_full", k)
        t0, p0 = planet_psi("all_newton", k)
        fit = op.fit_secular(t1, p1 - p0)
        out[k] = {"rate": fit["rate"], "rms_resid_arcsec": fit["rms_resid_arcsec"]}
    return out


# --------------------------------------------------------------------------
# The whole budget
# --------------------------------------------------------------------------


def analyse() -> dict:
    series = {r.name: run_series(r.name) for r in RUNS}
    t_days = nbody.load("all_newton")["t"]
    de = de440_series(t_days)

    res = {"definition": ("psi = turn of Mercury's eccentricity vector about its orbit normal "
                          "(frame-independent); rates in arcsec per Julian century at J2000 from "
                          "a quadratic fit over 1550-2650"),
           "de440": rates(de), "runs": {}}
    for name, s in series.items():
        res["runs"][name] = rates(s)

    two0 = series["two_body_newton"]
    contrib = {"planets": {}, "gr_two_body": {}}
    for k in PERTURBER_KEYS:
        contrib["planets"][k] = difference_rate(series[f"pair_{k}"], two0)
    contrib["planets_sum"] = sum(v["rate"] for v in contrib["planets"].values())
    contrib["all_newton"] = difference_rate(series["all_newton"], two0)
    contrib["nonadditivity"] = contrib["all_newton"]["rate"] - contrib["planets_sum"]
    for g in ("gr", "gr_full", "gr_potential"):
        contrib["gr_two_body"][g] = difference_rate(series[f"two_body_{g}"], two0)
    contrib["gr_in_system"] = {
        "gr_full": difference_rate(series["all_gr_full"], series["all_newton"]),
        "gr": difference_rate(series["all_gr"], series["all_newton"]),
    }
    contrib["j2"] = difference_rate(series["two_body_j2"], two0)
    contrib["lense_thirring"] = difference_rate(series["two_body_lt"], two0)
    contrib["j2_plus_lt_in_system"] = difference_rate(series["full_model"], series["all_gr_full"])
    contrib["full_model_total"] = res["runs"]["full_model"]["psi"]["rate"]
    contrib["de440_total"] = res["de440"]["psi"]["rate"]
    contrib["full_model_minus_de440"] = difference_rate(series["full_model"], de)
    contrib["numerical_floor"] = res["runs"]["two_body_newton"]["psi"]["rate"]
    res["contributions"] = contrib

    # The dogleg: varpi in the J2000 ecliptic vs the in-plane angle.
    for label, s in (("de440", de), ("all_newton", series["all_newton"])):
        r = res["de440"] if label == "de440" else res["runs"][label]
        res.setdefault("frames", {})[label] = {
            "psi": r["psi"]["rate"], "varpi_ecliptic_j2000": r["varpi_ecl"]["rate"],
            "varpi_icrf_equator": r["varpi_icrf"]["rate"], "Omega_ecliptic": r["Omega_ecl"]["rate"],
            "mean_inc_ecliptic_deg": float(np.degrees(s["inc_ecl"].mean())),
            "dogleg_predicted": float(r["Omega_ecl"]["rate"] * (1 - np.cos(s["inc_ecl"].mean()))),
        }

    # Epoch-by-epoch comparison with DE440.
    res["vs_de440"] = {}
    for name in ("all_newton", "all_gr", "all_gr_full", "full_model"):
        d = direction_offset(series[name], de) * ARCSEC_PER_RAD
        rec = nbody.load(name)
        r_de, _ = eph.heliocentric_state("mercury", J2000_JD + rec["t"])
        dpos = np.linalg.norm(rec["mercury"][:, :3].T - r_de, axis=0) * eph.au_km()
        ends = {"1550": float(d[0]), "2650": float(d[-1])}
        res["vs_de440"][name] = {"perihelion_offset_arcsec_at": ends,
                                 "max_abs_perihelion_offset_arcsec": float(np.abs(d).max()),
                                 "position_offset_km_at": {"1550": float(dpos[0]),
                                                           "2650": float(dpos[-1])},
                                 "max_position_offset_km": float(dpos.max())}
    res["analytic"] = analytic_expectations()
    res["gr_by_planet"] = gr_by_planet()
    return res


def whfast_test() -> dict:
    """Spurious precession from a symplectic integrator's step, against IAS15.

    Both series are fitted at *identical* times (the IAS15 run is sampled every
    2 days; a WHFast run with a longer step is compared on the IAS15 samples at
    its own output times), so that sampling cannot masquerade as integrator error.
    """
    ref = run_series("ias15_200yr")
    out = {"ias15_rate": op.fit_secular(ref["t_cy"], ref["psi"])["rate"], "runs": {}}
    for spec in WHFAST_RUNS:
        if spec.integrator != "whfast":
            continue
        s = run_series(spec.name)
        idx = np.searchsorted(ref["t_cy"], s["t_cy"])
        idx = np.clip(idx, 0, len(ref["t_cy"]) - 1)
        ok = np.abs(ref["t_cy"][idx] - s["t_cy"]) * CENTURY_DAYS < 1e-6
        rr = op.fit_secular(ref["t_cy"][idx[ok]], ref["psi"][idx[ok]])["rate"]
        rate = op.fit_secular(s["t_cy"][ok], s["psi"][ok])["rate"]
        out["runs"][f"{spec.dt_days:g}"] = {"dt_days": spec.dt_days, "rate": rate,
                                            "ias15_same_times": rr, "error": rate - rr,
                                            "n_compared": int(ok.sum())}
    return out
