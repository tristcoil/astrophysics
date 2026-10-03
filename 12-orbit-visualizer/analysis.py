"""Every measurement in the project, made from the cached Horizons data.

The analysis never touches the network (``horizons.load`` reads the cache and
nothing else), and it never looks up the answer it is trying to measure:
published values from ``config.PUBLISHED`` are only brought in at the end of
each section, to grade the measurement.

Five questions, one per section:

1. **Kepler's third law** -- does P**2 / a**3 come out the same for every
   planet, and what does the Sun weigh? (And why is it *not quite* the same?)
2. **Kepler's second law** -- does Halley's comet sweep equal areas in equal
   times while its speed changes sixty-fold?
3. **The interstellar visitors** -- how fast were they moving before the Sun
   pulled on them, where in the sky did they come from, and how do they move
   through the Galaxy?
4. **The escaping spacecraft** -- what did each planetary flyby do to them,
   and where are they going?
5. **Voyager 1's light-day** -- exactly when is it 24 light-hours from Earth?

Plus two "right now" pieces: where everything is on the snapshot date, and
what is in the sky tonight from the observing site.
"""

from __future__ import annotations

import json
import re
import warnings

import numpy as np

import horizons as hz
import orbit_physics as op
from config import (
    COMETS, DWARFS, FLYBYS, G_SI, G_SI_REL_ERR, INTERSTELLAR, OBSERVER, PLANETS,
    PUBLISHED, RESULTS_PATH, SNAPSHOT_DATE, SPACECRAFT, SUN_PECULIAR_UVW_KMS, body,
)

warnings.filterwarnings("ignore", message=".*dubious year.*")

# Earth-orientation data: use the tables that ship with astropy rather than
# downloading fresh ones, so the analysis stays offline. The bundled IERS-A
# predictions reach about a year past the package date, which covers the
# snapshot night. Astropy normally refuses predictions older than 30 days; a
# few months' drift in predicted UT1 is milliseconds of Earth rotation, far
# below anything that shows in a planet's altitude, so that check is relaxed.
from astropy.utils import iers as _iers

_iers.conf.auto_download = False
_iers.conf.auto_max_age = None

#: Only use states beyond this distance when reading off an interstellar
#: object's speed "at infinity": out here the whole solar system pulls like a
#: single point mass at the barycentre, so the two-body formula is exact enough.
FAR_FIELD_AU = 100.0

#: Days either side of a flyby at which the heliocentric orbit is read before
#: and after -- far enough out that the planet's own gravity no longer matters.
FLYBY_LEVER_DAYS = 40.0


def tdb_to_utc_iso(jd_tdb: float, *, seconds: bool = True, precision: int = 0) -> str:
    """A TDB Julian date as a UTC ISO string (TDB runs ~69.2 s ahead of UTC now)."""
    from astropy.time import Time

    t = Time(jd_tdb, format="jd", scale="tdb").utc
    t.precision = precision
    s = t.isot
    return s if seconds else s[:16]


def _constellation(ra_deg: float, dec_deg: float) -> str:
    import astropy.units as u
    from astropy.coordinates import SkyCoord, get_constellation

    name = get_constellation(SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg))
    return {"Ophiucus": "Ophiuchus"}.get(name, name)  # astropy's table misspells it


def _sky_separation_deg(ra1, dec1, ra2, dec2) -> float:
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    a = SkyCoord(ra=ra1 * u.deg, dec=dec1 * u.deg)
    b = SkyCoord(ra=ra2 * u.deg, dec=dec2 * u.deg)
    return float(a.separation(b).deg)


def body_gms() -> dict[str, float]:
    """Each planet system's GM (au^3/day^2), read from Horizons' own element headers.

    Horizons fits each planet's heliocentric osculating orbit with
    ``GM_sun + GM_planet``; the difference from the Sun's GM is the planet's.
    """
    out = {}
    for b in PLANETS + DWARFS:
        gm = hz.load(hz.elements_key(b))["meta"]["keplerian_gm_au3_d2"]
        out[b.key] = max(gm - op.GM_SUN, 0.0)
    return out


def total_gm() -> float:
    """The Sun plus every planet: the mass a distant object feels."""
    return op.GM_SUN + sum(body_gms().values())


# ==========================================================================
# 1. Kepler's third law
# ==========================================================================

def kepler_third_law() -> tuple[dict, dict]:
    gms = body_gms()
    bodies = PLANETS + DWARFS
    rows, per_orbit_gm = {}, {}
    for b in bodies:
        d = hz.load(hz.orbit_key(b))
        m = op.measure_orbit(d["jd"], d["r"], d["v"])
        gm_eff = float(op.gm_from_orbit(m["a_au"], m["period_days"]))
        gm_k = op.gm_from_orbit(m["a_per_orbit_au"], m["periods_days"])
        per_orbit_gm[b.key] = gm_k
        # Scatter of single-orbit values, divided down by the number of orbits.
        # Perturbations make P and a wander *together* (a stretched orbit is
        # also a slower one), so this is far smaller than the P and a scatters
        # would suggest separately.
        err = float(np.std(gm_k, ddof=1) / np.sqrt(len(gm_k))) if len(gm_k) >= 3 else None
        rows[b.key] = {
            "name": b.name, "n_orbits": m["n_orbits"], "span_years": m["span_days"] / 365.25,
            "period_days": m["period_days"], "period_scatter_days": m["period_scatter_days"],
            "a_au": m["a_au"], "a_scatter_au": m["a_scatter_au"], "e": m["e"],
            "gm_eff": gm_eff, "gm_eff_err": err,
            "excess_ppm": (gm_eff / op.GM_SUN - 1) * 1e6,
            "excess_err_ppm": err / op.GM_SUN * 1e6 if err is not None else None,
            "own_mass_ppm": gms[b.key] / op.GM_SUN * 1e6,
        }

    # The Newtonian prediction: the body's own mass, plus every other planet
    # smeared into a ring along its orbit.
    keys = [b.key for b in bodies]
    a_all = [rows[k]["a_au"] for k in keys]
    gm_all = [gms[k] for k in keys]
    for i, k in enumerate(keys):
        pred, contrib = op.effective_central_gm(i, a_all, gm_all)
        rows[k]["predicted_ppm"] = (pred / op.GM_SUN - 1) * 1e6
        rows[k]["ring_terms_ppm"] = {keys[j]: c / op.GM_SUN * 1e6 for j, c in contrib.items()
                                     if abs(c / op.GM_SUN * 1e6) >= 0.05}
        rows[k]["residual_vs_model_ppm"] = rows[k]["excess_ppm"] - rows[k]["predicted_ppm"]
        std = PUBLISHED["planet_a_au"]["values"].get(k)
        rows[k]["a_vs_jpl_mean_ppm"] = (rows[k]["a_au"] / std - 1) * 1e6 if std else None

    # Kepler's law as a straight line: log P against log a, planets and Pluto.
    fit_keys = [b.key for b in PLANETS] + ["pluto"]
    la = np.log10([rows[k]["a_au"] for k in fit_keys])
    lp = np.log10([rows[k]["period_days"] for k in fit_keys])
    (slope, intercept), cov = np.polyfit(la, lp, 1, cov=True)

    # The Sun's GM from the four inner planets, where perturbations are
    # smallest -- first raw, then with each planet's predicted offset removed.
    inner = ["mercury", "venus", "earth", "mars"]
    raw = np.array([rows[k]["gm_eff"] for k in inner])
    corrected = np.array([rows[k]["gm_eff"] / (1 + rows[k]["predicted_ppm"] * 1e-6)
                          for k in inner])
    gm_sun = float(np.mean(corrected))
    gm_sun_err = float(np.std(corrected, ddof=1) / np.sqrt(len(corrected)))
    gm_sun_si = gm_sun * op.AU_KM**3 / op.DAY_S**2 * 1e9  # m^3/s^2
    mass_sun = gm_sun_si / G_SI

    # Jupiter's mass: its orbit's excess, minus what the other planets explain.
    j = rows["jupiter"]
    others = sum(v for k, v in j["ring_terms_ppm"].items())
    jup_ppm = (j["gm_eff"] / gm_sun - 1) * 1e6 - others
    jup_err_ppm = j["excess_err_ppm"]

    halley = measure_halley_orbit()

    result = {
        "bodies": rows,
        "loglog_slope": float(slope), "loglog_slope_err": float(np.sqrt(cov[0, 0])),
        "gm_sun_raw_inner_spread_ppm": float((raw.max() - raw.min()) / op.GM_SUN * 1e6),
        "gm_sun_measured_au3_d2": gm_sun,
        "gm_sun_measured_err_ppm": gm_sun_err / gm_sun * 1e6,
        "gm_sun_measured_si": gm_sun_si,
        "gm_sun_vs_de440_ppm": (gm_sun / op.GM_SUN - 1) * 1e6,
        "sun_mass_kg": mass_sun,
        "sun_mass_rel_err": G_SI_REL_ERR,
        "jupiter_mass_ratio_measured": jup_ppm * 1e-6,
        "jupiter_mass_ratio_err": jup_err_ppm * 1e-6 if jup_err_ppm else None,
        "jupiter_sun_over_jupiter_measured": 1e6 / jup_ppm,
        "jupiter_sun_over_jupiter_jpl": op.GM_SUN / gms["jupiter"],
        "halley": halley,
    }
    art = {"per_orbit_gm": per_orbit_gm}
    return result, art


def measure_halley_orbit() -> dict:
    """Halley's one complete orbit in the data: perihelion 1986 to perihelion 2061.

    Period from perihelion to perihelion. Size from the orbit's long axis:
    ``a = (perihelion distance + aphelion distance) / 2``, with the two
    perihelion distances averaged. (The conic fit used for the planets is the
    wrong tool at e = 0.967: there ``a = p / (1 - e**2)``, so an error of one
    part in 10^4 in e moves a by 0.3%.) Unlike the planets' averages over many
    orbits this is a single orbit, and Halley's period changes by a year or
    more between apparitions as Jupiter and Saturn tug on it -- so it is a test
    point for the law, not a precision measurement.
    """
    d = hz.load(hz.track_key(body("halley")))
    jd, r, v = d["jd"], d["r"], d["v"]
    rn, rdot = op.radial_distance_and_rate(r, v)
    peri = op.extrema(jd, rn, rdot, "min")
    (t1, q1), (t2, q2) = peri[0], peri[-1]
    Q = max(op.extrema(jd, rn, rdot, "max"), key=lambda x: x[1])[1]
    a = 0.5 * (0.5 * (q1 + q2) + Q)
    period = t2 - t1
    gm = float(op.gm_from_orbit(a, period))
    return {"period_days": float(period), "a_au": float(a), "e": float((Q - q1) / (Q + q1)),
            "excess_ppm": (gm / op.GM_SUN - 1) * 1e6}


# ==========================================================================
# 2. Kepler's second law: Halley's comet
# ==========================================================================

def _swept(jd, r, v, t0, t1):
    """Area (au^2) and in-plane angle (deg) swept between two times."""
    m = (jd >= t0) & (jd <= t1)
    area = float(np.trapezoid(op.areal_velocity(r[m], v[m]), jd[m]))
    normal = op.unit(np.cross(r[m], v[m]).mean(axis=0))
    theta, _ = op.in_plane_angle(r[m], v[m], normal)
    return area, float(np.degrees(theta[-1] - theta[0]))


def kepler_second_law() -> tuple[dict, dict]:
    from astropy.timeseries import LombScargle

    d = hz.load(hz.track_key(body("halley")))
    jd, r, v = d["jd"], d["r"], d["v"]
    rn, rdot = op.radial_distance_and_rate(r, v)
    peri = op.extrema(jd, rn, rdot, "min")
    aph = op.extrema(jd, rn, rdot, "max")
    t_p1, q1 = peri[0]
    t_p2, q2 = peri[-1]
    t_a, Q = aph[0]
    speed = op.to_kms(op.norm(v))
    av = op.areal_velocity(r, v)
    orbit = (jd >= t_p1) & (jd <= t_p2)
    mean_av = float(av[orbit].mean())
    dev = av[orbit] / mean_av - 1

    # The ripple: away from perihelion, where the slow comet is most sensitive
    # to the Sun being tugged around by Jupiter.
    far = orbit & (np.abs(jd - t_p1) > 730) & (np.abs(jd - t_p2) > 730)
    yrs = (jd[far] - jd[far][0]) / 365.25
    freq, power = LombScargle(yrs, av[far] / mean_av - 1).autopower(
        minimum_frequency=1 / 40, maximum_frequency=1 / 2, samples_per_peak=50)
    ripple_period = float(1 / freq[np.argmax(power)])

    # How big should it be? Jupiter swings the Sun around the barycentre at
    # (m_J / M_sun) times Jupiter's orbital speed. Halley's angular momentum is
    # its sideways speed times distance, so a 12 m/s wobble in the reference
    # point is a fractional change of (12 m/s) / (sideways speed).
    gms = body_gms()
    snap_j = hz.load(hz.snapshot_key(body("jupiter")))
    v_sun_reflex_kms = float(gms["jupiter"] / op.GM_SUN * op.to_kms(op.norm(snap_j["v"][0])))
    v_side = op.to_kms(op.norm(np.cross(r[far], v[far])) / rn[far])
    expected_ripple = v_sun_reflex_kms / float(v_side.min())

    # Equal areas in equal times: one year centred on perihelion vs one year
    # centred on aphelion.
    half = 365.25 / 2
    area_p, ang_p = _swept(jd, r, v, t_p1 - half, t_p1 + half)
    area_a, ang_a = _swept(jd, r, v, t_a - half, t_a + half)

    snap = hz.load(hz.snapshot_key(body("halley")))
    result = {
        "perihelion_1_utc": tdb_to_utc_iso(t_p1, seconds=False), "q1_au": q1,
        "perihelion_2_utc": tdb_to_utc_iso(t_p2, seconds=False), "q2_au": q2,
        "aphelion_utc": tdb_to_utc_iso(t_a, seconds=False), "Q_au": Q,
        "perihelion_interval_years": (t_p2 - t_p1) / 365.25,
        "speed_min_kms": float(speed[orbit].min()), "speed_max_kms": float(speed[orbit].max()),
        "speed_ratio": float(speed[orbit].max() / speed[orbit].min()),
        "areal_velocity_au2_per_day": mean_av,
        "areal_velocity_au2_per_year": mean_av * 365.25,
        "areal_velocity_max_dev_percent": float(np.abs(dev).max() * 100),
        "areal_velocity_rms_dev_percent": float(np.sqrt(np.mean(dev**2)) * 100),
        "ripple_period_years": ripple_period,
        "sun_reflex_speed_ms": v_sun_reflex_kms * 1000,
        "expected_ripple_percent": expected_ripple * 100,
        "year_at_perihelion": {"area_au2": area_p, "angle_deg": ang_p},
        "year_at_aphelion": {"area_au2": area_a, "angle_deg": ang_a},
        "equal_area_ratio": area_p / area_a,
        "today_r_au": float(op.norm(snap["r"][0])),
        "today_speed_kms": float(op.to_kms(op.norm(snap["v"][0]))),
        "today_inbound": bool(op.dot(snap["r"][0], snap["v"][0]) < 0),
        "source": d["meta"],
    }
    art = {"jd": jd, "r": r, "v": v, "speed_kms": speed, "areal": av, "mean_areal": mean_av,
           "t_peri": (t_p1, t_p2), "t_aph": t_a}
    return result, art


# ==========================================================================
# 3. The interstellar visitors
# ==========================================================================

def solar_apex() -> dict:
    """Where the Sun is heading relative to the local stars (the LSR).

    Stars -- or rocks -- sitting still in the LSR would stream towards us from
    exactly this direction, so it is the natural yardstick for a radiant.
    """
    uvw = np.array(SUN_PECULIAR_UVW_KMS)
    gal = uvw / np.linalg.norm(uvw)
    icrs = gal @ op._icrs_to_galactic_matrix()  # inverse of a rotation is its transpose
    ra, dec = op.spherical_deg(icrs)
    l, b = op.spherical_deg(gal)
    return {"ra_deg": float(ra), "dec_deg": float(dec), "l_deg": float(l), "b_deg": float(b),
            "speed_kms": float(np.linalg.norm(uvw)),
            "constellation": _constellation(float(ra), float(dec))}


def interstellar() -> tuple[dict, dict]:
    gm_tot = total_gm()
    apex = solar_apex()
    out, art = {}, {}
    for b in INTERSTELLAR:
        bary = hz.load(f"bary_{b.key}")
        passage = hz.load(f"passage_{b.key}")
        track = hz.load(hz.track_key(b))

        # Perihelion: closest approach, from the 2-hourly passage track.
        rn, rdot = op.radial_distance_and_rate(passage["r"], passage["v"])
        tp, q = min(op.extrema(passage["jd"], rn, rdot, "min"), key=lambda x: x[1])
        i = int(np.argmin(np.abs(passage["jd"] - tp)))
        el = op.elements_from_state(passage["r"][i], passage["v"][i])

        # Arrival: barycentric states far outside the planets. The earliest
        # one (1900, hundreds of au out) defines the incoming asymptote.
        rb, vb, jb = bary["r"], bary["v"], bary["jd"]
        dist = op.norm(rb)
        vinf = op.to_kms(op.hyperbolic_excess_speed(rb, vb, gm_tot))
        inbound = (jb < tp) & (dist > FAR_FIELD_AU)
        outbound = (jb > tp) & (dist > FAR_FIELD_AU)
        radiant, _, v_in = op.asymptotes(rb[0], vb[0], gm_tot)
        _, heading, v_out = op.asymptotes(rb[-1], vb[-1], gm_tot)
        ra, dec = op.radec_from_ecliptic(radiant)
        ra_out, dec_out = op.radec_from_ecliptic(heading)
        lb = op.spherical_deg(op.icrs_to_galactic(op.ecliptic_to_icrs(radiant)))

        # Galactic velocity: the incoming velocity at infinity points *away*
        # from the radiant.
        uvw = op.galactic_uvw_kms(-float(v_in) * radiant)
        uvw_lsr = uvw + np.array(SUN_PECULIAR_UVW_KMS)

        pub = PUBLISHED[b.key]
        comp = {
            "v_inf_kms": {"ours": float(op.to_kms(v_in)), "published": pub["v_inf_kms"]},
            "radiant_separation_deg": _sky_separation_deg(
                float(ra), float(dec), pub["radiant_ra_deg"], pub["radiant_dec_deg"]),
        }
        if "uvw_kms" in pub:
            comp["uvw_diff_kms"] = [float(x) for x in uvw - np.array(pub["uvw_kms"])]
        if "e" in pub:
            comp["e"] = {"ours": float(el.e), "published": pub["e"]}
            comp["q_au"] = {"ours": float(q), "published": pub["q_au"]}

        snap = hz.load(hz.snapshot_key(b))
        out[b.key] = {
            "name": b.name,
            "perihelion_utc": tdb_to_utc_iso(tp, seconds=False), "q_au": float(q),
            "e": float(el.e), "inclination_deg": float(np.degrees(el.inc)),
            "speed_at_perihelion_kms": float(op.to_kms(op.norm(passage["v"][i]))),
            "escape_speed_at_perihelion_kms": float(op.to_kms(op.escape_speed(q))),
            "v_inf_in_kms": float(np.nanmean(vinf[inbound])),
            "v_inf_out_kms": float(np.nanmean(vinf[outbound])),
            "v_inf_change_ms": float((np.nanmean(vinf[outbound]) - np.nanmean(vinf[inbound])) * 1e3),
            "radiant_ra_deg": float(ra), "radiant_dec_deg": float(dec),
            "radiant_l_deg": float(lb[0]), "radiant_b_deg": float(lb[1]),
            "radiant_constellation": _constellation(float(ra), float(dec)),
            "radiant_from_apex_deg": _sky_separation_deg(float(ra), float(dec),
                                                         apex["ra_deg"], apex["dec_deg"]),
            "heading_ra_deg": float(ra_out), "heading_dec_deg": float(dec_out),
            "heading_constellation": _constellation(float(ra_out), float(dec_out)),
            "uvw_kms": [float(x) for x in uvw],
            "uvw_lsr_kms": [float(x) for x in uvw_lsr],
            "speed_vs_lsr_kms": float(np.linalg.norm(uvw_lsr)),
            "today_r_au": float(op.norm(snap["r"][0])),
            "today_speed_kms": float(op.to_kms(op.norm(snap["v"][0]))),
            "comparison": comp, "published_source": pub["source"],
            "orbit_solution": passage["meta"].get("target_source"),
            "n_obs": passage["meta"].get("n_obs"),
            "data_arc": passage["meta"].get("data_arc") or passage["meta"].get("obs_span"),
        }
        art[b.key] = {"track": track, "passage": passage, "bary_jd": jb, "bary_dist": dist,
                      "v_inf_t": vinf, "tp": tp}
    return {"objects": out, "solar_apex": apex}, art


# ==========================================================================
# 4. The escaping spacecraft
# ==========================================================================

def _planet_radius_km(text_meta_key: str) -> float | None:
    q = hz.plan_by_key()[text_meta_key]
    m = re.search(r"^Center radii\s*:\s*([0-9.]+)", q.path.read_text(), flags=re.MULTILINE)
    return float(m.group(1)) if m else None


def flybys() -> tuple[list, dict]:
    gms = body_gms()
    km3_s2 = op.AU_KM**3 / op.DAY_S**2
    rows, art = [], {}
    for craft, planet, _date in FLYBYS:
        key = hz.flyby_key(craft, planet)
        d = hz.load(key)
        jd, r, v = d["jd"], d["r"], d["v"]
        rn, rdot = op.radial_distance_and_rate(r, v)
        tca, rp = min(op.extrema(jd, rn, rdot, "min"), key=lambda x: x[1])
        mu = gms[planet]
        vinf = op.to_kms(op.hyperbolic_excess_speed(r, v, mu))
        radiant_in, _, _ = op.asymptotes(r[0], v[0], mu)
        _, heading_out, _ = op.asymptotes(r[-1], v[-1], mu)
        turn_obs = float(np.degrees(op.angle_between(-radiant_in, heading_out)))
        v_mean = 0.5 * (vinf[0] + vinf[-1])
        rp_km = rp * op.AU_KM
        turn_pred = float(np.degrees(op.flyby_turn_angle(v_mean, rp_km, mu * km3_s2)))

        # The heliocentric orbit before and after, from the cruise track.
        tr = hz.load(hz.track_key(body(craft)))
        i0 = np.searchsorted(tr["jd"], tca - FLYBY_LEVER_DAYS)
        i1 = np.searchsorted(tr["jd"], tca + FLYBY_LEVER_DAYS)
        before = op.elements_from_state(tr["r"][i0], tr["v"][i0])
        after = op.elements_from_state(tr["r"][i1], tr["v"][i1])
        radius = _planet_radius_km(key)
        rows.append({
            "craft": craft, "planet": planet,
            "closest_approach_utc": tdb_to_utc_iso(tca, seconds=False),
            "closest_approach_km": rp_km,
            "altitude_km": rp_km - radius if radius else None,
            "planet_radius_km": radius,
            "v_inf_in_kms": float(vinf[0]), "v_inf_out_kms": float(vinf[-1]),
            "turn_observed_deg": turn_obs, "turn_predicted_deg": turn_pred,
            "helio_speed_before_kms": float(op.to_kms(op.norm(tr["v"][i0]))),
            "helio_speed_after_kms": float(op.to_kms(op.norm(tr["v"][i1]))),
            "e_before": float(before.e), "e_after": float(after.e),
            "bound_before": bool(before.energy < 0), "bound_after": bool(after.energy < 0),
        })
        art[key] = {"jd": jd, "dist_km": rn * op.AU_KM, "vinf_kms": vinf, "tca": tca}
    return rows, art


def spacecraft() -> tuple[dict, dict]:
    gm_tot = total_gm()
    fly, fly_art = flybys()
    out, art = {}, {"flybys": fly_art}
    for b in SPACECRAFT:
        tr = hz.load(hz.track_key(b))
        snap = hz.load(hz.snapshot_key(b))
        geo = hz.load(hz.snapshot_geo_key(b))
        r, v = snap["r"][0], snap["v"][0]
        _, heading, vinf = op.asymptotes(tr["r"][-1], tr["v"][-1], gm_tot)
        ra, dec = op.radec_from_ecliptic(heading)
        gra, gdec = op.radec_from_ecliptic(geo["r"][0])
        energy = op.specific_energy(tr["r"], tr["v"])
        launched_unbound = bool(energy[0] > 0)
        escape_flyby = next((f for f in fly if f["craft"] == b.key and f["bound_before"]
                             and not f["bound_after"]), None)
        speed = float(op.to_kms(op.norm(v)))
        out[b.key] = {
            "name": b.name,
            "today_r_au": float(op.norm(r)),
            "today_speed_kms": speed,
            "today_au_per_year": speed * 365.25 * op.DAY_S / op.AU_KM,
            "today_ecliptic_lat_deg": float(op.spherical_deg(r)[1]),
            "today_earth_distance_au": float(geo["rg"][0]),
            "today_light_hours": float(geo["rg"][0] / op.C_AU_DAY * 24),
            "today_constellation": _constellation(float(gra), float(gdec)),
            "heading_ra_deg": float(ra), "heading_dec_deg": float(dec),
            "heading_constellation": _constellation(float(ra), float(dec)),
            "v_inf_kms": float(op.to_kms(vinf)),
            "launched_unbound": launched_unbound,
            "escape_by": escape_flyby["planet"] if escape_flyby else None,
            "trajectory_source": tr["meta"].get("target_source"),
            "track_end": tr["date"][-1],
        }
        art[b.key] = {"jd": tr["jd"], "r": tr["r"], "v": tr["v"]}

    # Where the Voyagers were when their own instruments saw them leave the
    # heliosphere: the trajectory against the published distances.
    helio = {}
    for key, h in PUBLISHED["heliopause"].items():
        from astropy.time import Time

        tr = art[key]
        t = Time(h["date"], scale="utc").tdb.jd
        pos = op.Track(tr["jd"], tr["r"], tr["v"]).position(t)[0]
        helio[key] = {"date": h["date"], "ours_au": float(op.norm(pos)),
                      "published_au": h["au"], "source": h["source"]}
    return {"craft": out, "flybys": fly, "heliopause": helio}, art


# ==========================================================================
# 5. Voyager 1 and the light-day
# ==========================================================================

def light_day() -> tuple[dict, dict]:
    from astropy.time import Time
    from scipy.optimize import brentq

    g = hz.load("lightday_v1_geo")
    jd, rg, rr = g["jd"], g["rg"], g["rr"]
    L = op.LIGHT_DAY_AU
    cross = op.crossings(jd, rg, rr, L)
    up = [t for t, s in cross if s > 0]
    t_cross = up[0]
    i = int(np.searchsorted(jd, t_cross))
    rate_kms = float(op.to_kms(np.interp(t_cross, jd, rr)))

    # After the crossing, does Earth's own orbit ever carry it back inside?
    later_min = [(t, x) for t, x in op.extrema(jd, rg, rr, "min") if t > t_cross]
    maxima = op.extrema(jd, rg, rr, "max")
    t_min, r_min = later_min[0]

    # Light-travel-time versions, solved from barycentric positions: a signal
    # *sent* from Earth takes 24 h to arrive earlier than the geometric moment
    # (Voyager keeps receding while it travels); a signal *received* at Earth
    # left Voyager 24 h earlier still.
    v1 = op.Track(*(hz.load("lightday_v1_ssb")[k] for k in ("jd", "r", "v")))
    ea = op.Track(*(hz.load("lightday_earth_ssb")[k] for k in ("jd", "r", "v")))

    def uplink(t):
        tau = 1.0
        for _ in range(8):
            tau = float(op.norm(v1.position(t + tau)[0] - ea.position(t)[0]) / op.C_AU_DAY)
        return tau

    def downlink(t):
        tau = 1.0
        for _ in range(8):
            tau = float(op.norm(ea.position(t)[0] - v1.position(t - tau)[0]) / op.C_AU_DAY)
        return tau

    lo, hi = t_cross - 10, t_cross + 10
    t_up = brentq(lambda t: uplink(t) - 1.0, lo, hi, xtol=1e-9)
    t_down = brentq(lambda t: downlink(t) - 1.0, lo, hi, xtol=1e-9)

    # And one light-day from the Sun, from the heliocentric track.
    tr = hz.load(hz.track_key(body("voyager1")))
    rn, rdot = op.radial_distance_and_rate(tr["r"], tr["v"])
    t_sun = [t for t, s in op.crossings(tr["jd"], rn, rdot, L) if s > 0][0]

    pub = PUBLISHED["voyager1_light_day"]
    t_pub = Time(pub["utc"], scale="utc").tdb.jd
    result = {
        "light_day_au": L,
        "crossing_utc": tdb_to_utc_iso(t_cross, precision=1),
        "published_utc": pub["utc"],
        "difference_s": float((t_cross - t_pub) * op.DAY_S),
        "range_rate_kms": rate_kms,
        "definition": "geometric distance from Earth's centre = c x 1 day",
        "uplink_utc": tdb_to_utc_iso(t_up, seconds=False),
        "downlink_utc": tdb_to_utc_iso(t_down, seconds=False),
        "next_minimum_utc": tdb_to_utc_iso(t_min, seconds=False),
        "next_minimum_au": r_min,
        "margin_au": r_min - L,
        "margin_light_minutes": (r_min - L) / op.C_AU_DAY * 1440,
        "annual_swing_au": float(maxima[0][1] - [x for t, x in op.extrema(jd, rg, rr, "min")
                                                  if t > maxima[0][0]][0]),
        "sun_light_day_utc": tdb_to_utc_iso(t_sun, seconds=False),
        "published_source": pub["source"],
    }
    art = {"jd": jd, "rg": rg, "t_cross": t_cross, "t_min": t_min, "r_min": r_min,
           "maxima": maxima, "minima": op.extrema(jd, rg, rr, "min")}
    return result, art


# ==========================================================================
# Right now: the snapshot, and the sky tonight
# ==========================================================================

def snapshot() -> dict:
    out = {}
    for b in PLANETS + DWARFS + COMETS + INTERSTELLAR + SPACECRAFT:
        s = hz.load(hz.snapshot_key(b))
        r, v = s["r"][0], s["v"][0]
        lon, lat = op.spherical_deg(r)
        el = op.elements_from_state(r, v)
        out[b.key] = {"name": b.name, "group": b.group, "r": [float(x) for x in r],
                      "dist_au": float(op.norm(r)), "ecl_lon_deg": float(lon),
                      "ecl_lat_deg": float(lat), "speed_kms": float(op.to_kms(op.norm(v))),
                      "escape_speed_kms": float(op.to_kms(op.escape_speed(op.norm(r)))),
                      "bound": bool(el.energy < 0), "source": s["meta"].get("target_source")}
    return {"date_tdb": SNAPSHOT_DATE, "bodies": out}


def elements_crosscheck() -> dict:
    """Our state-vector -> elements conversion against Horizons' own elements.

    Same epoch, same state, same GM (Horizons' 'Keplerian GM'): the two must
    agree to rounding, which tests the conversion on real orbits of every
    inclination and eccentricity in the set.
    """
    worst = {"a_rel": 0.0, "e_abs": 0.0, "inc_deg": 0.0}
    for b in PLANETS + DWARFS:
        el_h = hz.load(hz.elements_key(b))
        s = hz.load(hz.snapshot_key(b))
        mine = op.elements_from_state(s["r"][0], s["v"][0], el_h["meta"]["keplerian_gm_au3_d2"])
        worst["a_rel"] = max(worst["a_rel"], abs(float(mine.a) / el_h["a"][0] - 1))
        worst["e_abs"] = max(worst["e_abs"], abs(float(mine.e) - el_h["e"][0]))
        worst["inc_deg"] = max(worst["inc_deg"], abs(float(np.degrees(mine.inc)) - el_h["inc"][0]))
    return worst


def sky_tonight() -> tuple[dict, dict]:
    sun = hz.load(hz.sky_key("sun"))
    t, sun_el = sun["jd_utc"], sun["el"]
    off = OBSERVER["utc_offset_hours"] / 24

    def local(jd_utc):
        from astropy.time import Time

        return Time(jd_utc + off, format="jd", scale="utc").iso[11:16]

    def when(level, rising):
        s = np.sign(sun_el - level)
        idx = np.nonzero((s[:-1] < s[1:]) if rising else (s[:-1] > s[1:]))[0]
        if not len(idx):
            return None
        i = idx[0]
        f = (level - sun_el[i]) / (sun_el[i + 1] - sun_el[i])
        return float(t[i] + f * (t[i + 1] - t[i]))

    dark = sun_el < -12  # nautical twilight or darker
    twilight = {name: {"evening": local(when(lvl, False)), "morning": local(when(lvl, True))}
                for name, lvl in (("sunset", -0.833), ("civil", -6), ("nautical", -12),
                                  ("astronomical", -18))}
    planets, art = {}, {"jd_utc": t, "sun_el": sun_el, "offset": off}
    for key, name, _hid in hz.SKY_TARGETS[1:]:
        d = hz.load(hz.sky_key(key))
        el_dark = np.where(dark, d["el"], -90)
        imax = int(np.argmax(el_dark))
        up = dark & (d["el"] > 10)
        planets[key] = {
            "name": name, "magnitude": float(np.nanmedian(d["mag"])),
            "constellation": _constellation(float(d["ra"][imax]), float(d["dec"][imax])),
            "max_alt_dark_deg": float(el_dark[imax]),
            "max_alt_time_local": local(t[imax]) if el_dark[imax] > 0 else None,
            "hours_above_10deg_dark": float(up.sum() * 10 / 60),
            "visible": bool(up.any()),
        }
        art[key] = d

    # Cross-check against astropy's own ephemeris at local midnight: two
    # independent codes must put Saturn in the same place in the sky.
    check = _astropy_sky_check(art["saturn"])
    return {"site": OBSERVER, "twilight_local": twilight, "bodies": planets,
            "astropy_check": check}, art


def _astropy_sky_check(saturn) -> dict:
    import astropy.units as u
    from astropy.coordinates import AltAz, EarthLocation, get_body
    from astropy.time import Time

    i = len(saturn["jd_utc"]) // 2
    when = Time(saturn["jd_utc"][i], format="jd", scale="utc")
    site = EarthLocation(lat=OBSERVER["lat"] * u.deg, lon=OBSERVER["lon"] * u.deg,
                         height=OBSERVER["elevation_km"] * u.km)
    frame = AltAz(obstime=when, location=site, pressure=1010 * u.hPa,
                  temperature=10 * u.deg_C, relative_humidity=0.5, obswl=0.55 * u.micron)
    aa = get_body("saturn", when, site).transform_to(frame)
    return {"time_utc": when.iso[:16],
            "alt_diff_arcsec": float((aa.alt.deg - saturn["el"][i]) * 3600),
            "az_diff_arcsec": float(((aa.az.deg - saturn["az"][i] + 180) % 360 - 180) * 3600
                                    * np.cos(np.radians(saturn["el"][i])))}


# ==========================================================================
# Everything
# ==========================================================================

def run_all() -> tuple[dict, dict]:
    results, art = {}, {}
    results["kepler3"], art["kepler3"] = kepler_third_law()
    results["kepler2"], art["kepler2"] = kepler_second_law()
    results["interstellar"], art["interstellar"] = interstellar()
    results["spacecraft"], art["spacecraft"] = spacecraft()
    results["light_day"], art["light_day"] = light_day()
    results["snapshot"] = snapshot()
    results["sky"], art["sky"] = sky_tonight()
    results["elements_crosscheck"] = elements_crosscheck()

    nh = results["spacecraft"]["craft"]["new_horizons"]
    at = results["interstellar"]["objects"]["atlas"]
    results["spacecraft"]["new_horizons_vs_atlas_radiant_deg"] = _sky_separation_deg(
        nh["heading_ra_deg"], nh["heading_dec_deg"], at["radiant_ra_deg"], at["radiant_dec_deg"])
    return results, art


def _jsonable(x):
    if isinstance(x, dict):
        return {k: _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist())
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, float) and not np.isfinite(x):
        return None
    return x


def save_results(results: dict) -> None:
    RESULTS_PATH.write_text(json.dumps(_jsonable(results), indent=1))


if __name__ == "__main__":
    res, _ = run_all()
    save_results(res)
    print(json.dumps(_jsonable(res), indent=1)[:6000])
