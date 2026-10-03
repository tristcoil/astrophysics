"""Self-checks against cases with known answers.

Run with ``python test_physics.py``. No test framework and no network needed.

Four groups:

1. **Physics** -- every formula in ``neo_physics.py`` against a textbook case:
   H = 22 is 140 m, escape speed is 11.19 km/s, Chelyabinsk's light converts to
   441 kilotons, a numerically integrated flyby bends by exactly the angle the
   hyperbola formula says, and so on.
2. **Risk scales** -- the Torino scale on the historical cases (Apophis 2004
   rated 4, 2024 YR4 rated 3) and on every region of its chart; the Palermo
   scale on its own definition.
3. **Algorithms** -- the distance-index estimator and the near-miss-to-impact
   extrapolation, fed samples built so the right answer is known in advance.
   These are test fixtures, generated deterministically; no figure uses them.
4. **Cached data** -- JPL's numbers checked against each other and against our
   formulas; skipped when ``fetch_data.py`` has not been run.
"""

from __future__ import annotations

import sys
import warnings

import numpy as np
from scipy.integrate import solve_ivp

warnings.filterwarnings("ignore", message=".*dubious year.*")

import neo_physics as nphys
from config import AU_KM, DATA_DIR, EARTH_RADIUS_KM, GM_EARTH_KM3_S2

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float, unit: str = "") -> None:
    """Assert ``got`` is within ``tol`` (fractional) of ``expected``."""
    got, expected = float(got), float(expected)
    rel = abs(got - expected) / abs(expected) if expected else abs(got)
    status = "ok  " if rel <= tol else "FAIL"
    if rel > tol:
        FAILURES.append(label)
    print(f"  [{status}] {label:<56} {got:>12.6g} vs {expected:>11.6g} {unit:<5} "
          f"({rel * 100:.3f}% off, tol {tol * 100:g}%)")


def check_true(label: str, condition: bool) -> None:
    if not condition:
        FAILURES.append(label)
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label}")


# ---------------------------------------------------------------------------
# 1. Physics
# ---------------------------------------------------------------------------

def _integrate_flyby(b_km: float, v_inf: float) -> tuple[float, float]:
    """Two-body flyby, integrated numerically: closest approach and bend angle.

    Starts a billion km out -- far enough that the path there is straight to
    within 0.002 degrees -- with the speed energy conservation gives, pointed so
    the angular momentum is exactly b * v_inf (what "aimed to miss by b" means),
    and integrates Newton's law until it is as far away again. Closest approach
    is where the radial velocity changes sign. An independent check on the
    closed-form hyperbola formulas, which come from conservation laws."""
    r0 = 1.0e9
    big_r = np.hypot(r0, b_km)
    v0 = np.sqrt(v_inf**2 + 2 * GM_EARTH_KM3_S2 / big_r)
    theta = np.arcsin(b_km * v_inf / (v0 * big_r)) - np.arctan2(b_km, r0)
    state = [-r0, b_km, v0 * np.cos(theta), v0 * np.sin(theta)]

    def rhs(_t, s):
        x, y, vx, vy = s
        r3 = (x * x + y * y) ** 1.5
        return [vx, vy, -GM_EARTH_KM3_S2 * x / r3, -GM_EARTH_KM3_S2 * y / r3]

    def periapsis(_t, s):
        return s[0] * s[2] + s[1] * s[3]  # r . v, negative inbound

    periapsis.direction = 1
    sol = solve_ivp(rhs, (0, 2 * r0 / v_inf), state, method="DOP853", rtol=1e-12, atol=1e-8,
                    events=periapsis)
    q = float(np.hypot(*sol.y_events[0][0][:2]))
    v_in, v_out = sol.y[2:, 0], sol.y[2:, -1]
    bend = np.degrees(np.arccos(np.dot(v_in, v_out) / np.linalg.norm(v_in) / np.linalg.norm(v_out)))
    return q, bend


def physics() -> None:
    print("Size from brightness")
    check("H = 22, albedo 0.14 -> 141 m (the '140 m' convention)", nphys.diameter_km(22.0) * 1e3,
          141.40, 1e-4, "m")
    check("H = 17.75 -> 1 km", nphys.diameter_km(17.75), 1.0, 2e-3, "km")
    check("round trip D -> H -> D", nphys.diameter_km(nphys.absolute_magnitude(0.34)), 0.34, 1e-12, "km")
    check("albedo 0.05 makes it 1.67x bigger", nphys.diameter_km(22, 0.05) / nphys.diameter_km(22), 1.673, 1e-3)

    print("\nEncounter dynamics")
    check("escape speed at Earth's surface", nphys.escape_speed(), 11.186, 1e-4, "km/s")
    check("slowest possible impact = escape speed", nphys.impact_speed(0.0), 11.186, 1e-4, "km/s")
    check("capture radius at v_inf = v_esc is R*sqrt(2)",
          nphys.capture_radius_km(nphys.escape_speed()), EARTH_RADIUS_KM * np.sqrt(2), 1e-12, "km")
    check("a very fast asteroid is barely focused", nphys.capture_radius_km(1e4), EARTH_RADIUS_KM, 1e-6, "km")
    check("impact weight with Earth itself as reference", nphys.impact_weight(15.0, EARTH_RADIUS_KM), 1.0, 1e-12)
    # A 1 m sphere at 20 km/s: m = 2600 * pi/6 kg, E = m v^2 / 2
    check("1 m stony asteroid at 20 km/s", nphys.impact_energy_mt(1e-3, 20.0) * 1e3,
          0.5 * 2600 * np.pi / 6 * 4e8 / 4.184e12, 1e-12, "kt")
    check("diameter for energy inverts energy for diameter",
          nphys.diameter_for_energy_km(nphys.impact_energy_mt(0.05, 17.0), 17.0), 0.05, 1e-12, "km")

    for b_km, v_inf in ((50_000.0, 5.0), (20_000.0, 12.0)):
        q_num, bend_num = _integrate_flyby(b_km, v_inf)
        # invert b(q) numerically from the closed form
        qs = np.linspace(1_000, b_km, 200_001)
        q_formula = qs[np.argmin(np.abs(nphys.impact_parameter_km(qs, v_inf) - b_km))]
        check(f"flyby b = {b_km:,.0f} km, v = {v_inf:g} km/s: closest approach",
              q_num, q_formula, 2e-4, "km")
        check("  ... bend angle vs 2 arcsin(1/e)", bend_num, nphys.turning_angle_deg(q_num, v_inf), 2e-3, "deg")

    print("\nFireballs (Brown et al. 2002)")
    check("Chelyabinsk: 3.75e14 J of light -> impact energy", nphys.fireball_impact_energy_kt(3.75e14),
          441.0, 3e-3, "kt")
    check("impacts per year of at least 1 kt", nphys.brown2002_rate(1.0), 3.7, 1e-12, "/yr")
    check("Chelyabinsk-sized events: once per ~65 years", 1 / nphys.brown2002_rate(441.0), 64.83, 1e-4, "yr")

    print("\nOrbits")
    x, y, z = nphys.orbit_xyz(1.0, 0.0, 10.0, 30.0, 40.0, np.linspace(0, 360, 13))
    check("circular orbit stays at a", np.max(np.abs(np.sqrt(x**2 + y**2 + z**2) - 1.0)) + 1.0, 1.0, 1e-12, "au")
    x, y, z = nphys.orbit_xyz(2.0, 0.5, 0.0, 0.0, 0.0, [0.0, 180.0])
    check("e = 0.5: perihelion at a(1 - e)", np.hypot(x[0], y[0]), 1.0, 1e-12, "au")
    check("e = 0.5: aphelion at a(1 + e)", np.hypot(x[1], y[1]), 3.0, 1e-12, "au")

    # Yarkovsky: closed-form orbit average vs brute-force average of Gauss's equation
    a2, a_au, e = -2.902e-14, 0.9224, 0.1911
    mean_anomaly = np.linspace(0, 2 * np.pi, 200_001)[:-1]
    ecc = mean_anomaly.copy()
    for _ in range(50):  # Kepler's equation, by iteration
        ecc = mean_anomaly + e * np.sin(ecc)
    r = a_au * (1 - e * np.cos(ecc))
    p = a_au * (1 - e**2)
    h = np.sqrt(nphys.GM_SUN_AU3_D2 * p)
    dadt = (2 * a_au**2 / h) * (p / r) * a2 / r**2  # au/day, transverse push only
    brute = dadt.mean() * 365.25 * AU_KM * 1e3
    check("Yarkovsky drift: closed form vs orbit average", nphys.yarkovsky_drift_m_per_year(a2, a_au, e),
          brute, 1e-6, "m/yr")


# ---------------------------------------------------------------------------
# 2. Risk scales
# ---------------------------------------------------------------------------

def risk_scales() -> None:
    print("\nTorino scale (Binzel 2000, revised 2005)")
    cases = [
        ("Chelyabinsk (0.44 Mt), certain", 1.0, 0.44, 0),
        ("Tunguska-like 10 Mt, certain", 1.0, 10.0, 8),
        ("regional 1,000 Mt, certain", 1.0, 1e3, 9),
        ("Chicxulub-like 1e8 Mt, certain", 1.0, 1e8, 10),
        ("Apophis, Dec 2004: 2.7%, ~500 Mt (rated 4)", 0.027, 500.0, 4),
        ("2024 YR4, Feb 2025: 3.1%, ~8 Mt (rated 3)", 0.031, 8.0, 3),
        ("10% chance of a 7 m rock (2010 RF12-like)", 0.10, 0.0086, 0),
        ("1 in 10,000 of 1,000 Mt", 1e-4, 1e3, 1),
        ("1 in 2,000 of 10,000 Mt", 5e-4, 1e4, 2),
        ("50% of 1,000 Mt", 0.5, 1e3, 5),
        ("1 in 1,000 of 10 million Mt", 1e-3, 1e7, 6),
        ("10% of 1 million Mt", 0.1, 1e6, 7),
        ("1 in 100,000 of 100 Mt", 1e-5, 100.0, 0),
    ]
    for label, p, e, expected in cases:
        got = nphys.torino_scale(p, e)
        check_true(f"{label}: level {got} (expected {expected})", got == expected)
    lp, le = np.meshgrid(np.linspace(-9, 0, 400), np.linspace(-3, 9, 400))
    grid = nphys.torino_scale(10.0**lp, 10.0**le)
    check_true("every level 0-10 appears on the chart", set(np.unique(grid)) == set(range(11)))
    check_true("more probable never means a lower level", bool((np.diff(grid, axis=1) >= 0).all()))
    check_true("more energetic never means a lower level", bool((np.diff(grid, axis=0) >= 0).all()))

    print("\nPalermo scale (Chesley et al. 2002)")
    check("background frequency at 1 Mt", nphys.background_impact_rate(1.0), 0.03, 1e-12, "/yr")
    e, years = 50.0, 40.0
    p0 = nphys.background_impact_rate(e) * years
    check("PS = 0 when as likely as the background", nphys.palermo_scale(p0, e, years) + 1.0, 1.0, 1e-12)
    check("PS = -2 at 1% of the background", nphys.palermo_scale(p0 / 100, e, years), -2.0, 1e-12)
    # 1950 DA on 2026-09-29: Sentry lists p = 3.833e-4, E = 75,180 Mt, impact 2880-03-17, PS -0.92
    check("1950 DA (Sentry's own inputs) -> PS", nphys.palermo_scale(3.833e-4, 75_180.0, 2880.2077 - 2026.7448),
          -0.92, 6e-3)


# ---------------------------------------------------------------------------
# 3. Algorithms
# ---------------------------------------------------------------------------

def algorithms() -> None:
    import pandas as pd

    import analysis

    print("\nDistance index on samples with a known answer")
    rng = np.random.default_rng(1)
    d_max = 0.05
    # Uniform over a target disc: N(<d) grows as d^2 exactly
    disc = d_max * np.sqrt(rng.random(40_000))
    check("uniform disc (k = 2)", analysis.distance_index(disc, d_max), 2.0, 0.02)
    line = d_max * rng.random(40_000)  # uniform in d: N(<d) grows as d
    check("uniform in distance (k = 1)", analysis.distance_index(line, d_max), 1.0, 0.02)

    print("\nNear misses -> hits, on a synthetic stream with a known impact count")
    # One million asteroids aimed evenly over a disc of radius B in the plane of
    # aim, all at 8 km/s. Each one's true closest approach q follows from its
    # aim b by inverting b^2 = q^2 (1 + 2GM/(q v^2)); the impacts are those with
    # q below Earth's radius. The estimator sees only the approaches within a
    # reference distance and must predict the number of impacts.
    v = 8.0
    n = 1_000_000
    big_b = 60_000.0
    b = big_b * np.sqrt(rng.random(n))
    mu = GM_EARTH_KM3_S2 / v**2
    q = -mu + np.sqrt(mu**2 + b**2)  # the positive root of q^2 + 2 mu q - b^2 = 0
    impacts = int((q < EARTH_RADIUS_KM).sum())
    ref = 40_000.0  # every point with q < ref has b < big_b, so the sample is complete
    near = q < ref
    frame = pd.DataFrame({"des": np.arange(near.sum()).astype(str), "v_inf": v,
                          "dist_km": q[near]})
    predicted = analysis.extrapolation_weights(frame, ref).sum()
    check(f"predicted vs actual impacts ({impacts:,} actual)", predicted, impacts,
          3 / np.sqrt(impacts), "hits")

    print("\nBootstrap by asteroid")
    frame = pd.DataFrame({"des": ["a", "a", "b", "c"], "w": [1.0, 1.0, 2.0, 4.0]})
    (boot,) = analysis.bootstrap_object_sums(frame["des"], [frame["w"]], np.random.default_rng(0), 4000)
    # Each asteroid's total is 2, 2, 4: the resampled sum of three draws has mean 3 * 8/3 = 8
    check("bootstrap mean of summed weights", boot.mean(), 8.0, 0.01)

    print("\nParsing")
    import catalog
    days = catalog.parse_time_sigma_days(["< 00:01", "00:05", "3_19:55"])
    check("'< 00:01' -> 0", days[0] + 1, 1.0, 1e-12)
    check("'00:05' -> 5 minutes", days[1] * 1440, 5.0, 1e-9, "min")
    check("'3_19:55' -> 3 d 19 h 55 min", days[2], 3 + 19 / 24 + 55 / 1440, 1e-12, "d")
    check("decimal year of 2026-01-01", catalog.decimal_year("2026-01-01"), 2026.0, 1e-12)
    check("Sentry date '2880-03-17.00' (beyond pandas' range)",
          catalog.decimal_year("2880-03-17.00"), 2880 + (31 + 29 + 16) / 366, 1e-9)
    labels = list(catalog.size_class([1500.0, 999.0, 140.0, 139.0, 20.0, 5.0]))
    check_true(f"size classes {labels}",
               labels == ["≥ 1 km", "140 m – 1 km", "140 m – 1 km", "50 – 140 m", "20 – 50 m", "< 10 m"])


# ---------------------------------------------------------------------------
# 4. Cached data
# ---------------------------------------------------------------------------

def cached_data() -> None:
    print("\nCached data")
    if not (DATA_DIR / "manifest.json").exists():
        print("  (skipped: run `python fetch_data.py` first)")
        return
    import catalog
    from analysis import _closest_point

    cad = catalog.load_apophis_approaches()
    row = cad.loc[cad["dist_au"].idxmin()]
    check("Apophis 2029 miss distance (JPL close-approach table)", row["dist_km"], 38_011.5, 1e-4, "km")
    check("its speed at closest approach from v_inf and energy", nphys.speed_at_distance(row["dist_km"], row["v_inf"]),
          row["v_rel"], 1e-5, "km/s")

    fine = catalog.load_horizons_vectors("apophis_vectors_1")
    ca = _closest_point(fine)
    check("Horizons path: closest distance agrees with the table", ca["r_km"], row["dist_km"], 1e-5, "km")
    i = ca["index"]
    check("Horizons path: speed at closest approach", np.linalg.norm(fine[["vx", "vy", "vz"]].iloc[i]),
          row["v_rel"], 2e-4, "km/s")
    sigma = catalog.load_apophis_range_sigma()
    k = sigma["range_km"].idxmin()
    check("3σ: Horizons range uncertainty vs half the table's 3σ span",
          sigma["range_3sigma_km"].iloc[k], nphys.au_to_km(row["dist_max_au"] - row["dist_min_au"]) / 2,
          0.02, "km")

    rec = catalog.load_apophis_record()
    pars = {p["name"]: float(p["value"]) for p in rec["orbit"]["model_pars"] if p["kind"] == "EST"}
    els = {e["name"]: float(e["value"]) for e in rec["orbit"]["elements"]}
    drift = nphys.yarkovsky_drift_m_per_year(pars["A2"], els["a"], els["e"])
    check_true(f"Apophis Yarkovsky drift from JPL's A2: {drift:.0f} m/yr (inward, a few hundred m/yr)",
               -300 < drift < -100)

    fb = catalog.load_fireballs()
    ours = nphys.fireball_impact_energy_kt(fb["radiated_j"])
    rel = np.abs(ours / fb["impact_kt"] - 1)
    check_true(f"JPL's fireball impact energies are Brown et al.'s conversion "
               f"(median diff {np.median(rel):.2%}; differences are rounding)", np.median(rel) < 0.01)
    chel = fb.loc[fb["impact_kt"].idxmax()]
    check_true(f"largest fireball is Chelyabinsk: {chel['time']:%Y-%m-%d}, {chel['impact_kt']:.0f} kt",
               chel["time"].strftime("%Y-%m-%d") == "2013-02-15")
    tc3 = fb[fb["time"].dt.strftime("%Y-%m-%d") == "2008-10-07"]
    check_true("2008 TC3, the first asteroid seen before impact, is in the fireball record",
               len(tc3) == 1 and abs(tc3["lon"].iloc[0] - 31.4) < 1)

    vis = catalog.load_sentry_impacts()
    now = catalog.decimal_year(__import__("fetch_data").load_manifest()["sentry_impacts"]["retrieved_utc"][:19])
    ahead = vis["impact_year"] - now
    within = ahead <= 100
    ts = nphys.torino_scale(vis["ip"].to_numpy(), vis["energy"].to_numpy())
    check_true(f"Torino: ours equals Sentry's for all {within.sum():,} dates within 100 years",
               bool((ts[within] == vis.loc[within, "ts"]).all()))
    check_true("Sentry leaves exactly the dates beyond 100 years unrated",
               bool(((~within) == vis["ts"].isna()).all()))
    ps_ours = nphys.palermo_scale(vis["ip"], vis["energy"], np.maximum(ahead, 1.0))
    far = ahead > 2
    agree = float(((ps_ours - vis["ps"]).abs()[far] <= 0.02).mean())
    check_true(f"Palermo: {agree:.1%} of dates over 2 years ahead agree with Sentry to 0.02", agree > 0.99)

    ap = catalog.load_approaches()
    check_true("every approach found its asteroid in the catalogue", bool(ap["first_obs"].notna().all()))
    check_true("no approach listed twice", not ap.duplicated(["des", "jd"]).any())
    check_true("all approaches within 0.05 au", bool((ap["dist_au"] <= 0.05).all()))

    import json
    res = json.loads((DATA_DIR / "results.json").read_text()) if (DATA_DIR / "results.json").exists() else None
    if res:
        for label in ("≥ 1 km", "140 m – 1 km"):
            lo, hi = res["impact_rates"]["distance_index"][label]["k_2p5_97p5"]
            check_true(f"{label}: k = 2 inside the 95% range {lo:.2f}–{hi:.2f}", lo <= 2 <= hi)


def main() -> int:
    physics()
    risk_scales()
    algorithms()
    cached_data()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: {', '.join(FAILURES)}")
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
