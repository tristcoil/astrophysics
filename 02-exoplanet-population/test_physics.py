"""Self-checks against cases with known answers.

Run with ``python test_physics.py``. No test framework and no network needed.

Four groups:

1. **Physics** -- every formula in ``planet_physics.py`` against a textbook case
   (Earth's year gives 1 AU, Jupiter moves the Sun at 12.5 m/s, and so on).
2. **Measurement algorithms** -- the valley finder and the tilt fit, fed inputs
   constructed so the right answer is known in advance. These inputs are
   test fixtures, built deterministically; no figure ever uses them.
3. **Catalogue hygiene** -- the flags that keep the archive's calculated values
   out of the science, on a hand-made miniature table.
4. **Cached data** -- sanity checks on the downloaded data and model curves,
   skipped when ``fetch_data.py`` has not been run.
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

import planet_physics as pp

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float, unit: str = "") -> None:
    """Assert ``got`` is within ``tol`` (fractional) of ``expected``."""
    got, expected = float(got), float(expected)
    rel = abs(got - expected) / abs(expected) if expected else abs(got)
    status = "ok  " if rel <= tol else "FAIL"
    if rel > tol:
        FAILURES.append(label)
    print(f"  [{status}] {label:<50} {got:>11.5g} vs {expected:>10.5g} {unit:<6} "
          f"({rel * 100:.2f}% off, tol {tol * 100:g}%)")


def check_true(label: str, condition: bool) -> None:
    if not condition:
        FAILURES.append(label)
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label}")


def physics() -> None:
    print("Kepler's third law")
    check("a(365.256 d, 1 Msun)", pp.semi_major_axis_au(365.256, 1.0), 1.0, 1e-4, "AU")
    check("a(11.862 yr, 1 Msun) = Jupiter", pp.semi_major_axis_au(11.862 * 365.25, 1.0), 5.204, 2e-3, "AU")
    check("P(1 AU, 1 Msun)", pp.period_days(1.0, 1.0), 365.256, 1e-4, "d")
    check("round trip a -> P -> a", pp.semi_major_axis_au(pp.period_days(0.3, 0.6), 0.6), 0.3, 1e-12, "AU")

    print("\nStarlight")
    check("L(1 Rsun, 5772 K)", pp.luminosity_lsun(1.0, 5772.0), 1.0, 1e-12, "Lsun")
    check("Earth insolation", pp.insolation(1.0, 1.0), 1.0, 1e-12)
    check("Venus insolation (a = 0.723 AU)", pp.insolation(1.0, 0.723), 1.913, 1e-3)

    print("\nHabitable zone (Kopparapu et al. 2014 coefficients)")
    for limit, coeffs in pp.KOPPARAPU_2014.items():
        check(f"S_eff at 5780 K = S_eff_sun ({limit})", pp.hz_flux_limit(limit, 5780.0),
              coeffs[0], 1e-12)
    # The Sun's habitable zone as usually quoted from these papers:
    # conservative 0.95-1.67 AU, optimistic 0.75-1.77 AU.
    check("Sun: runaway greenhouse (inner, conservative)",
          pp.hz_distance_au("runaway_greenhouse", 5772, 1.0), 0.95, 0.01, "AU")
    check("Sun: maximum greenhouse (outer, conservative)",
          pp.hz_distance_au("maximum_greenhouse", 5772, 1.0), 1.67, 0.01, "AU")
    check("Sun: recent Venus (inner, optimistic)",
          pp.hz_distance_au("recent_venus", 5772, 1.0), 0.75, 0.01, "AU")
    check("Sun: early Mars (outer, optimistic)",
          pp.hz_distance_au("early_mars", 5772, 1.0), 1.77, 0.01, "AU")
    check_true("every limit sits at lower flux for a 3000 K star than the Sun",
               all(pp.hz_flux_limit(k, 3000.0) < pp.hz_flux_limit(k, 5780.0) for k in pp.KOPPARAPU_2014))
    check_true("below 2600 K the fit is held at its 2600 K value",
               pp.hz_flux_limit("early_mars", 2400.0) == pp.hz_flux_limit("early_mars", 2600.0))
    zones = pp.hz_zone([1.0, 0.431, 1.913, 1.5, 0.33, 0.1], [5772.0] * 6)
    check_true(f"zones: Earth, Mars conservative; Venus out; S=1.5, 0.33 optimistic; "
               f"0.1 out -> {list(zones)}",
               list(zones) == ["conservative", "conservative", "", "optimistic", "optimistic", ""])

    print("\nRadial velocity")
    k_jup = pp.rv_semi_amplitude(pp.SOLAR_SYSTEM["Jupiter"]["mass"], 4332.59, 1.0, ecc=0.0489)
    check("Jupiter moves the Sun at", k_jup, 12.5, 0.02, "m/s")
    k_earth = pp.rv_semi_amplitude(1.0, 365.256, 1.0)
    check("Earth moves the Sun at", k_earth, 0.0895, 0.02, "m/s")
    m_1ms = pp.rv_minimum_mass_mearth(1.0, 1.0, 1.0)
    check("mass for K = 1 m/s at 1 AU (Sun-like)", m_1ms, 11.2, 0.01, "Me")
    check("round trip: that mass gives K", pp.rv_semi_amplitude(m_1ms, pp.period_days(1.0, 1.0), 1.0),
          1.0, 1e-4, "m/s")

    print("\nTransits and microlensing")
    check("Earth transit probability", pp.transit_probability(1.0, 1.0) * 100, 0.465, 2e-3, "%")
    check("hot Jupiter at 0.05 AU", pp.transit_probability(1.0, 0.05) * 100, 9.30, 2e-3, "%")
    # Independent route: theta_E = sqrt(kappa M pi_rel), kappa = 8.144 mas/Msun,
    # and r_E = theta_E * D_L. For 0.5 Msun at 4 kpc, source at 8 kpc: 2.854 AU.
    theta_mas = np.sqrt(8.144 * 0.5 * (1000 / 4000 - 1000 / 8000))
    check("Einstein radius vs theta_E * D_L", pp.einstein_radius_au(0.5, 4000, 8000),
          theta_mas / 1000 * 4000, 1e-3, "AU")

    print("\nBulk properties and reference values")
    check("Earth density (mean radius 6371.0 km)",
          pp.bulk_density_gcc(1.0, 6371.0e3 / pp.R_EARTH), 5.513, 2e-3, "g/cc")
    check("Jupiter in Earth masses", pp.SOLAR_SYSTEM["Jupiter"]["mass"], 317.8, 1e-3, "Me")
    check("Jupiter in Earth radii (equatorial)", pp.SOLAR_SYSTEM["Jupiter"]["radius"], 11.21, 1e-3, "Re")
    check("13 Jupiter masses in Earth masses", pp.DEUTERIUM_LIMIT_MEARTH, 4132, 1e-3, "Me")
    check("Van Eylen et al. valley at 10 days", pp.valley_radius(10.0, -0.09, 0.37), 1.905, 1e-3, "Re")


def _normal_sample(centre_log: float, width_dex: float, n: int) -> np.ndarray:
    """A deterministic 'sample': evenly spaced quantiles of a normal distribution."""
    return centre_log + width_dex * norm.ppf((np.arange(n) + 0.5) / n)


def algorithms() -> None:
    import analysis

    print("\nValley finder")
    lo, hi = np.log10(1.4), np.log10(2.5)
    two_peaks = np.concatenate([_normal_sample(lo, 0.05, 2000), _normal_sample(hi, 0.05, 2000)])
    v = analysis.find_valley(two_peaks)
    # Equal peaks of equal width: the minimum lies midway between them in log R.
    check("two equal peaks: valley at the log midpoint", v["valley"], 10 ** ((lo + hi) / 2), 0.01, "Re")
    check("two equal peaks: low peak", v["peak_low"], 1.4, 0.02, "Re")
    check("two equal peaks: high peak", v["peak_high"], 2.5, 0.02, "Re")
    check_true(f"two equal peaks: a real dip (contrast {v['contrast']:.2f} < 0.5)", v["contrast"] < 0.5)
    one_peak = analysis.find_valley(_normal_sample(np.log10(2.0), 0.1, 3000))
    check_true("one peak only: no interior minimum, flagged on_edge", one_peak["on_edge"])

    print("\nTilt fit")
    # Two populations either side of a valley whose radius falls with period as
    # P**-0.12. The fit should find that slope by making the valley deepest.
    true_slope, n = -0.12, 4000
    log_p = np.linspace(0.0, 2.0, n)  # 1-100 days
    side = np.where(np.arange(n) % 2 == 0, -0.13, 0.13)
    jitter = 0.04 * norm.ppf((np.arange(n) * 7919 % n + 0.5) / n)  # deterministic scramble
    log_r = np.log10(1.8) + true_slope * (log_p - 1.0) + side + jitter
    slopes = np.round(np.arange(-0.30, 0.2001, 0.01), 3)
    best, _ = analysis.fit_valley_tilt(log_r, log_p - 1.0, slopes)
    check_true(f"recovers slope {true_slope:+.2f}: got {best:+.2f}", abs(best - true_slope) <= 0.02)


def catalogue() -> None:
    import catalog

    print("\nCatalogue hygiene (hand-made table)")
    calc = "<a refstr=CALCULATED_VALUE href=/docs/pscp_calc.html>Calculated Value</a>"
    paper = "<a refstr=SMITH_2020 href=x>Smith et al. 2020</a>"
    df = pd.DataFrame({
        "pl_rade": [1.5, 2.0, 1.0, np.nan],
        "pl_rade_reflink": [paper, calc, paper, np.nan],
        "pl_radelim": [0, 0, 1, np.nan],
        "pl_bmasse": [5.0, 10.0, 3.0, 300.0],
        "pl_bmasse_reflink": [paper, paper, calc, paper],
        "pl_bmassprov": ["Mass", "Msini", "M-R relationship", "Mass"],
        "pl_bmasselim": [0, 0, 0, 1],
    })
    rad = df["pl_rade"].gt(0) & ~catalog._is_calculated(df["pl_rade_reflink"]) & catalog._no_limit(df["pl_radelim"])
    check_true(f"radius kept only when measured and not a limit -> {rad.tolist()}",
               rad.tolist() == [True, False, False, False])
    rel = df["pl_bmassprov"].eq(catalog.MASS_FROM_RELATION) | catalog._is_calculated(df["pl_bmasse_reflink"])
    mass = df["pl_bmasse"].gt(0) & ~rel & catalog._no_limit(df["pl_bmasselim"])
    check_true(f"mass from the M-R relation and mass upper limits dropped -> {mass.tolist()}",
               mass.tolist() == [True, True, False, False])
    err = catalog._relative_error(pd.Series([2.0]), pd.Series([0.1]), pd.Series([-0.3]))
    check("symmetrised relative error (+0.1/-0.3 on 2.0)", err.iloc[0], 0.1, 1e-12)
    check_true("reference names stripped of HTML",
               catalog._plain(paper) == "Smith et al. 2020")


def cached_data() -> None:
    from config import CATALOG_PATH, MODEL_DIR, ZENG_MODELS

    print("\nCached data (skipped if not fetched)")
    if not (CATALOG_PATH.exists() and all((MODEL_DIR / f).exists() for f, _ in ZENG_MODELS.values())):
        print("  [skip] run `python fetch_data.py` first")
        return
    import catalog
    from analysis import ModelCurve
    from fetch_data import load_model

    check("Zeng Earth-like model: radius at 1 Earth mass",
          ModelCurve(*load_model("earthlike")).radius(1.0), 1.0, 0.01, "Re")
    check("Zeng cold-hydrogen model at Jupiter's mass",
          ModelCurve(*load_model("hydrogen")).radius(317.8), 11.2, 0.05, "Re")

    df = catalog.load_catalog()
    samples = catalog.all_samples(df)
    mr = samples["mass_radius"].frame
    check_true("mass-radius sample holds no calculated radius or mass",
               not (mr["radius_calculated"].any() or mr["mass_from_relation"].any()))
    check_true("mass-radius sample holds no minimum masses or mass limits",
               bool(mr["mass_is_true"].all()) and bool(mr["pl_bmasselim"].fillna(0).eq(0).all()))
    check_true("no controversial planets outside the census",
               not any(s.frame["controversial"].any() for n, s in samples.items() if n != "census"))
    # Agol et al. (2021) irradiations for TRAPPIST-1 e, f, g.
    for name, s in (("TRAPPIST-1 e", 0.646), ("TRAPPIST-1 f", 0.373), ("TRAPPIST-1 g", 0.252)):
        row = df[df["pl_name"] == name]
        if len(row):
            check(f"{name} insolation (Agol et al. 2021)", row["pl_insol"].iloc[0], s, 0.02)


def main() -> int:
    physics()
    algorithms()
    catalogue()
    cached_data()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: {', '.join(FAILURES)}")
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
