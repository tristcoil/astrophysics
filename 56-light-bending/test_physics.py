"""Checks of the physics, the geometry and the transcriptions.

    python test_physics.py
"""

from __future__ import annotations

import sys
import traceback

import numpy as np

import config as C


def approx(a, b, rel=None, abs_=None):
    tol = max(rel * abs(b) if rel else 0.0, abs_ if abs_ else 0.0) or 1e-12 * abs(b)
    return abs(a - b) <= tol


# ---------------------------------------------------------------------------
# Physics
# ---------------------------------------------------------------------------

def test_limb_deflection():
    """4GM/(c^2 R) for the IAU nominal Sun is 1.7512"."""
    assert approx(C.DEFLECTION_LIMB_GR, 1.7512, abs_=5e-4), C.DEFLECTION_LIMB_GR


def test_deflection_formula_limits():
    """cot form -> 4GM/(c^2 b) near the Sun; 4.07 mas at 90 deg; 0 at 180 deg."""
    import stars as S

    d_au = 1.0
    chi = np.deg2rad(0.3)
    exact = S.deflection_arcsec(chi, d_au)
    small_angle = 4 * C.GM_SUN / (C.C_LIGHT**2 * d_au * C.AU * chi) / C.ARCSEC
    assert approx(exact, small_angle, rel=1e-5)
    assert approx(S.deflection_arcsec(np.pi / 2, 1.0) * 1e3, 4.07, abs_=0.01)
    assert S.deflection_arcsec(np.pi - 1e-9, 1.0) < 1e-8


def test_shapiro_log_approximation():
    """Near conjunction the exact one-way delay matches (1+g)GM/c^3 ln(4 r1 r2 / b^2)."""
    import shapiro

    r1 = np.array([C.AU, 0.0, 0.0])
    b = 2 * C.R_SUN
    r2 = np.array([-7.4 * C.AU, b * 8.4 / 1.0, 0.0])   # line passing ~b from the Sun
    bb = shapiro.impact_parameter(r1, r2)
    exact = shapiro.one_way_delay(r1, r2)
    approx_ = 2 * shapiro.GM_C3 * np.log(4 * np.linalg.norm(r1) * np.linalg.norm(r2) / bb**2)
    assert approx(exact, approx_, rel=1e-4), (exact, approx_)


def test_cassini_geometry():
    """Closest approach 1.6-1.7 solar radii on 2002 June 21; ~0.26 ms round trip."""
    import shapiro

    s = shapiro.summary()
    assert s["closest_date"].startswith("2002-06-21")
    assert 1.55 < s["b_min_rsun"] < 1.75
    assert 240 < s["delay_2way_max_us"] < 280
    assert s["y_vs_approx_max_rel"] < 0.01


def test_refraction_at_zenith_is_isotropic():
    """At the zenith refraction is a pure, isotropic change of scale of -A."""
    import refraction as R

    # a field passing through the zenith at Greenwich-ish latitude at some date
    from astropy.time import Time
    import astropy.units as u
    from astropy.coordinates import AltAz, EarthLocation, SkyCoord

    t = Time("2002-06-21T00:00:00")
    loc = EarthLocation.from_geodetic(0 * u.deg, 51.5 * u.deg, 0 * u.m)
    zen = SkyCoord(alt=89.999 * u.deg, az=0 * u.deg, frame=AltAz(obstime=t, location=loc)).icrs
    r = R.plate_effects(zen.ra.deg, zen.dec.deg, t.isot, 0.0, 51.5, 0.0, temp_c=10, pressure_hpa=1013.25,
                        half_width_deg=0.5)
    A = np.deg2rad(r["refco_A_arcsec"] / 3600)
    ref = r["refraction"]
    assert approx(ref["a"], -A, rel=0.02), (ref["a"], -A)
    assert approx(ref["e"], -A, rel=0.02)
    assert abs(ref["b"]) < 0.02 * A and abs(ref["d"]) < 0.02 * A


def test_lsq_recovers_parameters():
    import eclipse1919 as E

    rng = np.random.default_rng(1)
    A = np.column_stack([np.ones(50), rng.normal(size=50), rng.normal(size=50)])
    beta = np.array([0.3, -1.2, 2.0])
    y = A @ beta + rng.normal(scale=0.01, size=50)
    f = E.lsq(A, y)
    assert np.all(np.abs(f["beta"] - beta) < 5 * np.sqrt(np.diag(f["cov"])))


# ---------------------------------------------------------------------------
# The 1919 star field
# ---------------------------------------------------------------------------

def test_table1_frame_matches_gaia():
    """Table I's coordinates are Gaia's, in units of 50' of the 1919.4 equinox."""
    import stars as S

    f = S.fit_table1_frame()
    assert f["rms_arcsec"] < 2.5, f["rms_arcsec"]
    assert abs(f["scale"] - 1) < 3 * f["scale_err"] + 2e-4
    assert abs(f["rotation_deg"]) < 0.05


def test_table1_predictions_match_gr():
    """The 1919 predicted displacements (1.75" at the limb) agree with ours."""
    import stars as S

    f = S.fit_table1_frame()
    t1 = S.table1().set_index("no")
    for site in ("sobral", "principe"):
        pat = S.predicted_pattern(f, site, limb_value=1.75)
        d = np.hypot(pat["dx"] - t1[f"{site}_dx"], pat["dy"] - t1[f"{site}_dy"])
        assert d.max() < 0.08, (site, d.max())


def test_contact_times():
    """JPL's Moon and Sun give totality where the expeditions saw it: Sobral's
    second contact 11:58:18 UT (assumed, p.300); Principe 14:13:05-14:18:07 (computed, p.314).
    No duration is given for Sobral, but the exposures, timed from the "go" called as the
    crescent disappeared (p.299), ran until +305 s (astrographic No. 19): totality lasted longer."""
    import pandas as pd

    import analysis

    c = analysis.step_contacts()

    def sec(hms):
        h, m, s = hms.split(":")
        return 3600 * int(h) + 60 * int(m) + float(s)

    assert abs(sec(c["sobral"]["second_contact_utc"]) - sec("11:58:18")) < 20, c["sobral"]
    assert abs(sec(c["principe"]["second_contact_utc"]) - sec("14:13:05")) < 20, c["principe"]
    assert abs(sec(c["principe"]["third_contact_utc"]) - sec("14:18:07")) < 20, c["principe"]
    ex = pd.read_csv(C.TRANSCRIBED_DIR / "dyson1920_sobral_exposures.csv", comment="#", skipinitialspace=True)
    last_end = float((ex["t_after_second_contact_s_derived"] + ex["exposure_s"]).max())
    assert last_end == 305 and c["sobral"]["duration_s"] > last_end, (last_end, c["sobral"]["duration_s"])


def test_plate1_identified_and_dated():
    """The seven dash pairs of Plate 1 are the seven 4-inch stars (one assignment fits
    to ~1 px, the next is > 20 px off); the Moon's disc has JPL's radius to 1%, and
    its centre falls on JPL's track during totality."""
    import analysis

    p = analysis.step_plate1()
    assert p["rms_px"] < 2.0 and p["next_best_rms_px"] > 20.0, (p["rms_px"], p["next_best_rms_px"])
    assert not p["mirrored"]
    assert abs(p["moon"]["r_arcmin"] / p["moon"]["r_arcmin_horizons"] - 1) < 0.01
    assert p["miss_px"] < 3.0
    assert "11:58:10" <= p["dated_utc"] <= "12:03:28", p["dated_utc"]


def test_gamma_history_conversions():
    """Probable errors become sigma (x 1.4826); the 1919 Sobral 4-inch 1.98 +- 0.12 (p.e.)
    is (1+gamma)/2 = 1.131 +- 0.102; Cassini's gamma - 1 = (2.1 +- 2.3)e-5 is sigma_gamma 2.3e-5."""
    import history as H

    t = H.timeline().set_index("id")
    r = t.loc["E1919-sobral-4inch"]
    assert approx(r["half_1pg"], 1.98 / C.DEFLECTION_LIMB_GR, rel=1e-9)
    assert approx(r["sigma_half"], 0.12 / 0.6745 / C.DEFLECTION_LIMB_GR, rel=1e-9)
    assert approx(t.loc["T2002-cassini", "sigma_gamma"], 2.3e-5, rel=1e-9)
    fr = H.frontier()
    assert fr["sigma_gamma"].is_monotonic_decreasing


def test_4inch_replication():
    """Tables III-VII redone from Table II (with the two values the 1920 reduction used):
    Table VII's means and the 1.98" at the limb (p.306, p.331); every plate's alpha within
    the documented 1920 slips (dyson1920_sobral4inch_notes.md, D8: <= 0.0035 r)."""
    import eclipse1919 as E

    rep = E.replicate_4inch()
    m = rep["table7_means"]
    for (coord, kind), printed in {("RA", "eclipse"): 0.120, ("RA", "comparison"): 0.015,
                                   ("Dec", "eclipse"): 0.129, ("Dec", "comparison"): 0.031}.items():
        assert abs(m.loc[coord, kind] - printed) < 0.0012, (coord, kind, m.loc[coord, kind])
    assert abs(rep["limb_1920"] - 1.98) < 0.005, rep["limb_1920"]
    cmp = E.compare_printed_4inch(rep)
    sol = cmp[cmp["quantity"] == "sol_alpha"]
    assert sol["diff"].abs().max() < 0.0036, sol["diff"].abs().max()


def test_table2_misprints():
    """Plate I star 2 Dx and plate 14_2a star 2 Dx: the printed values do not give the printed
    normal equations; the 1920 working values (-0.798, -0.346) give them to the last digit."""
    import eclipse1919 as E

    for printed, tol in ((False, 0.0015), (True, None)):
        cmp = E.compare_printed_4inch(E.replicate_4inch(E.load_4inch(as_printed=printed)))
        for table, plate in (("III", "I"), ("IV", "14_2a")):
            d = cmp[(cmp["table"] == table) & (cmp["plate"] == plate) & cmp["quantity"].str.startswith("rhs")]
            if tol:
                assert d["diff"].abs().max() < tol, (table, plate, d["diff"].abs().max())
            else:
                assert d["diff"].abs().max() > 0.008, (table, plate)


def test_principe_replication():
    """The printed Delta4 y give the printed dk on every comparison (p.323-325), and the four
    limb values average to 1.65", with -0.04" to 1.61" (p.327-328)."""
    import principe as PR

    rep = PR.replicate()
    per = rep["per_comparison"]
    assert (per["dk"] - per["dk_printed"]).abs().max() < 0.8, per[["dk", "dk_printed"]]
    assert abs(rep["mean"] - 1.65) < 0.005 and abs(rep["limb"] - 1.61) < 0.005, rep["limb"]


def test_wallal_replication():
    """Table 9 of the 1928 bulletin from its Table 6: beta and weight per plate pair (inputs
    rounded to 0.01", so beta within 0.08"), and the mean 1.82" +- 0.15" (p.e.)."""
    import wallal as W

    rep = W.replicate()
    printed = {"M15_M7": (2.38, 0.85), "M16_M11": (1.66, 0.96), "M17_M9": (1.36, 0.83),
               "N18_N11": (2.00, 0.73), "N19_N13": (2.54, 0.89), "N20_N9": (1.01, 0.88)}
    per = rep["per_pair"].set_index("pair")
    for p, (b, w) in printed.items():
        assert abs(per.loc[p, "beta"] - b) < 0.08, (p, per.loc[p, "beta"])
        assert abs(per.loc[p, "weight"] - w) < 0.012, (p, per.loc[p, "weight"])
    assert abs(rep["beta"] - 1.82) < 0.01 and abs(rep["pe"] - 0.15) < 0.006, (rep["beta"], rep["pe"])


def test_astrographic_replication():
    """p.311: each plate's alpha + 0.240 e from the printed Dy and Table IX's d agrees with the
    printed value (three documented 1920 slips aside); the means give 0.93" and 1.52" (p.312)."""
    import pandas as pd

    import eclipse1919 as E

    rep = E.replicate_astro()
    pr = pd.read_csv(C.TRANSCRIBED_DIR / "dyson1920_p311_alpha_plus_e.csv", comment="#", skipinitialspace=True)
    pr["eclipse_plate"] = pr["eclipse_plate"].astype(str)
    m = rep["per_plate"].merge(pr, left_on="plate", right_on="eclipse_plate")
    diff = (m["alpha_plus_ke"] - m["value_rev"]).abs()
    # dyson1920_astrographic_notes.md: D2-D4 are slips in the 1920 arithmetic. Plates 7 and 12
    # appear twice in Table IX and the d used for them is not stated (D6); we use the mean d.
    slips = {"2", "6", "15"}
    ok = ~m["plate"].isin(slips | {"7", "12"})
    assert diff[ok].max() < 0.0006, diff[ok].max()
    assert diff[m["plate"].isin({"7", "12"})].max() < 0.001
    assert diff.max() < 0.004
    # the printed means used alpha rounded to 0.024 r and 0.039 r: 0.93" and 1.52"
    assert abs(rep["limb_plates_scale"] - 0.93) < 0.03 and abs(rep["limb_refraction_scale"] - 1.52) < 0.04


# ---------------------------------------------------------------------------

def main() -> int:
    tests = [(k, v) for k, v in globals().items() if k.startswith("test_") and callable(v)]
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"PASS  {name}")
        except Exception:
            failed += 1
            print(f"FAIL  {name}")
            traceback.print_exc(limit=2)
    print(f"{len(tests) - failed}/{len(tests)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
