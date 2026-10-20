"""Self-checks against cases with known answers.

Run with ``python test_physics.py``. No test framework and no network needed.

Groups:

1. **Geometry** -- phase angles, the Lambert sphere, image orientation (north
   up, east right), sidereal time at J2000, sub-points, angular sizes.
2. **Photometry** -- a synthetic disc image measured back: centre, radius,
   histograms and stored means; the display transform's inverse.
3. **Periodicity** -- an injected 24-h signal recovered by the periodogram,
   the autocorrelation and the thinning experiment; a frame-level brightness
   jump flagged. These simulated inputs are test fixtures; no figure uses them.
4. **Surface and Moon** -- the Lambert weighting of the land share; gallery
   file-name times; the event finder on a hand-made Moon track.
5. **Cached data** -- sanity checks on the downloads, skipped if absent.
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

import geometry as G
import periodicity as P
import photometry as ph
import transform as T
from config import R_EARTH_EQ_KM, R_MOON_KM

FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    status = "ok  " if condition else "FAIL"
    print(f"  [{status}] {name}" + (f"  ({detail})" if detail else ""))
    if not condition:
        FAILURES.append(name)


# ---------------------------------------------------------------------------
# 1. Geometry
# ---------------------------------------------------------------------------

def test_phase_and_lambert():
    sun = np.array([[1.5e8, 0, 0]])
    dsc = np.array([[1.5e6 * np.cos(np.radians(8)), 1.5e6 * np.sin(np.radians(8)), 0]])
    a = G.phase_angle_deg(dsc, sun)[0]
    check("phase angle of a spacecraft 8 deg off the Sun line", abs(a - 8.0) < 1e-9, f"{a:.6f}")
    check("Lambert sphere: Phi(0) = 1, Phi(90) = 1/pi",
          abs(G.lambert_phase(np.array([0.0]))[0] - 1) < 1e-12
          and abs(G.lambert_phase(np.array([90.0]))[0] - 1 / np.pi) < 1e-12)
    phi = G.lambert_phase(np.linspace(0, 180, 181))
    check("Lambert phase function falls monotonically", np.all(np.diff(phi) <= 1e-12))
    check("sunlit fraction of the disc: 1 at full phase, 1/2 at quadrature",
          abs(G.sunlit_fraction_of_disc(np.array([0.0]))[0] - 1) < 1e-12
          and abs(G.sunlit_fraction_of_disc(np.array([90.0]))[0] - 0.5) < 1e-12)


def test_image_orientation():
    dsc = np.array([[1.5e6, 0.0, 0.0]])  # viewer on +x, looking back at Earth
    los, up, right = G.image_axes(dsc)
    check("line of sight points from DSCOVR to Earth", np.allclose(los, [[-1, 0, 0]]))
    check("image up is Earth's north", np.allclose(up, [[0, 0, 1]]))
    check("image right is east (+y on the facing hemisphere)", np.allclose(right, [[0, 1, 0]]))
    target = np.array([[0.0, 1000.0, 2000.0]])  # east and north of Earth's centre
    x, y = G.sky_offset_rad(dsc, target)
    check("a target east and north of Earth appears right and up", x[0] > 0 and y[0] > 0,
          f"x={x[0]:.2e}, y={y[0]:.2e}")
    r = G.angular_radius_rad(np.array([1.5e6]), R_EARTH_EQ_KM)[0]
    check("Earth from 1.5 million km is 0.487 deg across", abs(2 * np.degrees(r) - 0.4872) < 2e-4,
          f"{2 * np.degrees(r):.4f} deg")


def test_sidereal_time_and_subpoints():
    g = np.degrees(G.gmst_rad(pd.to_datetime(["2000-01-01 12:00:00"])))[0]
    check("GMST at J2000.0 = 280.46062 deg (IAU 1982)", abs(g - 280.46061837) < 1e-6, f"{g:.8f}")
    t1 = pd.to_datetime(["2020-03-01 00:00:00", "2020-03-02 00:00:00"])
    dg = (np.degrees(np.diff(G.gmst_rad(t1)))[0]) % 360
    check("Earth turns 0.9856 deg more than once per solar day", abs(dg - 0.98565) < 1e-4, f"{dg:.5f}")
    gm = np.array([0.0])
    lat, lon = G.sub_point(np.array([[1.0, 0.0, 1.0]]), gm)
    check("sub-point of (1, 0, 1) at GMST 0 is 45N 0E", abs(lat[0] - 45) < 1e-9 and abs(lon[0]) < 1e-9)


# ---------------------------------------------------------------------------
# 2. Photometry
# ---------------------------------------------------------------------------

def _synthetic_disc(n=400, cx=200.5, cy=199.0, r=150.0, level=(80, 100, 120)):
    yy, xx = np.mgrid[:n, :n]
    rgb = np.zeros((n, n, 3), np.uint8)
    inside = np.hypot(xx - cx, yy - cy) < r
    for c in range(3):
        rgb[..., c][inside] = level[c]
    # a brighter "cloud" patch
    patch = inside & (xx > cx) & (yy > cy)
    rgb[patch] = 200
    return rgb, inside, patch


def test_measure_synthetic_disc():
    rgb, inside, patch = _synthetic_disc()
    m = ph.measure(rgb)
    check("disc pixel count", m["n_disc"] == inside.sum(), f"{m['n_disc']} vs {inside.sum()}")
    check("disc radius from extent within 1 px", abs(m["r_extent"] - 150.0) < 1.0, f"{m['r_extent']:.2f}")
    check("histograms count every disc pixel", int(m["hist"][0, 0].sum()) == inside.sum())
    v = rgb[..., 0][inside] / 255.0
    expect = float(np.mean(v ** 1.8))
    got = m[ph.moment_column("disc", "red", 1.8)]
    check("stored disc mean of (V/255)^1.8 is exact", abs(got - expect) < 1e-12, f"{got:.6f}")
    means = ph.disc_means(m["hist"][0][None], (np.arange(256) / 255.0) ** 1.0)[0]
    check("disc_means from histograms = direct means",
          np.allclose(means, [rgb[..., c][inside].mean() / 255 for c in range(3)]))


def test_transform_inverse():
    v = np.arange(256)
    for g in (1.4, 1.8, 2.4):
        back = T.to_display(T.to_linear(v, g), g)
        if not np.allclose(back, v, atol=1e-9):
            check(f"display -> linear -> display at gamma {g}", False)
            return
    check("display -> linear -> display is the identity", True)


# ---------------------------------------------------------------------------
# 3. Periodicity (simulated fixtures)
# ---------------------------------------------------------------------------

def _fake_light_curve(rng, days=120, per_day=15, period_h=24.0, noise=0.03):
    t = np.sort(rng.uniform(0, days, int(days * per_day)))
    y = 1 + 0.1 * np.sin(2 * np.pi * t * 24 / period_h) + 0.04 * np.sin(4 * np.pi * t * 24 / period_h + 1)
    return t, y + noise * rng.standard_normal(t.size)


def test_periodogram_recovers_day():
    rng = np.random.default_rng(1)
    t, y = _fake_light_curve(rng, period_h=23.9345)  # a sidereal day, to test resolution
    f, p = P.periodogram(t, y, 0.5, 1.5)
    fpk, _ = P.peak_near(f, p, 1.0, 0.05)
    check("Lomb-Scargle finds an injected 23.9345 h period", abs(24 / fpk - 23.9345) < 0.005,
          f"{24 / fpk:.4f} h")
    months = np.floor(t / 30).astype(int)
    jk = P.jackknife_peak(t, y, months, 1.0, 0.05)
    check("jackknife error is finite and small", 0 < jk["period_err_h"] < 0.01, f"{jk['period_err_h']:.5f} h")


def test_acf_and_sampling():
    rng = np.random.default_rng(2)
    t, y = _fake_light_curve(rng, days=60)
    grid_t, grid = P.hourly_grid(t, y)
    r = P.acf(grid, 24 * 4)
    lag, rv = P.acf_peak(np.arange(len(r)).astype(float), r)
    check("autocorrelation peaks at 24 h", abs(lag - 24.0) < 0.5, f"{lag:.2f} h, r = {rv:.2f}")
    res = P.sampling_test(t, y, 6483.0, [1, 12], rng, mode="random", n_segments=40)
    check("thinning test: dense sampling always finds the day", res[0]["p_detect"] > 0.95,
          f"{res[0]['p_detect']:.2f}")


def test_frame_jump():
    import lightcurve as LC
    t = pd.date_range("2024-06-01", periods=40, freq="90min")
    base = 0.4 + 0.02 * np.sin(np.arange(40) / 3)
    df = pd.DataFrame({"t": t, "D_red": base, "D_green": base, "D_blue": base})
    df.loc[20, ["D_red", "D_green", "D_blue"]] *= 1.08
    j = LC.frame_jump(df, pd.Series(True, index=df.index))
    check("a frame 8% brighter than its neighbours is flagged", abs(j[20] - 0.08) < 0.01 and
          np.nanmax(np.abs(j.drop(index=[19, 20, 21]))) < 0.01, f"jump {j[20]:+.3f}")


# ---------------------------------------------------------------------------
# 4. Surface and Moon
# ---------------------------------------------------------------------------

def test_lambert_weights():
    # Observer and Sun both over (0N, 0E): the weight cos(theta0) cos(theta)
    # integrated over the visible hemisphere; the disc's land share of a mask
    # that is land everywhere must be 1, and the polar share of light ~ the
    # analytic value for |lat| > 60 deg on a sphere seen pole-on... here edge-on.
    from config import LANDMASK_PATH
    if not LANDMASK_PATH.exists():
        print("  [skip] land mask not downloaded")
        return
    import surface
    fr = surface.disc_fractions(np.array([0.0]), np.array([0.0]), np.array([0.0]), np.array([0.0]))
    # analytic polar share: integral of cos^2 over |lat| > 60 / over the hemisphere
    lat = np.radians(np.linspace(-90, 90, 20001))
    lon = np.radians(np.linspace(-90, 90, 2001))
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    w = (np.cos(LA) * np.cos(LO)) ** 2 * np.cos(LA)
    expect = w[np.abs(LA) > np.radians(60)].sum() / w.sum()
    check("polar share of a Lambert sphere seen at full phase", abs(fr["polar_flux"][0] - expect) < 0.003,
          f"{fr['polar_flux'][0]:.4f} vs {expect:.4f}")
    check("land shares lie in 0..1", 0 <= fr["land_flux"][0] <= 1 and 0 <= fr["land_area"][0] <= 1)


def test_moon_bits():
    import moon
    t1 = moon.frame_time("epic_1b_20210211014104_03.png")
    t2 = moon.frame_time("197_2015197220604.png")
    check("gallery file names: new style", t1 == pd.Timestamp("2021-02-11 01:41:04"))
    check("gallery file names: 2015 day-of-year style", t2 == pd.Timestamp("2015-07-16 22:06:04"))
    import ephemeris as E
    t = pd.date_range("2030-01-01", periods=48, freq="30min")
    r_e = np.full(48, np.arctan(R_EARTH_EQ_KM / 1.5e6))
    r_m = np.full(48, np.arctan(R_MOON_KM / 1.2e6))
    x = np.linspace(-0.02, 0.02, 48)
    view = pd.DataFrame({"t": t, "x_rad": x, "y_rad": 0.0, "sep_rad": np.abs(x), "r_moon_rad": r_m,
                         "r_earth_rad": r_e, "moon_phase_deg": 9.0, "in_front": True, "phase_deg": 9.0})
    view["in_fov"] = (np.abs(x) < E.HALF_FOV_RAD + r_m)
    view["overlaps_earth"] = view["sep_rad"] < r_e + r_m
    ev = E.find_events(view)
    check("event finder: one transit through the centre", len(ev) == 1 and ev["kind"].iloc[0] == "transit"
          and ev["min_sep_earth_radii"].iloc[0] < 0.2)


# ---------------------------------------------------------------------------
# 5. Cached data
# ---------------------------------------------------------------------------

def test_cached_data():
    from config import CACHE_DIR
    path = CACHE_DIR / "metadata.parquet"
    if not path.exists():
        print("  [skip] metadata not built yet")
        return
    meta = pd.read_parquet(path)
    check("metadata: > 50,000 images since 2015-06-13", len(meta) > 50000 and str(meta["t"].min())[:10] == "2015-06-13",
          f"{len(meta)} images")
    dsc = meta[["dsc_x", "dsc_y", "dsc_z"]].to_numpy()
    sun = meta[["sun_x", "sun_y", "sun_z"]].to_numpy()
    a = G.phase_angle_deg(dsc, sun)
    check("phase angle stays within 1-13 deg", a.min() > 1 and a.max() < 13, f"{a.min():.2f}-{a.max():.2f}")
    lat, lon = G.sub_point(dsc, G.gmst_rad(meta["t"]))
    dlon = (lon - meta["lon"].to_numpy() + 180) % 360 - 180
    check("our sub-DSCOVR point matches the API's disc centre (lat within 0.3 deg)",
          np.percentile(np.abs(lat - meta["lat"].to_numpy()), 99) < 0.3)
    check("... and longitude within 0.7 deg (J2000 vs date frame)", np.percentile(np.abs(dlon), 99) < 0.7,
          f"median {np.median(dlon):+.2f} deg")


def main() -> int:
    groups = [
        ("Geometry", [test_phase_and_lambert, test_image_orientation, test_sidereal_time_and_subpoints]),
        ("Photometry", [test_measure_synthetic_disc, test_transform_inverse]),
        ("Periodicity (simulated fixtures)", [test_periodogram_recovers_day, test_acf_and_sampling,
                                              test_frame_jump]),
        ("Surface and Moon", [test_lambert_weights, test_moon_bits]),
        ("Cached data", [test_cached_data]),
    ]
    for title, tests in groups:
        print(title)
        for test in tests:
            test()
    print(f"\n{'All checks passed' if not FAILURES else f'{len(FAILURES)} FAILED: ' + ', '.join(FAILURES)}")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
