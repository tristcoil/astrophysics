"""Checks of the physics and numerics against answers we know.

    python test_physics.py          (or: python -m pytest -q test_physics.py)

Most tests run from the cached data; the simulation checks read the cached runs
(data/derived) and are skipped if those do not exist yet.
"""

from __future__ import annotations

import numpy as np

import apollo
import config as C
import earth
import ephem
import gr_predictions as G
import literature as L
import llr_model


# --- data readers ------------------------------------------------------------

def test_apollo_parser_reproduces_documented_example():
    """The worked example on APOLLO's normal-point page, field by field."""
    r = apollo.parse_line("512008 330123950000000026170710379889370610207   317312B 72439 -5134 5320A 250A")
    assert (r["year"], r["month"], r["day"], r["hour"], r["minute"]) == (2008, 3, 30, 12, 39)
    assert r["sec_1e7"] == 500_000_000 and r["rt_1e13"] == 26_170_710_379_889
    assert r["reflector"] == 3 and r["station"] == 70610 and r["photons"] == 207
    assert r["sigma_01ps"] == 317 and r["snr_x10"] == 312 and r["grade"] == "B"
    assert r["p_001mbar"] == 72439 and r["t_01c"] == -51 and r["rh_pct"] == 34 and r["wavelength_A"] == 5320
    assert r["duration_s"] == 250


def test_apollo_files_complete():
    d = apollo.load()
    assert len(d) == 3846
    assert d.groupby("group").size().to_dict() == {"a": 941, "b": 514, "c": 361, "d": 760, "e": 915, "f": 355}
    assert 2.32 < d.rt_s.min() and d.rt_s.max() < 2.69


def test_ortho_eop_matches_iers_test_case():
    """IERS 2010 ORTHO_EOP.F documented test case (MJD 47100)."""
    dx, dy, du = earth.ortho_eop(47100.0)
    assert abs(dx[0] - (-162.8386373279636530)) < 1e-9
    assert abs(dy[0] - 117.7907525842668974) < 1e-9
    assert abs(du[0] - (-23.39092370609808214)) < 1e-9


def test_apollo_station_geodetic():
    r = earth.geocentric_to_xyz(C.APOLLO_GEOCENTRIC["radius_km"] * 1e3, C.APOLLO_GEOCENTRIC["lat_deg"],
                                C.APOLLO_GEOCENTRIC["lon_deg"])
    phi, lam, h = earth.geodetic(r)
    assert abs(np.degrees(phi) - 32.7803) < 1e-3 and abs(np.degrees(lam) + 105.8204) < 1e-3
    assert 2700 < h < 2900                                    # Apache Point is ~2,788 m up


# --- ephemeris evaluation ---------------------------------------------------------

def test_exact_time_evaluation_is_split_invariant():
    """(day, sec) and (day - 1, sec + 86400) must give the same Moon to < 1 micrometre."""
    e = ephem.default()
    p1, _ = e.posvel(ephem.MOON, 7000, 1234.5678, center=ephem.EARTH)
    p2, _ = e.posvel(ephem.MOON, 6999, 1234.5678 + 86400.0, center=ephem.EARTH)
    assert np.abs(p1 - p2).max() < 1e-9                    # km


def test_naif_spk_equals_jpl_ascii_de440():
    path = C.RAW_DIR / "jpl_ascii" / "ascp01950.440"
    if not path.exists():
        return
    recs = ephem.load_de_ascii(path)
    de = ephem.DEAscii(recs)
    e = ephem.default()
    day, sec = np.array([2300, 5000, 7300]), np.array([1234.5, 40000.0, 86000.0])
    a, _ = de("moon", day, sec)
    b, _ = e.posvel(ephem.MOON, day, sec, center=ephem.EARTH)
    assert np.abs(a - b).max() < 1e-7                      # km: 0.1 mm
    la, _ = de("librations", day, sec)
    lb, _ = e.moon_angles(day, sec)
    assert np.abs(la - lb).max() < 1e-12


def test_chebyshev_derivative():
    e = ephem.default()
    p1, _ = e.posvel(ephem.MOON, 6000, 4999.5, center=ephem.EARTH)
    p2, _ = e.posvel(ephem.MOON, 6000, 5000.5, center=ephem.EARTH)
    _, v = e.posvel(ephem.MOON, 6000, 5000.0, center=ephem.EARTH)
    assert np.abs((p2 - p1) - v).max() < 1e-8              # km/s


# --- range-model ingredients ----------------------------------------------------------

def test_zenith_delay_one_atmosphere():
    """W&B 2020: 'At the zenith and on the geoid, one atmosphere of air delays the range by 2.4 m.'"""
    zhd, zwd = llr_model.zenith_delay_m(np.radians(45.0), 0.0, 1013.25, 0.0, 0.532)
    assert 2.35 < zhd < 2.45 and zwd == 0.0


def test_fcula_mapping():
    assert abs(llr_model.fcula(np.radians(45.0), 0.0, 15.0, np.radians(90.0)) - 1.0) < 1e-12
    m30 = llr_model.fcula(np.radians(45.0), 0.0, 15.0, np.radians(30.0))
    assert 1.98 < m30 < 2.0                                # ~1/sin(30 deg), slightly less (curvature)


def test_water_vapour_pressure():
    e = llr_model.water_vapour_hpa(100.0, 20.0, 1013.25)
    assert 23.0 < e < 23.6                                 # saturation at 20 C: 23.4 hPa


def test_solar_shapiro_delay_for_lunar_ranging():
    """W&B 2020: the Sun's mean Shapiro delay for LLR is ~50 ns round trip (~7.5 m one way)."""
    rho = 384_400.0
    r = C.AU / 1e3
    dt = llr_model.shapiro_s(C.GM_SUN, r, r, rho)
    assert 7.3 < dt * C.C_LIGHT < 7.8


def test_lunar_tide_constant_part_magnitude():
    """The static tide raised by the Earth at Apollo 15 is ~0.4 m up (W&B 2020, Table A.2: 383 mm)."""
    site = np.array([C.REFLECTORS_PA["Apollo 15"]])
    earth_pa = np.array([[385_000e3, 0.0, 0.0]])
    d = llr_model.tide_on_body(site, [(earth_pa, C.GM_EARTH, True)], C.GM_MOON, C.R_MOON, C.MOON_LOVE)
    up = float((d @ site.T).item() / np.linalg.norm(site))
    assert 0.33 < up < 0.43


# --- GR's predictions -------------------------------------------------------------------

def test_classic_predictions():
    assert abs(G.mercury_perihelion()["value"] - 42.98) < 0.01
    assert abs(G.light_deflection_limb()["value"] - 1.7512) < 0.0005
    assert abs(G.lunar_geodetic_precession()["value"] - 19.2) < 0.05
    assert abs(G.redshift_tower()["value"] - 2.46e-15) < 0.01e-15
    assert abs(G.lageos_node_lense_thirring()["value"] - 30.7) < 0.2
    assert abs(G.redshift_sun()["seen_from_earth_kms"] * 1e3 - 633.3) < 0.5


def test_hulse_taylor_pk_parameters():
    """Weisberg & Huang (2016): masses 1.438 + 1.390 give Pb_dot(GR) = -2.40263e-12 and omega_dot 4.2266 deg/yr."""
    pk = G.pk_parameters(1.438, 1.390, 0.322997448918, 0.6171340)
    assert abs(pk["pbdot"] / -2.40263e-12 - 1) < 5e-4
    assert abs(pk["omega_dot_deg_yr"] - 4.226585) < 1e-3
    assert abs(pk["gamma_s"] - 0.004307) < 5e-6


def test_double_pulsar_masses_fit_all_pk_parameters():
    """Kramer et al. (2021) masses reproduce gamma_E, Pb_dot and s to their quoted precision."""
    pv = L.pulsar_value
    ma, mb = pv("dp_k21_mA")[0], pv("dp_k21_mB")[0]
    pk = G.pk_parameters(ma, mb, pv("dp_k21_pb")[0], pv("dp_k21_eT")[0], pv("dp_k21_x")[0])
    assert abs(pk["gamma_s"] * 1e3 - pv("dp_k21_gammaE")[0]) < 3 * pv("dp_k21_gammaE")[1]
    assert abs(pk["pbdot"] / pv("dp_k21_pbdot_gw")[0] - 1) < 2e-4
    assert abs(pk["s"] - pv("dp_k21_sini")[0]) < 5e-5


# --- literature conversions -------------------------------------------------------------

def test_one_sigma_conversions():
    v, e = L.one_sigma("wep", "bp1972")                    # 95% interval
    assert abs(e - 0.9e-12 / 1.96) < 1e-20
    v, e = L.one_sigma("wep", "microscope2022")            # stat and syst in quadrature
    assert abs(e - np.hypot(2.3e-15, 1.5e-15)) < 1e-22
    t = L.gamma_timeline().set_index("id")
    assert abs(t.loc["E1919-sobral-4inch", "sigma"] - 0.12 * 1.4826 / 1.7512) < 1e-6   # probable error


# --- simulations (cached) ------------------------------------------------------------------

def test_nordtvedt_coefficient():
    import nordtvedt
    if not nordtvedt.CACHE.exists():
        return
    a = nordtvedt.amplitude()
    assert 12.7 < a["cosD"] < 13.4                         # literature: 12.9-13.3 m per unit eta
    assert a["nonlinearity_max_m"] < 0.01                  # +eta and -eta runs mirror each other


def test_triple_pulsar_signature():
    import triple_pulsar
    if not triple_pulsar.CACHE.exists():
        return
    s = triple_pulsar.signature()
    assert 28 < s["amplitude_ns_per_1e-6"] < 38            # Archibald et al. 2018: 30 ns
    assert s["nonlinearity_ns_per_1e-6"] < 1.0


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
