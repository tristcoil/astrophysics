"""Checks of the physics and the parsers against known answers.

    python test_physics.py      (or: python -m pytest test_physics.py)
"""

from __future__ import annotations

import numpy as np

import blackcat
import gaia_orbits
import graveyard
import masses
import outbursts
from config import AU_M, DAY_S, GM_SUN


def test_mass_function_earth_sun():
    # Earth around the Sun: K = 29.78 km/s, P = 1 yr -> f(M) = M_sun^3/(M_sun + M_earth)^2 ~ 1
    f = masses.mass_function(365.25636, 29.7827)
    assert abs(f - 1.0) < 2e-3, f


def test_mass_function_v404_cyg():
    # Casares & Charles 1994: P = 6.4714 d, K2 = 208.5 km/s -> f(M) = 6.08 M_sun
    assert abs(masses.mass_function(6.4714, 208.5) - 6.08) < 0.01


def test_mass_round_trip():
    rng = np.random.default_rng(0)
    m1, m2, incl = rng.uniform(3, 30, 50), rng.uniform(0.2, 2, 50), rng.uniform(15, 89, 50)
    porb = rng.uniform(0.2, 30, 50)
    a = (GM_SUN * (m1 + m2) * (porb * DAY_S / (2 * np.pi)) ** 2) ** (1 / 3)
    k2 = 2 * np.pi * a * m1 / (m1 + m2) * np.sin(np.radians(incl)) / (porb * DAY_S) / 1e3
    fm = masses.mass_function(porb, k2)
    assert np.allclose(masses.m1_from_fm(fm, m2 / m1, incl), m1, rtol=1e-9)
    assert np.allclose(masses.solve_m_bh(fm, m2, incl), m1, rtol=1e-6)


def test_cyg_x1_central_values():
    # Miller-Jones et al. 2021 medians: K = 75.21, e = 0.0189, i = 27.51, M_opt = 40.6 -> ~21 M_sun
    f = masses.mass_function(5.599829, 75.21, 0.0189)
    assert 20.5 < masses.solve_m_bh(f, 40.6, 27.51) < 22.0


def test_eggleton():
    assert abs(masses.roche_lobe_fraction(1.0) - 0.3789) < 1e-3
    # a low-mass donor (q = 0.067, V404 Cyg) eclipses only above ~80 degrees
    assert 79 < masses.max_inclination_no_eclipse(0.067) < 81


def test_thiele_innes_round_trip():
    rng = np.random.default_rng(1)
    a0, i = rng.uniform(0.1, 5, 200), rng.uniform(1, 179, 200)
    w, W = rng.uniform(-179, 179, 200), rng.uniform(0, 179.9, 200)
    A, B, F, G = gaia_orbits.thiele_innes(a0, i, w, W)
    a0b, ib, wb, Wb = gaia_orbits.campbell(A, B, F, G)
    assert np.allclose(a0b, a0) and np.allclose(ib, i, atol=1e-7)
    assert np.allclose(np.cos(np.radians(wb - w)), 1) and np.allclose(Wb, W, atol=1e-7)


def test_astrometric_mass_function():
    # dark companion M2 around M1: the star's orbit is a M2/(M1+M2)
    m1, m2, p_yr, plx = 0.93, 9.6, 185.6 / 365.25, 2.1
    a_au = ((m1 + m2) * p_yr**2) ** (1 / 3)
    a0 = a_au * m2 / (m1 + m2) * plx
    f = masses.astrometric_mass_function(a0, plx, p_yr * 365.25)
    assert abs(f - m2**3 / (m1 + m2) ** 2) < 1e-9
    assert abs(masses.solve_m2_astrometric(f, m1) - m2) < 1e-6


def test_population_models_normalised():
    grid = np.linspace(0.01, 200, 400001)
    for model, theta in ((masses.PowerLaw, [5.0, 15.0, -2.0]), (masses.PowerLaw, [4.0, 12.0, -1.0]),
                         (masses.Exponential, [5.5, 1.5]), (masses.Gaussian, [7.8, 1.2])):
        t = np.array([theta])
        p = model.pdf(grid, t)[0]
        assert abs(np.trapezoid(p, grid) - 1) < 1e-3, model.name
        q01 = model.quantile(t, 0.01)[0]
        cdf = np.trapezoid(p[grid <= q01], grid[grid <= q01])
        assert abs(cdf - 0.01) < 1e-3, (model.name, cdf)


def test_population_cdfs_match_pdfs():
    g = np.linspace(0.3, 100, 200001)
    for model, theta in ((masses.PowerLaw, [5.0, 15.0, -2.0]), (masses.Exponential, [5.5, 1.5]),
                         (masses.Gaussian, [7.8, 1.2])):
        t = np.array([theta])
        p = model.pdf(g, t)[0]
        num = np.concatenate([[0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(g))])
        assert np.max(np.abs(num - model.cdf(g, t)[0])) < 2e-4, model.name


def test_population_likelihood_finite_for_spikes():
    # An exponential collapsing onto one mass must give the systems' densities
    # there (log N(6.9; 7, 0.3) + log N(6.9; 7.5, 0.4) = -0.90), not diverge.
    rng = np.random.default_rng(2)
    samples = {"a": rng.normal(7.0, 0.3, 400_000), "b": rng.normal(7.5, 0.4, 400_000)}
    like = masses.PopulationLikelihood(samples, masses.Exponential)
    val = like(np.array([[6.9, 1e-7]]))[0]
    assert abs(val - (-0.898)) < 0.05, val


def test_crab_energy_flux():
    # the canonical Crab gives ~2.2e-8 erg/cm2/s in 2-10 keV
    assert 2.0e-8 < outbursts.crab_energy_flux(2, 10) < 2.4e-8


def test_veiling_bias():
    assert masses.veiling_mass_bias(0.0) == 1.0
    assert abs(masses.veiling_mass_bias(0.3) - 0.7**-1.5) < 1e-12


def test_blackcat_value_parsing():
    p = blackcat.parse_value("26.79377±0.00007")
    assert p["kind"] == "pm" and p["value"] == 26.79377 and p["plus"] == 7e-5
    p = blackcat.parse_value("7.752340±2E-06")
    assert abs(p["plus"] - 2e-6) < 1e-15
    p = blackcat.parse_value(">67 – 81 (or 63±3 for the orientation angle of the radio jet)")
    assert p["kind"] == "range" and p["flag"] == ">" and (p["low"], p["high"]) == (67, 81)
    p = blackcat.parse_value("9 +0.2 -0.6")
    assert p["kind"] == "asym" and (p["plus"], p["minus"]) == (0.2, 0.6)
    assert blackcat.parse_value("≤0.053")["flag"] == "≤"
    assert blackcat.parse_value("≈0.03")["flag"] == "≈"
    # a literal '<' is a limit, not an HTML tag
    assert blackcat.clean_cell("<70<br>67&plusmn;4") == ["<70", "67±4"]
    assert blackcat.clean_cell("8.40&plusmn;<sup>0.77</sup><sub>0.70</sub>") == ["8.40 +0.77 -0.70"]


def test_freire_value_parsing():
    assert graveyard.parse_freire_value("2.01(4)") == (2.01, 0.04, 0.04)
    v, up, down = graveyard.parse_freire_value("1.4380(+12/&minus;11)")
    assert v == 1.438 and abs(up - 0.0012) < 1e-12 and abs(down - 0.0011) < 1e-12
    assert graveyard.parse_freire_value("1.68+(+17/&minus;15)")[0] == 1.68


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"{len(tests)} tests passed")
