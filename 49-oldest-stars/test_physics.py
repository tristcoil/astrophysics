"""Checks of the machinery: interpolation, extinction, cosmology, the fitter itself.

    python test_physics.py        # or: python -m pytest test_physics.py

The "injection" tests build a synthetic ridge line from an isochrone of known
age (a labelled simulation, not data) and check that the fitter recovers it.
"""

from __future__ import annotations

import numpy as np

import agefit as A
import catalog as K_
import config as C
import cosmology as CO
from fiducial import K, Ridge, turnoff_and_sgb_point
from isochrones import IsoGrid, Extinction, interpolate_in_age, load_grid


def _sgb_M(iso: dict) -> float:
    c = iso["F606W"] - iso["F814W"]
    ok = np.isfinite(c) & (iso["F606W"] < 6.5)
    return turnoff_and_sgb_point(c[ok], iso["F606W"][ok])[2]


def test_grid_nodes_reproduced():
    """Interpolating exactly at a grid node returns the node's own isochrone."""
    g = load_grid("MIST")
    i, j, k = list(g.feh).index(-2.0), list(g.alpha).index(0.4), 8
    age = 10 ** g.log_age[k] / 1e9
    iso = g.interpolate(-2.0, 0.4, [age])
    ref = g.data["F606W"][i, j, k]
    ok = np.isfinite(ref)
    assert np.nanmax(abs(iso["F606W"][0][ok] - ref[ok])) < 1e-4


def test_feh_interpolation_error():
    """Quadratic [Fe/H] interpolation over 0.5 dex steps, checked on a held-out MIST node.

    MIST's -1.75 node is predicted from the -2.5, -2.0 and -1.5 nodes only (the
    0.5 dex spacing MIST and DSED have where the metal-poor clusters sit). The
    subgiant point must agree to better than 0.02 mag (~0.3 Gyr).
    """
    g = load_grid("MIST")
    keep = [list(g.feh).index(f) for f in (-2.5, -2.0, -1.5)]
    sub = IsoGrid("sub", g.feh[keep], g.alpha, g.log_age, g.eep, {q: v[keep] for q, v in g.data.items()},
                  feh_order=2, eep_tip=g.eep_tip)
    lin = IsoGrid("lin", sub.feh, sub.alpha, sub.log_age, sub.eep, sub.data, feh_order=1, eep_tip=sub.eep_tip)
    for age in (11.0, 13.0):
        true = _sgb_M({q: v[0] for q, v in g.interpolate(-1.75, 0.4, [age]).items()})
        d = abs(_sgb_M({q: v[0] for q, v in sub.interpolate(-1.75, 0.4, [age]).items()}) - true)
        d_lin = abs(_sgb_M({q: v[0] for q, v in lin.interpolate(-1.75, 0.4, [age]).items()}) - true)
        print(f"      {age:g} Gyr: subgiant point off by {d:.4f} mag (quadratic), {d_lin:.4f} mag (linear)")
        assert d < 0.02, d


def test_age_interpolation_error():
    """Interpolating BaSTI between 12 and 14 Gyr reproduces its own 13 Gyr isochrone."""
    g = load_grid("BaSTI")
    ages = 10 ** g.log_age / 1e9
    comp = g.at_composition(-2.2, 0.4)
    pick = [int(np.argmin(abs(ages - a))) for a in (12.0, 14.0)]
    approx = interpolate_in_age({q: v[pick] for q, v in comp.items()}, g.log_age[pick], [13.0])
    true = interpolate_in_age(comp, g.log_age, [13.0])
    d = abs(_sgb_M({q: v[0] for q, v in approx.items()}) - _sgb_M({q: v[0] for q, v in true.items()}))
    print(f"      12-14 Gyr bracket: subgiant point off by {d:.4f} mag")
    assert d < 0.01, d


def test_extinction_ratios():
    """A_X / E(B-V) for a turnoff star, against standard R_V = 3.1 extinction laws.

    For ACS F606W the Cardelli et al. (1989) and Fitzpatrick (1999) laws give
    ~2.8-2.9 (Schlafly & Finkbeiner 2011 Table 6, rescaled to true E(B-V)), and
    ~1.6-1.9 for F814W. Also checks that a cool giant is reddened less in F606W.
    """
    ext = Extinction(-2.0, 0.4)
    a606, a814 = ext.ratio("F606W", 3.79, 4.2), ext.ratio("F814W", 3.79, 4.2)
    print(f"      turnoff star: A_F606W/E = {a606:.3f}, A_F814W/E = {a814:.3f}, colour {a606 - a814:.3f}")
    assert 2.7 < a606 < 3.0, a606
    assert 1.55 < a814 < 1.95, a814
    assert ext.ratio("F606W", 3.65, 2.0) < a606


def test_age_formula_matches_astropy():
    from astropy.cosmology import FlatLambdaCDM

    for h0, om in ((67.4, 0.315), (73.5, 0.332), (70.0, 1.0)):
        ref = FlatLambdaCDM(H0=h0, Om0=om, Tcmb0=0.0).age(0).value
        assert abs(CO.age_flat_lcdm(h0, om) - ref) < 1e-3, (h0, om)
    z = 11.0
    ref = FlatLambdaCDM(H0=67.4, Om0=0.315, Tcmb0=0.0).age(z).value
    assert abs(CO.age_at_redshift(z, 67.4, 0.315) - ref) < 1e-3


def test_external_tables():
    h = K_.harris()
    assert h["NGC 6341"]["feh"] == -2.31 and h["NGC 104"]["ebv"] == 0.04
    b = K_.bv21()
    assert abs(b["NGC 6341"]["mean"][0] - 8.501) < 1e-6
    assert b["NGC 7078"]["gaia_parallax"] is None


def test_quality_mask_keeps_bright_stars():
    """Regression: stars measured in one short exposure (RMS = 0) must survive the cuts."""
    cat = K_.read_hugs("ngc6341")
    ok = K_.quality_mask(cat)
    m = cat["m_F606W"]
    sgb = (m > 16) & (m < 18) & (cat["memb"] >= 90)
    assert ok[sgb].mean() > 0.5, ok[sgb].mean()


def _synthetic_ridge(grid_name: str, age: float, feh: float, dm: float, ebv: float, noise: float, seed: int):
    """A ridge line drawn from an isochrone of known age (simulation), with noise."""
    g = load_grid(grid_name)
    f = A.ClusterFitter.__new__(A.ClusterFitter)  # only to borrow models()
    f.grid, f.ext, f.thin = g, Extinction(feh, 0.4), 1
    f.ages, f.m_lo, f.m_hi = np.array([age]), -10, 40
    m, c, _ = A.ClusterFitter.models(f, A.Draw(feh, 0.4, dm, ebv))
    m, c = m[0], c[0]
    ok = np.isfinite(m) & np.isfinite(c)
    m, c = m[ok], c[ok]
    c_to, m_to, m_sgb = turnoff_and_sgb_point(c, m)
    keep = (m < m_to + 3.0) & (m > m_sgb - 2.5)
    m, c = m[keep], c[keep]
    # resample evenly in arc length, as the data ridge is
    x, y = K * c, m
    s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))])
    si = np.arange(0, s[-1], 0.06)
    xi, yi = np.interp(si, s, x), np.interp(si, s, y)
    dx, dy = np.gradient(xi), np.gradient(yi)
    n = np.hypot(dx, dy)
    nx, ny = -dy / n, dx / n
    rng = np.random.default_rng(seed)
    t = noise * rng.standard_normal(len(si))
    ci, mi = (xi + t * nx) / K, yi + t * ny
    i_to = int(np.argmin(ci))
    return Ridge(s=si, c=ci, m=mi, sig=np.full(len(si), noise), nx=nx, ny=ny, n=np.full(len(si), 100),
                 i_to=i_to, width=np.full(len(si), 0.02))


def test_injection_same_model():
    """A MIST-made ridge of 12.0 Gyr, fitted with MIST, returns 12.0 Gyr."""
    for true_age in (11.0, 13.5):
        r = _synthetic_ridge("MIST", true_age, -2.2, 14.6, 0.03, 0.003, seed=1)
        fit = A.ClusterFitter(load_grid("MIST"), r, feh0=-2.2)
        b = fit.best(A.Draw(-2.2, 0.4, 14.6, 0.03))
        print(f"      true {true_age} Gyr -> fitted {b['age']:.3f} Gyr, colour shift {b['dc']:+.4f}")
        assert abs(b["age"] - true_age) < 0.15, (true_age, b)
        assert abs(b["dc"]) < 0.005, b


def test_injection_colour_offset():
    """A colour error of the models (0.03 mag) is absorbed by the shift, not the age."""
    r = _synthetic_ridge("MIST", 12.5, -2.0, 14.6, 0.03, 0.003, seed=2)
    r = Ridge(r.s, r.c + 0.03, r.m, r.sig, r.nx, r.ny, r.n, r.i_to, r.width)
    b = A.ClusterFitter(load_grid("MIST"), r, feh0=-2.0).best(A.Draw(-2.0, 0.4, 14.6, 0.03))
    print(f"      true 12.5 Gyr, models 0.03 mag too blue -> {b['age']:.3f} Gyr, shift {b['dc']:+.4f}")
    assert abs(b["age"] - 12.5) < 0.25 and abs(b["dc"] - 0.03) < 0.01, b


def test_injection_distance_sensitivity():
    """Fitting with a distance modulus 0.1 mag too large makes the cluster ~1 Gyr younger."""
    r = _synthetic_ridge("MIST", 13.0, -2.2, 14.6, 0.03, 0.003, seed=3)
    f = A.ClusterFitter(load_grid("MIST"), r, feh0=-2.2)
    right = f.best(A.Draw(-2.2, 0.4, 14.6, 0.03))["age"]
    far = f.best(A.Draw(-2.2, 0.4, 14.7, 0.03))["age"]
    print(f"      distance modulus +0.1 mag: {right:.2f} -> {far:.2f} Gyr ({right - far:.2f} Gyr younger)")
    assert 0.8 < right - far < 1.8, (right, far)


def test_injection_cross_model():
    """The model-set systematic: a 13.0 Gyr DSED cluster, fitted with MIST and BaSTI.

    Not a pass/fail check of the code but a measurement recorded in the README:
    DSED's subgiant branch is brighter at a given age, so MIST and BaSTI read
    a DSED-made cluster as ~0.5 Gyr younger.
    """
    r = _synthetic_ridge("DSED", 13.0, -2.2, 14.6, 0.03, 0.003, seed=4)
    for model in ("MIST", "BaSTI"):
        age = A.ClusterFitter(load_grid(model), r, feh0=-2.2).best(A.Draw(-2.2, 0.4, 14.6, 0.03))["age"]
        print(f"      13.0 Gyr DSED cluster fitted with {model}: {age:.2f} Gyr")
        assert 12.2 < age < 12.8, (model, age)
    same = A.ClusterFitter(load_grid("DSED"), r, feh0=-2.2).best(A.Draw(-2.2, 0.4, 14.6, 0.03))["age"]
    assert abs(same - 13.0) < 0.1, same


if __name__ == "__main__":
    import sys
    import time

    failed = 0
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            t0 = time.time()
            try:
                fn()
                print(f"PASS  {name}  ({time.time() - t0:.1f} s)")
            except AssertionError as exc:
                failed += 1
                print(f"FAIL  {name}: {exc}")
    sys.exit(1 if failed else 0)
