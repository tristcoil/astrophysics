"""Checks of the physics and the machinery, independent of the survey data.

    python -m pytest test_physics.py -q     (or: python test_physics.py)

Each test checks one piece against something known exactly: an analytic
result, a limiting case, or a published number reproduced from the
published inputs.
"""

from __future__ import annotations

import numpy as np

import cosmo_physics as cp
import environment as env
import galaxies as gx
import literature as lit
import sfh


# ---------------------------------------------------------------------------
# Cosmic star-formation history
# ---------------------------------------------------------------------------


def test_peak_is_where_the_derivative_vanishes():
    p = lit.MD14_FIT
    zp = sfh.z_peak(p["B"], p["C"], p["D"])
    eps = 1e-4
    lo, mid, hi = (sfh.psi_md(zp + d, **p) for d in (-eps, 0.0, eps))
    assert mid > lo and mid > hi
    # MD14 quote z ~ 1.9; the formula gives 1.86
    assert abs(zp - 1.863) < 0.002


def test_integral_converges_and_matches_md14():
    """rho_*(0) from MD14's own curve: their text says 5.8e8 M_sun/Mpc^3."""
    h = sfh.history(lit.MD14_FIT)
    assert abs(h["rho_star"][0] / 5.8e8 - 1) < 0.01
    # independent quadrature in redshift
    from scipy.integrate import quad
    c = sfh.COSMO
    f = lambda z: sfh.psi_md(z, **lit.MD14_FIT) / ((1 + z) * c.H(z).to_value("km/(s Mpc)")) \
        * 977.79e9
    direct = (1 - 0.27) * quad(f, 0, 30, limit=400)[0]
    assert abs(h["rho_star"][0] / direct - 1) < 1e-4


def test_md14_milestones():
    s = sfh.summarise(lit.MD14_FIT)
    assert abs(s["z_half_mass"] - 1.3) < 0.03
    assert abs(s["z_last_quarter"] - 0.7) < 0.01
    assert s["frac_before_reionization_end"] < 0.01
    assert abs(s["t_peak_gyr"] - 3.5) < 0.1


def test_fit_recovers_known_parameters():
    """Points drawn exactly from a known curve are refitted to it."""
    p = {"A": 0.02, "B": 2.5, "C": 3.1, "D": 6.0}
    z = np.linspace(0.05, 7.5, 40)
    pts = {"z": z, "logpsi": np.log10(sfh.psi_md(z, **p)),
           "err_plus": np.full(40, 0.1), "err_minus": np.full(40, 0.1)}
    f = sfh.fit_md(pts)
    for k in "ABCD":
        assert abs(f[k] / p[k] - 1) < 1e-4


def test_efold_of_pure_exponential():
    """A curve that is exactly exponential in time returns its own tau."""
    z = np.linspace(0.01, 0.99, 30)
    t = sfh.cosmic_time_gyr(z)
    logpsi = np.log10(np.exp(-t / 3.0))
    assert abs(sfh.efold_time_gyr(z, logpsi)["tau_gyr"] - 3.0) < 1e-6


# ---------------------------------------------------------------------------
# Distances, K-corrections, V_max
# ---------------------------------------------------------------------------


def test_low_redshift_hubble_law():
    z = 0.01
    assert abs(cp.comoving_distance(z) / (cp.C_KMS * z / 70.0) - 1) < 0.01


def test_k_correction_vanishes_at_zero_redshift():
    assert np.allclose(cp.k_correction("r", 0.0, np.array([0.3, 0.7, 1.0])), 0.0)


def test_vmax_limits_bracket_the_observed_redshift():
    """A galaxy seen at z_obs with r inside the window is visible at z_obs."""
    z_obs = np.array([0.05, 0.08, 0.12])
    gr = np.array([0.6, 0.9, 0.8])
    r = np.array([16.0, 17.2, 17.5])
    M = cp.absolute_magnitude(r, z_obs, cp.k_correction("r", z_obs, gr))
    zmin, zmax = cp.vmax_limits(M, gr, (0.01, 0.25), 17.6, 14.5)
    assert np.all(zmin < z_obs) and np.all(zmax > z_obs)
    # and at z_max the galaxy sits exactly at the faint limit
    m_at = M + cp.distance_modulus(zmax) + cp.k_correction("r", zmax, gr)
    assert np.allclose(m_at, 17.6, atol=1e-6)


# ---------------------------------------------------------------------------
# BPT classification
# ---------------------------------------------------------------------------


def _lines(x, y, sn=50.0):
    """Line fluxes with the given BPT ratios and uniform S/N (before scaling)."""
    n = len(x)
    ha = np.full(n, 100.0)
    hb = np.full(n, 30.0)
    g = {"h_alpha_flux": ha, "h_beta_flux": hb,
         "nii_6584_flux": ha * 10 ** np.asarray(x), "oiii_5007_flux": hb * 10 ** np.asarray(y)}
    for k in list(g):
        g[k.replace("_flux", "_flux_err")] = g[k] / sn
    return g


def test_bpt_regions():
    # well inside each region of the diagram
    x = np.array([-0.8, -0.25, 0.0, 0.1])
    y = np.array([-0.2, 0.0, 0.9, 0.0])
    cls = gx.bpt_classify(_lines(x, y))["bpt"]
    assert list(cls) == [gx.BPT_SF, gx.BPT_COMPOSITE, gx.BPT_SEYFERT, gx.BPT_LINER]


def test_bpt_weak_and_none():
    g = _lines([-0.5, -0.5], [0.0, 0.0], sn=50.0)
    g["h_beta_flux_err"] = g["h_beta_flux"] / 1.0  # H-beta undetected
    g["h_alpha_flux_err"][1] = g["h_alpha_flux"][1] / 1.0  # nothing detected
    cls = gx.bpt_classify(g)["bpt"]
    assert list(cls) == [gx.BPT_WEAK, gx.BPT_NONE]


# ---------------------------------------------------------------------------
# Ridge line and valley
# ---------------------------------------------------------------------------


def test_ridge_line_recovers_a_planted_sequence():
    rng = np.random.default_rng(1)
    n = 200_000
    logm = rng.uniform(9.0, 11.0, n)
    sf = rng.random(n) < 0.6
    logsfr = np.where(sf, 0.7 * logm - 7.0 + rng.normal(0, 0.3, n),
                      0.7 * logm - 9.0 + rng.normal(0, 0.4, n))
    ms = gx.ridge_line(logm, logsfr, np.ones(n))
    assert abs(ms["slope"] - 0.7) < 0.03
    assert abs(ms["intercept"] + 7.0) < 0.3
    v = gx.valley(logsfr - gx.ms_sfr(logm, ms), np.ones(n))
    assert -1.6 < v["position"] < -0.7


# ---------------------------------------------------------------------------
# Environment machinery
# ---------------------------------------------------------------------------


def test_disc_coverage_full_and_half():
    import healpy as hp
    nside = 256
    mask = np.zeros(hp.nside2npix(nside), bool)
    ra, dec = hp.pix2ang(nside, np.arange(len(mask)), lonlat=True)
    mask[(ra > 150) & (ra < 200) & (dec > 0) & (dec < 40)] = True
    fp = {"mask": mask, "nside": nside}
    full = env.disc_coverage(fp, np.array([175.0]), np.array([20.0]), np.array([2.0]))
    half = env.disc_coverage(fp, np.array([175.0]), np.array([0.0]), np.array([2.0]))
    assert full[0] > 0.99
    assert abs(half[0] - 0.5) < 0.08


def test_target_sampling_rate_of_a_collided_pair():
    # two targets 30" apart, one observed; an isolated observed target
    tgt = {"ra": np.array([180.0, 180.0 + 30 / 3600, 190.0]),
           "dec": np.array([10.0, 10.0, 10.0]),
           "observed": np.array([True, False, True])}
    tsr = env.target_sampling_rate(tgt)
    assert np.allclose(tsr, [0.5, 0.5, 1.0])


def test_nth_neighbour_on_a_lattice():
    """Galaxies on a plane perpendicular to the line of sight, unit spacing."""
    d = 200.0
    gx_, gy_ = np.meshgrid(np.arange(-10, 11), np.arange(-10, 11))
    xyz = np.column_stack([np.full(gx_.size, d), gx_.ravel(), gy_.ravel()]).astype(float)
    z = np.full(len(xyz), float(cp.redshift_at_distance(d)))
    centre = np.array([[d, 0.0, 0.0]])
    res = env.nth_neighbour_density(centre, z[:1], xyz, np.ones(len(xyz)), n_th=5, k=50)
    # neighbours at distance 1 (x4), sqrt(2) (x4): the 5th lies at sqrt(2)
    assert abs(res["d_n"][0] - np.sqrt(2)) < 1e-3


if __name__ == "__main__":
    import inspect
    import sys

    fails = 0
    for name, fn in sorted(inspect.getmembers(sys.modules[__name__], inspect.isfunction)):
        if name.startswith("test_"):
            try:
                fn()
                print(f"PASS {name}")
            except AssertionError as exc:
                fails += 1
                print(f"FAIL {name}: {exc}")
    sys.exit(1 if fails else 0)
