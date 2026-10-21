"""Checks of the physics and the machinery, independent of the survey data.

    python test_physics.py      (or: python -m pytest test_physics.py -q)

Each test checks one piece against something known exactly: an analytic
result, a limiting case, or a published number rebuilt from published inputs.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import config as C
import distances as D
import hosts as H
import lightcurves as LC
import progenitors as P
import rates as R


# ---------------------------------------------------------------------------
# Distances and volumes
# ---------------------------------------------------------------------------

def test_distance_modulus_matches_astropy():
    z = np.array([0.001, 0.01, 0.05, 0.1, 0.3])
    assert np.allclose(D.distance_modulus(z), D.cosmology().distmod(z).value, atol=2e-4)
    assert np.allclose(D.redshift_at_modulus(D.distance_modulus(z)), z, rtol=1e-3)


def test_cmb_frame_shift_toward_the_apex():
    """A galaxy straight toward the dipole apex is blueshifted by our motion:
    its CMB-frame velocity is ~370 km/s larger; opposite it, smaller."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    apex = SkyCoord(l=C.CMB_DIPOLE_L_DEG * u.deg, b=C.CMB_DIPOLE_B_DEG * u.deg, frame="galactic").icrs
    z = 0.02
    z_toward = D.cmb_frame_redshift(z, apex.ra.deg, apex.dec.deg)
    z_away = D.cmb_frame_redshift(z, (apex.ra.deg + 180) % 360, -apex.dec.deg)
    # (1+z)/(1-beta) - 1 and (1+z)/(1+beta) - 1: both ~ (1+z) * 369.8 km/s
    assert abs((z_toward - z) * C.C_KMS - 377.7) < 0.5
    assert abs((z - z_away) * C.C_KMS - 376.8) < 0.5


def test_euclidean_vmax():
    """Perley et al.'s D_max: M = -19.3 at the 18.5 mag limit -> 363 Mpc."""
    v = R.vmax_euclidean(-19.3, 18.5)
    d = (v * 3 / (4 * np.pi)) ** (1 / 3)
    assert abs(d - 10 ** (0.2 * 37.8 - 5)) < 1e-6 * d


def test_cosmological_vmax_is_smaller_and_consistent():
    M = np.array([-15.0, -17.0, -19.3, -21.5])
    zmax = R.zmax_cosmological(M)
    # the peak at z_max is exactly at the limit (with the BTS K-term)
    m_at = M + D.distance_modulus(zmax) - 2.5 * np.log10(1 + zmax)
    assert np.allclose(m_at, C.M_LIM, atol=2e-3)
    # comoving volume with time dilation < Euclidean volume of luminosity distance
    assert np.all(R.vmax_cosmological(M) < R.vmax_euclidean(M))
    # and they agree in the nearby limit
    assert abs(R.vmax_cosmological(-13.0) / R.vmax_euclidean(-13.0) - 1) < 0.01


def test_rate_weight_units():
    """One supernova, V_max = 1 Gpc^3, one year, nothing lost: 1 per Gpc^3 per yr."""
    s = pd.DataFrame({"x": [0]})
    w = R.rate_weights(s, M=np.array([-19.0]), vmax=np.array([1e9]), fcl=np.array([1.0]),
                       T_yr=1.0, f_sky=1.0, f_ext=1.0, f_rec=1.0)
    assert abs(w[0] - 1.0) < 1e-12


def test_rate_ratio_cancels_survey_factors():
    rng = np.random.default_rng(1)
    s = pd.DataFrame({"w": np.tile([1.0, 2.0, 3.0, 4.0], 25), "a": np.tile([True, True, False, False], 25)})
    r1 = R.rate_ratio(s, s["a"], ~s["a"], rng, n_boot=50)["ratio"]
    s2 = s.assign(w=s["w"] * 7.3)    # any common factor (f_sky, T, f_ext ...)
    r2 = R.rate_ratio(s2, s2["a"], ~s2["a"], rng, n_boot=50)["ratio"]
    assert abs(r1 - 3 / 7) < 1e-12 and abs(r1 - r2) < 1e-12


# ---------------------------------------------------------------------------
# Light curves
# ---------------------------------------------------------------------------

def test_lyman_bolometric_correction():
    """Lyman et al. (2014) Type II fit, BC_g = 0.053 - 0.089 x - 0.736 x^2."""
    assert abs(LC.bolometric_correction_g(0.0) - 0.053) < 1e-12
    assert abs(LC.bolometric_correction_g(1.0) - (0.053 - 0.089 - 0.736)) < 1e-12
    # outside the validity range the colour is clipped, not extrapolated
    assert LC.bolometric_correction_g(3.0) == LC.bolometric_correction_g(1.3)


def test_solar_luminosity_roundtrip():
    """M_bol = 4.74 is one solar luminosity (g - r = 0 has BC_g = 0.053)."""
    L = LC.lbol_from_gr(4.74 - 0.053, 4.74 - 0.053)
    assert abs(L / LC.L_SUN - 1) < 1e-9


def test_nickel_mass_inverts_decay_power():
    t = 150.0
    L = 0.07 * (LC.NI_E1 * np.exp(-t / LC.NI_TAU1) + LC.NI_E2 * np.exp(-t / LC.NI_TAU2))
    assert abs(LC.nickel_mass(L, t) - 0.07) < 1e-12
    # on the tail cobalt dominates: the luminosity falls 0.98 mag per 100 days
    dm = 2.5 * np.log10(LC.nickel_mass(1.0, 250.0) / LC.nickel_mass(1.0, 150.0))
    assert abs(dm - 0.976) < 0.01


def test_goldberg_inversion_matches_their_equation():
    """Solving Goldberg et al.'s two forward scalings must give their published
    inversion (their Eq. 21) to the rounding of its coefficients."""
    for L50, tp, mni, R_ in ((1.2e42, 110.0, 0.03, 500.0), (4e41, 85.0, 0.01, 800.0),
                             (2.5e42, 130.0, 0.08, 350.0)):
        E, M = LC.explosion_properties(L50, tp, mni, R_)
        l, t, n, r = np.log10(L50 / 1e42), np.log10(tp / 100), np.log10(mni), np.log10(R_ / 500)
        logE = -0.728 + 2.148 * l - 0.280 * n + 2.091 * t - 1.632 * r
        logM = -0.947 + 1.474 * l - 0.518 * n + 3.867 * t - 1.120 * r
        assert abs(np.log10(E) - logE) < 0.01
        assert abs(np.log10(M / 10) - logM) < 0.01


def test_fermi_function_shape():
    t = np.array([0.0, 100.0, 200.0])
    m = LC.fermi_mag(t, -17.0, 0.0, 2.0, 100.0, 3.0)
    assert abs(m[0] + 17.0) < 1e-6 and abs(m[1] + 16.0) < 1e-12 and abs(m[2] + 15.0) < 1e-6


def test_plateau_fit_recovers_a_synthetic_light_curve():
    """A noiseless plateau ending at day 104 is found to a fraction of a day."""
    t = np.arange(22.0, 160.0, 3.0)
    m = LC.fermi_mag(t, -16.8, 0.008, 1.9, 104.0, 4.0)
    fit = LC.fit_plateau_end(t, m, np.full_like(t, 0.03))
    assert fit is not None and fit["accepted"]
    assert abs(fit["t_pt"] - 104.0) < 0.5


# ---------------------------------------------------------------------------
# Hosts
# ---------------------------------------------------------------------------

def test_dlr_distance_round_galaxy():
    """For a round galaxy the DLR distance is separation / effective radius."""
    cand = pd.DataFrame({"fracDeV_r": [0.0], "deVRad_r": [9.0], "deVAB_r": [1.0], "deVPhi_r": [0.0],
                         "expRad_r": [4.0], "expAB_r": [1.0], "expPhi_r": [30.0],
                         "ra": [150.0], "dec": [0.0], "sn_ra": [150.0], "sn_dec": [6 / 3600],
                         "sep_arcsec": [6.0]})
    assert abs(H.dlr_distance(cand)[0] - 1.5) < 1e-9


def test_dlr_distance_along_minor_axis():
    """An edge-on disc (b/a = 0.25, PA 0 = major axis north-south): a SN 2"
    east lies along the minor axis, where the radius is 0.25 r_e."""
    cand = pd.DataFrame({"fracDeV_r": [0.0], "expRad_r": [8.0], "expAB_r": [0.25], "expPhi_r": [0.0],
                         "deVRad_r": [1.0], "deVAB_r": [1.0], "deVPhi_r": [0.0],
                         "ra": [10.0], "dec": [0.0], "sn_ra": [10.0 + 2 / 3600], "sn_dec": [0.0],
                         "sep_arcsec": [2.0]})
    assert abs(H.dlr_distance(cand)[0] - 1.0) < 1e-6


def test_zmatched_weights_follow_host_redshifts():
    gal_z = np.array([0.02, 0.02, 0.05, 0.05, 0.05])
    w = H.zmatched_weights(gal_z, np.ones(5), np.array([0.021, 0.051, 0.052, 0.053]))
    assert abs(w[gal_z < 0.03].sum() - 0.25) < 1e-12


# ---------------------------------------------------------------------------
# Progenitors
# ---------------------------------------------------------------------------

def test_smartt_counting_argument():
    """Smartt (2015): 70% of 8-100 M_sun stars are below 18 M_sun; with 30
    progenitors below, 13 should lie above; Poisson p(0 or 1 | 13) = 3e-5."""
    c = P.counting_argument(30, 1)
    assert abs(c["f_8_18_of_8_100"] - 0.688) < 0.001
    assert 13 <= c["n_expected_above"] < 14
    from scipy.stats import poisson
    assert abs(poisson.cdf(1, 13) - 3.2e-5) < 0.2e-5


def test_luminosity_slope_from_salpeter():
    """Davies & Beasor (2020): Gamma_L = (1 - alpha + Gamma)/alpha = -1.675."""
    assert abs(P.salpeter_gamma_L(2.0) + 1.675) < 1e-12


def test_db18_mass_luminosity_relation():
    """Their STARS fit: log M = -1.2417 + 0.47868 log L (e.g. log L 5.24 -> ~18 M_sun)."""
    assert abs(P.mass_from_logL(5.24) - 18.3) < 0.2


def test_smartt_likelihood_prefers_range_containing_detections():
    lo, hi = np.array([7.0, 9.0]), np.array([14.0, 30.0])
    lnL = P.smartt_likelihood([10.0], [2.0], [2.0], [], lo, hi)
    # a range containing the star beats a range that excludes it
    assert lnL[1, 0] > np.log(1e-299)
    lnL2 = P.smartt_likelihood([10.0], [2.0], [2.0], [], np.array([11.0]), np.array([30.0]))
    assert lnL[1, 0] > lnL2[0, 0]


def test_imf_counts():
    """Stars per unit mass formed: Salpeter (0.1-100) between 8 and 50 M_sun is
    0.00703, the textbook number behind 'one supernova per ~140 M_sun'."""
    from analysis import imf_number_per_mass
    assert abs(imf_number_per_mass(8, 50, "salpeter") - 0.00703) < 0.00002
    assert imf_number_per_mass(8, 50, "kroupa") > imf_number_per_mass(8, 50, "salpeter")


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
