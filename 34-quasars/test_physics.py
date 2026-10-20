"""Checks of the physics and the machinery.

    python -m pytest test_physics.py -q     (or: python test_physics.py)

Each test checks one piece against something known exactly: an analytic
result, a limiting case, a planted signal, or a published number reproduced
from the published inputs. The last few need the cached data (run
fetch_data.py first) and are skipped without it.
"""

from __future__ import annotations

import gzip
import shutil
import tempfile
from pathlib import Path

import numpy as np

import composite
import fitsstream
import literature as lit
import masses
import photometry
import physics
import qlf
from config import COMPOSITE_GRID, DR7_TABLE, LIT_DIR, WUSHEN_TABLE


# ---------------------------------------------------------------------------
# Eddington-limited growth and cosmology
# ---------------------------------------------------------------------------

def test_eddington_luminosity_and_time():
    # 4 pi G m_p c / sigma_T = 1.26e38 erg/s per solar mass; sigma_T c / (4 pi G m_p) = 450 Myr
    assert abs(physics.eddington_per_msun() / 1.2571e38 - 1) < 1e-3
    assert abs(physics.eddington_time_myr() - 450.5) < 0.5
    # e-folding every ~50 Myr at 10% efficiency (the atlas spec: "~50 Myr")
    assert abs(physics.efold_time_myr(0.1) - 50.05) < 0.05


def test_seed_growth_round_trip():
    # 100 Msun -> 1e9 Msun takes ln(1e7) e-folds: 16.1 x 50 Myr = 807 Myr ("~0.8 Gyr")
    t = physics.time_to_grow_myr(100.0, 1e9)
    assert abs(t - 806.8) < 1.0
    assert abs(physics.seed_mass(1e9, t) - 100.0) < 1e-6
    # half the Eddington ratio, twice the time
    assert abs(physics.time_to_grow_myr(100.0, 1e9, lam=0.5) / t - 2.0) < 1e-9


def test_cosmology_ages():
    # D3's fit quotes 13.825 +/- 0.033 Gyr; astropy with D3's rounded (H0, Om) gives 13.804, inside
    # that error. The reference (0.7, 0.3) model without radiation gives 13.47 Gyr.
    assert abs(physics.age_myr(0.0) / 1e3 - 13.825) < 0.033
    assert abs(physics.age_myr(0.0, "ref") / 1e3 - 13.467) < 0.01
    # J0313-1806 at z = 7.6423: "~670 Myr after the Big Bang" (Wang et al. 2021)
    assert abs(physics.age_myr(7.6423) - 670.0) < 15.0
    # round trip
    assert abs(physics.z_at_age(physics.age_myr(3.0)) - 3.0) < 1e-3


def test_luminosity_2500_from_mi():
    # Richards et al. 2006 eq. 4: M_i(z=2) = -27.2 (3C 273) -> log L_nu(2500) ~ 31.1 erg/s/Hz
    l = physics.luminosity_2500(-27.2)
    expected = -0.4 * (-27.2 + 48.6 + 2.5 * np.log10(3)) + np.log10(4 * np.pi * (3.0857e19) ** 2)
    assert abs(l - expected) < 1e-3
    assert 30.9 < l < 31.3


# ---------------------------------------------------------------------------
# Synthetic photometry and K-corrections
# ---------------------------------------------------------------------------

def test_flat_fnu_spectrum_has_zero_colours():
    w = np.arange(2000.0, 12000.0, 1.0)
    flam = 1.0 / w**2  # f_nu constant
    mags = [photometry.ab_mag(w, flam, b) for b in "ugriz"]
    assert np.ptp(mags) < 1e-6


def test_power_law_k_correction_is_analytic():
    """For f_nu ~ nu^alpha the full K-correction is -2.5(1+alpha)log10((1+z)/3)."""
    w = np.arange(500.0, 12000.0, 0.5)
    for alpha in (-0.5, -1.0, 0.3):
        flam = w ** (-(alpha + 2.0))
        z = np.array([0.5, 1.0, 2.0, 3.0])
        k = photometry.kcorr_total(w, flam, z, "i", igm=False)
        assert np.max(np.abs(k - photometry.kcorr_continuum(z, alpha))) < 2e-3


def test_k_correction_pegged_at_two():
    assert abs(photometry.kcorr_continuum(2.0)) < 1e-12
    # Richards et al. eq. 1: M_i(z=0) - M_i(z=2) = 0.596 for alpha = -0.5
    assert abs(photometry.kcorr_continuum(0.0) - 0.596) < 1e-3


def test_madau_transmission_limits():
    z = 3.0
    lam_red = 1215.67 * (1 + z) * 1.01
    assert photometry.madau95_transmission(np.array([lam_red]), z)[0] == 1.0
    # between Ly beta and Ly alpha only the Ly alpha term acts
    lam = 1100.0 * (1 + z)
    t = photometry.madau95_transmission(np.array([lam]), z)[0]
    assert abs(t - np.exp(-3.6e-3 * (lam / 1215.67) ** 3.46)) < 1e-12
    # more absorption at higher redshift; below the Lyman limit the mean transmission (averaged
    # over sightlines, some of which miss every Lyman-limit system) keeps falling
    assert photometry.madau95_transmission(np.array([1100.0 * 6]), 5.0)[0] < t
    tl = photometry.madau95_transmission(np.array([900.0, 800.0, 700.0]) * 4, 3.0)
    assert np.all(np.diff(tl) < 0) and tl[-1] < 0.1


# ---------------------------------------------------------------------------
# Black hole masses
# ---------------------------------------------------------------------------

def test_virial_recipes():
    # VP06 H beta: log M = 6.91 + 2 log(FWHM/1000) + 0.5 log(L5100/1e44)
    assert abs(masses.virial_mass(2000.0, 44.0, "HB_VP06") - (6.91 + 2 * np.log10(2))) < 1e-9
    # VP06 C IV: zero point 6.66, slope 0.53
    assert abs(masses.virial_mass(1000.0, 45.0, "CIV_VP06") - (6.66 + 0.53)) < 1e-9
    # J0313-1806: Wang et al.'s FWHM and L3000 with the VO09 recipe give their 1.6e9
    w = lit.WANG21_J0313
    m = masses.virial_mass(w["fwhm_mgii"], np.log10(w["lambda_l3000"]), "MGII_VO09")
    assert abs(10 ** m / w["mbh"] - 1) < 0.01


def test_coatman_correction():
    # zero blueshift: FWHM / 0.62, zero point 6.71
    m = masses.coatman_corrected_mass(1000.0, 0.0, 44.0)
    assert abs(m - (6.71 + 2 * np.log10(1 / 0.62))) < 1e-9
    # a blueshift of 3000 km/s divides the width by 1.85: the mass drops by 0.53 dex vs 0 km/s
    d = masses.coatman_corrected_mass(1000.0, 3000.0, 44.0) - m
    assert abs(d - 2 * np.log10(0.62 / 1.85)) < 1e-9


def test_civ_recalibration_recovers_a_planted_relation():
    """Plant FWHM(C IV)/FWHM(Mg II) = 0.3 bs/1000 + 0.7 and C IV-free Mg II masses."""
    rng = np.random.default_rng(1)
    n = 20000
    bs = rng.uniform(-500, 3500, n)
    fwhm_mg = 10 ** rng.normal(3.5, 0.15, n)
    l1350 = rng.normal(45.5, 0.4, n)
    fwhm_c4 = fwhm_mg * (0.3 * bs / 1e3 + 0.7) * 10 ** rng.normal(0, 0.03, n)
    m_mg = 0.8 + 0.53 * (l1350 - 44) + 2 * np.log10(fwhm_mg)
    cal = masses.fit_civ_recalibration(bs, fwhm_c4, fwhm_mg, l1350, m_mg)
    assert abs(cal["slope"] - 0.3) < 0.01 and abs(cal["intercept"] - 0.7) < 0.01
    assert abs(cal["a"] - 0.8) < 0.01
    d = masses.recalibrated_civ_mass(bs, fwhm_c4, l1350, cal) - m_mg
    assert abs(np.median(d)) < 0.01 and np.polyfit(bs / 1e3, d, 1)[0] < 0.01


def test_eddington_ratio_definition():
    led = np.log10(physics.eddington_per_msun())
    assert abs(masses.eddington_ratio(led + 9.0, 9.0)) < 1e-12


# ---------------------------------------------------------------------------
# The luminosity function estimator
# ---------------------------------------------------------------------------

def test_page_carrera_recovers_a_planted_luminosity_function():
    """Quasars drawn from a known Phi(M) (a power law), uniform in comoving volume,
    observed through flux limits and a 50% random selection, in 20 realisations:
    the weighted 1/V_a estimator must return the input in bins wholly inside the
    limits. In bins only partly inside, Page & Carrera's method assumes Phi is flat
    across the bin; with Phi rising towards the faint, unseen part it reads a few
    per cent low (7% at 10% coverage here), which is why such bins are flagged."""
    k = qlf.KCorrection([0, 6], [0, 0])
    survey = qlf.Survey(area_deg2=500.0, kcorr=k, z_split=3.0, i_faint_lo=19.1, i_faint_hi=19.1,
                        i_bright=15.0, cosmology="ref")
    zlo, zhi, mlo, mhi = 1.0, 2.0, -29.0, -24.0
    phi0, slope = 1e-6, 0.4  # Phi(M) = phi0 * 10^(slope (M + 26)) per Mpc^3 per mag
    zz = np.linspace(zlo, zhi, 20001)
    dv = physics.dvdz_per_deg2(zz, "ref") * survey.area_deg2
    mg = np.linspace(mlo, mhi, 50001)
    phi = phi0 * 10 ** (slope * (mg + 26))
    cz = np.concatenate([[0], np.cumsum(0.5 * (dv[1:] + dv[:-1]) * np.diff(zz))])
    cm = np.concatenate([[0], np.cumsum(0.5 * (phi[1:] + phi[:-1]) * np.diff(mg))])
    medges = np.arange(-28.5, -24.49, 0.5)
    ratios, fills = [], None
    for seed in range(20):
        rng = np.random.default_rng(seed)
        ntot = rng.poisson(cz[-1] * cm[-1])
        z = np.interp(rng.random(ntot), cz / cz[-1], zz)
        m = np.interp(rng.random(ntot), cm / cm[-1], mg)
        i = m + physics.distance_modulus(z, "ref")
        keep = (i < 19.1) & (i > 15.0) & (rng.random(ntot) < 0.5)
        q = qlf.binned_qlf(z[keep], m[keep], np.full(keep.sum(), 2.0), survey, np.array([zlo, zhi]), medges)
        truth = np.array([np.mean(phi0 * 10 ** (slope * (np.linspace(a, b, 1001) + 26)))
                          for a, b in zip(q["m_lo"], q["m_hi"])])
        ratios.append((q["phi"] / truth).to_numpy())
        fills = q["fill"].to_numpy()
    r = np.array(ratios)
    mean, err = r.mean(axis=0), r.std(axis=0) / np.sqrt(len(r))
    full = fills > 0.999
    assert full.sum() >= 3
    assert np.all(np.abs(mean[full] - 1) < 3.5 * err[full] + 0.005), mean.round(3).tolist()
    assert np.all(mean[~full] < 1.02) and np.all(mean[~full] > 0.88), mean.round(3).tolist()


def test_accessible_volume_of_a_complete_bin():
    k = qlf.KCorrection([0, 6], [0, 0])
    survey = qlf.Survey(area_deg2=100.0, kcorr=k, z_split=3.0, cosmology="ref")
    va, fill = survey.accessible_volume(1.0, 1.2, -29.0, -28.5)
    zz = np.linspace(1.0, 1.2, 4001)
    exact = np.trapezoid(physics.dvdz_per_deg2(zz, "ref") * 100.0, zz) * 0.5
    assert fill == 1.0
    assert abs(va / exact - 1) < 1e-4


# ---------------------------------------------------------------------------
# The composite spectrum
# ---------------------------------------------------------------------------

def test_rebinning_conserves_flux():
    rng = np.random.default_rng(0)
    edges_in = np.cumsum(rng.uniform(0.5, 1.5, 3000)) + 1000.0
    vals = rng.uniform(0, 2, len(edges_in) - 1)
    edges_out = np.arange(1100.0, 3000.0, 1.0)
    out = composite.rebin_flux(edges_in, vals, edges_out)
    total_in = np.sum(vals * np.diff(edges_in)[: len(vals)]
                      * ((edges_in[:-1] >= 1100) & (edges_in[1:] <= 2999)))
    total_out = np.nansum(out * np.diff(edges_out))
    # the two windows differ only by the partial pixels at the ends
    assert abs(total_out - total_in) < 4.0


def test_stack_recovers_a_planted_template():
    """Spectra cut from one template at many redshifts, randomly scaled and noisy:
    the chained median composite must return the template's shape."""
    rng = np.random.default_rng(5)
    w = COMPOSITE_GRID
    template = w ** (-1.5) * (1 + 3 * np.exp(-0.5 * ((w - 1549) / 15) ** 2)
                              + 2 * np.exp(-0.5 * ((w - 2798) / 20) ** 2))
    n = 600
    zs = np.exp(rng.uniform(np.log1p(0.1), np.log1p(4.0), n)) - 1
    flux = np.full((n, len(w)), np.nan, np.float32)
    for k, z in enumerate(zs):
        cov = (w > 3650 / (1 + z)) & (w < 10200 / (1 + z))
        flux[k, cov] = (template[cov] * rng.uniform(0.2, 5.0) * (1 + rng.normal(0, 0.03, cov.sum())))
    res = composite.median_composite(flux, w, z=zs, n_iter=4)
    ok = res["n"] >= 30
    ratio = res["flux"][ok] / template[ok]
    assert np.std(ratio / np.median(ratio)) < 0.01
    # the power-law slope comes back: f_lambda ~ lambda^-1.5 is alpha_nu = -0.5
    a, _ = composite.power_law_fit(w, w ** (-1.5), composite.VB01_UV_WINDOWS)
    assert abs(a + 0.5) < 1e-3   # window means, not point values: exact to ~1e-5


def test_ccm89_at_v_band():
    # A_V / A_V = 1 at 5495 A (x = 1.82) for any R_V
    assert abs(composite.ccm89_alambda_over_av(np.array([1e4 / 1.82]))[0] - 1.0) < 1e-6


# ---------------------------------------------------------------------------
# The streaming FITS reader
# ---------------------------------------------------------------------------

def test_fitsstream_matches_astropy():
    from astropy.io import fits
    from astropy.table import Table
    rng = np.random.default_rng(1)
    n = 3000
    t = Table({"NAME": np.array([f"J{i:05d}" for i in range(n)]), "Z": rng.random(n),
               "ARR": rng.random((n, 6)), "PLATE": rng.integers(0, 12000, n).astype(np.int64),
               "U16": rng.integers(0, 65535, n).astype(np.uint16)})
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "t.fits"
        fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU(t)]).writeto(p)
        with open(p, "rb") as fi, gzip.open(str(p) + ".gz", "wb") as fo:
            shutil.copyfileobj(fi, fo)
        out = fitsstream.read_columns(str(p) + ".gz", {1: ["Z", "ARR", "PLATE", "U16", "NAME"]},
                                      chunk_rows=700)
    for c in ("Z", "ARR", "PLATE", "U16"):
        assert np.array_equal(out[1][c], np.asarray(t[c])), c
    assert (out[1]["NAME"] == np.asarray(t["NAME"])).all()


# ---------------------------------------------------------------------------
# Reproductions of published numbers (need the cached data)
# ---------------------------------------------------------------------------

def _have(path) -> bool:
    return Path(path).exists()


def test_selection_function_at_grid_nodes():
    if not _have(LIT_DIR / "J_AJ_131_2766_table1.tsv"):
        return
    import fetch_data
    t1 = fetch_data.load_vizier_table("J/AJ/131/2766/table1")
    sf = qlf.SelectionFunction(t1)
    row = t1[(t1["imag"] == 18.0) & (np.isclose(t1["z"], 2.8))].iloc[0]
    assert abs(sf(np.array([18.0]), np.array([2.8]), "point")[0] - row["point"]) < 1e-9


def test_richards06_luminosity_function_rebuilt():
    if not _have(LIT_DIR / "J_AJ_131_2766_table5.tsv"):
        return
    import analysis
    q = analysis.qlf_analysis(analysis.kcorr_r06())
    r = q["r06_reproduction"]
    assert r["n_matched"] == 101
    assert abs(r["median_dlog"]) < 0.01 and r["frac_within_0p05"] > 0.9
    s = q["sk12_comparison_without_image_factor"]
    assert abs(s["median_dlog"]) < 0.005 and s["frac_within_0p05"] > 0.8


def test_wushen_masses_reproduced():
    if not _have(WUSHEN_TABLE):
        return
    import pandas as pd
    cols = ["Z_SYS", "IF_BOSS_SDSS", "HBETA_BR_2", "HBETA_BR_4", "MGII_BR_2", "MGII_BR_4", "CIV_1", "CIV_2",
            "CIV_4", "LOGL5100", "LOGL3000", "LOGL1350", "LOGMBH_HB", "LOGMBH_MGII", "LOGMBH_CIV",
            "HBETA_BR_ERR_4", "MGII_BR_ERR_4", "CIV_ERR_4", "LOGL5100_ERR", "LOGL3000_ERR", "LOGL1350_ERR",
            "ZSYS_LINES_4", "ZSYS_LINES_ERR_4", "Z_FIT"]
    ws = pd.read_parquet(WUSHEN_TABLE, columns=cols).sample(50000, random_state=0)
    r = masses.recompute(ws)
    for ours, theirs in (("m_hb", "LOGMBH_HB"), ("m_mgii", "LOGMBH_MGII"), ("m_civ", "LOGMBH_CIV")):
        m = np.isfinite(r[ours]) & (ws[theirs] > 0)
        assert m.sum() > 1000
        assert np.max(np.abs(r.loc[m, ours] - ws.loc[m, theirs])) < 1e-6


def test_dr7_uniform_sample_size():
    if not _have(DR7_TABLE):
        return
    import pandas as pd
    d7 = pd.read_parquet(DR7_TABLE)
    u = d7[(d7["Uni"] == 1) & (d7["z"] >= 0.3) & (d7["z"] <= 5.0)]
    assert len(u) == lit.SK12["n_quasars"]   # Shen & Kelly 2012's 57,959


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
