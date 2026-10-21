"""Self-checks: the formulas and fitters against independent calculations.

    python test_physics.py

Simulated light curves and mock samples appear here only as tests (fit a
known input, get it back); no figure or result uses them.
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

import config
import cosmology as cosmo
import standardize as st

PASSED, FAILED = [], []


def check(name: str, ok: bool, detail: str = "") -> None:
    (PASSED if ok else FAILED).append(name)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))


def test_distances():
    from astropy.cosmology import FlatLambdaCDM

    z = np.array([0.001, 0.01, 0.05, 0.1, 0.3, 1.0, 2.2])
    ours = cosmo.distance_modulus(z, z, 0.332, 73.5)
    ref = FlatLambdaCDM(H0=73.5, Om0=0.332, Tcmb0=0).distmod(z).value
    err = np.max(np.abs(ours - ref))
    check("distance modulus matches astropy (flat LCDM, no radiation)", err < 1e-6, f"max {err:.1e} mag")
    # The heliocentric factor: d_L = (1 + z_hel) D_C(z_cmb).
    a = cosmo.distance_modulus([0.05], [0.051], 0.3, 70)[0]
    b = cosmo.distance_modulus([0.05], [0.05], 0.3, 70)[0]
    check("one (1 + z_hel) factor in d_L", abs((a - b) - 5 * np.log10(1.051 / 1.05)) < 1e-12)
    # Low-z limit: d_L -> (c z / H0) (1 + (1 - q0) z / 2), q0 = 3 Om / 2 - 1.
    z0, q0 = 1e-3, 1.5 * 0.3 - 1
    low = cosmo.distance_modulus([z0], [z0], 0.3, 70)[0]
    expect = 5 * np.log10(cosmo.C_KM_S * z0 / 70 * (1 + (1 - q0) * z0 / 2)) + 25
    check("Hubble law plus deceleration term at low z", abs(low - expect) < 1e-6,
          f"{low - expect:+.1e} mag")


def test_cmb_frame():
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    apex = SkyCoord(l=cosmo.CMB_DIPOLE_L * u.deg, b=cosmo.CMB_DIPOLE_B * u.deg, frame="galactic").icrs
    anti = SkyCoord(l=(cosmo.CMB_DIPOLE_L + 180) * u.deg, b=-cosmo.CMB_DIPOLE_B * u.deg,
                    frame="galactic").icrs
    beta = cosmo.CMB_DIPOLE_V / cosmo.C_KM_S
    z_apex = cosmo.z_hel_to_cmb([0.05], [apex.ra.deg], [apex.dec.deg])[0]
    z_anti = cosmo.z_hel_to_cmb([0.05], [anti.ra.deg], [anti.dec.deg])[0]
    check("towards the dipole apex z_cmb = (1 + z_hel)/(1 - v/c) - 1",
          abs(z_apex - (1.05 / (1 - beta) - 1)) < 1e-9, f"{z_apex:.6f}")
    check("away from the apex z_cmb < z_hel", z_anti < 0.05 and abs(z_anti - (1.05 / (1 + beta) - 1)) < 1e-9)
    # Pantheon+ lists both frames: our conversion should reproduce its zCMB.
    p = pd.read_csv(config.PANTHEON_DIR / "Pantheon+SH0ES.dat", sep=r"\s+")
    zc = cosmo.z_hel_to_cmb(p["zHEL"], p["RA"], p["DEC"])
    d = np.abs(zc - p["zCMB"])
    check("reproduces Pantheon+ zCMB from zHEL", np.median(d) < 2e-5 and np.percentile(d, 99) < 2e-4,
          f"median |dz| {np.median(d):.1e}, 99% {np.percentile(d, 99):.1e}")


def test_velocity_error():
    z = 0.01
    exact = cosmo.sigma_mu_from_velocity(z, 300.0)
    approx = 5 / np.log(10) * 300 / (cosmo.C_KM_S * z)
    check("velocity error -> 5/ln10 * v/(cz) at low z", abs(exact / approx - 1) < 0.01,
          f"{exact:.4f} vs {approx:.4f} mag")
    check("300 km/s at z = 0.02 is ~0.1 mag (the spec's rule of thumb)",
          0.09 < cosmo.sigma_mu_from_velocity(0.02, 300.0) < 0.12)


def test_d1_recipe():
    p = pd.read_csv(config.PANTHEON_DIR / "Pantheon+SH0ES.dat", sep=r"\s+")
    r = st.recover_d1_recipe(p)
    rebuilt = p["mB"] + r.alpha * p["x1"] - r.beta * p["c"] - p["biasCor_m_b"] + r.const
    err = np.max(np.abs(rebuilt - p["m_b_corr"]))
    check("D1's m_b_corr rebuilt from mB, x1, c and biasCor", err < 2e-4, f"max {err:.1e} mag")
    check("alpha, beta match Brout et al. 2022 Table 2",
          abs(r.alpha - 0.148) < 0.004 and abs(r.beta - 3.09) < 0.04, f"{r.alpha:.4f}, {r.beta:.4f}")
    mb = -2.5 * np.log10(p["x0"]) + 10.635
    check("mB = -2.5 log10(x0) + 10.635 (SNANA convention)", np.max(np.abs(mb - p["mB"])) < 1e-3)


def test_tripp_error():
    rng = np.random.default_rng(1)
    cov = np.array([[0.02 ** 2, 0.002, 0.0003], [0.002, 0.15 ** 2, 0.001], [0.0003, 0.001, 0.03 ** 2]])
    draws = rng.multivariate_normal([18.0, 0.5, 0.05], cov, 200000)
    a, b = 0.15, 3.1
    mc = np.std(draws[:, 0] + a * draws[:, 1] - b * draws[:, 2])
    row = pd.DataFrame({"mB_err": [0.02], "x1_err": [0.15], "c_err": [0.03], "cov_mB_x1": [0.002],
                        "cov_mB_c": [0.0003], "cov_x1_c": [0.001]})
    ana = float(st.tripp_error(row, a, b).iloc[0])
    check("distance error propagation matches Monte Carlo", abs(ana / mc - 1) < 0.01, f"{ana:.4f} vs {mc:.4f}")


def test_fit_tripp_recovers_mock():
    rng = np.random.default_rng(2)
    n = 900
    x1 = rng.normal(0, 1, n)
    c = rng.normal(0.03, 0.09, n)
    mass = rng.normal(10, 0.7, n)
    step = np.where(mass > 10, 0.5, -0.5)
    mu = rng.uniform(34, 37, n)
    true = {"alpha": 0.16, "beta": 3.05, "M": -19.3, "gamma": 0.12}
    mB = mu + true["M"] - true["alpha"] * x1 + true["beta"] * c - true["gamma"] * step \
        + rng.normal(0, 0.12, n)
    e = pd.DataFrame({"mB": mB + rng.normal(0, 0.03, n), "x1": x1 + rng.normal(0, 0.1, n),
                      "c": c + rng.normal(0, 0.02, n), "mB_err": 0.03, "x1_err": 0.1,
                      "c_err": 0.02, "cov_mB_x1": 0.0, "cov_mB_c": 0.0, "cov_x1_c": 0.0})
    fit = st.fit_tripp(e, mu, np.zeros(n), host_mass=mass)
    ok = all(abs(fit.params[k] - true[k]) < 3.5 * fit.errors[k] for k in true)
    check("Tripp fit recovers injected alpha, beta, gamma, M (mock sample)", ok,
          ", ".join(f"{k} {fit.params[k]:.3f}+-{fit.errors[k]:.3f}" for k in true)
          + f", sigma_int {fit.sigma_int:.3f} (true ~0.12)")


def test_salt_fit_recovers_simulated_lightcurve():
    import sncosmo

    import saltfit

    model = saltfit.get_model()
    true = {"z": 0.05, "t0": 59000.0, "x0": 2e-3, "x1": 0.8, "c": 0.05, "mwebv": 0.04}
    model.set(**true)
    rng = np.random.default_rng(3)
    rows = []
    for band in ("ztfg", "ztfr", "ztfi"):
        times = np.arange(58960, 59060, 2.0 if band != "ztfi" else 5.0) + rng.uniform(0, 0.3)
        flux = model.bandflux(band, times, zp=30, zpsys="ab")
        err = np.sqrt(np.clip(flux, 0, None) * 2.0 + 300.0 ** 2)
        rows.append(pd.DataFrame({"mjd": times, "band": band, "flux": flux + rng.normal(0, err),
                                  "fluxerr": err, "zp": 30.0, "good": True, "flag": 0}))
    lc = pd.concat(rows, ignore_index=True)
    res = saltfit.fit_lightcurve(lc, true["z"], true["mwebv"], 59003.0)
    ok = res["ok"] and all(abs(res[k] - true[k]) < 4 * res[k + "_err"] for k in ("t0", "x1", "c"))
    check("SALT3 fit recovers a simulated light curve", ok,
          f"t0 {res.get('t0', np.nan) - true['t0']:+.2f} d, x1 {res.get('x1', np.nan):.3f}+-{res.get('x1_err', np.nan):.3f}, "
          f"c {res.get('c', np.nan):.3f}+-{res.get('c_err', np.nan):.3f}")
    check("peak magnitude from x0 uses the 10.635 convention",
          abs(res["mB"] - (-2.5 * np.log10(res["x0"]) + 10.635)) < 1e-12)


def test_d1_model_conventions():
    import saltfit

    info = saltfit.b21_info()
    check("D1's model: MAG_OFFSET and rest-frame range read from SALT2.INFO",
          info["MAG_OFFSET"] == 0.27 and info["RESTLAMBDA_RANGE"] == [2800.0, 7000.0], str(info))
    lc = pd.DataFrame({"mjd": [1.0, 2.0, 3.0], "band": ["ztfg", "ztfr", "ztfi"], "flux": 1.0,
                       "fluxerr": 0.1, "zp": 30.0, "good": True})
    low = set(saltfit._restframe_bands(lc, 0.05, "salt2-b21")["band"])
    high = set(saltfit._restframe_bands(lc, 0.20, "salt2-b21")["band"])
    check("SNANA's band rule: ZTF i dropped at z = 0.05, kept at z = 0.2",
          low == {"ztfg", "ztfr"} and high == {"ztfg", "ztfr", "ztfi"}, f"{sorted(low)}, {sorted(high)}")
    # The two models agree on the B-band peak of the average supernova to ~0.01 mag,
    # so the 0.27 mag is a convention, not a difference between the models.
    m3, m2 = saltfit.get_model("salt3"), saltfit.get_model("salt2-b21")
    for m in (m3, m2):
        m.set(z=0.0, t0=0.0, x0=1.0, x1=0.0, c=0.0, mwebv=0.0)
    d = m3.source.peakmag("bessellb", "vega") - m2.source.peakmag("bessellb", "vega")
    check("SALT3 and D1's SALT2 share their amplitude scale (B peak at x0 = 1)", abs(d) < 0.02,
          f"{d:+.4f} mag")


def test_alerce_conversion():
    from lightcurves import alerce_forced_to_lightcurve, subtract_baseline, detected_phases

    rows = [{"mjd": 1.0, "fid": 2, "mag": 20.0, "e_mag": 0.1, "isdiffpos": 1,
             "extra_fields": {"magzpsci": 26.0, "procstatus": "0"}},
            {"mjd": 2.0, "fid": 1, "mag": 21.0, "e_mag": 0.5, "isdiffpos": -1,
             "extra_fields": {"magzpsci": 26.0, "procstatus": "0"}},
            {"mjd": 3.0, "fid": 3, "mag": 19.0, "e_mag": 0.05, "isdiffpos": 1,
             "extra_fields": {"magzpsci": 26.0, "procstatus": "56"}}]
    lc = alerce_forced_to_lightcurve(rows)
    f = dict(zip(lc["mjd"], lc["flux"]))
    check("ALeRCE magnitudes -> signed fluxes on zp = 30",
          abs(f[1.0] - 10 ** (0.4 * 10)) < 1e-6 and f[2.0] < 0 and not lc.loc[lc["mjd"] == 3.0, "good"].iloc[0])
    base = pd.DataFrame({"mjd": np.arange(0, 60.0), "band": "ztfr", "flux": 50.0, "fluxerr": 10.0,
                         "zp": 30.0, "good": True})
    out, off = subtract_baseline(base, 50.0, 0.0)
    check("baseline subtraction removes a constant offset", abs(out["flux"]).max() < 1e-9 and "ztfr" in off)
    lc2 = pd.DataFrame({"mjd": [100.1, 100.2, 101.1, 102.1, 99.0], "band": ["ztfg", "ztfg", "ztfg", "ztfr", "ztfr"],
                        "flux": [100.0, 100, 100, 100, 100], "fluxerr": [10.0] * 5, "zp": 30.0, "good": True})
    cov = detected_phases(lc2, 101.0, 0.0)
    check("same-night same-band detections count once", cov["n_phases"] == 4 and cov["n_before"] == 2,
          str(cov))


if __name__ == "__main__":
    for fn in [v for k, v in list(globals().items()) if k.startswith("test_")]:
        fn()
    print(f"\n{len(PASSED)} passed, {len(FAILED)} failed")
    sys.exit(1 if FAILED else 0)
