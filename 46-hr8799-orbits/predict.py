"""The prediction test: orbits fitted to data up to 2016 Sep 19 (Wang et al.'s
data set) predict where the planets will be, and are compared with positions
measured afterwards by other instruments (SPHERE 2017-2021, LBT 2020, GRAVITY
2018, JWST 2023, Keck 2025) - none of which the fit has seen.

Each comparison uses the measurement's errors, the spread of the predictions
(from posterior samples) and the instrument's calibration uncertainty (plate
scale and north angle, config.CALIBRATION), combined as 2x2 covariances.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import astrometry as A
import config as C
import kepler as K

#: everything measured after the training set ends, from instruments that
#: Wang et al. did not use
TEST_SOURCES = ("IRDIS-Z22", "IFS-Z22", "LUCI-Z22", "GRAVITY-L19", "NIRCam-B25", "AMI-N26", "NIRC2-N26")


def predicted_positions(model, samples: np.ndarray, epochs_mjd: np.ndarray) -> dict:
    """planet -> (ra, dec) arrays of shape (n_samples, n_epochs), in mas."""
    phys = model.physical(samples)
    t = np.asarray(epochs_mjd, float)[None, :]
    el = {p: tuple(np.asarray(v)[:, None] if np.ndim(v) else v for v in phys["elements"][p])
          for p in C.INNER_TO_OUTER}
    masses = ({p: np.asarray(phys["masses"][p])[:, None] for p in C.INNER_TO_OUTER}
              if model.fit_masses else None)
    sky = K.system_offsets(t, el, np.asarray(phys["plx"])[:, None], np.asarray(phys["mstar"])[:, None], masses)
    return sky


def test_rows() -> pd.DataFrame:
    df = A.load(TEST_SOURCES, planets=C.PLANETS)
    return df[df.mjd > C.PREDICTION_CUTOFF_MJD].reset_index(drop=True)


def compare(model, samples: np.ndarray, rows: pd.DataFrame | None = None, n: int = 4000,
            seed: int = C.RANDOM_SEED, jitter: dict | None = None) -> pd.DataFrame:
    """Predictions vs later measurements. ``jitter`` (source -> mas) adds an
    instrument's extra scatter in quadrature to its quoted errors - e.g. the
    jitters fitted to all the data (SPHERE's quoted 1-2 mas errors are too
    small: IFS and IRDIS disagree by up to ~10 mas on the same night)."""
    rows = test_rows() if rows is None else rows
    if len(samples) > n:
        samples = samples[np.random.default_rng(seed).choice(len(samples), n, replace=False)]
    epochs = rows.mjd.to_numpy()
    sky = predicted_positions(model, samples, epochs)
    out = []
    for k, r in rows.iterrows():
        ra_s, dec_s = sky[r.planet][0][:, k], sky[r.planet][1][:, k]
        pred = np.array([np.median(ra_s), np.median(dec_s)])
        cov_pred = np.cov(np.vstack([ra_s, dec_s]))
        obs = np.array([r.ra, r.dec])
        c = r["corr"] if np.isfinite(r["corr"]) else 0.0
        cov_obs = np.array([[r.ra_err ** 2, c * r.ra_err * r.dec_err], [c * r.ra_err * r.dec_err, r.dec_err ** 2]])
        if jitter and r.source in jitter:
            cov_obs = cov_obs + jitter[r.source] ** 2 * np.eye(2)
        # calibration: radial (scale) and tangential (north angle) components
        cal = C.CALIBRATION.get(r.source, (0.0, 0.0, ""))
        sep = np.hypot(*pred)
        u_r = pred / sep
        u_t = np.array([u_r[1], -u_r[0]])
        cov_cal = (cal[1] * sep) ** 2 * np.outer(u_r, u_r) + (np.radians(cal[0]) * sep) ** 2 * np.outer(u_t, u_t)
        cov = cov_obs + cov_pred + cov_cal
        res = obs - pred
        d2 = float(res @ np.linalg.solve(cov, res))
        out.append(dict(planet=r.planet, source=r.source, date=r.date, year=r.year,
                        obs_ra=r.ra, obs_dec=r.dec, pred_ra=pred[0], pred_dec=pred[1],
                        res_ra=res[0], res_dec=res[1], res_mas=float(np.hypot(*res)),
                        res_radial=float(res @ u_r), res_tangential=float(res @ u_t),
                        sigma_pred_mas=float(np.sqrt(np.trace(cov_pred) / 2)),
                        sigma_obs_mas=float(np.sqrt(np.trace(cov_obs) / 2)),
                        mahalanobis2=d2, p_value=float(np.exp(-0.5 * d2))))
    return pd.DataFrame(out)


def summary(table: pd.DataFrame) -> dict:
    """Headline numbers of a prediction test."""
    chi2 = table.mahalanobis2.sum()
    dof = 2 * len(table)
    from scipy.stats import chi2 as chi2dist
    return {"n": int(len(table)), "chi2": float(chi2), "dof": int(dof),
            "p_value": float(chi2dist.sf(chi2, dof)),
            "frac_within_2sigma": float(np.mean(table.p_value > np.exp(-0.5 * 6.18))),   # 2D 2-sigma (95.4%)
            "median_miss_mas": float(table.res_mas.median()),
            "worst": table.sort_values("mahalanobis2").iloc[-1][["planet", "source", "date", "res_mas"]].to_dict()}
