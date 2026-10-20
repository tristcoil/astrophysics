"""The Cosmicflows-4 Hubble diagram and the Hubble constant.

Each CF4 group has a distance measured without its redshift. Dividing the
expansion velocity (CMB-frame cz times CF4's curvature factor f) by that
distance gives a Hubble parameter H_i for every group. Nearby, peculiar
velocities of a few hundred km/s are a large fraction of cz, so H_i scatters
wildly; beyond 4,000 km/s they are a small fraction, and the mean of log H_i
measures the Hubble constant (Tully et al. 2023, Sect. 12). Averaging the
*logarithm* matters: distance errors are symmetric in log d, not in d.

The distances of the five methods that reach far (TF, FP, SBF, SN Ia, SN II)
were tied together by CF4 on an arbitrary scale and then shifted onto the
scale of Cepheids, the tip of the red giant branch and the maser galaxy NGC
4258 by one number, the zero point. H0 is therefore exactly as good as that
zero point: the first rungs of the distance ladder (project B4), the same rungs
that set D1's supernova H0.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import catalogs
import config
from config import C_KMS


def groups(version: str = "corrected") -> pd.DataFrame:
    if version == "corrected":
        import cf4fix

        return cf4fix.corrected_groups()
    return catalogs.cf4_groups()


def hubble_parameters(g: pd.DataFrame | None = None, version: str = "corrected") -> pd.DataFrame:
    """Per-group Hubble parameter recomputed from the distance modulus."""
    g = groups(version) if g is None else g
    out = g[["1PGC", "V3k", "Vls", "fV3k", "DMzp", "e_DMzp", "Dist", "logHi", "method", "RAdeg", "DEdeg",
             "GLON", "GLAT"]].copy()
    out["d_mpc"] = 10 ** ((out["DMzp"] - 25.0) / 5.0)
    from cosmology import cf4_curvature_factor

    out["fcz"] = cf4_curvature_factor(out["V3k"] / C_KMS) * out["V3k"]
    with np.errstate(divide="ignore", invalid="ignore"):
        out["logH"] = np.log10(out["fcz"] / out["d_mpc"])
    return out


def _mean_log_h(logh: np.ndarray) -> dict:
    n = len(logh)
    m = float(np.mean(logh))
    sd = float(np.std(logh, ddof=1))
    se = sd / np.sqrt(n)
    H0 = 10**m
    return {"n": n, "mean_logH": m, "rms_logH": sd, "H0": H0, "H0_stat": H0 * np.log(10) * se}


def hubble_constant(version: str = "corrected") -> dict:
    """H0 the way CF4 measured it, plus robustness checks and the zero point."""
    h = hubble_parameters(version=version)
    far = h[(h["V3k"] > config.HUBBLE_VMIN) & np.isfinite(h["logH"])]
    res = {"version": version, "baseline": _mean_log_h(far["logH"].to_numpy())}
    if version == "published":
        res["baseline"]["check_vs_table_logHi_max_abs"] = float(np.max(np.abs(far["logH"] - far["logHi"])))
    # The zero point: one magnitude offset between the arbitrary five-method
    # scale (DMav) and the calibrated one (DMzp). An error of dmu in it scales
    # every distance by 10^(dmu/5), hence H0 by 10^(-dmu/5).
    g = catalogs.cf4_groups()  # the zero point is read from the published table
    off = (g["DMzp"] - g["DMav"])[g["DMav"] > 0]
    zp = float(np.median(off))
    res["zero_point"] = {
        "offset_mag": zp, "n_groups_with_offset": int(np.isclose(off, zp, atol=1e-3).sum()),
        "paper_eq8_offset": -0.138, "paper_eq8_err": 0.022,
        "H0_if_paper_offset": res["baseline"]["H0"] * 10 ** ((zp - (-0.138)) / 5.0),
        "H0_zp_err": res["baseline"]["H0"] * np.log(10) / 5.0 * 0.022,
        "shift_to_planck_mag": float(5 * np.log10(res["baseline"]["H0"] / config.PLANCK18.H0)),
    }
    # Robustness: how CF4's figure averages, weighting, velocity ranges, methods, sky sectors.
    rob = {}
    edges = np.arange(config.HUBBLE_VMIN, far["V3k"].max() + config.HUBBLE_BIN_KMS, config.HUBBLE_BIN_KMS)
    idx = np.digitize(far["V3k"], edges)
    binmeans = far.groupby(idx)["logH"].mean()
    rob["mean_of_bin_means"] = float(10 ** binmeans.mean())
    w = 1.0 / far["e_DMzp"] ** 2
    rob["inverse_variance_weighted"] = float(10 ** np.average(far["logH"], weights=w))
    rob["median"] = float(10 ** np.median(far["logH"]))
    rob["velocity_ranges"] = {}
    for lo, hi in [(4000, 10000), (10000, 16000), (16000, 33000)]:
        s = far[(far["V3k"] > lo) & (far["V3k"] <= hi)]
        rob["velocity_ranges"][f"{lo}-{hi}"] = _mean_log_h(s["logH"].to_numpy())
    rob["methods"] = {m: _mean_log_h(s["logH"].to_numpy()) for m, s in far.groupby("method") if len(s) >= 5}
    north_eq = far["DEdeg"] > 0
    north_gal = far["GLAT"] > 0
    rob["sectors"] = {}
    for ne in (True, False):
        for ng in (True, False):
            s = far[(north_eq == ne) & (north_gal == ng)]
            key = f"{'N' if ne else 'S'} celestial / {'N' if ng else 'S'} galactic"
            rob["sectors"][key] = _mean_log_h(s["logH"].to_numpy())
    res["robustness"] = rob
    # The curve of log H against velocity (CF4's Fig. 18): bin means and errors.
    rob["sdss_vs_not"] = {
        "with_sdss_members": _mean_log_h(far.loc[far["1PGC"].isin(_sdss_groups(version)), "logH"].to_numpy()),
        "without": _mean_log_h(far.loc[~far["1PGC"].isin(_sdss_groups(version)), "logH"].to_numpy()),
    }
    allh = h[(h["V3k"] > 0) & np.isfinite(h["logH"])]
    edges_all = np.arange(0, 33000, config.HUBBLE_BIN_KMS)
    idx = np.digitize(allh["V3k"], edges_all)
    grp = allh.groupby(idx)["logH"]
    res["binned"] = {
        "v_mid": [float(edges_all[i - 1] + config.HUBBLE_BIN_KMS / 2) for i in grp.mean().index],
        "mean_logH": grp.mean().tolist(), "sem_logH": (grp.std() / np.sqrt(grp.count())).tolist(),
        "n": grp.count().tolist(),
    }
    res["published"] = config.PUBLISHED["cf4_H0"]
    return res


def _sdss_groups(version: str) -> set:
    import cf4fix

    g = cf4fix.corrected_groups()
    return set(g.loc[g["sdss_corrected"], "1PGC"])
