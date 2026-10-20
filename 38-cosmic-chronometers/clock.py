"""Inside a chronometer: reading the 4000 Å break as a clock.

Old stars are cool and their spectra are full of metal absorption lines blueward
of 4000 Å, so an old galaxy's spectrum steps up sharply at 4000 Å. The size of
that step, D4000_n (the ratio of the mean flux in 4000-4100 Å to 3850-3950 Å;
Balogh et al. 1999), grows steadily as a passive galaxy ages. Over the range
that matters here, stellar-population models make the growth linear,

    D4000_n = A(Z) * age + B,

with a slope A (per Gyr) that depends mostly on the stars' metallicity Z. The
cosmic-chronometer equation then becomes (Moresco et al. 2012)

    H(z) = -A(Z) / (1 + z) * dz / dD4000_n,

so the slope of the D4000-redshift relation of the oldest galaxies, observed,
gives H once A is taken from the models. This module re-reads Moresco et al.
(2012)'s two SDSS measurements (z = 0.18, 0.20) from the binned medians and
calibration printed in that paper, and repeats the binning on J3's SDSS
galaxies with simple passive cuts to show how much the selection matters.
"""

from __future__ import annotations

import numpy as np

import config as cfg
import fetch_data

#: D4000_n above which the "high D4000" calibration regime applies.
REGIME_SPLIT = 1.8


def a_of_z_relations(rows) -> dict[str, np.ndarray]:
    """Quadratic A(Z) through the three tabulated metallicities, for both D4000
    regimes, plus the curves through A - err and A + err (the spread from the
    star-formation-history assumption), as Moresco et al. (2012) did."""
    rows = np.asarray(rows)
    Z = rows[:, 0]
    out = {}
    for regime, (ia, ie) in {"high": (1, 2), "low": (3, 4)}.items():
        for tag, sign in (("", 0), ("_lo", -1), ("_hi", 1)):
            out[regime + tag] = np.polyfit(Z, rows[:, ia] + sign * rows[:, ie], 2)
    return out


def _a_value(relations, regime, Z, tag=""):
    return np.polyval(relations[regime + tag], Z)


def h_from_pair(row_i, row_j, relations) -> dict:
    """H(z) from two binned D4000 medians: rows are (z, D, sigma_D, Z, sigma_Z, ...).

    Statistical error from the medians' errors; systematic error from the
    spread of H over the minimum, median and maximum A allowed by the
    metallicity range (Z +- sigma_Z) and the star-formation-history spread.
    """
    z1, d1, s1, Z1, sZ1 = row_i[:5]
    z2, d2, s2, Z2, sZ2 = row_j[:5]
    zmid = 0.5 * (z1 + z2)
    dz, dd = z2 - z1, d2 - d1
    Z, sZ = 0.5 * (Z1 + Z2), 0.5 * (sZ1 + sZ2)
    regime = "high" if 0.5 * (d1 + d2) > REGIME_SPLIT else "low"
    a_med = _a_value(relations, regime, Z)
    a_min = _a_value(relations, regime, Z - sZ, "_lo")
    a_max = _a_value(relations, regime, Z + sZ, "_hi")
    to_h = cfg.HUBBLE_TIME_GYR * dz / abs(dd) / (1.0 + zmid)  # (km/s/Mpc) per (1/Gyr)
    h_vals = np.array([a_min, a_med, a_max]) * to_h
    sigma_dd = np.hypot(s1, s2)
    return {"z": zmid, "dz": dz, "dD": dd, "A": a_med, "Z": Z, "regime": regime,
            "H": float(h_vals[1]), "stat": float(h_vals[1] * sigma_dd / abs(dd)),
            "syst": float(0.5 * (h_vals.max() - h_vals.min())),
            "dt_gyr": abs(dd) / a_med}


def combine(values) -> dict:
    """Inverse-variance mean of independent H estimates (statistical errors as
    weights, as in Moresco et al. 2012); the systematic part is averaged."""
    H = np.array([v["H"] for v in values])
    s = np.array([v["stat"] for v in values])
    w = 1.0 / s**2
    return {"z": float(np.average([v["z"] for v in values], weights=w)),
            "H": float(np.average(H, weights=w)), "stat": float(1.0 / np.sqrt(w.sum())),
            "syst": float(np.average([v["syst"] for v in values], weights=w))}


def reread_moresco2012(n_step: int = 2) -> dict:
    """Moresco et al.'s SDSS H(z) from their own D4000 medians (their Table 2,
    first four bins of each mass range, z < 0.3) and their A(Z) (Table 3).
    Pairs bin i with bin i + 2, as they did, for both SPS models."""
    m12 = fetch_data.load_m12_d4000()
    out = {"pairs": {}, "combined": {}}
    for model, rows in m12["A"].items():
        rel = a_of_z_relations(rows)
        per_mass = {}
        for mass, table in m12["bins"].items():
            sdss = [r for r in table if r[0] < 0.3]
            per_mass[mass] = [h_from_pair(sdss[i], sdss[i + n_step], rel)
                              for i in range(len(sdss) - n_step)]
        out["pairs"][model] = per_mass
        n_pairs = min(len(v) for v in per_mass.values())
        out["combined"][model] = [combine([per_mass[m][k] for m in per_mass])
                                  for k in range(n_pairs)]
    out["medians"] = {mass: [r for r in t if r[0] < 0.3] for mass, t in m12["bins"].items()}
    out["A_rows"] = m12["A"]
    return out


def sdss_simple_cut_medians(gal=None, *, z_edges=cfg.SDSS_Z_EDGES,
                            mass_bins=cfg.SDSS_MASS_BINS, cuts=cfg.SDSS_SIMPLE_CUTS) -> dict:
    """Median D4000_n of J3's SDSS galaxies in Moresco et al.'s bins, selected
    with simple passive cuts (no SED-template fits, no 2MASS, no metallicities).

    Returns per mass bin: rows (z_median, D, sigma_D, Z, sigma_Z, log M, N) in the
    same layout as the published table, with Z/Zsun = 1.1 +- 0.1 assumed (the
    value Moresco et al. adopted where metallicities were not measured).
    """
    gal = fetch_data.load_j3_galaxies() if gal is None else gal
    if gal is None:
        return {}
    d = gal["d4000_n"]
    passive = ((gal["reliable"] == 1) & np.isfinite(d) & (d > 0)
               & (gal["h_alpha_eqw"] > cuts["ha_eqw_min"])
               & (gal["dered_u"] - gal["dered_r"] > cuts["u_minus_r_min"])
               & (gal["specsfr_tot_p50"] < cuts["log_ssfr_max"]))
    out = {}
    for lo, hi in mass_bins:
        rows = []
        in_mass = passive & (gal["lgm_tot_p50"] > lo) & (gal["lgm_tot_p50"] < hi)
        for z0, z1 in zip(z_edges[:-1], z_edges[1:]):
            b = in_mass & (gal["z"] >= z0) & (gal["z"] < z1)
            n = int(b.sum())
            med = float(np.median(d[b]))
            mad = 1.482 * float(np.median(np.abs(d[b] - med)))
            rows.append((float(np.median(gal["z"][b])), med, mad / np.sqrt(n), 1.1, 0.1,
                         float(np.median(gal["lgm_tot_p50"][b])), n))
        out[f"{lo:.2f}-{hi:.2f}"] = rows
    return out


def h_from_rows(rows_by_mass, model_rows, n_step=2) -> dict:
    """The same pairing and calibration applied to any binned medians."""
    rel = a_of_z_relations(model_rows)
    per_mass = {m: [h_from_pair(r[i], r[i + n_step], rel) for i in range(len(r) - n_step)]
                for m, r in rows_by_mass.items()}
    return per_mass
