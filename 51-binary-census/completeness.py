"""What fraction of companions could Gaia have found? The detection windows.

For a sample of primaries (real stars from the 100 pc census), we create mock
companions in narrow bins of orbital period (``mocks.make_mocks``) and ask, for
each, which of Gaia's channels would have caught it:

* astrometric orbit / acceleration / RUWE > 1.4 only: gaiamock runs DR3's
  astrometric cascade on the simulated epoch astrometry (``astrometric.py``);
* spectroscopic orbit: our radial-velocity model (``spectro.py``), using each
  star's own number of RVS transits and per-transit precision;
* resolved wide pair: the projected separation at a random moment of a
  random orientation, against the contrast curve measured in ``wide.py``.

The completeness of a channel in a period bin is the fraction of mock
companions it would catch, averaged over the assumed distributions of mass
ratio and eccentricity (Moe & Di Stefano 2017; varied as a systematic).
"""

from __future__ import annotations

import numpy as np

import astrometric
import binaries as B
import mocks
import stellar
import wide

#: log10(P/day) bins for the census: from 0.3 d to ~3e8 d (beyond 1 pc for Sun-like pairs).
LOGP_EDGES = np.arange(-0.5, 9.01, 0.25)
#: gaiamock is only needed where the light centre's wobble could matter.
LOGP_ASTROMETRY = (-0.5, 6.0)


def _f(col):
    return np.ma.filled(np.ma.asarray(col, dtype=float), np.nan)


def primary_arrays(t) -> dict:
    """The per-star quantities the detection models need."""
    return {
        "ra": _f(t["ra"]), "dec": _f(t["dec"]), "plx": _f(t["parallax"]), "pmra": _f(t["pmra"]),
        "pmdec": _f(t["pmdec"]), "g": _f(t["phot_g_mean_mag"]), "m1": np.asarray(t["mass_colour"], dtype=float),
        "dist": np.asarray(t["dist_pc"], dtype=float), "grvs": _f(t["grvs_mag"]),
        "n_rv": _f(t["rv_nb_transits"]), "rv_err": _f(t["radial_velocity_error"]),
        "teff_rv": _f(t["rv_template_teff"]), "bp_rp": np.asarray(t["bp_rp"], dtype=float),
        "source_id": np.asarray(t["source_id"], dtype=np.int64),
    }


def astrometric_outcomes(mock: mocks.MockSet, prim: dict, fratio_g, *, seed: int, workers: int = 6) -> np.ndarray:
    """gaiamock outcome for each mock (see astrometric.run_one); systems outside
    LOGP_ASTROMETRY are marked 1 (no astrometric signal)."""
    rng = np.random.default_rng(seed)
    out = np.zeros((len(mock), 4))
    out[:, 0] = 1
    out[:, 1:] = np.nan
    use = (mock.logp >= LOGP_ASTROMETRY[0]) & (mock.logp < LOGP_ASTROMETRY[1])
    idx = np.where(use)[0]
    tasks = []
    for i in idx:
        k = mock.star_index[i]
        g_sys = prim["g"][k] - 2.5 * np.log10(1 + fratio_g[i])
        tasks.append((int(rng.integers(0, 2**31 - 1)), float(prim["ra"][k]), float(prim["dec"][k]),
                      float(prim["plx"][k]), float(prim["pmra"][k]), float(prim["pmdec"][k]),
                      float(mock.m1[i]), float(mock.m2[i]), float(mock.period[i]), float(mock.t_peri[i]),
                      float(mock.ecc[i]), float(mock.node[i]), float(np.degrees(mock.inc[i])), float(mock.omega[i]),
                      float(g_sys), float(fratio_g[i])))
    if tasks:
        out[idx] = astrometric.run_many(tasks, workers=workers)
    return out


def resolved_probability(mock: mocks.MockSet, prim: dict, curve: dict, rng) -> np.ndarray:
    """Probability that El-Badry et al.'s catalogue contains the companion as a
    resolved pair: projected separation at a random moment, Gaia's contrast
    curve, and a parallax for the companion (G < 20.3)."""
    k = mock.star_index
    a = B.semi_major_axis_au(mock.period, mock.m1 + mock.m2)
    s = a * B.projected_separation_factor(rng, mock.ecc)
    theta = s * prim["plx"][k] / 1000.0  # arcsec
    mg1 = stellar.mg_from_mass(mock.m1)
    mg2 = stellar.mg_from_mass(mock.m2)
    dg = mg2 - mg1
    p = wide.detection_probability(theta, np.where(np.isfinite(dg), dg, 13.0), curve)
    g2 = prim["g"][k] + np.where(np.isfinite(dg), dg, 99.0)
    p = np.where(g2 < 20.3, p, 0.0)
    # El-Badry et al. searched out to 1 pc of projected separation
    return np.where(s < 206_265.0, p, 0.0)


def spectroscopic_outcomes(mock: mocks.MockSet, prim: dict, noise, grvs_coeffs, eff_table, rng) -> dict:
    """Spectroscopic verdict for every mock: the idealised DR3 selection
    (``spectro.simulate``) times the acceptance efficiency calibrated on SB9."""
    import spectro
    n = len(mock)
    outcome = np.full(n, "", dtype="U9")
    p_orbit = np.zeros(n)
    beta = stellar.flux_ratio_rvs(mock.m1, mock.m2, grvs_coeffs)
    k1 = B.k1_kms(mock.period, mock.m1, mock.m2, mock.ecc, np.sin(mock.inc))
    for i in range(n):
        k = mock.star_index[i]
        teff, grvs, nrv = prim["teff_rv"][k], prim["grvs"][k], prim["n_rv"][k]
        if not (spectro.TEFF_RANGE[0] < teff < spectro.TEFF_RANGE[1]):
            continue
        # the system's G_RVS includes the companion's light
        grvs_sys = grvs - 2.5 * np.log10(1 + beta[i]) if np.isfinite(grvs) else np.nan
        if not (spectro.GRVS_RANGE[0] <= grvs_sys <= spectro.GRVS_RANGE[1]) or not (nrv >= spectro.N_MIN):
            continue
        t = spectro.rvs_epochs(prim["ra"][k], prim["dec"][k], int(nrv), rng)
        res = spectro.simulate(t, float(noise(grvs_sys)), mock.period[i], mock.t_peri[i], mock.ecc[i],
                               mock.omega[i], k1[i], mock.q[i], beta[i], rng)
        outcome[i] = res
        if res.startswith("orbit"):
            p_orbit[i] = spectro.efficiency(mock.logp[i], mock.ecc[i], eff_table, sb2=res == "orbit_sb2")
    return {"outcome": outcome, "p_orbit": p_orbit, "beta_rvs": beta, "k1": k1}


def channel_completeness(mc: dict, spec: dict, *, efficiency: str = "sb9") -> dict:
    """Per log P bin: the fraction of mock companions each channel catches.

    ``efficiency``: "sb9" applies DR3's spectroscopic acceptance efficiency
    calibrated on SB9; "ideal" assumes DR3 published every orbit that passed
    its documented cuts (an upper bound on completeness).
    """
    mock, ast = mc["mock"], mc["ast"]
    edges = mc["edges"]
    idx = np.clip(np.digitize(mock.logp, edges) - 1, 0, len(edges) - 2)
    code = ast[:, 0]
    p_ast = (code == 12).astype(float)
    p_acc = astrometric.dr3_acceleration_published(code, ast[:, 2], ast[:, 3]).astype(float)
    # accelerations found in main processing but dropped by the final s > 20 cut
    # leave the star with RUWE > 1.4 and no solution
    acc_dropped = np.isin(code, (7, 9)) & (p_acc == 0)
    p_spec = spec["p_orbit"] if efficiency == "sb9" else np.char.startswith(spec["outcome"], "orbit").astype(float)
    p_orbit = 1 - (1 - p_spec) * (1 - p_ast)
    p_wide = mc["p_wide"]
    p_trend = (spec["outcome"] == "trend").astype(float)
    p_ruwe = ((code == -12) | (code == -99) | acc_dropped).astype(float) * (1 - p_orbit) * (1 - p_acc)
    p_any_period = 1 - (1 - p_orbit) * (1 - p_wide)
    out = {"edges": edges, "centres": 0.5 * (edges[1:] + edges[:-1])}
    for name, arr in (("spectroscopic", p_spec), ("astrometric", p_ast), ("orbit", p_orbit),
                      ("acceleration", p_acc), ("rv_trend", p_trend), ("wide", p_wide), ("ruwe_only", p_ruwe),
                      ("period_known", p_any_period)):
        out[name] = np.array([arr[idx == k].mean() if np.any(idx == k) else np.nan for k in range(len(edges) - 1)])
    out["per_mock"] = {"spec": p_spec, "ast": p_ast, "acc": p_acc, "wide": p_wide, "trend": p_trend,
                       "ruwe": p_ruwe, "orbit": p_orbit, "bin": idx}
    return out
