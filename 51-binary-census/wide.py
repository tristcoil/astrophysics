"""Resolved wide binaries: Gaia's sensitivity to close pairs, and mass ratios.

Gaia resolves two stars as separate sources only if they are far enough apart
on the sky for their brightness difference: near-equal pairs down to ~0.7",
a companion 4 mag fainter only beyond ~1" and 7 mag fainter beyond ~2"
(El-Badry et al. 2021, Sec. 4.1). We measure that sensitivity ourselves from
the catalogue: within a narrow range of brightness difference Delta G the true
number of pairs per log(angular separation) varies smoothly (a power law over
the 0.3-15" range), so the shortfall at small separations *is* the
incompleteness. Following El-Badry et al. (2019, Eq. D5) we model the
detection probability as

    S(theta) = 1 / (1 + (theta / theta_0)^-beta),

and fit theta_0 and beta in bins of Delta G.

Masses of both stars come from their absolute magnitudes on Mamajek's dwarf
sequence, giving mass ratios q = m2/m1 precise to ~2%: good enough to see
the excess of "twins" with q > 0.95 (El-Badry et al. 2019).
"""

from __future__ import annotations

import numpy as np
from scipy import optimize

import samples
import stellar

DG_BINS = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9.5, 11, 13])
LOGTHETA_FIT = (-0.6, 1.2)
LOGTHETA_STEP = 0.05


def sensitivity(theta, theta0, beta):
    theta = np.asarray(theta, dtype=float)
    return 1.0 / (1.0 + (np.maximum(theta, 1e-6) / theta0) ** (-beta))


def _fit_bin(logtheta: np.ndarray):
    edges = np.arange(LOGTHETA_FIT[0], LOGTHETA_FIT[1] + 1e-9, LOGTHETA_STEP)
    n, _ = np.histogram(logtheta, edges)
    x = 0.5 * (edges[1:] + edges[:-1])

    def model(p):
        log_a, gamma, log_t0, beta = p
        return np.exp(log_a + gamma * x * np.log(10)) * sensitivity(10 ** x, 10 ** log_t0, beta)

    def nll(p):
        if not (1.0 < p[3] < 40.0):
            return 1e12
        mu = model(p)
        return float(np.sum(mu - n * np.log(mu + 1e-300)))

    best = None
    for t0 in (-0.2, 0.0, 0.2, 0.4, 0.6):
        p0 = [np.log(max(n.max(), 1)), -0.3, t0, 8.0]
        r = optimize.minimize(nll, p0, method="Nelder-Mead",
                              options={"maxiter": 20000, "xatol": 1e-6, "fatol": 1e-6})
        if best is None or r.fun < best.fun:
            best = r
    # curvature errors on log theta0
    h = optimize.approx_fprime(best.x, nll, 1e-5)
    try:
        hess = np.array([optimize.approx_fprime(best.x, lambda p, i=i: optimize.approx_fprime(p, nll, 1e-5)[i],
                                                1e-5) for i in range(4)])
        cov = np.linalg.inv(0.5 * (hess + hess.T))
        err = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        err = np.full(4, np.nan)
    del h
    return best.x, err, (x, n, model(best.x))


def fit_contrast_curve(pairs=None) -> dict:
    """theta_0 and beta of the detection probability in bins of Delta G, plus a
    smooth interpolation used by the completeness model."""
    if pairs is None:
        pairs = samples.wide_pairs()
    lt = np.log10(np.asarray(pairs["theta_arcsec"]))
    dg = np.asarray(pairs["dg"])
    rows = []
    for lo, hi in zip(DG_BINS[:-1], DG_BINS[1:]):
        m = (dg >= lo) & (dg < hi) & np.isfinite(lt)
        p, e, curve = _fit_bin(lt[m])
        rows.append({"dg_lo": lo, "dg_hi": hi, "dg_mid": float(np.median(dg[m])), "n": int(m.sum()),
                     "log_theta0": p[2], "log_theta0_err": e[2], "beta": p[3], "beta_err": e[3],
                     "gamma": p[1], "curve": curve})
    dgm = np.array([r["dg_mid"] for r in rows])
    lt0 = np.array([r["log_theta0"] for r in rows])
    beta = np.array([r["beta"] for r in rows])
    return {"bins": rows, "dg_mid": dgm, "log_theta0": lt0, "beta": beta}


def detection_probability(theta_arcsec, delta_g, curve: dict):
    """Probability that Gaia resolves a pair at this separation and contrast."""
    dg = np.clip(np.asarray(delta_g, dtype=float), curve["dg_mid"][0], curve["dg_mid"][-1])
    t0 = 10 ** np.interp(dg, curve["dg_mid"], curve["log_theta0"])
    beta = np.interp(dg, curve["dg_mid"], curve["beta"])
    return sensitivity(theta_arcsec, t0, beta)


# ---------------------------------------------------------------------------
# Mass ratios of resolved main-sequence pairs
# ---------------------------------------------------------------------------

def ms_pairs(max_dist_pc: float = 200.0, min_poe: float = 20.0):
    """Main-sequence + main-sequence pairs with photometric masses for both stars.

    Both stars must lie in the main-sequence band about Mamajek's sequence (so
    neither is a giant, subgiant or white dwarf, nor an obviously unresolved
    binary far above the sequence), and have parallax/error > ``min_poe``.
    """
    w = samples.wide_pairs()
    d = 1000.0 / np.asarray(w["plx"])
    ridge1 = stellar.ms_ridge_mg(np.asarray(w["bprp1"]))
    ridge2 = stellar.ms_ridge_mg(np.asarray(w["bprp2"]))
    dg1 = np.asarray(w["abs_g1"]) - ridge1
    dg2 = np.asarray(w["abs_g2"]) - ridge2
    # El-Badry 2019 used parallax_over_error > 20 for both stars; our VizieR
    # extraction has parallax errors for both.
    m1 = stellar.mass_from_mg(np.asarray(w["abs_g1"]))
    m2 = stellar.mass_from_mg(np.asarray(w["abs_g2"]))
    keep = ((d < max_dist_pc) & (np.abs(dg1) < 0.75) & (np.abs(dg2) < 0.75)
            & np.isfinite(m1) & np.isfinite(m2) & (w["bin_type"] == "MSMS"))
    out = w[keep]
    out["m1"] = m1[keep]
    out["m2"] = m2[keep]
    out["q"] = out["m2"] / out["m1"]
    out["dist_pc"] = d[keep]
    return out


# ---------------------------------------------------------------------------
# The twin excess: El-Badry et al. (2019)'s model, fitted with selection
# ---------------------------------------------------------------------------

#: El-Badry et al. 2019's bins (their Table G1).
TWIN_SEP_BINS = [(50, 350), (350, 600), (600, 1000), (1000, 2500), (2500, 5000), (5000, 15000),
                 (15000, 50000)]
TWIN_MASS_BINS = [(0.1, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.2), (1.2, 2.5)]
Q_GRID = np.linspace(0.3, 1.0, 141)
G_FAINT = 20.3  # companions fainter than this rarely have parallax/error > 5 within 200 pc


def colour_subset(pairs):
    """Pairs in which both stars have BP-RP colours (needed to call them dwarfs)."""
    return pairs[np.isfinite(np.asarray(pairs["bprp1"])) & np.isfinite(np.asarray(pairs["bprp2"]))]


def selection_matrix(p, curve, q_grid=Q_GRID):
    """S[i, j]: probability that pair i, had its mass ratio been q_grid[j], would
    be resolved (contrast curve) with a companion bright enough for a parallax."""
    m1 = np.asarray(p["m1"])[:, None]
    m2 = m1 * q_grid[None, :]
    mg1 = stellar.mg_from_mass(m1)
    mg2 = stellar.mg_from_mass(m2)
    dg = mg2 - mg1
    theta = np.asarray(p["theta_arcsec"])[:, None] * np.ones_like(dg)
    s = detection_probability(theta, np.where(np.isfinite(dg), dg, 13.0), curve)
    g2 = mg2 + 5 * np.log10(np.asarray(p["dist_pc"])[:, None] / 10.0)
    return np.where(np.isfinite(g2) & (g2 < G_FAINT), s, 0.0)


def twin_model_pdf(q, gamma_lo, gamma_hi, f_twin, q_twin, q_break=0.5):
    """El-Badry et al. (2019): broken power law in q (break at q_break) plus a
    uniform excess on [q_twin, 1] holding a fraction F_twin of all q > 0.3
    binaries. Normalised on [0.3, 1]."""
    q = np.asarray(q, dtype=float)
    base = np.where(q < q_break, (q / q_break) ** gamma_lo, (q / q_break) ** gamma_hi)
    grid = Q_GRID
    gb = np.where(grid < q_break, (grid / q_break) ** gamma_lo, (grid / q_break) ** gamma_hi)
    norm = np.trapezoid(gb, grid)
    pdf = (1 - f_twin) * base / norm + f_twin * ((q >= q_twin) / (1.0 - q_twin))
    return pdf


def fit_twins(p, s_matrix, *, q_break=0.5, n_walkers=24, n_steps=900, seed=1):
    """Posterior of (gamma_lo, gamma_hi, F_twin, q_twin) by MCMC (emcee)."""
    import emcee
    q_obs = np.asarray(p["q"])
    keep = (q_obs >= Q_GRID[0]) & (q_obs <= 1.0)
    q_obs, S = q_obs[keep], s_matrix[keep]
    j = np.clip(np.searchsorted(Q_GRID, q_obs), 0, len(Q_GRID) - 1)
    s_at_obs = S[np.arange(len(q_obs)), j]
    good = s_at_obs > 0
    q_obs, S = q_obs[good], S[good]
    dq = np.gradient(Q_GRID)

    def lnprob(theta):
        g_lo, g_hi, ft, qt = theta
        if not (-6 < g_lo < 6 and -6 < g_hi < 6 and 0 <= ft < 0.6 and 0.93 <= qt < 0.995):
            return -np.inf
        pdf_obs = twin_model_pdf(q_obs, g_lo, g_hi, ft, qt, q_break)
        pdf_grid = twin_model_pdf(Q_GRID, g_lo, g_hi, ft, qt, q_break)
        norm = S @ (pdf_grid * dq)
        if np.any(pdf_obs <= 0) or np.any(norm <= 0):
            return -np.inf
        return float(np.sum(np.log(pdf_obs) - np.log(norm)))

    rng = np.random.default_rng(seed)
    p0 = np.column_stack([rng.normal(0, 0.3, n_walkers), rng.normal(-0.5, 0.3, n_walkers),
                          rng.uniform(0.01, 0.1, n_walkers), rng.uniform(0.94, 0.97, n_walkers)])
    sampler = emcee.EnsembleSampler(n_walkers, 4, lnprob)
    sampler.random_state = np.random.RandomState(seed).get_state()  # emcee copies numpy's global RNG
    sampler.run_mcmc(p0, n_steps, progress=False)
    chain = sampler.get_chain(discard=n_steps // 3, flat=True)
    ft = chain[:, 2]
    return {"n": int(len(q_obs)), "f_twin": float(np.median(ft)),
            "f_twin_lo": float(np.percentile(ft, 2.275)), "f_twin_hi": float(np.percentile(ft, 97.725)),
            "f_twin_16": float(np.percentile(ft, 16)), "f_twin_84": float(np.percentile(ft, 84)),
            "gamma_lo": float(np.median(chain[:, 0])), "gamma_hi": float(np.median(chain[:, 1])),
            "q_twin": float(np.median(chain[:, 3])), "chain": chain[::10]}
