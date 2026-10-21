"""Contact binaries by the hundred thousand: periods, the 0.22-day edge, colours.

The ZTF catalogue of periodic variables (Chen, Wang & Deng 2020) has 369,707
EW-type (W UMa) binaries. Their period distribution has a famous sharp edge
near 0.22 days (Rucinski 1992, 2007): almost no contact binaries orbit faster.

Before reading physics into the edge we test it against selection effects:
the same edge must appear at every apparent magnitude (the survey's depth
cannot carve a feature at a fixed period), in absolute-magnitude slices
(Gaia distances), and in Gaia's independent eclipsing-binary census.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

import config as C


def edge_model(logp, amp, logp_edge, width, slope, curv):
    """Smooth step (logistic in log P) times a curved power law: the density of
    contact binaries per unit log P near the short-period end."""
    step = 1.0 / (1.0 + np.exp(-(logp - logp_edge) / width))
    x = logp - (-0.55)
    return amp * step * np.exp(slope * x + curv * x * x)


def fit_edge(periods, logp_lo=-0.85, logp_hi=-0.35, nbin=100):
    """Fit ``edge_model`` to the histogram of log P (Poisson likelihood, via
    least squares on Anscombe-transformed counts). Returns the edge period,
    its width (in log P), the half-density point and bootstrap-free errors from
    the Jacobian."""
    lp = np.log10(periods)
    edges = np.linspace(logp_lo, logp_hi, nbin + 1)
    n, _ = np.histogram(lp, edges)
    mid = 0.5 * (edges[1:] + edges[:-1])
    y = 2 * np.sqrt(n + 3 / 8)

    def resid(x):
        m = edge_model(mid, *x) * (edges[1] - edges[0])
        return 2 * np.sqrt(np.maximum(m, 0) + 3 / 8) - y

    x0 = [n.max() / (edges[1] - edges[0]), np.log10(0.225), 0.01, 0.0, 0.0]
    r = least_squares(resid, x0, bounds=([0, -0.8, 1e-4, -50, -200], [np.inf, -0.5, 0.2, 50, 200]))
    J = r.jac
    cov = np.linalg.inv(J.T @ J)
    err = np.sqrt(np.diag(cov))
    return {"p_edge": 10 ** r.x[1], "logp_edge": r.x[1], "logp_edge_err": err[1], "width_dex": r.x[2],
            "width_err": err[2], "params": r.x, "hist": (mid, n), "chi2": float(np.sum(r.fun**2)), "dof": nbin - 5}


def quantiles(periods, qs=(0.001, 0.005, 0.01, 0.05, 0.5)):
    return {q: float(np.quantile(periods, q)) for q in qs}


def ew_sample(ztf: pd.DataFrame, amp_min: float = 0.08) -> pd.DataFrame:
    """EW binaries with both-band light curves; Chen et al.'s amplitude limit."""
    s = ztf[(ztf["type"] == "EW") & (ztf[["amp_g", "amp_r"]].max(axis=1) >= amp_min)].copy()
    s["logp"] = np.log10(s["period"])
    s["g_r"] = s["gmag"] - s["rmag"]
    return s


def mag_slices(ew: pd.DataFrame, edges=(12.5, 15.0, 16.5, 17.5, 18.5, 20.5), band="rmag"):
    """Edge position in slices of apparent magnitude."""
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = ew[(ew[band] >= lo) & (ew[band] < hi)]
        if len(s) < 2000:
            continue
        f = fit_edge(s["period"].to_numpy())
        out.append({"lo": lo, "hi": hi, "n": len(s), **{k: f[k] for k in ("p_edge", "logp_edge_err", "width_dex")},
                    "frac_below_022": float(np.mean(s["period"] < 0.22))})
    return out
