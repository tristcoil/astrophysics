"""sigma_8 from the CMB: project D3's own fit, run forward to today.

Project D3 (09-planck-cmb) fitted the Planck temperature power spectrum and
measured the universe's ingredients at z ~ 1090. The CMB fixes the amplitude
of the primordial fluctuations, A_s, and everything needed to grow them to
today *if* LCDM is right. Running that growth forward with CAMB gives the
sigma_8 that the early universe predicts for the late one: the number that
cluster counts test.

D3 fitted (Omega_b h^2, Omega_c h^2, 100 theta_MC, ln 10^10 A_s, n_s) with the
optical depth tau fixed at 0.0522. Temperature data constrain only
A_s exp(-2 tau), so here tau is drawn from Planck's TT+lowE constraint,
0.0522 +- 0.0080 (from large-scale polarisation, which D3 did not use), and
A_s rescaled to keep A_s exp(-2 tau) at D3's value. D3's parameter
covariance (its Gaussian approximation) gives the other draws. Same CAMB
settings as D3: one massive neutrino of 0.06 eV, N_eff = 3.046.
"""

from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import config as C
import literature as L

CACHE = C.CACHE_DIR / "d3_sigma8_samples.npz"
NAMES = ["ombh2", "omch2", "theta100", "logA", "ns"]


def d3_fit() -> dict:
    """D3's best fit and correlation matrix (live file if present, else the copy)."""
    if C.D3_FIT_PATH.exists():
        d = json.loads(C.D3_FIT_PATH.read_text())
        mean = np.array([d["parameters"][n]["value"] for n in NAMES])
        err = np.array([d["parameters"][n]["error"] for n in NAMES])
        corr = np.array(d["correlation"]["matrix"])
        return {"mean": mean, "err": err, "corr": corr, "tau_fixed": d["tau_fixed"], "source": str(C.D3_FIT_PATH)}
    d = L.D3_TT_FIT
    return {"mean": np.array(d["mean"]), "err": np.array(d["err"]), "corr": np.array(d["corr"]),
            "tau_fixed": d["tau_fixed"], "source": "literature.D3_TT_FIT (copy)"}


def _sigma8(p):
    import camb

    ombh2, omch2, theta, logA, ns, tau = p
    pars = camb.set_params(ombh2=ombh2, omch2=omch2, cosmomc_theta=theta / 100, tau=tau,
                           As=np.exp(logA) * 1e-10, ns=ns, mnu=0.06, nnu=3.046,
                           WantTransfer=True, redshifts=[0.0], kmax=10.0)
    pars.WantCls = False
    res = camb.get_results(pars)
    h = pars.H0 / 100
    om = (pars.ombh2 + pars.omch2 + pars.omnuh2) / h**2
    return float(res.get_sigma8_0()), float(om), float(h)


def sample(n: int = 300, seed: int = 31, processes: int = 6) -> dict:
    if CACHE.exists():
        return dict(np.load(CACHE))
    fit = d3_fit()
    rng = np.random.default_rng(seed)
    cov = fit["corr"] * np.outer(fit["err"], fit["err"])
    draws = rng.multivariate_normal(fit["mean"], cov, size=n)
    tau = rng.normal(L.PLANCK_TAU_TT_LOWE[0], L.PLANCK_TAU_TT_LOWE[1], size=n)
    draws[:, 3] += 2 * (tau - fit["tau_fixed"])  # keep A_s exp(-2 tau) fixed
    args = [tuple(d) + (t,) for d, t in zip(draws, tau)]
    import multiprocessing as mp

    with ProcessPoolExecutor(processes, mp_context=mp.get_context("fork")) as ex:
        out = np.array(list(ex.map(_sigma8, args)))
    best = _sigma8(tuple(fit["mean"]) + (fit["tau_fixed"],))
    res = {"sigma8": out[:, 0], "Om": out[:, 1], "h": out[:, 2], "tau": tau,
           "best_sigma8": best[0], "best_Om": best[1], "best_h": best[2]}
    np.savez(CACHE, **res)
    return res


def summary(s: dict) -> dict:
    S8 = s["sigma8"] * np.sqrt(s["Om"] / 0.3)
    q = lambda v: (float(np.median(v)), float(np.std(v)))  # noqa: E731
    return {"sigma8": q(s["sigma8"]), "Om": q(s["Om"]), "S8": q(S8),
            "corr_Om_sigma8": float(np.corrcoef(s["Om"], s["sigma8"])[0, 1]),
            "best_sigma8": float(s["best_sigma8"]), "n": int(s["sigma8"].size)}


def planck_check() -> float:
    """sigma_8 for Planck 2018's own best fit, through this same code path."""
    p = L.PLANCK2018_BESTFIT
    return _sigma8((p["ombh2"], p["omch2"], p["theta100"], np.log(p["As"] * 1e10), p["ns"], p["tau"]))[0]
