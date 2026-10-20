"""The first clock: deuterium from the first three minutes.

Deuterium is the "baryometer": the more baryons there were per photon, the
more completely deuterium burned into helium, so the leftover D/H falls
steeply with eta (roughly as eta^-1.6). We

1. run the BBN code PRyMordial (Burns, Tait & Valli 2024) over a grid of eta
   with both of its nuclear-rate compilations -- "PRIMAT" (theory-based
   d + d rates, Pitrou et al.) and "NACRE II" (data-driven), both with LUNA's
   measured d(p,gamma)3He rate -- because the choice moves the answer by
   about 2 sigma, the main systematic of the whole method;
2. run each 60 times with the 12 key reaction rates and the neutron lifetime
   drawn from their uncertainties, for the theory error on D/H;
3. combine with the measured primordial D/H (12 quasar absorbers, PDG
   2025/26 selection) to get eta and Omega_b h^2, with the theory error
   added in quadrature.

Runs go through ``prym_runner.py`` in isolated subprocesses (PRyMordial is
configured by module globals and reads its tables from the current
directory). Results are cached in ``data/derived/bbn/``.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy import interpolate, optimize

import baryons
import config as cfg
import literature as lit

OUT = cfg.DERIVED_DIR / "bbn"
OUT.mkdir(parents=True, exist_ok=True)
RATE_SETS = ("primat", "nacreii")
ETA_GRID = np.geomspace(4.0e-10, 8.0e-10, 21)
MC_ETA = 6.1e-10
MC_DRAWS = 60
MC_CHUNKS = 3


def _run(rates: str, mode_args: list[str], net: str = "full") -> list[dict]:
    cmd = [sys.executable, "-I", str(cfg.PROJECT_DIR / "prym_runner.py"),
           str(cfg.PRYMORDIAL_DIR), rates, net] + mode_args
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    if res.returncode != 0:
        raise RuntimeError(f"PRyMordial failed: {res.stderr[-2000:]}")
    return [json.loads(line) for line in res.stdout.splitlines() if line.startswith("{")]


def run_all(force: bool = False, workers: int = 6) -> dict:
    """Grid and Monte Carlo for both rate sets (cached as JSON lines)."""
    jobs = []
    for rates in RATE_SETS:
        path = OUT / f"grid_{rates}.jsonl"
        if force or not path.exists():
            halves = np.array_split(ETA_GRID, 2)
            for k, part in enumerate(halves):
                jobs.append((path, k, rates, ["grid"] + [f"{e:.6e}" for e in part]))
        path = OUT / f"mc_{rates}.jsonl"
        if force or not path.exists():
            per = MC_DRAWS // MC_CHUNKS
            for k in range(MC_CHUNKS):
                seed = 1000 * (1 + RATE_SETS.index(rates)) + k
                jobs.append((path, k, rates, ["mc", f"{MC_ETA:.6e}", str(per), str(seed)]))
    if jobs:
        print(f"  BBN: {len(jobs)} PRyMordial jobs ({workers} at a time, ~10 s per run)")
        with ThreadPoolExecutor(workers) as ex:
            outs = list(ex.map(lambda j: (j[0], j[1], _run(j[2], j[3])), jobs))
        by_path: dict = {}
        for path, k, rows in outs:
            by_path.setdefault(path, {})[k] = rows
        for path, parts in by_path.items():
            rows = [r for k in sorted(parts) for r in parts[k]]
            path.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    out = {}
    for rates in RATE_SETS:
        out[rates] = {
            "grid": [json.loads(x) for x in (OUT / f"grid_{rates}.jsonl").read_text().splitlines()],
            "mc": [json.loads(x) for x in (OUT / f"mc_{rates}.jsonl").read_text().splitlines()],
        }
    return out


class DeuteriumCurve:
    """D/H(eta) for one rate set: a spline in log-log, plus the fractional
    theory error from the Monte Carlo at eta = 6.1e-10."""

    def __init__(self, grid: list[dict], mc: list[dict]):
        eta = np.array([r["eta"] for r in grid])
        dh = np.array([r["DoH_e5"] for r in grid])
        o = np.argsort(eta)
        self.eta, self.dh = eta[o], dh[o]
        self._spl = interpolate.CubicSpline(np.log(self.eta), np.log(self.dh))
        self.yp = np.array([r["Yp_bbn"] for r in grid])[o]
        self.li7 = np.array([r["Li7oH_e10"] for r in grid])[o]
        self.he3 = np.array([r["He3oH_e5"] for r in grid])[o]
        d_mc = np.array([r["DoH_e5"] for r in mc])
        self.mc_dh = d_mc
        self.frac_theory = float(np.std(d_mc, ddof=1) / np.mean(d_mc))
        self.mc_mean = float(np.mean(d_mc))
        self.slope = float(self._spl(math.log(MC_ETA), 1))       # dlnDH/dln eta

    def __call__(self, eta):
        return np.exp(self._spl(np.log(eta)))

    def eta_from_dh(self, dh_obs: float, dh_err: float, n: int = 4001) -> dict:
        """Posterior for eta (flat prior in ln eta) from a Gaussian D/H
        measurement plus a log-normal theory error."""
        grid = np.geomspace(4.2e-10, 7.8e-10, n)
        pred = self(grid)
        sig_ln = math.hypot(dh_err / dh_obs, self.frac_theory)
        lnl = -0.5 * ((np.log(dh_obs) - np.log(pred)) / sig_ln) ** 2
        w = np.exp(lnl - lnl.max())
        w /= np.trapezoid(w, np.log(grid))
        cdf = np.concatenate([[0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(np.log(grid)))])
        q = np.interp([0.16, 0.5, 0.84], cdf, grid)
        best = float(optimize.brentq(lambda e: math.log(self(e)) - math.log(dh_obs), 4.3e-10, 7.7e-10))
        return {"eta_best": best, "eta_median": float(q[1]), "eta_lo": float(q[0]),
                "eta_hi": float(q[2]), "eta_err": float(0.5 * (q[2] - q[0])),
                "sigma_ln_total": sig_ln, "frac_theory": self.frac_theory}


def analyse(force: bool = False) -> dict:
    runs = run_all(force=force)
    curves = {r: DeuteriumCurve(runs[r]["grid"], runs[r]["mc"]) for r in RATE_SETS}
    dh_all = lit.dh_weighted_mean()
    dh_cooke = lit.dh_log_mean([r for r in lit.DH_SYSTEMS if r[5]])
    out = {"dh_measured": {"pdg2026_12": dh_all, "cooke2018_7": dh_cooke}, "rates": {}}
    for rates, cv in curves.items():
        res = {}
        for key, (m, e) in {"pdg2026_12": (dh_all["mean"], dh_all["err_scaled"]),
                            "cooke2018_7": (dh_cooke["mean"], dh_cooke["err"])}.items():
            r = cv.eta_from_dh(m, e)
            for k in ("eta_best", "eta_median", "eta_lo", "eta_hi", "eta_err"):
                r[k.replace("eta", "ombh2")] = float(baryons.ombh2_from_eta(r[k]))
            res[key] = r
        res["dh_at_planck_eta"] = float(cv(baryons.eta_from_ombh2(cfg.OMBH2_PLANCK)))
        res["frac_theory"] = cv.frac_theory
        res["slope_dlnDH_dlneta"] = cv.slope
        res["yp_at_planck_eta"] = float(np.interp(baryons.eta_from_ombh2(cfg.OMBH2_PLANCK),
                                                  cv.eta, cv.yp))
        res["li7_at_planck_eta"] = float(np.interp(baryons.eta_from_ombh2(cfg.OMBH2_PLANCK),
                                                   cv.eta, cv.li7))
        out["rates"][rates] = res
    return out, curves


if __name__ == "__main__":
    import pprint

    r, _ = analyse()
    pprint.pprint(r)
