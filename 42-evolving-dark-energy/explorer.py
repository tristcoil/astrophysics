"""An interactive explorer of the w0-wa plane (figures/dark_energy_explorer.html).

For every point of a grid in (w0, wa) and every headline data combination
(DESI + CMB, alone and with each supernova sample) we find the best values of
the other parameters (H0, omega_b, Omega_m) -- a *profile likelihood*. The page
then shows, for whichever (w0, wa) the reader picks:

* how much worse (or better) that universe fits than the best LCDM, per combination;
* the DESI distances and the binned supernova residuals against its prediction;
* the w(z) and dark energy density it implies.

The profile Delta chi^2 map is also a prior-free view of the evidence: unlike a
Bayesian posterior it does not depend on how much prior volume each region of
the plane has (one of the caveats debated after DESI's result).

    python explorer.py      # ~10-20 min on 6 cores; cached in data/derived/explorer_grid.npz
"""

from __future__ import annotations

import json
import multiprocessing as mp
import os
import time

import numpy as np
from scipy import optimize

import analysis as A
import background as bg
from config import C_KM_S, DERIVED_DIR, FIGURE_DIR, PROJECT_DIR

GRID_PATH = DERIVED_DIR / "explorer_grid.npz"
W0 = np.linspace(-1.4, 0.0, 36)
WA = np.linspace(-3.0, 1.2, 43)
COMBOS = [A.Combo(cmb="desi")] + [A.Combo(sn=s, cmb="desi") for s in A.SN_SAMPLES]
SCALE = np.array([1.0, 0.0002, 0.01])  # H0, ombh2, Om


def _profile_row(args):
    """Profile over (H0, ombh2, Om) along one row of wa, warm-starting each point."""
    ci, j, start = args
    from threadpoolctl import threadpool_limits

    combo = COMBOS[ci]
    with threadpool_limits(1):
        P = A.problem("w0wa", combo)
        out = np.full((W0.size, 3 + 3), np.nan)  # H0, ombh2, Om, chi2 bao/sn/cmb
        x = np.asarray(start, float)
        for i, w0 in enumerate(W0):
            wa = WA[j]
            if w0 + wa >= 0:
                continue

            def f(u):
                th = np.array([[*(u * SCALE), w0, wa]])
                return float(P.chi2_parts(th)["total"][0])

            res = optimize.minimize(f, x / SCALE, method="Nelder-Mead",
                                    options={"xatol": 1e-4, "fatol": 1e-4, "maxfev": 3000,
                                             "adaptive": True})
            if not np.isfinite(res.fun):
                continue
            x = res.x * SCALE
            parts = P.chi2_parts(np.array([[*x, w0, wa]]))
            out[i] = [*x, parts["bao"][0], parts["sn"][0] if "sn" in parts else 0.0,
                      parts["cmb"][0]]
        return ci, j, out


def build_grid(maps: dict, *, processes: int = 6, force: bool = False) -> dict:
    if GRID_PATH.exists() and not force:
        with np.load(GRID_PATH) as z:
            return {k: z[k] for k in z.files}
    t0 = time.monotonic()
    tasks = []
    for ci, c in enumerate(COMBOS):
        lc = maps[c.key]["lcdm"]
        for j in range(WA.size):
            tasks.append((ci, j, [lc["H0"], lc["ombh2"], lc["Om"]]))
    grid = np.full((len(COMBOS), WA.size, W0.size, 6), np.nan)
    with mp.get_context("fork").Pool(processes) as pool:
        for n, (ci, j, out) in enumerate(pool.imap_unordered(_profile_row, tasks)):
            grid[ci, j] = out
            if n % 40 == 0:
                print(f"    {n}/{len(tasks)} rows, {time.monotonic() - t0:.0f} s", flush=True)
    lcdm_chi2 = np.array([maps[c.key]["lcdm"]["chi2"] for c in COMBOS])
    out = {"grid": grid, "w0": W0, "wa": WA, "lcdm_chi2": lcdm_chi2}
    np.savez(GRID_PATH, **out)
    return out


def _binned_sn(s: str, ref: dict) -> dict:
    """Binned supernova residuals against the reference LCDM, with their covariance."""
    sn = A.sn_likelihood(s)
    c_ref = bg.Cosmo.from_Om(ref["Om"], ref["H0"] / 100, ref["ombh2"])
    g = bg.DistanceGrid(sn.z_hd, z_max=float(sn.z_hd.max()))
    dm = g.comoving(c_ref)[0]
    m_ref = sn.model_mu(dm[None, :])[0]
    if s == "Union3":
        cov = np.linalg.inv(sn.whiten.T @ sn.whiten)
        return {"z": sn.z_hd.tolist(), "d": (sn.mu - m_ref).tolist(), "cov": cov.tolist()}
    edges = np.array([0.01, 0.02, 0.03, 0.045, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
                      0.85, 1.0, 1.4, 2.4])
    which = np.digitize(sn.z_hd, edges) - 1
    used = [k for k in range(len(edges) - 1) if np.sum(which == k) > 0]
    B = np.zeros((sn.n, len(used)))
    for j, k in enumerate(used):
        B[which == k, j] = 1.0
    Aw = sn.whiten @ B
    cov = np.linalg.inv(Aw.T @ Aw)
    val = cov @ Aw.T @ (sn.whiten @ (sn.mu - m_ref))
    zb = [float(np.median(sn.z_hd[which == k])) for k in used]
    return {"z": zb, "d": val.tolist(), "cov": cov.tolist()}


def page_data(results: dict, grid: dict) -> dict:
    """Everything the page needs, as plain lists (a few hundred kB)."""
    import fetch_data as fd

    maps = results["maps"]
    ref = maps["desi_cmb-desi"]["lcdm"]
    data = fd.load_bao("dr2")
    early = A.early_universe()
    g = grid["grid"]
    # r_d at every grid point (depends on omega_b, omega_c only)
    H0, ob, Om = g[..., 0], g[..., 1], g[..., 2]
    ok = np.isfinite(H0)
    rd = np.full(H0.shape, np.nan)
    c = bg.Cosmo.from_Om(Om[ok], H0[ok] / 100, ob[ok])
    rd[ok] = early.rd(c)
    combos = []
    for ci, combo in enumerate(COMBOS):
        chi2 = g[ci, ..., 3:].sum(-1)
        combos.append({"label": " + ".join(["DESI", "CMB"] + ([A.SN_LABEL[combo.sn]] if combo.sn
                                                              else [])), "sn": combo.sn,
                       "dchi2": np.round(chi2 - grid["lcdm_chi2"][ci], 3).tolist(),
                       "H0": np.round(H0[ci], 4).tolist(), "ombh2": np.round(ob[ci], 6).tolist(),
                       "Om": np.round(Om[ci], 5).tolist(), "rd": np.round(rd[ci], 4).tolist(),
                       "lcdm_chi2": float(grid["lcdm_chi2"][ci]),
                       "w0wa_best": maps[combo.key]["w0wa"], "sigma": maps[combo.key]["sigma"],
                       "dchi2_best": maps[combo.key]["dchi2"]})
    c_ref = bg.Cosmo.from_Om(ref["Om"], ref["H0"] / 100, ref["ombh2"])
    sn = {s: _binned_sn(s, ref) for s in A.SN_SAMPLES}
    return {"w0": grid["w0"].tolist(), "wa": grid["wa"].tolist(), "combos": combos,
            "bao": {"z": data["z"].tolist(), "value": data["value"].tolist(),
                    "kind": data["kind"].tolist(),
                    "err": np.sqrt(np.diag(data["cov"])).tolist()},
            "ref": {**ref, "rd": float(early.rd(c_ref)[0])}, "sn": sn,
            "labels": {s: A.SN_LABEL[s] for s in A.SN_SAMPLES},
            "omega_r": bg.OMEGA_R, "omega_nu": bg.OMEGA_NU, "c": C_KM_S}


def write_page(payload: dict) -> None:
    template = PROJECT_DIR / "explorer_template.html"
    html = template.read_text().replace("/*__DATA__*/null", json.dumps(payload, separators=(",", ":")))
    out = FIGURE_DIR / "dark_energy_explorer.html"
    out.write_text(html)
    print(f"  wrote {out} ({out.stat().st_size / 1e3:.0f} kB)")


def main(processes: int = 6, force: bool = False) -> None:
    results = json.loads(A.RESULTS_PATH.read_text())
    print("Profile-likelihood grid over (w0, wa)")
    grid = build_grid(results["maps"], processes=processes, force=force)
    write_page(page_data(results, grid))


if __name__ == "__main__":
    main(processes=int(os.environ.get("D5_PROCESSES", "6")))
