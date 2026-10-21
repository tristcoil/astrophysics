"""The 1922 Wallal eclipse (Lick Observatory, Western Australia), replicated and refitted.

Campbell and Trumpler photographed the eclipsed Sun on 21 September 1922 with
two instruments: two 15-foot cameras (92 stars on four plate pairs; 1.72" +- 0.11",
probable error; Lick Obs. Bull. 11, 41, 1923) and a pair of five-foot cameras with
a wider field (145 stars on six plate pairs; 1.82" +- 0.15"; Lick Obs. Bull. 13,
130, 1928), combined to 1.75" +- 0.09".

Their displacements are *relative*: each plate's scale was fixed so that the
reference stars far from the Sun (> 5 deg for the five-foot cameras, > 2 deg for
the 15-foot ones) average to zero. A deflection beta R/r then shows up together
with a small residual scale term, and each plate pair is solved for both:

    D_r = beta R / r - alpha r        (R = the Sun's radius, 0.2661 deg; r in deg)

with each star's printed weight (weight 1 = a probable error of 0.33" for the
five-foot cameras). The deflection at the limb is beta.

Here: the 1928 Table 9 solutions redone from Table 6; our own errors (from the
weights and from the scatter); a test of the 1/r law (beta (R/r)^n); the 15-foot
cameras' stars (1923, Table II) fitted the same way; and the two combined.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import config as C

R_SUN_DEG = 0.2661                 # the Sun's radius adopted by Campbell & Trumpler (p.150)
PE_UNIT_5FT = 0.33                 # probable error of a radial component of weight 1, five-foot (p.138)
UNIT_15FT_DEG = 1.25222            # the 1923 coordinate unit: 100 mm = 4508" (p.50)
PE = 0.6745


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(C.TRANSCRIBED_DIR / name, comment="#", skipinitialspace=True)


def table6() -> pd.DataFrame:
    return _read("wallal1928_table6_radial.csv")


def pairs(t6: pd.DataFrame) -> list[str]:
    return [c[3:] for c in t6.columns if c.startswith("Dr_") and c != "Dr_mean"]


def solve(r_deg, dr, w, n: float = 1.0, scale: bool = True) -> dict:
    """Weighted least squares of D_r = beta (R/r)^n - alpha r. Returns beta, alpha, the
    weight of beta (1/variance in units of the unit weight), the residual variance of
    unit weight, and the number of stars."""
    ok = np.isfinite(dr) & np.isfinite(w) & (w > 0)
    r, y, w = np.asarray(r_deg)[ok], np.asarray(dr)[ok], np.asarray(w)[ok]
    cols = [(R_SUN_DEG / r) ** n]
    if scale:
        cols.append(-r)
    A = np.column_stack(cols)
    N = (A * w[:, None]).T @ A
    beta = np.linalg.solve(N, (A * w[:, None]).T @ y)
    res = y - A @ beta
    dof = len(y) - A.shape[1]
    cov_unit = np.linalg.inv(N)
    return {"beta": float(beta[0]), "alpha": float(beta[1]) if scale else 0.0, "weight": float(1 / cov_unit[0, 0]),
            "s2_unit": float(np.sum(w * res ** 2) / dof), "n_stars": int(len(y)), "dof": dof}


def replicate(t6: pd.DataFrame | None = None) -> dict:
    """Table 9 (p.151) from Table 6: beta, weight and alpha per plate pair; the weighted
    mean of beta and its probable error from the adopted unit weight (0.33")."""
    t6 = table6() if t6 is None else t6
    r = t6["distance_deg"].to_numpy(float)
    rows = []
    for p in pairs(t6):
        s = solve(r, t6[f"Dr_{p}"].to_numpy(float), t6[f"wt_{p}"].to_numpy(float))
        rows.append({"pair": p, **s})
    per = pd.DataFrame(rows)
    w = per["weight"].to_numpy()
    mean = float(np.sum(w * per["beta"]) / w.sum())
    return {"per_pair": per, "beta": mean, "pe": float(PE_UNIT_5FT / np.sqrt(w.sum())),
            "alpha": float(np.sum(w * per["alpha"]) / w.sum())}


def fit_5ft(t6: pd.DataFrame | None = None, n_free: bool = False) -> dict:
    """Our fit: the six plate pairs together, one beta (and n), each pair its own alpha;
    one-sigma errors both from the adopted weights and from the scatter (chi2/dof)."""
    from scipy.optimize import least_squares

    t6 = table6() if t6 is None else t6
    r = t6["distance_deg"].to_numpy(float)
    data = []
    for k, p in enumerate(pairs(t6)):
        d, w = t6[f"Dr_{p}"].to_numpy(float), t6[f"wt_{p}"].to_numpy(float)
        ok = np.isfinite(d) & np.isfinite(w) & (w > 0)
        data.append((k, r[ok], d[ok], w[ok]))
    npair = len(data)
    sig_unit = PE_UNIT_5FT / PE

    def resid(par):
        beta, n = par[0], (par[1] if n_free else 1.0)
        alphas = par[2:] if n_free else par[1:]
        out = []
        for k, rr, dd, ww in data:
            out.append(np.sqrt(ww) / sig_unit * (dd - beta * (R_SUN_DEG / rr) ** n + alphas[k] * rr))
        return np.concatenate(out)

    p0 = [1.75, 1.0] + [0.01] * npair if n_free else [1.75] + [0.01] * npair
    sol = least_squares(resid, p0, xtol=1e-12, ftol=1e-12, gtol=1e-12)
    J = sol.jac
    cov = np.linalg.inv(J.T @ J)
    chi2 = float(np.sum(sol.fun ** 2))
    dof = len(sol.fun) - len(p0)
    scale = max(1.0, np.sqrt(chi2 / dof))
    out = {"beta": float(sol.x[0]), "beta_err_weights": float(np.sqrt(cov[0, 0])),
           "beta_err": float(np.sqrt(cov[0, 0]) * scale), "chi2": chi2, "dof": dof,
           "n_measures": int(len(sol.fun))}
    if n_free:
        out.update({"n": float(sol.x[1]), "n_err": float(np.sqrt(cov[1, 1]) * scale),
                    "corr_beta_n": float(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1]))})
    return out


def fit_15ft() -> dict:
    """The 15-foot cameras' 92 stars (1923, Table II): the corrected radial components D2
    (mean of the four plate pairs and both observers), with their weights, fitted with
    beta R/r - alpha r. Errors from the scatter (the 1923 unit weight is not restated)."""
    t = _read("wallal1923_table2_displacements.csv")
    r = t["dist"].to_numpy(float) * UNIT_15FT_DEG
    out = {}
    for col in ("D1", "D2"):
        s = solve(r, t[col].to_numpy(float), t["weight"].to_numpy(float))
        err = float(np.sqrt(s["s2_unit"] / s["weight"]))
        out[col] = {"beta": s["beta"], "beta_err": err, "alpha": s["alpha"], "n_stars": s["n_stars"]}
    return out


def stars_for_plot(t6: pd.DataFrame | None = None, alpha: float | None = None) -> pd.DataFrame:
    """Each five-foot star's mean radial displacement with the mean scale term added back
    (so that it can be compared with beta R/r), its one-sigma error from its total weight."""
    t6 = table6() if t6 is None else t6
    a = replicate(t6)["alpha"] if alpha is None else alpha
    r = t6["distance_deg"].to_numpy(float)
    w = t6["wt_total"].to_numpy(float)
    return pd.DataFrame({"no": t6["no"], "r_deg": r, "inv_r": R_SUN_DEG / r,
                         "dr": t6["Dr_mean"].to_numpy(float) + a * r,
                         "err": (PE_UNIT_5FT / PE) / np.sqrt(w)})


def analyse() -> dict:
    t6 = table6()
    rep = replicate(t6)
    fit = fit_5ft(t6)
    pw = fit_5ft(t6, n_free=True)
    f15 = fit_15ft()
    # the two instruments combined, as in 1928 (where the 15-foot cameras had twice the 5-foot's
    # weight). The published 15-foot result (1.72") used the uncorrected radial components D1;
    # the check-star-corrected D2 gave 2.05" in 1923 (p.54), which the authors called an extrapolation.
    b5, e5 = fit["beta"], fit["beta_err"]
    b15, e15 = f15["D1"]["beta"], f15["D1"]["beta_err"]
    w5, w15 = 1 / e5 ** 2, 1 / e15 ** 2
    comb = {"beta": float((w5 * b5 + w15 * b15) / (w5 + w15)), "beta_err": float(1 / np.sqrt(w5 + w15)),
            "chi2": float(w5 * (b5 - b15) ** 2 * w15 / (w5 + w15))}
    return {"replication": {"per_pair": rep["per_pair"].to_dict("records"), "beta": rep["beta"], "pe": rep["pe"],
                            "alpha": rep["alpha"]},
            "fit": {"limb": fit["beta"], "limb_err": fit["beta_err"], "limb_err_weights": fit["beta_err_weights"],
                    "chi2": fit["chi2"], "dof": fit["dof"], "n_measures": fit["n_measures"]},
            "fit_power": pw, "fit_15ft": f15, "combined": comb,
            "stars": stars_for_plot(t6).to_dict("records"), "R_deg": R_SUN_DEG}


if __name__ == "__main__":
    import json

    a = analyse()
    print(pd.DataFrame(a["replication"]["per_pair"]).round(4).to_string(index=False))
    print(json.dumps({k: v for k, v in a.items() if k not in ("stars", "replication")}, indent=1))
    print("replication: beta", round(a["replication"]["beta"], 3), "pe", round(a["replication"]["pe"], 3),
          "alpha", round(a["replication"]["alpha"], 4))
