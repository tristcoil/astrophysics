"""The Principe plates of 1919, replicated and refitted.

Eddington and Cottingham's two usable eclipse plates (X and W, 5 stars each)
were measured (X at Principe a few days after the eclipse, W at Cambridge on
22-23 August; p.321) against comparison plates of the same field taken
at Oxford in January (G1, H1 for X; D1, I2 for W; plate U against K2 was
judged worthless). All numbers are in units of 0.003" (the fifth decimal of
the 5-mm reseau interval) and plate coordinates in reseau intervals (~5'); we
find the latter are exactly 10 x Table I's coordinates + 14.

The 1920 reduction (pp.321-328), per comparison:

1. remove the "theoretical plate constants": differential refraction and
   aberration, plus the change of scale measured on check plates of a field
   near Arcturus taken at both places (a x + b y in x, d x + e y in y, b = d);
2. find the orientation theta from the x measures, Delta2 x = c + theta y,
   after removing 3/4 of Einstein's x-displacement (a compromise: the x
   measures were to fix only the orientation);
3. from the y measures with the orientation removed, solve
   Delta4 y = df + dk alpha_y, where Delta4 y = Delta3 y - Ey has Einstein's
   full y-displacement (Ey = 184 alpha_y) taken out: kappa = 184 + dk is the
   deflection at 50' in units of 0.003", and 1.75" x kappa/184 at the limb.

Results: X-G1 1.94", X-H1 1.44", W-D1 1.55", W-I2 1.67"; mean 1.65"; with a
-0.04" correction for the 1/4 Ex neglected in step 2, 1.61" +- 0.30" (probable
error: 0.25" from the measures, 0.14" from the scale, common to all four).

Our refit solves each comparison's x and y measures together for the zero
points, the orientation and the deflection, with the exact GR pattern of
the five stars (Gaia positions, JPL Sun) and the check-plate scale, and
(as a variant) with the scale left free.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import config as C

COMPARISONS = ["X-G1", "X-H1", "W-D1", "W-I2"]          # used in 1920; U-K2 was set aside
ALL_COMPARISONS = COMPARISONS + ["U-K2"]
PLATE_OF = {"X-G1": "X", "X-H1": "X", "W-D1": "W", "W-I2": "W", "U-K2": "U"}
UNIT = 0.003                                              # arcsec per measuring unit
KAPPA_FULL = 184.0                                        # Einstein's 1.75" at the limb, at 50'
SCALE_PE_LIMB = 0.14                                      # p.320/327: scale error, common to all four


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(C.TRANSCRIBED_DIR / name, comment="#", skipinitialspace=True)


def load() -> dict:
    """Measures, star coordinates (reseau units), theoretical constants and the
    1920 Einstein displacements (Ex, Ey in 0.003" units)."""
    p = _read("dyson1920_principe_plates.csv")
    k = _read("dyson1920_principe_plate_constants.csv")
    t1 = _read("dyson1920_table1_stars.csv").set_index("no")
    stars = sorted(p["star"].unique())
    xy = pd.DataFrame({"x": 10 * t1.loc[stars, "x"] + 14, "y": 10 * t1.loc[stars, "y"] + 14}, index=stars)
    # Einstein's full displacements as used in 1920: Ey printed in every y table; Ex printed
    # in full only for U-K2, otherwise as 3/4 Ex (X-G1, W-D1)
    ey = p[p["coord"] == "y"].dropna(subset=["Ey"]).groupby("star")["Ey"].first()
    ex_full = p[(p["coord"] == "x")].dropna(subset=["Ex"]).groupby("star")["Ex"].first()
    ex_34 = p[(p["coord"] == "x")].dropna(subset=["three_quarter_Ex"]).groupby("star")["three_quarter_Ex"].first()
    ex = pd.Series({s: ex_full.get(s, np.nan) if s in ex_full.index else ex_34.get(s) * 4 / 3 for s in stars})
    const = {}
    for comp in ALL_COMPARISONS:
        sub = k[(k["context"] == comp)].set_index("quantity")["value"]
        const[comp] = {"a": float(sub["theoretical_a"]), "b": float(sub["theoretical_b=d"]),
                       "e": float(sub["theoretical_e"])}
    printed = {comp: float(k[(k["context"] == comp) & (k["quantity"] == "delta_kappa")]["value"].iloc[0])
               for comp in ALL_COMPARISONS}
    return {"plates": p, "xy": xy, "Ex": ex, "Ey": ey, "const": const, "printed_dk": printed}


def _comparison(d: dict, comp: str) -> pd.DataFrame:
    """One comparison's stars with x, y, Delta x, Delta y and Delta1 (theoretical constants removed)."""
    p = d["plates"][d["plates"]["comparison"] == comp]
    px = p[p["coord"] == "x"].set_index("star")
    py = p[p["coord"] == "y"].set_index("star")
    stars = [s for s in px.index if s in py.index]
    c = d["const"][comp]
    x, y = d["xy"].loc[stars, "x"].to_numpy(), d["xy"].loc[stars, "y"].to_numpy()
    dx, dy = px.loc[stars, "delta"].to_numpy(float), py.loc[stars, "delta"].to_numpy(float)
    return pd.DataFrame({"star": stars, "x": x, "y": y, "dx": dx, "dy": dy,
                         "d1x": dx - (c["a"] * x + c["b"] * y), "d1y": dy - (c["b"] * x + c["e"] * y),
                         "d1x_printed": px.loc[stars, "delta1"].to_numpy(float),
                         "d1y_printed": py.loc[stars, "delta1"].to_numpy(float),
                         "Ex": d["Ex"].loc[stars].to_numpy(float), "Ey": d["Ey"].loc[stars].to_numpy(float)})


def replicate(d: dict | None = None) -> dict:
    """The 1920 steps for each comparison. The deflection step is redone exactly from
    the printed Delta4 y (it reproduces the printed dk); the orientation step is also
    redone by plain least squares, which does not give the printed theta (see notes)."""
    d = load() if d is None else d
    k = _read("dyson1920_principe_plate_constants.csv")
    rows = []
    for comp in COMPARISONS:
        t = _comparison(d, comp)
        p = d["plates"][d["plates"]["comparison"] == comp]
        d4y_printed = p[p["coord"] == "y"].set_index("star").loc[t["star"], "delta4"].to_numpy(float)
        theta_printed = float(k[(k["context"] == comp) & (k["quantity"] == "orientation_theta")]["value"].iloc[0])
        ay = t["Ey"].to_numpy() / KAPPA_FULL
        A = np.column_stack([np.ones(len(t)), ay])
        (df, dk), *_ = np.linalg.lstsq(A, d4y_printed, rcond=None)
        # the orientation by plain least squares from Delta2 x = Delta1 x - 3/4 Ex = c + theta y
        d2x = t["d1x_printed"] - 0.75 * t["Ex"]
        (cx, theta_ls), *_ = np.linalg.lstsq(np.column_stack([np.ones(len(t)), t["y"]]), d2x, rcond=None)
        d4y_ls = t["d1y_printed"] + theta_ls * t["x"] - t["Ey"]
        (df2, dk_ls), *_ = np.linalg.lstsq(A, d4y_ls, rcond=None)
        rows.append({"comparison": comp, "theta_printed": theta_printed, "theta_ls": theta_ls,
                     "dk": dk, "dk_printed": d["printed_dk"][comp],
                     "limb": C.DEFLECTION_LIMB_1919 * (KAPPA_FULL + dk) / KAPPA_FULL,
                     "limb_theta_ls": C.DEFLECTION_LIMB_1919 * (KAPPA_FULL + dk_ls) / KAPPA_FULL,
                     "d1_max_abs_diff": float(max(np.abs(t["d1x"] - t["d1x_printed"]).max(),
                                                  np.abs(t["d1y"] - t["d1y_printed"]).max()))})
    per = pd.DataFrame(rows)
    mean = float(per["limb"].mean())
    return {"per_comparison": per, "mean": mean, "limb": mean - 0.04,
            "mean_theta_ls": float(per["limb_theta_ls"].mean())}


def gr_pattern(stars) -> pd.DataFrame:
    """Exact GR displacements (0.003" units) of the given stars at Principe's mid-totality."""
    import stars as S

    f = S.fit_table1_frame()
    pat = S.predicted_pattern(f, "principe")
    return pd.DataFrame({"Ex": pat.loc[stars, "dx"] / UNIT, "Ey": pat.loc[stars, "dy"] / UNIT}, index=stars)


def _design(t: pd.DataFrame, E: pd.DataFrame, free_scale: bool) -> tuple[np.ndarray, np.ndarray]:
    """Rows: x measures then y measures. Unknowns: c_x, c_y, theta, lambda (+ da, de)."""
    n = len(t)
    one, zero = np.ones(n), np.zeros(n)
    x, y = t["x"].to_numpy(), t["y"].to_numpy()
    ex, ey = E.loc[t["star"], "Ex"].to_numpy(), E.loc[t["star"], "Ey"].to_numpy()
    cols_x = [one, zero, y, ex]
    cols_y = [zero, one, -x, ey]
    if free_scale:               # corrections to the check-plate scale, separately in x and y
        cols_x += [x, zero]
        cols_y += [zero, y]
    A = np.vstack([np.column_stack(cols_x), np.column_stack(cols_y)])
    b = np.concatenate([t["d1x"].to_numpy(), t["d1y"].to_numpy()])
    return A, b


def fit(d: dict | None = None, free_scale: bool = False, comparisons=COMPARISONS, pattern: str = "gr") -> dict:
    """Least squares per comparison (both coordinates). lambda = deflection in units of
    the pattern's (GR's exact pattern: 1.7512" at the limb; '1920': Ex, Ey for 1.75").
    The measuring noise is one sigma per coordinate, pooled over the comparisons."""
    d = load() if d is None else d
    comps = {c: _comparison(d, c) for c in comparisons}
    stars = sorted({s for t in comps.values() for s in t["star"]})
    E = gr_pattern(stars) if pattern == "gr" else pd.DataFrame({"Ex": d["Ex"], "Ey": d["Ey"]}).loc[stars]
    limb_unit = C.DEFLECTION_LIMB_GR if pattern == "gr" else C.DEFLECTION_LIMB_1919
    sols, ssx, ssy, nx, ny, npar = {}, 0.0, 0.0, 0, 0, 0
    for c, t in comps.items():
        A, b = _design(t, E, free_scale)
        beta, *_ = np.linalg.lstsq(A, b, rcond=None)
        r = b - A @ beta
        n = len(t)
        ssx += float(np.sum(r[:n] ** 2))
        ssy += float(np.sum(r[n:] ** 2))
        nx += n
        ny += n
        npar += A.shape[1]
        sols[c] = {"A": A, "b": b, "beta": beta, "resid": r, "n": n}
    # pooled noise per coordinate (degrees of freedom shared out in proportion)
    dof = nx + ny - npar
    s2x = ssx / (dof * nx / (nx + ny))
    s2y = ssy / (dof * ny / (nx + ny))
    out = {"sigma_x_units": float(np.sqrt(s2x)), "sigma_y_units": float(np.sqrt(s2y)), "dof": dof,
           "per_comparison": {}}
    for c, s in sols.items():
        n = s["n"]
        w = np.concatenate([np.full(n, 1 / s2x), np.full(n, 1 / s2y)])
        cov = np.linalg.pinv((s["A"] * w[:, None]).T @ s["A"])
        beta = cov @ (s["A"] * w[:, None]).T @ s["b"]
        out["per_comparison"][c] = {"lambda": float(beta[3]), "lambda_err": float(np.sqrt(cov[3, 3])),
                                    "theta": float(beta[2]), "theta_err": float(np.sqrt(cov[2, 2])),
                                    "limb": float(beta[3] * limb_unit), "limb_err": float(np.sqrt(cov[3, 3]) * limb_unit)}
    # combine: the two comparisons of one eclipse plate share its measuring errors
    # (eclipse and comparison plates assumed equally noisy: var(mean of 2) = 0.75 var(one))
    plates = {}
    for plate in sorted({PLATE_OF[c] for c in comparisons}):
        vals = [out["per_comparison"][c] for c in comparisons if PLATE_OF[c] == plate]
        lam = np.mean([v["lambda"] for v in vals])
        err = np.mean([v["lambda_err"] for v in vals]) * (np.sqrt(0.75) if len(vals) == 2 else 1.0)
        plates[plate] = {"lambda": float(lam), "lambda_err": float(err)}
    lam = np.mean([v["lambda"] for v in plates.values()])
    err = np.sqrt(np.sum([v["lambda_err"] ** 2 for v in plates.values()])) / len(plates)
    out["plates"] = plates
    out["lambda"], out["lambda_err"] = float(lam), float(err)
    out["limb"], out["limb_err_measures"] = float(lam * limb_unit), float(err * limb_unit)
    # the check-plate scale is uncertain (p.e. 0.14" at the limb, common to all): add it
    out["limb_err_scale"] = SCALE_PE_LIMB / 0.6745 if not free_scale else 0.0
    out["limb_err"] = float(np.hypot(out["limb_err_measures"], out["limb_err_scale"]))
    return out


if __name__ == "__main__":
    rep = replicate()
    print(rep["per_comparison"].round(3).to_string(index=False))
    print(f"mean {rep['mean']:.3f}  -> with -0.04: {rep['limb']:.3f}  (printed 1.65, 1.61)")
    for kw in ({}, {"pattern": "1920"}, {"free_scale": True}):
        f = fit(**kw)
        print(kw, f"limb {f['limb']:.3f} +- {f['limb_err']:.3f} (measures {f['limb_err_measures']:.3f}); "
                  f"sigma x {f['sigma_x_units'] * UNIT:.3f}\", y {f['sigma_y_units'] * UNIT:.3f}\"",
              {c: round(v['limb'], 2) for c, v in f["per_comparison"].items()})
