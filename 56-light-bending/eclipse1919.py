"""The 1919 eclipse plates, refitted from the measurements printed in 1920.

Dyson, Eddington & Davidson measured, for every star on every plate, the
small difference between its position on that plate and on a reference
plate, with a micrometer. Their equations of condition (p.304) are

    Dx = a x + b y + c + alpha Ex          Dy = d x + e y + f + alpha Ey

for a star at standard coordinates (x, y) (units of 50'): the six "plate
constants" absorb the shift (c, f), the change of scale (a, e) and the
rotation (b, d) between the two plates, and alpha is the deflection at 50'
from the Sun's centre, multiplying the pattern (Ex, Ey) of a 1/r outward
displacement. GR predicts alpha = 1.7512" x (Sun's radius / 50').

Three data sets:

* **Sobral 4-inch** (7 eclipse + 7 comparison plates, 7 stars, both
  coordinates; each plate measured against an intermediate "scale plate");
* **Sobral astrographic** (16 eclipse plates, declinations only; we have the
  five bright stars of p.311 and the per-plate solutions of Table IX);
* **Principe** (plates X and W, each against two Oxford plates, 5 stars, with
  the scale fixed from check plates of a field around Arcturus).

For each: (1) an exact replication of the 1920 reduction from the printed
measures, (2) modern least-squares refits with honest standard errors, and
the variants (scale free or fixed, plates or stars left out) that the
century-long argument about these plates turns on.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import config as C

PE = 0.6745   # probable error = 0.6745 sigma (the 1920 convention)

ECLIPSE_4IN = ["I", "II", "III", "IV", "V", "VII", "VIII"]
COMPARISON_4IN = ["14_2a", "14_2b", "15_1", "15_2", "17_1", "17_2", "18_2"]
STARS_4IN = [11, 5, 4, 3, 6, 10, 2]


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(C.TRANSCRIBED_DIR / name, comment="#", skipinitialspace=True)


def gr_alpha_at_50(sun_radius_arcmin: float) -> float:
    """GR's deflection (arcsec) at 50' from the Sun's centre."""
    return C.DEFLECTION_LIMB_GR * sun_radius_arcmin / C.UNIT_TABLE1_ARCMIN


def lsq(A: np.ndarray, y: np.ndarray, w: np.ndarray | None = None) -> dict:
    """Weighted linear least squares with the formal covariance (unit-weight scaled
    by the residuals unless weights are absolute)."""
    w = np.ones_like(y) if w is None else w
    sw = np.sqrt(w)
    beta, *_ = np.linalg.lstsq(A * sw[:, None], y * sw, rcond=None)
    r = y - A @ beta
    dof = len(y) - np.linalg.matrix_rank(A)
    s2 = float(np.sum(w * r**2) / dof) if dof > 0 else np.nan
    cov_unit = np.linalg.pinv((A * w[:, None]).T @ A)
    return {"beta": beta, "resid": r, "dof": dof, "s2": s2, "cov_unit": cov_unit, "cov": cov_unit * s2}


# ===========================================================================
# Sobral 4-inch
# ===========================================================================

#: Two numbers in the printed Table II are not the ones the 1920 reduction used.
#: The arithmetic downstream proves it twice over: every one of the 13 numbers of
#: Table III (or IV) for that plate, and the star's residual in Table VIII,
#: follow from the second value (transcribed/dyson1920_sobral4inch_notes.md, D1
#: and D2; both are most likely typesetting slips in Table II). The transcription
#: keeps the printed values; by default the analysis uses the 1920 working ones.
#: (The plate VIII Dy offset, printed -1.322 for -1.332 (D3), is a constant per
#: plate and drops out of every fit.)
TABLE2_AS_REDUCED_1920 = {("I", 2, "Dx_rev"): -0.798, ("14_2a", 2, "Dx_rev"): -0.346}


def load_4inch(as_printed: bool = False) -> dict:
    """Table II measures (revolutions) with the condition coefficients. By default the
    two Table II values the 1920 reduction actually used replace the misprinted ones;
    ``as_printed=True`` keeps the table exactly as printed."""
    coeff = _read("dyson1920_4inch_condition_coeffs.csv").set_index("star").loc[STARS_4IN]
    t2 = _read("dyson1920_table2_4inch.csv")
    if not as_printed:
        for (plate, star, col), v in TABLE2_AS_REDUCED_1920.items():
            sel = (t2["plate"].astype(str) == plate) & (t2["star"] == star)
            assert sel.sum() == 1, (plate, star)
            t2.loc[sel, col] = v
    off = _read("dyson1920_table2_offsets.csv").set_index("plate")
    Dx = t2.pivot(index="star", columns="plate", values="Dx_rev").loc[STARS_4IN]
    Dy = t2.pivot(index="star", columns="plate", values="Dy_rev").loc[STARS_4IN]
    plates = ECLIPSE_4IN + COMPARISON_4IN
    Dx, Dy = Dx[plates], Dy[plates]
    # the constants "taken out ... for arithmetical convenience" (p.303 footnote)
    Dx_used = Dx - off.loc[plates, "Dx_offset"].to_numpy()[None, :]
    Dy_used = Dy - off.loc[plates, "Dy_offset"].to_numpy()[None, :]
    return {"coeff": coeff, "Dx": Dx, "Dy": Dy, "Dx_used": Dx_used, "Dy_used": Dy_used, "offsets": off}


def _design_4inch(coeff: pd.DataFrame, coord: str) -> np.ndarray:
    x, y = coeff["x"].to_numpy(float), coeff["y"].to_numpy(float)
    one = np.ones_like(x)
    if coord == "RA":     # unknowns in the paper's order: c, b, a, alpha
        return np.column_stack([one, y, x, coeff["Ex"].to_numpy(float)])
    return np.column_stack([one, x, y, coeff["Ey"].to_numpy(float)])   # f, d, e, alpha


def _gauss_steps(N: np.ndarray, r: np.ndarray) -> list[float]:
    """Right-hand sides printed by Dyson et al.: the 4 normal equations, then the
    3, 2 and 1 equations left after successive Gaussian elimination."""
    N, r = N.astype(float).copy(), r.astype(float).copy()
    out = list(r)
    for k in range(3):
        f = N[k + 1:, k] / N[k, k]
        N[k + 1:, k + 1:] -= np.outer(f, N[k, k + 1:])
        r[k + 1:] -= f * r[k]
        out += list(r[k + 1:])
    return out


def replicate_4inch(d: dict | None = None) -> dict:
    """Tables III-VII from the Table II measures, and the 1".98 result."""
    d = load_4inch() if d is None else d
    rows = []
    for coord, D in (("RA", d["Dx_used"]), ("Dec", d["Dy_used"])):
        A = _design_4inch(d["coeff"], coord)
        for plate in D.columns:
            y = D[plate].to_numpy(float)
            N, rhs = A.T @ A, A.T @ y
            beta = np.linalg.solve(N, rhs)
            steps = _gauss_steps(N, rhs)
            rows.append({"coord": coord, "plate": plate,
                         "kind": "eclipse" if plate in ECLIPSE_4IN else "comparison",
                         "alpha": beta[3], "scale": beta[2], "orient": beta[1], "zero": beta[0],
                         **{f"rhs{i + 1}": v for i, v in enumerate(steps)},
                         "N": N})
    sol = pd.DataFrame(rows)
    out = {"solutions": sol}
    means = sol.groupby(["coord", "kind"])["alpha"].mean().unstack()
    out["table7_means"] = means
    a_ra = means.loc["RA", "eclipse"] - means.loc["RA", "comparison"]
    a_dec = means.loc["Dec", "eclipse"] - means.loc["Dec", "comparison"]
    # "the declination result is of double the weight of that from the right ascensions"
    a_rev = (a_ra + 2 * a_dec) / 3
    out.update({"alpha_ra_rev": a_ra, "alpha_dec_rev": a_dec, "alpha_rev": a_rev,
                "alpha_arcsec": a_rev * C.REV_4INCH,
                "limb_1920": a_rev * C.REV_4INCH * C.UNIT_TABLE1_ARCMIN / C.SUN_RADIUS_1919_ARCMIN,
                "limb_ra_1920": a_ra * C.REV_4INCH * C.UNIT_TABLE1_ARCMIN / C.SUN_RADIUS_1919_ARCMIN,
                "limb_dec_1920": a_dec * C.REV_4INCH * C.UNIT_TABLE1_ARCMIN / C.SUN_RADIUS_1919_ARCMIN})
    # the 1920 probable error "judging from the accordance of the separate determinations" is
    # "about 6 per cent" (p.306; method not stated). Ours, from the scatter of the eclipse and of
    # the comparison plates' alphas in each coordinate (RA 7.9 %, Dec 6.0 %), combined 1:2 like
    # the result: 4.8 %.
    se = {}
    for coord in ("RA", "Dec"):
        g = sol[sol["coord"] == coord]
        e, c = g.loc[g["kind"] == "eclipse", "alpha"], g.loc[g["kind"] == "comparison", "alpha"]
        se[coord] = np.sqrt(e.var(ddof=1) / len(e) + c.var(ddof=1) / len(c))
    pe = PE * np.sqrt(se["RA"] ** 2 + 4 * se["Dec"] ** 2) / 3
    out["pe_fraction_1920"] = float(pe / a_rev)
    out["pe_fraction_by_coord"] = {"RA": float(PE * se["RA"] / a_ra), "Dec": float(PE * se["Dec"] / a_dec)}
    return out


def compare_printed_4inch(rep: dict) -> pd.DataFrame:
    """Our Tables III-VI against the printed ones (rhs1..rhs10 and solutions)."""
    pr = _read("dyson1920_tables3to6_normal_equations.csv")
    table = {("RA", "eclipse"): "III", ("RA", "comparison"): "IV",
             ("Dec", "eclipse"): "V", ("Dec", "comparison"): "VI"}
    rows = []
    for _, s in rep["solutions"].iterrows():
        t = table[(s["coord"], s["kind"])]
        p = pr[(pr["table"] == t) & (pr["plate"] == s["plate"])].set_index("quantity")["value"]
        for q in [f"rhs{i}" for i in range(1, 11)]:
            if q in p.index:
                rows.append({"table": t, "plate": s["plate"], "quantity": q, "printed": p[q], "ours": s[q]})
        sol_names = {"sol_alpha": "alpha", "sol_a": "scale", "sol_b": "orient", "sol_e": "scale", "sol_d": "orient"}
        for q, col in sol_names.items():
            if q in p.index:
                rows.append({"table": t, "plate": s["plate"], "quantity": q, "printed": p[q], "ours": s[col]})
    df = pd.DataFrame(rows)
    df["diff"] = df["ours"] - df["printed"]
    return df


def joint_fit_4inch(d: dict | None = None, coords=("RA", "Dec"), drop_stars=(), drop_plates=(),
                    alpha_free: bool = True, n_iter: int = 20) -> dict:
    """Modern least squares over all 14 plates at once.

    D[p, i] = c_p + a_p x_i + b_p y_i + s_i + alpha E_i [eclipse plates only] + noise

    s_i is the error of the intermediate scale plate at star i, common to every
    plate measured against it (made orthogonal to 1, x, y, which the plate
    constants already absorb). The noise variance is estimated separately for
    eclipse and comparison plates in each coordinate (iteratively reweighted,
    with degrees of freedom counted per group). Returns alpha in arcsec at 50'.
    """
    d = load_4inch() if d is None else d
    stars = [s for s in STARS_4IN if s not in drop_stars]
    plates = [p for p in ECLIPSE_4IN + COMPARISON_4IN if p not in drop_plates]
    coeff = d["coeff"].loc[stars]
    x, y = coeff["x"].to_numpy(float), coeff["y"].to_numpy(float)
    ns, npl = len(stars), len(plates)
    # orthonormal basis of star-patterns orthogonal to (1, x, y)
    Q, _ = np.linalg.qr(np.column_stack([np.ones(ns), x, y, np.eye(ns)]))
    P = Q[:, 3:ns]                                    # ns x (ns-3)
    blocks, yv, groups = [], [], []
    ncols_per_coord = 3 * npl + (ns - 3)
    ncol = ncols_per_coord * len(coords) + (1 if alpha_free else 0)
    for ci, coord in enumerate(coords):
        D = (d["Dx_used"] if coord == "RA" else d["Dy_used"]).loc[stars, plates] * C.REV_4INCH
        E = coeff["Ex" if coord == "RA" else "Ey"].to_numpy(float)
        off = ci * ncols_per_coord
        for pj, plate in enumerate(plates):
            A = np.zeros((ns, ncol))
            A[:, off + 3 * pj] = 1.0
            A[:, off + 3 * pj + 1] = x
            A[:, off + 3 * pj + 2] = y
            A[:, off + 3 * npl: off + 3 * npl + ns - 3] = P
            ecl = plate in ECLIPSE_4IN
            if ecl and alpha_free:
                A[:, -1] = E
            blocks.append(A)
            yv.append(D[plate].to_numpy(float))
            groups += [f"{coord}-{'ecl' if ecl else 'cmp'}"] * ns
    A, yv, groups = np.vstack(blocks), np.concatenate(yv), np.array(groups)
    gnames = sorted(set(groups))
    var = {g: 0.1**2 for g in gnames}
    for _ in range(n_iter):
        w = np.array([1.0 / var[g] for g in groups])
        fit = lsq(A, yv, w)
        # per-group variance with an approximate dof share (leverage-based)
        H = A @ fit["cov_unit"] @ (A * w[:, None]).T
        lev = np.diag(H)
        new = {g: float(np.sum(fit["resid"][groups == g] ** 2) / np.sum(1 - lev[groups == g])) for g in gnames}
        if all(abs(new[g] / var[g] - 1) < 1e-6 for g in gnames):
            var = new
            break
        var = new
    w = np.array([1.0 / var[g] for g in groups])
    fit = lsq(A, yv, w)
    chi2 = float(np.sum(w * fit["resid"] ** 2))
    out = {"sigma_by_group": {g: float(np.sqrt(v)) for g, v in var.items()}, "chi2": chi2, "dof": fit["dof"],
           "n_data": len(yv), "fit": fit, "groups": groups, "plates": plates, "stars": stars}
    if alpha_free:
        out["alpha"] = float(fit["beta"][-1])
        out["alpha_err"] = float(np.sqrt(fit["cov_unit"][-1, -1]))   # absolute weights: no rescaling
    return out


def jackknife_4inch(d: dict | None = None, coords=("RA", "Dec")) -> pd.DataFrame:
    """Leave one eclipse plate, comparison plate or star out at a time."""
    d = load_4inch() if d is None else d
    rows = [{"left_out": "none", **_alpha_row(joint_fit_4inch(d, coords))}]
    for p in ECLIPSE_4IN + COMPARISON_4IN:
        rows.append({"left_out": f"plate {p}", **_alpha_row(joint_fit_4inch(d, coords, drop_plates=(p,)))})
    for s in STARS_4IN:
        rows.append({"left_out": f"star {s}", **_alpha_row(joint_fit_4inch(d, coords, drop_stars=(s,)))})
    return pd.DataFrame(rows)


def _alpha_row(f: dict) -> dict:
    return {"alpha": f["alpha"], "alpha_err": f["alpha_err"], "chi2": f["chi2"], "dof": f["dof"]}


def to_limb(alpha_arcsec: float, sun_radius_arcmin: float) -> float:
    """Deflection at 50' -> at the limb (1/r law), and as a fraction of GR."""
    return alpha_arcsec * C.UNIT_TABLE1_ARCMIN / sun_radius_arcmin


# ===========================================================================
# Sobral astrographic (declinations only)
# ===========================================================================

STARS_ASTRO = [5, 4, 3, 6, 10]


def load_astro() -> dict:
    dy = _read("dyson1920_p311_astrographic_dy.csv")
    dy["eclipse_plate"] = dy["eclipse_plate"].astype(str)
    Dy = dy.pivot(index="eclipse_plate", columns="star", values="Dy_rev")[STARS_ASTRO]
    t9 = _read("dyson1920_table9_astrographic.csv")
    t9["eclipse_plate"] = t9["eclipse_plate"].astype(str)
    t9_mean = t9[t9["eclipse_plate"] == "Mean"]
    t9 = t9[t9["eclipse_plate"] != "Mean"].reset_index(drop=True)      # the printed means are a separate row
    coeff = _read("dyson1920_4inch_condition_coeffs.csv").set_index("star").loc[STARS_ASTRO]
    order = sorted(Dy.index, key=lambda s: int(s))
    return {"Dy": Dy.loc[order], "t9": t9, "t9_printed_mean": t9_mean, "coeff": coeff}


def replicate_astro(d: dict | None = None) -> dict:
    """p.311-312: the combination alpha + 0.240 e per plate, its mean, and
    the deflection for the two choices of the scale e (0.082 r, 0.020 r)."""
    d = load_astro() if d is None else d
    Dy, t9 = d["Dy"], d["t9"]
    # Table IX lists plates 7 and 12 twice (two comparison plates) but p.311 has one
    # row for each, and the paper does not say which d goes with it. We use the mean
    # of the two (the arithmetic of the printed values favours 18_3's d for plate 7
    # and the mean for plate 12; the choice moves the result by < 0.01").
    d_of = t9.groupby("eclipse_plate")["d_rev"].mean()
    rows = []
    for p, r in Dy.iterrows():
        if np.isfinite(r[10]):
            rhs = r[3] + r[4] - r[5] - 0.5 * (r[6] + r[10])
            k_e, k_d, k_a = 0.729, 0.118, 3.041
        else:   # plate 15: star 10 missing; eliminate f between (1), (4) and (2)+(3) instead
            rhs, k_a, k_e, k_d = _astro_combination_without_10(d["coeff"], r)
        val = (rhs - k_d * d_of.get(p, 0.0)) / k_a
        rows.append({"plate": p, "coef_e": k_e / k_a, "alpha_plus_ke": val})
    per = pd.DataFrame(rows)
    mean_coef = per["coef_e"].mean()
    mean_val = per["alpha_plus_ke"].mean()
    to_limb_rev = C.REV_ASTRO * C.UNIT_TABLE1_ARCMIN / C.SUN_RADIUS_1919_ARCMIN
    out = {"per_plate": per, "mean_coef_e": mean_coef, "mean_value_rev": mean_val}
    for label, e in (("plates_scale", 0.082), ("refraction_scale", 0.020)):
        a = mean_val - mean_coef * e
        out[f"alpha_rev_{label}"] = a
        out[f"limb_{label}"] = a * to_limb_rev
    out["alpha_vs_e"] = lambda e: (mean_val - mean_coef * np.asarray(e)) * to_limb_rev
    # per-plate least squares (f, d, e, alpha free), p.312: "mean result 0".99"
    A = np.column_stack([np.ones(5), d["coeff"]["x"], d["coeff"]["y"], d["coeff"]["Ey"]])
    a_pl = []
    for p, r in Dy.iterrows():
        ok = np.isfinite(r.to_numpy(float))
        if ok.sum() >= 5:
            a_pl.append(np.linalg.lstsq(A[ok], r.to_numpy(float)[ok], rcond=None)[0][3])
        else:   # 4 stars, 4 unknowns: exactly determined
            a_pl.append(np.linalg.solve(A[ok], r.to_numpy(float)[ok])[3])
    out["per_plate_lsq_alpha_rev"] = np.array(a_pl)
    out["per_plate_lsq_limb_mean"] = float(np.mean(a_pl) * to_limb_rev)
    t9u = t9.copy()
    out["table9_limb_mean"] = float(t9u["alpha_limb_arcsec"].mean())
    return out


def _astro_combination_without_10(coeff: pd.DataFrame, r: pd.Series):
    """For plate 15 (star 10 not shown): (2)+(3) minus (1) minus (4), the analogue
    of the p.311 combination, with its coefficients of alpha, e and d."""
    w = {3: 1.0, 4: 1.0, 5: -1.0, 6: -1.0}
    rhs = sum(w[s] * r[s] for s in w)
    k_a = sum(w[s] * coeff.loc[s, "Ey"] for s in w)
    k_e = sum(w[s] * coeff.loc[s, "y"] for s in w)
    k_d = sum(w[s] * coeff.loc[s, "x"] for s in w)
    return rhs, k_a, k_e, k_d


def fit_astro(d: dict | None = None, scale: str = "common", e_fixed: float | None = None) -> dict:
    """Joint least squares over the 16 plates (five bright stars, declinations):

        Dy[p, i] = f_p + d_p x_i + e y_i + alpha Ey_i

    scale = 'per_plate' (e_p free on each plate), 'common' (one e for all), or
    'fixed' (e = e_fixed, in revolutions). alpha in arcsec at 50'."""
    d = load_astro() if d is None else d
    Dy = d["Dy"]
    x, y, Ey = (d["coeff"][k].to_numpy(float) for k in ("x", "y", "Ey"))
    plates = list(Dy.index)
    npl = len(plates)
    n_e = npl if scale == "per_plate" else (1 if scale == "common" else 0)
    ncol = 2 * npl + n_e + 1
    rows, yv = [], []
    for j, p in enumerate(plates):
        vals = Dy.loc[p].to_numpy(float)
        for i in range(5):
            if not np.isfinite(vals[i]):
                continue
            a = np.zeros(ncol)
            a[2 * j], a[2 * j + 1] = 1.0, x[i]
            if scale == "per_plate":
                a[2 * npl + j] = y[i]
            elif scale == "common":
                a[2 * npl] = y[i]
            a[-1] = Ey[i]
            rows.append(a)
            yv.append(vals[i] - (e_fixed * y[i] if scale == "fixed" else 0.0))
    A, yv = np.array(rows), np.array(yv) * C.REV_ASTRO
    fit = lsq(A, yv)
    out = {"alpha": float(fit["beta"][-1]), "alpha_err": float(np.sqrt(fit["cov"][-1, -1])),
           "sigma": float(np.sqrt(fit["s2"])), "dof": fit["dof"], "fit": fit}
    if scale == "common":
        out["e_rev"] = float(fit["beta"][2 * npl] / C.REV_ASTRO)
        out["e_err_rev"] = float(np.sqrt(fit["cov"][2 * npl, 2 * npl]) / C.REV_ASTRO)
        out["corr_alpha_e"] = float(fit["cov"][-1, 2 * npl] / np.sqrt(fit["cov"][-1, -1] * fit["cov"][2 * npl, 2 * npl]))
    return out
