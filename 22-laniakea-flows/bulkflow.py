"""Bulk flows: is the local universe streaming too fast for Lambda-CDM?

The method follows Nusser & Tully (2026) exactly, so their numbers can be
reproduced, and then extends it to "spheres of growing radius".

1. Every CF4 galaxy is *placed* at its CMB-frame redshift distance r (h^-1 Mpc),
   and its measured distance is turned into a radial velocity with the linear
   log-distance estimator u = -K a (DM - DM_z), a = ln(10)/5 (K ~ cz nearby).
   The velocity error is K a e_DM, plus a 100 km/s floor in quadrature.
2. In a radial bin, a weighted least-squares fit of u = m + V . rhat gives a
   monopole m (an overall expansion or contraction, which also absorbs any
   error in the distance zero point) and a dipole V, the bulk flow, in
   supergalactic components.
3. What should Lambda-CDM give *for this exact set of galaxies and weights*?
   Linear theory predicts the covariance of the radial velocities of any two
   galaxies from the matter power spectrum (Gorski 1989). Pushing it through
   the same fit gives the covariance of (m, V) expected from cosmic variance;
   adding the measurement covariance gives the predictive covariance, and a
   chi-squared of the measured (m, V) against it gives a p-value.

Step 3 involves every pair of ~51,500 galaxies (2.7 billion pairs). Rather
than a 20 GB matrix, the code accumulates, in chunks, the covariance of each
thin shell's weighted sums b_s = X_s^T W_s u_s; the fit in any bin or sphere
is (sum F_s)^-1 sum b_s, so one pass serves every bin and every sphere.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
from scipy import stats

import catalogs
import config
import cosmology as cg
from config import C_KMS

#: The Hubble constant behind DM_z, the distance modulus a galaxy would have
#: with no peculiar velocity. It must match CF4's distance scale (H0 ~ 75):
#: any mismatch becomes a spurious monopole. Nusser & Tully do not quote the
#: value; 75 (CF4's own convention) reproduces their monopoles (see README).
DMZ_H0 = 75.0


#: Which CF4 table: "corrected" (SDSS sign fixed, ``cf4fix``) or "published" (as on VizieR).
VERSION = "corrected"


def sample(which: str = "all", H0: float = DMZ_H0, version: str | None = None) -> pd.DataFrame:
    """CF4 individual galaxies prepared for the fits.

    ``which``: "all" (the combined distance of every galaxy, N&T's "All
    Individual"), or one method's own distances: "TF", "FP", "SN" (SN Ia and
    SN II combined by inverse variance), "FP north" / "FP south" (the SDSS and
    6dFGS Fundamental Plane samples, separated on the sky).
    """
    version = version or VERSION
    if version == "corrected":
        import cf4fix

        gal = cf4fix.corrected_galaxies()
    else:
        gal = catalogs.cf4_galaxies()
    df = pd.DataFrame({"PGC": gal["PGC"], "z": gal["z_cmb"], "zh": gal["z_helio"], "DEdeg": gal["DEdeg"],
                       "GLAT": gal["GLAT"], "sgx": gal["sgx"], "sgy": gal["sgy"], "sgz": gal["sgz"]})
    if which == "all":
        dm, e = gal["DM"], gal["e_DM"]
    elif which in ("TF", "FP", "FP north", "FP south"):
        col = {"TF": "DMtf"}.get(which, "DMfp")
        dm, e = gal[col], gal["e_" + col]
    elif which == "SN":
        w1 = 1 / gal["e_DMsnIa"] ** 2
        w2 = 1 / gal["e_DMsnII"] ** 2
        num = (gal["DMsnIa"] * w1).fillna(0) + (gal["DMsnII"] * w2).fillna(0)
        den = w1.fillna(0) + w2.fillna(0)
        dm = (num / den).where(den > 0)
        e = (1 / np.sqrt(den)).where(den > 0)
    else:
        raise ValueError(which)
    df["DM"], df["e_DM"] = dm.astype(float), e.astype(float)
    if which == "FP north":
        df = df[(df["DEdeg"] > -3.0) & (df["GLAT"] > 0)]
    elif which == "FP south":
        df = df[df["DEdeg"] < 0.0]
    df = df[df["DM"].notna() & (df["e_DM"] > 0) & (df["z"] > 0)].copy()
    df["r"] = cg.comoving_distance(df["z"], config.PLANCK18.Om0)
    df["cz"] = df["z"] * C_KMS
    df["K"] = cg.nt_K(df["z"], df["zh"])
    df["DMz"] = cg.dm_from_redshift(df["z"], df["zh"], H0)
    df["dmu"] = df["DM"] - df["DMz"]
    df["u"] = -df["K"] * cg.A_LN10_5 * df["dmu"]
    df["sig_u"] = np.sqrt((df["K"] * cg.A_LN10_5 * df["e_DM"]) ** 2 + config.BULK_SIGMA_FLOOR**2)
    return df.reset_index(drop=True)


def fitted_objects(df: pd.DataFrame) -> pd.DataFrame:
    """The objects that enter any fit: beyond cz_CMB = 4,000 km/s and inside 300 h^-1 Mpc."""
    return df[(df["cz"] > config.BULK_CZ_MIN) & (df["r"] < config.SPHERE_RMAX)].reset_index(drop=True)


# --------------------------------------------------------------------------
# Shell sums and their Lambda-CDM covariance
# --------------------------------------------------------------------------

def shell_index(r: np.ndarray) -> np.ndarray:
    return np.floor(np.asarray(r) / config.SPHERE_SHELL).astype(int)


def shell_sums(df: pd.DataFrame) -> dict:
    """Per thin shell: b_s = X^T W u (4) and F_s = X^T W X (4x4), X = [1, rhat_SG]."""
    nsh = int(round(config.SPHERE_RMAX / config.SPHERE_SHELL))
    X = np.column_stack([np.ones(len(df)), df[["sgx", "sgy", "sgz"]].to_numpy()])
    w = 1.0 / df["sig_u"].to_numpy() ** 2
    u = df["u"].to_numpy()
    sh = shell_index(df["r"].to_numpy())
    b = np.zeros((nsh, 4))
    F = np.zeros((nsh, 4, 4))
    n = np.zeros(nsh, int)
    sw = np.zeros(nsh)
    sw2 = np.zeros(nsh)
    for s in range(nsh):
        m = sh == s
        Xs, ws = X[m], w[m]
        b[s] = Xs.T @ (ws * u[m])
        F[s] = (Xs * ws[:, None]).T @ Xs
        n[s] = m.sum()
        sw[s], sw2[s] = ws.sum(), (ws**2).sum()
    return {"b": b, "F": F, "n": n, "sw": sw, "sw2": sw2, "X": X, "w": w, "shell": sh}


def lcdm_shell_covariance(df: pd.DataFrame, *, chunk: int = 256, cache_key: str | None = None,
                          verbose: bool = True) -> np.ndarray:
    """Cosmic-variance covariance of the shell sums: M[s,a ; t,b] = sum_ij G_i,sa S_ij G_j,tb,
    with G_i = w_i X_i in the columns of object i's shell and S_ij the linear-theory
    covariance of the radial velocities of objects i and j (Gorski 1989).

    Returns an array of shape (nsh*4, nsh*4). Cached on disk by ``cache_key``.
    """
    nsh = int(round(config.SPHERE_RMAX / config.SPHERE_SHELL))
    path = config.DERIVED_DIR / f"lcdm_cov_{cache_key}.npz" if cache_key else None
    if path is not None and path.exists():
        d = np.load(path)
        if d["n"] == len(df):
            return d["M"]
    s_tab, par, perp = cg.velocity_correlations()
    ds = s_tab[1] - s_tab[0]
    par32, perp32 = par.astype(np.float32), perp.astype(np.float32)
    r = df["r"].to_numpy()
    rhat = df[["sgx", "sgy", "sgz"]].to_numpy()
    sums = shell_sums(df)
    X, w, sh = sums["X"], sums["w"], sums["shell"]
    G = np.zeros((len(df), nsh * 4))
    rows = np.arange(len(df))
    for a in range(4):
        G[rows, sh * 4 + a] = w * X[:, a]
    M = np.zeros((nsh * 4, nsh * 4))
    r32 = r.astype(np.float32)
    rhat32 = rhat.astype(np.float32)
    t0 = time.time()
    for i0 in range(0, len(df), chunk):
        i1 = min(i0 + chunk, len(df))
        cos = rhat32[i0:i1] @ rhat32.T
        ri = r32[i0:i1, None]
        rj = r32[None, :]
        s2 = np.maximum(ri * ri + rj * rj - 2.0 * ri * rj * cos, 0.0)
        s = np.sqrt(s2)
        x = s / ds
        k = np.minimum(x.astype(np.int32), len(s_tab) - 2)
        fr = x - k
        p_par = par32[k] * (1 - fr) + par32[k + 1] * fr
        p_perp = perp32[k] * (1 - fr) + perp32[k + 1] * fr
        with np.errstate(invalid="ignore", divide="ignore"):
            proj = np.where(s2 > 1e-8, (ri - rj * cos) * (ri * cos - rj) / s2, 0.0)
        S = (p_perp * cos + (p_par - p_perp) * proj).astype(np.float64)
        M += G[i0:i1].T @ (S @ G)
        if verbose and (i0 // chunk) % 40 == 0:
            print(f"    covariance rows {i1:,}/{len(df):,}  ({time.time() - t0:.0f} s)", flush=True)
    M = 0.5 * (M + M.T)
    if path is not None:
        np.savez(path, M=M, n=len(df))
    return M


# --------------------------------------------------------------------------
# Fits in bins and spheres
# --------------------------------------------------------------------------

PARAMS = {"joint": [0, 1, 2, 3], "dipole": [1, 2, 3], "monopole": [0]}


def combine(sums: dict, M: np.ndarray, shells: list[int], model: str = "joint") -> dict:
    """The fit over a set of thin shells, with measurement and Lambda-CDM covariances."""
    idx = PARAMS[model]
    F = sums["F"][shells].sum(axis=0)[np.ix_(idx, idx)]
    b = sums["b"][shells].sum(axis=0)[idx]
    Finv = np.linalg.inv(F)
    beta = Finv @ b
    cols = np.concatenate([np.array(idx) + 4 * s for s in shells])
    Msum = M[np.ix_(cols, cols)].reshape(len(shells), len(idx), len(shells), len(idx)).sum(axis=(0, 2))
    C_cv = Finv @ Msum @ Finv
    C_meas = Finv
    return {"beta": beta, "C_meas": C_meas, "C_cv": C_cv, "C_pred": C_meas + C_cv, "idx": idx}


def cross_cov(sums: dict, M: np.ndarray, shells_a: list[int], shells_b: list[int], model: str = "joint"):
    """Cosmic covariance between the fitted parameters of two disjoint sets of shells."""
    idx = PARAMS[model]
    Fa = np.linalg.inv(sums["F"][shells_a].sum(axis=0)[np.ix_(idx, idx)])
    Fb = np.linalg.inv(sums["F"][shells_b].sum(axis=0)[np.ix_(idx, idx)])
    ca = np.concatenate([np.array(idx) + 4 * s for s in shells_a])
    cb = np.concatenate([np.array(idx) + 4 * s for s in shells_b])
    Mab = M[np.ix_(ca, cb)].reshape(len(shells_a), len(idx), len(shells_b), len(idx)).sum(axis=(0, 2))
    return Fa @ Mab @ Fb.T


def pvalue(theta: np.ndarray, C: np.ndarray) -> tuple[float, float]:
    """Chi-squared p-value of theta against zero-mean covariance C, and its two-sided Gaussian sigma."""
    Q = float(theta @ np.linalg.solve(C, theta))
    p = float(stats.chi2.sf(Q, len(theta)))
    return p, float(stats.norm.isf(p / 2))


def _shells_for(lo: float, hi: float) -> list[int]:
    return list(range(int(round(lo / config.SPHERE_SHELL)), int(round(hi / config.SPHERE_SHELL))))


def describe(fit: dict, sums: dict, shells: list[int]) -> dict:
    """Summary numbers for a joint fit: (m, V) with errors, |V| and its ΛCDM tests."""
    beta, Cp, Cm = fit["beta"], fit["C_pred"], fit["C_meas"]
    V = beta[1:4]
    Vamp = float(np.linalg.norm(V))
    grad = V / Vamp
    sw, sw2 = sums["sw"][shells].sum(), sums["sw2"][shells].sum()
    l_sg, b_sg = catalogs.lonlat(V)
    l_g, b_g = catalogs.lonlat(catalogs.sg_to_gal(V))
    pm = pvalue(beta[:1], Cp[:1, :1])
    pV = pvalue(V, Cp[1:4, 1:4])
    pmV = pvalue(beta, Cp)
    return {
        "n": int(sums["n"][shells].sum()), "neff": float(sw**2 / sw2) if sw2 > 0 else 0.0,
        "m": float(beta[0]), "m_err_meas": float(np.sqrt(Cm[0, 0])), "m_err_pred": float(np.sqrt(Cp[0, 0])),
        "V_sg": V.tolist(), "V_err_meas": np.sqrt(np.diag(Cm)[1:4]).tolist(),
        "V_err_pred": np.sqrt(np.diag(Cp)[1:4]).tolist(),
        "Vamp": Vamp, "Vamp_err_meas": float(np.sqrt(grad @ Cm[1:4, 1:4] @ grad)),
        "dir_sg": [float(l_sg), float(b_sg)], "dir_gal": [float(l_g), float(b_g)],
        "sigma_dev": (beta / np.sqrt(np.diag(Cp))).tolist(),
        "pm": pm[0], "pm_sigma": pm[1], "pV": pV[0], "pV_sigma": pV[1], "pmV": pmV[0], "pmV_sigma": pmV[1],
        "lcdm_rms_V": float(np.sqrt(np.trace(fit["C_cv"][1:4, 1:4]))),
        "pred_rms_V": float(np.sqrt(np.trace(Cp[1:4, 1:4]))),
    }


def vamp_percentiles(C: np.ndarray, q=(16, 50, 84, 97.7), n: int = 200_000, seed: int = 1) -> list[float]:
    """Percentiles of |V| for V ~ N(0, C): the Lambda-CDM band for a window-matched dipole."""
    rng = np.random.default_rng(seed)
    v = rng.multivariate_normal(np.zeros(3), C, size=n)
    return np.percentile(np.linalg.norm(v, axis=1), q).tolist()


def analyse(which: str = "all", *, H0: float = DMZ_H0, version: str | None = None, verbose: bool = True) -> dict:
    """Bins (Nusser & Tully's) and cumulative spheres for one sample.

    The Lambda-CDM covariance depends only on positions and weights, which the
    SDSS correction does not change, so both catalogue versions share it.
    """
    version = version or VERSION
    df = fitted_objects(sample(which, H0, version))
    key = f"{which.replace(' ', '_')}_H{H0:g}"
    M = lcdm_shell_covariance(df, cache_key=key, verbose=verbose)
    sums = shell_sums(df)
    res = {"sample": which, "version": version, "H0_dmz": H0, "n_fitted": int(len(df)), "bins": {}, "spheres": [],
           "global": {}}
    used = []
    for lo, hi in config.BULK_BINS:
        shells = _shells_for(lo, hi)
        sw, sw2 = sums["sw"][shells].sum(), sums["sw2"][shells].sum()
        n = int(sums["n"][shells].sum())
        if n < 10 or sw2 == 0 or sw**2 / sw2 < config.BULK_NEFF_MIN:
            continue
        fit = combine(sums, M, shells)
        d = describe(fit, sums, shells)
        d["reff"] = float(np.sum(sums["w"][np.isin(sums["shell"], shells)] * df["r"].to_numpy()[np.isin(sums["shell"], shells)]) / sw)
        d["lcdm_band_V"] = vamp_percentiles(fit["C_pred"][1:4, 1:4])
        res["bins"][f"{lo:g}-{hi:g}"] = d
        used.append(shells)
    # Global tests across bins, with the full cross-bin covariance (N&T Table 4).
    for label, upto in (("40-160", 160.0), ("40-300", 300.0)):
        sel = [s for s in used if s[-1] * config.SPHERE_SHELL < upto]
        if len(sel) < 2:
            continue
        thetas, blocks = [], []
        for a in sel:
            thetas.append(combine(sums, M, a)["beta"])
        nb = len(sel)
        C = np.zeros((4 * nb, 4 * nb))
        for i, a in enumerate(sel):
            C[4 * i:4 * i + 4, 4 * i:4 * i + 4] = combine(sums, M, a)["C_pred"]
            for j, bsh in enumerate(sel):
                if j > i:
                    X = cross_cov(sums, M, a, bsh)
                    C[4 * i:4 * i + 4, 4 * j:4 * j + 4] = X
                    C[4 * j:4 * j + 4, 4 * i:4 * i + 4] = X.T
        th = np.concatenate(thetas)
        im = np.arange(0, 4 * nb, 4)
        iv = np.array([k for k in range(4 * nb) if k % 4])
        res["global"][label] = {
            "nbins": nb, "nobj": int(sum(sums["n"][s].sum() for s in sel)),
            "pm": pvalue(th[im], C[np.ix_(im, im)])[0], "pV": pvalue(th[iv], C[np.ix_(iv, iv)])[0],
            "pmV": pvalue(th, C)[0],
        }
        for k in ("pm", "pV", "pmV"):
            res["global"][label][k + "_sigma"] = float(stats.norm.isf(res["global"][label][k] / 2))
    # Cumulative spheres: every fitted object out to R.
    first = int(np.argmax(sums["n"] > 0))
    for R in np.arange(config.SPHERE_SHELL * (first + 2), config.SPHERE_RMAX + 1e-6, config.SPHERE_SHELL):
        shells = list(range(first, int(round(R / config.SPHERE_SHELL))))
        if sums["n"][shells].sum() < 50:
            continue
        fit = combine(sums, M, shells)
        d = describe(fit, sums, shells)
        d["R"] = float(R)
        mask = np.isin(sums["shell"], shells)
        d["reff"] = float(np.sum(sums["w"][mask] * df["r"].to_numpy()[mask]) / sums["sw"][shells].sum())
        d["lcdm_band_V"] = vamp_percentiles(fit["C_pred"][1:4, 1:4], n=50_000)
        res["spheres"].append(d)
    return res


def ideal_sphere_expectation() -> dict:
    """Lambda-CDM rms bulk flow of ideal, fully sampled top-hat spheres (no noise, no window)."""
    R = np.arange(10.0, 301.0, 5.0)
    return {"R": R.tolist(), "rms_tophat": cg.bulk_flow_rms(R, "tophat").tolist(),
            "rms_gaussian": cg.bulk_flow_rms(R, "gaussian").tolist()}
