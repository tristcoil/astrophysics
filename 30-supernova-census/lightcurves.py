"""Type II plateaus: from ZTF light curves to explosion energies and ejecta masses.

A Type II-P supernova is a red supergiant's hydrogen envelope thrown off at
several thousand km/s. As it expands and cools to ~6,000 K the hydrogen
recombines, and a recombination front eats inward through the envelope. While
it does, the light stays nearly constant: the plateau. When the front reaches
the bottom of the envelope the light drops by 1-3 mag within a few weeks onto
the radioactive tail, powered by cobalt-56 decay. So

* the plateau's **duration** grows with the envelope mass (more to recombine)
  and falls with the explosion energy (faster expansion, faster cooling);
* its **brightness** grows with energy and with the star's radius;
* the **tail** weighs the nickel-56.

Popov (1993) worked out the analytic scalings; Goldberg, Bildsten & Paxton
(2019) recalibrated them on 57 modern explosion models (MESA + STELLA):

    log L50  = 42.16 - 0.40 log M10 + 0.74 log E51 + 0.76 log R500
    log t_p  = 2.184 + 0.134 log M_Ni + 0.411 log M10 - 0.282 log E51

(L50 in erg/s at day 50; t_p in days; M10 = M_ej / 10 M_sun; E51 = E / 1e51
erg; R500 = R / 500 R_sun; M_Ni in M_sun). Two observables cannot fix three
unknowns: E and M_ej come out as functions of the progenitor radius R (their
Eq. 21), a degeneracy we show rather than hide.

Measurements here: explosion time from the last non-detection; plateau end
t_PT from a Fermi-function fit to the r-band fall (Valenti et al. 2016);
bolometric L50 from g and r with Lyman et al.'s (2014) Type II bolometric
correction; nickel from the tail with full gamma-ray trapping.
"""

from __future__ import annotations

import json
import warnings

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

import config as C

L_SUN = 3.828e33  # erg/s (IAU nominal)
# Nadyozhin (1994): energy deposition of 1 M_sun of 56Ni and its 56Co daughter
NI_E1, NI_TAU1 = 6.45e43, 8.8      # erg/s, days (56Ni)
NI_E2, NI_TAU2 = 1.45e43, 111.3    # erg/s, days (56Co)


# ---------------------------------------------------------------------------
# Reading ALeRCE light curves
# ---------------------------------------------------------------------------

def load_lc(oid: str) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    path = C.ALERCE_DIR / f"{oid}.json"
    if not path.exists():
        return None
    data = json.loads(path.read_text())
    if data.get("missing"):
        return None
    det = pd.DataFrame(data["detections"])
    non = pd.DataFrame(data["non_detections"])
    if det.empty:
        return None
    # positive subtractions only (the transient got brighter than the
    # reference), and not flagged dubious by the broker
    det = det[(det["isdiffpos"].astype(float) > 0) & ~det["dubious"].fillna(False).astype(bool)]
    det = det[np.isfinite(det["magpsf"]) & (det["fid"].isin([1, 2]))]
    det = det.sort_values("mjd").drop_duplicates(["mjd", "fid"]).reset_index(drop=True)
    if det.empty:
        return None
    if non.empty:
        non = pd.DataFrame(columns=["mjd", "fid", "diffmaglim"])
    non = non[non["fid"].isin([1, 2])].sort_values("mjd").reset_index(drop=True) if len(non) else non
    return det, non


def nightly(det: pd.DataFrame) -> pd.DataFrame:
    """Inverse-variance mean per night and filter (several exposures a night)."""
    w = 1.0 / np.maximum(det["sigmapsf"].to_numpy(), 0.02) ** 2
    d = pd.DataFrame({"night": np.floor(det["mjd"].to_numpy() - 0.3), "fid": det["fid"].to_numpy(),
                      "w": w, "wt": w * det["mjd"].to_numpy(), "wm": w * det["magpsf"].to_numpy()})
    g = d.groupby(["night", "fid"]).agg(w=("w", "sum"), wt=("wt", "sum"), wm=("wm", "sum"),
                                         n=("w", "size")).reset_index()
    out = pd.DataFrame({"night": g["night"], "fid": g["fid"], "mjd": g["wt"] / g["w"],
                        "mag": g["wm"] / g["w"], "err": g["w"] ** -0.5, "n": g["n"]})
    return out.sort_values("mjd").reset_index(drop=True)


def explosion_epoch(det: pd.DataFrame, non: pd.DataFrame, z: float) -> dict:
    """Midpoint of the last non-detection (deeper than the first detection)
    and the first detection, with half the gap as the uncertainty."""
    t1 = det["mjd"].min()
    m1 = det.loc[det["mjd"].idxmin(), "magpsf"]
    out = {"t_first": float(t1), "t_exp": np.nan, "t_exp_err": np.nan, "gap_rest": np.nan}
    if non is None or len(non) == 0:
        return out
    before = non[(non["mjd"] < t1) & (non["diffmaglim"] > m1)]
    if before.empty:
        return out
    t0 = before["mjd"].max()
    gap = (t1 - t0) / (1 + z)
    out["gap_rest"] = float(gap)
    if gap <= C.EXPLOSION_MAX_GAP_DAYS:
        out["t_exp"] = float(0.5 * (t0 + t1))
        out["t_exp_err"] = float(0.5 * (t1 - t0))
    return out


# ---------------------------------------------------------------------------
# The plateau fit
# ---------------------------------------------------------------------------

def fermi_mag(t, m0, p0, a0, t_pt, w0):
    """The fall from the plateau in magnitudes (Valenti et al. 2016, Eq. 1,
    written for magnitudes): a step of depth a0 at t_pt with width w0, on top
    of a linear decline p0 (mag/day) shared by plateau and tail."""
    return m0 + p0 * (t - t_pt) + a0 / (1.0 + np.exp(-(t - t_pt) / w0))


def fit_plateau_end(phase: np.ndarray, mag: np.ndarray, err: np.ndarray) -> dict | None:
    """Fit the r-band light curve from day 20 on and decide whether it really
    shows a plateau ending in a fall onto a tail.

    The light curve is cut at the first gap longer than 60 days (stray late
    detections, e.g. from a variable host, would otherwise drive the fit), the
    fit is started from a grid of end times (the steepest-decline guess can
    lock onto the post-peak decline instead), and a fit is accepted only if
    (a) the fall is deeper than 0.8 mag and much steeper than the plateau's
    own decline, (b) at least four nights cover >= 20 days of plateau before
    it, (c) at least two nights within 60 days after it show the tail, and
    (d) the end time is known to better than 7 days."""
    order = np.argsort(phase)
    t, m, e = phase[order], mag[order], np.maximum(err[order], 0.03)
    use = t > 20
    t, m, e = t[use], m[use], e[use]
    if len(t) < 8:
        return None
    gaps = np.flatnonzero(np.diff(t) > 60)
    if len(gaps):
        t, m, e = t[: gaps[0] + 1], m[: gaps[0] + 1], e[: gaps[0] + 1]
    if len(t) < 8:
        return None
    best = None
    for t_try in np.arange(35.0, max(t.max() - 5.0, 36.0), 7.0):
        before = m[t < t_try]
        m0 = np.median(before) if len(before) else m[0]
        guess = [m0, 0.01, 1.5, t_try, 3.0]
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                p, cov = curve_fit(fermi_mag, t, m, p0=guess, sigma=e, absolute_sigma=True,
                                   bounds=([-25, -0.02, 0.3, 30, 0.5], [-5, 0.05, 6, 300, 20]),
                                   maxfev=20000)
        except (RuntimeError, ValueError):
            continue
        chi2 = float(np.sum(((m - fermi_mag(t, *p)) / e) ** 2))
        if best is None or chi2 < best[0]:
            best = (chi2, p, cov)
    if best is None:
        return None
    chi2, p, cov = best
    perr = np.sqrt(np.diag(cov))
    m0, p0, a0, t_pt, w0 = p
    plateau = (t < t_pt - w0) & (t > t_pt - 70)
    tail = (t > t_pt + w0) & (t < t_pt + 60)
    n_before, n_after = int(plateau.sum()), int(tail.sum())
    span_before = float(t[plateau].max() - t[plateau].min()) if n_before else 0.0
    steep = a0 / (4 * w0)            # mag/day at the middle of the fall
    ok = (a0 > 0.8 and steep > 3 * abs(p0) and n_before >= 4 and span_before >= 20
          and n_after >= 2 and np.isfinite(perr[3]) and perr[3] < 7 and w0 < 15
          and chi2 / max(len(t) - 5, 1) < 15)
    return {"m0": m0, "p0": p0, "a0": a0, "t_pt": t_pt, "t_pt_err": float(perr[3]), "w0": w0,
            "n_before": n_before, "n_after": n_after, "span_before": span_before,
            "chi2_dof": chi2 / max(len(t) - 5, 1), "accepted": bool(ok)}


def interp_band(nl: pd.DataFrame, fid: int, phase_col: str, at: float, window: float = 10.0):
    """Magnitude at a phase, linearly interpolated between nightly points no
    more than `window` days away on each side (or the nearest single point
    within 3 days)."""
    b = nl[nl["fid"] == fid]
    if b.empty:
        return np.nan
    x, y = b[phase_col].to_numpy(), b["mag"].to_numpy()
    lo, hi = x[x <= at], x[x >= at]
    if len(lo) and len(hi) and (at - lo.max() <= window) and (hi.min() - at <= window):
        return float(np.interp(at, x, y))
    j = np.argmin(np.abs(x - at))
    return float(y[j]) if abs(x[j] - at) <= 3 else np.nan


def bolometric_correction_g(g_r) -> np.ndarray:
    """Lyman, Bersier & James (2014) for Type II supernovae: BC_g(g-r)."""
    c0, c1, c2 = C.LYMAN_BC_II
    x = np.clip(np.asarray(g_r, float), *C.LYMAN_BC_COLOR_RANGE)
    return c0 + c1 * x + c2 * x ** 2


def lbol_from_gr(M_g, M_r) -> float:
    m_bol = M_g + bolometric_correction_g(M_g - M_r)
    return float(L_SUN * 10 ** (-0.4 * (m_bol - C.M_BOL_SUN)))


def nickel_mass(L_tail: float, t_rest: float) -> float:
    """56Ni mass that powers L_tail at t_rest days after explosion, assuming
    every gamma ray is trapped (Nadyozhin 1994 decay energies)."""
    per_msun = NI_E1 * np.exp(-t_rest / NI_TAU1) + NI_E2 * np.exp(-t_rest / NI_TAU2)
    return L_tail / per_msun


# ---------------------------------------------------------------------------
# Explosion properties (Goldberg, Bildsten & Paxton 2019)
# ---------------------------------------------------------------------------

#: log L50(erg/s) = a + bM log M10 + bE log E51 + bR log R500  (their Eq. 10)
GOLDBERG_L50 = (42.16, -0.40, 0.74, 0.76)
#: log t_p(days) = a + bNi log M_Ni + bM log M10 + bE log E51  (their Eq. 13)
GOLDBERG_TP = (2.184, 0.134, 0.411, -0.282)


def explosion_properties(L50: float, t_p: float, m_ni: float, R_rsun: float) -> tuple[float, float]:
    """Solve Goldberg et al.'s two scalings for (E51, M_ej in M_sun)."""
    aL, bM, bE, bR = GOLDBERG_L50
    aT, bN, cM, cE = GOLDBERG_TP
    rhs1 = np.log10(L50) - aL - bR * np.log10(R_rsun / 500.0)
    rhs2 = np.log10(t_p) - aT - bN * np.log10(m_ni)
    A = np.array([[bM, bE], [cM, cE]])
    logM10, logE51 = np.linalg.solve(A, np.array([rhs1, rhs2]))
    return float(10 ** logE51), float(10 * 10 ** logM10)


# ---------------------------------------------------------------------------
# Running it over the sample
# ---------------------------------------------------------------------------

def analyse_one(row: pd.Series) -> dict | None:
    got = load_lc(row["ztf_id"])
    if got is None:
        return None
    det, non = got
    z = float(row["z_cmb"]) if np.isfinite(row["z_cmb"]) else float(row["redshift"])
    out = {"ztf_id": row["ztf_id"], "n_det": int(len(det))}
    ep = explosion_epoch(det, non, z)
    out.update(ep)
    nl = nightly(det)
    t0 = ep["t_exp"] if np.isfinite(ep["t_exp"]) else ep["t_first"]
    nl["phase"] = (nl["mjd"] - t0) / (1 + z)
    out["phase_ref"] = "explosion" if np.isfinite(ep["t_exp"]) else "first detection"
    # absolute magnitudes (same distance, extinction and K conventions as bts.py)
    k = 2.5 * np.log10(1 + z)
    a = {1: C.A_G_PER_AV * row["A_V"], 2: C.A_R_PER_AV * row["A_V"]}
    nl["M"] = nl["mag"] - row["mu"] - nl["fid"].map(a) + k
    out["t_last"] = float(nl["phase"].max())

    r = nl[nl["fid"] == 2]
    fit = fit_plateau_end(r["phase"].to_numpy(), r["M"].to_numpy(), r["err"].to_numpy()) \
        if len(r) >= 8 else None
    if fit:
        out.update({f"fit_{k}": v for k, v in fit.items()})
    # day-50 photometry and bolometric luminosity
    nl_abs = nl.assign(mag=nl["M"])
    Mg50 = interp_band(nl_abs, 1, "phase", 50.0)
    Mr50 = interp_band(nl_abs, 2, "phase", 50.0)
    out.update({"M_g50": Mg50, "M_r50": Mr50})
    if np.isfinite(Mg50) and np.isfinite(Mr50):
        out["L50"] = lbol_from_gr(Mg50, Mr50)
    # the tail: 25-60 days after the plateau ends
    if fit and fit["accepted"]:
        t_tail = fit["t_pt"] + 30.0
        Mg_t = interp_band(nl_abs, 1, "phase", t_tail, window=15)
        Mr_t = interp_band(nl_abs, 2, "phase", t_tail, window=15)
        if np.isfinite(Mr_t):
            if np.isfinite(Mg_t):
                L_tail = lbol_from_gr(Mg_t, Mr_t)
                out["tail_color"] = "g-r measured"
            else:
                # no g on the tail: assume the median tail colour of SNe II
                # that have it (set by the caller), here g - r = 1.2
                L_tail = lbol_from_gr(Mr_t + 1.2, Mr_t)
                out["tail_color"] = "g-r assumed 1.2"
            out["L_tail"] = L_tail
            out["t_tail"] = t_tail
            if out["phase_ref"] == "explosion":
                out["M_ni_tail"] = nickel_mass(L_tail, t_tail)
    return out


def analyse_sample(sample: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in sample.iterrows():
        res = analyse_one(row)
        if res is not None:
            rows.append(res)
    return pd.DataFrame(rows)


def finish_explosions(res: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Nickel for SNe without a measured tail from the L50-nickel correlation of
    those with one; then E and M_ej for the fiducial radius and its range."""
    good = res.get("fit_accepted", pd.Series(False, index=res.index)).fillna(False).astype(bool) \
        & (res["phase_ref"] == "explosion") & res["L50"].notna() \
        & (res["M_r50"] > C.LC_MAX_PLATEAU_MR)
    res = res.copy()
    res["too_luminous"] = res["M_r50"] <= C.LC_MAX_PLATEAU_MR
    res["plateau_ok"] = good
    tail = good & res.get("M_ni_tail", pd.Series(np.nan, index=res.index)).notna() \
        & (res.get("M_ni_tail", 0) > 0)
    x = np.log10(res.loc[tail, "L50"])
    y = np.log10(res.loc[tail, "M_ni_tail"])
    rel = {"n": int(tail.sum())}
    if tail.sum() >= 10:
        b, a = np.polyfit(x, y, 1)
        rel.update({"slope": float(b), "intercept": float(a),
                    "scatter_dex": float(np.std(y - (a + b * x)))})
        res["M_ni"] = np.where(tail, res.get("M_ni_tail"), 10 ** (a + b * np.log10(res["L50"])))
        res["M_ni_source"] = np.where(tail, "tail", "L50 relation")
    else:
        res["M_ni"] = np.where(tail, res.get("M_ni_tail"), 0.03)
        res["M_ni_source"] = np.where(tail, "tail", "assumed 0.03")
    for label, R in (("", C.R_PROGENITOR_RSUN), ("_Rlo", C.R_PROGENITOR_RANGE_RSUN[0]),
                     ("_Rhi", C.R_PROGENITOR_RANGE_RSUN[1])):
        E, M = [], []
        for _, r in res.iterrows():
            if r["plateau_ok"] and np.isfinite(r["M_ni"]) and r["M_ni"] > 0:
                e, m = explosion_properties(r["L50"], r["fit_t_pt"], r["M_ni"], R)
            else:
                e, m = np.nan, np.nan
            E.append(e)
            M.append(m)
        res["E51" + label] = E
        res["M_ej" + label] = M
    return res, rel
