"""KIC 9832227: the merger that was predicted for 2022 and did not happen.

Molnar et al. (2017) fitted an exponentially shrinking period to 17 years of
eclipse times and predicted a red-nova merger in 2022.2 +- 0.7. Socia et al.
(2018) showed that the earliest point, from the 1999 NSVS survey, was wrong
by half a day plus one orbit; without it, nothing is accelerating.

Here every eclipse time we can measure ourselves comes from the original
photometry, with one template and one reference ephemeris:

* NSVS 1999-2000 (VizieR II/287), whose time column is labelled
  "JD - 2400000" but holds MJD = JD - 2400000.5;
* SuperWASP 2004-2010; Kepler 2009-2013; ASAS-SN 2014-2018; TESS 2019-2024.

The published compilation of Kovacs et al. (2019) fills the gaps (Vulcan,
HAT, ASAS, AAVSO, MLO) and checks our values where both exist.
"""

from __future__ import annotations

import gzip
import json

import numpy as np
from astropy.io import fits

import config as C
import timing as T

#: Reference ephemeris of Socia et al. (2018) and Kovacs et al. (2019), so our
#: O - C values compare one-to-one with their tables.
T_REF = 2454953.949183   # BJD_TDB of a primary minimum
P_REF = 0.4579515        # days

ARCH_DIR = C.KEPLER_DIR.parent / "archival"
RA, DEC = C.KIC9832227["ra"], C.KIC9832227["dec"]


# ---------------------------------------------------------------------------
# Loaders: each returns dict(t = BJD_TDB, flux = relative flux, err, ...)
# ---------------------------------------------------------------------------


def _norm(flux, err):
    med = np.nanmedian(flux)
    return flux / med, err / med


def load_kepler() -> dict:
    """Kepler long cadence, Q0-Q17: PDCSAP flux, quality = 0, each quarter
    divided by its median."""
    ts, fs, es, qs = [], [], [], []
    for p in sorted(C.KEPLER_DIR.glob("kplr*_llc.fits")):
        with fits.open(p) as h:
            d = h[1].data
            q = h[0].header["QUARTER"]
            good = (d["SAP_QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"]) & np.isfinite(d["TIME"])
            f, e = _norm(d["PDCSAP_FLUX"][good], d["PDCSAP_FLUX_ERR"][good])
            ts.append(d["TIME"][good] + C.KEPLER_BJD_OFFSET)
            fs.append(f)
            es.append(e)
            qs.append(np.full(good.sum(), q))
    return {"t": np.concatenate(ts), "flux": np.concatenate(fs), "err": np.concatenate(es),
            "quarter": np.concatenate(qs)}


#: TESS products in order of preference per sector.
_TESS_PREF = ("tess", "tess-spoc", "qlp")


def load_tess() -> dict:
    """One light curve per TESS sector: SPOC 2-min if it exists, otherwise the
    TESS-SPOC full-frame-image product, otherwise QLP."""
    import re
    by_sector = {}
    for p in sorted(C.TESS_DIR.glob("*.fits")):
        name = p.name.lower()
        if name.startswith("hlsp_qlp"):
            kind, sec = "qlp", int(re.search(r"_s(\d{4})-", name).group(1))
        elif name.startswith("hlsp_tess-spoc"):
            kind, sec = "tess-spoc", int(re.search(r"[-_]s(\d{4})[-_]", name).group(1))
        elif name.startswith("tess"):
            kind, sec = "tess", int(re.search(r"-s(\d{4})-", name).group(1))
        else:
            continue
        by_sector.setdefault(sec, {})[kind] = p
    out = {"t": [], "flux": [], "err": [], "sector": [], "product": [], "cadence_s": []}
    for sec in sorted(by_sector):
        kind = next(k for k in _TESS_PREF if k in by_sector[sec])
        with fits.open(by_sector[sec][kind]) as h:
            d = h[1].data
            cols = d.columns.names
            if kind == "qlp":
                fcol = "KSPSAP_FLUX" if "KSPSAP_FLUX" in cols else "SAP_FLUX"
                ecol = "KSPSAP_FLUX_ERR" if "KSPSAP_FLUX_ERR" in cols else None
                qual = d["QUALITY"] == 0
            else:
                fcol, ecol = "PDCSAP_FLUX", "PDCSAP_FLUX_ERR"
                qual = d["QUALITY"] == 0
            good = qual & np.isfinite(d[fcol]) & np.isfinite(d["TIME"])
            f = d[fcol][good].astype(float)
            e = d[ecol][good].astype(float) if ecol else np.full(f.size, np.nan)
            if not np.isfinite(e).any():
                e = np.full(f.size, 1.4826 * np.median(np.abs(np.diff(f))) / np.sqrt(2))
            f, e = _norm(f, e)
            cad = float(np.median(np.diff(d["TIME"][good]))) * C.DAY_S
        out["t"].append(d["TIME"][good] + C.TESS_BJD_OFFSET)
        out["flux"].append(f)
        out["err"].append(e)
        out["sector"].append(np.full(f.size, sec))
        out["product"].append(np.full(f.size, kind))
        out["cadence_s"].append(np.full(f.size, cad))
    return {k: np.concatenate(v) for k, v in out.items()}


def nsvs_frames() -> dict:
    """NSVS frame table: frame ID (= line number) -> the 'JD' column, which
    holds MJD (JD - 2400000.5) despite the ReadMe's label 'JD-2400000'."""
    mjd = {}
    with gzip.open(ARCH_DIR / "II_287_frames.dat.gz", "rt") as fh:
        for i, ln in enumerate(fh, start=1):
            parts = ln.split("|")
            mjd[i] = (float(parts[3]), parts[5].strip(), float(parts[6]))
    return mjd


def load_nsvs(time_convention: str = "correct") -> dict:
    """NSVS 1999-2000 unfiltered photometry of both object IDs.

    ``time_convention``: 'correct' adds 2400000.5 to the frame times (they are
    MJD; checked against the UT date column); 'as_labelled' adds 2400000.0, as
    the ReadMe's label suggests, which makes every time half a day early.
    Good points: (flags & 32004) == 0, 0 < err < 0.4 mag, 5 < mag < 16
    (Wozniak et al. 2004, Table 3).
    """
    frames = nsvs_frames()
    offset = {"correct": 2400000.5, "as_labelled": 2400000.0}[time_convention]
    t, mag, err, ident, iso = [], [], [], [], []
    for oid in (5597755, 5620022):
        for ln in (ARCH_DIR / f"nsvs_{oid}.dat").read_text().splitlines():
            o, fr, _, _, mm, em, fl = (int(x) for x in ln.split())
            m, e = mm / 1000.0, em / 1000.0
            if (fl & 32004) or not (0 < e < 0.4) or not (5 < m < 16):
                continue
            mjd, date, exptime = frames[fr]
            t.append(mjd + offset)
            mag.append(m)
            err.append(e)
            ident.append(o)
            iso.append(date)
    t = np.array(t)
    bjd = T.jd_utc_to_bjd_tdb(t, RA, DEC)
    mag = np.array(mag)
    flux = 10 ** (-0.4 * (mag - np.median(mag)))
    e = np.array(err) * flux * np.log(10) / 2.5
    return {"t": bjd, "flux": flux, "err": e, "id": np.array(ident), "jd_utc": t, "date": np.array(iso)}


def load_wasp(source: str = "dr1") -> dict:
    """SuperWASP: DR1 (2004-2008, NASA Exoplanet Archive) or VeSPA (2004-2010).
    TMID is the mid-exposure HJD(UTC) in seconds after JD_REF. TAMFLUX2 is the
    detrended flux; each camera gets its own offset and scale in the fits."""
    name = {"dr1": "1SWASP_J192916.01+463719.9_DR1.fits",
            "vespa": "1SWASP_J192916.01+463719.9_vespa_v0.9.fits"}[source]
    with fits.open(ARCH_DIR / name) as h:
        d = h[1].data
        jd_ref = h[0].header["JD_REF"]
        hjd = jd_ref + d["TMID"].astype(float) / C.DAY_S
        f = d["TAMFLUX2"].astype(float)
        e = d["TAMFLUX2_ERR"].astype(float)
        cam = np.array([str(i)[:3] for i in d["IMAGEID"]])
    good = np.isfinite(f) & np.isfinite(e) & (f > 0) & (e > 0)
    hjd, f, e, cam = hjd[good], f[good], e[good], cam[good]
    # remove exact duplicate time stamps (two cameras exposing simultaneously keep both;
    # identical rows from the same camera do not)
    key = np.char.add(np.round(hjd * C.DAY_S).astype(np.int64).astype(str), cam)
    _, first = np.unique(key, return_index=True)
    first = np.sort(first)
    hjd, f, e, cam = hjd[first], f[first], e[first], cam[first]
    bjd = T.hjd_utc_to_bjd_tdb(hjd, RA, DEC)
    med = np.median(f)
    return {"t": bjd, "flux": f / med, "err": e / med, "camera": cam}


def load_asassn() -> dict:
    """ASAS-SN Variable Stars Database light curve: V band, 2014-2018, HJD(UTC)."""
    j = json.loads((ARCH_DIR / "asassn_variables_J192915.86+463719.9.json").read_text())
    r = j["results"]
    hjd = np.array([x["hjd"] for x in r], float)
    mag = np.array([x["mag"] for x in r], float)
    me = np.array([x["mag_err"] for x in r], float)
    cam = np.array([x["camera"] for x in r])
    o = np.argsort(hjd)
    hjd, mag, me, cam = hjd[o], mag[o], me[o], cam[o]
    good = np.isfinite(mag) & (me > 0) & (me < 0.2)
    bjd = T.hjd_utc_to_bjd_tdb(hjd[good], RA, DEC)
    flux = 10 ** (-0.4 * (mag[good] - np.median(mag[good])))
    return {"t": bjd, "flux": flux, "err": me[good] * flux * np.log(10) / 2.5, "camera": cam[good]}


# ---------------------------------------------------------------------------
# Template and timings
# ---------------------------------------------------------------------------


def kepler_template(kep: dict, order: int = 14, n_iter: int = 2):
    """Mean Kepler light curve as a Fourier template, with phase zero at the
    deeper minimum, built after removing each quarter's own O - C so that the
    four-year period drift does not smear it."""
    phase = (kep["t"] - T_REF) / P_REF
    tpl, _ = T.Template.fit(np.mod(phase, 1), kep["flux"], order=order)
    for _ in range(n_iter):
        shifts = np.zeros_like(phase)
        for q in np.unique(kep["quarter"]):
            s = kep["quarter"] == q
            r = T.fit_shift(kep["t"][s], kep["flux"][s], kep["err"][s], tpl, P_REF, T_REF, window=0.1, n_grid=201)
            shifts[s] = r["delta"]
        tpl, _ = T.Template.fit(np.mod(phase - shifts, 1), kep["flux"], order=order)
    return tpl


def segment_times(data: dict, tpl, groups_key=None, segment_days=None, split_gap=60.0, min_points=30,
                  nsig=4.0, window=0.25, predict=None, min_group=25, min_cover=0.8):
    """Template shifts for consecutive segments of a data set.

    Segments: observing seasons (separated by gaps > split_gap days), further
    cut into pieces of ``segment_days`` if given. ``predict`` (a function of
    time returning an expected O - C in days) centres each search window, so
    that the cycle count stays continuous: the two minima of KIC 9832227
    differ in depth by only ~2%, and with spots the template can otherwise
    lock onto the other minimum, half a cycle away. Groups (photometric
    systems) with fewer than ``min_group`` points in a segment are dropped,
    and a segment must cover at least ``min_cover`` of 20 phase bins.
    Returns a list of dicts with the mean time, O - C in days (relative to
    T_REF, P_REF) and its error.
    """
    t = data["t"]
    order = np.argsort(t)
    t_sorted = t[order]
    br = np.nonzero(np.diff(t_sorted) > split_gap)[0] + 1
    seasons = np.split(order, br)
    out = []
    for s in seasons:
        if segment_days:
            edges = np.arange(t[s].min(), t[s].max() + segment_days, segment_days)
            pieces = [s[(t[s] >= a) & (t[s] < b)] for a, b in zip(edges[:-1], edges[1:])]
        else:
            pieces = [s]
        for p in pieces:
            if groups_key:
                g = data[groups_key][p]
                labs, cnt = np.unique(g, return_counts=True)
                p = p[np.isin(g, labs[cnt >= min_group])]
            if p.size < min_points:
                continue
            # the segment must sample most of the orbit, or the shift is ill-defined
            cover = np.unique(np.floor(np.mod((t[p] - T_REF) / P_REF, 1.0) * 20)).size / 20
            if cover < min_cover:
                continue
            g = data[groups_key][p] if groups_key else None
            centre = 0.0 if predict is None else float(predict(np.mean(t[p]))) / P_REF
            r = T.sigma_clip_fit(t[p], data["flux"][p], data["err"][p], tpl, P_REF, T_REF + centre * P_REF,
                                 nsig=nsig, groups=g, window=window)
            if not np.isfinite(r["delta_err"]):
                continue
            out.append({"t": r["t_mean"], "t_min": r["t_min"], "t_max": r["t_max"], "n": int(r["keep"].sum()),
                        "oc_d": (centre + r["delta"]) * P_REF, "oc_err_d": r["delta_err"] * P_REF,
                        "chi2_red": r["chi2_red"]})
    return out


def symmetric_template(tpl):
    """The even-harmonic part of a template: identical minima half a cycle
    apart, insensitive to the odd harmonics that spots and the O'Connell
    effect produce (Molnar et al. 2017 timed with even Fourier orders)."""
    c = tpl.coef.copy()
    for k in range(1, tpl.order + 1):
        if k % 2:
            c[2 * k - 1] = c[2 * k] = 0.0
    return T.Template(c)


#: AAVSO bands kept for timing (visual estimates are too coarse).
AAVSO_BANDS = ("V", "B", "R", "I", "CV", "CR", "TG", "TB", "TR", "SG", "SR", "SI")


def load_aavso() -> dict:
    """AAVSO observers' CCD/DSLR photometry, 2015-2026. Each observer-band pair
    is a separate photometric system, so the fits give each its own offset
    and scale ('system' key). JD is UTC."""
    with gzip.open(ARCH_DIR / "aavso_kic9832227.json.gz") as fh:
        obs = json.load(fh)
    jd, mag, err, sysid = [], [], [], []
    for o in obs:
        band = o.get("Band", "")
        if band not in AAVSO_BANDS:
            continue
        try:
            m = float(o["Mag"])
            e = float(o["Uncertainty"]) if o.get("Uncertainty") not in (None, "") else 0.02
        except (TypeError, ValueError):
            continue
        if o.get("ValFlag") == "T" or not (8.0 < m < 16.0) or e <= 0 or e > 0.1:
            continue
        jd.append(float(o["JD"]))
        mag.append(m)
        err.append(max(e, 0.003))
        sysid.append(f"{o.get('ObsCode', '?')}:{band}")
    jd, mag, err, sysid = np.array(jd), np.array(mag), np.array(err), np.array(sysid)
    o = np.argsort(jd)
    jd, mag, err, sysid = jd[o], mag[o], err[o], sysid[o]
    bjd = T.jd_utc_to_bjd_tdb(jd, RA, DEC)
    # magnitudes -> flux relative to each system's median
    flux = np.empty_like(mag)
    for s in np.unique(sysid):
        k = sysid == s
        flux[k] = 10 ** (-0.4 * (mag[k] - np.median(mag[k])))
    return {"t": bjd, "flux": flux, "err": err * flux * np.log(10) / 2.5, "system": sysid}


# ---------------------------------------------------------------------------
# O - C models
# ---------------------------------------------------------------------------

def year_of(bjd):
    return 2000.0 + (np.asarray(bjd) - C.JD_J2000) / C.YEAR_D


def jd_of(year):
    return C.JD_J2000 + (np.asarray(year) - 2000.0) * C.YEAR_D


def oc_exponential(t, p0, tau, t0, n0):
    """O - C (days) of the exponential ephemeris P(t) = P0 exp(tau / (t - t0))
    (Tylenda et al. 2011; Molnar et al. 2017) against the reference ephemeris.

    The model's cycle count is n0 + integral of dt / P from T_REF to t, and
    O - C = (t - T_REF) - P_REF * cycles. Defined only for t < t0.
    """
    t = np.atleast_1d(np.asarray(t, float))
    lo = min(t.min(), T_REF) - 1.0
    hi = max(t.max(), T_REF) + 1.0
    if hi >= t0:
        return np.full(t.shape, np.nan)
    grid = np.linspace(lo, hi, 20001)
    inv_p = np.exp(-tau / (grid - t0)) / p0
    n = np.concatenate([[0.0], np.cumsum(0.5 * (inv_p[1:] + inv_p[:-1]) * np.diff(grid))])
    n -= np.interp(T_REF, grid, n)
    cyc = n0 + np.interp(t, grid, n)
    return (t - T_REF) - P_REF * cyc


def _jitter_nll(res, err, log_jit):
    s2 = err**2 + np.exp(2 * log_jit)
    return 0.5 * np.sum(res**2 / s2 + np.log(2 * np.pi * s2))


def fit_oc(t, oc, err, model: str, *, fix_jitter=None, x0=None):
    """Maximum-likelihood fit of an O - C model with a jitter term.

    ``model``: 'linear', 'quadratic', 'cubic', 'quad+ltte' (circular light-travel
    orbit) or 'exponential'. Times in BJD, O - C and errors in days. Returns
    best parameters, their errors (inverse Hessian), the jitter, -log L, BIC
    and a callable model.
    """
    from scipy.optimize import minimize
    t = np.asarray(t, float)
    oc = np.asarray(oc, float)
    err = np.asarray(err, float)
    tau_y = (t - T_REF) / C.YEAR_D

    if model in ("linear", "quadratic", "cubic"):
        deg = {"linear": 1, "quadratic": 2, "cubic": 3}[model]

        def f(p, tt):
            return np.polyval(p[::-1], (tt - T_REF) / C.YEAR_D)
        p0 = np.zeros(deg + 1)
        p0[:2] = np.polyfit(tau_y, oc, 1)[::-1]
    elif model == "quad+ltte":
        def f(p, tt):
            y = (tt - T_REF) / C.YEAR_D
            return p[0] + p[1] * y + p[2] * y * y + p[3] * np.sin(2 * np.pi * (y - p[5]) / p[4])
        q = np.polyfit(tau_y, oc, 2)[::-1]
        p0 = np.array([q[0], q[1], q[2], 0.001, 20.0, 0.0]) if x0 is None else np.asarray(x0, float)
    elif model == "exponential":
        def f(p, tt):
            return oc_exponential(tt, P_REF * (1 + p[0] * 1e-6), p[1], jd_of(p[2]), p[3])
        # parameters: (P0 / P_REF - 1) in ppm, tau (d), t0 (year), n0 (cycles)
        p0 = np.array([10.0, 0.1, 2022.0, 0.0]) if x0 is None else np.asarray(x0, float)
    else:
        raise ValueError(model)

    def nll(x):
        p, lj = (x, fix_jitter) if fix_jitter is not None else (x[:-1], x[-1])
        mdl = f(p, t)
        if not np.all(np.isfinite(mdl)):
            return 1e30
        return _jitter_nll(oc - mdl, err, np.log(lj) if fix_jitter is not None else lj)

    start = p0 if fix_jitter is not None else np.concatenate([p0, [np.log(np.median(err))]])
    best = None
    for method in ("Nelder-Mead", "Powell", "Nelder-Mead"):
        r = minimize(nll, start if best is None else best.x, method=method,
                     options={"maxiter": 40000, "maxfev": 40000, "xatol": 1e-10, "fatol": 1e-10})
        if best is None or r.fun < best.fun:
            best = r
    x = best.x
    # numerical Hessian for errors
    h = np.maximum(np.abs(x) * 1e-4, 1e-7)
    n = x.size
    H = np.zeros((n, n))
    f0 = nll(x)
    for i in range(n):
        for j in range(i, n):
            ei = np.zeros(n); ej = np.zeros(n)
            ei[i] = h[i]; ej[j] = h[j]
            H[i, j] = H[j, i] = (nll(x + ei + ej) - nll(x + ei - ej) - nll(x - ei + ej) + nll(x - ei - ej)) / (4 * h[i] * h[j])
    try:
        cov = np.linalg.inv(H)
        perr = np.sqrt(np.clip(np.diag(cov), 0, np.inf))
    except np.linalg.LinAlgError:
        cov, perr = None, np.full(n, np.nan)
    p = x if fix_jitter is not None else x[:-1]
    jit = fix_jitter if fix_jitter is not None else float(np.exp(x[-1]))
    k = n
    bic = 2 * f0 + k * np.log(t.size)
    return {"model": model, "params": p, "errors": perr[:p.size], "jitter_d": jit, "nll": float(f0), "bic": float(bic),
            "n": int(t.size), "k": int(k), "f": (lambda tt, p=p: f(p, np.asarray(tt, float))), "cov": cov}


def pdot_from_quadratic(fit):
    """Period derivative (dimensionless, days per day) of a quadratic O - C in years."""
    c2 = fit["params"][2]           # days / yr^2
    return 2 * c2 * P_REF / C.YEAR_D**2


# ---------------------------------------------------------------------------
# All eclipse timings in one table
# ---------------------------------------------------------------------------

#: Kovacs et al. (2019) data sets we have no photometry for; their template
#: values (OC4) fill the gaps.
LITERATURE_ONLY = ("Vulcan", "HAT-", "ASAS-IV", "WISE-", "MLO")


def measure_all(cache: bool = True):
    """Every O - C we measure ourselves, plus the literature-only epochs.

    Returns (table, template). The table (pandas) has one row per timing:
    source, kind ('ours' / 'literature'), BJD (mean, min, max), O - C and
    error in minutes, number of points. NSVS appears twice: with the correct
    time conversion and with the ReadMe's label taken literally.
    """
    import pandas as pd
    import pickle
    import literature as L
    path = C.DERIVED_DIR / "kic9832227_oc.csv"
    tpath = C.DERIVED_DIR / "kic9832227_template.pkl"
    kep = load_kepler()
    if cache and path.exists() and tpath.exists():
        return pd.read_csv(path), pickle.loads(tpath.read_bytes())
    tpl = kepler_template(kep)
    tpath.write_bytes(pickle.dumps(tpl))
    rows = []

    def add(source, segs, kind="ours"):
        for s in segs:
            rows.append({"source": source, "kind": kind, "t": s["t"], "t_min": s["t_min"], "t_max": s["t_max"],
                         "oc_min": s["oc_d"] * 1440, "err_min": s["oc_err_d"] * 1440, "n": s["n"]})

    # 1. the robust core: Kepler, NSVS, WASP, ASAS-SN (no cycle-count ambiguity)
    add("Kepler", segment_times(kep, tpl, segment_days=90))
    for conv, label in (("correct", "NSVS"), ("as_labelled", "NSVS (time as labelled)")):
        add(label, segment_times(load_nsvs(conv), tpl, groups_key="id", split_gap=1000, min_points=50))
    add("SuperWASP", segment_times(load_wasp("dr1"), tpl, groups_key="camera", min_points=50))
    vespa = load_wasp("vespa")
    late = vespa["t"] > 2454800.0  # VeSPA's 2009-2010 seasons are not in DR1
    add("SuperWASP (VeSPA)", segment_times({k: v[late] for k, v in vespa.items()}, tpl, groups_key="camera",
                                           min_points=50))
    add("ASAS-SN", segment_times(load_asassn(), tpl, split_gap=150, segment_days=365, min_points=40))
    core = pd.DataFrame(rows)
    core = core[core["source"] != "NSVS (time as labelled)"]
    q = np.polyfit((core["t"] - T_REF) / C.YEAR_D, core["oc_min"], 2)

    def pred(tt):
        return np.polyval(q, (np.asarray(tt) - T_REF) / C.YEAR_D) / 1440.0

    # 2. TESS, sector by sector, then AAVSO seasons centred on a cubic through everything so far
    tess = load_tess()
    add("TESS", segment_times(tess, tpl, groups_key="sector", split_gap=20, min_points=200, window=0.2, predict=pred))
    sofar = pd.DataFrame(rows)
    sofar = sofar[sofar["source"] != "NSVS (time as labelled)"]
    c3 = np.polyfit((sofar["t"] - T_REF) / C.YEAR_D, sofar["oc_min"], 3)

    def pred3(tt):
        return np.polyval(c3, (np.asarray(tt) - T_REF) / C.YEAR_D) / 1440.0

    add("AAVSO", segment_times(load_aavso(), tpl, groups_key="system", split_gap=60, segment_days=120,
                               min_points=60, window=0.2, predict=pred3))
    # 3. literature-only epochs (Kovacs et al. 2019, template values OC4)
    for r in L.kovacs_table1():
        if r["source"].startswith(LITERATURE_ONLY):
            rows.append({"source": r["source"] + " (Kovacs+19)", "kind": "literature", "t": r["t_mid"],
                         "t_min": r["t_mid"] - r["half_span"], "t_max": r["t_mid"] + r["half_span"],
                         "oc_min": r["oc"][3], "err_min": r["eoc"][3], "n": np.nan})
    df = pd.DataFrame(rows).sort_values("t").reset_index(drop=True)
    df["year"] = year_of(df["t"])
    df.to_csv(path, index=False)
    return df, tpl


def _wls_jitter(A, y, err, positive_last=False):
    """Linear least squares with a fitted jitter: maximise the Gaussian
    likelihood over the coefficients (closed form for a given jitter) and the
    jitter (1-D search). If ``positive_last``, the last coefficient is kept >= 0."""
    from scipy.optimize import minimize_scalar

    def solve(lj):
        s2 = err**2 + np.exp(2 * lj)
        w = 1 / np.sqrt(s2)
        coef, *_ = np.linalg.lstsq(A * w[:, None], y * w, rcond=None)
        if positive_last and coef[-1] < 0:
            coef2, *_ = np.linalg.lstsq(A[:, :-1] * w[:, None], y * w, rcond=None)
            coef = np.concatenate([coef2, [0.0]])
        r = y - A @ coef
        return 0.5 * np.sum(r**2 / s2 + np.log(2 * np.pi * s2)), coef

    res = minimize_scalar(lambda lj: solve(lj)[0], bounds=(np.log(1e-6), np.log(0.1)), method="bounded")
    nll, coef = solve(res.x)
    return nll, coef, float(np.exp(res.x))


def exponential_profile(t, oc, err, t0_years):
    """Profile likelihood of the exponential period law over its merger date t0.

    For tau << |t - t0| (true except in the last days), the law's O - C is
        O - C = a (t - T_REF) + b + c ln[(t0 - t) / (t0 - T_REF)],
    with c = P_REF tau / P0 >= 0 (a period that shrinks towards t0), linear in
    (a, b, c) for a fixed t0. For each trial t0 we maximise the likelihood over
    (a, b, c, jitter). Returns arrays of t0, -log L, and the coefficients.
    """
    t = np.asarray(t, float)
    out_nll, out_coef, out_jit = [], [], []
    for ty in t0_years:
        t0 = float(jd_of(ty))
        if t0 <= t.max() + 1:
            out_nll.append(np.inf); out_coef.append([np.nan] * 3); out_jit.append(np.nan)
            continue
        A = np.column_stack([(t - T_REF), np.ones_like(t), np.log((t0 - t) / (t0 - T_REF))])
        nll, coef, jit = _wls_jitter(A, oc, err, positive_last=True)
        out_nll.append(nll); out_coef.append(coef); out_jit.append(jit)
    return np.asarray(t0_years), np.array(out_nll), np.array(out_coef), np.array(out_jit)


def exponential_tau(coef, t0_year=None):
    """tau (days) from the linearised coefficients: c = P_REF tau / P0, P0 = P_REF / (1 - a)."""
    a, _, c = coef
    p0 = P_REF / (1 - a)
    return c * p0 / P_REF, p0


def kovacs_model(t_bjd):
    """Kovacs et al. (2019) Eq. 4 with their Table 3 least-squares coefficients:
    O - C [min] = c1 + c2 t + c3 t^2 + c4 sin(2 pi t / P2) + c5 cos(2 pi t / P2),
    t = HJD - 2455000 (their reference ephemeris is ours)."""
    t = np.asarray(t_bjd, float) - 2455000.0
    c1, c2, c3, c4, c5, p2 = -0.76257, -0.015618, -0.21632e-5, 10.659, -2.6196, 4925.0
    return c1 + c2 * t + c3 * t * t + c4 * np.sin(2 * np.pi * t / p2) + c5 * np.cos(2 * np.pi * t / p2)


def molnar_model(t_bjd, n0=0.0):
    """Molnar et al. (2017)'s exponential ephemeris (P0, tau, t0) = (0.4579615 d,
    0.114 d, 2022.2) as O - C in minutes (zero point n0 in cycles is free)."""
    import literature as L
    m = L.MOLNAR2017
    return oc_exponential(t_bjd, m["p0"], m["tau_d"], float(jd_of(m["t0_year"])), n0) * 1440.0


#: The data sets Molnar et al. (2017) had: NSVS, ASAS, WASP 2007-2008, Kepler,
#: then their own Calvin photometry 2013-2016, for which ASAS-SN and AAVSO
#: seasons of the same years stand in here.
MOLNAR_SOURCES = ("NSVS", "ASAS-IV (Kovacs+19)", "SuperWASP", "Kepler", "ASAS-SN", "AAVSO")
