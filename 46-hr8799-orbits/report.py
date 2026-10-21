"""Collect every number for the README into data/results.json, print a
summary, and draw the figures and the animation.

    python report.py              # everything (uses cached fits / runs)
    python report.py figures      # just the figures

Slow pieces (diagnostics, prediction tests) are cached as JSON in data/cache.
"""

from __future__ import annotations

import json
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

import analysis as AN
import astrometry as A
import config as C
import fits
import kepler as K
import predict
import published as P

FINAL_FIT = "ALL/near1248_coplanar"


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x


def cached(name: str, fn, *args, **kw):
    path = C.CACHE / f"report_{name}.json"
    if path.exists():
        return json.loads(path.read_text())
    val = _jsonable(fn(*args, **kw))
    path.write_text(json.dumps(val, indent=1))
    return val


# ---------------------------------------------------------------------------
# pieces
# ---------------------------------------------------------------------------

def data_summary() -> dict:
    df = A.load(include_excluded=True)
    out = {}
    for name in ("W18", "ALL"):
        d = A.dataset(name)
        out[name] = {"n": int(len(d)), "numbers": int(2 * len(d)), "first": float(d.year.min()),
                     "last": float(d.year.max()), "by_source": d.groupby("source").size().to_dict()}
    out["checks_only"] = df[df.source.isin(["AMI-N26"])].groupby("planet").size().to_dict()
    return out


def _orbitize_columns(post, labels) -> dict:
    return {lab.rstrip("1"): post[:, k] for k, lab in enumerate(labels)}


def _orbitize_chi2(col: dict, df: pd.DataFrame, tau_ref: float = C.TAU_REF_MJD) -> np.ndarray:
    """chi2 of each orbitize! sample against its planet's data, using orbitize!'s
    own Kepler solver (so no convention of ours enters)."""
    from orbitize import kepler as OK
    mtot = col["mtot"] if "mtot" in col else col["m0"]
    ra, dec, _ = OK.calc_orbit(df.mjd.to_numpy(), col["sma"], col["ecc"], col["inc"], col["aop"], col["pan"],
                               col["tau"], col["plx"], mtot, tau_ref_epoch=tau_ref)
    ra, dec = np.atleast_2d(ra), np.atleast_2d(dec)            # (epochs, samples)
    sr, sd, rho = df.ra_err.to_numpy()[:, None], df.dec_err.to_numpy()[:, None], df["corr"].fillna(0).to_numpy()[:, None]
    x, y = (ra - df.ra.to_numpy()[:, None]) / sr, (dec - df.dec.to_numpy()[:, None]) / sd
    return ((x ** 2 + y ** 2 - 2 * rho * x * y) / (1 - rho ** 2)).sum(axis=0)


def orbitize_vs_ours() -> dict:
    """orbitize! (OFTI and MCMC, one planet at a time) against our joint
    Unconstrained fit of the same 2004-2016 data."""
    t = fits.load("W18/unconstrained").table()
    data = A.dataset("W18")
    out = {}
    for method in ("ofti", "mcmc"):
        out[method] = {}
        for p in C.PLANETS:
            path = C.CACHE / f"{method}_W18_{p}.npz"
            if not path.exists():
                continue
            z = np.load(path)
            col = _orbitize_columns(z["post"], list(z["labels"]))
            rows = {}
            for ours, theirs, conv in (("a", "sma", 1.0), ("e", "ecc", 1.0), ("i", "inc", 180 / np.pi),
                                       ("Omega", "pan", 180 / np.pi)):
                x_o = col[theirs] * conv
                x_u = t[f"{ours}_{p}"].to_numpy()
                if ours == "Omega":
                    x_o, x_u = K.wrap180(x_o), K.wrap180(x_u)
                rows[ours] = {"orbitize": AN.q(x_o), "ours": AN.q(x_u)}
            dp = data[data.planet == p]
            rng = np.random.default_rng(C.RANDOM_SEED)
            k = rng.choice(len(z["post"]), min(2000, len(z["post"])), replace=False)
            chi2 = _orbitize_chi2({c: v[k] for c, v in col.items()}, dp)
            out[method][p] = {"n_samples": int(len(z["post"])), "params": rows, "n_points": int(2 * len(dp)),
                              "chi2_median": float(np.median(chi2)), "chi2_min": float(chi2.min())}
    return out


def betapic(r: dict | None = None) -> dict:
    """beta Pic b: our nested-sampling orbit, its 2017 passage, and orbitize!'s
    MCMC of the same astrometry as a cross-check."""
    import betapic as B
    r = r or B.analyse()
    out = {k: v for k, v in r.items() if not k.startswith("_")}
    path = C.CACHE / "mcmc_betapic.npz"
    if path.exists():
        z = np.load(path)
        col = _orbitize_columns(z["post"], list(z["labels"]))
        s = r["_samples"]
        P_o = K.period_days(col["sma"], col["mtot"]) / 365.25
        data = B.load_data()
        mjd = data.epoch.to_numpy(float)
        sep, pa, sep_e, pa_e = (data[c].to_numpy(float) for c in ("sep", "pa", "sep_err", "pa_err"))

        def chi2(a, e, i, w, O, tau, plx, m, tref):
            """chi2 against the astrometry, through orbitize!'s own Kepler solver."""
            from orbitize import kepler as OK
            ra, dec, _ = OK.calc_orbit(mjd, a, e, i, w, O, tau, plx, m, tau_ref_epoch=tref)
            msep, mpa = np.hypot(ra, dec), np.degrees(np.arctan2(ra, dec)) % 360
            dpa = (mpa - pa[:, None] + 180) % 360 - 180
            return (((msep - sep[:, None]) / sep_e[:, None]) ** 2 + (dpa / pa_e[:, None]) ** 2).sum(axis=0)

        k = np.random.default_rng(C.RANDOM_SEED).choice(len(z["post"]), 2000, replace=False)
        c_orb = chi2(*(col[n][k] for n in ("sma", "ecc", "inc", "aop", "pan", "tau", "plx", "mtot")), float(z["tau_ref"]))
        c_our = chi2(s["a"], s["e"], s["inc"], s["aop"], s["pan"], s["tau"], s["plx"], s["m"], C.TAU_REF_MJD)
        out["orbitize_mcmc"] = {
            "n_samples": int(len(z["post"])),
            "n_numbers": int(2 * len(data)),
            "chi2_orbitize": {"median": float(np.median(c_orb)), "min": float(c_orb.min())},
            "chi2_ours": {"median": float(np.median(c_our)), "min": float(c_our.min())},
            "a_au": AN.q(col["sma"]), "e": AN.q(col["ecc"]), "i_deg": AN.q(np.degrees(col["inc"])),
            "Omega_folded_deg": AN.q(K.wrap180(np.degrees(col["pan"]))), "P_yr": AN.q(P_o),
            "ours_Omega_folded_deg": AN.q(K.wrap180(np.degrees(s["pan"])))}
    return out


def consistency(prefix: str) -> dict:
    return AN.consistency(prefix)


# ---------------------------------------------------------------------------
# HGCA: the star's wobble
# ---------------------------------------------------------------------------

#: Brandt et al. (2021): m_c = m_d = m_e = 1.25 m_b (Wang et al.'s luminosity masses)
BRANDT_RATIOS = np.array([1.0, 1.0, 1.0, 0.8])            # e, d, c, b (relative to e)
BRANDT_JACOBIAN = {"ra": [2.00, -3.94, -3.01, -24.20], "dec": [1.14, 3.80, -13.74, -26.05]}   # b c d e


def _elements_from_table(t: pd.DataFrame):
    return {p: tuple(t[f"{q}_{p}"].to_numpy() * (np.pi / 180 if q in ("i", "omega", "Omega") else 1)
                     for q in ("a", "e", "i", "omega", "Omega", "tau")) for p in C.INNER_TO_OUTER}


def _mass_scale_posterior(D, mtot, ratios, d, cov):
    """Sum over orbit samples of the Gaussian likelihood of one mass scale."""
    W = np.linalg.inv(cov)
    G = np.einsum("nkj,k->nj", D, ratios) * C.MJUP_MSUN / mtot[:, None]
    a = np.einsum("ni,ij,nj->n", G, W, G)
    b = np.einsum("ni,ij,j->n", G, W, d)
    s, se = b / a, 1 / np.sqrt(a)
    res = d[None, :] - s[:, None] * G
    chi2 = np.einsum("ni,ij,nj->n", res, W, res)
    grid = np.linspace(0, 30, 3001)
    pdf = np.exp(-0.5 * ((grid[None] - s[:, None]) / se[:, None]) ** 2).sum(0)
    cdf = np.cumsum(pdf) / pdf.sum()
    lo, med, hi = np.interp([0.16, 0.5, 0.84], cdf, grid)
    p95 = float(np.interp(0.95, cdf, grid))
    return {"median": float(med), "plus": float(hi - med), "minus": float(med - lo), "p95": p95,
            "chi2_median": float(np.median(chi2))}


def hgca_tests() -> dict:
    import hgca as H
    import stability as St
    d4, cov4 = H.measured()
    d, cov = d4[:2], cov4[:2, :2]                      # Gaia minus long-term (as Brandt et al.)
    out = {"measured": {"delta_mu_masyr": d.tolist(), "cov": cov.tolist(),
                        "chi2_vs_zero": float(d @ np.linalg.solve(cov, d))}}
    # (1) Brandt et al.'s set-up: Wang's stable coplanar orbits, scaled to the EDR3 parallax
    t = fits.load("W18/near1248_coplanar").table()
    st = St.stable_systems("W18/near1248_coplanar")
    sub = t.iloc[st["sample"].to_numpy()].reset_index(drop=True)
    f = sub.plx.to_numpy() / C.DR3_PARALLAX[0]
    el = _elements_from_table(sub)
    el = {p: (el[p][0] * f,) + el[p][1:] for p in el}
    mst = sub.mstar.to_numpy() * f ** 3
    D = H.reflex_design(el, np.full(len(sub), C.DR3_PARALLAX[0]), mst)[:, :, :2]
    J = np.median(D * C.MJUP_MSUN / mst[:, None, None], axis=0) * 1e3
    out["brandt_reproduction"] = {
        "jacobian_uas_per_mjup": {"ra": J[[3, 2, 1, 0], 0].tolist(), "dec": J[[3, 2, 1, 0], 1].tolist()},
        "brandt_jacobian": BRANDT_JACOBIAN,
        "mass_e": _mass_scale_posterior(D, mst, BRANDT_RATIOS, d, cov),
        "brandt_mass_e": P.BRANDT2021_MASS_E}
    # (2) the final orbits, 1998-2025 (Jacobi geometry, masses as fitted are not used here)
    tf = fits.load(FINAL_FIT).table(n=3000)
    elf = _elements_from_table(tf)
    Df = H.reflex_design(elf, tf.plx.to_numpy(), tf.mstar.to_numpy())[:, :, :2]
    mtf = tf.mstar.to_numpy()
    out["final_orbits"] = {"fit": FINAL_FIT, "mass_e": _mass_scale_posterior(Df, mtf, BRANDT_RATIOS, d, cov)}
    # (3) can hot-start masses explain the acceleration?
    for label, (mb, mc) in {"wang_combined": (5.8, 7.2), "hotstart_42Myr": (6.0, 8.7)}.items():
        masses = np.array([mc, mc, mc, mb])                # e, d, c, b
        pred = np.einsum("nkj,k->nj", Df, masses * C.MJUP_MSUN) / mtf[:, None]
        pm = np.median(pred, axis=0)
        chi2 = float((d - pm) @ np.linalg.solve(cov, d - pm))
        out["final_orbits"][label] = {"masses_b_cde": [mb, mc], "predicted_masyr": pm.tolist(),
                                      "chi2_2dof": chi2, "amplitude_ratio_measured_over_predicted":
                                          float(np.linalg.norm(d) / np.linalg.norm(pm))}
    # (4) Gaia DR3 RUWE of HR 8799 against stars of similar brightness and colour
    g = pd.read_csv(C.RAW / "gaia_dr3_hr8799_betapic.csv")
    hr = g[g.source_id == C.GAIA_DR3_SOURCE_ID].iloc[0]
    comp = pd.read_csv(C.RAW / "gaia_dr3_bright_comparison.csv").dropna(subset=["ruwe"])
    out["ruwe"] = {"hr8799": float(hr.ruwe), "excess_noise_mas": float(hr.astrometric_excess_noise),
                   "n_comparison": int(len(comp)), "percentile": float(100 * np.mean(comp.ruwe < hr.ruwe)),
                   "comparison_median": float(comp.ruwe.median())}
    return out


# ---------------------------------------------------------------------------
# Prediction tests
# ---------------------------------------------------------------------------

def fitted_jitters(key: str = "ALL/unconstrained") -> dict:
    f = fits.load(key)
    m = f.model
    return {s: float(np.median(10 ** f.equal[:, m.labels.index(f"log10_jit_{s}")])) for s in m.data.jitter_sources}


def _stable_samples(key: str):
    import stability as St
    f = fits.load(key)
    st = St.stable_systems(key)
    return f.model, f.equal[st["sample"].to_numpy()]


def prediction_tests() -> dict:
    jit = fitted_jitters()
    rows = predict.test_rows()
    out = {"jitter_used_mas": jit, "n_measurements": int(len(rows)),
           "by_source_counts": rows.groupby("source").size().to_dict()}
    sets = {}
    for key in ("W18/unconstrained", "W18/near1248_coplanar", "W18R/unconstrained", "W18R/near1248_coplanar"):
        f = fits.load(key)
        sets[key] = (f.model, f.equal)
    sets["W18/near1248_coplanar (stable)"] = _stable_samples("W18/near1248_coplanar")
    tables = {}
    for name, (model, samples) in sets.items():
        tab = predict.compare(model, samples, rows, jitter=jit)
        tables[name] = tab
        summ = predict.summary(tab)
        by = tab.groupby("source").agg(n=("res_mas", "size"), median_miss_mas=("res_mas", "median"),
                                       median_tangential=("res_tangential", "median"),
                                       median_radial=("res_radial", "median"), chi2=("mahalanobis2", "sum"))
        g = tab[tab.source == "GRAVITY-L19"].iloc[0]
        out[name] = {"summary": summ, "by_source": by.round(3).to_dict(orient="index"),
                     "gravity_2018_e": {"predicted": [g.pred_ra, g.pred_dec], "measured": [g.obs_ra, g.obs_dec],
                                        "miss_mas": g.res_mas, "sigma_pred_mas": g.sigma_pred_mas,
                                        "radial": g.res_radial, "tangential": g.res_tangential}}
    out["gravity2019_wang_stable_prediction"] = P.GRAVITY2019_PREDICTION
    pd.concat({k: v for k, v in tables.items()}, names=["fit"]).to_csv(C.DERIVED / "prediction_tests.csv")
    return out


def witp_gravity_check() -> dict:
    """Final orbits vs the unpublished GRAVITY positions inside whereistheplanet's
    data file (never fitted), and vs its 2026 predictions."""
    import h5py
    from astropy.time import Time
    f = h5py.File(C.RAW / "witp_post_hr8799.hdf5", "r")
    data = np.asarray(f["data"])
    planet_of = {1: "b", 2: "c", 3: "d", 4: "e"}
    yr = Time(data["epoch"], format="mjd").decimalyear
    rows = []
    for r, y in zip(data, yr):
        if y > 2019 and r["quant_type"].decode() == "radec" and r["quant1_err"] < 2.0:
            rows.append(dict(planet=planet_of[int(r["object"])], mjd=float(r["epoch"]), year=float(y),
                             ra=r["quant1"], ra_err=r["quant1_err"], dec=r["quant2"], dec_err=r["quant2_err"],
                             corr=r["quant12_corr"], source="GRAVITY-unpublished", kind="radec", date=""))
    g = pd.DataFrame(rows)
    fit = fits.load(FINAL_FIT)
    tab = predict.compare(fit.model, fit.equal, g)
    out = {"n_points": int(len(g)), "median_miss_mas": float(tab.res_mas.median()),
           "max_miss_mas": float(tab.res_mas.max()),
           "median_sigma_pred_mas": float(tab.sigma_pred_mas.median()),
           "chi2": float(tab.mahalanobis2.sum()), "dof": int(2 * len(tab)),
           "per_planet_median_miss": tab.groupby("planet").res_mas.median().to_dict()}
    # 2026 Oct 8 positions: ours vs whereistheplanet's posterior (20 samples). Its
    # positions are computed exactly as whereistheplanet.py does (Keplerian orbit
    # about the star and the inner planets, plus the inner planets' pull on the
    # star) - this reproduces the file's own GRAVITY points to < 1.5 mas, whereas
    # orbitize! 3.4's System, reading the v2.0 file, misses them by up to 100 mas.
    t_now = Time("2026-10-08").mjd
    ours = predict.predicted_positions(fit.model, fit.equal[:3000], np.array([t_now]))
    post = np.asarray(f["post"])
    tref = float(f.attrs["tau_ref_epoch"])
    out["today"] = {}
    for k, p in enumerate(("b", "c", "d", "e")):
        ra_w, dec_w = _witp_predict(post, k, np.array([t_now]), tref)
        ra_o, dec_o = ours[p][0][:, 0], ours[p][1][:, 0]
        out["today"][p] = {"ours": [float(np.median(ra_o)), float(np.median(dec_o))],
                           "ours_sigma": [float(np.std(ra_o)), float(np.std(dec_o))],
                           "witp": [float(np.median(ra_w)), float(np.median(dec_w))],
                           "difference_mas": float(np.hypot(np.median(ra_o) - np.median(ra_w),
                                                            np.median(dec_o) - np.median(dec_w)))}
    return out


def _witp_predict(post, k, t, tref, n_pl=4):
    """whereistheplanet.print_prediction for planet k (0 = b) of a 4-planet posterior."""
    from orbitize import kepler as OK
    sma, ecc, inc, aop, pan, tau = post[:, 6 * k:6 * k + 6].T
    plx = post[:, 6 * n_pl]
    mp, ms = post[:, -1 - n_pl:-1], post[:, -1]
    smas = post[0, 0:6 * n_pl:6]
    within = np.where(smas <= smas[k])[0]
    ra, dec, _ = OK.calc_orbit(t, sma, ecc, inc, aop, pan, tau, plx, ms + mp[:, within].sum(axis=1),
                               tau_ref_epoch=tref)
    for j in within:
        if j == k:
            continue
        inner = np.where(smas < smas[j])[0]
        imtot = ms + mp[:, inner].sum(axis=1)
        a2, e2, i2, w2, O2, t2 = post[:, 6 * j:6 * j + 6].T
        r2, d2, _ = OK.calc_orbit(t, a2, e2, i2, w2, O2, t2, plx, imtot, tau_ref_epoch=tref)
        ra = ra + mp[:, j] / imtot * r2
        dec = dec + mp[:, j] / imtot * d2
    return np.ravel(ra), np.ravel(dec)


# ---------------------------------------------------------------------------
# Fits: reproduction of Wang et al., the other data sets, instrument offsets
# ---------------------------------------------------------------------------

def reproduction() -> dict:
    """Our four fits of Wang et al.'s data against their Table 2, with pulls."""
    res = AN.w18_reproduction("W18")
    for r in res.values():
        pl = AN.pulls(r["posterior"], r["wang"])
        r["pulls_vs_wang"] = pl
        r["max_abs_pull"] = float(max(abs(v) for v in pl.values()))
    return res


def other_fits() -> dict:
    """The four models fitted to the recalibrated-GPI set (W18R) and to all data
    1998-2025 (ALL); evidences, chi2_nu, and the final fit's masses and jitters."""
    out = {}
    for prefix in ("W18R", "ALL"):
        for kind in C.MODELS:
            key = f"{prefix}/{kind}"
            f = fits.load(key)
            g = AN.goodness_of_fit(f)
            out[key] = {"posterior": AN.posterior_summary(f), "logz": f.logz, "logzerr": float(f.raw["logzerr"]),
                        "chi2nu": AN.q(g["chi2nu"]), "n_data": g["n_data"], "k": g["k"],
                        "max_lnL": float(f.raw["logl"].max())}
        base = out[f"{prefix}/unconstrained"]["logz"]
        for kind in C.MODELS:
            out[f"{prefix}/{kind}"]["ln_bayes_vs_unconstrained"] = out[f"{prefix}/{kind}"]["logz"] - base
    f = fits.load("ALLGPIW18/coplanar")
    out["ALLGPIW18/coplanar"] = {"logz": f.logz, "logzerr": float(f.raw["logzerr"])}
    f = fits.load(FINAL_FIT)
    t = f.table()
    lab = f.model.labels
    out["final"] = {
        "fit": FINAL_FIT,
        "masses_mjup": {p: AN.q(t[f"m_{p}"]) for p in C.PLANETS},
        "mass_prior_mjup": list(f.model.mass_range),
        "jitter_mas": {s: AN.q(10 ** f.equal[:, lab.index(f"log10_jit_{s}")]) for s in f.model.data.jitter_sources},
        "plane": {"i": AN.q(t.i_e), "Omega_folded": AN.q(K.wrap180(t.Omega_e.to_numpy()))},
    }
    return out


def _ms(x) -> dict:
    x = np.asarray(x, float)
    return {"mean": float(np.mean(x)), "std": float(np.std(x))}


def calibration() -> dict:
    """Each instrument's north-angle and plate-scale offset, fitted with the orbits
    (coplanar model, all data). GPI twice: as published (2018) and recalibrated."""
    import figures as F
    out = {}
    for key in ("ALL/coplanar", "ALLGPIW18/coplanar"):
        o = F._offset_draws(fits.load(key))
        out[key] = {s: {"rotation_deg": _ms(np.degrees(v[:, 1])), "scale_permille": _ms(1e3 * v[:, 0])}
                    for s, v in o.items()}
    out["gpi_north_correction_deg"] = {
        "published_2018": -out["ALLGPIW18/coplanar"]["GPI-W18"]["rotation_deg"]["mean"],
        "published_2018_err": out["ALLGPIW18/coplanar"]["GPI-W18"]["rotation_deg"]["std"],
        "recalibrated_2020": -out["ALL/coplanar"]["GPI-D20"]["rotation_deg"]["mean"],
        "recalibrated_2020_err": out["ALL/coplanar"]["GPI-D20"]["rotation_deg"]["std"],
        "derosa2020_change": P.DEROSA2020_NORTH_CHANGE}
    return out


# ---------------------------------------------------------------------------
# Stability, resonances, twins, the future
# ---------------------------------------------------------------------------

def _resonance_frame(key: str, tag: str = "") -> tuple[pd.DataFrame, list]:
    path = C.CACHE / f"resonances_{key.replace('/', '_')}{tag}.npz"
    if not path.exists():
        return pd.DataFrame(), []
    recs = list(np.load(path, allow_pickle=True)["records"])
    return pd.DataFrame([{k: v for k, v in r.items() if k != "_series"} for r in recs]), recs


def resonance_summary(key: str, tag: str = "") -> dict:
    """How often each resonant angle librates in the stable systems (Wang Sec. 4.3)."""
    import stability as St
    df, _ = _resonance_frame(key, tag)
    if df.empty:
        return {}
    st = df[df.outcome == "stable"]
    out = {"n_replayed": int(len(df)), "n_stable_on_replay": int(len(st))}
    for name in St.ANGLES:
        f = st[f"lib_{name}"]
        heavy = st[st.m_cde > 6][f"lib_{name}"]
        out[name] = {"median_fraction": float(f.median()), "max_fraction": float(f.max()),
                     "frac_systems_ge90": float((f >= 0.9).mean()), "frac_systems_ge50": float((f >= 0.5).mean()),
                     "frac_systems_lt5": float((f < 0.05).mean()),
                     "n_mcde_gt6": int(len(heavy)),
                     "frac_systems_ge90_mcde_gt6": float((heavy >= 0.9).mean()) if len(heavy) else None}
    return out


def twins_summary(key: str, tag: str = "") -> dict:
    path = C.CACHE / f"twins_{key.replace('/', '_')}{tag}.csv"
    if not path.exists():
        return {}
    df = pd.read_csv(path).drop_duplicates("draw")
    return {"n": int(len(df)), "twin_also_stable": int((df.outcome == "stable").sum()),
            "twin_failed_median_yr": float(df[df.outcome != "stable"].t_end_yr.median())
            if (df.outcome != "stable").any() else None}


def future_summary(key: str, tag: str = "", years: float = 100e6) -> dict:
    import stability as St
    path = C.CACHE / f"future_{key.replace('/', '_')}{tag}_{years / 1e6:.0f}Myr.csv"
    if not path.exists():
        return {}
    df = pd.read_csv(path).drop_duplicates("draw").merge(St.stable_systems(key, tag)[["draw", "m_cde", "m_b"]],
                                                         on="draw")
    out = {"years": years, "n": int(len(df)), "survive": int((df.outcome == "stable").sum()),
           "t_end_percentiles_myr": {str(p): float(np.percentile(df.t_end_yr, p) / 1e6) for p in (10, 25, 50)}}
    bins = [4, 5, 6, 7, 11]
    df["bin"] = pd.cut(df.m_cde, bins)
    g = df.groupby("bin", observed=True)
    out["by_mcde"] = {str(k): {"n": int(len(v)), "survive": int((v.outcome == "stable").sum())} for k, v in g}
    return out


def stability_results() -> dict:
    import run
    import stability as St
    W = P.wang2018_models()
    out = {"runs": {}, "wang_stable_per_million": {m: W[m]["stable_per_million"] for m in C.MODELS}}
    for key, n, prior, tag in run.STABILITY_RUNS:
        if St.results_path(key, tag).exists():
            s = AN.stability_summary(key, tag)
            s.update(planned_draws=n, mass_prior=prior)
            out["runs"][key + tag] = s
    key, tag = run.MAIN_STABLE
    sp = AN.stable_posterior(key, tag)
    out["stable_posterior"] = sp
    out["wang_table3"] = P.wang2018_stable()
    out["wang_text"] = P.WANG2018_TEXT
    out["stable_pulls_vs_wang"] = AN.pulls({k: v for k, v in sp.items() if k != "masses"}, out["wang_table3"])
    out["resonances"] = resonance_summary(key, tag)
    out["resonances_hotstart"] = resonance_summary(key, "_hotstart")
    out["twins"] = twins_summary(key, tag)
    out["future"] = future_summary(key, tag)
    chaos = C.CACHE / "chaos_last_digit_experiment.json"
    if chaos.exists():                      # see stability.read_results
        c = json.loads(chaos.read_text())
        out["chaos_last_digit"] = {"n": c["n"], "n_broke_up": c["n_broke_up"], "what": c["what"],
                                   "t_end_broke_up_myr": sorted(round(x["t_end_perturbed_yr"] / 1e6, 2)
                                                                for x in c["systems"]
                                                                if x["outcome_perturbed"] != "stable")}
    hot = St.results_path(key, "_hotstart")
    if hot.exists():
        h = St.stable_systems(key, "_hotstart")
        out["hotstart_stable_masses"] = {"n": int(len(h)), "m_b": AN.q(h.m_b) if len(h) else None,
                                         "m_cde": AN.q(h.m_cde) if len(h) else None,
                                         "m_b_p5": float(np.percentile(h.m_b, 5)) if len(h) else None,
                                         "m_cde_p5": float(np.percentile(h.m_cde, 5)) if len(h) else None}
    return out


# ---------------------------------------------------------------------------
# Everything
# ---------------------------------------------------------------------------

def results() -> dict:
    """Every number in the README (cached pieces are reused)."""
    out = {
        "data": cached("data_summary", data_summary),
        "reproduction_w18": cached("reproduction", reproduction),
        "other_fits": cached("other_fits", other_fits),
        "goodness_all_with_nuisances": cached("goodness_all", lambda: {
            f"ALL/{m}": AN.goodness_with_nuisances(fits.load(f"ALL/{m}"), n=500) for m in C.MODELS}),
        "consistency": {"W18": cached("consistency_W18", consistency, "W18"),
                        "ALL": cached("consistency_ALL", consistency, "ALL")},
        "calibration": cached("calibration", calibration),
        "orbitize": orbitize_vs_ours(),
        "betapic": cached("betapic", betapic) if (C.CACHE / "report_betapic.json").exists() else betapic(),
        "prediction": cached("prediction_tests", prediction_tests),
        "witp_gravity": cached("witp_gravity_check", witp_gravity_check),
        "hgca": cached("hgca_tests", hgca_tests),
        "keplerian_vs_nbody": cached("keplerian_vs_nbody", AN.keplerian_vs_nbody, FINAL_FIT, 200),
        "stability": stability_results(),
    }
    if "orbitize_mcmc" in out["betapic"] and not (C.CACHE / "report_betapic.json").exists():
        (C.CACHE / "report_betapic.json").write_text(json.dumps(_jsonable(out["betapic"]), indent=1))
    return out


def stage_betapic():
    import betapic as B
    import figures as F
    r = B.analyse()
    print(F.betapic_figure(r)[0])
    print(json.dumps(_jsonable(betapic(r)), indent=1)[:2000])


def stage_prediction():
    import figures as F
    cached("prediction_tests", prediction_tests)
    cached("witp_gravity_check", witp_gravity_check)
    tab = pd.read_csv(C.DERIVED / "prediction_tests.csv", index_col=[0, 1])
    print(F.prediction_figure(tab)[0])


def stage_hgca():
    import figures as F
    import hgca as H
    h = cached("hgca_tests", hgca_tests)
    tf = fits.load(FINAL_FIT).table(n=1500)
    Df = H.reflex_design(_elements_from_table(tf), tf.plx.to_numpy(), tf.mstar.to_numpy())[:, :, :2]
    print(F.hgca_figure(h, Df, tf.mstar.to_numpy())[0])


def stage_figures():
    """Every static figure, from cached results only."""
    from astropy.time import Time
    import figures as F
    import run
    import stability as St
    res = cached("reproduction", reproduction)
    print(F.reproduction_figure(res)[0])
    print(F.coplanarity_figure(fits.load("W18/unconstrained"), fits.load("W18/coplanar"))[0])
    lnz = {"2004-2016 (Wang et al.'s data)": {m: (fits.load(f"W18/{m}").logz,
                                                   float(fits.load(f"W18/{m}").raw["logzerr"])) for m in C.MODELS},
           "1998-2025 (all data)": {m: (fits.load(f"ALL/{m}").logz,
                                        float(fits.load(f"ALL/{m}").raw["logzerr"])) for m in C.MODELS}}
    print(F.evidence_figure(lnz)[0])
    print(F.calibration_figure(fits.load("ALL/coplanar"), fits.load("ALLGPIW18/coplanar"))[0])
    stage_prediction()
    stage_hgca()
    stage_betapic()
    fu, fc = fits.load("W18/unconstrained"), fits.load("W18/near1248_coplanar")
    key, tag = run.MAIN_STABLE
    ms, ss = _stable_samples(key)
    sets = [("Unconstrained", fu.model, fu.equal), ("Near 1:2:4:8, coplanar, low e", fc.model, fc.equal),
            ("…and stable for 40 Myr (N-body)", ms, ss)]
    print(F.sky_orbits_figure(sets, A.dataset("W18"), stem="sky_orbits_w18",
                              title="Orbits allowed by the 2004-2016 data, as constraints are added "
                                    "(Wang et al. 2018, Fig. 5)")[0])
    summ = {m: AN.stability_summary(f"W18/{m}") for m in C.MODELS if St.results_path(f"W18/{m}").exists()}
    W = P.wang2018_models()
    extra = {k: AN.stability_summary(k) for k in ("ALL/coplanar", "ALL/near1248_coplanar")
             if St.results_path(k).exists()}
    hot = St.stable_systems(key, "_hotstart") if St.results_path(key, "_hotstart").exists() else None
    all_st = St.stable_systems("ALL/near1248_coplanar") if St.results_path("ALL/near1248_coplanar").exists() else None
    print(F.stability_figure(summ, W, St.stable_systems(key, tag), hot, extra=extra, all_masses=all_st)[0])
    df, recs = _resonance_frame(key, tag)
    if len(df):
        print(F.resonance_figure(df, recs)[0])
    fpath = C.CACHE / f"future_{key.replace('/', '_')}{tag}_100Myr.csv"
    tpath = C.CACHE / f"twins_{key.replace('/', '_')}{tag}.csv"
    if fpath.exists() and tpath.exists():
        fut = pd.read_csv(fpath).drop_duplicates("draw").merge(St.stable_systems(key, tag)[["draw", "m_cde", "m_b"]],
                                                                on="draw")
        print(F.future_figure(fut, pd.read_csv(tpath).drop_duplicates("draw"))[0])
    print(F.keeper_figure(fits.load(FINAL_FIT), A.dataset("ALL"), Time("2026-10-08").mjd)[0])


def stage_animation():
    import animate
    import interactive
    animate.main()
    print(interactive.build(FINAL_FIT))


def summary(r: dict) -> str:
    """A short plain-text digest of results.json."""
    L = []
    rep = r["reproduction_w18"]
    L.append("Wang et al. 2018 reproduction (2004-2016 data):")
    for m in C.MODELS:
        x = rep[m]
        L.append(f"  {m:18s} lnZ {x['logz']:8.2f}  lnB {x['ln_bayes_vs_unconstrained']:+6.1f}  "
                 f"chi2nu {x['chi2nu']['median']:.2f}  max|pull| {x['max_abs_pull']:.2f}")
    of = r["other_fits"]
    L.append("All data 1998-2025:")
    g2 = r.get("goodness_all_with_nuisances", {})
    for m in C.MODELS:
        x = of[f"ALL/{m}"]
        nu = g2.get(f"ALL/{m}", {}).get("chi2nu", {}).get("median", float("nan"))
        L.append(f"  {m:18s} lnZ {x['logz']:8.2f}  lnB {x['ln_bayes_vs_unconstrained']:+6.1f}  "
                 f"chi2nu {x['chi2nu']['median']:.2f} (with jitter and offsets {nu:.2f})")
    g = r["calibration"]["gpi_north_correction_deg"]
    L.append(f"GPI north correction asked by the planets: {g['published_2018']:+.3f} ± {g['published_2018_err']:.3f} deg "
             f"(De Rosa+2020: +{g['derosa2020_change']}); after recalibration {g['recalibrated_2020']:+.3f} ± "
             f"{g['recalibrated_2020_err']:.3f}")
    st = r["stability"]
    for k, x in st["runs"].items():
        s = x["stable"]
        L.append(f"  stability {k:32s} {s['k']:4d}/{s['n']:7d} = {s['per_million']:6.0f} per million")
    w = r["witp_gravity"]
    L.append(f"Unpublished GRAVITY 2019-2022 ({w['n_points']} points): median miss {w['median_miss_mas']:.1f} mas, "
             f"max {w['max_miss_mas']:.1f}")
    h = r["hgca"]
    L.append(f"HGCA mass of e (final orbits): {h['final_orbits']['mass_e']['median']:.1f} "
             f"+{h['final_orbits']['mass_e']['plus']:.1f} -{h['final_orbits']['mass_e']['minus']:.1f} MJup")
    b = r["betapic"]
    L.append(f"beta Pic b: conjunction {b['conj2017_iso']}, miss {b['conj2017_sep_mas'][0]:.2f} mas "
             f"= {b['conj2017_sep_stellar_radii'][0]:.0f} stellar radii")
    return "\n".join(L)


def write_all():
    r = results()
    path = C.DATA / "results.json"
    path.write_text(json.dumps(_jsonable(r), indent=1))
    print(f"wrote {path}")
    stage_figures()
    stage_animation()
    print(summary(json.loads(path.read_text())))


if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[f"stage_{sys.argv[1]}"]()
    else:
        write_all()
