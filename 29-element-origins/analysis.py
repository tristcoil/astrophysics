"""Every measurement of L4, part by part, with results cached.

    knee      the alpha-knee of the thick disc and of Gaia-Sausage-Enceladus
    process   prompt vs delayed share of 30 elements (two-process model)
    ncapture  barium and europium: trends, GSE vs disc, scatter, DR4's Eu
    kilonova  AT2017gfo: ejecta masses by nested sampling
    budget    europium: one kilonova vs the whole Galaxy

Slow steps (the disc mixture fits, ~5 min; nested sampling, ~10 min) are
cached in data/cache/ and reused unless --refit is given.
"""

from __future__ import annotations

import json
import pickle
import re

import numpy as np

import galah
from config import CACHE_DIR, ELEMENTS, REFERENCE, RESULTS_PATH

FEH_RANGE_DISC = (-1.6, 0.4)


def _cached(name, fn, refit=False):
    path = CACHE_DIR / f"{name}.pkl"
    if path.exists() and not refit:
        with open(path, "rb") as f:
            return pickle.load(f)
    out = fn()
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(out, f)
    return out


def _clean(o):
    """Make results JSON-safe (numpy -> python, NaN -> None)."""
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


class Context:
    """The GALAH sample and the selections every part shares."""

    def __init__(self):
        import gse_h3

        d = galah.load()
        self.d = d
        self.eu = galah.load_europium(d)
        self.disc = galah.disc_kinematics(d)
        self.giants = galah.giants(d)
        self.dwarfs = galah.dwarfs(d)
        h3 = gse_h3.load()
        sid = d["gaiadr3_source_id"]
        self.gse = {
            "GALAH orbits (Feuillet+21)": galah.gse_dynamical(d),
            "H3: Koppelman+19b box": np.isin(sid, h3["source_id"][h3["koppelman2019b"]]),
            "H3: Helmi+18 box": np.isin(sid, h3["source_id"][h3["helmi2018"]]),
        }
        self.h3_counts = dict(helmi=int(h3["n_helmi"]), koppelman=int(h3["n_koppelman"]))

    def cls(self, name):
        return self.giants if name == "giants" else self.dwarfs


# ---------------------------------------------------------------------------
# Part 1: the alpha-knee
# ---------------------------------------------------------------------------


def part_knee(ctx: Context, refit=False):
    import knee

    d = ctx.d
    res, art = {}, {}

    def disc_fit(cls):
        m = (ctx.cls(cls) & ctx.disc & np.isfinite(d["mg_fe"])
             & (d["fe_h"] > FEH_RANGE_DISC[0]) & (d["fe_h"] < FEH_RANGE_DISC[1]))
        return knee.two_sequence_fit(d["fe_h"][m], d["mg_fe"][m], d["e_mg_fe"][m])

    for cls in ("giants", "dwarfs"):
        fit = _cached(f"disc_knee_{cls}", lambda c=cls: disc_fit(c), refit)
        m = ctx.cls(cls) & ctx.disc & np.isfinite(d["mg_fe"])
        ridge = knee.high_alpha_ridge(d["fe_h"][m], d["mg_fe"][m])
        res[f"disc_{cls}"] = {k: v for k, v in fit.items() if not k.startswith(("grid", "profile"))}
        res[f"disc_{cls}"]["ridge_method"] = {
            k: v for k, v in knee.knee_from_ridge(ridge).items() if k in ("knee", "knee_lo", "knee_hi", "plateau", "slope")}
        art[f"disc_{cls}"] = dict(fit=fit, ridge=ridge)

    for name, sel in ctx.gse.items():
        for cls in ("giants", "all"):
            m = sel & np.isfinite(d["mg_fe"]) & (ctx.giants if cls == "giants" else True)
            tag = re.sub(r"[^a-z0-9]+", "_", f"gse_knee_{name}_{cls}".lower())
            f = _cached(tag, lambda m=m: knee.bootstrap_knee(
                d["fe_h"][m].astype(float), d["mg_fe"][m].astype(float), d["e_mg_fe"][m].astype(float),
                n_boot=300), refit)
            key = f"gse[{name}][{cls}]"
            res[key] = {k: v for k, v in f.items() if k not in ("grid", "lnl", "boot_knees")}
            res[key]["median_feh"] = float(np.median(d["fe_h"][m]))
            art[key] = f
    res["h3_counts"] = ctx.h3_counts
    res["gse_minus_disc_knee"] = (res["gse[GALAH orbits (Feuillet+21)][giants]"]["knee"]
                                  - res["disc_giants"]["knee"])
    res.update(_apogee_check(ctx))
    return res, art


def _apogee_check(ctx: Context) -> dict:
    """APOGEE DR17 against GALAH for stars both surveys observed, and the
    disc knee from APOGEE giants that GALAH puts on disc orbits."""
    from config import APOGEE_DIR

    path = APOGEE_DIR / "apogee_dr17_mgfe.fits"
    if not path.exists():
        return {"apogee": "not fetched"}
    from astropy.table import Table

    import knee

    a = Table.read(path)
    ok = ((a["f_[Fe/H]"] == 0) & (a["f_[Mg/Fe]"] == 0) & ((a["AFlag"] & (1 << 23)) == 0)
          & np.isfinite(a["[Mg/Fe]"]) & np.isfinite(a["[Fe/H]"]) & (a["GaiaEDR3"] > 0))
    a = a[ok]
    sid = np.asarray(a["GaiaEDR3"], dtype=np.int64)
    order = np.argsort(sid)
    d = ctx.d
    pos = np.clip(np.searchsorted(sid[order], d["gaiadr3_source_id"]), 0, sid.size - 1)
    hit = sid[order][pos] == d["gaiadr3_source_id"]
    a_fe = np.where(hit, np.asarray(a["[Fe/H]"])[order][pos], np.nan)
    a_mg = np.where(hit, np.asarray(a["[Mg/Fe]"])[order][pos], np.nan)
    out = {"apogee_n_clean": int(len(a))}
    for cls in ("giants", "dwarfs"):
        m = hit & ctx.cls(cls) & np.isfinite(d["mg_fe"])
        dfe, dmg = d["fe_h"][m] - a_fe[m], d["mg_fe"][m] - a_mg[m]
        out[f"apogee_vs_galah_{cls}"] = dict(n=int(m.sum()), d_feh_median=float(np.median(dfe)),
                                             d_feh_sigma=float(1.4826 * np.median(np.abs(dfe - np.median(dfe)))),
                                             d_mgfe_median=float(np.median(dmg)),
                                             d_mgfe_sigma=float(1.4826 * np.median(np.abs(dmg - np.median(dmg)))))
    # The disc knee from APOGEE alone: its giants, with a kinematic disc cut
    # from Gaia astrometry + APOGEE velocities (Toomre |v - v_LSR| < 180 km/s;
    # Schoenrich et al. 2010 solar motion), parallaxes better than 20%.
    a_ok = ((np.asarray(a["logg"]) > 1.0) & (np.asarray(a["logg"]) < 3.2)
            & (np.asarray(a["Teff"]) > 4000) & (np.asarray(a["Teff"]) < 5500)
            & (np.asarray(a["plx"]) > 5 * np.asarray(a["e_plx"])) & np.isfinite(np.asarray(a["HRV"])))
    a = a[a_ok]
    u, v, w = _uvw(np.asarray(a["RAJ2000"]), np.asarray(a["DEJ2000"]), np.asarray(a["plx"]),
                   np.asarray(a["pmRA"]), np.asarray(a["pmDE"]), np.asarray(a["HRV"]))
    disc = np.sqrt((u + 11.1) ** 2 + (v + 12.24) ** 2 + (w + 7.25) ** 2) < 180.0
    fe, mg = np.asarray(a["[Fe/H]"], float), np.asarray(a["[Mg/Fe]"], float)
    m = disc & (fe > -1.6) & (fe < 0.4)
    out["apogee_disc_giants"] = int(m.sum())
    f = _cached("disc_knee_apogee", lambda: knee.two_sequence_fit(fe[m], mg[m], np.full(m.sum(), 0.03)))
    out["apogee_disc_knee"] = {k: v for k, v in f.items() if not k.startswith(("grid", "profile"))}
    return out


def _uvw(ra, dec, plx_mas, pmra, pmdec, rv):
    """Heliocentric Galactic velocities U (to the Galactic Centre), V, W in km/s."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    c = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, distance=(1.0 / plx_mas) * u.kpc,
                 pm_ra_cosdec=pmra * u.mas / u.yr, pm_dec=pmdec * u.mas / u.yr,
                 radial_velocity=rv * u.km / u.s, frame="icrs").galactic
    vel = c.velocity
    return vel.d_x.to_value(u.km / u.s), vel.d_y.to_value(u.km / u.s), vel.d_z.to_value(u.km / u.s)


# ---------------------------------------------------------------------------
# Part 2: prompt and delayed shares (two-process model)
# ---------------------------------------------------------------------------


def part_process(ctx: Context, knee_res: dict):
    import knee
    import twoprocess as tp

    d = ctx.d
    els = [e.lower() for e in ELEMENTS if e not in ("Fe", "Mg", "Eu")]
    extra = {"eu": (ctx.eu["eu_golden"], ctx.eu["cov_cnn"]),
             "eu_dr3": (ctx.eu["eu_dr3"], np.isfinite(ctx.eu["eu_dr3"])),
             "eu_dr4": (d["eu_fe"], np.ones(d["fe_h"].size, bool))}
    out = {}
    for cls in ("giants", "dwarfs"):
        fit = knee_res[f"disc_{cls}"]
        sel = ctx.cls(cls) & ctx.disc & np.isfinite(d["mg_fe"])
        p_high = knee.high_alpha_probability(fit, d["fe_h"], d["mg_fe"], d["e_mg_fe"])
        out[f"{cls}_ours"] = tp.decompose(d, sel, p_high > 0.5, els, fit["plateau"], n_boot=100,
                                          extra=extra)
    # Griffith et al. (2022)'s own recipe on DR4 dwarfs: their dividing line
    # and their assumed plateau of 0.3 (a like-for-like reproduction).
    sel = ctx.dwarfs & ctx.disc & np.isfinite(d["mg_fe"])
    line = np.where(d["fe_h"] < 0, d["mg_fe"] > 0.12 - 0.13 * d["fe_h"], d["mg_fe"] > 0.12)
    zp = tp.zero_points(d, sel, sel & ~line, [])["mg_fe"]
    out["dwarfs_griffith_recipe"] = tp.decompose(d, sel, line, els, REFERENCE["griffith2022_plateau"] + zp,
                                                 n_boot=100, extra=extra)
    return out


# ---------------------------------------------------------------------------
# Part 3: neutron-capture elements
# ---------------------------------------------------------------------------


def part_ncapture(ctx: Context):
    import ncapture as nc

    d, eu = ctx.d, ctx.eu
    g = ctx.giants
    gse = ctx.gse["GALAH orbits (Feuillet+21)"]
    disc = ctx.disc
    edges = np.round(np.arange(-1.7, 0.51, 0.1), 2)
    quantities = {
        "mg_fe": (d["mg_fe"], d["e_mg_fe"]), "ba_fe": (d["ba_fe"], d["e_ba_fe"]),
        "eu_fe": (eu["eu_best"], None), "eu_mg": (eu["eu_best"] - d["mg_fe"], None),
        "ba_eu": (d["ba_fe"] - eu["eu_best"], None),
    }
    res = {"trends": {}, "gse_minus_disc": {}, "gse_minus_disc_cnn_regime": {}}
    for q, (y, e) in quantities.items():
        res["trends"][q] = {
            "disc": nc.binned(d["fe_h"][g & disc], y[g & disc], edges,
                              err=None if e is None else e[g & disc]),
            "gse": nc.binned(d["fe_h"][g & gse], y[g & gse], np.round(np.arange(-2.0, -0.39, 0.2), 2),
                             min_n=8),
        }
        for key, lo in (("gse_minus_disc", -1.7), ("gse_minus_disc_cnn_regime", -1.0)):
            off = nc.offset_at_fixed_feh(d["fe_h"][g & gse], y[g & gse], d["fe_h"][g & disc],
                                         y[g & disc], np.round(np.arange(lo, -0.69, 0.1), 2))
            res[key][q] = dict(offset=off[0], err=off[1], n=off[2])
    # Every element, GSE minus disc at fixed [Fe/H] (giants): GSE's fingerprint.
    res["fingerprint"] = {}
    for el in ELEMENTS:
        if el == "Fe":
            continue
        y = eu["eu_best"] if el == "Eu" else d[f"{el.lower()}_fe"]
        off = nc.offset_at_fixed_feh(d["fe_h"][g & gse], y[g & gse], d["fe_h"][g & disc],
                                     y[g & disc], np.round(np.arange(-1.5, -0.69, 0.1), 2), n_boot=200)
        res["fingerprint"][el] = dict(offset=off[0], err=off[1], n=off[2])
    res["gse_median"] = {q: float(np.nanmedian(y[g & gse])) for q, (y, _) in quantities.items()}
    res["gse_median"]["fe_h"] = float(np.median(d["fe_h"][g & gse]))
    res["buder2022"] = REFERENCE["buder2022_gse_dynamical"]
    res["ba_eu_pure_r"] = nc.pure_r_ba_eu()
    res["ba_eu_pure_s"] = nc.pure_s_ba_eu()
    res["s_fraction"] = {"Ba": nc.bisterzo_s_fraction("Ba"), "Eu": nc.bisterzo_s_fraction("Eu")}

    # Scatter: GALAH giants (all orbits) and JINAbase metal-poor stars.
    sc_edges = np.round(np.arange(-1.7, 0.31, 0.2), 2)
    res["scatter_galah"] = {q: nc.binned(d["fe_h"][g], quantities[q][0][g], sc_edges,
                                         err=None if quantities[q][1] is None else quantities[q][1][g])
                            for q in ("mg_fe", "ba_fe", "eu_fe")}
    j = nc.load_jinabase()
    j_edges = np.round(np.arange(-4.0, -0.49, 0.5), 2)
    res["scatter_jinabase"] = {}
    for el in ("mg", "ba", "eu"):
        det = ~j[f"{el}_upper"]
        res["scatter_jinabase"][f"{el}_fe"] = nc.binned(j["feh"][det], j[f"{el}_fe"][det], j_edges)
        res["scatter_jinabase"][f"{el}_upper_fraction"] = [
            float(np.mean(j[f"{el}_upper"][(j["feh"] >= lo) & (j["feh"] < hi) & np.isfinite(j[f"{el}_fe"])]))
            if np.isfinite(j[f"{el}_fe"][(j["feh"] >= lo) & (j["feh"] < hi)]).sum() else None
            for lo, hi in zip(j_edges[:-1], j_edges[1:])]
    res["jinabase_n"] = int(len(j["feh"]))

    # DR4's own europium against DR3's measurement of the same stars.
    both = np.isfinite(d["eu_fe"]) & np.isfinite(eu["eu_dr3"])
    # The CNN catalogue holds only stars DR3 did not measure; compare it with
    # DR4's own values for those stars instead.
    cnn = np.isfinite(eu["eu_cnn"]) & np.isfinite(d["eu_fe"])
    res["eu_dr4_vs_dr3"] = dict(
        n=int(both.sum()), r=float(np.corrcoef(d["eu_fe"][both], eu["eu_dr3"][both])[0, 1]),
        sigma=nc.robust_sigma(d["eu_fe"][both] - eu["eu_dr3"][both]),
        median=float(np.median(d["eu_fe"][both] - eu["eu_dr3"][both])))
    res["eu_dr4_vs_cnn"] = dict(
        n=int(cnn.sum()), r=float(np.corrcoef(d["eu_fe"][cnn], eu["eu_cnn"][cnn])[0, 1]),
        sigma=nc.robust_sigma(d["eu_fe"][cnn] - eu["eu_cnn"][cnn]),
        median=float(np.median(d["eu_fe"][cnn] - eu["eu_cnn"][cnn])))
    res["eu_counts"] = {k: int(np.isfinite(v).sum()) for k, v in eu.items() if k.startswith("eu_")}
    return res, j


# ---------------------------------------------------------------------------
# Part 4: the kilonova
# ---------------------------------------------------------------------------


def part_kilonova(refit=False):
    import kilonova as K

    kn = K.Kilonova()
    ns = _cached("kilonova_nested", lambda: K.nested(kn, nlive=400), refit)
    s = ns["samples"]
    frac_red = s[:, 6] / (s[:, 3] + s[:, 6])
    modes = {"all": np.ones(len(s), bool), "purple-dominated": frac_red < 0.5,
             "red-dominated": frac_red >= 0.5}
    res = dict(logz=ns["logz"], logz_err=ns["logz_err"], n_calls=ns["n_calls"],
               max_loglike=ns["max_loglike"], n_points=int(len(kn.mag)),
               n_upper=int(kn.upper.sum()), n_filters=len(kn.bands.ids),
               distance_mpc=kn.distance_mpc, modes={})
    for name, m in modes.items():
        if m.sum() < 20:
            continue
        q = {p: np.percentile(s[m, i], [16, 50, 84]).tolist() for i, p in enumerate(K.PARAMS)}
        q["m_total"] = np.percentile(s[m, 0] + s[m, 3] + s[m, 6], [16, 50, 84]).tolist()
        q["m_purple_red"] = np.percentile(s[m, 3] + s[m, 6], [16, 50, 84]).tolist()
        res["modes"][name] = dict(weight=float(m.mean()), quantiles=q)
    best = ns["best"]
    m = kn.model_mag(best)
    det = ~kn.upper
    res["best"] = dict(zip(K.PARAMS, best.tolist()))
    res["best_rms"] = float(np.std((kn.mag - m)[det]))
    res["villar_rms_our_model"] = None
    v = REFERENCE["villar2017"]
    c = K.C / K.KM
    vil = np.array([v["blue"]["m"], v["blue"]["v"] * c, v["blue"]["T"], v["purple"]["m"],
                    v["purple"]["v"] * c, v["purple"]["T"], v["red"]["m"], v["red"]["v"] * c,
                    v["red"]["T"], v["sigma"]])
    res["villar_loglike_our_model"] = float(kn.log_likelihood(vil))
    res["villar_rms_our_model"] = float(np.std((kn.mag - kn.model_mag(vil))[det]))
    res["best_loglike"] = float(kn.log_likelihood(best))
    return res, dict(kn=kn, samples=s, best=best, villar=vil)


# ---------------------------------------------------------------------------
# Part 5: the europium budget
# ---------------------------------------------------------------------------


def part_budget(ctx: Context, kn_samples):
    import budget as B

    sh = B.europium_shares()
    eu = ctx.eu["eu_best"] + ctx.d["fe_h"]          # [Eu/H]
    sel = ctx.giants & ctx.disc & np.isfinite(eu)
    mean_lin = float(np.mean(10 ** eu[sel]))
    gal = B.galaxy_europium(mean_lin)
    s = kn_samples
    per_event = {k: np.percentile(v, [16, 50, 84]).tolist()
                 for k, v in B.kilonova_europium(s[:, 0], s[:, 3], s[:, 6], sh).items()}
    rates = {"low": REFERENCE["gwtc4_bns_rate"][0],
             "mid": float(np.sqrt(np.prod(REFERENCE["gwtc4_bns_rate"]))),
             "high": REFERENCE["gwtc4_bns_rate"][1]}
    budgets = {}
    for reading, (lo, mid, hi) in per_event.items():
        budgets[reading] = B.merger_budget(mid, gal["m_eu"], list(rates.values()))
    f_hist, tt, rel = B.history_factor()
    return dict(shares=sh, mean_eu_h_linear_disc_giants=mean_lin, n_disc_giants_eu=int(sel.sum()),
                galaxy=gal, per_event_eu=per_event, rates=rates, budgets=budgets,
                history_factor=f_hist, villar_check=B.villar_2017_check(sh),
                history_curve=dict(t_gyr=tt[::20].tolist(), rate_rel=rel[::20].tolist()))


# ---------------------------------------------------------------------------


def run(parts=("knee", "process", "ncapture", "kilonova", "budget"), refit=False):
    results, art = {}, {}
    ctx = Context() if set(parts) & {"knee", "process", "ncapture", "budget"} else None
    if ctx is not None:
        results["sample"] = dict(ctx.d["_counts"], n_giants=int(ctx.giants.sum()),
                                 n_dwarfs=int(ctx.dwarfs.sum()), n_disc=int(ctx.disc.sum()),
                                 n_gse={k: int(v.sum()) for k, v in ctx.gse.items()})
    if "knee" in parts or "process" in parts:
        results["knee"], art["knee"] = part_knee(ctx, refit)
    if "process" in parts:
        results["process"] = part_process(ctx, results["knee"])
    if "ncapture" in parts:
        results["ncapture"], art["jinabase"] = part_ncapture(ctx)
    if "kilonova" in parts or "budget" in parts:
        results["kilonova"], art["kilonova"] = part_kilonova(refit)
    if "budget" in parts:
        results["budget"] = part_budget(ctx, art["kilonova"]["samples"])
    art["ctx"] = ctx
    return results, art


def save_results(results) -> None:
    old = json.loads(RESULTS_PATH.read_text()) if RESULTS_PATH.exists() else {}
    old.update(_clean(results))
    RESULTS_PATH.write_text(json.dumps(old, indent=1))
