"""The Higgs search from start to finish, and what each choice does to the answer.

1. **Replica check.** Our re-implementation of the outreach analysis is
   compared bin by bin with the twelve plots CMS published (record 12360).
2. **What the published analysis does not do.** Its 2e2mu data are the sum of
   the DoubleMuParked and DoubleElectron streams, and a collision that fired
   both triggers is in both (same run, luminosity block and event number):
   counted twice. Its same-flavour pairing keeps the best mass in an integer.
3. **Significance, step by step**: the published histogram as released, the
   same with duplicates removed, our own selection, and our selection with the
   misidentified-lepton (Z+X) background included.
4. **Mass and signal strength** with statistical and systematic uncertainties,
   after the lepton-scale calibration; the Z -> 4l peak checks that
   calibration in the four-lepton events themselves.
5. **Pseudo-experiments**: the distribution of the test statistic, the
   expected significance, and the look-elsewhere effect.
"""

from __future__ import annotations

import json

import numpy as np
from iminuit import Minuit
from scipy import stats

import calibration
import config as cfg
import events
import fourlepton as fl
import higgsfit as hf
import outreach_plots

CACHE = cfg.CACHE_DIR / "higgs_search.json"
VARIABLES = {"Higgs_mass": "m4l", "Z1_mass": "mz1", "Z2_mass": "mz2"}
STATES = ("4mu", "4e", "2e2mu")
DATA_STREAMS = {"4mu": ("Run2012B_DoubleMuParked", "Run2012C_DoubleMuParked"),
                "4e": ("Run2012B_DoubleElectron", "Run2012C_DoubleElectron"),
                "2e2mu": cfg.DATA_4L}
MC_FOR = {"4mu": "ZZTo4mu", "4e": "ZZTo4e", "2e2mu": "ZZTo2e2mu"}
#: Resolution in data relative to simulation, from the Z peaks (calibration.py):
#: muons ~1.0, electrons ~1.1; 2e2mu in between.
SEED = 20120704          # the day the discovery was announced


def _w(name):
    """Event weight exactly as skim.cxx computes it (single precision)."""
    f = np.float32
    w = f(cfg.XSEC_PB[name]) / f(cfg.N_GENERATED[name])
    if name != "SMHiggsToZZTo4L":
        w = w * f(cfg.ZZ_SCALE_FACTOR)
    return float(w * (f(11.58) * f(1000.0)))


# ---------------------------------------------------------------------------
# 1-2. The replica against the published plots
# ---------------------------------------------------------------------------


def replica_histograms(tag: str = "", dedupe: bool = False) -> dict:
    """{state: {variable: {'data', 'signal', 'background'}}} from our replica."""
    reps = {n: events.load_replica(n) for n in (*cfg.SIMULATION, *cfg.DATA_4L)}
    out = {}
    for st in STATES:
        out[st] = {}
        for var, key in VARIABLES.items():
            lo, hi = outreach_plots.RANGES[var]
            edges = np.linspace(lo, hi, outreach_plots.NBINS + 1)
            vals, seen = [], set()
            for n in DATA_STREAMS[st]:
                r = reps[n][st + tag]
                for run, lumi, ev, v in zip(r["run"], r["lumi"], r["event"], r[key]):
                    k = (int(run), int(lumi), int(ev))
                    if dedupe and k in seen:
                        continue
                    seen.add(k)
                    vals.append(v)
            data = np.histogram(vals, edges)[0]
            sig = np.histogram(reps["SMHiggsToZZTo4L"][st + tag][key], edges)[0] * _w(
                "SMHiggsToZZTo4L")
            bkg = np.histogram(reps[MC_FOR[st]][st + tag][key], edges)[0] * _w(MC_FOR[st])
            out[st][var] = {"edges": edges, "data": data, "signal": sig, "background": bkg}
    return out


def compare_with_published() -> dict:
    pub = outreach_plots.load()
    res = {}
    for tag in ("", "_fixed"):
        ours = replica_histograms(tag)
        rows = {}
        for st in STATES:
            for var in VARIABLES:
                p, o = pub[st][var], ours[st][var]
                rows[f"{st}/{var}"] = {
                    "data_identical": bool(np.array_equal(np.array(p["data"]), o["data"])),
                    "data_bins_differing": int(np.sum(np.array(p["data"]) != o["data"])),
                    "n_data": int(o["data"].sum()),
                    "max_abs_diff_signal": float(np.max(np.abs(np.array(p["signal"]) - o["signal"]))),
                    "max_abs_diff_background": float(np.max(np.abs(np.array(p["background"])
                                                                   - o["background"])))}
        res["as_published" if tag == "" else "integer_pairing_fixed"] = rows
    return res


def duplicates() -> dict:
    """2e2mu events present in both the DoubleMuParked and DoubleElectron streams."""
    seen = {}
    for n in cfg.DATA_4L:
        r = events.load_replica(n)["2e2mu"]
        for run, lumi, ev, m in zip(r["run"], r["lumi"], r["event"], r["m4l"]):
            seen.setdefault((int(run), int(lumi), int(ev)), []).append((n, float(m)))
    dup = {k: v for k, v in seen.items() if len(v) > 1}
    in_plot = sorted(v[0][1] for v in dup.values() if 70 < v[0][1] < 180)
    return {"entries": int(sum(len(v) for v in seen.values())), "unique": len(seen),
            "duplicated_events": len(dup), "duplicated_in_70_180": len(in_plot),
            "duplicated_masses_70_180": in_plot,
            "duplicated_in_peak_bins": sum(121.94 < m < 128.06 for m in in_plot),
            "all_from_both_streams": all({x[0].split("_")[1] for x in v} ==
                                         {"DoubleMuParked", "DoubleElectron"}
                                         for v in dup.values())}


def binned_significance(h: dict) -> dict:
    """q0 from the published-style histogram (combined channels, 70-180 GeV):
    Poisson likelihood with the simulated Higgs template scaled by mu >= 0."""
    n = sum(h[st]["Higgs_mass"]["data"] for st in STATES).astype(float)
    s = sum(h[st]["Higgs_mass"]["signal"] for st in STATES)
    b = sum(h[st]["Higgs_mass"]["background"] for st in STATES)

    def nll(mu):
        nu = np.maximum(mu * s + b, 1e-12)
        return 2 * np.sum(nu - n * np.log(nu))

    m = Minuit(nll, mu=1.0)
    m.errordef = 1.0
    m.limits["mu"] = (0, 50)
    m.migrad()
    q0 = max(nll(0.0) - m.fval, 0.0)
    peak = (h["4mu"]["Higgs_mass"]["edges"][:-1] > 121) & (h["4mu"]["Higgs_mass"]["edges"][:-1] < 128)
    return {"mu_hat": m.values["mu"], "q0": q0, "z_local": float(np.sqrt(q0)),
            "n_data_70_180": int(n.sum()), "n_peak_bins": int(n[peak].sum()),
            "b_peak_bins": float(b[peak].sum()), "s_peak_bins": float(s[peak].sum())}


# ---------------------------------------------------------------------------
# 3-5. Our selection: calibration, fits, toys
# ---------------------------------------------------------------------------


def z4l_scale_check(data: dict, mc: dict) -> dict:
    """Fit the data's Z -> 4l peak (80-100 GeV) with the simulated ZZ shape
    stretched by (1 + delta): delta is the data/simulation 4l mass scale."""
    lo, hi = 80.0, 100.0
    m = data["m4l"][(data["m4l"] > lo) & (data["m4l"] < hi)]
    zz = mc["ZZ"]
    sel = (zz["m4l"] > lo - 10) & (zz["m4l"] < hi + 10)
    grid = np.linspace(lo - 5, hi + 5, 2001)
    dens = hf.kde(zz["m4l"][sel], zz["weight"][sel], 0.4, grid)

    def nll(delta):
        g = grid * (1 + delta)
        d = dens / np.trapezoid(dens, g)
        x = np.linspace(lo, hi, 801)
        norm = np.trapezoid(np.interp(x, g, d), x)
        return -2 * np.sum(np.log(np.interp(m, g, d) / norm + 1e-300))

    ds = np.linspace(-0.03, 0.03, 601)
    prof = np.array([nll(x) for x in ds])
    best, lo_, hi_ = hf.interval(ds, prof)
    return {"n_events": int(len(m)), "delta": float(best), "delta_lo": float(lo_),
            "delta_hi": float(hi_)}


def build(scale=None, *, resolution_ratio=None, os_ss_ratio=1.0, zz_sigma=0.10,
          zx_sigma=0.50):
    data = fl.data_candidates(scale=scale)
    mc = fl.mc_candidates()
    ssd = fl.data_candidates(scale=scale, same_sign_z2=True)
    ssm = fl.mc_candidates(same_sign_z2=True)
    model = hf.Model(mc, ssd, ssm, resolution_ratio=resolution_ratio, zz_sigma=zz_sigma,
                     zx_sigma=zx_sigma, os_ss_ratio=os_ss_ratio)
    return data, mc, ssd, ssm, model


def fit_summary(L: hf.Likelihood, mh_known: float, *, profiles: bool = True) -> dict:
    best = L.fit()
    b0 = L.fit(mu=0.0, mh=125.0)
    q0 = max(b0["nll"] - best["nll"], 0.0)
    at_known = L.fit(mh=mh_known)
    q0_known = max(b0["nll"] - at_known["nll"], 0.0)
    if not profiles:
        return {"mh": float(best["mh"]), "mu": float(best["mu"]), "q0": float(q0),
                "z_local_best": float(np.sqrt(q0)), "z_local_known_mass": float(np.sqrt(q0_known))}
    masses = np.arange(120.0, 130.0 + 1e-9, 0.05)
    prof_m = L.profile_mass(masses)
    mh, mlo, mhi = hf.interval(masses, prof_m)
    mus = np.linspace(0.0, 3.0, 121)
    prof_mu = np.array([L.fit(mu=x, mh=best["mh"])["nll"] for x in mus])
    mu, mulo, muhi = hf.interval(mus, prof_mu)
    return {"mh": float(best["mh"]), "mh_lo": float(mlo), "mh_hi": float(mhi),
            "mu": float(best["mu"]), "mu_lo": float(mulo), "mu_hi": float(muhi),
            "theta_zz": float(best["tz"]), "theta_zx": float(best["tx"]),
            "q0": float(q0), "z_local_best": float(np.sqrt(q0)),
            "q0_at_known_mass": float(q0_known), "z_local_known_mass": float(np.sqrt(q0_known)),
            "mu_at_known_mass": float(at_known["mu"]),
            "profile_mass": {"m": masses.tolist(), "nll": (prof_m - prof_m.min()).tolist()},
            "profile_mu": {"mu": mus.tolist(), "nll": (prof_mu - prof_mu.min()).tolist()}}


def run(*, n_toys: int = 20000, refresh: bool = False) -> dict:
    if CACHE.exists() and not refresh:
        return json.loads(CACHE.read_text())
    out = {"replica": compare_with_published(), "duplicates": duplicates()}

    # significance of the published histogram, with and without the duplicates
    out["binned"] = {"as_published": binned_significance(replica_histograms("")),
                     "duplicates_removed": binned_significance(replica_histograms("", dedupe=True))}

    # lepton-scale calibration from the dilepton Z peaks
    cal = calibration.calibrate()
    out["calibration"] = cal
    scale = calibration.scale_function(cal)
    rr = {"4mu": 1.0, "4e": cal["e"]["barrel"]["data"]["resolution"]
          / cal["e"]["barrel"]["sim"]["resolution"]}
    rr["2e2mu"] = float(np.sqrt(0.5 * (1 + rr["4e"] ** 2)))
    out["resolution_ratio"] = rr
    mh_known = cfg.M_HIGGS[0]

    variants = {}
    # our selection, before calibration, with and without Z+X
    d0, mc, ssd0, ssm, m0 = build(None, resolution_ratio={c: 1.0 for c in fl.CHANNELS})
    out["z4l_check_uncalibrated"] = z4l_scale_check(d0, mc)
    m0.X = {c: 0.0 for c in fl.CHANNELS}
    variants["ours_no_zx_uncalibrated"] = fit_summary(hf.Likelihood(m0, d0), mh_known)
    # calibrated, nominal model
    data, mc, ssd, ssm, model = build(scale, resolution_ratio=rr)
    out["z4l_check_calibrated"] = z4l_scale_check(data, mc)
    variants["ours_calibrated_with_zx"] = fit_summary(hf.Likelihood(model, data), mh_known)
    nozx = build(scale, resolution_ratio=rr)[4]
    nozx.X = {c: 0.0 for c in fl.CHANNELS}
    variants["ours_calibrated_no_zx"] = fit_summary(hf.Likelihood(nozx, data), mh_known)
    out["variants"] = variants
    out["expected_yields"] = model.expected()
    out["signal_shapes"] = model.sig_par
    out["zx_shape"] = model.zx_shape
    out["zx"] = {"n_same_sign_data": model.n_ss, "n_same_sign_sim": model.n_ss_mc,
                 "os_ss_ratio": model.os_ss_ratio,
                 "total_all_masses": float(sum(max(model.n_ss[c] - model.n_ss_mc[c], 0.0)
                                               for c in fl.CHANNELS) * model.os_ss_ratio),
                 "in_fit_range": model.X}
    out["data_events"] = [{"channel": str(c), "m4l": float(m), "mz1": float(a), "mz2": float(b),
                           "run": int(r), "event": int(e), "dataset": str(ds)}
                          for c, m, a, b, r, e, ds in zip(data["channel"], data["m4l"], data["mz1"],
                                                          data["mz2"], data["run"], data["event"],
                                                          data["dataset"])]
    out["same_sign_events"] = [{"channel": str(c), "m4l": float(m)}
                               for c, m in zip(ssd["channel"], ssd["m4l"])]
    lo, hi = 70.0, 180.0
    out["yields_70_180"] = {"data": fl.yields(data, lo, hi),
                            "ZZ": fl.yields(mc["ZZ"], lo, hi, mc["ZZ"]["weight"]),
                            "H": fl.yields(mc["H"], lo, hi, mc["H"]["weight"])}
    out["yields_onshell_ZZ"] = {"data": fl.yields(data, 180.0, 1000.0),
                                "ZZ": fl.yields(mc["ZZ"], 180.0, 1000.0, mc["ZZ"]["weight"])}

    # systematic variations of the mass and significance
    nominal = variants["ours_calibrated_with_zx"]
    syst = {}
    for flav in ("mu", "e"):
        for sign in (+1, -1):
            c2 = json.loads(json.dumps(cal))
            for reg in c2[flav]:
                k = c2[flav][reg]
                # calibration statistics + 0.1% for the line-shape difference
                # between Drell-Yan (data) and ZZ (simulation) Z bosons
                k["scale"] = k["scale"] * (1 + sign * np.hypot(k["scale_err"], 0.001))
            d2 = fl.data_candidates(scale=calibration.scale_function(c2))
            syst[f"scale_{flav}_{'up' if sign > 0 else 'down'}"] = fit_summary(
                hf.Likelihood(model, d2), mh_known, profiles=False)
    m_rr1 = build(scale, resolution_ratio={c: 1.0 for c in fl.CHANNELS})[4]
    syst["resolution_as_simulated"] = fit_summary(hf.Likelihood(m_rr1, data), mh_known,
                                                  profiles=False)
    m_os = build(scale, resolution_ratio=rr, os_ss_ratio=1.5)[4]
    syst["zx_os_ss_1.5"] = fit_summary(hf.Likelihood(m_os, data), mh_known, profiles=False)
    out["systematics"] = {k: {kk: v[kk] for kk in ("mh", "mu", "z_local_best",
                                                   "z_local_known_mass")}
                          for k, v in syst.items()}
    dm = []
    for flav in ("mu", "e"):
        up = syst[f"scale_{flav}_up"]["mh"] - nominal["mh"]
        dn = syst[f"scale_{flav}_down"]["mh"] - nominal["mh"]
        dm.append(0.5 * (abs(up) + abs(dn)))
    dm.append(abs(syst["resolution_as_simulated"]["mh"] - nominal["mh"]))
    out["mass_systematic"] = {"scale_mu": dm[0], "scale_e": dm[1], "resolution": dm[2],
                              "total": float(np.sqrt(np.sum(np.square(dm))))}

    # pseudo-experiments
    rng_seed = SEED
    bg = hf.toy_q0_scans(model, n_toys, seed=rng_seed)
    i_known = int(np.argmin(np.abs(hf.SCAN - mh_known)))
    q_obs = nominal["q0_at_known_mass"]
    p_toy = float(np.mean(bg[:, i_known] >= q_obs))
    sig = hf.toy_q0_scans(model, max(n_toys // 10, 500), seed=rng_seed + 1, mu_sig=1.0,
                          mh_sig=mh_known)
    out["toys"] = {
        "n_background": n_toys, "n_signal": int(sig.shape[0]),
        "p_local_known_mass_toys": p_toy,
        "z_local_known_mass_toys": float(stats.norm.isf(max(p_toy, 0.5 / n_toys))),
        "p_local_known_mass_asymptotic": float(0.5 * stats.chi2.sf(q_obs, 1)),
        "expected_q0_median": float(np.median(sig[:, i_known])),
        "expected_z_median": float(np.sqrt(np.median(sig[:, i_known]))),
        "lee": hf.gross_vitells(nominal["q0"], bg, u0=1.0),
        "q0_bg_hist": np.histogram(bg[:, i_known], bins=np.linspace(0, 16, 65))[0].tolist(),
        "max_q0_bg_hist": np.histogram(bg.max(axis=1), bins=np.linspace(0, 20, 81))[0].tolist(),
        "scan_masses": hf.SCAN.tolist(),
    }
    # the observed local-significance scan (nominal model)
    L = hf.Likelihood(model, data)
    out["scan_q0"] = L.scan_q0().tolist()
    CACHE.write_text(json.dumps(out, indent=1, default=_json))
    return out


def _json(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    raise TypeError(type(o))


if __name__ == "__main__":
    r = run(refresh=True)
    print(json.dumps({k: r[k] for k in ("duplicates", "binned", "mass_systematic")}, indent=1))
    for k, v in r["variants"].items():
        print(k, {kk: round(v[kk], 3) for kk in ("mh", "mh_lo", "mh_hi", "mu", "z_local_best",
                                                  "z_local_known_mass")})
    print("toys", {k: v for k, v in r["toys"].items() if not k.endswith("hist")
                   and k != "scan_masses"})
