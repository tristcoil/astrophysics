"""Counting and weighing neutrinos with the CMB and BAO: our own fits.

The data combinations form a ladder, each rung adding information:

* **tt** — Planck's binned temperature spectrum (ℓ = 30–2508) with independent
  errors, exactly project D3's likelihood, plus a Gaussian prior on τ from
  Planck's low-ℓ polarisation (D3 held τ fixed instead).
* **tt_bao** — the same plus DESI DR2's distances (the spec's "D3 + D5").
* **plik_tt_tau** — the same temperature data through Planck's own
  foreground-marginalised likelihood (plik_lite, with the full covariance
  between bins), still with the τ prior: isolates D3's shortcut.
* **plik_tt** — + Planck's low-ℓ likelihoods (TT and EE, ℓ = 2–29) instead of
  the τ prior: our version of Planck's "TT+lowE".
* **plik** — temperature *and* polarisation (TT, TE, EE) + low-ℓ: "TT,TE,EE+lowE".
* **plik_lens18** — + Planck 2018's lensing reconstruction.
* **planck18** — + the BAO data Planck 2018 used (6dF, MGS, BOSS DR12): our
  stand-in for Planck 2018's "TT,TE,EE+lowE+lensing+BAO", the validation rung.
* **cmb_desi** — Planck + ACT DR6/Planck PR4 lensing + DESI DR2 BAO: our
  stand-in for DESI DR2's "DESI + CMB" (plik variant).
* **cmb_desi_sn_w0wa** — the same + DES-SN5YR supernovae with dark energy free
  to evolve (w0, wa), for how much that loosens the neutrino mass limit.

For each rung we fit ΛCDM, then ΛCDM + N_eff (one 0.06 eV neutrino, as in
Planck's ``base_nnu``), then ΛCDM + Σm_ν (three degenerate neutrinos, as in
Planck's and DESI's ``base_mnu``), each with a profile likelihood.
"""

from __future__ import annotations

import json

import numpy as np

import cmb_model as cm
import config as cfg
import fetch_data as fd
import fitting as ft
import likelihoods as lk

FIT_DIR = cfg.DERIVED_DIR / "fits"
FIT_DIR.mkdir(parents=True, exist_ok=True)

RUNGS = ("tt", "tt_bao", "plik_tt_tau", "plik_tt", "plik", "plik_lens18", "planck18",
         "cmb_desi", "cmb_desi_sn_w0wa")
#: Which rungs get which extension.
NEFF_RUNGS = ("tt", "tt_bao", "plik_tt_tau", "plik_tt", "plik", "planck18", "cmb_desi")
MNU_RUNGS = ("tt", "tt_bao", "plik_tt", "plik", "plik_lens18", "planck18", "cmb_desi",
             "cmb_desi_sn_w0wa")
W0WA_RUNGS = ("cmb_desi_sn_w0wa",)
#: The published result each rung is graded against (keys of config tables).
COMPARE = {
    "tt": {"neff": "Planck TT+lowE", "mnu": "Planck TT+lowE"},
    "tt_bao": {"mnu": "Planck TT+lowE+BAO"},
    "plik_tt_tau": {"neff": "Planck TT+lowE"},
    "plik_tt": {"neff": "Planck TT+lowE", "mnu": "Planck TT+lowE"},
    "plik": {"neff": "Planck TT,TE,EE+lowE", "mnu": "Planck TT,TE,EE+lowE"},
    "plik_lens18": {"mnu": "Planck TT,TE,EE+lowE+lensing"},
    "planck18": {"neff": "Planck TT,TE,EE+lowE+lensing+BAO",
                 "mnu": "Planck TT,TE,EE+lowE+lensing+BAO"},
    "cmb_desi": {"neff": "DESI DR2 + CMB", "mnu": "DESI DR2 + CMB (plik)"},
    "cmb_desi_sn_w0wa": {"mnu": "DESI DR2 + CMB + DESY5, w0wa"},
}


def _tau_prior():
    m, s = cfg.TAU_PRIOR
    return lk.GaussianPrior("tau", m, s, name="tau prior (Planck lowE)")


_LIKES: dict = {}


def _like(name: str):
    """Each likelihood is built once, on first use (ACT's needs the reduced matrices)."""
    if name not in _LIKES:
        make = {
            "tt": lambda: lk.BinnedTT(fd.load_planck_binned_tt()),
            "tau": _tau_prior,
            "lowl": lambda: lk.PlanckLowlTT(fd.load_lowl_tt()),
            "lowe": lambda: lk.PlanckLowlEE(fd.load_lowl_ee()),
            "bao": lambda: lk.BAO(fd.load_bao()),
            "bao18": lambda: lk.BAO2018(fd.load_bao2018()),
            "plik": lambda: lk.PlikLite(fd.load_plik_lite()),
            "plik_tt": lambda: lk.PlikLite(fd.load_plik_lite(), use=("TT",)),
            "lens18": lambda: lk.PlanckLensing2018(fd.load_planck_lensing()),
            "act": lambda: lk.ACTLensing(fd.load_act_lensing()),
            "desy5": lambda: lk.Supernovae(fd.load_sn("DESY5")),
        }
        _LIKES[name] = make[name]()
    return _LIKES[name]


RUNG_SPEC = {
    "tt": ("TT (D3's likelihood) + τ prior", ("tt", "tau")),
    "tt_bao": ("TT (D3) + τ prior + DESI BAO", ("tt", "tau", "bao")),
    "plik_tt_tau": ("Planck TT + τ prior (no low-ℓ data)", ("plik_tt", "tau")),
    "plik_tt": ("Planck TT + low-ℓ", ("plik_tt", "lowl", "lowe")),
    "plik": ("Planck TT,TE,EE + low-ℓ", ("plik", "lowl", "lowe")),
    "plik_lens18": ("… + Planck lensing", ("plik", "lowl", "lowe", "lens18")),
    "planck18": ("… + lensing + BAO (Planck 2018's data)",
                 ("plik", "lowl", "lowe", "lens18", "bao18")),
    "cmb_desi": ("Planck + ACT lensing + DESI DR2 BAO", ("plik", "lowl", "lowe", "act", "bao")),
    "cmb_desi_sn_w0wa": ("… + DES-SN5YR, w0wa dark energy",
                         ("plik", "lowl", "lowe", "act", "bao", "desy5")),
}


def build_problems(rungs=RUNGS) -> dict:
    """{rung: (label, likelihood list)} for the requested rungs."""
    return {k: (RUNG_SPEC[k][0], [_like(n) for n in RUNG_SPEC[k][1]]) for k in rungs}


def _free(likes, extra=()):
    free = list(cfg.LCDM)
    for lkh in likes:
        free += [n for n in getattr(lkh, "nuisance", ()) if n not in free]
    return free + list(extra)


def _neff_grid(rung: str):
    return {"tt": np.arange(1.6, 5.01, 0.4), "tt_bao": np.arange(2.0, 4.61, 0.3),
            "plik_tt_tau": np.arange(1.8, 4.61, 0.35), "plik_tt": np.arange(1.8, 4.61, 0.35),
            "plik": np.arange(2.0, 4.01, 0.25), "planck18": np.arange(2.3, 3.81, 0.15),
            "cmb_desi": np.arange(2.5, 4.01, 0.15)}[rung].round(3)


def _mnu_grid(rung: str):
    grids = {
        "tt": [0.0, 0.06, 0.15, 0.3, 0.45, 0.6, 0.8, 1.0, 1.3, 1.6],
        "tt_bao": [0.0, 0.03, 0.06, 0.1, 0.15, 0.2, 0.3, 0.4, 0.55, 0.7],
        "plik_tt": [0.0, 0.06, 0.15, 0.3, 0.45, 0.6, 0.8, 1.0, 1.3, 1.6],
        "plik": [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8],
        "plik_lens18": [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8],
        "planck18": [0.0, 0.03, 0.06, 0.09, 0.12, 0.16, 0.2, 0.25, 0.32, 0.4],
        "cmb_desi": [0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.13, 0.16, 0.2, 0.25],
        "cmb_desi_sn_w0wa": [0.0, 0.03, 0.06, 0.1, 0.15, 0.2, 0.27, 0.35, 0.45, 0.6],
    }
    return np.array(grids[rung])


def _load(name):
    path = FIT_DIR / name
    return json.loads(path.read_text()) if path.exists() else None


def _save(name, obj):
    (FIT_DIR / name).write_text(json.dumps(obj, indent=1))


def _start(rung: str, base: dict | None = None) -> dict:
    """A good starting point: D3's fit for the TT rungs, Planck's best fit otherwise."""
    p = {"ombh2": 0.02237, "omch2": 0.1200, "theta100": 1.04092, "logA": 3.044, "ns": 0.9649,
         "tau": cfg.TAU_PRIOR[0], "A_planck": 1.0}
    if base:
        p.update(base)
    return p


def _extra(key: str) -> tuple:
    return ("w0", "wa") if key in W0WA_RUNGS else ()


def run(force: bool = False, rungs=RUNGS) -> dict:
    problems = build_problems(rungs)
    out = {"rungs": {k: problems[k][0] for k in rungs}, "lcdm": {}, "neff": {}, "mnu": {}}
    with cm.TheoryPool(cfg.N_WORKERS) as pool:
        # ---- LCDM on every rung (one massive 0.06 eV neutrino) ----
        for key in rungs:
            label, likes = problems[key]
            name = f"lcdm_{key}.json"
            res = None if force else _load(name)
            if res is None:
                print(f"  fitting LCDM to {label}")
                prob = ft.Problem(label, likes, _free(likes, _extra(key)), nu_model="one")
                start = _start(key)
                if key in W0WA_RUNGS:
                    start.update({"w0": -0.75, "wa": -0.8})
                res = ft.fit(prob, start, pool)
                _save(name, res)
            out["lcdm"][key] = res
            print(f"  LCDM  {label:48s} chi2 = {res['chi2']:8.2f} for {res['n_data']} data, "
                  f"H0 = {res['derived']['H0']:.2f}")
        # ---- N_eff ----
        for key in rungs:
            if key not in NEFF_RUNGS:
                continue
            label, likes = problems[key]
            name = f"neff_{key}.json"
            res = None if force else _load(name)
            if res is None:
                print(f"  N_eff free: {label}")
                prob = ft.Problem(label, likes, _free(likes, ("nnu",)), nu_model="one")
                start = _start(key, out["lcdm"][key]["best"])
                start["nnu"] = cfg.NEFF_SM
                fit = ft.fit(prob, start, pool)
                prof = ft.profile(prob, fit, "nnu", _neff_grid(key), pool)
                res = {"fit": fit, "profile": prof}
                _save(name, res)
            out["neff"][key] = summarise_neff(res, label, key)
        # ---- Sum of neutrino masses ----
        for key in rungs:
            if key not in MNU_RUNGS:
                continue
            label, likes = problems[key]
            name = f"mnu_{key}.json"
            res = None if force else _load(name)
            if res is None:
                print(f"  Sum m_nu free: {label}")
                prob = ft.Problem(label, likes, _free(likes, _extra(key) + ("mnu",)),
                                  nu_model="degenerate")
                start = _start(key, out["lcdm"][key]["best"])
                start["mnu"] = 0.06
                fit = ft.fit(prob, start, pool)
                prof = ft.profile(prob, fit, "mnu", _mnu_grid(key), pool)
                res = {"fit": fit, "profile": prof}
                _save(name, res)
            out["mnu"][key] = summarise_mnu(res, label, key)
    return out


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------


def profile_with_best(res: dict, name: str) -> tuple[np.ndarray, np.ndarray]:
    """The profile's grid points plus the full fit itself (the profile's minimum).

    The grid alone is too coarse to locate the minimum: a shape-preserving
    interpolation through it can only place the minimum on a grid node.
    """
    fit, prof = res["fit"], res["profile"]
    pts = [(p["value"], p["chi2"]) for p in prof["points"]] + [(fit["best"][name], fit["chi2"])]
    pts.sort()
    x, c = [], []
    for v, ch in pts:
        if x and abs(v - x[-1]) < 1e-9:
            c[-1] = min(c[-1], ch)
        else:
            x.append(v)
            c.append(ch)
    c = np.array(c)
    return np.array(x), c - c.min()


def summarise_neff(res: dict, label: str, key: str = "") -> dict:
    fit, prof = res["fit"], res["profile"]
    x, y = profile_with_best(res, "nnu")
    x, y = x.tolist(), y.tolist()
    par = ft.parabola_fit(x, y)
    iv = ft.gaussian_interval(x, y)
    ref = COMPARE.get(key, {}).get("neff")
    return {"label": label, "laplace": {"mean": fit["best"]["nnu"], "sigma": fit["errors"]["nnu"]},
            "profile": {"values": x, "dchi2": y}, "parabola": par, "interval68": iv,
            "chi2_min": fit["chi2"], "H0": fit["derived"]["H0"],
            "H0_at": {f"{p['value']:.3f}": p["derived"]["H0"] for p in prof["points"]},
            "compare": {"label": ref, "value": cfg.NEFF_PUBLISHED.get(ref)} if ref else None}


def summarise_mnu(res: dict, label: str, key: str = "") -> dict:
    fit, prof = res["fit"], res["profile"]
    x, y = profile_with_best(res, "mnu")
    par = ft.parabola_fit(x, y, max_dchi2=9.0)
    out = {"label": label, "profile": {"values": x.tolist(), "dchi2": y.tolist()},
           "parabola": par, "best_fit_mnu": fit["best"]["mnu"]}
    for tag, lower in (("flat_ge_0", 0.0), ("flat_ge_NO", cfg.MIN_SUM_NO),
                       ("flat_ge_IO", cfg.MIN_SUM_IO)):
        try:
            out[f"limit95_{tag}"] = ft.bayes_upper_limit(x, y, lower=lower)
        except ValueError as err:
            out[f"limit95_{tag}"] = None
            out[f"limit95_{tag}_error"] = str(err)
        out[f"limit95_{tag}_gauss"] = ft.truncated_gaussian_limit(par["mu"], par["sigma"],
                                                                  lower=lower)
    # how strongly the data disfavour the inverted ordering's minimum relative to normal's
    f = np.interp
    out["dchi2_IO_min_vs_NO_min"] = float(f(cfg.MIN_SUM_IO, x, y) - f(cfg.MIN_SUM_NO, x, y))
    out["sigma8_at"] = {f"{p['value']:.3f}": p["derived"]["sigma8"] for p in prof["points"]}
    out["H0_at"] = {f"{p['value']:.3f}": p["derived"]["H0"] for p in prof["points"]}
    ref = COMPARE.get(key, {}).get("mnu")
    out["compare"] = {"label": ref, "value": cfg.MNU_PUBLISHED.get(ref)} if ref else None
    return out


# ---------------------------------------------------------------------------
# DESI's own chains (the answer key)
# ---------------------------------------------------------------------------


def weighted_upper_limit(x, w, level=0.95):
    order = np.argsort(x)
    c = np.cumsum(w[order])
    return float(np.interp(level * c[-1], c, x[order]))


def weighted_mean_std(x, w):
    m = float(np.average(x, weights=w))
    return m, float(np.sqrt(np.average((x - m) ** 2, weights=w)))


def answer_keys() -> dict:
    """The published limits, recomputed from the collaborations' own chains.

    DESI DR2 (Cobaya chains; 30% burn-in as in D5), Planck 2018 (CosmoMC
    chains, burn-in already removed) and ACT DR6 (Cobaya, 30% burn-in), plus
    the getdist summaries ("margestats") shipped with each run.
    """
    out = {"desi": {}, "desi_margestats": {}, "planck": {}, "planck_margestats": {}, "act": {}}
    for run in fd.DESI_CHAIN_RUNS:
        try:
            ch = fd.load_desi_chain(run)
        except FileNotFoundError:
            continue
        w = ch["weight"]
        out["desi"][run] = {"mnu_95": weighted_upper_limit(ch["mnu"], w),
                            "mnu_mean_std": weighted_mean_std(ch["mnu"], w),
                            "H0": weighted_mean_std(ch["H0"], w),
                            "sigma8": weighted_mean_std(ch["sigma8"], w),
                            "n_samples": int(w.size)}
    for model in fd.DESI_SUMMARY_MODELS:
        for path in sorted((fd.DESI_DIR / "cobaya" / model).glob("*/chain.margestats")):
            run = f"{model}/{path.parent.name}"
            try:
                ms = fd.load_desi_margestats(run)
            except (FileNotFoundError, StopIteration):
                continue
            if "mnu" in ms:
                out["desi_margestats"][run] = {"mnu": ms["mnu"], **{k: ms[k] for k in
                                                                    ("w", "wa") if k in ms}}
    for key in fd.PLANCK_GRID_RUNS:
        try:
            ms = fd.load_planck_margestats(key)
        except FileNotFoundError:
            continue
        par = "nnu" if key.startswith("nnu") else "mnu"
        out["planck_margestats"][key] = ms.get(par)
    for key in fd.PLANCK_GRID_CHAINS:
        try:
            ch = fd.load_planck_chain(key)
        except FileNotFoundError:
            continue
        w = ch["weight"]
        par = "nnu" if key.startswith("nnu") else "mnu"
        row = {"mean_std": weighted_mean_std(ch[par], w), "H0": weighted_mean_std(ch["H0*"]
               if "H0*" in ch else ch["H0"], w)}
        if par == "mnu":
            row["mnu_95"] = weighted_upper_limit(ch["mnu"], w)
        out["planck"][key] = row
    try:
        ch = fd.load_act_chain()
        w = ch["weight"]
        out["act"]["P-ACT-LB nnu"] = {"mean_std": weighted_mean_std(ch["nnu"], w),
                                      "H0": weighted_mean_std(ch["H0"], w)}
    except (FileNotFoundError, KeyError):
        pass
    return out


def desi_answer_key() -> dict:
    return answer_keys()
