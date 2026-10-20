"""Every fit and derived number in the project.

The steps, in the order the README tells them:

1. **Validation** -- our chi^2 at each of DESI's published best-fit points,
   likelihood by likelihood; the compressed CMB prior recomputed from the
   Planck chains (and the neutrino question it raises); D3's own CMB prior.
2. **LCDM** -- BAO alone gives Omega_m and h*r_d; compare with the CMB and with
   each supernova sample (D1's Omega_m among them).
3. **w0waCDM** -- best fits (MAP) for every data combination and their
   Delta chi^2 against LCDM: the significance of evolving dark energy, and how
   it moves when the supernova sample (or the CMB prior) is swapped.
4. **Posteriors** -- emcee chains for the headline combinations, compared with
   DESI's published chains.

Results are cached in ``data/derived/``; ``python run.py`` drives all of it.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import os
import time
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

import background as bg
import fetch_data as fd
import fitting as F
import likelihoods as L
from config import (D1_RESULTS, D3_FIT, DERIVED_DIR, DESI_DCHI2, PLANCK_EARLY_BURN,
                    UNCALIBRATED_H0)
from early_universe import EarlyUniverse

CHAIN_DIR = DERIVED_DIR / "chains"
CHAIN_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_PATH = DERIVED_DIR / "results.json"

SN_SAMPLES = fd.SN_SAMPLES
SN_LABEL = fd.SN_LABELS
CMB_VARIANTS = ("desi", "consistent", "d3")
CMB_LABEL = {"desi": "CMB", "consistent": "CMB'", "d3": "CMB(D3)"}

# ---------------------------------------------------------------------------
# Data combinations
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Combo:
    bao: bool = True
    sn: str | None = None
    cmb: str | None = None

    @property
    def label(self) -> str:
        parts = (["DESI"] if self.bao else []) + ([CMB_LABEL[self.cmb]] if self.cmb else []) + \
            ([SN_LABEL[self.sn]] if self.sn else [])
        return "+".join(parts)

    @property
    def key(self) -> str:
        return "_".join((["desi"] if self.bao else []) + ([f"cmb-{self.cmb}"] if self.cmb else [])
                        + ([self.sn.lower()] if self.sn else []))

    def desi_label(self) -> str | None:
        """The matching label of DESI's published results (CMB 'as run' only)."""
        if self.cmb not in (None, "desi") or not self.bao or self.sn == "DES-Dovekie":
            return None
        parts = ["DESI"] + (["CMBc"] if self.cmb else []) + ([SN_LABEL[self.sn]] if self.sn else [])
        return "+".join(parts)


def all_combos() -> list[Combo]:
    combos = [Combo()]
    combos += [Combo(bao=False, sn=s) for s in SN_SAMPLES]
    combos += [Combo(sn=s) for s in SN_SAMPLES]
    combos += [Combo(bao=False, cmb=c) for c in CMB_VARIANTS]  # LCDM from the CMB prior alone
    for c in CMB_VARIANTS:
        combos += [Combo(cmb=c)] + [Combo(sn=s, cmb=c) for s in SN_SAMPLES]
    return combos


# ---------------------------------------------------------------------------
# Building likelihoods (cached per process)
# ---------------------------------------------------------------------------


@lru_cache(maxsize=None)
def bao_likelihood(release: str = "dr2") -> L.BAOLikelihood:
    return L.BAOLikelihood(**fd.load_bao(release), name=f"DESI {release.upper()}")


@lru_cache(maxsize=None)
def sn_likelihood(name: str) -> L.SNLikelihood:
    return L.SNLikelihood.from_sample(fd.load_sn(name))


@lru_cache(maxsize=None)
def early_universe() -> EarlyUniverse:
    return EarlyUniverse()


@lru_cache(maxsize=None)
def planck_early_chains() -> dict:
    return fd.load_planck_early_chains(PLANCK_EARLY_BURN)


@lru_cache(maxsize=None)
def cmb_prior(variant: str) -> L.CMBPrior:
    if variant == "desi":
        return L.cmb_prior_desi()
    if variant == "consistent":
        return L.cmb_prior_from_chains(planck_early_chains(), use_omegam_mean=False)
    if variant == "desi_from_chains":
        return L.cmb_prior_from_chains(planck_early_chains(), use_omegam_mean=True)
    if variant == "d3":
        if not D3_FIT.exists():
            raise FileNotFoundError(f"{D3_FIT} missing: run project D3 first")
        return L.cmb_prior_from_d3(json.loads(D3_FIT.read_text()), early_universe())
    raise KeyError(variant)


def problem(model: str, combo: Combo) -> F.Problem:
    return F.Problem(model,
                     bao=bao_likelihood() if combo.bao else None,
                     sn=sn_likelihood(combo.sn) if combo.sn else None,
                     cmb=cmb_prior(combo.cmb) if combo.cmb else None,
                     early=early_universe() if combo.cmb else None,
                     label=combo.label)


# ---------------------------------------------------------------------------
# 1. Validation against DESI's published best fits
# ---------------------------------------------------------------------------


def desi_point(P: F.Problem, row: dict) -> np.ndarray:
    vals = {"Om": row["omm"], "hrd": row.get("hrdrag"), "H0": row.get("H0"),
            "ombh2": row.get("ombh2"), "w0": row.get("w"), "wa": row.get("wa")}
    return np.array([[vals[n] for n in P.names]])


def validate_against_desi() -> dict:
    """Our chi^2, likelihood by likelihood, at every DESI best-fit point."""
    bestfits = fd.load_desi_bestfits()
    rows = []
    sn_of = {"Pantheon+": "PantheonPlus", "Union3": "Union3", "DESY5": "DESY5"}
    for (model, label), row in sorted(bestfits.items()):
        if "full" in label:
            continue
        sn = next((v for k, v in sn_of.items() if label.endswith(k)), None)
        combo = Combo(sn=sn, cmb="desi" if "CMBc" in label else None)
        m = "lcdm" if model == "base" else "w0wa"
        P = problem(m, combo)
        theta = desi_point(P, row)
        parts = P.chi2_parts(theta)
        rec = {"model": m, "label": label, "file": row["_file"], "params": {
            n: float(v) for n, v in zip(P.names, theta[0])}}
        for k, dk in (("bao", "chi2__BAO"), ("sn", "chi2__SN"), ("cmb", "chi2__CMB_compressed")):
            if k in parts and dk in row:
                rec[f"chi2_{k}"] = float(parts[k][0])
                rec[f"desi_chi2_{k}"] = float(row[dk])
        if "cmb" in parts:  # the neutrino question: which omega_bc does DESI compare?
            c = P.cosmology(theta)
            prior = cmb_prior("desi")
            rec["chi2_cmb_if_bcnu"] = float(prior.chi2(parts["theta_star"], c.ombh2,
                                                       c.omega_cb + bg.OMEGA_NU)[0])
            rec["theta_star"] = float(parts["theta_star"][0])
            rec["omega_bc"] = float(c.omega_cb[0])
            rec["ombh2"] = float(c.ombh2[0])
        rec["desi_chi2"] = float(row["chi2"])
        rec["chi2"] = float(parts["total"][0])
        rows.append(rec)
    # One constant relative offset in theta_* (CAMB 1.5.4 vs our 2.0.4-based value)
    # should explain all remaining CMB differences if nothing else differs.
    cmb_rows = [r for r in rows if "chi2_cmb" in r]
    prior = cmb_prior("desi")

    def mismatch(eps):
        return sum((prior.chi2(np.array([r["theta_star"] * (1 + eps)]), np.array([r["ombh2"]]),
                               np.array([r["omega_bc"]]))[0] - r["desi_chi2_cmb"]) ** 2
                   for r in cmb_rows)

    eps_grid = np.linspace(-4e-5, 2e-5, 601)
    eps = float(eps_grid[np.argmin([mismatch(e) for e in eps_grid])])
    for r in cmb_rows:
        r["chi2_cmb_theta_shifted"] = float(prior.chi2(
            np.array([r["theta_star"] * (1 + eps)]), np.array([r["ombh2"]]),
            np.array([r["omega_bc"]]))[0])
    return {"rows": rows, "theta_offset_rel": eps}


def cmb_prior_report() -> dict:
    """The three CMB priors side by side, and DESI's numbers recomputed from chains."""
    out = {}
    for v in ("desi", "desi_from_chains", "consistent", "d3"):
        try:
            p = cmb_prior(v)
        except FileNotFoundError as err:
            out[v] = {"error": str(err)}
            continue
        out[v] = {"name": p.name, "mean": p.mean.tolist(), "sigma": p.sigma.tolist(),
                  "corr": (p.cov / np.outer(p.sigma, p.sigma)).tolist(), "note": p.note}
    ch = planck_early_chains()
    w = ch["weight"]
    out["chain_means"] = {
        "ombh2+omch2": float(np.average(ch["ombh2"] + ch["omch2"], weights=w)),
        "omegamh2": float(np.average(ch["omegamh2"], weights=w)),
        "omnuh2_implied": float(np.average(ch["omegamh2"] - ch["ombh2"] - ch["omch2"], weights=w)),
        "n_rows": int(w.size), "sum_weights": float(w.sum())}
    return out


# ---------------------------------------------------------------------------
# 2-3. Best fits for every combination
# ---------------------------------------------------------------------------


MAPS_PATH = DERIVED_DIR / "maps.json"


def map_fits(combos=None, models=("lcdm", "w0wa"), verbose: bool = True,
             cache: bool = True) -> dict:
    """Best fits for every combination; each is saved as soon as it is found."""
    out = json.loads(MAPS_PATH.read_text()) if cache and MAPS_PATH.exists() else {}
    for combo in combos or all_combos():
        if combo.cmb == "d3" and not D3_FIT.exists():
            continue
        if combo.key in out:
            continue
        if not combo.bao and combo.sn is None and combo.cmb is not None:
            models_here = ("lcdm",)  # the CMB prior alone cannot constrain w0, wa
        else:
            models_here = models
        rec = {"label": combo.label, "combo": combo.__dict__}
        for m in models_here:
            P = problem(m, combo)
            starts = []
            if m == "w0wa" and "lcdm" in rec:  # also start from the LCDM best fit
                lc = rec["lcdm"]
                starts.append([lc[n] if n in lc else {"w0": -1.0, "wa": 0.0}[n]
                               for n in P.free_names("hrd" in P.names)])
            res = F.find_map(P, starts)
            rec[m] = res.as_dict()
            rec[m]["n_data"] = int((P.bao.n if P.bao else 0) + (P.sn.n if P.sn else 0)
                                   + (3 if P.cmb else 0))
            rec[m]["n_params"] = P.ndim + (1 if P.sn else 0)
        if "lcdm" in rec and "w0wa" in rec:
            d = rec["lcdm"]["chi2"] - rec["w0wa"]["chi2"]
            rec["dchi2"] = d
            rec["sigma"] = F.significance(d)
        dl = combo.desi_label()
        if dl:
            rec["desi_label"] = dl
        if verbose:
            s = f"  {combo.label:<26}"
            for m in models_here:
                s += f" {m}: chi2 = {rec[m]['chi2']:9.3f}"
            if "dchi2" in rec:
                s += f"   Δχ² = {rec['dchi2']:6.2f}  ({rec['sigma']:.2f}σ)"
            print(s, flush=True)
        out[combo.key] = rec
        if cache:
            MAPS_PATH.write_text(json.dumps(to_jsonable(out), indent=1))
    return out


def desi_map_comparison(maps: dict) -> list[dict]:
    """Our Delta chi^2 against DESI's: from the paper and from their best-fit files."""
    bestfits = fd.load_desi_bestfits()
    rows = []
    for key, rec in maps.items():
        dl = rec.get("desi_label")
        if not dl or "dchi2" not in rec:
            continue
        row = {"label": rec["label"], "desi_label": dl, "dchi2": rec["dchi2"],
               "sigma": rec["sigma"]}
        lc, ww = bestfits.get(("base", dl)), bestfits.get(("base_w_wa", dl))
        if lc and ww:
            row["desi_files_dchi2"] = lc["chi2"] - ww["chi2"]
            row["desi_files_sigma"] = F.significance(row["desi_files_dchi2"])
            row["desi_chi2_lcdm"], row["desi_chi2_w0wa"] = lc["chi2"], ww["chi2"]
            row["chi2_lcdm"], row["chi2_w0wa"] = rec["lcdm"]["chi2"], rec["w0wa"]["chi2"]
        if dl in DESI_DCHI2:
            row["paper_dchi2"], row["paper_sigma"] = -DESI_DCHI2[dl][0], DESI_DCHI2[dl][1]
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# 4. Posteriors
# ---------------------------------------------------------------------------

#: (model, combo, nsteps): what we sample. BAO-only w0wa is long and banana-
#: shaped (cut by the priors), so it gets the longest run; it is also the cheapest.
MCMC_RUNS = (
    [("lcdm", Combo(), 6000), ("w0wa", Combo(), 30000), ("lcdm", Combo(bao=False, cmb="desi"), 6000),
     ("lcdm", Combo(cmb="desi"), 6000), ("w0wa", Combo(cmb="desi"), 12000),
     ("w0wa", Combo(cmb="consistent"), 12000)]
    + [("lcdm", Combo(bao=False, sn=s), 6000) for s in SN_SAMPLES]
    + [("w0wa", Combo(bao=False, sn=s), 12000) for s in SN_SAMPLES]
    + [("w0wa", Combo(sn=s), 8000) for s in SN_SAMPLES]
    + [("w0wa", Combo(sn=s, cmb="desi"), 8000) for s in SN_SAMPLES]
    + [("w0wa", Combo(sn=s, cmb="consistent"), 8000) for s in SN_SAMPLES]
    + [("w0wa", Combo(cmb="d3"), 12000)]
    + [("w0wa", Combo(sn=s, cmb="d3"), 8000) for s in SN_SAMPLES]
)


def chain_path(model: str, combo: Combo):
    return CHAIN_DIR / f"{model}__{combo.key}.npz"


def _mcmc_worker(args):
    model, combo, nsteps, start, seed = args
    from threadpoolctl import threadpool_limits

    with threadpool_limits(int(os.environ.get("D5_BLAS_THREADS", "2"))):
        P = problem(model, combo)
        out = F.run_mcmc(P, np.asarray(start, float), nsteps=nsteps, seed=seed)
        # derived quantities stored alongside the samples
        th = out["samples"]
        keep = np.linspace(0, th.shape[0] - 1, min(th.shape[0], 60000)).astype(int)
        th = th[keep]
        parts = {}
        for i in range(0, th.shape[0], 2000):
            p = P.chi2_parts(th[i:i + 2000])
            for k, v in p.items():
                parts.setdefault(k, []).append(v)
        parts = {k: np.concatenate(v) for k, v in parts.items()}
        c = P.cosmology(th)
        derived = {"Om": c.Om, "rd": parts.get("rd", np.full(th.shape[0], np.nan)),
                   "H0": 100 * c.h if P.calibrated else np.full(th.shape[0], np.nan),
                   "omch2": c.omch2, "chi2": parts["total"]}
        if not P.calibrated and "hrd" in P.names:
            derived["hrd"] = th[:, P.names.index("hrd")]
        elif P.calibrated:
            derived["hrd"] = derived["rd"] * c.h
        np.savez(chain_path(model, combo), samples=th, names=np.array(P.names),
                 log_prob=out["log_prob"][keep], tau=out["tau"],
                 acceptance=out["acceptance"], nsteps=nsteps, nwalkers=out["nwalkers"],
                 burn=out["burn"], seconds=out["seconds"], n_eff=out["n_eff"],
                 **{f"derived_{k}": v for k, v in derived.items()})
        return model, combo.key, out["seconds"], float(np.max(out["tau"])), out["n_eff"]


def run_all_mcmc(maps: dict, *, processes: int = 6, force: bool = False) -> None:
    todo = []
    for i, (model, combo, nsteps) in enumerate(MCMC_RUNS):
        if combo.cmb == "d3" and not D3_FIT.exists():
            continue
        if chain_path(model, combo).exists() and not force:
            continue
        rec = maps[combo.key][model]
        start = [rec[n] for n in problem(model, combo).names]
        todo.append((model, combo, nsteps, start, 1000 + i))
    if not todo:
        return
    print(f"  sampling {len(todo)} posteriors with {processes} processes ...", flush=True)
    t0 = time.monotonic()
    ctx = mp.get_context("fork")
    with ctx.Pool(processes) as pool:
        for model, key, sec, tau, neff in pool.imap_unordered(_mcmc_worker, todo):
            print(f"    {model:5s} {key:<34} {sec:6.0f} s   tau_max = {tau:6.0f}"
                  f"   ~{neff:,.0f} independent samples", flush=True)
    print(f"  done in {(time.monotonic() - t0) / 60:.1f} min")


def load_chain(model: str, combo: Combo) -> dict:
    with np.load(chain_path(model, combo)) as z:
        d = {k: z[k] for k in z.files}
    d["names"] = [str(n) for n in d["names"]]
    return d


def chain_column(ch: dict, name: str) -> np.ndarray:
    if name in ch["names"]:
        return ch["samples"][:, ch["names"].index(name)]
    return ch[f"derived_{name}"]


def posterior_summaries() -> dict:
    out = {}
    for model, combo, _ in MCMC_RUNS:
        if not chain_path(model, combo).exists():
            continue
        ch = load_chain(model, combo)
        rec = {"label": combo.label, "n_eff": float(ch["n_eff"]),
               "tau_max": float(np.max(ch["tau"])), "acceptance": float(ch["acceptance"])}
        for n in ("Om", "H0", "hrd", "rd", "w0", "wa", "ombh2"):
            try:
                x = chain_column(ch, n)
            except KeyError:
                continue
            if np.all(np.isfinite(x)):
                rec[n] = F.summarize(x)
        out[f"{model}__{combo.key}"] = rec
    return out


def desi_chain_summaries() -> dict:
    """DESI's own posteriors from their published chains (burn-in 30%)."""
    out = {}
    for path in sorted((fd.DESI_DIR / "cobaya").glob("*/*")):
        try:
            ch = fd.load_desi_chain(path.parent.name, path.name)
        except FileNotFoundError:
            continue
        w = ch["weight"]
        rec = {"n_rows": int(w.size)}
        for ours, theirs in (("Om", "omegam"), ("H0", "H0"), ("hrd", "hrdrag"), ("w0", "w"),
                             ("wa", "wa"), ("rd", "rdrag"), ("ombh2", "ombh2")):
            if theirs in ch:
                rec[ours] = F.summarize(ch[theirs], w)
        if "H0" not in ch and "hrdrag" in ch:  # uncalibrated runs fix H0
            rec.pop("rd", None)
        out[f"{path.parent.name}/{path.name}"] = rec
    return out


# ---------------------------------------------------------------------------
# Derived curves for the figures
# ---------------------------------------------------------------------------


def w_band(model: str, combo: Combo, z=np.linspace(0, 2.5, 126), n: int = 4000) -> dict:
    """Posterior bands of w(z) and of the dark energy density rho_de(z)/rho_de(0)."""
    ch = load_chain(model, combo)
    rng = np.random.default_rng(7)
    idx = rng.choice(ch["samples"].shape[0], size=min(n, ch["samples"].shape[0]), replace=False)
    w0, wa = chain_column(ch, "w0")[idx], chain_column(ch, "wa")[idx]
    w = bg.w_of_z(z, w0, wa)
    rho = bg.f_de(z, w0, wa)
    q = lambda a: np.percentile(a, [2.275, 15.865, 50, 84.135, 97.725], axis=0)
    return {"z": z, "w": q(w), "rho": q(rho)}


def d1_flat_om() -> tuple[float, float] | None:
    """Project D1's own flat-LCDM Omega_m from Pantheon+ (if D1 has been run)."""
    if not D1_RESULTS.exists():
        return None
    r = json.loads(D1_RESULTS.read_text())
    try:
        s = r["flat_lcdm"]["statsys"]["Om"]
        return float(s["mean"]), float(s["std"])
    except (KeyError, TypeError):
        return None


def to_jsonable(x):
    if isinstance(x, dict):
        return {str(k): to_jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [to_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    return x


def run_core(*, processes: int = 6, force_mcmc: bool = False) -> dict:
    """Validation + all best fits + all posteriors; returns the results dict."""
    t0 = time.monotonic()
    print("Recombination grid (CAMB) and CMB priors")
    early_universe()
    priors = cmb_prior_report()
    print("Validation: our chi^2 at DESI's best-fit points")
    val = validate_against_desi()
    print("Best fits (MAP) for every data combination")
    maps = map_fits()
    comparison = desi_map_comparison(maps)
    print("Posteriors (emcee)")
    run_all_mcmc(maps, processes=processes, force=force_mcmc)
    results = {"cmb_priors": priors, "validation": val, "maps": maps,
               "desi_comparison": comparison, "posteriors": posterior_summaries(),
               "desi_posteriors": desi_chain_summaries(), "d1_flat_om": d1_flat_om(),
               "uncalibrated_H0": UNCALIBRATED_H0,
               "seconds": time.monotonic() - t0}
    RESULTS_PATH.write_text(json.dumps(to_jsonable(results), indent=1))
    return results
