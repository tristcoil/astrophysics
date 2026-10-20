"""Every analysis of the project, each cached as data/derived/results/<name>.json.

Run them all with run.py; each function here can also be called alone. The
heavy steps (sample building with Monte Carlo errors, model libraries, MCMC)
are cached, so a second run takes minutes.
"""

from __future__ import annotations

import json
import time

import numpy as np

from config import DERIVED_DIR, RANDOM_SEED
import fit as F
import forward as FW
import samples as S

RESULTS = DERIVED_DIR / "results"
RESULTS.mkdir(exist_ok=True)


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist())
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x


def cached(name: str):
    """Decorator: run once, store the returned dict as results/<name>.json."""
    def deco(fn):
        def wrapper(*args, rerun: bool = False, **kw):
            path = RESULTS / f"{name}.json"
            if path.exists() and not rerun:
                return json.loads(path.read_text())
            t0 = time.time()
            out = _jsonable(fn(*args, **kw))
            out["_seconds"] = round(time.time() - t0, 1)
            path.write_text(json.dumps(out, indent=1))
            return out
        wrapper.__name__ = fn.__name__
        wrapper.result_name = name
        return wrapper
    return deco


# --------------------------------------------------------------------------
# Samples (cached as npz: they carry Monte Carlo errors that take minutes)
# --------------------------------------------------------------------------

def sample(name: str, *, rebuild: bool = False, **kw) -> tuple[dict, list]:
    path = DERIVED_DIR / f"sample_{name}.npz"
    fpath = DERIVED_DIR / f"sample_{name}_flow.json"
    if path.exists() and fpath.exists() and not rebuild:
        z = np.load(path, allow_pickle=True)
        return {k: z[k] for k in z.files}, json.loads(fpath.read_text())
    if name.startswith("banik2024"):
        s, flow0 = S.ps23_table(False)
        sub, flow = S.banik2024(s, flow0, **kw)
    elif name.startswith("cookson2026"):
        s, flow0 = S.ps23_table(False)
        sub, flow = S.cookson2026(s, flow0, **kw)
    elif name.startswith("shifted_banik2024"):
        s, flow0 = S.ps23_table(True)
        sub, flow = S.banik2024(s, flow0, **kw)
    else:
        raise ValueError(name)
    np.savez(path, **sub)
    fpath.write_text(json.dumps(flow))
    return sub, flow


# --------------------------------------------------------------------------
# The forward-model wide binary test (Banik et al. 2024 style)
# --------------------------------------------------------------------------

def forward_fit(tag: str, sub: dict, *, binning: FW.Binning | None = None, fixed: dict | None = None,
                prior=None, mcmc: bool = True, n_steps: int = 2500, n_walkers: int = 32,
                los_kind: str = "banik", shifted=None, lib_kw: dict | None = None, quiet: bool = False,
                gauss: dict | None = None):
    """Max-likelihood at alpha = 0 and 1, alpha free, then MCMC; returns a dict."""
    binning = binning or FW.Binning()
    lib = F.library_for(tag, sub["m1"], sub["m2"], binning=binning, quiet=quiet, **(lib_kw or {}))
    counts = F.observed_counts(sub["rsky_kau"], sub["vt"], binning)
    los = FW.los_template(binning, los_kind, shifted)
    prior = F.PRIOR if prior is None else prior
    out = dict(tag=tag, n_pairs_in_region=int(counts.sum()), fixed=fixed or {})
    best = {}
    for label, fx in (("newton", {"alpha": 0.0}), ("mond", {"alpha": 1.0}), ("free", {})):
        ft = F.Fitter(lib, counts, los, prior=prior, fixed={**(fixed or {}), **fx}, gauss=gauss)
        th, lnp = ft.max_like()
        best[label] = dict(theta=dict(zip(F.NAMES, th)), lnP=lnp)
        if not quiet:
            print(f"    {tag}: {label:6s} lnP = {lnp:9.2f}  " +
                  " ".join(f"{n}={v:.3g}" for n, v in zip(F.NAMES, th)), flush=True)
    out["best"] = best
    dl = best["newton"]["lnP"] - best["mond"]["lnP"]
    out["delta_lnP_newton_minus_mond"] = dl
    out["sigma_newton_over_mond"] = float(np.sign(dl) * np.sqrt(2 * abs(dl)))
    # pixel-level comparison of the two best fits
    pn = FW.Model(lib, los).pixel_probabilities(np.array(list(best["newton"]["theta"].values())))
    pm = FW.Model(lib, los).pixel_probabilities(np.array(list(best["mond"]["theta"].values())))
    out["pixels"] = dict(counts=counts, p_newton=pn, p_mond=pm, r_edges=binning.r_edges,
                         v_edges=binning.v_edges)
    if mcmc:
        ft = F.Fitter(lib, counts, los, prior=prior, fixed=fixed or {}, gauss=gauss)
        x0 = np.array(list(best["free"]["theta"].values()))
        chain, lnp, tau, acc = ft.run_mcmc(x0, n_walkers=n_walkers, n_steps=n_steps,
                                           seed=RANDOM_SEED % 1000, burn=n_steps // 3)
        names = [F.NAMES[i] for i in ft.free]
        out["posterior"] = F.summarize(chain, names)
        out["mcmc"] = dict(n_walkers=n_walkers, n_steps=n_steps, acceptance=acc,
                           tau=dict(zip(names, tau)), names=names,
                           chain_thin=chain[:: max(1, len(chain) // 4000)])
        if not quiet:
            a = out["posterior"].get("alpha")
            if a:
                print(f"    {tag}: alpha = {a['median']:.3f} (+{a['hi'] - a['median']:.3f} "
                      f"-{a['median'] - a['lo']:.3f}), mode {a['mode']:.3f}; acceptance {acc:.2f}", flush=True)
    return out


@cached("forward_banik2024")
def forward_banik2024():
    sub, flow = sample("banik2024")
    out = forward_fit("banik2024", sub)
    out["cut_flow"] = flow
    return out


# --------------------------------------------------------------------------
# Chae's method, and injection tests through it
# --------------------------------------------------------------------------

@cached("gamma_profile")
def gamma_profile(gammas=(1.0, 1.32, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)):
    """Best Newtonian and best MOND fits to the Banik-style sample with the
    eccentricity index gamma held fixed: how much of Newton's advantage rests
    on the eccentricity law."""
    sub, _ = sample("banik2024")
    binning = FW.Binning()
    lib = F.library_for("banik2024", sub["m1"], sub["m2"], binning=binning, quiet=True)
    counts = F.observed_counts(sub["rsky_kau"], sub["vt"], binning)
    los = FW.los_template(binning, "banik", None)
    out = dict(gamma=list(gammas), newton=[], mond=[], theta_newton=[], theta_mond=[])
    for g in gammas:
        for label, a in (("newton", 0.0), ("mond", 1.0)):
            th, lnp = F.Fitter(lib, counts, los, prior=F.PRIOR, fixed={"gamma": g, "alpha": a}).max_like()
            out[label].append(lnp)
            out["theta_" + label].append(dict(zip(F.NAMES, th)))
        print(f"    gamma {g:.2f}: lnP newton {out['newton'][-1]:9.2f}  mond {out['mond'][-1]:9.2f}  "
              f"Newton - MOND {out['newton'][-1] - out['mond'][-1]:7.2f}", flush=True)
    return out


def _vr_cos(s):
    """|cos psi| between the projected separation and the perspective-corrected
    relative sky velocity of each pair (both in the plane perpendicular to the
    pair's mean line of sight)."""
    import kinematics as KIN
    v, _ = KIN.sky_velocity(s["ra1"], s["dec1"], s["pmra1"], s["pmdec1"], s["ra2"], s["dec2"], s["pmra2"],
                            s["pmdec2"], s["dist"], s["dist"], np.nan_to_num(s["rv_sys"]), perspective=True)
    n1, n2 = KIN.unit(s["ra1"], s["dec1"])[0], KIN.unit(s["ra2"], s["dec2"])[0]
    nm = n1 + n2
    nm /= np.linalg.norm(nm, axis=-1, keepdims=True)
    dr = n1 - n2
    dr = dr - np.sum(dr * nm, axis=-1, keepdims=True) * nm
    return np.abs(np.sum(v * dr, -1)) / (np.linalg.norm(v, axis=-1) * np.linalg.norm(dr, axis=-1))


@cached("vr_angle_test")
def vr_angle_test(n: int = 600_000, vmax: float = 2.0):
    """The angle psi between each pair's separation and relative velocity
    carries the eccentricity information the (r_sky, v~) forward model ignores
    (Hwang et al. 2022). For the best Newtonian and MOND fits at each fixed
    gamma (gamma_profile), simulate the population with its companions and
    the data's velocity errors, and score the observed |cos psi| distribution
    in three separation ranges (multinomial log-likelihood, conditional on
    r_sky range, v~ < vmax)."""
    sub, _ = sample("banik2024")
    gp = gamma_profile()
    cos_d = _vr_cos(sub)
    r_d, vt_d = sub["rsky_kau"], sub["vt"]
    ranges = ((2.0, 5.0), (5.0, 12.0), (12.0, 30.0))
    edges = np.linspace(0, 1, 11)
    counts = np.array([np.histogram(cos_d[(r_d > lo) & (r_d < hi) & (vt_d < vmax)], edges)[0] for lo, hi in ranges])
    noise = empirical_noise_sampler(r_d, sub["vt_err"])
    out = dict(gamma=gp["gamma"], ranges=ranges, edges=edges, counts=counts, newton=[], mond=[], p_newton=[],
               p_mond=[])
    for i, g in enumerate(gp["gamma"]):
        for label in ("newton", "mond"):
            th = np.array([gp["theta_" + label][i][k] for k in F.NAMES])
            rng = np.random.default_rng(4000 + i)
            sim = FW.simulate(th, sub["m1"], sub["m2"], n, rng, a_lo=0.5, a_hi=100.0)
            ok = sim["valid"]
            v = sim["vt_vec"][ok] + noise(sim["r_sky_kau"][ok], rng)[:, None] * rng.standard_normal((ok.sum(), 2))
            r = sim["r_sky_kau"][ok]
            vt = np.linalg.norm(v, axis=1)
            cs = np.abs(np.sum(v * sim["sep_dir"][ok], 1)) / vt
            p = []
            for lo, hi in ranges:
                h = np.histogram(cs[(r > lo) & (r < hi) & (vt < vmax)], edges)[0].astype(float) + 0.5
                h /= h.sum()
                f_los = th[5]                       # chance alignments: isotropic, |cos psi| uniform
                p.append((1 - f_los) * h + f_los / (len(edges) - 1))
            p = np.array(p)
            out[label].append(float(np.sum(counts * np.log(p))))
            out["p_" + label].append(p)
        print(f"    v-r angle test, gamma {g:.2f}: lnL newton {out['newton'][-1]:9.2f}  mond {out['mond'][-1]:9.2f}",
              flush=True)
    return out


@cached("alpha_at_fixed_gamma")
def alpha_at_fixed_gamma(gammas=(1.5, 2.0, 2.5)):
    """Best-fitting gravity dial alpha on the Banik-style sample with the
    eccentricity index held at values the v-r angles allow (maximum
    likelihood; no MCMC)."""
    sub, _ = sample("banik2024")
    binning = FW.Binning()
    lib = F.library_for("banik2024", sub["m1"], sub["m2"], binning=binning, quiet=True)
    counts = F.observed_counts(sub["rsky_kau"], sub["vt"], binning)
    los = FW.los_template(binning, "banik", None)
    out = dict(gamma=list(gammas), alpha=[], lnP=[])
    for g in gammas:
        th, lnp = F.Fitter(lib, counts, los, prior=F.PRIOR, fixed={"gamma": g}).max_like()
        out["alpha"].append(float(th[F.NAMES.index("alpha")]))
        out["lnP"].append(lnp)
        print(f"    alpha at fixed gamma {g:.2f}: alpha = {out['alpha'][-1]:.3f}  lnP = {lnp:.2f}", flush=True)
    return out


@cached("alpha_at_fixed_gamma_chae")
def alpha_at_fixed_gamma_chae(gammas=(1.5, 2.0), tag: str = "chae2023_2to30_relative"):
    """The same on Chae's sample (2-30 kAU, relative companions, its noise and
    selection modelled): best-fitting alpha with the eccentricity index held
    fixed, against the Newtonian fit at the same gamma."""
    import chae as CH
    import injection as INJ
    d = CH.chae_sample(CH.load_chae_base())
    sub = chae_for_forward(d)
    binning = FW.Binning()
    # the same library as forward_fit_chae builds (relative companions, noise, selection)
    lib = F.library_for(tag, sub["m1"], sub["m2"], binning=binning, quiet=True, cb_mode="relative",
                        noise_sampler=empirical_noise_sampler(sub["rsky_kau"], sub["sig_vt"]), vt_cap=elbadry_cap)
    counts = F.observed_counts(sub["rsky_kau"], sub["vt"], binning)
    los = FW.los_template(binning, "banik")
    acc = acceptance_map(binning, *INJ.r_acceptance_table())
    out = dict(tag=tag, gamma=list(gammas), alpha=[], lnP_free=[], lnP_newton=[])
    for g in gammas:
        thf, lf = F.Fitter(lib, counts, los, prior=F.PRIOR, fixed={"gamma": g}, acceptance=acc).max_like()
        _, ln_ = F.Fitter(lib, counts, los, prior=F.PRIOR, fixed={"gamma": g, "alpha": 0.0}, acceptance=acc).max_like()
        out["alpha"].append(float(thf[F.NAMES.index("alpha")]))
        out["lnP_free"].append(lf)
        out["lnP_newton"].append(ln_)
        print(f"    Chae sample, gamma fixed at {g:.2f}: alpha = {out['alpha'][-1]:.3f}, lnP {lf:.2f} "
              f"(alpha = 0: {ln_:.2f})", flush=True)
    return out


@cached("vr_angle_noise")
def vr_angle_noise(n: int = 600_000, vmax: float = 2.0, seeds=(4100, 4200, 4300)):
    """Monte Carlo noise of vr_angle_test's log-likelihoods: the best Newtonian
    model at gamma = 2 re-simulated with other seeds."""
    sub, _ = sample("banik2024")
    gp = gamma_profile()
    vt0 = vr_angle_test()
    i = gp["gamma"].index(2.0)
    th = np.array([gp["theta_newton"][i][k] for k in F.NAMES])
    counts = np.array(vt0["counts"])
    edges, ranges = np.array(vt0["edges"]), vt0["ranges"]
    noise = empirical_noise_sampler(sub["rsky_kau"], sub["vt_err"])
    out = []
    for sd in seeds:
        rng = np.random.default_rng(sd)
        sim = FW.simulate(th, sub["m1"], sub["m2"], n, rng, a_lo=0.5, a_hi=100.0)
        ok = sim["valid"]
        v = sim["vt_vec"][ok] + noise(sim["r_sky_kau"][ok], rng)[:, None] * rng.standard_normal((ok.sum(), 2))
        r, vt = sim["r_sky_kau"][ok], np.linalg.norm(v, axis=1)
        cs = np.abs(np.sum(v * sim["sep_dir"][ok], 1)) / vt
        p = []
        for lo, hi in ranges:
            h = np.histogram(cs[(r > lo) & (r < hi) & (vt < vmax)], edges)[0].astype(float) + 0.5
            p.append((1 - th[5]) * h / h.sum() + th[5] / (len(edges) - 1))
        out.append(float(np.sum(counts * np.log(np.array(p)))))
        print(f"    v-r angle noise check, seed {sd}: lnL {out[-1]:.2f}", flush=True)
    return dict(seeds=list(seeds), lnL=out, original=vt0["newton"][i], sd=float(np.std(out + [vt0["newton"][i]])))


@cached("injection_chae")
def injection_chae(n_mc: int = 200):
    import chae as CH
    import injection as INJ
    d = CH.chae_sample(CH.load_chae_base())
    fb = forward_banik2024()
    thN = np.array([fb["best"]["newton"]["theta"][n] for n in F.NAMES])
    out = dict(n_pairs=len(d["s"]), theta_newton=dict(zip(F.NAMES, thN)), runs={})
    out["profile_data"] = INJ.vtilde_profile(d)
    corrected = dict(mock_noise=True, mock_persp=True)
    mocks = {
        "newton_full": dict(theta=thN, noise=True, perspective=True),
        "newton_cb_only": dict(theta=thN, noise=False, perspective=False),
        "newton_noise_only": dict(theta=thN, noise=True, perspective=False),
        "newton_persp_only": dict(theta=thN, noise=False, perspective=True),
        "newton_no_cb": dict(theta=np.r_[thN[:3], 0.0, thN[4:]], noise=True, perspective=True),
        "mond_full": dict(theta=np.r_[thN[:6], 1.0], noise=True, perspective=True),
        "mond_cb_only": dict(theta=np.r_[thN[:6], 1.0], noise=False, perspective=False),
    }
    plan = [("data", None, {}), ("data", None, corrected)]
    plan += [(k, k, {}) for k in mocks]
    plan += [("newton_full", "newton_full", corrected), ("mond_full", "mond_full", corrected)]
    cache = {}
    for i, (label, mk, kw) in enumerate(plan):
        if mk is None:
            sample_ = d
        else:
            if mk not in cache:
                cfg = mocks[mk]
                cache[mk] = INJ.mock_sample(d, cfg["theta"], np.random.default_rng(1000 + i), n_pool=3_000_000,
                                            noise=cfg["noise"], perspective=cfg["perspective"])
                out.setdefault("profiles", {})[mk] = INJ.vtilde_profile(cache[mk])
            sample_ = cache[mk]
        r = INJ.chae_on(sample_, n_mc=n_mc, seed=31 + i, **kw)
        key = f"{label}|{'corrected' if kw else 'original'}"
        out["runs"][key] = r
        print(f"    {key:32s} f_multi {r['f_multi']:.3f}  g_obs/g_pred = {np.round(r['gamma_g'], 3)}  "
              f"delta = {np.round(r['delta'], 4)} +- {np.round(r['delta_sd'], 4)}", flush=True)
    return out


if __name__ == "__main__":
    r = forward_banik2024()
    print(json.dumps({k: v for k, v in r.items() if k in ("best", "posterior", "sigma_newton_over_mond")},
                     indent=1)[:3000])


# --------------------------------------------------------------------------
# The forward model on Chae's own sample
# --------------------------------------------------------------------------

def chae_for_forward(d):
    """Chae's sample in the forward model's format (his masses, his v_p)."""
    import chae as CH
    import kinematics as KIN
    vp = CH.observed_vp(d)
    mt = d["M1"] + d["M2"]
    vc = np.sqrt(KIN.GM_KMS2_AU * mt / (d["s"] * 1000))
    sig = (KIN.K_PM * np.hypot(np.hypot(d["mu1ra_err"], d["mu2ra_err"]), np.hypot(d["mu1dec_err"], d["mu2dec_err"]))
           / np.sqrt(2) * d["dM"] / 1000 / vc)
    return dict(rsky_kau=d["s"], vt=vp / vc, m1=d["M1"], m2=d["M2"], sig_vt=sig)


def empirical_noise_sampler(rsky, sig, n_bins: int = 24):
    """sigma_v~ drawn from the data pairs in the same log-separation bin."""
    edges = np.linspace(np.log10(rsky.min()) - 1e-6, np.log10(rsky.max()) + 1e-6, n_bins + 1)
    ib = np.clip(np.digitize(np.log10(rsky), edges) - 1, 0, n_bins - 1)
    groups = [np.sort(sig[ib == k]) for k in range(n_bins)]

    def sampler(r_model, rng):
        k = np.clip(np.digitize(np.log10(np.maximum(r_model, 1e-6)), edges) - 1, 0, n_bins - 1)
        out = np.empty(len(r_model))
        for b in range(n_bins):
            m = k == b
            g = groups[b] if len(groups[b]) else groups[np.argmin([abs(i - b) if len(x) else 99 for i, x in enumerate(groups)])]
            out[m] = g[rng.integers(0, len(g), m.sum())]
        return out
    return sampler


def elbadry_cap(m_est, r_sky_kau, sig_vt):
    """El-Badry et al.'s cut in v~ units: sqrt(5 Msun / M) plus 2 sigma."""
    return np.sqrt(5.0 / m_est) + 2.0 * np.sqrt(2.0) * sig_vt


def acceptance_map(binning, ls_edges, vt_edges, acc):
    rc = 0.5 * (binning.r_edges[1:] + binning.r_edges[:-1])
    vc = 0.5 * (binning.v_edges[1:] + binning.v_edges[:-1])
    il = np.clip(np.searchsorted(ls_edges, np.log10(rc)) - 1, 0, len(ls_edges) - 2)
    iv = np.clip(np.searchsorted(vt_edges, vc) - 1, 0, len(vt_edges) - 2)
    return acc[il][:, iv]


def forward_fit_chae(tag: str, *, r_lo: float = 2.0, cb_mode: str = "relative", noise: bool = True,
                     selection: bool = True, mcmc: bool = True, n_steps: int = 2500, prior=None,
                     fixed: dict | None = None, d=None, start_from: str | None = None):
    """The forward model on Chae's sample. ``start_from``: a previous result
    whose best fits seed a single, capped optimisation per gravity law (used
    for the 0.5-30 kAU library, whose likelihood is ~10x slower)."""
    import chae as CH
    import injection as INJ
    d = CH.chae_sample(CH.load_chae_base()) if d is None else d
    sub = chae_for_forward(d)
    if r_lo < 2.0:
        r_edges = np.concatenate([np.array([r_lo, 0.75, 1.0, 1.4]), FW.R_EDGES_BANIK])
        r_edges = r_edges[r_edges >= r_lo]
        binning = FW.Binning(r_edges=r_edges, a_lo_kau=r_lo / 2.0)
    else:
        binning = FW.Binning()
    lib_kw = dict(cb_mode=cb_mode)
    if cb_mode == "absolute":
        lib_kw["log_k"] = np.arange(0.0, 3.001, 0.25)     # log10 A_max / AU
    if noise:
        lib_kw["noise_sampler"] = empirical_noise_sampler(sub["rsky_kau"], sub["sig_vt"])
    if selection:
        lib_kw["vt_cap"] = elbadry_cap
    lib = F.library_for(tag, sub["m1"], sub["m2"], binning=binning, **lib_kw)
    counts = F.observed_counts(sub["rsky_kau"], sub["vt"], binning)
    los = FW.los_template(binning, "banik")
    acc = acceptance_map(binning, *INJ.r_acceptance_table()) if selection else None
    out = dict(tag=tag, n_pairs_in_region=int(counts.sum()), r_lo=r_lo, cb_mode=cb_mode, noise=noise,
               selection=selection, quick_fit_from=start_from)
    best = {}
    prior = F.PRIOR if prior is None else prior
    seed_fit = json.loads((RESULTS / f"{start_from}.json").read_text()) if start_from else None
    for label, fx in (("newton", {"alpha": 0.0}), ("mond", {"alpha": 1.0}), ("free", {})):
        ft = F.Fitter(lib, counts, los, prior=prior, fixed={**(fixed or {}), **fx}, acceptance=acc)
        if seed_fit:
            x0 = np.array([seed_fit["best"][label]["theta"][n] for n in F.NAMES])
            x0[5] = max(x0[5], 1e-4)             # move f_LOS off its boundary
            th, lnp = ft.max_like(starts=[x0[ft.free]], nm_maxiter=1500, powell_maxiter=3)
        else:
            th, lnp = ft.max_like()
        best[label] = dict(theta=dict(zip(F.NAMES, th)), lnP=lnp)
        print(f"    {tag}: {label:6s} lnP = {lnp:9.2f}  " + " ".join(f"{n}={v:.3g}" for n, v in zip(F.NAMES, th)),
              flush=True)
    out["best"] = best
    dl = best["newton"]["lnP"] - best["mond"]["lnP"]
    out["delta_lnP_newton_minus_mond"] = dl
    out["sigma_newton_over_mond"] = float(np.sign(dl) * np.sqrt(2 * abs(dl)))
    m = FW.Model(lib, los, acc)
    out["pixels"] = dict(counts=counts, p_newton=m.pixel_probabilities(np.array(list(best["newton"]["theta"].values()))),
                         p_mond=m.pixel_probabilities(np.array(list(best["mond"]["theta"].values()))),
                         r_edges=binning.r_edges, v_edges=binning.v_edges)
    if mcmc:
        ft = F.Fitter(lib, counts, los, prior=prior, fixed=fixed or {}, acceptance=acc)
        chain, lnp, tau, accf = ft.run_mcmc(np.array(list(best["free"]["theta"].values())), n_steps=n_steps,
                                            seed=RANDOM_SEED % 1000, burn=n_steps // 3)
        names = [F.NAMES[i] for i in ft.free]
        out["posterior"] = F.summarize(chain, names)
        out["mcmc"] = dict(acceptance=accf, tau=dict(zip(names, tau)), names=names,
                           chain_thin=chain[:: max(1, len(chain) // 4000)])
        a = out["posterior"]["alpha"]
        print(f"    {tag}: alpha = {a['median']:.3f} (+{a['hi'] - a['median']:.3f} -{a['median'] - a['lo']:.3f})",
              flush=True)
    return out


@cached("forward_chae2023")
def forward_chae2023():
    return forward_fit_chae("chae2023_2to30")


# --------------------------------------------------------------------------
# Median v~ against r_sky / r_M (Banik et al. 2024 fig. 9; Cookson et al. 2026)
# --------------------------------------------------------------------------
RM1_KAU = FW.RM1_AU / 1000.0


def _median_bins(x, y, edges, n_boot: int = 300, seed: int = 0):
    rng = np.random.default_rng(seed)
    med, err, xm, n = [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = (x >= lo) & (x < hi)
        yy = y[k]
        n.append(int(k.sum()))
        if k.sum() < 15:
            med.append(np.nan); err.append(np.nan); xm.append(np.sqrt(lo * hi)); continue
        med.append(float(np.median(yy)))
        bs = np.median(yy[rng.integers(0, len(yy), (n_boot, len(yy)))], axis=1)
        err.append(float(np.std(bs)))
        xm.append(float(np.exp(np.mean(np.log(x[k])))))
    return dict(x=xm, median=med, err=err, n=n)


def profile_for(sub, theta_newton, *, vt_cuts=(1.5, 2.0, 2.5, 5.0), r_range=(2.0, 30.0), n_bins: int = 8,
                n_sim: int = 1_500_000, seed: int = 5, noise_sampler=None):
    """Data and model medians of v~ in bins of r_sky / r_M (equal-count bins
    in the data). Models: the fitted Newtonian population with alpha = 0 and
    the same population with alpha = 1 (only gravity changes)."""
    x = sub["rsky_kau"] / (RM1_KAU * np.sqrt(sub["m1"] + sub["m2"]))
    inr = (sub["rsky_kau"] > r_range[0]) & (sub["rsky_kau"] < r_range[1])
    qs = np.quantile(x[inr], np.linspace(0, 1, n_bins + 1))
    qs[0] *= 0.999
    qs[-1] *= 1.001
    out = dict(edges=qs, cuts={})
    sims = {}
    for label, al in (("newton", 0.0), ("mond", 1.0)):
        rng = np.random.default_rng(seed)
        sm = FW.simulate(np.r_[np.asarray(theta_newton)[:6], al], sub["m1"], sub["m2"], n_sim, rng,
                         a_lo=max(r_range[0] / 2.0, 0.05), a_hi=100.0)
        vt = np.hypot(sm["vt_vec"][:, 0], sm["vt_vec"][:, 1])
        if noise_sampler is not None:
            sg = noise_sampler(sm["r_sky_kau"], np.random.default_rng(seed + 1))
            vt = np.hypot(sm["vt_vec"][:, 0] + sg * rng.standard_normal(n_sim),
                          sm["vt_vec"][:, 1] + sg * rng.standard_normal(n_sim))
        ok = sm["valid"] & (sm["r_sky_kau"] > r_range[0]) & (sm["r_sky_kau"] < r_range[1])
        sims[label] = (sm["r_sky_kau"][ok] / (RM1_KAU * np.sqrt(sm["m_est"][ok])), vt[ok])
    for c in vt_cuts:
        k = inr & (sub["vt"] < c)
        d = _median_bins(x[k], sub["vt"][k], qs)
        mods = {}
        for label, (xs, vs) in sims.items():
            kk = vs < c
            mods[label] = _median_bins(xs[kk], vs[kk], qs, n_boot=20)
        out["cuts"][str(c)] = dict(data=d, **mods)
    return out


@cached("median_profiles")
def median_profiles():
    fb = forward_banik2024()
    thN = np.array([fb["best"]["newton"]["theta"][n] for n in F.NAMES])
    out = {}
    sub, _ = sample("banik2024")
    out["banik2024"] = profile_for(sub, thN)
    sub_c, _ = sample("cookson2026")
    out["cookson2026"] = profile_for(sub_c, thN, r_range=(1.0, 30.0), vt_cuts=(1.5, 2.0, 2.5), n_bins=6)
    import chae as CH
    d = CH.chae_sample(CH.load_chae_base())
    sc = chae_for_forward(d)
    out["chae2023"] = profile_for(sc, thN, r_range=(0.2, 30.0), n_bins=8,
                                  noise_sampler=empirical_noise_sampler(sc["rsky_kau"], sc["sig_vt"]))
    precise = sc["sig_vt"] < 0.05
    out["chae2023_precise"] = profile_for({k: v[precise] for k, v in sc.items()}, thN, r_range=(0.2, 30.0),
                                          n_bins=6)
    return out


# --------------------------------------------------------------------------
# Which choices move the answer? Variants of the forward-model test
# --------------------------------------------------------------------------

def _variant_sample(name):
    """Samples for the variants (cached like the others)."""
    if name == "banik2024":
        return sample("banik2024")[0]
    base = sample("banik2024")[0]
    if name == "mass_1_2":
        k = (base["mtot"] > 1.0) & (base["mtot"] < 2.0)
    elif name == "mstar_gt_05":
        k = (base["m1"] > 0.5) & (base["m2"] > 0.5)
    elif name == "d_lt_150":
        k = base["dist"] < 150.0
    elif name == "ruwe_lt_125":
        k = (base["ruwe1"] < 1.25) & (base["ruwe2"] < 1.25)
    elif name == "both_rv":
        k = np.isfinite(base["radial_velocity1"]) & np.isfinite(base["radial_velocity2"])
    elif name == "no_perspective":
        sub, _ = sample("banik2024_noperspective", perspective=False)
        return sub
    else:
        raise ValueError(name)
    return {kk: (v[k] if isinstance(v, np.ndarray) and v.shape[:1] == k.shape else v) for kk, v in base.items()}


VARIANTS = {
    # name: (sample, forward_fit keyword arguments, label)
    "nominal": ("banik2024", {}, "nominal (our rebuild of Banik et al.)"),
    "fcb_03": ("banik2024", dict(fixed={"f_CB": 0.30}), "f_CB fixed at 0.30"),
    "ecc_hwang": ("banik2024", dict(gauss={"gamma": (1.32, 0.09)}), "eccentricity prior gamma = 1.32 +- 0.09 (Hwang)"),
    "noise": ("banik2024", dict(noise=True), "measurement noise modelled"),
    "absolute_cb": ("banik2024", dict(lib_kw=dict(cb_mode="absolute")), "companion orbits independent of the wide orbit"),
    "vt_lt_25": ("banik2024", dict(vmax=2.5), "only v~ < 2.5 fitted"),
    "mass_1_2": ("mass_1_2", {}, "total mass 1-2 Msun"),
    "mstar_gt_05": ("mstar_gt_05", {}, "both stars > 0.5 Msun"),
    "d_lt_150": ("d_lt_150", {}, "distance < 150 pc"),
    "ruwe_lt_125": ("ruwe_lt_125", {}, "RUWE < 1.25 for both stars"),
    "both_rv": ("both_rv", {}, "both stars have Gaia RVs"),
    "no_perspective": ("no_perspective", {}, "no perspective correction"),
}


def run_variant(name: str, n_steps: int = 1500):
    path = RESULTS / f"variant_{name}.json"
    if path.exists():
        return json.loads(path.read_text())
    sname, kw, label = VARIANTS[name]
    sub = _variant_sample(sname)
    kw = dict(kw)
    binning = FW.Binning()
    lib_kw = dict(kw.pop("lib_kw", {}))
    if kw.pop("noise", False):
        lib_kw["noise_sampler"] = empirical_noise_sampler(sub["rsky_kau"], sub["vt_err"])  # per-component sigma ~ sigma of the magnitude
    vmax = kw.pop("vmax", None)
    if vmax is not None:
        binning = FW.Binning(v_edges=FW.V_EDGES_BANIK[FW.V_EDGES_BANIK <= vmax + 1e-9])
    if lib_kw.get("cb_mode") == "absolute":
        lib_kw["log_k"] = np.arange(0.0, 3.001, 0.25)      # log10 A_max / AU
        prior = F.PRIOR.copy()
        prior[4] = [0.0, 3.0]
        kw["prior"] = prior
    tag = "banik2024" if (sname == "banik2024" and not lib_kw and vmax is None) else f"variant_{name}"
    t0 = time.time()
    out = forward_fit(tag, sub, binning=binning, n_steps=n_steps, lib_kw=lib_kw, quiet=True, **kw)
    out["label"] = label
    out["n_pairs"] = int(len(sub["vt"]))
    out["_seconds"] = round(time.time() - t0, 1)
    out.pop("pixels", None)
    out["mcmc"].pop("chain_thin", None)
    path.write_text(json.dumps(_jsonable(out), indent=1))
    a = out["posterior"]["alpha"]
    print(f"    variant {name:15s} alpha = {a['median']:+.3f} (+{a['hi'] - a['median']:.3f} -{a['median'] - a['lo']:.3f})"
          f"  Newton vs MOND {out['sigma_newton_over_mond']:+.1f} sigma  N = {out['n_pairs']}", flush=True)
    return out


# --------------------------------------------------------------------------
# The 3D test (36 pairs with precise RV differences)
# --------------------------------------------------------------------------

@cached("three_d")
def three_d(n: int = 40000):
    import threed as T
    obs = T.observables(T.load())
    v3 = np.sqrt(obs["v_par"] ** 2 + obs["v_perp"] ** 2 + obs["v_z"] ** 2)
    import kinematics as KIN
    vesc = np.sqrt(2 * KIN.GM_KMS2_AU * (obs["m1"] + obs["m2"]) / obs["r_perp_au"])
    ratio = v3 / vesc
    order = np.argsort(-ratio)
    out = dict(chae_id=obs["chae_id"], v_over_vesc=ratio, r_perp_kau=obs["r_perp_au"] / 1000, v3d=v3,
               runs={})
    for pa in ("opik", "bpl"):
        for k in (0, 2, 4, 6):
            drop = set(order[:k].tolist()) if k else None
            r = T.gamma_posterior(pa, n=n, drop=drop)
            out["runs"][f"{pa}_drop{k}"] = dict(median=r["median"], q16=r["q16"], q84=r["q84"],
                                                gammas=r["gammas"], posterior_log10=r["posterior_log10"],
                                                n_pairs=r["n_pairs"])
            print(f"    3D {pa:4s} without the {k} fastest pairs: gamma = {r['median']:.3f} "
                  f"(+{r['q84'] - r['median']:.3f} -{r['median'] - r['q16']:.3f})", flush=True)
    st = np.load(T.RAW_DIR / "chae2026_3d" / "posterior_samples_saad_ting.npz")
    out["saad_ting"] = {k: dict(median=float(np.median(st[k])), q16=float(np.percentile(st[k], 16)),
                                q84=float(np.percentile(st[k], 84))) for k in st.files}
    return out


@cached("three_d_robustness")
def three_d_robustness(n: int = 40000):
    """The 3D likelihood with one error-model choice changed at a time (broken
    power-law P(a)), the outlier-mixture version, and injection tests on
    synthetic 36-pair samples with known gamma."""
    import threed as T
    out = dict(runs={})
    for name, kw in (("systemic_velocity_free", dict(vsys_sigma=30.0)),
                     ("perspective_ignored", dict(perspective=False)),
                     ("student_t_rvs_and_pm_jitter", dict(rv_student_nu=5, pm_jitter_mas=0.05)),
                     ("vr_angle_ignored", dict(ignore_vr_angle=True))):
        r = T.gamma_posterior("bpl", n=n, **kw)
        out["runs"][name] = dict(median=r["median"], q16=r["q16"], q84=r["q84"])
        print(f"    3D robustness {name:28s} gamma = {r['median']:.3f} (+{r['q84'] - r['median']:.3f} "
              f"-{r['median'] - r['q16']:.3f})", flush=True)
    m = T.gamma_posterior_mixture("bpl", n=n)
    out["mixture"] = dict(median=m["median"], q16=m["q16"], q84=m["q84"], f_mode=m["f_mode"],
                          p_outlier=m["p_outlier"], chae_id=m["chae_id"])
    print(f"    3D robustness outlier mixture              gamma = {m['median']:.3f} (+{m['q84'] - m['median']:.3f} "
          f"-{m['median'] - m['q16']:.3f}), outlier fraction at the mode {m['f_mode']:.2f}", flush=True)
    for g in (1.0, 1.6):
        meds = T.injection_test(g)
        out[f"injection_{g:.1f}"] = dict(recovered_medians=meds, mean=float(np.mean(meds)), sd=float(np.std(meds)))
        print(f"    3D injection gamma {g:.1f}: recovered {np.round(meds, 3)} -> mean {np.mean(meds):.3f} "
              f"+- {np.std(meds):.3f}", flush=True)
    return out


@cached("three_d_saad_ting_like")
def three_d_saad_ting_like(n: int = 40000):
    """Saad & Ting's (2026) error model inside our likelihood: each pair's
    systemic velocity left free (30 km/s prior, so its perspective terms shift
    both the sky-velocity and the RV difference), Student-t RV errors and a
    proper-motion jitter, broken power-law P(a)."""
    import threed as T
    out = {}
    for name, kw in (("systemic_velocity_free_in_rv_and_sky", dict(vsys_sigma=30.0, vsys_in_rv=True)),
                     ("saad_ting_error_model", dict(vsys_sigma=30.0, vsys_in_rv=True, rv_student_nu=5,
                                                    pm_jitter_mas=0.05))):
        r = T.gamma_posterior("bpl", n=n, **kw)
        out[name] = dict(median=r["median"], q16=r["q16"], q84=r["q84"])
        print(f"    3D {name:38s} gamma = {r['median']:.3f} (+{r['q84'] - r['median']:.3f} "
              f"-{r['median'] - r['q16']:.3f})", flush=True)
    return out


# --------------------------------------------------------------------------
# Forward fits of Chae's sample
# --------------------------------------------------------------------------
CHAE_FORWARD = {
    "chae2023_2to30_relative": dict(r_lo=2.0, cb_mode="relative"),
    "chae2023_2to30_absolute": dict(r_lo=2.0, cb_mode="absolute"),
    # the 0.5-30 kAU likelihood is ~10x slower: capped fits seeded from the 2-30 kAU solution, no MCMC
    "chae2023_05to30_absolute": dict(r_lo=0.5, cb_mode="absolute", mcmc=False,
                                     start_from="forward_chae2023_2to30_absolute"),
}


def forward_chae_all(n_steps: int = 1500):
    out = {}
    for tag, kw in CHAE_FORWARD.items():
        path = RESULTS / f"forward_{tag}.json"
        if path.exists():
            out[tag] = json.loads(path.read_text())
            continue
        kw = dict(kw)
        if kw["cb_mode"] == "absolute":
            prior = F.PRIOR.copy()
            prior[4] = [0.0, 3.0]
            kw["prior"] = prior
        r = forward_fit_chae(tag, n_steps=n_steps, **kw)
        path.write_text(json.dumps(_jsonable(r), indent=1))
        out[tag] = r
    return out


# --------------------------------------------------------------------------
# Headline numbers
# --------------------------------------------------------------------------

def write_summary():
    from config import DATA_DIR
    res = {p.stem: json.loads(p.read_text()) for p in RESULTS.glob("*.json")}
    s = {}
    fb = res.get("forward_banik2024")
    if fb:
        s["forward_banik2024"] = dict(n_pairs=fb["n_pairs_in_region"], alpha=fb["posterior"]["alpha"],
                                      sigma_newton_over_mond=fb["sigma_newton_over_mond"],
                                      delta_lnP=fb["delta_lnP_newton_minus_mond"],
                                      newton_best=fb["best"]["newton"]["theta"], cut_flow=fb.get("cut_flow"))
    inj = res.get("injection_chae")
    if inj:
        s["chae_method"] = {k: dict(gamma_g=v["gamma_g"], delta=v["delta"], delta_sd=v["delta_sd"],
                                    f_multi=v["f_multi"]) for k, v in inj["runs"].items()}
    for k, v in res.items():
        if k.startswith("variant_"):
            s.setdefault("variants", {})[k[8:]] = dict(label=v["label"], alpha=v["posterior"]["alpha"],
                                                       n_pairs=v["n_pairs"],
                                                       sigma_newton_over_mond=v["sigma_newton_over_mond"])
        if k.startswith("forward_chae2023"):
            s.setdefault("forward_chae2023", {})[k] = dict(alpha=v.get("posterior", {}).get("alpha"),
                                                           sigma_newton_over_mond=v["sigma_newton_over_mond"],
                                                           n_pairs=v["n_pairs_in_region"])
    td = res.get("three_d")
    if td:
        s["three_d"] = {k: dict(median=v["median"], q16=v["q16"], q84=v["q84"], n_pairs=v["n_pairs"])
                        for k, v in td["runs"].items()}
        s["three_d"]["saad_ting"] = td["saad_ting"]
    if "three_d_robustness" in res:
        r3 = res["three_d_robustness"]
        s.setdefault("three_d", {})["robustness"] = dict(
            runs=r3["runs"], mixture={k: r3["mixture"][k] for k in ("median", "q16", "q84", "f_mode")},
            injections={k: r3[k] for k in r3 if k.startswith("injection_")})
    if "three_d_saad_ting_like" in res:
        s.setdefault("three_d", {})["saad_ting_error_model"] = res["three_d_saad_ting_like"]

    def chae_runs(runs):
        return {k: dict(gamma_g=v["gamma_g"], delta=v["delta"], delta_sd=v["delta_sd"], f_multi=v["f_multi"])
                for k, v in runs.items() if isinstance(v, dict) and "gamma_g" in v}
    for name in ("chae_variants", "injection_controls", "injection_mechanism", "injection_mechanism_opik",
                 "injection_subsets", "injection_chae_v2"):
        if name in res:
            s[name] = chae_runs(res[name].get("runs", res[name]))
    if "gamma_profile" in res:
        g = res["gamma_profile"]
        s["gamma_profile"] = dict(gamma=g["gamma"], lnP_newton=g["newton"], lnP_mond=g["mond"])
    if "vr_angle_test" in res:
        g = res["vr_angle_test"]
        s["vr_angle_test"] = dict(gamma=g["gamma"], lnL_newton=g["newton"], lnL_mond=g["mond"])
    if "profile_tests" in res:
        s["profile_tests"] = res["profile_tests"]
    (DATA_DIR / "results.json").write_text(json.dumps(_jsonable(s), indent=1))
    return s


@cached("profile_tests")
def profile_tests():
    """Cookson et al.'s model-light test on our samples: a flat line and MOND's
    step (MOND over Newton in our forward model), each with a free normalisation,
    fitted to the binned median v~; and, for Chae's sample, the ratio of the
    median v~ of wide (s/r_M > 0.5) to tight (s/r_M < 0.1) pairs by distance
    and by velocity precision, with bootstrap errors."""
    import chae as CH
    mp = median_profiles()
    out = dict(chi2={}, ratios={})
    for k, cuts in mp.items():
        if not isinstance(cuts, dict) or "cuts" not in cuts:
            continue
        for cut, c in cuts["cuts"].items():
            m, e = np.array(c["data"]["median"]), np.array(c["data"]["err"])
            step = np.array(c["mond"]["median"]) / np.array(c["newton"]["median"])
            res_ = {}
            for name, shape in (("flat", np.ones_like(m)), ("mond_step", step)):
                A_ = np.sum(m * shape / e**2) / np.sum(shape**2 / e**2)
                res_[name] = dict(chi2=float(np.sum(((m - A_ * shape) / e) ** 2)), norm=float(A_))
            res_["n_bins"] = len(m)
            out["chi2"][f"{k}|vt<{cut}"] = res_
    d = CH.chae_sample(CH.load_chae_base())
    sc = chae_for_forward(d)
    vt, sig = sc["vt"], sc["sig_vt"]
    x = d["s"] / (7.029 * np.sqrt(d["M1"] + d["M2"]))
    rng = np.random.default_rng(0)
    groups = {"all": np.ones(len(x), bool), "d < 80 pc": d["dM"] < 80, "80-140 pc": (d["dM"] > 80) & (d["dM"] < 140),
              "d > 140 pc": d["dM"] > 140, "sigma < 0.05": sig < 0.05, "sigma 0.05-0.1": (sig > 0.05) & (sig < 0.1),
              "sigma > 0.1": sig > 0.1}
    for name, msk in groups.items():
        k = msk & (vt < 2.0)
        w, t = vt[k & (x > 0.5)], vt[k & (x < 0.1)]
        bs = [np.median(rng.choice(w, len(w))) / np.median(rng.choice(t, len(t))) for _ in range(400)]
        out["ratios"][name] = dict(ratio=float(np.median(w) / np.median(t)), err=float(np.std(bs)),
                                   n_wide=int(len(w)), n_tight=int(len(t)))
    return out


@cached("chae_variants")
def chae_variants(n_mc: int = 200):
    """Chae's method on the real data with one choice changed at a time."""
    import chae as CH
    import injection as INJ
    base = CH.load_chae_base()
    d = CH.chae_sample(base)
    sc = chae_for_forward(d)
    out = {}
    plan = {
        "published": (d, {}),
        "statistical_eccentricities": (d, dict(echoice="s")),
        "precise_only": ({k: (v[sc["sig_vt"] < 0.05] if isinstance(v, np.ndarray) and v.shape[:1] == sc["sig_vt"].shape
                              else v) for k, v in d.items()}, {}),
        "within_80pc": (CH.chae_sample(base, dmax=80.0, ext_corr=False), {}),
        "noise_and_perspective_in_expectation": (d, dict(mock_noise=True, mock_persp=True)),
    }
    for i, (name, (dd, kw)) in enumerate(plan.items()):
        r = INJ.chae_on(dd, n_mc=n_mc, seed=71 + i, **kw)
        r["n_pairs"] = int(len(dd["s"]))
        out[name] = r
        print(f"    Chae method, {name:38s} N = {r['n_pairs']:6d}  g_obs/g_pred = {np.round(r['gamma_g'], 3)}  "
              f"f_multi {r['f_multi']:.2f}", flush=True)
    return out


@cached("injection_chae_v2")
def injection_chae_v2(n_mc: int = 200, source: str = "forward_chae2023_05to30_absolute"):
    """Second-generation injection: mocks drawn from the Newtonian (and MOND)
    fit of the forward model to Chae's own sample, so that the mock matches his
    data at all separations; validated by v~ quantiles before use."""
    import chae as CH
    import injection as INJ
    d = CH.chae_sample(CH.load_chae_base())
    fr = json.loads((RESULTS / f"{source}.json").read_text())
    cb_mode = fr["cb_mode"]
    out = dict(source=source, runs={}, profiles={}, profile_data=INJ.vtilde_profile(d))
    for label, th_key, alpha in (("newton", "newton", 0.0), ("mond_same_population", "newton", 1.0),
                                 ("mond_refit", "mond", 1.0)):
        th = np.array([fr["best"][th_key]["theta"][n] for n in F.NAMES])
        th[-1] = alpha
        m = INJ.mock_sample(d, th, np.random.default_rng(5000 + len(out["runs"])), n_pool=3_000_000,
                            cb_mode=cb_mode)
        out["profiles"][label] = INJ.vtilde_profile(m)
        for kw, tagkw in (({}, "original"), (dict(mock_noise=True, mock_persp=True), "corrected")):
            r = INJ.chae_on(m, n_mc=n_mc, seed=91 + len(out["runs"]), **kw)
            out["runs"][f"{label}|{tagkw}"] = r
            print(f"    v2 mock {label:22s} {tagkw:9s} g_obs/g_pred = {np.round(r['gamma_g'], 3)}  "
                  f"f_multi {r['f_multi']:.2f}", flush=True)
    return out


@cached("injection_controls")
def injection_controls(n_mc: int = 200):
    """Control mocks for Chae's pipeline, one ingredient at a time, all without
    measurement noise or perspective effects:
    c1 pure Newtonian binaries, eccentricities from Hwang et al.'s prior alpha(s)
       (the pipeline's own assumptions: it should report 1.00);
    c2 the same with superthermal eccentricities (gamma = 2.1, our fit);
    c3 hidden companions as fitted, eccentricities from Hwang's alpha(s)."""
    import chae as CH
    import hwang as HW
    import injection as INJ
    d = CH.chae_sample(CH.load_chae_base())
    fb = forward_banik2024()
    thN = np.array([fb["best"]["newton"]["theta"][n] for n in F.NAMES])
    no_cb = np.r_[thN[:3], 0.0, thN[4:]]
    plan = {"c1_pure_hwang_e": (no_cb, HW.alpha_fit), "c2_pure_superthermal_e": (no_cb, None),
            "c3_companions_hwang_e": (thN, HW.alpha_fit)}
    out = dict(runs={}, profiles={})
    for i, (k, (th, efn)) in enumerate(plan.items()):
        m = INJ.mock_sample(d, th, np.random.default_rng(7000 + i), n_pool=3_000_000, noise=False,
                            perspective=False, e_alpha_fn=efn)
        out["profiles"][k] = INJ.vtilde_profile(m)
        r = INJ.chae_on(m, n_mc=n_mc, seed=131 + i)
        out["runs"][k] = r
        print(f"    control {k:26s} g_obs/g_pred = {np.round(r['gamma_g'], 3)}  f_multi {r['f_multi']:.2f}", flush=True)
    return out


@cached("injection_mechanism")
def injection_mechanism(n_mc: int = 200):
    """Why does a pure Newtonian mock come out 'boosted'? The control mock c1
    (no companions, noise or perspective; eccentricities from Hwang et al.'s
    prior) analysed three ways: with Hwang-style per-pair eccentricities
    (Chae's standard choice), with each pair's TRUE eccentricity (known in a
    mock), and with Chae's statistical option (eccentricities drawn from the
    separation-dependent population law, ignoring each pair's v-r angle)."""
    import chae as CH
    import hwang as HW
    import injection as INJ
    d = CH.chae_sample(CH.load_chae_base())
    fb = forward_banik2024()
    thN = np.array([fb["best"]["newton"]["theta"][n] for n in F.NAMES])
    no_cb = np.r_[thN[:3], 0.0, thN[4:]]
    m = INJ.mock_sample(d, no_cb, np.random.default_rng(7000), n_pool=3_000_000, noise=False,
                        perspective=False, e_alpha_fn=HW.alpha_fit)
    out = {}
    out["hwang_individual"] = INJ.chae_on(m, n_mc=n_mc, seed=151)
    mt = dict(m)
    te = np.clip(m["mock_true_e"], 0.001, 0.999)
    mt["e"], mt["e0"], mt["e1"] = te, te, te          # the true value, no spread
    out["true_eccentricity"] = INJ.chae_on(mt, n_mc=n_mc, seed=152)
    out["statistical_eccentricity"] = INJ.chae_on(m, n_mc=n_mc, seed=153, echoice="s")
    for k, r in out.items():
        print(f"    c1 mock, {k:26s} g_obs/g_pred = {np.round(r['gamma_g'], 3)}  f_multi {r['f_multi']:.2f}", flush=True)
    return out


@cached("injection_subsets")
def injection_subsets(n_mc: int = 200):
    """Newtonian and MOND mocks (the injection_chae population) of two of the
    chae_variants subsets, through Chae's method: within 80 pc (small proper-
    motion errors) and the precise half (sigma(v~) < 0.05)."""
    import chae as CH
    import injection as INJ
    base = CH.load_chae_base()
    d = CH.chae_sample(base)
    sc = chae_for_forward(d)
    keep = sc["sig_vt"] < 0.05
    sub = lambda x: {k: (v[keep] if isinstance(v, np.ndarray) and v.shape[:1] == keep.shape else v)
                     for k, v in x.items()}
    fb = forward_banik2024()
    thN = np.array([fb["best"]["newton"]["theta"][n] for n in F.NAMES])
    d80 = CH.chae_sample(base, dmax=80.0, ext_corr=False)
    out = {}
    for i, (label, th) in enumerate((("newton", thN), ("mond", np.r_[thN[:6], 1.0]))):
        m80 = INJ.mock_sample(d80, th, np.random.default_rng(8000 + i), n_pool=3_000_000)
        out[f"within_80pc|{label}"] = INJ.chae_on(m80, n_mc=n_mc, seed=181 + i)
        m = INJ.mock_sample(d, th, np.random.default_rng(8100 + i), n_pool=3_000_000)
        out[f"precise_only|{label}"] = INJ.chae_on(sub(m), n_mc=n_mc, seed=191 + i)
    for k, r in out.items():
        print(f"    subset mock {k:22s} g_obs/g_pred = {np.round(r['gamma_g'], 3)}  f_multi {r['f_multi']:.2f}",
              flush=True)
    return out


@cached("injection_mechanism_opik")
def injection_mechanism_opik(n_mc: int = 200):
    """Two follow-ups to injection_mechanism.
    (1) The same control mock but with a scale-free (Opik, dN/da ~ 1/a)
        separation law, the population Chae's Newtonian reference implicitly
        assumes when it deprojects each observed s at a random orbital phase:
        with true eccentricities the method should then return 1.00.
    (2) The full Newtonian mock of injection_chae (same seed) analysed with
        Chae's statistical-eccentricity option, to compare with the data's
        value under that option."""
    import chae as CH
    import hwang as HW
    import injection as INJ
    d = CH.chae_sample(CH.load_chae_base())
    fb = forward_banik2024()
    thN = np.array([fb["best"]["newton"]["theta"][n] for n in F.NAMES])
    no_cb = np.r_[thN[:3], 0.0, thN[4:]]
    orig = FW.sample_broken_power_law
    FW.sample_broken_power_law = lambda rng, n, a_lo, a_hi, a_break, beta, slope_in=-1.6: \
        orig(rng, n, a_lo, a_hi, 1e9, -1.0, slope_in=-1.0)
    try:
        m = INJ.mock_sample(d, no_cb, np.random.default_rng(7001), n_pool=3_000_000, noise=False,
                            perspective=False, e_alpha_fn=HW.alpha_fit)
    finally:
        FW.sample_broken_power_law = orig
    out = {}
    mt = dict(m)
    te = np.clip(m["mock_true_e"], 0.001, 0.999)
    mt["e"], mt["e0"], mt["e1"] = te, te, te
    out["opik_true_eccentricity"] = INJ.chae_on(mt, n_mc=n_mc, seed=161)
    out["opik_hwang_individual"] = INJ.chae_on(m, n_mc=n_mc, seed=162)
    # injection_chae drew its newton_full mock with default_rng(1000 + 2)
    mf = INJ.mock_sample(d, thN, np.random.default_rng(1002), n_pool=3_000_000, noise=True, perspective=True)
    out["newton_full_statistical"] = INJ.chae_on(mf, n_mc=n_mc, seed=163, echoice="s")
    for k, r in out.items():
        print(f"    mechanism-2 {k:26s} g_obs/g_pred = {np.round(r['gamma_g'], 3)}  f_multi {r['f_multi']:.2f}",
              flush=True)
    return out


@cached("noise_inflation")
def noise_inflation(scales=(1.0, 1.3, 1.5, 2.0)):
    """Median v~ profiles of Newtonian mocks of Chae's sample whose proper-motion
    noise is the quoted Gaia error times a factor (the catalogue cuts still use
    the quoted errors), against the data's profile."""
    import chae as CH
    import injection as INJ
    d = CH.chae_sample(CH.load_chae_base())
    fb = forward_banik2024()
    thN = np.array([fb["best"]["newton"]["theta"][n] for n in F.NAMES])
    out = dict(data=INJ.vtilde_profile(d), mocks={})
    for es in scales:
        m = INJ.mock_sample(d, thN, np.random.default_rng(42), n_pool=2_500_000, err_scale=es)
        out["mocks"][str(es)] = INJ.vtilde_profile(m)
        print(f"    noise x{es}: median v~ {np.round(out['mocks'][str(es)]['median'], 3)}", flush=True)
    return out
