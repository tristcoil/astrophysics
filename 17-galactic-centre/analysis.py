"""Every measurement in the project, from the cached data to ``results.json``.

The fits, in the order the README tells the story:

======================  ==========================  ===========================
fit                     data                        physics
======================  ==========================  ===========================
``keck_newton``         Keck positions + all RVs    Kepler ellipse (as Do+19)
``keck_gr``             same                        GR: precession, redshift,
                                                    Roemer delay
``keck_redshift``       same                        GR, redshift factor f free
``vlt_kepler``          VLT positions + VLT RVs     Kepler ellipse (as G+17)
``vlt_gr``              same                        GR
``joint_gr``            both teams, two frames      GR  <- headline M and R0
``joint_newton``        same                        Kepler ellipse
``joint_redshift``      same                        f free  <- Einstein test
``joint_precession``    same                        precession factor free
======================  ==========================  ===========================

Four of them repeat, on the same data, what each team published, so our
code is graded before it is trusted. Each fit is a best-fit search followed by
an emcee run; chains are cached in ``data/chains/`` and reused unless the model
definition changes or ``--resample`` is given.
"""

from __future__ import annotations

import hashlib
import json
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import black_hole_case as bh
import fitting as ft
import orbit_physics as op
import s2data as sd
from config import CHAIN_DIR, PUBLISHED, RESULTS_PATH

#: Bump when the likelihood or the physics changes, to invalidate cached chains.
MODEL_VERSION = "2026-09-30.1"

PHYSICS = {
    #   f (redshift), f_sp (precession), Roemer delay, free relativity factors
    "newton": ({"f": 0.0, "f_sp": 0.0}, False, []),
    "gr": ({"f": 1.0, "f_sp": 1.0}, True, []),
    "redshift": ({"f_sp": 1.0}, True, ["f"]),
    "precession": ({"f": 1.0}, True, ["f_sp"]),
}

FIT_SPECS = {
    "keck_newton": {"frames": [("keck", "correlated")], "rv": "keck", "physics": "newton",
                    "label": "Keck, Newtonian", "compare": "Newton"},
    "keck_gr": {"frames": [("keck", "correlated")], "rv": "keck", "physics": "gr",
                "label": "Keck, general relativity", "compare": "GR"},
    "keck_redshift": {"frames": [("keck", "correlated")], "rv": "keck", "physics": "redshift",
                      "label": "Keck, redshift free", "compare": "redshift"},
    "vlt_kepler": {"frames": [("vlt", "fixed")], "rv": "vlt", "physics": "newton",
                   "label": "VLT, Keplerian", "compare": "g17_s2_vlt_nopriors"},
    "vlt_gr": {"frames": [("vlt", "fixed")], "rv": "vlt", "physics": "gr",
               "label": "VLT, general relativity"},
    "joint_gr": {"frames": [("keck", "correlated"), ("vlt", "fixed")], "rv": "keck",
                 "physics": "gr", "label": "Keck + VLT, general relativity"},
    "joint_newton": {"frames": [("keck", "correlated"), ("vlt", "fixed")], "rv": "keck",
                     "physics": "newton", "label": "Keck + VLT, Newtonian"},
    "joint_redshift": {"frames": [("keck", "correlated"), ("vlt", "fixed")], "rv": "keck",
                       "physics": "redshift", "label": "Keck + VLT, redshift free"},
    "joint_precession": {"frames": [("keck", "correlated"), ("vlt", "fixed")], "rv": "keck",
                         "physics": "precession", "label": "Keck + VLT, precession free"},
}

#: Reference epochs of each team's frame offsets: J2000 for Keck (as Do et al.),
#: 2009.0 for the VLT (as Gillessen et al.), so the offsets compare directly.
T_REF = {"keck": 2000.0, "vlt": 2009.0}

#: Where the searches start: roughly Do et al.'s Newtonian solution, with every
#: frame at zero. The optimiser moves well away from here when it needs to.
START = {"gm6": 4.15, "r0": 8.1, "period": 16.045, "t0": 2018.38, "e": 0.885, "inc": 133.8,
         "omega": 65.3, "Omega": 226.1, "f": 1.0, "f_sp": 1.0, "vz0": 0.0,
         "nirc2_offset": 70.0, "corr_length_keck": 25.0, "corr_mix_keck": 0.5}
START_FRAME = {"keck": (-0.4, -0.8, 0.02, 0.22), "vlt": (1.0, -1.9, 0.09, 0.0)}

MCMC_SETTINGS = {"n_burn": 5_000, "n_steps": 20_000, "thin": 20}


# --------------------------------------------------------------------------
# Building the models
# --------------------------------------------------------------------------


def load_data() -> dict:
    return {
        "keck_astro": sd.load_keck_astrometry(),
        "keck_rv": sd.load_keck_rv(),
        "vlt_astro": sd.load_vlt_astrometry(),
        "vlt_rv": sd.load_vlt_rv(),
    }


def build_model(name: str, data: dict | None = None) -> ft.Model:
    spec = FIT_SPECS[name]
    data = data or load_data()
    frames = []
    free = ["gm6", "r0", "period", "t0", "e", "inc", "omega", "Omega"]
    fixed, roemer, extra = PHYSICS[spec["physics"]]
    free += extra
    for key, noise in spec["frames"]:
        astro = data[f"{key}_astro"].used()
        frames.append(ft.Frame(key, astro, T_REF[key], noise=noise))
        free += [f"x0_{key}", f"y0_{key}", f"vx0_{key}", f"vy0_{key}"]
        if noise == "correlated":
            free += [f"corr_length_{key}", f"corr_mix_{key}"]
        elif noise == "scaled":
            free += [f"lnscale_astro_{key}"]
    rv = data[f"{spec['rv']}_rv"]
    free += ["vz0"]
    if np.any(rv.instrument == "NIRC2"):
        free += ["nirc2_offset"]
    return ft.Model(name, frames, rv, free, dict(fixed), roemer=roemer,
                    description=spec["label"])


def start_values(model: ft.Model) -> dict:
    s = dict(START)
    for fr in model.frames:
        for n, v in zip(("x0", "y0", "vx0", "vy0"), START_FRAME[fr.key]):
            s[f"{n}_{fr.key}"] = v
    return s


def _signature(name: str) -> str:
    blob = json.dumps({"spec": FIT_SPECS[name], "version": MODEL_VERSION, "mcmc": MCMC_SETTINGS,
                       "t_ref": T_REF}, sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()[:16]


# --------------------------------------------------------------------------
# Running (or loading) one fit
# --------------------------------------------------------------------------


def run_fit(name: str, *, resample: bool = False, quiet: bool = False) -> dict:
    """Best fit + MCMC for one model, cached as ``data/chains/<name>.npz``."""
    path = CHAIN_DIR / f"{name}.npz"
    sig = _signature(name)
    if path.exists() and not resample:
        z = np.load(path, allow_pickle=False)
        if str(z["signature"]) == sig:
            return {k: z[k] for k in z.files} | {"names": list(z["names"])}

    t_start = time.time()
    model = build_model(name)
    best = ft.best_fit(model, start_values(model))
    # A second search from the first answer, in case the start was far away.
    best = ft.best_fit(model, best["values"])
    mc = ft.run_mcmc(model, best["theta"], **MCMC_SETTINGS)
    # If the chain found a better point than the optimiser, polish from there.
    k = int(np.argmax(mc["log_prob"]))
    if mc["log_prob"][k] > best["log_prob"] + 1e-6:
        best = ft.best_fit(model, dict(zip(model.free, mc["chain"][k])))
    out = {
        "names": np.array(model.free),
        "chain": mc["chain"],
        "log_prob": mc["log_prob"],
        "tau": mc["tau"],
        "acceptance": mc["acceptance"],
        "n_walkers": mc["n_walkers"],
        "n_steps": mc["n_steps"],
        "best_theta": best["theta"],
        "best_log_prob": best["log_prob"],
        "signature": sig,
        "seconds": time.time() - t_start,
    }
    np.savez_compressed(path, **out)
    if not quiet:
        print(f"  {name}: {out['seconds']:.0f} s, acceptance {mc['acceptance']:.2f}, "
              f"max tau {np.nanmax(mc['tau']):.0f} steps")
    out["names"] = list(model.free)
    return out


def run_all_fits(*, resample: bool = False, processes: int = 4) -> dict:
    """Run every fit, the uncached ones in parallel worker processes."""
    todo = []
    for name in FIT_SPECS:
        path = CHAIN_DIR / f"{name}.npz"
        cached = path.exists() and str(np.load(path)["signature"]) == _signature(name)
        if resample or not cached:
            todo.append(name)
    if todo:
        print(f"Sampling {len(todo)} posterior(s): {', '.join(todo)} "
              f"(a few minutes each; cached afterwards)")
        with ProcessPoolExecutor(max_workers=min(processes, len(todo))) as pool:
            list(pool.map(_run_fit_worker, [(n, resample) for n in todo]))
    return {name: run_fit(name) for name in FIT_SPECS}


def _run_fit_worker(args):
    name, resample = args
    return run_fit(name, resample=resample)["seconds"]


# --------------------------------------------------------------------------
# Turning chains into numbers
# --------------------------------------------------------------------------


def fit_summary(name: str, fit: dict, model: ft.Model) -> dict:
    names = list(fit["names"])
    samples = fit["chain"]
    summ = ft.summarize(samples, names)
    # Frame offsets moved to J2000, the epoch Gillessen et al.'s Table 1 turns out to
    # use for them (matching their numbers shows it; the text says 2009.0 only for the
    # Plewa et al. priors).
    for key in ("vlt", "keck"):
        if f"x0_{key}" in names:
            dt = 2000.0 - T_REF[key]
            for c in ("x", "y"):
                v = samples[:, names.index(f"{c}0_{key}")] + dt * samples[:, names.index(f"v{c}0_{key}")]
                summ[f"{c}0_{key}_j2000"] = ft.summarize(v[:, None], ["v"])["v"]
    best = dict(zip(names, fit["best_theta"]))
    chi = model.chi2(fit["best_theta"][None, :])
    n_params = len(names)
    k = names.index("gm6"), names.index("r0")
    lm, lr = np.log(samples[:, k[0]]), np.log(samples[:, k[1]])
    cov = np.cov(lm, lr)
    derived = derived_quantities(samples, names, model)
    out = {
        "label": FIT_SPECS[name]["label"],
        "n_data": model.n_data,
        "n_params": n_params,
        "parameters": summ,
        "best_fit": best,
        "chi2": chi,
        "chi2_total": float(sum(v for k2, v in chi.items() if not k2.startswith("n_"))),
        "log_likelihood_max": float(fit["best_log_prob"]),
        "bic": float(-2 * fit["best_log_prob"] + n_params * np.log(model.n_data)),
        "mass_distance_exponent": float(cov[0, 1] / cov[1, 1]),
        "mass_distance_correlation": float(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])),
        "derived": derived,
        "mcmc": {"walkers": int(fit["n_walkers"]), "steps": int(fit["n_steps"]),
                 "samples": int(samples.shape[0]), "acceptance": float(fit["acceptance"]),
                 "max_autocorr_steps": float(np.nanmax(fit["tau"])),
                 "steps_per_autocorr": float(fit["n_steps"] / np.nanmax(fit["tau"]))},
    }
    return out


def derived_quantities(samples: np.ndarray, names: list[str], model: ft.Model) -> dict:
    """Physical numbers implied by each posterior sample, summarised."""
    p = {n: samples[:, k] for k, n in enumerate(names)}
    for n, v in model.fixed.items():
        p.setdefault(n, np.full(samples.shape[0], v))
    a_au = op.semi_major_axis_au(p["gm6"], p["period"])
    rp = a_au * (1 - p["e"])
    vp = op.pericentre_speed_kms(p["gm6"], p["period"], p["e"])
    rs = op.schwarzschild_radius_au(p["gm6"])
    rho = bh.density_bound_msun_pc3(p["gm6"], rp)
    quantities = {
        "mass_1e6_msun": p["gm6"],  # GM / GM_sun = M / M_sun by definition
        "semi_major_axis_au": a_au,
        "semi_major_axis_mas": a_au / p["r0"],
        "pericentre_au": rp,
        "pericentre_light_hours": rp * op.AU_KM / op.C_KMS / 3600.0,
        "pericentre_speed_kms": vp,
        "pericentre_speed_frac_c": vp / op.C_KMS,
        "apocentre_au": a_au * (1 + p["e"]),
        "schwarzschild_radius_au": rs,
        "pericentre_in_rs": rp / rs,
        "precession_arcmin_per_orbit": np.degrees(op.precession_per_orbit(p["gm6"], a_au, p["e"])) * 60,
        "redshift_at_pericentre_kms": op.relativistic_redshift_kms(p["gm6"], rp, vp),
        "density_bound_msun_pc3": rho,
        "log10_density_bound": np.log10(rho),
        "theta_g_uas": op.gravitational_angle_uas(p["gm6"], p["r0"]),
        "shadow_diameter_uas": op.shadow_diameter_uas(p["gm6"], p["r0"]),
        "mass_over_distance": p["gm6"] / p["r0"],
    }
    return {k: ft.summarize(v[:, None], [k])[k] for k, v in quantities.items()}


def compare_with_do2019(ours: dict, which: str) -> dict:
    """Our posterior against the authors' own chain on the same data and model.

    In their relativistic chains the period and pericentre angle are osculating
    elements at J2000 of a numerically integrated orbit, while ours refer to the
    2018 pericentre, so those two are only compared for the Newtonian model.
    """
    chain = sd.load_do2019_chain(which)
    w = chain["weight"]
    skip = {"period", "omega"} if which in ("GR", "redshift") else set()
    out = {}
    for n in ("gm6", "r0", "f", "period", "e", "inc", "omega", "Omega", "x0_keck", "y0_keck",
              "vx0_keck", "vy0_keck", "vz0", "nirc2_offset", "corr_length_keck", "corr_mix_keck"):
        if n not in ours["parameters"] or n in skip:
            continue
        theirs = chain[n.replace("_keck", "") if n.replace("_keck", "") in chain else
                       {"corr_length_keck": "corr_length", "corr_mix_keck": "corr_mix"}[n]]
        m = float(np.sum(w * theirs))
        s = float(np.sqrt(np.sum(w * (theirs - m) ** 2)))
        o = ours["parameters"][n]
        out[n] = {"ours_mean": o["mean"], "ours_sd": o["sd"], "theirs_mean": m, "theirs_sd": s,
                  "difference_in_sigma": (o["mean"] - m) / s if s > 0 else None,
                  "sd_ratio": o["sd"] / s if s > 0 else None}
    return out


def compare_with_published(ours: dict, key: str) -> dict:
    pub = PUBLISHED[key]
    out = {"ref": pub["ref"]}
    mapping = {"gm6": "gm6", "r0": "r0", "f": "f", "f_sp": "f_sp",
               "vz_kms": "vz0", "alpha_mas": "x0_vlt_j2000", "delta_mas": "y0_vlt_j2000",
               "valpha_uas_yr": "vx0_vlt", "vdelta_uas_yr": "vy0_vlt"}
    for k, v in pub.items():
        if k not in mapping or mapping[k] not in ours["parameters"]:
            continue
        o = ours["parameters"][mapping[k]]
        scale = 1e3 if k.endswith("uas_yr") else 1.0
        val, err = v[0], float(np.hypot(*v[1:])) if len(v) > 2 else v[1]
        out[k] = {"published": val, "published_err": err, "ours": o["median"] * scale,
                  "ours_err": 0.5 * (o["hi68"] - o["lo68"]) * scale,
                  "difference_in_sigma": (o["median"] * scale - val) / float(np.hypot(err, 0.5 * (o["hi68"] - o["lo68"]) * scale))}
    return out


# --------------------------------------------------------------------------
# The one-shot "frames mixed naively" demonstration
# --------------------------------------------------------------------------


def naive_mixing_test(data: dict) -> dict:
    """What goes wrong if both teams' positions are treated as one frame.

    Both variants use the published errors with no correlated noise, so the
    only difference is one shared frame versus one per team. Cached.
    """
    path = CHAIN_DIR / "naive_frames.json"
    if path.exists():
        cached = json.loads(path.read_text())
        if cached.get("signature") == MODEL_VERSION:
            return cached["result"]
    out = _naive_mixing_fits(data)
    path.write_text(json.dumps({"signature": MODEL_VERSION, "result": _jsonable(out)}, indent=1))
    return out


def _naive_mixing_fits(data: dict) -> dict:
    ka, va = data["keck_astro"].used(), data["vlt_astro"]
    merged = sd.Astrometry("Keck+VLT as one frame", np.concatenate([ka.t, va.t]),
                           np.concatenate([ka.x, va.x]), np.concatenate([ka.y, va.y]),
                           np.concatenate([ka.x_err, va.x_err]), np.concatenate([ka.y_err, va.y_err]),
                           np.concatenate([ka.instrument, va.instrument]))
    fixed = {"f": 1.0, "f_sp": 1.0}
    base = ["gm6", "r0", "period", "t0", "e", "inc", "omega", "Omega"]
    frame_names = lambda k: [f"x0_{k}", f"y0_{k}", f"vx0_{k}", f"vy0_{k}"]  # noqa: E731
    one = ft.Model("naive", [ft.Frame("keck", merged, 2000.0, noise="fixed")], data["keck_rv"],
                   base + frame_names("keck") + ["vz0", "nirc2_offset"], fixed)
    two = ft.Model("two_frames", [ft.Frame("keck", ka, 2000.0, noise="fixed"),
                                  ft.Frame("vlt", va, 2009.0, noise="fixed")], data["keck_rv"],
                   base + frame_names("keck") + frame_names("vlt") + ["vz0", "nirc2_offset"], fixed)
    out = {}
    for label, model in (("one_shared_frame", one), ("separate_frames", two)):
        s = start_values(model)
        best = ft.best_fit(model, s)
        best = ft.best_fit(model, best["values"])
        chi = model.chi2(best["theta"][None, :])
        chi_tot = float(sum(v for k, v in chi.items() if not k.startswith("n_")))
        out[label] = {"gm6": best["values"]["gm6"], "r0": best["values"]["r0"], "chi2": chi_tot,
                      "n_data": model.n_data, "n_params": model.ndim,
                      "chi2_pieces": chi, "best": best["values"]}
    out["delta_chi2"] = out["one_shared_frame"]["chi2"] - out["separate_frames"]["chi2"]
    b = out["separate_frames"]["best"]
    out["frame_difference_vlt_minus_keck"] = {
        "x_mas_at_2009": (b["x0_vlt"]) - (b["x0_keck"] + b["vx0_keck"] * 9.0),
        "y_mas_at_2009": (b["y0_vlt"]) - (b["y0_keck"] + b["vy0_keck"] * 9.0),
        "vx_mas_yr": b["vx0_vlt"] - b["vx0_keck"],
        "vy_mas_yr": b["vy0_vlt"] - b["vy0_keck"],
    }
    return out


# --------------------------------------------------------------------------
# Systematics: how far do M and R0 move under other reasonable choices?
# --------------------------------------------------------------------------

#: Variants of the joint GR fit. Each changes one analysis choice.
SYSTEMATICS = {
    "baseline": "as the headline fit",
    "vlt_rv_zero_point": "a separate velocity zero point for the VLT's SINFONI spectra",
    "keck_independent_errors": "Keck positions with independent errors x a fitted factor, "
                               "instead of Do et al.'s correlated-noise model",
    "vlt_errors_scaled": "VLT positions with their errors x a fitted factor",
    "no_2002_vlt_positions": "drop the ten 2002 NTT/VLT positions taken while S2 was blended near pericentre",
    "no_nirc2_velocities": "drop the six NIRC2 velocities instead of fitting their offset",
    "keck_ao_positions_only": "drop the 19 Keck speckle positions (1995-2005), keep adaptive optics only",
    "vz0_prior": "Sgr A* at rest along the line of sight: prior vz0 = 0 +- 5 km/s (as Gillessen et al. 2017)",
    "vlt_frame_priors": "the VLT frame's priors from Plewa et al. 2015: mass within 0.2 mas of the origin in "
                        "2009.0, drifting < 0.1 mas/yr (as Gillessen et al. 2017)",
    "both_gillessen_priors": "both of the above",
}

#: Plewa et al. (2015) as used by Gillessen et al. (2017): the VLT frame's origin is
#: Sgr A*'s radio position to 0.2 mas at 2009.0, and it drifts by less than 0.1 mas/yr.
VLT_FRAME_PRIORS = {"x0_vlt": (0.0, 0.2), "y0_vlt": (0.0, 0.2), "vx0_vlt": (0.0, 0.1), "vy0_vlt": (0.0, 0.1)}
VZ0_PRIOR = {"vz0": (0.0, 5.0)}


#: Each variant is applied to the GR fit (for M and R0) and to the free-redshift fit (for f).
SYSTEMATICS_BASES = ("joint_gr", "joint_redshift")


def build_variant(variant: str, data: dict, base_name: str = "joint_gr") -> ft.Model:
    base = build_model(base_name, data)
    frames, rv, free = [], base.rv, list(base.free)
    for fr in base.frames:
        astro, noise = fr.data, fr.noise
        if fr.key == "keck" and variant == "keck_independent_errors":
            noise = "scaled"
            free = [n for n in free if n not in ("corr_length_keck", "corr_mix_keck")] + ["lnscale_astro_keck"]
        if fr.key == "vlt" and variant == "vlt_errors_scaled":
            noise = "scaled"
            free.append("lnscale_astro_vlt")
        if fr.key == "vlt" and variant == "no_2002_vlt_positions":
            keep = ~((astro.t > 2002.0) & (astro.t < 2003.0))
            astro = sd.Astrometry(astro.name, astro.t[keep], astro.x[keep], astro.y[keep],
                                  astro.x_err[keep], astro.y_err[keep], astro.instrument[keep])
        if fr.key == "keck" and variant == "keck_ao_positions_only":
            keep = astro.instrument == "Keck AO"
            astro = sd.Astrometry(astro.name, astro.t[keep], astro.x[keep], astro.y[keep],
                                  astro.x_err[keep], astro.y_err[keep], astro.instrument[keep])
        frames.append(ft.Frame(fr.key, astro, fr.t_ref, noise=noise))
    if variant == "vlt_rv_zero_point":
        free.append("vlt_rv_offset")
    if variant == "no_nirc2_velocities":
        rv = rv.subset(rv.instrument != "NIRC2")
        free.remove("nirc2_offset")
    priors = {}
    if variant in ("vz0_prior", "both_gillessen_priors"):
        priors |= VZ0_PRIOR
    if variant in ("vlt_frame_priors", "both_gillessen_priors"):
        priors |= VLT_FRAME_PRIORS
    return ft.Model(f"{base_name}/{variant}", frames, rv, free, dict(base.fixed), roemer=base.roemer,
                    gaussian_priors=priors)


def _systematic_worker(job: tuple[str, str]) -> dict:
    variant, base_name = job
    data = load_data()
    model = build_variant(variant, data, base_name)
    start = start_values(model) | {"vlt_rv_offset": 0.0, "lnscale_astro_keck": 0.0, "lnscale_astro_vlt": 0.0}
    cached = CHAIN_DIR / f"{base_name}.npz"
    if cached.exists():  # start from that model's optimum when we have it
        z = np.load(cached)
        start |= dict(zip([str(n) for n in z["names"]], z["best_theta"]))
    best = ft.best_fit(model, start)
    best = ft.best_fit(model, best["values"])
    cov = ft.numerical_hessian_cov(model, best["theta"])
    sd_ = np.sqrt(np.diag(cov))
    names = model.free
    out = {"description": SYSTEMATICS[variant], "n_data": model.n_data, "log_prob": best["log_prob"]}
    for n in ("gm6", "r0", "f", "vlt_rv_offset", "lnscale_astro_keck", "lnscale_astro_vlt", "vz0",
              "nirc2_offset", "vy0_vlt"):
        if n in names:
            out[n] = float(best["values"][n])
            out[f"{n}_err"] = float(sd_[names.index(n)])
    return out


def systematics_table(*, processes: int = 4) -> dict:
    """Best fit (with Laplace errors) of every variant of both bases; cached in data/chains/."""
    path = CHAIN_DIR / "systematics.json"
    sig = hashlib.sha256((json.dumps([SYSTEMATICS, SYSTEMATICS_BASES], sort_keys=True)
                          + MODEL_VERSION).encode()).hexdigest()[:16]
    if path.exists():
        cached = json.loads(path.read_text())
        if cached.get("signature") == sig:
            return cached["variants"]
    jobs = [(v, b) for b in SYSTEMATICS_BASES for v in SYSTEMATICS]
    print(f"Re-fitting {len(jobs)} analysis variants for the systematics table (a few minutes) ...")
    with ProcessPoolExecutor(max_workers=processes) as pool:
        out = list(pool.map(_systematic_worker, jobs))
    rows = {b: {} for b in SYSTEMATICS_BASES}
    for (v, b), row in zip(jobs, out):
        rows[b][v] = row
    path.write_text(json.dumps({"signature": sig, "variants": rows}, indent=1))
    return rows


def summarize_systematics(rows: dict) -> dict:
    out = {}
    for b, table in rows.items():
        base = table["baseline"]
        keys = [k for k in ("gm6", "r0", "f") if k in base]
        shifts = {v: {f"d_{k}": row[k] - base[k] for k in keys} for v, row in table.items() if v != "baseline"}
        summary = {"shifts": shifts}
        for k in keys:
            d = np.array([sh[f"d_{k}"] for sh in shifts.values()])
            summary[f"{k}_systematic_rms"] = float(np.sqrt(np.mean(d**2)))
            summary[f"{k}_max_shift"] = float(np.max(np.abs(d)))
            summary[f"{k}_range"] = [float(min(r[k] for r in table.values())), float(max(r[k] for r in table.values()))]
        out[b] = summary
    return out


# --------------------------------------------------------------------------
# Can these data see the Schwarzschild precession? A profile likelihood.
# --------------------------------------------------------------------------

#: Precession factors at which everything else is re-optimised.
PROFILE_FSP = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 0.0, -1.0, -2.0]
AO_INSTRUMENTS = ("Keck AO", "VLT NACO")


def _ao_only(model: ft.Model) -> ft.Model:
    frames = []
    for fr in model.frames:
        a = fr.data
        keep = np.isin(a.instrument, AO_INSTRUMENTS)
        frames.append(ft.Frame(fr.key, sd.Astrometry(a.name, a.t[keep], a.x[keep], a.y[keep], a.x_err[keep],
                                                      a.y_err[keep], a.instrument[keep]), fr.t_ref, noise=fr.noise))
    return ft.Model(model.name + "/ao_only", frames, model.rv, model.free, dict(model.fixed), roemer=model.roemer)


def _profile_worker(config: str) -> dict:
    """Best log-likelihood at each fixed f_SP, warm-starting along the grid."""
    data = load_data()
    base = build_model("joint_gr", data)
    if config == "ao_only":
        base = _ao_only(base)
    z = np.load(CHAIN_DIR / "joint_gr.npz")
    start = start_values(base) | dict(zip([str(n) for n in z["names"]], z["best_theta"]))
    out, cur = {}, dict(start)
    for fsp in PROFILE_FSP:
        if fsp == 0.0:  # second leg: walk down from f_SP = 1 again
            cur = dict(out_values_at_one)
        model = ft.Model(f"profile/{config}/{fsp}", base.frames, base.rv, list(base.free),
                         {"f": 1.0, "f_sp": fsp}, roemer=True)
        best = ft.best_fit(model, cur)
        best = ft.best_fit(model, best["values"])
        cur = best["values"]
        if fsp == 1.0:
            out_values_at_one = dict(cur)
        out[fsp] = {"log_prob": best["log_prob"], "gm6": cur["gm6"], "r0": cur["r0"]}
    keys = sorted(out)
    lp = np.array([out[k]["log_prob"] for k in keys])
    return {"f_sp": keys, "log_prob": lp.tolist(), "delta_chi2": (2 * (lp.max() - lp)).tolist(),
            "gm6": [out[k]["gm6"] for k in keys], "r0": [out[k]["r0"] for k in keys],
            "n_positions": {fr.key: fr.data.n for fr in base.frames},
            "delta_chi2_gr_vs_best": float(2 * (lp.max() - out[1.0]["log_prob"])),
            "delta_chi2_newton_vs_gr": float(2 * (out[1.0]["log_prob"] - out[0.0]["log_prob"]))}


def precession_profile() -> dict:
    path = CHAIN_DIR / "precession_profile.json"
    sig = hashlib.sha256((json.dumps(PROFILE_FSP) + MODEL_VERSION).encode()).hexdigest()[:16]
    if path.exists():
        cached = json.loads(path.read_text())
        if cached.get("signature") == sig:
            return cached["profiles"]
    print("Profiling the likelihood in the precession factor (all positions, and AO-era only) ...")
    configs = ["all_positions", "ao_only"]
    with ProcessPoolExecutor(max_workers=2) as pool:
        profiles = dict(zip(configs, pool.map(_profile_worker, configs)))
    path.write_text(json.dumps({"signature": sig, "profiles": profiles}, indent=1))
    return profiles


def mass_at_distance(fit: dict, r0_target: float) -> dict:
    """The mass our posterior implies *if* the distance were ``r0_target``.

    Along the banana, ln M = a + alpha ln R0 with small scatter; this is the
    comparison that is immune to the distance disagreement between teams.
    """
    names = list(fit["names"])
    lm = np.log(fit["chain"][:, names.index("gm6")])
    lr = np.log(fit["chain"][:, names.index("r0")])
    alpha, a = np.polyfit(lr, lm, 1)
    scatter = float(np.std(lm - (a + alpha * lr)))
    m = float(np.exp(a + alpha * np.log(r0_target)))
    return {"r0": r0_target, "gm6": m, "gm6_err": m * scatter, "alpha": float(alpha)}


# --------------------------------------------------------------------------
# Kepler's third law by hand
# --------------------------------------------------------------------------


def kepler_by_hand(joint: dict) -> dict:
    """The back-of-envelope version of the whole project, from our numbers."""
    par = joint["parameters"]
    a_mas = joint["derived"]["semi_major_axis_mas"]["median"]
    period = par["period"]["median"]
    r0 = par["r0"]["median"]
    a_au = a_mas * r0
    return {
        "a_mas": a_mas, "period_yr": period, "r0_kpc": r0, "a_au": a_au,
        "mass_msun": float(op.mass_from_orbit(a_au, period)),
        "mass_if_r0_were_7_kpc": float(op.mass_from_orbit(a_mas * 7.0, period)),
        "mass_if_r0_were_9_kpc": float(op.mass_from_orbit(a_mas * 9.0, period)),
    }


def data_summary(data: dict) -> dict:
    out = {}
    for key in ("keck_astro", "vlt_astro"):
        d = data[key]
        u = d.used()
        per_inst = {}
        for inst in np.unique(u.instrument):
            m = u.instrument == inst
            per_inst[str(inst)] = {"n": int(m.sum()), "first": float(u.t[m].min()),
                                   "last": float(u.t[m].max()),
                                   "median_error_mas": float(np.median(np.hypot(u.x_err[m], u.y_err[m]) / np.sqrt(2)))}
        out[key] = {"n_total": d.n, "n_used": u.n, "excluded_epochs": d.t[d.excluded].round(3).tolist(),
                    "first": float(d.t.min()), "last": float(d.t.max()), "by_instrument": per_inst}
    for key in ("keck_rv", "vlt_rv"):
        d = data[key]
        per_inst = {}
        for inst in np.unique(d.instrument):
            m = d.instrument == inst
            per_inst[str(inst)] = {"n": int(m.sum()), "first": float(d.t[m].min()), "last": float(d.t[m].max()),
                                   "median_error_kms": float(np.median(d.v_err[m]))}
        out[key] = {"n": d.n, "first": float(d.t.min()), "last": float(d.t.max()), "by_instrument": per_inst}
    out["rv_overlap"] = sd.rv_overlap_report(data["keck_rv"], data["vlt_rv"])
    return out


# --------------------------------------------------------------------------
# Everything
# --------------------------------------------------------------------------


def run_all(*, resample: bool = False) -> tuple[dict, dict]:
    data = load_data()
    fits = run_all_fits(resample=resample)
    models = {name: build_model(name, data) for name in FIT_SPECS}
    results = {"data": data_summary(data), "fits": {}}
    for name, fit in fits.items():
        results["fits"][name] = fit_summary(name, fit, models[name])
        cmp = FIT_SPECS[name].get("compare")
        if cmp in ("Newton", "GR", "redshift"):
            results["fits"][name]["vs_do2019_chain"] = compare_with_do2019(results["fits"][name], cmp)
        elif cmp:
            results["fits"][name]["vs_published"] = compare_with_published(results["fits"][name], cmp)

    joint = results["fits"]["joint_gr"]
    results["headline"] = {
        "vs_gravity2022": compare_with_published(joint, "gravity2022"),
        "vs_gillessen2017_multistar": compare_with_published(joint, "g17_multistar"),
        "vs_do2019_gr": compare_with_published(results["fits"]["keck_gr"], "do2019_gr"),
    }
    g22 = PUBLISHED["gravity2022"]
    results["headline"]["mass_at_gravity_distance"] = {
        name: mass_at_distance(fits[name], g22["r0"][0]) | {"gravity2022_gm6": g22["gm6"][0],
                                                              "gravity2022_gm6_err": g22["gm6"][1]}
        for name in ("keck_gr", "vlt_gr", "joint_gr")}
    results["relativity"] = relativity_tests(results)
    rows = systematics_table()
    results["systematics"] = {"variants": rows, "summary": summarize_systematics(rows)}
    results["precession_profile"] = precession_profile()
    results["naive_frames"] = naive_mixing_test(data)
    results["kepler_by_hand"] = kepler_by_hand(joint)
    results["black_hole_case"] = bh.case_for_black_hole(fits["joint_gr"], list(fits["joint_gr"]["names"]))
    results["roemer_delay_days"] = roemer_range(joint)
    return results, {"data": data, "fits": fits, "models": models}


def relativity_tests(results: dict) -> dict:
    f = results["fits"]
    out = {}
    for key in ("keck_redshift", "joint_redshift"):
        par = f[key]["parameters"]["f"]
        out[key] = {"f_median": par["median"], "f_lo68": par["lo68"], "f_hi68": par["hi68"],
                    "f_sd": par["sd"], "sigma_from_newton": par["mean"] / par["sd"],
                    "sigma_from_einstein": (par["mean"] - 1.0) / par["sd"]}
    par = f["joint_precession"]["parameters"]["f_sp"]
    out["joint_precession"] = {"f_sp_median": par["median"], "f_sp_lo68": par["lo68"],
                               "f_sp_hi68": par["hi68"], "f_sp_sd": par["sd"],
                               "sigma_from_newton": par["mean"] / par["sd"],
                               "sigma_from_einstein": (par["mean"] - 1.0) / par["sd"]}
    for a, b in (("keck_gr", "keck_newton"), ("joint_gr", "joint_newton")):
        out[f"{a}_vs_{b}"] = {
            "delta_chi2": f[b]["chi2_total"] - f[a]["chi2_total"],
            "delta_log_likelihood": f[a]["log_likelihood_max"] - f[b]["log_likelihood_max"],
            "delta_bic": f[b]["bic"] - f[a]["bic"],
            "mass_shift": f[b]["parameters"]["gm6"]["median"] - f[a]["parameters"]["gm6"]["median"],
            "r0_shift": f[b]["parameters"]["r0"]["median"] - f[a]["parameters"]["r0"]["median"],
        }
    return out


def roemer_range(joint: dict) -> dict:
    b = joint["best_fit"]
    t = np.linspace(b["t0"] - 8.5, b["t0"] + 8.5, 20_001)
    z = op.orbit_state(t, b["gm6"], b["period"], b["t0"], b["e"], b["inc"], b["omega"], b["Omega"],
                       f_sp=1.0)[2]
    delay = z / op.C_AUYR * 365.25
    return {"min_days": float(delay.min()), "max_days": float(delay.max()),
            "at_pericentre_days": float(np.interp(b["t0"], t, delay))}


def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    return obj


def save_results(results: dict) -> None:
    RESULTS_PATH.write_text(json.dumps(_jsonable(results), indent=1))


def load_cached() -> tuple[dict, dict]:
    """``results.json`` plus the cached chains and models: everything the figures
    need, without re-running a single fit (for iterating on figures)."""
    res = json.loads(RESULTS_PATH.read_text())
    data = load_data()
    fits = {name: run_fit(name) for name in FIT_SPECS}
    models = {name: build_model(name, data) for name in FIT_SPECS}
    return res, {"data": data, "fits": fits, "models": models}


if __name__ == "__main__":
    res, _ = run_all()
    save_results(res)
    print(json.dumps(_jsonable(res["headline"]), indent=1))
