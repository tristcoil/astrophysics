"""Weighing the universe: what the peaks measure, and a fit of LCDM to them.

The model spectra come from CAMB (Lewis, Challinor & Lasenby 2000), the
Boltzmann code the Planck team itself uses. Given six numbers, it evolves the
photons, baryons, dark matter and neutrinos from the early universe through
recombination and outputs the power spectrum they would leave on the sky.

The six numbers of the standard model ("LCDM"):

====================  ==========================================================
``ombh2``             physical baryon density, Omega_b h^2
``omch2``             physical cold-dark-matter density, Omega_c h^2
``theta100``          100 x the angular size of the sound horizon (CosmoMC's
                      approximation to it, theta_MC) -- where the peaks sit
``logA``              ln(10^10 A_s): amplitude of the primordial fluctuations
``ns``                spectral index: how that amplitude changes with scale
``tau``               optical depth to reionisation: the fraction of CMB
                      photons re-scattered once the first stars ionised the gas
====================  ==========================================================

Everything else follows from these (plus flatness and standard neutrinos):
the Hubble constant, the age, the matter fraction.

**What the fit is, and is not.** We fit Planck's binned TT spectrum
(83 bands, l = 30-2508) by least squares with independent Gaussian errors.
Planck's real likelihood also models foregrounds and calibration, uses the
correlations between bands, and treats the lowest multipoles with a
non-Gaussian likelihood. Using the published *CMB-only* spectrum takes care of
the foregrounds (this is what Planck's own "Plik-lite" likelihood does); the
other simplifications cost some accuracy in the error bars. And tau is fixed,
not fitted -- see ``config.FIT_TAU``.
"""

from __future__ import annotations

import json
import time

import numpy as np

import cmb_physics as cp
from config import (CAMB_LMAX, DERIVED_DIR, FIT_LMAX, FIT_LMIN, FIT_PARAMETERS, FIT_START,
                    FIT_TAU, PLANCK_TT_LOWE)

FIT_CACHE = DERIVED_DIR / "tt_fit.json"
RESPONSE_CACHE = DERIVED_DIR / "parameter_responses.npz"

#: Finite-difference steps for the Jacobian: a quarter of Planck's TT+lowE
#: uncertainty on each parameter. Big enough to rise above CAMB's numerical
#: noise, small enough that the model is linear over the step.
STEPS = {k: 0.25 * PLANCK_TT_LOWE[k][1] for k in FIT_PARAMETERS}


def _camb():
    import camb
    return camb


def camb_params(p: dict):
    """Build a CAMB parameter set in Planck's baseline conventions.

    Curvature is fixed by ``omk`` (default flat). Either ``theta100`` or ``H0``
    must be given: with theta the expansion rate follows from the requirement
    that the sound horizon subtend that angle, which keeps the peaks in place
    while other parameters change.
    """
    camb = _camb()
    kw = dict(ombh2=p["ombh2"], omch2=p["omch2"], tau=p.get("tau", FIT_TAU),
              As=np.exp(p["logA"]) * 1e-10, ns=p["ns"], mnu=0.06, nnu=3.046,
              omk=p.get("omk", 0.0), lmax=CAMB_LMAX, lens_potential_accuracy=1)
    if "theta100" in p:
        kw["cosmomc_theta"] = p["theta100"] / 100.0
    else:
        kw["H0"] = p["H0"]
    return camb.set_params(**kw)


def model_tt(p: dict, *, derived: bool = False):
    """Lensed TT spectrum D_l [uK^2] for l = 0..FIT_LMAX; optionally derived parameters."""
    camb = _camb()
    pars = camb_params(p)
    res = camb.get_results(pars)
    dl = res.get_cmb_power_spectra(pars, CMB_unit="muK", lmax=FIT_LMAX)["total"][:, 0]
    if not derived:
        return dl
    d = res.get_derived_params()
    om = (pars.ombh2 + pars.omch2 + pars.omnuh2) / (pars.H0 / 100) ** 2
    der = {
        "H0": pars.H0,
        "omegam": om,
        "omegal": 1.0 - om - pars.omk,
        "omegamh2": pars.ombh2 + pars.omch2 + pars.omnuh2,
        "age": d["age"],
        "zstar": d["zstar"],
        "rstar": d["rstar"],
        "thetastar100": d["thetastar"],
        "DMstar_gpc": d["DAstar"],
        "zdrag": d["zdrag"],
        "rdrag": d["rdrag"],
        "zeq": d["zeq"],
        "Ase2tau": np.exp(p["logA"]) * 1e-10 * np.exp(-2 * p.get("tau", FIT_TAU)) * 1e9,
    }
    return dl, der


def _binned_model(p: dict, *, derived: bool = False):
    out = model_tt(p, derived=derived)
    dl, der = out if derived else (out, None)
    ell = np.arange(len(dl))
    _, db = cp.planck_bin(ell, dl, lmin=FIT_LMIN, lmax=FIT_LMAX)
    return (db, der) if derived else db


# --------------------------------------------------------------------------
# Checks and illustrations
# --------------------------------------------------------------------------


def reproduce_bestfit(bestfit: dict, theory) -> dict:
    """CAMB, fed Planck's best-fit parameters, against Planck's best-fit spectrum.

    Before trusting CAMB to fit anything, check that it reproduces the curve
    the Planck team computed for the same parameters.
    """
    p = {"ombh2": bestfit["omegabh2"], "omch2": bestfit["omegach2"],
         "theta100": bestfit["theta"], "logA": bestfit["logA"], "ns": bestfit["ns"],
         "tau": bestfit["tau"]}
    dl, der = model_tt(p, derived=True)
    m = theory.ell >= 30
    rel = dl[theory.ell] / theory.tt - 1.0
    return {"rms_rel_diff": float(np.sqrt(np.mean(rel[m] ** 2))),
            "max_rel_diff": float(np.max(np.abs(rel[m]))),
            "H0_camb": der["H0"], "H0_planck": bestfit["H0"],
            "age_camb": der["age"], "age_planck": bestfit["age"]}


#: Each panel of the "what the peaks measure" figure varies one parameter away
#: from Planck's best fit, holding the rest fixed. Baryons, dark matter and tilt
#: are varied at fixed acoustic angle, so the peaks stay put and only their
#: heights respond; curvature is varied at fixed H0, because curvature acts
#: precisely by changing the angle.
RESPONSES = {
    "ombh2": {"values": (0.70, 1.30), "relative": True, "hold": "theta100"},
    "omch2": {"values": (0.70, 1.30), "relative": True, "hold": "theta100"},
    "omk": {"values": (-0.05, 0.05), "relative": False, "hold": "H0"},
    "ns": {"values": (0.90, 1.03), "relative": False, "hold": "theta100"},
}


def parameter_responses(bestfit: dict, *, force: bool = False) -> dict:
    """Model spectra with one parameter changed at a time (for the explainer figure)."""
    if RESPONSE_CACHE.exists() and not force:
        with np.load(RESPONSE_CACHE, allow_pickle=False) as z:
            return {k: z[k] for k in z.files}
    base = {"ombh2": bestfit["omegabh2"], "omch2": bestfit["omegach2"],
            "theta100": bestfit["theta"], "logA": bestfit["logA"], "ns": bestfit["ns"],
            "tau": bestfit["tau"]}
    out = {"ell": np.arange(FIT_LMAX + 1), "base": model_tt(base)}
    print("  computing parameter-response spectra with CAMB (cached afterwards)")
    for name, spec in RESPONSES.items():
        for i, v in enumerate(spec["values"]):
            p = dict(base)
            if spec["hold"] == "H0":
                p.pop("theta100")
                p["H0"] = bestfit["H0"]
            p[name] = base[name] * v if spec["relative"] else v
            out[f"{name}_{i}"] = model_tt(p)
            out[f"{name}_{i}_value"] = np.array(p[name])
    np.savez(RESPONSE_CACHE, **out)
    return out


# --------------------------------------------------------------------------
# The fit
# --------------------------------------------------------------------------


def fit_tt(binned, *, force: bool = False, max_iter: int = 8, verbose: bool = True) -> dict:
    """Fit five LCDM parameters to Planck's binned TT spectrum by Gauss-Newton.

    Gauss-Newton linearises the model around the current guess, solves the
    resulting linear least-squares problem, and repeats. Near the answer the
    CMB spectrum is very nearly linear in these parameters, so it converges in
    a few steps. The same linearisation gives the parameter covariance,
    ``(J^T J)^-1``, with J the Jacobian of the error-normalised residuals.
    """
    if FIT_CACHE.exists() and not force:
        return json.loads(FIT_CACHE.read_text())

    y, sig = binned.dl, binned.err
    names = list(FIT_PARAMETERS)
    x = np.array([FIT_START[k] for k in names], float)
    t0 = time.monotonic()

    def as_dict(v):
        return dict(zip(names, v))

    def resid(v):
        return (_binned_model(as_dict(v)) - y) / sig

    def jacobian(v):
        cols = []
        for i, k in enumerate(names):
            h = STEPS[k]
            up, dn = v.copy(), v.copy()
            up[i] += h
            dn[i] -= h
            cols.append((resid(up) - resid(dn)) / (2 * h))
        return np.column_stack(cols)

    r = resid(x)
    chi2 = float(r @ r)
    history = [chi2]
    if verbose:
        print(f"  fitting LCDM to the TT spectrum with CAMB: start chi2 = {chi2:.1f}")
    for it in range(max_iter):
        j = jacobian(x)
        step = np.linalg.solve(j.T @ j, -j.T @ r)
        lam = 1.0
        while True:  # halve the step until chi2 improves
            trial = x + lam * step
            r_trial = resid(trial)
            chi2_trial = float(r_trial @ r_trial)
            if chi2_trial <= chi2 or lam < 1 / 16:
                break
            lam /= 2
        x, r, improvement = trial, r_trial, chi2 - chi2_trial
        chi2 = chi2_trial
        history.append(chi2)
        if verbose:
            print(f"    iteration {it + 1}: chi2 = {chi2:.3f}  ({time.monotonic() - t0:.0f} s)")
        if abs(improvement) < 0.01:
            break

    j = jacobian(x)
    cov = np.linalg.inv(j.T @ j)
    err = np.sqrt(np.diag(cov))

    # Derived parameters and their errors, by propagating the covariance through
    # a Jacobian of each derived quantity with respect to the fitted ones.
    best_db, best_der = _binned_model(as_dict(x), derived=True)
    grads = {k: np.zeros(len(names)) for k in best_der}
    for i, k in enumerate(names):
        h = STEPS[k]
        up, dn = x.copy(), x.copy()
        up[i] += h
        dn[i] -= h
        _, der_up = _binned_model(as_dict(up), derived=True)
        _, der_dn = _binned_model(as_dict(dn), derived=True)
        for q in best_der:
            grads[q][i] = (der_up[q] - der_dn[q]) / (2 * h)
    derived = {q: {"value": float(v), "error": float(np.sqrt(grads[q] @ cov @ grads[q]))}
               for q, v in best_der.items()}

    corr = cov / np.outer(err, err)
    result = {
        "parameters": {k: {"value": float(v), "error": float(e)} for k, v, e in zip(names, x, err)},
        "derived": derived,
        "tau_fixed": FIT_TAU,
        "chi2": chi2,
        "n_bins": int(len(y)),
        "dof": int(len(y) - len(names)),
        "chi2_history": history,
        "correlation": {"names": names, "matrix": corr.tolist()},
        "model_binned": best_db.tolist(),
        "l_eff": binned.ell.tolist(),
        "seconds": time.monotonic() - t0,
    }
    FIT_CACHE.write_text(json.dumps(result, indent=2))
    return result


def compare_to_planck(fit: dict) -> list[dict]:
    """Our values against Planck's TT+lowE posterior, in units of Planck's error."""
    rows = []
    for group in ("parameters", "derived"):
        for k, v in fit[group].items():
            if k not in PLANCK_TT_LOWE:
                continue
            ref, ref_err = PLANCK_TT_LOWE[k]
            rows.append({"name": k, "ours": v["value"], "ours_err": v["error"],
                         "planck": ref, "planck_err": ref_err,
                         "shift_sigma": (v["value"] - ref) / ref_err,
                         "error_ratio": v["error"] / ref_err})
    return rows
