"""Labelled calculations with CAMB: what N_eff and Σm_ν do to the spectra.

These are not fits. Each one changes a single ingredient of a reference
universe (Planck 2018's best fit) and records what the CMB or the matter
power spectrum does, so the figures can show *why* the data constrain
neutrinos.

**N_eff in the CMB.** More radiation means a faster early expansion. Most of
that is absorbed by other parameters: the acoustic angle θ_s is measured so
precisely that we hold it fixed (CAMB adjusts H0), and the epoch of matter–
radiation equality z_eq is pinned by the peak heights, so we can hold it too
by scaling the matter densities with (1 + 0.2271 N_eff). What remains is
the damping tail: photons random-walk a distance ∝ H^(−1/2) while the sound
horizon scales as H^(−1), so a faster expansion makes the damping scale
*larger relative to* the acoustic scale and the small-scale peaks fall
(Hou et al. 2013, Sect. II.B, who hold ρ_b, z_eq and θ_s fixed as we do). Holding the damping angle fixed too (by adjusting the
helium fraction) leaves the subtle phase shift of the peaks caused by
free-streaming neutrinos (Bashinsky & Seljak 2004; Follin et al. 2015).

**Σm_ν in the matter power spectrum.** Massive neutrinos are hot dark matter:
they stream out of structures smaller than their free-streaming length and
do not fall into them, so on small scales a fraction f_ν = Ω_ν/Ω_m of the
matter does not cluster, and the rest grows more slowly. The net small-scale
suppression of the linear power spectrum at z = 0 is about 8 f_ν (Hu,
Eisenstein & Tegmark 1998, eq. 2; Lesgourgues & Pastor 2006, eq. 142; at
fixed Ω_m, h, ω_b and A_s), and non-linear growth deepens it to about 10 f_ν
at k ~ 1 h/Mpc (Brandbyge et al. 2008; Bird, Viel & Haehnelt 2012).
"""

from __future__ import annotations

import json

import numpy as np

import cmb_model as cm
import config as cfg

#: Reference universe: Planck 2018 TT,TE,EE+lowE+lensing best fit (Planck 2018 VI,
#: Table 2, first column), with N_eff = 3.044 and one 0.06 eV neutrino.
REFERENCE = {"ombh2": 0.022383, "omch2": 0.12011, "theta100": 1.040909, "logA": 3.0448,
             "ns": 0.96605, "tau": 0.0543, "nnu": cfg.NEFF_SM, "mnu": 0.06}

ZEQ_FACTOR = 0.2271  # rho_nu/rho_gamma per neutrino species: (7/8)(4/11)^(4/3)


def _cls(p: dict, nu_model: str = "one", yhe: float | None = None) -> dict:
    """Lensed D_ℓ (µK²) and derived numbers for one universe (no cache)."""
    import camb

    if yhe is not None:  # a fixed helium fraction instead of BBN's
        p = {**p, "YHe": yhe}
    pars = cm.camb_params(p, nu_model)
    res = camb.get_results(pars)
    cl = res.get_cmb_power_spectra(pars, CMB_unit="muK", lmax=cm.LMAX)
    der = res.get_derived_params()
    return {"tt": cl["total"][:, 0], "ee": cl["total"][:, 1], "te": cl["total"][:, 3],
            "unlensed_tt": cl["unlensed_scalar"][:, 0],
            "H0": pars.H0, "YHe": pars.YHe, "zeq": der["zeq"], "thetad": der["thetad"],
            "thetastar": der["thetastar"], "rstar": der["rstar"], "zstar": der["zstar"]}


def _scaled_for_zeq(neff: float, omnuh2: float) -> dict:
    """ω_c raised so z_eq is unchanged: ω_m ∝ (1 + 0.2271 N_eff), ω_b held fixed.

    Holding ω_b fixed keeps the electron density (and the baryon loading of
    the sound waves) unchanged, so what remains is the expansion rate itself.
    """
    s = (1 + ZEQ_FACTOR * neff) / (1 + ZEQ_FACTOR * cfg.NEFF_SM)
    om_ref = REFERENCE["ombh2"] + REFERENCE["omch2"] + omnuh2
    return {"omch2": s * om_ref - REFERENCE["ombh2"] - omnuh2}


def _yhe_for_thetad(neff: float, target_thetad: float, base: dict) -> float:
    """Helium fraction that restores the damping angle (bisection on CAMB runs)."""
    lo, hi = 0.15, 0.40
    for _ in range(18):
        mid = 0.5 * (lo + hi)
        td = _cls({**base, "nnu": neff}, yhe=mid)["thetad"]
        # more helium -> fewer free electrons -> longer mean free path -> larger damping angle
        if td > target_thetad:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def neff_spectra(neffs=(2.0, cfg.NEFF_SM, 4.0), *, force: bool = False) -> dict:
    """TT and EE at three N_eff values under three successively stricter controls."""
    path = cfg.DERIVED_DIR / "neff_spectra.npz"
    if path.exists() and not force:
        with np.load(path, allow_pickle=False) as z:
            return {k: z[k] for k in z.files}
    out = {"ell": np.arange(cm.LMAX + 1), "neffs": np.array(neffs)}
    ref = _cls(REFERENCE)
    omnuh2 = cm.camb_params(REFERENCE, "one", want_cls=False).omnuh2
    # Hold the exact acoustic angle theta_* (not CosmoMC's approximation theta_MC,
    # which drifts by ~0.1% when the helium fraction changes: as large as the
    # phase shift we want to see).
    base = {k: v for k, v in REFERENCE.items() if k != "theta100"}
    base["thetastar100"] = ref["thetastar"]
    for i, n in enumerate(neffs):
        # A: only the acoustic angle held (CAMB adjusts H0); helium from BBN
        a = _cls({**base, "nnu": n})
        # B: acoustic angle and matter-radiation equality held; helium held too
        b_par = {**base, **_scaled_for_zeq(n, omnuh2), "nnu": n}
        b = _cls(b_par, yhe=ref["YHe"])
        # C: ... and the damping angle held, by adjusting the helium fraction
        same = abs(n - cfg.NEFF_SM) < 1e-9
        yhe = ref["YHe"] if same else _yhe_for_thetad(n, ref["thetad"], b_par)
        c = _cls(b_par, yhe=yhe)
        for tag, s in (("theta", a), ("zeq", b), ("thetad", c)):
            out[f"{tag}_tt_{i}"] = s["tt"]
            out[f"{tag}_ee_{i}"] = s["ee"]
            out[f"{tag}_utt_{i}"] = s["unlensed_tt"]
            for q in ("H0", "YHe", "zeq", "thetad", "thetastar", "rstar"):
                out[f"{tag}_{q}_{i}"] = np.array(s[q])
    np.savez(path, **out)
    return out


def matter_power(masses=(0.0, 0.06, 0.12, 0.24, 0.48), *, force: bool = False) -> dict:
    """Linear and non-linear P(k) at z = 0 for several Σm_ν (three degenerate states).

    Two families: at fixed ω_b, ω_m = ω_b + ω_c + ω_ν, h and A_s (the textbook
    comparison behind "−8 f_ν"), and at fixed θ_s instead of h (what the CMB
    actually pins).
    """
    import camb

    path = cfg.DERIVED_DIR / "matter_power.npz"
    if path.exists() and not force:
        with np.load(path, allow_pickle=False) as z:
            return {k: z[k] for k in z.files}
    ref_pars = cm.camb_params(REFERENCE, "one", want_cls=False)
    h_ref = ref_pars.H0 / 100.0
    omegam_h2 = REFERENCE["ombh2"] + REFERENCE["omch2"] + ref_pars.omnuh2
    kh = np.geomspace(1e-4, 10.0, 400)
    out = {"kh": kh, "masses": np.array(masses), "h_ref": np.array(h_ref)}
    for fam in ("fixed_h", "fixed_theta"):
        for i, m in enumerate(masses):
            p = dict(REFERENCE)
            p["mnu"] = m
            if fam == "fixed_h":
                p.pop("theta100")
                p["H0"] = 100 * h_ref
            pars = cm.camb_params(p, "degenerate", want_cls=False)
            # keep total omega_m fixed: omega_c = omega_m - omega_b - omega_nu(m)
            pars.set_cosmology(ombh2=p["ombh2"], omch2=omegam_h2 - p["ombh2"] - pars.omnuh2,
                               mnu=m, nnu=cfg.NEFF_SM, num_massive_neutrinos=3 if m > 0 else 0,
                               neutrino_hierarchy="degenerate", tau=p["tau"],
                               **({"H0": 100 * h_ref} if fam == "fixed_h"
                                  else {"cosmomc_theta": p["theta100"] / 100}))
            pars.set_matter_power(redshifts=[0.0], kmax=20.0, nonlinear=False)
            pars.NonLinear = camb.model.NonLinear_none
            res = camb.get_results(pars)
            k_lin, _, pk_lin = res.get_matter_power_spectrum(minkh=1e-4, maxkh=10.0, npoints=400)
            pars.NonLinear = camb.model.NonLinear_pk
            pars.NonLinearModel.set_params(halofit_version="mead2020")
            res_nl = camb.get_results(pars)
            k_nl, _, pk_nl = res_nl.get_matter_power_spectrum(minkh=1e-4, maxkh=10.0, npoints=400)
            omm = (pars.ombh2 + pars.omch2 + pars.omnuh2) / (pars.H0 / 100) ** 2
            out[f"{fam}_lin_{i}"] = np.interp(kh, k_lin, pk_lin[0])
            out[f"{fam}_nl_{i}"] = np.interp(kh, k_nl, pk_nl[0])
            out[f"{fam}_fnu_{i}"] = np.array(pars.omnuh2 / (omm * (pars.H0 / 100) ** 2))
            out[f"{fam}_h_{i}"] = np.array(pars.H0 / 100)
            out[f"{fam}_sigma8_{i}"] = np.array(res.get_sigma8_0())
    np.savez(path, **out)
    return out


def save_json(obj, name: str):
    (cfg.DERIVED_DIR / name).write_text(json.dumps(obj, indent=2))
