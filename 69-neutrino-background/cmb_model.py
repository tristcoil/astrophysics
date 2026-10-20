"""Theory: CAMB turns cosmological parameters into every prediction we compare
with data.

One call computes, for one universe:

* the lensed CMB spectra TT, TE, EE, BB (D_ℓ in µK², ℓ ≤ 3000) and the lensing
  potential spectrum C_L^φφ, for the Planck and lensing likelihoods;
* distances at the DESI redshifts and the sound horizon at the drag epoch r_d,
  for the BAO likelihood, and D_M on a fine redshift grid for supernovae;
* derived numbers: H0, Ω_m, σ8, the age, the helium fraction, Ω_ν h², ...

Neutrinos follow the conventions of the published analyses we compare with:

* ``nu_model="one"``: the Planck/DESI baseline, one massive eigenstate of
  Σm_ν (0.06 eV by default) plus massless ones, N_eff = 3.044 in total;
* ``nu_model="degenerate"``: three equal masses Σm_ν/3, as in Planck's and
  DESI's neutrino-mass runs (DESI's ``base_mnu``: ``num_massive_neutrinos: 3``);
* Σm_ν = 0: three massless neutrinos.

The helium fraction follows from BBN given ω_b and N_eff (CAMB's PArthENoPE
table, as DESI), so changing N_eff also changes the helium, as in Planck's
``base_nnu`` runs. Dark energy is a cosmological constant unless w0/wa are
given (CAMB's PPF model, as DESI).

Every result is cached on disk under a hash of its parameters, so a second
run costs nothing. :class:`TheoryPool` evaluates many universes in parallel,
one OpenMP thread per worker (CAMB would otherwise start 16 threads in each).
"""

from __future__ import annotations

import hashlib
import json
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np

import config as cfg

#: Bump to invalidate every cached theory vector (e.g. after changing settings).
THEORY_VERSION = 4

#: Highest multipole of the lensed CMB spectra we keep, and of C_L^φφ.
LMAX = 3000
#: CAMB accuracy: lensing potential accuracy 4 and a 500-multipole margin above
#: LMAX for the lensed spectra (DESI's runs used 4 and 1250 for ℓ ≤ 4000; we
#: stop at 3000). Checked in test_physics.py against a higher-accuracy run.
LENS_POTENTIAL_ACCURACY = 4
LENS_MARGIN = 500
#: Non-linear correction to the lensing potential (as DESI's runs).
HALOFIT = "mead2016"

#: Redshift grid on which D_M(z) is stored for the supernova likelihoods.
SN_Z_GRID = np.concatenate([[0.0], np.geomspace(1e-3, 2.6, 700)])

#: Default (Planck 2018 TT,TE,EE+lowE+lensing+BAO-like) starting universe.
DEFAULTS = {
    "ombh2": 0.02237, "omch2": 0.1200, "theta100": 1.04092, "logA": 3.044, "ns": 0.9649,
    "tau": 0.0544, "nnu": cfg.NEFF_SM, "mnu": 0.06, "w0": -1.0, "wa": 0.0,
}
CAMB_KEYS = ("ombh2", "omch2", "theta100", "thetastar100", "H0", "YHe", "logA", "ns", "tau", "nnu", "mnu", "w0", "wa")


def _camb():
    import camb
    return camb


def theory_key(p: dict, nu_model: str) -> str:
    """Stable hash of everything that changes CAMB's output."""
    items = {k: round(float(p[k]), 12) for k in CAMB_KEYS if k in p}
    items["nu_model"] = nu_model
    items["v"] = THEORY_VERSION
    blob = json.dumps(items, sort_keys=True).encode()
    return hashlib.sha256(blob).hexdigest()[:24]


def camb_params(p: dict, nu_model: str = "one", *, lmax: int = LMAX, want_cls: bool = True):
    """A CAMBparams object for parameter dict ``p`` (missing keys take DEFAULTS)."""
    camb = _camb()
    q = {**DEFAULTS, **p}
    mnu = float(q["mnu"])
    if mnu < 0:
        raise ValueError("CAMB cannot take a negative neutrino mass")
    if mnu == 0.0:
        nu = {"mnu": 0.0, "num_massive_neutrinos": 0}
    elif nu_model == "one":
        nu = {"mnu": mnu, "num_massive_neutrinos": 1}
    elif nu_model == "degenerate":
        nu = {"mnu": mnu, "num_massive_neutrinos": 3, "neutrino_hierarchy": "degenerate"}
    else:
        raise ValueError(nu_model)
    kw = dict(ombh2=q["ombh2"], omch2=q["omch2"], tau=q["tau"], As=np.exp(q["logA"]) * 1e-10,
              ns=q["ns"], nnu=q["nnu"], omk=0.0, halofit_version=HALOFIT, **nu)
    if "YHe" in p:  # a fixed helium fraction (illustrations only)
        kw["YHe"] = p["YHe"]
    else:  # helium from BBN given omega_b and N_eff
        kw["bbn_predictor"] = cfg.BBN_PREDICTOR
    if "H0" in p:
        kw["H0"] = p["H0"]
    elif "thetastar100" in p:  # the exact acoustic angle (illustrations only)
        kw["thetastar"] = p["thetastar100"] / 100.0
        kw["theta_H0_range"] = [20, 120]
    else:
        kw["cosmomc_theta"] = q["theta100"] / 100.0
        kw["theta_H0_range"] = [20, 120]
    if q["w0"] != -1.0 or q["wa"] != 0.0:
        kw.update(dark_energy_model="ppf", w=q["w0"], wa=q["wa"])
    pars = camb.set_params(**kw)
    pars.set_matter_power(redshifts=[0.0], kmax=2.0)  # for sigma8
    if want_cls:
        pars.set_for_lmax(lmax, lens_potential_accuracy=LENS_POTENTIAL_ACCURACY,
                          lens_output_margin=LENS_MARGIN)
    else:
        pars.WantCls = False
    return pars


def compute_theory(p: dict, nu_model: str = "one") -> dict:
    """Run CAMB once and collect every prediction (see module docstring)."""
    camb = _camb()
    pars = camb_params(p, nu_model)
    res = camb.get_results(pars)
    cl = res.get_cmb_power_spectra(pars, CMB_unit="muK", lmax=LMAX, raw_cl=False)
    tot = cl["total"]
    pp = cl["lens_potential"][:, 0]  # [L(L+1)]^2 C_L^phiphi / 2 pi
    der = res.get_derived_params()
    h = pars.H0 / 100.0
    zb = np.asarray(cfg.DESI_Z_UNIQUE, float)
    out = {
        "ell": np.arange(LMAX + 1),
        "tt": tot[:, 0], "ee": tot[:, 1], "bb": tot[:, 2], "te": tot[:, 3], "pp": pp,
        "bao_z": zb,
        "bao_dm": res.comoving_radial_distance(zb),  # flat: D_M = chi
        "bao_hz": res.hubble_parameter(zb),
        "sn_z": SN_Z_GRID,
        "sn_dm": res.comoving_radial_distance(SN_Z_GRID),
        "sn_hz": res.hubble_parameter(SN_Z_GRID),
    }
    omm = (pars.ombh2 + pars.omch2 + pars.omnuh2) / h**2
    sig8 = float(res.get_sigma8_0())
    derived = {
        "H0": pars.H0, "omegam": omm, "omegamh2": omm * h**2, "omnuh2": pars.omnuh2,
        "sigma8": sig8, "S8": sig8 * np.sqrt(omm / 0.3), "YHe": pars.YHe,
        "age": der["age"], "zstar": der["zstar"], "rstar": der["rstar"],
        "thetastar100": der["thetastar"], "zdrag": der["zdrag"], "rdrag": der["rdrag"],
        "zeq": der["zeq"], "thetad100": der["thetad"], "thetaeq100": der["thetaeq"],
        "kd": der["kd"],
        "Ase2tau": pars.InitPower.As * 1e9 * np.exp(-2.0 * pars.Reion.optical_depth),
    }
    out["derived"] = np.array([derived[k] for k in DERIVED_KEYS], float)
    return out


DERIVED_KEYS = ("H0", "omegam", "omegamh2", "omnuh2", "sigma8", "S8", "YHe", "age", "zstar",
                "rstar", "thetastar100", "zdrag", "rdrag", "zeq", "thetad100", "thetaeq100",
                "kd", "Ase2tau")


def derived_dict(theory: dict) -> dict:
    return dict(zip(DERIVED_KEYS, (float(v) for v in theory["derived"])))


# ---------------------------------------------------------------------------
# Disk cache
# ---------------------------------------------------------------------------


def _cache_path(key: str) -> Path:
    return cfg.THEORY_CACHE / key[:2] / f"{key}.npz"


def load_cached(p: dict, nu_model: str) -> dict | None:
    path = _cache_path(theory_key(p, nu_model))
    if not path.exists():
        return None
    try:
        with np.load(path, allow_pickle=False) as z:
            return {k: z[k] for k in z.files}
    except (OSError, ValueError, EOFError):  # a half-written file from a killed run
        return None


def save_cached(p: dict, nu_model: str, theory: dict) -> None:
    path = _cache_path(theory_key(p, nu_model))
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".part.npz")
    np.savez(tmp, **theory)
    tmp.replace(path)


def theory(p: dict, nu_model: str = "one") -> dict:
    """Cached single evaluation (in this process, with all CAMB threads)."""
    t = load_cached(p, nu_model)
    if t is None:
        t = compute_theory(p, nu_model)
        save_cached(p, nu_model, t)
    return t


# ---------------------------------------------------------------------------
# Parallel evaluation
# ---------------------------------------------------------------------------

_LIMITS = None


def _init_worker():
    global _LIMITS
    from threadpoolctl import threadpool_limits
    _LIMITS = threadpool_limits(1)  # one OpenMP (CAMB) and one BLAS thread per worker


def _work(args):
    p, nu_model = args
    t = load_cached(p, nu_model)
    if t is None:
        try:
            t = compute_theory(p, nu_model)
        except Exception as err:  # noqa: BLE001 - report which universe failed
            return ("error", f"{type(err).__name__}: {err} for {p}")
        save_cached(p, nu_model, t)
    return ("ok", None)


class TheoryPool:
    """Evaluate many universes at once; results come back from the disk cache."""

    def __init__(self, n_workers: int = cfg.N_WORKERS):
        self.n_workers = n_workers
        self._pool = None
        self.n_computed = 0
        self.seconds = 0.0

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        if self._pool is not None:
            self._pool.close()
            self._pool.join()
            self._pool = None

    def _ensure_pool(self):
        if self._pool is None:
            ctx = mp.get_context("fork")
            self._pool = ctx.Pool(self.n_workers, initializer=_init_worker)

    def map(self, plist: list[dict], nu_model: str = "one") -> list[dict]:
        missing = [p for p in plist if load_cached(p, nu_model) is None]
        if missing:
            # de-duplicate by key, keep order
            seen, todo = set(), []
            for p in missing:
                k = theory_key(p, nu_model)
                if k not in seen:
                    seen.add(k)
                    todo.append(p)
            t0 = time.monotonic()
            if len(todo) == 1 or self.n_workers <= 1:
                status = [_work((p, nu_model)) for p in todo]
            else:
                self._ensure_pool()
                status = self._pool.map(_work, [(p, nu_model) for p in todo], chunksize=1)
            errors = [msg for s, msg in status if s == "error"]
            if errors:
                raise RuntimeError("CAMB failed: " + "; ".join(errors[:3]))
            self.n_computed += len(todo)
            self.seconds += time.monotonic() - t0
        out = []
        for p in plist:
            t = load_cached(p, nu_model)
            if t is None:
                raise RuntimeError(f"theory missing after evaluation: {p}")
            out.append(t)
        return out
