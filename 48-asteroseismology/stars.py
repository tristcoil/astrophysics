"""The per-star pipeline: light curve -> spectrum -> numax, Dnu, epsilon, Delta Pi_1.

``analyse_star(kic)`` runs every step for one red giant and caches two files
in ``data/cache/stars/``:

* ``<kic>.json``: every measured number (and the settings that produced it);
* ``<kic>.npz``: the arrays the figures need (the log-binned spectrum, the
  native-resolution spectrum around numax divided by the background, and
  posterior samples of numax).

The steps:

1. prepare the Kepler long-cadence light curve (``lightcurves``);
2. Parseval-calibrated power density spectrum (``spectrum``);
3. blind maximum-likelihood fit of background + Gaussian envelope from 24
   starting values of numax, then an emcee posterior (``background``);
4. Dnu from the autocorrelation of the background-divided spectrum, refined
   with the universal pattern, which also gives epsilon (``deltanu``);
5. the radial modes in each order and a straight-line fit through them;
6. the gravity-mode period spacing Delta Pi_1 (``period_spacing``).
"""

from __future__ import annotations

import json
import time
import traceback

import numpy as np

import background as B
import config as C
import deltanu as D
import lightcurves as L
import period_spacing as P
import spectrum as S

#: Effective number of modes per radial order (Bedding et al. 2010; Stello et
#: al. 2011), used for the amplitude per radial mode (Kjeldsen et al. 2008).
MODES_PER_ORDER = 3.04


def _paths(kic: int):
    return C.STAR_CACHE_DIR / f"{kic}.json", C.STAR_CACHE_DIR / f"{kic}.npz"


def load_star(kic: int) -> tuple[dict, dict] | None:
    pj, pn = _paths(kic)
    if not (pj.exists() and pn.exists()):
        return None
    with np.load(pn) as z:
        arrays = {k: z[k] for k in z.files}
    return json.loads(pj.read_text()), arrays


def global_parameters(ts: L.TimeSeries, *, seed: int, nyquist_note: str = "",
                      search=C.NUMAX_SEARCH, fit_nu_min=C.FIT_NU_MIN_UHZ,
                      fit_nu_max=None, n_perturb=C.DNU_PERTURBATIONS,
                      period_spacing: bool = True, eps_width: float = 0.15,
                      oversample: int = 1, dnu_method: str = "pattern") -> tuple[dict, dict]:
    """Steps 2-6 on a prepared light curve (shared with 16 Cyg and the Sun).

    ``dnu_method``: "pattern" (red giants) takes the autocorrelation Dnu
    refined with the red-giant universal pattern; "acf" (main-sequence stars
    and the Sun, whose spectra do not follow that pattern) keeps the
    autocorrelation value, the pattern fit then only places epsilon for the
    echelle diagrams.
    """
    out: dict = {"label": ts.label, "n_points": int(len(ts.time)),
                 "span_days": ts.span_days, "duty_cycle": ts.duty_cycle,
                 "rms_ppm": float(np.std(ts.flux)), "cadence_s": ts.cadence_s}
    t0 = time.time()
    sp = S.power_density(ts.time, ts.flux, cadence_s=ts.cadence_s, oversample=oversample)
    out["parseval_ratio"] = float(np.sum(sp.psd) * sp.df / np.var(ts.flux))
    out["resolution_uhz"] = sp.resolution
    out["nyquist_uhz"] = sp.nyquist
    bs = B.bin_for_fit(sp, trend_sigma_days=ts.trend_sigma_days, nu_min=fit_nu_min,
                       nu_max=fit_nu_max)
    fit = B.fit_envelope(bs, search=search)
    fit.samples = B.sample_posterior(fit, bs, seed=seed)
    summ = fit.summary()
    out["fit"] = summ
    post = summ["posterior"]
    numax, numax_lo, numax_hi = post["numax"]
    out["numax"] = numax
    out["numax_err"] = 0.5 * (numax_lo + numax_hi)
    sig = post["sigma"][0]
    out["envelope_fwhm"] = 2.3548 * sig
    out["delta_lnL_envelope"] = summ["delta_lnL_envelope"]
    out["t_fit_s"] = time.time() - t0

    th = fit.theta.copy()
    th[B.IDX["numax"]] = numax
    th[B.IDX["sigma"]] = sig
    eta2 = S.apodization(sp.nu, sp.nyquist)
    tr = L.filter_transfer(sp.nu, ts.trend_sigma_days)
    bg = B.background_only(th, sp.nu, eta2, tr)
    snr = sp.psd / bg

    acf = D.measure_dnu(sp.nu, sp.psd, bg, numax, sig, seed=seed, n_perturb=n_perturb)
    out["dnu_acf"] = acf["dnu"]
    out["dnu_acf_err"] = acf["dnu_err"]
    out["dnu_expected"] = acf["dnu_expected"]
    ref = D.refine_dnu(sp.nu, sp.psd, bg, numax, sig, acf["dnu"], seed=seed, eps_width=eps_width)
    out["dnu_pattern"] = ref["dnu"]
    out["dnu_pattern_err"] = ref["dnu_err"]
    out["dnu_method"] = dnu_method
    if dnu_method == "pattern":
        out["dnu"], out["dnu_err"] = ref["dnu"], ref["dnu_err"]
    else:
        out["dnu"], out["dnu_err"] = acf["dnu"], acf["dnu_err"]
    out["epsilon"] = ref["epsilon"]
    rad = D.radial_modes(sp.nu, snr, ref["dnu"], ref["epsilon"], numax, sig)
    out["radial"] = {k: v for k, v in rad.items()}

    # Height and granulation in the convention of Yu et al. / the SYD
    # pipeline, whose power density is twice ours (a mode's peak has area A^2
    # rather than A^2/2); see README.
    H = post["H"][0]
    out["env_height"] = H
    out["amp_per_radial_mode"] = float(np.sqrt(2.0 * H * out["dnu"] / MODES_PER_ORDER))
    out["gran_at_numax"] = float(2.0 * B.granulation_at(th, numax))

    if period_spacing:
        t1 = time.time()
        ps = P.measure_dpi1(sp.nu, snr, out["dnu"], numax, out["epsilon"], sig)
        power = ps.pop("power", None)
        out["dpi1"] = ps
        ob = P.observed_period_spacing(sp.nu, snr, out["dnu"], numax, out["epsilon"])
        out["dp_obs"] = {"dp_obs": ob["dp_obs"], "height": ob["height"]}
        out["t_dpi1_s"] = time.time() - t1
    else:
        power = None

    # arrays for figures: native spectrum within 5 envelope widths of numax
    w = (sp.nu > max(sp.nu[0], numax - 5 * sig)) & (sp.nu < min(sp.nyquist, numax + 5 * sig))
    arrays = {
        "bin_nu": bs.nu.astype(np.float32), "bin_psd": bs.p.astype(np.float32),
        "bin_n": bs.n.astype(np.float32),
        "nu": sp.nu[w].astype(np.float64), "psd": sp.psd[w].astype(np.float32),
        "bg": bg[w].astype(np.float32),
        "numax_samples": fit.samples[::8, B.IDX["numax"]].astype(np.float32),
        "theta_ml": fit.theta, "theta_post": np.array([post[n][0] for n in B.PARAM_NAMES]),
        "acf_lags": acf["lags"][: int(1.8 * acf["dnu_expected"] / sp.df)].astype(np.float32),
        "acf": acf["acf"][: int(1.8 * acf["dnu_expected"] / sp.df)].astype(np.float32),
        "dnu_trials": ref["trials"].astype(np.float32), "dnu_score": ref["score"].astype(np.float32),
    }
    if power is not None:
        arrays["dpi1_power"] = power
    return out, arrays


def analyse_star(kic: int, *, force: bool = False) -> dict:
    kic = int(kic)
    pj, pn = _paths(kic)
    if pj.exists() and pn.exists() and not force:
        return json.loads(pj.read_text())
    t0 = time.time()
    try:
        ts = L.kepler_long_cadence(kic)
        out, arrays = global_parameters(ts, seed=kic % (2 ** 31))
        out["kic"] = kic
        out["ok"] = True
    except Exception as e:  # recorded, so one bad star does not stop the run
        out = {"kic": kic, "ok": False, "error": f"{type(e).__name__}: {e}",
               "traceback": traceback.format_exc()}
        arrays = {}
    out["t_total_s"] = time.time() - t0
    np.savez_compressed(pn, **arrays)
    pj.write_text(json.dumps(out, indent=1, default=float))
    return out


def _worker(kic):
    from threadpoolctl import threadpool_limits

    with threadpool_limits(1):
        r = analyse_star(kic)
    return kic, r.get("ok", False), r.get("t_total_s", 0.0)


def analyse_all(kics, *, workers: int = 6, quiet: bool = False) -> None:
    """Run the pipeline on every star not yet cached (resumable)."""
    import multiprocessing as mp

    todo = [int(k) for k in kics if load_star(int(k)) is None]
    if not todo:
        return
    if not quiet:
        print(f"  analysing {len(todo)} stars on {workers} processes "
              f"(~{len(todo) * 40 / workers / 60:.0f} min)...")
    ctx = mp.get_context("fork")
    t0 = time.time()
    with ctx.Pool(workers) as pool:
        for i, (kic, ok, dt) in enumerate(pool.imap_unordered(_worker, todo), 1):
            if not quiet:
                flag = "" if ok else "  FAILED"
                print(f"    {i}/{len(todo)} KIC {kic} {dt:.0f} s{flag} ({time.time() - t0:.0f} s)",
                      flush=True)
