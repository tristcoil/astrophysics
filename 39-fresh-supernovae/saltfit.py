"""SALT3 light-curve fits with sncosmo.

SALT models a Type Ia supernova's spectrum at every phase as

    F(phase, wavelength) = x0 * [M0 + x1 * M1] * exp(c * CL(wavelength))

M0 is the average supernova, M1 how it changes with "stretch" x1 (broader,
slower light curves have larger x1), and CL a colour law (positive c = redder).
x0 is the overall amplitude, from which the peak magnitude mB follows. The
model spectrum is redshifted, dimmed by Milky Way dust and integrated through
the ZTF g, r and i passbands, and (t0, x0, x1, c) are adjusted until it matches
the photometry. Redshift and Milky Way dust are fixed, not fitted.

The procedure follows the ZTF DR2 papers (Rigault et al. 2025, Sect. 4.7): a
first fit within -15..+50 rest-frame days of a peak-time guess pins down t0,
then a second fit uses only -10..+40 days, where SALT is well trained. The
model's own uncertainty is included in the fit (sncosmo's ``modelcov=True``);
without it the colour is biased and its error is underestimated (checked on
ZTF18aaqfziz: c = 0.164 +- 0.011 without, 0.184 +- 0.025 with, against the
ZTF team's 0.190 +- 0.028).
"""

from __future__ import annotations

import hashlib
import json
import os
import time
import warnings

import numpy as np
import pandas as pd
from scipy import stats

import config

# Point sncosmo at the project's own cache before anything else touches it,
# so models and passbands always come from data/sncosmo (and runs work offline).
import sncosmo as _sncosmo  # noqa: E402

_sncosmo.conf.data_dir = str(config.SNCOSMO_DIR)

#: First-pass and final rest-frame phase windows [days relative to peak].
WINDOW_FIRST = (-15.0, 50.0)
WINDOW_FINAL = (-10.0, 40.0)

#: Fit bounds. Wide on purpose: the quality cuts come afterwards, so a
#: peculiar supernova shows up as an extreme value instead of a failed fit.
BOUNDS_X1 = (-5.0, 5.0)
BOUNDS_C = (-1.0, 2.0)
T0_SEARCH = 15.0  # days either side of the guess, first pass

#: SNANA/Pantheon+ convention linking amplitude and peak magnitude:
#: mB = -2.5 log10(x0) + 10.635 holds exactly in the Pantheon+ table.
MB_OFFSET = 10.635

#: The two models we fit with.
#:   "salt3"     - SALT3 (Kenworthy et al. 2021), the model the spec asks for;
#:                 every ZTF-internal result uses it.
#:   "salt2-b21" - D1's model: the SALT2 surfaces the Pantheon+ team retrained
#:                 (Brout et al. 2022), from the Pantheon+ release, used the way
#:                 SNANA uses them, so D1's alpha, beta and M apply exactly.
MODELS = ("salt3", "salt2-b21")

_MODELS: dict = {}


def b21_info() -> dict:
    """Settings in the B21 model's SALT2.INFO that SNANA applies and sncosmo
    does not: the rest-frame wavelength range of usable bands, and a global
    magnitude offset added to every model magnitude (MAG_OFFSET)."""
    info = {}
    for line in (config.SALT2_B21_DIR / "SALT2.INFO").read_text().splitlines():
        key, _, val = line.partition(":")
        if key.strip() in ("RESTLAMBDA_RANGE", "MAG_OFFSET"):
            nums = [float(v) for v in val.split("#")[0].split()]
            info[key.strip()] = nums if len(nums) > 1 else nums[0]
    return info


def fit_fingerprint(kind: str = "salt3") -> str:
    """Hash of everything that changes the fits, to invalidate stale caches."""
    import sncosmo

    spec = {
        "source": config.SALT_SOURCE, "version": config.SALT_VERSION,
        "first": WINDOW_FIRST, "final": WINDOW_FINAL, "x1": BOUNDS_X1, "c": BOUNDS_C,
        "t0": T0_SEARCH, "rv": config.MW_RV, "sfd": config.SFD_SCALING,
        "sncosmo": sncosmo.__version__, "dust": "F99", "modelcov": True,
    }
    if kind != "salt3":
        spec.update(source=kind, version="Pantheon+ MODEL000", restframe=b21_info())
    return hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:16]


def get_model(kind: str = "salt3"):
    """The SALT model with Milky Way dust, built once per process."""
    if kind not in _MODELS:
        import sncosmo

        sncosmo.conf.data_dir = str(config.SNCOSMO_DIR)
        if kind == "salt3":
            source = sncosmo.get_source(config.SALT_SOURCE, version=config.SALT_VERSION)
        elif kind == "salt2-b21":
            source = sncosmo.SALT2Source(modeldir=str(config.SALT2_B21_DIR),
                                         name="salt2-b21", version="Pantheon+")
        else:
            raise ValueError(kind)
        _MODELS[kind] = sncosmo.Model(source=source,
                                      effects=[sncosmo.F99Dust(r_v=config.MW_RV)],
                                      effect_names=["mw"], effect_frames=["obs"])
    return _MODELS[kind]


def _restframe_bands(lc: pd.DataFrame, z: float, kind: str) -> pd.DataFrame:
    """SNANA drops bands whose rest-frame mean wavelength is outside the
    model's RESTLAMBDA_RANGE (2800-7000 A for B21): ZTF i below z ~ 0.13."""
    if kind == "salt3":
        return lc
    import sncosmo

    lo, hi = b21_info()["RESTLAMBDA_RANGE"]
    keep = {b for b in lc["band"].unique()
            if lo <= sncosmo.get_bandpass(b).wave_eff / (1 + z) <= hi}
    return lc[lc["band"].isin(keep)]


def _table(lc: pd.DataFrame):
    from astropy.table import Table

    return Table({"time": lc["mjd"].to_numpy(), "band": lc["band"].to_numpy(),
                  "flux": lc["flux"].to_numpy(), "fluxerr": lc["fluxerr"].to_numpy(),
                  "zp": lc["zp"].to_numpy(), "zpsys": np.full(len(lc), "ab")})


def _select(lc: pd.DataFrame, t0: float, z: float, window) -> pd.DataFrame:
    phase = (lc["mjd"] - t0) / (1 + z)
    return lc[lc["good"] & (phase >= window[0]) & (phase <= window[1])]


def _guess_t0(lc: pd.DataFrame, z: float) -> float:
    """Time of the highest significant flux, smoothed over 3 nights."""
    det = lc[lc["good"] & (lc["flux"] / lc["fluxerr"] > 5)]
    if det.empty:
        return np.nan
    # Normalise each band by its own peak so a bright band does not dominate.
    norm = det["flux"] / det.groupby("band")["flux"].transform("max")
    t = det["mjd"].to_numpy()
    score = [norm[np.abs(t - ti) < 1.5].mean() for ti in t]
    return float(t[int(np.argmax(score))])


def fit_lightcurve(lc: pd.DataFrame, z: float, mwebv: float, t0_guess: float | None = None,
                   name: str = "", kind: str = "salt3") -> dict:
    """Two-pass SALT fit. Returns a flat dict (one row of the results table).

    Never raises for bad data: failures come back with ``ok = False`` and a
    ``fail`` reason, so one odd supernova cannot stop a run of thousands.
    """
    import sncosmo

    out = {"name": name, "z_fit": z, "mwebv_used": mwebv, "ok": False, "fail": ""}
    t_start = time.time()
    if t0_guess is None or not np.isfinite(t0_guess):
        t0_guess = _guess_t0(lc, z)
    if not np.isfinite(t0_guess):
        out["fail"] = "no detection to guess the peak time"
        return out
    out["t0_guess"] = t0_guess

    model = get_model(kind)
    model.set(z=z, mwebv=mwebv, t0=t0_guess, x0=1e-3, x1=0.0, c=0.0)
    lc_all = lc
    lc = _restframe_bands(lc, z, kind)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # Pass 1: locate the peak.
            data = _select(lc, t0_guess, z, WINDOW_FIRST)
            if len(data) < 5:
                out["fail"] = f"only {len(data)} good points near the guessed peak"
                return out
            res, fitted = sncosmo.fit_lc(
                _table(data), model, ["t0", "x0", "x1", "c"],
                bounds={"t0": (t0_guess - T0_SEARCH, t0_guess + T0_SEARCH),
                        "x1": BOUNDS_X1, "c": BOUNDS_C},
                modelcov=False, guess_t0=False)
            t0_1 = fitted.get("t0")
            # Pass 2 (and a third if the peak moved by more than 2 days):
            # the final window, with the model covariance.
            for _ in range(2):
                data = _select(lc, t0_1, z, WINDOW_FINAL)
                if len(data) < 5:
                    out["fail"] = f"only {len(data)} good points in -10..+40 d"
                    return out
                fitted.set(t0=t0_1)
                res, fitted = sncosmo.fit_lc(
                    _table(data), fitted, ["t0", "x0", "x1", "c"],
                    bounds={"t0": (t0_1 - 10, t0_1 + 10), "x1": BOUNDS_X1, "c": BOUNDS_C},
                    modelcov=True, guess_amplitude=False, guess_t0=False)
                moved = abs(fitted.get("t0") - t0_1)
                t0_1 = fitted.get("t0")
                if moved < 2.0:
                    break
    except Exception as err:  # sncosmo raises a zoo of error types
        out["fail"] = f"{type(err).__name__}: {err}"[:200]
        return out

    names = list(res.vparam_names)
    cov = res.covariance
    if cov is None or not np.all(np.isfinite(cov)):
        out["fail"] = "no covariance (fit did not converge)"
        return out
    par = dict(zip(res.param_names, res.parameters))
    for i, k in enumerate(names):
        out[k] = float(par[k])
        out[k + "_err"] = float(np.sqrt(cov[i, i]))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if j > i:
                out[f"cov_{a}_{b}"] = float(cov[i, j])
    out["chisq"] = float(res.chisq)
    out["ndof"] = int(res.ndof)
    out["fitprob"] = float(stats.chi2.sf(res.chisq, res.ndof)) if res.ndof > 0 else np.nan
    out["n_points"] = int(len(data))
    out["bands_fitted"] = "".join(sorted({b[-1] for b in data["band"]}))
    at_bound = (abs(out["x1"] - BOUNDS_X1[0]) < 1e-3 or abs(out["x1"] - BOUNDS_X1[1]) < 1e-3
                or abs(out["c"] - BOUNDS_C[0]) < 1e-3 or abs(out["c"] - BOUNDS_C[1]) < 1e-3)
    out["at_bound"] = bool(at_bound)

    # Peak magnitude and its covariances (linearised: dmB/dx0 = -2.5/(ln10 x0)).
    if out["x0"] > 0:
        k = -2.5 / np.log(10) / out["x0"]
        out["mB"] = -2.5 * np.log10(out["x0"]) + MB_OFFSET
        out["mB_err"] = abs(k) * out["x0_err"]
        out["cov_mB_x1"] = k * out["cov_x0_x1"]
        out["cov_mB_c"] = k * out["cov_x0_c"]
        out["ok"] = True
    else:
        out["fail"] = "non-positive amplitude"

    from lightcurves import detected_phases, good_sampling

    covg = detected_phases(lc_all, out["t0"], z, WINDOW_FINAL)
    out.update({("det_bands" if k == "bands" else f"lc_{k}"): v for k, v in covg.items()})
    out["good_sampling"] = good_sampling(covg)
    out["fit_seconds"] = round(time.time() - t_start, 3)
    return out


# ---------------------------------------------------------------------------
# Many supernovae in parallel
# ---------------------------------------------------------------------------

def _init_worker(kind: str = "salt3"):
    # One BLAS thread per worker: several forked workers each spawning 16
    # OpenBLAS threads only fight over the cores.
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    try:
        from threadpoolctl import threadpool_limits

        threadpool_limits(1)
    except ImportError:
        pass
    get_model(kind)


def _fit_dr2_task(task: dict) -> dict:
    from lightcurves import read_dr2_lightcurve

    name = task["name"]
    try:
        lc, header = read_dr2_lightcurve(name)
    except FileNotFoundError:
        return {"name": name, "ok": False, "fail": "no light curve file in the release"}
    guess = header.get("mjd_max_init", np.nan)
    return fit_lightcurve(lc, task["z"], task["mwebv"], guess, name=name,
                          kind=task.get("kind", "salt3"))


def fit_dr2_sample(tasks: list[dict], n_workers: int = 6, cache_path=None,
                   progress_every: int = 200, kind: str = "salt3") -> pd.DataFrame:
    """Fit every DR2 supernova in `tasks` (dicts with name, z, mwebv).

    Results are cached in a parquet file tagged with the fit fingerprint; a
    cache made with other settings is ignored and rebuilt. While running,
    finished fits are checkpointed every `progress_every` supernovae to a
    ``.partial.parquet`` file, and a restarted run resumes from it.
    """
    import multiprocessing as mp

    fingerprint = fit_fingerprint(kind)
    tasks = [{**t, "kind": kind} for t in tasks]
    names = {t["name"] for t in tasks}
    if cache_path is not None and cache_path.exists():
        cached = pd.read_parquet(cache_path)
        if "fingerprint" in cached and (cached["fingerprint"] == fingerprint).all() \
                and names <= set(cached["name"]):
            return cached
    partial_path = None if cache_path is None else cache_path.with_suffix(".partial.parquet")
    results = []
    if partial_path is not None and partial_path.exists():
        part = pd.read_parquet(partial_path)
        if "fingerprint" in part and (part["fingerprint"] == fingerprint).all():
            results = [r for r in part.drop(columns="fingerprint").to_dict("records")
                       if r["name"] in names]
            print(f"    resuming: {len(results)} fits already done", flush=True)
    done = {r["name"] for r in results}
    todo = [t for t in tasks if t["name"] not in done]

    def checkpoint(path):
        df = _parquet_safe(pd.DataFrame(results))
        df["fingerprint"] = fingerprint
        tmp = path.with_suffix(".tmp")
        df.to_parquet(tmp)
        tmp.replace(path)
        return df

    t_start = time.time()
    if todo:
        ctx = mp.get_context("fork")
        with ctx.Pool(n_workers, initializer=_init_worker, initargs=(kind,)) as pool:
            for i, r in enumerate(pool.imap_unordered(_fit_dr2_task, todo, chunksize=4), 1):
                results.append(r)
                if i % progress_every == 0 or i == len(todo):
                    rate = i / (time.time() - t_start)
                    print(f"    {len(results)}/{len(tasks)} fitted ({rate:.1f}/s, "
                          f"{sum(not x['ok'] for x in results)} failed)", flush=True)
                    if partial_path is not None:
                        checkpoint(partial_path)
    if cache_path is None:
        df = _parquet_safe(pd.DataFrame(results))
        df["fingerprint"] = fingerprint
        return df
    df = checkpoint(cache_path)
    if partial_path.exists():
        partial_path.unlink()
    return df


def _parquet_safe(df: pd.DataFrame) -> pd.DataFrame:
    """Failed fits lack most keys, leaving object columns that mix values with
    NaN. Give those columns one type: nullable booleans, or strings."""
    for col in df.columns[df.dtypes == object]:
        vals = df[col].dropna()
        if len(vals) and all(isinstance(v, (bool, np.bool_)) for v in vals):
            df[col] = df[col].astype("boolean")
        else:
            df[col] = df[col].where(df[col].notna(), None).astype("string")
    return df
