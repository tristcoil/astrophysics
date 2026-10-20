"""Every fit and number in the project; slow steps are cached in data/cache/.

    python analysis.py fit [group ...]   # run the missing nested-sampling fits
    python analysis.py images            # our RML images + the EHT pipeline images
    python analysis.py summarize         # derived numbers -> data/results.json

Groups: m87_2017 (fiducial ring model, 8 data sets), xsring (Paper VI's crescent,
8 data sets), comparison (13 models on April 11), m87_2018, sgra. ``run.py``
calls all of these in order.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import subprocess
import sys
import time
import zlib
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")

import numpy as np  # noqa: E402

import fitting as F  # noqa: E402
import models as M  # noqa: E402
import uvdata as U  # noqa: E402
from config import (C, CACHE_DIR, G, H1_FALLBACK, H1_RESULTS_PATH, M87_2017_BANDS, M87_2017_DAYS,  # noqa: E402
                    M87_2018_BANDS, M_SUN, MPC, PAPER6_XSRING, PAPER6_XSRINGAUSS, PC, PUBLISHED, RAD_PER_UAS,
                    RAW_DIR, RESULTS_PATH, SGRA_DAYS, SHADOW_FACTOR_SCHWARZSCHILD, ensure_dirs)

FIT_DIR = CACHE_DIR / "fits"
IMG_DIR = CACHE_DIR / "images"
EHTIM_DIR = CACHE_DIR / "ehtim"
NJOBS = int(os.environ.get("I5_NJOBS", "6"))

# ---------------------------------------------------------------------------
# Observations
# ---------------------------------------------------------------------------
M87_2017_UV = RAW_DIR / "m87_2017" / "EHTC_FirstM87Results_Apr2019" / "uvfits"
M87_2018_UV = RAW_DIR / "m87_2018" / "uvfits"
SGRA_UV = RAW_DIR / "sgra_2017" / "EHTC_FirstSgrAResults_May2022" / "uvfits_norm"
M87_2018_USED_DAYS = {"111": "April 21", "115": "April 25"}  # April 22 has 4-5 stations: too few closures


def m87_2017_path(day: str, band: str) -> Path:
    return M87_2017_UV / f"SR1_M87_2017_{day}_{band}_hops_netcal_StokesI.uvfits"


def m87_2018_path(day: str, band: str, pipeline: str = "hops") -> Path:
    return M87_2018_UV / f"L2V1_M87_2018_{day}_{band}_{pipeline}_netcal_10s_StokesI.uvfits"


def sgra_path(day: str, band: str, pipeline: str = "hops") -> Path:
    return SGRA_UV / f"ER6_SGRA_2017_{day}_{band}_{pipeline}_netcal-LMTcal-norm_StokesI.uvfits"


OBS: dict[str, dict] = {}
for _day, _date in M87_2017_DAYS.items():
    for _band in M87_2017_BANDS:
        OBS[f"m87_2017_{_day}_{_band}"] = dict(kind="m87", year=2017, day=_day, band=_band,
                                                label=f"M87 2017 {_date} {_band.upper()}",
                                                paths=[m87_2017_path(_day, _band)])
for _day, _date in M87_2018_USED_DAYS.items():
    for _band in ("b3", "b4"):  # 227 and 229 GHz, the 2017 frequencies (b1/b2 are 213/215 GHz)
        OBS[f"m87_2018_{_day}_{_band}"] = dict(kind="m87", year=2018, day=_day, band=_band,
                                                label=f"M87 2018 {_date} {_band}",
                                                paths=[m87_2018_path(_day, _band)])
for _band in ("lo", "hi"):
    OBS[f"sgra_2017_{_band}"] = dict(kind="sgra", year=2017, band=_band, label=f"Sgr A* 2017 April 6+7 {_band.upper()}",
                                     paths=[sgra_path(day, _band) for day in SGRA_DAYS])
    OBS[f"sgra_2017_{_band}_casa"] = dict(kind="sgra", year=2017, band=_band,
                                          label=f"Sgr A* 2017 April 6+7 {_band.upper()} (CASA pipeline)",
                                          paths=[sgra_path(day, _band, "casa") for day in SGRA_DAYS])


def load_obs(key: str):
    """M87: scan-averaged, 1% systematic, closures. Sgr A*: the same with thermal errors only
    and no deblurring (the fitted noise model and the scattered model do the rest)."""
    o = OBS[key]
    if o["kind"] == "sgra":
        return U.prepare_sgra_raw(o["paths"], o["label"])
    return U.prepare(o["paths"][0], o["label"])


def load_obs_deblurred(key: str):
    """Sgr A* for plots: deblurred, with Paper IV's measured noise budget added (before deblurring)."""
    o = OBS[key]
    return U.prepare_sgra(o["paths"], o["label"], a4=PUBLISHED["sgra_var_a4"], b=PUBLISHED["sgra_var_b"])


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------
MODELS = {
    "gaussian": lambda: M.SymmetricShape("gaussian"),
    "disk": lambda: M.SymmetricShape("disk"),
    "ring": lambda: M.SymmetricShape("ring"),
    "two_gauss": lambda: M.GaussianSum(),
    "crescent": lambda: M.Crescent(),
    "mg1": lambda: M.MGRing(1),
    "mg2": lambda: M.MGRing(2),
    "mg3": lambda: M.MGRing(3),
    "mg1g2": lambda: M.MGRing(1, n_gauss=2),
    "mg2g1": lambda: M.MGRing(2, n_gauss=1),
    "mg2g2": lambda: M.MGRing(2, n_gauss=2),
    "mg2g3": lambda: M.MGRing(2, n_gauss=3),
    "mg3g2": lambda: M.MGRing(3, n_gauss=2),
    "xs2": lambda: M.XSRing(n_gauss=2),
    "xs3": lambda: M.XSRing(n_gauss=3),
    # Sgr A*: Paper IV's mG-ring (ring + central Gaussian floor, its priors), orders 1-3,
    # with the noise budget fitted alongside (fractional, refractive floor, variability).
    # No large-scale component: Sgr A* has no jet, and the light-curve-normalised data put
    # exactly the total flux on the shortest baselines (a spare component there let the fit
    # decouple the long baselines from the total flux and pick a 40-uas crescent).
    **{f"sg{m}": (lambda m=m: F.NoisyModel(M.MGRing(m, large=False, w_max=40.0, floor_flux_max=4.0,
                                                    floor_fwhm_range=(20.0, 200.0)))) for m in (1, 2, 3)},
    # Sgr A* from its amplitudes beyond 2 Glambda (the nulls), flux scale free, no floor needed there
    # (the fractional term absorbs residual station-gain errors; the data want ~20%, so it may reach 30%)
    **{f"sa{m}": (lambda m=m: F.AmplitudeModel(M.MGRing(m, large=False, floor=False, w_max=40.0, total_flux=True),
                                               f_max=0.3))
       for m in (1, 2)},
    # M87 from its amplitudes beyond 2 Glambda: no scattering; the LMT's 2017 gains were too unstable
    "ma1": lambda: F.AmplitudeModel(M.MGRing(1, large=False, floor=False, w_max=40.0, total_flux=True),
                                    scattered=False, exclude_stations=("LM",)),
}
FIDUCIAL = "mg2g2"  # chosen by Bayesian evidence on April 11 (see the 'comparison' group)
COMPARISON_OBS = "m87_2017_101_lo"
COMPARISON_MODELS = ["gaussian", "disk", "ring", "two_gauss", "crescent", "mg2", "mg3", "mg1g2", "mg2g1",
                     "mg2g2", "mg3g2", "xs2"]

GROUPS = {
    "m87_2017": [(k, FIDUCIAL) for k in OBS if k.startswith("m87_2017")],
    "xsring": [(k, "xs2") for k in OBS if k.startswith("m87_2017")],
    # (the crescent fit to April 11 comes from the 'xsring' group; the summary reads it from the cache)
    "comparison": [(COMPARISON_OBS, m) for m in COMPARISON_MODELS if m not in ("xs2", FIDUCIAL)],
    "m87_2018": [(k, FIDUCIAL) for k in OBS if k.startswith("m87_2018")],
    "sgra": [("sgra_2017_lo", "sg2"), ("sgra_2017_lo", "sg3"), ("sgra_2017_hi", "sg2"), ("sgra_2017_hi", "sg3"),
             ("sgra_2017_lo", "sg1"), ("sgra_2017_hi", "sg1"), ("sgra_2017_lo_casa", "sg2"),
             ("sgra_2017_hi_casa", "sg2")],
    "sgra_amp": [("sgra_2017_lo", "sa1"), ("sgra_2017_hi", "sa1"), ("sgra_2017_lo", "sa2"), ("sgra_2017_hi", "sa2"),
                 ("sgra_2017_lo_casa", "sa1"), ("sgra_2017_hi_casa", "sa1")],
    "m87_amp": [(k, "ma1") for k in OBS if k.startswith("m87_2017")],
}


def fit_path(obs_key: str, model_key: str) -> Path:
    return FIT_DIR / f"{obs_key}__{model_key}"


def run_fit(obs_key: str, model_key: str, *, refit: bool = False, nlive: int | None = None,
            seed_offset: int = 0) -> dict:
    """Nested sampling with random-slice proposals.

    Settings chosen by experiment (README, 'Sampling'): on April 11 the default
    random-walk sampler got stuck in local modes (best chi^2 276-280, posteriors
    0.05 uas wide), while slice sampling with 1000 live points found chi^2 222 and
    reproduced ln Z to 0.3 and the diameter to 0.1 uas between independent seeds.
    """
    path = fit_path(obs_key, model_key)
    cached = None if refit else F.load_fit(path)
    if cached is not None:
        return cached
    FIT_DIR.mkdir(parents=True, exist_ok=True)
    model = MODELS[model_key]()
    prep = load_obs(obs_key)
    t0 = time.time()
    seed = zlib.crc32(f"{obs_key}/{model_key}".encode()) % 10_000 + seed_offset  # reproducible across runs
    # nuisance Gaussians make the posterior rugged and multimodal: give those fits more live points
    nlive = nlive or (1000 if getattr(model, "n_gauss", 0) else 500)
    fit = F.run_nested(prep, model, nlive=nlive, dlogz=0.1, njobs=NJOBS, seed=seed, sample="rslice")
    F.save_fit(path, fit, model)
    s = fit.summary()
    print(f"  {obs_key:22s} {model_key:9s} lnZ {fit.logz:8.1f}  chi2/N "
          + ", ".join(f"{k} {v:.2f}" for k, v in s["reduced_chi2"].items())
          + f"  ({time.time() - t0:.0f} s)", flush=True)
    return F.load_fit(path)


def run_group(group: str, *, refit: bool = False) -> None:
    for obs_key, model_key in GROUPS[group]:
        run_fit(obs_key, model_key, refit=refit)


# ---------------------------------------------------------------------------
# Images
# ---------------------------------------------------------------------------
IMAGER_SETTINGS = dict(npix=64, fov_uas=128.0, flux=0.6, prior_fwhm_uas=40.0)
IMAGER_WEIGHTS = dict(cphase=1.0, lcamp=1.0, amp=1.0, mem=1.0, tsv=1.0, tv=0.3, centroid=0.1)


def our_image(day: str, *, weights=None, refit: bool = False) -> np.ndarray:
    """Our own RML image of one 2017 night (both bands), cached."""
    from imaging import ClosureImager

    tag = "" if weights is None else "_" + "_".join(f"{k}{v:g}" for k, v in sorted(weights.items()))
    path = IMG_DIR / f"m87_2017_{day}{tag}.npy"
    if path.exists() and not refit:
        return np.load(path)
    IMG_DIR.mkdir(parents=True, exist_ok=True)
    preps = [U.prepare(m87_2017_path(day, b), exclude_intra_site=True) for b in M87_2017_BANDS]
    im = ClosureImager(preps, amplitudes=True, **IMAGER_SETTINGS)
    w = dict(IMAGER_WEIGHTS)
    w.update(weights or {})
    img = im.image(weights=w)
    np.save(path, img)
    (path.with_suffix(".json")).write_text(json.dumps({"reduced_chi2": im.history[-1], "weights": w}))
    return img


def ehtim_pipeline_image(day: str) -> np.ndarray | None:
    """Run the EHT's own eht-imaging script (2019-D01-02) on one night, both bands."""
    from astropy.io import fits

    out = EHTIM_DIR / f"m87_2017_{day}.fits"
    if not out.exists():
        EHTIM_DIR.mkdir(parents=True, exist_ok=True)
        src = (RAW_DIR / "eht-imaging_pipeline.py").read_text()
        old = "ttype     = 'nfft'              # Type of Fourier transform ('direct', 'nfft', or 'fast')"
        if old not in src:
            raise RuntimeError("the EHT pipeline script changed; cannot apply the one-line NFFT patch")
        patched = EHTIM_DIR / "eht-imaging_pipeline_direct.py"
        patched.write_text(src.replace(old, "ttype     = 'direct'            # PATCHED (I5): NFFT is not installed"))
        env = dict(os.environ, OPENBLAS_NUM_THREADS="4", OMP_NUM_THREADS="4")
        cmd = [sys.executable, str(patched), "-i", str(m87_2017_path(day, "lo")), "-i2", str(m87_2017_path(day, "hi")),
               "-o", str(out)]
        with open(EHTIM_DIR / f"pipeline_{day}.log", "w") as log:
            subprocess.run(cmd, cwd=EHTIM_DIR, stdout=log, stderr=subprocess.STDOUT, env=env, check=True)
    with fits.open(out) as h:
        return h[0].data.astype(float)[::-1, :]  # FITS rows run south -> north; ours north -> south


def ring_features(img: np.ndarray, fov_uas: float = 128.0, *, r_range=(5.0, 35.0), search: float = 8.0) -> dict:
    """Measure a ring in an image the way the EHT's image-domain feature extraction does.

    For trial centres, take radial profiles at 360 angles and find each one's peak;
    the centre is where the peak radii scatter least (REx, Chael 2019; M87 Paper IV
    sec. 9). Diameter = twice the mean peak radius; width = mean FWHM of the radial
    profiles; brightness PA = direction of the intensity-weighted mean of the ring.
    """
    from scipy.ndimage import map_coordinates

    n = img.shape[0]
    ps = fov_uas / n
    # pixel coordinates: column 0 = east (+x), row 0 = north (+y)
    def to_pix(x, y):
        return (0.5 * fov_uas - 0.5 * ps - y) / ps, (0.5 * fov_uas - 0.5 * ps - x) / ps

    from scipy.ndimage import gaussian_filter

    X = (0.5 * fov_uas - 0.5 * ps) - np.arange(n) * ps
    # start the centre search at the ring's dark hole: the faintest point of the smoothed image
    # inside the bright region (the light-weighted centroid is pulled towards the bright side)
    sm = gaussian_filter(img, 5.0 / ps)
    tot = img.sum()
    cx0 = np.sum(img.sum(axis=0) * X) / tot
    cy0 = np.sum(img.sum(axis=1) * X) / tot
    XX, YY = np.meshgrid(X, X)
    near = np.hypot(XX - cx0, YY - cy0) < 15.0
    k = np.argmin(np.where(near, sm, np.inf))
    cx, cy = XX.ravel()[k], YY.ravel()[k]
    th = np.radians(np.arange(360))
    rr = np.linspace(0, r_range[1] + 10, 200)

    def profiles(x0, y0):
        xs = x0 + np.outer(np.sin(th), rr)
        ys = y0 + np.outer(np.cos(th), rr)
        r_, c_ = to_pix(xs, ys)
        return map_coordinates(img, [r_.ravel(), c_.ravel()], order=1).reshape(xs.shape)

    def score(x0, y0):
        P = profiles(x0, y0)
        sel = (rr >= r_range[0]) & (rr <= r_range[1])
        rpk = rr[sel][np.argmax(P[:, sel], axis=1)]
        return np.std(rpk), rpk, P

    best = None
    for x0 in np.arange(cx - search, cx + search + 0.1, 1.0):
        for y0 in np.arange(cy - search, cy + search + 0.1, 1.0):
            s, rpk, P = score(x0, y0)
            if best is None or s < best[0]:
                best = (s, x0, y0, rpk, P)
    _, x0, y0, rpk, P = best
    widths = []
    for k in range(len(th)):
        prof = P[k]
        ipk = np.argmin(abs(rr - rpk[k]))
        half = prof[ipk] / 2
        lo = np.flatnonzero(prof[:ipk] < half)
        hi = np.flatnonzero(prof[ipk:] < half)
        if len(lo) and len(hi):
            widths.append(rr[ipk + hi[0]] - rr[lo[-1]])
    peak_I = P[np.arange(len(th)), np.argmin(abs(rr[None, :] - rpk[:, None]), axis=1)]
    pa = np.degrees(np.angle(np.sum(peak_I * np.exp(1j * th)))) % 360
    # ratio of the mean brightness at the centre to that on the ring
    centre = map_coordinates(img, [[to_pix(x0, y0)[0]], [to_pix(x0, y0)[1]]], order=1)[0]
    return {"d": float(2 * np.mean(rpk)), "d_scatter": float(2 * np.std(rpk)), "width": float(np.median(widths)),
            "pa_deg": float(pa), "x0": float(x0), "y0": float(y0), "contrast": float(np.mean(peak_I) / max(centre, 1e-12))}


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------
def _fit_summary(obs_key: str, model_key: str) -> dict | None:
    f = F.load_fit(fit_path(obs_key, model_key))
    if f is None:
        return None
    s = f["summary"]
    out = {"logz": s["logz"], "logz_err": s["logz_err"], "reduced_chi2": s["reduced_chi2"], "chi2": s["chi2"],
           "n_data": s["n_data"], "n_params": len(f["names"])}
    if s.get("posterior_method"):
        out["posterior_method"] = s["posterior_method"]
    # derived quantities recomputed from the samples with the current model code
    model = MODELS[model_key]()
    if hasattr(model, "derived"):
        f["derived"] = F.derived_samples(model, f["samples"])
    for k, v in f["derived"].items():
        out[k] = F.summarize_angle_deg(v) if k == "pa_deg" else F.summarize(v)
    der = f["derived"]
    if "d" in der and "w" in der:
        # peak-brightness diameter of the blurred ring (Sgr A* Paper IV eq. 49), comparable with image rings
        out["d_peak"] = F.summarize(der["d"] - der["w"] ** 2 / (4 * np.log(2) * der["d"]))
    return out


def _combine(values: list[float], errors: list[float]) -> dict:
    v, e = np.asarray(values), np.asarray(errors)
    w = 1 / e**2
    mean = float(np.sum(w * v) / np.sum(w))
    return {"weighted_mean": mean, "stat_err": float(np.sqrt(1 / np.sum(w))), "mean": float(np.mean(v)),
            "rms_scatter": float(np.std(v, ddof=1)) if len(v) > 1 else 0.0, "n": int(len(v)),
            "error_of_mean": float(np.std(v, ddof=1) / np.sqrt(len(v))) if len(v) > 1 else float(e[0])}


def theta_g_to_mass(theta_uas, distance_mpc):
    return np.asarray(theta_uas) * RAD_PER_UAS * np.asarray(distance_mpc) * MPC * C**2 / G / M_SUN


def mass_to_theta_g(mass_msun, distance_mpc):
    return G * np.asarray(mass_msun) * M_SUN / C**2 / (np.asarray(distance_mpc) * MPC) / RAD_PER_UAS


def h1_prediction() -> dict:
    """H1's (17-galactic-centre) Sgr A* mass and distance, from its results file if present."""
    if H1_RESULTS_PATH.exists():
        r = json.loads(H1_RESULTS_PATH.read_text())["fits"]["joint_gr"]
        der, par = r["derived"], r["parameters"]
        return {"source": "17-galactic-centre/data/results.json (joint Keck+VLT GR fit)",
                "mass_1e6": (par["gm6"]["median"], par["gm6"]["sd"]), "r0_kpc": (par["r0"]["median"], par["r0"]["sd"]),
                "theta_g_uas": (der["theta_g_uas"]["median"], der["theta_g_uas"]["sd"]),
                "shadow_diameter_uas": (der["shadow_diameter_uas"]["median"], der["shadow_diameter_uas"]["sd"])}
    out = dict(H1_FALLBACK)
    out["source"] = "H1 README (results file not found on this machine)"
    return out


def summarize() -> dict:
    res: dict = {"published": PUBLISHED}
    rng = np.random.default_rng(42)

    # --- data volumes ---
    data = {}
    for key, o in OBS.items():
        if o["kind"] != "m87" or o["year"] != 2017:
            continue
        prep = U.prepare(o["paths"][0], o["label"])
        data[key] = {"n_10s": len(prep.raw), "n_scan_avg": len(prep.avg), "n_cphase": len(prep.cphase),
                     "n_lcamp": len(prep.lcamp), "stations": prep.raw.stations,
                     "max_baseline_glambda": float(prep.avg.uvdist.max() / 1e9),
                     "n_scans": int(len(np.unique(prep.avg.scan)))}
        triv = U.trivial_closure_phases(prep.avg)
        z = np.angle(np.exp(1j * triv[:, 0])) / triv[:, 1]
        data[key].update({"trivial_cphase_n": int(len(triv)), "trivial_cphase_chi2_per_n": float(np.mean(z**2)),
                          "trivial_cphase_median_abs_deg": float(np.degrees(np.median(np.abs(triv[:, 0]))))})
    res["data_2017"] = data

    # --- M87 2017 fits ---
    per = {}
    for key, o in OBS.items():
        if not key.startswith("m87_2017"):
            continue
        entry = {"label": o["label"], "fiducial": _fit_summary(key, FIDUCIAL), "xsring": _fit_summary(key, "xs2"),
                 "amplitudes": _fit_summary(key, "ma1")}
        entry["paper6_xsring"] = PAPER6_XSRING[(o["day"], o["band"])]
        entry["paper6_xsringauss"] = PAPER6_XSRINGAUSS[(o["day"], o["band"])]
        per[key] = entry
    res["m87_2017_fits"] = per
    fid = [v["fiducial"] for v in per.values() if v["fiducial"]]
    if fid:
        res["m87_2017_ring"] = {
            "d": _combine([f["d"]["median"] for f in fid], [f["d"]["sd"] for f in fid]),
            "w": _combine([f["w"]["median"] for f in fid], [max(f["w"]["sd"], 0.1) for f in fid]),
            "asymmetry": _combine([f["asymmetry"]["median"] for f in fid], [f["asymmetry"]["sd"] for f in fid]),
            "pa_by_day": {d: [per[f"m87_2017_{d}_{b}"]["fiducial"]["pa_deg"]["median"] for b in M87_2017_BANDS
                              if per[f"m87_2017_{d}_{b}"]["fiducial"]] for d in M87_2017_DAYS},
            "reduced_chi2_cphase": [f["reduced_chi2"]["cphase"] for f in fid],
            "reduced_chi2_lcamp": [f["reduced_chi2"]["lcamp"] for f in fid],
        }
    am = [v["amplitudes"] for v in per.values() if v["amplitudes"]]
    if am:
        res["m87_2017_amplitude_ring"] = {"d": _combine([f["d"]["median"] for f in am], [f["d"]["sd"] for f in am])}
    xs = [v["xsring"] for v in per.values() if v["xsring"]]
    if xs:
        ours = [f["d"]["median"] for f in xs]
        theirs = [v["paper6_xsring"]["d"][0] for v in per.values() if v["xsring"]]
        res["m87_2017_xsring"] = {"d": _combine(ours, [f["d"]["sd"] for f in xs]),
                                   "paper6_mean": float(np.mean(theirs)),
                                   "differences_uas": [float(a - b) for a, b in zip(ours, theirs)]}

    # --- model comparison ---
    comp = {}
    for m in COMPARISON_MODELS:
        s = _fit_summary(COMPARISON_OBS, m)
        if s:
            comp[m] = {"name": MODELS[m]().name, "logz": s["logz"], "logz_err": s["logz_err"],
                       "reduced_chi2": s["reduced_chi2"], "n_params": s["n_params"],
                       "d": s.get("d", {}).get("median") if isinstance(s.get("d"), dict) else None}
    if comp:
        best = max(v["logz"] for v in comp.values())
        for v in comp.values():
            v["delta_logz"] = v["logz"] - best
    res["model_comparison"] = comp

    # --- mass ---
    if xs and fid:
        res["m87_mass"] = mass_budget(res, rng)

    # --- 2018 ---
    per18 = {k: {"label": o["label"], "fiducial": _fit_summary(k, FIDUCIAL)} for k, o in OBS.items()
             if k.startswith("m87_2018")}
    res["m87_2018_fits"] = per18
    f18 = [v["fiducial"] for v in per18.values() if v["fiducial"]]
    if f18:
        res["m87_2018_ring"] = {"d": _combine([f["d"]["median"] for f in f18], [f["d"]["sd"] for f in f18]),
                                "pa": [f["pa_deg"]["median"] for f in f18]}

    # --- Sgr A* ---
    sg = {}
    for b in ("lo", "hi"):
        for m in (1, 2, 3):
            s = _fit_summary(f"sgra_2017_{b}", f"sg{m}")
            if s:
                sg[f"{b}_m{m}"] = s
        s = _fit_summary(f"sgra_2017_{b}_casa", "sg2")
        if s:
            sg[f"{b}_casa_m2"] = s
    res["sgra_fits"] = sg
    sa = {}
    for b in ("lo", "hi"):
        for m in (1, 2):
            s = _fit_summary(f"sgra_2017_{b}", f"sa{m}")
            if s:
                sa[f"{b}_m{m}"] = s
        s = _fit_summary(f"sgra_2017_{b}_casa", "sa1")
        if s:
            sa[f"{b}_casa_m1"] = s
    res["sgra_amp_fits"] = sa
    res["h1"] = h1_prediction()
    if sg or sa:
        res["sgra_test"] = sgra_test(sg, res["h1"], sa)

    # --- images ---
    feats = {}
    for day in M87_2017_DAYS:
        p_ours = IMG_DIR / f"m87_2017_{day}.npy"
        p_eht = EHTIM_DIR / f"m87_2017_{day}.fits"
        e = {}
        if p_ours.exists():
            e["ours"] = ring_features(np.load(p_ours))
            e["ours"]["reduced_chi2"] = json.loads(p_ours.with_suffix(".json").read_text())["reduced_chi2"]
        if p_eht.exists():
            e["ehtim"] = ring_features(ehtim_pipeline_image(day))
        feats[day] = e
    res["image_features"] = feats
    RESULTS_PATH.write_text(json.dumps(_jsonable(res), indent=1))
    return res


def mass_budget(res: dict, rng) -> dict:
    """Ring diameter -> angular gravitational radius -> mass, under two interpretations."""
    n = 200_000
    D = rng.normal(*PUBLISHED["m87_distance_mpc"], n)
    out = {}
    # (1) Paper VI's GRMHD calibration, applied to the same crescent model it was derived for
    xs = res["m87_2017_xsring"]["d"]
    d_err = max(xs["rms_scatter"], xs["stat_err"])  # day-to-day scatter dominates
    d_xs = rng.normal(xs["weighted_mean"], d_err, n)
    alpha = rng.normal(PUBLISHED["alpha_xsring"], PUBLISHED["alpha_xsring"] * PUBLISHED["alpha_theory_frac"], n)
    theta = d_xs / alpha
    out["grmhd_calibrated"] = {
        "diameter_uas": float(xs["weighted_mean"]), "diameter_err": float(d_err), "alpha": PUBLISHED["alpha_xsring"],
        "theta_g": F.summarize(theta), "mass_1e9": F.summarize(theta_g_to_mass(theta, D) / 1e9),
        "mass_1e9_no_alpha_err": F.summarize(theta_g_to_mass(d_xs / PUBLISHED["alpha_xsring"], D) / 1e9),
    }
    # (2) the ring as the shadow's edge (pure general relativity, no emission physics)
    fd = res["m87_2017_ring"]["d"]
    d_err2 = max(fd["rms_scatter"], fd["stat_err"])
    d_ring = rng.normal(fd["weighted_mean"], d_err2, n)
    for label, a in (("shadow_schwarzschild", SHADOW_FACTOR_SCHWARZSCHILD), ("shadow_kerr_min", 9.6)):
        th = d_ring / a
        out[label] = {"diameter_uas": float(fd["weighted_mean"]), "diameter_err": float(d_err2), "alpha": a,
                      "theta_g": F.summarize(th), "mass_1e9": F.summarize(theta_g_to_mass(th, D) / 1e9)}
    # what each earlier mass predicts for the ring (distance-free: theta_dyn from Paper VI Table 9)
    pred = {}
    for name, key in (("stars (Gebhardt+11)", "theta_dyn_stars_uas"), ("gas (Walsh+13)", "theta_dyn_gas_uas")):
        th, up, lo = PUBLISHED[key]
        pred[name] = {"theta_g": th, "theta_up": up, "theta_lo": lo,
                      "crescent_diameter_uas": th * PUBLISHED["alpha_xsring"],
                      "shadow_diameter_uas": th * SHADOW_FACTOR_SCHWARZSCHILD}
    for name, key in (("stars (Liepold+23)", "mass_stars_liepold2023_1e9"), ("stars (Simon+24)", "mass_stars_simon2024_1e9")):
        m, up, lo = PUBLISHED[key]
        th = float(mass_to_theta_g(m * 1e9, 16.8))
        pred[name] = {"theta_g": th, "theta_up": th * up / m, "theta_lo": th * lo / m,
                      "crescent_diameter_uas": th * PUBLISHED["alpha_xsring"],
                      "shadow_diameter_uas": th * SHADOW_FACTOR_SCHWARZSCHILD}
    out["predictions"] = pred
    # which dynamical measurement does the ring support? (in theta_g, so distance drops out)
    tg = out["grmhd_calibrated"]["theta_g"]
    sig = 0.5 * (tg["hi68"] - tg["lo68"])
    out["tension_sigma"] = {}
    for name, p in pred.items():
        side = p["theta_up"] if tg["median"] > p["theta_g"] else p["theta_lo"]
        out["tension_sigma"][name] = float((tg["median"] - p["theta_g"]) / np.hypot(sig, side))
    return out


def sgra_test(sg: dict, h1: dict, sa: dict | None = None) -> dict:
    """Compare the Sgr A* ring with H1's prediction from the orbit of S2.

    Primary measurement: the amplitude fits (the nulls; `sa`). The closure-only fits
    (`sg`) are reported alongside: for Sgr A* they are degenerate with a broad halo.
    """
    th_h1, th_err = h1["theta_g_uas"]
    alpha_eht = PUBLISHED["sgra_ring_diameter_uas"][0] / PUBLISHED["sgra_theta_g_uas"][0]
    out = {"theta_g_h1": th_h1, "theta_g_h1_err": th_err,
           "shadow_prediction_uas": th_h1 * SHADOW_FACTOR_SCHWARZSCHILD,
           "shadow_prediction_err": th_err * SHADOW_FACTOR_SCHWARZSCHILD,
           "alpha_eht_sgra": alpha_eht, "ring_prediction_uas": alpha_eht * th_h1}

    def family(fits: dict, prefix: str) -> dict:
        res = {"best_by_band": {}}
        for b in ("lo", "hi"):
            cands = {k: v for k, v in fits.items() if k.startswith(f"{b}_m")}
            if cands:
                k = max(cands, key=lambda kk: cands[kk]["logz"])
                c = cands[k]
                res["best_by_band"][b] = {"model": f"{prefix}{k}", "d": c["d"], "d_peak": c["d_peak"], "w": c.get("w"),
                                          "pa": c.get("pa_deg"), "reduced_chi2": c["reduced_chi2"],
                                          "noise": {n: c[n]["median"] for n in ("f", "sigma_ref", "a4") if n in c}}
        allp = [v["d_peak"]["median"] for v in fits.values() if "d_peak" in v]
        res["d_peak_all_fits"] = {"min": float(min(allp)), "max": float(max(allp)), "n": len(allp)} if allp else None
        best = res["best_by_band"]
        if best:
            comb = _combine([v["d_peak"]["median"] for v in best.values()], [v["d_peak"]["sd"] for v in best.values()])
            err = max(comb["rms_scatter"], comb["stat_err"])
            if res["d_peak_all_fits"]:  # half the spread over orders and pipelines, as a systematic
                err = float(np.hypot(err, 0.5 * (res["d_peak_all_fits"]["max"] - res["d_peak_all_fits"]["min"])))
            res["d_peak"] = {"value": comb["weighted_mean"], "err": err, **comb}
            r = comb["weighted_mean"] / th_h1
            res["ratio_d_over_theta_g"] = {"value": r, "err": r * np.hypot(err / comb["weighted_mean"], th_err / th_h1)}
            res["measured_over_predicted"] = {"value": comb["weighted_mean"] / (alpha_eht * th_h1),
                                              "err": err / (alpha_eht * th_h1)}
        return res

    if sa:
        out["amplitudes"] = family(sa, "amplitudes ")
    if sg:
        out["closures"] = family(sg, "closures ")
    return out


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k) if not isinstance(k, str) else k: _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist())
    return x


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("what", choices=["fit", "fitjobs", "verify", "images", "summarize"])
    parser.add_argument("groups", nargs="*", default=list(GROUPS))
    parser.add_argument("--refit", action="store_true")
    args = parser.parse_args(argv)
    ensure_dirs()
    if args.what == "fit":
        for g in args.groups:
            print(f"[{g}]", flush=True)
            run_group(g, refit=args.refit)
    elif args.what == "fitjobs":  # explicit obs:model pairs, in the order given
        for job in args.groups:
            obs_key, model_key = job.split(":")
            run_fit(obs_key, model_key, refit=args.refit)
    elif args.what == "verify":
        for e in verify_fits(models=tuple(args.groups) if args.groups != list(GROUPS) else None):
            a = e["attempts"][-1]
            print(f"  {e['obs']:22s} {e['model']:6s} {'ok ' if e['ok'] else 'BAD'} sampler {a['sampler_max_lnl']:.1f}, "
                  f"best polished {a['best_lnl']:.1f} from {a['best_from']} ({len(e['attempts'])} attempt(s))", flush=True)
    elif args.what == "images":
        for day in M87_2017_DAYS:
            with contextlib.redirect_stdout(io.StringIO()):
                our_image(day, refit=args.refit)
            ehtim_pipeline_image(day)
            print(f"  images for {M87_2017_DAYS[day]} done", flush=True)
    else:
        summarize()
        print(f"wrote {RESULTS_PATH}")
    return 0



# ---------------------------------------------------------------------------
# Mode check: did nested sampling find the best mode?
# ---------------------------------------------------------------------------
def _polish(like, theta0, model):
    """Local Powell polish of a starting point, inside the prior box."""
    from scipy.optimize import minimize

    lo = np.array([p.lo for p in model.params])
    hi = np.array([p.hi for p in model.params])
    ang = np.array([p.kind == "angle" for p in model.params])

    def neg(t):
        t = np.where(ang, (t + np.pi) % (2 * np.pi) - np.pi, t)
        if np.any((t < lo) & ~ang) or np.any((t > hi) & ~ang):
            return 1e30
        return -like(t)

    r = minimize(neg, np.asarray(theta0, float), method="Powell", options={"maxfev": 20000, "xtol": 1e-4, "ftol": 1e-7})
    t = np.where(ang, (r.x + np.pi) % (2 * np.pi) - np.pi, r.x)
    return float(-r.fun), t


def mode_check(obs_key: str, model_key: str, peers: list[str]) -> dict:
    """Polish from this fit's best point and from the best points of `peers` (same model, other data).

    Returns the sampler's best ln L, the best polished ln L and where it came from.
    """
    f = F.load_fit(fit_path(obs_key, model_key))
    model = MODELS[model_key]()
    prep = load_obs(obs_key)
    like = F.ClosureLikelihood(prep, model) if not isinstance(model, F.NoisyModel) else F.NoiseModelLikelihood(prep, model)
    sampler_best = np.array([f["summary"]["best"][n] for n in model.names])
    out = {"sampler_max_lnl": float(like(sampler_best)), "polished": []}
    starts = [(obs_key, sampler_best)]
    for pk in peers:
        g = F.load_fit(fit_path(pk, model_key))
        if g is not None and pk != obs_key:
            starts.append((pk, np.array([g["summary"]["best"][n] for n in model.names])))
    best = (out["sampler_max_lnl"], obs_key, sampler_best)
    for src, th in starts:
        lnl, t = _polish(like, th, model)
        out["polished"].append({"from": src, "lnl": lnl})
        if lnl > best[0]:
            best = (lnl, src, t)
    out.update({"best_lnl": best[0], "best_from": best[1], "best_theta": best[2].tolist(),
                "improvement": best[0] - out["sampler_max_lnl"]})
    return out


def verify_fits(*, threshold: float = 3.0, models: tuple | None = None) -> list[dict]:
    """Mode check for every cached fit; redo those that another data set's solution beats.

    Nested sampling can miss the dominant mode of these rugged posteriors (it did for
    M87 on 10 April 2017, low band). For each fit we polish its own best point and the
    best points of every other data set fitted with the same model. If a polished
    solution started elsewhere beats the sampler's best by more than `threshold` in
    ln L, the fit is redone with a new seed and 2000 live points (up to twice).
    The outcome is stored next to the fit as <fit>.modecheck.json.
    """
    from collections import defaultdict

    by_model = defaultdict(list)
    for jobs in GROUPS.values():
        for obs_key, model_key in jobs:
            if F.load_fit(fit_path(obs_key, model_key)) is not None and obs_key not in by_model[model_key]:
                by_model[model_key].append(obs_key)
    report = []
    for model_key, obs_keys in by_model.items():
        if models and model_key not in models:
            continue
        if isinstance(MODELS[model_key](), F.NoisyModel):
            continue  # Sgr A*: two data sets with different noise; checked by the m-order and variant fits
        for obs_key in obs_keys:
            path = fit_path(obs_key, model_key).with_suffix(".modecheck.json")
            if path.exists():
                report.append(json.loads(path.read_text()))
                continue
            entry = {"obs": obs_key, "model": model_key, "attempts": []}
            for attempt in range(3):
                mc = mode_check(obs_key, model_key, obs_keys)
                entry["attempts"].append({k: mc[k] for k in ("sampler_max_lnl", "best_lnl", "best_from", "improvement")})
                if mc["best_from"] == obs_key or mc["improvement"] <= threshold or attempt == 2:
                    break
                print(f"  mode check: {obs_key} {model_key}: a solution from {mc['best_from']} is better by "
                      f"{mc['improvement']:.1f} in ln L; refitting with 2000 live points", flush=True)
                run_fit(obs_key, model_key, refit=True, nlive=2000, seed_offset=101 * (attempt + 1))
            last = entry["attempts"][-1]
            entry["ok"] = last["best_from"] == obs_key or last["improvement"] <= threshold
            if not entry["ok"]:
                print(f"  mode check: {obs_key} {model_key}: nested sampling missed the better mode {len(entry['attempts'])}"
                      f" times; sampling it locally with emcee instead", flush=True)
                summ = local_posterior(obs_key, model_key, mc["best_theta"])
                entry["fallback"] = summ["posterior_method"]
                mc2 = mode_check(obs_key, model_key, obs_keys)
                entry["attempts"].append({k: mc2[k] for k in ("sampler_max_lnl", "best_lnl", "best_from", "improvement")})
                entry["ok"] = mc2["best_from"] == obs_key or mc2["improvement"] <= threshold
            path.write_text(json.dumps(entry, indent=1))
            report.append(entry)
    return report



def local_posterior(obs_key: str, model_key: str, theta0, *, nwalkers: int = 64, nsteps: int = 6000,
                    burn: int = 3000, seed: int = 1) -> dict:
    """Fallback when nested sampling keeps missing a narrow, better mode: emcee around the polished solution.

    Walkers start in a tiny ball around `theta0` and explore the mode with the same
    likelihood and priors (uniform; log-uniform where declared; angles periodic). The
    result is that mode's posterior; there is no evidence for it, so ln Z is stored as NaN.
    """
    import emcee

    model = MODELS[model_key]()
    prep = load_obs(obs_key)
    like = F.ClosureLikelihood(prep, model)
    lo = np.array([p.lo for p in model.params])
    hi = np.array([p.hi for p in model.params])
    ang = np.array([p.kind == "angle" for p in model.params])
    logk = np.array([p.kind == "log" for p in model.params])

    def logp(t):
        t = np.where(ang, (t + np.pi) % (2 * np.pi) - np.pi, t)
        if np.any(((t < lo) | (t > hi)) & ~ang):
            return -np.inf
        return like(t) - np.sum(np.log(t[logk]))

    rng = np.random.default_rng(seed)
    scale = 1e-3 * (hi - lo)
    p0 = np.array(theta0, float) + scale * rng.normal(size=(nwalkers, model.ndim))
    p0 = np.where(ang, p0, np.clip(p0, lo + 1e-9, hi - 1e-9))
    t0 = time.time()
    sampler = emcee.EnsembleSampler(nwalkers, model.ndim, logp)
    sampler.random_state = np.random.RandomState(seed).get_state()  # emcee copies numpy's global RNG otherwise
    sampler.run_mcmc(p0, nsteps, progress=False)
    chain = sampler.get_chain(discard=burn, flat=True)
    chain[:, ang] = (chain[:, ang] + np.pi) % (2 * np.pi) - np.pi
    lp = sampler.get_log_prob(discard=burn, flat=True)
    best = chain[np.argmax(lp)]
    fit = F.FitResult(prep.label, model.name, model.names, chain[rng.choice(len(chain), 20000, replace=False)],
                      float("nan"), float("nan"), float(like(best)), best, like.chi2_parts(best),
                      {"cphase": len(prep.cphase), "lcamp": len(prep.lcamp)}, nwalkers * nsteps, time.time() - t0)
    path = fit_path(obs_key, model_key)
    F.save_fit(path, fit, model)
    summ = json.loads(path.with_suffix(".json").read_text())
    summ["posterior_method"] = ("emcee around the polished best mode (nested sampling repeatedly missed it); "
                                f"acceptance {np.mean(sampler.acceptance_fraction):.2f}")
    path.with_suffix(".json").write_text(json.dumps(summ, indent=1))
    return summ


if __name__ == "__main__":
    sys.exit(main())
