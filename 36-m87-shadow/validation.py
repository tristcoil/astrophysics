"""Imaging validation on SIMULATED data (clearly labelled as such in every figure).

Question: could the ring be an artefact of sparse sampling plus our imaging
choices? Test: give the imager data from sources that are *not* rings, sampled
exactly like the real April 11 observation (same baselines, times and thermal
noise, plus random station gains as large as the real calibration allows), and
see what comes out. The EHT ran the same kind of test (M87 Paper IV, section 5).
"""

from __future__ import annotations

import json
from dataclasses import replace

import numpy as np

import models as M
import uvdata as U
from config import CACHE_DIR, M87_2017_BANDS, SYSTEMATIC_FRACTION
from imaging import ClosureImager, gain_sigma

VAL_DIR = CACHE_DIR / "validation"
COMPACT_FLUX = 0.6  # Jy, as assumed for M87's compact emission


def _truths():
    """Name -> (model, theta) for the simulated sources (compact flux scaled to 0.6 Jy)."""
    disk = M.SymmetricShape("disk", large=False)
    gauss = M.SymmetricShape("gaussian", large=False)
    two = M.GaussianSum(large=False)
    ring = M.MGRing(1, large=False, floor=False)
    return {
        "disk": ("a filled disk, 45 µas", disk, np.array([45.0])),
        "gaussian": ("a Gaussian blob, 40 µas FWHM", gauss, np.array([40.0])),
        "double": ("two blobs", two, np.array([28.0, 0.8, 0.5, 0.6, 26.0, -18.0, 16.0, 0.7, 1.2])),
        "ring": ("a lopsided ring, 44 µas", ring, np.array([44.0, 10.0, 0.25, np.radians(-170.0)])),
    }


def simulate(raw: U.VisData, model, theta, rng) -> U.VisData:
    """Synthetic 10-s visibilities with the real (u, v), noise and scan structure.

    Each station gets a random amplitude gain per scan (log-normal with the a-priori
    calibration error of that station) and a random phase per scan (the real data
    were phase-stabilised within scans by the HOPS pipeline). Intra-site baselines
    see an extra 0.6 Jy of arcsecond-scale emission, as M87's do.
    """
    V = model.vis(theta, raw.u, raw.v)
    V0 = model.vis(theta, np.array([1e3]), np.array([0.0]))[0]
    V = COMPACT_FLUX * V / V0
    V = V + np.where(raw.intra_site(), COMPACT_FLUX, 0.0)
    gains = {}
    for s in np.unique(raw.scan):
        for st in raw.stations:
            gains[(s, st)] = np.exp(rng.normal(0, gain_sigma(st)) + 1j * rng.uniform(-np.pi, np.pi))
    g = np.array([gains[(s, a)] * np.conj(gains[(s, b)]) for s, a, b in zip(raw.scan, raw.t1, raw.t2)])
    noise = raw.sigma * (rng.normal(size=len(V)) + 1j * rng.normal(size=len(V)))
    return replace(raw, vis=V * g + noise)


def prepared_from(raw: U.VisData) -> U.Prepared:
    avg = U.with_systematics(U.scan_average(raw), SYSTEMATIC_FRACTION)
    return U.Prepared(raw, avg, U.closure_phases(avg), U.log_closure_amplitudes(avg, exclude_intra_site=True))


def run(day: str = "101", *, settings: dict, weights: dict, seed: int = 5, refit: bool = False) -> dict:
    """Image each simulated source; cache truth and reconstruction."""
    import analysis as A

    VAL_DIR.mkdir(parents=True, exist_ok=True)
    out = {}
    rng = np.random.default_rng(seed)
    raws = [U.read_uvfits(A.m87_2017_path(day, b)) for b in M87_2017_BANDS]
    for key, (label, model, theta) in _truths().items():
        path = VAL_DIR / f"{key}_{day}.npz"
        if path.exists() and not refit:
            z = np.load(path)
            out[key] = {"label": label, "truth": z["truth"], "recon": z["recon"],
                        "reduced_chi2": json.loads(str(z["chi2"]))}
            continue
        preps = [prepared_from(simulate(r, model, theta, rng)) for r in raws]
        im = ClosureImager(preps, amplitudes=True, **settings)
        recon = im.image(weights=dict(weights))
        truth = model.image(theta, im.X, im.Y)
        truth = COMPACT_FLUX * truth / truth.sum()
        np.savez_compressed(path, truth=truth, recon=recon, chi2=json.dumps(im.history[-1]))
        out[key] = {"label": label, "truth": truth, "recon": recon, "reduced_chi2": im.history[-1]}
    return out
