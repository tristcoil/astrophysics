"""Stacking the CMB on supervoids and superclusters.

The compensated top-hat filter
------------------------------
For a structure at position p and a radius R, the filtered temperature is the
mean temperature inside R minus the mean in the surrounding ring from R to
sqrt(2) R (same area as the disc). The ring subtraction removes fluctuations
larger than the structure -- the primary CMB has most of its power on large
scales, so this matters enormously (Granett et al. 2008).

Linear algebra instead of loops
-------------------------------
A filter is a fixed weighted sum of pixels: +1/N_in on the disc, -1/N_ring on
the ring, 0 elsewhere. Writing every filter we will ever need (every structure,
every radius) as a row of one sparse matrix S turns "stack this sky" into the
single product S @ map. That is what makes 1,000 simulated skies cheap. Big
filters are evaluated on a coarser HEALPix grid (at least ~8 pixels across the
radius) so the matrix stays small.

Blinding
--------
``real_values`` refuses to touch the real CMB for a catalogue whose analysis
choices are not frozen in data/preregistration.json (see config.PREREGISTRATION).
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone

import numpy as np
from scipy import sparse

import cmb
from config import (MIN_UNMASKED, NSIDE, PREREG_PATH, PREREGISTRATION, prereg_hash)

#: Coarsest resolution used for big filters, and pixels required across R.
NSIDES = (512, 256, 128)
PIX_PER_RADIUS = 8.0


def _hp():
    import healpy as hp
    return hp


def pick_nside(radius_deg: float) -> int:
    hp = _hp()
    for ns in NSIDES[::-1]:  # coarsest first
        if radius_deg / np.degrees(hp.nside2resol(ns)) >= PIX_PER_RADIUS:
            return ns
    return NSIDES[0]


def multires(m, mask) -> tuple[dict, dict]:
    """The map and mask at every resolution in NSIDES (masked-mean degrading)."""
    maps, masks = {NSIDES[0]: m}, {NSIDES[0]: mask}
    for ns in NSIDES[1:]:
        maps[ns], masks[ns] = cmb.degrade(m, mask, ns, 0.5)
    return maps, masks


class FastDegrader:
    """``multires`` for many maps sharing one mask: same result, ~20x faster.

    In HEALPix's NESTED ordering the 4^k children of a coarse pixel are
    consecutive, so degrading is a reshape and a sum. The reorderings and the
    mask coverage fractions are computed once.
    """

    def __init__(self, mask):
        hp = _hp()
        n0 = NSIDES[0]
        self.mask = mask
        self.to_nest = hp.nest2ring(n0, np.arange(hp.nside2npix(n0)))  # ring index of nest pixel
        mnest = mask[self.to_nest].astype(np.float64)
        self.levels = []
        for ns in NSIDES[1:]:
            f = (n0 // ns) ** 2
            frac = mnest.reshape(-1, f).mean(axis=1)
            keep = frac >= 0.5 - 1e-9
            back = hp.ring2nest(ns, np.arange(hp.nside2npix(ns)))  # nest index of ring pixel
            self.levels.append((ns, f, frac, keep, back))
        self.masks = {n0: mask, **{ns: keep[back] for ns, _, _, keep, back in self.levels}}

    def __call__(self, m) -> dict:
        mnest = np.where(self.mask, m, 0.0)[self.to_nest].astype(np.float64)
        maps = {NSIDES[0]: m}
        for ns, f, frac, keep, back in self.levels:
            s = mnest.reshape(-1, f).mean(axis=1)
            out = np.zeros_like(s)
            out[keep] = s[keep] / frac[keep]
            maps[ns] = out[back].astype(np.float32)
        return maps


@dataclass
class FilterBank:
    """Many compensated top-hat filters, one sparse matrix per resolution."""

    lon: np.ndarray  # Galactic l (deg) of each filter's centre
    lat: np.ndarray
    radius: np.ndarray  # deg
    nside: np.ndarray
    valid: np.ndarray
    mats: dict = field(default_factory=dict)  # nside -> (row indices, csr matrix)

    def __len__(self):
        return len(self.radius)

    @classmethod
    def build(cls, lon, lat, radius, masks: dict) -> "FilterBank":
        hp = _hp()
        lon, lat, radius = (np.asarray(x, float) for x in (lon, lat, radius))
        n = len(radius)
        nside = np.array([pick_nside(r) for r in radius])
        valid = np.zeros(n, bool)
        vec = cmb.gal_to_vec(lon, lat)
        bank = cls(lon, lat, radius, nside, valid)
        for ns in np.unique(nside):
            rows = np.flatnonzero(nside == ns)
            mask = masks[ns]
            indptr, indices, data = [0], [], []
            for r in rows:
                inner = hp.query_disc(ns, vec[r], np.radians(radius[r]))
                outer = hp.query_disc(ns, vec[r], np.radians(radius[r] * np.sqrt(2)))
                ring = np.setdiff1d(outer, inner, assume_unique=True)
                ki, kr = inner[mask[inner]], ring[mask[ring]]
                ok = (len(ki) >= MIN_UNMASKED * len(inner) and len(kr) >= MIN_UNMASKED * len(ring)
                      and len(ki) > 0 and len(kr) > 0)
                if ok:
                    valid[r] = True
                    indices.extend((ki, kr))
                    data.extend((np.full(len(ki), 1.0 / len(ki), np.float32),
                                 np.full(len(kr), -1.0 / len(kr), np.float32)))
                    indptr.append(indptr[-1] + len(ki) + len(kr))
                else:
                    indptr.append(indptr[-1])
            ind = np.concatenate(indices) if indices else np.zeros(0, int)
            dat = np.concatenate(data) if data else np.zeros(0, np.float32)
            mat = sparse.csr_matrix((dat, ind, np.array(indptr)),
                                    shape=(len(rows), hp.nside2npix(ns)))
            bank.mats[int(ns)] = (rows, mat)
        return bank

    def apply(self, maps: dict) -> np.ndarray:
        """Filtered temperature of every filter (NaN where invalid)."""
        out = np.full(len(self), np.nan, np.float64)
        for ns, (rows, mat) in self.mats.items():
            out[rows] = mat @ maps[ns].astype(np.float64)
        out[~self.valid] = np.nan
        return out

    @property
    def nnz(self) -> int:
        return sum(m.nnz for _, m in self.mats.values())


# --------------------------------------------------------------------------
# Pre-registration and blinding
# --------------------------------------------------------------------------


def frozen_entries() -> dict:
    if PREREG_PATH.exists():
        return json.loads(PREREG_PATH.read_text())
    return {}


def freeze(keys: list[str]) -> dict:
    """Write the named pre-registration entries (with hash and UTC time) to disk.

    An entry already frozen is never overwritten: if config.py has changed since,
    that is reported and the old record kept, so the record shows what was
    decided before the data were seen.
    """
    frozen = frozen_entries()
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    for k in keys:
        entry = PREREGISTRATION[k]
        h = prereg_hash(entry)
        if k in frozen:
            if frozen[k]["sha256"] != h:
                raise RuntimeError(f"pre-registration '{k}' was frozen on {frozen[k]['frozen_utc']} "
                                   f"and config.py now differs from it; refusing to change it")
            continue
        frozen[k] = {"frozen_utc": now, "sha256": h, "entry": entry}
    PREREG_PATH.write_text(json.dumps(frozen, indent=1))
    return frozen


def assert_unblinded(*keys: str) -> None:
    frozen = frozen_entries()
    for k in keys:
        if k not in frozen:
            raise PermissionError(f"analysis '{k}' is still blinded: freeze its choices first")
        if frozen[k]["sha256"] != prereg_hash(PREREGISTRATION[k]):
            raise PermissionError(f"analysis '{k}': config.py differs from the frozen choices")


# --------------------------------------------------------------------------
# Running the bank on data and simulations
# --------------------------------------------------------------------------


def data_maps(which: str = "smica"):
    d = cmb.prepare_maps()
    return multires(d[which], d["mask"])


def run_sims(bank: FilterBank, n_sims: int, *, extra=None, log_every: int = 50) -> np.ndarray:
    """Filtered values for simulated skies 0..n_sims-1: array (n_sims, n_filters).

    ``extra`` is an optional callable(i, alm, maps512) returning a 1-D array of
    additional statistics per sky; its results are returned as a second array.
    """
    d = cmb.prepare_maps()
    mask = d["mask"]
    cl = cmb.theory_cl()
    deg = FastDegrader(mask)
    out = np.empty((n_sims, len(bank)), np.float32)
    extras = []
    t0 = time.monotonic()
    for i in range(n_sims):
        alm = cmb.sim_alm(i, cl)
        m = cmb.alm_to_map(alm, NSIDE)
        out[i] = bank.apply(deg(m))
        if extra is not None:
            extras.append(extra(i, alm, m))
        if log_every and (i + 1) % log_every == 0:
            el = time.monotonic() - t0
            print(f"    sim {i + 1}/{n_sims}  ({el:.0f} s, {el / (i + 1):.2f} s/sky)", flush=True)
    if extra is not None:
        return out, np.asarray(extras)
    return out


# --------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------


def stack_mean(values: np.ndarray, rows: np.ndarray) -> np.ndarray:
    """Mean over the given filter rows, ignoring invalid (NaN) filters.

    Works on a 1-D data vector or a 2-D (sims, filters) array.
    """
    v = values[..., rows]
    return np.nanmean(v, axis=-1)


def summarize(obs: float, sims: np.ndarray, sign: int = 1) -> dict:
    """Observed value against its simulated null distribution.

    ``sign`` is the direction LCDM predicts (+1 hot, -1 cold). The p-value is
    one-sided in that direction, counted from the simulations; the Gaussian
    significance is obs / sigma.
    """
    sims = np.asarray(sims, float)
    mu, sd = float(np.mean(sims)), float(np.std(sims, ddof=1))
    p = (1 + np.sum(sign * sims >= sign * obs)) / (len(sims) + 1)
    return {"value": float(obs), "sigma": sd, "null_mean": mu, "snr": float(obs / sd),
            "p_one_sided": float(p), "n_sims": int(len(sims))}


# --------------------------------------------------------------------------
# Stacked images
# --------------------------------------------------------------------------


def gnomonic_patch(m, mask, lon, lat, size_deg, npix, nside):
    """A square gnomonic (tangent-plane) patch of a HEALPix map; masked pixels NaN."""
    hp = _hp()
    half = np.tan(np.radians(size_deg / 2))
    x = np.linspace(-half, half, npix)
    xx, yy = np.meshgrid(x, x)  # xx: east (increasing l), yy: north
    # tangent-plane point -> unit vector, centre at (0, 0, 1) first
    v = np.stack([xx, yy, np.ones_like(xx)], axis=-1)
    v /= np.linalg.norm(v, axis=-1, keepdims=True)
    lo, la = np.radians(lon), np.radians(lat)
    # basis at the centre: e_r (radial), e_l (east), e_b (north)
    er = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    el = np.array([-np.sin(lo), np.cos(lo), 0.0])
    eb = np.cross(er, el)
    world = v[..., 0, None] * el + v[..., 1, None] * eb + v[..., 2, None] * er
    pix = hp.vec2pix(nside, world[..., 0], world[..., 1], world[..., 2])
    out = m[pix].astype(np.float64)
    out[~mask[pix]] = np.nan
    return out


def stack_image(m, mask, lon, lat, signs, *, size_deg, npix, scale=None, nside=NSIDE):
    """Mean of patches centred on each structure (voids sign-flipped if signs = -1).

    With ``scale`` (per-structure angular radius), each patch spans
    size_deg * scale, i.e. the image axis is in units of the structure's radius.
    """
    acc = np.zeros((npix, npix))
    cnt = np.zeros((npix, npix))
    for i in range(len(lon)):
        sz = size_deg * (scale[i] if scale is not None else 1.0)
        p = gnomonic_patch(m, mask, lon[i], lat[i], sz, npix, nside) * signs[i]
        ok = np.isfinite(p)
        acc[ok] += p[ok]
        cnt[ok] += 1
    with np.errstate(invalid="ignore"):
        return acc / cnt
