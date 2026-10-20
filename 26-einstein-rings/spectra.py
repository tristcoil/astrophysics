"""Two galaxies in one fibre: how SLACS found its lenses.

Each SLACS lens was first noticed in its SDSS spectrum. A 3"-wide fibre placed
on a massive elliptical galaxy also caught light from a galaxy far behind it,
lensed into an arc or ring that falls partly inside the fibre. The spectrum then
shows the elliptical's starlight at one redshift and, on top, the emission lines
of the distant star-forming galaxy -- [O II], H-beta, [O III], sometimes
H-alpha -- at a second, higher redshift (Bolton et al. 2004).

We repeat that search blind. The SDSS pipeline's best-fit galaxy model is
subtracted from each spectrum, slow template mismatch is filtered out, and
every trial background redshift is scored by how strongly the residuals show
the emission-line pattern there. The catalogued source redshift is consulted
only afterwards, to grade the result.

Wavelengths are vacuum, as in SDSS spectra.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter

from config import SDSS_DIR

#: Emission lines of a star-forming galaxy (vacuum wavelengths, Angstrom), as
#: the search uses them: each entry is one template of (wavelength, relative
#: strength). The [O II] pair is unresolved-to-marginal at SDSS resolution and
#: its ratio depends on gas density, so it enters 1:1. The [O III] pair shares
#: its upper level, so its ratio is fixed by atomic physics: 4959/5007 = 1/2.98.
EMISSION_GROUPS = {
    "[O II]": ((3727.092, 1.0), (3729.875, 1.0)),
    "H-beta": ((4862.683, 1.0),),
    "[O III]": ((4960.295, 1 / 2.98), (5008.240, 1.0)),
    "H-alpha": ((6564.614, 1.0),),
}

#: The same lines one by one, for labelling and for per-line signal-to-noise.
EMISSION_LINES = {
    "[O II]": (3727.092, 3729.875),
    "H-beta": (4862.683,),
    "[O III] 4959": (4960.295,),
    "[O III] 5007": (5008.240,),
    "H-alpha": (6564.614,),
}

#: Lines the *lens* galaxy itself may emit: the low-ionisation nuclear emission
#: common in massive ellipticals ([O II], H-beta, [O III], [O I], [N II],
#: H-alpha, [S II]). Their residuals are masked at the lens redshift, so the
#: search for a second galaxy is not fooled by the first. Weaker lines are left
#: alone: masking them would hide some sources' [O II] for no gain.
FOREGROUND_LINES = (3727.092, 3729.875, 4862.68, 4960.30, 5008.24, 6302.05, 6365.54,
                    6549.86, 6564.61, 6585.27, 6718.29, 6732.67)
FOREGROUND_HALF_WIDTH_KMS = 350.0

#: Absorption features of the lens galaxy's old stars, for labelling.
ABSORPTION_LINES = {
    "Ca K": 3934.78, "Ca H": 3969.59, "G band": 4305.6, "H-beta": 4862.68,
    "Mg b": 5176.7, "Na D": 5895.6,
}

#: Bright night-sky lines whose subtraction residuals mimic emission lines.
SKY_LINES = (5578.5, 5891.6, 5897.6, 6301.7, 6365.5, 7246.0)
SKY_HALF_WIDTH = 6.0  # Angstrom

DLOG = 1e-4  # SDSS pixels are uniform in log10(wavelength)
LINE_SIGMA_KMS = 110.0  # line width: source-galaxy motions plus instrument (~70 km/s)
C_KMS = 299792.458


@dataclass
class Spectrum:
    name: str
    wave: np.ndarray  # vacuum Angstrom
    flux: np.ndarray  # 1e-17 erg/s/cm^2/A
    ivar: np.ndarray
    model: np.ndarray  # SDSS pipeline's best-fit galaxy model at the lens redshift
    z_pipeline: float
    vdisp: float
    vdisp_err: float
    sn_median: float

    @property
    def residual(self) -> np.ndarray:
        return self.flux - self.model


@lru_cache(maxsize=None)
def read_spectrum(name: str) -> Spectrum:
    with fits.open(SDSS_DIR / f"{name}.fits") as hdul:
        co = hdul["COADD"].data
        so = hdul["SPECOBJ"].data
        return Spectrum(
            name=name,
            wave=10 ** co["loglam"].astype(float),
            flux=co["flux"].astype(float),
            ivar=co["ivar"].astype(float),
            model=co["model"].astype(float),
            z_pipeline=float(so["Z"][0]),
            vdisp=float(so["VDISP"][0]),
            vdisp_err=float(so["VDISP_ERR"][0]),
            sn_median=float(so["SN_MEDIAN_ALL"][0]),
        )


def clean_residual(spec: Spectrum, window: int = 101, z_lens: float | None = None
                   ) -> tuple[np.ndarray, np.ndarray]:
    """Residual with slow template mismatch removed, and honest weights.

    A running median (101 pixels, ~7000 km/s) takes out broad continuum
    mismatch while leaving lines a few pixels wide untouched. Pixels near
    bright sky lines, or without data, get zero weight. The weights are then
    rescaled so the residuals have unit robust scatter: the pipeline's errors
    do not include template mismatch."""
    w = spec.ivar.copy()
    for sky in SKY_LINES:
        w[np.abs(spec.wave - sky) < SKY_HALF_WIDTH] = 0.0
    if z_lens is not None:
        for lam in FOREGROUND_LINES:
            c = lam * (1 + z_lens)
            w[np.abs(spec.wave - c) < c * FOREGROUND_HALF_WIDTH_KMS / C_KMS + 2.0] = 0.0
    r = spec.residual
    r_hp = r - median_filter(np.where(w > 0, r, 0.0), size=window, mode="nearest")
    good = w > 0
    chi = r_hp[good] * np.sqrt(w[good])
    scale = 1.4826 * np.median(np.abs(chi - np.median(chi)))
    if np.isfinite(scale) and scale > 0:
        w = w / scale**2
    return r_hp, w


def _matched_filter(r, w, centers_pix, offsets_pix, sigma_pix, weights=None):
    """Amplitude S/N of a Gaussian line (or a fixed-ratio group of lines) at
    each trial centre: sum(r t w) / sqrt(sum(t^2 w))."""
    weights = weights or [1.0] * len(offsets_pix)
    half = int(np.ceil(5 * sigma_pix + np.max(np.abs(offsets_pix))))
    idx = np.arange(-half, half + 1)
    base = np.floor(centers_pix).astype(int)
    frac = centers_pix - base
    pix = base[:, None] + idx[None, :]
    inside = (pix >= 0) & (pix < len(r))
    pixc = np.clip(pix, 0, len(r) - 1)
    t = np.zeros(pix.shape)
    for off, wt in zip(offsets_pix, weights):
        t += wt * np.exp(-0.5 * ((idx[None, :] - frac[:, None] - off) / sigma_pix) ** 2)
    t *= inside
    num = np.sum(r[pixc] * w[pixc] * t, axis=1)
    den = np.sum(t * t * w[pixc], axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        sn = np.where(den > 0, num / np.sqrt(den), 0.0)
    # A line needs most of its profile observed to count.
    coverage = np.sum((w[pixc] > 0) * t, axis=1) / np.maximum(np.sum(t, axis=1), 1e-9)
    return np.where(coverage > 0.6, sn, 0.0)


@dataclass
class SearchResult:
    name: str
    z_grid: np.ndarray
    score: np.ndarray  # sum of squared positive line S/N
    line_sn: dict
    z_best: float
    score_best: float
    lines_at_best: dict


def search_background(spec: Spectrum, z_lens: float, z_max: float = 1.45,
                      min_dz: float = 0.05) -> SearchResult:
    """Score every trial background redshift by its emission-line evidence:
    the sum of the squared (positive) S/N of each line group."""
    r, w = clean_residual(spec, z_lens=z_lens)
    z_lo = z_lens + min_dz
    step = 2e-5  # in log10(1+z), ~14 km/s
    lz = np.arange(np.log10(1 + z_lo), np.log10(1 + z_max), step)
    z_grid = 10**lz - 1
    group_sn = _scan(spec, r, w, lz, EMISSION_GROUPS)
    score = sum(np.clip(sn, 0, None) ** 2 for sn in group_sn.values())
    i = int(np.argmax(score))
    per_line = _scan(spec, r, w, lz[i:i + 1],
                     {k: tuple((lam, 1.0) for lam in v) for k, v in EMISSION_LINES.items()})
    return SearchResult(spec.name, z_grid, score, group_sn, float(z_grid[i]), float(score[i]),
                        {k: float(v[0]) for k, v in per_line.items()})


def _scan(spec, r, w, lz, groups) -> dict:
    """S/N of each line group at every trial log10(1+z) in ``lz``."""
    loglam0 = np.log10(spec.wave[0])
    sigma_pix = LINE_SIGMA_KMS / C_KMS / np.log(10) / DLOG
    out = {}
    for name, lines in groups.items():
        ref = lines[-1][0]  # the strongest line anchors the template
        pix = (np.log10(ref) + lz - loglam0) / DLOG
        offsets = [(np.log10(lam) - np.log10(ref)) / DLOG for lam, _ in lines]
        weights = [wt for _, wt in lines]
        ok = (pix + min(offsets) > 5) & (pix + max(offsets) < len(r) - 6)
        sn = np.zeros(len(lz))
        if ok.any():
            sn[ok] = _matched_filter(r, w, pix[ok], offsets, sigma_pix, weights)
        out[name] = sn
    return out


def line_sn_at(spec: Spectrum, z: float, z_lens: float | None = None) -> dict:
    """S/N of each emission line, one by one, at a given redshift."""
    r, w = clean_residual(spec, z_lens=z_lens)
    per = _scan(spec, r, w, np.array([np.log10(1 + z)]),
                {k: tuple((lam, 1.0) for lam in v) for k, v in EMISSION_LINES.items()})
    return {k: float(v[0]) for k, v in per.items()}
