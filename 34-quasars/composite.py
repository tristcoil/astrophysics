"""A composite quasar spectrum, stacked the way Vanden Berk et al. (2001) did it.

Each SDSS spectrum shows a quasar through a fixed observed window (3,600-10,400 A
for BOSS). Shifted to the quasar's rest frame, that window lands on a different
stretch of the spectrum for every redshift: low-redshift quasars show the
optical lines, high-redshift ones the far ultraviolet. Chaining a few thousand
spectra spread evenly in redshift covers ~700-9,000 A in one average spectrum.

Steps (Vanden Berk et al. 2001, section 3):
1. Correct each spectrum for Milky Way dust (Cardelli, Clayton & Mathis 1989,
   R_V = 3.1) and shift it to the rest frame with its systemic redshift.
2. Rebin onto 1 A pixels, conserving flux.
3. Scale every spectrum to the growing composite where they overlap, ignoring
   the Lyman-alpha forest (it varies wildly from quasar to quasar); iterate.
4. Median composite: the median of the scaled spectra in each pixel (keeps line
   ratios). Geometric-mean composite: spectra covering 3020-3100 A, normalised
   there, averaged in log (keeps the mean power-law slope).
"""

from __future__ import annotations

import concurrent.futures as cf

import numpy as np
import pandas as pd
from astropy.io import fits

from config import BOSS_TRUSTED, COMPOSITE_GRID

#: BOSS pixel mask bits treated as bad when set in every exposure (and_mask):
#: FULLREJECT (18), NOSKY (22), BRIGHTSKY (23), NODATA (24).
BAD_BITS = (1 << 18) | (1 << 22) | (1 << 23) | (1 << 24)
#: Scaling ignores the forest: Vanden Berk et al. used no flux blueward of Ly alpha.
SCALE_MIN_WAVE = 1230.0
GM_WINDOW = (3020.0, 3100.0)
#: Continuum windows of Vanden Berk et al. 2001 (section 4.1).
VB01_UV_WINDOWS = ((1350.0, 1365.0), (4200.0, 4230.0))
VB01_RED_WINDOWS = ((6005.0, 6035.0), (7160.0, 7180.0))


def ccm89_alambda_over_av(wave_angstrom, r_v: float = 3.1):
    """Cardelli, Clayton & Mathis (1989) extinction curve A_lambda / A_V, 0.3 <= x <= 3.3 /um."""
    x = 1e4 / np.asarray(wave_angstrom, float)
    a = np.zeros_like(x)
    b = np.zeros_like(x)
    ir = (x >= 0.3) & (x < 1.1)
    a[ir] = 0.574 * x[ir] ** 1.61
    b[ir] = -0.527 * x[ir] ** 1.61
    opt = (x >= 1.1) & (x <= 3.3)
    y = x[opt] - 1.82
    a[opt] = (1 + 0.17699 * y - 0.50447 * y**2 - 0.02427 * y**3 + 0.72085 * y**4
              + 0.01979 * y**5 - 0.77530 * y**6 + 0.32999 * y**7)
    b[opt] = (1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3 - 5.38434 * y**4
              - 0.62251 * y**5 + 5.30260 * y**6 - 2.09002 * y**7)
    return a + b / r_v


def rebin_flux(edges_in, values, edges_out):
    """Flux-conserving rebinning of a density sampled in bins ``edges_in`` onto ``edges_out``."""
    cum = np.concatenate([[0.0], np.cumsum(values * np.diff(edges_in))])
    c_out = np.interp(edges_out, edges_in, cum, left=np.nan, right=np.nan)
    return np.diff(c_out) / np.diff(edges_out)


def read_boss_spectrum(path):
    """Observed wavelength (A), f_lambda (1e-17 erg/s/cm2/A), inverse variance, good-pixel mask."""
    with fits.open(path, memmap=False) as h:
        d = h[1].data
        wave = 10.0 ** d["loglam"].astype(float)
        flux = d["flux"].astype(float)
        ivar = d["ivar"].astype(float)
        and_mask = d["and_mask"].astype(np.int64)
    good = (ivar > 0) & ((and_mask & BAD_BITS) == 0) & np.isfinite(flux)
    good &= (wave >= BOSS_TRUSTED[0]) & (wave <= BOSS_TRUSTED[1])
    return wave, flux, ivar, good


def to_rest_grid(path, z: float, ebv: float, grid=COMPOSITE_GRID):
    """One spectrum dereddened, shifted to the rest frame and rebinned onto ``grid``.

    Returns (flux, ivar) on the grid, NaN where fewer than 90% of the input
    pixels were good or the spectrum does not reach.
    """
    wave, flux, ivar, good = read_boss_spectrum(path)
    ext = 10.0 ** (0.4 * 3.1 * ebv * ccm89_alambda_over_av(wave))
    flux = flux * ext
    var = np.where(good, 1.0 / np.where(ivar > 0, ivar, 1.0), 0.0) * ext**2
    # bin edges in the observed frame (BOSS pixels are uniform in log lambda)
    loglam = np.log10(wave)
    mid = 0.5 * (loglam[1:] + loglam[:-1])
    edges_obs = 10.0 ** np.concatenate([[loglam[0] - (mid[0] - loglam[0])], mid,
                                        [loglam[-1] + (loglam[-1] - mid[-1])]])
    edges_rest = edges_obs / (1.0 + z)
    g_edges = np.concatenate([[grid[0] - 0.5], 0.5 * (grid[1:] + grid[:-1]), [grid[-1] + 0.5]])
    gf = good.astype(float)
    cover = rebin_flux(edges_rest, gf, g_edges)
    fsum = rebin_flux(edges_rest, np.where(good, flux, 0.0), g_edges)
    # Variance of the overlap-weighted mean, sum_j var_j o_j^2 / (sum_j o_j)^2, with
    # o_j^2 approximated by w_j o_j (exact for input pixels wholly inside an output
    # pixel, up to 2x high for split ones): integrating the density var_j * w_j
    # gives sum_j var_j w_j o_j. Only used for diagnostics; the stack is unweighted.
    w_in = np.diff(edges_rest)
    cum_v = np.concatenate([[0.0], np.cumsum(var * w_in**2)])
    v_sum = np.diff(np.interp(g_edges, edges_rest, cum_v, left=np.nan, right=np.nan))
    out_w = np.diff(g_edges)
    with np.errstate(invalid="ignore", divide="ignore"):
        f_out = fsum / cover
        var_out = v_sum / (cover * out_w) ** 2
    ok = np.isfinite(cover) & (cover > 0.9)
    f_out = np.where(ok, f_out, np.nan)
    ivar_out = np.where(ok & (var_out > 0), 1.0 / var_out, np.nan)
    return f_out.astype(np.float32), ivar_out.astype(np.float32)


def load_stack(sample: pd.DataFrame, paths, workers: int = 8):
    """Rest-frame fluxes of every spectrum, shape (n_spectra, n_grid)."""
    n, m = len(sample), len(COMPOSITE_GRID)
    flux = np.full((n, m), np.nan, np.float32)
    ivar = np.full((n, m), np.nan, np.float32)
    args = [(p, float(z), float(e)) for p, z, e in zip(paths, sample["Z_SYS"], sample["EBV"])]

    def work(k):
        p, z, e = args[k]
        try:
            return k, to_rest_grid(p, z, e)
        except (OSError, ValueError, KeyError):
            return k, None

    bad = []
    with cf.ThreadPoolExecutor(max_workers=workers) as pool:
        for k, res in pool.map(work, range(n)):
            if res is None:
                bad.append(k)
                continue
            flux[k], ivar[k] = res
    return flux, ivar, bad


def _scale_to(spec, comp, wave, min_wave=SCALE_MIN_WAVE, min_pix: int = 200):
    m = np.isfinite(spec) & np.isfinite(comp) & (wave >= min_wave) & (comp > 0)
    if m.sum() < min_pix:
        return np.nan
    return float(np.median(spec[m] / comp[m]))


def median_composite(flux, wave=COMPOSITE_GRID, z=None, n_iter: int = 6):
    """Iteratively scaled median composite.

    Start: spectra covering 3020-3100 A normalised there (that window is seen by
    most of the sample); spectra that miss it are then scaled to the provisional
    composite in their overlap, nearest redshifts first. Then every spectrum is
    rescaled to the full composite and the median recomputed, a few times over.
    """
    n = flux.shape[0]
    scale = np.full(n, np.nan)
    win = (wave >= GM_WINDOW[0]) & (wave <= GM_WINDOW[1])
    cov = np.isfinite(flux[:, win]).mean(axis=1) > 0.9
    scale[cov] = np.nanmean(flux[cov][:, win], axis=1)
    scale[scale <= 0] = np.nan
    order = np.argsort(np.abs(np.log1p(z) - np.median(np.log1p(z[cov])))) if z is not None else np.arange(n)
    with np.errstate(all="ignore"):
        comp = np.nanmedian(flux[np.isfinite(scale)] / scale[np.isfinite(scale), None], axis=0)
    for k in order:  # chain the rest on, nearest redshifts first
        if np.isfinite(scale[k]):
            continue
        s = _scale_to(flux[k], comp, wave)
        if np.isfinite(s) and s > 0:
            scale[k] = s
            with np.errstate(all="ignore"):
                comp = np.where(np.isfinite(comp), comp, flux[k] / s)
    history = []
    for _ in range(n_iter):
        with np.errstate(all="ignore"):
            comp = np.nanmedian(flux / scale[:, None], axis=0)
        new = np.array([_scale_to(flux[k], comp, wave) for k in range(n)])
        new[new <= 0] = np.nan
        change = np.nanmedian(np.abs(np.log(new / scale)))
        history.append(float(change))
        scale = new
    norm = flux / scale[:, None]
    with np.errstate(all="ignore"):
        comp = np.nanmedian(norm, axis=0)
        q16, q84 = np.nanpercentile(norm, [16, 84], axis=0)
    nspec = np.isfinite(norm).sum(axis=0)
    err = 0.5 * (q84 - q16) / np.sqrt(np.maximum(nspec, 1))
    return dict(flux=comp, err=err, n=nspec, scale=scale, history=history)


def geometric_mean_composite(flux, wave=COMPOSITE_GRID):
    """Spectra covering 3020-3100 A, normalised to unit mean there, averaged in log."""
    win = (wave >= GM_WINDOW[0]) & (wave <= GM_WINDOW[1])
    cov = np.isfinite(flux[:, win]).mean(axis=1) > 0.9
    sub = flux[cov]
    norm = np.nanmean(sub[:, win], axis=1)
    keep = norm > 0
    sub = sub[keep] / norm[keep, None]
    with np.errstate(all="ignore"):
        logf = np.where(sub > 0, np.log(sub), np.nan)
        gm = np.exp(np.nanmean(logf, axis=0))
    n = np.isfinite(logf).sum(axis=0)
    nonpos = float(np.mean(sub[np.isfinite(sub)] <= 0))
    return dict(flux=gm, n=n, n_spectra=int(keep.sum()), frac_nonpositive=nonpos)


def power_law_fit(wave, flux, windows):
    """alpha_nu of a power law through the mean flux of two continuum windows.

    f_lambda ~ lambda^alpha_lambda with alpha_nu = -(alpha_lambda + 2).
    Returns (alpha_nu, normalisation function).
    """
    pts = []
    for lo, hi in windows:
        m = (wave >= lo) & (wave <= hi) & np.isfinite(flux)
        pts.append((np.log10(np.mean(wave[m])), np.log10(np.mean(flux[m]))))
    (x1, y1), (x2, y2) = pts
    alpha_lambda = (y2 - y1) / (x2 - x1)
    alpha_nu = -(alpha_lambda + 2.0)

    def model(w):
        return 10.0 ** (y1 + alpha_lambda * (np.log10(w) - x1))

    return float(alpha_nu), model
