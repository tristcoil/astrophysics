"""Turn archive photometry into clean relative-brightness time series (ppm).

Three kinds of input share one recipe:

* Kepler long cadence (29.4 min) of the red giants, one FITS file per quarter;
* Kepler short cadence (58.8 s) of 16 Cyg A and B, one file per month;
* the Sun: SOHO/VIRGO's green-channel sun photometer (1-minute samples).

The recipe, per continuous piece of data (a Kepler quarter or month):

1. keep cadences Kepler did not flag as bad (lightkurve's "default" quality
   bitmask: attitude tweaks, safe modes, coarse pointing, Earth pointing,
   desaturations, manual excludes, detector anomalies, thruster firings);
2. divide by the median: Kepler rolled 90 degrees every quarter and each
   quarter landed the star on a different detector, so raw levels jump;
3. divide by a smooth trend (a Gaussian-weighted running mean, sigma =
   ``TREND_SIGMA_DAYS``) to remove slow instrumental drifts and starspots,
   computed as a *normalised* convolution so that gaps do not drag it down;
4. clip points more than ``CLIP_SIGMA`` robust standard deviations out.

The trend filter is a high-pass filter with a known transfer function,
``filter_transfer``, which the background model includes.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d

import config as C

#: lightkurve's ``KeplerQualityFlags.DEFAULT_BITMASK``: attitude tweak, safe
#: mode, coarse point, Earth point, desaturation, manual exclude, detector
#: anomaly, no data, thruster firing.
QUALITY_BITMASK = 1130799


@dataclass
class TimeSeries:
    """A prepared light curve: time in days, relative flux in ppm."""

    time: np.ndarray  # days
    flux: np.ndarray  # ppm
    cadence_s: float
    trend_sigma_days: float
    label: str
    n_segments: int

    @property
    def span_days(self) -> float:
        return float(self.time.max() - self.time.min())

    @property
    def duty_cycle(self) -> float:
        n_possible = self.span_days * 86400.0 / self.cadence_s + 1
        return float(len(self.time) / n_possible)


def filter_transfer(nu_uhz, sigma_days: float) -> np.ndarray:
    """Power transfer of 'divide by a Gaussian running mean' at frequency nu.

    Subtracting a Gaussian-smoothed copy multiplies the signal's Fourier
    amplitude by 1 - exp(-2 pi^2 sigma^2 nu^2); power goes as the square.
    """
    s = sigma_days * 86400.0
    nu = np.asarray(nu_uhz) * 1e-6
    return (1.0 - np.exp(-2.0 * np.pi ** 2 * s ** 2 * nu ** 2)) ** 2


def _detrend_segment(cad: np.ndarray, rel: np.ndarray, sigma_cad: float) -> np.ndarray:
    """Normalised Gaussian convolution on the regular cadence grid."""
    idx = cad - cad.min()
    n = int(idx.max()) + 1
    y = np.zeros(n)
    w = np.zeros(n)
    y[idx] = rel
    w[idx] = 1.0
    num = gaussian_filter1d(y, sigma_cad, mode="constant", truncate=4.0)
    den = gaussian_filter1d(w, sigma_cad, mode="constant", truncate=4.0)
    return num[idx] / den[idx]


def _robust_clip(flux: np.ndarray, nsig: float, iterations: int = 3) -> np.ndarray:
    keep = np.isfinite(flux)
    for _ in range(iterations):
        med = np.median(flux[keep])
        mad = 1.4826 * np.median(np.abs(flux[keep] - med))
        new = keep & (np.abs(flux - med) < nsig * mad)
        if new.sum() == keep.sum():
            break
        keep = new
    return keep


def prepare_segments(segments, *, cadence_s: float, trend_sigma_days: float,
                     clip_sigma: float = C.CLIP_SIGMA, label: str = "") -> TimeSeries:
    """Steps 2-4 of the recipe for a list of (time_days, flux, cadence_index)."""
    sigma_cad = trend_sigma_days * 86400.0 / cadence_s
    ts, fs = [], []
    for t, f, cad in segments:
        ok = np.isfinite(t) & np.isfinite(f) & (f > 0)
        t, f, cad = t[ok], f[ok], cad[ok]
        if len(t) < 50:
            continue
        rel = f / np.median(f)
        trend = _detrend_segment(cad, rel, sigma_cad)
        ppm = (rel / trend - 1.0) * 1e6
        ts.append(t)
        fs.append(ppm)
    t = np.concatenate(ts)
    f = np.concatenate(fs)
    order = np.argsort(t)
    t, f = t[order], f[order]
    keep = _robust_clip(f, clip_sigma)
    f = f[keep] - np.mean(f[keep])
    return TimeSeries(time=t[keep], flux=f, cadence_s=cadence_s,
                      trend_sigma_days=trend_sigma_days, label=label, n_segments=len(ts))


# ---------------------------------------------------------------------------
# Kepler
# ---------------------------------------------------------------------------


def read_kepler_file(path: Path):
    """(time [BJD-2454833], PDCSAP flux, cadence number, quarter) of one file."""
    from astropy.io import fits

    with fits.open(path, memmap=False) as h:
        d = h[1].data
        quarter = int(h[0].header.get("QUARTER", -1))
        t = np.asarray(d["TIME"], float)
        f = np.asarray(d["PDCSAP_FLUX"], float)
        q = np.asarray(d["SAP_QUALITY"], np.int64)
        cad = np.asarray(d["CADENCENO"], np.int64)
    good = (q & QUALITY_BITMASK) == 0
    return t[good], f[good], cad[good], quarter


def kepler_long_cadence(kic: int) -> TimeSeries:
    """All long-cadence quarters of one star, prepared and stitched."""
    from fetch_data import kepler_dir

    files = sorted(kepler_dir(kic).glob("*_llc.fits"))
    if not files:
        raise FileNotFoundError(f"no Kepler light curves cached for KIC {kic}")
    segs = []
    for p in files:
        t, f, cad, _ = read_kepler_file(p)
        segs.append((t, f, cad))
    return prepare_segments(segs, cadence_s=C.KEPLER_LC_CADENCE_S,
                            trend_sigma_days=C.TREND_SIGMA_DAYS_LC, label=f"KIC {kic}")


def kepler_short_cadence(name: str) -> TimeSeries:
    """16 Cyg A or B: monthly short-cadence files from the custom apertures."""
    folder = C.KEPLER_DIR / name.replace(" ", "_").lower()
    files = sorted(folder.glob("*_slc.fits"))
    if not files:
        raise FileNotFoundError(f"no short-cadence files cached for {name}")
    segs = []
    for p in files:
        t, f, cad, _ = read_kepler_file(p)
        segs.append((t, f, cad))
    return prepare_segments(segs, cadence_s=C.KEPLER_SC_CADENCE_S,
                            trend_sigma_days=C.TREND_SIGMA_DAYS_SC, label=name)


# ---------------------------------------------------------------------------
# The Sun (SOHO/VIRGO SPM green channel)
# ---------------------------------------------------------------------------


def read_virgo(path: Path):
    """Time (days since 1996-01-01 TAI-ish, as in the file) and flux of VIRGO.

    Returns the raw columns; ``virgo_segments`` cuts them into pieces.
    """
    from astropy.io import fits

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with fits.open(path, memmap=True) as h:
            hdu = h[1] if len(h) > 1 and h[1].data is not None else h[0]
            data = hdu.data
            cols = [c.lower() for c in data.columns.names] if hasattr(data, "columns") else []
            if cols:
                return {c: np.asarray(data[n], float) for c, n in zip(cols, data.columns.names)}, \
                    dict(hdu.header), dict(h[0].header)
            return {"flux": np.asarray(data, float)}, dict(hdu.header), dict(h[0].header)
