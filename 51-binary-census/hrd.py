"""Gaia's binaries on the Hertzsprung-Russell diagram (the 100 pc sample).

An unresolved binary appears as one dot with the summed light of both stars.
A faint companion (q < ~0.6) adds little: the dot stays on the main sequence.
An equal-mass twin doubles the light: 2.5 log10(2) = 0.753 mag brighter at the
same colour. B1 found that twin sequence among K and early-M dwarfs; here we ask
which of Gaia's detection channels populate it.

Each channel has its own photometric signature, because each is sensitive to a
different mass ratio:

* SB2 (both stars' lines visible) needs a bright companion: high q, elevated;
* SB1 and astrometric orbits need a faint companion (a twin's light centre and
  blended lines hardly move): low q, near the single-star ridge;
* RUWE > 1.4 without a published solution: everything in between.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

import binaries as B
import samples


def empirical_ridge(colour, abs_g, single_like, *, edges=None, smooth_mag=0.05):
    """Most common M_G in narrow colour bins among single-looking stars,
    smoothed across colour: the single-star main sequence of this sample."""
    if edges is None:
        edges = np.arange(0.3, 4.21, 0.05)
    centres = 0.5 * (edges[1:] + edges[:-1])
    ridge = np.full(len(centres), np.nan)
    mag_edges = np.arange(-2, 17, 0.02)
    for k in range(len(centres)):
        m = single_like & (colour >= edges[k]) & (colour < edges[k + 1])
        if m.sum() < 30:
            continue
        # look near Mamajek's sequence (+-1 mag) for the densest point
        h, _ = np.histogram(abs_g[m], mag_edges)
        hs = ndimage.gaussian_filter1d(h.astype(float), smooth_mag / 0.02)
        ridge[k] = 0.5 * (mag_edges[np.argmax(hs)] + mag_edges[np.argmax(hs) + 1])
    ok = np.isfinite(ridge)
    # light smoothing across colour bins (3-bin running median)
    r = ridge.copy()
    for k in np.where(ok)[0]:
        lo, hi = max(k - 1, 0), min(k + 2, len(r))
        r[k] = np.nanmedian(ridge[lo:hi])
    return centres[ok], r[ok]


def offsets(t, centres, ridge):
    """M_G minus the ridge at each star's colour (negative = brighter)."""
    return t["abs_g"] - np.interp(t["bp_rp"], centres, ridge, left=np.nan, right=np.nan)


def twin_peak(delta, *, window=(-1.1, -0.5), bin_width=0.02, smooth=0.04):
    """Position of the secondary (twin) peak on the bright side of the ridge."""
    edges = np.arange(-2.0, 1.0 + bin_width, bin_width)
    h, _ = np.histogram(delta[np.isfinite(delta)], edges)
    hs = ndimage.gaussian_filter1d(h.astype(float), smooth / bin_width)
    c = 0.5 * (edges[1:] + edges[:-1])
    inside = np.where((c >= window[0]) & (c <= window[1]))[0]
    k = inside[np.argmax(hs[inside])]
    return float(c[k]), edges, h


def channel_labels(t) -> np.ndarray:
    """One label per star for the HR diagram, in order of precedence."""
    lab = np.full(len(t), "none", dtype="U14")
    ruwe = np.ma.filled(np.ma.asarray(t["ruwe"], dtype=float), np.nan)
    lab[ruwe > 1.4] = "ruwe"
    lab[(t["accel_type"] != "") | (t["trend_type"] != "")] = "acceleration"
    lab[np.isin(t["nss_type"], ("Orbital", "OrbitalTargetedSearch", "OrbitalTargetedSearchValidated",
                                "OrbitalAlternative", "OrbitalAlternativeValidated"))] = "astrometric"
    lab[np.isin(t["nss_type"], ("SB1", "SB1C", "AstroSpectroSB1"))] = "sb1"
    lab[np.isin(t["nss_type"], ("SB2", "SB2C"))] = "sb2"
    lab[np.isin(t["nss_type"], ("EclipsingBinary", "EclipsingSpectro"))] = "eclipsing"
    return lab


def twin_offset_predicted() -> float:
    return -2.5 * np.log10(2.0)


def mass_ratio_from_offset(delta, m1_grid=None):
    """Not used for science: a reminder of how little a low-mass companion adds.
    Returns the brightening (mag) of a 1 Msun star for q = 0.1 ... 1."""
    import stellar
    q = np.linspace(0.1, 1.0, 10)
    f = stellar.flux_ratio_g(1.0, q)
    return q, -2.5 * np.log10(1 + f)


__all__ = ["empirical_ridge", "offsets", "twin_peak", "channel_labels", "samples", "B"]
