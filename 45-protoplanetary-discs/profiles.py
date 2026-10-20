"""Radial brightness profiles, and the gaps and rings in them.

Once a disc is deprojected, a ring is a circle, so averaging the brightness around each
circle beats the noise down and turns the image into a one-dimensional profile I(r).
Gaps are its local minima ("D" for dark, in DSHARP's naming, followed by the radius in
au), rings its local maxima ("B" for bright).

Widths and depths follow DSHARP II section 3.2 exactly. For a gap D followed outward by
its ring B, I_mean = (I_D + I_B)/2; the gap's outer edge is where the profile climbs
through I_mean between D and B, its inner edge the last crossing inside D, and the ring's
outer edge the first crossing outside B. The gap width over the outer edge's radius is
DSHARP VII's Delta = (r_out - r_in)/r_out, the number that becomes a planet mass.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage, signal

from alma_images import Image
from geometry import DiscGeometry, sky_from_polar


@dataclass
class Profile:
    r: np.ndarray  # bin centres, arcsec
    intensity: np.ndarray  # mean, Jy/beam
    error: np.ndarray  # scatter / sqrt(number of beams around the arc), Jy/beam
    scatter: np.ndarray  # standard deviation of the pixels in the bin
    npix: np.ndarray
    nbeam: np.ndarray
    bin_width: float  # arcsec
    wedge_deg: float | None = None
    exclude: tuple = ()

    def at(self, r):
        return np.interp(r, self.r, self.intensity)


def azimuth_mask(theta_deg: np.ndarray, wedge_deg: float | None, exclude: tuple) -> np.ndarray:
    """Pixels kept for the average: all azimuths, or only those within +/- wedge of the
    major axis (theta = +/-90 deg); minus excluded ranges."""
    keep = np.ones(theta_deg.shape, dtype=bool)
    if wedge_deg is not None:
        keep = (np.abs(theta_deg - 90.0) <= wedge_deg) | (np.abs(theta_deg + 90.0) <= wedge_deg)
    for lo, hi in exclude:
        keep &= ~((theta_deg > lo) & (theta_deg < hi))
    return keep


def _arc_length(r: float, geom: DiscGeometry, wedge_deg, exclude, n: int = 1441) -> float:
    """Length on the sky (arcsec) of the part of the projected circle that is averaged."""
    th = np.linspace(-180.0, 180.0, n)
    x, y = sky_from_polar(r, th, 0.0, 0.0, geom.inc, geom.pa)
    seg = np.hypot(np.diff(x), np.diff(y))
    mid = 0.5 * (th[1:] + th[:-1])
    return float(np.sum(seg[azimuth_mask(mid, wedge_deg, exclude)]))


def radial_profile(img: Image, geom: DiscGeometry, bin_width: float, r_max: float,
                   wedge_deg: float | None = None, exclude: tuple = ()) -> Profile:
    xx, yy = img.grids()
    r, th = geom.polar(xx, yy)
    keep = azimuth_mask(th, wedge_deg, exclude) & (r < r_max)
    edges = np.arange(0.0, r_max + bin_width, bin_width)
    idx = np.digitize(r[keep], edges) - 1
    vals = img.data[keep]
    nb = edges.size - 1
    cnt = np.bincount(idx, minlength=nb)[:nb]
    s1 = np.bincount(idx, weights=vals, minlength=nb)[:nb]
    s2 = np.bincount(idx, weights=vals * vals, minlength=nb)[:nb]
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = s1 / cnt
        sd = np.sqrt(np.maximum(s2 / cnt - mean**2, 0.0))
    centres = 0.5 * (edges[1:] + edges[:-1])
    nbeam = np.array([max(_arc_length(rc, geom, wedge_deg, exclude) / img.beam_mean, 1.0) for rc in centres])
    return Profile(centres, mean, sd / np.sqrt(nbeam), sd, cnt, nbeam, bin_width, wedge_deg, exclude)


# ---------------------------------------------------------------------------
# Finding gaps and rings
# ---------------------------------------------------------------------------
@dataclass
class Extremum:
    kind: str  # "D" gap (minimum) or "B" ring (maximum)
    r: float  # arcsec, refined with a parabola through the three bins around the extremum
    intensity: float
    prominence: float  # in units of the profile error
    contrast: float  # prominence as a fraction of the local intensity


def find_extrema(prof: Profile, beam: float, r_min: float, min_sigma: float = 5.0,
                 min_contrast: float = 0.01, noise: float | None = None) -> list[Extremum]:
    """Local minima and maxima of the profile that stand out by > min_sigma times its
    error (and by > min_contrast of the local brightness), at least half a beam apart."""
    i = prof.intensity.copy()
    good = np.isfinite(i)
    i[~good] = 0.0
    err = np.where(np.isfinite(prof.error) & (prof.error > 0), prof.error, np.inf)
    if noise is not None:  # floor: the image noise averaged over the beams of the arc
        err = np.maximum(err, noise / np.sqrt(prof.nbeam))
    dist = max(1, int(0.5 * beam / prof.bin_width))
    out = []
    for kind, sgn in (("B", 1.0), ("D", -1.0)):
        peaks, props = signal.find_peaks(sgn * i, prominence=0.0, distance=dist)
        for p, prom in zip(peaks, props["prominences"]):
            if prof.r[p] < r_min or not good[p] or p == 0 or p == i.size - 1:
                continue
            s = prom / err[p]
            c = prom / max(abs(i[p]), 1e-30)
            if s < min_sigma or c < min_contrast:
                continue
            y1, y2, y3 = i[p - 1], i[p], i[p + 1]
            den = y1 - 2 * y2 + y3
            shift = 0.5 * (y1 - y3) / den if den != 0 else 0.0
            out.append(Extremum(kind, float(prof.r[p] + np.clip(shift, -1, 1) * prof.bin_width), float(i[p]),
                                float(s), float(c)))
    return sorted(out, key=lambda e: e.r)


# ---------------------------------------------------------------------------
# Widths and depths (DSHARP II sec. 3.2 and Appendix B; DSHARP VII sec. 3.3.2)
# ---------------------------------------------------------------------------
@dataclass
class GapRing:
    gap: str  # labels, e.g. "D61"
    ring: str
    gap_r: float  # arcsec
    ring_r: float
    i_gap: float
    i_ring: float
    i_mean: float
    r_gap_in: float  # inner edge of the gap
    r_gap_out: float  # boundary between gap and ring
    r_ring_out: float  # outer edge of the ring
    depth: float  # I_D / I_B
    gap_width: float  # arcsec
    ring_width: float
    delta: float  # DSHARP VII: (r_out - r_in) / r_out
    case: str = "1"


def _crossings(r, i, level):
    """Radii where the (finely interpolated) profile crosses ``level``. A sample lying
    exactly on the level counts as above it, so such a crossing is not lost."""
    d = np.asarray(i, dtype=float) - level
    s = np.where(d >= 0, 1.0, -1.0)
    k = np.where(s[:-1] * s[1:] < 0)[0]
    return r[k] + (level - i[k]) * (r[k + 1] - r[k]) / (i[k + 1] - i[k])


def pair_measures(prof: Profile, features: list[tuple[str, str, float]], step: float) -> list[GapRing]:
    """Widths and depths of every gap/ring pair of a disc.

    ``features`` is the disc's ordered list of (label, kind, radius) with kind "D" or "B".
    A pair is a gap and the ring immediately outside it. Special cases (DSHARP II App. B):
    case 2, the ring inside the gap peaks below I_mean -> the gap's inner edge is that
    ring's outer edge; case 3, the gap outside the ring stays above I_mean -> the ring's
    outer edge is that gap's inner edge.
    """
    feats = sorted(features, key=lambda f: f[2])
    rf = np.arange(prof.r[0], prof.r[-1], step)
    ok = np.isfinite(prof.intensity)
    i_f = np.interp(rf, prof.r[ok], prof.intensity[ok])
    binval = lambda r: float(prof.intensity[np.argmin(np.abs(prof.r - r))])  # noqa: E731

    pairs = []
    for k in range(len(feats) - 1):
        (lab_d, kd, rd), (lab_b, kb, rb) = feats[k], feats[k + 1]
        if kd != "D" or kb != "B":
            continue
        i_d, i_b = binval(rd), binval(rb)
        if not i_b > i_d:
            continue
        lev = 0.5 * (i_d + i_b)
        cr = _crossings(rf, i_f, lev)
        between = cr[(cr > rd) & (cr < rb)]
        if between.size == 0:
            continue
        inner = feats[k - 1] if k >= 1 else None
        outer = feats[k + 2] if k + 2 < len(feats) else None
        lo_lim = inner[2] if inner is not None else 0.0
        hi_lim = outer[2] if outer is not None else np.inf
        inside = cr[(cr < rd) & (cr > lo_lim)]
        outside = cr[(cr > rb) & (cr < hi_lim)]
        pairs.append(dict(k=k, gap=lab_d, ring=lab_b, rd=rd, rb=rb, i_d=i_d, i_b=i_b, lev=lev, r_do=float(between[0]),
                          r_di=float(inside[-1]) if inside.size else np.nan,
                          r_bo=float(outside[0]) if outside.size else np.nan, case="1",
                          inner=inner, outer=outer))
    by_k = {p["k"]: p for p in pairs}
    for p in pairs:  # case 2
        inner = p["inner"]
        if inner is not None and inner[1] == "B" and binval(inner[2]) < p["lev"]:
            prev = by_k.get(p["k"] - 2)
            p["r_di"] = prev["r_bo"] if prev is not None else np.nan
            p["case"] = "2"
    for p in pairs:  # case 3
        outer = p["outer"]
        if outer is not None and outer[1] == "D" and binval(outer[2]) > p["lev"]:
            nxt = by_k.get(p["k"] + 2)
            p["r_bo"] = nxt["r_di"] if (nxt is not None and nxt["case"] == "1") else np.nan
            p["case"] = "3" if p["case"] == "1" else p["case"] + "+3"
    out = []
    for p in pairs:
        gw = p["r_do"] - p["r_di"]
        rw = p["r_bo"] - p["r_do"]
        out.append(GapRing(p["gap"], p["ring"], p["rd"], p["rb"], p["i_d"], p["i_b"], p["lev"], p["r_di"], p["r_do"],
                           p["r_bo"], p["i_d"] / p["i_b"], gw, rw, gw / p["r_do"], p["case"]))
    return out


def correlated_noise(n: int, sigma_bins: float, rng: np.random.Generator, size: int) -> np.ndarray:
    """Unit-variance Gaussian noise correlated over ``sigma_bins`` bins (the beam), to
    perturb a profile realistically in Monte Carlo error estimates."""
    w = rng.standard_normal((size, n))
    if sigma_bins > 0.3:
        w = ndimage.gaussian_filter1d(w, sigma_bins, axis=1, mode="nearest")
        w /= w.std(axis=1, keepdims=True)
    return w


def pair_measures_with_errors(prof: Profile, features, step: float, beam: float, rng: np.random.Generator,
                              n_mc: int = 300) -> tuple[list[GapRing], dict]:
    """:func:`pair_measures` plus 1-sigma errors from Monte Carlo perturbations of the
    profile by its own errors, correlated over a beam."""
    best = pair_measures(prof, features, step)
    sig = beam / 2.3548 / prof.bin_width
    noise = correlated_noise(prof.r.size, sig, rng, n_mc)
    keys = ("gap_width", "ring_width", "depth", "delta", "r_gap_out", "r_gap_in", "r_ring_out")
    samples = {(m.gap, m.ring): {k: [] for k in keys} for m in best}
    err = np.nan_to_num(prof.error)
    for j in range(n_mc):
        p2 = Profile(prof.r, prof.intensity + noise[j] * err, prof.error, prof.scatter, prof.npix, prof.nbeam,
                     prof.bin_width)
        for m in pair_measures(p2, features, step):
            if (m.gap, m.ring) in samples:
                for k in keys:
                    samples[(m.gap, m.ring)][k].append(getattr(m, k))
    errors = {}
    for pair, d in samples.items():
        errors[pair] = {k: (float(np.nanstd(v)) if np.sum(np.isfinite(v)) > 20 else np.nan) for k, v in d.items()}
    return best, errors
