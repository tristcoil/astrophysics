"""Part 1: GD-1, a star cluster torn into an 80-degree-long stream.

The pipeline is told only where GD-1 lies on the sky (the stream frame of
Koposov et al. 2010, which follows the stream's discovered track). Everything
else is found in Gaia DR3:

1. **The template.** The globular cluster M92 is as old and metal-poor as
   GD-1's parent must have been, and it lies at the same distance, so its
   colour-magnitude diagram -- measured by Gaia with Gaia's own noise at the
   same apparent magnitudes -- is the filter. This is how GD-1 was found in
   2006 (Grillmair & Dionatos used a globular cluster's diagram to build a
   "matched filter" for SDSS photometry).
2. **The motion.** Stars on the stream track and in control bands beside it
   are compared in proper-motion space, piece by piece along the stream. The
   stream shows up as a knot on the track that the control bands lack. A
   smooth curve through the knots is the stream's proper-motion track.
3. **The distance.** Stars moving with the stream, minus the control bands'
   stars, form a colour-magnitude sequence. Sliding M92's sequence up and down
   to match it measures the distance modulus, segment by segment.
4. **Members.** Stars matching the motion, the colour-magnitude sequence at
   that distance, and a parallax consistent with being far away. Then: the
   stream's track, its density along its length (gaps), and stars off the
   main track (the "spur").
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from astropy.table import Table
from scipy.interpolate import UnivariateSpline
from scipy.ndimage import convolve

import galaxy as gx
import photometry as ph
from config import (
    GD1_OFF_STREAM_PHI2, GD1_ON_STREAM_PHI2, GD1_PHI1_RANGE, GD1_PHI2_RANGE,
    GD1_PM_BIN, GD1_PM_WINDOW, GD1_SEARCH_SEGMENT_DEG, M92, REFERENCE, RUWE_MAX,
)
from fetch_data import GD1_COUNTS_PATH, GD1_PATH, HUANG_PATH, M92_PATH

# --------------------------------------------------------------------------
# Tunable numbers of the method (documented in the README)
# --------------------------------------------------------------------------

#: The template's colour-magnitude sequence: magnitude bins with at least this
#: many M92 stars each, but no wider than MAX_BIN_MAG.
TEMPLATE_MIN_STARS = 40
TEMPLATE_MIN_BIN_MAG = 0.08
TEMPLATE_MAX_BIN_MAG = 0.5
#: Width floor added in quadrature to M92's measured sequence width: a real
#: population's spread in age/metallicity, and GD-1's small depth along the
#: line of sight.
TEMPLATE_WIDTH_FLOOR = 0.02
#: Blue horizontal branch: M92 stars bluer than this colour, and brighter than
#: HB_MAX_ABSMAG, are horizontal-branch stars.
HB_MAX_COLOUR = 0.35
HB_MAX_ABSMAG = 2.5
HB_HALF_WIDTH_MAG = 0.3
#: CMD match: within this many sequence widths of the template.
CMD_NSIGMA = 2.0
#: Parallax: GD-1 is ~8 kpc away (parallax ~0.12 mas). A star is dropped only
#: if its parallax says it is closer than 4 kpc at 3 sigma.
PARALLAX_FAR_MAS = 0.25
#: Gaia DR3 parallax systematics correlated on degree scales (Lindegren et al.
#: 2021), mas: a floor under any average over one patch of sky.
PARALLAX_SYSTEMATIC_MAS = 0.010
#: Distance moduli searched: 5.8 to 16 kpc.
DM_GRID = np.arange(13.8, 16.02, 0.02)
#: Broad CMD pre-filter used while the distance is still unknown: the template
#: at any distance modulus in this range.
DM_PREFILTER = (14.0, 15.6)
#: Proper-motion search: counts in a disc of this radius (mas/yr) on the track
#: minus the scaled control, divided by the Poisson noise.
PM_APERTURE = 0.75
#: A segment counts as a detection above this signal-to-noise.
PM_DETECT_SNR = 6.0
#: Second pass: search within this radius (mas/yr) of the first-pass track.
PM_REFINE_RADIUS = 3.0
#: Proper-motion membership: intrinsic spread of the stream around its track
#: (plus the track model's own error), mas/yr, added to each star's errors;
#: accept chi^2 < 5.99 (95% for two degrees of freedom).
PM_INTRINSIC = 0.5
PM_CHI2_MAX = 5.99
#: Degree of the polynomial tracks in phi1.
TRACK_DEGREE = 3
#: Density profile: counts in windows this wide, stepped by DENSITY_STEP.
DENSITY_WINDOW = 3.0
DENSITY_STEP = 1.0
#: Control bands for background subtraction, in |phi2 - track| (degrees):
#: beyond the spur (+0.9..1.7) and blob, to keep them out of the background.
BACKGROUND_BAND = (2.5, 5.5)


# --------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------


@dataclass
class Stars:
    """A column-oriented star table with the derived stream quantities."""

    cols: dict

    def __getitem__(self, key):
        return self.cols[key]

    def __len__(self):
        return len(next(iter(self.cols.values())))

    def subset(self, mask):
        return Stars({k: v[mask] for k, v in self.cols.items()})


def _read(path) -> dict:
    t = Table.read(path, format="fits")
    out = {}
    for c in t.colnames:
        col = t[c]
        data = np.ma.getdata(col)
        if data.dtype.kind == "f":
            values = np.asarray(data, float).copy()
            values[np.ma.getmaskarray(col)] = np.nan
        elif data.dtype.kind in "iu":
            values = np.asarray(data).astype(np.int64)
        else:
            values = np.asarray(np.ma.filled(col, "")).astype(str)
        out[c] = values
    return out


def prepare(d: dict) -> Stars:
    """Stream coordinates, dereddened photometry and quality flags."""
    phi1, phi2, pm1, pm2, m = gx.gd1_coordinates(d["ra"], d["dec"], d["pmra"], d["pmdec"])
    e1, e2, corr = gx.rotate_pm_covariance(m, d["pmra_error"], d["pmdec_error"],
                                           d["pmra_pmdec_corr"])
    ebv = ph.ebv_sfd(d["ra"], d["dec"])
    g0, c0 = ph.deredden(d["phot_g_mean_mag"], d["bp_rp"], ebv)
    good = (ph.clean_colour(d["phot_bp_rp_excess_factor"], d["bp_rp"], d["phot_g_mean_mag"])
            & (d["ruwe"] < RUWE_MAX) & np.isfinite(pm1))
    cols = dict(d)
    cols.update(phi1=phi1, phi2=phi2, pm1=pm1, pm2=pm2, pm1_err=e1, pm2_err=e2, pm_corr=corr,
                ebv=ebv, g0=g0, c0=c0, good=good)
    return Stars(cols)


def load_field() -> Stars:
    return prepare(_read(GD1_PATH))


def load_countmap():
    """The unfiltered star counts: (phi1 edges, phi2 edges, counts)."""
    from config import GD1_COUNT_MAP_BIN_DEG as b

    t = _read(GD1_COUNTS_PATH)
    e1 = np.arange(GD1_PHI1_RANGE[0], GD1_PHI1_RANGE[1] + b / 2, b)
    e2 = np.arange(GD1_PHI2_RANGE[0], GD1_PHI2_RANGE[1] + b / 2, b)
    H = np.zeros((len(e1) - 1, len(e2) - 1))
    i1 = (t["i1"] - round(e1[0] / b)).astype(int)
    i2 = (t["i2"] - round(e2[0] / b)).astype(int)
    ok = (i1 >= 0) & (i1 < H.shape[0]) & (i2 >= 0) & (i2 < H.shape[1])
    np.add.at(H, (i1[ok], i2[ok]), t["n"][ok])
    return e1, e2, H


# --------------------------------------------------------------------------
# 1. The M92 template
# --------------------------------------------------------------------------


@dataclass
class Template:
    """M92's colour-magnitude sequence in absolute magnitude."""

    absmag: np.ndarray          # bin centres of the main sequence / giant branch
    colour: np.ndarray          # median dereddened colour in each bin
    width: np.ndarray           # robust colour spread in each bin
    n: np.ndarray
    hb_colour: np.ndarray       # blue horizontal branch: colour bin centres
    hb_absmag: np.ndarray       # and the median absolute magnitude there
    members_c0: np.ndarray
    members_absmag: np.ndarray
    pm_knot: tuple
    dm: float
    n_members: int
    _fid: object = field(default=None, repr=False)
    _wid: object = field(default=None, repr=False)

    def __post_init__(self):
        # Smoothing splines weighted by the statistical error of each bin's
        # median and spread, with s = number of bins: the curve follows the
        # data to within their noise (chi^2 ~ N) and no further.
        err_med = np.maximum(1.2533 * self.width / np.sqrt(self.n), 0.002)
        self._fid = UnivariateSpline(self.absmag, self.colour, w=1 / err_med, k=3,
                                     s=len(self.absmag))
        # Width: mostly photometric noise, growing smoothly towards faint
        # magnitudes, plus the subgiant branch's tilt. A cubic in log(width),
        # weighted by each bin's star count, keeps the few giants per bin at
        # the bright end from making it wiggle.
        self._wid_coef = np.polyfit(self.absmag, np.log(np.maximum(self.width, 0.005)), 3,
                                    w=np.sqrt(self.n))
        self._wid = lambda m: np.exp(np.polyval(self._wid_coef, m))

    @property
    def absmag_range(self):
        return float(self.absmag[0]), float(self.absmag[-1])

    def colour_at(self, absmag):
        return self._fid(absmag)

    def width_at(self, absmag):
        w = np.maximum(self._wid(np.clip(absmag, *self.absmag_range)), 0.01)
        return np.sqrt(w**2 + TEMPLATE_WIDTH_FLOOR**2)

    def hb_absmag_at(self, colour):
        return np.interp(colour, self.hb_colour, self.hb_absmag, left=np.nan, right=np.nan)

    def match(self, c0, absmag, nsigma=CMD_NSIGMA):
        """True for stars on the template's sequence (giant branch, subgiant
        branch, main sequence) or its blue horizontal branch."""
        lo, hi = self.absmag_range
        inside = (absmag > lo) & (absmag < hi)
        on_seq = inside & (np.abs(c0 - self.colour_at(absmag)) < nsigma * self.width_at(absmag))
        hb = self.hb_absmag_at(c0)
        on_hb = (c0 < HB_MAX_COLOUR) & np.isfinite(hb) & (np.abs(absmag - hb) < HB_HALF_WIDTH_MAG)
        return on_seq | on_hb

    def weight(self, c0, absmag):
        """Normalised template density in colour at fixed magnitude (1/mag):
        the probability of this colour for a template star of this absolute
        magnitude. Zero off the template's magnitude range."""
        lo, hi = self.absmag_range
        w = self.width_at(absmag)
        dens = np.exp(-0.5 * ((c0 - self.colour_at(absmag)) / w) ** 2) / (np.sqrt(2 * np.pi) * w)
        return np.where((absmag > lo) & (absmag < hi), dens, 0.0)


def build_template(g_max: float) -> Template:
    d = _read(M92_PATH)
    dx = (d["ra"] - M92.ra_deg) * np.cos(np.radians(M92.dec_deg))
    r = np.hypot(dx, d["dec"] - M92.dec_deg)
    # The cluster's proper motion, found as the densest knot near its centre.
    inner = (r < 0.1) & np.isfinite(d["pmra"])
    bins = np.arange(-15, 15.001, 0.1)
    H, xe, ye = np.histogram2d(d["pmra"][inner], d["pmdec"][inner], bins=[bins, bins])
    i, j = np.unravel_index(np.argmax(H), H.shape)
    knot = (xe[i] + 0.05, ye[j] + 0.05)
    member = ((r > M92.core_exclude_deg) & (r < M92.radius_deg)
              & (np.hypot(d["pmra"] - knot[0], d["pmdec"] - knot[1]) < 1.0)
              & ph.clean_colour(d["phot_bp_rp_excess_factor"], d["bp_rp"], d["phot_g_mean_mag"])
              & (d["phot_g_mean_mag"] < g_max) & np.isfinite(d["bp_rp"]))
    ebv = ph.ebv_sfd(d["ra"][member], d["dec"][member])
    g0, c0 = ph.deredden(d["phot_g_mean_mag"][member], d["bp_rp"][member], ebv)
    dm = float(ph.distance_modulus(M92.distance_kpc))
    absmag = g0 - dm

    hb = (c0 < HB_MAX_COLOUR) & (absmag < HB_MAX_ABSMAG)
    seq = ~hb & (c0 > 0.3)
    # Sequence: equal-count bins in absolute magnitude.
    order = np.argsort(absmag[seq])
    m_s, c_s = absmag[seq][order], c0[seq][order]
    centres, colours, widths, ns = [], [], [], []
    start = 0
    while start < len(m_s):
        # Grow the bin until it has enough stars AND is wide enough, but never
        # wider than the maximum (the sparse upper giant branch).
        stop = start + 1
        while stop < len(m_s) and (stop - start < TEMPLATE_MIN_STARS
                                   or m_s[stop - 1] - m_s[start] < TEMPLATE_MIN_BIN_MAG):
            if m_s[stop] - m_s[start] > TEMPLATE_MAX_BIN_MAG:
                break
            stop += 1
        if stop - start >= 8:
            cc = c_s[start:stop]
            med = np.median(cc)
            # Iterative 3-sigma clip against stray non-members.
            for _ in range(3):
                mad = 1.4826 * np.median(np.abs(cc - med))
                keep = np.abs(cc - med) < 3 * max(mad, 0.01)
                med = np.median(cc[keep])
            centres.append(np.median(m_s[start:stop]))
            colours.append(med)
            widths.append(1.4826 * np.median(np.abs(cc[keep] - med)))
            ns.append(keep.sum())
        start = stop
    # Horizontal branch: median magnitude in colour bins.
    hb_c, hb_m = [], []
    for lo in np.arange(-0.3, HB_MAX_COLOUR, 0.05):
        s = hb & (c0 >= lo) & (c0 < lo + 0.05)
        if s.sum() >= 3:
            hb_c.append(np.median(c0[s]))
            hb_m.append(np.median(absmag[s]))
    return Template(np.array(centres), np.array(colours), np.array(widths), np.array(ns),
                    np.array(hb_c), np.array(hb_m), c0, absmag, knot, dm, int(member.sum()))


# --------------------------------------------------------------------------
# 2. The proper-motion track
# --------------------------------------------------------------------------


def _disc_kernel(radius_bins: float) -> np.ndarray:
    n = int(np.ceil(radius_bins))
    y, x = np.mgrid[-n:n + 1, -n:n + 1]
    return (np.hypot(x, y) <= radius_bins).astype(float)


def _bands(phi2, track_phi2):
    d = phi2 - track_phi2
    on = np.abs(d) < GD1_ON_STREAM_PHI2
    off = (np.abs(d) > GD1_OFF_STREAM_PHI2[0]) & (np.abs(d) < GD1_OFF_STREAM_PHI2[1])
    # Area ratio: the off bands on both sides, clipped to the field's edges.
    return on, off


def _band_area_ratio(track_phi2: float) -> float:
    lo, hi = GD1_PHI2_RANGE
    on = min(track_phi2 + GD1_ON_STREAM_PHI2, hi) - max(track_phi2 - GD1_ON_STREAM_PHI2, lo)
    a, b = GD1_OFF_STREAM_PHI2
    up = max(0.0, min(track_phi2 + b, hi) - max(track_phi2 + a, lo))
    down = max(0.0, min(track_phi2 - a, hi) - max(track_phi2 - b, lo))
    return on / (up + down)


def pm_significance_map(pm1_on, pm2_on, pm1_off, pm2_off, alpha, window=GD1_PM_WINDOW):
    """Signal-to-noise map of (on - alpha * off) counts in a disc aperture."""
    b = GD1_PM_BIN
    e1 = np.arange(window[0][0], window[0][1] + b / 2, b)
    e2 = np.arange(window[1][0], window[1][1] + b / 2, b)
    h_on, _, _ = np.histogram2d(pm1_on, pm2_on, bins=[e1, e2])
    h_off, _, _ = np.histogram2d(pm1_off, pm2_off, bins=[e1, e2])
    k = _disc_kernel(PM_APERTURE / b)
    n_on = convolve(h_on, k, mode="constant")
    n_off = convolve(h_off, k, mode="constant")
    excess = n_on - alpha * n_off
    snr = excess / np.sqrt(alpha * (1 + alpha) * n_off + 1.0)
    return 0.5 * (e1[1:] + e1[:-1]), 0.5 * (e2[1:] + e2[:-1]), snr, excess


def _cmd_prefilter(stars: Stars, tpl: Template) -> np.ndarray:
    """Stars that match the template at *some* distance in DM_PREFILTER."""
    ok = np.zeros(len(stars), bool)
    for dm in np.arange(DM_PREFILTER[0], DM_PREFILTER[1] + 1e-9, 0.1):
        ok |= tpl.match(stars["c0"], stars["g0"] - dm, nsigma=1.5)
    return ok


def far_enough(stars: Stars) -> np.ndarray:
    """Parallax consistent with lying beyond ~4 kpc (3 sigma)."""
    return stars["parallax"] - 3 * stars["parallax_error"] < PARALLAX_FAR_MAS


@dataclass
class PMTrack:
    knots: list            # one dict per searched segment
    coef1: np.ndarray      # polynomial in phi1 for pm_phi1_cosphi2
    coef2: np.ndarray      # and for pm_phi2
    maps: dict             # a few significance maps, for figures

    def at(self, phi1):
        return np.polyval(self.coef1, phi1), np.polyval(self.coef2, phi1)


def find_pm_track(stars: Stars, tpl: Template, phi2_track=None) -> PMTrack:
    """Find the stream's proper motion along its length, without being told it.

    Pass 1 searches the whole window in each segment; a smooth curve through
    the confident detections then guides pass 2, which searches near the curve
    to pick up the stream where it is fainter."""
    base = stars["good"] & far_enough(stars) & _cmd_prefilter(stars, tpl)
    lo, hi = GD1_PHI1_RANGE
    seg = GD1_SEARCH_SEGMENT_DEG
    edges = np.arange(lo, hi + 1e-9, seg / 2)  # segments overlap by half
    segments = [(a, a + seg) for a in edges if a + seg <= hi + 1e-9]
    track2 = phi2_track if phi2_track is not None else (lambda p: np.zeros_like(p))

    def search(a, b, window):
        s = base & (stars["phi1"] >= a) & (stars["phi1"] < b)
        t2 = float(np.mean(track2(np.array([0.5 * (a + b)]))))
        on, off = _bands(stars["phi2"][s], t2)
        alpha = _band_area_ratio(t2)
        x, y, snr, excess = pm_significance_map(
            stars["pm1"][s][on], stars["pm2"][s][on], stars["pm1"][s][off], stars["pm2"][s][off],
            alpha, window)
        i, j = np.unravel_index(np.nanargmax(snr), snr.shape)
        return dict(phi1=0.5 * (a + b), pm1=float(x[i]), pm2=float(y[j]), snr=float(snr[i, j]),
                    excess=float(excess[i, j]), n_on=int(on.sum())), (x, y, snr)

    knots, maps = [], {}
    for a, b in segments:
        k, m = search(a, b, GD1_PM_WINDOW)
        k["pass"] = 1
        knots.append(k)
        maps[(a, b)] = m
    good = [k for k in knots if k["snr"] > PM_DETECT_SNR]
    if len(good) < TRACK_DEGREE + 2:
        raise RuntimeError("GD-1's proper motion was not found in enough segments")
    # Robust first-pass track: drop detections far from a straight line first.
    p = np.array([k["phi1"] for k in good])
    m1 = np.array([k["pm1"] for k in good])
    m2 = np.array([k["pm2"] for k in good])
    c1, c2 = np.polyfit(p, m1, 1), np.polyfit(p, m2, 1)
    keep = np.hypot(m1 - np.polyval(c1, p), m2 - np.polyval(c2, p)) < 2.0
    c1 = np.polyfit(p[keep], m1[keep], min(TRACK_DEGREE, keep.sum() - 1))
    c2 = np.polyfit(p[keep], m2[keep], min(TRACK_DEGREE, keep.sum() - 1))

    # Pass 2: a narrow search around the first-pass track, every segment.
    knots2 = []
    for a, b in segments:
        c = 0.5 * (a + b)
        p1, p2 = np.polyval(c1, c), np.polyval(c2, c)
        w = ((p1 - PM_REFINE_RADIUS, p1 + PM_REFINE_RADIUS),
             (p2 - PM_REFINE_RADIUS, p2 + PM_REFINE_RADIUS))
        k, _ = search(a, b, w)
        k["pass"] = 2
        k["detected"] = k["snr"] > PM_DETECT_SNR
        knots2.append(k)
    det = [k for k in knots2 if k["detected"]]
    p = np.array([k["phi1"] for k in det])
    wts = np.sqrt(np.array([k["snr"] for k in det]))
    c1 = np.polyfit(p, [k["pm1"] for k in det], TRACK_DEGREE, w=wts)
    c2 = np.polyfit(p, [k["pm2"] for k in det], TRACK_DEGREE, w=wts)
    return PMTrack(knots + knots2, c1, c2, maps)


def pm_chi2(stars: Stars, track: PMTrack) -> np.ndarray:
    """Chi-squared distance of each star's proper motion from the track."""
    t1, t2 = track.at(stars["phi1"])
    d1, d2 = stars["pm1"] - t1, stars["pm2"] - t2
    s1 = stars["pm1_err"] ** 2 + PM_INTRINSIC**2
    s2 = stars["pm2_err"] ** 2 + PM_INTRINSIC**2
    c12 = stars["pm_corr"] * stars["pm1_err"] * stars["pm2_err"]
    det = s1 * s2 - c12**2
    return (s2 * d1**2 - 2 * c12 * d1 * d2 + s1 * d2**2) / det


# --------------------------------------------------------------------------
# 3. The distance along the stream
# --------------------------------------------------------------------------


def _background_density(c0, g0):
    """Colour density of the field at each magnitude (1/mag), for weighting."""
    ce = np.arange(-0.5, 1.61, 0.05)
    ge = np.arange(12.0, 21.01, 0.25)
    H, _, _ = np.histogram2d(g0, c0, bins=[ge, ce])
    H = H + 0.5  # never zero
    H /= H.sum(axis=1, keepdims=True) * 0.05
    ig = np.clip(np.digitize(g0, ge) - 1, 0, len(ge) - 2)
    ic = np.clip(np.digitize(c0, ce) - 1, 0, len(ce) - 2)
    return lambda cc, gg: H[np.clip(np.digitize(gg, ge) - 1, 0, len(ge) - 2),
                            np.clip(np.digitize(cc, ce) - 1, 0, len(ce) - 2)]


def matched_filter_signal(c0_on, g0_on, c0_off, g0_off, alpha, tpl: Template, bg):
    """Background-subtracted matched-filter signal as a function of distance
    modulus: sum over on-stream stars of template/background, minus the same
    for the scaled control bands."""
    out = np.empty(len(DM_GRID))
    b_on, b_off = bg(c0_on, g0_on), bg(c0_off, g0_off)
    for i, dm in enumerate(DM_GRID):
        s_on = np.sum(tpl.weight(c0_on, g0_on - dm) / b_on)
        s_off = np.sum(tpl.weight(c0_off, g0_off - dm) / b_off)
        out[i] = s_on - alpha * s_off
    return out


def peak_of(signal, smooth: float = 0.06):
    """Distance modulus at the peak of a signal curve: smoothed by a Gaussian
    of ``smooth`` mag, then a parabola through the three points at the top."""
    from scipy.ndimage import gaussian_filter1d

    step = DM_GRID[1] - DM_GRID[0]
    s = gaussian_filter1d(signal, smooth / step, mode="nearest")
    i = int(np.clip(np.argmax(s), 1, len(s) - 2))
    y0, y1, y2 = s[i - 1], s[i], s[i + 1]
    denom = y0 - 2 * y1 + y2
    shift = 0.5 * (y0 - y2) / denom if denom < 0 else 0.0
    return float(DM_GRID[i] + np.clip(shift, -1, 1) * step)


def fit_distance(stars: Stars, tpl: Template, track: PMTrack, phi2_track,
                 n_boot: int = 200, seed: int = 1):
    """Distance modulus in segments along the stream, with bootstrap errors."""
    rng = np.random.default_rng(seed)
    base = stars["good"] & far_enough(stars) & (pm_chi2(stars, track) < PM_CHI2_MAX)
    bg = _background_density(stars["c0"][base], stars["g0"][base])
    lo, hi = GD1_PHI1_RANGE
    seg = GD1_SEARCH_SEGMENT_DEG
    out = []
    for a in np.arange(lo, hi - seg + 1e-9, seg):
        s = base & (stars["phi1"] >= a) & (stars["phi1"] < a + seg)
        t2 = phi2_track(stars["phi1"][s])
        d = stars["phi2"][s] - t2
        on = np.abs(d) < GD1_ON_STREAM_PHI2
        off = (np.abs(d) > BACKGROUND_BAND[0]) & (np.abs(d) < BACKGROUND_BAND[1])
        alpha = _band_area_ratio_general(float(np.mean(t2)) if len(t2) else 0.0,
                                         GD1_ON_STREAM_PHI2, BACKGROUND_BAND)
        c_on, g_on = stars["c0"][s][on], stars["g0"][s][on]
        c_off, g_off = stars["c0"][s][off], stars["g0"][s][off]
        sig = matched_filter_signal(c_on, g_on, c_off, g_off, alpha, tpl, bg)
        best = peak_of(sig)
        boots = []
        for _ in range(n_boot):
            i_on = rng.integers(0, len(c_on), len(c_on))
            i_off = rng.integers(0, len(c_off), len(c_off))
            sb = matched_filter_signal(c_on[i_on], g_on[i_on], c_off[i_off], g_off[i_off],
                                       alpha, tpl, bg)
            boots.append(peak_of(sb))
        boots = np.array(boots)
        excess = len(c_on) - alpha * len(c_off)
        out.append(dict(phi1=a + seg / 2, dm=float(best),
                        dm_err=float(0.5 * (np.percentile(boots, 84) - np.percentile(boots, 16))),
                        signal=sig, n_on=int(len(c_on)), excess=float(excess)))
    return out


def _band_area_ratio_general(track_phi2, half_width, band):
    lo, hi = GD1_PHI2_RANGE
    on = min(track_phi2 + half_width, hi) - max(track_phi2 - half_width, lo)
    a, b = band
    up = max(0.0, min(track_phi2 + b, hi) - max(track_phi2 + a, lo))
    down = max(0.0, min(track_phi2 - a, hi) - max(track_phi2 - b, lo))
    return on / (up + down)


# --------------------------------------------------------------------------
# 4. Members, track, density
# --------------------------------------------------------------------------


def select_members(stars: Stars, tpl: Template, track: PMTrack, dm_of_phi1, extent) -> dict:
    """The selection, step by step (each a boolean over the field).

    The last step keeps only the stretch of sky where the stream's motion was
    actually detected: beyond it the proper-motion track would be an
    extrapolated polynomial, and near phi1 > 5 deg it runs into the field
    stars' own motion, letting in thousands of them."""
    steps = {}
    steps["quality"] = stars["good"]
    steps["parallax"] = steps["quality"] & far_enough(stars)
    steps["proper_motion"] = steps["parallax"] & (pm_chi2(stars, track) < PM_CHI2_MAX)
    absmag = stars["g0"] - dm_of_phi1(stars["phi1"])
    steps["cmd"] = steps["proper_motion"] & tpl.match(stars["c0"], absmag)
    steps["extent"] = steps["cmd"] & (stars["phi1"] >= extent[0]) & (stars["phi1"] <= extent[1])
    return steps


@dataclass
class Midline:
    """The stream's midline phi2(phi1): a smoothing spline through medians in
    sliding windows (PWB2018 used 4-degree windows shifted by 2 degrees), which
    follows the data to within their errors. Held flat beyond its ends."""

    centres: np.ndarray
    median: np.ndarray
    err: np.ndarray

    def __post_init__(self):
        self._spl = UnivariateSpline(self.centres, self.median, w=1 / self.err, k=3,
                                     s=len(self.centres))

    def __call__(self, phi1):
        return self._spl(np.clip(phi1, self.centres[0], self.centres[-1]))


def fit_phi2_track(phi1, phi2, window=4.0, step=2.0, max_offset=2.0, guide=None) -> Midline:
    centres, med, err = [], [], []
    for a in np.arange(GD1_PHI1_RANGE[0], GD1_PHI1_RANGE[1] - window + 1e-9, step):
        c = a + window / 2
        ref = float(guide(np.array([c]))[0]) if guide is not None else 0.0
        s = (phi1 >= a) & (phi1 < a + window) & (np.abs(phi2 - ref) < max_offset)
        if s.sum() >= 15:
            # Median about a local reference; iterate once so the window
            # centres on the stream rather than on phi2 = 0.
            m = np.median(phi2[s])
            s = (phi1 >= a) & (phi1 < a + window) & (np.abs(phi2 - m) < max_offset)
            centres.append(c)
            med.append(np.median(phi2[s]))
            err.append(max(1.2533 * np.std(phi2[s]) / np.sqrt(s.sum()), 0.02))
    return Midline(*map(np.array, (centres, med, err)))


def density_profile(phi1, phi2, phi2_track, ranges=GD1_PHI1_RANGE):
    """Background-subtracted number of member stars per degree along phi1."""
    d = phi2 - phi2_track(phi1)
    on = np.abs(d) < REFERENCE["pwb2018"]["stream_half_width_deg"]
    off = (np.abs(d) > BACKGROUND_BAND[0]) & (np.abs(d) < BACKGROUND_BAND[1])
    centres = np.arange(ranges[0] + DENSITY_WINDOW / 2, ranges[1] - DENSITY_WINDOW / 2 + 1e-9,
                        DENSITY_STEP)
    dens, err, bgd = [], [], []
    for c in centres:
        s = np.abs(phi1 - c) < DENSITY_WINDOW / 2
        t2 = float(phi2_track(np.array([c]))[0])
        alpha = _band_area_ratio_general(t2, REFERENCE["pwb2018"]["stream_half_width_deg"],
                                         BACKGROUND_BAND)
        n_on, n_off = np.sum(s & on), np.sum(s & off)
        dens.append((n_on - alpha * n_off) / DENSITY_WINDOW)
        err.append(np.sqrt(n_on + alpha**2 * n_off) / DENSITY_WINDOW)
        bgd.append(alpha * n_off / DENSITY_WINDOW)
    return centres, np.array(dens), np.array(err), np.array(bgd)


def find_gaps(centres, dens, err, reach_deg=10.0, min_sigma=2.0, min_depth=0.25):
    """Gaps: local minima of the density profile that are significantly lower
    than the highest density on *both* sides within ``reach_deg`` (their
    "prominence", the standard way to measure how deep a valley is).

    depth = 1 - n_gap / n_shoulder, where n_shoulder is the lower of the two
    shoulders; sigma = (n_shoulder - n_gap) / combined Poisson error. Peaks are
    found the same way, upside down. Neighbouring windows overlap, so features
    closer than 3 degrees are merged, keeping the most significant."""
    step = centres[1] - centres[0]
    w = int(round(reach_deg / step))

    def extrema(sign):
        out = []
        y = sign * dens
        for i in range(1, len(y) - 1):
            if not (y[i] <= y[i - 1] and y[i] <= y[i + 1]):
                continue
            left, right = y[max(0, i - w):i], y[i + 1:i + 1 + w]
            if len(left) == 0 or len(right) == 0:
                continue
            il, ir = max(0, i - w) + int(np.argmax(left)), i + 1 + int(np.argmax(right))
            j = il if y[il] < y[ir] else ir  # the lower shoulder
            shoulder = sign * y[j]
            sigma = abs(shoulder - dens[i]) / np.hypot(err[i], err[j])
            depth = 1 - dens[i] / shoulder if sign > 0 else 1 - shoulder / dens[i]
            if sigma > min_sigma and depth > min_depth:
                out.append(dict(phi1=float(centres[i]), density=float(dens[i]),
                                shoulder=float(shoulder), shoulder_phi1=float(centres[j]),
                                depth_fraction=float(depth), sigma=float(sigma),
                                significant=bool(sigma > 3.0)))
        merged = []
        for e in sorted(out, key=lambda e: -e["sigma"]):
            if all(abs(e["phi1"] - m["phi1"]) > 3.0 for m in merged):
                merged.append(e)
        return sorted(merged, key=lambda e: e["phi1"])

    return extrema(+1), extrema(-1)


def off_track_test(phi1, phi2, phi2_track, phi1_range, offset_range):
    """PWB2018-style test of an off-track feature: stars in a box beside the
    stream versus the mirror-image box on the other side (same area)."""
    s = (phi1 > phi1_range[0]) & (phi1 < phi1_range[1])
    d = phi2[s] - phi2_track(phi1[s])
    lo, hi = offset_range
    n_feature = int(np.sum((d > lo) & (d < hi)))
    n_mirror = int(np.sum((d > -hi) & (d < -lo)))
    # Significance of the excess, treating the mirror box as the background.
    from scipy.stats import poisson

    p = float(poisson.sf(n_feature - 1, max(n_mirror, 0.5)))
    return dict(phi1=phi1_range, offset=offset_range, n_feature=n_feature, n_mirror=n_mirror,
                p_value=p)


# --------------------------------------------------------------------------
# Checks against the literature
# --------------------------------------------------------------------------


def load_huang2019():
    rows = [ln.split("\t") for ln in HUANG_PATH.read_text().splitlines()
            if ln and not ln.startswith("#")]
    header, data = rows[0], rows[3:]
    idx = {h.strip(): i for i, h in enumerate(header)}
    def col(name, cast=float):
        return np.array([cast(r[idx[name]].strip()) for r in data])
    return dict(ra=col("_RA"), dec=col("_DE"), g=col("Gmag"), rv=col("vLOS"),
                rv_err=col("e_vLOS"), feh=col("[Fe/H]"))


def crossmatch(ra1, dec1, ra2, dec2, radius_arcsec=1.5):
    """Nearest neighbour of each (ra1, dec1) among (ra2, dec2)."""
    from scipy.spatial import cKDTree

    v2 = gx.unit_vector(ra2, dec2).T
    v1 = gx.unit_vector(ra1, dec1).T
    tree = cKDTree(v2)
    dist, idx = tree.query(v1)
    sep = np.degrees(2 * np.arcsin(dist / 2)) * 3600
    return idx, sep < radius_arcsec, sep


def zero_point_corrected_parallax(stars: Stars) -> np.ndarray:
    """Parallax minus the Lindegren et al. (2021) zero-point (NaN where the
    recipe does not apply)."""
    import warnings

    from zero_point import zpt

    from galaxy import rotate_sky

    zpt.load_tables()
    # Ecliptic latitude from the mean obliquity at J2000 (23.4393 deg).
    eps = np.radians(23.4392911)
    rot = np.array([[1, 0, 0], [0, np.cos(eps), np.sin(eps)], [0, -np.sin(eps), np.cos(eps)]])
    _, ecl_lat = rotate_sky(rot, stars["ra"], stars["dec"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        zp = zpt.get_zpt(stars["phot_g_mean_mag"], stars["nu_eff_used_in_astrometry"],
                         stars["pseudocolour"], ecl_lat, stars["astrometric_params_solved"],
                         _warnings=False)
    return stars["parallax"] - zp


def analyse(g_max: float, n_boot: int = 200) -> tuple[dict, dict]:
    """The whole GD-1 pipeline. Returns (results for results.json, artefacts
    for the figures)."""
    tpl = build_template(g_max)
    stars = load_field()

    # Motion: first with bands centred on phi2 = 0, then again with bands that
    # follow the stream's measured midline.
    track0 = find_pm_track(stars, tpl)
    pre = (stars["good"] & far_enough(stars) & (pm_chi2(stars, track0) < PM_CHI2_MAX)
           & _cmd_prefilter(stars, tpl))
    mid0 = fit_phi2_track(stars["phi1"][pre], stars["phi2"][pre])
    track = find_pm_track(stars, tpl, mid0)

    # Where the stream is detected in proper motion (second pass).
    det = [k for k in track.knots if k["pass"] == 2 and k["detected"]]
    seg = GD1_SEARCH_SEGMENT_DEG
    extent = (min(k["phi1"] for k in det) - seg / 2, max(k["phi1"] for k in det) + seg / 2)

    # Distance: segments, then a smooth model over the detected extent.
    dist = fit_distance(stars, tpl, track, mid0, n_boot=n_boot)
    use = [d for d in dist if extent[0] <= d["phi1"] <= extent[1] and d["dm_err"] < 0.3]
    p = np.array([d["phi1"] for d in use])
    dm = np.array([d["dm"] for d in use])
    # Bootstrap errors are statistical only; M92's own distance (0.018 mag)
    # and the template's colour calibration (~0.05 mag) add a floor.
    dm_err = np.hypot([d["dm_err"] for d in use], 0.03)
    dm_coef, dm_cov = np.polyfit(p, dm, 2, w=1 / dm_err, cov="unscaled")
    lin_coef, lin_cov = np.polyfit(p, dm, 1, w=1 / dm_err, cov="unscaled")
    fit_lo, fit_hi = p.min(), p.max()

    def dm_of_phi1(x):
        return np.polyval(dm_coef, np.clip(x, fit_lo, fit_hi))

    steps = select_members(stars, tpl, track, dm_of_phi1, extent)
    members = steps["extent"]
    phi2_track = fit_phi2_track(stars["phi1"][members], stars["phi2"][members], guide=mid0)

    m_phi1, m_phi2 = stars["phi1"][members], stars["phi2"][members]
    dens_c, dens, dens_err, dens_bg = density_profile(m_phi1, m_phi2, phi2_track)
    in_ext = (dens_c >= extent[0] + DENSITY_WINDOW / 2) & (dens_c <= extent[1] - DENSITY_WINDOW / 2)
    gaps, peaks = find_gaps(dens_c[in_ext], dens[in_ext], dens_err[in_ext])

    ref = REFERENCE["pwb2018"]
    spur = off_track_test(m_phi1, m_phi2, phi2_track, ref["spur_phi1"], ref["spur_phi2"])
    spur_control = off_track_test(m_phi1, m_phi2, phi2_track, ref["control_phi1"],
                                  ref["spur_phi2"])
    blob = off_track_test(m_phi1, m_phi2, phi2_track, ref["blob_phi1"],
                          (-ref["spur_phi2"][1], -ref["spur_phi2"][0]))
    blob_control = off_track_test(m_phi1, m_phi2, phi2_track, ref["control_phi1"],
                                  (-ref["spur_phi2"][1], -ref["spur_phi2"][0]))

    # Stars within the stream band: the "high-confidence" sample.
    on_band = members & (np.abs(stars["phi2"] - phi2_track(stars["phi1"]))
                         < ref["stream_half_width_deg"])
    in_extent = (stars["phi1"] >= extent[0]) & (stars["phi1"] <= extent[1])
    n_band = int(np.sum(on_band & in_extent))
    # Expected field contamination in the band, from the control bands.
    off_band = members & in_extent & (
        (np.abs(stars["phi2"] - phi2_track(stars["phi1"])) > BACKGROUND_BAND[0])
        & (np.abs(stars["phi2"] - phi2_track(stars["phi1"])) < BACKGROUND_BAND[1]))
    alpha_band = _band_area_ratio_general(0.0, ref["stream_half_width_deg"], BACKGROUND_BAND)
    n_contam = float(alpha_band * np.sum(off_band))

    # ---- Checks against independent data -------------------------------
    huang = load_huang2019()
    idx, ok, sep = crossmatch(huang["ra"], huang["dec"], stars["ra"], stars["dec"])
    h_phi1, h_phi2 = gx.gd1_coordinates(huang["ra"], huang["dec"])
    in_field = ((h_phi1 >= GD1_PHI1_RANGE[0]) & (h_phi1 <= GD1_PHI1_RANGE[1])
                & (h_phi2 >= GD1_PHI2_RANGE[0]) & (h_phi2 <= GD1_PHI2_RANGE[1]))
    fate = []
    for i in range(len(huang["ra"])):
        if not in_field[i]:
            fate.append("outside our field")
        elif not ok[i]:
            fate.append("not in our download (cuts at the archive)")
        else:
            j = idx[i]
            if not steps["quality"][j]:
                fate.append("failed quality (RUWE or colour excess)")
            elif not steps["parallax"][j]:
                fate.append("parallax says nearby")
            elif not steps["proper_motion"][j]:
                fate.append("proper motion off the track")
            elif not steps["cmd"][j]:
                fate.append("colour-magnitude off the template")
            elif not steps["extent"][j]:
                fate.append("beyond the detected stretch")
            else:
                fate.append("member")
    fate = np.array(fate)

    # Mean parallax of band members versus the photometric distance. The band
    # holds ~12% field stars that passed every cut; being nearer, they pull the
    # mean parallax up. So the control bands' stars are subtracted, with the
    # same weights, exactly as for the star counts.
    #
    # Errors come from bootstrapping the stars, not from their individual
    # parallax errors: one bright field star at 0.5 mas carries a large weight
    # and moves a segment's mean far more than the formal error admits.
    # Segments where the field stars' weights rival the stream's (at the far
    # end, where the stream's motion approaches the field's) are skipped: the
    # subtraction is unstable there.
    plx = zero_point_corrected_parallax(stars)
    rng = np.random.default_rng(3)

    def bg_mean(p_on, w_on, p_off, w_off):
        sw = np.sum(w_on) - alpha_band * np.sum(w_off)
        return (np.sum(w_on * p_on) - alpha_band * np.sum(w_off * p_off)) / sw, sw

    par_rows, par_skipped = [], []
    for a in np.arange(extent[0], extent[1] - 9.99, 10.0):
        sl = (stars["phi1"] >= a) & (stars["phi1"] < a + 10) & np.isfinite(plx)
        s_on, s_off = on_band & sl, off_band & sl
        if s_on.sum() < 20:
            continue
        p_on, p_off = plx[s_on], plx[s_off]
        w_on = 1 / stars["parallax_error"][s_on] ** 2
        w_off = 1 / stars["parallax_error"][s_off] ** 2
        mean, sw = bg_mean(p_on, w_on, p_off, w_off)
        c = a + 5
        if sw < 0.5 * np.sum(w_on):
            par_skipped.append(dict(phi1=c, n=int(s_on.sum()), n_control=int(s_off.sum()),
                                    stream_weight_fraction=float(sw / np.sum(w_on))))
            continue
        boots = []
        for _ in range(300):
            i = rng.integers(0, len(p_on), len(p_on))
            j = rng.integers(0, len(p_off), len(p_off))
            mb, swb = bg_mean(p_on[i], w_on[i], p_off[j], w_off[j])
            if swb > 0:
                boots.append(mb)
        boots = np.array(boots)
        par_rows.append(dict(
            phi1=c, n=int(s_on.sum()), n_control=int(s_off.sum()), parallax=float(mean),
            parallax_err=float(0.5 * (np.percentile(boots, 84) - np.percentile(boots, 16))),
            parallax_err_formal=float(np.sqrt(np.sum(w_on) + alpha_band**2 * np.sum(w_off)) / sw),
            parallax_raw=float(np.sum(w_on * p_on) / np.sum(w_on)),
            stream_weight_fraction=float(sw / np.sum(w_on)),
            photometric_parallax=float(1 / ph.distance_from_modulus(dm_of_phi1(c)))))

    # Overall agreement of parallax and photometric distance. Gaia DR3
    # parallaxes carry systematics of ~0.01 mas correlated over degree scales
    # (Lindegren et al. 2021; Maiz Apellaniz et al. 2021), which no number of
    # stars averages away; the chi^2 is shown with and without that floor.
    diff = np.array([r["parallax"] - r["photometric_parallax"] for r in par_rows])
    e_stat = np.array([r["parallax_err"] for r in par_rows])
    w = 1 / e_stat**2
    e_sys = np.hypot(e_stat, PARALLAX_SYSTEMATIC_MAS)
    parallax_summary = dict(
        mean_offset_mas=float(np.sum(w * diff) / np.sum(w)),
        mean_offset_err_mas=float(1 / np.sqrt(np.sum(w))),
        chi2_stat=float(np.sum((diff / e_stat) ** 2)),
        chi2_with_systematic=float(np.sum((diff / e_sys) ** 2)),
        dof=int(len(diff)),
        mean_parallax_distance_ratio=float(np.sum(w * np.array([r["parallax"] / r["photometric_parallax"]
                                                               for r in par_rows])) / np.sum(w)),
    )

    # Gaia's own radial velocities for the brightest members, against the
    # published spectroscopic velocities (Koposov 2010 + Huang 2019), which
    # trace a smooth curve along the stream.
    has_rv = on_band & np.isfinite(stars["radial_velocity"])
    kop = REFERENCE["koposov2010"]
    lit_p = np.concatenate([kop["rv_phi1"], h_phi1])
    lit_v = np.concatenate([kop["rv_kms"], huang["rv"]])
    lit_e = np.hypot(np.concatenate([kop["rv_err_kms"], huang["rv_err"]]), 10.0)
    rv_coef = np.polyfit(lit_p, lit_v, 3, w=1 / lit_e)
    rv_resid = stars["radial_velocity"][has_rv] - np.polyval(rv_coef, stars["phi1"][has_rv])
    rv_tol = 3 * np.hypot(stars["radial_velocity_error"][has_rv], 15.0)
    rv_check = dict(n=int(has_rv.sum()), n_consistent=int(np.sum(np.abs(rv_resid) < rv_tol)),
                    residuals_kms=[float(x) for x in rv_resid],
                    g_mag=[float(x) for x in stars["phot_g_mean_mag"][has_rv]],
                    literature_track_coef=rv_coef.tolist())

    # The stream's proper motion against Price-Whelan & Bonaca's measurement:
    # their mean for phi1 in (-18, -10) deg at (RA, Dec) = (177.01, 53.99).
    pwb = REFERENCE["pwb2018"]
    ra_p, dec_p = pwb["radec_at_m14"]
    p1_p, _, m1_p, m2_p, _ = gx.gd1_coordinates(np.array([ra_p]), np.array([dec_p]),
                                                np.array([pwb["pm_radec_at_m14"][0]]),
                                                np.array([pwb["pm_radec_at_m14"][1]]))
    t1_p, t2_p = track.at(p1_p)
    pm_check = dict(phi1=float(p1_p[0]), published=[float(m1_p[0]), float(m2_p[0])],
                    ours=[float(t1_p[0]), float(t2_p[0])],
                    difference=float(np.hypot(t1_p[0] - m1_p[0], t2_p[0] - m2_p[0])))

    results = dict(
        pm_check=pm_check,
        template=dict(cluster=M92.name, n_members=tpl.n_members,
                      pm_knot=[float(x) for x in tpl.pm_knot], distance_modulus=tpl.dm,
                      distance_kpc=M92.distance_kpc, turnoff_absmag=float(
                          tpl.absmag[np.argmin(tpl.colour)]),
                      turnoff_colour=float(np.min(tpl.colour))),
        field=dict(n_stars=len(stars), n_good=int(stars["good"].sum())),
        selection={k: int(v.sum()) for k, v in steps.items()},
        n_band=n_band, n_band_contamination=n_contam,
        extent_phi1=[float(x) for x in extent],
        pm_knots=[{k: v for k, v in kk.items()} for kk in track.knots],
        pm_track=dict(coef_pm1=track.coef1.tolist(), coef_pm2=track.coef2.tolist()),
        distance=[{k: v for k, v in d.items() if k != "signal"} for d in dist],
        distance_model=dict(quadratic=dm_coef.tolist(), linear=lin_coef.tolist(),
                            linear_err=np.sqrt(np.diag(lin_cov)).tolist(),
                            phi1_range=[float(fit_lo), float(fit_hi)]),
        phi2_track=dict(centres=phi2_track.centres.tolist(), median=phi2_track.median.tolist(),
                        err=phi2_track.err.tolist()),
        density=dict(phi1=dens_c.tolist(), per_degree=dens.tolist(), err=dens_err.tolist(),
                     background=dens_bg.tolist()),
        gaps=gaps, peaks=peaks,
        spur=spur, spur_control=spur_control, blob=blob, blob_control=blob_control,
        huang2019=dict(n=int(len(fate)), fate={f: int(np.sum(fate == f)) for f in np.unique(fate)},
                       median_feh=float(np.median(huang["feh"]))),
        parallax_check=par_rows, parallax_skipped=par_skipped,
        parallax_summary=parallax_summary,
        n_with_gaia_rv=int(has_rv.sum()), rv_check=rv_check,
    )
    artefacts = dict(tpl=tpl, stars=stars, steps=steps, track=track, track0=track0,
                     dist=dist, dm_of_phi1=dm_of_phi1, phi2_track=phi2_track, members=members,
                     on_band=on_band, dens_in_ext=in_ext, huang=huang, rv_coef=rv_coef,
                     extent=extent,
                     huang_fate=fate, huang_phi1=h_phi1, huang_phi2=h_phi2, has_rv=has_rv,
                     parallax_corrected=plx)
    return results, artefacts
