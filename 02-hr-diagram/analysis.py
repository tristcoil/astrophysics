"""The analysis: from raw Gaia rows to a clean HR diagram, and what we read off it.

The pipeline, in order:

1. **Derive** absolute magnitude, tangential velocity and the photometric
   quality metric C* for every candidate.
2. **Cut**, one criterion at a time, recording for every rejected star the
   *first* cut it failed. That attribution is what lets the figures show what
   each cut removes and where on the sky and the diagram those stars were.
3. **Map the diagram**: the main-sequence ridge line, and rough regions (main
   sequence, giants, white dwarfs) for counting.
4. **Measure** features with known published values -- the binary sequence, the
   red clump, the Jao gap, the white dwarf branches -- so the diagram is tested,
   not just drawn.
5. **Kinematics**: the same diagram split by how fast stars move.
6. **Sky map** of the whole catalogue's star counts.

Nothing here plots; ``figures.py`` draws from what this module returns.
"""

from __future__ import annotations

import json
from dataclasses import dataclass

import numpy as np
from scipy import ndimage

import hr_physics as hp
from config import (
    EXCESS_G_MIN,
    EXCESS_NSIGMA,
    MAX_BP_MAG,
    MAX_RUWE,
    MIN_PARALLAX_OVER_ERROR,
    PARALLAX_MIN_MAS,
    RESULTS_PATH,
    SUN_ABS_G,
    SUN_BP_RP,
)


# --------------------------------------------------------------------------
# 1. Derived quantities
# --------------------------------------------------------------------------

FLOAT_COLUMNS = [
    "ra", "dec", "l", "b", "parallax", "parallax_error", "parallax_over_error",
    "pmra", "pmdec", "ruwe", "phot_g_mean_mag", "phot_bp_mean_mag",
    "phot_rp_mean_mag", "bp_rp", "phot_g_mean_flux_over_error",
    "phot_bp_mean_flux_over_error", "phot_rp_mean_flux_over_error",
    "phot_bp_rp_excess_factor", "radial_velocity",
]


def derive(table) -> dict[str, np.ndarray]:
    """Plain arrays for every column the analysis uses, plus derived ones."""
    d = {name: np.asarray(table[name], dtype=float) for name in FLOAT_COLUMNS}
    d["source_id"] = np.asarray(table["source_id"], dtype=np.int64)
    d["params_solved"] = np.asarray(table["astrometric_params_solved"], dtype=int)
    d["abs_g"] = hp.absolute_magnitude(d["phot_g_mean_mag"], d["parallax"])
    d["distance_pc"] = hp.distance_pc(d["parallax"])
    d["v_tan"] = hp.tangential_velocity(d["pmra"], d["pmdec"], d["parallax"])
    d["c_star"] = hp.flux_excess_corrected(d["bp_rp"], d["phot_bp_rp_excess_factor"])
    d["c_star_sigma"] = hp.flux_excess_sigma(d["phot_g_mean_mag"])
    return d


# --------------------------------------------------------------------------
# 2. Quality cuts
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Cut:
    key: str
    #: Short label naming what the cut *removes*, for figures and tables.
    label: str
    #: One sentence on why the rejected stars are untrustworthy.
    why: str
    passes: np.ndarray


def quality_cuts(d) -> list[Cut]:
    """The cuts, in the order they are applied. Each keeps a star if True."""
    has_colour = np.isfinite(d["bp_rp"]) & np.isfinite(d["phot_g_mean_mag"])
    g = d["phot_g_mean_mag"]
    # C* is compared against its expected scatter at the star's magnitude. The
    # published scatter relation does not hold for saturated stars (G < 4), so
    # those pass this particular cut by construction.
    excess_ok = has_colour & (
        (g < EXCESS_G_MIN)
        | (np.abs(d["c_star"]) < EXCESS_NSIGMA * d["c_star_sigma"])
    )
    return [
        Cut("parallax", f"parallax measured to worse than {100 / MIN_PARALLAX_OVER_ERROR:.0f}%",
            "the distance, and so the absolute magnitude, is too uncertain",
            d["parallax_over_error"] > MIN_PARALLAX_OVER_ERROR),
        Cut("ruwe", f"RUWE ≥ {MAX_RUWE}",
            "the motion on the sky does not fit a single star: usually an "
            "unresolved binary, sometimes a spurious solution",
            d["ruwe"] < MAX_RUWE),
        Cut("colour", "no BP−RP colour",
            "without a colour a star cannot be placed on the diagram",
            has_colour),
        Cut("excess", f"BP+RP flux inconsistent with G (|C*| ≥ {EXCESS_NSIGMA:.0f}σ)",
            "a neighbour's light has leaked into the BP/RP measurement, "
            "corrupting the colour",
            excess_ok),
        Cut("faint_bp", f"G_BP fainter than {MAX_BP_MAG}",
            "BP fluxes this faint are biased bright, making the colour too blue",
            d["phot_bp_mean_mag"] < MAX_BP_MAG),
    ]


def apply_cuts(d, cuts: list[Cut]):
    """Apply cuts in order. Returns (kept mask, first failed cut index, flow).

    ``first_fail`` is -1 for stars that pass everything, otherwise the index of
    the first cut they failed. Attributing each star to one cut makes the flow
    add up exactly, at the price of order-dependence, which the flow records.
    """
    n = len(d["parallax"])
    alive = np.ones(n, dtype=bool)
    first_fail = np.full(n, -1, dtype=int)
    flow = [{"cut": f"parallax > {PARALLAX_MIN_MAS:.0f} mas (the archive query)",
             "removed": 0, "remaining": n, "why": "naive distance under 100 pc"}]
    for i, cut in enumerate(cuts):
        failing = alive & ~cut.passes
        first_fail[failing] = i
        alive &= cut.passes
        flow.append({"cut": cut.label, "removed": int(failing.sum()),
                     "remaining": int(alive.sum()), "why": cut.why})
    return alive, first_fail, flow


# --------------------------------------------------------------------------
# 3. The main-sequence ridge line, and regions of the diagram
# --------------------------------------------------------------------------

#: Fewest stars in a colour window for its ridge point to be trusted.
RIDGE_MIN_STARS = 50


def ridge_line(colour, absmag, *, lo=0.2, hi=4.4, step=0.05, width=0.10,
               mag_bin=0.02, smooth_mag=0.06):
    """The main sequence's ridge: the most common M_G at each colour.

    The *mode*, not the median. Unresolved binaries sit above the main sequence
    and would drag a median upwards; the densest point is the single-star
    sequence itself. Within each colour window the M_G distribution is
    histogrammed finely, lightly smoothed, and its peak taken.

    Returns (colour grid, ridge M_G, number of stars per window).
    """
    grid = np.arange(lo, hi + 1e-9, step)
    ridge = np.full(grid.size, np.nan)
    count = np.zeros(grid.size, dtype=int)
    edges = np.arange(-3.0, 20.0 + mag_bin, mag_bin)
    sigma_bins = smooth_mag / mag_bin
    for i, c in enumerate(grid):
        sel = np.abs(colour - c) < width / 2
        count[i] = int(sel.sum())
        if count[i] < RIDGE_MIN_STARS:
            continue
        h, _ = np.histogram(absmag[sel], bins=edges)
        hs = ndimage.gaussian_filter1d(h.astype(float), sigma_bins)
        ridge[i] = 0.5 * (edges[:-1] + edges[1:])[np.argmax(hs)]
    return grid, ridge, count


def ridge_at(colour, grid, ridge):
    """Interpolate the ridge line to arbitrary colours (NaN outside it)."""
    ok = np.isfinite(ridge)
    if ok.sum() < 2:
        return np.full(np.shape(colour), np.nan)
    return np.interp(colour, grid[ok], ridge[ok], left=np.nan, right=np.nan)


# Regions, for counting only. Boundaries are drawn through empty gaps.

#: White dwarfs lie below the line M_G = WD_LINE[0] + WD_LINE[1] * (BP-RP),
#: drawn through the nearly empty gap between them and the main sequence.
WD_LINE = (9.0, 3.0)
#: Giants and subgiants: redder than this and more than GIANT_DELTA mag above
#: the main-sequence ridge. Bluer than it, evolved stars and the upper main
#: sequence overlap too much to separate by position alone.
GIANT_MIN_COLOUR = 0.75
GIANT_DELTA = -1.4
#: Stars more than this far *below* the ridge but above the white-dwarf line.
BELOW_MS_DELTA = 1.0


def regions(colour, absmag, delta):
    wd = absmag > WD_LINE[0] + WD_LINE[1] * colour
    giant = ~wd & (colour > GIANT_MIN_COLOUR) & (delta < GIANT_DELTA)
    below = ~wd & ~giant & (delta > BELOW_MS_DELTA)
    ms = ~wd & ~giant & ~below
    return {"main_sequence": ms, "white_dwarfs": wd, "giants": giant, "other": below}


# --------------------------------------------------------------------------
# 4. Measurements with published counterparts
# --------------------------------------------------------------------------

def _bootstrap(fn, *arrays, n=200, seed=0):
    """Standard deviation of ``fn`` over bootstrap resamples of the stars."""
    rng = np.random.default_rng(seed)
    size = len(arrays[0])
    values = []
    for _ in range(n):
        idx = rng.integers(0, size, size)
        v = fn(*(a[idx] for a in arrays))
        if v is not None and np.isfinite(v):
            values.append(v)
    return float(np.std(values)) if len(values) > n // 2 else None


def _smoothed_peak(values, lo, hi, *, bin_width=0.02, smooth=0.04):
    """Location of a genuine local maximum of a smoothed histogram in [lo, hi].

    Returns None when the highest point in the window sits on its edge: that
    is the flank of something outside the window, not a peak of its own.
    """
    edges = np.arange(lo - 0.2, hi + 0.2 + bin_width / 2, bin_width)
    h, _ = np.histogram(values, bins=edges)
    hs = ndimage.gaussian_filter1d(h.astype(float), smooth / bin_width)
    centres = 0.5 * (edges[:-1] + edges[1:])
    inside = np.where((centres >= lo) & (centres <= hi))[0]
    k = inside[np.argmax(hs[inside])]
    if k in (inside[0], inside[-1]):
        return None
    return float(centres[k])


def binary_sequence(colour, delta, rejected_delta=None, *, colour_range=(1.0, 1.8),
                    window=(-1.1, -0.5)):
    """Offset of the unresolved-twin sequence above the single-star ridge.

    ``delta`` is M_G minus the ridge at each star's colour. Two near-identical
    stars too close for Gaia to separate appear as one dot at almost the same
    colour, 2.5 log10(2) = 0.753 mag brighter: twice the light. We look for
    that secondary peak on the bright side of the ridge.

    The colour range is where it can be seen. Among K and early-M dwarfs the
    single-star sequence is narrow (about 0.1 mag), so the twins stand clear of
    it; among later M dwarfs it is several times broader -- a spread of
    metallicities -- and the twin peak drowns in its wing.

    ``rejected_delta`` (stars removed by the RUWE cut, at the same colours) is
    carried along for the figure: that cut removes most unresolved binaries,
    so their distribution shows the binary sequence the clean sample lacks.
    """
    sel = (colour > colour_range[0]) & (colour < colour_range[1]) & np.isfinite(delta)
    dd = delta[sel]
    if sel.sum() < 500:
        return {"ok": False, "n_stars": int(sel.sum())}
    offset = _smoothed_peak(dd, *window)
    if offset is None:
        return {"ok": False, "n_stars": int(sel.sum()),
                "summary": "no distinct twin-binary peak found"}
    err = _bootstrap(lambda v: _smoothed_peak(v, *window), dd)
    predicted = float(hp.combined_magnitude_offset(1.0))
    edges = np.arange(-2.0, 1.0 + 0.02, 0.02)
    hist, _ = np.histogram(dd, bins=edges)
    out = {
        "ok": True, "offset": offset, "offset_err": err, "predicted": predicted,
        "colour_range": list(colour_range), "n_stars": int(sel.sum()),
        "hist_edges": edges.tolist(), "hist": hist.tolist(),
        "summary": (f"twin-binary sequence {abs(offset):.2f} ± {err:.2f} mag above the "
                    f"single-star ridge (BP−RP {colour_range[0]}–{colour_range[1]}); two "
                    f"identical stars predict {abs(predicted):.3f}"),
    }
    if rejected_delta is not None:
        rej = rejected_delta[np.isfinite(rejected_delta)]
        out["rejected_hist"] = np.histogram(rej, bins=edges)[0].tolist()
        out["n_rejected"] = int(rej.size)
        # How much of each group sits on the binary side (0.4-1.1 mag above the ridge)?
        bright = lambda v: float(np.mean((v < -0.4) & (v > -1.1)))
        out["binary_side_fraction_kept"] = bright(dd)
        out["binary_side_fraction_rejected"] = bright(rej)
    return out


def _gauss_line(x, a, mu, sigma, c0, c1):
    return a * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + c0 + c1 * (x - 0.5)


def _fit_clump(m):
    """Gaussian + straight-line fit to the M_G histogram of red giants."""
    from scipy.optimize import curve_fit

    edges = np.arange(-1.5, 2.5 + 0.1, 0.1)
    h, _ = np.histogram(m, bins=edges)
    c = 0.5 * (edges[:-1] + edges[1:])
    try:
        p, _ = curve_fit(_gauss_line, c, h, p0=[h.max(), 0.45, 0.2, np.median(h), 0.0],
                         sigma=np.sqrt(np.maximum(h, 1)), maxfev=20000)
    except RuntimeError:
        return None
    return p


def red_clump(colour, absmag, *, box=((0.95, 1.55), (-1.5, 2.5))):
    """The red clump's absolute magnitude, a classic standard candle.

    Red-clump stars burn helium in their cores at a nearly fixed luminosity,
    set by the mass of the helium core at ignition. We fit a Gaussian plus a
    linear background (the red giant branch passing through) to the M_G
    distribution of giants in the clump's colour range.
    """
    (c_lo, c_hi), (m_lo, m_hi) = box
    sel = (colour > c_lo) & (colour < c_hi) & (absmag > m_lo) & (absmag < m_hi)
    m, c = absmag[sel], colour[sel]
    if sel.sum() < 40:   # a fit to a handful of stars would just return its guess
        return {"ok": False, "n_giants": int(sel.sum())}
    p = _fit_clump(m)
    if p is None or not (m_lo < p[1] < m_hi):
        return {"ok": False, "n_giants": int(sel.sum())}
    err = _bootstrap(lambda v: (lambda q: q[1] if q is not None else None)(_fit_clump(v)), m)
    in_clump = np.abs(m - p[1]) < abs(p[2])
    # Two published calibrations, which differ from each other by 0.055 mag:
    # the clump's brightness depends on the age and metallicity mix of the
    # stars, and on how extinction was treated.
    literature = [
        {"ref": "Ruiz-Dern et al. 2018", "abs_g": 0.495, "err": 0.009,
         "note": "2,482 low-extinction TGAS giants, no dust correction (A&A 609, A116)"},
        {"ref": "Hawkins et al. 2017", "abs_g": 0.44, "err": 0.01,
         "note": "972 TGAS giants, corrected with a 3D dust map (MNRAS 471, 722)"},
    ]
    return {
        "ok": True, "abs_g": float(p[1]), "abs_g_err": err, "width": abs(float(p[2])),
        "bp_rp": float(np.median(c[in_clump])), "n_giants": int(sel.sum()),
        "n_in_clump": int(in_clump.sum()), "fit": [float(v) for v in p],
        "literature": literature,
        "summary": (f"red clump at M_G = {p[1]:.2f} ± {err:.2f} (spread {abs(p[2]):.2f}), "
                    f"{int(in_clump.sum())} stars within 1σ; published "
                    + ", ".join(f"{r['abs_g']:+.3f} ± {r['err']:.3f} ({r['ref']})"
                                for r in literature)),
    }


def _masked_smooth(h, sigma_bins, mask):
    """Gaussian smoothing that ignores masked bins (normalised convolution)."""
    w = (~mask).astype(float)
    num = ndimage.gaussian_filter1d(h * w, sigma_bins, mode="nearest")
    den = ndimage.gaussian_filter1d(w, sigma_bins, mode="nearest")
    return num / np.maximum(den, 1e-9)


def find_dip(values, *, bin_width=0.01, baseline=0.15, gap_halfwidth=0.025):
    """The most significant narrow deficit in a histogram, against a smooth baseline.

    The baseline is the histogram smoothed on a scale much wider than any gap.
    Smoothing *across* a gap would drag the baseline down into it and shrink the
    measured deficit, so the search is iterated: find the deepest dip, mask it,
    rebuild the baseline without it, and measure the dip against that. Only the
    well-populated core of the distribution is searched -- its thinning flanks
    would otherwise masquerade as deficits.
    """
    edges = np.arange(np.percentile(values, 1), np.percentile(values, 99) + bin_width,
                      bin_width)
    h = np.histogram(values, bins=edges)[0].astype(float)
    centres = 0.5 * (edges[:-1] + edges[1:])
    sig = baseline / bin_width
    half = max(1, int(round(gap_halfwidth / bin_width)))
    idx = np.arange(len(h))
    mask = np.zeros_like(h, dtype=bool)
    for _ in range(3):
        base = _masked_smooth(h, sig, mask)
        local = ndimage.uniform_filter1d(h, 2 * half + 1, mode="nearest")
        z = (local - base) / np.sqrt(np.maximum(base / (2 * half + 1), 1e-9))
        z[base < 0.5 * base.max()] = np.inf
        k = int(np.argmin(z))
        mask = np.abs(idx - k) <= 2 * half
    base = _masked_smooth(h, sig, mask)
    window = np.abs(idx - k) <= half
    observed, expected = h[window].sum(), base[window].sum()
    return {
        "centre": float(centres[k]),
        "deficit": float(1 - observed / expected),
        # Deficit in Poisson sigmas: positive means missing stars.
        "significance": float((expected - observed) / np.sqrt(expected)),
        "centres": centres.tolist(), "hist": h.tolist(), "baseline": base.tolist(),
    }


#: Jao et al. (2018, ApJL 861, L11) quote the DR2 gap at M_G = 10.09 for
#: BP-RP 2.35-2.40 and 10.24 for 2.55-2.60: a line of slope 0.75 through
#: (2.375, 10.09). We adopt their slope and measure the gap's position.
JAO_SLOPE, JAO_PIVOT, JAO_DR2_MAG = 0.75, 2.375, 10.09
JAO_RANGE = (2.30, 2.60)
#: Control colour ranges with no known gap: the same search there shows how
#: deep a dip noise alone produces.
JAO_CONTROLS = [(1.90, 2.20), (2.00, 2.30), (2.65, 2.95), (2.80, 3.10), (3.00, 3.30)]


def jao_dr2_gap(colour):
    """The gap's M_G in Gaia DR2 as a function of colour (Jao et al. 2018)."""
    return JAO_DR2_MAG + JAO_SLOPE * (np.asarray(colour) - JAO_PIVOT)


def jao_gap(colour, absmag, delta, *, band=0.4):
    """Measure the Jao gap: a thin deficit of stars across the M dwarfs.

    It marks the mass (~0.35 M_sun) at which red dwarfs become fully
    convective. The gap is only ~0.05 mag wide and slightly tilted, so we
    straighten the diagram along it -- m' = M_G - 0.75 (BP-RP - 2.375), constant
    along the published DR2 gap -- and search the m' distribution of
    main-sequence stars (within ``band`` of the ridge) for its deepest dip.
    Using the published slope makes this a confirmation test of a known feature,
    with only the gap's position left free; the control ranges calibrate how
    significant a dip noise alone can produce.
    """
    def aligned(lo, hi, pivot):
        sel = (colour >= lo) & (colour < hi) & (np.abs(delta) < band)
        return absmag[sel] - JAO_SLOPE * (colour[sel] - pivot), int(sel.sum())

    mp, n = aligned(*JAO_RANGE, JAO_PIVOT)
    if n < 2000:
        return {"ok": False, "n_stars": n}
    dip = find_dip(mp)
    err = _bootstrap(lambda v: find_dip(v)["centre"], mp, n=100)
    controls = []
    for lo, hi in JAO_CONTROLS:
        v, nc = aligned(lo, hi, lo + 0.075)
        cd = find_dip(v)
        controls.append({"colour_lo": lo, "colour_hi": hi, "n_stars": nc,
                         "deficit": cd["deficit"], "significance": cd["significance"]})
    worst = max(c["significance"] for c in controls)
    ok = dip["significance"] > worst + 1.0
    dip.update({
        "ok": ok, "n_stars": n, "centre_err": err, "colour_range": list(JAO_RANGE),
        "slope": JAO_SLOPE, "pivot": JAO_PIVOT, "dr2_mag": JAO_DR2_MAG,
        "controls": controls, "max_control_significance": worst,
        "literature": "Jao et al. 2018 (ApJL 861, L11), DR2: M_G = 10.09 at BP−RP "
                      "2.35–2.40, 10.24 at 2.55–2.60; ~0.05 mag wide, up to 28% deficit",
        "summary": (f"Jao gap at M_G = {dip['centre']:.2f} ± {err:.2f} (at BP−RP "
                    f"{JAO_PIVOT}; DR2: {JAO_DR2_MAG}), {dip['deficit']:.0%} fewer stars, "
                    f"{dip['significance']:.1f}σ; the same search in control ranges "
                    f"finds at most {worst:.1f}σ"),
    })
    return dip


def white_dwarf_branches(colour, absmag, *, lo=12.0, hi=14.0, step=0.2):
    """The split of the white dwarf cooling sequence into two branches.

    Gaia Collaboration, Babusiaux et al. (2018) found the white dwarf sequence
    divides into a main branch A, on the 0.6 M_sun hydrogen-atmosphere cooling
    track, and a second branch B just below it (bluer at a given M_G), which
    their colour-colour analysis links to helium atmospheres. In each narrow
    slice of M_G we histogram the colours and look for two peaks.
    """
    edges = np.arange(-0.4, 0.9 + 0.01, 0.02)
    centres = 0.5 * (edges[:-1] + edges[1:])
    rows = []
    for m0 in np.arange(lo, hi - 1e-9, step):
        sel = (absmag >= m0) & (absmag < m0 + step)
        hs = ndimage.gaussian_filter1d(np.histogram(colour[sel], bins=edges)[0].astype(float), 1.0)
        peaks = [i for i in range(1, len(hs) - 1)
                 if hs[i] > hs[i - 1] and hs[i] >= hs[i + 1] and hs[i] > 0.25 * hs.max()]
        peaks = sorted(sorted(peaks, key=lambda i: hs[i])[-2:])
        row = {"abs_g_lo": float(m0), "abs_g_hi": float(m0 + step), "n": int(sel.sum())}
        if len(peaks) == 2:
            b, a = peaks   # bluer peak = branch B, redder = branch A
            row.update({"b_colour": float(centres[b]), "a_colour": float(centres[a]),
                        "separation": float(centres[a] - centres[b]),
                        "a_over_b": float(hs[a] / hs[b])})
        rows.append(row)
    split = [r for r in rows if "separation" in r]
    ok = len(split) >= 0.8 * len(rows)
    sep = float(np.median([r["separation"] for r in split])) if split else None
    ratio = float(np.median([r["a_over_b"] for r in split])) if split else None
    return {
        "ok": ok, "slices": rows, "n_split": len(split), "n_slices": len(rows),
        "median_separation": sep, "median_a_over_b": ratio,
        "literature": "Gaia Collaboration, Babusiaux et al. 2018 (A&A 616, A10), Sect. 4.4",
        "summary": (f"white dwarf sequence splits in two in {len(split)} of {len(rows)} "
                    f"slices of M_G {lo:g}–{hi:g}: branches {sep:.2f} mag apart in "
                    f"BP−RP, the redder (A) {ratio:.1f}× as populous at its peak"
                    if split else "no white dwarf bifurcation found"),
    }


# --------------------------------------------------------------------------
# 5. Kinematics
# --------------------------------------------------------------------------

def tangential_velocity_vectors(d, sel):
    """Heliocentric tangential velocities as 3-vectors in Galactic coordinates.

    Returns (p, rhat): p[i] is star i's velocity across the sky in km/s, as a
    vector (x towards the Galactic centre, y in the direction of rotation, z to
    the north Galactic pole); rhat[i] is the unit vector towards the star.
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    c = SkyCoord(ra=d["ra"][sel] * u.deg, dec=d["dec"][sel] * u.deg,
                 pm_ra_cosdec=d["pmra"][sel] * u.mas / u.yr,
                 pm_dec=d["pmdec"][sel] * u.mas / u.yr, frame="icrs").galactic
    l = c.l.radian
    b = c.b.radian
    k = hp.KM_S_PER_AU_PER_YR / d["parallax"][sel]
    v_l = k * c.pm_l_cosb.to_value(u.mas / u.yr)
    v_b = k * c.pm_b.to_value(u.mas / u.yr)
    e_l = np.stack([-np.sin(l), np.cos(l), np.zeros_like(l)], axis=1)
    e_b = np.stack([-np.cos(l) * np.sin(b), -np.sin(l) * np.sin(b), np.cos(b)], axis=1)
    rhat = np.stack([np.cos(l) * np.cos(b), np.sin(l) * np.cos(b), np.sin(b)], axis=1)
    return v_l[:, None] * e_l + v_b[:, None] * e_b, rhat


def proper_motion_dispersion(p, rhat, *, clip=4.0, iterations=5):
    """Mean space velocity and velocity dispersion S from proper motions alone.

    The method of Dehnen & Binney (1998, MNRAS 298, 387), used by Binney,
    Dehnen & Bertelli (2000) to measure Parenago's discontinuity. Each star's
    tangential velocity p is the projection of its space velocity onto the sky,
    p = A v with A = I - rhat rhat^T. For stars spread over the whole sky the
    mean velocity follows by least squares, <v> = (sum A)^-1 sum p, and

        S^2 = < |p - A <v>|^2 >

    measures the random motions left over. The missing radial component is why
    this works only statistically, over many directions -- for one star it
    cannot.

    Stars more than ``clip`` S away are dropped and the fit repeated: a handful
    of halo stars moving at 300 km/s would otherwise dominate S for the whole
    bin.
    """
    keep = np.ones(len(p), dtype=bool)
    for _ in range(iterations):
        A = np.eye(3)[None, :, :] - rhat[keep, :, None] * rhat[keep, None, :]
        v_mean = np.linalg.solve(A.sum(axis=0), p[keep].sum(axis=0))
        proj = v_mean - (rhat @ v_mean)[:, None] * rhat   # A_i <v> for every star
        resid = np.linalg.norm(p - proj, axis=1)
        s = float(np.sqrt(np.mean(resid[keep] ** 2)))
        new_keep = resid < clip * s
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep
    return {"v_mean": v_mean.tolist(), "S": s, "n": int(keep.sum()),
            "n_clipped": int((~keep).sum())}


def kinematic_classes(v_tan):
    from config import VTAN_HALO_MIN, VTAN_THICK_DISK, VTAN_THIN_DISK_MAX

    return {
        f"thin-disk-like (v_tan < {VTAN_THIN_DISK_MAX:.0f} km/s)": v_tan < VTAN_THIN_DISK_MAX,
        f"thick-disk-like ({VTAN_THICK_DISK[0]:.0f}–{VTAN_THICK_DISK[1]:.0f} km/s)":
            (v_tan > VTAN_THICK_DISK[0]) & (v_tan < VTAN_THICK_DISK[1]),
        f"halo-like (v_tan > {VTAN_HALO_MIN:.0f} km/s)": v_tan > VTAN_HALO_MIN,
    }


def parenago_curve(d, sel, colour, *, lo=-0.1, hi=2.6, step=0.1, min_stars=150):
    """Velocity dispersion S of main-sequence stars against colour.

    Bluer main-sequence stars are more massive and burn out sooner, so in any
    colour bin bluer than the colour of the oldest surviving main-sequence
    stars, everyone is young. Young stars have had little time to be stirred up
    by spiral arms and molecular clouds, so they move slowly. Redward of that
    colour the main sequence holds stars of every age and the dispersion stops
    rising: Parenago's discontinuity.
    """
    p, rhat = tangential_velocity_vectors(d, sel)
    c = colour[sel]
    rows = []
    for c0 in np.arange(lo, hi, step):
        in_bin = (c >= c0) & (c < c0 + step)
        if in_bin.sum() < min_stars:
            continue
        fit = proper_motion_dispersion(p[in_bin], rhat[in_bin])
        err = fit["S"] / np.sqrt(2 * fit["n"])   # large-N error on an rms
        rows.append({"colour": float(np.median(c[in_bin])), "S": fit["S"], "S_err": err,
                     "n": fit["n"], "n_clipped": fit["n_clipped"]})
    return rows


def _break_colour(rows):
    """Colour where the S(colour) curve changes slope: best continuous two-line fit.

    Returns (break colour, its 1-sigma range from the chi-square profile, the
    fit coefficients, chi-square, degrees of freedom). Two straight lines are
    only an approximation to the real curve, so chi2/dof exceeds 1; the errors
    are rescaled until it equals 1 (delta-chi2 = chi2_min/dof instead of 1),
    which widens the range honestly instead of quoting a spuriously tight one.
    """
    x = np.array([r["colour"] for r in rows])
    y = np.array([r["S"] for r in rows])
    w = 1 / np.array([r["S_err"] for r in rows]) ** 2
    trials = np.arange(x[2], x[-3], 0.005)
    chi2, coefs = [], []
    for xb in trials:
        # y = a + b (x - xb) below the break, a + c (x - xb) above: continuous at xb.
        X = np.stack([np.ones_like(x), np.minimum(x - xb, 0), np.maximum(x - xb, 0)], axis=1)
        coef, *_ = np.linalg.lstsq(X * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)
        chi2.append(float(np.sum(w * (y - X @ coef) ** 2)))
        coefs.append(coef)
    chi2 = np.array(chi2)
    best = int(np.argmin(chi2))
    dof = len(x) - 4
    within = trials[chi2 <= chi2[best] + max(1.0, chi2[best] / dof)]
    return float(trials[best]), (float(within.min()), float(within.max())), coefs[best], \
        float(chi2[best]), dof


# --------------------------------------------------------------------------
# 6. The all-sky star-count map
# --------------------------------------------------------------------------

def sky_density(counts, meta, *, level=None, width=1600):
    """Stars per square degree on a regular Galactic-coordinate grid.

    ``counts`` holds, for each HEALPix pixel (ICRS, nested ordering -- the
    scheme built into Gaia's source_id), how many sources of the random
    subsample fell in it. Scaling by (catalogue size / subsample size) and
    dividing by the pixel area gives the full catalogue's density.

    ``level`` optionally coarsens the map: in the nested scheme the parent of
    pixel i one level up is i // 4, so summing to a lower level is exact. That
    trades resolution for less sampling noise where the sky is sparse.

    The grid's x coordinate is -l, so Galactic longitude increases to the
    *left*, the convention for maps of the sky seen from inside.
    """
    import astropy.units as u
    from astropy.coordinates import ICRS, SkyCoord
    from astropy_healpix import HEALPix

    fetched_level = int(meta["LEVEL"])
    level = fetched_level if level is None else int(level)
    if level > fetched_level:
        raise ValueError(f"cannot refine level {fetched_level} counts to level {level}")
    n_total, n_sample = int(meta["NTOTAL"]), int(meta["NSAMPLE"])
    hpx = HEALPix(nside=2 ** level, order="nested", frame=ICRS())
    per_pixel = np.zeros(hpx.npix)
    parent = np.asarray(counts["hpx"], dtype=np.int64) // 4 ** (fetched_level - level)
    np.add.at(per_pixel, parent, np.asarray(counts["n"], dtype=float))
    area_deg2 = hpx.pixel_area.to_value(u.deg ** 2)
    scale = n_total / n_sample
    density = per_pixel * scale / area_deg2

    # Where on the sky the stars are: Galactic latitude of every pixel centre.
    centres = hpx.healpix_to_skycoord(np.arange(hpx.npix)).galactic
    abs_b = np.abs(centres.b.deg)
    order = np.argsort(abs_b)
    cumulative = np.cumsum(per_pixel[order]) / per_pixel.sum()
    half_within_b = float(abs_b[order][np.searchsorted(cumulative, 0.5)])
    plane_area_fraction = float(np.mean(abs_b < half_within_b))

    x_edges = np.linspace(-180.0, 180.0, width + 1)
    y_edges = np.linspace(-90.0, 90.0, width // 2 + 1)
    xc = 0.5 * (x_edges[:-1] + x_edges[1:])
    yc = 0.5 * (y_edges[:-1] + y_edges[1:])
    X, Y = np.meshgrid(xc, yc)
    grid_coords = SkyCoord(l=np.mod(-X, 360.0) * u.deg, b=Y * u.deg, frame="galactic")
    grid = density[hpx.skycoord_to_healpix(grid_coords)]

    # Poisson noise of the subsample, per pixel, at the sparsest sky.
    sparse = np.percentile(per_pixel, 1)
    return {
        "x_edges": x_edges, "y_edges": y_edges, "density": grid,
        "stats": {
            "healpix_level": level,
            "pixel_area_deg2": area_deg2,
            "n_catalogue": n_total,
            "n_sampled": n_sample,
            "sample_fraction": n_sample / n_total,
            "density_min_per_deg2": float(density[density > 0].min()),
            "density_p01_per_deg2": float(np.percentile(density, 1)),
            "density_median_per_deg2": float(np.median(density)),
            "density_p99_per_deg2": float(np.percentile(density, 99)),
            "density_max_per_deg2": float(density.max()),
            "half_of_stars_within_abs_b_deg": half_within_b,
            "sky_fraction_within_that_band": plane_area_fraction,
            "sampled_count_p01": float(sparse),
            "poisson_noise_at_p01": float(1 / np.sqrt(sparse)) if sparse > 0 else None,
        },
    }


# --------------------------------------------------------------------------
# 7. Famous stars
# --------------------------------------------------------------------------

def locate_named_stars(d, kept, first_fail, cuts, named):
    """Where each named star landed -- or which cut removed it, and why."""
    order = np.argsort(d["source_id"])
    sorted_ids = d["source_id"][order]
    out = {}
    for name, info in named["stars"].items():
        sid = info["source_id"]
        rec = {"source_id": sid, "simbad": info["simbad"]}
        j = np.searchsorted(sorted_ids, sid) if sid is not None else None
        if sid is None or j >= len(sorted_ids) or sorted_ids[j] != sid:
            rec["status"] = "not among the parallax > 10 mas candidates"
            out[name] = rec
            continue
        i = order[j]
        rec.update({
            "g_mag": float(d["phot_g_mean_mag"][i]),
            "bp_rp": float(d["bp_rp"][i]),
            "abs_g": float(d["abs_g"][i]),
            "distance_pc": float(d["distance_pc"][i]),
            "v_tan_km_s": float(d["v_tan"][i]),
            # Gaia's spectrometer measures the line-of-sight speed for bright
            # stars; with it the full speed relative to the Sun follows.
            "radial_velocity_km_s": (float(d["radial_velocity"][i])
                                     if np.isfinite(d["radial_velocity"][i]) else None),
            "total_speed_km_s": (float(np.hypot(d["v_tan"][i], d["radial_velocity"][i]))
                                 if np.isfinite(d["radial_velocity"][i]) else None),
            "ruwe": float(d["ruwe"][i]),
            "kept": bool(kept[i]),
            "status": "on the diagram" if kept[i]
            else f"removed: {cuts[first_fail[i]].label}",
        })
        out[name] = rec
    return out


# --------------------------------------------------------------------------
# 8. The whole pipeline
# --------------------------------------------------------------------------

def analyse(candidates, dwarf, named, sky_counts=None):
    """Run the full analysis. Returns (results for JSON, artefacts for figures)."""
    d = derive(candidates)
    cuts = quality_cuts(d)
    kept, first_fail, flow = apply_cuts(d, cuts)

    colour, absmag = d["bp_rp"], d["abs_g"]
    # Ridge line from stars that cannot be white dwarfs, so the white dwarf
    # sequence -- which crosses the same colours -- never captures the mode.
    not_wd = kept & ~(absmag > WD_LINE[0] + WD_LINE[1] * colour)
    grid, ridge, ridge_n = ridge_line(colour[not_wd], absmag[not_wd])
    delta = absmag - ridge_at(colour, grid, ridge)

    reg = regions(colour, absmag, delta)
    composition = {k: int((m & kept).sum()) for k, m in reg.items()}
    n_final = int(kept.sum())

    fainter = int((absmag[kept] > SUN_ABS_G).sum())
    sun = {"abs_g": SUN_ABS_G, "bp_rp": SUN_BP_RP, "n_fainter": fainter,
           "fraction_fainter": fainter / n_final}

    k = kept
    # Stars removed *only* by the RUWE cut whose photometry is otherwise clean,
    # for comparison with the kept stars in the binary measurement.
    z = np.abs(d["c_star"] / d["c_star_sigma"])
    ruwe_only = ((first_fail == 1) & np.isfinite(colour)
                 & ((d["phot_g_mean_mag"] < EXCESS_G_MIN) | (z < EXCESS_NSIGMA))
                 & (d["phot_bp_mean_mag"] < MAX_BP_MAG))
    bcr = (1.0, 1.8)
    in_bcr = (colour > bcr[0]) & (colour < bcr[1])
    features = {
        "binary_sequence": binary_sequence(colour[k], delta[k], delta[ruwe_only & in_bcr],
                                           colour_range=bcr),
        "red_clump": red_clump(colour[k], absmag[k]),
        "jao_gap": jao_gap(colour[k], absmag[k], delta[k]),
        "white_dwarf_branches": white_dwarf_branches(colour[k & reg["white_dwarfs"]],
                                                     absmag[k & reg["white_dwarfs"]]),
    }

    classes = kinematic_classes(d["v_tan"])
    kin = {"classes": {name: {"n": int((m & kept).sum()),
                              "fraction": float((m & kept).sum() / n_final)}
                       for name, m in classes.items()}}
    ms_single = kept & reg["main_sequence"] & (np.abs(delta) < 0.5)
    rows = parenago_curve(d, ms_single, colour)
    kin["parenago_rows"] = rows
    br = _break_colour(rows) if len(rows) > 8 else None
    okb = np.isfinite(dwarf["B-V"]) & np.isfinite(dwarf["Bp-Rp"])
    bv06 = float(np.interp(0.6, dwarf["B-V"][okb], dwarf["Bp-Rp"][okb]))
    if br:
        xb, (lo, hi), coef, chi2, dof = br
        red = [r["S"] for r in rows if r["colour"] > xb]
        kin["parenago"] = {
            "break_colour": xb, "break_range": [lo, hi],
            "slope_blue": float(coef[1]), "slope_red": float(coef[2]),
            "S_at_break": float(coef[0]), "S_bluest": rows[0]["S"],
            "S_red_mean": float(np.mean(red)), "chi2": chi2, "dof": dof,
            "literature_bp_rp": bv06,
            "literature": "Binney, Dehnen & Bertelli 2000 (MNRAS 318, 658): "
                          "discontinuity near B−V = 0.6",
            "summary": (f"Parenago's discontinuity: dispersion rises from "
                        f"{rows[0]['S']:.0f} to {coef[0]:.0f} km/s, then stops rising at "
                        f"BP−RP = {xb:.2f} (1σ range {lo:.2f}–{hi:.2f}); B−V = 0.6 "
                        f"(Binney et al. 2000) corresponds to BP−RP = {bv06:.2f}"),
        }

    named_stars = locate_named_stars(d, kept, first_fail, cuts, named)

    # Completeness from geometry alone: doubling the distance multiplies the
    # volume by 8, so a complete sample should grow eightfold. Where it grows
    # less, faint stars have started to drop out.
    dist = d["distance_pc"][kept]
    within = {r: int((dist < r).sum()) for r in (12.5, 25, 50, 100)}
    completeness = {
        "within_pc": {f"{r:g}": n for r, n in within.items()},
        "growth_per_doubling": {f"{r:g}->{2 * r:g} pc": (within[2 * r] / within[r]
                                                  if within[r] else None)
                                for r in (12.5, 25, 50)},
        "expected_growth": 8.0,
        "dr2_comparison": "Gaia Collaboration, Babusiaux et al. 2018 (DR2): 3,724 stars "
                          "within 25 pc, 29,683 within 50 pc, 212,728 within 100 pc",
    }

    results = {
        "retrieved_utc": candidates.meta.get("RETRIEVD"),
        "query": candidates.meta.get("QUERY"),
        "n_candidates": len(d["parallax"]),
        "cut_flow": flow,
        "cuts": {"min_parallax_over_error": MIN_PARALLAX_OVER_ERROR, "max_ruwe": MAX_RUWE,
                 "excess_nsigma": EXCESS_NSIGMA, "max_bp_mag": MAX_BP_MAG},
        "n_final": n_final,
        "composition": composition,
        "sun": sun,
        "completeness": completeness,
        "features": features,
        "kinematics": kin,
        "named_stars": named_stars,
        "absent_bright_stars": named.get("absent_bright_stars", {}),
        "systematics": {
            "max_abs_g_error_from_parallax": float(hp.absolute_magnitude_error(
                d["parallax_over_error"][kept]).max()),
            "median_abs_g_error_from_parallax": float(np.median(hp.absolute_magnitude_error(
                d["parallax_over_error"][kept]))),
            "max_zero_point_shift_mag": float(hp.parallax_zero_point_shift(
                d["parallax"][kept]).max()),
        },
    }
    art = {"d": d, "kept": kept, "first_fail": first_fail, "cuts": cuts,
           "delta": delta, "regions": reg, "ridge": (grid, ridge, ridge_n),
           "classes": classes, "dwarf": dwarf}

    if sky_counts is not None:
        sky = sky_density(sky_counts, sky_counts.meta)
        results["sky"] = sky["stats"]
        art["sky"] = sky
    return results, art


def save_results(results: dict) -> None:
    RESULTS_PATH.write_text(json.dumps(results, indent=2, ensure_ascii=False))
