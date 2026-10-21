"""The measurements.

1. **Census** -- how the 1.57 million orbits divide between the belt's regions.
2. **The gaps, found blind** -- a search for deficits in the semi-major axis
   histogram that knows nothing about Jupiter. Only afterwards are the gaps
   compared with the resonance locations Kepler's third law predicts.
3. **Are they real?** -- the gaps at every asteroid size, and as the
   catalogue grew from Kirkwood's 91 asteroids to today's.
4. **Resonance lock** -- the Hildas and Trojans, which live *in* resonances
   Jupiter protects instead of clearing, measured directly from their orbits.
5. **Families** -- how much sharper collisional families become in proper
   elements than in osculating ones.
6. **The simulation** (if it has been run) -- where test particles survive.

Every sample cut is logged, and every number that reaches the README is saved
to ``data/results.json``.
"""

from __future__ import annotations

import json
from fractions import Fraction

import numpy as np
from scipy.optimize import least_squares
from scipy.special import erf

import orbit_physics as op
from config import (
    A_RANGE_BELT,
    EARTH_APHELION_AU,
    H_COMPLETE,
    HILDA_A,
    HILDA_E_MAX,
    HILDA_I_MAX_DEG,
    HORIZONS_PATH,
    JUPITER_PERIHELION_AU,
    MARS_APHELION_AU,
    MAX_UNCERTAINTY_U,
    RESULTS_PATH,
    TROJAN_A,
)

# --------------------------------------------------------------------------
# Samples
# --------------------------------------------------------------------------


def quality_mask(d: dict) -> tuple[np.ndarray, list]:
    """Orbits good enough to measure narrow features with, and the cut log.

    The MPC's uncertainty parameter U estimates how far the asteroid's
    predicted position could drift in ten years: U = 0 means under a second
    of arc, U = 9 means the orbit is barely constrained. One-opposition orbits
    may have the eccentricity *assumed* rather than fitted (flag E), which
    makes their semi-major axes guesses.
    """
    log = [("all orbits in MPCORB", int(len(d["a"])))]
    m = d["U"] >= 0
    log.append(("U is a digit (not an E/D/F one-opposition flag)", int(m.sum())))
    m &= d["U"] <= MAX_UNCERTAINTY_U
    log.append((f"U <= {MAX_UNCERTAINTY_U} (position known to ~a few arcsec over a decade)",
                int(m.sum())))
    m &= d["e"] < 1
    log.append(("bound orbit (e < 1)", int(m.sum())))
    return m, log


def load_planets() -> dict:
    return json.loads(HORIZONS_PATH.read_text())


def mean_longitudes(d: dict, t_jd: float) -> dict:
    """Two-body propagate every orbit to ``t_jd``; return lambda, varpi, M."""
    M = op.propagate_mean_anomaly(d["M"], d["n"], t_jd - d["epoch_jd"])
    varpi = op.wrap360(d["node"] + d["peri"])
    lam = op.wrap360(varpi + M)
    return {"M": M, "varpi": varpi, "lam": lam}


def jupiter_longitudes(planets: dict) -> dict:
    j = planets["bodies"]["Jupiter"]
    varpi = op.wrap360(j["node_deg"] + j["peri_deg"])
    return {"lam": float(op.wrap360(varpi + j["mean_anomaly_deg"])), "varpi": float(varpi),
            "theta": float(np.degrees(np.arctan2(j["xyz_au"][1], j["xyz_au"][0])))}


# --------------------------------------------------------------------------
# 1. Census
# --------------------------------------------------------------------------

#: Region boundaries, in semi-major axis, drawn at the resonances that
#: separate them. This is the conventional division of the belt.
REGIONS = [
    ("Hungaria zone", 1.78, 2.00),
    ("inner belt (4:1 to 3:1)", 2.00, 2.50),
    ("middle belt (3:1 to 5:2)", 2.50, 2.82),
    ("outer belt (5:2 to 2:1)", 2.82, 3.28),
    ("Cybele zone (2:1 to 5:3)", 3.28, 3.70),
    ("Hilda zone (around 3:2)", 3.70, 4.20),
    ("Jupiter Trojans (1:1)", TROJAN_A[0], TROJAN_A[1]),
]


def census(d: dict, quality: np.ndarray) -> dict:
    a, e = d["a"], d["e"]
    q = a * (1 - e)
    out = {"total": int(len(a)), "numbered": int((d["number"] > 0).sum()),
           "sections": {"numbered": int((d["section"] == 0).sum()),
                        "unnumbered multi-opposition": int((d["section"] == 1).sum()),
                        "one-opposition": int((d["section"] == 2).sum())},
           "quality_sample": int(quality.sum()), "regions": {}}
    for name, lo, hi in REGIONS:
        sel = (a >= lo) & (a < hi) & (q > EARTH_APHELION_AU)
        out["regions"][name] = int(sel.sum())
    out["near_earth"] = int((q < 1.3).sum())
    out["mars_crossing_not_neo"] = int(((q >= 1.3) & (q < MARS_APHELION_AU)).sum())
    out["beyond_jupiter"] = int((a > 5.5).sum())
    return out


# --------------------------------------------------------------------------
# 2. The gaps, found blind
# --------------------------------------------------------------------------

GAP_WIDTHS = (0.004, 0.006, 0.008, 0.012, 0.016, 0.024, 0.032, 0.048, 0.064)


def detect_gaps(a_sorted: np.ndarray, lo: float, hi: float, *, step: float = 0.0005,
                widths=GAP_WIDTHS, min_sigma: float = 5.0, min_wall: int = 25) -> list:
    """Find gaps without knowing where to look.

    Slide a window of width w along the semi-major axis. A *gap* is a window
    holding fewer asteroids than **both** neighbouring windows of the same
    width would predict -- a dip with a wall on each side. Requiring both
    walls is what keeps the edge of a population (say, the inner edge of the
    belt) from masquerading as a gap. The deficit is measured against the
    *lower* wall, so a gap next to a dense family is not exaggerated.

    Every width in ``widths`` is tried at every position; then the most
    significant (position, width) is kept, every window that overlaps it is
    discarded -- a wide window that merely contains part of a gap already
    found is an echo, not a new gap -- and the next one taken, until nothing
    is left above ``min_sigma``. On pure noise this finds nothing (checked in
    test_algorithms.py).
    """
    xs = np.arange(lo, hi, step)
    sig = np.full((len(widths), len(xs)), -np.inf)
    depth = np.zeros_like(sig)
    for iw, w in enumerate(widths):
        i_lo = np.searchsorted(a_sorted, xs - w / 2)
        i_hi = np.searchsorted(a_sorted, xs + w / 2)
        n_in = i_hi - i_lo
        n_left = i_lo - np.searchsorted(a_sorted, xs - 1.5 * w)
        n_right = np.searchsorted(a_sorted, xs + 1.5 * w) - i_hi
        wall = np.minimum(n_left, n_right).astype(float)
        ok = wall >= min_wall
        sig[iw, ok] = (wall[ok] - n_in[ok]) / np.sqrt(wall[ok])
        depth[iw, ok] = 1 - n_in[ok] / wall[ok]
    found = []
    work = sig.copy()
    while True:
        iw, ix = np.unravel_index(np.argmax(work), work.shape)
        if work[iw, ix] < min_sigma:
            break
        found.append({"a_detect": float(xs[ix]), "width_detect": float(widths[iw]),
                      "sigma": float(sig[iw, ix]), "depth_detect": float(depth[iw, ix])})
        for jw, w in enumerate(widths):
            work[jw, np.abs(xs - xs[ix]) < (widths[iw] + w) / 2] = -np.inf
    return sorted(found, key=lambda g: g["a_detect"])


def _dip(x, mu, h, s):
    """A box of half-width h centred on mu, with edges blurred over s.
    Equals ~1 inside the box when h >> s, and a Gaussian-like dip if h -> 0."""
    r2 = np.sqrt(2.0) * s
    return 0.5 * (erf((x - mu + h) / r2) - erf((x - mu - h) / r2))


def _gap_model(p, x, x0):
    b0, b1, b2, depth, mu, h, s = p
    return (b0 + b1 * (x - x0) + b2 * (x - x0) ** 2) * (1.0 - depth * _dip(x, mu, h, s))


def fit_gap(a_sorted: np.ndarray, x0: float, w0: float, *, bin_au: float = 0.001,
            span: float = 2.5, quadratic: bool = True) -> dict:
    """Refine a gap: centre, width and depth from a fit to the local histogram.

    Model: a smooth (quadratic) background times (1 - depth x smoothed box).
    Fitted by Poisson maximum likelihood -- the counts are counts. The gap's
    width is the full width at half depth of the fitted dip.
    """
    half = span * w0 + 0.01
    edges = np.arange(x0 - half, x0 + half + bin_au / 2, bin_au)
    counts = np.diff(np.searchsorted(a_sorted, edges)).astype(float)
    x = 0.5 * (edges[1:] + edges[:-1])
    bg = np.median(counts[np.abs(x - x0) > w0]) if (np.abs(x - x0) > w0).any() else counts.mean()
    p0 = [max(bg, 1.0), 0.0, 0.0, 0.8, x0, w0 / 2, w0 / 8]
    lo_b = [0.0, -np.inf, -np.inf if quadratic else -1e-12, 0.0, x0 - w0, 2e-4, 2e-4]
    hi_b = [np.inf, np.inf, np.inf if quadratic else 1e-12, 1.0, x0 + w0, 2 * w0, w0]

    def resid(p):
        m = np.maximum(_gap_model(p, x, x0), 1e-9)
        dev = 2 * (m - counts + np.where(counts > 0, counts * np.log(np.maximum(counts, 1e-300) / m), 0))
        return np.sign(counts - m) * np.sqrt(np.maximum(dev, 0))

    fit = least_squares(resid, p0, bounds=(lo_b, hi_b), x_scale="jac")
    b0, b1, b2, depth, mu, h, s = fit.x
    fine = np.linspace(mu - 3 * w0, mu + 3 * w0, 6001)
    prof = depth * _dip(fine, mu, h, s)
    inside = fine[prof >= prof.max() / 2]
    fwhm = float(inside.max() - inside.min()) if len(inside) else float("nan")
    # Formal error on the centre from the Jacobian (statistical only).
    try:
        J = fit.jac
        cov = np.linalg.pinv(J.T @ J)
        mu_err = float(np.sqrt(cov[4, 4]))
    except np.linalg.LinAlgError:
        mu_err = float("nan")
    return {"center": float(mu), "center_err_stat": mu_err, "fwhm": fwhm,
            "depth": float(depth * _dip(mu, mu, h, s)), "edge_blur": float(s),
            "n_bins": int(len(x)), "success": bool(fit.success)}


def fit_gap_robust(a_sorted: np.ndarray, x0: float, w0: float) -> dict:
    """``fit_gap`` repeated over reasonable analysis choices -- window size
    and background shape -- so the quoted spread reflects how much the centre
    depends on them. That spread, not the tiny formal error, is the honest
    uncertainty."""
    base = fit_gap(a_sorted, x0, w0)
    centers = [base["center"]]
    for span in (1.8, 3.5):
        for quad in (True, False):
            try:
                centers.append(fit_gap(a_sorted, x0, w0, span=span, quadratic=quad)["center"])
            except (ValueError, np.linalg.LinAlgError):
                pass
    base["center_spread"] = float(np.ptp(centers) / 2)
    base["center_err"] = float(max(base["center_spread"], base["center_err_stat"]))
    return base


def refine_gap(a_sorted: np.ndarray, x0: float, w0: float, *, max_iter: int = 6) -> dict:
    """Iterate ``fit_gap`` until the centre stops moving, each time centring
    the window on the last fit and sizing it from the fitted width. This
    removes the dependence on where the search happened to start."""
    g = fit_gap(a_sorted, x0, w0)
    for _ in range(max_iter):
        w = float(np.clip(g["fwhm"], 0.004, 0.06)) if np.isfinite(g["fwhm"]) else w0
        new = fit_gap(a_sorted, g["center"], w)
        moved = abs(new["center"] - g["center"])
        g = new
        if moved < 2e-4:
            break
    w = float(np.clip(g["fwhm"], 0.004, 0.06))
    out = fit_gap_robust(a_sorted, g["center"], w)
    out["walls"] = _walls(a_sorted, out["center"], max(out["fwhm"], 0.004))
    # Two estimators of the centre -- the model fit and the model-free
    # deepest point. We quote their mean, and count their disagreement as
    # part of the uncertainty.
    dp = deepest_point(a_sorted, out["center"], w)
    out["center_fit"] = out["center"]
    out["center_deepest"] = dp["center"]
    out["center"] = 0.5 * (out["center_fit"] + dp["center"])
    out["center_err"] = float(max(out["center_err"], abs(out["center_fit"] - dp["center"]) / 2))
    return out


def deepest_point(a_sorted: np.ndarray, center: float, w: float,
                  *, bin_au: float = 0.0005) -> dict:
    """A second, model-free estimate of a gap's centre: its deepest point.

    The histogram is smoothed with a Gaussian a quarter of the gap's width
    (0.001-0.005 AU) and its minimum located. A gap's shoulders -- such as the
    long, gradual decline on the inner side of the 3:1 -- barely move this,
    whereas they can drag a symmetric model fit sideways.
    """
    sigma = float(np.clip(w / 4, 0.001, 0.005))
    half = 2 * w + 4 * sigma
    edges = np.arange(center - half, center + half + bin_au / 2, bin_au)
    x = 0.5 * (edges[1:] + edges[:-1])
    n = np.diff(np.searchsorted(a_sorted, edges)).astype(float)
    k = np.arange(-int(4 * sigma / bin_au), int(4 * sigma / bin_au) + 1)
    kern = np.exp(-0.5 * (k * bin_au / sigma) ** 2)
    sm = np.convolve(n, kern / kern.sum(), mode="same")
    core = np.flatnonzero(np.abs(x - center) <= w)
    j = core[np.argmin(sm[core])]
    # Sub-bin position from a parabola through the minimum and its neighbours.
    if 0 < j < len(x) - 1:
        y0, y1, y2 = sm[j - 1], sm[j], sm[j + 1]
        den = y0 - 2 * y1 + y2
        shift = float(np.clip(0.5 * (y0 - y2) / den, -0.5, 0.5)) if den > 0 else 0.0
    else:
        shift = 0.0
    return {"center": float(x[j] + shift * bin_au), "smoothing_au": sigma}


def _walls(a_sorted: np.ndarray, center: float, w: float) -> dict:
    """Density just outside the gap on each side. A gap the population
    recovers from has walls of comparable height; an *edge* -- where the
    population simply ends -- has one wall far lower than the other."""
    n_l = (np.searchsorted(a_sorted, center - w) - np.searchsorted(a_sorted, center - 2 * w))
    n_r = (np.searchsorted(a_sorted, center + 2 * w) - np.searchsorted(a_sorted, center + w))
    ratio = min(n_l, n_r) / max(n_l, n_r) if max(n_l, n_r) else 0.0
    return {"left": int(n_l), "right": int(n_r), "ratio": float(ratio),
            "is_edge": bool(ratio < EDGE_WALL_RATIO)}


#: Below this ratio of lower to higher wall, a dip is called an edge.
EDGE_WALL_RATIO = 0.25


def identify(center: float, *, max_q: int = 7) -> dict:
    """Name a gap after the simplest period ratio close to it.

    The period ratio P_Jupiter / P_asteroid at the gap is approximated by a
    fraction p/q with q <= max_q; the gap is then compared with the exact
    location of that resonance.
    """
    ratio = float(op.period_ratio(center))
    frac = Fraction(ratio).limit_denominator(max_q)
    p, q = frac.numerator, frac.denominator
    a_res = op.resonance_a(p, q)
    return {"ratio": ratio, "p": p, "q": q, "order": abs(p - q), "a_resonance": a_res,
            "delta_a": center - a_res, "label": f"{p}:{q}"}


def _spike(a_sorted: np.ndarray, center: float, *, bin_au: float = 0.0005,
           around: float = 0.01) -> dict:
    """Peak count in 0.0005 AU bins within +-1 mAU of ``center``, against the
    median bin in the surrounding +-0.01 AU."""
    edges = np.arange(center - around, center + around + bin_au / 2, bin_au)
    n = np.diff(np.searchsorted(a_sorted, edges))
    x = 0.5 * (edges[1:] + edges[:-1])
    near = np.abs(x - center) <= 0.001
    k = np.flatnonzero(near)[np.argmax(n[near])]
    return {"peak_count": int(n[k]), "peak_a": float(x[k]),
            "median_count": float(np.median(n)), "ratio": float(n[k] / max(np.median(n), 1))}


def other_resonances_near(center: float, window: float = 0.012) -> list[dict]:
    """Resonances that are *not* two-body resonances with Jupiter, near a
    feature: two-body resonances with Mars, and three-body resonances with
    Jupiter and Saturn (m_J n_J + m_S n_S + m n = 0, small integers).
    Listed for features the Jupiter resonances cannot explain."""
    out = []
    for p, q in ((1, 2), (2, 5), (1, 3), (3, 7), (2, 3), (3, 5), (4, 7)):
        a = float(op.semi_major_axis(MARS_PERIOD_DAYS * q / p))
        if abs(a - center) < window:
            out.append({"label": f"Mars {p}:{q}", "a": a, "order": abs(p - q)})
    from math import gcd
    for m in (1, 2, 3):
        for mj in range(-7, 8):
            for ms in range(-5, 6):
                if mj == 0 or ms == 0 or gcd(gcd(abs(mj), abs(ms)), m) != 1:
                    continue
                # Asteroid mean motion the resonance requires; must be positive.
                n = (mj * 360.0 / op.JUPITER_PERIOD_DAYS + ms * 360.0 / op.SATURN_PERIOD_DAYS) / m
                if n <= 0:
                    continue
                a = float(op.semi_major_axis(360.0 / n))
                order = abs(mj + ms - m)
                if np.isfinite(a) and abs(a - center) < window and order <= 3:
                    out.append({"label": f"{mj}J{ms:+d}S{-m:+d}", "a": a, "order": order})
    return sorted(out, key=lambda r: abs(r["a"] - center))


#: Mars' mean sidereal period, days (NASA Mars fact sheet).
MARS_PERIOD_DAYS = 686.980


def chance_of_alignment(delta_a: list[float], *, a_range=A_RANGE_BELT, max_q: int = 3,
                        trials: int = 20_000_000, chunk: int = 1_000_000,
                        seed: int = 0) -> dict:
    """How likely is it that gaps placed at random would sit this close to
    low-order resonances?

    Draws ``len(delta_a)`` random positions in ``a_range`` and asks how often
    each lands at least as close (in AU) to a resonance with q <= ``max_q``
    as the observed gaps do, matched from the closest outward.
    """
    fracs = sorted({Fraction(p, q) for q in range(1, max_q + 1) for p in range(q, 6 * q)
                    if Fraction(p, q) >= 1})
    res_a = np.array(sorted(op.resonance_a(f.numerator, f.denominator) for f in fracs))
    rng = np.random.default_rng(seed)
    obs = np.sort(np.abs(np.asarray(delta_a)))
    n_hits = 0
    for _ in range(trials // chunk):
        x = rng.uniform(*a_range, size=(chunk, len(delta_a)))
        dist = np.min(np.abs(x[..., None] - res_a[None, None, :]), axis=-1)
        dist.sort(axis=1)
        n_hits += int(np.all(dist <= obs[None, :], axis=1).sum())
    p = n_hits / trials
    # The same question answered for each gap alone, for intuition: the
    # fraction of the range lying that close to some resonance.
    grid = np.linspace(*a_range, 200_001)
    gdist = np.min(np.abs(grid[:, None] - res_a[None, :]), axis=1)
    single = [float((gdist <= o).mean()) for o in obs]
    return {"probability": float(p), "hits": n_hits, "trials": trials,
            "n_gaps": len(delta_a), "single_gap_chances": single,
            "upper_limit_if_zero": float(3.0 / trials) if p == 0 else None,
            "resonances_considered": [f"{f.numerator}:{f.denominator}" for f in fracs
                                      if a_range[0] <= op.resonance_a(f.numerator, f.denominator)
                                      <= a_range[1]]}


#: Low-order resonances measured individually (where each is looked for is
#: predicted, but the fit finds its own centre).
NAMED_RESONANCES = [(4, 1), (7, 2), (3, 1), (8, 3), (5, 2), (7, 3), (9, 4), (2, 1)]


def measure_gaps(d: dict, quality: np.ndarray, proper: dict | None) -> tuple[dict, dict]:
    a_osc = np.sort(d["a"][quality])
    detections = detect_gaps(a_osc, 1.85, 3.6)
    for g in detections:
        g.update(refine_gap(a_osc, g["a_detect"], g["width_detect"]))
        g.update(identify(g["center"]))
        g["kind"] = "edge" if g["walls"]["is_edge"] else "gap"
        # A Jupiter resonance is a credible explanation only if it is of
        # modest order and lies inside the fitted gap.
        g["credible"] = bool(g["order"] <= 6
                             and abs(g["delta_a"]) < max(0.005, g["fwhm"] / 2))
        if not g["credible"]:
            g["other_candidates"] = other_resonances_near(g["center"])
    # For a feature no Jupiter resonance explains, look for independent
    # evidence at each candidate: in proper elements, asteroids *captured* in
    # a resonance all get the resonance's semi-major axis, so a real
    # resonance shows up as a narrow spike there.
    a_prop_all = np.sort(proper["a"][~proper["secres"]]) if proper is not None else None
    for g in detections:
        for cand in g.get("other_candidates", []):
            if a_prop_all is not None:
                cand["proper_spike"] = _spike(a_prop_all, cand["a"])
    # Two search windows can converge on the same gap once refined; keep the
    # more significant detection of each.
    unique: list[dict] = []
    for g in sorted(detections, key=lambda g: -g["sigma"]):
        if not any(abs(g["center"] - u["center"]) < max(u["fwhm"], g["fwhm"]) / 2
                   for u in unique):
            unique.append(g)
    detections = sorted(unique, key=lambda g: g["center"])

    # The four most significant interior gaps (edges excluded) are the
    # classic Kirkwood gaps. How improbable is it that four gaps placed at
    # random would sit this close to resonances with q <= 3?
    interior = [g for g in detections if g["kind"] == "gap"]
    top = sorted(interior, key=lambda g: -g["sigma"])[:4]
    chance = chance_of_alignment([g["center"] - op.resonance_a(g["p"], g["q"]) for g in top])

    # Each named resonance, fitted in osculating and in proper elements.
    named = []
    a_prop = None
    if proper is not None:
        a_prop = np.sort(proper["a"][~proper["secres"]])
    for p, q in NAMED_RESONANCES:
        a_res = op.resonance_a(p, q)
        w0 = 0.03 if abs(p - q) <= 3 else 0.012
        row = {"label": f"{p}:{q}", "order": abs(p - q), "a_predicted": a_res,
               "osculating": refine_gap(a_osc, a_res, w0)}
        if a_prop is not None:
            row["proper"] = refine_gap(a_prop, a_res, w0)
        named.append(row)

    results = {"detections": detections, "classic_four": [g["label"] for g in top],
               "chance_alignment": chance, "named": named,
               "jupiter_period_days": float(op.JUPITER_PERIOD_DAYS),
               "jupiter_a_mean": op.jupiter_semi_major_axis()}
    arrays = {"a_osc_sorted": a_osc, "a_proper_sorted": a_prop}
    return results, arrays


# --------------------------------------------------------------------------
# 3. Are the gaps real?
# --------------------------------------------------------------------------

#: Absolute-magnitude slices. Bright asteroids (low H) are big; H = 13 is
#: roughly 10 km across at a typical albedo, H = 18 roughly 1 km.
H_BINS = [(-5.0, 13.0), (13.0, 15.0), (15.0, 16.5), (16.5, 18.0), (18.0, 30.0)]


def _core_and_walls(a_sorted: np.ndarray, center: float, fwhm: float) -> tuple[int, float]:
    """Count in the gap's core -- its central half, |a - centre| < FWHM/4 -- and
    the count expected there from the walls just outside the gap (the mean of
    two windows of the same width, one each side, 1-1.5 FWHM from the centre).
    """
    w = fwhm / 2
    n_in = int(np.searchsorted(a_sorted, center + w / 2) - np.searchsorted(a_sorted, center - w / 2))
    n_l = np.searchsorted(a_sorted, center - fwhm) - np.searchsorted(a_sorted, center - fwhm - w)
    n_r = np.searchsorted(a_sorted, center + fwhm + w) - np.searchsorted(a_sorted, center + fwhm)
    return n_in, 0.5 * float(n_l + n_r)


def _depth_at(a_sorted: np.ndarray, center: float, fwhm: float, min_expected: float = 20):
    """1 - (count in the gap's core) / (count the walls predict there)."""
    n_in, expected = _core_and_walls(a_sorted, center, fwhm)
    return float(1 - n_in / expected) if expected >= min_expected else float("nan")


def gaps_by_size_and_era(d: dict, quality: np.ndarray, gaps: dict) -> dict:
    """The main gaps' depths in every brightness slice, and in the catalogue
    as it stood at five points in its history (numbering order)."""
    classic = {row["label"]: row["osculating"] for row in gaps["named"]
               if row["label"] in ("3:1", "5:2", "7:3", "9:4")}
    by_h = []
    H = d["H"]
    for lo, hi in H_BINS:
        sel = quality & (H >= lo) & (H < hi)
        a_s = np.sort(d["a"][sel])
        row = {"H_range": [lo, hi], "n": int(sel.sum()),
               "n_belt": int(((a_s > 2.0) & (a_s < 3.3)).sum()),
               "diameter_km_at_albedo_0.14": [float(op.diameter_km(hi, 0.14)),
                                              float(op.diameter_km(lo, 0.14))] if lo > -5 else
                                             [float(op.diameter_km(hi, 0.14)), None]}
        for lab, g in classic.items():
            row[lab] = _depth_at(a_s, g["center"], g["fwhm"])
        by_h.append(row)

    numbered = d["number"] > 0
    order = np.argsort(d["number"][numbered])
    a_num = d["a"][numbered][order]
    eras = []
    for n_first in (91, 1_000, 10_000, 100_000, int(numbered.sum())):
        a_s = np.sort(a_num[:n_first])
        row = {"first_n": n_first}
        for lab, g in classic.items():
            c, core = g["center"], g["fwhm"] / 2
            row[lab + "_inside"] = int(np.searchsorted(a_s, c + core / 2)
                                      - np.searchsorted(a_s, c - core / 2))
            # What a gap-free belt would put in the core: the density in a
            # band 0.1 AU wide on each side, just outside the gap itself.
            half = g["fwhm"] / 2
            flank = ((np.searchsorted(a_s, c - half) - np.searchsorted(a_s, c - half - 0.1))
                     + (np.searchsorted(a_s, c + half + 0.1) - np.searchsorted(a_s, c + half)))
            row[lab + "_expected"] = float(flank * core / 0.2)
        eras.append(row)

    # Where discovery stops: the faintest H reached, as a function of a. Only
    # orbits that never cross Mars' (q > 1.666 AU) count here: objects that
    # swing close to Earth get discovered at far fainter H, which would mix a
    # second, different selection effect into the curve.
    edges = np.arange(1.8, 5.5 + 1e-9, 0.1)
    faint = []
    q_ok = d["a"] * (1 - d["e"]) > MARS_APHELION_AU
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (quality & q_ok & (d["a"] >= lo) & (d["a"] < hi) & np.isfinite(H) & (H < 30))
        if sel.sum() >= 200:
            faint.append([float((lo + hi) / 2), float(np.percentile(H[sel], 95)), int(sel.sum())])
    return {"by_H": by_h, "by_era": eras, "faint_limit_95pct": faint,
            "h_complete_assumed": H_COMPLETE}


# --------------------------------------------------------------------------
# 4. Resonance lock: Hildas and Trojans
# --------------------------------------------------------------------------

def resonant_populations(d: dict, quality: np.ndarray, planets: dict) -> tuple[dict, dict]:
    t = planets["epoch_jd_tdb"]
    L = mean_longitudes(d, t)
    J = jupiter_longitudes(planets)
    a, e, inc, H = d["a"], d["e"], d["incl"], d["H"]

    hilda = ((a > HILDA_A[0]) & (a < HILDA_A[1]) & (e < HILDA_E_MAX)
             & (inc < HILDA_I_MAX_DEG) & quality)
    sigma = op.resonant_angle(3, 2, L["lam"], J["lam"], L["varpi"])
    # Control: the Cybele zone just inside, which is *not* in the 3:2.
    control = (a > 3.30) & (a < 3.60) & quality
    hist_edges = np.linspace(-180, 180, 25)

    # True positions, for the triangle: where each Hilda actually is, measured
    # as an angle from Jupiter around the Sun.
    x, y, _ = op.elements_to_xyz(a[hilda], e[hilda], inc[hilda], d["node"][hilda],
                                 d["peri"][hilda], L["M"][hilda])
    theta = op.wrap180(np.degrees(np.arctan2(y, x)) - J["theta"])

    trojan = (a > TROJAN_A[0]) & (a < TROJAN_A[1]) & quality
    dlam = op.longitude_from_jupiter(L["lam"], J["lam"])
    ratios = []
    for hmax in (99.0, 16.0, 15.0, 14.0, 13.0, 12.0):
        sel = trojan & (H < hmax)
        n4, n5 = int((sel & (dlam > 0)).sum()), int((sel & (dlam < 0)).sum())
        r = n4 / n5 if n5 else float("nan")
        # Binomial error on the ratio, from the fraction f = n4 / (n4 + n5).
        n = n4 + n5
        f = n4 / n if n else float("nan")
        err = float(np.sqrt(f * (1 - f) / n) / (1 - f) ** 2) if n else float("nan")
        ratios.append({"H_max": hmax, "L4": n4, "L5": n5, "ratio": r, "ratio_err": err})

    res = {
        "hildas": {
            "n": int(hilda.sum()),
            "frac_sigma_within_90": float((np.abs(sigma[hilda]) < 90).mean()),
            "median_abs_sigma": float(np.median(np.abs(sigma[hilda]))),
            "control_n": int(control.sum()),
            "control_frac_sigma_within_90": float((np.abs(sigma[control]) < 90).mean()),
            "theta_peaks_deg": _three_peaks(theta),
        },
        "trojans": {
            "n": int(trojan.sum()),
            "L4_L5_by_H": ratios,
            "median_abs_dlam": float(np.median(np.abs(dlam[trojan]))),
            "literature_ratio": 1.45, "literature_ratio_err": 0.05,
            "literature": "Vokrouhlicky et al. 2024 (AJ 167, 138), Catalina Sky Survey, H < 15",
        },
        "jupiter": J,
    }
    arrays = {"hilda_sigma": sigma[hilda], "control_sigma": sigma[control],
              "hilda_theta": theta, "trojan_dlam": dlam[trojan], "trojan_H": H[trojan],
              "hist_edges": hist_edges}
    return res, arrays


def _three_peaks(theta: np.ndarray) -> list[float]:
    """Directions (deg from Jupiter) of the three densest 30-degree sectors,
    one per 120-degree third of the circle centred on the expected vertices."""
    peaks = []
    for centre in (-60.0, 60.0, 180.0):
        rel = op.wrap180(theta - centre)
        sel = np.abs(rel) < 60
        h, edges = np.histogram(rel[sel], bins=24, range=(-60, 60))
        k = np.argmax(np.convolve(h, np.ones(5), mode="same"))
        peaks.append(float(op.wrap180(centre + 0.5 * (edges[k] + edges[k + 1]))))
    return peaks


# --------------------------------------------------------------------------
# 5. Families: osculating vs proper
# --------------------------------------------------------------------------

def _iqr(x: np.ndarray) -> float:
    q75, q25 = np.nanpercentile(x, [75, 25])
    return float(q75 - q25)


def families(d: dict, proper: dict, famtab, n_top: int = 12) -> dict:
    """The largest families, and how much tighter they are in proper elements.

    For each family, numbered members are matched to their MPC osculating
    orbits by number, and the spread (interquartile range) of e and sin i is
    compared between the two element sets. A ratio well above 1 means the
    osculating elements scatter the family that the proper elements hold
    together.
    """
    fam = proper["family"]
    num = proper["number"]
    in_belt = (proper["a"] > 1.78) & (proper["a"] < 3.7)
    idx_of_number = np.full(int(d["number"].max()) + 1, -1, dtype=np.int64)
    idx_of_number[d["number"][d["number"] > 0]] = np.flatnonzero(d["number"] > 0)

    rows = []
    for _, r in famtab.head(40).iterrows():
        fid = int(r["family"])
        members = (fam == fid) & in_belt
        if members.sum() < 100:
            continue
        numbered = members & (num > 0) & (num < len(idx_of_number))
        j = idx_of_number[num[numbered]]
        ok = j >= 0
        sini_osc = np.sin(np.radians(d["incl"][j[ok]]))
        e_osc = d["e"][j[ok]]
        sini_prop = proper["sini"][numbered][ok]
        e_prop = proper["e"][numbered][ok]
        rows.append({
            "id": fid, "name": r["name"], "n_members": int(members.sum()),
            "a_median": float(np.median(proper["a"][members])),
            "e_median": float(np.nanmedian(proper["e"][members])),
            "sini_median": float(np.median(proper["sini"][members])),
            "sini_iqr_osc": _iqr(sini_osc), "sini_iqr_proper": _iqr(sini_prop),
            "e_iqr_osc": _iqr(e_osc), "e_iqr_proper": _iqr(e_prop),
            "secular_resonant": bool(r["e_min"] > 0.6),
        })
        if len(rows) >= n_top:
            break
    for row in rows:
        row["sini_sharpening"] = row["sini_iqr_osc"] / row["sini_iqr_proper"]
        row["e_sharpening"] = (row["e_iqr_osc"] / row["e_iqr_proper"]
                               if row["e_iqr_proper"] > 0 else float("nan"))
    belt = in_belt & (proper["a"] > 2.0) & (proper["a"] < 3.3)
    return {"top": rows,
            "fraction_in_families": float((fam[belt] > 0).mean()),
            "n_families_in_table": int(len(famtab)),
            "n_proper": int(len(proper["a"]))}


# --------------------------------------------------------------------------
# Snapshot
# --------------------------------------------------------------------------

def snapshot(d: dict, planets: dict) -> dict:
    """Heliocentric ecliptic positions of every object on the snapshot date.

    Each orbit is carried from its epoch (for 99.4% of objects, 2026 June 9)
    to the snapshot by two-body motion, ignoring the planets' pull over those
    112 days -- an error far below a pixel on the map.
    """
    t = planets["epoch_jd_tdb"]
    M = op.propagate_mean_anomaly(d["M"], d["n"], t - d["epoch_jd"])
    x, y, z = op.elements_to_xyz(d["a"], d["e"], d["incl"], d["node"], d["peri"], M)
    return {"x": x.astype(np.float32), "y": y.astype(np.float32), "z": z.astype(np.float32),
            "t_jd": t, "dt_days_median": float(np.median(t - d["epoch_jd"]))}


# --------------------------------------------------------------------------
# 6. Simulation summary
# --------------------------------------------------------------------------

def simulation_summary(sim: dict | None) -> dict | None:
    if sim is None:
        return None
    a0 = sim["ic_a"]
    removed = np.isfinite(sim["t_removed"])
    out = {"n_particles": int(len(a0)), "years": float(sim["t_years"][-1]),
           "removed": int(removed.sum()),
           "removed_mars": int((sim["cause"] == 1).sum()),
           "removed_jupiter": int((sim["cause"] == 2).sum()),
           "removed_unbound": int((sim["cause"] == 3).sum()),
           "energy_drift_max": float(sim["energy_drift_max"]),
           "wall_minutes": float(sim["wall_seconds"]) / 60, "zones": []}
    for p, q in ((4, 1), (3, 1), (5, 2), (7, 3), (2, 1)):
        c = op.resonance_a(p, q)
        zone = np.abs(a0 - c) < 0.015
        if zone.sum():
            t_rem = sim["t_removed"][zone]
            out["zones"].append({
                "label": f"{p}:{q}", "a": c, "n": int(zone.sum()),
                "removed_frac": float(np.isfinite(t_rem).mean()),
                "median_removal_kyr": (float(np.nanmedian(t_rem)) / 1e3
                                       if np.isfinite(t_rem).any() else None)})
    far = np.ones(len(a0), bool)
    for p, q in ((4, 1), (7, 2), (3, 1), (8, 3), (5, 2), (7, 3), (9, 4), (2, 1)):
        far &= np.abs(a0 - op.resonance_a(p, q)) > 0.03
    # The innermost belt is a special case: close enough to Mars that a
    # modest eccentricity makes an orbit Mars-crossing, and home to the nu6
    # secular resonance. Reported separately so it doesn't blur the contrast
    # between resonant and non-resonant particles elsewhere.
    inner = a0 < 2.3
    out["inner_edge"] = {"a_max": 2.3, "n": int(inner.sum()),
                         "removed_frac": float(removed[inner].mean())}
    out["away_from_resonances"] = {"a_min": 2.3, "n": int((far & ~inner).sum()),
                                   "removed_frac": float(removed[far & ~inner].mean())}
    return out


# --------------------------------------------------------------------------
# Everything
# --------------------------------------------------------------------------

def analyse(d: dict, proper: dict | None, famtab, sim: dict | None = None,
            quiet: bool = False) -> tuple[dict, dict]:
    planets = load_planets()
    quality, cut_log = quality_mask(d)
    results = {"cuts": cut_log, "snapshot_jd": planets["epoch_jd_tdb"]}
    arrays: dict = {"quality": quality}

    if not quiet:
        print("  census")
    results["census"] = census(d, quality)

    if not quiet:
        print("  searching for gaps")
    results["gaps"], arr = measure_gaps(d, quality, proper)
    arrays.update(arr)

    if not quiet:
        print("  gaps by size and by catalogue era")
    results["reality"] = gaps_by_size_and_era(d, quality, results["gaps"])

    if not quiet:
        print("  Hildas and Trojans")
    results["resonant"], arr = resonant_populations(d, quality, planets)
    arrays.update(arr)

    if proper is not None and famtab is not None:
        if not quiet:
            print("  families")
        results["families"] = families(d, proper, famtab)

    arrays["snapshot"] = snapshot(d, planets)
    arrays["planets"] = planets
    results["simulation"] = simulation_summary(sim)
    return results, arrays


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, float) and not np.isfinite(o):
        return None
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def save_results(results: dict) -> None:
    RESULTS_PATH.write_text(json.dumps(_jsonable(results), indent=2))
