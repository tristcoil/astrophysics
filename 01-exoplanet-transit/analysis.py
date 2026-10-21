"""The transit-detection pipeline.

This is where the science happens. Given a cached Kepler light curve we:

1. **Flatten** it, removing slow stellar variability while protecting the
   transits themselves from the filter.
2. **Search** for a periodic box-shaped dip with a Box Least Squares
   periodogram. The period is something we *find*, never something we look up.
3. **Fold** the light curve on the recovered period so every transit stacks.
4. **Measure** the depth and duration from the folded curve.
5. **Validate** the detection against the ways a transit search fools itself.
6. **Derive** the planet's physical parameters and compare them to the
   literature.

The output is a JSON results file plus in-memory arrays the figure code reuses,
so plotting never repeats the analysis.
"""

from __future__ import annotations

import json
import warnings
from dataclasses import asdict, dataclass, field

import numpy as np

import transit_physics as tp
from config import (
    ASSUMED_ALBEDO,
    LIMB_DARKENING_U1,
    LIMB_DARKENING_U2,
    RESULTS_PATH,
    Target,
)

warnings.filterwarnings("ignore", category=UserWarning)

#: Kepler's long-cadence sampling interval in days (29.4244 minutes).
LONG_CADENCE_DAYS = 29.4244 / 60.0 / 24.0


# --------------------------------------------------------------------------
# Step 1 -- detrending
# --------------------------------------------------------------------------

def _window_in_cadences(window_days: float, cadence_days: float) -> int:
    """Convert a detrending window from days to the odd cadence count SG needs.

    ``LightCurve.flatten`` takes ``window_length`` in **cadences**, not days, and
    a Savitzky-Golay filter requires an odd window. Getting this wrong is the
    single easiest way to silently destroy the signal: a window comparable to the
    transit duration lets the filter follow the transit down and flatten it away.
    """
    n = int(round(window_days / cadence_days))
    if n % 2 == 0:
        n += 1
    return max(n, 11)


def flatten_light_curve(lc, target: Target, transit_mask=None):
    """Remove slow variability, optionally protecting known in-transit points.

    Stars are not constant. Spots rotating in and out of view, pulsations, and
    residual instrumental drift all produce variability on timescales of days --
    much longer than a transit, which lasts hours. A Savitzky-Golay filter fitted
    over a window several times the transit duration tracks the slow part and
    divides it out.

    The ``transit_mask`` argument is what keeps this honest. On the second pass,
    once we know where the transits are, we exclude those cadences from the
    fitting so the filter interpolates *across* each transit instead of bending
    down into it. Without this the measured depth comes out systematically
    shallow.
    """
    window = _window_in_cadences(target.flatten_window_days, LONG_CADENCE_DAYS)
    flat, trend = lc.flatten(
        window_length=window,
        polyorder=2,
        break_tolerance=5,
        mask=transit_mask,
        return_trend=True,
    )
    return flat.remove_nans(), trend, window


# --------------------------------------------------------------------------
# Step 2 -- period search
# --------------------------------------------------------------------------

def bls_period_grid(
    period_min: float,
    period_max: float,
    baseline_days: float,
    min_duration_days: float,
    oversample: float = 4.0,
) -> np.ndarray:
    """Build a period grid fine enough not to smear the transit away.

    This is the step that decides whether a search works at all, and it is easy
    to get silently wrong. If two adjacent trial periods differ by ``dP``, then
    over a baseline ``B`` the predicted mid-transit time drifts by ``(B/P) * dP``.
    Once that drift approaches the transit duration, the transits no longer stack
    when folded and the signal washes out. Requiring the drift to stay below
    ``duration / oversample`` gives::

        dP < P * duration / (oversample * B)

    Because the requirement scales with ``P``, the grid is built uniform in
    *frequency* rather than period, which spaces periods as ``dP ~ P**2`` and so
    satisfies the condition everywhere with no wasted samples at short periods.

    For context on why this matters: Kepler's 1470-day baseline covers about 1760
    orbits of a 0.84-day planet, so a period error of one part in 10,000 already
    displaces the last transit by nearly a full duration.

    Note that astropy's own ``autoperiod`` is stricter still -- it would ask for
    ~1.3e8 trial periods here, which is not computable on a laptop. That is the
    motivation for the two-stage search in :func:`run_bls`.
    """
    df = min_duration_days / (oversample * baseline_days * period_max)
    frequencies = np.arange(1.0 / period_max, 1.0 / period_min, df)
    return (1.0 / frequencies)[::-1]


def run_bls(lc, target: Target, *, discovery_baseline_days: float = 300.0,
            quiet: bool = True):
    """Two-stage Box Least Squares search.

    BLS slides a box-shaped dip of trial duration through the light curve folded
    at every trial period and reports how well a box fits at each. A real
    transiting planet produces a sharp peak at its orbital period. Unlike a
    Lomb-Scargle periodogram, which looks for sinusoids, BLS is matched to the
    flat-bottomed shape a transit actually has, which is why it is the standard
    tool for this job.

    Searching the full period range at the resolution the full 1470-day baseline
    demands is computationally out of reach, so the search is hierarchical, which
    is what real transit pipelines do too:

    **Stage 1 (discovery)** searches the entire period range, but using only the
    first ``discovery_baseline_days`` of data. A shorter baseline needs a much
    coarser grid -- the resolution requirement scales as ``1/B`` -- which brings
    the cost down by more than two orders of magnitude. Roughly 300 days still
    contains hundreds of transits of a short-period planet, so the signal is
    easily strong enough to localise the period.

    **Stage 2 (refinement)** then searches a narrow window around the stage-1
    candidate using *all* the data, which is what delivers the final precision.

    Returns a dict holding both periodograms and the refined parameters.
    """
    time_arr = np.asarray(lc.time.value, dtype=float)
    full_baseline = float(time_arr.max() - time_arr.min())
    min_duration = min(target.durations_days)
    durations = list(target.durations_days)

    # --- Stage 1: coarse search over the whole period range -----------------
    sub = lc[time_arr - time_arr.min() < discovery_baseline_days]
    sub_time = np.asarray(sub.time.value, dtype=float)
    sub_baseline = float(sub_time.max() - sub_time.min())

    coarse_grid = bls_period_grid(
        target.period_min_days, target.period_max_days, sub_baseline, min_duration
    )
    if not quiet:
        print(f"  BLS stage 1: {len(coarse_grid):,} periods over "
              f"{sub_baseline:.0f} d of data")
    coarse = sub.to_periodogram(method="bls", period=coarse_grid, duration=durations)
    candidate = float(coarse.period_at_max_power.value)

    # --- Stage 2: refine around the candidate using every cadence ------------
    # The window is +/-0.3% of the candidate, comfortably wider than stage 1's
    # own grid spacing there, so the true period cannot fall outside it.
    step = candidate * min_duration / (8.0 * full_baseline)
    fine_grid = np.arange(candidate * 0.997, candidate * 1.003, step)
    if not quiet:
        print(f"  BLS stage 2: {len(fine_grid):,} periods over "
              f"{full_baseline:.0f} d of data")
    fine = lc.to_periodogram(method="bls", period=fine_grid, duration=durations)

    return {
        "coarse": coarse,
        "fine": fine,
        "candidate_period": candidate,
        "period": float(fine.period_at_max_power.value),
        "t0": float(fine.transit_time_at_max_power.value),
        "duration": float(fine.duration_at_max_power.value),
        "depth": float(fine.depth_at_max_power.value),
        "discovery_baseline_days": sub_baseline,
        "n_periods_coarse": int(len(coarse_grid)),
        "n_periods_fine": int(len(fine_grid)),
    }


def harmonic_check(bls, best_period: float) -> dict:
    """Compare BLS power at the best period against its harmonics and aliases.

    BLS frequently peaks at an integer multiple or fraction of the true period.
    Folding on 2P still stacks every transit (they just land in two groups), and
    folding on P/2 puts real transits on top of empty phase. Both can look
    convincing.

    We report the power at P/2, P and 2P. If a harmonic carries comparable power
    the folded light curve must be inspected: the tell-tale sign of a 2P error is
    two transits per fold, and of a P/2 error is alternating deep and empty
    transits, which the odd/even test below also catches.
    """
    period_grid = bls.period.value
    power = bls.power.value

    def power_near(p: float) -> float | None:
        if p < period_grid.min() or p > period_grid.max():
            return None
        return float(power[np.argmin(np.abs(period_grid - p))])

    best = power_near(best_period)
    out = {"power_at_best": best}
    for name, factor in (("half", 0.5), ("double", 2.0), ("third", 1 / 3), ("triple", 3.0)):
        p = power_near(best_period * factor)
        out[f"power_at_{name}"] = p
        out[f"ratio_{name}"] = (p / best) if (p is not None and best) else None
    return out


def signal_detection_efficiency(bls, best_period: float) -> float:
    """SDE -- how far the BLS peak stands above the periodogram's own noise.

    ``SDE = (peak_power - mean_power) / std_power``

    computed over the whole periodogram. This is more meaningful than raw BLS
    power, which has arbitrary units that scale with the size of the dataset.

    A warning about the familiar "SDE > 7 is a detection" rule of thumb: it is
    not safe on real Kepler sampling. Block-shuffling this kind of light curve to
    destroy any signal while preserving its noise structure still produces
    maximum SDE values around 7, occasionally approaching 9, purely from noise
    interacting with the mission's data gaps. A genuine detection should stand
    far above that -- both our targets land above 60 -- and marginal candidates
    near the threshold need an empirical null, not folklore.

    We exclude a window around the peak itself from the mean and standard
    deviation, so a very strong signal does not inflate its own baseline.
    """
    period_grid = bls.period.value
    power = np.asarray(bls.power.value, dtype=float)
    peak = float(np.nanmax(power))

    # Mask +/-5% around the peak period so the peak does not bias the baseline.
    away = np.abs(period_grid - best_period) > 0.05 * best_period
    baseline = power[away & np.isfinite(power)]
    if baseline.size < 10:
        return float("nan")
    return float((peak - np.nanmean(baseline)) / np.nanstd(baseline))


# --------------------------------------------------------------------------
# Step 3/4 -- fold and measure
# --------------------------------------------------------------------------

def measure_depth(folded, duration_days: float, period_days: float) -> dict:
    """Measure transit depth from the folded light curve.

    We deliberately do not trust the depth BLS reports, because it depends on the
    trial duration grid. Instead we measure directly:

    * **in-transit** points are those within +/-35% of the transit duration from
      mid-transit. Staying inside the flat bottom rather than spanning the full
      duration avoids diluting the depth with the partial ingress and egress,
      where the planet is only part-way onto the disc.
    * **out-of-transit** baseline uses points between 1.5 and 5 transit durations
      from mid-transit -- far enough to be clear of the wings, close enough to
      avoid any residual long-term trend.

    The depth is the difference of the medians. Medians rather than means, so a
    handful of outliers (cosmic rays, flares) cannot move the answer.
    """
    phase_days = np.asarray(folded.time.value, dtype=float)
    flux = np.asarray(folded.flux.value, dtype=float)
    good = np.isfinite(phase_days) & np.isfinite(flux)
    phase_days, flux = phase_days[good], flux[good]

    half = duration_days / 2.0
    in_transit = np.abs(phase_days) < 0.7 * half
    out_transit = (np.abs(phase_days) > 1.5 * half) & (np.abs(phase_days) < 5.0 * half)

    if in_transit.sum() < 3 or out_transit.sum() < 10:
        return {
            "depth": None, "depth_err": None, "depth_ppm": None,
            "n_in_transit": int(in_transit.sum()),
            "n_out_transit": int(out_transit.sum()),
            "scatter": None, "snr": None, "baseline": None,
        }

    baseline = float(np.median(flux[out_transit]))
    in_level = float(np.median(flux[in_transit]))
    depth = baseline - in_level

    # Scatter of the out-of-transit baseline, via the median absolute deviation
    # scaled to a Gaussian sigma. Robust to the outliers a plain std would chase.
    mad = float(np.median(np.abs(flux[out_transit] - baseline)))
    scatter = 1.4826 * mad

    n_in = int(in_transit.sum())
    # Standard error on the median is ~1.253x that on the mean.
    depth_err = 1.253 * scatter / np.sqrt(n_in) if n_in else None
    snr = tp.transit_snr(depth, scatter, n_in)

    return {
        "depth": depth,
        "depth_err": depth_err,
        "depth_ppm": depth * 1e6,
        "depth_err_ppm": depth_err * 1e6 if depth_err else None,
        "n_in_transit": n_in,
        "n_out_transit": int(out_transit.sum()),
        "scatter": scatter,
        "scatter_ppm": scatter * 1e6,
        "snr": snr,
        "baseline": baseline,
    }


def trapezoid_model(t, t0: float, depth: float, t14: float, t23: float):
    """A trapezoidal transit: flat bottom, linear ingress and egress.

    This is the standard simple description of a transit shape. ``t14`` is the
    total duration from first to fourth contact -- the moment the planet's disc
    first touches the stellar limb, to the moment it fully leaves. ``t23`` is the
    flat-bottomed portion, from second to third contact, when the planet is
    entirely in front of the star. The ingress and egress ramps each last
    ``(t14 - t23) / 2``.

    A real transit is not exactly a trapezoid: limb darkening rounds the bottom,
    because the planet blocks brighter light near the disc centre than near the
    limb. The trapezoid still recovers the two durations well, which is what we
    want it for.
    """
    t = np.asarray(t, dtype=float)
    x = np.abs(t - t0)
    half14 = t14 / 2.0
    half23 = min(t23, t14) / 2.0
    ingress = half14 - half23

    flux = np.ones_like(x)
    flux = np.where(x <= half23, 1.0 - depth, flux)
    if ingress > 1e-12:
        ramp = (x > half23) & (x < half14)
        flux = np.where(ramp, 1.0 - depth * (half14 - x) / ingress, flux)
    return flux


def fit_trapezoid(folded, duration_guess_days: float, depth_guess: float) -> dict:
    """Least-squares trapezoid fit to the folded transit.

    Fitting a shape, rather than reading a width off the curve, is what lets us
    quote the same T14 that published catalogues quote. Measuring the width at
    half depth instead would give roughly ``(T14 + T23) / 2``, which is
    systematically shorter than T14 by half the ingress duration -- a bias of
    10-20% for these systems, and exactly the kind of definitional mismatch that
    makes a comparison against the literature look worse than it is.

    Fitted on binned medians rather than raw cadences: binning first makes the
    fit robust to outliers and keeps the residual scatter roughly uniform.
    """
    from scipy.optimize import least_squares

    phase = np.asarray(folded.time.value, dtype=float)
    flux = np.asarray(folded.flux.value, dtype=float)
    good = np.isfinite(phase) & np.isfinite(flux)
    phase, flux = phase[good], flux[good]

    window = 3.0 * duration_guess_days
    sel = np.abs(phase) < window
    if sel.sum() < 40:
        return {"ok": False, "reason": "too few points near transit"}
    phase, flux = phase[sel], flux[sel]

    # Bin to ~150 points across the window for a stable fit.
    nbins = 150
    edges = np.linspace(-window, window, nbins + 1)
    idx = np.digitize(phase, edges) - 1
    bc, bf, bw = [], [], []
    for b in range(nbins):
        m = idx == b
        n = int(m.sum())
        if n >= 3:
            bc.append(0.5 * (edges[b] + edges[b + 1]))
            bf.append(float(np.median(flux[m])))
            bw.append(np.sqrt(n))
    if len(bc) < 30:
        return {"ok": False, "reason": "too few filled bins"}
    bc = np.asarray(bc)
    bf = np.asarray(bf)
    bw = np.asarray(bw)

    def residual(p):
        t0, depth, t14, frac = p
        t23 = frac * t14
        return (trapezoid_model(bc, t0, depth, t14, t23) - bf) * bw

    # Parametrised by the flat fraction t23/t14 in [0, 1], which keeps the
    # geometry valid for any parameter the optimiser tries.
    p0 = [0.0, max(depth_guess, 1e-6), duration_guess_days, 0.6]
    lo = [-0.25 * duration_guess_days, 1e-7, 0.2 * duration_guess_days, 0.0]
    hi = [0.25 * duration_guess_days, 0.5, 4.0 * duration_guess_days, 1.0]

    try:
        fit = least_squares(residual, p0, bounds=(lo, hi), method="trf", max_nfev=4000)
    except Exception as exc:  # pragma: no cover - defensive
        return {"ok": False, "reason": str(exc)}

    t0, depth, t14, frac = fit.x
    t23 = frac * t14
    resid = residual(fit.x) / bw
    rms = float(np.sqrt(np.mean(resid**2)))

    # --- Correct for the finite exposure time -------------------------------
    # Each Kepler long-cadence sample is a 29.4-minute *integration*, not an
    # instantaneous reading. The measured light curve is therefore the true one
    # convolved with a box of that width, which stretches the transit. For a
    # trapezoid convolved with a box of width w (with w shorter than the flat
    # part) the result is exactly another trapezoid with
    #
    #     T14_observed = T14_true + w        T23_observed = T23_true - w
    #
    # so the correction is a simple subtraction. This matters a great deal for
    # short transits: Kepler-10 b's transit lasts under two hours, so a 29-minute
    # smear inflates its apparent duration by roughly a quarter.
    w = LONG_CADENCE_DAYS
    t14_true = t14 - w
    t23_true = t23 + w
    smearing_valid = bool(t23 > w and t14_true > 0)

    return {
        "ok": bool(fit.success),
        "t0_offset_days": float(t0),
        "depth": float(depth),
        "depth_ppm": float(depth * 1e6),
        "t14_days": float(t14),
        "t14_hours": float(t14 * 24.0),
        "t23_days": float(t23),
        "t23_hours": float(t23 * 24.0),
        "ingress_hours": float((t14 - t23) / 2.0 * 24.0),
        "flat_fraction": float(frac),
        "rms_residual_ppm": rms * 1e6,
        "n_bins": int(len(bc)),
        # Deconvolved values -- what an instantaneous camera would have seen.
        "cadence_hours": float(w * 24.0),
        "smearing_correction_valid": smearing_valid,
        "t14_corrected_days": float(t14_true) if smearing_valid else None,
        "t14_corrected_hours": float(t14_true * 24.0) if smearing_valid else None,
        "t23_corrected_hours": float(t23_true * 24.0) if smearing_valid else None,
        "ingress_corrected_minutes": (
            float((t14_true - t23_true) / 2.0 * 24.0 * 60.0) if smearing_valid else None
        ),
    }


def measure_duration(folded, duration_guess_days: float) -> float | None:
    """Width at half depth -- a model-free cross-check on the trapezoid fit.

    We bin the folded curve finely, find the minimum, and walk outward on each
    side until the flux recovers past the halfway point between the in-transit
    floor and the out-of-transit baseline. For a trapezoid this returns
    ``(T14 + T23) / 2``, so it is deliberately NOT the quantity catalogues call
    the transit duration; it is here as an independent sanity check that the
    fitted shape is not wildly off.
    """
    phase = np.asarray(folded.time.value, dtype=float)
    flux = np.asarray(folded.flux.value, dtype=float)
    good = np.isfinite(phase) & np.isfinite(flux)
    phase, flux = phase[good], flux[good]

    window = 4.0 * duration_guess_days
    sel = np.abs(phase) < window
    if sel.sum() < 20:
        return None
    phase, flux = phase[sel], flux[sel]

    nbins = 120
    edges = np.linspace(-window, window, nbins + 1)
    idx = np.digitize(phase, edges) - 1
    centres, levels = [], []
    for b in range(nbins):
        m = idx == b
        if m.sum() >= 3:
            centres.append(0.5 * (edges[b] + edges[b + 1]))
            levels.append(np.median(flux[m]))
    if len(centres) < 20:
        return None
    centres = np.asarray(centres)
    levels = np.asarray(levels)

    baseline = float(np.median(levels[np.abs(centres) > 1.5 * duration_guess_days]))
    floor = float(np.min(levels))
    if not np.isfinite(baseline) or not np.isfinite(floor) or baseline <= floor:
        return None
    half_level = 0.5 * (baseline + floor)

    below = levels < half_level
    if not below.any():
        return None
    first, last = np.argmax(below), len(below) - 1 - np.argmax(below[::-1])
    return float(centres[last] - centres[first])


# --------------------------------------------------------------------------
# Step 5 -- validation
# --------------------------------------------------------------------------

def per_transit_depths(lc_flat, period: float, t0: float, duration: float) -> dict:
    """Measure the depth of every individual transit separately.

    Doing this once buys two things that a single folded measurement cannot
    give. First, the *empirical* scatter of transit depths, which is the honest
    uncertainty to use when comparing subsets: formal photon-noise errors ignore
    residual systematics and can overstate a depth difference several-fold,
    turning a perfectly ordinary planet into an apparent false positive. Second,
    a check that the depth is stable over the mission rather than drifting with
    detector or quarter, which is how instrumental artefacts betray themselves.
    """
    time = np.asarray(lc_flat.time.value, dtype=float)
    flux = np.asarray(lc_flat.flux.value, dtype=float)
    good = np.isfinite(time) & np.isfinite(flux)
    time, flux = time[good], flux[good]

    epoch = np.round((time - t0) / period).astype(int)
    phase = time - (t0 + epoch * period)
    half = duration / 2.0
    in_t = np.abs(phase) < 0.7 * half
    out_t = (np.abs(phase) > 1.5 * half) & (np.abs(phase) < 5.0 * half)

    epochs, depths, times = [], [], []
    for e in np.unique(epoch):
        sel = epoch == e
        si, so = sel & in_t, sel & out_t
        # Require enough coverage on both sides for the depth to mean anything.
        if si.sum() < 2 or so.sum() < 6:
            continue
        base = float(np.median(flux[so]))
        epochs.append(int(e))
        depths.append(base - float(np.median(flux[si])))
        times.append(float(np.median(time[sel])))

    return {
        "epochs": np.asarray(epochs, dtype=int),
        "depths": np.asarray(depths, dtype=float),
        "times": np.asarray(times, dtype=float),
        "n_usable_transits": len(epochs),
    }


def depth_stability_test(per_transit: dict, n_chunks: int = 6) -> dict:
    """Is the transit depth constant across the mission?

    A real planet blocks the same fraction of light every orbit. Instrumental
    artefacts -- rolling-band noise, changing apertures after each quarterly
    roll, detector drift -- vary with time and CCD. Splitting the transits into
    consecutive chunks and comparing their depths is a cheap way to catch that.

    We compare the scatter of the chunk depths against what the per-transit
    scatter predicts. A ratio near 1 means the variation is just noise.
    """
    depths = per_transit["depths"]
    times = per_transit["times"]
    if depths.size < 4 * n_chunks:
        return {"ok": False, "reason": "too few usable transits"}

    order = np.argsort(times)
    depths, times = depths[order], times[order]
    chunks = np.array_split(depths, n_chunks)

    means = np.array([float(np.median(c)) for c in chunks])
    overall = float(np.median(depths))
    per_transit_scatter = 1.4826 * float(np.median(np.abs(depths - overall)))
    expected = per_transit_scatter / np.sqrt(len(depths) / n_chunks)
    observed = float(np.std(means))

    return {
        "ok": True,
        "chunk_depths_ppm": [m * 1e6 for m in means],
        "mean_depth_ppm": overall * 1e6,
        "observed_scatter_ppm": observed * 1e6,
        "expected_scatter_ppm": expected * 1e6,
        "ratio": float(observed / expected) if expected else None,
        "consistent": bool(expected and observed / expected < 3.0),
    }


def plausibility_gate(radius_ratio: float, stellar_radius_rsun: float | None) -> dict:
    """Is the implied companion small enough to actually be a planet?

    This turns out to be the most decisive discriminator against eclipsing
    binaries, and it is nothing more than a size check. Objects above roughly two
    Jupiter radii are brown dwarfs or stars, not planets -- and crucially,
    degenerate objects stop growing with mass, so anything from about one Jupiter
    radius up to the hydrogen-burning limit sits near the same size. A deep
    eclipse implying several Jupiter radii is therefore a stellar companion,
    whatever the light-curve shape looks like.

    Real Kepler eclipsing binaries produce eclipses of 8% and 30%, implying 2.8
    and 5.5 Jupiter radii respectively -- both comfortably rejected here, whereas
    the V-shape and odd/even tests do not reliably separate them.
    """
    if not stellar_radius_rsun:
        return {"ok": False, "reason": "no stellar radius available"}
    rp_rjup = radius_ratio * stellar_radius_rsun * tp.RSUN_PER_RJUP
    if rp_rjup > 2.0:
        verdict = "TOO LARGE to be a planet -- almost certainly a stellar companion"
    elif rp_rjup > 1.5:
        verdict = "large; possible but at the top of the planetary range"
    else:
        verdict = "within the planetary size range"
    return {
        "ok": True,
        "implied_radius_rjup": float(rp_rjup),
        "threshold_rjup": 2.0,
        "passes": bool(rp_rjup <= 2.0),
        "verdict": verdict,
    }


def odd_even_test(lc_flat, period: float, t0: float, duration: float,
                  per_transit: dict | None = None) -> dict:
    """Compare the depths of odd- and even-numbered transits.

    This is the classic eclipsing-binary discriminator. If what we found is
    really two stars orbiting each other, the "transits" are alternating primary
    and secondary eclipses of unequal depth, and folding on the true period makes
    every other event shallower. A planet produces identical depths.

    It also catches a period error: if BLS locked onto half the true period, the
    odd and even sets contain different physical events and their depths differ.

    We report the depth difference in units of its combined uncertainty. A
    significance above ~3 sigma is a red flag.
    """
    if per_transit is None:
        per_transit = per_transit_depths(lc_flat, period, t0, duration)

    epochs = per_transit["epochs"]
    depths = per_transit["depths"]

    results = {}
    if depths.size < 6:
        for name in ("odd", "even"):
            results[name] = {"depth_ppm": None, "err_ppm": None, "n": 0}
    else:
        for name, parity in (("odd", 1), ("even", 0)):
            sel = (np.abs(epochs) % 2) == parity
            n = int(sel.sum())
            if n < 3:
                results[name] = {"depth_ppm": None, "err_ppm": None, "n": n}
                continue
            d = depths[sel]
            median = float(np.median(d))
            # Empirical scatter of individual transit depths, not a photon-noise
            # estimate. This is the difference between a sane answer and a false
            # alarm: formal errors ignore transit-to-transit systematics and can
            # overstate the odd/even significance by roughly sevenfold.
            scatter = 1.4826 * float(np.median(np.abs(d - median)))
            results[name] = {
                "depth_ppm": median * 1e6,
                "err_ppm": 1.253 * scatter / np.sqrt(n) * 1e6,
                "per_transit_scatter_ppm": scatter * 1e6,
                "n": n,
            }

    odd, even = results["odd"], results["even"]
    if odd["depth_ppm"] is not None and even["depth_ppm"] is not None:
        diff = abs(odd["depth_ppm"] - even["depth_ppm"])
        sigma = np.hypot(odd["err_ppm"], even["err_ppm"])
        mean_depth = 0.5 * (odd["depth_ppm"] + even["depth_ppm"])
        frac = diff / mean_depth if mean_depth else None

        results["difference_ppm"] = diff
        results["significance_sigma"] = float(diff / sigma) if sigma else None
        results["fractional_difference"] = float(frac) if frac is not None else None

        # Two criteria, because either alone misleads. On a very bright target
        # the folded depth is known to ~1 ppm, so a physically trivial 0.1%
        # difference registers as several sigma; conversely on a faint target a
        # genuinely alarming 30% difference may sit below 3 sigma. A real
        # eclipsing binary shows a LARGE fractional difference. We treat the
        # detection as planet-like unless the depths differ both significantly
        # and substantially.
        results["consistent"] = bool(
            sigma
            and frac is not None
            and (diff / sigma < 3.0 or frac < 0.05)
        )
        results["criterion"] = (
            "planet-like unless depths differ by >3 sigma AND >5% fractionally"
        )
    else:
        results["difference_ppm"] = None
        results["significance_sigma"] = None
        results["fractional_difference"] = None
        results["consistent"] = None
    return results


def secondary_eclipse_test(folded, period: float, duration: float,
                           primary_depth: float | None = None) -> dict:
    """Look for a dip at phase 0.5, halfway between transits.

    When the companion passes *behind* the star, we lose whatever light the
    companion itself was contributing. What matters is not whether a secondary
    exists but how deep it is **relative to the primary transit**:

    * A secondary comparable to the primary (tens of percent of it) means the
      companion is itself a star. That is an eclipsing binary, not a planet.
    * A secondary of order 0.1-1% of the primary is exactly what a very hot
      planet produces, from its own thermal glow plus reflected starlight. It is
      a detection *of the planet's atmosphere*, not a red flag. Kepler found
      these for several ultra-hot Jupiters.

    So the verdict here keys on the ratio, not on statistical significance alone:
    with millions of folded cadences even a 7 ppm dip is highly significant.

    Absence of a secondary is consistent with a planet but proves nothing on its
    own -- a sufficiently faint stellar companion also hides.
    """
    phase = np.asarray(folded.time.value, dtype=float)
    flux = np.asarray(folded.flux.value, dtype=float)
    good = np.isfinite(phase) & np.isfinite(flux)
    phase, flux = phase[good], flux[good]

    half_period = period / 2.0
    half_dur = duration / 2.0

    at_secondary = np.abs(np.abs(phase) - half_period) < 0.7 * half_dur
    nearby = (np.abs(np.abs(phase) - half_period) > 1.5 * half_dur) & (
        np.abs(np.abs(phase) - half_period) < 6.0 * half_dur
    )
    if at_secondary.sum() < 3 or nearby.sum() < 10:
        return {"depth_ppm": None, "significance_sigma": None, "detected": None,
                "binary_like": None}

    base = float(np.median(flux[nearby]))
    depth = base - float(np.median(flux[at_secondary]))
    mad = float(np.median(np.abs(flux[nearby] - base)))
    err = 1.253 * 1.4826 * mad / np.sqrt(at_secondary.sum())
    sig = float(depth / err) if err else None
    detected = bool(sig is not None and sig > 5.0)

    ratio = None
    binary_like = None
    interpretation = "no significant secondary eclipse"
    if primary_depth and depth > 0:
        ratio = depth / primary_depth
        # A stellar companion produces a secondary that is a substantial
        # fraction of the primary. 10% is a deliberately conservative line.
        binary_like = bool(detected and ratio > 0.10)
        if detected and binary_like:
            interpretation = (
                f"secondary is {ratio:.0%} of the primary -- too deep for a "
                "planet, suspect an eclipsing binary"
            )
        elif detected:
            interpretation = (
                f"secondary detected at {ratio:.2%} of the primary -- consistent "
                "with thermal emission and reflected light from a hot planet, "
                "not with a stellar companion"
            )
    elif detected:
        interpretation = "secondary detected; no primary depth available for comparison"

    return {
        "depth_ppm": depth * 1e6,
        "err_ppm": err * 1e6,
        "significance_sigma": sig,
        "n": int(at_secondary.sum()),
        "detected": detected,
        "ratio_to_primary": ratio,
        "binary_like": binary_like,
        "interpretation": interpretation,
    }


def transit_shape_ratio(folded, duration: float) -> dict:
    """Ratio of flat-bottom width to total width -- a V-versus-U shape metric.

    A planet much smaller than its star produces a flat-bottomed (U-shaped)
    transit: once the planet is fully on the disc, the blocked fraction stops
    changing. A grazing eclipse, or one by a body comparable in size to the star,
    never has a fully-on phase and produces a V.

    **This test does not actually work at long cadence, and we report it only for
    completeness.** Kepler's 29.4-minute integration smears ingress and egress
    over a large fraction of a short transit, rounding every event toward a V
    regardless of what caused it. Tested against real Kepler eclipsing binaries,
    the metric gives values that overlap with confirmed planets, so it cannot
    separate them. The size gate in :func:`plausibility_gate` is the check that
    genuinely discriminates.
    """
    phase = np.asarray(folded.time.value, dtype=float)
    flux = np.asarray(folded.flux.value, dtype=float)
    good = np.isfinite(phase) & np.isfinite(flux)
    phase, flux = phase[good], flux[good]

    window = 2.0 * duration
    sel = np.abs(phase) < window
    if sel.sum() < 30:
        return {"flat_fraction": None}
    phase, flux = phase[sel], flux[sel]

    nbins = 80
    edges = np.linspace(-window, window, nbins + 1)
    idx = np.digitize(phase, edges) - 1
    centres, levels = [], []
    for b in range(nbins):
        m = idx == b
        if m.sum() >= 3:
            centres.append(0.5 * (edges[b] + edges[b + 1]))
            levels.append(np.median(flux[m]))
    if len(centres) < 20:
        return {"flat_fraction": None}
    centres, levels = np.asarray(centres), np.asarray(levels)

    base = float(np.median(levels[np.abs(centres) > 1.2 * duration]))
    floor = float(np.min(levels))
    if base <= floor:
        return {"flat_fraction": None}

    def width_at(frac: float) -> float | None:
        level = base - frac * (base - floor)
        below = levels < level
        if not below.any():
            return None
        first = np.argmax(below)
        last = len(below) - 1 - np.argmax(below[::-1])
        return float(centres[last] - centres[first])

    w90 = width_at(0.9)   # near the floor -- the flat part
    w50 = width_at(0.5)   # half depth -- the total width
    if not w50 or not w90:
        return {"flat_fraction": None}
    return {
        "width_90pct_days": w90,
        "width_50pct_days": w50,
        "flat_fraction": w90 / w50,
        "shape": "U (flat-bottomed, planet-like)" if w90 / w50 > 0.5 else "V (grazing or binary-like)",
    }


# --------------------------------------------------------------------------
# The pipeline
# --------------------------------------------------------------------------

@dataclass
class TransitResult:
    """Everything the pipeline measured and derived for one target."""

    target: str
    planet: str
    kic: str
    difficulty: str
    n_cadences: int
    baseline_days: float
    n_quarters: int | None
    detrend_window_cadences: int
    detrend_window_days: float

    detected_period_days: float
    detected_t0_bkjd: float
    bls_duration_days: float
    bls_depth: float
    sde: float
    harmonics: dict

    measured_depth_ppm: float | None
    measured_depth_err_ppm: float | None
    measured_duration_hours: float | None
    scatter_ppm: float | None
    snr: float | None
    n_in_transit: int
    n_transits_observed: int

    search_meta: dict = field(default_factory=dict)
    validation: dict = field(default_factory=dict)
    derived: dict = field(default_factory=dict)
    published: dict = field(default_factory=dict)
    comparison: dict = field(default_factory=dict)


def analyse(target: Target, lc, published: dict | None = None, *, quiet: bool = False):
    """Run the whole pipeline for one target.

    Returns ``(result, artefacts)`` where ``artefacts`` holds the intermediate
    light-curve objects the figure code needs, so nothing is recomputed later.
    """
    if not quiet:
        print(f"\n{target.planet}  [{target.difficulty}]")
        print(f"  {len(lc)} cadences")

    # --- Pass 1: detrend blind, search, learn where the transits are ----------
    flat1, trend1, window = flatten_light_curve(lc, target)
    if not quiet:
        print(f"  detrending window {window} cadences ({window * LONG_CADENCE_DAYS:.2f} d)")

    search1 = run_bls(flat1, target, quiet=quiet)
    if not quiet:
        print(f"  pass 1 period = {search1['period']:.6f} d")

    # --- Pass 2: re-detrend with the transits masked out ---------------------
    # Now that we know where the transits fall, we hide them from the filter so
    # it cannot bend down into them and shave the depth. This is why the search
    # runs twice: the first pass exists only to locate the transits.
    fine1 = search1["fine"]
    mask = fine1.get_transit_mask(
        period=fine1.period_at_max_power,
        transit_time=fine1.transit_time_at_max_power,
        duration=fine1.duration_at_max_power,
    )
    flat, trend, _ = flatten_light_curve(lc, target, transit_mask=mask)
    search = run_bls(flat, target, quiet=quiet)

    bls = search["coarse"]      # full period range -- what the periodogram figure shows
    bls_fine = search["fine"]   # the narrow, high-resolution refinement
    period = search["period"]
    t0 = search["t0"]
    bls_duration = search["duration"]
    bls_depth = search["depth"]

    # Significance and harmonics are judged on the coarse periodogram, since that
    # is the one that spans enough period range for "stands above the noise" and
    # "compare against 2P" to mean anything.
    sde = signal_detection_efficiency(bls, search["candidate_period"])
    harmonics = harmonic_check(bls, search["candidate_period"])

    if not quiet:
        print(f"  pass 2 period = {period:.7f} d, SDE = {sde:.1f}")

    # --- Fold and measure ----------------------------------------------------
    folded = flat.fold(period=period, epoch_time=t0)
    depth_info = measure_depth(folded, bls_duration, period)
    fwhm_days = measure_duration(folded, bls_duration)
    trap = fit_trapezoid(folded, bls_duration, depth_info["depth"] or bls_depth)

    # T14 from the trapezoid fit is what catalogues mean by "transit duration",
    # and we use the exposure-time-corrected value so we are comparing like with
    # like against published values derived from properly-integrated models.
    if trap.get("ok"):
        duration_days = trap.get("t14_corrected_days") or trap["t14_days"]
    else:
        duration_days = fwhm_days

    time_arr = np.asarray(flat.time.value, dtype=float)
    span = float(time_arr.max() - time_arr.min())
    n_transits = int(np.floor(span / period)) + 1

    if not quiet and depth_info["depth_ppm"]:
        print(
            f"  depth = {depth_info['depth_ppm']:.1f} +/- "
            f"{depth_info['depth_err_ppm']:.1f} ppm, SNR = {depth_info['snr']:.1f}"
        )

    # --- Validation ----------------------------------------------------------
    per_transit = per_transit_depths(flat, period, t0, bls_duration)
    ratio = float(np.sqrt(depth_info["depth"])) if depth_info["depth"] else None

    validation = {
        "odd_even": odd_even_test(flat, period, t0, bls_duration, per_transit),
        "secondary_eclipse": secondary_eclipse_test(
            folded, period, bls_duration, primary_depth=depth_info["depth"]
        ),
        "shape": transit_shape_ratio(folded, bls_duration),
        "depth_stability": depth_stability_test(per_transit),
        "plausibility": (
            plausibility_gate(ratio, (published or {}).get("st_rad"))
            if ratio else {"ok": False, "reason": "no depth measured"}
        ),
        "n_usable_transits": per_transit["n_usable_transits"],
        "trapezoid": trap,
        "fwhm_hours": (fwhm_days * 24.0) if fwhm_days else None,
        "sde": sde,
        "harmonics": harmonics,
    }

    # --- Derive physical parameters -----------------------------------------
    derived: dict = {}
    comparison: dict = {}
    pub = published or {}
    rstar = pub.get("st_rad")
    mstar = pub.get("st_mass")
    teff = pub.get("st_teff")

    if depth_info["depth"] and rstar and mstar and teff:
        d = tp.derive_all(
            depth=depth_info["depth"],
            period_days=period,
            stellar_radius_rsun=rstar,
            stellar_mass_msun=mstar,
            stellar_teff_k=teff,
            observed_duration_hours=(duration_days * 24.0) if duration_days else None,
            depth_err=depth_info["depth_err"],
            stellar_radius_err_rsun=pub.get("st_raderr1"),
            albedo=ASSUMED_ALBEDO,
        )
        derived = asdict(d)
        derived["limb_darkening"] = tp.limb_darkening_note(
            depth_info["depth"], u1=LIMB_DARKENING_U1, u2=LIMB_DARKENING_U2
        )
        derived["stellar_inputs"] = {
            "st_rad": rstar, "st_mass": mstar, "st_teff": teff,
            "source": "NASA Exoplanet Archive pscomppars",
        }

        # Grade ourselves against the literature.
        if pub.get("pl_orbper"):
            comparison["period"] = {
                "ours": period,
                "published": pub["pl_orbper"],
                "percent_diff": 100.0 * (period - pub["pl_orbper"]) / pub["pl_orbper"],
            }
        if pub.get("pl_rade"):
            comparison["radius_rearth"] = {
                "ours": d.planet_radius_rearth,
                "published": pub["pl_rade"],
                "percent_diff": 100.0
                * (d.planet_radius_rearth - pub["pl_rade"]) / pub["pl_rade"],
            }
            # The limb-darkening-corrected radius assumes a central transit, so
            # it over-corrects for a real, non-central one. Quoting both values
            # brackets the truth rather than pretending to a precision we do not
            # have from a depth measurement alone.
            comparison["radius_rearth_ld_corrected"] = {
                "ours": d.planet_radius_rearth_ld_corrected,
                "published": pub["pl_rade"],
                "percent_diff": 100.0
                * (d.planet_radius_rearth_ld_corrected - pub["pl_rade"]) / pub["pl_rade"],
            }
            lo = min(d.planet_radius_rearth_ld_corrected, d.planet_radius_rearth)
            hi = max(d.planet_radius_rearth_ld_corrected, d.planet_radius_rearth)
            comparison["radius_bracket_rearth"] = {
                "low": lo, "high": hi,
                "published": pub["pl_rade"],
                "published_inside": bool(lo <= pub["pl_rade"] <= hi),
            }
        if pub.get("depth_ppm"):
            comparison["depth_ppm"] = {
                "ours": depth_info["depth_ppm"],
                "published": pub["depth_ppm"],
                "percent_diff": 100.0
                * (depth_info["depth_ppm"] - pub["depth_ppm"]) / pub["depth_ppm"],
            }
        if pub.get("pl_trandur") and duration_days:
            comparison["duration_hours"] = {
                "ours": duration_days * 24.0,
                "published": pub["pl_trandur"],
                "percent_diff": 100.0
                * (duration_days * 24.0 - pub["pl_trandur"]) / pub["pl_trandur"],
            }

    result = TransitResult(
        target=target.name,
        planet=target.planet,
        kic=target.kic,
        difficulty=target.difficulty,
        n_cadences=len(lc),
        baseline_days=span,
        n_quarters=None,
        detrend_window_cadences=window,
        detrend_window_days=window * LONG_CADENCE_DAYS,
        detected_period_days=period,
        detected_t0_bkjd=t0,
        bls_duration_days=bls_duration,
        bls_depth=bls_depth,
        sde=sde,
        harmonics=harmonics,
        search_meta={
            "candidate_period_days": search["candidate_period"],
            "discovery_baseline_days": search["discovery_baseline_days"],
            "n_periods_coarse": search["n_periods_coarse"],
            "n_periods_fine": search["n_periods_fine"],
            "period_min_days": target.period_min_days,
            "period_max_days": target.period_max_days,
            "trial_durations_days": list(target.durations_days),
        },
        measured_depth_ppm=depth_info["depth_ppm"],
        measured_depth_err_ppm=depth_info.get("depth_err_ppm"),
        measured_duration_hours=(duration_days * 24.0) if duration_days else None,
        scatter_ppm=depth_info.get("scatter_ppm"),
        snr=depth_info.get("snr"),
        n_in_transit=depth_info["n_in_transit"],
        n_transits_observed=n_transits,
        validation=validation,
        derived=derived,
        published=pub,
        comparison=comparison,
    )

    artefacts = {
        "raw": lc,
        "flat": flat,
        "trend": trend,
        "bls": bls,
        "bls_fine": bls_fine,
        "search": search,
        "folded": folded,
        "depth_info": depth_info,
        "mask": mask,
    }
    return result, artefacts


def save_results(results: list[TransitResult]) -> None:
    payload = {
        "_generated_utc": __import__("datetime").datetime.now(
            __import__("datetime").timezone.utc
        ).isoformat(timespec="seconds"),
        "results": [asdict(r) for r in results],
    }
    RESULTS_PATH.write_text(json.dumps(payload, indent=2, default=float))
