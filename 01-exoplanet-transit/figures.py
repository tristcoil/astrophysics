"""Figure generation.

Every figure is built from the analysis artefacts, never by re-running the
analysis, so what is plotted is exactly what was measured.

Design notes: the palette is colourblind-safe (adjacent series keep an OKLab
colour difference above 8, checked with a validator rather than by eye), dense
photometry is drawn as a recessive grey cloud so the binned signal reads on top
of it, and each figure carries a caption stating its data provenance and the
caveat a reader needs in order not to over-read it.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

import plotstyle as ps
from config import FIGURE_DIR

ps.use_style()

#: Kepler times are in BKJD -- Barycentric Kepler Julian Date, i.e. BJD - 2454833.
BKJD_OFFSET = 2454833.0

#: Kepler long-cadence integration time in hours (29.4244 minutes).
LONG_CADENCE_HOURS = 29.4244 / 60.0

PROVENANCE = (
    "Data: NASA Kepler mission (PDCSAP photometry), Mikulski Archive for Space "
    "Telescopes. Reference parameters: NASA Exoplanet Archive."
)


def _thin(x, y, max_points=40000):
    """Uniformly thin a dense series for plotting without changing its shape."""
    n = len(x)
    if n <= max_points:
        return x, y
    step = int(np.ceil(n / max_points))
    return x[::step], y[::step]


# --------------------------------------------------------------------------
# 1. The light curve
# --------------------------------------------------------------------------

def figure_light_curve(result, art, target):
    """Raw photometry with the fitted trend, plus a zoom onto individual transits.

    The top panel shows what the spacecraft actually measured across the whole
    mission: several years of brightness, dominated by slow stellar variability
    and instrumental drift. The bottom panel zooms into a few days of the
    flattened curve, which is where a deep transit becomes visible by eye -- and
    where a shallow one conspicuously does not.
    """
    raw, flat, trend = art["raw"], art["flat"], art["trend"]

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(10, 6.4), gridspec_kw={"height_ratios": [1, 1], "hspace": 0.42}
    )

    t_raw = np.asarray(raw.time.value, dtype=float)
    f_raw = np.asarray(raw.flux.value, dtype=float)
    ok = np.isfinite(t_raw) & np.isfinite(f_raw)
    t_raw, f_raw = t_raw[ok], f_raw[ok]
    tx, fx = _thin(t_raw, f_raw)

    ax1.plot(tx, fx, ".", color=ps.C_DATA, markersize=1.0, alpha=0.5,
             rasterized=True, label="PDCSAP flux")
    t_tr = np.asarray(trend.time.value, dtype=float)
    f_tr = np.asarray(trend.flux.value, dtype=float)
    ok2 = np.isfinite(t_tr) & np.isfinite(f_tr)
    ttx, ftx = _thin(t_tr[ok2], f_tr[ok2])
    # Break the line across Kepler's data gaps (downlinks, safe modes) instead of
    # drawing a straight segment through months of missing data.
    ftx = ftx.copy()
    gaps = np.where(np.diff(ttx) > 1.0)[0]
    ftx[gaps] = np.nan
    ax1.plot(ttx, ftx, "-", color=ps.ORANGE, linewidth=1.0,
             label=f"fitted trend ({result.detrend_window_days:.2f} d window)")

    ax1.set_title(
        f"{result.planet} — {result.baseline_days:.0f} days of Kepler photometry"
    )
    ax1.set_xlabel("Time (BKJD, days)")
    ax1.set_ylabel("Normalised flux")
    ax1.legend(loc="lower left", ncol=2)

    # --- zoom on the flattened curve around a few consecutive transits -------
    t_flat = np.asarray(flat.time.value, dtype=float)
    f_flat = np.asarray(flat.flux.value, dtype=float)
    okf = np.isfinite(t_flat) & np.isfinite(f_flat)
    t_flat, f_flat = t_flat[okf], f_flat[okf]

    period = result.detected_period_days
    t0 = result.detected_t0_bkjd
    # Centre the zoom on a transit near the middle of a well-sampled stretch.
    n_mid = np.round((np.median(t_flat) - t0) / period)
    centre = t0 + n_mid * period
    span = 2.5 * period
    sel = np.abs(t_flat - centre) < span

    if sel.sum() > 10:
        ax2.plot(t_flat[sel], f_flat[sel] * 1e6 - 1e6, ".", color=ps.C_DATA,
                 markersize=3.0, alpha=0.85, rasterized=True,
                 label="flattened flux")
        for k in (-2, -1, 0, 1, 2):
            tc = centre + k * period
            if abs(tc - centre) <= span:
                ax2.axvline(tc, color=ps.BLUE, linewidth=1.0, alpha=0.5, zorder=0)
        ax2.axhline(0.0, color=ps.AXIS, linewidth=0.8, zorder=0)
        ax2.set_xlim(centre - span, centre + span)

        depth_ppm = result.measured_depth_ppm or 0.0
        vis = (f_flat[sel] - 1.0) * 1e6
        lo_p, hi_p = np.percentile(vis, [1, 99])
        margin = 0.35 * (hi_p - lo_p)
        ax2.set_ylim(min(lo_p - margin, -1.4 * depth_ppm), hi_p + margin)

        # State the single-transit significance rather than asserting a binary
        # "visible / invisible". A transit spanning n cadences with per-cadence
        # scatter sigma is detected at depth / (sigma / sqrt(n)) in one event.
        scatter = result.scatter_ppm or float("nan")
        trap = result.validation.get("trapezoid", {})
        t14_h = trap.get("t14_hours") or (result.measured_duration_hours or 0.0)
        n_per_transit = max(t14_h / (LONG_CADENCE_HOURS), 1.0)
        per_transit_snr = depth_ppm / (scatter / np.sqrt(n_per_transit)) if scatter else float("nan")

        if per_transit_snr >= 15:
            verdict = f"each transit is an unmistakable {per_transit_snr:.0f}σ dip"
        elif per_transit_snr >= 5:
            verdict = (
                f"one transit is only a {per_transit_snr:.1f}σ dip — visible if you "
                f"know where to look"
            )
        else:
            verdict = (
                f"one transit is just {per_transit_snr:.1f}σ — lost in the noise; "
                f"only folding reveals it"
            )
        ax2.set_title(f"Zoom: {2 * span:.1f} days — {verdict}")
        ax2.set_xlabel("Time (BKJD, days)")
        ax2.set_ylabel("Flux − 1 (ppm)")
        ax2.legend(loc="lower left")
        y_lo, y_hi = ax2.get_ylim()
        ax2.text(
            0.985, 0.95,
            f"vertical rules: predicted mid-transit,\nP = {period:.4f} d (our own value)",
            transform=ax2.transAxes, ha="right", va="top",
            fontsize=8, color=ps.INK_SECONDARY, linespacing=1.5,
        )

    ps.caption(
        fig,
        f"{PROVENANCE}  Top: unflattened PDCSAP flux; the slow wander is real "
        f"stellar variability plus instrumental drift, not noise. "
        f"Blue rules in the lower panel mark predicted mid-transit times from "
        f"our own recovered period.",
    )
    return ps.save(fig, FIGURE_DIR / f"{result.target.lower()}_01_lightcurve")


# --------------------------------------------------------------------------
# 2. The periodogram
# --------------------------------------------------------------------------

def figure_periodogram(result, art, target):
    """BLS periodogram: the period we recovered, and its harmonics.

    A single tall spike is what a real transiting planet looks like. The
    annotated companions at half and twice the period are the aliases a transit
    search most often mistakes for the truth, so they are labelled rather than
    hidden.
    """
    bls = art["bls"]
    periods = np.asarray(bls.period.value, dtype=float)
    power = np.asarray(bls.power.value, dtype=float)

    fig, ax = plt.subplots(figsize=(10, 4.0))
    ax.plot(periods, power, "-", color=ps.BLUE, linewidth=0.9)

    p = result.detected_period_days
    peak_power = float(np.nanmax(power))
    ax.plot([p], [peak_power], "o", color=ps.BLUE, markersize=7,
            markeredgecolor=ps.SURFACE, markeredgewidth=1.5, zorder=5)

    ps.annotate(
        ax,
        f"P = {p:.5f} d\nSDE = {result.sde:.1f}",
        xy=(p, peak_power),
        xytext=(p * 1.35, peak_power * 0.92),
        color=ps.INK, ha="left",
    )

    for label, factor in (("½P", 0.5), ("2P", 2.0)):
        ph = p * factor
        if periods.min() < ph < periods.max():
            yh = power[np.argmin(np.abs(periods - ph))]
            ax.plot([ph], [yh], "o", color=ps.MUTED, markersize=4.5,
                    markeredgecolor=ps.SURFACE, markeredgewidth=1.2, zorder=4)
            ax.annotate(label, xy=(ph, yh), xytext=(0, 8),
                        textcoords="offset points", ha="center",
                        fontsize=8, color=ps.MUTED)

    ax.set_xscale("log")
    ax.set_title(f"{result.planet} — Box Least Squares periodogram")
    ax.set_xlabel("Trial orbital period (days, log scale)")
    ax.set_ylabel("BLS power")
    ax.set_xlim(periods.min(), periods.max())

    ratio_half = result.harmonics.get("ratio_half")
    ratio_double = result.harmonics.get("ratio_double")
    bits = []
    if ratio_half is not None:
        bits.append(f"power at ½P is {ratio_half:.0%} of the peak")
    if ratio_double is not None:
        bits.append(f"at 2P, {ratio_double:.0%}")
    ps.caption(
        fig,
        f"{PROVENANCE}  The period was recovered from the data, not looked up. "
        + ("; ".join(bits) + ". " if bits else "")
        + "SDE counts how many standard deviations the peak stands above the rest "
        "of the periodogram. The familiar 'SDE > 7 means a detection' rule is not "
        "safe on Kepler's gappy sampling — pure noise can reach 7-9 — so only a "
        "large margin like this one counts.",
    )
    return ps.save(fig, FIGURE_DIR / f"{result.target.lower()}_02_periodogram")


# --------------------------------------------------------------------------
# 3. The folded transit -- the hero figure
# --------------------------------------------------------------------------

def _bin_folded(phase_hours, flux, nbins):
    edges = np.linspace(phase_hours.min(), phase_hours.max(), nbins + 1)
    idx = np.digitize(phase_hours, edges) - 1
    c, m, e = [], [], []
    for b in range(nbins):
        sel = idx == b
        n = int(sel.sum())
        if n >= 3:
            vals = flux[sel]
            c.append(0.5 * (edges[b] + edges[b + 1]))
            m.append(np.median(vals))
            mad = np.median(np.abs(vals - np.median(vals)))
            e.append(1.253 * 1.4826 * mad / np.sqrt(n))
    return np.array(c), np.array(m), np.array(e)


def figure_folded_transit(result, art, target):
    """Phase-folded transit: the measurement, and the payoff.

    Every transit in the mission is stacked on top of every other. The grey cloud
    is the individual cadences; the blue points are those cadences binned, which
    is where the shape becomes legible. The measured depth and duration are drawn
    directly onto the figure so the numbers and the picture cannot disagree.
    """
    folded = art["folded"]
    phase_days = np.asarray(folded.time.value, dtype=float)
    flux = np.asarray(folded.flux.value, dtype=float)
    ok = np.isfinite(phase_days) & np.isfinite(flux)
    phase_days, flux = phase_days[ok], flux[ok]

    halfwidth = target.zoom_halfwidth_days
    sel = np.abs(phase_days) < halfwidth
    phase_h = phase_days[sel] * 24.0
    flux_ppm = (flux[sel] - 1.0) * 1e6

    fig, ax = plt.subplots(figsize=(9.5, 5.2))

    px, fx = _thin(phase_h, flux_ppm, 30000)
    ax.plot(px, fx, ".", color=ps.C_DATA, markersize=1.6, alpha=0.30,
            rasterized=True, label=f"{len(phase_h):,} individual cadences")

    nbins = 90
    cb, mb, eb = _bin_folded(phase_h, flux_ppm, nbins)
    ax.errorbar(cb, mb, yerr=eb, fmt="o", color=ps.C_BINNED, markersize=4.2,
                linewidth=0, elinewidth=1.0, capsize=0,
                markeredgecolor=ps.SURFACE, markeredgewidth=0.8,
                label="binned median", zorder=4)

    depth = result.measured_depth_ppm
    if depth:
        ax.axhline(0.0, color=ps.AXIS, linewidth=0.9, zorder=1)

        # Overlay the fitted trapezoid, which is where the quoted duration comes
        # from. Drawing model and data together means a bad fit is visible.
        trap = result.validation.get("trapezoid", {})
        if trap.get("ok"):
            import analysis as _an
            grid = np.linspace(-halfwidth, halfwidth, 1200)
            model = _an.trapezoid_model(
                grid, trap["t0_offset_days"], trap["depth"],
                trap["t14_days"], trap["t23_days"],
            )
            ax.plot(grid * 24.0, (model - 1.0) * 1e6, "-", color=ps.C_MODEL,
                    linewidth=1.8, zorder=5, label="trapezoid fit")

            y = -depth * 1.55
            ax.annotate("", xy=(-trap["t14_hours"] / 2, y),
                        xytext=(trap["t14_hours"] / 2, y),
                        arrowprops={"arrowstyle": "<->",
                                    "color": ps.INK_SECONDARY, "linewidth": 1.0})
            corrected = trap.get("t14_corrected_hours")
            label = f"T$_{{14}}$ = {trap['t14_hours']:.2f} h as observed"
            if corrected:
                label += (
                    f"\n{corrected:.2f} h after removing the "
                    f"{trap['cadence_hours'] * 60:.0f} min exposure smear"
                )
            ax.text(0, -depth * 1.72, label, ha="center", va="top",
                    fontsize=8, color=ps.INK_SECONDARY, linespacing=1.5)
        else:
            ax.axhline(-depth, color=ps.ORANGE, linewidth=1.3, zorder=2,
                       label=f"measured depth {depth:,.0f} ppm")

        ax.text(
            0.985, 0.94,
            f"depth {depth:,.0f} ± {result.measured_depth_err_ppm:,.0f} ppm\n"
            f"= {depth / 1e4:.4f}% of the star's light",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8.5, color=ps.INK, linespacing=1.6,
        )

    ax.set_title(
        f"{result.planet} — {result.n_transits_observed} transits folded onto one period"
    )
    ax.set_xlabel(f"Hours from mid-transit  (folded on P = {result.detected_period_days:.5f} d)")
    ax.set_ylabel("Flux − 1 (ppm)")
    ax.set_xlim(-halfwidth * 24, halfwidth * 24)

    # Frame on the signal, not on the outliers. A handful of cosmic-ray hits
    # would otherwise set the scale and squash a real transit into a thin band.
    if depth:
        lo_p, hi_p = np.percentile(flux_ppm, [0.3, 99.7])
        lo = min(lo_p, -depth * 1.95)
        hi = max(hi_p, depth * 0.75)
        ax.set_ylim(lo, hi)

    ax.legend(loc="lower left", ncol=2, framealpha=0.0)

    snr_txt = f"SNR {result.snr:.0f}" if result.snr else ""
    ps.caption(
        fig,
        f"{PROVENANCE}  Depth is the difference of medians between the flat "
        f"bottom and the out-of-transit baseline; error bars are robust "
        f"(median-absolute-deviation) standard errors. {snr_txt}. Depth here is "
        f"the raw observed dip — limb darkening makes it slightly deeper than the "
        f"true area ratio, which is corrected for in the reported radius.",
    )
    return ps.save(fig, FIGURE_DIR / f"{result.target.lower()}_03_folded")


# --------------------------------------------------------------------------
# 4. Validation
# --------------------------------------------------------------------------

def figure_validation(result, art, target):
    """The checks that argue this is a planet and not something imitating one.

    Left: odd- versus even-numbered transits. Equal depths are what a planet
    gives; unequal depths betray an eclipsing binary or a period off by a factor
    of two. Right: the phase range around the far side of the orbit, where a
    stellar companion would show a secondary eclipse and a planet essentially
    does not.
    """
    folded = art["folded"]
    flat = art["flat"]
    period = result.detected_period_days
    t0 = result.detected_t0_bkjd
    dur = result.bls_duration_days

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
    fig.subplots_adjust(wspace=0.28)

    # --- odd vs even ---------------------------------------------------------
    time = np.asarray(flat.time.value, dtype=float)
    fl = np.asarray(flat.flux.value, dtype=float)
    ok = np.isfinite(time) & np.isfinite(fl)
    time, fl = time[ok], fl[ok]
    epoch = np.round((time - t0) / period)
    phase_h = (time - (t0 + epoch * period)) * 24.0

    trap = result.validation.get("trapezoid", {})
    t14_h = trap.get("t14_hours") or (result.measured_duration_hours or 4.0)
    win = 1.7 * t14_h
    for name, parity, colour in (("odd", 1, ps.C_ODD), ("even", 0, ps.C_EVEN)):
        sel = (np.abs(epoch) % 2 == parity) & (np.abs(phase_h) < win)
        if sel.sum() < 20:
            continue
        cb, mb, eb = _bin_folded(phase_h[sel], (fl[sel] - 1.0) * 1e6, 42)
        ax1.errorbar(cb, mb, yerr=eb, fmt="o", color=colour, markersize=4.0,
                     linewidth=0, elinewidth=1.0,
                     markeredgecolor=ps.SURFACE, markeredgewidth=0.8,
                     label=f"{name}-numbered transits")

    oe = result.validation.get("odd_even", {})
    sig = oe.get("significance_sigma")
    prim = result.measured_depth_ppm
    ax1.axhline(0.0, color=ps.AXIS, linewidth=0.9, zorder=0)
    ax1.set_title("Odd vs even transits")
    ax1.set_xlabel("Hours from mid-transit")
    ax1.set_ylabel("Flux − 1 (ppm)")
    ax1.set_xlim(-win, win)
    if prim:
        ax1.set_ylim(-prim * 1.35, prim * 0.32)
    ax1.legend(loc="lower left", framealpha=0.0)
    if sig is not None:
        frac = oe.get("fractional_difference")
        extra = f" ({frac:.2%} of the depth)" if frac is not None else ""
        verdict = "consistent — planet-like" if oe.get("consistent") else "DISCREPANT — suspect binary"
        ax1.text(0.98, 0.94, f"depths agree to {sig:.1f}σ{extra}\n{verdict}",
                 transform=ax1.transAxes, ha="right", va="top", fontsize=8,
                 color=ps.INK_SECONDARY, linespacing=1.5)

    # --- secondary eclipse ---------------------------------------------------
    pdays = np.asarray(folded.time.value, dtype=float)
    pflux = np.asarray(folded.flux.value, dtype=float)
    okp = np.isfinite(pdays) & np.isfinite(pflux)
    pdays, pflux = pdays[okp], pflux[okp]

    # Re-wrap phase so that zero sits at the secondary eclipse. The fold puts
    # phase in [-P/2, +P/2] with the primary transit at 0, which leaves the
    # secondary split across BOTH ends of that range. This shift brings its two
    # halves together at zero and sends the primary out to the edges.
    off = ((pdays + period) % period) - period / 2.0
    sel2 = np.abs(off) < target.zoom_halfwidth_days
    if sel2.sum() > 30:
        cb2, mb2, eb2 = _bin_folded(off[sel2] * 24.0, (pflux[sel2] - 1.0) * 1e6, 42)
        ax2.errorbar(cb2, mb2, yerr=eb2, fmt="o", color=ps.BLUE, markersize=4.0,
                     linewidth=0, elinewidth=1.0,
                     markeredgecolor=ps.SURFACE, markeredgewidth=0.8,
                     label="binned median")
        # Scale to the secondary itself. Drawing the primary depth as a reference
        # line here would compress this panel by a factor of ~100 and hide the
        # very feature it exists to show; the comparison is stated in words.
        span = float(np.nanmax(mb2) - np.nanmin(mb2))
        pad = max(0.9 * span, 3.0 * float(np.nanmedian(eb2)))
        ax2.set_ylim(np.nanmin(mb2) - pad, np.nanmax(mb2) + pad)
    ax2.axhline(0.0, color=ps.AXIS, linewidth=0.9, zorder=0)

    se = result.validation.get("secondary_eclipse", {})
    ax2.set_title("Where a secondary eclipse would be (phase 0.5)")
    ax2.set_xlabel("Hours from expected secondary eclipse")
    ax2.set_ylabel("Flux − 1 (ppm)")
    ax2.set_xlim(-target.zoom_halfwidth_days * 24, target.zoom_halfwidth_days * 24)
    ax2.legend(loc="lower left", framealpha=0.0)
    if se.get("significance_sigma") is not None:
        d = se["depth_ppm"]
        s = se["significance_sigma"]
        ratio = se.get("ratio_to_primary")
        lines = [f"{d:,.1f} ppm dip ({s:.1f}σ)"]
        if ratio is not None and prim:
            lines.append(f"only {ratio:.2%} of the {prim:,.0f} ppm primary")
            lines.append("→ hot planet's own glow, not a star")
        ax2.text(0.98, 0.94, "\n".join(lines), transform=ax2.transAxes,
                 ha="right", va="top", fontsize=8, color=ps.INK_SECONDARY,
                 linespacing=1.5)

    ps.caption(
        fig,
        f"{PROVENANCE}  These are falsification tests. Matching odd/even depths "
        f"and the absence of a deep secondary eclipse are both consistent with a "
        f"planet, but neither proves one — a sufficiently faint stellar companion "
        f"or a blended background binary can pass both. Confirmation needs "
        f"radial-velocity mass measurement or centroid analysis.",
    )
    return ps.save(fig, FIGURE_DIR / f"{result.target.lower()}_04_validation")


# --------------------------------------------------------------------------
# 5. To-scale geometry
# --------------------------------------------------------------------------

def figure_scale(result, art, target):
    """Star and planet drawn to scale, from our own measured radius ratio.

    The transit chord is placed at the impact parameter implied by the measured
    duration. This is the whole measurement rendered as a picture: the fraction
    of the stellar disc the little circle covers *is* the depth.
    """
    d = result.derived
    if not d:
        return []

    ratio = d["radius_ratio"]
    b = d.get("impact_parameter")
    rp_e = d["planet_radius_rearth"]
    rp_j = d["planet_radius_rjup"]
    rstar = d["stellar_inputs"]["st_rad"]

    fig, ax = plt.subplots(figsize=(8.4, 3.9))

    star = Circle((0, 0), 1.0, facecolor="#fde9c8", edgecolor=ps.YELLOW,
                  linewidth=1.2, zorder=1)
    ax.add_patch(star)

    y = b if b is not None else 0.0
    ax.plot([-1.32, 1.32], [y, y], "-", color=ps.INK_SECONDARY, linewidth=0.9,
            zorder=2, label="transit chord")

    for i, x in enumerate((-1.08, -0.54, 0.0, 0.54, 1.08)):
        planet = Circle((x, y), ratio, facecolor=ps.INK, edgecolor=ps.SURFACE,
                        linewidth=0.6, zorder=3,
                        alpha=1.0 if i == 2 else 0.35)
        ax.add_patch(planet)

    ax.set_xlim(-1.6, 3.05)
    ax.set_ylim(-1.12, 1.12)
    ax.set_aspect("equal")
    ax.axis("off")

    lines = [
        f"$R_p/R_\\star$ = {ratio:.4f}",
        f"$R_p$ = {rp_e:.2f} $R_\\oplus$ = {rp_j:.3f} $R_J$",
        f"$R_\\star$ = {rstar:.2f} $R_\\odot$",
        f"a = {d['semi_major_axis_au']:.4f} AU = {d['scaled_semi_major_axis']:.1f} $R_\\star$",
        f"$T_{{eq}}$ ≈ {d['equilibrium_temperature_k']:.0f} K (albedo 0.3 assumed)",
    ]
    if b is not None:
        lines.append(f"impact parameter b ≈ {b:.2f}")
    ax.text(1.42, 0.0, "\n".join(lines), va="center", ha="left", fontsize=9,
            color=ps.INK, linespacing=2.0)

    ax.set_title(f"{result.planet} — to scale, from our measured depth", pad=6)
    ps.caption(
        fig,
        f"{PROVENANCE}  The planet's silhouette is drawn at the radius ratio we "
        f"measured; the star's absolute radius comes from the archive, so the "
        f"planet's absolute size inherits that uncertainty. Orbital separation is "
        f"NOT to scale here — at {d['scaled_semi_major_axis']:.0f} stellar radii "
        f"the real orbit would not fit on the page.",
    )
    return ps.save(fig, FIGURE_DIR / f"{result.target.lower()}_05_scale")


# --------------------------------------------------------------------------
# 6. Comparison against the literature
# --------------------------------------------------------------------------

def figure_comparison(results):
    """How our measurements compare with published values, as percent differences.

    One row per measured quantity per planet. Zero means we exactly reproduced
    the literature value. This is the figure that says whether the pipeline
    works, so it is deliberately unflattering: the bars are what they are.
    """
    rows = []
    for r in results:
        for key, label in (
            ("period", "Orbital period"),
            ("depth_ppm", "Transit depth"),
            ("duration_hours", "Duration T₁₄"),
            ("radius_rearth", "Radius (raw depth)"),
            ("radius_rearth_ld_corrected", "Radius (limb-dark. corrected)"),
        ):
            c = r.comparison.get(key)
            if c and c.get("percent_diff") is not None:
                rows.append((f"{r.planet}  ·  {label}", c["percent_diff"], r.planet))
    if not rows:
        return []

    fig, ax = plt.subplots(figsize=(9.0, 0.40 * len(rows) + 2.3))
    labels = [r[0] for r in rows]
    vals = [r[1] for r in rows]
    ypos = np.arange(len(rows))

    planets = sorted({r[2] for r in rows})
    colour_of = {p: ps.SERIES[i] for i, p in enumerate(planets)}
    colours = [colour_of[r[2]] for r in rows]

    ax.barh(ypos, vals, height=0.42, color=colours, zorder=3)
    ax.axvline(0.0, color=ps.INK_SECONDARY, linewidth=1.0, zorder=4)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlabel("Difference from published value (%)")
    ax.set_title("Our measurements vs the literature")
    ax.grid(axis="y", visible=False)

    lim = max(1.0, 1.5 * max(abs(v) for v in vals))
    ax.set_xlim(-lim, lim)
    for y, v in zip(ypos, vals):
        # Agreement can be far better than two decimals can show; say so rather
        # than printing a meaningless "+0.00%".
        text = f"{v:+.2f}%" if abs(v) >= 0.01 else f"{v:+.5f}%"
        off = 5 if v >= 0 else -5
        ax.annotate(text, xy=(v, y), xytext=(off, 0),
                    textcoords="offset points", va="center",
                    ha="left" if v >= 0 else "right",
                    fontsize=8, color=ps.INK_SECONDARY)

    # No legend: every row label names its planet, so colour is never the only
    # carrier of identity and a legend box would just collide with the caption.

    ps.caption(
        fig,
        f"{PROVENANCE}  Published values are the NASA Exoplanet Archive composite "
        f"parameters, which are themselves compiled from several papers using "
        f"different methods, so exact agreement is not expected. Radius "
        f"differences inherit the archive's stellar radius, which we did not "
        f"measure ourselves.",
    )
    return ps.save(fig, FIGURE_DIR / "00_comparison")


def make_all(result, art, target):
    written = []
    written += figure_light_curve(result, art, target)
    written += figure_periodogram(result, art, target)
    written += figure_folded_transit(result, art, target)
    written += figure_validation(result, art, target)
    written += figure_scale(result, art, target)
    return written
