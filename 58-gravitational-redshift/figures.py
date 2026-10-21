"""Every figure of the project, drawn from cached results (see run.py)."""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config as C
import plotstyle as S

S.use_style()
FIG = C.FIGURE_DIR


# ---------------------------------------------------------------------------
# GPS clock budget
# ---------------------------------------------------------------------------

def gps_budget(df: pd.DataFrame, summary: dict) -> list:
    import gps_budget as GB
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11.5, 4.8), gridspec_kw={"width_ratios": [1.45, 1]})
    alt = np.linspace(0, 40_000, 800)
    rates = GB.circular_rates((alt + 6_371.0) * 1e3)
    ax.axhline(0, color=S.AXIS, lw=0.8)
    ax.plot(alt, rates["gravity"], color=S.GRAVITY, lw=2, label="higher up: runs fast (gravity)")
    ax.plot(alt, rates["speed"], color=S.SPEED, lw=2, label="moving: runs slow (speed)")
    ax.plot(alt, rates["total"], color=S.INK, lw=2.2, label="net, vs. a clock on the ground")
    x0 = summary["crossover_altitude_km"]
    ax.plot([x0], [0], "o", color=S.INK, ms=6, mec=S.SURFACE, zorder=6)
    S.annotate(ax, f"no net effect at {x0:,.0f} km up\n(1.5 Earth radii from the centre)", (x0, 0),
               (7_500, 16), va="top")
    # E14/E18 sweep between perigee and apogee, drawn under the points
    rp, ra = 23_450.0, 32_500.0
    sweep = GB.circular_rates(np.linspace(rp, ra, 50) * 1e3)["total"]
    ax.plot(np.linspace(rp, ra, 50) - 6371, sweep, color=S.BLUE_300, lw=6, alpha=0.6,
            solid_capstyle="round", zorder=3)
    # labels stacked below the curve, ordered so that leader lines do not cross
    label_pos = {"ISS": (1500, -30), "GLONASS": (27_000, 6), "GPS": (27_000, 12),
                 "BeiDou MEO": (27_000, 18), "Galileo": (27_000, 24), "geostationary": (33_000, 37)}
    for name, r_km in GB.REFERENCE_ORBITS.items():
        a = r_km - 6_371.0
        y = GB.circular_rates(np.array([r_km * 1e3]))["total"][0]
        ax.plot([a], [y], "o", color=S.INK, ms=5, mec=S.SURFACE, zorder=7)
        S.annotate(ax, f"{name} {y:+.1f} µs/day", (a, y), label_pos[name], fontsize=7.5)
    S.annotate(ax, "E14/E18: swings along this band\nonce per orbit", (25_500, sweep[25]), (27_000, 30),
               fontsize=7.5)
    ax.set_xlim(0, 40_000)
    ax.set_ylim(-36, 62)
    ax.set_xlabel("height above the ground (km), circular orbit")
    ax.set_ylabel("rate relative to a ground clock (µs per day)")
    ax.set_title("How much faster a clock in orbit runs")
    ax.legend(loc="lower right", fontsize=7.5)

    per_sat = df.groupby("prn").agg(drift=("published_drift_us_day", "median"),
                                    gain=("predicted_gain_us_day", "median")).reset_index()
    y = np.arange(len(per_sat))
    bx.axvline(0, color=S.AXIS, lw=0.8)
    bx.scatter(per_sat["drift"], y, s=22, color=S.BLUE, marker="o", zorder=3,
               label="published drift (offset applied)")
    bx.scatter(per_sat["gain"], y, s=30, color=S.INK, marker="|", zorder=3,
               label="GR's gain for its orbit (no offset)")
    bx.set_yticks(y)
    bx.set_yticklabels(per_sat["prn"], fontsize=6.5)
    bx.set_xlim(-3, 42)
    bx.set_ylim(-1, len(y) + 3.5)
    bx.set_xlabel("clock rate relative to the ground (µs per day)")
    bx.set_title("The 32 GPS clocks, 2016–2026")
    bx.legend(loc="upper center", fontsize=7.5, ncol=1, bbox_to_anchor=(0.5, 1.0))
    bx.annotate(f"factory offset\n{summary['factory_offset_us_day']:.3f} µs/day", (19, len(y) * 0.42),
                fontsize=7.5, color=S.INK_SECONDARY, ha="center")
    bx.annotate("", (1.4, len(y) * 0.36), (37.2, len(y) * 0.36),
                arrowprops={"arrowstyle": "->", "color": S.MUTED, "lw": 0.9})
    S.caption(fig, "Left: (W0 − GM/r)/c² and −GM/(2rc²) for circular orbits, W0 the geoid potential. "
              f"Right: {summary['n_satellite_days']:,} GPS satellite-days from ESA's final products (first day of each month): "
              "each clock's median rate against the ground reference, and the gain GR predicts from its actual orbit (J2 included).",
              y=-0.04)
    return S.save(fig, FIG / "gps_clock_budget")


# ---------------------------------------------------------------------------
# The Sun
# ---------------------------------------------------------------------------

def solar(lines: pd.DataFrame, windows: pd.DataFrame, res: dict) -> list:
    import solar as SO
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11.5, 4.6), gridspec_kw={"width_ratios": [1.5, 1]})
    weak = lines[~((lines["ew_mA"] > 150) & (lines["ew_mA"] < 550))]
    strong = lines[(lines["ew_mA"] > 150) & (lines["ew_mA"] < 550)]
    ax.scatter(weak["depth"], weak["v_core"], s=5, color=S.CONTEXT, label=f"{len(weak):,} Fe I lines")
    ax.scatter(strong["depth"], strong["v_core"], s=14, color=S.BLUE, marker="o", edgecolor=S.SURFACE,
               lw=0.5, label=f"{len(strong)} strongest (equivalent width 150–550 mÅ)")
    b = SO.binned(lines)
    ax.errorbar(b["center"], b["median"], yerr=b["sem"], color=S.INK, lw=1.4, marker="s", ms=3.5,
                label="median in depth bins")
    ax.axhline(SO.V_GRS_EARTH, color=S.INK, ls="--", lw=1.2)
    ax.text(0.02, SO.V_GRS_EARTH + 18, f"general relativity: {SO.V_GRS_EARTH:.1f} m/s", fontsize=8, color=S.INK)
    ax.axhline(0, color=S.AXIS, lw=0.8)
    ax.annotate("", (0.2, 260), (0.2, SO.V_GRS_EARTH - 8),
                arrowprops={"arrowstyle": "->", "color": S.MUTED, "lw": 1})
    ax.text(0.215, 420, "convective\nblueshift", fontsize=7.5, color=S.INK_SECONDARY)
    ax.set_xlim(0, 1)
    ax.set_ylim(-150, 1000)
    ax.set_xlabel("line depth (0 = no line, 1 = black core)")
    ax.set_ylabel("line-core shift (m/s, + = redshift)")
    ax.set_title("Iron lines in sunlight: the third signature of granulation")
    ax.legend(loc="upper left", fontsize=7.5)

    bx.axhline(SO.V_GRS_EARTH, color=S.INK, ls="--", lw=1.2, label="general relativity (633.1 m/s)")
    bx.errorbar(windows["half_window_kms"], windows["mean"], yerr=windows["sem"], color=S.BLUE, marker="o",
                ms=5, lw=1.4, capsize=2,
                label=f"IAG atlas, strongest lines (this work; {res['n_strong']} after 4σ clip)")
    bx.errorbar([1.56], [res["gh_lines_iag"]], yerr=[res["gh_lines_iag_err"]], color=S.ORANGE, marker="s",
                ms=5, capsize=2, ls="none", label="IAG atlas, their 15 lines (this work)")
    bx.errorbar([1.66], [639], yerr=[14], color=S.MUTED, marker="D", ms=5, capsize=2, ls="none",
                label="González Hernández+2020, HARPS-LFC\n(15 lines; five-pixel core fit ≈ ±1.6 km/s)")
    bx.set_xlim(0.8, 2.2)
    bx.set_ylim(540, 700)
    bx.set_xlabel("half-width of the core fit (km/s)")
    bx.set_ylabel("mean core shift of the strongest lines (m/s)")
    bx.set_title("Measured redshift vs. how much core is fitted")
    bx.legend(loc="lower left", fontsize=7)
    S.caption(fig, "IAG solar flux atlas (Reiners et al. 2016) with Nave et al. (1994) laboratory wavelengths on the NIST scale "
              "(+6.7e-8). Left: parabola over ±1.6 km/s of each core; depths and equivalent widths from Reiners et al.'s Table A1. "
              "Right: wider fits reach the convectively blueshifted flanks.", y=-0.06)
    return S.save(fig, FIG / "solar_redshift")


# ---------------------------------------------------------------------------
# Galileo: the clocks themselves
# ---------------------------------------------------------------------------

def galileo_clocks(start: str = "2020-03-02", ndays: int = 3) -> list:
    """Three days of E14 and a nominal Galileo clock: published, with the relativity term added
    back, and what is left after subtracting general relativity's prediction."""
    import clockfit as F
    import prepare
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 8.2), sharex=True,
                             gridspec_kw={"height_ratios": [0.8, 1.6, 1.1]})
    for prn, style in (("E14", S.SAT_STYLE["E14"]), ("E26", S.SAT_STYLE["control"])):
        d = prepare.load(prn)
        idx = list(d["date"]).index(start)
        for k in range(ndays):
            i = idx + k
            hours = (np.arange(289) * 300.0) / 3600 + 24 * k
            y, tpl = F.series_for(d, i, "full")
            f = F.fit_day(y, {}, keep_resid=True)
            phys = d["clk"][i] + d["kc"][i]
            gr = d["mono"][i] + d["j2"][i] + d["tide"][i] + d["kin"][i]
            ok = f.mask
            X = np.vstack([np.ones(ok.sum()), F.TC[ok]]).T
            detr = lambda z: z[ok] - X @ np.linalg.lstsq(X, z[ok], rcond=None)[0]
            lab = (lambda s: s if k == 0 else None)
            axes[0].plot(hours, (d["r"][i] - 6_371e3) / 1e3, color=style["color"], lw=1.6,
                         label=lab(f"{prn}" + (" (eccentric)" if prn == "E14" else " (nominal orbit)")))
            axes[1].plot(hours[ok], detr(phys) * 1e9, color=style["color"], lw=0 if prn == "E14" else 1.4,
                         marker="o" if prn == "E14" else None, ms=2.2,
                         label=lab(f"{prn} clock, relativity term added back"))
            if prn == "E14":
                axes[1].plot(hours[ok], detr(gr) * 1e9, color=S.INK, lw=1.2,
                             label=lab("general relativity's prediction"))
            axes[2].plot(hours[ok], f.resid[ok] * 1e12, color=style["color"], lw=1.0,
                         label=lab(f"{prn}: clock minus prediction"))
    axes[0].set_ylabel("height (km)")
    axes[0].set_title("Watch a satellite clock feel gravity: E14, 2–4 March 2020")
    axes[0].legend(loc="lower right", fontsize=7.5, ncol=2)
    axes[1].set_ylabel("clock offset (ns)\n(daily offset and rate removed)")
    axes[1].set_ylim(-600, 640)
    axes[1].legend(loc="upper center", fontsize=7.5, ncol=3)
    axes[2].set_ylabel("residual (ps)")
    axes[2].legend(loc="lower left", fontsize=7.5, ncol=2)
    axes[2].set_xlabel("hours from 2 March 2020, 00:00 GPS time")
    axes[2].set_xlim(0, 24 * ndays)
    for ax in axes:
        for k in range(1, ndays):
            ax.axvline(24 * k, color=S.AXIS, lw=0.8)
    S.annotate(axes[1], "deep and fast through perigee, the clock loses time; a quarter orbit later it is\n"
               "0.4 µs behind, and high and slow around apogee it catches up", (15.6, -440), (22.5, -560),
               va="bottom")
    S.caption(fig, "ESA/ESOC final products (5-min orbits and clocks). Published clocks follow the IGS convention: "
              "the once-per-orbit term −2r·v/c² is left for the user to add, which we do. Each day is a separate "
              "solution, so offset and rate are removed day by day. Note the different units: ns above, ps below.",
              y=-0.02)
    return S.save(fig, FIG / "galileo_clocks")


def galileo_fold(res: dict) -> list:
    import redshift_test as R
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    for prn in ("E14",):
        f = R.folded(prn)
        x = np.degrees(f["centers"])
        ax.plot(x, f["phys"] * 1e9, color=S.SAT_STYLE[prn]["color"], marker="o", ms=3, lw=0,
                label=f"{prn} clock, {int(f['n'].max()):,} days folded")
        ax.plot(x, f["tpl"] * 1e9 * 2, color=S.INK, lw=1.2, label="general relativity (gravity + speed)")
        ax.plot(x, f["tpl"] * 1e9, color=S.GRAVITY, lw=1.0, ls="--", label="its gravity half (what α scales)")
        bx.fill_between(x, (f["resid"] - f["resid_sem"]) * 1e12, (f["resid"] + f["resid_sem"]) * 1e12,
                        color=S.BLUE_150, lw=0, alpha=0.8)
        bx.plot(x, f["resid"] * 1e12, color=S.SAT_STYLE[prn]["color"], marker="o", ms=3, lw=1,
                label=f"{prn}: clock minus general relativity (±1 s.e. band)")
        for a, ls, lab in ((1.4e-4, ":", "a violation α = ±1.4×10⁻⁴ (GP-A's limit)"),
                           (2.5e-5, "--", "α = ±2.5×10⁻⁵ (Delva et al. 2018's error)")):
            bx.plot(x, a * f["tpl"] * 1e12, color=S.INK_SECONDARY, lw=1, ls=ls, label=lab)
            bx.plot(x, -a * f["tpl"] * 1e12, color=S.INK_SECONDARY, lw=1, ls=ls)
    for a in (ax, bx):
        a.set_xlim(-180, 180)
        a.set_xticks(np.arange(-180, 181, 90))
        a.set_xlabel("orbital phase from perigee (mean anomaly, degrees)")
    ax.set_ylabel("clock offset (ns)")
    ax.set_title("Folded over 7.4 years: the swing")
    ax.legend(loc="lower left", fontsize=7.5)
    bx.set_ylabel("clock minus GR (ps)")
    bx.set_title("…and what is left of it")
    bx.set_ylim(-40, 40)
    bx.legend(loc="upper right", fontsize=7, ncol=1)
    S.caption(fig, "E14 (GSAT-0202) passive hydrogen maser, Aug 2016 – Dec 2023: each day's clock (offset and rate "
              "removed) averaged in 5° bins of mean anomaly, every day weighted once. No Sun-locked signals are "
              "modelled here; the global fit's nuisance terms remove part of what remains.", y=-0.04)
    return S.save(fig, FIG / "galileo_fold")


def galileo_years(res: dict) -> list:
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    for prn, dx in (("E14", -0.12), ("E18", 0.12)):
        st = S.SAT_STYLE[prn]
        yr = pd.DataFrame(res["per_satellite"][prn]["yearly"])
        lo, hi = res["segments"][prn]
        used = (yr["year"] >= int(lo[:4])) & (yr["year"] <= int(hi[:4]))
        for mask, alpha_, lab in ((used, 1.0, f"{st['label']}, one year at a time"), (~used, 0.3, None)):
            ax.errorbar(yr["year"][mask] + dx, yr["alpha"][mask] * 1e5, yerr=yr["err_boot"][mask] * 1e5,
                        color=st["color"], marker=st["marker"], ms=5, ls="none", capsize=2, alpha=alpha_,
                        label=lab)
        if "alpha_orbit" in yr:
            m = yr["alpha_orbit"].notna()
            ax.errorbar(yr["year"][m] + dx * 2.2, yr["alpha_orbit"][m] * 1e5, yerr=yr["alpha_orbit_err"][m] * 1e5,
                        color=S.CONTEXT, marker=st["marker"], mfc="none", ms=4.5, ls="none", capsize=1.5,
                        label="what laser-ranged orbit errors alone would give (○ E14, □ E18)" if prn == "E14" else None)
    c = res["combined"]
    ax.axhspan((c["alpha"] - c["total"]) * 1e5, (c["alpha"] + c["total"]) * 1e5, color=S.BLUE_150, alpha=0.5, lw=0,
               label=f"this work, all years with Sun-locked signals modelled: ({c['alpha'] * 1e5:+.1f} ± {c['total'] * 1e5:.1f})×10⁻⁵")
    ax.axhline(0, color=S.INK, lw=1)
    ax.errorbar([2016.0], [6.75], yerr=[5.62], color=S.MUTED, marker="D", ms=5, capsize=2, ls="none",
                label="Delva et al. 2018, GSAT-0202 = E14 (2015–17)")
    ax.errorbar([2017.0], [4.5], yerr=[3.1], color=S.MUTED, marker="^", ms=6, capsize=2, ls="none",
                label="Herrmann et al. 2018, E14's maser (2016–17), α for gravity alone")
    ax.set_ylim(-60, 60)
    ax.set_xlabel("year")
    ax.set_ylabel("α: deviation from GR's redshift (×10⁻⁵)")
    ax.set_title("Year by year, each eccentric clock on its own")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), fontsize=7.5, ncol=2)
    S.caption(fig, "Points: daily fits combined per calendar year (no Sun-locked signals modelled), block-bootstrap errors; "
              "faded points are clock periods excluded from the test (E18's non-nominal maser in 2016-17, E14's degrading "
              "maser from 2024). Hollow grey: the same template fitted to ESA's SLR pass residuals. E18's 2016 value "
              "(−288, its non-nominal maser) lies beyond the plotted range.", y=-0.2)
    return S.save(fig, FIG / "galileo_alpha_by_year")


def galileo_sun_vs_perigee() -> list:
    import redshift_test as R
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("years", [S.BLUE_150, S.BLUE_550])
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.4), sharex=True, sharey=True)
    for row, prn in enumerate(("E18", "E14")):
        df = R.once_per_rev(prn)
        lo, hi = R.SEGMENTS[prn]
        df = df[(df["date"] >= lo) & (df["date"] <= hi) & (df["q"] < 2e-24)]
        yrs = pd.to_datetime(df["date"]).dt.year + pd.to_datetime(df["date"]).dt.dayofyear / 365.25
        for col, (key, title) in enumerate((("phase_sun_deg", "phase measured from the Sun"),
                                            ("phase_perigee_deg", "phase measured from perigee"))):
            ax = axes[row, col]
            big = df["amp_ps"] > np.percentile(df["amp_ps"], 50)
            sc = ax.scatter(df["beta_deg"][big], df[key][big], c=yrs[big], cmap=cmap, s=4, lw=0,
                            vmin=2016.5, vmax=2026.8)
            amp = f" (median amplitude {df['amp_ps'].median():.0f} ps)" if col == 0 else ""
            ax.set_title(f"{prn}: {title}{amp}", fontsize=10)
            ax.set_ylim(-180, 180)
            ax.set_yticks(np.arange(-180, 181, 90))
            if col == 0:
                ax.set_ylabel("phase of the once-per-orbit\nclock signal (degrees)")
            if row == 1:
                ax.set_xlabel("Sun's elevation above the orbit plane, β (degrees)")
    cb = fig.colorbar(sc, ax=axes, shrink=0.6, pad=0.02)
    cb.set_label("year")
    S.caption(fig, "Each dot is one day (the larger half of the daily once-per-orbit signals). E18's signal keeps a fixed "
              "phase relative to the Sun (it flips around 2020-21) while its phase from perigee wanders: a thermal or "
              "radiation-pressure effect, not gravity, which would be locked to perigee. E14's daily signals are mostly "
              "noise.", y=-0.01)
    return S.save(fig, FIG / "galileo_sun_vs_perigee")


def galileo_history() -> list:
    import redshift_test as R
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(10.5, 6.2), sharex=True)
    for prn in ("E14", "E18"):
        st = S.SAT_STYLE[prn]
        h = R.clock_history(prn)
        t = pd.to_datetime(h["date"])
        ax.plot(t, h["rate_us_day"], ls="none", marker=st["marker"], ms=1.6, color=st["color"], label=st["label"])
        bx.semilogy(t, np.sqrt(h["q"] * 300.0) * 1e12, ls="none", marker=st["marker"], ms=1.6, color=st["color"],
                    label=st["label"])
        lo, hi = R.SEGMENTS[prn]
        y0 = 0.93 if prn == "E14" else 0.86
        bx.annotate("", (pd.Timestamp(lo), y0), (pd.Timestamp(hi), y0), xycoords=("data", "axes fraction"),
                    arrowprops={"arrowstyle": "<->", "color": st["color"], "lw": 1.2})
        bx.text(pd.Timestamp(lo), y0 + 0.015, f"{prn} days used", transform=bx.get_xaxis_transform(), fontsize=7.5,
                color=S.INK_SECONDARY)
    ax.set_ylabel("clock rate vs. ground (µs/day)")
    ax.set_ylim(-45, 50)
    ax.set_title("Ten years of two clocks: frequency re-settings and an ageing maser")
    ax.legend(loc="upper right", markerscale=4, fontsize=7.5)
    bx.set_ylabel("random-walk step per 5 min (ps)")
    bx.set_ylim(3, 3e4)
    bx.set_xlabel("date")
    S.caption(fig, "Top: each day's mean clock rate against ESA's ground reference (it includes the ~40 µs/day relativistic "
              "gain; the steps are the operators' frequency settings). Bottom: the clock noise per 5-minute step from each "
              "day's fit. E14's maser is quiet until 2024, then its noise grows several hundredfold; the 2025-26 levels look "
              "like a rubidium clock.", y=-0.02)
    return S.save(fig, FIG / "galileo_clock_history")


def galileo_budget_and_j2(res: dict) -> list:
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11.5, 4.6), gridspec_kw={"width_ratios": [1.25, 1]})
    rows = [("Gravity Probe A (1976 rocket)", 0.0, 14.0, S.MUTED, "^"),
            ("Herrmann et al. 2018 (E14, gravity part)", 4.5, 3.1, S.MUTED, "^"),
            ("Delva et al. 2018 (E14 + E18)", 0.19, 2.48, S.MUTED, "D")]
    b14, b18, c = res["per_satellite"]["E14"]["budget"], res["per_satellite"]["E18"]["budget"], res["combined"]
    rows += [("this work: E18, 2017–2026", b18["alpha"] * 1e5, b18["total"] * 1e5, S.ORANGE, "s"),
             ("this work: E14, 2016–2023", b14["alpha"] * 1e5, b14["total"] * 1e5, S.BLUE, "o"),
             ("this work: combined", c["alpha"] * 1e5, c["total"] * 1e5, S.INK, "*")]
    for k, (name, a, e, col, mk) in enumerate(rows):
        ax.errorbar([a], [k], xerr=[e], color=col, marker=mk, ms=7 if mk != "*" else 11, capsize=3, lw=1.6)
        ax.text(-23.5, k + 0.18, name, fontsize=7.8, color=S.INK, va="bottom")
    ax.text(-23.5, -0.95, "GP-A measured the total shift; its α limit is ±14×10⁻⁵", fontsize=7, color=S.MUTED)
    ax.axvline(0, color=S.INK, lw=1)
    ax.set_yticks([])
    ax.set_ylim(-1.2, len(rows) - 0.3)
    ax.set_xlim(-24, 24)
    ax.set_xlabel("α: fractional deviation from GR's redshift (×10⁻⁵, 1σ)")
    ax.set_title("Is the redshift exactly as Einstein said?")
    j2 = pd.DataFrame(res["j2_tide"])
    j2 = j2[j2["beyond_err"] < 0.05].reset_index(drop=True)   # E11: only 42 usable days
    yy = np.arange(len(j2))
    for k, r in j2.iterrows():
        st = S.SAT_STYLE.get(r["prn"], S.SAT_STYLE["control"])
        bx.errorbar([r["beyond"]], [k], xerr=[r["beyond_err"]], color=st["color"], marker=st["marker"], ms=6,
                    capsize=2, lw=1.4)
    bx.axvline(1, color=S.INK, lw=1.2, ls="--")
    bx.axvline(0, color=S.AXIS, lw=0.8)
    bx.set_yticks(yy)
    bx.set_yticklabels(j2["prn"])
    bx.set_ylim(-1.3, len(yy) - 0.5)
    bx.set_xlim(-0.2, 1.8)
    bx.set_xlabel("size of the beyond-Kepler term in the clocks (GR = 1)")
    bx.set_title("Every clock feels the Earth's flattening")
    bx.text(1.03, -0.85, "general relativity", fontsize=7.5, color=S.INK_SECONDARY, va="center")
    bx.text(0.03, -0.85, "no effect", fontsize=7.5, color=S.INK_SECONDARY, va="center")
    S.caption(fig, "Left: 1σ totals including systematics (ours: clock noise calibrated by null tests, nuisance-model "
              "choice, laser-ranging orbit check, and Delva et al.'s magnetic and thermal limits). Right: the twice-per-orbit "
              "relativistic signal (110-150 ps peak to peak) the published clocks keep because the IGS convention removes only "
              "the Keplerian term: J2 bending the orbit, J2's own potential and the tides, fitted in every clock.", y=-0.04)
    return S.save(fig, FIG / "galileo_result")
