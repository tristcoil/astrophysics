"""Every figure of the project (PNG + SVG in figures/).

    python figures.py            # all figures whose inputs exist
    python figures.py proxima    # one figure by name
"""

from __future__ import annotations

import json
import sys

import matplotlib.pyplot as plt
import numpy as np

import config as C
import plotstyle as PS

PS.use_style()
FIG = C.FIGURES


def _save(fig, name):
    return PS.save(fig, str(FIG / name))


# ---------------------------------------------------------------------------
# 1. Architecture of triples: the Multiple Star Catalog
# ---------------------------------------------------------------------------

SHOWCASE = {
    # wds: (label, inner, outer, label position in the zoom panel (log P_in, log P_out))
    "21044+5114": ("TIC 290061484", "Aa,Ab", "A,B", (0.02, 1.27)),
    "04007+1229": ("λ Tau", "Aa,Ab", "Aab,Ac", (0.85, 1.42)),
    "19499+4107": ("KOI-126", "Aa,Ab", "Ac,Aab", (-0.15, 1.78)),
    "19171+4116": ("HD 181068", "Ab1,Ab2", "Aa,Ab", (-0.42, 1.95)),
    "03082+4057": ("Algol", "Aa,Ab", "A,B", (0.15, 3.05)),
    "22047+4626": ("SS Lac", "Aa1,Aa2", "Aa,Ab", (1.20, 3.05)),
    "10246-1906": ("HS Hya", "Aa,Ab", "Aab,Ac", (-0.15, 2.55)),
    "17569-3445": ("V907 Sco", "Aa1,Aa2", "Aa,Ab", (0.95, 1.85)),
}


def fig_architecture():
    import architecture as A
    comp, sys_, orb, notes, pairs = A.msc_pairs()
    st = A.stability_table(pairs)
    fig, axs = plt.subplots(1, 3, figsize=(15.5, 5.4), gridspec_kw=dict(width_ratios=[1.15, 1, 1]))
    ax, zx, bx = axs
    both = pairs["inner_from_orbit"] & pairs["outer_from_orbit"]
    ok = np.isfinite(pairs["P_in"]) & np.isfinite(pairs["P_out"])
    est = pairs[ok & ~both]
    orbs = pairs[ok & both]
    for a in (ax, zx):
        a.scatter(np.log10(est["P_in"]), np.log10(est["P_out"]), s=5, color=PS.CONTEXT, lw=0,
                  label=f"a period estimated from a separation ({len(est):,})", rasterized=True)
        a.scatter(np.log10(orbs["P_in"]), np.log10(orbs["P_out"]), s=9 if a is ax else 16, color=PS.BLUE,
                  lw=0, alpha=0.85, label=f"both periods from orbits ({len(orbs):,})")
    x = np.linspace(-1.5, 9.8, 50)
    r0 = float(C.ma01_period_ratio(2.0, 1.0, 0.0))
    r6 = float(C.ma01_period_ratio(2.0, 1.0, 0.6))
    for a, xt in ((ax, 6.2), (zx, 0.9)):
        a.plot(x, x, color=PS.INK_SECONDARY, lw=0.8)
        a.plot(x, x + np.log10(r0), color=PS.INK, lw=1.0)
        a.plot(x, x + np.log10(r6), color=PS.INK, lw=1.0)
    ax.text(6.25, 5.85, "P_out = P_in", color=PS.INK_SECONDARY, fontsize=7.5, rotation=36)
    ax.text(5.05, 5.95, "stability limit,\ne_out = 0", color=PS.INK, fontsize=7.5, rotation=36)
    ax.text(3.15, 5.25, "e_out = 0.6", color=PS.INK, fontsize=7.5, rotation=36)
    ac = pairs[(pairs["wds"] == "14396-6050") & (pairs["inner"] == "A,B")]
    xx, yy = np.log10(ac["P_in"].iloc[0]), np.log10(ac["P_out"].iloc[0])
    ax.scatter([xx], [yy], s=46, facecolor="none", edgecolor=PS.INK, lw=1.2, zorder=5)
    PS.annotate(ax, "α Cen + Proxima", (xx, yy), (xx - 2.6, yy + 0.45), relpos=(1, 0.5))
    # the zoom box
    zr = (-0.5, 1.45, 1.15, 3.2)
    ax.plot([zr[0], zr[1], zr[1], zr[0], zr[0]], [zr[2], zr[2], zr[3], zr[3], zr[2]], color=PS.INK, lw=0.8)
    PS.annotate(ax, "compact triples\n(next panel)", (zr[1], 0.5 * (zr[2] + zr[3])), (2.0, 1.3), relpos=(0, 0.5))
    ax.set_xlim(-1.0, 7.5)
    ax.set_ylim(0.5, 9.6)
    ax.set_xlabel("inner period, log10 (days)")
    ax.set_ylabel("outer period, log10 (days)")
    ax.set_title("Every triple in the Multiple Star Catalog")
    ax.legend(loc="lower right", fontsize=7.5)
    for wds, (lab, inn, out, pos) in SHOWCASE.items():
        row = pairs[(pairs["wds"] == wds) & (pairs["inner"] == inn) & (pairs["outer"] == out)]
        if row.empty:
            continue
        xx, yy = np.log10(row["P_in"].iloc[0]), np.log10(row["P_out"].iloc[0])
        zx.scatter([xx], [yy], s=60, facecolor="none", edgecolor=PS.INK, lw=1.3, zorder=5)
        PS.annotate(zx, lab, (xx, yy), pos, relpos=(0.5, 0.5), fontsize=8, color=PS.INK)
    zx.text(1.02, 1.17, "P_out = P_in", color=PS.INK_SECONDARY, fontsize=7.5, rotation=31)
    zx.text(0.98, 1.58, "limit, e_out = 0", color=PS.INK, fontsize=7.5, rotation=31)
    zx.text(-0.33, 1.32, "limit, e_out = 0.6", color=PS.INK, fontsize=7.5, rotation=31)
    zx.set_xlim(zr[0], zr[1])
    zx.set_ylim(zr[2], zr[3])
    zx.set_xlabel("inner period, log10 (days)")
    zx.set_title("The compact end: periods of days and weeks")
    # right: margin above the limit
    so = st[st["both_orbits"]]
    m = np.log10(so["margin"])
    bins = np.arange(-1.5, 4.6, 0.1)
    bx.hist(m, bins=bins, color=PS.BLUE, edgecolor=PS.SURFACE, lw=0.6)
    bx.axvline(0, color=PS.INK, lw=1.0)
    n_below = int((so["margin"] < 1).sum())
    n_ecc = int(((so["margin"] < 1) & (so["e_out"] > 0.7)).sum())
    bx.text(0.03, 0.70, f"{n_below} of {len(so)} sit below the limit;\n{n_ecc} of them have outer e > 0.7\n"
                        "(long, poorly covered orbits)", transform=bx.transAxes, va="top", fontsize=7.5,
            color=PS.INK_SECONDARY)
    bx.text(0.02, 0.12, "unstable\nside", transform=bx.transAxes, fontsize=7.5, color=PS.INK_SECONDARY)
    tic = so[(so["wds"] == "21044+5114") & (so["inner"] == "Aa,Ab")]
    if not tic.empty:
        v = float(np.log10(tic["margin"].iloc[0]))
        top = bx.get_ylim()[1]
        bx.set_ylim(0, top * 1.28)
        bx.plot([v, v], [0, top * 1.08], color=PS.INK, lw=0.8)
        bx.text(v + 0.05, top * 1.08, f"TIC 290061484 ({tic['margin'].iloc[0]:.1f}× the limit)",
                fontsize=7.5, ha="left", va="bottom")
    bx.set_xlabel("period ratio / critical ratio, log10")
    bx.set_ylabel("couples with both orbits known")
    bx.set_title("How far from falling apart?")
    PS.caption(fig, "Tokovinin's Multiple Star Catalog (release 2026-06-19): every inner/outer pair of "
                    "subsystems. Lines: Mardling & Aarseth (2001) limit for equal masses. Right: each couple's "
                    "own masses and outer eccentricity (orbits assumed coplanar).", y=-0.03)
    fig.tight_layout()
    return _save(fig, "triples_architecture")


# ---------------------------------------------------------------------------
# 2. Orbit alignment
# ---------------------------------------------------------------------------

def fig_alignment():
    import architecture as A
    t17 = A.tokovinin2017_reproduction()
    comp, sys_, orb, notes, pairs = A.msc_pairs()
    v26 = A.msc2026_two_orbits(pairs, orb)
    fig, axs = plt.subplots(1, 3, figsize=(14.5, 4.6))
    ax = axs[0]
    run = np.array(t17["running"])
    ax.fill_between(run[:, 0], run[:, 1] - run[:, 2], run[:, 1] + run[:, 2], color=PS.BLUE, alpha=0.18, lw=0)
    ax.plot(run[:, 0], run[:, 1], color=PS.BLUE, lw=1.6)
    s = np.logspace(0.5, 4.3, 50)
    ax.plot(s, 1.31 - 0.41 * np.log10(s), color=PS.VIOLET, lw=1.2)
    PS.annotate(ax, "Tokovinin (2017):\nC = 1.31 − 0.41 log s", (300, 1.31 - 0.41 * np.log10(300)), (500, 0.75),
                color=PS.VIOLET, relpos=(0, 0.5))
    PS.annotate(ax, "this work: running mean\nof 24 triples", (run[-1, 0], run[-1, 1]),
                (2.5e3, 0.05), color=PS.BLUE, relpos=(0, 0.5))
    ax.axhline(0, color=PS.AXIS, lw=0.8)
    ax.set_xscale("log")
    ax.set_ylim(-0.6, 1.2)
    ax.set_xlabel("projected outer separation s (au)")
    ax.set_ylabel("sign correlation C  (1 = all co-rotating)")
    ax.set_title("Tight triples spin the same way")
    ax = axs[1]
    labels = []
    for k, d in enumerate(t17["by_mass"][:3]):
        ax.errorbar([k], [d["C"]], [d["C_err"]], fmt="o", color=PS.INK, ms=6, capsize=0, lw=1.2)
    pub = [(0.61, 0.10), (0.23, 0.11), (0.18, 0.12)]
    for k, (c, e) in enumerate(pub):
        ax.errorbar([k + 0.15], [c], [e], fmt="s", color=PS.VIOLET, ms=5, capsize=0, lw=1.2,
                    markerfacecolor=PS.SURFACE)
    ax.set_xticks([0.07, 1.07, 2.07], ["< 1 Msun", "1-2 Msun", "> 2 Msun"])
    ax.set_xlim(-0.5, 2.6)
    ax.set_ylim(-0.1, 0.9)
    ax.text(1.6, 0.8, "● this work", color=PS.INK, fontsize=8)
    ax.text(1.6, 0.73, "■ Tokovinin 2017", color=PS.VIOLET, fontsize=8)
    ax.set_ylabel("sign correlation C")
    ax.set_xlabel("mass of the inner primary")
    ax.set_title("…and low-mass ones the most")
    ax = axs[2]
    for phis, col, lab in ((np.array(t17["two_orbit_phis"]).ravel(), PS.INK, f"Tokovinin 2017, {t17['n_two_orbits']} triples"),
                           (v26[["phi1", "phi2"]].to_numpy().ravel(), PS.BLUE, f"MSC 2026, {len(v26)} triples")):
        xs = np.sort(phis)
        ax.plot(xs, np.arange(1, len(xs) + 1) / len(xs), color=col, lw=1.6, drawstyle="steps-post")
        labels.append((xs, col, lab))
    rng = np.random.default_rng(1)
    sim_r = A.simulate_phi_pairs(0.0, 180, 200000, rng).ravel()
    sim_a = A.simulate_phi_pairs(0.8, 70, 200000, rng).ravel()
    for xs, col, lab in ((np.sort(sim_r), PS.MUTED, "random orientations"),
                         (np.sort(sim_a), PS.VIOLET, "80% aligned within 70° (Tokovinin)")):
        ax.plot(xs, np.arange(1, len(xs) + 1) / len(xs), color=col, lw=1.2)
    ax.text(100, 0.33, "random orientations", color=PS.MUTED, fontsize=8)
    ax.text(75, 0.83, "80% within 70°\n(Tokovinin's model)", color=PS.VIOLET, fontsize=8)
    ax.text(5, 0.88, labels[0][2], color=PS.INK, fontsize=8)
    ax.text(5, 0.80, labels[1][2], color=PS.BLUE, fontsize=8)
    ax.set_xlim(0, 180)
    ax.set_ylim(0, 1)
    ax.set_xlabel("angle between the two orbits, Φ (deg; both node solutions)")
    ax.set_ylabel("cumulative fraction")
    ax.set_title("Mutual inclinations from two visual orbits")
    PS.caption(fig, "Left, middle: Tokovinin (2017) Table 1, 216 triples whose outer pair's sense of rotation is "
                    "known, recomputed. Right: triples with both orbits resolved; each gives two possible Φ (the "
                    "node of a visual orbit is ambiguous), both plotted.", y=-0.03)
    fig.tight_layout()
    return _save(fig, "orbit_alignment")


# ---------------------------------------------------------------------------
# 3. Tertiary companions of close binaries
# ---------------------------------------------------------------------------

def fig_tertiary():
    import architecture as A
    t06 = A.tokovinin2006_reproduction()
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    for q in t06["quartiles"]:
        xm = np.sqrt(q["P_lo"] * q["P_hi"])
        ax.errorbar([xm], [q["f_raw"]], [q["f_err"]], xerr=[[xm - q["P_lo"]], [q["P_hi"] - xm]], fmt="o",
                    color=PS.INK, ms=6, lw=1.2, capsize=0)
    # published: raw and corrected (Table 6), 96% below 2.9 d
    ax.errorbar([np.sqrt(1 * 7)], [0.80], [0.06], xerr=[[np.sqrt(7) - 1], [7 - np.sqrt(7)]], fmt="s",
                color=PS.VIOLET, ms=5, markerfacecolor=PS.SURFACE, lw=1.2)
    ax.errorbar([np.sqrt(7 * 30)], [0.40], [0.06], xerr=[[np.sqrt(210) - 7], [30 - np.sqrt(210)]], fmt="s",
                color=PS.VIOLET, ms=5, markerfacecolor=PS.SURFACE, lw=1.2)
    ax.plot([1.0, 2.9], [0.96, 0.96], color=PS.VIOLET, lw=1.6)
    PS.annotate(ax, "96% below 2.9 d\n(corrected for missed companions)", (2.0, 0.96), (3.6, 1.0),
                color=PS.VIOLET, relpos=(0, 0.5))
    PS.annotate(ax, "corrected, Tokovinin et al. 2006", (np.sqrt(210), 0.40), (9, 0.18), color=PS.VIOLET,
                relpos=(0, 0.5))
    PS.annotate(ax, "observed fraction, this work\n(165 binaries, in quarters)", (1.8, 0.79), (1.05, 0.52),
                color=PS.INK, relpos=(0, 0.5))
    ax.set_xscale("log")
    ax.set_xlim(0.9, 35)
    ax.set_ylim(0, 1.1)
    ax.set_xlabel("period of the close binary (days)")
    ax.set_ylabel("fraction with a third star")
    ax.set_title("The closest binaries almost always have a third star")
    PS.caption(fig, "Tokovinin, Thomas, Sterzik & Udry (2006, A&A 450, 681): 165 solar-type spectroscopic "
                    "binaries searched for tertiaries. Our counts from their Table 5; their maximum-likelihood "
                    "correction for undetected companions raises the fractions (violet).", y=-0.04)
    fig.tight_layout()
    return _save(fig, "tertiary_fraction")


# ---------------------------------------------------------------------------
# 8. Proxima's orbit, to scale
# ---------------------------------------------------------------------------

def fig_proxima():
    import proxima as X
    k17 = X.kervella2017(n=4000)
    upd = X.update_gaia_dr3(n=4000)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(13, 5.4), gridspec_kw=dict(width_ratios=[1, 1.25]))
    # left: the orbit on the sky, in gnomonic (tangent-plane) coordinates
    # centred on alpha Cen A+B; each point is scaled by its own distance
    d_ab_au = 1000 / 750.81 * 206264.806

    def sky_track(el, k, n=400):
        a, e, i, Om, w = (np.atleast_1d(el[q])[k] for q in ("a", "e", "i", "Omega", "omega"))
        E = np.linspace(0, 2 * np.pi, n)
        nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
        r = a * (1 - e * np.cos(E))
        u = nu + np.radians(w)
        Om, i = np.radians(Om), np.radians(i)
        north = r * (np.cos(u) * np.cos(Om) - np.sin(u) * np.sin(Om) * np.cos(i))
        east = r * (np.cos(u) * np.sin(Om) + np.sin(u) * np.cos(Om) * np.cos(i))
        away = r * np.sin(u) * np.sin(i)
        dist = d_ab_au + away
        return np.degrees(east / dist), np.degrees(north / dist)
    rng = np.random.default_rng(3)
    for k in rng.choice(len(upd["mc"]["a"]), 60, replace=False):
        if upd["mc"]["a"][k] > 0:
            e_, n_ = sky_track(upd["mc"], k)
            ax.plot(e_, n_, color=PS.BLUE, lw=0.4, alpha=0.35)
    e_u, n_u = sky_track(upd["nominal"], 0)
    ax.plot(e_u, n_u, color=PS.BLUE, lw=1.6)
    e_, n_ = sky_track(k17["nominal"], 0)
    ax.plot(e_, n_, color=PS.VIOLET, lw=1.2)
    ax.scatter([0], [0], s=60, color=PS.YELLOW, edgecolor=PS.INK, lw=0.8, zorder=5)
    PS.annotate(ax, "α Cen A+B", (0, 0), (1.3, -0.5), relpos=(0, 0.5))
    # Proxima today (Gaia DR3), standard coordinates about A+B (Akeson et al., J2016.0)
    import pandas as pd
    g3 = pd.read_csv(C.RAW / "gaia_dr3_proxima.csv").iloc[0]
    a21 = C.SYSTEMS["alpha Cen"]["a21"]
    ra0 = np.radians(a21["ra_deg"] + a21["pmra"] * (2016.0 - 2019.5) / 3.6e6 / np.cos(np.radians(a21["dec_deg"])))
    de0 = np.radians(a21["dec_deg"] + a21["pmdec"] * (2016.0 - 2019.5) / 3.6e6)
    ra1, de1 = np.radians(g3["ra"]), np.radians(g3["dec"])
    D = np.sin(de0) * np.sin(de1) + np.cos(de0) * np.cos(de1) * np.cos(ra1 - ra0)
    dE = np.degrees(np.cos(de1) * np.sin(ra1 - ra0) / D)
    dN = np.degrees((np.cos(de0) * np.sin(de1) - np.sin(de0) * np.cos(de1) * np.cos(ra1 - ra0)) / D)
    ax.scatter([dE], [dN], s=40, color=PS.RED, edgecolor=PS.INK, lw=0.6, zorder=6)
    PS.annotate(ax, "Proxima today", (dE, dN), (dE - 0.15, dN - 0.55), relpos=(0, 0.5))
    k_top = int(np.argmax(n_u))
    PS.annotate(ax, "orbit from Gaia DR3 + Akeson et al.\n2021 (this work); thin: 60 draws",
                (e_u[k_top], n_u[k_top]), (1.4, 1.85), color=PS.BLUE, relpos=(0, 0.5))
    ax.text(0.03, 0.04, "violet: from Kervella et al. (2017)'s inputs", transform=ax.transAxes,
            color=PS.VIOLET, fontsize=8)
    ax.set_aspect("equal")
    ax.set_xlim(1.5, -2.1)          # east to the left, as on the sky
    ax.set_ylim(-2.75, 2.15)
    ax.set_xlabel("east of α Cen (degrees on the sky)")
    ax.set_ylabel("north of α Cen (degrees)")
    ax.set_title("Proxima's orbit on our sky")
    # right: powers of ten
    s = summarize_proxima(upd)
    rows = [
        ("Earth–Sun", 1.0, 1.0, PS.INK),
        ("Neptune", 30.1, 30.1, PS.INK),
        ("Kuiper belt", 30, 50, PS.MUTED),
        ("α Cen A–B orbit", 0.4805 * 23.3, 1.5195 * 23.3, PS.YELLOW),
        ("Voyager 1 (2026)", 170, 170, PS.INK),
        ("Proxima's orbit", s["peri"], s["apo"], PS.BLUE),
        ("Oort cloud (inner to outer)", 2000, 100000, PS.MUTED),
        ("distance Sun – α Cen", 274_000, 274_000, PS.INK),
    ]
    for k, (lab, lo, hi, col) in enumerate(rows):
        y = len(rows) - k
        if hi > lo:
            bx.plot([lo, hi], [y, y], color=col, lw=7, solid_capstyle="butt")
        else:
            bx.scatter([lo], [y], s=40, color=col, zorder=4)
        bx.text(hi * 1.6, y, lab, va="center", fontsize=8.5, color=PS.INK)
    bx.set_xscale("log")
    bx.set_xlim(0.5, 3e6)
    bx.set_ylim(0.3, len(rows) + 0.7)
    bx.set_yticks([])
    bx.spines["left"].set_visible(False)
    bx.set_xlabel("distance from the central star(s) (au, logarithmic)")
    bx.set_title("Powers of ten: Proxima stays beyond ~4,000 au")
    PS.caption(fig, "Left: Keplerian orbit of Proxima relative to the A+B barycentre from 3-D positions and "
                    "velocities (radial velocities corrected for gravitational redshift, as in Kervella et al. "
                    "2017); thin lines: 60 Monte Carlo draws. The orbit is too slow (0.28 km/s) to follow on the "
                    "sky; the arc shown spans ~500,000 years. Right: Proxima's periastron-apastron range "
                    "against the Solar System's scales.", y=-0.03)
    fig.tight_layout()
    return _save(fig, "proxima_orbit")


def summarize_proxima(res):
    mc = res["mc"]
    b = mc["bound"]
    return {k: float(np.median(mc[k][b])) for k in ("a", "e", "peri", "apo", "P_yr")}


# ---------------------------------------------------------------------------
# 5. Triply eclipsing triples: a gallery
# ---------------------------------------------------------------------------

def _kepler(name):
    import glob
    from astropy.io import fits
    t, f = [], []
    for path in sorted(glob.glob(str(C.RAW / "kepler" / name / "*.fits"))):
        d = fits.open(path)[1].data
        ok = (d["SAP_QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"])
        ff = d["PDCSAP_FLUX"][ok]
        t.append(d["TIME"][ok] + 2454833.0)
        f.append(ff / np.nanmedian(ff))
    t, f = np.concatenate(t), np.concatenate(f)
    o = np.argsort(t)
    return t[o], f[o]


def fig_gallery():
    import etv
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
    # HD 181068: the pair crossing the red giant (Kepler, Q1-Q6)
    t, f = _kepler("hd_181068")
    tc = 2454833.0 + 190.12
    sel = np.abs(t - tc) < 3.0
    axs[0].plot(t[sel] - tc, f[sel], ".", ms=3, color=PS.INK)
    axs[0].set_title("HD 181068: two red dwarfs cross a red giant")
    axs[0].set_xlabel("days from mid-event (BJD 2455023.1)")
    axs[0].set_ylabel("relative brightness (Kepler)")
    axs[0].annotate("the red-dwarf pair passes in front\nof the giant for ~1.5 days", xy=(-0.8, 0.9765),
                    xytext=(0.30, 0.88), textcoords="axes fraction", fontsize=8, color=PS.INK_SECONDARY,
                    arrowprops=dict(arrowstyle="-", color=PS.INK_SECONDARY, lw=0.8))
    # KOI-126: the M-dwarf pair transiting the F star
    t, f = _kepler("koi_126")
    tc = 2454833.0 + 337.48
    sel = np.abs(t - tc) < 0.8
    axs[1].plot(t[sel] - tc, f[sel], ".", ms=3, color=PS.INK)
    axs[1].set_title("KOI-126: two M dwarfs transit an F star")
    axs[1].set_xlabel("days from mid-event (BJD 2455170.5)")
    axs[1].annotate("one M dwarf in front,\nthen both at once", xy=(0.0, 0.9727), xytext=(0.58, 0.25),
                    textcoords="axes fraction", fontsize=8, color=PS.INK_SECONDARY,
                    arrowprops=dict(arrowstyle="-", color=PS.INK_SECONDARY, lw=0.8))
    # TIC 290061484: our TESS photometry around a third-body eclipse
    lcs = etv.light_curves()
    lc = [x for x in lcs if x["sector"] == 75][0]
    tc = 2460344.80
    sel = np.abs(lc["jd"] - tc) < 3.0
    axs[2].plot(lc["jd"][sel] - tc, lc["flux"][sel], ".", ms=2, color=PS.INK, alpha=0.7)
    axs[2].set_title("TIC 290061484: three B stars, all eclipsing")
    axs[2].set_xlabel("days from BJD 2460344.80 (TESS sector 75)")
    axs[2].annotate("extra dip: the third star\nand the pair line up", xy=(0.02, 0.918), xytext=(0.56, 0.08),
                    textcoords="axes fraction", fontsize=8, color=PS.INK_SECONDARY,
                    arrowprops=dict(arrowstyle="-", color=PS.INK_SECONDARY, lw=0.8))
    axs[2].annotate("the 1.8-day pair's\nown eclipses", xy=(-2.22, 0.86), xytext=(0.03, 0.08),
                    textcoords="axes fraction", fontsize=8, color=PS.INK_SECONDARY,
                    arrowprops=dict(arrowstyle="-", color=PS.INK_SECONDARY, lw=0.8))
    PS.caption(fig, "Kepler long-cadence PDCSAP photometry (HD 181068, KOI-126; Derekas et al. 2011, "
                    "Carter et al. 2011) and our aperture photometry of TESS full-frame images "
                    "(TIC 290061484; Kostov et al. 2024).", y=-0.04)
    fig.tight_layout()
    return _save(fig, "triply_eclipsing_gallery")


# ---------------------------------------------------------------------------
# 6. Kozai-Lidov cycles in Algol, and what stops them
# ---------------------------------------------------------------------------

def fig_kozai():
    import kozai as K
    ts = K.timescales()
    runs = {}
    for lab in ("newton", "gr", "tides", "tides_Bsoft", "tides_Aonly"):
        path = K.OUT / f"algol_{lab}.json"
        if path.exists():
            runs[lab] = json.loads(path.read_text())
    a = C.SYSTEMS["Algol"]
    e_touch = 1 - (a["R_A"] + a["R_B"]) * C.R_SUN_AU / ts["a_in_au"]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw=dict(width_ratios=[1.45, 1]))
    col = {"newton": PS.RED, "gr": PS.YELLOW, "tides": PS.BLUE, "tides_Bsoft": PS.VIOLET,
           "tides_Aonly": PS.MAGENTA}
    for lab, r in runs.items():
        ax.plot(r["t_yr"], r["e_in"], color=col[lab], lw=1.8 if lab in ("newton", "gr", "tides") else 1.3,
                zorder=3 if lab != "tides" else 2)
        if r["touched_yr"]:
            ax.scatter([r["t_yr"][-1]], [r["e_in"][-1]], s=40, color=col[lab], edgecolor=PS.INK, lw=0.6, zorder=5)
    if "newton" in runs:
        t_n = runs["newton"]["touched_yr"]
        t_g = runs.get("gr", {}).get("touched_yr")
        txt = f"point masses: the stars touch after {t_n:.0f} yr"
        if t_g:
            txt += f"\nadding general relativity: after {t_g:.0f} yr\n(the red and yellow curves overlap)"
        ax.annotate(txt, xy=(t_n, e_touch), xytext=(t_n + 160, 0.40), fontsize=8, va="top",
                    arrowprops=dict(arrowstyle="-", color=PS.INK_SECONDARY, lw=0.8))
    if "tides_Aonly" in runs:
        r = runs["tides_Aonly"]
        k = int(np.argmax(r["e_in"]))
        ax.annotate(f"bulge of A only (B rigid): e peaks at {r['e_max']:.2f};\nno contact in 3,000 yr",
                    xy=(r["t_yr"][k], r["e_in"][k]), xytext=(r["t_yr"][k] + 260, r["e_in"][k] - 0.02),
                    fontsize=8, va="top", arrowprops=dict(arrowstyle="-", color=col["tides_Aonly"], lw=0.9))
    if "tides" in runs:
        emax = max(runs[lab]["e_max"] for lab in ("tides", "tides_Bsoft") if lab in runs)
        ax.text(2980, 0.035, "tidal and rotational bulges of both stars: e stays below "
                f"{emax:.3f}\n(blue; violet: B's bulge 10x weaker, just as flat)", fontsize=8, ha="right",
                va="bottom")
    ax.axhline(e_touch, color=PS.INK_SECONDARY, lw=0.8)
    ax.text(3050, e_touch + 0.012, f"periastron = R_A + R_B (e = {e_touch:.2f}): the stars collide",
            fontsize=7.5, color=PS.INK_SECONDARY, ha="right")
    ax.set_xlim(0, 3100)
    ax.set_ylim(-0.01, 0.62)
    ax.set_xlabel("years from now (simulation)")
    ax.set_ylabel("eccentricity of the inner orbit (A-B, 2.87 d)")
    ax.set_title("Algol's perpendicular orbits: what saves the eclipsing pair?")
    # right: the competing clocks
    rows = [("Kozai-Lidov cycle", ts["t_KL"], PS.RED), ("GR periastron advance", ts["t_GR"], PS.YELLOW),
            ("bulges of A and B", ts["t_tide"], PS.BLUE),
            ("bulges, B's 10x weaker", K.timescales(K.VARIANTS["tides_Bsoft"])["t_tide"], PS.VIOLET),
            ("bulge of A only", K.timescales(K.VARIANTS["tides_Aonly"])["t_tide"], PS.MAGENTA)]
    ys = np.arange(len(rows))[::-1]
    bx.barh(ys, [r[1] for r in rows], color=[r[2] for r in rows], height=0.55)
    for y, (_, v, _) in zip(ys, rows):
        bx.text(v * 1.25, y, f"{v:,.0f} yr" if v > 10 else f"{v:.1f} yr", va="center", fontsize=8.5)
    bx.set_yticks(ys, [r[0] for r in rows])
    bx.axvline(ts["t_KL"], color=PS.RED, lw=0.8, ls=(0, (3, 3)))
    bx.set_xscale("log")
    bx.set_xlim(0.5, 3e4)
    bx.set_xlabel("timescale (years, logarithmic)")
    bx.set_title("The race: precession faster than Kozai wins")
    PS.caption(fig, "Labelled simulations (REBOUND IAS15; REBOUNDx gr_full and tides_spin, conservative) "
                    "started from Baron et al. (2012)'s orbits: masses 3.17, 0.70, 1.76 Msun, mutual "
                    "inclination 90.2 deg. Bulges: apsidal constants 0.005 (A) and 0.05 (B, the "
                    "Roche-lobe-filling subgiant), synchronous spins. Bars: Kozai-Lidov timescale "
                    "and periastron-precession periods, standard formulae.", y=-0.04)
    fig.tight_layout()
    return _save(fig, "kozai_algol")


# ---------------------------------------------------------------------------
# 4. Eclipses that switch off and on
# ---------------------------------------------------------------------------

SWITCH_SYSTEMS = ("SS Lac", "HS Hya", "V907 Sco")
SWITCH_YEARS = {"SS Lac": (1880, 2300), "HS Hya": (1880, 2300), "V907 Sco": (1880, 2090)}


def posterior_draws(name, n=200, seed=0):
    """Equal-weight posterior draws: the final file if present, else the live
    points of the running nested sampler (for previews)."""
    import switch_fit as F
    path = F.posterior_path(name, F.FINAL_TAG[name])
    rng = np.random.default_rng(seed)
    if path.exists():
        post = np.load(path, allow_pickle=True)["posterior"]
        return post[rng.choice(len(post), min(n, len(post)), replace=False)], True
    import __main__
    import dynesty
    __main__.GlobalFit = F.GlobalFit
    smp = dynesty.NestedSampler.restore(str(F.CACHE / f"dynesty_{F._slug(name)}.save"))
    fit = F.GlobalFit(name)
    live = np.array([fit.prior_transform(u) for u in smp.live_u])
    return live[rng.choice(len(live), min(n, len(live)), replace=False)], False


def folded_tilt(i_deg):
    """Tilt folded to 0-90 deg (90 = edge-on): eclipses cannot tell i from 180 - i."""
    return 90 - np.abs(90 - np.asarray(i_deg))


def published_tilt(name, years):
    """Published models of the tilt (deg, folded), for comparison."""
    yrs = np.asarray(years, float)
    if name == "SS Lac":
        pub = C.SYSTEMS[name]["pub"]
        I, eps, t90, Pn = np.radians(84.0), np.radians(15.0), 1905.8, 600.0
        Om90 = np.arccos(1 / np.tan(I) / np.tan(eps))
        cosi = np.cos(I) * np.cos(eps) - np.sin(I) * np.sin(eps) * np.cos(Om90 + 2 * np.pi / Pn * (yrs - t90))
        return {"Torres 2001 (nodal model)": folded_tilt(np.degrees(np.arccos(cosi)))}
    if name == "HS Hya":
        # Davenport et al. 2021: linear tilt change through i = 71.3 deg in 2021.1,
        # one 'tumble' per 1,194 yr
        i_lin = 71.3 + 360.0 / 1194.0 * (2021.1 - yrs)
        return {"Davenport et al. 2021 (linear)": folded_tilt(i_lin)}
    if name == "V907 Sco":
        # Zasche et al. 2023 secular part: J = 26.2, I = 71, P = 63.5 yr,
        # anchored at i = 76.61 deg (TESS S39b, 2021.43) and increasing
        import switching as S
        sp = S.SPECS[name]
        eps, _, _, _, _, _ = S.secular(sp.m1, sp.m2, 1.062, sp.P, sp.P_out, 0.0, sp.e_out, np.radians(26.2))
        I = np.radians(71.0)
        ci = np.cos(np.radians(76.61))
        ch = (np.cos(I) * np.cos(eps) - ci) / (np.sin(I) * np.sin(eps))
        h0 = -np.arccos(np.clip(ch, -1, 1))          # increasing tilt: h moving toward 0
        h = h0 + 2 * np.pi / 63.5 * (yrs - 2021.43)
        cosi = np.cos(I) * np.cos(eps) - np.sin(I) * np.sin(eps) * np.cos(h)
        return {"Zasche et al. 2023 (secular part)": folded_tilt(np.degrees(np.arccos(cosi)))}
    return {}


def fig_switching():
    import switch_fit as F
    import switching as S
    fig, axs = plt.subplots(3, 2, figsize=(14, 12.5), gridspec_kw=dict(width_ratios=[2.2, 1]))
    for row, name in enumerate(SWITCH_SYSTEMS):
        ax, bx = axs[row]
        col = PS.SYSTEM_COLOR[name]
        draws, final = posterior_draws(name)
        m = S.Model(name)
        y0, y1 = SWITCH_YEARS[name]
        yrs = np.linspace(y0, y1, 600)
        jd = C.year_to_jd(yrs)
        tilts = []
        rss = []
        for th in draws:
            g = m.geometry(m.unpack(th))
            tilts.append(folded_tilt(np.degrees(m.inclination(g, jd))))
            rss.append(g["r1"] + g["r2"])
        tilts = np.array(tilts)
        rs = float(np.median(rss))
        i_edge = np.degrees(np.arccos(rs))          # eclipses need tilt above this (circular orbit)
        ax.axhspan(i_edge, 90, color=col, alpha=0.10, lw=0)
        ax.text(y1 - 2, 89.2, "eclipses happen above this line", ha="right", va="top", fontsize=7.5,
                color=PS.INK_SECONDARY)
        ax.axhline(i_edge, color=col, lw=0.8)
        q = np.percentile(tilts, [16, 50, 84], axis=0)
        ax.fill_between(yrs, q[0], q[2], color=col, alpha=0.25, lw=0)
        ax.plot(yrs, q[1], color=col, lw=1.8)
        # the best retrograde solution (a fit limited to J > 90 deg), if one was run
        alt = F.alternative_solution(name) if final else None
        if alt:
            lnB = alt["lnB"]
            rpost = np.load(F.posterior_path(name, alt["tag"]), allow_pickle=True)["posterior"]
            rsel = rpost[np.random.default_rng(1).choice(len(rpost), min(200, len(rpost)), replace=False)]
            rt = np.median([folded_tilt(np.degrees(m.inclination(m.geometry(m.unpack(th)), jd))) for th in rsel],
                           axis=0)
            ax.plot(yrs, rt, color=col, lw=1.2, ls=(0, (4, 3)))
            r2 = float(f"{np.exp(abs(lnB)):.2g}")          # two significant figures
            odds = f"odds 1:{r2:,.0f}" if lnB > 0 else f"odds {r2:,.0f}:1"
            vis = np.where((rt > 64) & (rt < 87) & (yrs > y0 + 0.62 * (y1 - y0)))[0]
            if len(vis):
                k = int(vis[len(vis) // 3])
                sense = "retrograde" if "retro" in alt["tag"] else "prograde"
                ax.text(yrs[k], rt[k] + 1.0, f"dashed: best {sense} solution ({odds})", color=PS.INK_SECONDARY,
                        fontsize=7.5, ha="center", va="bottom", zorder=6,
                        bbox=dict(facecolor=PS.SURFACE, edgecolor="none", alpha=0.85, pad=1.5))
        for lab, curve in published_tilt(name, yrs).items():
            ax.plot(yrs, curve, color=PS.VIOLET, lw=1.1)
            # label where the curve is on the plot, late in the range
            vis = np.where((curve > 63) & (curve < 87) & (yrs > y0 + 0.45 * (y1 - y0)))[0]
            k = int(vis[len(vis) // 2]) if len(vis) else int(0.3 * len(yrs))
            ax.text(yrs[k], min(88, max(62, curve[k] + 1.0)), lab, color=PS.VIOLET, fontsize=7.5,
                    ha="center", va="bottom", clip_on=True)
        # plates, era by era (global fit's best draw)
        best = draws[0]
        for e in F.era_inclinations(name, best):
            yc = e["year"]          # the era's tilt refers to its median plate date
            if e["dlnl_vs_none"] > 2.0:
                ax.errorbar([yc], [e["i_best"]], [[e["i_best"] - e["i_lo"]], [e["i_hi"] - e["i_best"]]],
                            xerr=[[yc - e["lo"]], [e["hi"] - yc]], fmt="o", color=PS.INK, ms=4, lw=1.0)
            else:
                ax.errorbar([yc], [e["i_hi95"]], xerr=[[yc - e["lo"]], [e["hi"] - yc]], fmt="v", color=PS.INK,
                            ms=5, lw=1.0, mfc=PS.SURFACE)
        if name == "V907 Sco":
            # Zasche et al. 2023, Table 2: their own tilt measurements
            for mjd, inc, err, src in C.SYSTEMS[name]["pub"]["inclinations"]:
                yy = float(C.jd_to_year(mjd + 2400000.0))
                ax.errorbar([yy + 0.6], [inc], [err], fmt="D", color=PS.VIOLET, ms=3.5, lw=0.9,
                            mfc=PS.SURFACE)
            ax.text(1962, 88.6, "◇ Zasche et al. 2023, Table 2", color=PS.VIOLET, fontsize=7.5)
        # TESS
        for pr in F.profiles(name):
            yt = C.jd_to_year(pr["jd_mid"])
            x = np.array(pr["x"])
            c = (np.array(pr["chi2"]) - pr["chi2_min"]) / pr["scale"]
            ok = c < 1.0
            ti = lambda xx: np.degrees(np.arccos(np.clip(np.asarray(xx) * rs, 0, 1)))
            if pr["detected"]:
                lo, hi = ti(x[ok].max()), ti(x[ok].min())
                ax.errorbar([yt], [ti(pr["x_best"])], [[ti(pr["x_best"]) - lo], [hi - ti(pr["x_best"])]],
                            fmt="s", color=PS.INK, ms=4, lw=1.0)
            else:
                # no eclipse: tilt below the value where chi2 rises by 4
                lim = ti(x[c < 4.0].min())
                ax.errorbar([yt], [lim], fmt="v", color=PS.INK, ms=5, mfc=PS.SURFACE)
        # Hipparcos and ASAS-3 (in the final fit for HS Hya; checks for the others)
        gpath = F.CACHE / f"ground_profiles_{F._slug(name)}.json"
        for pr in (json.loads(gpath.read_text()) if gpath.exists() else []):
            yt = C.jd_to_year(pr["jd_mid"])
            x = np.array(pr["x"])
            c = (np.array(pr["chi2"]) - pr["chi2_min"]) / pr["scale"]
            ok = c < 1.0
            ti = lambda xx: np.degrees(np.arccos(np.clip(np.asarray(xx) * rs, 0, 1)))
            if pr["detected"] or x[ok].max() < x.max() - 0.05:
                lo, hi = ti(x[ok].max()), ti(x[ok].min())
                ax.errorbar([yt], [ti(pr["x_best"])], [[ti(pr["x_best"]) - lo], [hi - ti(pr["x_best"])]],
                            fmt="D", color=PS.INK, ms=3.6, lw=0.9, mfc=PS.SURFACE)
            else:
                ax.errorbar([yt], [ti(x[c < 4.0].min())], fmt="v", color=PS.INK, ms=4, mfc=PS.SURFACE, alpha=0.7)
        ax.axvline(2026.8, color=PS.INK_SECONDARY, lw=0.8)
        ax.text(2028, 61.3, "now", fontsize=7.5, color=PS.INK_SECONDARY)
        ax.set_xlim(y0, y1)
        ax.set_ylim(60, 90.5)
        ax.set_ylabel("tilt of the eclipsing pair's orbit\n(90° = edge-on)")
        ax.set_title(f"{name}: " + {"SS Lac": "eclipses stopped around 1950",
                                    "HS Hya": "eclipses appeared ~1920, stopped by 2021",
                                    "V907 Sco": "on, off, and on again since 2019"}[name])
        if row == 2:
            ax.set_xlabel("year")
        if not final:
            ax.text(0.01, 0.03, "PREVIEW: sampler still running", transform=ax.transAxes, color=PS.RED)
        # right: folded plates, eclipsing era vs silent era
        d = S.dasch_data(name)
        g = m.geometry(m.unpack(best))
        fit = F.GlobalFit(name)
        eras = {"SS Lac": ((1890, 1920), (1955, 1990)), "HS Hya": ((1890, 1918), (1940, 1956)),
                "V907 Sco": ((1900, 1920), (1925, 1945))}[name]
        for k, (lo, hi) in enumerate(eras):
            sel = (d["year"] >= lo) & (d["year"] < hi)
            t = d["date_jd"].to_numpy()[sel]
            # phase from the global model's clock: time since the nearest primary conjunction
            tc = np.array([fit.conjunction_time(g, tt) for tt in t])
            P = 2 * np.pi / g["n_lam"]
            ph = ((t - tc) / P + 0.25) % 1.0 - 0.25
            off = 1.4 * k
            mag = d["mag"].to_numpy()[sel] - np.median(d["mag"]) + off
            bx.plot(ph, mag, ".", ms=2.0, color=PS.CONTEXT)
            edges = np.linspace(-0.25, 0.75, 41)
            idx = np.digitize(ph, edges) - 1
            med = [np.median(mag[idx == j]) if np.sum(idx == j) >= 4 else np.nan for j in range(40)]
            bx.plot(0.5 * (edges[1:] + edges[:-1]), med, color=col, lw=1.6)
            bx.text(-0.24, off - 0.45, f"{lo}-{hi} ({sel.sum()} plates)", fontsize=7.5, ha="left")
        bx.set_ylim(2.4, -0.75)
        bx.set_xlim(-0.25, 0.75)
        bx.set_xlabel("orbital phase (primary eclipse at 0)" if row == 2 else "")
        bx.set_ylabel("plate magnitude (offset)")
        bx.set_title("Harvard plates, folded")
    PS.caption(fig, "Tilt of each eclipsing pair's orbit to the line of sight. Shaded band: our fit (68% of "
                    "posterior draws; quadrupole secular precession about the total angular momentum). Circles: "
                    "plate eras at their median plate date (DASCH; bars: era span; triangles: no eclipse, upper "
                    "limit); squares: TESS sectors; open diamonds: "
                    "Hipparcos (1990-93) and ASAS-3 seasons (used in HS Hya's fit, checks for the others). Violet: "
                    "published models. Right: plates folded on the fitted eclipse clock; line: binned median.",
               y=0.015)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.96, bottom=0.07, hspace=0.32, wspace=0.18)
    return _save(fig, "eclipses_switching")


# ---------------------------------------------------------------------------
# 5. TIC 290061484: the outer orbit in the eclipse times
# ---------------------------------------------------------------------------

def fig_tic290_etv():
    import tic290_nbody as N
    obs = N.observed_times()
    eph = json.loads((C.DERIVED / "tic290061484" / "ephemeris.json").read_text())
    fitr = N.run_fit(variant="triple")
    best = dict(fitr["base"], **fitr["best"])
    obs = obs.assign(oc_min=(obs["t_obs"] - (eph["T0"] + obs["cycle"] * eph["P"])) * 1440)
    assert np.allclose(obs["t_obs"], fitr["t_obs"])
    obs = obs.assign(res_min=(obs["t_obs"] - np.array(fitr["t_model"])) * 1440)
    wins = [(obs.loc[obs["sector"].isin(secs), "t_obs"].min() - 0.6,
             obs.loc[obs["sector"].isin(secs), "t_obs"].max() + 0.6) for secs in N.SEASONS.values()]

    def oc_curve(params):
        tm, km = N.model_minima(params, wins)
        cyc = np.round((tm - eph["T0"]) / eph["P"] * 2) / 2
        return tm, km, (tm - (eph["T0"] + cyc * eph["P"])) * 1440
    tm, km, oc_m = oc_curve(best)
    tk, kk, oc_k = oc_curve(fitr["base"])
    scan = N.etv_period_scan(obs)
    per = N.outer_mean_periods(best)
    pub = C.SYSTEMS["TIC 290061484"]["pub"]
    fig = plt.figure(figsize=(13, 8.6))
    gs = fig.add_gridspec(3, 4, height_ratios=[2.2, 1.0, 1.5], hspace=0.42, wspace=0.08)
    kinds = (("pri", 1.0, PS.BLUE, "primary eclipses (Ab in front of Aa)"),
             ("sec", -1.0, PS.ORANGE, "secondary eclipses (Aa in front of Ab)"))
    lo = min(obs["oc_min"].min(), oc_m.min(), oc_k.min()) - 3
    hi = max(obs["oc_min"].max(), oc_m.max(), oc_k.max()) + 9
    for c, (season, secs) in enumerate(N.SEASONS.items()):
        ax = fig.add_subplot(gs[0, c])
        rx = fig.add_subplot(gs[1, c], sharex=ax)
        g = obs[obs["sector"].isin(secs)]
        t1, t2 = wins[c]
        x0 = np.floor(t1)
        for kind, kv, col, lab in kinds:
            sel = (km == kv) & (tm > t1) & (tm < t2)
            ax.plot(tm[sel] - x0, oc_m[sel], color=col, lw=1.3, zorder=2)
            sk = (kk == kv) & (tk > t1) & (tk < t2)
            ax.plot(tk[sk] - x0, oc_k[sk], color=col, lw=0.9, ls=(0, (2, 2)), alpha=0.9, zorder=1,
                    label=("Kostov et al. 2024 solution (inner period calibrated)" if (c == 0 and kind == "pri")
                           else None))
            h = g[g["kind"] == kind]
            ax.errorbar(h["t_obs"] - x0, h["oc_min"], yerr=h["err"] * 1440, fmt="o", ms=3.2, color=col,
                        mfc=PS.SURFACE if kind == "sec" else col, mew=0.9, elinewidth=0.7, capsize=0,
                        zorder=3, label=lab if c == 0 else None)
            rx.errorbar(h["t_obs"] - x0, h["res_min"], yerr=h["err"] * 1440, fmt="o", ms=2.6, color=col,
                        mfc=PS.SURFACE if kind == "sec" else col, mew=0.8, elinewidth=0.6, capsize=0)
        rx.axhline(0, color=PS.INK_SECONDARY, lw=0.8)
        ax.set_ylim(lo, hi)
        rx.set_ylim(-11, 11)
        y_, m_ = C.jd_to_year(np.array([t1]))[0], None
        ax.set_title(f"TESS sectors {secs[0]}-{secs[1]} ({season[:4]})", fontsize=10)
        rx.set_xlabel(f"days from BJD {x0:.0f}")
        if c == 0:
            ax.set_ylabel("O - C (minutes)")
            rx.set_ylabel("data - model\n(minutes)")
        else:
            ax.tick_params(labelleft=False)
            rx.tick_params(labelleft=False)
        ax.tick_params(labelbottom=False)
        ax.text(0.03, 0.04, f"rms {np.sqrt(np.mean(g['res_min'] ** 2)):.1f} min", transform=rx.transAxes,
                fontsize=7.5, color=PS.INK_SECONDARY)
    fig.axes[0].legend(loc="upper left", fontsize=7.2, frameon=False)
    # bottom: the period in the timing variations
    bx = fig.add_subplot(gs[2, :])
    d = (scan["chi2"] - scan["chi2"].min()) / scan["chi2_scale"]
    bx.plot(scan["periods"], d, color=PS.INK, lw=1.4)
    bx.axvspan(scan["P_lo"], scan["P_hi"], color=PS.CONTEXT, lw=0)
    bx.axvline(pub["P_out_mean"], color=PS.VIOLET, lw=1.0, ls=(0, (4, 3)))
    bx.text(pub["P_out_mean"] + 0.4, 0.82 * min(d.max(), 400), f"Kostov et al. 2024:\nmean {pub['P_out_mean']:.2f} d",
            color=PS.VIOLET, fontsize=8, va="top")
    bx.text(scan["P_best"] - 0.6, 0.82 * min(d.max(), 400),
            f"timing alone:\n{scan['P_best']:.2f} d ({scan['P_lo']:.2f}-{scan['P_hi']:.2f})", fontsize=8,
            ha="right", va="top")
    bx.set_ylim(-5, min(d.max(), 400) * 1.02)
    bx.set_xlim(scan["periods"].min(), scan["periods"].max())
    bx.set_xlabel("trial period of the timing wave (days)")
    bx.set_ylabel("Δχ² (rescaled)")
    bx.set_title("Which period makes the eclipses run early and late? (each season fitted separately: "
                 "trend + two harmonics)", fontsize=10, loc="left")
    fig.suptitle("TIC 290061484: the inner pair's eclipses keep time with a 24.5-day outer orbit",
                 x=0.07, ha="left", fontsize=13, fontweight="bold", y=0.985)
    PS.caption(fig, f"Our eclipse times from TESS full-frame images ({len(obs)} eclipses; template fits), relative "
                    f"to a linear ephemeris (P = {eph['P']:.5f} d); eclipses near third-star events left out. Solid "
                    f"lines: our four-body REBOUND fit to these times (fourth star held at Kostov et al.'s orbit; "
                    f"chi2 {fitr['chi2']:.0f} for {fitr['dof']} degrees of freedom; outer orbit "
                    f"{per['conjunction_B_in_front']['P']:.3f} d between conjunctions, e = {best['A_B.e']:.3f}). "
                    f"Dashed: Kostov et al.'s solution with only the inner period and phase adjusted "
                    f"(chi2 {fitr['chi2_start']:.0f}).", y=0.012)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.92, bottom=0.11)
    return _save(fig, "tic290061484_etv")


# ---------------------------------------------------------------------------
# 6. Stability: chaos maps and a triple pushed past the edge
# ---------------------------------------------------------------------------

def _edges(c):
    c = np.asarray(c, float)
    m = 0.5 * (c[1:] + c[:-1])
    return np.concatenate([[c[0] - (m[0] - c[0])], m, [c[-1] + (c[-1] - m[-1])]])


def fig_stability():
    import matplotlib.colors as mcolors
    import stability as St
    sysT = C.SYSTEMS["TIC 290061484"]
    q_out = sysT["masses"]["B"] / (sysT["masses"]["Aa"] + sysT["masses"]["Ab"])
    ratio_real = sysT["A_B"]["a_rsun"] / sysT["Aa_Ab"]["a_rsun"]
    e_real, i_real = sysT["A_B"]["e"], sysT["pub"]["i_mut"]
    cls_cols = [PS.BLUE_RAMP[0], PS.BLUE_RAMP[5], PS.BLUE_RAMP[12]]
    cmap = mcolors.ListedColormap(cls_cols)
    fig = plt.figure(figsize=(13, 8.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], hspace=0.38, wspace=0.2)
    maps = (("ratio_e", "outer eccentricity", np.linspace(0, 0.8, 200),
             lambda y: C.ma01_ratio(q_out, y, i_real), e_real,
             f"Coplanar-ish (mutual inclination {i_real:.0f}°): eccentric outer orbits need more room"),
            ("ratio_i", "mutual inclination (degrees)", np.linspace(0, 180, 200),
             lambda y: C.ma01_ratio(q_out, e_real, y), i_real,
             f"Outer eccentricity {e_real:.2f}: tilted orbits are wilder than the criterion says"))
    for col, (kind, ylabel, ygrid, crit, y_real, title) in enumerate(maps):
        d = np.load(St.OUT / f"megno_{kind}.npz")
        cls = np.where(d["broken"], 2, np.where(d["megno"] >= 2.5, 1, 0))
        ax = fig.add_subplot(gs[0, col])
        ax.pcolormesh(_edges(d["ratios"]), _edges(d["ys"]), cls, cmap=cmap, vmin=-0.5, vmax=2.5,
                      shading="flat", rasterized=True)
        xc = np.array([crit(y) for y in ygrid])
        ax.plot(xc, ygrid, color=PS.ORANGE, lw=2.0)
        ax.scatter([ratio_real], [y_real], marker="*", s=220, color=PS.YELLOW, edgecolor=PS.INK, lw=0.8, zorder=5)
        ax.annotate("TIC 290061484", xy=(ratio_real, y_real),
                    xytext=(ratio_real + 0.3, y_real - 0.075 if kind == "ratio_e" else y_real + 14),
                    fontsize=8.5, color=PS.INK, fontweight="bold")
        if kind == "ratio_e":
            ax.scatter([2.5], [e_real], marker="X", s=90, color=PS.SURFACE, edgecolor=PS.INK, lw=1.0, zorder=5)
            ax.annotate("start of the 400\nbreak-up runs (below)", xy=(2.5, e_real), xytext=(2.15, 0.52),
                        fontsize=7.5, color=PS.SURFACE, arrowprops=dict(arrowstyle="-", color=PS.SURFACE, lw=0.7))
            ax.text(8.95, 0.36, "Mardling & Aarseth\n(2001) limit", color=PS.INK, fontsize=8,
                    ha="right", va="top")
        else:
            ax.text(4.95, 16, "Mardling & Aarseth\n(2001) limit", color=PS.INK, fontsize=8,
                    ha="left", va="center")
            ax.text(8.95, 116, "near 90°, Kozai-Lidov cycles\npump the inner orbit's eccentricity", fontsize=8,
                    color=PS.INK, ha="right", va="bottom")
            ax.set_yticks([0, 30, 60, 90, 120, 150, 180])
        ax.set_xlim(d["ratios"][0] - 0.1, d["ratios"][-1] + 0.1)
        ax.set_ylim(d["ys"][0] - 0.5 * (d["ys"][1] - d["ys"][0]), d["ys"][-1] + 0.5 * (d["ys"][1] - d["ys"][0]))
        ax.set_xlabel("outer orbit size / inner orbit size  (a_out / a_in)")
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=10, loc="left")
        ax.grid(False)
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=c, edgecolor=PS.AXIS, label=l) for c, l in
               zip(cls_cols, ("regular (MEGNO ≈ 2)", "chaotic (MEGNO > 2.5), still intact",
                              "broke up within 300 outer orbits"))]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, 0.955))
    # bottom left: who gets thrown out
    r = json.loads((St.OUT / "disrupt_ratio2.50.json").read_text())
    res = r["results"]
    ej = [x for x in res if x["ejected"]]
    theo = St.ejection_theory()
    bx = fig.add_subplot(gs[1, 0])
    names = ["Aa", "Ab", "B"]
    m = sysT["masses"]
    xs = np.arange(3)
    fr = np.array([sum(1 for x in ej if x["ejected"] == n) / len(ej) for n in names])
    er = np.sqrt(fr * (1 - fr) / len(ej))
    bx.bar(xs, fr, width=0.55, color=[PS.STAR_COLOR[n] for n in names])
    bx.errorbar(xs, fr, yerr=er, fmt="none", ecolor=PS.INK, elinewidth=1.0, capsize=3)
    bx.scatter(xs, [theo[n] for n in names], marker="D", s=46, color=PS.SURFACE, edgecolor=PS.INK, lw=1.2,
               zorder=5, label="statistical theory: P ∝ m⁻³")
    for x, n, f in zip(xs, names, fr):
        bx.text(x + 0.31, f, f"{100 * f:.0f}%", va="center", fontsize=8.5)
    bx.set_xticks(xs, [f"{n}\n({m[n]:.2f} M☉)" for n in names])
    bx.set_ylabel("fraction of break-ups")
    bx.set_ylim(0, 0.62)
    bx.legend(loc="upper left", frameon=False, fontsize=8)
    bx.set_title(f"Pushed past the edge (a_out/a_in = 2.5): who is thrown out? ({len(ej)} of {len(res)} runs)",
                 fontsize=10, loc="left")
    # bottom right: how long it takes
    cx = fig.add_subplot(gs[1, 1])
    t = np.sort([x["t_out_orbits"] for x in ej])
    cx.step(t, np.arange(1, len(t) + 1) / len(res), where="post", color=PS.INK, lw=1.6)
    med = float(np.median(t))
    cx.axvline(med, color=PS.INK_SECONDARY, lw=0.8, ls=(0, (3, 3)))
    cx.text(med * 1.08, 0.30, f"half have broken up\nafter {med:.0f} outer orbits", fontsize=8)
    surv = len(res) - len(ej)
    cx.text(0.98, 0.03, f"{surv} of {len(res)} still triple after {r['t_max_out']} outer orbits", transform=cx.transAxes,
            fontsize=8, ha="right", color=PS.INK_SECONDARY)
    ab = np.median([x["a_bin"] for x in ej])
    eb = np.median([x["e_bin"] for x in ej])
    cx.text(0.02, 0.97, f"the pair left behind is tighter\n(median a = {ab:.2f} x the original)\n"
                        f"and eccentric (median e = {eb:.2f})", transform=cx.transAxes, fontsize=8, va="top")
    cx.set_xscale("log")
    cx.set_xlim(5, r["t_max_out"])
    cx.set_ylim(0, 1.0)
    cx.set_xlabel("time until one star escapes (outer orbits, log)")
    cx.set_ylabel("fraction of runs broken up")
    cx.set_title("How long a doomed triple lasts", fontsize=10, loc="left")
    fig.suptitle("When does a triple fall apart?", x=0.07, ha="left", fontsize=13, fontweight="bold", y=0.995)
    PS.caption(fig, "Labelled simulations (REBOUND IAS15) with TIC 290061484's masses (Aa 6.85, Ab 6.11, B 7.90 Msun). "
                    "Maps: one random orbital phase per cell, 300 outer orbits; MEGNO from REBOUND's variational "
                    "equations ('broke up': a star escaped, or the inner orbit changed size by 2x or reached e > 0.99). "
                    "Break-ups: 400 random phases at e_out = 0.20 and 9° mutual inclination, up to 2,000 outer orbits.",
               y=0.012)
    fig.subplots_adjust(left=0.065, right=0.985, top=0.885, bottom=0.1)
    return _save(fig, "stability")


FIGURES = {"architecture": fig_architecture, "alignment": fig_alignment, "tertiary": fig_tertiary,
           "proxima": fig_proxima, "gallery": fig_gallery, "kozai": fig_kozai,
           "switching": fig_switching, "tic290_etv": fig_tic290_etv,
           "stability": fig_stability}

if __name__ == "__main__":
    names = sys.argv[1:] or list(FIGURES)
    for n in names:
        print(n, FIGURES[n]())
