"""Every figure, drawn from the cached products of ``analysis.py``."""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D

import config as cfg
import galaxy
import outer
import plotstyle as ps
import reference
import survey
import terminal

ps.use_style()
FIG = cfg.FIGURE_DIR
P = cfg.PRODUCTS_DIR


def _load(name: str) -> dict:
    with np.load(P / f"{name}.npz", allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def _results() -> dict:
    return json.loads(cfg.RESULTS_PATH.read_text())


def _wrap(l_deg):
    """Longitude on (-180, 180], for plots centred on the Galactic centre."""
    return np.where(np.asarray(l_deg) > 180, np.asarray(l_deg) - 360, l_deg)


def _break_at_wrap(x, y, jump=90.0):
    """Insert NaNs where a wrapped longitude track jumps across the plot."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    cut = np.flatnonzero(np.abs(np.diff(x)) > jump) + 1
    return np.insert(x, cut, np.nan), np.insert(y, cut, np.nan)


def _lon_axis(ax, which="x"):
    ticks = np.arange(-180, 181, 30)
    labels = [f"{t % 360:d}°" for t in ticks]
    if which == "x":
        ax.set_xticks(ticks, labels)
        ax.set_xlim(180, -180)
    ax.set_xlabel("Galactic longitude l")


# --------------------------------------------------------------------------
# 1. The longitude-velocity diagram
# --------------------------------------------------------------------------

#: Arm labels in the longitude-velocity diagram: (text, wrapped l, v, leader
#: target or None). Placed where the diagram is empty (beyond the terminal
#: velocity, or between well-separated tracks) so no label hides data.
LV_LABELS = [
    ("Scutum", 36, 128, (31.6, 104.0)),
    ("Sagittarius", 53, 97, (48.5, 66.0)),
    ("Norma", -40, -140, (-32.3, -103.0)),
    ("Centaurus", -55, -70, (-53.7, -45.0)),
    ("Carina", -77, -38, (-74.0, -7.0)),
    ("Local", 87, 20, None),
    ("Perseus", -150, 34, None),
    ("Outer", -140, 67, None),
    ("Outer Scutum-Centaurus", 40, -118, (38.0, -70.0)),
]


def fig_lv():
    lv = survey.lv_diagram(0.0, 0.0)
    nhi = survey.column_density_map()
    mas = _load("masers")
    obs = galaxy.Observer.reid19()

    l = _wrap(lv["l"])
    order = np.argsort(l)
    l, T = l[order], lv["T"][order]
    v = lv["v"]
    lb = _wrap(nhi["l"])
    ob = np.argsort(lb)

    fig = plt.figure(figsize=(13, 7.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.2, 1], width_ratios=[1, 0.018],
                          hspace=0.09, wspace=0.015)
    ax0 = fig.add_subplot(gs[0, 0])
    cax0 = fig.add_subplot(gs[0, 1])
    ax = fig.add_subplot(gs[1, 0], sharex=ax0)
    cax = fig.add_subplot(gs[1, 1])

    # (a) the hydrogen sky
    N = nhi["NHI"][:, ob]
    im0 = ax0.imshow(N, origin="lower", extent=[lb[ob][0], lb[ob][-1], nhi["b"][0], nhi["b"][-1]],
                     cmap=ps.CMAP_GREY, norm=LogNorm(vmin=2e20, vmax=2e22), aspect="equal",
                     interpolation="nearest")
    ax0.set_xlim(180, -180)
    ax0.set_ylim(-11, 11)
    ax0.set_yticks([-10, 0, 10], ["−10°", "0°", "+10°"])
    ax0.set_ylabel("b")
    ax0.grid(False)
    ax0.tick_params(labelbottom=False)
    ax0.set_title("The Milky Way's hydrogen, seen from inside: where it is (top) and how it moves "
                  "(bottom)", loc="left")
    cb0 = fig.colorbar(im0, cax=cax0)
    cb0.set_label("N$_{HI}$ (cm$^{-2}$)", fontsize=7.5)
    cb0.ax.tick_params(labelsize=7)

    # (b) longitude-velocity at b = 0
    stretch = lambda x: np.arcsinh(np.clip(x, 0, None) / 4.0)
    im = ax.imshow(stretch(T.T), origin="lower", aspect="auto",
                   extent=[l[0], l[-1], v[0], v[-1]], cmap=ps.CMAP_GREY,
                   vmin=0, vmax=stretch(130.0), interpolation="nearest")
    cb = fig.colorbar(im, cax=cax)
    ticks = [0, 2, 5, 10, 20, 50, 100]
    cb.set_ticks([stretch(t) for t in ticks], labels=[str(t) for t in ticks])
    cb.set_label("brightness temperature T$_b$ at b = 0 (K)", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)

    dashed = dict(color=ps.C_ARM, lw=0.9, ls=(0, (3, 2)))
    for arm in galaxy.REID19_ARMS:
        la, va, d, ind = galaxy.arm_lv(arm, galaxy.theta_reid19, obs)
        x, y = _break_at_wrap(_wrap(la), va)
        xin, yin = _break_at_wrap(_wrap(np.where(ind, la, np.nan)), np.where(ind, va, np.nan))
        ax.plot(x, y, **dashed)
        ax.plot(xin, yin, color=ps.C_ARM, lw=1.5)
    sct = next(a for a in galaxy.REID19_ARMS if a.name == "Scutum-Centaurus")
    beta = np.linspace(*galaxy.OSC_BETA, 400)
    xo, yo = sct.xy(beta)
    lo_, do_ = galaxy.ld_from_xy(xo, yo, obs.r0)
    ax.plot(*_break_at_wrap(_wrap(lo_), galaxy.v_lsr(lo_, 0.0, do_, galaxy.theta_reid19, obs)),
            **dashed)
    for name, (v0, slope) in galaxy.THREE_KPC_ARMS.items():
        L = np.linspace(*galaxy.THREE_KPC_L_RANGE, 20)
        ax.plot(L, v0 + slope * L, color=ps.C_ARM, lw=1.5)
    box = dict(boxstyle="round,pad=0.15", fc=ps.SURFACE, ec="none", alpha=0.8)
    for text, lx, vy, target in LV_LABELS + [("near 3-kpc arm", 9, -140, (-6.0, -78.0)),
                                             ("far 3-kpc arm", -19, 38, (-9.0, 19.3))]:
        if target is None:
            ax.text(lx, vy, text, fontsize=8, fontweight="bold", ha="center", va="center",
                    color=ps.INK, bbox=box)
        else:
            ax.annotate(text, xy=target, xytext=(lx, vy), fontsize=8, fontweight="bold",
                        ha="center", va="center", color=ps.INK, bbox=box,
                        arrowprops=dict(arrowstyle="-", color=ps.INK_SECONDARY, lw=0.7,
                                        shrinkA=0, shrinkB=2))

    sel = mas["plx"] > 0
    ax.scatter(_wrap(mas["l"][sel]), mas["vlsr"][sel], s=9, color=ps.C_MASER, lw=0.4,
               edgecolor=ps.SURFACE, zorder=5)

    # annotations of the pattern itself
    kw = dict(color=ps.INK_SECONDARY)
    ps.annotate(ax, "terminal velocity:\nthe tangent points", (61, 40), (80, 150), **kw)
    ps.annotate(ax, "terminal velocity", (-24, -125), (-68, -185), **kw)
    ps.annotate(ax, "gas beyond the solar circle\n(outer Galaxy)", (140, -48), (140, -160), **kw)
    ps.annotate(ax, "gas beyond the solar circle", (-150, 55), (-120, 165), **kw)
    ps.annotate(ax, "Galactic centre:\nnuclear disc and bar", (1.5, 200), (32, 222), **kw)
    ax.set_ylim(-250, 250)
    ax.set_ylabel("line-of-sight velocity v$_{LSR}$ (km/s)")
    _lon_axis(ax)
    ax.grid(False)
    handles = [Line2D([], [], color=ps.C_ARM, lw=1.5, label="spiral arms: Reid et al. (2019) "
                      "maser model, projected with circular rotation (dashed: extrapolated)"),
               Line2D([], [], marker="o", ls="", color=ps.C_MASER, ms=4,
                      label="masers in young massive-star regions (Reid et al. 2019)")]
    ax.legend(handles=handles, loc="lower left", fontsize=7.5, frameon=True, facecolor=ps.SURFACE,
              edgecolor="none", framealpha=0.85)
    ps.caption(fig, "HI4PI 21-cm survey (HI4PI Collaboration 2016): 16' resolution, 1.3 km/s "
               "channels. Top: hydrogen column density |b| ≤ 11°, all velocities within ±350 km/s. "
               "Bottom: brightness at b = 0 (arcsinh stretch); the pale stripe at l = 0 is hydrogen "
               "absorbing the radio glow of Sgr A. The 3-kpc arms are drawn from their measured CO "
               "tracks (Dame & Thaddeus 2008), since they are not in circular rotation.")
    return ps.save(fig, FIG / "01_lv_diagram")


# --------------------------------------------------------------------------
# 2. How a terminal velocity is measured
# --------------------------------------------------------------------------

def fig_terminal_method():
    lv = survey.lv_diagram(-cfg.TERMINAL_B_HALFWIDTH_DEG, cfg.TERMINAL_B_HALFWIDTH_DEG)
    tv = _load("terminal")
    mgd = reference.mgd16_q1()
    v = lv["v"]
    examples = [30.0, 50.0, 320.0]
    fig = plt.figure(figsize=(12.5, 7.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15], hspace=0.38, wspace=0.28)
    for k, l0 in enumerate(examples):
        ax = fig.add_subplot(gs[0, k])
        i = int(np.argmin(np.abs(lv["l"] - l0)))
        sign = int(terminal.quadrant_sign(l0))
        fit = terminal.fit_edge(v, lv["T"][i], sign)
        spec = lv["T"][i]
        ax.plot(v, spec, color=ps.C_OBSERVED, lw=1.3, label="HI4PI, mean over |b| ≤ 0.5°")
        u_thr = sign * fit.v_threshold
        uu = np.linspace(u_thr - cfg.TERMINAL_W_IN_KMS, u_thr + cfg.TERMINAL_W_OUT_KMS, 300)
        model = terminal.two_erfc(uu, sign * fit.v_terminal, fit.amplitudes[0], fit.widths[0],
                                  fit.amplitudes[1], fit.widths[1])
        ax.plot(sign * uu, model, color=ps.ORANGE, lw=2.0, label="two-error-function edge fit")
        ax.axhline(cfg.TERMINAL_THRESHOLD_K, color=ps.MUTED, lw=0.8)
        ax.axvline(fit.v_terminal, color=ps.ORANGE, lw=0.9)
        top = np.nanmax(spec) * 1.15
        ax.set_ylim(-2, top)
        side = "left" if sign > 0 else "right"
        ax.text(fit.v_terminal + sign * 3, 0.8 * top, f"V$_t$ = {fit.v_terminal:+.1f} km/s",
                color=ps.INK, fontsize=8, ha=side, va="center")
        if sign > 0:
            ax.set_xlim(fit.v_terminal - 90, fit.v_terminal + 45)
        else:
            ax.set_xlim(fit.v_terminal - 45, fit.v_terminal + 90)
        ax.text(ax.get_xlim()[1 if sign > 0 else 0] - sign * 3, cfg.TERMINAL_THRESHOLD_K + 0.03 * top,
                "10 K threshold", color=ps.MUTED, fontsize=7.5,
                ha="right" if sign > 0 else "left", va="bottom")
        ax.set_xlabel("v$_{LSR}$ (km/s)")
        if k == 0:
            ax.set_ylabel("brightness temperature (K)")
            ax.legend(loc="lower left", fontsize=7)
        ax.set_title(f"l = {l0:.0f}°" + ("  (first quadrant)" if sign > 0 else "  (fourth quadrant)"),
                     fontsize=9.5)
        # MGD16 value for the Q1 examples
        if sign > 0:
            j = int(np.argmin(np.abs(mgd["l"] - l0)))
            ax.plot([mgd["v"][j]], [0.6 * top], marker="v", color=ps.C_PUBLISHED, ms=6)
            ax.text(mgd["v"][j] + 3, 0.6 * top, f"VLA survey: {mgd['v'][j]:+.1f}\n(MGD 2016)",
                    color=ps.C_PUBLISHED, fontsize=7, va="center", ha="left")

    # bottom: first-quadrant l-v zoom with the measured edge
    ax = fig.add_subplot(gs[1, :])
    sel = (lv["l"] >= 12) & (lv["l"] <= 78)
    stretch = lambda x: np.arcsinh(np.clip(x, 0, None) / 4.0)
    im = ax.imshow(stretch(lv["T"][sel].T), origin="lower", aspect="auto",
                   extent=[lv["l"][sel][0], lv["l"][sel][-1], v[0], v[-1]], cmap=ps.CMAP_GREY,
                   vmin=0, vmax=stretch(110), interpolation="nearest")
    ax.plot(mgd["l"], mgd["v"], "o", color=ps.C_PUBLISHED, ms=3.2, mew=0, zorder=3,
            label="McClure-Griffiths & Dickey 2016 (VLA + Green Bank, 1')")
    ax.plot(tv["Q1_l"], np.where(tv["Q1_reliable"], tv["Q1_v_terminal"], np.nan), color=ps.ORANGE,
            lw=1.1, zorder=4, label="terminal velocity, this project (HI4PI, 16')")
    ax.set_xlim(78, 12)
    ax.set_ylim(-40, 170)
    ax.set_xlabel("Galactic longitude l (deg)")
    ax.set_ylabel("v$_{LSR}$ (km/s)")
    ax.grid(False)
    ax.legend(loc="upper right", fontsize=7.5, frameon=True, facecolor=ps.SURFACE,
              edgecolor="none", framealpha=0.85)
    ax.set_title("First quadrant, |b| ≤ 0.5°: the emission's upper edge is the terminal "
                 "velocity", fontsize=9.5)
    ps.caption(fig, "The edge is located where the latitude-averaged brightness last exceeds "
               "10 K, then fitted over −25 to +45 km/s around it with two complementary error "
               "functions sharing one centre (McClure-Griffiths & Dickey 2007): the centre, where "
               "emission has fallen to half its level inside the edge, is the terminal velocity. "
               "Grey scale: arcsinh stretch of T_b.")
    return ps.save(fig, FIG / "02_measuring_terminal_velocity")


# --------------------------------------------------------------------------
# 3. Terminal velocities: both quadrants, against the literature
# --------------------------------------------------------------------------

def fig_terminal_curve():
    tv = _load("terminal")
    ns = _load("ns_asymmetry")
    res = _results()["terminal"]
    mgd = reference.mgd16_q1()
    fig = plt.figure(figsize=(12.5, 8.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], hspace=0.34, wspace=0.22)

    ax = fig.add_subplot(gs[0, :])
    for quad, color, name in (("Q1", ps.C_NORTH, "first quadrant (north), this project"),
                              ("Q4", ps.C_SOUTH, "fourth quadrant (south), this project")):
        l = tv[f"{quad}_l"]
        s = np.abs(np.sin(np.radians(l)))
        v_ok = np.where(tv[f"{quad}_reliable"], np.abs(tv[f"{quad}_v_terminal"]), np.nan)
        ax.plot(s, v_ok, color=color, lw=1.4, label=name)
    ms = np.abs(np.sin(np.radians(mgd["l"])))
    ax.plot(ms, mgd["v"], ".", color=ps.C_PUBLISHED, ms=2.5,
            label="first quadrant, McClure-Griffiths & Dickey 2016 (VLA)")
    s07 = np.abs(np.sin(np.radians(339.4)))
    ax.plot([s07], [135.72], marker="D", color=ps.C_PUBLISHED, ms=4, ls="",
            label="fourth quadrant, McClure-Griffiths & Dickey 2007 (10 published values, mean)")
    ax.set_xlabel("|sin l|   (= tangent-point radius R / R$_0$)")
    ax.set_ylabel("|terminal velocity| (km/s)")
    ax.set_xlim(0.25, 0.97)
    ax.set_ylim(0, 170)
    ax.legend(loc="upper right", fontsize=7.5)
    ax.set_title("Terminal velocities in both halves of the inner Galaxy", loc="left")
    q = res["q1_vs_mgd16"]
    ax.text(0.27, 12, f"first quadrant vs the VLA survey, {q['n']} longitudes: median offset "
            f"{q['median']:+.2f} km/s, rms {q['rms']:.2f} km/s", fontsize=7.5,
            color=ps.INK_SECONDARY)

    # residuals against MGD16
    ax = fig.add_subplot(gs[1, 0])
    ours = np.interp(mgd["l"], tv["Q1_l"], tv["Q1_v_terminal"])
    ax.axhline(0, color=ps.AXIS, lw=0.8)
    ax.plot(mgd["l"], ours - mgd["v"], ".", color=ps.C_NORTH, ms=2.5)
    ax.set_xlim(67.5, 18)
    ax.set_ylim(-4, 4)
    ax.set_xlabel("Galactic longitude l (deg)")
    ax.set_ylabel("this project − VLA survey (km/s)")
    ax.set_title("Point by point, first quadrant", fontsize=9.5)

    # north-south asymmetry
    ax = fig.add_subplot(gs[1, 1])
    ax.axhline(0, color=ps.AXIS, lw=0.8)
    ax.plot(ns["l"], ns["diff"], color=ps.C_CONTEXT, lw=0.8)
    ax.plot(ns["l"], ns["diff_smooth"], color=ps.C_OBSERVED, lw=1.5,
            label="this project (1° smoothing)")
    pub = res["ns_asymmetry"]["published"]
    ax.plot([25, 40, 65], [pub["at_25"], pub["at_40"], pub["at_65"]], "D", color=ps.C_PUBLISHED,
            ms=5, label="McClure-Griffiths & Dickey 2016")
    ax.plot([55, 67], [pub["mean_55_67"]] * 2, color=ps.C_PUBLISHED, lw=2)
    L = ns["l"]
    for pi0 in (7, 14):
        ax.plot(L, 2 * pi0 * np.cos(np.radians(L)), color=ps.MUTED, lw=0.9, ls=(0, (3, 2)))
        ax.text(19, 2 * pi0 * np.cos(np.radians(19)) + 1, f"LSR outflow {pi0} km/s",
                fontsize=7, color=ps.MUTED, ha="left", va="bottom")
    ax.set_xlim(18, 67.5)
    ax.set_ylim(-20, 33)
    ax.set_xlabel("|l| (deg)")
    ax.set_ylabel("north − south |V$_t$| (km/s)")
    ax.set_title("The north-south asymmetry", fontsize=9.5)
    ax.legend(loc="upper right", fontsize=7)
    nflag = sum(len(res["reliability"][q]["flagged_l"]) for q in ("Q1", "Q4"))
    ps.caption(fig, "Top: at the same |sin l| the two sightlines graze the same radius, so a "
               "symmetric galaxy would give identical curves. Gaps: the "
               f"{nflag} of 1442 longitudes where our 12 method variants disagree by more than 3 km/s "
               "(an ambiguous edge). Bottom right: an outward drift of "
               "the local standard of rest would add 2Π₀cos l (dashed); the measured difference "
               "has the wrong shape for that (best fit Π₀ = "
               f"{res['ns_asymmetry']['lsr_outflow_fit_Pi0']:.1f} km/s).")
    return ps.save(fig, FIG / "03_terminal_velocities")


# --------------------------------------------------------------------------
# 4. The rotation curve
# --------------------------------------------------------------------------

def fig_rotation_curve():
    rot = _load("rotation")
    mas = _load("masers")
    eil = reference.eilers19()
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw={"width_ratios": [1.5, 1]})

    for quad, color, name in (("Q1", ps.C_NORTH, "HI tangent points, north (first quadrant)"),
                              ("Q4", ps.C_SOUTH, "HI tangent points, south (fourth quadrant)")):
        good = rot[f"{quad}_good"]
        for a in (ax, ax2):
            order = np.argsort(rot[f"{quad}_R"])
            Rq = rot[f"{quad}_R"][order]
            tq = np.where(good[order], rot[f"{quad}_theta"][order], np.nan)
            a.plot(Rq, tq, color=color, lw=1.0, label=name if a is ax else None)
        c = rot[f"{quad}_bin_centres"]
        ax2.fill_between(c, rot[f"{quad}_bin_sys_lo"], rot[f"{quad}_bin_sys_hi"], color=color,
                         alpha=0.15, lw=0)
    # masers
    g = mas["good"]
    ax.errorbar(mas["R"][g], mas["vrot"][g], xerr=mas["R_err"][g], yerr=mas["vrot_err"][g],
                fmt="o", ms=3, color=ps.C_MASER, ecolor=ps.C_MASER, elinewidth=0.6, alpha=0.8,
                mec=ps.SURFACE, mew=0.4,
                label=f"{g.sum()} masers: parallax + proper motion (Reid et al. 2019 data)")
    # Eilers
    ax.fill_between(eil["R"], eil["vc"] - eil["err_lo"], eil["vc"] + eil["err_hi"],
                    color=ps.C_PUBLISHED, alpha=0.25, lw=0)
    ax.plot(eil["R"], eil["vc"], color=ps.C_PUBLISHED, lw=1.4,
            label="23,000 red giants (Eilers et al. 2019)")
    R = np.linspace(0.5, 26, 300)
    ax.plot(R, galaxy.theta_reid19(R), color=ps.C_MODEL, lw=1.0, ls=(0, (3, 2)),
            label="Reid et al. (2019) fitted rotation curve")
    for a in (ax, ax2):
        a.axvline(cfg.R0_KPC, color=ps.AXIS, lw=0.8)
    ax.text(cfg.R0_KPC + 0.25, 280, "the Sun", fontsize=7.5, color=ps.MUTED, va="top")
    ax.set_xlim(0, 26)
    ax.set_ylim(150, 285)
    ax.set_xlabel("distance from the Galactic centre R (kpc)")
    ax.set_ylabel("circular speed Θ (km/s)")
    ax.set_title("The Milky Way's rotation curve", loc="left")
    ax.legend(loc="lower right", fontsize=7.2)

    # zoom on the HI with constant-uncertainty band and literature fits
    R2 = np.linspace(2.0, 8.2, 100)
    obs = galaxy.Observer.reid19()
    # MGD16 joint fit re-expressed with our constants: their model terminal
    # velocity V_t(l) = 220 (k0 sin l + k1 - sin l), converted as our data are
    s = R2 / obs.r0
    vt_model = cfg.THETA0_IAU_KMS * (0.171 * s + 0.889 - s)
    th_q1 = vt_model + obs.v * s + obs.u * np.sqrt(1 - s ** 2)
    ax2.plot(R2, th_q1, color=ps.C_PUBLISHED, lw=1.4,
             label="McClure-Griffiths & Dickey 2016 joint fit\n(re-expressed with our constants)")
    ax2.plot(R2, galaxy.theta_clemens85(R2 * cfg.R0_IAU_KPC / obs.r0) - cfg.THETA0_IAU_KMS * s
             + obs.v * s, color=ps.YELLOW, lw=1.2,
             label="Clemens 1985 CO polynomial (same conversion)")
    ax2.set_xlim(2.0, 8.2)
    ax2.set_ylim(195, 262)
    ax2.set_xlabel("R (kpc)")
    ax2.set_title("Inner Galaxy, tangent points only", fontsize=9.5)
    ax2.legend(loc="lower right", fontsize=7)
    ps.caption(fig, "HI curves: Θ(R) = V_t + (Sun's orbital motion − the LSR convention) along "
               "the line of sight, at R = R₀|sin l|, with R₀ = 8.15 kpc and the Sun's angular "
               "velocity 30.32 km/s/kpc (Reid et al. 2019). Shaded bands (right): the effect of "
               "those constants' uncertainties. Eilers et al. assumed R₀ = 8.122 kpc.")
    return ps.save(fig, FIG / "04_rotation_curve")


# --------------------------------------------------------------------------
# 5. The Milky Way among other galaxies
# --------------------------------------------------------------------------

def fig_context():
    rot = _load("rotation")
    mas = _load("masers")
    eil = reference.eilers19()
    keep = 0.5 * (eil["err_lo"] + eil["err_hi"]) < 5.0  # drop the sparse outermost bins
    eil = {k: v[keep] for k, v in eil.items()}
    sp = reference.sparc()
    res = _results()["context"]
    fig, (ax, axh, axm) = plt.subplots(1, 3, figsize=(13.5, 4.8),
                                       gridspec_kw={"width_ratios": [1.5, 0.8, 1]})
    for name, g in sp.items():
        if g["quality"] == 3:
            continue
        ax.plot(g["R"], g["Vobs"], color=ps.C_CONTEXT, lw=0.6, zorder=1)
    # the Milky Way
    for quad in ("Q1", "Q4"):
        c, m = rot[f"{quad}_bin_centres"], rot[f"{quad}_bin_median"]
        ax.plot(c, m, color=ps.C_OBSERVED, lw=2.0, zorder=3)
    g = mas["good"]
    edges = np.arange(4, 17, 1.5)
    idx = np.digitize(mas["R"][g], edges)
    mc = [np.median(mas["R"][g][idx == i]) for i in range(1, edges.size) if (idx == i).sum() >= 3]
    mv = [np.median(mas["vrot"][g][idx == i]) for i in range(1, edges.size) if (idx == i).sum() >= 3]
    ax.plot(mc, mv, color=ps.C_OBSERVED, lw=2.0, zorder=3)
    ax.plot(eil["R"], eil["vc"], color=ps.C_OBSERVED, lw=2.0, zorder=3)
    ax.text(17, 245, "Milky Way", fontsize=9, fontweight="bold", color=ps.INK)
    ax.text(38, 60, "175 SPARC galaxies\n(Lelli et al. 2016)", fontsize=8, color=ps.MUTED)
    ax.set_xlim(0, 50)
    ax.set_ylim(0, 330)
    ax.set_xlabel("R (kpc)")
    ax.set_ylabel("rotation speed (km/s)")
    ax.set_title("Our galaxy's rotation among 175 others", loc="left")

    vflat = np.array([g["v_flat"] for g in sp.values()])
    vflat = vflat[np.isfinite(vflat)]
    axh.hist(vflat, bins=np.arange(0, 340, 20), color=ps.C_CONTEXT, edgecolor=ps.SURFACE, lw=1.5)
    axh.axvline(cfg.THETA0_KMS, color=ps.C_OBSERVED, lw=1.6)
    axh.text(cfg.THETA0_KMS + 5, axh.get_ylim()[1] * 0.92,
             f"Milky Way\n{cfg.THETA0_KMS:.0f} km/s:\nfaster than\n"
             f"{100 * (1 - res['sparc_fraction_faster_than_mw']):.0f}%", fontsize=7.5,
             color=ps.INK, va="top")
    axh.yaxis.get_major_locator().set_params(integer=True)
    axh.set_xlabel("flat rotation speed V$_{flat}$ (km/s)")
    axh.set_ylabel("galaxies")
    axh.set_title("How fast is fast?", fontsize=9.5)

    # enclosed mass
    Rm = np.linspace(2.0, 25.0, 200)
    for quad in ("Q1", "Q4"):
        c, m = rot[f"{quad}_bin_centres"], rot[f"{quad}_bin_median"]
        axm.plot(c, galaxy.enclosed_mass(c, m) / 1e10, color=ps.C_OBSERVED, lw=1.6)
    axm.plot(eil["R"], galaxy.enclosed_mass(eil["R"], eil["vc"]) / 1e10, color=ps.C_PUBLISHED,
             lw=1.6, label="from Eilers et al. 2019")
    axm.plot(mc, galaxy.enclosed_mass(np.array(mc), np.array(mv)) / 1e10, "o", color=ps.C_MASER,
             ms=4, label="from masers")
    mr0 = res["enclosed_mass_Msun"]["R0_reid_theta0"] / 1e10
    axm.plot([cfg.R0_KPC], [mr0], marker="*", color=ps.INK, ms=9, ls="")
    axm.text(cfg.R0_KPC + 0.6, mr0 - 1.5, f"inside the Sun's orbit:\n{mr0 / 10:.2f} × 10¹¹ M$_\\odot$",
             fontsize=7.5, va="top")
    axm.set_xlabel("R (kpc)")
    axm.set_ylabel("mass inside R, V²R/G (10¹⁰ M$_\\odot$)")
    axm.set_title("Weighing the Galaxy", fontsize=9.5)
    axm.legend(loc="upper left", fontsize=7)
    axm.set_xlim(0, 26)
    axm.set_ylim(0, 25)
    ps.caption(fig, "Left: Milky Way in black (HI tangent points 2-8 kpc, masers 4-16 kpc, red giants "
               "5-21 kpc, where their errors stay below 5 km/s) over SPARC's rotation curves "
               "(the 12 quality-3 curves omitted). Right: the "
               "mass that a spherical distribution would need to hold each orbit; a flattened "
               "disc needs up to ~30% less.")
    return ps.save(fig, FIG / "05_milky_way_among_galaxies")


# --------------------------------------------------------------------------
# 6-7. The outer disc: face-on map and warp
# --------------------------------------------------------------------------

def _galaxy_frame(ax, r0=cfg.R0_KPC, radii=(10, 15, 20, 25, 30)):
    th = np.linspace(0, 2 * np.pi, 400)
    for r in radii:
        ax.plot(r * np.sin(th), r * np.cos(th), color=ps.AXIS, lw=0.5, zorder=2)
        ax.text(r * np.sin(np.radians(-135)), r * np.cos(np.radians(-135)), f"{r} kpc",
                fontsize=6.5, color=ps.MUTED, ha="center", va="center")
    ax.plot(0, r0, marker="*", color=ps.INK, ms=9, zorder=6)
    ax.text(0.8, r0 + 0.8, "Sun", fontsize=7.5)
    ax.plot(0, 0, marker="+", color=ps.INK, ms=8, zorder=6)
    ax.text(0.8, -1.2, "Galactic\ncentre", fontsize=7, va="top")
    ax.set_aspect("equal")
    ax.set_xlabel("x (kpc)   →  toward l = 90°")
    ax.set_ylabel("y (kpc)")
    ax.grid(False)


def _polar_mesh(vs):
    rc, pc = vs["r_centres"], vs["phi_centres"]
    dr = np.median(np.diff(rc))
    dp = np.median(np.diff(pc))
    re = np.concatenate([rc - dr / 2, [rc[-1] + dr / 2]])
    pe = np.radians(np.concatenate([pc - dp / 2, [pc[-1] + dp / 2]]))
    PP, RR = np.meshgrid(pe, re)
    # LBH06 azimuth: phi = 0 beyond the centre, 90 deg on the l = 90 side
    return RR * np.sin(PP), -RR * np.cos(PP)


def fig_outer_faceon():
    from scipy.ndimage import gaussian_filter

    prod = _load("outer_default")
    circ = _load("outer_circular")
    mas = _load("masers")
    res = _results()
    ell = res["outer_elliptic_fit"]
    obs = galaxy.Observer.reid19()
    e = prod["xy_edges"]
    sig = prod["xy_sigma"]
    x = 0.5 * (e[1:] + e[:-1])
    X, Y = np.meshgrid(x, x)
    R = np.hypot(X, Y)

    fig = plt.figure(figsize=(14.5, 6.3))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.62], height_ratios=[1, 1],
                          wspace=0.25, hspace=0.55)
    ax = fig.add_subplot(gs[:, 0])
    ax2 = fig.add_subplot(gs[:, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[1, 2])

    im = ax.imshow(np.where(sig > 0, sig, np.nan), origin="lower",
                   extent=[e[0], e[-1], e[0], e[-1]], cmap=ps.CMAP_BLUE,
                   norm=LogNorm(vmin=0.03, vmax=15), interpolation="nearest")
    cb = fig.colorbar(im, ax=ax, shrink=0.72, pad=0.09, location="bottom")
    cb.set_label("HI surface density (M$_\\odot$/pc², hydrogen only)")
    _galaxy_frame(ax)
    ax.set_xlabel("")
    ax.set_xlim(-30, 30)
    ax.set_ylim(-30, 30)
    ax.set_title("The outer disc seen from above", loc="left")

    # surface density relative to the fitted exponential decline
    valid = (sig != 0) & (R > obs.r0 + 1) & (R <= 25)
    w = valid.astype(float)
    smooth = gaussian_filter(np.where(valid, sig, 0.0), 1.0) / np.maximum(gaussian_filter(w, 1.0), 1e-9)
    h = res["outer"]["default"]["sigma_scale_length_kpc"]
    prof = res["outer"]["default"]["sigma_R"]
    rr = np.array([float(k) for k in prof])
    ss = np.array([np.nan if v is None else v for v in prof.values()])
    ok = (rr >= 12) & (rr <= 25) & np.isfinite(ss) & (ss > 0)
    k, c = np.polyfit(rr[ok], np.log(ss[ok]), 1)
    rel = np.where(valid, smooth / np.exp(c + k * R), np.nan)
    im2 = ax2.imshow(np.log2(rel), origin="lower", extent=[e[0], e[-1], e[0], e[-1]],
                     cmap=ps.CMAP_DIVERGING, norm=TwoSlopeNorm(0, -1.5, 1.5),
                     interpolation="nearest")
    cb2 = fig.colorbar(im2, ax=ax2, shrink=0.72, pad=0.09, location="bottom")
    cb2.set_ticks([-1, 0, 1], labels=["½×", "typical for its radius", "2×"])
    cb2.set_label(f"surface density relative to an exponential decline (scale length {h:.1f} kpc)")
    _galaxy_frame(ax2, radii=(10, 15, 20, 25))
    ax2.set_xlabel("")
    for arm in galaxy.REID19_ARMS:
        if arm.name not in ("Perseus", "Outer"):
            continue
        beta = np.linspace(*arm.beta_draw, 2000)
        xa, ya = arm.xy(beta)
        ind = (beta >= arm.beta_data[0]) & (beta <= arm.beta_data[1])
        outside = np.hypot(xa, ya) > obs.r0
        ax2.plot(np.where(outside, xa, np.nan), np.where(outside, ya, np.nan), color=ps.INK,
                 lw=0.9, ls=(0, (3, 2)))
        ax2.plot(np.where(ind & outside, xa, np.nan), np.where(ind & outside, ya, np.nan),
                 color=ps.INK, lw=1.5)
    box = dict(boxstyle="round,pad=0.15", fc=ps.SURFACE, ec="none", alpha=0.8)
    ax2.text(-12.0, 11.3, "Perseus", fontsize=8, fontweight="bold", ha="center", bbox=box)
    ax2.text(-14.0, 17.2, "Outer", fontsize=8, fontweight="bold", ha="center", bbox=box)
    sel = np.isin(mas["arm"], ["Per", "Out"]) & (mas["e_plx"] / mas["plx"] < 0.2)
    xm, ym = galaxy.xy_from_ld(mas["l"][sel], 1 / mas["plx"][sel], obs.r0)
    ax2.scatter(xm, ym, s=12, color=ps.C_MASER, edgecolor=ps.SURFACE, lw=0.4, zorder=7,
                label="Perseus and Outer arm masers (parallax distances)")
    ax2.plot([], [], color=ps.INK, lw=1.5, label="arms fitted to them (Reid et al. 2019)")
    ax2.legend(loc="lower left", fontsize=7)
    ax2.set_xlim(-26, 26)
    ax2.set_ylim(-26, 26)
    ax2.set_title("Where the gas is denser or thinner than usual", fontsize=9.5)

    # the anticentre matching that fixes the elliptical-orbit correction
    rb = np.arange(10.0, 30.01, 1.0)
    rc = 0.5 * (rb[1:] + rb[:-1])
    for a, p, title in ((ax3, circ, "Either side of the anticentre:\ncircular orbits"),
                        (ax4, prod, f"with radial motions fitted\n(α = {ell['ours']['alpha']:.1f}, "
                                    f"β = {ell['ours']['beta']:.2f} km/s)")):
        s165, s195 = outer.anticentre_profiles(p, obs.r0, rb)
        a.semilogy(rc, s165, color=ps.C_NORTH, lw=1.5, label="l = 155-165°")
        a.semilogy(rc, s195, color=ps.C_SOUTH, lw=1.5, label="l = 195-205°")
        a.set_xlim(10, 30)
        a.set_ylim(0.01, 20)
        a.set_ylabel("Σ (M$_\\odot$/pc²)")
        a.set_title(title, fontsize=8.5)
    ax3.legend(loc="upper right", fontsize=7)
    ax4.set_xlabel("R (kpc)")
    ps.caption(fig, "Each HI4PI brightness measurement placed where Reid et al.'s rotation curve puts "
               "gas of that velocity, outside the solar circle only; white wedges: within 15° of the "
               "Galactic centre and anticentre, where velocities are too small to give distances. "
               "Radial streaks are velocity noise stretched along sightlines. Circular orbits leave a "
               "jump in density across l = 180° (top right); a small radial motion v_R = v_π(R) cos φ, "
               "fitted as Levine, Blitz & Heiles (2006) did, removes it and is used throughout.")
    return ps.save(fig, FIG / "06_outer_disc")


def fig_warp():
    vs = _load("outer_default_vertical")
    res = _results()["outer"]
    fig = plt.figure(figsize=(13.5, 6.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1], hspace=0.45, wspace=0.2)
    ax = fig.add_subplot(gs[:, 0])
    X, Y = _polar_mesh(vs)
    zb = np.where((vs["sigma"] > 0.02) & (vs["r_centres"][:, None] <= 30), vs["zbar"], np.nan)
    m = ax.pcolormesh(X, Y, zb, cmap=ps.CMAP_DIVERGING, norm=TwoSlopeNorm(0, -3, 3),
                      shading="flat", rasterized=True)
    cb = fig.colorbar(m, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label("mean height of the hydrogen above the plane (kpc)")
    _galaxy_frame(ax)
    ax.set_xlim(-30, 30)
    ax.set_ylim(-30, 30)
    ax.set_title("The warp: the gas disc bends up on one side, down on the other", loc="left")

    # warp amplitudes vs R, like for like with LBH06
    ax2 = fig.add_subplot(gs[0, 1])
    Rl = np.linspace(8, 30, 200)
    lbh = outer.lbh06_warp_amplitudes(Rl)
    names = {"W0": "W0: whole disc lifted", "W1": "W1: tilt, up one side, down the other",
             "W2": "W2: both sides up (the asymmetry)"}
    for key, color in (("W0", ps.MUTED), ("W1", ps.C_NORTH), ("W2", ps.C_SOUTH)):
        rr = [r["R"] for r in res["lbh06"]["rings"]]
        ww = [r[key] for r in res["lbh06"]["rings"]]
        ax2.plot(rr, np.array(ww) / 1e3, "o", color=color, ms=3.5)
        ax2.plot(Rl, lbh[key] / 1e3, color=color, lw=1.3)
        ax2.text(Rl[-1] + 0.4, lbh[key][-1] / 1e3, key, color=color, fontsize=8, va="center",
                 fontweight="bold")
    ax2.set_xlim(10, 32)
    ax2.set_xlabel("R (kpc), on LBH06's scale (R$_0$ = 8.5 kpc, flat 220 km/s)")
    ax2.set_ylabel("warp amplitude (kpc)")
    ax2.set_title("Warp modes: this project (dots) vs Levine, Blitz & Heiles 2006 (lines)",
                  fontsize=9.5)

    # flaring: dispersion of the gas about its mean height
    ax3 = fig.add_subplot(gs[1, 1])
    rr = vs["r_centres"]
    disp = np.nanmedian(vs["disp"], axis=1)
    fwhm = 2 * np.sqrt(2 * np.log(2)) * disp
    ax3.plot(rr, fwhm, color=ps.INK, lw=1.5)
    ax3.set_xlim(9, 32)
    ax3.set_ylim(0, None)
    ax3.set_xlabel("R (kpc)")
    ax3.set_ylabel("disc thickness, FWHM (kpc)")
    ax3.set_title("The disc flares: it thickens outward", fontsize=9.5)
    ps.caption(fig, "Left: mean height of the hydrogen after dropping gas more than two dispersions "
               "off the disc (LBH06's filter), cells with Σ > 0.02 M☉/pc². Top right: the warp fitted "
               "ring by ring as W0 + W1 sin(φ−φ1) + W2 sin(2φ−φ2), rerun with LBH06's own "
               "Galactic constants and elliptical-orbit correction so the radii are on their scale. "
               "Bottom right: median over azimuth.")
    return ps.save(fig, FIG / "07_warp")


def fig_comparison():
    """Our numbers against the published ones, in two panels: velocities
    (km/s differences) and structure (% differences)."""
    r = _results()
    t = r["terminal"]
    mz = r["masers"]
    el = r["outer_elliptic_fit"]
    rings = {ring["R"]: ring for ring in r["outer"]["lbh06"]["rings"]}
    lbh = outer.lbh06_warp_amplitudes(np.array([20.5]))
    mgd = r["rotation"]["mgd_style_fit"]
    th5 = lambda k0, k1: cfg.THETA0_IAU_KMS * (k0 * 5.0 / cfg.R0_IAU_KPC + k1)
    vel = [  # (label, ours - published, published uncertainty or None)
        ("terminal velocities, first quadrant, median\nof 748 longitudes (VLA, MGD 2016)",
         t["q1_vs_mgd16"]["median"], None),
        ("terminal velocities, fourth quadrant,\n10 printed values (MGD 2007)",
         t["q4_vs_mgd07_table"]["mean"], None),
        ("linear fit v_circ, first quadrant", t["linear_fit"]["Q1"]["v_circ"] - 172.5, None),
        ("linear fit v_circ, fourth quadrant", t["linear_fit"]["Q4"]["v_circ"] - 177.7, None),
        ("linear fit v', first quadrant", t["linear_fit"]["Q1"]["v_prime"] - 20.7, None),
        ("linear fit v', fourth quadrant", t["linear_fit"]["Q4"]["v_prime"] - 18.4, None),
        ("north − south, mean over |l| = 55-67°", t["ns_asymmetry"]["mean_55_67"] + 5.2, None),
        ("joint rotation fit at R = 5 kpc (IAU constants)",
         th5(mgd["k0"], mgd["k1"]) - th5(0.171, 0.889), None),
        ("Oort A·R₀, fourth quadrant, |sin l| 0.35-0.95", t["AR0"]["Q4_0.35-0.95"]["AR0"] - 89, 4),
        ("Oort A·R₀, fourth quadrant, |sin l| 0.80-0.95", t["AR0"]["Q4_0.80-0.95"]["AR0"] - 110, 2),
        ("maser lag behind rotation (Reid et al. 2019)",
         -mz["vrot_minus_reid_urc"]["median"] - 4.3, 5.6),
        ("maser inward drift (Reid et al. 2019)", -mz["vrad_median_outward"] - 6.0, 1.4),
        ("elliptical-orbit α, LBH06 constants", el["lbh06_constants"]["alpha"] - 8.67, None),
        ("elliptical-orbit β, LBH06 constants", el["lbh06_constants"]["beta"] + 1.08, None),
    ]
    struct = [
        ("warp W0 at R = 20.5 kpc (LBH06)", 100 * (rings[20.5]["W0"] / lbh["W0"][0] - 1)),
        ("warp W1 at R = 20.5 kpc (LBH06)", 100 * (rings[20.5]["W1"] / lbh["W1"][0] - 1)),
        ("warp W2 at R = 20.5 kpc (LBH06)", 100 * (rings[20.5]["W2"] / lbh["W2"][0] - 1)),
        ("HI scale length, same constants, co-rotating\ngas (Kalberla & Dedes 2008: 4.0 kpc)",
         100 * (r["outer"]["lbh06"]["sigma_scale_length_kpc"] / 4.0 - 1)),
    ]
    fig = plt.figure(figsize=(12.5, 8.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[len(vel), len(struct) + 0.6], hspace=0.28)
    ax = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    for a, rows, unit, lim in ((ax, vel, "km/s", 10), (ax2, [(l, d, None) for l, d in struct], "%", 30)):
        y = np.arange(len(rows))[::-1]
        for yi, (label, d, err) in zip(y, rows):
            if err is not None:
                a.barh(yi, 2 * err, left=-err, height=0.8, color=ps.C_PUBLISHED, alpha=0.15, lw=0)
            a.barh(yi, d, height=0.45, color=ps.ORANGE)
            a.text(d + (0.25 if unit == "km/s" else 0.8) * np.sign(d if d else 1), yi,
                   f"{d:+.2f}" if unit == "km/s" else f"{d:+.0f}%", va="center",
                   ha="left" if d >= 0 else "right", fontsize=7.5, color=ps.INK_SECONDARY)
        a.axvline(0, color=ps.INK_SECONDARY, lw=1)
        a.set_yticks(y, [row[0] for row in rows], fontsize=8)
        a.set_xlim(-lim, lim)
        a.grid(axis="y", visible=False)
        a.tick_params(axis="y", length=0)
    ax.set_xlabel("this project − published (km/s)")
    ax2.set_xlabel("this project relative to published (%)")
    ax.set_title("Our numbers against the literature", loc="left", pad=22)
    ax.text(0, 1.01, "violet band: the published uncertainty, where one is quoted",
            transform=ax.transAxes, fontsize=8, color=ps.INK_SECONDARY, va="bottom")
    ps.caption(fig, "Published: McClure-Griffiths & Dickey 2007 (ApJ 671, 427) and 2016 (ApJ 831, "
               "124); Reid et al. 2019 (ApJ 885, 131); Levine, Blitz & Heiles 2006 (ApJ 643, 881); "
               "Kalberla & Dedes 2008 (A&A 487, 951). The one clear miss, Oort A·R₀ fitted over "
               "|sin l| = 0.80-0.95, swings from 48 to 133 km/s in our data when that range shifts "
               "by 0.05: it measures local streaming more than rotation, as MGD 2007 warned.")
    return ps.save(fig, FIG / "00_comparison")


def make_all() -> list:
    written = []
    for f in (fig_comparison, fig_lv, fig_terminal_method, fig_terminal_curve, fig_rotation_curve,
              fig_context, fig_outer_faceon, fig_warp):
        written += f()
    return written


if __name__ == "__main__":
    for p in make_all():
        print(p)
