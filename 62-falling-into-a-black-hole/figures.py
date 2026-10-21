"""Every figure of the project (PNG + SVG in figures/).

The calculations live in infall.py, raymap.py, render.py and kerr.py; this
module only draws them. Views of the sky are rendered at modest resolution
here (the movie has the full-resolution ones).
"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Polygon

import config as C
import infall
import kerr
import plotstyle as P
import raymap
import render as RD

P.use_style()
OUT = C.FIGURE_DIR


def _fmt_pow(x: float, digits: int = 1) -> str:
    """1.3x10^-4 with Unicode superscripts."""
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10**e
    sup = str(e).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{m:.{digits}f}×10{sup}" if digits else f"10{sup}"


def _fmt_mass(m: float) -> str:
    if m >= 1e9:
        return f"{m / 1e9:.1f} billion M⊙"
    if m >= 1e6:
        return f"{m / 1e6:.2f} million M⊙"
    return f"{m:.3g} M⊙"


def _fmt_time(s: float) -> str:
    if s >= 3600:
        return f"{s / 3600:.1f} h"
    if s >= 60:
        return f"{s / 60:.1f} min"
    if s >= 1:
        return f"{s:.1f} s"
    if s >= 0.1:
        return f"{s:.2f} s"
    if s >= 0.01:
        return f"{s * 1e3:.0f} ms"
    if s >= 1e-3:
        return f"{s * 1e3:.2f} ms"
    return f"{s * 1e6:.0f} µs"


# ---------------------------------------------------------------------------
# 1. Two clocks
# ---------------------------------------------------------------------------

def clocks(bhs: dict):
    tu = bhs["sgra"]["gm_c3_s"]
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), gridspec_kw={"width_ratios": [1.15, 1, 1.05]})

    # (a) radius against the faller's clock and against a far-away clock, Sgr A*
    ax = axs[0]
    r0 = 20.0                                          # 10 r_s
    r = np.concatenate([np.geomspace(r0, 2.0 + 1e-9, 3000), np.linspace(2.0, 0.0, 400)[1:]])
    tau = (infall.rain_tau_left(r0) - infall.rain_tau_left(r)) * tu
    out = r > 2.0
    t = (infall.rain_schwarzschild_t(r[out]) - infall.rain_schwarzschild_t(r0)) * tu
    ax.plot(tau / 60, r / 2, color=P.BLUE, label="faller's own clock")
    ax.plot(t / 60, r[out] / 2, color=P.ORANGE, label="far-away clock (Schwarzschild t)")
    # dropped from rest near the horizon: released at 1.5 r_s (the photon sphere)
    rd0 = 3.0
    eta = np.linspace(0, math.pi, 800)
    rd, taud = infall.drop_r_tau(eta, rd0)
    # released from a ship hovering at 1.5 r_s just as the first faller passes: the two
    # clocks start from that event's readings (own: the first faller's watch; far: its t)
    tau_rel = float(infall.rain_tau_left(r0) - infall.rain_tau_left(rd0)) * tu
    t_rel = float(infall.rain_schwarzschild_t(rd0) - infall.rain_schwarzschild_t(r0)) * tu
    ax.plot((tau_rel + taud * tu) / 60, rd / 2, color=P.BLUE, ls="--", lw=1.2,
            label="dropped from rest at 1.5 rₛ: own clock")
    od = rd > 2.0
    ax.plot((t_rel + infall.drop_schwarzschild_t(eta[od], rd0) * tu) / 60, rd[od] / 2, color=P.ORANGE, ls="--", lw=1.2,
            label="dropped from rest at 1.5 rₛ: far-away clock")
    ax.axhline(1.0, color=P.INK, ls="--", lw=0.8)
    ax.text(0.2, 1.05, "horizon", fontsize=8, color=P.INK)
    ax.set_xlim(0, 22)
    ax.set_ylim(0, 10.2)
    ax.set_xlabel("time since passing 10 rₛ (minutes, Sgr A*)")
    ax.set_ylabel("distance from the centre (rₛ)")
    ax.set_title("a. Two clocks watch the same fall")
    ax.legend(loc="upper right", fontsize=7.5)
    P.annotate(ax, f"own clock: horizon after {infall.rain_tau_left(r0) * tu / 60 - infall.rain_tau_left(2.0) * tu / 60:.1f} min,\n"
               f"centre {infall.rain_tau_left(2.0) * tu:.0f} s later", xy=(tau[np.argmin(np.abs(r - 2))] / 60, 1.0),
               xytext=(1.5, 4.0))
    P.annotate(ax, "far-away clock: never\nsees the crossing", xy=(19.5, 1.07), xytext=(14.5, 3.4))

    # (b) time inside against release height
    ax = axs[1]
    r0s = np.geomspace(2.0 + 1e-7, 2e5, 400)
    ins = np.array([infall.drop_tau_inside(x) for x in r0s])
    ax.semilogx(r0s / 2, ins, color=P.BLUE)
    ax.axhline(math.pi, color=P.INK, ls=":", lw=0.8)
    ax.axhline(4 / 3, color=P.INK, ls=":", lw=0.8)
    ax.text(1.3e4, math.pi + 0.05, "π GM/c³: the most anyone gets", fontsize=8, ha="right")
    ax.text(1.3e4, 4 / 3 + 0.05, "4/3 GM/c³: falling from far away", fontsize=8, ha="right")
    ax.set_xlabel("released from rest at r₀ (rₛ)")
    ax.set_ylabel("time from horizon to singularity (GM/c³)")
    ax.set_ylim(0, 3.6)
    sec = ax.secondary_yaxis("right", functions=(lambda y: y * tu, lambda s: s / tu))
    sec.set_ylabel("for Sgr A* (seconds)")
    ax.set_title("b. Inside, the singularity is a moment, not a place")

    # (c) the real black holes
    ax = axs[2]
    keys = list(C.BH_ORDER)
    y = np.arange(len(keys))[::-1]
    for yi, k in zip(y, keys):
        b = bhs[k]
        st = P.BH_STYLE[k]
        ax.barh(yi + 0.18, b["horizon_to_singularity_max_s"], height=0.34, color=st["color"], alpha=0.35)
        ax.barh(yi - 0.18, b["horizon_to_singularity_rain_s"], height=0.34, color=st["color"])
        ax.text(b["horizon_to_singularity_max_s"] * 1.4, yi + 0.18, _fmt_time(b["horizon_to_singularity_max_s"]),
                va="center", fontsize=7.5, color=P.INK_SECONDARY)
        ax.text(b["horizon_to_singularity_rain_s"] * 1.4, yi - 0.18, _fmt_time(b["horizon_to_singularity_rain_s"]),
                va="center", fontsize=7.5, color=P.INK)
    ax.set_xscale("log")
    ax.set_yticks(y, [f"{bhs[k]['name']}\n{_fmt_mass(bhs[k]['mass_msun'])}" for k in keys], fontsize=8)
    ax.set_xlim(1e-5, 3e6)
    ax.set_xlabel("time from horizon to singularity (s)")
    ax.set_title("c. The last moments, by your own clock")
    ax.text(0.98, 0.03, "solid: falling from far away\npale: longest possible (π GM/c³)", transform=ax.transAxes,
            ha="right", fontsize=7.5, color=P.INK_SECONDARY)
    ax.grid(axis="y", visible=False)
    P.caption(fig, "Schwarzschild (non-spinning) black holes, radial fall. Masses: Sgr A* from S2's orbit (H1), M87* from its "
              "EHT ring (I5), Cyg X-1 from its orbit (I4), GW250114's remnant (LVK 2025), Gaia BH1 (El-Badry et al. 2023).")
    fig.tight_layout()
    return P.save(fig, OUT / "01_clocks")


# ---------------------------------------------------------------------------
# 2. Tides
# ---------------------------------------------------------------------------

def tides(bhs: dict, uni: dict):
    fig, axs = plt.subplots(1, 2, figsize=(13.5, 5.0))
    L = C.BODY_LENGTH_M
    ax = axs[0]
    offsets = {"m87": 3.0, "sgra": 3.0, "gw250114": 0.06, "cygx1": 1.0, "gaiabh1": 15.0}
    for k in C.BH_ORDER:
        b = bhs[k]
        st = P.BH_STYLE[k]
        rs_km = b["rs_km"]
        r_km = np.geomspace(rs_km * 1e-3, rs_km * 1e5, 400)
        a = infall.tidal_stretch_si(r_km * 1e3, b["mass_msun"], L, C.GM_SUN) / C.G0
        ax.loglog(r_km, a, color=st["color"], ls=st["ls"], lw=1.5)
        ax.plot(rs_km, b["stretch_at_horizon_g"], "o", color=st["color"], ms=6, mec=P.SURFACE)
        ax.text(rs_km * 2.2, b["stretch_at_horizon_g"] * offsets[k],
                f"{st['label']}: {_fmt_pow(b['stretch_at_horizon_g'])} g at the horizon", color=st["color"], fontsize=7.5,
                va="center")
    ax.axhline(C.SURVIVAL_G, color=P.INK, lw=0.9, ls="--")
    ax.fill_between([1e-3, 1e18], C.SURVIVAL_G, 1e14, color=P.RED, alpha=0.06, lw=0)
    ax.text(2e9, C.SURVIVAL_G * 3, "10 g head to foot: spaghettification", fontsize=8)
    ax.set_xlim(1e0, 1e16)
    ax.set_ylim(1e-14, 1e12)
    ax.set_xlabel("distance from the centre (km)")
    ax.set_ylabel("stretch across a 2 m body (g)")
    ax.set_title("a. Tides: dots mark each horizon")
    P.annotate(ax, f"Cyg X-1 tears a body apart at {bhs['cygx1']['r_10g_km']:,.0f} km,\n"
               f"{bhs['cygx1']['r_10g_rs']:.0f} times its horizon radius",
               xy=(bhs["cygx1"]["r_10g_km"], C.SURVIVAL_G), xytext=(3e9, 1e6))

    # (b) against time left: one universal curve, Delta a = (4/9) L / tau^2
    ax = axs[1]
    tau = np.geomspace(1e-5, 1e7, 400)
    ax.loglog(tau, 4.0 / 9.0 * L / tau**2 / C.G0, color=P.INK, lw=1.4)
    for k in C.BH_ORDER:
        b = bhs[k]
        st = P.BH_STYLE[k]
        th = b["horizon_to_singularity_rain_s"]
        ax.plot(th, 4.0 / 9.0 * L / th**2 / C.G0, "o", color=st["color"], ms=7, mec=P.SURFACE)
        ax.text(th * 1.6, 4.0 / 9.0 * L / th**2 / C.G0 * 1.6, st["label"] + " horizon", color=st["color"], fontsize=8)
    t10 = bhs["sgra"]["time_from_10g_to_singularity_s"]
    ax.axhline(C.SURVIVAL_G, color=P.INK, lw=0.9, ls="--")
    ax.axvline(t10, color=P.RED, lw=0.9)
    ax.text(t10 * 1.3, 1e-9, f"{t10:.3f} s before the end,\nfor every black hole", color=P.RED, fontsize=8)
    ax.set_xlabel("your own time left before the singularity (s)")
    ax.set_ylabel("stretch across a 2 m body (g)")
    ax.set_title("b. Against your own clock: one curve for all")
    ax.set_xlim(1e-5, 1e7)
    ax.set_ylim(1e-14, 1e12)
    P.caption(fig, "Falling from rest far away, the stretch is 2GML/r³ (exact in GR for radial falls; see symbolic.py) and "
              "the time left is (2/3)√(r³/2GM), so stretch = (4/9) L / τ² whatever the mass. Mass only "
              f"decides where the horizon falls on that curve: above ~{round(uni['critical_mass_10g_at_horizon_msun'], -2):,.0f} M⊙ "
              "you cross it intact.")
    fig.tight_layout()
    return P.save(fig, OUT / "02_tides")


# ---------------------------------------------------------------------------
# 3. What a far-away friend sees
# ---------------------------------------------------------------------------

def friend(bhs: dict):
    tu = bhs["sgra"]["gm_c3_s"]
    r_f = 2000.0
    tau, t, z1, r = infall.rain_received(r_f, 20.0, n=20000)
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.6))
    ax = axs[0]
    ax.semilogy(t * tu / 60, z1, color=P.ORANGE, label="1 + z of the faller's light")
    ax.semilogy(t * tu / 60, z1**-4 * z1[0] ** 4, color=P.BLUE, ls="--",
                label="its brightness (relative): the (1+z)⁻⁴ factor")
    efold = infall.late_time_efold(t, z1)
    sel = z1 > 50
    ax.semilogy(t[sel] * tu / 60, np.exp(t[sel] / efold) * z1[sel][0] / np.exp(t[sel][0] / efold), color=P.INK, lw=0.8, ls=":")
    ax.set_ylim(1e-8, 1e6)
    ax.set_xlim(0, 45)
    ax.set_xlabel("friend's clock (minutes since the faller passed 10 rₛ, Sgr A*)")
    ax.set_ylabel("redshift factor 1 + z,  relative brightness")
    ax.set_title("a. The faller reddens and fades")
    ax.legend(loc="lower left", fontsize=7.5)
    ax.text(0.03, 0.62, f"late on, the redshift grows ×e every\n{efold:.3f} GM/c³ = {efold * tu:.0f} s (dotted), and the\n"
            f"brightness falls ×e every GM/c³ = {tu:.0f} s", transform=ax.transAxes, fontsize=8)
    ax = axs[1]
    ax.plot(t * tu / 60, tau * tu / 60, color=P.BLUE)
    cross = (infall.rain_tau_left(20.0) - infall.rain_tau_left(2.0)) * tu / 60
    ax.axhline(cross, color=P.INK, ls="--", lw=0.8)
    ax.text(0.5, cross + 0.3, f"the faller's watch freezes at {cross:.2f} min: the horizon crossing", fontsize=8)
    ax.set_xlabel("friend's clock (minutes)")
    ax.set_ylabel("faller's watch as the friend sees it (minutes)")
    ax.set_title("b. The faller's watch, seen through a telescope")
    ax.set_ylim(0, cross * 1.25)
    ax.set_xlim(0, 45)
    P.caption(fig, f"A friend hovering at {r_f / 2:,.0f} rₛ directly above receives light sent straight up. In the far-away "
              "friend's time the faller never crosses; in their own time they cross after 13 minutes and are gone 26 s later. "
              "The brightness curve keeps only the (1+z)⁻⁴ factor (photon energy, arrival rate and beaming), not the "
              "image's apparent size, which stops changing near the horizon. Light can also circle the photon sphere "
              "before escaping, so faint delayed echoes outlast this direct image; for a "
              "collapsing star they dominate the last light, which fades as exp(−t/3√3 GM/c³) (Ames & Thorne 1968).")
    fig.tight_layout()
    return P.save(fig, OUT / "03_friend_view")


# ---------------------------------------------------------------------------
# 4. Penrose diagram
# ---------------------------------------------------------------------------

PENROSE_P = 0.5       # compactification exponent (see infall.penrose)


def penrose(bhs: dict, p: float = PENROSE_P):
    tu = bhs["sgra"]["gm_c3_s"]
    fig, ax = plt.subplots(figsize=(8.2, 6.8))
    h = math.pi / 2
    ext = Polygon([(0, 0), (h, -h), (2 * h, 0), (h, h)], closed=True, fc=P.WASH, ec="none")
    bh = Polygon([(0, 0), (h, h), (-h, h)], closed=True, fc="#dfe7f4", ec="none")
    ax.add_patch(ext)
    ax.add_patch(bh)
    ax.plot([0, h], [0, h], color=P.INK, lw=1.4)
    xs = np.linspace(-h, h, 400)
    ax.plot(xs, h + 0.025 * np.sin(xs * 40), color=P.INK, lw=1.2)
    ax.plot([h, 2 * h, h], [h, 0, -h], color=P.INK_SECONDARY, lw=0.8)
    ax.plot([0, h], [0, -h], color=P.INK_SECONDARY, lw=0.8, ls=":")
    ax.text(-0.2, h + 0.1, "singularity, r = 0: a moment in the future of everything inside", ha="center", fontsize=8.5)
    ax.text(h * 0.22 - 0.06, h * 0.22 + 0.06, "event horizon", rotation=45, fontsize=8.5, ha="center", va="bottom")
    ax.text(1.5 * h + 0.12, 0.5 * h + 0.06, "future null infinity", rotation=-45, fontsize=8, ha="center",
            color=P.INK_SECONDARY)
    ax.text(2 * h + 0.05, 0, "i⁰", fontsize=10, va="center")
    ax.text(h, h + 0.07, "i⁺", fontsize=10, ha="center")
    ax.text(h, -h - 0.16, "i⁻", fontsize=10, ha="center")
    ax.text(1.42 * h, -0.28, "outside\n(region I)", fontsize=9, ha="center", color=P.INK_SECONDARY)
    ax.text(-0.42 * h, 0.66 * h, "inside\n(region II)", fontsize=9, ha="center", color=P.INK_SECONDARY)

    r = np.concatenate([np.geomspace(5e4, 2.0001, 4000), np.linspace(2.0, 1e-5, 600)])
    T = -infall.rain_tau_left(r)
    v = infall.rain_ef_v(r, T)
    U, V = infall.kruskal_from_v_r(v, r)
    x, y = infall.penrose(U, V, p)
    ax.plot(x, y, color=P.BLUE, lw=2.0, label="faller (from rest far away)", zorder=4)
    friend_r = 10.0
    rs_f = float(infall.tortoise(friend_r))
    t_f = np.linspace(-3000, 3000, 6000)
    Uf, Vf = infall.static_world_line(friend_r, t_f)
    xf, yf = infall.penrose(Uf, Vf, p)
    ax.plot(xf, yf, color=P.ORANGE, lw=1.8, label=f"friend hovering at {friend_r / 2:.0f} rₛ", zorder=4)
    tau_left = infall.rain_tau_left(r) * tu
    i_h = int(np.argmin(np.abs(r - 2.0)))
    ax.plot(x[i_h], y[i_h], "o", ms=9, mfc="none", mec=P.INK, mew=1.2, zorder=6)
    P.annotate(ax, f"crossing the horizon,\n{tau_left[i_h]:.0f} s before the end", xy=(x[i_h], y[i_h]),
               xytext=(x[i_h] + 0.75, y[i_h] - 0.2))
    for m in (600, 300, 120, 60, 40, 13.0):
        i = int(np.argmin(np.abs(tau_left - m)))
        Ui, Vi = U[i], V[i]
        ax.plot(x[i], y[i], "o", color=P.BLUE, ms=4, zorder=5)
        lab = f"{m / 60:.0f} min" if m >= 60 else f"{m:.0f} s"
        ax.text(x[i] - 0.05, y[i], lab, fontsize=7, ha="right", va="center", color=P.BLUE)
        if r[i] >= friend_r:
            continue            # still above the friend: light sent outward never reaches them
        if Ui < 0:
            t_hit = rs_f - 4.0 * math.log(-Ui)
            Vh = math.exp((t_hit + rs_f) / 4.0)
            Vs = np.geomspace(Vi, Vh, 80)
            ls = "-"
        else:
            Vs = np.linspace(Vi, 1.0 / Ui, 80)
            ls = "--"
        px, py = infall.penrose(np.full_like(Vs, Ui), Vs, p)
        ax.plot(px, py, color=P.MUTED, lw=0.7, ls=ls, zorder=3)
    ax.set_xlim(-h - 0.1, 2 * h + 0.25)
    ax.set_ylim(-h - 0.25, h + 0.3)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.legend(loc="lower left", fontsize=8)
    ax.set_title("Penrose diagram: falling into Sgr A*", loc="left")
    P.caption(fig, "Light travels at 45°. Dots: the faller's own time left before the singularity; grey lines: light "
              "they send outward once below the friend. Signals sent before the horizon reach the friend later and later; "
              "those sent after it (dashed) end on the singularity. Outside (I) and black-hole (II) regions of Schwarzschild spacetime, "
              f"compactified as arctan(|U|^{p:g}), arctan(|V|^{p:g}) of the Kruskal-Szekeres coordinates.")
    fig.tight_layout()
    return P.save(fig, OUT / "04_penrose")


# ---------------------------------------------------------------------------
# 5. The river model
# ---------------------------------------------------------------------------

def river():
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1, 1.15]})
    ax = axs[0]
    r = np.linspace(0.08, 6.0, 600)            # in r_s
    b = 1.0 / np.sqrt(r)                       # sqrt(2M/r) with r_s = 2M
    ax.plot(r, -b, color=P.BLUE, lw=2, label="the river (space), and a raindrop in it")
    ax.plot(r, 1 - b, color=P.ORANGE, lw=1.6, label="light aimed outward")
    ax.plot(r, -1 - b, color=P.ORANGE, lw=1.6, ls="--", label="light aimed inward")
    ax.axhline(0, color=P.AXIS, lw=0.8)
    ax.axvline(1.0, color=P.INK, ls="--", lw=0.8)
    ax.text(1.05, 0.55, "horizon: the river\nflows at c", fontsize=8)
    ax.set_ylim(-4.5, 1.0)
    ax.set_xlabel("distance from the centre (rₛ)")
    ax.set_ylabel("speed dr/dT (units of c; negative = inward)")
    ax.set_title("a. Speeds in the river (Gullstrand-Painlevé time)")
    ax.legend(loc="lower right", fontsize=7.5)

    ax = axs[1]
    ax.axvspan(0, 1.0, color="#101010")
    ax.axvline(1.0, color=P.INK, lw=1.0, ls="--")
    dt = 0.32
    rows = (2.8, 1.8, 1.3, 1.0, 0.7)
    for j, r0 in enumerate(rows):
        yc = 2.0 - j * 1.0
        beta = 1.0 / math.sqrt(r0)
        drift = beta * dt
        inside = r0 < 1.0
        ax.annotate("", xy=(r0 - drift, yc), xytext=(r0, yc),
                    arrowprops={"arrowstyle": "->", "color": P.BLUE_300 if not inside else P.BLUE_150, "lw": 1.4})
        ax.add_patch(Circle((r0 - drift, yc), dt, fc="none", ec=P.ORANGE, lw=1.4))
        ax.plot([r0], [yc], "o", color=P.ORANGE if inside else P.INK, ms=4, mec="white" if inside else P.INK)
        out_edge = r0 - drift + dt
        ax.plot([out_edge], [yc], ">", color=P.ORANGE, ms=6)
        msg = ("outward light gains ground" if out_edge > r0 + 1e-9 else
               "outward light stands still" if abs(out_edge - r0) < 1e-9 else "outward light is swept in")
        ax.text(3.35, yc, f"{r0:g} rₛ: {msg}", fontsize=8, va="center",
                color=P.INK if not inside else P.INK_SECONDARY)
    ax.text(0.5, 2.6, "inside", color="white", fontsize=9, ha="center")
    ax.set_xlim(0, 5.6)
    ax.set_ylim(-2.6, 2.9)
    ax.set_aspect("equal")
    ax.set_xticks([0, 1, 2, 3])
    ax.set_yticks([])
    ax.set_xlabel("distance from the centre (rₛ)")
    ax.grid(False)
    ax.set_title("b. A flash of light, a moment later, carried by the river", loc="left")
    P.caption(fig, "Hamilton & Lisle 2008: in Gullstrand-Painlevé coordinates space is flat (verified in symbolic.py) and "
              "flows inward at the Newtonian escape speed √(2GM/r). Light moves at c relative to the flowing space "
              "(circles, radius cΔt, centres carried inward by the flow, blue arrows): at the horizon the flow reaches "
              "c and outward-aimed light (triangles) stands still; inside, it too is swept toward the centre.")
    fig.tight_layout()
    return P.save(fig, OUT / "05_river")


# ---------------------------------------------------------------------------
# 6. How big the dark region looks, and the colour shifts
# ---------------------------------------------------------------------------

def shadow_size():
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
    ax = axs[0]
    r = np.concatenate([np.geomspace(2000.0, 2.0005, 300), np.geomspace(1.9995, 0.01, 200)])
    fall = np.degrees([raymap.shadow_edge(x) for x in r])
    out = r > 2.0
    hov = []
    for x in r[out]:
        s = math.sqrt(27.0) * math.sqrt(1 - 2 / x) / x
        hov.append(math.degrees(math.asin(min(s, 1.0)) if x >= 3 else math.pi - math.asin(s)))
    ax.semilogx(r / 2, fall, color=P.FALLING, lw=2, label="falling from rest far away")
    ax.semilogx(r[out] / 2, hov, color=P.HOVERING, lw=2, ls="--", label="hovering on rockets")
    ax.axvline(1.0, color=P.INK, ls="--", lw=0.8)
    ax.axvline(1.5, color=P.INK_SECONDARY, ls=":", lw=0.8)
    ax.text(1.62, 112, "photon sphere", fontsize=8, rotation=90, ha="right")
    ax.text(0.93, 130, "horizon", fontsize=8, rotation=90, ha="left")
    e = math.degrees(raymap.shadow_edge(2.0))
    ax.plot(1.0, e, "o", color=P.FALLING)
    P.annotate(ax, f"{e:.1f}° as you cross", xy=(1.0, e), xytext=(3.0, 70))
    ax.set_xlabel("distance from the centre (rₛ)")
    ax.set_ylabel("angular radius of the dark region (degrees)")
    ax.set_ylim(0, 180)
    ax.set_yticks([0, 30, 60, 90, 120, 150, 180])
    ax.invert_xaxis()
    ax.set_title("a. How much of the sky is black")
    ax.legend(loc="upper left", fontsize=8)
    ax = axs[1]
    for rr, col in ((2000.0, P.MUTED), (20.0, P.BLUE_150), (6.0, P.BLUE_300), (2.0, P.BLUE), (0.2, P.BLUE_550)):
        tab = raymap.RayTable(rr, 0.0)
        ax.semilogy(np.degrees(tab.chi), tab.g, color=col, lw=1.6, label=f"{rr / 2:g} rₛ")
    ax.axhline(1, color=P.AXIS, lw=0.8)
    ax.set_xlim(0, 180)
    ax.set_xlabel("viewing angle from straight ahead (degrees; 180 = straight back)")
    ax.set_ylabel("frequency shift g (>1: bluer, brighter)")
    ax.set_title("b. Colour shift of the sky, falling from far away")
    ax.legend(title="distance", fontsize=7.5)
    P.caption(fig, "Hovering observers (Synge 1966: sin ψ = √27 GM/c² √(1−2GM/rc²)/r) see the dark region "
              "grow to fill the whole sky at the horizon. A falling observer moves toward the light from below, which "
              "aberration compresses forward: the hole looks smaller, and the sky ahead is blueshifted (up to g = 3.9 at "
              "the horizon) while the sky behind is redshifted (g = 1/2 straight back at the horizon). Light arriving "
              "exactly sideways is never shifted (g = 1/(1 − β cos 90°) = 1, at every radius).")
    fig.tight_layout()
    return P.save(fig, OUT / "06_shadow_size_and_shift")


# ---------------------------------------------------------------------------
# 7-9. Views of the real sky
# ---------------------------------------------------------------------------

def _view(scene, r, kind="ahead", size=640, fov=150.0, v_rel=0.0):
    zhat, up = RD.default_frame()
    table = raymap.RayTable(r, v_rel)
    if kind == "ahead":
        cam = RD.Camera(size, size, forward=-zhat, up=up, fov_deg=fov)
    elif kind == "behind":
        cam = RD.Camera(size, size, forward=zhat, up=up, fov_deg=fov)
    elif kind == "whole":
        cam = RD.Camera(size, size, forward=-zhat, up=up, projection="equidistant")
    elif kind == "side":
        zh = RD.unit(zhat)
        cam = RD.Camera(size, size, forward=RD.unit(np.cross(zh, up)), up=zh, fov_deg=fov)
    return RD.tone_map(RD.render_view(scene, cam, table, zhat))


def views_sequence(scene, bhs):
    tu = bhs["sgra"]["gm_c3_s"]
    radii = [2000.0, 60.0, 12.0, 4.0, 2.0, 0.5]
    fig, axs = plt.subplots(2, 3, figsize=(15, 10.6), facecolor="black")
    for ax, r in zip(axs.flat, radii):
        ax.imshow(_view(scene, r, "ahead", 720))
        ax.axis("off")
        left = infall.rain_tau_left(r) * tu
        tl = f"{left / 86400:.1f} days" if left > 86400 * 2 else _fmt_time(left)
        title = f"{r / 2:g} rₛ" + ("  (the horizon)" if r == 2.0 else "") + f"   — {tl} left (Sgr A*)"
        ax.set_title(title, color="white", fontsize=12, loc="left")
    fig.suptitle("Falling toward Sgr A*, looking ahead (150° field, real Gaia/Hipparcos sky)", color="white",
                 fontsize=14, x=0.01, ha="left")
    fig.text(0.01, 0.005, "SIMULATION. The Milky Way (Gaia DR3) is bent into Einstein rings around the shadow; as you fall, "
             "the sky crowds forward into a bright, blueshifted ring. Sky as seen from the Sun; no gas around the hole.",
             color="#bbbbbb", fontsize=9)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    out = P.save(fig, OUT / "07_views_ahead", formats=("jpg",))
    return out


def hover_vs_fall(scene):
    fig, axs = plt.subplots(2, 3, figsize=(15, 10.4), facecolor="black")
    for j, r in enumerate((6.0, 3.0, 2.1)):
        axs[0, j].imshow(_view(scene, r, "whole", 620, v_rel=raymap.velocity_rel_rain(r, "static")))
        axs[1, j].imshow(_view(scene, r, "whole", 620, v_rel=0.0))
        axs[0, j].set_title(f"hovering at {r / 2:g} rₛ", color="white", loc="left")
        axs[1, j].set_title(f"falling past {r / 2:g} rₛ", color="white", loc="left")
        for ax in axs[:, j]:
            ax.axis("off")
    fig.suptitle("The whole sky (centre: toward the hole; rim: straight away): same place, different motion",
                 color="white", fontsize=14, x=0.01, ha="left")
    fig.text(0.01, 0.005, "A hovering observer's sky shrinks into a small, blueshifted window overhead (the rim here) as "
             "the horizon nears; a falling observer's does not (aberration). Hamilton & Polhemus (2010) make the same point.",
             color="#bbbbbb", fontsize=9)
    fig.tight_layout(rect=(0, 0.02, 1, 0.96))
    return P.save(fig, OUT / "08_hovering_vs_falling", formats=("jpg",))


def near_singularity(scene, bhs):
    tu = bhs["sgra"]["gm_c3_s"]
    fig, axs = plt.subplots(1, 3, figsize=(15, 5.4), facecolor="black")
    for ax, r in zip(axs, (1.0, 0.2, 0.04)):
        ax.imshow(_view(scene, r, "side", 640))
        ax.axis("off")
        ax.set_title(f"{r / 2:g} rₛ, {_fmt_time(infall.rain_tau_left(r) * tu)} left: looking sideways",
                     color="white", loc="left")
    fig.text(0.01, 0.02, "Up = away from the hole. The tide that stretches the faller lengthwise squeezes their view of the "
             "outside universe into a thin blueshifted band around the waist; above is dim and red, below is the dark "
             "region (Hamilton's \"band around your waist\").", color="#bbbbbb", fontsize=9)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    return P.save(fig, OUT / "09_near_the_singularity", formats=("jpg",))


# ---------------------------------------------------------------------------
# 10. Spinning black holes
# ---------------------------------------------------------------------------

def kerr_figure(bhs):
    fig = plt.figure(figsize=(15, 5.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.15])
    for j, key in enumerate(("gw250114", "cygx1")):
        b = bhs[key]
        a = float(b["spin"])
        ax = fig.add_subplot(gs[0, j])
        rp, rm = kerr.horizons(a)
        ax.add_patch(Circle((0, 0), 2.0, fc=P.WASH, ec=P.INK_SECONDARY, lw=0.8, ls=":"))
        ax.add_patch(Circle((0, 0), rp, fc="#111111", ec="none"))
        d = kerr.plunge_bl(a, 1.0, 0.0, 9.0)
        ax.plot(d["x"], d["y"], color=P.BLUE, lw=1.6, label="dropped with zero angular momentum")
        E, L = kerr.circular_EL(kerr.isco(a), a)
        q = kerr.plunge_bl(a, E, L, kerr.isco(a) - 0.1)
        ax.plot(q["x"], q["y"], color=P.ORANGE, lw=1.1, label="gas leaving the innermost stable orbit")
        ri = kerr.isco(a)
        th = np.linspace(0, 2 * np.pi, 300)
        ax.plot(ri * np.cos(th), ri * np.sin(th), color=P.ORANGE, lw=0.6, ls=":")
        ax.set_aspect("equal")
        lim = 4.6
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("x (GM/c²)")
        ax.set_ylabel("y (GM/c²)")
        turns = d["phi"][-1] / (2 * math.pi)
        ax.set_title(f"{'ab'[j]}. {b['name']}: spin a = {a:g}", loc="left")
        ax.text(0.02, 0.02, f"dropped body: swept round {turns:.1f} turns before r₊ + 0.002 GM/c²\n"
                f"ISCO {b['kerr']['r_isco_km']:.0f} km; ergosphere (shaded) to {b['kerr']['ergosphere_equator_km']:.0f} km",
                transform=ax.transAxes, fontsize=7.5)
        if j == 0:
            ax.legend(loc="upper left", fontsize=7.5)
    ax = fig.add_subplot(gs[0, 2])
    b = bhs["gw250114"]
    a = float(b["spin"])
    rp, _ = kerr.horizons(a)
    r = np.linspace(rp + 1e-4, 8, 600)
    lo, hi = kerr.allowed_omega(r, a)
    om = kerr.frame_dragging_omega(r, a)
    ax.fill_between(r, lo, hi, color=P.BLUE_150, alpha=0.6, lw=0, label="possible angular velocities")
    ax.plot(r, om, color=P.BLUE_550, lw=1.4, label="zero-angular-momentum observers")
    ax.axhline(0, color=P.INK, lw=0.8)
    ax.axvline(2.0, color=P.INK_SECONDARY, ls=":", lw=0.9)
    ax.axvline(rp, color=P.INK, ls="--", lw=0.9)
    ax.text(2.07, -0.3, "edge of the\nergosphere", fontsize=8)
    ax.text(rp + 0.06, 0.32, "horizon", fontsize=8, rotation=90)
    oh = a / (rp**2 + a * a)
    P.annotate(ax, f"at the horizon everything turns with it:\n{b['kerr']['omega_horizon_hz']:.0f} times a second for GW250114's remnant",
               xy=(rp, oh), xytext=(3.0, 0.33))
    ax.set_xlabel("distance r (GM/c², Boyer-Lindquist)")
    ax.set_ylabel("dφ/dt (c³/GM)")
    ax.set_ylim(-0.4, 0.45)
    ax.set_title("c. Inside the ergosphere nothing can stand still", loc="left")
    ax.legend(loc="lower right", fontsize=7.5)
    P.caption(fig, "Kerr equatorial geodesics in Boyer-Lindquist coordinates, whose angle is measured against the distant "
              "stars; drawn to just outside the horizon (r₊, black), which in these coordinates is reached only as "
              "t → ∞, the body co-rotating with the horizon. By its own clock it crosses after a finite twist. "
              "Spins: GW250114's remnant 0.68 (LVK 2025); Cyg X-1 > 0.9985 (Zhao et al. 2021, contested).")
    fig.tight_layout()
    return P.save(fig, OUT / "10_spinning_black_holes")


def poster_frame():
    """A still of the movie at the horizon crossing for the README."""
    from PIL import Image
    import movie
    fr = movie.schedule()
    idx = next(i for i, f in enumerate(fr) if f["act"] == "inside")
    src = C.FRAME_DIR / f"frame_{idx:05d}.png"
    if src.exists():
        Image.open(src).convert("RGB").resize((1280, 720), Image.LANCZOS).save(OUT / "movie_poster.jpg", quality=88)
        return [str(OUT / "movie_poster.jpg")]
    return []
