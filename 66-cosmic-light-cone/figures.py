"""Every figure, from data/results.json plus the expansion history (recomputed: < 1 s).

    python figures.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # works under python -I too

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import config as C  # noqa: E402
import plotstyle as ps  # noqa: E402
from cosmology import Universe  # noqa: E402

FIG = C.FIG_DIR


def _gyr_grid(t_max: float, n: int = 3000) -> np.ndarray:
    return np.unique(np.concatenate([np.geomspace(1e-6, 0.6, 600), np.linspace(0.6, t_max, n)]))


def spacetime_curves(u: Universe, t_max: float, *, conformal: bool = False) -> dict:
    """Light cone, horizons and Hubble sphere, in proper and comoving distance (Gly).
    With ``conformal`` the grid runs in conformal time all the way to t = infinity."""
    if conformal:
        eta_end = float(u.eta_of_a(u.a_max * 0.999))
        eta = np.unique(np.concatenate([np.geomspace(1e-5, 1.0, 300), np.linspace(1.0, eta_end, 3000)]))
        a = u.a_of_eta(eta)
        t = u.t_of_eta(eta)
    else:
        t = _gyr_grid(t_max)
        a = u.a_of_t(t)
    eta = u.eta_of_t(t) if not conformal else eta
    H = u.H_of_a(a)
    past = t <= u.age
    out = {"t": t, "a": a, "eta": eta,
           "lc_chi": np.where(past, u.eta0 - eta, np.nan),
           "eh_chi": u.eta_inf - eta, "ph_chi": eta, "hs_chi": 1.0 / (a * H)}
    for k in ("lc", "eh", "ph", "hs"):
        out[f"{k}_D"] = out[f"{k}_chi"] * a
    return out


def _scale_factor_axis(ax, u: Universe, kind: str, ticks):
    """Label the right-hand side with the scale factor a (the same time axis, relabelled)."""
    if kind == "t":
        ys = [float(u.t_of_a(a)) for a in ticks]
    else:
        ys = [float(u.eta_of_a(a)) if np.isfinite(a) else u.eta_inf for a in ticks]
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks(ys)
    ax2.set_yticklabels(["∞" if not np.isfinite(a) else f"{a:g}" for a in ticks])
    ax2.grid(False)
    ax2.spines["right"].set_visible(True)
    ax2.tick_params(colors=ps.MUTED, labelsize=7.5)
    ax2.set_ylabel("scale factor a  (1 = today)", color=ps.INK_SECONDARY, fontsize=8)
    return ax2


def _mirror(ax, x, y, **kw):
    label = kw.pop("label", None)
    ax.plot(x, y, label=label, **kw)
    ax.plot(-np.asarray(x), y, **kw)


# ---------------------------------------------------------------------------
# Figure 1: the three spacetime diagrams (Davis & Lineweaver 2004, Fig. 1, for Planck 2018)
# ---------------------------------------------------------------------------


def fig_spacetime(u: Universe, results: dict, path, *, title: str, landmarks: bool = True,
                  t_max: float = 28.0, xlim: float = 65.0, worldline_z=C.WORLDLINE_Z,
                  extra_chi=(60.0,)):
    cur_t = spacetime_curves(u, t_max)
    cur_c = spacetime_curves(u, t_max, conformal=True)
    t = cur_t["t"]
    fig, axes = plt.subplots(3, 1, figsize=(8.6, 13.6), gridspec_kw={"hspace": 0.32})
    fig.subplots_adjust(bottom=0.1)
    marks = results["cone"]["landmarks"] if landmarks else []
    show = [m for m in marks if m["name"] in ("CMB (last scattering)", "MoM-z14", "GN-z11",
                                             "J0313-1806", "SN 2023adsy", "GW190403_051519")]
    for k, ax in enumerate(axes):
        proper = k == 0
        conformal = k == 2
        cur = cur_c if conformal else cur_t
        if conformal:
            y = cur["eta"]
            ymax = u.eta_inf if np.isfinite(u.eta_inf) else 70.0
            y_now = u.eta0
        else:
            y = t
            ymax = t_max
            y_now = u.age
        key = "D" if proper else "chi"
        # Superluminal region: beyond the Hubble sphere (recession faster than light).
        hs = cur[f"hs_{key}"]
        ax.fill_betweenx(y, hs, xlim * 1.2, color=ps.WASH_BLUE, lw=0, zorder=0)
        ax.fill_betweenx(y, -xlim * 1.2, -hs, color=ps.WASH_BLUE, lw=0, zorder=0)
        # Worldlines of comoving galaxies (labelled by today's redshift of their light).
        for z in worldline_z:
            chi = float(u.chi_of_z(z))
            xw = cur["a"] * chi if proper else np.full_like(y, chi)
            _mirror(ax, xw, y, color=ps.C_WORLD, lw=0.7, zorder=1)
            # label where the line leaves the panel (top or side)
            inside = np.abs(xw) < xlim
            j = np.flatnonzero(inside)[-1]
            xl, yl = xw[j], y[j]
            lab = f"z = {z:g}"
            if yl > ymax * 0.97:
                ax.text(xl, ymax * 0.985, lab, ha="center", va="top", fontsize=7, color=ps.INK_SECONDARY,
                        bbox={"fc": ps.SURFACE, "ec": "none", "pad": 0.6})
            else:
                ax.text(xlim * 0.985, yl, lab, ha="right", va="bottom", fontsize=7,
                        color=ps.INK_SECONDARY, bbox={"fc": ps.SURFACE, "ec": "none", "pad": 0.6})
        for chi in extra_chi:  # a galaxy we cannot see yet
            xw = cur["a"] * chi if proper else np.full_like(y, chi)
            _mirror(ax, xw, y, color=ps.C_WORLD, lw=0.5, ls=(0, (1, 2)), zorder=1)
        _mirror(ax, [0], [0], color=ps.INK)
        ax.axvline(0, color=ps.INK_SECONDARY, lw=0.8, zorder=2)
        # Horizons and the light cone.
        _mirror(ax, cur[f"ph_{key}"], y, color=ps.C_PARTICLE, lw=1.6, zorder=3)
        if np.isfinite(u.eta_inf):
            _mirror(ax, cur[f"eh_{key}"], y, color=ps.C_EVENT, lw=1.6, zorder=3)
        _mirror(ax, hs, y, color=ps.C_HUBBLE, lw=1.6, zorder=3)
        _mirror(ax, cur[f"lc_{key}"], y, color=ps.C_CONE, lw=2.0, zorder=4)
        ax.axhline(y_now, color=ps.INK_SECONDARY, lw=0.7, zorder=2)
        ax.text(-xlim * 0.985, y_now, f"now ({u.age:.2f} Gyr)", ha="left", va="bottom",
                fontsize=7.5, color=ps.INK_SECONDARY)
        # Real objects where their light set out (right half only).
        for m in show:
            ye = float(u.eta_of_z(m["z"])) if conformal else m["t_emit_gyr"]
            xe = m["d_emit_gly"] if proper else m["chi_gly"]
            ax.plot(xe, ye, "o", ms=4.2, mfc=ps.INK, mec=ps.SURFACE, mew=1.0, zorder=6)
            if conformal:
                name = m["name"].replace("CMB (last scattering)", "CMB").replace("_051519", "")
                ax.annotate(name, xy=(xe, ye), xytext=(7, 0), textcoords="offset points",
                            fontsize=6.8, color=ps.INK_SECONDARY, va="center", zorder=7,
                            bbox={"fc": ps.SURFACE, "ec": "none", "pad": 0.4, "alpha": 0.85})
        ax.set_xlim(-xlim, xlim)
        ax.set_ylim(0, ymax)
        ax.grid(False)
        if proper:
            ax.set_xlabel("proper distance D (billion light-years) — the distance a ruler would measure at that moment")
            ax.set_ylabel("cosmic time t (billion years)")
            ax.set_title("a  Proper distance: how far away things were at each moment")
            _scale_factor_axis(ax, u, "t", [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5])
            j = np.argmin(np.abs(t - 8.0))
            ps.annotate(ax, "our past light cone\n(everything we see today)", (cur["lc_D"][j], 8.0),
                        (12, 6.5))
            j = np.argmin(np.abs(t - 17.0))
            ps.annotate(ax, "Hubble sphere\n(recession speed = c)", (cur["hs_D"][j], 17.0), (24, 16.0))
            j = np.argmin(np.abs(t - 24.0))
            ps.annotate(ax, "event horizon", (cur["eh_D"][j], 24.0), (27, 22.8))
            j = np.argmin(np.abs(t - 9.0))
            ps.annotate(ax, "particle horizon", (cur["ph_D"][j], 9.0), (50, 7.5))
            ax.text(-44, 3.0, "receding faster than light", fontsize=7.5, color=ps.BLUE, ha="center")
        elif not conformal:
            ax.set_xlabel("comoving distance χ (billion light-years) — proper distance today")
            ax.set_ylabel("cosmic time t (billion years)")
            ax.set_title("b  Comoving distance: galaxies stand still, the grid expands")
            _scale_factor_axis(ax, u, "t", [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5])
            j = np.argmin(np.abs(t - 3.0))
            ps.annotate(ax, "our past light cone", (-cur["lc_chi"][j], 3.0), (-38, 9.0))
            j = np.argmin(np.abs(t - 7.6))
            ps.annotate(ax, "Hubble sphere at its widest:\nthe expansion starts to accelerate",
                        (-cur["hs_chi"][j], 7.6), (-31.5, 17.5))
            j = np.argmin(np.abs(t - 22.0))
            ps.annotate(ax, "particle horizon", (-cur["ph_chi"][j], 22.0), (-62, 25.5))
        else:
            ax.set_xlabel("comoving distance χ (billion light-years)")
            ax.set_ylabel("conformal time η (billion years): light moves at 45°")
            ax.set_title("c  Conformal time: light rays are straight lines; infinite time fits on the page")
            _scale_factor_axis(ax, u, "eta", [0.001, 0.01, 0.1, 0.3, 0.6, 1.0, 2.0, 5.0, np.inf])
            ps.annotate(ax, "our past light cone", (-20.0, u.eta0 - 20.0), (-44, 34))
            ps.annotate(ax, "event horizon\n(our light cone at t = ∞)", (-16.0, u.eta_inf - 16.0),
                        (-46, 52))
            ps.annotate(ax, "particle horizon", (-50.0, 50.0), (-63, 58))
            yy = cur["eta"]
            j = np.argmin(np.abs(yy - 26.0))
            ps.annotate(ax, "Hubble sphere", (-cur["hs_chi"][j], 26.0), (-32, 18))
        ax.tick_params(axis="both", which="both", length=3)
    fig.suptitle(title, x=0.06, ha="left", y=0.925, fontsize=12)
    from matplotlib.lines import Line2D

    handles = [Line2D([], [], color=ps.C_CONE, lw=2.0, label="our past light cone"),
               Line2D([], [], color=ps.C_HUBBLE, lw=1.6, label="Hubble sphere (v = c)"),
               Line2D([], [], color=ps.C_EVENT, lw=1.6, label="event horizon"),
               Line2D([], [], color=ps.C_PARTICLE, lw=1.6, label="particle horizon"),
               Line2D([], [], color=ps.C_WORLD, lw=0.7, label="galaxies (worldlines)")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.06, 0.912), ncol=5,
               fontsize=7.6, handlelength=1.8, columnspacing=1.2)
    ps.caption(fig, "Planck 2018 cosmology (H0 = 67.66 km/s/Mpc, Ωm = 0.311, flat). Grey lines are "
                    "galaxies labelled by the redshift we see them at today; the dotted grey line is a galaxy "
                    "60 billion light-years away that no light has yet reached us from. Shaded: regions "
                    "receding faster than light. Dots: the CMB, MoM-z14, GN-z11, quasar J0313-1806, "
                    "SN 2023adsy and the farthest gravitational-wave merger, where their light set out. "
                    "After Davis & Lineweaver (2004), Fig. 1.", y=0.062)
    return ps.save(fig, path)


# ---------------------------------------------------------------------------
# Figure 2: where and when the light we see set out
# ---------------------------------------------------------------------------

POP_ORDER = [("sdss_galaxies", "SDSS galaxies"), ("gw_events", "Gravitational-wave mergers"),
             ("supernovae", "Type Ia supernovae"), ("quasars", "SDSS quasars"),
             ("jades", "JWST galaxies (JADES)"), ("highz_quasars", "Quasars beyond z = 5.6")]


def _t_of_z_axis(ax, u: Universe, zs=(0.1, 0.5, 1, 2, 3, 5, 10)):
    top = ax.secondary_xaxis("top")
    top.set_xticks([float(u.t_of_z(z)) for z in zs])
    top.set_xticklabels([f"{z:g}" for z in zs])
    top.tick_params(colors=ps.MUTED, labelsize=7.5, length=3)
    top.set_xlabel("redshift we see today", color=ps.INK_SECONDARY, fontsize=8)
    return top


def fig_objects(u: Universe, results: dict, arrays, path):
    cone = results["cone"]
    marks = {m["name"]: m for m in cone["landmarks"]}
    t_hs = results["headline"]["hubble_sphere_now"]["t_emit_gyr"]
    t_turn = results["headline"]["light_cone_turn"]["t_gyr"]
    t_eh = results["headline"]["event_horizon_now"]["t_emit_gyr"]
    fig = plt.figure(figsize=(8.6, 10.4))
    gs = fig.add_gridspec(1 + len(POP_ORDER), 1, height_ratios=[3.2] + [0.62] * len(POP_ORDER),
                          hspace=0.12, top=0.9, bottom=0.115, left=0.24, right=0.8)
    ax = fig.add_subplot(gs[0])
    t = np.linspace(1e-6, u.age, 3000)
    a, eta, H = u.a_of_t(t), u.eta_of_t(t), u.H_of_t(t)
    d_lc = a * (u.eta0 - eta)
    d_hs = 1.0 / H
    ax.fill_between(t, d_hs, 9, color=ps.WASH_BLUE, lw=0, zorder=0)
    ax.plot(t, d_hs, color=ps.C_HUBBLE, lw=1.6, zorder=2)
    ax.plot(t, d_lc, color=ps.C_CONE, lw=2.0, zorder=3)
    ps.annotate(ax, "Hubble sphere (recession = c):\nabove it, faster than light",
                (5.52, float(np.interp(5.52, t, d_hs))), (6.3, 7.45))
    ps.annotate(ax, "our past light cone: where every photon\nwe receive today was, at each moment",
                (9.0, float(np.interp(9.0, t, d_lc))), (9.6, 4.6))
    ax.plot(t_turn, results["headline"]["light_cone_turn"]["width_gly"], "o", ms=5, mfc=ps.SURFACE,
            mec=ps.INK, mew=1.2, zorder=5)
    ps.annotate(ax, f"t = {t_turn:.2f} Gyr: the photons stop receding and\nstart to approach us, "
                    f"{results['headline']['light_cone_turn']['width_gly']:.2f} billion ly away",
                (t_turn, results["headline"]["light_cone_turn"]["width_gly"]), (6.6, 6.35))
    label_pos = {"CMB (last scattering)": (0.55, 0.9), "MoM-z14": (1.0, 1.25),
                 "GN-z11": (1.15, 1.95), "J0313-1806": (1.3, 2.65), "GRB 250314A + SN": (1.45, 3.35),
                 "SN 2023adsy": (1.15, 5.55), "GW190403_051519": (6.3, 3.4),
                 "GW150914": (10.6, 0.6), "GW170817": (12.0, 1.95)}
    for name, (xt, yt) in label_pos.items():
        m = marks.get(name)
        if m is None:
            continue
        ax.plot(m["t_emit_gyr"], m["d_emit_gly"], "o", ms=4.4, mfc=ps.INK, mec=ps.SURFACE, mew=1.0,
                zorder=6)
        nice = {"CMB (last scattering)": "CMB: 42 million ly away",
                "GW190403_051519": "farthest GW merger (GW190403)",
                "GRB 250314A + SN": "GRB 250314A + supernova"}.get(name, name)
        an = ps.annotate(ax, nice, (m["t_emit_gyr"], m["d_emit_gly"]), (xt, yt), fontsize=7)
        an.set_bbox({"fc": ps.SURFACE, "ec": "none", "pad": 0.5, "alpha": 0.85})
    for tv in (t_eh, t_turn, t_hs):
        ax.axvline(tv, color=ps.MUTED, lw=0.7, zorder=1)
    ax.set_xlim(0, u.age)
    ax.set_ylim(0, 8.2)
    ax.set_ylabel("proper distance when the light\nset out (billion light-years)")
    ax.tick_params(labelbottom=False)
    _t_of_z_axis(ax, u)
    ax.set_title("Where and when the light we see today set out", pad=34)
    edges = arrays["t_edges"]
    centers = 0.5 * (edges[1:] + edges[:-1])
    pops = cone["populations"]
    for i, (key, label) in enumerate(POP_ORDER):
        axl = fig.add_subplot(gs[i + 1], sharex=ax)
        h = arrays[f"{key}_t_hist"].astype(float)
        h = h / h.max()
        axl.fill_between(centers, 0, h, step="mid", color=ps.lighten(ps.BLUE, 0.55), lw=0)
        axl.step(centers, h, where="mid", color=ps.BLUE, lw=0.9)
        for tv in (t_eh, t_turn, t_hs):
            axl.axvline(tv, color=ps.MUTED, lw=0.7)
        axl.set_ylim(0, 1.15)
        axl.set_yticks([])
        axl.grid(False)
        axl.spines["left"].set_visible(False)
        p = pops[key]
        axl.text(-0.02, 0.5, f"{label}\n{p['n']:,}", transform=axl.transAxes, ha="right", va="center",
                 fontsize=7.6, color=ps.INK)
        if p["n_superluminal_now"] and p["n_superluminal_now"] < 0.01 * p["n"]:
            txt = (f"{p['n_superluminal_now']} recede faster than light now\n"
                   f"{p['n_beyond_event_horizon']} beyond our event horizon")
        elif p["n_superluminal_now"]:
            txt = (f"{100 * p['n_superluminal_now'] / p['n']:.0f}% recede faster than light now\n"
                   f"{100 * (1 - p['frac_reachable']):.0f}% beyond our event horizon")
        else:
            txt = "all recede slower than light\nall still reachable"
        axl.text(1.02, 0.5, txt, transform=axl.transAxes, ha="left", va="center", fontsize=6.9,
                 color=ps.INK_SECONDARY)
        if i < len(POP_ORDER) - 1:
            axl.tick_params(labelbottom=False)
    axl.set_xlabel("cosmic time when the light set out (billion years after the Big Bang)")
    ax.text(t_eh - 0.08, 0.2, "← earlier: beyond our\nevent horizon now", ha="right", va="bottom",
            fontsize=6.8, color=ps.INK_SECONDARY)
    ax.text(t_hs + 0.08, 4.0, "later: receding slower\nthan light today →", ha="left", va="center",
            fontsize=6.8, color=ps.INK_SECONDARY)
    ps.caption(fig, "Top: proper distance between us and each photon now arriving, against the cosmic "
                    "time it was there (Planck 2018). Every object we see lies on this curve, at the moment "
                    "its light set out. Lanes: when the light of each catalogue set out (counts per 0.1 Gyr, "
                    "each lane scaled to its own peak). Vertical lines at "
                    f"{t_eh:.2f}, {t_turn:.2f} and {t_hs:.2f} Gyr: light set out earlier than these comes from "
                    "galaxies now beyond our event horizon (z > 1.85), from galaxies receding faster than "
                    "light when it left (z > 1.59), and from galaxies receding faster than light today "
                    "(z > 1.48).", y=0.06)
    return ps.save(fig, path)


# ---------------------------------------------------------------------------
# Figure 3: recession speeds, now and when the light left
# ---------------------------------------------------------------------------


def fig_speeds(u: Universe, results: dict, path):
    from analysis import beta_sr

    desi = results["desi"]
    sn = results["sn_velocities"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 5.2), gridspec_kw={"wspace": 0.22,
                                                                          "width_ratios": [1.05, 1]})
    z = np.geomspace(1e-3, 3000, 1500)
    vn, vt = u.v_now(z), u.v_then(z)
    for ax in (ax1, ax2):
        ax.axhline(1.0, color=ps.INK, lw=0.9, zorder=1)
        ax.plot(z, vn, color=ps.C_NOW, lw=2.0, zorder=3)
        ax.plot(z, vt, color=ps.C_THEN, lw=2.0, zorder=3)
        ax.plot(z, beta_sr(z), color=ps.MUTED, lw=1.3, ls=(0, (5, 3)), zorder=2)
        ax.plot(z, z, color=ps.MUTED, lw=1.0, ls=(0, (1, 2)), zorder=2)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(1e-2, 2000)
    ax1.set_ylim(1e-2, 150)
    ax1.set_xlabel("redshift z")
    ax1.set_ylabel("recession speed / speed of light")
    ax1.set_title("a  Galaxies we see recede faster than light, now and then")
    cmb = results["headline"]["cmb"]
    ax1.plot([cmb["z"], cmb["z"]], [cmb["v_now_c"], cmb["v_then_c"]], "o", ms=5, mfc=ps.INK,
             mec=ps.SURFACE, mew=1, zorder=6)
    ps.annotate(ax1, f"CMB: {cmb['v_then_c']:.0f}c when it set out", (cmb["z"], cmb["v_then_c"]),
                (25, 75))
    ps.annotate(ax1, f"CMB: {cmb['v_now_c']:.2f}c today", (cmb["z"], cmb["v_now_c"]), (60, 1.45))
    mom = [m for m in results["cone"]["landmarks"] if m["name"] == "MoM-z14"][0]
    ax1.plot([mom["z"], mom["z"]], [mom["v_now_c"], mom["v_then_c"]], "o", ms=4.5, mfc=ps.INK,
             mec=ps.SURFACE, mew=1, zorder=6)
    ps.annotate(ax1, f"MoM-z14: {mom['v_then_c']:.1f}c then,\n{mom['v_now_c']:.2f}c now",
                (mom["z"], mom["v_then_c"]), (2.2, 18))
    ax1.text(4.2e-2, 1.12, "speed of light", fontsize=7.5, color=ps.INK, va="bottom")
    ax1.text(14, 0.62, "special-relativistic Doppler reading:\nnever reaches c", fontsize=7,
             color=ps.INK_SECONDARY, ha="left", va="top")
    ax1.text(0.03, 0.06, "v = cz", fontsize=7, color=ps.INK_SECONDARY, rotation=33)
    ax1.text(250, 4.1, "today", fontsize=7.8, color=ps.INK, ha="center", va="bottom")
    ax1.text(55, 40, "when the light set out", fontsize=7.8, color=ps.INK, ha="right", va="bottom")
    # zoom with data
    ax2.set_xlim(0, 2.6)
    ax2.set_ylim(0, 1.62)
    ax2.set_xlabel("redshift z")
    ax2.set_title("b  Measured: DESI (when the light left) and supernovae (today)")
    rows = desi["rows"]
    ax2.errorbar([r["z"] for r in rows], [r["v_then_c"] for r in rows],
                 yerr=[r["v_then_sigma"] for r in rows], fmt="o", ms=6, mfc=ps.C_THEN, mec=ps.SURFACE,
                 mew=1.2, ecolor=ps.C_THEN, elinewidth=1.2, capsize=0, zorder=6)
    for r in rows:
        ax2.annotate(r["tracer"].replace("Ly-alpha", "Lyα"), (r["z"], r["v_then_c"]), xytext=(-6, 8),
                     textcoords="offset points", fontsize=6.6, color=ps.INK_SECONDARY, ha="right")
    b = sn["bins"]
    ax2.errorbar([x["z_eff"] for x in b], [x["v_now_c"] for x in b], yerr=[x["v_now_sigma"] for x in b],
                 fmt="s", ms=5.5, mfc=ps.C_NOW, mec=ps.SURFACE, mew=1.2, ecolor=ps.C_NOW, elinewidth=1.2,
                 capsize=0, zorder=6)
    s1 = sn["beyond_hubble_sphere"]
    ax2.errorbar([x["z"] for x in s1], [x["v_now_c"] for x in s1], yerr=[x["v_now_sigma"] for x in s1],
                 fmt="s", ms=3.5, mfc=ps.SURFACE, mec=ps.C_NOW, mew=0.9, ecolor=ps.lighten(ps.C_NOW, 0.45),
                 elinewidth=0.8, capsize=0, zorder=5)
    last = b[-1]
    ps.annotate(ax2, f"the 7 supernovae beyond z = 1.48, together:\n{last['v_now_c']:.3f} ± "
                     f"{last['v_now_sigma']:.3f} c today ({last['sigma_from_c']:.1f}σ above c)",
                (last["z_eff"], last["v_now_c"] - last["v_now_sigma"]), (1.93, 0.905), va="center")
    ly = rows[-1]
    ps.annotate(ax2, f"Lyman-α forest, z = 2.33: {ly['v_then_c']:.3f} ± {ly['v_then_sigma']:.3f} c\n"
                     f"when the light left ({ly['sigma_from_c']:.0f}σ above c)",
                (ly["z"] - 0.02, ly["v_then_c"] + 0.02), (1.08, 1.555), va="top")
    cr = desi["crossing"]
    ax2.errorbar([cr["z_cross_median"]], [1.0], xerr=[[cr["z_cross_median"] - cr["z_cross_lo68"]],
                                                      [cr["z_cross_hi68"] - cr["z_cross_median"]]],
                 fmt="D", ms=5, mfc=ps.SURFACE, mec=ps.C_THEN, mew=1.4, ecolor=ps.C_THEN, elinewidth=1.4,
                 capsize=3, zorder=7)
    ps.annotate(ax2, f"DESI alone: v = c at emission for\nz = {cr['z_cross_median']:.2f} "
                     f"(+{cr['z_cross_hi68'] - cr['z_cross_median']:.2f}/−"
                     f"{cr['z_cross_median'] - cr['z_cross_lo68']:.2f}); Planck: 1.59",
                (cr["z_cross_median"], 0.985), (1.62, 0.5), va="center")
    from matplotlib.lines import Line2D

    handles = [Line2D([], [], color=ps.C_NOW, lw=2, label="today (Planck 2018)"),
               Line2D([], [], color=ps.C_THEN, lw=2, label="when the light left (Planck 2018)"),
               Line2D([], [], color=ps.C_THEN, marker="o", ls="", ms=6, label="DESI DR2 BAO: D_M / ((1+z) D_H)"),
               Line2D([], [], color=ps.C_NOW, marker="s", ls="", ms=5.5, label="Pantheon+ supernovae (bins; 7 singles)"),
               Line2D([], [], color=ps.MUTED, lw=1.3, ls=(0, (5, 3)), label="special-relativistic Doppler"),
               Line2D([], [], color=ps.MUTED, lw=1.0, ls=(0, (1, 2)), label="v = cz")]
    ax2.legend(handles=handles, loc="upper left", fontsize=7.0, bbox_to_anchor=(0.0, 0.93))
    ps.caption(fig, "Recession speed = H × proper distance, the rate at which the proper distance to a "
                    "galaxy grows (Hubble's law holds at every distance in an expanding universe). Nothing "
                    "outruns light locally; space between us and the galaxy is what grows. DESI's ratio "
                    "D_M/((1+z) D_H) equals H(z) D_A / c, the recession speed at emission in units of c, "
                    "without the size of the BAO ruler or H0. Supernovae measure H0 D_L, calibrated on the "
                    "nearby Hubble flow, which gives today's speed H0 D_M / c. Both assume a spatially flat "
                    "universe.", y=-0.02)
    return ps.save(fig, path)


# ---------------------------------------------------------------------------
# Figure 4: the angular-size turnover, measured by DESI
# ---------------------------------------------------------------------------


def fig_angular(u: Universe, results: dict, path):
    desi = results["desi"]
    rd = desi["r_drag_mpc"] * C.GLY_PER_MPC  # Gly
    rows = desi["rows"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 5.0), gridspec_kw={"wspace": 0.24})
    z = np.geomspace(0.01, 1500, 1500)
    mpc = C.GLY_PER_MPC
    arcmin = 180 / np.pi * 60
    for uu, lab, kw in [(u, "Planck 2018", dict(color=ps.INK, lw=2.0)),
                        (Universe(C.EDS), "matter only (Einstein–de Sitter)", dict(color=ps.MUTED, lw=1.2,
                                                                                 ls=(0, (5, 3)))),
                        (Universe(C.MILNE), "empty (Milne): no turnover", dict(color=ps.MUTED, lw=1.0,
                                                                               ls=(0, (1, 2))))]:
        ax1.plot(z, mpc / uu.da_of_z(z) * arcmin, label=lab, **kw)
    zmin = results["headline"]["summary"]["z_min_angular_size"]
    th_min = mpc / float(u.da_of_z(zmin)) * arcmin
    ax1.plot(zmin, th_min, "o", ms=5.5, mfc=ps.SURFACE, mec=ps.INK, mew=1.3, zorder=6)
    ps.annotate(ax1, f"smallest at z = {zmin:.2f}: {th_min:.2f}′\n(farther objects look bigger again)",
                (zmin, th_min * 0.97), (12.0, 1.13))
    ax1.errorbar([r["z"] for r in rows], [r["theta_1mpc_arcmin"] for r in rows],
                 yerr=[r["theta_1mpc_arcmin_sigma"] for r in rows], fmt="o", ms=5.5, mfc=ps.BLUE,
                 mec=ps.SURFACE, mew=1.1, ecolor=ps.BLUE, elinewidth=1.1, label="DESI DR2 (r_d from CAMB)",
                 zorder=7)
    zs = results["headline"]["z_star"]
    th_cmb = mpc / float(u.da_of_z(zs)) * arcmin
    ax1.plot(zs, th_cmb, "o", ms=5, mfc=ps.INK, mec=ps.SURFACE, mew=1, zorder=6)
    ps.annotate(ax1, f"at the CMB a 1 Mpc object\nspans {th_cmb / 60:.1f}°", (zs, th_cmb), (30, 2.9))
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(0.01, 1500)
    ax1.set_ylim(0.95, 400)
    ax1.set_xlabel("redshift z")
    ax1.set_ylabel("angle spanned by an object 1 Mpc across (arcminutes)")
    ax1.set_title("a  How big a fixed-size object looks")
    ax1.legend(loc="upper center", fontsize=7.4, bbox_to_anchor=(0.42, 1.0))
    # (b) D_A / r_d with DESI, and DESI's measured slope at each redshift
    zz = np.linspace(0.01, 2.7, 600)
    ax2.plot(zz, u.da_of_z(zz) / rd, color=ps.INK, lw=2.0, label="Planck 2018 (r_d = %.2f Mpc)" %
             desi["r_drag_mpc"])
    for r in rows:
        zc, y0 = r["z"], r["DA_over_rd"]
        slope = r["DH_over_rd"] / (1 + zc) * (1 - r["v_then_c"])     # dD_A/dz in r_d units
        s_sig = r["DH_over_rd"] / (1 + zc) * r["v_then_sigma"]
        dz = 0.085
        xs = np.array([zc - dz, zc + dz])
        ax2.fill_between(xs, y0 + (slope - s_sig) * (xs - zc), y0 + (slope + s_sig) * (xs - zc),
                         color=ps.lighten(ps.BLUE, 0.7), lw=0, zorder=4)
        ax2.plot(xs, y0 + slope * (xs - zc), color=ps.BLUE, lw=1.8, zorder=5)
    ax2.errorbar([r["z"] for r in rows], [r["DA_over_rd"] for r in rows],
                 yerr=[r["DA_over_rd_sigma"] for r in rows], fmt="o", ms=5.5, mfc=ps.BLUE, mec=ps.SURFACE,
                 mew=1.1, ecolor=ps.BLUE, elinewidth=1.1, zorder=7, label="DESI DR2: D_A / r_d and its slope")
    ly = rows[-1]
    ps.annotate(ax2, f"z = 2.33: D_A shrinks with redshift\n(slope from D_M and D_H: "
                     f"{ly['sigma_from_c']:.0f}σ below zero)", (ly["z"], ly["DA_over_rd"] - 0.25), (1.55, 9.0))
    ax2.axvline(zmin, color=ps.MUTED, lw=0.7)
    ax2.text(zmin + 0.03, 13.35, f"peak z = {zmin:.2f}", fontsize=7, color=ps.INK_SECONDARY, va="top")
    ax2.set_xlim(0, 2.7)
    ax2.set_ylim(6, 13.6)
    ax2.set_xlabel("redshift z")
    ax2.set_ylabel("angular-diameter distance D_A / r_d")
    ax2.set_title("b  DESI measures the turnover, and its slope, directly")
    ax2.legend(loc="lower right", fontsize=7.4)
    ps.caption(fig, "In flat space D_A is also the proper distance to the source when the light left, so "
                    "the turnover is the widest point of our teardrop-shaped light cone. Short blue lines: "
                    "the slope dD_A/dz = D_H (1 − v_then/c)/(1+z) measured by DESI at each redshift from its "
                    "D_M and D_H (pale wedge ±1σ), independent of the ruler's size. The 1 Mpc angles in "
                    "panel a use the Planck sound horizon r_d = 147.21 Mpc to turn D_A/r_d into megaparsecs.",
               y=-0.02)
    return ps.save(fig, path)


# ---------------------------------------------------------------------------
# Figure 5: horizons over time, and what we can still reach
# ---------------------------------------------------------------------------


def fig_reach(u: Universe, results: dict, path):
    reach = results["reachability"]
    pops = results["cone"]["populations"]
    fig = plt.figure(figsize=(12.0, 5.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1], height_ratios=[1, 1], wspace=0.28, hspace=0.55)
    ax = fig.add_subplot(gs[:, 0])
    t = np.geomspace(1e-4, 300, 3000)
    eta = u.eta_of_t(t)
    a, H = u.a_of_t(t), u.H_of_t(t)
    ax.plot(t, eta, color=ps.C_PARTICLE, lw=1.8, label="particle horizon: how far light has got since the Big Bang")
    ax.plot(t, u.eta_inf - eta, color=ps.C_EVENT, lw=1.8, label="event horizon: how far light sent now will ever get")
    ax.plot(t, 1 / (a * H), color=ps.C_HUBBLE, lw=1.6, label="Hubble sphere (recession = c)")
    ax.axhline(u.eta_inf, color=ps.MUTED, lw=0.7)
    ax.text(1.2e-4, u.eta_inf + 0.8, f"{u.eta_inf:.1f} billion ly: nothing farther will ever be seen",
            fontsize=7, color=ps.INK_SECONDARY)
    marks = {m["name"]: m for m in results["cone"]["landmarks"]}
    box = ["Light we will ever receive from them left", "no later than these ages:"]
    for name, lab in (("CMB (last scattering)", "gas that emitted the CMB"), ("MoM-z14", "MoM-z14"),
                      ("J0313-1806", "quasar J0313-1806")):
        m = marks[name]
        chi, tl = m["chi_gly"], m["t_last_seen_gyr"]
        ax.plot([m["t_emit_gyr"], tl], [chi, chi], color=ps.INK_SECONDARY, lw=0.9)
        ax.plot(m["t_emit_gyr"], chi, "|", ms=7, color=ps.INK_SECONDARY)
        ax.plot(tl, chi, "o", ms=4.5, mfc=ps.INK, mec=ps.SURFACE, mew=1, zorder=6)
        if m["t_emit_gyr"] > 0.01:
            ax.text(m["t_emit_gyr"] / 1.3, chi, lab, ha="right", va="center", fontsize=6.9,
                    color=ps.INK_SECONDARY)
        else:
            ax.text(m["t_emit_gyr"] * 1.2, chi + 0.9, lab, ha="left", va="bottom", fontsize=6.9,
                    color=ps.INK_SECONDARY)
        box.append(f"  {lab}: {tl:.2f} Gyr")
    ax.text(25, 30, "\n".join(box), fontsize=6.9, color=ps.INK_SECONDARY, va="center", ha="left",
            bbox={"fc": ps.SURFACE, "ec": ps.GRID, "pad": 4})
    ax.axvline(u.age, color=ps.INK_SECONDARY, lw=0.7)
    ax.text(u.age * 1.05, 2, "now", fontsize=7.5, color=ps.INK_SECONDARY)
    eh_now = u.eta_inf - u.eta0
    ax.plot(u.age, eh_now, "o", ms=5, mfc=ps.SURFACE, mec=ps.C_EVENT, mew=1.4, zorder=6)
    ax.annotate(f"{eh_now:.1f} billion ly today\n(galaxies at z = 1.85)", (u.age, eh_now), xytext=(9, 2),
                textcoords="offset points", fontsize=7, color=ps.INK_SECONDARY)
    ax.plot(u.age, u.eta0, "o", ms=5, mfc=ps.SURFACE, mec=ps.C_PARTICLE, mew=1.4, zorder=6)
    ax.annotate(f"{u.eta0:.1f} billion ly today", (u.age, u.eta0), xytext=(9, -10),
                textcoords="offset points", fontsize=7, color=ps.INK_SECONDARY)
    ax.set_xscale("log")
    ax.set_xlim(1e-4, 300)
    ax.set_ylim(0, 68)
    ax.set_xlabel("cosmic time (billion years)")
    ax.set_ylabel("comoving distance (billion light-years, measured today)")
    ax.set_title("a  The horizons in comoving distance")
    ax.legend(loc="upper left", fontsize=7, bbox_to_anchor=(0.0, 0.92))
    # (b) reachable fraction of the observable universe
    axb = fig.add_subplot(gs[0, 1])
    frac = np.minimum(1, ((u.eta_inf - eta) / eta) ** 3)
    axb.plot(t, 100 * frac, color=ps.C_EVENT, lw=1.8)
    axb.plot(u.age, 100 * reach["fraction_now"], "o", ms=5, mfc=ps.SURFACE, mec=ps.C_EVENT, mew=1.4)
    axb.annotate(f"today: {100 * reach['fraction_now']:.1f}%", (u.age, 100 * reach["fraction_now"]),
                 xytext=(8, 2), textcoords="offset points", fontsize=7.2, color=ps.INK_SECONDARY)
    axb.set_xscale("log")
    axb.set_yscale("log")
    axb.set_xlim(0.5, 100)
    axb.set_ylim(0.01, 150)
    axb.set_xlabel("cosmic time (billion years)")
    axb.set_ylabel("% still reachable")
    axb.set_title("b  Share of the observable universe a signal sent then could reach")
    axb.axvline(reach["t_horizons_equal_gyr"], color=ps.MUTED, lw=0.7)
    axb.text(reach["t_horizons_equal_gyr"] * 0.93, 0.03, f"{reach['t_horizons_equal_gyr']:.1f} Gyr:\n"
             "horizons equal", fontsize=6.8, color=ps.INK_SECONDARY, ha="right")
    # (c) catalogues still reachable
    axc = fig.add_subplot(gs[1, 1])
    order = [("sdss_galaxies", "SDSS galaxies"), ("gw_events", "GW mergers"), ("supernovae", "supernovae"),
             ("quasars", "SDSS quasars"), ("jades", "JWST galaxies"), ("highz_quasars", "quasars z > 5.6")]
    vals = [100 * pops[k]["frac_reachable"] for k, _ in order]
    yy = np.arange(len(order))[::-1]
    axc.barh(yy, vals, height=0.55, color=ps.BLUE)
    for y, v in zip(yy, vals):
        axc.text(v + 1.5, y, f"{v:.1f}%" if 0 < v < 99.95 else f"{v:.0f}%", va="center", fontsize=7.2,
                 color=ps.INK)
    axc.set_yticks(yy)
    axc.set_yticklabels([lab for _, lab in order], fontsize=7.4, color=ps.INK)
    axc.set_xlim(0, 118)
    axc.set_xlabel("% of the catalogue a light signal sent today would still reach")
    axc.set_title("c  Which of the objects we see we could still reach")
    axc.grid(axis="y", visible=False)
    ps.caption(fig, "Planck 2018. Comoving distances are proper distances today; a galaxy sits at a fixed "
                    "comoving distance while the universe expands. A galaxy beyond the event horizon "
                    "(orange) is cut off: light it emits from then on never reaches us, so we only ever see "
                    "its history up to that moment, frozen and ever more redshifted. Panel c counts objects "
                    "with z < 1.852 (Planck), the galaxies still inside our event horizon today.", y=-0.03)
    return ps.save(fig, path)


# ---------------------------------------------------------------------------
# Figure 6: is redshift a Doppler shift? Davis & Lineweaver's supernova test
# ---------------------------------------------------------------------------


def fig_sr(u: Universe, results: dict, path):
    from analysis import hubble_free_dl

    sr = results["sr_test"]
    rep = sr["dl2004_reproduction"]
    pts = rep["points"]
    zp, mp, ep = (np.array(pts[k]) for k in ("z", "m", "err"))
    u_dl = Universe(C.DL2004)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 5.0), gridspec_kw={"wspace": 0.24})
    z = np.linspace(0.01, 1.0, 300)
    ref = 5 * np.log10(hubble_free_dl("frw", z, u_dl))
    refp = 5 * np.log10(hubble_free_dl("frw", zp, u_dl))
    ax1.errorbar(zp, mp - refp - 5 * np.log10(C.C_KM_S) - C.DL2004_TEXT["script_M"][0], yerr=ep, fmt="o",
                 ms=4, mfc=ps.INK_SECONDARY, mec=ps.SURFACE, mew=0.8, ecolor=ps.lighten(ps.INK_SECONDARY, 0.4),
                 elinewidth=0.9, zorder=3, label="Perlmutter et al. 1999, as plotted by D&L")
    curves = [("lcdm", "frw", u_dl, dict(color=ps.INK, lw=1.8), "ΛCDM (0.3, 0.7)"),
              ("SR_as_plotted", "SR_as_plotted", None, dict(color=ps.ORANGE, lw=1.8),
               "SR, as drawn in their Fig. 5 (D_L = D)"),
              ("SR", "SR", None, dict(color=ps.ORANGE, lw=1.3, ls=(0, (5, 3))),
               "SR, as their text describes (D_L = D(1+z))"),
              ("linear_as_plotted", "linear_as_plotted", None, dict(color=ps.BLUE, lw=1.6),
               "v = cz, as drawn (D_L = cz/H)"),
              ("linear", "linear", None, dict(color=ps.BLUE, lw=1.2, ls=(0, (5, 3))),
               "v = cz, as described (= de Sitter)")]
    for key, model, uu, kw, lab in curves:
        ax1.plot(z, 5 * np.log10(hubble_free_dl(model, z, uu)) - ref, label=lab, **kw)
    ax1.axhline(0, color=ps.MUTED, lw=0.6)
    ax1.set_xlim(0, 1.0)
    ax1.set_ylim(-2.1, 1.88)
    ax1.set_xlabel("redshift z")
    ax1.set_ylabel("magnitude relative to ΛCDM (0.3, 0.7)\n(fainter ↑, M = −3.45 fixed as in D&L)")
    ax1.set_title("a  Davis & Lineweaver's test, rebuilt from their own figure")
    ax1.legend(loc="lower left", fontsize=6.9)
    rr = rep["rows"]
    ax1.text(0.98, 0.97, "Δχ² against ΛCDM (their errors scaled ×%.2f to give their χ² = 38.9):\n"
             "SR as drawn %.0f (%.1fσ; they quote 530, 23σ)\n"
             "v = cz as drawn %.0f (%.1fσ; they quote 154, 12σ)\n"
             "SR as described %.0f (%.1fσ); with M free %.1f (%.1fσ)"
             % (rep["error_scale_for_their_chi2"],
                rr["SR_as_plotted"]["dchi2_fixedM_scaled_err"], rr["SR_as_plotted"]["sigma_fixedM_scaled_err"],
                rr["linear_as_plotted"]["dchi2_fixedM_scaled_err"],
                rr["linear_as_plotted"]["sigma_fixedM_scaled_err"],
                rr["SR"]["dchi2_fixedM_scaled_err"], rr["SR"]["sigma_fixedM_scaled_err"],
                rr["SR"]["dchi2_freeM_scaled_err"], rr["SR"]["sigma_freeM_scaled_err"]),
             transform=ax1.transAxes, ha="right", va="top", fontsize=6.8, color=ps.INK_SECONDARY,
             bbox={"fc": ps.SURFACE, "ec": "none", "pad": 1.5})
    # (b) Pantheon+
    pan = sr["pantheon"]
    bn = pan["binned_vs_milne"]
    zb, dm, sg = (np.array(bn[k]) for k in ("z_eff", "dmag", "sigma"))
    ax2.errorbar(zb, dm, yerr=sg, fmt="o", ms=5, mfc=ps.INK, mec=ps.SURFACE, mew=1, ecolor=ps.INK,
                 elinewidth=1, zorder=6, label="Pantheon+ (1,590 light curves, binned)")
    z2 = np.geomspace(0.01, 2.4, 400)
    mods = pan["models"]
    mil = 5 * np.log10(hubble_free_dl("milne", z2)) + mods["milne"]["offset"]
    u_best = Universe(C.PLANCK18.but(Om=pan["Om_best_flat"], T_cmb=0.0, N_eff=0.0, m_nu=()))
    for key, model, uu, kw, lab in [
            ("flat_lcdm_best", "frw", u_best, dict(color=ps.INK, lw=1.8), "flat ΛCDM, best fit"),
            ("milne", "milne", None, dict(color=ps.AQUA, lw=1.6), "empty universe (Milne) = SR done consistently"),
            ("SR", "SR", None, dict(color=ps.ORANGE, lw=1.6, ls=(0, (5, 3))), "SR, D&L's text recipe"),
            ("eds", "frw", Universe(C.EDS), dict(color=ps.MUTED, lw=1.2, ls=(0, (1, 2))), "matter only"),
            ("de_sitter", "frw", Universe(C.DE_SITTER), dict(color=ps.BLUE, lw=1.2, ls=(0, (5, 3))),
             "v = cz (= Λ only)")]:
        y = 5 * np.log10(hubble_free_dl(model, z2, uu)) + mods[key]["offset"] - mil
        lab2 = lab + ("" if key == "flat_lcdm_best" else f": Δχ² = {mods[key]['dchi2']:.0f} "
                                                        f"({mods[key]['sigma']:.1f}σ)")
        ax2.plot(z2, y, label=lab2, **kw)
    ax2.set_xscale("log")
    ax2.set_xlim(0.01, 2.4)
    ax2.set_ylim(-0.75, 0.55)
    ax2.set_xlabel("redshift z")
    ax2.set_ylabel("magnitude relative to the empty universe\n(each model with its best offset)")
    ax2.set_title("b  The same question asked of 1,590 supernovae")
    ax2.legend(loc="lower left", fontsize=6.8)
    ps.caption(fig, "Left: the 42 supernovae of D&L's Fig. 5, read back from their EPS file, with the "
                    "curves they drew. Their 'special-relativistic' and 'v = cz' curves match D_L = D and "
                    "D_L = cz/H to 0.005 mag, not the D_L = D(1+z) and cz(1+z)/H their text gives; their "
                    "23σ and 12σ follow from the drawn curves. Right: Pantheon+ with its full covariance and "
                    "the absolute magnitude free. Treating redshift as a Doppler shift through static space "
                    "is cleanly excluded either way; the self-consistent special-relativistic picture "
                    "(an empty, coasting universe) fails at 7σ, because the expansion accelerates.", y=-0.02)
    return ps.save(fig, path)


# ---------------------------------------------------------------------------
# Figure 7: our recomputation against Davis & Lineweaver's own curves
# ---------------------------------------------------------------------------


def fig_dl_check(results: dict, path):
    import dl2004
    from analysis import _model_x

    u = Universe(C.DL2004)
    src = C.LIT_DIR / "dl2004_src"
    fig, axes = plt.subplots(3, 2, figsize=(12.0, 11.0), gridspec_kw={"width_ratios": [1.35, 1],
                                                                       "hspace": 0.42, "wspace": 0.22})
    cols = {"light_cone": ps.C_CONE, "hubble_sphere": ps.C_HUBBLE, "event_horizon": ps.C_EVENT,
            "particle_horizon": ps.C_PARTICLE}
    names = {"light_cone": "light cone", "hubble_sphere": "Hubble sphere", "event_horizon": "event horizon",
             "particle_horizon": "particle horizon"}
    for row, (fname, panel, title) in enumerate([("fig1a", "a", "proper distance"),
                                                 ("fig1b", "b", "comoving distance"),
                                                 ("fig1c", "c", "conformal time")]):
        cur = dl2004.figure1_curves(src / f"{fname}.eps")
        ax, axr = axes[row]
        ymax = cur["frame"][3]
        for w in cur["worldlines"]:
            ax.plot(w[:, 0], w[:, 1], color=ps.lighten(ps.MUTED, 0.3), lw=2.6, solid_capstyle="butt")
        for key, c in cols.items():
            for w in cur[key]:
                ax.plot(w[:, 0], w[:, 1], color=ps.lighten(c, 0.6), lw=4.0, solid_capstyle="butt", zorder=2)
            # ours on top, thin
            y = np.linspace(0.02, ymax * 0.999, 600)
            x = _model_x(u, panel, key, y)
            if key == "light_cone":
                ok = y < (u.eta0 if panel == "c" else u.age)
                y, x = y[ok], x[ok]
            ax.plot(x, y, color=c, lw=1.0, zorder=4)
            ax.plot(-x, y, color=c, lw=1.0, zorder=4)
            # residuals
            for w in cur[key]:
                xx, yy = np.abs(w[:, 0]), w[:, 1]
                keep = (yy > 0.2) & (xx < cur["frame"][2] - 0.5) & (yy < ymax - 0.05)
                if key == "light_cone":
                    keep &= yy < (u.eta0 if panel == "c" else u.age) - 0.05
                r = xx[keep] - _model_x(u, panel, key, yy[keep])
                axr.plot(yy[keep], r, ".", ms=1.6, color=c, alpha=0.6)
        # worldlines of ours (proper panel: a(t) chi; others vertical)
        for z in (1.0, 3.0, 10.0, 1000.0):
            chi = float(u.chi_of_z(z))
            yy = np.linspace(0.02, ymax, 300)
            xw = u.a_of_t(yy) * chi if panel == "a" else np.full_like(yy, chi)
            ax.plot(xw, yy, color=ps.INK_SECONDARY, lw=0.6, zorder=3)
            ax.plot(-xw, yy, color=ps.INK_SECONDARY, lw=0.6, zorder=3)
        ax.set_xlim(cur["frame"][0], cur["frame"][2])
        ax.set_ylim(0, ymax)
        ax.grid(False)
        ax.set_title(f"{panel}  Fig. 1{panel}, {title}: pale = theirs, thin = ours")
        ax.set_xlabel("distance (billion light-years)")
        ax.set_ylabel("conformal time (Gyr)" if panel == "c" else "cosmic time (Gyr)")
        axr.axhline(0, color=ps.MUTED, lw=0.6)
        rec = results["dl2004"]["panels"][fname]["curves"]
        txt = "\n".join(f"{names[k]}: median {rec[k]['dl2004']['median_gly']:+.3f}, "
                        f"max {rec[k]['dl2004']['max_abs_gly']:.3f}" for k in cols)
        rk = "light_cone" if panel == "c" else "particle_horizon"
        rad = rec[rk]["dl2004_rad"]["median_gly"]
        axr.text(0.02, 0.03, txt + f"\n(with radiation included, the {names[rk]} would be off by {rad:+.2f})",
                 transform=axr.transAxes, fontsize=6.6, color=ps.INK_SECONDARY, va="bottom")
        axr.set_ylim(-0.32, 0.12)
        axr.set_xlabel(axes[row][0].get_ylabel())
        axr.set_ylabel("theirs − ours (billion ly)")
        axr.set_title(f"{panel}  Theirs minus ours")
    from matplotlib.lines import Line2D

    fig.legend(handles=[Line2D([], [], color=cols[k], lw=1.6, label=names[k]) for k in cols]
               + [Line2D([], [], color=ps.INK_SECONDARY, lw=0.6, label="worldlines z = 1, 3, 10, 1000")],
               loc="upper left", bbox_to_anchor=(0.12, 0.935), ncol=5, fontsize=7.6)
    ps.caption(fig, "Davis & Lineweaver (2004), Fig. 1, read back from the vector EPS in their arXiv source "
                    "(astro-ph/0310808 v2) and redrawn under our recomputation for their universe: Ωm = 0.3, "
                    "ΩΛ = 0.7, H0 = 70 km/s/Mpc, no radiation. Light cone, Hubble sphere and event horizon "
                    "agree to a few hundredths of a billion light-years (the page resolution is 0.011); their "
                    "conformal times run 0.05 (today) to 0.16 (t = ∞) billion years short of the exact "
                    "integrals, which shifts their particle horizon by ~0.1.", y=0.04)
    return ps.save(fig, path)


# ---------------------------------------------------------------------------
# Figure 8: what "distance" means
# ---------------------------------------------------------------------------


def fig_distances(u: Universe, results: dict, path):
    desi = results["desi"]
    rd = desi["r_drag_mpc"] * C.GLY_PER_MPC
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    z = np.geomspace(0.01, 1500, 1200)
    curves = [(u.chi_of_z(z), "comoving distance: where the source is now", dict(color=ps.INK, lw=2.0)),
              (u.dl_of_z(z), "luminosity distance (how faint it looks)", dict(color=ps.ORANGE, lw=1.5)),
              (u.light_travel_distance(z), "light-travel distance (c × time in flight)",
               dict(color=ps.AQUA, lw=1.5)),
              (u.emission_distance(z), "proper distance when the light left (= D_A)", dict(color=ps.BLUE, lw=1.5))]
    for y, lab, kw in curves:
        ax.plot(z, y, label=lab, **kw)
    ax.plot(z, z * u.hubble_distance, color=ps.MUTED, lw=1.0, ls=(0, (1, 2)), label="naive cz / H0")
    rows = desi["rows"]
    ax.errorbar([r["z"] for r in rows], [r["DA_over_rd"] * (1 + r["z"]) * rd for r in rows],
                yerr=[r["DA_over_rd_sigma"] * (1 + r["z"]) * rd for r in rows], fmt="o", ms=5,
                mfc=ps.INK, mec=ps.SURFACE, mew=1, ecolor=ps.INK, label="DESI D_M (× Planck r_d)", zorder=6)
    ax.errorbar([r["z"] for r in rows], [r["DA_over_rd"] * rd for r in rows],
                yerr=[r["DA_over_rd_sigma"] * rd for r in rows], fmt="o", ms=5, mfc=ps.BLUE,
                mec=ps.SURFACE, mew=1, ecolor=ps.BLUE, label="DESI D_A (× Planck r_d)", zorder=6)
    sb = results["sn_velocities"]["bins"]
    ax.errorbar([b["z_eff"] for b in sb], [b["v_now_c"] * (1 + b["z_eff"]) * u.hubble_distance for b in sb],
                yerr=[b["v_now_sigma"] * (1 + b["z_eff"]) * u.hubble_distance for b in sb], fmt="s", ms=4.5,
                mfc=ps.ORANGE, mec=ps.SURFACE, mew=1, ecolor=ps.ORANGE, zorder=6,
                label="Pantheon+ D_L (× Planck c/H0)")
    cmb = results["headline"]["cmb"]
    for key, lab in (("chi_gly", f"{cmb['chi_gly']:.1f} billion"),
                     ("light_travel_gly", f"{cmb['light_travel_gly']:.1f} billion"),
                     ("d_emit_gly", f"{cmb['d_emit_gly'] * 1000:.0f} million")):
        ax.plot(cmb["z"], cmb[key], "o", ms=4.5, mfc=ps.SURFACE, mec=ps.INK, mew=1.1, zorder=7)
        ax.annotate(f"CMB: {lab} ly", (cmb["z"], cmb[key]), xytext=(-8, 5), textcoords="offset points",
                    ha="right", fontsize=7, color=ps.INK_SECONDARY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.01, 1500)
    ax.set_ylim(0.03, 200)
    ax.set_xlabel("redshift z")
    ax.set_ylabel("distance (billion light-years)")
    ax.set_title("Five meanings of 'how far away', for the same source")
    ax.legend(loc="upper left", fontsize=7.2)
    ps.caption(fig, "Planck 2018. The CMB's light has travelled for 13.8 billion years; the gas that emitted "
                    "it is now 45.3 billion light-years away, and was only 42 million light-years away when "
                    "the light set out. DESI measures D_M and D_A in units of the sound horizon r_d (here "
                    "multiplied by Planck's r_d = 147.21 Mpc); supernovae measure H0 D_L (here multiplied by "
                    "Planck's c/H0 = 14.45 billion ly).", y=-0.02)
    return ps.save(fig, path)


# ---------------------------------------------------------------------------
# Figure 9: how much the answers depend on the cosmology
# ---------------------------------------------------------------------------


def fig_sensitivity(results: dict, path):
    rows = results["sensitivity"]
    unc = results["uncertainties"]
    err = {"planck18": unc["planck18"], "d3_tt": unc.get("d3_tt")}
    labels = {"planck18": "Planck 2018 (+BAO)", "dl2004": "Davis & Lineweaver 2004 (0.3, 0.7, 70)",
              "dl2004_rad": "same, with radiation", "d3_tt": "D3: our fit to Planck TT",
              "d5_lcdm": "D5: DESI DR2 + CMB", "d1_shoes": "D1: SH0ES H0 = 73.6 + Pantheon+",
              "d5_w0wa_frozen": "D5: DESI + CMB + DESY5 w0wa, w frozen after today",
              "d5_w0wa_cpl": "same, w(a) extrapolated: dark energy fades"}
    order = ["planck18", "d3_tt", "d5_lcdm", "d1_shoes", "dl2004", "dl2004_rad", "d5_w0wa_frozen", "d5_w0wa_cpl"]
    by = {r["name"]: r for r in rows}
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.9), sharey=True, gridspec_kw={"wspace": 0.08})
    qs = [("particle_horizon_gly", "particle horizon today\n(billion light-years)", (40, 47.5)),
          ("event_horizon_gly", "event horizon today\n(billion light-years)", (14.5, 28)),
          ("z_hubble_sphere", "redshift receding at c today", (1.42, 1.56))]
    yy = np.arange(len(order))[::-1]
    for ax, (q, lab, lim) in zip(axes, qs):
        for y, key in zip(yy, order):
            v = by[key][q]
            if isinstance(v, str) or not np.isfinite(v):
                ax.text(lim[1] - 0.01 * (lim[1] - lim[0]), y, "no event horizon →", ha="right", va="center",
                        fontsize=7, color=ps.INK_SECONDARY)
                continue
            e = err.get(key)
            col = ps.INK if key == "planck18" else ps.INK_SECONDARY
            if e is not None and q in e:
                ax.errorbar(v, y, xerr=e[q]["sigma"], fmt="o", ms=5.5, color=col, mfc=col, mec=ps.SURFACE,
                            elinewidth=1.2, capsize=0, zorder=5)
            else:
                ax.plot(v, y, "o", ms=5.5, mfc=col, mec=ps.SURFACE, mew=1, zorder=5)
            ax.annotate(f"{v:.2f}" if q != "z_hubble_sphere" else f"{v:.3f}", (v, y), xytext=(0, 6),
                        textcoords="offset points", ha="center", fontsize=6.6, color=ps.INK_SECONDARY)
        ax.set_xlim(*lim)
        ax.set_xlabel(lab)
        ax.grid(axis="y", visible=False)
    axes[0].set_yticks(yy)
    axes[0].set_yticklabels([labels[k] for k in order], fontsize=7.4, color=ps.INK)
    axes[0].set_ylim(-0.7, len(order) - 0.4)
    axes[0].set_title("How the horizons depend on the cosmology (error bars: Planck 2018 and D3's fit, 1σ)",
                      loc="left")
    return ps.save(fig, path)


def make_all(results: dict | None = None) -> list:
    ps.use_style()
    if results is None:
        results = json.loads(C.RESULTS_PATH.read_text())
    FIG.mkdir(parents=True, exist_ok=True)
    u = Universe(C.PLANCK18)
    arrays = np.load(C.DERIVED_DIR / "cone_objects.npz")
    out = []
    out += fig_spacetime(u, results, FIG / "01_spacetime_diagrams",
                         title="Our light cone and horizons, Planck 2018")
    out += fig_objects(u, results, arrays, FIG / "02_objects_on_the_light_cone")
    out += fig_speeds(u, results, FIG / "03_recession_speeds")
    out += fig_angular(u, results, FIG / "04_angular_size_turnover")
    out += fig_reach(u, results, FIG / "05_horizons_and_reach")
    out += fig_sr(u, results, FIG / "06_doppler_or_expansion")
    out += fig_dl_check(results, FIG / "07_davis_lineweaver_check")
    out += fig_distances(u, results, FIG / "08_five_distances")
    out += fig_sensitivity(results, FIG / "09_cosmology_dependence")
    return out


if __name__ == "__main__":
    for p in make_all():
        print(p)
