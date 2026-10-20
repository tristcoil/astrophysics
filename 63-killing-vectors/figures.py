"""The figures of project N6, drawn from the cached analysis results."""

from __future__ import annotations

import pickle

import matplotlib.pyplot as plt
import numpy as np

import config as C
import literature as lit
import plotstyle as ps

ps.use_style()
FIG = C.FIGURE_DIR


def load(name):
    return pickle.loads((C.CACHE_DIR / f"{name}.pkl").read_bytes())


# ---------------------------------------------------------------------------
# 1. How many symmetries?
# ---------------------------------------------------------------------------

def killing_counts(sym):
    vb = sym["vector_bounds"]
    explicit = {"Minkowski": "Minkowski", "Schwarzschild": "Schwarzschild", "Kerr": "Kerr",
                "FLRW, matter era": "FLRW, matter era", "de Sitter (FLRW form)": "de Sitter (conformal time)"}
    names = list(vb)
    fig, axs = plt.subplots(3, 5, figsize=(13, 7.4), sharex=True, sharey=True)
    for ax, name in zip(axs.flat, names):
        d = vb[name]
        seq = d["sequence"]
        N = np.arange(1, len(seq) + 1)
        ax.step(N, seq, where="mid", color=ps.BLUE, lw=1.6, zorder=3)
        ax.plot(N, seq, "o", color=ps.BLUE, ms=5, mec=ps.SURFACE, mew=1.0, zorder=4)
        known = d["known"] if d["known"] is not None else C.KILLING_COUNTS.get(name)
        if known is not None:
            ax.axhline(known, color=ps.MUTED, lw=0.9, ls=(0, (3, 2)), zorder=2)
        ex = explicit.get(name)
        if ex and ex in sym["explicit"]:
            nf = sym["explicit"][ex]["n_found"]
            ax.plot([N[-1] + 0.45], [nf], marker="D", color=ps.AQUA, ms=6, mec=ps.SURFACE, mew=1.0, zorder=5,
                    clip_on=False)
        ax.set_title(f"{name}\n{d['final']} Killing vectors", fontsize=8.5, pad=4)
        ax.set_ylim(-0.5, 11)
        ax.set_xlim(0.5, 6.9)
        ax.set_yticks([0, 2, 4, 6, 8, 10])
    for ax in axs.flat[len(names):]:
        ax.axis("off")
    ax = axs.flat[len(names)]
    ax.text(0.0, 0.95, "How to read the panels", fontsize=9, fontweight="bold", va="top", transform=ax.transAxes)
    ax.text(0.0, 0.80,
            "Blue: the most symmetries the Taylor\nexpansion of Killing's equation allows,\n"
            "at expansion order N (exact arithmetic).\nIt starts at 10 and drops as higher\n"
            "orders feel the curvature.\n\nGrey dashes: the textbook answer.\n"
            "Aqua diamond: independent Killing\nvectors found explicitly and checked\nby substitution.",
            fontsize=7.8, va="top", color=ps.INK_SECONDARY, transform=ax.transAxes)
    for ax in axs[-1]:
        ax.set_xlabel("expansion order N")
    for ax in axs[:, 0]:
        ax.set_ylabel("Killing vectors (upper bound)")
    fig.suptitle("How many symmetries does each spacetime have?", x=0.01, ha="left", y=1.0)
    fig.tight_layout()
    ps.caption(fig, "Upper bound = dimension of the order-N solutions of Killing's equation projected onto a vector "
                    "field's value and first derivatives at one generic point (Taylor jets, exact rationals). "
                    "Where the bound meets an explicit set of verified solutions, the count is exact.", y=0.0)
    return ps.save(fig, FIG / "01_symmetry_counts")


# ---------------------------------------------------------------------------
# 2. What Killing vectors look like
# ---------------------------------------------------------------------------

def _flow(ax, f, xlim, ylim, density=1.1, color=ps.CONTEXT, mask=None):
    X, Y = np.meshgrid(np.linspace(*xlim, 300), np.linspace(*ylim, 300))
    U, V = f(X, Y)
    if mask is not None:
        bad = mask(X, Y)
        U, V = np.ma.array(U, mask=bad), np.ma.array(V, mask=bad)
    ax.streamplot(X, Y, U, V, density=density, color=color, linewidth=0.7, arrowsize=0.7)


def killing_fields(sym):
    fig, axs = plt.subplots(1, 4, figsize=(15, 4.3))
    # (a) a boost in flat spacetime
    ax = axs[0]
    lim = (-2, 2)
    X, T = np.meshgrid(np.linspace(*lim, 400), np.linspace(*lim, 400))
    norm = -X**2 + T**2  # g(K, K) for K = x d_t + t d_x
    ax.contourf(X, T, (norm < 0).astype(float), levels=[-0.5, 0.5, 1.5], colors=[ps.SURFACE, ps.BLUE_100], zorder=0)
    _flow(ax, lambda x_, t_: (t_, x_), lim, lim)
    ax.plot([-2, 2], [-2, 2], color=ps.INK, lw=1.0)
    ax.plot([-2, 2], [2, -2], color=ps.INK, lw=1.0)
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_aspect("equal")
    ax.set_xlabel("x (space)")
    ax.set_ylabel("t (time)")
    ax.set_title("A boost in flat spacetime")
    ax.text(1.2, 0.05, "timelike\n(Rindler wedge)", ha="center", fontsize=7.5, color=ps.BLUE_550)
    ax.text(0.0, 1.45, "spacelike", ha="center", fontsize=7.5, color=ps.INK_SECONDARY)
    ps.annotate(ax, "light cone: boost is null", (-1.3, 1.3), (-1.9, 1.75))
    # (b) Schwarzschild's time translation in Kruskal coordinates
    ax = axs[1]
    lim = (-2, 2)
    Xk, Tk = np.meshgrid(np.linspace(*lim, 400), np.linspace(*lim, 400))
    region = (Xk**2 - Tk**2 > 0)
    sing = Tk**2 - Xk**2 > 1
    ax.contourf(Xk, Tk, region.astype(float), levels=[-0.5, 0.5, 1.5], colors=[ps.SURFACE, ps.BLUE_100], zorder=0)
    ax.contourf(Xk, Tk, sing.astype(float), levels=[0.5, 1.5], colors=[ps.WASH], zorder=1)
    _flow(ax, lambda x_, t_: (t_, x_), lim, lim, mask=lambda x_, t_: t_**2 - x_**2 > 1)  # d_t = (T d_X + X d_T)/4M
    ax.plot([-2, 2], [-2, 2], color=ps.INK, lw=1.2)
    ax.plot([-2, 2], [2, -2], color=ps.INK, lw=1.2)
    xs = np.linspace(-2, 2, 400)
    for sgn in (1, -1):
        ax.plot(xs, sgn * np.sqrt(1 + xs**2), color=ps.INK_SECONDARY, lw=1.6, ls=(0, (1, 1.2)))
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_aspect("equal")
    ax.set_xlabel("Kruskal X")
    ax.set_ylabel("Kruskal T")
    ax.set_title("Schwarzschild's 'time' symmetry")
    ax.text(1.3, -0.05, "outside:\ntimelike", ha="center", fontsize=7.5, color=ps.BLUE_550)
    ax.text(0, 0.6, "inside:\nspacelike", ha="center", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.text(0, 1.62, "singularity r = 0", ha="center", fontsize=7.5, color=ps.INK_SECONDARY)
    ps.annotate(ax, "horizon: the Killing\nvector is null", (-0.9, -0.9), (-1.95, -1.75))
    # (c) Kerr: the ergosphere, where d_t is spacelike
    ax = axs[2]
    a_ = 0.95
    th = np.linspace(0, np.pi, 400)
    rp = 1 + np.sqrt(1 - a_**2)
    rE = 1 + np.sqrt(1 - a_**2 * np.cos(th) ** 2)
    for rad, kw in ((rE, {"color": ps.ORANGE, "lw": 1.6}), (np.full_like(th, rp), {"color": ps.INK, "lw": 1.6})):
        xx = np.sqrt(rad**2 + a_**2) * np.sin(th)
        zz = rad * np.cos(th)
        ax.plot(xx, zz, **kw)
        ax.plot(-xx, zz, **kw)
    xo, zo = np.sqrt(rE**2 + a_**2) * np.sin(th), rE * np.cos(th)
    xi, zi = np.sqrt(rp**2 + a_**2) * np.sin(th), rp * np.cos(th)
    for sgn in (1, -1):
        ax.fill(np.concatenate([sgn * xo, sgn * xi[::-1]]), np.concatenate([zo, zi[::-1]]), color=ps.ORANGE,
                alpha=0.18, lw=0)
        ax.fill(sgn * xi, zi, color=ps.INK, alpha=0.85, lw=0)
    import orbits as ob
    risco = float(ob.isco_bpt(a_))
    for sgn in (1, -1):
        ax.plot([sgn * np.sqrt(risco**2 + a_**2)], [0], "o", color=ps.BLUE, ms=6, mec=ps.SURFACE)
    ax.set_xlim(-3.2, 3.2)
    ax.set_ylim(-2.6, 2.6)
    ax.set_aspect("equal")
    ax.set_xlabel("distance from axis (GM/c²)")
    ax.set_ylabel("along spin axis (GM/c²)")
    ax.set_title("Kerr (a = 0.95): the ergosphere")
    ps.annotate(ax, "static limit: d/dt\nturns spacelike", (-1.9, 0.4), (-3.1, 1.9))
    ps.annotate(ax, "horizon", (0.2, -1.25), (0.9, -2.2))
    ps.annotate(ax, "prograde ISCO", (np.sqrt(risco**2 + a_**2), 0), (1.9, -1.0))
    # (d) the three rotations of a sphere
    ax = axs[3]
    t_ = np.linspace(0, 2 * np.pi, 300)
    tilt = np.radians(20)
    colors = [ps.BLUE, ps.ORANGE, ps.AQUA]
    labels = ["J_x", "J_y", "J_z"]
    # orthographic view direction
    view = np.array([np.cos(np.radians(25)) * np.cos(np.radians(-35)), np.cos(np.radians(25)) * np.sin(np.radians(-35)),
                     np.sin(np.radians(25))])
    e1 = np.cross([0, 0, 1], view)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(view, e1)
    ax.add_patch(plt.Circle((0, 0), 1, fill=False, color=ps.AXIS, lw=1.0))
    for k_, (col, lab) in enumerate(zip(colors, labels)):
        axis = np.eye(3)[k_]
        u_ = np.eye(3)[(k_ + 1) % 3]
        w_ = np.eye(3)[(k_ + 2) % 3]
        for h in (-0.75, -0.35, 0.0, 0.35, 0.75):
            rr = np.sqrt(1 - h**2)
            P = h * axis[:, None] + rr * (np.cos(t_) * u_[:, None] + np.sin(t_) * w_[:, None])
            front = P.T @ view >= 0
            xx, yy = P.T @ e1, P.T @ e2
            xs_ = np.where(front, xx, np.nan)
            ax.plot(xs_, np.where(front, yy, np.nan), color=col, lw=1.1)
            ax.plot(np.where(~front, xx, np.nan), np.where(~front, yy, np.nan), color=col, lw=0.5, alpha=0.35)
        ax.plot([], [], color=col, lw=1.4, label=f"rotation {lab}")
    ax.set_aspect("equal")
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.axis("off")
    ax.legend(loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.12), fontsize=7.5, handlelength=1.2)
    ax.set_title("Three rotations: angular momentum")
    fig.tight_layout()
    ps.caption(fig, "Flow lines show the direction in which each symmetry moves points. Shading: where the Killing "
                    "vector is timelike. A horizon is where a time-translation Killing vector becomes null "
                    "(a 'Killing horizon'); in Kruskal coordinates Schwarzschild's d/dt is literally a boost, with "
                    "rate kappa = c^4/4GM, the surface gravity that sets Hawking's temperature.", y=0.0)
    return ps.save(fig, FIG / "02_killing_fields")


# ---------------------------------------------------------------------------
# 3. Effective potentials
# ---------------------------------------------------------------------------

def potentials(orb):
    import orbits as ob
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
    ax = axs[0]
    r = np.linspace(2.05, 30, 2000)
    Ls = [3.2, np.sqrt(12), 3.8, 4.3]
    cols = [ps.BLUE_300, ps.BLUE, ps.BLUE_550, ps.BLUE_700]
    for L, c_ in zip(Ls, cols):
        V = np.sqrt(ob.V2_schwarzschild(r, L))
        lab = "L = 2√3 (ISCO)" if abs(L - np.sqrt(12)) < 1e-9 else f"L = {L:.1f}"
        ax.plot(r, V, color=c_, lw=1.6 if abs(L - np.sqrt(12)) > 1e-9 else 2.2, label=lab)
    ax.axhline(1, color=ps.MUTED, lw=0.8, ls=(0, (3, 2)))
    E6, L6 = ob.circular_schwarzschild(6.0)
    ax.plot([6], [E6], "o", color=ps.INK, ms=6, mec=ps.SURFACE, zorder=5)
    ps.annotate(ax, "ISCO, r = 6GM/c²\nE = √(8/9): 5.7% of\nthe rest energy released", (6, E6), (9.5, 0.905))
    ax.set_ylim(0.88, 1.06)
    ax.set_xlabel("r (GM/c²)")
    ax.set_ylabel("effective potential V(r)  (energy per unit rest mass)")
    ax.set_title("Schwarzschild: massive particles")
    ax.legend(loc="lower right", fontsize=7.5)
    ax = axs[1]
    r = np.linspace(2.02, 20, 2000)
    b_eff = r / np.sqrt(1 - 2 / r)
    ax.plot(r, b_eff, color=ps.BLUE, lw=1.8)
    ax.axhline(np.sqrt(27), color=ps.MUTED, lw=0.8, ls=(0, (3, 2)))
    ax.plot([3], [np.sqrt(27)], "o", color=ps.INK, ms=6, mec=ps.SURFACE, zorder=5)
    ps.annotate(ax, "photon sphere r = 3GM/c²\ncritical b = √27 GM/c²", (3, np.sqrt(27)), (6.5, 3.6))
    ax.fill_between(r, 0, np.minimum(b_eff, np.sqrt(27)), color=ps.WASH, zorder=0)
    ax.text(12, 2.4, "rays with b < √27 GM/c²\nare captured: the shadow", ha="center", fontsize=7.5,
            color=ps.INK_SECONDARY)
    ax.set_ylim(0, 14)
    ax.set_xlabel("r (GM/c²)")
    ax.set_ylabel("turning-point impact parameter b(r) (GM/c²)")
    ax.set_title("Schwarzschild: light")
    ax = axs[2]
    a_ = 0.9
    r = np.geomspace(1 + np.sqrt(1 - a_**2) + 0.01, 30, 4000)
    pro = ob.isco_kerr_derived(a_)
    ret = ob.isco_kerr_derived(-a_)
    Vp = ob.V_eff_kerr(r, pro["L"], a_)
    Vr = ob.V_eff_kerr(r, abs(ret["L"]), -a_)
    ax.plot(r, Vp, color=ps.PRO, lw=1.8, label=f"prograde, ISCO at {pro['r_isco']:.2f}")
    ax.plot(r, Vr, color=ps.RETRO, lw=1.8, ls=(0, (5, 2)), label=f"retrograde, ISCO at {ret['r_isco']:.2f}")
    ax.plot([pro["r_isco"]], [pro["E"]], "o", color=ps.PRO, ms=6, mec=ps.SURFACE, zorder=5)
    ax.plot([ret["r_isco"]], [ret["E"]], "s", color=ps.RETRO, ms=6, mec=ps.SURFACE, zorder=5)
    ax.axhline(1, color=ps.MUTED, lw=0.8, ls=(0, (3, 2)))
    ax.axvline(1 + np.sqrt(1 - a_**2), color=ps.INK, lw=1.0)
    ax.text(1 + np.sqrt(1 - a_**2) + 0.2, 0.71, "horizon", fontsize=7.5, color=ps.INK_SECONDARY, rotation=90)
    ax.set_ylim(0.7, 1.05)
    ax.set_xscale("log")
    ax.set_xticks([1.5, 2, 3, 5, 10, 20, 30])
    ax.set_xticklabels(["1.5", "2", "3", "5", "10", "20", "30"])
    ax.minorticks_off()
    ax.set_xlabel("r (GM/c²), log scale")
    ax.set_ylabel("effective potential (equatorial)")
    ax.set_title("Kerr, a = 0.9: ISCO-energy orbits")
    ax.legend(loc="lower right", fontsize=7.5)
    fig.tight_layout()
    ps.caption(fig, "Each curve uses only the conserved energy and angular momentum from the time-translation and "
                    "rotation Killing vectors. Prograde orbits (with the spin) reach closer and release more energy.",
               y=0.0)
    return ps.save(fig, FIG / "03_effective_potentials")


# ---------------------------------------------------------------------------
# 4. Spin: special radii, efficiency, rotational energy
# ---------------------------------------------------------------------------

def spin_figure(orb, th):
    import orbits as ob
    tab = orb["table"]
    a = np.array(tab["a"])
    fig, axs = plt.subplots(1, 3, figsize=(14.5, 4.4))
    ax = axs[0]
    for key_p, key_r, col, lab in (("r_isco_pro", "r_isco_retro", ps.BLUE, "innermost stable orbit"),
                                   ("r_mb_pro", "r_mb_retro", ps.AQUA, "marginally bound orbit"),
                                   ("r_ph_pro", "r_ph_retro", ps.ORANGE, "photon orbit")):
        ax.plot(a, tab[key_p], color=col, lw=1.8, label=lab)
        ax.plot(a, tab[key_r], color=col, lw=1.5, ls=(0, (5, 2)))
    ax.plot(a, tab["r_plus"], color=ps.INK, lw=1.8, label="horizon")
    ax.text(0.03, 0.45, "dashed: retrograde orbits (against the spin)", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_xlabel("spin a = cJ/GM²")
    ax.set_ylabel("radius (GM/c²)")
    ax.set_ylim(0, 9.5)
    ax.set_title("Special orbits against spin")
    ax.legend(loc="upper left", fontsize=7.5)
    # (b) efficiency with real spins
    ax = axs[1]
    af = np.linspace(0, 0.9999, 600)
    ax.plot(af, 100 * ob.efficiency_from_spin(af), color=ps.PRO, lw=1.8, label="prograde disc")
    ax.plot(af, 100 * ob.efficiency_from_spin(-af), color=ps.RETRO, lw=1.6, ls=(0, (5, 2)), label="retrograde disc")
    for s in lit.XRB_SPINS:
        eff = 100 * float(ob.efficiency_from_spin(s["a"]))
        mk, col = ("o", ps.BLUE) if s["method"] == "continuum" else ("s", ps.ORANGE)
        if s["kind"] == "value":
            lo, hi = max(s["a"] - s["lo"], 0), min(s["a"] + s["hi"], 0.9999)
            seg = np.linspace(lo, hi, 60)
            ax.plot(seg, 100 * ob.efficiency_from_spin(seg), color=col, lw=3.0, alpha=0.25, solid_capstyle="round")
        else:
            ax.annotate("", xy=(min(s["a"] + 0.05, 1.0), eff), xytext=(s["a"], eff),
                        arrowprops={"arrowstyle": "->", "color": col, "lw": 0.8})
        ax.plot([s["a"]], [eff], mk, color=col, ms=5.5, mec=ps.SURFACE, mew=0.8, zorder=4)
    gw = th["census"]["gw250114"]
    for name, chi in (("GW250114", gw["chi_final"]["median"]),):
        ax.plot([chi], [100 * ob.efficiency_from_spin(chi)], "D", color=ps.AQUA, ms=6, mec=ps.SURFACE, zorder=5)
        ps.annotate(ax, "GW250114 remnant", (chi, 100 * float(ob.efficiency_from_spin(chi))), (0.25, 17))
    ps.annotate(ax, "Cyg X-1 (a* > 0.9985)", (0.9985, 100 * float(ob.efficiency_from_spin(0.9985))), (0.55, 28))
    ax.plot([], [], "o", color=ps.BLUE, label="X-ray binary, continuum fitting")
    ax.plot([], [], "s", color=ps.ORANGE, label="X-ray binary, iron-line reflection")
    ax.plot([], [], "D", color=ps.AQUA, label="merger remnant (LVK)")
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 45)
    ax.set_xlabel("spin a")
    ax.set_ylabel("accretion efficiency 1 − E_ISCO (%)")
    ax.set_title("How much light a disc can make")
    ax.legend(loc="upper left", fontsize=7)
    # (c) rotational energy
    ax = axs[2]
    m_irr = np.sqrt((1 + np.sqrt(1 - af**2)) / 2)
    ax.plot(af, 100 * (1 - m_irr), color=ps.BLUE, lw=1.8)
    ax.axhline(100 * (1 - 1 / np.sqrt(2)), color=ps.MUTED, lw=0.8, ls=(0, (3, 2)))
    ax.text(0.02, 100 * (1 - 1 / np.sqrt(2)) + 0.6, "29.3%: the most any black hole can hold", fontsize=7.5,
            color=ps.INK_SECONDARY)
    for r_ in th["census"]["rows"]:
        if r_["chi"] is None:
            continue
        mk, col = ("D", ps.AQUA) if r_["group"] == "merger remnant" else ("o", ps.BLUE)
        ax.plot([r_["chi"]], [100 * r_["E_rot_frac"]], mk, color=col, ms=5.5, mec=ps.SURFACE, mew=0.8, zorder=4)
    ps.annotate(ax, "GW250114 remnant: 6.8%\n(4 suns' worth)", (gw["chi_final"]["median"],
                                                             100 * gw["E_rot_final_frac"]["median"]), (0.12, 9.5))
    ax.plot([], [], "o", color=ps.BLUE, label="X-ray binary (published spin, I4 mass)")
    ax.plot([], [], "D", color=ps.AQUA, label="merger remnant")
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 32)
    ax.set_xlabel("spin a")
    ax.set_ylabel("extractable rotational energy 1 − M_irr/M (%)")
    ax.set_title("Spin energy (Penrose process)")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.92), fontsize=7)
    fig.tight_layout()
    ps.caption(fig, "Radii from solving the Killing-reduced radial equation numerically (they match Bardeen, Press & "
                    "Teukolsky 1972 to machine precision). Spins: published X-ray values (literature.py; strongly "
                    "model-dependent, see README) and LVK's GW250114 posterior. M_irr = sqrt(A/16pi) is the mass "
                    "the area theorem does not let anything take away.", y=0.0)
    return ps.save(fig, FIG / "04_spin")


# ---------------------------------------------------------------------------
# 5. S2: what is conserved, season by season
# ---------------------------------------------------------------------------

def s2_conservation(s2r):
    se = s2r["seasons"]
    rows = se["seasons"]
    fig, axs = plt.subplots(2, 2, figsize=(12.5, 7.6))
    t = np.array([r_["t_center"] for r_ in rows])
    peri = np.array([r_["near_pericentre"] for r_ in rows])
    h = np.array([r_["h_best"] for r_ in rows])
    herr = np.array([np.hypot(r_["h_stat_err"], r_["h_frame_sd"]) for r_ in rows])
    mh = se["model_h"]
    ax = axs[0, 0]
    ax.axhspan(mh["median"] - mh["sd"], mh["median"] + mh["sd"], color=ps.BLUE_150, lw=0, zorder=0)
    ax.axhline(mh["median"], color=ps.BLUE, lw=1.4, label="whole-orbit fit (H1 posterior)")
    ax.errorbar(t[~peri], h[~peri], herr[~peri], fmt="o", color=ps.INK, ms=5, mec=ps.SURFACE, capsize=0, lw=1,
                label="one season's positions only")
    ax.errorbar(t[peri], h[peri], herr[peri], fmt="o", color=ps.INK, mfc=ps.SURFACE, ms=5, capsize=0, lw=1,
                label="pericentre season (S2 turns too fast)")
    c = se["h_consistency"]
    ax.set_title("Equal areas in equal times: S2's angular momentum")
    ax.set_ylabel("sky-projected angular momentum\nX dY/dt − Y dX/dt (mas²/yr)")
    ax.text(0.02, 0.04, f"{c['n']} seasons: weighted mean {c['weighted_mean']:.0f} ± {c['weighted_mean_err']:.0f}, "
                        f"fit {c['model_median']:.0f};  χ² = {c['chi2_vs_model']:.1f} for {c['n']}",
            transform=ax.transAxes, fontsize=7.5, color=ps.INK_SECONDARY)
    ax.legend(loc="upper left", fontsize=7.2)
    ax.set_ylim(min(h.min() - 400, 800), 3200)
    ax = axs[1, 0]
    v = np.array([r_["sky_speed_kms"] for r_ in rows])
    ax.plot(t, v, "o-", color=ps.INK_SECONDARY, ms=5, mec=ps.SURFACE, lw=0.8)
    ax.set_yscale("log")
    ax.set_ylabel("speed across the sky (km/s)")
    ax.set_xlabel("year")
    ax.set_title("...while its speed changes 13-fold")
    ax.set_xlim(axs[0, 0].get_xlim())
    ax = axs[0, 1]
    E = np.array([r_.get("E_best", np.nan) for r_ in rows]) / 1e6
    Eerr = np.array([np.hypot(r_.get("E_stat_err", np.nan), r_.get("E_frame_sd", np.nan)) for r_ in rows]) / 1e6
    mE = se["model_E"]
    ok = np.isfinite(E) & ~peri
    ax.axhspan((mE["median"] - mE["sd"]) / 1e6, (mE["median"] + mE["sd"]) / 1e6, color=ps.BLUE_150, lw=0)
    ax.axhline(mE["median"] / 1e6, color=ps.BLUE, lw=1.4, label="whole-orbit fit: −GM/2a")
    ax.errorbar(t[ok], E[ok], Eerr[ok], fmt="o", color=ps.INK, ms=5, mec=ps.SURFACE, lw=1, capsize=0,
                label="one season: v²/2 − GM/r")
    ce = se["E_consistency"]
    ax.text(0.02, 0.04, f"{ce['n']} seasons: χ² = {ce['chi2_vs_model']:.1f} for {ce['n']}",
            transform=ax.transAxes, fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_ylabel("energy per unit mass (10⁶ km²/s²)")
    ax.set_title("S2's orbital energy")
    ax.set_ylim(-2.1, -1.3)
    ax.legend(loc="upper left", fontsize=7.2)
    ax = axs[1, 1]
    g = s2r["geodesic"]
    tt = np.array(g["t_yr"])
    sel = np.abs(tt) < 1.0
    ax.plot(tt[sel], 100 * np.array(g["E_newton_rel"])[sel], color=ps.NEWTON, lw=1.8,
            label="Newtonian energy v²/2 − GM/r (% of binding energy)")
    ax.plot(tt[sel], 100 * np.array(g["L_newton_rel"])[sel], color=ps.AQUA, lw=1.8, ls=(0, (5, 2)),
            label="Newtonian angular momentum r² dφ/dt (%)")
    ax.axhline(0, color=ps.GR, lw=1.6, label="Killing E and L (exactly constant)")
    ax.set_xlabel("years from pericentre")
    ax.set_ylabel("change along the orbit (%)")
    ax.set_title("General relativity vs Newton along S2's orbit")
    ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.5), fontsize=7.2)
    ax.text(0.55, 0.06, f"Newtonian E swings by {100 * g['E_newton_swing_frac']:.2f}% of the\nbinding energy and L by "
                        f"{100 * g['L_newton_swing_frac']:.3f}%, all within\nweeks of pericentre: below one season's\n"
                        "errors (3-10%). What H1 measured there is\nthe photons' Killing energy: the 200 km/s redshift.",
            transform=ax.transAxes, fontsize=7.0, color=ps.INK_SECONDARY)
    fig.tight_layout()
    ps.caption(fig, "Data: H1 (17-galactic-centre), Keck 1995-2018 (Do et al. 2019) and NTT/VLT 1992-2016 "
                    "(Gillessen et al. 2017). Each season fitted on its own with a quadratic in light-travel-corrected "
                    "time, relative to the black hole as located by H1's fit (each team's frame offset and drift "
                    "included; error bars add the spread over 400 posterior draws). Energy also uses that season's "
                    "spectra and the orbital plane.", y=0.0)
    return ps.save(fig, FIG / "05_s2_conservation")


# ---------------------------------------------------------------------------
# 6. S2 in its effective potential; precession from E and L
# ---------------------------------------------------------------------------

def s2_potential(s2r):
    import orbits as ob
    sm = s2r["summary"]
    g = s2r["geodesic"]
    p, e, rg = g["p"], g["e"], g["rg_au"]
    E, L = ob.bound_orbit_EL(p, e)
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [1.5, 1]})
    ax = axs[0]
    r = np.geomspace(2.05, 1e5, 4000)               # in GM/c^2
    Vgr = 0.5 * (ob.V2_schwarzschild(r, L) - 1) * 1e5  # (V^2 - 1)/2 -> Newtonian limit; units 1e-5 c^2
    Vn = (-1 / r + L**2 / (2 * r**2)) * 1e5
    ax.plot(r * rg, Vgr, color=ps.GR, lw=1.8, label="general relativity: (V² − 1)/2")
    ax.plot(r * rg, Vn, color=ps.NEWTON, lw=1.4, ls=(0, (5, 2)), label="Newton: −GM/r + L²/2r²")
    eps = 0.5 * (E**2 - 1) * 1e5
    ax.axhline(eps, color=ps.INK, lw=0.9)
    rp, ra = p / (1 + e) * rg, p / (1 - e) * rg
    ax.plot([rp, ra], [eps, eps], "o", color=ps.INK, ms=6, mec=ps.SURFACE, zorder=5)
    ps.annotate(ax, f"pericentre {rp:.0f} AU", (rp, eps), (rp * 0.12, eps - 1.0))
    ps.annotate(ax, f"apocentre {ra:.0f} AU", (ra, eps), (ra * 0.35, eps + 1.6))
    ax.axvline(2 * rg, color=ps.INK, lw=1.2)
    ax.text(2 * rg * 1.15, -3.45, "horizon\n0.078 AU", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_xscale("log")
    ax.set_ylim(-3.6, 2.2)
    ax.set_xlim(0.05, 4e3)
    ax.set_xlabel("distance from Sgr A* (AU), log scale")
    ax.set_ylabel("energy per unit rest mass (10⁻⁵ c²)")
    ax.set_title(f"S2 in its effective potential: L = {sm['L_gr_geom']['median']:.1f} GM/c")
    ax.legend(loc="upper right", fontsize=7.5)
    ins = ax.inset_axes([0.2, 0.08, 0.3, 0.36])
    rr = np.linspace(2.02, 40, 2000)
    ins.plot(rr, 0.5 * (ob.V2_schwarzschild(rr, L) - 1), color=ps.GR, lw=1.4)
    ins.plot(rr, -1 / rr + L**2 / (2 * rr**2), color=ps.NEWTON, lw=1.2, ls=(0, (5, 2)))
    ins.set_xlim(2, 40)
    ins.set_ylim(-5, 160)
    ins.set_xlabel("r (GM/c²)", fontsize=6.5)
    ins.set_ylabel("c²", fontsize=6.5)
    ins.tick_params(labelsize=6)
    ins.set_title("near the hole: GR's barrier turns over", fontsize=6.8, loc="left", pad=3)
    ax.text(0.1, 0.97, f"E = 1 − {sm['binding_gr']['median']:.4e} c² (GR) vs binding {sm['binding_newton']['median']:.4e} "
                        f"(Newton)\nThe two differ by {abs(sm['binding_gr_over_newton_minus_1']['median']) * 1e5:.1f} "
                        f"x 10⁻⁵ of the binding energy; L by {sm['L_gr_over_newton_minus_1']['median'] * 1e4:.1f} x 10⁻⁴.",
            transform=ax.transAxes, fontsize=7.2, color=ps.INK_SECONDARY, va="top")
    ax = axs[1]
    pr = np.array(s2r["samples"]["precession_exact_arcmin"])
    lo, med, hi = np.percentile(pr, [16, 50, 84])
    ax.axvspan(lo, hi, color=ps.BLUE_150, lw=0, label=f"from S2's E and L: {med:.2f}′ ± {(hi - lo) / 2:.2f}′ (H1 posterior)")
    ax.axvline(med, color=ps.BLUE, lw=1.6)
    ax.set_ylim(0, 1)
    ymax = 1.0
    ys = np.linspace(0.85, 0.35, len(lit.GRAVITY_FSP)) * ymax
    for gv, yy in zip(lit.GRAVITY_FSP, ys):
        ref = 12.10
        ax.errorbar([gv["f"] * ref], [yy], xerr=[gv["err"] * ref], fmt="s", color=ps.ORANGE, ms=5, mec=ps.SURFACE,
                    lw=1, capsize=0)
        ax.text(gv["f"] * ref + gv["err"] * ref + 0.3, yy, gv["label"], fontsize=7, va="center",
                color=ps.INK_SECONDARY)
    ax.plot([], [], "s", color=ps.ORANGE, label="measured by GRAVITY (f_SP x 12.1′)")
    ax.set_xlim(7, 22)
    ax.set_xlabel("pericentre advance per orbit (arcmin)")
    ax.set_yticks([])
    ax.set_title("Precession computed from E and L")
    ax.legend(loc="upper right", fontsize=7)
    fig.tight_layout()
    ps.caption(fig, "Posterior over H1's 80,000 samples (joint Keck+VLT GR fit). Exact Schwarzschild orbit integral "
                    f"{sm['precession_exact_arcmin']['median']:.2f}′ ± {sm['precession_exact_arcmin']['sd']:.2f}′ "
                    f"per orbit (6πGM/c²a(1−e²): {sm['precession_6pi_arcmin']['median']:.2f}′). GRAVITY's prediction "
                    "is 12.1′ because their mass and distance are larger (Keck vs GRAVITY distance tension, see H1).",
               y=0.0)
    return ps.save(fig, FIG / "06_s2_potential_precession")


# ---------------------------------------------------------------------------
# 7. Shadows
# ---------------------------------------------------------------------------

def shadows_figure(shd):
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.9))
    cols = {"0.0": ps.INK, "0.5": ps.BLUE, "0.94": ps.ORANGE, "0.998": ps.AQUA}
    for ax, tgt in zip(axs[:2], ("sgra", "m87")):
        if tgt == "sgra":
            tg = shd["sgra"]["theta_g_uas"][0]
            ring, rerr = shd["sgra"]["ring_uas"]
            title = f"Sgr A*: H1's mass and distance (i = {shd['sgra']['inclination_assumed_deg']:.0f}°)"
        else:
            tg = shd["m87"]["by_mass"]["stars, Gebhardt+11"]["theta_g_uas"][0]
            ring, rerr = shd["m87"]["ring_uas"][:2]
            title = f"M87*: stellar-dynamics mass (i = {shd['m87']['inclination_assumed_deg']:.0f}°)"
        th = np.linspace(0, 2 * np.pi, 400)
        for rr in (ring / 2 - rerr / 2, ring / 2 + rerr / 2):
            pass
        ax.fill(np.concatenate([(ring / 2 + rerr / 2) * np.cos(th), (ring / 2 - rerr / 2) * np.cos(th[::-1])]),
                np.concatenate([(ring / 2 + rerr / 2) * np.sin(th), (ring / 2 - rerr / 2) * np.sin(th[::-1])]),
                color=ps.YELLOW, alpha=0.35, lw=0, label=f"ring measured in I5: {ring:.1f} ± {rerr:.1f} µas")
        for s, cv in shd["curves"][tgt].items():
            al, be = np.array(cv["alpha"]), np.array(cv["beta"])
            ac = 0.5 * (al.max() + al.min())
            ax.plot((al - ac) * tg, be * tg, color=cols[s], lw=1.5, label=f"critical curve, a = {s}")
        ax.set_aspect("equal")
        lim = 32 if tgt == "sgra" else 26
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("µas")
        ax.set_ylabel("µas")
        ax.set_title(title, fontsize=9.5)
        ax.legend(loc="lower center", fontsize=6.8, bbox_to_anchor=(0.5, -0.47), ncol=2)
        if tgt == "m87":
            b = shd["m87"]["by_mass"]["stars, Gebhardt+11"]["theta_g_uas"]
            ax.text(0.02, 0.98, f"θg = {b[0]:.2f} +{b[1]:.2f}/−{b[2]:.2f} µas: the non-spinning critical curve\n"
                                f"spans {2 * 27**0.5 * (b[0] - b[2]):.0f}-{2 * 27**0.5 * (b[0] + b[1]):.0f} µas within the mass errors",
                    transform=ax.transAxes, fontsize=6.8, color=ps.INK_SECONDARY, va="top")
    ax = axs[2]
    grid = shd["grid"]
    spins = np.array(grid["spins"])
    D = np.array(grid["mean_diameter_M"])
    incs = grid["incs"]
    ramp = [ps.BLUE_300, ps.BLUE, ps.BLUE_550, ps.BLUE_700]
    for j, inc in enumerate(incs):
        if inc not in (17, 30, 60, 90):
            continue
        ax.plot(spins, D[:, j], color=ramp[[17, 30, 60, 90].index(inc)], lw=1.6, label=f"critical curve, i = {inc}°")
    sg = shd["sgra"]["ring_over_theta_g"]
    ax.errorbar([1.05], [sg[0]], [sg[1]], fmt="o", color=ps.INK, ms=6, mec=ps.SURFACE, capsize=0, clip_on=False)
    ax.text(1.08, sg[0], "Sgr A* ring\n(I5, H1 mass)", fontsize=7, va="center", color=ps.INK_SECONDARY)
    for k_, mk, xx in (("stars, Gebhardt+11", "s", 1.25), ("stars, Liepold+23", "^", 1.45)):
        b = shd["m87"]["by_mass"][k_]
        ax.errorbar([xx], [b["ring_over_theta_g"]], [[b["ring_over_theta_g_lo"]], [b["ring_over_theta_g_up"]]],
                    fmt=mk, color=ps.INK, ms=6, mec=ps.SURFACE, capsize=0, clip_on=False)
        ax.text(xx + 0.03, b["ring_over_theta_g"], "M87 ring\n" + k_.split(", ")[1], fontsize=7, va="center",
                color=ps.INK_SECONDARY)
    ax.axhspan(10.8, 11.6, color=ps.WASH, lw=0, zorder=0)
    ax.text(0.03, 11.2, "EHT's simulated emission rings (α ≈ 10.8-11.6)", fontsize=7, color=ps.INK_SECONDARY,
            va="center")
    ax.set_xlim(0, 1.0)
    ax.set_ylim(9.0, 15.5)
    ax.set_xlabel("spin a")
    ax.set_ylabel("mean diameter / (GM/c²D)")
    ax.set_title("Shadow size barely depends on spin")
    ax.legend(loc="upper left", fontsize=7)
    fig.tight_layout()
    ps.caption(fig, "Critical curves from Bardeen's (1973) photon-orbit constants: L/E from the rotation Killing "
                    "vector, Q/E² from Carter's Killing tensor. Mean diameter 9.66-10.39 GM/c² over all spins <= 0.998 "
                    "and inclinations (EHT: 9.6-10.4). The bright ring sits a few per cent outside the critical curve.",
               y=0.02)
    return ps.save(fig, FIG / "07_shadows")


# ---------------------------------------------------------------------------
# 8. The hidden symmetry: regular versus chaotic orbits
# ---------------------------------------------------------------------------

def hidden_symmetry(ch):
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, n_, col in ((axs[0], "Kerr", ps.BLUE), (axs[1], "ZV", ps.ORANGE)):
        mass = 1.0 if n_ == "Kerr" else 2.0
        secs = ch["sections"][n_]
        for k_, s in enumerate(secs):
            P = np.array(s["points"])
            shade = ps.SEQ_BLUE[3 + (k_ * 9) // max(len(secs), 1)] if n_ == "Kerr" else \
                ["#f7c6b0", "#f1a07d", "#eb6834", "#c4501f", "#8f3714"][k_ * 5 // max(len(secs), 1)]
            ax.plot(P[:, 0] / mass, P[:, 1], ".", ms=1.6, color=shade)
        w = ch["wells"][n_]
        ax.set_xlabel("radius where the orbit crosses the equator (GM/c²)")
        ax.set_ylabel("radial momentum p_r there")
        ax.set_title("Kerr (a = 0.9): has Carter's constant" if n_ == "Kerr" else
                     "Zipoy-Voorhees (delta = 2): no Killing tensor", fontsize=10)
        ax.text(0.97, 0.04, f"L = {w['L'] / mass:.2f} GM/c, E = {w['E']:.4f}", transform=ax.transAxes, fontsize=7.5,
                color=ps.INK_SECONDARY, va="bottom", ha="right")
    ax = axs[2]
    for n_, col, ls in (("Kerr", ps.BLUE, "-"), ("ZV", ps.ORANGE, (0, (5, 2)))):
        for k_, d in enumerate(ch["divergence"][n_]):
            lam, dist = np.array(d["lam"]), np.array(d["dist"])
            per = ch["wells"][n_]["lam_per_crossing"]
            ax.plot(lam / per, np.maximum.accumulate(dist), color=col, lw=1.6 if k_ == 0 else 1.0, ls=ls,
                    label=f"{n_}" + (" (near the barrier)" if k_ == 0 else " (mid-well)"))
    ax.set_yscale("log")
    ax.set_xlabel("time (radial periods)")
    ax.set_ylabel("largest separation so far (orbits started 10⁻⁸ apart)")
    ax.text(0.03, 0.92, "chaos would mean exponential growth: a straight line rising\nsteeply on this log scale. "
                        "Both stay within ~10x of the start.", transform=ax.transAxes, fontsize=7,
            color=ps.INK_SECONDARY, va="top")
    ax.set_title("Do neighbouring orbits fly apart?")
    ax.legend(fontsize=7)
    fig.tight_layout()
    ps.caption(fig, "Poincare sections: one dot per upward crossing of the equatorial plane; orbits with the same "
                    "barrier height (0.995 c² per unit mass) in both spacetimes. Kerr's curves are guaranteed by Carter's "
                    "constant. Zipoy-Voorhees has no Killing tensor (exact count: 4 tensors, all built from g, d_t, d_phi), "
                    "yet these orbits still look regular and separate only linearly: they stay outside ~3 GM/c², where it is "
                    "nearly integrable, and by the KAM theorem chaos then occupies layers too thin to see here.", y=0.0)
    return ps.save(fig, FIG / "08_hidden_symmetry")


# ---------------------------------------------------------------------------
# 9. Noether's theorem, numerically
# ---------------------------------------------------------------------------

def noether_figure(no):
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.4))
    for ax, key, title in ((axs[0], "schwarzschild", "Schwarzschild: an inclined orbit (p = 12, e = 0.5)"),
                           (axs[1], "kerr", "Kerr (a = 0.9): orbit between r = 4 and 9, tilted 40°")):
        d = no[key]
        tau = np.array(d["tau"])
        cols = [ps.BLUE, ps.ORANGE, ps.AQUA, ps.YELLOW, ps.MAGENTA]
        k_ = 0
        for name, vals in d["values"].items():
            if name == "norm":
                continue
            v = np.array(vals)
            dev = np.abs(v - v[0]) / max(np.max(np.abs(v)), 1e-30)
            dev = np.maximum(dev, 1e-17)
            style = (0, (5, 2)) if ("not Killing" in name or "not a symmetry" in name) else "-"
            ax.plot(tau, dev, color=cols[k_ % len(cols)], lw=1.3, ls=style, label=name)
            k_ += 1
        ax.set_yscale("log")
        ax.set_ylim(1e-17, 10)
        ax.set_xlabel("proper time τ (GM/c³)")
        ax.set_ylabel("|Q(τ) − Q(0)| / max |Q|")
        ax.set_title(title, fontsize=9.5)
        ax.legend(loc="center right", fontsize=7)
    fig.tight_layout()
    ps.caption(fig, "Orbits integrated from the geodesic equation alone (Christoffel symbols from SymPy, DOP853, "
                    "rtol 1e-12). Every Killing quantity, and Carter's K from the Killing tensor, stays fixed to the "
                    "integrator's precision. Dashed: quantities that are not symmetries change by order unity: r d/dr, "
                    "and in Kerr even Schwarzschild's total angular momentum L² (a spinning hole is not spherical).", y=0.0)
    return ps.save(fig, FIG / "09_noether_numerically")


# ---------------------------------------------------------------------------
# 10. Thermodynamics
# ---------------------------------------------------------------------------

def thermodynamics(th):
    import thermo
    cen = th["census"]
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), gridspec_kw={"width_ratios": [1.5, 1, 1]})
    ax = axs[0]
    M = np.geomspace(0.5, 1e11, 400)
    ax.plot(M, thermo.kerr_thermo(M, 0.0)["T_K"], color=ps.INK_SECONDARY, lw=1.2, label="non-spinning: T ∝ 1/M")
    ax.axhline(C.T_CMB, color=ps.ORANGE, lw=1.4, ls=(0, (5, 2)), label="today's microwave background, 2.7 K")
    groups = {"merger remnant": ("D", ps.AQUA), "X-ray binary": ("o", ps.BLUE),
              "dormant (Gaia)": ("s", ps.BLUE_550), "supermassive": ("*", ps.INK)}
    for grp, (mk, col) in groups.items():
        rr = [r_ for r_ in cen["rows"] if r_["group"] == grp]
        ax.plot([r_["M_msun"] for r_ in rr], [r_["T_K"] for r_ in rr], mk, color=col, ms=6 if mk != "*" else 9,
                mec=ps.SURFACE, mew=0.6, ls="none", label=grp, zorder=4)
    for name, dx, dy in (("Sgr A*", 1.6, 1.0), ("M87*", 0.12, 3.0), ("Cyg X-1", 1.8, 0.35), ("GW250114 remnant", 1.6, 2.2)):
        r_ = next(r_ for r_ in cen["rows"] if r_["name"] == name)
        ax.annotate(name, (r_["M_msun"], r_["T_K"]), (r_["M_msun"] * dx, r_["T_K"] * dy), fontsize=7,
                    color=ps.INK_SECONDARY)
    q = cen["quasars"]
    if q:
        edges = 10 ** np.array(q["hist_edges"])
        cnt = np.array(q["hist_counts"], float)
        mid = np.sqrt(edges[1:] * edges[:-1])
        Tq = thermo.kerr_thermo(mid, 0.0)["T_K"]
        w = 1 + 6 * cnt / cnt.max()
        ax.fill_between(mid, Tq / w, Tq * w, color=ps.BLUE_150, lw=0, zorder=1, alpha=0.9,
                        label=f"{q['n']:,} SDSS quasars (I3): width = how many")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("mass (solar masses)")
    ax.set_ylabel("Hawking temperature (K)")
    ax.set_title("Every black hole in the atlas is far colder than space")
    ax.legend(loc="lower left", fontsize=6.8)
    ax.set_ylim(1e-20, 1e3)
    # (b) GW250114 areas
    ax = axs[1]
    gw = cen["gw250114"]
    import h5py
    with h5py.File(C.RAW_DIR / C.GW250114["samples_file"], "r") as f:
        ps_ = f["bilby-NRSur7dq4_prod-reweighted/posterior_samples"][:]
    t1 = thermo.kerr_thermo(ps_["mass_1_source"], ps_["a_1"])
    t2 = thermo.kerr_thermo(ps_["mass_2_source"], ps_["a_2"])
    tf = thermo.kerr_thermo(ps_["final_mass_source"], ps_["final_spin"])
    Ai = (t1["area_km2"] + t2["area_km2"]) / 1e5
    Af = tf["area_km2"] / 1e5
    bins = np.linspace(2.1, 4.0, 70)
    ax.hist(Ai, bins=bins, color=ps.BLUE_300, edgecolor=ps.SURFACE, lw=0.3, label="two black holes before")
    ax.hist(Af, bins=bins, color=ps.ORANGE, alpha=0.75, edgecolor=ps.SURFACE, lw=0.3, label="one black hole after")
    ax.set_xlabel("total horizon area (10⁵ km²)")
    ax.set_yticks([])
    ax.set_ylabel("posterior samples")
    ax.set_title("GW250114: the area (entropy) grows")
    ax.legend(loc="upper center", fontsize=7)
    ax.text(0.5, 0.62, f"entropy {gw['S_initial_kB']['median']:.2e} → {gw['S_final_kB']['median']:.2e} k_B\n"
                        f"area ratio {gw['area_ratio']['median']:.2f}; radiated {100 * gw['E_radiated_frac']['median']:.1f}%"
                        f" of the mass\n(the area theorem allows up to {100 * gw['max_radiated_by_area_theorem_frac']['median']:.0f}%)",
            transform=ax.transAxes, fontsize=7, color=ps.INK_SECONDARY, ha="center")
    # (c) lifetimes vs time until the CMB is colder
    ax = axs[2]
    M = np.geomspace(1, 1e11, 300)
    ax.plot(M, thermo.evaporation_time_yr(M), color=ps.BLUE, lw=1.6, label="evaporation time (Page 1976)")
    ff = th["far_future"]["black_holes"]
    names = list(ff)
    Ms = []
    for n_ in names:
        Ms.append(1.0 if n_ == "1 solar mass" else next(r_["M_msun"] for r_ in cen["rows"] if r_["name"] == n_))
    ax.plot(Ms, [ff[n_]["years_from_now"] for n_ in names], "o", color=ps.ORANGE, ms=6, mec=ps.SURFACE,
            label="wait until the CMB is colder")
    ax.axhline(13.8e9, color=ps.MUTED, lw=0.9, ls=(0, (3, 2)))
    ax.text(2, 13.8e9 / 300, "age of the universe (13.8 billion years)", fontsize=7, color=ps.INK_SECONDARY)
    ax.text(2, 4e14, "the CMB cools below each one's temperature\nin 0.3-0.7 trillion years (Planck 2018 LCDM)",
            fontsize=7, color=ps.INK_SECONDARY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("mass (solar masses)")
    ax.set_ylabel("years")
    ax.set_title("When could they start to evaporate?")
    ax.legend(loc="center right", fontsize=7)
    fig.tight_layout()
    ps.caption(fig, "Masses from earlier projects: I4 (X-ray binaries, Gaia BH1-3), H1 (Sgr A*), I5 (M87*), I1 "
                    "(Kormendy & Ho extremes), I3 (quasars), C1 (GW150914), and LVK's GW250114 samples. Spins only where "
                    "published (T falls with spin). Hawking radiation is a prediction: never observed from any black hole.",
               y=0.0)
    return ps.save(fig, FIG / "10_thermodynamics")


# ---------------------------------------------------------------------------
# 11. The expanding universe
# ---------------------------------------------------------------------------

def cosmology_figure(cs):
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
    ax = axs[0]
    tc = cs["tcmb"]
    P = tc["points"]
    z, T = np.array(P["z"]), np.array(P["T"])
    ep, em, es = np.array(P["ep"]), np.array(P["em"]), np.array(P["es"])
    meth = np.array(P["method"])
    zz = np.linspace(0, 7, 200)
    ax.plot(zz, C.T_CMB * (1 + zz), color=ps.GR, lw=1.6, label="conformal Killing vector: T = T₀(1 + z)")
    ax.plot(zz, np.full_like(zz, C.T_CMB), color=ps.NEWTON, lw=1.4, ls=(0, (5, 2)), label="if photon energy were conserved")
    sz = np.array(["SZ" in m_ or "Sunyaev" in m_ for m_ in meth])
    for sel, mk, lab in ((sz, "o", "galaxy clusters (SZ effect)"), (~sz, "s", "molecules and atoms in quasar spectra")):
        ax.errorbar(z[sel], T[sel], [np.hypot(em[sel], es[sel]), np.hypot(ep[sel], es[sel])], fmt=mk, color=ps.INK,
                    ms=4.5, mec=ps.SURFACE, lw=0.8, capsize=0, label=lab)
    ax.set_xlabel("redshift z")
    ax.set_ylabel("CMB temperature (K)")
    ax.set_title("Light cools as space expands (D7's 46 measurements)")
    ax.legend(loc="upper left", fontsize=7)
    ax.text(0.55, 0.12, f"T = T₀(1+z)^(1−β): β = {tc['beta']:.4f} +{tc['plus']:.4f}/−{tc['minus']:.4f}\n"
                        f"χ² = {tc['chi2_expansion']:.1f} (expansion) vs {tc['chi2_static']:.0f} (static), {tc['n']} points",
            transform=ax.transAxes, fontsize=7.2, color=ps.INK_SECONDARY)
    ax = axs[1]
    g = cs["geodesics"]
    ph, ms = g["photon"], g["massive"]
    a_ph, E_ph = np.array(ph["a"]), np.array(ph["E"])
    ax.plot(a_ph, E_ph / E_ph[0], color=ps.BLUE, lw=1.8, label="photon energy E (from recombination)")
    ax.plot(a_ph, np.array(ph["aE"]) / ph["aE"][0], color=ps.BLUE, lw=1.2, ls=(0, (1, 1.5)),
            label="photon a·E: constant (conformal Killing vector)")
    a_m, p_m = np.array(ms["a"]), np.array(ms["p"])
    ax.plot(a_m, p_m / p_m[0], color=ps.ORANGE, lw=1.8, label="galaxy's peculiar momentum p (from z = 10)")
    ax.plot(a_m, np.array(ms["aE"]) / ms["aE"][0], color=ps.ORANGE, lw=1.2, ls=(0, (1, 1.5)),
            label="galaxy's a·E: NOT constant (massive)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("scale factor a (1 = today)")
    ax.set_ylabel("relative to its starting value")
    ax.set_title("Integrated geodesics in Planck 2018 LCDM")
    ax.legend(loc="lower left", fontsize=7)
    ax.text(0.45, 0.85, f"photon: E fell by {ph['E_ratio_start_end']:.0f}x = 1 + z\n"
                        f"galaxy: 0.6c at z = 10 → {ms['v_today'] * 299792:.0f} km/s today",
            transform=ax.transAxes, fontsize=7.2, color=ps.INK_SECONDARY)
    fig.tight_layout()
    ps.caption(fig, "An expanding universe has six Killing vectors (translations and rotations of space) and no time "
                    "translation, so nothing conserves a photon's energy; its conformal Killing vector a d/dt conserves "
                    "a·E for light only. CMB temperatures: D7's compilation (37-expansion-tests), fitted with D7's "
                    "own statistic.", y=0.0)
    return ps.save(fig, FIG / "11_cosmology")


def make_all() -> list:
    written = []
    sym, orb, s2r = load("symbolic"), load("orbits"), load("s2")
    shd, th, cs, no = load("shadows"), load("thermo"), load("cosmology"), load("noether")
    ch = load("chaos")
    written += killing_counts(sym)
    written += killing_fields(sym)
    written += potentials(orb)
    written += spin_figure(orb, th)
    written += s2_conservation(s2r)
    written += s2_potential(s2r)
    written += shadows_figure(shd)
    written += hidden_symmetry(ch)
    written += noether_figure(no)
    written += thermodynamics(th)
    written += cosmology_figure(cs)
    return written
