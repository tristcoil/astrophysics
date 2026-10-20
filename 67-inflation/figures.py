"""Every figure, drawn from the cached analysis in data/derived (run analysis.py first).

    python figures.py            # all figures
    python figures.py 05 08      # only figures whose names start with these
"""

from __future__ import annotations

import json
import math
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

import config as cfg  # noqa: E402
import plotstyle as ps  # noqa: E402

LN10 = math.log(10.0)

#: Diverging blue <-> red with a neutral midpoint, for CMB temperatures.
CMB_CMAP = LinearSegmentedColormap.from_list(
    "cmb", ["#0d366b", "#2a78d6", "#9ec5f4", "#f0efec", "#f3b3b2", "#e34948", "#8f1f1f"])

#: The three graded data combinations: plane key -> (label, hue).
GRADED = {
    "planck18_bk18_ours": ("Planck 2018 + BK18 (our BK18 likelihood)", ps.C_PLANCK),
    "pactlb_bk18": ("ACT DR6: P-ACT-LB + BK18", ps.C_ACT),
    "cmbspa_desi_x_bk18": ("SPT-3G: CMB-SPA + DESI, with our BK18 r", ps.C_SPT),
}
MODEL_SHORT = {"phi23": "φ^(2/3)", "phi1": "φ", "phi2": "φ²", "phi4": "φ⁴",
               "natural": "natural", "starobinsky": "Starobinsky R²",
               "hilltop4": "hilltop φ⁴", "tmodel": "α-attractor T"}


def _load(name: str) -> tuple[dict, dict]:
    with np.load(cfg.DERIVED_DIR / f"{name}.npz") as z:
        arrays = {k: z[k] for k in z.files}
    return arrays, json.loads((cfg.DERIVED_DIR / f"{name}.json").read_text())


def _out(stem: str) -> str:
    return str(cfg.FIGURE_DIR / stem)


# ---------------------------------------------------------------------------
# 01 The horizon problem on the real CMB sky
# ---------------------------------------------------------------------------


def _rgba(img: np.ndarray, vmin=-300, vmax=300) -> np.ndarray:
    out = CMB_CMAP(Normalize(vmin, vmax)(np.nan_to_num(img, nan=0.0, posinf=0.0, neginf=0.0)))
    outside = ~np.isfinite(img) & (img < 0) if img.dtype.kind == "f" else np.zeros(img.shape, bool)
    out[np.isnan(img)] = (0.86, 0.86, 0.84, 1.0)            # Planck's mask: light grey
    out[np.isinf(img)] = (0, 0, 0, 0)                        # outside the Mollweide ellipse
    out[outside & np.isinf(img)] = (0, 0, 0, 0)
    return out


def _circle_xy(lon0, lat0, rad_deg, centre, n=72):
    import horizons as hz

    lat1, lon1, d = map(math.radians, (lat0, lon0, rad_deg))
    b = np.linspace(0, 2 * math.pi, n)
    lat2 = np.arcsin(math.sin(lat1) * math.cos(d) + math.cos(lat1) * math.sin(d) * np.cos(b))
    lon2 = lon1 + np.arctan2(np.sin(b) * math.sin(d) * math.cos(lat1),
                             math.cos(d) - math.sin(lat1) * np.sin(lat2))
    return hz.gnomonic_forward(np.degrees(lon2), np.degrees(lat2), *centre)


def fig_01_cmb_horizons():
    import horizons as hz

    a, s = _load("horizons")
    hzn, sky = s["horizon"], s["sky"]
    theta = hzn["theta_H_deg"]
    fig = plt.figure(figsize=(12.4, 5.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.72, 1], wspace=0.06)
    ax = fig.add_subplot(gs[0])
    ext = a["sky_moll_extent"]
    ax.imshow(_rgba(a["sky_mollweide"]), origin="lower", extent=ext, interpolation="nearest")
    ox, oy = a["sky_moll_outline"]
    ax.plot(ox, oy, color=ps.INK, lw=1.1)
    ax.annotate("zoom →", (np.max(ox), np.mean(oy)), xytext=(8, 0), textcoords="offset points",
                fontsize=8, color=ps.INK_SECONDARY, va="center")
    ax.set_axis_off()
    ax.set_title("The CMB sky (Planck 2018): uniform to 1 part in 30,000", loc="left")
    # the zoom, with the patches that could have exchanged light without inflation
    az = fig.add_subplot(gs[1])
    half = sky["zoom_deg"] / 2
    centre = tuple(a["sky_zoom_center"])
    az.imshow(_rgba(a["sky_zoom"]), origin="lower", extent=(-half, half, -half, half),
              interpolation="bilinear")
    rad = theta / 2
    step = 2 * rad
    rows = int(np.ceil(2 * half / (step * math.sqrt(3) / 2))) + 2
    cols = int(np.ceil(2 * half / step)) + 2
    for i in range(-rows // 2, rows // 2 + 1):
        for j in range(-cols // 2, cols // 2 + 1):
            xc = j * step + (step / 2 if i % 2 else 0.0)
            yc = i * step * math.sqrt(3) / 2
            if abs(xc) > half + rad or abs(yc) > half + rad:
                continue
            lon, lat = hz.gnomonic_inverse(np.array(xc), np.array(yc), *centre)
            cx, cy = _circle_xy(float(lon), float(lat), rad, centre)
            az.plot(cx, cy, color="white", lw=0.45, alpha=0.85)
    # one patch picked out, with the larger disc of 'no common past at all'
    lon, lat = hz.gnomonic_inverse(np.array(0.0), np.array(0.0), *centre)
    cx, cy = _circle_xy(float(lon), float(lat), rad, centre)
    az.plot(cx, cy, color=ps.INK, lw=1.6)
    cx, cy = _circle_xy(float(lon), float(lat), theta, centre)
    az.plot(cx, cy, color=ps.INK, lw=1.0, ls="--")
    # the full Moon, for scale
    mx, my = -half + 1.4, -half + 1.4
    lon_m, lat_m = hz.gnomonic_inverse(np.array(mx), np.array(my), *centre)
    cx, cy = _circle_xy(float(lon_m), float(lat_m), 0.26, centre)
    az.fill(cx, cy, color="white", alpha=0.95, lw=0)
    az.text(mx, my + 0.75, "full Moon", ha="center", va="bottom", fontsize=7, color=ps.INK,
            bbox=dict(fc="white", ec="none", alpha=0.8, pad=1))
    az.set_xlim(half, -half)  # longitude increases to the left, as on the sky
    az.set_ylim(-half, half)
    az.set_xticks([])
    az.set_yticks([])
    for sp in az.spines.values():
        sp.set_visible(True)
        sp.set_color(ps.AXIS)
    az.set_title(f"{sky['zoom_deg']:.0f}° zoom: {theta:.2f}° causal patches", loc="left")
    az.text(0.02, 0.98,
            "solid circle: every pair of points inside could have\n"
            "exchanged light by the time the CMB was released\n"
            "dashed: points outside share no past at all",
            transform=az.transAxes, va="top", fontsize=7, color=ps.INK,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=2))
    cax = fig.add_axes([0.13, 0.17, 0.32, 0.025])
    sm = plt.cm.ScalarMappable(norm=Normalize(-300, 300), cmap=CMB_CMAP)
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("temperature − 2.7255 K  [µK]", fontsize=8, color=ps.INK_SECONDARY)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=7, colors=ps.MUTED)
    ps.caption(fig,
               f"Without inflation, light could have travelled only η* = {hzn['eta_star_mpc']:.0f} Mpc "
               f"(comoving) before the CMB was released, at a distance of {hzn['D_M_star_mpc']:,.0f} Mpc: "
               f"θ_H = {theta:.2f}° on the sky. The sky holds ≈{hzn['patches_exchange']:,.0f} such patches "
               f"(≈{hzn['patches_common_past']:,.0f} discs of radius θ_H with no shared past). "
               "Planck SMICA, mean of the half-mission maps, 15′ smoothing, Galactic coordinates; "
               "grey: Planck common mask. Horizon from CAMB with Planck 2018 parameters.", y=0.05)
    return ps.save(fig, _out("01_cmb_horizons"))


# ---------------------------------------------------------------------------
# 02 The shrinking Hubble sphere and the flatness problem
# ---------------------------------------------------------------------------


def fig_02_hubble_sphere():
    a, s = _load("horizons")
    hist = s["history"]
    la = a["ln_a"] / LN10
    lr = a["ln_inv_aH_mpc"] / LN10
    ph = a["phase"]
    n_show = 70.0  # e-folds of inflation shown before its end
    la_end = math.log10(hist["a_end"])
    keep = la >= la_end - n_show / LN10
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10.2, 8.0), sharex=True,
                                   gridspec_kw={"hspace": 0.07})
    om = 0.3153
    a_eq = hist["omega_r"] / om
    a_acc = (om / (2 * (1 - om - hist["omega_r"]))) ** (1 / 3)
    x0, x1 = la_end - n_show / LN10, math.log10(30)
    bands = (("inflation", x0, la_end, ps.lighten(ps.BLUE, 0.85)),
             ("radiation", la_end, math.log10(a_eq), ps.lighten(ps.YELLOW, 0.88)),
             ("matter", math.log10(a_eq), math.log10(a_acc), ps.lighten(ps.NEUTRAL, 0.8)),
             ("dark\nenergy", math.log10(a_acc), x1, ps.lighten(ps.VIOLET, 0.86)))
    for axx in (ax1, ax2):
        for _, lo, hi, col in bands:
            axx.axvspan(lo, hi, color=col, lw=0, zorder=0)
    for name, lo, hi, _ in bands:
        if name.startswith("dark"):
            ax1.text(hi - 0.2, 25.5, "dark energy", ha="right", va="top", fontsize=8,
                     color=ps.INK_SECONDARY, rotation=90)
            continue
        ax1.text(0.5 * (lo + hi), 25.5, name, ha="center", va="top", fontsize=8,
                 color=ps.INK_SECONDARY)
    # top: comoving Hubble radius
    inf_ = keep & (ph == 0)
    ax1.plot(la[inf_], lr[inf_], color=ps.INK, lw=1.8)
    ax1.plot(la[ph > 0], lr[ph > 0], color=ps.INK, lw=1.8, label="our model universe")
    ax1.plot(a["noinf_ln_a"] / LN10, a["noinf_ln_inv_aH_mpc"] / LN10, color=ps.NEUTRAL, lw=1.4,
             ls="--", label="hot Big Bang without inflation (back to the Planck temperature)")
    xn, yn = a["noinf_ln_a"][0] / LN10, a["noinf_ln_inv_aH_mpc"][0] / LN10
    ps.annotate(ax1, "without inflation, at the Planck temperature\n"
                     f"the Hubble sphere was {10 ** yn:.0e} Mpc (comoving)",
                (xn, yn), (xn - 21, yn - 6), fontsize=7.5)
    h_today = math.log10(s["horizon"]["hubble_radius_today_mpc"])
    pivot = math.log10(1 / cfg.K_PIVOT)
    for y, txt, dy, va in ((h_today, "today's Hubble radius, 4.45 Gpc", 0.6, "bottom"),
                           (pivot, "the pivot scale, 1/k* = 20 Mpc", -0.6, "top")):
        ax1.axhline(y, color=ps.INK_SECONDARY, lw=0.7)
        ax1.text(la_end + 1.0, y + dy, txt, fontsize=7.5, color=ps.INK_SECONDARY, va=va)
    n_star = hist["n_star"]
    x_exit = la_end - n_star / LN10
    ax1.plot([x_exit], [pivot], "o", ms=6, mfc=ps.SURFACE, mec=ps.INK, mew=1.2, zorder=5)
    ps.annotate(ax1, f"the pivot leaves the Hubble sphere\n{n_star:.1f} e-folds before inflation ends",
                (x_exit, pivot), (x_exit + 3, pivot - 9), fontsize=7.5)
    hot = ph >= 1
    i_re = np.argmin(np.abs(lr[hot] - pivot) + 1000 * (la[hot] < la_end + 1))
    x_re = la[hot][i_re]
    ax1.plot([x_re], [pivot], "o", ms=6, mfc=ps.INK, mec=ps.INK, zorder=5)
    ps.annotate(ax1, "... and comes back in", (x_re, pivot), (x_re - 9, pivot + 6), fontsize=7.5)
    n_h = hist["n_min_horizon"]
    x_h = la_end - n_h / LN10
    ax1.plot([x_h], [h_today], "s", ms=5, color=ps.INK, zorder=5)
    ps.annotate(ax1, f"everything we see today was inside one Hubble sphere\n"
                     f"{n_h:.1f} e-folds before the end (the horizon problem, solved)",
                (x_h, h_today), (x_h + 4.5, h_today + 7), fontsize=7.5)
    ax1.set_ylabel("comoving Hubble radius 1/(aH)  [log₁₀ Mpc]")
    ax1.set_ylim(-36, 26)
    ax1.legend(loc="lower right", fontsize=7.5)
    ax1.set_title("Inflation shrinks the Hubble sphere; everything since has grown it", loc="left")
    # bottom: |Omega - 1|
    od = a["omega_dev"]
    ok = keep & (ph == 0) & (od <= 1.0)
    ax2.plot(la[ok], np.log10(od[ok]), color=ps.INK, lw=1.8)
    ax2.plot(la[ph >= 1], np.log10(od[ph >= 1]), color=ps.INK, lw=1.8)
    ax2.plot(a["noinf_ln_a"] / LN10, np.log10(a["noinf_omega_dev"]), color=ps.NEUTRAL, lw=1.4,
             ls="--")
    bound = hist["omega_k_bound"]
    ax2.plot([0.0], [math.log10(bound)], "o", ms=6, color=ps.INK, zorder=5)
    ps.annotate(ax2, f"today: |Ω − 1| < {bound:.4f} (95%, O1: DESI BAO + CMB)",
                (0.0, math.log10(bound)), (-30, 0.5), fontsize=7.5)
    t_pl = 1.220890e19
    pl = s["omega_dev_planck_temperature"]
    xpl = a["noinf_ln_a"][0] / LN10
    ax2.plot([xpl], [math.log10(pl)], "o", ms=5, mfc=ps.SURFACE, mec=ps.NEUTRAL, mew=1.2, zorder=5)
    ps.annotate(ax2, f"at the Planck temperature ({t_pl:.1e} GeV): {pl:.0e}\n"
                     "— fine-tuned to 1 part in 10⁶², if no inflation", (xpl, math.log10(pl)),
                (xpl - 26, math.log10(pl) + 8), fontsize=7.5)
    temps = a["temperature_gev"]
    for T, lab in ((160.0, "electroweak, 160 GeV"), (0.1565, "quarks → hadrons"),
                   (1e-3, "nucleosynthesis, 1 MeV")):
        i = np.nanargmin(np.abs(np.log(np.where(np.isfinite(temps), temps, np.inf)) - math.log(T)))
        ax2.plot([la[i]], [math.log10(od[i])], "o", ms=4, color=ps.INK_SECONDARY, zorder=5)
        ax2.annotate(f"{lab}: {od[i]:.0e}", (la[i], math.log10(od[i])), xytext=(4, -10),
                     textcoords="offset points", fontsize=7, color=ps.INK_SECONDARY)
    n_f = hist["n_min_flatness"]
    x_f = la_end - n_f / LN10
    ax2.plot([x_f], [0.0], "s", ms=5, color=ps.INK, zorder=5)
    ps.annotate(ax2, f"|Ω − 1| ≈ 1 here is enough:\n{n_f:.1f} e-folds of inflation drive it to 10⁻⁵⁵",
                (x_f, 0.0), (x_f + 4, -16), fontsize=7.5)
    ax2.set_ylim(-68, 6)
    ax2.set_xlim(x0, x1)
    ax2.set_ylabel("|Ω − 1|  [log₁₀]")
    ax2.set_xlabel("scale factor a  [log₁₀; a = 1 today]")
    ax2.set_title("The flatness problem: |Ω − 1| grows as (aH)⁻² whenever the expansion decelerates",
                  loc="left")
    ps.caption(fig,
               "A labelled calculation. Inflation: the Starobinsky model solved exactly, with instant "
               f"reheating at T = {hist['T_reh_gev']:.1e} GeV (N* = {n_star:.1f}); afterwards the Friedmann "
               "equation with Planck 2018 parameters and the Standard Model's g*(T) (free-particle gas with a "
               "smooth quark-hadron crossover). |Ω − 1| = |Ω_k|(a₀H₀/aH)², anchored to today's 95% bound.",
               y=0.02)
    return ps.save(fig, _out("02_hubble_sphere"))


# ---------------------------------------------------------------------------
# 03 The potentials and where the field sits
# ---------------------------------------------------------------------------

REPRESENTATIVE = {"natural": 7.6, "hilltop4": 20.0, "tmodel": 1.0}


def fig_03_potentials():
    import inflation as inf

    mod, ms = _load("models")
    fig, axes = plt.subplots(2, 4, figsize=(12.6, 6.0))
    for ax, name in zip(axes.ravel(), cfg.MODELS):
        spec = cfg.MODELS[name]
        params = mod[f"{name}__param"]
        if "param" in spec:
            i = int(np.argmin(np.abs(np.log(params / REPRESENTATIVE[name]))))
            par = float(params[i])
        else:
            i, par = 0, None
        m = inf.model_from_spec(spec, par)
        phi_s = mod[f"{name}__phi_star"][i]
        phi_e = mod[f"{name}__phi_end"][i][0]
        nst = mod[f"{name}__n_star"][i]
        lo = min(phi_e, phi_s.min())
        hi = max(phi_e, phi_s.max())
        span = hi - lo
        if spec["kind"] == "natural":
            grid = np.linspace(0, math.pi * par, 400)
        elif spec["kind"] == "hilltop":
            grid = np.linspace(0, par * 0.999, 400)
        else:
            grid = np.linspace(max(0.0, lo - 0.6 * span), hi + 0.35 * span, 400)
        v = np.array([m.v(x)[0] for x in grid])
        vmax = v.max()
        ax.plot(grid, v / vmax, color=ps.INK, lw=1.6)
        j60, j50 = int(np.argmin(np.abs(nst - 60))), int(np.argmin(np.abs(nst - 50)))
        for j, mk, lab in ((j60, "o", "N = 60"), (j50, "o", "N = 50")):
            x = phi_s[j]
            ax.plot([x], [m.v(x)[0] / vmax], mk, ms=6, mfc=ps.SURFACE if j == j60 else ps.INK,
                    mec=ps.INK, mew=1.2, zorder=5)
        ax.plot([phi_e], [m.v(phi_e)[0] / vmax], "x", ms=7, mew=1.6, color=ps.INK, zorder=5)
        title = MODEL_SHORT[name]
        if par is not None:
            title += f"  ({spec['param'].replace('alpha', 'α').replace('mu', 'μ')} = {par:.3g}" + \
                     (" M_Pl)" if spec["param"] != "alpha" else ")")
        ax.set_title(title, loc="left", fontsize=9.5)
        dphi = abs(phi_s[j60] - phi_e)
        # the empty corner: upper left where V rises to the right, lower left where it falls
        rising = spec["kind"] not in ("natural", "hilltop")
        ax.text(0.03, 0.95 if rising else 0.05, f"Δφ (N* = 60 → end) = {dphi:.2f} M_Pl",
                transform=ax.transAxes, ha="left", va="top" if rising else "bottom", fontsize=7,
                color=ps.INK_SECONDARY)
        ax.set_ylim(-0.04, 1.08)
        ax.set_xlabel("φ  [M_Pl]", fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel("V / V_max")
    handles = [Line2D([], [], marker="o", ls="", mfc=ps.SURFACE, mec=ps.INK, ms=6,
                      label="field when the pivot scale left the horizon, N* = 60"),
               Line2D([], [], marker="o", ls="", mfc=ps.INK, mec=ps.INK, ms=6, label="... N* = 50"),
               Line2D([], [], marker="x", ls="", color=ps.INK, ms=7, mew=1.6,
                      label="end of inflation (ε = 1)")]
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Eight shapes for the inflaton's hill: the field rolls toward the minimum",
                 x=0.06, ha="left", y=1.0)
    fig.tight_layout(rect=(0, 0.04, 1, 0.98))
    return ps.save(fig, _out("03_potentials"))


# ---------------------------------------------------------------------------
# 04 Quantum fluctuations frozen in: the Mukhanov-Sasaki equation
# ---------------------------------------------------------------------------


def fig_04_mukhanov_sasaki():
    md, msum = _load("modes")
    mod, ms = _load("models")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.4, 4.8), gridspec_kw={"width_ratios": [1.25, 1]})
    x = md["starobinsky__n_since_exit"]
    amp = md["starobinsky__sqrt_PR"] * 1e5 / math.sqrt(msum["starobinsky"]["final_PR"])
    # normalise so the frozen amplitude equals the measured sqrt(A_s) ~ 4.6e-5
    a_s = math.exp(cfg.LN10_10_AS) * 1e-10
    scale = math.sqrt(a_s) / math.sqrt(msum["starobinsky"]["final_PR"])
    re = np.abs(md["starobinsky__re_sqrt_PR"] * scale)
    ab = md["starobinsky__sqrt_PR"] * scale
    at = md["starobinsky__sqrt_Pt"] * scale
    ax1.plot(x, re * 1e5, color=ps.lighten(ps.BLUE, 0.4), lw=0.7,
             label="curvature perturbation: |real part| (oscillation)")
    ax1.plot(x, ab * 1e5, color=ps.BLUE, lw=1.9, label="curvature perturbation: amplitude")
    ax1.plot(x, at * 1e5, color=ps.INK, lw=1.4, ls="--", label="gravitational wave: amplitude")
    ax1.axvline(0, color=ps.INK_SECONDARY, lw=0.8)
    ax1.set_yscale("log")
    ax1.set_ylim(0.05, 3e3)
    ax1.set_xlim(x[0], x[-1])
    ax1.set_xlabel("e-folds since the mode left the Hubble sphere (k = aH)")
    ax1.set_ylabel("√P per mode  [×10⁻⁵]")
    ax1.text(-5.4, 0.3, "inside the horizon: a quantum oscillation,\n"
                        "its amplitude falling as 1/a", fontsize=7.5, color=ps.INK_SECONDARY)
    ax1.annotate(f"outside: frozen at √A_s = {math.sqrt(a_s) * 1e5:.1f}×10⁻⁵ —\n"
                 "the seed of the CMB's ripples", (3.0, ab[-1] * 1e5), xytext=(1.2, 120),
                 textcoords="data", fontsize=7.5, color=ps.INK_SECONDARY,
                 arrowprops={"arrowstyle": "-", "color": ps.INK_SECONDARY, "lw": 0.7})
    ax1.annotate(f"gravitational waves frozen at √(r A_s): r = "
                 f"{(at[-1] / ab[-1]) ** 2:.4f}", (3.0, at[-1] * 1e5), xytext=(0.6, 0.09),
                 textcoords="data", fontsize=7.5, color=ps.INK_SECONDARY,
                 arrowprops={"arrowstyle": "-", "color": ps.INK_SECONDARY, "lw": 0.7})
    ax1.legend(loc="upper right", fontsize=7.3)
    ax1.set_title("One fluctuation, followed exactly (Starobinsky, pivot at N* = 55)", loc="left")
    # right: how good is slow roll?
    names = list(cfg.MODELS)
    rows = []
    for name in names:
        spec = cfg.MODELS[name]
        params = mod[f"{name}__param"]
        if "param" in spec:
            i = int(np.argmin(np.abs(np.log(params / REPRESENTATIVE[name]))))
        else:
            i = 0
        nst = mod[f"{name}__n_star"][i]
        j = int(np.argmin(np.abs(nst - 55)))
        ex = mod[f"{name}__ns"][i, j]
        rows.append((MODEL_SHORT[name], (mod[f"{name}__ns_sr1"][i, j] - ex) * 1e3,
                     (mod[f"{name}__ns_sr2"][i, j] - ex) * 1e3))
    y = np.arange(len(rows))[::-1]
    sig = 0.0033e3
    ax2.axvline(sig, color=ps.INK_SECONDARY, lw=0.9)
    ax2.text(sig - 0.05, len(rows) - 0.75, "today's best 1σ\nerror on n_s\n(P-ACT-LB)",
             ha="right", va="top", fontsize=7, color=ps.INK_SECONDARY)
    ax2.axvline(0, color=ps.AXIS, lw=0.8)
    for yy, (lab, d1, d2) in zip(y, rows):
        ax2.plot([d1], [yy], "o", ms=7, mfc=ps.SURFACE, mec=ps.INK, mew=1.3)
        ax2.plot([d2], [yy], "s", ms=5, color=ps.INK)
    ax2.set_yticks(y)
    ax2.set_yticklabels([r[0] for r in rows], fontsize=8, color=ps.INK)
    ax2.set_xlabel("slow-roll n_s − exact n_s   [×10⁻³]  (N* = 55)")
    ax2.legend(handles=[Line2D([], [], marker="o", ls="", mfc=ps.SURFACE, mec=ps.INK, ms=7,
                               label="first-order slow roll (textbook formulas)"),
                        Line2D([], [], marker="s", ls="", color=ps.INK, ms=5,
                               label="second-order (Hubble-flow) slow roll")],
               loc="lower right", fontsize=7.2)
    ax2.set_xlim(-0.4, 3.6)
    ax2.set_ylim(-1.6, len(rows) - 0.4)
    ax2.set_title("How good are the slow-roll formulas?", loc="left")
    fig.tight_layout()
    return ps.save(fig, _out("04_mukhanov_sasaki"))


# ---------------------------------------------------------------------------
# 05 The keeper: the n_s - r plane
# ---------------------------------------------------------------------------


def _plane(c: dict, key: str):
    return c["plane_x"], c["plane_y"], c[f"plane__{key}__p"].astype(float), \
        c[f"plane__{key}__levels"]


def _draw_region(ax, x, y, p, levels, color, *, fill=True, lw=1.0, ls="-", z=2, alpha=0.75,
                 log_extend=False):
    if log_extend:  # continue the density below the grid's first r row (flat in r there)
        ylo = np.geomspace(1e-6, y[0], 40, endpoint=False)
        y = np.concatenate([ylo, y])
        p = np.concatenate([np.repeat(p[:, :1], ylo.size, axis=1), p], axis=1)
    light, dark = ps.RAMPS.get(color, (ps.lighten(color, 0.7), ps.lighten(color, 0.35)))
    l95, l68 = levels[1], levels[0]
    if fill:
        ax.contourf(x, y, p.T, levels=[l95, l68, p.max() * 2], colors=[light, dark], alpha=alpha,
                    zorder=z, antialiased=True)
    ax.contour(x, y, p.T, levels=[l95, l68], colors=[color, color], linewidths=[lw, lw],
               linestyles=[ls, ls], zorder=z + 0.5)


def _family_band(mod, name):
    ns, r = mod[f"{name}__ns"], mod[f"{name}__r"]
    bx = np.concatenate([ns[:, 0], ns[::-1, -1]])
    by = np.concatenate([r[:, 0], r[::-1, -1]])
    return bx, by


def _draw_models(ax, mod, *, logy=False, labels=True):
    styles = {"natural": dict(ls="-.", hatch=None), "hilltop4": dict(ls="--", hatch=None),
              "tmodel": dict(ls=":", hatch=None)}
    for name, st in styles.items():
        bx, by = _family_band(mod, name)
        ax.fill(bx, by, color=ps.NEUTRAL_LIGHT, alpha=0.45, lw=0, zorder=3)
        ax.plot(np.append(bx, bx[0]), np.append(by, by[0]), color=ps.INK_SECONDARY, lw=0.9,
                ls=st["ls"], zorder=3.2)
    for name in ("phi23", "phi1", "phi2", "phi4", "starobinsky"):
        ns, r = mod[f"{name}__ns"][0], mod[f"{name}__r"][0]
        ax.plot(ns, r, color=ps.INK, lw=2.0, zorder=6, solid_capstyle="round")
        ax.plot([ns[0]], [r[0]], "o", ms=5.5, mfc=ps.SURFACE, mec=ps.INK, mew=1.2, zorder=7)
        ax.plot([ns[-1]], [r[-1]], "o", ms=5.5, mfc=ps.INK, mec=ps.INK, zorder=7)


def fig_05_ns_r_plane():
    c, cs = _load("constraints")
    mod, ms = _load("models")
    fig, (ax, axl) = plt.subplots(1, 2, figsize=(13.2, 6.2), gridspec_kw={"width_ratios": [1.25, 1]})
    for axx, logy in ((ax, False), (axl, True)):
        x, y, p, lev = _plane(c, "planck18")
        _draw_region(axx, x, y, p, lev, ps.NEUTRAL, fill=False, lw=1.0, ls="-", z=1.5,
                     log_extend=logy)
        for z, (key, (label, col)) in enumerate(GRADED.items()):
            x, y, p, lev = _plane(c, key)
            # SPT's plane is a product approximation (and aqua is the low-contrast hue):
            # an outline only, so the two full joint fits stay readable underneath
            outline = key == "cmbspa_desi_x_bk18"
            _draw_region(axx, x, y, p, lev, col, z=2 + z, alpha=0.62, log_extend=logy,
                         fill=not outline, lw=1.6 if outline else 1.0)
        _draw_models(axx, mod, logy=logy)
        d3 = cs["d3"]
        yy = 0.006 if not logy else 1.5e-5
        axx.errorbar([d3["value"]], [yy], xerr=[d3["error"]], fmt="D", ms=5, color=ps.INK,
                     mfc=ps.SURFACE, mew=1.2, capsize=0, lw=1.2, zorder=8)
    ax.set_xlim(0.940, 0.990)
    ax.set_ylim(0.0, 0.20)
    ax.set_xlabel("spectral index n_s  (tilt of the primordial ripples; 1 = scale-invariant)")
    ax.set_ylabel("tensor-to-scalar ratio r  (gravitational waves; at k = 0.05 Mpc⁻¹)")
    axl.set_yscale("log")
    axl.set_xlim(0.950, 0.985)
    axl.set_ylim(1e-5, 0.3)
    axl.set_xlabel("spectral index n_s")
    axl.set_ylabel("r  (log scale)")
    # direct labels: models
    lab = dict(fontsize=7.5, color=ps.INK, zorder=9,
               bbox=dict(fc=ps.SURFACE, ec="none", alpha=0.75, pad=0.6))
    for name, (dx, dy) in {"phi4": (0.0, 0.006), "phi2": (0.0016, 0.004), "phi1": (0.0016, 0.0),
                           "phi23": (0.0016, -0.002)}.items():
        ns, r = mod[f"{name}__ns"][0], mod[f"{name}__r"][0]
        if ns[-1] + dx < 0.989 and r[-1] + dy < 0.2:
            ax.text(ns[-1] + dx, r[-1] + dy, MODEL_SHORT[name], **lab)
    ns, r = mod["starobinsky__ns"][0], mod["starobinsky__r"][0]
    ps.annotate(ax, "Starobinsky R²\n(and Higgs inflation)", (ns[-1], r[-1]), (0.9815, 0.010),
                fontsize=7.5, color=ps.INK)
    n4 = mod["phi4__ns"][0]
    ax.annotate(f"φ⁴ is off the top: r = {mod['phi4__r'][0].min():.2f}–{mod['phi4__r'][0].max():.2f}",
                (0.5 * (n4.min() + n4.max()), 0.199), xytext=(0.9405, 0.188), fontsize=7.5,
                color=ps.INK, arrowprops={"arrowstyle": "->", "color": ps.INK, "lw": 0.8})
    ps.annotate(ax, "natural inflation", (0.952, 0.07), (0.9425, 0.105), fontsize=7.5, color=ps.INK)
    ps.annotate(ax, "hilltop φ⁴", (0.9735, 0.055), (0.981, 0.085), fontsize=7.5, color=ps.INK)
    ps.annotate(ax, "α-attractors (T)", (0.9605, 0.045), (0.9445, 0.03), fontsize=7.5, color=ps.INK)
    d3v = cs["d3"]["value"]
    ps.annotate(ax, "D3: our own n_s\n(Planck TT, r = 0)", (d3v - 0.003, 0.006), (0.9435, 0.012),
                fontsize=7, color=ps.INK_SECONDARY)
    ns, r = mod["starobinsky__ns"][0], mod["starobinsky__r"][0]
    axl.text(ns[-1] + 0.0008, r[-1] * 1.25, "Starobinsky", **lab)
    axl.text(0.9612, 4.5e-5, "α-attractors", **lab)
    axl.text(0.9700, 0.12, "φ", **lab)
    axl.text(0.9790, 0.05, "φ^(2/3)", **lab)
    ps.annotate(axl, "hilltop φ⁴", (0.9565, 0.0016), (0.9505, 0.0045), fontsize=7.5, color=ps.INK)
    ps.annotate(axl, "natural", (0.9605, 0.055), (0.951, 0.1), fontsize=7.5, color=ps.INK)
    handles = [Patch(fc=ps.RAMPS[col][1], ec=col, label=label) if key != "cmbspa_desi_x_bk18"
               else Patch(fc="none", ec=col, lw=1.6, label=label + " (outline)")
               for key, (label, col) in GRADED.items()]
    handles.insert(0, Patch(fc="none", ec=ps.NEUTRAL, label="Planck 2018 alone (no B-mode data)"))
    handles += [Line2D([], [], marker="o", ls="", mfc=ps.SURFACE, mec=ps.INK, ms=6,
                       label="model at N* = 50 e-folds"),
                Line2D([], [], marker="o", ls="", mfc=ps.INK, mec=ps.INK, ms=6, label="... at N* = 60"),
                Patch(fc=ps.NEUTRAL_LIGHT, ec=ps.INK_SECONDARY, label="model families, N* = 50–60"),
                Line2D([], [], marker="D", ls="", mfc=ps.SURFACE, mec=ps.INK, ms=5,
                       label="D3: our own n_s (Planck TT, τ fixed)")]
    ax.legend(handles=handles, loc="upper right", fontsize=7.4, framealpha=0.9, frameon=True,
              facecolor=ps.SURFACE, edgecolor="none")
    ax.set_title("Inflation on trial: models against the CMB (68% and 95% regions)", loc="left")
    axl.set_title("The same on a log scale: where plateau models live", loc="left")
    ps.caption(fig,
               "Contours from the published MCMC chains: Planck 2018 (PLA grid); Planck 2018 + BK18 = Planck's "
               "r-free chain weighted by our own BK18 likelihood; ACT DR6.02 P-ACT-LB + BK18 (LAMBDA); SPT-3G "
               "D1 CMB-SPA + DESI DR2 (LAMBDA, 2026 bug-fix chains; no r fit, so its n_s posterior × our BK18 r "
               "likelihood, an approximation checked on ACT). Models solved exactly (Mukhanov-Sasaki) for "
               "N* = 50-60 e-folds; families span natural f = 3-100 M_Pl, hilltop μ = 5-100 M_Pl, α = 0.01-30.",
               y=0.0)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    return ps.save(fig, _out("05_ns_r_plane"))


# ---------------------------------------------------------------------------
# 06 n_s across data sets
# ---------------------------------------------------------------------------

NS_ROWS = (  # (chain key or special, label, colour key)
    ("planck18_tt_lcdm", "Planck 2018 TT + lowE", None),
    ("planck18_lcdm", "Planck 2018 TT,TE,EE + lowE + lensing", None),
    ("ours", "Planck 2018 + BK18, r free (ours)", "planck18_bk18_ours"),
    ("act_planck_lcdm_camb", "Planck (ACT's run: plik-lite, Sroll2 low-ℓ EE)", None),
    ("act_actlite_lcdm_camb", "ACT DR6 alone", None),
    ("act_p-actlite_lcdm_camb", "P-ACT (Planck + ACT)", None),
    ("act_p-actlite-l-b_lcdm_camb", "P-ACT-LB (+ lensing + DESI DR1 BAO)", "pactlb_bk18"),
    ("spt3g", "SPT-3G D1 alone", None),
    ("cmbspa", "CMB-SPA (SPT + Planck + ACT)", None),
    ("cmbspa_desi", "CMB-SPA + DESI DR2", "cmbspa_desi_x_bk18"),
    ("d3", "D3: our own fit (Planck TT power spectrum, τ fixed)", None),
)


def fig_06_ns_datasets():
    c, cs = _load("constraints")
    mod, ms = _load("models")
    fig, ax = plt.subplots(figsize=(10.4, 6.0))
    staro = ms["models"]["starobinsky"]
    lo, hi = staro["table"]["N50"]["ns"], staro["table"]["N60"]["ns"]
    ax.axvspan(lo, hi, color=ps.NEUTRAL_LIGHT, alpha=0.6, lw=0)
    ax.text(0.5 * (lo + hi), len(NS_ROWS) - 0.35, "Starobinsky\nN* = 50–60", ha="center",
            va="bottom", fontsize=7.5, color=ps.INK_SECONDARY)
    inst = staro["reheating"]["instant"]
    ax.axvline(inst["ns"], color=ps.INK_SECONDARY, lw=0.8)
    ax.text(inst["ns"] - 0.0003, -0.95, f"instant reheating\n(N* = {inst['n_star']:.1f})",
            ha="right", va="bottom", fontsize=7, color=ps.INK_SECONDARY)
    labels = []
    for k, (key, label, gk) in enumerate(NS_ROWS):
        yy = len(NS_ROWS) - 1 - k
        labels.append((yy, label))
        if key == "ours":
            s = cs["ours"]["ns"]
            tau = cs["ours"]["tau"]["mean"]
        elif key == "d3":
            s = {"mean": cs["d3"]["value"], "std": cs["d3"]["error"]}
            tau = 0.0522
        else:
            ch = cs["chains"][key]
            s, tau = ch["ns"], ch.get("tau", {}).get("mean")
        col = GRADED[gk][1] if gk else ps.INK
        m, sd = s["mean"], s["std"]
        ax.plot([m - 2 * sd, m + 2 * sd], [yy, yy], color=col, lw=1.0, solid_capstyle="round")
        ax.plot([m - sd, m + sd], [yy, yy], color=col, lw=3.2, solid_capstyle="round")
        ax.plot([m], [yy], "o", ms=6, mfc=ps.SURFACE, mec=col, mew=1.6, zorder=5)
        txt = f"{m:.4f} ± {sd:.4f}" + (f"   τ = {tau:.3f}" if tau else "")
        if key == "d3":
            txt = f"{m:.4f} ± {sd:.4f}   τ fixed: 0.0522"
        ax.text(1.01, yy, txt, va="center", ha="left", fontsize=7.5, color=ps.INK_SECONDARY,
                transform=ax.get_yaxis_transform())
    ax.set_xlim(0.930, 0.985)
    ax.set_xticks(np.arange(0.93, 0.986, 0.01))
    ax.set_ylim(-1.2, len(NS_ROWS) + 0.4)
    ax.set_yticks([yy for yy, _ in labels])
    ax.set_yticklabels([lab for _, lab in labels], fontsize=8, color=ps.INK)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("spectral index n_s (LCDM; thick: 68%, thin: 95%)")
    ax.set_title("n_s has moved by ~2σ between data combinations", loc="left")
    ps.caption(fig,
               "Every value from the published chains (burn-in handled per release), except D3 (atlas project D3's "
               "own CAMB fit) and 'ours' (Planck 2018's r-free chain weighted by our BK18 likelihood). Coloured "
               "rows are the three combinations the models are graded against. The ACT runs use Planck's Sroll2 "
               "low-ℓ EE likelihood, which prefers a higher optical depth τ; n_s and τ are correlated.", y=0.0)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    return ps.save(fig, _out("06_ns_datasets"))


# ---------------------------------------------------------------------------
# 07 Our BK18 likelihood
# ---------------------------------------------------------------------------


def fig_07_bk18():
    b, bs = _load("bk18")
    mc, mcs = _load("bk18_mcmc")
    v, vs = _load("bk18_validation")
    c, cs = _load("constraints")
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.3), gridspec_kw={"width_ratios": [1.15, 1.15, 1]})
    ax = axes[0]
    r = mc["chain"][:, 0]
    h, e = np.histogram(r, bins=np.linspace(0, 0.1, 51))
    ax.bar(0.5 * (e[1:] + e[:-1]), h / h.max() * 0.98 * b["r_like"].max(), width=np.diff(e) * 0.9,
           color=ps.NEUTRAL_LIGHT, lw=0, label="our MCMC samples")
    ax.fill_between(b["r_grid"], b["r_like"] - b["r_like_boot_std"], b["r_like"] + b["r_like_boot_std"],
                    color=ps.lighten(ps.BLUE, 0.6), lw=0)
    ax.plot(b["r_grid"], b["r_like"], color=ps.BLUE, lw=1.8, label="ours (Rao-Blackwell)")
    ax.plot(b["r_pub"], b["l_pub"], color=ps.INK, lw=1.1, ls="--", label="BK18's own (published)")
    ax.axvline(bs["r95"], color=ps.BLUE, lw=0.8)
    ax.axvline(0.036, color=ps.INK, lw=0.8, ls="--")
    ax.text(bs["r95"] + 0.002, 0.9, f"95%: r < {bs['r95']:.4f}\n(BK18: < 0.036)", fontsize=7.5,
            color=ps.INK_SECONDARY, va="top")
    ax.set_xlim(0, 0.08)
    ax.set_ylim(0, 1.08)
    ax.set_xlabel("tensor-to-scalar ratio r (k = 0.05 Mpc⁻¹)")
    ax.set_ylabel("likelihood (peak = 1)")
    ax.legend(loc="upper right", fontsize=7.2)
    ax.set_title("BK18 alone: our likelihood vs theirs", loc="left")
    ax = axes[1]
    g = c["r1d_grid"]
    for key, col, ls, lab in (("planck18", ps.NEUTRAL, "-", "Planck 2018 alone"),
                              ("planck18_bk18_ours", ps.BLUE, "-", "Planck 2018 + BK18 (ours)"),
                              ("planck_bk18_act", ps.INK, "--", "Planck + BK18 (ACT's joint run)"),
                              ("pactlb_bk18", ps.ORANGE, "-", "P-ACT-LB + BK18 (ACT)")):
        ax.plot(g, c[f"r1d__{key}"], color=col, lw=1.6 if ls == "-" else 1.1, ls=ls, label=lab)
    ax.set_xlim(0, 0.2)
    ax.set_ylim(0, 1.08)
    ax.set_xlabel("r")
    ax.set_title("What B-modes add to the CMB temperature", loc="left")
    ax.legend(loc="upper right", fontsize=7.2)
    ax.text(0.98, 0.42,
            f"95% limits\nPlanck alone: {cs['chains']['planck18_r']['r95']:.3f}\n"
            f"+ BK18 (ours): {cs['ours']['r95']:.4f}\n"
            f"ACT's joint run: {cs['chains']['act_planck-bk_r_camb']['r95']:.4f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.2, color=ps.INK_SECONDARY)
    ax = axes[2]
    d = v["ours"] - v["cobaya"]
    ax.axhspan(-0.01, 0.01, color=ps.lighten(ps.BLUE, 0.85), lw=0)
    ax.axhline(0, color=ps.AXIS, lw=0.8)
    ax.plot(v["cobaya"], d, "o", ms=4, mfc=ps.SURFACE, mec=ps.BLUE, mew=1.0)
    ax.set_ylim(-0.06, 0.06)
    ax.set_xlabel("Cobaya's BK18 χ² for the sample")
    ax.set_ylabel("ours − Cobaya's")
    ax.set_title("Sample by sample: our χ² vs Cobaya's", loc="left")
    ax.text(0.03, 0.95, f"{vs['n']} posterior samples of ACT's Planck + BK18 run,\n"
                        f"each at its own cosmology (CAMB)\nrms difference {vs['rms_diff']:.3f}, "
                        f"max {vs['max_abs_diff']:.3f}\n(χ² itself ≈ 540 for 594 bandpowers)",
            transform=ax.transAxes, va="top", fontsize=7.2, color=ps.INK_SECONDARY)
    fig.tight_layout()
    return ps.save(fig, _out("07_bk18_likelihood"))


# ---------------------------------------------------------------------------
# 08 B-modes: where the gravitational waves hide
# ---------------------------------------------------------------------------


def fig_08_bmodes():
    a, s = _load("bmodes")
    fig, ax = plt.subplots(figsize=(10.4, 6.0))
    ell = a["ell"]
    sel = (ell >= 2) & (ell <= 1500)
    ax.plot(ell[sel], a["lens_bb"][sel], color=ps.INK, lw=1.8, label="gravitational lensing of E-modes (LCDM)")
    for r, lw in zip(s["r_show"], (1.6, 1.2, 1.0, 0.8)):
        t = a[f"tens_bb_r{r}"]
        lt = np.arange(t.size)
        st = lt >= 2
        ax.plot(lt[st], t[st], color=ps.INK, lw=lw, ls="--")
        # labels just under each curve's peak: above it, r = 0.036's label sat on the
        # lensing curve, which crosses the tensor curves near l = 85
        i = 90
        ax.text(lt[i] * 1.03, t[i] / 1.25, f"r = {r:g}", fontsize=7.5, color=ps.INK, va="top")
    for name, col, lab in (("dust150", ps.ORANGE, "polarized dust at 150 GHz (our BK18 fit, 68%)"),
                           ("dust95", ps.lighten(ps.ORANGE, 0.4), "... at 95 GHz (BICEP3)")):
        q = a[f"{name}_q"]
        ax.fill_between(a[f"{name}_l"], q[1], q[3], color=col, alpha=0.45, lw=0, label=lab)
    # BK18's measured 150 GHz auto-spectrum: CMB and dust together (noise-debiased,
    # corrected for filtering); points whose 1-sigma error reaches zero are drawn without
    # a lower bar on the log axis
    lb, v150, e150 = a["bk18_bp_ell"], a["bk18_bp_150"], a["bk18_bp_150_err"]
    lo_err = np.where(v150 - e150 > 0, e150, v150 * 0.999)
    ax.errorbar(lb * 1.02, v150, yerr=[lo_err, e150], fmt="s", ms=4.5, color=ps.NEUTRAL,
                mfc=ps.SURFACE, mec=ps.NEUTRAL, mew=1.1, elinewidth=1.0, capsize=0, zorder=5,
                label="BK18 measured at 150 GHz (CMB + dust, ±1σ)")
    comp = a["bk18_comp"]
    lc = comp[:, 1]
    cmb, lo68, hi68, lo95, hi95 = comp[:, 3], comp[:, 4], comp[:, 5], comp[:, 6], comp[:, 7]
    ax.errorbar(lc, cmb, yerr=[cmb - lo95, hi95 - cmb], fmt="none", ecolor=ps.BLUE, elinewidth=0.8,
                capsize=0, zorder=6)
    ax.errorbar(lc, cmb, yerr=[cmb - lo68, hi68 - cmb], fmt="o", ms=5, color=ps.BLUE, mfc=ps.SURFACE,
                mec=ps.BLUE, mew=1.3, elinewidth=2.2, capsize=0, zorder=7,
                label="BK18's CMB B-mode power, foregrounds removed (68/95%)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(20, 1500)
    ax.set_ylim(2e-5, 0.3)
    ax.set_xlabel("multipole ℓ  (≈ 180° / angular scale)")
    ax.set_ylabel("B-mode power ℓ(ℓ+1)C_ℓ / 2π  [µK²]")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_title("Where primordial gravitational waves would hide: the B-mode spectrum", loc="left")
    ax.text(0.47, 0.97, f"At ℓ = 80, dust at 150 GHz is {s['ratio_dust150_to_r0.036_at_80']:.1f}× "
                        f"brighter than r = 0.036,\nand lensing {s['ratio_lens_to_r0.036_at_80']:.2f}× "
                        "its strength. Tensors peak near ℓ ≈ 90 (2°).",
            transform=ax.transAxes, ha="left", va="top", fontsize=7.5, color=ps.INK_SECONDARY)
    ps.caption(fig, "Theory curves: CAMB with Planck 2018 parameters (tensors with n_t = −r/8). Dust: our MCMC "
                    "of the BK18 likelihood, A_d (ℓ/80)^α_d scaled from 353 GHz through each map's bandpass. "
                    "Points: BK18's released per-bin decomposition (BK-XIII Fig. 16).", y=0.0)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    return ps.save(fig, _out("08_bmodes"))


# ---------------------------------------------------------------------------
# 09 What r means: the energy scale and the field excursion
# ---------------------------------------------------------------------------


def fig_09_energy_scale():
    mod, ms = _load("models")
    c, cs = _load("constraints")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.4, 4.8))
    a_s = math.exp(cfg.LN10_10_AS) * 1e-10
    rr = np.geomspace(1e-5, 0.5, 200)
    ax1.plot(rr, (1.5 * math.pi ** 2 * a_s * rr) ** 0.25 * cfg.M_PL_GEV, color=ps.INK_SECONDARY,
             lw=1.0, label="slow roll: V^¼ = (3π² A_s r / 2)^¼ M_Pl")
    for name in cfg.MODELS:
        r = mod[f"{name}__r"]
        v = mod[f"{name}__V_star_quarter_gev"]
        d = mod[f"{name}__delta_phi"]
        mk = "o" if "param" not in cfg.MODELS[name] else "."
        ax1.plot(r.ravel(), v.ravel(), mk, ms=4 if mk == "o" else 3, color=ps.INK, alpha=0.9)
        ax2.plot(r.ravel(), d.ravel(), mk, ms=4 if mk == "o" else 3, color=ps.INK, alpha=0.9)
    for axx in (ax1, ax2):
        axx.set_xscale("log")
        axx.axvline(cs["ours"]["r95"], color=ps.BLUE, lw=1.2)
        axx.axvspan(cs["ours"]["r95"], 0.5, color=ps.lighten(ps.BLUE, 0.85), lw=0)
        axx.axvline(0.003, color=ps.INK_SECONDARY, lw=0.8, ls="--")
        axx.axvline(0.001, color=ps.INK_SECONDARY, lw=0.8, ls=":")
        axx.set_xlim(1e-5, 0.5)
        axx.set_xlabel("tensor-to-scalar ratio r")
    ax1.set_yscale("log")
    ax1.set_ylabel("energy scale of inflation V*^¼  [GeV]")
    ax1.text(0.13, 3.0e15, "excluded at 95%\nby Planck 2018\n+ BK18 (ours)", fontsize=7.5,
             color=ps.INK_SECONDARY, ha="center")
    for axx in (ax1, ax2):
        axx.text(0.0031, 0.97, "Simons Observatory goal: σ(r) ≈ 0.003", rotation=90, fontsize=7,
                 color=ps.INK_SECONDARY, transform=axx.get_xaxis_transform(), va="top")
        axx.text(0.00103, 0.97, "LiteBIRD goal: δr < 0.001", rotation=90, fontsize=7,
                 color=ps.INK_SECONDARY, transform=axx.get_xaxis_transform(), va="top")
    ax1.axhline(1e16, color=ps.AXIS, lw=0.6)
    ax1.legend(loc="lower right", fontsize=7.4)
    ax1.set_title("r measures the energy of inflation directly", loc="left")
    n = np.array([50.0, 60.0])
    for nn, ls in zip(n, ("-", "--")):
        ax2.plot(rr, nn * np.sqrt(rr / 8), color=ps.INK_SECONDARY, lw=0.9, ls=ls,
                 label=f"Δφ = {nn:.0f} √(r/8): ε frozen at its N* = {nn:.0f} value")
    ax2.plot(rr, 4.6 * np.sqrt(rr / 8), color=ps.INK_SECONDARY, lw=0.9, ls=":",
             label="the CMB's own range only (~4.6 e-folds; Lyth 1997)")
    ax2.axhline(1.0, color=ps.AXIS, lw=0.8)
    ax2.text(1.2e-5, 1.08, "Δφ = M_Pl", fontsize=7, color=ps.INK_SECONDARY)
    ax2.set_yscale("log")
    ax2.set_ylabel("field travelled from N* to the end, Δφ  [M_Pl]")
    ax2.legend(loc="lower right", fontsize=7.2)
    ax2.set_title("... and how far the field had to roll (the Lyth bound)", loc="left")
    ps.caption(fig, "Dots: our exact solutions of every model for N* = 50–60 (large dots: single-shape models; "
                    "small: the families). The energy scale uses each model's own P_R normalised to Planck's "
                    "A_s. Blue: r excluded at 95% by Planck 2018 + BK18 (ours). If ε grows toward the end of "
                    "inflation, Δφ ≥ N*√(r/8): detecting r ≳ 0.002–0.003 would mean a field excursion beyond "
                    "M_Pl and inflation at ~10¹⁶ GeV.", y=0.0)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    return ps.save(fig, _out("09_energy_scale"))


# ---------------------------------------------------------------------------
# 10 Starobinsky against ACT: what reheating would it need?
# ---------------------------------------------------------------------------


def fig_10_reheating():
    mod, ms = _load("models")
    c, cs = _load("constraints")
    fig, ax = plt.subplots(figsize=(10.4, 5.6))
    ns_combo = {"planck18_bk18_ours": cs["ours"]["ns"],
                "pactlb_bk18": cs["chains"]["act_p-actlite-bk-l-b_r_camb"]["ns"],
                "cmbspa_desi_x_bk18": cs["chains"]["cmbspa_desi"]["ns"]}
    n = mod["staro_scan_n"]
    nsv = mod["staro_scan_ns"]
    for key, (label, col) in GRADED.items():
        s = ns_combo[key]
        ax.axhspan(s["mean"] - s["std"], s["mean"] + s["std"], color=ps.RAMPS[col][0], alpha=0.75, lw=0)
        ax.axhline(s["mean"], color=col, lw=1.0)
        need = float(np.interp(s["mean"], nsv, n))
        x_lab = 41.0 if key != "planck18_bk18_ours" else 62.5
        ax.text(x_lab, s["mean"] + 0.00015, f"{label}: n_s = {s['mean']:.4f} ± {s['std']:.4f}"
                + (f" → needs N* ≈ {need:.0f}" if need < n.max() else ""),
                fontsize=7.5, color=ps.INK, va="bottom")
    ax.plot(n, nsv, color=ps.INK, lw=2.0, label="Starobinsky R², exact n_s(N*)")
    rng = ms["starobinsky_nstar_range"]
    y0 = 0.9505
    for k, (lab, key) in enumerate((("matter-like reheating (w = 0), T_reh from instant down to 4 MeV", "w=0.00"),
                                    ("radiation-like (w = 1/3): any T_reh", "w=0.33"),
                                    ("stiff (w = 2/3)", "w=0.67"),
                                    ("kination (w = 1), down to 4 MeV", "w=1.00"))):
        lo, hi = rng[key]
        yy = y0 + 0.0011 * k
        ax.plot([lo, max(hi, lo + 0.25)], [yy, yy], color=ps.INK_SECONDARY, lw=3.0, solid_capstyle="butt")
        ax.text(max(hi, lo + 0.25) + 0.5, yy, f"{lab}: N* = {lo:.1f}–{hi:.1f}" if hi - lo > 0.1
                else f"{lab}: N* = {lo:.1f}", va="center", fontsize=7, color=ps.INK_SECONDARY)
    inst = ms["models"]["starobinsky"]["reheating"]["instant"]
    ax.plot([inst["n_star"]], [inst["ns"]], "o", ms=7, mfc=ps.SURFACE, mec=ps.INK, mew=1.4, zorder=6)
    ps.annotate(ax, f"instant reheating: N* = {inst['n_star']:.1f}, n_s = {inst['ns']:.4f}",
                (inst["n_star"], inst["ns"]), (57.5, 0.9585), fontsize=7.5)
    ax.axvspan(50, 60, color=ps.NEUTRAL_LIGHT, alpha=0.35, lw=0)
    ax.text(55, 0.9805, "the usual band\nN* = 50–60", ha="center", fontsize=7.5, color=ps.INK_SECONDARY)
    ax.set_xlim(40, 84)
    ax.set_ylim(0.949, 0.982)
    ax.set_xlabel("N*: e-folds between the pivot scale leaving the horizon and the end of inflation")
    ax.set_ylabel("spectral index n_s")
    ax.set_title("Starobinsky inflation vs the 2025 data: what would reheating have to do?", loc="left")
    ax.legend(loc="lower right", fontsize=7.5)
    ps.caption(fig, "Bands: n_s ± 1σ of the three graded combinations (r free where available). Bars: N* allowed by "
                    "reheating histories with equation of state w, from instant reheating down to T_reh = 4 MeV "
                    "(nucleosynthesis), each solved self-consistently with the model's own energy scale. Only w > 1/3 "
                    "(a 'stiff' phase) raises N* above the instant-reheating value.", y=0.0)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    return ps.save(fig, _out("10_reheating"))


# ---------------------------------------------------------------------------
# 11 The verdicts
# ---------------------------------------------------------------------------

STATUS = (("allowed", 2.0, "#0ca30c", "✓"), ("under pressure", 3.0, "#fab219", "!"),
          ("ruled out", 99.0, "#d03b3b", "✗"))


def fig_11_verdicts():
    c, cs = _load("constraints")
    ver = cs["verdicts"]
    names = list(cfg.MODELS)
    cols = list(GRADED)
    cw, x0 = 1.45, 2.0          # cell width, where the cells start
    fig, ax = plt.subplots(figsize=(13.4, 6.6))
    ax.set_xlim(0, x0 + cw * len(cols) + 0.05)
    ax.set_ylim(-0.6, len(names) + 0.55)
    ax.axis("off")
    for j, key in enumerate(cols):
        label, col = GRADED[key]
        xc = x0 + cw * (j + 0.5)
        ax.text(xc, len(names) - 0.25, label.replace(": ", ":\n").replace(" (our", "\n(our")
                .replace(", with", ",\nwith"), ha="center", va="bottom", fontsize=8, color=ps.INK)
        ax.plot([xc - 0.45 * cw, xc + 0.45 * cw], [len(names) - 0.4] * 2, color=col, lw=3)

    def param_txt(name, v):
        if v["param"] is None:
            return ""
        p = cfg.MODELS[name]["param"].replace("alpha", "α").replace("mu", "μ")
        return f", {p} = {v['param']:.3g}"

    for i, name in enumerate(names):
        y = len(names) - 1 - i
        ax.text(0.02, y, cfg.MODELS[name]["label"], va="center", fontsize=8.6, color=ps.INK)
        for j, key in enumerate(cols):
            v = ver[name][key]["best_in_band"]
            w = ver[name][key]["instant_reheating"]
            sig = v["sigma"]
            lab, _, colr, icon = next(st for st in STATUS if sig < st[1])
            xl = x0 + cw * j + 0.04 * cw
            ax.add_patch(plt.Rectangle((x0 + cw * j + 0.02 * cw, y - 0.42), 0.96 * cw, 0.84,
                                       fc=ps.lighten(colr, 0.82), ec="none"))
            sig_txt = f"> {sig:.1f}σ" if v["beyond_resolution"] else f"{sig:.1f}σ"
            wsig = f"> {w['sigma']:.1f}σ" if w["beyond_resolution"] else f"{w['sigma']:.1f}σ"
            ax.text(xl, y + 0.2, f"{icon} {lab}", va="center", fontsize=8, color=ps.INK,
                    fontweight="bold")
            ax.text(xl, y - 0.04, f"best in band: {sig_txt} (N* = {v['n_star']:.0f}"
                    f"{param_txt(name, v)})", va="center", fontsize=6.9, color=ps.INK_SECONDARY)
            ax.text(xl, y - 0.25, f"instant reheating: {wsig} (N* = {w['n_star']:.1f})",
                    va="center", fontsize=6.9, color=ps.INK_SECONDARY)
    ax.set_title("The verdicts: each model at its best within N* = 50–60 and its parameter range",
                 loc="left", pad=14)
    ps.caption(fig, "σ-equivalent of the smallest 2D credible region in (n_s, r) that contains the model's point: "
                    "1σ ↔ 68.3%, 2σ ↔ 95.4%. Allowed < 2σ, under pressure 2–3σ, ruled out > 3σ. 'Instant "
                    "reheating' uses the N* that the model's own energy scale implies if the universe "
                    "thermalised at once. Beyond 99.9% (3.3σ) the chains run out of samples ('>').", y=0.03)
    return ps.save(fig, _out("11_verdicts"))


FIGURES = [fig_01_cmb_horizons, fig_02_hubble_sphere, fig_03_potentials, fig_04_mukhanov_sasaki,
           fig_05_ns_r_plane, fig_06_ns_datasets, fig_07_bk18, fig_08_bmodes, fig_09_energy_scale,
           fig_10_reheating, fig_11_verdicts]


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    ps.use_style()
    for fn in FIGURES:
        tag = fn.__name__.split("_")[1]
        if argv and not any(tag.startswith(a) for a in argv):
            continue
        out = fn()
        print("  " + ", ".join(str(o).rsplit("/", 1)[-1] for o in out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
