"""The project's figures, made from cached data products and data/results.json."""

from __future__ import annotations

import json
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

import analysis as A  # noqa: E402
import fitting as F  # noqa: E402
import models as M  # noqa: E402
import plotstyle as S  # noqa: E402
import uvdata as U  # noqa: E402
from config import (FIGURE_DIR, M87_2017_BANDS, M87_2017_DAYS, PAPER6_XSRING, PAPER6_XSRINGAUSS,  # noqa: E402
                    PUBLISHED, RAD_PER_UAS, RESULTS_PATH, SHADOW_FACTOR_SCHWARZSCHILD, STATION_NAMES)

warnings.filterwarnings("ignore")
S.use_style()
GL = 1e9  # one billion wavelengths


def _res() -> dict:
    return json.loads(RESULTS_PATH.read_text())


def orientation(u, v) -> np.ndarray:
    """'ew', 'ns' or 'diag' from the baseline's position angle in the (u, v) plane."""
    ang = np.degrees(np.arctan2(u, v)) % 180.0
    out = np.full(len(u), "diag", dtype=object)
    out[(ang > 60) & (ang < 120)] = "ew"
    out[(ang < 30) | (ang > 150)] = "ns"
    return out


ORIENT_COLOR = {"ew": S.C_EW, "ns": S.C_NS, "diag": S.C_DIAG}
ORIENT_LABEL = {"ew": "east-west baselines", "ns": "north-south", "diag": "diagonal"}


def _prep(day, band):
    return U.prepare(A.m87_2017_path(day, band), f"{M87_2017_DAYS[day]} {band}")


def _best_model(obs_key, model_key=A.FIDUCIAL):
    f = F.load_fit(A.fit_path(obs_key, model_key))
    if f is None:
        return None, None
    model = A.MODELS[model_key]()
    best = np.array([f["summary"]["best"][n] for n in model.names])
    return model, best


def _bl_label(a, b):
    return f"{STATION_NAMES.get(a, a)}-{STATION_NAMES.get(b, b)}"


# ---------------------------------------------------------------------------
# Figure 1: what an interferometer measures
# ---------------------------------------------------------------------------
def fig_uv_and_amplitudes():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.6), gridspec_kw={"width_ratios": [1, 1.25]})
    prep = {b: _prep("101", b) for b in M87_2017_BANDS}
    model, best = _best_model("m87_2017_101_lo")
    # background: |V| of the fitted April 11 model over the (u, v) plane
    if model is not None:
        g = np.linspace(-9.5, 9.5, 400) * GL
        UU, VV = np.meshgrid(g, g)
        amp = np.abs(model.compact_vis(best, UU.ravel(), VV.ravel())).reshape(UU.shape)
        amp /= np.abs(model.compact_vis(best, np.array([1e3]), np.array([0.0])))[0]
        ax1.imshow(np.log10(np.clip(amp, 1e-3, 1)), extent=[9.5, -9.5, -9.5, 9.5], origin="lower",
                   cmap="Greys", vmin=-3, vmax=0.3, alpha=0.55, interpolation="bilinear")
    for b in M87_2017_BANDS:
        d = prep[b].raw
        o = orientation(d.u, d.v)
        intra = d.intra_site()
        for key in ("ew", "ns", "diag"):
            m = (o == key) & ~intra
            for sgn in (1, -1):
                ax1.plot(sgn * d.u[m] / GL, sgn * d.v[m] / GL, ".", ms=1.2, color=ORIENT_COLOR[key],
                         label=ORIENT_LABEL[key] if (sgn == 1 and b == "lo") else None)
    ax1.set_xlim(9.5, -9.5)
    ax1.set_ylim(-9.5, 9.5)
    ax1.set_aspect("equal")
    ax1.grid(False)
    ax1.set_xlabel("u (billion wavelengths, east)")
    ax1.set_ylabel("v (billion wavelengths, north)")
    ax1.set_title("1a  Where the EHT sampled M87 on 11 April 2017")
    ax1.legend(loc="lower left", markerscale=6, fontsize=7.5)
    ax1.text(0.02, 0.98, "grey: |V| of our fitted ring model\n(each track samples it along a line)",
             transform=ax1.transAxes, va="top", fontsize=7.5, color=S.INK_SECONDARY)
    for circ_gl in (3.67, 8.4):
        ax1.add_patch(Circle((0, 0), circ_gl, fill=False, ls=":", lw=0.7, color=S.MUTED))

    # amplitudes vs baseline length, all four nights
    for day, col in S.DAY_COLORS.items():
        for b in M87_2017_BANDS:
            d = _prep(day, b).avg
            a = U.debiased_amplitude(np.abs(d.vis), d.sigma)
            ok = ~d.intra_site()
            ax2.errorbar(d.uvdist[ok] / GL, a[ok], yerr=d.sigma[ok], fmt="o", ms=2.3, lw=0.5, color=col,
                         alpha=0.75, label=M87_2017_DAYS[day] if b == "lo" else None)
    rho = np.linspace(0.01, 9.5, 600) * GL
    F0 = 0.6
    ring = np.abs(M.SymmetricShape("ring", large=False).vis(np.array([43.0]), rho, 0 * rho)) * F0
    disk = np.abs(M.SymmetricShape("disk", large=False).vis(np.array([43.0]), rho, 0 * rho)) * F0
    gau = np.abs(M.SymmetricShape("gaussian", large=False).vis(np.array([43.0 / 1.7]), rho, 0 * rho)) * F0
    ax2.plot(rho / GL, ring, color=S.C_OURS, lw=1.3, label="thin ring, 43 µas")
    ax2.plot(rho / GL, disk, color=S.MUTED, lw=1.1, ls="--", label="filled disk, 43 µas")
    ax2.plot(rho / GL, gau, color=S.MUTED, lw=1.1, ls=":", label="Gaussian blob")
    ax2.axvline(3.4, color=S.C_GR, lw=0.8, alpha=0.8)
    ax2.text(3.3, 0.0042, "first null ≈ 3.4 Gλ:\na thin ring 0.766/3.4 Gλ\n= 46 µas across", fontsize=7.5,
             color=S.C_GR, va="bottom", ha="right")
    ax2.annotate("second minimum ≈ 8.3 Gλ", xy=(8.25, 0.01), xytext=(4.6, 0.0042), fontsize=7.5, color=S.INK_SECONDARY,
                 arrowprops={"arrowstyle": "-", "color": S.INK_SECONDARY, "lw": 0.7})
    ax2.set_yscale("log")
    ax2.set_ylim(3e-3, 2.5)
    ax2.set_xlim(0, 9.5)
    ax2.set_xlabel("baseline length (billion wavelengths)")
    ax2.set_ylabel("visibility amplitude (Jy)")
    ax2.set_title("1b  The ring's signature: a deep null, then a second peak")
    ax2.legend(loc="upper right", fontsize=7.2, ncol=3)
    S.caption(fig, "EHT 2017 M87 data release 2019-D01-01 (HOPS pipeline, both bands), scan-averaged; intra-site "
                   "baselines omitted on the right. Error bars: thermal noise plus 1% systematic. Model curves are "
                   "shapes with the compact flux of 0.6 Jy, not fits.")
    return S.save(fig, FIGURE_DIR / "01_uv_and_amplitudes")


# ---------------------------------------------------------------------------
# Figure 2: closure phases
# ---------------------------------------------------------------------------
TRIANGLES = [("AA", "LM", "SM"), ("AA", "AZ", "LM"), ("AA", "AZ", "SM"), ("AZ", "LM", "SM")]


def _closure_series(prep, tri):
    """All closure phases of one triangle in a prepared data set (not only the minimal set)."""
    d = U.canonical_orientation(prep.avg)
    a, b, c = sorted(tri)
    idx = {(s, x, y): i for i, (s, x, y) in enumerate(zip(d.scan, d.t1, d.t2))}
    t, psi, sig = [], [], []
    for s in np.unique(d.scan):
        if (s, a, b) in idx and (s, b, c) in idx and (s, a, c) in idx:
            i, j, k = idx[(s, a, b)], idx[(s, b, c)], idx[(s, a, c)]
            t.append(d.time[i])
            psi.append(np.angle(d.vis[i] * d.vis[j] * np.conj(d.vis[k])))
            sig.append(np.sqrt(sum((d.sigma[q] / abs(d.vis[q])) ** 2 for q in (i, j, k))))
    return np.array(t), np.degrees(psi), np.degrees(sig), (i, j, k) if t else None


def gmst_hours(mjd: int, ut_hours):
    """Greenwich mean sidereal time: equal GMST = equal (u, v) points on different nights."""
    from astropy.time import Time

    return Time(mjd + np.asarray(ut_hours) / 24.0, format="mjd").sidereal_time("mean", "greenwich").hour


def fig_closure_phases():
    fig, axs = plt.subplots(2, 2, figsize=(12, 7.2), sharex=True)
    for ax, tri in zip(axs.ravel(), TRIANGLES):
        for day, col in S.DAY_COLORS.items():
            prep = _prep(day, "lo")
            t, psi, sig, _ = _closure_series(prep, tri)
            if len(t) == 0:
                continue
            ax.errorbar(gmst_hours(prep.raw.mjd, t), psi, yerr=sig, fmt="o", ms=3, lw=0.7, color=col,
                        label=M87_2017_DAYS[day])
            model, best = _best_model(f"m87_2017_{day}_lo")
            if model is not None:
                tt = np.linspace(t.min() - 0.1, t.max() + 0.1, 300)
                # interpolate (u, v) of each baseline along the night
                d = U.canonical_orientation(prep.avg)
                a, b, c = sorted(tri)
                uv = {}
                for x, y in ((a, b), (b, c), (a, c)):
                    m = (d.t1 == x) & (d.t2 == y)
                    uv[(x, y)] = (np.interp(tt, d.time[m], d.u[m]), np.interp(tt, d.time[m], d.v[m]))
                ph = {k: np.angle(model.vis(best, *v)) for k, v in uv.items()}
                psi_m = np.degrees(np.angle(np.exp(1j * (ph[(a, b)] + ph[(b, c)] - ph[(a, c)]))))
                jump = np.abs(np.diff(psi_m)) > 180  # don't draw the wrap from +180 to -180
                psi_m[1:][jump] = np.nan
                ax.plot(gmst_hours(prep.raw.mjd, tt), psi_m, color=col, lw=0.9, alpha=0.8)
        ax.axhline(0, color=S.MUTED, lw=0.6)
        ax.set_ylim(-180, 180)
        ax.set_yticks([-180, -90, 0, 90, 180])
        ax.set_title(" - ".join(STATION_NAMES[s] for s in tri), fontsize=10)
        ax.set_ylabel("closure phase (deg)")
    for ax in axs[1]:
        ax.set_xlabel("Greenwich sidereal time (hours): equal times sample the same Fourier points")
    axs[0, 0].legend(loc="lower left", fontsize=7.5, ncol=2)
    fig.suptitle("2  Closure phases: nonzero, so M87 is lopsided; changing, so it evolved over six days", x=0.01,
                 ha="left")
    S.caption(fig, "Low band. A point-symmetric source (ring, disk, Gaussian) would give only 0 or ±180 deg. Lines: our "
                   "fitted ring model for each night, sampled along that night's baseline tracks. Station phase errors "
                   "cancel exactly in these sums.")
    return S.save(fig, FIGURE_DIR / "02_closure_phases")




# ---------------------------------------------------------------------------
# Figure 3: the ring measured night by night, against Paper VI
# ---------------------------------------------------------------------------
def fig_ring_fits():
    res = _res()
    per = res["m87_2017_fits"]
    keys = list(per)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.6))
    x = np.arange(len(keys))
    labels = [f"{M87_2017_DAYS[k.split('_')[2]].replace('April ', 'Apr ')}\n{k.split('_')[3].upper()}" for k in keys]
    for i, k in enumerate(keys):
        e = per[k]
        p6, p6g = e["paper6_xsring"], e["paper6_xsringauss"]
        ax1.errorbar(i - 0.24, p6["d"][0], yerr=[[p6["d"][2]], [p6["d"][1]]], fmt="s", ms=4, color=S.C_PUBLISHED,
                     mfc="white", label="Paper VI crescent (closures)" if i == 0 else None)
        ax1.errorbar(i - 0.08, p6g["d"][0], yerr=[[p6g["d"][2]], [p6g["d"][1]]], fmt="D", ms=3.5, color=S.C_PUBLISHED,
                     label="Paper VI crescent (amplitudes)" if i == 0 else None)
        if e.get("xsring"):
            q = e["xsring"]["d"]
            ax1.errorbar(i + 0.08, q["median"], yerr=[[q["median"] - q["lo68"]], [q["hi68"] - q["median"]]], fmt="s",
                         ms=4, color=S.C_OURS, label="ours, same crescent model" if i == 0 else None)
        if e.get("fiducial"):
            q = e["fiducial"]["d"]
            ax1.errorbar(i + 0.24, q["median"], yerr=[[q["median"] - q["lo68"]], [q["hi68"] - q["median"]]], fmt="o",
                         ms=4.5, color=S.C_GR, label="ours, blurred ring (δ-ring diameter)" if i == 0 else None)
        if e.get("amplitudes"):
            q = e["amplitudes"]["d"]
            ax1.errorbar(i + 0.38, q["median"], yerr=[[q["median"] - q["lo68"]], [q["hi68"] - q["median"]]], fmt="^",
                         ms=4.5, color=S.C_GR, mfc="white", label="ours, from the amplitude nulls alone" if i == 0 else None)
        for ax, key in ((ax2, "pa_deg"),):
            if e.get("fiducial"):
                q = e["fiducial"]["pa_deg"]
                ax.errorbar(i + 0.15, q["median"], yerr=[[q["median"] - q["lo68"]], [q["hi68"] - q["median"]]],
                            fmt="o", ms=4.5, color=S.C_GR, label="ours, brightest side of the ring" if i == 0 else None)
            if e.get("xsring"):
                q = e["xsring"]["pa_deg"]
                ax.errorbar(i, q["median"], yerr=[[q["median"] - q["lo68"]], [q["hi68"] - q["median"]]], fmt="s",
                            ms=4, color=S.C_OURS, label="ours, crescent (where its light is)" if i == 0 else None)
            ax.errorbar(i - 0.15, p6["pa"][0], yerr=[[p6["pa"][2]], [p6["pa"][1]]], fmt="s", ms=4, mfc="white",
                        color=S.C_PUBLISHED, label="Paper VI crescent" if i == 0 else None)
    lo, hi = PUBLISHED["m87_2017_2021_mean_diameter_uas"][0] - PUBLISHED["m87_2017_2021_mean_diameter_uas"][1], \
        PUBLISHED["m87_2017_2021_mean_diameter_uas"][0] + PUBLISHED["m87_2017_2021_mean_diameter_uas"][1]
    ax1.axhspan(42.0, 46.4, color=S.C_GR, alpha=0.08, lw=0)
    ax1.text(len(keys) - 0.6, 46.2, "EHT 2025: δ-ring diameter, 2017 (all methods)", fontsize=7, color=S.C_GR,
             ha="right", va="top")
    ax1.set_xticks(x, labels, fontsize=8)
    ax1.set_ylabel("ring diameter (µas)")
    ax1.set_title("3a  Diameter, each night and band")
    ax1.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, columnspacing=1.0)
    ax2.set_xticks(x, labels, fontsize=8)
    ax2.set_ylabel("position angle (deg east of north)")
    ax2.set_title("3b  The bright side: south-east, then south-west")
    ax2.axhline(180, color=S.MUTED, lw=0.6)
    ax2.text(-0.4, 181, "due south", fontsize=7, color=S.MUTED)
    ax2.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, columnspacing=1.0)
    fig.subplots_adjust(bottom=0.3)
    S.caption(fig, "Our fits use closure phases and log closure amplitudes only, with the full covariance, like Paper "
                   "VI's dynesty fits. Two ways to define 'the diameter' differ by a few µas: Paper VI's crescent "
                   "averages its inner and outer edges, while a Gaussian-blurred thin ring keeps the null of the thin "
                   "ring.")
    return S.save(fig, FIGURE_DIR / "03_ring_fits")


# ---------------------------------------------------------------------------
# Figure 4: is it really a ring? Model comparison by Bayesian evidence
# ---------------------------------------------------------------------------
def fig_model_comparison():
    res = _res()
    comp = res["model_comparison"]
    order = sorted(comp, key=lambda k: comp[k]["logz"])  # best at the top
    ring_like = {"ring", "crescent", "mg1", "mg2", "mg3", "mg1g2", "mg2g1", "mg2g2", "mg2g3", "mg3g2", "xs2"}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.4), sharey=True, gridspec_kw={"width_ratios": [1.3, 1]})
    y = np.arange(len(order))
    best = max(comp[k]["logz"] for k in order)
    for i, k in enumerate(order):
        v = comp[k]
        dz = v["logz"] - best
        col = S.C_GR if k in ring_like else S.MUTED
        ax1.barh(i, max(-dz, 0.3), color=col, height=0.62)
        ax1.text(max(-dz, 0.3) * 1.15, i, f"{dz:,.0f}" if dz < -0.5 else "best", va="center", fontsize=7.5,
                 color=S.INK_SECONDARY)
        red = v["reduced_chi2"]
        ax2.plot(red.get("cphase", np.nan), i, "o", color=S.C_EW, ms=5.5, label="closure phases" if i == 0 else None)
        ax2.plot(red.get("lcamp", np.nan), i, "s", color=S.C_NS, ms=5, label="log closure amplitudes" if i == 0 else None)
    ax1.set_yticks(y, [f"{comp[k]['name']}  ({comp[k]['n_params']})" for k in order], fontsize=8)
    ax1.set_xscale("log")
    ax1.set_xlim(0.2, None)
    ax1.set_xlabel("how much worse than the best model:  −Δ ln(evidence)   (5 is already decisive)")
    ax1.set_title("4a  Rings win, by hundreds to thousands of log-units")
    ax2.set_xscale("log")
    ax2.axvline(1, color=S.INK, lw=0.8)
    ax2.set_xlabel("reduced χ² of the best fit  (1 = fits down to the noise)")
    ax2.set_title("4b  How well each model fits")
    ax2.legend(fontsize=7.5, loc="lower right")
    S.caption(fig, "M87, 11 April 2017, low band; number of parameters in brackets. Every model is fitted to the same "
                   "106 closure phases and 108 log closure amplitudes by nested sampling. The Bayesian evidence "
                   "rewards fit quality and penalises needless parameters. Aqua: models with a dark centre "
                   "(rings, crescents); grey: models without one.")
    return S.save(fig, FIGURE_DIR / "04_model_comparison")


# ---------------------------------------------------------------------------
# Figure 5: images, ours and the EHT pipeline's
# ---------------------------------------------------------------------------
def _show_image(ax, img, fov=128.0, title=None, feats=None, vmax=None):
    ax.imshow(img, extent=[fov / 2, -fov / 2, -fov / 2, fov / 2], cmap=S.IMAGE_CMAP, vmin=0,
              vmax=vmax if vmax is not None else img.max(), interpolation="bicubic")
    if feats:
        t = np.linspace(0, 2 * np.pi, 200)
        r = feats["d"] / 2
        ax.plot(feats["x0"] + r * np.sin(t), feats["y0"] + r * np.cos(t), ls=":", lw=0.8, color="#7fe3ff")
    ax.set_xticks([40, 0, -40])
    ax.set_yticks([-40, 0, 40])
    ax.tick_params(labelsize=7)
    ax.grid(False)
    if title:
        ax.set_title(title, fontsize=8.5)


def fig_images():
    res = _res()
    fig, axs = plt.subplots(2, 4, figsize=(13, 7.0))
    for j, day in enumerate(M87_2017_DAYS):
        f = res["image_features"].get(day, {})
        ours = np.load(A.IMG_DIR / f"m87_2017_{day}.npy")
        eht = A.ehtim_pipeline_image(day)
        fo, fe = f.get("ours"), f.get("ehtim")
        _show_image(axs[0, j], ours, title=f"{M87_2017_DAYS[day]}: our imager\nring {fo['d']:.1f} µas, "
                                             f"brightest at {fo['pa_deg']:.0f}°", feats=fo)
        _show_image(axs[1, j], eht, title=f"EHT's eht-imaging script\nring {fe['d']:.1f} µas, "
                                            f"brightest at {fe['pa_deg']:.0f}°", feats=fe)
    for ax in axs[:, 0]:
        ax.set_ylabel("µas north")
    for ax in axs[1]:
        ax.set_xlabel("µas east")
    fig.suptitle("5  Images from the same data: our own imager (top) and the EHT's published pipeline (bottom)",
                 x=0.01, ha="left")
    S.caption(fig, "Both bands of each night. Top: our regularised-maximum-likelihood imager (closure phases, log "
                   "closure amplitudes and amplitudes with each station's gain free within its calibration budget; "
                   "maximum entropy towards a 40 µas Gaussian, plus smoothness). Bottom: the EHT's own eht-imaging "
                   "script (2019-D01-02), run unchanged except for using the exact Fourier transform. Dotted circles: "
                   "rings found by our image-domain feature extraction. Field of view 128 µas, east to the left.")
    return S.save(fig, FIGURE_DIR / "05_images")


def fig_validation():
    import validation as V

    out = V.run("101", settings=A.IMAGER_SETTINGS, weights=A.IMAGER_WEIGHTS)
    fig, axs = plt.subplots(2, len(out), figsize=(13, 6.6))
    for j, (k, r) in enumerate(out.items()):
        _show_image(axs[0, j], r["truth"], title=f"SIMULATED truth:\n{r['label']}")
        red = r["reduced_chi2"]
        _show_image(axs[1, j], r["recon"], title=f"reconstruction\nχ²/N: phases {red['cphase']:.1f}, "
                                                    f"amps {red['amp']:.1f}")
    axs[0, 0].set_ylabel("truth (simulation)")
    axs[1, 0].set_ylabel("what our imager returns")
    fig.suptitle("6  Validation on simulated data: does the imager invent rings?  No.", x=0.01, ha="left")
    S.caption(fig, "Simulated sources sampled exactly like the real 11 April 2017 observation: the same baselines, "
                   "times and thermal noise, random station gains as large as the real calibration allows, and 0.6 Jy "
                   "of extra flux on the shortest baselines. Imaged with the same settings as Figure 5. Only the "
                   "ring comes back as a ring. Not real data.")
    return S.save(fig, FIGURE_DIR / "06_imaging_validation")


# ---------------------------------------------------------------------------
# Figure 7: from ring to mass
# ---------------------------------------------------------------------------
def fig_mass():
    res = _res()
    mb = res["m87_mass"]
    D = PUBLISHED["m87_distance_mpc"][0]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.4), gridspec_kw={"width_ratios": [1.45, 1]})
    rows = []  # (label, median, lo, hi, color, kind)
    g = mb["grmhd_calibrated"]["theta_g"]
    rows.append(("this project: ring + GRMHD calibration\n(crescent d = %.1f µas, α = %.2f)" %
                 (mb["grmhd_calibrated"]["diameter_uas"], mb["grmhd_calibrated"]["alpha"]),
                 g["median"], g["lo68"], g["hi68"], S.C_OURS, "ours"))
    s = mb["shadow_schwarzschild"]["theta_g"]
    k = mb["shadow_kerr_min"]["theta_g"]
    rows.append(("this project: if the ring were the shadow's edge\n(δ-ring d = %.1f µas, α = 10.39 to 9.6)" %
                 mb["shadow_schwarzschild"]["diameter_uas"], s["median"], s["lo68"], k["hi68"], S.C_GR, "ours"))
    th, er = PUBLISHED["m87_theta_g_uas"]
    rows.append(("EHT 2019 (Paper VI)", th, th - er, th + er, S.C_PUBLISHED, "pub"))
    for name, p in mb["predictions"].items():
        col = S.C_GAS if "gas" in name else S.C_STARS
        rows.append((f"{name} dynamics" if "dynamics" not in name else name, p["theta_g"], p["theta_g"] - p["theta_lo"],
                     p["theta_g"] + p["theta_up"], col, "dyn"))
    ys = []
    for i, (lab, med, lo, hi, col, kind) in enumerate(rows):
        y = len(rows) - 1 - i
        ys.append(y)
        ax1.errorbar(med, y, xerr=[[med - lo], [hi - med]], fmt="o" if kind != "pub" else "s", ms=6, color=col,
                     mfc=col if kind != "pub" else "white", lw=1.6, capsize=0)
    ax1.set_xlim(0.5, 8.0)
    ax1.set_yticks(ys, [r[0] for r in rows], fontsize=7.6)
    for tl, r in zip(ax1.get_yticklabels(), rows):
        tl.set_color(S.INK if r[5] != "dyn" else r[4])
    ax1.set_xlabel("angular size of one gravitational radius, GM/c²D (µas)")
    sec = ax1.secondary_xaxis("top", functions=(lambda t: np.asarray(A.theta_g_to_mass(t, D)) / 1e9,
                                                 lambda m: np.asarray(A.mass_to_theta_g(np.asarray(m) * 1e9, D))))
    sec.set_xlabel(f"mass at {D:g} Mpc (billion suns)")
    ax1.set_title("7a  Weighing M87*: the ring sides with the stars", pad=28)
    # right: the April 11 image with the ring sizes each mass predicts
    img = A.ehtim_pipeline_image("101")
    _show_image(ax2, img)
    t = np.linspace(0, 2 * np.pi, 300)
    f = res["image_features"]["101"]["ehtim"]
    for name, col, ls in (("gas (Walsh+13)", S.C_GAS, "--"), ("stars (Gebhardt+11)", S.C_STARS, "-")):
        dpred = mb["predictions"][name]["crescent_diameter_uas"]
        ax2.plot(f["x0"] + dpred / 2 * np.sin(t), f["y0"] + dpred / 2 * np.cos(t), ls=ls, color=col, lw=1.3,
                 label=f"ring predicted from {'stellar' if name.startswith('stars') else 'gas'} dynamics: {dpred:.0f} µas")
    ax2.legend(loc="upper left", fontsize=6.8, frameon=True, framealpha=0.85, facecolor="white", edgecolor="none")
    ax2.set_title("7b  Predicted rings on the 11 April image")
    S.caption(fig, "Right: the EHT's own imaging script re-run on the public data. Error bars 68%. Our "
                   "GRMHD-calibrated value uses Paper VI's α for the same crescent model and "
                   "carries its ±11% spread over simulations; the night-to-night scatter of the diameter is included. "
                   "Dynamical values are mass/distance ratios (Paper VI Table 9; Liepold et al. 2023 and Simon et al. "
                   "2024 converted at 16.8 Mpc), so this comparison needs no distance. Predicted rings: θ × α.")
    return S.save(fig, FIGURE_DIR / "07_mass")


# ---------------------------------------------------------------------------
# Figure 8: Sgr A* and H1's prediction
# ---------------------------------------------------------------------------
def fig_sgra():
    res = _res()
    st = res["sgra_test"]
    fig = plt.figure(figsize=(13.5, 5.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.45, 0.8, 1.1], wspace=0.28)
    ax1, ax3, ax2 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
    for key, col, lab in (("sgra_2017_lo", S.BLUE, "227 GHz"), ("sgra_2017_hi", S.ORANGE, "229 GHz")):
        for k, prep in enumerate(A.load_obs_deblurred(key)):
            d = prep.avg
            a = U.debiased_amplitude(np.abs(d.vis), d.sigma)
            ax1.errorbar(d.uvdist / GL, a, yerr=d.sigma, fmt="o", ms=2.0, lw=0.4, color=col, alpha=0.5,
                         label=f"{lab}, April 6 + 7" if k == 0 else None)
    best = st["amplitudes"]["best_by_band"].get("lo")
    if best:
        mkey = "sa" + best["model"].split("_m")[-1]
        model, theta = _best_model("sgra_2017_lo", mkey)
        src, th = model.source, theta[: model.source.ndim]
        rho = np.linspace(0.02, 9, 500) * GL
        for j, ang in enumerate(np.radians([0, 45, 90, 135])):
            V = src.vis(th, rho * np.sin(ang), rho * np.cos(ang))
            ax1.plot(rho / GL, np.abs(V), color=S.INK, lw=0.8, alpha=0.75,
                     label="ring fitted beyond 2 Gλ (4 directions)" if j == 0 else None)
        _show_image(ax3, posterior_mean_image("sgra_2017_lo", mkey), title="fitted ring (no scattering)")
    ax1.axvspan(0, 2.0, color=S.MUTED, alpha=0.08, lw=0)
    ax1.text(1.0, 0.0025, "not fitted:\nshaped by the\nnormalisation and\nLMT calibration", ha="center",
             fontsize=6.8, color=S.INK_SECONDARY)
    ax1.set_yscale("log")
    ax1.set_ylim(2e-3, 1.6)
    ax1.set_xlim(0, 9)
    ax1.set_xlabel("baseline length (billion wavelengths)")
    ax1.set_ylabel("amplitude / total flux, scattering removed")
    ax1.set_title("8a  Sgr A*: a ring's nulls at 2.9 and 6.5 Gλ")
    ax1.legend(fontsize=7.0, loc="upper right")
    # right: the diameters
    pred, perr = st["shadow_prediction_uas"], st["shadow_prediction_err"]
    rp = st["ring_prediction_uas"]
    ax2.axvspan(pred - perr, pred + perr, color=S.C_GR, alpha=0.3, lw=0)
    ax2.axvline(pred, color=S.C_GR, lw=1.2)
    ax2.axvline(rp, color=S.C_GR, lw=1.0, ls="--")
    ax2.text(pred - 0.3, 6.85, f"H1 (star S2):\nshadow {pred:.1f} µas", color=S.C_GR, fontsize=7.2, ha="right",
             va="top")
    ax2.text(rp + 0.3, 6.85, f"ring expected\n(EHT calibration)\n{rp:.1f} µas", color=S.C_GR, fontsize=7.2, va="top")
    y = 5.6
    rows = [("amplitudes", k, v) for k, v in res.get("sgra_amp_fits", {}).items()] + \
           [("closures", k, v) for k, v in res.get("sgra_fits", {}).items()]
    for fam, key, v in rows:
        q = v["d_peak"]
        col = S.BLUE if key.startswith("lo") else S.ORANGE
        filled = fam == "amplitudes"
        ax2.errorbar(q["median"], y, xerr=[[q["median"] - q["lo68"]], [q["hi68"] - q["median"]]],
                     fmt="D" if "casa" in key else "o", color=col, ms=5, mfc=col if filled else "white", capsize=0)
        lab = key.replace("lo", "227").replace("hi", "229").replace("_casa", " CASA").replace("_m", " m=")
        ax2.text(36.3, y, f"{'amps' if filled else 'closures'} {lab}", fontsize=6.4, va="center",
                 color=S.INK_SECONDARY)
        y -= 0.42
    eht = PUBLISHED["sgra_ring_diameter_uas"]
    ax2.errorbar(eht[0], y - 0.15, xerr=eht[1], fmt="s", mfc="white", color=S.C_PUBLISHED, ms=6, capsize=0)
    ax2.text(36.3, y - 0.15, "EHT 2022, all methods", fontsize=6.4, va="center", color=S.C_PUBLISHED)
    ax2.set_ylim(y - 0.7, 7.1)
    ax2.set_yticks([])
    ax2.set_xlim(36, 62)
    ax2.set_xlabel("ring diameter at peak brightness (µas)")
    ax2.set_title("8b  Measured ring vs the orbit's prediction")
    S.caption(fig, "EHT 2022-D02-01, light-curve-normalised, both nights; left: scan-averaged amplitudes divided by "
                   "the Johnson et al. (2018) scattering kernel. Filled: rings fitted to the amplitudes beyond 2 Gλ "
                   "(the nulls), scattering applied to the model, noise budget fitted alongside. Hollow: fits to "
                   "closure quantities alone, which a broad halo makes degenerate for Sgr A*. Diamonds: the "
                   "independent CASA calibration. The prediction uses only H1's mass and distance from the star S2.")
    return S.save(fig, FIGURE_DIR / "08_sgra")


# ---------------------------------------------------------------------------
# Figure 9: the same ring a year later
# ---------------------------------------------------------------------------
def posterior_mean_image(obs_key: str, model_key: str = A.FIDUCIAL, *, n: int = 120, npix: int = 96,
                         fov: float = 128.0, compact_only: bool = True) -> np.ndarray:
    """Average of model images over posterior samples (what the fit 'looks like')."""
    f = F.load_fit(A.fit_path(obs_key, model_key))
    model = A.MODELS[model_key]()
    from imaging import pixel_grid

    X, Y, _ = pixel_grid(npix, fov)
    rng = np.random.default_rng(0)
    samp = f["samples"][rng.choice(len(f["samples"]), size=min(n, len(f["samples"])), replace=False)]
    src = getattr(model, "source", model)
    img = np.zeros_like(X)
    for th in samp:
        th = th[: src.ndim]
        im = src.image(th, X, Y)
        img += im / im.sum()
    return img / len(samp)


def fig_two_years():
    res = _res()
    fig = plt.figure(figsize=(13, 7.4))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.05], hspace=0.38, wspace=0.25)
    ax_d = fig.add_subplot(gs[0, :2])
    ax_p = fig.add_subplot(gs[0, 2:])
    groups = [("2017", res["m87_2017_fits"], S.BLUE), ("2018", res.get("m87_2018_fits", {}), S.ORANGE)]
    x = 0
    ticks, labels, group_pos = [], [], []
    for year, fits, col in groups:
        for k, e in fits.items():
            f = e.get("fiducial")
            if not f:
                continue
            for ax, key in ((ax_d, "d"), (ax_p, "pa_deg")):
                q = f[key]
                ax.errorbar(x, q["median"], yerr=[[q["median"] - q["lo68"]], [q["hi68"] - q["median"]]], fmt="o",
                            color=col, ms=5)
            ticks.append(x)
            lab = e["label"].replace("M87 ", "").replace(year + " ", "").replace("April ", "")
            labels.append(lab.replace(" ", "\n", 1))
            x += 1
        group_pos.append((year, np.mean(ticks[-len([e for e in fits.values() if e.get("fiducial")]):])))
        x += 0.6
    for ax in (ax_d, ax_p):
        ax.set_xticks(ticks, labels, fontsize=6.5)
        for year, xc in group_pos:
            ax.text(xc, -0.2, f"April {year}", transform=ax.get_xaxis_transform(), ha="center", fontsize=8,
                    color=S.INK_SECONDARY)
    ax_d.axhspan(PUBLISHED["m87_2017_2021_mean_diameter_uas"][0] - PUBLISHED["m87_2017_2021_mean_diameter_uas"][1],
                 PUBLISHED["m87_2017_2021_mean_diameter_uas"][0] + PUBLISHED["m87_2017_2021_mean_diameter_uas"][1],
                 color=S.C_PUBLISHED, alpha=0.15, lw=0)
    ax_d.text(ticks[-1], PUBLISHED["m87_2017_2021_mean_diameter_uas"][0] + 0.7, "EHT 2025: 43.9 ± 0.6 µas, 2017-2021",
              ha="right", fontsize=7, color=S.C_PUBLISHED)
    ax_d.set_ylabel("ring diameter, δ-ring (µas)")
    ax_d.set_title("9a  The ring's size did not change")
    ax_p.set_ylabel("brightest side (deg east of north)")
    ax_p.set_title("9b  Its bright side moved")
    ax_p.axhspan(*PUBLISHED["m87_pa_2018_deg"], color=S.ORANGE, alpha=0.12, lw=0)
    ax_p.axhspan(*PUBLISHED["m87_pa_2017_deg"], color=S.BLUE, alpha=0.08, lw=0)
    shows = [("m87_2017_096_lo", "6 Apr 2017"), ("m87_2017_101_lo", "11 Apr 2017"),
             ("m87_2018_111_b3", "21 Apr 2018"), ("m87_2018_115_b3", "25 Apr 2018")]
    for j, (key, lab) in enumerate(shows):
        ax = fig.add_subplot(gs[1, j])
        if F.load_fit(A.fit_path(key, A.FIDUCIAL)) is None:
            ax.axis("off")
            continue
        img = posterior_mean_image(key)
        _show_image(ax, img, title=f"fitted model, {lab}")
    S.caption(fig, "Same model (blurred two-mode ring + central floor + two Gaussians) fitted to each data set's closure "
                   "quantities. 2018: EHT 2024-D01-01, HOPS, 227 and 229 GHz, with the Greenland Telescope. Bottom: "
                   "posterior-mean model images. Shaded: EHT 2025 ranges of the brightness position angle for 2017 "
                   "(blue) and 2018 (orange).")
    return S.save(fig, FIGURE_DIR / "09_two_years")


def main(which=None):
    figs = {"uv": fig_uv_and_amplitudes, "cphase": fig_closure_phases, "fits": fig_ring_fits,
            "comparison": fig_model_comparison, "images": fig_images, "validation": fig_validation,
            "mass": fig_mass, "sgra": fig_sgra, "years": fig_two_years}
    for name, fn in figs.items():
        if which and name not in which:
            continue
        try:
            print("  ", fn()[0], flush=True)
        except (KeyError, FileNotFoundError, TypeError) as err:  # a result not computed yet
            print(f"   skipped {name}: {type(err).__name__} {err}", flush=True)


if __name__ == "__main__":
    import sys

    main(sys.argv[1:] or None)
