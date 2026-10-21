"""The keeper: HR 8799's four planets moving on their fitted orbits, 1998-2050.

    python animate.py            # MP4 + GIF in figures/

Each frame shows the planets where the fitted orbits put them (median of 200
posterior orbits, with the spread as a faint cloud), a fading trail, and every
measured position appearing on the date it was taken. After the last
measurement (Keck, 2025 Oct 9) the motion is a prediction, and the frame says
so. Orbits: the ALL data set (1998-2025) - see run.py for which model.
East is to the left, as on the sky.
"""

from __future__ import annotations

import warnings

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter, PillowWriter
from astropy.time import Time

import astrometry as A
import config as C
import plotstyle as S
import predict

warnings.filterwarnings("ignore", message=".*dubious year.*")
S.use_style()

YEARS = (1998.0, 2050.0)
STEP_YR = 0.1
N_ORBITS = 200
TRAIL_YR = 6.0


def _mjd(year):
    return Time(np.asarray(year, float), format="decimalyear").mjd


def track(model, samples, years):
    """Median and 1-sigma spread of each planet's position at the given years."""
    sky = predict.predicted_positions(model, samples, _mjd(years))
    out = {}
    for p in C.PLANETS:
        ra, dec = sky[p]
        out[p] = dict(ra=np.median(ra, axis=0), dec=np.median(dec, axis=0), ra_all=ra, dec_all=dec)
    return out


def full_orbits(model, samples, n=40):
    """n complete orbits per planet (over one period each) for the faint ellipses."""
    phys = model.physical(samples[:n])
    out = {}
    for p in C.PLANETS:
        a, e, inc, aop, pan, tau = phys["elements"][p]
        P = np.median(np.asarray(_period_days(a, phys, p)))
        t = _mjd(2000.0) + np.linspace(0, P, 400)[None, :]
        out[p] = t
    sky = {}
    for p in C.PLANETS:
        sk = predict.predicted_positions(model, samples[:n], out[p][0])
        sky[p] = sk[p]
    return sky


def _period_days(a, phys, p):
    import kepler as K
    masses = phys["masses"]
    m = phys["mstar"] + sum(masses[q] for q in C.INNER_TO_OUTER[: C.INNER_TO_OUTER.index(p) + 1])
    return K.period_days(a, m)


def make(model, samples, stem="hr8799_orbits_1998_2050", fps=30, label="", data_sources=C.ALL_SOURCES,
         years_range=YEARS, step=STEP_YR, gif=True):
    rng = np.random.default_rng(C.RANDOM_SEED)
    s = samples[rng.choice(len(samples), min(N_ORBITS, len(samples)), replace=False)]
    years = np.arange(years_range[0], years_range[1] + 1e-9, step)
    tr = track(model, s, years)
    ell = full_orbits(model, s, n=30)
    data = A.load(data_sources)
    last_data = data.year.max()

    fig = plt.figure(figsize=(7.2, 7.2), dpi=150)
    ax = fig.add_axes([0.07, 0.07, 0.86, 0.80])
    ax.set_aspect("equal")
    ext = np.array([[ell[p][0].min(), ell[p][0].max(), ell[p][1].min(), ell[p][1].max()] for p in C.PLANETS])
    xlo, xhi, ylo, yhi = ext[:, 0].min(), ext[:, 1].max(), ext[:, 2].min(), ext[:, 3].max()
    pad = 0.04 * max(xhi - xlo, yhi - ylo)
    xlo, xhi, ylo, yhi = xlo - pad, xhi + pad, ylo - pad, yhi + pad
    ax.set_xlim(xhi, xlo)                     # east to the left
    ax.set_ylim(ylo, yhi)
    ax.set_xlabel("ΔRA (mas, east to the left)")
    ax.set_ylabel("ΔDec (mas, north up)")
    for p in C.PLANETS:
        ra, dec = ell[p]
        for k in range(ra.shape[0]):
            ax.plot(ra[k], dec[k], color=S.PLANET_COLOR[p], lw=0.5, alpha=0.06, zorder=1)
    ax.scatter([0], [0], marker="*", s=220, color=S.INK, zorder=6, edgecolors=S.SURFACE, linewidths=1)
    # scale bar: 20 au at the Gaia distance
    plx = C.DR3_PARALLAX[0]
    x0, y0 = xlo + 0.04 * (xhi - xlo), ylo + 0.04 * (yhi - ylo)
    ax.plot([x0, x0 + 20 * plx], [y0, y0], color=S.INK_SECONDARY, lw=2)
    ax.text(x0 + 10 * plx, y0 + 0.015 * (yhi - ylo), "20 au", ha="center", va="bottom", fontsize=8,
            color=S.INK_SECONDARY)
    title = fig.text(0.07, 0.95, "HR 8799: four planets orbiting their star", fontsize=13, fontweight="bold")
    sub = fig.text(0.07, 0.932, label, fontsize=8, color=S.INK_SECONDARY, va="top", linespacing=1.4)
    year_txt = fig.text(0.93, 0.95, "", fontsize=16, fontweight="bold", ha="right")
    mode_txt = fig.text(0.93, 0.915, "", fontsize=8.5, ha="right", color=S.INK_SECONDARY)

    dots, trails, clouds, labels = {}, {}, {}, {}
    for p in C.PLANETS:
        c = S.PLANET_COLOR[p]
        clouds[p] = ax.scatter([], [], s=2, color=c, alpha=0.15, lw=0, zorder=2)
        trails[p], = ax.plot([], [], color=c, lw=2, alpha=0.9, zorder=3, solid_capstyle="round")
        dots[p] = ax.scatter([], [], s=90, color=c, edgecolors=S.SURFACE, linewidths=2, zorder=7)
        labels[p] = ax.text(0, 0, p, fontsize=11, fontweight="bold", color=S.INK, zorder=8)
    seen = ax.scatter([], [], s=9, color=S.INK, edgecolors=S.SURFACE, linewidths=0.5, zorder=5)
    ax.text(0.01, 0.01, "dots: every measured position, shown from the date it was taken\n"
            "colour cloud: spread of 200 orbits allowed by the data", transform=ax.transAxes,
            fontsize=7.5, color=S.MUTED, va="bottom")

    n_trail = int(round(TRAIL_YR / step))

    def draw(i):
        yr = years[i]
        for p in C.PLANETS:
            ra, dec = tr[p]["ra"], tr[p]["dec"]
            j0 = max(0, i - n_trail)
            trails[p].set_data(ra[j0:i + 1], dec[j0:i + 1])
            dots[p].set_offsets([[ra[i], dec[i]]])
            clouds[p].set_offsets(np.column_stack([tr[p]["ra_all"][:, i], tr[p]["dec_all"][:, i]]))
            r = np.hypot(ra[i], dec[i])
            labels[p].set_position((ra[i] * (1 + 70 / r), dec[i] * (1 + 70 / r)))
        d = data[data.year <= yr]
        seen.set_offsets(np.column_stack([d.ra, d.dec]) if len(d) else np.empty((0, 2)))
        year_txt.set_text(f"{yr:.1f}")
        mode_txt.set_text("fitted to measurements" if yr <= last_data else "prediction")
        return []

    out = []
    w = FFMpegWriter(fps=fps, bitrate=6000, codec="libx264", extra_args=["-pix_fmt", "yuv420p"])
    path = C.FIGURES / f"{stem}.mp4"
    with w.saving(fig, str(path), dpi=150):
        for i in range(len(years)):
            draw(i)
            w.grab_frame()
    out.append(path)
    if gif:                                   # a lighter GIF: every 3rd frame, smaller
        gw = PillowWriter(fps=12)
        gpath = C.FIGURES / f"{stem}.gif"
        with gw.saving(fig, str(gpath), dpi=70):
            for i in range(0, len(years), 3):
                draw(i)
                gw.grab_frame()
        out.append(gpath)
    # one still for the README (2026)
    i26 = int(np.argmin(abs(years - 2026.0)))
    draw(i26)
    fig.savefig(C.FIGURES / f"{stem}_2026.png", dpi=150)
    plt.close(fig)
    return out


def main():
    import fits
    import report
    f = fits.load(report.FINAL_FIT)
    n = len(A.load(C.ALL_SOURCES))
    label = (f"Orbits fitted to {n} positions measured 1998-2025\n"
             f"coplanar, near 1:2:4:8 periods; distance {1000 / C.DR3_PARALLAX[0]:.1f} pc (Gaia DR3)")
    for p in make(f.model, f.equal, label=label):
        print(p)


if __name__ == "__main__":
    main()
