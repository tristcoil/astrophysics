"""Keeper animation: one day of E14 flying its eccentric orbit while its clock runs slow and fast.

    python animation.py      # writes figures/galileo_clock.mp4 (~1 min)

Left: the orbit to scale in its own plane, Earth at the centre, a nominal Galileo
orbit for comparison. Middle: how fast the clock ticks relative to a ground clock
at each moment, split into gravity and speed. Right: the clock's accumulated
offset, real ESA data (relativity term added back, daily offset and rate removed)
against general relativity's prediction.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation

import clockfit as F
import config as C
import gnss
import plotstyle as S
import prepare

S.use_style()
DAY = "2020-03-02"


def frames_data(prn: str = "E14", day: str = DAY) -> dict:
    d = prepare.load(prn)
    i = list(d["date"]).index(day)
    path = next(p for p in gnss.sp3_files() if p.name.split("_")[1].startswith(
        f"{day[:4]}{np.datetime64(day).astype('datetime64[D]').item().timetuple().tm_yday:03d}"))
    sp3 = gnss.read_sp3(path)
    od = gnss.orbit_day(sp3, prn)
    od26 = gnss.orbit_day(sp3, "E26")
    # perifocal coordinates: x towards perigee, y along the motion
    el = gnss.kepler_elements(od.r, od.v)
    h = el["h"].mean(axis=0)
    h /= np.linalg.norm(h)
    k = np.argmin(el["r"])
    px = od.r[k] / np.linalg.norm(od.r[k])
    py = np.cross(h, px)
    xy = np.stack([od.r @ px, od.r @ py], axis=1) / 1e3
    h26 = np.cross(od26.r[0], od26.v[0])
    h26 /= np.linalg.norm(h26)
    e1 = od26.r[0] / np.linalg.norm(od26.r[0])
    xy26 = np.stack([od26.r @ e1, od26.r @ np.cross(h26, e1)], axis=1) / 1e3
    y, tpl = F.series_for(d, i, "full")
    f = F.fit_day(y, {}, keep_resid=True)
    ok = f.mask
    X = np.vstack([np.ones(ok.sum()), F.TC[ok]]).T
    detr = lambda z: z - np.vstack([np.ones_like(F.TC), F.TC]).T @ np.linalg.lstsq(X, z[ok], rcond=None)[0]
    phys = detr(d["clk"][i] + d["kc"][i])
    gr = detr(d["mono"][i] + d["j2"][i] + d["tide"][i] + d["kin"][i])
    T = od.terms
    c2 = C.C_LIGHT**2
    grav_rate = (C.W0_GEOID / c2 - (T["rate_mono"] + T["rate_j2"])) * 86_400e6
    speed_rate = -T["rate_kin"] * 86_400e6
    return {"xy": xy, "xy26": xy26, "hours": np.arange(289) * 300 / 3600, "phys": phys * 1e9, "gr": gr * 1e9,
            "ok": ok, "grav": grav_rate, "speed": speed_rate, "alt": (el["r"] - 6_371e3) / 1e3}


def make(out=C.FIGURE_DIR / "galileo_clock.mp4", fps: int = 24) -> str:
    D = frames_data()
    fig = plt.figure(figsize=(12.8, 5.2))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.05, 1, 1], height_ratios=[1, 1], wspace=0.35, hspace=0.45)
    ax = fig.add_subplot(gs[:, 0])
    rx = fig.add_subplot(gs[0, 1:])
    cx = fig.add_subplot(gs[1, 1:], sharex=rx)
    # orbit panel
    ax.add_patch(plt.Circle((0, 0), 6_371, color=S.BLUE_300, zorder=2))
    ax.plot(D["xy26"][:, 0], D["xy26"][:, 1], color=S.AQUA, lw=1, ls="--", label="nominal Galileo orbit (E26)")
    ax.plot(D["xy"][:, 0], D["xy"][:, 1], color=S.CONTEXT, lw=1.2, label="E14's orbit")
    sat, = ax.plot([], [], "o", color=S.BLUE, ms=9, mec=S.SURFACE, zorder=5)
    trail, = ax.plot([], [], color=S.BLUE, lw=2.2, zorder=4)
    ax.set_aspect("equal")
    lim = 34_000
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim + 7_000)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    kp = int(np.argmin(np.hypot(D["xy"][:, 0], D["xy"][:, 1])))
    ax.annotate("perigee", (D["xy"][kp, 0], D["xy"][kp, 1]), (D["xy"][kp, 0] + 2500, D["xy"][kp, 1] + 7000),
                fontsize=8, color=S.INK_SECONDARY,
                arrowprops={"arrowstyle": "-", "color": S.MUTED, "lw": 0.7})
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.06), fontsize=7.5)
    ax.set_title("E14 (GSAT-0202), 2 March 2020", loc="center")
    alt_txt = ax.text(0.02, 0.97, "", fontsize=9, color=S.INK, ha="left", va="top", zorder=6,
                      transform=ax.transAxes)
    # rate panel
    hrs = D["hours"]
    rx.plot(hrs, D["grav"], color=S.GRAVITY, lw=1.2, alpha=0.25)
    rx.plot(hrs, D["speed"], color=S.SPEED, lw=1.2, alpha=0.25)
    rx.plot(hrs, D["grav"] + D["speed"], color=S.INK, lw=1.2, alpha=0.25)
    g_line, = rx.plot([], [], color=S.GRAVITY, lw=2, label="gravity: higher = faster")
    s_line, = rx.plot([], [], color=S.SPEED, lw=2, label="speed: faster = slower")
    t_line, = rx.plot([], [], color=S.INK, lw=2, label="net")
    rx.set_ylabel("µs gained per day\nvs. a ground clock")
    rx.set_ylim(-12, 56)
    rx.legend(loc="center", fontsize=7, ncol=3, bbox_to_anchor=(0.5, 0.42))
    rx.set_title("How fast the clock is ticking", loc="left")
    # clock panel
    ok = D["ok"]
    cx.plot(hrs, D["gr"], color=S.INK, lw=1.2, alpha=0.25)
    gr_line, = cx.plot([], [], color=S.INK, lw=1.6, label="general relativity")
    pts, = cx.plot([], [], "o", color=S.BLUE, ms=2.5, label="E14's clock (ESA data)")
    cx.set_ylabel("clock ahead (+) or\nbehind (−), ns")
    cx.set_xlabel("hours (GPS time)")
    cx.set_xlim(0, 24)
    cx.set_ylim(-520, 520)
    cx.legend(loc="lower right", fontsize=7)
    cx.set_title("Where the clock is, after daily offset and rate", loc="left")
    fig.text(0.5, 0.005, "Data: ESA/ESOC final orbits and clocks. The clock's once-per-orbit swing (±0.4 µs) is "
             "half gravity, half speed; its daily average gain (~40 µs/day) is removed with the day's rate.",
             ha="center", fontsize=7.5, color=S.MUTED)
    n = len(hrs)

    def update(j):
        j = min(j, n - 1)
        sat.set_data([D["xy"][j, 0]], [D["xy"][j, 1]])
        trail.set_data(D["xy"][max(0, j - 30):j + 1, 0], D["xy"][max(0, j - 30):j + 1, 1])
        alt_txt.set_text(f"{hrs[j]:5.2f} h   {D['alt'][j]:,.0f} km up\n"
                         f"clock gains {D['grav'][j] + D['speed'][j]:+.1f} µs/day")
        g_line.set_data(hrs[:j + 1], D["grav"][:j + 1])
        s_line.set_data(hrs[:j + 1], D["speed"][:j + 1])
        t_line.set_data(hrs[:j + 1], (D["grav"] + D["speed"])[:j + 1])
        gr_line.set_data(hrs[:j + 1], D["gr"][:j + 1])
        m = ok[:j + 1]
        pts.set_data(hrs[:j + 1][m], D["phys"][:j + 1][m])
        return sat, trail, alt_txt, g_line, s_line, t_line, gr_line, pts

    anim = FuncAnimation(fig, update, frames=n + fps * 2, blit=False)
    anim.save(str(out), writer=FFMpegWriter(fps=fps, bitrate=1800, codec="libx264",
                                            extra_args=["-pix_fmt", "yuv420p", "-crf", "28"]), dpi=110)
    plt.close(fig)
    # a still of the last frame for the README
    return str(out)


if __name__ == "__main__":
    print(make())
