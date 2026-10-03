"""The week as a time-lapse: the Sun turning in 193 A, with the X-rays underneath.

Every AIA 193 A image of the week (one every 3 hours) becomes one frame, all on
the same brightness scale, so what changes is the Sun: regions carried across the
disc by rotation, loops brightening and fading, the coronal hole crossing the
centre. Under the image, the GOES 1-8 A flux for the week with a cursor at the
frame's time: every spike is a flare.

Writes ``sun_week_193.mp4`` (for watching) and a smaller ``sun_week_193.gif`` (for
the README, where video does not play).
"""

from __future__ import annotations

import warnings

import matplotlib

matplotlib.use("Agg")

import astropy.units as u  # noqa: E402
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from astropy.visualization import ImageNormalize, LogStretch  # noqa: E402
from matplotlib import animation  # noqa: E402

import config as C  # noqa: E402
import fetch_data  # noqa: E402
import limb  # noqa: E402
import plotstyle as ps  # noqa: E402
import sdo  # noqa: E402

warnings.filterwarnings("ignore")

import sunpy.visualization.colormaps  # noqa: E402,F401

FPS = 8


def week_frames(paths: C.EpochPaths) -> list:
    manifest = fetch_data.load_manifest(paths)
    rels = [v for _, v in sorted(manifest["week"]["aia"].items()) if v]
    return [sdo.load_aia(paths.data / r) for r in rels]


def common_norm(maps) -> ImageNormalize:
    """One scale for every frame, from the middle frame, so brightness is comparable."""
    m = maps[len(maps) // 2]
    d = np.asarray(m.data, dtype=float)
    r = limb.fast_radius(m)
    vmax = float(np.nanpercentile(d[r < 0.98], 99.9))
    vmin = max(float(np.nanpercentile(d[(r > 1.15) & (r < 1.25)], 20)), vmax * 1e-3)
    return ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch(1000), clip=True)


def render(paths: C.EpochPaths, xrs: pd.DataFrame, flare_peak: str | None, out_dir,
           *, size_in: float = 6.4, dpi: int = 120, gif_dpi: int = 48) -> list:
    maps = week_frames(paths)
    norm = common_norm(maps)
    lim = 1180
    fig = plt.figure(figsize=(size_in, size_in * 1.27), facecolor=ps.DARK_PAGE)
    ax = fig.add_axes([0, 0.24, 1, 0.76])
    ax.set_facecolor(ps.DARK_PAGE)
    ax.set_axis_off()
    m0 = maps[0]
    bl, tr = m0.bottom_left_coord, m0.top_right_coord
    ext = [bl.Tx.to_value(u.arcsec), tr.Tx.to_value(u.arcsec),
           bl.Ty.to_value(u.arcsec), tr.Ty.to_value(u.arcsec)]
    im = ax.imshow(maps[0].data, origin="lower", cmap="sdoaia193", norm=norm, extent=ext,
                   interpolation="lanczos")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    title = ax.text(0.03, 0.97, "", transform=ax.transAxes, color=ps.DARK_INK, fontsize=11,
                    fontweight="bold", va="top")
    ax.text(0.03, 0.915, "SDO/AIA 193 Å (Fe XII, 1.6 MK corona) · one frame every 3 h",
            transform=ax.transAxes, color=ps.DARK_INK_SECONDARY, fontsize=7.5, va="top")
    ax.text(0.97, 0.03, "solar north up · the Sun turns left to right", transform=ax.transAxes,
            color=ps.DARK_MUTED, fontsize=7, ha="right", va="bottom")

    axg = fig.add_axes([0.1, 0.05, 0.86, 0.15])
    axg.set_facecolor(ps.DARK_PAGE)
    axg.plot(xrs.index, xrs["long"], color=ps.DARK_BLUE, lw=0.9)
    axg.set_yscale("log")
    axg.set_ylim(1e-7, 1e-5)
    axg.set_xlim(xrs.index[0], xrs.index[-1])
    for sp_ in ("top", "right"):
        axg.spines[sp_].set_visible(False)
    for sp_ in ("left", "bottom"):
        axg.spines[sp_].set_color(ps.DARK_MUTED)
    axg.tick_params(colors=ps.DARK_MUTED, labelsize=7)
    axg.xaxis.set_major_locator(mdates.DayLocator())
    axg.xaxis.set_major_formatter(mdates.DateFormatter("%-d %b"))
    axg.set_yticks([1e-7, 1e-6, 1e-5])
    axg.set_yticklabels(["B", "C", "M"])
    axg.grid(False)
    axg.text(0.0, 1.08, "GOES X-rays, 1–8 Å: each spike is a flare", transform=axg.transAxes,
             color=ps.DARK_INK_SECONDARY, fontsize=7.5)
    if flare_peak:
        axg.axvline(pd.Timestamp(flare_peak), color=ps.DARK_MUTED, lw=0.8, ls=":")
    cursor = axg.axvline(xrs.index[0], color=ps.DARK_INK, lw=1.2)

    def draw(k):
        m = maps[k]
        im.set_data(m.data)
        t = pd.Timestamp(m.date.isot, tz="UTC")
        title.set_text(f"{t:%a %-d %B %Y, %H:%M} UTC")
        cursor.set_xdata([t, t])
        return im, title, cursor

    anim = animation.FuncAnimation(fig, draw, frames=len(maps), blit=False)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    mp4 = out_dir / "sun_week_193.mp4"
    anim.save(mp4, writer=animation.FFMpegWriter(fps=FPS, bitrate=2400,
                                                 extra_args=["-pix_fmt", "yuv420p"]),
              dpi=dpi, savefig_kwargs={"facecolor": ps.DARK_PAGE})
    written.append(str(mp4))
    gif = out_dir / "sun_week_193.gif"
    anim.save(gif, writer=animation.PillowWriter(fps=FPS), dpi=gif_dpi,
              savefig_kwargs={"facecolor": ps.DARK_PAGE})
    written.append(str(gif))
    plt.close(fig)
    return written
