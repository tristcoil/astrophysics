"""Animations: a year of Earth with its "dot", and the Moon crossing Earth's face.

* ``earth_2024.mp4``: every 2024 natural-colour image (~5,500 frames, 30 per
  second, ~3 minutes). Beside the planet: the same moment as a single point of
  light, coloured as EPIC's camera would record it, and the dot's brightness in
  three colours over the previous five days.
* ``earth_one_day.gif``: one summer day of the same, for pages where video
  does not play.
* ``moon_transit_2015.gif``, ``_2016.gif``, ``_2021.gif``: the lunar transits
  from the EPIC team's galleries.

Frames are composed with PIL and piped to ffmpeg (H.264). The light-curve
strip is drawn once with matplotlib (antialiased) and cropped per frame.
"""

from __future__ import annotations

import subprocess

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from PIL import Image, ImageDraw, ImageFont  # noqa: E402

import epicdata  # noqa: E402
from config import FIGURE_DIR, GALLERY_DIR  # noqa: E402
from figures import dot_colour  # noqa: E402

W, H = 960, 540
EARTH_PX = 540
SURFACE = (26, 26, 25)  # dark surface of the house palette
INK = (255, 255, 255)
INK2 = (195, 194, 183)
MUTED = (137, 135, 129)
GRID = (60, 60, 58)
#: Channel colours, dark-mode steps of the validated palette.
CHAN_HEX = {"red": "#e66767", "green": "#2f9e2f", "blue": "#3987e5"}
FPS = 30
DAYS_SHOWN = 5.0
STRIP_PX_PER_DAY = 70
STRIP_H = 190
YLIM = 25.0  # percent


def _font(size: int, bold: bool = False):
    prop = font_manager.FontProperties(family="DejaVu Sans", weight="bold" if bold else "normal")
    return ImageFont.truetype(font_manager.findfont(prop), size)


def _earth(image: str, size: int) -> Image.Image:
    with Image.open(epicdata.image_path("jpg", image)) as im:
        return im.convert("RGB").resize((size, size), Image.LANCZOS)


def _strip(frames: pd.DataFrame, t0: pd.Timestamp, n_days: float) -> np.ndarray:
    """The whole period's light curve as one long antialiased image."""
    width = int(n_days * STRIP_PX_PER_DAY)
    dpi = 100
    fig = plt.figure(figsize=(width / dpi, STRIP_H / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    fig.patch.set_facecolor(np.array(SURFACE) / 255)
    ax.set_facecolor(np.array(SURFACE) / 255)
    x = (frames["t"] - t0).dt.total_seconds().to_numpy() / 86400.0
    ok = frames["ok"].to_numpy()
    for d in range(int(n_days) + 1):
        ax.axvline(d, color=np.array(GRID) / 255, lw=1.0)
    for yv in (-20, -10, 0, 10, 20):
        ax.axhline(yv, color=np.array(GRID) / 255, lw=0.8 if yv else 1.2)
    for ch in ("red", "green", "blue"):
        y = 100 * (frames[f"F_{ch}"] / frames.loc[ok, f"F_{ch}"].mean() - 1).to_numpy()
        yy = np.where(ok, y, np.nan)
        gaps = np.diff(x, prepend=x[0]) > 3.0 / 24
        yy[gaps] = np.nan  # break the line across gaps
        ax.plot(x, yy, color=CHAN_HEX[ch], lw=1.6, solid_capstyle="round")
        ax.plot(x[ok], y[ok], ls="none", marker="o", ms=2.2, color=CHAN_HEX[ch])
    ax.set_xlim(0, n_days)
    ax.set_ylim(-YLIM, YLIM)
    ax.set_axis_off()
    fig.canvas.draw()
    img = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
    plt.close(fig)
    return img


class Composer:
    def __init__(self, frames: pd.DataFrame, title: str):
        self.frames = frames.sort_values("t").reset_index(drop=True)
        self.t0 = self.frames["t"].iloc[0].normalize() - pd.Timedelta(days=DAYS_SHOWN)
        n_days = (self.frames["t"].iloc[-1] - self.t0).total_seconds() / 86400 + 1
        self.strip = _strip(self.frames, self.t0, n_days)
        ok = self.frames["ok"]
        self.mean = {ch: self.frames.loc[ok, f"F_{ch}"].mean() for ch in ("red", "green", "blue")}
        self.title = title
        self.f_title, self.f_big = _font(17, True), _font(16)
        self.f_mid, self.f_small, self.f_tiny = _font(13), _font(11), _font(10)
        # a Gaussian point-spread function for the dot
        yy, xx = np.mgrid[-60:61, -60:61]
        self.psf = np.exp(-(xx**2 + yy**2) / (2 * 14.0**2))

    def frame(self, k: int) -> Image.Image:
        r = self.frames.iloc[k]
        img = Image.new("RGB", (W, H), SURFACE)
        img.paste(_earth(r["image"], EARTH_PX), (0, 0))
        d = ImageDraw.Draw(img)
        x0 = EARTH_PX + 18
        d.text((x0, 14), self.title, font=self.f_title, fill=INK)
        d.text((x0, 42), f"{r['t']:%Y-%m-%d  %H:%M} UTC", font=self.f_big, fill=INK)
        lon = r["lon"]
        lon_s = f"{abs(lon):.0f}°{'E' if lon >= 0 else 'W'}"
        lat_s = f"{abs(r['lat']):.0f}°{'N' if r['lat'] >= 0 else 'S'}"
        d.text((x0, 66), f"facing {lon_s}, {lat_s} · phase angle {r['phase_deg']:.1f}°",
               font=self.f_mid, fill=INK2)
        # the dot
        box = 150
        bx, by = x0, 96
        dot = np.zeros((box, box, 3))
        col = np.array(dot_colour(r["F_red"], r["F_green"], r["F_blue"]))
        c = box // 2
        dot[c - 60:c + 61, c - 60:c + 61] = self.psf[..., None] * col
        img.paste(Image.fromarray((255 * np.clip(dot, 0, 1)).astype(np.uint8)), (bx, by))
        tx = bx + box + 14
        d.text((tx, by + 4), "Earth as a dot", font=self.f_mid, fill=INK)
        d.text((tx, by + 24), "the whole planet in", font=self.f_small, fill=INK2)
        d.text((tx, by + 39), "one point of light", font=self.f_small, fill=INK2)
        if r["ok"]:
            br = 100 * (np.mean([r[f"F_{ch}"] / self.mean[ch] for ch in ("red", "green", "blue")]) - 1)
            d.text((tx, by + 66), f"brightness {br:+.0f}%", font=self.f_small, fill=INK2)
            d.text((tx, by + 82), f"blue/red {r['blue_red']:.2f}", font=self.f_small, fill=INK2)
        elif r["t"].strftime("%Y-%m-%d") == "2024-04-08":
            d.text((tx, by + 66), "total solar eclipse:", font=self.f_small, fill=INK2)
            d.text((tx, by + 82), "the Moon's shadow", font=self.f_small, fill=INK2)
        else:
            d.text((tx, by + 66), "frame left out", font=self.f_small, fill=MUTED)
            d.text((tx, by + 82), "of the curve", font=self.f_small, fill=MUTED)
        # light curve: the last DAYS_SHOWN days up to now
        lx, ly = x0, 290
        lw = int(DAYS_SHOWN * STRIP_PX_PER_DAY)
        xnow = int((r["t"] - self.t0).total_seconds() / 86400 * STRIP_PX_PER_DAY)
        crop = self.strip[:, max(0, xnow - lw):xnow]
        if crop.shape[1] < lw:
            crop = np.pad(crop, ((0, 0), (lw - crop.shape[1], 0), (0, 0)), constant_values=SURFACE[0])
        img.paste(Image.fromarray(crop), (lx, ly))
        d.line([(lx + lw - 1, ly), (lx + lw - 1, ly + STRIP_H)], fill=INK, width=1)
        d.text((lx, ly - 22), "brightness, % from the year's mean (last 5 days)", font=self.f_small, fill=INK2)
        for yv in (-20, 0, 20):
            yp = ly + STRIP_H / 2 - yv / YLIM * STRIP_H / 2
            d.text((lx + lw + 6, yp - 7), f"{yv:+d}", font=self.f_tiny, fill=MUTED)
        # day labels under the strip
        day_now = r["t"].normalize()
        for back in range(0, int(DAYS_SHOWN) + 1):
            day = day_now - pd.Timedelta(days=back)
            xp = lx + lw - (r["t"] - day).total_seconds() / 86400 * STRIP_PX_PER_DAY
            if lx <= xp <= lx + lw:
                d.text((xp + 2, ly + STRIP_H + 4), f"{day:%d %b}", font=self.f_tiny, fill=MUTED)
        ky = ly + STRIP_H + 24
        for i, (ch, lab) in enumerate((("red", "680 nm"), ("green", "551 nm"), ("blue", "443 nm"))):
            kx = lx + 8 + i * 110
            d.line([(kx, ky + 7), (kx + 22, ky + 7)], fill=CHAN_HEX[ch], width=3)
            d.text((kx + 28, ky), lab, font=self.f_small, fill=INK2)
        d.text((x0, H - 18), "DSCOVR/EPIC natural colour · NASA EPIC Team · Astrophysics Atlas A4",
               font=self.f_tiny, fill=MUTED)
        return img


def write_mp4(composer: Composer, path, *, crf: int = 30, scale: float = 1.0) -> str:
    w, h = int(W * scale) // 2 * 2, int(H * scale) // 2 * 2
    cmd = ["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24",
           "-s", f"{w}x{h}", "-r", str(FPS), "-i", "-", "-c:v", "libx264", "-preset", "slow",
           "-crf", str(crf), "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(path)]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    for k in range(len(composer.frames)):
        im = composer.frame(k)
        if scale != 1.0:
            im = im.resize((w, h), Image.LANCZOS)
        proc.stdin.write(im.tobytes())
    proc.stdin.close()
    if proc.wait() != 0:
        raise RuntimeError("ffmpeg failed")
    return str(path)


def write_gif(images: list, path, fps: float) -> str:
    pal = [im.convert("P", palette=Image.ADAPTIVE, colors=255) for im in images]
    pal[0].save(path, save_all=True, append_images=pal[1:], duration=int(1000 / fps), loop=0,
                optimize=True)
    return str(path)


def transit_gif(key: str, folder: str, size: int = 420, fps: float = 3.0) -> str:
    files = sorted((GALLERY_DIR / folder).glob("*.png"))
    from moon import frame_time

    f_big, f_small = _font(15, True), _font(11)
    frames = []
    for p in files:
        with Image.open(p) as im:
            im = im.convert("RGB").resize((size, size), Image.LANCZOS)
        canvas = Image.new("RGB", (size, size + 44), SURFACE)
        canvas.paste(im, (0, 0))
        d = ImageDraw.Draw(canvas)
        d.text((10, size + 4), f"The Moon's far side crossing Earth, {frame_time(p.name):%d %b %Y %H:%M} UTC",
               font=f_small, fill=INK)
        d.text((10, size + 22), "DSCOVR/EPIC from 1.5 million km · NASA EPIC Team gallery", font=f_small,
               fill=MUTED)
        frames.append(canvas)
    return write_gif(frames, FIGURE_DIR / f"moon_transit_{key}.gif", fps)


def make_all(results, art) -> list:
    jpg = art["jpg"]
    year = jpg[jpg["t"].dt.year == 2024]
    year = year[np.isfinite(year["F_red"])]
    comp = Composer(year, "EARTH AS AN EXOPLANET · 2024")
    written = [write_mp4(comp, FIGURE_DIR / "earth_2024.mp4")]
    # one day as a GIF (half size)
    day = year[year["t"].dt.date == pd.Timestamp("2024-06-17").date()]
    idx = [comp.frames.index[comp.frames["image"] == im][0] for im in day["image"]]
    imgs = [comp.frame(k).resize((W // 2, H // 2), Image.LANCZOS) for k in idx]
    written.append(write_gif(imgs, FIGURE_DIR / "earth_one_day.gif", fps=5))
    for key, folder in (("2015", "lunar_transit_2015"), ("2016", "lunar_transit_2016"),
                        ("2021", "lunar_transit_2021")):
        written.append(transit_gif(key, folder))
    return written
