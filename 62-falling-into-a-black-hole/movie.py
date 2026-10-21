"""The movie: falling into Sgr A* from 1,000 Schwarzschild radii to the singularity.

    python movie.py              # render missing frames (resumable), then encode
    python movie.py --encode     # only encode existing frames

The faller drops from rest at infinity (a "raindrop") straight toward Sgr A*,
coming from the Sun's direction, and looks ahead. What they see depends only
on r / (GM/c^2); the black hole's mass sets only the clock. With H1's mass
(3.96 million Suns) GM/c^3 = 19.5 s:

* act 1 (20 s): 1,000 -> 10 Schwarzschild radii (9.5 days of the faller's time),
* act 2 (22 s): 10 -> 1 Schwarzschild radius (13 minutes),
* act 3 (26 s): from the horizon to 0.01 Schwarzschild radii, in real time.

Frame layout (1920 x 1080): ahead (stereographic, 150 deg wide) on the left; behind
(away from the hole) and the whole sky (centre ahead, rim behind) on the right;
the numbers of the moment in the right-hand column.
"""

from __future__ import annotations

import json
import math
import multiprocessing as mp
import os
import subprocess
import sys
import time
from functools import lru_cache

import numpy as np
from PIL import Image, ImageDraw, ImageFont

import config as C
import infall
import raymap
import render as RD

W, H = 1920, 1080
MAIN = 1080
SIDE = 540
COL_X = MAIN + SIDE          # info column starts here (300 px wide)
FPS = C.MOVIE_FPS
TITLE_S, END_S = 3.0, 5.0
ACT1_S, ACT2_S = 20.0, 22.0
R_START, R_ACT2, R_END = 2000.0, 20.0, 0.022     # in units of GM/c^2 (r_s = 2)
MOVIE_PATH = C.FIGURE_DIR / "falling_into_sgr_a.mp4"
FRAME_DIR = C.FRAME_DIR


@lru_cache(maxsize=1)
def bh():
    return C.black_holes()["sgra"]


def mass():
    return bh()["mass"][0]


def time_unit():
    """GM/c^3 in seconds for Sgr A*."""
    return C.gm_over_c3_s(mass())


def schedule():
    """List of frame descriptors: dict(r, act, pace) in units of M (= GM/c^2)."""
    frames = []
    for _ in range(int(TITLE_S * FPS)):
        frames.append({"r": R_START, "act": "title"})
    n1 = int(ACT1_S * FPS)
    for i in range(n1):
        x = i / n1
        frames.append({"r": R_START * (R_ACT2 / R_START) ** x, "act": "approach"})
    n2 = int(ACT2_S * FPS)
    for i in range(n2):
        x = i / n2
        frames.append({"r": R_ACT2 * (2.0 / R_ACT2) ** x, "act": "final approach"})
    # act 3: real time for Sgr A*, from the horizon until r = R_END
    tu = time_unit()
    tau_h = float(infall.rain_tau_left(2.0))
    tau_end = float(infall.rain_tau_left(R_END))
    n3 = int(round((tau_h - tau_end) * tu * FPS))
    for i in range(n3 + 1):
        tau_left = tau_h - i / FPS / tu
        frames.append({"r": float(infall.solve_r_at_tau(max(tau_left, tau_end))), "act": "inside"})
    last = frames[-1]["r"]
    for _ in range(int(END_S * FPS)):
        frames.append({"r": last, "act": "end"})
    # pace: faller seconds per movie second, from neighbouring frames of the same act
    tau = np.array([float(infall.rain_tau_left(f["r"])) for f in frames]) * tu
    for i, f in enumerate(frames):
        j = i + 1 if i + 1 < len(frames) and frames[i + 1]["act"] == f["act"] else i
        k = i - 1 if i > 0 and frames[i - 1]["act"] == f["act"] else i
        f["pace"] = (tau[k] - tau[j]) / ((j - k) / FPS) if j > k else 0.0
        f["tau_left_s"] = tau[i]
    return frames


# ---------------------------------------------------------------------------
# Numbers shown on screen
# ---------------------------------------------------------------------------

def fmt_duration(s: float) -> str:
    if s >= 86400 * 2:
        return f"{s / 86400:.1f} days"
    if s >= 3600:
        h = int(s // 3600)
        m = int((s - 3600 * h) // 60)
        return f"{h} h {m:02d} min"
    if s >= 60:
        m = int(s // 60)
        return f"{m} min {int(s - 60 * m):02d} s"
    if s >= 10:
        return f"{s:.1f} s"
    return f"{s:.2f} s"


def fmt_sci(x: float, digits: int = 2) -> str:
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    if -2 <= e <= 4:
        return f"{x:,.{max(0, digits - 1 - e)}f}"
    m = x / 10**e
    sup = str(e).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{m:.{digits - 1}f}×10{sup}"


def frame_numbers(r: float) -> dict:
    m = mass()
    gm = C.GM_SUN * m
    r_m = r * C.gm_over_c2_m(m)
    beta = math.sqrt(2.0 / r)
    stretch = infall.tidal_stretch_si(r_m, m, C.BODY_LENGTH_M, C.GM_SUN) / C.G0
    out = {"r_rs": r / 2.0, "r_km": r_m / 1e3, "tau_left_s": float(infall.rain_tau_left(r)) * time_unit(),
           "stretch_g": float(stretch), "beta": beta}
    if r > 2.0:
        out["one_plus_z"] = 1.0 / (1.0 - beta)
    return out


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------

@lru_cache(maxsize=4)
def font(size: int, bold: bool = False):
    import matplotlib.font_manager as fm
    path = fm.findfont(fm.FontProperties(family="DejaVu Sans", weight="bold" if bold else "normal"))
    return ImageFont.truetype(path, size)


INK = (236, 236, 232)
MUTE = (150, 150, 146)
ACCENT = (140, 180, 255)


def radius_bar(draw, x0, y0, h, r):
    """Vertical log-radius bar from 1,000 r_s (top) to 0.01 r_s (bottom) with the faller's dot."""
    lo, hi = -2.0, 3.0
    def ypos(r_rs):
        return y0 + (hi - math.log10(r_rs)) / (hi - lo) * h
    draw.line([(x0, y0), (x0, y0 + h)], fill=(90, 90, 88), width=2)
    for r_rs, label, dy in ((1000, "1,000 rₛ", 0), (100, "100", 0), (10, "10", 0),
                            (1.5, "photon sphere", -10), (1, "horizon", 6), (0.1, "0.1", 0), (0.01, "0.01", 0)):
        y = ypos(r_rs)
        col = ACCENT if label == "horizon" else MUTE
        draw.line([(x0 - 6, y), (x0 + 6, y)], fill=col, width=2 if label == "horizon" else 1)
        draw.text((x0 + 12, y - 8 + dy), label, font=font(14), fill=col)
    y = ypos(max(r / 2.0, 0.0101))
    draw.ellipse([x0 - 7, y - 7, x0 + 7, y + 7], fill=(255, 200, 90))


def info_panel(img: Image.Image, f: dict):
    d = ImageDraw.Draw(img)
    x = COL_X + 18
    n = frame_numbers(f["r"])
    y = 24
    d.text((x, y), "Falling into", font=font(22), fill=INK); y += 28
    d.text((x, y), "Sgr A*", font=font(34, True), fill=INK); y += 46
    d.text((x, y), f"{mass() / 1e6:.2f} million Suns (H1)", font=font(15), fill=MUTE); y += 34

    def row(label, value, col=INK, sub=None):
        nonlocal y
        d.text((x, y), label, font=font(15), fill=MUTE); y += 20
        d.text((x, y), value, font=font(22, True), fill=col); y += 30
        if sub:
            d.text((x, y), sub, font=font(15), fill=MUTE); y += 20
        y += 8

    row("distance from the centre", f"{fmt_sci(n['r_rs'], 3)} rₛ", sub=f"= {fmt_sci(n['r_km'], 3)} km")
    row("your own clock: time left", fmt_duration(n["tau_left_s"]), (255, 200, 90))
    if f["r"] > 2.0:
        row("speed past hovering observers", f"{n['beta']:.3f} c")
        row("light you send outward", f"arrives {fmt_sci(n['one_plus_z'], 3)}× redder")
    else:
        row("speed past hovering observers", "no one can hover")
        row("light you send outward", "can no longer escape", ACCENT)
    row("head-to-foot stretch (2 m)", f"{fmt_sci(n['stretch_g'], 2)} g")
    pace = f.get("pace", 0.0)
    if f["act"] in ("title", "end"):
        ptxt = "paused"
    elif abs(pace - 1.0) < 0.02:
        ptxt = "real time"
    else:
        ptxt = f"1 s here = {fmt_duration(pace)}"
    row("movie speed", ptxt)
    radius_bar(d, x + 8, y + 10, 230, f["r"])
    yb = H - 92
    for line in ("SIMULATION (general relativity,", "Schwarzschild). Sky: Gaia DR3 +",
                 "Hipparcos, as seen from the Sun.", "No gas: a clean view."):
        d.text((x, yb), line, font=font(13), fill=MUTE); yb += 17


def labels(img: Image.Image, f: dict):
    d = ImageDraw.Draw(img)
    d.text((16, 14), "AHEAD: toward the black hole", font=font(18, True), fill=INK)
    d.text((MAIN + 12, 10), "BEHIND you", font=font(15, True), fill=INK)
    d.text((MAIN + 12, SIDE + 10), "WHOLE SKY", font=font(15, True), fill=INK)
    d.text((MAIN + 12, SIDE + 30), "centre: ahead, rim: behind", font=font(13), fill=MUTE)
    if f["act"] == "title":
        msg = ["What would you see falling into the black hole", "at the centre of the Milky Way?"]
        for i, line in enumerate(msg):
            d.text((60, 860 + 34 * i), line, font=font(28, True), fill=INK)
    if f["act"] == "inside":
        d.text((16, 40), "Inside the horizon. These last seconds run in real time.", font=font(16), fill=ACCENT)
    if f["act"] == "end":
        left = fmt_duration(frame_numbers(f["r"])["tau_left_s"])
        msg = [f"{left} before the singularity, the outside universe is squeezed into a thin,",
               "blueshifted band around the faller's waist, while tides stretch them",
               "lengthwise and crush them sideways. At r = 0 general relativity stops working."]
        for i, line in enumerate(msg):
            d.text((40, 930 + 28 * i), line, font=font(20), fill=INK)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

_SCENE = None


def views(scene, r):
    zhat, up = RD.default_frame()
    table = raymap.RayTable(r, 0.0)
    ahead = RD.Camera(MAIN, MAIN, forward=-zhat, up=up, fov_deg=150)
    behind = RD.Camera(SIDE, SIDE, forward=zhat, up=up, fov_deg=150)
    whole = RD.Camera(SIDE, SIDE, forward=-zhat, up=up, projection="equidistant")
    return [RD.tone_map(RD.render_view(scene, cam, table, zhat)) for cam in (ahead, behind, whole)]


def render_frame(i: int, f: dict, cache: dict | None = None):
    path = FRAME_DIR / f"frame_{i:05d}.png"
    if path.exists():
        return path
    key = round(f["r"], 12)
    if cache is not None and key in cache:
        a, b, c = cache[key]
    else:
        a, b, c = views(_SCENE, f["r"])
        if cache is not None:
            cache.clear()
            cache[key] = (a, b, c)
    canvas = np.zeros((H, W, 3), np.uint8)
    canvas[:, :MAIN] = a
    canvas[:SIDE, MAIN:MAIN + SIDE] = b
    canvas[SIDE:, MAIN:MAIN + SIDE] = c
    canvas[:, MAIN:MAIN + 1] = 40
    canvas[SIDE, MAIN:MAIN + SIDE] = 40
    img = Image.fromarray(canvas)
    info_panel(img, f)
    labels(img, f)
    tmp = path.with_suffix(".tmp.png")
    img.save(tmp, compress_level=1)
    os.replace(tmp, path)
    return path


def _worker(args):
    chunk = args
    cache: dict = {}
    for i, f in chunk:
        render_frame(i, f, cache)
    return len(chunk)


def _init_worker():
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except ImportError:
        pass


def render_all(n_workers: int = 6):
    global _SCENE
    frames = schedule()
    todo = [(i, f) for i, f in enumerate(frames) if not (FRAME_DIR / f"frame_{i:05d}.png").exists()]
    print(f"{len(frames)} frames, {len(todo)} to render with {n_workers} workers", flush=True)
    if not todo:
        return frames
    _SCENE = RD.Scene()
    zhat, _ = RD.default_frame()
    _SCENE.star_geometry(RD.unit(zhat))     # computed once, shared by the forked workers
    # contiguous chunks so that repeated still frames (title, end) reuse one render
    size = max(1, len(todo) // (n_workers * 8))
    chunks = [todo[k:k + size] for k in range(0, len(todo), size)]
    ctx = mp.get_context("fork")
    t0 = time.time()
    done = 0
    with ctx.Pool(n_workers, initializer=_init_worker) as pool:
        for n in pool.imap_unordered(_worker, chunks):
            done += n
            el = time.time() - t0
            print(f"  {done}/{len(todo)} frames, {el / 60:.1f} min, ~{el / done * (len(todo) - done) / 60:.0f} min left",
                  flush=True)
    return frames


def encode(crf: int = 26):
    frames = schedule()
    n = len(frames)
    missing = [i for i in range(n) if not (FRAME_DIR / f"frame_{i:05d}.png").exists()]
    if missing:
        raise RuntimeError(f"{len(missing)} frames missing, e.g. {missing[:5]}")
    cmd = ["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", str(FRAME_DIR / "frame_%05d.png"),
           "-frames:v", str(n), "-c:v", "libx264", "-preset", "slow", "-crf", str(crf), "-pix_fmt", "yuv420p",
           "-movflags", "+faststart", str(MOVIE_PATH)]
    subprocess.run(cmd, check=True)
    meta = {"frames": n, "fps": FPS, "seconds": n / FPS, "bytes": MOVIE_PATH.stat().st_size,
            "mass_msun": mass(), "gm_c3_s": time_unit(), "crf": crf}
    (C.CACHE_DIR / "movie_meta.json").write_text(json.dumps(meta, indent=1))
    print(f"wrote {MOVIE_PATH} ({meta['bytes'] / 1e6:.1f} MB, {meta['seconds']:.1f} s)")
    return meta


if __name__ == "__main__":
    if "--encode" not in sys.argv:
        nw = int(os.environ.get("N5_WORKERS", "6"))
        render_all(nw)
    encode()
