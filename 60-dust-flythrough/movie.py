"""Record the viewer's guided tour as a movie: figures/dust_flythrough.mp4.

The frames come from the same page as the interactive viewer
(figures/dust_flythrough.html), rendered by headless Chromium on the GPU
(ANGLE + Vulkan; tools/render_frames.js refuses a software renderer, which
would take ~25 s a frame) and encoded with ffmpeg (H.264). Render mode: Solid
(dust also absorbs what lies behind it; the opacity scale is a display
choice, stated on screen and in the README).
"""

from __future__ import annotations

import shutil
import subprocess

import config as C

FPS = 30
SIZE = "1280x720"
MODE = 2
FRAMES = C.CACHE_DIR / "movie_frames"
OUT = C.FIGURE_DIR / "dust_flythrough.mp4"
POSTER = C.FIGURE_DIR / "dust_flythrough_poster.png"


def make(force: bool = False) -> str:
    if OUT.exists() and not force:
        return str(OUT)
    if FRAMES.exists():
        shutil.rmtree(FRAMES)
    page = C.FIGURE_DIR / "dust_flythrough.html"
    subprocess.run(["node", str(C.PROJECT_DIR / "tools" / "render_frames.js"), str(page), str(FRAMES),
                    "--mode", str(MODE), "--fps", str(FPS), "--from", "0", "--to", "100", "--size", SIZE, "--jpeg"],
                   check=True)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", str(FRAMES / "frame_%05d.jpg"),
                    "-c:v", "libx264", "-preset", "slow", "-crf", "22", "-pix_fmt", "yuv420p",
                    "-movflags", "+faststart", str(OUT)], check=True)
    poster = sorted(FRAMES.glob("frame_*.jpg"))[int(68 * FPS)]   # the face-on view of the wave
    from PIL import Image
    Image.open(poster).save(POSTER)
    shutil.rmtree(FRAMES)
    return str(OUT)


if __name__ == "__main__":
    import sys
    print(make(force="--force" in sys.argv))
