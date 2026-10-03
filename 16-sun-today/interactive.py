"""A self-contained web page for exploring the Sun of the epoch.

Two parts, both working offline from one HTML file (all images embedded):

1. **Twelve views of the Sun**: every AIA channel, the magnetogram and a
   three-colour composite, picked from a temperature ladder, each with a note on
   what it shows and why.
2. **The week in 193 A**: a time-lapse to scrub with a slider, play, or jump
   through by clicking the X-ray curve underneath.
"""

from __future__ import annotations

import base64
import io
import json
import warnings
from datetime import datetime

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

import config as C
import figures
import limb
import movie

warnings.filterwarnings("ignore")

VIEW_PX = 720
WEEK_PX = 440
LIM_ARCSEC = 1150

#: What each view shows, in plain words. {coronal_hole} etc. are filled from the data.
NOTES = {
    "hmi": "The magnetic field at the visible surface, measured from the Zeeman splitting "
           "of an iron line (HMI). White: field pointing toward us; black: away. Active "
           "regions are pairs of opposite polarity; the gas in every other view is shaped by "
           "this field.",
    4500: "White light from the photosphere, the 5,800 K layer we call the surface: the "
          "deepest we can see. Sunspots are the small dark dots; the disc darkens toward its "
          "edge because there we see higher, cooler layers. The faint diagonal streaks are "
          "stray light inside the telescope, not the Sun.",
    1700: "Ultraviolet light from the top of the photosphere, where the atmosphere is "
          "coolest (about 4,500 K). Bright patches (plage) mark magnetic field around active "
          "regions; the fine mesh is the magnetic network at the edges of convection cells "
          "some 30,000 km across.",
    1600: "Ultraviolet from the upper photosphere plus carbon (C IV) emission from 100,000 K "
          "gas. It looks much like 1700 Å but flares light it up first, where their energy "
          "slams into the lower atmosphere.",
    304: "Helium (He II) at about 50,000 K: the chromosphere and transition region, a layer "
         "only a few thousand km thick between the surface and the corona. Prominences, cool "
         "dense gas held up by the field, arch above the limb; seen on the disc they are dark "
         "filaments.",
    131: "Iron (Fe VIII) at 0.4 MK on the quiet Sun, but Fe XXI at 10 MK wherever a flare "
         "heats gas, so this is the channel that sees flares first and hottest.",
    171: "Iron (Fe IX) at 0.6 MK: the classic picture of coronal loops, gas trapped on "
         "magnetic field lines that arch from one polarity to the other.",
    193: "Iron (Fe XII) at 1.6 MK: the corona proper. The dark patch {hole} is a coronal "
         "hole, where field lines open into space and the fast solar wind escapes; it is dark "
         "because the gas there is thin, not because it is cold.",
    211: "Iron (Fe XIV) at 2 MK: the hotter corona above active regions. The coronal hole is "
         "darker still at this temperature.",
    335: "Iron (Fe XVI) at 2.5 MK: the hottest long-lived gas above active regions. On the "
         "quiet Sun cooler lines of magnesium, aluminium and silicon supply much of this "
         "channel's faint signal.",
    94: "Iron (Fe XVIII) at 6 MK: flaring loops and the hottest active-region cores. The "
        "faint quiet-Sun glow in this channel is mostly Fe X at 1 MK (O'Dwyer et al. 2010).",
    "rgb": "211, 193 and 171 Å as red, green and blue: 2, 1.6 and 0.6 MK gas in one picture. "
           "Cooler loops come out blue, hotter active regions yellow-red, and the coronal hole "
           "dark in all three.",
}


def _jpeg(rgb: np.ndarray, px: int, quality: int) -> str:
    img = Image.fromarray((np.clip(rgb, 0, 1) * 255).astype(np.uint8))
    img = img.resize((px, px), Image.LANCZOS)
    buf = io.BytesIO()
    img.save(buf, format="JPEG", quality=quality, optimize=True, progressive=True)
    return "data:image/jpeg;base64," + base64.b64encode(buf.getvalue()).decode()


def _crop(rgb: np.ndarray, m, lim: float = LIM_ARCSEC, fill: float = 0.0) -> np.ndarray:
    """The +-lim arcsec square around disc centre, north up, row 0 at the top.

    Images with finer pixels (HMI's 2.0" against AIA's 2.4") do not reach that far,
    so the square is a canvas of sky colour with the image copied into it; slicing
    past an array's edge would silently wrap around instead.
    """
    from astropy.coordinates import SkyCoord

    c = m.world_to_pixel(SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=m.coordinate_frame))
    half = int(round(lim / m.scale.axis1.to_value(u.arcsec / u.pix)))
    x0, y0 = int(round(c.x.value)), int(round(c.y.value))
    ny, nx = rgb.shape[:2]
    out = np.full((2 * half, 2 * half, 3), fill, dtype=float)
    ya, yb = max(0, y0 - half), min(ny, y0 + half)
    xa, xb = max(0, x0 - half), min(nx, x0 + half)
    out[ya - (y0 - half):yb - (y0 - half), xa - (x0 - half):xb - (x0 - half)] = rgb[ya:yb, xa:xb]
    return out[::-1]


def _render_channel(m, key) -> np.ndarray:
    """RGB image of one view with the same scaling as figure 1, north up."""
    d = np.asarray(m.data, dtype=float)
    if key == "hmi":
        norm = figures.ImageNormalize(vmin=-250, vmax=250, clip=True)
        rgb = plt.get_cmap("gray")(np.clip(norm(np.nan_to_num(d)), 0, 1))[..., :3]
        rgb[~np.isfinite(d)] = 0.0
    else:
        norm = figures.channel_norm(m, key)
        rgb = plt.get_cmap(C.CHANNEL_BY_WAVELENGTH[key].cmap)(np.clip(norm(d), 0, 1))[..., :3]
    return _crop(rgb, m)


def _coronal_hole_phrase(maps) -> str:
    """Where the largest dark 193 A region sits, in words (for the 193 note)."""
    from scipy import ndimage

    m = maps[193]
    d = np.asarray(m.data, dtype=float)
    r = limb.fast_radius(m)
    on = r < 0.9
    smooth = ndimage.gaussian_filter(d, 6)
    dark = on & (smooth < 0.35 * np.median(d[on]))
    lab, n = ndimage.label(dark)
    if n == 0:
        return "(if any)"
    sizes = ndimage.sum(dark, lab, range(1, n + 1))
    yy, xx = np.nonzero(lab == (int(np.argmax(sizes)) + 1))
    cy, cx = np.mean(yy), np.mean(xx)
    c = m.pixel_to_world(cx * u.pix, cy * u.pix)
    x, y = c.Tx.to_value(u.arcsec), c.Ty.to_value(u.arcsec)
    ew = "east" if x < -150 else ("west" if x > 150 else "")
    ns = "north" if y > 150 else ("south" if y < -150 else "")
    where = "-".join(w for w in (ns, ew) if w)
    return f"{where} of the disc centre" if where else "near the disc centre"


def views(products, results) -> list[dict]:
    maps, hmi = products["maps"], products["hmi"]
    out = []
    order = ["hmi", 4500, 1700, 1600, 304, 131, 171, 193, 211, 335, 94, "rgb"]
    hole = _coronal_hole_phrase(maps)
    for key in order:
        if key == "hmi":
            rgb, title, sub = _render_channel(hmi, "hmi"), "Magnetic field", "HMI · surface"
            log_t = None
        elif key == "rgb":
            rgb = _crop(figures.rgb_composite(maps), maps[171])
            title, sub, log_t = "Composite", "211 + 193 + 171 Å", None
        else:
            rgb = _render_channel(maps[key], key)
            ch = C.CHANNEL_BY_WAVELENGTH[key]
            title = f"{key} Å"
            sub = f"{ch.ions} · " + " & ".join(figures._t_label(10 ** t) for t in ch.log_t)
            log_t = ch.log_t[0]
        out.append({"key": str(key), "title": title, "sub": sub, "logT": log_t,
                    "note": NOTES[key].format(hole=hole),
                    "src": _jpeg(rgb, VIEW_PX, 86)})
    return out


def week(products, results) -> dict:
    paths = products["epoch_paths"]
    maps = movie.week_frames(paths)
    norm = movie.common_norm(maps)
    cmap = plt.get_cmap("sdoaia193")
    frames, times = [], []
    for m in maps:
        rgb = cmap(np.clip(norm(np.asarray(m.data, dtype=float)), 0, 1))[..., :3]
        frames.append(_jpeg(_crop(rgb, m), WEEK_PX, 78))
        times.append(pd.Timestamp(m.date.isot, tz="UTC"))
    xrs = products["xrs"]["long"].resample("10min").mean()
    t0, t1 = xrs.index[0], xrs.index[-1]
    span = (t1 - t0).total_seconds()
    # the X-ray curve as SVG path data in a 1000 x 160 box, log flux 1e-7..1e-5
    pts, path, pen = [], [], "M"
    for t, v in xrs.items():
        x = 1000 * (t - t0).total_seconds() / span
        if not (np.isfinite(v) and v > 0):
            pen = "M"
            continue
        y = 160 * (1 - (np.log10(v) + 7) / 2)
        path.append(f"{pen}{x:.1f},{y:.1f}")
        pen = "L"
        pts.append([round(x, 1), round(float(v), 10), t.strftime("%-d %b %H:%M")])
    flare = results["flare"]["selection"]
    fx = 1000 * (pd.Timestamp(flare["peak"]) - t0).total_seconds() / span if flare else None
    return {
        "frames": frames,
        "times": [t.strftime("%a %-d %b %Y, %H:%M UTC") for t in times],
        "x": [round(1000 * (t - t0).total_seconds() / span, 2) for t in times],
        "path": " ".join(path),
        "points": pts[::3],
        "flare_x": fx,
        "days": [{"x": round(1000 * (d - t0).total_seconds() / span, 1), "label": d.strftime("%-d %b")}
                 for d in pd.date_range(t0.ceil("1D"), t1, freq="1D")],
    }


PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>The Sun Today</title>
<style>
:root { color-scheme: dark; --page:#0d0d0d; --surface:#1a1a19; --ink:#ffffff; --ink2:#c3c2b7;
  --muted:#898781; --line:#383835; --accent:#3987e5; }
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--ink);
  font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1160px; margin: 0 auto; padding: 24px 16px 48px; }
h1 { font-size: 26px; margin: 0 0 4px; }
h2 { font-size: 19px; margin: 40px 0 6px; }
p.lead { color: var(--ink2); margin: 0 0 18px; max-width: 900px; }
.views { display: grid; grid-template-columns: minmax(0, 720px) minmax(260px, 1fr); gap: 20px; }
.frame { width: 100%; aspect-ratio: 1; background: #000; border-radius: 6px; overflow: hidden; }
.frame img { width: 100%; height: 100%; display: block; }
.ladder { display: flex; flex-direction: column; gap: 4px; }
.ladder button { all: unset; cursor: pointer; display: grid; grid-template-columns: 92px 1fr;
  gap: 8px; padding: 6px 10px; border-radius: 6px; color: var(--ink2); background: var(--surface);
  border: 1px solid transparent; }
.ladder button:hover { border-color: var(--line); color: var(--ink); }
.ladder button[aria-pressed="true"] { border-color: var(--accent); color: var(--ink); }
.ladder button:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; }
.ladder b { font-weight: 600; color: var(--ink); }
.ladder .t { font-size: 13px; }
.note { margin-top: 12px; color: var(--ink2); min-height: 6em; }
.note b { color: var(--ink); }
.hint { color: var(--muted); font-size: 13px; }
.week { display: grid; grid-template-columns: minmax(0, 440px) 1fr; gap: 20px; align-items: start; }
.controls { display: flex; gap: 10px; align-items: center; margin: 10px 0; }
.controls button { all: unset; cursor: pointer; padding: 5px 14px; border-radius: 6px;
  background: var(--surface); border: 1px solid var(--line); color: var(--ink); }
.controls button:focus-visible { outline: 2px solid var(--accent); }
input[type=range] { width: 100%; accent-color: var(--accent); }
.stamp { font-weight: 600; }
svg.xr { width: 100%; height: auto; background: var(--surface); border-radius: 6px; cursor: crosshair; }
svg.xr text { fill: var(--muted); font-size: 22px; }
.tip { position: fixed; pointer-events: none; background: #000; color: var(--ink);
  border: 1px solid var(--line); border-radius: 4px; padding: 3px 7px; font-size: 13px; display: none; }
footer { margin-top: 40px; color: var(--muted); font-size: 13px; }
@media (max-width: 820px) {
  .views, .week { grid-template-columns: minmax(0, 1fr); }
}
</style>
</head>
<body>
<main>
<h1>The Sun on __DATE__</h1>
<p class="lead">Real images from NASA's Solar Dynamics Observatory, taken at __TIME__ UTC. Each
view shows gas at a different temperature, from the 5,800 K surface to 6-million-degree
flare loops. Pick one from the ladder (or use the ← → keys).</p>
<section class="views">
  <div>
    <div class="frame"><img id="view" alt=""></div>
    <div class="hint">Solar north is up; the Sun rotates from left (east) to right (west).</div>
  </div>
  <div>
    <div class="ladder" id="ladder" role="group" aria-label="Views, coolest to hottest"></div>
    <div class="note" id="note"></div>
  </div>
</section>

<h2>The week before, in 193 Å</h2>
<p class="lead">One image every 3 hours for 7 days. Watch active regions ride across the disc
as the Sun turns (once every 27 days as seen from Earth). Below, the X-ray output measured by
GOES: every spike is a flare. Click the curve to jump to that moment.</p>
<section class="week">
  <div>
    <div class="frame"><img id="wk" alt="The Sun in 193 angstroms"></div>
    <div class="controls">
      <button id="play" aria-label="Play or pause">Play</button>
      <span class="stamp" id="stamp"></span>
    </div>
    <input type="range" id="slider" min="0" value="0" aria-label="Time">
  </div>
  <div>
    <svg class="xr" id="xr" viewBox="-60 -20 1080 230" role="img"
         aria-label="GOES X-ray flux over the week">
      <line x1="0" x2="1000" y1="0" y2="0" stroke="#383835"/>
      <line x1="0" x2="1000" y1="80" y2="80" stroke="#383835"/>
      <line x1="0" x2="1000" y1="160" y2="160" stroke="#383835"/>
      <text x="-50" y="8">M</text><text x="-50" y="88">C</text><text x="-50" y="168">B</text>
      <g id="days"></g>
      <line id="flare" y1="-10" y2="170" stroke="#898781" stroke-dasharray="4 4"/>
      <path id="curve" fill="none" stroke="#3987e5" stroke-width="2"/>
      <line id="cursor" y1="-10" y2="170" stroke="#ffffff" stroke-width="2.5"/>
    </svg>
    <p class="hint">GOES-18, 1–8 Å, 10-minute means. Dotted line: the week's most energetic
    flare, __FLARE__. Letters are flare classes (each ten times the one below).</p>
  </div>
</section>
<footer>Data: SDO/AIA and SDO/HMI (NASA; JSOC, Stanford), GOES-18 XRS (NOAA SWPC and NCEI).
Built by the Astrophysics Project Atlas, project F2.</footer>
</main>
<div class="tip" id="tip"></div>
<script>
const VIEWS = __VIEWS__;
const WEEK = __WEEK__;
const view = document.getElementById("view"), note = document.getElementById("note"),
      ladder = document.getElementById("ladder");
let current = 0;
function show(i) {
  current = (i + VIEWS.length) % VIEWS.length;
  const v = VIEWS[current];
  view.src = v.src;
  view.alt = v.title + ": " + v.sub;
  note.innerHTML = "<b>" + v.title + "</b> · " + v.sub + "<br>" + v.note;
  ladder.querySelectorAll("button").forEach((b, k) => b.setAttribute("aria-pressed", k === current));
}
VIEWS.forEach((v, k) => {
  const b = document.createElement("button");
  b.innerHTML = "<b>" + v.title + "</b><span class='t'>" + v.sub + "</span>";
  b.addEventListener("click", () => show(k));
  ladder.appendChild(b);
});
document.addEventListener("keydown", e => {
  if (e.target.tagName === "INPUT") return;
  let step = 0;
  if (e.key === "ArrowRight" || e.key === "ArrowDown") step = 1;
  if (e.key === "ArrowLeft" || e.key === "ArrowUp") step = -1;
  if (!step) return;
  show(current + step);
  if (ladder.contains(document.activeElement)) ladder.children[current].focus();
  e.preventDefault();
});
show(0);

const wk = document.getElementById("wk"), slider = document.getElementById("slider"),
      stamp = document.getElementById("stamp"), cursor = document.getElementById("cursor"),
      play = document.getElementById("play"), svg = document.getElementById("xr"),
      tip = document.getElementById("tip");
slider.max = WEEK.frames.length - 1;
document.getElementById("curve").setAttribute("d", WEEK.path);
const flare = document.getElementById("flare");
if (WEEK.flare_x === null) flare.remove(); else { flare.setAttribute("x1", WEEK.flare_x); flare.setAttribute("x2", WEEK.flare_x); }
const days = document.getElementById("days");
WEEK.days.forEach(d => {
  days.insertAdjacentHTML("beforeend",
    `<line x1="${d.x}" x2="${d.x}" y1="160" y2="172" stroke="#898781"/>` +
    `<text x="${d.x}" y="198" text-anchor="middle">${d.label}</text>`);
});
function frame(k) {
  k = Math.max(0, Math.min(WEEK.frames.length - 1, k));
  slider.value = k;
  wk.src = WEEK.frames[k];
  stamp.textContent = WEEK.times[k];
  cursor.setAttribute("x1", WEEK.x[k]);
  cursor.setAttribute("x2", WEEK.x[k]);
}
slider.addEventListener("input", () => { stop(); frame(+slider.value); });
let timer = null;
function stop() { if (timer) { clearInterval(timer); timer = null; play.textContent = "Play"; } }
play.addEventListener("click", () => {
  if (timer) { stop(); return; }
  play.textContent = "Pause";
  timer = setInterval(() => {
    const k = (+slider.value + 1) % WEEK.frames.length;
    frame(k);
  }, 125);
});
function svgX(evt) {
  const pt = svg.createSVGPoint(); pt.x = evt.clientX; pt.y = evt.clientY;
  return pt.matrixTransform(svg.getScreenCTM().inverse()).x;
}
function nearest(arr, x) {
  let best = 0;
  for (let k = 1; k < arr.length; k++) if (Math.abs(arr[k] - x) < Math.abs(arr[best] - x)) best = k;
  return best;
}
svg.addEventListener("click", e => { stop(); frame(nearest(WEEK.x, svgX(e))); });
function cls(f) {
  const L = [["A", 1e-8], ["B", 1e-7], ["C", 1e-6], ["M", 1e-5], ["X", 1e-4]];
  let c = L[0]; for (const l of L) if (f >= l[1]) c = l;
  return c[0] + (Math.floor(f / c[1] * 10 + 1e-6) / 10).toFixed(1);
}
svg.addEventListener("mousemove", e => {
  const x = svgX(e);
  const k = nearest(WEEK.points.map(p => p[0]), x), p = WEEK.points[k];
  tip.style.display = "block";
  tip.style.left = (e.clientX + 12) + "px"; tip.style.top = (e.clientY - 28) + "px";
  tip.textContent = p[2] + " UTC · " + cls(p[1]) + " level";
});
svg.addEventListener("mouseleave", () => { tip.style.display = "none"; });
frame(WEEK.frames.length - 1);
</script>
</body>
</html>
"""


def build(products, results, out_dir) -> str:
    epoch = datetime.fromisoformat(results["epoch"])
    sel = results["flare"]["selection"]
    reg = (results["flare"].get("location") or {}).get("noaa") or {}
    flare_txt = (f"{sel['class']} from NOAA {reg.get('region', '?')} at "
                 f"{pd.Timestamp(sel['peak']):%-d %b %H:%M} UTC") if sel else "none"
    html = (PAGE.replace("__DATE__", f"{epoch:%-d %B %Y}")
                .replace("__TIME__", f"{epoch:%H:%M}")
                .replace("__FLARE__", flare_txt)
                .replace("__VIEWS__", json.dumps(views(products, results)))
                .replace("__WEEK__", json.dumps(week(products, results))))
    out = out_dir / "sun_today.html"
    out.write_text(html)
    return str(out)
