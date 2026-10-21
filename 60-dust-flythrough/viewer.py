"""Build the interactive 3D fly-through (WebGL 2, no libraries) from the cached grids.

Outputs, both opened by double-click or served as they are:

* ``figures/dust_flythrough.html``: one self-contained file (data embedded,
  ~28 MB), for offline use;
* ``figures/dust_viewer/``: the same page with the volumes as separate files
  (``index.html`` + base64 text parts of at most 10 MB, since the publishing
  host serves text but not raw binary), which is what gets published online.

The dust is stored as three nested volumes of 8-bit codes on a log scale
(code 0 = below 2e-5 ZGR23 mag/pc, i.e. n_H < 0.03 cm^-3, treated as empty):

* inner: the main map on 2.5 pc cells, |x|, |y| < 400 pc, |z| < 250 pc;
* middle: the main map on 5 pc cells, |x|, |y| < 1250 pc, |z| < 400 pc, used
  within 1240 pc of the Sun;
* outer: the 2 kpc validation map on 10 pc cells, |x|, |y| < 2000 pc,
  |z| < 500 pc, used beyond (the authors call it less reliable).

The GPU's trilinear filtering of the log codes is a geometric interpolation,
the same convention as the map's own tools.
"""

from __future__ import annotations

import base64
import json
import shutil
import zlib

import healpy as hp
import numpy as np

import catalogs
import config as C
import grids
from deredden import coefficient

LO, HI = np.log10(2e-5), np.log10(0.2)
TEMPLATE = C.PROJECT_DIR / "viewer_template.html"
OUT_HTML = C.FIGURE_DIR / "dust_flythrough.html"
OUT_DIR = C.FIGURE_DIR / "dust_viewer"


def quantise(a: np.ndarray) -> np.ndarray:
    """Density (ZGR23 mag/pc) -> uint8 codes, 0 = below 10**LO."""
    with np.errstate(divide="ignore", invalid="ignore"):
        lg = np.log10(np.where(a > 0, a, 1e-30))
    code = np.round((np.clip(lg, LO, HI) - LO) / (HI - LO) * 254 + 1)
    code[~(a >= 10**LO)] = 0
    return code.astype(np.uint8)


def levels() -> tuple[list, dict]:
    g = grids.build_main()
    w = grids.build_wide()
    inner_xs, inner_zs = g["inner_xs"], g["inner_zs"]
    ii = np.flatnonzero(np.abs(inner_xs) < 400)
    kk = np.flatnonzero(np.abs(inner_zs) < 250)
    inner = g["inner"][np.ix_(ii, ii, kk)]
    v5zs = g["vol5_zs"]
    km = np.flatnonzero(np.abs(v5zs) < 400)
    middle = g["vol5"][:, :, km]
    outer = w["vol10"]
    specs = [("inner", inner, 2.5, (-400, -400, -250), 1240.0, 1.0),
             ("middle", middle, 5.0, (-1250, -1250, -400), 1240.0, 2.0),
             ("outer", outer, 10.0, (-2000, -2000, -500), 1985.0, 4.0)]
    meta, blobs = [], {}
    for name, arr, cell, lo_corner, rmax, step in specs:
        codes = quantise(arr)
        raw = np.ascontiguousarray(np.transpose(codes)).tobytes()     # x fastest, then y, then z
        blobs[name] = zlib.compress(raw, 9)
        shape = list(arr.shape)
        meta.append({"name": name, "shape": shape, "cell": cell, "min": list(lo_corner),
                     "max": [lo_corner[i] + cell * shape[i] for i in range(3)], "rmax": rmax, "step": step})
    return meta, blobs


def overlays(wave_curve: np.ndarray | None = None) -> dict:
    hr = catalogs.hunt_reffert()
    young = hr[hr["reliable"] & (hr["age_myr"] < 30) & (hr["dist50"] < 2000)]
    clusters = [[round(r.x, 1), round(r.y, 1), round(r.z, 1), round(r.age_myr, 1), int(r.N)]
                for r in young.itertuples()]
    if wave_curve is None:
        import wave
        x, y, z, _ = wave.kon_curve(wave.KON_PUBLISHED["clouds"], 400)
        wave_curve = np.c_[x, y, z]
    import bubble
    shell = bubble.shell_model()
    d = shell["d_peak"].to_numpy().reshape(-1, 64)        # NESTED: 64 children of each Nside-32 pixel
    d32 = np.nanmedian(np.where(np.isfinite(d), d, np.nan), axis=1)
    l32, b32 = hp.pix2ang(32, np.arange(hp.nside2npix(32)), nest=True, lonlat=True)
    ok = np.isfinite(d32)
    from dustmap import lbd_to_xyz
    sx, sy, sz = lbd_to_xyz(l32[ok], b32[ok], d32[ok])
    shell_pts = np.round(np.c_[sx, sy, sz], 1).tolist()
    labs = catalogs.cloud_labels()
    keep = {"Taurus": 600, "Perseus": 900, "California": 1100, "Orion": 1500, "Mon R2": 2200, "Mon OB1": 2200,
            "CMa OB1": 3500, "Vela": 2500, "Chamaeleon": 900, "Lupus": 800, "Ophiuchus": 800,
            "Corona Australis": 800, "Aquila Rift": 1200, "Serpens": 1800, "Cepheus": 1500,
            "North America": 2500, "Cygnus X": 3500, "Coalsack": 900, "IC5146": 2500, "IC1396": 2500,
            "Rosette": 3000}
    labels = [{"name": r.complex, "pos": [round(r.x, 1), round(r.y, 1), round(r.z, 1)], "range": keep[r.complex]}
              for r in labs.itertuples() if r.complex in keep]
    labels.append({"name": "Radcliffe Wave", "pos": [float(v) for v in wave_curve[int(0.62 * len(wave_curve))]],
                   "cls": "wave big", "range": 6000})
    labels.append({"name": "Local Bubble wall", "pos": [-30.0, -150.0, 160.0], "cls": "shell", "range": 2500})
    labels.append({"name": "Sun", "pos": [0.0, 0.0, 0.0], "cls": "big", "range": 5000})
    labels.append({"name": "to the Galactic Centre", "pos": [1900.0, 0.0, 0.0], "range": 9000})
    return {"clusters": clusters, "wave": np.round(wave_curve, 1).tolist(), "shell": shell_pts,
            "rings": [250, 500, 1000, 1500, 2000], "labels": labels}


def views(curve: np.ndarray) -> dict:
    geo = wave_geometry(curve)
    mid, n = geo["mid"], geo["n"]
    return {
        "sun": {"pos": [0.0, 0.0, 0.0], "target": [-141.0, 19.0, -36.0]},
        "taurus": {"pos": [-40.0, -70.0, 30.0], "target": [-150.0, 20.0, -40.0]},
        "orion": {"pos": [-150.0, -260.0, -20.0], "target": [-380.0, -170.0, -120.0]},
        "above": {"pos": [350.0, -1150.0, 1900.0], "target": [-150.0, 100.0, 0.0]},
        "wave": {"pos": [float(mid[0] - 2600 * n[0]), float(mid[1] - 2600 * n[1]), 120.0],
                 "target": [float(mid[0]), float(mid[1]), 0.0]},
    }


def wave_geometry(curve: np.ndarray) -> dict:
    """Chord, in-plane normal (pointing away from the anticentre side) and midpoint of the wave model."""
    a, b = curve[0, :2], curve[-1, :2]
    d = (b - a) / np.linalg.norm(b - a)
    n = np.array([d[1], -d[0]])
    if n[0] < 0:
        n = -n
    mid = 0.5 * (a + b)
    return {"dir": d, "n": n, "mid": mid}


def tour(curve: np.ndarray) -> dict:
    """Keyframes (camera position and look-at point, pc) and captions, in seconds."""
    geo = wave_geometry(curve)
    mid, n, d = geo["mid"], geo["n"], geo["dir"]
    m3 = [float(mid[0]), float(mid[1]), 0.0]
    side = lambda dist, along=0.0, z=60.0: [float(mid[0] - dist * n[0] + along * d[0]),
                                           float(mid[1] - dist * n[1] + along * d[1]), z]
    s0 = float(mid @ n)
    a0, a1 = float(curve[0, :2] @ d), float(curve[-1, :2] @ d)
    slab = {"n": [float(n[0]), float(n[1]), 0.0], "lo": s0 - 350.0, "hi": s0 + 350.0,
            "d": [float(d[0]), float(d[1]), 0.0], "lo2": min(a0, a1) - 250.0, "hi2": max(a0, a1) + 250.0}
    above = [float(mid[0] + 700), float(mid[1] - 1300), 2600.0]
    k = [
        (0, [0, 0, 0], [-141, 19, -36], {"clusters": True, "wave": False, "shell": False, "labels": True, "grid": False}, None),
        (5, [-25, 5, -6], [-141, 19, -36], None, None),
        (13, [-95, -25, -20], [-160, 30, -45], None, None),
        (19, [-150, -40, -40], [-255, 95, -90], None, None),
        (25, [-215, -10, -70], [-300, 110, -95], None, None),
        (31, [-250, -60, -80], [-375, -175, -120], None, None),
        (38, [-300, -140, -95], [-400, -180, -125], None, None),
        (43, [-180, -330, 60], [-330, -170, -100], None, None),
        (50, [150, -650, 520], [-80, 0, -10], {"clusters": True, "wave": False, "shell": True, "labels": True, "grid": True}, None),
        (57, [420, -560, 380], [-40, 0, 0], None, None),
        (64, above, m3, {"clusters": True, "wave": True, "shell": False, "labels": True, "grid": True}, None),
        (71, [above[0] - 900, above[1] - 300, 2300.0], m3, None, None),
        (78, side(2600, 0, 400), m3, None, slab),
        (84, side(2600, 0, 60), m3, None, None),
        (91, side(2500, 350, 90), [m3[0] + 250 * d[0], m3[1] + 250 * d[1], 0.0], None, None),
        (96, [700, -700, 350], [0, 0, 0], {"clusters": True, "wave": True, "shell": False, "labels": True, "grid": True}, "off"),
        (100, [350, -350, 180], [0, 0, 0], None, None),
    ]
    keys = []
    for t, pos, target, layers, clip in k:
        key = {"t": t, "pos": [float(v) for v in pos], "target": [float(v) for v in target]}
        if layers:
            key["layers"] = layers
        if clip == "off":
            key["clip"] = None
        elif clip is not None:
            key["clip"] = clip
        keys.append(key)
    captions = [
        {"t0": 0.3, "t1": 6.5, "text": "We start at the Sun.",
         "note": "The map has no data within 69 pc: we sit in the Local Bubble, a cavity carved by supernovae."},
        {"t0": 7, "t1": 17, "text": "Taurus, 140 pc away",
         "note": "The nearest large star-forming cloud. Blue dots are star clusters younger than 30 million years."},
        {"t0": 18, "t1": 27, "text": "Perseus, about 290 pc", "note": "Young clusters IC 348 and NGC 1333 sit in its dust."},
        {"t0": 28.5, "t1": 41, "text": "Orion, about 400 pc",
         "note": "The nearest giant star factory: Orion A, Orion B and the lambda Orionis ring."},
        {"t0": 45, "t1": 61, "text": "The Local Bubble's dusty wall",
         "note": "Green dots: where dust first piles up along each line of sight from the Sun. Many young clusters near the Sun sit on it."},
        {"t0": 63.5, "t1": 76, "text": "The Radcliffe Wave",
         "note": "A 2.7 kpc chain of clouds from Canis Major to Cygnus. Orange: the wave model fitted to cloud distances."},
        {"t0": 78, "t1": 93, "text": "Seen from the side, it rises and falls",
         "note": "Cutaway: only dust within 350 pc of the Wave's axis is drawn. Its clouds swing up to ~200 pc above and below the plane."},
        {"t0": 94.5, "t1": 100, "text": "Back toward the Sun", "note": "The detailed map reaches 1.25 kpc; the coarser one beyond it, 2 kpc."},
    ]
    return {"duration": 100.0, "keys": keys, "captions": captions}


def build(wave_curve: np.ndarray | None = None) -> dict:
    lv, blobs = levels()
    ov = overlays(wave_curve)
    curve = np.asarray(ov["wave"])
    texts = {name: base64.b64encode(blob).decode() for name, blob in blobs.items()}
    part = 10_000_000
    files = {name: [f"dust_{name}_{k}.txt" for k in range(-(-len(txt) // part))] for name, txt in texts.items()}
    meta = {"files": files, "levels": lv, "rmin": float(C.R_INNER), "lo": float(LO), "hi": float(HI),
            "kappa_rgb": [coefficient(650), coefficient(550), coefficient(450)], "solid_opacity": 4.0,
            "solid_fade": [6e-5, 5e-4],
            "zgr23_to_av": C.ZGR23_TO_AV, "glow": [0.02, 12.0], "overlays": ov, "views": views(curve),
            "tour": tour(curve)}
    template = TEMPLATE.read_text()
    page = template.replace("/*__META__*/", json.dumps(meta, separators=(",", ":")))
    # online build: page + the volumes as base64 text parts
    if OUT_DIR.exists():
        shutil.rmtree(OUT_DIR)
    OUT_DIR.mkdir()
    (OUT_DIR / "index.html").write_text(page.replace("<!--__DATA__-->", ""), encoding="utf-8")
    for name, txt in texts.items():
        for k, fname in enumerate(files[name]):
            (OUT_DIR / fname).write_text(txt[k * part:(k + 1) * part])
    # offline build: everything in one file
    data = "\n".join(f'<script type="application/octet-stream" id="data-{name}">{txt}</script>'
                     for name, txt in texts.items())
    # a standalone file needs its own document head (the publishing skeleton provides one online)
    head = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n')
    body = page.replace("<!--__DATA__-->", data)
    t_end = body.index("</style>") + len("</style>")
    OUT_HTML.write_text(head + body[:t_end] + "\n</head>\n<body>\n" + body[t_end:] + "\n</body>\n</html>\n",
                        encoding="utf-8")
    sizes = {name: len(blob) for name, blob in blobs.items()}
    return {"html_bytes": OUT_HTML.stat().st_size, "volume_bytes": sizes,
            "index_bytes": (OUT_DIR / "index.html").stat().st_size}


if __name__ == "__main__":
    print(json.dumps(build(), indent=1))
