"""The interactive "build a universe" page, written from the analysis caches.

    python explorer.py

writes ``figures/build_a_universe.html``: one self-contained file (no network
needed; web fonts are optional and fall back to system fonts) in which the
Friedmann equation is solved live in JavaScript, with the same method as
``friedmann.py``. The page carries our credible regions for supernovae, BAO,
the CMB and all three combined, a coarse grid of the supernova chi^2 for the
"what the data say" readout, the 32 cosmic chronometers and DESI's
line-of-sight BAO points for its speedometer. It also writes the same page
without the document skeleton to ``data/derived/`` (the form the atlas
publishes as a web page).
"""

from __future__ import annotations

import json
from datetime import date

import contourpy
import numpy as np

import config as cfg
import fetch_data as fd
import friedmann as fr
import likelihoods as lk

TEMPLATE = cfg.PROJECT_DIR / "explorer_template.html"
OUT = cfg.FIGURE_DIR / "build_a_universe.html"
FRAGMENT = cfg.DERIVED_DIR / "build_a_universe_fragment.html"


def _load(name: str):
    with np.load(cfg.DERIVED_DIR / f"{name}.npz", allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def contours(X, Y, p, max_pts: int = 180) -> dict:
    """68% and 95% highest-density polygons of a grid posterior, in (Om, OL)."""
    t68, t95, _ = lk.credible_levels(p)
    gen = contourpy.contour_generator(np.asarray(X, float), np.asarray(Y, float),
                                      np.asarray(p, float))
    out = {}
    for name, level in (("p95", t95), ("p68", t68)):
        polys = []
        for line in gen.lines(level):
            if len(line) < 4:
                continue
            step = max(1, len(line) // max_pts)
            pts = line[::step]
            if not np.allclose(pts[-1], line[-1]):
                pts = np.vstack([pts, line[-1]])
            polys.append(np.round(pts, 4).tolist())
        out[name] = polys
    return out


def build_data() -> dict:
    res = json.loads(cfg.RESULTS_PATH.read_text())
    sn, bao, cmb = _load("sn"), _load("bao"), _load("cmb_band_pr4")
    allc = _load("combo_sn+bao+cmb")
    OMs, OLs = np.meshgrid(sn["Om"], sn["OL"], indexing="ij")
    OMb = np.repeat(bao["Om"][:, None], bao["OL"].shape[1], axis=1)
    OMc = np.repeat(cmb["Om"][:, None], cmb["OL"].shape[1], axis=1)
    OMa = np.repeat(allc["Om"][:, None], allc["Ok"].size, axis=1)
    data = {"contours": {"sn": contours(OMs, OLs, sn["p"]),
                         "bao": contours(OMb, bao["OL"], bao["p"]),
                         "cmb": contours(OMc, cmb["OL"], cmb["p"]),
                         "all": contours(OMa, allc["OL"], allc["p"])}}
    # supernova delta-chi^2 on a coarse grid (every 4th point of D1's grid)
    d = (sn["chi2"] - sn["chi2"].min())[::4, ::4]
    data["sngrid"] = {"om0": float(sn["Om"][0]), "dom": float(sn["Om"][4] - sn["Om"][0]),
                      "nom": int(d.shape[0]), "ol0": float(sn["OL"][0]),
                      "dol": float(sn["OL"][4] - sn["OL"][0]), "nol": int(d.shape[1]),
                      "d": np.round(np.minimum(d, 999.0), 1).ravel().tolist()}
    cc = fd.load_chronometers()
    data["cc"] = [[float(z), float(h), float(e)] for z, h, e in zip(cc["z"], cc["H"], cc["err"])]
    b = fd.load_bao()
    rd = res["ours"]["epochs"]["r_drag"]
    dh = b["kind"] == "DH_over_rs"
    data["desi"] = [[float(z), float(cfg.C_KM_S / (v * rd)),
                     float(cfg.C_KM_S / (v * rd) * np.sqrt(c) / v)]
                    for z, v, c in zip(b["z"][dh], b["value"][dh], np.diag(b["cov"])[dh])]
    p = res["ours"]["params"]
    data["ours"] = {"Om": round(p["Om"], 4), "OL": round(p["OL"], 4), "Or": p["Or"],
                    "H0": round(p["H0"], 2)}
    # presets: the gallery, in its order, plus a universe balanced on the edge
    presets = []
    Or_cmb = fr.omega_r_today(p["H0"])
    for u in res["gallery"]["universes"]:
        Or = "cmb" if abs(u["Or"] - Or_cmb) < 1e-9 else float(u["Or"])
        name = {"ours": "Our universe"}.get(u["key"], u["name"])
        presets.append({"key": u["key"], "name": name, "Om": round(u["Om"], 4),
                        "OL": round(u["OL"], 4), "Or": Or, "w": u["w"], "story": u["story"]})
    edge = float(fr.no_big_bang_boundary([0.3], Or_cmb)[0])
    presets.insert(10, {"key": "edge", "name": "On the edge (Eddington–Lemaître)", "Om": 0.3,
                        "OL": round(edge - 2e-4, 4), "Or": "cmb", "w": -1.0,
                        "story": "just below the no-Big-Bang edge: it lingers for tens of billions of "
                                 "years near Einstein's static size before expanding again"})
    data["presets"] = presets
    data["hubbleTimeGyrTimesH0"] = cfg.HUBBLE_TIME_GYR_X_H0
    data["omegaR_h2"] = fr.OMEGA_R_TOY
    data["credits"] = (f"Built on {date.today().isoformat()} from atlas project O1 "
                       "(65-friedmann-universes). Everything here is AI-generated; the "
                       "numbers are compared with published values in the project's README.")
    return data


def main() -> int:
    data = build_data()
    blob = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    tpl = TEMPLATE.read_text()
    page = tpl.replace("/*__DATA__*/null", blob)
    if page == tpl:
        raise RuntimeError("data placeholder not found in the template")
    head, body = page.split("<!--BODY-->")
    FRAGMENT.write_text(head + body)
    OUT.write_text("<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
                   "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, "
                   "viewport-fit=cover\">\n" + head +
                   "<style>body{margin:0}[hidden]{display:none!important}"
                   ":root{padding-top:env(safe-area-inset-top,0px);"
                   "padding-bottom:env(safe-area-inset-bottom,0px)}</style>\n</head>\n<body>\n"
                   + body + "\n</body>\n</html>\n")
    print(f"  wrote {OUT.relative_to(cfg.PROJECT_DIR)} ({OUT.stat().st_size / 1e3:.0f} kB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
