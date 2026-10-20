"""Build the interactive inflation explorer: one self-contained HTML file, no libraries.

    python explorer.py      # writes figures/inflation_explorer.html (~0.5 MB)

It also writes data/derived/inflation_explorer_artifact.html, the same page as a
fragment for a viewer that supplies the document skeleton and the light/dark
theme (so it has no theme buttons of its own).

The page carries, from the cached analysis: the 68% / 95% regions of the four
data planes as polygons (linear and log-r views), each plane's credibility on a
grid (so any point in the plane can be graded live), and every model's exact
predictions on a grid of N* (and of its parameter). Nothing is computed from
scratch in the browser except interpolation, so it works offline.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # works under python -I too

import contourpy  # noqa: E402
import numpy as np  # noqa: E402

import config as cfg  # noqa: E402

TEMPLATE = cfg.PROJECT_DIR / "explorer_template.html"
OUT = cfg.FIGURE_DIR / "inflation_explorer.html"
OUT_ARTIFACT = cfg.DERIVED_DIR / "inflation_explorer_artifact.html"
THEME_BUTTONS = """      <div class="seg themebtn" role="group" aria-label="Theme">
        <button type="button" data-theme-set="light" aria-pressed="false">Light</button>
        <button type="button" data-theme-set="dark" aria-pressed="false">Dark</button>
      </div>"""

PLANES = {  # key -> (label, colour slot: 0 neutral, 1-3 data hues, outline only?)
    "planck18": ("Planck 2018 alone (no B-mode data)", 0, True),
    "planck18_bk18_ours": ("Planck 2018 + BK18 (our BK18 likelihood)", 1, False),
    "pactlb_bk18": ("ACT DR6: P-ACT-LB + BK18", 2, False),
    "cmbspa_desi_x_bk18": ("SPT-3G: CMB-SPA + DESI, with our BK18 r", 3, True),
}


def _load(name: str):
    with np.load(cfg.DERIVED_DIR / f"{name}.npz") as z:
        arrays = {k: z[k] for k in z.files}
    return arrays, json.loads((cfg.DERIVED_DIR / f"{name}.json").read_text())


def _g(v: float, sig: int = 4) -> float:
    return float(f"{v:.{sig}g}")


def _rings(x, y, p, level, max_pts=260) -> list:
    gen = contourpy.contour_generator(x, y, p.T, fill_type=contourpy.FillType.OuterOffset)
    points, offsets = gen.filled(level, float(p.max()) * 2)
    rings = []
    for pts, offs in zip(points, offsets):
        for a, b in zip(offs[:-1], offs[1:]):
            ring = pts[a:b]
            if len(ring) < 4:
                continue
            step = max(1, int(math.ceil(len(ring) / max_pts)))
            ring = np.vstack([ring[::step], ring[:1]])
            rings.append([[round(float(u), 5), _g(float(v))] for u, v in ring])
    return rings


def _credibility_grid(p: np.ndarray) -> np.ndarray:
    """Mass of the plane at densities above each cell's (the HPD credibility)."""
    flat = p.ravel()
    order = np.argsort(-flat, kind="stable")
    above = np.empty_like(flat)
    csum = np.cumsum(flat[order])
    above[order] = csum - flat[order]
    return (above / flat.sum()).reshape(p.shape)


def build() -> Path:
    c, cs = _load("constraints")
    eg, egs = _load("explorer_grid")
    x, y = c["plane_x"], c["plane_y"]
    planes = []
    for key, (label, slot, outline) in PLANES.items():
        p = c[f"plane__{key}__p"].astype(float)
        lev = c[f"plane__{key}__levels"]
        ylo = np.geomspace(1e-6, y[0], 40, endpoint=False)
        y_ext = np.concatenate([ylo, y])
        p_ext = np.concatenate([np.repeat(p[:, :1], ylo.size, axis=1), p], axis=1)
        cred = _credibility_grid(p)[::3, ::3]
        planes.append({
            "key": key, "label": label, "slot": slot, "outline": outline,
            "lin": {"r95": _rings(x, y, p, lev[1]), "r68": _rings(x, y, p, lev[0])},
            "log": {"r95": _rings(x, y_ext, p_ext, lev[1]), "r68": _rings(x, y_ext, p_ext, lev[0])},
            "cred": (np.round(cred * 1000)).astype(int).tolist(),
        })
    grid = {"x0": float(x[::3][0]), "dx": float(x[3] - x[0]), "nx": int(x[::3].size),
            "y0": float(y[::3][0]), "dy": float(y[3] - y[0]), "ny": int(y[::3].size)}
    models = []
    n_star = [float(v) for v in eg["n_star"]]
    scen = egs["reheating_scenarios"]
    for name, spec in cfg.MODELS.items():
        params = eg[f"{name}__param"]
        entry = {"key": name, "label": spec["label"], "kind": spec["kind"],
                 "param": None if "param" not in spec else {
                     "name": spec["param"], "range": list(spec["range"]),
                     "values": [_g(v, 5) for v in params]}}
        for k, sig in (("ns", 6), ("r", 4), ("alpha_s", 3), ("ns_sr1", 6), ("r_sr1", 4),
                       ("V_star_quarter_gev", 4), ("delta_phi", 4)):
            entry[k] = [[_g(float(v), sig) for v in row] for row in eg[f"{name}__{k}"]]
        reh = eg[f"{name}__reheating"]
        entry["reheating"] = [[None if not np.isfinite(v[0]) else round(float(v[0]), 2) for v in row]
                              for row in reh]
        models.append(entry)
    d3 = cs["d3"]
    data = {
        "planes": planes, "grid": grid, "models": models, "n_star": n_star,
        "reheating_labels": scen, "d3": {"ns": d3["value"], "err": d3["error"]},
        "r95": {k: (cs["ours"]["r95"] if k == "planck18_bk18_ours" else cs["planes"][k]["r95_from_plane"])
                for k in PLANES},
        "footnote": (
            "Regions: 68% and 95% credible regions from the published MCMC chains — Planck 2018 (Planck "
            "Legacy Archive grid, r free); Planck 2018 + BK18 = Planck's chain weighted by this project's "
            "own re-implementation of the BICEP/Keck BK18 likelihood; ACT DR6.02 P-ACT-LB + BK18 (LAMBDA); "
            "SPT-3G D1 CMB-SPA + DESI DR2 (LAMBDA, 2026 bug-fix chains; fitted without r, so its n_s "
            "posterior is multiplied by our BK18 r likelihood — an approximation checked on ACT). Models: "
            "each potential solved exactly (inflaton + Friedmann equations, Mukhanov–Sasaki modes); n_s and "
            "r at k = 0.05 Mpc⁻¹, interpolated between the computed grid points (N* every 2.5 e-folds, 24 "
            "parameter values). Verdicts: σ-equivalent of the smallest 2D credible region containing the "
            "point (1σ ↔ 68.3%, 2σ ↔ 95.4%). Part of project P1 of the astrophysics atlas: 67-inflation/."),
    }
    fragment = TEMPLATE.read_text().replace("__DATA__", json.dumps(data, separators=(",", ":")))
    page = fragment.replace("<!--THEME_BUTTONS-->", THEME_BUTTONS)
    head, body = page.split("<main>", 1)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
                   '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
                   + head + "</head>\n<body>\n<main>" + body + "</body>\n</html>\n")
    OUT_ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    OUT_ARTIFACT.write_text(fragment.replace("<!--THEME_BUTTONS-->\n", "").replace("<!--THEME_BUTTONS-->", ""))
    return OUT


def main() -> int:
    p = build()
    print(f"  {p.name} ({p.stat().st_size / 1e3:.0f} kB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
