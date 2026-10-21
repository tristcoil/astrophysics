"""Read Cohen, De Rujula & Glashow's (1998) published curves from their figures.

Their arXiv source (astro-ph/9707087) ships the figures as vector PostScript:
three made with XMgr (absolute page coordinates, with labelled tick marks we
calibrate against) and one with Mathematica (normalised coordinates with an
explicit ``MathScale`` line). The curves are polylines, so they can be read
back exactly -- no digitising by eye. They are the answer key for our own
calculation of the annihilation rate and the gamma-ray signal in
``domains.py``:

* Fig. 3 (``curr.eps``): J(y), annihilations per cm^2 per s at a domain
  boundary; solid = their numerical fluid solution, dashed = the
  "hydrodynamic" estimate J = (3/4)^4 n c_s with their numerical temperature.
* Fig. 2 (``temp.eps``): plasma temperature T(y) in eV; solid = inside the
  annihilation electrons' range, dashed = UV-photon heating only, and the
  conventional (no-annihilation) universe.
* Fig. 4 (``flux.eps``): predicted diffuse gamma-ray spectrum today for
  domain sizes d0 = 20 Mpc (upper) and 1000 Mpc (lower), photons
  cm^-2 s^-1 MeV^-1 sr^-1.

Their fiducial cosmology: Omega = 1, Lambda = 0, h = 0.75, eta = 2e-10, and
annihilations counted only for 1100 > y > 20 (y = 1 + z).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

import config as cfg

SRC = cfg.ARXIV_DIR / "astro-ph_9707087"

#: Their fiducial parameters (Section 1).
FIDUCIAL = {"h": 0.75, "omega_m": 1.0, "omega_l": 0.0, "eta": 2e-10, "y_r": 1100.0,
            "y_s": 20.0, "gbar": 3.8, "a_dedl": 10.0}


@dataclass
class Path_:
    color: tuple
    width: float
    dash: str
    xy: list = field(default_factory=list)


def _xmgr_paths(text: str) -> list[Path_]:
    """Stroked polylines of an XMgr EPS, with their colour, width and dash."""
    paths, cur = [], None
    color, width, dash = (0.0, 0.0, 0.0), 1.0, "[]"
    for line in text.splitlines():
        s = line.strip()
        m = re.fullmatch(r"(-?\d+) (-?\d+) ([ml])", s)
        if m:
            x, y, op = float(m.group(1)), float(m.group(2)), m.group(3)
            if op == "m" or cur is None:
                if cur is not None and len(cur.xy) > 1:
                    paths.append(cur)
                cur = Path_(color, width, dash)
            cur.xy.append((x, y))
            continue
        if s == "stroke":
            if cur is not None and len(cur.xy) > 1:
                paths.append(cur)
            cur = None
        elif s.endswith("setrgbcolor"):
            color = tuple(float(v) for v in s.split()[:3])
        elif s.endswith("setlinewidth"):
            width = float(s.split()[0])
        elif s.endswith("setdash"):
            dash = s.rsplit(" ", 2)[0]
    return paths


def _xmgr_ticks(text: str) -> tuple[dict, dict]:
    """Tick-label positions: x labels '(100) CS' and y labels '(10)^n'."""
    xt, yt = {}, {}
    for m in re.finditer(r"(\d+) (\d+) translate\s+0 rotate\s+0 -20\s+m\s+\(([\d.]+)\) CS", text):
        xt[float(m.group(3))] = float(m.group(1))
    # y labels are right-justified on their full text, then drawn: a power of
    # ten as '(10n) RJ (10) show ... (n) show', a plain number as
    # '(0.1) RJ (0.1) show'.
    for m in re.finditer(r"(\d+) (\d+) translate\s+0 rotate\s+0 -20\s+m\s+\(([\-\d.]+)\) RJ"
                         r"\s+\(([\-\d.]+)\) show", text):
        label, shown = m.group(3), m.group(4)
        if label == shown:
            yt[float(label)] = float(m.group(2))
        elif shown == "10" and label.startswith("10"):
            yt[10.0 ** int(label[2:])] = float(m.group(2))
    return xt, yt


def _log_calibration(ticks: dict):
    vals = np.array(sorted(ticks))
    pos = np.array([ticks[v] for v in vals])
    slope, icpt = np.polyfit(np.log10(vals), pos, 1)
    resid = pos - (icpt + slope * np.log10(vals))
    assert np.max(np.abs(resid)) < 1.5, "tick marks not on a log scale"
    return lambda p: 10.0 ** ((np.asarray(p) - icpt) / slope)


def read_xmgr(name: str):
    text = (SRC / name).read_text(errors="replace")
    xt, yt = _xmgr_ticks(text)
    fx, fy = _log_calibration(xt), _log_calibration(yt)
    out = []
    for p in _xmgr_paths(text):
        xy = np.array(p.xy)
        out.append({"color": p.color, "width": p.width, "dash": p.dash,
                    "x": fx(xy[:, 0]), "y": fy(xy[:, 1]), "n": len(xy)})
    return out


def _data_curves(paths: list, min_points: int = 30) -> list:
    """Long polylines only (axes, ticks and legend swatches are short)."""
    return [p for p in paths if p["n"] >= min_points]


def annihilation_rate() -> dict:
    """Fig. 3: J(y) in cm^-2 s^-1, numerical (solid) and hydrodynamic (dashed)."""
    curves = _data_curves(read_xmgr("cdg98_J.eps"))
    out = {}
    for c in curves:
        key = "numerical" if c["dash"] == "[]" else "hydrodynamic_estimate"
        order = np.argsort(c["x"])
        out[key] = {"y": c["x"][order], "J_cm2_s": c["y"][order]}
    assert set(out) == {"numerical", "hydrodynamic_estimate"}, out.keys()
    return out


def temperatures() -> dict:
    """Fig. 2: T(y) in eV for the three cases."""
    curves = _data_curves(read_xmgr("cdg98_T.eps"))
    out = {}
    for c in curves:
        r, g, b = c["color"]
        if b > 0.5:
            key = "conventional"
        elif r > 0.5:
            key = "uv_heating"
        else:
            key = "electron_heating"
        order = np.argsort(c["x"])
        x, y = c["x"][order], c["y"][order]
        if key in out:                                   # a curve split into pieces
            x = np.concatenate([out[key]["y"], x])
            y = np.concatenate([out[key]["T_eV"], y])
            o = np.argsort(x)
            x, y = x[o], y[o]
        out[key] = {"y": x, "T_eV": y}
    return out


def gamma_spectra() -> dict:
    """Fig. 4: predicted CDG spectra (photons cm^-2 s^-1 MeV^-1 sr^-1).

    The figure was assembled in a drawing program (MacDraw-style PostScript,
    which draws the visible axes, labels and data points) from two
    Mathematica pictures. Both pictures hold the *same* curve, computed for
    d0 = 20 Mpc; the second is simply placed 75 pt lower on the page to stand
    for d0 = 1000 Mpc. We therefore return:

    * ``computed``: the curve in Mathematica's own coordinates (its scaling
      line reads x = a + b log10(E/MeV), y = c + d log10(F), confirmed by its
      tick labels) -- what the authors calculated, for d0 = 20 Mpc; the signal
      is exactly proportional to 1/d0;
    * ``as_drawn``: both curves read against the drawn axes, using the page
      matrices Ghostscript reports for each picture (cached in
      ``data/derived/cdg98_flux_ctm.json``). The upper curve sits on the
      drawn axes to 0.01 dex in flux but 0.05 dex right in energy; the lower
      one is 1.58 decades below it (x 1/38), not the x 1/50 the text implies.
    """
    text = (SRC / "cdg98_flux.eps").read_text(errors="replace")
    starts = [m.start() for m in re.finditer(r"\nMathPictureStart\n%% Graphics", text)]
    pictures = []
    for start in starts:
        body = text[start:text.index("MathPictureEnd", start)]
        curves = _mathematica_curves(body)
        assert len(curves) == 1
        pictures.append(curves[0])
    assert len(pictures) == 2
    (e1, f1, xy1), (e2, f2, xy2) = pictures
    # same curve (the lower copy is cut at 39.5 MeV to end on the frame)
    grid = np.log(np.geomspace(0.3, 39.0, 50))
    same = np.interp(grid, np.log(e1), np.log10(f1)) - np.interp(grid, np.log(e2), np.log10(f2))
    assert np.max(np.abs(same)) < 1e-4, "the two pictures are expected to hold one curve"
    ctms = _picture_matrices()
    drawn = [_as_drawn(xy, m, ctms["default"]) for xy, m in zip((xy1, xy2), ctms["pictures"])]
    drawn.sort(key=lambda t: -np.median(t[1]))
    return {"computed": {"d0_mpc": 20.0, "E_MeV": e1, "flux": f1},
            "as_drawn": {"d0_20Mpc": {"E_MeV": drawn[0][0], "flux": drawn[0][1]},
                         "d0_1000Mpc": {"E_MeV": drawn[1][0], "flux": drawn[1][1]}}}


#: The drawing program's page matrix (``/mT``) and its axes, in its own
#: units (1/300 inch, y down), read from the tick marks in the file.
MD_MATRIX = (0.24, 29.999, 760.769)          # x = 0.24 u + 29.999, y = 760.769 - 0.24 v
MD_X_TICKS = {1.0: 675.0, 10.0: 1406.0}      # E (MeV) -> u
MD_Y_TICKS = {1.0: 140.0, 1e-6: 1331.0}      # F -> v


def _picture_matrices() -> dict:
    """Current transformation matrices after each picture's MathScale.

    Ghostscript runs a copy of the figure with a print statement inserted
    after each ``MathScale`` (null device, 72 dpi, -dSAFER). The result is
    cached, so Ghostscript is only needed once.
    """
    import json
    import shutil
    import subprocess
    import tempfile

    cache = cfg.DERIVED_DIR / "cdg98_flux_ctm.json"
    if cache.exists():
        return json.loads(cache.read_text())
    if shutil.which("gs") is None:
        raise RuntimeError("Ghostscript is needed once to place CDG's Fig. 4 curves")
    lines = (SRC / "cdg98_flux.eps").read_text(errors="replace").splitlines()
    out = []
    for i, line in enumerate(lines):
        out.append(line)
        if i == 4:
            out.append("(DEFAULT ) print matrix currentmatrix == flush")
        if line.strip().endswith("] MathScale"):
            out.append("(MATHCTM ) print matrix currentmatrix == flush")
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "flux_ctm.eps"
        p.write_text("\n".join(out))
        res = subprocess.run(["gs", "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=nullpage",
                              "-r72", str(p)], capture_output=True, text=True, timeout=120)
    mats = {"default": None, "pictures": []}
    for line in res.stdout.splitlines():
        nums = [float(v) for v in re.findall(r"-?[\d.]+(?:e-?\d+)?", line.split("[", 1)[-1])]
        if line.startswith("DEFAULT"):
            mats["default"] = nums
        elif line.startswith("MATHCTM"):
            mats["pictures"].append(nums)
    assert mats["default"] and len(mats["pictures"]) == 2, res.stdout + res.stderr
    cache.write_text(json.dumps(mats, indent=1))
    return mats


def _as_drawn(xy: np.ndarray, m: list, default: list):
    """Normalised Mathematica coordinates -> the drawn axes' (E, F)."""
    a, b, c, d, e, f = m
    xd = a * xy[:, 0] + c * xy[:, 1] + e                # device space
    yd = b * xy[:, 0] + d * xy[:, 1] + f
    A, B, C, D, E, F = default                          # device -> page points
    det = A * D - B * C
    xp = (D * (xd - E) - C * (yd - F)) / det
    yp = (-B * (xd - E) + A * (yd - F)) / det
    s, x0, y0 = MD_MATRIX
    u, v = (xp - x0) / s, (y0 - yp) / s
    (e1, u1), (e2, u2) = sorted(MD_X_TICKS.items())
    (f1, v1), (f2, v2) = sorted(MD_Y_TICKS.items())
    energy = 10 ** (np.log10(e1) + (u - u1) * np.log10(e2 / e1) / (u2 - u1))
    flux = 10 ** (np.log10(f1) + (v - v1) * np.log10(f2 / f1) / (v2 - v1))
    return energy, flux


def _mathematica_curves(body: str) -> list:
    a, b, c, d = (float(v) for v in re.search(
        r"% Scaling calculations\s+([\d.\-]+) ([\d.\-]+) ([\d.\-]+) ([\d.\-]+) \[", body).groups())
    # the tick labels confirm the scaling: '1' at x = a, '10.' at a + b,
    # '0.1' at y = c - d
    assert abs(float(re.search(r"\[\(1\)\] ([\d.]+)", body).group(1)) - a) < 1e-4
    assert abs(float(re.search(r"\[\(10\.\)\] ([\d.]+)", body).group(1)) - (a + b)) < 1e-3
    assert abs(float(re.search(r"\[\(0\.1\)\] [\d.]+ ([\d.]+)", body).group(1)) - (c - d)) < 1e-4
    paths, cur, width = [], None, None
    for line in body.splitlines():
        s = line.strip()
        m = re.fullmatch(r"([\d.\-]+) ([\d.\-]+) ([mL])", s)
        if m:
            x, y = float(m.group(1)), float(m.group(2))
            if m.group(3) == "m":
                if cur and len(cur[1]) > 1:
                    paths.append(cur)
                cur = (width, [(x, y)])
            elif cur is not None:
                cur[1].append((x, y))
            continue
        mw = re.fullmatch(r"([\d.]+) w", s)
        if mw:
            width = float(mw.group(1))
        if s == "s" and cur is not None:
            if len(cur[1]) > 1:
                paths.append(cur)
            cur = None
    thick = [np.array(xy) for w, xy in paths if w is not None and w >= 0.0035 and len(xy) > 5]
    # join segments that continue one another into whole curves
    curves = []
    for seg in thick:
        for cv in curves:
            if np.allclose(cv[-1], seg[0], atol=1e-4):         # seg continues cv
                cv.extend(seg[1:].tolist())
                break
            if np.allclose(seg[-1], cv[0], atol=1e-4):         # seg precedes cv
                cv[:0] = seg[:-1].tolist()
                break
        else:
            curves.append(seg.tolist())
    out = []
    for cv in curves:
        xy = np.array(cv)
        e = 10.0 ** ((xy[:, 0] - a) / b)
        f = 10.0 ** ((xy[:, 1] - c) / d)
        out.append((e, f, xy))
    return out


def scales() -> list:
    """Fig. 1 (dist.eps): horizon, look-back size of a 20 Mpc domain, the
    reheated (electron range), depletion and annihilation zones, in Mpc."""
    return [c for c in _data_curves(read_xmgr("cdg98_dist.eps"), 10)]


if __name__ == "__main__":
    j = annihilation_rate()
    for k, v in j.items():
        print(k, len(v["y"]), v["y"][[0, -1]], v["J_cm2_s"][[0, -1]])
    t = temperatures()
    for k, v in t.items():
        print(k, len(v["y"]), v["y"][[0, -1]], v["T_eV"][[0, -1]], v["T_eV"].max())
    g = gamma_spectra()
    c = g["computed"]
    print("computed (20 Mpc):", np.interp([1, 3, 10, 30], c["E_MeV"], c["flux"]))
    for k, v in g["as_drawn"].items():
        print(k, "as drawn:", np.interp([1, 3, 10, 30], v["E_MeV"], v["flux"]))
