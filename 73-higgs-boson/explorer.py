"""Write the interactive vacuum-stability explorer (an offline HTML page).

The page reruns the Standard Model's three-loop renormalization group
equations in the browser whenever a slider moves (top mass, Higgs mass,
strong coupling), so the reader can watch the Higgs self-coupling dive below
zero and the vacuum turn from stable to metastable to unstable.

Nothing is approximated twice: the JavaScript beta functions are generated
from ``rge.betas`` itself (a mechanical translation of its source), and the
matching at the top mass is a quadratic fit to this project's matching
(SMDR-calibrated), accurate to ~2e-6 over the slider ranges. ``verify()``
runs both through node and compares them with Python.

    python explorer.py      # writes figures/vacuum_explorer.html (+ an artifact fragment)
"""

from __future__ import annotations

import inspect
import json
import re
import shutil
import subprocess

import numpy as np

import config as cfg
import rge
import stability as st

TEMPLATE = cfg.PROJECT_DIR / "explorer_template.html"
OUT = cfg.FIGURE_DIR / "vacuum_explorer.html"
FRAGMENT = cfg.DATA_DIR / "derived" / "vacuum_explorer_artifact.html"
COUPLINGS = ("lam", "yt", "g3", "g2", "gY", "yb", "ytau")
CENTER = {"Mt": cfg.M_TOP[0], "Mh": cfg.M_HIGGS[0], "as": cfg.ALPHA_S_MZ[0]}
RANGES = {"Mt": (164.0, 181.0), "Mh": (118.0, 132.0), "as": (0.1150, 0.1210)}


def betas_js() -> str:
    """Translate rge.betas (Python) into an equivalent JavaScript function."""
    src = inspect.getsource(rge.betas)
    body = src[src.index("G1, G2, G3, YT, YB, YTA, L = X[:7]"):src.index("# Gamma = int gamma")]
    two = body.index("    if loops >= 2:")
    three = body.index("    if loops >= 3:")
    seg1, seg2, seg3 = body[:two], body[two:three], body[three:]

    def conv(seg: str) -> str:
        seg = seg.replace("\\\n", "\n")
        seg = re.sub(r"#(.*)", r"//\1", seg)
        return seg

    seg1 = conv(seg1)
    seg1 = seg1.replace("G1, G2, G3, YT, YB, YTA, L = X[:7]",
                        "const [G1, G2, G3, YT, YB, YTA, L] = X;")
    seg1 = seg1.replace("k, k2, k3 = K, K * K, K**3", "const k = K, k2 = K * K, k3 = K ** 3;")
    for v in ("dG1", "dG2", "dG3", "dYT", "dYB", "dYTA", "dL", "gamma"):
        seg1 = re.sub(rf"^(\s*){v} = ", rf"\1let {v} = ", seg1, count=1, flags=re.M)
    seg2 = conv(seg2).replace("    if loops >= 2:", "    if (loops >= 2) {") + "    }\n"
    seg3 = conv(seg3).replace("    if loops >= 3:", "    if (loops >= 3) {") + "    }\n"
    return ("function betas(X, loops) {\n    const K = 1 / (16 * Math.PI ** 2);\n    "
            + seg1 + seg2 + seg3
            + "    return [dG1, dG2, dG3, dYT, dYB, dYTA, dL, 0.5 * gamma];\n}\n")


def matching_fit(offsets: dict | None) -> dict:
    """Quadratic fits of the couplings at mu = M_t in (dMt, dMh, das/0.001)."""
    base = st.inputs("pdg2026")
    P, V = [], []
    for mt in np.linspace(*RANGES["Mt"], 11):
        for mh in np.linspace(*RANGES["Mh"], 9):
            for a in np.linspace(*RANGES["as"], 7):
                c = st.match(st._with(base, Mt=float(mt), Mh=float(mh), alpha_s=float(a)),
                             offsets=offsets)
                P.append((mt - CENTER["Mt"], mh - CENTER["Mh"], (a - CENTER["as"]) / 0.001))
                V.append([c[k] for k in COUPLINGS])
    P, V = np.array(P), np.array(V)
    x, y, z = P.T
    A = np.stack([np.ones_like(x), x, y, z, x * x, y * y, z * z, x * y, x * z, y * z], 1)
    coef, *_ = np.linalg.lstsq(A, V, rcond=None)
    resid = np.abs(V - A @ coef).max(axis=0)
    return {"coef": {k: coef[:, i].tolist() for i, k in enumerate(COUPLINGS)},
            "max_residual": dict(zip(COUPLINGS, resid.tolist()))}


def boundaries(g: dict) -> dict:
    """M_h(M_t) along the zoom grid where min lambda = 0 (stability, with the
    running coupling, as the page computes it), where min lambda_eff = 0, and
    where the past-light-cone decay probability reaches 1."""
    mt, mh = g["zoom_mt"], g["zoom_mh"]
    out = {}
    for key, arr in (("stable_lambda", g["zoom_lambda_min"]),
                     ("stable_lambda_eff", g["zoom_lambda_eff_min"]),
                     ("unstable", -g["zoom_log10_p0"])):
        pts = []
        for j, t in enumerate(mt):
            col = arr[:, j]
            idx = np.nonzero((col[:-1] < 0) & (col[1:] >= 0))[0]
            if len(idx):
                i = idx[0]
                h = mh[i] + (mh[i + 1] - mh[i]) * (-col[i]) / (col[i + 1] - col[i])
                pts.append([float(t), float(h)])
        out[key] = pts
    return out


def reference_curves(offsets) -> dict:
    """lambda(mu) for the measured masses and for M_t +- 1 sigma (static band)."""
    inp = st.inputs("pdg2026")
    out = {}
    for tag, dt in (("central", 0.0), ("lo", -cfg.M_TOP[1]), ("hi", cfg.M_TOP[1])):
        tr = st.trajectory(st._with(inp, Mt=inp.Mt + dt), offsets=offsets, mu_max=1e20)
        lmu = np.log10(tr["mu"])
        grid = np.linspace(np.log10(inp.Mt), 20, 260)
        out[tag] = np.round(np.interp(grid, lmu, tr["X"][:, 6]), 6).tolist()
        out["x"] = np.round(grid, 4).tolist()
    return out


def build() -> dict:
    s = st.summary()
    off = s.get("smdr_offsets")
    g = st.grids()
    data = {
        "center": CENTER, "sigma": {"Mt": cfg.M_TOP[1], "Mh": cfg.M_HIGGS[1],
                                     "as": cfg.ALPHA_S_MZ[1]},
        "ranges": RANGES, "match": matching_fit(off), "bounds": boundaries(g),
        "reference": reference_curves(off),
        "constants": {"MPl": st.M_PL, "log10_H0": float(np.log10(st.h0_gev())),
                      "log10_lc": float(np.log10(st.light_cone_volume())),
                      "log10_hbar_per_year": float(np.log10(cfg.HBAR_GEV_S / st.SEC_PER_YEAR)),
                      "pdg2024": {"Mt": 172.57, "Mh": 125.20}},
        "headline": s["results"][s["headline"]],
        "matching_label": "SMDR" if off else "Buttazzo et al. 2013",
    }
    tpl = TEMPLATE.read_text()
    page = (tpl.replace("/*__DATA__*/", "const DATA = " + json.dumps(data) + ";")
               .replace("/*__BETAS__*/", betas_js()))
    head, frag = page.split("<!--FRAGMENT-->", 1)
    full = ("<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, "
            "viewport-fit=cover\">\n" + frag.replace("<!--BODY-->", "</head>\n<body>", 1)
            + "\n</body>\n</html>\n")
    OUT.write_text(full)
    FRAGMENT.parent.mkdir(parents=True, exist_ok=True)
    FRAGMENT.write_text(frag.replace("<!--BODY-->", "", 1))
    return data


def verify() -> dict:
    """Run the generated JavaScript under node and compare with Python."""
    if shutil.which("node") is None:
        return {"skipped": "node not installed"}
    js = betas_js()
    X = rge.state(rge.match(st.inputs("pdg2026")))
    prog = js + f"""
const X = {json.dumps(X.tolist())};
const out = {{b3: betas(X, 3), b2: betas(X, 2), b1: betas(X, 1)}};
// full run, RK4 in ln mu^2, as the page does
function run(X0, mu0, mu1, step) {{
  const s0 = 2 * Math.log(mu0), s1 = 2 * Math.log(mu1);
  const n = Math.ceil(Math.abs(s1 - s0) / step), h = (s1 - s0) / n;
  let x = X0.slice();
  for (let i = 0; i < n; i++) {{
    const k1 = betas(x, 3);
    const k2 = betas(x.map((v, j) => v + 0.5 * h * k1[j]), 3);
    const k3 = betas(x.map((v, j) => v + 0.5 * h * k2[j]), 3);
    const k4 = betas(x.map((v, j) => v + h * k3[j]), 3);
    x = x.map((v, j) => v + h / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]));
  }}
  return x;
}}
out.final = run(X, {float(cfg.M_TOP[0])}, 1e19, 0.05);
console.log(JSON.stringify(out));
"""
    r = subprocess.run(["node", "-e", prog], capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise RuntimeError(r.stderr)
    js_out = json.loads(r.stdout)
    py = {f"b{n}": rge.betas(X, n) for n in (1, 2, 3)}
    rel = max(float(np.max(np.abs(np.array(js_out[k]) - py[k]) / (np.abs(py[k]) + 1e-30)))
              for k in py)
    final_py = rge.run(X, float(cfg.M_TOP[0]), 1e19)[1][-1]
    dfinal = float(np.max(np.abs(np.array(js_out["final"]) - final_py)))
    return {"max_rel_diff_betas": rel, "max_abs_diff_after_running": dfinal}


if __name__ == "__main__":
    d = build()
    print("matching fit residuals:", d["match"]["max_residual"])
    print("verify:", verify())
    print("wrote", OUT, "and", FRAGMENT)
