"""Build the interactive page, figures/wimp_freeze_out_lab.html, from this project's results.

    python explorer.py

The page (explorer_template.html) is a body fragment holding its own copy of the freeze-out solver
(the implicit stepper of freezeout.yield_grid, in JavaScript). This script injects the data the
solver needs - Saikawa & Shirai's g*(T), the clock t(T), the thermal-relic curves for four particle
types and s/p-wave, and our Fermi-LAT dwarf limits - and writes two copies: a standalone document
for opening from disk, and the bare fragment for publishing as an Artifact
(data/derived/wimp_freeze_out_lab_artifact.html).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # works under python -I too

import numpy as np  # noqa: E402

import config as C  # noqa: E402
import freezeout as fo  # noqa: E402
import thermo  # noqa: E402

TEMPLATE = C.HERE / "explorer_template.html"
OUT = C.FIG_DIR / "wimp_freeze_out_lab.html"
OUT_ARTIFACT = C.DERIVED_DIR / "wimp_freeze_out_lab_artifact.html"


def _r(a, digits=7):
    return [float(f"{v:.{digits}g}") for v in np.asarray(a, dtype=float)]


def page_data() -> dict:
    F = dict(np.load(C.DERIVED_DIR / "figdata.npz"))
    e = thermo.load("ss2020")
    lnT = np.linspace(np.log(1e-9), np.log(1e6), 2400)
    T = np.exp(lnT)
    Tt, t = e.time_of(np.logspace(-6, 5.5, 1200))
    thermal = {}
    curves = {"2s": "th_main", "4s": "th_type_dirac", "1s": "th_type_real_scalar", "3s": "th_type_vector",
              "2p": "th_pwave", "4p": "th_pwave_dirac", "1p": "th_pwave_real_scalar", "3p": "th_pwave_vector"}
    for key, name in curves.items():
        if name in F:
            thermal[key] = {"m": _r(F["th_m"], 6), "sv": _r(F[name], 6)}
    for key in curves:                     # fall back to the Majorana curve of the same wave
        thermal.setdefault(key, thermal["2" + key[1]])
    fermi = {}
    for ch in ("bb", "tau"):
        m, ul = F[f"dw_{ch}_pppc_benchmark_m"], F[f"dw_{ch}_pppc_benchmark_ul"]
        ok = np.isfinite(ul)
        fermi[ch] = {"m": _r(m[ok], 6), "ul": _r(ul[ok], 6)}
    fermi["mmin"], fermi["mmax"] = float(C.GAMMA_MASS_MIN), 1.0e4
    return {
        "consts": {"M_PL": C.M_PL, "GEV2_TO_CM3S": C.GEV2_TO_CM3S, "OMEGA_TARGET": C.OMEGA_C_H2,
                   "OMEGA_PER_MY": fo.omega_h2_from_Y(1.0, 1.0, e.g_s0)},
        "eos": {"lnT": _r(lnT, 8), "gs": _r(e.g_s(T), 7), "gh": _r(e.gstar_half(T), 7),
                "time_mlnT": _r(-np.log(Tt), 8), "time_lnt": _r(np.log(t), 8)},
        "thermal": thermal,
        "fermi": fermi,
    }


def build(results: dict | None = None) -> Path:
    frag = TEMPLATE.read_text().replace("__DATA__", json.dumps(page_data(), separators=(",", ":")))
    OUT_ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    OUT_ARTIFACT.write_text(frag)
    i = frag.index("<main")
    head, body = frag[:i], frag[i:]
    doc = ("<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
           "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, viewport-fit=cover\">\n"
           + head + "</head>\n<body>\n" + body + "</body>\n</html>\n")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(doc)
    return OUT


if __name__ == "__main__":
    print("  wrote", build().relative_to(C.HERE))
