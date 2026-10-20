"""Build the interactive page figures/freezeout_thermometer.html.

The page runs the same hadron resonance gas in the browser: the ~4,200 mass
points of the PDG 2025 list (Breit-Wigner nodes included), the 64-point
momentum quadrature, the conservation-law solver and each dataset's
feed-down matrix are embedded as data, so moving a slider recomputes every
yield from scratch. ``check()`` compares the page's numbers with Python's.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import config as cfg  # noqa: E402
import datasets as D  # noqa: E402
import hrg  # noqa: E402
import thermal_fit as F  # noqa: E402
import universe as U  # noqa: E402

TEMPLATE = cfg.PROJECT_DIR / "explorer_template.html"
OUT = cfg.FIGURE_DIR / "freezeout_thermometer.html"


def _g(x: float, digits: int = 10) -> float:
    return float(f"{x:.{digits}g}")


def _dataset_block(ds, best: dict) -> dict:
    model = F.YieldModel(ds, F.FitSettings())
    rows = []
    for o in range(len(ds.yields)):
        nz = np.flatnonzero(model.A[o])
        rows.append([[int(i), _g(model.A[o, i])] for i in nz])
    pred = F.model_yields(ds, best["T"], best["muB"], best["V"])
    return {"name": ds.name, "title": ds.title, "qb": ds.q_over_b,
            "sqrt_s": ds.sqrt_s, "charm": [_g(c) for c in model.c],
            "yields": [{"key": y.key, "html": y.html, "value": y.value, "err": y.err,
                        "fitted": y.fitted, "note": y.note} for y in ds.yields],
            "A": rows, "best": best, "python_model": [float(v) for v in pred]}


def payload() -> dict:
    gas = hrg.model()
    hl = gas.hl
    main = json.loads((cfg.DERIVED_DIR / "main_fit.json").read_text())
    star = json.loads((cfg.DERIVED_DIR / "star.json").read_text())
    phase = json.loads((cfg.DERIVED_DIR / "phase_diagram.json").read_text())
    p = main["params"]
    sets = [_dataset_block(D.alice("0-10"), {"T": p["T"]["value"], "muB": p["muB"]["value"],
                                              "V": p["V"]["value"], "T_err": p["T"]["error"],
                                              "muB_err": p["muB"]["error"]})]
    for e in reversed(D.STAR_ENERGIES):
        r = star["ours"][f"{e}"]["eq"]
        sets.append(_dataset_block(D.star(e), {"T": r["T"], "muB": r["muB"], "V": r["V"],
                                               "T_err": r["T_err"], "muB_err": r["muB_err"]}))
    cal = U.p2_calendar()
    T_cal = np.logspace(0, np.log10(3000), 160)
    t_cal = U.age_at(T_cal, cal)
    lat = phase["lattice"]["hotqcd"]
    keep = [i for i, m in enumerate(lat["mu_B"]) if m <= 300]
    up = phase["universe_path"]
    return {
        "quad": {"x": [float(v) for v in gas._x], "w": [float(v) for v in gas._w]},
        "pts": {"mass": [float(v) for v in gas.pt_mass], "pref": [float(v) for v in gas.pt_pref],
                "stat": [int(v) for v in gas.pt_stat], "sp": [int(v) for v in gas.pt_species]},
        "species": {"B": hl.B.tolist(), "Q": hl.Q.tolist(), "S": hl.S.tolist()},
        "datasets": sets,
        "phase": {
            "lattice": {k: [lat[k][i] for i in keep] for k in ("mu_B", "T", "lo", "hi")},
            "star": [{"label": f"STAR {float(k):g} GeV (published)", **v}
                     for k, v in phase["star_hepdata"].items()],
            "andronic": phase["published"]["andronic"],
            "other": sum((phase["published"].get(g, []) for g in ("other_a", "other_b", "other_c")), []),
            "ours": phase["ours"],
            "universe": {"T": up["T_MeV"], "muB": up["muB_MeV"]},
        },
        "calendar": {"T": [float(v) for v in T_cal], "t_us": [float(v) * 1e6 for v in t_cal]},
        "n_species": len(hl), "n_points": len(gas.pt_mass),
    }


def main() -> Path:
    data = json.dumps(payload(), separators=(",", ":"))
    html = TEMPLATE.read_text().replace("__PAYLOAD__", data)
    OUT.write_text(html)
    print(f"  wrote {OUT.relative_to(cfg.PROJECT_DIR)} ({OUT.stat().st_size / 1e3:.0f} kB)")
    return OUT


if __name__ == "__main__":
    main()
