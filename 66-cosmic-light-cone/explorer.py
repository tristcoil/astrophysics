"""Build the interactive light-cone explorer: one self-contained HTML file, no libraries.

    python explorer.py      # writes figures/light_cone_explorer.html (~100 kB)

It also writes data/derived/light_cone_explorer_artifact.html, the same page as
a fragment for publishing in a viewer that supplies the document skeleton and
the light/dark theme (so it has no theme buttons of its own).

The page carries the Planck 2018 expansion history on a grid (t, a, eta, H) and
the objects from results.json; everything else (light cones seen at any
moment, horizons, the photon's journey) is computed in the browser from that
grid, so it works offline.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # works under python -I too

import numpy as np  # noqa: E402

import config as C  # noqa: E402
from analysis import object_row  # noqa: E402
from cosmology import Universe  # noqa: E402

TEMPLATE = C.HERE / "explorer_template.html"
OUT = C.FIG_DIR / "light_cone_explorer.html"
OUT_ARTIFACT = C.DERIVED_DIR / "light_cone_explorer_artifact.html"
THEME_BUTTONS = """      <div class="seg themebtn" role="group" aria-label="Theme">
        <button type="button" data-theme-set="light" aria-pressed="false">Light</button>
        <button type="button" data-theme-set="dark" aria-pressed="false">Dark</button>
      </div>"""

KIND_LABEL = {"cmb": "cosmic microwave background", "galaxy": "galaxy", "quasar": "quasar",
              "grb": "gamma-ray burst", "supernova": "supernova", "gw": "gravitational-wave merger",
              "hypothetical": "an example galaxy (not a catalogued object)"}

SHOW = ["CMB (last scattering)", "MoM-z14", "JADES-GS-z14-0", "GN-z11", "GRB 090423", "J0313-1806",
        "GRB 250314A + SN", "SN 2023adsy", "SN colfax", "GW190403_051519", "GW150914", "GW170817"]


def build(results: dict | None = None) -> Path:
    if results is None:
        results = json.loads(C.RESULTS_PATH.read_text())
    u = Universe(C.PLANCK18)
    t = np.geomspace(1e-6, 300.0, 1500)
    grid = {"t": t, "a": u.a_of_t(t), "eta": u.eta_of_t(t), "H": u.H_of_t(t)}
    grid = {k: [float(f"{x:.7g}") for x in v] for k, v in grid.items()}
    marks = {m["name"]: m for m in results["cone"]["landmarks"]}
    objects = []
    for name in SHOW:
        m = marks.get(name)
        if m is None:
            continue
        nice = {"CMB (last scattering)": "CMB (gas at last scattering)", "SN colfax": "SN Colfax",
                "GW190403_051519": "GW190403 (farthest GW merger)",
                "GRB 250314A + SN": "GRB 250314A and its supernova"}.get(name, name)
        objects.append({"name": nice, "kind_label": KIND_LABEL.get(m["kind"], m["kind"]), "ref": m["ref"],
                        "note": m.get("note", ""), "z": m["z"], "z_label": f"{m['z']:.4g}",
                        "chi": m["chi_gly"], "t_emit": m["t_emit_gyr"], "d_emit": m["d_emit_gly"],
                        "v_then": m["v_then_c"], "v_now": m["v_now_c"], "reachable": m["reachable"],
                        "t_last": m["t_last_seen_gyr"] if isinstance(m["t_last_seen_gyr"], float) else None})
    s = results["headline"]["summary"]
    for z, label, note in [
            (s["z_hubble_sphere"], "a galaxy receding at exactly c today",
             "It sits on today's Hubble sphere, c/H0 = 14.45 billion ly away."),
            (s["z_light_cone_meets_hubble_sphere"], "a galaxy at the light cone's widest point",
             "Its light set out exactly on the Hubble sphere, so the photon hovered at a fixed distance "
             "before approaching. Anything farther emitted light that first moved away from us."),
            (s["z_event_horizon"], "a galaxy crossing our event horizon today",
             "Light it emits from now on will never reach us; a signal we send it now will never arrive."),
            (1.0, "a galaxy at z = 1", "")]:
        r = object_row(u, z, name=label, kind="hypothetical", ref="")
        objects.append({"name": label, "kind_label": KIND_LABEL["hypothetical"], "ref": "",
                        "note": note, "z": z, "z_label": f"{z:.4g}", "chi": r["chi_gly"],
                        "t_emit": r["t_emit_gyr"], "d_emit": r["d_emit_gly"], "v_then": r["v_then_c"],
                        "v_now": r["v_now_c"], "reachable": r["reachable"], "t_last": r["t_last_seen_gyr"]})
    objects.sort(key=lambda o: -o["z"])
    default = next(i for i, o in enumerate(objects) if o["name"] == "GN-z11")
    data = {
        "grid": grid, "eta_inf": u.eta_inf, "age": u.age, "eta0": u.eta0,
        "worldlines": [{"z": z, "chi": float(u.chi_of_z(z))} for z in (1, 3, 10, 1090)],
        "objects": objects, "default_index": default,
        "footnote": ("Planck 2018 cosmology (H0 = 67.66 km/s/Mpc, Ωm = 0.311, flat, with photons, neutrinos "
                     "and one 0.06 eV neutrino mass), integrated by this project's own code (checked against "
                     "astropy to 2e-6 and CAMB to 1.5e-5). Distances in billions of light-years; 1 Gly per Gyr "
                     "is the speed of light. Redshifts: JADES DR4, Naidu et al. 2026 (MoM-z14), Wang et al. "
                     "2021 (J0313-1806), Tanvir et al. 2009 (GRB 090423), Levan et al. 2025 (GRB 250314A), "
                     "Pierel et al. 2024 (SN 2023adsy), Pantheon+ (SN Colfax), GWOSC (mergers; redshift from "
                     "the measured distance). Part of project O2 of the astrophysics atlas: "
                     "66-cosmic-light-cone/."),
    }
    fragment = TEMPLATE.read_text().replace("__DATA__", json.dumps(data, separators=(",", ":")))
    # the standalone page: own skeleton and theme buttons
    page = fragment.replace("<!--THEME_BUTTONS-->", THEME_BUTTONS)
    head, body = page.split("<main>", 1)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
                   '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
                   + head + "</head>\n<body>\n<main>" + body + "</body>\n</html>\n")
    OUT_ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    OUT_ARTIFACT.write_text(fragment.replace("<!--THEME_BUTTONS-->\n", "").replace("<!--THEME_BUTTONS-->", ""))
    return OUT


if __name__ == "__main__":
    p = build()
    print(p, f"{p.stat().st_size / 1e3:.0f} kB")
