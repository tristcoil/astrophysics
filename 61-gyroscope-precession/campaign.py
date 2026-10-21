"""All the orbit-fitting runs of the project, in priority order, resumable.

    python campaign.py            # everything not yet cached
    python campaign.py acs        # only the analysis-centre comparison

Each configuration re-fits every weekly orbit of LAGEOS and LARES 2 (about 437
arcs, ~10 s each) into its own cache folder (``data/cache/arcs/<sat><tag>/<ac>``):

* pass 1 (tag ""): ASI orbits, GOCO06s as published; its LAGEOS - LARES 2
  difference measures the Earth's J2 week by week (framedrag.j2_from_difference);
* pass 2 (tag "_c20"): ASI, with that J2 correction, C_R fitted weekly;
* baseline (tag "_c20_cr"): every analysis centre, J2 corrected, the radiation-pressure
  scale C_R held at one value per satellite (a weekly C_R absorbs part of frame
  dragging, see global_cr());
* systematics on ASI's orbits: frame dragging switched on in the model (injection
  tests: mu must drop by exactly the gain), C_R at the nominal ILRS values, CSR's
  monthly C20 instead of our own J2, the GGM05S gravity field, no ocean tides, no de Sitter;
* pass 1 with the 18.6-year equilibrium ocean tide in phase with the solid tide
  (tag "_om1"), to check the J2 the pair measures against CSR's.
"""

from __future__ import annotations

import json
import sys

import numpy as np

import config as C
import framedrag as FD
import nodes
import gravity
import frames as F

J2_TABLE = C.CACHE_DIR / "dc20_pass1_asi.json"


def j2_table() -> tuple:
    """Correction to the normalised C20 measured by the pair (ASI, pass 1), cached."""
    if not J2_TABLE.exists():
        ser = FD.series(C.PRIMARY_AC, "")
        j = FD.j2_from_difference(ser)
        mjd = np.concatenate([[j["mjd"][0] - 30], j["mjd"], [j["mjd"][-1] + 30]])
        dc = np.concatenate([[j["dcbar20"][0]], j["dcbar20"], [j["dcbar20"][-1]]])
        J2_TABLE.write_text(json.dumps({"mjd": mjd.tolist(), "dcbar20": dc.tolist(),
                                        "beta": j["beta"].tolist(),
                                        "beta_smooth": j["beta_smooth"].tolist(),
                                        "mjd_raw": j["mjd"].tolist()}))
    t = json.loads(J2_TABLE.read_text())
    return (t["mjd"], t["dcbar20"])


def csr_table() -> tuple:
    """CSR's monthly SLR C20 expressed as a correction to GOCO06s (time-variable model)."""
    c = np.loadtxt(C.GRAVITY_DIR / "C20_RL06.txt", comments="#", usecols=(0, 1))
    mjd = 51_544.5 + (c[:, 0] - 2000.0) * 365.25
    sel = (mjd > 59_500) & (mjd < 61_400)
    goco = np.array([gravity.load(C.GRAVITY_MODEL, 2, mjd=m).cbar[2, 0] for m in mjd[sel]])
    return (mjd[sel].tolist(), (c[sel, 1] - goco).tolist())


def global_cr() -> dict:
    """One radiation-pressure scale per satellite: the median of ASI's weekly free fits.

    Fitting C_R every week absorbs part of frame dragging (in eclipse seasons the Earth's
    shadow lets radiation pressure turn the orbital plane): on synthetic orbits the weekly
    jumps came out at 81-100% of the true frame dragging with C_R free and 99.7-99.9% with it
    fixed (test_physics.test_jump_gain). So the adopted fits hold C_R at one value per satellite.
    """
    path = C.CACHE_DIR / "global_cr.json"
    if not path.exists():
        med = {}
        for sat in ("lageos1", "lares2"):
            med[sat] = float(np.median([a["prm"][0] for a in nodes.load_arcs(sat, C.PRIMARY_AC, "_c20")]))
        path.write_text(json.dumps(med, indent=1))
    return json.loads(path.read_text())


def configurations() -> list:
    """(name, settings), in priority order. Tags ending in _cr hold C_R at the global value."""
    acs = [C.PRIMARY_AC] + [a for a in C.ACS if a != C.PRIMARY_AC]
    asi = (C.PRIMARY_AC,)
    own = {"dc20_table": "own"}
    fix = {"fix_cr": True, "cr": "global"}
    out = [("pass1", dict(acs=asi, tag="", model_kw={})),
           ("free:ac:asi", dict(acs=asi, tag="_c20", model_kw=own)),
           ("ac:asi", dict(acs=asi, tag="_c20_cr", model_kw=own, fit_kw=fix)),
           ("free:injection", dict(acs=asi, tag="_c20_mu1", model_kw={**own, "mu_lt": 1.0})),
           ("crnominal", dict(acs=asi, tag="_c20_crfix", model_kw=own, fit_kw={"fix_cr": True})),
           ("injection", dict(acs=asi, tag="_c20_cr_mu1", model_kw={**own, "mu_lt": 1.0}, fit_kw=fix))]
    for ac in acs[1:]:
        out.append((f"ac:{ac}", dict(acs=(ac,), tag="_c20_cr", model_kw=own, fit_kw=fix)))
    out += [
        ("csr", dict(acs=asi, tag="_csr_cr", model_kw={"dc20_table": "csr"}, fit_kw=fix)),
        ("ggm05s", dict(acs=asi, tag="_c20_cr_ggm05s", model_kw={**own, "gravity_model": "GGM05S"}, fit_kw=fix)),
        ("noocean", dict(acs=asi, tag="_c20_cr_noocean", model_kw={**own, "ocean": False}, fit_kw=fix)),
        ("nods", dict(acs=asi, tag="_c20_cr_nods", model_kw={**own, "mu_ds": 0.0}, fit_kw=fix)),
        # pass 1 again with the 18.6-yr ocean tide in phase with the solid tide (J2 check)
        ("om1", dict(acs=asi, tag="_om1", model_kw={"om1_flip": True})),
        # the J2 comparison with C_R held fixed (a weekly C_R also absorbs part of a J2 error)
        ("pass1cr", dict(acs=asi, tag="_cr", model_kw={}, fit_kw=fix)),
        ("om1cr", dict(acs=asi, tag="_om1_cr", model_kw={"om1_flip": True}, fit_kw=fix)),
    ]
    return out


def run(only: str | None = None, workers: int = 8) -> None:
    for name, cfg in configurations():
        if only and not name.startswith(only):
            continue
        kw = dict(cfg["model_kw"])
        if kw.get("dc20_table") == "own":
            kw["dc20_table"] = j2_table()
        elif kw.get("dc20_table") == "csr":
            kw["dc20_table"] = csr_table()
        fkw = dict(cfg.get("fit_kw") or {})
        if fkw.get("cr") == "global":
            fkw["cr"] = global_cr()
        print(f"== {name}")
        nodes.run_all(sats=("lageos1", "lares2"), acs=cfg["acs"], model_kw=kw, tag=cfg["tag"],
                      workers=workers, fit_kw=fkw or None)


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else None)
