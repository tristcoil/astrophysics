"""The interactive page: hide an antimatter domain from the MeV sky.

    python explorer.py

Writes ``figures/antimatter_domain_explorer.html`` (a complete offline page)
and ``data/derived/antimatter_domain_explorer_fragment.html`` (the same page
without the document skeleton, for publishing as an Artifact). All data are
embedded: predicted spectra for 3 heating models x 4 end-of-annihilation
epochs, the COMPTEL, EGRET, Fermi-LAT and Swift-BAT measurements.
"""

from __future__ import annotations

import json
import pickle

import numpy as np

import config as cfg
import domains as dm
import fetch_data as fd
import literature as lit

TEMPLATE = cfg.PROJECT_DIR / "explorer_template.html"
OUT = cfg.FIGURE_DIR / "antimatter_domain_explorer.html"
FRAGMENT = cfg.DERIVED_DIR / "antimatter_domain_explorer_fragment.html"


def _band_average(E, F, lo, hi):
    out = []
    for a, b in zip(lo, hi):
        g = np.geomspace(a, b, 60)
        out.append(np.trapezoid(np.exp(np.interp(np.log(g), np.log(E), np.log(np.maximum(F, 1e-300)))), g)
                   / (b - a))
    return np.array(out)


def build_data() -> dict:
    u = dm.planck2018()
    ph = dm.AnnihilationPhotons(lit.PBARP_CHANNELS, n_events=40000)
    norm = lit.PBARP_MULTIPLICITY["n_gamma"][0] / ph.n_gamma

    class Photons:
        def __call__(self, e):
            return norm * ph(e)

    runs = {k: pickle.load(open(cfg.DERIVED_DIR / "domains" / f"sol_{name}.pkl", "rb"))
            for k, name in (("none", "planck_no_heating"), ("cdg", "planck_cdg_heating"),
                            ("full", "planck_full_heating"))}
    E = np.geomspace(0.05, 120, 110)
    comptel = {"w99": lit.comptel("weidenspointner"), "k98": lit.comptel("kappadath")}
    spectra, ratio1, frac1, bound = {}, {}, {}, {}
    for heat, sol in runs.items():
        yy_all = np.geomspace(10, 1100, 3000)
        jj = np.exp(np.interp(np.log(yy_all), np.log(sol.y[::-1]), np.log(sol.J[::-1])))
        for ymin in (10, 20, 30, 50):
            key = f"{heat}|{ymin}"
            F = dm.diffuse_spectrum(u, sol.y, sol.J, Photons(), E, d0_mpc=1000.0, y_min=float(ymin))
            spectra[key] = (E**2 * F).tolist()
            ratio1[key] = {}
            for a, c in comptel.items():
                pred = _band_average(E, F, c["E_lo"], c["E_hi"])
                ratio1[key][a] = (pred / (c["I"] + 2 * c["err"])).tolist()
            bound[key] = max(ratio1[key]["w99"])
            m = yy_all >= ymin
            frac1[key] = float(np.trapezoid(jj[m] / (cfg.GPC * yy_all[m] ** 3 * u.hubble(yy_all[m]) * u.n_b(1.0)),
                                            yy_all[m]))
    fl = {k: fd.load_fermi_igrb(k) for k in "ABC"}
    fe = fl["A"]["E_mid"]
    lo = np.min([fl[k]["igrb"] - fl[k]["igrb_err_lo"] for k in "ABC"], axis=0)
    hi = np.max([fl[k]["igrb"] + fl[k]["igrb_err_up"] for k in "ABC"], axis=0)
    eg = lit.egret()
    ok = ~eg["is_ul"]
    bat = lit.swift_bat()
    dist = {str(y): u.comoving_distance(float(y)) / cfg.GPC for y in (10, 20, 30, 50)}
    return {
        "E": E.tolist(), "spectra": spectra, "ratio1": ratio1, "frac1": frac1, "bound": bound,
        "comptel": {a: {"lo": c["E_lo"].tolist(), "hi": c["E_hi"].tolist(), "mid": c["E_mid"].tolist(),
                        "I": c["I"].tolist(), "err": c["err"].tolist()} for a, c in comptel.items()},
        "egret": {"mid": eg["E_mid"][ok].tolist(), "I": eg["I"][ok].tolist(), "err": eg["err"][ok].tolist()},
        "fermi": {"mid": fe.tolist(), "lo": (fe**2 * np.maximum(lo, 1e-30)).tolist(), "hi": (fe**2 * hi).tolist(),
                  "elo": fl["A"]["E_lo"].tolist(), "ehi": fl["A"]["E_hi"].tolist()},
        "bat": {"mid": bat["E_mid"].tolist(), "E2I": bat["E2I"].tolist(), "err": bat["E2I_err"].tolist()},
        "dist_cmb": u.comoving_distance(1100.0) / cfg.GPC, "dist_zend": dist,
    }


def main() -> None:
    data = build_data()
    html = TEMPLATE.read_text().replace("/*DATA*/null", json.dumps(data, separators=(",", ":")))
    FRAGMENT.write_text(html)
    head = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
            '</head>\n<body>\n')
    OUT.write_text(head + html + "\n</body>\n</html>\n")
    print("   ", OUT, f"({OUT.stat().st_size / 1e3:.0f} kB)")


if __name__ == "__main__":
    main()
