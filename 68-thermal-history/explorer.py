"""Build the offline explorer page, figures/cosmic_thermometer.html, from this project's curves.

    python explorer.py

The page (explorer_template.html) is a body fragment; this script injects the data and writes
two copies: a standalone document for opening from disk, and the bare fragment for publishing
as an Artifact (data/derived/cosmic_thermometer_artifact.html).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

import config as C  # noqa: E402

TEMPLATE = C.HERE / "explorer_template.html"
OUT = C.FIG_DIR / "cosmic_thermometer.html"
OUT_ARTIFACT = C.DERIVED_DIR / "cosmic_thermometer_artifact.html"

GROUP_KEYS = [  # (label, parts in the equation of state) - same order and colours as figure 2
    ("Photons", ("photon",)),
    ("Neutrinos", ("nu_e", "nu_mu", "nu_tau")),
    ("Electrons & positrons", ("electron", "qed")),
    ("Muons & taus", ("muon", "tau")),
    ("W, Z & Higgs bosons", ("W", "Z", "higgs")),
    ("Quarks & gluons / hadrons", ("strong",)),
    ("Matter (dark + atoms)", ()),
    ("Dark energy", ()),
]

#: Plain-language notes for each event (shown on the page).
NOTES = {
    "ew_crossover": "The Higgs field settles towards its present value, and W and Z bosons, quarks and charged "
                    "leptons take on their masses. In the Standard Model this is a smooth crossover, not a "
                    "sudden phase transition (lattice result: 159.5 ± 1.5 GeV).",
    "top_annihilation": "Collisions are no longer energetic enough to make top quarks, the heaviest known "
                        "particle (172.6 GeV). The remaining ones annihilate or decay.",
    "W_annihilation": "W and Z bosons and the Higgs boson disappear the same way.",
    "bottom_annihilation": "Bottom quarks annihilate.",
    "tau_annihilation": "Tau leptons annihilate.",
    "charm_annihilation": "Charm quarks annihilate (as free quarks; by now the strong force is reshaping everything).",
    "qcd_crossover": "Quarks and gluons bind into hadrons: protons, neutrons and short-lived pions. Lattice QCD puts "
                     "this crossover at 156.5 ± 1.5 MeV. The number of particle species drops from about 60 to 17, "
                     "and the pressure briefly sags.",
    "pion_annihilation": "Pions decay and annihilate. Of the hadrons only a trace of protons and neutrons is left, "
                         "about one for every billion photons.",
    "muon_annihilation": "Muons annihilate. The universe is now photons, electrons and positrons, three kinds of "
                         "neutrino and that trace of nucleons.",
    "nu_mu_decoupling": "Neutrinos stop sharing energy with everything else (muon and tau neutrinos first, electron "
                        "neutrinos at 1.5 MeV). They still fill space today: 336 per cubic centimetre at 1.95 K.",
    "nu_e_decoupling": "Electron neutrinos decouple too.",
    "np_freeze_out": "Weak reactions can no longer turn neutrons into protons and back fast enough. The ratio stops "
                     "near one neutron for six protons; free neutrons then slowly decay.",
    "electron_annihilation": "Electrons and positrons annihilate. Their energy heats the photons but not the "
                             "neutrinos, which end up 1.4 times colder. One electron in a billion survives.",
    "nuclei_form": "Deuterium nuclei now survive the photons that used to break them apart, and within minutes "
                   "nearly every neutron ends up in helium-4, about a quarter of all the mass in atoms.",
    "equality": "Matter now carries more energy than radiation. From here on, gravity can grow the small "
                "lumps that later become galaxies.",
    "last_scattering": "Electrons and nuclei combine into neutral atoms. The universe becomes transparent and the "
                       "light we now see as the cosmic microwave background sets out.",
    "today": "The photons have cooled to 2.7255 K. The relic neutrinos are at about 1.95 K.",
}

SHORT = {
    "ew_crossover": "The Higgs field switches on", "top_annihilation": "Top quarks disappear",
    "W_annihilation": "W, Z and Higgs bosons disappear", "bottom_annihilation": "Bottom quarks disappear",
    "tau_annihilation": "Tau leptons disappear", "charm_annihilation": "Charm quarks disappear",
    "qcd_crossover": "Quarks bind into protons and neutrons", "pion_annihilation": "Pions disappear",
    "muon_annihilation": "Muons disappear", "nu_mu_decoupling": "Neutrinos decouple",
    "nu_e_decoupling": "Electron neutrinos decouple", "np_freeze_out": "The neutron-to-proton ratio freezes",
    "electron_annihilation": "Electrons and positrons annihilate", "nuclei_form": "The first nuclei form",
    "equality": "Matter outweighs radiation", "last_scattering": "Atoms form; the CMB is released",
    "today": "Today",
}


def _data(res: dict) -> dict:
    d = dict(np.load(C.DERIVED_DIR / "curves.npz", allow_pickle=False))
    cT, ct = d["cal_T"][::-1], d["cal_t"][::-1]  # ascending in time
    n = 900
    lt_grid = np.linspace(np.log10(ct[0]), np.log10(ct[-1]), n)
    lT = np.interp(lt_grid, np.log10(ct), np.log10(cT))
    T = 10 ** lT

    def on_t(arr):
        return np.interp(lt_grid, np.log10(ct), arr[::-1])

    rho_bath, rho_m = on_t(d["cal_rho_bath"]), on_t(d["cal_rho_m"])
    rho_L = res["calendar"]["late_universe"]["rho_L"]
    tot = rho_bath + rho_m + rho_L
    # thermal-bath composition: the equation of state above 10 keV, photons + neutrinos below
    names = [str(x) for x in d["parts_names"]]
    parts = d["parts_g_rho"]
    eT = d["T"]
    gtot = parts.sum(axis=0)
    comp = []
    photon_bath = (np.pi ** 2 / 15 * T ** 4) / rho_bath
    for label, keys in GROUP_KEYS[:6]:
        g = sum(parts[names.index(k)] for k in keys) / gtot
        f = np.interp(np.log(T), np.log(eT), g)
        if label == "Photons":
            f = np.where(T < eT[0], photon_bath, f)
        elif label == "Neutrinos":
            f = np.where(T < eT[0], 1 - photon_bath, f)
        else:
            f = np.where(T < eT[0], 0.0, f)
        comp.append(f * rho_bath / tot)
    comp.append(rho_m / tot)
    comp.append(np.full_like(T, rho_L) / tot)
    comp = np.array(comp)
    comp = comp / comp.sum(axis=0)
    g_rho = np.interp(np.log(T), np.log(d["cal_T"]), d["cal_g_rho"])
    g_s = np.interp(np.log(T), np.log(eT), d["g_s"], left=float(d["g_s"][0]), right=float(d["g_s"][-1]))

    def r(x, nd=4):
        return [float(f"{v:.{nd}g}") for v in x]

    ev = []
    for e in res["calendar"]["events"]:
        if e["key"] in ("W_annihilation", "nu_e_decoupling"):
            continue
        ev.append({"key": e["key"], "title": SHORT.get(e["key"], e["label"]), "note": NOTES.get(e["key"], ""),
                   "t": e["t_s"], "T": e["T_GeV"], "source": e["source"]})
    le = res["calendar"]["late_epochs"]
    t_acc = le["acceleration_begins"]["t_Gyr"] * 1e9 * C.YEAR_S
    T_acc = float(np.interp(np.log10(t_acc), lt_grid, lT))
    ev.append({"key": "acceleration", "title": "The expansion starts to speed up",
               "note": "Dark energy's push overtakes matter's pull (at redshift 0.64); three billion years later "
                       "dark energy also dominates the energy budget.",
               "t": t_acc, "T": 10 ** T_acc, "source": "this work (Planck 2018 densities)"})
    ev.sort(key=lambda e: e["t"])
    labs = [{"key": l["key"], "label": l["label"], "t": l["t_s"], "T": l["T_GeV"]} for l in res["calendar"]["labs"]]
    E = {e["key"]: e["t"] for e in ev}
    eras = [
        (0, E["ew_crossover"], "Before the Higgs field switched on",
         "Every Standard Model particle is a nearly massless member of one hot soup."),
        (E["ew_crossover"], E["qcd_crossover"], "Quark–gluon plasma",
         "Quarks and gluons move freely; heavy particles drop out as the soup cools."),
        (E["qcd_crossover"], E["muon_annihilation"], "Hadrons",
         "Quarks are locked inside hadrons; pions and muons are still plentiful."),
        (E["muon_annihilation"], E["nu_mu_decoupling"], "Leptons",
         "Photons, electrons, positrons and neutrinos, all in balance."),
        (E["nu_mu_decoupling"], E["nuclei_form"], "Neutrinos fly free",
         "Neutrinos have decoupled; electrons and positrons annihilate."),
        (E["nuclei_form"], 1200.0, "Nucleosynthesis", "Deuterium, helium and a little lithium are made."),
        (1200.0, E["equality"], "Radiation era", "A hot, opaque plasma of photons, electrons and nuclei."),
        (E["equality"], E["last_scattering"], "Matter era (still opaque)",
         "Dark matter starts to clump; the plasma is still a fog."),
        (E["last_scattering"], le["matter_lambda_equality"]["t_Gyr"] * 1e9 * C.YEAR_S, "Neutral atoms",
         "Transparent and dark at first; later the first stars and galaxies light up."),
        (le["matter_lambda_equality"]["t_Gyr"] * 1e9 * C.YEAR_S, 1e30, "Dark energy era",
         "Dark energy dominates and the expansion accelerates."),
    ]
    return {
        "lt": r(lt_grid, 6), "lT": r(lT, 6), "g": r(g_rho), "gs": r(g_s),
        "comp": [r(100 * c, 3) for c in comp], "groups": [g[0] for g in GROUP_KEYS],
        "lRH": r(np.log10(on_t(d["cal_RH"])), 5), "lMH": r(np.log10(on_t(d["cal_MH"])), 5),
        "lf": r(np.log10(on_t(d["cal_f"])), 5),
        "events": ev, "labs": labs,
        "eras": [{"t0": a, "t1": b, "name": nm, "desc": ds} for a, b, nm, ds in eras],
        "kB": C.K_B_GEV_PER_K, "T0K": C.T0_K,
        "Neff": res["neutrinos"]["N_eff_ours"], "age": res["calendar"]["age_today_Gyr"],
        "pta": res["calendar"]["pta_band_T_GeV"],
    }


def build(res: dict | None = None) -> Path:
    if res is None:
        res = json.loads(C.RESULTS_PATH.read_text())
    data = _data(res)
    frag = TEMPLATE.read_text().replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":")))
    C.DERIVED_DIR.mkdir(parents=True, exist_ok=True)
    OUT_ARTIFACT.write_text(frag)
    head, body = frag.split("<!--BODY-->", 1)
    doc = ("<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
           "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, viewport-fit=cover\">\n"
           + head + "</head>\n<body>\n" + body + "\n</body>\n</html>\n")
    OUT.write_text(doc)
    return OUT


if __name__ == "__main__":
    print("wrote", build())
