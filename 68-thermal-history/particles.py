"""The particle content of the hot universe, from the PDG tables shipped with `particle`.

Two lists:

* ``standard_model()``: the 17 elementary particles (plus antiparticles, spins and colours),
  with their PDG masses. g counts every internal state: a quark flavour has 3 colours x 2 spins x
  (quark + antiquark) = 12, the gluons 8 colours x 2 polarisations = 16, the W 2 charges x 3
  polarisations = 6, a neutrino flavour 1 helicity x (nu + anti-nu) = 2.
* ``hadrons()``: every hadron in the PDG table below a mass cut, each charge state and
  antiparticle with its 2J+1 spin states: the hadron-resonance gas.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache

import config as C
from thermo import Species

#: Pole masses of the charm and bottom quarks (PDG quark-mass review: 1.67 +- 0.07 and
#: 4.78 +- 0.06 GeV), an alternative to the PDG's MS-bar masses for a free-gas count.
POLE_MASS = {4: 1.67, 5: 4.78}


def _pdg():
    from particle import Particle

    return Particle


def pdg_mass(pdgid: int) -> float:
    """PDG mass in GeV."""
    return _pdg().from_pdgid(pdgid).mass / 1000.0


@lru_cache(maxsize=None)
def standard_model(quark_masses: str = "pdg") -> tuple[Species, ...]:
    """The Standard Model's species. ``quark_masses``: "pdg" (MS-bar for c, b) or "pole"."""
    sp = []
    quark_names = {1: ("down", "d"), 2: ("up", "u"), 3: ("strange", "s"),
                   4: ("charm", "c"), 5: ("bottom", "b"), 6: ("top", "t")}
    for pid, (name, label) in quark_names.items():
        m = POLE_MASS[pid] if (quark_masses == "pole" and pid in POLE_MASS) else pdg_mass(pid)
        sp.append(Species(name, label, "quark", m, 12, -1, pid))
    sp.append(Species("gluon", "g", "gluon", 0.0, 16, +1, 21))
    sp.append(Species("photon", "γ", "photon", 0.0, 2, +1, 22))
    for pid, name, label in ((11, "electron", "e"), (13, "muon", "μ"), (15, "tau", "τ")):
        sp.append(Species(name, label, "lepton", pdg_mass(pid), 4, -1, pid))
    for pid, name, label in ((12, "nu_e", "νe"), (14, "nu_mu", "νμ"), (16, "nu_tau", "ντ")):
        sp.append(Species(name, label, "neutrino", 0.0, 2, -1, pid))
    sp.append(Species("W", "W", "gauge", pdg_mass(24), 6, +1, 24))
    sp.append(Species("Z", "Z", "gauge", pdg_mass(23), 3, +1, 23))
    sp.append(Species("higgs", "H", "higgs", pdg_mass(25), 1, +1, 25))
    return tuple(sp)


def by_name(species, name: str) -> Species:
    for s in species:
        if s.name == name:
            return s
    raise KeyError(name)


def _spin(p) -> float | None:
    j = p.J
    if j is None:
        return None
    return float(j)


@lru_cache(maxsize=None)
def hadrons(mass_max: float = C.HRG_MASS_MAX, keep_broad_scalars: bool = False,
            statuses: tuple = ("Common", "Rare", "Unsure")) -> tuple[Species, ...]:
    """Every PDG hadron below ``mass_max`` (GeV), one Species per charge state and antiparticle.

    K_S and K_L are dropped (K0 and anti-K0 are already listed). f0(500) and K0*(700) are left
    out unless ``keep_broad_scalars``. States without a mass or spin in the table are skipped.
    """
    P = _pdg()
    out = []
    for p in P.all():
        pid = int(p.pdgid)
        if not p.pdgid.is_hadron or abs(pid) >= 1_000_000_000:
            continue
        if pid in C.HRG_DUPLICATES:
            continue
        if not keep_broad_scalars and pid in C.HRG_EXCLUDE:
            continue
        if p.status.name not in statuses:
            continue
        j = _spin(p)
        if p.mass is None or j is None:
            continue
        m = p.mass / 1000.0
        if m > mass_max:
            continue
        stat = -1 if abs(j - round(j)) > 0.25 else +1
        out.append(Species(p.name, p.name, "baryon" if p.pdgid.is_baryon else "meson",
                           m, 2 * j + 1, stat, pid))
    return tuple(sorted(out, key=lambda s: (s.mass, s.pdgid)))


def hadron_summary(hs) -> dict:
    """Counts used in the README: particles, states, mesons, baryons, the lightest few."""
    n_states = sum(s.g for s in hs)
    c = Counter(s.group for s in hs)
    return {"entries": len(hs), "states": n_states, "mesons": c["meson"], "baryons": c["baryon"],
            "states_below_1GeV": sum(s.g for s in hs if s.mass < 1.0),
            "states_below_2GeV": sum(s.g for s in hs if s.mass < 2.0)}


def hagedorn_spectrum(hs, m_grid):
    """Cumulative number of hadron states (spin x charge x antiparticle) lighter than m."""
    import numpy as np

    masses = np.array([s.mass for s in hs])
    g = np.array([s.g for s in hs])
    order = np.argsort(masses)
    cum = np.cumsum(g[order])
    idx = np.searchsorted(masses[order], m_grid, side="right")
    return np.where(idx > 0, cum[np.maximum(idx - 1, 0)], 0.0)
