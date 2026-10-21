"""The hadron list: every known hadron (and light nucleus) with its mass,
width, spin degeneracy, quantum numbers and decays.

The lists are Thermal-FIST's text files (PDG 2025 by default, PDG 2014 as a
variant), which tabulate the Particle Data Group's established states. We
read only the particles; antiparticles are generated here, and so are the
decays of antiparticles (the charge-conjugate channels).

Feed-down is a matrix ``F``: ``F[i, j]`` is how many particles ``j`` come out,
on average, at the end of the decay chain started by one particle ``i``
(for a particle that does not decay, ``F[i, i] = 1``). Final yields are then
``F.T @ primary``. Which particles count as "decaying" is a convention that
must match the experiment (see ``FEEDDOWN``):

* ``"strong"``: strong and electromagnetic decays happen; weakly decaying
  hadrons (K0S, Lambda, Sigma+-, Xi, Omega) are kept. This is ALICE's
  definition of a primary particle: its proton yield excludes protons from
  Lambda decays, its Lambda yield includes Sigma0 -> Lambda gamma.
* ``"weak"``: the strange baryons' weak decays happen too. STAR's beam-energy
  scan protons are inclusive (Lambda -> p pi- included).
* ``"none"``: primordial, before any decay.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np

import config as cfg

#: Weakly decaying hadrons whose decays a "weak" feed-down convention adds.
WEAK_HYPERONS = (3122, 3222, 3112, 3322, 3312, 3334)


@dataclass
class HadronList:
    name: str
    pdg: np.ndarray            # PDG codes, antiparticles negative
    names: list
    mass: np.ndarray           # GeV
    width: np.ndarray          # GeV
    threshold: np.ndarray      # GeV (lowest decay threshold)
    g: np.ndarray              # spin (x isospin-member) degeneracy
    stat: np.ndarray           # +1 fermion, -1 boson
    B: np.ndarray
    Q: np.ndarray
    S: np.ndarray
    C: np.ndarray
    absS: np.ndarray           # number of strange + antistrange valence quarks
    stable: np.ndarray         # True if the list marks it stable (no strong/EM decay)
    decays: dict = field(default_factory=dict)  # pdg -> [(BR, [daughter pdgs]), ...]
    index: dict = field(default_factory=dict)   # pdg -> row

    def __len__(self) -> int:
        return len(self.pdg)

    def idx(self, pdg: int) -> int:
        return self.index[pdg]

    @property
    def is_nucleus(self) -> np.ndarray:
        return np.abs(self.pdg) > 1_000_000_000

    def label(self, pdg: int) -> str:
        return self.names[self.index[pdg]]


def _read_list(path: Path) -> list:
    rows = []
    for line in path.read_text().splitlines():
        line = line.split("#")[0].strip()
        if not line:
            continue
        p = line.split()
        rows.append({"pdg": int(p[0]), "name": p[1], "stable": int(p[2]) == 1,
                     "mass": float(p[3]), "g": float(p[4]), "stat": int(p[5]),
                     "B": int(p[6]), "Q": int(p[7]), "S": int(p[8]), "C": int(p[9]),
                     "absS": float(p[10]), "absC": float(p[11]),
                     "width": float(p[12]), "threshold": float(p[13])})
    return rows


def _read_decays(path: Path) -> dict:
    """Thermal-FIST decay file: pdg, number of channels, then BR + daughters."""
    tokens = []
    for line in path.read_text().splitlines():
        line = line.split("#")[0].strip()
        if line:
            tokens.append(line.split())
    decays, i = {}, 0
    while i < len(tokens):
        pdg = int(tokens[i][0])
        n = int(tokens[i + 1][0])
        chans = []
        for row in tokens[i + 2:i + 2 + n]:
            chans.append((float(row[0]), [int(x) for x in row[1:]]))
        decays[pdg] = chans
        i += 2 + n
    return decays


@lru_cache(maxsize=None)
def load(list_name: str = cfg.DEFAULT_LIST, *, with_nuclei: bool = True) -> HadronList:
    """Read one of Thermal-FIST's lists, add antiparticles and their decays."""
    lst_rel, dec_rel = cfg.HADRON_LISTS[list_name]
    rows = _read_list(cfg.LISTS_DIR / lst_rel)
    if not with_nuclei:
        rows = [r for r in rows if abs(r["pdg"]) < 1_000_000_000]
    decays = _read_decays(cfg.LISTS_DIR / dec_rel)

    full = []
    for r in rows:
        if r["g"] <= 0:          # e.g. f0(500): listed with zero degeneracy
            continue
        full.append(r)
        if r["B"] != 0 or r["Q"] != 0 or r["S"] != 0 or r["C"] != 0:
            a = dict(r)
            a["pdg"] = -r["pdg"]
            a["name"] = "anti-" + r["name"]
            for q in ("B", "Q", "S", "C"):
                a[q] = -r[q]
            full.append(a)
    pdgs = [r["pdg"] for r in full]
    self_conj = {r["pdg"] for r in full if -r["pdg"] not in pdgs}

    def conj(d: int) -> int:
        if d in self_conj or d == 22:
            return d
        return -d

    dec_all = {}
    for p in pdgs:
        src = decays.get(abs(p))
        if src is None:
            continue
        if p > 0 or abs(p) in self_conj:
            dec_all[p] = [(br, list(ds)) for br, ds in src]
        else:
            dec_all[p] = [(br, [conj(d) for d in ds]) for br, ds in src]

    arr = lambda k, dt=float: np.array([r[k] for r in full], dtype=dt)  # noqa: E731
    hl = HadronList(
        name=list_name, pdg=arr("pdg", np.int64), names=[r["name"] for r in full],
        mass=arr("mass"), width=arr("width"), threshold=arr("threshold"), g=arr("g"),
        stat=arr("stat", int), B=arr("B", int), Q=arr("Q", int), S=arr("S", int),
        C=arr("C", int), absS=arr("absS"), stable=arr("stable", bool), decays=dec_all)
    hl.index = {int(p): i for i, p in enumerate(hl.pdg)}
    return hl


def decaying_mask(hl: HadronList, convention: str) -> np.ndarray:
    """Which particles decay under a feed-down convention."""
    if convention == "none":
        return np.zeros(len(hl), bool)
    mask = ~hl.stable
    if convention == "strong":
        return mask
    if convention == "weak":
        mask = mask.copy()
        for p in WEAK_HYPERONS:
            for s in (p, -p):
                if s in hl.index:
                    mask[hl.index[s]] = True
        return mask
    raise ValueError(convention)


@lru_cache(maxsize=None)
def feeddown_matrix(list_name: str, convention: str, with_nuclei: bool = True) -> np.ndarray:
    """F[i, j]: mean number of j at the end of i's decay chain (memoised recursion).

    Branching ratios are normalised to one per particle; daughters that are
    not hadrons in the list (photons, leptons) are dropped.
    """
    hl = load(list_name, with_nuclei=with_nuclei)
    n = len(hl)
    dec = decaying_mask(hl, convention)
    rows: dict = {}

    def row(i: int, depth: int = 0) -> np.ndarray:
        if i in rows:
            return rows[i]
        if depth > 50:
            raise RuntimeError(f"decay chain loop at {hl.names[i]}")
        out = np.zeros(n)
        chans = hl.decays.get(int(hl.pdg[i]))
        if not dec[i] or not chans:
            out[i] = 1.0
        else:
            total = sum(br for br, _ in chans)
            for br, ds in chans:
                for d in ds:
                    j = hl.index.get(d)
                    if j is not None:
                        out += (br / total) * row(j, depth + 1)
        rows[i] = out
        return out

    return np.array([row(i) for i in range(n)])


@lru_cache(maxsize=None)
def feeddown_into(list_name: str, convention: str, target: int,
                  with_nuclei: bool = True) -> np.ndarray:
    """Vector v[i]: how many ``target`` particles the decay chain of one ``i``
    ends with, with the target itself kept intact even if it can decay (so a
    measured phi or K* counts phis from heavier resonances, but is not
    decayed into kaons). For a stable target this is F[:, target]."""
    hl = load(list_name, with_nuclei=with_nuclei)
    t = hl.index[target]
    if not decaying_mask(hl, convention)[t]:
        return feeddown_matrix(list_name, convention, with_nuclei)[:, t].copy()
    dec = decaying_mask(hl, convention).copy()
    dec[t] = False
    memo: dict = {}

    def count(i: int, depth: int = 0) -> float:
        if i == t:
            return 1.0
        if i in memo:
            return memo[i]
        if depth > 50:
            raise RuntimeError("decay chain loop")
        chans = hl.decays.get(int(hl.pdg[i]))
        val = 0.0
        if dec[i] and chans:
            total = sum(br for br, _ in chans)
            for br, ds in chans:
                for d in ds:
                    j = hl.index.get(d)
                    if j is not None:
                        val += (br / total) * count(j, depth + 1)
        memo[i] = val
        return val

    return np.array([count(i) for i in range(len(hl))])


@lru_cache(maxsize=None)
def charm_feeddown_per_pair(T: float = 0.1565) -> dict:
    """Light hadrons from the decays of one charm-anticharm pair.

    Charm quarks are made in the first hard collisions, not thermally, but
    they hadronise like everything else: the D, D*, D_s, Lambda_c ... are
    populated in thermal proportions (Boltzmann weights at T), decay
    strongly to the ground states and then weakly, within ~0.1 mm of the
    vertex, so ALICE counts their pions and kaons as primary.
    Returns {pdg: mean multiplicity per c-cbar pair} for pi, K, K0, p.
    """
    from scipy.special import kn

    hl = load("PDG2025charm", with_nuclei=False)
    dec = (~hl.stable) | (hl.C != 0)
    n = len(hl)
    memo: dict = {}

    def row(i: int) -> np.ndarray:
        if i in memo:
            return memo[i]
        out = np.zeros(n)
        chans = hl.decays.get(int(hl.pdg[i]))
        if not dec[i] or not chans:
            out[i] = 1.0
        else:
            total = sum(br for br, _ in chans)
            for br, ds in chans:
                for d in ds:
                    j = hl.index.get(d)
                    if j is not None:
                        out += (br / total) * row(j)
        memo[i] = out
        return out

    w = np.zeros(n)
    for i in np.flatnonzero(np.abs(hl.C) == 1):
        m = hl.mass[i]
        w[i] = hl.g[i] * m * m * kn(2, m / T)
    pair = np.zeros(n)
    for sign in (1, -1):                      # one charm and one anticharm quark
        ws = np.where(hl.C == sign, w, 0.0)
        pair += (ws / ws.sum()) @ np.array([row(i) for i in range(n)])
    return {p: float(pair[hl.index[p]]) for p in (211, -211, 321, -321, 311, -311, 2212, -2212)}


def branching_ratio_sums(list_name: str = cfg.DEFAULT_LIST) -> dict:
    """Sum of listed branching ratios per unstable particle (should be ~1)."""
    hl = load(list_name)
    out = {}
    for i in np.flatnonzero(~hl.stable):
        chans = hl.decays.get(int(hl.pdg[i]))
        out[hl.names[i]] = sum(br for br, _ in chans) if chans else 0.0
    return out


def conserved_in_decays(list_name: str = cfg.DEFAULT_LIST) -> list:
    """Channels that violate B, Q or S (none should, except weak decays changing S)."""
    hl = load(list_name)
    bad = []
    for p, chans in hl.decays.items():
        i = hl.index[p]
        for br, ds in chans:
            js = [hl.index[d] for d in ds if d in hl.index]
            leptons_q = sum({11: -1, -11: 1, 13: -1, -13: 1, 15: -1, -15: 1}.get(d, 0) for d in ds)
            dB = hl.B[i] - sum(hl.B[j] for j in js)
            dQ = hl.Q[i] - sum(hl.Q[j] for j in js) - leptons_q
            if dB != 0 or dQ != 0:
                bad.append((hl.names[i], br, ds, int(dB), int(dQ)))
    return bad
