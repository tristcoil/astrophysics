"""Our own four-lepton selection, run on the good-lepton skims of ``events.py``.

It follows the CMS 2012 discovery analysis (Phys. Lett. B 716, 30) more
closely than the outreach code does, and fixes what we found in the latter:

* **Any number of leptons.** The outreach code demands *exactly* four muons
  (or four electrons) and rejects the event if any lepton fails a cut; here
  every event with at least four good leptons is considered.
* **Z1 first.** Among all opposite-charge, same-flavour pairs, Z1 is the one
  closest to the Z mass (91.1876 GeV, compared without truncation), with
  40 < m(Z1) < 120 GeV; Z2 is the remaining such pair with the largest summed
  transverse momentum, 12 < m(Z2) < 120 GeV.
* **Kinematics.** One lepton above 20 GeV and another above 10 GeV; all
  leptons separated by Delta R > 0.02; every opposite-charge pair heavier than
  4 GeV (no J/psi or Upsilon decays); electrons within Delta R < 0.05 of a
  muon are dropped (a muon track reconstructed twice).
* **No double counting.** 4mu events are taken from the DoubleMuParked data,
  4e from DoubleElectron, 2e2mu from both with duplicates (same run,
  luminosity block and event number) removed.
* **A control region** for backgrounds the outreach code neglects (Z + jets,
  top pairs: events where some "leptons" are misidentified or come from
  heavy-quark decays): Z2 built from a *same-sign* pair, which no real ZZ or
  Higgs decay produces.

Lepton momenta can be rescaled before selection (``scale``), for the
momentum-scale calibration of ``calibration.py``.
"""

from __future__ import annotations

import math
from collections import defaultdict

import numpy as np

import config as cfg
import events

MZ_PDG = 91.1876
CHANNELS = ("4mu", "4e", "2e2mu")


def _p4(pt, eta, phi, m):
    px, py, pz = pt * math.cos(phi), pt * math.sin(phi), pt * math.sinh(eta)
    return (math.sqrt(px * px + py * py + pz * pz + m * m), px, py, pz)


def _m(*vs):
    e = sum(v[0] for v in vs)
    px = sum(v[1] for v in vs)
    py = sum(v[2] for v in vs)
    pz = sum(v[3] for v in vs)
    m2 = e * e - px * px - py * py - pz * pz
    return math.copysign(math.sqrt(abs(m2)), m2)


def _dr(a, b):
    dphi = (b[2] - a[2] + math.pi) % (2 * math.pi) - math.pi
    return math.hypot(b[1] - a[1], dphi)


def _event_leptons(skim: dict, scale=None):
    """Yield (index, leptons) per skimmed event; a lepton is
    (flavour, pt, eta, phi, mass, charge)."""
    n_mu, n_e = skim["n_mu"], skim["n_e"]
    off_mu = np.concatenate([[0], np.cumsum(n_mu)])
    off_e = np.concatenate([[0], np.cumsum(n_e)])
    cols = {f: skim[f].tolist() for f in skim if f.startswith(("mu_", "e_"))}
    for i in range(len(n_mu)):
        leps = []
        for flav, off in (("mu", off_mu), ("e", off_e)):
            for k in range(off[i], off[i + 1]):
                pt = cols[f"{flav}_pt"][k]
                eta = cols[f"{flav}_eta"][k]
                if scale is not None:
                    pt *= scale(flav, eta)
                leps.append((flav, pt, eta, cols[f"{flav}_phi"][k], cols[f"{flav}_mass"][k],
                             cols[f"{flav}_charge"][k]))
        yield i, leps


def best_candidate(leps, *, same_sign_z2: bool = False):
    """The ZZ candidate of one event, or None. Returns a dict with masses,
    the channel and the four leptons."""
    # electrons sharing a track with a muon
    mus = [l for l in leps if l[0] == "mu"]
    leps = [l for l in leps if l[0] == "mu" or all(_dr(l, m) >= 0.05 for m in mus)]
    if len(leps) < 4:
        return None
    vec = [_p4(l[1], l[2], l[3], l[4]) for l in leps]
    n = len(leps)
    best = None
    for i in range(n):
        for j in range(i + 1, n):
            if leps[i][0] != leps[j][0] or leps[i][5] == leps[j][5]:
                continue
            m = _m(vec[i], vec[j])
            if 40 < m < 120 and (best is None or abs(m - MZ_PDG) < abs(best[0] - MZ_PDG)):
                best = (m, i, j)
    if best is None:
        return None
    mz1, i1, j1 = best
    z2 = None
    for i in range(n):
        for j in range(i + 1, n):
            if len({i, j, i1, j1}) < 4 or leps[i][0] != leps[j][0]:
                continue
            opposite = leps[i][5] != leps[j][5]
            if opposite == same_sign_z2:
                continue
            m = _m(vec[i], vec[j])
            if 12 < m < 120:
                ptsum = leps[i][1] + leps[j][1]
                if z2 is None or ptsum > z2[0]:
                    z2 = (ptsum, m, i, j)
    if z2 is None:
        return None
    _, mz2, i2, j2 = z2
    four = [i1, j1, i2, j2]
    pts = sorted((leps[k][1] for k in four), reverse=True)
    if not (pts[0] > 20 and pts[1] > 10):
        return None
    for a in range(4):
        for b in range(a + 1, 4):
            la, lb = leps[four[a]], leps[four[b]]
            if _dr(la, lb) < 0.02:
                return None
            if la[5] != lb[5] and _m(vec[four[a]], vec[four[b]]) < 4:
                return None
    nmu = sum(leps[k][0] == "mu" for k in four)
    channel = {4: "4mu", 0: "4e", 2: "2e2mu"}[nmu]
    return {"m4l": _m(*(vec[k] for k in four)), "mz1": mz1, "mz2": mz2, "channel": channel,
            "leptons": [leps[k] for k in four]}


def select(skim: dict, *, scale=None, same_sign_z2: bool = False) -> dict:
    """Run the selection over one skim; returns per-candidate arrays."""
    out = defaultdict(list)
    for i, leps in _event_leptons(skim, scale):
        c = best_candidate(leps, same_sign_z2=same_sign_z2)
        if c is None:
            continue
        out["run"].append(int(skim["run"][i]))
        out["lumi"].append(int(skim["lumi"][i]))
        out["event"].append(int(skim["event"][i]))
        out["m4l"].append(c["m4l"])
        out["mz1"].append(c["mz1"])
        out["mz2"].append(c["mz2"])
        out["channel"].append(c["channel"])
        out["lep_pt"].append([l[1] for l in c["leptons"]])
        out["lep_eta"].append([l[2] for l in c["leptons"]])
        out["lep_flav"].append([l[0] for l in c["leptons"]])
    res = {k: np.array(v) for k, v in out.items()}
    if not res:
        res = {k: np.zeros(0) for k in ("run", "lumi", "event", "m4l", "mz1", "mz2")}
        res["channel"] = np.zeros(0, dtype="<U5")
    return res


def mc_weight(name: str) -> float:
    """The outreach normalisation: luminosity x cross-section / generated events
    (x 1.386 for ZZ)."""
    w = cfg.LUMINOSITY_PB * cfg.XSEC_PB[name] / cfg.N_GENERATED[name]
    return w * (cfg.ZZ_SCALE_FACTOR if name != "SMHiggsToZZTo4L" else 1.0)


#: Which simulated sample describes each channel (as in the outreach analysis).
MC_FOR_CHANNEL = {"4mu": "ZZTo4mu", "4e": "ZZTo4e", "2e2mu": "ZZTo2e2mu"}
#: Which data stream each channel is taken from (trigger logic).
DATA_FOR_CHANNEL = {"4mu": ("Run2012B_DoubleMuParked", "Run2012C_DoubleMuParked"),
                    "4e": ("Run2012B_DoubleElectron", "Run2012C_DoubleElectron"),
                    "2e2mu": cfg.DATA_4L}


def data_candidates(*, scale=None, same_sign_z2: bool = False) -> dict:
    """Selected data events, each collision once."""
    per_file = {n: select(events.load_skim(n), scale=scale, same_sign_z2=same_sign_z2)
                for n in cfg.DATA_4L}
    rows, seen, dup = [], set(), 0
    for ch in CHANNELS:
        for n in DATA_FOR_CHANNEL[ch]:
            r = per_file[n]
            for k in np.nonzero(r["channel"] == ch)[0]:
                key = (r["run"][k], r["lumi"][k], r["event"][k])
                if key in seen:
                    dup += 1
                    continue
                seen.add(key)
                rows.append({f: r[f][k] for f in r})
                rows[-1]["dataset"] = n
    out = {f: np.array([row[f] for row in rows]) for f in (rows[0] if rows else {})}
    out["n_duplicates_removed"] = dup
    return out


def mc_candidates(*, scale=None, same_sign_z2: bool = False) -> dict:
    """{'H': {...}, 'ZZ': {...}} with per-candidate weights, channel-matched."""
    res = {}
    h = select(events.load_skim("SMHiggsToZZTo4L"), scale=scale, same_sign_z2=same_sign_z2)
    h["weight"] = np.full(len(h["m4l"]), mc_weight("SMHiggsToZZTo4L"))
    res["H"] = h
    parts = []
    for ch, name in MC_FOR_CHANNEL.items():
        r = select(events.load_skim(name), scale=scale, same_sign_z2=same_sign_z2)
        keep = r["channel"] == ch
        parts.append({k: v[keep] for k, v in r.items()} | {
            "weight": np.full(int(keep.sum()), mc_weight(name)), "sample": name})
    res["ZZ"] = {k: np.concatenate([p[k] for p in parts]) for k in ("m4l", "mz1", "mz2",
                                                                  "channel", "weight")}
    return res


def yields(cands: dict, lo: float, hi: float, weight=None) -> dict:
    m = cands["m4l"]
    w = np.ones(len(m)) if weight is None else weight
    out = {}
    for ch in CHANNELS:
        sel = (cands["channel"] == ch) & (m > lo) & (m < hi)
        out[ch] = float(np.sum(w[sel]))
    out["all"] = sum(out[ch] for ch in CHANNELS)
    return out
