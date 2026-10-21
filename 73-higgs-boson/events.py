"""Read the CMS events and reduce them to what the analyses need.

Three passes over the ROOT files, each cached in ``data/cache``:

1. **Dimuon pass** (61.5 million events, 2.2 GB): every event with exactly two
   muons of opposite charge, as in the CMS dimuon-spectrum example (record
   12342), kept as compact per-event arrays: the pair's invariant mass, both
   muons' pseudorapidities, the pair's transverse momentum and rapidity, and
   the Collins-Soper decay angle.
2. **Outreach replica** (record 12360, ``skim.cxx``): the four-lepton
   selection exactly as published -- including two quirks of the original
   C++ (see ``pair_closest_to_z``) -- so that our histograms can be compared
   bin by bin with the plots CMS released.
3. **Good-lepton skim**: every event with at least four leptons passing the
   identification cuts (any charges), with those leptons' kinematics, for our
   own selection in ``fourlepton.py``.

Kinematics: a particle with transverse momentum pT, pseudorapidity eta,
azimuth phi and mass m has px = pT cos(phi), py = pT sin(phi),
pz = pT sinh(eta) and E = sqrt(pT^2 cosh^2(eta) + m^2). The invariant mass of
a system is sqrt((sum E)^2 - |sum p|^2); computed here in double precision
from the stored single-precision values, as ROOT's PtEtaPhiMVector does.
"""

from __future__ import annotations

import multiprocessing as mp
import time
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

import config as cfg

STEP = 1_000_000          # events per chunk read from a ROOT file
Z_MASS = cfg.Z_MASS_NOMINAL


# ---------------------------------------------------------------------------
# Kinematics
# ---------------------------------------------------------------------------


def p4(pt, eta, phi, m):
    """(E, px, py, pz) in float64 from pT, eta, phi, m arrays."""
    pt, eta, phi, m = (np.asarray(x, dtype=np.float64) for x in (pt, eta, phi, m))
    px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    e = np.sqrt(px * px + py * py + pz * pz + m * m)
    return np.stack([e, px, py, pz])


def mass(v):
    """Invariant mass of four-vector(s) v = (E, px, py, pz); negative m^2 -> -sqrt."""
    m2 = v[0] ** 2 - v[1] ** 2 - v[2] ** 2 - v[3] ** 2
    return np.sign(m2) * np.sqrt(np.abs(m2))


def delta_r(eta1, phi1, eta2, phi2):
    """ROOT::VecOps::DeltaR: sqrt(deta^2 + dphi^2), dphi wrapped into [-pi, pi]."""
    dphi = np.mod(phi2 - phi1 + np.pi, 2 * np.pi) - np.pi
    return np.sqrt((eta2 - eta1) ** 2 + dphi**2)


def _sip3d(dxy, dz, dxyerr, dzerr):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.sqrt(dxy * dxy + dz * dz) / np.sqrt(dxyerr * dxyerr + dzerr * dzerr)


# ---------------------------------------------------------------------------
# 1. The dimuon pass
# ---------------------------------------------------------------------------

DIMUON_CACHE = cfg.CACHE_DIR / "dimuon_events.npz"
DIMUON_BRANCHES = ["nMuon", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge"]


def collins_soper_cos(lep, alep):
    """cos(theta*) in the Collins-Soper frame, from the negative lepton's and the
    positive lepton's four-vectors (Collins & Soper 1977), signed by the pair's
    longitudinal direction as is usual for pp collisions."""
    q = lep + alep
    qm = mass(q)
    qt2 = q[1] ** 2 + q[2] ** 2
    lp, lm = (lep[0] + lep[3]) / np.sqrt(2), (lep[0] - lep[3]) / np.sqrt(2)
    ap, am = (alep[0] + alep[3]) / np.sqrt(2), (alep[0] - alep[3]) / np.sqrt(2)
    cos = 2 * (lp * am - lm * ap) / (qm * np.sqrt(qm * qm + qt2))
    return cos * np.sign(q[3])


def dimuon_pass(path: Path | None = None, *, refresh: bool = False) -> dict:
    """Opposite-charge pairs from events with exactly two muons (record 12342)."""
    if DIMUON_CACHE.exists() and not refresh:
        return dict(np.load(DIMUON_CACHE))
    path = path or cfg.RAW_DIR / f"{cfg.DIMUON_FILE}.root"
    tree = uproot.open(path)["Events"]
    out = {k: [] for k in ("mass", "eta1", "eta2", "pt", "y", "cos_cs")}
    counts = {"events": 0, "two_muons": 0, "opposite_charge": 0}
    t0 = time.time()
    for arr in tree.iterate(DIMUON_BRANCHES, step_size=STEP * 2, library="ak"):
        counts["events"] += len(arr)
        two = arr[arr["nMuon"] == 2]
        counts["two_muons"] += len(two)
        q = ak.to_numpy(two["Muon_charge"]).reshape(-1, 2)
        os_ = q[:, 0] != q[:, 1]
        counts["opposite_charge"] += int(os_.sum())
        cols = {b: ak.to_numpy(two[b]).reshape(-1, 2)[os_] for b in DIMUON_BRANCHES[1:]}
        v1 = p4(cols["Muon_pt"][:, 0], cols["Muon_eta"][:, 0], cols["Muon_phi"][:, 0],
                cols["Muon_mass"][:, 0])
        v2 = p4(cols["Muon_pt"][:, 1], cols["Muon_eta"][:, 1], cols["Muon_phi"][:, 1],
                cols["Muon_mass"][:, 1])
        v = v1 + v2
        m = mass(v)
        neg_first = cols["Muon_charge"][:, 0] < 0
        lep = np.where(neg_first, v1, v2)
        alep = np.where(neg_first, v2, v1)
        out["mass"].append(m.astype(np.float32))
        out["eta1"].append(cols["Muon_eta"][:, 0].astype(np.float16))
        out["eta2"].append(cols["Muon_eta"][:, 1].astype(np.float16))
        out["pt"].append(np.hypot(v[1], v[2]).astype(np.float32))
        with np.errstate(divide="ignore", invalid="ignore"):
            out["y"].append((0.5 * np.log((v[0] + v[3]) / (v[0] - v[3]))).astype(np.float16))
            out["cos_cs"].append(collins_soper_cos(lep, alep).astype(np.float16))
        print(f"    dimuon pass: {counts['events'] / 1e6:5.1f} M events "
              f"({time.time() - t0:4.0f} s)", flush=True)
    res = {k: np.concatenate(v) for k, v in out.items()}
    res.update({f"n_{k}": np.int64(v) for k, v in counts.items()})
    np.savez(DIMUON_CACHE, **res)
    return res


DIELECTRON_CACHE = cfg.CACHE_DIR / "dielectron_events.npz"


def dielectron_pass(*, refresh: bool = False) -> dict:
    """Z -> ee calibration sample: DoubleElectron events with exactly two good
    electrons of opposite charge (the electron counterpart of the dimuon pass)."""
    if DIELECTRON_CACHE.exists() and not refresh:
        return dict(np.load(DIELECTRON_CACHE))
    branches = list(cfg.ELECTRON_BRANCHES)
    out = {k: [] for k in ("mass", "eta1", "eta2")}
    for name in ("Run2012B_DoubleElectron", "Run2012C_DoubleElectron"):
        tree = uproot.open(cfg.RAW_DIR / f"{name}.root")["Events"]
        for arr in tree.iterate(branches, step_size=STEP * 2, library="ak"):
            el = good_leptons(arr, "e")
            two = ak.to_numpy(ak.num(el["pt"]) == 2)
            c = {f: ak.to_numpy(el[f][two]).reshape(-1, 2) for f in ("pt", "eta", "phi", "mass",
                                                                     "charge")}
            os_ = (c["charge"][:, 0] != c["charge"][:, 1]) & (np.max(c["pt"], axis=1) > 20) & (
                np.min(c["pt"], axis=1) > 10)
            v = p4(c["pt"][os_, 0], c["eta"][os_, 0], c["phi"][os_, 0], c["mass"][os_, 0]) + \
                p4(c["pt"][os_, 1], c["eta"][os_, 1], c["phi"][os_, 1], c["mass"][os_, 1])
            out["mass"].append(mass(v).astype(np.float32))
            out["eta1"].append(c["eta"][os_, 0].astype(np.float16))
            out["eta2"].append(c["eta"][os_, 1].astype(np.float16))
    res = {k: np.concatenate(v) for k, v in out.items()}
    np.savez(DIELECTRON_CACHE, **res)
    return res


# ---------------------------------------------------------------------------
# 2. The outreach replica (record 12360, skim.cxx)
# ---------------------------------------------------------------------------

LEPTON_BRANCHES = list(cfg.MUON_BRANCHES + cfg.ELECTRON_BRANCHES + cfg.EVENT_ID_BRANCHES)
COMBINATIONS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]   # ROOT's Combinations(v, 2)


def pair_closest_to_z(v, charge, *, emulate_int: bool = True):
    """Indices (i1, i2) of the opposite-charge pair closest to 91.2 GeV.

    ``v``: (4, N, 4) four-vectors of the four leptons; ``charge``: (N, 4).
    ``skim.cxx`` keeps the best mass so far in ``auto best_mass = -1;`` -- an
    int -- so each accepted pair's mass is truncated to whole GeV before it is
    compared with the next pair. With ``emulate_int`` we reproduce that
    exactly; without it, the comparison uses the true mass.
    """
    n = charge.shape[0]
    best = np.full(n, -1.0)
    i1 = np.zeros(n, dtype=int)
    i2 = np.zeros(n, dtype=int)
    for a, b in COMBINATIONS:
        m = mass(v[:, :, a] + v[:, :, b])
        upd = (charge[:, a] != charge[:, b]) & (np.abs(Z_MASS - m) < np.abs(Z_MASS - best))
        best = np.where(upd, np.trunc(m) if emulate_int else m, best)
        i1 = np.where(upd, a, i1)
        i2 = np.where(upd, b, i2)
    return i1, i2


def _take(arr2d, idx):
    return np.take_along_axis(arr2d, idx[:, None], axis=1)[:, 0]


def _z_pair_masses(za, zb):
    """Order two Z candidates as skim.cxx does: the one closer to 91.2 GeV first."""
    ma, mb = mass(za), mass(zb)
    first = np.abs(ma - Z_MASS) < np.abs(mb - Z_MASS)
    z1 = np.where(first, za, zb)
    z2 = np.where(first, zb, za)
    return z1, z2


def _regular(arr, prefix, n, fields):
    """Jagged lepton branches of events with exactly/at least n entries -> (N, n) float64."""
    return {f: ak.to_numpy(ak.pad_none(arr[f"{prefix}_{f}"], n, clip=True)).astype(np.float64)
            for f in fields}


def replica_samekind(arr, flavour: str, *, emulate_int: bool = True) -> dict:
    """skim.cxx MinimalSelection + ReconstructHiggs for 4mu or 4e."""
    p = "Muon" if flavour == "mu" else "Electron"
    iso = "pfRelIso04_all" if flavour == "mu" else "pfRelIso03_all"
    ptmin, etamax = (5.0, 2.4) if flavour == "mu" else (7.0, 2.5)
    nlep = arr[f"n{p}"]
    sip = _sip3d(arr[f"{p}_dxy"], arr[f"{p}_dz"], arr[f"{p}_dxyErr"], arr[f"{p}_dzErr"])
    sel = ((nlep >= 4)
           & ak.all(np.abs(arr[f"{p}_{iso}"]) < 0.40, axis=1)
           & ak.all(arr[f"{p}_pt"] > ptmin, axis=1) & ak.all(np.abs(arr[f"{p}_eta"]) < etamax, axis=1)
           & ak.all(sip < 4, axis=1) & ak.all(np.abs(arr[f"{p}_dxy"]) < 0.5, axis=1)
           & ak.all(np.abs(arr[f"{p}_dz"]) < 1.0, axis=1)
           & (nlep == 4) & (ak.sum(arr[f"{p}_charge"] == 1, axis=1) == 2)
           & (ak.sum(arr[f"{p}_charge"] == -1, axis=1) == 2))
    ev = arr[ak.to_numpy(sel)]
    if len(ev) == 0:
        return _empty_result()
    c = _regular(ev, p, 4, ["pt", "eta", "phi", "mass", "charge"])
    v = p4(c["pt"], c["eta"], c["phi"], c["mass"])            # (4, N, 4)
    i1, i2 = pair_closest_to_z(v, c["charge"], emulate_int=emulate_int)
    # second Z: the leptons not in the first pair, in index order. If no pair
    # qualified (every opposite-charge pair above 183.4 GeV), skim.cxx keeps
    # i1 = i2 = 0, and its Delta-R filter then rejects the event (dR = 0).
    used = (np.arange(4) == i1[:, None]) | (np.arange(4) == i2[:, None])
    rest = np.argsort(used, axis=1, kind="stable")[:, :2]
    j1, j2 = rest[:, 0], rest[:, 1]
    ok = i1 != i2
    for a, b in ((i1, i2), (j1, j2)):
        dr = delta_r(_take(c["eta"], a), _take(c["phi"], a), _take(c["eta"], b), _take(c["phi"], b))
        ok &= ~(dr < 0.02)
    rows = np.arange(len(ev))
    za = v[:, rows, i1] + v[:, rows, i2]
    zb = v[:, rows, j1] + v[:, rows, j2]
    z1, z2 = _z_pair_masses(za, zb)
    return _finish(ev, ok, z1, z2)


def replica_2e2mu(arr) -> dict:
    """skim.cxx MinimalSelection + ReconstructHiggs for 2e2mu."""
    pre = ak.to_numpy((arr["nElectron"] >= 2) & (arr["nMuon"] >= 2))
    ev = arr[pre]
    if len(ev) == 0:
        return _empty_result()
    ok = ak.to_numpy(ak.all(np.abs(ev["Electron_eta"]) < 2.5, axis=1)
                     & ak.all(np.abs(ev["Muon_eta"]) < 2.4, axis=1))
    mu_sorted = ak.to_numpy(ak.pad_none(ak.sort(ev["Muon_pt"], ascending=False), 2, clip=True))
    el_sorted = ak.to_numpy(ak.pad_none(ak.sort(ev["Electron_pt"], ascending=False), 2, clip=True))
    ok &= ((mu_sorted[:, 0] > 20) & (mu_sorted[:, 1] > 10)) | (
        (el_sorted[:, 0] > 20) & (el_sorted[:, 1] > 10))
    mu = _regular(ev, "Muon", 2, ["pt", "eta", "phi", "mass"])
    el = _regular(ev, "Electron", 2, ["pt", "eta", "phi", "mass"])
    mu_dr = delta_r(mu["eta"][:, 0], mu["phi"][:, 0], mu["eta"][:, 1], mu["phi"][:, 1])
    el_dr = delta_r(el["eta"][:, 0], el["phi"][:, 0], el["eta"][:, 1], el["phi"][:, 1])
    ok &= ~((mu_dr < 0.02) | (el_dr < 0.02))
    ok &= ak.to_numpy(ak.all(np.abs(ev["Electron_pfRelIso03_all"]) < 0.40, axis=1)
                      & ak.all(np.abs(ev["Muon_pfRelIso04_all"]) < 0.40, axis=1))
    for p in ("Electron", "Muon"):
        sip = _sip3d(ev[f"{p}_dxy"], ev[f"{p}_dz"], ev[f"{p}_dxyErr"], ev[f"{p}_dzErr"])
        ok &= ak.to_numpy(ak.all(sip < 4, axis=1) & ak.all(np.abs(ev[f"{p}_dxy"]) < 0.5, axis=1)
                          & ak.all(np.abs(ev[f"{p}_dz"]) < 1.0, axis=1))
    ok &= ak.to_numpy((ak.sum(ev["Electron_charge"], axis=1) == 0)
                      & (ak.sum(ev["Muon_charge"], axis=1) == 0))
    za = p4(mu["pt"][:, 0], mu["eta"][:, 0], mu["phi"][:, 0], mu["mass"][:, 0]) + \
        p4(mu["pt"][:, 1], mu["eta"][:, 1], mu["phi"][:, 1], mu["mass"][:, 1])
    zb = p4(el["pt"][:, 0], el["eta"][:, 0], el["phi"][:, 0], el["mass"][:, 0]) + \
        p4(el["pt"][:, 1], el["eta"][:, 1], el["phi"][:, 1], el["mass"][:, 1])
    z1, z2 = _z_pair_masses(za, zb)
    return _finish(ev, ok, z1, z2)


def _empty_result() -> dict:
    return {k: np.zeros(0, dtype=t) for k, t in
            (("run", np.int64), ("lumi", np.int64), ("event", np.uint64),
             ("m4l", np.float64), ("mz1", np.float64), ("mz2", np.float64))}


def _finish(ev, ok, z1, z2) -> dict:
    m1, m2 = mass(z1), mass(z2)
    ok = ok & (m1 > 40) & (m1 < 120) & (m2 > 12) & (m2 < 120)
    h = z1 + z2
    return {"run": ak.to_numpy(ev["run"])[ok].astype(np.int64),
            "lumi": ak.to_numpy(ev["luminosityBlock"])[ok].astype(np.int64),
            "event": ak.to_numpy(ev["event"])[ok].astype(np.uint64),
            "m4l": mass(h)[ok], "mz1": m1[ok], "mz2": m2[ok]}


# ---------------------------------------------------------------------------
# 3. The good-lepton skim (input to our own selection)
# ---------------------------------------------------------------------------

#: Lepton identification, as in the outreach analysis and the CMS 2012 paper:
#: kinematic acceptance, relative isolation < 0.4, and a track compatible
#: with the primary vertex (3D impact-parameter significance < 4).
GOOD = {"mu": {"pt": 5.0, "eta": 2.4, "iso": "pfRelIso04_all"},
        "e": {"pt": 7.0, "eta": 2.5, "iso": "pfRelIso03_all"}}
SKIM_FIELDS = ("pt", "eta", "phi", "mass", "charge", "iso", "sip3d", "dxy", "dz")


def good_leptons(arr, flavour: str):
    p = "Muon" if flavour == "mu" else "Electron"
    g = GOOD[flavour]
    sip = _sip3d(arr[f"{p}_dxy"], arr[f"{p}_dz"], arr[f"{p}_dxyErr"], arr[f"{p}_dzErr"])
    good = ((arr[f"{p}_pt"] > g["pt"]) & (np.abs(arr[f"{p}_eta"]) < g["eta"])
            & (np.abs(arr[f"{p}_{g['iso']}"]) < 0.40) & (sip < 4)
            & (np.abs(arr[f"{p}_dxy"]) < 0.5) & (np.abs(arr[f"{p}_dz"]) < 1.0))
    rec = {"pt": arr[f"{p}_pt"], "eta": arr[f"{p}_eta"], "phi": arr[f"{p}_phi"],
           "mass": arr[f"{p}_mass"], "charge": arr[f"{p}_charge"],
           "iso": arr[f"{p}_{g['iso']}"], "sip3d": sip, "dxy": arr[f"{p}_dxy"],
           "dz": arr[f"{p}_dz"]}
    return {k: v[good] for k, v in rec.items()}


def skim_chunk(arr) -> dict:
    mu, el = good_leptons(arr, "mu"), good_leptons(arr, "e")
    nmu, nel = ak.num(mu["pt"]), ak.num(el["pt"])
    keep = ak.to_numpy((nmu + nel) >= 4)
    out = {"run": ak.to_numpy(arr["run"])[keep].astype(np.int64),
           "lumi": ak.to_numpy(arr["luminosityBlock"])[keep].astype(np.int64),
           "event": ak.to_numpy(arr["event"])[keep].astype(np.uint64),
           "n_mu": ak.to_numpy(nmu)[keep].astype(np.int16),
           "n_e": ak.to_numpy(nel)[keep].astype(np.int16)}
    for flav, rec in (("mu", mu), ("e", el)):
        for f in SKIM_FIELDS:
            out[f"{flav}_{f}"] = ak.to_numpy(ak.flatten(rec[f][keep])).astype(
                np.int8 if f == "charge" else np.float32)
    return out


# ---------------------------------------------------------------------------
# Driver over the four-lepton files
# ---------------------------------------------------------------------------


def replica_cache(name: str) -> Path:
    return cfg.CACHE_DIR / f"replica_{name}.npz"


def skim_cache(name: str) -> Path:
    return cfg.CACHE_DIR / f"skim_{name}.npz"


def process_file(name: str, refresh: bool = False) -> str:
    """Outreach replica (both int emulations) + good-lepton skim of one file."""
    if replica_cache(name).exists() and skim_cache(name).exists() and not refresh:
        return name
    tree = uproot.open(cfg.RAW_DIR / f"{name}.root")["Events"]
    states = cfg.FINAL_STATES[name]
    rep = {f"{s}{tag}": [] for s in states for tag in ("", "_fixed")}
    skim = []
    n_events = 0
    t0 = time.time()
    for arr in tree.iterate(LEPTON_BRANCHES, step_size=STEP, library="ak"):
        n_events += len(arr)
        for s in states:
            if s == "2e2mu":
                r = replica_2e2mu(arr)
                rep[s].append(r)
                rep[f"{s}_fixed"].append(r)        # no pairing ambiguity in 2e2mu
            else:
                flav = "mu" if s == "4mu" else "e"
                rep[s].append(replica_samekind(arr, flav, emulate_int=True))
                rep[f"{s}_fixed"].append(replica_samekind(arr, flav, emulate_int=False))
        skim.append(skim_chunk(arr))
        print(f"    {name}: {n_events / 1e6:5.1f} M events ({time.time() - t0:4.0f} s)",
              flush=True)
    out = {"n_events": np.int64(n_events)}
    for key, parts in rep.items():
        for f in parts[0]:
            out[f"{key}__{f}"] = np.concatenate([p[f] for p in parts])
    np.savez(replica_cache(name), **out)
    sk = {f: np.concatenate([p[f] for p in skim]) for f in skim[0]}
    sk["n_events"] = np.int64(n_events)
    np.savez(skim_cache(name), **sk)
    return name


def _worker(args):
    import threadpoolctl
    with threadpoolctl.threadpool_limits(1):
        return process_file(*args)


def process_all(names=None, *, workers: int = 3, refresh: bool = False) -> list[str]:
    names = list(names or (*cfg.SIMULATION, *cfg.DATA_4L))
    todo = [n for n in names if refresh or not (replica_cache(n).exists()
                                                and skim_cache(n).exists())]
    todo = [n for n in todo if (cfg.RAW_DIR / f"{n}.root").exists()]
    if not todo:
        return names
    if workers <= 1 or len(todo) == 1:
        for n in todo:
            process_file(n, refresh)
    else:
        ctx = mp.get_context("fork")
        with ctx.Pool(min(workers, len(todo))) as pool:
            for done in pool.imap_unordered(_worker, [(n, refresh) for n in todo]):
                print(f"  processed {done}", flush=True)
    return names


def load_replica(name: str) -> dict:
    """{state: {field: array}} for one file's replica (keys like '4mu', '4mu_fixed')."""
    z = np.load(replica_cache(name))
    out = {}
    for k in z.files:
        if "__" in k:
            state, f = k.split("__")
            out.setdefault(state, {})[f] = z[k]
    out["n_events"] = int(z["n_events"])
    return out


def load_skim(name: str) -> dict:
    return dict(np.load(skim_cache(name)))


if __name__ == "__main__":
    import sys
    if "dimuon" in sys.argv[1:]:
        d = dimuon_pass(refresh=True)
        print({k: int(v) for k, v in d.items() if k.startswith("n_")})
    else:
        process_all(sys.argv[1:] or None)
