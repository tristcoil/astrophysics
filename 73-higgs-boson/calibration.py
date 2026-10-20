"""Lepton momentum scale: is a 45 GeV muon in the data measured as 45 GeV?

The outreach files carry the muon momenta and electron energies as
reconstructed in 2012, without the calibrations CMS applied in its own
analyses. Their scale shows up directly in the position of the Z peak. We fit
the Z -> mu mu and Z -> e e peaks (Breit-Wigner x Crystal Ball, as in
``dimuon.fit_z``) in the data -- millions of pairs -- and in the simulation,
where the same Z bosons appear inside the ZZ samples, separately for leptons
in the central (barrel) and forward (endcap) parts of the detector. The ratio
simulation/data is the factor by which the data momenta must be multiplied so
that both are on the same scale: the Higgs signal model comes from the
simulation, so this, not the PDG Z mass, is what the Higgs mass needs.
"""

from __future__ import annotations

import numpy as np

import dimuon
import events

#: |eta| regions: muon barrel / endcap (CMS muon system), ECAL barrel / endcap.
REGIONS = {"mu": {"barrel": (0.0, 0.9), "endcap": (1.2, 2.4)},
           "e": {"barrel": (0.0, 1.479), "endcap": (1.566, 2.5)}}
WINDOW = (80.0, 102.0)


def _in(eta, region):
    a = np.abs(eta)
    return (a >= region[0]) & (a < region[1])


def mc_pairs(flavour: str) -> dict:
    """All opposite-charge same-flavour pairs of good leptons in the ZZ
    simulation (true Z pairs and a smooth combinatorial background)."""
    name = "ZZTo4mu" if flavour == "mu" else "ZZTo4e"
    s = events.load_skim(name)
    n = s[f"n_{flavour}"]
    off = np.concatenate([[0], np.cumsum(n)])
    pt, eta, phi, m, q = (s[f"{flavour}_{f}"] for f in ("pt", "eta", "phi", "mass", "charge"))
    masses, e1, e2 = [], [], []
    for i in range(len(n)):
        a, b = off[i], off[i + 1]
        for j in range(a, b):
            for k in range(j + 1, b):
                if q[j] != q[k]:
                    v = events.p4(pt[[j, k]], eta[[j, k]], phi[[j, k]], m[[j, k]])
                    masses.append(events.mass(v[:, 0] + v[:, 1]))
                    e1.append(eta[j])
                    e2.append(eta[k])
    return {"mass": np.array(masses), "eta1": np.array(e1), "eta2": np.array(e2)}


def data_pairs(flavour: str) -> dict:
    return events.dimuon_pass() if flavour == "mu" else events.dielectron_pass()


def fit_peaks(pairs: dict, flavour: str) -> dict:
    out = {}
    e1 = pairs["eta1"].astype(np.float64)
    e2 = pairs["eta2"].astype(np.float64)
    m = pairs["mass"].astype(np.float64)
    for reg, rng in REGIONS[flavour].items():
        sel = _in(e1, rng) & _in(e2, rng) & (m > WINDOW[0]) & (m < WINDOW[1])
        r = dimuon.fit_z(m[sel], WINDOW)
        out[reg] = {"mass": r["mass"], "mass_err": r["mass_err"], "resolution": r["resolution"],
                    "n": int(sel.sum()), "chi2": r["chi2"], "ndof": r["ndof"],
                    "converged": r["converged"]}
    return out


def calibrate() -> dict:
    """Peak positions in data and simulation; scale = M_sim / M_data per region."""
    res = {}
    for flav in ("mu", "e"):
        d = fit_peaks(data_pairs(flav), flav)
        s = fit_peaks(mc_pairs(flav), flav)
        res[flav] = {}
        for reg in REGIONS[flav]:
            k = s[reg]["mass"] / d[reg]["mass"]
            kerr = k * np.hypot(s[reg]["mass_err"] / s[reg]["mass"],
                                d[reg]["mass_err"] / d[reg]["mass"])
            res[flav][reg] = {"data": d[reg], "sim": s[reg], "scale": k, "scale_err": kerr}
    return res


def scale_function(cal: dict):
    """Per-lepton momentum factor for data: barrel/endcap values, linear in |eta|
    across the gap between them."""
    def f(flav, eta):
        reg = REGIONS[flav]
        kb, ke = cal[flav]["barrel"]["scale"], cal[flav]["endcap"]["scale"]
        a = abs(eta)
        b_hi, e_lo = reg["barrel"][1], reg["endcap"][0]
        if a <= b_hi:
            return kb
        if a >= e_lo:
            return ke
        t = (a - b_hi) / (e_lo - b_hi)
        return (1 - t) * kb + t * ke
    return f
