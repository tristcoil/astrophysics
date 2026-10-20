"""The particle zoo in one plot: the dimuon mass spectrum, and the masses we read from it.

Any particle that can decay into a muon and an antimuon leaves a bump in the
distribution of the pair's invariant mass at the particle's own mass. From
61.5 million 2012 events with two muons we build the spectrum CMS published as
an outreach example (record 12342), then fit each bump to measure the masses
of nine particles -- from the eta meson (0.548 GeV) to the Z boson (91.19 GeV)
-- and compare them with the Particle Data Group's 2026 values.

How well a peak is measured depends on the detector: muons crossing the
central "barrel" (|eta| < 0.9) are measured more precisely than those in the
endcaps, which is what separates the three Upsilon states cleanly.
"""

from __future__ import annotations

import numpy as np
from iminuit import Minuit
from particle import Particle
from scipy.special import erf, wofz

import config as cfg
import events

#: Peaks we fit: label, PDG id(s), fit window (GeV), signal model.
RESONANCES = [
    {"key": "eta", "label": "η", "pdgid": 221, "window": (0.47, 0.63), "model": "gauss"},
    {"key": "omega", "label": "ρ/ω", "pdgid": 223, "window": (0.66, 0.90), "model": "rho_omega"},
    {"key": "phi", "label": "φ", "pdgid": 333, "window": (0.90, 1.14), "model": "voigt"},
    {"key": "jpsi", "label": "J/ψ", "pdgid": 443, "window": (2.80, 3.40), "model": "cb"},
    {"key": "psi2s", "label": "ψ(2S)", "pdgid": 100443, "window": (3.40, 3.98), "model": "cb"},
    {"key": "upsilon", "label": "ϒ(1S,2S,3S)", "pdgid": (553, 100553, 200553),
     "window": (8.6, 11.2), "model": "upsilon"},
    {"key": "z", "label": "Z", "pdgid": 23, "window": (70.0, 110.0), "model": "z"},
]

BARREL = 0.9           # |eta| below which a muon is in the barrel


def pdg_mass(pdgid: int) -> tuple[float, float, float]:
    """(mass, uncertainty, width) in GeV from the PDG 2026 table in `particle`."""
    p = Particle.from_pdgid(pdgid)
    width = p.width if p.width is not None else 0.0
    return p.mass / 1000, max(p.mass_upper or 0, p.mass_lower or 0) / 1000, width / 1000


# ---------------------------------------------------------------------------
# Spectrum
# ---------------------------------------------------------------------------


def load() -> dict:
    return events.dimuon_pass()


def categories(d: dict) -> dict[str, np.ndarray]:
    """Event classes by where the two muons crossed the detector."""
    a1, a2 = np.abs(d["eta1"].astype(np.float32)), np.abs(d["eta2"].astype(np.float32))
    return {"all": np.ones(len(a1), dtype=bool),
            "barrel": (a1 < BARREL) & (a2 < BARREL),
            "endcap": (a1 > 1.2) & (a2 > 1.2)}


def spectrum(d: dict, lo: float = 0.25, hi: float = 300.0, nbins: int = 1500) -> dict:
    """Log-binned spectrum: counts and events per GeV."""
    edges = np.geomspace(lo, hi, nbins + 1)
    n = np.histogram(d["mass"], edges)[0]
    return {"edges": edges, "counts": n, "per_gev": n / np.diff(edges)}


def outreach_histogram(d: dict) -> np.ndarray:
    """The CMS example's own histogram: 30,000 linear bins from 0.25 to 300 GeV."""
    return np.histogram(d["mass"], np.linspace(0.25, 300.0, 30001))[0]


# ---------------------------------------------------------------------------
# Line shapes (all unnormalised; normalised numerically over the fit window)
# ---------------------------------------------------------------------------


def gauss(x, mu, s):
    return np.exp(-0.5 * ((x - mu) / s) ** 2)


def crystal_ball(x, mu, s, alpha, n):
    """Gaussian core with a power-law tail on the low side (photon radiation)."""
    t = (x - mu) / s
    a = np.abs(alpha)
    A = (n / a) ** n * np.exp(-0.5 * a * a)
    B = n / a - a
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        tail = A * np.power(np.maximum(B - t, 1e-12), -n)
    return np.where(t > -a, np.exp(-0.5 * t * t), tail)


def voigt(x, mu, s, gamma):
    """Breit-Wigner of full width gamma convolved with a Gaussian of width s."""
    z = ((x - mu) + 0.5j * gamma) / (s * np.sqrt(2))
    return np.real(wofz(z))


def breit_wigner_rel(x, m, gamma):
    """Relativistic Breit-Wigner in m (as used for the broad rho)."""
    return x * x * m * m * gamma * gamma / ((x * x - m * m) ** 2 + m * m * gamma * gamma)


def cheb(x, lo, hi, coeffs):
    """exp(polynomial) background in the scaled variable u in [-1, 1]."""
    u = 2 * (x - lo) / (hi - lo) - 1
    return np.exp(np.polynomial.chebyshev.chebval(u, np.concatenate([[0.0], coeffs])))


# ---------------------------------------------------------------------------
# Binned Poisson fits
# ---------------------------------------------------------------------------


class PeakFit:
    """Extended binned Poisson likelihood for one window: sum of shapes x yields."""

    def __init__(self, mass: np.ndarray, window: tuple, bin_width: float):
        lo, hi = window
        self.lo, self.hi = lo, hi
        nb = int(round((hi - lo) / bin_width))
        self.edges = np.linspace(lo, hi, nb + 1)
        self.n = np.histogram(mass, self.edges)[0].astype(float)
        e = self.edges
        # Simpson nodes per bin
        self.x3 = np.stack([e[:-1], 0.5 * (e[:-1] + e[1:]), e[1:]])
        self.w = np.diff(e)
        self.grid = np.linspace(lo, hi, 4001)

    def integrate_bins(self, f):
        y = f(self.x3)
        return (y[0] + 4 * y[1] + y[2]) / 6 * self.w

    def norm(self, f):
        return np.trapezoid(f(self.grid), self.grid)

    def nll(self, mu):
        mu = np.maximum(mu, 1e-300)
        return 2 * np.sum(mu - self.n * np.log(mu))


def _expected(fit: PeakFit, comps):
    """comps: list of (yield, shape function)."""
    total = np.zeros_like(fit.n)
    for y, f in comps:
        total += y * fit.integrate_bins(f) / fit.norm(f)
    return total


def fit_resonance(mass: np.ndarray, res: dict) -> dict:
    """Fit one peak (or a family of peaks) in its window."""
    lo, hi = res["window"]
    model = res["model"]
    if model == "z":
        fit = PeakFit(mass, (lo, hi), 0.1)
    elif model == "upsilon":
        fit = PeakFit(mass, (lo, hi), 0.01)
    elif model in ("cb",):
        fit = PeakFit(mass, (lo, hi), 0.004)
    else:
        fit = PeakFit(mass, (lo, hi), 0.002)
    nsum = fit.n.sum()

    if model == "upsilon":
        m1, m2, m3 = (pdg_mass(p)[0] for p in res["pdgid"])

        def comps(p):
            (n1, n2, n3, nb, mu1, mu2, mu3, s1, f2, r2, alpha, nn, c1, c2, c3) = p
            out = [(nb, lambda x: cheb(x, lo, hi, [c1, c2, c3]))]
            for n_, mu_ in ((n1, mu1), (n2, mu2), (n3, mu3)):
                s = s1 * mu_ / mu1          # resolution scales with mass
                out.append((n_, lambda x, mu_=mu_, s=s: (1 - f2) * crystal_ball(x, mu_, s, alpha, nn)
                            + f2 * gauss(x, mu_, r2 * s)))
            return out

        names = ["n1", "n2", "n3", "nb", "mu1", "mu2", "mu3", "s1", "f2", "r2", "alpha", "n",
                 "c1", "c2", "c3"]
        start = [0.1 * nsum, 0.04 * nsum, 0.03 * nsum, 0.8 * nsum, m1, m2, m3, 0.08, 0.3, 2.0,
                 1.5, 3.0, -0.5, 0.0, 0.0]
        limits = {"s1": (0.02, 0.3), "f2": (0.0, 1.0), "r2": (1.0, 5.0), "alpha": (0.3, 5),
                  "n": (1.01, 50)}
        peak_params = [("upsilon1s", "mu1", 553), ("upsilon2s", "mu2", 100553),
                       ("upsilon3s", "mu3", 200553)]
    elif model == "z":
        mz, _, gz = pdg_mass(23)

        def comps(p):
            ns, nb, mu, s, gamma, alpha, nn, frac, c1, c2 = p
            sig = lambda x: (1 - frac) * voigt(x, mu, s, gamma) + frac * crystal_ball(
                x, mu, s, alpha, nn) * voigt(mu, mu, s, gamma)
            return [(ns, sig), (nb, lambda x: cheb(x, lo, hi, [c1, c2]))]

        names = ["ns", "nb", "mu", "s", "gamma", "alpha", "n", "frac", "c1", "c2"]
        start = [0.9 * nsum, 0.1 * nsum, mz, 1.2, gz, 1.0, 3.0, 0.2, -1.0, 0.0]
        limits = {"s": (0.3, 4.0), "gamma": (0.5, 6.0), "alpha": (0.2, 5), "n": (1.01, 50),
                  "frac": (0.0, 1.0)}
        peak_params = [("z", "mu", 23)]
    elif model == "cb":
        m0, _, _ = pdg_mass(res["pdgid"])

        def comps(p):
            ns, nb, mu, s, f2, r2, alpha, nn, c1, c2 = p
            sig = lambda x: (1 - f2) * crystal_ball(x, mu, s, alpha, nn) + f2 * gauss(x, mu, r2 * s)
            return [(ns, sig), (nb, lambda x: cheb(x, lo, hi, [c1, c2]))]

        names = ["ns", "nb", "mu", "s", "f2", "r2", "alpha", "n", "c1", "c2"]
        start = [0.3 * nsum, 0.7 * nsum, m0, 0.025 * m0 / 3.1, 0.3, 2.0, 1.5, 3.0, -0.5, 0.0]
        limits = {"s": (0.003, 0.2), "f2": (0.0, 1.0), "r2": (1.0, 5.0), "alpha": (0.3, 5),
                  "n": (1.01, 50)}
        peak_params = [(res["key"], "mu", res["pdgid"])]
    elif model in ("gauss", "voigt"):
        m0, _, w0 = pdg_mass(res["pdgid"])
        width = w0 if model == "voigt" else 0.0

        def comps(p):
            ns, nb, mu, s, c1, c2, c3 = p
            sig = (lambda x: voigt(x, mu, s, width)) if width > 0 else (lambda x: gauss(x, mu, s))
            return [(ns, sig), (nb, lambda x: cheb(x, lo, hi, [c1, c2, c3]))]

        names = ["ns", "nb", "mu", "s", "c1", "c2", "c3"]
        start = [0.05 * nsum, 0.95 * nsum, m0, 0.008, -0.5, 0.0, 0.0]
        limits = {"s": (0.001, 0.05)}
        peak_params = [(res["key"], "mu", res["pdgid"])]
    elif model == "rho_omega":
        mo, _, go = pdg_mass(223)
        mr, _, gr = pdg_mass(113)

        def comps(p):
            no, nr, nb, mu, s, c1, c2, c3 = p
            return [(no, lambda x: voigt(x, mu, s, go)),
                    (nr, lambda x: breit_wigner_rel(x, mr, gr)),
                    (nb, lambda x: cheb(x, lo, hi, [c1, c2, c3]))]

        names = ["no", "nr", "nb", "mu", "s", "c1", "c2", "c3"]
        start = [0.05 * nsum, 0.05 * nsum, 0.9 * nsum, mo, 0.008, -0.5, 0.0, 0.0]
        limits = {"s": (0.001, 0.05), "nr": (0, None)}
        peak_params = [("omega", "mu", 223)]
    else:
        raise ValueError(model)

    def cost(*p):
        return fit.nll(_expected(fit, comps(p)))

    m = Minuit(cost, *start, name=names)
    m.errordef = Minuit.LIKELIHOOD / 0.5       # cost is 2 * NLL -> Delta = 1
    for k, v in limits.items():
        m.limits[k] = v
    for k in names:
        if k.startswith("n") and k not in ("n",):
            m.limits[k] = (0, None)
    m.strategy = 1
    m.migrad(ncall=20000)
    m.hesse()
    best = np.array(m.values)
    expected = _expected(fit, comps(best))
    chi2 = float(np.sum((fit.n - expected) ** 2 / np.maximum(expected, 1)))
    out = {"window": [lo, hi], "converged": bool(m.valid), "chi2": chi2,
           "ndof": int(len(fit.n) - len(names)), "params": dict(zip(names, best.tolist())),
           "errors": dict(zip(names, np.array(m.errors).tolist())),
           "edges": fit.edges.tolist(), "counts": fit.n.tolist(),
           "expected": expected.tolist(), "peaks": {}}
    for key, par, pid in peak_params:
        m0, m0_err, w0 = pdg_mass(pid)
        mu, err = m.values[par], m.errors[par]
        out["peaks"][key] = {"mass": mu, "mass_err": err, "pdg_mass": m0, "pdg_err": m0_err,
                             "pdg_width": w0, "ratio_minus_1": mu / m0 - 1,
                             "ratio_err": err / m0}
    if model == "z":
        out["peaks"]["z"].update({"width": m.values["gamma"], "width_err": m.errors["gamma"],
                                  "resolution": m.values["s"]})
    return out


def fit_z(mass: np.ndarray, window=(80.0, 102.0), *, gamma_free: bool = False,
          bin_width: float = 0.1) -> dict:
    """Z line shape: relativistic Breit-Wigner convolved with a Crystal Ball
    resolution function (Gaussian core, low-side tail from photon radiation),
    over an exponential continuum (Drell-Yan photon exchange + interference).

    With the width fixed to PDG's 2.4955 GeV this measures the peak mass; with
    it free, the fit trades width against resolution tails, so a free width
    is only indicative.
    """
    lo, hi = window
    fit = PeakFit(mass, window, bin_width)
    mz0, _, gz0 = pdg_mass(23)
    step = 0.01
    grid = np.arange(lo - 15, hi + 15 + step / 2, step)
    kx = np.arange(-15, 15 + step / 2, step)
    nsum = fit.n.sum()

    def signal_on_grid(mu, gamma, s, alpha, nn):
        bw = grid**2 / ((grid**2 - mu**2) ** 2 + mu**2 * gamma**2)
        ker = crystal_ball(kx, 0.0, s, alpha, nn)
        ker /= ker.sum()
        # (BW * r)(m) = sum_k BW(x_k) r(m - x_k): the kernel enters unreversed
        return np.convolve(bw, ker, mode="same")

    def comps(p):
        ns, nb, mu, s, alpha, nn, c1, gamma = p
        sg = signal_on_grid(mu, gamma, s, alpha, nn)
        return [(ns, lambda x: np.interp(x, grid, sg)),
                (nb, lambda x: cheb(x, lo, hi, [c1]))]

    names = ["ns", "nb", "mu", "s", "alpha", "n", "c1", "gamma"]
    start = [0.9 * nsum, 0.1 * nsum, mz0, 1.2, 1.2, 3.0, -1.0, gz0]

    def cost(*p):
        return fit.nll(_expected(fit, comps(p)))

    m = Minuit(cost, *start, name=names)
    m.errordef = 1.0
    m.limits["s"] = (0.3, 5.0)
    m.limits["alpha"] = (0.2, 5.0)
    m.limits["n"] = (1.01, 100.0)
    m.limits["ns"] = (0, None)
    m.limits["nb"] = (0, None)
    m.limits["gamma"] = (0.5, 6.0)
    if not gamma_free:
        m.fixed["gamma"] = True
    m.migrad(ncall=20000)
    m.hesse()
    best = np.array(m.values)
    expected = _expected(fit, comps(best))
    chi2 = float(np.sum((fit.n - expected) ** 2 / np.maximum(expected, 1)))
    return {"window": list(window), "converged": bool(m.valid), "chi2": chi2,
            "ndof": int(len(fit.n) - (len(names) - (0 if gamma_free else 1))),
            "params": dict(zip(names, best.tolist())),
            "errors": dict(zip(names, np.array(m.errors).tolist())),
            "edges": fit.edges.tolist(), "counts": fit.n.tolist(), "expected": expected.tolist(),
            "mass": m.values["mu"], "mass_err": m.errors["mu"], "pdg_mass": mz0,
            "ratio_minus_1": m.values["mu"] / mz0 - 1, "ratio_err": m.errors["mu"] / mz0,
            "width": m.values["gamma"], "width_err": m.errors["gamma"] if gamma_free else 0.0,
            "resolution": m.values["s"]}


def fit_all(d: dict, category: str = "all") -> dict:
    mask = categories(d)[category]
    mass = d["mass"][mask].astype(np.float64)
    results = {}
    for res in RESONANCES:
        lo, hi = res["window"]
        sel = mass[(mass > lo) & (mass < hi)]
        if res["model"] == "z":
            r = fit_z(sel)
            r["peaks"] = {"z": {k: r[k] for k in ("mass", "mass_err", "pdg_mass", "ratio_minus_1",
                                                  "ratio_err", "width", "resolution")}}
            r["peaks"]["z"]["pdg_err"] = pdg_mass(23)[1]
            results["z"] = r
        else:
            results[res["key"]] = fit_resonance(sel, res)
        results[res["key"]]["label"] = res["label"]
    return results


# ---------------------------------------------------------------------------
# Forward-backward asymmetry: the photon and the Z interfere
# ---------------------------------------------------------------------------


def afb(d: dict, edges=None, ymin: float = 0.0) -> dict:
    """A_FB(m) = (N_F - N_B)/(N_F + N_B) in the Collins-Soper frame.

    Because the beams are both protons, the quark's direction is guessed from
    the pair's boost (the quark usually carries more momentum than the
    antiquark); the guess is more often right at large rapidity, where the
    measured asymmetry is less diluted.
    """
    edges = np.geomspace(40, 300, 41) if edges is None else edges
    m = d["mass"]
    c = d["cos_cs"].astype(np.float32)
    y = np.abs(d["y"].astype(np.float32))
    sel = (y >= ymin) & np.isfinite(c)
    nf = np.histogram(m[sel & (c > 0)], edges)[0].astype(float)
    nb = np.histogram(m[sel & (c < 0)], edges)[0].astype(float)
    tot = nf + nb
    a = np.where(tot > 0, (nf - nb) / np.maximum(tot, 1), np.nan)
    err = np.where(tot > 0, np.sqrt(np.maximum(1 - a * a, 0) / np.maximum(tot, 1)), np.nan)
    return {"edges": edges, "afb": a, "err": err, "n": tot}


def afb_quark_level(m, *, s2w: float = 0.23155, quark_weights=None, MZ: float = 91.1876,
                    GZ: float = 2.4955, GF: float = 1.1663785e-5, alpha_mz: float = 1 / 127.95):
    """Leading-order A_FB of q qbar -> gamma*/Z -> mu mu at the quark level
    (Ellis, Stirling & Webber, QCD and Collider Physics, sec. 9.3), mixed over
    u and d quarks. In proton collisions the measured asymmetry is this times
    a dilution factor (the quark direction is guessed from the pair's boost)."""
    w = quark_weights or {"u": 2.0, "d": 1.0}
    s = np.asarray(m, dtype=float) ** 2
    kappa = np.sqrt(2) * GF * MZ**2 / (4 * np.pi * alpha_mz)
    den = (s - MZ**2) ** 2 + GZ**2 * MZ**2
    chi1 = kappa * s * (s - MZ**2) / den
    chi2 = kappa**2 * s * s / den
    vl, al = -0.5 + 2 * s2w, -0.5
    num = np.zeros_like(s)
    tot = np.zeros_like(s)
    for q, wq in w.items():
        Q, T3 = (2 / 3, 0.5) if q == "u" else (-1 / 3, -0.5)
        vq, aq = T3 - 2 * Q * s2w, T3
        A0 = Q * Q - 2 * Q * vl * vq * chi1 + (al * al + vl * vl) * (aq * aq + vq * vq) * chi2
        A1 = -4 * Q * al * aq * chi1 + 8 * al * vl * aq * vq * chi2
        num += wq * A1
        tot += wq * A0
    return 3 / 8 * num / tot


def afb_with_dilution(r: dict, *, resolution: float = 1.5, m_min: float = 110.0) -> dict:
    """Fit one dilution factor D to A_FB(m) = D x quark-level prediction
    (smeared by the mass resolution) above m_min only. Below the peak, Z
    bosons that radiated a photon migrate to lower reconstructed masses and
    drag the measured asymmetry towards its value at the peak; without a
    model of that radiation the region below ~100 GeV is not fitted."""
    edges = r["edges"]
    mid = np.sqrt(edges[:-1] * edges[1:])
    fine = np.linspace(20, 400, 7601)
    pred = afb_quark_level(fine)
    ker = np.exp(-0.5 * ((fine[:, None] - fine[None, :]) / resolution) ** 2)
    smeared = (ker @ pred) / ker.sum(axis=1)
    p_mid = np.interp(mid, fine, smeared)
    ok = np.isfinite(r["afb"]) & (r["err"] > 0) & (mid > m_min)
    wts = 1 / r["err"][ok] ** 2
    D = float(np.sum(wts * r["afb"][ok] * p_mid[ok]) / np.sum(wts * p_mid[ok] ** 2))
    D_err = float(1 / np.sqrt(np.sum(wts * p_mid[ok] ** 2)))
    chi2 = float(np.sum(wts * (r["afb"][ok] - D * p_mid[ok]) ** 2))
    m0 = float(fine[np.nonzero((smeared[:-1] < 0) & (smeared[1:] >= 0) & (fine[:-1] > 60))[0][0]])
    return {"D": D, "D_err": D_err, "chi2": chi2, "ndof": int(ok.sum() - 1),
            "curve_m": fine, "curve": D * smeared, "zero_crossing_prediction": m0}
