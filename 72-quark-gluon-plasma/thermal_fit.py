"""Fit the hadron resonance gas to measured yields: chemical freeze-out.

The model yield of each measured particle is

    N = V * sum_i A[obs, i] * n_i(T, mu_B, mu_Q, mu_S)

where n_i are the primordial thermal densities of all ~460 hadrons and
nuclei, A is the feed-down matrix for that measurement's convention, and V is
the volume per unit of rapidity. The free parameters are T, mu_B and V
(optionally gamma_S); mu_S and mu_Q follow from the conservation laws. The
chi^2 uses the measurement errors only (statistical and systematic in
quadrature), as in the published fits, so model systematics are studied
separately by refitting with other model choices.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass

import numpy as np
from iminuit import Minuit

import config as cfg
import datasets
import hadrons
import hrg


@dataclass(frozen=True)
class FitSettings:
    list_name: str = cfg.DEFAULT_LIST
    widths: str = cfg.DEFAULT_WIDTHS
    stats: str = "quantum"
    free_gammaS: bool = False
    free_muS: bool = False          # STAR's grand-canonical fits leave mu_S free (and mu_Q = 0)
    fix_muB: float | None = None    # MeV
    exclude: tuple = ()             # hadron states removed from the gas
    charm_pairs: float | None = None  # c-cbar pairs per unit y (None: the dataset's estimate)
    label: str = "default"


def andronic_mu_b(sqrt_s_gev: float) -> float:
    """Andronic et al. 2018 parametrisation of mu_B(sqrt s) in MeV (start values only)."""
    return cfg.ANDRONIC_2018["param_a"] / (1 + 0.288 * sqrt_s_gev)


class YieldModel:
    def __init__(self, ds: datasets.Dataset, settings: FitSettings = FitSettings()):
        self.ds, self.settings = ds, settings
        self.gas = hrg.model(settings.list_name, settings.widths, settings.stats,
                             exclude=tuple(settings.exclude))
        hl = self.gas.hl
        A = np.zeros((len(ds.yields), len(hl)))
        D = np.zeros_like(A)             # primordial (direct) part only
        for o, y in enumerate(ds.yields):
            for pdg, coef in y.targets.items():
                A[o] += coef * hadrons.feeddown_into(settings.list_name, y.feeddown, pdg)
                D[o, hl.index[pdg]] += coef
        self.A, self.D = A, D
        # Charm-hadron decays add a fixed number of pions and kaons (charm is
        # made in hard scatterings: it does not scale with T or V)
        self.c = np.zeros(len(ds.yields))
        self.charm_pairs = ds.charm_pairs if settings.charm_pairs is None else settings.charm_pairs
        if self.charm_pairs:
            per_pair = hadrons.charm_feeddown_per_pair()
            for o, y in enumerate(ds.yields):
                if y.group == "meson":
                    self.c[o] = self.charm_pairs * sum(
                        coef * per_pair.get(pdg, 0.0) for pdg, coef in y.targets.items())
        self._mu = (0.0, 0.0)            # warm start for mu_S, mu_Q

    def state(self, T: float, muB: float, gammaS: float = 1.0, muS: float | None = None):
        """Thermal state (GeV units) with the conservation laws imposed."""
        q_b = self.ds.q_over_b
        if muS is not None:
            return self.gas.state(T, muB, 0.0, muS, gammaS)
        st = self.gas.solve_mu(T, muB, q_b, gammaS, muS0=self._mu[0], muQ0=self._mu[1])
        self._mu = (st.muS, st.muQ)
        return st

    def densities(self, T: float, muB: float, gammaS: float = 1.0,
                  muS: float | None = None) -> tuple:
        st = self.state(T, muB, gammaS, muS)
        return self.A @ st.n, st


def _profile_volume(y, s, m) -> float:
    """Best volume for fixed shape: minimises sum ((y - V m)/s)^2."""
    return float(np.sum(y * m / s ** 2) / np.sum(m * m / s ** 2))


def fit(ds: datasets.Dataset, settings: FitSettings = FitSettings(), *,
        minos: bool = True, keys: tuple | None = None) -> dict:
    """Fit (T, mu_B, V [, gamma_S, mu_S]) to the dataset's fitted yields.

    ``keys`` restricts the fit to a subset of yields (e.g. nuclei only).
    """
    model = YieldModel(ds, settings)
    use = np.array([(y.fitted if keys is None else y.key in keys) for y in ds.yields])
    yv = np.array([y.value for y in ds.yields])
    sv = np.array([y.err for y in ds.yields])

    def chi2(T, muB, V, gammaS=1.0, muS=None):
        dens, _ = model.densities(T / 1e3, muB / 1e3, gammaS,
                                  None if muS is None else muS / 1e3)
        r = (yv[use] - V * dens[use] - cv[use]) / sv[use]
        return float(r @ r)

    cv = model.c
    T0 = 155.0
    muB0 = 0.0 if ds.sqrt_s > 1000 else andronic_mu_b(ds.sqrt_s)
    if settings.fix_muB is not None:
        muB0 = settings.fix_muB
    dens0, st0 = model.densities(T0 / 1e3, muB0 / 1e3)
    V0 = _profile_volume(yv[use] - cv[use], sv[use], dens0[use])

    names = ["T", "muB", "V"]
    start = [T0, muB0, V0]
    if settings.free_gammaS:
        names.append("gammaS")
        start.append(1.0)
    if settings.free_muS:
        names.append("muS")
        start.append(st0.muS * 1e3 if ds.mu_q == "zero" else muB0 / 4)

    def f(*p):
        kw = dict(zip(names, p))
        return chi2(kw["T"], kw["muB"], kw["V"], kw.get("gammaS", 1.0), kw.get("muS"))

    f.errordef = Minuit.LEAST_SQUARES
    m = Minuit(f, *start, name=names)
    m.limits["T"] = (90.0, 220.0)
    m.limits["V"] = (1e-3, None)
    m.limits["muB"] = (-200.0, 800.0)
    if settings.free_gammaS:
        m.limits["gammaS"] = (0.05, 3.0)
    if settings.fix_muB is not None:
        m.fixed["muB"] = True
    m.migrad()
    m.hesse()
    if minos and m.valid:
        try:
            m.minos(*[n for n in ("T", "muB") if not m.fixed[n]])
        except RuntimeError:
            pass
    best = dict(zip(names, [float(v) for v in m.values]))
    T, muB, V = best["T"], best["muB"], best["V"]
    dens, st = model.densities(T / 1e3, muB / 1e3, best.get("gammaS", 1.0),
                               None if "muS" not in best else best["muS"] / 1e3)
    pred = V * dens + cv
    direct = V * (model.D @ st.n)
    tot = model.gas.totals(st)
    ev = model.gas.excluded_volume(st)
    ev_factor = float(np.average(ev["suppression"][st.n > 0], weights=st.n[st.n > 0]))
    nfit = int(use.sum())
    nfree = len(names) - sum(1 for n in names if m.fixed[n])
    out = {
        "dataset": ds.name, "settings": asdict(settings), "valid": bool(m.valid),
        "chi2": float(m.fval), "n_yields": nfit, "ndf": nfit - nfree,
        "params": {}, "corr": {}, "muS": st.muS * 1e3, "muQ": st.muQ * 1e3,
        "charm_pairs": model.charm_pairs,
        "yields": [], "thermo": {
            "p_GeV_fm3": tot["p"], "e_GeV_fm3": tot["e"], "s_fm3": tot["s"],
            "n_fm3": tot["n"], "nB_fm3": tot["nB"],
            "e_over_n_GeV": tot["e"] / tot["n"], "s_over_T3": tot["s"] * cfg.HBARC ** 3 / (st.T ** 3),
        },
        "V_excluded_volume": V / ev_factor, "excluded_volume_factor": ev_factor,
    }
    for n in names:
        e = float(m.errors[n])
        rec = {"value": best[n], "error": e}
        if minos and n in m.merrors:
            me = m.merrors[n]
            rec.update({"minos_lo": float(me.lower), "minos_hi": float(me.upper)})
        out["params"][n] = rec
    if m.covariance is not None:
        cov = np.array(m.covariance)
        for i, a in enumerate(names):
            for j, b in enumerate(names):
                if j > i and cov[i, i] > 0 and cov[j, j] > 0:
                    out["corr"][f"{a}-{b}"] = float(cov[i, j] / math.sqrt(cov[i, i] * cov[j, j]))
    out["V_excluded_volume_err"] = out["params"]["V"]["error"] / ev_factor
    for o, y in enumerate(ds.yields):
        out["yields"].append({
            "key": y.key, "label": y.label, "html": y.html, "value": y.value, "err": y.err,
            "stat": y.stat, "sys": y.sys, "model": float(pred[o]), "direct": float(direct[o]),
            "decay_fraction": float(1 - direct[o] / pred[o]) if pred[o] > 0 else float("nan"),
            "charm": float(cv[o]),
            "fitted": bool(use[o]), "pull": float((y.value - pred[o]) / y.err),
            "ratio": float(y.value / pred[o]), "group": y.group, "mass": _mass(model, y),
            "g": _degeneracy(model, y), "source": y.source, "note": y.note,
        })
    return out


def _mass(model: YieldModel, y: datasets.Yield) -> float:
    hl = model.gas.hl
    return float(hl.mass[hl.index[next(iter(y.targets))]])


def _degeneracy(model: YieldModel, y: datasets.Yield) -> float:
    hl = model.gas.hl
    return float(hl.g[hl.index[next(iter(y.targets))]])


def profile_T(ds: datasets.Dataset, settings: FitSettings = FitSettings(),
              Ts: np.ndarray | None = None, keys: tuple | None = None) -> dict:
    """chi^2 minimised over mu_B and V on a grid of T (the fit's T profile)."""
    model = YieldModel(ds, settings)
    use = np.array([(y.fitted if keys is None else y.key in keys) for y in ds.yields])
    yv = np.array([y.value for y in ds.yields])[use] - model.c[use]
    sv = np.array([y.err for y in ds.yields])[use]
    Ts = np.arange(140.0, 175.01, 0.5) if Ts is None else Ts
    out_chi2, out_mu = [], []
    muB = 0.0 if ds.sqrt_s > 1000 else andronic_mu_b(ds.sqrt_s)
    for T in Ts:
        def f(mb):
            dens, _ = model.densities(T / 1e3, mb / 1e3)
            m_ = dens[use]
            V = _profile_volume(yv, sv, m_)
            r = (yv - V * m_) / sv
            return float(r @ r)
        f.errordef = Minuit.LEAST_SQUARES
        m = Minuit(f, muB, name=["muB"])
        if settings.fix_muB is not None:
            m.values["muB"] = settings.fix_muB
            m.fixed["muB"] = True
        m.migrad()
        muB = float(m.values["muB"])
        out_chi2.append(float(m.fval))
        out_mu.append(muB)
    return {"T": [float(t) for t in Ts], "chi2": out_chi2, "muB": out_mu}


def model_yields(ds: datasets.Dataset, T: float, muB: float, V: float,
                 settings: FitSettings = FitSettings()) -> np.ndarray:
    """Model yields for every entry of a dataset at given (T, mu_B in MeV, V)."""
    model = YieldModel(ds, settings)
    dens, _ = model.densities(T / 1e3, muB / 1e3)
    return V * dens + model.c
