"""Fit isochrone ages to a cluster's ridge line, with a Monte Carlo over the inputs.

For one set of inputs (distance modulus, reddening, [Fe/H], [alpha/Fe]) the
fit compares the cluster's ridge line (fiducial.py) with model isochrones of
every age from 8 to 16 Gyr. Each ridge point is compared along its own normal
direction: on the main and giant branches that is a colour difference at fixed
brightness, on the flat subgiant branch a brightness difference at fixed colour.

* Distance and reddening are *fixed* to the external values, never fitted: a
  turnoff that looks fainter must be older, not farther away.
* One free colour shift (delta c) is allowed. Model colours depend on the
  colour-temperature transformation and the convection treatment at the ~0.01-
  0.03 mag level, and the cluster's reddening is uncertain at a similar level.
  Absorbing that into one shift leaves the age to be set by brightness and
  shape: chiefly the subgiant branch's brightness, which is what a turnoff age
  means.
* A model "floor" (0.01 mag in colour, 0.02 mag in brightness) is added to the
  ridge's statistical errors: with tens of thousands of stars the ridge is
  known to ~0.002 mag, far better than any model reproduces a whole CMD.

The Monte Carlo repeats the fit for N draws of (distance, reddening, [Fe/H],
[alpha/Fe], ridge noise) and records the best age of each: the spread of those
ages is the age uncertainty, and running it with one input varied at a time
gives the error budget.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

import config as C
from fiducial import K, Ridge, turnoff_and_sgb_point
from isochrones import Extinction, IsoGrid, interpolate_in_age

#: Model floor added to the ridge's statistical errors (mag).
FLOOR_COLOUR = 0.01
FLOOR_MAG = 0.02
#: Region of the ridge used in the fit, relative to the turnoff (TO) and the
#: subgiant point 0.05 mag redder than it (SGB): from MS_DEPTH mag below the
#: turnoff to RGB_HEIGHT mag above the subgiant point.
MS_DEPTH = 2.0
RGB_HEIGHT = 1.5
#: Largest colour shift allowed (mag).
DC_MAX = 0.08
#: Residual assigned to a ridge point the model never crosses (scaled units).
MISS_PENALTY = 0.3


def fit_region(ridge: Ridge, ms_depth: float = MS_DEPTH, rgb_height: float = RGB_HEIGHT) -> np.ndarray:
    _, m_to, m_sgb = turnoff_and_sgb_point(ridge.c, ridge.m)
    idx = np.arange(len(ridge.m))
    before = idx <= ridge.i_to
    keep = np.where(before, ridge.m <= m_to + ms_depth, ridge.m >= m_sgb - rgb_height)
    return keep


@dataclass
class Draw:
    feh: float
    alpha: float
    dm: float
    ebv: float


class ClusterFitter:
    """Ages from one model set for one cluster's ridge line."""

    def __init__(self, grid: IsoGrid, ridge: Ridge, *, feh0: float, ages=None,
                 ms_depth: float = MS_DEPTH, rgb_height: float = RGB_HEIGHT,
                 colour_shift: bool = True, floor_colour: float = FLOOR_COLOUR,
                 floor_mag: float = FLOOR_MAG):
        self.grid = grid
        self.ext = Extinction(feh0, 0.4)
        a0, a1, da = C.AGE_GRID
        ages = np.round(np.arange(a0, a1 + 1e-9, da), 4) if ages is None else np.asarray(ages)
        # never extrapolate past the oldest isochrone of the grid (DSED stops at 15 Gyr)
        self.age_max = float(min(10 ** grid.log_age.max() / 1e9, a1))
        self.ages = ages[ages <= self.age_max + 1e-6]
        self.ages_coarse = np.arange(self.ages[0], self.age_max + 1e-6, 0.4)
        keep = fit_region(ridge, ms_depth, rgb_height)
        self.full_ridge = ridge
        self.ridge = ridge.subset(keep)
        self.colour_shift = colour_shift
        #: BaSTI samples its isochrones ~3x more densely than MIST along the
        #: sequence; every second point is plenty for the crossings below.
        self.thin = 2 if grid.name.startswith("BaSTI") else 1
        r = self.ridge
        self.px, self.py = K * r.c, r.m
        self.nx, self.ny = r.nx, r.ny
        self.tx, self.ty = r.ny, -r.nx  # tangent: the normal rotated back by -90 degrees
        self.sig2 = r.sig**2 + (K * floor_colour * r.nx) ** 2 + (floor_mag * r.ny) ** 2
        self.m_lo, self.m_hi = r.m.min() - 1.2, r.m.max() + 1.2

    # -- models ----------------------------------------------------------------
    def models(self, d: Draw, ages=None):
        """Apparent (F606W, F606W - F814W) of every model point: arrays (n_age, n_eep)."""
        comp = self.grid.at_composition(d.feh, d.alpha)
        if self.thin > 1:
            comp = {k: v[:, ::self.thin] for k, v in comp.items()}
        iso = interpolate_in_age(comp, self.grid.log_age, self.ages if ages is None else ages)
        av = C.R_V * d.ebv
        good = np.isfinite(iso["F606W"]) & np.isfinite(iso["F814W"])
        lt = np.where(good, iso["logTe"], 3.8)
        lg = np.where(good, iso["logg"], 4.0)
        a606 = self.ext("F606W", lt, lg, av)
        a814 = self.ext("F814W", lt, lg, av)
        m = iso["F606W"] + d.dm + a606
        c = (iso["F606W"] + a606) - (iso["F814W"] + a814)
        bad = ~good | (m < self.m_lo) | (m > self.m_hi)
        m = np.where(bad, np.nan, m)
        c = np.where(bad, np.nan, c)
        # drop EEP columns that are empty at every age
        cols = np.any(np.isfinite(m), axis=0)
        return m[:, cols], c[:, cols], {k: v[:, cols] for k, v in iso.items()}

    # -- residuals ---------------------------------------------------------------
    def residuals(self, m: np.ndarray, c: np.ndarray, dc: np.ndarray) -> np.ndarray:
        """Data-minus-model offsets along each ridge normal: (n_age, n_ridge), scaled units."""
        qx = K * (c + dc[:, None])  # (n_age, n_eep)
        qy = m
        dx = qx[:, None, :] - self.px[None, :, None]  # (n_age, n_ridge, n_eep)
        dy = qy[:, None, :] - self.py[None, :, None]
        u = dx * self.tx[None, :, None] + dy * self.ty[None, :, None]
        tau = dx * self.nx[None, :, None] + dy * self.ny[None, :, None]
        u0, u1 = u[..., :-1], u[..., 1:]
        cross = (u0 * u1 <= 0) & np.isfinite(u0) & np.isfinite(u1) & (u0 != u1)
        f = np.where(cross, u0 / np.where(cross, u0 - u1, 1.0), 0.0)
        tau_x = tau[..., :-1] + f * (tau[..., 1:] - tau[..., :-1])
        tau_x = np.where(cross, tau_x, np.inf)
        j = np.argmin(abs(tau_x), axis=2)
        best = np.take_along_axis(tau_x, j[..., None], axis=2)[..., 0]
        r = -best
        return np.where(np.isfinite(r), r, MISS_PENALTY)

    def chi2_curve(self, d: Draw, ridge_shift: np.ndarray | None = None, ages=None, models=None):
        """chi^2(age) with the colour shift optimised at each age; also the shifts."""
        if models is None:
            models = self.models(d, ages)
        m, c = models[0], models[1]
        w = 1.0 / self.sig2
        dc = np.zeros(m.shape[0])
        r = self.residuals(m, c, dc)
        if ridge_shift is not None:
            r = r + ridge_shift[None, :]
        if self.colour_shift:
            g = -K * self.nx  # d r / d(delta c) for a model parallel to the ridge
            for _ in range(2):
                step = -np.sum(w * r * g, axis=1) / np.sum(w * g * g)
                dc = np.clip(dc + step, -DC_MAX, DC_MAX)
                r = self.residuals(m, c, dc)
                if ridge_shift is not None:
                    r = r + ridge_shift[None, :]
        chi2 = np.sum(w * r**2, axis=1)
        return chi2, dc, r

    def best(self, d: Draw, ridge_shift=None) -> dict:
        """Best age: a coarse scan (0.4 Gyr), then a fine one (0.05 Gyr) around its minimum."""
        coarse = self.ages_coarse
        chi2, _, _ = self.chi2_curve(d, ridge_shift, ages=coarse)
        k = int(np.argmin(chi2))
        fine = np.arange(coarse[k] - 0.5, coarse[k] + 0.5001, 0.05)
        fine = fine[(fine >= self.ages[0] - 1e-9) & (fine <= self.age_max + 1e-9)]
        chi2, dc, _ = self.chi2_curve(d, ridge_shift, ages=fine)
        k = int(np.argmin(chi2))
        age = fine[k]
        if 0 < k < len(fine) - 1:  # parabola through the minimum
            y0, y1, y2 = chi2[k - 1:k + 2]
            den = y0 - 2 * y1 + y2
            if den > 0:
                age = fine[k] + 0.5 * (y0 - y2) / den * (fine[1] - fine[0])
        edge = (age <= self.ages[0] + 0.05) or (age >= self.age_max - 0.05)
        return {"age": float(age), "dc": float(dc[k]), "chi2": float(chi2[k]), "edge": bool(edge),
                "n": len(self.px)}

    # -- the classic subgiant-point age (a cross-check) -----------------------------
    def sgb_point_age(self, d: Draw) -> float:
        """Age from the subgiant point alone (Chaboyer et al. 1996 style)."""
        _, _, m_sgb_obs = turnoff_and_sgb_point(self.full_ridge.c, self.full_ridge.m)
        m, c, _ = self.models(d)
        m_sgb = np.array([turnoff_and_sgb_point(c[a], m[a])[2] if np.isfinite(c[a]).sum() > 5 else np.nan
                          for a in range(len(self.ages))])
        ok = np.isfinite(m_sgb)
        # the subgiant point fades monotonically with age
        return float(np.interp(m_sgb_obs, m_sgb[ok], self.ages[ok]))


def draw_inputs(inputs: dict, n: int, rng: np.random.Generator, *, vary=("dm", "ebv", "feh", "alpha"),
                dm_key: str = "dm") -> list[Draw]:
    """Monte Carlo draws of the external inputs; inputs not in ``vary`` stay at their centres.

    ``dm_key`` picks the distance: "dm" (Baumgardt & Vasiliev mean, with the
    systematic floor) or "dm_geo" (geometric distances only).
    """
    def gauss(key: str, err_key: str, label: str) -> np.ndarray:
        centre = inputs[key]
        if label not in vary:
            return np.full(n, centre)
        return centre + inputs[err_key] * rng.standard_normal(n)

    dm = gauss(dm_key, dm_key + "_err", "dm")
    ebv = np.clip(gauss("ebv", "ebv_err", "ebv"), 0.0, None)
    feh = gauss("feh", "feh_err", "feh")
    alpha = gauss("alpha", "alpha_err", "alpha")
    return [Draw(float(f), float(a), float(m), float(e)) for f, a, m, e in zip(feh, alpha, dm, ebv)]


def monte_carlo(fitter: ClusterFitter, inputs: dict, n: int, rng: np.random.Generator, *,
                vary=("dm", "ebv", "feh", "alpha", "ridge"), dm_key: str = "dm") -> dict:
    draws = draw_inputs(inputs, n, rng, vary=vary, dm_key=dm_key)
    ages, dcs, chi2s, edges = [], [], [], 0
    for d in draws:
        shift = None
        if "ridge" in vary:  # statistical noise of each ridge point, along its normal
            shift = np.sqrt(fitter.ridge.sig**2) * rng.standard_normal(len(fitter.px))
        b = fitter.best(d, shift)
        ages.append(b["age"])
        dcs.append(b["dc"])
        chi2s.append(b["chi2"])
        edges += b["edge"]
    ages = np.array(ages)
    return {
        "ages": ages,
        "dc": np.array(dcs),
        "chi2": np.array(chi2s),
        "n_edge": edges,
        "draws": np.array([[d.feh, d.alpha, d.dm, d.ebv] for d in draws]),
        "summary": summarize(ages),
    }


def summarize(x: np.ndarray) -> dict:
    p16, p50, p84 = np.percentile(x, [15.865, 50, 84.135])
    return {"median": float(p50), "minus": float(p50 - p16), "plus": float(p84 - p50),
            "mean": float(np.mean(x)), "std": float(np.std(x)),
            "p2.5": float(np.percentile(x, 2.5)), "p97.5": float(np.percentile(x, 97.5))}
