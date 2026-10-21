"""Bolometric light curves from multi-band photometry.

At each epoch: interpolate every band's light curve (GP), remove the dust, turn
magnitudes into flux densities at the filters' pivot wavelengths, and
integrate over wavelength with the trapezoid rule. The result is kept in two
parts, because they are known to very different accuracy:

* ``L_uvoi``: ultraviolet + optical, from the bluest band up to the reddest
  optical band (I or z, < 1 micron), plus an ultraviolet tail below the bluest
  band. The tail is a straight line down to zero flux at 1600 A (with Swift
  UV data), 2000 A (from U) or 3000 A (from B, when U is missing), because
  line blanketing suppresses the UV of supernovae far below a blackbody. For
  Type II supernovae in their first
  weeks, when the photosphere is hotter than 12,000 K, a blackbody (Wien) tail
  is used instead.
* ``L_nir``: near-infrared, from the reddest optical band through J, H, Ks
  plus a Rayleigh-Jeans tail beyond Ks. Only where JHKs were observed.

Where the near-infrared was not observed, ``analysis.py`` adds it using the
infrared fraction measured at other epochs (see ``nir_fraction``). We do not
extrapolate blackbodies into the infrared: in the nebular phase a supernova's
spectrum is emission lines, not a blackbody, and the test against SN 1987A's
published UVOIR light curve showed blackbody tails off by 10-40%.

Uncertainties come from re-drawing every band's magnitude within its GP error
(Monte Carlo); the distance error is a separate, common scale factor.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import config as C
import extinction as E
import photometry as P
from gp import LightCurveGP

UV_HOT_K = 12000.0
UV_ZERO_WITH_UV = 1600.0
UV_ZERO_OPTICAL = 2000.0
#: Without U (only from B redward), the line goes to zero at 3000 A: nebular
#: spectra have little flux below ~3500 A (e.g. Dimitriadis et al. 2017).
UV_ZERO_NO_U = 3000.0
NIR_MIN_WAVE = 10000.0  # bands redder than this are "near-infrared"
T_GRID = np.geomspace(2500.0, 60000.0, 400)
H_C_OVER_K = C.H_CGS * C.C_CGS / C.K_CGS * 1e8  # hc/k in Angstrom * K


def planck_lam(wave, temp):
    """Blackbody B_lambda shape (arbitrary normalisation), wave in A; broadcasts over T."""
    x = H_C_OVER_K / (np.asarray(wave)[..., None] * np.asarray(temp))
    with np.errstate(over="ignore"):
        return 1.0 / (np.asarray(wave)[..., None] ** 5 * np.expm1(x))


def bb_integral_below(lam0, temp):
    """Int_0^{lam0} B_lambda d lambda / B_lambda(lam0), in A (vectorised in T)."""
    grid = np.linspace(0.02 * lam0, lam0, 400)
    shape = planck_lam(grid, temp) / planck_lam(np.array([lam0]), temp)[0]
    return np.trapezoid(shape, grid, axis=0)


def fit_bb_temperature(wave, flam):
    """Grid-search blackbody temperature for SEDs ``flam`` (n_samples x n_bands), log-space LSQ."""
    logf = np.log(np.clip(flam, 1e-300, None))
    logb = np.log(planck_lam(wave, T_GRID))  # (n_bands, n_T)
    best = np.empty(flam.shape[0])
    for i in range(flam.shape[0]):
        r = logf[i][:, None] - logb
        r -= r.mean(axis=0)
        best[i] = T_GRID[np.argmin(np.sum(r * r, axis=0))]
    return best


class Bolometric:
    """GP-interpolated multi-band photometry of one supernova, ready to integrate."""

    def __init__(self, sn: str, phot: pd.DataFrame, bands=None, min_points: int = 4,
                 hot_uv_days: float = 0.0, host_dust: tuple[float, float] | None = None):
        self.sn = sn
        self.cfg = dict(C.SUPERNOVAE[sn])
        if host_dust is not None:  # (E(B-V), R_V) to test another dust law
            self.cfg["ebv_host"], self.cfg["rv_host"] = host_dust
        self.hot_uv_days = hot_uv_days
        bands = bands or sorted(phot["band"].unique(), key=P.pivot)
        self.gps, self.system = {}, {}
        for b in bands:
            d = phot[(phot["band"] == b) & (phot["phase"] > 0.3)]
            if len(d) < min_points:
                continue
            self.gps[b] = LightCurveGP(d["phase"].values, d["mag"].values, d["err"].values)
            self.system[b] = d["system"].iloc[0]
        self.bands = sorted(self.gps, key=P.pivot)
        self.wave = np.array([P.pivot(b) for b in self.bands])
        self.a_band = np.array([self.extinction(w) for w in self.wave])
        self.dist_cm = 10 ** (0.2 * self.cfg["mu"] + 1.0) * C.PC_CM

    def extinction(self, wave):
        cfg = self.cfg
        return float(E.a_lambda(wave, cfg["ebv_mw"], 3.1)[0]
                     + E.a_lambda(wave, cfg["ebv_host"], cfg["rv_host"])[0])

    def _covered(self, band, t):
        # Near-infrared light curves are sparser: allow wider gaps.
        if P.pivot(band) > NIR_MIN_WAVE:
            return self.gps[band].covered(t, frac=0.2, min_days=3.0)[0]
        return self.gps[band].covered(t)[0]

    def luminosity(self, phases, n_mc: int = 200, required=("B", "V"), seed: int = 1,
                   allowed_bands=None) -> pd.DataFrame:
        """UV+optical and near-infrared luminosity (erg/s) with Monte Carlo errors.

        ``required`` bands must be covered at an epoch for it to be used;
        ``allowed_bands`` optionally restricts which bands enter.
        """
        rng = np.random.default_rng(seed)
        four_pi_d2 = 4.0 * np.pi * self.dist_cm**2
        rows = []
        for t in np.atleast_1d(phases):
            use = [i for i, b in enumerate(self.bands)
                   if self._covered(b, t) and (allowed_bands is None or b in allowed_bands)]
            names = [self.bands[i] for i in use]
            opt = [i for i in use if self.wave[i] < NIR_MIN_WAVE]
            if any(r not in names for r in required) or len(opt) < 3:
                continue
            mu = np.empty(len(use))
            sd = np.empty(len(use))
            for k, i in enumerate(use):
                m, s = self.gps[self.bands[i]].predict(t)
                mu[k], sd[k] = m[0], s[0]
            mags = mu[None, :] + sd[None, :] * rng.standard_normal((n_mc, len(use)))
            mags = np.vstack([mu, mags])  # row 0: the best estimate
            flam = np.column_stack([P.mag_to_flam(mags[:, k] - self.a_band[i], self.bands[i],
                                                  self.system[self.bands[i]])
                                    for k, i in enumerate(use)])
            wave = self.wave[use]
            out = self._integrate(wave, flam, hot_uv=t <= self.hot_uv_days)
            row = {"phase": float(t), "bands": ",".join(names),
                   "T_bb": float(out["temp"][0]), "has_uv": bool(wave.min() < 3000),
                   "has_nir": bool(out["has_nir"])}
            for key in ("uvoi", "nir", "uv_tail"):
                vals = out[key] * four_pi_d2
                row[f"L_{key}"] = float(vals[0])
                row[f"L_{key}_err"] = float(np.std(vals[1:]))
            tot = (out["uvoi"] + out["nir"]) * four_pi_d2
            row["L_tot_obs"] = float(tot[0]) if out["has_nir"] else np.nan
            row["L_tot_obs_err"] = float(np.std(tot[1:])) if out["has_nir"] else np.nan
            rows.append(row)
        return pd.DataFrame(rows)

    @staticmethod
    def _integrate(wave, flam, hot_uv: bool = False) -> dict:
        """Integrate SEDs (n_samples x n_bands, f_lambda in erg/s/cm^2/A) over wavelength."""
        opt = wave < NIR_MIN_WAVE
        nir = ~opt
        wo, fo = wave[opt], flam[:, opt]
        uvoi = np.trapezoid(fo, wo, axis=1)
        red = (wo >= 4500.0)
        temp = fit_bb_temperature(wo[red], fo[:, red]) if red.sum() >= 3 else np.full(len(flam), np.nan)

        lam_blue = wo[0]
        if hot_uv and np.isfinite(temp[0]) and temp[0] >= UV_HOT_K:
            uv = fo[:, 0] * bb_integral_below(lam_blue, temp)
        else:
            if lam_blue < 3000.0:
                lam_zero = UV_ZERO_WITH_UV
            elif lam_blue < 4000.0:
                lam_zero = UV_ZERO_OPTICAL
            else:
                lam_zero = UV_ZERO_NO_U
            uv = 0.5 * fo[:, 0] * max(lam_blue - lam_zero, 0.0)

        has_nir = bool(nir.any() and wave[nir].max() > 20000.0)
        if has_nir:
            wn = np.concatenate([[wo[-1]], wave[nir]])
            fn = np.column_stack([fo[:, -1], flam[:, nir]])
            nir_l = np.trapezoid(fn, wn, axis=1) + flam[:, -1] * wave[-1] / 3.0  # Rayleigh-Jeans tail
        else:
            nir_l = np.full(len(flam), np.nan)
        return {"uvoi": uvoi + uv, "uv_tail": uv, "nir": nir_l, "temp": temp, "has_nir": has_nir}


def nir_fraction(table: pd.DataFrame) -> pd.DataFrame:
    """Measured near-infrared fraction L_nir / (L_uvoi + L_nir) where JHKs were observed."""
    t = table[table["has_nir"]].copy()
    t["f_nir"] = t["L_nir"] / (t["L_uvoi"] + t["L_nir"])
    return t[["phase", "f_nir"]]


def apply_nir_fraction(table: pd.DataFrame, f_phase, f_value, f_err) -> pd.DataFrame:
    """Total luminosity: measured where JHKs exist, else L_uvoi / (1 - f_nir(t)).

    ``f_phase``/``f_value`` describe the adopted infrared fraction (linear in
    time between anchors, constant beyond them); ``f_err`` is its absolute
    1-sigma uncertainty, propagated into ``L_err``.
    """
    out = table.copy()
    f = np.interp(out["phase"], f_phase, f_value)
    fe = np.interp(out["phase"], f_phase, f_err)
    est = out["L_uvoi"] / (1.0 - f)
    est_err = np.hypot(out["L_uvoi_err"] / (1.0 - f), out["L_uvoi"] * fe / (1.0 - f) ** 2)
    measured = out["has_nir"].values
    out["L"] = np.where(measured, out["L_tot_obs"], est)
    out["L_err"] = np.where(measured, out["L_tot_obs_err"], est_err)
    out["nir_source"] = np.where(measured, "measured", "adopted fraction")
    out["f_nir_used"] = np.where(measured, out["L_nir"] / out["L_tot_obs"], f)
    return out
