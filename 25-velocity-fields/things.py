"""From a THINGS data cube to velocity fields.

A data cube is a stack of images of the sky, one per radio frequency
(channel). Each pixel therefore holds a spectrum of the 21-cm hydrogen line:
gas moving away from us shows up at lower frequency (higher velocity). To
build a *velocity field* -- one velocity per pixel -- you must decide which
velocity represents a line profile that is often lopsided: gas in spiral arms
and gas between them, or gas at different heights in a thick disc, adds up
along the line of sight. de Blok et al. (2008) compared five recipes and chose
a Gauss-Hermite fit with a skewness term (h3), because it follows the
dominant peak of a skewed profile instead of being dragged toward its tail.
We make four velocity fields, as they did:

* **Gauss-Hermite (h3)**: fit A exp(-w^2/2) [1 + h3 H3(w)], w = (v - v0)/sigma.
  We keep both the centre parameter v0 and the velocity where the fitted
  profile peaks; the two differ by ~sqrt(3) h3 sigma.
* **Single Gaussian**: the classical 1970s-80s recipe.
* **Peak**: the channel with the most emission, refined by a parabola.
* **First moment**: the intensity-weighted mean, as published by THINGS.

The fits use a vectorized Levenberg-Marquardt solver, so the ~10^5-10^6
spectra of a galaxy fit in a minute or two, and are cached in
``data/derived/``.
"""

from __future__ import annotations

import re
import time
from dataclasses import dataclass
from functools import cached_property

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from config import (C_KMS, DERIVED_DIR, HI_REST_HZ, MOM0_SNR_MIN, NOISE_EDGE_CHANNELS,
                    PEAK_SNR_MIN, THINGS_DIR, Galaxy)

SQ2, SQ3, SQ6 = np.sqrt(2.0), np.sqrt(3.0), np.sqrt(6.0)


# --------------------------------------------------------------------------
# Reading cubes
# --------------------------------------------------------------------------

def channel_velocities(header: fits.Header) -> np.ndarray:
    """Heliocentric velocity (optical convention, km/s) of every channel.

    THINGS cubes come in two AIPS flavours. 'FELO-HEL' (NGC 2403, NGC 2903,
    DDO 154) is an axis linear in frequency but labelled in heliocentric
    optical velocity: v = c (f0/f - 1) through each channel's frequency.
    A plain 'FREQ' axis (NGC 3198) holds the frequencies *as observed*, not in
    the Sun's frame; AIPS then records the heliocentric velocity of one
    channel (ALTRVAL at pixel ALTRPIX), and the others follow from Doppler
    ratios, (c + v) / (c + v_alt) = f_alt / f. Ignoring this shifts NGC 3198's
    velocities by 24.5 km/s, the Earth's motion on the observing date.
    Checked against the channel velocities of the THINGS global spectra.
    """
    n = header["NAXIS3"]
    pix = np.arange(1, n + 1, dtype=float)
    ctype = header["CTYPE3"].strip()
    rest = header.get("RESTFREQ", HI_REST_HZ)
    if header.get("VELREF", 2) > 256:
        raise ValueError("radio-convention velocities are not handled")
    if ctype.startswith("FREQ"):
        freq = header["CRVAL3"] + (pix - header["CRPIX3"]) * header["CDELT3"]
        if "ALTRVAL" in header:
            v_alt = header["ALTRVAL"] / 1e3
            f_alt = header["CRVAL3"] + (header["ALTRPIX"] - header["CRPIX3"]) * header["CDELT3"]
            return (C_KMS + v_alt) * f_alt / freq - C_KMS
    elif ctype.startswith("FELO"):
        v_ref = header["CRVAL3"] / 1e3  # m/s -> km/s
        dv_ref = header["CDELT3"] / 1e3
        f_ref = rest / (1 + v_ref / C_KMS)
        df = -dv_ref * f_ref / (C_KMS + v_ref)  # the frequency step at the reference pixel
        freq = f_ref + (pix - header["CRPIX3"]) * df
    else:
        raise ValueError(f"unsupported spectral axis {ctype!r}")
    return C_KMS * (rest / freq - 1.0)


def clean_beam(header: fits.Header) -> tuple[float, float, float]:
    """(BMAJ, BMIN in arcsec, BPA in deg) of the restoring beam.

    AIPS writes it into the history ('CLEAN BMAJ= ...'); the last one wins.
    """
    if "BMAJ" in header:
        return header["BMAJ"] * 3600, header["BMIN"] * 3600, header["BPA"]
    found = None
    for card in header.get("HISTORY", []):
        m = re.search(r"BMAJ=\s*([\d.E+-]+)\s+BMIN=\s*([\d.E+-]+)\s+BPA=\s*([\d.E+-]+)", str(card))
        if m:
            found = m
    if not found:
        raise KeyError("no clean beam in header")
    return float(found.group(1)) * 3600, float(found.group(2)) * 3600, float(found.group(3))


@dataclass
class Cube:
    galaxy: Galaxy
    path: str

    @cached_property
    def _hdul(self):
        return fits.open(self.path, memmap=True, do_not_scale_image_data=True)

    @cached_property
    def header(self) -> fits.Header:
        return self._hdul[0].header

    @cached_property
    def data(self) -> np.ndarray:
        """(channel, y, x) in Jy/beam, memory-mapped; ascending velocity."""
        d = self._hdul[0].data
        d = d.reshape(d.shape[-3:])
        return d[::-1] if self._descending else d

    @cached_property
    def _descending(self) -> bool:
        v = channel_velocities(self.header)
        return bool(v[1] < v[0])

    @cached_property
    def velocity(self) -> np.ndarray:
        v = channel_velocities(self.header)
        return v[::-1] if self._descending else v

    @property
    def channel_width(self) -> float:
        return float(np.median(np.abs(np.diff(self.velocity))))

    @cached_property
    def wcs(self) -> WCS:
        return WCS(self.header).celestial

    @property
    def pixel_arcsec(self) -> float:
        return abs(self.header["CDELT2"]) * 3600

    @cached_property
    def beam(self) -> tuple[float, float, float]:
        return clean_beam(self.header)

    @property
    def beam_fwhm_arcsec(self) -> float:
        """Geometric-mean beam FWHM."""
        return float(np.sqrt(self.beam[0] * self.beam[1]))

    @property
    def beam_area_pixels(self) -> float:
        return float(np.pi / (4 * np.log(2)) * self.beam[0] * self.beam[1] / self.pixel_arcsec ** 2)

    @cached_property
    def noise(self) -> float:
        """Channel noise (Jy/beam): robust scatter of the line-free end channels.

        The cubes have no primary-beam correction, so the noise is the same
        everywhere; 1.4826 x median absolute deviation estimates it while
        ignoring any stray emission.
        """
        k = NOISE_EDGE_CHANNELS // 2
        edge = np.concatenate([np.asarray(self.data[:k]).ravel(), np.asarray(self.data[-k:]).ravel()])
        edge = edge[np.isfinite(edge) & (edge != 0)]
        return float(1.4826 * np.median(np.abs(edge - np.median(edge))))

    def channel_noise_profile(self) -> np.ndarray:
        """Robust noise of every channel, from pixels far from the galaxy."""
        mom0 = load_moment(self.galaxy, 0)
        empty = mom0 == 0
        out = []
        for k in range(self.data.shape[0]):
            plane = np.asarray(self.data[k])[empty]
            plane = plane[np.isfinite(plane) & (plane != 0)]
            out.append(1.4826 * np.median(np.abs(plane - np.median(plane))))
        return np.array(out)


def load_cube(galaxy: Galaxy) -> Cube:
    return Cube(galaxy, str(THINGS_DIR / galaxy.cube_file()))


def load_moment(galaxy: Galaxy, moment: int) -> np.ndarray:
    """A THINGS moment map, 2-D. Moment 0 in Jy/beam km/s, 1 and 2 in km/s."""
    with fits.open(THINGS_DIR / galaxy.moment_file(moment)) as hdul:
        data = np.array(hdul[0].data, dtype=float).reshape(hdul[0].data.shape[-2:])
        unit = hdul[0].header.get("BUNIT", "").upper()
    if "M/S" in unit or "METR/SEC" in unit:
        data = data / 1e3
    return data


# --------------------------------------------------------------------------
# Profile fitting
# --------------------------------------------------------------------------

def gauss_hermite(v: np.ndarray, amp, v0, sigma, h3) -> np.ndarray:
    """Gauss-Hermite series truncated at h3 (van der Marel & Franx 1993)."""
    w = (v - v0) / sigma
    h3poly = (2 * SQ2 * w ** 3 - 3 * SQ2 * w) / SQ6
    return amp * np.exp(-0.5 * w * w) * (1 + h3 * h3poly)


def _model_and_jacobian(v: np.ndarray, p: np.ndarray, with_h3: bool):
    """Model (n, m) and Jacobian (n, m, k) for n spectra on m channels."""
    amp, v0, sig = p[:, 0:1], p[:, 1:2], p[:, 2:3]
    h3 = p[:, 3:4] if with_h3 else 0.0
    w = (v[None, :] - v0) / sig
    e = np.exp(-0.5 * w * w)
    poly = (2 * SQ2 * w ** 3 - 3 * SQ2 * w) / SQ6
    dpoly = (6 * SQ2 * w ** 2 - 3 * SQ2) / SQ6
    shape = 1 + h3 * poly
    model = amp * e * shape
    j_amp = e * shape
    j_v0 = amp * e * (w * shape - h3 * dpoly) / sig
    j_sig = amp * e * (w * w * shape - h3 * dpoly * w) / sig
    cols = [j_amp, j_v0, j_sig] + ([amp * e * poly] if with_h3 else [])
    return model, np.stack(cols, axis=-1)


def fit_profiles(v: np.ndarray, spectra: np.ndarray, p0: np.ndarray, *, with_h3: bool,
                 noise: float, n_iter: int = 60) -> dict[str, np.ndarray]:
    """Least-squares fit of many profiles at once (Levenberg-Marquardt).

    ``spectra`` is (n, m); ``p0`` holds starting values (amp, v0, sigma[, h3]).
    Bounds keep the fit physical: sigma between 0.3 channel and half the band,
    |h3| <= 0.6, v0 inside the band. Returns parameters, their 1-sigma errors
    (from the curvature matrix and the channel noise), chi^2 and a flag for
    fits that ended on a bound.
    """
    k = 4 if with_h3 else 3
    p = p0[:, :k].astype(float).copy()
    dv = float(np.median(np.abs(np.diff(v))))
    lo = np.array([0.0, v.min(), 0.3 * dv, -0.6][:k])
    hi = np.array([np.inf, v.max(), 0.5 * (v.max() - v.min()), 0.6][:k])
    lam = np.full(len(p), 1e-2)
    model, jac = _model_and_jacobian(v, p, with_h3)
    chi2 = np.sum((spectra - model) ** 2, axis=1)
    eye = np.eye(k)
    for _ in range(n_iter):
        r = spectra - model
        jtj = np.einsum("nmi,nmj->nij", jac, jac)
        jtr = np.einsum("nmi,nm->ni", jac, r)
        diag = np.einsum("nii->ni", jtj)
        a = jtj + lam[:, None, None] * diag[:, :, None] * eye
        try:
            step = np.linalg.solve(a, jtr[..., None])[..., 0]
        except np.linalg.LinAlgError:
            a += 1e-12 * eye
            step = np.linalg.lstsq(a.reshape(-1, k), jtr.ravel(), rcond=None)[0].reshape(-1, k)
        trial = np.clip(p + step, lo, hi)
        m_trial, j_trial = _model_and_jacobian(v, trial, with_h3)
        chi2_trial = np.sum((spectra - m_trial) ** 2, axis=1)
        better = chi2_trial < chi2
        p[better], model[better], jac[better], chi2[better] = (trial[better], m_trial[better],
                                                               j_trial[better], chi2_trial[better])
        lam = np.where(better, lam / 3, lam * 4)
        lam = np.clip(lam, 1e-7, 1e7)
    jtj = np.einsum("nmi,nmj->nij", jac, jac)
    try:
        cov = np.linalg.inv(jtj + 1e-30 * eye)
        err = np.sqrt(np.clip(np.einsum("nii->ni", cov), 0, None)) * noise
    except np.linalg.LinAlgError:
        err = np.full_like(p, np.nan)
    at_bound = np.any(np.isclose(p, lo, rtol=0, atol=1e-6) | np.isclose(p, hi, rtol=0, atol=1e-6), axis=1)
    return {"params": p, "errors": err, "chi2": chi2 / noise ** 2, "at_bound": at_bound}


def profile_peak(v0, sigma, h3, amp) -> np.ndarray:
    """Velocity where a fitted Gauss-Hermite profile peaks (grid + parabola)."""
    t = np.linspace(-2.5, 2.5, 251)
    w = t[None, :]
    poly = (2 * SQ2 * w ** 3 - 3 * SQ2 * w) / SQ6
    shape = np.exp(-0.5 * w * w) * (1 + h3[:, None] * poly)
    i = np.clip(np.argmax(shape, axis=1), 1, len(t) - 2)
    rows = np.arange(len(i))
    y0, y1, y2 = shape[rows, i - 1], shape[rows, i], shape[rows, i + 1]
    denom = y0 - 2 * y1 + y2
    frac = np.where(denom != 0, 0.5 * (y0 - y2) / denom, 0.0)
    return v0 + sigma * (t[i] + frac * (t[1] - t[0]))


def _initial_guesses(v: np.ndarray, spectra: np.ndarray) -> np.ndarray:
    """Peak, centroid and width of the strongest feature in each spectrum."""
    dv = float(np.median(np.diff(v)))
    i = np.argmax(spectra, axis=1)
    amp = spectra[np.arange(len(i)), i]
    window = np.abs(np.arange(spectra.shape[1])[None, :] - i[:, None]) <= 4
    s = np.where(window, np.clip(spectra, 0, None), 0.0)
    tot = s.sum(axis=1) + 1e-30
    cen = (s * v[None, :]).sum(axis=1) / tot
    var = (s * (v[None, :] - cen[:, None]) ** 2).sum(axis=1) / tot
    sig = np.clip(np.sqrt(np.clip(var, 0, None)), dv, 4 * dv)
    return np.column_stack([amp, v[i], sig, np.zeros_like(amp)])


def peak_velocity(v: np.ndarray, spectra: np.ndarray) -> np.ndarray:
    """Velocity of the brightest channel, refined by a parabola through it."""
    i = np.clip(np.argmax(spectra, axis=1), 1, spectra.shape[1] - 2)
    rows = np.arange(len(i))
    y0, y1, y2 = spectra[rows, i - 1], spectra[rows, i], spectra[rows, i + 1]
    denom = y0 - 2 * y1 + y2
    frac = np.where(denom < 0, 0.5 * (y0 - y2) / np.where(denom == 0, 1, denom), 0.0)
    return v[i] + np.clip(frac, -0.5, 0.5) * (v[1] - v[0])


# --------------------------------------------------------------------------
# Whole-galaxy velocity fields
# --------------------------------------------------------------------------

FIELDS = ("gh_amp", "gh_v0", "gh_sigma", "gh_h3", "gh_vpeak", "gh_v0_err", "gh_chi2",
          "gauss_amp", "gauss_v", "gauss_sigma", "gauss_v_err", "peak_v", "peak_snr",
          "mom1_things", "mom0_things", "n_line_channels")


def fit_cube(data, v: np.ndarray, noise: float, candidates: np.ndarray) -> dict[str, np.ndarray]:
    """Fit every candidate spectrum of a (channel, y, x) cube, in row blocks.

    Each spectrum gets a single-Gaussian fit, then a Gauss-Hermite fit started
    from it. Works on a memory-mapped file or an in-memory array (the model
    cubes of the beam-smearing test go through exactly this code).
    """
    dv = float(np.median(np.abs(np.diff(v))))
    ny, nx = candidates.shape
    out = {k: np.full((ny, nx), np.nan, dtype=np.float32) for k in FIELDS}
    block = 32 if nx <= 1024 else 16
    for y0 in range(0, ny, block):
        rows = candidates[y0:y0 + block]
        if not rows.any():
            continue
        slab = np.asarray(data[:, y0:y0 + block, :], dtype=np.float64)  # (m, b, nx)
        yy, xx = np.nonzero(rows)
        spectra = np.nan_to_num(slab[:, yy, xx].T)  # (n, m)
        p0 = _initial_guesses(v, spectra)
        g = fit_profiles(v, spectra, p0, with_h3=False, noise=noise)
        p0h = np.column_stack([g["params"], np.zeros(len(spectra))])
        h = fit_profiles(v, spectra, p0h, with_h3=True, noise=noise)
        gp, hp = g["params"], h["params"]
        Y, X = yy + y0, xx
        out["gauss_amp"][Y, X], out["gauss_v"][Y, X], out["gauss_sigma"][Y, X] = gp.T
        out["gauss_v_err"][Y, X] = g["errors"][:, 1]
        out["gh_amp"][Y, X], out["gh_v0"][Y, X], out["gh_sigma"][Y, X], out["gh_h3"][Y, X] = hp.T
        out["gh_v0_err"][Y, X] = h["errors"][:, 1]
        out["gh_chi2"][Y, X] = h["chi2"] / max(len(v) - 4, 1)
        out["gh_vpeak"][Y, X] = profile_peak(hp[:, 1], hp[:, 2], hp[:, 3], hp[:, 0])
        bad = h["at_bound"]
        out["gh_vpeak"][Y[bad], X[bad]] = np.nan
        out["gh_v0"][Y[bad], X[bad]] = np.nan
        bad_g = g["at_bound"]
        out["gauss_v"][Y[bad_g], X[bad_g]] = np.nan
        out["peak_v"][Y, X] = peak_velocity(v, spectra)
        out["peak_snr"][Y, X] = spectra.max(axis=1) / noise
        # channels within +-2 sigma of the fitted line: the "N" of the moment-0 noise
        out["n_line_channels"][Y, X] = np.clip(np.round(4 * hp[:, 2] / dv) + 1, 3, len(v))
    return out


def build_velocity_fields(galaxy: Galaxy, *, quiet: bool = False) -> dict[str, np.ndarray]:
    """Fit every spectrum inside the THINGS emission region; return maps.

    Candidates are pixels where the survey's own masked moment-0 map has
    emission. Quality masks are applied later (``velocity_field``), so the
    raw fits can be inspected.
    """
    cube = load_cube(galaxy)
    mom0 = load_moment(galaxy, 0)
    t0 = time.time()
    # Fit in chunks of rows, checkpointing each, so an interrupted run (the
    # biggest cube takes ~50 min) resumes where it stopped.
    candidates = mom0 > 0
    part = DERIVED_DIR / f"{galaxy.key}_profile_fits.partial.npz"
    chunk = 256
    if part.exists():
        with np.load(part) as f:
            out = {k: f[k] for k in FIELDS}
            done_rows = int(f["done_rows"])
    else:
        out = {k: np.full(mom0.shape, np.nan, dtype=np.float32) for k in FIELDS}
        done_rows = 0
    for y0 in range(done_rows, mom0.shape[0], chunk):
        mask = np.zeros_like(candidates)
        mask[y0:y0 + chunk] = candidates[y0:y0 + chunk]
        if mask.any():
            got = fit_cube(cube.data, cube.velocity, cube.noise, mask)
            for k in FIELDS:
                out[k][y0:y0 + chunk] = got[k][y0:y0 + chunk]
        np.savez(part, done_rows=min(y0 + chunk, mom0.shape[0]), **out)
    out["mom1_things"][:] = load_moment(galaxy, 1)
    out["mom0_things"][:] = mom0
    if not quiet:
        print(f"    {galaxy.name}: fitted {int((mom0 > 0).sum()):,} spectra in {time.time() - t0:.0f} s "
              f"(noise {cube.noise * 1e3:.3f} mJy/beam, channels {cube.channel_width:.2f} km/s)", flush=True)
    return out


def velocity_fields(galaxy: Galaxy, *, rebuild: bool = False, quiet: bool = False) -> dict[str, np.ndarray]:
    """Cached profile fits for one galaxy (see ``build_velocity_fields``)."""
    path = DERIVED_DIR / f"{galaxy.key}_profile_fits.npz"
    if path.exists() and not rebuild:
        with np.load(path) as f:
            return {k: f[k] for k in f.files}
    maps = build_velocity_fields(galaxy, quiet=quiet)
    np.savez_compressed(path, **maps)
    (DERIVED_DIR / f"{galaxy.key}_profile_fits.partial.npz").unlink(missing_ok=True)
    return maps


def quality_mask(maps: dict[str, np.ndarray], cube, *, kind: str = "gh") -> np.ndarray:
    """de Blok et al. (2008, Sect. 3.4) criteria for a usable fitted profile.

    * fitted peak intensity above 3 x the channel noise;
    * fitted dispersion larger than the channel spacing;
    * integrated HI (THINGS moment 0) above 3 sigma_N, sigma_N = sqrt(N) x
      channel noise x channel width, N the number of channels the line covers.
    """
    noise, dv = cube.noise, cube.channel_width
    amp = maps["gh_amp"] if kind == "gh" else maps["gauss_amp"]
    sig = maps["gh_sigma"] if kind == "gh" else maps["gauss_sigma"]
    sigma_n = np.sqrt(maps["n_line_channels"]) * noise * dv
    ok = np.isfinite(amp) & (amp > PEAK_SNR_MIN * noise) & (sig > dv)
    ok &= maps["mom0_things"] > MOM0_SNR_MIN * sigma_n
    return ok


def velocity_field(galaxy: Galaxy, kind: str = "gh_vpeak", maps=None, cube=None) -> np.ndarray:
    """One velocity field (km/s, NaN where unusable).

    ``kind``: 'gh_vpeak' (Gauss-Hermite profile peak, our default), 'gh_v0'
    (its centre parameter), 'gauss_v', 'peak_v', or 'mom1_things'.
    """
    maps = maps if maps is not None else velocity_fields(galaxy)
    cube = cube if cube is not None else load_cube(galaxy)
    vf = maps[kind].astype(float).copy()
    if kind == "mom1_things":
        return vf  # the survey's own map, with its own mask
    mask = quality_mask(maps, cube, kind="gauss" if kind == "gauss_v" else "gh")
    if kind == "peak_v":
        mask &= maps["peak_snr"] > PEAK_SNR_MIN
    vf[~mask] = np.nan
    return vf
