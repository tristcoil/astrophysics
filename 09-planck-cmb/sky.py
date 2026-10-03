"""The Planck sky maps: from 50 million pixels to a power spectrum.

We use Planck's SMICA CMB map -- the microwave sky with the Galaxy and other
foregrounds removed by combining Planck's nine frequency channels -- in its two
**half-mission** versions: the same sky, made from the first and second halves
of Planck's observations. The sky is identical in both; the instrument noise is
independent. That one fact does a lot of work:

* ``(hm1 + hm2) / 2`` is the sky map.
* ``(hm1 - hm2) / 2`` cancels the sky exactly and leaves a map of pure noise,
  with the same noise level as the sum map.
* Correlating hm1 with hm2 (a *cross*-spectrum) measures the sky's power with
  no contribution from noise at all, because the two noises do not correlate.

The pipeline:

1. Read both maps; build the analysis mask (Planck's common confidence mask,
   which removes the Galactic plane and bright point sources, times the mask of
   pixels unobserved in either half) and taper its edges.
2. Remove the monopole and dipole fitted outside the mask -- the dipole is our
   own motion through the CMB, not a property of the early universe.
3. Spherical-harmonic transform each masked half map; form the cross, sum and
   difference spectra.
4. Undo the cut sky -- not just its area but the way it mixes neighbouring
   multipoles, via correlation functions (``cmb_physics.deconvolve_mask``) --
   and the smoothing by the beam and the pixels (divide by B_l^2 p_l^2).
5. For contrast, also compute the naive spectrum: the whole sky, no mask, no
   corrections.
6. Check the estimator on a *simulated* sky with a known spectrum, pushed
   through the identical mask, beam and code. (Labelled as a simulation
   wherever it appears; it tests the method, it is not data.)

Everything expensive is cached in ``data/derived/`` so figures regenerate in
seconds.
"""

from __future__ import annotations

import time

import numpy as np

import cmb_physics as cp
from config import (APODIZATION_FWHM_DEG, BEAM_PATH, DERIVED_DIR, DISPLAY_NSIDE, LMAX, NSIDE,
                    RANDOM_SEED, T_CMB_K, product)

K_TO_UK = 1e6
SPECTRA_CACHE = DERIVED_DIR / "map_spectra.npz"
DISPLAY_CACHE = DERIVED_DIR / "display_maps.npz"
SIM_CACHE = DERIVED_DIR / "estimator_validation.npz"
#: Lowest multipole used when fitting the point-source excess.
POISSON_FIT_LMIN = 1000

#: Centre of the zoom-in panel, in Galactic coordinates (l, b) degrees. Chosen
#: far from the Galactic plane, clear of the Magellanic Clouds, in a patch the
#: mask leaves untouched.
ZOOM_CENTER = (230.0, -62.0)
ZOOM_SIZE_DEG = 20.0
ZOOM_PIXELS = 800


def _hp():
    import healpy as hp
    return hp


def _log(msg: str, t0: float | None = None) -> None:
    extra = f"  ({time.monotonic() - t0:.1f} s)" if t0 is not None else ""
    print(f"    {msg}{extra}", flush=True)


# --------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------


def read_half_missions():
    """Both half-mission temperature maps in microkelvin, RING ordering.

    Planck stores maps in NESTED pixel ordering (good for pixel neighbours);
    healpy's harmonic transforms want RING (pixels in rings of constant
    latitude), so ``read_map`` reorders on the fly.
    """
    hp = _hp()
    maps = []
    for key in ("hm1", "hm2"):
        m = hp.read_map(product(key).path, field=0, dtype=np.float32)
        maps.append(m * np.float32(K_TO_UK))
    return maps


def read_beam(lmax: int = LMAX) -> np.ndarray:
    """Total transfer function from sky to map: SMICA's beam times the pixel window.

    Two effects smooth the map. The first is SMICA's effective beam B_l, read
    from Planck's own file -- within 0.1% of a 5-arcmin Gaussian below l = 2000.
    The second is the HEALPix pixel window p_l: each pixel reports an average
    over its 1.7-arcmin area rather than a value at a point.

    Whether the second one applies is a genuine trap. The FITS header note says
    the beam "includes the pixel window function", which would mean dividing by
    B_l alone. But the values in the file are a bare 5' Gaussian, the same note
    writes the relation as C_l(map) = C_l(sky) B_l^2 p_l^2, and the 2018 product
    description is explicit: the map "is not deconvolved from the pixel window
    function". Our own data settle it independently: without p_l, the recovered
    spectrum falls short of Planck's by 2.2% at l = 1000 -- which is exactly
    1 - p_l^2 there. So both are divided out.
    """
    from astropy.io import fits

    with fits.open(BEAM_PATH) as hdul:
        b = np.asarray(hdul[1].data["INT_BEAM"], float)
    return b[: lmax + 1] * _hp().pixwin(NSIDE, lmax=lmax)


def binary_mask() -> np.ndarray:
    """Pixels we trust: common confidence mask AND observed in both halves."""
    hp = _hp()
    common = hp.read_map(product("mask_common").path, dtype=np.float32)
    observed = hp.read_map(product("mask_hm").path, dtype=np.float32)
    return ((common > 0.5) & (observed > 0.5)).astype(np.float32)


def apodize(mask: np.ndarray, fwhm_deg: float = APODIZATION_FWHM_DEG) -> np.ndarray:
    """Taper a binary mask smoothly to zero at its edges.

    A hard-edged mask is a step function, and a step has power at every
    multipole; multiplying the sky by it smears power from large scales into
    small ones. We therefore (1) erode the mask -- smooth it with a Gaussian and
    keep only pixels still above 0.99, i.e. more than ~2.3 sigma from any edge
    -- and then (2) smooth the eroded mask with the same Gaussian. The result
    rises from ~0 at the original edge to 1 about a beam-width inside it, so the
    taper sits entirely within sky we trust. Multiplying by the original mask
    guarantees exact zeros everywhere we do not.
    """
    hp = _hp()
    lmax_s = int(4.5 * np.sqrt(8 * np.log(2)) / np.radians(fwhm_deg))  # beam < 1e-8 beyond
    fwhm = np.radians(fwhm_deg)
    smooth = hp.smoothing(mask.astype(np.float64), fwhm=fwhm, lmax=lmax_s, iter=0)
    eroded = (smooth > 0.99).astype(np.float64)
    del smooth
    tapered = hp.smoothing(eroded, fwhm=fwhm, lmax=lmax_s, iter=0)
    return (np.clip(tapered, 0.0, 1.0) * mask).astype(np.float32)


def remove_monopole_dipole(m: np.ndarray, mask: np.ndarray) -> tuple[float, np.ndarray]:
    """Fit monopole + dipole on the unmasked sky and subtract them everywhere.

    Returns the monopole and the dipole vector (uK) that were removed. The SMICA
    maps already had the large kinematic dipole (3.4 mK, our motion through the
    CMB) taken out; what remains is a small residual, but a residual dipole
    outside a mask leaks into every multipole of the masked spectrum, so it has
    to go.
    """
    hp = _hp()
    fit_map = np.where(mask > 0.5, m, hp.UNSEEN)
    mono, dip = hp.fit_dipole(fit_map, bad=hp.UNSEEN)
    del fit_map
    npix = m.size
    chunk = npix // 32
    for start in range(0, npix, chunk):
        ipix = np.arange(start, min(start + chunk, npix))
        x, y, z = hp.pix2vec(NSIDE, ipix)
        m[ipix] -= np.float32(mono) + (dip[0] * x + dip[1] * y + dip[2] * z).astype(np.float32)
    return float(mono), np.asarray(dip, float)


# --------------------------------------------------------------------------
# Spectra
# --------------------------------------------------------------------------


def mask_power_spectrum(apod: np.ndarray) -> np.ndarray:
    """Power spectrum of the mask itself, to twice the analysis lmax.

    The mask's correlation function must be known more finely than the map's,
    because it multiplies the sky: a product of two functions band-limited at
    lmax contains multipoles up to 2 lmax.
    """
    hp = _hp()
    return hp.alm2cl(hp.map2alm(apod.astype(np.float64), lmax=2 * LMAX, iter=3))


def _debeam(cl: np.ndarray, beam: np.ndarray) -> np.ndarray:
    """Divide out the beam, leaving the monopole and dipole (where B_0 = 0) at zero."""
    out = np.zeros_like(cl)
    ok = beam > 0
    out[ok] = cl[ok] / beam[ok] ** 2
    out[:2] = 0.0
    return out


def compute_map_spectra(*, force: bool = False) -> dict:
    """Cross, auto, noise and naive spectra of the SMICA half-mission maps.

    All spectra are returned as C_l in uK^2 for l = 0..LMAX. The corrected
    ones have had the cut sky and the beam removed; ``naive_raw`` has had
    nothing done to it at all.
    """
    if SPECTRA_CACHE.exists() and not force:
        with np.load(SPECTRA_CACHE, allow_pickle=False) as z:
            return {k: z[k] for k in z.files}

    hp = _hp()
    print("  computing spectra from the sky maps (cached afterwards)")
    t0 = time.monotonic()
    hm1, hm2 = read_half_missions()
    _log("read two half-mission maps", t0)

    # The naive spectrum first, from the untouched sum map: whole sky, raw.
    t = time.monotonic()
    alm_naive = hp.map2alm(0.5 * (hm1.astype(np.float64) + hm2), lmax=LMAX, iter=3)
    naive_raw = hp.alm2cl(alm_naive)
    del alm_naive
    _log("naive full-sky spectrum", t)

    t = time.monotonic()
    mask = binary_mask()
    apod = apodize(mask)
    moments_binary = cp.mask_moments(mask)
    moments = cp.mask_moments(apod)
    _log(f"mask: f_sky {moments_binary['fsky']:.3f} binary, "
         f"w2 {moments['w2']:.3f} after tapering", t)

    t = time.monotonic()
    mono1, dip1 = remove_monopole_dipole(hm1, mask)
    mono2, dip2 = remove_monopole_dipole(hm2, mask)
    _log(f"removed monopole/dipole: residual dipole amplitudes "
         f"{np.linalg.norm(dip1):.2f} and {np.linalg.norm(dip2):.2f} uK", t)

    t = time.monotonic()
    alm1 = hp.map2alm(apod * hm1.astype(np.float64), lmax=LMAX, iter=3)
    alm2 = hp.map2alm(apod * hm2.astype(np.float64), lmax=LMAX, iter=3)
    _log("masked harmonic transforms", t)

    beam = read_beam()
    pseudo_cross = hp.alm2cl(alm1, alm2)
    pseudo_sum = hp.alm2cl(0.5 * (alm1 + alm2))
    pseudo_diff = hp.alm2cl(0.5 * (alm1 - alm2))
    del alm1, alm2

    t = time.monotonic()
    mask_cl = mask_power_spectrum(apod)
    cross, auto, noise = cp.deconvolve_mask([pseudo_cross, pseudo_sum, pseudo_diff], mask_cl,
                                            lmax=LMAX)
    _log("mask coupling removed (correlation-function method)", t)

    # Map statistics on the trusted sky, at full resolution.
    good = mask > 0.5
    sum_map = 0.5 * (hm1 + hm2)
    diff_map = 0.5 * (hm1 - hm2)
    out = {
        "ell": np.arange(LMAX + 1),
        "beam": beam,
        "cross": _debeam(cross, beam),
        "auto": _debeam(auto, beam),
        "noise": _debeam(noise, beam),
        # The simple f_sky correction, kept to show what the full one fixes.
        "cross_fsky": _debeam(pseudo_cross / moments["w2"], beam),
        "naive_raw": naive_raw,
        "pseudo_cross": pseudo_cross,
        "mask_cl": mask_cl,
        "w1": moments["w1"], "w2": moments["w2"], "w4": moments["w4"],
        "fsky_eff": moments["fsky_eff"],
        "fsky_binary": moments_binary["fsky"],
        "fsky_common": float(np.mean(hp.read_map(product("mask_common").path,
                                                  dtype=np.float32) > 0.5)),
        "monopole_removed": np.array([mono1, mono2]),
        "dipole_removed": np.vstack([dip1, dip2]),
        "rms_sum_uK": float(np.std(sum_map[good])),
        "rms_diff_uK": float(np.std(diff_map[good])),
        "rms_sum_fullsky_uK": float(np.std(sum_map)),
        "min_sum_uK": float(sum_map[good].min()),
        "max_sum_uK": float(sum_map[good].max()),
        "max_abs_fullsky_uK": float(np.max(np.abs(sum_map))),
    }
    np.savez(SPECTRA_CACHE, **out)
    _log(f"cached to {SPECTRA_CACHE.relative_to(DERIVED_DIR.parent.parent)}", t0)
    return out


def compute_display_maps(*, force: bool = False) -> dict:
    """Low-resolution maps for the whole-sky figures, and a full-resolution zoom.

    The whole-sky maps are averaged down to DISPLAY_NSIDE (6.9 arcmin pixels) --
    about one pixel of the printed figure -- so the image shows what is there
    rather than aliased pixel noise. The zoom keeps Planck's full resolution.
    """
    if DISPLAY_CACHE.exists() and not force:
        with np.load(DISPLAY_CACHE, allow_pickle=False) as z:
            return {k: z[k] for k in z.files}

    hp = _hp()
    print("  building display maps (cached afterwards)")
    t0 = time.monotonic()
    hm1, hm2 = read_half_missions()
    mask = binary_mask()
    remove_monopole_dipole(hm1, mask)
    remove_monopole_dipole(hm2, mask)
    sum_map = 0.5 * (hm1 + hm2)
    diff_map = 0.5 * (hm1 - hm2)

    proj = hp.projector.GnomonicProj(rot=(ZOOM_CENTER[0], ZOOM_CENTER[1], 0.0),
                                     xsize=ZOOM_PIXELS, ysize=ZOOM_PIXELS,
                                     reso=ZOOM_SIZE_DEG * 60.0 / ZOOM_PIXELS)
    to_pix = lambda x, y, z: hp.vec2pix(NSIDE, x, y, z)  # noqa: E731
    zoom_sum = proj.projmap(sum_map, to_pix)
    zoom_diff = proj.projmap(diff_map, to_pix)
    zoom_mask = proj.projmap(mask, to_pix)

    out = {
        "sum": hp.ud_grade(sum_map, DISPLAY_NSIDE).astype(np.float32),
        "diff": hp.ud_grade(diff_map, DISPLAY_NSIDE).astype(np.float32),
        "mask": hp.ud_grade(mask, DISPLAY_NSIDE).astype(np.float32),
        "zoom_sum": np.asarray(zoom_sum, np.float32),
        "zoom_diff": np.asarray(zoom_diff, np.float32),
        "zoom_mask_fraction": float(np.mean(zoom_mask)),
        "zoom_center": np.array(ZOOM_CENTER),
        "zoom_size_deg": ZOOM_SIZE_DEG,
    }
    np.savez(DISPLAY_CACHE, **out)
    _log("display maps cached", t0)
    return out


def validate_estimator(theory_ell, theory_dl, *, force: bool = False) -> dict:
    """Run the estimator on a simulated sky whose spectrum we know. (SIMULATION.)

    A Gaussian random sky is drawn from Planck's best-fit spectrum, blurred by
    the SMICA beam, and pushed through exactly the same mask, taper and
    corrections as the real maps. We compare the result with the full-sky
    spectrum *of the same simulated sky*, so that cosmic variance -- which would
    otherwise dominate a single realisation -- largely cancels, and what is left
    is the error of our cut-sky correction itself.

    The simulated sky contains no noise and no foregrounds: it isolates the
    geometry and the harmonic transforms, which is what is being tested.
    """
    if SIM_CACHE.exists() and not force:
        with np.load(SIM_CACHE, allow_pickle=False) as z:
            return {k: z[k] for k in z.files}

    hp = _hp()
    print("  validating the estimator on a simulated sky (cached afterwards)")
    t0 = time.monotonic()
    # Extend the model past l = 2508 with a smooth exponential continuation so
    # the simulated sky has the damping tail the real one has.
    ell = np.arange(LMAX + 1)
    dl = np.zeros(LMAX + 1)
    dl[theory_ell] = theory_dl
    tail = (theory_ell > 2300)
    slope = np.polyfit(theory_ell[tail], np.log(theory_dl[tail]), 1)
    ext = ell > theory_ell.max()
    dl[ext] = np.exp(np.polyval(slope, ell[ext]))
    cl_in = cp.cl_from_dl(ell, dl)

    np.random.seed(RANDOM_SEED)
    alm = hp.synalm(cl_in, lmax=LMAX, new=True)
    beam = read_beam()
    sim = hp.alm2map(hp.almxfl(alm, beam), NSIDE, lmax=LMAX)
    cl_full = hp.alm2cl(alm)

    mask = binary_mask()
    apod = apodize(mask)
    w2 = cp.mask_moments(apod)["w2"]
    pseudo = hp.alm2cl(hp.map2alm(apod * sim, lmax=LMAX, iter=3))
    mask_cl = mask_power_spectrum(apod)
    corrected = cp.deconvolve_mask([pseudo], mask_cl, lmax=LMAX)[0]
    out = {"ell": ell, "cl_full": cl_full, "cl_input": cl_in,
           "cl_masked": _debeam(corrected, beam),
           "cl_masked_fsky": _debeam(pseudo / w2, beam)}
    np.savez(SIM_CACHE, **out)
    _log("simulation done", t0)
    return out


# --------------------------------------------------------------------------
# Summaries
# --------------------------------------------------------------------------


def bin_cl(ell, cl, *, lmin: int = 30, lmax: int = 2508):
    """Bin a C_l spectrum with Planck's scheme; returns (l_eff, D_b)."""
    return cp.planck_bin(ell[: lmax + 1], cp.dl_from_cl(ell, cl)[: lmax + 1], lmin=lmin, lmax=lmax)


def binned_errors(l_eff, cl_signal_b, cl_noise_b, *, fsky_eff: float, width: int = 30):
    """Knox errors on binned D_l for the cross-spectrum of two half-mission maps.

    ``cl_noise_b`` is the noise of the *sum* map, as measured from the
    difference map; each half map is twice as noisy.
    """
    n_half = 2.0 * cl_noise_b
    err_cl = cp.knox_error(l_eff, cl_signal_b, n_half, fsky=fsky_eff, delta_ell=width,
                           cross_nl=n_half)
    return cp.dl_from_cl(l_eff, err_cl)


def compare_to_planck(spec: dict, planck_binned, planck_bestfit_binned) -> dict:
    """How well does our cross-spectrum reproduce Planck's published one?"""
    ell = spec["ell"]
    l_eff, cross_b = bin_cl(ell, spec["cross"])
    _, auto_b = bin_cl(ell, spec["auto"])
    _, noise_b = bin_cl(ell, spec["noise"])
    model_cl_b = cp.cl_from_dl(l_eff, planck_bestfit_binned)
    err_b = binned_errors(l_eff, model_cl_b, cp.cl_from_dl(l_eff, noise_b),
                          fsky_eff=float(spec["fsky_eff"]))
    ratio = cross_b / planck_binned.dl
    out = {"l_eff": l_eff, "cross": cross_b, "auto": auto_b, "noise": noise_b,
           "err": err_b, "ratio_to_planck": ratio}
    for lo, hi in ((30, 500), (500, 1000), (1000, 1500), (1500, 2000), (2000, 2508)):
        sel = (l_eff >= lo) & (l_eff < hi)
        out[f"mean_ratio_{lo}_{hi}"] = float(np.mean(ratio[sel]))
        resid = (cross_b[sel] - planck_bestfit_binned[sel]) / err_b[sel]
        out[f"chi2_{lo}_{hi}"] = float(np.sum(resid**2))
        out[f"n_{lo}_{hi}"] = int(sel.sum())
    # Where does the noise overtake the signal in the auto-spectrum?
    above = np.nonzero(noise_b > planck_bestfit_binned)[0]
    out["noise_equals_signal_ell"] = float(l_eff[above[0]]) if len(above) else float("nan")

    # The excess at high l, tested against the shape of unresolved point sources:
    # randomly placed sources have a white (l-independent) C_l, i.e. D_l ~ l^2.
    # Amplitude quoted as D_l at l = 3000, the convention of Planck's likelihood.
    sel = l_eff >= POISSON_FIT_LMIN
    shape = (l_eff[sel] / 3000.0) ** 2
    excess = cross_b[sel] - planck_bestfit_binned[sel]
    w = 1.0 / err_b[sel] ** 2
    amp = float(np.sum(w * excess * shape) / np.sum(w * shape**2))
    out["poisson_amp"] = amp
    out["poisson_amp_err"] = float(1.0 / np.sqrt(np.sum(w * shape**2)))
    out["poisson_chi2"] = float(np.sum(w * (excess - amp * shape) ** 2))
    out["poisson_dof"] = int(sel.sum() - 1)
    return out



def predicted_rms(theory_ell, theory_dl, beam) -> float:
    """RMS temperature a sky with Planck's best-fit spectrum would have in this map."""
    ell = np.asarray(theory_ell)
    cl = cp.cl_from_dl(ell, theory_dl)
    return cp.map_rms_from_cl(ell, cl, beam[ell])


def delta_t_over_t(rms_uk: float) -> float:
    return rms_uk / (T_CMB_K * 1e6)
