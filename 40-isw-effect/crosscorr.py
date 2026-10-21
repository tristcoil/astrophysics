"""Galaxy-CMB cross-correlations: the ISW effect over the whole sky.

Stacking looks at a few hundred special places. A cross-correlation uses every
pixel: if dark energy decays the potentials, the CMB should be a little hotter
wherever there are more galaxies, on angular scales of a few degrees. The
amplitude of that correlation is predicted by linear theory (CAMB), once we
know how strongly the galaxies trace matter (their bias, measured from the
galaxies' own clustering).

Tracers
-------
* **NVSS**: 1.4 million radio sources over 80% of the sky. Radio galaxies are
  seen out to z ~ 3 (mean z ~ 1.2), right where dark energy starts to act; NVSS
  has been the single most powerful ISW tracer for 20 years. Its redshift
  distribution is not measured source by source; we use Planck 2015 XXI's.
* **BOSS**, in two redshift slices (LOWZ and CMASS): spectroscopic redshifts,
  so the redshift distribution is exact.
* **The BOSS ISW template**: the predicted ISW sky itself (boss.py), fitted to
  the CMB with an inverse-variance filter.

The cut sky
-----------
Multiplying a map by a mask mixes neighbouring multipoles. In angle space the
mixing is simply a product: the masked correlation function is the true one
times the mask's own correlation function (Szapudi et al. 2001). So we predict
exactly what our estimator (pseudo-C_l divided by the mean of the mask product)
measures for a given theory: transform the theory to a correlation function,
multiply by the mask's, transform back. Data and theory are compared in the
same "pseudo" space and no inversion is needed.
"""

from __future__ import annotations

import gzip
import time

import numpy as np

import boss
import cmb
import cosmology
from config import (BOSS_SLICES, CC_DL, CC_LMAX, CC_LMIN, CC_NSIDE, DERIVED_DIR, LITERATURE,
                    N_ROTATIONS, NVSS_BMIN, NVSS_DEC_BAND_DEG, NVSS_DEC_MIN, NVSS_SMIN_MJY, product)

NVSS_CACHE = DERIVED_DIR / "nvss.npz"
TRACER_CACHE = DERIVED_DIR / "cc_tracers.npz"
CC_SIMS = DERIVED_DIR / "sims_cc.npz"
LMAX = 3 * CC_NSIDE - 1
EDGES = np.arange(CC_LMIN, CC_LMAX + 1, CC_DL)  # 10, 20, ..., 100
TEMPLATE_TAGS = ("fiducial", "raw20", "smooth5", "smooth20", "nokaiser", "bias_lo", "bias_hi")


def _hp():
    import healpy as hp
    return hp


def _log(msg, t0=None):
    extra = f"  ({time.monotonic() - t0:.1f} s)" if t0 is not None else ""
    print(f"    {msg}{extra}", flush=True)


def bin_cl(cl):
    ell = np.arange(len(cl))
    return np.array([cl[(ell >= lo) & (ell < hi)].mean() for lo, hi in zip(EDGES[:-1], EDGES[1:])])


ELL_BIN = 0.5 * (EDGES[:-1] + EDGES[1:]) - 0.5


# --------------------------------------------------------------------------
# Tracer maps
# --------------------------------------------------------------------------


def nvss_catalog() -> dict:
    """RA, Dec (deg) and 1.4 GHz flux (mJy) of all 1,773,484 NVSS components."""
    if NVSS_CACHE.exists():
        d = np.load(NVSS_CACHE)
        return {k: d[k] for k in d.files}
    ra, dec, flux = [], [], []
    with gzip.open(product("nvss").path, "rt") as fh:
        for line in fh:
            if len(line) < 83:
                continue
            ra.append((int(line[40:42]) + int(line[43:45]) / 60 + float(line[46:51]) / 3600) * 15)
            sgn = -1.0 if line[52] == "-" else 1.0
            dec.append(sgn * (int(line[53:55]) + int(line[56:58]) / 60 + float(line[59:63]) / 3600))
            flux.append(float(line[75:83]))
    out = {"ra": np.array(ra), "dec": np.array(dec), "flux": np.array(flux)}
    np.savez(NVSS_CACHE, **out)
    return out


def _pixel_dec(nside):
    """Declination of every Galactic pixel centre."""
    hp = _hp()
    th, ph = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    r = hp.Rotator(coord=["G", "C"])
    th_c, _ = r(th, ph)
    return 90 - np.degrees(th_c)


def nvss_map():
    """NVSS overdensity map (Galactic, Nside 64) and its mask.

    Sources above 2.5 mJy; the survey limit Dec >= -40 deg (pixels must lie
    entirely north of it); |b| >= 10 deg. NVSS is known to have declination-
    dependent source densities (different VLA configurations and calibrations,
    Blake & Wall 2002), so the mean density is computed in 1-deg declination
    bands and each pixel compared with its own band.
    """
    hp = _hp()
    cat = nvss_catalog()
    sel = (cat["flux"] >= NVSS_SMIN_MJY) & (cat["dec"] >= NVSS_DEC_MIN)
    l, b = cmb.equatorial_to_galactic(cat["ra"][sel], cat["dec"][sel])
    npix = hp.nside2npix(CC_NSIDE)
    counts = np.bincount(hp.ang2pix(CC_NSIDE, l, b, lonlat=True), minlength=npix).astype(float)
    pdec = _pixel_dec(CC_NSIDE)
    _, pb = hp.pix2ang(CC_NSIDE, np.arange(npix), lonlat=True)
    margin = np.degrees(hp.max_pixrad(CC_NSIDE))
    mask = (pdec >= NVSS_DEC_MIN + margin) & (np.abs(pb) >= NVSS_BMIN) & (counts > 0)
    band = np.floor(pdec / NVSS_DEC_BAND_DEG).astype(int)
    delta = np.zeros(npix)
    for bb in np.unique(band[mask]):
        s = mask & (band == bb)
        delta[s] = counts[s] / counts[s].mean() - 1.0
    nbar_sr = counts[mask].mean() / hp.nside2pixarea(CC_NSIDE)
    return delta, mask, {"n_sources": int(counts[mask].sum()), "nbar_sr": float(nbar_sr),
                         "fsky": float(mask.mean())}


def boss_slice_maps():
    """BOSS overdensity maps (Galactic, Nside 64) in the two redshift slices.

    Expected counts come from the rasterised footprint: within each n(z) region,
    a pixel expects (its completeness) x (the region's galaxies) / (the region's
    total completeness). A Nside-64 pixel is used if at least half of it lies in
    the footprint.
    """
    hp = _hp()
    gal, foot = boss.galaxies(), boss.footprint()
    comp = foot["comp"].max(axis=0)
    grp = foot["group"]
    nb = boss.BOSS_MASK_NSIDE
    ra, dec = hp.pix2ang(nb, np.arange(hp.nside2npix(nb)), lonlat=True)
    l, b = cmb.equatorial_to_galactic(ra, dec)
    fine_to_coarse = hp.ang2pix(CC_NSIDE, l, b, lonlat=True)
    npix = hp.nside2npix(CC_NSIDE)
    cover = np.bincount(fine_to_coarse, (comp > 0).astype(float), npix) / (nb // CC_NSIDE) ** 2
    mask = cover >= 0.5
    gl, gb = cmb.equatorial_to_galactic(gal["ra"], gal["dec"])
    gpix = hp.ang2pix(CC_NSIDE, gl, gb, lonlat=True)
    maps, info = [], []
    for z1, z2 in BOSS_SLICES:
        s = (gal["z"] >= z1) & (gal["z"] < z2)
        expected_fine = np.zeros(len(comp))
        for ig in range(len(boss.GROUPS)):
            gsel = grp == ig
            n_g = gal["w"][s & (gal["group"] == ig)].sum()
            expected_fine[gsel] = comp[gsel] * n_g / comp[gsel].sum()
        expected = np.bincount(fine_to_coarse, expected_fine, npix)
        counts = np.bincount(gpix[s], gal["w"][s], npix)
        delta = np.zeros(npix)
        delta[mask] = counts[mask] / expected[mask] - 1.0
        kappa = float(np.sum(gal["w"][s] ** 2) / np.sum(gal["w"][s]))
        nbar_sr = float(counts[mask].sum() / (mask.sum() * hp.nside2pixarea(CC_NSIDE)) / kappa)
        hist, e = np.histogram(gal["z"][s], np.arange(z1, z2 + 1e-9, 0.01), weights=gal["w"][s])
        maps.append(delta)
        info.append({"z": 0.5 * (e[1:] + e[:-1]), "dndz": hist, "nbar_sr": nbar_sr,
                     "n_gal": int(s.sum()), "fsky": float(mask.mean())})
    return maps, mask, info


def tracers(*, force: bool = False) -> dict:
    """All tracer maps, masks and metadata (cached)."""
    if TRACER_CACHE.exists() and not force:
        d = np.load(TRACER_CACHE, allow_pickle=True)
        return d["tracers"].item()
    nv, nvmask, nvinfo = nvss_map()
    bmaps, bmask, binfo = boss_slice_maps()
    out = {"nvss": {"map": nv, "mask": nvmask, **nvinfo}}
    for (z1, z2), m, inf in zip(BOSS_SLICES, bmaps, binfo):
        out[f"boss_{z1:.2f}_{z2:.2f}"] = {"map": m, "mask": bmask, **inf}
    np.savez(TRACER_CACHE, tracers=np.array(out, dtype=object))
    return out


# --------------------------------------------------------------------------
# Spectra, forward model, theory
# --------------------------------------------------------------------------


def pseudo_cl(a, wa, b=None, wb=None):
    """Pseudo-spectrum of masked maps divided by the mean mask product."""
    hp = _hp()
    if b is None:
        b, wb = a, wa
    c = hp.anafast(a * wa, b * wb, lmax=LMAX, use_pixel_weights=False)
    return c / np.mean(wa * wb)


def forward(cl_theory, wa, wb, n_nodes: int = 600):
    """Exact expected pseudo-C_l (same normalisation as pseudo_cl) for a theory C_l.

    xi_masked(theta) = xi_true(theta) * xi_W(theta): multiply the correlation
    functions of the theory and of the two masks, project back to multipoles.
    """
    from numpy.polynomial import legendre as npl
    from scipy.special import roots_legendre

    hp = _hp()
    cw = hp.anafast(wa.astype(float), wb.astype(float), lmax=3 * CC_NSIDE - 1)
    mu, wts = roots_legendre(n_nodes)
    lw = np.arange(len(cw))
    xi_w = npl.legval(mu, (2 * lw + 1) / (4 * np.pi) * cw)
    lt = np.arange(len(cl_theory))
    xi_t = npl.legval(mu, (2 * lt + 1) / (4 * np.pi) * cl_theory)
    prod = xi_t * xi_w
    out = np.empty(LMAX + 1)
    pl_prev, pl = np.ones_like(mu), mu.copy()
    out[0] = 2 * np.pi * np.sum(wts * prod)
    out[1] = 2 * np.pi * np.sum(wts * prod * mu)
    for ell in range(2, LMAX + 1):
        pl_prev, pl = pl, ((2 * ell - 1) * mu * pl - (ell - 1) * pl_prev) / ell
        out[ell] = 2 * np.pi * np.sum(wts * prod * pl)
    return out / np.mean(wa * wb)


def nvss_window(bias_amp: float = 1.0):
    z = np.linspace(0.005, 6.0, 1200)
    p = LITERATURE["planck2015_xxi"]["nvss_dndz"]
    dndz = (z / p["z0"]) ** p["alpha"] * np.exp(-p["alpha"] * z / p["z0"])
    bias = bias_amp * 0.90 * (1 + 0.54 * (1 + z) ** 2)
    return {"z": z, "dndz": dndz, "bias": bias}


def boss_window(info: dict, b_eff: float, z_eff: float):
    L = cosmology.LCDM
    z = np.concatenate([[info["z"][0] - 0.01], info["z"], [info["z"][-1] + 0.01]])
    dndz = np.concatenate([[0.0], info["dndz"], [0.0]])
    bias = b_eff * L.D(z_eff) / L.D(z)
    return {"z": z, "dndz": dndz, "bias": bias}


def theory(window, *, lmax: int = LMAX, cosmo_overrides=None):
    """CAMB C_l^{gT} (uK) and C_l^{gg} for one tracer, times the Nside-64 pixel window^2."""
    hp = _hp()
    cls = cosmology.source_cls([window], lmax=lmax, cosmo_overrides=cosmo_overrides)
    pw = hp.pixwin(CC_NSIDE, lmax=lmax)[: lmax + 1] ** 2
    return {"gT": cls["TxW1"][: lmax + 1] * cmb.T_CMB_UK * pw,
            "gg": cls["W1xW1"][: lmax + 1] * pw,
            "TT_isw_free": None}


def fit_bias(name: str, tr: dict, *, b0: float):
    """Bias amplitude from the auto-spectrum over l = 10-100 (shot noise subtracted)."""
    hp = _hp()
    data = bin_cl(pseudo_cl(tr["map"], tr["mask"]) - 1.0 / tr["nbar_sr"])
    win = nvss_window(1.0) if name == "nvss" else None
    if win is None:
        zeff = float(np.sum(tr["z"] * tr["dndz"]) / np.sum(tr["dndz"]))
        win = boss_window(tr, b0, zeff)
    th = theory(win)
    model = bin_cl(forward(th["gg"], tr["mask"], tr["mask"]))
    # Gaussian errors for the auto-spectrum binned over Delta l = 10
    fsky = tr["mask"].mean()
    err = np.sqrt(2.0 / ((2 * ELL_BIN + 1) * CC_DL * fsky)) * (model + 1.0 / tr["nbar_sr"])
    # amplitude^2 multiplies the density term; RSD makes it slightly non-quadratic,
    # so iterate once with the updated bias.
    amp2 = float(np.sum(data * model / err**2) / np.sum(model**2 / err**2))
    return {"amp": float(np.sqrt(max(amp2, 0))), "data": data, "model_unit": model, "err": err}


# --------------------------------------------------------------------------
# Measurement, simulations, amplitude fits
# --------------------------------------------------------------------------


def models(*, force: bool = False) -> dict:
    """Bias fits and the forward-modelled pseudo-C_l^{gT} (A_ISW = 1) for every tracer."""
    path = DERIVED_DIR / "cc_models.npz"
    if path.exists() and not force:
        return np.load(path, allow_pickle=True)["models"].item()
    tr = tracers()
    b3d = boss.measure_bias()
    cmbmask = cmb.prepare_maps()["mask64"].astype(float)
    out = {}
    for name, t in tr.items():
        if name == "nvss":
            fb = fit_bias(name, t, b0=1.0)
            win = nvss_window(fb["amp"])
            b_label = f"{fb['amp']:.3f} x Planck 2015 XXI b(z)"
        else:
            b0 = b3d["lo"]["b"] if name.startswith("boss_0.15") else b3d["hi"]["b"]
            fb = fit_bias(name, t, b0=b0)
            zeff = float(np.sum(t["z"] * t["dndz"]) / np.sum(t["dndz"]))
            win = boss_window(t, fb["amp"] * b0, zeff)
            b_label = f"b = {fb['amp'] * b0:.3f} at z = {zeff:.2f} (3D P(k): {b0:.3f})"
        th = theory(win)
        gt = bin_cl(forward(th["gT"], t["mask"].astype(float), cmbmask))
        gg = bin_cl(forward(th["gg"], t["mask"].astype(float), t["mask"].astype(float)))
        out[name] = {"gT_model": gt, "gg_model": gg, "gg_data": fb["data"], "gg_err": fb["err"],
                     "bias_amp": fb["amp"], "bias_label": b_label, "cl_gT_full": th["gT"]}
        _log(f"{name}: {b_label}")
    np.savez(path, models=np.array(out, dtype=object))
    return out


def template_filters(tags=TEMPLATE_TAGS):
    """For each BOSS ISW template t: u = W F t and norm = <u, t>, so A(d) = <u, d> / norm.

    F keeps l = 10-100 and weights each multipole by 1 / C_l^TT (the CMB is the
    noise); W is the CMB mask. Any fixed filter gives an unbiased A; this one is
    close to optimal.
    """
    import healpy as hp

    W = cmb.prepare_maps()["mask64"].astype(float)
    cl = cmb.theory_cl(LMAX) * cmb.beam(LMAX, CC_NSIDE) ** 2
    fl = np.zeros(LMAX + 1)
    sel = (np.arange(LMAX + 1) >= CC_LMIN) & (np.arange(LMAX + 1) <= CC_LMAX)
    fl[sel] = 1.0 / cl[sel]
    out = {}
    for tag in tags:
        m = boss.predicted_isw(tag)
        t64 = hp.ud_grade(m, CC_NSIDE)
        ft = hp.alm2map(hp.almxfl(hp.map2alm(t64 * W, lmax=LMAX), fl), CC_NSIDE, lmax=LMAX)
        u = W * ft
        out[tag] = (u, float(np.sum(u * t64)), t64)
    return out


def measure(d64, mask64, tr, tfilters):
    """Binned cross-spectra with every tracer and template amplitudes for one CMB map."""
    import healpy as hp

    alm_t = hp.map2alm(d64 * mask64, lmax=LMAX)
    res = {}
    for name, t in tr.items():
        cl = hp.alm2cl(alm_t, t["alm"]) / np.mean(mask64 * t["mask"])
        res[name] = bin_cl(cl)
    for tag, (u, norm, _) in tfilters.items():
        res[f"template_{tag}"] = float(np.sum(u * d64) / norm)
    return res


def _prepared_tracers():
    import healpy as hp

    tr = tracers()
    for t in tr.values():
        t["alm"] = hp.map2alm(t["map"] * t["mask"], lmax=LMAX)
    return tr


def run_cc_sims(n_sims: int, *, force: bool = False) -> dict:
    """Null distributions: simulated CMB skies against the real galaxy maps."""
    if CC_SIMS.exists() and not force:
        d = np.load(CC_SIMS, allow_pickle=True)
        return d["sims"].item()
    import stacking as S

    d = cmb.prepare_maps()
    mask = d["mask"]
    tr = _prepared_tracers()
    tf = template_filters()
    cl = cmb.theory_cl()
    acc: dict = {}
    t0 = time.monotonic()
    for i in range(n_sims):
        m512 = cmb.alm_to_map(cmb.sim_alm(i, cl), 512)
        m64, mask64 = cmb.degrade(np.where(mask, m512, 0), mask, CC_NSIDE, cmb.LOWRES_KEEP)
        for k, v in measure(m64, mask64.astype(float), tr, tf).items():
            acc.setdefault(k, []).append(v)
        if (i + 1) % 100 == 0:
            _log(f"cc sim {i + 1}/{n_sims}", t0)
    out = {k: np.asarray(v) for k, v in acc.items()}
    np.savez(CC_SIMS, sims=np.array(out, dtype=object))
    return out


def rotated_measurements() -> dict:
    """The real CMB rotated in Galactic longitude by 10, 20, ..., 350 deg."""
    import healpy as hp

    d = cmb.prepare_maps()
    tr = _prepared_tracers()
    tf = template_filters()
    acc: dict = {}
    for k in range(1, N_ROTATIONS + 1):
        rot = hp.Rotator(rot=(10.0 * k, 0.0, 0.0), coord="G")
        m = rot.rotate_map_pixel(d["smica64"].astype(float))
        mk = rot.rotate_map_pixel(d["mask64"].astype(float)) > 0.5
        for key, v in measure(m * mk, mk.astype(float), tr, tf).items():
            acc.setdefault(key, []).append(v)
    return {k: np.asarray(v) for k, v in acc.items()}


def fit_amplitude(data, model, sims):
    """A = (d C^-1 t) / (t C^-1 t) with C from the simulations (Hartlap-corrected)."""
    n, p = sims.shape
    cov = np.cov(sims, rowvar=False)
    icov = np.linalg.inv(cov) * (n - p - 2) / (n - 1)
    norm = model @ icov @ model
    a = float(model @ icov @ data / norm)
    sigma = float(1 / np.sqrt(norm))
    a_sims = (sims - sims.mean(0)) @ icov @ model / norm
    chi2 = float((data - a * model) @ icov @ (data - a * model))
    return {"A": a, "sigma": sigma, "sigma_sims": float(np.std(a_sims, ddof=1)), "snr": a / sigma,
            "chi2": chi2, "dof": int(p - 1), "A_sims": a_sims}
