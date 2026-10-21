"""Every THINGS measurement in the project, and its grading against published work.

For each galaxy: velocity fields from four recipes (plus the survey's own
moment-1 map), our tilted-ring fit to each, the dynamical centre and systemic
velocity, the rotation curve with approaching and receding halves, the
residual velocity field and its harmonic decomposition, and the gas.
For DDO 154 additionally the cusp-core test: density profiles, halo fits, and
model galaxies run through the identical pipeline. The MaNGA half lives in
``manga.py``.

Heavy products are cached in ``data/derived/`` (pickles), so figures can be
redrawn without refitting. ``run.py`` calls everything and prints the summary.
"""

from __future__ import annotations

import pickle
import time
from dataclasses import dataclass, field

import numpy as np

import galmod as GM
import harmonics as H
import literature as L
import massmodel as MM
import things
import tiltedring as TR
from config import (DERIVED_DIR, GALAXIES, GAS_PER_HI, MODEL_SCALE_HEIGHT_KPC, RINGS_PER_BEAM,
                    Galaxy, galaxy_by_name)

#: The velocity-field recipes compared (first = the one every result uses).
VELOCITY_KINDS = {
    "gh_vpeak": "Gauss-Hermite h3, profile peak",
    "gh_v0": "Gauss-Hermite h3, centre parameter",
    "gauss_v": "single Gaussian",
    "peak_v": "brightest channel",
    "mom1_things": "first moment (THINGS map)",
}
PRIMARY = "gh_vpeak"


# --------------------------------------------------------------------------
# Coordinates
# --------------------------------------------------------------------------

@dataclass
class Frame:
    """Arcsecond offsets for a galaxy's maps: x west (with the columns), y north."""

    pixel: float
    cref: float
    rref: float
    shape: tuple[int, int]
    wcs: object

    def xy(self, rows, cols):
        return (np.asarray(cols) - self.cref) * self.pixel, (np.asarray(rows) - self.rref) * self.pixel

    def grid(self):
        ny, nx = self.shape
        return np.meshgrid((np.arange(nx) - self.cref) * self.pixel, (np.arange(ny) - self.rref) * self.pixel)

    def radec(self, x, y):
        ra, dec = self.wcs.pixel_to_world_values(x / self.pixel + self.cref, y / self.pixel + self.rref)
        return float(ra), float(dec)

    def from_radec(self, ra, dec):
        col, row = self.wcs.world_to_pixel_values(ra, dec)
        return float((col - self.cref) * self.pixel), float((row - self.rref) * self.pixel)


def frame_of(cube: things.Cube) -> Frame:
    h = cube.header
    return Frame(cube.pixel_arcsec, h["CRPIX1"] - 1, h["CRPIX2"] - 1, (h["NAXIS2"], h["NAXIS1"]), cube.wcs)


def ring_step(cube: things.Cube) -> float:
    """Half the beam, rounded to 0.5 arcsec: two rings per resolution element."""
    return round(cube.beam_fwhm_arcsec / RINGS_PER_BEAM * 2) / 2


def offset_arcsec(ra1, dec1, ra2, dec2) -> float:
    return float(np.hypot((ra1 - ra2) * np.cos(np.radians(dec2)), dec1 - dec2) * 3600)


# --------------------------------------------------------------------------
# One galaxy
# --------------------------------------------------------------------------

@dataclass
class GalaxyResult:
    name: str
    key: str
    distance_mpc: float
    kpc_per_arcsec: float
    beam: tuple
    noise_mjy: float
    channel_kms: float
    pixel_arcsec: float
    step_arcsec: float
    n_fitted: int
    n_good: dict
    tr: dict  # kind -> TiltedRingResult
    centre_radec: tuple
    deblok: L.PublishedCurve
    sparc: L.PublishedCurve
    comparisons: dict
    harmonics: H.Harmonics  # around the adopted (smoothed) rings: the residual maps' view
    harmonic_summary: dict
    gas: dict
    flux: dict
    extra: dict = field(default_factory=dict)

    @property
    def primary(self) -> TR.TiltedRingResult:
        return self.tr[PRIMARY]

    def r_kpc(self, kind=PRIMARY) -> np.ndarray:
        return self.tr[kind].radius * self.kpc_per_arcsec


def compare_curves(r_ours, v_ours, e_ours, pub: L.PublishedCurve, inner_kpc=None) -> dict:
    """Our curve interpolated at a published curve's radii, where both exist."""
    ok = np.isfinite(v_ours)
    vi = np.interp(pub.r_kpc, r_ours[ok], v_ours[ok], left=np.nan, right=np.nan)
    ei = np.interp(pub.r_kpc, r_ours[ok], e_ours[ok], left=np.nan, right=np.nan)
    d = vi - pub.v
    m = np.isfinite(d)
    out = {"n": int(m.sum()), "median_diff": float(np.median(d[m])), "rms_diff": float(np.sqrt(np.mean(d[m] ** 2))),
           "mean_abs_frac": float(np.mean(np.abs(d[m]) / pub.v[m])),
           "within_1sigma": float(np.mean(np.abs(d[m]) <= np.hypot(ei[m], pub.err[m]))),
           "r_max_kpc": float(pub.r_kpc[m].max())}
    if inner_kpc:
        inner = m & (pub.r_kpc <= inner_kpc)
        out["inner_kpc"] = inner_kpc
        out["inner_median_diff"] = float(np.median(d[inner])) if inner.any() else None
    return out


def profile_midpoint(cube: things.Cube, mask: np.ndarray, level: float = 0.5) -> float:
    """Midpoint of the global HI profile at ``level`` of its peak (W50 by default).

    An independent check on the fitted systemic velocity. Scanning inward from both ends of the spectrum to the first channel above
    the level, with linear interpolation between channels.
    """
    spec = _global_spectrum(cube, mask)
    v = cube.velocity
    thr = level * spec.max()
    i = int(np.argmax(spec > thr))
    j = len(spec) - 1 - int(np.argmax(spec[::-1] > thr))
    lo = v[i - 1] + (thr - spec[i - 1]) * (v[i] - v[i - 1]) / (spec[i] - spec[i - 1]) if i > 0 else v[i]
    hi = v[j] + (thr - spec[j]) * (v[j + 1] - v[j]) / (spec[j + 1] - spec[j]) if j < len(v) - 1 else v[j]
    return float(0.5 * (lo + hi))


def _pixels(galaxy: Galaxy, kind: str, maps, cube, frame: Frame):
    vf = things.velocity_field(galaxy, kind, maps, cube)
    rows, cols = np.nonzero(np.isfinite(vf))
    x, y = frame.xy(rows, cols)
    return x, y, vf[rows, cols], maps["mom0_things"][rows, cols]


def _fit(galaxy: Galaxy, kind: str, maps, cube, frame: Frame, step: float) -> TR.TiltedRingResult:
    x, y, v, w = _pixels(galaxy, kind, maps, cube, frame)
    return TR.analyse(x, y, v, w, ring_step=step,
                      beam_fwhm=cube.beam_fwhm_arcsec, beam_area_arcsec2=cube.beam_area_pixels * cube.pixel_arcsec ** 2,
                      inner_fixed_arcsec=galaxy.geometry_fixed_inside_arcsec)


def save_maps(galaxy: Galaxy, res: TR.TiltedRingResult, maps, cube, frame: Frame) -> None:
    """Velocity, model and residual maps of the primary fit, cached for the figures."""
    vf = things.velocity_field(galaxy, PRIMARY, maps, cube)
    rows, cols = np.nonzero(np.isfinite(vf))
    x, y = frame.xy(rows, cols)
    model, _ = TR.model_field(x, y, res)
    resid = np.full(vf.shape, np.nan, dtype=np.float32)
    resid[rows, cols] = vf[rows, cols] - model
    model_map = np.full(vf.shape, np.nan, dtype=np.float32)
    model_map[rows, cols] = model
    np.savez_compressed(DERIVED_DIR / f"{galaxy.key}_maps.npz", velocity=vf.astype(np.float32),
                        model=model_map, residual=resid, mom0=maps["mom0_things"].astype(np.float32),
                        mom1=maps["mom1_things"].astype(np.float32),
                        gh_sigma=np.where(np.isfinite(vf), maps["gh_sigma"], np.nan).astype(np.float32),
                        gh_h3=np.where(np.isfinite(vf), maps["gh_h3"], np.nan).astype(np.float32))


def analyse_galaxy(galaxy: Galaxy, *, quiet=False) -> GalaxyResult:
    t0 = time.time()
    cube = things.load_cube(galaxy)
    maps = things.velocity_fields(galaxy, quiet=quiet)
    frame = frame_of(cube)
    d_mpc = L.walter_distance_mpc(galaxy.key)
    kpa = d_mpc * 1e3 * np.pi / 180 / 3600
    step = ring_step(cube)

    vsys_w50 = profile_midpoint(cube, maps["mom0_things"] > 0)
    tr, n_good = {}, {}
    for kind in VELOCITY_KINDS:
        tr[kind] = _fit(galaxy, kind, maps, cube, frame, step)
        n_good[kind] = int(np.isfinite(things.velocity_field(galaxy, kind, maps, cube)).sum())
    res = tr[PRIMARY]

    # centre and geometry against de Blok et al.'s Table 3
    ra, dec = frame.radec(*res.centre)
    geo = L.deblok_geometry()[galaxy.key]
    deblok = L.deblok_rotation_curve(galaxy.key)
    sparc = L.sparc_rotation_curve(galaxy.sparc, d_mpc)
    inner = 0.15 * deblok.r_kpc.max()
    comparisons = {"geometry": {"centre_offset_arcsec": offset_arcsec(ra, dec, geo["ra_deg"], geo["dec_deg"]),
                                "vsys": res.vsys, "vsys_deblok": geo["vsys"],
                                "incl_mean": float(np.nanmean(res.incl)), "incl_deblok": geo["incl"],
                                "pa_mean": float(np.nanmean(res.pa)), "pa_deblok": geo["pa"],
                                "step_arcsec": step, "step_deblok": geo["ring_step_arcsec"],
                                "vsys_global_w50": vsys_w50,
                                },
                   "deblok": {}, "sparc": {}}
    for kind, t in tr.items():
        r = t.radius * kpa
        comparisons["deblok"][kind] = compare_curves(r, t.vrot, t.vrot_err, deblok, inner_kpc=inner)
        comparisons["deblok"][kind]["vsys"] = t.vsys
    comparisons["sparc"] = compare_curves(res.radius * kpa, res.vrot, res.vrot_err, sparc)
    comparisons["sparc"]["incl_sparc"] = L.sparc_galaxy_table()[galaxy.sparc]["incl"]
    comparisons["sparc"]["distance_sparc"] = L.sparc_galaxy_table()[galaxy.sparc]["distance_mpc"]
    comparisons["sparc"]["refs"] = L.sparc_galaxy_table()[galaxy.sparc]["refs"]

    # harmonic decomposition on the primary velocity field
    vf = things.velocity_field(galaxy, PRIMARY, maps, cube)
    rows, cols = np.nonzero(np.isfinite(vf))
    x, y = frame.xy(rows, cols)
    harm = H.decompose(x, y, vf[rows, cols], res)
    # like-for-like with Trachternach et al.: rings refitted without weights, as RESWRI does
    harm_reswri = H.decompose(x, y, vf[rows, cols], res, geometry="reswri")
    trach = L.trachternach_harmonics().get(galaxy.key, {})
    hsum = H.summary(harm, kpa)
    hsum_trach = H.summary(harm_reswri, kpa, r_max_arcsec=trach.get("r_max_arcsec"))

    save_maps(galaxy, res, maps, cube, frame)

    # gas: face-on surface density around the rings, and its own rotation curve
    X, Y = frame.grid()
    mom0 = maps["mom0_things"]
    ring_mom0 = MM.ring_profile(mom0.ravel(), X.ravel(), Y.ravel(), res)
    sigma_hi = MM.hi_surface_density(ring_mom0, cube.beam) * np.cos(np.radians(res.incl))
    r_kpc = res.radius * kpa
    good = np.isfinite(sigma_hi)
    sigma_gas = GAS_PER_HI * sigma_hi
    v2_gas = MM.disc_rotation_curve(r_kpc, r_kpc[good], sigma_gas[good], z0_kpc=MODEL_SCALE_HEIGHT_KPC)
    total_flux = float(np.sum(mom0[mom0 > 0]) / cube.beam_area_pixels)  # Jy km/s
    tv, tf = L.walter_global_spectrum(galaxy.key)
    walter_flux = float(np.sum(tf) * np.median(np.abs(np.diff(tv))))
    sparc_mhi = float(_sparc_mhi(galaxy.sparc) * (d_mpc / L.sparc_galaxy_table()[galaxy.sparc]["distance_mpc"]) ** 2)
    # our own global spectrum: every pixel with emission in the THINGS mask, per channel
    spec = _global_spectrum(cube, mom0 > 0)
    gas = {"r_kpc": r_kpc, "sigma_hi": sigma_hi, "sigma_gas": sigma_gas,
           "v2_gas": v2_gas, "m_hi": MM.hi_mass_msun(total_flux, d_mpc)}
    flux = {"mom0_flux_jykms": total_flux, "walter_flux_jykms": walter_flux,
            "m_hi_msun": MM.hi_mass_msun(total_flux, d_mpc), "m_hi_walter_msun": MM.hi_mass_msun(walter_flux, d_mpc),
            "m_hi_sparc_msun": sparc_mhi * 1e9, "spectrum_v": cube.velocity, "spectrum_ours": spec,
            "spectrum_walter_v": tv, "spectrum_walter": tf}

    if not quiet:
        c = comparisons["deblok"][PRIMARY]
        print(f"    {galaxy.name}: {len(res.radius)} rings to {r_kpc.max():.1f} kpc; vs de Blok "
              f"{c['median_diff']:+.2f} km/s median, {c['rms_diff']:.2f} rms ({time.time() - t0:.0f} s)", flush=True)
    return GalaxyResult(galaxy.name, galaxy.key, d_mpc, kpa, cube.beam, cube.noise * 1e3, cube.channel_width,
                        cube.pixel_arcsec, step, int((maps["mom0_things"] > 0).sum()), n_good, tr, (ra, dec),
                        deblok, sparc, comparisons, harm,
                        {"adopted_rings": hsum, "reswri_rings": hsum_trach, "trachternach": trach},
                        gas, flux, extra={"harmonics_reswri": harm_reswri})


def _sparc_mhi(name: str) -> float:
    """SPARC HI mass (1e9 Msun) at SPARC's distance."""
    for ln in (L.REF_DIR / "SPARC_Lelli2016c.mrt").read_text(encoding="latin-1").splitlines():
        parts = ln.split()
        if parts and parts[0] == name and len(parts) >= 19:
            return float(parts[13])
    return float("nan")


def _global_spectrum(cube: things.Cube, mask: np.ndarray) -> np.ndarray:
    """Flux density (Jy) per channel, summed over the masked pixels."""
    rows = np.nonzero(mask.any(axis=1))[0]
    out = np.zeros(len(cube.velocity))
    for y0 in range(rows.min(), rows.max() + 1, 64):
        slab = np.asarray(cube.data[:, y0:y0 + 64, :], dtype=np.float64)
        m = mask[y0:y0 + 64]
        out += np.nansum(slab[:, m], axis=1)
    return out / cube.beam_area_pixels


def galaxy_results(*, rebuild=False, quiet=False) -> dict[str, GalaxyResult]:
    path = DERIVED_DIR / "galaxy_results.pkl"
    if path.exists() and not rebuild:
        with open(path, "rb") as fh:
            return pickle.load(fh)
    out = {g.key: analyse_galaxy(g, quiet=quiet) for g in GALAXIES}
    with open(path, "wb") as fh:
        pickle.dump(out, fh)
    return out


# --------------------------------------------------------------------------
# DDO 154: the cusp-core test
# --------------------------------------------------------------------------

#: Oh et al. (2011) fitted DDO 154's inner slope inside 0.68 kpc (their five
#: innermost points, read from their figure); we use the same break radius.
CUSP_BREAK_KPC = 0.68
DENSITY_WINDOW = 5


def cusp_test(gr: GalaxyResult, *, rebuild=False, quiet=False) -> dict:
    """Density slopes, halo fits and the beam-smearing simulation for DDO 154."""
    path = DERIVED_DIR / "cusp_test.pkl"
    if path.exists() and not rebuild:
        with open(path, "rb") as fh:
            return pickle.load(fh)
    t0 = time.time()
    res = gr.primary
    r, v, e = gr.r_kpc(), res.vrot, res.vrot_err
    ok = np.isfinite(v)
    r, v, e = r[ok], v[ok], e[ok]
    db = gr.deblok
    out: dict = {"break_kpc": CUSP_BREAK_KPC, "window": DENSITY_WINDOW}

    # density profiles and inner slopes, ours and de Blok's curve through the same code
    for label, rr, vv in (("ours", r, v), ("deblok", db.r_kpc, db.v)):
        rho = MM.density_smoothed(rr, vv, DENSITY_WINDOW)
        out[f"rho_{label}"] = (rr, rho)
        out[f"alpha_{label}"] = MM.slope_within(rr, rho, CUSP_BREAK_KPC)
        out[f"vslope_{label}"] = MM.velocity_power_slope(rr, vv, 0.7)
        out[f"vslope_{label}_nofirst"] = MM.velocity_power_slope(rr, vv, 0.7, r_min_kpc=0.2)
    # baryons: subtract the gas ("minimum disc + gas")
    v2g = np.interp(r, gr.gas["r_kpc"], gr.gas["v2_gas"])
    v_dm = np.sqrt(np.clip(v ** 2 - v2g, 0, None))
    rho_dm = MM.density_smoothed(r, v_dm, DENSITY_WINDOW)
    out["rho_ours_dm"] = (r, rho_dm)
    out["alpha_ours_dm"] = MM.slope_within(r, rho_dm, CUSP_BREAK_KPC)

    # halo fits, minimum disc (as Oh et al.) and minimum disc + gas
    out["halo"] = {}
    for label, rr, vv, ee, vb2 in (("ours", r, v, e, None), ("deblok", db.r_kpc, db.v, db.err, None),
                                   ("ours_gas", r, v, e, v2g)):
        out["halo"][label] = {k: MM.fit_halo(k, rr, vv, ee, vb2) for k in ("ISO", "NFW")}
    lcdm = MM.fit_lcdm_nfw(r, v, e, 3.0)
    out["lcdm_nfw"] = lcdm

    # beam smearing: model galaxies through the identical pipeline (labelled simulation)
    galaxy = galaxy_by_name(gr.name)
    cube = things.load_cube(galaxy)
    maps = things.velocity_fields(galaxy)
    frame = frame_of(cube)
    m0 = maps["mom0_things"]
    X, Y = frame.grid()
    prof = MM.ring_profile(m0.ravel(), X.ravel(), Y.ravel(), res)
    incl, pa = float(np.nanmedian(res.incl)), float(np.nanmedian(res.pa))
    s0 = np.nan_to_num(prof * np.cos(np.radians(incl)))
    q = things.quality_mask(maps, cube)
    disp = float(np.nanmedian(maps["gh_sigma"][q]))
    rows, cols = np.nonzero(m0 > 0)
    mar = int(3 * cube.beam_fwhm_arcsec / cube.pixel_arcsec)
    r0, r1 = max(rows.min() - mar, 0), min(rows.max() + mar, frame.shape[0])
    c0, c1 = max(cols.min() - mar, 0), min(cols.max() + mar, frame.shape[1])
    xcol = (np.arange(frame.shape[1]) - frame.cref) * frame.pixel
    yrow = (np.arange(frame.shape[0]) - frame.rref) * frame.pixel
    vin_obs = np.where(np.isfinite(res.vrot), res.vrot, np.interp(res.radius, res.radius[np.isfinite(res.vrot)],
                                                                  res.vrot[np.isfinite(res.vrot)]))
    models = {"observed curve": vin_obs,
              "LCDM cusp (NFW)": MM.v_nfw(gr.r_kpc(), lcdm["c"], lcdm["v200"])}
    rng = np.random.default_rng(154)
    out["models"] = {}
    for name, vin in models.items():
        spec = GM.ModelSpec(res.centre, res.vsys, pa, incl, res.radius, vin, s0, disp,
                            MODEL_SCALE_HEIGHT_KPC / gr.kpc_per_arcsec)
        model = GM.make_cube(spec, (r1 - r0, c1 - c0), xcol[c0:c1], yrow[r0:r1], cube.velocity)
        clean = GM.observe(model, cube.beam, cube.pixel_arcsec, 0.0, rng)
        obs = GM.observe(model, cube.beam, cube.pixel_arcsec, cube.noise, rng)
        rec = GM.recover(obs, cube.velocity, cube.noise, xcol[c0:c1], yrow[r0:r1], spec, res.width,
                         mom0_noiseless=clean.sum(axis=0) * cube.channel_width)
        rk = rec.radius * gr.kpc_per_arcsec
        rho_rec = MM.density_smoothed(rk, rec.v_recovered, DENSITY_WINDOW)
        rho_in = MM.density_smoothed(rk, rec.v_input, DENSITY_WINDOW)
        out["models"][name] = {"r_kpc": rk, "v_input": rec.v_input, "v_recovered": rec.v_recovered,
                               "alpha_input": MM.slope_within(rk, rho_in, CUSP_BREAK_KPC),
                               "alpha_recovered": MM.slope_within(rk, rho_rec, CUSP_BREAK_KPC),
                               "vslope_recovered": MM.velocity_power_slope(rk, rec.v_recovered, 0.7),
                               "vslope_recovered_nofirst": MM.velocity_power_slope(rk, rec.v_recovered, 0.7, 0.2),
                               "dispersion": disp, "z0_kpc": MODEL_SCALE_HEIGHT_KPC, "n_pixels": rec.n_pixels}
        if not quiet:
            print(f"    model '{name}' observed and refitted ({time.time() - t0:.0f} s)", flush=True)
    # how far the smeared cusp sits above the real galaxy inside 1 kpc
    cusp = out["models"]["LCDM cusp (NFW)"]
    inner = (r <= 1.0)
    v_cusp = np.interp(r[inner], cusp["r_kpc"], cusp["v_recovered"])
    out["cusp_excess"] = {"r_kpc": r[inner], "excess_kms": v_cusp - v[inner], "err_kms": e[inner],
                          "mean_excess_kms": float(np.mean(v_cusp - v[inner])),
                          "mean_ratio": float(np.mean(v_cusp / v[inner])),
                          "chi2": float(np.sum(((v_cusp - v[inner]) / e[inner]) ** 2)),
                          "n": int(inner.sum())}
    with open(path, "wb") as fh:
        pickle.dump(out, fh)
    return out
