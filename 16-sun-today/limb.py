"""What the edge of the Sun says about the layers it is made of.

* **Visible light: limb darkening.** The photosphere is opaque: we see down to
  roughly one optical depth along our line of sight. At disc centre that depth lies
  deep and hot; near the limb the line of sight runs obliquely, reaches the same
  optical depth higher up, where it is cooler, and the disc looks darker. The shape of
  I(mu) therefore maps temperature against depth (Eddington-Barbier).
* **Extreme UV: limb brightening.** The corona is transparent to its own EUV lines,
  so an image adds up all the emission along the line of sight. Near the limb that
  path through the corona is longest, and the disc edge looks brighter.
* **Above the limb: the corona's scale height.** In hydrostatic balance, hotter gas
  is supported to greater heights (H = kT / mu m_H g). How fast each channel's
  emission fades above the limb measures the temperature of the gas it sees.
"""

from __future__ import annotations

import warnings

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord

import config as C
import sdo
import solar_physics as sp

warnings.filterwarnings("ignore")


def fast_radius(m) -> np.ndarray:
    """r / R_sun for every pixel from the disc centre's pixel position and the plate
    scale. At these small angles the TAN projection is linear to 1e-5."""
    centre = m.world_to_pixel(SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=m.coordinate_frame))
    ny, nx = m.data.shape
    y, x = np.ogrid[0:ny, 0:nx]
    scale = m.scale.axis1.to_value(u.arcsec / u.pix)
    rho = np.hypot(x - centre.x.value, y - centre.y.value) * scale
    return rho / m.rsun_obs.to_value(u.arcsec)


# --- Visible light ---------------------------------------------------------------------------


def continuum_limb_darkening(path) -> dict:
    """Measure I(mu)/I(1) on an HMI continuum image and interpret it.

    The median in each mu bin ignores sunspots (too few pixels to move a median);
    I(1) is the median within mu > 0.995 (the central 10% of the radius).
    """
    m = sdo.load_map(path)
    data = np.asarray(m.data, dtype=float)
    r = fast_radius(m)
    mu = sdo.mu_from_radius(r)
    on = (r < 1.0) & np.isfinite(data)
    centre = float(np.nanmedian(data[on & (mu > 0.995)]))

    edges = np.linspace(0.0, 1.0, C.LIMB_N_BINS + 1)
    mu_c, ratio, n = [], [], []
    d, mm = data[on], mu[on]
    order = np.argsort(mm)
    d, mm = d[order], mm[order]
    cuts = np.searchsorted(mm, edges)
    for k in range(C.LIMB_N_BINS):
        seg = slice(cuts[k], cuts[k + 1])
        if cuts[k + 1] - cuts[k] < 50:
            continue
        mu_c.append(float(np.median(mm[seg])))
        ratio.append(float(np.median(d[seg])) / centre)
        n.append(int(cuts[k + 1] - cuts[k]))
    mu_c, ratio, n = np.array(mu_c), np.array(ratio), np.array(n)

    # Fit Allen's form 1 - u - v + u mu + v mu^2 over the trusted range
    fit = mu_c >= C.LIMB_MU_MIN
    x = np.column_stack([mu_c[fit] - 1.0, mu_c[fit] ** 2 - 1.0])
    (u_fit, v_fit), *_ = np.linalg.lstsq(x, ratio[fit] - 1.0, rcond=None)
    model_fit = 1 + u_fit * (mu_c - 1) + v_fit * (mu_c**2 - 1)
    u_allen, v_allen = sp.allen_limb_coefficients(C.HMI_WAVELENGTH_A)
    allen = sp.allen_limb_darkening(mu_c, C.HMI_WAVELENGTH_A)
    eddington = sp.eddington_limb_darkening(mu_c)
    resid_allen = ratio[fit] / allen[fit] - 1

    # Temperature against optical depth, tau = mu, anchored at T(2/3) = T_eff
    wl = C.HMI_WAVELENGTH_A * 1e-10
    ratio_23 = float(1 + u_fit * (2 / 3 - 1) + v_fit * ((2 / 3) ** 2 - 1))
    t_tau = sp.eddington_barbier_temperature(mu_c[fit], ratio[fit], wl, ratio_at_anchor=ratio_23)
    t_centre = float(sp.eddington_barbier_temperature(1.0, 1.0, wl, ratio_at_anchor=ratio_23))
    brightness_t_limb = float(sp.eddington_barbier_temperature(C.LIMB_MU_MIN,
                                                               np.interp(C.LIMB_MU_MIN, mu_c, ratio),
                                                               wl, ratio_at_anchor=ratio_23))

    # Scattered light: the median just outside the limb, as a fraction of centre
    off = (r > 1.01) & (r < 1.03) & np.isfinite(data)
    scattered = float(np.nanmedian(data[off]) / centre) if off.any() else float("nan")

    return {
        "date": m.date.isot,
        "wavelength_A": C.HMI_WAVELENGTH_A,
        "mu": mu_c.tolist(),
        "ratio": ratio.tolist(),
        "n_pixels": n.tolist(),
        "fit_mu_min": C.LIMB_MU_MIN,
        "u_fit": float(u_fit),
        "v_fit": float(v_fit),
        "u_allen": u_allen,
        "v_allen": v_allen,
        "limb_ratio_fit": float(1 - u_fit - v_fit),
        "limb_ratio_allen": float(1 - u_allen - v_allen),
        "limb_ratio_eddington": 0.4,
        "max_abs_resid_allen": float(np.max(np.abs(resid_allen))),
        "rms_resid_allen": float(np.sqrt(np.mean(resid_allen**2))),
        "rms_resid_eddington": float(np.sqrt(np.mean((ratio[fit] / eddington[fit] - 1) ** 2))),
        "ratio_at_mu_0p1": float(np.interp(0.1, mu_c, ratio)),
        "allen_at_mu_0p1": float(sp.allen_limb_darkening(0.1, C.HMI_WAVELENGTH_A)),
        "model_fit": model_fit.tolist(),
        "tau": mu_c[fit].tolist(),
        "t_tau": t_tau.tolist(),
        "t_tau_grey": sp.grey_temperature(mu_c[fit]).tolist(),
        "t_at_tau_1": t_centre,
        "t_at_tau_0p1": brightness_t_limb,
        "scattered_light_fraction": scattered,
    }


# --- Extreme ultraviolet ----------------------------------------------------------------------


def radial_profiles(maps: dict, r_edges=None) -> dict:
    """Median intensity against r for each channel, normalised to disc centre
    (r < 0.3), from disc centre out into the corona."""
    if r_edges is None:
        r_edges = np.concatenate([np.linspace(0, 0.9, 19), np.linspace(0.9, 1.0, 21)[1:],
                                  np.linspace(1.0, 1.5, 51)[1:]])
    out = {"r": (0.5 * (r_edges[1:] + r_edges[:-1])).tolist()}
    for wl, m in maps.items():
        d = np.asarray(m.data, dtype=float)
        r = fast_radius(m)
        centre = float(np.nanmedian(d[r < 0.3]))
        prof = []
        for lo, hi in zip(r_edges[:-1], r_edges[1:]):
            sel = (r >= lo) & (r < hi)
            prof.append(float(np.nanmedian(d[sel])) / centre if sel.any() else float("nan"))
        out[str(wl)] = prof
        rr = np.array(out["r"])
        inside = rr < 1.0
        out[f"{wl}_limb_peak"] = float(np.nanmax(np.array(prof)[inside]))
    return out


def corona_scale_heights(maps: dict) -> dict:
    """Hydrostatic temperatures from how fast each channel fades above the limb.

    The limb is cut into ``CORONA_N_SECTORS`` position-angle sectors. In each, the
    median intensity is measured in 2-Mm height bins and fitted with the spherical,
    isothermal hydrostatic model between ``CORONA_FIT_MIN_MM`` and
    ``CORONA_FIT_MAX_MM``. Sectors crossed by an active region (any of 171/193/211 in
    its top quarter at 15 Mm) or a coronal hole (193 in its bottom fifth) are set
    aside; the rest, the quiet corona, give the headline temperatures. The same
    sectors are used for every channel.
    """
    ref = maps[193]
    r = fast_radius(ref)
    pa = sdo.position_angle_map(ref)
    rsun_mm = C.R_SUN_REF_M / 1e6  # the radius the header's RSUN_OBS was computed with
    height = (r - 1.0) * rsun_mm
    h_edges = np.arange(0.0, 202.0, 2.0)
    h_mid = 0.5 * (h_edges[1:] + h_edges[:-1])
    width = 360.0 / C.CORONA_N_SECTORS
    sector = np.floor(pa / width).astype(int)
    hbin = np.digitize(height, h_edges) - 1
    valid = (hbin >= 0) & (hbin < len(h_mid))

    profiles = {}
    for wl in C.CORONA_CHANNELS:
        d = np.asarray(maps[wl].data, dtype=float)
        prof = np.full((C.CORONA_N_SECTORS, len(h_mid)), np.nan)
        for s in range(C.CORONA_N_SECTORS):
            sel = valid & (sector == s)
            idx = hbin[sel]
            vals = d[sel]
            for b in np.unique(idx):
                prof[s, b] = np.median(vals[idx == b])
        profiles[wl] = prof

    fit_rng = (h_mid >= C.CORONA_FIT_MIN_MM) & (h_mid <= C.CORONA_FIT_MAX_MM)
    base_bin = int(np.argmin(np.abs(h_mid - 15.0)))
    base = {wl: profiles[wl][:, base_bin] for wl in C.CORONA_CHANNELS}
    active = np.zeros(C.CORONA_N_SECTORS, bool)
    for wl in (171, 193, 211):
        active |= base[wl] > np.nanpercentile(base[wl], 75)
    hole = base[193] < np.nanpercentile(base[193], 20)
    quiet = ~active & ~hole

    out = {"height_mm": h_mid.tolist(), "fit_range_mm": [C.CORONA_FIT_MIN_MM, C.CORONA_FIT_MAX_MM],
           "sector_width_deg": width, "quiet_sectors": np.flatnonzero(quiet).tolist(),
           "active_sectors": np.flatnonzero(active).tolist(),
           "hole_sectors": np.flatnonzero(hole & ~active).tolist(), "channels": {}}
    for wl in C.CORONA_CHANNELS:
        temps, rms = np.full(C.CORONA_N_SECTORS, np.nan), np.full(C.CORONA_N_SECTORS, np.nan)
        for s in range(C.CORONA_N_SECTORS):
            y = profiles[wl][s, fit_rng]
            if np.sum(np.isfinite(y) & (y > 0)) < 10:
                continue
            ok = np.isfinite(y) & (y > 0)
            temps[s], rms[s], _ = sp.fit_offlimb_temperature(h_mid[fit_rng][ok], y[ok])
        # the quiet corona's median profile, fitted as a whole
        med = np.nanmedian(profiles[wl][quiet], axis=0)
        ok = fit_rng & np.isfinite(med) & (med > 0)
        t_med, rms_med, model = sp.fit_offlimb_temperature(h_mid[ok], med[ok])
        lam_em = sp.emission_scale_height_mm(h_mid[ok], med[ok])
        ch = C.CHANNEL_BY_WAVELENGTH[wl]
        tq = temps[quiet]
        out["channels"][str(wl)] = {
            "ion": ch.ions,
            "log_t_formation": ch.log_t[0],
            "t_formation_MK": 10 ** ch.log_t[0] / 1e6,
            "t_sector_MK": (temps / 1e6).tolist(),
            "rms_sector_dex": rms.tolist(),
            "t_quiet_median_MK": float(np.nanmedian(tq) / 1e6),
            "t_quiet_p16_MK": float(np.nanpercentile(tq, 16) / 1e6),
            "t_quiet_p84_MK": float(np.nanpercentile(tq, 84) / 1e6),
            "t_active_median_MK": float(np.nanmedian(temps[active]) / 1e6),
            "n_quiet": int(np.sum(np.isfinite(tq))),
            "t_quiet_profile_MK": t_med / 1e6,
            "rms_quiet_profile_dex": rms_med,
            "quiet_profile": med.tolist(),
            "quiet_model": model.tolist(),
            "quiet_model_heights": h_mid[ok].tolist(),
            "emission_scale_height_mm": lam_em,
            "t_plain_exponential_MK": float(sp.temperature_from_scale_height(2 * lam_em * 1e6) / 1e6),
            "density_scale_height_mm": float(sp.pressure_scale_height_m(t_med) / 1e6),
        }
    return out
