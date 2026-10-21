"""The problems inflation solves, in numbers (a labelled calculation).

1. **Horizons.** Without inflation, light emitted at the Big Bang has travelled
   a comoving distance eta_* (the conformal time) by the time the CMB is
   released. Seen from today, that distance subtends theta_H = eta_* / D_M(z_*)
   on the sky: points farther apart could never have exchanged a signal, yet
   they have the same temperature to 1 part in 10^5.
2. **The shrinking Hubble sphere.** The comoving Hubble radius 1/(aH) grows in
   any decelerating universe; during inflation it shrinks by e^N. We follow it
   from inside an inflationary model (Starobinsky, solved exactly, instant
   reheating) through the hot Big Bang to today's dark-energy era.
3. **Flatness.** |Omega - 1| = |Omega_k| (a0 H0 / aH)^2 grows with 1/(aH)^2, so
   today's bound on curvature (O1's) traced back gives how finely balanced the
   early universe had to be -- unless inflation drove it there.

The hot Big Bang part needs the number of relativistic species g*(T): here a
free-particle gas of the Standard Model with a smooth quark-hadron crossover
at 156.5 MeV (HotQCD 2019) -- adequate for a plot spanning 60 orders of
magnitude; project P2 builds g*(T) properly with lattice QCD.
"""

from __future__ import annotations

import math
from functools import lru_cache

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad

import config as cfg

# Standard Model species: name, mass [GeV] (PDG 2024), internal dof, fermion?, phase
SPECIES = (
    ("photon", 0.0, 2, False, "always"),
    ("e", 0.000511, 4, True, "always"),
    ("mu", 0.10566, 4, True, "always"),
    ("tau", 1.77693, 4, True, "always"),
    ("nu", 0.0, 6, True, "neutrino"),
    ("W", 80.369, 6, False, "always"),
    ("Z", 91.188, 3, False, "always"),
    ("H", 125.20, 1, False, "always"),
    ("gluon", 0.0, 16, False, "qgp"),
    ("u", 0.00216, 12, True, "qgp"),
    ("d", 0.00470, 12, True, "qgp"),
    ("s", 0.0935, 12, True, "qgp"),
    ("c", 1.2730, 12, True, "qgp"),
    ("b", 4.183, 12, True, "qgp"),
    ("t", 172.57, 12, True, "qgp"),
    ("pi+-", 0.13957, 2, False, "hadron"),
    ("pi0", 0.13498, 1, False, "hadron"),
    ("K", 0.4937, 4, False, "hadron"),
    ("eta", 0.5479, 1, False, "hadron"),
    ("rho", 0.7753, 9, False, "hadron"),
    ("omega", 0.7827, 3, False, "hadron"),
    ("nucleon", 0.9389, 8, True, "hadron"),
)
T_QCD = 0.1565       # GeV, lattice QCD pseudo-critical temperature (HotQCD 2019)
DT_QCD = 0.010       # GeV, width of the smooth switch (an approximation)
T_NU_DEC = 0.0015    # GeV, neutrino decoupling taken as instantaneous (P3 does better)


_X_TAB = np.geomspace(1e-4, 60.0, 1500)


@lru_cache(maxsize=None)
def _tables(fermion: bool) -> tuple[np.ndarray, np.ndarray]:
    """ln(rho/T^4) and ln(P/T^4) per internal degree of freedom on a grid of x = m/T."""
    sgn = 1.0 if fermion else -1.0
    rho, prs = [], []
    for x in _X_TAB:
        f_rho = lambda u: u * u * math.sqrt(max(u * u - x * x, 0.0)) / (math.exp(u) + sgn)
        f_p = lambda u: (max(u * u - x * x, 0.0)) ** 1.5 / (math.exp(u) + sgn)
        rho.append(quad(f_rho, x, x + 80, limit=200)[0] / (2 * math.pi ** 2))
        prs.append(quad(f_p, x, x + 80, limit=200)[0] / (6 * math.pi ** 2))
    return np.log(rho), np.log(prs)


def _free_gas(x: float, fermion: bool) -> tuple[float, float]:
    """rho/T^4 and P/T^4 for one internal degree of freedom of mass x = m/T."""
    if x >= _X_TAB[-1]:
        return 0.0, 0.0
    lx = math.log(max(x, _X_TAB[0]))
    lr, lp = _tables(fermion)
    grid = np.log(_X_TAB)
    return float(np.exp(np.interp(lx, grid, lr))), float(np.exp(np.interp(lx, grid, lp)))


def _plasma(T: float) -> tuple[float, float, float, float]:
    """(g*_rho, g*_s) of everything except neutrinos, and of photons + e+- alone."""
    sw = 0.5 * (1 + math.tanh((T - T_QCD) / DT_QCD))
    rho = s = rho_ge = s_ge = 0.0
    for name, m, g, ferm, phase in SPECIES:
        if phase == "neutrino":
            continue
        weight = sw if phase == "qgp" else (1 - sw) if phase == "hadron" else 1.0
        if weight < 1e-12:
            continue
        r, p = _free_gas(m / T, ferm)
        rho += weight * g * r
        s += weight * g * (r + p)
        if name in ("photon", "e"):
            rho_ge += g * r
            s_ge += g * (r + p)
    k_rho, k_s = 30 / math.pi ** 2, 45 / (2 * math.pi ** 2)
    return rho * k_rho, s * k_s, rho_ge * k_rho, s_ge * k_s


def g_star(T: float) -> tuple[float, float]:
    """(g*_rho, g*_s) at photon temperature T [GeV], neutrinos included."""
    g_rho, g_s, _, gs_ge = _plasma(T)
    if T >= T_NU_DEC:
        tnu = 1.0
    else:  # e+e- annihilation heats the photons but not the decoupled neutrinos; the
        # textbook instantaneous decoupling with relativistic e+e- (g_s = 11/2) gives
        # T_nu/T -> (4/11)^(1/3) today (P3 treats decoupling properly: N_eff = 3.044)
        tnu = min(1.0, (gs_ge / 5.5) ** (1 / 3))
    return g_rho + 7 / 8 * 6 * tnu ** 4, g_s + 7 / 8 * 6 * tnu ** 3


# ---------------------------------------------------------------------------
# 1. The horizon at the release of the CMB
# ---------------------------------------------------------------------------


def cmb_horizon() -> dict:
    """Conformal time at last scattering and the angle it subtends (Planck 2018, CAMB)."""
    import camb

    pars = camb.set_params(H0=cfg.H0_FID, ombh2=cfg.OMBH2_FID, omch2=cfg.OMCH2_FID,
                           mnu=cfg.MNU_FID, tau=cfg.TAU_FID, ns=cfg.NS_FID, As=2.1e-9)
    res = camb.get_background(pars)
    der = res.get_derived_params()
    z_star = der["zstar"]
    eta_star = float(res.conformal_time(z_star))
    eta0 = float(res.tau0)
    dm = eta0 - eta_star
    theta = eta_star / dm
    patch = lambda rad: 4 * math.pi / (2 * math.pi * (1 - math.cos(rad)))
    return {
        "z_star": z_star, "eta_star_mpc": eta_star, "eta0_mpc": eta0, "D_M_star_mpc": dm,
        "theta_H_deg": math.degrees(theta),
        # never exchanged a signal: separation > theta_H; every pair inside a disc of
        # diameter theta_H could have -> about 4 pi / (pi (theta_H/2)^2) such discs
        "patches_exchange": patch(theta / 2),
        # no common causal past at all: separation > 2 theta_H -> discs of radius theta_H
        "patches_common_past": patch(theta),
        "sound_horizon_mpc": der["rstar"], "age_gyr": der["age"],
        "hubble_radius_today_mpc": cfg.C_KM_S / cfg.H0_FID,
    }


# ---------------------------------------------------------------------------
# 2-3. The comoving Hubble radius and |Omega - 1| through the whole history
# ---------------------------------------------------------------------------


def history(model_res: dict, bg, v0: float, *, omega_k_bound: float,
            omega_m: float = 0.3153) -> dict:
    """1/(aH) [comoving Mpc] and |Omega - 1| from inside inflation to the far future.

    ``bg`` is an exact inflationary background (V0 = 1 units) and ``v0`` its scale
    fixed by A_s; reheating is instant, into a Standard-Model plasma.
    """
    mpl = cfg.M_PL_GEV
    gev_to_mpc_inv = 1.0 / cfg.MPC_INV_GEV
    h0 = cfg.H0_FID / cfg.C_KM_S                            # Mpc^-1
    t0 = cfg.T_CMB_K * cfg.K_B_GEV_PER_K                    # GeV
    g0_rho, g0_s = g_star(1e-13)                            # today's (3.36, 3.91)
    # end of inflation -> instant reheating
    phi_e, pi_e = bg.state(bg.n_end)
    v_end = bg.model.v(phi_e)[0]
    rho_end = 1.5 * v0 * v_end * mpl ** 4                   # 3 H^2 M^2 with eps = 1
    g_reh = g_star(1e4)[0]
    t_reh = (30 * rho_end / (math.pi ** 2 * g_reh)) ** 0.25
    a_end = (g0_s / g_star(t_reh)[1]) ** (1 / 3) * t0 / t_reh
    # inflation: a = a_end e^(N - N_end), H from the exact solution
    n = np.linspace(0.0, bg.n_end, 4000)
    ln_h = np.array([bg.ln_h(x) for x in n]) + 0.5 * math.log(v0) + math.log(mpl)  # ln H [GeV]
    ln_a_inf = math.log(a_end) + (n - bg.n_end)
    ln_aH_inf = ln_a_inf + ln_h + math.log(gev_to_mpc_inv)
    # hot Big Bang: temperature grid down to today, then the future (T < T0)
    temps = np.geomspace(t_reh, t0, 1400)
    gs = np.array([g_star(float(T)) for T in temps])
    a_hot = (g0_s / gs[:, 1]) ** (1 / 3) * t0 / temps
    omega_r = (math.pi ** 2 / 30 * g0_rho * t0 ** 4) / (3 * mpl ** 2) / (h0 / gev_to_mpc_inv) ** 2
    rho_r = math.pi ** 2 / 30 * gs[:, 0] * temps ** 4                       # GeV^4
    h2 = rho_r / (3 * mpl ** 2) * gev_to_mpc_inv ** 2 + h0 ** 2 * (
        omega_m * a_hot ** -3 + (1 - omega_m - omega_r))
    ln_aH_hot = np.log(a_hot) + 0.5 * np.log(h2)
    a_fut = np.geomspace(1.0, 30.0, 200)[1:]
    h2f = h0 ** 2 * (omega_r * a_fut ** -4 + omega_m * a_fut ** -3 + 1 - omega_m - omega_r)
    ln_aH_fut = np.log(a_fut) + 0.5 * np.log(h2f)
    # without inflation: the hot Big Bang extrapolated back to the Planck temperature
    t_pl = 1.220890e19
    temps_ni = np.geomspace(t_pl, t_reh, 300)
    gs_ni = np.array([g_star(float(T)) for T in temps_ni])
    a_ni = (g0_s / gs_ni[:, 1]) ** (1 / 3) * t0 / temps_ni
    h_ni = np.sqrt(math.pi ** 2 / 30 * gs_ni[:, 0] * temps_ni ** 4 / (3 * mpl ** 2)) * gev_to_mpc_inv
    ln_aH_ni = np.log(a_ni) + np.log(h_ni)
    ln_a = np.concatenate([ln_a_inf, np.log(a_hot), np.log(a_fut)])
    ln_aH = np.concatenate([ln_aH_inf, ln_aH_hot, ln_aH_fut])
    phase = np.concatenate([np.zeros(n.size), np.ones(temps.size), 2 * np.ones(a_fut.size)])
    temp = np.concatenate([np.full(n.size, np.nan), temps, np.full(a_fut.size, np.nan)])
    omega_dev = omega_k_bound * np.exp(2 * (math.log(h0) - ln_aH))
    # e-folds needed. aH grows monotonically during inflation (the comoving Hubble
    # sphere shrinks); find when it passed today's a0 H0 (horizon problem: all we
    # see today was once inside one Hubble sphere) and a0 H0 sqrt(|Omega_k|)
    # (flatness: |Omega - 1| was still ~1 then)
    def efolds_before_end(target):
        if not (ln_aH_inf[0] <= target <= ln_aH_inf[-1]):
            return None
        return float(bg.n_end - np.interp(target, ln_aH_inf, n))

    n_min_h = efolds_before_end(math.log(h0))
    n_min_f = efolds_before_end(math.log(h0) + 0.5 * math.log(omega_k_bound))
    # pivot scale check: k_* = aH at N_end - N_*
    n_star = model_res["n_star"]
    ln_k_star = float(np.interp(bg.n_end - n_star, n, ln_aH_inf))

    def at_temp(T):
        return float(np.interp(math.log(T), np.log(temps[::-1]), np.log(omega_dev[n.size:n.size + temps.size][::-1])))

    return {
        "arrays": {"ln_a": ln_a, "ln_inv_aH_mpc": -ln_aH, "omega_dev": omega_dev,
                   "phase": phase, "temperature_gev": temp,
                   "g_temps": temps, "g_rho": gs[:, 0], "g_s": gs[:, 1],
                   "noinf_ln_a": np.log(a_ni), "noinf_ln_inv_aH_mpc": -ln_aH_ni,
                   "noinf_omega_dev": omega_k_bound * np.exp(2 * (math.log(h0) - ln_aH_ni)),
                   "noinf_temperature_gev": temps_ni},
        "summary": {
            "T_reh_gev": t_reh, "a_end": a_end, "g_reh": g_reh, "rho_end_gev4": rho_end,
            "g_star_today": [g0_rho, g0_s], "omega_r": omega_r,
            "k_star_check_mpc_inv": math.exp(ln_k_star), "n_star": n_star,
            "inv_aH_end_mpc": math.exp(-ln_aH_inf[-1]),
            "n_min_horizon": n_min_h, "n_min_flatness": n_min_f,
            "omega_dev_at": {"BBN_1MeV": math.exp(at_temp(1e-3)),
                             "QCD_156MeV": math.exp(at_temp(0.1565)),
                             "EW_160GeV": math.exp(at_temp(160.0)),
                             "reheating": float(omega_dev[n.size])},
            "omega_k_bound": omega_k_bound,
        },
    }


def omega_dev_planck_temperature(omega_k_bound: float, omega_m: float = 0.3153) -> float:
    """|Omega - 1| extrapolated back to T = 1.22e19 GeV (the Planck energy) assuming a
    hot Big Bang all the way, with the Standard Model's g* there."""
    mpl = cfg.M_PL_GEV
    gev_to_mpc_inv = 1.0 / cfg.MPC_INV_GEV
    h0 = cfg.H0_FID / cfg.C_KM_S
    t0 = cfg.T_CMB_K * cfg.K_B_GEV_PER_K
    tp = 1.220890e19
    g0_rho, g0_s = g_star(1e-13)
    g_rho, g_s = g_star(tp)
    a = (g0_s / g_s) ** (1 / 3) * t0 / tp
    h = math.sqrt(math.pi ** 2 / 30 * g_rho * tp ** 4 / (3 * mpl ** 2)) * gev_to_mpc_inv
    return omega_k_bound * (h0 / (a * h)) ** 2


# ---------------------------------------------------------------------------
# The CMB sky, for the picture
# ---------------------------------------------------------------------------


def gnomonic_forward(lon, lat, lon0, lat0):
    """Gnomonic projection (degrees in, tangent-plane coordinates in degrees of arc at
    the centre out): x toward increasing longitude, y toward increasing latitude."""
    lam, phi = np.radians(lon), np.radians(lat)
    lam0, phi0 = math.radians(lon0), math.radians(lat0)
    cosc = math.sin(phi0) * np.sin(phi) + math.cos(phi0) * np.cos(phi) * np.cos(lam - lam0)
    x = np.cos(phi) * np.sin(lam - lam0) / cosc
    y = (math.cos(phi0) * np.sin(phi) - math.sin(phi0) * np.cos(phi) * np.cos(lam - lam0)) / cosc
    return np.degrees(x), np.degrees(y)


def gnomonic_inverse(x, y, lon0, lat0):
    """Inverse of gnomonic_forward (x, y in degrees of arc at the centre)."""
    xr, yr = np.radians(x), np.radians(y)
    phi0, lam0 = math.radians(lat0), math.radians(lon0)
    rho = np.hypot(xr, yr)
    c = np.arctan(rho)
    with np.errstate(invalid="ignore", divide="ignore"):
        lat = np.arcsin(np.cos(c) * math.sin(phi0) + np.where(rho > 0, yr * np.sin(c) * math.cos(phi0) / rho, 0))
        lon = lam0 + np.arctan2(xr * np.sin(c), rho * math.cos(phi0) * np.cos(c) - yr * math.sin(phi0) * np.sin(c))
    return np.degrees(lon), np.degrees(lat)


def cmb_sky(nside: int = 1024, fwhm_arcmin: float = 15.0, zoom_center=(240.0, -58.0),
            zoom_deg: float = 20.0, zoom_px: int = 800) -> dict:
    """Planck SMICA (mean of the two half-mission maps), smoothed: a Mollweide image
    (Galactic coordinates; Planck's common mask set to NaN) and a gnomonic zoom on
    a clean high-latitude patch, sampled on our own tangent-plane grid."""
    import healpy as hp

    maps = [hp.read_map(cfg.SKY_DIR / f, field=0, dtype=np.float32) for f in cfg.SMICA_FILES]
    m = 0.5 * (maps[0] + maps[1]) * 1e6                    # K -> uK
    del maps
    mask = hp.read_map(cfg.SKY_DIR / cfg.MASK_FILE, field=0, dtype=np.float32)
    m = hp.ud_grade(m, nside)
    mask = hp.ud_grade(mask, nside)
    m = hp.remove_monopole(m, gal_cut=30)
    alm = hp.map2alm(m, lmax=3 * nside - 1, use_pixel_weights=False)
    hp.almxfl(alm, hp.gauss_beam(math.radians(fwhm_arcmin / 60), lmax=3 * nside - 1), inplace=True)
    sm = hp.alm2map(alm, nside)
    masked = np.where(mask > 0.5, sm, np.nan)
    moll = hp.projector.MollweideProj(xsize=1600)
    img = moll.projmap(masked, lambda x, y, z: hp.vec2pix(nside, x, y, z))
    half = zoom_deg / 2
    g = np.linspace(-half, half, zoom_px)
    xx, yy = np.meshgrid(g, g)
    lon, lat = gnomonic_inverse(xx, yy, *zoom_center)
    zoom = hp.get_interp_val(sm, lon, lat, lonlat=True)
    zmask = hp.get_interp_val(mask, lon, lat, lonlat=True)
    # the zoom patch's outline in the Mollweide image's own coordinates
    edge = np.concatenate([np.linspace(-half, half, 50)])
    ex = np.concatenate([edge, np.full(50, half), edge[::-1], np.full(50, -half)])
    ey = np.concatenate([np.full(50, -half), edge, np.full(50, half), edge[::-1]])
    olon, olat = gnomonic_inverse(ex, ey, *zoom_center)
    mx, my = moll.ang2xy(olon, olat, lonlat=True)
    return {"mollweide": img.astype(np.float32), "zoom": zoom.astype(np.float32),
            "zoom_center": np.array(zoom_center), "zoom_deg": zoom_deg,
            "zoom_mask_fraction": float(np.mean(zmask < 0.5)),
            "moll_outline": np.vstack([mx, my]), "moll_extent": np.array(moll.get_extent()),
            "rms_uK": float(np.nanstd(masked)), "fwhm_arcmin": fwhm_arcmin,
            "masked_fraction": float(np.mean(mask < 0.5)), "nside": nside}
