"""The pipeline: from 175 rotation curves to the case for dark matter.

It works at two levels.

**One galaxy at a time** (``analyse_galaxy``). Where does the observed rotation
part company with what the visible matter predicts, and by how much? Could the
stars simply be heavier than assumed (``maximum_disk``)? What shape of dark
halo fills the gap, and is it the shape cold-dark-matter simulations predict
(``fit_halo``)?

**The whole sample** (``analyse_population``): the 153 galaxies that pass the
quality cuts. The radial acceleration relation -- at every radius of every
galaxy, the observed acceleration against the visible matter's. The baryonic
Tully-Fisher relation -- each galaxy's total visible mass against its flat
rotation speed. And the maximum-disc test, repeated for every galaxy.

Published values are loaded only in ``compare_with_literature``, after our own
numbers exist, and never feed into them.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, is_dataclass
from datetime import datetime, timezone

import numpy as np
from scipy import optimize, stats
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from scipy.spatial import cKDTree

import rotation_physics as rp
import sparc
from config import (
    H0,
    HI_MASS_ERROR,
    MAX_DISK_TOLERANCE_SIGMA,
    MAX_FRACTIONAL_VELOCITY_ERROR,
    MAX_QUALITY_FLAG,
    MIN_INCLINATION_DEG,
    RESULTS_PATH,
    UPSILON_BULGE,
    UPSILON_DISK,
    UPSILON_PLAUSIBLE,
    UPSILON_SCATTER_DEX,
    Showcase,
)

#: Bootstrap resamples for uncertainties that should respect the fact that the
#: points within one galaxy are not independent. Fixed seed: reproducible.
N_BOOTSTRAP = 1000
RNG_SEED = 20160920

#: The published numbers we grade ourselves against. Used only in
#: ``compare_with_literature``.
PUBLISHED = {
    "rar": {
        "reference": "McGaugh, Lelli & Schombert 2016, PRL 117, 201101",
        "g_dagger": 1.20e-10,
        "g_dagger_random": 0.02e-10,
        "g_dagger_systematic": 0.24e-10,
        "rms_dex": 0.13,
        "gaussian_width_dex": 0.11,
        "n_points": 2693,
        "n_galaxies": 153,
    },
    "btfr": {
        "reference": "Lelli et al. 2019, MNRAS 484, 3267",
        "slope": 3.85,
        "slope_err": 0.09,
    },
    "max_disk": {"reference": "Starkman et al. 2018, MNRAS 480, 2292"},
}


# --------------------------------------------------------------------------
# One galaxy
# --------------------------------------------------------------------------

def baryon_components(g: sparc.Galaxy, upsilon_disk: float = UPSILON_DISK,
                      upsilon_bulge: float = UPSILON_BULGE) -> dict[str, np.ndarray]:
    """Each visible component's rotation speed at the given mass-to-light ratios.

    Speeds are signed (see ``rotation_physics.signed_square``); ``baryons_v2``
    is the signed square of the total, which is what adds up.
    """
    c = g.curve
    v2 = rp.baryonic_velocity_squared(c.v_gas, c.v_disk, c.v_bulge, upsilon_disk, upsilon_bulge)
    return {
        "gas": c.v_gas,
        "disk": np.sqrt(upsilon_disk) * c.v_disk,
        "bulge": np.sqrt(upsilon_bulge) * c.v_bulge,
        "baryons": rp.signed_sqrt(v2),
        "baryons_v2": v2,
    }


def dark_matter_required(g: sparc.Galaxy) -> dict[str, np.ndarray]:
    """What the unseen mass must contribute at each radius, straight from the data.

    Speeds add in quadrature, so whatever holds the observed orbit beyond the
    visible matter must supply ``V_dm**2 = V_obs**2 - V_bar**2``. No halo model
    is involved: this is the gap itself, in velocity units. The error bar maps
    the observed speed's +/-1 sigma through the same subtraction, which keeps
    it honest where the gap is small (and can dip below zero, where the visible
    matter alone already over-predicts).
    """
    c = g.curve
    v2_bar = baryon_components(g)["baryons_v2"]
    v2 = c.v_obs**2 - v2_bar
    lo = (c.v_obs - c.e_v_obs) ** 2 - v2_bar
    hi = (c.v_obs + c.e_v_obs) ** 2 - v2_bar
    return {"v": rp.signed_sqrt(v2), "v_lo": rp.signed_sqrt(lo), "v_hi": rp.signed_sqrt(hi)}


@dataclass
class MaxDisk:
    """Could heavier stars do the job without dark matter? One galaxy's answer.

    Stellar mass-to-light ratios are expressed as a multiple ``scale`` of the
    fiducial ones (Y_disk = 0.5, Y_bulge = 0.7), so disc and bulge scale
    together.
    """

    #: The heaviest the stars can be before the visible matter alone spins the
    #: galaxy faster than observed (+ tolerance) somewhere: the "maximum disc".
    scale_max: float
    upsilon_disk_max: float
    #: The radius where that limit bites.
    binding_radius_kpc: float
    #: How heavy the stars must be for visible matter alone to reach the observed
    #: speed (- tolerance) at the outermost measured point.
    scale_needed: float
    upsilon_disk_needed: float
    #: ``scale_needed / scale_max``. Above 1, *no* single mass-to-light ratio,
    #: however implausible, lets stars and gas explain both the inner and the
    #: outer rotation curve.
    overshoot: float
    #: With the stars heavy enough to explain the edge, how much faster than
    #: observed the visible matter would spin the galaxy somewhere inside.
    inner_excess: float
    #: Whether the needed mass-to-light ratio is one stars can actually have.
    needed_is_plausible: bool
    tolerance_sigma: float


def maximum_disk(g: sparc.Galaxy, tolerance_sigma: float = MAX_DISK_TOLERANCE_SIGMA) -> MaxDisk | None:
    """The maximum-disc test for one galaxy.

    The single largest uncertainty in the visible-matter prediction is the
    mass-to-light ratio of the stars. So give the no-dark-matter hypothesis
    every chance: scale the stars up as far as the *inner* rotation curve
    allows, and ask whether that is enough to hold the *outer* curve.

    Two robustness notes. The tolerance runs in the hypothesis's favour at both
    ends (the ceiling is the observed speed plus ``tolerance_sigma`` errors, the
    target at the edge is the observed speed minus them). And because an error
    in the galaxy's distance or inclination rescales the stellar term or the
    observed speeds uniformly at every radius, it moves ``scale_max`` and
    ``scale_needed`` by the same factor: ``overshoot`` compares the *shapes* of
    the visible and observed curves, so it is immune to both.
    """
    c = g.curve
    stars2 = UPSILON_DISK * rp.signed_square(c.v_disk) + UPSILON_BULGE * rp.signed_square(c.v_bulge)
    gas2 = rp.signed_square(c.v_gas)
    ok = stars2 > 0
    if not ok.any():
        return None

    ceiling = (c.v_obs + tolerance_sigma * c.e_v_obs) ** 2 - gas2
    ratios = ceiling[ok] / stars2[ok]
    k = int(np.argmin(ratios))
    scale_max = max(float(ratios[k]), 0.0)

    target = (max(c.v_obs[-1] - tolerance_sigma * c.e_v_obs[-1], 0.0)) ** 2 - gas2[-1]
    if stars2[-1] > 0:
        scale_needed = max(float(target / stars2[-1]), 0.0)
    else:
        scale_needed = float("inf") if target > 0 else 0.0

    if np.isfinite(scale_needed):
        v_needed = rp.signed_sqrt(gas2 + scale_needed * stars2)
        inner_excess = float(np.max(v_needed / c.v_obs))
    else:
        inner_excess = float("inf")
    y_needed = UPSILON_DISK * scale_needed
    return MaxDisk(
        scale_max=scale_max,
        upsilon_disk_max=UPSILON_DISK * scale_max,
        binding_radius_kpc=float(c.r[ok][k]),
        scale_needed=scale_needed,
        upsilon_disk_needed=y_needed,
        overshoot=(scale_needed / scale_max) if scale_max > 0 else float("inf"),
        inner_excess=inner_excess,
        needed_is_plausible=bool(UPSILON_PLAUSIBLE[0] <= y_needed <= UPSILON_PLAUSIBLE[1]),
        tolerance_sigma=tolerance_sigma,
    )


# --- dark matter haloes ----------------------------------------------------

#: Parameters are fitted as log10 values: it keeps them positive, and puts
#: parameters that span orders of magnitude on an even footing. The start
#: points are a coarse grid; the best of all the local fits is kept.
HALO_MODELS = {
    "NFW": {
        "params": ("v200_kms", "concentration"),
        "bounds": ((1.0, 0.0), (3.0, 2.0)),        # V200 10-1000 km/s, c 1-100
        "starts": [(np.log10(v), np.log10(cc)) for v in (50, 120, 250) for cc in (3, 10, 30)],
    },
    "ISO": {
        "params": ("v_inf_kms", "r_core_kpc"),
        "bounds": ((0.5, -2.0), (3.0, 2.5)),       # V_inf 3-1000 km/s, r_c 0.01-300 kpc
        "starts": [(np.log10(v), np.log10(rc)) for v in (40, 100, 250) for rc in (0.3, 1, 3, 10)],
    },
}


def halo_velocity(model: str, r_kpc, params: dict) -> np.ndarray:
    if model == "NFW":
        return rp.v_nfw(r_kpc, params["v200_kms"], params["concentration"], H0)
    if model == "ISO":
        return rp.v_iso(r_kpc, params["v_inf_kms"], params["r_core_kpc"])
    raise KeyError(model)


@dataclass
class HaloFit:
    model: str
    params: dict
    #: 1-sigma errors from the curvature of chi^2, inflated by sqrt(chi2_nu)
    #: when the fit is worse than the error bars allow. Halo parameters are
    #: strongly correlated, so treat these as indicative only.
    errors: dict
    chi2: float
    dof: int
    chi2_nu: float
    #: True if a parameter ended on the edge of its allowed range -- a sign the
    #: model wants to be somewhere unphysical, or is unconstrained.
    at_bound: bool
    #: The misfit in plain km/s. chi2 can be large simply because a galaxy's
    #: error bars are tiny (SPARC's include random errors only, not the few-km/s
    #: wiggles non-circular gas motions add), so this is the better guide to
    #: whether a fit is actually good.
    rms_residual_kms: float
    max_abs_residual_kms: float
    derived: dict


def fit_halo(g: sparc.Galaxy, model: str) -> HaloFit:
    """Fit a dark halo to fill the gap, with the stars fixed at Y = 0.5.

    ``V_model**2 = V_bar**2 + V_halo**2`` is fitted to the observed curve,
    weighted by the published errors. The stellar mass-to-light ratio is held
    at its population-synthesis value rather than fitted: letting it float
    trades disc mass against halo mass (the "disc-halo degeneracy") and would
    tell us more about that degeneracy than about the halo.
    """
    spec = HALO_MODELS[model]
    c = g.curve
    v2_bar = baryon_components(g)["baryons_v2"]
    lo, hi = (np.array(b, dtype=float) for b in spec["bounds"])

    def to_params(logp):
        return dict(zip(spec["params"], 10.0 ** np.asarray(logp)))

    def residuals(logp):
        vh = halo_velocity(model, c.r, to_params(logp))
        v_model = np.sqrt(np.clip(v2_bar + vh**2, 0.0, None))
        return (v_model - c.v_obs) / c.e_v_obs

    best = None
    for start in spec["starts"]:
        x0 = np.clip(start, lo + 1e-6, hi - 1e-6)
        sol = optimize.least_squares(residuals, x0, bounds=(lo, hi), x_scale="jac")
        if best is None or sol.cost < best.cost:
            best = sol

    chi2 = float(2.0 * best.cost)
    dof = max(len(c) - len(spec["params"]), 1)
    chi2_nu = chi2 / dof
    params = to_params(best.x)
    try:
        cov = np.linalg.inv(best.jac.T @ best.jac) * max(1.0, chi2_nu)
        sig_log = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    except np.linalg.LinAlgError:
        sig_log = np.full(len(best.x), np.nan)
    errors = {k: float(params[k] * np.log(10.0) * s) for k, s in zip(spec["params"], sig_log)}
    at_bound = bool(np.any(np.abs(best.x - lo) < 1e-3) or np.any(np.abs(best.x - hi) < 1e-3))

    derived: dict = {}
    if model == "NFW":
        v200, conc = params["v200_kms"], params["concentration"]
        m200 = float(rp.m200_msun(v200, H0))
        c_lcdm = float(rp.nfw_concentration_lcdm(m200))
        derived = {
            "r200_kpc": float(rp.r200_kpc(v200, H0)),
            "r_s_kpc": float(rp.r200_kpc(v200, H0) / conc),
            "m200_msun": m200,
            "concentration_lcdm": c_lcdm,
            #: How far the fitted concentration sits from the simulations'
            #: median for this mass, in dex and in units of their 0.11 dex scatter.
            "concentration_offset_dex": float(np.log10(conc / c_lcdm)),
            "concentration_offset_sigma": float(np.log10(conc / c_lcdm) / 0.11),
        }
    else:
        derived = {"rho0_msun_pc3": float(rp.iso_central_density(params["v_inf_kms"],
                                                                  params["r_core_kpc"]))}

    vh_last = float(halo_velocity(model, c.r[-1:], params)[0])
    v2_tot_last = v2_bar[-1] + vh_last**2
    derived["halo_fraction_last"] = float(vh_last**2 / v2_tot_last) if v2_tot_last > 0 else None
    misfit = residuals(best.x) * c.e_v_obs

    return HaloFit(model=model, params={k: float(v) for k, v in params.items()}, errors=errors,
                   chi2=chi2, dof=dof, chi2_nu=chi2_nu, at_bound=at_bound,
                   rms_residual_kms=float(np.sqrt(np.mean(misfit**2))),
                   max_abs_residual_kms=float(np.max(np.abs(misfit))), derived=derived)


# --- putting one galaxy together -------------------------------------------

@dataclass
class GalaxyResult:
    name: str
    label: str
    kind: str
    why: str
    hubble_type: str
    distance_mpc: float
    distance_method: str
    inclination_deg: float
    quality: int
    n_points: int
    r_last_kpc: float
    v_obs_last: float
    e_v_obs_last: float
    #: What the visible matter predicts at the outermost point.
    v_bar_last: float
    v_flat: float | None
    #: Catalogue totals, Msun: stars at the fiducial Y, gas = 1.33 x HI.
    stellar_mass: float
    gas_mass: float
    baryonic_mass: float
    gas_fraction: float
    sb_disk0: float
    #: Mass gravity requires inside the outermost point (spherical-equivalent,
    #: V_obs**2 r / G), Msun -- and how many times the galaxy's *entire* visible
    #: mass that is. Comparing with the catalogue total rather than with
    #: V_bar**2 r / G sidesteps disc geometry: a flat disc's spherical-equivalent
    #: mass overstates what it actually contains (for NGC 3198 it exceeds the
    #: galaxy's whole visible mass), whereas a total is a total.
    dynamical_mass_last: float
    dynamical_over_all_visible_last: float
    #: Radius beyond which the orbits need more mass than every star and gas
    #: cloud in the galaxy put together; ``None`` if that never happens.
    needs_more_than_all_visible_from_kpc: float | None
    #: (V_obs / V_bar)**2 at the outermost point: how much more gravity there is
    #: than the visible matter can account for.
    mass_ratio_last: float
    #: 1 - (V_bar / V_obs)**2 at the outermost point: the fraction of the pull
    #: there that the visible matter does not supply.
    dark_fraction_last: float
    #: Radius beyond which the unseen mass supplies most of the pull at every
    #: measured point; ``None`` if that never happens.
    dark_dominates_from_kpc: float | None
    dark_dominates_everywhere: bool
    vbar_peak: float
    r_vbar_peak_kpc: float
    #: V_bar / V_obs where the visible matter peaks. Near 0.85-0.9 a disc is
    #: called "maximal"; far below that, dark matter dominates even the centre.
    vbar_over_vobs_at_peak: float
    max_disk: dict
    halo_fits: dict
    better_halo: str


def analyse_galaxy(g: sparc.Galaxy, showcase: Showcase | None = None) -> tuple[GalaxyResult, dict]:
    """Every per-galaxy result, plus the arrays the figures are drawn from."""
    c = g.curve
    comp = baryon_components(g)
    v2_bar = comp["baryons_v2"]
    frac_visible = v2_bar / c.v_obs**2

    dominated = frac_visible < 0.5
    dark_from = None
    everywhere = bool(dominated.all())
    if dominated[-1]:
        i = len(dominated) - 1
        while i > 0 and dominated[i - 1]:
            i -= 1
        dark_from = float(c.r[i])

    k_peak = int(np.argmax(comp["baryons"]))
    stellar = UPSILON_DISK * g.luminosity_disk + UPSILON_BULGE * g.luminosity_bulge
    gas = rp.HELIUM_FACTOR * g.m_hi
    m_bar = rp.baryonic_mass(g.luminosity_disk, g.luminosity_bulge, g.m_hi, UPSILON_DISK, UPSILON_BULGE)

    md = maximum_disk(g)
    fits = {m: fit_halo(g, m) for m in HALO_MODELS}
    better = min(fits, key=lambda m: fits[m].chi2)

    m_dyn = rp.enclosed_mass(c.v_obs, c.r)
    exceeds = m_dyn > m_bar
    beyond_all_visible = None
    if exceeds[-1]:
        i = len(exceeds) - 1
        while i > 0 and exceeds[i - 1]:
            i -= 1
        beyond_all_visible = float(c.r[i])

    result = GalaxyResult(
        name=g.name,
        label=g.label,
        kind=showcase.kind if showcase else "",
        why=showcase.why if showcase else "",
        hubble_type=g.type_name,
        distance_mpc=g.distance_mpc,
        distance_method=g.distance_method_name,
        inclination_deg=g.inclination_deg,
        quality=g.quality,
        n_points=len(c),
        r_last_kpc=float(c.r[-1]),
        v_obs_last=float(c.v_obs[-1]),
        e_v_obs_last=float(c.e_v_obs[-1]),
        v_bar_last=float(comp["baryons"][-1]),
        v_flat=g.v_flat,
        stellar_mass=float(stellar),
        gas_mass=float(gas),
        baryonic_mass=float(m_bar),
        gas_fraction=float(gas / m_bar),
        sb_disk0=g.sb_disk0,
        dynamical_mass_last=float(m_dyn[-1]),
        dynamical_over_all_visible_last=float(m_dyn[-1] / m_bar),
        needs_more_than_all_visible_from_kpc=beyond_all_visible,
        mass_ratio_last=float(1.0 / frac_visible[-1]),
        dark_fraction_last=float(1.0 - frac_visible[-1]),
        dark_dominates_from_kpc=dark_from,
        dark_dominates_everywhere=everywhere,
        vbar_peak=float(comp["baryons"][k_peak]),
        r_vbar_peak_kpc=float(c.r[k_peak]),
        vbar_over_vobs_at_peak=float(comp["baryons"][k_peak] / c.v_obs[k_peak]),
        max_disk=asdict(md) if md else {},
        halo_fits={m: asdict(f) for m, f in fits.items()},
        better_halo=better,
    )

    # --- arrays for the figures -----------------------------------------
    r_fine = np.linspace(0.0, c.r[-1], 400)[1:]
    v2_bar_fine = np.interp(r_fine, c.r, v2_bar)
    halo_curves = {}
    for m, f in fits.items():
        vh = halo_velocity(m, r_fine, f.params)
        halo_curves[m] = {"halo": vh, "total": np.sqrt(np.clip(v2_bar_fine + vh**2, 0, None))}
    r_kep = np.linspace(0.25 * c.r[-1], 1.25 * c.r[-1], 200)
    art = {
        "galaxy": g,
        "components": comp,
        "dark_required": dark_matter_required(g),
        "r_fine": r_fine,
        "halo_curves": halo_curves,
        "r_kepler": r_kep,
        "v_kepler": rp.keplerian_velocity(m_bar, r_kep),
        "mass_dynamical": m_dyn,
    }
    if md:
        stars2 = UPSILON_DISK * rp.signed_square(c.v_disk) + UPSILON_BULGE * rp.signed_square(c.v_bulge)
        gas2 = rp.signed_square(c.v_gas)
        art["max_disk_curves"] = {
            "max": rp.signed_sqrt(gas2 + md.scale_max * stars2),
            "needed": (rp.signed_sqrt(gas2 + md.scale_needed * stars2)
                       if np.isfinite(md.scale_needed) else None),
        }
    return result, art


# --------------------------------------------------------------------------
# The whole sample
# --------------------------------------------------------------------------

def select_sample(galaxies: dict[str, sparc.Galaxy]) -> list[sparc.Galaxy]:
    """The galaxies passing the McGaugh et al. (2016) cuts: 153 of 175."""
    return [g for g in galaxies.values()
            if g.quality <= MAX_QUALITY_FLAG and g.inclination_deg >= MIN_INCLINATION_DEG]


@dataclass
class RARData:
    """Every usable (g_bar, g_obs) point, with the galaxy it came from."""

    galaxy: np.ndarray
    r_kpc: np.ndarray
    log_gbar: np.ndarray
    log_gobs: np.ndarray
    e_log_gbar: np.ndarray
    e_log_gobs: np.ndarray

    def __len__(self) -> int:
        return len(self.log_gbar)

    def subset(self, mask) -> RARData:
        return RARData(*(getattr(self, f)[mask] for f in (
            "galaxy", "r_kpc", "log_gbar", "log_gobs", "e_log_gbar", "e_log_gobs")))


def rar_data(sample: list[sparc.Galaxy], upsilon_scale: float = 1.0) -> RARData:
    """Observed and baryonic accelerations at every point that passes the cuts.

    ``upsilon_scale`` multiplies both stellar mass-to-light ratios, for the
    systematic test. Points are kept if their velocity is known to 10% and the
    visible matter's net pull is inward (a handful of innermost points in
    galaxies with a central gas hole have a net outward pull and no logarithm).
    """
    cols: dict[str, list] = {k: [] for k in ("galaxy", "r", "x", "y", "ex", "ey")}
    ud, ub = UPSILON_DISK * upsilon_scale, UPSILON_BULGE * upsilon_scale
    for g in sample:
        c = g.curve
        v2_bar = rp.baryonic_velocity_squared(c.v_gas, c.v_disk, c.v_bulge, ud, ub)
        keep = (c.e_v_obs / c.v_obs < MAX_FRACTIONAL_VELOCITY_ERROR) & (v2_bar > 0) & (c.r > 0)
        if not keep.any():
            continue
        stars2 = ud * rp.signed_square(c.v_disk) + ub * rp.signed_square(c.v_bulge)
        cols["galaxy"] += [g.name] * int(keep.sum())
        cols["r"].append(c.r[keep])
        cols["x"].append(np.log10(rp.acceleration(np.sqrt(v2_bar[keep]), c.r[keep])))
        cols["y"].append(np.log10(rp.acceleration(c.v_obs[keep], c.r[keep])))
        cols["ex"].append(rp.log_gbar_error(stars2[keep] / v2_bar[keep], UPSILON_SCATTER_DEX,
                                            HI_MASS_ERROR))
        cols["ey"].append(rp.log_gobs_error(c.v_obs[keep], c.e_v_obs[keep], g.inclination_deg,
                                            g.inclination_err_deg, g.distance_mpc,
                                            g.distance_err_mpc))
    return RARData(np.array(cols["galaxy"]), *(np.concatenate(cols[k]) for k in
                                               ("r", "x", "y", "ex", "ey")))


def _rar_model(log_gbar, log_gd):
    return np.log10(rp.rar_mcgaugh(10.0**log_gbar, 10.0**log_gd))


def fit_g_dagger(x, y, ex, ey) -> dict:
    """Fit the acceleration scale g_dagger, with errors on both axes.

    McGaugh et al. used orthogonal distance regression. We use its standard
    first-order equivalent, the *effective variance* method: each point's
    x-error is converted to a y-error through the local slope of the relation,
    and ``chi2 = sum (y - f(x))**2 / (ey**2 + slope**2 * ex**2)`` is minimised.
    (SciPy's ODR wrapper is deprecated; on the published data the two methods
    agree to 0.5%.)
    """
    def chi2(lg):
        s = rp.rar_log_slope(10.0**x, 10.0**lg)
        return float(np.sum((y - _rar_model(x, lg)) ** 2 / (ey**2 + s**2 * ex**2)))

    sol = optimize.minimize_scalar(chi2, bounds=(-11.0, -9.0), method="bounded",
                                   options={"xatol": 1e-7})
    lg = float(sol.x)
    n = len(x)
    chi2_nu = sol.fun / (n - 1)
    # Formal error from the curvature of chi^2 (a parabola near the minimum),
    # inflated by sqrt(chi2_nu) since the scatter exceeds the error bars.
    h = 1e-3
    curv = (chi2(lg + h) - 2 * sol.fun + chi2(lg - h)) / h**2
    sig_lg = np.sqrt(2.0 / curv) * np.sqrt(max(chi2_nu, 1.0))
    return {"log_g_dagger": lg, "g_dagger": 10.0**lg, "chi2_nu": float(chi2_nu),
            "g_dagger_err_formal": float(10.0**lg * np.log(10.0) * sig_lg)}


def fit_rar(data: RARData, sample: list[sparc.Galaxy]) -> dict:
    """The radial acceleration relation: fit, scatter, and honest error bars."""
    main = fit_g_dagger(data.log_gbar, data.log_gobs, data.e_log_gbar, data.e_log_gobs)
    lg = main["log_g_dagger"]
    resid = data.log_gobs - _rar_model(data.log_gbar, lg)

    # Unweighted least squares in log g_obs, for reference: no error model at all.
    ols = optimize.minimize_scalar(
        lambda t: float(np.sum((data.log_gobs - _rar_model(data.log_gbar, t)) ** 2)),
        bounds=(-11.0, -9.0), method="bounded")

    # Bootstrap over *galaxies*. The formal error treats 2700 points as
    # independent, but points in one galaxy share its distance, inclination and
    # stellar mass-to-light ratio. Resampling whole galaxies respects that.
    rng = np.random.default_rng(RNG_SEED)
    names = np.unique(data.galaxy)
    index = {n: np.flatnonzero(data.galaxy == n) for n in names}
    boot = []
    for _ in range(N_BOOTSTRAP):
        pick = np.concatenate([index[n] for n in rng.choice(names, size=len(names))])
        boot.append(fit_g_dagger(data.log_gbar[pick], data.log_gobs[pick],
                                 data.e_log_gbar[pick], data.e_log_gobs[pick])["g_dagger"])
    boot = np.array(boot)

    # Systematic: the mass-to-light normalisation is uncertain by ~20%.
    systematic = {}
    for scale in (0.8, 1.2):
        d = rar_data(sample, upsilon_scale=scale)
        systematic[f"{scale:.1f}"] = fit_g_dagger(d.log_gbar, d.log_gobs, d.e_log_gbar,
                                                  d.e_log_gobs)["g_dagger"]

    mad = np.median(np.abs(resid - np.median(resid)))
    return {
        **main,
        "g_dagger_err_bootstrap": float(np.std(boot)),
        "g_dagger_bootstrap_16_84": [float(v) for v in np.percentile(boot, [16, 84])],
        "g_dagger_upsilon_minus20": systematic["0.8"],
        "g_dagger_upsilon_plus20": systematic["1.2"],
        "g_dagger_ols": float(10.0**ols.x),
        "n_points": len(data),
        "n_galaxies": int(len(names)),
        "residual_mean_dex": float(np.mean(resid)),
        "residual_rms_dex": float(np.sqrt(np.mean(resid**2))),
        "residual_gaussian_width_dex": float(1.4826 * mad),
        "median_error_gbar_dex": float(np.median(data.e_log_gbar)),
        "median_error_gobs_dex": float(np.median(data.e_log_gobs)),
        # Where does the discrepancy set in? The observed pull exceeds the
        # visible matter's by 10% where f(g_bar) = 1.1 g_bar, and doubles where
        # f(g_bar) = 2 g_bar.
        "log_g_bar_where_discrepancy_10pct": float(optimize.brentq(
            lambda t: _rar_model(t, lg) - t - np.log10(1.1), -14.0, -7.0)),
        "log_g_bar_where_discrepancy_2x": float(optimize.brentq(
            lambda t: _rar_model(t, lg) - t - np.log10(2.0), -14.0, -7.0)),
    }


#: Bin edges matching the published binned relation: 14 bins over 3 dex.
RAR_BIN_EDGES = np.linspace(-11.8, -8.8, 15)


def bin_rar(data: RARData) -> dict:
    """Binned relation: log of the mean accelerations in each bin, and the spread."""
    idx = np.digitize(data.log_gbar, RAR_BIN_EDGES) - 1
    out = {"centre": [], "log_gbar": [], "log_gobs": [], "sd": [], "n": []}
    for b in range(len(RAR_BIN_EDGES) - 1):
        m = idx == b
        out["centre"].append(float(0.5 * (RAR_BIN_EDGES[b] + RAR_BIN_EDGES[b + 1])))
        out["n"].append(int(m.sum()))
        if m.sum() == 0:
            out["log_gbar"].append(None), out["log_gobs"].append(None), out["sd"].append(None)
            continue
        out["log_gbar"].append(float(np.log10(np.mean(10.0 ** data.log_gbar[m]))))
        out["log_gobs"].append(float(np.log10(np.mean(10.0 ** data.log_gobs[m]))))
        out["sd"].append(float(np.std(data.log_gobs[m])))
    return out


def max_disk_population(sample: list[sparc.Galaxy], tolerance_sigma: float) -> dict:
    """The maximum-disc test for every galaxy in the sample."""
    rows = []
    for g in sample:
        md = maximum_disk(g, tolerance_sigma)
        if md is None:
            continue
        rows.append({"name": g.name, "sb_disk0": g.sb_disk0, "n_points": len(g.curve),
                     "r_last_over_r_disk": float(g.curve.r[-1] / g.r_disk_kpc), **asdict(md)})
    over = np.array([r["overshoot"] for r in rows])
    exceptions = [r for r in rows if r["overshoot"] <= 1.0]
    return {
        "tolerance_sigma": tolerance_sigma,
        "n_galaxies": len(rows),
        #: Galaxies where no mass-to-light ratio at all makes visible matter work.
        "n_no_single_upsilon_works": int((over > 1.0).sum()),
        #: Galaxies where some single ratio works -- and the ratios that takes.
        "exceptions": [{k: r[k] for k in ("name", "upsilon_disk_needed", "upsilon_disk_max",
                                          "n_points", "r_last_over_r_disk")} for r in exceptions],
        "n_explained_with_plausible_upsilon": int(sum(
            1 for r in exceptions if r["needed_is_plausible"])),
        "overshoot_median": float(np.median(over)),
        "upsilon_needed_median": float(np.median([r["upsilon_disk_needed"] for r in rows])),
        "rows": rows,
    }


# --- baryonic Tully-Fisher ---------------------------------------------------

@dataclass
class BTFRData:
    name: np.ndarray
    log_vflat: np.ndarray
    e_log_vflat: np.ndarray
    log_mbar: np.ndarray
    e_log_mbar: np.ndarray
    log_mstar: np.ndarray
    e_log_mstar: np.ndarray
    gas_fraction: np.ndarray


def btfr_data(sample: list[sparc.Galaxy]) -> BTFRData:
    """Total visible mass and flat rotation speed for every galaxy that has one.

    Mass errors combine the distance (masses scale as D**2), the 3.6um
    luminosity, the mass-to-light scatter and the HI flux calibration, each
    weighted by the share of the mass it affects. Velocity errors combine the
    published error with the inclination uncertainty.
    """
    rows = []
    for g in sample:
        if g.v_flat is None:
            continue
        m_star = UPSILON_DISK * g.luminosity_disk + UPSILON_BULGE * g.luminosity_bulge
        m_gas = rp.HELIUM_FACTOR * g.m_hi
        m_bar = m_star + m_gas
        f_star, f_gas = m_star / m_bar, m_gas / m_bar
        ln10 = np.log(10.0)
        e_distance = 2.0 * g.distance_err_mpc / g.distance_mpc / ln10
        e_light = g.luminosity_err / g.luminosity / ln10
        e_log_m = np.sqrt(e_distance**2 + (f_star * e_light) ** 2
                          + (f_star * UPSILON_SCATTER_DEX) ** 2
                          + (f_gas * np.log10(1.0 + HI_MASS_ERROR)) ** 2)
        e_log_mstar = np.sqrt(e_distance**2 + e_light**2 + UPSILON_SCATTER_DEX**2)
        inc = np.radians(g.inclination_deg)
        e_log_v = np.sqrt((g.v_flat_err / g.v_flat) ** 2
                          + (np.radians(g.inclination_err_deg) / np.tan(inc)) ** 2) / ln10
        rows.append((g.name, np.log10(g.v_flat), e_log_v, np.log10(m_bar), e_log_m,
                     np.log10(m_star), e_log_mstar, f_gas))
    cols = list(zip(*rows))
    return BTFRData(np.array(cols[0]), *(np.array(c, dtype=float) for c in cols[1:]))


def fit_line_intrinsic(x, y, ex, ey, pivot: float, scatter: str = "orthogonal") -> dict:
    """Straight line with errors on both axes plus intrinsic scatter (max. likelihood).

    The relation is allowed real scatter beyond the error bars, and the fit
    finds the slope, zero point and that scatter together. *Where* the scatter
    is assumed to lie changes the slope, so both conventions are available:

    * ``"orthogonal"`` -- scatter perpendicular to the line, the convention of
      Lelli et al. (2019) and their BayesLineFit code. Each point's orthogonal
      distance ``(y - a x - b) / sqrt(1 + a**2)`` has variance
      ``(a**2 ex**2 + ey**2) / (1 + a**2) + sigma_perp**2``.
    * ``"vertical"`` -- scatter in ``y`` only, variance
      ``ey**2 + a**2 ex**2 + sigma_y**2``. Its likelihood carries an extra
      ``ln(1 + a**2)`` per point relative to the orthogonal one, which pulls
      steep slopes down: the same data give a shallower line.
    """
    def nll(theta):
        a, b, ln_s = theta
        resid = y - a * (x - pivot) - b
        if scatter == "orthogonal":
            var = (a**2 * ex**2 + ey**2) / (1.0 + a**2) + np.exp(2 * ln_s)
            return 0.5 * np.sum(resid**2 / (1.0 + a**2) / var + np.log(var))
        var = ey**2 + a**2 * ex**2 + np.exp(2 * ln_s)
        return 0.5 * np.sum(resid**2 / var + np.log(var))

    a0, b0 = np.polyfit(x - pivot, y, 1)
    sol = optimize.minimize(nll, [a0, b0, np.log(0.05)], method="Nelder-Mead",
                            options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 40000})
    a, b, ln_s = sol.x
    return {"slope": float(a), "zero_point_at_pivot": float(b), "pivot": pivot,
            "intercept": float(b - a * pivot), "scatter_convention": scatter,
            "sigma_int_dex": float(np.exp(ln_s)),
            "rms_vertical_dex": float(np.sqrt(np.mean((y - a * (x - pivot) - b) ** 2)))}


def fit_btfr(d: BTFRData) -> dict:
    """The baryonic Tully-Fisher relation, fitted the way Lelli et al. (2019) did."""
    pivot = 2.0  # log10(100 km/s): puts the zero point in the middle of the data
    main = fit_line_intrinsic(d.log_vflat, d.log_mbar, d.e_log_vflat, d.e_log_mbar, pivot)
    rng = np.random.default_rng(RNG_SEED)
    boot = []
    for _ in range(N_BOOTSTRAP):
        k = rng.integers(0, len(d.name), len(d.name))
        boot.append(fit_line_intrinsic(d.log_vflat[k], d.log_mbar[k], d.e_log_vflat[k],
                                       d.e_log_mbar[k], pivot)["slope"])
    vertical = fit_line_intrinsic(d.log_vflat, d.log_mbar, d.e_log_vflat, d.e_log_mbar, pivot,
                                  scatter="vertical")
    # The same fit with stars only: gas-rich dwarfs fall off this one.
    stellar = fit_line_intrinsic(d.log_vflat, d.log_mstar, d.e_log_vflat, d.e_log_mstar, pivot)
    return {**main, "slope_err_bootstrap": float(np.std(boot)), "n_galaxies": len(d.name),
            "vertical_scatter_fit": vertical, "stellar_only": stellar}


# --- a units check on the published components -----------------------------

def freeman_check(sample: list[sparc.Galaxy]) -> dict:
    """Do SPARC's disc speeds mean what we think they mean?

    SPARC says its stellar-disc curves are for a mass-to-light ratio of exactly
    1. If that is right, then for bulgeless galaxies the peak of the published
    disc curve should match an ideal exponential disc (Freeman 1970) with the
    catalogued total luminosity and scale length. Real discs are only roughly
    exponential and have some thickness, so expect agreement to ~10%, with the
    median close to 1. A factor of sqrt(2) off would mean the column was really
    for Y = 0.5 -- and every conclusion here would be wrong by a factor of two.
    """
    ratios = []
    for g in sample:
        c = g.curve
        if g.luminosity_bulge > 0 or g.r_disk_kpc <= 0 or c.r[-1] < 3.0 * g.r_disk_kpc:
            continue
        r = np.linspace(0.01, 10.0 * g.r_disk_kpc, 4000)
        v_peak = rp.freeman_disk_velocity(r, g.luminosity, g.r_disk_kpc).max()
        ratios.append(c.v_disk.max() / v_peak)
    ratios = np.array(ratios)
    return {"n_galaxies": len(ratios), "median_ratio": float(np.median(ratios)),
            "ratio_16_84": [float(v) for v in np.percentile(ratios, [16, 84])]}


def analyse_population(galaxies: dict[str, sparc.Galaxy]) -> tuple[dict, dict]:
    """All the sample-wide results, plus arrays for the figures."""
    sample = select_sample(galaxies)
    data = rar_data(sample)
    rar = fit_rar(data, sample)
    btfr_d = btfr_data(sample)
    result = {
        "sample": {
            "n_catalogue": len(galaxies),
            "n_rejected_quality": sum(g.quality > MAX_QUALITY_FLAG for g in galaxies.values()),
            "n_rejected_inclination": sum(g.quality <= MAX_QUALITY_FLAG
                                          and g.inclination_deg < MIN_INCLINATION_DEG
                                          for g in galaxies.values()),
            "n_galaxies": len(sample),
            "n_points_total": int(sum(len(g.curve) for g in sample)),
        },
        "rar": rar,
        "rar_bins": bin_rar(data),
        "max_disk": max_disk_population(sample, MAX_DISK_TOLERANCE_SIGMA),
        # The same test with half the slack, to show how much the tolerance matters.
        "max_disk_1sigma": {k: v for k, v in max_disk_population(sample, 1.0).items()
                            if k != "rows"},
        "btfr": fit_btfr(btfr_d),
        "freeman_check": freeman_check(sample),
    }
    return result, {"rar": data, "btfr": btfr_d, "sample": sample}


# --------------------------------------------------------------------------
# Grading ourselves against the literature
# --------------------------------------------------------------------------

def match_published_points(ours: RARData, pub: sparc.PublishedRAR, tol: float) -> int:
    """How many published points we reproduce, one-to-one, within ``tol`` dex.

    The published file lists points without galaxy names, in a different
    order, rounded to 0.01 dex. A maximum bipartite matching pairs each of our
    points with at most one published point within the tolerance on both axes,
    so no point is counted twice.
    """
    a = np.column_stack([ours.log_gbar, ours.log_gobs])
    b = np.column_stack([pub.log_gbar, pub.log_gobs])
    neighbours = cKDTree(b).query_ball_point(a, r=tol, p=np.inf)
    rows = np.repeat(np.arange(len(a)), [len(n) for n in neighbours])
    cols = np.concatenate([np.asarray(n, dtype=int) for n in neighbours])
    graph = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(a), len(b)))
    match = maximum_bipartite_matching(graph, perm_type="column")
    return int((match >= 0).sum())


def compare_with_literature(population: dict, pop_art: dict,
                            galaxy_results: list[GalaxyResult]) -> dict:
    pub_rar, pub_bins = sparc.load_published_rar()
    pub_btfr = sparc.load_published_btfr()
    pub_md = sparc.load_published_max_disk()
    rar = population["rar"]
    ours = pop_art["rar"]

    # The published data refitted with our method: does the fitting method, or
    # the data, account for any difference from the quoted g_dagger?
    refit = fit_g_dagger(pub_rar.log_gbar, pub_rar.log_gobs, pub_rar.e_log_gbar,
                         pub_rar.e_log_gobs)

    bins = population["rar_bins"]
    bin_diffs = [o - p for o, p, n in zip(bins["log_gobs"], pub_bins.log_gobs, bins["n"])
                 if o is not None and n > 0]

    btfr_d = pop_art["btfr"]
    mb_diff = np.array([btfr_d.log_mbar[i] - pub_btfr[n]["log_mbar"]
                        for i, n in enumerate(btfr_d.name) if n in pub_btfr])
    # Their published table through our fitting code: if this returns their
    # slope, any difference in ours comes from the data, not the method.
    with_vf = [n for n, v in pub_btfr.items() if v["v_flat"]]
    btfr_refit = fit_line_intrinsic(
        np.log10([pub_btfr[n]["v_flat"] for n in with_vf]),
        np.array([pub_btfr[n]["log_mbar"] for n in with_vf]),
        np.array([pub_btfr[n]["e_v_flat"] / pub_btfr[n]["v_flat"] for n in with_vf]) / np.log(10.0),
        np.array([pub_btfr[n]["e_log_mbar"] for n in with_vf]),
        pivot=2.0,
    )

    md_rows = population["max_disk"]["rows"]
    ours_y = np.array([r["upsilon_disk_max"] for r in md_rows])
    theirs_y = np.array([pub_md[r["name"]]["upsilon_disk"] for r in md_rows])
    good = theirs_y > 0

    p_rar, p_btfr = PUBLISHED["rar"], PUBLISHED["btfr"]
    return {
        "rar": {
            "reference": p_rar["reference"],
            "g_dagger_ours": rar["g_dagger"],
            "g_dagger_published": p_rar["g_dagger"],
            "percent_diff": 100.0 * (rar["g_dagger"] / p_rar["g_dagger"] - 1.0),
            "g_dagger_published_data_refit": refit["g_dagger"],
            "rms_ours": rar["residual_rms_dex"],
            "rms_published": p_rar["rms_dex"],
            "gaussian_width_ours": rar["residual_gaussian_width_dex"],
            "gaussian_width_published": p_rar["gaussian_width_dex"],
            "n_points_ours": rar["n_points"],
            "n_points_published": p_rar["n_points"],
            "n_galaxies_ours": population["sample"]["n_galaxies"],
            "n_galaxies_published": p_rar["n_galaxies"],
            "points_matched_within_rounding": match_published_points(ours, pub_rar, 0.0051),
            "points_matched_within_0p015": match_published_points(ours, pub_rar, 0.015),
            "binned_diff_dex": [float(v) for v in bin_diffs],
            "binned_max_abs_diff_dex": float(np.max(np.abs(bin_diffs))),
            "binned_n_ours": bins["n"],
            "binned_n_published": [int(v) for v in pub_bins.n],
            "median_error_gbar_published": float(np.median(pub_rar.e_log_gbar)),
            "median_error_gobs_published": float(np.median(pub_rar.e_log_gobs)),
        },
        "btfr": {
            "reference": p_btfr["reference"],
            "slope_ours": population["btfr"]["slope"],
            "slope_err_ours": population["btfr"]["slope_err_bootstrap"],
            "slope_published": p_btfr["slope"],
            "slope_err_published": p_btfr["slope_err"],
            "slope_published_data_refit": btfr_refit["slope"],
            "sigma_int_published_data_refit": btfr_refit["sigma_int_dex"],
            "sigma_int_ours": population["btfr"]["sigma_int_dex"],
            "n_galaxies_ours": population["btfr"]["n_galaxies"],
            "n_galaxies_published_with_vflat": sum(1 for v in pub_btfr.values() if v["v_flat"]),
            "log_mbar_median_diff_dex": float(np.median(mb_diff)),
            "log_mbar_rms_diff_dex": float(np.sqrt(np.mean(mb_diff**2))),
        },
        "max_disk": {
            "reference": PUBLISHED["max_disk"]["reference"],
            "n_galaxies": int(good.sum()),
            "spearman": float(stats.spearmanr(ours_y[good], theirs_y[good])[0]),
            "median_ratio_ours_over_published": float(np.median(ours_y[good] / theirs_y[good])),
            "ratio_16_84": [float(v) for v in np.percentile(ours_y[good] / theirs_y[good], [16, 84])],
            "showcase": {r.name: {"ours": r.max_disk.get("upsilon_disk_max"),
                                  "published": pub_md.get(r.name, {}).get("upsilon_disk")}
                         for r in galaxy_results},
            "points": {"ours": ours_y.tolist(), "published": theirs_y.tolist(),
                       "sb_disk0": [r["sb_disk0"] for r in md_rows]},
        },
    }


# --------------------------------------------------------------------------
# Saving
# --------------------------------------------------------------------------

def _jsonable(obj):
    """Recursively turn results into plain JSON (non-finite floats become null)."""
    if is_dataclass(obj):
        return _jsonable(asdict(obj))
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _jsonable(obj.tolist())
    if isinstance(obj, (np.floating, float)):
        return float(obj) if np.isfinite(obj) else None
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def save_results(galaxy_results: list[GalaxyResult], population: dict, comparison: dict) -> None:
    payload = {
        "_generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "_assumptions": {"upsilon_disk": UPSILON_DISK, "upsilon_bulge": UPSILON_BULGE, "H0": H0},
        "galaxies": galaxy_results,
        "population": population,
        "comparison": comparison,
    }
    RESULTS_PATH.write_text(json.dumps(_jsonable(payload), indent=2))
