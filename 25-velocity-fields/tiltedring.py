"""Tilted-ring fits: from a velocity field to a rotation curve.

A thin disc in circular rotation, seen at inclination i, shows a line-of-sight
velocity at each point

    v = V_sys + V_rot(R) sin(i) cos(theta),

theta being the angle around the disc from the receding half of the major
axis. Lines of equal velocity fan out from the centre: the "spider diagram".
Real discs warp, so the model is a set of concentric rings, each with its own
centre, systemic velocity, position angle and inclination (Rogstad, Lockhart &
Wright 1974; Begeman 1987, whose GIPSY task ROTCUR de Blok et al. used). We
reimplement ROTCUR's scheme with scipy's least squares and follow de Blok et
al.'s (2008) sequence:

1. all parameters free, to find the dynamical centre and systemic velocity;
2. centre and V_sys fixed, position angle and inclination free per ring --
   this is where warps appear;
3. those profiles smoothed (box-car), then held fixed;
4. V_rot alone, for both sides of the disc together and for the approaching
   and receding halves separately.

Pixels near the minor axis carry little rotation signal, so each is weighted
by |cos theta|. Coordinates are in arcseconds on the sky, x increasing to the
west (with the pixel columns) and y to the north; position angles run from
north through east to the receding major axis.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.optimize import least_squares

from config import GEOMETRY_SMOOTH_RINGS, MIN_BEAMS_PER_RING, WEIGHT_POWER

PARAMS = ("x0", "y0", "vsys", "pa", "incl", "vrot")


def disc_coordinates(x, y, x0, y0, pa_deg, incl_deg):
    """Radius in the disc plane and the angle theta from the receding major axis."""
    pa, inc = np.radians(pa_deg), np.radians(incl_deg)
    dx, dy = x - x0, y - y0
    along = -dx * np.sin(pa) + dy * np.cos(pa)
    across = (-dx * np.cos(pa) - dy * np.sin(pa)) / np.cos(inc)
    r = np.hypot(along, across)
    with np.errstate(invalid="ignore", divide="ignore"):
        cos_t = np.where(r > 0, along / r, 1.0)
        sin_t = np.where(r > 0, across / r, 0.0)
    return r, cos_t, sin_t


def los_velocity(x, y, p: dict) -> np.ndarray:
    _, cos_t, _ = disc_coordinates(x, y, p["x0"], p["y0"], p["pa"], p["incl"])
    return p["vsys"] + p["vrot"] * np.sin(np.radians(p["incl"])) * cos_t


@dataclass
class Ring:
    radius: float  # arcsec
    width: float
    p: dict  # x0, y0, vsys, pa, incl, vrot
    err: dict = field(default_factory=dict)
    npix: int = 0
    rms: float = np.nan  # rms of line-of-sight residuals in the ring, km/s
    side: str = "both"
    ok: bool = True


def _select(x, y, v, p, radius, width, side):
    r, cos_t, _ = disc_coordinates(x, y, p["x0"], p["y0"], p["pa"], p["incl"])
    sel = np.abs(r - radius) <= 0.5 * width
    if side == "receding":
        sel &= cos_t > 0
    elif side == "approaching":
        sel &= cos_t < 0
    return sel


def fit_ring(x, y, v, radius, width, start: dict, free: tuple[str, ...], *, side="both",
             min_pixels=20, max_iter=8, weight_power=WEIGHT_POWER) -> Ring:
    """Fit one ring, re-selecting its pixels as the geometry settles."""
    p = dict(start)
    sel = _select(x, y, v, p, radius, width, side)
    ring = Ring(radius, width, dict(p), side=side)
    for _ in range(max_iter):
        if sel.sum() < min_pixels:
            ring.ok = False
            ring.npix = int(sel.sum())
            return ring
        xs, ys, vs = x[sel], y[sel], v[sel]

        def resid(q):
            pp = dict(p)
            pp.update(zip(free, q))
            _, cos_t, _ = disc_coordinates(xs, ys, pp["x0"], pp["y0"], pp["pa"], pp["incl"])
            model = pp["vsys"] + pp["vrot"] * np.sin(np.radians(pp["incl"])) * cos_t
            return np.abs(cos_t) ** (weight_power / 2) * (vs - model)

        lo = {"x0": p["x0"] - 60, "y0": p["y0"] - 60, "vsys": p["vsys"] - 100, "pa": -720,
              "incl": 5.0, "vrot": 0.0}
        hi = {"x0": p["x0"] + 60, "y0": p["y0"] + 60, "vsys": p["vsys"] + 100, "pa": 720,
              "incl": 88.0, "vrot": 600.0}
        q0 = np.array([np.clip(p[k], lo[k] + 1e-6, hi[k] - 1e-6) for k in free])
        sol = least_squares(resid, q0, bounds=([lo[k] for k in free], [hi[k] for k in free]),
                            x_scale="jac")
        p.update(zip(free, sol.x))
        new_sel = _select(x, y, v, p, radius, width, side)
        converged = np.array_equal(new_sel, sel)
        sel = new_sel
        if converged:
            break
    # formal errors from the Jacobian, scaled by the residual scatter
    dof = max(len(sol.fun) - len(free), 1)
    s2 = float(np.sum(sol.fun ** 2) / dof)
    try:
        cov = np.linalg.inv(sol.jac.T @ sol.jac) * s2
        err = dict(zip(free, np.sqrt(np.clip(np.diag(cov), 0, None))))
    except np.linalg.LinAlgError:
        err = {k: np.nan for k in free}
    p["pa"] = p["pa"] % 360.0
    ring.p, ring.err, ring.npix = p, err, int(sel.sum())
    ring.rms = float(np.sqrt(np.mean((v[sel] - los_velocity(x[sel], y[sel], p)) ** 2)))
    return ring


def fit_vrot_fixed(x, y, v, ring_p: dict, radius, width, side="both"):
    """V_rot (and its spread) for a ring of fixed geometry: a weighted mean.

    With geometry fixed the model is linear in V_rot, so the least-squares
    answer is closed-form: V = sum w (v - V_sys) s / sum w s^2, with
    s = sin(i) cos(theta) and w = |cos theta|. Returns (vrot, formal error,
    rms line-of-sight residual, npix).
    """
    sel = _select(x, y, v, ring_p, radius, width, side)
    if sel.sum() < 3:
        return np.nan, np.nan, np.nan, int(sel.sum())
    _, cos_t, _ = disc_coordinates(x[sel], y[sel], ring_p["x0"], ring_p["y0"], ring_p["pa"],
                                   ring_p["incl"])
    s = np.sin(np.radians(ring_p["incl"])) * cos_t
    w = np.abs(cos_t) ** WEIGHT_POWER
    dv = v[sel] - ring_p["vsys"]
    vrot = float(np.sum(w * dv * s) / np.sum(w * s * s))
    res = dv - vrot * s
    # ROTCUR's "dispersion along the ring": sqrt(chi^2 / (N - 1)) with the
    # |cos theta| weights left unnormalized. It is ~0.8x the plain rms
    # (<|cos theta|> = 2/pi); with it our error bars match de Blok et al.'s to
    # ~20%, against ~60% for the plain deprojected rms.
    dispersion = float(np.sqrt(np.sum(w * res ** 2) / max(sel.sum() - 1, 1)))
    err = float(np.sqrt(np.sum(w * res ** 2) / max(sel.sum() - 1, 1) / np.sum(w * s * s)))
    return vrot, err, dispersion, int(sel.sum())


# --------------------------------------------------------------------------
# The full sequence
# --------------------------------------------------------------------------

@dataclass
class TiltedRingResult:
    radius: np.ndarray  # arcsec
    width: float
    centre: tuple[float, float]  # arcsec (x, y) in the map frame
    vsys: float
    free_pa: np.ndarray  # stage 2, per ring
    free_pa_err: np.ndarray
    free_incl: np.ndarray
    free_incl_err: np.ndarray
    pa: np.ndarray  # adopted (smoothed)
    incl: np.ndarray
    vrot: np.ndarray  # stage 4, both sides
    vrot_formal_err: np.ndarray
    vrot_app: np.ndarray
    vrot_rec: np.ndarray
    ring_rms: np.ndarray  # ROTCUR-style dispersion along each ring, km/s (line of sight)
    npix: np.ndarray
    vrot_err: np.ndarray  # de Blok et al.'s definition, see below
    stage1: list = field(default_factory=list)
    notes: list = field(default_factory=list)

    def ring_params(self, k: int) -> dict:
        return {"x0": self.centre[0], "y0": self.centre[1], "vsys": self.vsys,
                "pa": self.pa[k], "incl": self.incl[k], "vrot": self.vrot[k]}


def initial_geometry(x, y, v, intensity, centre=None) -> dict:
    """Rough starting values from the maps themselves.

    Centre: the given one (the pointing centre of the observation, i.e. the
    galaxy's catalogued optical centre), else the intensity-weighted centroid
    of the gas -- which a lopsided gas disc can put ~50" off. V_sys: midway
    between the 5th and 95th percentile velocities. Position angle: the
    direction in which velocity rises fastest. Inclination: from the axis
    ratio of the gas distribution's second moments, cos i = b/a.
    """
    w = np.clip(intensity, 0, None)
    if centre is None:
        x0, y0 = np.sum(w * x) / np.sum(w), np.sum(w * y) / np.sum(w)
    else:
        x0, y0 = centre
    lo, hi = np.percentile(v, [5, 95])
    vsys = 0.5 * (lo + hi)
    gx, gy = np.sum((v - vsys) * (x - x0)), np.sum((v - vsys) * (y - y0))
    pa = np.degrees(np.arctan2(-gx, gy)) % 360
    cxx = np.sum(w * (x - x0) ** 2) / np.sum(w)
    cyy = np.sum(w * (y - y0) ** 2) / np.sum(w)
    cxy = np.sum(w * (x - x0) * (y - y0)) / np.sum(w)
    ev = np.linalg.eigvalsh(np.array([[cxx, cxy], [cxy, cyy]]))
    incl = np.degrees(np.arccos(np.clip(np.sqrt(ev[0] / ev[1]), 0.1, 0.99)))
    return {"x0": x0, "y0": y0, "vsys": vsys, "pa": pa, "incl": incl, "vrot": 0.5 * (hi - lo) / 2}


def boxcar(values: np.ndarray, width: int) -> np.ndarray:
    """Running mean ignoring NaNs, window shrinking at the ends."""
    out = np.full_like(values, np.nan, dtype=float)
    h = width // 2
    for k in range(len(values)):
        seg = values[max(0, k - h): k + h + 1]
        seg = seg[np.isfinite(seg)]
        if len(seg):
            out[k] = seg.mean()
    return out


def _unwrap_pa(pa: np.ndarray, ref: float) -> np.ndarray:
    return ref + (pa - ref + 180) % 360 - 180


def _min_pixels_rule(x, y, ring_step, incl_deg, beam_area_arcsec2):
    pixel_area = _pixel_area(x, y)

    def min_pixels(r: float) -> int:
        """A ring needs 4 beams' worth of pixels, or half of its full area if smaller."""
        full = 2 * np.pi * r * ring_step * np.cos(np.radians(incl_deg)) / pixel_area
        return int(max(10, min(MIN_BEAMS_PER_RING * beam_area_arcsec2 / pixel_area, 0.5 * full)))

    return min_pixels


def stage_one(x, y, v, start, *, ring_step, beam_fwhm, min_pixels, free=PARAMS, self_referencing=False,
              max_rings=400) -> list[Ring]:
    """Rings fitted outwards one by one, each from the running median of those before.

    One wild ring cannot drag the next ones along. A leash keeps wild rings
    out: the strict one keeps a ring only if its centre stays within two
    beams of the starting guess (and its V_sys within 30 km/s); the
    self-referencing one allows six beams until three rings agree, then two
    beams around their running median, for when the starting guess is poor
    (a lopsided gas disc can put its centroid ~50" from its centre).
    """
    accepted, misses = [], 0
    for r in ring_step * np.arange(1, max_rings + 1):
        if accepted:
            recent = accepted[-5:]
            p = {k: float(np.median([rg.p[k] for rg in accepted])) for k in ("x0", "y0", "vsys")}
            p["pa"] = float(np.median([_unwrap_pa(rg.p["pa"], recent[-1].p["pa"]) for rg in recent]))
            p["incl"] = float(np.median([rg.p["incl"] for rg in recent]))
            p["vrot"] = recent[-1].p["vrot"]
        else:
            p = dict(start)
        ring = fit_ring(x, y, v, r, ring_step, p, free, min_pixels=min_pixels(r))
        if self_referencing and len(accepted) >= 3:
            ref, leash = p, 2 * beam_fwhm
        elif self_referencing:
            ref, leash = start, 6 * beam_fwhm
        else:
            ref, leash = start, 2 * beam_fwhm
        sane = ring.ok and np.hypot(ring.p["x0"] - ref["x0"], ring.p["y0"] - ref["y0"]) < leash \
            and abs(ring.p["vsys"] - ref["vsys"]) < 30 and 6 < ring.p["incl"] < 87
        if not sane:
            misses += bool(accepted) and not ring.ok
            if misses >= 3:
                break
            continue
        misses = 0
        accepted.append(ring)
    return accepted


def _centre_rings(rings: list[Ring], beam_fwhm, centre_rings_max=None) -> np.ndarray:
    """Rings between 2 beams and half the disc: resolved, and with many pixels."""
    r1 = np.array([rg.radius for rg in rings])
    use = (r1 >= 2 * beam_fwhm) & (r1 <= (centre_rings_max or 0.5 * r1.max()))
    if use.sum() < 3:
        use = r1 <= (centre_rings_max or 0.5 * r1.max())
    return use


def analyse(x, y, v, intensity, *, ring_step, beam_fwhm, beam_area_arcsec2, start=None,
            centre_rings_max=None, inner_fixed_beams=2.0, inner_fixed_arcsec=None,
            max_rings=400) -> TiltedRingResult:
    """Run the whole de Blok et al. (2008) sequence on one velocity field.

    ``x, y, v, intensity`` are 1-D arrays of the usable pixels (arcsec, km/s).
    """
    start = start or initial_geometry(x, y, v, intensity)
    notes = []
    min_pixels = _min_pixels_rule(x, y, ring_step, start["incl"], beam_area_arcsec2)

    # Stage 1: every parameter free, rings fitted outwards; the rings find the
    # centre themselves (the self-referencing leash of ``stage_one``).
    stage1 = stage_one(x, y, v, start, ring_step=ring_step, beam_fwhm=beam_fwhm, min_pixels=min_pixels,
                       self_referencing=True, max_rings=max_rings)
    if not stage1:
        raise RuntimeError("no ring could be fitted")
    r1 = np.array([rg.radius for rg in stage1])
    r_last = r1.max()
    # Centre and V_sys: median over rings between 2 beams and half the disc,
    # where the centre is both resolved and constrained by many pixels.
    use = _centre_rings(stage1, beam_fwhm, centre_rings_max)
    x0 = float(np.median([rg.p["x0"] for rg, u in zip(stage1, use) if u]))
    y0 = float(np.median([rg.p["y0"] for rg, u in zip(stage1, use) if u]))
    vsys = float(np.median([rg.p["vsys"] for rg, u in zip(stage1, use) if u]))
    pa_ref = float(np.median([rg.p["pa"] for rg, u in zip(stage1, use) if u]))

    # Stage 2: centre and V_sys fixed; PA, inclination and V_rot free, on the
    # regular grid of radii out to the last ring stage 1 could fit. Each ring
    # starts from the nearest stage-1 ring.
    radius = ring_step * np.arange(1, int(round(r_last / ring_step)) + 1)
    free_pa, free_pa_err, free_incl, free_incl_err = (np.full(len(radius), np.nan) for _ in range(4))
    p = {**start, "x0": x0, "y0": y0, "vsys": vsys, "pa": pa_ref}
    for k, r in enumerate(radius):
        rg1 = stage1[int(np.argmin(np.abs(r1 - r)))]
        p_start = {**p, "pa": _unwrap_pa(rg1.p["pa"], pa_ref), "incl": rg1.p["incl"],
                   "vrot": rg1.p["vrot"]}
        ring = fit_ring(x, y, v, r, ring_step, p_start, ("pa", "incl", "vrot"),
                        min_pixels=min_pixels(r))
        if ring.ok and 6 < ring.p["incl"] < 87:
            free_pa[k] = _unwrap_pa(ring.p["pa"], pa_ref)
            free_incl[k] = ring.p["incl"]
            free_pa_err[k] = ring.err.get("pa", np.nan)
            free_incl_err[k] = ring.err.get("incl", np.nan)

    # Stage 3: smooth; inside `inner_fixed_beams` beams the geometry is not
    # resolved, so those rings take the value of the first resolved ring.
    pa_s = boxcar(free_pa, GEOMETRY_SMOOTH_RINGS)
    in_s = boxcar(free_incl, GEOMETRY_SMOOTH_RINGS)
    r_fix = inner_fixed_arcsec if inner_fixed_arcsec else inner_fixed_beams * beam_fwhm
    inner = radius < r_fix
    if inner.any() and (~inner).any():
        k0 = np.argmax(~inner)
        pa_s[inner], in_s[inner] = pa_s[k0], in_s[k0]
        notes.append(f"PA and inclination inside R = {r_fix:.0f}\" fixed to the ring at {radius[k0]:.0f}\"")
    for arr in (pa_s, in_s):  # fill any gaps by interpolation
        bad = ~np.isfinite(arr)
        if bad.any() and (~bad).any():
            arr[bad] = np.interp(radius[bad], radius[~bad], arr[~bad])

    # Stage 4: V_rot only, both sides and each side.
    n = len(radius)
    vrot, vrot_ferr, rms, npix, v_app, v_rec = (np.full(n, np.nan) for _ in range(6))
    for k in range(n):
        geo = {"x0": x0, "y0": y0, "vsys": vsys, "pa": pa_s[k], "incl": in_s[k]}
        vrot[k], vrot_ferr[k], rms[k], npix[k] = fit_vrot_fixed(x, y, v, geo, radius[k], ring_step)
        v_app[k] = fit_vrot_fixed(x, y, v, geo, radius[k], ring_step, "approaching")[0]
        v_rec[k] = fit_vrot_fixed(x, y, v, geo, radius[k], ring_step, "receding")[0]

    # de Blok et al.'s uncertainty: the dispersion along the ring (ROTCUR's,
    # see fit_vrot_fixed) and a quarter of the approaching-receding
    # difference, added in quadrature.
    asym = np.abs(v_app - v_rec) / 4
    vrot_err = np.sqrt(rms ** 2 + np.nan_to_num(asym) ** 2)

    return TiltedRingResult(radius=radius, width=ring_step, centre=(x0, y0), vsys=vsys,
                            free_pa=free_pa, free_pa_err=free_pa_err, free_incl=free_incl,
                            free_incl_err=free_incl_err, pa=pa_s, incl=in_s, vrot=vrot,
                            vrot_formal_err=vrot_ferr, vrot_app=v_app, vrot_rec=v_rec,
                            ring_rms=rms, npix=npix, vrot_err=vrot_err, stage1=stage1, notes=notes)


def _pixel_area(x, y) -> float:
    """Area of one pixel (arcsec^2), from the spacing of the coordinate grid."""
    ux = np.unique(np.round(x, 4))
    step = np.median(np.diff(ux)) if len(ux) > 1 else 1.0
    return float(step ** 2)


# --------------------------------------------------------------------------
# Model velocity fields
# --------------------------------------------------------------------------

def model_field(x, y, res: TiltedRingResult, *, vrot=None) -> tuple[np.ndarray, np.ndarray]:
    """Model line-of-sight velocity at each (x, y), and the ring index used.

    Each pixel is assigned to the ring whose own geometry places it closest to
    that ring's radius (a warped disc has no single radius per pixel); where
    the geometry changes between neighbouring rings, a pixel can fall between
    them, and then goes to the nearest within one ring width. Pixels beyond
    the last ring stay NaN.
    """
    vrot = res.vrot if vrot is None else vrot
    best = np.full(x.shape, np.inf)
    out = np.full(x.shape, np.nan)
    idx = np.full(x.shape, -1)
    for k, r in enumerate(res.radius):
        if not np.isfinite(vrot[k]):
            continue
        rr, cos_t, _ = disc_coordinates(x, y, res.centre[0], res.centre[1], res.pa[k], res.incl[k])
        d = np.abs(rr - r) / res.width
        take = (d <= 1.0) & (d < best)
        best[take] = d[take]
        out[take] = res.vsys + vrot[k] * np.sin(np.radians(res.incl[k])) * cos_t[take]
        idx[take] = k
    # the central half-ring: interpolate linearly from V = 0 at R = 0
    rr, cos_t, _ = disc_coordinates(x, y, res.centre[0], res.centre[1], res.pa[0], res.incl[0])
    inside = (rr < res.radius[0] - 0.5 * res.width) & ~np.isfinite(out)
    out[inside] = res.vsys + vrot[0] * rr[inside] / res.radius[0] * np.sin(np.radians(res.incl[0])) * cos_t[inside]
    idx[inside] = 0
    return out, idx


def with_vrot(res: TiltedRingResult, vrot: np.ndarray) -> TiltedRingResult:
    return replace(res, vrot=np.asarray(vrot, dtype=float))
