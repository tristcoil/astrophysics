"""Self-checks against cases with known answers.

Run with ``python test_physics.py``. No test framework and no network: these
check the methods, not the data.

Each step that turns raw numbers into a rotation curve is pinned here to a
synthetic case whose answer is known exactly: the spectral axis of the cubes,
the line-profile fits, the tilted rings (on a warped model disc), the
harmonic decomposition, the gas disc's own pull (Freeman's exponential disc),
the density formula (NFW and isothermal haloes), the kinematic position
angles (against pafit itself), and the reader that recovers published curves
from PostScript figures. The last check runs a whole model galaxy -- cube,
beam, noise, profile fits, rings -- through the pipeline.
"""

from __future__ import annotations

import sys

import numpy as np
from astropy.io import fits

import galmod as GM
import harmonics as H
import kinpa
import literature as L
import massmodel as MM
import things
import tiltedring as TR
from config import C_KMS, HI_REST_HZ

FAILURES: list[str] = []


def check(label, got, expected, tol, unit=""):
    """``got`` within ``tol`` (absolute) of ``expected``."""
    got = float(got)
    ok = abs(got - expected) <= tol
    if not ok:
        FAILURES.append(label)
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<58} {got:>11.5g} vs {expected:>10.5g} {unit:<7} (tol {tol:g})")


def check_true(label, condition):
    if not condition:
        FAILURES.append(label)
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label}")


def test_spectral_axis():
    print("Spectral axes of the cubes")
    # FELO-HEL: linear in frequency, labelled in optical velocity
    h = fits.Header({"NAXIS3": 5, "CTYPE3": "FELO-HEL", "CRVAL3": 375000.0, "CDELT3": -2583.0,
                     "CRPIX3": 3, "RESTFREQ": HI_REST_HZ, "VELREF": 2})
    v = things.channel_velocities(h)
    check("FELO: reference pixel gets CRVAL3", v[2], 375.0, 1e-9, "km/s")
    f = HI_REST_HZ / (1 + v / C_KMS)
    check_true("FELO: frequencies equally spaced", np.allclose(np.diff(f), np.diff(f)[0], rtol=1e-9))
    check("FELO: increment at the reference pixel", v[3] - v[2], -2.583, 2e-4, "km/s")
    # FREQ with an alternate velocity reference (AIPS): the alternate pixel gets ALTRVAL
    h = fits.Header({"NAXIS3": 72, "CTYPE3": "FREQ", "CRVAL3": 1.41583745194e9, "CDELT3": 24414.0625,
                     "CRPIX3": -18.0, "RESTFREQ": HI_REST_HZ, "VELREF": 2, "ALTRVAL": 782087.0, "ALTRPIX": 13.0})
    v = things.channel_velocities(h)
    check("FREQ: ALTRPIX gets ALTRVAL (NGC 3198)", v[12], 782.087, 1e-6, "km/s")
    check("FREQ: channel spacing near 660 km/s", abs(np.median(np.diff(v))), 5.176, 0.002, "km/s")


def test_profiles():
    print("\nLine-profile fits")
    rng = np.random.default_rng(1)
    v = np.arange(300, 500, 5.2)
    truth = (8e-3, 402.0, 11.0, 0.18)
    spec = things.gauss_hermite(v, *truth)[None, :] + rng.normal(0, 3e-4, (200, len(v)))
    p0 = things._initial_guesses(v, spec)
    g = things.fit_profiles(v, spec, p0, with_h3=False, noise=3e-4)
    h = things.fit_profiles(v, spec, np.column_stack([g["params"], np.zeros(len(spec))]), with_h3=True, noise=3e-4)
    med = np.median(h["params"], axis=0)
    check("Gauss-Hermite centre recovered (S/N 27)", med[1], truth[1], 0.3, "km/s")
    check("Gauss-Hermite dispersion recovered", med[2], truth[2], 0.3, "km/s")
    check("Gauss-Hermite h3 recovered", med[3], truth[3], 0.02)
    sd = np.std(h["params"][:, 1])
    check("quoted centre error matches the scatter (ratio)", np.median(h["errors"][:, 1]) / sd, 1.0, 0.25)
    fine = np.linspace(350, 450, 200001)
    true_peak = fine[np.argmax(things.gauss_hermite(fine, *truth))]
    check("profile peak of a skewed Gauss-Hermite", things.profile_peak(np.array([402.0]), np.array([11.0]),
                                                                       np.array([0.18]), np.array([8e-3]))[0],
          true_peak, 0.01, "km/s")
    check_true("skewness moves the peak away from v0 (h3 > 0: peak below v0)", true_peak < truth[1] - 2)


def _model_disc(n=181, pix=4.0, vsys=500.0, warp=True):
    """A disc with a mild warp; returns coordinates, velocities and the truth."""
    xs = (np.arange(n) - n // 2) * pix
    X, Y = np.meshgrid(xs, xs)
    r_true = np.linspace(0, 400, 401)
    v_true = 150 * (1 - np.exp(-r_true / 60))
    pa_true = 40 + (8 * np.clip((r_true - 200) / 200, 0, 1) if warp else 0 * r_true)
    incl_true = 60 + (6 * np.clip((r_true - 200) / 200, 0, 1) if warp else 0 * r_true)
    v = np.full(X.shape, np.nan)
    # assign each pixel through the ring whose geometry puts it at its own radius
    best = np.full(X.shape, np.inf)
    for k in range(5, 400, 2):
        rr, ct, _ = TR.disc_coordinates(X, Y, 3.0, -2.0, pa_true[k], incl_true[k])
        d = np.abs(rr - r_true[k])
        take = (d < 1.0) & (d < best)
        best[take] = d[take]
        v[take] = vsys + v_true[k] * np.sin(np.radians(incl_true[k])) * ct[take]
    ok = np.isfinite(v)
    return X[ok], Y[ok], v[ok], (r_true, v_true, pa_true, incl_true)


def test_tilted_rings():
    print("\nTilted rings on a warped model disc (noise-free)")
    x, y, v, (rt, vt, pat, it) = _model_disc()
    res = TR.analyse(x, y, v, np.ones_like(v), ring_step=10.0, beam_fwhm=10.0, beam_area_arcsec2=100.0,
                     start={"x0": 0.0, "y0": 0.0, "vsys": 495.0, "pa": 50.0, "incl": 50.0, "vrot": 100.0})
    check("centre x", res.centre[0], 3.0, 0.3, "arcsec")
    check("centre y", res.centre[1], -2.0, 0.3, "arcsec")
    check("systemic velocity", res.vsys, 500.0, 0.2, "km/s")
    sel = (res.radius > 40) & (res.radius < 360)
    dv = res.vrot[sel] - np.interp(res.radius[sel], rt, vt)
    check("rotation curve, rms error 40-360\"", np.sqrt(np.nanmean(dv ** 2)), 0.0, 0.8, "km/s")
    dpa = res.free_pa[sel] - np.interp(res.radius[sel], rt, pat)
    check("free position angles follow the warp (rms)", np.sqrt(np.nanmean(dpa ** 2)), 0.0, 0.6, "deg")
    di = res.free_incl[sel] - np.interp(res.radius[sel], rt, it)
    check("free inclinations follow the warp (rms)", np.sqrt(np.nanmean(di ** 2)), 0.0, 1.0, "deg")
    check("approaching = receding for a symmetric disc",
          np.nanmedian(np.abs(res.vrot_app - res.vrot_rec)[sel]), 0.0, 0.5, "km/s")


def test_harmonics():
    print("\nHarmonic decomposition recovers injected non-circular terms")
    x, y, v, _ = _model_disc(warp=False)
    res = TR.analyse(x, y, v, np.ones_like(v), ring_step=10.0, beam_fwhm=10.0, beam_area_arcsec2=100.0,
                     start={"x0": 3.0, "y0": -2.0, "vsys": 500.0, "pa": 40.0, "incl": 60.0, "vrot": 100.0})
    _, ct, st = TR.disc_coordinates(x, y, res.centre[0], res.centre[1], 40.0, 60.0)
    th = np.arctan2(st, ct)
    v2 = v + 4.0 * np.sin(th) + 3.0 * np.cos(2 * th) + 2.0 * np.sin(3 * th)
    h = H.decompose(x, y, v2, res)
    sel = (res.radius > 60) & (res.radius < 340)
    sin_i = np.sin(np.radians(60))
    check("s1 (radial inflow/outflow signature)", np.nanmedian(h.term("s", 1, deproject=False)[sel]), 4.0, 0.15, "km/s")
    check("c2", np.nanmedian(h.term("c", 2, deproject=False)[sel]), 3.0, 0.15, "km/s")
    check("s3", np.nanmedian(h.term("s", 3, deproject=False)[sel]), 2.0, 0.15, "km/s")
    check("A_r = sqrt(4^2+3^2+2^2) / sin i", np.nanmedian(h.a_r[sel]), np.sqrt(29) / sin_i, 0.2, "km/s")


def test_mass_models():
    print("\nThe gas disc's own pull, and density from rotation")
    r = np.linspace(0.05, 40, 800)
    sig = 50 * np.exp(-r / 4.0)
    re = np.array([1.0, 4.0, 8.6, 15.0, 25.0])
    v2 = MM.disc_rotation_curve(re, r, sig)
    ex = MM.freeman_disc_v2(re, 50, 4.0)
    check("thin disc vs Freeman (1970), worst of 5 radii", np.max(np.abs(np.sqrt(v2 / ex) - 1)), 0.0, 0.005)
    check("Freeman disc peaks at 2.15 scale lengths",
          np.linspace(0.1, 20, 4000)[np.argmax(MM.freeman_disc_v2(np.linspace(0.1, 20, 4000), 50, 4.0))] / 4.0,
          2.15, 0.01)
    check("21-cm: 1 Jy km/s at 1 Mpc is 2.356e5 Msun of HI", MM.hi_mass_msun(1.0, 1.0), 2.356e5, 1)
    check("brightness temperature of 1 Jy in a 10\" x 10\" beam", MM.brightness_temperature_k(1.0, 10, 10),
          1.222e6 / (1.420405752 ** 2 * 100), 0.01, "K")
    rr = np.logspace(-2, 1, 400)
    v_iso = MM.v_iso(rr, 0.03, 1.5)
    rho = MM.density_from_curve(rr, v_iso)
    check("isothermal sphere: central density from its curve", rho[5], 0.03 / (1 + (rr[5] / 1.5) ** 2), 1e-4,
          "Msun/pc^3")
    v_n = MM.v_nfw(rr, 10.0, 50.0)
    rho_n = MM.density_from_curve(rr, v_n)
    inner = rr < 0.05
    slope = np.polyfit(np.log10(rr[inner]), np.log10(rho_n[inner]), 1)[0]
    check("NFW: logarithmic slope near the centre", slope, -1.0, 0.03)
    check("rotation-curve slope: V ~ R^0.5 means alpha = -1",
          MM.velocity_power_slope(rr, 3 * rr ** 0.5, 1.0)["alpha"], -1.0, 1e-6)
    check("rotation-curve slope: solid body means alpha = 0",
          MM.velocity_power_slope(rr, 20 * rr, 1.0)["alpha"], 0.0, 1e-6)
    m200 = MM.m200_from_v200(200.0)
    check("M200 for V200 = 200 km/s (H0 = 70)", np.log10(m200), np.log10(200 ** 3 / (10 * 4.300917e-6 * 0.07)), 1e-9)


def test_kinematic_pa():
    print("\nKinematic position angles")
    rng = np.random.default_rng(3)
    x, y = rng.uniform(-15, 15, (2, 400))
    pa_true, inc = 125.0, np.radians(55)
    t = np.radians(pa_true)
    along = -x * np.sin(t) + y * np.cos(t)
    across = (-x * np.cos(t) - y * np.sin(t)) / np.cos(inc)
    rr = np.hypot(along, across)
    v = 180 * np.tanh(rr / 4) * np.sin(inc) * along / rr + rng.normal(0, 8, x.size)
    axis, err, vsys, _ = kinpa.fit_kinematic_pa(x, y, v, 8.0)
    check("axis of a rotating disc (mod 180)", axis, pa_true % 180, 1.5, "deg")
    check("receding side recovered", kinpa.receding_pa(x, y, v, axis), pa_true, 1.5, "deg")
    check("counter-rotation flips the receding side", kinpa.receding_pa(x, y, -v, axis), (pa_true + 180) % 360,
          1.5, "deg")
    check("angle between 10 and 350 deg", kinpa.angle_between(10, 350), 20, 1e-9, "deg")
    try:
        from pafit.fit_kinematic_pa import fit_kinematic_pa as pafit
        a2, e2, s2 = pafit(x, y, v, dvel=8.0, plot=False, quiet=True)
        check("same axis as pafit itself", axis, a2, 1e-9, "deg")
        check("same error as pafit itself", err, e2, 1e-9, "deg")
    except ImportError:
        print("  (pafit not installed: comparison skipped)")


def test_postscript_reader():
    print("\nReading plotted points back from SuperMongo PostScript")
    # a frame 1000 x 800 units, x ticks at 0/10/20 (100 + 40 units each), y ticks at 0/50/100
    ps = ["%!PS", "20 dict begin", " 80 600 translate", " 0 0 0 C",
          "100 100 M 1000 0 l 0 800 l 100 900 L 100 100 L"]
    for k, xv in enumerate([100, 500, 900]):
        ps.append(f"{xv} 100 M 0 40 l")
    for k, yv in enumerate([100, 500, 900]):
        ps.append(f"100 {yv} M 40 0 l")
    ps.append(" 150 150 150 C")
    pts = [(300, 400, 30), (700, 650, 50)]  # device x, y, half error bar
    for (px, py, e) in pts:
        ps.append(f"{px} {py} M B")
        ps.append(f"6 0 m {px + 6} {py + 6} L {px} {py + 8} L {px - 6} {py + 6} L {px - 8} {py} L "
                  f"{px - 6} {py - 6} L {px} {py - 8} L {px + 6} {py - 6} L {px + 8} {py} L")
        ps.append("CF M")
        ps.append(f"{px} {py - e} M 0 {2 * e} l")
    ps.append(" 0 0 0 C")
    plot = L.read_sm_postscript("\n".join(ps))
    xa, ya = L.calibrate_axes(plot, [0, 10, 20], [0, 50, 100])
    centres = np.array(sorted(s.centre for s in plot.filled((150, 150, 150))))
    check("first point, x", xa(centres[0, 0]), 5.0, 1e-9)
    check("first point, y", ya(centres[0, 1]), 37.5, 1e-9)
    check("second point, y", ya(centres[1, 1]), 68.75, 1e-9)
    bars = [s for s in plot.strokes((150, 150, 150))]
    check_true("error bars found as grey vertical strokes", len(bars) >= 2)


def test_pipeline_end_to_end():
    print("\nA model galaxy through the whole pipeline (labelled simulation)")
    rng = np.random.default_rng(7)
    pix = 3.0
    xs = (np.arange(121) - 60) * pix
    vel = np.arange(300.0, 700.0, 5.2)
    r = np.arange(6.0, 180.0, 6.0)
    vin = 140 * (1 - np.exp(-r / 40))
    spec = GM.ModelSpec((0.0, 0.0), 500.0, 70.0, 60.0, r, vin, 2.0 * np.exp(-r / 120), 8.0, 0.0)
    model = GM.make_cube(spec, (121, 121), xs, xs, vel, sub=3, n_z=1)
    beam = (6.0, 6.0, 0.0)
    obs = GM.observe(model, beam, pix, 3e-4, rng)
    clean = GM.observe(model, beam, pix, 0.0, rng)
    rec = GM.recover(obs, vel, 3e-4, xs, xs, spec, 6.0, mom0_noiseless=clean.sum(axis=0) * 5.2)
    outer = rec.radius >= 24
    check("recovered curve beyond 4 beams (rms)", np.sqrt(np.nanmean((rec.v_recovered - rec.v_input)[outer] ** 2)),
          0.0, 1.5, "km/s")
    check_true("beam smearing lowers the innermost point", rec.v_recovered[0] < rec.v_input[0] + 1.0)


def main() -> int:
    test_spectral_axis()
    test_profiles()
    test_tilted_rings()
    test_harmonics()
    test_mass_models()
    test_kinematic_pa()
    test_postscript_reader()
    test_pipeline_end_to_end()
    print(f"\n{'All checks passed.' if not FAILURES else f'{len(FAILURES)} FAILED: {FAILURES}'}")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
