"""Self-checks for every piece of the pipeline that a wrong sign or factor could break.

    python test_physics.py          # ~1 minute; a few checks need data/raw from fetch_data.py

Each test prints PASS/FAIL with the number it checked. Synthetic inputs (a ring image, a
Gaussian disc, a spiral, a spectral cube) are labelled as such: they test the methods
against known answers and never enter the results.
"""

from __future__ import annotations

import os
import sys
import warnings

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
from astropy.io import fits  # noqa: E402
from scipy import ndimage  # noqa: E402

import alma_images as A  # noqa: E402
import geometry as G  # noqa: E402
import kinematics as K  # noqa: E402
import literature as L  # noqa: E402
import planets as P  # noqa: E402
import profiles as PR  # noqa: E402
import spirals as S  # noqa: E402
from config import AU_CGS, IMAGE_DIR, L_SUN_CGS, PROFILE_DIR, SIGMA_SB_CGS  # noqa: E402

RESULTS: list[tuple[str, bool, str]] = []
RNG = np.random.default_rng(20261006)


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, bool(ok), detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}: {detail}")


def synthetic_image(r_ring=0.5, width=0.02, inc=40.0, pa=120.0, beam=(0.04, 0.035, 30.0), n=301, half=0.9,
                    noise=0.0, x0=0.01, y0=-0.02) -> A.Image:
    """A thin Gaussian ring seen at (inc, PA), convolved with an elliptical beam (synthetic)."""
    pix = 2 * half / (n - 1)
    x = np.linspace(half, -half, n)  # east to the left, as in the FITS images (CDELT1 < 0)
    y = np.linspace(-half, half, n)
    xx, yy = np.meshgrid(x, y)
    r, _ = G.disc_polar(xx, yy, x0, y0, inc, pa)
    sky = np.exp(-0.5 * ((r - r_ring) / width) ** 2)
    # beam convolution through an explicit kernel (pixel offsets: +column = west)
    bmaj, bmin, bpa = beam
    k = int(3 * bmaj / pix)
    ky, kx = np.mgrid[-k:k + 1, -k:k + 1]
    dx_east, dy_north = -kx * pix, ky * pix
    p = np.radians(bpa)
    u = dx_east * np.sin(p) + dy_north * np.cos(p)
    v = dx_east * np.cos(p) - dy_north * np.sin(p)
    ker = np.exp(-0.5 * ((u / (bmaj / 2.3548)) ** 2 + (v / (bmin / 2.3548)) ** 2))
    ker /= ker.sum()
    img = ndimage.convolve(sky, ker, mode="constant")
    img = img / img.max() * 1e-3 + noise * RNG.standard_normal(img.shape)
    return A.Image("synthetic", img, x, y, pix, bmaj, bmin, bpa, 2.4e11, 0.0, 0.0, fits.Header())


# ---------------------------------------------------------------------------
def test_geometry_conventions():
    x, y = G.sky_from_polar(1.0, 0.0, 0, 0, 0.0, 0.0)
    check("DSHARP azimuth convention", np.allclose([x, y], [1, 0]) and np.allclose(G.sky_from_polar(1.0, 90.0, 0, 0, 0, 0), [0, 1]),
          "PA = 0: theta = 0 points east, theta = 90 north (DSHARP II App. A)")
    r = RNG.uniform(0.05, 2, 1000)
    th = RNG.uniform(-180, 180, 1000)
    xs, ys = G.sky_from_polar(r, th, 0.03, -0.01, 57.0, 141.0)
    r2, th2 = G.disc_polar(xs, ys, 0.03, -0.01, 57.0, 141.0)
    err = max(np.abs(r2 - r).max(), np.abs(G.angle_diff_180(2 * th2, 2 * th) / 2).max())
    check("sky <-> disc round trip", err < 1e-9, f"max error {err:.1e}")
    check("theta = 90 lies along the PA", abs(G.sky_pa_of_azimuth(90.0, 141.0, 57.0) - 141.0) < 1e-9,
          f"{G.sky_pa_of_azimuth(90.0, 141.0, 57.0):.6f} deg")


def test_sampson_distance():
    x0, y0, a, inc, pa = 0.02, -0.01, 0.6, 52.0, 33.0
    th = np.linspace(-180, 180, 73)[:-1]
    off = RNG.normal(0, 0.01, th.size)
    xs, ys = G.sky_from_polar(a, th, x0, y0, inc, pa)
    # push each point off the ellipse along the sky-plane normal by a known distance
    eps = 1e-5
    xa, ya = G.sky_from_polar(a, th + eps, x0, y0, inc, pa)
    tx, ty = xa - xs, ya - ys
    nrm = np.hypot(tx, ty)
    nx, ny = ty / nrm, -tx / nrm
    d = G.ellipse_distance((x0, y0, a, np.cos(np.radians(inc)), pa), xs + off * nx, ys + off * ny)
    rel = np.max(np.abs(np.abs(d) - np.abs(off)) / np.abs(off))
    check("first-order ellipse distance", rel < 0.03, f"max relative error {rel * 100:.2f}% for 10-mas offsets")


def test_ring_recovery():
    img = synthetic_image(noise=5e-6)
    guess = G.DiscGeometry(0.0, 0.0, 30.0, 110.0)
    pts = G.extract_ridge(img, guess, "B", 0.42, 0.58, 360 / 60, 5e-6)
    f1 = G.fit_ellipse(pts, guess, 0.5, (0, 180), RNG)
    g1 = G.DiscGeometry(f1.x0, f1.y0, f1.inc, f1.pa)
    pts = G.extract_ridge(img, g1, "B", 0.42, 0.58, 360 / 60, 5e-6)
    f = G.fit_ellipse(pts, g1, 0.5, (0, 180), RNG)
    ok = abs(f.inc - 40) < 0.5 and abs(G.angle_diff_180(f.pa, 120)) < 1.0 and abs(f.a - 0.5) < 0.003
    check("ellipse fit recovers a synthetic ring", ok,
          f"i {f.inc:.2f} (40), PA {f.pa:.2f} (120), a {f.a * 1e3:.1f} mas (500), centre ({f.x0 * 1e3:.1f}, {f.y0 * 1e3:.1f}) (10, -20)")


def test_gaussian_fit():
    n, half = 301, 1.2
    x = np.linspace(half, -half, n)
    y = np.linspace(-half, half, n)
    xx, yy = np.meshgrid(x, y)
    inc, pa = 50.0, 30.0
    fw = 0.4
    beam = (0.08, 0.05, 70.0)
    c = G._cov(fw / G.FWHM, fw * np.cos(np.radians(inc)) / G.FWHM, pa) + G._cov(beam[0] / G.FWHM, beam[1] / G.FWHM, beam[2])
    ci = np.linalg.inv(c)
    dx, dy = xx - 0.05, yy + 0.03
    data = np.exp(-0.5 * (ci[0, 0] * dx * dx + 2 * ci[0, 1] * dx * dy + ci[1, 1] * dy * dy))
    img = A.Image("synthetic", data, x, y, x[0] - x[1], *beam, 2.4e11, 0, 0, fits.Header())
    g = G.fit_gaussian(img, 0.0, 0.0, 0.9)
    check("Gaussian fit deconvolves the beam", abs(g.inc - inc) < 0.3 and abs(G.angle_diff_180(g.pa, pa)) < 0.3,
          f"i {g.inc:.2f} (50), PA {g.pa:.2f} (30), FWHM {g.fwhm_maj:.3f} (0.400)")


def test_profile_and_widths():
    img = synthetic_image(r_ring=0.5, width=0.03, noise=0.0)
    prof = PR.radial_profile(img, G.DiscGeometry(0.01, -0.02, 40.0, 120.0), 0.005, 0.85)
    rpk = prof.r[np.nanargmax(prof.intensity)]
    check("profile of a synthetic ring peaks at its radius", abs(rpk - 0.5) <= 0.005, f"{rpk:.3f}\" (0.500)")
    # widths on a piecewise-linear profile with known crossings
    r = np.arange(0.0, 5.0, 0.01) + 0.005
    i = np.interp(r, [0, 1, 2, 3, 4, 5], [0.5, 1.0, 0.0, 1.0, 0.0, 0.0])
    p = PR.Profile(r, i, np.full(r.size, 0.01), np.zeros(r.size), np.ones(r.size), np.ones(r.size), 0.01)
    m = PR.pair_measures(p, [("B1", "B", 1.0), ("D2", "D", 2.0), ("B3", "B", 3.0), ("D4", "D", 4.0)], 0.001)[0]
    ok = abs(m.gap_width - 1.0) < 0.02 and abs(m.ring_width - 1.0) < 0.02 and abs(m.delta - 0.4) < 0.01
    check("DSHARP II width definitions", ok, f"gap {m.gap_width:.3f} (1), ring {m.ring_width:.3f} (1), "
          f"Delta {m.delta:.3f} (0.4), depth {m.depth:.3f} (0)")


def test_zhang_relation():
    stars = L.dsharp1_stars()
    devs = []
    for p in [p for p in L.dsharp7_planets() if p.from_fit]:
        for (model, alpha), m_pub in p.masses.items():
            stu = dict(zip(("DSD1", "1mm", "DSD2"), p.st_used))[model]
            m, _, _ = P.planet_mass_mjup(p.width, p.h_over_r, alpha, stu, model, p.mstar)
            devs.append(abs(np.log10(m / m_pub)))
        hr = P.aspect_ratio(p.r_gap_au, stars[p.key].lstar, stars[p.key].mstar)
        devs_hr = abs(hr - p.h_over_r)
        if devs_hr > 0.006:
            check(f"h/r of {p.key} {p.r_gap_au:.0f} au", False, f"{hr:.3f} vs {p.h_over_r}")
    devs = np.array(devs)
    n_ok = int(np.sum(devs < 0.05))
    check("DSHARP VII Table 3 from its own inputs", n_ok >= len(devs) - 5,
          f"{n_ok} of {len(devs)} masses within 0.05 dex (the rest: Elias 24, whose 3 entries imply Delta = 0.334 "
          "not the printed 0.32, and masses rounded to 0.01 M_J)")
    t1 = P.midplane_temperature(1.0, 1.0)
    t_ref = (0.02 * L_SUN_CGS / (8 * np.pi * AU_CGS**2 * SIGMA_SB_CGS)) ** 0.25
    hr = P.aspect_ratio(100.0, 1.0, 1.0)
    check("passive-disc temperature and thickness", abs(t1 - t_ref) < 1e-9 and abs(hr - 0.07) < 0.002,
          f"T(1 au) = {t1:.1f} K; h/r(100 au, Sun-like star) = {hr:.4f} (DSHARP II: phi = 0.02 'corresponds to h/r ~ 0.07')")


def test_brightness_temperatures():
    t = np.array([5.0, 20.0, 80.0])
    nu = 2.4e11
    back = A.brightness_temperature(A.planck(nu, t), nu)
    ok1 = np.allclose(back, t, rtol=1e-10)
    f = PROFILE_DIR / "AS209.profile.txt"
    if not f.exists():
        check("Planck inversion", ok1, "round trip (DSHARP profile file not downloaded)")
        return
    rel = np.loadtxt(f)
    with fits.open(IMAGE_DIR / "AS209_continuum.fits") as h:
        hdr = h[0].header
    img = A.Image("AS209", np.zeros((2, 2)), np.zeros(2), np.zeros(2), 1.0, hdr["BMAJ"] * 3600, hdr["BMIN"] * 3600,
                  hdr["BPA"], float(hdr["CRVAL3"]), 0, 0, hdr)
    i_cgs = A.jy_per_beam_to_cgs(rel[:, 2], img)
    tb_rj = A.brightness_temperature(i_cgs, img.freq_hz, rayleigh_jeans=True)
    tb_pl = A.brightness_temperature(i_cgs, img.freq_hz)
    ok = rel[:, 2] > 0
    e1 = np.max(np.abs(tb_rj[ok] / rel[ok, 4] - 1))
    e2 = np.nanmax(np.abs(tb_pl[ok] / rel[ok, 6] - 1))
    check("brightness temperatures match DSHARP's profile files", ok1 and e1 < 2e-3 and e2 < 2e-3,
          f"RJ max deviation {e1 * 100:.3f}%, Planck {e2 * 100:.3f}% (beam, frequency from the FITS header)")


def test_kinematic_model():
    m = K.DiscModel(0.0, 0.0, 45.0, 300.0, 2.0, 5.0, 0.0, 1.0, 0.0, 4.0, distance=100.0, near=1)
    pa = np.radians(300.0)
    v, _ = m.velocity(np.array([1.0 * np.sin(pa)]), np.array([1.0 * np.cos(pa)]))
    vk = np.sqrt(6.67430e-8 * 2.0 * 1.98840987e33 / (100 * AU_CGS)) / 1e5
    check("thin Keplerian disc on its red major axis", abs(v[0] - (5.0 + vk * np.sin(np.radians(45)))) < 1e-9,
          f"{v[0]:.4f} km/s = v_sys + v_K sin i")
    for near in (1, -1):
        ms = K.DiscModel(0.01, -0.02, 46.7, 312.8, 2.0, 5.76, 0.28, 1.2, 0.006, 3.8, distance=101.0, near=near)
        r = RNG.uniform(0.3, 3.0, 500)
        ph = RNG.uniform(-np.pi, np.pi, 500)
        xs, ys = ms.sky(r, ph, ms.height(r))
        r2, ph2, _ = ms.disc_coords(xs, ys, n_iter=30)
        err = max(np.abs(r2 - r).max(), np.abs(np.angle(np.exp(1j * (ph2 - ph)))).max())
        check(f"emission-surface round trip (near = {near:+d})", err < 1e-6, f"max error {err:.1e}")
    a = K.DiscModel(0, 0, 46.7, 312.8, 2.0, 5.76, 0.28, 1.2, 0, 4, near=1)
    b = K.DiscModel(0, 0, 46.7, 312.8, 2.0, 5.76, 0.28, 1.2, 0, 4, near=-1)
    check("the near-side flag flips the near side", abs(((a.near_side_pa - b.near_side_pa) % 360) - 180) < 1e-9,
          f"{a.near_side_pa:.1f} vs {b.near_side_pa:.1f} deg")


def test_line_centres():
    vel = np.arange(-5, 5.01, 0.32)
    ny, nx = 20, 30
    true = RNG.uniform(-2, 2, (ny, nx))
    data = np.exp(-0.5 * ((vel[:, None, None] - true[None]) / 0.35) ** 2)
    rms = 0.002
    data = data + rms * RNG.standard_normal(data.shape)
    cube = A.Cube("synthetic", data.astype(np.float32), vel, np.arange(nx), np.arange(ny), 1.0, np.ones(vel.size),
                  np.ones(vel.size), np.zeros(vel.size), 2.3e11)
    lm = K.line_centres(cube, rms, -5, 5)
    err = np.nanmax(np.abs(lm.v0 - true))
    check("line centres of synthetic Gaussian lines", err < 0.03, f"max error {err * 1e3:.0f} m/s "
          "(quadratic interpolation is biased by < 0.1 channel for lines this narrow)")


def test_spirals():
    # a synthetic two-armed log spiral, r = R0 exp(b theta), b = 0.25 (pitch 14.0 deg), face-on
    n, half = 401, 1.2
    x = np.linspace(half, -half, n)
    y = np.linspace(-half, half, n)
    xx, yy = np.meshgrid(x, y)
    geom = G.DiscGeometry(0.0, 0.0, 30.0, 150.0)
    r, th = geom.polar(xx, yy)
    b = 0.25
    phase = 2 * (np.radians(th) - np.log(np.maximum(r, 1e-3) / 0.3) / b)
    data = (1 + np.cos(phase)) * np.exp(-r / 0.5) * (r > 0.1)
    img = A.Image("synthetic", data, x, y, x[0] - x[1], 0.03, 0.03, 0.0, 2.4e11, 0, 0, fits.Header())
    pol = S.polar_image(img, geom, 1.0, 0.0075, 1.0)
    p, amp = S.log_spiral_spectrum(pol, 0.2, 0.9, m=2)
    _, bb, mu = S.pitch_from_spectrum(p, amp)
    check("log-spiral Fourier pitch of a synthetic spiral", abs(mu - 14.04) < 0.6 and bb > 0,
          f"{mu:.2f} deg (14.04), b {bb:+.3f} (+0.25)")
    thd = np.linspace(-150, 150, 40)
    arm = S.ArmPoints(thd, 0.3 * np.exp(0.25 * np.radians(thd)) * (1 + 0.003 * RNG.standard_normal(40)), np.full(40, 0.003))
    f = S.fit_log_spiral(arm, np.full(40, 0.003))
    check("log-spiral fit of traced points", abs(f.b - 0.25) < 0.005, f"b {f.b:.4f} (0.25), pitch {f.pitch:.2f} deg")


def test_literature_tables():
    st, im, fe, ge, pl = L.dsharp1_stars(), L.dsharp1_images(), L.dsharp2_features(), L.dsharp2_geometry(), L.dsharp7_planets()
    counts = (len(st), len(im), len(fe), len(ge), sum(p.from_fit for p in pl), sum(not p.from_fit for p in pl))
    check("table sizes", counts == (20, 20, 107, 18, 12, 7), f"{counts} (20, 20, 107, 18, 12, 7)")
    f = {x.name: x for x in fe}
    ok = (f["AS209:D61"].r_mas.value == 503 and f["AS209:D61"].width_au.value == 15.5 and
          f["HD163296:B100"].inc.value == 46.63 and abs(st["Elias27"].distance.err_hi - 19) < 1e-9 and
          abs(st["GWLup"].log_lstar.value + 0.48) < 1e-9)
    check("spot values (multi-line rows, asymmetric errors, signs outside math mode)", ok,
          "AS 209 D61 503 mas / 15.5 au; HD 163296 B100 i = 46.63; Elias 27 d +19; GW Lup log L -0.48")


def test_distances():
    from distances import gaia_distances

    try:
        g = gaia_distances()
    except FileNotFoundError:
        check("distances", True, "skipped (Gaia table not downloaded)")
        return
    ok = all(s["zero_point"] < 0 and s["parallax_corr"] > s["parallax"] for k, s in g.items() if "parallax" in s)
    check("zero point makes parallaxes larger", ok, "all corrected parallaxes exceed the catalogue values")
    check("Elias 27 takes its cloud's distance", g["Elias27"]["distance_source"].startswith("mean of"),
          f"{g['Elias27']['distance']:.1f} pc ({g['Elias27']['distance_source'][:40]}...)")


def main() -> int:
    tests = [test_geometry_conventions, test_sampson_distance, test_ring_recovery, test_gaussian_fit,
             test_profile_and_widths, test_zhang_relation, test_brightness_temperatures, test_kinematic_model,
             test_line_centres, test_spirals, test_literature_tables, test_distances]
    for t in tests:
        print(t.__name__)
        t()
    n_fail = sum(not ok for _, ok, _ in RESULTS)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
