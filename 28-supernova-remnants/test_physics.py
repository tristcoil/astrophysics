"""Check every method against an answer known independently.

    python test_physics.py

The inputs here are synthetic -- made-up star fields, textures and vectors
with a known answer -- and exist only to test the code; no figure or result
uses them. Each check prints PASS or FAIL with the numbers compared, and the
script exits non-zero if any check fails. Runs in about a minute, offline.
"""

from __future__ import annotations

import sys
from datetime import datetime

import numpy as np
from astropy.io import fits
from skimage.feature import match_template

import astrometry as A
import casa_xray as X
import config
import kicks as K
import motions as M
import remnants as R

RESULTS: list[tuple[str, bool, str]] = []
rng = np.random.default_rng(12345)


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, bool(ok), detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else ""))


# --------------------------------------------------------------------------

def test_dates_and_projection():
    print("Dates and the tangent plane")
    from astropy.time import Time

    jy = A.julian_year(datetime(1951, 11, 7, 10, 25))
    ref = Time("1951-11-07T10:25:00", scale="tt").jyear
    check("Julian year of 1951-11-07 10:25 matches astropy", abs(jy - ref) < 1e-9, f"{jy:.6f}")
    h = fits.Header({"DATE-OBS": "1951-11-07T09:60:00", "EXPOSURE": 50.0})
    check("DSS date with '60' minutes = 10:00 + half of 50 min",
          abs(A.dss_epoch(h) - A.julian_year(datetime(1951, 11, 7, 10, 25))) < 1e-9)
    ra0, dec0 = 83.633, 22.0145
    ra = ra0 + rng.uniform(-0.1, 0.1, 50)
    dec = dec0 + rng.uniform(-0.1, 0.1, 50)
    xi, eta = A.to_tangent(ra, dec, ra0, dec0)
    ra2, dec2 = A.from_tangent(xi, eta, ra0, dec0)
    err = np.max(np.hypot((ra2 - ra) * np.cos(np.radians(dec)), dec2 - dec)) * 3.6e6
    check("gnomonic projection round trip", err < 1e-3, f"max {err:.1e} mas")
    # 1 mas/yr at 1 kpc = 4.7405 km/s (1 AU / 1 yr).
    check("4.74 km/s per mas/yr per kpc (1 mas at 1 kpc is 1 AU)",
          abs(config.KMS_PER_MASYR_KPC - 1.495978707e8 / 3.15576e7) < 1e-4)


def test_plate_solution():
    print("Plate solution from a synthetic star field")
    ny = nx = 600
    # A made-up plate: 1.0"/px, rotated 0.3 deg, with a quadratic distortion.
    x = rng.uniform(20, nx - 20, 300)
    y = rng.uniform(20, ny - 20, 300)
    u, v = (x - 300) / 300, (y - 300) / 300
    rot = np.radians(0.3)
    xi = -(np.cos(rot) * (x - 300) - np.sin(rot) * (y - 300)) * 1.0 + 0.8 * u * u
    eta = (np.sin(rot) * (x - 300) + np.cos(rot) * (y - 300)) * 1.0 - 0.6 * u * v
    ra, dec = A.from_tangent(xi, eta, 83.633, 22.0145)
    img = np.zeros((ny, nx)) + 100.0
    yy, xx = np.mgrid[:ny, :nx]
    for xs, ys in zip(x, y):
        sl = (slice(int(ys) - 6, int(ys) + 7), slice(int(xs) - 6, int(xs) + 7))
        img[sl] += 500 * np.exp(-0.5 * ((xx[sl] - xs) ** 2 + (yy[sl] - ys) ** 2) / 1.3**2)
    img += rng.normal(0, 2, img.shape)
    gaia = A.GaiaStars(np.arange(len(ra)), ra, dec, np.zeros_like(ra), np.zeros_like(ra), np.zeros_like(ra),
                       np.full(len(ra), 16.0), np.zeros_like(ra), np.ones_like(ra))
    from astropy.wcs import WCS

    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [83.633, 22.0145]
    w.wcs.crpix = [301.5, 301.5]  # deliberately off by ~1.5 px from the truth
    w.wcs.cdelt = [-1 / 3600, 1 / 3600]
    sol = A.solve_plate("synthetic", img, w, 2016.0, gaia, 83.633, 22.0145, star_radius_px=1.3)
    gx, gy = np.meshgrid(np.linspace(30, 570, 10), np.linspace(30, 570, 10))
    gu, gv = (gx - 300) / 300, (gy - 300) / 300
    txi = -(np.cos(rot) * (gx - 300) - np.sin(rot) * (gy - 300)) + 0.8 * gu * gu
    teta = (np.sin(rot) * (gx - 300) + np.cos(rot) * (gy - 300)) - 0.6 * gu * gv
    pxi, peta = sol.to_sky(gx.ravel(), gy.ravel())
    err = np.max(np.hypot(pxi - txi.ravel(), peta - teta.ravel()))
    check("cubic plate solution recovers a distorted, rotated plate", err < 0.03 and sol.n_stars > 200,
          f"{sol.n_stars} stars, rms {sol.rms_arcsec:.3f}\", max map error {err:.3f}\"")
    bx, by = sol.to_pixel(pxi, peta)
    check("inverse map returns the pixels", np.max(np.hypot(bx - gx.ravel(), by - gy.ravel())) < 0.02)


def test_correlation_and_tracking():
    print("Correlation and patch tracking")
    a = rng.normal(size=(60, 60))
    tpl = a[20:36, 22:38]
    cc = match_template(a, tpl)
    mc = M.masked_ncc(tpl, np.ones_like(tpl, bool), a, np.ones_like(a, bool), min_overlap=0.0)
    check("masked NCC equals skimage's match_template with no mask", np.allclose(cc, mc, atol=1e-6),
          f"max diff {np.max(np.abs(cc - mc)):.1e}")
    # A smooth random texture expanding about a centre by 7% between two
    # "epochs" 60 yr apart: the fit must return the centre and T = t1 - 60/0.07.
    grid = A.Grid(83.633, 22.0145, 0.5, 600)
    from scipy import ndimage

    tex = ndimage.gaussian_filter(rng.normal(size=(600, 600)), 4.0)
    s, c = 1.07, (12.0, -5.0)  # stretch, centre (xi, eta) in arcsec
    rr, cc_ = np.mgrid[:600, :600].astype(float)
    cpx = 299.5
    c_col, c_row = cpx - c[0] / 0.5, cpx + c[1] / 0.5
    img2 = ndimage.map_coordinates(tex, [c_row + (rr - c_row) / s, c_col + (cc_ - c_col) / s], order=3)
    region = np.zeros((600, 600), bool)
    region[150:450, 150:450] = True
    ep1 = np.full((600, 600), 1951.0)
    ep2 = np.full((600, 600), 2011.0)
    ok = np.ones((600, 600), bool)
    p1 = M.Prepared("a", ep1, tex, tex, ~ok, ok)
    p2 = M.Prepared("b", ep2, img2, img2, ~ok, ok)
    v = M.track(p1, p2, grid, region, patch=20.0, step=10.0, search=20.0)
    v = v.subset(M.quality_cut(v))
    fit = M.fit_expansion(v)
    T_true = 1951.0 - 60.0 / (s - 1.0)
    check("tracking + fit recover a known convergence date", abs(fit.T - T_true) < 3.0,
          f"T = {fit.T:.1f} vs {T_true:.1f} from {fit.n} vectors")
    check("... and the centre", np.hypot(fit.xi_c - c[0], fit.eta_c - c[1]) < 1.0,
          f"({fit.xi_c:.2f}, {fit.eta_c:.2f}) vs {c}")


def test_stretch_fit():
    print("Poisson stretch fit (the Chandra method)")
    from scipy import ndimage

    n = 300
    yy, xx = np.mgrid[:n, :n].astype(float)
    r = np.hypot(xx - 150, yy - 150)
    ref = 50.0 * np.exp(-0.5 * ((r - 90) / 4.0) ** 2) + 5.0 * ndimage.gaussian_filter(rng.random((n, n)), 2) + 0.5
    s_true, d_true = 1.012, (0.4, -0.3)
    model = X.StretchModel(ref * 20, (150.0, 150.0))  # a "deep" reference, 20x the counts
    truth = X.StretchModel(ref, (150.0, 150.0))
    sel = np.ones((n, n), bool)
    expect = truth.predict(s_true, *d_true, sel).reshape(n, n)
    counts = rng.poisson(expect).astype(float)
    ep = X.Epoch("test", (0,), 2010.0, 5e4, counts, np.ones((n, n), bool))
    g = X.fit_global(model, ep, r < 120)
    check("global stretch recovered", abs(g["s"] - s_true) < 5 * g["s_err"] + 2e-4,
          f"s = {g['s']:.5f} +- {g['s_err']:.5f} vs {s_true}")
    check("pointing shift recovered", abs(g["dx"] - d_true[0]) < 0.05 and abs(g["dy"] - d_true[1]) < 0.05,
          f"({g['dx']:.3f}, {g['dy']:.3f}) vs {d_true}")
    rate = X.rate_fit(np.array([2000.0, 2010.0, 2020.0]), np.array([0.98, 1.0, 1.02]), np.full(3, 1e-4), 2010.0)
    check("rate from a straight line of stretches", abs(rate["a"] - 0.002) < 1e-9, f"{rate['a']:.6f} /yr")


def test_blast_waves():
    print("Blast waves")
    # Sedov-Taylor: R = 1.15 (E t^2 / rho)^(1/5); for E = 1e51 erg, n = 1,
    # mu = 1.4, t = 1000 yr: rho = 2.342e-24 g/cm^3, R = 1.536e19 cm = 4.98 pc.
    rho = 1.4 * 1.6726e-24
    expect = 1.15 * (1e51 * (1000 * 3.15576e7) ** 2 / rho) ** 0.2 / 3.0857e18
    got = R.sedov_radius_pc(1000.0)
    check("Sedov-Taylor radius at 1000 yr", abs(got - expect) < 1e-9 and abs(got - 4.98) < 0.01, f"{got:.3f} pc")
    check("R grows as t^0.4", abs(np.log(R.sedov_radius_pc(2000) / R.sedov_radius_pc(1000)) / np.log(2) - 0.4) < 1e-9)
    check("coasting 5000 km/s for 1000 yr = 5.11 pc", abs(R.free_radius_pc(1000, 5000) - 5.114) < 0.001)


def test_galactic_rotation():
    print("Galactic rotation (Hobbs et al. 2005 model)")
    from astropy import units as u
    from astropy.coordinates import CartesianDifferential, Galactocentric, SkyCoord

    # Compare our flat-rotation-curve proper motion with astropy's frame
    # transformation for a star moving on a circular orbit at 225 km/s.
    sun_pec = config.HOBBS_SUN_PECULIAR
    l0, b0 = np.radians(sun_pec[1]), np.radians(sun_pec[2])
    vpec = sun_pec[0] * np.array([np.cos(b0) * np.cos(l0), np.cos(b0) * np.sin(l0), np.sin(b0)])
    gc_frame = Galactocentric(galcen_distance=config.HOBBS_R0_KPC * u.kpc, z_sun=0 * u.pc, roll=0 * u.deg,
                              galcen_v_sun=CartesianDifferential([vpec[0], config.HOBBS_V0_KMS + vpec[1], vpec[2]]
                                                                 * u.km / u.s))
    worst = 0.0
    for l_deg, b_deg, d in ((120.0, 2.0, 3.0), (300.0, -5.0, 2.0), (30.0, 10.0, 1.0)):
        c = SkyCoord(l=l_deg * u.deg, b=b_deg * u.deg, distance=d * u.kpc, frame="galactic")
        g = c.transform_to(gc_frame)
        X_, Y_ = g.x.to_value(u.kpc), g.y.to_value(u.kpc)
        R_ = np.hypot(X_, Y_)
        # astropy's Galactocentric: Sun at x = -R0, moving towards +y, so a
        # circular orbit has v = V0 (Y, -X) / R (at the Sun: (0, +V0)).
        vx, vy = config.HOBBS_V0_KMS * Y_ / R_, -config.HOBBS_V0_KMS * X_ / R_
        star = SkyCoord(x=g.x, y=g.y, z=g.z, v_x=vx * u.km / u.s, v_y=vy * u.km / u.s, v_z=0 * u.km / u.s,
                        frame=gc_frame).transform_to("galactic")
        ml, mb = K.galactic_rotation_pm(l_deg, b_deg, d)
        worst = max(worst, abs(ml - star.pm_l_cosb.to_value(u.mas / u.yr)), abs(mb - star.pm_b.to_value(u.mas / u.yr)))
    check("rotation proper motions agree with astropy's Galactocentric frame", worst < 0.05, f"worst {worst:.3f} mas/yr")
    sig, _ = K.rayleigh_sigma(np.abs(rng.normal(0, 265, (20000, 2))).dot([0, 0]) + np.hypot(*rng.normal(0, 265, (2, 20000))))
    check("Rayleigh sigma recovers 265 km/s from projected Maxwellian speeds", abs(sig - 265) < 5, f"{sig:.1f}")


def main() -> int:
    for t in (test_dates_and_projection, test_plate_solution, test_correlation_and_tracking, test_stretch_fit,
              test_blast_waves, test_galactic_rotation):
        t()
    n_fail = sum(1 for _, ok, _ in RESULTS if not ok)
    print(f"\n{len(RESULTS) - n_fail} of {len(RESULTS)} checks passed.")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
