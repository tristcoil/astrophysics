"""Self-checks for every transform and estimator, against known answers.

    python test_physics.py

No network and no test framework needed (the McMillan 2017 grid must have
been built once, which run.py or fetch-free `python -c "import galaxy;
galaxy.mcmillan17()"` does).

Everything in this project flows through these functions. A sign error in a
velocity transform, or a spiral-dating estimator that prefers some age for no
reason, would produce figures that look plausible and are wrong -- so each is
checked against an independent implementation (astropy, gala, galpy, a direct
orbit integration) or against a case whose answer is known.

The estimator checks (anisotropy, spiral dating) use simulated stars. That is
the one place this project uses synthetic data: to prove the code recovers a
*known* answer. No figure is ever made from it.
"""

from __future__ import annotations

import re
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")

FAILURES: list[str] = []


def check(label, got, expected, tol, unit="", absolute=False):
    got, expected = float(got), float(expected)
    err = abs(got - expected) if absolute else abs(got - expected) / abs(expected)
    ok = err <= tol
    if not ok:
        FAILURES.append(label)
    shown = f"{err:.2g} off" if absolute else f"{err * 100:.3g}% off"
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<58} {got:>12.6g} vs {expected:>12.6g} {unit:<9}"
          f" ({shown})")


def claim(label, condition):
    if not condition:
        FAILURES.append(label)
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label}")


def adql_to_numpy(expr: str, cols: dict):
    """Evaluate one of our ADQL expressions with numpy, to prove the server-side
    cuts compute what the local code computes."""
    py = expr
    for fn, np_fn in [("RADIANS", "np.radians"), ("DEGREES", "np.degrees"), ("SQRT", "np.sqrt"),
                      ("POWER", "np.power"), ("ATAN2", "np.arctan2"), ("ASIN", "np.arcsin"),
                      ("COS", "np.cos"), ("SIN", "np.sin")]:
        py = re.sub(rf"\b{fn}\(", f"{np_fn}(", py)
    py = py.replace("g.", "")
    return eval(py, {"np": np}, cols)


def random_stars(n=2000, seed=0):
    rng = np.random.default_rng(seed)
    ra = rng.uniform(0, 360, n)
    dec = np.degrees(np.arcsin(rng.uniform(-1, 1, n)))
    return dict(ra=ra, dec=dec, parallax=rng.uniform(0.2, 5, n), pmra=rng.normal(0, 30, n),
                pmdec=rng.normal(0, 30, n), radial_velocity=rng.normal(0, 80, n))


def main() -> int:
    import astropy.coordinates as coord
    import astropy.units as u

    import galaxy as gx
    import photometry as ph
    from config import R0_KPC, V_SUN_KMS, Z_SUN_KPC

    s = random_stars()
    dist = 1 / s["parallax"]

    print("GD-1 frame (against gala's GD1Koposov10)")
    import gala.coordinates as gc

    c = coord.SkyCoord(ra=s["ra"] * u.deg, dec=s["dec"] * u.deg,
                       pm_ra_cosdec=s["pmra"] * u.mas / u.yr, pm_dec=s["pmdec"] * u.mas / u.yr)
    g1 = c.transform_to(gc.GD1Koposov10())
    p1, p2, m1, m2, _ = gx.gd1_coordinates(s["ra"], s["dec"], s["pmra"], s["pmdec"])
    dphi1 = np.abs((p1 - g1.phi1.wrap_at(180 * u.deg).deg + 180) % 360 - 180)
    check("max |phi1 - gala| (deg)", dphi1.max(), 0, 1e-8, "deg", absolute=True)
    # Koposov et al. give the matrix to 10 digits, so it is orthonormal only to
    # ~1e-10; gala and we round differently at the 0.1-milliarcsecond level.
    check("max |phi2 - gala| (deg)", np.abs(p2 - g1.phi2.deg).max(), 0, 1e-6, "deg", absolute=True)
    check("max |pm_phi1 - gala|", np.abs(m1 - g1.pm_phi1_cosphi2.value).max(), 0, 1e-8,
          "mas/yr", absolute=True)
    check("max |pm_phi2 - gala|", np.abs(m2 - g1.pm_phi2.value).max(), 0, 1e-8,
          "mas/yr", absolute=True)
    back = gx.gd1_to_icrs(p1, p2)
    check("GD-1 -> ICRS round trip (deg)", np.abs(((back[0] - s["ra"]) + 180) % 360 - 180).max(),
          0, 1e-6, "deg", absolute=True)
    # The PWB2018 progenitor field: (RA, Dec) = (177.01, 53.99) lies near phi1 = -14.
    q1, q2 = gx.gd1_coordinates(np.array([177.01]), np.array([53.99]))
    claim(f"PWB2018's densest field (177.01, 53.99) maps to phi1 = {q1[0]:.1f}, "
          f"phi2 = {q2[0]:.2f} (their -18..-10)", -18 < q1[0] < -10 and abs(q2[0]) < 1)

    print("\nGalactocentric transform (against astropy)")
    frame = coord.Galactocentric(galcen_distance=R0_KPC * u.kpc, z_sun=Z_SUN_KPC * u.kpc,
                                 galcen_v_sun=coord.CartesianDifferential(list(V_SUN_KMS) * u.km / u.s))
    c6 = coord.SkyCoord(ra=s["ra"] * u.deg, dec=s["dec"] * u.deg, distance=dist * u.kpc,
                        pm_ra_cosdec=s["pmra"] * u.mas / u.yr, pm_dec=s["pmdec"] * u.mas / u.yr,
                        radial_velocity=s["radial_velocity"] * u.km / u.s).transform_to(frame)
    g = gx.galactocentric(s["ra"], s["dec"], dist, s["pmra"], s["pmdec"], s["radial_velocity"])
    check("max position difference (kpc)",
          np.abs(np.array([g.x, g.y, g.z]) - c6.cartesian.xyz.to_value(u.kpc)).max(), 0, 1e-9,
          "kpc", absolute=True)
    v_ast = c6.velocity.d_xyz.to_value(u.km / u.s)
    check("max velocity difference (km/s)", np.abs(np.array([g.vx, g.vy, g.vz]) - v_ast).max(),
          0, 1e-8, "km/s", absolute=True)
    sun = gx.galactocentric(np.array([0.0]), np.array([0.0]), np.array([1e-9]), np.array([0.0]),
                            np.array([0.0]), np.array([0.0]))
    check("the Sun's v_phi (prograde positive)", sun.v_phi[0], 245.6, 0.001, "km/s")
    check("the Sun's L_z = R0 x 245.6", sun.L_z[0], R0_KPC * 245.6, 0.002, "kpc km/s")
    claim("the Sun sits at azimuth 180 deg", abs(sun.phi_deg[0] - 180) < 0.01)

    print("\nServer-side cuts (the ADQL strings, evaluated with numpy)")
    import fetch_data as fd

    uvw = gx.heliocentric_uvw(s["ra"], s["dec"], dist, s["pmra"], s["pmdec"], s["radial_velocity"])
    got = np.array([adql_to_numpy(e, s) for e in fd.adql_uvw()])
    check("ADQL (U,V,W) vs numpy, max |diff|", np.abs(got - uvw).max(), 0, 1e-8, "km/s",
          absolute=True)
    check("ADQL |v - v_LSR| vs numpy, max |diff|",
          np.abs(adql_to_numpy(fd.adql_speed_vs_lsr(), s) - gx.speed_vs_lsr(*uvw)).max(),
          0, 1e-8, "km/s", absolute=True)
    ast_uvw = coord.SkyCoord(ra=s["ra"] * u.deg, dec=s["dec"] * u.deg, distance=dist * u.kpc,
                             pm_ra_cosdec=s["pmra"] * u.mas / u.yr, pm_dec=s["pmdec"] * u.mas / u.yr,
                             radial_velocity=s["radial_velocity"] * u.km / u.s).galactic.velocity.d_xyz
    check("numpy (U,V,W) vs astropy Galactic, max |diff|",
          np.abs(uvw - ast_uvw.to_value(u.km / u.s)).max(), 0, 1e-8, "km/s", absolute=True)
    ph1, ph2 = (adql_to_numpy(e, s) for e in fd.adql_gd1())
    check("ADQL phi1 vs numpy, max |diff|", np.abs(ph1 - p1).max(), 0, 1e-8, "deg", absolute=True)
    check("ADQL phi2 vs numpy, max |diff|", np.abs(ph2 - p2).max(), 0, 1e-8, "deg", absolute=True)
    l, b = gx.galactic_lb(s["ra"], s["dec"])
    R_adql = adql_to_numpy(fd.adql_galactocentric_radius(), dict(s, l=l, b=b))
    check("ADQL ring radius vs full transform, max |diff|", np.abs(R_adql - g.R).max(), 0, 0.02,
          "kpc", absolute=True)
    near = dist < 1.5
    claim(f"...within 5 pc for stars < 1.5 kpc (max {np.abs(R_adql - g.R)[near].max() * 1000:.1f} pc;"
          f" the fetch margin)", np.abs(R_adql - g.R)[near].max() < 0.005)
    from gaia_archive import healpix_source_id_range

    check("HEALPix level-1 pixel 1 starts at 2**57", healpix_source_id_range(1, 1)[0], 2.0**57, 0)

    print("\nSolar reflex correction (against gala.coordinates.reflex_correct)")
    cc = coord.SkyCoord(ra=s["ra"] * u.deg, dec=s["dec"] * u.deg, distance=dist * u.kpc,
                        pm_ra_cosdec=s["pmra"] * u.mas / u.yr, pm_dec=s["pmdec"] * u.mas / u.yr,
                        radial_velocity=s["radial_velocity"] * u.km / u.s)
    rc = gc.reflex_correct(cc, frame)
    a, bb = gx.reflex_correct_pm(s["ra"], s["dec"], dist, s["pmra"], s["pmdec"])
    check("max |pmra difference|", np.abs(a - rc.pm_ra_cosdec.value).max(), 0, 1e-8, "mas/yr",
          absolute=True)
    check("max |pmdec difference|", np.abs(bb - rc.pm_dec.value).max(), 0, 1e-8, "mas/yr",
          absolute=True)

    print("\nPotentials")
    from galpy.potential import (HernquistPotential, MiyamotoNagaiPotential, NFWPotential,
                                 evaluatePotentials, vcirc)

    from config import HELMI17

    h = gx.AnalyticHelmi17()
    ro, vo = 8.0, 220.0
    G = 4.300917270038e-06
    amp_mass = lambda M: M * G / (vo**2 * ro)  # galpy natural mass units
    mn = MiyamotoNagaiPotential(amp=amp_mass(HELMI17.disc_mass_msun), a=HELMI17.disc_a_kpc / ro,
                                b=HELMI17.disc_b_kpc / ro)
    # galpy's Hernquist amp is 2 x mass (its density normalisation).
    hq = HernquistPotential(amp=2 * amp_mass(HELMI17.bulge_mass_msun), a=HELMI17.bulge_c_kpc / ro)
    cc_ = HELMI17.halo_concentration
    nfw = NFWPotential(amp=amp_mass(HELMI17.halo_mass_msun) / (np.log1p(cc_) - cc_ / (1 + cc_)),
                       a=HELMI17.halo_rs_kpc / ro)
    pts = [(8.122, 0.0), (4.0, 1.0), (15.0, 5.0), (30.0, -10.0)]
    for R, z in pts:
        check(f"Helmi17 Phi({R:g}, {z:g}) vs galpy", h.phi(R, z),
              evaluatePotentials([mn, hq, nfw], R / ro, z / ro, use_physical=False) * vo**2,
              1e-10, "km2/s2")
    check("Helmi17 v_c(R0) vs galpy", h.vcirc(R0_KPC),
          vcirc([mn, hq, nfw], R0_KPC / ro, use_physical=False) * vo, 1e-10, "km/s")
    m17 = gx.mcmillan17()
    from galpy.potential.mwpotentials import McMillan17

    rng = np.random.default_rng(5)
    Rs, zs = rng.uniform(0.3, 40, 12), rng.uniform(-8, 8, 12)
    direct = np.array([evaluatePotentials(McMillan17, r / 8.21, zz / 8.21, use_physical=False)
                       for r, zz in zip(Rs, zs)]) * 233.1**2
    check("McMillan17 grid vs direct galpy, max rel. error", np.max(np.abs(m17.phi(Rs, zs) / direct - 1)),
          0, 1e-6, "", absolute=True)
    check("McMillan17 v_c(8.21 kpc) (McMillan 2017: 233.1)", m17.vcirc(8.21), 233.1, 0.001, "km/s")
    # Circularity: a circular orbit has circularity +1 (prograde) or -1.
    R = np.array([5.0, 8.0, 12.0])
    E = m17.phi(R, np.zeros(3)) + 0.5 * m17.vcirc(R) ** 2
    circ = gx.circularity(m17, E, R * m17.vcirc(R))
    check("circularity of circular orbits (max |1 - circ|)", np.max(np.abs(circ - 1)), 0, 1e-3, "",
          absolute=True)

    print("\nVertical orbits (against galpy's own vertical force, integrated directly)")
    from galpy.potential import evaluatez2derivs, evaluatezforces
    from scipy.integrate import solve_ivp

    import spiral as sp

    vo_ = sp.VerticalOrbits(m17)
    to_kms2_per_kpc = 233.1**2 / 8.21  # galpy natural force units -> (km/s)^2 / kpc

    def fz(z):
        return evaluatezforces(McMillan17, R0_KPC / 8.21, z / 8.21, use_physical=False) \
            * to_kms2_per_kpc

    for amp in (0.1, 0.5, 1.0):
        vz0 = np.sqrt(2 * vo_.vertical_energy(amp))
        # z in kpc, v_z in km/s, time in kpc/(km/s).
        sol = solve_ivp(lambda t, y: [y[1], fz(y[0])], (0, 2.0), [0.0, vz0],
                        events=lambda t, y: y[0], dense_output=True, rtol=1e-8, atol=1e-10)
        crossings = sol.t_events[0][sol.t_events[0] > 1e-6]
        period = crossings[1] * sp.KPC_PER_KMS_MYR  # z = 0 again, moving up
        check(f"vertical period at amplitude {amp} kpc", 2 * np.pi / vo_.frequency(amp), period,
              0.002, "Myr")
        ts = np.linspace(0.02, 0.95, 40) * crossings[1]
        zz, vv = sol.sol(ts)
        th, om, _ = vo_.angle_frequency(zz, vv)
        pred = np.mod(om * ts * sp.KPC_PER_KMS_MYR, 2 * np.pi)
        diff = np.abs((th - pred + np.pi) % (2 * np.pi) - np.pi)
        check(f"angle = Omega t along that orbit (max |diff|, rad)", diff.max(), 0, 0.01, "rad",
              absolute=True)
    nu_galpy = np.sqrt(evaluatez2derivs(McMillan17, R0_KPC / 8.21, 0.0, use_physical=False)) \
        * 233.1 / 8.21 / sp.KPC_PER_KMS_MYR
    check("harmonic frequency vs galpy's d2Phi/dz2 at (R0, 0)", vo_.harmonic_frequency(),
          nu_galpy, 0.005, "rad/Myr")

    print("\nSpiral dating recovers a known age (simulated stars, test only)")
    rng = np.random.default_rng(7)
    n = 400_000
    amp = rng.uniform(0.2, 1.0, n)
    om = vo_.frequency(amp)
    for t_true in (300.0, 550.0):
        # Half the stars phase-mixed, half kicked together t_true ago.
        theta = np.where(rng.random(n) < 0.5, rng.uniform(0, 2 * np.pi, n),
                         np.mod(rng.normal(0.3, 0.6, n) + om * t_true, 2 * np.pi))
        res = sp.date_spiral(theta, om, np.ones(n), m=1, n_boot=5)
        check(f"recovered age for t = {t_true:.0f} Myr", res["t_myr"], t_true, 0.02, "Myr")
    theta = rng.uniform(0, 2 * np.pi, n)
    amp_null, *_ = sp.phase_stack(theta, om, np.ones(n))
    claim(f"no spiral -> no peak (max stacked amplitude {amp_null.max():.4f} < 0.01)",
          amp_null.max() < 0.01)

    print("\nVelocity anisotropy (simulated stars, test only)")
    import halo as hl

    for beta_true in (0.0, 0.9):
        s_r = 150.0
        s_t = s_r * np.sqrt(1 - beta_true)
        v = rng.normal(0, [s_r, s_t, s_t], (6000, 3))
        e = hl.velocity_ellipsoid(v[:, 0], v[:, 1], v[:, 2], max_disc=0)
        check(f"beta = {beta_true} recovered (halo only)", e["beta"], beta_true, 0.05, "",
              absolute=True)
    # A halo hidden under a rotating disc three times its size: the mixture
    # must still recover the halo's shape.
    halo_v = rng.normal(0, [150.0, 150 * np.sqrt(0.1), 150 * np.sqrt(0.1)], (3000, 3)) + [0, 0, 20]
    disc_v = rng.normal(0, [45.0, 25.0, 30.0], (9000, 3)) + [0, 0, 200]
    v = np.vstack([halo_v, disc_v])
    e = hl.velocity_ellipsoid(v[:, 0], v[:, 1], v[:, 2])
    check("beta = 0.9 recovered under a 3x larger disc", e["beta"], 0.9, 0.05, "", absolute=True)
    check("... and the halo's share (0.25)", e["weight"], 0.25, 0.04, "", absolute=True)

    print("\nPhotometry")
    check("distance modulus at 8.501 kpc (M92)", ph.distance_modulus(8.501), 14.647, 1e-4, "mag")
    k = ph.extinction_coefficient("G", 0.7, 0.1)
    claim(f"k_G for a turn-off star ({k:.3f}) between 0.8 and 0.9", 0.8 < k < 0.9)
    # Babusiaux et al. (2018) Eq. 1 by hand for (BP-RP)_0 = 0.75, A_0 = 0.155:
    # k_BP = 1.071892, k_RP = 0.594332, so E(BP-RP)/E(B-V) = 3.1 x 0.477560 = 1.48044.
    c_true, ebv = 0.75, 0.05
    a0 = 3.1 * ebv
    kbp = ph.extinction_coefficient("BP", c_true, a0)
    krp = ph.extinction_coefficient("RP", c_true, a0)
    check("E(BP-RP)/E(B-V) at (BP-RP)0 = 0.75 (hand value 1.48044)", (kbp - krp) * 3.1, 1.48044, 1e-5)
    kg = ph.extinction_coefficient("G", c_true, a0)
    g0, c0 = ph.deredden(np.array([18.0 + kg * a0]), np.array([c_true + (kbp - krp) * a0]),
                         np.array([ebv]))
    check("deredden() recovers a reddened star's colour", c0[0], c_true, 1e-6, "mag", absolute=True)
    check("... and its magnitude", g0[0], 18.0, 1e-6, "mag", absolute=True)

    print()
    if FAILURES:
        print(f"{len(FAILURES)} check(s) FAILED: {FAILURES}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
