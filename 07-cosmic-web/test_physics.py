"""Self-checks for the physics and the clustering machinery.

Run with ``python test_physics.py``. No test framework and no network needed:
these check the formulas and the estimators against cases with known answers,
not the data.

Every result in this project flows through these functions. A slip in a
distance, a K-correction or the pair counter would produce a wrong map or a
wrong correlation length that still looks entirely plausible on a plot, so
each gets checked against something independent.
"""

from __future__ import annotations

import sys

import numpy as np
from astropy import units as u

import clustering as cl
import cosmo_physics as cp

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float, unit: str = "",
          *, absolute: bool = False) -> None:
    """Assert ``got`` is within ``tol`` of ``expected`` (fractional unless ``absolute``)."""
    err = abs(got - expected) if (absolute or expected == 0) else abs(got - expected) / abs(expected)
    ok = err <= tol
    if not ok:
        FAILURES.append(label)
    kind = "abs" if (absolute or expected == 0) else "rel"
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<52} {got:>12.6g} vs {expected:>12.6g} "
          f"{unit:<6} ({kind} err {err:.2g}, tol {tol:g})")


def test_distances() -> None:
    print("Distances")
    # Nearby, redshift is just velocity: Hubble's law d = cz / H0.
    z = 0.001
    hubble = cp.C_KMS * z / cp.COSMO.H0.to_value(u.km / u.s / u.Mpc)
    check("Hubble's law limit, D_C(z=0.001)", float(cp.comoving_distance(z)), hubble, 1e-3, "Mpc")
    # The interpolation table must reproduce astropy's own integral.
    for z in (0.02, 0.1, 0.3):
        exact = cp.COSMO.comoving_distance(z).to_value(u.Mpc)
        check(f"interpolated D_C(z={z}) vs astropy", float(cp.comoving_distance(z)), exact, 1e-7, "Mpc")
    check("distance modulus z=0.1 vs astropy", float(cp.distance_modulus(0.1)),
          cp.COSMO.distmod(0.1).value, 1e-6, "mag")
    # In h^-1 Mpc, D_C(z = 0.1) is close to the ~ 293 h^-1 Mpc of any Om ~ 0.3 model.
    check("D_C(z=0.1) in h^-1 Mpc", float(cp.comoving_distance_h(0.1)), 292.9, 3e-3, "Mpc/h")
    check("inverse: z at D_C(0.0734)", float(cp.redshift_at_distance(cp.comoving_distance(0.0734))),
          0.0734, 1e-6)
    # A +/-1000 km/s window at low z spans ~ +/-10 h^-1 Mpc (Hubble's law, H0 = 100 h).
    check("LOS depth of 1000 km/s at z=0.001", float(cp.los_window(1000.0, 0.001)), 10.0, 3e-3,
          "Mpc/h")


def test_k_corrections() -> None:
    print("\nK-corrections (Chilingarian et al. 2010, 2012)")
    # The worked example in the authors' calc_kcor.py documentation.
    check("K_g(z=0.2, g-r=1.1), their documented example",
          float(cp.k_correction("g", 0.2, 1.1)), 0.5209713976, 1e-8, "mag")
    # Values from running the authors' calc_kcor.py itself, guarding against a
    # slip when the r-band coefficient table was transcribed. (Their IDL
    # header's g-band example, 0.481612 at z = 0.15, predates the 2012
    # coefficients: their own current code, and ours, give 0.429621.)
    for z, c, ref in ((0.05, 0.6, 0.041234262625), (0.1, 0.9, 0.0840588705),
                      (0.2, 1.2, 0.191570128)):
        check(f"K_r(z={z}, g-r={c}) vs authors' code",
              float(cp.k_correction("r", z, c)), ref, 1e-8, "mag")
    # No redshift, nothing to correct.
    check("K_r(z=0) for any colour", float(cp.k_correction("r", 0.0, 0.9)), 0.0, 1e-12, "mag")
    # Red galaxies at z = 0.1: K_r ~ 0.1 mag, the familiar size of the correction.
    k = float(cp.k_correction("r", 0.1, 0.95))
    check("K_r(z=0.1, red g-r=0.95) is ~0.1 mag", k, 0.1, 0.05, "mag", absolute=True)


def test_magnitudes() -> None:
    print("\nMagnitudes")
    m, z, k = 17.77, 0.1, 0.12
    big_m = float(cp.absolute_magnitude(m, z, k))
    check("round trip m -> M -> m", big_m + float(cp.distance_modulus(z)) + k, m, 1e-12, "mag")
    # Quoting for h = 1 brightens nothing physically; it shifts by -5 log10(h).
    shift = float(cp.to_literature_magnitude(0.0, 0.1, 0.0, 0.0))
    check("h = 1 convention shift, -5 log10(0.6766)", shift, 0.8483, 1e-3, "mag")
    # Evolution: seen at z = 0.02, a galaxy was brighter at z = 0.1 by Q x 0.08.
    ev = float(cp.to_literature_magnitude(0.0, 0.02, 0.0, 1.62)) - shift
    check("evolution to z=0.1 from z=0.02", ev, -0.1296, 1e-3, "mag")


def test_geometry() -> None:
    print("\nLine-of-sight geometry")
    # A pair lined up along the line of sight: all separation is pi.
    rp, pi = cp.rp_pi([[0, 0, 100.0]], [[0, 0, 110.0]])
    check("radial pair: pi", float(pi[0]), 10.0, 1e-12)
    check("radial pair: r_p", float(rp[0]), 0.0, 1e-9, absolute=True)
    # A pair side by side at the same distance: all separation is r_p.
    rp, pi = cp.rp_pi([[-5.0, 0, 100.0]], [[5.0, 0, 100.0]])
    check("transverse pair: r_p", float(rp[0]), 10.0, 1e-12)
    check("transverse pair: pi", float(pi[0]), 0.0, 1e-12, absolute=True)
    # Round trip through Cartesian coordinates.
    xyz = cp.cartesian([123.4], [-1.5], [250.0])
    check("|cartesian(ra, dec, d)| = d", float(np.linalg.norm(xyz[0])), 250.0, 1e-12)


def test_power_law_projection() -> None:
    print("\nProjected correlation function of a power law")
    # Integrate xi(r) = (r/r0)^-gamma along the line of sight numerically and
    # compare with the closed form used for fitting.
    r0, g = 5.5, 1.8
    for rp in (0.3, 3.0, 15.0):
        pi = np.logspace(-6, 6, 200_001)
        xi = (np.sqrt(rp**2 + pi**2) / r0) ** (-g)
        numeric = 2.0 * np.trapezoid(xi, pi)
        check(f"w_p({rp}) numeric vs closed form", numeric, float(cp.wp_power_law(rp, r0, g)),
              2e-3)


def test_velocity_dispersion() -> None:
    print("\nVelocity dispersion")
    # Members at z0 +/- dz: the rest-frame dispersion is c dz / (1 + z0),
    # times the sqrt(n / (n - 1)) of the unbiased standard deviation.
    z0, dz, n = 0.08, 0.003, 10
    z = np.array([z0 - dz, z0 + dz] * (n // 2))
    sigma, _ = cp.velocity_dispersion(z)
    check("sigma_v of a symmetric pair set", sigma,
          cp.C_KMS * dz / (1 + z0) * np.sqrt(n / (n - 1)), 1e-12, "km/s")


def test_cluster_dispersion() -> None:
    print("\nCluster velocity dispersion with 3-sigma clipping")
    from scipy.spatial import cKDTree

    import analysis as an

    rng = np.random.default_rng(4)
    z0 = 0.05

    def measure(v):
        z = z0 + (1 + z0) * v / cp.C_KMS
        ra = 180.0 + rng.normal(0.0, 0.05, len(z))
        dec = rng.normal(0.0, 0.05, len(z))
        g = {"ra": ra, "dec": dec, "z": z}
        tree = cKDTree(cp.cartesian(ra, dec, np.ones(len(ra))))
        return an.cluster_dispersion(g, tree, 180.0, 0.0, z0)["sigma_v"]

    # A clean cluster: clipping plus the truncation correction must return the
    # true dispersion (statistical precision ~1% with 5,000 members).
    clean = rng.normal(0.0, 800.0, 5000)
    check("clean cluster, true sigma 800 km/s", measure(clean), 800.0, 0.03, "km/s")
    # With 10% unrelated galaxies spread over +/-3000 km/s, clipping sheds those
    # beyond 3 sigma but not those inside the window, which still inflate the
    # answer. A known limitation: the result must land between the truth and
    # the unclipped standard deviation.
    mixed = np.concatenate([rng.normal(0.0, 800.0, 3000), rng.uniform(-3000, 3000, 300)])
    got = measure(mixed)
    naive = float(np.std(mixed, ddof=1))
    ok = 800.0 < got < naive
    if not ok:
        FAILURES.append("clipping with interlopers")
    print(f"  [{'ok  ' if ok else 'FAIL'}] 10% interlopers: clipped {got:.0f} km/s lies between "
          f"the true 800 and the unclipped {naive:.0f}")


def test_randoms() -> None:
    print("\nRandom catalogues")
    rng = np.random.default_rng(1)
    d = cl.random_volume_distances(400_000, 0.0, 100.0, rng)
    # Half the radius holds an eighth of the volume.
    check("uniform-in-volume: fraction inside d/2", float((d < 50.0).mean()), 0.125, 0.01)
    from astropy_healpix import HEALPix

    hp = HEALPix(nside=64, order="nested")
    everything = np.ones(hp.npix, dtype=bool)
    full = cl.Mask(nside=64, counts=everything.astype(int), surveyed=everything,
                   interior=everything)
    check("full-sky mask area", full.area_deg2(), 41252.96, 1e-6, "deg2")
    ra, dec = cl.random_sky(400_000, full, rng)
    # Uniform on a sphere: the band |dec| < 30 deg holds exactly half the area
    # (sin 30 = 1/2). Sampling uniformly in dec instead would put 2/3 there.
    check("random_sky: fraction within |dec| < 30", float((np.abs(dec) < 30).mean()), 0.5, 0.005)
    # HEALPix pixels have equal areas, so the fraction of uniform points landing
    # in any set of pixels equals that set's share of pixels.
    lon, lat = hp.healpix_to_lonlat(np.arange(hp.npix))
    band = np.abs(lat.to_value("deg")) < 30.0
    in_band = band[full.pixel(ra, dec)].mean()
    check("equal-area pixels: point share = pixel share", float(in_band), float(band.mean()), 0.005)
    mask = cl.Mask(nside=64, counts=band.astype(int), surveyed=band, interior=band)
    ra, dec = cl.random_sky(50_000, mask, rng)
    check("random_sky: every point inside a partial mask",
          float(mask.contains(ra, dec).mean()), 1.0, 1e-12)


def test_pair_counter() -> None:
    print("\nPair counter against brute force")
    rng = np.random.default_rng(2)
    # Points in a shell segment far from the origin, like a survey.
    n = 1500
    ra = rng.uniform(170, 190, n)
    dec = np.degrees(np.arcsin(rng.uniform(np.sin(np.radians(-10)), np.sin(np.radians(10)), n)))
    d = cl.random_volume_distances(n, 150.0, 250.0, rng)
    p = cp.cartesian(ra, dec, d)
    reg = (ra > 180).astype(int)
    b = cl.Binning(rp_edges=10 ** np.arange(-0.9, 1.51, 0.2), pi_max=40.0, pi_bin=2.0,
                   map_bin=2.0, map_max=20.0)
    pc = cl.count_pairs(p, reg, b, n_regions=2, workers=1)
    # Brute force every pair with the explicit 3D formula.
    i, j = np.triu_indices(n, k=1)
    rp, pi = cp.rp_pi(p[i], p[j])
    in_bins = (pi < b.pi_max) & (rp >= b.rp_edges[0]) & (rp < b.rp_edges[-1])
    h, _, _ = np.histogram2d(rp[in_bins], pi[in_bins],
                             bins=[b.rp_edges, np.arange(0, b.pi_max + 1e-9, b.pi_bin)])
    check("auto-pair counts, total in bins", float(pc.wp.sum()), float(h.sum()), 0, "pairs")
    check("auto-pair counts, largest bin mismatch", float(np.abs(pc.wp - h).max()), 0.0, 0,
          absolute=True)
    # Jackknife bookkeeping: pairs touching region 0 = all pairs minus those
    # entirely inside region 1.
    both1 = (reg[i] == 1) & (reg[j] == 1) & in_bins
    check("pairs touching region 0", float(pc.wp_region[0].sum()),
          float(in_bins.sum() - both1.sum()), 0, "pairs")

    # Cross pairs.
    q = cp.cartesian(rng.uniform(170, 190, 900), rng.uniform(-10, 10, 900),
                     cl.random_volume_distances(900, 150.0, 250.0, rng))
    pcx = cl.count_pairs(p, reg, b, n_regions=2, p2=q, reg2=np.zeros(900, int), workers=1)
    ii, jj = np.meshgrid(np.arange(n), np.arange(900), indexing="ij")
    rp, pi = cp.rp_pi(p[ii.ravel()], q[jj.ravel()])
    sel = (pi < b.pi_max) & (rp >= b.rp_edges[0]) & (rp < b.rp_edges[-1])
    check("cross-pair counts, total in bins", float(pcx.wp.sum()), float(sel.sum()), 0, "pairs")
    # Parallel and serial must agree exactly.
    pcp = cl.count_pairs(p, reg, b, n_regions=2, workers=4, target_pairs_per_task=2e4)
    check("4 workers vs 1 worker", float(np.abs(pcp.wp - pc.wp).max()), 0.0, 0, absolute=True)


def test_null_clustering() -> None:
    print("\nNull test: random points are not clustered")
    # Treat one random catalogue as 'data' and measure it against another. The
    # estimator must return w_p consistent with zero -- if it does not, the
    # machinery manufactures clustering on its own.
    rng = np.random.default_rng(3)

    def shell(n):
        ra = rng.uniform(150, 210, n)
        dec = np.degrees(np.arcsin(rng.uniform(np.sin(np.radians(-15)), np.sin(np.radians(15)), n)))
        return cp.cartesian(ra, dec, cl.random_volume_distances(n, 60.0, 250.0, rng)), ra

    fake, fra = shell(12_000)
    rand, rra = shell(36_000)
    b = cl.Binning(rp_edges=10 ** np.arange(-0.1, 1.31, 0.2), pi_max=40.0, pi_bin=2.0,
                   map_bin=2.0, map_max=10.0)
    nreg = 6
    fr = np.minimum(((fra - 150) / 10).astype(int), nreg - 1)
    rr_ = np.minimum(((rra - 150) / 10).astype(int), nreg - 1)
    dd = cl.count_pairs(fake, fr, b, nreg, workers=4)
    dr = cl.count_pairs(fake, fr, b, nreg, p2=rand, reg2=rr_, workers=4)
    rr = cl.count_pairs(rand, rr_, b, nreg, workers=4)
    res = cl.measure(dd, dr, rr, np.bincount(fr, minlength=nreg).astype(float),
                     np.bincount(rr_, minlength=nreg).astype(float), b)
    chi2 = float(np.sum((res.wp / res.wp_err) ** 2))
    print(f"         w_p values: {np.array2string(res.wp, precision=2)}")
    check("chi^2 of w_p against zero, per bin", chi2 / len(res.wp), 1.0, 1.5, absolute=True)


def main() -> int:
    test_distances()
    test_k_corrections()
    test_magnitudes()
    test_geometry()
    test_power_law_projection()
    test_velocity_dispersion()
    test_cluster_dispersion()
    test_randoms()
    test_pair_counter()
    test_null_clustering()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILED: {FAILURES}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
