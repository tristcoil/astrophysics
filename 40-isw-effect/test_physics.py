"""Self-checks against cases with known answers.

Run with ``python test_physics.py`` (no test framework). The first group needs
no data at all; the second uses cached products when they exist and skips
otherwise.
"""

from __future__ import annotations

import sys

import numpy as np
from scipy.special import erf

import cosmology as C

FAILURES: list[str] = []


def check(label, got, expected, tol, *, absolute=False):
    err = abs(got - expected) if absolute else abs(got - expected) / max(abs(expected), 1e-300)
    ok = err <= tol
    if not ok:
        FAILURES.append(label)
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<58} {got:>12.5g} vs {expected:>12.5g}")


def test_background():
    print("Background and growth")
    from config import COSMO, COSMO_OMEGA_M

    h = COSMO["H0"] / 100
    om = (COSMO["ombh2"] + COSMO["omch2"] + COSMO["mnu"] / C.NU_EV) / h**2
    check("Omega_m from D3's densities vs D3's CAMB value", om, COSMO_OMEGA_M, 1e-4)
    # Exact Lambda-CDM growth: f is close to Omega_m(z)^0.55 (Linder 2005)
    for z in (0.0, 0.5, 1.0):
        L = C.LCDM
        omz = L.omega_m * (1 + z) ** 3 / L.E(z) ** 2
        check(f"growth rate f(z={z}) vs Omega_m(z)^0.55", float(L.f(z)), float(omz**0.55), 6e-3)
    for z in (0.3, 1.0, 3.0):
        check(f"Einstein-de Sitter: ISW kernel at z={z} is zero", float(C.EDS.isw_kernel(z)), 0.0,
              1e-15, absolute=True)
        check(f"Einstein-de Sitter: D(z={z}) = 1/(1+z)", float(C.EDS.D(z)), 1 / (1 + z), 1e-4)


def test_growth_vs_camb():
    print("Growth factor against CAMB")
    zs = [0.3, 0.6, 1.0]
    k, z, pk, s8 = C.linear_power(zs)
    i = np.argmin(np.abs(k - 0.02))
    i0 = int(np.argmin(np.abs(z - 0.0)))
    for zz in zs:
        j = int(np.argmin(np.abs(z - zz)))
        d_camb = np.sqrt(pk[j, i] / pk[i0, i])
        # massive neutrinos make CAMB's growth slightly scale dependent: 0.3% budget
        check(f"D({zz}) ours vs CAMB at k = 0.02 h/Mpc", float(C.LCDM.D(zz)), float(d_camb), 3e-3)


def test_spherical_potential():
    print("Spherical Poisson solver against an analytic Gaussian blob")
    d0, r0 = -0.3, 50.0
    r = np.linspace(0, 1500, 150001)
    delta = d0 * np.exp(-(r / r0) ** 2)
    phi = C.spherical_potential(r, delta)
    x = r / r0
    mass = d0 * r0**3 * (np.sqrt(np.pi) / 4 * erf(x) - x * np.exp(-x**2) / 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        analytic = -C.LCDM.poisson_factor * (np.where(r > 0, mass / r, 0) + d0 * r0**2 / 2 * np.exp(-x**2))
    for rr in (0.0, 50.0, 150.0, 500.0):
        i = int(np.argmin(np.abs(r - rr)))
        check(f"Phi_0(r = {rr:.0f}) of a Gaussian void", float(phi[i]), float(analytic[i]), 2e-3)
    check("a void is a potential hill (Phi > 0)", float(np.sign(phi[0])), 1.0, 0, absolute=True)


def test_compensated_profile():
    print("Compensated (Finelli) void profile")
    r = np.linspace(0, 3000, 300001)
    d = C.finelli_profile(r, 0.25, 195.0)
    mass = np.trapezoid(d * r**2, r)
    check("total mass deficit integrates to zero", float(mass), 0.0, 1e-3 * 195.0**3, absolute=True)


def test_supervoid_isw():
    print("Supervoid ISW against the literature")
    th = np.array([0.0])
    r = np.linspace(0, 1500, 30001)
    dz = C.finelli_profile(r, 0.25, 195.0)  # Finelli et al. quote the depth today
    t0 = float(C.spherical_isw_profile(th, 0.155, r, dz)[0])
    # Kovacs et al. 2022: "the expected central ISW imprint is of order -20 uK"
    check("Eridanus supervoid central ISW (Finelli model) ~ -20 uK", t0, -20.0, 0.15)


def test_filters():
    print("Compensated top hat and forward model")
    import healpy as hp

    import stacking as S

    nside = 128
    mask = np.ones(hp.nside2npix(nside), bool)
    masks = {512: np.ones(hp.nside2npix(512), bool), 256: np.ones(hp.nside2npix(256), bool),
             128: mask}
    bank = S.FilterBank.build([30.0, 100.0], [40.0, -50.0], [4.0, 6.0], masks)
    const = {ns: np.full(hp.nside2npix(ns), 7.0, np.float32) for ns in masks}
    v = bank.apply(const)
    check("CTH of a constant map is zero", float(np.max(np.abs(v))), 0.0, 1e-5, absolute=True)
    # a uniform disc of +10 uK exactly matching the inner radius
    disc = {}
    for ns in masks:
        m = np.zeros(hp.nside2npix(ns), np.float32)
        m[hp.query_disc(ns, hp.ang2vec(30.0, 40.0, lonlat=True), np.radians(4.0))] = 10.0
        disc[ns] = m
    check("CTH of a matching +10 uK disc is +10 uK", float(bank.apply(disc)[0]), 10.0, 1e-6)

    import crosscorr as X

    ell = np.arange(X.LMAX + 1)
    cl = 1e-3 / (ell + 10.0) ** 2
    full = np.ones(hp.nside2npix(X.CC_NSIDE))
    fw = X.forward(cl, full, full)
    # the 600-node Gauss-Legendre projection is exact to ~1.5e-6
    check("forward model with no mask returns the input C_l (l=30)", float(fw[30]), float(cl[30]), 1e-4)


def test_cached():
    print("Checks on cached products (skipped if absent)")
    from config import PREREG_PATH, PREREGISTRATION, prereg_hash

    if PREREG_PATH.exists():
        import json

        frozen = json.loads(PREREG_PATH.read_text())
        for k, v in frozen.items():
            check(f"pre-registration '{k}' unchanged since {v['frozen_utc'][:16]}",
                  float(v["sha256"] == prereg_hash(PREREGISTRATION[k])), 1.0, 0, absolute=True)
    try:
        import boss
        from astropy.io import fits

        poly = fits.open(boss.product("boss_mask_north").path)[1].data
        p = boss._mangle_simple_pixel(poly["RAMID"], poly["DECMID"], 6)
        check("mangle simple pixel of BOSS polygon midpoints (fraction right)",
              float(np.mean(p == poly["PIXEL"])), 1.0, 5e-4)
    except FileNotFoundError:
        print("  [skip] BOSS mask not downloaded")
    try:
        import cmb
        import stacking as S

        d = cmb.prepare_maps()
        m = cmb.sim_map(5)
        slow, _ = S.multires(np.where(d["mask"], m, 0).astype(np.float32), d["mask"])
        fast = S.FastDegrader(d["mask"])(m)
        check("fast degrader equals healpy ud_grade (Nside 128)",
              float(np.max(np.abs(fast[128] - slow[128]))), 0.0, 1e-3, absolute=True)
    except FileNotFoundError:
        print("  [skip] Planck maps not downloaded")


def main() -> int:
    for t in (test_background, test_growth_vs_camb, test_spherical_potential,
              test_compensated_profile, test_supervoid_isw, test_filters, test_cached):
        t()
    print(f"\n{len(FAILURES)} failure(s)" + (": " + ", ".join(FAILURES) if FAILURES else ""))
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
