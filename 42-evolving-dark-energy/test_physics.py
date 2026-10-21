"""Self-checks: our physics against CAMB, closed forms and brute force.

    python test_physics.py

Each check prints PASS/FAIL with the measured deviation and its tolerance.
CAMB (in the shared venv) is the reference for the expansion history and the
sound horizons; our likelihoods are checked against direct transcriptions of
Cobaya's formulas, and the C pair counter against a brute-force count.
"""

from __future__ import annotations

import sys

import numpy as np

import background as bg
import early_universe as eu

RESULTS = []


def check(name: str, value: float, tol: float) -> None:
    ok = bool(np.isfinite(value) and abs(value) <= tol)
    RESULTS.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<62} {value: .2e}  (tol {tol:.0e})")


# Three universes: Planck LCDM, DESI's best fit for DESI+CMB+Pantheon+, a strong w0wa.
CASES = [(0.6736, 0.02237, 0.1200, -1.0, 0.0),
         (0.67495, 0.022314, 0.118672, -0.8532, -0.5216),
         (0.637, 0.02223, 0.1199, -0.43, -1.72)]


def test_against_camb() -> None:
    print("Expansion history and rulers vs CAMB 2.0.4 (same settings as DESI)")
    z = np.array([0.0, 0.1, 0.5, 1.0, 2.33, 10.0, 100.0, 1000.0, 1090.0, 1e4, 1e6])
    zt = np.array([0.295, 0.51, 0.706, 0.934, 1.321, 1.484, 2.33])
    sh = bg.SoundHorizonGrid()
    for h, ob, oc, w0, wa in CASES:
        d = eu.camb_derived(ob, oc, 100 * h, w0, wa)
        res = d["_results"]
        c = bg.Cosmo(h=h, ombh2=ob, omch2=oc, w0=w0, wa=wa)
        tag = f"(w0={w0:+.2f}, wa={wa:+.2f})"
        check(f"H(z), z = 0 to 1e6 {tag}",
              np.max(np.abs(bg.hubble(z, c)[0] / res.hubble_parameter(z) - 1)), 1e-6)
        g = bg.DistanceGrid(zt, z_max=1400.0)
        check(f"D_M at the DESI redshifts {tag}",
              np.max(np.abs(g.comoving(c)[0] / res.comoving_radial_distance(zt) - 1)), 2e-6)
        dms = g.comoving_at(np.array([d["zstar"]]), c)[0]
        check(f"D_M(z*) {tag}", dms / d["DMstar"] - 1, 2e-6)
        rs = sh.r_s(np.array([d["zstar"], d["zdrag"]]), c.take(np.array([0, 0])))
        check(f"r_* at CAMB's z* {tag}", rs[0] / d["rstar"] - 1, 1e-5)
        check(f"r_d at CAMB's z_d {tag}", rs[1] / d["rdrag"] - 1, 1e-5)
        check(f"theta_* {tag}", (rs[0] / dms) / d["thetastar"] - 1, 1e-5)


def test_recombination_grid() -> None:
    print("Interpolated recombination redshifts and rulers vs fresh CAMB runs")
    E = eu.EarlyUniverse()
    rng = np.random.default_rng(5)
    dz = {"zstar": [], "zdrag": []}
    dr = {"rstar": [], "rdrag": []}
    for ob, oc in zip(rng.uniform(0.0205, 0.0245, 12), rng.uniform(0.10, 0.14, 12)):
        d = eu.camb_derived(ob, oc)
        for k in dz:
            dz[k].append(E(k, ob, oc)[0] - d[k])
        for k in dr:
            dr[k].append(E(k, ob, oc)[0] / d[k] - 1)
    check("z* (absolute)", np.max(np.abs(dz["zstar"])), 0.01)
    check("z_d (absolute)", np.max(np.abs(dz["zdrag"])), 0.02)
    check("r_* (relative)", np.max(np.abs(dr["rstar"])), 1e-5)
    check("r_d (relative)", np.max(np.abs(dr["rdrag"])), 2e-5)
    a, b = eu.camb_derived(0.0222, 0.12, 55.0), eu.camb_derived(0.0222, 0.12, 80.0)
    check("r_d independent of H0 at fixed omega_b, omega_c", b["rdrag"] / a["rdrag"] - 1, 1e-6)


def test_neutrinos() -> None:
    print("Massive neutrino density")
    check("relativistic limit: rho -> N_eff/3 of a massless species",
          bg._fd_energy_ratio(np.array([1e-6]))[0] - 1, 1e-10)
    y = np.array([500.0])
    # m n (1 + <p^2>/2m^2): int x^4/(e^x+1) = (15/16) 24 zeta(5), int x^2/(e^x+1) = (3/2) zeta(3)
    i4 = 15.0 / 16.0 * 24.0 * 1.0369277551433699
    nr = y * bg._FD_NR / bg._FD_F0 * (1 + i4 / bg._FD_NR / (2 * y**2))
    check("non-relativistic limit: rho -> m n (1 + <p^2>/2m^2)",
          bg._fd_energy_ratio(y)[0] / nr[0] - 1, 1e-8)
    check("today's density equals CAMB's omnuh2",
          bg.omega_nu_of_z(np.array([0.0]))[0] / bg.OMEGA_NU - 1, 1e-4)


def test_sn_likelihood() -> None:
    print("Supernova chi^2 vs a literal transcription of Cobaya's formula")
    import fetch_data as fd
    from likelihoods import SNLikelihood

    for name in ("PantheonPlus", "DES-Dovekie"):
        s = fd.load_sn(name)
        lik = SNLikelihood.from_sample(s)
        inv = np.linalg.inv(s["cov"]) if "cov" in s else s["invcov"]
        d = inv @ np.ones(len(s["mu"]))
        inv_m = inv - np.outer(d, d) / d.sum()  # Cobaya's _marginalize_abs_mag
        c = bg.Cosmo(h=0.68, ombh2=0.0224, omch2=0.12, w0=-0.9, wa=-0.3)
        g = bg.DistanceGrid(lik.z_hd, z_max=float(lik.z_hd.max()))
        dm = g.comoving(c)
        mine, _ = lik.chi2(dm)
        r = s["mu"] - lik.model_mu(dm)[0]
        check(f"{name}: profiled chi^2", mine[0] / (r @ inv_m @ r) - 1, 1e-9)


def test_bao_profile() -> None:
    print("BAO: analytic profiling of h*r_d vs brute-force minimisation")
    import fetch_data as fd
    from likelihoods import BAOLikelihood
    from scipy.optimize import minimize_scalar

    bao = BAOLikelihood(**fd.load_bao("dr2"))
    c = bg.Cosmo(h=0.6736, ombh2=0.02237, omch2=0.12, w0=-1.0, wa=0.0)
    g = bg.DistanceGrid(bao.unique_z, z_max=3.0)
    un = bao.model_unscaled(g.comoving(c), bg.hubble(bao.unique_z, c))
    chi2, rd = bao.chi2_profile_scale(un)
    best = minimize_scalar(lambda x: bao.chi2(un / x)[0], bounds=(100, 200), method="bounded",
                           options={"xatol": 1e-9})
    check("profiled chi^2", chi2[0] - best.fun, 1e-6)
    check("profiled r_d", rd[0] / best.x - 1, 1e-6)


def test_pair_counter() -> None:
    print("C pair counter vs brute force (2,000 + 1,500 points)")
    import bao_measurement as B

    rng = np.random.default_rng(0)

    def cat(n):
        return {"pos": rng.uniform(-250, 250, (n, 3)) + np.array([1400.0, 150.0, 250.0]),
                "w": rng.uniform(0.5, 1.5, n), "reg": rng.integers(0, 4, n).astype(np.int32),
                "z": rng.uniform(0.4, 0.6, n)}

    a, b = cat(2000), cat(1500)
    for auto in (True, False):
        c2 = a if auto else b
        dv = c2["pos"][None, :, :] - a["pos"][:, None, :]
        s = np.sqrt((dv**2).sum(-1))
        los = a["pos"][:, None, :] + c2["pos"][None, :, :]
        mu = np.abs((dv * los).sum(-1)) / (s * np.sqrt((los**2).sum(-1)) + 1e-300)
        ww = a["w"][:, None] * c2["w"][None, :]
        m = (s < B.S_MAX) & (s > 0)
        if auto:
            m &= np.triu(np.ones_like(m), 1).astype(bool)
        H, _, _ = np.histogram2d(s[m], mu[m], bins=[B.NS, B.NMU], range=[[0, B.S_MAX], [0, 1]],
                                 weights=ww[m])
        res = B.pair_counts(a, None if auto else b, nreg=4)
        check(f"{'auto' if auto else 'cross'}-pair counts", np.abs(res["total"] - H).max() / H.max(),
              1e-12)


def test_reconstruction() -> None:
    print("Reconstruction solver on analytic fields")
    import scipy.fft as sfft

    import reconstruction as R

    rng = np.random.default_rng(1)
    pts = rng.uniform(0, 1, (1000, 3)) * 1024 + np.array([1e6, 0, 0])  # far away: r^ ~ x
    mesh = R.Mesh(pts, cell=8.0, pad=0.0)
    x = mesh.origin[0] + (np.arange(mesh.shape[0]) + 0.5) * mesh.cell
    k0, A = 2 * np.pi * 4 / (mesh.shape[0] * mesh.cell), 0.1
    delta = np.broadcast_to((A * np.sin(k0 * x))[:, None, None], mesh.shape).astype(np.float32)
    dk = sfft.rfftn(delta)
    kvec = mesh.kvectors()
    psi = R.displacement(mesh, dk, kvec, workers=4)
    check("Zel'dovich displacement of a plane wave (del.Psi = -delta)",
          np.abs(psi[0][:, 3, 3] - (A / k0) * np.cos(k0 * x)).max() / (A / k0), 1e-5)
    f, b = R.GROWTH, R.BIAS
    rhat = mesh.unit_los()
    p = R.displacement(mesh, dk / b, kvec, 4)
    for _ in range(R.NITER - 1):
        proj = p[0] * rhat[0] + p[1] * rhat[1] + p[2] * rhat[2]
        div_k = sum(1j * kc * sfft.rfftn(proj * rc) for kc, rc in zip(kvec, rhat))
        p = R.displacement(mesh, (dk + f * div_k) / b, kvec, 4)
    amp = A / (k0 * b * (1 + f / b))
    check("redshift-space term along the line of sight (amplitude 1/(1+f/b))",
          np.abs(p[0][:, 3, 3] - amp * np.cos(k0 * x)).max() / amp, 2e-3)
    pos = rng.uniform(100, 900, (100000, 3)) + np.array([1e6, 0, 0])
    m2 = R.Mesh(pos, cell=8.0, pad=50.0)
    check("cloud-in-cell painting conserves weight",
          m2.paint(pos, np.ones(len(pos))).sum() / len(pos) - 1, 1e-6)


if __name__ == "__main__":
    test_neutrinos()
    test_against_camb()
    test_recombination_grid()
    test_sn_likelihood()
    test_bao_profile()
    test_pair_counter()
    test_reconstruction()
    n_fail = RESULTS.count(False)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
    sys.exit(1 if n_fail else 0)
