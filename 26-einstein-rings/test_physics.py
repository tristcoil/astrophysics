"""Check every formula in lensing.py against an answer known independently.

    python test_physics.py

Each check prints PASS or FAIL with the numbers compared. The script exits
non-zero if any check fails. If ``lenstronomy`` happens to be importable, the
SIE deflections are also compared with its implementation (the project itself
does not depend on it).
"""

from __future__ import annotations

import sys

import numpy as np
from scipy import integrate, special

import lensing as L

RESULTS: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, bool(ok), detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else ""))


def close(a, b, rtol=1e-6, atol=0.0):
    return bool(np.all(np.isclose(a, b, rtol=rtol, atol=atol)))


# --------------------------------------------------------------------------

def test_einstein_radius_and_mass():
    print("Einstein radius, critical density, mass")
    # A 250 km/s isothermal sphere with D_ls/D_s = 0.6:
    # theta_E = 4 pi (250/c)^2 * 0.6 rad = 1.0815"
    theta = L.theta_e_sis(250.0, 1.0, 0.6)
    expect = 4 * np.pi * (250 / 299792.458) ** 2 * 0.6 * (180 * 3600 / np.pi)
    check("SIS Einstein radius for 250 km/s, D_ls/D_s = 0.6", close(theta, expect, 1e-9),
          f"{theta:.5f}\" vs {expect:.5f}\"")
    check("sigma_SIE inverts theta_E", close(L.sigma_sie(theta, 1.0, 0.6), 250.0, 1e-12))

    # For an SIS the projected mass is M(<R) = pi sigma^2 R / G; the lensing mass
    # pi Sigma_cr R_E^2 must agree at R = R_E.
    d = L.distances(0.2, 0.6, 70.0, 0.3)
    th = L.theta_e_sis(260.0, d["d_s"], d["d_ls"])
    m_lens = L.einstein_mass(th, d["d_l"], d["d_s"], d["d_ls"])
    r_e_kpc = th * L.kpc_per_arcsec(d["d_l"])
    m_sis = np.pi * 260.0**2 * r_e_kpc / L.G_KPC
    check("pi Sigma_cr R_E^2 = pi sigma^2 R_E / G for an SIS", close(m_lens, m_sis, 1e-9),
          f"{m_lens:.4e} vs {m_sis:.4e} Msun")

    # Flat universe: D_ls = (D_M(z_s) - D_M(z_l)) / (1 + z_s).
    cosmo = L.cosmology(70.0, 0.3)
    dm_l = cosmo.comoving_transverse_distance(0.2).value
    dm_s = cosmo.comoving_transverse_distance(0.6).value
    check("D_ls from comoving distances", close(d["d_ls"], (dm_s - dm_l) / 1.6, 1e-9),
          f"{d['d_ls']:.2f} Mpc")
    # Sigma_cr in familiar units: c^2/(4 pi G) = 1.663e18 Msun/Mpc.
    k = L.C_KMS**2 / (4 * np.pi * L.G_KPC) * 1e3
    check("c^2 / 4 pi G = 1.663e18 Msun/Mpc", close(k, 1.6625e18, 1e-3), f"{k:.4e}")

    # Jorgensen: sigma_e8 / sigma_e2 = 4^0.04 (Bolton et al. 2008b: log ratio 0.0241)
    r = L.jorgensen_sigma(250.0, 1.5, 2.0, 1 / 8) / L.jorgensen_sigma(250.0, 1.5, 2.0, 1 / 2)
    check("Jorgensen sigma_e8/sigma_e2 = 10^0.0241", close(np.log10(r), 0.0241, 1e-2),
          f"log ratio {np.log10(r):.4f}")


def test_sie():
    print("Singular isothermal ellipsoid")
    rng = np.random.default_rng(1)
    x, y = rng.uniform(-3, 3, 200), rng.uniform(-3, 3, 200)
    ax, ay = L.sie_deflection(x, y, 1.2, 0.999999)
    r = np.hypot(x, y)
    check("q -> 1 gives the isothermal sphere, alpha = b r_hat",
          close(ax, 1.2 * x / r, 1e-4) and close(ay, 1.2 * y / r, 1e-4))

    # Divergence of the deflection = 2 kappa; curl = 0 (it is a gradient).
    b, q, h = 1.3, 0.6, 1e-5
    pts = [(0.7, 0.4), (-1.1, 0.9), (0.2, -1.6), (2.0, 0.1)]
    div_err, curl_err = [], []
    for px, py in pts:
        axp, _ = L.sie_deflection(px + h, py, b, q)
        axm, _ = L.sie_deflection(px - h, py, b, q)
        _, ayp = L.sie_deflection(px, py + h, b, q)
        _, aym = L.sie_deflection(px, py - h, b, q)
        div = (axp - axm) / (2 * h) + (ayp - aym) / (2 * h)
        div_err.append(abs(div - 2 * L.sie_convergence(px, py, b, q)))
        _, ay_xp = L.sie_deflection(px + h, py, b, q)
        _, ay_xm = L.sie_deflection(px - h, py, b, q)
        ax_yp, _ = L.sie_deflection(px, py + h, b, q)
        ax_ym, _ = L.sie_deflection(px, py - h, b, q)
        curl_err.append(abs((ay_xp - ay_xm) / (2 * h) - (ax_yp - ax_ym) / (2 * h)))
    check("div alpha = 2 kappa", max(div_err) < 1e-6, f"max error {max(div_err):.1e}")
    check("curl alpha = 0", max(curl_err) < 1e-6, f"max error {max(curl_err):.1e}")

    # Gauss: the flux of alpha through the critical curve equals 2 * mass / Sigma_cr.
    # Intermediate-axis normalisation: mass inside the critical curve = pi b^2.
    t = np.linspace(0, 2 * np.pi, 20001)
    xe, ye = b / np.sqrt(q) * np.cos(t), b * np.sqrt(q) * np.sin(t)
    dx, dy = np.gradient(xe, t), np.gradient(ye, t)
    ax, ay = L.sie_deflection(xe, ye, b, q)
    flux = np.trapezoid(ax * dy - ay * dx, t)
    check("mass inside the critical curve is pi b^2 (any q)", close(flux / 2, np.pi * b**2, 1e-5),
          f"{flux / 2:.6f} vs {np.pi * b**2:.6f}")

    # The critical curve is where det(A) = 0.
    xc, yc = L.sie_critical_curve(b, q, 0.0)
    dets = []
    for px, py in zip(xc[::90], yc[::90]):
        axp, ayp = L.sie_deflection(px + h, py, b, q)
        axm, aym = L.sie_deflection(px - h, py, b, q)
        axq, ayq = L.sie_deflection(px, py + h, b, q)
        axn, ayn = L.sie_deflection(px, py - h, b, q)
        a11 = 1 - (axp - axm) / (2 * h)
        a22 = 1 - (ayq - ayn) / (2 * h)
        a12 = -(axq - axn) / (2 * h)
        a21 = -(ayp - aym) / (2 * h)
        dets.append(abs(a11 * a22 - a12 * a21))
    check("det(A) = 0 on the critical curve", max(dets) < 1e-6, f"max |det| {max(dets):.1e}")

    # Rotating the lens rotates the deflection field.
    phi = 0.7
    xr, yr = L.rotate(np.array([0.9]), np.array([-0.3]), -phi)  # point in the rotated frame
    axs, ays = L.sie_deflection_sky(xr, yr, b, q, phi)
    ax0, ay0 = L.sie_deflection(np.array([0.9]), np.array([-0.3]), b, q)
    axr, ayr = L.rotate(axs, ays, phi)
    check("rotating the lens rotates its deflection", close(axr, ax0, 1e-10) and close(ayr, ay0, 1e-10))

    try:  # optional cross-check against an independent implementation
        from lenstronomy.LensModel.lens_model import LensModel
        from lenstronomy.Util.param_util import phi_q2_ellipticity
    except ImportError:
        print("  (lenstronomy not installed: independent SIE cross-check skipped)")
        return
    lm = LensModel(["SIE"])
    phi_l, q_l = 0.4, 0.55
    e1, e2 = phi_q2_ellipticity(phi_l, q_l)
    xs, ys = rng.uniform(-2.5, 2.5, 500), rng.uniform(-2.5, 2.5, 500)
    ax_l, ay_l = lm.alpha(xs, ys, [{"theta_E": 1.1, "e1": e1, "e2": e2, "center_x": 0, "center_y": 0}])
    ax_o, ay_o = L.sie_deflection_sky(xs, ys, 1.1, q_l, phi_l)
    err = np.max(np.hypot(ax_l - ax_o, ay_l - ay_o))
    check("SIE deflection matches lenstronomy (theta_E = b, intermediate axis)", err < 1e-6,
          f"max |delta alpha| {err:.1e}\"")


def test_profiles():
    print("Light and mass profiles")
    check("de Vaucouleurs b_4 = 7.669", close(L.sersic_b(4.0), 7.66925, 1e-5),
          f"{float(L.sersic_b(4.0)):.5f}")
    check("half the light inside R_e", close(L.sersic_light_fraction(2.0, 2.0), 0.5, 1e-12))
    f2 = L.sersic_light_fraction(1.0, 2.0)
    check("32.0% of de Vaucouleurs light inside R_e/2 (Bolton et al. 2008b)",
          close(f2, 0.320, 2e-3), f"{f2:.4f}")

    a = 1.0
    half = L.hernquist_mass2d(1.8153 * a, a)
    check("Hernquist: half the projected mass inside R = 1.8153 a", close(half, 0.5, 1e-4),
          f"{half:.5f}")
    num = integrate.quad(lambda R: 2 * np.pi * R * L.hernquist_surface(R, a), 0, 2.3)[0]
    check("Hernquist: projected mass = integral of surface density",
          close(num, L.hernquist_mass2d(2.3, a), 1e-6))
    check("Hernquist: M_3D(a) = 1/4", close(L.hernquist_mass3d(a, a), 0.25, 1e-12))

    check("isothermal: M_3D(R)/M_2D(R) = 2/pi",
          close(L.powerlaw_mass3d(1.0, 2.0, 1.0), 2 / np.pi, 1e-12))
    # gamma = 2.3: projected mass of the 3D profile, integrated numerically.
    g = 2.3
    k3 = L.powerlaw_mass3d(1.0, g, 1.0) * (3 - g) / (4 * np.pi)  # = rho0 r0^g

    def sigma_pl(R):
        return 2 * integrate.quad(lambda z: k3 * (R * R + z * z) ** (-g / 2), 0, np.inf)[0]
    m2 = integrate.quad(lambda R: 2 * np.pi * R * sigma_pl(R), 0, 1.0)[0]
    check("power law gamma = 2.3: projected mass inside theta_E = 1", close(m2, 1.0, 1e-5),
          f"{m2:.6f}")

    rs = 3.0
    def sigma_nfw(R):
        return 2 * integrate.quad(lambda z: 1 / (4 * np.pi * (np.hypot(R, z) / rs)
                                                 * (1 + np.hypot(R, z) / rs) ** 2 * rs**3),
                                  0, np.inf, limit=200)[0]
    for R in (0.5, 3.0, 7.0):
        m_num = integrate.quad(lambda x: 2 * np.pi * x * sigma_nfw(x), 0, R, limit=200)[0]
        m_ana = L.nfw_mass2d(np.array([R]), rs)[0]
        check(f"NFW projected mass at R = {R/rs:.2f} r_s", close(m_num, m_ana, 1e-4),
              f"{m_num:.5f} vs {m_ana:.5f}")


def test_jeans():
    print("Jeans equation and the fibre")
    w = L.fibre_weight(np.array([0.5, 1.4, 1.6, 3.0]), 1.5, 0.05)
    check("sharp seeing: the fibre takes all light inside 1.5\" and none outside",
          close(w, [1, 1, 0, 0], atol=1e-3), str(np.round(w, 4)))
    R = np.linspace(0, 8, 4000)
    total = np.trapezoid(L.fibre_weight(R, 1.5, 1.5) * 2 * np.pi * R, R)
    check("seeing moves light around but keeps the fibre's area", close(total, np.pi * 1.5**2, 1e-3),
          f"{total:.4f} vs {np.pi * 2.25:.4f}")

    # Hernquist's own isotropic sigma_r(r) (Hernquist 1990, eq. 10), G M = a = 1.
    jm = L.JeansModel(r_eff=1.8153, aperture_radius=1.5, seeing_fwhm=1.5, n_r=3000)

    def v2_hern(r):
        return L.hernquist_mass3d(r, 1.0) / r
    nsr2 = jm.nu_sigma_r2(v2_hern)
    r = np.array([0.1, 0.5, 1.0, 3.0])
    sr2 = np.interp(r, jm.r, nsr2) / L.hernquist_density(r, 1.0)
    exact = (1 / 12) * (12 * r * (r + 1) ** 3 * np.log((r + 1) / r)
                        - r / (r + 1) * (25 + 52 * r + 42 * r**2 + 12 * r**3))
    check("self-consistent Hernquist sigma_r(r) (Hernquist 1990, eq. 10)", close(sr2, exact, 2e-3),
          f"max rel. error {np.max(np.abs(sr2 / exact - 1)):.1e}")

    # Power-law tracer nu ~ r^-3.5 in a flat rotation curve v_c = 300 km/s:
    # sigma_r^2 = v_c^2 / (delta - 2 beta), sigma_p^2 = sigma_r^2 (1 - beta (delta-1)/delta)
    # at every radius, so any fibre average returns it.
    delta, vc = 3.5, 300.0
    tracer = (lambda r: r ** (-delta),
              lambda R: R ** (1 - delta) * np.sqrt(np.pi) * special.gamma((delta - 1) / 2)
              / special.gamma(delta / 2))
    for beta in (0.0, 0.3, -0.3):
        jm = L.JeansModel(r_eff=2.0, aperture_radius=1.5, seeing_fwhm=1.5, beta=beta,
                          tracer=tracer, n_r=2000, n_u=600)
        s = jm.sigma_fibre(lambda r: np.full_like(r, vc**2))
        exact = vc * np.sqrt((1 - beta * (delta - 1) / delta) / (delta - 2 * beta))
        check(f"power-law tracer in a flat rotation curve, beta = {beta:+.1f}",
              close(s, exact, 3e-3), f"{s:.2f} vs {exact:.2f} km/s")

    v2 = L.powerlaw_v2(2.0, 1.2, 250.0)
    check("isothermal: G M(r)/r = 2 sigma_SIS^2 at every radius",
          close(v2(np.array([0.1, 1.0, 10.0])), 2 * 250.0**2, 1e-12))
    # two-component model with all mass in a Hernquist profile: projected mass
    # inside theta_E must be the lensing mass.
    v2s = L.two_component_v2(1.0, 1.2, 250.0, 0.8, 30.0)
    m3 = v2s(np.array([1e4]))[0] * 1e4 / (np.pi * 250.0**2 * 1.2)  # total stellar mass / M_E
    check("two-component model: stars normalised to the lensing mass",
          close(m3 * L.hernquist_mass2d(1.2, 0.8), 1.0, 1e-3), f"{m3 * L.hernquist_mass2d(1.2, 0.8):.5f}")


def test_spectral_search():
    """Plant a background galaxy's lines at a known redshift in pure noise."""
    print("Blind search for a second redshift (injected lines)")
    import spectra
    rng = np.random.default_rng(3)
    loglam = np.arange(np.log10(3800), np.log10(9200), 1e-4)
    wave = 10**loglam
    z_true = 0.6123
    flux = rng.normal(0, 1.0, wave.size)
    sig = 1.5  # Angstrom, observed line width
    for lam, amp in ((3727.092, 3.0), (3729.875, 3.0), (4862.683, 2.0), (4960.295, 1.0), (5008.240, 3.0)):
        flux += amp * np.exp(-0.5 * ((wave - lam * (1 + z_true)) / sig) ** 2)
    # Lens at z = 0.25: its own possible emission lines (masked by the search) then
    # avoid the planted ones. At z = 0.2 its [O III] 5007 would sit on the planted [O II].
    spec = spectra.Spectrum("injected", wave, flux, np.ones_like(wave), np.zeros_like(wave),
                            0.25, 250.0, 10.0, 10.0)
    res = spectra.search_background(spec, 0.25)
    check("injected lines at z = 0.6123 found", abs(res.z_best - z_true) < 3e-4,
          f"found {res.z_best:.4f}, [O III] 5007 at {res.lines_at_best['[O III] 5007']:.0f} sigma")


def test_lens_fitter(quick: bool = False):
    """Make a fake Hubble image with a known lens, then fit it blind."""
    if quick:
        print("Lens fitter injection test skipped (--quick)")
        return
    print("Lens fitter (injected lens, about a minute)")
    import lensmodel as LM
    n = 160
    idx = (np.arange(n) - n // 2) * 0.04
    x, y = np.meshgrid(idx, idx)
    truth = dict(b=1.21, q_m=0.78, phi_m=0.6, xs=0.12, ys=-0.05, r_s=0.12, n_s=1.2, q_s=0.7, phi_s=1.9)
    var = np.full((n, n), 0.004**2)
    cut = LM.Cutout("injected", np.zeros((n, n)), var, x, y, np.pi / 2, 1, 0.004, 2000.0, (0, 0))
    fitter = LM.LensFitter(cut)
    # lens galaxy: a de Vaucouleurs-like profile; source lensed by the true SIE
    r_ell = np.sqrt(0.85 * x**2 + y**2 / 0.85)
    galaxy = 3.0 * np.exp(-7.669 * ((np.maximum(r_ell, 0.01) / 1.6) ** 0.25 - 1)) * 0.01
    from scipy.signal import fftconvolve
    galaxy = fftconvolve(galaxy, fitter.psf, mode="same")
    fitter.light_shape = (0.0, 0.0, 0.85, 0.0)
    ring = fitter.lensed_source(**truth)
    rng = np.random.default_rng(7)
    sci = galaxy + 0.25 * ring / ring.max() + rng.normal(0, 0.004, (n, n))
    cut.sci[:] = sci
    fitter = LM.LensFitter(cut)
    fitter.fit_light()
    res = fitter.fit()
    check("Einstein radius of an injected lens recovered to 1%", abs(res["b"] / truth["b"] - 1) < 0.01,
          f"{res['b']:.3f}\" vs {truth['b']}\"")
    check("its flattening recovered to 0.05", abs(res["q"] - truth["q_m"]) < 0.05,
          f"q = {res['q']:.3f} vs {truth['q_m']}")


def main() -> int:
    quick = "--quick" in sys.argv
    test_einstein_radius_and_mass()
    test_sie()
    test_profiles()
    test_jeans()
    test_spectral_search()
    test_lens_fitter(quick)
    n_fail = sum(not ok for _, ok, _ in RESULTS)
    print(f"\n{len(RESULTS) - n_fail} of {len(RESULTS)} checks passed.")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
