"""Checks of the physics and the numerics against results known exactly or from the literature.

    python test_physics.py        (each test prints 'ok'; the SymPy one takes ~1 min)
"""

from __future__ import annotations

import math

import mpmath as mp
import numpy as np

import colour
import geodesic
import infall
import kerr
import raymap as R

B_C = math.sqrt(27.0)


# ---------------------------------------------------------------------------
# The ray map
# ---------------------------------------------------------------------------

def _breaks(lo, hi, centre, width):
    """Quadrature breakpoints clustered geometrically around a narrow peak."""
    pts = {lo, hi}
    for k in range(0, 12):
        for sgn in (-1, 1):
            x = centre + sgn * width * mp.mpf(10) ** k
            if lo < x < hi:
                pts.add(x)
    if lo < centre < hi:
        pts.add(centre)
    return sorted(pts)


def _phi_mpmath(r, chi):
    """The orbit integral at 40 digits. Tanh-sinh quadrature copes with end-point singularities;
    the narrow peaks near the photon sphere get breakpoints scaled to their width."""
    mp.mp.dps = 40
    u_obs = mp.mpf(1) / r
    bt = mp.sqrt(2 / mp.mpf(r))
    c, s = mp.cos(chi), mp.sin(chi)
    E = 1 - bt * c
    b = r * s / E
    k = 1 / b**2
    f = lambda u: 1 / mp.sqrt(2 * u**3 - u**2 + k)
    third = mp.mpf(1) / 3
    if bt - c > 0:          # backward ray moves out the whole way: P(u) ~ delta + (u - 1/3)^2 near 1/3
        width = mp.sqrt(abs(k - mp.mpf(1) / 27)) + mp.mpf(10) ** -38
        return float(mp.re(mp.quad(f, _breaks(mp.mpf(0), u_obs, third, width))))
    u1, up, u3 = sorted([x.real for x in mp.polyroots([2, -1, 0, k], extraprec=80) if abs(x.imag) < 1e-25])
    # u = up - s^2 and P = 2 (u - u1)(u - up)(u - u3): the factor s cancels exactly, leaving
    # 2 / sqrt(2 (up - s^2 - u1)(u3 - up + s^2)), peaked at s = 0 with width sqrt(u3 - up)
    g = lambda s_: 2 / mp.sqrt(2 * (up - s_**2 - u1) * (u3 - up + s_**2))
    w = mp.sqrt(u3 - up)
    val = (mp.quad(g, _breaks(mp.mpf(0), mp.sqrt(up), mp.mpf(0), w))
           + mp.quad(g, _breaks(mp.mpf(0), mp.sqrt(up - u_obs), mp.mpf(0), w)))
    return float(mp.re(val))


def test_carlson_against_40_digit_quadrature():
    """Our double-precision closed form against 40-digit quadrature, from the shadow's edge
    (phi -> infinity: rays winding round the photon sphere) to straight overhead.

    Within ~1e-7 rad of the edge the answer is ill-conditioned: b - b_c, computed in double
    precision from the viewing angle, carries round-off that shifts phi by ~1e-8. For an
    observer exactly on the photon sphere (r = 3) the edge is a maximum of b(chi), so
    b_c - b is quadratically small there and the same round-off moves phi by ~1e-3;
    those two points are excluded (the renderer averages that region anyway)."""
    worst = {}
    for r in (1000.0, 50.0, 6.0, 3.5, 3.01, 3.0, 2.99, 2.5, 2.0, 1.0, 0.3, 0.05):
        chs = R.shadow_edge(r)
        for d in (1e-7, 1e-4, 1e-2, 0.3, 1.0):
            if d == 1e-7 and abs(r - 3.0) < 0.02:
                continue
            chi = min(chs + d * (math.pi - chs), math.pi - 1e-6)
            ours = float(R.deflection(r, np.cos([chi]))[0][0])
            ref = _phi_mpmath(r, chi)
            tol = 1e-7 if d == 1e-7 else 1e-8
            assert abs(ours - ref) < tol, (r, d, ours, ref)
            worst[d] = max(worst.get(d, 0.0), abs(ours - ref))
    return worst


def test_carlson_against_brute_force_geodesics():
    """Christoffel symbols from SymPy, momentum from the observer's own frame, DOP853 steps."""
    for r, deg in ((1000.0, 0.6), (50.0, 25.0), (6.0, 60.0), (6.0, 120.0), (3.2, 70.0),
                   (2.5, 100.0), (2.0, 80.0), (1.0, 120.0), (0.5, 95.0)):
        chi = math.radians(deg)
        esc, phi_num, _ = geodesic.trace(r, chi)
        phi = float(R.deflection(r, np.cos([chi]))[0][0])
        assert esc and abs(phi_num - phi) < 1e-8, (r, deg, phi_num, phi)
    # and a ray inside the shadow really does not come from the sky
    esc, _, _ = geodesic.trace(6.0, math.radians(10.0))
    assert not esc


def test_weak_field_deflection():
    """Light passing at impact parameter b >> M bends by 4M/b + (15 pi/4) M^2/b^2 + (128/3) M^3/b^3."""
    for b in (30.0, 100.0, 1000.0):
        k = 1.0 / b**2
        u1, u2, u3 = R._roots_real(np.array([k]))
        half = R._carlson_segment(u2, np.zeros(1), u1, u2, u3)[0]
        alpha = 2 * half - math.pi
        series = 4 / b + 15 * math.pi / 4 / b**2 + 128 / 3 / b**3
        assert abs(alpha - series) < 2e3 / b**4, (b, alpha, series)


def test_strong_deflection_limit():
    """Near the photon sphere alpha = -ln(b/b_c - 1) + b_bar, b_bar = -pi + 2 ln(6(2 - sqrt 3)) + ln 6
    = -0.40023 (Bozza 2002, PRD 66, 103001): light can circle the hole any number of times."""
    b_bar = -math.pi + 2 * math.log(6 * (2 - math.sqrt(3))) + math.log(6)
    for eps in (1e-5, 1e-7, 1e-9):
        b = B_C * (1 + eps)
        k = 1.0 / b**2
        u1, u2, u3 = R._roots_real(np.array([k]))
        alpha = 2 * R._carlson_segment(u2, np.zeros(1), u1, u2, u3)[0] - math.pi
        assert abs(alpha - (-math.log(eps) + b_bar)) < 20 * math.sqrt(eps) + 1e-6, (eps, alpha)


def test_shadow_sizes():
    # falling from rest at infinity, at the horizon: tan(chi/2) = 2 / sqrt 27 -> 42.103 degrees
    assert abs(R.shadow_edge(2.0) - 2 * math.atan(2 / B_C)) < 1e-12
    # hovering observers: sin(psi) = sqrt 27 sqrt(1 - 2/r) / r (Synge 1966), > 90 deg inside r = 3
    for r in (100.0, 6.0, 4.0, 3.0, 2.5, 2.1):
        tab = R.RayTable(r, R.velocity_rel_rain(r, "static"))
        s = B_C * math.sqrt(1 - 2 / r) / r
        expect = math.asin(min(s, 1.0)) if r >= 3 else math.pi - math.asin(s)
        assert abs(tab.chi_sh - expect) < 1e-9, (r, tab.chi_sh, expect)
    assert abs(math.degrees(R.RayTable(6.0, R.velocity_rel_rain(6.0, "static")).chi_sh) - 45.0) < 1e-8


def test_frequency_shifts():
    # a hovering observer sees the whole sky blueshifted by 1/sqrt(1 - 2/r), whatever the direction
    for r in (20.0, 3.0, 2.2):
        tab = R.RayTable(r, R.velocity_rel_rain(r, "static"))
        assert np.allclose(tab.g, 1 / math.sqrt(1 - 2 / r), rtol=1e-12)
    # the raindrop at the horizon: light from straight above arrives at half its frequency,
    # light at the shadow's edge at 1/(1 - cos 42.1 deg) = 3.875
    _, _, _, g = R.photon_constants(2.0, np.array([-1.0, math.cos(R.shadow_edge(2.0))]))
    assert abs(g[0] - 0.5) < 1e-12 and abs(g[1] - 1 / (1 - math.cos(R.shadow_edge(2.0)))) < 1e-12
    assert abs(g[1] - 3.875) < 1e-3


def test_observer_velocities():
    # a body dropped from rest is, at that moment, a hovering observer; far below its release
    # point it approaches the raindrop's speed
    assert abs(R.velocity_rel_rain(10.0, "drop", r0=10.0 + 1e-12) - math.sqrt(0.2)) < 1e-5
    v = R.velocity_rel_rain(2.0, "drop", r0=1e8)
    assert 0 < v < 1e-3
    # aberration round trip
    c = np.linspace(-1, 1, 7)
    cr, _ = R.boost_to_rain(c, 0.6)
    back, _ = R.boost_to_rain(cr, -0.6)
    assert np.allclose(back, c)


# ---------------------------------------------------------------------------
# Falling: clocks, tides, the friend's view
# ---------------------------------------------------------------------------

def test_proper_times():
    assert abs(infall.rain_tau_inside() - 4 / 3) < 1e-14
    # dropped from rest just above the horizon: the longest possible stay inside, pi M
    assert abs(infall.drop_tau_inside(2.0 + 1e-9) - math.pi) < 1e-3
    assert infall.drop_tau_inside(3.0) < math.pi and infall.drop_tau_inside(1e6) > 4 / 3
    assert abs(infall.drop_tau_inside(1e8) - 4 / 3) < 1e-3
    # Newton's free-fall time pi sqrt(r0^3 / 8), checked by quadrature of dtau/dr
    r0 = 37.0
    num, _ = __import__("scipy").integrate.quad(lambda x: 1 / math.sqrt(2 / x - 2 / r0), 0, r0)
    assert abs(num - infall.drop_tau_total(r0)) < 1e-8


def test_schwarzschild_time_of_a_dropped_body():
    """MTW eq. 31.10 against direct quadrature."""
    r0 = 12.0
    for eta in (0.3, 1.0, 1.6, 1.9):
        r, _ = infall.drop_r_tau(eta, r0)
        if r <= 2.0:
            continue
        assert abs(infall.drop_schwarzschild_t(eta, r0) - infall.drop_schwarzschild_t_numeric(float(r), r0)) < 1e-6


def test_raindrop_schwarzschild_time():
    """Closed form against quadrature of dt/dr = -sqrt(r/2) r/(r - 2)."""
    from scipy.integrate import quad
    num, _ = quad(lambda r: math.sqrt(r / 2) * r / (r - 2), 5.0, 40.0)
    assert abs((infall.rain_schwarzschild_t(5.0) - infall.rain_schwarzschild_t(40.0)) - num) < 1e-9


def test_friend_sees_exponential_redshift():
    """Late-time redshift of a faller grows as exp(t / 4M) (e.g. Frolov & Novikov 1998)."""
    tau, t, z1, r = infall.rain_received(1e4, 50.0, n=20000)
    assert abs(infall.late_time_efold(t, z1) - 4.0) < 0.01
    # and the faller's clock, as received, freezes at its horizon-crossing reading
    assert abs(tau[-1] - (infall.rain_tau_left(50.0) - 4 / 3)) < 1e-5


def test_kruskal_singularity_and_horizon():
    U, V = infall.kruskal_from_v_r(np.array([1.0, 5.0]), np.array([0.0, 2.0]))
    assert abs(U[0] * V[0] - 1.0) < 1e-12 and abs(U[1]) < 1e-12
    x, y = infall.penrose(U[0], V[0])
    assert abs(y - math.pi / 2) < 1e-12        # r = 0 is a horizontal line: a moment, not a place


def test_tidal_numbers():
    gm = 1.32712440041e20
    c = 299_792_458.0
    # 10 Msun: ~10^7 g across 2 m at the horizon (the atlas spec); Sgr A* (4e6): ~1e-4 g
    rs10 = 2 * gm * 10 / c**2
    assert 1e7 < infall.tidal_stretch_si(rs10, 10, 2.0, gm) / 9.80665 < 3e7
    rs = 2 * gm * 4e6 / c**2
    assert 5e-5 < infall.tidal_stretch_si(rs, 4e6, 2.0, gm) / 9.80665 < 3e-4
    m_crit = infall.critical_mass_for_horizon_stretch(2.0, 10 * 9.80665, gm, c)
    assert abs(infall.tidal_stretch_si(2 * gm * m_crit / c**2, m_crit, 2.0, gm) - 98.0665) < 1e-6


# ---------------------------------------------------------------------------
# Colour
# ---------------------------------------------------------------------------

def test_blackbody_shift_and_limits():
    lt = math.log10(5772.0)
    assert abs(float(colour.luminance_ratio(lt, 1.0)) - 1.0) < 1e-12
    # Rayleigh-Jeans: at enormous blueshift the visible brightness grows only like g
    r1, r2 = float(colour.luminance_ratio(lt, 1e4)), float(colour.luminance_ratio(lt, 2e4))
    assert abs(r2 / r1 - 2.0) < 0.01
    # the colour of an infinitely hot blackbody: about (148, 177, 255) in sRGB (Charity's table: 155, 176, 255)
    sw = colour.blackbody_srgb_swatch(1e9) * 255
    assert abs(sw[0] - 150) < 10 and abs(sw[1] - 178) < 6 and sw[2] > 250


# ---------------------------------------------------------------------------
# Kerr
# ---------------------------------------------------------------------------

def test_kerr_standard_numbers():
    assert abs(kerr.isco(0.0) - 6.0) < 1e-12
    assert abs(kerr.efficiency(0.0) - (1 - math.sqrt(8 / 9))) < 1e-12          # 5.72%
    assert abs(kerr.isco(0.998) - 1.2369) < 1e-3                                  # Thorne limit
    assert abs(kerr.efficiency(0.998) - 0.321) < 2e-3                             # Thorne 1974: 32%
    assert abs(kerr.efficiency(1 - 1e-12) - (1 - 1 / math.sqrt(3))) < 1e-3       # 42.3% as a -> 1
    assert abs(kerr.isco(0.0, prograde=False) - 6.0) < 1e-12
    assert abs(kerr.isco(1.0 - 1e-12, prograde=False) - 9.0) < 1e-3
    rp, rm = kerr.horizons(0.68)
    assert abs(rp * rm - 0.68**2) < 1e-12
    # inside the ergosphere nothing can stand still; at the horizon everything turns with Omega_H
    lo, hi = kerr.allowed_omega(np.array([1.9, 2.5]), 0.68)
    assert lo[0] > 0 and lo[1] < 0
    lo, hi = kerr.allowed_omega(np.array([rp + 1e-9]), 0.68)
    assert abs(lo[0] - 0.68 / (rp**2 + 0.68**2)) < 1e-4 and abs(hi[0] - lo[0]) < 1e-3


def test_kerr_plunge_reduces_to_schwarzschild():
    """For a = 0 a body with E = 1, L = 0 takes (2/3)(r^1.5 - 2^1.5)/sqrt 2 to reach the horizon."""
    out = kerr.plunge(1e-9, 1.0, 0.0, 20.0, 1.5)
    expect = math.sqrt(2) / 3 * (20.0**1.5 - 2.0**1.5)
    assert abs(out["tau_to_horizon"] - expect) < 2e-2


def test_symbolic_geometry():
    import symbolic
    res = symbolic.derive()
    for key in ("gp_equals_schwarzschild", "spatial_metric_flat", "raindrop_is_geodesic",
                "kretschmann_ok", "tidal_rain_ok", "tidal_drop_ok"):
        assert res[key] is True, (key, res)


if __name__ == "__main__":
    import sys
    only = sys.argv[1:]
    for name, fn in list(globals().items()):
        if name.startswith("test_") and (not only or any(o in name for o in only)):
            fn()
            print("ok", name, flush=True)
