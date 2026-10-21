"""Checks of the physics and the numerics against textbook results.

    python -m pytest -q test_physics.py      (or: python test_physics.py)

Each test isolates one ingredient: the numerical orbit model must reproduce the
analytic Lense-Thirring, de Sitter and J2 node rates; the frame rotation must
turn SP3 velocities consistently; the GP-B integrator must reproduce the
circular-orbit formulas; the pulsar formula must reproduce the published GR
rates; the tide model must respond to the Moon's 18.6-year cycle in the right
sense.
"""

from __future__ import annotations

import numpy as np

import config as C
import literature as L

MAS = np.pi / 180 / 3600e3


def _circular_state(a, inc, node, u=0.0):
    v = np.sqrt(C.GM_EARTH / a)
    r = a * np.array([np.cos(u), np.sin(u), 0.0])
    vv = v * np.array([-np.sin(u), np.cos(u), 0.0])
    ci, si, cO, sO = np.cos(inc), np.sin(inc), np.cos(node), np.sin(node)
    R = np.array([[cO, -sO * ci, sO * si], [sO, cO * ci, -cO * si], [0, si, ci]])
    return np.concatenate([R @ r, R @ vv])


def _node_change(model, y0, days, **kw):
    import frames as F
    T, Y = model.run(y0, [1.1, 0, 0, 0, 0], 0.0, days * 86400.0, stride=10, **kw)
    el = F.kepler_elements(Y[:, :3], Y[:, 3:6])
    node = np.unwrap(el["node"])
    # mean over the last and first revolution-ish windows (2 days) to beat short-period terms
    n = len(T)
    w = int(2 * 86400 / (T[1] - T[0]))
    return (node[-w:].mean() - node[:w].mean()), (T[-w:].mean() - T[:w].mean())


def test_lense_thirring_node_rate():
    """mu_LT = 1 minus mu_LT = 0 must turn a LAGEOS orbit's node by 2GJ/(c^2 a^3) per unit time."""
    import arcmodel as AM
    import framedrag as FD
    m = AM.build(60_400.0, 9.0, "lageos1", ocean=False, solid=False, pole=False)
    y0 = _circular_state(12_270e3, np.radians(109.84), np.radians(40.0))
    d1, dt = _node_change(m, y0, 8.0, mu_lt=1.0)
    d0, _ = _node_change(m, y0, 8.0, mu_lt=0.0)
    rate = (d1 - d0) / dt * C.JULIAN_YEAR_S / MAS
    expect = FD.lt_rate_mas_yr(12_270e3, 0.0)
    assert abs(rate / expect - 1) < 0.01, (rate, expect)


def test_de_sitter_node_rate():
    """The de Sitter term turns every orbit about the ecliptic pole: ~17.6 mas/yr on the node."""
    import arcmodel as AM
    m = AM.build(60_400.0, 9.0, "lares2", ocean=False, solid=False, pole=False)
    y0 = _circular_state(12_266e3, np.radians(70.16), np.radians(40.0))
    d1, dt = _node_change(m, y0, 8.0, mu_ds=1.0)
    d0, _ = _node_change(m, y0, 8.0, mu_ds=0.0)
    rate = (d1 - d0) / dt * C.JULIAN_YEAR_S / MAS
    # rate depends on the node's position relative to the ecliptic pole direction: 19.2 x (cos eps
    # + sin eps * cot i * cos(...)) -> between ~15 and ~20 mas/yr for these orbits
    assert 14.0 < rate < 21.5, rate


def test_j2_node_rate():
    """With a J2-only Earth (no other harmonics, no Sun or Moon, no radiation pressure) the
    numerical node drift matches the first-order formula -3/2 n J2 (R/a)^2 cos i / (1-e^2)^2
    evaluated with the trajectory's mean elements, to the size of second-order terms (~3e-4)."""
    import arcmodel as AM
    import framedrag as FD
    import frames as F
    m = AM.build(60_400.0, 9.0, "lageos1", ocean=False, solid=False, pole=False, srp=False)
    c20 = m.Cs[2, 0]
    m.Cs[:] = 0.0
    m.Ss[:] = 0.0
    m.Cs[0, 0] = 1.0
    m.Cs[2, 0] = c20
    m.gmb[:] = 0.0
    y0 = _circular_state(12_270e3, np.radians(109.84), np.radians(40.0))
    T, Y = m.run(y0, [0, 0, 0, 0, 0], 0.0, 8 * 86400.0, stride=1, mu_ds=0.0)
    # elements about the CIP (the model's J2 axis), not the GCRS pole
    import erfa
    x, yy = erfa.xy06(2_400_000.5, 60_404.0)
    Q = erfa.c2ixys(x, yy, erfa.s06(2_400_000.5, 60_404.0, x, yy))
    el = F.kepler_elements(Y[:, :3] @ Q.T, Y[:, 3:6] @ Q.T)
    node = np.unwrap(el["node"])
    rate = np.polyfit(T, node, 1)[0] * C.JULIAN_YEAR_S / MAS
    expect = FD.j2_node_rate_mas_yr(el["a"].mean(), el["e"].mean(), el["i"].mean(), j2=-c20)
    assert abs(rate / expect - 1) < 1e-3, (rate / 3.6e6, expect / 3.6e6)


def test_frame_velocity_consistency():
    """GCRS velocities from rotating SP3 velocities agree with differentiated GCRS positions."""
    import frames as F
    t = 60_400.0 + np.arange(0, 3600, 1.0) / 86400.0
    # a synthetic Earth-fixed track: a point fixed in the GCRS seen from the rotating Earth
    R = F.itrs_to_gcrs_matrices(t)
    r_g = np.tile([7.0e6, 1.0e6, 2.0e6], (t.size, 1))
    r_i = np.einsum("nji,nj->ni", R, r_g)
    v_i = np.gradient(r_i, 1.0, axis=0)
    _, v_g = F.itrs_to_gcrs(t, r_i, v_i)
    assert np.max(np.linalg.norm(v_g[5:-5], axis=1)) < 1e-3     # a fixed point: ~zero velocity


def test_gpb_circular_formulas():
    """With a point-mass Earth the integrated geodetic rate equals 3/2 (GM)^1.5 / (c^2 a^2.5)."""
    import gpb
    p = gpb.predict(n_days=10.0, step_s=20.0, lunisolar=False)
    r = p["mean_r_km"] * 1e3
    circ = 1.5 * C.GM_EARTH ** 1.5 / (C.C_LIGHT ** 2 * r ** 2.5) * C.JULIAN_YEAR_S / MAS
    assert abs(-p["geodetic_point_mass"]["ns_mas_yr"] / circ - 1) < 2e-3, (p["geodetic_point_mass"], circ)
    # frame dragging averaged over a polar orbit: G J / (2 c^2 a^3), projected by cos(dec)
    fd = C.G_NEWTON * C.EARTH_SPIN_PER_MASS * C.GM_EARTH / C.G_NEWTON / (2 * C.C_LIGHT ** 2 * r ** 3)
    fd = fd * C.JULIAN_YEAR_S / MAS * np.cos(np.radians(L.IM_PEG["dec_deg"]))
    assert abs(-p["frame_dragging"]["we_mas_yr"] / fd - 1) < 0.01, (p["frame_dragging"], fd)


def test_solar_geodetic_projection():
    """19.2 mas/yr about the ecliptic pole, projected at IM Peg: +7.3 NS, -16.2 WE (Everitt 2011)."""
    import gpb
    gs = gpb.guide_star_unit()
    ns, we = gpb.ns_we_basis(gs)
    eps = np.radians(23.4393)
    pole = np.array([0.0, -np.sin(eps), np.cos(eps)])
    drift = 19.2 * np.cross(pole, gs)
    assert abs(drift @ ns - 7.3) < 0.15 and abs(drift @ we + 16.2) < 0.15, (drift @ ns, drift @ we)


def test_pulsar_rates():
    """Barker-O'Connell with the published masses reproduces the published GR rates."""
    import pulsars
    pred = pulsars.predictions()
    for m in L.PULSAR_PRECESSION:
        assert abs(pred[m["system"]] / m["gr"] - 1) < 0.01, (m["system"], pred[m["system"]], m["gr"])


def test_solid_tide_nodal_phase():
    """At the 2025 major lunar standstill the Moon reaches high declinations and the
    tidal bulge is less equatorial: the solid tide's 18.6-yr term raises C20."""
    import arcmodel as AM
    import tides
    mjd = np.arange(55_000.0, 55_000.0 + 19 * 365.25, 1.0)
    beta = tides.doodson_arguments(mjd, mjd)
    pb = AM.body_positions(mjd)
    sc, _ = tides.solid_tide_coefficients(pb[:, 0], pb[:, 1], beta)
    N = beta[:, 4]
    A = np.stack([np.ones_like(N), np.cos(N), np.sin(N)], 1)
    c, *_ = np.linalg.lstsq(A, sc[:, 2, 0], rcond=None)
    assert 3.0e-10 < c[1] < 4.5e-10, c


def test_jump_gain(weeks: int = 3, fit_kw: dict | None = None, start_mjd: float = 60_400.0,
                   node_deg: float = 40.0):
    """Synthetic data with frame dragging, fitted week by week without it: every node jump
    must equal GR's Lense-Thirring rate x 7 days (the method passes the signal at full gain).

    With the radiation-pressure scale C_R held fixed (the adopted setting) the gain is
    0.997-0.999 at four seasons; with C_R fitted every week it is 0.81-1.00 (mean ~0.92):
    in eclipse seasons radiation pressure can turn the orbital plane, so a free C_R
    absorbs part of frame dragging. Pass fit_kw={} to see that.
    """
    if fit_kw is None:
        fit_kw = {"fix_cr": True, "cr": 1.13}
    import arcmodel as AM
    import framedrag as FD
    import frames as F
    import nodes
    import odfit
    mjd0 = start_mjd + 69.184 / 86400.0            # an arc starting at 00:00 UTC, in TT
    span = weeks * 7.0
    fmjd = mjd0 + span / 2                         # one gravity-field epoch for truth and fits
    truth_model = AM.build(mjd0, span + 1.0, "lageos1", mu_lt=1.0, field_mjd=fmjd)
    y0 = _circular_state(12_270e3, np.radians(109.84), np.radians(node_deg))
    T, Y = truth_model.run(y0, [1.13, 0, 0, 0, 0], 0.0, span * 86400.0, stride=2)   # every 120 s
    mjd_tt = mjd0 + T / 86400.0
    mjd_utc = F.tt_mjd_to_utc_mjd(mjd_tt)
    r_itrs, v_itrs = F.gcrs_to_itrs(mjd_utc, Y[:, :3], Y[:, 3:6])
    jumps = []
    prev = None
    for w in range(weeks):
        sel = (T >= w * 7 * 86400.0 - 1e-6) & (T <= (w + 1) * 7 * 86400.0 + 1e-6)
        orbit = F.Orbit(sat="lageos1", ac="synthetic", week=f"w{w}", path=None,
                        mjd_utc=mjd_utc[sel], r_itrs=r_itrs[sel], v_itrs=v_itrs[sel], header="")
        fit, model = odfit.fit_arc(orbit, extra_days=1.2,
                                   model_kw={"mu_lt": 0.0, "t_before": 50_000.0, "field_mjd": fmjd},
                                   **fit_kw)
        Tf, Yf = model.run(fit.y0, fit.prm, 0.0, 7.5 * 86400.0, stride=2)
        Tb, Yb = model.run(fit.y0, fit.prm, 0.0, -0.5 * 86400.0, stride=2)
        # one-day windows centred on the week boundaries, as in nodes.process_arc
        keep_f = Tf >= 6.5 * 86400.0 - 1e-6
        early = Tf <= 0.5 * 86400.0 + 1e-6
        Tb = np.concatenate([Tb[::-1][:-1], Tf[early]])
        Yb = np.concatenate([Yb[::-1][:-1], Yf[early]])
        cur = {"mjd0": fit.mjd0, "fwd": (Tf[keep_f], Yf[keep_f]), "bwd": (Tb, Yb)}
        if prev is not None:
            ta = prev["mjd0"] * 86400.0 + prev["fwd"][0]
            tb = cur["mjd0"] * 86400.0 + cur["bwd"][0]
            common, ia, ib = np.intersect1d(np.round(ta, 2), np.round(tb, 2), return_indices=True)
            ea = F.kepler_elements(prev["fwd"][1][ia, :3], prev["fwd"][1][ia, 3:6])
            eb = F.kepler_elements(cur["bwd"][1][ib, :3], cur["bwd"][1][ib, 3:6])
            jumps.append(np.mean((eb["node"] - ea["node"] + np.pi) % (2 * np.pi) - np.pi) / MAS)
        prev = cur
    expect = FD.lt_rate_mas_yr(12_270e3, 0.0) * 7 / 365.25
    gain = np.mean(jumps) / expect
    assert abs(gain - 1) < 0.02, (jumps, expect)
    return gain, jumps


if __name__ == "__main__":
    import inspect
    import sys
    failed = 0
    for name, fn in sorted(inspect.getmembers(sys.modules[__name__], inspect.isfunction)):
        if name.startswith("test_"):
            try:
                fn()
                print(f"ok    {name}")
            except AssertionError as exc:
                failed += 1
                print(f"FAIL  {name}: {exc}")
    sys.exit(1 if failed else 0)
