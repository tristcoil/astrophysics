"""Orbital elements, precession rates, and the formulas they are checked against.

Three ways to say "the perihelion turned"
=========================================

An orbit's *eccentricity vector* points from the Sun to the perihelion. Its
direction is what precesses. There are three common ways to put a number on that,
and for Mercury they differ by up to 14 arcseconds per century:

1. **The in-plane angle** ``psi``: how far the eccentricity vector has turned
   about the orbit's own normal, relative to the distant stars. It does not
   depend on any coordinate system, and it is what the textbook numbers mean:
   general relativity's 43"/cy, the planets' ~532"/cy and MESSENGER's
   575.31"/cy total (Park et al. 2017). We compute it from successive samples as
   ``d psi = (e_prev x e_next) . h`` (unit vectors), which equals the familiar
   ``d omega + cos(i) d Omega`` in any fixed frame.
2. **The longitude of perihelion** ``varpi = Omega + omega`` in the fixed ecliptic
   of J2000 -- a "dogleg" angle, measured along the ecliptic to the node and then
   along the orbit. While Mercury's orbit plane slides along the ecliptic
   (``Omega`` falls by ~452"/cy), the dogleg changes by ``(1 - cos i) dOmega``,
   so ``varpi`` advances ~3.4"/cy *less* than ``psi``.
3. The same ``varpi`` on the ICRF equator, where Mercury's orbit is tilted by
   28.6 degrees: ~14"/cy less again.

Conventions: positions in au, velocities in au/day, GM in au^3/day^2, angles in
radians internally and arcseconds per Julian century (``"/cy``) when reported.
"""

from __future__ import annotations

import numpy as np

from config import ARCSEC_PER_RAD, CENTURY_DAYS

# --------------------------------------------------------------------------
# Elements
# --------------------------------------------------------------------------


def elements(r: np.ndarray, v: np.ndarray, mu: float) -> dict:
    """Osculating elements of relative states ``r``, ``v`` with shape (3, N).

    Returns a, e, inc, Omega, omega, varpi (radians; Omega/omega/varpi unwrapped
    along the series), plus the unit eccentricity and angular-momentum vectors.
    """
    r, v = np.asarray(r, float), np.asarray(v, float)
    if r.ndim == 1:
        r, v = r[:, None], v[:, None]
    rn = np.linalg.norm(r, axis=0)
    h = np.cross(r.T, v.T).T
    hn = np.linalg.norm(h, axis=0)
    evec = np.cross(v.T, h.T).T / mu - r / rn
    e = np.linalg.norm(evec, axis=0)
    a = 1.0 / (2.0 / rn - (v * v).sum(0) / mu)
    inc = np.arccos(np.clip(h[2] / hn, -1.0, 1.0))
    Omega = np.arctan2(h[0], -h[1])
    node = np.array([np.cos(Omega), np.sin(Omega), np.zeros_like(Omega)])
    hhat = h / hn
    perp = np.cross(hhat.T, node.T).T  # in the plane, 90 deg ahead of the node
    omega = np.arctan2((evec * perp).sum(0), (evec * node).sum(0))
    Omega_u = np.unwrap(Omega)
    omega_u = np.unwrap(omega)
    with np.errstate(invalid="ignore", divide="ignore"):  # e = 0: no perihelion direction
        ehat = evec / e
    return {"a": a, "e": e, "inc": inc, "Omega": Omega_u, "omega": omega_u,
            "varpi": Omega_u + omega_u, "ehat": ehat, "hhat": hhat}


def in_plane_angle(ehat: np.ndarray, hhat: np.ndarray) -> np.ndarray:
    """Cumulative turn of the eccentricity vector about the orbit normal (rad).

    Frame-independent: each step's rotation is ``atan2((e1 x e2).h, e1.e2)``
    with ``h`` the mean of the two normals. Starts at zero.
    """
    e1, e2 = ehat[:, :-1], ehat[:, 1:]
    hm = hhat[:, :-1] + hhat[:, 1:]
    hm = hm / np.linalg.norm(hm, axis=0)
    cross = np.cross(e1.T, e2.T).T
    step = np.arctan2((cross * hm).sum(0), (e1 * e2).sum(0))
    return np.concatenate([[0.0], np.cumsum(step)])


def angle_series(r: np.ndarray, v: np.ndarray, mu: float, rot: np.ndarray | None = None) -> dict:
    """Elements plus the in-plane angle ``psi`` for a time series of states.

    ``rot`` rotates the states into the frame in which Omega/omega/varpi are
    wanted (``psi`` does not care).
    """
    if rot is not None:
        r, v = rot @ r, rot @ v
    el = elements(r, v, mu)
    el["psi"] = in_plane_angle(el["ehat"], el["hhat"])
    return el


# --------------------------------------------------------------------------
# Measuring a secular rate
# --------------------------------------------------------------------------


def fit_secular(t_cy: np.ndarray, angle: np.ndarray, deg: int = 2) -> dict:
    """Least-squares polynomial in time (centuries from J2000) to an angle series.

    The coefficient of t is the rate *at J2000*; with ``deg=2`` the slow change
    of the rate itself (planetary orbits evolve) does not bias it even on an
    asymmetric span. Rates in "/cy, curvature as d(rate)/dt in "/cy^2.
    """
    A = np.vander(t_cy, deg + 1, increasing=True)
    coef, *_ = np.linalg.lstsq(A, angle, rcond=None)
    resid = angle - A @ coef
    out = {"rate": coef[1] * ARCSEC_PER_RAD,
           "rms_resid_arcsec": float(resid.std() * ARCSEC_PER_RAD),
           "max_resid_arcsec": float(np.abs(resid).max() * ARCSEC_PER_RAD),
           "span_cy": [float(t_cy.min()), float(t_cy.max())], "deg": deg, "n": int(len(t_cy))}
    if deg >= 2:
        out["curvature"] = 2.0 * coef[2] * ARCSEC_PER_RAD
    return out


def rate_stability(t_cy: np.ndarray, angle: np.ndarray, halves=(1.0, 2.0, 3.0, 4.0, 4.495)) -> dict:
    """The J2000 rate from symmetric spans of several lengths (linear fits).

    A symmetric span makes a linear fit unbiased by the slow curvature; the
    spread over span lengths measures how much the periodic terms that the fit
    does not model can move the answer.
    """
    rates = {}
    for h in halves:
        m = np.abs(t_cy) <= h + 1e-12
        if t_cy[m].min() > -h + 0.01 or t_cy[m].max() < h - 0.01:
            continue
        rates[f"+-{h * 100:.0f}yr"] = fit_secular(t_cy[m], angle[m], deg=1)["rate"]
    vals = np.array(list(rates.values()))
    return {"rates": rates, "spread": float(vals.max() - vals.min()) if len(vals) else None}


# --------------------------------------------------------------------------
# Analytic predictions
# --------------------------------------------------------------------------


def mean_motion(mu: float, a: float) -> float:
    """rad/day."""
    return np.sqrt(mu / a**3)


def per_century(rate_rad_per_day: float) -> float:
    return rate_rad_per_day * CENTURY_DAYS * ARCSEC_PER_RAD


def schwarzschild_per_orbit(gm: float, a: float, e: float, c: float, *, gamma=1.0, beta=1.0):
    """Perihelion advance per orbit, radians: 6 pi GM / (c^2 a (1-e^2)) x PPN factor.

    The PPN factor (2 + 2 gamma - beta) / 3 is 1 in general relativity;
    Newton (no relativistic terms at all) corresponds to zero precession.
    """
    return 6.0 * np.pi * gm / (c**2 * a * (1.0 - e**2)) * (2.0 + 2.0 * gamma - beta) / 3.0


def schwarzschild_rate(gm: float, a: float, e: float, c: float, **ppn) -> float:
    """The same in "/cy, using the Keplerian period of the orbit."""
    n = mean_motion(gm, a)
    return per_century(schwarzschild_per_orbit(gm, a, e, c, **ppn) * n / (2.0 * np.pi))


def compactness(gm: float, a: float, e: float, c: float) -> float:
    """GM / (c^2 a (1 - e^2)): how deep in gravity an orbit is (dimensionless)."""
    return gm / (c**2 * a * (1.0 - e**2))


def j2_in_plane_rate(gm: float, a: float, e: float, j2: float, radius: float, inc_eq: float) -> float:
    """Orbit-averaged in-plane apsidal rate from the central body's J2, "/cy.

    ``inc_eq`` is the orbit's tilt to the body's equator. From the classical
    rates dOmega/dt = -3/2 n J2 (R/p)^2 cos i and domega/dt = 3/4 n J2 (R/p)^2
    (5 cos^2 i - 1): domega + cos i dOmega = 3/4 n J2 (R/p)^2 (3 cos^2 i - 1).
    """
    n = mean_motion(gm, a)
    p = a * (1.0 - e**2)
    return per_century(0.75 * n * j2 * (radius / p) ** 2 * (3.0 * np.cos(inc_eq) ** 2 - 1.0))


def lense_thirring_in_plane_rate(g_spin: float, a: float, e: float, c: float, inc_eq: float) -> float:
    """In-plane apsidal rate from the central body's spin (gravitomagnetism), "/cy.

    ``g_spin`` is G x the angular momentum (au^5/day^3). The node advances by
    2 G S / (c^2 a^3 (1-e^2)^1.5) and the perigee by -6 G S cos i / (...), so the
    in-plane combination is -4 G S cos i / (c^2 a^3 (1-e^2)^1.5).
    """
    return per_century(-4.0 * g_spin * np.cos(inc_eq) / (c**2 * a**3 * (1.0 - e**2) ** 1.5))


def laplace_coefficient(s: float, j: int, alpha: float) -> float:
    """b_s^(j)(alpha) = (1/pi) int_0^{2pi} cos(j psi) / (1 - 2 alpha cos psi + alpha^2)^s dpsi."""
    # The integrand is smooth and periodic: the trapezoid rule on an even grid
    # converges exponentially (here to ~1e-15 with 4096 points for alpha < 0.8).
    psi = np.linspace(0.0, 2.0 * np.pi, 4096, endpoint=False)
    f = np.cos(j * psi) / (1.0 - 2.0 * alpha * np.cos(psi) + alpha**2) ** s
    return 2.0 * f.mean()


def laplace_lagrange_rate(gm_central: float, gm_inner: float, a_inner: float, e_inner: float,
                          varpi_inner: float, gm_outer: float, a_outer: float, e_outer: float,
                          varpi_outer: float) -> dict:
    """First-order secular theory (Laplace-Lagrange) for an inner orbit's perihelion.

    The classical calculation Le Verrier and Newcomb did by hand: average the
    outer planet's pull over both orbits, keep terms to first order in the
    masses and second order in the eccentricities (Murray & Dermott 1999, sect.
    7.2). The inner perihelion turns at

        dvarpi/dt = A11 + A12 (e_outer / e_inner) cos(varpi_inner - varpi_outer)

    with A11 = n/4 (m'/M) alpha^2 b_{3/2}^{(1)}(alpha) and
    A12 = -n/4 (m'/M) alpha^2 b_{3/2}^{(2)}(alpha), alpha = a/a'. Returns "/cy.
    """
    alpha = a_inner / a_outer
    n = mean_motion(gm_central + gm_inner, a_inner)
    mass = gm_outer / (gm_central + gm_inner)
    b1 = laplace_coefficient(1.5, 1, alpha)
    b2 = laplace_coefficient(1.5, 2, alpha)
    a11 = 0.25 * n * mass * alpha**2 * b1
    a12 = -0.25 * n * mass * alpha**2 * b2
    coupling = a12 * (e_outer / e_inner) * np.cos(varpi_inner - varpi_outer)
    return {"free": per_century(a11), "coupling": per_century(coupling),
            "total": per_century(a11 + coupling), "alpha": alpha}


def ppn_precession_factor(gamma: float, beta: float) -> float:
    """(2 + 2 gamma - beta) / 3: the perihelion advance relative to general relativity."""
    return (2.0 + 2.0 * gamma - beta) / 3.0


# --------------------------------------------------------------------------
# Gauss's method: the exact first-order secular rate
# --------------------------------------------------------------------------


def kepler_points(a: float, e: float, inc: float, Omega: float, omega: float, mu: float,
                  n: int) -> tuple[np.ndarray, np.ndarray]:
    """Positions and velocities (3, n) at n equally spaced mean anomalies."""
    M = 2.0 * np.pi * (np.arange(n) + 0.5) / n
    E = M.copy()
    for _ in range(50):
        E = E - (E - e * np.sin(E) - M) / (1.0 - e * np.cos(E))
    cosE, sinE = np.cos(E), np.sin(E)
    b = a * np.sqrt(1.0 - e**2)
    x, y = a * (cosE - e), b * sinE
    nmot = np.sqrt(mu / a**3)
    edot = nmot / (1.0 - e * cosE)
    vx, vy = -a * sinE * edot, b * cosE * edot
    cO, sO, ci, si, cw, sw = (np.cos(Omega), np.sin(Omega), np.cos(inc), np.sin(inc),
                              np.cos(omega), np.sin(omega))
    P = np.array([cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si])
    Q = np.array([-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si])
    r = P[:, None] * x + Q[:, None] * y
    v = P[:, None] * vx + Q[:, None] * vy
    return r, v


def elements_angles(r: np.ndarray, v: np.ndarray, mu: float) -> dict:
    """Scalar osculating elements (a, e, inc, Omega, omega) of one state."""
    el = elements(r, v, mu)
    return {k: float(np.ravel(el[k])[0]) for k in ("a", "e", "inc", "Omega", "omega")}


def gauss_secular_rate(mu_c: float, r_in: np.ndarray, v_in: np.ndarray, gm_out: float,
                       r_out: np.ndarray, v_out: np.ndarray, n: int = 720) -> dict:
    """First-order secular perihelion rate of an inner orbit, exact in e and i ("/cy).

    Gauss (1818): spread the outer planet along its ellipse with a density
    proportional to the time it spends there (a "Gauss ring"), compute the
    ring's pull at every point of the inner orbit, and average the resulting
    change of the inner eccentricity vector over the inner orbit, weighting by
    time again. The planet's indirect pull on the Sun averages to zero over its
    orbit. ``mu_c`` is G(M_sun + m_inner); states are heliocentric, any frame.
    The in-plane rate is <de/dt> . (h x e) / e^2.
    """
    ei = elements_angles(r_in, v_in, mu_c)
    eo = elements_angles(r_out, v_out, mu_c)  # mu only shapes the outer ring's timing
    ri, vi = kepler_points(ei["a"], ei["e"], ei["inc"], ei["Omega"], ei["omega"], mu_c, n)
    ro, _ = kepler_points(eo["a"], eo["e"], eo["inc"], eo["Omega"], eo["omega"], mu_c, n)
    F = np.zeros_like(ri)
    for k in range(0, n, 120):  # chunks keep memory small
        d = ro[:, None, :] - ri[:, k:k + 120, None]
        F[:, k:k + 120] = gm_out * (d / np.linalg.norm(d, axis=0) ** 3).mean(axis=2)
    h = np.cross(ri.T, vi.T).T
    dedt = (np.cross(F.T, h.T).T + np.cross(vi.T, np.cross(ri.T, F.T)).T) / mu_c
    evec = np.cross(vi.T, h.T).T / mu_c - ri / np.linalg.norm(ri, axis=0)
    e0 = evec.mean(axis=1)  # constant along a Kepler orbit
    h0 = h.mean(axis=1)
    q = np.cross(h0, e0) / np.linalg.norm(np.cross(h0, e0))
    rate = dedt.mean(axis=1) @ q / np.linalg.norm(e0)
    return {"rate": per_century(rate), "alpha": ei["a"] / eo["a"]}
