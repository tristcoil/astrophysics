"""A simple lens model, to see why delays alone cannot weigh the lens.

The arrival time of light from a source at angular position beta, seen at
image position theta, is (Schneider 1985; Blandford & Narayan 1986)

    t(theta) = (D_dt / c) * tau(theta),   tau = |theta - beta|^2 / 2 - psi(theta),

the "Fermat potential": a geometric term (the bent path is longer) minus the
lensing potential psi (time dilation in the lens's gravity well). Images form
where tau is stationary (minima, saddles, maxima), and the delay between two
images is D_dt/c times their Fermat-potential difference. Measure the delay,
know tau, and you have D_dt, hence H0.

The mass model here is an elliptical power law (EPL) plus external shear.
The EPL has convergence (surface density in units of the critical density)
kappa = (2 - t)/2 * (b / R)^t with elliptical radius R = sqrt(q^2 x^2 + y^2)
and t = gamma' - 1, where gamma' is the 3-D logarithmic density slope
(gamma' = 2 is isothermal, i.e. a flat rotation curve). Its deflection is
the series of Tessore & Metcalf (2015, A&A 580, A79).

``mst`` applies the mass-sheet transformation (Falco, Gorenstein & Shapiro
1985): kappa -> lambda kappa + (1 - lambda). The images stay exactly where
they were (the source moves to lambda beta and changes size, which nobody can
see), but every Fermat-potential difference, and so every predicted delay,
is multiplied by lambda.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

import config as C


def _omega(phi, t, q, nmax=400, tol=1e-15):
    """Tessore & Metcalf (2015) angular series Omega(phi)."""
    f = (1 - q) / (1 + q)
    if f < 1e-12:
        return np.exp(1j * phi)
    n_it = min(nmax, int(np.log(tol) / np.log(f)) + 2)
    term = np.exp(1j * phi)
    fact = -f * np.exp(2j * phi)
    out = np.zeros_like(term)
    for n in range(1, n_it):
        out += term
        term = term * (2 * n - (2 - t)) / (2 * n + (2 - t)) * fact
    return out + term


@dataclass
class EPLShear:
    """EPL (Einstein-radius scale b, axis ratio q, angle phi_q, slope gamma')
    centred at (x0, y0), plus external shear (g1, g2), plus an optional
    uniform mass sheet kappa_s (used for the mass-sheet transformation)."""

    b: float
    q: float
    phi_q: float
    gamma: float  # 3-D density slope gamma'
    x0: float
    y0: float
    g1: float = 0.0
    g2: float = 0.0
    lam: float = 1.0  # mass-sheet transformation parameter

    def _local(self, x, y):
        c, s = np.cos(self.phi_q), np.sin(self.phi_q)
        dx, dy = x - self.x0, y - self.y0
        return c * dx + s * dy, -s * dx + c * dy, c, s

    def _epl(self, x, y):
        """EPL deflection (alpha_x, alpha_y) and potential in the sky frame."""
        t = self.gamma - 1.0
        xl, yl, c, s = self._local(x, y)
        z = self.q * xl + 1j * yl
        R = np.abs(z)
        R = np.where(R == 0, 1e-12, R)
        phi = np.angle(z)
        a = 2 * self.b / (1 + self.q) * (self.b / R) ** t * R / self.b * _omega(phi, t, self.q)
        axl, ayl = a.real, a.imag
        psi = (xl * axl + yl * ayl) / (2 - t)
        return c * axl - s * ayl, s * axl + c * ayl, psi

    def deflection(self, x, y):
        ax, ay, _ = self._epl(x, y)
        ax = ax + self.g1 * x + self.g2 * y
        ay = ay + self.g2 * x - self.g1 * y
        # MST: alpha -> lam * alpha + (1 - lam) * theta
        return self.lam * ax + (1 - self.lam) * x, self.lam * ay + (1 - self.lam) * y

    def potential(self, x, y):
        _, _, psi = self._epl(x, y)
        psi = psi + 0.5 * self.g1 * (x * x - y * y) + self.g2 * x * y
        return self.lam * psi + (1 - self.lam) * 0.5 * (x * x + y * y)

    def convergence(self, x, y):
        t = self.gamma - 1.0
        xl, yl, _, _ = self._local(x, y)
        R = np.sqrt((self.q * xl) ** 2 + yl**2)
        return self.lam * (2 - t) / 2 * (self.b / np.maximum(R, 1e-9)) ** t + (1 - self.lam)

    def source(self, x, y):
        ax, ay = self.deflection(x, y)
        return x - ax, y - ay

    def fermat(self, x, y, bx, by):
        """Fermat potential tau (arcsec^2) for source position (bx, by)."""
        return 0.5 * ((x - bx) ** 2 + (y - by) ** 2) - self.potential(x, y)

    def jacobian(self, x, y, h=1e-5):
        """d(beta)/d(theta) by central differences."""
        bxp, byp = self.source(x + h, y)
        bxm, bym = self.source(x - h, y)
        bxq, byq = self.source(x, y + h)
        bxr, byr = self.source(x, y - h)
        return np.array([[(bxp - bxm) / (2 * h), (bxq - bxr) / (2 * h)],
                         [(byp - bym) / (2 * h), (byq - byr) / (2 * h)]])

    def mst(self, lam: float) -> "EPLShear":
        return EPLShear(self.b, self.q, self.phi_q, self.gamma, self.x0, self.y0, self.g1, self.g2, lam * self.lam)


# ---------------------------------------------------------------------------
# HE 0435-1223: the observed configuration
# ---------------------------------------------------------------------------


def he0435_positions():
    """Image and lens positions (arcsec) in a frame with x = west, y = north, origin at the lens.

    Courbin et al. (2011) give (Delta alpha, Delta delta) relative to image A,
    with Delta alpha positive to the east; x = -Delta alpha makes the frame
    right-handed as seen on the sky (it is only a mirror flip either way and
    changes no delay).
    """
    ast = C.LENSES["HE0435"]["astrometry"]
    gx, gy = -ast["G"][0], ast["G"][1]
    img = {k: (-ast[k][0] - gx, ast[k][1] - gy) for k in "ABCD"}
    return img


def fit_positions(images: dict, gamma: float, sigma: float = 0.003, centre_free: bool = False,
                  start: EPLShear | None = None, delay_ratios: dict | None = None,
                  ratio_ref: tuple = ("A", "D")) -> tuple[EPLShear, float]:
    """Fit an EPL + shear of fixed slope to point-image positions.

    Image-plane chi^2 via the local Jacobian: an image-position error dtheta
    maps to a source-plane offset A dtheta, so the image-plane residual of
    each back-projected image relative to the best common source is
    A^-1 (beta_i - beta). The best source position is solved for linearly.
    The lens centre is fixed at the observed light centroid (0, 0) unless
    ``centre_free``.

    ``delay_ratios``: optional {X: (value, error)} for (t_X - t_A) / (t_D - t_A)
    (with ``ratio_ref`` = ("A", "D")). Delay *ratios* do not depend on H0 or
    on a mass sheet, so they are fair lens-model constraints.
    """
    names = list(images)
    xy = np.array([images[k] for k in names])

    def model(p):
        b, e1, e2, g1, g2 = p[:5]
        x0, y0 = (p[5], p[6]) if centre_free else (0.0, 0.0)
        e = np.hypot(e1, e2)
        q = (1 - e) / (1 + e)
        phi = 0.5 * np.arctan2(e2, e1)
        return EPLShear(b, max(q, 0.05), phi, gamma, x0, y0, g1, g2)

    def resid(p):
        m = model(p)
        bx, by = m.source(xy[:, 0], xy[:, 1])
        Ainv = [np.linalg.inv(m.jacobian(x, y)) for x, y in xy]
        # best source: minimize sum |Ainv_i (beta_i - beta)|^2 (linear in beta)
        M = sum(a.T @ a for a in Ainv)
        v = sum(a.T @ a @ np.array([u, w]) for a, u, w in zip(Ainv, bx, by))
        beta = np.linalg.solve(M, v)
        r = [a @ (np.array([u, w]) - beta) for a, u, w in zip(Ainv, bx, by)]
        out = np.concatenate(r) / sigma
        if delay_ratios:
            tau = {k: float(m.fermat(x, y, beta[0], beta[1])) for k, (x, y) in zip(names, xy)}
            r0, r1 = ratio_ref
            den = tau[r1] - tau[r0]
            extra = [((tau[k] - tau[r0]) / den - v) / e for k, (v, e) in delay_ratios.items()]
            out = np.concatenate([out, extra])
        return out

    if start is not None:
        e = (1 - start.q) / (1 + start.q)
        p0 = [start.b, e * np.cos(2 * start.phi_q), e * np.sin(2 * start.phi_q), start.g1, start.g2]
    else:
        p0 = [1.2, 0.05, 0.0, 0.0, 0.0]
    if centre_free:
        p0 += [0.0, 0.0]
    best = None
    rng = np.random.default_rng(3)
    for trial in range(12):
        x0 = np.array(p0) + (0 if trial == 0 else rng.normal(0, [0.05, 0.08, 0.08, 0.04, 0.04] + ([0.02, 0.02] if centre_free else [])))
        try:
            res = least_squares(resid, x0, method="lm", x_scale="jac", max_nfev=4000)
        except (np.linalg.LinAlgError, ValueError):
            continue
        if best is None or res.cost < best.cost:
            best = res
    m = model(best.x)
    return m, float(2 * best.cost)


def source_of(m: EPLShear, images: dict) -> tuple[float, float]:
    xy = np.array(list(images.values()))
    bx, by = m.source(xy[:, 0], xy[:, 1])
    return float(np.mean(bx)), float(np.mean(by))


def fermat_differences(m: EPLShear, images: dict, ref: str = "A") -> dict:
    """tau(ref) - tau(X) for each image X: positive when X arrives after ref... (arcsec^2).

    Arrival time is D_dt/c * tau, so t_X - t_ref = D_dt/c * (tau_X - tau_ref).
    """
    bx, by = source_of(m, images)
    tau = {k: float(m.fermat(x, y, bx, by)) for k, (x, y) in images.items()}
    return {k: tau[k] - tau[ref] for k in images if k != ref}


def ddt_from(delay_days: float, dphi_arcsec2: float) -> float:
    """Time-delay distance (Mpc) from a delay (days) and a Fermat-potential difference (arcsec^2)."""
    c_mpc_day = C.C_KM_S * C.DAY_S / C.MPC_KM
    return c_mpc_day * delay_days / (dphi_arcsec2 * C.ARCSEC_RAD**2)
