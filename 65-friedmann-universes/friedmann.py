"""Solving the Friedmann equations: how a universe evolves, given what it contains.

A universe that is the same everywhere and in every direction is described by a
single function of time, the **scale factor** a(t): the size of the universe
relative to today (a = 1 now; a = 1/2 when everything was twice as close). Its
growth rate H = (da/dt)/a is the Hubble parameter, H0 today. Einstein's
equations reduce to the **Friedmann equation**

    H(a)^2 = H0^2 E(a)^2,   E(a)^2 = Or a^-4 + Om a^-3 + Ok a^-2 + OL a^(-3(1+w)),

in which each ingredient dilutes in its own way as space grows: radiation (Or)
as a^-4 (fewer photons per volume *and* each one stretched), matter (Om) as
a^-3, the curvature term (Ok = 1 - Or - Om - OL) as a^-2, and dark energy with
equation of state w = pressure / (density c^2) as a^(-3(1+w)) -- not at all for
a cosmological constant (w = -1). The densities are today's, as fractions of the
critical density 3 H0^2 / (8 pi G), so they sum to one.

Multiplying by a^4 gives the polynomial at the heart of this module,

    P(a) = a^4 E(a)^2 = Or + Om a + Ok a^2 + OL a^(1-3w),     (da/dt)^2 = H0^2 P(a) / a^2,

which must stay positive for the universe to exist at that size. Today P(1) = 1.
Everything follows from where P has roots:

* a root at a_min < 1 is a **bounce**: running time backwards the universe stops
  shrinking at a_min and re-expands. It never had a Big Bang.
* a root at a_max > 1 is a **turnaround**: the expansion stops and reverses
  into a **Big Crunch**.
* no roots: a Big Bang (a -> 0) in the past and expansion forever.
* a double root (P = P' = 0) is Einstein's static universe: a universe on that
  boundary approaches it forever (Eddington-Lemaitre), and one just beside it
  *loiters* there for a long time (Lemaitre's "hesitating" universe).

Time follows by integration, dt = a da / (H0 sqrt(P(a))). The square-root
singularities at turning points are removed by substitution (a = r + s^2), so
the integrals are smooth and quadrature is accurate to ~1e-12. The second-order
**acceleration equation**, d^2a/dt^2 = H0^2 (a P'(a) - 2 P(a)) / (2 a^3), is
integrated as an independent check (``integrate_ode``): it passes through
turning points without any special treatment.

Units: times are in units of the Hubble time 1/H0 unless a function says Gyr.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy import integrate, optimize

import config as cfg

# --- Radiation today ----------------------------------------------------------

_G = 6.67430e-11
_SIGMA_SB = 5.670374419e-8
_C = 2.99792458e8
_MPC_M = 3.085677581e22
#: Photon density omega_gamma = Omega_gamma h^2 of a 2.7255 K blackbody.
OMEGA_GAMMA = (4.0 * _SIGMA_SB * cfg.T_CMB**4 / _C**3) * 8.0 * np.pi * _G / (
    3.0 * (100e3 / _MPC_M) ** 2)
#: Photons plus N_eff neutrino species, all treated as relativistic (the toy
#: universes; the data fits in ``background.py`` give one neutrino its mass).
OMEGA_R_TOY = OMEGA_GAMMA * (1.0 + cfg.N_EFF * 7.0 / 8.0 * (4.0 / 11.0) ** (4.0 / 3.0))


def omega_r_today(H0: float) -> float:
    """Omega_r for a universe expanding at H0 today: ~9.2e-5 at 67.4 km/s/Mpc."""
    return OMEGA_R_TOY / (H0 / 100.0) ** 2


OR_GALLERY = omega_r_today(cfg.H0_GALLERY)


def hubble_time_gyr(H0: float) -> float:
    """1/H0 in Gyr: the age a universe would have had if it always expanded at H0."""
    return cfg.HUBBLE_TIME_GYR_X_H0 / H0


# --- One universe ---------------------------------------------------------------


@dataclass(frozen=True)
class Universe:
    """A Friedmann universe: today's density parameters, dark energy's w, and H0."""

    Om: float
    OL: float
    Or: float = 0.0
    w: float = -1.0
    H0: float = cfg.H0_GALLERY
    name: str = ""

    @property
    def Ok(self) -> float:
        return 1.0 - self.Om - self.OL - self.Or

    @property
    def n(self) -> float:
        """Power of a carried by dark energy in P(a): 1 - 3w (4 for Lambda)."""
        return 1.0 - 3.0 * self.w

    @property
    def t_hubble_gyr(self) -> float:
        return hubble_time_gyr(self.H0)

    def with_(self, **kw) -> "Universe":
        return replace(self, **kw)

    # -- The polynomial and its derivative ------------------------------------

    def P(self, a):
        a = np.asarray(a, float)
        return self.Or + self.Om * a + self.Ok * a**2 + self.OL * a**self.n

    def dP(self, a):
        a = np.asarray(a, float)
        return self.Om + 2.0 * self.Ok * a + self.n * self.OL * a ** (self.n - 1.0)

    def E2(self, a):
        """(H(a)/H0)^2."""
        a = np.asarray(a, float)
        return self.P(a) / a**4

    def hubble(self, z):
        """H(z) in km/s/Mpc (NaN where the universe never reached that redshift)."""
        e2 = self.E2(1.0 / (1.0 + np.asarray(z, float)))
        with np.errstate(invalid="ignore"):
            return self.H0 * np.sqrt(np.where(e2 > 0, e2, np.nan))

    def acceleration(self, a):
        """d^2a/dt^2 in units of H0^2: matter and radiation pull, Lambda pushes."""
        a = np.asarray(a, float)
        return (a * self.dP(a) - 2.0 * self.P(a)) / (2.0 * a**3)

    def q0(self) -> float:
        """Today's deceleration parameter q0 = -a a'' / a'^2 (negative: speeding up)."""
        return self.Or + 0.5 * self.Om + 0.5 * (1.0 + 3.0 * self.w) * self.OL

    # -- Turning points -----------------------------------------------------------

    def roots(self) -> np.ndarray:
        """Every a > 0 where P(a) = 0 (sorted), simple and double roots alike."""
        if self.w == -1.0:
            r = np.roots([self.OL, 0.0, self.Ok, self.Om, self.Or])
            real = r[np.abs(r.imag) <= 1e-9 * np.maximum(1.0, np.abs(r.real))].real
            real = real[real > 0]
            # polish (np.roots loses digits on near-double roots)
            out = []
            for x in np.sort(real):
                try:
                    x = optimize.newton(self.P, x, fprime=self.dP, tol=1e-15, maxiter=50)
                except (RuntimeError, OverflowError):
                    pass
                out.append(float(x))
            return np.unique(np.round(out, 12))
        # general w: scan log a, bracket sign changes, catch tangencies
        x = np.logspace(-10, 10, 40001)
        p = self.P(x)
        out = [optimize.brentq(self.P, x[i], x[i + 1], xtol=1e-15, rtol=1e-14)
               for i in np.flatnonzero(np.sign(p[:-1]) * np.sign(p[1:]) < 0)]
        mins = np.flatnonzero((p[1:-1] < p[:-2]) & (p[1:-1] < p[2:]) & (p[1:-1] > 0)
                              & (p[1:-1] < 1e-10))
        out += [float(x[i + 1]) for i in mins]
        return np.array(sorted(out))

    def bounce_a(self) -> float | None:
        """The bounce: the largest root of P below a = 1, if there is one."""
        r = self.roots()
        r = r[r < 1.0]
        return float(r.max()) if r.size else None

    def turnaround_a(self) -> float | None:
        """The maximum size before recollapse: the smallest root above a = 1."""
        r = self.roots()
        r = r[r > 1.0]
        return float(r.min()) if r.size else None

    def _is_double(self, r: float) -> bool:
        return abs(self.dP(r)) < 1e-7 * max(1.0, abs(self.Om) + abs(self.Ok) + abs(self.OL))

    # -- Time integrals --------------------------------------------------------

    def _t_from_zero(self, a: float) -> float:
        """Time since a = 0 (the Big Bang), in 1/H0, with a = u^2: smooth integrand."""
        if a <= 0:
            return 0.0

        def g(u):
            return 2.0 * u**3 / np.sqrt(self.P(u * u))

        # diverges only if P vanishes at least as fast as a^4 at a -> 0 (no
        # radiation, matter or curvature): an eternal past, as in de Sitter space
        powers = [p for p, c in ((0, self.Or), (1, self.Om), (2, self.Ok), (self.n, self.OL))
                  if c != 0]
        if min(powers) >= 4:
            return np.inf
        val, _ = integrate.quad(g, 0.0, np.sqrt(a), epsabs=0.0, epsrel=1e-12, limit=400)
        return float(val)

    def _P_over_distance(self, r: float, x):
        """P(x) / (x - r) for a root r of P, without the cancellation that ruins
        P(x) itself close to r: a Taylor series there, the direct ratio further out."""
        x = np.asarray(x, float)
        d = x - r
        n = self.n
        p1 = self.dP(r)
        p2 = 2.0 * self.Ok + n * (n - 1.0) * self.OL * r ** (n - 2.0)
        p3 = n * (n - 1.0) * (n - 2.0) * self.OL * r ** (n - 3.0)
        p4 = n * (n - 1.0) * (n - 2.0) * (n - 3.0) * self.OL * r ** (n - 4.0)
        series = p1 + p2 * d / 2.0 + p3 * d**2 / 6.0 + p4 * d**3 / 24.0
        with np.errstate(divide="ignore", invalid="ignore"):
            direct = self.P(x) / d
        return np.where(np.abs(d) < 1e-3 * max(r, 1e-12), series, direct)

    def _t_from_root(self, r: float, a: float) -> float:
        """Time to go from a turning point r to a, in 1/H0: a = r +- s^2 removes the
        1/sqrt singularity (P ~ P'(r) (a - r) there)."""
        if self._is_double(r):
            return np.inf
        sign = 1.0 if a >= r else -1.0
        smax = np.sqrt(abs(a - r))

        def g(s):
            x = r + sign * s * s
            # 2 s x / sqrt(|P(x)|) with |P(x)| = s^2 |P(x) / (x - r)|
            return 2.0 * x / np.sqrt(np.abs(self._P_over_distance(r, x)))

        val, _ = integrate.quad(g, 0.0, smax, epsabs=0.0, epsrel=1e-12, limit=400)
        return float(val)

    def _t_between(self, a1: float, a2: float) -> float:
        """Time between two sizes on a stretch of expansion with no turning point."""
        val, _ = integrate.quad(lambda a: a / np.sqrt(self.P(a)), a1, a2, epsabs=0.0,
                                epsrel=1e-12, limit=400)
        return float(val)

    # -- The story of this universe ------------------------------------------------

    def fate(self) -> dict:
        """Beginning and end: Big Bang or bounce; Big Crunch, Big Rip or forever.

        Returns a dict with times in Gyr *relative to today* (negative = past)
        plus the classification and the key sizes.
        """
        tH = self.t_hubble_gyr
        out = {"name": self.name, "Om": self.Om, "OL": self.OL, "Or": self.Or, "Ok": self.Ok,
               "w": self.w, "H0": self.H0, "q0": self.q0()}
        a_min = self.bounce_a()
        if a_min is not None:
            out["past"] = "no Big Bang: bounce" if not self._is_double(a_min) else \
                "emerged from a static universe (Eddington-Lemaitre)"
            out["a_min"] = a_min
            out["z_max"] = 1.0 / a_min - 1.0
            out["t_bounce_gyr"] = -self._t_from_root(a_min, 1.0) * tH
            out["age_gyr"] = np.inf
        else:
            age = self._t_from_zero(1.0)
            out["past"] = "Big Bang" if np.isfinite(age) else "no beginning (eternal past)"
            out["age_gyr"] = age * tH
            out["t_bang_gyr"] = -age * tH
        a_max = self.turnaround_a()
        if a_max is not None:
            if self._is_double(a_max):
                out["future"] = "approaches a static universe forever"
            else:
                t_turn = self._t_from_root(a_max, 1.0) * tH
                out["future"] = "Big Crunch"
                out["a_max"] = a_max
                out["t_turnaround_gyr"] = t_turn
                # by time symmetry the collapse from a_max back to a -> 0 takes
                # as long as the expansion from a -> 0 to a_max did
                t_up = out["age_gyr"] + t_turn
                out["t_crunch_gyr"] = t_turn + t_up
                out["lifetime_gyr"] = 2.0 * t_up
        else:
            rip = self.time_to_rip_gyr()
            if np.isfinite(rip):
                out["future"] = "Big Rip"
                out["t_rip_gyr"] = rip
            elif self.OL > 0 and self.w < -1.0 / 3.0:
                out["future"] = "expands forever, accelerating"
            else:
                out["future"] = "expands forever, slowing down"
        return out

    def age_gyr(self) -> float:
        """Time since the Big Bang (inf if there was none)."""
        if self.bounce_a() is not None:
            return np.inf
        return self._t_from_zero(1.0) * self.t_hubble_gyr

    def time_to_rip_gyr(self) -> float:
        """Phantom dark energy (w < -1) grows as space expands and tears it apart in
        a finite time: t_rip - t0 = int_1^inf a da / sqrt(P). Infinite otherwise."""
        if self.w >= -1.0 or self.OL <= 0 or self.turnaround_a() is not None:
            return np.inf

        def g(v):  # a = 1/v, da = -dv / v^2
            return 1.0 / (v**3 * np.sqrt(self.P(1.0 / v)))

        val, _ = integrate.quad(g, 0.0, 1.0, epsabs=0.0, epsrel=1e-11, limit=400)
        return float(val) * self.t_hubble_gyr

    def lookback_gyr(self, z) -> np.ndarray:
        """How long ago light now arriving from redshift z set out (Big Bang universes)."""
        z = np.atleast_1d(np.asarray(z, float))
        return np.array([self._t_between(1.0 / (1.0 + zi), 1.0) for zi in z]) * self.t_hubble_gyr

    def time_of_a(self, a) -> np.ndarray:
        """Cosmic time (since the Big Bang) at scale factor a <= 1, in Gyr."""
        a = np.atleast_1d(np.asarray(a, float))
        return np.array([self._t_from_zero(x) for x in a]) * self.t_hubble_gyr

    # -- The whole history, for plotting ---------------------------------------

    def history(self, t_past_gyr: float = -40.0, t_future_gyr: float = 60.0,
                n: int = 600) -> list[tuple[np.ndarray, np.ndarray]]:
        """a(t) over a window of time around today, as a list of (t_gyr, a) pieces.

        Built from the time integrals, branch by branch: the expanding branch
        through today; for a bounce, its mirror image (the contracting phase
        before the bounce); for a recollapse, the mirror image after turnaround.
        The window is in Gyr relative to today.
        """
        tH = self.t_hubble_gyr
        fate = self.fate()
        a_min = fate.get("a_min")
        a_max = fate.get("a_max")
        eternal_past = a_min is None and not np.isfinite(fate.get("age_gyr", np.inf))
        m = max(n // 2, 20)
        # sizes on the expanding branch: past part (start -> 1), future part (1 -> end)
        if a_min is not None:  # bounce: dense near the turning point
            past = a_min + np.linspace(0.0, np.sqrt(1.0 - a_min), m) ** 2
        elif eternal_past:  # de Sitter-like: a -> 0 only as t -> -infinity
            past = np.geomspace(1e-6, 1.0, m)
        else:
            past = np.concatenate([[0.0], np.geomspace(1e-6, 1.0, m)])
        if a_max is not None:  # turnaround: dense near it
            fut = a_max - np.linspace(np.sqrt(a_max - 1.0), 0.0, m) ** 2
        else:
            fut = np.geomspace(1.0, self._size_at_future_time(t_future_gyr / tH), m)
        past[-1], fut[0] = 1.0, 1.0

        def seg(a1, a2):  # time to grow from a1 to a2, each piece singularity-free
            if a1 == 0.0:
                return self._t_from_zero(a2)
            if a_min is not None and a1 == a_min:
                return self._t_from_root(a_min, a2)
            if a_max is not None and a2 == a_max:
                return self._t_from_root(a_max, a1)
            return self._t_between(a1, a2)

        dt_past = np.array([seg(past[i], past[i + 1]) for i in range(past.size - 1)])
        dt_fut = np.array([seg(fut[i], fut[i + 1]) for i in range(fut.size - 1)])
        t_past = -np.concatenate([np.cumsum(dt_past[::-1])[::-1], [0.0]])
        t_fut = np.concatenate([[0.0], np.cumsum(dt_fut)])
        x = np.concatenate([past, fut[1:]])
        t_gyr = np.concatenate([t_past, t_fut[1:]]) * tH
        pieces = [(t_gyr, x)]
        if a_min is not None:  # contracting before the bounce: mirror about t_bounce
            tb = t_gyr[0]
            pieces.insert(0, ((2 * tb - t_gyr)[::-1], x[::-1]))
        if a_max is not None:  # collapsing after turnaround: mirror about t_turn
            tt = t_gyr[-1]
            pieces.append((2 * tt - t_gyr[::-1], x[::-1]))
        clipped = []
        for tg, aa in pieces:
            keep = (tg >= t_past_gyr - 1e-9) & (tg <= t_future_gyr + 1e-9)
            if keep.sum() > 1:
                clipped.append((tg[keep], aa[keep]))
        return clipped

    def _size_at_future_time(self, t_rel: float) -> float:
        """Scale factor reached t_rel (in 1/H0) after today, expanding forever."""
        if t_rel <= 0:
            return 1.0
        rip = self.time_to_rip_gyr() / self.t_hubble_gyr
        if np.isfinite(rip) and t_rel >= rip:
            return 1e6

        def f(lna):
            return self._t_between(1.0, np.exp(lna)) - t_rel

        hi = 1.0
        while f(hi) < 0:
            hi *= 2.0
            if hi > 200:
                return float(np.exp(hi))
        return float(np.exp(optimize.brentq(f, 0.0, hi, xtol=1e-10)))


# --- The acceleration equation, integrated directly (an independent check) -------


def integrate_ode(u: Universe, t_span_hubble: tuple[float, float], *, a_floor: float = 1e-4,
                  rtol: float = 1e-11, atol: float = 1e-13):
    """Integrate a'' = (a P' - 2P) / (2 a^3) from today (a = 1, a' = 1) in units
    of 1/H0. Stops if a drops below ``a_floor`` (a Big Bang or Big Crunch).

    Returns the scipy solution with dense output, for t in ``t_span_hubble``.
    """
    def rhs(_t, y):
        return [y[1], float(u.acceleration(y[0]))]

    def hit_floor(_t, y):
        return y[0] - a_floor

    hit_floor.terminal = True
    hit_floor.direction = -1
    return integrate.solve_ivp(rhs, t_span_hubble, [1.0, 1.0], method="DOP853", rtol=rtol,
                               atol=atol, dense_output=True, events=hit_floor)


# --- Closed-form solutions, for the textbook universes and for testing ----------


def age_einstein_de_sitter() -> float:
    """Flat, matter only: a = (t/t0)^(2/3), t0 = 2/(3 H0)."""
    return 2.0 / 3.0


def age_flat_lcdm(Om: float) -> float:
    """Flat matter + Lambda, no radiation: t0 H0 = 2/(3 sqrt(OL)) asinh(sqrt(OL/Om))."""
    OL = 1.0 - Om
    return 2.0 / (3.0 * np.sqrt(OL)) * np.arcsinh(np.sqrt(OL / Om))


def age_open_matter(Om: float) -> float:
    """Open, matter only (Om < 1): the hyperbolic 'cycloid'."""
    k = 1.0 - Om
    return 1.0 / k - Om / (2.0 * k**1.5) * np.arccosh(2.0 / Om - 1.0)


def closed_matter_parametric(Om: float, eta):
    """Closed, matter only (Om > 1): a cycloid, from Big Bang (eta = 0) to Big
    Crunch (eta = 2 pi). Returns (t H0, a) at development angle eta."""
    k = Om - 1.0
    eta = np.asarray(eta, float)
    return Om / (2.0 * k**1.5) * (eta - np.sin(eta)), Om / (2.0 * k) * (1.0 - np.cos(eta))


def age_closed_matter(Om: float) -> float:
    eta0 = np.arccos(1.0 - 2.0 * (Om - 1.0) / Om)
    return float(closed_matter_parametric(Om, eta0)[0])


def lifetime_closed_matter(Om: float) -> float:
    """Big Bang to Big Crunch: pi Om / (Om - 1)^(3/2), in 1/H0."""
    return np.pi * Om / (Om - 1.0) ** 1.5


def t_radiation_matter(a, Or: float, Om: float):
    """Flat radiation + matter (Or + Om = 1): the exact solution
    t H0 = 2 / (3 Om^2) [(Om a - 2 Or) sqrt(Or + Om a) + 2 Or^(3/2)]."""
    a = np.asarray(a, float)
    return 2.0 / (3.0 * Om**2) * ((Om * a - 2.0 * Or) * np.sqrt(Or + Om * a) + 2.0 * Or**1.5)


def big_rip_approx_gyr(w: float, Om: float, H0: float) -> float:
    """Caldwell, Kamionkowski & Weinberg (2003): t_rip - t0 ~ (2/3) / (|1+w| H0 sqrt(1-Om))."""
    return 2.0 / 3.0 / (abs(1.0 + w) * np.sqrt(1.0 - Om)) * hubble_time_gyr(H0)


# --- The geography of the (Omega_m, Omega_Lambda) plane --------------------------


def double_root_curve(a_star, Or: float = 0.0, w: float = -1.0):
    """Universes balanced on a static state at size a_star: P(a*) = P'(a*) = 0.

    Both conditions are linear in (Om, OL) once Ok = 1 - Om - OL - Or is
    substituted, so the boundary is an exact parametric curve:

    * a_star < 1: the edge of the "no Big Bang" region. Above it universes
      bounce; on it they emerged from Einstein's static universe infinitely long
      ago; just below it they loiter at a ~ a_star before expanding again.
    * a_star > 1: the edge of the recollapsing region (for Om > 1): below it
      the expansion stops and reverses.

    For Lambda and no radiation this reproduces the closed forms of Felten &
    Isaacman (1986) and Carroll, Press & Turner (1992) (checked in the tests).
    """
    a = np.asarray(a_star, float)
    n = 1.0 - 3.0 * w
    # Om (a - a^2) + OL (a^n - a^2) = -Or (1 - a^2) - a^2
    # Om (1 - 2a) + OL (n a^(n-1) - 2a) = -2a (1 - Or)
    a11, a12 = a - a**2, a**n - a**2
    a21, a22 = 1.0 - 2.0 * a, n * a ** (n - 1.0) - 2.0 * a
    b1, b2 = -Or * (1.0 - a**2) - a**2, -2.0 * a * (1.0 - Or)
    det = a11 * a22 - a12 * a21
    with np.errstate(divide="ignore", invalid="ignore"):
        Om = (b1 * a22 - a12 * b2) / det
        OL = (a11 * b2 - a21 * b1) / det
    return Om, OL


def no_big_bang_boundary(Om_grid, Or: float = 0.0, w: float = -1.0) -> np.ndarray:
    """Omega_Lambda above which a universe with this Omega_m has no Big Bang.

    Inverts the a_star < 1 branch of ``double_root_curve`` (Om rises
    monotonically from 0 to infinity as a_star goes from its smallest physical
    value up to 1). With radiation the curve starts at Omega_m = 0 with
    Omega_Lambda ~ 1 + 2 sqrt(Or) and a_star ~ Or^(1/4): see ``max_bounce_redshift``.
    """
    a_lo = smallest_physical_static_a(Or, w)
    # dense near a_lo, where Om rises from 0 steeply, and up to a* ~ 0.95 (Om ~ 400)
    a = a_lo + (0.95 - a_lo) * np.linspace(0.0, 1.0, 200001) ** 2
    Om_c, OL_c = double_root_curve(a, Or, w)
    Om_c = np.where(np.abs(Om_c) < 1e-12, 0.0, Om_c)
    good = np.isfinite(Om_c) & (Om_c >= 0)
    order = np.argsort(Om_c[good])
    return np.interp(np.asarray(Om_grid, float), Om_c[good][order], OL_c[good][order],
                     left=np.nan, right=np.nan)


def recollapse_boundary(Om_grid, Or: float = 0.0, w: float = -1.0) -> np.ndarray:
    """Omega_Lambda below which a universe recollapses into a Big Crunch.

    Zero for Om <= 1 - Or (any negative Lambda eventually wins); above that,
    the a_star > 1 branch of ``double_root_curve``.
    """
    Om_grid = np.asarray(Om_grid, float)
    a = 1.0 + np.geomspace(1e-5, 1e6, 200001)
    Om_c, OL_c = double_root_curve(a, Or, w)
    good = np.isfinite(Om_c) & (Om_c > 0)
    order = np.argsort(Om_c[good])
    high = np.interp(Om_grid, Om_c[good][order], OL_c[good][order], left=0.0, right=np.nan)
    return np.where(Om_grid <= 1.0 - Or, 0.0, high)


def smallest_physical_static_a(Or: float, w: float = -1.0) -> float:
    """The static size a_star at which the no-Big-Bang boundary reaches Om = 0.

    For Lambda: Om(a*) = 0 when Or (1/a*^2 - 1)^2 = 1, i.e. a* = (1 + Or^-1/2)^-1/2.
    Without radiation it is 0.
    """
    if Or <= 0:
        return 1e-6
    if w == -1.0:
        return float((1.0 + Or**-0.5) ** -0.5)
    return float(optimize.brentq(lambda a: double_root_curve(a, Or, w)[0], 1e-6, 0.999))


def max_bounce_redshift(Or: float, Om_min: float = 0.0, w: float = -1.0) -> float:
    """The highest redshift any bouncing universe can ever show us.

    A universe that bounced at a_min never was smaller than a_min, so nothing in
    it can appear more redshifted than z = 1/a_min - 1. Over all bouncing
    universes with Omega_m >= Om_min, a_min is smallest on the boundary itself,
    where it equals a_star. So the answer is 1/a_star(Om_min) - 1.
    """
    if Om_min <= 0:
        a = smallest_physical_static_a(Or, w)
    else:
        a = optimize.brentq(lambda x: double_root_curve(x, Or, w)[0] - Om_min,
                            smallest_physical_static_a(Or, w) + 1e-12, 0.999999)
    return 1.0 / a - 1.0


def classify_grid(Om, OL, Or: float = 0.0) -> np.ndarray:
    """Vectorised fate of every (Om, OL) on a grid (Lambda): codes
    0 = expands forever and is accelerating today, 1 = expands forever but is
    decelerating today, 2 = recollapses (Big Crunch), 3 = no Big Bang (bounce).
    Checked point by point against ``Universe.fate`` in the tests."""
    Om, OL = np.broadcast_arrays(np.asarray(Om, float), np.asarray(OL, float))
    code = np.where(OL > 0.5 * Om + Or, 0, 1)
    code = np.where(OL < recollapse_boundary(Om, Or), 2, code)
    nbb = no_big_bang_boundary(Om, Or)
    code = np.where(np.isfinite(nbb) & (OL > nbb), 3, code)
    return code


def bounce_size_grid(Om, OL, Or: float = 0.0) -> np.ndarray:
    """a_min of every bouncing universe on a grid (NaN elsewhere), vectorised.

    For Lambda the bounce is the largest root below 1 of the quartic P(a);
    found by bisection on [a_star, 1], where P changes sign exactly once.
    """
    Om, OL = np.broadcast_arrays(np.asarray(Om, float), np.asarray(OL, float))
    Ok = 1.0 - Om - OL - Or
    nbb = no_big_bang_boundary(Om, Or)
    bounce = np.isfinite(nbb) & (OL > nbb)
    lo = np.zeros(Om.shape)
    hi = np.ones(Om.shape)
    # P(1) = 1 > 0, and P < 0 just above the bounce's lower partner root. Start
    # the bracket at the minimum of P on (0, 1), which lies between the two roots.
    P = lambda a: Or + Om * a + Ok * a**2 + OL * a**4  # noqa: E731
    with np.errstate(invalid="ignore", divide="ignore"):
        grid = np.linspace(1e-4, 1.0, 400)
        vals = np.stack([P(g) for g in grid])
        lo = grid[np.argmin(vals, axis=0)]
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        neg = P(mid) < 0
        lo = np.where(neg, mid, lo)
        hi = np.where(neg, hi, mid)
    return np.where(bounce, hi, np.nan)


def age_grid(Om, OL, Or: float = 0.0, n_nodes: int = 400) -> np.ndarray:
    """t0 H0 for every Big Bang universe on a grid (NaN for bounces), vectorised.

    Gauss-Legendre in u = sqrt(a), t0 H0 = int_0^1 2u^3 du / sqrt(P(u^2)). Accurate
    to ~1e-8 away from the loitering edge, where the integrand spikes and the
    age diverges (there the value is a lower bound).
    """
    Om, OL = np.broadcast_arrays(np.asarray(Om, float), np.asarray(OL, float))
    Ok = 1.0 - Om - OL - Or
    x, wts = np.polynomial.legendre.leggauss(n_nodes)
    u = 0.5 * (x + 1.0)
    wts = 0.5 * wts
    out = np.zeros(Om.shape)
    bad = np.zeros(Om.shape, bool)
    for ui, wi in zip(u, wts):
        a = ui * ui
        p = Or + Om * a + Ok * a * a + OL * a**4
        bad |= p <= 0
        with np.errstate(invalid="ignore", divide="ignore"):
            out += wi * 2.0 * ui**3 / np.sqrt(np.where(p > 0, p, np.nan))
    out[bad] = np.nan
    return out
