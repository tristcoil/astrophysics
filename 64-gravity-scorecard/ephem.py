"""JPL DE440 positions and the Moon's orientation, evaluated at picosecond-exact times.

Laser ranging needs the Earth and Moon at three instants (fire, bounce, return)
whose differences must be good to ~1 ps. A Julian date in one double resolves
only ~40 us, and jplephem's evaluator forms "seconds since segment start" as a
difference of two ~5e9 s numbers (~0.6 us of rounding, i.e. ~1-2 cm for the
Earth's 30 km/s). So times here are (day, seconds) pairs counted from J2000 TDB,
and the Chebyshev offset is built with exact integer arithmetic before the
seconds are added. The coefficients themselves are read with jplephem.

Positions are in km and velocities in km/s, ICRF axes, barycentric (SSB) or as
requested. The lunar libration angles (phi, theta, psi) of the DE440
principal-axis frame give r_ICRF = Rz(-phi) Rx(-theta) Rz(-psi) r_PA
(Park et al. 2021, eq. 8, with frame-rotation matrices).
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from jplephem.pck import PCK
from jplephem.spk import SPK

import config as C

SECONDS_PER_DAY = 86_400.0


class Chebyshev:
    """One SPK/PCK type-2 segment, evaluated with exact (day, second) times."""

    def __init__(self, segment):
        init, intlen, coefficients = segment._data if segment._data is not None else segment._load()
        segment._data = (init, intlen, coefficients)
        self.init = float(init)          # seconds past J2000 (an integer number of seconds)
        self.intlen = float(intlen)
        self.coef = coefficients         # (n_coef, n_comp, n_intervals); coef[-1] is the constant term
        self.n = coefficients.shape[2]
        if self.init != round(self.init) or self.intlen != round(self.intlen):
            raise ValueError("segment start or interval length is not an integer number of seconds")

    def __call__(self, day, sec, derivative: bool = True):
        """Value (and d/dt per second) at t = day*86400 + sec seconds past J2000 TDB."""
        day = np.atleast_1d(np.asarray(day, dtype=np.int64))
        sec = np.atleast_1d(np.asarray(sec, dtype=float))
        base = day * 86_400 - np.int64(self.init)                # exact integer seconds
        index = base // np.int64(self.intlen)
        offset = (base - index * np.int64(self.intlen)).astype(float) + sec
        shift = np.floor(offset / self.intlen).astype(np.int64)
        index = index + shift
        offset = offset - shift * self.intlen
        at_end = index == self.n
        index[at_end] -= 1
        offset[at_end] += self.intlen
        if (index < 0).any() or (index >= self.n).any():
            raise ValueError("time outside the ephemeris segment")
        cf = self.coef[:, :, index]
        s = 2.0 * offset / self.intlen - 1.0
        s2 = 2.0 * s
        w0 = w1 = dw0 = dw1 = 0.0
        for c in cf[:-1]:
            w2, w1 = w1, w0
            w0 = c + (s2 * w1 - w2)
            if derivative:
                dw2, dw1 = dw1, dw0
                dw0 = 2.0 * w1 + dw1 * s2 - dw2
        value = cf[-1] + (s * w0 - w1)
        if not derivative:
            return value
        rate = (w0 + s * dw0 - dw1) * (2.0 / self.intlen)
        return value, rate


class Ephemeris:
    def __init__(self, spk_path=None, pck_path=None):
        self.spk = SPK.open(str(spk_path or C.NAIF_DIR / "de440s.bsp"))
        self.pck = PCK.open(str(pck_path or C.NAIF_DIR / "moon_pa_de440_200625.bpc"))
        self._seg = {}
        for seg in self.spk.segments:
            self._seg[(seg.center, seg.target)] = Chebyshev(seg)
        self._lib = [Chebyshev(seg) for seg in self.pck.segments]
        self._lib_range = [(s.init, s.init + s.intlen * s.n) for s in self._lib]

    def _chain(self, target: int):
        """SSB -> target as a list of (center, target) segment keys."""
        if target == 0:
            return []
        for (c, t) in self._seg:
            if t == target:
                return self._chain(c) + [(c, t)]
        raise KeyError(target)

    def posvel(self, target: int, day, sec, center: int = 0):
        """Position (km) and velocity (km/s) of target relative to center (NAIF ids)."""
        def ssb(body):
            p = v = 0.0
            for key in self._chain(body):
                pp, vv = self._seg[key](day, sec)
                p = p + pp
                v = v + vv
            return np.asarray(p), np.asarray(v)
        pt, vt = ssb(target)
        if center == 0:
            return pt, vt
        pc, vc = ssb(center)
        return pt - pc, vt - vc

    def moon_angles(self, day, sec):
        """Libration Euler angles (phi, theta, psi) in rad and their rates (rad/s)."""
        t = np.atleast_1d(day).astype(float) * SECONDS_PER_DAY + np.atleast_1d(sec)
        for seg, (a, b) in zip(self._lib, self._lib_range):
            if np.all((t >= a) & (t <= b)):
                return seg(day, sec)
        raise ValueError("times span more than one libration segment")

    def moon_pa_to_icrf(self, day, sec):
        """Rotation matrices (N, 3, 3) taking DE440 PA-frame vectors to ICRF."""
        ang = self.moon_angles(day, sec)[0]
        phi, theta, psi = ang
        return rz(-phi) @ rx(-theta) @ rz(-psi)


def rx(a):
    """Frame-rotation matrix about x (Park et al. 2021 eq. 9), stacked over a."""
    a = np.atleast_1d(a)
    c, s = np.cos(a), np.sin(a)
    m = np.zeros(a.shape + (3, 3))
    m[..., 0, 0] = 1
    m[..., 1, 1] = c
    m[..., 1, 2] = s
    m[..., 2, 1] = -s
    m[..., 2, 2] = c
    return m


def rz(a):
    a = np.atleast_1d(a)
    c, s = np.cos(a), np.sin(a)
    m = np.zeros(a.shape + (3, 3))
    m[..., 0, 0] = c
    m[..., 0, 1] = s
    m[..., 1, 0] = -s
    m[..., 1, 1] = c
    m[..., 2, 2] = 1
    return m


@lru_cache(maxsize=1)
def default() -> Ephemeris:
    return Ephemeris()


# ---------------------------------------------------------------------------
# EPM2021 (IAA RAS): an independent lunar ephemeris for cross-checks
# ---------------------------------------------------------------------------

EPM_SLICE = C.RAW_DIR / "epm2021" / "epm2021_de_ascii_2005-2021.part"
EPM_CACHE = C.DERIVED_DIR / "epm2021_records_2005_2021.npy"


def load_de_ascii(path) -> np.ndarray:
    """Records (each 1018 numbers: JD start, JD end, coefficients) from a JPL-format ASCII
    ephemeris file or a byte-range slice of one (partial records at the edges are dropped)."""
    recs, cur, n = [], None, None
    for line in open(path, errors="replace"):
        t = line.split()
        if len(t) == 2 and all(x.isdigit() for x in t):
            if cur is not None and len(cur) >= n:
                recs.append(np.array(cur[:n], float))
            n, cur = int(t[1]), []
        elif cur is not None:
            cur += [float(x.replace("D", "E")) for x in t]
    if cur is not None and len(cur) >= n:
        recs.append(np.array(cur[:n], float))
    return np.array(recs)


class DEAscii:
    """Geocentric Moon (item 10) and lunar librations (item 13) from JPL-format records,
    evaluated with exact (day, second) arithmetic like ``Chebyshev``."""

    ITEMS = {"moon": (441, 13, 8, 3), "librations": (899, 10, 4, 3)}   # pointer, ncoef, nsub, ncomp

    def __init__(self, recs: np.ndarray):
        self.recs = recs
        self.start = np.round((recs[:, 0] - 2451545.0) * SECONDS_PER_DAY).astype(np.int64)
        self.span = int(round((recs[0, 1] - recs[0, 0]) * SECONDS_PER_DAY))

    def __call__(self, what: str, day, sec):
        ptr, ncoef, nsub, ncomp = self.ITEMS[what]
        day = np.atleast_1d(np.asarray(day, dtype=np.int64))
        sec = np.atleast_1d(np.asarray(sec, dtype=float))
        t_int = day * 86_400 + np.floor(sec).astype(np.int64)
        idx = np.searchsorted(self.start, t_int, side="right") - 1
        if (idx < 0).any() or (t_int - self.start[idx] > self.span).any():
            raise ValueError("time outside the EPM2021 slice")
        off = (day * 86_400 - self.start[idx]).astype(float) + sec
        sub = self.span / nsub
        k = np.minimum(np.floor(off / sub).astype(int), nsub - 1)
        s = 2.0 * (off - k * sub) / sub - 1.0
        rows = self.recs[idx]                                           # (N, 1018)
        base = ptr - 1 + (k[:, None] * ncomp + np.arange(ncomp)[None, :]) * ncoef     # (N, ncomp)
        cf = np.take_along_axis(rows[:, None, :].repeat(ncomp, axis=1),
                                base[:, :, None] + np.arange(ncoef)[None, None, :], axis=2)  # (N, ncomp, ncoef)
        # Chebyshev values and derivatives by recurrence
        tk = np.zeros((ncoef,) + s.shape)
        dk = np.zeros_like(tk)
        tk[0] = 1.0
        if ncoef > 1:
            tk[1] = s
            dk[1] = 1.0
        for j in range(2, ncoef):
            tk[j] = 2 * s * tk[j - 1] - tk[j - 2]
            dk[j] = 2 * tk[j - 1] + 2 * s * dk[j - 1] - dk[j - 2]
        val = np.einsum("nck,kn->cn", cf, tk)
        der = np.einsum("nck,kn->cn", cf, dk) * (2.0 / sub)
        return val, der


class MixedEphemeris(Ephemeris):
    """DE440 for the Sun and planets, EPM2021 for the geocentric Moon and its librations.

    Only the Earth-Moon vector and the Moon's orientation matter for laser ranging,
    so this isolates the lunar solution. (EPM2021's principal-axis frame differs from
    DE440's, so reflector coordinates must be refitted when using it.)
    """

    def __init__(self, recs: np.ndarray | None = None):
        super().__init__()
        if recs is None:
            if EPM_CACHE.exists():
                recs = np.load(EPM_CACHE)
            else:
                recs = load_de_ascii(EPM_SLICE)
                np.save(EPM_CACHE, recs)
        self.epm = DEAscii(recs)

    def posvel(self, target: int, day, sec, center: int = 0):
        if target == MOON:
            pe, ve = super().posvel(EARTH, day, sec)
            gm, gv = self.epm("moon", day, sec)
            pm, vm = pe + gm, ve + gv
            if center == 0:
                return pm, vm
            pc, vc = super().posvel(center, day, sec) if center != MOON else (pm, vm)
            return pm - pc, vm - vc
        if center == MOON:
            pt, vt = self.posvel(target, day, sec)
            pm, vm = self.posvel(MOON, day, sec)
            return pt - pm, vt - vm
        return super().posvel(target, day, sec, center)

    def moon_angles(self, day, sec):
        return self.epm("librations", day, sec)


@lru_cache(maxsize=1)
def epm2021() -> MixedEphemeris:
    return MixedEphemeris()


SUN, EMB, EARTH, MOON = 10, 3, 399, 301
