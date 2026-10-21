"""The Cassini navigation team's SCE1 short-arc trajectory (SPK "SM727D").

The SCE1 archive (CORS_0021, sce1_ancillary/spk) carries the navigation
team's own trajectory for the experiment, "SM727D_SHORTARC_SCE1.BSP"
(J. Bordi, 5 August 2002): fitted to DSS 25's plasma-calibrated two-way
X-band Doppler of 2002 June 9 - July 5, estimating the spacecraft state, the
RTG thermal acceleration and solar pressure. It is the trajectory the
radio-science analysis started from, so it is the natural alternative to
Horizons' reconstruction (``cassini_doppler.Ephemeris(cassini="nav")``).

The file is a SPICE transfer file (ASCII "DAFETF") holding SPK type 1
segments (modified difference arrays) of Cassini (-82) about the Sun (10) in
J2000 = ICRF axes, TDB. This module decodes the transfer format and
evaluates type 1 records with the algorithm of SPICE's SPKE01, with no SPICE
dependency.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

import config as C

SPK_PATH = C.RAW_DIR / "cassini_sce1" / "spk" / "c32easc2002_160_2002_186.spk"
#: J2000 TDB epoch (JD 2451545.0) to the project's T0 (2002-06-04 00:00 TDB), s
T0_FROM_J2000 = (2452429.5 - 2451545.0) * 86400.0
#: The SPK's fitted arc ends with the data on July 5; on July 6 (day 32.6
#: after T0) its records turn short and the trajectory departs from Horizons'
#: (a modelled event): use it only up to here (days after T0).
NAV_END_DAY = 32.5


def _decode(tok: str) -> float:
    """A DAF transfer-format number: signed hex mantissa (an implied leading
    '0.') and a hex exponent of 16, e.g. '-2EA91B454B6242^-7'; plain integers
    are hex too ('-52' = -82)."""
    tok = tok.strip().strip("'")
    if "^" not in tok:
        return float(int(tok, 16))
    m, e = tok.split("^")
    sign = -1.0 if m.startswith("-") else 1.0
    m = m.lstrip("-")
    return sign * int(m, 16) * 16.0 ** (int(e, 16) - len(m))


def read_transfer(path: Path = SPK_PATH) -> list[dict]:
    """The segments of an SPK transfer file: summary and data values."""
    lines = Path(path).read_text().splitlines()
    segs, i = [], 0
    while i < len(lines):
        m = re.match(r"BEGIN_ARRAY\s+(\d+)\s+(\d+)", lines[i])
        if not m:
            i += 1
            continue
        n_values = int(m.group(2))
        name = lines[i + 1].strip().strip("'").strip()
        summ = [_decode(lines[i + 2 + k]) for k in range(6)]
        j = i + 8
        data = []
        while not lines[j].startswith("END_ARRAY"):
            s = lines[j].strip()
            if s.startswith("'"):
                data.append(_decode(s))
            j += 1             # block counts (bare integers) are skipped
        if len(data) != n_values:
            raise ValueError(f"segment {m.group(1)}: {len(data)} values, expected {n_values}")
        segs.append({"name": name, "et0": summ[0], "et1": summ[1], "target": int(summ[2]),
                     "center": int(summ[3]), "frame": int(summ[4]), "type": int(summ[5]),
                     "data": np.array(data)})
        i = j + 1
    return segs


class Type1Segment:
    """SPK type 1 (modified difference arrays), evaluated as in SPICE SPKE01."""

    def __init__(self, seg: dict):
        d = seg["data"]
        n = int(round(d[-1]))
        self.n = n
        self.records = d[:71 * n].reshape(n, 71)
        self.epochs = d[71 * n:72 * n]
        self.et0, self.et1 = seg["et0"], seg["et1"]

    def state(self, et: float) -> np.ndarray:
        i = int(np.searchsorted(self.epochs, et))
        i = min(i, self.n - 1)
        rec = self.records[i]
        tl = rec[0]
        g = rec[1:16]
        refpos = rec[[16, 18, 20]]
        refvel = rec[[17, 19, 21]]
        dt = rec[22:67].reshape(3, 15).T               # DT(15, 3), Fortran order
        kqmax1 = int(rec[67])
        kq = rec[68:71].astype(int)
        delta = et - tl
        tp = delta
        mq2 = kqmax1 - 2
        ks = kqmax1 - 1
        fc = np.zeros(17)
        wc = np.zeros(17)
        w = np.zeros(20)
        fc[0] = 1.0
        for j in range(1, mq2 + 1):                     # Fortran J = 1..MQ2
            fc[j] = tp / g[j - 1]                       # FC(J+1)
            wc[j - 1] = delta / g[j - 1]                # WC(J)
            tp = delta + g[j - 1]
        for j in range(1, kqmax1 + 1):
            w[j - 1] = 1.0 / j                          # W(J)
        jx = 0
        ks1 = ks - 1
        while ks >= 2:
            jx += 1
            for j in range(1, jx + 1):
                # W(J+KS) = FC(J+1)*W(J+KS-1) - WC(J)*W(J+KS)
                w[j + ks - 1] = fc[j] * w[j + ks - 2] - wc[j - 1] * w[j + ks - 1]
            ks = ks1
            ks1 -= 1
        out = np.zeros(6)
        for i3 in range(3):
            s = 0.0
            for j in range(kq[i3], 0, -1):
                s += dt[j - 1, i3] * w[j + ks - 1]
            out[i3] = refpos[i3] + delta * (refvel[i3] + delta * s)
        for j in range(1, jx + 1):
            w[j + ks - 1] = fc[j] * w[j + ks - 2] - wc[j - 1] * w[j + ks - 1]
        ks -= 1
        for i3 in range(3):
            s = 0.0
            for j in range(kq[i3], 0, -1):
                s += dt[j - 1, i3] * w[j + ks - 1]
            out[3 + i3] = refvel[i3] + delta * s
        return out


class NavTrajectory:
    """Cassini about the Sun (m, m/s; J2000/ICRF; TDB) from the SM727D SPK,
    tabulated every 60 s (positions: Hermite; velocities: cubic spline of
    the exact velocities); t = TDB s since T0."""

    STEP = 60.0

    def __init__(self, path: Path = SPK_PATH):
        from scipy.interpolate import CubicHermiteSpline, CubicSpline

        segs = [Type1Segment(s) for s in read_transfer(path) if s["type"] == 1 and s["target"] == -82]
        if not segs:
            raise ValueError("no Cassini type-1 segment")
        # the Sun-centred segment: the short arc and its prediction; stop
        # before the modelled maneuver of 2002 July 6 (many short records)
        et0 = min(s.et0 for s in segs)
        et1 = min(max(s.et1 for s in segs), T0_FROM_J2000 + NAV_END_DAY * 86400.0)
        t = np.arange(np.ceil((et0 - T0_FROM_J2000) / self.STEP) * self.STEP,
                      (et1 - T0_FROM_J2000), self.STEP)
        st = np.zeros((len(t), 6))
        for k, tk in enumerate(t):
            et = tk + T0_FROM_J2000
            seg = next(s for s in reversed(segs) if s.et0 <= et <= s.et1)
            st[k] = seg.state(et)
        self.t, self.states = t, st * 1e3                 # km -> m
        self.span = (t[0], t[-1])
        self._s = CubicHermiteSpline(t, self.states[:, :3], self.states[:, 3:], axis=0)
        self._v = CubicSpline(t, self.states[:, 3:], axis=0)

    def _mask(self, t, x):
        t = np.asarray(t)
        out = (t < self.span[0]) | (t > self.span[1])
        if np.any(out):
            x = np.array(x, dtype=float)
            x[out] = np.nan
        return x

    def pos(self, t) -> np.ndarray:
        """NaN outside the SPK's span (2002 June 9 15:50 - July 6)."""
        return self._mask(t, self._s(t))

    def vel(self, t) -> np.ndarray:
        return self._mask(t, self._v(t))


if __name__ == "__main__":
    for s in read_transfer():
        print(s["name"], s["target"], s["center"], s["frame"], s["type"], s["et0"], s["et1"], len(s["data"]))
