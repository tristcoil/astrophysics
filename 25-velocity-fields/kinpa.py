"""Kinematic position angles: which way is a velocity field rotating?

The standard tool is ``fit_kinematic_pa`` (Krajnovic et al. 2006, Appendix C;
Cappellari's ``pafit`` package), which Zhou et al. (2022) used on MaNGA. For
each trial axis it builds the bi-antisymmetric version of the map a rotating
disc would have -- reflect it about the axis and its perpendicular, flipping
the sign where rotation demands -- and keeps the axis where the map agrees
best with its own symmetrized copy.

``pafit`` re-triangulates the bins for each of its 361 trial angles. A
Delaunay triangulation does not change when the points are rotated, so here
it is built once, and the reflected points are mapped back to the original
frame instead. The arithmetic is otherwise ``pafit``'s, line for line
(including its error: the half-range of angles with chi^2 within
9 + 3 sqrt(2N) of the minimum, which ``pafit`` labels a 3-sigma error), and
``test_physics.py`` checks that both give the same answers.

``pafit`` returns an *axis* (0-180 deg): the sense of rotation is left to the
user. To tell co- from counter-rotation we add it: the receding side is the
end of the axis where the velocities are positive.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay


def _rotate(x, y, angle_deg):
    t = np.radians(angle_deg)
    c, s = np.cos(t), np.sin(t)
    return x * c - y * s, x * s + y * c


def symmetrized(interp, x, y, vel, pa):
    """pafit's bi-antisymmetric (sym=1) velocity field for axis ``pa``."""
    xr, yr = _rotate(x, y, -(pa + 90))
    out = [vel]
    for sx, sy, sign in ((-1, 1, -1.0), (1, -1, 1.0), (-1, -1, -1.0)):
        bx, by = _rotate(sx * xr, sy * yr, pa + 90)  # reflected point, back in the original frame
        out.append(sign * interp(bx, by))
    return np.nanmean(np.vstack(out), axis=0)


def fit_kinematic_pa(x, y, vel, dvel, nsteps=361) -> tuple[float, float, float, np.ndarray]:
    """Kinematic axis (deg, 0-180), its error (pafit's convention) and V_sys offset.

    ``vel`` must already have an estimate of the systemic velocity removed.
    Also returns chi^2 at every trial angle.
    """
    x, y, vel = (np.asarray(a, float).ravel() for a in (x, y, vel))
    dvel = np.broadcast_to(np.asarray(dvel, float), vel.shape)
    interp = LinearNDInterpolator(Delaunay(np.column_stack([x, y])), vel)
    angles = np.linspace(0, 180, nsteps)
    chi2 = np.array([np.sum(((vel - symmetrized(interp, x, y, vel, a)) / dvel) ** 2) for a in angles])
    k = int(np.argmin(chi2))
    best = float(angles[k])
    f = chi2 - chi2[k] <= 9 + 3 * np.sqrt(2 * vel.size)
    min_err = max(0.5, 45.0 / (nsteps - 1))
    if f.sum() > 1:
        err = (angles[f].max() - angles[f].min()) / 2.0
        if err >= 45:
            good = np.degrees(np.arctan(np.tan(np.radians(angles[f]))))
            err = (good.max() - good.min()) / 2.0
        err = max(err, min_err)
    else:
        err = min_err
    vsyst = float(np.median(vel - symmetrized(interp, x, y, vel, best)))
    return best, float(err), vsyst, chi2


def receding_pa(x, y, vel, axis_pa) -> float:
    """Turn an axis into the position angle of the receding side (0-360).

    pafit draws its major axis at (-sin PA, cos PA); the receding end is the
    one along which the velocities are positive.
    """
    t = np.radians(axis_pa)
    along = -np.asarray(x) * np.sin(t) + np.asarray(y) * np.cos(t)
    return float(axis_pa % 360 if np.sum(vel * along) >= 0 else (axis_pa + 180) % 360)


def angle_between(pa1, pa2) -> float:
    """Smallest angle between two directions, 0-180 deg."""
    d = abs((pa1 - pa2) % 360)
    return float(min(d, 360 - d))
