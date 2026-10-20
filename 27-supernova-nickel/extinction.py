"""Interstellar dust extinction: Cardelli, Clayton & Mathis (1989) with O'Donnell (1994).

A_lambda / A_V = a(x) + b(x) / R_V, with x = 1/lambda in inverse microns. R_V
sets how steeply extinction rises toward the blue: 3.1 for the average Milky
Way, but much lower (1.4) for the dust in front of SN 2014J. Valid from 0.1 to
3.3 microns.
"""

from __future__ import annotations

import numpy as np


def ccm89(wave_angstrom, r_v: float = 3.1):
    """A_lambda / A_V for the CCM89 law with O'Donnell's (1994) optical coefficients."""
    x = 1e4 / np.atleast_1d(np.asarray(wave_angstrom, dtype=float))
    a = np.zeros_like(x)
    b = np.zeros_like(x)

    ir = (x >= 0.3) & (x < 1.1)
    a[ir] = 0.574 * x[ir] ** 1.61
    b[ir] = -0.527 * x[ir] ** 1.61

    opt = (x >= 1.1) & (x < 3.3)
    y = x[opt] - 1.82
    a[opt] = np.polyval([-0.505, 1.647, -0.827, -1.718, 1.137, 0.701, -0.609, 0.104, 1.0], y)
    b[opt] = np.polyval([3.347, -10.805, 5.491, 11.102, -7.985, -3.989, 2.908, 1.952, 0.0], y)

    uv = (x >= 3.3) & (x < 8.0)
    xu = x[uv]
    fa = np.where(xu >= 5.9, -0.04473 * (xu - 5.9) ** 2 - 0.009779 * (xu - 5.9) ** 3, 0.0)
    fb = np.where(xu >= 5.9, 0.2130 * (xu - 5.9) ** 2 + 0.1207 * (xu - 5.9) ** 3, 0.0)
    a[uv] = 1.752 - 0.316 * xu - 0.104 / ((xu - 4.67) ** 2 + 0.341) + fa
    b[uv] = -3.090 + 1.825 * xu + 1.206 / ((xu - 4.62) ** 2 + 0.263) + fb

    fuv = (x >= 8.0) & (x <= 11.0)
    xf = x[fuv] - 8.0
    a[fuv] = -1.073 - 0.628 * xf + 0.137 * xf**2 - 0.070 * xf**3
    b[fuv] = 13.670 + 4.257 * xf - 0.420 * xf**2 + 0.374 * xf**3

    if np.any((x < 0.3) | (x > 11.0)):
        raise ValueError("CCM89 is defined for 0.091-3.33 micron only")
    return a + b / r_v


def a_lambda(wave_angstrom, ebv: float, r_v: float = 3.1):
    """Extinction in magnitudes at the given wavelengths for E(B-V) and R_V."""
    return ebv * r_v * ccm89(wave_angstrom, r_v)
