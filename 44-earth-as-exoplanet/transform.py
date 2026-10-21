"""From display brightness back to light: what we know about the EPIC colour pipeline.

The natural-colour images are made from EPIC's calibrated 680, 551 and 443 nm
reflectances, then "color and brightness adjusted to represent what a
conventional camera would produce" (EPIC site). The adjustment is not
published, so we model it as a per-channel power law,

    V / 255 = (c_ch * R)^(1 / gamma),

and test gamma six ways (``run.py`` re-derives each; numbers from the final
run, 20-31 July 2016 unless stated):

1. Partly cloudy pixels are linear mixtures (in light) of ocean and cloud, so
   in display space their blue-vs-red and green-vs-red loci bend by an amount
   set by gamma alone: gamma = 2.36 (blue) and 2.29 (green), residuals 0.5-1
   display units; the implied clear-ocean colours (blue/red 7.7x cloud's,
   green/red 2.8x) are close to Rayleigh scattering's (~5.5 and ~2.4).
2. Ocean vs cloud colour (blue/red), with Rayleigh physics for the ocean
   (reflectance ratio 4.5-7): gamma = 1.4-2.0.
3. Day-to-day relative scatter of the disc mean against Yang et al. (2018),
   who used calibrated L1B data for the same 12 days: gamma = 1.67 (680 nm)
   and 1.74 (443 nm).
4. Their peak-to-peak swings: gamma = 2.4 (680 nm); 443 nm needs > 2.6.
5. The darkest ocean pixels (top-of-atmosphere reflectance ~0.02 at 680 nm
   and ~0.12 at 443 nm) against the disc mean: gamma = 1.5-1.8 and 1.2-1.3.
6. The Moon in the lunar transits against EPIC's calibrated far-side
   reflectance (``moon.py``): gamma = 1.1-1.2 (443 nm), below 1 (680 nm).

The tests disagree beyond their errors, so the curve is not one power law;
the frame-to-frame tone adjustments (``lightcurve.frame_jump``) show it also
varies from image to image. What this project reports is mostly disc-
integrated brightness, for which test 3 is the direct calibration: we adopt
gamma = 1.8 and carry 1.4-2.4 as the systematic range (the pixel-level
curvature of test 1 at the top, the ocean anchors at the bottom). The Moon
(test 6) falls outside it; see the README.

The per-channel scale c_ch is tied to Yang et al.'s calibrated disc-mean
reflectances for the same July 2016 days (0.27 at 443 nm, 0.21 at 680 nm; the
551 nm band has no published value, so green stays in relative units).
"""

from __future__ import annotations

import numpy as np

GAMMA_NOMINAL = 1.8
GAMMA_RANGE = (1.4, 2.4)
#: Exponents stored per image, so every result can be recomputed at the ends
#: of the systematic range without re-reading images.
GAMMA_GRID = (1.0, 1.4, 1.8, 2.0, 2.2, 2.4, 2.6)

#: Yang, Marshak, Varnai & Knyazikhin 2018, Remote Sensing 10, 254: simple
#: disc mean of calibrated L1B reflectance (their <R>_1, pixels with solar
#: zenith angle < 78 deg), 20-31 July 2016, and its day-to-day variability.
YANG2018 = {
    "window": ("2016-07-20", "2016-07-31"),
    "mean": {"blue": 0.27, "red": 0.21},
    "daily_rel_std": {"blue": 0.04, "red": 0.08},
    "daily_peak_to_peak": {"blue": 0.23, "red": 0.33},
    # Same quantities with the cos(view zenith) weighting, their <R>_2.
    "mean_R2": {"blue": 0.28, "red": 0.22},
}

#: Clear-ocean top-of-atmosphere reflectance near the disc centre, from
#: single Rayleigh scattering (optical depth 0.236 at 443 nm, 0.037 at 680 nm,
#: phase function 1.5 at backscatter) plus water and aerosol: an anchor good to
#: ~30%.
OCEAN_REFLECTANCE = {"blue": 0.12, "red": 0.02}


def levels() -> np.ndarray:
    """Display values 0..255 as fractions."""
    return np.arange(256) / 255.0


def to_linear(v, gamma: float = GAMMA_NOMINAL):
    """Display value(s) in 0..255 -> relative linear light (0..1)."""
    return (np.asarray(v, float) / 255.0) ** gamma


def to_display(lin, gamma: float = GAMMA_NOMINAL):
    """Relative linear light (0..1) -> display value 0..255."""
    return 255.0 * np.clip(np.asarray(lin, float), 0.0, 1.0) ** (1.0 / gamma)


def gamma_index(gamma: float) -> int:
    return GAMMA_GRID.index(gamma)
