"""Extinction to any star from the Edenhofer et al. (2024) 3D dust map: a dereddening tool.

    python deredden.py 184.6 -5.8 140            # l, b in degrees, distance in pc
    python deredden.py --radec 83.82 -5.39 400   # RA, Dec (ICRS, degrees), distance in pc
    python deredden.py --parallax 56.67 24.12 7.35   # RA, Dec, parallax in mas

From Python (vectorised; this is what B1 and B8 can import):

    from deredden import Extinction
    ext = Extinction()                          # loads the map (1.6 GB, a few seconds)
    ext.av(l, b, d)                             # A_V at 540 nm, magnitudes
    ext.bands(l, b, d)                          # dict of A at several wavelengths

What it computes: the map's voxels are summed from the Sun to the star (each
voxel constant between its radial boundaries, plus the map's own integral
inside 69 pc), and the four nearest HEALPix sightlines are combined with
healpy's bilinear weights. Uncertainty: a range from the map's per-voxel
standard deviations, added in quadrature along the line of sight (as if
voxels were independent: a lower bound) and added linearly (as if fully
correlated: an upper bound). In the released posterior samples neighbouring
voxels correlate by +0.3 to +0.5, so the truth lies in between.

Units: the map's "ZGR23" extinction E times the Zhang, Green & Rix (2023)
extinction curve gives the monochromatic extinction at any wavelength; we
quote A_V at 540 nm (2.80 E, the map authors' convention), A_G (2.04 E, at
672 nm, as O'Neill et al. 2024), E(B-V) ~ A(440 nm) - A(550 nm) = 0.97 E,
E(BP-RP) ~ A(517 nm) - A(782 nm) = 1.36 E, and the 2MASS J, H, Ks values the
curve tabulates. Broad-band extinctions depend on the star's own spectrum by
up to ~10%: treat these as monochromatic approximations.

Limits: 69 pc to 1.25 kpc (inside 69 pc the inner integral is scaled
linearly; beyond 1.25 kpc the value is the map's full column and marked).
"""

from __future__ import annotations

import argparse
from functools import cached_property

import numpy as np

import config as C
from dustmap import DustSphere


def zgr23_curve():
    """The Zhang, Green & Rix (2023) extinction curve: wavelength (nm) -> A_lambda / E."""
    data = np.loadtxt(C.CATALOG_DIR / "zgr23" / "extinction_curve.txt", skiprows=1)
    return data[:, 0], data[:, 1]


def coefficient(wavelength_nm: float) -> float:
    lam, curve = zgr23_curve()
    return float(np.interp(wavelength_nm, lam, curve))


class Extinction:
    """Integrated extinction from the Sun along (l, b) to distance d."""

    def __init__(self, which: str = "main"):
        self.sphere = DustSphere(which)

    @cached_property
    def _std_cumulative(self) -> np.ndarray:
        """Cumulative linear sum of voxel standard deviations (fully correlated voxels)."""
        s = self.sphere
        cum = np.empty((s.nshell + 1, s.npix), dtype=np.float32)
        acc = s.inner_std.astype(np.float64)
        cum[0] = acc
        widths = np.diff(s.bounds)
        for i in range(s.nshell):
            acc += widths[i] * s.std[i]
            cum[i + 1] = acc
        return cum

    @cached_property
    def _var_cumulative(self) -> np.ndarray:
        """Cumulative variance to each radial boundary (voxels independent)."""
        s = self.sphere
        cum = np.empty((s.nshell + 1, s.npix), dtype=np.float32)
        acc = s.inner_std.astype(np.float64) ** 2
        cum[0] = acc
        widths = np.diff(s.bounds)
        for i in range(s.nshell):
            acc += (widths[i] * s.std[i]) ** 2
            cum[i + 1] = acc
        return cum

    def zgr23(self, l, b, d) -> np.ndarray:
        return self.sphere.extinction(l, b, d)

    def zgr23_err(self, l, b, d) -> tuple[np.ndarray, np.ndarray]:
        """1-sigma error range: (independent voxels, fully correlated voxels)."""
        import healpy as hp
        s = self.sphere
        l, b, d = (np.atleast_1d(np.asarray(a, dtype=float)) for a in np.broadcast_arrays(l, b, d))
        pix, w = hp.get_interp_weights(s.nside, l, b, nest=True, lonlat=True)
        k = np.clip(np.searchsorted(s.bounds, d, side="right") - 1, 0, s.nshell - 1)
        part = np.clip(d - s.bounds[k], 0, None)
        var = self._var_cumulative[k[None, :], pix] + (s.std[k[None, :], pix] * part) ** 2
        lin = self._std_cumulative[k[None, :], pix] + s.std[k[None, :], pix] * part
        return np.sqrt((var * w).sum(0)), (lin * w).sum(0)

    def av(self, l, b, d) -> np.ndarray:
        return C.ZGR23_TO_AV * self.zgr23(l, b, d)

    def bands(self, l, b, d) -> dict:
        e = self.zgr23(l, b, d)
        return {"E_zgr23": e, "A_V(540nm)": C.ZGR23_TO_AV * e, "A_G(672nm)": C.ZGR23_TO_AG * e,
                "E(B-V)~A440-A550": (coefficient(440) - coefficient(550)) * e,
                "E(BP-RP)~A517-A782": (coefficient(517) - coefficient(782)) * e,
                "A_J(1.24um)": coefficient(1235) * e, "A_H(1.66um)": coefficient(1662) * e,
                "A_Ks(2.16um)": coefficient(2159) * e}


def _parse():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("a", type=float, help="l (deg), or RA with --radec/--parallax")
    ap.add_argument("b", type=float, help="b (deg), or Dec with --radec/--parallax")
    ap.add_argument("d", type=float, help="distance (pc), or parallax (mas) with --parallax")
    ap.add_argument("--radec", action="store_true", help="input is ICRS RA, Dec")
    ap.add_argument("--parallax", action="store_true", help="input is RA, Dec, parallax (mas)")
    return ap.parse_args()


def main() -> None:
    args = _parse()
    l, b, d = args.a, args.b, args.d
    if args.radec or args.parallax:
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        g = SkyCoord(ra=l * u.deg, dec=b * u.deg, frame="icrs").galactic
        l, b = g.l.deg, g.b.deg
    if args.parallax:
        d = 1000.0 / d
    ext = Extinction()
    vals = ext.bands(l, b, d)
    lo_err, hi_err = (e[0] for e in ext.zgr23_err(l, b, d))
    print(f"l = {l:.3f} deg, b = {b:.3f} deg, d = {d:.1f} pc")
    if d > C.R_OUTER_MAIN:
        print("  (beyond the map's 1.25 kpc: this is the map's whole column, a lower limit)")
    print("  (errors: voxels independent ... fully correlated)")
    for key, v in vals.items():
        val = float(np.ravel(v)[0])
        scale = val / max(float(np.ravel(vals["E_zgr23"])[0]), 1e-12)
        print(f"  {key:20s} {val:.3f} +- {scale * lo_err:.3f} ... {scale * hi_err:.3f} mag")


if __name__ == "__main__":
    main()
