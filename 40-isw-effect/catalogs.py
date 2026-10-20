"""The superstructure catalogues, read into one common format.

Three independent catalogues, all found with the same family of algorithm
(ZOBOV, Neyrinck 2008), which grows voids around density minima of the Voronoi
tessellation of the galaxies (and superclusters around maxima, by inverting the
density):

* **Granett et al. 2008** -- 50 supervoids and 50 superclusters in ~1.1 million
  photometric SDSS luminous red galaxies (0.4 < z < 0.75). The catalogue behind
  the famous 4-sigma claim. Positions, redshifts and sky radii; volumes in
  (h^-1 Mpc)^3.
* **Nadathur & Hotchkiss 2014** -- voids and superclusters in the SDSS DR7
  spectroscopic samples: luminous red galaxies (lrgdim 0.16 < z < 0.36,
  lrgbright 0.36 < z < 0.44) and four volume-limited main-galaxy samples
  (z < 0.2). Type1 voids must start from a zone below 0.3 of the mean density.
* **Mao et al. 2017** -- 1,228 voids in BOSS DR12 (LOWZ and CMASS, 0.2 < z <
  0.7) after quality cuts. Effective radius = radius of a sphere of the same
  volume.

Every structure gets an angular radius ``theta`` (deg): the effective radius
divided by the comoving distance. Filters scale with it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

import cosmology
from config import product


@dataclass
class Catalog:
    name: str
    kind: str  # "void" or "cluster"
    ra: np.ndarray
    dec: np.ndarray
    z: np.ndarray
    r_eff: np.ndarray  # h^-1 Mpc
    theta: np.ndarray  # deg, angular effective radius
    extra: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.ra)

    def subset(self, sel, name: str | None = None) -> "Catalog":
        sel = np.asarray(sel)
        return Catalog(name or self.name, self.kind, self.ra[sel], self.dec[sel], self.z[sel],
                       self.r_eff[sel], self.theta[sel],
                       {k: np.asarray(v)[sel] for k, v in self.extra.items()})

    @property
    def sign(self) -> int:
        """Expected sign of the ISW imprint: voids cold, clusters hot."""
        return -1 if self.kind == "void" else 1


def _theta_from_r(r_eff, z):
    return np.degrees(r_eff / cosmology.LCDM.chi(z))


# --------------------------------------------------------------------------


def granett(kind: str) -> Catalog:
    """Granett, Neyrinck & Szapudi (2008), as archived on their (now offline) page.

    Columns: z ra dec mean_r_sky max_r_sky vol d_all d_neg|d_pos d_min|d_max
    denscontrast prob. ``mean_r_sky`` is the mean angular distance (deg) of
    member galaxies from the centre: our ``theta``. Rows are ordered by density
    contrast, most significant first.
    """
    key = "granett_voids" if kind == "void" else "granett_clusters"
    a = np.loadtxt(product(key).path, skiprows=1)
    vol = a[:, 5]
    r_eff = (3 * vol / (4 * np.pi)) ** (1 / 3)
    return Catalog(f"Granett+08 {'supervoids' if kind == 'void' else 'superclusters'}", kind,
                   a[:, 1], a[:, 2], a[:, 0], r_eff, a[:, 3],
                   {"max_r_sky": a[:, 4], "denscontrast": a[:, 9], "prob": a[:, 10],
                    "rank": np.arange(len(a))})


def mao() -> Catalog:
    """Mao et al. (2017) BOSS DR12 voids, fixed-width columns per the CDS ReadMe."""
    rows = []
    with open(product("mao_voids").path) as fh:
        for line in fh:
            if not line.strip():
                continue
            rows.append((line[0:11].strip(), int(line[12:17]), float(line[18:25]),
                         float(line[26:32]), float(line[33:38]), int(line[39:45]),
                         float(line[46:55]), float(line[56:63]), float(line[64:73]),
                         float(line[74:80]), float(line[81:86]), float(line[87:96]),
                         float(line[97:104])))
    sample = np.array([r[0] for r in rows])
    num = np.array([r[1:] for r in rows], float)
    ra, dec, z, ngal, vol, reff, nmin, delmin, ratio, prob, dbound = num[:, 1:].T
    return Catalog("Mao+17 BOSS voids", "void", ra, dec, z, reff, _theta_from_r(reff, z),
                   {"sample": sample, "ngal": ngal, "delmin": delmin, "ratio": ratio,
                    "prob": prob, "dbound": dbound})


NH14_SAMPLES = {
    "lrgdim": ("lrgdt1v", "lrgdt1c"),
    "lrgbright": ("lrgbt1v", "lrgbt1c"),
    "dim1": ("dim1t1v", "dim1t1c"),
    "dim2": ("dim2t1v", "dim2t1c"),
    "bright1": ("bri1t1v", "bri1t1c"),
    "bright2": ("bri2t1v", "bri2t1c"),
}


def nh14(sample: str, kind: str) -> Catalog:
    """Nadathur & Hotchkiss (2014) Type1 voids or superclusters of one DR7 sample.

    Columns (CDS ReadMe): zone, RA, Dec, z, R_eff (h^-1 Mpc), theta_eff (deg),
    average density, min (voids) or max (clusters) density, edge flag, density
    ratio. ``theta_eff`` is used as published, as Hotchkiss et al. (2015) did.
    """
    stem = NH14_SAMPLES[sample][0 if kind == "void" else 1]
    a = np.atleast_2d(np.loadtxt(product(f"nh14_{stem}").path))
    label = "Type1 voids" if kind == "void" else "superclusters"
    return Catalog(f"NH14 {sample} {label}", kind, a[:, 1], a[:, 2], a[:, 3], a[:, 4], a[:, 5],
                   {"zone": a[:, 0], "avg_dens": a[:, 6], "ext_dens": a[:, 7],
                    "edge": a[:, 8].astype(int), "ratio": a[:, 9]})


def galactic(cat: Catalog):
    """Galactic (l, b) in degrees."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    c = SkyCoord(ra=cat.ra * u.deg, dec=cat.dec * u.deg, frame="icrs").galactic
    return c.l.deg, c.b.deg


def all_catalogs() -> dict[str, Catalog]:
    cats = {"granett_voids": granett("void"), "granett_clusters": granett("cluster"),
            "mao_voids": mao()}
    for s in ("lrgdim", "lrgbright"):
        cats[f"nh14_{s}_voids"] = nh14(s, "void")
        cats[f"nh14_{s}_clusters"] = nh14(s, "cluster")
    return cats
