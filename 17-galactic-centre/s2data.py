"""Read the two teams' S2 measurements into one common form.

Conventions after loading, for both teams:

* times are decimal Julian years (J2000.0 = 2000.0);
* ``x`` is the offset in right ascension, positive to the **East**, and ``y``
  the offset in declination, positive North, both in milliarcseconds from
  the origin of that team's reference frame (their estimate of where the
  radio source Sgr A* sits);
* radial velocities are in km/s, corrected to the Local Standard of Rest.

A note on signs, because it cost us an afternoon: in the UCLA release the
``x`` column is East-positive (the position angles in the same file confirm
it), while the printed Table S5 of Do et al. (2019) and their fitting code use
West-positive offsets. Their fitted frame offsets ``x0`` and ``vx0`` therefore
change sign when compared with ours; the orbital elements do not.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field

import numpy as np

from config import KECK_EXCLUDED_EPOCHS, RAW_DIR
from orbit_physics import mjd_to_year

DO2019_DIR = RAW_DIR / "do2019"


@dataclass
class Astrometry:
    """Sky positions of the star from one team, in one reference frame."""

    name: str
    t: np.ndarray
    x: np.ndarray
    y: np.ndarray
    x_err: np.ndarray
    y_err: np.ndarray
    instrument: np.ndarray
    #: Points the team itself excluded (kept for plotting, never fitted).
    excluded: np.ndarray = field(default=None)

    def __post_init__(self):
        if self.excluded is None:
            self.excluded = np.zeros(self.t.size, dtype=bool)

    def used(self) -> "Astrometry":
        keep = ~self.excluded
        return Astrometry(self.name, self.t[keep], self.x[keep], self.y[keep],
                          self.x_err[keep], self.y_err[keep], self.instrument[keep])

    @property
    def n(self) -> int:
        return int(self.t.size)


@dataclass
class RadialVelocities:
    name: str
    t: np.ndarray
    v: np.ndarray
    v_err: np.ndarray
    instrument: np.ndarray

    @property
    def n(self) -> int:
        return int(self.t.size)

    def subset(self, mask) -> "RadialVelocities":
        return RadialVelocities(self.name, self.t[mask], self.v[mask], self.v_err[mask],
                                self.instrument[mask])


# --------------------------------------------------------------------------
# UCLA / Keck (Do et al. 2019)
# --------------------------------------------------------------------------


def load_keck_astrometry() -> Astrometry:
    """46 Keck positions (27 speckle 1995-2005, then adaptive optics), in mas."""
    rows = list(csv.DictReader(open(DO2019_DIR / "S0-2_astro.csv")))
    t = mjd_to_year([float(r["mjd"]) for r in rows])
    x = np.array([float(r["x"]) for r in rows]) * 1e3
    y = np.array([float(r["y"]) for r in rows]) * 1e3
    xe = np.array([float(r["xerr"]) for r in rows]) * 1e3
    ye = np.array([float(r["yerr"]) for r in rows]) * 1e3
    # Speckle imaging on Keck I ran to mid-2005, then adaptive optics on Keck II
    # (Do et al. 2019, S1.2). This cut reproduces their count of 27 speckle epochs:
    # 19 in the file plus the 8 they excluded (1998 x5, 2002 x3).
    inst = np.where(t < 2005.5, "Keck speckle", "Keck AO")
    excluded = np.zeros(t.size, dtype=bool)
    for epoch in KECK_EXCLUDED_EPOCHS:
        excluded |= np.abs(t - epoch) < 0.002
    return Astrometry("Keck (UCLA)", t, x, y, xe, ye, inst, excluded)


def load_keck_rv() -> RadialVelocities:
    """116 LSR radial velocities: Keck, Subaru, Gemini, and the VLT's own values.

    The VLT (SINFONI) velocities in this file are the MPE team's published ones,
    so a joint fit must not add Gillessen et al.'s copies a second time.
    """
    rows = list(csv.DictReader(open(DO2019_DIR / "S0-2_rv.csv")))
    t = mjd_to_year([float(r["mjd"]) for r in rows])
    v = np.array([float(r["vlsr"]) for r in rows])
    ve = np.array([float(r["vzerr"]) for r in rows])
    inst = np.array([r["source"].strip() for r in rows])
    return RadialVelocities("Do et al. 2019", t, v, ve, inst)


def load_do2019_chain(which: str) -> dict:
    """The authors' posterior samples: 'Newton', 'GR' or 'redshift'.

    MultiNest output: column 0 is each sample's weight. GM is in units of
    10^6 GM_sun. Their x0 and vx0 are West-positive (see the module note); we
    flip them here so they compare directly with our East-positive fits.
    """
    names = ["weight", "gm6", "r0", "f", "x0", "y0", "vx0", "vy0", "vz0", "period", "t0",
             "e", "inc", "omega", "Omega", "nirc2_offset", "corr_length", "corr_mix"]
    raw = np.loadtxt(DO2019_DIR / f"chains_{which}.txt")
    chain = {n: raw[:, k] for k, n in enumerate(names)}
    chain["weight"] = chain["weight"] / chain["weight"].sum()
    chain["x0"] = -chain["x0"]
    chain["vx0"] = -chain["vx0"]
    return chain


# --------------------------------------------------------------------------
# MPE / VLT (Gillessen et al. 2017)
# --------------------------------------------------------------------------


def _g17_s2_rows():
    """S2's columns from the fixed-width CDS table (bytes 47-75, 1-based)."""
    for line in open(RAW_DIR / "gillessen2017_table5.dat"):
        line = line.rstrip("\n").ljust(557)
        chunk = line[46:75]
        if not chunk.strip():
            continue
        yield {
            "t": float(line[0:8]),
            "flag": line[9:11].strip(),
            "x": line[46:51], "xe": line[52:55], "y": line[56:61], "ye": line[62:65],
            "rv": line[66:71], "rve": line[72:75],
            "inst": line[551:557].strip(),
        }


def load_vlt_astrometry() -> Astrometry:
    """145 positions: NTT speckle (1992-2002), one Gemini epoch, VLT/NACO AO (2002-2016)."""
    rows = [r for r in _g17_s2_rows() if r["flag"] == "a"]
    t = np.array([r["t"] for r in rows])
    x = np.array([float(r["x"]) for r in rows])
    y = np.array([float(r["y"]) for r in rows])
    xe = np.array([float(r["xe"]) for r in rows])
    ye = np.array([float(r["ye"]) for r in rows])
    inst = np.array([{"VLT": "VLT NACO", "NTT": "NTT speckle", "GEMINI": "Gemini"}[r["inst"]]
                     for r in rows])
    return Astrometry("VLT (MPE)", t, x, y, xe, ye, inst)


def load_vlt_rv() -> RadialVelocities:
    """44 LSR velocities: 41 from VLT/SINFONI and 3 early Keck ones (2000, 2002)."""
    rows = [r for r in _g17_s2_rows() if r["flag"] == "rv"]
    t = np.array([r["t"] for r in rows])
    v = np.array([float(r["rv"]) for r in rows])
    ve = np.array([float(r["rve"]) for r in rows])
    inst = np.array(["VLT" if r["inst"] == "VLT" else "Keck (G17)" for r in rows])
    return RadialVelocities("Gillessen et al. 2017", t, v, ve, inst)


def load_s_star_elements() -> list[dict]:
    """Gillessen et al. 2017 Table 3: orbital elements of 40 S-stars."""
    stars = []
    for line in open(RAW_DIR / "gillessen2017_table3.dat"):
        line = line.rstrip("\n").ljust(118)
        per = line[90:97].strip()
        stars.append({
            "star": line[0:4].strip(),
            "a_arcsec": float(line[7:15]), "a_err": float(line[16:22]),
            "e": float(line[23:29]), "e_err": float(line[30:36]),
            "inc": float(line[37:43]), "Omega": float(line[49:55]), "omega": float(line[62:68]),
            "tp": float(line[75:82]),
            "period": float(per) if per else np.nan,
            "period_err": float(line[98:105]) if line[98:105].strip() else np.nan,
            "sptype": line[106:107].strip(),
            "kmag": float(line[108:113]),
        })
    return stars


def rv_overlap_report(keck_rv: RadialVelocities, vlt_rv: RadialVelocities) -> dict:
    """How Gillessen et al.'s velocities relate to the UCLA file.

    Their 41 VLT/SINFONI velocities should appear in the UCLA file unchanged;
    their three early Keck velocities (2000, 2002) appear there re-reduced.
    """
    vlt = vlt_rv.instrument == "VLT"
    identical = 0
    for t, v, ve in zip(vlt_rv.t[vlt], vlt_rv.v[vlt], vlt_rv.v_err[vlt]):
        d = np.abs(keck_rv.t - t)
        k = int(np.argmin(d))
        if d[k] < 0.004 and keck_rv.instrument[k] == "VLT" and abs(keck_rv.v[k] - v) < 0.5 \
                and abs(keck_rv.v_err[k] - ve) < 0.5:
            identical += 1
    changes = []
    for t, v in zip(vlt_rv.t[~vlt], vlt_rv.v[~vlt]):
        k = int(np.argmin(np.abs(keck_rv.t - t)))
        changes.append(float(keck_rv.v[k] - v))
    return {"n_vlt_velocities": int(vlt.sum()), "n_identical_in_ucla_file": identical,
            "n_early_keck_velocities": int((~vlt).sum()), "early_keck_revisions_kms": changes}
