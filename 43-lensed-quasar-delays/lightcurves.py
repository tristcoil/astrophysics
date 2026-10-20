"""Read the COSMOGRAIL light curves into one common structure.

Every data set is a table of epochs. At each epoch all images of the lensed
quasar were measured on the same frames by simultaneous deconvolution
photometry (Magain, Courbin & Sohy 1998), so the images share their dates.
Magnitudes are *relative* (arbitrary zero point per image and telescope):
only changes matter for time delays.

Time stamps are heliocentric Julian dates minus 2,400,000 (the CDS tables use
HJD - 2400000 or heliocentric MJD; the half-day difference between the two
conventions is irrelevant because every delay is measured within one data set).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

import config as C


@dataclass
class LightCurves:
    """Light curves of the images of one lensed quasar."""

    key: str  # config.LENSES key, or the COSMOGRAIL file stem for the gallery
    name: str
    images: list[str]
    t: np.ndarray  # (N,) days, HJD - 2400000
    mag: dict[str, np.ndarray]  # image -> (N,) relative magnitude
    err: dict[str, np.ndarray]  # image -> (N,) 1-sigma error
    tel: np.ndarray  # (N,) telescope / instrument label
    source: str = ""
    notes: list[str] = field(default_factory=list)

    def __post_init__(self):
        order = np.argsort(self.t, kind="stable")
        self.t = np.asarray(self.t, float)[order]
        self.tel = np.asarray(self.tel).astype(str)[order]
        self.mag = {k: np.asarray(v, float)[order] for k, v in self.mag.items()}
        self.err = {k: np.asarray(v, float)[order] for k, v in self.err.items()}

    @property
    def n(self) -> int:
        return len(self.t)

    @property
    def telescopes(self) -> list[str]:
        """Telescopes in order of first appearance."""
        _, idx = np.unique(self.tel, return_index=True)
        return [self.tel[i] for i in sorted(idx)]

    def select(self, mask) -> "LightCurves":
        mask = np.asarray(mask, bool)
        return LightCurves(self.key, self.name, list(self.images), self.t[mask],
                           {k: v[mask] for k, v in self.mag.items()},
                           {k: v[mask] for k, v in self.err.items()},
                           self.tel[mask], self.source, list(self.notes))

    def seasons(self, gap: float = 60.0) -> np.ndarray:
        """Season index of each epoch: a new season starts after a gap > ``gap`` days."""
        return np.concatenate([[0], np.cumsum(np.diff(self.t) > gap)])


def decimal_year(t) -> np.ndarray:
    """HJD - 2400000 -> calendar year (J2000.0 = JD 2451545.0)."""
    return 2000.0 + (np.asarray(t, float) + 2_400_000.0 - 2_451_545.0) / 365.25


# ---------------------------------------------------------------------------
# File readers
# ---------------------------------------------------------------------------


def _read_whitespace(path: Path, labels: list[str], has_tel: bool = True):
    """COSMOGRAIL XIX layout: time, (mag, err) per image, telescope."""
    t, mags, errs, tels = [], {k: [] for k in labels}, {k: [] for k in labels}, []
    for line in path.read_text().splitlines():
        tok = line.split()
        if not tok:
            continue
        t.append(float(tok[0]))
        vals = tok[1:-1] if has_tel else tok[1:]
        if len(vals) != 2 * len(labels):
            raise ValueError(f"{path.name}: expected {2 * len(labels)} numbers, got {len(vals)}: {line!r}")
        for i, k in enumerate(labels):
            mags[k].append(float(vals[2 * i]))
            errs[k].append(float(vals[2 * i + 1]))
        tels.append(tok[-1] if has_tel else path.stem)
    return np.array(t), {k: np.array(v) for k, v in mags.items()}, {k: np.array(v) for k, v in errs.items()}, np.array(tels)


def _read_cds(folder: str, filename: str):
    """A CDS table with its ReadMe (astropy matches the two by file name)."""
    from astropy.io import ascii

    d = C.DATA / folder
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # unit-parsing chatter about 'mag'
        return ascii.read(str(d / filename), format="cds", readme=str(d / "ReadMe"))


def _flux_sum_mag(m1, e1, m2, e2):
    """Magnitude and error of the summed flux of two images."""
    f1, f2 = 10 ** (-0.4 * m1), 10 ** (-0.4 * m2)
    f = f1 + f2
    ef = np.hypot(f1 * e1, f2 * e2) * 0.4 * np.log(10)
    return -2.5 * np.log10(f), ef / f / (0.4 * np.log(10))


#: COSMOGRAIL XIX sub-folder -> image labels, in column order.
CG19_LABELS = {"lcab": ["A", "B"], "lcabcd": ["A", "B", "C", "D"], "lcapd": ["A+D", "B", "C"], "lcapb": ["A+B", "C", "D"]}


def read_cosmograil19(relpath: str, key: str | None = None, name: str | None = None) -> LightCurves:
    sub, fname = relpath.split("/")
    labels = list(CG19_LABELS[sub])
    if fname.startswith("WFI2026"):
        labels = ["A", "B", "C"]  # no D column (ReadMe note 1)
    t, mag, err, tel = _read_whitespace(C.DATA / "cosmograil19" / sub / fname, labels)
    return LightCurves(key or Path(fname).stem, name or Path(fname).stem, labels, t, mag, err, tel,
                       source="Millon et al. 2020a (COSMOGRAIL XIX), VizieR J/A+A/640/A105")


def read_lens(key: str, which: str = "lc") -> LightCurves:
    """Light curves of one of the lenses in ``config.LENSES``.

    ``which="lc"`` gives the data set we measure delays on; ``"lc_h0licow"``
    (RX J1131, HE 0435 only) gives the shorter data set behind the published
    lens model's delays.
    """
    L = C.LENSES[key]
    folder, files = L[which]
    name = L["name"]

    if folder == "cosmograil19":
        lc = read_cosmograil19(files, key, name)
        # Image labels in the file match the lens's own labels.
        assert lc.images == L["images"], (lc.images, L["images"])
        return lc

    if key == "RXJ1131" and folder == "tewes2013_rxj1131":
        tab = _read_cds(folder, files)
        ims = L["images"]
        return LightCurves(key, name, ims, np.asarray(tab["HJD"]) - 2_400_000.0,
                           {k: np.asarray(tab[f"{k}mag"]) for k in ims}, {k: np.asarray(tab[f"e_{k}mag"]) for k in ims},
                           np.asarray(tab["Tel"]), source="Tewes et al. 2013b, VizieR J/A+A/556/A22")

    if key == "HE0435" and folder == "bonvin2017_he0435":
        tab = _read_cds(folder, files)
        ims = L["images"]
        return LightCurves(key, name, ims, np.asarray(tab["MHJD"]),
                           {k: np.asarray(tab[f"mag{k}"]) for k in ims}, {k: np.asarray(tab[f"e_mag{k}"]) for k in ims},
                           np.full(len(tab), "all"), source="Bonvin et al. 2017, VizieR J/MNRAS/465/4914")

    if key == "PG1115":
        tab = _read_cds(folder, files)
        ims = L["images"]
        return LightCurves(key, name, ims, np.asarray(tab["MHJD"]),
                           {k: np.asarray(tab[f"mag{k}"]) for k in ims}, {k: np.asarray(tab[f"e_mag{k}"]) for k in ims},
                           np.full(len(tab), "WFI"), source="Bonvin et al. 2018, VizieR J/A+A/616/A183",
                           notes=["A is the unresolved merging pair A1 + A2"])

    if key == "WFI2033":
        # Euler (C2, ECAM) and SMARTS measure A = A1 + A2 together; the 2.2 m WFI
        # resolves A1 and A2, so its two fluxes are summed to make A as well.
        t, tel = [], []
        mag = {k: [] for k in ("A", "B", "C")}
        err = {k: [] for k in ("A", "B", "C")}
        for fname in files:
            tab = _read_cds(folder, fname)
            if fname == "wfi.dat":
                mA, eA = _flux_sum_mag(np.asarray(tab["magA1"]), np.asarray(tab["e_magA1"]),
                                       np.asarray(tab["magA2"]), np.asarray(tab["e_magA2"]))
            else:
                mA, eA = np.asarray(tab["magA"]), np.asarray(tab["e_magA"])
            t.append(np.asarray(tab["MHJD"]))
            mag["A"].append(mA)
            err["A"].append(eA)
            for k in ("B", "C"):
                mag[k].append(np.asarray(tab[f"mag{k}"]))
                err[k].append(np.asarray(tab[f"e_mag{k}"]))
            tel.append(np.full(len(tab), {"ecam.dat": "ECAM", "c2.dat": "EulerC2", "smarts.dat": "SMARTS", "wfi.dat": "WFI"}[fname]))
        cat = np.concatenate
        return LightCurves(key, name, ["A", "B", "C"], cat(t), {k: cat(v) for k, v in mag.items()},
                           {k: cat(v) for k, v in err.items()}, cat(tel),
                           source="Bonvin et al. 2019, VizieR J/A+A/629/A97",
                           notes=["A is the merging pair A1 + A2 (summed flux)"])

    if key == "J1206":
        tab = _read_cds(folder, files)
        return LightCurves(key, name, ["A", "B"], np.asarray(tab["MJD"]),
                           {k: np.asarray(tab[f"R{k}mag"]) for k in "AB"}, {k: np.asarray(tab[f"e_R{k}mag"]) for k in "AB"},
                           np.asarray(tab["Tel"]), source="Eulaers et al. 2013, VizieR J/A+A/553/A121")

    if key == "DES0408":
        tab = _read_cds(folder, files)
        ims = L["images"]
        return LightCurves(key, name, ims, np.asarray(tab["MHJD"]),
                           {k: np.asarray(tab[f"mag{k}"]) for k in ims}, {k: np.asarray(tab[f"e_mag{k}"]) for k in ims},
                           np.asarray(tab["Tel"]), source="Courbin et al. 2018, VizieR J/A+A/609/A71")

    raise KeyError(f"no light curves for {key} ({which})")


def read_gallery() -> list[LightCurves]:
    """All 23 COSMOGRAIL XIX systems, in the order of the CDS star list."""
    out = []
    for line in (C.DATA / "cosmograil19" / "stars.dat").read_text().splitlines():
        name, relpath = line[0:18].strip(), line[43:67].strip()
        if relpath:
            out.append(read_cosmograil19(relpath, key=Path(relpath).stem, name=name))
    return out


if __name__ == "__main__":
    for key in C.MEASURED:
        lc = read_lens(key)
        span = decimal_year([lc.t.min(), lc.t.max()])
        print(f"{key:8s} {lc.n:5d} epochs  {span[0]:.1f}-{span[1]:.1f}  images {lc.images}  "
              f"telescopes {dict(zip(*np.unique(lc.tel, return_counts=True)))}")
    for key in ("RXJ1131", "HE0435"):
        lc = read_lens(key, "lc_h0licow")
        span = decimal_year([lc.t.min(), lc.t.max()])
        print(f"{key:8s} (H0LiCOW data) {lc.n:5d} epochs  {span[0]:.1f}-{span[1]:.1f}")
    print(len(read_gallery()), "COSMOGRAIL XIX systems")
