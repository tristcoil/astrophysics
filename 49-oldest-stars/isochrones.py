"""Theoretical isochrones from three independent stellar-evolution codes, on one footing.

An isochrone is the locus of stars of one age and composition in the
colour-magnitude diagram. Each model set tabulates them on its own grid of age,
iron abundance [Fe/H] and alpha-element enhancement [alpha/Fe]. To fit a cluster
we need isochrones *between* the grid points, so every set is read into the same
structure and interpolated the same way:

* along the isochrone, points are matched by evolutionary phase, not by mass
  ("equivalent evolutionary points", EEPs: MIST and DSED label them, BaSTI's
  2,100 points per isochrone follow fixed phases, so its row number is one);
* between grid points, each quantity is interpolated at fixed EEP: linearly in
  log(age) and [alpha/Fe], and in [Fe/H] linearly (MIST, BaSTI: 0.25 dex steps
  or finer) or quadratically (DSED: 0.5 dex steps; see test_physics.py for the
  size of the interpolation error).

Magnitudes here are absolute and unreddened. ``Extinction`` turns the MIST
bolometric-correction tables into the dimming of each model star by dust,
computed with that star's own spectrum (hot turnoff stars and cool giants are
dimmed by slightly different amounts in the same filter).

Built grids are cached in ``data/derived/isogrid_<set>.npz``; building MIST the
first time streams the 1.1 GB tarball (8 GB uncompressed) once, ~3 minutes.
"""

from __future__ import annotations

import io
import re
import tarfile
from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import RegularGridInterpolator

import config as C

QUANTS = ("mass", "logTe", "logg", "logL", "F606W", "F814W")
#: Ages kept from each grid (Gyr): wide enough for every cluster and for the
#: interpolation brackets at the ends of the fitted age range.
KEEP_AGES = (6.5, 18.0)
#: Faintest absolute F606W kept on the main sequence (the fits use stars at
#: most ~3 mag below the turnoff, M_F606W < 8).
FAINT_LIMIT = 9.0

LOGG_SUN = 4.438  # IAU 2015 nominal solar values
LOGTE_SUN = np.log10(5772.0)


@dataclass
class IsoGrid:
    """Isochrones on a regular (feh, alpha, log_age, eep) grid, NaN where absent."""

    name: str
    feh: np.ndarray
    alpha: np.ndarray
    log_age: np.ndarray
    eep: np.ndarray
    data: dict = field(repr=False)
    feh_order: int = 1  # 1 = linear, 2 = quadratic interpolation in [Fe/H]
    #: EEP (or row) of the main-sequence turnoff region and of the red-giant-branch
    #: tip: the fits use points between the lower main sequence and the tip.
    eep_tip: int = 0

    # -- persistence ----------------------------------------------------------
    def save(self, path) -> None:
        np.savez_compressed(
            path, name=self.name, feh=self.feh, alpha=self.alpha, log_age=self.log_age,
            eep=self.eep, feh_order=self.feh_order, eep_tip=self.eep_tip,
            **{f"q_{k}": v for k, v in self.data.items()})

    @classmethod
    def load(cls, path) -> "IsoGrid":
        z = np.load(path)
        return cls(name=str(z["name"]), feh=z["feh"], alpha=z["alpha"], log_age=z["log_age"],
                   eep=z["eep"], data={k: z[f"q_{k}"] for k in QUANTS},
                   feh_order=int(z["feh_order"]), eep_tip=int(z["eep_tip"]))

    # -- interpolation -------------------------------------------------------
    def _feh_weights(self, feh: float) -> tuple[np.ndarray, np.ndarray]:
        nodes = self.feh
        if self.feh_order == 2 and len(nodes) >= 3:
            # three consecutive nodes centred on the one nearest to feh (Lagrange)
            k = int(np.clip(np.argmin(abs(nodes - feh)), 1, len(nodes) - 2))
            idx = np.array([k - 1, k, k + 1])
            x = nodes[idx]
            w = np.array([
                (feh - x[1]) * (feh - x[2]) / ((x[0] - x[1]) * (x[0] - x[2])),
                (feh - x[0]) * (feh - x[2]) / ((x[1] - x[0]) * (x[1] - x[2])),
                (feh - x[0]) * (feh - x[1]) / ((x[2] - x[0]) * (x[2] - x[1])),
            ])
            return idx, w
        return _linear_weights(nodes, feh)

    def at_composition(self, feh: float, alpha: float) -> dict:
        """All native ages at one ([Fe/H], [alpha/Fe]): arrays (n_age, n_eep)."""
        fi, fw = self._feh_weights(feh)
        ai, aw = _linear_weights(self.alpha, alpha)
        out = {}
        for q, arr in self.data.items():
            acc = None
            for i, wi in zip(fi, fw):
                for j, wj in zip(ai, aw):
                    term = (wi * wj) * arr[i, j].astype(np.float64)
                    acc = term if acc is None else acc + term
            out[q] = acc
        return out

    def interpolate(self, feh: float, alpha: float, ages_gyr) -> dict:
        """Isochrones at the requested ages (Gyr): arrays (n_ages, n_eep)."""
        comp = self.at_composition(feh, alpha)
        return interpolate_in_age(comp, self.log_age, ages_gyr)


def _linear_weights(nodes: np.ndarray, x: float) -> tuple[np.ndarray, np.ndarray]:
    if len(nodes) == 1:
        return np.array([0]), np.array([1.0])
    j = int(np.clip(np.searchsorted(nodes, x), 1, len(nodes) - 1))
    t = (x - nodes[j - 1]) / (nodes[j] - nodes[j - 1])  # extrapolates linearly off the grid
    return np.array([j - 1, j]), np.array([1.0 - t, t])


def interpolate_in_age(comp: dict, log_age_nodes: np.ndarray, ages_gyr) -> dict:
    """Linear interpolation in log(age) at fixed EEP."""
    la = np.log10(np.atleast_1d(np.asarray(ages_gyr, dtype=float)) * 1e9)
    j = np.clip(np.searchsorted(log_age_nodes, la), 1, len(log_age_nodes) - 1)
    t = ((la - log_age_nodes[j - 1]) / (log_age_nodes[j] - log_age_nodes[j - 1]))[:, None]
    return {q: (1.0 - t) * v[j - 1] + t * v[j] for q, v in comp.items()}


def logg_from(mass, logL, logTe):
    """Surface gravity from mass, luminosity and temperature (solar units)."""
    return LOGG_SUN + np.log10(mass) + 4.0 * (logTe - LOGTE_SUN) - logL


# ---------------------------------------------------------------------------
# MIST v2.5
# ---------------------------------------------------------------------------

_MIST_NAME = re.compile(r"feh_([mp])(\d{3})_afe_([mp])(\d)_(vvcrit[\d.]+)_full\.iso\.HST_ACS_WFC$")
MIST_EEP_TIP = 605  # red-giant-branch tip


def _parse_mist_member(text: str) -> dict:
    """Parse one MIST .iso file, keeping the ages in KEEP_AGES."""
    lines = text.splitlines()
    out = {}
    i, n = 0, len(lines)
    lo, hi = np.log10(np.array(KEEP_AGES) * 1e9)
    while i < n:
        line = lines[i]
        if line.startswith("# number of EEPs"):
            neep = int(line.split("=")[1].split()[0])
            header = lines[i + 2].lstrip("#").split()
            first = lines[i + 3].split()
            log_age = float(first[1])
            if lo - 0.06 <= log_age <= hi + 0.06:
                block = np.array([ln.split() for ln in lines[i + 3:i + 3 + neep]], dtype=float)
                col = {name: k for k, name in enumerate(header)}
                out[round(log_age, 4)] = {
                    "eep": block[:, col["EEP"]].astype(int),
                    "mass": block[:, col["initial_mass"]],
                    "logTe": block[:, col["log_Teff"]],
                    "logg": block[:, col["log_g"]],
                    "logL": block[:, col["log_L"]],
                    "F606W": block[:, col["ACS_WFC_F606W"]],
                    "F814W": block[:, col["ACS_WFC_F814W"]],
                }
            i += 3 + neep
        else:
            i += 1
    return out


def build_mist(feh_max: float = 0.0, feh_min: float = -3.0) -> IsoGrid:
    path = C.MIST_DIR / "HST_ACS_WFC.txz"
    blocks = {}
    print(f"  building the MIST grid from {path.name} (streams 8 GB, a few minutes)")
    with tarfile.open(path, mode="r|xz") as tar:
        for member in tar:
            m = _MIST_NAME.search(member.name)
            if not m or m.group(5) != C.MIST_ROTATION:
                continue
            feh = (-1 if m.group(1) == "m" else 1) * int(m.group(2)) / 100
            afe = (-1 if m.group(3) == "m" else 1) * int(m.group(4)) / 10
            if not (feh_min - 1e-6 <= feh <= feh_max + 1e-6):
                continue
            blocks[(feh, afe)] = _parse_mist_member(tar.extractfile(member).read().decode())
            print(f"    [Fe/H] = {feh:+.2f}  [alpha/Fe] = {afe:+.1f}: {len(blocks[(feh, afe)])} ages")
    return _assemble(blocks, name="MIST", eep_tip=MIST_EEP_TIP, feh_order=FEH_ORDER["MIST"])


def _assemble(blocks: dict, *, name: str, eep_tip: int, feh_order: int) -> IsoGrid:
    fehs = np.array(sorted({k[0] for k in blocks}))
    alphas = np.array(sorted({k[1] for k in blocks}))
    ages = np.array(sorted({a for v in blocks.values() for a in v}))
    # EEPs on the kept range: from the faintest point brighter than FAINT_LIMIT
    # in *any* isochrone, to the RGB tip.
    eep_min = min(int(b["eep"][np.argmax(b["F606W"] < FAINT_LIMIT)])
                  for v in blocks.values() for b in v.values())
    eeps = np.arange(eep_min, eep_tip + 1)
    data = {q: np.full((len(fehs), len(alphas), len(ages), len(eeps)), np.nan, np.float32) for q in QUANTS}
    for (feh, afe), v in blocks.items():
        i, j = int(np.searchsorted(fehs, feh)), int(np.searchsorted(alphas, afe))
        for la, b in v.items():
            k = int(np.searchsorted(ages, la))
            sel = (b["eep"] >= eep_min) & (b["eep"] <= eep_tip)
            pos = b["eep"][sel] - eep_min
            for q in QUANTS:
                data[q][i, j, k, pos] = b[q][sel]
    return IsoGrid(name=name, feh=fehs, alpha=alphas, log_age=ages, eep=eeps, data=data,
                   feh_order=feh_order, eep_tip=eep_tip)


# ---------------------------------------------------------------------------
# BaSTI-IAC
# ---------------------------------------------------------------------------

#: BaSTI rows kept: 150 (M_F606W ~ 8.5-9 on the lower main sequence) to the
#: last red-giant point before the helium flash (row 1290).
BASTI_ROWS = (150, 1290)
BASTI_ROW_TIP = 1290


def _parse_basti_tarball(path) -> dict:
    out = {}
    lo, hi = KEEP_AGES
    with tarfile.open(path) as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue
            age_myr = float(member.name.split("/")[-1].split("z")[0])
            if not (lo * 1e3 - 150 <= age_myr <= hi * 1e3 + 150):
                continue
            a = np.loadtxt(tar.extractfile(member), comments="#")
            r0, r1 = BASTI_ROWS
            a = a[r0:r1 + 1]
            mass, mfin, logL, logTe = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
            out[round(np.log10(age_myr * 1e6), 5)] = {
                "eep": np.arange(r0, r1 + 1),
                "mass": mass,
                "logTe": logTe,
                "logg": logg_from(mfin, logL, logTe),
                "logL": logL,
                "F606W": a[:, 7],
                "F814W": a[:, 10],
            }
    return out


def build_basti(grids=("P04O1D1E1Y247", "P00O1D1E1Y247"), name="BaSTI") -> IsoGrid:
    blocks = {}
    for grid in grids:
        afe = 0.4 if grid.startswith("P04") else 0.0
        for feh_label, node in C.BASTI_FEH_NODES.items():
            path = C.BASTI_DIR / grid / f"{node}.tar.gz"
            if not path.exists():
                continue
            # Nodes are labelled by their nominal [Fe/H]; the alpha-enhanced "-1.70"
            # node is really -1.71, a difference far below every error here.
            blocks[(feh_label, afe)] = _parse_basti_tarball(path)
            print(f"    BaSTI {grid} [Fe/H] = {feh_label:+.2f}: {len(blocks[(feh_label, afe)])} ages")
    return _assemble_rows(blocks, name=name)


def _assemble_rows(blocks: dict, *, name: str) -> IsoGrid:
    fehs = np.array(sorted({k[0] for k in blocks}))
    alphas = np.array(sorted({k[1] for k in blocks}))
    ages = np.array(sorted({a for v in blocks.values() for a in v}))
    r0, r1 = BASTI_ROWS
    eeps = np.arange(r0, r1 + 1)
    data = {q: np.full((len(fehs), len(alphas), len(ages), len(eeps)), np.nan, np.float32) for q in QUANTS}
    for (feh, afe), v in blocks.items():
        i, j = int(np.searchsorted(fehs, feh)), int(np.searchsorted(alphas, afe))
        for la, b in v.items():
            k = int(np.searchsorted(ages, la))
            for q in QUANTS:
                data[q][i, j, k] = b[q]
    return IsoGrid(name=name, feh=fehs, alpha=alphas, log_age=ages, eep=eeps, data=data,
                   feh_order=1, eep_tip=BASTI_ROW_TIP)


# ---------------------------------------------------------------------------
# Dartmouth DSED
# ---------------------------------------------------------------------------

_DSED_NAME = re.compile(r"feh([mp])(\d\d)afe([mp])(\d)\.HST_ACSWF$")
#: DSED isochrones stop at the RGB tip; their EEP labels run to ~280.
DSED_EEP_MAX = 290


def _parse_dsed_member(text: str) -> dict:
    lines = text.splitlines()
    out = {}
    i, n = 0, len(lines)
    lo, hi = KEEP_AGES
    while i < n:
        line = lines[i]
        if line.startswith("#AGE="):
            m = re.match(r"#AGE=\s*([\d.]+)\s+EEPS=\s*(\d+)", line)
            age, neep = float(m.group(1)), int(m.group(2))
            header = lines[i + 1].lstrip("#").split()
            if lo - 0.6 <= age <= hi + 0.6:
                block = np.array([ln.split() for ln in lines[i + 2:i + 2 + neep]], dtype=float)
                col = {name: k for k, name in enumerate(header)}
                out[round(np.log10(age * 1e9), 5)] = {
                    "eep": block[:, col["EEP"]].astype(int),
                    "mass": block[:, col["M/Mo"]],
                    "logTe": block[:, col["LogTeff"]],
                    "logg": block[:, col["LogG"]],
                    "logL": block[:, col["LogL/Lo"]],
                    "F606W": block[:, col["F606W"]],
                    "F814W": block[:, col["F814W"]],
                }
            i += 2 + neep
        else:
            i += 1
    return out


def build_dsed() -> IsoGrid:
    path = C.DSED_DIR / "HST_ACSWF.tgz"
    blocks = {}
    with tarfile.open(path) as tar:
        for member in tar.getmembers():
            m = _DSED_NAME.search(member.name)
            if not m:  # skips the helium-rich (y33, y40) and young (_2) files
                continue
            # Nominal composition from the file name (the headers quote [Fe/H]
            # recomputed from Z, which drifts by up to 0.02 dex with [alpha/Fe]).
            feh = (-1 if m.group(1) == "m" else 1) * int(m.group(2)) / 10
            afe = (-1 if m.group(3) == "m" else 1) * int(m.group(4)) / 10
            # [Fe/H] > 0 is not needed; [Fe/H] <= -3 exists only at [alpha/Fe] = +0.4
            # and only for 10-14 Gyr, so it would punch holes in the grid.
            if feh > 0.01 or feh < -2.51:
                continue
            blocks[(feh, afe)] = _parse_dsed_member(tar.extractfile(member).read().decode())
    # [Fe/H] < -2.5 exists only at [alpha/Fe] = +0.4: fill the other alpha values
    # by carrying over the alpha effect at the nearest complete node (-2.5), where
    # metals barely matter. Done at fixed EEP, before assembling.
    full = sorted({f for f, a in blocks if all((f, x) in blocks for x in (-0.2, 0.0, 0.2, 0.4, 0.6, 0.8))})
    alphas = sorted({a for _, a in blocks})
    for feh in sorted({f for f, _ in blocks}):
        if feh in full:
            continue
        ref = min(full, key=lambda f: abs(f - feh))
        for afe in alphas:
            if (feh, afe) in blocks:
                continue
            blocks[(feh, afe)] = _shift_alpha(blocks[(feh, 0.4)], blocks[(ref, afe)], blocks[(ref, 0.4)])
    return _assemble(blocks, name="DSED", eep_tip=DSED_EEP_MAX, feh_order=FEH_ORDER["DSED"])


def _shift_alpha(base: dict, ref_a: dict, ref_04: dict) -> dict:
    out = {}
    for la, b in base.items():
        if la not in ref_a or la not in ref_04:
            continue
        common = np.intersect1d(np.intersect1d(b["eep"], ref_a[la]["eep"]), ref_04[la]["eep"])
        new = {"eep": common}
        for q in QUANTS:
            pick = lambda d: d[q][np.searchsorted(d["eep"], common)]  # noqa: E731
            new[q] = pick(b) + (pick(ref_a[la]) - pick(ref_04[la]))
        out[la] = new
    return out


# ---------------------------------------------------------------------------
# Loading with a cache
# ---------------------------------------------------------------------------

#: Interpolation order in [Fe/H]. MIST's packaged ACS isochrones skip [Fe/H] =
#: -2.75 and -2.25 (the website lists them for the tracks only), so MIST and DSED
#: both have 0.5 dex steps where the metal-poor clusters sit: quadratic there.
FEH_ORDER = {"MIST": 2, "DSED": 2, "BaSTI": 1, "BaSTI_Y275": 1}

BUILDERS = {
    "MIST": build_mist,
    "BaSTI": build_basti,
    "DSED": build_dsed,
    "BaSTI_Y275": lambda: build_basti(grids=("P04O1D1E1Y275",), name="BaSTI_Y275"),
}


def load_grid(name: str) -> IsoGrid:
    path = C.DERIVED_DIR / f"isogrid_{name}.npz"
    if path.exists():
        grid = IsoGrid.load(path)
    else:
        grid = BUILDERS[name]()
        grid.save(path)
    grid.feh_order = FEH_ORDER[name]
    return grid


# ---------------------------------------------------------------------------
# Extinction from the MIST bolometric-correction tables
# ---------------------------------------------------------------------------


class Extinction:
    """A_X(Teff, log g, A_V) = BC_X(A_V = 0) - BC_X(A_V) at R_V = 3.1.

    Uses the table for the composition nearest the cluster's; the dependence of
    the *extinction* on composition is far below the other uncertainties.
    """

    def __init__(self, feh: float, alpha: float = 0.4, system: str = "HST_ACS_WFC"):
        path = C.MIST_DIR / f"BC_{system}.txz"
        fnodes = np.arange(-3.0, 0.51, 0.25)
        anodes = np.array([-0.2, 0.0, 0.2, 0.4, 0.6])
        f = fnodes[np.argmin(abs(fnodes - feh))]
        a = anodes[np.argmin(abs(anodes - alpha))]
        member = f"feh{f:+.2f}_afe{a:+.1f}.{system}"
        cache = C.DERIVED_DIR / f"bc_{member}.npz"
        if not cache.exists():
            with tarfile.open(path) as tar:
                text = tar.extractfile(member).read().decode()
            lines = text.splitlines()
            header = lines[3].lstrip("#").split()
            arr = np.loadtxt(io.StringIO("\n".join(lines[4:])))
            np.savez_compressed(cache, header=np.array(header), arr=arr)
        z = np.load(cache)
        header, arr = list(z["header"]), z["arr"]
        lt, lg, av = (np.unique(arr[:, header.index(k)]) for k in ("lgTef", "logg", "Av"))
        self.filters = [h for h in header if h.startswith(("ACS_", "WFC3_"))]
        shape = (len(lt), len(lg), len(av))
        # the table is ordered Teff, logg, Av (fastest); check rather than assume
        order = np.lexsort((arr[:, header.index("Av")], arr[:, header.index("logg")],
                            arr[:, header.index("lgTef")]))
        arr = arr[order]
        self._interp = {}
        for fname in self.filters:
            bc = arr[:, header.index(fname)].reshape(shape)
            dimming = bc[:, :, :1] - bc  # A_X at each A_V
            self._interp[fname] = RegularGridInterpolator((lt, lg, av), dimming,
                                                          bounds_error=False, fill_value=None)
        self.composition = (float(f), float(a))

    def __call__(self, band: str, logTe, logg, A_V: float):
        key = band if band in self._interp else next(k for k in self._interp if k.endswith(band))
        pts = np.column_stack([np.ravel(logTe), np.ravel(logg), np.full(np.size(logTe), A_V)])
        return self._interp[key](pts).reshape(np.shape(logTe))

    def ratio(self, band: str, logTe: float, logg: float, ebv: float = 0.05) -> float:
        """A_X / E(B-V) for a star of this temperature and gravity."""
        return float(self(band, np.array([logTe]), np.array([logg]), C.R_V * ebv)[0] / ebv)
