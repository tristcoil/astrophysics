"""Published results, read back from the papers' own source files.

Used only to grade our measurements at the end; the analysis never starts from
them.

**de Blok et al. (2008)** published their THINGS rotation curves only as
figures. Those figures are vector PostScript written by SuperMongo (SM), in
which every plotted symbol is a small polygon drawn at its exact position in
units of 1/300 inch. So instead of digitizing by eye, we interpret the drawing
commands: the grey filled decagons are the THINGS points, the grey vertical
strokes through them their error bars, and the black open symbols the older
curves they were compared with. The axes are calibrated from the plot's own
tick marks. On a 180 km/s axis 2,000 units tall, positions are recovered to
~0.1 km/s, far finer than any error bar.

**Trachternach et al. (2008)**: harmonic decompositions of the same velocity
fields, from the LaTeX table. **de Blok et al. (2008)**: the adopted centres,
systemic velocities, mean inclinations and position angles, from the LaTeX
table. **Walter et al. (2008)**, via VizieR: distances, and global HI spectra.
**SPARC** (Lelli et al. 2016), as in project D4.
"""

from __future__ import annotations

import re
import tarfile
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np

from config import REF_DIR

DEBLOK_TAR = REF_DIR / "deblok2008_arxiv0810.2100.tar.gz"
TRACHTERNACH_TAR = REF_DIR / "trachternach2008_arxiv0810.2116.tar.gz"
OH2011_TAR = REF_DIR / "oh2011_arxiv1011.0899.tar.gz"


# --------------------------------------------------------------------------
# A minimal interpreter for SuperMongo's PostScript
# --------------------------------------------------------------------------

@dataclass
class Shape:
    """One drawn path: a filled polygon or a stroked polyline."""

    kind: str  # "fill" or "stroke"
    color: tuple[int, int, int]  # 0-255
    points: np.ndarray  # (n, 2) device coordinates

    @property
    def centre(self) -> tuple[float, float]:
        lo, hi = self.points.min(axis=0), self.points.max(axis=0)
        return float((lo[0] + hi[0]) / 2), float((lo[1] + hi[1]) / 2)

    @property
    def size(self) -> float:
        lo, hi = self.points.min(axis=0), self.points.max(axis=0)
        return float(max(hi - lo))


@dataclass
class SMPlot:
    shapes: list[Shape] = field(default_factory=list)

    def filled(self, color) -> list[Shape]:
        return [s for s in self.shapes if s.kind == "fill" and s.color == tuple(color)]

    def strokes(self, color) -> list[Shape]:
        return [s for s in self.shapes if s.kind == "stroke" and s.color == tuple(color)]


_TOKEN = re.compile(r"\(.*?\)|\[|\]|[^\s\[\]]+")


def read_sm_postscript(text: str) -> SMPlot:
    """Interpret the path-drawing subset of an SM PostScript file.

    SM's prologue defines one-letter procedures (``M`` moveto, ``m`` rmoveto,
    ``L`` lineto, ``l`` rlineto, ``B`` stroke-and-restart, ``CF`` fill,
    ``C`` set colour, ...). We replay them, tracking the current point, the
    current path and the colour, and record every path that gets stroked or
    filled. Text is drawn by SM as strokes too; it is recorded but unused.
    """
    body = text[text.index("translate") + len("translate"):]
    plot = SMPlot()
    stack: list = []
    cur = np.zeros(2)
    color = (0, 0, 0)
    subpaths: list[list[np.ndarray]] = []

    def new_subpath(at: np.ndarray) -> None:
        subpaths.append([at.copy()])

    def emit(kind: str) -> None:
        for sp in subpaths:
            if len(sp) > 1:
                plot.shapes.append(Shape(kind, color, np.array(sp)))
        subpaths.clear()

    def pop_point() -> np.ndarray:
        y = stack.pop()
        x = stack.pop()
        return np.array([x, y], dtype=float)

    for tok in _TOKEN.findall(body):
        try:
            stack.append(float(tok))
            continue
        except ValueError:
            pass
        if tok == "M":
            cur = pop_point()
            new_subpath(cur)
        elif tok == "m":
            cur = cur + pop_point()
            new_subpath(cur)
        elif tok in ("L", "l"):
            p = pop_point()
            cur = p if tok == "L" else cur + p
            if not subpaths:
                new_subpath(cur)
            subpaths[-1].append(cur.copy())
        elif tok == "B":  # stroke what is there, restart the path here
            emit("stroke")
            new_subpath(cur)
        elif tok == "CF":  # fill, leaving the current point on the stack
            emit("fill")
            stack.extend(cur.tolist())
        elif tok == "CS":
            emit("stroke")
            stack.extend(cur.tolist())
        elif tok == "C":  # stroke with the old colour, then switch colour
            b, g, r = stack.pop(), stack.pop(), stack.pop()
            emit("stroke")
            new_subpath(cur)
            color = (int(r), int(g), int(b))
        elif tok == "stroke":
            emit("stroke")
        elif tok == "[":
            stack.append("[")
        elif tok == "]":
            while stack and stack[-1] != "[":
                stack.pop()
            if stack:
                stack.pop()
            stack.append("array")
        elif tok == "setdash":
            del stack[-2:]
        elif tok in ("setlinewidth", "P"):
            if stack:
                stack.pop()
        elif tok.startswith("("):
            stack.append(tok)
        elif tok == "T":
            if stack:
                stack.pop()
        # gsave, grestore, showpage, end, newpath: no effect on what we record
    emit("stroke")
    return plot


@dataclass(frozen=True)
class Axis:
    """A linear map from device units to data units along one axis."""

    slope: float
    offset: float
    tick_rms: float  # how far the major ticks sit from the fitted map, data units

    def __call__(self, device):
        return self.slope * np.asarray(device, dtype=float) + self.offset


def frames(plot: SMPlot, min_size: float = 400.0) -> list[Shape]:
    """The plot frames (large black rectangles), left to right, top to bottom."""
    def is_rectangle(s: Shape) -> bool:
        lo, hi = s.points.min(axis=0), s.points.max(axis=0)
        on_edge = (np.isclose(s.points[:, 0], lo[0], atol=1) | np.isclose(s.points[:, 0], hi[0], atol=1)
                   | np.isclose(s.points[:, 1], lo[1], atol=1) | np.isclose(s.points[:, 1], hi[1], atol=1))
        return bool(on_edge.all())

    out = [s for s in plot.strokes((0, 0, 0)) if len(s.points) >= 4 and np.ptp(s.points[:, 0]) > min_size
           and np.ptp(s.points[:, 1]) > min_size and is_rectangle(s)]
    return sorted(out, key=lambda s: (-s.points[:, 1].min(), s.points[:, 0].min()))


def calibrate_axes(plot: SMPlot, x_ticks: list[float], y_ticks: list[float],
                   frame: Shape | None = None) -> tuple[Axis, Axis]:
    """Find the frame and its major tick marks; map them to the given values.

    The frame is the largest black rectangle. Major ticks are the longer of
    the short black strokes rising from its bottom edge (x) or running from
    its left edge (y); they are matched in order to ``x_ticks``/``y_ticks``,
    the values printed under them (read off the rendered figure once and
    written into ``DEBLOK_FIGURES``).
    """
    black = plot.strokes((0, 0, 0))
    if frame is None:
        rects = [s for s in black if len(s.points) >= 4]
        frame = max(rects, key=lambda s: np.ptp(s.points[:, 0]) * np.ptp(s.points[:, 1]))
    x0, y0 = frame.points.min(axis=0)
    x1, y1 = frame.points.max(axis=0)

    def ticks(axis: int) -> np.ndarray:
        found = []
        for s in black:
            if len(s.points) != 2:
                continue
            (ax, ay), (bx, by) = s.points
            if axis == 0 and abs(ax - bx) < 0.5 and min(ay, by) == y0 and 0 < abs(by - ay) < 0.1 * (y1 - y0):
                found.append((ax, abs(by - ay)))
            if axis == 1 and abs(ay - by) < 0.5 and min(ax, bx) == x0 and 0 < abs(bx - ax) < 0.1 * (x1 - x0):
                found.append((ay, abs(bx - ax)))
        found = [f for f in found if (x0 - 1 <= f[0] <= x1 + 1 if axis == 0 else y0 - 1 <= f[0] <= y1 + 1)]
        found = np.array(sorted(set(found)))
        major = found[found[:, 1] >= 0.8 * found[:, 1].max(), 0]
        return np.unique(major)

    axes = []
    for axis, values in ((0, x_ticks), (1, y_ticks)):
        pos = ticks(axis)[: len(values)]
        if len(pos) != len(values):
            raise ValueError(f"found {len(pos)} major ticks on axis {axis}, expected {len(values)}")
        slope, offset = np.polyfit(pos, values, 1)
        rms = float(np.sqrt(np.mean((slope * pos + offset - np.array(values)) ** 2)))
        axes.append(Axis(float(slope), float(offset), rms))
    return axes[0], axes[1]


# --------------------------------------------------------------------------
# de Blok et al. (2008) rotation curves
# --------------------------------------------------------------------------

#: Which figure holds each galaxy's THINGS curve, and the tick values printed
#: on its axes (major ticks, in order from the origin).
DEBLOK_FIGURES = {
    "NGC2403": dict(member="fig7.ps", x_ticks=[0, 5, 10, 15], y_ticks=[0, 50, 100, 150],
                    others={"open circle": "Begeman (1987)", "open square": "Sicking (1997)",
                            "star": "Fraternali et al. (2002)"}),
    "NGC2903": dict(member="fig9.ps", x_ticks=[0, 10, 20, 30], y_ticks=[0, 50, 100, 150, 200, 250],
                    others={"open circle": "Begeman (1987)"}),
    "NGC3198": dict(member="fig12.ps", x_ticks=[0, 10, 20, 30], y_ticks=[0, 50, 100, 150],
                    others={"open circle": "Begeman (1987)", "open square": "Sicking (1997)"}),
    "DDO154": dict(member="fig15.ps", x_ticks=[0, 2, 4, 6, 8], y_ticks=[0, 20, 40, 60],
                   others={"open circle": "Carignan & Purton (1998)"}),
}

GREY = (150, 150, 150)


@dataclass
class PublishedCurve:
    galaxy: str
    source: str
    r_kpc: np.ndarray
    v: np.ndarray
    err_up: np.ndarray
    err_down: np.ndarray
    distance_mpc: float | None = None
    calibration_rms: tuple[float, float] = (0.0, 0.0)

    @property
    def err(self) -> np.ndarray:
        return 0.5 * (self.err_up + self.err_down)


def _tar_member(tar_path: Path, member: str) -> str:
    with tarfile.open(tar_path) as tar:
        return tar.extractfile(member).read().decode("latin-1")


@lru_cache(maxsize=None)
def deblok_plot(galaxy_key: str) -> SMPlot:
    return read_sm_postscript(_tar_member(DEBLOK_TAR, DEBLOK_FIGURES[galaxy_key]["member"]))


def deblok_rotation_curve(galaxy_key: str) -> PublishedCurve:
    """The THINGS rotation curve of de Blok et al. (2008), from its vector figure."""
    spec = DEBLOK_FIGURES[galaxy_key]
    plot = deblok_plot(galaxy_key)
    xa, ya = calibrate_axes(plot, spec["x_ticks"], spec["y_ticks"])

    # The legend repeats one grey symbol; it sits right of / below the data
    # cloud's own bars, so keep only symbols that carry an error bar.
    centres = np.array([s.centre for s in plot.filled(GREY)])
    bars = []
    for s in plot.strokes(GREY):
        for a, b in zip(s.points[:-1], s.points[1:]):
            if abs(a[0] - b[0]) < 0.5 and abs(a[1] - b[1]) > 0:
                bars.append((a[0], min(a[1], b[1]), max(a[1], b[1])))
    bars = np.array(bars) if bars else np.zeros((0, 3))

    rows = []
    for cx, cy in centres:
        near = bars[(np.abs(bars[:, 0] - cx) <= 1.0) & (bars[:, 1] <= cy + 1) & (bars[:, 2] >= cy - 1)]
        if len(near) == 0:
            continue
        top, bottom = near[:, 2].max(), near[:, 1].min()
        rows.append((cx, cy, top - cy, cy - bottom))
    rows = np.array(sorted(rows))
    r, v = xa(rows[:, 0]), ya(rows[:, 1])
    up, down = rows[:, 2] * ya.slope, rows[:, 3] * ya.slope
    return PublishedCurve(galaxy_key, "de Blok et al. 2008 (THINGS)", r, v, up, down,
                          calibration_rms=(xa.tick_rms, ya.tick_rms))


def deblok_comparison_curves(galaxy_key: str) -> dict[str, np.ndarray]:
    """The older curves de Blok et al. plotted for comparison (black open symbols).

    Returns {label: array of (r_kpc, v)} -- symbol centres only.
    """
    spec = DEBLOK_FIGURES[galaxy_key]
    plot = deblok_plot(galaxy_key)
    xa, ya = calibrate_axes(plot, spec["x_ticks"], spec["y_ticks"])
    out: dict[str, list] = {}
    for s in plot.strokes((0, 0, 0)):
        n = len(s.points)
        closed = np.allclose(s.points[0], s.points[-1], atol=1.0)
        if not (8 <= s.size <= 45):
            continue
        if n in (10, 11) and closed:
            kind = "open circle"
        elif n == 5 and closed:
            kind = "open square"
        else:
            continue
        if kind in spec["others"]:
            out.setdefault(spec["others"][kind], []).append(s.centre)
    result = {}
    for label, pts in out.items():
        pts = np.array(sorted(pts))
        result[label] = np.column_stack([xa(pts[:, 0]), ya(pts[:, 1])])
    return result


# --------------------------------------------------------------------------
# Tables from the LaTeX sources
# --------------------------------------------------------------------------

def _latex_rows(tex: str, label: str) -> list[list[str]]:
    """Rows of the deluxetable carrying ``\\label{label}``, split into cells."""
    at = tex.index(f"\\label{{{label}}}")
    env_start = tex.rindex("\\begin{deluxetable", 0, at)
    env_end = tex.index("\\end{deluxetable", at)
    start = tex.index("\\startdata", env_start, env_end)
    stop = tex.index("\\enddata", start, env_end)
    block = tex[start + len("\\startdata"):stop]
    rows = []
    for line in re.split(r"\\\\(?:\[[^\]]*\])?", block):
        line = re.sub(r"%.*", "", line).strip()
        if not line or line.startswith("\\hline"):
            continue
        rows.append([c.strip() for c in line.split("&")])
    return rows


def _num(cell: str) -> float:
    cell = cell.replace("--", "-").replace("$", "").replace("{\\bf", "").replace("}", "")
    m = re.search(r"-?\d+(\.\d+)?", cell.replace(" ", ""))
    return float(m.group()) if m else float("nan")


def _sexagesimal(cell: str, hours: bool) -> float:
    parts = [float(p) for p in cell.replace("--", "-").split()]
    sign = -1.0 if cell.strip().startswith("-") else 1.0
    value = abs(parts[0]) + parts[1] / 60 + parts[2] / 3600
    return sign * value * (15.0 if hours else 1.0)


@lru_cache(maxsize=None)
def deblok_geometry() -> dict[str, dict]:
    """de Blok et al. (2008) Table 3: centre, D, ring step, V_sys, <i>, <PA>."""
    tex = _tar_member(DEBLOK_TAR, "deblok_astroph.tex")
    out = {}
    for cells in _latex_rows(tex, "bigtable"):
        name = cells[0].replace(" ", "")
        out[name] = dict(ra_deg=_sexagesimal(cells[1], True), dec_deg=_sexagesimal(cells[2], False),
                         distance_mpc=_num(cells[3]), ring_step_arcsec=_num(cells[4]),
                         vsys=_num(cells[5]), incl=_num(cells[6]), pa=_num(cells[7]))
    return out


@lru_cache(maxsize=None)
def trachternach_harmonics() -> dict[str, dict]:
    """Trachternach et al. (2008) Table 3: median non-circular amplitudes.

    A_r is the median over rings of sqrt(s1^2 + c2^2 + s2^2 + c3^2 + s3^2)
    (km/s), whole disc and inner 1 kpc, with lower/upper quartiles; also the
    potential's elongation and the median |residual| after the decomposition.
    """
    with tarfile.open(TRACHTERNACH_TAR) as tar:
        tex = tar.extractfile("trachternach.tex").read().decode("latin-1")
    out = {}
    for cells in _latex_rows(tex, "table:harm-decomp"):
        name = cells[0].replace(" ", "").replace("{\\bf", "").strip()
        if not re.match(r"(NGC|DDO|IC)\d+", name):
            continue

        def amp(cell):
            m = re.match(r"\$?([\d.]+)_\{-([\d.]+)\}\^\{\+([\d.]+)\}", cell.replace(" ", ""))
            if m:
                return tuple(float(g) for g in m.groups())
            return (_num(cell), float("nan"), float("nan"))

        a_all, a_1kpc = amp(cells[1]), amp(re.sub(r"\\footnote\{.*", "", cells[2]))
        eps = re.findall(r"-?[\d.]+", cells[4].replace("$-$", "-"))
        out[name] = dict(A_r=a_all[0], A_r_q=(a_all[0] - a_all[1], a_all[0] + a_all[2]),
                         A_r_1kpc=a_1kpc[0], A_r_over_vmax_pct=_num(cells[3]),
                         eps_pot=float(eps[0]) if eps else float("nan"),
                         eps_pot_err=float(eps[1]) if len(eps) > 1 else float("nan"),
                         median_residual=_num(cells[5]), r_max_arcsec=_num(cells[6]))
    return out


def oh2011_ddo154_density() -> np.ndarray:
    """The five inner density points Oh et al. (2011) fitted for DDO 154.

    Read from the "Mass density profile" panel of their Fig. A.27 (SM vector
    PostScript, log-log axes): the grey filled symbols are the points their
    slope alpha = -0.29 +- 0.15 was fitted to. Returns (R kpc, rho Msun/pc^3).
    """
    with tarfile.open(OH2011_TAR) as tar:
        text = tar.extractfile("FigureA27.ps").read().decode("latin-1")
    plot = read_sm_postscript(text)
    panel = max(frames(plot), key=lambda s: s.points[:, 0].min())  # the right-hand panel
    lo, hi = panel.points.min(axis=0), panel.points.max(axis=0)
    inside = SMPlot([s for s in plot.shapes if np.all((s.points >= lo - 150) & (s.points <= hi + 5))])
    xa, ya = calibrate_axes(inside, [-1, 0], [-3, -2, -1], frame=panel)
    pts = np.array(sorted(s.centre for s in inside.filled((170, 170, 170))))
    return np.column_stack([10 ** xa(pts[:, 0]), 10 ** ya(pts[:, 1])])


#: Oh et al. (2011) printed DDO 154's inner density slope inside a figure panel
#: (their Fig. A.27, "Mass density profile"), measured from de Blok et al.'s
#: Hermite rotation curve assuming minimum disk; the sample mean is in the text.
OH2011_DDO154_ALPHA = (-0.29, 0.15)
OH2011_DWARF_MEAN_ALPHA = (-0.29, 0.07)


# --------------------------------------------------------------------------
# VizieR tables (Walter et al. 2008) and SPARC
# --------------------------------------------------------------------------

def _read_vizier_tsv(path: Path) -> list[dict]:
    lines = [ln for ln in path.read_text().splitlines() if ln and not ln.startswith("#")]
    header = lines[0].split("\t")
    rows = []
    for ln in lines[3:]:  # names, units, dashes
        cells = ln.split("\t")
        if len(cells) == len(header):
            rows.append({h: c.strip() for h, c in zip(header, cells)})
    return rows


@lru_cache(maxsize=None)
def walter_table1() -> dict[str, dict]:
    out = {}
    for row in _read_vizier_tsv(REF_DIR / "walter2008_table1.tsv"):
        out[row["Name"].replace(" ", "")] = row
    return out


def walter_distance_mpc(galaxy_key: str) -> float:
    return float(walter_table1()[galaxy_key]["Dist"])


def walter_global_spectrum(galaxy_key: str) -> tuple[np.ndarray, np.ndarray]:
    """THINGS global HI spectrum: velocity (km/s) and flux density (Jy)."""
    vel, flux = [], []
    for row in _read_vizier_tsv(REF_DIR / "walter2008_table4.tsv"):
        if row["Name"].replace(" ", "") == galaxy_key and row["SHI"]:
            vel.append(float(row["Vel"]))
            flux.append(float(row["SHI"]))
    order = np.argsort(vel)
    return np.array(vel)[order], np.array(flux)[order]


@lru_cache(maxsize=None)
def sparc_galaxy_table() -> dict[str, dict]:
    """SPARC Table 1: distance (Mpc), inclination (deg) and data references."""
    out = {}
    lines = (REF_DIR / "SPARC_Lelli2016c.mrt").read_text(encoding="latin-1").splitlines()
    for ln in lines:
        parts = ln.split()
        if len(parts) >= 19 and re.match(r"^[A-Z]", parts[0]) and parts[1].isdigit():
            try:
                out[parts[0]] = dict(distance_mpc=float(parts[2]), incl=float(parts[5]),
                                     incl_err=float(parts[6]), refs=parts[-1])
            except ValueError:
                continue
    return out


def sparc_rotation_curve(sparc_name: str, distance_mpc: float | None = None) -> PublishedCurve:
    """SPARC's rotation curve, optionally rescaled to another distance.

    Radii scale with distance; velocities do not.
    """
    rows = []
    for ln in (REF_DIR / "MassModels_Lelli2016c.mrt").read_text(encoding="latin-1").splitlines():
        parts = ln.split()
        if len(parts) >= 8 and parts[0] == sparc_name:
            rows.append([float(p) for p in parts[1:5]])  # D, R, Vobs, errV
    rows = np.array(rows)
    d_sparc = rows[0, 0]
    scale = (distance_mpc / d_sparc) if distance_mpc else 1.0
    info = sparc_galaxy_table().get(sparc_name, {})
    return PublishedCurve(sparc_name, f"SPARC ({info.get('refs', '')})", rows[:, 1] * scale,
                          rows[:, 2], rows[:, 3], rows[:, 3],
                          distance_mpc=distance_mpc or d_sparc)
