"""Geometry on the sphere: equal-area map projections, polygon areas, and which
country a point lies in.

Why equal-area projections? Every map in this project is a map of *density* --
meteorites per unit area. On an equal-area map, equal areas on the ground are
equal areas on paper, so a patch of dots means the same thing wherever it
sits. On the familiar Mercator map Antarctica is inflated without limit, which
would exaggerate the very collection bias we are trying to measure. We use:

* **Equal Earth** (Savric, Patterson & Jenny 2018, Int. J. Geogr. Inf. Sci. 33,
  454) for world maps -- a modern equal-area projection with a closed-form
  formula, so no mapping library is needed;
* **Lambert azimuthal equal-area**, south polar aspect, for Antarctica.

Both are exactly equal-area: the area of any region on the unit sphere equals
its area on the projected plane. ``test_physics.py`` checks this numerically,
and ``polygon_area_km2`` relies on it to measure countries.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from matplotlib.path import Path as MplPath

from config import EARTH_RADIUS_KM

# --- Equal Earth -----------------------------------------------------------------

#: Polynomial coefficients from Savric, Patterson & Jenny (2018), Eq. 2.
_A1, _A2, _A3, _A4 = 1.340264, -0.081106, 0.000893, 0.003796
_M = np.sqrt(3.0) / 2.0


def _ee_y(theta):
    t2 = theta * theta
    t6 = t2 * t2 * t2
    return theta * (_A1 + _A2 * t2 + t6 * (_A3 + _A4 * t2))


def _ee_dy(theta):
    """d y / d theta: the denominator of the x formula."""
    t2 = theta * theta
    t6 = t2 * t2 * t2
    return _A1 + 3 * _A2 * t2 + t6 * (7 * _A3 + 9 * _A4 * t2)


def equal_earth(lon_deg, lat_deg):
    """Project longitude/latitude (degrees) to Equal Earth x, y on the unit sphere.

    The map spans x in [-2.7064, 2.7064] and y in [-1.3173, 1.3173]; its total
    area is 4 pi, the area of the unit sphere.
    """
    lam = np.radians(np.asarray(lon_deg, dtype=float))
    phi = np.radians(np.asarray(lat_deg, dtype=float))
    theta = np.arcsin(_M * np.sin(phi))
    x = 2.0 * np.sqrt(3.0) * lam * np.cos(theta) / (3.0 * _ee_dy(theta))
    y = _ee_y(theta)
    return x, y


def equal_earth_inverse(x, y, *, iterations: int = 12):
    """Invert ``equal_earth`` by Newton's method on the parametric latitude."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    theta = y / _A1
    for _ in range(iterations):
        theta = theta - (_ee_y(theta) - y) / _ee_dy(theta)
    lam = 3.0 * x * _ee_dy(theta) / (2.0 * np.sqrt(3.0) * np.cos(theta))
    phi = np.arcsin(np.clip(np.sin(theta) / _M, -1.0, 1.0))
    return np.degrees(lam), np.degrees(phi)


# --- Lambert azimuthal equal-area, south polar aspect ------------------------------

def south_polar(lon_deg, lat_deg):
    """Project to Lambert azimuthal equal-area centred on the South Pole (unit sphere).

    The Greenwich meridian points up and 90 E to the right, the usual orientation
    for maps of Antarctica. ``rho`` is the projected distance from the pole:
    2 sin(c / 2) for a point at angular distance c from the pole.
    """
    lam = np.radians(np.asarray(lon_deg, dtype=float))
    phi = np.radians(np.asarray(lat_deg, dtype=float))
    rho = 2.0 * np.sin((np.pi / 2 + phi) / 2.0)
    return rho * np.sin(lam), rho * np.cos(lam)


# --- Areas ----------------------------------------------------------------------------

def shoelace(x, y) -> float:
    """Signed area of a closed ring in the plane (positive if counter-clockwise)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def ring_area_km2(lon, lat) -> float:
    """Area of a lon/lat ring on the Earth (km^2), measured on the Equal Earth plane.

    Because the projection preserves area, the planar area of the projected ring
    *is* the spherical area (to the accuracy with which straight projected edges
    follow the true edges -- negligible for Natural Earth's closely spaced vertices).
    """
    x, y = equal_earth(lon, lat)
    return abs(shoelace(x, y)) * EARTH_RADIUS_KM**2


def band_area_km2(lat1: float, lat2: float, lon_span_deg: float = 360.0) -> float:
    """Exact area of the latitude band lat1..lat2 over ``lon_span_deg`` of longitude."""
    s1, s2 = np.sin(np.radians([lat1, lat2]))
    return float(EARTH_RADIUS_KM**2 * np.radians(lon_span_deg) * abs(s2 - s1))


EARTH_AREA_KM2 = 4.0 * np.pi * EARTH_RADIUS_KM**2


# --- Countries --------------------------------------------------------------------------

@dataclass
class Country:
    """One Natural Earth admin-0 unit: identifiers plus its polygons."""

    name: str
    adm0_a3: str
    iso3: str                      # ISO 3166 alpha-3, as the World Bank uses it
    continent: str
    pop_est: float
    #: Each polygon is (exterior ring, [hole rings]); rings are (N, 2) lon/lat arrays.
    polygons: list = field(repr=False)
    area_km2: float = 0.0
    bbox: tuple = (0.0, 0.0, 0.0, 0.0)  # lon_min, lat_min, lon_max, lat_max

    def contains(self, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
        """Boolean mask: which of the points lie inside this country."""
        lon = np.asarray(lon, dtype=float)
        lat = np.asarray(lat, dtype=float)
        inside = np.zeros(lon.shape, dtype=bool)
        x0, y0, x1, y1 = self.bbox
        cand = (lon >= x0) & (lon <= x1) & (lat >= y0) & (lat <= y1)
        if not cand.any():
            return inside
        pts = np.column_stack([lon[cand], lat[cand]])
        hit = np.zeros(len(pts), dtype=bool)
        for exterior, holes in self.polygons:
            in_poly = MplPath(exterior).contains_points(pts)
            for hole in holes:
                in_poly &= ~MplPath(hole).contains_points(pts)
            hit |= in_poly
        inside[np.flatnonzero(cand)] = hit
        return inside


def _polygons_from_geometry(geom: dict) -> list:
    if geom["type"] == "Polygon":
        parts = [geom["coordinates"]]
    elif geom["type"] == "MultiPolygon":
        parts = geom["coordinates"]
    else:
        raise ValueError(f"unexpected geometry type {geom['type']}")
    polys = []
    for rings in parts:
        arrs = [np.asarray(r, dtype=float)[:, :2] for r in rings]
        polys.append((arrs[0], arrs[1:]))
    return polys


def load_countries(path: Path) -> list[Country]:
    """Read a Natural Earth admin-0 GeoJSON file into ``Country`` objects."""
    features = json.loads(Path(path).read_text())["features"]
    out = []
    for f in features:
        p = f["properties"]
        polys = _polygons_from_geometry(f["geometry"])
        area = sum(ring_area_km2(*ext.T) - sum(ring_area_km2(*h.T) for h in holes)
                   for ext, holes in polys)
        allpts = np.vstack([ext for ext, _ in polys])
        iso3 = p.get("ISO_A3_EH") or p.get("ISO_A3")
        if not iso3 or iso3 == "-99":
            iso3 = p["ADM0_A3"]
        out.append(Country(
            name=p["NAME"], adm0_a3=p["ADM0_A3"], iso3=iso3,
            continent=p["CONTINENT"], pop_est=float(p.get("POP_EST") or 0.0),
            polygons=polys, area_km2=area,
            bbox=(allpts[:, 0].min(), allpts[:, 1].min(), allpts[:, 0].max(), allpts[:, 1].max()),
        ))
    return out


def destination(lon_deg, lat_deg, bearing_deg, distance_km):
    """The point reached by travelling ``distance_km`` along a great circle."""
    lam1, phi1 = np.radians(lon_deg), np.radians(lat_deg)
    b = np.radians(bearing_deg)
    d = distance_km / EARTH_RADIUS_KM
    phi2 = np.arcsin(np.sin(phi1) * np.cos(d) + np.cos(phi1) * np.sin(d) * np.cos(b))
    lam2 = lam1 + np.arctan2(np.sin(b) * np.sin(d) * np.cos(phi1),
                             np.cos(d) - np.sin(phi1) * np.sin(phi2))
    lon2 = (np.degrees(lam2) + 540.0) % 360.0 - 180.0
    return lon2, np.degrees(phi2)


def assign_countries(lon, lat, countries: list[Country], *,
                     snap_km: tuple = (2.0, 5.0, 10.0, 25.0), directions: int = 16):
    """Index into ``countries`` for each point, or -1 if it is out at sea.

    Points are first tested against the polygons as given. A point that lands
    just offshore -- a meteorite on a beach, or coordinates rounded to 0.01
    degree, against a 1:50m coastline -- is snapped to the nearest country found
    on rings of ``snap_km`` radius around it. Returns (index, snap distance in km;
    0 for a direct hit, NaN if never matched).
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    idx = np.full(lon.shape, -1, dtype=int)
    for i, c in enumerate(countries):
        todo = idx < 0
        if not todo.any():
            break
        hit = c.contains(lon[todo], lat[todo])
        idx[np.flatnonzero(todo)[hit]] = i
    snapped = np.where(idx >= 0, 0.0, np.nan)

    bearings = np.arange(directions) * 360.0 / directions
    for r in snap_km:
        todo = np.flatnonzero(idx < 0)
        if len(todo) == 0:
            break
        # Test every point on the ring; take the country hit most often.
        votes = np.full((len(todo), directions), -1, dtype=int)
        for j, b in enumerate(bearings):
            lo, la = destination(lon[todo], lat[todo], b, r)
            sub = np.full(len(todo), -1, dtype=int)
            for i, c in enumerate(countries):
                free = sub < 0
                if not free.any():
                    break
                hit = c.contains(lo[free], la[free])
                sub[np.flatnonzero(free)[hit]] = i
            votes[:, j] = sub
        for k, row in enumerate(votes):
            got = row[row >= 0]
            if len(got):
                idx[todo[k]] = np.bincount(got).argmax()
                snapped[todo[k]] = r
    return idx, snapped


# --- Map furniture -----------------------------------------------------------------------

def graticule(step_deg: float = 30.0, *, n: int = 181):
    """Meridians and parallels as lists of lon/lat lines (for drawing)."""
    lines = []
    lat = np.linspace(-90, 90, n)
    for lon0 in np.arange(-180 + step_deg, 180, step_deg):
        lines.append((np.full_like(lat, lon0), lat))
    lon = np.linspace(-180, 180, 2 * n)
    for lat0 in np.arange(-90 + step_deg, 90, step_deg):
        lines.append((lon, np.full_like(lon, lat0)))
    return lines


def world_outline(n: int = 361):
    """The Equal Earth map's edge (the +-180 meridians), as a closed lon/lat ring."""
    lat = np.linspace(-90, 90, n)
    lon = np.concatenate([np.full(n, 180.0), np.full(n, -180.0)])
    return lon, np.concatenate([lat, lat[::-1]])
