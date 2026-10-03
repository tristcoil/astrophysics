"""The formulas: distances, magnitudes, coordinates, temperatures and star colours.

Kept apart from the pipeline so each can be checked against textbook cases
(see test_physics.py). A unit slip in any of these would put every star in the
wrong place while still producing a perfectly plausible-looking map.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

#: Light years per parsec (IAU: 1 pc = 648000/pi au; 1 ly = 63241.077 au).
LY_PER_PC = 3.261563777

# --------------------------------------------------------------------------
# Distances and magnitudes
# --------------------------------------------------------------------------


def distance_pc(parallax_mas):
    """Distance in parsecs from parallax in milliarcseconds.

    This is the definition of the parsec: the distance at which one astronomical
    unit subtends one arcsecond. It is an unbiased distance estimate only when
    the parallax is precise; the sample requires 10% precision, where the bias
    of 1/parallax is below 1%.
    """
    return 1000.0 / np.asarray(parallax_mas, dtype=float)


def absolute_magnitude(apparent_mag, parallax_mas):
    """Absolute magnitude: how bright the star would look from 10 pc.

    M = m - 5 log10(d / 10 pc) = m + 5 log10(parallax_mas) - 10.
    No extinction term: within 100 pc we sit inside the Local Bubble, a cavity
    almost empty of dust, where reddening is typically below 0.02 mag.
    """
    return np.asarray(apparent_mag, dtype=float) + 5.0 * np.log10(np.asarray(parallax_mas, dtype=float)) - 10.0


def limiting_distance_pc(absolute_mag, limiting_apparent_mag):
    """How far away a star of this absolute magnitude can be and still be seen.

    Inverts the distance modulus: d = 10^((m_lim - M + 5) / 5) pc. This single
    line is the whole story of selection effects in a magnitude-limited survey.
    """
    return 10.0 ** ((limiting_apparent_mag - np.asarray(absolute_mag, dtype=float) + 5.0) / 5.0)


def corrected_excess_factor(excess_factor, bp_rp, g_mag):
    """Gaia's corrected BP/RP flux excess C* and its expected scatter (Riello et al. 2021).

    Gaia measures BP and RP by dispersing light over a small window. The summed
    BP + RP flux should match G; when it exceeds G, extra light from a
    neighbouring source or the background has leaked into the BP/RP window, and
    the colour is contaminated -- typically too blue for a faint red star. C*
    removes the colour dependence of the raw excess factor (their eq. 6), and
    sigma_C* is its expected 1-sigma scatter at each magnitude (their eq. 18).
    """
    x = np.asarray(bp_rp, dtype=float)
    f = np.where(x < 0.5, 1.154360 + 0.033772 * x + 0.032277 * x**2,
                 np.where(x < 4.0, 1.162004 + 0.011464 * x + 0.049255 * x**2 - 0.005879 * x**3,
                          1.057572 + 0.140537 * x))
    sigma = 0.0059898 + 8.817481e-12 * np.asarray(g_mag, dtype=float) ** 7.618399
    return np.asarray(excess_factor, dtype=float) - f, sigma


# --------------------------------------------------------------------------
# Coordinates
# --------------------------------------------------------------------------


def galactic_xyz(ra_deg, dec_deg, dist_pc):
    """Heliocentric Galactic Cartesian coordinates, in parsecs.

    X points to the Galactic centre, Y along the direction the Sun orbits the
    Galaxy (l = 90 deg), Z to the North Galactic Pole. The origin is the Sun, not
    the Galactic mid-plane, which lies some 20 pc below us.

    Returns (x, y, z, l_deg, b_deg).
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    c = SkyCoord(ra=np.asarray(ra_deg) * u.deg, dec=np.asarray(dec_deg) * u.deg,
                 distance=np.asarray(dist_pc) * u.pc, frame="icrs").galactic
    cart = c.cartesian
    return (cart.x.to_value(u.pc), cart.y.to_value(u.pc), cart.z.to_value(u.pc),
            c.l.deg, c.b.deg)


def propagate_position(ra_deg, dec_deg, pmra_masyr, pmdec_masyr, years):
    """Move a position along its proper motion, to first order.

    pmra is Gaia/Hipparcos convention: already multiplied by cos(dec). Good to
    well under an arcsecond over the 25 years between Hipparcos and Gaia even
    for the fastest star in the sky, which is all the cross-match needs.
    """
    dec = np.asarray(dec_deg, dtype=float)
    ddec = np.asarray(pmdec_masyr) * years / 3.6e6
    dra = np.asarray(pmra_masyr) * years / 3.6e6 / np.cos(np.radians(dec))
    return np.asarray(ra_deg) + dra, dec + ddec


# --------------------------------------------------------------------------
# The dwarf sequence: colour <-> temperature <-> spectral type <-> luminosity
# --------------------------------------------------------------------------

SPECTRAL_CLASSES = ["O", "B", "A", "F", "G", "K", "M", "L", "T", "Y"]


@dataclass
class DwarfSequence:
    """Mamajek's mean main-sequence table, as interpolation-ready arrays.

    Every row is one spectral subtype of a main-sequence (dwarf) star. Missing
    entries are NaN. The table is empirical: each colour is the mean observed
    for well-classified stars of that type.
    """

    spt: np.ndarray
    teff: np.ndarray
    b_v: np.ndarray
    g_v: np.ndarray
    bp_rp: np.ndarray
    m_g: np.ndarray

    # -- parsing ---------------------------------------------------------

    @classmethod
    def from_text(cls, text: str) -> "DwarfSequence":
        """Parse the table body: everything between the first two '#SpT' headers."""
        lines = text.splitlines()
        header_idx = [i for i, ln in enumerate(lines) if ln.startswith("#SpT")]
        header = lines[header_idx[0]].lstrip("#").split()
        body = lines[header_idx[0] + 1: header_idx[1]]

        def number(v: str) -> float:
            v = v.rstrip(":")  # a trailing colon marks an uncertain value
            try:
                return float(v)
            except ValueError:
                return np.nan

        rows = [ln.split() for ln in body if ln.strip() and not ln.startswith("#")]
        col = {name: i for i, name in enumerate(header)}

        def column(name):
            return np.array([number(r[col[name]]) for r in rows])

        return cls(
            spt=np.array([r[0] for r in rows]),
            teff=column("Teff"), b_v=column("B-V"), g_v=column("G-V"),
            bp_rp=column("Bp-Rp"), m_g=column("M_G"),
        )

    # -- helpers ---------------------------------------------------------

    def _monotonic(self, x: np.ndarray, *ys: np.ndarray):
        """Rows usable for interpolating ys against x, with x strictly rising.

        Walks down the spectral sequence and keeps a row only if its x exceeds
        every x seen so far. That drops the few late-M rows where the tabulated
        Bp-Rp *falls* again -- Gaia's BP band runs out of photons for the
        coolest dwarfs -- instead of letting them break the interpolation.
        """
        ok = np.isfinite(x)
        for y in ys:
            ok &= np.isfinite(y)
        keep, best = [], -np.inf
        for i in np.flatnonzero(ok):
            if x[i] > best:
                keep.append(i)
                best = x[i]
        keep = np.array(keep)
        return (x[keep],) + tuple(y[keep] for y in ys)

    # -- the relations ---------------------------------------------------

    def teff_from_bp_rp(self, bp_rp):
        x, t = self._monotonic(self.bp_rp, self.teff)
        return np.interp(bp_rp, x, t, left=np.nan, right=np.nan)

    def teff_from_b_v(self, b_v):
        x, t = self._monotonic(self.b_v, self.teff)
        return np.interp(b_v, x, t, left=np.nan, right=np.nan)

    def bp_rp_from_b_v(self, b_v):
        x, y = self._monotonic(self.b_v, self.bp_rp)
        return np.interp(b_v, x, y, left=np.nan, right=np.nan)

    def g_minus_v_from_b_v(self, b_v):
        x, y = self._monotonic(self.b_v, self.g_v)
        return np.interp(b_v, x, y, left=np.nan, right=np.nan)

    def g_minus_v_from_bp_rp(self, bp_rp):
        x, y = self._monotonic(self.bp_rp, self.g_v)
        return np.interp(bp_rp, x, y, left=np.nan, right=np.nan)

    def main_sequence_mg(self, bp_rp):
        """Absolute G magnitude of a dwarf of this colour: the main-sequence locus."""
        x, y = self._monotonic(self.bp_rp, self.m_g)
        return np.interp(bp_rp, x, y, left=np.nan, right=np.nan)

    def main_sequence_mg_from_teff(self, teff):
        """Absolute G magnitude of a dwarf of this temperature.

        The same locus as ``main_sequence_mg``, but indexed by temperature so it
        also works for stars whose temperature came from B-V (Hipparcos) rather
        than from a Gaia colour. Interpolated in log T, which is smooth.
        """
        x, y = self._monotonic(-np.log10(self.teff), self.m_g)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.interp(-np.log10(np.asarray(teff, dtype=float)), x, y, left=np.nan, right=np.nan)

    def teff_from_mg(self, m_g):
        """Temperature of the dwarf with this absolute magnitude.

        Used only for the coolest objects (L/T dwarfs), where Gaia's own colour is
        too noisy to trust but the absolute magnitude still is.
        """
        x, t = self._monotonic(self.m_g, self.teff)
        return np.interp(m_g, x, t, left=np.nan, right=np.nan)

    @property
    def teff_range_bp_rp(self) -> tuple[float, float]:
        x, _ = self._monotonic(self.bp_rp, self.teff)
        return float(x[0]), float(x[-1])

    def class_boundaries(self) -> dict[str, tuple[float, float]]:
        """Temperature range of each spectral class, derived from the table.

        A boundary sits midway between the coolest subtype of one class and the
        hottest subtype of the next (e.g. K/M between K9V at 3930 K and M0V at
        3850 K), so no boundary is hand-picked.
        """
        letters = np.array([s[0] for s in self.spt])
        bounds = {}
        present = [c for c in SPECTRAL_CLASSES if c in letters]
        for i, c in enumerate(present):
            hot = np.nanmax(self.teff[letters == c])
            cool = np.nanmin(self.teff[letters == c])
            upper = np.inf if i == 0 else 0.5 * (hot + np.nanmin(self.teff[letters == present[i - 1]]))
            lower = 0.0 if i == len(present) - 1 else 0.5 * (cool + np.nanmax(self.teff[letters == present[i + 1]]))
            bounds[c] = (lower, upper)
        return bounds

    def spectral_class(self, teff):
        """Spectral class letter for each temperature (dwarf scale)."""
        teff = np.asarray(teff, dtype=float)
        out = np.full(teff.shape, "", dtype="<U1")
        for c, (lo, hi) in self.class_boundaries().items():
            out[(teff >= lo) & (teff < hi)] = c
        return out

    def spectral_subtype(self, teff):
        """Nearest tabulated subtype (e.g. 'M3.5V') for each temperature."""
        teff = np.asarray(teff, dtype=float)
        ok = np.isfinite(self.teff)
        table_t, table_s = self.teff[ok], self.spt[ok]
        out = np.full(teff.shape, "", dtype="<U8")
        good = np.isfinite(teff)
        idx = np.abs(np.log(teff[good, None]) - np.log(table_t[None, :])).argmin(axis=1)
        out[good] = table_s[idx]
        return out


# --------------------------------------------------------------------------
# What colour is a star? Blackbody -> human eye -> screen
# --------------------------------------------------------------------------

H = 6.62607015e-34   # Planck constant, J s
C_LIGHT = 2.99792458e8  # m/s
K_B = 1.380649e-23   # Boltzmann constant, J/K
WIEN_B = 2.897771955e-3  # Wien displacement constant, m K


def planck_lambda(wavelength_nm, teff):
    """Blackbody spectral radiance B_lambda(T), arbitrary overall scale.

    Only the *shape* matters for colour, so constants that scale every
    wavelength equally are kept but the absolute units are irrelevant.
    """
    lam = np.asarray(wavelength_nm, dtype=float) * 1e-9
    t = np.asarray(teff, dtype=float)
    with np.errstate(over="ignore"):
        return (2 * H * C_LIGHT**2 / lam**5) / np.expm1(H * C_LIGHT / (lam * K_B * t))


def _lobe(lam, mu, s1, s2):
    s = np.where(lam < mu, s1, s2)
    return np.exp(-0.5 * ((lam - mu) / s) ** 2)


def cie_xyz_bar(wavelength_nm):
    """CIE 1931 2-degree colour-matching functions (x-bar, y-bar, z-bar).

    Uses the multi-lobe analytic fit of Wyman, Sloan & Shirley (2013, JCGT 2,
    1), which reproduces the tabulated CIE functions to about 1% -- far better
    than the colour differences between stars. These three curves are how the
    average human eye responds to light; y-bar is the eye's brightness response.
    """
    lam = np.asarray(wavelength_nm, dtype=float)
    x = (1.056 * _lobe(lam, 599.8, 37.9, 31.0) + 0.362 * _lobe(lam, 442.0, 16.0, 26.7)
         - 0.065 * _lobe(lam, 501.1, 20.4, 26.2))
    y = 0.821 * _lobe(lam, 568.8, 46.9, 40.5) + 0.286 * _lobe(lam, 530.9, 16.3, 31.1)
    z = 1.217 * _lobe(lam, 437.0, 11.8, 36.0) + 0.681 * _lobe(lam, 459.0, 26.0, 13.8)
    return x, y, z


#: Linear sRGB from CIE XYZ, D65 white point (IEC 61966-2-1).
XYZ_TO_LINEAR_SRGB = np.array([
    [3.2404542, -1.5371385, -0.4985314],
    [-0.9692660, 1.8760108, 0.0415560],
    [0.0556434, -0.2040259, 1.0572252],
])


def _srgb_encode(linear):
    linear = np.clip(linear, 0.0, 1.0)
    return np.where(linear <= 0.0031308, 12.92 * linear, 1.055 * linear ** (1 / 2.4) - 0.055)


def blackbody_srgb(teff):
    """The colour of a blackbody at temperature T, as sRGB in [0, 1].

    Integrates the Planck spectrum against the eye's colour-matching functions,
    converts to screen RGB, and scales so the brightest channel is 1 -- i.e.
    this is the star's *colour*, deliberately separated from its brightness. The
    reddest temperatures fall outside what a screen can show; those are pulled
    toward white just enough to fit, which keeps their hue.
    """
    teff = np.atleast_1d(np.asarray(teff, dtype=float))
    lam = np.arange(380.0, 781.0, 1.0)
    xb, yb, zb = cie_xyz_bar(lam)
    spec = planck_lambda(lam[None, :], teff[:, None])
    xyz = np.stack([(spec * xb).sum(1), (spec * yb).sum(1), (spec * zb).sum(1)], axis=1)
    rgb = xyz @ XYZ_TO_LINEAR_SRGB.T
    rgb -= np.minimum(rgb.min(axis=1, keepdims=True), 0.0)  # desaturate out-of-gamut
    rgb /= rgb.max(axis=1, keepdims=True)
    return _srgb_encode(rgb)


def srgb_to_hex(rgb) -> list[str]:
    rgb = np.clip(np.round(np.atleast_2d(rgb) * 255), 0, 255).astype(int)
    return ["#{:02x}{:02x}{:02x}".format(*row) for row in rgb]


class StarColours:
    """Fast temperature -> colour lookup for hundreds of thousands of stars.

    Computing a spectrum per star would be wasteful; colour varies smoothly with
    temperature, so a table on a fine log-temperature grid is interpolated.
    """

    def __init__(self, t_min: float = 1000.0, t_max: float = 40000.0, n: int = 400):
        self.log_t = np.linspace(np.log10(t_min), np.log10(t_max), n)
        self.rgb = blackbody_srgb(10 ** self.log_t)

    def __call__(self, teff):
        lt = np.log10(np.clip(np.asarray(teff, dtype=float), 10 ** self.log_t[0], 10 ** self.log_t[-1]))
        return np.stack([np.interp(lt, self.log_t, self.rgb[:, k]) for k in range(3)], axis=-1)


# --------------------------------------------------------------------------
# Names
# --------------------------------------------------------------------------

_GREEK = {
    "alf": "α", "bet": "β", "gam": "γ", "del": "δ", "eps": "ε", "zet": "ζ", "eta": "η",
    "tet": "θ", "iot": "ι", "kap": "κ", "lam": "λ", "mu.": "μ", "nu.": "ν", "ksi": "ξ",
    "omi": "ο", "pi.": "π", "rho": "ρ", "sig": "σ", "tau": "τ", "ups": "υ", "phi": "φ",
    "chi": "χ", "psi": "ψ", "ome": "ω",
}


_GREEK_ID = re.compile(r"^(" + "|".join(re.escape(k) for k in _GREEK) + r")(\d*)$")


def tidy_simbad_name(name: str) -> str:
    """Turn SIMBAD's catalogue-style identifiers into something readable.

    '* alf CMa' -> 'α CMa';  '* ome02 Tau' -> 'ω2 Tau';  'NAME Barnard's star'
    -> 'Barnard's star';  'Wolf  359' -> 'Wolf 359'. SIMBAD pads identifiers
    with spaces and marks stars with a leading '*' or 'V*'; neither helps a
    human reader.
    """
    name = name.strip()
    for prefix in ("NAME ", "* ", "V* "):
        if name.startswith(prefix):
            name = name[len(prefix):]
    parts = name.split()
    if parts:
        m = _GREEK_ID.match(parts[0])
        if m:
            parts[0] = _GREEK[m.group(1)] + m.group(2).lstrip("0")
    return " ".join(parts)
