"""The quasar luminosity function: how many quasars of each luminosity per unit volume.

A flux-limited survey sees faint quasars only nearby, so counting them needs the
volume each could have been seen in. We use the Page & Carrera (2000) form of
the 1/V_a estimator, as Richards et al. (2006) did: in a bin of redshift and
absolute magnitude,

    Phi = sum_i w_i / V_a,   V_a = Omega * integral over the bin's accessible part of dV/dz dz dM,

where "accessible" means a quasar at that (z, M) would fall inside the survey's
apparent-magnitude limits, and w_i = 1 / (selection probability) corrects for
quasars the colour selection missed. Phi is per Mpc^3 per magnitude.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

import physics
from config import I_BRIGHT_LIMIT, I_LIMIT_HIGHZ, I_LIMIT_LOWZ, SELFUNC_FLOOR


class KCorrection:
    """K(z) from a table (Richards et al. 2006 Table 4, or our own), linearly interpolated."""

    def __init__(self, z, k, label: str = ""):
        order = np.argsort(z)
        self.z = np.asarray(z, float)[order]
        self.k = np.asarray(k, float)[order]
        self.label = label

    def __call__(self, z):
        return np.interp(np.asarray(z, float), self.z, self.k)


class SelectionFunction:
    """Richards et al. 2006 Table 1: completeness of the colour/radio selection vs (i, z)."""

    KINDS = ("point", "radio", "extend")

    def __init__(self, table: pd.DataFrame):
        imag = np.sort(table["imag"].unique())
        zz = np.sort(table["z"].unique())
        self.imag, self.z = imag, zz
        self._interp = {}
        for kind in self.KINDS:
            grid = (table.pivot(index="imag", columns="z", values=kind)
                    .reindex(index=imag, columns=zz).to_numpy())
            self._interp[kind] = RegularGridInterpolator((imag, zz), grid, bounds_error=False,
                                                         fill_value=None)

    def __call__(self, i, z, kind="point"):
        i = np.clip(np.asarray(i, float), self.imag[0], self.imag[-1])
        z = np.clip(np.asarray(z, float), self.z[0], self.z[-1])
        kind = np.broadcast_to(np.asarray(kind), i.shape)
        out = np.empty(i.shape)
        for k in self.KINDS:
            m = kind == k
            if np.any(m):
                out[m] = self._interp[k](np.column_stack([i[m], z[m]]))
        return np.clip(out, 0.0, 1.0)

    def marginal(self, z, i_max: float = I_LIMIT_LOWZ, kind="point"):
        """Selection function averaged over i (uniformly in magnitude, 15 <= i <= i_max)."""
        ii = np.arange(15.0, i_max + 1e-9, 0.1)
        return np.array([np.mean(self(ii, np.full_like(ii, zz), kind)) for zz in np.atleast_1d(z)])


@dataclass
class Survey:
    """Geometry and flux limits of a uniformly selected quasar sample."""

    area_deg2: float
    kcorr: KCorrection
    z_split: float                    # faint limit switches from i_faint_lo to i_faint_hi here
    i_faint_lo: float = I_LIMIT_LOWZ
    i_faint_hi: float = I_LIMIT_HIGHZ
    i_bright: float = I_BRIGHT_LIMIT
    cosmology: str = "ref"            # "ref" (h=0.7, Om=0.3) or "d3"
    name: str = ""
    meta: dict = field(default_factory=dict)

    def i_faint(self, z):
        z = np.asarray(z, float)
        return np.where(z < self.z_split, self.i_faint_lo, self.i_faint_hi)

    def m_faint(self, z):
        """Faintest accessible M_i(z=2) at each redshift."""
        return physics.absolute_mag(self.i_faint(z), z, self.kcorr, self.cosmology)

    def m_bright(self, z):
        return physics.absolute_mag(np.full_like(np.asarray(z, float), self.i_bright), z,
                                    self.kcorr, self.cosmology)

    def accessible_volume(self, z_lo, z_hi, m_lo, m_hi, nz: int = 600) -> tuple[float, float]:
        """V_a (Mpc^3 mag) of the bin, and the accessible fraction of its (z, M) area.

        m_lo < m_hi are absolute magnitudes (m_lo is the brighter edge).
        """
        edges = np.linspace(z_lo, z_hi, nz + 1)
        zc = 0.5 * (edges[1:] + edges[:-1])
        dz = np.diff(edges)
        faint = np.minimum(self.m_faint(zc), m_hi)
        bright = np.maximum(self.m_bright(zc), m_lo)
        width = np.clip(faint - bright, 0.0, None)
        dv = physics.dvdz_per_deg2(zc, self.cosmology) * self.area_deg2
        va = float(np.sum(dv * dz * width))
        frac = float(np.sum(dz * width) / ((z_hi - z_lo) * (m_hi - m_lo)))
        return va, frac


def weights_from_selection(p, floor: float = SELFUNC_FLOOR, image_completeness: float = 0.95):
    """1 / completeness. Richards et al. 2006 floor the selection function at 1/3 (their
    section 3.3) and multiply by a global 95% image-quality completeness."""
    return 1.0 / (np.maximum(np.asarray(p, float), floor) * image_completeness)


def binned_qlf(z, mag, weights, survey: Survey, z_edges, m_edges, *, min_fill: float = 0.0,
               nz: int = 600) -> pd.DataFrame:
    """Page & Carrera (2000) binned luminosity function. Poisson errors as in Richards et al."""
    z = np.asarray(z, float)
    mag = np.asarray(mag, float)
    w = np.asarray(weights, float)
    rows = []
    zi = np.digitize(z, z_edges) - 1
    mi = np.digitize(mag, m_edges) - 1
    for a in range(len(z_edges) - 1):
        for b in range(len(m_edges) - 1):
            sel = (zi == a) & (mi == b)
            n = int(sel.sum())
            if n == 0:
                continue
            va, fill = survey.accessible_volume(z_edges[a], z_edges[a + 1], m_edges[b],
                                                m_edges[b + 1], nz=nz)
            if va <= 0 or fill < min_fill:
                continue
            nw = float(w[sel].sum())
            rows.append(dict(z_lo=z_edges[a], z_hi=z_edges[a + 1], z_mid=0.5 * (z_edges[a] + z_edges[a + 1]),
                             m_lo=m_edges[b], m_hi=m_edges[b + 1], m_mid=0.5 * (m_edges[b] + m_edges[b + 1]),
                             n=n, n_weighted=nw, v_a=va, fill=fill,
                             phi=nw / va, phi_err=float(np.sqrt(np.sum(w[sel] ** 2))) / va,
                             z_mean=float(z[sel].mean())))
    return pd.DataFrame(rows)


def cumulative_density(qlf: pd.DataFrame, m_cut: float) -> pd.DataFrame:
    """Number density of quasars brighter than m_cut in each redshift bin: sum of Phi dM.

    Only bins entirely brighter than m_cut are summed (choose m_cut on a bin edge).
    The error adds the bins' Poisson errors in quadrature.
    """
    rows = []
    for (zlo, zhi), g in qlf.groupby(["z_lo", "z_hi"]):
        g = g[g["m_hi"] <= m_cut + 1e-6]
        if g.empty:
            continue
        dm = (g["m_hi"] - g["m_lo"]).to_numpy()
        rows.append(dict(z_lo=zlo, z_hi=zhi, z_mid=0.5 * (zlo + zhi),
                         z_mean=float(np.average(g["z_mean"], weights=g["n"])),
                         density=float(np.sum(g["phi"] * dm)),
                         err=float(np.sqrt(np.sum((g["phi_err"] * dm) ** 2))),
                         n=int(g["n"].sum()), min_fill=float(g["fill"].min())))
    return pd.DataFrame(rows)


def complete_redshift_range(survey: Survey, m_cut: float, m_brightest: float, zgrid=None):
    """Redshifts where every M between m_brightest and m_cut is inside the flux limits."""
    zgrid = np.linspace(0.1, 6.0, 5901) if zgrid is None else zgrid
    ok = (survey.m_faint(zgrid) >= m_cut) & (survey.m_bright(zgrid) <= m_brightest)
    return zgrid[ok]
