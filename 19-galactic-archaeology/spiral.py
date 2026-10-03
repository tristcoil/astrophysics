"""Part 3: the phase spiral -- the disc is still ringing.

Every disc star bobs up and down through the Galactic plane as it orbits. Plot
height z against vertical velocity v_z and each star traces a loop: fast
through the plane, stopping at the top and bottom of its swing. Stars with
bigger swings take longer (the disc's pull weakens with height), so if some
event once pushed a group of stars in the same direction at the same time,
the group shears into a spiral -- a "snail shell" in the (z, v_z) plane that
winds tighter with time (Antoja et al. 2018, from Gaia DR2).

Two things are measured here:

1. **The spiral itself** in the ~2.9 million DR3 stars with full 3D motions in a
   ring 0.2 kpc wide around the Sun's Galactocentric radius: counts, and the
   median azimuthal and radial velocity across the (z, v_z) plane.
2. **Its age.** In "angle-frequency" coordinates -- where each star is on its
   vertical loop (the angle theta) and how fast it goes round (the frequency
   Omega) -- a phase spiral is a set of straight stripes, theta = theta0 +
   Omega t. The slope t is the time since the disturbance. We find the t at
   which the pattern lines up best.

The frequencies come from a model of the Galaxy's vertical gravity (McMillan
2017), so the age is model-dependent: other published potentials are tried as a
systematic check.
"""

from __future__ import annotations

import numpy as np

import galaxy as gx
from config import R0_KPC, RING_HALF_WIDTH_KPC, RUWE_MAX
from fetch_data import RING_PATH
from gd1 import _read

#: kpc / (km/s) in Myr.
KPC_PER_KMS_MYR = 977.792221
#: Orbits used for dating: vertical amplitudes where the sample is well
#: populated on every part of the loop (Gaia's parallax cut thins out stars
#: far from the Sun, i.e. at large |z|).
DATING_ZMAX_KPC = (0.15, 1.0)
#: Trial times, Myr.
T_GRID_MYR = np.arange(0.0, 1500.1, 2.0)
#: Frequency bins for the phase-stacking estimator.
N_OMEGA_BINS = 40


def load_ring() -> dict:
    d = _read(RING_PATH)
    dist = 1.0 / d["parallax"]
    g = gx.galactocentric(d["ra"], d["dec"], dist, d["pmra"], d["pmdec"], d["radial_velocity"])
    d.update(distance=dist, R=g.R, z=g.z, v_R=g.v_R, v_phi=g.v_phi, v_z=g.vz, L_z=g.L_z,
             phi=g.phi_deg)
    d["in_ring"] = np.abs(g.R - R0_KPC) < RING_HALF_WIDTH_KPC
    d["good"] = d["ruwe"] < RUWE_MAX
    return d


# --------------------------------------------------------------------------
# Maps of the (z, v_z) plane
# --------------------------------------------------------------------------


def zvz_maps(z, vz, v_phi, v_R, *, z_lim=1.2, vz_lim=60.0, dz=0.02, dvz=1.0, min_count=5):
    """Counts and median v_phi, v_R in (z, v_z) bins (Antoja et al. 2018 Fig. 1)."""
    from scipy.stats import binned_statistic_2d

    ez = np.arange(-z_lim, z_lim + dz / 2, dz)
    ev = np.arange(-vz_lim, vz_lim + dvz / 2, dvz)
    counts, _, _ = np.histogram2d(z, vz, bins=[ez, ev])
    med_phi = binned_statistic_2d(z, vz, v_phi, "median", bins=[ez, ev]).statistic
    med_R = binned_statistic_2d(z, vz, v_R, "median", bins=[ez, ev]).statistic
    med_phi[counts < min_count] = np.nan
    med_R[counts < min_count] = np.nan
    return ez, ev, counts, med_phi, med_R


def relative_contrast(counts, smooth_bins=(6.0, 6.0)):
    """Counts divided by a smoothed version of themselves, minus one: small-
    scale structure (the spiral) stands out from the smooth overall shape."""
    from scipy.ndimage import gaussian_filter

    smooth = gaussian_filter(counts.astype(float), smooth_bins, mode="nearest")
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(smooth > 3, counts / smooth - 1.0, np.nan)


# --------------------------------------------------------------------------
# Vertical angle-frequency coordinates
# --------------------------------------------------------------------------


class VerticalOrbits:
    """Oscillation in the vertical potential Phi(R0, z) - Phi(R0, 0).

    The standard local approximation (Antoja et al. 2018 use it too): each star
    oscillates in the vertical potential at the Sun's radius, independent of
    its in-plane motion."""

    def __init__(self, pot, R=R0_KPC, z_max=3.0, n_z=3000, n_amp=400, n_psi=400):
        self.pot = pot
        zz = np.linspace(0, z_max, n_z)
        self._z = zz
        self._phi = pot.vertical_profile(R, zz)
        # Tabulate the time from z = 0 to height z on the way up, for a grid of
        # amplitudes, with z = A sin(psi) so the integrand stays finite.
        amps = np.geomspace(0.002, z_max * 0.98, n_amp)
        psi = np.linspace(0, np.pi / 2, n_psi)
        A, P = np.meshgrid(amps, psi, indexing="ij")
        Ez = np.interp(amps, zz, self._phi)[:, None]
        with np.errstate(invalid="ignore", divide="ignore"):
            ke = 2 * (Ez - np.interp(A * np.sin(P), zz, self._phi))
            integrand = A * np.cos(P) / np.sqrt(np.maximum(ke, 0))
        # At psi = pi/2 the ratio has a finite limit; take it from the harmonic
        # form near the turning point (evaluated one step earlier).
        integrand[:, -1] = integrand[:, -2]
        tau = np.concatenate([np.zeros((n_amp, 1)),
                              np.cumsum(0.5 * (integrand[:, 1:] + integrand[:, :-1])
                                        * np.diff(psi)[None, :], axis=1)], axis=1)
        self._amps, self._psi, self._tau = amps, psi, tau
        self._quarter = tau[:, -1]

    def vertical_energy(self, amp):
        """Phi(amp) - Phi(0): the vertical energy of an orbit of this amplitude."""
        return np.interp(amp, self._z, self._phi)

    def amplitude(self, z, vz):
        Ez = 0.5 * vz**2 + np.interp(np.abs(z), self._z, self._phi)
        return np.interp(Ez, self._phi, self._z)

    def frequency(self, amp):
        """Vertical angular frequency, rad/Myr."""
        quarter = np.interp(amp, self._amps, self._quarter)
        return 2 * np.pi / (4 * quarter * KPC_PER_KMS_MYR)

    def angle_frequency(self, z, vz):
        """(theta, Omega, amplitude): theta in [0, 2 pi), 0 when crossing the
        plane upwards, pi/2 at the top of the swing."""
        from scipy.interpolate import RegularGridInterpolator

        amp = np.clip(self.amplitude(z, vz), self._amps[0], self._amps[-1])
        psi = np.arcsin(np.clip(np.abs(z) / amp, 0, 1))
        interp = RegularGridInterpolator((np.log(self._amps), self._psi), self._tau,
                                         bounds_error=False, fill_value=None)
        tau = interp(np.column_stack([np.log(amp), psi]))
        quarter = np.interp(amp, self._amps, self._quarter)
        frac = tau / quarter * (np.pi / 2)  # angle from the upward crossing
        up = vz >= 0
        theta = np.where(z >= 0, np.where(up, frac, np.pi - frac),
                         np.where(up, 2 * np.pi - frac, np.pi + frac))
        return np.mod(theta, 2 * np.pi), self.frequency(amp), amp

    def harmonic_frequency(self):
        """Small-amplitude limit sqrt(d2Phi/dz2) at z = 0, rad/Myr: fit
        Phi = a z^2 + b |z|^3 over the first 20 pc (a disc whose density is
        exponential in |z| has exactly this form near the midplane)."""
        s = (self._z > 0) & (self._z <= 0.02)
        A = np.column_stack([self._z[s] ** 2, self._z[s] ** 3])
        a, _ = np.linalg.lstsq(A, self._phi[s], rcond=None)[0]
        return np.sqrt(2 * a) / KPC_PER_KMS_MYR


def phase_stack(theta, omega, weights, m=1, t_grid=T_GRID_MYR, n_bins=N_OMEGA_BINS):
    """Amplitude of an m-armed winding spiral as a function of trial age.

    Stars are grouped in narrow bins of frequency. In each bin the complex
    amplitude c_k = <w exp(i m theta)> measures the lopsidedness of the bin's
    stars around their loop. A spiral that has wound for time t has
    arg(c_k) = const + m Omega_k t; so rotating each bin back by m Omega_k t
    and adding them up gives the largest total at the true t."""
    edges = np.quantile(omega, np.linspace(0, 1, n_bins + 1))
    k = np.clip(np.digitize(omega, edges) - 1, 0, n_bins - 1)
    n = np.maximum(np.bincount(k, minlength=n_bins), 1)
    c = (np.bincount(k, weights * np.cos(m * theta), minlength=n_bins)
         + 1j * np.bincount(k, weights * np.sin(m * theta), minlength=n_bins)) / n
    om = np.bincount(k, omega, minlength=n_bins) / n
    amp = np.abs(np.exp(-1j * m * np.outer(t_grid, om)) @ c) / n_bins
    return amp, om, c


def date_spiral(theta, omega, weights, m=1, n_boot=100, seed=2):
    rng = np.random.default_rng(seed)
    amp, om, c = phase_stack(theta, omega, weights, m)
    # Skip the region near t = 0, where any non-winding pattern (e.g. the
    # sample's own selection in height) adds up coherently.
    min_t = 150.0
    ok = T_GRID_MYR >= min_t
    t_best = float(T_GRID_MYR[ok][np.argmax(amp[ok])])
    boots = []
    for _ in range(n_boot):
        i = rng.integers(0, len(theta), len(theta))
        a, _, _ = phase_stack(theta[i], omega[i], weights[i], m)
        boots.append(T_GRID_MYR[ok][np.argmax(a[ok])])
    boots = np.array(boots)
    return dict(t_myr=t_best, t_lo=float(np.percentile(boots, 16)),
                t_hi=float(np.percentile(boots, 84)), amplitude=amp.tolist(),
                bin_omega=om.tolist(), bin_phase=np.angle(c).tolist(), m=m)


#: The baseline amplitude range for dating, and alternatives that probe how
#: much the answer depends on that choice. Below ~0.3 kpc a small error in the
#: Sun's assumed height or vertical velocity (tens of pc, ~1 km/s) makes a
#: lopsided pattern that does not wind and swamps the spiral's phase.
DATING_RANGES_KPC = ((0.3, 1.0), (0.25, 0.9), (0.3, 0.8), (0.35, 1.1), (0.4, 1.2), (0.15, 1.0))
#: The local disc's surface density within 1.1 kpc of the plane, measured from
#: the vertical motions of stars (Bovy & Rix 2013: 68 +/- 4 Msun/pc^2).
SIGMA_11_MEASURED = (68.0, 4.0)
N_AMP_BINS = 80


def phase_mixed_residuals(theta, amp, *values, n_amp=N_AMP_BINS, n_theta=72):
    """What a fully phase-mixed disc would look like, subtracted star by star.

    Stars are grouped by vertical amplitude. Within a group, a phase-mixed
    population is spread evenly around its loop (uniform in angle) and every
    velocity has the same average at every angle. So the star-count contrast
    is the count in each (amplitude, angle) cell relative to the group's mean,
    and each velocity's residual is the value minus the group's mean. What is
    left is the spiral -- plus whatever the sample's selection does to the
    angle distribution (fewer stars far from the plane)."""
    ab = np.quantile(amp, np.linspace(0, 1, n_amp + 1))
    ia = np.clip(np.digitize(amp, ab) - 1, 0, n_amp - 1)
    tb = np.linspace(0, 2 * np.pi, n_theta + 1)
    it = np.clip(np.digitize(theta, tb) - 1, 0, n_theta - 1)
    H = np.bincount(ia * n_theta + it, minlength=n_amp * n_theta).reshape(n_amp, n_theta)
    contrast = (H / H.mean(axis=1, keepdims=True) - 1)[ia, it]
    n = np.bincount(ia, minlength=n_amp)
    resid = [v - (np.bincount(ia, v, minlength=n_amp) / n)[ia] for v in values]
    return contrast, resid


def binned_map(z, vz, value, ez, ev, min_count=20):
    from scipy.stats import binned_statistic_2d

    st = binned_statistic_2d(z, vz, value, "mean", bins=[ez, ev]).statistic
    n = binned_statistic_2d(z, vz, None, "count", bins=[ez, ev]).statistic
    return np.where(n >= min_count, st, np.nan)


def surface_density_11(pot, R=R0_KPC):
    """Sigma(|z| < 1.1 kpc) = K_z(1.1 kpc) / (2 pi G), Msun/pc^2."""
    from config import G_KPC_KMS2_MSUN

    zz = np.array([1.099, 1.101])
    kz = np.diff(pot.vertical_profile(R, zz))[0] / 0.002
    return float(kz / (2 * np.pi * G_KPC_KMS2_MSUN) / 1e6)


def analyse(n_boot: int = 100) -> tuple[dict, dict]:
    d = load_ring()
    s = d["in_ring"] & d["good"]
    ref = dict(n_fetched=int(len(d["source_id"])), n_ring=int(np.sum(d["in_ring"])),
               n_ring_good=int(s.sum()))
    z, vz, vphi, vR = d["z"][s], d["v_z"][s], d["v_phi"][s], d["v_R"][s]
    res = dict(ref)
    res["median_errors_kms"] = dict(radial_velocity=float(np.median(d["radial_velocity_error"][s])))

    pots = {"McMillan (2017)": gx.mcmillan17(), "Helmi et al. (2017)": gx.AnalyticHelmi17()}
    res["surface_density_11"] = {k: surface_density_11(p) for k, p in pots.items()}
    res["surface_density_11_measured"] = list(SIGMA_11_MEASURED)

    # Angle-frequency coordinates in the baseline potential.
    vo = VerticalOrbits(pots["McMillan (2017)"])
    theta, omega, amp = vo.angle_frequency(z, vz)
    contrast, (dR, dphi) = phase_mixed_residuals(theta, amp, vR, vphi)
    ez = np.arange(-1.0, 1.0001, 0.025)
    ev = np.arange(-50.0, 50.001, 1.25)
    maps = dict(contrast=binned_map(z, vz, contrast, ez, ev),
                dR=binned_map(z, vz, dR, ez, ev), dphi=binned_map(z, vz, dphi, ez, ev),
                counts=np.histogram2d(z, vz, bins=[ez, ev])[0])
    res["harmonic_period_myr"] = float(2 * np.pi / vo.harmonic_frequency())

    # Dating: three tracers x several amplitude ranges (systematics), plus
    # bootstrap errors for the baseline range.
    tracers = {"counts": None, "v_phi": dphi, "v_R": dR}
    table = []
    for lo, hi in DATING_RANGES_KPC:
        u = (amp > lo) & (amp < hi)
        for name, w in tracers.items():
            weights = np.ones(u.sum()) if w is None else w[u]
            nb = n_boot if (lo, hi) == DATING_RANGES_KPC[0] else 0
            r = date_spiral(theta[u], omega[u], weights, 1, n_boot=nb) if nb else None
            if r is None:
                a, _, _ = phase_stack(theta[u], omega[u], weights, 1)
                ok = T_GRID_MYR >= 150.0
                r = dict(t_myr=float(T_GRID_MYR[ok][np.argmax(a[ok])]))
            table.append(dict(amp_lo=lo, amp_hi=hi, tracer=name, **{
                k: v for k, v in r.items() if k in ("t_myr", "t_lo", "t_hi")}))
    res["dating_table"] = table
    ages = np.array([row["t_myr"] for row in table])
    base = [row for row in table if (row["amp_lo"], row["amp_hi"]) == DATING_RANGES_KPC[0]]
    res["age_summary"] = dict(median_myr=float(np.median(ages)), min_myr=float(ages.min()),
                              max_myr=float(ages.max()),
                              baseline={row["tracer"]: row for row in base})
    # Full stacking curves for the baseline range (figures).
    u = (amp > DATING_RANGES_KPC[0][0]) & (amp < DATING_RANGES_KPC[0][1])
    curves = {}
    for name, w in tracers.items():
        weights = np.ones(u.sum()) if w is None else w[u]
        a, om, c = phase_stack(theta[u], omega[u], weights, 1)
        curves[name] = dict(amplitude=a.tolist(), bin_omega=om.tolist(),
                            bin_phase=np.angle(c).tolist())
    res["baseline_curves"] = curves
    # The same in the other potential (baseline range, counts and v_phi).
    vo_h = VerticalOrbits(pots["Helmi et al. (2017)"])
    th_h, om_h, amp_h = vo_h.angle_frequency(z, vz)
    _, (dR_h, dphi_h) = phase_mixed_residuals(th_h, amp_h, vR, vphi)
    uh = (amp_h > DATING_RANGES_KPC[0][0]) & (amp_h < DATING_RANGES_KPC[0][1])
    alt = {}
    for name, w in (("counts", None), ("v_phi", dphi_h), ("v_R", dR_h)):
        weights = np.ones(uh.sum()) if w is None else w[uh]
        a, _, _ = phase_stack(th_h[uh], om_h[uh], weights, 1)
        ok = T_GRID_MYR >= 150.0
        alt[name] = float(T_GRID_MYR[ok][np.argmax(a[ok])])
    res["age_other_potential"] = {"Helmi et al. (2017)": alt}

    artefacts = dict(d=d, s=s, ez=ez, ev=ev, maps=maps, theta=theta, omega=omega, amp=amp,
                     use=u, dphi=dphi, dR=dR, contrast=contrast)
    return res, artefacts
