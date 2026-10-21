"""Sirius B in Gaia DR3: what a straight-line astrometric model makes of a curving orbit.

Sirius A is far too bright for Gaia, but Sirius B (G = 8.52, 10.6" away in
2016) is in DR3 with a parallax of 374.49 +- 0.23 mas, 4.4 mas below the
Hipparcos-based value Bond et al. adopted (378.9 +- 1.4). Because the mass
scales as parallax^-3, the choice matters at the 3.5% level.

DR3's five-parameter model assumes uniform motion, but during the 34 months
of DR3 data Sirius B was accelerating along its orbit at ~29 mas/yr^2; the
curvature alone departs from a straight line by up to 19 mas. Gaia observed
the star in only 26 field-of-view transits, grouped into 10 visibility
periods, so a 5-parameter fit can trade some of the curvature against the
parallax and proper motion. This module rebuilds that situation:

1. Gaia's DR3 commanded scanning law (ESA, 10-s sampling) gives every
   transit of Sirius B's position across either field of view, with the scan
   direction; DPAC's DR3 astrometric data gaps are removed.
2. Sirius B's true path is simulated from our orbit and mass ratio (plus
   parallax and the barycentre's uniform motion), observed along-scan only,
   and fitted with Gaia's five-parameter model.
3. The difference between the fitted and the true parallax (and proper
   motion) is the bias. It is compared with the DR3 - Hipparcos difference.

This is a forward model with the nominal scanning law, not DPAC's pipeline:
the real solution had per-CCD weights, outlier rejection and calibration
terms we cannot reproduce. Its purpose is to show the size and sign of the
effect.
"""

from __future__ import annotations

import gzip

import numpy as np
import pandas as pd

import config as C
import orbits

SCANLAW = C.GAIA_DIR / "CommandedScanLaw_001.csv.gz"
GAPS = C.GAIA_DIR / "dr3_astrometry_gaps.csv"
#: jd_time in the scan law is TCB days since 2010-01-01T00:00 (JD 2455197.5).
JD0_SCANLAW = 2455197.5
#: Half-widths of a field of view: across-scan ~0.35 deg (7 CCD rows of 0.1
#: deg); along-scan, one 10-s sample sweeps 0.167 deg, so the transit sample is
#: the one with |AL| <= half of that.
AC_HALF_DEG = 0.35
AL_STEP_DEG = 60.0 * 10.0 / 3600.0
#: DR3 reference epoch and data span (Lindegren et al. 2021).
DR3_EPOCH = 2016.0


def transits(ra_deg, dec_deg, path=SCANLAW, chunksize=2_000_000):
    """Field-of-view transits of a fixed position in the DR3 commanded scanning law.

    Returns a DataFrame: time (Julian year, barycentric), fov (1 or 2), scan
    position angle psi (deg, direction of the along-scan motion, east of north),
    and the along/across-scan offsets of the sample used.
    """
    ra0, dec0 = np.radians(ra_deg), np.radians(dec_deg)
    v0 = np.array([np.cos(dec0) * np.cos(ra0), np.cos(dec0) * np.sin(ra0), np.sin(dec0)])
    rows = []
    cols = ["jd_time", "bjd_fov1", "bjd_fov2", "ra_fov1", "dec_fov1", "scan_angle_fov1",
            "ra_fov2", "dec_fov2", "scan_angle_fov2"]
    with gzip.open(path, "rt") as fh:
        for chunk in pd.read_csv(fh, usecols=cols, chunksize=chunksize):
            for k in (1, 2):
                ra = np.radians(chunk[f"ra_fov{k}"].to_numpy())
                dec = np.radians(chunk[f"dec_fov{k}"].to_numpy())
                cosd = np.sin(dec) * v0[2] + np.cos(dec) * (np.cos(ra) * v0[0] + np.sin(ra) * v0[1])
                near = cosd > np.cos(np.radians(1.0))
                if not near.any():
                    continue
                c = chunk[near]
                ra, dec = ra[near], dec[near]
                # gnomonic offsets of the target from the FOV centre (east, north), deg
                dra = ra0 - ra
                den = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(dra)
                xi = np.degrees(np.cos(dec0) * np.sin(dra) / den)
                eta = np.degrees((np.sin(dec0) * np.cos(dec) - np.cos(dec0) * np.sin(dec) * np.cos(dra)) / den)
                psi = np.radians(c[f"scan_angle_fov{k}"].to_numpy())
                al = xi * np.sin(psi) + eta * np.cos(psi)
                ac = -xi * np.cos(psi) + eta * np.sin(psi)
                hit = (np.abs(ac) < AC_HALF_DEG) & (np.abs(al) <= 0.5 * AL_STEP_DEG + 1e-9)
                for j in np.where(hit)[0]:
                    rows.append({"jd_time": float(c["jd_time"].iloc[j]), "bjd": float(c[f"bjd_fov{k}"].iloc[j]),
                                 "fov": k, "psi": float(np.degrees(psi[j]) % 360.0),
                                 "al": float(al[j]), "ac": float(ac[j])})
    t = pd.DataFrame(rows).sort_values("jd_time").reset_index(drop=True)
    # one sample per transit: drop duplicates closer than 1 hour within a FOV
    keep = np.ones(len(t), bool)
    for k in (1, 2):
        idx = np.where(t["fov"].to_numpy() == k)[0]
        tt = t["jd_time"].to_numpy()[idx]
        dup = np.r_[False, np.diff(tt) < 1.0 / 24.0]
        keep[idx[dup]] = False
    t = t[keep].reset_index(drop=True)
    t["year"] = orbits.jd_to_year(t["bjd"] + JD0_SCANLAW)
    return t


def remove_gaps(t: pd.DataFrame, gaps: pd.DataFrame) -> pd.DataFrame:
    """Drop transits inside DPAC's DR3 astrometric data gaps (OBMT revolutions)."""
    # The gap list is in OBMT revolutions; convert transit TCB to OBMT with the
    # relation used by gaiaunlimited (Lindegren et al. 2021, eq. 1 inverse):
    # tcb_jd2010 = 1717.6256 + (obmt_rev - 1717.6256 ... ) approximately linear.
    obmt = tcb_to_obmt(t["jd_time"].to_numpy())
    bad = np.zeros(len(t), bool)
    for s, e in zip(gaps["start"].to_numpy(), gaps["end"].to_numpy()):
        bad |= (obmt > s) & (obmt < e)
    return t[~bad].reset_index(drop=True)


def tcb_to_obmt(tcb_jd2010):
    """TCB (days since 2010.0) -> OBMT (6-h revolutions), Lindegren et al. 2021 Eq. (1) inverted.

    OBMT runs at 4 revolutions per day, offset so that rev 1717.6256 corresponds
    to TCB 2015.0 (JD 2457023.5); the small periodic term (< 2 ms) is ignored.
    """
    return 1717.6256 + 4.0 * (np.asarray(tcb_jd2010) - (2457023.5 - JD0_SCANLAW))


def observer_xyz(year):
    """Gaia's barycentric position (au, ICRS): the Earth's, pushed 1% outward to L2."""
    from astropy.coordinates import get_body_barycentric
    from astropy.time import Time
    tm = Time(orbits.year_to_jd(year), format="jd", scale="tdb")
    e = get_body_barycentric("earth", tm)
    s = get_body_barycentric("sun", tm)
    xyz = np.vstack([e.x.to_value("au"), e.y.to_value("au"), e.z.to_value("au")])
    sun = np.vstack([s.x.to_value("au"), s.y.to_value("au"), s.z.to_value("au")])
    return sun + (xyz - sun) * 1.0100


def parallax_factors(year, ra_deg, dec_deg):
    """(f_ra*, f_dec): displacement per unit parallax, east and north."""
    X, Y, Z = observer_xyz(year)
    a, d = np.radians(ra_deg), np.radians(dec_deg)
    f_ra = X * np.sin(a) - Y * np.cos(a)
    f_de = X * np.cos(a) * np.sin(d) + Y * np.sin(a) * np.sin(d) - Z * np.cos(d)
    return f_ra, f_de


def design(t, psi_deg, f_ra, f_de, epoch=DR3_EPOCH):
    """Along-scan design matrix of the 5-parameter model (dra*, ddec, plx, pmra*, pmdec)."""
    s, c = np.sin(np.radians(psi_deg)), np.cos(np.radians(psi_deg))
    dt = t - epoch
    return np.column_stack([s, c, f_ra * s + f_de * c, dt * s, dt * c])


def simulate_fit(tr: pd.DataFrame, ra_deg, dec_deg, orbit_theta, frac_B, plx_true_mas,
                 sigma_al_mas=1.0, n_ccd=9, rng=None, noise=False):
    """Fit Gaia's 5-parameter model to Sirius B's simulated along-scan positions.

    orbit_theta: relative orbit (P, T, e, a, i, Omega, omega); frac_B = M_A/M,
    the share of the relative orbit that is B's motion about the barycentre.
    Returns the fitted parameters, the 'truth' they should be compared with
    (parallax; B's velocity at the reference epoch), residuals and RUWE-like
    statistics.
    """
    t = tr["year"].to_numpy()
    psi = tr["psi"].to_numpy()
    if "f_ra" in tr:
        f_ra, f_de = tr["f_ra"].to_numpy(), tr["f_de"].to_numpy()
    else:
        f_ra, f_de = parallax_factors(t, ra_deg, dec_deg)
    xn, ye = orbits.relative_xy(t, *orbit_theta)                # north, east (arcsec)
    b_e, b_n = frac_B * ye * 1e3, frac_B * xn * 1e3              # B about the barycentre, mas
    s, c = np.sin(np.radians(psi)), np.cos(np.radians(psi))
    w = (b_e * s + b_n * c) + plx_true_mas * (f_ra * s + f_de * c)
    if noise:
        rng = rng or np.random.default_rng(C.SEED)
        w = w + rng.normal(0, sigma_al_mas / np.sqrt(n_ccd), len(w))
    A = design(t, psi, f_ra, f_de)
    coef, *_ = np.linalg.lstsq(A, w, rcond=None)
    res = w - A @ coef
    vn, ve = orbits.relative_velocity(DR3_EPOCH, *orbit_theta)
    xn0, ye0 = orbits.relative_xy(DR3_EPOCH, *orbit_theta)
    truth = {"plx": plx_true_mas, "pmra_orb": frac_B * ve * 1e3, "pmdec_orb": frac_B * vn * 1e3,
             "dra_orb": frac_B * ye0 * 1e3, "ddec_orb": frac_B * xn0 * 1e3}
    cov = np.linalg.inv(A.T @ A) * (sigma_al_mas**2 / n_ccd)
    return {"coef": coef, "truth": truth, "res": res, "rms_res": float(np.sqrt(np.mean(res**2))),
            "plx_bias": float(coef[2] - plx_true_mas),
            "pmra_bias": float(coef[3] - truth["pmra_orb"]), "pmdec_bias": float(coef[4] - truth["pmdec_orb"]),
            "formal_err": np.sqrt(np.diag(cov)), "n": len(t)}


# ---------------------------------------------------------------------------
# What DR3 should have measured, from Hipparcos + our orbit
# ---------------------------------------------------------------------------

#: Hipparcos observing window (Perryman et al. 1997): 1989.85 - 1993.21.
HIP_WINDOW = (1989.85, 1993.21)


def hipparcos_barycentre(orbit_theta, aA, hip=None, dmso=None):
    """The barycentre's proper motion and parallax implied by Hipparcos.

    Hipparcos solved Sirius with Sirius A's photocentre orbit fixed at Gatewood
    & Gatewood's (1978) elements (DMSA/O), so its proper motion refers to the
    barycentre *if* that orbit was right. Any difference between GG78's orbit
    of A and the true one leaks, over the 3.4-yr window, into the fitted
    proper motion; we remove that leak with our orbit and Bond et al.'s a_A.
    Returns (pmra*, pmdec) in mas/yr at J1991.25 and the leak itself.
    """
    import literature as L
    hip = hip or L.HIP2_SIRIUS
    d = dmso or L.HIP1_DMSA_O_SIRIUS
    t = np.linspace(*HIP_WINDOW, 400)
    P_gg = d["P_d"] / C.JULIAN_YEAR
    T_gg = orbits.jd_to_year(d["T_jd2440000"] + 2440000.0)
    # A about the barycentre: GG78 (photocentre elements, omega of A) vs ours (-a_A/a x relative)
    xa_gg, ya_gg = orbits.relative_xy(t, P_gg, T_gg, d["e"], d["a0_mas"] / 1e3, d["i"], d["Omega"], d["omega"])
    xr, yr = orbits.relative_xy(t, *orbit_theta)
    k = aA / orbit_theta[3]
    dn, de = (-k * xr - xa_gg) * 1e3, (-k * yr - ya_gg) * 1e3     # mas, north / east
    slope_n = np.polyfit(t - 1991.25, dn, 1)[0]
    slope_e = np.polyfit(t - 1991.25, de, 1)[0]
    pmra = hip["pmra"][0] - slope_e
    pmdec = hip["pmdec"][0] - slope_n
    return {"pmra": pmra, "pmdec": pmdec, "leak_pmra": slope_e, "leak_pmdec": slope_n,
            "pmra_err": hip["pmra"][1], "pmdec_err": hip["pmdec"][1]}


def perspective_pm_change(pmra, pmdec, plx_mas, rv_kms, dt_yr):
    """Change of proper motion from perspective acceleration over dt (mas/yr).

    d(mu)/dt = -2 mu v_r / d: an approaching star's proper motion grows.
    """
    rate = -2.0 * rv_kms * (plx_mas / 1e3) / (C.AU_PER_YR_KMS * 206264.806)   # per year
    return pmra * rate * dt_yr, pmdec * rate * dt_yr


def expected_dr3(orbit_theta, aA, plx_true_mas, rv_kms=-8.47):
    """Sirius B's true parallax and mean proper motion at J2016.0 (what an ideal DR3 would give)."""
    bary = hipparcos_barycentre(orbit_theta, aA)
    dpa, dpd = perspective_pm_change(bary["pmra"], bary["pmdec"], plx_true_mas, rv_kms, DR3_EPOCH - 1991.25)
    vn, ve = orbits.relative_velocity(DR3_EPOCH, *orbit_theta)
    fA = 1.0 - aA / orbit_theta[3]            # B's share of the relative motion = M_A / M
    return {"plx": plx_true_mas, "pmra": bary["pmra"] + dpa + fA * ve * 1e3,
            "pmdec": bary["pmdec"] + dpd + fA * vn * 1e3, "bary": bary, "persp": (dpa, dpd),
            "orbital_pm": (fA * ve * 1e3, fA * vn * 1e3)}


def visibility_periods(tr: pd.DataFrame, gap_days=4.0) -> np.ndarray:
    """Group transits into visibility periods (gaps > 4 days, as in Gaia DR3)."""
    t = tr["year"].to_numpy()
    return np.r_[0, np.cumsum(np.diff(t) * C.JULIAN_YEAR > gap_days)]


def subset_monte_carlo(tr, ra_deg, dec_deg, orbit_theta, frac_B, plx_true_mas, n_vp=10, n_transits=26,
                       n_draw=20000, seed=C.SEED):
    """Biases of the 5-parameter fit for random choices of which visibility periods DR3 used.

    DR3 used 26 transits in 10 visibility periods for Sirius B (matched_transits,
    visibility_periods_used); the scanning law offers 55 transits in 19 periods,
    one of which (June 2016) holds 24 consecutive transits. Each draw picks 10
    periods at random and keeps their transits, trimming the largest period at
    random so that the total does not exceed 26.
    """
    rng = np.random.default_rng(seed)
    tr = tr.copy()
    tr["f_ra"], tr["f_de"] = parallax_factors(tr["year"].to_numpy(), ra_deg, dec_deg)
    vp = visibility_periods(tr)
    groups = [np.where(vp == g)[0] for g in range(vp.max() + 1)]
    out = np.full((n_draw, 4), np.nan)
    for k in range(n_draw):
        pick = rng.choice(len(groups), n_vp, replace=False)
        idx = [groups[g] for g in pick]
        sizes = np.array([len(i) for i in idx])
        excess = sizes.sum() - n_transits
        if excess > 0:                       # trim the biggest period
            j = int(np.argmax(sizes))
            keep = max(1, sizes[j] - excess)
            idx[j] = np.sort(rng.choice(idx[j], keep, replace=False))
        sel = np.concatenate(idx)
        r = simulate_fit(tr.iloc[sel], ra_deg, dec_deg, orbit_theta, frac_B, plx_true_mas)
        out[k] = (r["plx_bias"], r["pmra_bias"], r["pmdec_bias"], r["rms_res"])
    return pd.DataFrame(out, columns=["plx_bias", "pmra_bias", "pmdec_bias", "rms_res"])


def hipparcos_leak_mc(orbit_theta, aA, n_epochs=47, n_draw=2000, seed=C.SEED):
    """How GG78's outdated orbit of Sirius A leaks into the Hipparcos parallax and proper motion.

    Hipparcos fitted position, parallax and proper motion after subtracting
    GG78's photocentre orbit. The residual orbit (ours minus GG78's) is
    'observed' along random scan directions at n_epochs random times in the
    mission window (Sirius had 47 field transits in the new reduction) and
    fitted with the 5-parameter model; the fitted parallax and proper motion are
    the leak. Returns a DataFrame of leaks (mas, mas/yr) over random scan
    histories: the real scan history is not reproduced, only its spread.
    """
    import literature as L
    d = L.HIP1_DMSA_O_SIRIUS
    rng = np.random.default_rng(seed)
    P_gg = d["P_d"] / C.JULIAN_YEAR
    T_gg = orbits.jd_to_year(d["T_jd2440000"] + 2440000.0)
    k = aA / orbit_theta[3]
    ra, dec = L.HIP2_SIRIUS["ra"], L.HIP2_SIRIUS["dec"]
    rows = []
    for _ in range(n_draw):
        t = np.sort(rng.uniform(*HIP_WINDOW, n_epochs))
        psi = rng.uniform(0, 360, n_epochs)
        xa_gg, ya_gg = orbits.relative_xy(t, P_gg, T_gg, d["e"], d["a0_mas"] / 1e3, d["i"], d["Omega"], d["omega"])
        xr, yr = orbits.relative_xy(t, *orbit_theta)
        dn, de = (-k * xr - xa_gg) * 1e3, (-k * yr - ya_gg) * 1e3
        s, c = np.sin(np.radians(psi)), np.cos(np.radians(psi))
        w = de * s + dn * c
        f_ra, f_de = parallax_factors(t, ra, dec)
        A = design(t, psi, f_ra, f_de, epoch=1991.25)
        coef, *_ = np.linalg.lstsq(A, w, rcond=None)
        rows.append(coef[2:])
    return pd.DataFrame(rows, columns=["plx_leak", "pmra_leak", "pmdec_leak"])
