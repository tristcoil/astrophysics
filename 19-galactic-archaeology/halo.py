"""Part 2: Gaia-Sausage-Enceladus, the galaxy the Milky Way ate.

About 10 billion years ago a dwarf galaxy of a few billion suns fell into the
young Milky Way on a nearly radial orbit and was shredded. Its stars still
move on the elongated orbits it arrived on. Near the Sun they are a few per
cent of the stars at most, but they dominate the *stellar halo*, and they show
up two ways:

* **In velocity** (Belokurov et al. 2018): halo stars near the Sun do not
  orbit the Galaxy like the disc does. Their radial velocities (towards and
  away from the Galactic Centre) are huge while their rotation is ~zero: in the
  plane of radial versus azimuthal velocity they form a long thin "sausage".
* **In orbital energy and angular momentum** (Helmi et al. 2018): stars from
  one progenitor keep similar integrals of motion long after they have spread
  around the sky. GSE is a large blob at L_z ~ 0.

The trap, as the atlas warns: the usual way to pick halo stars is by velocity
(|v - v_LSR| > 210 km/s). That cut shapes what we see. So the same question is
asked twice: once with the velocity cut, as the discovery papers did, and once
with Gaia's own spectroscopic metallicities and *no* cut on velocity at all.
"""

from __future__ import annotations

import numpy as np

import galaxy as gx
import photometry as ph
from config import (
    FEH_BINS, GSPSPEC_ALPHA_CALIB, GSPSPEC_ALPHA_LOGG_RANGE, GSPSPEC_MEDIUM_MAX,
    GSPSPEC_MH_CALIB, HALO_MAX_DISTANCE_KPC, HALO_SPEED_KMS, N_DR3_SOURCES,
    RANDOM_FRACTION, REFERENCE, RUWE_MAX,
)
from fetch_data import GSPSPEC_PATH, HALO_PATH
from gd1 import _read

#: The chemistry-only samples use stars within this distance (Belokurov et al.
#: 2018 reach ~10 kpc with SDSS dwarfs; Gaia's RVS stars are bright and near).
CHEM_MAX_DISTANCE_KPC = 3.0
#: GSP-Spec metallicities are trusted only for cool giants. Hot stars have
#: few metal lines in the RVS window: among "metal-poor" ([M/H] < -1.66)
#: medium-quality stars, half sit on cold thin-disc orbits, and those have a
#: median Teff of 7,500 K -- disc A/F stars with metallicities that are
#: simply wrong. Giants in this range put only 6% there, as a real metal-weak
#: disc would.
CHEM_TEFF_RANGE = (4000.0, 5500.0)
CHEM_LOGG_MAX = 3.5
#: Parallaxes good to 10%: at |v| ~ 300 km/s a 20% distance error would add
#: ~60 km/s of scatter to the tangential velocity, inflating the dispersions.
CHEM_PARALLAX_OVER_ERROR = 10.0
N_BOOT = 40
#: Overlapping windows for the finer look (and the interactive page).
FINE_WINDOW_HALF = 0.15
FINE_CENTRES = np.round(np.arange(-2.45, -0.849, 0.1), 2)
#: Coarser metallicity bins shown in the sausage figure: below Belokurov et
#: al.'s transition, GSE's own metallicity, and where the thick disc takes over.
DISPLAY_BINS = ((-3.0, -1.66), (-1.66, -1.0), (-1.0, -0.7))
#: Metallicity windows for comparing [alpha/Fe] at fixed [M/H].
ALPHA_BINS = ((-1.6, -1.2), (-1.2, -0.9), (-0.9, -0.6))


def calibrated_mh(d):
    logg = d["logg_gspspec"]
    return d["mh_gspspec"] + sum(p * logg**k for k, p in enumerate(GSPSPEC_MH_CALIB))


def calibrated_alpha(d):
    logg = d["logg_gspspec"]
    a = d["alphafe_gspspec"] + sum(p * logg**k for k, p in enumerate(GSPSPEC_ALPHA_CALIB))
    lo, hi = GSPSPEC_ALPHA_LOGG_RANGE
    return np.where((logg > lo) & (logg < hi), a, np.nan)


def gspspec_medium(d):
    """Recio-Blanco et al. (2023) medium-quality sample (first 13 flags)."""
    flags = d["flags_gspspec"]
    limit = np.array(GSPSPEC_MEDIUM_MAX)
    out = np.zeros(len(flags), bool)
    for i, f in enumerate(flags):
        if len(f) >= 13 and f[:13].isdigit():
            out[i] = np.all(np.frombuffer(f[:13].encode(), np.uint8) - 48 <= limit)
    return out & np.isfinite(d["mh_gspspec"])


def enrich(d: dict, pot_h17, pot_m17) -> dict:
    """Add distances, Galactocentric kinematics, integrals of motion and
    calibrated chemistry to a star table (in place, returned)."""
    dist = 1.0 / d["parallax"]
    g = gx.galactocentric(d["ra"], d["dec"], dist, d["pmra"], d["pmdec"], d["radial_velocity"])
    U, V, W = gx.heliocentric_uvw(d["ra"], d["dec"], dist, d["pmra"], d["pmdec"],
                                  d["radial_velocity"])
    v_r, v_theta, _ = g.spherical_velocities()
    d.update(
        distance=dist, x=g.x, y=g.y, z=g.z, R=g.R, v_R=g.v_R, v_phi=g.v_phi, v_z=g.vz,
        v_r=v_r, v_theta=v_theta, L_z=g.L_z, L_perp=g.L_perp, U=U, V=V, W=W,
        speed_lsr=gx.speed_vs_lsr(U, V, W),
        E_h17=gx.energy(pot_h17, g), E_m17=gx.energy(pot_m17, g),
    )
    d["circ_m17"] = gx.circularity(pot_m17, d["E_m17"], d["L_z"])
    d["good"] = d["ruwe"] < RUWE_MAX
    if "mh_gspspec" in d:
        d["mh"] = calibrated_mh(d)
        d["alpha"] = calibrated_alpha(d)
        d["chem_ok"] = gspspec_medium(d)
    return d


#: No part of the analysis looks beyond this distance (E-L_z uses 5 kpc, as
#: Helmi et al. did), so farther stars are dropped on loading. They are more
#: than half of the fast file: distant inner-disc stars move fast relative to
#: the Sun's neighbourhood simply because the disc rotates differentially.
LOAD_MAX_DISTANCE_KPC = 5.0


def _read_near(path):
    d = _read(path)
    keep = (d["parallax"] > 1.0 / LOAD_MAX_DISTANCE_KPC)
    return {k: v[keep] for k, v in d.items()}, int(keep.size)


def load_samples():
    pot_h17, pot_m17 = gx.AnalyticHelmi17(), gx.mcmillan17()
    raw, n_file = _read_near(HALO_PATH)
    fast = enrich(raw, pot_h17, pot_m17)
    fast["_n_file"] = n_file
    fast["is_random"] = fast["random_index"] < int(round(RANDOM_FRACTION * N_DR3_SOURCES))
    mp = enrich(_read(GSPSPEC_PATH), pot_h17, pot_m17)
    return fast, mp, pot_h17, pot_m17


def trusted_giants(d):
    """Stars whose GSP-Spec metallicity we trust, with good velocities: cool
    giants, medium-quality flags, parallax better than 10%, within 3 kpc."""
    return (d["good"] & d["chem_ok"] & (d["distance"] < CHEM_MAX_DISTANCE_KPC)
            & (d["teff_gspspec"] > CHEM_TEFF_RANGE[0]) & (d["teff_gspspec"] < CHEM_TEFF_RANGE[1])
            & (d["logg_gspspec"] < CHEM_LOGG_MAX)
            & (d["parallax"] / d["parallax_error"] > CHEM_PARALLAX_OVER_ERROR))


# --------------------------------------------------------------------------
# Velocity ellipsoids
# --------------------------------------------------------------------------


def _log_gauss(X, mean, cov):
    d = X - mean
    L = np.linalg.cholesky(cov)
    z = np.linalg.solve(L, d.T)
    return -0.5 * np.sum(z**2, axis=0) - np.sum(np.log(np.diag(L))) - 1.5 * np.log(2 * np.pi)


def halo_mixture(v_r, v_theta, v_phi, n_disc: int, n_iter: int = 500, tol: float = 1e-8):
    """Gaussian mixture in (v_r, v_theta, v_phi) with one *halo* component and
    ``n_disc`` disc components, fitted by expectation-maximisation.

    The halo is constrained the way Belokurov et al. (2018) model it: no net
    radial or polar motion (mean v_r = v_theta = 0), free rotation, and
    independent dispersions along the three spherical axes -- so its shape
    is the velocity ellipsoid whose anisotropy we want. Disc components are
    unconstrained 3D Gaussians. Returns (components, log-likelihood, number
    of free parameters)."""
    X = np.column_stack([v_r, v_theta, v_phi]).astype(float)
    n = len(X)
    comps = [dict(kind="halo", w=0.5 if n_disc else 1.0, mean=np.array([0.0, 0.0, 0.0]),
                  cov=np.diag([150.0**2, 100.0**2, 100.0**2]))]
    for k in range(n_disc):
        comps.append(dict(kind="disc", w=0.5 / n_disc, mean=np.array([0.0, 0.0, 210.0 - 50 * k]),
                          cov=np.diag([(40.0 + 30 * k) ** 2, (25.0 + 20 * k) ** 2,
                                       (30.0 + 25 * k) ** 2])))
    prev = -np.inf
    for _ in range(n_iter):
        logp = np.array([np.log(c["w"]) + _log_gauss(X, c["mean"], c["cov"]) for c in comps])
        mx = logp.max(axis=0)
        ll = float(np.sum(mx + np.log(np.sum(np.exp(logp - mx), axis=0))))
        resp = np.exp(logp - mx - np.log(np.sum(np.exp(logp - mx), axis=0)))
        for c, r in zip(comps, resp):
            nk = r.sum()
            c["w"] = max(nk / n, 1e-6)
            if c["kind"] == "halo":
                mean = np.array([0.0, 0.0, np.sum(r * X[:, 2]) / nk])
                var = np.sum(r[:, None] * (X - mean) ** 2, axis=0) / nk
                c["mean"], c["cov"] = mean, np.diag(np.maximum(var, 25.0))
            else:
                mean = np.sum(r[:, None] * X, axis=0) / nk
                d = X - mean
                c["mean"] = mean
                c["cov"] = (r[:, None] * d).T @ d / nk + np.eye(3) * 4.0
        if ll - prev < tol * abs(ll):
            break
        prev = ll
    n_par = 4 + 1 + n_disc * 10  # halo: v_phi mean + 3 variances (+weight); disc: 3+6 (+weight)
    return comps, ll, n_par - 1


def velocity_ellipsoid(v_r, v_theta, v_phi, max_disc: int = 2, return_probability=False):
    """The halo's velocity ellipsoid from the best mixture (by BIC) and its
    anisotropy beta = 1 - (sigma_theta^2 + sigma_phi^2) / (2 sigma_r^2).

    beta = 0: isotropic; beta -> 1: purely radial orbits (a sausage);
    beta < 0: tangentially biased (a rotating disc)."""
    n = len(v_r)
    best, best_bic = None, np.inf
    for k in range(0, max_disc + 1):
        comps, ll, p = halo_mixture(v_r, v_theta, v_phi, k)
        bic = -2 * ll + p * np.log(n)
        if bic < best_bic:
            best, best_bic = (comps, k), bic
    comps, k = best
    h = comps[0]
    s_r, s_t, s_p = np.sqrt(np.diag(h["cov"]))
    out = dict(n_disc_components=int(k), weight=float(h["w"]), mean_v_phi=float(h["mean"][2]),
               sigma=[float(s_r), float(s_t), float(s_p)],
               beta=float(1 - (s_t**2 + s_p**2) / (2 * s_r**2)))
    if return_probability:
        X = np.column_stack([v_r, v_theta, v_phi]).astype(float)
        logp = np.array([np.log(c["w"]) + _log_gauss(X, c["mean"], c["cov"]) for c in comps])
        mx = logp.max(axis=0)
        out["_p_halo"] = np.exp(logp[0] - mx - np.log(np.sum(np.exp(logp - mx), axis=0)))
    return out


def moments_beta(v_r, v_theta, v_phi):
    """beta from plain velocity moments (no mixture): for a sample that has
    already been cut to the halo by velocity."""
    s_r, s_t, s_p = np.std(v_r), np.std(v_theta), np.std(v_phi)
    return dict(mean_v_phi=float(np.mean(v_phi)), sigma=[float(s_r), float(s_t), float(s_p)],
                beta=float(1 - (s_t**2 + s_p**2) / (2 * s_r**2)))


def fine_windows(mh, v_r, v_theta, v_phi, passing, centres=None, half=FINE_WINDOW_HALF):
    """beta, the halo's share, and how much of the halo a velocity cut keeps,
    in overlapping metallicity windows (0.3 dex wide, every 0.1 dex)."""
    centres = FINE_CENTRES if centres is None else centres
    rows = []
    for c in centres:
        w = (mh >= c - half) & (mh < c + half)
        row = {"c": float(c), "n": int(w.sum())}
        if w.sum() >= 60:
            e = velocity_ellipsoid(v_r[w], v_theta[w], v_phi[w], return_probability=True)
            p = e.pop("_p_halo")
            row.update(beta=round(e["beta"], 3), share=round(e["weight"], 3),
                       sig=[round(x) for x in e["sigma"]],
                       comp=round(float(np.sum(p[passing[w]]) / np.sum(p)), 3))
        rows.append(row)
    return rows


def beta_with_errors(v_r, v_theta, v_phi, n_boot=N_BOOT, seed=1, fit=velocity_ellipsoid):
    rng = np.random.default_rng(seed)
    main = fit(v_r, v_theta, v_phi)
    betas = []
    for _ in range(n_boot):
        i = rng.integers(0, len(v_r), len(v_r))
        betas.append(fit(v_r[i], v_theta[i], v_phi[i])["beta"])
    betas = np.array(betas)
    main["beta_err"] = float(0.5 * (np.percentile(betas, 84) - np.percentile(betas, 16)))
    return main


# --------------------------------------------------------------------------
# Selections in integrals of motion
# --------------------------------------------------------------------------


def in_box(x, rng):
    return (x > rng[0]) & (x < rng[1])


def substructure_masks(d):
    """Published selection boxes, each in the potential it was defined in."""
    h18 = REFERENCE["helmi2018"]
    k19b = REFERENCE["koppelman2019b"]
    k19 = REFERENCE["koppelman2019"]
    e5 = d["E_m17"] / 1e5
    return {
        "GSE (Helmi+18)": in_box(d["L_z"], h18["gse_lz"]) & (d["E_h17"] > h18["gse_energy_min"]),
        "GSE (Koppelman+19b)": in_box(e5, k19b["gse"]["energy"]) & in_box(d["circ_m17"], k19b["gse"]["circ"]),
        "Sequoia": in_box(e5, k19b["sequoia"]["energy"]) & in_box(d["circ_m17"], k19b["sequoia"]["circ"]),
        "Thamnos 1": in_box(e5, k19b["thamnos1"]["energy"]) & in_box(d["circ_m17"], k19b["thamnos1"]["circ"]),
        "Thamnos 2": in_box(e5, k19b["thamnos2"]["energy"]) & in_box(d["circ_m17"], k19b["thamnos2"]["circ"]),
        "Helmi streams": in_box(d["L_perp"], k19["hs_lperp"]) & in_box(d["L_z"], k19["hs_lz"]),
    }


# --------------------------------------------------------------------------
# The analysis
# --------------------------------------------------------------------------


def analyse() -> tuple[dict, dict]:
    fast, mp, pot_h17, pot_m17 = load_samples()
    res: dict = dict(n_fast_file=fast.pop("_n_file"), n_fast_within_5kpc=int(len(fast["ra"])))

    # ---- The kinematic route (Helmi et al. 2018) ---------------------------
    near = fast["good"] & (fast["distance"] < HALO_MAX_DISTANCE_KPC)
    halo = near & (fast["speed_lsr"] > HALO_SPEED_KMS)
    rnd = fast["good"] & fast["is_random"]
    rnd_near = rnd & (fast["distance"] < HALO_MAX_DISTANCE_KPC)
    res["n_halo_2p5kpc"] = int(halo.sum())
    res["halo_fraction_random"] = float(np.mean(fast["speed_lsr"][rnd_near] > HALO_SPEED_KMS))
    res["n_random_2p5kpc"] = int(rnd_near.sum())
    # The halo's velocity ellipsoid, with the velocity cut: the sausage.
    res["halo_kinematic_ellipsoid"] = beta_with_errors(
        fast["v_r"][halo], fast["v_theta"][halo], fast["v_phi"][halo], fit=moments_beta)
    res["halo_kinematic_moments"] = dict(
        mean_v_phi=float(np.mean(fast["v_phi"][halo])),
        sigma_v_R=float(np.std(fast["v_R"][halo])), sigma_v_phi=float(np.std(fast["v_phi"][halo])),
        sigma_v_z=float(np.std(fast["v_z"][halo])),
        frac_retrograde=float(np.mean(fast["v_phi"][halo] < 0)),
    )
    # Integrals of motion within 5 kpc (Helmi et al. 2018, Extended Data Fig. 1).
    halo5 = fast["good"] & (fast["distance"] < 5.0) & (fast["speed_lsr"] > HALO_SPEED_KMS) \
        & (fast["phot_bp_rp_excess_factor"] < 1.27)
    masks = substructure_masks(fast)
    res["n_halo_5kpc"] = int(halo5.sum())
    res["substructure_counts_5kpc"] = {k: int(np.sum(m & halo5)) for k, m in masks.items()}
    res["substructure_counts_2p5kpc"] = {k: int(np.sum(m & halo)) for k, m in masks.items()}
    gse = masks["GSE (Helmi+18)"] & halo
    res["gse_fraction_of_halo_2p5kpc"] = float(gse.sum() / halo.sum())

    # ---- The chemical route (no velocity cut) ------------------------------
    # Stars with trustworthy Gaia RVS metallicities, from both files (the
    # metal-poor file has every such star below [M/H] = -0.7; the random 2%
    # supplies the metal-rich disc for reference). No cut on velocity.
    chem_sample = trusted_giants
    chem_mp = chem_sample(mp)
    chem_rnd = rnd & chem_sample(fast)
    # The diagnosis that motivated the giant-only cut, kept as a result.
    loose = mp["good"] & mp["chem_ok"] & (mp["distance"] < CHEM_MAX_DISTANCE_KPC) & (mp["mh"] < -1.66)
    cold = (mp["v_phi"] > 180) & (np.abs(mp["v_R"]) < 60) & (np.abs(mp["v_z"]) < 50)
    res["gspspec_diagnosis"] = dict(
        n_metal_poor_medium=int(loose.sum()),
        cold_disc_fraction=float(np.mean(cold[loose])),
        cold_disc_median_teff=float(np.median(mp["teff_gspspec"][loose & cold])),
        others_median_teff=float(np.median(mp["teff_gspspec"][loose & ~cold])),
        giants_cold_disc_fraction=float(np.mean(cold[chem_mp & (mp["mh"] < -1.66)])),
    )
    rows = []
    edges = np.array(FEH_BINS)
    for lo, hi in zip(edges[:-1], edges[1:]):
        if hi <= -0.7:
            src, s = mp, chem_mp & (mp["mh"] >= lo) & (mp["mh"] < hi)
        else:
            src, s = fast, chem_rnd & (fast["mh"] >= lo) & (fast["mh"] < hi)
        n = int(s.sum())
        row = dict(feh_lo=float(lo), feh_hi=float(hi), n=n, source="metal-poor file"
                   if src is mp else "random 2%")
        if n >= 60:
            row.update(beta_with_errors(src["v_r"][s], src["v_theta"][s], src["v_phi"][s]))
            row["frac_would_pass_speed_cut"] = float(np.mean(src["speed_lsr"][s] > HALO_SPEED_KMS))
            # How well would the usual velocity cut have found this bin's
            # halo? The mixture gives each star a probability of belonging to
            # the halo: sum them over the stars the cut keeps (completeness)
            # and compare with how many stars it keeps (purity).
            p = velocity_ellipsoid(src["v_r"][s], src["v_theta"][s], src["v_phi"][s],
                                   return_probability=True)["_p_halo"]
            passing = src["speed_lsr"][s] > HALO_SPEED_KMS
            row["cut_completeness"] = float(np.sum(p[passing]) / np.sum(p))
            row["cut_purity"] = float(np.sum(p[passing]) / max(passing.sum(), 1))
        rows.append(row)
    res["beta_by_metallicity"] = rows
    # Coarser bins, for the sausage figure.
    disp = []
    for lo, hi in DISPLAY_BINS:
        s = chem_mp & (mp["mh"] >= lo) & (mp["mh"] < hi)
        disp.append(dict(feh_lo=lo, feh_hi=hi, n=int(s.sum()),
                         **beta_with_errors(mp["v_r"][s], mp["v_theta"][s], mp["v_phi"][s])))
    res["beta_display_bins"] = disp
    res["beta_fine_windows"] = fine_windows(
        mp["mh"][chem_mp], mp["v_r"][chem_mp], mp["v_theta"][chem_mp], mp["v_phi"][chem_mp],
        mp["speed_lsr"][chem_mp] > HALO_SPEED_KMS)

    # ---- Chemistry of the debris --------------------------------------------
    # Every star in the fast file has its GSP-Spec values (if any), whatever
    # its metallicity, so kinematically selected GSE stars get an unbiased
    # metallicity distribution. Same trusted giants as above. Two GSE
    # selections: Helmi et al.'s box in (L_z, E) among |v - v_LSR| > 210 km/s
    # stars, which reaches slightly prograde orbits where heated thick-disc
    # stars ("the Splash") also live; and Koppelman et al.'s (2019b)
    # near-radial orbits (|circularity| < 0.2), which avoid them.
    giants = chem_sample(fast) & np.isfinite(fast["alpha"])
    gse_sel = {
        "Koppelman+19b": masks["GSE (Koppelman+19b)"] & giants,
        "Helmi+18": masks["GSE (Helmi+18)"] & (fast["speed_lsr"] > HALO_SPEED_KMS) & giants,
    }
    # Chemistry of every published box (halo stars within 5 kpc): accreted
    # debris should be metal-poor; a box full of heated thick-disc stars will
    # sit near the thick disc's [M/H] ~ -0.5.
    res["box_chemistry"] = {}
    for k, m in masks.items():
        sel = m & halo5 & giants
        res["box_chemistry"][k] = dict(
            n_giants=int(sel.sum()),
            mh_median=float(np.median(fast["mh"][sel])) if sel.sum() >= 10 else None)
    res["gse_chemistry"] = {
        k: dict(n=int(m.sum()), mh_median=float(np.median(fast["mh"][m])),
                mh_quartiles=np.percentile(fast["mh"][m], [25, 75]).tolist())
        for k, m in gse_sel.items()}
    # Alpha sequences at fixed metallicity: GSE versus thick-disc giants
    # (prograde, 80-180 km/s from the LSR: the thick disc's velocities).
    thick = giants & (fast["v_phi"] > 120) & (fast["speed_lsr"] > 80) & (fast["speed_lsr"] < 180)
    rng = np.random.default_rng(4)

    def med_err(x):
        b = [np.median(rng.choice(x, len(x))) for _ in range(300)]
        return float(np.median(x)), float(np.std(b))

    rows = []
    for lo, hi in ALPHA_BINS:
        a = gse_sel["Koppelman+19b"] & in_box(fast["mh"], (lo, hi))
        t = thick & in_box(fast["mh"], (lo, hi))
        if a.sum() >= 20 and t.sum() >= 20:
            ma, ea = med_err(fast["alpha"][a])
            mt, et = med_err(fast["alpha"][t])
            rows.append(dict(mh_lo=lo, mh_hi=hi, n_gse=int(a.sum()), n_thick=int(t.sum()),
                             alpha_gse=ma, alpha_gse_err=ea, alpha_thick=mt, alpha_thick_err=et,
                             difference=ma - mt, difference_err=float(np.hypot(ea, et))))
    res["alpha_by_metallicity"] = rows
    d = np.array([r["difference"] for r in rows])
    e = np.array([r["difference_err"] for r in rows])
    res["alpha_difference_combined"] = dict(value=float(np.sum(d / e**2) / np.sum(1 / e**2)),
                                            err=float(1 / np.sqrt(np.sum(1 / e**2))))
    gse = masks["GSE (Helmi+18)"] & halo

    # ---- The halo's colour-magnitude diagram (Helmi et al. 2018, Fig. 2c) ---
    # Halo stars with precise distances and little dust: Gaia DR2 found they
    # form two parallel sequences, the bluer made of accreted (GSE) stars.
    ebv = ph.ebv_sfd(fast["ra"], fast["dec"])
    hrd = (halo & (fast["parallax"] / fast["parallax_error"] > 10) & (ebv < 0.015)
           & ph.clean_colour(fast["phot_bp_rp_excess_factor"], fast["bp_rp"],
                             fast["phot_g_mean_mag"]))
    res["n_hrd"] = int(hrd.sum())
    fast["abs_g"] = fast["phot_g_mean_mag"] + 5 * np.log10(fast["parallax"] / 100.0)

    artefacts = dict(fast=fast, mp=mp, halo=halo, halo5=halo5, rnd=rnd, rnd_near=rnd_near,
                     masks=masks, chem_mp=chem_mp, chem_rnd=chem_rnd, gse=gse, gse_sel=gse_sel,
                     thick=thick, giants=giants, hrd=hrd, pot_h17=pot_h17, pot_m17=pot_m17)
    return res, artefacts
