"""The wide-binary samples of the published tests, rebuilt from public data.

* ``chae2023``: Chae (2023)'s main sample from El-Badry et al. (2021):
  MS+MS, R < 0.01, distances agreeing within 3 sigma, every proper-motion
  component measured to better than 1% (or 0.5%), 0.2 < s < 30 kAU,
  4 < M_G < 14 for both stars, within 200 pc (26,615 pairs in his paper).
* ``banik2024``: Banik et al. (2024)'s 8,611 pairs, from our rebuild of
  Pittordis & Sutherland's (2023) search (pairs.py) plus their extra cuts.
* ``cookson2026``: Cookson et al. (2026)'s high-purity sample within 130 pc
  (1,421 pairs).

Each builder returns a dict of per-pair columns (suffix 1/2 for the two
stars) and a cut-flow list of (step, number left), so every published count
can be compared with ours.
"""

from __future__ import annotations

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from config import RAW_DIR
import kinematics as KIN
import masses as ML

# --------------------------------------------------------------------------
# Dust: Schlegel, Finkbeiner & Davis (1998)
# --------------------------------------------------------------------------
_SFD = {}


def sfd_ebv(l_deg, b_deg):
    """SFD E(B-V) (mag) at Galactic (l, b), nearest pixel of the 4096^2 maps."""
    l_deg, b_deg = np.atleast_1d(l_deg).astype(float), np.atleast_1d(b_deg).astype(float)
    out = np.empty(len(l_deg))
    for hemi, sel in (("ngp", b_deg >= 0), ("sgp", b_deg < 0)):
        if not np.any(sel):
            continue
        if hemi not in _SFD:
            with fits.open(RAW_DIR / "sfd" / f"SFD_dust_4096_{hemi}.fits") as h:
                _SFD[hemi] = (h[0].data.astype(np.float32), WCS(h[0].header))
        data, w = _SFD[hemi]
        x, y = w.wcs_world2pix(l_deg[sel], b_deg[sel], 0)
        xi = np.clip(np.round(x).astype(int), 0, data.shape[1] - 1)
        yi = np.clip(np.round(y).astype(int), 0, data.shape[0] - 1)
        out[sel] = data[yi, xi]
    return out


def dust_layer_AG(l_deg, b_deg, d_pc, h_pc: float = 125.0):
    """A_G to distance d: the SFD column scaled by the fraction of an
    exponential dust layer (scale height h) in front of the star. An
    approximation, for the |b| > 15 deg, d < 250 pc stars where it matters at
    the 0.01-0.1 mag level; A_G = 0.789 A_V = 0.789 * 3.1 E(B-V)."""
    ebv = sfd_ebv(l_deg, b_deg)
    s = np.abs(np.sin(np.radians(b_deg)))
    frac = (1 - np.exp(-d_pc * np.maximum(s, 0.05) / h_pc))
    return 0.789 * 3.1 * ebv * frac


def mg(G, dist_pc):
    return G - 5.0 * np.log10(dist_pc / 10.0)


# --------------------------------------------------------------------------
# Chae (2023) from El-Badry et al. (2021)
# --------------------------------------------------------------------------

def chae2023(eb, *, dmax_pc: float = 200.0, pm_rel: float = 0.01, mg_range=(4.0, 14.0),
             s_range_kau=(0.2, 30.0), dust: bool = True, r_max: float = 0.01,
             bayestar_footprint: bool = True) -> tuple[dict, list]:
    """``eb``: the El-Badry 'near' table (elbadry.load('near')).

    Chae corrected magnitudes beyond 80 pc with the Bayestar 3D dust map,
    which covers only dec > -28 deg (his fig. 3), so his 200 pc sample has no
    pairs farther than 80 pc south of that ("because of ... the availability
    of dust extinction correction", Chae 2024a). ``bayestar_footprint``
    reproduces that; our own correction is the SFD-scaled dust layer."""
    g = lambda c: np.asarray(eb[c])
    flow = [("El-Badry et al. pairs with parallax1 > 3 mas", len(eb))]
    bt = np.char.strip(g("binary_type").astype(str))
    m = bt == "MSMS"
    flow.append(("both main sequence (MSMS)", int(m.sum())))
    m &= g("R_chance_align") < r_max
    flow.append((f"R < {r_max}", int(m.sum())))
    d1, d2 = 1000.0 / g("parallax1"), 1000.0 / g("parallax2")
    s1, s2 = d1 * g("parallax_error1") / g("parallax1"), d2 * g("parallax_error2") / g("parallax2")
    m &= np.abs(d1 - d2) < 3.0 * np.hypot(s1, s2)
    flow.append(("distances agree within 3 sigma", int(m.sum())))
    rel = np.maximum.reduce([np.abs(g(f"{c}_error{w}") / g(f"{c}{w}")) for c in ("pmra", "pmdec")
                             for w in (1, 2)])
    m &= rel < pm_rel
    flow.append((f"all PM components to < {pm_rel:g} relative error", int(m.sum())))
    s = g("sep_AU") / 1000.0
    m &= (s > s_range_kau[0]) & (s < s_range_kau[1])
    flow.append((f"{s_range_kau[0]} < s < {s_range_kau[1]} kAU", int(m.sum())))
    dM = KIN.weighted_distance(g("parallax1"), g("parallax_error1"), g("parallax2"), g("parallax_error2"))
    m &= dM < dmax_pc
    flow.append((f"d_M < {dmax_pc:g} pc", int(m.sum())))
    if bayestar_footprint:
        m &= (dM < 80.0) | (g("dec1") > -28.0)
        flow.append(("beyond 80 pc only dec > -28 deg (Bayestar footprint)", int(m.sum())))
    AG = np.zeros(len(eb))
    if dust:
        AG[m] = dust_layer_AG(g("l1")[m], g("b1")[m], dM[m])
        AG = np.where(dM > 80.0, AG, 0.0)       # Chae: dust correction beyond 80 pc
    MG1 = mg(g("phot_g_mean_mag1"), dM) - AG
    MG2 = mg(g("phot_g_mean_mag2"), dM) - AG
    m &= (MG1 > mg_range[0]) & (MG1 < mg_range[1]) & (MG2 > mg_range[0]) & (MG2 < mg_range[1])
    flow.append((f"{mg_range[0]} < M_G < {mg_range[1]} for both", int(m.sum())))
    out = {c: g(c)[m] for c in eb.colnames}
    out.update(dist=dM[m], MG1=MG1[m], MG2=MG2[m], AG=AG[m])
    return out, flow


if __name__ == "__main__":
    import elbadry
    eb = elbadry.load("near")
    for kw in (dict(), dict(dust=False), dict(pm_rel=0.005)):
        out, flow = chae2023(eb, **kw)
        print(kw)
        for step, nleft in flow:
            print(f"   {nleft:>9,}  {step}")


# --------------------------------------------------------------------------
# PS23-style pairs joined to their full Gaia DR3 rows
# --------------------------------------------------------------------------
STAR_FIELDS = ("source_id", "ra", "dec", "l", "b", "parallax", "parallax_error", "pmra", "pmra_error",
               "pmdec", "pmdec_error", "parallax_pmra_corr", "parallax_pmdec_corr", "pmra_pmdec_corr",
               "phot_g_mean_mag", "phot_g_mean_flux_over_error", "bp_rp", "ruwe", "ipd_frac_multi_peak",
               "astrometric_chi2_al", "astrometric_n_good_obs_al", "radial_velocity",
               "radial_velocity_error", "non_single_star", "ipd_gof_harmonic_amplitude", "duplicated_source")


def ps23_table(shifted: bool = False):
    """PS23 candidate pairs (after steps 1-5) with full Gaia columns per star;
    star 1 is the brighter in G. For the shifted catalogue, star 2's sky
    position is its shifted one (so r_sky and the perspective geometry are
    those of the chance alignment) while its other columns are its own."""
    from astropy.table import Table
    from config import GAIA_DIR
    import pairs as PR
    pr, flow = PR.cached_pairs(shifted)
    keep = pr["keep"]
    full = Table.read(GAIA_DIR / "ps23_stars_full.fits")
    sid = np.asarray(full["source_id"])
    order = np.argsort(sid)

    def rows(ids):
        j = order[np.searchsorted(sid[order], ids)]
        assert np.array_equal(sid[j], ids)
        return j

    ja, jb = rows(pr["source_id_a"][keep]), rows(pr["source_id_b"][keep])
    cols = {}
    for f in STAR_FIELDS:
        va, vb = np.asarray(full[f])[ja], np.asarray(full[f])[jb]
        cols[f + "_a"], cols[f + "_b"] = va, vb
    if shifted:
        rb, db = PR.shifted_positions(cols["ra_b"], cols["dec_b"])
        cols["ra_b"], cols["dec_b"] = rb, db
    # star 1 = brighter
    swap = cols["phot_g_mean_mag_a"] > cols["phot_g_mean_mag_b"]
    out = {}
    for f in STAR_FIELDS:
        a, b = cols[f + "_a"], cols[f + "_b"]
        out[f + "1"] = np.where(swap, b, a)
        out[f + "2"] = np.where(swap, a, b)
    out["rsky_ps23_au"] = pr["rsky_au"][keep]
    out["dv_ps23_kms"] = pr["dv_kms"][keep]
    return out, flow


def arenou_ok(chi2, ngood, G):
    """Arenou et al. (2018) eq. 1, as used by PS23 (their eq. 1)."""
    u = np.sqrt(chi2 / np.maximum(ngood - 5, 1))
    return u <= 1.2 * np.maximum(1.0, np.exp(-0.2 * (G - 19.5)))


def wd_below(MG, bprp):
    """El-Badry et al. (2021) white-dwarf side of the CMD."""
    return np.isfinite(bprp) & (MG > 3.25 * bprp + 9.625)


def pair_kinematics(s: dict, mass_relation: str = "banik", perspective: bool = True,
                    mass_override=None) -> dict:
    """Nominal r_sky (kAU), v~, masses for each pair in s (suffix 1/2)."""
    dist = KIN.weighted_distance(s["parallax1"], s["parallax_error1"], s["parallax2"], s["parallax_error2"])
    rv, erv = KIN.systemic_rv(s["radial_velocity1"], s["radial_velocity_error1"], s["radial_velocity2"],
                              s["radial_velocity_error2"])
    vsky, rsky = KIN.sky_velocity(s["ra1"], s["dec1"], s["pmra1"], s["pmdec1"], s["ra2"], s["dec2"],
                                  s["pmra2"], s["pmdec2"], dist, dist, np.nan_to_num(rv),
                                  perspective=perspective)
    MG1, MG2 = mg(s["phot_g_mean_mag1"], dist), mg(s["phot_g_mean_mag2"], dist)
    if mass_override is not None:
        m1, m2 = mass_override
    else:
        m1, m2 = ML.mass_from_MG(MG1, mass_relation), ML.mass_from_MG(MG2, mass_relation)
    vs = np.linalg.norm(vsky, axis=-1)
    return dict(dist=dist, rv_sys=rv, rv_sys_err=erv, rsky_kau=rsky / 1000.0, vsky=vs, MG1=MG1, MG2=MG2,
                m1=m1, m2=m2, mtot=m1 + m2, vt=KIN.vtilde(vs, rsky, m1 + m2),
                vc=np.sqrt(KIN.GM_KMS2_AU * (m1 + m2) / rsky))


def _take(d: dict, m):
    return {k: (v[m] if isinstance(v, np.ndarray) and v.shape[:1] == m.shape else v) for k, v in d.items()}


def banik2024(s: dict, flow0, *, plx_min: float = 4.0, rsky_range=(2.0, 30.0), vt_max: float = 5.0,
              mass_range=(0.464, 4.31), av_max: float = 0.5, ipd_max: int = 2, need_rv: bool = True,
              rv_match: bool = True, err_cut: bool = True, n_trials: int = 512, perspective: bool = True,
              plx_err_scale: float = 1.0, seed: int = 1):
    """Banik et al. (2024) Sections 2.1-2.4 on our PS23 rebuild."""
    flow = [(k, v) for k, v in flow0.items()]
    m = np.ones(len(s["ra1"]), bool)
    m &= arenou_ok(s["astrometric_chi2_al1"], s["astrometric_n_good_obs_al1"], s["phot_g_mean_mag1"])
    m &= arenou_ok(s["astrometric_chi2_al2"], s["astrometric_n_good_obs_al2"], s["phot_g_mean_mag2"])
    flow.append(("Arenou et al. (2018) astrometric cut on both stars (PS23 CWB-EDR3: 73,159)", int(m.sum())))
    m &= (s["parallax1"] > plx_min) & (s["parallax2"] > plx_min)
    flow.append((f"parallax > {plx_min:g} mas for both", int(m.sum())))
    k = pair_kinematics(s, "banik", perspective=perspective)
    m &= (k["rsky_kau"] > rsky_range[0]) & (k["rsky_kau"] < rsky_range[1])
    flow.append((f"{rsky_range[0]:g} < r_sky < {rsky_range[1]:g} kAU", int(m.sum())))
    m &= k["vt"] < vt_max
    flow.append((f"v~ < {vt_max:g}", int(m.sum())))
    m &= (k["mtot"] > mass_range[0]) & (k["mtot"] < mass_range[1])
    flow.append((f"{mass_range[0]} < M_tot < {mass_range[1]} Msun", int(m.sum())))
    av1 = 3.1 * sfd_ebv(s["l1"], s["b1"])
    av2 = 3.1 * sfd_ebv(s["l2"], s["b2"])
    m &= (av1 < av_max) & (av2 < av_max)
    flow.append((f"SFD A_V < {av_max:g} for both (Banik: -8.5%)", int(m.sum())))
    m &= ~wd_below(k["MG1"], s["bp_rp1"]) & ~wd_below(k["MG2"], s["bp_rp2"])
    flow.append(("no white dwarf (Banik: -66)", int(m.sum())))
    m &= (np.nan_to_num(s["ipd_frac_multi_peak1"], nan=99) <= ipd_max) & \
         (np.nan_to_num(s["ipd_frac_multi_peak2"], nan=99) <= ipd_max)
    flow.append((f"ipd_frac_multi_peak <= {ipd_max} for both", int(m.sum())))
    has1 = np.isfinite(s["radial_velocity1"])
    has2 = np.isfinite(s["radial_velocity2"])
    if need_rv:
        m &= has1 | has2
        flow.append(("radial velocity for at least one star", int(m.sum())))
    if rv_match:
        both = has1 & has2
        drv = np.abs(np.nan_to_num(s["radial_velocity1"]) - np.nan_to_num(s["radial_velocity2"]))
        lim = np.maximum(3.0 * np.hypot(np.nan_to_num(s["radial_velocity_error1"]),
                                        np.nan_to_num(s["radial_velocity_error2"])), 3.0 * k["vc"])
        m &= ~both | (drv <= lim)
        flow.append(("RVs agree when both known (Banik: -308)", int(m.sum())))
    sub = _take({**s, **k}, m)
    if err_cut:
        med, std = KIN.mc_vtilde(sub, lambda x: ML.mass_from_MG(x, "banik"), n_trials=n_trials, seed=seed,
                                 plx_err_scale=plx_err_scale)
        ok = std <= 0.1 * np.maximum(1.0, sub["vt"] / 2.0)
        sub = _take(sub, ok)
        sub["vt_err"] = std[ok]
        flow.append(("sigma(v~) <= 0.1 max(1, v~/2) (Banik: final 8,611)", int(ok.sum())))
    return sub, flow


# --------------------------------------------------------------------------
# Cookson et al. (2026): the "quality framework" sample
# --------------------------------------------------------------------------
MG_SUN_G = 4.67          # Willmer (2018), the Sun's absolute G magnitude


def cookson_mass(MG):
    """Cookson et al. (2026) eq. 12 (+0.05 Msun below 0.7 Msun); used for
    every star (they use FLAME masses where Gaia has them, ~19%)."""
    m = 10.0 ** (0.0725 * (MG_SUN_G - MG))
    return np.where(m < 0.7, m + 0.05, m)


def hr_parallelogram(MG, bprp, ridge_c, half_height: float = 0.4):
    """0.5 < BP-RP < 3.5, 4 < M_G < 14, within +-0.4 mag of a line of slope
    3.3 (their 0.8-mag-high parallelogram). The intercept is not published;
    ``ridge_c`` is either one number (a straight band through the sample's
    main-sequence ridge) or a callable c(BP-RP) for a band that follows the
    ridge's curvature."""
    c = ridge_c(bprp) if callable(ridge_c) else ridge_c
    return ((bprp > 0.5) & (bprp < 3.5) & (MG > 4) & (MG < 14)
            & (np.abs(MG - (3.3 * bprp + c)) < half_height))


def ms_ridge_offset(MG, bprp, m):
    """c(BP-RP): running median of M_G - 3.3 (BP-RP) along the main sequence."""
    x = bprp[m]
    y = (MG - 3.3 * bprp)[m]
    edges = np.linspace(0.5, 3.5, 13)
    mid = 0.5 * (edges[1:] + edges[:-1])
    med = np.array([np.median(y[(x >= a) & (x < b)]) if np.sum((x >= a) & (x < b)) > 20 else np.nan
                    for a, b in zip(edges[:-1], edges[1:])])
    ok = np.isfinite(med)
    return lambda c: np.interp(c, mid[ok], med[ok])


def cookson2026(s: dict, flow0, *, dmax_pc: float = 130.0, rsky_range=(1.0, 30.0), vt_max: float = 2.5,
                ruwe_max: float = 1.25, drv_max: float = 10.0, plx_snr: float = 40.0, ipd_zero: bool = False,
                vt_err_max: float = 0.1, n_trials: int = 512, seed: int = 2, hr_band: str = "ridge"):
    flow = [(k, v) for k, v in flow0.items()]
    m = (s["parallax1"] / s["parallax_error1"] > plx_snr) & (s["parallax2"] / s["parallax_error2"] > plx_snr)
    flow.append((f"parallax S/N > {plx_snr:g} for both (their catalogue C25)", int(m.sum())))
    k = pair_kinematics(s, mass_override=None, mass_relation="banik")
    MG1, MG2 = k["MG1"], k["MG2"]
    m1, m2 = cookson_mass(MG1), cookson_mass(MG2)
    k.update(m1=m1, m2=m2, mtot=m1 + m2)
    k["vt"] = KIN.vtilde(k["vsky"], k["rsky_kau"] * 1000.0, k["mtot"])
    k["vc"] = np.sqrt(KIN.GM_KMS2_AU * k["mtot"] / (k["rsky_kau"] * 1000.0))
    m &= k["dist"] < dmax_pc
    flow.append((f"distance < {dmax_pc:g} pc", int(m.sum())))
    m &= (k["rsky_kau"] > rsky_range[0]) & (k["rsky_kau"] < rsky_range[1])
    flow.append((f"{rsky_range[0]:g} < r_sky < {rsky_range[1]:g} kAU", int(m.sum())))
    m &= (s["ruwe1"] < ruwe_max) & (s["ruwe2"] < ruwe_max)
    flow.append((f"RUWE < {ruwe_max:g} for both", int(m.sum())))
    both = np.isfinite(s["radial_velocity1"]) & np.isfinite(s["radial_velocity2"])
    m &= both & (np.abs(s["radial_velocity1"] - s["radial_velocity2"]) < drv_max)
    flow.append((f"both RVs, |dRV| < {drv_max:g} km/s", int(m.sum())))
    if hr_band == "ridge":
        allMG = np.concatenate([MG1[m], MG2[m]])
        allc = np.concatenate([s["bp_rp1"][m], s["bp_rp2"][m]])
        c = ms_ridge_offset(allMG, allc, np.isfinite(allc))
        label = "following the main-sequence ridge"
    else:
        c = float(np.median((MG1 - 3.3 * s["bp_rp1"])[m & (s["bp_rp1"] > 0.5) & (s["bp_rp1"] < 3.5)]))
        label = f"straight, intercept {c:.2f}"
    m &= hr_parallelogram(MG1, s["bp_rp1"], c) & hr_parallelogram(MG2, s["bp_rp2"], c)
    flow.append((f"HR band 0.8 mag high, slope 3.3, {label}, for both", int(m.sum())))
    if ipd_zero:
        m &= (s["ipd_frac_multi_peak1"] == 0) & (s["ipd_frac_multi_peak2"] == 0)
        flow.append(("ipd_frac_multi_peak = 0 for both", int(m.sum())))
    m &= k["vt"] < vt_max
    flow.append((f"v~ < {vt_max:g}", int(m.sum())))
    sub = _take({**s, **k}, m)
    med, std = KIN.mc_vtilde(sub, cookson_mass, n_trials=n_trials, seed=seed)
    ok = std < vt_err_max
    sub = _take(sub, ok)
    sub["vt_err"] = std[ok]
    flow.append((f"sigma(v~) < {vt_err_max:g} (Cookson et al.: 1,421)", int(ok.sum())))
    return sub, flow
