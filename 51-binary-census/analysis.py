"""The binary-star census: every measurement, from cached data.

``run()`` returns one results dict that ``figures.py`` draws from and
``summarize()`` turns into data/results.json. Expensive steps (the gaiamock
simulations, the twin MCMC) are cached under data/cache/ and reused when their
inputs have not changed.
"""

from __future__ import annotations

import hashlib
import json
import pickle
import zlib

import numpy as np

import config as C
import literature as L
import pe
import raghavan
import samples
import sb9
import wide

RNG = np.random.default_rng(C.SEED)


def _cached(name: str, key: dict, fn):
    """Run ``fn()`` unless a cache with the same key exists."""
    h = hashlib.sha1(json.dumps(key, sort_keys=True, default=str).encode()).hexdigest()[:12]
    path = C.CACHE_DIR / f"{name}_{h}.pkl"
    if path.exists():
        return pickle.loads(path.read_bytes())
    out = fn()
    path.write_bytes(pickle.dumps(out))
    return out


# ---------------------------------------------------------------------------
# 1. Period-eccentricity diagrams and tidal circularization
# ---------------------------------------------------------------------------

def sb9_circularization() -> dict:
    t = sb9.load_sb9()
    ms = sb9.main_sequence_mask(t)
    fgk = ms & np.isin(t["sp_class"], ["F", "G", "K"])
    p, e = np.asarray(t["period"]), np.asarray(t["ecc"])
    ee = np.asarray(t["ecc_err"])
    fixed = np.asarray(t["ecc_fixed"])
    sel = fgk & np.isfinite(p) & np.isfinite(e)
    mm = pe.fit_meibom_mathieu(p[sel], e[sel])
    boot = pe.bootstrap(lambda pp, ee_: pe.fit_meibom_mathieu(pp, ee_, p_grid=np.arange(0.5, 40, 0.05))["p_circ"],
                        [p[sel], e[sel]], n_boot=200)
    circ = fixed | pe.consistent_with_circular(e, np.where(np.isfinite(ee), ee, 0.0))
    cf = pe.circular_fraction_curve(p[sel], circ[sel])
    giants = (t["lum_group"] == "giant") & np.isfinite(p) & np.isfinite(e)
    return {"table": t, "fgk_dwarf": sel, "giant": np.asarray(giants), "circular": circ,
            "mm": {k: v for k, v in mm.items() if k not in ("cost", "p_grid")},
            "mm_boot": boot, "p_circ": mm["p_circ"], "p_circ_err": float(np.std(boot)),
            "circ_frac": cf, "n_fgk": int(sel.sum())}


#: Bashi, Mazeh & Faigler (2023) Teff bins (K), for a like-for-like comparison.
BASHI_TEFF_BINS = [(5300, 5500), (5500, 5600), (5600, 5700), (5700, 5800), (5800, 5900), (5900, 6000),
                   (6000, 6100), (6100, 6200), (6200, 6300), (6300, 6400), (6400, 6500), (6500, 6600),
                   (6600, 6700), (6700, 6800), (6800, 6900), (6900, 7100), (7100, 7300), (7300, 7700)]
BASHI_PCUT = [7.00, 7.19, 6.57, 6.16, 6.52, 5.90, 5.44, 4.58, 4.52, 3.95, 4.22, 3.39, 3.19, 2.55, 2.90, 2.94, 3.46,
              3.28]
BASHI_PCUT_ERR = [0.77, 0.86, 0.68, 0.49, 0.63, 0.47, 0.41, 0.32, 0.34, 0.40, 0.43, 0.44, 0.46, 0.41, 0.63, 0.38,
                  0.43, 0.38]
#: Meibom & Mathieu (2005) Table 3: (population, age in Gyr, circularization period, +err, -err).
MM05_TABLE3 = [("PMS", 10 ** -2.5, 7.1, 1.2, 1.2), ("Pleiades", 0.1, 7.2, 1.8, 1.9), ("M35", 0.15, 10.2, 1.0, 1.5),
               ("Hyades/Praesepe", 0.63, 3.2, 1.2, 1.2), ("M67", 10 ** 0.6, 12.1, 1.0, 1.5),
               ("NGC 188", 10 ** 0.8, 14.5, 1.4, 2.2), ("Field", 9.0, 10.3, 1.5, 3.1), ("Halo", 10.0, 15.6, 2.3, 3.2)]


SIG_MIN_MM = 20.0
SIG_MIN_ENVELOPE = 40.0


def _pcirc_with_errors(period, ecc, n_boot=150, seed=1):
    r = pe.fit_meibom_mathieu(period, ecc, p_grid=np.arange(0.5, 40, 0.02))
    if r["n"] < 10:
        return {"p_circ": np.nan, "lo": np.nan, "hi": np.nan, "n": r["n"]}
    boot = pe.bootstrap(lambda a, b: pe.fit_meibom_mathieu(a, b, p_grid=np.arange(0.5, 40, 0.05))["p_circ"],
                        [np.asarray(period), np.asarray(ecc)], n_boot=n_boot, seed=seed)
    return {"p_circ": r["p_circ"], "lo": float(np.percentile(boot, 16)), "hi": float(np.percentile(boot, 84)),
            "n": r["n"]}


def gaia_circularization() -> dict:
    """Tidal circularization in DR3's spectroscopic orbits of dwarfs within
    500 pc: by temperature (Bashi et al.'s test), by kinematic population and by
    open-cluster age. Eccentricities consistent with zero (Lucy & Sweeney's
    e < 2.45 sigma_e) are set to zero before fitting, as Meibom & Mathieu's
    estimator assumes measured circular orbits are circular."""
    def compute():
        g = pe.gaia_spectroscopic_orbits()
        poe = np.asarray(g["plx_best"]) / np.ma.filled(np.ma.asarray(g["parallax_error"], dtype=float), np.nan)
        d = np.asarray(g["dist_pc"])
        sig = np.ma.filled(np.ma.asarray(g["significance"], dtype=float), np.nan)
        # Gosset et al. (2025) advise significance K/sigma_K of 20-50 for science;
        # below ~40, DR3's SB1s include spurious short-period, highly eccentric
        # orbits (Bashi et al. 2022). The Meibom-Mathieu estimate barely moves
        # between cuts of 0, 20 and 40 (9.5, 9.1, 9.1 d); the envelope estimator
        # needs 40.
        # Double-lined (SB2) solutions are left out: 63% of DR3's SB2 orbits
        # shorter than 3 days have "significant" eccentricities beyond the
        # envelope where the stars would nearly touch at periastron (13% of
        # SB1s), so their eccentricities are not usable here (Bashi et al.
        # 2023 also used SB1 orbits only).
        single_lined = np.isin(np.asarray(g["nss_solution_type"]), ("SB1", "SB1C", "AstroSpectroSB1"))
        ms = ((np.asarray(g["dg_ms"]) > -1.25) & (np.asarray(g["dg_ms"]) < 0.75) & (poe > 10) & (d < 500)
              & np.isfinite(np.asarray(g["period"])) & (sig > SIG_MIN_MM) & single_lined)
        P = np.asarray(g["period"], dtype=float)
        e = np.asarray(g["ecc"], dtype=float)
        ee = np.asarray(g["ecc_err"], dtype=float)
        circ = np.asarray(g["circular_model"]) | pe.consistent_with_circular(e, ee)
        e_ls = np.where(circ, 0.0, e)
        teff = np.asarray(g["teff"], dtype=float)
        pop = pe.population(g)
        short = P < 100
        out = {"P": P[ms], "e": e[ms], "e_err": ee[ms], "circ": circ[ms], "teff": teff[ms], "pop": pop[ms],
               "type": np.asarray(g["nss_solution_type"])[ms], "n": int(ms.sum())}
        import binaries as B
        base = ((np.asarray(g["dg_ms"]) > -1.25) & (np.asarray(g["dg_ms"]) < 0.75) & (poe > 10) & (d < 500)
                & (sig > SIG_MIN_MM) & (P < 3))
        above = ~circ & (e > B.e_max_envelope(P) + 0.05)
        out["short_period_eccentricity_check"] = {
            t: {"n": int((base & (np.asarray(g["nss_solution_type"]) == t)).sum()),
                "frac_above_envelope": float(np.mean(above[base & (np.asarray(g["nss_solution_type"]) == t)]))}
            for t in ("SB1", "SB2", "AstroSpectroSB1")}
        out["all"] = _pcirc_with_errors(P[ms & short], e_ls[ms & short])
        out["circ_frac"] = pe.circular_fraction_curve(P[ms], circ[ms])
        tb = []
        for (lo, hi) in BASHI_TEFF_BINS:
            m = ms & short & (teff >= lo) & (teff < hi)
            r = _pcirc_with_errors(P[m], e_ls[m])
            r["teff_lo"], r["teff_hi"] = lo, hi
            tb.append(r)
        out["teff_bins"] = tb
        # Bashi et al.'s own estimator on the same bins: SB1 only, significantly
        # eccentric orbits (e > 2.45 sigma_e), P < 30 d
        sb1 = np.asarray(g["nss_solution_type"]) == "SB1"
        ecc_sig = ~circ & np.isfinite(e)
        tb_env = []
        for (lo, hi) in BASHI_TEFF_BINS:
            m = ms & sb1 & ecc_sig & (P < 30) & (teff >= lo) & (teff < hi) & (sig > SIG_MIN_ENVELOPE)
            r = pe.fit_envelope(P[m], e[m], tau_grid=np.arange(1.0, 6.01, 0.25))
            coarse = {"p_grid": np.arange(1.0, 12.0, 0.1), "tau_grid": np.arange(1.0, 6.01, 0.5)}
            boot = pe.bootstrap(lambda a, b: pe.fit_envelope(a, b, **coarse)["p_cut"], [P[m], e[m]], n_boot=40,
                                seed=3)
            r.update({"teff_lo": lo, "teff_hi": hi, "lo": float(np.nanpercentile(boot, 16)),
                      "hi": float(np.nanpercentile(boot, 84))})
            tb_env.append(r)
        out["teff_bins_envelope"] = tb_env
        pops = {}
        solar_like = (teff >= 5300) & (teff < 6300)
        for name in pe.POPULATIONS:
            m = ms & short & solar_like & (pop == name)
            pops[name] = _pcirc_with_errors(P[m], e_ls[m])
            pops[name]["circ_frac"] = pe.circular_fraction_curve(P[ms & solar_like & (pop == name)],
                                                                 circ[ms & solar_like & (pop == name)])
        out["populations"] = pops
        # temperature-matched comparison: thin-disc binaries reweighted to the
        # thick disc's Teff distribution (P_circ depends on Teff, and old stars
        # of a given colour differ), by resampling within 100 K bins
        rng = np.random.default_rng(5)
        thick = ms & short & solar_like & (pop == "thick disc")
        thin = ms & short & solar_like & (pop == "thin disc")
        edges_t = np.arange(5300, 6301, 100)
        it, ik = np.digitize(teff[thin], edges_t), np.digitize(teff[thick], edges_t)
        idx_thin = np.where(thin)[0]
        draws = []
        for b in np.unique(ik):
            cand = idx_thin[it == b]
            k = int(np.sum(ik == b))
            if len(cand) and k:
                draws.append(rng.choice(cand, size=4 * k, replace=True))
        matched = np.concatenate(draws)
        pops["thin disc (Teff-matched)"] = _pcirc_with_errors(P[matched], e_ls[matched])
        pops["thin disc (Teff-matched)"]["circ_frac"] = pe.circular_fraction_curve(P[matched], circ[matched])
        for name in list(pops):
            if name == "thin disc (Teff-matched)":
                m_all = matched
            else:
                m_all = np.where(ms & solar_like & (pop == name))[0]
            bs = []
            for _ in range(60):
                i = rng.choice(m_all, size=len(m_all), replace=True)
                bs.append(pe.circular_fraction_curve(P[i], circ[i])["p50"])
            pops[name]["p50"] = pe.circular_fraction_curve(P[m_all], circ[m_all])["p50"]
            pops[name]["p50_lo"], pops[name]["p50_hi"] = (float(np.percentile(bs, 16)), float(np.percentile(bs, 84)))
            pops[name]["n_all"] = int(len(m_all))
            pops[name]["median_teff"] = float(np.median(teff[m_all]))
        # open clusters (Hunt & Reffert 2023 members, matched by source_id)
        from astropy.table import Table
        mem = Table.read(C.VIZIER_DIR / "hunt23_members_nss.fits")
        same = np.asarray(mem["source_id"]).astype(np.int64) == np.asarray(mem["GaiaDR3"]).astype(np.int64)
        mem = mem[same & (np.ma.filled(np.ma.asarray(mem["Prob"], dtype=float), 0) > 0.5)]
        cl = Table.read(C.VIZIER_DIR / "hunt23_clusters.fits")
        age = dict(zip(np.char.strip(np.asarray(cl["Name"]).astype(str)),
                       np.ma.filled(np.ma.asarray(cl["logAge50"], dtype=float), np.nan)))
        mem_names = np.char.strip(np.asarray(mem["Name"]).astype(str))
        mage = {int(s): age.get(n, np.nan) for s, n in zip(np.asarray(mem["source_id"]), mem_names)}
        gid = np.asarray(g["source_id"]).astype(np.int64)
        la = np.array([mage.get(int(x), np.nan) for x in gid])
        in_cl = (np.isfinite(la) & (np.asarray(g["dg_ms"]) > -1.25) & (np.asarray(g["dg_ms"]) < 0.75) & np.isfinite(P)
                 & (sig > SIG_MIN_MM) & single_lined)
        clusters = {}
        # Cluster binaries bright enough for Gaia's spectroscopy are mostly hot
        # (median 6,700 K in clusters < 300 Myr), and hotter stars circularize
        # less far; so the comparison with the field uses the same 5300-6300 K
        # window as the kinematic populations. Older clusters have too few.
        for name, (lo, hi), tlim in (("< 1 Gyr, 5300-6300 K", (0, 9.0), (5300, 6300)),
                                     ("< 1 Gyr, all temperatures", (0, 9.0), (0, 99999)),
                                     ("> 1 Gyr, all temperatures", (9.0, 11), (0, 99999))):
            m = in_cl & (la >= lo) & (la < hi) & (P < 100) & (teff >= tlim[0]) & (teff < tlim[1])
            r = _pcirc_with_errors(P[m], e_ls[m])
            r["median_age_gyr"] = float(10 ** np.median(la[m]) / 1e9) if m.any() else np.nan
            r["median_teff"] = float(np.median(teff[m])) if m.any() else np.nan
            r["n_all"] = int((in_cl & (la >= lo) & (la < hi)).sum())
            clusters[name] = r
        out["clusters"] = clusters
        out["cluster_members"] = {"P": P[in_cl], "e": e[in_cl], "logage": la[in_cl], "circ": circ[in_cl]}
        return out
    return _cached("gaia_circ", {"v": 7, "sig": [SIG_MIN_MM, SIG_MIN_ENVELOPE]}, compute)


# ---------------------------------------------------------------------------
# 2. Wide pairs: Gaia's resolution, and the twin excess
# ---------------------------------------------------------------------------

def contrast_curves() -> dict:
    pairs = samples.wide_pairs()
    return {"all": wide.fit_contrast_curve(pairs), "colour": wide.fit_contrast_curve(wide.colour_subset(pairs))}


TWIN_MIN_THETA = 3.5  # arcsec: beyond Gaia's colour-blending zone for Delta G < 8


def twin_grid(curves: dict) -> dict:
    def compute():
        p = wide.ms_pairs()
        p = p[np.asarray(p["theta_arcsec"]) >= TWIN_MIN_THETA]
        out = {}
        for i, (m_lo, m_hi) in enumerate(wide.TWIN_MASS_BINS):
            for j, (s_lo, s_hi) in enumerate(wide.TWIN_SEP_BINS):
                m = ((p["m1"] >= m_lo) & (p["m1"] < m_hi) & (p["sep_au"] >= s_lo) & (p["sep_au"] < s_hi))
                sub = p[m]
                if len(sub) < 40:
                    out[(i, j)] = {"n": int(len(sub))}
                    continue
                S = wide.selection_matrix(sub, curves["colour"])
                q_break = 0.3 if m_lo >= 1.2 else 0.5  # as El-Badry et al. 2019
                r = wide.fit_twins(sub, S, q_break=q_break, seed=1000 + 10 * i + j)
                r.pop("chain", None)
                out[(i, j)] = r
        # q histograms for the figure: solar-type primaries, three separation ranges
        hists = {}
        for (s_lo, s_hi) in [(50, 600), (600, 2500), (2500, 15000)]:
            m = (p["m1"] >= 0.8) & (p["m1"] < 1.2) & (p["sep_au"] >= s_lo) & (p["sep_au"] < s_hi)
            h, e = np.histogram(np.asarray(p["q"][m]), bins=np.arange(0.3, 1.0001, 0.025))
            hists[(s_lo, s_hi)] = (h, e)
        return {"fits": out, "hists": hists, "n_pairs": int(len(p))}
    return _cached("twins", {"v": 3, "theta": TWIN_MIN_THETA}, compute)


# ---------------------------------------------------------------------------
# 3. The 100 pc census: samples and detection windows
# ---------------------------------------------------------------------------

def solar_sample():
    """Raghavan-like F6-K3 dwarfs within 100 pc, one row per system: stars
    whose brighter resolved companion is also in the sample are dropped (as
    Raghavan et al. dropped nine)."""
    t = samples.load_100pc()
    sol = samples.solar_type_mask(t)
    ids_in = set(np.asarray(t["source_id"])[sol].tolist())
    pairs = samples.wide_pairs()
    sec_of_member = set(int(b) for a, b in zip(np.asarray(pairs["source_id1"]), np.asarray(pairs["source_id2"]))
                        if int(a) in ids_in)
    keep = sol & ~np.isin(np.asarray(t["source_id"]), list(sec_of_member))
    return t[keep]


N_PER_BIN = 300


def mock_census(sample_name: str, prim_table, *, q_model="mds", e_model="mds", n_per_bin=N_PER_BIN,
                workers=6, edges=None) -> dict:
    """Mock companions around a sample's stars, with every channel's verdict."""
    import completeness as CP
    import mocks
    if edges is None:
        edges = CP.LOGP_EDGES

    def compute():
        rng = np.random.default_rng(C.SEED + zlib.crc32(sample_name.encode()) % 1000)  # stable across runs
        prim = CP.primary_arrays(prim_table)
        mk = mocks.make_mocks(rng, n_per_bin, edges, prim["m1"], q_model=q_model, e_model=e_model)
        f_g = mocks.flux_ratios(mk, np.zeros(5))[0]
        print(f"    {sample_name}: {len(mk):,} mock companions; gaiamock on log P < {CP.LOGP_ASTROMETRY[1]}",
              flush=True)
        ast = CP.astrometric_outcomes(mk, prim, f_g, seed=C.SEED, workers=workers)
        curve = contrast_curves()["all"]
        p_wide = CP.resolved_probability(mk, prim, curve, rng)
        return {"mock": mk, "ast": ast, "p_wide": p_wide, "f_g": f_g, "edges": edges,
                "n_primaries": len(prim_table)}

    key = {"v": 3, "sample": sample_name, "n": len(prim_table), "q": q_model, "e": e_model, "npb": n_per_bin,
           "gaiamock": C.GAIAMOCK_COMMIT, "edges": [float(edges[0]), float(edges[-1]), len(edges)]}
    if edges is CP.LOGP_EDGES or np.array_equal(edges, CP.LOGP_EDGES):
        key.pop("edges")  # keep the cache key of the original full-range runs
    return _cached(f"mocks_{sample_name}_{q_model}_{e_model}", key, compute)


#: Solution types whose acceptance rules the completeness model reproduces.
#: Targeted-search and "alternative" astrometric orbits used different, relaxed
#: rules for pre-selected candidates and are left out of the corrected census.
CENSUS_SPECTRO = ("SB1", "SB1C", "SB2", "SB2C", "AstroSpectroSB1", "EclipsingSpectro")
CENSUS_ASTRO = ("Orbital",)


def observed_companions(sample) -> dict:
    """Every companion Gaia reveals around the sample's stars, by channel.

    Orbits carry their DR3 period. Resolved pairs get one the way Raghavan et
    al. estimated theirs: a = 10^0.13 x projected separation (Duquennoy & Mayor
    1991), then Kepler's law with both masses (from absolute magnitudes on the
    dwarf sequence; a white-dwarf companion is given 0.6 Msun).
    """
    import binaries as B
    import stellar
    ids = np.asarray(sample["source_id"], dtype=np.int64)
    m1 = np.asarray(sample["mass_colour"], dtype=float)
    types = np.asarray(sample["nss_type"])
    per = np.asarray(sample["nss_period"], dtype=float)
    spec = np.isin(types, CENSUS_SPECTRO)
    ast = np.isin(types, CENSUS_ASTRO)
    ecl = types == "EclipsingBinary"
    acc = (np.asarray(sample["accel_type"]) != "") | (np.asarray(sample["trend_type"]) != "")
    ruwe = np.ma.filled(np.ma.asarray(sample["ruwe"], dtype=float), np.nan)
    other_nss = (types != "") & ~spec & ~ast & ~ecl
    ruwe_only = (ruwe > 1.4) & ~spec & ~ast & ~ecl & ~acc & ~other_nss

    pairs = samples.wide_pairs()
    index = {s: i for i, s in enumerate(ids)}
    w_logp, w_star, w_sep, w_q = [], [], [], []
    for k in range(len(pairs)):
        for me, other_absg, other_bt in ((pairs["source_id1"][k], pairs["abs_g2"][k], 2),
                                         (pairs["source_id2"][k], pairs["abs_g1"][k], 1)):
            i = index.get(int(me))
            if i is None:
                continue
            bt = str(pairs["bin_type"][k])
            other_is_wd = (bt[2:4] == "WD") if other_bt == 2 else (bt[:2] == "WD")
            m2 = 0.6 if other_is_wd else float(stellar.mass_from_mg(other_absg))
            if not np.isfinite(m2):
                m2 = 0.08  # fainter than the M9.5 end of the sequence: a brown dwarf or very late M
            a = L.RAGHAVAN_A_OVER_S * float(pairs["sep_au"][k])
            w_logp.append(np.log10(B.period_days_from_a(a, m1[i] + m2)))
            w_star.append(i)
            w_sep.append(float(pairs["sep_au"][k]))
            w_q.append(m2 / m1[i])
    return {
        "n_primaries": len(sample),
        "spectro_logp": np.log10(per[spec]), "astro_logp": np.log10(per[ast]), "ecl_logp": np.log10(per[ecl]),
        "n_spectro": int(spec.sum()), "n_astro": int(ast.sum()), "n_ecl": int(ecl.sum()),
        "n_accel": int((acc & ~spec & ~ast & ~ecl).sum()), "n_ruwe_only": int(ruwe_only.sum()),
        "n_other_nss": int(other_nss.sum()),
        "wide_logp": np.array(w_logp), "wide_star": np.array(w_star, int), "wide_sep": np.array(w_sep),
        "wide_q": np.array(w_q),
        "any_flag": int(np.sum(spec | ast | ecl | acc | other_nss | ruwe_only |
                               np.isin(np.arange(len(ids)), np.array(w_star, int)))),
    }


def rv_models() -> dict:
    """The spectroscopic detection model's ingredients: per-transit noise,
    the RVS-band colour relation, and the efficiency calibrated on SB9."""
    import spectro
    t = samples.load_100pc()
    f = lambda c: np.ma.filled(np.ma.asarray(t[c], dtype=float), np.nan)
    noise = spectro.NoiseModel(f("grvs_mag"), f("radial_velocity_error"), f("rv_nb_transits"),
                               f("rv_renormalised_gof"), f("rv_chisq_pvalue"), f("rv_template_teff"))
    coeffs = samples.fit_grvs_colour(t)

    def validate():
        return spectro.validate_on_sb9(sb9.load_sb9(), samples.gaia_table("nss_orbits"), noise, n_sim=40)
    val = _cached("sb9_validation", {"v": 2}, validate)
    eff = spectro.efficiency_from_validation(val, np.asarray(sb9.main_sequence_mask(val["table"])))
    return {"noise": noise, "grvs_coeffs": coeffs, "validation": val, "efficiency": eff}


def intrinsic_models(logp, m1=1.0) -> dict:
    """Companions per primary per decade of period under the published fits."""
    import mocks
    gauss = lambda x, mu, s: np.exp(-0.5 * ((x - mu) / s) ** 2) / (s * np.sqrt(2 * np.pi))
    # Moe & Di Stefano's q > 0.3 fit scaled to q > 0.1 with their own q-shape
    rng = np.random.default_rng(1)
    q_ratio = []
    for lp in logp:
        qs = np.array([mocks._draw_q(rng, float(L.mds_gamma_smallq(np.array([lp]), m1)[0]),
                                     float(L.mds_gamma_largeq(np.array([lp]), m1)[0]),
                                     float(L.mds_ftwin(np.array([lp]), m1)[0])) for _ in range(4000)])
        q_ratio.append(1.0 / np.mean(qs > 0.3))
    mds = np.nan_to_num(L.mds_flogp_q03(logp, m1)) * np.array(q_ratio)
    return {
        "raghavan": L.RAGHAVAN_N_COMPANIONS_FIG13 / L.RAGHAVAN_N_PRIMARIES
        * gauss(logp, L.RAGHAVAN_LOGP_MEAN, L.RAGHAVAN_LOGP_SIGMA),
        "tokovinin": L.TOKOVININ_ALL["eps"] * gauss(logp, L.TOKOVININ_ALL["x0"], L.TOKOVININ_ALL["sigma"]),
        "mds": mds,
    }


def solar_census(q_model="mds", e_model="mds") -> dict:
    """The period distribution of companions to Sun-like stars within 100 pc."""
    import completeness as CP
    sol = solar_sample()
    rv = rv_models()
    mc = mock_census("solar100", sol, q_model=q_model, e_model=e_model)
    prim = CP.primary_arrays(sol)

    def spec_fn():
        return CP.spectroscopic_outcomes(mc["mock"], prim, rv["noise"], rv["grvs_coeffs"], rv["efficiency"],
                                         np.random.default_rng(C.SEED + 11))
    spec = _cached(f"spec_solar100_{q_model}_{e_model}", {"v": 2, "n": len(sol), "eff": rv["efficiency"]["eff"].tolist(),
                                                          "sb2": rv["efficiency"]["eff_sb2"]}, spec_fn)
    comp = CP.channel_completeness(mc, spec, efficiency="sb9")
    comp_ideal = CP.channel_completeness(mc, spec, efficiency="ideal")
    obs = observed_companions(sol)
    edges = comp["edges"]
    width = np.diff(edges)
    n = obs["n_primaries"]
    # observed companions with a period, by bin (wide pairs: only q > 0.1, as in the model)
    keep_w = obs["wide_q"] >= 0.1
    h_orbit = np.histogram(np.concatenate([obs["spectro_logp"], obs["astro_logp"], obs["ecl_logp"]]), edges)[0]
    h_spec = np.histogram(obs["spectro_logp"], edges)[0]
    h_ast = np.histogram(obs["astro_logp"], edges)[0]
    h_wide = np.histogram(obs["wide_logp"][keep_w], edges)[0]
    h_all = h_orbit + h_wide
    def corrected(cmp):
        c = cmp["period_known"]
        ok = c >= 0.1
        f = np.where(ok, h_all / (n * np.where(ok, c, 1) * width), np.nan)
        e = np.where(ok, np.sqrt(np.maximum(h_all, 1)) / (n * np.where(ok, c, 1) * width), np.nan)
        return f, e
    f_corr, f_err = corrected(comp)
    f_ideal, f_ideal_err = corrected(comp_ideal)
    models = intrinsic_models(comp["centres"])
    # forward model: predicted counts per channel under each published distribution
    pred, pred_ideal = {}, {}
    for name, f in models.items():
        ff = np.where(comp["centres"] < 8.0, f, 0.0) if name == "mds" else f
        for dst, cmp in ((pred, comp), (pred_ideal, comp_ideal)):
            dst[name] = {ch: float(np.nansum(n * ff * width * cmp[ch]))
                         for ch in ("orbit", "spectroscopic", "astrometric", "acceleration", "rv_trend", "wide",
                                    "ruwe_only")}
    observed_counts = {"orbit": int(h_orbit.sum()), "spectroscopic": int(h_spec.sum()), "astrometric": int(h_ast.sum()),
                       "wide": int(h_wide.sum()), "acceleration": None, "ruwe_only": obs["n_ruwe_only"]}
    acc_ast = int(np.sum((np.asarray(sol["accel_type"]) != "") & (np.asarray(sol["nss_type"]) == "")))
    trend = int(np.sum((np.asarray(sol["trend_type"]) != "") & (np.asarray(sol["nss_type"]) == "")
                       & (np.asarray(sol["accel_type"]) == "")))
    observed_counts["acceleration"] = acc_ast
    observed_counts["rv_trend"] = trend
    return {"n_primaries": n, "completeness": {k: v for k, v in comp.items() if k != "per_mock"},
            "completeness_ideal": {k: v for k, v in comp_ideal.items() if k != "per_mock"},
            "observed": obs, "h_orbit": h_orbit, "h_spec": h_spec, "h_ast": h_ast, "h_wide": h_wide,
            "f_corr": f_corr, "f_err": f_err, "f_ideal": f_ideal, "f_ideal_err": f_ideal_err,
            "models": models, "predicted": pred, "predicted_ideal": pred_ideal,
            "observed_counts": observed_counts, "q_model": q_model, "e_model": e_model}


#: Close companions (P < 1000 d) against primary mass, within 100 pc.
CLOSE_MASS_BINS = [("M0-M4", 0.25, 0.60), ("K", 0.60, 0.80), ("G", 0.80, 1.05), ("F", 1.05, 1.50),
                   ("A", 1.50, 2.50)]
CLOSE_EDGES = np.arange(-0.5, 3.51, 0.25)


def mass_bin_sample(m_lo, m_hi):
    t = samples.load_100pc()
    keep = (t["good_astrometry"] & t["ms_band"] & (t["wide_role"] != "secondary")
            & (np.asarray(t["mass_colour"]) >= m_lo) & (np.asarray(t["mass_colour"]) < m_hi))
    return t[keep]


def close_census() -> list[dict]:
    """Companions per star with 1 d < P < 1000 d (q > 0.1), by primary mass:
    Gaia's orbits corrected with the same detection model as the solar census."""
    import completeness as CP
    rv = rv_models()
    out = []
    for name, lo, hi in CLOSE_MASS_BINS:
        smp = mass_bin_sample(lo, hi)
        mc = mock_census(f"close_{name}", smp, n_per_bin=200, edges=CLOSE_EDGES)
        prim = CP.primary_arrays(smp)

        def spec_fn():
            return CP.spectroscopic_outcomes(mc["mock"], prim, rv["noise"], rv["grvs_coeffs"], rv["efficiency"],
                                             np.random.default_rng(C.SEED + 21))
        spec = _cached(f"spec_close_{name}", {"v": 2, "n": len(smp), "eff": rv["efficiency"]["eff"].tolist(),
                                              "sb2": rv["efficiency"]["eff_sb2"]}, spec_fn)
        types = np.asarray(smp["nss_type"])
        per = np.asarray(smp["nss_period"], dtype=float)
        is_orbit = np.isin(types, CENSUS_SPECTRO + CENSUS_ASTRO + ("EclipsingBinary",))
        row = {"name": name, "m_lo": lo, "m_hi": hi, "n_primaries": len(smp),
               "m_med": float(np.median(np.asarray(smp["mass_colour"])))}
        for eff_mode in ("sb9", "ideal"):
            comp = CP.channel_completeness(mc, spec, efficiency=eff_mode)
            edges = comp["edges"]
            h = np.histogram(np.log10(per[is_orbit]), edges)[0]
            sel = (comp["centres"] > 0.0) & (comp["centres"] < 3.0)
            c = comp["orbit"]
            f = h[sel] / (len(smp) * np.maximum(c[sel], 1e-3))
            row[f"close_{eff_mode}"] = float(f.sum())
            row[f"close_{eff_mode}_err"] = float(np.sqrt(np.sum(h[sel] / (len(smp) * np.maximum(c[sel], 1e-3)) ** 2)))
            row[f"completeness_{eff_mode}"] = float(np.mean(c[sel]))
            row["n_orbits"] = int(h[sel].sum())
        x = np.linspace(0.2, 3.0, 400)
        q_fac = intrinsic_models(np.array([1.5]), m1=max(row["m_med"], 0.8))["mds"][0] / \
            max(float(L.mds_flogp_q03(np.array([1.5]), max(row["m_med"], 0.8))[0]), 1e-9)
        row["mds_close"] = float(np.trapezoid(np.nan_to_num(L.mds_flogp_q03(x, max(row["m_med"], 0.8))), x) * q_fac)
        out.append(row)
    return out


# ---------------------------------------------------------------------------
# 4. Gaia's binaries on the HR diagram
# ---------------------------------------------------------------------------

HRD_CHANNELS = ("none", "ruwe", "acceleration", "astrometric", "sb1", "sb2")


def hr_diagram() -> dict:
    import hrd
    t = samples.load_100pc()
    ruwe = np.ma.filled(np.ma.asarray(t["ruwe"], dtype=float), np.nan)
    good = t["good_astrometry"] & np.isfinite(t["bp_rp"]) & np.isfinite(t["abs_g"])
    single = good & (ruwe < 1.4) & (t["channel"] == "") & (t["n_wide"] == 0) & t["ms_band"]
    centres, ridge = hrd.empirical_ridge(t["bp_rp"], t["abs_g"], single)
    delta = hrd.offsets(t, centres, ridge)
    lab = hrd.channel_labels(t)
    out = {"bp_rp": np.asarray(t["bp_rp"])[good], "abs_g": np.asarray(t["abs_g"])[good],
           "label": lab[good], "delta": np.asarray(delta)[good], "ms_band": np.asarray(t["ms_band"])[good],
           "ridge_c": centres, "ridge_g": ridge, "stats": {}}
    for name, (lo, hi) in (("FGK", (0.6, 1.3)), ("KM", (1.0, 1.8))):
        sel = out["ms_band"] & (out["bp_rp"] >= lo) & (out["bp_rp"] < hi)
        st = {}
        for ch in HRD_CHANNELS:
            m = sel & (out["label"] == ch)
            d = out["delta"][m]
            d = d[np.isfinite(d)]
            st[ch] = {"n": int(m.sum()), "median": float(np.median(d)) if len(d) else np.nan,
                      "frac_twin_zone": float(np.mean((d > -0.9) & (d < -0.6))) if len(d) else np.nan,
                      "frac_elevated": float(np.mean(d < -0.3)) if len(d) else np.nan}
        tz = sel & (out["delta"] > -0.9) & (out["delta"] < -0.6)
        st["twin_zone_composition"] = {ch: float(np.mean(out["label"][tz] == ch)) for ch in HRD_CHANNELS}
        st["twin_zone_n"] = int(tz.sum())
        out["stats"][name] = st
    return out


# ---------------------------------------------------------------------------
# 5. DR3's non-single-star catalogue as a whole
# ---------------------------------------------------------------------------

#: Arenou et al. (2023), Table 1: solutions per type in Gaia DR3.
ARENOU23_COUNTS = {"Acceleration7": 246_947, "Acceleration9": 91_268, "Orbital": 134_598,
                   "AstroSpectroSB1": 33_467, "SB1 or SB2": 186_905, "EclipsingSpectro": 155,
                   "EclipsingBinary": 86_918, "FirstDegreeTrendSB1": 24_083, "SecondDegreeTrendSB1": 32_725,
                   "VIMF": 870}


def nss_counts() -> dict:
    """Solutions per type in our downloads, against Arenou et al.'s Table 1."""
    out = {}
    for name in ("nss_orbits", "nss_accel", "nss_trend", "nss_vim"):
        t = samples.gaia_table(name)
        types, n = np.unique(np.char.strip(np.asarray(t["nss_solution_type"]).astype(str)), return_counts=True)
        out.update(dict(zip(types.tolist(), n.tolist())))
    ours = {k: out.get(k, 0) for k in ARENOU23_COUNTS if k != "SB1 or SB2"}
    ours["SB1 or SB2"] = sum(out.get(k, 0) for k in ("SB1", "SB1C", "SB2", "SB2C"))
    sources = set()
    for name in ("nss_orbits", "nss_accel", "nss_trend", "nss_vim"):
        sources.update(np.asarray(samples.gaia_table(name)["source_id"]).tolist())
    return {"by_type": out, "compared": {k: (ours.get(k), v) for k, v in ARENOU23_COUNTS.items()},
            "unique_sources": len(sources)}


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def run() -> dict:
    import multiplicity
    res = {}
    print("  DR3 non-single-star counts", flush=True)
    res["counts"] = nss_counts()
    print("  Raghavan et al. 2010 rebuilt from their Table 13", flush=True)
    res["raghavan"] = raghavan.refit()
    print("  SB9: period-eccentricity and circularization", flush=True)
    res["sb9"] = sb9_circularization()
    print("  Gaia DR3: circularization by temperature, population and cluster age", flush=True)
    res["gaia_circ"] = gaia_circularization()
    print("  Gaia's resolution for close pairs (contrast curves)", flush=True)
    res["contrast"] = contrast_curves()
    print("  Twin excess among wide pairs (MCMC per bin)", flush=True)
    res["twins"] = twin_grid(res["contrast"])
    print("  HR diagram of the 100 pc sample", flush=True)
    res["hrd"] = hr_diagram()
    print("  Radial-velocity detection model, calibrated on SB9", flush=True)
    rv = rv_models()
    res["rv_validation"] = {k: v for k, v in rv["validation"].items() if k != "table"}
    res["rv_efficiency"] = rv["efficiency"]
    res["rv_noise"] = {"grvs": rv["noise"].grvs, "sigma": rv["noise"].sigma}
    print("  The census of Sun-like stars within 100 pc (gaiamock + RV model)", flush=True)
    res["solar"] = solar_census()
    print("  ... systematics: flat mass ratios, thermal eccentricities", flush=True)
    res["solar_flatq"] = solar_census(q_model="flat")
    res["solar_thermal"] = solar_census(e_model="thermal")
    print("  Close companions against primary mass", flush=True)
    res["close"] = close_census()
    print("  Wide-binary fraction against primary mass", flush=True)
    res["wbf"] = multiplicity.wide_binary_fraction(np.random.default_rng(C.SEED + 3), res["contrast"]["all"])
    return res


def _integral(sc, key, lo, hi):
    """Companions per star between log P = lo and hi from a census' corrected bins."""
    c = sc["completeness"]["centres"]
    w = np.diff(sc["completeness"]["edges"])
    m = (c > lo) & (c < hi)
    f = np.asarray(sc[key])[m]
    if not np.all(np.isfinite(f)):
        return float("nan")
    return float(np.sum(f * w[m]))


def summarize(res: dict) -> dict:
    out = {}
    out["dr3_counts"] = res["counts"]["compared"]
    out["dr3_unique_nss_sources"] = res["counts"]["unique_sources"]
    r = res["raghavan"]
    out["raghavan_refit"] = {"hist_mean_logP": r["hist_mean"], "hist_sigma": r["hist_sigma"],
                             "ml_mean_logP": r["ml_mean"], "ml_sigma": r["ml_sigma"],
                             "n_companions": r["n_companions"], "published": [L.RAGHAVAN_LOGP_MEAN,
                                                                             L.RAGHAVAN_LOGP_SIGMA]}
    s = res["sb9"]
    out["sb9_fgk_dwarfs"] = {"n": s["n_fgk"], "p_circ_meibom_mathieu": s["p_circ"], "p_circ_err": s["p_circ_err"],
                             "half_circular_period": s["circ_frac"]["p50"]}
    g = res["gaia_circ"]
    out["gaia_circularization"] = {
        "n": g["n"], "p_circ_all": g["all"], "half_circular_period": g["circ_frac"]["p50"],
        "teff_bins": [{k: v for k, v in b.items()} for b in g["teff_bins"]],
        "teff_bins_envelope": [{k: v for k, v in b.items()} for b in g["teff_bins_envelope"]],
        "populations": {k: {kk: vv for kk, vv in v.items() if kk != "circ_frac"} for k, v in g["populations"].items()},
        "clusters": g["clusters"]}
    out["twins"] = {f"M{wide.TWIN_MASS_BINS[i]}_s{wide.TWIN_SEP_BINS[j]}": f for (i, j), f in res["twins"]["fits"].items()}
    out["hrd"] = res["hrd"]["stats"]
    out["rv_validation_sb9"] = res["rv_validation"]
    out["rv_efficiency"] = {"eff": res["rv_efficiency"]["eff"], "n": res["rv_efficiency"]["n"],
                            "eff_sb2": res["rv_efficiency"]["eff_sb2"]}
    sc = res["solar"]
    obs = sc["observed"]
    out["solar_census"] = {
        "n_primaries": sc["n_primaries"], "any_gaia_sign": obs["any_flag"],
        "any_gaia_sign_fraction": obs["any_flag"] / sc["n_primaries"],
        "observed_counts": sc["observed_counts"], "predicted": sc["predicted"], "predicted_ideal": sc["predicted_ideal"],
        "companions_per_star_logP_0_3": [_integral(sc, "f_corr", 0, 3), _integral(sc, "f_ideal", 0, 3)],
        "companions_per_star_logP_2_3": [_integral(sc, "f_corr", 2, 3), _integral(sc, "f_ideal", 2, 3)],
        "companions_per_star_logP_5_9": _integral(sc, "f_corr", 5.0, 9.0),
        "systematics_flat_q": {"logP_0_3": _integral(res["solar_flatq"], "f_corr", 0, 3),
                               "logP_5_9": _integral(res["solar_flatq"], "f_corr", 5.0, 9.0),
                               "predicted_mds": res["solar_flatq"]["predicted"]["mds"]},
        "systematics_thermal_e": {"logP_0_3": _integral(res["solar_thermal"], "f_corr", 0, 3),
                                  "logP_5_9": _integral(res["solar_thermal"], "f_corr", 5.0, 9.0),
                                  "predicted_mds": res["solar_thermal"]["predicted"]["mds"]},
    }
    # companions per star over all periods: measured where Gaia measures periods,
    # and Moe & Di Stefano's value in the gap (where the forward model of
    # accelerations and RUWE excess is consistent with it)
    x = sc["completeness"]["centres"]
    w = np.diff(sc["completeness"]["edges"])
    measured = np.isfinite(sc["f_corr"])
    gap = ~measured & (x > 3) & (x < 5)
    for key in ("f_corr", "f_ideal"):
        meas = float(np.nansum((np.asarray(sc[key]) * w)[measured]))
        out["solar_census"][f"total_{key}"] = {"measured": meas,
                                               "gap_from_mds": float(np.sum((sc["models"]["mds"] * w)[gap])),
                                               "gap_logP": [float(x[gap].min() - 0.125), float(x[gap].max() + 0.125)]}
    for name, f in sc["models"].items():
        w = np.diff(sc["completeness"]["edges"])
        out["solar_census"][f"model_{name}_logP_0_3"] = float(np.sum((f * w)[(x > 0) & (x < 3)]))
        out["solar_census"][f"model_{name}_logP_5_9"] = float(np.sum((f * w)[(x > 5) & (x < 9)]))
        out["solar_census"][f"model_{name}_total"] = float(np.sum((f * w)[(x > -0.5) & (x < 9)]))
    out["close_vs_mass"] = res["close"]
    out["wide_binary_fraction"] = res["wbf"]
    return out


def print_summary(res: dict) -> None:
    s = summarize(res)
    sc = s["solar_census"]
    print(f"Sun-like stars within 100 pc: {sc['n_primaries']:,}; any Gaia sign of a companion: "
          f"{100 * sc['any_gaia_sign_fraction']:.1f}%")
    print("Observed vs predicted (Moe & Di Stefano):",
          {k: (v, round(sc['predicted']['mds'].get(k, np.nan))) for k, v in sc["observed_counts"].items()})
    print(f"Companions per star, 1-1000 d: {sc['companions_per_star_logP_0_3'][0]:.3f} (SB9-calibrated) / "
          f"{sc['companions_per_star_logP_0_3'][1]:.3f} (documented cuts); Moe & Di Stefano "
          f"{sc['model_mds_logP_0_3']:.3f}, Raghavan {sc['model_raghavan_logP_0_3']:.3f}")
    g = s["gaia_circularization"]
    print(f"Circularization: Gaia dwarfs P_circ = {g['p_circ_all']['p_circ']:.2f} d; SB9 FGK dwarfs "
          f"{s['sb9_fgk_dwarfs']['p_circ_meibom_mathieu']:.2f} d")
