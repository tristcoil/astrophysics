"""The science: one function per question, each returning (results, artefacts).

``results`` holds plain numbers for results.json and the README; ``artefacts``
holds the arrays the figures need. Nothing here draws anything.
"""

from __future__ import annotations

import json
import re
import warnings
from dataclasses import replace

import numpy as np
import pandas as pd

import rv_fit as fit
import rv_physics as phys
import rvdata
from config import (
    DISCOVERY_ANNOUNCED_BJD, FAP_THRESHOLD, KP_51PEG, MCMC_BURN, MCMC_STEPS,
    MCMC_THIN, MCMC_WALKERS_PER_DIM, MQ95, OVERSAMPLE,
    PERIOD_MIN_D, RESULTS_PATH, STARS, STELLAR_MASS,
)
from fetch_data import NASA_HOSTS_PATH, NASA_PS_PATH

warnings.filterwarnings("ignore", category=RuntimeWarning)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _plain(value) -> str:
    """Strip the HTML the archive wraps around reference names."""
    import html as _html
    if not isinstance(value, str):
        return ""
    return _html.unescape(re.sub(r"<[^>]+>", "", value)).strip()


def check_stellar_masses() -> dict:
    """Confirm every adopted stellar mass is a row of the archive's table."""
    hosts = pd.read_csv(NASA_HOSTS_PATH)
    hosts["ref"] = hosts["st_refname"].map(_plain)
    out = {}
    for key, (mass, err, ref) in STELLAR_MASS.items():
        rows = hosts[(hosts["hostname"] == STARS[key].archive_host) & (hosts["ref"] == ref)]
        ok = bool(len(rows)) and np.isclose(rows["st_mass"].iloc[0], mass) and \
            np.isclose(rows["st_masserr1"].iloc[0], err)
        if not ok:
            raise RuntimeError(f"adopted mass for {key} ({mass} +- {err}, {ref}) not in stellarhosts")
        others = hosts[(hosts["hostname"] == STARS[key].archive_host) & hosts["st_mass"].notna()]
        out[key] = {"mass": mass, "err": err, "ref": ref,
                    "archive_range": [float(others["st_mass"].min()), float(others["st_mass"].max())],
                    "n_published": int(len(others))}
    return out


def mass_draws(key: str, rng, n: int) -> np.ndarray:
    mass, err, _ = STELLAR_MASS[key]
    return np.abs(rng.normal(mass, err, n))


def published(planet: str) -> pd.DataFrame:
    """Every published solution for one planet in the archive's ps table."""
    ps = pd.read_csv(NASA_PS_PATH)
    ps = ps[ps["pl_name"] == planet].copy()
    ps["ref"] = ps["pl_refname"].map(_plain)
    year = ps["ref"].str.extract(r"(\d{4})")[0].astype(float)
    ps["year"] = year
    return ps.sort_values("year")


def literature_rows(planet: str, cols=("pl_orbper", "pl_rvamp", "pl_orbeccen", "pl_msinij", "pl_msinie")) -> list[dict]:
    rows = []
    for _, r in published(planet).iterrows():
        entry = {"ref": r["ref"], "year": None if np.isnan(r["year"]) else int(r["year"]),
                 "default": bool(r["default_flag"])}
        for c in cols + ("pl_orbpererr1", "pl_rvamperr1", "pl_msinieerr1", "st_mass"):
            v = r.get(c)
            entry[c] = None if v is None or (isinstance(v, float) and np.isnan(v)) else float(v)
        rows.append(entry)
    return rows


def run_posterior(res: fit.FitResult, tag: str, quiet: bool, *, steps=MCMC_STEPS, burn=MCMC_BURN) -> fit.Posterior:
    return fit.run_mcmc(res, tag=tag, n_steps=steps, burn=burn, thin=MCMC_THIN,
                        walkers_per_dim=MCMC_WALKERS_PER_DIM, seed=_seed(tag), quiet=quiet)


def _seed(tag: str) -> int:
    import hashlib
    return int(hashlib.sha256(tag.encode()).hexdigest()[:8], 16)


def convergence(post: fit.Posterior) -> dict:
    tau = post.autocorr
    return {
        "samples": int(len(post.samples)),
        "acceptance": round(post.acceptance, 3),
        "max_autocorr_steps": None if tau is None else float(np.max(tau)),
        "steps_over_tau": None if tau is None else float(post.n_steps / np.max(tau)),
    }


def planet_posteriors(post: fit.Posterior, star_key: str, rng, labels: list[str]) -> list[dict]:
    """Summaries of every planet's parameters, including m sin i with the
    stellar mass drawn from its uncertainty."""
    P, K, e, w, tc = post.planet_samples()
    mstar = mass_draws(star_key, rng, len(P))
    out = []
    for j, label in enumerate(labels):
        mj = phys.msini_mjup(P[:, j], K[:, j], e[:, j], mstar)
        me = phys.msini_mearth(P[:, j], K[:, j], e[:, j], mstar)
        a = phys.semi_major_axis_au(P[:, j], mstar * phys.GM_SUN, me * phys.GM_EARTH)
        entry = {
            "label": label,
            "P": fit.summarize(P[:, j]), "K": fit.summarize(K[:, j]),
            "tc": fit.summarize(tc[:, j]),
            "msini_mjup": fit.summarize(mj), "msini_mearth": fit.summarize(me),
            "a_au": fit.summarize(a),
        }
        if f"h{j + 1}" in post.model.names:
            entry["e"] = fit.summarize(e[:, j], upper_limit=0.95)
            entry["omega_deg"] = fit.summarize(np.degrees(w[:, j]))
        out.append(entry)
    return out


def lucy_sweeney(e_samples) -> dict:
    """Lucy & Sweeney (1971): the probability that a truly circular orbit,
    measured with this precision, would show an eccentricity at least this
    large. Measurement noise always pushes a fitted e above zero, so a small
    e is only real if this probability is small (< 5%)."""
    e_hat = float(np.median(e_samples))
    sigma = float(np.std(e_samples))
    p = float(np.exp(-e_hat ** 2 / (2 * sigma ** 2)))
    return {"e_hat": e_hat, "sigma_e": sigma, "p_circular": p, "significant": p < 0.05}


def best_peak(data, freqs, power, y, err, multi_offset=True) -> dict:
    top = phys.top_peaks(freqs, power, n=6)
    i = top[0]
    df = freqs[1] - freqs[0]
    f_ref, p_ref = phys.refine_peak(data.t, y, err, freqs[i], df, data.inst if multi_offset else None)
    fap = phys.fap_baluev(p_ref, data.t, err, freqs.max(), n_offsets=data.n_inst if multi_offset else 1)
    return {"period": float(1 / f_ref), "power": float(p_ref), "fap": float(fap),
            "top": [{"period": float(1 / freqs[j]), "power": float(power[j])} for j in top]}


def _start_from(res: fit.FitResult, model: fit.KeplerModel) -> np.ndarray:
    """A starting vector for ``model`` taking every shared parameter from a
    previous fit (e.g. eccentric from circular)."""
    th = model.x0.copy()
    for i, n in enumerate(model.names):
        if n in res.model.names:
            th[i] = res.theta[res.model.names.index(n)]
    return th


def fit_eccentric_from(circ: fit.FitResult, *, trend=None, activity=None, which=None) -> fit.FitResult:
    """Refit with eccentric orbits (all planets, or those in ``which``),
    starting from a circular solution and from several eccentricity vectors;
    the nearly circular start guarantees the result is at least as good."""
    base = circ.planets()
    specs = circ.model.planets

    def factory(start):
        planets = []
        for j, (p, q) in enumerate(zip(specs, base)):
            ecc = which is None or j in which
            h, k = start if ecc else (0.0, 0.0)
            planets.append(replace(p, P=q["P"], K=q["K"], tc=q["tc"], h=h, k=k, circular=not ecc))
        model = fit.KeplerModel(circ.model.data, planets, trend=circ.model.trend if trend is None else trend,
                                activity=circ.model.activity if activity is None else activity,
                                t_ref=circ.model.t_ref, fixed_jitter=circ.model.fixed_jitter,
                                night_term=circ.model.night_insts)
        return model, _start_from(circ, model)

    starts = [(1e-3, 1e-3), (0.15, 0.0), (-0.15, 0.0), (0.0, 0.15), (0.0, -0.15)]
    return fit.fit_best_of(factory, starts)


def instrument_table(data: fit.RVData, res: fit.FitResult) -> list[dict]:
    """Per-instrument counts, time span, fitted zero point and noise."""
    lin, lin_err = res.model.linear_params(res.theta)
    p = res.params()
    resid = res.model.residuals(res.theta)
    from astropy.time import Time
    out = []
    for k, name in enumerate(data.inst_names):
        m = data.inst == k
        years = Time(data.t[m], format="jd").decimalyear
        out.append({
            "instrument": name, "n": int(m.sum()),
            "nights": int(np.unique(np.floor(data.t[m] - 0.3)).size),
            "first": float(years.min()), "last": float(years.max()),
            "median_error": float(np.median(data.err[m])),
            "offset": float(lin[f"gamma[{name}]"]), "offset_err": float(lin_err[f"gamma[{name}]"]),
            "jitter": float(p.get(f"jit[{name}]", np.nan)),
            "night_noise": float(p[f"night[{name}]"]) if f"night[{name}]" in p else None,
            "rms": float(np.std(resid[m])),
        })
    return out


def alias_safe_bounds(P: float, fallback=None) -> tuple[float, float]:
    """Prior range for a period: its own periodogram peak, but not the
    one-year aliases on either side (+-1/365 cycles per day away). Without
    this a sampler can hop between a weak planet's period and its yearly
    alias, which the data only mildly disfavour. For periods so long that the
    rule is undefined, ``fallback`` (or a factor of two) is used."""
    df = 0.5 / phys.YEAR_D
    f = 1.0 / P
    if f - df <= 0:
        return fallback if fallback else (P / 2, P * 2)
    return (1.0 / (f + df), 1.0 / (f - df))


def aliases_of(period: float) -> dict:
    f = 1 / period
    f_sid = 1 / phys.SIDEREAL_DAY_D
    f_yr = 1 / phys.YEAR_D
    return {
        "sidereal_day": [float(1 / abs(f_sid - f)), float(1 / (f_sid + f))],
        "solar_day": [float(1 / abs(1.0 - f)), float(1 / (1.0 + f))],
        "year": [float(1 / (f + f_yr)), float(1 / abs(f - f_yr))],
    }


# ---------------------------------------------------------------------------
# 51 Pegasi b
# ---------------------------------------------------------------------------

def analyse_51peg(rng, *, quiet: bool = False):
    """The planet that started it all, from four spectrographs over 24 years."""
    if not quiet:
        print("51 Pegasi")
    data, log = rvdata.build("51peg")
    res: dict = {"cuts": log.as_dict(), "counts": data.counts(), "n": data.n,
                 "baseline_days": data.baseline}
    art: dict = {"data": data}

    # --- 1. What Mayor & Queloz had: ELODIE before the announcement ---------
    elo = data.only("ELODIE")
    disc = elo.subset(elo.t < DISCOVERY_ANNOUNCED_BJD)
    f_disc = phys.frequency_grid(disc.t, PERIOD_MIN_D, 30.0, 20)
    flat = np.ones(disc.n)  # the archive gives no per-spectrum errors: equal weights
    p_disc = phys.gls(disc.t, disc.y, flat, f_disc)
    peak_disc = best_peak(disc, f_disc, p_disc, disc.y, flat, multi_offset=False)
    pw, amp, ph = phys.gls(disc.t, disc.y, flat, np.array([1 / peak_disc["period"]]), return_fit=True)
    spec = fit.planet_from_sinusoid(peak_disc["period"], amp[0], ph[0], np.median(disc.t))
    spec.P_bounds = (peak_disc["period"] * 0.97, peak_disc["period"] * 1.03)
    disc_fit = fit.fit_map(fit.KeplerModel(disc, [spec], night_term=[]))
    disc_err, _ = fit.laplace_errors(disc_fit.model, disc_fit.theta)
    dq = disc_fit.planets()[0]
    from astropy.time import Time
    res["discovery_season"] = {
        "n": disc.n, "first": Time(disc.t.min(), format="jd").iso[:10],
        "last": Time(disc.t.max(), format="jd").iso[:10],
        "nights": int(np.unique(np.floor(disc.t - 0.3)).size),
        "peak": peak_disc,
        "P": dq["P"], "P_err": float(disc_err[disc_fit.model.names.index("P1")]),
        "K": dq["K"], "K_err": float(disc_err[disc_fit.model.names.index("K1")]),
        "rms": disc_fit.rms(), "mayor_queloz_1995": MQ95,
    }
    art["discovery"] = {"data": disc, "freqs": f_disc, "power": p_disc, "fit": disc_fit}
    if not quiet:
        print(f"  discovery season: {disc.n} ELODIE spectra -> P = {dq['P']:.4f} d, "
              f"K = {dq['K']:.1f} m/s, rms {disc_fit.rms():.1f} m/s")

    # --- 2. Blind search of everything ------------------------------------------
    # ELODIE's archival velocities carry a seasonal systematic that tracks the
    # barycentric correction (see README); one ELODIE-only regression term
    # absorbs it. It is zero for every other instrument.
    elo_k = data.inst_names.index("ELODIE")
    data.activity["elodie_berv"] = np.where(data.inst == elo_k, data.activity["BERV"], 0.0)
    freqs = phys.frequency_grid(data.t, PERIOD_MIN_D, None, OVERSAMPLE)
    steps, circ = fit.search_planets(data, freqs, n_max=1, fap_threshold=FAP_THRESHOLD,
                                     activity=("elodie_berv",), quiet=quiet)
    first = steps[0]
    res["search"] = [{"planets_before": s.n_planets_before, "peak_period": s.peak_period,
                      "peak_power": s.peak_power, "fap": s.fap, "rms_before": s.rms_before,
                      "top": list(zip(s.top_periods, s.top_powers))} for s in steps]
    res["aliases_predicted"] = aliases_of(first.peak_period)
    window = phys.spectral_window(data.t, freqs)
    art["search"] = steps
    art["window"] = window
    art["freqs"] = freqs
    # Is the next-highest residual peak another planet? Its significance with
    # nights counted once, and whether every instrument sees the same signal.
    nxt = steps[-1]
    res["next_signal"] = {"period": nxt.peak_period, "power": nxt.peak_power, "fap": nxt.fap,
                          "significant": nxt.fap < FAP_THRESHOLD,
                          "by_instrument": fit.signal_by_instrument(circ, nxt.peak_period)}
    no_berv = fit.fit_map(fit.KeplerModel(data, circ.model.planets, night_term=circ.model.night_insts),
                          circ.theta)
    lin, lin_err = circ.model.linear_params(circ.theta)
    res["elodie_berv_term"] = {"slope_ms_per_kms": lin["c[elodie_berv]"],
                               "slope_err": lin_err["c[elodie_berv]"],
                               "delta_log_like": circ.log_like - no_berv.log_like,
                               "elodie_jitter_with": circ.params()["jit[ELODIE]"],
                               "elodie_jitter_without": no_berv.params()["jit[ELODIE]"]}

    # --- 3. Circular and eccentric orbits, MAP and posterior --------------------
    ecc = fit_eccentric_from(circ)
    post_c = run_posterior(circ, "51peg_circ", quiet)
    post_e = run_posterior(ecc, "51peg_ecc", quiet)
    _, _, e_s, w_s, _ = post_e.planet_samples()
    ls = lucy_sweeney(e_s[:, 0])
    res["fits"] = {
        "circular": {"log_like": circ.log_like, "bic": circ.bic, "rms": circ.rms(),
                     "planet": planet_posteriors(post_c, "51peg", rng, ["51 Peg b"])[0],
                     "convergence": convergence(post_c)},
        "eccentric": {"log_like": ecc.log_like, "bic": ecc.bic, "rms": ecc.rms(),
                      "planet": planet_posteriors(post_e, "51peg", rng, ["51 Peg b"])[0],
                      "convergence": convergence(post_e)},
        "delta_log_like": ecc.log_like - circ.log_like,
        "delta_bic": ecc.bic - circ.bic,
        "lucy_sweeney": ls,
        "adopted": "eccentric" if ls["significant"] else "circular",
    }
    adopted = ecc if ls["significant"] else circ
    post = post_e if ls["significant"] else post_c
    res["instruments"] = instrument_table(data, adopted)
    art["fit"] = adopted
    art["post"] = post
    art["fit_ecc"] = ecc

    # --- 4. The true mass, from the planet's own orbital velocity --------------
    P, K, e, w, tc = post.planet_samples()
    mstar = mass_draws("51peg", rng, len(P))
    true_mass = {}
    for ref, (kp, up, dn) in KP_51PEG.items():
        kp_s = rng.normal(kp, (up + dn) / 2, len(P)) * 1e3  # km/s -> m/s
        q = K[:, 0] / kp_s  # mass ratio m_p / M_star
        m_true = q * mstar * phys.GM_SUN / phys.GM_JUP
        v_orb = phys.orbital_speed(P[:, 0], mstar * phys.GM_SUN * (1 + q))
        sin_i = np.clip((K[:, 0] + kp_s) / v_orb, None, 1.0)
        incl = np.degrees(np.arcsin(sin_i))
        true_mass[ref] = {"kp_kms": kp, "mass_mjup": fit.summarize(m_true),
                          "sin_i": fit.summarize(sin_i), "incl_deg": fit.summarize(incl),
                          "incl_deg_95_lower": float(np.percentile(incl, 5))}
    res["true_mass"] = true_mass

    # --- 5. Diagnostics -----------------------------------------------------------
    res["residual_search"] = res["search"][-1]
    res["lick_time_convention"] = _lick_time_check(data, adopted)
    res["keck_talor_check"] = _talor_check(adopted)
    res["literature"] = literature_rows("51 Peg b")
    if not quiet:
        pl = res["fits"][res["fits"]["adopted"]]["planet"]
        print(f"  P = {pl['P']['median']:.6f} +- {pl['P']['plus']:.6f} d, K = {pl['K']['median']:.2f} "
              f"+- {pl['K']['plus']:.2f} m/s, m sin i = {pl['msini_mjup']['median']:.3f} M_J; "
              f"e: Lucy-Sweeney p = {ls['p_circular']:.2f}")
    return res, art


def _lick_time_check(data: fit.RVData, res: fit.FitResult) -> dict:
    """Are Fischer et al.'s Lick "JD"s already barycentric? Refit with them
    taken as they are and compare likelihoods: the wrong convention shifts
    each time by up to +-8 minutes over the year, smearing the orbit's phase."""
    alt, _ = rvdata.build("51peg", lick_convert_time=False)
    elo_k = alt.inst_names.index("ELODIE")
    alt.activity["elodie_berv"] = np.where(alt.inst == elo_k, alt.activity["BERV"], 0.0)
    conv = data.t[np.isin(data.inst, [k for k, n in enumerate(data.inst_names) if n.startswith("Lick")])]
    raw = alt.t[np.isin(alt.inst, [k for k, n in enumerate(alt.inst_names) if n.startswith("Lick")])]
    model = fit.KeplerModel(alt, res.model.planets, activity=res.model.activity,
                            night_term=res.model.night_insts)
    alt_fit = fit.fit_map(model, _start_from(res, model))
    return {"log_like_converted": res.log_like, "log_like_as_is": alt_fit.log_like,
            "delta_log_like": res.log_like - alt_fit.log_like,
            "max_shift_minutes": float(np.max(np.abs(np.sort(conv) - np.sort(raw))) * 1440)}


def _talor_check(res: fit.FitResult) -> dict:
    """Keck velocities from Butler et al. 2017 vs Tal-Or et al.'s 2019
    zero-point-corrected version: same orbit, less jitter?"""
    out = {}
    for version in ("butler", "talor"):
        d, _ = rvdata.build("51peg", sources=("keck",), keck_version=version)
        pl = res.planets()[0]
        spec = fit.PlanetSpec(P=pl["P"], K=pl["K"], tc=pl["tc"], circular=True,
                              P_bounds=(pl["P"] * 0.999, pl["P"] * 1.001))
        f = fit.fit_map(fit.KeplerModel(d, [spec]))
        p = f.params()
        out[version] = {"n": d.n, "K": f.planets()[0]["K"],
                        "jitter": {k: v for k, v in p.items() if k.startswith(("jit", "night"))},
                        "rms": f.rms()}
    return out


# ---------------------------------------------------------------------------

def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def save_results(results: dict) -> None:
    RESULTS_PATH.write_text(json.dumps(_jsonable(results), indent=2))


# ---------------------------------------------------------------------------
# 55 Cancri: five planets, an alias and a transit
# ---------------------------------------------------------------------------

#: Published periods, used only to *name* what the blind search finds.
KNOWN_55CNC = {"b": 14.65, "c": 44.4, "d": 5000.0, "e": 0.7365, "f": 260.0}
CYCLE_BOUNDS_D = (3000.0, 4600.0)  # prior range for the magnetic-cycle term
D_BOUNDS_D = (4600.0, 9000.0)  # ...and planet d's, so the two cannot swap


def _name_55cnc(P: float) -> str | None:
    for name, P0 in KNOWN_55CNC.items():
        if abs(np.log(P / P0)) < np.log(1.25 if name != "d" else 1.5):
            return name
    return None


def _noise_draw(model: fit.KeplerModel, theta, rng):
    """One realisation of the fitted noise model: white noise (quoted error
    and jitter) plus a shared offset for each night of night-term instruments."""
    d = model.data
    white = np.sqrt(d.err ** 2 + model.jitter(theta)[0, d.inst] ** 2)
    noise = rng.normal(0, white)
    if model._g_starts.size:
        c = model._night_variance_per_group(theta)[0]
        sizes = np.diff(np.r_[model._g_starts, model._g_order.size])
        noise[model._g_order] += np.repeat(rng.normal(0, np.sqrt(c)), sizes)
    return noise


def alias_test(full: fit.FitResult, j: int, periods: tuple[float, float], rng, n_sim: int = 500) -> dict:
    """Which of two alias periods is the real one? (after Dawson & Fabrycky 2010)

    The velocities with every other signal removed (the full fit's
    residuals plus signal j) are the data. Each candidate period is taken in
    turn as the truth: a sinusoid with the amplitude and phase the data have
    at that period is added to fresh noise drawn from the fitted noise
    model, sampled at the real times, and the periodogram heights at both
    candidates are recorded. The true period must reproduce the heights
    actually observed. Returns, for each hypothesis, how often simulated data
    show the observed ordering of the two peaks, and how often a height ratio
    at least as extreme as the observed one.
    """
    m = full.model
    d = m.data
    r = m.residuals(full.theta) + m.keplerian(full.theta, planets=[j])[0]
    err = m.effective_errors(full.theta)
    fit_no_e = full  # the noise model to draw from
    T = d.baseline
    windows = [np.linspace(1 / P - 2.5 / T, 1 / P + 2.5 / T, 161) for P in periods]
    grid = np.concatenate(windows)

    def heights(y):
        p = phys.gls(d.t, y, err, grid, d.inst)
        return np.array([p[:161].max(), p[161:].max()])

    obs = heights(r)
    out = {"periods": list(periods), "observed_power": obs.tolist(),
           "observed_ratio": float(obs[0] / obs[1]), "hypotheses": {}, "residual": r, "err": err}
    for P_true in periods:
        pw, amp, ph = phys.gls(d.t, r, err, windows[periods.index(P_true)], d.inst, return_fit=True)
        i = int(np.argmax(pw))
        f_true = windows[periods.index(P_true)][i]
        sims = np.empty((n_sim, 2))
        for k in range(n_sim):
            y = amp[i] * np.cos(2 * np.pi * f_true * d.t - ph[i]) + _noise_draw(m, fit_no_e.theta, rng)
            sims[k] = heights(y)
        ratio = sims[:, 0] / sims[:, 1]
        same_order = np.mean((ratio > 1) == (out["observed_ratio"] > 1))
        as_extreme = (np.mean(ratio >= out["observed_ratio"]) if out["observed_ratio"] > 1
                      else np.mean(ratio <= out["observed_ratio"]))
        out["hypotheses"][f"{P_true:.4f}"] = {
            "amplitude": float(amp[i]), "fraction_same_order": float(same_order),
            "fraction_as_extreme": float(as_extreme),
            "ratio_median": float(np.median(ratio)), "ratio_16_84": np.percentile(ratio, [16, 84]).tolist(),
            "sim_ratios": ratio,
        }
    return out


def analyse_55cnc(rng, *, quiet: bool = False):
    """Planets one at a time; the 0.74-day planet hiding behind its alias."""
    if not quiet:
        print("55 Cancri")
    data, log = rvdata.build("55cnc", sources=("lick", "keck"))
    res: dict = {"cuts": log.as_dict(), "counts": data.counts(), "n": data.n,
                 "baseline_days": data.baseline}
    art: dict = {"data": data}
    freqs = phys.frequency_grid(data.t, PERIOD_MIN_D, None, OVERSAMPLE)

    # --- 1. Blind search: strongest signal first, refit everything each time --
    steps, search_fit = fit.search_planets(data, freqs, n_max=5, fap_threshold=FAP_THRESHOLD, quiet=quiet)
    names = [_name_55cnc(s.peak_period) for s in steps]
    res["search"] = [{"planets_before": s.n_planets_before, "peak_period": s.peak_period,
                      "peak_power": s.peak_power, "fap": s.fap, "rms_before": s.rms_before,
                      "identified_as": nm, "top": list(zip(s.top_periods, s.top_powers))}
                     for s, nm in zip(steps, names)]
    art["search"] = steps
    art["freqs"] = freqs
    found = [nm for nm in names[:len(search_fit.planets())]]
    if sorted(n for n in found if n) != sorted(KNOWN_55CNC):
        raise RuntimeError(f"blind search found {found}, expected all of b-f")

    # --- 2. The alias: the step where e appeared ---------------------------------
    k_e = names.index("e")
    se = steps[k_e]
    over_1d = (1 / freqs) > 1.0
    j = int(np.argmax(np.where(over_1d, se.power, 0)))
    res["e_step"] = {
        "step": k_e, "top": list(zip(se.top_periods, se.top_powers)),
        "best_period_search_below_1d": se.peak_period,
        "best_period_if_search_stops_at_1d": float(1 / freqs[j]), "power_at_that": float(se.power[j]),
        "butler_2017_table2": 2.81738,
        "aliases_predicted": aliases_of(KNOWN_55CNC["e"]),
    }

    # --- 3. The fiducial model ------------------------------------------------------
    sp = {nm: q for nm, q in zip(found, search_fit.planets())}
    order = ["b", "c", "d", "e", "f"]
    specs = []
    for nm in order:
        q = sp[nm]
        bounds = D_BOUNDS_D if nm == "d" else alias_safe_bounds(q["P"])
        specs.append(fit.PlanetSpec(P=q["P"], K=q["K"], tc=q["tc"], circular=True, P_bounds=bounds, label=nm))
    five = fit.fit_map(fit.KeplerModel(data, specs))
    # the star's magnetic cycle: a sinusoid seeded at the S-index period
    S_keck = data.activity["S"]
    keck = np.isin(data.inst, [i for i, n in enumerate(data.inst_names) if n.startswith("HIRES")])
    fS = phys.frequency_grid(data.t[keck], 2.0, None, OVERSAMPLE)
    pS = phys.gls(data.t[keck], S_keck[keck], np.full(keck.sum(), 0.01), fS)
    topS = phys.top_peaks(fS, pS, 4)
    P_cycle_S = float(1 / fS[topS[0]])
    res["s_index"] = {"n": int(keck.sum()), "peak_period": P_cycle_S, "peak_power": float(pS[topS[0]]),
                      "top": [(float(1 / fS[i]), float(pS[i])) for i in topS]}
    art["s_index"] = {"t": data.t[keck], "S": S_keck[keck], "freqs": fS, "power": pS}
    # The cycle's period is anchored to the S-index (+-0.5 d); amplitude and
    # phase are free. Left free as well, the period slides along a ridge
    # with planet d's (tested below as an alternative) and the model is not
    # identifiable.
    pw, amp, ph = fit.residual_periodogram(five, np.array([1 / P_cycle_S]))
    cyc = fit.planet_from_sinusoid(P_cycle_S, max(amp[0], 3.0), ph[0], five.model.t_ref)
    cyc.P_bounds = (P_cycle_S - 0.5, P_cycle_S + 0.5)
    cyc.label = "cycle"
    base = [replace(p, P=q["P"], K=q["K"], tc=q["tc"]) for p, q in zip(specs, five.planets())]
    circ6 = fit.fit_map(fit.KeplerModel(data, base + [cyc], t_ref=five.model.t_ref))
    fid = fit_eccentric_from(circ6, which=[0, 1, 2, 4])  # e (index 3) and the cycle stay circular
    post = run_posterior(fid, "55cnc_fiducial", quiet, steps=30000, burn=10000)
    labels = ["55 Cnc b", "55 Cnc c", "55 Cnc d", "55 Cnc e", "55 Cnc f", "magnetic cycle"]
    planets = planet_posteriors(post, "55cnc", rng, labels)
    for pl in planets[:5]:
        if "e" in pl:
            P, K, e_, w_, tc_ = post.planet_samples()
            j = labels.index(pl["label"])
            pl["lucy_sweeney"] = lucy_sweeney(e_[:, j])
    # the same masses with the other published stellar mass
    from config import STELLAR_MASS_ALTERNATIVE_55CNC as ALT
    scale = (ALT[0] / STELLAR_MASS["55cnc"][0]) ** (2 / 3)
    res["fiducial"] = {"log_like": fid.log_like, "bic": fid.bic, "rms": fid.rms(),
                       "planets": planets, "convergence": convergence(post),
                       "noise": {k: v for k, v in fid.params().items() if k.startswith(("jit", "night"))},
                       "mass_scale_if_alternative_mstar": scale, "alternative_mstar": ALT}
    res["instruments"] = instrument_table(data, fid)
    art["fit"] = fid
    art["post"] = post
    art["five_circ"] = five
    if not quiet:
        for pl in planets:
            print(f"  {pl['label']:14s} P = {pl['P']['median']:.5f} d, K = {pl['K']['median']:.2f} "
                  f"+- {pl['K']['plus']:.2f} m/s, m sin i = {pl['msini_mearth']['median']:.1f} M_E")

    # --- 4. Planet d depends on how the long-period power is modelled ---------------
    alts = {}
    five_ecc = fit_eccentric_from(five, which=[0, 1, 2, 4])
    alts["five planets, no cycle term"] = five_ecc
    pw, amp, ph = fit.residual_periodogram(five_ecc, freqs)
    top = phys.top_peaks(freqs, pw, 3)
    P_extra = float(1 / freqs[top[0]])
    extra = fit.planet_from_sinusoid(P_extra, amp[top[0]], ph[top[0]], five.model.t_ref)
    extra.P_bounds = (1000.0, 3000.0)
    base_e = [replace(p, P=q["P"], K=q["K"], tc=q["tc"], h=np.sqrt(q["e"]) * np.cos(q["omega"]),
                      k=np.sqrt(q["e"]) * np.sin(q["omega"])) for p, q in zip(five_ecc.model.planets, five_ecc.planets())]
    m_x = fit.KeplerModel(data, base_e + [extra], t_ref=five.model.t_ref)
    alts[f"five planets + {P_extra:.0f}-d sinusoid"] = fit.fit_map(m_x, _start_from(five_ecc, m_x))
    free = replace(fid.model.planets[5], P=fid.planets()[5]["P"], K=fid.planets()[5]["K"],
                   tc=fid.planets()[5]["tc"], P_bounds=CYCLE_BOUNDS_D)
    base_f = [replace(p, P=q["P"], K=q["K"], tc=q["tc"], h=np.sqrt(q["e"]) * np.cos(q["omega"]),
                      k=np.sqrt(q["e"]) * np.sin(q["omega"])) for p, q in zip(fid.model.planets[:5], fid.planets()[:5])]
    m_f = fit.KeplerModel(data, base_f + [free], t_ref=fid.model.t_ref)
    alts["five planets + cycle term, period free (3,000-4,600 d)"] = fit.fit_map(m_f, _start_from(fid, m_f))
    res["alternatives"] = {}
    for label, f_ in alts.items():
        pd_ = f_.planets()[2]
        res["alternatives"][label] = {"log_like": f_.log_like, "bic": f_.bic, "rms": f_.rms(),
                                      "d": pd_, "all": f_.planets()}
    res["alternatives"]["fiducial: five planets + cycle term at the S-index period"] = {
        "log_like": fid.log_like, "bic": fid.bic, "rms": fid.rms(), "d": fid.planets()[2], "all": fid.planets()}
    res["leftover_after_five"] = {"period": P_extra, "power": float(pw[top[0]]),
                                  "fap": float(phys.fap_baluev(pw[top[0]], data.t,
                                                               five_ecc.model.effective_errors(five_ecc.theta),
                                                               freqs.max(), data.n_inst,
                                                               five_ecc.model.n_effective())),
                                  "by_instrument": fit.signal_by_instrument(five_ecc, P_extra)}
    art["alternatives"] = alts
    pw_f, _, _ = fit.residual_periodogram(fid, freqs)
    topf = phys.top_peaks(freqs, pw_f, 4)
    res["residual_after_fiducial"] = {
        "top": [(float(1 / freqs[i]), float(pw_f[i])) for i in topf],
        "fap": float(phys.fap_baluev(pw_f[topf[0]], data.t, fid.model.effective_errors(fid.theta),
                                     freqs.max(), data.n_inst, fid.model.n_effective()))}
    art["residual_fiducial"] = pw_f

    # --- 5. The alias test: everything but e removed, noise model as fitted --------
    p_alias = res["e_step"]["best_period_if_search_stops_at_1d"]
    at = alias_test(fid, 3, (KNOWN_55CNC["e"], p_alias), rng)
    art["alias_test"] = at
    res["alias_test"] = {k: v for k, v in at.items() if k not in ("hypotheses", "residual", "err")}
    res["alias_test"]["hypotheses"] = {k: {kk: vv for kk, vv in v.items() if kk != "sim_ratios"}
                                      for k, v in at["hypotheses"].items()}
    art["no_e_power"] = phys.gls(data.t, at["residual"], at["err"], freqs, data.inst)
    if not quiet:
        hy = res["alias_test"]["hypotheses"]
        print(f"  alias test: observed power ratio {at['observed_ratio']:.2f}; "
              + "; ".join(f"if {k} d were true, same ordering in {v['fraction_same_order']:.1%}"
                          for k, v in hy.items()))

    # --- 6. HARPS: two nights in 2012 that the fit never saw ------------------------
    harps, hlog = rvdata.build("55cnc", sources=("harps",))
    th = fid.theta
    all_pl = fid.model.keplerian(th, harps.t)[0]
    only_e = fid.model.keplerian(th, harps.t, planets=[3])[0]
    nights = np.floor(harps.t - 0.3)
    out_n = []
    for with_e in (True, False):
        pred = all_pl if with_e else all_pl - only_e
        r = harps.y - pred
        for u in np.unique(nights):
            r[nights == u] -= np.average(r[nights == u], weights=1 / harps.err[nights == u] ** 2)
        out_n.append(float(np.std(r)))
    res["harps_check"] = {"n": harps.n, "nights": int(np.unique(nights).size),
                          "rms_with_e": out_n[0], "rms_without_e": out_n[1],
                          "median_error": float(np.median(harps.err))}
    art["harps"] = {"data": harps, "all": all_pl, "e": only_e}

    # --- 7. TESS -------------------------------------------------------------------------
    res["tess"], art["tess"] = _tess_55cnc(post, fid, rng, quiet)

    # --- 8. Literature -----------------------------------------------------------------
    res["literature"] = {nm: literature_rows(f"55 Cnc {nm}") for nm in order}
    res["alias_history"] = [r for r in literature_rows("55 Cnc e")
                            if r["pl_orbper"] is not None and r["pl_orbper"] > 1.5]
    return res, art


def _tess_55cnc(post: fit.Posterior, fid: fit.FitResult, rng, quiet: bool):
    """55 Cnc e's transit in TESS: which period, how big, and when."""
    import transit as tr
    from config import STELLAR_RADIUS_55CNC, TEFF_LOGG_FEH_55CNC

    t, f, e, sec = tr.load_light_curves()
    first = tr.flatten(t, f)
    ok = np.isfinite(first)
    bls = tr.bls_search(t[ok], first[ok])
    # re-flatten with the transits masked, so they do not bias the trend
    dt = np.mod(t - bls["t0"] + 0.5 * bls["best_period"], bls["best_period"]) - 0.5 * bls["best_period"]
    flat = tr.flatten(t, f, mask=np.abs(dt) < 0.06)
    ok = np.isfinite(flat)
    t, flat, e = t[ok], flat[ok], e[ok]

    u1, u2 = tr.claret_quadratic(*TEFF_LOGG_FEH_55CNC[:3])
    mstar, _, _ = STELLAR_MASS["55cnc"]
    rstar, rstar_err, _ = STELLAR_RADIUS_55CNC
    rho = mstar / rstar ** 3
    tc_mid = bls["t0"] + bls["best_period"] * np.round((np.median(t) - bls["t0"]) / bls["best_period"])
    a_rs0 = tr.a_over_r_from_density(bls["best_period"], rho)
    tf = tr.fit_transit(t, flat, e, bls["best_period"], tc_mid, a_rs0, u1, u2)
    p = tf["params"]; pe = tf["errors"]
    # radius, inclination, true mass and density, with uncertainties by sampling
    n = 20000
    x = rng.multivariate_normal([p[k] for k in ("tc", "P", "rp_rs", "b", "a_rs", "f0")], tf["cov"], n)
    rs = rng.normal(rstar, rstar_err, n)
    rp_earth = x[:, 2] * rs * tr.R_SUN_M / tr.R_EARTH_M
    inc = np.degrees(np.arccos(np.clip(x[:, 3] / x[:, 4], 0, 1)))
    P_, K_, e_, w_, tc_ = post.planet_samples()
    je = 3
    mst = mass_draws("55cnc", rng, len(P_))
    msini = phys.msini_mearth(P_[:, je], K_[:, je], e_[:, je], mst)
    idx = rng.integers(0, len(msini), n)
    mtrue = msini[idx] / np.sin(np.radians(inc))
    rho_p = phys.RHO_EARTH * mtrue / rp_earth ** 3
    # stellar density from the transit shape alone, as a consistency check
    rho_star_transit = 3 * np.pi * x[:, 4] ** 3 / (6.674e-11 * (x[:, 1] * 86400) ** 2) / 1408.0
    # the RV orbit's prediction of the transit time, 6-9 years after the last velocity
    n_orb = np.round((p["tc"] - tc_[:, je]) / P_[:, je])
    tc_pred = tc_[:, je] + n_orb * P_[:, je]
    offset_min = (np.median(tc_pred) - p["tc"]) * 1440
    folds = {}
    for label, P in (("true", p["P"]), ("alias", 2.81738)):
        tc0 = p["tc"] if label == "true" else np.median(tc_pred)
        ph, mean, se, cnt = tr.fold_binned(t, flat, P, tc0, bin_minutes=10.0)
        folds[label] = {"P": P, "phase": ph, "flux": mean, "err": se}
    ph_grid = np.linspace(-0.12, 0.12, 400)
    model_curve = p["f0"] * tr.transit_model(p["tc"] + ph_grid, p["tc"], p["P"], p["rp_rs"], p["b"],
                                              p["a_rs"], u1, u2)
    res = {
        "n_points": int(t.size), "sectors": sorted(set(int(s) for s in sec)),
        "bls": {k: v for k, v in bls.items() if k not in ("periods", "power")},
        "limb_darkening": [u1, u2],
        "fit": {k: {"value": p[k], "err": pe[k]} for k in p},
        "depth_ppm": float((1 - tr.transit_model(np.array([p["tc"]]), p["tc"], p["P"], p["rp_rs"], p["b"],
                                                  p["a_rs"], u1, u2)[0]) * 1e6),
        "radius_earth": fit.summarize(rp_earth), "inclination_deg": fit.summarize(inc),
        "true_mass_earth": fit.summarize(mtrue), "density_gcc": fit.summarize(rho_p),
        "stellar_density_transit_solar": fit.summarize(rho_star_transit),
        "stellar_density_adopted_solar": rho,
        "rv_predicted_transit_offset_minutes": float(offset_min),
        "rv_predicted_transit_sigma_minutes": float(np.std(tc_pred) * 1440),
        "orbits_extrapolated": int(np.median(n_orb)),
    }
    art = {"t": t, "flat": flat, "bls": bls, "folds": folds, "model_phase": ph_grid,
           "model_flux": model_curve, "rp_earth": rp_earth, "mtrue": mtrue, "u": (u1, u2)}
    if not quiet:
        print(f"  TESS: BLS {bls['best_period']:.5f} d; Rp = {res['radius_earth']['median']:.2f} "
              f"+- {res['radius_earth']['plus']:.2f} R_E, M = {res['true_mass_earth']['median']:.2f} M_E, "
              f"density {res['density_gcc']['median']:.1f} g/cc; RV-predicted transit "
              f"{offset_min:+.0f} +- {res['rv_predicted_transit_sigma_minutes']:.0f} min")
    return res, art


# ---------------------------------------------------------------------------
# Activity: shared tools
# ---------------------------------------------------------------------------

ACTIVITY_LABELS = {
    "Halpha": "H-alpha index (chromosphere)",
    "NaD1": "Na D1 index (chromosphere)",
    "RHK": "Ca II H&K, R'HK (chromosphere)",
    "CRX": "chromatic index (colour-dependent velocity)",
    "dLW": "differential line width",
    "FWHM": "cross-correlation width (FWHM)",
    "BIS": "bisector span (line asymmetry)",
}


def indicator_periodograms(data: fit.RVData, names, *, period_min=2.0, planet_periods=()):
    """Periodogram of each activity indicator (zero point per instrument),
    with its strongest peaks and its power at the given planet periods."""
    out, art = {}, {}
    for col in names:
        v = data.activity[col]
        ok = np.isfinite(v)
        if ok.sum() < 30:
            continue
        sub = data.subset(ok)
        fr = phys.frequency_grid(sub.t, period_min, None, OVERSAMPLE)
        p = phys.gls(sub.t, v[ok], np.ones(ok.sum()), fr, sub.inst)
        top = phys.top_peaks(fr, p, 5)
        fap = float(phys.fap_baluev(p[top[0]], sub.t, np.ones(ok.sum()), fr.max(), sub.n_inst))
        at = {}
        for P in planet_periods:
            j = int(np.argmin(np.abs(fr - 1 / P)))
            sl = slice(max(j - 8, 0), j + 9)
            at[f"{P:.4g}"] = float(np.max(p[sl]))
        out[col] = {"n": int(ok.sum()), "top": [(float(1 / fr[i]), float(p[i])) for i in top],
                    "fap_top": fap, "power_at_planets": at}
        art[col] = {"freqs": fr, "power": p}
    return out, art


def stacked_power(t, r, err, P, rng, amp_vec=None, n_sim=300, step=5, n0=20):
    """How a signal at period P builds up as the data accumulate, in time
    order (after Mortier & Collier Cameron 2017, who do this with a Bayesian
    periodogram). The quantity tracked is the chi-square improvement from
    fitting a sinusoid at P, which for a strictly periodic signal grows in
    proportion to the number of measurements; the band is what the same
    signal, injected into white noise with the data's errors and sampling,
    produces. (The normalised periodogram power would not do: it is a
    fraction of the variance and levels off for any steady signal.) A planet
    stays in the band; activity, whose spots are born and die, need not."""
    f0 = np.array([1 / P])
    N = np.arange(n0, t.size + 1, step)
    w = 1 / err ** 2

    def dchi2(y, n):
        yn, wn = y[:n], w[:n]
        chi2_0 = np.sum(wn * (yn - np.average(yn, weights=wn)) ** 2)
        return phys.gls(t[:n], yn, err[:n], f0)[0] * chi2_0

    obs = np.array([dchi2(r, n) for n in N])
    if amp_vec is None:
        ph = 2 * np.pi * t / P
        X = np.column_stack([np.cos(ph), np.sin(ph), np.ones(t.size)])
        amp_vec = np.linalg.solve(X.T @ (w[:, None] * X), X.T @ (w * r))[:2]
    sims = np.empty((n_sim, N.size))
    for k in range(n_sim):
        y = amp_vec[0] * np.cos(2 * np.pi * t / P) + amp_vec[1] * np.sin(2 * np.pi * t / P) + rng.normal(0, err)
        sims[k] = [dchi2(y, n) for n in N]
    lo, med, hi = np.percentile(sims, [5, 50, 95], axis=0)
    return {"N": N, "dchi2": obs, "band_lo": lo, "band_med": med, "band_hi": hi,
            "fraction_outside": float(np.mean((obs < lo) | (obs > hi)))}


def epoch_consistency(t, r, err, P, edges_year) -> dict:
    """Amplitude and phase of a sinusoid at P fitted separately in each
    epoch, and a chi-square test that they agree (a planet's must)."""
    from astropy.time import Time
    from scipy import stats

    yr = Time(t, format="jd").decimalyear
    rows, vecs, covs = [], [], []
    for a, b in zip(edges_year[:-1], edges_year[1:]):
        m = (yr >= a) & (yr < b)
        ph = 2 * np.pi * t[m] / P
        X = np.column_stack([np.cos(ph), np.sin(ph), np.ones(m.sum())])
        W = 1 / err[m] ** 2
        cov = np.linalg.inv(X.T @ (W[:, None] * X))
        beta = cov @ (X.T @ (W * r[m]))
        rows.append({"from": a, "to": b, "n": int(m.sum()), "amplitude": float(np.hypot(*beta[:2])),
                     "amplitude_err": float(np.sqrt(np.mean(np.diag(cov[:2, :2])))),
                     "phase_deg": float(np.degrees(np.arctan2(beta[1], beta[0])))})
        vecs.append(beta[:2]); covs.append(cov[:2, :2])
    inv = [np.linalg.inv(c) for c in covs]
    mean = np.linalg.solve(sum(inv), sum(i @ v for i, v in zip(inv, vecs)))
    chi2 = float(sum((v - mean) @ i @ (v - mean) for v, i in zip(vecs, inv)))
    dof = 2 * (len(vecs) - 1)
    return {"epochs": rows, "chi2": chi2, "dof": dof, "p_consistent": float(stats.chi2.sf(chi2, dof))}


def _sinusoid_term(fit_res: fit.FitResult, P: float, bounds, label: str) -> fit.PlanetSpec:
    """A circular Keplerian (i.e. a sinusoid) seeded from the residuals at P,
    used to model a stellar-activity signal alongside the planets."""
    pw, amp, ph = fit.residual_periodogram(fit_res, np.array([1 / P]))
    spec = fit.planet_from_sinusoid(P, max(amp[0], 0.3), ph[0], fit_res.model.t_ref)
    spec.P_bounds = bounds
    spec.label = label
    return spec


# ---------------------------------------------------------------------------
# Proxima Centauri b: hard mode
# ---------------------------------------------------------------------------

PROXIMA_INDICATORS = ("Halpha", "NaD1", "RHK", "CRX", "dLW", "FWHM", "BIS")
PROXIMA_ROT_BOUNDS = (70.0, 110.0)


def analyse_proxima(rng, *, quiet: bool = False):
    """A 1.2 m/s planet around a flaring red dwarf."""
    from config import PROXIMA_ROTATION_D
    if not quiet:
        print("Proxima Centauri")
    data, log = rvdata.build("proxima", sources=("harps",), flare_cut=True)
    res: dict = {"cuts": log.as_dict(), "counts": data.counts(), "n": data.n,
                 "baseline_days": data.baseline}
    art: dict = {"data": data}
    freqs = phys.frequency_grid(data.t, PERIOD_MIN_D, None, OVERSAMPLE)

    # 1. blind search: the first peak, then what is left
    steps, one = fit.search_planets(data, freqs, n_max=1, fap_threshold=FAP_THRESHOLD, quiet=quiet)
    res["search"] = [{"planets_before": s.n_planets_before, "peak_period": s.peak_period,
                      "peak_power": s.peak_power, "fap": s.fap,
                      "top": list(zip(s.top_periods, s.top_powers))} for s in steps]
    art["search"] = steps
    art["freqs"] = freqs
    P_b = steps[0].peak_period
    if abs(P_b - 11.19) > 0.05:
        raise RuntimeError(f"Proxima b not the strongest signal (found {P_b:.3f} d)")
    nxt = steps[1]
    f_sid = 1 / phys.SIDEREAL_DAY_D
    res["next_signal"] = {"period": nxt.peak_period, "fap": nxt.fap,
                          "daily_alias_of": float(1 / abs(1 / nxt.peak_period - f_sid))}

    # 2. activity indicators: rotation, and nothing at the planet
    ind, ind_art = indicator_periodograms(data, PROXIMA_INDICATORS, planet_periods=(P_b,))
    res["indicators"] = ind
    art["indicators"] = ind_art
    res["rotation_literature_d"] = PROXIMA_ROTATION_D

    # 3. fiducial: planet b + a sinusoid for rotation-modulated activity
    pw, _, _ = fit.residual_periodogram(one, freqs)
    band = (1 / freqs >= PROXIMA_ROT_BOUNDS[0]) & (1 / freqs <= PROXIMA_ROT_BOUNDS[1])
    P_rot_rv = float(1 / freqs[band][np.argmax(pw[band])])
    rot = _sinusoid_term(one, P_rot_rv, PROXIMA_ROT_BOUNDS, "rotation")
    b_spec = replace(one.model.planets[0], P=one.planets()[0]["P"], K=one.planets()[0]["K"],
                     tc=one.planets()[0]["tc"], P_bounds=alias_safe_bounds(P_b), label="b")
    m_fid = fit.KeplerModel(data, [b_spec, rot], t_ref=one.model.t_ref)
    fid = fit.fit_map(m_fid, _start_from(one, m_fid))
    post = run_posterior(fid, "proxima_fiducial", quiet)
    pl = planet_posteriors(post, "proxima", rng, ["Proxima Cen b", "rotation signal"])
    res["fiducial"] = {"log_like": fid.log_like, "bic": fid.bic, "rms": fid.rms(),
                       "planet": pl[0], "rotation_term": {k: pl[1][k] for k in ("P", "K")},
                       "convergence": convergence(post),
                       "noise": {k: v for k, v in fid.params().items() if k.startswith(("jit", "night"))}}
    art["fit"] = fid
    art["post"] = post
    # K_b under other treatments of activity (a systematic on the mass)
    alt = {"planet only": one, "planet + rotation (fiducial)": fid}
    q = fid.planets()
    half = _sinusoid_term(fid, q[1]["P"] / 2, (35.0, 55.0), "half rotation")
    specs = [replace(s, P=x["P"], K=x["K"], tc=x["tc"]) for s, x in zip(fid.model.planets, q)]
    m_h = fit.KeplerModel(data, specs + [half], t_ref=fid.model.t_ref)
    alt["planet + rotation + half-rotation"] = fit.fit_map(m_h, _start_from(fid, m_h))
    for acts in (("dLW",), ("dLW", "Halpha", "CRX")):
        m_a = fit.KeplerModel(data, [b_spec], activity=acts, t_ref=one.model.t_ref)
        alt["planet + linear " + " + ".join(acts)] = fit.fit_map(m_a, _start_from(one, m_a))
    mstar = STELLAR_MASS["proxima"][0]
    res["activity_models"] = {k: {"K_b": f.planets()[0]["K"], "P_b": f.planets()[0]["P"],
                                  "msini_mearth": float(phys.msini_mearth(f.planets()[0]["P"], f.planets()[0]["K"],
                                                                          0.0, mstar)),
                                  "log_like": f.log_like, "bic": f.bic} for k, f in alt.items()}
    # 4. does the planet's signal build up like a planet's should?
    r = fid.model.residuals(fid.theta) + fid.model.keplerian(fid.theta, planets=[0])[0]
    err = fid.model.effective_errors(fid.theta)
    sp = stacked_power(data.t, r - np.average(r, weights=1 / err ** 2), err, q[0]["P"], rng)
    res["stacked_power_b"] = {"fraction_outside_band": sp["fraction_outside"], "final_dchi2": float(sp["dchi2"][-1])}
    art["stacked_b"] = sp
    res["literature"] = literature_rows("Proxima Cen b")
    if not quiet:
        p0 = pl[0]
        print(f"  Proxima b: P = {p0['P']['median']:.4f} d, K = {p0['K']['median']:.2f} +- {p0['K']['plus']:.2f} "
              f"m/s, m sin i = {p0['msini_mearth']['median']:.2f} M_E; rotation term "
              f"{pl[1]['P']['median']:.1f} d")
    return res, art


# ---------------------------------------------------------------------------
# Gliese 581 d: a planet that was stellar activity
# ---------------------------------------------------------------------------

def analyse_gj581(rng, *, quiet: bool = False):
    from config import GJ581_ROTATION_D, GJ581D_PUBLISHED
    from astropy.time import Time
    if not quiet:
        print("Gliese 581")
    data, log = rvdata.build("gj581", sources=("harps",))
    res: dict = {"cuts": log.as_dict(), "counts": data.counts(), "n": data.n,
                 "baseline_days": data.baseline}
    art: dict = {"data": data}
    freqs = phys.frequency_grid(data.t, PERIOD_MIN_D, None, OVERSAMPLE)

    # 1. blind search, exactly as for any other star
    steps, four = fit.search_planets(data, freqs, n_max=4, fap_threshold=FAP_THRESHOLD, quiet=quiet)
    res["search"] = [{"planets_before": s.n_planets_before, "peak_period": s.peak_period,
                      "peak_power": s.peak_power, "fap": s.fap,
                      "top": list(zip(s.top_periods, s.top_powers))} for s in steps]
    art["search"] = steps
    art["freqs"] = freqs
    periods = [p["P"] for p in four.planets()]
    i_d = int(np.argmin(np.abs(np.array(periods) - 66.8)))
    if abs(periods[i_d] - 66.8) > 2:
        raise RuntimeError(f"no ~66.8-d signal among {periods}")
    res["signal_d"] = {"P": periods[i_d], "K": four.planets()[i_d]["K"], "published": GJ581D_PUBLISHED}

    # 2. the star's rotation in H-alpha, and its half-period harmonic
    ind, ind_art = indicator_periodograms(data, ("Halpha", "NaD1", "NaD2", "CRX", "dLW"),
                                          planet_periods=(periods[i_d], 5.3687, 12.914, 3.149))
    res["indicators"] = ind
    art["indicators"] = ind_art
    ha = data.activity["Halpha"]
    ha_err = data.activity["Halpha_err"]
    fr = ind_art["Halpha"]["freqs"]
    p1, a1, ph1 = phys.gls(data.t, ha, ha_err, fr, return_fit=True)
    i1 = int(np.argmax(p1))
    ha_res = ha - (a1[i1] * np.cos(2 * np.pi * fr[i1] * data.t - ph1[i1]))
    p2 = phys.gls(data.t, ha_res, ha_err, fr)
    top2 = phys.top_peaks(fr, p2, 5)
    near = (1 / fr > 60) & (1 / fr < 76)
    j = int(np.argmax(np.where(near, p2, 0)))
    res["halpha_rotation"] = {"period": float(1 / fr[i1]), "power": float(p1[i1]),
                              "literature_rotation": GJ581_ROTATION_D,
                              "after_removing_rotation_top": [(float(1 / fr[k]), float(p2[k])) for k in top2],
                              "half_rotation_peak": {"period": float(1 / fr[j]), "power": float(p2[j]),
                                                     "fap": float(phys.fap_baluev(p2[j], data.t, ha_err, fr.max()))}}
    art["halpha"] = {"freqs": fr, "power": p1, "power_after": p2}

    # 3. velocities against H-alpha, with the three accepted planets removed
    keep = [k for k in range(len(periods)) if k != i_d]
    specs3 = [replace(four.model.planets[k], P=periods[k], K=four.planets()[k]["K"], tc=four.planets()[k]["tc"],
                      P_bounds=alias_safe_bounds(periods[k])) for k in keep]
    m3 = fit.KeplerModel(data, specs3, t_ref=four.model.t_ref)
    three = fit.fit_map(m3, _start_from(four, m3))
    r3 = three.model.residuals(three.theta)
    err3 = three.model.effective_errors(three.theta)
    yr = Time(data.t, format="jd").decimalyear
    season = np.floor(yr - 0.3).astype(int)
    corr = {"all": float(np.corrcoef(r3, ha)[0, 1])}
    for s in np.unique(season):
        m = season == s
        if m.sum() >= 10:
            corr[str(s)] = float(np.corrcoef(r3[m], ha[m])[0, 1])
    res["rv_halpha_correlation"] = corr
    art["three"] = three
    art["season"] = season

    # 4. models: does linear H-alpha decorrelation remove the 66.8-d signal?
    data.activity["Ha"] = ha
    for s in np.unique(season):
        sel = season == s
        data.activity[f"Ha_{s}"] = np.where(sel, ha - ha[sel].mean(), 0.0)
    per_season = tuple(f"Ha_{s}" for s in np.unique(season) if (season == s).sum() >= 10)
    models = {"b, c, e": three, "b, c, e, d": four}
    for label, acts in (("b, c, e + H-alpha", ("Ha",)), ("b, c, e + H-alpha by season", per_season)):
        m = fit.KeplerModel(data, specs3, activity=acts, t_ref=four.model.t_ref)
        models[label] = fit.fit_map(m, _start_from(three, m))
    res["models"] = {}
    for label, f in models.items():
        pw, _, _ = fit.residual_periodogram(f, freqs)
        k = int(np.argmin(np.abs(freqs - 1 / periods[i_d])))
        sl = slice(k - 30, k + 31)
        kk = sl.start + int(np.argmax(pw[sl]))
        entry = {"log_like": f.log_like, "bic": f.bic, "rms": f.rms(),
                 "power_near_66.8d": float(pw[kk]),
                 "fap_near_66.8d": float(phys.fap_baluev(pw[kk], data.t, f.model.effective_errors(f.theta),
                                                         freqs.max(), 1, f.model.n_effective()))}
        if "H-alpha" in label and "season" not in label:
            lin, le = f.model.linear_params(f.theta)
            entry["halpha_coefficient"] = [lin["c[Ha]"], le["c[Ha]"]]
        res["models"][label] = entry

    # 5. is the 66.8-d signal coherent, as a planet's must be?
    res["epochs"] = epoch_consistency(data.t, r3, err3, periods[i_d], [2004, 2008.9, 2013])
    sp = stacked_power(data.t, r3, err3, periods[i_d], rng)
    res["stacked_power_d"] = {"fraction_outside_band": sp["fraction_outside"]}
    art["stacked_d"] = sp

    # 6. the accepted system, with the 66.8-d activity signal modelled as a sinusoid
    act = _sinusoid_term(three, periods[i_d], (60.0, 75.0), "activity (half rotation)")
    base = [replace(s, P=x["P"], K=x["K"], tc=x["tc"]) for s, x in zip(three.model.planets, three.planets())]
    m_f = fit.KeplerModel(data, base + [act], t_ref=three.model.t_ref)
    circ = fit.fit_map(m_f, _start_from(three, m_f))
    # b and c get eccentric orbits; e (K ~ 1.7 m/s) cannot constrain one and
    # stays circular, as tides would also make it at 0.03 AU.
    i_e = int(np.argmin([abs(s.P - 3.149) for s in circ.model.planets[:3]]))
    fid = fit_eccentric_from(circ, which=[k for k in range(3) if k != i_e])
    post = run_posterior(fid, "gj581_fiducial", quiet)
    order = np.argsort([s.P for s in fid.model.planets[:3]])
    names = {5.37: "GJ 581 b", 12.9: "GJ 581 c", 3.15: "GJ 581 e"}
    labels = []
    for s in fid.model.planets[:3]:
        labels.append(min(names.items(), key=lambda kv: abs(kv[0] - s.P))[1])
    labels.append("activity (half rotation)")
    pls = planet_posteriors(post, "gj581", rng, labels)
    P_, K_, e_, w_, tc_ = post.planet_samples()
    for j, pl in enumerate(pls[:3]):
        if "e" in pl:
            pl["lucy_sweeney"] = lucy_sweeney(e_[:, j])
    res["fiducial"] = {"log_like": fid.log_like, "bic": fid.bic, "rms": fid.rms(), "planets": pls,
                       "convergence": convergence(post)}
    art["fit"] = fid
    art["post"] = post
    res["literature"] = {lab: literature_rows(lab) for lab in labels[:3]}
    if not quiet:
        for pl in pls:
            print(f"  {pl['label']:26s} P = {pl['P']['median']:.4f} d, K = {pl['K']['median']:.2f} "
                  f"+- {pl['K']['plus']:.2f} m/s, m sin i = {pl['msini_mearth']['median']:.2f} M_E")
    return res, art
