"""The analysis: every measurement the figures and the README rely on.

Everything here works from the cached files and never touches the network. It
produces the cleaned, labelled catalogue (the figures draw from it) and a
results dictionary holding every number quoted in the README, saved to
``data/results.json`` so the prose can be checked against the code.

The organising idea: meteorites *fall* almost uniformly over the Earth, so any
pattern in the catalogue that is not uniform tells us about how meteorites are
*found*. Each function below measures one such pattern.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import binomtest, chi2

import catalog
import geo
from config import (
    ANTARCTIC_LATITUDE,
    DESERT_TRADE_NAMES,
    EVATT_2016,
    FLUX_REFERENCES,
    POPULATION_JSON,
    POPULATION_YEARS,
    RESULTS_PATH,
    WITNESS_MIN_MASS_G,
    WITNESS_REFERENCE_RATE,
)

ENVIRONMENTS = ("fall", "antarctic", "desert", "chance")
ENV_NAMES = {
    "fall": "Witnessed falls",
    "antarctic": "Antarctic ice",
    "desert": "Hot-desert collection areas",
    "chance": "Individual finds elsewhere",
}

RNG_SEED = 20130514  # the date of the NASA snapshot's last update

#: Population-density bins (people per km^2) for pooling countries in the witness test.
DENSITY_BIN_EDGES = (0.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1e5)


# --- Small statistics helpers ----------------------------------------------------------

def proportion(k: int, n: int) -> dict:
    """A fraction with its 95% Wilson score interval."""
    if n == 0:
        return {"k": 0, "n": 0, "fraction": None, "ci95": [None, None]}
    ci = binomtest(int(k), int(n)).proportion_ci(confidence_level=0.95, method="wilson")
    return {"k": int(k), "n": int(n), "fraction": k / n, "ci95": [ci.low, ci.high]}


def poisson_ci(k: int, conf: float = 0.95) -> tuple[float, float]:
    """Exact (Garwood) confidence interval on a Poisson mean, given k events."""
    a = 1 - conf
    lo = 0.0 if k == 0 else chi2.ppf(a / 2, 2 * k) / 2
    hi = chi2.ppf(1 - a / 2, 2 * k + 2) / 2
    return float(lo), float(hi)


def bootstrap_median(x, *, n_boot: int = 2000, seed: int = RNG_SEED) -> dict:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    rng = np.random.default_rng(seed)
    meds = np.median(rng.choice(x, size=(n_boot, len(x)), replace=True), axis=1)
    return {"median": float(np.median(x)), "ci95": [float(np.percentile(meds, 2.5)),
                                                   float(np.percentile(meds, 97.5))],
            "n": int(len(x))}


# --- Census and data quality ------------------------------------------------------------

def census(df: pd.DataFrame, ledger: dict) -> dict:
    v = df[df["valid"]]
    by_env = v["environment"].value_counts()
    return {
        "rows": ledger["rows"],
        "valid": int(len(v)),
        "relict": ledger["relict"],
        "falls": int((v["environment"] == "fall").sum()),
        "finds": int(v["fall"].eq("Found").sum()),
        "by_environment": {e: int(by_env.get(e, 0)) for e in (*ENVIRONMENTS, "unknown", "mars")},
        "located": int(v["located"].sum()),
        "located_by_environment": {e: int((v["located"] & v["environment"].eq(e)).sum())
                                   for e in ENVIRONMENTS},
        "year_range": [int(v["year"].min()), int(v["year"].max())],
        "n_classes": int(df["recclass"].nunique()),
    }


def data_quality(df: pd.DataFrame, ledger: dict) -> dict:
    early = df[df["year_calendar_fixed"]].sort_values("year")
    q = dict(ledger)
    q["calendar_fixed_examples"] = early[["name", "year_raw", "year"]].to_dict("records")
    q["no_location"] = int((~df["located"] & df["valid"]).sum())
    q["no_location_by_area"] = (df[~df["located"] & df["valid"] & df["area"].notna()]
                                ["area"].value_counts().head(8).to_dict())
    nwa = df["area"].eq("Northwest Africa")
    q["northwest_africa"] = {"n": int(nwa.sum()), "located": int((nwa & df["located"]).sum())}
    return q


def naming_rule_check(df: pd.DataFrame) -> dict:
    """How reliable is 'season-year numbering means Antarctica'? Test it on the
    meteorites whose coordinates we know."""
    known = df["located"] & df["valid"] & df["area"].notna() & ~df["fell"]
    rule = known & df["antarctic_number"]
    ant = df["lat"] < ANTARCTIC_LATITUDE
    areas = df[rule].groupby("area")["lat"].median()
    return {
        "meteorites_following_rule": int(rule.sum()),
        "of_which_antarctic": proportion(int((rule & ant).sum()), int(rule.sum())),
        "areas_following_rule": int(len(areas)),
        "areas_not_antarctic": sorted(areas[areas >= ANTARCTIC_LATITUDE].index),
        "desert_meteorites_following_rule": int((known & ~ant & df["antarctic_number"]).sum()),
        "applied_to": (df[df["environment"].eq("antarctic") & df["area_lat"].isna()]
                       ["area"].value_counts().to_dict()),
    }


# --- Where the finds come from ----------------------------------------------------------

def collection_areas(df: pd.DataFrame, countries: list) -> dict:
    v = df[df["valid"] & ~df["fell"] & df["environment"].isin(["antarctic", "desert"])]
    g = (v.groupby("area")
         .agg(n=("id", "size"), located=("located", "sum"),
              environment=("environment", "first"),
              lat=("lat", lambda s: s[v.loc[s.index, "located"]].median()),
              lon=("lon", lambda s: s[v.loc[s.index, "located"]].median()),
              first_year=("year", "min"), last_year=("year", "max"))
         .sort_values("n", ascending=False))
    finds = df[df["valid"] & ~df["fell"] & df["environment"].isin(["antarctic", "desert", "chance"])]
    top = g.head(20).reset_index()
    top["share_of_finds"] = top["n"] / len(finds)

    located_finds = finds[finds["located"]]
    ant_area = next(c.area_km2 for c in countries if c.name == "Antarctica")
    land_area = sum(c.area_km2 for c in countries)
    n_ant = int((located_finds["lat"] < ANTARCTIC_LATITUDE).sum())

    # Finds per million km^2, country by country.
    per_country = (located_finds[located_finds["country_idx"] >= 0]
                   .groupby("country_idx").size())
    rows = []
    for i, n in per_country.items():
        c = countries[i]
        rows.append({"country": c.name, "finds": int(n), "area_km2": c.area_km2,
                     "per_million_km2": n / c.area_km2 * 1e6})
    dens = pd.DataFrame(rows).sort_values("per_million_km2", ascending=False)
    with_finds = len(dens)
    zero = [c.name for i, c in enumerate(countries)
            if i not in per_country.index and c.area_km2 > 3e5]
    return {
        "n_areas": int(len(g)),
        "top": top.to_dict("records"),
        "top10_share_of_finds": float(g["n"].head(10).sum() / len(finds)),
        "n_finds": int(len(finds)),
        "located_finds": int(len(located_finds)),
        "antarctic_share_of_located_finds": proportion(n_ant, len(located_finds)),
        "antarctica_area_km2": ant_area,
        "antarctica_share_of_earth": ant_area / geo.EARTH_AREA_KM2,
        "antarctica_share_of_land": ant_area / land_area,
        "land_area_km2": land_area,
        "density_antarctica": n_ant / ant_area * 1e6,
        "density_rest_of_land": (len(located_finds) - n_ant) / (land_area - ant_area) * 1e6,
        "densest_countries": dens.head(12).to_dict("records"),
        "countries_with_finds": with_finds,
        "large_countries_without_finds": zero,
        "all_country_density": dens.to_dict("records"),
    }


def ocean(df: pd.DataFrame, countries: list) -> dict:
    """The 71% of the Earth nobody searches."""
    land = sum(c.area_km2 for c in countries)
    off = df[df["valid"] & df["located"] & df["country_idx"].lt(0)]
    return {
        "land_share_of_earth": land / geo.EARTH_AREA_KM2,
        "offshore": off[["name", "fall", "environment", "year", "lat", "lon", "recclass",
                         "mass_g"]].to_dict("records"),
        "snapped_to_coast": int((df["valid"] & df["located"] & df["snap_km"].gt(0)).sum()),
    }


# --- What falls vs what is found ---------------------------------------------------------

#: The samples compared in the composition figure. Desert finds are split in two:
#: meteorites bought through the trade (named for a whole region, find site
#: unknown) and meteorites from desert collection areas whose find site was
#: recorded -- by search teams or by collectors with a GPS.
SAMPLES = ("fall", "antarctic", "desert_searched", "desert_trade", "chance")
SAMPLE_NAMES = {
    "fall": "Witnessed falls",
    "antarctic": "Antarctic ice",
    "desert_searched": "Hot deserts, find site recorded",
    "desert_trade": "Hot deserts, bought (site unknown)",
    "chance": "Individual finds elsewhere",
}


def sample_masks(df: pd.DataFrame) -> dict:
    trade = df["area"].isin(DESERT_TRADE_NAMES)
    desert = df["environment"].eq("desert")
    return {
        "fall": df["environment"].eq("fall"),
        "antarctic": df["environment"].eq("antarctic"),
        "desert_searched": desert & ~trade,
        "desert_trade": desert & trade,
        "chance": df["environment"].eq("chance"),
    }


def composition(df: pd.DataFrame) -> dict:
    """Family shares in each sample (classified meteorites only)."""
    out = {}
    v = df["valid"] & df["family"].ne("unclassified")
    for key, m in sample_masks(df).items():
        s = df[v & m]
        fam = s["family"].value_counts()
        h, l_, ll = (int(s["group"].eq(g).sum()) for g in ("H", "L", "LL"))
        out[key] = {
            "n": int(len(s)),
            "families": {f: proportion(int(fam.get(f, 0)), len(s)) for f in catalog.FAMILIES},
            "iron_or_stony_iron": proportion(int(s["family"].isin(["iron", "stony_iron"]).sum()),
                                             len(s)),
            "martian": int(s["group"].eq("Martian").sum()),
            "lunar": int(s["group"].eq("Lunar").sum()),
            "h_fraction_of_h_plus_l": proportion(h, h + l_),
            "ordinary_groups": {"H": h, "L": l_, "LL": ll},
        }
    return out


def evatt_comparison(df: pd.DataFrame) -> dict:
    """Reproduce Evatt et al. (2016) Table 1 from the 2013 snapshot.

    Their counts include every named meteorite (classified or not); so do ours.
    """
    v = df[df["valid"]]
    iron = v["family"].isin(["iron", "stony_iron"])
    ant = v["environment"].eq("antarctic")
    finds_rest = v["fall"].eq("Found") & ~ant & v["environment"].ne("mars")
    rows = {
        "LaPaz Icefield": v["area"].eq("LaPaz Icefield"),
        "Frontier Mountain": v["area"].eq("Frontier Mountain"),
        "All Antarctic finds": ant,
        "Falls (rest of world)": v["fell"],
        "Finds (rest of world)": finds_rest,
        # Our split of the last row:
        "  of which hot-desert areas": v["environment"].eq("desert"),
        "  of which individual finds": v["environment"].eq("chance"),
    }
    out = {}
    for label, m in rows.items():
        ours = proportion(int((m & iron).sum()), int(m.sum()))
        pub = EVATT_2016.get(label)
        out[label] = {"ours": ours,
                      "published": proportion(pub["k"], pub["n"]) if pub else None}
    return out


# --- Masses ------------------------------------------------------------------------------

def masses(df: pd.DataFrame) -> dict:
    v = df[df["valid"] & df["environment"].isin(ENVIRONMENTS)]
    out = {"median_g": {}, "quartiles_g": {}}
    for env in ENVIRONMENTS:
        m = v.loc[v["environment"] == env, "mass_g"].dropna()
        out["median_g"][env] = bootstrap_median(m)
        out["quartiles_g"][env] = [float(m.quantile(0.25)), float(m.quantile(0.75))]
    heavy = v.dropna(subset=["mass_g"]).sort_values("mass_g", ascending=False)
    out["heaviest"] = heavy.head(12)[["name", "recclass", "family", "mass_g", "fall", "year",
                                      "country"]].to_dict("records")
    stones = heavy[~heavy["family"].isin(["iron", "stony_iron"])]
    out["heaviest_stones"] = stones.head(3)[["name", "recclass", "mass_g", "fall", "year",
                                             "country"]].to_dict("records")
    total = v["mass_g"].sum()
    by_fam = v.groupby("family")["mass_g"].sum()
    out["total_mass_kg"] = float(total / 1e3)
    out["mass_share_by_family"] = {f: float(by_fam.get(f, 0) / total) for f in catalog.FAMILIES}
    out["count_share_by_family"] = {f: float((v["family"] == f).mean()) for f in catalog.FAMILIES}
    top10 = heavy.head(10)["mass_g"].sum()
    out["top10_share_of_mass"] = float(top10 / total)
    return out


# --- The witness test ---------------------------------------------------------------------

def load_population(countries: list) -> dict:
    """Mean World Bank population over POPULATION_YEARS, keyed by ISO3 code."""
    meta, rows = json.loads(POPULATION_JSON.read_text())
    pop = pd.DataFrame([(r["countryiso3code"], int(r["date"]), r["value"]) for r in rows],
                       columns=["iso3", "year", "pop"]).dropna()
    y0, y1 = POPULATION_YEARS
    pop = pop[(pop["year"] >= y0) & (pop["year"] <= y1)]
    full = pop.groupby("iso3")["year"].nunique() == (y1 - y0 + 1)
    mean = pop.groupby("iso3")["pop"].mean()
    wanted = {c.iso3 for c in countries}
    return {k: float(v) for k, v in mean.items() if k in wanted and full.get(k, False)}


def world_population() -> dict:
    """World Bank world totals (code WLD) for the first and last year of the window."""
    meta, rows = json.loads(POPULATION_JSON.read_text())
    w = {int(r["date"]): r["value"] for r in rows if r["countryiso3code"] == "WLD"}
    y0, y1 = POPULATION_YEARS
    return {"first_year": y0, "last_year": y1, "first": float(w[y0]), "last": float(w[y1]),
            "ratio": float(w[y1] / w[y0])}


def fit_efficiency(n, expected, density, *, d_ref: float = 100.0) -> dict:
    """Poisson fit of N_i ~ Poisson(E_i * eps0 * (D_i / d_ref)^gamma).

    E_i is the number of falls country i would have recorded if every fall were
    found; eps0 is the fraction recovered at population density d_ref per km^2;
    gamma says how recovery scales with density (0 = not at all, 1 = in direct
    proportion to the number of people).
    """
    n = np.asarray(n, dtype=float)
    logE = np.log(np.asarray(expected, dtype=float))
    x = np.log(np.asarray(density, dtype=float) / d_ref)

    def nll(p):
        mu = np.exp(logE + p[0] + p[1] * x)
        return float(np.sum(mu - n * (logE + p[0] + p[1] * x)))

    def grad(p):
        mu = np.exp(logE + p[0] + p[1] * x)
        r = mu - n
        return np.array([r.sum(), (r * x).sum()])

    res = minimize(nll, x0=[np.log(max(n.sum(), 1) / np.exp(logE).sum()), 0.5], jac=grad,
                   method="BFGS")
    p = res.x
    mu = np.exp(logE + p[0] + p[1] * x)
    # Fisher information for the log-link Poisson model.
    info = np.array([[mu.sum(), (mu * x).sum()], [(mu * x).sum(), (mu * x * x).sum()]])
    cov = np.linalg.inv(info)
    dof = len(n) - 2
    pearson = float(np.sum((n - mu) ** 2 / mu) / dof)
    # Likelihood-ratio test against gamma = 0 (recovery independent of people).
    p0 = np.log(n.sum() / np.exp(logE).sum())
    nll0 = nll([p0, 0.0])
    return {
        "log_eps0": float(p[0]), "gamma": float(p[1]),
        "eps0": float(np.exp(p[0])),
        "gamma_se_poisson": float(np.sqrt(cov[1, 1])),
        "gamma_se_overdispersed": float(np.sqrt(cov[1, 1] * max(pearson, 1.0))),
        "pearson_dispersion": pearson,
        "lr_stat_gamma0": float(2 * (nll0 - res.fun)),
        "n_countries": int(len(n)), "n_falls": int(n.sum()), "d_ref": d_ref,
    }


def witness(df: pd.DataFrame, countries: list) -> tuple[dict, pd.DataFrame]:
    """Do witnessed falls follow people? Country by country, 1960-2009."""
    y0, y1 = POPULATION_YEARS
    years = y1 - y0 + 1
    pops = load_population(countries)
    falls = df[df["valid"] & df["fell"] & df["year"].between(y0, y1)
               & df["mass_g"].ge(WITNESS_MIN_MASS_G)]
    flux = WITNESS_REFERENCE_RATE / 1e6  # per km^2 per year

    counts = falls[falls["country_idx"] >= 0].groupby("country_idx").size()
    rows = []
    for i, c in enumerate(countries):
        if c.name == "Antarctica":
            continue
        rows.append({
            "country": c.name, "iso3": c.iso3, "continent": c.continent,
            "area_km2": c.area_km2, "population": pops.get(c.iso3, np.nan),
            "falls": int(counts.get(i, 0)),
        })
    t = pd.DataFrame(rows)
    t["expected_if_all_found"] = flux * t["area_km2"] * years
    t["density"] = t["population"] / t["area_km2"]
    t["efficiency"] = t["falls"] / t["expected_if_all_found"]
    ci = np.array([poisson_ci(k) for k in t["falls"]])
    t["eff_lo"] = ci[:, 0] / t["expected_if_all_found"]
    t["eff_hi"] = ci[:, 1] / t["expected_if_all_found"]

    usable = t["population"].notna() & (t["area_km2"] > 0)
    fit = fit_efficiency(t.loc[usable, "falls"], t.loc[usable, "expected_if_all_found"],
                         t.loc[usable, "density"])

    # Bootstrap over countries: robust to countries differing for reasons the
    # model does not know about (overdispersion).
    rng = np.random.default_rng(RNG_SEED)
    u = t[usable].reset_index(drop=True)
    gammas, eps = [], []
    for _ in range(1000):
        s = u.iloc[rng.integers(0, len(u), len(u))]
        f = fit_efficiency(s["falls"], s["expected_if_all_found"], s["density"])
        gammas.append(f["gamma"])
        eps.append(f["eps0"])
    fit["gamma_ci95_bootstrap"] = [float(np.percentile(gammas, 2.5)),
                                   float(np.percentile(gammas, 97.5))]
    fit["eps0_ci95_bootstrap"] = [float(np.percentile(eps, 2.5)), float(np.percentile(eps, 97.5))]

    # Countries pooled in bins of population density: fewer, sturdier points.
    binned = []
    for lo, hi in zip(DENSITY_BIN_EDGES[:-1], DENSITY_BIN_EDGES[1:]):
        b = u[(u["density"] >= lo) & (u["density"] < hi)]
        k, e = int(b["falls"].sum()), float(b["expected_if_all_found"].sum())
        clo, chi = poisson_ci(k)
        w_d = np.exp(np.average(np.log(b["density"]), weights=b["area_km2"]))
        binned.append({"density_lo": lo, "density_hi": hi, "density_typical": float(w_d),
                       "countries": int(len(b)), "falls": k, "expected": e,
                       "efficiency": k / e, "ci95": [clo / e, chi / e]})
    fit["binned"] = binned

    # Global: every recorded fall of >= 100 g, located or not, against all land
    # outside Antarctica (where no fall has ever been witnessed).
    land = sum(c.area_km2 for c in countries if c.name != "Antarctica")
    n_all = int(len(falls))
    expected_all = flux * land * years
    no_pop_falls = t[~usable & (t["falls"] > 0)][["country", "falls"]].to_dict("records")
    enough = t[usable & (t["expected_if_all_found"] >= 20)].sort_values("efficiency",
                                                                        ascending=False)
    result = {
        "years": [y0, y1], "min_mass_g": WITNESS_MIN_MASS_G,
        "reference_rate_per_million_km2_yr": WITNESS_REFERENCE_RATE,
        "falls_in_window": n_all,
        "falls_in_window_unlocated": int((~falls["located"]).sum()),
        "falls_per_year": n_all / years,
        "land_area_ex_antarctica_km2": land,
        "expected_per_year_on_land": flux * land,
        "global_efficiency": n_all / expected_all,
        "years_per_km2_per_fall": 1.0 / flux,
        "fit": fit,
        "falls_without_population_data": no_pop_falls,
        "best_countries": enough.head(8)[["country", "falls", "expected_if_all_found",
                                          "efficiency", "density"]].to_dict("records"),
        "zero_fall_countries_big": t[usable & (t["falls"] == 0)
                                     & (t["expected_if_all_found"] >= 100)]
        .sort_values("expected_if_all_found", ascending=False)["country"].head(10).tolist(),
    }
    return result, t


def recovery_vs_mass(df: pd.DataFrame, countries: list) -> dict:
    """Recorded falls per year above a mass, against the published flux on land."""
    y0, y1 = POPULATION_YEARS
    years = y1 - y0 + 1
    falls = df[df["valid"] & df["fell"] & df["year"].between(y0, y1)]
    m = np.sort(falls["mass_g"].dropna().to_numpy())
    land = sum(c.area_km2 for c in countries if c.name != "Antarctica")
    refs = []
    for r in FLUX_REFERENCES:
        recorded = int((m >= r["mass_g"]).sum())
        expected = r["rate"] * land / 1e6
        refs.append({**r, "recorded_per_year": recorded / years,
                     "expected_per_year_on_land": expected,
                     "efficiency": recorded / years / expected})
    return {
        "years": [y0, y1],
        "falls_with_mass": int(len(m)),
        "falls_without_mass": int(falls["mass_g"].isna().sum()),
        "masses_g": m.tolist(),
        "land_area_ex_antarctica_km2": land,
        "references": refs,
    }


# --- History ---------------------------------------------------------------------------------

def history(df: pd.DataFrame) -> dict:
    v = df[df["valid"] & df["environment"].isin(ENVIRONMENTS) & df["year"].notna()]
    per_year = (v.groupby(["year", "environment"]).size().unstack(fill_value=0)
                .reindex(columns=list(ENVIRONMENTS), fill_value=0))
    per_year.index = per_year.index.astype(int)
    falls = v[v["environment"] == "fall"]
    decades = (falls["year"] // 10 * 10).astype(int).value_counts().sort_index()
    finds = v[v["environment"] != "fall"]
    ant = v[v["environment"] == "antarctic"]
    first_ant = ant.sort_values("year").iloc[0]
    since_1969 = (finds["year"] >= 1969).mean()
    peak = per_year.drop(columns="fall").sum(axis=1)
    return {
        "per_year": {int(y): {e: int(r[e]) for e in ENVIRONMENTS} for y, r in per_year.iterrows()},
        "falls_per_decade": {int(k): int(n) for k, n in decades.items()},
        "falls_per_year_1900_2009": float(((falls["year"] >= 1900) & (falls["year"] <= 2009)).sum()
                                          / 110),
        "falls_per_year_1800_1899": float(((falls["year"] >= 1800) & (falls["year"] <= 1899)).sum()
                                          / 100),
        "undated": int((df["valid"] & df["environment"].isin(ENVIRONMENTS)
                        & df["year"].isna()).sum()),
        "first_antarctic_find": {"name": first_ant["name"], "year": int(first_ant["year"])},
        "share_of_finds_since_1969": float(since_1969),
        "peak_find_year": {"year": int(peak.idxmax()), "n": int(peak.max())},
        "world_population": world_population(),
        "oldest_falls": falls.sort_values("year").head(5)[["name", "year", "year_raw", "recclass",
                                                          "country"]].to_dict("records"),
    }


# --- Everything -------------------------------------------------------------------------------

def analyse():
    """Load, clean and measure. Returns (catalogue, country witness table, results)."""
    df, ledger, countries = catalog.load_catalog()
    wit, wtable = witness(df, countries)
    results = {
        "census": census(df, ledger),
        "data_quality": data_quality(df, ledger),
        "naming_rule": naming_rule_check(df),
        "collection_areas": collection_areas(df, countries),
        "ocean": ocean(df, countries),
        "composition": composition(df),
        "evatt": evatt_comparison(df),
        "masses": masses(df),
        "witness": wit,
        "recovery_vs_mass": recovery_vs_mass(df, countries),
        "history": history(df),
    }
    return df, wtable, countries, results


def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        return None if not np.isfinite(obj) else float(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def save_results(results: dict, path=RESULTS_PATH) -> None:
    slim = json.loads(json.dumps(_jsonable(results)))
    slim["recovery_vs_mass"].pop("masses_g", None)
    slim["collection_areas"].pop("all_country_density", None)
    path.write_text(json.dumps(slim, indent=2, ensure_ascii=False) + "\n")
