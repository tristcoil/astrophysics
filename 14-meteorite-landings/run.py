"""Run the whole project end to end.

    python run.py

Downloads the data if it is not already cached, cleans and labels the
catalogue, measures everything, writes the figures and the interactive map, and
prints every number the README quotes. Once the cache is warm this runs offline.
"""

from __future__ import annotations

import argparse
import sys
import warnings

warnings.filterwarnings("ignore")

import analysis
from config import FIGURE_DIR
from fetch_data import fetch_all, load_manifest


def _pct(p: dict, digits: int = 1) -> str:
    lo, hi = p["ci95"]
    return (f"{p['fraction']:.{digits}%} ({p['k']:,}/{p['n']:,}; 95% CI "
            f"{lo:.{digits}%}-{hi:.{digits}%})")


def print_summary(r: dict) -> None:
    m = load_manifest()["meteorites"]
    c = r["census"]
    q = r["data_quality"]
    print("\n" + "=" * 78)
    print(f"Meteoritical Society catalogue, 2013 snapshot (NASA; retrieved {m['retrieved_utc'][:10]})")
    print("=" * 78)
    print(f"  {c['rows']:,} rows: {c['valid']:,} meteorites + {c['relict']} relicts; "
          f"{c['falls']:,} falls, {c['finds']:,} finds; {c['n_classes']} class labels")
    e = c["by_environment"]
    print(f"  finds by environment: Antarctic ice {e['antarctic']:,}, hot deserts {e['desert']:,}, "
          f"individual finds {e['chance']:,}, unknown {e['unknown']}, on Mars {e['mars']}")

    print("\nDATA REPAIRS AND FLAGS")
    print(f"  no usable location: {q['no_location']:,} ({q['coords_missing']:,} blank, "
          f"{q['coords_placeholder_0_0']:,} at 0,0, {q['coords_zeroed_latitude']} Yamato with "
          f"latitude zeroed, {q['coords_mars']} on Mars)")
    print(f"  Northwest Africa: {q['northwest_africa']['n']:,} meteorites, "
          f"{q['northwest_africa']['located']} with a find site")
    print(f"  years before 1583 shifted by one: {q['year_calendar_fixed']} "
          f"(e.g. {', '.join(f'{x['name']} {x['year_raw']:.0f}->{x['year']:.0f}' for x in q['calendar_fixed_examples'][:3])})")
    print(f"  years after the snapshot: {[x['name'] for x in q['year_after_snapshot']]}")
    nr = r["naming_rule"]
    print(f"  season-year numbering rule: {_pct(nr['of_which_antarctic'], 2)} of located "
          f"meteorites following it are Antarctic; applied to "
          f"{sum(nr['applied_to'].values()):,} without coordinates")
    o = r["ocean"]
    print(f"  offshore after snapping to coasts: {len(o['offshore'])} "
          f"({', '.join(x['name'] for x in o['offshore'])})")

    print("\nWHERE FINDS COME FROM")
    ca = r["collection_areas"]
    print(f"  Antarctica: {_pct(ca['antarctic_share_of_located_finds'])} of located finds, "
          f"from {ca['antarctica_share_of_land']:.1%} of the land "
          f"({ca['antarctica_area_km2'] / 1e6:.2f} million km^2, "
          f"{ca['antarctica_share_of_earth']:.1%} of Earth)")
    print(f"  finds per million km^2: Antarctica {ca['density_antarctica']:,.0f}, "
          f"rest of land {ca['density_rest_of_land']:,.0f}")
    for d in ca["densest_countries"][:6]:
        print(f"    {d['country']:<26}{d['finds']:>7,} finds  {d['per_million_km2']:>9,.0f} per M km^2")
    print(f"  {ca['n_areas']} collection areas; the top 10 hold {ca['top10_share_of_finds']:.0%} "
          f"of {ca['n_finds']:,} finds")

    print("\nWITNESSED FALLS FOLLOW PEOPLE")
    w = r["witness"]
    f = w["fit"]
    print(f"  {w['falls_in_window']} falls >= {w['min_mass_g']:.0f} g in {w['years'][0]}-{w['years'][1]}"
          f" = {w['falls_per_year']:.2f}/yr, vs {w['expected_per_year_on_land']:,.0f}/yr landing on "
          f"land (ex-Antarctica): 1 in {1 / w['global_efficiency']:,.0f} recorded")
    print(f"  recovery ~ density^gamma: gamma = {f['gamma']:.2f} +- {f['gamma_se_overdispersed']:.2f} "
          f"(bootstrap 95% {f['gamma_ci95_bootstrap'][0]:.2f}-{f['gamma_ci95_bootstrap'][1]:.2f}); "
          f"at 100 people/km^2: {f['eps0']:.3%}; LR test vs gamma=0: {f['lr_stat_gamma0']:.0f}; "
          f"dispersion {f['pearson_dispersion']:.2f}; {f['n_countries']} countries")
    print(f"  one km^2 receives a >= 100 g fall once every {w['years_per_km2_per_fall']:,.0f} years")
    print("  recovery by mass (recorded vs published flux on land):")
    for ref in r["recovery_vs_mass"]["references"]:
        print(f"    >= {ref['mass_g']:>5} g: {ref['recorded_per_year']:.2f}/yr vs "
              f"{ref['expected_per_year_on_land']:>7,.0f}/yr -> 1 in {1 / ref['efficiency']:,.0f}"
              f"   [{ref['source'].split(',')[0]}]")

    print("\nWHAT FALLS VS WHAT IS FOUND (irons + stony-irons)")
    for label, row in r["evatt"].items():
        pub = row["published"]
        pubs = f"   Evatt+2016: {pub['fraction']:.2%} ({pub['k']}/{pub['n']:,})" if pub else ""
        print(f"  {label:<30}{_pct(row['ours'], 2)}{pubs}")
    comp = r["composition"]
    for s in analysis.SAMPLES:
        fams = comp[s]["families"]
        print(f"  {analysis.SAMPLE_NAMES[s]:<38} n={comp[s]['n']:>6,}  " +
              "  ".join(f"{k[:5]} {v['fraction']:5.1%}" for k, v in fams.items()))

    print("\nMASSES")
    ms = r["masses"]
    for env, v in ms["median_g"].items():
        print(f"  median {analysis.ENV_NAMES[env]:<30}{v['median']:>10,.1f} g "
              f"(95% CI {v['ci95'][0]:,.1f}-{v['ci95'][1]:,.1f})")
    print(f"  irons: {ms['mass_share_by_family']['iron']:.0%} of {ms['total_mass_kg'] / 1e3:,.0f} t "
          f"total; ten heaviest: {ms['top10_share_of_mass']:.0%}; heaviest "
          f"{ms['heaviest'][0]['name']} {ms['heaviest'][0]['mass_g'] / 1e6:g} t; heaviest stone "
          f"{ms['heaviest_stones'][0]['name']} {ms['heaviest_stones'][0]['mass_g'] / 1e6:g} t")

    print("\nHISTORY")
    h = r["history"]
    print(f"  falls per year: {h['falls_per_year_1800_1899']:.1f} in the 1800s, "
          f"{h['falls_per_year_1900_2009']:.1f} in 1900-2009; oldest: "
          f"{h['oldest_falls'][0]['name']} {h['oldest_falls'][0]['year']:.0f}")
    print(f"  {h['share_of_finds_since_1969']:.0%} of finds made since 1969; busiest year "
          f"{h['peak_find_year']['year']} ({h['peak_find_year']['n']:,} finds); first Antarctic "
          f"find {h['first_antarctic_find']['name']} {h['first_antarctic_find']['year']}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--no-figures", action="store_true", help="skip figure generation")
    args = parser.parse_args(argv)

    fetch_all(quiet=True)
    df, wtable, countries, results = analysis.analyse()
    analysis.save_results(results)

    if not args.no_figures:
        import figures
        import interactive

        written = figures.make_all(df, wtable, countries, results)
        written.append(interactive.write_map(df, results))
        print(f"wrote {len(written)} figure files to {FIGURE_DIR}")

    print_summary(results)
    return 0


if __name__ == "__main__":
    sys.exit(main())
