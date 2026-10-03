"""Download and cache every Horizons query the project needs.

    python fetch_data.py            # fetch whatever is missing
    python fetch_data.py --list     # show the plan without fetching
    python fetch_data.py --refresh snapshot_   # re-fetch keys starting "snapshot_"

JPL Horizons is free and needs no account. It is also a shared service that
JPL's own navigators rely on, so requests go one at a time with a pause
between them, and nothing already on disk is fetched again. The full plan is
83 requests and about 75 MB of text, and takes about four minutes the first time.
After that the analysis runs offline.
"""

from __future__ import annotations

import argparse
import sys
import time

import horizons as hz
from config import HORIZONS_DIR, MANIFEST_PATH

#: Seconds between requests. Horizons has no published rate limit, but it is
#: shared infrastructure and there is no hurry.
PAUSE_S = 0.5


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--list", action="store_true", help="print the query plan and exit")
    parser.add_argument("--refresh", action="append", default=[], metavar="PREFIX",
                        help="re-download queries whose key starts with PREFIX (repeatable)")
    args = parser.parse_args(argv)

    plan = hz.query_plan()
    if args.list:
        for q in plan:
            cached = "cached " if q.path.exists() else "missing"
            span = f"{q.start} .. {q.stop} step {q.step}"
            print(f"  [{cached}] {q.key:<32} {q.kind:<9} {span}")
        print(f"\n{len(plan)} queries")
        return 0

    n_new = n_cached = 0
    failures = []
    t0 = time.time()
    for i, q in enumerate(plan, 1):
        refresh = any(q.key.startswith(p) for p in args.refresh)
        try:
            text, downloaded = hz.fetch(q, refresh=refresh)
        except hz.HorizonsError as exc:
            failures.append(str(exc))
            print(f"  [{i:>3}/{len(plan)}] FAILED {q.key}")
            continue
        if downloaded:
            n_new += 1
            rows = text.count("\n", text.index("$$SOE"), text.index("$$EOE")) - 1
            print(f"  [{i:>3}/{len(plan)}] {q.key:<32} {rows:>7,} rows   ({q.purpose})")
            time.sleep(PAUSE_S)
        else:
            n_cached += 1

    size_mb = sum(p.stat().st_size for p in HORIZONS_DIR.glob("*.txt")) / 1e6
    print(f"\n{n_new} downloaded, {n_cached} already cached, {len(failures)} failed "
          f"({time.time() - t0:.0f} s). Cache: {size_mb:.1f} MB in {HORIZONS_DIR.name}/")
    print(f"Provenance of every query: {MANIFEST_PATH.relative_to(MANIFEST_PATH.parents[1])}")
    for f in failures:
        print("\n" + f, file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
