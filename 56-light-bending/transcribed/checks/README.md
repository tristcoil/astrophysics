# Checks of the transcribed tables

These scripts check the CSV files one folder up against the papers' own arithmetic: sums, means,
normal equations, solutions, residuals and weighted means. Each prints a report, and the notes
files quote those reports. They read only the CSVs, so they run from anywhere:

    python transcribed/checks/principe.py

| script | what it checks | notes file | the notes call it |
|---|---|---|---|
| `sobral4inch_tables1to6.py` | Sobral 4-inch: Table I's predicted displacements (one Sun, one 1/r law), the p.304 coefficients, Table II, the normal equations and solutions of Tables III–VI | `dyson1920_sobral4inch_notes.md` | `data/scratch/T1/verify_t1.py` |
| `sobral4inch_results.py` | Sobral 4-inch: Table VII, p.307, Table VIII, p.309, p.331, the exposures | `dyson1920_sobral4inch_results_notes.md` | `data/scratch/T2/verify_t2.py` |
| `astrographic.py`, `astrographic_lsq.py` | Sobral astrographic: Table IX, the p.311 Dy table and α + 0.24e list, the least-squares means | `dyson1920_astrographic_notes.md` | `data/scratch/T3-audit/scripts/checks.py`, `check_lsq.py` |
| `harvey1979.py` | Harvey 1979, Tables I–III | `harvey1979_summary.md` | `data/scratch/T3-audit/scripts/check_harvey.py` |
| `principe.py` | Príncipe: Tables XI–XIII, the plate reductions, p.328 | `dyson1920_principe_notes.md` | `data/scratch/T4/tools/verify.py` |
| `wallal.py`, `wallal_tables2_3.py` | Wallal: the measures, weights, means and solutions of the 1928 and 1923 tables; the 1928 plate constants (Tables 2–3) | `wallal1928_notes.md` | `data/scratch/T5/scripts/verify_all.py`, `verify_t23.py` |

The scratch folders named in the notes are git-ignored. The scripts here are copies of them, with
the path to the CSVs made relative, and they print the same reports. The second, blind readings
of every table are not kept in the repository. The notes give their counts and results.
