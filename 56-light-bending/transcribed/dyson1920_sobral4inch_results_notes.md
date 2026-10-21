# Sobral 4-inch results, exposures and conditions in Dyson, Eddington & Davidson (1920): transcription notes

These notes cover the T2 tables. They are Table VII and the text results (p.306), the two tables of p.307,
Table VIII (p.308), the observed and calculated displacements (p.309), the radial displacements (p.331), the
exposure tables (p.300) and the observing conditions (pp.295-302, 309-312, 330-331).

Source: Dyson, F. W., Eddington, A. S. & Davidson, C. 1920, *Phil. Trans. R. Soc. A* 220, 291-333. Scans:
archive.org item `philtrans06337895`, leaf N = printed page 291 + N. Every table digit was read on the
full-resolution 5120 × 6600 JPEG 2000 leaves (`06337895_00NN.jp2`): leaves 9, 10, 15, 16, 17, 18 and 40, with
leaf 2 for the Table I star positions used in check 8. Each block was read once from 0.5× crops and then re-read in
full a second time. Any cell that the arithmetic or OCR questioned was viewed again at 2-3× zoom. None turned out to
be a misreading.

## Files

| File | Printed source | Content |
|---|---|---|
| `dyson1920_table7_alpha.csv` | Table VII, p.306 | alpha per plate (r): 7 rows × 4 columns plus the Mean row; the text results of pp.306-307 (and their restatement on pp.330-331) quoted in the header |
| `dyson1920_p307_refraction_aberration.csv` | p.307, first table | a, e, b, d (r) for the eclipse plates, the scale plate and the 7 comparison plates |
| `dyson1920_p307_scale_orientation.csv` | p.307, second table | scale from x and y, orientation from x and y, adopted scale and orientation (r), 14 plates |
| `dyson1920_table8_residuals.csv` | Table VIII, p.308 | 196 plate values plus 28 printed means (arcsec), long format |
| `dyson1920_p309_obs_calc.csv` | p.309 | observed and calculated R.A. and Decl. displacements of the 7 stars (arcsec) |
| `dyson1920_p331_radial.csv` | p.331 | "Radial Displacement of Individual Stars": calculation and observation (arcsec) |
| `dyson1920_sobral_exposures.csv` | p.300 | 19 astrographic and 8 4-inch eclipse exposures, as printed, plus clearly marked derived columns |
| `dyson1920_sobral_conditions.md` | pp.295-302, 306, 309-312, 330-331 | station, instruments, weather, timing, Dr. Morize's temperatures, the July comparison exposures (p.301 table), plate selection, the astrographic focus change and the decision; with quotations |

Conventions:

- Values are the printed strings: explicit sign and the printed number of decimals. Zeros printed without a sign
  stay unsigned: "0.00" twice in Table VIII, and "0.000" for the adopted scale of the seven eclipse plates on p.307.
- Old typesetting is written plainly (0r·105 becomes 0.105). Plain ASCII is used in all CSV cells.
- Plate labels follow T1's convention: eclipse I II III IV V VII VIII; comparison 14_2a 14_2b 15_1 15_2 17_1 17_2
  18_2 (printed as subscripts); the scale plate is 18_1.
- Star numbers are those of Table I, in the printed row order: 11 5 4 3 6 10 2 (pp.308-309) and 11 10 6 5 4 2 3
  (p.331).
- Axes as in Table I: x = R.A. (east +), y = Decl. (north +), unit 50′. 1 r = 6″.25 (p.302).
- **Added, not printed**: the plate labels in the Table VII CSV (Table VII prints none; see check 2); the
  `station_instrument` column of the p.331 CSV (the table has one group only; "Sobral 4-inch" comes from the sentence
  introducing it); and the `*_derived`, `plate_type`, `gmt_astronomical` (fields carried down) and
  `plate_label_tables2to8` columns of the exposures CSV.

## Checks

The script is `data/scratch/T2/verify_t2.py` (git-ignored); its output is in `data/scratch/T2/verify_t2_output.txt`.
It reads only the final CSVs, plus T1's Table I, p.304, Table II and Tables III-VI CSVs (read-only) for
cross-checks. The Horizons and altitude computations are in `data/scratch/T2/work/` (`contacts_horizons.py`,
`altitude_check2.py`).

Result codes:

- **PASS**: within the stated tolerance.
- **ORIGINAL**: the deviation is in the 1920 printing or arithmetic. The value has been re-read, is kept as printed,
  and is explained under Discrepancies.
- **INFO**: a sanity check with no tolerance.

| # | Check | n | max dev | Result |
|---|---|---|---|---|
| 0 | Format: sign and printed decimals in every numeric cell; row counts 8/9/14/224/7/7/27 | 418 | - | PASS |
| 1 | Table VII Mean row = mean of the 7 values in each column (+0.11957, +0.01457, +0.12886, +0.03114) | 4 | 0.0004 r | PASS |
| 2 | Table VII = printed alpha solutions of Tables III, IV, V, VI (T1 CSV), plates in order I..VIII and 14_2a..18_2 | 28 | 0 | PASS: identical; this fixes the plate labels of Table VII |
| 3 | p.306 arithmetic: 0.120 - 0.015 = 0.105; 0.129 - 0.031 = 0.098; (0.105 + 2 × 0.098)/3 = 0.1003 → 0.100; 0.100 × 6.25 = 0.625″; 0.625″ × 50/15.8 = 1.978″ → 1″.98 | 5 | 0.002″ | PASS |
| 3b | The same results restated in §39 (p.330): 0r.098 → 1″.938 (printed 1″.94); 0r.105 → 2″.077 (printed 2″.06) | 2 | 0.02″ | ORIGINAL (R5) |
| 3c | Probable error from the scatter of Table VII (7 + 7 plates, eclipse and comparison scatter combined, weights 1:2): R.A. 8.1 %, Decl. 5.9 %, mean 4.8 %. Paper: "about 6 per cent" (p.306); ±0″.12 is 6.1 % of 1″.98 (p.331) | 28 | - | INFO: same size; the paper's method is not stated |
| 4 | p.307 scale/orientation from the solutions of Tables III-VI plus the refraction and aberration difference. Scale from x = a + Δa, scale from y = e + Δe, orientation from x = b + Δb, orientation from y = -(d + Δd), where Δ = row of the first p.307 table minus the scale-plate row | 56 | 0.010 r | ORIGINAL: 55/56 exact; 14_2b scale from x (R3) |
| 4b | Adopted orientation = mean of orientation from x and from y | 14 | 0.0005 r | PASS |
| 4c | Adopted scale = mean of the printed scale values of a group: July 14 plates +0.0130 (adopted +0.013); July 15-18 plates -0.0844 (adopted -0.084); eclipse plates -0.0026 (adopted 0.000) | 3 groups | - | PASS for the comparison groups; eclipse 0.000 is not explained (R8) |
| 4d | First p.307 table: b = d on all 9 rows; 14_2a = 14_2b; a, e, b never increase as the field altitude printed on p.301 rises (scale plate 33.6° to 18_2 at 39.2°) | 9 | - | PASS |
| 5 | Table VIII printed Mean = mean of the 7 plate values (4 blocks × 7 stars) | 28 | 0.0057″ | ORIGINAL: 27/28 within 0.005″; eclipse Decl. star 11 (R2) |
| 5b | Table VIII column sums: R.A. -0.03 to +0.02; Decl. eclipse +3.49 to +3.53, comparison +1.29 to +1.34. A misread digit would push one column out of line | 32 | - | PASS |
| 6 | p.309 observed = Table VIII eclipse Mean - comparison Mean (this is how the page defines it, "Subtracting the results of the comparison plates") | 14 | 0.20″ | ORIGINAL: 13/14 exact; star 6 R.A. sign (R1) |
| 7 | p.309 calculated = Table I Sobral displacements (T1 CSV) | 14 | 0.01″ | ORIGINAL: 13/14 identical; star 4 Decl. +0.74 vs +0.73 (R7) |
| 7b | p.309 calculated = K × (Ex, Ey) of p.304 with one K | 14 | - | PASS: all 14 for K = 0.5501-0.5513″, i.e. 1″.75 at a solar radius of 15′.72-15′.75. K = 0.5530 (radius 15′.8) gives 12/14 |
| 8 | p.331 = p.309 components projected on the unit vector from the Sun to the star. Star x, y from Table I (leaf 2). The Sun is fitted to the 25 regular Table I Sobral displacements: x = +0.4364, y = -0.2627, K = 0.5504″ | 14 | calc 0.005″, obs 0.017″ | ORIGINAL: calculation 7/7 within rounding, observation 6/7; star 6 needs R.A. -0.10 (R1). Same result with the Sun at (0.44, -0.26) |
| 8b | Alternative for the Calculation column: K/r instead of projecting the rounded p.309 values | 7 | 0.008″ | Calculation is the projection: K/r misses stars 10 and 4 by 0.006-0.008″, and star 3 sits on the rounding edge |
| 9 | Table VIII rebuilt from Table II (T1 CSV) with the p.307 adopted scale, orientation and refraction; zero point per column = mean | 196 | 0.44″ | ORIGINAL: 165/196 within 0.0105″, 174/196 within 0.02″. Every deviation above 0.02″ is in column 17_1 or in a column containing one of the six cells of R4 (their wrong value also shifts the column zero point) |
| 9b | Same with 3 free constants per column (zero point, x and y terms) | 196 | 0.22″ | ORIGINAL: 22/28 columns fit all 7 stars within 0.0105″. Each of the other 6 has one outlying star; without it the other six fit within 0.0082″ (R4) |
| 10 | Exposures (p.300): 19 astrographic exposures alternating 5 and 10 s, 8 4-inch exposures of 28 s (as p.299 says), no overlaps (gaps of 4-11 s). The first exposures start 3 s (4-inch) and 5 s (astrographic) after the assumed second contact; the last end 302 s and 305 s after it | 27 | - | PASS |
| 11 | Totality at the station (JPL Horizons, topocentric, 5 s steps; 40.35° W, latitude as printed): 11:58:10-12:03:28 UT (5 m 18 s), so every exposure (11:58:21-12:03:23) falls inside. The paper assumes second contact at 11:58:18 | - | - | PASS (8 s later than the smooth-limb prediction; the paper: contacts were "somewhat before those calculated") |
| 12 | Longitude implied by the 25 field altitudes printed on p.301 (field centre = Sun or 7-star centroid; time = printed G.M.T. or 1 min later): 2h 42m 48s to 2h 44m 06s W, rms 0.07° (4-inch) and 0.29° (astrographic) | 25 | - | INFO (R6) |
| 13 | tesseract 5.5 second opinion on the full-resolution table crops: 260/276 digit strings found. The 16 misses are OCR garbles ('O' for 0, '4' for +, merged tokens, a dropped line); the cells involved were re-checked at 2-3× zoom | 276 | - | PASS |

## Discrepancies

Printed values are kept in every CSV. None of these is a reading error.

**R1. p.309, star 6, R.A. "Observed": printed +0.10; the arithmetic gives -0.10.**

- Table VIII gives -0.03 (eclipse Mean) - (+0.07) (comparison Mean) = -0.10. The other 13 "observed" values are
  exactly these differences.
- The p.331 radial observation for star 6, 0.56, requires -0.10: projected, -0.10 gives 0.556 and +0.10 gives 0.578.
- The R.A. Means of Table VIII sum to 0.00 (eclipse) and +0.02 (comparison), so the p.309 R.A. column should sum to
  -0.02. It does with -0.10, but sums to +0.18 with +0.10.
- The plus sign is clear at full resolution (`crops/fr/z_p309_star6.png`).
- One sign change restores agreement with Table VIII and p.331, so this is a sign misprint on p.309. Alpha (from
  Table VII) and p.331 are unaffected; anyone plotting p.309 should use -0.10.

**R2. Table VIII, eclipse Decl., star 11, Mean: printed +0.03; the 7 printed values average +0.0357, i.e. +0.04.**

- The values are 0.00, -0.08, -0.03, +0.02, +0.17, +0.16, +0.01, all clear at full resolution.
- This is 0.0007″ beyond strict rounding. The unrounded entries could in principle average as little as +0.031, but
  the Table II rebuilds (checks 9 and 9b) give +0.036 to +0.037.
- p.309 uses the printed +0.03: Decl. observed for star 11 = +0.03 - (-0.13) = +0.16, as printed. p.331 (0.20) is
  insensitive to it (the projection factor is 0.06).
- Probably a 0.01″ slip in the 1920 mean. It cannot be resolved; kept.

**R3. p.307, plate 14_2b, "Scale Value. From x.": printed +0.008; the arithmetic gives +0.018.**

- sol_a (Table IV) +0.032 + a(14_2b) 0.409 - a(scale plate) 0.423 = +0.018. The other 55 values of the table follow
  this rule exactly. T1 noticed the same.
- The first table is not at fault. 14_2a and 14_2b are two exposures on one plate, both rows print a = 0.409, and
  0.409 fits the altitude sequence (check 4d).
- The adopted scale +0.013 for 14_2a/14_2b is the mean of +0.010, +0.004, +0.008 and +0.030, so +0.008 was used
  downstream. With +0.018 the adopted scale would be +0.0155.
- This is therefore a 1920 computing slip, not only a typesetting error.
- Effect on Table VIII for 14_2a and 14_2b: at most 6.25 × 0.0025 × 1.32 = 0.02″ in the scale terms. No effect on
  alpha.

**R4. Table VIII against Table II (cross-table, using T1's Table II).** Six single cells and one column do not
follow from Table II with the printed p.307 constants. The six cells and the 17_1 column were re-read at full
resolution, the doubtful ones at 2-3× zoom, and are printed unambiguously. These are 1920 differences between two printed tables. They do not affect Table VII's alpha, which
comes from the normal equations, but they do enter p.309 and p.331, which are built from Table VIII.

| | Cell | Printed | Table II implies (plate constants fixed by the other six stars) | Explanation |
|---|---|---|---|---|
| a | eclipse R.A., plate I, star 2 | +0.75 | +0.81 | Table II prints -0.789. T1 (D1) shows the normal equations used -0.798, a transposition; with -0.798 Table VIII is reproduced (+0.752). A **Table II misprint**; Table VIII is consistent. |
| b | comparison R.A., 14_2a, star 2 | +0.31 | -0.00 (+0.04 if the zero point is the mean of all seven) | Table II prints -0.396. T1 (D2) shows -0.346 was used; with it Table VIII is reproduced (+0.306). A **Table II misprint**. |
| c | comparison R.A., 14_2b, star 5 | +0.16 | -0.35 | The other six stars fit within 0.005″. Table IV agrees with Table II for this plate, so only the Table VIII computation used a different value, equivalent to Dx about +1.095 instead of +1.013. The column sums to -0.02 with +0.16, so the zero point was formed with it. A **1920 slip in Table VIII**. With the zero point re-formed as the column mean, the corrected 14_2b values would be -0.28 for star 5 and about +0.07 higher for the other six. The comparison R.A. Mean of star 5 would then be about -0.32 instead of -0.25, and its p.331 radial value about 0.50 instead of 0.54. |
| d | eclipse R.A., II, star 4 | +0.11 | +0.092 | Others fit within 0.002″; the column sums to 0.00 with +0.11. A 0.018″ slip (0.003 r in Dx). |
| e | eclipse R.A., III, star 6 | -0.09 | -0.127 | Others fit within 0.003″; the column sums to 0.00. A 0.037″ slip (0.006 r). |
| f | eclipse R.A., VIII, star 3 | -0.01 | +0.016 | Others fit within 0.004″. The column sums to -0.03; it would sum to -0.01 with +0.01 and 0.00 with +0.02. Either a sign misprint (-0.01 for +0.01) or a slip. Not resolvable; the printed Mean -0.09 fits both. |
| g | comparison R.A., whole column 17_1 | - | - | The residuals follow Table II with an x coefficient of -0.0694 r instead of the -0.060 r implied by p.307 (adopted scale -0.084 and a(17_1) - a(scale) = -0.024). The y coefficient agrees (-0.140). This gives deviations from -0.05″ to +0.09″, linear in x. It is as if a = 0.409 (the 14_2 value) or a scale of -0.093 had been used. Effect on the star means is at most 0.012″. |

The Decl. blocks follow Table II with the adopted constants within 0.013″ (eclipse) and 0.018″ (comparison, 18_2
star 11), i.e. within one to two rounding steps.

**R5. §39, p.330: "From right ascensions . . . 2″.06"; 0r.105 gives 2″.08.**

- With the paper's own factors (6″.25 per r, 50′/15′.8), 0r.105 gives 2″.077 and 0r.098 gives 1″.938, which matches
  the printed 1″.94.
- No single solar radius gives both printed values. 1″.94 needs 15′.75-15′.83; 2″.06 needs 15′.89-15′.97.
- 2″.06 would need alpha = 0.1040-0.1045 r. Solving the printed normal equations exactly gives 0.1046 r, i.e.
  2″.07.
- The mean 1″.98 follows either way: (2.06 + 2 × 1.94)/3 = 1.980, (2.077 + 2 × 1.938)/3 = 1.984, and 0.100 r gives
  1.978.
- Kept as printed: an unexplained 0″.02 (1 %) difference in the R.A.-only value.

**R6. p.296: "longitude 2h. 47m. 25s. west" is inconsistent with the paper's own timing.**

The figure is printed unambiguously ("47m", `crops/fr/z_p296_longitude.png`).

- JPL Horizons (topocentric Moon and Sun, 5 s steps, latitude as printed):
  - At 2h 47m 25s W (41.854° W; the Horizons run used 41.857° W, which changes the contacts by under 1 s, see
    Audit): second contact 11:55:54 UT, third contact 12:00:21, totality 4 m 27 s.
  - At the project's station, 40.35° W (2h 41m 24s): 11:58:10-12:03:28, 5 m 18 s.
- The paper assumes second contact at 11:58:18 civil G.M.T. Its 4-inch plates IV, V, VII and VIII, exposed between
  12:00:17 and 12:03:20, show the stars. With the printed longitude those plates would have been taken after totality.
- A one-digit misprint, "47" for "41", gives 2h 41m 25s = 40.354° W, which is the station.
- The field altitudes printed on p.301 point the same way. At the printed longitude the computed altitudes average
  0.99° below the printed ones (4-inch; 0.94° for the astrographic). At 2h 41m 25s they average 0.30° above (0.35°).
  A zero offset needs 2h 42m 48s to 2h 44m 06s W, depending on the assumed field centre and time tag (check 12).
- Kept as printed in the quotation in `dyson1920_sobral_conditions.md`. Use the 40.35° W station for any computation.

**R7. p.309 "Calculated" Decl. of star 4 is +0.74; Table I prints +0.73.**

- Both are roundings of about 0.735-0.737: K × Ey = (0.5501-0.5513) × 1.336 = 0.7349-0.7365 (check 7b), or 0.735
  from the Table I pattern (T1).
- A rounding-edge difference between two printings, not an error.

**R8. p.307: the adopted scale of the eclipse plates is 0.000, but their 14 printed scale values average -0.0026 r.**

- Unlike the comparison groups, it is not the group mean. No reason is printed. The difference is within the scatter
  of the individual values (±0.03 r).
- Note in passing: the eclipse plates' printed scales (mean -0.003) and the July 14 plates' (+0.013) differ by about
  0.016 r, and the July 15-18 plates differ from both by about 0.08-0.10 r. This is the "difference of scale value
  from July 15-18" that p.307 mentions.

**R9. Table VIII Decl. zero point.**

- The R.A. columns are zero-sum (zero point = column mean). The Decl. columns are not: they sum to about +3.50″
  (eclipse) and +1.32″ (comparison), nearly constant from plate to plate within each block (check 5b).
- The convention is not stated on pp.307-308. It cancels in p.309 (eclipse minus comparison).
- The constancy of these sums served as a check on the transcription.

## Ambiguities not resolved

- **R2** (eclipse Decl. star 11 Mean) and **R4f** (VIII, star 3, R.A.): the printed digits are certain, but the
  arithmetic cannot say whether a misprint or a 1920 slip is involved.
- **R5**: the source of 2″.06 is unknown.
- **R6**: the evidence that the printed longitude is wrong is strong. Whether the error was a typesetting slip ("47"
  for "41") or the expedition's own, the station is at about 2h 41m 24s W.
- Inferred labels: Table VII prints no plate labels (fixed by check 2). The p.331 table has a single, unlabelled
  group. Mapping 4-inch Ref. Nos. 1-8 to plates I-VIII rests on Table II-VIII omitting VI and on "No. 6 ... taken
  through cloud" (pp.300, 302).
- The absolute refraction and aberration coefficients of p.307 were not recomputed: the pressures and temperatures
  used in 1920 are not printed. Only their internal consistency (checks 4 and 4d) was tested.

## Audit (second reading)

Second, independent transcriber (blind double-keying, T2-audit, 2026-10-20). Working files are in
`data/scratch/T2-audit/` (git-ignored): `blind/` (own reading), `tools/compare.py`, `tools/checks.py`, `work/*.py`,
and their `*_output.txt` / `checks_*.txt`.

**Blind reading.** Every printed cell of the seven CSVs was read again from full-resolution JPEG 2000 crops (leaves 9,
15, 16, 17, 18 and 40; leaf 2 for the star x, y used in check H) before the first transcriber's value columns were
opened. Cell-by-cell comparison at the printed precision:

| File | Printed cells compared | Disagreements |
|---|---|---|
| `dyson1920_table7_alpha.csv` | 32 | 0 |
| `dyson1920_p307_refraction_aberration.csv` | 36 | 0 |
| `dyson1920_p307_scale_orientation.csv` | 84 | 0 |
| `dyson1920_table8_residuals.csv` | 224 | 0 |
| `dyson1920_p309_obs_calc.csv` | 28 | 0 |
| `dyson1920_p331_radial.csv` | 14 | 0 |
| `dyson1920_sobral_exposures.csv` (d, h, m, s, exposure, plate; printed blanks included) | 162 | 0 |
| **Total** | **580** | **0** |

The added columns also agree: the Table VII plate labels, the p.331 station column, and the derived exposure columns
(carried-down time, plate type, civil UT, plate label; the time after second contact is equal in value, "+5" vs "5").
The only difference was a convention: the blind copy carried d, h, m down, while the CSV keeps the printed blanks and
puts the carried-down time in `gmt_astronomical`. Both are faithful. No CSV value was changed.

The cells behind R1-R4 were viewed again by the second reader in tight full-resolution crops: p.309 star 6 R.A.
"+0·10"; Table VIII eclipse Decl. star 11 (0·00, -0·08, -0·03, +0·02, +0·17, +0·16, +0·01; Mean +0·03); p.307 14_2b
scale from x "+0·008"; Table VIII cells R4a-f (+0·75, +0·31, +0·16, +0·11, -0·09, -0·01) and the 17_1 R.A. column.
All are printed as transcribed; the discrepancies are in the 1920 printing or arithmetic.

**Checks re-run** with an independent script (`tools/checks.py`) on the final CSVs, and identically on the blind copy:

| Check | n | Max dev | Result |
|---|---|---|---|
| Format: printed decimals, sign on every non-zero value, ASCII | 418 | - | PASS |
| Table VII Mean row = mean of 7 rows | 4 | 0.00043 r | PASS |
| Table VII = printed alpha of Tables III-VI (T1 CSV), order I..VIII / 14_2a..18_2 | 28 | 0 | PASS |
| p.306: 0.120-0.015, 0.129-0.031, (0.105+2x0.098)/3 = 0.1003, x6.25 = 0.625″, x50/15.8 = 1.978″ | 5 | 0.0022″ | PASS |
| p.330: 0r.098 → 1.938″ (1″.94), 0r.105 → 2.077″ (2″.06) | 2 | 0.017″ | ORIGINAL (R5) |
| p.307 scale/orientation = solution + (row - scale plate) coefficient | 56 | 0.010 r | 55/56; 14_2b scale from x (R3) |
| Adopted orientation = mean of the two | 14 | 0.0005 r | PASS |
| Adopted scale = group means (+0.0130, -0.0844); eclipse 0.000 vs -0.0026 | 3 | - | PASS / R8 |
| p.307 first table: b = d, 14_2a = 14_2b | 9 | - | PASS |
| Table VIII Mean = mean of 7 plates | 28 | 0.0057″ | 27/28 (R2) |
| p.309 observed = eclipse Mean - comparison Mean | 14 | 0.20″ | 13/14 (R1) |
| p.309 calculated = Table I Sobral (T1 CSV) | 14 | 0.01″ | 13/14 (R7) |
| p.331 = p.309 projected away from the Sun (Sun fitted to all 26 Table I Sobral components: x = +0.4381, y = -0.2625, K = 0.5516″; also Sun at +0.44, -0.26) | 14 | calc 0.005″, obs 0.017″ | calculation 7/7, observation 6/7; star 6 gives 0.577 with +0.10 and 0.556 with -0.10 (R1) |
| Exposures: 19 x 5/10 s alternating, 8 x 28 s, no overlap (gaps 4-10 s and 10-11 s), +3/+5 s to +302/+305 s, derived columns, carry-down | 27 | - | PASS |

Extra checks (`work/extra_checks.py`, `work/rebuild_t8.py`) reproduce 3c (8.1 %, 5.9 %, 4.8 %), 7b (K = 0.5501-0.5513″,
r0 = 15′.72-15′.75; K = 0.553 gives 12/14), 8b (K/r misses stars 10 and 4 by +0.0057″ and -0.0080″) and check 9:
Table VIII rebuilt from Table II gives 165/196 cells within 0.0105″ with column-mean zero points (179/196 with median
zero points). Every larger deviation is one of R4a-g, or a Decl. cell of at most 0.020″ (18_2 star 11, I star 2).

The derived timing was recomputed independently with skyfield 1.55 and DE440s (smooth limbs, UT1 with ΔT = 21.2 s,
latitude as printed). At 40.35° W totality runs 11:58:11.5 to 12:03:26.5 UT (5 m 15 s); Horizons gave 11:58:10 to
12:03:28. At the printed 41.854° W it runs 11:55:56.0 to 12:00:20.5 (4 m 24.5 s); Horizons gave 11:55:54 to 12:00:21.
Every exposure (11:58:21 to 12:03:23) is inside totality at 40.35° W in both computations. The p.301 field altitudes
(field centre = the Sun at eclipse) give zero mean offset at 2h 42m 50s W (4-inch, printed time) to 2h 44m 04s W
(astrographic, time + 1 min). The mean offsets are -0.99° and -0.94° at the printed longitude and +0.30° and +0.36° at
2h 41m 25s W, with rms 0.08° and 0.29°. R6 and check 12 are reproduced.

**Corrections made (text only):**

1. R6, and §1 of `dyson1920_sobral_conditions.md`: 2h 47m 25s W is 41.854° W, not 41.857° W. The Horizons run used
   E-longitude 318.14306° (41.857° W, a 10″ slip in the conversion). This moves the contacts by under 1 s, so the
   times stand.
2. R7: "0.735-0.739" corrected to "0.735-0.737" (K × Ey = 0.7349-0.7365).
3. `dyson1920_sobral_conditions.md`: the station coordinates are in `fetch_data.py`, not `config.py`. "Father Cortie's
   lens" is replaced by the paper's own words (p.295, p.296). The local mean time of 08h civil G.M.T. at 40.35° W is
   05h 19m, not 05h 25m. Leaf 4 (p.295) is added to the list of scans.
4. Headers of `dyson1920_table7_alpha.csv` and `dyson1920_p307_refraction_aberration.csv`: the sentence "The satisfactory
   accordance ... calculated from the times and dates of exposure." is on p.306 (it was cited as pp.306-307). Comment
   lines only; the data rows are byte-identical.

**Noted, not changed:** p.312 prints "α + 243e = + 0r·0435" with no decimal point before 243 (checked at full
resolution). The paper's own numbers need 0.243: 0.0435 - 0.243 × 0.082 = 0.0236 → 0r.024, and 0.0435 - 0.243 ×
0.020 = 0.0386 → 0r.039. The quotation in `dyson1920_sobral_conditions.md` stays as printed.
