# Sobral 4-inch tables of Dyson, Eddington & Davidson (1920): transcription notes

Source: Dyson, F. W., Eddington, A. S. & Davidson, C. 1920, *Phil. Trans. R. Soc. A* 220, 291-333,
"A Determination of the Deflection of Light by the Sun's Gravitational Field, from Observations made at
the Total Eclipse of May 29, 1919". Scans: archive.org item `philtrans06337895`. Leaf N is printed page
291 + N. Every number was read on the full-resolution 5120 x 6600 JPEG 2000 leaves
(`06337895_jp2.zip/06337895_jp2/06337895_00NN.jp2`): leaf 2 (p.293), 12 (p.303, printed sideways),
13 (p.304) and 14 (p.305). Leaves 15-17 (pp.306-308) were read only to cross-check.

## Files

| File | Printed source | Content |
|---|---|---|
| `dyson1920_table1_stars.csv` | Table I, p.293 | 13 stars: photographic magnitude, x, y (unit 50'), predicted displacements at Sobral and Principe (arcsec) |
| `dyson1920_4inch_condition_coeffs.csv` | section 19, p.304 | x, y, Ex, Ey of the 7 stars on the 4-inch plates (left-hand sides of the equations of condition) |
| `dyson1920_table2_4inch.csv` | Table II, p.303 | Dx, Dy (revolutions) of 7 stars on 7 eclipse and 7 comparison plates, long format (98 rows) |
| `dyson1920_table2_offsets.csv` | Table II, p.303 | the number printed below the line in each column (14 plates), with notes |
| `dyson1920_tables3to6_normal_equations.csv` | Tables III-VI, pp.304-305 | right-hand sides rhs1-rhs10 and the 3 printed solutions, per table and plate (364 rows) |
| `dyson1920_tables3to6_lhs.csv` | Tables III-VI, pp.304-305 | the printed left-hand-side coefficients, once per table (40 rows) |

Conventions:

- Every value is exactly as printed: explicit sign, printed number of decimals. Old typesetting is
  written plainly (`1".98` becomes 1.98, the raised decimal point becomes `.`). No printed number has
  been changed. Where the paper's own arithmetic shows a misprint or a slip, the printed value is kept
  and the case is documented below and in the CSV header.
- Units: Table I x, y in units of 50 arcmin, x increasing to the east (R.A.) and y to the north
  (Decl.). The displacements are in arcsec on the same axes. Table II and Tables III-VI are in
  micrometer revolutions ("r"); 1 r = 6.25 arcsec ("The value of the micrometer screws (both in R.A.
  and Decl.) is 6".25", p.302). The orientation was confirmed by the Gaia fit (check A6): Table I x, y
  match tangent-plane coordinates with x east and y north, rotation +0.016 deg and scale 1.000. For
  example, kappa1 Tau lies 5.6' north of and 0.7' west of kappa2 Tau, which Table I gives as y +0.112
  and x -0.014 units.
- Plate labels: eclipse I II III IV V VII VIII (plate VI showed only prominences, p.302). The
  comparison plates are printed with subscripts 14_{2a} 14_{2b} 15_1 15_2 17_1 17_2 18_2 and are
  written `14_2a` ... `18_2`. Star rows follow the printed order 11 5 4 3 6 10 2.
- Model (p.304): `ax + by + c + alpha Ex (= Dx)` and `dx + ey + f + alpha Ey (= Dy)`. On p.304 the R.A.
  rows are printed in the order c, b (times y), a (times x), alpha, and the Decl. rows in the order
  f, d (times x), e (times y), alpha. Tables III-VI keep these orders, so the unknowns are (c, b, a,
  alpha) for R.A. and (f, d, e, alpha) for Decl.
- `quantity` in the normal-equations file:
  - rhs1-rhs4: the 4 normal equations.
  - rhs5-rhs7: the reduced equations after eliminating c or f.
  - rhs8-rhs9: after eliminating b or d.
  - rhs10: the final equation in alpha.
  - The solutions in printed order: sol_alpha, sol_a, sol_b (III, IV) and sol_alpha, sol_e, sol_d
    (V, VI).
  - Row numbers in the LHS file use the same numbering.
- Tables V and VI print no column headings. Their columns are in the same plate order as Tables III
  and IV, as shown by rebuilding each column from Table II and by Table VII (p.306).
- Offset convention, confirmed from the numbers: value used = printed value - offset. The footnote,
  quoted exactly (the minus signs are long dashes): "* The numbers -1.500, -0.554, &c., given below
  the line, were taken out to make the values of Dx, Dy small and positive for arithmetical
  convenience."
  - The asterisk is printed after the first offset of each half of the table (plate I Dx -1.500* and
    plate 14_2a Dx -0.552*).
  - Under this convention the printed rhs1 equals the column sum (section C below).
  - The offsets are:
    - the column minimum for 12 of the 14 Dy columns. 14_2a is 0.002 lower; VIII is 0.010 higher,
      see D3.
    - the column minimum for 5 of the 7 comparison Dx columns. 14_2a is 0.008 lower and 15_2 is
      0.002 lower.
    - round numbers below the minimum for the eclipse Dx columns (-1.500, +0.500, +0.400, -1.300).

## Checks run

Scripts are in `data/scratch/T1/` (git-ignored). `verify_t1.py` reads only the six CSV files. The
companion scripts are `gaia_xy_check.py`, `sun_check.py`, `ocr_check.py` and `rereduce.py`.

Results:

- **PASS** means within the stated tolerance.
- **ORIGINAL** means the deviation comes from the 1920 printing or arithmetic. The cause is shown and
  the printed value is kept.

| # | Check | n | max dev | Result |
|---|---|---|---|---|
| 0 | Format: explicit sign and printed precision in every numeric cell; row counts 13/7/98/14/364/40 | 787 | - | PASS |
| A1 | Table I Sobral displacements = k (s - s0)/\|s - s0\|^2 (fit x0, y0, k), all 26 values | 26 | 0.059" | rms 0.012"; one outlier (star 13 dx) |
| A2 | same without star 13 dx: x0 +0.4364, y0 -0.2627, k 0.5504" | 25 | 0.0048" | PASS (rms 0.0031"; rounding alone gives 0.0029") |
| A3 | Table I Principe, all 26 values | 26 | 0.020" | rms 0.0050"; one outlier (star 8 dx) |
| A4 | same without star 8 dx: x0 +0.5524, y0 -0.2456, k 0.5500" | 25 | 0.0073" | PASS (rms 0.0032") |
| A5 | all 13 displacement vectors point away from the fitted Sun (min cos 0.991 / 0.999) | 26 | - | PASS |
| A6 | Table I x, y vs Gaia DR3 propagated to 1919.4 and projected (tangent point, rotation, scale fitted; project's `stars.py`) | 26 | 0.0009 units (2.9") | PASS (rms 1.66"; every coordinate within the 0.001-unit printing step of Gaia) |
| A7 | coordinates sane: all \|x\|, \|y\| <= 1.321 units (66'), no star outside the field | 26 | - | PASS |
| A8 | Sun shift Sobral to Principe from the two fits: +5.80' east, +0.86' north. JPL Horizons, Sun at each mid-totality in the same frame: +5.12', +0.90' (`sun_check.py`). Fitted Suns vs Horizons: Sobral -0.29', -0.08'; Principe +0.39', -0.12' | - | - | plausible (direction and size agree; 1919 Sun places good to about 0.4') |
| B1 | x, y on p.304 = Table I (both printings per star agree) | 14 | 0 | PASS |
| B2 | (Ex, Ey) = (s - s0)/\|s - s0\|^2 with one Sun position x0 +0.4361, y0 -0.2627 (Table I Sobral fit: +0.4364, -0.2627) | 14 | 0.0014 | PASS |
| B3 | k_Sobral x (Ex, Ey) vs Table I Sobral displacements of the 7 stars | 14 | 0.0053" | PASS |
| C1 | value - printed offset >= 0 ("small and positive") | 196 | - | ORIGINAL: one exception, plate VIII Dy star 5 = -0.010 (offset misprint, D3) |
| C2 | printed rhs1 = sum(Table II value - offset), exactly | 28 columns | - | 24/28 exact; the other 4 are D1, D2, D3 and D9 |
| D | LHS of Tables III-VI rebuilt from p.304 (normal equations and Gauss reduction) | 80 | 0.042 (3.683) | ORIGINAL: D6, D7; all others within 0.0015 |
| E1 | rhs1-rhs10 and 3 solutions rebuilt from Table II as printed | 364 | 0.30 | ORIGINAL: 287/364 within 0.002 |
| E2 | same with the 3 Table II values the 1920 reduction used (D1-D3) | 364 | 0.30 (D4 sign) | 314/364 within 0.002; all 50 others traced (D4, D5, D8) |
| F | each printed rhs5-rhs10 and solution recomputed from the printed lines above it | 252 | 0.30 | 246/252 within 0.002 (exceptions: D4 [4 values], D5, and IV 15_2 sol_b 0.0022) |
| G | search over single-digit changes, adjacent transpositions and sign flips (7 values, offset, rhs1-4) for every column failing E1 | 14 columns | - | finds D1 and D2 uniquely; see D8 |
| H | rhs2-rhs4 as a hand computer forms them (each product rounded to 3 decimals), using the D1-D3 values | 84 | 0.009 | 79/84 within 0.001 (45/84 exact vs 28/84 for the rounded exact sum); the 5 others are the rhs4 slips D8. With Table II exactly as printed: 71/84 (41 vs 25) |
| I1 | printed alpha of Tables III-VI = Table VII (p.306), and its 4 means +0.120 +0.015 +0.129 +0.031 | 32 | 0 | PASS |
| I2 | final result from the printed alphas: R.A. +0.1050 r, Decl. +0.0977 r, weights 1:2 give +0.1001 r = 0.626" at 50', or 1.981" at the limb (paper: +0.105, +0.098, +0.100 r = +0.625", 1".98) | - | - | PASS |
| I3 | p.307 scale value and orientation from x and y = printed a, e, b and -d plus the p.307 refraction/aberration differences | 56 | 0.010 | PASS except one p.307 value (plate 14_2b scale from x printed +0.008; the arithmetic gives +0.018) |
| I4 | Table VIII (p.308) R.A. residuals rebuilt from Table II with adopted scale/orientation and refraction, zero point = mean (12 plates) | 84 | 0.032" | PASS (rms 0.0066"; Table VIII is rounded to 0.01") |
| J | tesseract second opinion on all numeric blocks | 666 tokens | - | 626 agree in digits. The other 40 are in 3 rows where tesseract dropped or garbled a number, plus 5 misread digits (re-read visually: -1.432, -0.317, +1.319, -0.104, '+ .107'). The 40 sign disagreements where the digits agree are all printed '+' read as '-'; the signs are confirmed by checks C2, E and F |

Re-reduction (`rereduce.py`): recomputing all 28 plate solutions exactly from Table II gives R.A.
+0.1048 r, Decl. +0.0977 r, combined +0.1001 r, and **1.980" at the limb** when the three values used
in 1920 (D1-D3) are applied. The paper prints 1.98". With Table II exactly as printed the R.A. value is
+0.1097 r and the limb deflection **2.012"**: only D1 and D2 matter, since the offset (D3) is absorbed
by f. The 1920 arithmetic slips below change single plates by up to 0.0035 r but cancel in the means.

## Discrepancies and their resolution

All the values below were re-read at full resolution. Each is printed unambiguously, so none is a
misreading on our side.

**D1. Plate I, star 2 (Piazzi IV 82), Dx: printed -0.789; the reduction used -0.798.** This is the
"2.168 vs 2.159" case.

- Table II (p.303) prints -0.789. The plate I Dx column sums to -8.332, and after the -1.500 offset
  rhs1 = 2.168. Table III prints rhs1 = +2.159.
- Rebuilt from the printed value, the whole of Table III column I deviates (printed - rebuilt): rhs1
  -0.0090, rhs2 +0.0037, rhs3 -0.0092, rhs4 -0.0137, ..., alpha -0.0044.
- These deviations match -0.009 x (1, y, x, Ex) of star 2 = (-0.009, +0.003, -0.010, -0.014) to
  within 0.0008. Swapping the last two digits
  (-0.789 to -0.798) brings all 13 printed numbers of column I within 0.0013 (rhs1-rhs4 within
  0.0008). No other single change does this.
- Table VIII (p.308) is a separate computation from the same measures. Its plate I R.A. residual for
  star 2 is +0.75". Rebuilt, it is +0.800" with -0.789 and +0.752" with -0.798.
- So both downstream computations used -0.798. The printed Table II value is most likely a
  typesetting transposition. It is not a misread here, and the normal equations themselves were
  formed correctly. The CSV keeps -0.789.
- Effect: alpha(I, R.A.) = +0.0972 r with -0.798 (printed +0.098), or +0.1024 r with -0.789.

**D2. Plate 14_2a, star 2, Dx: printed -0.396; the reduction used -0.346.**

- Table IV column 14_2a deviates from Table II as printed by up to 0.078 (rhs1 +0.051).
- One digit changed (-0.396 to -0.346) brings all 13 numbers within 0.0010. The alternative -0.345
  fits rhs1 exactly but rhs4 only to 0.0019.
- Table VIII's plate 14_2a star 2 residual of +0.31" is reproduced by -0.346 (+0.306") and not by
  -0.396 (+0.038").
- The CSV keeps -0.396.

**D3. Plate VIII, Dy offset: printed -1.322; -1.332 was used.**

- The Dy column's minimum is star 5, at -1.332. Every other Dy offset except 14_2a (0.002 lower)
  equals its column minimum.
- With -1.322, star 5 becomes -0.010, which contradicts "small and positive".
- Table V rhs1 = +1.906 = 1.836 + 7 x 0.010. With -1.332, rhs1-rhs3 agree to 0.0017 and all 13 numbers
  to 0.0036.
- Probable misprint of the offset (a single digit). It does not affect alpha. Kept as printed, with a
  note in the offsets file.

**D4. Table VI, plate 17_2, rhs5: printed +0.151; should be -0.151.**

- From the printed rhs2 and rhs1: 0.041 - (1.787/7) x 0.752 = -0.151. Table II gives -0.1507.
- Every later printed line follows -0.151: rhs8 -0.219 (+0.151 would give -0.357), rhs9 -0.019
  (-0.101), sol_d -0.010 (+0.073). The p.307 orientation from y, +0.029, uses d = -0.010.
- This is a sign misprint. Kept as printed.

**D5. Table III, plate V, sol_b: printed +0.042; the printed lines give +0.066.**

- Back-substituting (0.533 - 1.666 x 0.128 - 0.281 x 0.140)/4.271 gives 0.0657. Table II gives
  0.0662.
- The printed +0.042 follows exactly if 0.433 is used for rhs5 (printed and correct: +0.533). This
  is a slip in the 1920 back-substitution.
- p.307 (orientation from x, +0.008 = 0.042 - 0.034) carries it. b does not enter alpha.

**D6. Tables III and IV, LHS row 6 coef3: printed +3.683; should be 3.638.**

- 4.094 - 1.787^2/7 = 3.638. The next printed row uses 3.638: 3.638 - 1.666^2/4.271 = 2.988, which is
  printed.
- Tables V and VI print the same element as "+3.638d". This is a digit transposition, present in both
  R.A. tables.

**D7. Smaller LHS inconsistencies, all kept as printed.**

- Sum of x^2 is printed 4.094 in all four tables, but the p.304 x values give 4.0968. Sum of x =
  1.787 is exact, so no single x value explains it. The reduced rows follow 4.094, so the 1920
  computers used it.
- Sum of y Ex is printed 0.335 (exact 0.3343). Sum of x Ex is +2.534 in Table III but +2.535 in Table
  IV (exact 2.5344). Sum of y Ey is 3.694 (3.6920) and sum of x Ey is 1.840 (1.8392).
- The Decl. reduction rounds slightly at each step, ending at row 10 = 2.224 against 2.2286 from the
  exact sums.
- Effect on alpha: below 0.001 r.

**D8. Slips in rhs4 of the Decl. tables (1920 arithmetic).**

- The cases (printed minus the sum of Ey x value):

  | Column | Printed rhs4 | Sum of Ey x value | Difference | Effect on alpha |
  |---|---|---|---|---|
  | V I | +2.657 | 2.666 | -0.009 | -0.0034 r |
  | V II | +1.681 | 1.684 | -0.003 | -0.0006 r |
  | VI 15_1 | +0.220 | 0.229 | -0.009 | -0.0034 r |
  | VI 17_2 | +0.296 | 0.301 | -0.005 | -0.0022 r |
  | V VIII (after D3) | +1.707 | 1.709 | -0.002 | - |

- In each case rhs1-rhs3 agree with Table II, and the later printed lines follow the printed rhs4
  (stepwise check F). So no Table II value can cause the deviation. A wrong Ey would also break the
  printed sum of Ey = 3.316, which is exact.
- These are independent slips in 1920 (all negative), kept as printed.
- Accumulated rounding of 0.002-0.005, consistent with hand computation (products rounded to 3
  decimals, check H, then each reduction step rounded):
  - III VII: rhs3 +0.0027; with rounded products +0.001. The reduction carries it to alpha -0.0035.
  - III VIII: rhs5, rhs8 0.002.
  - III V: rhs4 -0.0024 (with rounded products -0.001), carried to rhs10 -0.0030; its sol_b is D5
    (added in the audit).
  - IV 15_2: rhs3 +0.0021.
  - IV 17_1: rhs3 +0.0018, rhs6 +0.0027.
  - V IV: rhs4 -0.0026.
  - VI 18_2: rhs8 -0.0021.
  - IV 15_2 sol_b: 0.0022 from its own lines.

**D9. Table VI, plate 14_2a, rhs1: printed +0.446; Table II gives 0.448.**

- With hand-rounded products, rhs2-rhs4 agree to 0.001.
- Several single changes to Table II would also bring rhs1-rhs4 within 0.001, for example star 2
  Dy -0.182 to -0.184. The printed digits are unambiguous, though, and the Table VIII Decl. pattern
  slightly favours -0.182.
- The later lines cannot tell 0.446 from 0.448. Only rhs7 = 0.168 matches 0.448 slightly better.
- Either a 0.002 slip or a misprint in rhs1, or a 0.002 error in one Table II value. Not resolvable.
  Effect on alpha: about 0.0004 r.

**D10. Table III, plate III, sol_alpha: printed "+ .107" without the leading zero.**

- Written +0.107 in the CSV. Table VII prints +0.107.

**D11. Table I, star 13, sobral_dx: printed -0.17; the radial pattern of the other 25 Sobral values
gives -0.110.**

- The residual is -0.060", 20 times the rounding scatter. Principe gives -0.14, and the Sobral dy
  (+0.40) fits.
- A single digit (-0.11) would fit, so this is a misprint or 1920 slip. Kept as printed. Star 13 is not
  on the 4-inch plates.

**D12. Table I, star 8, principe_dx: printed -0.29; the pattern gives -0.310 (residual +0.020").**

- Star 8's dy also deviates by -0.007".
- Most likely a small 1920 computing error. Kept as printed. Star 8 is not on the 4-inch plates.

**D13. Name: Table I star 11 is printed "66 Tauri" (kept as "66 Tau").**

- The star at these co-ordinates is 56 Tau = HD 27309; 66 Tau is about 12 degrees away. This is the
  project's Gaia identification in `config.py` (`TABLE1_GAIA`), confirmed by A6.

## Observations on other pages (not part of T1, for whoever transcribes them)

- p.307: plate 14_2b "Scale value from x" is printed +0.008, but sol_a +0.032 + 0.409 - 0.423 =
  +0.018. The other 55 values of that table follow the printed solutions exactly.
- Table VIII (p.308), R.A., plate 14_2b, star 5: printed +0.16" while Table II implies about -0.36".
  The printed mean -0.25 was computed with +0.16.
- Table VIII, R.A., plate 17_1: the residuals deviate from a rebuild linearly in x (-0.08" to +0.05"),
  as if a different scale had been applied.
- Table VIII, R.A., plate III, star 6: +0.03".
- Table VIII, R.A., plate VIII, star 3: printed -0.01", rebuilt +0.017" (-0.027"; added in the audit).
- Table VIII Decl.: the zero point is not the column mean. The columns sum to +3.5" (eclipse) and
  +1.3" (comparison), but the star-to-star pattern matches Table II to about 0.01" (0.02" at worst).

## Ambiguities not resolved

- For D1 and D2 the arithmetic proves which value the 1920 reduction used, but not which number was
  actually measured. A typesetting error in Table II is the simpler explanation, because two
  independent downstream computations agree on the alternative.
- For D9 we cannot tell whether 0.446 is a misprint or slip in rhs1, or whether one Table II value is
  0.002 off; the effect is negligible.
- For D11 and D12 the correct 1920 value cannot be recovered beyond the radial-pattern prediction.

## Audit (second reading)

An independent second transcriber double-keyed these tables blind on 2026-10-20. The working files are
in `data/scratch/T1-audit/` (git-ignored): `blind/` holds the second reading, `tools/` the scripts and
`crops/` the image crops.

**Method.**

- The full-resolution JPEG 2000 leaves 2, 12, 13 and 14 were downloaded again.
- Every numeric cell of the six CSV files was read from 1:1 crops into `blind/`, magnified 3-5x
  wherever a glyph was in doubt.
- This reading was completed before the value columns of the first transcriber's files were opened.
  Only their `#` headers and key columns were used, to align rows.
- `tools/compare.py` compared the two readings cell by cell.
- `tools/audit_checks.py` and `tools/audit_checks2.py` re-ran every internal check using only the six
  CSV files. They are independent of the scripts in `data/scratch/T1/`.

**Comparison.**

- 787 numeric cells were compared: Table I 91, p.304 coefficients 28, Table II 196, offsets 28,
  Tables III-VI 364 and LHS 80. The 80 empty LHS cells and the 13 names were compared as well. Row sets
  and row order agree.
- Result: **0 disagreements** and 0 format-only differences. No data cell was changed.
- These cells were zoomed at full resolution because they carry a documented discrepancy or a print
  blemish:
  - Table II: plate I star 2 Dx -0.789 (D1), 14_2a star 2 Dx -0.396 (D2), the VIII Dy offset -1.322
    (D3), and the 14_2a Dy column with its offset (D9).
  - Table II, 17_2 Dy star 2, -0.328: specks beside the minus sign, but at 3x it is a full-length dash.
  - LHS: +3.683 in Tables III and IV (D6), +2.534 (III) against +2.535 (IV), and +4.094.
  - Table I, star 11 magnitude: broken type. The second digit has lost its stem, but its flat top bar
    is that of a 5 (a 3 in this fount has a round top), so it reads 5.5, as transcribed. No arithmetic
    constrains the magnitudes.

**Checks re-run on the transcribed files.** Every number quoted in the sections above was reproduced,
except the tesseract tallies of check J, which were not re-run. Two values came out slightly
different:

- A6: the maximum deviation is 2.8" (the table above gives 2.9").
- A8: Sobral gives -0.07' instead of -0.08' (see the A8 item below).

- A1-A5:
  - Sobral, all 26 values: rms 0.0120", max 0.059".
  - Sobral without star 13 dx: x0 +0.4364, y0 -0.2627, k 0.5504", rms 0.0031". The model value for
    star 13 dx is -0.110.
  - Principe, all 26 values: rms 0.0050".
  - Principe without star 8 dx: x0 +0.5524, y0 -0.2456, k 0.5500", rms 0.0032". The model value for
    star 8 dx is -0.310, and star 8 dy deviates by -0.007".
  - Minimum cos: 0.991 (Sobral) and 0.999 (Principe).
- A6, done independently:
  - Inputs: an own Gaia DR3 TAP query of the 13 `config.py` source_ids, moved to epoch 1919.41 with
    astropy, transformed to FK5 equinox J1919.41, and projected gnomonically with the tangent point,
    rotation and scale fitted.
  - Result: rotation +0.016 deg, scale 1.0000, rms 1.66", max 0.0009 units.
  - The kappa1/kappa2 Tau example (5.6' north, 0.7' west) is confirmed.
- A8, done independently with JPL Horizons (topocentric astrometric Sun at mid-totality, projected in
  the Gaia-fitted frame):
  - Sobral: p.300 gives second contact as "May 28, 23h. 58m. 18s. G.M.T.", i.e. May 29 11:58:18 UT.
    Mid-totality was taken as 12:00:55 UT; this time is the auditor's assumption.
  - Principe: p.314 gives "2h. 13m. 5s. to 2h. 18m. 7s. G.M.T.", so mid-totality is 14:15:36 UT.
  - Shift from Sobral to Principe: +5.12', +0.90'.
  - Fitted minus Horizons: Sobral -0.29', -0.07'; Principe +0.39', -0.12'.
- B1-B3 as stated. In addition, the printed sums of Ex (0.226) and Ey (3.316) are exact. Sum Ex^2 =
  3.1417 and sum Ey^2 = 5.7845 round to the printed 3.142 and 5.784.
- C1 and C2 as stated: rhs1 is exact in 24 of 28 columns, and the other four are D1, D2, D3 and D9.
- D (LHS):
  - Rebuilt from the p.304 coefficients, only the D6 and D7 elements and their carried consequences
    deviate by more than 0.0015.
  - Rebuilt stepwise from the printed rows, only two values deviate: +3.683 (D6), and row 10 of
    Tables V and VI (3.941 - 2.453^2/3.508 = 2.2257 against the printed 2.224).
- E1 287/364, E2 314/364, F 246/252.
  - Of the 252 stepwise values, 199 are within 0.0005, 237 within 0.0010 and 245 within 0.0015.
  - All 50 E2 failures fall under D4, D5 and the D8 list. III V rhs4 and rhs10 have now been added to
    D8.
- G: the single-change search again finds only D1 (-0.798) and D2 (-0.346) for the large deviations.
  - Columns that fail by only 0.002-0.004 also accept trivial last-digit changes, for example V II
    star 4 Dy -0.944 to -0.946.
  - These are not evidence of misprints: rhs1 of those columns agrees exactly with Table II as
    printed.
- H: 79/84 with the D1-D3 values (45 exact against 28). With Table II exactly as printed it is 71/84;
  the table row now says this.
- D1-D13: the numbers are as stated. The effect of D7 on alpha is at most 0.0002 r.
- D13, checked in SIMBAD:
  - 66 Tau is r Tau, at Dec +09 27', about 12.3 deg from 56 Tau (+21 46').
  - 56 Tau's J2000 offset from kappa1 Tau is -80.2', -30.9'; Table I gives -79.75', -31.60'. The
    difference is the 0.44 deg rotation of the equinox-1919.41 frame plus proper motion.
  - In the 1919.41 Gaia fit above, star 11 fits within 0.0008 units.
- I1-I4 and the re-reduction:
  - Table VII and the p.306 text were read by the auditor: alpha = +0r.105 (R.A.), +0r.098 (Decl.),
    +0r.100 = +0".625; the sun's radius 15'.8; deflection at the limb 1".98.
  - Re-reduction: 2.012" with Table II as printed, 1.980" with the D1-D3 values.
  - p.307: 55 of 56 values agree.
  - Table VIII R.A.: 84 values, rms 0.0066", max 0.032".
- Table VIII Decl. (Observations, and the remark in D9), with p.308 read by the auditor:
  - The columns sum to +3.49" to +3.53" (eclipse) and +1.29" to +1.34" (comparison).
  - After removing a constant from each column, the star-to-star pattern matches a rebuild (adopted
    scale and orientation, refraction differences, VIII offset -1.332) with rms 0.0054" and max
    0.018".
  - For 14_2a star 2, -0.182 fits better than -0.184: pattern rms 0.0048" against 0.0058".
- Quotes checked against the page, all faithful:
  - the p.303 footnote;
  - the Table I text on p.293;
  - p.302: the micrometer value 6".25, sections 16-17 and plate No. 6;
  - p.304: the section 19 model and the printed rows;
  - the table titles on pp.304-305.

**Edits made by the audit.** These are documentation only; no value was changed.

1. `dyson1920_table2_4inch.csv` header: it said "Dx, Dy = (plate measured) - (scale plate) ...
   (sections 16-18, p.302)", but the paper does not state the sense. The header now marks it as
   inferred from alpha > 0 on the eclipse plates.
2. This file, check H: now states that the counts use the D1-D3 values, and gives the counts for
   Table II as printed.
3. This file, D8: added III V (rhs4 -0.0024, rhs10 -0.0030).
4. This file, Observations: added Table VIII R.A. plate VIII star 3 (-0.027").
