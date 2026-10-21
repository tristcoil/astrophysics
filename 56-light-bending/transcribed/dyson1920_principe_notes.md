# Dyson, Eddington & Davidson 1920: the Principe expedition, transcription notes

Source: F. W. Dyson, A. S. Eddington & C. Davidson 1920, "A Determination of the Deflection of Light by the
Sun's Gravitational Field, from Observations made at the Total Eclipse of May 29, 1919",
*Phil. Trans. R. Soc. A* 220, 291-333. Section IV, "The Expedition to Principe" (paras 26-38, printed
pages 312-330) and the General Conclusions (para 39, which begins on p.330; its Principe sentences are on
p.331). Scans: archive.org item `philtrans06337895`, leaf N = printed page 291 + N.

## Files

| file | content | printed page |
|---|---|---|
| `dyson1920_principe_exposures.csv` | List of Plates: check field (10 rows) and eclipse field (21 rows) | 315, 316 |
| `dyson1920_table11_check_measures.csv` | Table XI, check-plate measures Delta x, Delta y (11 stars x 5 pairs) | 318 |
| `dyson1920_table12_check_residuals.csv` | Table XII, check-plate residuals (arcsec) + Mean columns | 318 |
| `dyson1920_table13_check_constants.csv` | Table XIII, check-plate constants, uncorrected and corrected, + Mean | 319 |
| `dyson1920_principe_plates.csv` | The detailed reductions X-G1, X-H1, W-D1, W-I2, U-K2 (50 table rows, long format) | 323-326 |
| `dyson1920_principe_plate_constants.csv` | Every number printed in the running text: refraction, scale, theoretical constants, orientations, normal equations, delta-kappa, deflections, error budget, final result (132 rows, with the printed sentences) | 318-331 |
| `dyson1920_p328_principe_residuals.csv` | The untitled residual table of p.328 (arcsec) | 328 |

No printed value was changed anywhere. Every apparent misprint is kept as printed and listed under
"Discrepancies" below.

## How the numbers were read

* Every table was read cell by cell on enlarged crops of the 5120 x 6600 JPEG-2000 scans (leaves 24, 25, 27,
  28, 30-37), with tesseract as a second opinion. Table XI: 129 of 132 cells agree with tesseract; Table XII: all
  cells agree except eight of the star-8 row that tesseract garbled; the reduction tables of pp.323-326 and the
  p.328 table agree except tesseract's typical 3/8 and 1/7 confusions ("4814" for 4314, "4852" for 4352, "4889"
  for 4389, "187" for 137, "-15" for -75, "T" for 7). Every disagreement was settled on 2-4x zooms and is
  confirmed by the arithmetic below.
* The reduction tables are checked by their own arithmetic: nearly every printed number takes part in at least
  one identity that holds exactly (or to the 1-unit rounding the authors used), so a misread digit would show.
  Table XI is checked digit by digit against Table XII through the Table XIII constants (rms 0.0044"), so no tens
  or higher digit of Table XI can be misread. This agrees with the earlier reader's finding that Table XI is
  consistent with a second reading.
* The text numbers (constants, normal equations, error budget) were read on 0.4x crops of the full-resolution
  scan, which are perfectly legible, and checked against the OCR text layer.

## The reduction method, in our words

### Units (p.318, p.322)

* Plate coordinates x, y are in units of **5 mm, about 5'** (one reseau interval). The plate centre is near
  x = 14, y = 14. The eclipse-star coordinates are exactly x = 14 + 10 x(Table I), y = 14 + 10 y(Table I),
  because Table I's unit is 50' = 10 reseau intervals (checked for all 12 values).
* The measured differences Delta x, Delta y and all intermediate columns are in **units of the fifth decimal
  place** of that coordinate unit: 1e-5 x 5' = 1e-5 x 300" = **0".003**. The plate constants a, b, d, e, the
  orientation theta and kappa are expressed in the same 1e-5 units ("in units of the fifth place of decimals").
* **kappa** is the deflection at 50' (10 reseau intervals) from the Sun's centre, in units of 0".003. With the
  Sun's radius 15'.78, the full Einstein deflection 1".75 at the limb is 1".75 x 15.78/50 = 0".552 at 50' =
  184.1 units, so **kappa = 184** is the full Einstein value (p.322). A star's displacement in y is
  kappa alpha_y, with alpha_y = 10 (y - y_sun)/r^2 (r in plate units), and likewise Ex = 184 alpha_x. The printed
  alpha_y and Ex values fit this 1/r law with a single Sun position, x = 19.53, y = 11.53, to within their
  rounding (12 values, max 0.006 in alpha), and 0.003 Ey, 0.003 Ex agree with Table I's Principe column
  (p.293) to within its rounding (max 0.008").
* A change delta_kappa gives a limb deflection of 1".75 (184 + delta_kappa)/184.

### Step 1: the check field fixes the scale (paras 31-33, pp.317-320)

Because the eclipse happened in the afternoon, the eclipse field could not be photographed at night from
Principe at the same position in the sky. Instead a **check field** around Arcturus (R.A. 14h 12m 47s,
Dec. +20 deg 30', nearly the eclipse field's declination) was photographed at Oxford in January 1919 (plates
a1-e1, direct) and at Principe in May (q1, r1, s2, v1, w1, by reflection in the coelostat) at the altitude of
the eclipse (Z.D. 43-46 deg), with the same object-glass and 8-inch stop (p.314); the Oxford plates were taken
direct, without the coelostat, and at Principe the same part of the coelostat sector was used for the check field
as for the eclipse (p.313, p.317). An Oxford and a Principe plate were
laid film to film (the mirror image makes the stars nearly coincide) and the small differences Delta x,
Delta y (Principe - Oxford) of 11 check stars (the check field's own numbering, not Table I's) were measured
by Eddington on the Cambridge machine, 8 settings per image, half with the field reversed (Table XI).

For each of the five pairs a linear plate model Delta x = a x + b y + c, Delta y = d x + e y + f was fitted
(Table XIII "uncorrected") and the residuals listed in arc (Table XII; mean residual 0".21, p.e. of one
measure +-0".22). The differential refraction (with the different pressure and the ~40 deg F temperature
difference between Oxford and Principe) and the differential aberration of each plate were then removed to
give the **corrected constants a', b', d', e'** (Table XIII): the purely instrumental change of the plate
between Oxford and Principe. The Principe images are larger by a' ~ +31 and e' ~ +37 parts in 1e5 (focus set
about 1.2 mm further out; the lens's focal length grows with temperature, p.320). b' - d' is only the
arbitrary relative orientation of the two plates laid face to face; b' + d' should vanish (mean -1.3).
The adopted scale terms are the separate means a' = +31.2, e' = +37.3 (not one common scale) and
b = d = -0.6 (half of -1.3) (p.322). From the scatter of e', the p.e. of the mean is +-2.1 units, worth
+-0".14 at the limb for all comparisons alike (p.320, p.327).

### Step 2: the theoretical plate constants of each eclipse comparison (pp.321-326)

The plate constants between an eclipse plate and an Oxford plate are predicted, not fitted:
(differential refraction of the eclipse plate: a = -46.5, b = d = +8.2, e = -27.0, aberration zero because
the field is centred on the Sun) - (differential refraction + aberration of the Oxford plate; G1:
-19.1, +0.7, -28.3) + (the check-plate scale: +31.2, -0.6, +37.3). For X-G1 this gives a = +3.8,
b = d = +6.9, e = +38.6 (both sums check exactly). For the other comparisons only the totals are printed;
subtracting the eclipse and scale terms gives the implied Oxford-plate terms H1 (-19.1, -0.7, -28.3),
D1 (-20.2, +1.1, -29.4), I2 (-19.3, -1.5, -28.5), K2 (-18.1, -1.3, -27.4). Their b = d term is antisymmetric in
hour angle as refraction requires: G1 (HA -0.25h) +0.7, H1 (+0.25h) -0.7, D1 (-0.36h) +1.1, I2 (+0.49h) -1.5,
K2 (+0.43h) -1.3, i.e. -b/HA = 2.8-3.1 for all five (b and HA have opposite signs; hour angles recomputed from
the exposure table).

Because the predicted scale is used rather than one fitted from the five eclipse stars, the deflection is
measured over a larger range (p.329: the greatest relative deflection is 1".2 with the scale known, 0".67 with
it eliminated), and the weight is "more than 50 per cent. greater" than from the 12 Sobral stars with unknown
scale (p.329).

### Step 3: the eclipse-plate tables (para 35, pp.322-326)

For each comparison (eclipse plate minus Oxford plate, both in the sense Principe - Oxford):

1. **Delta1** = Delta - (a x + b y) for x, Delta - (d x + e y) for y, with the theoretical constants
   (constant terms c, f left in). Products are rounded to whole units first, so Delta1 agrees with the
   unrounded formula only to 1 unit.
2. **Delta2 x** = Delta1 x - 3/4 Ex: three quarters of the full Einstein x-displacement removed, a compromise
   "between the two hypotheses" (half or full deflection), which divides the possible bias by 4 (p.322).
3. **Orientation theta** by least squares from Delta2 x = c + theta y (p.322). A rotation gives
   Delta x = +theta y, Delta y = -theta x. **Delta3 x** = Delta2 x - theta y, **Delta3 y** = Delta1 y + theta x
   (column "-theta x" subtracted).
4. **Delta4 y** = Delta3 y - Ey, with Ey = 184 alpha_y (the full Einstein y-displacement);
   **Resid** = Delta4 y - mean (the mean rounded to a whole unit). The y-residuals therefore refer to the full
   theoretical deflection 1".75.
5. Normal equations of the conditions **residual = delta_f + alpha_y delta_kappa**, eliminated for
   delta_kappa; deflection = 1".75 (184 + delta_kappa)/184.
6. The x-columns are completed with **Delta4 x** = Delta3 x - 1/4 Ex and residuals from the mean.

X-G1 is printed in full (with the separate a x, b y, d x, e y columns, the coordinates and alpha_y); X-H1 and
W-I2 are "condensed" (no coordinates, no 3/4 Ex column), W-D1 prints coordinates, 3/4 Ex and alpha_y.

Plate **U** (p.326) lacks star 5, without which the orientation is badly determined, so a preliminary
theta = 240 was applied, the full Ex and Ey removed, and a joint least-squares solution of the x and y
residuals made for delta_c, delta_f, delta_theta, delta_kappa: delta_theta = +2, delta_kappa = +121, deflection
2".90 +- 0".87, "practically worthless" and ignored afterwards.

### Step 4: results and error budget (paras 36-37, pp.327-328)

The four determinations are X-G 1".94, X-H 1".44, W-D 1".55, W-I 1".67, mean 1".65. The p.e. of one y
measure on the eclipse plates is +-0".21 (the check plates gave +-0".22); the weight of delta_kappa is ~3 (3.23 for X,
2.87 for W), so +-0".12 on kappa, i.e. +-0".38 at the limb per plate; the four values from two plates give
about +-0".25; adding the common scale error +-0".14 and the correction -0".04 for having removed only 3/4 of
Ex when fitting the orientation gives **1".61 +- 0".30** (probable error). As a standard error this is about
+-0".44 (p.e. = 0.6745 sigma).

## Quotations (error budget, result, conclusions)

p.320: "From the internal discordances of our determination of e' (the most important of these constants) the
probable error of the mean is +-2.1. This, as shown later, will cause a probable error of our final
determination of the deflection, reduced to the limb of the sun, of amount +-0".14, affecting all
determinations systematically. Errors in the other constants have much smaller influence."

p.326-327 (plate U): "The probable error is, however, +-0".87, so that the result is practically worthless.
Further, it is much more likely to be affected by systematic error than the previous results. The large
probable error is partly due to the large residuals which are greater than in the previous measures; in
particular star 3 is unduly faint. If the same accuracy had been obtained, the theoretical weight would have
been half that of plates W and X; but having regard to possible systematic error, probably a quarter weight
would more nearly represent the true value. This determination is ignored in the subsequent discussion."

p.327 (para 36): "It is easy to calculate the effects of any errors in the adopted scale, orientation, &c., on
the final result (deflection at the limb). We give some illustrations. An error in the adopted scale of y of
10 units (in the fifth place of decimals) would lead to an error 0".68 in the result from either plate. Thus
the probable error +-2.1 in the determination of e' gives a probable error +-0".14 in the final result; or, if
we adopted the largest (rather discordant) value found for e' instead of the mean, we should reduce the result
by 0".52. An error of 10 units in the orientation gives an error in the result of 0".45 for plate X, and 0".22
for Plate W. It is therefore of less importance, and further it is not likely to be systematic. Errors in the
measurement of x only affect the result through the orientation. For Plate X, a probable error of +-0".20 in
the x-measures would give an error +-4.0 in the orientation, leading to an error +-0".18 in the result; whereas
an error of the same magnitude in the y measures gives directly an error +-0".35 in the result. For Plate W, the
probable error of +-0".20 in x gives an error +-3.5 in the orientation and +-0".08 in the result, compared with
+-0".38 for similar inaccuracy in y. It is particularly fortunate that the x-measures are so unimportant for
Plate W, because, as already mentioned, the images trailed on that plate. Finally, it will be remembered that
in order not to commit ourselves to the Einstein hypothesis prematurely we neglected the correction 1/4 Ex in
determining the orientation. This will make a difference of 0".029 in the results from Plate W and 0".092 from
Plate X. The effect is that the deduced deflection needs to be decreased, and the mean correction -0".06 should
be applied to the mean result obtained, or rather, to make the adopted deflection for x consistent with the
deduced value from y, the correction needed is -0".04."

p.327-328 (para 37): "The four determinations from the two eclipse plates are X - G 1".94, X - H 1".44,
W - D 1".55, W - I 1".67, giving a mean of 1".65. They evidently agree with EINSTEIN's predicted value 1".75.
The residuals* in the separate comparisons reduced to arc are as follows. They do not appear to show any special
peculiarities. [table] The average y-residual is +-0".22, which gives a probable error for y of +-0".21. It is
satisfactory that this agrees so nearly with the probable error (+-0".22) of the check plates, showing that the
images are of about the same degree of difficulty and therefore presumably comparable. The probable error of x
is +-0".25, but we are not so much concerned with this. The weight of the determination of delta-kappa is about
3 (strictly 3.23 for Plate X and 2.87 for Plate W). The probable error of kappa is therefore +-0".12, which
corresponds to a probable error of +-0".38 in the final values of the deflection. As the four determinations
involve only two eclipse plates and are not wholly independent, and further small accidental errors may arise
through inaccurate determination of the orientation, the probable error of our mean result will be about
+-0".25. There is further the error of +-0".14 affecting all four results equally, arising from the
determination of scale. Taking this into account, and including the small correction -0".04 previously
mentioned, our result may be written 1".61 +- 0".30. It will be seen that the error deduced in this way from
the residuals is considerably larger than at first seemed likely from the accordance of the four results.
Nevertheless the accuracy seems sufficient to give a fairly trustworthy confirmation of EINSTEIN's theory, and
to render the half-deflection at least very improbable."
Footnote *: "The residuals refer to the theoretical deflection 1".75, not the deduced deflections."

p.331 (para 39, General Conclusions): "The Principe observations were generally interfered with by cloud. The
unfavourable circumstances were perhaps partly compensated by the advantage of the extremely uniform
temperature of the island. The deflection obtained was 1".61. The probable error is about +-0".30, so that the
result has much less weight than the preceding. Both of these point to the full deflection 1".75 of EINSTEIN's
generalised relativity theory, the Sobral results definitely, and the Principe results perhaps with some
uncertainty."

## Checks performed

All checks use only the transcribed CSVs (scratch script
`56-light-bending/data/scratch/T4/tools/verify.py`, git-ignored; the formulas are given here so they can be
re-run). "n" = number of values compared. Tolerances are the printed rounding (0.005" for 2-decimal arcsec,
0.5 unit for whole units, 1 unit where two rounded products are subtracted).

### Exposure table (pp.315-316)

| check | n | max deviation | result |
|---|---|---|---|
| Eclipse exposures: start/end from mid-time and exposure give one uniform gap ("assumed uniform") | 15 gaps | all 11.5 s +- 0.5 s (whole-second rounding) | PASS |
| 16 exposures (123 s) + 15 gaps + 1.5 s before and 5 s after = calculated totality 2h13m5s-2h18m7s (302 s) | 1 | 0 s | PASS |
| Approx. Z.D. recomputed for the field centres (Oxford: printed L.S.T.; Principe: printed time as **astronomical** G.M.T.) | 31 | 0.63 deg | PASS |
| Same with the Principe times read as civil G.M.T. | 5 check rows | Z.D. 133-135 deg (field below the horizon) | rejected: the times are astronomical G.M.T. |

### Tables XI-XIII (pp.318-319)

| check | n | max deviation | result |
|---|---|---|---|
| Table XII Mean columns = mean of the printed residuals of the row | 22 | 0.004" | PASS |
| Blank cells of Table XII = blank cells of Table XI | 10 columns | identical | PASS |
| p.e. of one measure 0".22 = Peters 0.8453 sum abs(v) / sqrt(n(n-m)), n = 101, m = 30 | 1 | 0.221" | PASS |
| "mean residual without regard to sign" 0".21 | 1 | recomputed 0.219" | FLAG (rounding level, see D10) |
| Table XII = 0.003 x (Table XI - a x - b y - mean) with the Table XIII constants | 101 | 0.019"; rms 0.0044"; 20 values > 0.005" | PASS (to ~1.5 units) |
| Constants that turn Table XI into Table XII (LS fit of Delta - Resid/0.003) = Table XIII | 20 | 0.42 (w1-b1 e), all others <= 0.2 | PASS except w1-b1 e (D8) |
| Exact least-squares refit of Table XI = Table XIII constants | 20 | 4.7 units; rms 2.3 | FLAG (D9) |
| Table XII = residuals of the exact LS refit | 101 | 0.18"; rms 0.05" | FLAG (D9) |
| b'+d' = b' + d' | 5 | 0 | PASS |
| Means of a', e', b'+d' | 3 | 0.04 | PASS |
| b' - b = d' - d (one refraction/aberration correction for b and d) | 5 | 0 | PASS |
| p.e. of mean e' = +-2.1: Peters 0.8453 sum abs(v) / (n sqrt(n-1)) | 1 | 2.13 | PASS (Gauss's formula would give 1.83) |
| Corrections a'-a, b'-b, e'-e vs a rough refraction (k ~ P/T, tuned to the eclipse values) + annual-aberration model | 15 | ~3 units | consistent (model too rough for a strict test) |

### Theoretical plate constants (pp.321-326)

| check | n | max deviation | result |
|---|---|---|---|
| X-G1: (-46.5, +8.2, -27.0) - (-19.1, +0.7, -28.3) = (-27.4, +7.5, +1.3) | 3 | 0 | PASS |
| X-G1: (-27.4, +7.5, +1.3) + (+31.2, -0.6, +37.3) = (+3.8, +6.9, +38.6) | 3 | 0 | PASS |
| Check-scale b, d = -0.6 = (Table XIII mean b'+d')/2 = -0.65 | 1 | 0.05 | PASS |
| Implied comparison-plate b = d antisymmetric in hour angle (-b/HA = 2.8-3.1 for G1, H1, D1, I2, K2) | 5 | - | PASS |

### Eclipse-plate reductions (pp.323-326)

| check | n | max deviation | result |
|---|---|---|---|
| Coordinates = 14 + 10 x Table I | 12 | 0 | PASS |
| alpha_y, Ey, 3/4 Ex, Ex identical wherever the same star is printed twice | all | 0 | PASS |
| Ey = 184 alpha_y | 6 | 0.48 unit | PASS |
| 3/4 Ex = 0.75 Ex (stars printed with both) | 5 | 0.25 unit | PASS |
| alpha_y and Ex/184 follow 10 (r - r_sun)/r^2 for one Sun position (19.53, 11.53) | 12 | 0.006 | PASS |
| 0.003 Ey, 0.003 Ex vs Table I Principe column | 12 | 0.008" | PASS |
| p.321 relative y-displacements 5-3 = 1".2, 11-3 = 0".8 | 2 | 0.003" | PASS |
| X-G1 printed terms 3.8x, 6.9y, 6.9x, 38.6y | 20 | 0.59 | PASS (2 truncated, D6) |
| X-G1 Delta1 = Delta - printed terms | 10 | 0 | PASS |
| Delta1 x = Delta x - (a x + b y), Delta1 y = Delta y - (d x + e y), all comparisons | 50 | 0.94 unit | PASS |
| Orientation columns = theta y, -theta x | 50 | 0.5 unit, except star 6 in X-G1 (1.6) and W-D1 (0.9) | PASS (D5) |
| Delta2 x = Delta1 x - 3/4 Ex | 20 | 0 | PASS |
| Delta3 x = Delta2 x - theta y | 20 | 0 | PASS |
| Delta4 x = Delta3 x - 1/4 Ex | 20 | 0.33 unit | PASS |
| Delta3 y = Delta1 y + theta x | 20 | 0 | PASS |
| Delta4 y = Delta3 y - Ey (U-K2: Delta1 y + 240 x - Ey; Delta4 x = Delta1 x - 240 y - Ex) | 30 | 1 unit (U-K2 star 3) | PASS except D2 |
| X-G1 mean Delta4 x = -5942 (printed in the text) | 1 | 0 | PASS |
| Resid = Delta4 - round(mean Delta4) | 50 | 0, except W-I2 x star 10 (9 units) | PASS except D1 |
| Normal-equation coefficients 5, sum alpha (2.83, 3.42), sum alpha^2 (4.83, 5.21), eliminated (3.23, 2.87) | 12 | 0.004 | PASS |
| First right-hand side = sum of the printed residuals (-1, +1, -1, +2) | 4 | 0 | PASS |
| Second right-hand side = sum alpha x Resid (64, -105, -62, -24) | 4 | 0.81 unit (computed 64.09, -104.43, -62.81, -23.31) | rounding level (D11) |
| Eliminated right-hand side from the printed normal equations (64, -105, -61, -25) | 4 | 0.57 unit (computed 64.57, -105.57, -61.32, -25.37) | rounding level (D11) |
| delta_kappa = eliminated rhs / coefficient (20, -33, -21, -9) | 4 | 0.49 unit (19.8, -32.5, -21.3, -8.7) | PASS |
| Limb deflection from the printed eliminated equations (1.94, 1.44, 1.55, 1.67) | 4 | 0.003" (1.938, 1.441, 1.548, 1.667) | PASS |
| Limb deflection by exact least squares on the printed residuals | 4 | 0.006" (1.941, 1.440, **1.544**, 1.668) | PASS except W-D1 1.54 vs 1.55 (D11) |
| Steps "20 on 184 = 0".19", "0".31", "0".20", "0".08" = 1".75 delta_kappa/184 | 4 | 0.006" (0.190, 0.314, 0.200, 0.086 from the printed delta_kappa; W-I2's 0".08 follows from the unrounded -8.71: 0.083) | PASS |
| Mean of the four = 1".65; 1".65 - 0".04 = 1".61 | 2 | 0 | PASS |
| U-K2 joint LS (delta_c, delta_theta, delta_f, delta_kappa) from the printed residuals | 4 | delta_theta +1.6 (+2), delta_kappa +115.9 (+121), 2".85 (2".90), p.e. 0".84 Peters / 0".80 Gauss (0".87) | FLAG (D12) |

### p.328 table and error budget (pp.327-328)

| check | n | max deviation | result |
|---|---|---|---|
| p.328 residuals = 0.003 x Resid of pp.323-325 | 40 | 0.005" (7 half-way cases all rounded toward zero) | PASS |
| p.328 Mean = mean of the four | 8 | 0.005" | PASS |
| p.e. of x 0".25: Peters with n = 20, m = 4 on the printed x residuals | 1 | 0.2486" from the p.328 values (0.2504" from the unit residuals of pp.323-325) | PASS |
| average y-residual 0".22 | 1 | recomputed 0.2147" from the unit residuals x 0.003 (0.2135" from the rounded p.328 values) | FLAG (rounding level, D10) |
| p.e. of y 0".21 = 0.8453 x 0".22 x sqrt(20/16) | 1 | 0.208" | PASS (from the printed 0".22) |
| p.e. of kappa 0".21/sqrt(3) = 0".12; x 50/15.78 = 0".38; sqrt(0.38^2/2 + 0.14^2) = 0".30 | 3 | 0.004" | PASS |
| 10 units of e -> 0".68 (X 0.687, W 0.683) | 2 | 0.007" | PASS (one rounded value "for either plate") |
| 10 units of theta -> 0".45 (X), 0".22 (W) | 2 | 0.002" | PASS |
| p.e. 0".20 in x -> theta +-4.0 (X), +-3.5 (W) -> 0".18, 0".08 | 4 | 0.04 unit, 0.003" | PASS |
| p.e. 0".20 in y -> 0".35 (X), 0".38 (W) | 2 | X 0.353, W 0.374 | PASS / W rounding level |
| e' p.e. 2.1 -> 0".14; largest e' 45.0 instead of 37.3 -> 0".52 | 2 | 0.004" | PASS |
| neglected 1/4 Ex in the orientation -> 0".029 (W), 0".092 (X) | 2 | W 0.0291 (PASS); X 0.049 | FLAG for X (D13) |
| mean correction (0.029 + 0.092)/2 = 0".06; self-consistent correction for a deflection of 1.61/1.75 of Einstein's: 0.0605 x (0.92 - 0.75)/0.25 = 0".041 | 2 | 0.001" | PASS (explains the -0".04) |

Totals: about 790 printed values compared in 148 checks. Not checkable from the printed numbers alone: the
plate types, places and dates, the barometer and thermometer readings (used only in the rough refraction model),
the 1.2 mm focus shift and the 0".4 change from the provisional reduction.

## Discrepancies (all printed values kept)

**D1. W-I2 x residual of star 10 (p.325).** Printed Resid +90. Delta4 x - mean = 5219 - 5120 = +99 (the
other four residuals agree exactly with mean 5120; the printed residuals sum to -10, with +99 to -1). The single
change +90 -> +99 restores agreement; Delta4 x = 5219 itself is confirmed by Delta3 x - 1/4 Ex = 5225 - 5.75.
The slip was carried into p.328 (x, I, star 10 = +0".27 = 0.003 x 90; +99 would give +0".30), so it is a
computing slip of the original, not a typesetting error of p.325 alone. It does not affect any deflection
(x residuals are not used). Zoom confirms "+ 90" and "5219".

**D2. U-K2 y, star 3 (p.326).** Printed Delta4 y = 9089, but Delta1 y + 240 x - Ey = 5165 + 4195 - 272 = 9088.
The unrounded Delta1 y is 5985 - 8.9 x 17.48 - 37.7 x 17.60 = 5165.9, so 9089 follows from the unrounded
value; the printed 5165 comes from rounding the two products separately (156 + 664). 1 unit, kept.

**D3. U-K2 y table heading (p.326).** The coordinate column of the y table is headed "x." but holds the
y coordinates (12.40, 18.72, 17.60, 24.99, 27.21, identical to the y columns of X-G1 and W-D1). Heading
misprint; values kept, documented in the CSV header.

**D4. "absence of star 8" (p.326).** The text says "Since Plate U shows some good images it has been examined,
although owing to the absence of star 8 the weight is small." Star 8 (B.D. +20 deg 740, photographic mag 7.0)
appears on no Principe plate (p.320: "nothing fainter than 5'.5 is shown", a magnitude printed with a prime;
the stars listed for the eclipse plates are 3, 4, 5, 6, 10, 11, 12). What plate U lacks is **star 5**: its
stars are 11, 4, 3, 6, 10 (p.320: "Plate U. 6, 10 very bright; 11 faint"), and p.321 says that "the
possibility of a determination of deflection practically depends
on the appearance of star 5 ... the orientation is badly determined in the absence of star 5", which is exactly
the reason given for the low weight of U. "star 8" is therefore a misprint for "star 5". Kept as printed in the
quotation.

**D5. Orientation term of star 6 in X-G1 and W-D1 (pp.323-324).** Printed "+163y" = 4075 and "+91y" = 2275,
whereas 163 x 24.99 = 4073.4 and 91 x 24.99 = 2274.1; both correspond to y = 25.0. Delta3 x follows the printed
terms. 1-2 units in x only, no effect on the result. Kept.

**D6. X-G1 product columns (p.323).** "3.8x" for star 6 printed 75 (3.8 x 19.87 = 75.5) and "38.6y" for star 4
printed 722 (38.6 x 18.72 = 722.6): truncated instead of rounded. Delta1 follows the printed terms. Kept.

**D7. Barometer unit (pp.315-316).** The unit under "Bar." is printed "m." on both pages, but the readings
(29.45-30.48) are inches of mercury and the note speaks of a wave of "about 0.05 in."; "m." is a misprint for
"in.". The column is named barometer_inches. (The D1 reading "29.65" carries a printer's speck between the 9 and
the point.)

**D8. Table XII w1-b1 y residuals vs Table XIII (p.318-319).** The w1-b1 Delta1 y residuals of Table XII are
reproduced to the rounding (max 0.006") with d = 26.0, e = 54.0, but Table XIII prints d = +25.9, e = +53.6
(max difference 0.019", star 6). Either Table XIII's e is misprinted or the residuals were computed with a
slightly different e. With e = 54.0 the corrected e' would be 31.8 instead of 31.4, and the mean e' would still
round to 37.3, so the adopted scale is unaffected. Kept.

**D9. Table XIII constants are not the least-squares solution of Table XI.** The Table XII residuals are,
to their rounding, the residuals of the Table XIII constants (D8 aside) with c, f = the mean: sum v = 0 in every column
but sum v (x - xbar) and sum v (y - ybar) are far from zero (up to 1670 units). Refitting Delta = a x + b y + c
to the printed Table XI gives (Table XIII in brackets):

| pair | a | b | d | e |
|---|---|---|---|---|
| q1-a1 | +34.4 (+32.7) | +101.0 (+101.0) | -90.4 (-87.8) | +57.9 (+58.2) |
| w1-b1 | +24.2 (+26.2) | -18.6 (-16.0) | +28.7 (+25.9) | +58.3 (+53.6) |
| s2-c1 | +29.7 (+31.5) | +189.2 (+192.5) | -174.3 (-173.5) | +64.5 (+64.8) |
| r1-d1 | +24.0 (+28.2) | +162.6 (+165.0) | -145.8 (-146.8) | +69.3 (+69.8) |
| v1-e1 | +21.6 (+21.6) | -77.6 (-76.2) | +73.6 (+70.6) | +62.6 (+61.4) |

Explanations tested and rejected (none reproduces Table XIII for all pairs): coordinates rounded to 1, 0.5 or
0.1 unit; measures rounded to 10 or 100 units; normal equations rounded to 3, 4 or 5 significant figures;
leaving out any one of stars 1, 5, 6, 7, 8, 11, 13, or both faint stars 5 and 8; weights by the brightness classes
of p.319; every subset of 4 or more stars per column (only one column, r1-d1 a/b, is matched by one subset, by
chance); every grouping of the 11 stars into 3 groups (method of averages) for q1-a1; coordinates referred to the
plate centre (14, 14) without fitting c; a joint fit with a common offset per star; a partly converged
iterative (Gauss-Seidel) solution of the normal equations, starting with a or with b (1-10 iterations; it only
approaches the least-squares values). Unpublished extra stars in the fit (check stars 3 and 9 are missing from
the table) are implausible: sum v over the printed stars is ~0 (|sum| <= 10 units) in every column while
sum v (x - xbar) reaches 1670 units, which would need hidden residuals of about +-100 units (0".3) cancelling
each other in every column. Misreading of Table XI is excluded, because Table XII is reproduced from the
printed Table XI. The "a few units" found by the earlier
reader is therefore real and comes from how the 1919 computers obtained a, b, d, e (an approximate or
non-least-squares solution, or constants from an earlier state of the measures), which the printed numbers
cannot reveal. Consequence: with the least-squares constants the mean corrections would be a' = 30.0,
e' = 38.2, (b'+d')/2 = -1.25 instead of 31.2, 37.3, -0.65; re-reducing the four eclipse comparisons with these
(and a least-squares theta) changes the mean deflection by only -0".003.

**D10. Mean residuals quoted in the text.** p.318: "mean residual without regard to sign is +-0".21", but the
101 printed residuals give 0.219"; the printed p.e. 0".22 follows from 0.219" (Peters, m = 30), not from 0.21"
(-> 0.212"). p.328: "The average y-residual is +-0".22", but the 20 y residuals give 0.2147" (pp.323-325
units x 0.003; 0.2135" from the rounded p.328 entries), i.e. 0".21; the printed p.e. 0".21 follows from 0.22
with Peters' formula for n = 20, m = 4. Rounding-level (0.005-0.009"), kept.

**D11. Normal-equation right-hand sides.** The second right-hand sides 64, -105, -62, -24 differ from
sum alpha_y x Resid of the printed residuals (64.09, -104.43, -62.81, -23.31) by up to 0.8 unit, and the
eliminated right-hand sides 64, -105 (X) are 64.57, -105.57 from the printed normal equations. The
delta_kappa values (20, -33, -21, -9) and the deflections follow from the printed equations; an exact
least-squares solution of the printed residuals gives 1.941, 1.440, 1.544, 1.668 (W-D1: 1".54 vs printed 1".55,
because the printed rhs -62 should be -62.8). 0.006" at most; kept.

**D12. Plate U solution (p.326).** The paper gives delta_theta = +2, delta_kappa = +121, 2".90 +- 0".87, but a
least-squares solution of the printed U-K2 residuals for (delta_c, delta_theta, delta_f, delta_kappa) gives
delta_theta = +1.6, delta_kappa = +115.9, 2".85, p.e. +-0".84 (Peters) or +-0".80 (Gauss). Variants tried (x
weights 0.25-3, alpha rounded to 2 decimals, no alpha_x term in the x equations, theta fixed at +2, y equations
alone) give delta_kappa between 93 and 129 (the opposite sign of theta in y gives 28-59), but never 121
together with delta_theta = +2; the 1-unit slip D2 changes delta_kappa by only 0.4. No normal equations are printed, so the
paper's exact procedure cannot be reconstructed. The difference (0".05) is far inside the quoted p.e. and the
plate was discarded by the authors. Kept.

**D13. Neglected 1/4 Ex for plate X (p.327).** The text says the neglect of 1/4 Ex in fitting the orientation
changes the result by 0".029 for W and 0".092 for X. Recomputing (the slope of 1/4 Ex against y shifts theta by
1.30 units for W and 1.09 for X; 10 units of theta = 0".22 and 0".45): W 0.029" (reproduced exactly) but X
0.049", not 0.092". The same values follow from a full re-reduction (theta fitted after removing all of Ex:
X-G1 and X-H1 change by -0.049", W-D1 and W-I2 by -0.029"). The paper's -0".06 mean and the adopted -0".04
(0.0605 x (1.61/1.75 - 0.75)/0.25 = 0.041) both use 0.092; with 0.049 the mean correction would be -0".039 and
the self-consistent one about -0".03, final result ~1".62. Kept as printed.

**D14. Printed orientations are not the least-squares solutions of the printed Delta2 x.** p.322 says theta is
found "by least squares" from Delta2 x = c + theta y. The printed Delta2 x and y give theta = +165.6 (printed
+163), +14.3 (+10), +93.3 (+91), -28.8 (-30); all four printed values are 1.2-4.3 units lower. The printed
Delta3 x columns keep a visible trend with y. Not explained by leaving out one star, rounding y to whole units
(163.9, 13.7, 93.8, -29.3) or a joint x + y solution (163.8, 17.1, 93.2, -31.7); removing the full Ex before the
fit reproduces W-I2's -30 but gives 164.5, 13.2, 92.0 for the others. A possible cause, not provable from the
printed numbers: theta carried over from
the provisional reduction (p.321 footnote: preliminary results with a provisional scale, then recalculated
refraction and aberration), since b and theta both multiply y in the x-equation. Effect: with the least-squares
theta (keeping 3/4 Ex as the paper did) the four deflections become 2.05, 1.63, 1.59, 1.69 (mean 1".74 instead of
1".65); removing all of Ex before fitting theta (which makes the -0".04 correction unnecessary) gives 2.00,
1.58, 1.56, 1.66, mean 1".70. Printed values kept.

**D15. "0".67" (p.329).** "if the scale is unknown and must be eliminated, this [greatest relative deflection]
is reduced to 0".67". For the five Principe stars of X or W, removing f + e y (or f + d x + e y) from Ey leaves a
range of 0.74-0.78". The definition behind 0".67 is not stated (it may refer to the 12 Sobral stars); not
resolved.

## Time scale of the exposure list (pp.315-316)

The column is headed "Loc. Sid. T." for the Oxford plates; a rule and the sub-heading "G.M.T." precede the
Principe plates. The paper does not say astronomical or civil, but the numbers settle it: read as **astronomical
G.M.T.** (day beginning at noon, the 1919 astronomical convention), all 31 recomputed zenith distances agree with
the printed "Approx. Z.D." within 0.6 deg; read as civil time, the Arcturus field would be 133-135 deg from the
zenith for q1-w1. The check plates q1, r1, s2 were therefore taken around 00h25m-00h38m U.T. on 23 May (civil),
v1, w1 around 00h25m U.T. on 26 May, and the eclipse plates at 14h13m-14h18m U.T. on 29 May. The eclipse times
are themselves "deduced from the calculated time of totality" (2h. 13m. 5s. to 2h. 18m. 7s. G.M.T., p.314), not
observed.

## Unresolved

* How the Table XIII uncorrected constants were obtained (D9); whether w1-b1 e is 53.6 or 54.0 (D8).
* Why the printed orientations differ from least squares by 1-4 units (D14), the largest effect found (+0".10
  on the mean if the least-squares orientations are used).
* The U-K2 solution (D12) and the plate-X value 0".092 (D13) cannot be reproduced.
* The 0".67 of p.329 (D15).

## Audit (second reading)

Independent blind double-keying of all seven CSVs of this group (second transcriber, 2026-10-20). The blind files
(same names and columns), the reading log and the audit scripts are in the git-ignored scratch folder
`56-light-bending/data/scratch/T4-audit/` (`blind/`, `blind/readings_log.txt`, `tools/compare.py`,
`tools/check_*.py`, `tools/rereduce.py`; full output of the rerun checks in `checks_on_transcribed.txt`).

**Protocol.** Until the blind reading was complete, only the '#' header lines, the CSV header rows and the key
columns (field/ref, star, comparison/coord/star, context/quantity) of the first transcriber's files were looked at.
Every value was read from the 5120 x 6600 JPEG-2000 scans (table cells on 1.0-2.0x crops, running-text numbers on
0.42-0.5x crops of the full-resolution page), then compared cell by cell by script (numeric equality at the printed
precision, empty cells included).

**Comparison.**

| file | cells compared | printed numbers among them | disagreements |
|---|---|---|---|
| `dyson1920_principe_exposures.csv` | 372 | 217 (+31 page numbers) | 62, all representation only (see below) |
| `dyson1920_table11_check_measures.csv` | 132 | 123 | 0 |
| `dyson1920_table12_check_residuals.csv` | 132 | 123 | 0 |
| `dyson1920_table13_check_constants.csv` | 54 | 48 | 0 |
| `dyson1920_principe_plates.csv` | 800 | 410 (+50 page numbers, 50 printed headings) | 0 |
| `dyson1920_principe_plate_constants.csv` | 264 (value, page) | 132 | 3 sign conventions + 1 form (see below) |
| `dyson1920_p328_principe_residuals.csv` | 60 | 48 | 0 |
| total | 1814 | 1101 | no printed digit or sign read differently |

Rows and empty (dash) cells were identical in every file. The 65 flagged differences are all representation:
* Exposures: the first transcriber writes the date as "1919 January 16" and the time scale as the printed heading
  "Loc. Sid. T." / "G.M.T."; the blind file had "January 16" and "LST"/"GMT". Zooms confirm that "1919." is printed
  once at the head of the Date column and that the "G.M.T." sub-heading stands below a rule above the Principe rows,
  on both pages. The first transcriber's form is the more faithful one and was kept.
* Plate constants, `limb_deflection_change`: the paper prints these numbers unsigned ("an increase of 0".19",
  "1".75 - 0".31 = 1".44", "- 0".20", "- 0".08"). The blind file had 0.19, 0.31, 0.20, 0.08, the first transcriber
  +0.19, -0.31, -0.20, -0.08 (sign taken from the wording). The signed values were kept, and the convention is now
  stated in the CSV header.

**Corrections to transcribed data cells: none.**

**Cells re-checked on tight full-resolution zooms** (all confirmed as transcribed): W-I2 x Resid star 10 "+ 90" and
Delta4 x "5219" (D1); U-K2 y Delta4 y star 3 "9089" (D2); p.326 "absence of star 8" (D4); X-G1 x row of star 6
(3.8x = 75, +163y = 4075; D5, D6); X-G1 y row of star 4 (38.6y = 722; D6); W-D1 x row of star 6 (+91y = 2275; D5);
the Table XIII rows (w1-b1 e "+53.6", r1-d1 e "+69.8" with a speck after the 9, the Mean row with dashes under b', d');
Table XI w1-b1 Delta y star 6 "5140" and r1-d1 Delta y star 13 "4120"; Table XII w1-b1 Delta1 y star 6 "+0.27"
and r1-d1 Delta1 y star 13 "-0.38" (the two largest deviations of the Table XI -> XII check); "1919." and "m." at the
head of the exposure tables; the "x." heading of the U-K2 y table (D3).

**Checks rerun with independent scripts on the final files.** Exposures: 15 gaps 11-12 s, 1.5 s + 123 s + 15 gaps
+ 5 s = 302 s; Approx. Z.D. reproduced to 0.63 deg with L.S.T. (Oxford) and astronomical G.M.T. (Principe), civil
G.M.T. gives 132.6-135.3 deg for the check field. Tables XI-XIII: Mean columns 0.004"; XI -> XII with the Table XIII
constants rms 0.0044", max 0.019", 20 of 101 above 0.005"; implied constants within 0.20 except w1-b1 e (54.02);
exact LS refit differs by up to 4.74 (rms 2.26), residuals by up to 0.176" (rms 0.052) - the D9 table is reproduced
digit for digit; b'+d' and b'-b = d'-d exact; means +31.24, +37.26, -1.26; p.e. of mean e' 2.13 (Peters);
mean |v| 0.219", Peters (n = 101, m = 30) 0.221". Reductions (320 identities): all exact or within rounding except
D1, D2, D5, D6; Ey = 184 alpha_y within 0.48; one Sun position (19.534, 11.533) fits alpha_y and Ex/184 to 0.0056;
coordinates = 14 + 10 x Table I exactly. Normal equations, delta-kappa, deflections, the U-K2 joint solution
(+1.62, +115.95, 2".853, p.e. 0".84 Peters / 0".80 Gauss) and every number of the para 36-37 error budget reproduce
the values given in the tables above. p.328: all 40 entries = 0.003 x Resid within 0.005" (7 half-way cases, all
rounded toward zero); means within 0.005". So every check passes or fails only at a documented misprint (D1, D2, D5,
D6) or a documented analysis question (D8-D15).

**New observations from the audit.**
* Delta1 is reproduced *exactly* (not just to 1 unit) for all 50 rows as Delta minus the two products rounded
  separately, provided the products 3.8x = 75 (star 6) and 38.6y = 722 (star 4), printed truncated in X-G1 (D6),
  are also used in X-H1, which has the same a = 3.8 and e = 38.6: X-H1 x star 6 7429 - 75 - 207 = 7147 and X-H1 y
  star 4 1881 - 144 - 722 = 1015, both as printed (with rounded products they would be 7146 and 1014). So D6 carries
  over into X-H1, and every Delta and Delta1 digit of the five comparisons is confirmed at the units level.
* w1-b1 y (D8): with e = 53.6 the deviations of the printed residuals trend with y (+3.8, +3.2, -6.3 units for
  stars 11, 12, 6 at y = 0.93, 6.23, 27.33); with e = 54.0 all are within 1.9 units. This is a constant slightly
  different from the printed one, not one misread or misprinted cell (5140 and +0.27 confirmed on zooms).
* r1-d1 y star 13: printed -0".38 against -0".392 recomputed (+4.1 units, the second-largest deviation); the other
  stars of that column deviate by -2.4 to +0.4 units and the printed column sums to -0".02, so this is the
  computers' rounding scatter, not a misprint.
* Plate U: agree with D4 - the U tables list stars 11, 4, 3, 6, 10, so the missing star is 5; "star 8" on p.326 is
  a misprint for star 5.
* D9 consequence reproduced: re-reducing with the least-squares-based scale constants (a -1.2, b = d -0.6, e +0.9)
  and least-squares orientations changes the mean deflection by -0".0028.

**Fixes made to the summaries in this audit.**
* This file: "-b/HA = 2.8-3.1" (b and the hour angle have opposite signs; the text said "b/HA = 2.8-3.1"), twice;
  the eliminated right-hand side of W-I2 from the printed normal equations is -25.37 (was given as -24.68, which is
  the value from the residuals themselves); the exact least-squares W-I2 deflection is 1.668 (was 1.667), twice;
  the p.e. of x by Peters' formula is 0.2486" from the p.328 values, 0.2504" from the unit residuals (was 0.2496");
  the average y-residual 0.2147" is from the unit residuals, 0.2135" from the rounded p.328 values (made explicit,
  twice); Step 1 said the check field was taken at Oxford and Principe "with the same lens, 8-inch stop and
  coelostat", but the Oxford plates were taken direct, without the coelostat (p.317), now corrected; para 39 begins
  on p.330 (its Principe sentences are on p.331).
* `dyson1920_principe_plate_constants.csv` header: "paras 33-37" -> "paras 32-38" (the rows run from p.318,
  para 32, to p.329, para 38); added the sign convention of `limb_deflection_change`. No data cell changed.
* All quotations in this file and in the CSV headers and note columns were checked word by word against the
  scans of pp.312-331 (p.314 totality and plate list, p.315 notes, p.316, p.317, p.318, p.319, p.320, p.321
  including both footnotes, p.322, pp.323-326 text, pp.326-327 plate U, p.327 para 36 and 37, p.328, p.329, p.331)
  and found faithful apart from the points above.
