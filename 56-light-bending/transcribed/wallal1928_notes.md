# Wallal 1922 (Lick): transcription notes, checks and discrepancies

This file covers the transcription of Campbell & Trumpler 1928, *Lick Obs. Bull.* 13, 130–160 (five-foot cameras), and of Table II of Campbell & Trumpler 1923, *Lick Obs. Bull.* 11, 41–54 (15-foot cameras). All page numbers are printed bulletin pages. The PDF page is the printed page − 129 for 1928 and − 40 for 1923.

**In short.** About 10,400 printed numbers were transcribed into 15 CSV files. Every table was checked against the paper's own arithmetic and against the other tables that are derived from it. No transcription error remains that the checks can detect. A number of genuine inconsistencies in the printed bulletin are documented in §3. Every one of them is **kept as printed**, with what the arithmetic implies.

## 1. How the transcription was made

- **Images.** I worked from the 600-dpi, 1-bit JBIG2 page images embedded in the ADS PDFs (extracted losslessly with `pdfimages`). That is three times the resolution of the 200-dpi PNGs. Skewed pages were rotated by 0.1°–0.7° so that rows run horizontally. Table rows were located from the ink profile of the star-number column. Each block of 11–20 rows was read at about 300 dpi equivalent. Every doubtful cell was re-read at full 600-dpi resolution.
- **Second opinion.** I ran tesseract on every single cell (digits only): about 8,200 cells in Tables 4, 5 and 6 and 1,100 in the 1923 Table II. Every disagreement in digits was either decided by the arithmetic checks below or re-viewed at full resolution. In every case the OCR was wrong: a dropped or truncated digit, ″ read as 7, or a 3 read as 8.
- **Scripts.** The checks are scripts in the git-ignored scratch folder `data/scratch/T5/scripts/` (`verify_all.py`, `verify_t23.py`, `make_csvs.py`, `make_results.py`). `verify_all.py` reads the final CSVs exactly as `pandas.read_csv(path, comment='#', skipinitialspace=True)`.
- **CSV conventions.**
  - Numbers are copied at printed precision.
  - The leading zero is restored: `.55` is written 0.55 and `+″70` is written +0.70.
  - Explicit `+` signs are kept as printed.
  - An empty cell means nothing is printed.
  - A value printed in parentheses (rejected by the authors) is kept, with flag `rej_<pair> = 1`. Such values are excluded from every mean, as in the paper.

## 2. Checks performed (all on the final CSV files)

| Check | Values | Result |
|---|---|---|
| **Table 6** (pp.147–149): set of stars | 145 rows | Table 4 numbers 1–147 without 66 and 90 (dropped, p.132); distances non-decreasing as printed ✓ |
| Table 6: `wt_total` = Σ plate weights | 145 | exact in all rows ✓ |
| Table 6: `Dr_mean` = Σ w·D_r / Σ w | 145 | max deviation **0.0060″**. 11 rows lie between 0.005 and 0.0060 (10 of them above 0.0051), which is expected from inputs rounded to 0.01 ✓ |
| Table 6: plate D_r vs D_x·x/r + D_y·y/r from Table 4 | 828 | 819 agree within the rounding tolerance (≈0.010–0.015″, depending on x, y and the tangential part). **9 disagree by 0.014–0.037″** (§3.1). Every per-plate weight in Table 6 equals Table 4 ✓. Rejected Table 4 measures have no D_r in Table 6 ✓ |
| Table 6: `distance_deg` vs gnomonic distance atan(r·27084″) from Table 4 x, y | 145 | mean +0.0010°, rms 0.0052°, max 0.0143° ✓. This confirms the unit is degrees and the coordinate unit is 27084″ |
| **Table 4** (pp.139–144): `wt_mean` = Σ non-rejected weights | 145 | exact ✓; also equals Table 6 `wt_total` for every star ✓ |
| Table 4: `Dx_mean`, `Dy_mean` = weighted means | 290 | 287 within 0.0101″. Of these, 27 lie between 0.005 and 0.0101; the authors evidently averaged unrounded plate values. **3 are off by 0.016–0.017″** (stars 67, 102, 123; §3.2) |
| Table 4: radial projection of the mean (Dx, Dy) vs Table 6 `Dr_mean` | 145 | within rounding except stars 102 and 123 (0.012″; same cause as above) |
| Table 4: RA/Dec 1900 vs x, y (precessed to 1922.0 FK4, 4-parameter gnomonic fit) | 145 | rms 0.0013 / 0.0009 units in x / y, max 0.0070 / 0.0036. The fitted centre is RA 11h50.56m, Dec +1°01′, the Sun's position (1923 paper: 11h50.6m, +1°1′) ✓. **Only star 129 fails** (−0.125 unit in y; §3.3) |
| **Table 5** (pp.143–145): weight sums / vx, vy means | 75 / 150 | exact / max 0.0058″ (1 value > 0.005) ✓ |
| **Table 7** (pp.149–150) recomputed from Table 6: weights | 72 | exact ✓; the group sizes 5, 6, 9, 25, 20, 19, 19, 20, 22 are fixed by this check |
| Table 7: radial displacement | 72 | max 0.0063″ ✓ |
| Table 7: group distance | 72 | weighted means reproduce only to ≤ 0.031° (39 values beyond 0.005). The rule the authors used is unknown, and no single value stands out (§3.6) |
| Table 7, all plates: stars, weight, D_r, O−T | 9 groups | exact; exact; ≤ 0.0043″ from the per-plate values (≤ 0.006″ from the star means); exact ✓. "Einstein −″.009r" at the printed distance is within 0.006″. The p.e. equals 0″.33/√weight to rounding ✓ |
| **Table 8** (p.150) from the Table 5 means sorted by distance | 7 groups | the weights partition the 75 stars exactly into groups of 12, 10, 10, 11, 10, 10, 12 ✓; D_r within 0.0057″ ✓; distances only to ≤ 0.05° (§3.6) |
| **Table 9** (p.151) | 6 plates | weighted means give 1.7975 (printed 1″.80), **1.8438 (printed 1″.85)** and 1.8200 (1″.82); α 0.0095 (+.010). p.e. 0.33/√5.14 = 0.146 (±″.15). β-scatter p.e. of unit weight = 0.373 (±″.37) ✓ |
| Table 9: refit −αr + βR/r = D_r to the Table 6 columns (R = 0°.2661) | 6 | weights within 0.007, α within 0.0012, β within 0.074 (N18–N11: 2.074 vs 2.00), 0.042 (M16–M11), ≤ 0.032 otherwise. All are below β's p.e. (≈ 0″.35–0″.39). Fits to the Table 7 group means do worse, so the paper fitted star by star |
| **Table 10** (pp.152–153) vs 1923 Table II and 1928 Table 6 | 103 rows | weight = w₁₅ + w₅/3 within 0.0033 ✓. D_r = (w₁₅(D₁ + 0.060 r) + (w₅/3)·Dr₅)/w within 0.0091″ ✓. Distance within 0.0064° ✓. Einstein column within 0.0055″ ✓. All 65 five-foot stars ≤ 3°.45 and all 92 fifteen-foot stars present ✓ |
| **Table 11** (p.153) | 11 groups | groups 2–7 from consecutive Table 10 rows (17, 14, 12, 14, 18, 18 stars) and groups 8–11 = Table 7 groups 6–9 with weight/3: all ✓. **Group 1 weight 8.92 vs 13.25** (§3.4). Einstein column ≤ 0.0049″ ✓. Courvoisier column ≤ 0.0063″ (one value just over rounding). O−E and O−C exact ✓. The p.e. equals 0″.20/√weight ✓ |
| **Table 12** (p.155) | 15 stars | O − T = obs − theory exactly (30/30) ✓. Theory (1″.75R/r − 0″.009r along the radius) within 0.006″ ✓. Observed from 1923 (Ax, Ay with the +0″.060/° radial correction) and 1928 Table 4 means within 0.0074″ ✓ |
| **Areas** (pp.146–147) | 8 | with a central square \|x\|, \|y\| < 0.5 unit, the star counts 16, 15, 18, 17 / 18, 21, 24, 16 are reproduced exactly. 7 of 8 areas reproduce to ≤ 0.006″; **area IV Dx** does not (§3.5). corrected − observed = (+.04, +.07) for all outer areas and (−.04, −.06) for all central areas, i.e. uniform ✓ |
| **Table 2 / Table 3** (pp.134–137): det → corr → adopt → corr chain | 180 relations | as printed: 23 relations off by > 0.015″. 20 involve g or h′ (no g_a or h′_a row is printed), 2 involve i (the i_a sign) and 1 is M17–M9 e′ (§3.7). Modelling g_a = i_a and h′_a = k′_a and setting the i_a of "M7 and N9 Check F." to +.01: 179 of 180 within 0.0101″, 1 off by 0.02 (M17–M9 e′) |
| Table 3 corr columns vs the p.138 check−eclipse scale table | 8 | +0.75, +1.24, +0.99, +1.41, +1.14, +1.31, +1.16, +1.21 ✓; mean +1.15 ✓ |
| Text claims | | "more than half" of the Table 7 O−T below p.e. (5 of 9) ✓; Table 8 "four of the seven … smaller" ✓; Table 11 "6 being smaller, 5 larger" ✓; "+″.035" mean of four residuals ✓; Courvoisier "5–6 times the probable error" (−.17/.03, +.20/.04) ✓; "5 per cent" (0.09/1.75) ✓; (2·1.72 + 1.82)/3 = 1.753 ✓; 26 fifteen-foot reference stars also measured with the 5-foot cameras ✓ (§3.8) |
| **Table 1** (p.132): simultaneous M and N exposures must agree | 10 line pairs | 3 differences (§3.9) |
| **1923 Table II** (pp.49–50): dist = √(X² + Y²) | 92 | ≤ 0.0098 ✓ |
| 1923 Table II: D₁ = radial projection of (Ax, Ay) | 92 | ≤ 0.0175″, all within rounding tolerance ✓ |
| 1923 Table II: theory = 1″.75·0°.2661/(1.2522 d) − 0″.079 d | 92 | ≤ 0.0126″. 12 values lie between 0.005 and 0.013; these are the nearest stars, where rounding d to 0.01 explains the difference ✓ |
| 1923 Table III recomputed from Table II | 9 groups | weights exact for 8 of 9. D₁, D₂ and theory means within rounding except group 3 D₂. Distances within 0.008° except group 7, which is off by 0.055° (§3.10) |
| 1923 Table IV: (1.80 + 1.48 + 1.85 + 0.9·1.76)/3.9 | 1 | 1.7215 → 1.72 ✓. Plate means are the means of the two observers ✓ |

## 3. Discrepancies (printed numbers kept; nothing was changed to force a check)

### 3.1 Table 6 per-plate D_r vs Table 4 (9 cells)

The weights agree in all nine cells. Every printed value was re-read at 600 dpi in both tables, and the per-cell OCR agrees with my readings. For each case the table lists the single changes that would restore agreement. The last column says what the tables' own mean columns support. In the Table 6 row, the printed `Dr_mean` is recomputed (i) with the printed cell and (ii) with the value implied by Table 4.

| Star | Pair | Table 4 Dx, Dy (x, y) | D_r implied | Table 6 printed | Δ | Single changes that restore agreement | Support from mean columns |
|---|---|---|---|---|---|---|---|
| 57 | N20–N9 | +.10, −.49 (−.137, −.400) | +0.431 | +.41 | 0.021 | T6 → +.42…+.44; or T4 Dx +.13…+.16 or Dy −.46…−.48 | T6 mean fits +.43 slightly better (dev 0.000 vs 0.003) |
| 61 | M17–M9 | −1.27, −.57 (−.112, +.160) | +0.261 | +.28 | 0.019 | T6 → +.25…+.27; T4 Dy −.53…−.56; x −.113…−.115 | T6 mean fits +.26 slightly better (0.001 vs 0.004) |
| 62 | M17–M9 | −.04, −1.11 (−.076, −.045) | +0.600 | +.58 | 0.020 | T6 → +.59; T4 Dx −.01…−.03, Dy −1.10; y −.041…−.044 | T6 mean supports printed +.58 (0.001 vs 0.010) |
| 71 | N18–N11 | +.85, +.16 (+.015, −.152) | −0.076 | −.09 | 0.014 | T6 → −.07/−.08; T4 Dx +.75…+.84, Dy +.17/+.18 | neutral |
| 73 | M17–M9 | −1.23, −1.57 (+.017, +.098) | −1.757 | −1.74 | 0.017 | T6 → −1.75; T4 Dy −1.54…−1.56 | T6 mean supports printed (0.003 vs 0.009) |
| 80 | N20–N9 | +.25, −.47 (+.081, +.243) | −0.367 | −.35 | 0.017 | T6 → −.36/−.37; T4 Dx +.27, Dy −.44…−.46 | T6 mean fits −.37 slightly better |
| 96 | N20–N9 | −1.83, −.58 (+.230, −.126) | −1.326 | −1.31 | 0.016 | T6 → −1.32/−1.33; T4 Dx −1.80…−1.82, Dy −.59 | T6 mean supports printed (0.002 vs 0.005) |
| 101 | N19–N13 | −.68, +.18 (+.246, −.067) | −0.703 | −.72 | 0.017 | T6 → −.70/−.71; T4 Dx −.69 | T6 mean fits −.70 slightly better |
| 128 | N20–N9 | −.10, −.33 (+.565, +.909) | −0.333 | −.37 | 0.037 | T4 Dy −.36…−.38 | **both** mean columns support their own table: T6 mean −.02 needs −.37 (−.33 → dev 0.012); T4 Dy mean +.07 needs −.33 (−.37 gives +.063). The slip happened when D_r was formed from D_x, D_y, not in typesetting |

**Resolution.** All 18 values are kept as printed. The effect on any fit is negligible (≤ 0.037″ in single plate values with p.e. ≈ 0.3–0.6″).

### 3.2 Table 4 mean columns that disagree with the plate values (3 stars)

- **Star 67, Dx_mean.** Printed −.20; the weighted mean of the printed plate Dx is −0.216. The radial checks cannot locate the cause because x ≈ −0.02, so D_r hardly depends on Dx. Many single changes would restore agreement, for example M16–M11 Dx −.18 → −.13…−.10 or N20–N9 Dx −.76 → −.66…−.73. **Unresolved, kept.**
- **Star 102, Dy_mean.** Printed +.15; computed +0.133. Dx_mean is printed +.84 against a computed +0.847, which is a rounding edge. Table 6 Dr_mean is +.85, consistent with its own plate values, while the radial component of the Table 4 means is +0.84. No single change fixes all relations. **Unresolved, kept.**
- **Star 123, Dx_mean.** Printed +.01; computed −0.007 (Δ 0.017). Table 6 Dr_mean is +.02, which is consistent with its plate values but not with +.01 as Dx_mean (radial +0.033). Changing only the sign (−.01), or reading .00, restores **every** check. This is probably a sign misprint in Table 4. **Kept as printed (+.01).**

### 3.3 Star 129: catalogue position inconsistent with x, y

Table 4 prints star 129 as **B.D. +3°2612, α₁₉₀₀ 12h6m.7, δ₁₉₀₀ +3°41′**, with no magnitude (the only blank magnitude), and **x = +.570, y = +.465**. The printed BD number and declination agree with each other: BD +3°2612 has 1855 coordinates 12h04m20.6s, +3°56.1′, which precess to +3°41′04″ in 1900. The printed y, however, implies δ ≈ +4°36′. All other 144 stars fit their x, y to 0.0036 unit in y (1.6′); star 129 is off by −0.125 unit (−57′).

The x, y position matches **BD +4°2583** (1855: 12h04m16.4s, +4°51.5′; VizieR I/122) to (+0.0004, −0.0018) unit. The reduction used x, y: Table 6 gives 5°.52 for this star, which equals the gnomonic distance from x, y. The star actually measured was therefore very probably BD +4°2583, and the identification line of Table 4 contains an error. **Kept as printed.** For any use of the catalogue columns, trust x, y. The displacement data are unaffected.

### 3.4 Table 11, group 1 weight

The printed weight is **8.92**. The ten Table 10 stars of the group (0°.55–0°.93, n = 10 as printed) have weights summing to **13.25**. With those weights the group mean distance (0.674 vs 0°.675) and displacement (+0.632 vs +″.63) are reproduced, so the group's membership and weighting are as in Table 10. The printed p.e. ±″.07 equals 0″.20/√8.92 = 0.067, whereas √13.25 would give ±0.055. The 8.92 was therefore used for the p.e. No single digit change explains it, and no rule found here produces 8.92. **Unresolved, kept as printed.**

### 3.5 Area IV (p.146)

Dx is printed **.00**. The Table 4 means of the 17 area-IV stars give −0.031″; the count and Dy (−.08) match. The corrected value on p.147 (+.04 = .00 + 0.04, the same shift as the other outer areas) shows the authors worked with .00. The other seven areas reproduce. **Kept; cause unknown.** The magnitude-equation difference "bright-faint Δx = −″.08, Δy = −″.02" (p.146) could not be reproduced either. Weighted means of the Table 4 stars brighter vs fainter than 9m.0 give −0.071 and −0.054. The method is not described.

### 3.6 Group distances (Tables 7 and 8)

The weights and radial means of every group reproduce, but the printed mean distances do not follow exactly from the weighted means of the Table 6 (or Table 5) distances:

- per-plate Table 7: up to 0.031°;
- per-objective Table 7: up to 0.014°;
- all-plates Table 7: up to 0.008°;
- Table 8: up to 0.05°.

Unweighted, harmonic and differently weighted means, and an offset of the origin, all fit worse. No single entry stands out, so this is a property of the authors' arithmetic and not a misprint. **Kept.** Related small effects at the rounding level:

- Table 7 all-plates group 4 prints "Einstein −″.009r" +.20, while the formula at the printed 2°.18 gives +0.194; at 2°.17 it would give +0.195.
- Table 9 "Mean of obj. N" prints 1″.85; the weighted mean of the printed plate values is 1.8438.
- Table 11 Courvoisier column at 6°.56 prints −.07; the formula gives −0.064.

### 3.7 Tables 2 and 3: aberration terms

- **Missing rows.** The printed aberration formula contains g_a x³ (in A_x) and h′_a x²y (in A_y), but Table 2 has no g_a or h′_a row. With those terms set to zero, 20 Table 3 relations for g and h′ fail by 0.02–0.03″. In every printed column but one, the other aberration terms obey b_a = c′_a, c_a = −b′_a, f_a = −d_a and i_a = k′_a. Assuming the analogous g_a = i_a and h′_a = k′_a fixes 18 of the 20. The remaining two (g of M15–M7 and N20–N9) are fixed, together with the two i relations, by the i_a sign below. The unprinted rows are noted in the CSV header; nothing was added to the CSV.
- **One sign.** i_a of "M7 and N9, Check F." is printed **−.01**, while its k′_a is +.01; all other columns have i_a = k′_a. The Table 3 relations for g and i of M15–M7 and N20–N9 (the pairs with the M7/N9 check field) reproduce exactly with **+.01** and fail by 0.02″ with −.01. This is probably a sign misprint in Table 2. **Kept as printed (−.01).**
- **M17–M9 e′.** "Check F. det." is printed −2.48; deducting refraction, aberration and the inclination term (d_i = e′_i = +.08, confirmed by the d row) gives −2.25, while "corr." is printed −2.23. "Eclipse F. adopt." (−2.67) agrees with the printed −2.23. So either the det value is a misprint for −2.46 or the authors made a slip of 0.02. **Kept.**

### 3.8 Text p.152: "+0″.15 … 2°.52"

The "26" stars (15-foot reference stars beyond 2° that were also measured with the 5-foot cameras) are reproduced exactly. For them, Table 6 gives a mean 5-foot D_r of +0.188″ unweighted (+0.174 weighted with the 15-foot weights). The 5-foot minus 15-foot difference, weighted, is +0.153″. The mean distance is 2°.54. The printed +0″.15, which yields the +″.060/° scale correction, therefore corresponds to the inter-instrument difference. That is the quantity that matters, since the 15-foot reference stars average +0.02″ on their own system. Table 10 consistently uses +″.060. No discrepancy is claimed.

### 3.9 Table 1, simultaneous exposures (p.132)

M and N plates of one pair were exposed together, so their lines must agree. Three do not:

- **M16 hour angle 1h39m.4 vs N19 1h39m.2.** The 1.9-minute offset from the first exposure (M15/N18, 1h30m.7, H.A. 1h37m.3) implies 1h39m.2, so M16 is probably misprinted.
- **M9 check-field local time 15h27m.7 vs N11 15h27m.5.** Both have H.A. 1h51m.1, which the eclipse-field time predicts for 15h27m.6; undecidable.
- **N9 eclipse barometer 30.43 vs M7 30.34.** Probably a transposition; the other Tahiti readings are 30.24–30.34.

All are kept as printed.

### 3.10 1923 Table III (p.51), recomputed from Table II

- **Group 7 mean distance.** Printed **2.22**; weighted mean 2°.165 (2.17), unweighted 2°.19. The Theory column of the same row (+.08) corresponds to 2°.17 (the formula gives +0.070 at 2°.22). The weight 21.32 shows the membership is right. Probably a misprint. **Kept.**
- **Group 9 weight.** Printed **22.78**; the 14 member weights sum to 22.81. Every one of these weights is independently confirmed by the 1928 Table 10 (w₁₅ = w − w₅/3). **Kept.**
- **Group 3 D₂ mean.** Printed +.39; computed +0.397. The D₂ entries of stars 26–28 are crossed by a pen stroke in the ADS scan but are legible (+.24, +.32, +.09). **Kept.**
- **Observer means (Table IV).** "Mean for each observer" 1″.60 and 1″.78 are not the 1/p.e.²-weighted means (1.594 and 1.795). The paper used "the weights of the separate least squares solutions" (p.52), which are not printed, so this cannot be checked.
- **Quadrant means.** The "Mean of equatorial/polar quadrants" values (1.63 ± .15, 1.76 ± .18) are not weighted means of the quadrant values (1.649 ± 0.186 and 1.744 ± 0.176). They are probably joint solutions over 53 and 39 stars. Not checkable.

### 3.11 Minor textual points

- **p.138.** The eclipse-field p.e. "from the observed mean star displacements (last column of Table 5)" should read Table 4; Table 5 is the check field.
- **"55 reference stars … more than 5°" (p.136).** Table 6 has 56 stars beyond 5° (5°.04–10°.33). Not every star is measured on every plate. Their Σr²/Σr is 7.25°, matching the "mean distance of 7°.3" (p.137) for a scale fitted to the reference stars.
- **"exceed 0″.1 only in four cases" (p.132).** The four are stars 55, 58, 69 and 115 ✓.
- **p.130.** In the ADS scan the words "designed by Dr. F. E. Ross" (after "wide-angle quadruplet objectives of four inches aperture") are crossed out by hand. This is a reader's mark, not print.

## 4. Things a re-analysis must keep in mind

1. **All displacements are relative.** The eclipse-plate scale was set so that the stars beyond 5° (1928) or beyond 2° (1923) have zero mean displacement. Fit D_r = β R/r − α r (the paper's equation of condition) with α free; do not fit β R/r alone. In 1928, R = 0°.2661 and r is in degrees. The 1923 coordinates are in units of 100 mm = 4508″ = 1°.2522.
2. **Errors are probable errors.** σ = p.e. × 1.4826. Weight 1 corresponds to a p.e. of ±0″.33 per radial component on the five-foot plates (p.138; ±0″.37 from the β scatter, p.151) and ±0″.20 on the 15-foot weight system (p.152).
3. **Rejected measures.** Parenthesised values (17 eclipse-field, 7 check-field) were rejected by the authors. They are flagged, not deleted.
4. **Combining the instruments.** Table 10 already combines the two instruments for stars within 3°.5. Do not add Table 6 and Table 10 together for those stars.
5. **Mean columns.** Table 6 `Dr_mean`/`wt_total` and Table 4 means are the authors' own, rounded to 0.01″. The per-plate columns allow the plate-by-plate fits of Table 9.
6. **Table 9 refit.** Refitting Table 6 per plate with the printed weights reproduces Table 9's weights (≤ 0.007) and α (≤ 0.0012). β agrees to ≤ 0.07, well inside its p.e. ≈ ±0″.35–0″.39.
7. **Times.** The Table 1 (and 1923 Table I) local mean times are in astronomical reckoning, counted from noon. The Wallal eclipse exposures (1h30m.7–1h34m.7 on Sept. 21) were at about 5h28m–5h32m UT (Wallal at 8h2m.7 E). This is consistent with the printed solar hour angle 1h37m.3 W, which needs the Sun's 1922 right ascension of 11h50m.6.

## 5. Unresolved

- §3.2 star 67 Dx_mean and star 102 Dy_mean: which plate value or mean is wrong cannot be decided.
- §3.4 Table 11 group 1 weight 8.92 (vs 13.25).
- §3.5 area IV Dx (.00 vs −0.03) and the magnitude-equation numbers.
- §3.6 the rule behind the printed group distances of Tables 7 and 8.
- §3.9 M9/N11 check-field time (15h27m.7 vs 15h27m.5).
- §3.10 the 1923 observer and quadrant means. They are reproducible only with the unpublished least-squares weights.

## 6. Audit (second reading)

A second transcriber (T5-audit, 20 Oct 2026) re-read every numeric cell of the 15 CSV files of this group blind, from the same 600-dpi JBIG2 page images (extracted with `pdfimages`; crops viewed at about 200–600 dpi, every doubtful digit at full resolution). The readings are in the git-ignored scratch folder `data/scratch/T5-audit/blind/` (raw keyed text in `raw/`, scripts in `scripts/`). Before the reading was complete only the `#` header lines, the CSV header rows and the key columns of the first transcription were looked at. The key-column listing also showed, incidentally, the first value column of three small tables (Table 7 "Mean of all plates" `dist_deg`, Table 8 `dist_deg`, Table 9 `beta`) and the `n15`/`n5` columns of Tables 10 and 12. These cells were nevertheless re-read from the scans, and the arithmetic checks below fix them independently.

**Comparison.** 12,517 cells were compared cell by cell (10,843 non-empty printed values, plus empty cells and rejection flags; key columns excluded), aligned by key (by row for Tables 1, 10, 11 and 12). Result: **0 disagreements**. No row is missing or extra; printed precision and signs agree; the parenthesised (rejected) measures agree (17 eclipse-field measures in 11 stars, 7 check-field). A test with deliberately altered cells showed that the comparison script catches both value and precision differences. **No CSV data cell was changed by the audit.**

**Checks re-run** with independent scripts (`data/scratch/T5-audit/scripts/run_all_checks.sh`, output `checks_final.txt`) on the final files. Every number in the table of §2 was reproduced: Table 6 weight sums 145/145 exact, means ≤ 0.0060″; Table 4 means 287/290 within 0.0101″; Table 4 → Table 6 per-plate D_r 819/828 within a per-cell rounding bound; the RA/Dec → x, y fit (rms 0.0013/0.0009 unit, star 129 off by −0.125); Table 5; Tables 7–12 and the areas; the Table 9 refit; the Table 2/3 chain (23 relations off as printed, 5 with g_a = i_a and h′_a = k′_a); 1923 Tables II and III. Every printed inconsistency of §3 was re-viewed at full resolution in all tables involved (the nine cells of §3.1 in Tables 4 and 6, stars 67, 102 and 123, star 129, Table 11 group 1, area IV, the i_a sign and M17–M9 e′, the three Table 1 pairs, the pen-marked 1923 cells). All are printed as transcribed.

**New external checks** (B.D. positions from VizieR I/122 and the southern B.D. I/119, put through the fitted gnomonic model):

- **Table 4.** 142 of the 143 B.D.-numbered stars reproduce their printed x, y within 0.01 unit. The printed α, δ 1900 agree with the B.D. positions precessed to 1900 (sd 0.03 min and 0.6′; max 0.07 min and 1.6′). This confirms every B.D. number and coordinate reading. The one exception is star 129 (§3.3). The 1855 positions quoted in §3.3 are those of I/122, and BD +4°2583 fits the printed x, y to (+0.0004, −0.0018) unit, as stated there.
- **1923 Table II, star 101 (new printed inconsistency).** 89 of the 90 B.D.-numbered stars reproduce their X, Y. Star 101 is printed **"+1 2639"**, but its X, Y (+1.43, −.04) fit **BD +1°2640** (1855: 11h54m21.2s, +1°20.8′). BD +1°2639 (11h54m14.1s, +1°25.2′) lies 0.065 unit (4.9′) farther north. The same star is No. 99 of the 1928 Table 4 (paired with 101 in Table 10), printed "+1 2640" at x, y = +.240, −.007, which fits BD +1°2640 to 0.001 unit. Both prints were checked at full resolution. **Kept as printed** ("+1 2639"); the identification line of the 1923 table is wrong, not the measures.
- **Table 5 check-field x, y.** All 75 stars match a B.D. star after a 6-parameter affine fit (rms 0.0006/0.0011 unit): 74 within 0.003 unit, and 61c is BD +0°3936 (SAO 123497, V = 5.4; its B.D. magnitude field is coded 30.0). This confirms the x, y readings, which no internal arithmetic checks.
- **Stars common to both papers** (the 54 pairs of Table 10). The magnitudes are identical, as expected from p.132 ("for the stars observed with the 15-foot cameras they were taken from Table II of the former publication"). The positions agree to rms 0.0008 unit after a similarity fit. The B.D. numbers are identical except 1923 No. 101 / 1928 No. 99 (above).

**Markdown and quote review.** Every quotation in `wallal1928_method_and_results.md`, `wallal1923_summary.md` and the `quote` column of `wallal1928_results.csv` was located in the page text, the numbers inside each quote were compared one by one, and doubtful places were viewed on the 600-dpi images. Three errors were corrected:

1. **Courvoisier scale correction.** The p.155 text prints "A scale correction of ″.073 per 1°", not ″.0735. The value ″.0735 appears only in the Table 11 column heading "−″.0735r" (p.153). Fixed in `wallal1928_method_and_results.md` §10 (quote) and §13 (page now 153, 155), and in `wallal1928_results.csv` (quote, page 153;155, comment). The value 0.0735 is kept.
2. **`wallal1923_summary.md`.** The cross-reference "checks are in `wallal1928_notes.md` §9" pointed to the Table 1 section; it now reads "§2 and §3.10".
3. **This file, §2, Table 6 row.** "11 rows lie between 0.0051 and 0.0060" now reads "between 0.005 and 0.0060 (10 of them above 0.0051)".

Everything else that was checked agrees with the pages. This includes the remaining quotes; the 1923 Tables I, III, IV and V, the quadrant solutions and the probable errors; the 1928 screw values, image weights, probable errors, scale-difference tables (pp.138, 145), Table 9 and the combination table (p.151); the star counts (53/66/26 by distance, 17 rejected measures in 11 stars, the four stars with P.M. + Par. > 0″.1); and §3.8 (26 reference stars; the +0.153″ difference uses harmonic weights 1/(1/w₁₅ + 3/w₅); with the 15-foot weights alone it is +0.156″).
