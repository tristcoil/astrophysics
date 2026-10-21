# Dyson, Eddington & Davidson (1920): the Sobral astrographic reduction (pp. 309–312), transcription notes

Source: F. W. Dyson, A. S. Eddington & C. Davidson 1920, *Phil. Trans. R. Soc. A* **220**, 291–333, section III (Sobral), §§ 23–25. These notes cover printed pp. 309–312 (archive.org `philtrans06337895`, leaves 018–021).

| File | Content | Rows |
|---|---|---|
| `dyson1920_table9_astrographic.csv` | Table IX, p. 310: d, e, α in revolutions at 50′ and α at the limb, one row per eclipse/comparison pair, plus the printed Mean row | 18 + Mean |
| `dyson1920_p311_astrographic_dy.csv` | The untitled p. 311 table of measured Dy for stars 5, 4, 3, 6, 10, in long format | 80 (16 plates × 5; plate 15 / star 10 empty) |
| `dyson1920_p311_alpha_plus_e.csv` | The p. 311 list "α + 0·240e = …", one value per eclipse plate. Its header quotes the five condition equations, their eliminations (p. 311) and the p. 312 results | 16 |

## How the pages were read

- I read every number by eye on enlarged crops of the full-resolution JPEG-2000 scans (5120 × 6600): `06337895_0019.jp2` (p. 310), `_0020` (p. 311) and `_0021` (p. 312). Leaf 018 (p. 309) and leaves 013, 015 and 040 (pp. 304, 306 and 331) were used only for context and definitions.
- Tesseract was used as a second opinion. It agreed on every digit except where it is known to confuse 3 and 8, at Dy₅ of plate 10, α of plate 4 and α of plate 11. I re-read those three at higher zoom. They are +0·033, +0·033 and +0·032.
- Typography: the raised point is the decimal point (0·051 = 0.051). A superscript r means revolutions of the screw, and ″ means arcseconds. The CSVs keep the printed precision and signs, and values printed without a sign (0·000, 0·00) stay unsigned.
- Legibility: in the Table IX row for eclipse plate 15, the subscript of the comparison plate is partly broken in the print. It reads as 3, the same as the next row (15₃). It is written `15_3`.

## What the quantities are (from the paper)

- The astrographic object glass was "stopped down to 8 inches, mounted in a steel tube and fed by a 16-inch cœlostat" (p. 309). The images were "diffused and apparently out of focus", and the instrument "had returned to focus when the comparison plates were taken in July" (p. 309). The authors attribute this to "the effect of the sun's heat on the mirror".
- The plates "were measured in the astrographic duplex micrometer, the eclipse photographs being directly compared with the comparison plates taken in July" (p. 309). They "were only measured in declination, as the right ascensions were of little weight" (p. 310).
- Table IX (p. 310) gives "the value of α, the constant of the gravitational displacement, as calculated from the measures; the apparent difference of scale e between the eclipse and comparison plates; d the difference of orientation of the plates given by the measures of y". These values are "in Revolutions at 50′ Distance", with 1ʳ = 12″·3.
- The form of the equations is defined on p. 304 for the 4-inch plates: "dx + ey + f + αE_y (= Dy)". Here f is a zero correction, and "The quantity α denotes the deflection at unit distance (*i.e.*, 50′ from the sun's centre)". The coefficients of the five astrographic condition equations on p. 311 are exactly the declination coefficients for stars 5, 4, 3, 6 and 10 on p. 304, which I checked on the page.
- Conversion to the limb uses the same factor as the 4-inch result. On p. 306 the paper says "This is at a distance 50′ from the sun's centre. At the time of the eclipse the sun's radius was 15′·8". So α at the limb = α[r] × 12″·3 × 50/15·8 = α[r] × 38″·924.

## Checks performed

All checks were rerun from the CSV files as written, read back with `pandas.read_csv(path, comment='#', skipinitialspace=True)`.

| # | Check | Values | Result |
|---|---|---|---|
| 1 | Table IX: α at limb = α[r] × 12.3 × 50/15.8 | 18 rows | **Pass.** Seventeen rows reproduce exactly after rounding. Plate 11 gives 0.032 × 38.924 = 1.2456, printed 1.24, which is within the ±0.0005 r rounding of α (±0.019″). Max \|printed − computed\| = 0.0056″. A single factor between 38.878 and 38.906 ″/r reproduces all 18 rows exactly, so the limb values were probably converted from unrounded α, or with a radius slightly above 15′·80. |
| 2 | Table IX Mean row | 3 means | **Pass.** e: 1.480/18 = 0.0822 → +0.082. α: 0.398/18 = 0.0221 → +0.022. Limb: 15.48/18 = 0.8600 → +0.86. These are plain means of all **18 rows**, counting plates 7 and 12 twice. Averaging per plate (16) would give 0.077, 0.023 and 0.88, none of which is printed. No mean is printed for d. Text check, p. 310: "the mean value of α obtained from all the astrographic plates is 0″·86". |
| 3 | No. of Stars against the text (p. 310) | 18 rows | **Consistent.** Every plate has ≥ 7 stars except plate 15 (6). This matches "With the exception of plates Nos. 15 and 16, taken through clouds, the stars numbered 3, 4, 5, 6, 10, 11 and 12 are shown on all the plates". The maximum of 12 (plate 7 against 14₃) means stars 2–13, since No. 1 "is always drowned in the corona". This agrees with Harvey (1979) Table I, "Maximum no. of images in the field 12". |
| 4 | Elimination algebra on p. 311 | 9 coefficients | **Pass.** ½[(4)+(5)] + (1) gives 0.5635d (printed 0.564), 0.1030e and −0.1315α (printed −0.131). (2)+(3) gives 0.682d, 0.832e and 2.910α exactly. Their difference gives 3.0415α (printed 3.041; the rounded intermediates give 3.041 exactly), 0.729e and 0.1185d (printed 0.118). Dividing by the α coefficient, e/α = 0.729/3.041 = 0.2397, printed 0.240. |
| 5 | Plate-15 variant, star 10 missing (my derivation) | 1 | **Pass.** Use (4) alone in place of the mean of (4) and (5): (2)+(3) − (1) − (4) gives **2.973α + 0.840e + 0.255d = Dy₃ + Dy₄ − Dy₅ − Dy₆**, so α + 0.2825e + 0.0858d. The printed coefficient is "α + 0·283e". |
| 6 | Each "α + 0·240e" recomputed as (RHS − 0.118d)/3.041, with plate 15 as (RHS′ − 0.255d)/2.973 and d from Table IX | 16 plates | **11/16 within ±0.0005 r**, taking plate 7 with the d of 18₃ and plate 12 with the mean d (see D6). Plates 3 and 17 are just outside rounding (by 0.00006 and 0.00001 r) and plates 6, 15 and 2 are clearly off. See the table below and D2–D5. |
| 7 | Mean of the 16 printed values (p. 312: "α + 243e = + 0ʳ·0435") | 16 | **Pass, exactly.** The sum is 0.696 and 0.696/16 = 0.04350. The mean coefficient is (15 × 0.240 + 0.283)/16 = 0.24269 → 0.243. |
| 8 | α = 0.0435 − 0.243 × 0.082 | 1 | **Pass.** 0.02357 → "+ 0ʳ·024". Then 0.024 × 38.924 = 0.934 → "0″·93". The conversion was made from the rounded 0.024; the unrounded 0.02357 would give 0″·92. |
| 9 | α with e = 0ʳ·020 (refraction and aberration) | 1 | **Pass.** 0.0435 − 0.243 × 0.020 = 0.03864 → "+ 0ʳ·039". Then 0.039 × 38.924 = 1.518 → "+ 1″·52"; from the unrounded value it would be 1″·50. |
| 10 | Per-plate least-squares mean "0″·99" (p. 312) | 16 plates | **Approximately reproduced; the procedure is not specified.** I solved each plate's five equations (plate 15: four) for α, e and f, with d fixed at Table IX. The mean α is 0.0254 r = **0.99″** with the first-listed comparison for plates 7 and 12 (14₃, 11₁), 0.97″ with their mean d and 0.96″ with the second (18₃, 17₂). Solving for d as well gives 0.96″. |

### Check 6 in detail (p. 311 values recomputed from the p. 311 Dy and Table IX d)

RHS = Dy₃ + Dy₄ − Dy₅ − ½(Dy₆ + Dy₁₀). For plate 15, RHS′ = Dy₃ + Dy₄ − Dy₅ − Dy₆. Δ = printed − recomputed, in revolutions. One revolution is 38.9″ at the limb, so Δ = 0.001 r is 0.04″ at the limb.

| Plate | Comparison (d) | RHS | Recomputed | Printed | Δ | Verdict |
|---|---|---|---|---|---|---|
| 1 | 18₄ (+0.051) | 0.1770 | 0.05623 | 0.056 | −0.00023 | ok |
| 2 | 18₄ (−0.009) | 0.1370 | 0.04540 | 0.049 | **+0.00360** | **discrepant (D2)** |
| 3 | 18₄ (−0.074) | 0.1325 | 0.04644 | 0.047 | +0.00056 | marginal (D5) |
| 4 | 18₄ (−0.168) | 0.1590 | 0.05880 | 0.059 | +0.00020 | ok |
| 5 | 11₃ (+0.094) | 0.1630 | 0.04995 | 0.050 | +0.00005 | ok |
| 6 | 11₃ (+0.186) | 0.2050 | 0.06019 | 0.059 | **−0.00119** | **discrepant (D4)** |
| 7 | 14₃ (+0.006) | 0.1045 | 0.03413 | 0.036 | +0.00187 | not this comparison |
| 7 | 18₃ (−0.054) | 0.1045 | 0.03646 | 0.036 | −0.00046 | ok |
| 7 | mean (−0.024) | 0.1045 | 0.03529 | 0.036 | +0.00071 | marginal |
| 8 | 14₃ (+0.093) | 0.1495 | 0.04555 | 0.046 | +0.00045 | ok |
| 9 | 17₄ (−0.096) | 0.0965 | 0.03546 | 0.035 | −0.00046 | ok |
| 10 | 17₄ (+0.090) | 0.1560 | 0.04781 | 0.048 | +0.00019 | ok |
| 11 | 11₁ (+0.073) | 0.1465 | 0.04534 | 0.045 | −0.00034 | ok |
| 12 | 11₁ (−0.009) | 0.1725 | 0.05707 | 0.059 | +0.00193 | not this comparison |
| 12 | 17₂ (−0.102) | 0.1725 | 0.06068 | 0.059 | −0.00168 | not this comparison |
| 12 | mean (−0.0555) | 0.1725 | 0.05888 | 0.059 | +0.00012 | ok |
| 15 | 15₃ (+0.111), variant | 0.1110 | 0.02782 | 0.026 | **−0.00182** | **discrepant (D3)** |
| 16 | 15₃ (−0.002) | 0.0720 | 0.02375 | 0.024 | +0.00025 | ok |
| 17 | 17₂ (−0.022) | 0.0810 | 0.02749 | 0.028 | +0.00051 | marginal (D5) |
| 18 | 17₂ (+0.045) | 0.0935 | 0.02900 | 0.029 | −0.00000 | ok |

**The printed values are the ones the authors used.** The 16 printed values sum to exactly 16 × 0.0435, the printed mean. So a discrepancy cannot be a misprint of the result alone: replacing 0.049 by the recomputed 0.045 for plate 2 would change the mean to 0.0433. A discrepancy therefore means that a printed input (Dy or d) differs from the authors' working value, or that the authors slipped in their arithmetic. The printed values are kept everywhere.

**Effect on the result: small.** The mean of the 16 recomputed values is 0.04341 r instead of 0.0435 r. Then α = 0.04341 − 0.243 × 0.082 = 0.0235 r = 0.91″ instead of 0.93″, and with e = 0.020, α = 0.0386 r = 1.50″ instead of 1.52″. Taking the mean d for both plates 7 and 12 gives 0.91″ as well.

## Discrepancies and how they were resolved

**D1. p. 312 prints "α + 243e = + 0ʳ·0435": a misprint, quoted as printed.** The coefficient lacks "0·" and should be 0·243. The mean of the 16 coefficients is (15 × 0.240 + 0.283)/16 = 0.2427, and the next line of the paper uses it as 0.243: 0.0435 − 0.243 × 0.082 = 0.0236 → printed 0ʳ·024. The full-resolution scan shows clearly that nothing stands before "243". It is a typesetting omission, not a different number.

**D2. Plate 2: printed +0.049, recomputed +0.0454 (Δ = +0.0036 r ≈ +0.14″ at the limb). Kept as printed; the source of the difference cannot be identified.** All five Dy and d = −0·009 were re-read at high zoom and are certain as printed. Any one of these single changes restores agreement:
- d = −0.089, −0.099 or −0.109 instead of −0.009. The transposition −0.090 gives 0.0485 → 0.049.
- Dy₅ = 0.548 (printed 0.558), Dy₄ = 0.666 (0.656), Dy₃ = 0.734 (0.724) or Dy₆ = 0.648 (0.668).

Remarkably, plate 12 against 17₂ has d = −0.102, and that d also reproduces 0.049 exactly. The source cannot be identified from the paper.

**D3. Plate 15: printed +0.026 (with 0.283e), recomputed +0.0278 (Δ = −0.0018 r). Kept as printed.** The coefficient 0.283 confirms the variant equation of check 5, and the d coefficient follows from the same equation. These single changes restore agreement:
- Dy₅ = −0.032 to −0.034 (printed −0.038) or Dy₄ = 0.132–0.134 (printed 0.138). In both, an 8 → 3 change would do it, and 3 and 8 are easily confused in this type.
- Dy₃ = 0.101–0.103 (printed 0.107) or Dy₆ = 0.176–0.178 (printed 0.172).
- d = +0.131 (printed +0.111).

None of them can be singled out.

**D4. Plate 6: printed +0.059, recomputed +0.0602 (Δ = −0.0012 r). Kept as printed.** These single changes restore agreement:
- Dy₄ = 0.483–0.485 (printed 0.488) or Dy₃ = 0.473–0.475 (printed 0.478). The 8 → 3 change works in both.
- Dy₅ = 0.167–0.169 (printed 0.164), Dy₆ = 0.567 (printed 0.557) or Dy₁₀ = 0.647 (printed 0.637).

No change in d (+0·186) can do it, because the required d would be 0.204–0.230. The cause is unidentified.

**D5. Plates 3 and 17: Δ = +0.00056 and +0.00051 r, just beyond strict rounding (by 0.00006 and 0.00001 r). Treated as agreement within the authors' computational rounding, not as misprints.** Two different intermediate roundings reproduce the printed values:
- Plate 3: rounding the RHS 0.1325 up to 0.133 gives 0.0466 → 0.047.
- Plate 17: rounding RHS/3.041 = 0.0266 → 0.027 and 0.0388d = 0.0009 → 0.001 separately gives 0.028.

**D6. Plates 7 and 12, which comparison plate? Unresolved, so the CSV does not assign one.** Table IX lists each of these eclipse plates twice: plate 7 against 14₃ and 18₃, plate 12 against 11₁ and 17₂. The p. 311 Dy table and the α + 0·240e list have one row each, and p. 312 speaks of "the 16 photographs". Arithmetically:
- **Plate 7.** The printed 0.036 is reproduced by the d of **18₃** (0.0365). The mean of the two d's gives 0.0353, which is 0.0007 off and comparable with D5. 14₃ gives 0.0341, which is 0.0019 off.
- **Plate 12.** Neither comparison alone reproduces 0.059: 11₁ gives 0.0571 and 17₂ gives 0.0607. **The mean of the two d's (−0.0555) gives 0.0589 → 0.059.**

The most economical reading is that for these two plates the p. 311 Dy row and d are averages over the two comparison plates. Plate 12 needs this. Plate 7 is then 0.0007 off, or else it uses 18₃ alone. The paper does not say, so the choice is left to the analysis.

## Other notes

- The plate-2 Dy values (+0.56 to +0.72) are all much larger than those of the other plates. That is a large zero point f, which cancels in the combination used on p. 311. It is not an error.
- There are no eclipse plates 13 and 14 in Table IX. The 16 plates are 1–12 and 15–18, matching "the means of the 16 photographs" (p. 312) and Harvey (1979) Table I (16 useful astrographic eclipse plates).
- The paper does not explain the comparison-plate references (18₄, 11₃, 14₃, 18₃, 17₄, 11₁, 17₂, 15₃).
- For context, the conclusion on p. 331 calls the astrographic deflection "0″·93 discordant by an amount much beyond the limits of its accidental error" and attaches "not much weight" to it.

## Ambiguities left unresolved

1. Which comparison plate the single p. 311 rows for eclipse plates 7 and 12 refer to (D6).
2. Whether D2–D4 come from misprinted inputs (Dy or d) or from slips in the authors' computation. The candidate single-digit changes are listed above, but none can be singled out.
3. The exact per-plate least-squares procedure behind the p. 312 mean of 0″·99 (check 10).

## Audit (second reading)

A second, independent transcription (blind double-keying) was made on 2026-10-20.

- **Method.** I downloaded the full-resolution JPEG-2000 leaves `06337895_0019`–`_0021` separately and read every cell of the three CSVs on crops at 1.0–2.0×: each Table IX column, the Dy table in blocks of four rows, and both columns of the α + 0·240e list. My reading was written to `data/scratch/T3-audit/blind/` before I opened the value columns of the files here. Blindness was partial in one respect. The protocol let me see these files' `#` header lines and key columns (comparison plate, coefficient of e), and the headers quote a few values: plate 11's α and limb value, the Table IX means and the list of recomputed Δ.
- **Comparison: no disagreement.** All 421 cells of the three CSVs, keys included, are identical as strings, with the same signs and printed precision. That covers the 93 Table IX numbers, the 18 comparison-plate labels, the 79 Dy values plus the plate 15 / star 10 dash, and the 16 coefficients and 16 values. No CSV cell was changed.
- **Ambiguous print.** The only ambiguous print is the broken subscript of the plate-15 comparison plate in Table IX. Both readings give 15₃. No printed arithmetic involves the label, so it stays as read.
- **Cells behind the discrepancies.** Every input of plates 2, 3, 6, 15 and 17 (five Dy values and d) and their printed results were re-read in tight 2× crops of the full-resolution scan. The type is crisp, so D2–D5 are discrepancies in the paper, not misreadings. A script that tries every single-digit change of each input reproduces exactly the alternatives listed in D2, D3 and D4.
- **Tesseract, own run** (`--psm 6` on full-resolution crops). It misread or garbled six cells, mostly by confusing 3 and 8:
  - Dy₁₀ of plate 9: +0·180 read as "130"
  - Dy₄ of plate 15: +0·138 read as "188"
  - the plate-7 value: +0·036 read as "086"
  - the plate-15 coefficient: 0·283 read as "288"
  - α of plate 11: +0·032 read as "082"
  - the plate-1 limb value +1·28 (garbled)

  The eye reading of each is certain at 2×, and the arithmetic agrees with it in every case. With 0·130, plate 9 would come out at 0.044 instead of 0.035. A value of 0·086 would break the sum 0.696 = 16 × 0.0435. 0·288 ≠ 0.840/2.973. And 0·082 would give 3·19″ at the limb, not 1·24″.
- **Checks rerun.** I reran the checks with my own scripts (`data/scratch/T3-audit/scripts/checks.py` and `check_lsq.py`) on the files as they stand, and they reproduce every number in the checks table above:
  - limb factor range 38.878–38.906; means 0.0822, 0.0221 and 0.8600
  - elimination: 0.5635/0.103/−0.1315, 0.682/0.832/2.910 and 0.1185/0.729/3.0415
  - plate-15 variant: 2.973α + 0.840e + 0.255d, giving 0.2825
  - the 16 recomputed values and Δ: 11 of 16 within rounding, with 18₃ for plate 7 and the mean d for plate 12
  - sum 0.696 and mean coefficient 0.2427
  - α = 0.02357 r and 0.03864 r: 0″·92 and 1″·50 from the unrounded values, 0″·93 and 1″·52 from 0.024 and 0.039
  - per-plate least-squares means: 0.987″ (first-listed comparison), 0.974″ (mean d), 0.962″ (second-listed) and 0.958″ (d free)

  Every failure is one of the discrepancies D1–D5.
- **Notes reviewed.** Every quote was checked against pp. 304, 306, 309–312 and 331, and every number in checks 1–10 and D1–D6 was recomputed. One numerical wording was corrected: check 6 and D5 said plates 3 and 17 were "0.0001" outside rounding. They are 0.00006 r and 0.00001 r outside, and both places now say so.
