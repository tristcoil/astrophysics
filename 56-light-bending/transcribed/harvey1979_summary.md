# Harvey (1979): remeasurement of the 1919 Sobral plates, summary and checks

**Source.** G. M. Harvey (Royal Greenwich Observatory) 1979, "Gravitational deflection of light. A re-examination of the observations of the solar eclipse of 1919", *The Observatory* **99**, 195–198 (1979 December). ADS `1979Obs....99..195H`; file `data/raw/harvey1979/1979Obs....99..195H.pdf`.

I read the text and tables on the PDF's native page images (3456 × 5400, about 490 dpi, 1-bit, extracted with `pdfimages`), with high-zoom crops of every table entry. The paper is set in old-style figures, where 1 looks like a small-capital I and 0 like o. Tesseract is unreliable on them: it read the plate number 7 as "4" and 0·91 as "0·01". Every value below was checked by eye. The article begins in the lower half of p. 195; the upper half is the end of a conference discussion on supernova remnants. It ends on p. 198, where R. F. Griffin's "Spectroscopic binary orbits … Paper 29: HD 220007" begins, which is why p. 198 carries Griffin's running head.

**Transcribed tables.** Typography in all three files: 1″·75 = 1.75″, and "±" is written as separate value and error columns.

| File | Content |
|---|---|
| `harvey1979_table1.csv` | Table I (p. 196), the plates taken by the Sobral expedition |
| `harvey1979_table2.csv` | Table II (p. 197), λ ± standard error for each plate, and the weighted means |
| `harvey1979_table3.csv` | Table III (p. 198), displacement at the limb for the five determinations |

## Why the plates were remeasured

> "Most of the plates taken by this expedition have remained in the plate archives at the Royal Greenwich Observatory, and it was decided that, to mark the centenary of Einstein's birth, this important observational material should be re-examined by means of modern measuring and reduction techniques." (p. 196)

Only the **Sobral** plates were remeasured, from both the 4-inch and the 13-inch Astrographic. The Príncipe plates are not discussed.

## Plates (p. 196, Table I)

| Object glass | Size of plates | Scale ″/mm | Max. images in field | Useful plates: eclipse | Useful plates: comparison |
|---|---|---|---|---|---|
| 13-inch | 16 × 16 cm | 59·9 | 12 | 16 | 17 |
| 4-inch | 10 × 8 inches | 35·9 | 7 | 7 | 7 |

> "Of the seven eclipse plates from the 4-inch, one is now missing and another broken. One of the comparison plates also is missing. The Astrographic plates have all survived intact, though several have suffered considerable discolouration, and several of the comparison plates show some variability in background density." (p. 196)

## How the plates were measured (p. 196)

- **Machine.** The plates were measured on the **Zeiss Ascorecord**, by **Mr E. D. Clements**, not on the RGO's automatic machine GALAXY:

  > "The images on many of the plates, both eclipse and comparison, from both instruments, were non-circular, some grossly so, making them more suitable for measurement on the Zeiss Ascorecord, where the informed judgement of the operator can be exercised, than on GALAXY. A few plates were measured on GALAXY, as a trial, and the results confirmed this view. All the plates that appeared of reasonable quality were measured on the Zeiss, by Mr E. D. Clements."

- **Which plates were measured.**

  > "Five eclipse and five comparison plates taken with the 4-inch, and 16 eclipse and 12 comparison plates taken with the Astrographic, were selected."

- **Positions per plate instead of plate pairs.**

  > "In contrast to the method employed by Dyson, Eddington and Davidson, who made differential measurements using pairs of plates, it was decided to produce a set of image positions for each plate so that those from the comparison plates could be combined in a single reference catalogue against which each eclipse plate could be compared."

## Reduction (pp. 196–197)

- **Software.**

  > "The measurements were processed with existing computer programmes, which had been written for the reduction of astrometric plates measured on GALAXY; the sets from the two telescopes were treated separately throughout."

- **Plates used in the reduction.**

  > "All ten of the 4-inch plates were reduced, but of the Astrographic plates only nine eclipse and nine comparison plates were selected, these having at least the ten images most commonly appearing." (p. 196)

  The ten 4-inch plates are the 5 eclipse and 5 comparison plates. Table II lists 9 Astrographic eclipse plates (2, 3, 4, 6, 7, 8, 10, 11, 12) and 5 4-inch eclipse plates (1, 4, 5, 7, 8).

- **Standard coordinates and the reference catalogue.**

  > "The final stage of reduction calculated standard coördinates (expressed in seconds of arc) of each image, from each plate, referred to the centre of one particular comparison plate. One of the Astrographic comparison plates showed a very high standard error of unit weight, so the solution was repeated with this plate excluded." (pp. 196–197)

  > "Average standard coördinates of the star images from all the comparison plates were calculated and compared with the standard coördinates from each eclipse plate." (p. 197)

  *Inference, not stated: the final Astrographic reference catalogue therefore came from 8 comparison plates.*

- **Model.** "The following equations were solved simultaneously by least squares" (p. 197):

  > λ . δx + ax + by + c = ΔX
  > λ . δy + dx + ey + f = ΔY ,

  > "where, for each star, δx, δy are the displacements in the direction of right ascension and declination, respectively, predicted by Dyson *et al.* on the assumption of a radial displacement outwards of 1″·75 at the Sun's limb, the amount indicated by Einstein's theory; ΔX, ΔY are the corresponding measured displacements; x, y are the standard coördinates referred to the centre of the Sun. λ is the linear coefficient of the predicted displacements and a, b, c, d, e and f are linear plate constants." (p. 197)

  - **Both coordinates** (RA and Dec) were used. In 1920 the astrographic plates had been measured in declination only (Dyson et al. p. 310).
  - **Six linear plate constants per eclipse plate**, fitted together with λ: zero points c and f, and the four linear terms a, b, d and e. The scale difference between each eclipse plate and the reference catalogue is therefore a **free parameter of each plate, in each coordinate**. This contrasts with the 1920 astrographic result of 0″·93, where one common scale value (e = 0ʳ·082, the Table IX mean) was imposed on all plates (Dyson et al. p. 312). That is my comparison; Harvey does not discuss the scale.
  - λ = 1 corresponds to Einstein's 1″·75 at the limb, and the deflection is assumed radial and proportional to 1/r, as in the predicted δx and δy.
- **Rejected images.**

  > "On three of the Astrographic plates (numbers 2, 6 and 11) the measured displacements of the image nearest to the Sun were clearly in error, reflecting the difficulty in identifying and measuring a faint image within the corona. This star was excluded from the solutions for those three plates." (p. 197)

- **Weighting.**

  > "Table II shows the values of λ thus obtained and their standard errors. Each value was assigned a weight equal to the reciprocal of the square of its standard error; the weighted mean values of λ and their standard errors for each of the instruments are also given in Table II. It is evident from the errors that values of λ derived from the 4-inch plates are much better determined than those from the Astrographic, as was pointed out by Dyson *et al.*" (p. 197)

## Results

**Table II (p. 197), λ ± standard error.**

| 4-inch plate | λ | | 13-inch plate | λ |
|---|---|---|---|---|
| 1 | 1·08 ± 0·16 | | 2 | 0·77 ± 0·76 |
| 4 | 1·15 ± 0·25 | | 3 | 1·40 ± 0·40 |
| 5 | 1·11 ± 0·09 | | 4 | 0·59 ± 0·57 |
| 7 | 0·85 ± 0·21 | | 6 | 1·22 ± 0·63 |
| 8 | 1·12 ± 0·11 | | 7 | 1·06 ± 0·64 |
| | | | 8 | 1·33 ± 0·81 |
| | | | 10 | 0·82 ± 0·53 |
| | | | 11 | 1·51 ± 0·91 |
| | | | 12 | −0·03 ± 0·44 |
| **Mean** | **1·09 ± 0·06** | | **Mean** | **0·88 ± 0·19** |

> "*Results.* In the equations, λ was the coefficient of the predicted displacements calculated by Dyson *et al.* The mean values of λ must therefore be multiplied by 1″·75 to produce results in seconds of arc for displacements at the Sun's limb. These are given in Table III." (p. 197)

**Table III (p. 198), gravitational displacement at the Sun's limb, in seconds of arc.**

| Determination | Displacement | Error type |
|---|---|---|
| Predicted from Einstein's theory | 1·75 | — |
| 4-inch plates reduced by Dyson *et al.* | 1·98 ± 0·18 | Matches Dyson et al.'s **probable** error ±0″·12 (p. 331) converted to a **standard** error: 0.12/0.6745 = 0.178. The conversion is not stated by Harvey. |
| 4-inch plates measured on the Zeiss | 1·90 ± 0·11 | Standard error |
| Astrographic plates reduced by Dyson *et al.* | 0·93: | No error printed; the colon marks an uncertain value. The text says "with an unspecified, large error". |
| Astrographic plates measured on the Zeiss | 1·55 ± 0·34 | Standard error |

**Conclusion (p. 198, verbatim).**

> "For the 4-inch plates there is no great difference between the value obtained by Dyson *et al.* and that from the new measurements, but the error has been considerably reduced. For the Astrographic plates, however, a significant improvement has been achieved by the new measurements. Where the previous reduction yielded a value of 0″·93 with an unspecified, large error, the new determination is 1″·55±0″·34. This is still a weak result, but does provide support for that from the 4-inch plates. Combining the two fresh determinations, weighted according to their standard errors, gives 1″·87±0″·13, a result which is just within one standard error of the predicted value."

> "*Acknowledgements.* The author acknowledges with thanks the guidance of Mr C. A. Murray throughout the course of this work, and the skill of Mr E. D. Clements who willingly undertook the task of measuring the plates." (p. 198)

**Error type throughout.** Harvey quotes **standard errors**: Table II says "their standard errors", the text says "weighted according to their standard errors" and "within one standard error". Dyson et al. 1920 quote **probable errors**, for example 1″·98 "with a probable error of about ±0″·12" (p. 331). For a Gaussian, the standard error is the probable error divided by 0.6745. Convert before comparing.

## Checks performed

| # | Check | Result |
|---|---|---|
| H1 | Table II 4-inch weighted mean (5 plates, weights 1/σ²) | 1.0903 ± 0.0594, printed **1·09 ± 0·06**. **Pass.** χ² = 1.49 for 4 degrees of freedom. |
| H2 | Table II 13-inch weighted mean (9 plates) | 0.8798 ± 0.1908, printed **0·88 ± 0·19**. **Pass.** χ² = 7.42 for 8 degrees of freedom. |
| H3 | Table III 4-inch Zeiss = 1.75 × Table II mean | 1.75 × 1.0903 = 1.908 and 1.75 × 0.0594 = 0.104, which would round to 1.91 ± 0.10. Printed **1·90 ± 0·11**. **Pass within rounding, at its edge.** Each Table II entry carries ±0.005 of rounding. Searching that range moves the mean × 1.75 over at least 1.897–1.919 and the error over 0.099–0.108, so 1.90 ± 0.11 is attainable. Harvey evidently converted unrounded values; this is not a misprint. |
| H4 | Table III Astrographic Zeiss = 1.75 × Table II mean | 1.75 × 0.8798 = 1.540 and 1.75 × 0.1908 = 0.334 (1.54 ± 0.33). Printed **1·55 ± 0·34**. **Pass within rounding, at its edge**, by the same argument: the attainable ranges are at least 1.522–1.558 and 0.331–0.337. |
| H5 | Table III "4-inch plates reduced by Dyson et al. 1·98 ± 0·18" against Dyson et al. 1920 p. 331, "1″·98 with a probable error of about ±0″·12" | 0.12/0.6745 = 0.178 → 0.18. **Pass**: a probable error converted to a standard error. |
| H6 | Table III "Astrographic plates reduced by Dyson et al. 0·93:" against Dyson et al. pp. 312 and 331 | 0″·93. **Pass.** |
| H7 | Combined value 1″·87 ± 0″·13, inverse-variance mean of 1.90 ± 0.11 and 1.55 ± 0.34 | Mean 1.8668 → **1·87, pass.** Standard error 1/√(1/0.11² + 1/0.34²) = **0.105**, against the printed **0·13**: **not reproduced** (see D1). |
| H8 | Text consistency of "just within one standard error of the predicted value" | (1.87 − 1.75)/0.13 = 0.92σ, consistent with the printed 0·13. With the recomputed 0.105 the difference would be 1.12σ, just outside one standard error. |
| H9 | Plate counts | 4-inch: 7 eclipse plates (Table I), minus one missing and one broken, leaves 5 measured and reduced, matching the 5 rows of Table II. 5 + 5 comparison = "All ten of the 4-inch plates". Astrographic: 16 eclipse (Table I) and 16 measured; 9 reduced, matching the 9 rows of Table II. **Consistent.** |
| H10 | Table I against Dyson et al. 1920 | 16 astrographic eclipse plates ("the 16 photographs", p. 312), at most 12 stars on an astrographic plate (Table IX), and 7 stars on 7 eclipse plates for the 4-inch (pp. 304–306). **Consistent.** |

## Discrepancies and notes

**D1. The combined standard error ±0″·13 does not follow from the stated weighting. Kept as printed.** Inverse-variance weighting of 1″·90 ± 0″·11 and 1″·55 ± 0″·34 reproduces the mean 1″·87 but gives a standard error of 0″·105. The same holds for every other standard route I tried:
- combining the λ means first: 1.0717 ± 0.0567, i.e. 1″·88 ± 0″·10;
- pooling all 14 plates: 1″·88 ± 0″·10;
- scaling by the external scatter: 0″·10 for the two Table III values (χ² = 0.96 for 1 degree of freedom), or 0″·09 if all 14 plates are pooled (χ² = 10.0 for 13).

The text's "just within one standard error" agrees with the printed 0″·13, so 0″·13 is Harvey's own number and not a typesetting error. One formula that reproduces it from the printed numbers is the weighted mean of the two errors, (Σ σᵢ/σᵢ²)/(Σ 1/σᵢ²) = 0.132. That is not the standard error of a weighted mean, and attributing it to Harvey is my conjecture only. For re-analysis, the conventional combination is 1″·87 ± 0″·10.

**D2. Last-digit differences in Table III against Table II × 1.75.** All four numbers (1·90, 0·11, 1·55, 0·34) differ by 0.01 from 1.75 times the recomputed Table II means, once those are rounded (1.908 → 1.91, 0.104 → 0.10, 1.540 → 1.54, 0.334 → 0.33; H3, H4). Against 1.75 times the printed, rounded means (1.9075, 0.105, 1.540, 0.3325), 0·11 is reproduced (0.105 rounded half up), and 1·90, 1·55 and 0·34 differ by 0.0075, 0.010 and 0.0075. All four are within the rounding of Table II (H3, H4), so they were kept as printed and not treated as misprints.

**D3. Reference.** Harvey cites Dyson, Eddington & Davidson as "*Mem. R.A.S.*, **62**, 291, 1919". The version transcribed in this project is *Phil. Trans. R. Soc. A* **220**, 291–333 (1920). I did not check the Memoirs printing or its year, so this is recorded only.

**D4. Plate numbering.** Harvey's plate numbers are not cross-referenced to Dyson et al. His 13-inch numbers (2–12) lie within Dyson et al.'s Table IX eclipse-plate numbers. His 4-inch numbers 1, 4, 5, 7, 8 would fit Dyson's Roman-numbered 4-inch plates I, IV, V, VII, VIII; the paper's plate VI does not appear in its tables. That would make II and III the missing and broken ones. This correspondence is a plausible inference, not stated in the paper.

**Context (not in Harvey).** Dyson et al. themselves (p. 312) noted that if only the refraction and aberration scale change (0ʳ·020) is used, instead of the mean e of 0ʳ·082, the astrographic plates give 1″·52 at the limb. That is close to Harvey's 1″·55 ± 0″·34, which fits the plate scale freely for each plate in both coordinates. Harvey does not mention the 1″·52.

## Audit (second reading)

A second, independent transcription (blind double-keying) was made on 2026-10-20.

- **Method.** I extracted the native page images separately (`pdfimages -png`: 3456 × 5400 pixels, 1-bit) and read Tables I–III on crops at 1.0–1.6×. My reading was written to `data/scratch/T3-audit/blind/` before I opened the value columns here. Blindness was partial:
  - A key-column listing showed me the Table III displacement column and the Table I plate sizes before my reading.
  - The first transcriber's report quoted the Table II means and the Table III values.

  The 28 per-plate Table II entries and the numbers in Table I were read blind.
- **Comparison: no numerical disagreement.** All 96 cells of the three CSVs, keys included, were compared. Every number is identical: 8 in Table I, 32 in Table II, and 5 values and 3 errors in Table III. The only difference is the wording of the free-text note on the row "Astrographic plates reduced by Dyson et al.". Both readings record "0·93:" with no error printed, so the existing note was kept. No CSV cell was changed.
- **Checks rerun.** I reran the checks with my own script (`data/scratch/T3-audit/scripts/check_harvey.py`):
  - H1: 1.0903 ± 0.0594, χ² = 1.49 for 4 degrees of freedom
  - H2: 0.8798 ± 0.1908, χ² = 7.42 for 8
  - H3 and H4: attainable ranges 1.897–1.919 and 0.099–0.108 for the 4-inch plates, 1.522–1.558 and 0.331–0.337 for the Astrographic
  - H5: 0.12/0.6745 = 0.178
  - H7: 1.8668 ± 0.1047
  - H8: 0.92σ and 1.12σ
  - D1: 1.0717 ± 0.0567, i.e. 1″·88 ± 0″·10; the weighted mean of the two errors is 0.132

  All are reproduced. The only failure is the printed ±0″·13 (D1), which is Harvey's own number.
- **Quotes.** Every quotation above was checked against pp. 195–198 and is verbatim. The Table I–III values in this summary match the CSVs.
- **Corrections to this summary.**
  1. D2 said that all four Table III numbers differ by 0.01 from 1.75 times the printed Table II means. That holds for the recomputed means (H3, H4). Against the printed means it does not: 1.75 × 0·06 = 0.105 rounds to 0·11. D2 now says so.
  2. D1's "external scatter: 0″·10" now names the combination it applies to: the two Table III values. Pooling all 14 plates gives 0″·09.
- **Gloss.** "The RGO's automatic machine GALAXY" (under Machine) is a gloss, not Harvey's words. Harvey names GALAXY without describing it.
