# N2 literature notes: bending starlight from Soldner to Cassini

Checked on 2026-10-20 for project N2 ("Bending starlight: from the 1919 eclipse to Cassini"). The companion table `literature_reanalyses.csv` has one row per published number: 35 rows covering the 1919 results, their re-analyses, and the Cassini analyses.

**How this note relates to N7.** N7 (`64-gravity-scorecard/literature/checks/gamma.csv` and `gamma_notes.md`) already verified 58 deflection and time-delay results. This note does not repeat that work. It adds:
- the re-analyses of the 1919 data;
- the facts the N2 spec relies on (Soldner, Einstein, Cassini, Gaia);
- a "check against" table, with page and table references.

Other agents are transcribing Dyson, Eddington & Davidson (1920), Harvey (1979) and the Wallal bulletins, so their tables are not copied here.

**Conventions.**
- **p.e.** is a probable error, the old 50% half-width: p.e. = 0.6745σ, so σ = 1.4826 × p.e. Dyson et al. (1920), Campbell & Trumpler, Russell (1920) and Mikhailov (1959) all quote probable errors.
- **s.e.** is a standard error (1σ). Harvey (1979) uses standard errors.
- 1″·98 in old typesetting means 1.98″.
- **"Derived"** marks a number I computed myself. These are not published values. Geometry checks use the project's own Horizons file, `data/raw/horizons/cassini_earth_2002.csv` (heliocentric ICRF vectors, TDB).
- The papers I downloaded are in `data/raw/literature/<name>/`, which git ignores. Papers N7 already holds were read from `64-gravity-scorecard/data/literature/papers/`.

---

## 0. The key numbers at a glance

| what | published value | error type | where |
|---|---|---|---|
| Soldner's Newtonian deflection | ω = 0″,84 | — | Berliner Astron. Jahrbuch für 1804 (printed Berlin 1801), p. 170 |
| Einstein 1911 | "4·10⁻⁶ = 0,83 Bogensekunden" | — | Ann. Phys. 35, 898, p. 908 |
| Einstein 1915 | "1.7″ (statt 0.85″)" | — | Sitzungsber. Preuss. Akad. 1915, 831, p. 834 |
| GR at the limb, modern constants | 1.7512″ (1.7504″ if R☉ = 696 000 km) | — | derived: 4GM☉/(c²R☉) with R☉ = 695 700 km (IAU nominal) |
| Sobral 4-inch, 1919 | 1.98 ± 0.12″ (RA 2.06, Dec 1.94) | p.e. | Dyson et al. 1920, pp. 330–331 |
| Príncipe, 1919 | 1.61 ± 0.30″ | p.e. | Dyson et al. 1920, pp. 328, 331 |
| Sobral astrographic, 1919 | 0.93″ (1.52″ with the scale from refraction and aberration only) | none given | Dyson et al. 1920, p. 312 |
| Sobral plates re-measured (1978, published 1979) | 4-inch 1.90 ± 0.11″; astrographic 1.55 ± 0.34″; combined 1.87 ± 0.13″ | s.e. | Harvey 1979, Table III, p. 198 |
| Least-squares refit of the published 1919 measures | 4-inch 1.99″; astrographic 0.95″ (dispersion 0.90″); Príncipe 1.61″ | — | Gilmore & Tausch-Pebody 2022, Table A1 |
| Wallal 1922 (Lick) | 15-ft 1.72 ± 0.11″; 5-ft 1.82 ± 0.15″; combined 1.75 ± 0.09″ | p.e. | already in N7 (Lick Obs. Bull. 11, 41 and 13, 130) |
| VLBI | γ = 0.99992 ± 0.00012 | formal | Lambert & Le Poncin-Lafitte 2011, A&A 529, A70 (already in N7) |
| Cassini 2002 | γ − 1 = (2.1 ± 2.3) × 10⁻⁵ | 1σ ("at 1-σ", Ashby & Bertotti 2010) | Bertotti, Iess & Tortora 2003, abstract |
| Cassini, independent analysis of the closed-loop data | γ − 1 = (−1.3 ± 5.2) × 10⁻⁵ | 1σ, random | Anderson, Lau & Giampieri 2004, Table I |
| Solar deflection at 90° from the Sun (observer at 1 AU) | 4.07 mas | — | Ni 2016 (arXiv:1611.06025, p. 17); formula in Klioner 2003, eq. (68) |

---

## 1. Re-analyses of the 1919 results

### 1.1 What the 1920 report printed, and the variants published before it

**Dyson, Eddington & Davidson 1920** (Phil. Trans. A 220, 291–333; read from the archive.org scan in `data/raw/dyson1920`).

The final results, on p. 331:
- Sobral 4-inch: "the mean result is 1″·98 with a probable error of about ±0″·12".
- Príncipe: "1″·61. The probable error is about ±0″·30".
- Sobral astrographic: "0″·93 discordant by an amount much beyond the limits of its accidental error. For the reasons already described at length not much weight is attached to this determination."

The 4-inch result, separated by coordinate (p. 330): "From declinations 1″·94, From right ascensions 2″·06". On p. 306 the coefficient α = +0ʳ·100 (in micrometer screw revolutions; 1ʳ = 6″·25 for this instrument) = +0″·625 at 50′ from the Sun's centre, and the Sun's radius was 15′·8, so the limb value is 1″·98. Converting the RA coefficient on p. 306 (0ʳ·105) the same way gives 2.08″ (derived), not the printed 2.06.

The astrographic alternatives are on p. 312:
- "or with the value of the scale 0ʳ·082 from the previous table α = + 0ʳ·024 = 0″·93 at the limb." Here e = 0ʳ·082 is the scale change fitted from the plates themselves (Table IX).
- "It may be noticed that the change of scale arising from differences of refraction and aberration is 0ʳ·020. If this value of e be taken instead of 0ʳ·082 we obtain α = + 0ʳ·039 = + 1″·52 at the sun's limb." The superscript r means screw revolutions; the wording was checked on the page image.
- Fitting each plate by least squares gives a mean of 0″·99.
- The plates "were only measured in declination, as the right ascensions were of little weight" (p. 310).

There is also a mean of 0.86″. Page 310 says: "Thus the mean value of α obtained from all the astrographic plates is 0″·86, a figure considerably less than that obtained from the 4-inch photographs." This is the plain mean of the 18 rows of Table IX (16 plates, two of them measured against two comparison plates each; confirmed in the parallel transcription, `dyson1920_table9_astrographic.csv`). Note that the astrographic screw scale is 1r = 12″·3, not the 4-inch's 6″·25. Earman & Glymour (p. 74) and Gilmore & Tausch-Pebody (Table A1) quote the same 0.86″.

Príncipe (p. 328):
- The accidental p.e. is ±0″·25.
- "There is further the error of ±0″·14 affecting all four results equally, arising from the determination of scale."
- Plate U was reduced at Príncipe and gave 2″·90 ± 0″·87 (p.e.), "practically worthless … ignored" (pp. 326–327).
- Eddington called the Sobral astrographic results "largely vitiated by systematic errors" (pp. 328–329).

Plate counts (pp. 299–300): "19 plates being exposed in the astrographic telescope … and eight in the 4-inch camera". "Of the eight plates taken with the 4-inch lens, seven show seven stars, but No. 6, which was taken through cloud, does not show any."

Radial displacement of the individual 4-inch stars, mean of all plates (p. 331; read from the page image):

| star | 11 | 10 | 6 | 5 | 4 | 2 | 3 |
|---|---|---|---|---|---|---|---|
| calculated (″) | 0.32 | 0.33 | 0.40 | 0.53 | 0.75 | 0.85 | 0.88 |
| observed (″) | 0.20 | 0.32 | 0.56 | 0.54 | 0.84 | 0.97 | 1.02 |

**What was said before the report.** The joint RS/RAS meeting of 6 Nov 1919 is reported in Observatory 42, 389–398 (ADS scan in `observatory1919/`).
- Dyson (pp. 390–391): "The astrographic plates gave 0″·97 for the displacement at the limb when the scale-value was determined from the plates themselves, and 1″·40 when the scale-value was assumed from the check-plates. But the much better plates gave for the displacement at the limb 1″·98, Einstein's predicted value being 1·75."
- Dyson (p. 390) also said "Mr. Davidson secured 15 out of 18 photographs showing the required stellar images" and "Dr. Crommelin with the other lens had 7 successful plates out of 8". The report itself says 19 astrographic plates were exposed.
- Eddington (p. 392): "1″·6 ± 0″·3. The p.e. was determined from the residuals of individual stars."

Two later versions are known only from secondary sources:
- Crommelin, Nature 104, 281 (1919), as quoted by Kennefick 2007: "1.52″".
- Dyson, Nature 106, 786 (1921), as quoted identically by Kennefick 2009 and by Gilmore & Tausch-Pebody: "0.90″ … 1.56″".

**So four pairs of astrographic numbers exist:** 0.97/1.40 (Nov 1919), 0.93/1.52 (the report), 0.90/1.56 (1921), and the Table IX mean of 0.86. None carries an error.

### 1.2 Russell 1920: a least-squares re-solution of the 4-inch data

Russell, H. N. 1920, "Note on the Sobral eclipse photographs", MNRAS 81, 154–164 (ADS scan, OCR plus a visual check of Table IV).

- **Data:** the published 4-inch displacements, plus second-order refraction terms (Turner's formulae, Table I). He shows these terms are negligible, 0″·011 in x and 0″·003 in y.
- **Measurement errors:** the mean error of a measured coordinate difference is 0″·263, a p.e. of ±0″·177 per plate pair (p. 157).
- **Four least-squares solutions:**
  1. theoretical deflection, no distortion;
  2. deflection fitted, no distortion (four constants);
  3. theory plus a linear field distortion;
  4. deflection plus distortion (six constants, as in the report).
- **Table IV** gives the correction to the theoretical 0″·553 at 50′:
  - solution 2: +0″·101 (≈ 2.07″ at the limb, derived; not printed);
  - solution 4: +0″·080;
  - x alone: +0″·118; y alone: +0″·061.
- p. 159: "The deflection at the Sun's limb comes out 2″·00 in solution 4 as against 1″·98 in the Report,—an agreement far within the probable error."
- **Distortion** (p. 161): "S = 0″·26 ± 0″·06 at unit distance (50′)" at position angle 62° ± 7°. These are probable errors ("The probable errors are computed from the weights …"). "The existence of a distortion of the eclipse field appears to be proved" (p. 160). He attributes it to the coelostat mirror (p. 161).
- **Table V, sums of squared residuals** (observed against expected from the measurement errors):
  - theory, no distortion: 0.1623 (0.0990)
  - deflection fitted, no distortion: 0.1303 (0.0891)
  - theory plus distortion: 0.0479 (0.0792)
  - both fitted: 0.0244 (0.0693)
- Russell is also the source for Gilmore & Tausch-Pebody's statement that "an independent confirmation of the report analysis and its accuracy" existed in 1920.
- **Relevance to N2:** a plate model with only shift, scale and rotation (four constants) leaves this distortion in the 4-inch residuals.

### 1.3 Freundlich, von Klüber & von Brunn 1931

Z. Astrophys. 3, 171–198 (read from N7's OCR, `gamma-freundlich1931/`).

- They argue the 1919 4-inch data favour a deflection "noch oberhalb des damals errechneten Wertes E = 2″,0 näher an 2″,2" (still above the value of 2″·0 calculated at the time, closer to 2″·2; pp. 193–194).
- Summary on p. 198: "Englische Expedition 1919, 6 m-Kamera, E = 2″,0 bzw. 2″,2" (2″·0 or 2″·2). No error is given.
- Their re-reduction of the Lick 1922 data gives 2.2″.
- Gilmore & Tausch-Pebody (App. B) note these re-analyses "were not considered credible at the time, or later".

### 1.4 Danjon 1932

"Le déplacement apparent des étoiles autour du soleil éclipsé", J. Phys. Radium 3, 281–301 (HAL jpa-00233102).

- **Data:** only the Sobral 4-inch mean displacements of the 7 stars. Danjon treated the Príncipe plates as too poor (p. 288).
- **Method:** a straight-line fit of displacement against 1/r, with a free scale-like parameter q. The equations are images, so only the text could be read.
- **Results** (p. 289):
  - the radial displacements give slope 1.18 × 1.75″, i.e. 2.06″ (the authors' own value was 1″,98 "avec q = 0");
  - RA alone: 2″,11 with q = −0″,012;
  - Dec alone: 2″,00 with q = −0″,001;
  - the mean of the 7 single-star values is 2″,05 with an "écart moyen" of 0″,20.
- **Abstract:** Greenwich 1919 2″,06; Lick 15-ft 2″,05; Lick 5-ft 2″,07; Potsdam 1929 2″,06. "L'auteur adopte la valeur 2″,06."
- Gilmore & Tausch-Pebody (App. B) describe the free scale parameter as a way to "force the observations to match the 1/r law exactly".

### 1.5 Mikhailov 1959, and Merat et al. 1974

**Mikhailov** (MNRAS 119, 593–608, p. 596; N7's copy).
- **Data:** the deflections as revised by Dyson & Woolley (1937, *Eclipses of the Sun and Moon*) for second-order refraction.
- **Fit:** least squares on Δr = A/r + Br: "A = 2″·07 ± 0″·09 p.e. and the scale correction B = −0″·0158".
- With A fixed at 1.75″, B = +0″·0068.
- **Sums of squared residuals** (in units of (0″·01)²): 360 (A free), 690 (A = 1.75″), 402 (a straight line a + br).
- He concludes the deflection "is larger than Einstein's value if the theoretical hyperbolic law of inverse distance is assumed".

**Merat, Pecker, Vigier & Yourgrau 1974** (A&A 32, 471, Table 1, p. 472).
- 1919: ratio to GR 1.13 ± 0.09 (observer) and 1.18 ± 0.08 (Mikhailov's reduction), B = −0.0158. The error convention is not stated; 0.12″ p.e. converted to 1σ would give 0.10.
- **Their claim** (Table 3): an excess over GR of +0.139 ± 0.033″ at a mean distance of 3.43 R☉, and +0.081 ± 0.024″ at 6.14 R☉. Beyond about 5 R☉ the data agree with GR. They interpret the excess as a "dispersive medium" near the Sun.
- This is a fringe interpretation. Later radio and VLBI results are consistent with GR at small elongations (N7).

### 1.6 Harvey 1979: the plates re-measured (the "1979 remeasurement")

Harvey, G. M. 1979, The Observatory 99, 195–198. The tables are being transcribed by another agent; the method is summarized here from the page images.

- **When and who:**
  - Measured in **1978** on the Zeiss Ascorecord at the RGO (Herstmonceux) by E. D. Clements.
  - Software by A. Murray; reduction by Harvey (Kennefick 2007).
  - Published December 1979, for the Einstein centenary.
- **Plates** (p. 196):
  - Of the seven 4-inch eclipse plates, "one is now missing and another broken. One of the comparison plates also is missing."
  - Selected: 5 eclipse + 5 comparison plates from the 4-inch; 16 eclipse + 12 comparison plates from the astrographic. All ten 4-inch plates were reduced.
  - Of the astrographic plates, 9 eclipse + 9 comparison plates were used, those showing "at least the ten images most commonly appearing".
  - The Príncipe plates appear to be lost (Kennefick 2007).
- **Method** (p. 197):
  - A position catalogue was made for each plate; the comparison plates were averaged into one reference set.
  - For each eclipse plate the fit is λ·δx + ax + by + c = ΔX and λ·δy + dx + ey + f = ΔY, by least squares. δx and δy are Dyson et al.'s predicted 1.75″ displacements, and there are six linear plate constants per plate.
  - Weights are 1/s.e.², i.e. **standard errors**.
  - One astrographic comparison plate "showed a very high standard error of unit weight" and was dropped. The star nearest the Sun was dropped on astrographic plates 2, 6 and 11.
- **Results:**
  - Table II: λ = 1.09 ± 0.06 (4-inch) and 0.88 ± 0.19 (astrographic).
  - Table III: 1.90 ± 0.11″ and 1.55 ± 0.34″.
  - Combined: "1″·87 ± 0″·13, a result which is just within one standard error of the predicted value" (p. 198).
- **Derived context:** 1.55 ± 0.34 is 0.6σ below 1.75″ and 2.0σ above 0.87″. The spec's "also favour Einstein" is fair if phrased as "consistent with Einstein's value and 2σ from the half value".
- **Murray on why 0.93 came out low** (e-mail of 22 Nov 2003, quoted by Kennefick 2007): "the problem with the 1919 measurements was not so much in the quality of the images, but rather in the reduction method, which relied very heavily on the experimental determination of the scale constant e. The Herstmonceux plate reductions of course included both co-ordinates on each plate, giving a much better separation of the plate scale from the deflection".

### 1.7 Earman & Glymour 1980: where the bias accusation comes from

"Relativity and eclipses: the British eclipse expeditions of 1919 and their predecessors", Hist. Stud. Phys. Sci. 11(1), 49–85. Read from the CMU KiltHub scan (doi:10.1184/R1/6492578.v1), OCR'd here. Page numbers are the journal's.

- **Their numbers** (p. 75):
  - "The standard deviations in seconds of arc are approximately: Principe astrographic, 0.444″; Sobral astrographic, 0.48″; Sobral 4-inch, 0.178″."
  - The first and last are the probable errors × 1.4826. The 0.48″ is the scatter of the per-plate astrographic values ("a standard deviation can be calculated from the value of the deflection obtained from each of its plates").
  - p. 74 gives the astrographic as "a mean deflection of 0.86″ … from 18 of the astrographic eclipse plates". It is really 16 plates and 18 comparisons.
- **Key sentences:**
  - p. 75: "The dispersion of the measurements from the Principe astrographic is about the same as the dispersion from the Sobral astrographic. … In all these sets of measurement seem of about equal weight, and it is hard to see decisive grounds for dismissing one set but not the other."
  - p. 75: "the very quality of the Sobral 4-inch results prevent them from constituting an unequivocal confirmation of Einstein's predicted deflection of 1.74″: the mean value is too high and the dispersion too small. The mean of the 4-inch measurements is about 1.3 standard deviations from the Einstein value: if Einstein's value is the true one and the errors are random, about one chance in ten existed of obtaining a mean value as high as that given by the measurements." (One-sided; Gilmore & Tausch-Pebody's 0.20 is two-sided.)
  - p. 76: "If one kept the data from all three instruments, the best estimate of the deflection would have to be somewhere between the Newtonian value and the Einstein value. If one kept only the results of the Sobral 4-inch instrument, the best estimate of the deflection would be 1.98″, significantly above even Einstein's value."
  - p. 77: "Dyson spoke as the voice of both expeditions, and it is striking that he described the Principe expedition but did not mention its results."
  - p. 79: "Eddington won the argument by the power of the reference work."
  - p. 84: "The bias showed in his treatment of the evidence: he repeatedly posed a false trichotomy for the deflection results, claimed the superiority of the qualitatively inferior Principe data, and suppressed reference to the negative Sobral results."
  - p. 85: "Now the eclipse expeditions confirmed the theory only if part of the observations were thrown out and the discrepancies in the remainder ignored; Dyson and Eddington, who presented the results to the scientific world, threw out a good part of the data and ignored the discrepancies. This curious sequence of reasons might be cause enough for despair on the part of those who see in science a model of objectivity and rationality."
- **No refit.** They did no least-squares refit. Their "between the Newtonian value and the Einstein value" is a statement, not a calculation. Gilmore & Tausch-Pebody's Table A4 calculation contradicts it (§1.12).

### 1.8 Collins & Pinch, *The Golem* (1993; 2nd ed. 1998)

**Not accessed.** Archive.org has only lending copies, and search inside the book is unavailable. Two quotes are known only through Gilmore & Tausch-Pebody:
- "Thus Eddington and the Astronomer Royal did their own throwing out and ignoring of discrepancies."
- p. 45 (edition not specified): "In interpreting the observations this way, Eddington seemed to confirm not only Einstein's prediction about the actual displacement, but also his method of deriving the prediction from his theory—something that no experiment can do."

Kennefick (2007) gives a fairer summary of their position: "Collins and Pinch are still very careful in not relying on the claim that Eddington fudged any data, their argument being essentially that there is no such thing as a definitive experiment which resolves all doubts". They build on Earman & Glymour and add no new data analysis.

### 1.9 Kennefick 2007/2012, 2009 and 2019

**"Not only because of theory"** (arXiv:0709.0685, 2007; book chapter in *Einstein Studies* 12, 2012). N7's copy was read. Its archive findings, from the RGO papers in Cambridge University Library:
- Dyson's hand appears in the Sobral reduction; Eddington's does not.
- Dyson was sceptical of relativity.
- The 1.52″ alternative appears both in the report and in Crommelin's Nature summary.
- Harvey's 1978 work was set in motion by Andrew Murray.

**Physics Today 62(3), 37–42 (2009).** PDF from aip.brightspotcdn.com, read in full.
- p. 40: "Accordingly, we must believe that it was Dyson who made the decision to ignore the results from the astrographic plates".
- p. 41: "there are good grounds for believing that Dyson made the scientifically correct decision in choosing to ignore the astrographic data."
- p. 42 table: 4-inch "1.98″ ± 0.18″" (the p.e. converted) against 1979 "1.90″ ± 0.11″"; astrographic 0.93″ against "1.55″ ± 0.34″".
- p. 42: "the reanalysis provides after-the-fact justification for the view that the real problem with the Sobral astrographic data was the difficulty, with the limited means available in 1919, of separating the scale change from the light deflection."
- Errors in this article are listed in §5.

**No Shadow of a Doubt (Princeton UP, 2019).** Not accessed. The publisher's interview page is behind Cloudflare. A WebFetch summary, not verbatim, says Kennefick calls the bias complaints "unconvincing", says the decisions were taken "especially [by] Dyson", and says the 1978 re-analysis "completely vindicated" the original team.

### 1.10 Almassi 2009

"Trust in expert testimony: Eddington's 1919 eclipse expedition and the British response to general relativity", Stud. Hist. Phil. Mod. Phys. 40, 57–67 (doi:10.1016/j.shpsb.2008.08.003). Paywalled, with no abstract in Crossref. It is not a data analysis. Gilmore & Tausch-Pebody (note 13) quote it:

"if we look carefully at the reasons Eddington and his colleagues gave in defence of their data analysis, where and to whom these reasons were presented, and the actual critical response from contemporary observers, I suggest, the case for this allegation is seriously undermined."

### 1.11 Will 2015

"The 1919 measurement of the deflection of light", CQG 32, 124001 (arXiv:1409.7812, read in full). This is a review with no new analysis. It summarises Harvey ("errors reduced by 40 percent"; "1.55±0.34 arcseconds … consistent with general relativity") and Kennefick ("no credible evidence of bias"). It has several factual slips, listed in §5. Its useful line for the Gaia point is (p. 8): "at 90° from the Sun, the deflection is still 4 milliarcseconds".

### 1.12 Gilmore & Tausch-Pebody 2020/2022: the modern statistical re-analysis

"The 1919 eclipse results that verified general relativity and their later detractors: a story re-told", Notes Rec. R. Soc. 76(1), 155–180. Print date 20 Mar 2022 per Crossref. Online 21 Oct 2020 per the arXiv comment and the PDF's © 2020 notice; Crossref lists 21 Oct 2021. doi:10.1098/rsnr.2020.0040; arXiv:2010.13744 (title there: "…results which verified…"). The arXiv file is the published online version. Page numbers below are its pages 1–26, which presumably map to pp. 155–180.

**Data.** "All data used in this analysis are those published in the original eclipse report":
- The 4-inch differential positions (14 equations per plate).
- The Príncipe differential positions (2 eclipse plates × 2 comparison plates × 5 stars).
- For the astrographic, only Table IX summaries plus "measured relative displacements in declination … for five stars on the 16 plates". The full astrographic measures were never published.

**Fit.** Dx = ax + by + c + αEx and Dy = dx + ey + f + αEy, by least squares (MATLAB `lsqr`, a 2019 undergraduate project by I. Webber).

**Table A1** (p. 20): least-squares results against the report.

| instrument | report | least squares |
|---|---|---|
| 4-inch | 1.98 | 1.99 |
| astrographic, five stars | 0.93 / 0.99 | 0.95 |
| Príncipe | 1.61 | 1.61 |

**Table A2** (p. 22): the 4-inch least-squares dispersion is 0.33″, i.e. 0.13″ on the mean. The authors keep 0.18″ (the p.e. converted). Príncipe: 0.45″. Sobral astrographic: dispersion 0.90″, "no model-independent way … to derive a standard error on the mean".

**Tables 1 and A3** (two-sided normal test):

| data set | value | vs Einstein | vs Newton | vs zero deflection |
|---|---|---|---|---|
| 4-inch | 1.98 ± 0.18 | z = 1.28, P = 0.20 | z = 6.17 | z = 11.0 |
| Príncipe | 1.61 ± 0.45 | P = 0.76 | P = 0.10 | z = 3.58 |
| Sobral astrographic | 0.95 ± 0.9 | P = 0.37 | P = 0.93 | P = 0.29 |

**Table A4** (p. 22), combinations the authors call "indicative, not robust":

| combination | result | from Einstein | from Newton |
|---|---|---|---|
| 4-inch + Príncipe | 1.93 ± 0.17 | 1.06σ | 6.2σ |
| all three | 1.90 ± 0.16 | 0.94σ | 6.4σ |
| all three, astrographic error = Earman & Glymour's 0.48 | 1.82 ± 0.16 | 0.44σ | 5.9σ |

I checked the arithmetic of these combinations (derived).

**Conclusion** (p. 15): "The original 1919 analysis is statistically robust, with conclusions validly derived, supporting Einstein's prediction. The rejected third data set is indeed of such low weight that its suppression or inclusion has no effect on the final result for the light deflection, though the very large and poorly quantified systematic errors justify its rejection." Their diagnosis of Earman & Glymour (p. 13): they "failed to understand the difference between the dispersion of a set of measurements and an uncertainty, random plus systematic, on the value of the parameter being measured."

**Also useful:**
- Everitt's (1980) "0.93 ± 0.05 arcsec" is an invented error (note 30, and N7).
- Danjon 1932 (§1.4).
- Russell 1920 (§1.2).
- Hawking's *Brief History of Time* remark, corrected by Wayman & Murray 1989 (Observatory 109, 189).
- The Wallal confirmation.

### 1.13 Searched for and not found

- I found no re-reduction of the 1919 plates against Gaia positions.
- Lemos 2019 (arXiv:1912.05587) is a historical review with no new fit.
- **Not read:**
  - von Klüber 1960 (Vistas Astron. 3, 47; the ADS file in `data/raw/wallal` is a 243-byte failed download);
  - Ludendorff 1932 (AN 244, 321); Trumpler 1932 (Z. Astrophys. 4, 208); Jackson 1931 (Observatory 54, 292);
  - Moyer 1979;
  - Hentschel 1994.

### 1.14 A fair one-paragraph summary for the README

The 1919 data reject the half deflection: every combination is ≥ 5.9σ from 0.87″ (Gilmore & Tausch-Pebody 2022). They do not measure the full deflection better than about 10%. The 4-inch result sits 1.3σ above 1.75″, and every independent re-reduction of it gives 1.9–2.1″:
- Russell 1920: 2.00
- Danjon 1932: 2.06
- Mikhailov 1959: 2.07 ± 0.09 p.e.
- Gilmore & Tausch-Pebody: 1.99
- Harvey 1979: 1.90 ± 0.11 s.e.

The rejected astrographic plates give anything from 0.86″ to 1.56″ depending on how the scale is treated. Re-measured in both coordinates, they give 1.55 ± 0.34″ (s.e.). The rejection was decided at Greenwich under Dyson, not by Eddington. The bias claim (Earman & Glymour 1980, popularized by Collins & Pinch) rests on treating a scatter as an uncertainty.

---

## 2. Facts used in the N2 spec

### 2.1 Soldner: year and value

Source: Wikisource transcription, proofread, of "Ueber die Ablenkung eines Lichtstrals von seiner geradlinigen Bewegung, durch die Attraktion eines Weltkörpers, an welchem er nahe vorbei geht". Saved in `soldner1801/`.

- **Where and when:** *Astronomisches Jahrbuch für das Jahr 1804*, pp. 161–172, dated at the top "Berlin, im März 1801".
- **Printing year:** the volume's title page (archive.org `bub_gb_PN43AAAAMAAJ`) reads "Berlin, 1801". Wikimedia Commons also dates it 1801. Jaki 1978's title says "published in 1801". Treder & Jackisch 1981 say "In March 1801 BODE accepted a paper of him".
- **The value** (p. 170): "Wenn man in der Formel für tang ω die Beschleunigung der Schwere auf der Oberfläche der Sonne substituirt, und den Halbmesser dieses Körpers für die Einheit annimmt, so findet man ω=0″,84." (Substituting the Sun's surface gravity into the formula for tan ω, with the Sun's radius as the unit, one finds ω = 0″,84.) His formula is tan ω = 2g/(v√(v² − 4g)), p. 169.
- **His conclusion** (p. 171): the effect can be neglected in current practice.
- **The factor of 2.** Soldner's ω is formally the angle for one arm of the hyperbola; he doubles it explicitly only for the Moon case (p. 170). Treder & Jackisch 1981 (AN 302, 275; ADS scan in `treder1981/`) resolve this: Soldner's g is the old "g = s·t⁻²" (half the modern acceleration), and there are two misprints, including a missing factor 2 on p. 170. Their verdict: "it is indisputable that SOLDNER obtained the Newtonian value of the deflection of light, which with respect to the constants of his times amounts to 0″.84, and not to EINSTEIN's value." The modern-constant Newtonian value is 2GM☉/(c²R☉) = 0.8756″ (derived, R☉ = 695 700 km).
- **Context:**
  - Lenard reprinted the paper in Ann. Phys. (4) 65, 593 (1921) (Treder & Jackisch; Will 2015).
  - Cavendish's unpublished calculation of about 1784 was printed in 1921 in his *Scientific Papers*, vol. II, p. 437, edited by Dyson (Treder & Jackisch). Will 1988 (AJP 56, 413) was not read.
- **Verdict:** "Soldner 1804" is the volume year, not the year of writing or printing. Write "Soldner 1801".
  - Will 2014/2015's "1803" is wrong.
  - Kennefick 2009's "published in 1804" follows the volume year.

### 2.2 Einstein 1911 and 1915

**1911: "Über den Einfluß der Schwerkraft auf die Ausbreitung des Lichtes", Ann. Phys. (4) 35, 898–908.** PDF from Zenodo record 1424213, with p. 908 checked on the page image.
- The formula is α = 2kM/(c²Δ).
- The sentence: "Ein an der Sonne vorbeigehender Lichtstrahl erlitte demnach eine Ablenkung vom Betrage 4·10⁻⁶ = 0,83 Bogensekunden." (A light ray passing the Sun would accordingly be deflected by 4·10⁻⁶ = 0.83 arcseconds.) Dated "Prag, Juni 1911".
- **0.83″ is confirmed.** It comes from rounding the angle to 4 × 10⁻⁶ rad. Will 2015's "0.875" and Kennefick 2009's "0.87″" are modern-constant values, not what Einstein printed.

**1915: "Erklärung der Perihelbewegung des Merkur aus der allgemeinen Relativitätstheorie", Sitzungsber. Preuss. Akad. Wiss. 1915, 831–839, Gesamtsitzung vom 18. November 1915.** Archive.org item `SitzungsberichtePreussAkad19152`, leaf 361, checked on the image.
- p. 834: "Ein an der Oberfläche der Sonne vorbeigehender Lichtstrahl soll eine Ablenkung von 1.7″ (statt 0.85″) erleiden." (A light ray grazing the Sun's surface should be deflected by 1.7″ instead of 0.85″.)
- The formula is 2α/Δ, against the earlier α/Δ.
- p. 831 announces "eine stärkere (doppelt so starke) Lichtstrahlenkrümmung" (a stronger, twice as strong, bending of light rays).
- **1.7″ is confirmed.** Note that in 1915 Einstein recalled his old value as 0.85″, not the 0.83″ he had printed. The doubled prediction dates from 18 Nov 1915, not 1916 as Kennefick 2009 has it.

### 2.3 VLBI: Lambert & Le Poncin-Lafitte 2011

A&A 529, A70 (Research Note). Already verified by N7 from HAL hal-03785812, and re-read here.
- **The value:** γ = 0.99992 ± 0.00012, a formal error ("the formal error in γ decreases when additional sessions are processed").
- **The data are not "quasars passing near the Sun".** They are a global geodetic-VLBI solution: 5055 sessions, 3 Aug 1979 – 30 Aug 2010, "more than 7.3 millions ionosphere-free group delay measurements at 8 GHz", with γ a global parameter. A minority of VLBA Calibrator Survey sessions (1994–2007) "contain observations as close as 1.4° to the Sun", while routine sessions "stopped observing at less than 15° from the Sun after 2002". The deflection is measurable over the whole sky; it is still about 4 mas at 90°.
- **Systematics** (derived from their Table 1): the same 1979–2008 sessions give 0.99984 ± 0.00015 with the 2009 models and 1.00004 ± 0.00015 with the 2011 models. The modelling change alone moves γ by 2 × 10⁻⁴, matching their 2009 conclusion that γ "cannot be estimated at better than 2 × 10⁻⁴". N7 already recommends treating about 2 × 10⁻⁴ as the realistic error.

### 2.4 Cassini: Bertotti, Iess & Tortora 2003 (Nature 425, 374–376)

**What I could read.** The Nature text is paywalled and has no open copy (Unpaywall and Semantic Scholar both report it closed). I read:
- the abstract;
- the three **figure captions** and **figure images** from Nature's public figure pages;
- the three **supplementary files** (Figs S1–S2 and their legends);
- two later papers by the same authors: Bertotti, Ashby & Iess 2008, CQG 25, 045013 (arXiv:0709.1512v5), and Ashby & Bertotti 2010, CQG 27, 145013 (NIST copy);
- the independent JPL analysis by Anderson et al.;
- the PDS documentation of the data set.

All of these are saved in `bertotti2003/`, `ashby2010/`, `anderson2004/` and `pds_cors0028_docs/`.

**Result and error convention.**
- Abstract: "Our result, γ = 1 + (2.1 ± 2.3) × 10⁻⁵, agrees with the predictions of standard general relativity".
- Ashby & Bertotti 2010, eq. (1): "the fit (at 1-σ)". This now states explicitly what N7 inferred from "formal accuracy" in Bertotti et al. 2008.

**Dates.**
- Conjunction on **21 June 2002** (Fig. 1 caption).
- The PDS data set CO-SS-RSS-1-SCE1-V1.0 runs from START_TIME 2002-06-06 to STOP_TIME 2002-07-05.
- Anderson et al. 2004: "carried out on reaction wheels for an interval of 30 days, with no thrusting by the ACS".
- The 18 passes used, with centres read from the residual clusters in Fig. 3 (days from conjunction): −11.6, −10.6, −8.6, −7.6, −6.6, −5.6, −4.6, −3.6, +1.3 (only a handful of points), +2.3, +3.3, +4.3, +5.3, +6.3, +8.3, +9.3, +11.3, +13.3. Supplementary Fig. S1 shows the same 18 passes at day numbers −12 … +13. That puts them on about 9–17 June and 22 June – 4 July (derived).
- No data were taken in the three days before conjunction: "Unfortunately, no data could be acquired for three days just before conjunction owing to a failure of the ground transmitter" (Fig. 2 caption). Anderson et al. 2004 count four passes: the residuals "would be symmetric about conjunction if the Ka transmitter at DSS25 had not failed for four scheduled passes just prior to conjunction".

**Number of passes.** Fig. 2 caption: "the available 18 passages, each lasting about 8 h". Fig. 3: "1,094 data points (compressed at a 300-s integration time)".

**Impact parameters.**
- Fig. 1 caption: "The conjunction—at which the spacecraft (at a geocentric distance of 8.43 AU), the Sun and the Earth were almost aligned, in this order—occurred on 21 June 2002, with a minimum impact parameter b_min = 1.6 R_S, and no occultation."
- Ashby & Bertotti 2010, p. 3: "the minimum value of the impact parameter of the beam was 1.6 R⊙, but in effect only 18 passages have been used, with a minimum impact parameter of ≈6 R⊙".
- Derived from Horizons:
  - the straight-line b_min is 1.667 R☉ at 2002-06-21 13:10 TDB;
  - Cassini was 7.422 AU from the Sun, and Earth–Cassini was 8.438 AU;
  - the ray sweeps about 3.04 R☉ per day, as in Fig. 1;
  - at the pass centres b ≈ 4.4 R☉ (+1.3 d), 7.3 R☉ (+2.3 d), 11.3 R☉ (−3.6 d), and 36–41 R☉ at the ends.
- **So the data cover b ≈ 4 to 41 R☉, not 1.6 R☉.** Anderson et al. 2004: "data within one day of conjunction are useless, even with the Ka-Band uplink".

**Links used.**
- The three observables are "X/X, X/Ka, Ka/Ka" (Supplementary Fig. S1 legend). Frequencies from Anderson et al. 2004:
  - X-band: 7175 MHz up, 8425 MHz down;
  - Ka-band: 34 316 MHz up, 32 028 MHz down;
  - the Ka/Ka link uses the Ka-band translator built by the Italian Space Agency (Anderson et al. 2003 v1).
- All from DSS 25 at Goldstone, "the only station of the DSN … equipped with a Ka-Band … transmitter".
- With three links, the non-dispersive part y_nd and the up- and downlink plasma terms "are separately determined" (S1 legend).
- **Troposphere.** Anderson et al. 2004 used the DSS 25 Advanced Water Vapor Radiometer data, which are the PD1/PD2 files in the PDS set. The Nature S2 legend confirms a tropospheric calibration without saying how it was done: "after the plasma and tropospheric calibrations have been applied, the residuals are approximately white in frequency".

**Noise level of the residuals.**
- Fig. 3: "The r.m.s. value of the 1,094 data points … is 1.2 × 10⁻⁴ Hz, corresponding to a one-way range rate of 2.2 × 10⁻⁶ m s⁻¹". On the 8.4-GHz carrier that is 1.4 × 10⁻¹⁴ fractional at 300 s (derived).
- Bertotti et al. 2008: σ_y = 0.7 × 10⁻¹⁴ at 1000 s, "with a remarkably Gaussian distribution of the frequency residuals".
- 1.4 × 10⁻¹⁴ × √(300/1000) = 0.78 × 10⁻¹⁴, consistent with white noise (derived).
- Supplementary Fig. S1 shows the per-pass Allan deviation of y_nd at 1000 s at about 1–2 × 10⁻¹⁴, and about 4 × 10⁻¹³ on the +1-day pass.

**How the trajectory was handled.**
- Fig. 3 caption: "obtained from the final 9-parameter fit on the calibrated data".
- Bertotti et al. 2008: "The data have been analyzed with NASA's Orbit Determination Program (ODP) in the Barycentric Celestial Reference System". The ODP models the delay as differenced range.
- **The nine fitted parameters are not listed in anything I could read.** For comparison, the independent JPL fit used 8: the 6 state components, one constant acceleration toward Earth, and γ.

**What Fig. 2 shows.** It is titled "The gravitational signal":
- the predicted two-way relativistic frequency shift y_gr against days from conjunction (−15 to +15);
- the 18 passes drawn as thick segments on the curve, which is antisymmetric: negative before conjunction, positive after;
- the drawn extremes are about ±6 × 10⁻¹⁰ at about ±0.5 day.

The caption also says the data near closest approach were "particularly noisy", and that "A much larger plasma noise was detected in some passes after conjunction, and it was fully removed by the multi-link technique. Remarkably, during this time period, SOHO observations revealed large coronal mass ejections traversing the radio beam." Kopeikin et al. 2007 also quote the peak as 6 × 10⁻¹⁰.

**Our amplitude check (derived; unresolved).**
- I differentiated the standard two-way delay, 2 × (1+γ)GM☉/c³ × ln[(r₁ + r₂ + r₁₂)/(r₁ + r₂ − r₁₂)] with γ = 1, along the Horizons geometry. This gives extremes of ±4.25 × 10⁻¹⁰ with Earth and Cassini taken at the same instant, and ±4.23 × 10⁻¹⁰ with light time, as `shapiro.py` computes it.
- At the pass times it gives 3.0 × 10⁻¹⁰ (+1.3 d), 1.9 × 10⁻¹⁰ (+2.3 d) and 3.4–3.9 × 10⁻¹¹ (the end passes).
- The curve digitized from Fig. 2 is 1.34–1.44 times this at every epoch I checked (days −8 … +12). The far tails, which depend only on GM/c³ and time, are equally high.
- A curve scaled by 1.4 reproduces the drawing almost exactly. The overlay is saved as `data/raw/literature/bertotti2003/fig2_overlay_check.png`: red is the standard two-way signal with v = 3.03 R☉/day and b₀ = 1.67 R☉, blue is the same × 1.4.
- Bertotti et al. 2008 say the signal "attained the value 3 × 10⁻¹⁰ in the best passage", which matches **our** +1.3-day value.
- **Known equation errors:**
  - The Nature equations are acknowledged to be wrong by the authors: "Eq. (1) of [4] refers to the round-trip and has an additional factor 2; in eq. (2) an obvious factor 2 is missing" (Bertotti et al. 2008).
  - "As stated in the supplementary material of [6], in equation (2) the two terms on the right-hand side should obviously be multiplied by a factor 2. This error, of course, had no consequence on the computer fit." (Ashby & Bertotti 2010, footnote 4)
  - Neither factor explains a ratio of 1.4.
- **Recommendation:** compute y_gr from first principles for the N2 figure, and do not scale it to the published Fig. 2.

**Delay at the passes (derived).** The two-way excess delay peaks at 258.4 µs with light time, as `shapiro.py` computes it (258.6 µs if Earth and Cassini are taken at the same instant, which puts both legs at b_min). It is 221 µs at the +1.3-day pass and about 132–137 µs at the end passes. Will 2015's eq. (4) gives 260 µs at b_min. Doppler data measure only the rate of change, never the absolute delay.

**Independent analysis of the same data: Anderson, Lau & Giampieri.**
- **2003 preprint** (arXiv:gr-qc/0308010v1):
  - γ − 1 = (−4.8 ± 5.7) × 10⁻⁵ from 1569 closed-loop points;
  - "withdrawn at the recommendation of the Cassini Radio Science Team" (v2, 24 Oct 2003).
- **2004 Texas Symposium paper** (eConf C041213, paper 0305; Wayback copy):
  - **data:** JPL navigation's closed-loop Doppler, the same kind of data as the PDS ODFs this project downloaded; 1697 points after plasma calibration, 1575 in 17 passes after AWVR troposphere calibration;
  - **calibration:** the dual-band plasma combination y = 1.0575901 y_Ka − 0.0575901 y_X; the tropospheric noise floor falls from 4.2 × 10⁻¹⁴ to about 2 × 10⁻¹⁴ at 300 s;
  - **fit:** 8 parameters, with the state epoch at 6 June 2002 20:42:04.1848 ET; four batches with a stochastic acceleration ATAR;
  - **Table I:** γ − 1 = (−1.3 ± 5.2) × 10⁻⁵, ATAR = −61.38 ± 1.57 × 10⁻¹³ km s⁻², χ²/(N−8) = 0.966 for N = 1575;
  - residual σ is 3.95 (pre) and 3.98 (post) µm/s;
  - they "deleted data from the first pass immediately after conjunction";
  - they also give ω_BD > 9000 (95%).
- **This is the closest published analogue to a refit from the ODFs.** Its error is 2.3 times larger than Nature's.

**A published dispute.** Kopeikin et al. 2007 (Phys. Lett. A 367, 276; N7's copy) argued that the Sun's barycentric motion could bias γ by ±1.2 × 10⁻⁴. Bertotti, Ashby & Iess 2008 reply that the correction is "several orders of magnitude below our experimental accuracy", and that a two-way link largely cancels it.

### 2.5 Iess et al. 1999: the Doppler method

Iess, Giampieri, Anderson & Bertotti, "Doppler measurement of the solar gravitational deflection", CQG 16, 1487–1502. IOP blocks automated downloads, so only abstract fragments were obtained (via WebFetch):
- γ "can be substantially improved by measuring the Doppler frequency shift of a microwave beam transponded back to the ground by an interplanetary spacecraft near solar conjunction";
- the coronal dispersion "can be essentially overcome using skilful combinations of carriers with different frequencies";
- Cassini's "Ka-band link at 32–34 GHz … makes this possible";
- the noise budget covers "the solar corona, the non-gravitational accelerations, and thermal noise due to solar radio emission";
- "We estimate that an accuracy in gamma of about 10⁻⁵ is achievable".

The physics of the method, from the authors' later papers:
- y = d(Δt)/dt. Near conjunction the one-way shift is −2(1+γ)(GM/c³)(1/b)(db/dt) (Ashby & Bertotti 2010, eq. 5), and twice that for the two-way link.
- "this frequency shift is of order αv, where α is the deflection angle and v a linear combination of the velocities of the end points" (Bertotti et al. 2008).
- "As pointed out by Iess et al., the coefficient of −db/dt is just the angular deflection of the ray path … When thought of as a series of frequency measurements, a Doppler experiment has more in common with light deflection than it does with Shapiro time delay. … the correlation of γ with the spacecraft trajectory is negligible." (Anderson et al. 2003 v1)
- di Stefano, Cappuccio & Iess 2022 (arXiv:2201.05107, eqs. 4–5) print the one-way formula but call it "two-way". The factor-2 ambiguity runs through this literature.

### 2.6 Gaia's light-deflection model

- **The formula:** Klioner 2003, AJ 125, 1580 (arXiv:astro-ph/0107457v4), eq. (68): δ_pN = (1 + γ) GM_A/(c²|r_oA|) cot(ψ_A/2). PPN γ is explicit, and β too elsewhere in the model.
- **Klioner's Table 1:** the Sun's δ_pN = 1.75 × 10⁶ µas at the limb, and ψ_max(1 µas) = 180°. The Sun's deflection exceeds 1 µas anywhere on the sky.
- **4.07 mas at 90°.** Ni 2016 (arXiv:1611.06025, p. 17) states: "it is 4.07 mas at right angles to the solar direction for an observer at 1 AU from the Sun according to GR." From eq. (68): 2GM☉/(c² AU) = 4.072 mas (derived).
- **Gaia's range** (derived). Gaia observes 45°–135° from the Sun (ESA information sheet by de Bruijne 2009: a "45°-radius zone of avoidance … centered on the Sun"; and "The Sun … contributes significantly to light bending even 180° away from its center"). The solar deflection Gaia meets is therefore 9.8 mas at 45° to 1.7 mas at 135°, and about 4.03 mas at 90° from L2.
- **γ in the astrometric solution:** Lindegren et al. 2012 (A&A 538, A78; arXiv:1112.4139):
  - "The formulation adopted for Gaia (Klioner 2003, 2004) is based on the parametrized post-Newtonian (PPN) version of the relativistic framework";
  - γ enters as an optional global parameter, γ = 1 + g₀, with "g₀ = 0 … the baseline case of General Relativity";
  - "the inclusion of g₀ representing the PPN parameter γ considerably slows down the convergence … caused by the relatively strong correlation between the gravitational light deflection by the Sun (proportional to 1 + γ …) and trigonometric parallax";
  - "gravitational light deflection by the Sun is ∼2 × 10⁻⁸ rad".
- **No published Gaia γ.** As N7 found, none exists as of October 2026; the releases assume GR. I did not re-check the DR3 papers on this.

### 2.7 The astrographic "focus change" story

Confirmed in the primary sources. The cause was the coelostat mirror heated by the Sun, not the lens or the focus setting.

- **Dyson et al. 1920, p. 309:** "The images were diffused and apparently out of focus, although on the night of May 27 the focus was good. Worse still, this change was temporary, for without any change in the adjustments, the instrument had returned to focus when the comparison plates were taken in July. These changes must be attributed to the effect of the sun's heat on the mirror, but it is difficult to say whether this caused a real change of scale in the resulting photographs or merely blurred the images."
- **The footnote on the same page** quotes the field note: "May 30, 3 a.m., four of the astrographic plates were developed … there had been a serious change of focus, so that, while the stars were shown, the definition was spoilt. This change of focus can only be attributed to the unequal expansion of the mirror through the sun's heat. The readings of the focussing scale were checked next day, but were found unaltered at 11·0 mm. It seems doubtful whether much can be got from these plates."
- **Before the eclipse** (p. 298), focusing on Arcturus had already shown "serious astigmatism due to the figure of the mirror of the 16-inch coelostat", reduced with an 8-inch stop. The 4-inch images were "not quite perfect" too, so July comparison plates were needed for both instruments.
- **The report's advice for future eclipses** (p. 332): "It will probably be best to discard the use of coelostat mirrors."
- **Crommelin at the 6 Nov 1919 meeting** (Observatory 42, 391–392): "we can only ascribe this to a change of curvature of the cœlostat mirror, due to the Sun's heat; for the focus was good on the stars two days before the eclipse and again on the check-plates taken during July. The small cœlostat used with the 4-inch lens did not suffer from deformation".
- **Dyson's first report,** at the RAS meeting of Friday 1919 June 13 (Observatory 42, 261–262; ADS scan in `observatory1919june/`): "The images are reported to be diffuse. This is probably due to the mirror or the cœlostat in some way. The proper silvering of the mirror under the conditions prevailing before the expedition started proved very difficult." The heat explanation came later. Kennefick 2009's "streaky" is not in this report.
- **Why it matters.** A defocus can change the scale between the eclipse and comparison plates, and scale is nearly degenerate with deflection:
  - Russell 1920 (p. 159): "If the stars were all equidistant from the Sun, the gravitational displacement would obviously be indistinguishable from a certain change in scale; and, in fact, the changes in α and a nearly neutralise one another at the mean distance of the stars from the Sun, 55′."
  - the report: e and α are combined as α + 243e.
  - Fitting the scale from the blurred plates gave 0.93; fixing it gave 1.52.

### 2.8 Wallal 1922

N7's `gamma.csv` has all three, each verified with page references:
- `E1922-wallal-lick-15ft`: 1.72 ± 0.11 p.e. (Lick Obs. Bull. 11, 41, Table IV, p. 51);
- `E1922-wallal-lick-5ft`: 1.82 ± 0.15 p.e. (LOB 13, 130, p. 160 and Table 9, p. 151);
- `E1922-wallal-lick-combined`: 1.75 ± 0.09 p.e. (2:1 weighted, p. 160; "agrees exactly with Einstein's prediction").

The 1928 bulletin is the second instrument, not a re-analysis. The spec's 1.72 ± 0.11″ is a probable error. Campbell & Trumpler also obtained 2″·05 with check-star corrections, but left it open (1923, p. 54; N7). The Freundlich et al. (1931) re-reduction of Lick gives 2.2″, and Danjon gives 2.05″ (15-ft) and 2.07″ (5-ft).

### 2.9 Smaller checks

- **GR at the limb:** 4GM☉/(c²R☉) = 1.7512″ with GM☉ = 1.32712440041 × 10²⁰ m³ s⁻² and R☉ = 695 700 km (the IAU 2015 nominal radius; derived). With R☉ = 696 000 km it is 1.7504″, which is N7's "1.7505". The spec's "1.75″ × (R☉/b)" is right.
- **The two-way Shapiro delay** at the 2002 conjunction peaks at 258.4 µs (derived, with light time; 258.6 µs with both legs at b_min at the same instant). "About a quarter of a millisecond" is right, but no data were taken that close (§2.4).
- **"Sobral had 7":** 7 plates × 7 stars. The 4-inch was a 4-inch object glass of 19-ft focal length, lent by Father Cortie (Dyson et al. 1920, p. 298; Observatory 42, 262 says "20 feet").
- **Plate scales** (Harvey 1979, Table I): 35.9″/mm (4-inch), 59.9″/mm (astrographic).
- **Comparison plates:** Sobral's were taken "on nights between July 14 and July 18" (Dyson et al. 1920, p. 306), about 7 weeks after the eclipse. Príncipe's were taken at Oxford before the expedition.

---

## 3. "Check against" table for the README

These are the published numbers an N2 refit should be graded against. All errors are as published, with the 1σ equivalent in brackets.

| N2 output | grade against | type | reference |
|---|---|---|---|
| Sobral 4-inch refit from Dyson et al.'s per-plate measures | 1.98 ± 0.12″ (σ 0.18); RA 2.06, Dec 1.94 | p.e. | Dyson et al. 1920, pp. 330–331; derivation p. 306 |
| same, other least-squares re-solutions of the same numbers | 2.00″ (6-constant); 2.07″ (4-constant, derived from Table IV); 1.99″; 2.07 ± 0.09″ (one scale term, Dyson & Woolley 1937 data) | —; p.e. | Russell 1920, p. 159; Gilmore & Tausch-Pebody, Table A1; Mikhailov 1959, p. 596 |
| same, from re-measured plates (not reproducible from the 1920 tables) | 1.90 ± 0.11″ | s.e. | Harvey 1979, Table III |
| field distortion on the 4-inch plates | S = 0.26 ± 0.06″ at 50′, PA 62° ± 7° | p.e. | Russell 1920, p. 161 |
| per-star radial displacements, observed against calculated | the p. 331 table (§1.1) | — | Dyson et al. 1920, p. 331 |
| Sobral astrographic, 5-star declination refit | 0.93″ (scale fitted); 1.52″ (scale fixed); 0.99″ (per-plate least squares); 0.95″, dispersion 0.90″ (Gilmore & Tausch-Pebody) | none | Dyson et al. 1920, p. 312; Gilmore & Tausch-Pebody, Tables A1–A3 |
| Sobral astrographic, re-measured | 1.55 ± 0.34″ | s.e. | Harvey 1979, Table III |
| Príncipe refit | 1.61 ± 0.30″ (σ 0.44), incl. ±0.14″ scale | p.e. | Dyson et al. 1920, p. 328; least squares 1.61 in Gilmore & Tausch-Pebody, Table A1 |
| combined 1919 | 1.93 ± 0.17 (4-inch + Príncipe); 1.90 ± 0.16 (all three); 1.87 ± 0.13 (Harvey's two Sobral instruments) | 1σ; s.e. | Gilmore & Tausch-Pebody, Table A4; Harvey 1979, p. 198 |
| Wallal refit | 1.72 ± 0.11 (15-ft); 1.82 ± 0.15 (5-ft); 1.75 ± 0.09 (combined) | p.e. | Lick Obs. Bull. 11, 41 (Table IV, p. 51); 13, 130 (p. 160) |
| Cassini γ from the Nature analysis | γ − 1 = (2.1 ± 2.3) × 10⁻⁵ | 1σ | Bertotti et al. 2003, abstract; 1σ per Ashby & Bertotti 2010 |
| Cassini γ from closed-loop data, i.e. ODFs (best analogue for an N2 refit) | γ − 1 = (−1.3 ± 5.2) × 10⁻⁵; 1575 points in 17 passes; residual σ ≈ 4 µm/s | 1σ | Anderson, Lau & Giampieri 2004, Table I |
| Cassini residual noise | 1.2 × 10⁻⁴ Hz rms at 300 s on 8.4 GHz (2.2 µm/s one-way); σ_y = 0.7 × 10⁻¹⁴ at 1000 s | rms | Nature Fig. 3 caption; Bertotti et al. 2008 |
| Cassini passes | 18 passes of ~8 h; b from ~4 to ~41 R☉; conjunction 21 June 2002; b_min 1.6 R☉ (no data there) | — | Nature Figs 1–3 captions; Ashby & Bertotti 2010 |
| Cassini signal size | compute it: two-way peak ±4.2 × 10⁻¹⁰ (derived); ≈ 3 × 10⁻¹⁰ in the best pass (Bertotti et al. 2008); Fig. 2 draws ±6 × 10⁻¹⁰ | — | §2.4 |
| VLBI γ (history plot) | 0.99992 ± 0.00012 | formal (realistic about 2 × 10⁻⁴) | Lambert & Le Poncin-Lafitte 2011, abstract and Table 1 (S5) |
| Gaia deflection at 90° | 4.07 mas (1 AU) | — | Ni 2016, p. 17; Klioner 2003, eq. (68) |

---

## 4. Proposed corrections to the N2 spec

1. **Soldner.** Change "(Soldner 1804)" to "(Soldner 1801: 0.84″)". The paper was written in March 1801 and printed in 1801 in the *Berliner Astronomisches Jahrbuch für 1804*. A footnote could add that his 0.84″ is the full Newtonian value computed with 1801 constants (Treder & Jackisch 1981).
2. **Error types for 1919.** Say "1.98 ± 0.12″ and 1.61 ± 0.30″ (probable errors; 1σ ≈ 0.18″ and 0.44″)". For "Sobral 4-inch lens (7 stars)", write "(7 plates, 7 stars)".
3. **The astrographic numbers.** Replace "0.93″, which Dyson and Eddington set aside as faulty" with: "0.93″ (or 1.52″ with the scale change taken from refraction and aberration alone; no error quoted). The Sun's heat had distorted the coelostat mirror and blurred the images, and the report treated the result as 'largely vitiated by systematic errors'." Kennefick shows the decision was Dyson's at Greenwich, not Eddington's.
4. **The 1979 re-analysis.** Make it "a 1978 re-measurement (published 1979)". Replace "found the astrographic plates, properly reduced, also favour Einstein" with: "measured the surviving Sobral plates in both coordinates and found 1.90 ± 0.11″ (4-inch) and 1.55 ± 0.34″ (astrographic; standard errors). The astrographic value is consistent with 1.75″ and 2σ from the half value." Cite Harvey 1979, Kennefick 2009, and Gilmore & Tausch-Pebody 2022 (the least-squares refit of the published numbers: 1.93 ± 0.17″, ≥ 5.9σ from 0.87″ however combined).
5. **Wallal.** Add "(probable error)". Optionally add "the second camera pair gave 1.82 ± 0.15″, combined 1.75 ± 0.09″ (Campbell & Trumpler 1928)".
6. **VLBI.** Replace "VLBI of quasars passing near the Sun" with "geodetic VLBI of thousands of radio sources all over the sky, 1979–2010". Add "(formal error; modelling choices move it by about 2 × 10⁻⁴)".
7. **Cassini.** Replace "radio signals to the Cassini spacecraft grazing the Sun in 2002" with: "radio signals to the Cassini spacecraft, 8.4 AU from Earth on the far side of the Sun, in June 2002 (closest approach of the ray 1.6 R☉; the usable data span 4–41 R☉)". Add "(1σ)" to γ − 1 = (2.1 ± 2.3) × 10⁻⁵. Note that Cassini measured the rate of change of the delay (a frequency shift), not the delay itself.
8. **Method.** "a linear plate model (shift, scale, rotation)": also fit the six-constant affine model, as Dyson et al. and Harvey did. Russell (1920) found a significant linear distortion on the 4-inch plates, which a four-constant model leaves in the residuals. Report both.
9. **Data.** Dyson et al. 1920 print full per-plate measures only for the 4-inch and Príncipe. For the astrographic there are only per-plate summaries (Table IX) and the declination displacements of 5 stars (Gilmore & Tausch-Pebody, App. A). The full 12-star astrographic refit cannot be done from the report; only Harvey's re-measurement exists.
10. **Watch out for.** Change "the astrographic's focus shifted in the heat" to: "the astrographic's coelostat mirror warped in the Sun's heat. The images went out of focus during the eclipse, while the focusing scale stayed at 11.0 mm, and the focus was back by the July comparison plates (Dyson et al. 1920, p. 309)."
11. **Cassini figure.** Generate the predicted y_gr curve from first principles; do not copy Nature Fig. 2's amplitude (§2.4). The SCE1 pass list (18 passes) and the 300-s compression come from Bertotti et al. 2003. The PDS ODFs are closed-loop data, comparable to Anderson et al. 2004, not to the Nature analysis, which used open-loop Radio Science receiver data (as Anderson et al. 2004 state).
12. **Gaia.** "4 mas at 90°" is fine. Optionally: "(4.07 mas; Gaia sees 9.8–1.7 mas at its 45°–135° elongations)".

---

## 5. Errors found in secondary sources

So they are not propagated:

- **Will 2015** (CQG 32, 124001):
  - "In 1911, he determined that the deflection … should be 0.875 arcseconds" — Einstein printed 0.83″.
  - "Johann von Soldner in 1803" — 1801.
  - Príncipe "1.60 ± 0.31 arcseconds, or 0.91 ± 0.18" — Dyson et al. printed 1.61 ± 0.30 (p.e.).
  - "eight usable plates showing at least seven stars each" — 7 of the 8 plates.
  - "The nineteen plates taken on a second telescope turned out to be worthless because the telescope apparently changed its focal length just before totality of the eclipse, possibly as a result of heating by the Sun" — 19 were exposed and 16 measured, giving 0.93 or 1.52, "of much less weight", not worthless. The focusing scale was found unaltered; what changed was the figure of the coelostat mirror, heated by the Sun (Dyson et al. 1920, p. 309).
  - "with the spacecraft at 8.43 astronomical units from the Sun" — 8.43 AU is geocentric (Nature Fig. 1); heliocentric is 7.42 AU.
  - Lebach 1995 "(−8 ± 34) × 10⁻⁴" — should be (−4 ± 17) × 10⁻⁴ (already found by N7).
- **Will 2014 LRR:** the same "8.43 astronomical units from the Sun", and "von Soldner 1803". N7's `gamma.csv` row `T2002-cassini` copies "spacecraft at 8.43 AU … per Will 2014" and should say geocentric. Its error type can now be cited as explicit 1σ (Ashby & Bertotti 2010). Both are now corrected in N7's file and in N2's copy, `gamma_history.csv`.
- **Kennefick 2009** (Physics Today):
  - Einstein's 1911 value "0.87 arcsecond" — he printed 0.83″.
  - "In 1916 … Einstein realized" — 18 Nov 1915.
  - "The stars were noticeably streaky, a problem reported by Dyson at a meeting of the Royal Astronomical Society as early as 13 June" — that report says "diffuse" and blames the mirror silvering or coelostat, not heat.
- **Kennefick 2007** (arXiv:0709.0685): "the Greenwich team decided to measure the star positions on the Astrographic plates only in right ascension, and not in declination" — reversed. Dyson et al. 1920, p. 310 says declination only, and the Murray e-mail quoted in the same paper says "(declination)".
- **Gilmore & Tausch-Pebody:** "Eighteen exposures were obtained with the Sobral astrographic". The report (p. 299, and the exposure table) says 19, but Dyson at the Nov 1919 meeting said "15 out of 18 photographs", so the primary sources themselves disagree.
- **Earman & Glymour 1980:**
  - "18 of the astrographic eclipse plates" — 16 plates, 18 comparisons;
  - they use 1.74″ for Einstein's value;
  - their "one chance in ten" is one-sided.
- **Merat et al. 1974:** gives 1919's ratio error as 0.09; the 0.12″ p.e. converts to 0.10 at 1σ.
- **Bertotti et al. 2003** (Nature): the factor-2 errors in eqs. (1) and (2), acknowledged in Bertotti et al. 2008 and Ashby & Bertotti 2010. Fig. 2's curve is about 1.4 times the standard two-way signal (our check, unresolved).
- **di Stefano, Cappuccio & Iess 2022:** "the Cassini experiment was carried out when the spacecraft was at 7.2 Astronomical Units" (Horizons: 7.42 AU). They also label a one-way formula "two-way".

---

## 6. Not verified, or not accessible

- **Bertotti et al. 2003 main text** (paywalled). Missing: the actual eqs. (1) and (2); the list of the nine fitted parameters; the exact pass times and the conjunction epoch. The figure captions, figure images, abstract and supplementary files were read.
- **Iess et al. 1999, full text.** Only abstract fragments; IOP blocks automated access. The Doppler formula above comes from the authors' later papers.
- **Collins & Pinch, *The Golem*.** Quotes only via Gilmore & Tausch-Pebody.
- **Almassi 2009.** Quote only via Gilmore & Tausch-Pebody.
- **Kennefick 2019 book.** Not accessed; the publisher's interview is summarized, not verbatim.
- **Crommelin 1919 (Nature 104, 280) and Dyson 1921 (Nature 106, 786).** Values via Kennefick and via Gilmore & Tausch-Pebody.
- **Dyson & Woolley 1937.** Its revised deflections are known only via Mikhailov.
- **Not read at all:** von Klüber 1960 (Vistas Astron. 3, 47); Moyer 1979; Ludendorff 1932; Trumpler 1932; Jackson 1931; Hentschel 1994; Will 1988 (AJP, on Cavendish); Jaki 1978 (English translation of Soldner); the Tangherlini–Treder correspondence.
- **Danjon 1932.** His equations are images; the values were read from the running text.
- **Gaia DR3.** I did not check the DR3 papers' statement that γ = 1 was fixed (N7: no Gaia γ has been published).

## 7. Files downloaded

All under `56-light-bending/data/raw/literature/`, which git ignores.

| folder | contents |
|---|---|
| `gilmore2022/` | arXiv abstract and PDF, text |
| `will2015/` | arXiv PDF, text |
| `kennefick2009/` | Physics Today PDF; letters page (Nov 2009) |
| `earman1980/` | CMU scan, 300-dpi pages, `cmu_ocr.txt` |
| `russell1920/` | ADS scan, page OCR |
| `danjon1932/` | HAL PDF, text |
| `merat1974/` | ADS scan, OCR |
| `treder1981/` | ADS scan, page images |
| `soldner1801/` | Wikisource page wikitext; archive.org OCR of the 1804 Jahrbuch volume (title page "Berlin, 1801") |
| `einstein1911/` | Zenodo PDF of Ann. Phys. 35, 898 |
| `einstein1915/` | archive.org OCR of Sitzungsberichte 1915 part 2; page image of p. 834 |
| `observatory1919/` | Observatory 42, 389 (Nov 1919 meeting), OCR |
| `observatory1919june/` | Observatory 42, 261 (13 June 1919 meeting), OCR |
| `bertotti2003/` | Nature figure images; captions text; supplementary Figs S1, S2 and legends; `fig2_overlay_check.png` |
| `anderson2004/` | withdrawn arXiv v1 (2003); 22nd Texas Symposium paper (Wayback copy) |
| `ashby2010/` | NIST copy of CQG 27, 145013 |
| `pds_cors0028_docs/` | SCE1 DATASET.CAT, errata, readme |
| `klioner2003/` | arXiv PDF |
| `lindegren2012/` | arXiv PDF |
| `ni2016/` | arXiv PDF |
| `esa_gaia_lightbending/` | ESA information sheet |
| `distefano2022/` | arXiv PDF |
| `arxiv1912.05587/` | Lemos 2019 |
| `iess1999/` | Semantic Scholar metadata only |

Papers read from N7's folders (`64-gravity-scorecard/data/literature/papers/`):
- `gamma-sec-kennefick2007`
- `gamma-bertotti2008`
- `gamma-kopeikin2007`
- `gamma-lambert2011`
- `gamma-mikhailov1959`
- `gamma-freundlich1931`
- `gamma-bertotti2003`
- Will 2014 source (`will2014/x/article.tex`)
