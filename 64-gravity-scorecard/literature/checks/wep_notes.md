# Weak equivalence principle (WEP) timeline: verification notes

Companion to `wep.csv` (39 rows) for atlas project N7. Each number was checked against the
primary paper's own text (abstract, table or equation) wherever that text could be reached;
the exceptions are flagged as `secondary_only` (two rows) or described below. The downloaded
sources are in `data/literature/papers/wep-*/`.

Checked October 2026.

## 1. Conventions used in the CSV

- **Eötvös ratio.** η = 2(a1 − a2)/(a1 + a2) for two bodies falling toward the same attractor.
  Papers quote it signed, and for a plot of precision only the size of the uncertainty
  matters. Older papers use equivalent forms:
  - Eötvös et al. use κ, defined by f′ = f(1 + κ), so κ − κ_Pt ≈ η(X, Pt).
  - Adelberger 1990 uses differences of m_i/m_g.
  - LLR papers use (m_G/m_I)_Earth − (m_G/m_I)_Moon.
- **Values are as published.** Some are in physical units: µGal (1 µGal = 10⁻⁸ m s⁻²) or
  cm s⁻². The `notes` field gives the dimensionless conversion, marked "derived".
- **Superscripts.** In quotes, superscripts are written as `^`, for example `10^-6`.
- **New err_type value.** `3sigma` appears once, for INPOP17a: Viswanathan et al. 2018 label
  that uncertainty "3-σ" in their Table 6. It is not in the original vocabulary list.
- **err_syst.** This column holds a separately quoted second error:
  - systematic error for MICROSCOPE, Asenbaum, Baeßler and Adelberger 2001;
  - simulation error for ALPHA-g.

## 2. Recommended series for the "precision of η vs year" plot

**What to plot.** Plot the total 1σ uncertainty: statistical and systematic errors combined in
quadrature, which is how MICROSCOPE combine their own. Where only a limit is published, plot
the limit and mark it as a limit.

**Caveat on the early points.** The early "limits" (Newton to Renner) are sensitivity
statements, not statistically defined bounds. Footnote this. A factor of 2–3 is invisible on a
plot spanning 12 decades, but say so.

### A. Laboratory and space tests of material WEP (the main curve; record holders in bold)

| year | experiment | plot value | basis |
|---|---|---|---|
| **1687** | Newton, pendulums | 1e-3 | stated sensitivity ("thousandth part") |
| **1832** | Bessel, pendulums | 1.7e-5 | stated limit 1/60 000 |
| **1890** | Eötvös, torsion balance | 5e-8 | stated limit 1/20 000 000 |
| **1922** | Eötvös–Pekár–Fekete (measured 1906–08) | 5e-9 | their stated bound; individual 1σ errors are 1–3e-9 |
| 1935 | Renner | 5e-10 claimed | plot hollow and flagged; errors probably ×2.4–3 too small, giving about 1.5e-9 |
| **1964** | Roll, Krotkov & Dicke (Sun) | 1.0e-11 | ±1.0e-11, called a "probable error" but effectively 1σ; 95% bound 3e-11 |
| **1972** | Braginsky & Panov (Sun) | 0.9e-12 | 95% interval; about 0.4e-12 as 1σ (inferred) |
| 1979 | Keiser & Faller | 4e-11 | secondary only, conflicting; optional |
| 1987 | Niebauer et al. (free fall) | 5.1e-10 | 1σ: 0.50 µGal / g |
| 1990 | Adelberger et al. (Eöt-Wash) | 1.0e-11 | 1σ |
| 1994 | Su et al. (Eöt-Wash) | 1.9e-12 | 1σ, Be vs Al/Cu, Earth source |
| **2008** | Schlamminger et al. (Be–Ti) | 1.8e-13 | 1σ |
| **2012** | Wagner et al. (Be–Al) | 1.3e-13 | 1σ; best laboratory result |
| **2017** | MICROSCOPE, first results | 1.3e-14 | 9 ⊕ 9, as the authors' CQG 2019 states |
| **2022** | MICROSCOPE, final | 2.7e-15 | 2.3 ⊕ 1.5, as the authors' CQG 2022 states; current record |

**Record holders.** Braginsky & Panov 1972 stayed the best η until the 2000s. Tino et al.
2020 credit the Eöt-Wash 2001 Sun test with first going below it; Schlamminger 2008 is the
first clean single-pair result below it. Baeßler 1999 and Adelberger 2001 are Earth–Moon
scaled numbers (see section 3), so do not plot them as test-body η on this curve.

### B. Atom interferometry (separate colour)

| year | experiment | plot value | note |
|---|---|---|---|
| 1999 | Peters, Chung & Chu (Cs vs glass corner cube) | 7e-9 | gravimeter comparison |
| 2004 | Fray et al. (⁸⁵Rb/⁸⁷Rb) | 1.7e-7 | |
| 2010 | Merlet et al. (⁸⁷Rb vs corner cube) | 6.5e-9 | gravimeter comparison; Will's "matter waves" point |
| 2014 | Schlippert et al. (⁸⁷Rb/³⁹K) | 5.4e-7 | first test with two chemical elements |
| 2015 | Zhou et al. (⁸⁵Rb/⁸⁷Rb) | 3.0e-8 | |
| 2020 | Asenbaum et al. (⁸⁵Rb/⁸⁷Rb) | 3.8e-12 | 1.8 ⊕ 3.4; still the atom record |
| 2021 | Zhou et al. (⁸⁵Rb/⁸⁷Rb, selected states) | 6.7e-10 | not a record |

### C. Lunar laser ranging (separate symbol)

These are Earth–Moon *total* equivalence-principle tests: composition (WEP) plus
self-energy (SEP).

| year | analysis | plot value |
|---|---|---|
| 1976 | Shapiro et al. | about 7e-12 (derived from η_N = −0.001 ± 0.015) |
| 2004 | Williams et al. | 1.4e-13 |
| 2012 | Williams et al. | 1.3e-13 |
| 2018 | Viswanathan et al. (INPOP17a) | 7.1e-14, labelled 3σ |
| 2018 | Hofmann & Müller | 5e-14 |
| 2021 | Biskupek et al. | 2.4e-14 |

### D. Antimatter

ALPHA-g 2023 gives a_ḡ/g = 0.75 ± 0.21 (errors in quadrature). Put it in an inset or an
annotation, not on the log-precision curve. It shows the *sign* of antimatter gravity, not a
precision test. As an Eötvös ratio it is 0.29 ± 0.27 (derived).

## 3. Entry-by-entry caveats

### Newton (1687)

Verified in the primary text (Motte's 1729 translation, Book III Prop. VI): "a difference of
matter less than the thousandth part of the whole". This is a sensitivity estimate.

### Bessel (1832)

The paper says "no difference … amounting to the sixty-thousandth part of the whole" (p. 414,
page image checked), which is 1.67e-5. The "2×10⁻⁵" in reviews is this number rounded.
Potter (1923) remarks that Bessel's water result was a possible exception.

### Eötvös 1890 (added)

From the Hungarian original (Akadémiai Értesítő 1, p. 110): for brass, glass, antimonite and
cork the difference is "certainly smaller than one twenty-millionth" (5e-8). For brass vs air
it is below 1e-5.

### Eötvös, Pekár & Fekete (1922)

**Units.** The results are in units of **10⁻⁶**. The summary limit is
"κ … den Wert von 0,005·10⁻⁶ nicht erreicht" = 5×10⁻⁹ (p. 65). The errors are
"mittlere Fehler", Gauss's mean-square error, i.e. standard errors. They are not probable
errors.

**Full table (p. 65): κ − κ_Pt in units of 10⁻⁶, Earth's attraction.**

| material | κ − κ_Pt |
|---|---|
| magnalium | +0.004 ± 0.001 |
| snakewood | −0.001 ± 0.002 |
| copper | +0.004 ± 0.002 |
| water | −0.006 ± 0.003 |
| crystalline CuSO₄ | −0.001 ± 0.003 |
| CuSO₄ solution | −0.003 ± 0.003 |
| asbestos | +0.001 ± 0.003 |
| tallow | −0.002 ± 0.003 |

**Other results in the paper:**
- Chemical reactions: the Landolt Ag₂SO₄ + FeSO₄ reaction gives 0.000 ± 0.002; dissolving
  CuSO₄ gives 0.002 ± 0.002 (both ×10⁻⁶).
- Radium bromide: κ(RaBr₂) − κ_Pt = −0.25 ± 0.50 ×10⁻⁶.
- Sun (meridian) series: κ_magnalium − κ_Pt = +0.006×10⁻⁶, with no error quoted.

**The magnalium result.** It is 4× its own error. EPF dismissed it because the data showed
long runs of nearly constant deviations, i.e. systematics. In 1986 Fischbach et al. (PRL 56,
3) found the EPF differences correlate with baryon number per unit mass. That launched the
"fifth force" episode; the later Eöt-Wash and free-fall experiments in this file ruled it out.

**Dates.** Measurements 1906–1908; the meridian series ran to August 1908. Göttingen
prize 1909; the 1909 faculty report quotes 1/200 000 000. Published only in 1922, after
Eötvös's death.

**Other values in the literature:**
- 3e-9: Braginsky & Panov 1972; Tino et al. 2020.
- 1e-8: Peebles 2016; Nobili & Anselmi 2018; Eötvös's own 1909 London report, per Bod et al.
  1991.
- Will 2014 Fig. 1 plots about 5e-9, read off the figure.

### Renner (1935): corrected

Renner's own German abstract is reprinted verbatim in Bod, Fischbach, Marx & Náray-Ziegler
1991, Appendix II. In it he claims:
- accuracy 1 : 2 000 000 000 = **5e-10** in general, and 1 : 5×10⁹ = 2e-10 in one case;
- mean errors averaging ±0.52×10⁻⁹.

So the draft's "~2e-9" is not his claim. Roll, Krotkov & Dicke (1964) showed that his
time-interpolated readings are not statistically independent. His errors are probably too
small by about ×3 (×2.4 according to Király), and his results lie "suspiciously near zero".
Hence reviews quote about 1–2e-9; Will's Fig. 1 plots about 1.3e-9.

The Hungarian journal (Mat. Természettud. Értesítő 53, 542) was not reachable. Bod et al.
print the volume as "13", which looks like a typo for 53.

### Roll, Krotkov & Dicke (1964)

η(Au, Al) = (1.3 ± 1.0)×10⁻¹¹, and "with 95% confidence |η| < 3×10⁻¹¹".

The abstract calls ±1.0 a **"probable error"**. Numerically it is a 1σ standard error:
- With σ = 1.0, P(|η| < 3) = 95.5%, which matches their 95% statement.
- A true probable error (σ = 1.48) would make the 95% bound about 3.8e-11.

I kept `probable_error` because that is the published wording. Treat it as 1σ when plotting.
Only the abstract was available: the Elsevier text, via INSPIRE.

### Braginsky & Panov (1972)

η = (−0.3 ± 0.9)×10⁻¹² at **95%** confidence. The interval is a Student-t interval from 7
daily film records taken in February 1971. The Russian original appeared in 1971. As 1σ it
would be about 0.37e-12 (my inference).

### Keiser & Faller (1979): secondary only

The primary sources (a 1979 BAPS abstract and the 1982 Marcel Grossmann proceedings) are not
online. Will's figure has 4e-11 at about 1979. Tino et al. 2020 give 1e-10 for 1979 and 4e-11
for 1982. Unresolved; plot it only with a caveat.

### Niebauer, McHugh & Faller (1987)

Read from a CC0 copy of the paper on Zenodo. The table gives the composition term as
+0.13 ± 0.50 µGal, with the 1σ convention stated. That is (1.3 ± 5.1)×10⁻¹⁰ of g, so the
draft's numbers are right. Note that the paper prints only "less than 5 parts in 10¹⁰"
(abstract) and the µGal values.

### Eöt-Wash, 1989–2001

- **Heckel et al. 1989** (PRL 63, 2705) was the first rotating-balance result:
  Δa⊥(Be/Cu) = (0.9 ± 1.7)×10⁻¹¹ cm s⁻². It is not in the CSV.
- **Adelberger 1990.** Δa⊥ values for Be/Al and Be/Cu, and the WEP statement
  (m_i/m_g)(Cu) − (m_i/m_g)(Be) = (0.2 ± 1.0)×10⁻¹¹ at 1σ.
- **Su 1994.** Earth-source η are 1σ, as stated. The Sun-source results are published as Δa
  in cm s⁻². Dividing by the Sun's pull at 1 AU (0.593 cm s⁻²) gives η☉ of order
  (5–10)×10⁻¹².
- **Baeßler 1999.** Δa_CD/a_s = (+0.1 ± 2.7 ± 1.7)×10⁻¹³ (draft numbers correct). It is
  labelled "composition-dependent **Earth-Moon** acceleration": the test-body difference
  scaled by the Earth–Moon core-fraction difference (0.382 − 0.101 = 0.281; Williams et al.
  2009, p. 38). The raw test-body precision is about 3.6× worse.
- **Adelberger 2001** (CQG 18, 2397). The improved Earth–Moon value is
  (1.0 ± 1.4 ± 0.2)×10⁻¹³. The primary is paywalled and its abstract has no numbers. The
  value comes from Williams et al. 2009 and Viswanathan et al. 2018. It is the "WEP part"
  that every LLR WEP/SEP split uses.
- **Schlamminger 2008 and Wagner 2012.** Wagner's Table 3 repeats Schlamminger's Be–Ti
  result unchanged and adds Be–Al (Earth −0.7 ± 1.3; Sun −5.2 ± 4.0; ×10⁻¹³). In both,
  statistical and systematic errors are added in quadrature at 1σ.
- **Smith et al. 2000** (PRD 61, 022001) was a short-range test toward a 3-ton ²³⁸U
  attractor. Its abstract gives a_Cu − a_Pb = (1.0 ± 2.8)×10⁻¹³ cm s⁻² against a
  gravitational pull of 9.2×10⁻⁵ cm s⁻², so η ≈ (1 ± 3)×10⁻⁹ (derived). Not a timeline point.

### Atom interferometry

- **Peters 1999 and Merlet 2010** compare two absolute gravimeters. They are not co-located
  differential drops, so each instrument's absolute accuracy enters. Merlet et al. did not
  present their result as a WEP test; Will (2014) does.
- **Merlet sign.** Merlet's abstract gives the magnitude "(4.3 ± 6.4) µGal". Their Table 1
  gives the signed value CAG − FG5 = −4.3 ± 6.4 µGal at k = 1.
- **Peters 2001.** The full Peters account (Metrologia 38, 25, 2001) is quoted by Merlet as an
  atom–FG5 difference of (7 ± 7) µGal.
- **Conventions not stated.** Fray, Schlippert and the two Zhou papers do not state a
  confidence level. Their uncertainties are statistical or combined Type A + B values,
  effectively 1σ.
- **Search for a newer record.** arXiv was searched to October 2026. Zhou et al. 2021 exists
  but at 6.7e-10 is not a record; Asenbaum 2020 (3.8e-12) remains the atom record.
- **Newer preprints:**
  - hyperfine-state test of ⁸⁷Rb, 2.9e-11 (arXiv 2210.08533, no journal reference);
  - first in-orbit atom test on the China Space Station, (−3.1 ± 4.6)×10⁻⁷ (arXiv 2603.22981).

  Neither is a new record.

### Lunar laser ranging

All the LLR values are the **total** Earth–Moon difference: composition (WEP) plus
gravitational self-energy (SEP, Nordtvedt effect).

**Separated values exist only where a paper subtracts the laboratory WEP part**
(Adelberger 2001, (1.0 ± 1.4)×10⁻¹³):
- Williams 2004/2009: SEP = (−2.0 ± 2.0)×10⁻¹³ and η_N = (4.4 ± 4.5)×10⁻⁴.
- Viswanathan 2018, Table 7: SEP = (−1.4 ± 1.6)×10⁻¹³ and η_SEP = (3.1 ± 3.6)×10⁻⁴. That
  quadrature mixes a "3-σ" LLR error with an unspecified-convention laboratory error.
- Hofmann & Müller 2018: η_N = (−0.2 ± 1.1)×10⁻⁴. This treats the whole signal as SEP.

**Uncertainty conventions differ between groups:**
- The JPL papers do not state one.
- INPOP17a is labelled "3-σ".
- The Hannover results (Hofmann & Müller 2018, Biskupek 2021) state none in the text I could
  read. Viswanathan's Table 6 tags earlier Hannover EP results as "3-σ".

So the LLR uncertainties may not be on the same footing as the laboratory 1σ values. The
Hannover/INPOP ones may be conservative by up to ×3.

**Other notes:**
- **Viswanathan published version.** Viswanathan was read from arXiv v1; the published MNRAS
  version was blocked.
- **Biskupek citation slip.** Biskupek et al. 2021 cite MICROSCOPE as "Δm_g/m_i = (4 ± 12)×10⁻¹⁵".
  That matches neither Ti–Pt result. It looks like the Pt–Pt *reference-sensor* control,
  [+4 ± 4 ± 8]×10⁻¹⁵, with errors added linearly. Do not copy it.
- **Sept 2026 preprint** (arXiv 2609.15303, not peer-reviewed): (0.90 ± 3.18)×10⁻¹⁴ with
  EPM21 and (0.50 ± 3.62)×10⁻¹⁴ with INPOP21.

### MICROSCOPE

**2017.** δ(Ti, Pt) = [−1 ± 9 (stat) ± 9 (syst)]×10⁻¹⁵ from one 120-orbit session. The
statistical error is 1σ; the systematic is an "upper limit evaluation". The CQG 2019 version
repeats the numbers, calls both "1σ", and gives the combination as about 1.3×10⁻¹⁴.

**2022 final.** η(Ti, Pt) = [−1.5 ± 2.3 (stat) ± 1.5 (syst)]×10⁻¹⁵.
- The statistical error is 1σ; the systematic is the "overall systematics upper bound".
- CQG 39, 204009 gives the quadrature combination, 2.7×10⁻¹⁵.
- This is 4.6× better than 2017.
- The Pt–Pt control gave [0.0 ± 1.1 ± 2.3]×10⁻¹⁵.
- The test masses are alloys: TiAlV (90/6/4) and PtRh (90/10).

**Critique.** Nobili & Anselmi (Celest. Mech. Dyn. Astron. 137, 3, 2025; arXiv 2302.06400)
dispute how glitches were handled. Gaps were filled with reconstructed data, up to 35–40% of
some sessions. They also report thermal-noise anomalies. They give no alternative η, but the
issue deserves a footnote.

### ALPHA-g (2023)

a_ḡ = (0.75 ± 0.13 (statistical + systematic) ± 0.16 (simulation)) g, with g = 9.81 m s⁻².
- No confidence level is stated, though the component uncertainties are standard deviations.
- The probability that gravity does not act on antihydrogen is 2.9×10⁻⁴.
- Repulsive gravity is excluded (probability < 10⁻¹⁵).
- The authors stress this is not yet a precision measurement.

## 4. Discrepancies found in secondary sources

| quantity | what secondary sources say | what the primary says |
|---|---|---|
| Eötvös–Pekár–Fekete | 1e-8, 5e-9 or 3e-9, depending on the review | 5e-9 (0.005·10⁻⁶) |
| Bessel | "2e-5" | 1/60 000 = 1.7e-5 |
| Renner | ~2e-9 | claimed 5e-10. A web-search summary even turned 1/2 000 000 000 into "2×10⁻⁹" (a mis-conversion) |
| Potter 1923 | "3×10⁻⁶" (Tino et al. 2020) | agreement at the 1e-5 level, e.g. lead vs brass within 1/120 000, steel 1/50 000 (Proc. R. Soc. A 104, 588; read from a scan). A 1927 follow-up (Proc. R. Soc. A 113, 731) was not checked. Potter is not in the CSV. |
| Keiser & Faller | 4e-11 (1979) per Will's figure | Tino et al.: 1e-10 (1979) and 4e-11 (1982). No primary available. |
| MICROSCOPE | Biskupek 2021: (4 ± 12)×10⁻¹⁵ | see above; matches neither Ti–Pt result |
| Roll, Krotkov & Dicke | — | calls ±1.0e-11 a "probable error", but it behaves as 1σ |

Will 2014 Fig. 1 is schematic. Approximate values read from the figure:

| label | read-off value |
|---|---|
| Eötvös | 4.8e-9 |
| Renner | 1.3e-9 |
| Princeton | 1.0e-11 |
| Moscow | 1.0e-12 |
| Boulder | 3.7e-11 |
| Free-fall | 7.5e-10 |
| Eöt-Wash | 1e-11, 1e-12, 4.6e-13, 2.9e-13, 2.5e-13 |
| Matter waves | 4.6e-9 |

The x-axis spacing is not uniform. Use the CSV values, not the figure.

## 5. Not verified, or verified only partly

- **Keiser & Faller:** no primary source (secondary only).
- **Adelberger 2001:** number from secondary sources only.
- **Renner:** read through a verbatim reprint of his abstract; the journal itself was not
  reachable.
- **Abstract only** (full text paywalled or blocked), so error conventions beyond the abstract
  are unknown:
  - Roll, Krotkov & Dicke 1964;
  - Adelberger 1990;
  - Su 1994;
  - Baeßler 1999;
  - Hofmann & Müller 2018.
- **Viswanathan 2018:** accepted arXiv version only.
- **Will's TEGP Table 2.1** (the classic table of these numbers) was not available, so the
  comparison with Will rests on his 2014 Figure 1.

## 6. Where the sources are

Everything is in `data/literature/papers/wep-*`. Each folder holds a PDF or HTML file plus
extracted text. Download routes:
- arXiv PDFs: about 20 papers;
- APS abstract pages;
- the JETP archive (Braginsky & Panov);
- the MTA REAL repository (EPF 1922 scan, Göttingen 1909 report, Eötvös 1890);
- archive.org scans (Bessel 1832, Potter 1923);
- Wikisource (Newton);
- nature.com (Peters 1999, ALPHA-g);
- INSPIRE and Crossref records (Roll, Krotkov & Dicke; Hofmann & Müller; Biskupek);
- Zenodo (Niebauer);
- tudtor.kfki.hu (Bod et al. 1991 with appendices).
