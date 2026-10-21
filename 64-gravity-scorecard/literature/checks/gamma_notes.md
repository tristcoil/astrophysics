# PPN γ from light deflection and time delay, 1919–2026: notes on the literature check

These notes go with `gamma.csv`, which has 58 rows. Each row was checked in October 2026 against the primary paper, where I could reach it. The downloaded copies are in `data/literature/papers/gamma-*/`; folders named `gamma-sec-*` hold secondary sources. Section 6 has the plotting recommendation.

## 0. How to read `gamma.csv`

- **year** is the year the observations were made (eclipse, occultation or conjunction). For data sets that span many years (geodetic VLBI, Hipparcos, planetary ephemerides) it is the year of publication. Will's figure is not consistent about this either (see §5).
- **value and errors** are exactly as published, in the quantity the paper used:
  - `deflection_limb_arcsec` is the deflection at the solar limb in arcseconds.
  - `deflection_ratio_to_GR` is the measured deflection divided by the GR prediction, which equals (1+γ)/2.
  - `(1+gamma)/2`, `gamma` and `gamma-1` mean what they say.
  - Old typesetting: 1″·98 means 1.98″, and the German "2″,24" means 2.24″.
- **status** has four values:
  - `verified`: the draft's number matches the paper.
  - `corrected`: the draft was wrong; the notes column says how.
  - `added`: the draft had no entry, or named the paper without giving a number.
  - `secondary_only`: I could only read a secondary source, which is named in the row.
- **Vocabulary I had to extend:**
  - err_type `3sigma`: used for Pitjeva & Pitjev 2013 and the 99.7% intervals of INPOP20a.
  - err_type `interval`: used for INPOP-style "acceptable violation intervals". These give the range of γ over which the ephemeris post-fit residuals change by less than 5, 10 or 25%. They are **not** σ.
  - method `time_delay_pulsar`: used for PSR 1937+21.
  - method `galaxy_lensing`: used for Bolton et al. 2006, which Will includes in his Fig. 5.
- **err_type conventions:**
  - `standard_error` covers "standard error", "s.e.", "mean error (m.e.)" and the German "mittlerer Fehler". All of these are 1σ.
  - `realistic` is the authors' enlarged ("estimated" or "realistic") error. The formal error is given in the notes column.
  - `formal` is the 1σ error from the least-squares fit.

## 1. Probable error versus σ

Until about 1950, astronomers quoted the **probable error** (p.e.). This is the half-width of the interval that contains 50% of the probability, so p.e. = 0.6745σ, or σ = 1.483 × p.e.

- **Results quoted as probable errors:**
  - Dyson et al. 1920 (Sobral ±0″·12, Principe ±0″·30)
  - Campbell & Trumpler 1923 and 1928
  - Dodwell & Davidson 1924
  - Mikhailov 1959 (his re-reductions)
- **Results already at 1σ, under other names:**
  - Freundlich et al. 1931 ("mittlerer Fehler")
  - van Biesbroeck 1950 and 1953 ("m.e.")
  - Harvey 1979 ("standard error")
  - Texas team 1976 ("1σ")
- **The "1.98 ± 0.12 or 0.16?" question.** The paper prints only ±0″·12, and it is a probable error. As 1σ that becomes ±0.18, which is how Harvey (1979) lists it. The "±0.16" found in some secondary texts (for example Goldoni & Stefanini 2020) is not in the paper.
- **Will 2014, Fig. 5 draws the published errors without converting them.** Examples: Principe ±0.17 (= 0.30″/1.75), Lick 1922 ±0.06 and ±0.09. His early optical error bars are therefore about 1.5 times shorter than 1σ bars would be.
- **Reviews can also convert in the opposite direction.** Mikhailov (1959) turned van Biesbroeck's 1σ errors into probable errors (±0.27 m.e. became "±0.18 p.e.", ±0.10 became "±0.07 p.e."), yet still called them "the observer's result as printed".

The optical results converted to (1+γ)/2 ± 1σ are below. These are my derived numbers, not published values. I divided by 1.75″; the GR limb value is 1.7505″, the Lick papers used 1.745″ and Bruns used 1.751″, and the difference does not matter at this precision.

| eclipse / row id | as published | (1+γ)/2 ± 1σ |
|---|---|---|
| 1919 Sobral 4-inch (`E1919-sobral-4inch`) | 1.98 ± 0.12 p.e. | 1.13 ± 0.10 |
| 1919 Principe | 1.61 ± 0.30 p.e. | 0.92 ± 0.25 |
| 1919 Sobral astrographic (Dyson) | 0.93, no error | 0.53 |
| 1919 Sobral 4-inch, Harvey 1979 | 1.90 ± 0.11 s.e. | 1.09 ± 0.06 |
| 1919 Sobral astrographic, Harvey 1979 | 1.55 ± 0.34 s.e. | 0.89 ± 0.19 |
| 1919 Sobral, both instruments, Harvey 1979 | 1.87 ± 0.13 s.e. | 1.07 ± 0.07 |
| 1922 Wallal, Lick 15-ft | 1.72 ± 0.11 p.e. | 0.98 ± 0.09 |
| 1922 Wallal, Lick 5-ft | 1.82 ± 0.15 p.e. | 1.04 ± 0.13 |
| 1922 Wallal, Lick combined | 1.75 ± 0.09 p.e. | 1.00 ± 0.08 |
| 1922 Cordillo Downs | 1.77 ± 0.3 p.e. | 1.01 ± 0.25 |
| 1922 Wallal, Canadian | 1.42 / 1.75 / 2.16, no error | 0.81 – 1.23 |
| 1929 Takengon | 2.24 ± 0.10 (1σ) | 1.28 ± 0.06 |
| 1936 USSR, 1959 re-reduction | 2.68 ± 0.37 p.e. | 1.53 ± 0.31 |
| 1936 Japan (secondary only) | 1.71, no error | 0.98 |
| 1947 Bocaiuva | 2.01 ± 0.27 (1σ) | 1.15 ± 0.15 |
| 1952 Khartoum | 1.70 ± 0.10 (1σ) | 0.97 ± 0.06 |
| 1973 Chinguetti | 0.95 ± 0.11 (1σ, ratio) | 0.95 ± 0.11 |
| 2017 Wyoming (Bruns) | 1.752, ±3.4%, level not stated | 1.00 ± 0.03 |

## 2. The 1919 data and the controversy over which plates were used

Dyson, Eddington & Davidson (1920, Phil. Trans. A 220, 291) report three instruments:

| instrument | result at the limb | data | notes |
|---|---|---|---|
| Sobral 4-inch lens | 1″·98 ± 0″·12 (p.e.) | 7 plates | good images (pp. 306, 331) |
| Sobral 13-inch astrographic | 0″·93, no error given | 16 plates | diffuse, distorted images, measured in declination only |
| Principe (Eddington) | 1″·61 ± 0″·30 (p.e.) | 2 plates through cloud, 5 stars | plate scale from check plates taken at Oxford and Principe, not solved from the eclipse plates (pp. 328, 331) |

For the astrographic plates, the 0″·93 depends on solving for the plate scale from the eclipse plates themselves. If the scale change is taken only from refraction and aberration, the same plates give 1″·52 (p. 312).

The authors set the astrographic result aside as "largely vitiated by systematic errors". They concluded for Einstein's 1″·75 and against the "Newtonian" half deflection of 0″·87.

**The controversy.**

- In 1980 Earman & Glymour, and separately C. W. F. Everitt, argued that the astrographic plates were dropped mainly because they disagreed with Einstein, while the weak Principe result was kept.
- The story reached a wide audience through Collins & Pinch, *The Golem*, and Hawking, *A Brief History of Time*. (Source for this paragraph: Kennefick 2007.)
- Kennefick (2007, arXiv:0709.0685) used the Royal Greenwich Observatory (RGO) archives to argue that the decision was reasonable:
  - The rejection is recorded on the data sheets, apparently in Dyson's hand. Eddington was probably not involved, and Dyson himself was sceptical of relativity.
  - Dyson's first draft even averaged the two Sobral instruments to 1″·83, "very close to the value required by Einstein's theory". Eddington objected that combining a result that disagrees with one that agrees was arbitrary (letter of 21 Oct 1919, quoted by Kennefick).
  - The underlying problem with the astrographic plates was the plate-scale solution.
  - Will cites Kennefick's 2009 Physics Today article, which I did not read; the 2007 arXiv paper is the longer treatment.
- **Harvey (1979)** is the decisive modern evidence. In 1978 the RGO re-measured the surviving Sobral plates in both coordinates and reduced them with modern methods (the Principe plates appear to be lost). The results, all standard errors, were:
  - 4-inch: 1.90 ± 0.11
  - astrographic: 1.55 ± 0.34
  - both combined: 1.87 ± 0.13

  So the "Newtonian" 0.93 came from the original reduction, not from the plates.

**Fair summary.** The 1919 data rule out the half deflection. But they are 10–25% measurements, and neither the original numbers nor the 1979 ones confirm GR to better than about 1σ. Mikhailov's own re-reduction of the 4-inch data gives 2.07 ± 0.09 (p.e.), which shows how much these numbers depend on how the plates were reduced.

## 3. Conflicting values and re-analyses

### 1922, Wallal

- Campbell & Trumpler had two instruments:
  - The 1923 Lick Bulletin covers the 15-ft cameras: 1.72 ± 0.11.
  - The 1928 Bulletin is **not** a re-analysis. It covers the second instrument, the 5-ft cameras (1.82 ± 0.15), and then combines both into 1.75 ± 0.09, with weights 2:1.
  - All three numbers are probable errors.
  - Plot either the two instruments separately or the combination, never all three.
- In 1923 they also worked out 2″·05 using check-star "corrections" (p. 54), but left it open rather than adopting it.
- Freundlich et al. (1931) argued that the Lick data really implied about 2.2; Mikhailov (1959) re-reduced them to 1.83 ± 0.11 (p.e.).
- Two smaller results from the same eclipse:
  - Chant & Young (Canada). Published in 1923 in Publ. Dominion Astrophys. Obs. Victoria 2, 275, not 1924 and not a Toronto publication. Three solutions, 1.42 to 2.16, with no error given.
  - Dodwell & Davidson (Cordillo Downs, MNRAS 84, 150): 1.77 ± 0.3 (p.e.).

### 1929, Takengon (Potsdam expedition)

- 2.24 ± 0.10 (1σ), with the plate scale set by a collimator. This is clearly high, and the authors argued that the excess was real.
- Mikhailov re-reduced the plates, solving for the scale from the stars, and got 1.96 ± 0.08 (p.e.).

### 1936 eclipse

- **USSR.** I could not find the draft's 2.73 ± 0.31 anywhere I had access to.
  - In his own 1959 English-language account, Mikhailov lists 2″·70 as "given by observer" and gives a re-reduction of 2″·68 ± 0″·37 (p.e.).
  - The original Soviet publication was not available.
- **Japan.** Matukuma's 1.71 comes from a secondary source only.

### 1947 and 1952, van Biesbroeck

- Published: 2.01 ± 0.27 and 1.70 ± 0.10 ("m.e.", i.e. 1σ).
- Mikhailov re-reduced the published measurements to 2.20 ± 0.38 and 1.43 ± 0.18 (p.e.). He called the published 1952 value "quite inadmissible" given the measurements.
- His Table I gives a weighted mean over all six eclipses of 1″·93 ± 0″·05 (p.e., from his re-reductions) and a simple mean of 2″·03 ± 0″·10. Eclipse results as a group tend to come out high.

### 1973, Chinguetti (Mauritania)

- Published in two forms:
  - (0.95 ± 0.11) × the GR value (1σ, in both abstracts).
  - 1.66 ± 0.18″ (the final result in Jones 1976).
- 0.18/1.75 = 0.103, so the ratio's error was rounded up to 0.11. This is one measurement, so plot only one of the two rows.

### 2017, Wyoming (Bruns 2018, CQG 35, 075009)

- 1.752″, uncertainty 3.4% (1.69–1.81″), using amateur-class equipment.
- More precise than any earlier eclipse result, and absent from Will's figure.

### Radio, 1974–75 (Fomalont & Sramek)

- Three numbers come from two experiments:
  - 1974: 1.015 ± 0.011 (ApJ 199, 749)
  - 1975: 0.999 ± 0.011 (given in the 1977 review)
  - mean of the two: 1.007 ± 0.009 (PRL 36, 1475)
- The 1977 review restates the mean as γ = 1.014 ± 0.018.
- The draft's "1.015 ± 0.011" belongs to the ApJ paper, not to the PRL.
- Table I of Fomalont & Sramek (1977, Comments Astrophys. 7, 19) is the best single list of the 1969–75 radio results. Two of its entries come from sources I could not read:
  - Goldstone–Haystack VLBI 1970: 1.03 ± 0.11, from Weinberg's 1972 textbook.
  - NRAO 1971: 0.97 ± 0.08, from Sramek's 1972 Varenna lecture.

  Both are plotted in Will's figure, so I added them as `secondary_only`.

### Lebach et al. 1995: an error in Will 2014

- The paper gives γ = 0.9996 ± 0.0017 ("estimated standard error"). That equals a deflection ratio of 0.9998 ± 0.0008, or γ−1 = (−4 ± 17)×10⁻⁴.
- Will 2014 §3.4.1, and the draft, give γ−1 = (−8 ± 34)×10⁻⁴.
- Will's numbers are exactly twice the published ones. That is what happens if 0.9996 ± 0.0017 is read as (1+γ)/2.
- Fomalont et al. (2009) and Froeschlé et al. (1997) quote the paper correctly.

### Shapiro et al. 2004

- Published: γ = 0.99983 ± 0.00045. The abstract prints the last digits as subscripts (0.9998₃ ± 0.0004₅). This is an "estimated standard error".
- The formal error is 0.00026, according to Lambert & Le Poncin-Lafitte 2009.
- Will converts it correctly: (1+γ)/2 = 0.99992 ± 0.00023.

### Lambert & Le Poncin-Lafitte 2009 and 2011

- The 2009 abstract gives no single value. It says only that γ "cannot be estimated at better than 2 × 10⁻⁴".
- 0.99984 ± 0.00015 (Will's γ−1 = (−1.6 ± 1.5)×10⁻⁴) is one row of its Table 1. Other time spans in that table give anything from 0.99940 to 1.00031.
- The 2011 value, 0.99992 ± 0.00012, also has a formal error. Going by the 2009 paper's conclusion, a realistic error is nearer 2×10⁻⁴.
- The 2004, 2009 and 2011 VLBI analyses largely reuse the same data, so they are not independent.

### Titov et al. 2018

- A single dedicated session in 2017 gives γ−1 ≈ (0.5 ± 1.1)×10⁻⁴ (formal).
- Their full global solution gives γ−1 = (2.72 ± 0.92)×10⁻⁴, 3σ away from GR. They attribute this to systematics at the 2–3σ level.
- I found no later VLBI value of γ.

### Time delay

- **Shapiro et al. 1971:** λ = (1+γ)/2 = 1.02 ± 0.05; the formal error is 0.02. Will plots this point at 1.00.
- **Viking (Reasenberg et al. 1979):** published as γ = 1.000 ± 0.002, which is about twice the formal σ. This equals the draft's (1+γ)/2 = 1.000 ± 0.001, but it is not what the paper printed.
- **Voyager 2:** γ = 1 ± 0.03 (realistic); the formal error is 0.01.
- **Cassini:** γ−1 = (2.1 ± 2.3)×10⁻⁵.
  - The Nature abstract does not state the error convention, and I could not get the full text.
  - The same group (Bertotti, Ashby & Iess 2008) call 2.3×10⁻⁵ "the formal accuracy", and Kopeikin et al. (2007) call it "the quoted standard deviation". So it is a formal 1σ error.
- **Shapiro et al. 1968, the first detection:** the letter has no abstract and is paywalled. Will's error bar and Pössel (2021) both point to γ = 0.8 ± 0.4, i.e. (1+γ)/2 = 0.9 ± 0.2.

### Planetary ephemerides

These are a different kind of "measurement" from the rest of the table.

- **Intervals are not σ.** Most INPOP numbers are intervals, not σ. One paper can give intervals that differ by a factor of about 14, depending on the threshold chosen. For example, Verma et al. 2014 give ±2.5, ±1.1 or ±0.18 × 10⁻⁵.
- **Least-squares errors are admitted to be "very optimistic".** Verma's least-squares fit gives (4.53 ± 0.54)×10⁻⁵, which is nominally 8σ from GR.
- **Pitjeva & Pitjev 2013 contradict themselves on the error level.** Equation (3) labels γ−1 = (4 ± 6)×10⁻⁵ as "(3σ)", but the Conclusion labels the same numbers "(σ)".
- **A probable transcription error.** Fienga et al. 2015 (Table 2 of the arXiv version) and the 2024 Living Review (Table 5) list INPOP10a as γ−1 = (−6.2 ± 8.1)×10⁻⁵. That is in fact INPOP10a's β−1, from the run with γ held fixed. INPOP10a's actual γ−1 is (4.5 ± 7.5)×10⁻⁵ (Fienga et al. 2011, Table 4; Verma et al. 2014; Pitjeva & Pitjev 2013).
- **The MESSENGER papers did not measure γ.** Park et al. 2017 (AJ 153, 121) and Genova et al. 2018 (Nat. Commun. 9, 289) both took the Cassini value. Genova et al.: "γ was thus fixed to 1, while still considering the uncertainty obtained by Cassini."
- **The Living Review's (γ−1) = (0.55 ± 1.35)×10⁻⁵** is the overlap of published intervals, not a measurement.
- **Dilaton limits are not PPN measurements.** Bernus et al.'s limits from the dilaton model (|γ−1| of order 10⁻⁸, listed in Fienga et al. 2022) depend on that model.

## 4. Searched for and not found

- **Gaia:** no measurement of γ has been published as of October 2026. Published Gaia work on γ is simulations and forecasts. For example, Butkevich et al. 2022 (arXiv:2204.02789) show how an error in γ would bias the parallaxes.
- **BepiColombo (MORE radio-science experiment):** only simulations so far (di Stefano et al., arXiv:2201.05107; Cappuccio et al., arXiv:2201.05092). I found no published measured γ.
- **VLBI after Titov et al. 2018:** nothing new. The 2019 conference paper (arXiv:1908.00973) repeats the 2018 results. A 2026 forecast for SKA-VLBI (arXiv:2606.25549) is not a measurement.
- **Sources I could not read:**
  - the full text of Shapiro et al. 1968
  - Will's book *Theory and Experiment in Gravitational Physics* (TEGP), §7.1–7.2 and its tables
  - the original Mikhailov and Matukuma papers on the 1936 eclipse
  - Sramek's 1972 Varenna lecture
  - the primary paper for Mariner 9; probably Anderson et al. 1978, Acta Astronautica 5, 43
  - the primary paper for PSR 1937+21
  - the Nature full texts of Robertson & Carter 1984, Robertson et al. 1991 and Bertotti et al. 2003 (I read their abstracts)
  - the published journal versions of Fienga et al. 2015 and Bruns 2018 (I read the arXiv versions)

## 5. What Will's Figure 5 actually plots

I read the figure from its vector file (`data/literature/will2014/x/PLOTS/gammabound.eps`), parsing the Illustrator path coordinates and calibrating them on the axis tick marks. The values below are where the points are drawn, as (1+γ)/2.

| drawn at (year) | drawn value | matching row |
|---|---|---|
| 1919 | 0.92, bar to 1.09 | E1919-principe (p.e. used as the bar) |
| 1919 | arrow, off scale | E1919-sobral-4inch (1.13) |
| 1922 | 0.98, bar 0.92–1.04 | E1922-wallal-lick-15ft |
| 1922 | 1.04, bar 0.95 to off scale | E1922-wallal-lick-5ft |
| 1929, 1936, ~1950 | arrows | Takengon 1929; USSR 1936; Bocaiuva 1947 (drawn at about 1950) |
| 1952 | 0.97 ± 0.06 | E1952-khartoum |
| 1973 | 0.95, bar to 1.06 (lower end clipped) | E1973-chinguetti |
| 1969 | 1.04; 1.01 | R1969-goldstone; R1969-owens-valley |
| 1970 | 1.07; 1.03; 0.90 | R1970-cambridge-hill; R1970-goldstone-haystack-vlbi; R1970-nrao-sramek |
| 1971 | 0.97 | R1971-nrao-sramek |
| 1972 | 1.04; 0.99; 0.96 | R1972-cambridge-riley; R1972-vlbi-counselman; R1972-westerbork |
| 1973 | 1.04 | R1973-westerbork-dual |
| 1974 | 1.01 | R1974-nrao-fomalont |
| 1975 | 1.000 ± 0.012 | R1975-nrao-fomalont (1975 alone, 0.999 ± 0.011) |
| 1984 | 1.000 ± ~0.005 | V1984-robertson-carter (the paper gives 1.004 ± 0.0025) |
| 1990 | arrows, bars at 0.998/1.002 | V1991-robertson-carter-dillinger (paper: ±0.001) |
| 1996 | arrows | V1987-lebach |
| 2004 | arrows labelled "2×10⁻⁴" | V2004-shapiro |
| ~2008 | arrows labelled "10⁻⁴" | most likely V2009-lambert (γ−1 = (−1.6 ± 1.5)×10⁻⁴ in Will's text); V2005-vlba-fomalont would carry the same label |
| ~1998 | arrows, bars 0.996/1.001 | A1997-hipparcos |
| 2006 | 0.99, bar 0.95–1.03 | L2006-slacs-bolton |
| 1968 | bar only, ~0.94 (clipped) to 1.10 | T1967-shapiro-radar (0.9 ± 0.2) |
| 1971 | 1.00 ± 0.05 | T1970-shapiro-radar (the paper gives 1.02) |
| 1975 | 1.00 ± 0.03 | T1970-mariner6-7 |
| 1977 | 1.00 ± 0.02 | T1972-mariner9 (secondary) |
| 1979 | arrows | T1976-viking |
| 1986 | square, 1.02 ± 0.03 | T1985-psr1937 (secondary) |
| 1992 | 1.000 ± 0.015 | T1985-voyager2 |
| 2003 | arrows, "1×10⁻⁵" | T2002-cassini |

Things to notice about the figure:

- The time axis changes scale at 1960. Before that, one tick is 20 years; after it, one tick is 10 years.
- Some points are placed at the publication year rather than the year of the experiment: Bocaiuva 1947 at about 1950, Voyager 2 (1985) at about 1992, Lebach (1987) at 1996.
- The early optical error bars are probable errors, not σ.
- A few drawn values differ from the papers:
  - Robertson & Carter 1984 is drawn at 1.000; the paper gives 1.004.
  - The 1971 radar point is drawn at 1.00; the paper gives 1.02.
- Not shown at all: Lambert 2011, the Harvey 1979 re-reductions, Chant & Young, Dodwell & Davidson, Matukuma, Bruns 2017, and all ephemeris results. Fomalont et al. 2009 (VLBA) has no point of its own unless it is the ~2008 arrow pair.

## 6. Recommendation for the scorecard plot

**General approach.**

- Plot (1+γ)/2 (GR = 1) against the year of the experiment, with every error converted to 1σ (multiply a probable error by 1.483). Mark the converted points, for example with open symbols.
- Two panels in Will's style work well (deflection; time delay). So does a single panel with a log-scale inset of σ[(1+γ)/2] against year. The inset shows an improvement of more than four orders of magnitude, from ±0.1–0.25 in 1919 to ±1.2×10⁻⁵ with Cassini.
- Draw the points that fall off the scale (1929, 1936, 1947) as arrows, as Will does, or extend the y-axis to about 1.6.

**Recommended set: one point per independent experiment.** Values are (1+γ)/2 ± 1σ, derived from the CSV.

- **Optical eclipses:**
  - 1919 Sobral 4-inch: 1.13 ± 0.10 (from p.e.)
  - 1919 Principe: 0.92 ± 0.25 (from p.e.)
  - 1922 Lick 15-ft: 0.98 ± 0.09 (from p.e.)
  - 1922 Lick 5-ft: 1.04 ± 0.13 (from p.e.)
  - 1922 Cordillo Downs: 1.01 ± 0.25 (from p.e.)
  - 1929 Takengon: 1.28 ± 0.06
  - 1936 USSR: 1.53 ± 0.31 (1959 re-reduction, from p.e.)
  - 1947 Bocaiuva: 1.15 ± 0.15
  - 1952 Khartoum: 0.97 ± 0.06
  - 1973 Chinguetti: 0.95 ± 0.11
  - 2017 Wyoming: 1.00 ± 0.03 (confidence level not stated)
- **Radio interferometry:**
  - 1969: 1.01 ± 0.11; 1.04 (+0.15/−0.10)
  - 1970: 1.03 ± 0.11 (secondary); 0.90 ± 0.05; 1.07 ± 0.17
  - 1971: 0.97 ± 0.08 (secondary)
  - 1972: 1.04 ± 0.08; 0.96 ± 0.05; 0.99 ± 0.03
  - 1973: 1.038 ± 0.033
  - 1974: 1.015 ± 0.011
  - 1975: 0.999 ± 0.011
- **VLBI** (convert γ to (1+γ)/2 by halving both the deviation from 1 and the error):
  - 1984: 1.004 ± 0.0025 (formal)
  - 1987: 0.9998 ± 0.0009
  - 1991: 1.0001 ± 0.0010
  - 2004: 0.99992 ± 0.00023
  - 2005 (VLBA): 0.99990 ± 0.00015
  - 2011: 0.99996 ± 0.00006 (formal; the authors' 2009 conclusion suggests a realistic error of about 1×10⁻⁴)
- **Space astrometry:** Hipparcos 0.9985 ± 0.0015.
- **Time delay:**
  - 1967 radar: 0.9 ± 0.2 (secondary)
  - 1970 radar: 1.02 ± 0.05
  - 1970 Mariner 6/7: 1.00 ± 0.03
  - 1972 Mariner 9: 1.00 ± 0.02 (secondary)
  - 1976 Viking: 1.000 ± 0.001
  - 1985 PSR 1937+21: 1.02 ± 0.03 (secondary)
  - 1985 Voyager 2: 1.000 ± 0.015
  - 2002 Cassini: 1 + (1.05 ± 1.15)×10⁻⁵

**Plot, but flag:**

- The Harvey 1979 re-reductions of 1919. Use open symbols at 1919, or join each to its 1919 original. They are the clearest way to show the controversy.
- The `secondary_only` rows. Use grey or dashed symbols.

**Leave out, or only annotate:**

- Sobral astrographic 0.93: no error was given, and Harvey supersedes it.
- Chant & Young: a range of solutions, no error.
- Mikhailov's original 2.70: no error.
- Matukuma: secondary source only, no error.
- Rows that duplicate a measurement:
  - Jones 1976 (same as the Texas row)
  - the Lick combined value
  - the Fomalont & Sramek 1974+75 mean
  - Titov's global solution
  - Lambert 2009, if Lambert 2011 is plotted
  - Shapiro 2004, if you want strictly independent data; it shares data with Lambert 2011
- Ephemeris results. Their statistics differ (intervals, not σ). Show at most one as an illustration, for example INPOP20a's 99.7% interval, and never the Living Review's interval overlap.
- Bolton et al. 2006. It is not a solar-system test; include it only if you are reproducing Will's figure.
