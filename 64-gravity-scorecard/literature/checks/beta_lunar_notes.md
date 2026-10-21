# β, lunar laser ranging and Mercury: notes on the literature check

Companion to `beta_lunar.csv` (124 rows). It covers the PPN parameter β, the Nordtvedt parameter η_N = 4β − γ − 3 from lunar laser ranging (LLR), Mercury's perihelion, the other LLR tests (geodetic precession, Ġ/G), APOLLO, and Will's 2014 table of PPN limits.

Every number was read in the source itself: the paper's text or tables, the publisher's abstract page, or Will's TeX. Search-engine summaries were not used for any number. The sources are saved under `data/literature/papers/<name>/` (list at the end). Rows marked `secondary_only` name the secondary source in the `reference` and `location` columns.

## How to read the CSV

* **value / err_plus / err_minus** are plain numbers in absolute units. A table entry printed as "(−0.41 ± 0.78)" in a "×10⁴" column becomes −4.1e-5 ± 7.8e-5. The `quote` column keeps the published wording.
  * Where a paper publishes only an interval, the value is the interval's midpoint and the errors are its half-width. The notes column says so. This applies to Pitjeva 2013 and 2021 (Ġ/G).
  * Where a paper gives no central value (Chandler 1989; Shapiro 1988), the value is left blank and only the error is filled.
* **err_type = realistic** means the authors deliberately inflated the least-squares error. The factor is written in the notes column:
  * **IfE Hannover (Müller group):** 3 × formal. The exception is Hofmann & Müller 2018's Ġ/G, which is 4 × formal.
  * **IAA RAS / EPM:** the authors write "3σ", meaning 3 × the formal least-squares σ.
  * **INPOP17a** (Viswanathan 2018): labelled "3-σ".
  * **INPOP20a:** 99.7 % acceptance intervals.
  * **Ke et al. 2026:** 3 × formal combined in quadrature with 0.5 mm.

  The allowed vocabulary has no "3sigma" entry, so these cases are tagged `realistic` and explained in the notes. **Do not** compare them one-to-one with JPL numbers. JPL never states its convention (`not_stated`) and is believed to quote ~1σ-type "realistic" errors.
* **status**:
  * `verified`: a draft number or a Will-table transcription that is confirmed.
  * `corrected`: the draft number or attribution was wrong.
  * `added`: a new primary result.
  * `secondary_only`: the primary could not be read.

## A. Will (2014) Table 4: what each row really is

The transcription in the draft matches `article.tex` exactly. Will gives **no confidence levels** for any entry, and the rows are not like-for-like:

| Row | Limit | Primary behind it | What the primary says | Comment |
|---|---|---|---|---|
| γ−1, time delay | 2.3e-5 | Bertotti, Iess & Tortora 2003 | γ − 1 = (2.1 ± 2.3)e-5 | The limit is the error bar, not value + error. Confidence level is not in the abstract. |
| γ−1, VLBI | 2e-4 | Will's text cites Lambert & Le Poncin-Lafitte 2009: (−1.6 ± 1.5)e-4 | — | 2e-4 = \|value\| + error of their 2011 result (−0.8 ± 1.2)e-4, which Will does not cite. Primaries not read (γ is another session's topic). |
| β−1, perihelion | 8e-5 | **INPOP10a** (Fienga et al. 2011), with γ fixed to Cassini: (−4.1 ± 7.8)e-5 | Interval over which the post-fit residuals change by < 5 %; not a σ | Will's text calls these "data from the Messenger spacecraft that orbited Mercury", but INPOP10a used only the 3 MESSENGER **flybys**. The "J₂ = (2.2 ± 0.1)e-7" remark is helioseismology; INPOP10a itself had J₂ = 2.40e-7 and did not refit it in the PPN runs. |
| β−1, Nordtvedt | 2.3e-4 | Williams, Turyshev & Boggs 2004: β − 1 = (1.2 ± 1.1)e-4 | 1.2 + 1.1 = 2.3 | Value + error, whereas the perihelion row is the error alone. |
| ξ | 4e-9 | Shao & Wex 2013 | \|ξ̂\| < 3.9e-9, 95 % CL | Strong-field analogue. |
| α₁, LLR | 1e-4 | Not cited in Will 2014. Will 2006 cites Müller, Nordtvedt & Vokrouhlický 1996 | α₁ = (−8 ± 9)e-5 "realistic" [(−8 ± 4)e-5 at 90 % C.L.] | Later LLR: (−7 ± 9)e-5 (Müller et al. 2008). |
| α₁, PSR J1738+0333 | 7e-5 | Shao & Wex 2012 | α̂₁ = −0.4 (+3.7/−3.1)e-5, 95 % CL | **Conflict.** 7e-5 is the full width of the 95 % interval; the real bound is \|α̂₁\| ≲ 3.5e-5. |
| α₂ | 2e-9 | Shao et al. 2013 | \|α̂₂\| < 1.6e-9, 95 % CL | Rounded up. |
| α₃ | 4e-20 | Stairs et al. 2005 (method of Bell & Damour 1996) | \|α̂₃\| < 4.0e-20, 95 % | |
| ζ₁ | 2e-2 | none (TEGP combination) | — | `secondary_only`. |
| ζ₂ | 4e-5 | Will 1992, ApJ 393, L59 | \|ζ₂\| < 4e-5, built on the 95 % bound on the second derivative of the pulse period, P̈ | |
| ζ₃ | 1e-8 | Bartlett & van Buren 1986 (active/passive mass of Fe vs Al < 4e-12), converted in TEGP | — | Singh et al. 2023 PRL improved the underlying limit to 3.9e-14 (~100×) but give no ζ₃ number. |

Will 2014 Table 5 (Ġ/G) is also transcribed correctly. The LLR row (4 ± 9)e-13 /yr is Williams et al. 2004. The "Mars ephemeris 0.1 ± 1.6" row (Konopliv et al. 2011) is, according to two other ephemeris papers, a bound on d(GM☉)/dt / GM☉ and not on Ġ/G. Will's text concedes this.

## B. Mercury's perihelion over time

| Year | Who | Number (as published) | Notes |
|---|---|---|---|
| 1859 | Le Verrier, Comptes rendus 49, 380 | Excess of **38″ per century** ("augmenter de 38 secondes le mouvement séculaire du périhélie") | No error given. The often-quoted 38.3″ is from the full memoir, which was not checked. |
| 1882 | Newcomb, Astron. Papers 1, p. 473 | p = **+42.95″**, "greater by 43″" | No error given at that point. |
| 1895 | Newcomb, *Elements of the Four Inner Planets*, p. 109 | e·Δϖ̇ = **+8.48 ± 0.43″/cy** (mean error) | **Corrected.** Dividing by e = 0.2056 gives **41.24 ± 2.09″/cy** (my conversion), not ~43″. Newcomb used a "provisional value 40″.7". |
| 1947 | Clemence, RMP 19, 361 | 42.56 ± 0.94″/cy | `secondary_only` (via Pireaux & Rozelot 2003). Clemence's GR prediction of 43.03″ rested on non-standard AU and c values (Nobili & Will 1986). |
| 1972 | Shapiro et al., PRL 28, 1594 | (2 + 2γ − β)/3 ≃ **1.005 ± 0.007** (statistical) | The authors put the realistic uncertainty at **≈ 0.02**. J₂ was assumed negligible. With γ = 1 this is β − 1 ≈ −0.015 ± 0.021 (stat) or ± 0.06 (total), by my conversion. |
| "1976" | — | The draft's "Shapiro et al. 1976 radar" was **not found** | The 1976 Shapiro paper is the LLR paper (Shapiro, Counselman & King). It gives β = 1.003 ± 0.005 from combining LLR with other results "as if independent". The compilation value "(2+2γ−β)/3 = 1.003 ± 0.005" looks like that β mislabelled. Do not use. |
| 1986 | Nobili & Will, Nature 320, 39 | GR prediction **42.98″/cy**; radar agrees with GR to 0.5 % (±0.2″/cy) | This is the primary for the 42.98″ that Will uses (Will 2014 Eq. E34). |
| 1990 | Shapiro (GR12 proceedings) | \|2γ − β − 1\| < 3e-3 | `secondary_only` (Will 2006). This was Will's perihelion β row until 2006. |
| 2017 | Park et al., AJ 153, 121, Table 3 | GR (gravitoelectric) term **42.9799 ± 0.0009″/cy**; total precession 575.3100 ± 0.0015″/cy | Lense–Thirring −0.0020″, solar J₂ +0.0286″. |

Main caveat: the perihelion advance alone measures only the combination (2 + 2γ − β)/3 + 3e-4 (J₂/1e-7). Every "β from Mercury" result therefore depends on how γ and J₂ were handled. That is why the `fixed_parameters` column matters.

## C. β determinations: what was fixed

| Result | β − 1 (×1e-5) | γ | J₂ | Error type |
|---|---|---|---|---|
| INPOP10a, 2011 | −4.1 ± 7.8 (−6.2 ± 8.1 with γ = 1) | fixed: Cassini (or 1) | fitted once (2.40e-7), not refit in the PPN runs | < 5 % residual change |
| EPM2011, Pitjeva & Pitjev 2013 | −2 ± 3 (γ − 1 = 4 ± 6) | fitted together | fitted (2.0 ± 0.2)e-7 | Eq. (3) says "3σ" but the Conclusions say "(σ)"; the companion SSR paper confirms 3 × formal |
| Verma et al. 2014 (INPOP13a + MESSENGER orbit) | 0.2 ± 2.5 at 25 % (−0.15 ± 0.70 at 10 %; 0.02 ± 0.12 at 5 %) | varied together | refit | Residual-change thresholds. The authors call their least-squares 1σ (1.34 ± 0.043) "very optimistic". |
| Fienga et al. 2015 | arXiv v1: −0.25 ± 6.70 (1σ of a Monte Carlo sample). Published: −6.7 ± 6.9 (LS "3-σ", secondary via Fienga 2021) | varied/fitted | varied | Two different versions exist. Use the published one with caution. |
| Park et al. 2017 (MESSENGER orbit, JPL) | **−2.7 ± 3.9** (abstract) but −2.6 ± 3.9 in Table 1 | considered at Cassini | fitted (2.25 ± 0.09)e-7 | Not stated |
| Genova et al. 2018 (MESSENGER) | **−1.6 ± 1.8** | fixed = 1 (Cassini error considered) | fitted (2.246 ± 0.022)e-7 | Formal. **Imposes** η = 4β − γ − 3, so with γ = 1, η = 4(β − 1): not independent of η. Without the constraint β − 1 = (−14.3 ± 14.7)e-5. |
| INPOP20a, Fienga et al. 2021 | −1.12 ± 7.16 (WRSS); −1.9 ± 6.28 (LS) | varied (−1.69 ± 7.49) | refit | 99.7 % (3σ) |
| LLR, Williams et al. 2004 | 12 ± 11 | Cassini | — | Via η (WEP subtracted) |
| LLR, Müller et al. 2008 | 15 ± 18 | Cassini | — | Via η, WEP assumed; 3 × formal |
| LLR, Hofmann & Müller 2018 | **−4.5 ± 5.6** via η + Cassini; **−8.7 ± 9.0** direct EIH fit (thesis) | Cassini / fitted | — | 3 × formal. Biskupek 2021 quotes the second number as "Hofmann & Müller", the CQG abstract the first. |
| LLR, Biskupek et al. 2021 | 6.2 ± 7.2 (direct EIH fit) | fitted (17 ± 16) | — | Not stated (IfE practice is 3 × formal) |
| Review, Fienga & Minazzoli 2024 | −0.45 ± 1.75 | — | — | Overlap of intervals, not a measurement (`secondary_only`) |

No separate INPOP19a β result was found: the INPOP19a data set was used for INPOP20a. EPM2021 and Pitjeva & Pitjev 2014 (CeMDA 119, 237) β values were not located or not accessible. Pavlov & Dolgakov 2024 (Universe 10, 413) is about the Sun's pole, J₂ and μ̇ and gives no β in its abstract. BepiColombo reaches Mercury orbit only at the end of 2026, so no MORE results exist yet; INPOP20a + BepiColombo *simulations* project ±5e-5 (WRSS) or ±1.1e-5 (least squares) for β − 1. These are forecasts, not data.

## D. The Nordtvedt effect: WEP vs SEP

LLR measures the **total** Earth–Moon difference Δ(m_G/m_I)_EM = Δ_WEP (composition: iron core vs silicate Moon) + η × [(U/Mc²)_E − (U/Mc²)_M]. The self-energy factor is −4.45e-10 (WTB 2009 Eq. 7). The two pieces are separated by the Eöt-Wash torsion-balance test with test bodies that mimic Earth-core and Moon-mantle compositions:

* **Baessler et al. 1999:** (+0.1 ± 2.7 ± 1.7)e-13.
* **Adelberger 2001:** (1.0 ± 1.4 ± 0.2)e-13. Everyone uses the combined (1.0 ± 1.4)e-13. This row is `secondary_only` because it was read via WTB 2009.

So η is published in two different ways. Do not mix them on one plot without saying so.

1. **"Clean" SEP (lab WEP subtracted).**
   * WTB 2004 / 2009: η = (4.4 ± 4.5)e-4.
   * Hofmann 2017 thesis: η_SEP = (2.9 ± 3.4)e-4.
   * Viswanathan 2018: η_SEP = (3.1 ± 3.6)e-4.

   These have stalled at an error of ~3.5–4.5e-4 since 2004, because the **lab WEP error (1.4e-13) now dominates**.
2. **WEP assumed** (the whole signal attributed to SEP).
   * Hofmann & Müller 2018: η = (−0.2 ± 1.1)e-4.
   * INPOP17a: η_metric = (0.85 ± 1.59)e-4.
   * Ke et al. 2026 (EPM21): η = (−0.20 ± 0.72)e-4.
   * Biskupek 2021 gives only Δ = (−2.1 ± 2.4)e-14 and no η. Ke et al. convert it to η ≈ (0.47 ± 0.54)e-4.

   Theorists justify this route for *metric* theories, which satisfy the WEP automatically.

The LLR-only total Δ improved roughly 6×:

| Year | Analysis | Δ (×1e-13) |
|---|---|---|
| 2004 | WTB | −1.0 ± 1.4 |
| 2012 | WTB | −0.8 ± 1.3 |
| 2018 | HM | −0.30 ± 0.50 |
| 2018 | INPOP17a | −0.38 ± 0.71 |
| 2021 | Biskupek | −0.21 ± 0.24 |
| 2026 | Ke | 0.09 ± 0.32 |

The JPL numbers (2004, 2012) have no stated error convention; the others are about 3 × formal.

Converting η to β: β − 1 = [η + (γ − 1)]/4, with Cassini γ. This assumes ξ = α₁ = α₂ = ζ₁ = ζ₂ = 0, because the full expression is η_N = 4β − γ − 3 − 10ξ/3 − α₁ + 2α₂/3 − 2ζ₁/3 − ζ₂/3.

**Genova 2018's η = (−6.6 ± 7.2)e-5** is a different test. It comes from the Sun's self-energy shifting the solar-system barycentre, seen in Mercury's orbit; it is not the Earth–Moon test.

**WTB 2012 CQG gives only the EP number (−0.8 ± 1.3)e-13.** It has no η, SEP or β; the draft implied otherwise. Genova 2018 Table 1 attributes β − 1 = (1.2 ± 1.1)e-4 to it, which is a misattribution: that value is WTB 2004/2009.

## E. Size of the Nordtvedt signal (coefficient)

| Source | Coefficient | Status |
|---|---|---|
| Nordtvedt 1968, Phys. Rev. 170, 1186 (abstract) | "range oscillation of amplitude about **12 m**" if m_g/m_i − 1 ≈ the Earth's self-energy fraction | Abstract only. Damour & Vokrouhlický (D&V) say their full result is "approximately 60 % larger than the lowest-order estimate first derived by Nordtvedt in 1968". That makes the first-order value ≈ 13.10/1.6 ≈ 8.2 m (my arithmetic). Why the abstract says 12 m could not be checked. |
| "Standard estimate" until ~1995 (Nordtvedt 1973; Will TEGP 1981) | ≈ **9.3 η m** | `secondary_only`, via D&V 1996 |
| Nordtvedt 1995, Icarus 114, 51 | S ≈ −2.9e13 mm. WTB 2009 says 12.8 m; WTB 2012 Table 3 says 12.90 m | Conflicting secondaries. D&V could not extract a precise full value from that paper. |
| Damour & Vokrouhlický 1996 | **2.9427e12 δ̄₁₂ cm** (Moon-minus-Earth convention, = S of −2.9427e13 mm in JPL's Earth-minus-Moon form) and **13.10 η m** | Verified. First order × 1.62201 (tidal amplification). |
| Nordtvedt & Vokrouhlický 1997 | −2.943e13 mm, 13.10 m | `secondary_only` (WTB 2012) |
| JPL numerical, WTB 2009 / 2012 | −2.92e13 to −2.99e13 mm → 13.0–13.3 m; −2.992e13 mm → 13.31 m | Verified. Agrees with analytic theory to a few %. |
| Nordtvedt 2003 review | 2.9e12 δ cos D cm; feedback factor ≈ 1.75 | |
| Will 2001 and 2014 (Eq. E36) | **13.1 η_N m** | Verified. WTB 2004 round it to "13 meter". |

Self-energies (WTB 2009 Eqs. 5–7, originally from Williams, Newhall & Dickey 1996):

* **Earth: U/Mc² = −4.64e-10.** Two Earth models agree to 0.1 %; a uniform Earth would be 10 % smaller.
* **Moon: −1.90e-11.** Within 1 % of a uniform Moon. Williams et al. 1996 and D&V used −0.19e-10.
* **Difference: −4.45e-10.**
* Will 2014 rounds these to 4.6e-10 and 0.2e-10.

The draft numbers are all confirmed.

## F. Geodetic (de Sitter) precession of the lunar orbit

The GR value is **19.2 mas/yr** (≈ 1.92″/century). This is confirmed in Bertotti, Ciufolini & Bender 1987, Dickey et al. 1989 and Williams, Newhall & Dickey 1996.

| Year | Analysis | Result |
|---|---|---|
| 1987 | Bertotti et al. | ~10 % |
| 1988 | Shapiro et al. | σ = 0.04″/cy = 2 %, systematics included |
| 1989 | Dickey et al. | 2 % |
| 1996 | WND | "within 0.7 %" |
| 2004 | WTB | **K_gp = −0.0019 ± 0.0064** |
| 2008 | Müller et al. | (6 ± 10)e-3 ″/cy |
| 2018 | Hofmann & Müller | "confirmed within 0.09 %"; thesis h = (−5.6 ± 8.5)e-4, 3 × formal |

**Correction:** the draft put K_gp = −0.0019 ± 0.0064 under Williams, Newhall & Dickey 1996. It belongs to WTB 2004. Will 2014 summarises the LLR result as "about 0.6 percent".

## G. Ġ/G

| Source | Ġ/G (per yr) | Error convention / note |
|---|---|---|
| WND 1996 | \|Ġ/G\| ≤ 8e-12 | upper limit |
| WTB 2004 | (4 ± 9)e-13 | Will's Table 5 LLR row |
| Müller 2008 | (6 ± 8)e-13 | 3 × formal |
| Pitjeva & Pitjev 2013 | −7.0e-14 to +7.8e-14 | 95 % |
| Fienga 2015, arXiv v1 | (0.09 ± 0.91)e-13 | |
| Hofmann & Müller 2018 | **(7.1 ± 7.6)e-14** | **4** × formal; core rotation fixed, station coordinates constrained |
| Genova 2018 | < 4e-14 | |
| Pitjeva et al. 2021 (EPM2019) | −2.9e-14 to +4.6e-14 | 3σ |
| Biskupek 2021 | **(−5.0 ± 9.6)e-15**; (0.2 ± 1.3)e-14 when fitted jointly with G̈ | |

Planetary-ephemeris values are really measurements of d(GM☉)/dt. Turning them into Ġ/G requires subtracting a modelled solar mass loss of ≈ −(0.9 to 1.1)e-13 /yr. The result is therefore model-dependent; LLR is less affected by this.

## H. APOLLO and the precision of LLR

**Precision vs accuracy vs model fit.** APOLLO normal points:

* **Median uncertainty: 3.4 mm** raw, **4.2 mm** including systematic timing errors (Battat et al. 2023).
* **Median nightly uncertainty: 1.7 mm** (1.3 mm without the systematic term), and **1.0 mm since Sept 2016**.
* The **absolute calibration system** revealed a ~0.4 % hardware calibration error. After correction, systematics are **< 1 mm** (Colmenares et al. 2023).
* The early "one-millimeter range precision" claim (Murphy et al. 2008) refers to statistical precision for well-populated normal points.

**Model residuals are an order of magnitude larger:**

* JPL: ~18 mm weighted RMS for APOLLO and the Côte d'Azur station (OCA) (Murphy 2013).
* JPL 2005–2011, all stations: 1.8 cm.
* INPOP17a: **1.15–1.95 cm** for the four APOLLO groups.
* **DE440: "about 1.3 cm" for recent LLR data overall.** DE440 gives no APOLLO-only number, so the draft's request cannot be met from DE440 itself. DE440 used 2452 APOLLO normal points from 2006–2017.

**Time span:** APOLLO archive 2006-04-06 to 2020-12-27 (4296 normal points, 35 % of all LLR normal points in that period).

## I. Which entries to use for plots

**β timeline.** Use one marker style per method and keep the published error conventions. Do **not** rescale the "realistic" errors to 1σ: they are deliberately inflated for systematics.

* Radar/perihelion:
  * 1972 Shapiro et al. (convert from (2 + 2γ − β)/3; use the ±0.02 total error).
  * 1989 Chandler: σ_β = 0.02 with J₂ fixed, no central value.
  * 1990 Shapiro (`secondary_only`).
* LLR/Nordtvedt (+ Cassini γ):
  * 1976 Williams et al.: \|β − 1\| < 0.01.
  * 1976 Shapiro, Counselman & King: β − 1 = 0.003 ± 0.005, combined.
  * 1996 WND: < 1.4e-3 (6e-4 if WEP assumed).
  * 2004 WTB: (1.2 ± 1.1)e-4.
  * 2008 Müller: (1.5 ± 1.8)e-4.
  * 2018 Hofmann & Müller: (−4.5 ± 5.6)e-5.
* LLR direct (EIH):
  * 2018 Hofmann & Müller: (−8.7 ± 9.0)e-5.
  * 2021 Biskupek: (6.2 ± 7.2)e-5.
* Planetary ephemerides:
  * 2011 INPOP10a: (−4.1 ± 7.8)e-5.
  * 2013 EPM2011: (−2 ± 3)e-5.
  * 2014 Verma: (0.2 ± 2.5)e-5.
  * 2015 Fienga: published (−6.7 ± 6.9)e-5 (`secondary_only`) or arXiv (−0.25 ± 6.70)e-5.
  * 2017 Park: (−2.7 ± 3.9)e-5.
  * 2018 Genova: (−1.6 ± 1.8)e-5, flagged "Nordtvedt-constrained, γ = 1".
  * 2021 INPOP20a: (−1.12 ± 7.16)e-5.

**γ–β plane (Will-style).**

* **Cassini band:** γ − 1 = (2.1 ± 2.3)e-5.
* **Nordtvedt band:** 4(β − 1) − (γ − 1) = η.
  * Use a *clean-SEP* η for a model-independent plot: WTB 2004 (4.4 ± 4.5)e-4, the Will-era value; or Viswanathan 2018 (3.1 ± 3.6)e-4; or Hofmann 2017 (2.9 ± 3.4)e-4.
  * Or a WEP-assumed η, labelled as such: Hofmann & Müller 2018 (−0.2 ± 1.1)e-4; Ke 2026 EPM21 (−0.20 ± 0.72)e-4.
* **Perihelion band:** 2(γ − 1) − (β − 1) = 3[(2 + 2γ − β)/3 − 1].
  * Historic: Shapiro 1972, 3 × (0.005 ± 0.02); or \|2γ − β − 1\| < 3e-3 (Shapiro 1990, `secondary_only`).
  * Modern joint fits are better drawn as boxes, not bands:
    * EPM2011 (β − 1, γ − 1) = (−2 ± 3, 4 ± 6)e-5 (3 × formal).
    * Verma 2014 25 % box (0.2 ± 2.5, −0.3 ± 2.5)e-5.
    * INPOP20a 99.7 % box (−1.12 ± 7.16, −1.69 ± 7.49)e-5. Its least-squares β–γ correlation is ~25 %.
* **Do not** draw Park 2017, Genova 2018, INPOP10a (Cassini γ) or Will's 8e-5 as independent β bands. They all assumed γ, and Genova also assumed the Nordtvedt relation.

## J. Conflicts and traps found

1. Will 2014 Table 4, α₁ (PSR J1738+0333) = 7e-5 vs primary 95 % interval [−3.5, +3.3]e-5.
2. Will's two β rows use different definitions: error only (8e-5) vs value + error (2.3e-4).
3. Pitjeva & Pitjev 2013: "3σ" in Eq. (3) but "(σ)" in the Conclusions. This is also in the published version.
4. Park et al. 2017: −2.7 ± 3.9 in the abstract vs −2.6 ± 3.9 in Table 1.
5. Hofmann & Müller 2018 have two β values: −4.5 ± 5.6 (abstract, via η) and −8.7 ± 9.0 (direct; quoted by Biskupek 2021).
6. Hofmann & Müller 2018's η assumes WEP; their clean-SEP η is (2.9 ± 3.4)e-4.
7. Nordtvedt 1995 coefficient: 12.8 m (WTB 2009) vs 12.90 m (WTB 2012). The 1968 abstract says "about 12 m", but per D&V the first-order formula gives ≈ 8 m (their full result is ~60 % larger).
8. The draft's K_gp attribution (WND 1996 → actually WTB 2004).
9. The draft's "Newcomb 1895 ~43″" → Newcomb 1895 gives 41.24 ± 2.09″/cy after dividing his e·Δϖ̇ by e.
10. "Shapiro et al. 1976 radar" — no such paper; the compilation number is probably β from the LLR paper.
11. Fienga 2015: the arXiv v1 and the published CeMDA paper differ in title and numbers.
12. Transcription slips in Fienga & Minazzoli 2024 LRR Table 5:
    * INPOP10a's β (γ = 1) printed in the γ column.
    * The INPOP20a least-squares β printed as +1.9 instead of −1.9.
    * Park 2017 labelled "DE440".
    * Verma 2014 labelled "INPOP13c".
13. Genova 2018 Table 1 cites WTB 2012 for β − 1 = (1.2 ± 1.1)e-4; it is WTB 2004.
14. Konopliv 2011's "Ġ/G = 0.1 ± 1.6" in Will's Table 5 is a d(GM☉)/dt bound.
15. Ġ/G from Hofmann & Müller 2018 uses 4 × formal, not 3 ×.
16. The published Pireaux & Rozelot 2003 lists two authors; the arXiv preprint adds S. Godier.

## Sources saved (data/literature/papers/)

**arXiv PDFs and text:**

* williams2004prl, williams2004ijmpd, williams2009ijmpd, williams2012cqg
* damour1996prd, nordtvedt2003, muller2008llr
* fienga2011inpop10a, verma2014aa, pitjeva2013mnras, pitjeva2013ssr, fienga2015cemda
* viswanathan2018mnras, biskupek2021universe, fienga2021inpop20a, fienga2022moriond, fienga2024lrr, pitjeva2021aa, mariani2024aa, ke2026prd
* murphy2008pasp, murphy2013rpp, adelberger2017cqg, battat2023apollo, colmenares2023apollo
* shao2012cqg_alpha1, shao2013cqg_alpha2, shao2013cqg_xi, stairs2005apj
* will2001lrr, will2006lrr, pireaux2003apss

**Publisher or repository copies:**

* biskupek2021universe_published, singh2023prl_published (Hannover repository)
* genova2018natcomm (Nature open-access PDF)
* park2017aj, park2021de440 (IOP open-access HTML)
* hofmann2018cqg (IOP abstract page only; the PDF is behind a bot check)
* murphy2012cqg (abstract page)
* hofmann2017thesis (DGK C-797 PDF)
* aps_abstracts (APS abstract pages), muller1996prd_alpha1, bertotti2003nature, nobili1986nature (abstract pages)
* ntrs_abstracts (NASA NTRS records holding the published abstracts)
* chandler1989ntrs (NTRS PDF)

**Historical scans (archive.org and ADS):**

* leverrier1859cr (Comptes rendus scan + OCR)
* newcomb1882apae, newcomb1895elements (scan + OCR; Newcomb's p. 109 table checked on the image)
* will1992apjl (ADS scan)
