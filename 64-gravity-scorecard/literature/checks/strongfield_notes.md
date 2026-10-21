# Strong-field, gyroscope and large-scale GR tests: literature check

Companion to `strongfield.csv` (43 rows) for atlas project N7, the "Is Einstein still right?" scorecard. Checked on 2026-10-19. Every number in the CSV was read in the primary paper: the arXiv PDF or TeX source, the publisher's abstract page, or a free Extended Data table or peer-review file. No number comes from a search snippet or a review. Downloaded sources are in `data/literature/papers/sf-*/`.

## 1. How to read the CSV

- **Values are as published.** `value`, `err_plus` and `err_minus` keep the paper's own convention: ratio, deviation parameter, mas/yr, eV/c², and so on. `err_syst` is filled only when the paper quotes a separate systematic error. `err_type` gives the interval convention.
- **Two-sided bounds** (`gw170817_speed`, `psaltis2020_alpha13`): `value` holds the lower bound and `err_plus` holds the upper bound. In these two rows `err_plus` is a bound, not an error bar.
- **Upper limits** (graviton mass rows): `value` is the 90% upper limit. Plot these as arrows, not as points.
- **Rows that are not deviation measurements:**
  - `lares2_2023_first_results`: no number was published.
  - `gw150914_imr`: a GR quantile in %.
  - `gw250114_arealaw`: a significance in σ.
  - `eht_m87_2019_delta_gas`: illustrates prior dependence.

  The notes column of each says "NOT PLOTTABLE", "NOT A DEVIATION VALUE" or "DO NOT PLOT AS A GR TEST".
- **Quote notation.** Text quotes are verbatim. Superscripts and subscripts that the PDF splits across lines are written as `^{+a}_{−b}`. The Planck rows quote the TeX table row verbatim.
- **Converting to 1σ before plotting:**
  - LVK numbers are medians with symmetric 90% credible intervals. For a Gaussian, 1σ ≈ (90% half-width)/1.645.
  - 95% intervals: divide by 1.96.
  - The 68% EHT, Collett, Planck and DESI intervals are already about 1σ.
  - The EHT Sgr A* and DESI central values are posterior means. Most others are medians or best fits.
- **Statuses:** 18 verified, 4 corrected, 21 added, 0 secondary_only.

## 2. Corrections to the draft

1. **LAGEOS 2004 (Ciufolini & Pavlis, Nature 431, 958).** The published result is 0.99 ± 0.05: "99 ± 5 per cent", including all known errors. The authors then "allow for a total ± 10 per cent uncertainty to include underestimated and unknown sources of error". The draft's ±0.10 is that allowance, not the quoted measurement error. Use ±0.10 if you want a conservative error bar.
2. **GWTC-3 graviton mass.** The draft's m_g ≤ 1.27 × 10⁻²³ eV/c² is arXiv v1 (2021). The paper was finally published as Phys. Rev. D 112, 084080 (2025), and the published and v3 bound is **m_g ≤ 2.42 × 10⁻²³ eV/c² (90%)**. The v3 arXiv comment reads "fixed graviton mass bound in abstract". The same revision changed the quoted GWTC-2 bound from 1.76 to 3.09 × 10⁻²³. It also changed the comparison from "2.5 times better than" to "1.31 times better than the Solar System bound of 3.16 × 10⁻²³ eV/c²".
3. **Psaltis et al. 2020.** δ = −0.01 ± 0.17 is not in this paper. It is in EHT M87 Paper VI (ApJL 875, L6, §9.5). Psaltis et al. take the "∼17% at the 68-percentile level" agreement and turn it into bounds on 2PN-type metric parameters: −3.6 < α₁₃ < 5.9 (Johannsen–Psaltis metric) and −5.0 < γ₁,₂ < 4.9 (MGBK metric).
4. **GW150914 IMR consistency.** The paper quotes no fractional precision. It gives only Fig. 4 and the statement that GR lies on the contour enclosing 28% of the posterior. Its Table I does give a parameterized-inspiral bound, δφ̂₀ = −0.1 ± 0.1 (90%), and the paper gives a graviton bound of 1.2 × 10⁻²² eV/c² (90%). Both are added as rows.

Everything else in the draft checked out, including all the S2 numbers, both EHT Sgr A* δ values, GP-B, LARES 2016/2019, GW170817, Reyes, Collett and Nguyen et al. Some draft entries had no number; those values are now added.

## 3. Notes by group

### 3.1 Gyroscopes: Gravity Probe B
- Everitt et al. 2011 (PRL 106, 221101), Table II. Geodetic precession is −6601.8 ± 18.3 mas/yr against GR's −6606.1, a ratio of 0.99935 ± 0.00277 (computed here). Frame dragging is −37.2 ± 7.2 mas/yr against −39.2, a ratio of 0.949 ± 0.184. Both errors are 1σ totals, statistical plus systematic (Table III).
- The final CQG paper (Everitt et al. 2015, CQG 32, 224001) repeats exactly the same numbers. Cite it for the full analysis, but nothing changed.
- Caveats:
  - Large electrostatic "patch-effect" torques had to be modelled, so the frame-dragging error is far above the ~0.5 mas/yr design aim.
  - Will 2014 summarises the result as 0.3% (geodetic) and 20% (frame dragging).
  - GR's frame-dragging value depends on the guide-star direction (IM Pegasi).

### 3.2 Frame dragging with laser-ranged satellites (LAGEOS, LARES, LARES-2)
- **History.**

  | Year | Result | Source |
  |---|---|---|
  | 2004 | 0.99 ± 0.05, with ±0.10 allowed for unknown errors | LAGEOS + LAGEOS 2 |
  | 2016 | (0.994 ± 0.002) ± 0.05 | |
  | 2019 | (0.9910 ± 0.0006) ± 0.02 (0.03–0.04 if the gravity-field errors are inflated ×2–3) | |
  | 2020 | μ − 1 = (1.5 ± 7.4) × 10⁻³ ± 16 × 10⁻³, with the statistical error at 95% | Lucchesi et al., Universe 6, 139: an independent team (LARASE) using the same three satellites |
  | 2023 | First LARES 2 paper: no number, only a figure | |
  | 2026 | Nature result (below) | |

  - In the 2016 and 2019 Ciufolini results the small error is the formal fit error. The large one is an *estimated* systematic, which is not a formal σ.
  - The 2004 ±5% is said to include "all known random and systematic errors".
  - In 2026 the systematic part was estimated by Monte Carlo and an error budget.
- **The 2026 Nature paper.**
  - Reference: Ciufolini, Paolozzi, Pavlis, Ries et al., *Nature* 655, 332–335, doi:10.1038/s41586-026-10715-0, online 8 July 2026.
  - Result (Extended Data Table 2c): **μ = 1.001**, with precision/statistical error **0.001** and systematic error **0.002**. Total in quadrature ≈ 0.0022 (computed here).
  - The README's "~0.2%" is confirmed. It is the RSS of five 0.1% systematic terms: non-zonal harmonics, tides, albedo, thermal drag and eclipses, and orbital-parameter errors (ED Table 2b). The abstract's "one-part-in-a-thousand" refers to the statistical precision.
  - The authors cross-checked with a Monte Carlo (1.0001 ± 0.0019) and a covariance analysis (1.0007 ± 0.0019). In the peer-review file they say the total uncertainty "may be interpreted as a 1-sigma uncertainty".
  - The experiment sums the nodes of LAGEOS (a = 12,270,020 m, i = 109.85°) and LARES-2 (a = 12,264,567 m, i = 70.16°). Their inclinations are supplementary to 0.009° on average, so the huge Newtonian nodal drift from the Earth's oblateness cancels in the sum, while the ≈30.7 mas/yr frame-dragging drifts add.
  - The fit uses about one 1050-day period of the K₁ tide, the largest residual, whose amplitude is fitted jointly.
  - The main text is paywalled. The numbers come from the free Extended Data tables and the CC-BY peer-review file, both saved in `sf-ciufolini2026-nature/`.
- **Controversy, worth stating on the scorecard.**
  - One Nature referee (rounds 1–2) argued that relativistic perturbations missing from the equations of motion mean "the actual uncertainty in the reported measurement likely exceeds several percent". These were the GCRS/BCRS frame choice, third-body couplings and the relativistic quadrupole. The final-round referees were supportive.
  - L. Iorio has disputed the error budgets at each stage:
    - The 2004 LAGEOS result: Space Sci. Rev. 148, 363, which says the systematic "may be up to 3 to 4 times larger" than the published 5–10%.
    - The 2016 LARES result: EPJC 77, 73 (2017).
    - LARES 2: Universe 9, 211 (2023) argues the summed J₂ precessions do not cancel well enough ("still ≃5000 times larger" than the relativistic rates), and that eccentricities and inclinations would need "3-4 orders of magnitude" better accuracy. EPJC 85, 255 (2025) finds a J₂-driven bias of "nearly 1.4%" when the calibrated GGM05S J₂ error is used, and calls the claimed orbital errors "mere statistical errors of the fitting procedure".
  - Ciufolini et al. replied in EPJC 78, 880 (2018) and EPJC 84, 998 (2024).
  - **Recommendation:** plot the 2026 point as 1.001 ± 0.0022 (1σ, per the authors), labelled "Ciufolini et al. 2026, error budget disputed". Keep these as history:

    | Result | Precision |
    |---|---|
    | 2004 | 5–10% |
    | 2016 | 5% |
    | 2019 | 2% |
    | Independent LARASE 2020 | ≈1.6% systematic, with a 0.74% statistical error at 95% |

### 3.3 Gravitational waves
- **Speed of gravity (GW170817).** −3 × 10⁻¹⁵ ≤ Δv/v ≤ +7 × 10⁻¹⁶.
  - This is a deterministic, *conservative* bound, not a credible interval. It assumes a distance of 26 Mpc (the lower end of the 90% distance interval) and that the gamma rays left 0–10 s after the GW peak.
  - It tests propagation over ~40 Mpc of weak field, not strong-field dynamics. Place it separately on a potential-axis plot.
- **Graviton mass (dispersion).**

  | Data set | Bound (90%) | Source |
  |---|---|---|
  | GW150914 | 1.2 × 10⁻²² eV/c² | |
  | GWTC-3 | 2.42 × 10⁻²³ (published) | |
  | GWTC-4.0 | **1.92 × 10⁻²³ eV/c²** (newest) | arXiv:2603.19020 v2, July 2026, submitted to ApJL |

  - GWTC-4.0 quotes GWTC-3 as 2.23 × 10⁻²³ after re-weighting to a prior flat in m_g, so prior choices shift these bounds by ~10–20%.
  - GWTC-4.0 itself notes that planetary ephemerides give 1.01 × 10⁻²⁴ eV/c² (99.7%), ten times better.
  - The GWTC-5.0 tests-of-GR paper (arXiv:2607.19293, July 2026) does not repeat this test.
- **Inspiral–merger–ringdown (IMR) consistency.** GW250114 alone gives ΔM_f/M̄_f = 0.02 ^{+0.07}_{−0.06} (90%). The GWTC-4.0 hierarchical combination gives 0.00 ^{+0.07}_{−0.06}. GW150914 gave no numerical width.
  - Caveat: the GWTC-5.0 tests-of-GR paper dropped this test. Its reason: Narayan et al. "showed that this test can exhibit spurious biases even for unlensed, quasicircular binary signals in certain portions of the parameter space".
- **Ringdown / black-hole spectroscopy.**
  - GW250114, full-signal pSEOBNR: δf̂₂₂₀ = 0.02 ± 0.02 and δτ̂₂₂₀ = −0.01 ^{+0.10}_{−0.09} (90%). This is "within a few percent" for the dominant frequency.
  - This uses NR-calibrated amplitudes and pre-merger information. The more agnostic ringdown-only tests are weaker. The overtone gives δf₂₂₁ = 0.1 ± 0.3 (the "±30%" of the area-law paper). Post-peak pyRing KerrPostmerger gives δf̂₂₂₀ = 0.09 ^{+0.34}_{−0.22}.
- **The one persistent oddity: the catalog damping time.**
  - Combined over the GWTC-4.0 events, δτ̂₂₂₀ sits off zero: joint 0.17 ± 0.11 (90%), GR quantile 98.6% (bootstrap 98.6 ^{+1.4}_{−9.4}%).
  - In GWTC-5.0 (O1–O4b) it is 0.07 ^{+0.06}_{−0.05} joint, with GR "slightly outside the 90% credible region" (98.1%). Hierarchically it is 0.09 ^{+0.13}_{−0.10} (GR quantile 92.5%).
  - The LVK lists possible causes rather than a diagnosis: non-Gaussian or non-stationary noise, parameter correlations, waveform-model systematics, intrinsic variance from the limited number of events, and selection effects. It says there is "no strong evidence for a GR deviation". The loud GW250114 alone sits on GR.
  - Plot the hierarchical value (LVK's "main result") and mention the joint one.
- **Area theorem (GW250114).** The final horizon area exceeds the initial one at **4.4σ** (t< = −40 t_M). It is ≥3.4σ for all pre-merger cut times > −250 t_M and >5σ for cuts ≥ −10 t_M; the two-mode ringdown version gives 3.6σ.
  - The significance is a Gaussian approximation, (μ_f − μ_i)/√(σ_f² + σ_i²).
  - The test discards the merger itself.
  - The area law assumes GR plus the null energy condition and cosmic censorship.
  - This is a qualitative, pass/fail point, not a deviation parameter.
- **Inspiral PN coefficients.** For the 0PN (quadrupole-formula) term, GW250114 gives δφ̂₀ = 0.00 ± 0.03 (90%). GW150914 gave −0.1 ± 0.1. The double pulsar constrains the same term far better in the weak field.
- **Publication status.** GW250114 area law: PRL 135, 111403 (2025). GW250114 spectroscopy: PRL 136, 041403 (published 29 Jan 2026). GWTC-4.0 tests of GR I–III (arXiv March 2026) and GWTC-5.0 (July 2026) are still preprints, so numbers can change: the GWTC-3 graviton bound changed between arXiv v1 and the journal.

### 3.4 Black-hole shadows (EHT)
- **Two different δ definitions.**
  - M87* Paper VI uses δ ≡ θ_g/θ_dyn − 1: EHT-inferred against dynamically measured mass-to-distance.
  - Sgr A* Paper VI uses δ ≡ d_sh/d_sh,Sch − 1, the inferred shadow diameter against a Schwarzschild shadow of 6√3 θ_g.
  - Sgr A* Paper VI notes that "values of δ in the range [−0.075, 0] are consistent with the Kerr predictions" (spin/inclination). The GR "value" is therefore a band, not exactly 0.
- **Calibration (the "α" issue).** The EHT measures a bright *ring*, not the shadow.
  - For M87*, d = α θ_g, with α ≈ 11.5 calibrated on Kerr GRMHD images (9.6–10.4 for the photon ring itself). The M87* test is therefore a consistency (null-hypothesis) test that assumes Kerr in its own calibration.
  - For Sgr A*, α_c = α₁ × α₂: theory bias times measurement bias. α₁ was calibrated on Kerr GRMHD plus analytic Kerr and *non-Kerr* models. Across imaging methods and calibrations, Table 2 gives δ from −0.17 to −0.02. The fiducial values are eht-imaging with GRMHD.
- **Mass priors.**
  - M87*: −0.01 ± 0.17 (68%) with the stellar-dynamics prior, but 0.78 ± 0.3 with the gas-dynamics prior. The EHT mass, (6.5 ± 0.2 stat ± 0.7 sys) × 10⁹ M☉, agrees with the stellar value only.
  - Sgr A*: −0.08 ± 0.09 (VLTI/GRAVITY prior, θ_g = 5.125 μas) or −0.04 ^{+0.09}_{−0.10} (Keck prior, θ_g = 4.92 μas). The two priors differ by ~4%. Plot one, e.g. VLTI, since both use the same EHT data. The Sgr A* intervals are 68% around the posterior mean.
- **Controversy.** Gralla 2021 (PRD 103, 024023, "Can the EHT M87 results be used to test general relativity?") answers "No": a GR violation cannot be distinguished from failed "gastrophysical" assumptions. The EHT papers frame their results as null-hypothesis tests. Psaltis et al. 2020 argue the shadow size probes the metric at 2PN order, which is inaccessible to Solar-System tests.
- **Newer EHT data.** There is no newer EHT δ measurement.
  - M87* 2018 Paper II (A&A 693, A265, 2025) reports a ring "consistent with the 2017 observation" and focuses on accretion models.
  - The 2017–2021 M87* paper (arXiv:2509.24593) reports "a persistent diameter of 43.9 ± 0.6 µas" but no metric test.
  - Recommendation: use Sgr A* (VLTI prior) as the main shadow point and M87* as the second.

### 3.5 S2 and Sgr A*
- **Redshift.** f (GRAVITY) and Υ (Do et al., Keck/UCLA) scale the *combined* gravitational redshift plus special-relativistic transverse Doppler shift (~200 km/s at pericentre). Neither is a pure gravitational-redshift test.

  | Paper | Value |
  |---|---|
  | GRAVITY 2018 | f = 0.90 ± 0.09 (stat) ± 0.15 (sys) |
  | Do et al. 2019 | Υ = 0.88 ± 0.17, independent data |
  | GRAVITY 2019 | **f = 1.04 ± 0.05** (20σ; newest dedicated value) |
  | GRAVITY 2020 | f_RS = 0.965 ± 0.042 when fitted jointly with f_SP |

- **Schwarzschild precession.** f_SP scales the 1PN terms; GR gives 12.1′ per orbit, prograde.

  | Paper | f_SP | Data |
  |---|---|---|
  | GRAVITY 2020 | 1.10 ± 0.19 | S2 alone |
  | GRAVITY 2022 | 0.997 ± 0.144 | S2, S29, S38, S55 |
  | GRAVITY 2024 (A&A 692, A242) | **1.135 ± 0.110** (~10σ) | same 4 stars, S2 data to 2022.7 |

  - The newest value is GRAVITY 2024. The 2026 GRAVITY+ papers fix the precession at its GR value rather than refitting f_SP: these are the S301 discovery and the Yukawa limits with S301.
  - Caveats:
    - An extended mass (stars, remnants, dark matter) inside S2's orbit causes *retrograde* Newtonian precession, which is degenerate with f_SP. GRAVITY 2024 limits it to ≲1200 M☉ (1σ).
    - The NACO reference-frame zero points are partly degenerate with f_SP, and the results depend on the priors used for them. S2 alone with the 2024 prior gives 0.918 ± 0.128.

### 3.6 Cosmological and galactic scales
- **E_G (Reyes et al. 2010).** E_G = 0.39 ± 0.06 (1σ) at z ≈ 0.32, on 10–50 h⁻¹ Mpc.
  - GR's value is not a pure number: 0.408 ± 0.029 for their Ωm,0 = 0.2565 and f ≈ Ωm(z)^0.55.
  - As a ratio: 0.956 ± 0.162 (computed here, with the GR error included).
- **γ at ~2 kpc (Collett et al. 2018).** 0.97 ± 0.09 (68%) from ESO 325-G004.
  - The error is dominated by a systematic: velocity dispersions differ by 2.2% between stellar template libraries, which becomes 4.4% in dynamical mass and 8.8% in γ.
  - The lensing mass scales with the assumed Planck H₀: its 1.3% uncertainty becomes 2.6% in γ.
  - The analysis assumes γ is constant over the lens scale.
- **Planck 2018 μ₀ − 1 and Σ₀ − 1 (Table 7, 68%).**
  - Headline combination, Planck TT,TE,EE+lowE+lensing + BAO/RSD + DES-Y1 WL: μ₀ − 1 = −0.07 ^{+0.19}_{−0.32} and Σ₀ − 1 = 0.018 ^{+0.059}_{−0.048}.
  - Without CMB lensing, Planck alone gives Σ₀ − 1 = 0.27 ^{+0.15}_{−0.13}. Planck calls this "largely another reflection of the preference for A_L > 1", the lensing anomaly.
  - Parametrisation: μ(z) = 1 + E₁₁ Ω_DE(z), late-time only, ΛCDM background.
- **DESI (Ishak et al. 2025, JCAP 09, 053; DR1).** μ₀ = 0.05 ± 0.22 and Σ₀ = 0.008 ± 0.045 (means, 68%).
  - Data: DESI full-shape + BAO, Planck **PR4** (LoLLiPoP+HiLLiPoP) *without* CMB lensing, DES Y3 and DES Y5 SNe, ΛCDM background.
  - DESI's μ₀ is μ(today) − 1, the same quantity as Planck's μ₀ − 1.
  - With Planck PR3 the Σ₀ tension exceeds 3σ (Σ₀ = 0.388 ^{+0.11}_{−0.086}). It "goes away" with the PR4 likelihoods.
  - No DESI DR2 full-shape modified-gravity paper was found as of October 2026. Prefer DESI over Planck for Σ₀.
- **Growth index (Nguyen, Huterer & Wen 2023).** γ = 0.633 ^{+0.025}_{−0.024} against GR's 0.55, a **3.7σ** preference for suppressed growth. Planck PR3 + fσ8 gives 4.2σ.
  - This is the scorecard's one real "tension" point.
  - It is tied to the S8 tension and uses Planck PR3, the same likelihood family behind the Σ₀ anomaly that disappears with PR4. Whether γ survives PR4 was not checked here.
  - Table I does not state its interval level; the figures use 68%/95%.

## 4. Placing points by strength of gravity (approximate, computed here)

Order-of-magnitude dimensionless potential Φ/c² ≈ GM/(rc²) for each test, for a "potential vs precision" chart. These are my estimates, not published values. For a published compilation of potential and curvature by test, see Baker, Psaltis & Skordis 2015 (ApJ 802, 63); it was not checked here.

| Test | Φ/c² (approx.) | basis |
|---|---|---|
| LAGEOS / LARES-2 | 3.6 × 10⁻¹⁰ | GM⊕/c² = 4.435 mm, a = 12,270 km (ED Table 1) |
| Gravity Probe B | 6.3 × 10⁻¹⁰ | r ≈ 6371 + 642 km |
| Cosmology (E_G, μ₀, Σ₀, γ growth) | ~10⁻⁵ | large-scale-structure potentials (order of magnitude) |
| Collett lens | ~10⁻⁶ | galaxy potential (order of magnitude) |
| S2 at pericentre | 3.6 × 10⁻⁴ | 120 AU ≈ 1400 R_S → 1/2800; (v/c)² = 6.5 × 10⁻⁴ at 7650 km/s |
| EHT shadows | ~0.3 | photon orbit r ≈ 3GM/c² |
| GW inspiral to merger / ringdown | ~0.1–0.3 | ISCO (1/6) to light ring (1/3) |
| GW170817 speed | n/a | propagation test over ~40 Mpc; plot separately |

## 5. Suggested one-point-per-test choices (1σ after conversion)

| Test | Row | Point to plot | Precision |
|---|---|---|---|
| Geodetic precession | gpb2011_geodetic | ratio 0.99935 ± 0.00277 | 0.28% |
| Frame dragging, gyroscope | gpb2011_framedrag | ratio 0.949 ± 0.184 | 19% |
| Frame dragging, orbits | lares2_2026_nature | μ = 1.001 ± 0.0022 (disputed; see §3.2) | 0.2% (claimed) |
| GW propagation speed | gw170817_speed | bound, ≲10⁻¹⁵ | – |
| Graviton mass | gwtc4_graviton | < 1.92 × 10⁻²³ eV/c² (90%), arrow | – |
| GW inspiral (0PN) | gw250114_phi0 | 0.00 ± 0.018 | ~2% |
| IMR consistency | gw250114_imr_dMf | 0.02 ^{+0.043}_{−0.036} | ~4% |
| Ringdown frequency | gw250114_df220 (single event) or gwtc5_ringdown_df220_hier (catalog) | 0.02 ± 0.012 or 0.01 ± 0.018 | 1–2% |
| Ringdown damping (catalog) | gwtc5_ringdown_dtau220_hier | 0.09 ^{+0.079}_{−0.061} | mild tail; see §3.3 |
| Area theorem | gw250114_arealaw | 4.4σ pass (no deviation value) | – |
| Shadow, Sgr A* | eht_sgra_2022_delta_vlti | −0.08 ± 0.09 (Kerr band −0.075 to 0) | ~9% |
| Shadow, M87* | eht_m87_2019_delta_stellar | −0.01 ± 0.17 | 17% |
| S2 redshift | gravity2019_f_redshift | 1.04 ± 0.05 | 5% |
| S2 precession | gravity2024_fsp | 1.135 ± 0.110 | 11% |
| γ at 2 kpc | collett2018_gamma | 0.97 ± 0.09 | 9% |
| E_G | reyes2010_EG | ratio 0.956 ± 0.162 | 16% |
| Σ₀ (lensing) | desi2024_Sigma0 | 0.008 ± 0.045 | 4.5% |
| μ₀ (growth) | desi2024_mu0 | 0.05 ± 0.22 | 22% |
| Growth index | nguyen2023_growth_index | ratio 1.151 ± 0.045 (3.7σ off GR) | 4% |

## 6. Limits of this check

- **Pound–Rebka and Pound–Snider** (optional) were not redone. The APS full texts are paywalled, and another project has already verified them.
- **LARES-2 2026.** The main text, including its Eq. (2), is paywalled. The values come from Extended Data Table 2 and the peer-review file. The "≈3 yr" span is inferred from the 1050-day fits.
- **A&A and MDPI** block automated access. The EHT M87* 2018 Paper II (A&A 693, A265) was read only as an abstract, via Crossref and Semantic Scholar. Paper I of that series was not checked. The Lucchesi 2020 central value is from the arXiv version; its published abstract confirms the errors.
- **Preprints.** The GWTC-4.0 and GWTC-5.0 tests-of-GR papers are preprints, and their numbers may change on publication.
- **Not covered here:** pulsar tests (double pulsar, J1141−6545 frame dragging, pulsar spin precession) and lunar laser ranging belong to the other checks.
