# Binary and triple pulsars — literature check notes (atlas N7, gravity scorecard)

Companion to `pulsars.csv` (361 rows). Checked 2026-10-19.

## 1. How the numbers were checked

- Every value in the CSV was read from the primary paper's own TeX source (arXiv), in the table or sentence quoted. The PDF text was used to confirm table, equation and figure numbers. The sources are stored under `data/literature/papers/<short-name>/`.
- A script checked each `quote` field: after collapsing whitespace, every quote occurs verbatim in its source file. Two papers have no arXiv version: Taylor & Weisberg 1982 and 1989. Their quotes were transcribed by eye from the ADS page scans, because the ADS OCR text garbles the ± signs.
- Parentheses follow the usual convention: digits in parentheses are the error in the last digits, e.g. `4.226585(4)` = 4.226585 ± 0.000004.
- **Voisin et al. 2020/2025 use a different notation.** `1.62940(06)_{-30}^{+30}` means the value is 1.6294006 and the errors (−30/+30) apply to the digits in parentheses, i.e. ±0.0000030. The CSV gives the decoded number in `value`/`err_plus`/`err_minus` and the original TeX cell in `quote`. I checked the decoding against their own statement that the A18 vs V20 masses differ by 2.9σ, 2.6σ and 2.4σ. With the decoding, (1.4401 − 1.4359)/0.0015 = 2.8σ, (0.19780 − 0.19730)/0.00019 = 2.6σ and (0.41058 − 0.40962)/0.00040 = 2.4σ.
- Error conventions differ between papers, and the `err_type` column records each one:
  - WH16: formal TEMPO 1σ errors.
  - WNT10: "estimated" errors of 2–6 × formal.
  - Kramer+2006: twice the TEMPO errors, quoted as 1σ.
  - Kramer+2021 and Hu+2022: 1σ from Monte Carlo.
  - Ransom+2014 and Archibald+2018: 1σ from MCMC.
  - Voisin+2020/2025: 68% intervals, with the headline Δ given at 95%.
  - Taylor & Weisberg 1982: "conservative" adopted errors.
  - `not_stated` means the paper gives a ± without saying what confidence level it is.
- Row status:
  - `verified`: the draft value was confirmed.
  - `corrected`: the draft value was wrong (details in `notes`).
  - `added`: a value not in the draft, checked in the primary source.
  - `secondary_only`: found only in a review, or not printed in the primary paper.

## 2. Corrections to the draft

| Draft item | What the primary source says |
|---|---|
| WH16 Galactic correction Ṗb(gal) ≈ −0.0128e-12 | **−(0.025 ± 0.004)e-12**. The value −(0.0128 ± 0.0050)e-12 belongs to Weisberg & Taylor 2005, which gave the ratio 1.0013 ± 0.0021. |
| K21 Pb = 0.10225156248(5) d | That is the **Kramer 2006** value. K21 gives **0.1022515592973(10) d**. |
| K21 e = 0.0877775(9) | That is also the 2006 value. K21 gives **e_T = 0.087777023(61)** (the DD "time eccentricity"). |
| Ransom 2014 outer Pb = 327.257541(66) d | arXiv TeX table and PDF: **327.257541(7) d**. The Nature typeset version could not be opened (paywall), so this rests on the arXiv v1 source; I found no "(66)" anywhere. |
| Voisin 2020 \|Δ\| < 2.1e-6 | Published value **2.05e-6** (rounds to 2.1). Marked verified. |
| WH16 γ "= 0.004307(4) s" | The value is correct, but the paper's table header says "(ms)". The number is in **seconds** (4.307 ms): the public tempo file has GAMMA = 0.004306748 s, and WNT10 gave 4.2992(8) ms. Bhat+2008 makes the same slip for J1141−6545: "γ (ms) 0.000773(11)" is 0.773 ms. |

Everything else in the draft was confirmed:
- **B1913+16:** WH16 Pb, e, x, ω̇ and Ṗb(obs); masses 1.438/1.390; Ṗb(GR) −2.40263(5)e-12; ratio 0.9983 ± 0.0016; WNT10 0.997 ± 0.002.
- **Double pulsar (K21):** x_A, ω̇, γ_E, Ṗb, R = 1.0714(11), m_A, m_B, ratio 0.999963(63), and "1.3e-4 at 95%". The (63) is the 1σ error; the 95% statement is the paper's own.
- **Double pulsar (K06):** s ratio 0.99987 ± 0.00050.
- **Double pulsar (Breton+2008):** Ω_B = 4.77 +0.66/−0.65 °/yr (68%), GR value 5.0734 ± 0.0007.
- **Triple:** R14 masses and inner Pb; A18 \|Δ\| < 2.6e-6 (95%); V20 Δ = (+0.5 ± 1.8)e-6 (95%).
- **Pulsar–white dwarf binaries:** Antoniadis 1.05 ± 0.18; Zhu \|Δ\| < 0.002 (95%); Bhat 1.04 ± 0.06.

Other things a careful reader should know:
- **WH16 TOA count.** The abstract and Sect. 2 say 9257 TOAs, but the Conclusions say 9309. The public tempo file has NTOA 9257.
- **Hu+2022 mass of B.** Table 3 gives m_B = 1.248886(5), but the same table has M − m_A = 2.587050 − 1.338186 = 1.248864. One of the numbers must be a slip; use the Table 2 masses.
- **TW89 typo.** Sect. Vc writes "m2 = 1.442 ± 0.003" when substituting; the abstract has m1 = 1.442 and m2 = 1.386.
- **WT05 sign.** WT05 prints Ṗb(corrected) = (2.4056 ± 0.0051)e-12 without the minus sign.
- **V20 signal size.** Voisin 2020, Sect. 3.2 says \|Δ\| ~ 1e-6 gives "a signal amplitude of ~0.1 s via the Rømer delay". The same wording is in the final-version source on Zenodo. It must be a typo for ~0.1 µs (see §5.4).

## 3. B1913+16 (Hulse–Taylor)

- **Use WH16** (ApJ 829, 55): Table 1 for spin and astrometry, Table 2 for orbital parameters, plus the text values.
  - PK parameters: ω̇ = 4.226585(4) °/yr, γ = 4.307(4) ms, Ṗb(obs) = −2.423(1)e-12, s = 0.68 +0.10/−0.06, r = 9.6 +2.7/−3.5 µs. This is the first Shapiro-delay measurement in this system; s and r are 68% intervals from a χ² grid.
  - δθ(obs) = 4.0(25)e-6 cannot be used as a test, because it is mixed with an unknown aberration term.
- **Masses assume GR.** m1 = 1.438(1) and m2 = 1.390(1) M☉ come from ω̇ + γ with the GR equations. The Ṗb test is a three-parameter (ω̇–γ–Ṗb) test.
- **Ṗb correction.** WNT10, WH16 and Deller+2018 all use the Damour & Taylor (1991) "galactic" term. Deller+2018 Eqs. 2–3 write it out: planar differential Galactic rotation plus the Shklovskii term (Pb/c)·μ²d. The vertical acceleration is neglected because b = 2.1°.
  - WH16: R0 = 8.34 ± 0.16 kpc and Θ0 = 240 ± 8 km/s (Reid+2014), DM distance 9.8 ± 3.1 kpc → Ṗb^gal = −(0.025 ± 0.004)e-12, Ṗb^intr = −(2.398 ± 0.004)e-12.
  - Ratio to GR: **0.9983 ± 0.0016**.
- **The test is distance-limited.** Deller+2018 (ApJ 862, 139) measured a VLBI parallax of 0.24 +0.06/−0.08 mas, giving d = 4.1 +2.0/−0.7 kpc. Combined with the WH16 timing, this gives Ṗb^gal = −(0.008 +0.009/−0.003)e-12 and a ratio of **1.005 +0.001/−0.003**. That is "a formal discrepancy from general relativity of 1.6σ" and, in their words, "similar precision and lower accuracy" than WH16. The honest scorecard statement is agreement at the few-tenths-of-a-percent level, limited by the distance, not by the timing.
- **History of the ratio** (all rows in the CSV):

  | Paper | Ratio Ṗb(obs)/Ṗb(GR) | Galactic correction applied? |
  |---|---|---|
  | TW1982 | −2.30(22)e-12 vs −2.403(5)e-12, i.e. ≈ 0.96 ± 0.09 (our arithmetic; no ratio printed) | none |
  | TW1989 | 1.010 ± 0.011 (abstract: 1.01 ± 0.01) | none |
  | WT2005 | 1.0013 ± 0.0021 | yes |
  | WNT2010 | 0.997 ± 0.002 | yes |
  | WH2016 | 0.9983 ± 0.0016 | yes |
  | Deller+2018 | 1.005 +0.001/−0.003 | yes (VLBI distance) |

- **Data are public.** WH16 put the TEMPO input file and all 9257 TOAs on Zenodo (10.5281/zenodo.54764, CC0). The same files are in the arXiv source, already downloaded to `data/literature/papers/weisberg-huang-2016/src/`.
  - The file gives more digits than the paper, e.g. GAMMA 0.004306748 ± 0.000003466 s and OMDOT 4.2265854 ± 0.0000044.
  - The time scale is TDB. The README warns that the file will not work natively in tempo2.
  - Its SINI and M2 lines are single-fit values, not the paper's grid-refined s and r.

## 4. Double pulsar J0737−3039A/B

- **Use K21** (PRX 11, 041050; arXiv v2 TeX): Table IV for parameters, Table V for test ratios, Eqs. 35–48 for masses and the Ṗb budget.
- **Timing model and time scale.**
  - The timing model is DDS, so the fitted Shapiro shape parameter is z_s = −ln(1 − s) = 9.65(15), and s = 0.999936 +9/−10 is derived from it.
  - The eccentricity is the DD "time eccentricity" e_T.
  - All values are in TDB. Times must be divided by (1 − 1.550519768e-8) for TCB/tempo2, and T0 needs the IAU 2006 conversion (paper footnote).
- **Data set.** 916,648 30-s TOAs over MJD 52760–58640 with 26.0 µs rms; also 199,913 4-min TOAs with 13.0 µs rms.
- **Masses assume GR plus a neutron-star equation-of-state prior.**
  - m_A = 1.338185 +12/−14 and m_B = 1.248868 +13/−11 M☉ come from k (≡ ω̇) and s.
  - k includes a 2PN term (~35σ significant) and a Lense–Thirring term from A's spin, ω̇^LT,A = −4.83 +0.29/−0.35 × 10⁻⁴ °/yr, which is more than 30× the measurement error. The LT term needs the moment of inertia I_A, which was taken from GW170817 radius limits.
  - Masses in M☉ use the IAU nominal GM☉ = 1.3271244e26 cm³ s⁻². They are "Doppler-shifted" masses (true masses × 1/D), which does not matter for the tests.
  - **Practical warning for the mass–mass diagram:** a plain 1PN ω̇ curve will miss the intersection of the other curves by many σ at K21 precision. Either include the 2PN and LT terms (K21 Sect. V), or plot the ω̇ band with a realistic I_A range, as K21 Fig. 13 does.
- **Ṗb budget (K21 Eqs. 39–48).**

  | Term | Value (s/s) |
  |---|---|
  | Observed | −1.247920(78)e-12 |
  | Shklovskii | +1.72(15)e-16 |
  | Galactic (McMillan 2017 potential) | −3.40 +0.25/−0.24 × 10⁻¹⁶ |
  | Shklovskii + Galactic | −1.68 +0.11/−0.10 × 10⁻¹⁶ |
  | Intrinsic | −1.247752(79)e-12 |
  | Mass loss from A's spin-down | +2.9(2)e-17 |
  | **Compare with GR:** Ṗb^GW | **−1.247782(79)e-12** |
  | GR quadrupole (2.5PN) | −1.247810 +6/−7 × 10⁻¹² |
  | 3.5PN correction | −1.75e-17 |
  | GR total | −1.247827 +6/−7 × 10⁻¹² |
  | **Ratio** | **0.999963(63)**, "1.3e-4 with 95% confidence" |

  The Shklovskii and Galactic terms have opposite signs, so distance errors nearly cancel. With a wide uniform I_A prior the ratio is 0.999958(64).
- **Other tests (Table V).**
  - s: 1.00009(18), from the ω̇–γ_E–s test; "4e-4 at 95%".
  - r: 1.0016(34); "7e-3 at 95%".
  - γ_E: 1.00012(25).
  - ω̇: 1.000015(26).
  - δθ: 1.3(13). This is δθ = 13(13)e-6 against a GR value of 12.60889(4)e-6.
  - Ω_B: 0.94(13), taken from Breton+2008.
  - λ_LT: 0.7(9).
  - q_NLO: 1.15(13) total, 1.26(24) from deflection only, 1.32(24) from retardation only.
  - Moment of inertia: I_A < 3.0e45 g cm² (90%).
- **Other double-pulsar inputs.**
  - The mass ratio R = 1.0714(11) is the 2006 measurement; it has not been improved because pulsar B vanished in 2008.
  - Breton+2008 measured B's spin precession from eclipse modelling: 4.77 +0.66/−0.65 °/yr (68%) against a GR value of 5.0734(7). They also give c²σ_B/𝒢 = 3.38 +0.49/−0.46.
- **Newer result (added).** Hu et al. 2022 (A&A 667, A149) timed A with MeerKAT for 3 years.
  - Shapiro shape twice as precise: z_s = 9.669(77), s = 0.9999369 +46/−51. Range: r = 6.163(16) µs.
  - Masses: m_A = 1.338186(10), m_B = 1.248866(7) M☉, consistent with K21.
  - q_NLO = 0.999(79).
  - γ_E and Ṗb were *fixed* to the K21 values, so this paper does not give a new radiative test. Its x differs from K21 by ~3σ for that reason.
  - It is a confirmation and refinement, not a replacement of K21. No combined 16-year + MeerKAT paper was found up to October 2026.
- **Data availability.** Neither the K21 nor the Hu+2022 source contains a data-availability statement or a public TOA link (I grepped for "available", "zenodo", "github" and "request"). I found no public double-pulsar TOA set.

## 5. Triple system PSR J0337+1715

### 5.1 Parameter sets: use one paper's set, never mix

| Paper | Reference epoch | Equations of motion | Time scale / ephemeris | Notes |
|---|---|---|---|---|
| Ransom+2014 Table 1 | MJD 55920.0 | Newtonian three-body + special-relativistic transverse Doppler; Einstein and Shapiro delays added afterwards | DE405; time scale not stated (TEMPO2 used for barycentring) | 14 orbital parameters + spin; derived masses and geometry |
| Archibald+2018 Tables E.1–E.3 | MJD 55920.0 | 1PN; PPN β, γ with Solar-System priors; Δ | **TCB**, DE435 | pulsar proper time *not* rescaled |
| Voisin+2020 Table 2 | **MJD 56492** (position epoch 57205) | 1PN modified-EIH; Damour–Deruelle post-Keplerian elements mapped to initial conditions | not printed (a DE430 table line is commented out of the TeX); time scale not stated | includes Kopeikin, Shklovskii and aberration delays |
| Voisin+2025 Table E.1 | MJD 56492 | as V20 (runs with Δ ≡ 0) | DE430 | **ω and Ω conventions differ from V20**, e.g. ω_p = −7.06° vs V20 ω_I = 97.0° |

In every set the elements are *osculating* values at the reference epoch. The inner set describes the pulsar about the inner-binary centre of mass; the outer set describes the inner-binary centre of mass about the system centre of mass. In R14 the longitude of the inner ascending node is set to zero.

They are not comparable across papers. For example, P_b,I is 1.629401788(5) d in R14, 1.6293932(6) d in A18 and 1.6294006(30) d in V20. The R14–A18 gap of 8.6e-6 d is about 14 times A18's quoted error and about 1,700 times R14's. Voisin+2020 says direct comparison "is not straightforward because of i) slightly different definition of the parameters and ii) the fact that most of the parameters are not constants of motion". Only the masses and Δ can be compared directly.

**For the simulation in (b), use R14 Table 1 with Newtonian equations of motion.** It is the most complete and self-describing set: fixed, fitted and derived values, including q_I, δΩ, the inclinations and all three masses. Add Δ as G_pi = G_po = G(1+Δ), G_io = G. That is Voisin's definition, and at Newtonian order it is the same as giving the pulsar m_G = (1+Δ)m_I with the white dwarfs obeying the strong equivalence principle, as in A18.

### 5.2 Definition of Δ

- **A18:** "Following Damour and Schäfer we define Δ = m_G/m_I − 1, the fractional difference between the pulsar's inertial (m_I) and gravitational (m_G) masses." The table note describes the same quantity as "the fractional difference in acceleration between the pulsar and the inner white dwarf". The white dwarfs are treated as weak-field bodies. In quasi-Brans–Dicke theories, \|Δ\| = \|α_o(α_p − α_i)\|.
- **V20 and V25:** G_ab = G_N(1 + Δ_ab), where Δ_ab is a symmetric "relative GWEP parameter" between bodies a and b. Only the pulsar is strongly self-gravitating, so G_pi ≃ G_po ≃ G_N(1+Δ) and G_io ≃ G_N, with Δ ≡ Δ_pb(s_b = 0).
  - In the weak-field limit Δ_ab ≃ Δ_a + Δ_b, with Δ_a = (m_G/m)_a − 1.
  - In Bergmann–Wagoner scalar–tensor theories Δ ≃ −2ζ s_p. V20 gives ζ s_p = (−0.2 ± 0.9)e-6 (95%).
  - V20 also argues that m_G/m_I is not a meaningful per-body quantity in the strong field because of the s_a s_b cross term. Freire & Wex 2024 make the same point.
- **So Δ ≈ (m_G/m_I)_pulsar − (m_G/m_I)_WD** in both definitions. The white dwarfs' own Δ is roughly η·ε_grav ≲ 1e-5 × 1e-4, which is negligible.
- **Zhu+2019 (J1713+0747):** the same Δ, but the pulsar and white dwarf fall in the Galaxy's field (about 2e-10 m s⁻²) instead of the outer white dwarf's (about 1.7e-3 m s⁻²). That is why the limit is about 1000× weaker: \|Δ\| < 0.002 (95%).

### 5.3 Results and their disagreement

| Paper | Data | Δ | 95% limit |
|---|---|---|---|
| A18 | 27,110 TOAs, AO + GBT + WSRT, MJD 55956.7–57866.9 | (−1.09 ± 0.74)e-6 (1σ, including systematics) | \|Δ\| < 2.6e-6 |
| V20 | 9,303 Nançay TOAs, MJD 56492–58761 | (+0.5 ± 1.8)e-6 at 95%; 4.8 +9.5/−9.4 × 10⁻⁷ at 68% | \|Δ\| < 2.05e-6 |
| V25 (newer, added) | 12,474 Nançay TOAs to MJD 59480 | red-noise model: 0.97 +1.58/−1.48 × 10⁻⁶; planet model: 0.51 +1.14/−1.28 × 10⁻⁶ | 2.29e-6 (red noise); 1.46e-6 (planet) |

- **A18's error budget.** The formal (statistical) result was Δ = −1.1e-6 ± 0.2e-6, about 5σ from zero. The authors attributed the excess to unmodelled systematics such as solar wind and scattering, which can put correlated signal at the signature frequency. They modelled those systematics empirically, as random sinusoids at k·f_I + l·f_O whose 1σ scatter is 22 ns at the Δ-signature component. That inflated the error to 0.74e-6, so the result is "systematics-limited".
- **V20 is "statistics-limited."** Its noise is absorbed by an EFAC of 1.3153, which accounts for about 8% of the error, and it identifies red noise in the residuals.
- **Other reasons V20 gives for the difference (V20 Sect. 6.1):**
  1. It uses an independent, uniformly sampled single-telescope data set.
  2. Its timing model includes Kopeikin and Shklovskii delays and aberration. A18 used a local linear approximation for astrometry, and its timing astrometry disagrees with Gaia.
  3. V20 treats Δ like every other MCMC parameter and argues that A18's systematics model "heavily depends on the modelling choices for which no physical justification is currently available".
  4. The masses disagree: A18 m_p = 1.4359(3) vs V20 1.4401(15) (2.9σ); m_i 2.6σ; m_o 2.4σ. V20 suggests A18's unreported systematic error on the masses partly relieves this.
  5. V20 measures the outer node longitude Ω_O for the first time.
- Both results are consistent with Δ = 0 at 95%.
- **V25** finds the strong low-frequency residual is either exceptionally strong achromatic red noise or a ~0.5 lunar-mass planet. Chromatic, dispersion-measure noise is ruled out. The SEP limit depends on which model is right, by about 1.6×, so a scorecard should quote the range: \|Δ\| < (1.5–2.3)e-6 at 95%, with "model-dependent" noted.
- The Freire & Wex 2024 review still quotes V20 (Δ_p0 = (0.5 ± 1.8)e-6) and ω_BD ≳ 150,000. V20's own abstract says > 140,000.

### 5.4 Expected signal size and shape, for simulation (b)

- **A18, measured:** the Δ signature, i.e. what is left after refitting all other parameters with Δ ≡ 0, "is a sinusoidal variation of pulse arrival times whose amplitude is 30 ns when Δ is 10⁻⁶", at frequency **2f_inner − f_outer**.
  - It scales linearly with Δ.
  - For the best-fit Δ the data show 33 ns of signature against 22 ns of systematics (1σ).
  - For Δ = 2.6e-6 the largest component is 78 ns.
  - V20 also places the main signature at 2f_i − f_o.
- **Our own cross-check (not from the papers).** The differential acceleration Δ·g_o, with g_o = G·m_o/a_O,rel² ≈ 1.75e-3 m s⁻², rotates with the outer orbit. This forces an inner eccentricity e_F ≈ (3/2)·Δ·g_o/(n_I·a_R·n_O) ≈ 5.5e-8 per Δ = 1e-6. Here a_R = 15.9 lt-s is the inner relative semimajor axis, n_I = 2π/1.6294 d and n_O = 2π/327.26 d.
  - The pulsar's displacement is e_F·a_p ≈ 0.1 µs.
  - In the Rømer delay this appears as an eccentricity term (x·e_F/2) at 2n_I − n_O, about **34 ns per 1e-6**. That agrees with A18's 30 ns and confirms both the "0.1 µs" reading of V20's "0.1 s" and the signal frequency.
  - The term is a steady forced oscillation: it does not grow with time. The error on Δ shrinks roughly as σ_TOA/√N, and with coverage of many outer orbits, unless red noise dominates.
  - Rule of thumb (ours): with N TOAs of white-noise error σ, the error on the amplitude is ≈ σ·√(2/N). For A18 (σ ≈ 1.2 µs, N = 27,110) that is ≈ 10 ns, i.e. σ_Δ ≈ 0.35e-6. A18's formal statistical error was 0.2e-6, which is helped by the 0.6 µs Arecibo TOAs.
- **Without refitting** the other parameters, a change in Δ also rescales the effective masses, and the orbital phase drift grows to milliseconds over years. That part is completely degenerate with the masses and periods. **The simulation must refit the other parameters, or project out their partial derivatives, before comparing with the 30-ns scale.** V20 makes the same point: the correlations reduce what can be detected to "~50 ns".
- **Real timing precision.**
  - R14: 26,280 TOAs, weighted rms 1.34 µs; best TOAs 0.8 µs in 10 s.
  - A18: weighted rms 1.2 µs overall; Arecibo 0.6 µs (5,355 TOAs), GBT 1.5 µs (18,487), WSRT 1.8 µs (3,268).
  - V20: 20-min TOAs in four 128-MHz subbands; mean error 2.15 µs, median 1.89 µs; EFAC 1.3153.
  - V25: computed by us from the public `.tim` file: 12,474 active TOAs, mean error 2.18 µs, median 1.88 µs; EFAC 1.1177 in the GR run.

### 5.5 Data availability (triple system)

- **Archibald+2018**
  - Data: the paper says "Data is available from the authors upon request."
  - Code: on GitHub at https://github.com/aarchiba/triplesystem (no license file; last push 2018-05-30). The repository contains early TOA files from June 2013 to September 2014 (`code-v1/*.tim`), not the final 27,110-TOA set; the README describes the raw data as on the group's private servers.
- **Voisin+2020: public.** Zenodo **10.5281/zenodo.3778978** (CC-BY-4.0) holds:
  - `Data_and_Results.tar.bz2` (30.8 MB): Nançay TOAs and MCMC results;
  - `nutimo.tar.bz2`: the timing code (built with tempo2 2018.02.1, Minuit2 and Boost 1.55);
  - the final article sources, which I downloaded. They are identical to the arXiv TeX.
- **Voisin+2025: public.** Zenodo **10.5281/zenodo.13899771** holds:
  - `Data.tar.bz2` (0.49 MB): two `.tim` files, the full set and the one used in the paper with MJD 58631–58780 excised. I downloaded this and counted 12,474 active TOAs over MJD 56492.305–59480.160, matching the paper.
  - `nutimo` code, `Analysis.tar.bz2` (247 MB), and the paper sources. No license is given on the record.
  - The files are in `data/literature/papers/voisin-2025-zenodo/`.
- **Ransom+2014:** no data statement in the arXiv source.

## 6. Pulsar–white dwarf tests (dipole radiation, Ġ, universality of free fall)

| System | Ṗb corrections | Result | Masses |
|---|---|---|---|
| J1738+0333 (Freire+2012) | Shklovskii 8.2 +0.6/−0.5 and Galactic 0.58 +0.16/−0.14 fs/s, from timing parallax 0.68(5) mas and proper motion. Ṗb(obs) = −17.0(3.1) → Ṗb^Int = −25.9(3.2) fs/s, vs GR −27.7 +1.5/−1.9 fs/s. | **Ṗb^xs = +2.0 +3.7/−3.6 fs/s**; κ_D = (−0.8 ± 1.6)e-4 (alone, s_p = 0.15); with J0437: Ġ/G = (−0.6 ± 1.6)e-12 /yr and κ_D = (−0.3 ± 2.0)e-4 (1σ) | Optical spectroscopy (M_c = 0.181 +0.008/−0.007, q = 8.1(2)), so at leading order the GR prediction does not depend on GR-derived masses |
| J0348+0432 (Antoniadis+2013) | Kinematic + Ġ correction 0.016 ± 0.003 × 10⁻¹³ s/s judged negligible and not applied | Ṗb = −0.273(45)e-12 vs GR −0.258 +0.008/−0.011 × 10⁻¹²; **ratio 1.05 ± 0.18** | Optical (M_WD = 0.172(3), M_PSR = 2.01(4)) |
| J1141−6545 (Bhat+2008) | Shklovskii (6.71 ± 1.73)e-15 (v_T from scintillation, d = 3.7 kpc) and Galactic (−5.05 ± 0.44)e-15 → intrinsic (−4.01 ± 0.25)e-13 | **Ratio 1.04 ± 0.06**; the GR Ṗb value itself is not printed | **Assume GR** (DDGR fit: M = 2.2892(3), m_c = 1.02(1), m_p = 1.27(1)) |
| J1141−6545 (Venkatraman Krishnan+2020) | — | ẋ = (1.7 ± 0.3)e-13 explained by the white dwarf's rotational quadrupole plus Lense–Thirring; white-dwarf spin period 116 +17/−70 s or 397 +78/−242 s, depending on the prior (68%) | **Assume GR** (DDGR). This is a consistency demonstration, not a ratio test, so no single number belongs on the scorecard. |
| J1713+0747 (Zhu+2019) | — | **\|Δ\| < 0.002** (95%), with −0.0007 < Δ < 0.0023; Ġ/G = (−0.1 ± 0.9)e-12 /yr (confidence level not stated; presumably 1σ); κ_D = (−0.7 ± 2.2)e-4; −3e-20 < α̂₃ < 4e-20 (95%) | — |

Freire+2012 does **not** print a ratio Ṗb^Int/Ṗb^GR. The 0.94 ± 0.13 often quoted, for example in WH16 Table 3, is secondary; it is −25.9/−27.7. The paper's own test quantity is Ṗb^xs.

## 7. Will 2014 Living Review cross-check

- **Table 6 (B1913+16)**
  - The orbital and post-Keplerian entries equal WNT10 Table 3, except that ω0 is printed with error (4) where WNT10 has (6).
  - The RA, Dec and pulsar period are **not** from WNT10, despite the caption; I could not identify their source.
  - Ṗ = 8.62713(8)e-18 matches Taylor & Weisberg 1989 Table 3 (J2000 column).
  - The text values (masses, Ṗb^GR, Ṗb^gal, ratio 0.997 ± 0.002) equal WNT10. All are superseded by WH16.
- **Table 7, double-pulsar column**
  - x_A, e, Pb and ω̇, and the text masses 1.337/1.250, come from **Lyne+2004**.
  - x_B is printed 1.513(3), but Lyne has 1.513(4).
  - s is printed 0.9995(4), dropping Lyne's asymmetric −0.0032.
  - γ′ = 0.382(5), Ṗb = −1.21(6) and r = 6.2(5) match neither Lyne+2004 nor Kramer+2006; their source is unknown.
  - The whole column is superseded by K21.
- **Table 7, J1738+0333 and J1141−6545 columns:** they match Freire+2012 and Bhat+2008. The Ṗb values are the corrected intrinsic ones, as the caption says.
- **Text on J0337, J0348 and J1141:** these values match R14, A13 and Bhat+2008.

## 8. Newer results, 2022–2026

- **Added:**
  - Hu+2022 (MeerKAT double pulsar, §4).
  - Voisin+2025 (J0337, §5.3).
  - Deller+2018 (B1913+16 VLBI distance; older than 2022, but it changes the B1913+16 test).
- **Secondary only:** the Freire & Wex 2024 Living Review.
- **Searches.** arXiv API abstract searches for "J0737-3039", "J0337+1715" and "B1913+16" over 2022-01-01 to 2026-10-19, plus web searches. They found no new timing analysis that supersedes K21/Hu+2022, V20/V25 or WH16/Deller+2018. One paper is out of scope: arXiv 2404.09864 (Universe 2024), "Measuring a gravitomagnetic effect with the triple pulsar PSR J0337+1715" (only the title was checked; not a timing result).
- **New triple system.** PSR J0435+3233 (Yang et al. 2026, ApJL 1010, L9; arXiv 2608.01227, abstract only): a γ-ray pulsar with an inner 8.0-d helium white dwarf and a G-subgiant tertiary on a 73.5-yr orbit. No strong-equivalence-principle test is reported, and with that outer period it is not competitive with J0337.

## 9. Practical pitfalls for the computations

1. **Use intrinsic Ṗb, not observed.** B1913+16: WH16 −2.398(4)e-12, or Deller's version. Double pulsar: K21 Ṗb^GW = −1.247782(79)e-12. Compare with the quadrupole formula, plus the 3.5PN term for the double pulsar.
2. **Solar-mass time constant.** Use T☉ = GM☉/c³ = 4.925490947 µs, as K06, WNT10 and F12 do.
3. **Time scales.** TDB for WH16 (README), K21 and Hu+2022. TCB for A18 and F12. Not stated for R14, V20 and V25. For dimensionful parameters the difference is a factor 1 − 1.55e-8, which is negligible for the plots but not for re-using T0 in tempo2.
4. **Double-pulsar ω̇** must include the 2PN and Lense–Thirring terms; see §4.
5. **J0337:** never combine elements from different papers, and refit before quoting a signal amplitude; see §5.
