# Cassini 2002: the Shapiro frequency shift and γ from the raw Doppler data

*Analysis notes for project N2 (the Cassini part). Every number below comes from `data/cache/cassini_fit.json`, `data/cache/cassini_inventory.csv` and the code in `cassini_*.py`; `cassini_doppler.run_cassini()` rebuilds them.*

## 1. Result

Main analysis: the plasma-free combination of the simultaneous X/X and Ka/Ka two-way links at DSS 25 (Goldstone), modelled with the navigation team's own SCE1 trajectory plus a fitted correction; 967 points of 600 s from 21 daily tracks, 2002-06-09 – 2002-07-05, impact parameter 1.8–42 solar radii:

- λ − 1 = (−0.27 ± 0.84 stat ± 0.65 syst ± 0.78 trajectory) × 10⁻⁴, where λ = (1 + γ)/2 scales GR's Shapiro rate;
- **γ − 1 = (−0.5 ± 1.7 stat ± 1.3 syst) × 10⁻⁴**, and **(−0.5 ± 2.6) × 10⁻⁴** once the choice of starting trajectory (± 1.6 × 10⁻⁴) is counted.

The statistical error is the largest of the formal error (0.45 × 10⁻⁴ in λ, from a noise model fitted to the data), a delete-one-track jackknife (0.77 × 10⁻⁴) and a bootstrap over whole tracks (0.84 × 10⁻⁴, 200 resamples). The systematic error is the rms spread about the main value of twelve alternative analysis choices (§7). The trajectory term is half the difference from the same fit started from Horizons' reconstructed trajectory instead of the navigation team's (−1.8 ± 3.1 × 10⁻⁴ in λ − 1): the two agree within their errors, and the navigation trajectory fits the data better and gives a more stable result (§6), but the difference is the largest single change we found, so we count it.

GR (γ = 1) is well inside the errors (0.2 σ with all terms; 0.3 σ with statistical and systematic only), and so is Bertotti, Iess & Tortora's γ − 1 = (2.1 ± 2.3) × 10⁻⁵ (Nature 425, 374, 2003) from the same experiment. Our uncertainty is 9–11 times larger than theirs, for reasons set out in §8. The Shapiro signal itself is detected at overwhelming significance: λ is measured to 0.01 %.

The frequency shift Cassini measured, with the fitted trajectory correction removed, is in `cassini_measured_series.parquet` (§10), beside GR's prediction: it follows the curve, from +0.4 to +1.2 × 10⁻¹⁰ before conjunction (June 9–17; there is no Ka/Ka closer in) and from about −4 × 10⁻¹⁰ on June 21 back to −0.3 × 10⁻¹⁰ on July 5 (y > 0 when the delay grows, §10; figure 09 plots Δν/ν = −y), with residuals of 1.6 × 10⁻¹⁴ rms per 10-minute point beyond 20 R☉.

## 2. What the archive holds (inventory)

Closed-loop Doppler of the first Solar Conjunction Experiment, PDS volumes CORS_0021–0028 (`cassini_inventory.csv` lists every data set per day, station and link). Two-way data, good-quality 1-s points:

| station | link | source | days | points | hours |
|---|---|---|---|---|---|
| DSS 15 | X/X | ODF | 2 | 31555 | 8.8 |
| DSS 25 | Ka/Ka | ODF | 22 | 615743 | 179.0 |
| DSS 25 | Ka/Ka | TDF only | 7 | 64457 | 18.4 |
| DSS 25 | X/X | ODF | 29 | 898017 | 251.9 |
| DSS 34 | X/X | ODF | 3 | 14896 | 4.1 |
| DSS 45 | X/X | ODF | 24 | 161371 | 134.5 |
| DSS 54 | X/X | ODF | 5 | 93940 | 26.7 |
| DSS 65 | X/X | ODF | 22 | 468891 | 153.2 |

Also present but not used: 777,693 three-way points (uplink and downlink at different complexes, two station clocks) and 98,547 one-way points (spacecraft oscillator).

What matters for the plasma calibration:

- **DSS 25 carried X/X and Ka/Ka together on 23 of the 29 days** with two-way data (the main fit uses those from June 9, where the navigation trajectory starts). There is no Ka/Ka on June 18, 19 and 20, the days closest to conjunction before it (the Ka uplink was off), nor on July 1 and only 12 minutes on June 28. The closest Ka/Ka track is June 21 (b = 1.7–2.3 R☉, 4.9 h).
- **The X-up/Ka-down link (X/Ka) is not in the closed-loop data at all.** DSS 25 had one Ka-band closed-loop channel and it tracked the Ka/Ka carrier (for about an hour on June 8 it held the X/Ka carrier instead; those records are labelled Ka/Ka and are screened out). X/Ka exists only as open-loop Radio Science Receiver recordings (`rsr/…x25k25d.1b1`, ~120 MB per pass), which this analysis does not process. Bertotti et al.'s three-link plasma cancellation is therefore not reproducible from the closed-loop archive; the two-link combination is.
- The ODFs lack two complete Ka/Ka passes (June 14 and 17) that the Archival Tracking Data Files (TDF) hold; we read them from the TDFs. Their Ka uplink frequency is in no ramp table, so it is estimated from the simultaneous X/X data and those passes carry a free bias in every fit.

Usable DSS 25 tracks after cleaning (elevation ≥ 15° on both legs), hours of 60-s normal points per link, and each track's share of the Fisher information on λ in the main fit, after the trajectory correction and the free biases are projected out (with the fitted noise model; the prior on the acceleration holds the other 4.9%):

| track | start (UTC) | b (R☉) | X/X h | Ka/Ka h | combined h | info on λ | mean residual (10⁻¹⁴) |
|---|---|---|---|---|---|---|---|
| 25_2002-06-06T17 | 2002-06-06 17:08 | 44.6–45.7 | 9.3 | 2.0 | 2.0 | 0.0% |  |
| 25_2002-06-07T17 | 2002-06-07 17:05 | 41.5–42.7 | 9.3 | 4.0 | 4.0 | 0.0% |  |
| 25_2002-06-08T17 | 2002-06-08 17:02 | 38.5–39.7 | 8.4 | 3.6 | 3.6 | 0.0% |  |
| 25_2002-06-09T16 | 2002-06-09 16:59 | 35.4–36.6 | 9.3 | 9.3 | 9.3 | 0.8% | −0.1 |
| 25_2002-06-10T16 | 2002-06-10 16:56 | 32.4–33.6 | 9.3 | 9.3 | 9.3 | 1.6% | +0.0 |
| 25_2002-06-12T18 | 2002-06-12 19:00 | 26.3–27.2 | 7.1 | 7.1 | 7.1 | 0.2% | +0.1 |
| 25_2002-06-13T16 | 2002-06-13 16:46 | 23.2–24.4 | 9.3 | 9.3 | 9.2 | 1.8% | −0.1 |
| 25_2002-06-14T16 | 2002-06-14 16:42 | 20.1–21.3 | 9.3 | 9.3 | 9.3 | 0.0% | +0.0 |
| 25_2002-06-15T16 | 2002-06-15 16:39 | 17.0–18.2 | 9.3 | 9.3 | 9.3 | 2.4% | −0.3 |
| 25_2002-06-16T16 | 2002-06-16 16:36 | 14.0–15.1 | 9.3 | 9.2 | 9.3 | 14.2% | −0.5 |
| 25_2002-06-17T16 | 2002-06-17 16:33 | 10.9–12.1 | 6.7 | 6.9 | 4.4 | 0.0% | +0.0 |
| 25_2002-06-18T16 | 2002-06-18 16:29 | 7.9–9.0 | 9.3 | 0.0 | 0.0 | 0.0% |  |
| 25_2002-06-19T16 | 2002-06-19 16:26 | 4.9–6.0 | 9.3 | 0.0 | 0.0 | 0.0% |  |
| 25_2002-06-20T16 | 2002-06-20 16:23 | 2.2–3.1 | 9.2 | 0.0 | 0.0 | 0.0% |  |
| 25_2002-06-21T16 | 2002-06-21 16:20 | 1.7–2.3 | 9.2 | 4.9 | 4.8 | 0.0% | +128.1 |
| 25_2002-06-22T16 | 2002-06-22 16:16 | 3.9–5.0 | 9.3 | 9.0 | 9.0 | 1.0% | +28.6 |
| 25_2002-06-23T16 | 2002-06-23 16:13 | 6.8–8.0 | 9.3 | 9.3 | 9.1 | 32.2% | +0.6 |
| 25_2002-06-24T16 | 2002-06-24 16:10 | 9.9–11.0 | 9.0 | 8.2 | 8.2 | 0.5% | −1.2 |
| 25_2002-06-25T16 | 2002-06-25 16:07 | 12.9–14.1 | 9.3 | 9.2 | 9.2 | 23.0% | −0.8 |
| 25_2002-06-26T16 | 2002-06-26 16:03 | 16.0–17.2 | 9.3 | 5.8 | 5.8 | 2.1% | +0.7 |
| 25_2002-06-27T15 | 2002-06-27 16:00 | 19.1–20.3 | 9.2 | 9.2 | 9.2 | 1.1% | −0.1 |
| 25_2002-06-29T01 | 2002-06-29 01:02 | 23.3–23.3 | 0.2 | 0.2 | 0.2 | 0.0% | +0.0 |
| 25_2002-06-29T15 | 2002-06-29 15:53 | 25.2–26.4 | 9.3 | 9.3 | 9.2 | 4.8% | +0.1 |
| 25_2002-06-30T15 | 2002-06-30 15:50 | 28.3–29.5 | 9.3 | 9.2 | 9.2 | 3.2% | +0.1 |
| 25_2002-07-01T20 | 2002-07-01 20:26 | 32.0–32.0 | 0.7 | 0.0 | 0.0 | 0.0% |  |
| 25_2002-07-02T15 | 2002-07-02 15:43 | 34.4–35.4 | 7.6 | 7.5 | 7.5 | 0.7% | −0.1 |
| 25_2002-07-03T21 | 2002-07-03 21:48 | 38.3–38.7 | 3.2 | 3.0 | 3.0 | 3.1% | −0.5 |
| 25_2002-07-04T15 | 2002-07-04 15:37 | 40.6–41.7 | 9.0 | 9.0 | 9.0 | 2.4% | +0.6 |

June 14 and 17 (and the 12-minute track 25_2002-06-29T01) carry none: with a free bias, only the signal's small change within the pass counts (§9). June 21 carries none because of its plasma noise (§5).

## 3. The model

Everything is computed for each 1-s observable, in the BCRS with TDB as time, and nothing is fitted at this stage:

- **Observable.** The ODF Doppler observable is F = K·f_ref − f_R, with f_ref the record's reference frequency, K = 880/221 for an X-band downlink and 3344/221 for Ka band. Verified: K·f_ref − F is continuous when f_ref changes in mid-pass, and the X/X and Ka/Ka received frequencies divided by their turnaround ratios (880/749 and 14/15) times their uplink frequencies agree to 10⁻⁹ (plasma). The TDF counts give the identical observable (to 5 × 10⁻⁷ Hz on 74,000 points).
- **Time.** ODF time tags are UTC calendar seconds since 1950, at the middle of the 1-s count. TT = UTC + 64.184 s; TDB − TT from ERFA's series plus the topocentric term v⊕·x_station/c².
- **Ephemerides.** JPL Horizons' SSB vectors of the Earth and the Sun (2-minute steps; positions by cubic Hermite interpolation, velocities by a cubic spline of the printed velocities). For Cassini, the navigation team's SCE1 short-arc trajectory "SM727D" (J. Bordi, August 2002: fitted to DSS 25's plasma-calibrated X-band Doppler of June 9 – July 5, estimating the state, the RTG acceleration and solar pressure), read from the archive's SPK (SPICE transfer format, type 1, decoded here) — the base of the γ fit; and Horizons' reconstruction, which covers all of SCE1 (used for the inventory-wide products and as a variant).
- **Stations.** DSN 810-005 coordinates (ITRF93), solid-Earth tides (IERS 2010, degree 2 and 3, up to 0.3 m), rotation to the GCRS with IERS Earth orientation (astropy), then to the BCRS with Moyer's relativistic transformation.
- **Light time.** Solved per leg (up: station at t1 → Cassini at t2; down: t2 → station at t3) including the Sun's Shapiro delay (1PN, plus the 2PN terms of Klioner & Zschocke 2010, with the Sun where it was at the ray's closest approach as Bertotti, Ashby & Iess 2008 prescribe), the troposphere (zenith wet delay from DSS 25's water-vapour radiometers, dry from surface pressure, Niell mapping) and the Earth's ionosphere (DSN GPS-based calibrations, scaled as 1/f²).
- **Doppler.** Differentiating the two light-time equations exactly: dt2/dt3 = (1 − n_d·v_S3/c − dD_d/dt3)/(1 − n_d·v_C/c) and dt1/dt2 = (1 − n_u·v_C/c − dD_u/dt2)/(1 − n_u·v_S1/c), with the delay rates from full solutions 30 s either side. The received frequency is f_R = M₂ f_T(t1) (dt1/dt3) (dτ/dt)₁/(dτ/dt)₃, the last factor being the station clock's rate at transmission and reception.
- **Signal.** The model is formed with the Shapiro rate on and off; GR's predicted frequency shift is y_GR = 1 − F_on/F_off, and the residual y = (f_R,model − f_R,observed)/f_R,model (Shapiro off) is what relativity, trajectory errors, plasma and noise must explain. Away from the Sun this residual is a smooth offset of 1.1–1.4 × 10⁻¹⁰ (a 16–20 mm/s line-of-sight velocity difference from Horizons' reconstruction) plus GR's signal.

## 4. Checks, and the bugs they found

On a quiet pass (June 6, b = 45 R☉) the residual of 60-s points after a linear trend has rms 4.8 × 10⁻¹³ (X/X) and 0.98 × 10⁻¹³ (Ka/Ka) with every model term in; removing terms shows each is needed:

| change | X/X rms | Ka/Ka rms |
|---|---|---|
| full model | 4.8 × 10⁻¹³ | 0.98 × 10⁻¹³ |
| station clock-rate term with the wrong sign | 1.0 × 10⁻¹¹ | 3.0 × 10⁻¹³ |
| no station clock-rate term | 5.0 × 10⁻¹² | 1.7 × 10⁻¹³ |
| DSN troposphere instead of the radiometers | 4.8 × 10⁻¹³ | 1.35 × 10⁻¹³ |
| no troposphere or ionosphere | 1.6 × 10⁻¹² (bias +1.5 × 10⁻¹² below 25°) | 1.35 × 10⁻¹³ |

Time tags shifted by ±0.5 s triple the 1-s residual (6.7 → 20 × 10⁻¹²); by 32 s (TAI-like tags) it becomes 10⁻⁹–10⁻⁸: the UTC reading is right. X/X and Ka/Ka give the same residual offset (1.0927 vs 1.0935 × 10⁻¹⁰), which checks both turnaround ratios. The trajectory partials were checked against numerically perturbed ephemerides (agreement to 0.1%; the line-of-sight position partial, which is tiny, to 6%). GR's y_GR agrees with the project's independent geocentric code (`shapiro.py`) to 3 × 10⁻⁴ when both use the geocentre; at DSS 25 it is ~1% smaller because every pass was in daytime, when the station's rotation opposes the Earth's orbital motion and slows the sweep of the ray across the Sun.

Bugs found and fixed in the earlier model (λ = 0.953 ± 0.007 from X/X alone):

1. **Horizons epochs rounded to 0.86 ms.** The CSV's Julian dates carry 8 decimals; at 30 km/s the Earth is misplaced by up to 13 m at a node, and a Hermite spline on those times made 2-minute Doppler wiggles of ~3 × 10⁻¹⁰ — the dominant "noise" of the old residuals (37 × 10⁻¹²). The 2-minute grid is now rebuilt exactly.
2. **Cassini's velocity from rounded positions.** Horizons prints Cassini's coordinates (10⁹ km) to 1 cm; differentiating the position spline turned that into velocity errors of ~5 × 10⁻⁵ m/s inside each 2-minute interval (3 × 10⁻¹³ in y, 4 × 10⁻¹⁴ in 600-s averages). Velocities now come from the printed velocities (good to 10⁻¹¹ m/s). This halved the noise of the Ka/Ka residuals on quiet days.
3. **Shapiro delay about the wrong centre.** It was computed with SSB-centred vectors as if heliocentric; the Sun was 1.17 R☉ from the SSB (mostly along the line of sight, 0.03 R☉ across it). The Sun's own vectors are now fetched from Horizons.
4. **Modelling 60-s averages instead of 1-s points.** The curvature of the diurnal Doppler over a minute is ~2 × 10⁻¹² (more with gaps); the model is now evaluated at every second and averaged afterwards.
5. **No solid-Earth tides.** The 0.3-m tidal motion of the antenna adds ~2 × 10⁻¹³ to the two-way Doppler and shifts the daily means by up to 10⁻¹³ as the Moon's phase drifts through the month. Without the tides λ − 1 is −3.6 ± 4.0 × 10⁻⁴ (variant `no_solid_tides`): they move λ by +3.4 × 10⁻⁴ and cut the jackknife error 5-fold.
6. **ODF reader.** The data-summary group (key 105) was not recognised, so its records were read as ramp records of the last station.
7. **No plasma calibration** (the main cause of λ = 0.953 with the bugs above): X/X alone gives λ − 1 = −207 × 10⁻⁴ here too (§7).
8. **The least-squares solver.** The first version of the fit inverted the normal matrix; with columns spanning 20 orders of magnitude and nearly degenerate trajectory directions, that silently discarded information and made λ depend on the starting trajectory by up to 5 × 10⁻⁴. The fits now solve the whitened system by SVD with normalised columns (checked against an independent least-squares solve to 10⁻¹⁰).

Data problems found and screened out: a Ka/Ka stretch of June 7 (19:40–21:03) offset by −5.7 × 10⁻⁶ with 10⁻⁸ scatter (a non-coherent downlink), the June 8 X/Ka stretch labelled Ka/Ka (−1.57 × 10⁻⁴), and an X/X stretch of June 17 (16:33–18:36) offset by +3.2 × 10⁻⁴; plus uplink acquisition sweeps, points below 15° elevation and 1-s outliers beyond 5 robust σ.

## 5. Plasma: the two-link combination

The corona and the Earth's ionosphere add A/f² to the fractional frequency on each leg (A ∝ the rate of change of the electron column): y_link = y₀ + A_u/f_u² + A_d/f_d². With P = 1/f_u² + 1/f_d² (P_X = 0.0335 GHz⁻² for X/X at 7.2/8.4 GHz, P_Ka = 0.00182 GHz⁻² for Ka/Ka at 34.3/32.0 GHz: X/X carries 18.4 times more plasma),

  y₀ = (P_X·y_KaKa − P_Ka·y_XX)/(P_X − P_Ka) = 1.0576 y_KaKa − 0.0576 y_XX,

exact when A_u = A_d. That holds for the slowly changing corona that the line of sight sweeps through, the part that could mimic the Shapiro rate; for turbulence that differs between the legs the residue is 2.2 × 10⁻⁴ (A_d − A_u), about 4 times less than Ka/Ka's own and 77 times less than X/X's. The combination is formed second by second, so its noise comes from its own scatter. With the X/Ka link the separation would be exact (Bertotti, Comoretto & Iess 1993): y₀ = 1.0457 y_KaKa + 0.0287 y_XKa − 0.0744 y_XX.

What the data show (median white noise of a 60-s point at DSS 25):

| b (R☉) | Ka/Ka | X/X | X/X+Ka/Ka |
|---|---|---|---|
| (0, 3] | 1.1e-11 | 1.5e-10 | 1.8e-11 |
| (3, 5] | 7.6e-13 | 2.2e-11 | 1.3e-12 |
| (5, 8] | 1.9e-13 | 6.1e-12 | 1.5e-13 |
| (8, 12] | 1.4e-13 | 3.5e-12 | 8.6e-14 |
| (12, 20] | 5.9e-14 | 8.5e-13 | 4.2e-14 |
| (20, 30] | 4.0e-14 | 4.7e-13 | 3.4e-14 |
| (30, 50] | 4.3e-14 | 3.5e-13 | 3.9e-14 |

Away from the Sun the combination is quieter than Ka/Ka alone, and both are far quieter than X/X. Inside 5 R☉ it is noisier than Ka/Ka: the closed-loop X/X data there carry noise that does not scale as 1/f² (the receiver tracking a strongly scintillating carrier), which the −0.0576 weight passes on. So the closest days, where the Shapiro signal is 3–4 × 10⁻¹⁰, add little: June 21 and 22 together carry about 1% of the information on λ (§2).

Plasma by link (λ − 1, × 10⁻⁴): X/X alone −207; Ka/Ka alone −10.1; the combination −0.27. The mean corona lowers the apparent Shapiro rate (its phase advance grows as the ray approaches the Sun, opposite to the delay), as the earlier λ = 0.953 showed. A fitted coronal template (density ∝ r^−2.3, each leg's column ∝ b^−1.3) brings Ka/Ka alone back to +2.5 ± 7.3 with an error 4–9 times the combination's, but cannot rescue X/X alone (−223 ± 395: the real corona is not a smooth power law); added to the combination it changes nothing (+0.5 ± 1.3).

## 6. The fit

y = λ·y_GR + trajectory correction + free biases of the passes with an estimated Ka uplink + noise, by generalized least squares on 600-s averages:

- **Trajectory correction.** Both reconstructions were themselves fitted to Doppler with γ = 1, and they differ by 12–15 mm/s along the line of sight. The fit corrects Cassini's position and velocity at June 21.5 and adds a constant acceleration (RTG heat, solar pressure): nine parameters through the exact partials of the two-way Doppler.
- **Which trajectory to start from.** Horizons' and the navigation team's trajectories differ by a function that the nine-parameter correction absorbs to 3 × 10⁻¹⁵ rms, so in a linear world they would give the same λ. They nearly do: with identical fixed weights (b > 5 R☉) λ − 1 is −1.0 × 10⁻⁴ from Horizons and −0.2 × 10⁻⁴ from the navigation SPK, and with the data-driven noise model −1.8 ± 3.1 × 10⁻⁴ (Horizons, same span) against the main −0.27 — the small unabsorbed remainder lines up a little with the Shapiro shape once the weights emphasise a few quiet tracks. With the same nine parameters and weights the navigation trajectory fits better (χ² lower by 40 for 884 points), needs a smaller correction (the Horizons-based fit wants accelerations of ~4 × 10⁻⁸ m/s² across the line of sight, which a prior must hold back), and its λ moves less with the noise model. So the navigation trajectory is the base; the Horizons results are kept as variants, and half the difference is counted as a separate systematic (§1). A Gaussian prior of 10⁻⁸ m/s² per acceleration component keeps it physical (left free, the acceleration runs to 3 × 10⁻⁸ m/s², almost all of it across the line of sight (2 × 10⁻⁹ m/s² along it), nine times the RTG recoil; λ moves by −0.1 × 10⁻⁴). The fitted acceleration is 1.4 × 10⁻⁸ m/s² in all, 5.2 × 10⁻¹⁰ m/s² along the line of sight. The Shapiro signal is sharp near conjunction and antisymmetric about it; the correction is smooth; their largest correlation with λ is 0.64.
- **Noise.** Each pass has an excess noise above the white error of its points, with AR(1) correlation between successive bins; each track has a random offset common to its points, of variance s₀² + (κ·excess)². All are re-estimated from the residuals until stable: away from the Sun the excess is 7.1 × 10⁻¹⁵ per 600-s point (median), near it up to 2 × 10⁻¹¹; the offsets come out at s₀ = 0, κ = 0.25 (offsets of ~2.5 × 10⁻¹⁵ per track away from the Sun). χ²/dof = 0.93.
- **Uncertainty.** Away from the Sun the day-to-day errors are small: beyond 6 R☉ the post-fit daily mean residuals (§2) scatter by 0.5 × 10⁻¹⁴ rms, no more than the averaged within-pass noise and those small offsets predict, and none exceeds 1.2 × 10⁻¹⁴. Inside 5 R☉ residual plasma leaves 2.9 × 10⁻¹³ (June 22) and 1.3 × 10⁻¹² (June 21). But with only 21 tracks the fitted offset variance is itself uncertain, so the quoted statistical error is the largest of the formal error and two resampling estimates over whole tracks (here 1.8 times the formal one).

## 7. Variants

All variants re-fit λ with everything else as in the main analysis. "jk" is the delete-one-track jackknife. The error on γ − 1 in the last column is twice the larger of the formal and jackknife errors; for the main fit the bootstrap (§1) raises it from 1.5 to 1.7.

| variant | what changes | points | tracks | χ²/dof | (λ−1)·10⁴ | formal | jackknife | (γ−1)·10⁴ ± 2·max(formal, jk) |
|---|---|---|---|---|---|---|---|---|
| main | X/X + Ka/Ka plasma-free combination at DSS 25, navigation SPK trajectory + state/acceleration correction, 600-s bins, AR(1) noise per pass + random track offsets | 967 | 21 | 0.93 | −0.27 | 0.45 | 0.77 | −0.5 ± 1.5 |
| kaka_only | Ka/Ka link alone (no plasma calibration) | 983 | 21 | 1.00 | −10.06 | 2.48 | 5.10 | −20.1 ± 10.2 |
| kaka_only_plasma_template | Ka/Ka alone + fitted coronal template (b^-1.3) | 983 | 21 | 1.01 | +2.46 | 4.25 | 7.33 | +4.9 ± 14.7 |
| xx_only | X/X link alone at DSS 25 (no plasma calibration) | 1208 | 25 | 1.02 | −207.50 | 43.00 | 51.95 | −415.0 ± 103.9 |
| xx_only_plasma_template | X/X alone + fitted coronal template | 1208 | 25 | 1.02 | −223.35 | 67.05 | 395.43 | −446.7 ± 790.9 |
| combo_plasma_template | combination + template for its plasma residue | 967 | 21 | 0.93 | +0.52 | 0.55 | 1.30 | +1.0 ± 2.6 |
| combo_plus_xx_no_ka_template | combination + DSS 25 X/X where Ka was off (June 18-20), template | 1206 | 25 | 0.96 | −0.42 | 0.35 | 0.89 | −0.8 ± 1.8 |
| combo_plus_xx_other_stations | combination + X/X of DSS 34/45/54/65/15 with b > 15 Rsun (station biases) | 1640 | 52 | 0.84 | −0.28 | 0.54 | 0.73 | −0.6 ± 1.5 |
| bin_300s | 300-s bins | 1938 | 21 | 0.97 | −0.37 | 0.35 | 0.61 | −0.7 ± 1.2 |
| bin_1800s | 1800-s bins | 324 | 20 | 0.95 | −0.48 | 0.60 | 0.89 | −1.0 ± 1.8 |
| noise_white | independent errors (per-pass excess, no AR(1), no track offsets) | 967 | 21 | 0.94 | −0.38 | 0.31 | 0.71 | −0.8 ± 1.4 |
| noise_ar1_no_offsets | AR(1) per pass, no random track offsets | 967 | 21 | 0.90 | −0.41 | 0.35 | 0.69 | −0.8 ± 1.4 |
| noise_common_offset | AR(1) + random track offsets of one common variance | 967 | 21 | 0.94 | −0.34 | 0.33 | 0.77 | −0.7 ± 1.5 |
| trajectory_no_acceleration | trajectory: state only | 967 | 21 | 0.92 | −0.56 | 0.51 | 1.21 | −1.1 ± 2.4 |
| trajectory_acceleration_free | trajectory: state + free acceleration | 967 | 21 | 0.93 | −0.40 | 0.46 | 0.77 | −0.8 ± 1.5 |
| trajectory_poly3 | line-of-sight velocity as a cubic in time | 967 | 21 | 0.52 | −5.90 | 2.54 | 4.84 | −11.8 ± 9.7 |
| trajectory_poly5 | line-of-sight velocity as a quintic in time | 967 | 21 | 0.51 | −8.78 | 3.42 | 6.85 | −17.6 ± 13.7 |
| track_offsets_2e-14 | per-track offset parameters with a 2e-14 Gaussian prior (instead of random offsets) | 967 | 21 | 0.84 | −2.17 | 1.09 | 1.32 | −4.3 ± 2.6 |
| track_offsets_free | free per-track offsets (within-pass shape only) | 967 | 21 | 0.92 | −19.28 | 9.05 | 8.22 | −38.6 ± 18.1 |
| b_gt_5 | impact parameter > 5 Rsun | 884 | 19 | 0.92 | −0.18 | 0.36 | 0.76 | −0.4 ± 1.5 |
| b_gt_8 | impact parameter > 8 Rsun | 828 | 18 | 0.92 | +0.11 | 0.42 | 1.17 | +0.2 ± 2.3 |
| b_gt_12 | impact parameter > 12 Rsun | 752 | 16 | 0.91 | +0.14 | 0.48 | 1.24 | +0.3 ± 2.5 |
| b_lt_25 | impact parameter < 25 Rsun | 584 | 13 | 0.97 | −0.74 | 0.58 | 1.35 | −1.5 ± 2.7 |
| without_june_21_22 | the two closest tracks removed | 884 | 19 | 0.92 | −0.18 | 0.36 | 0.76 | −0.4 ± 1.5 |
| ingress_only | before conjunction only (June 9-17) | 402 | 8 | 0.91 | −3.14 | 11.22 | 5.80 | −6.3 ± 22.4 |
| egress_only | after conjunction only (June 21 - July 4) | 565 | 13 | 0.91 | −4.87 | 2.68 | 3.16 | −9.7 ± 6.3 |
| elevation_10 | elevation cut 10 deg (AWVR data only above 15) | 1032 | 21 | 0.95 | −0.51 | 0.39 | 0.47 | −1.0 ± 0.9 |
| elevation_20 | elevation cut 20 deg | 890 | 20 | 0.94 | −0.36 | 0.44 | 0.66 | −0.7 ± 1.3 |
| elevation_25 | elevation cut 25 deg | 812 | 20 | 0.94 | −0.24 | 0.40 | 0.47 | −0.5 ± 0.9 |
| tropo_dsn | DSN troposphere calibration instead of the AWVRs | 968 | 21 | 1.01 | −0.94 | 0.44 | 0.46 | −1.9 ± 0.9 |
| awvr2 | water vapour radiometer unit 2 instead of unit 1 | 967 | 21 | 0.92 | −0.27 | 0.46 | 0.68 | −0.5 ± 1.4 |
| no_ionosphere | no Earth-ionosphere calibration | 967 | 21 | 0.93 | −0.26 | 0.46 | 0.78 | −0.5 ± 1.6 |
| no_solid_tides | station not displaced by the solid-Earth tides (a check: the tides are real) | 967 | 21 | 0.73 | −3.64 | 1.20 | 4.02 | −7.3 ± 8.0 |
| no_2pn | the 2PN light-time terms left out of the model (a check) | 967 | 21 | 0.93 | −0.46 | 0.47 | 0.82 | −0.9 ± 1.6 |
| horizons_trajectory | Horizons' reconstructed trajectory instead of the navigation SPK (June 6 - July 5) | 1025 | 24 | 0.92 | −2.63 | 0.61 | 1.60 | −5.3 ± 3.2 |
| horizons_trajectory_span | Horizons' trajectory, same span as the main fit (June 9 - July 5) | 967 | 21 | 0.94 | −1.83 | 0.55 | 3.12 | −3.7 ± 6.2 |
| horizons_simple_weights | Horizons' trajectory, main span, b > 5 Rsun, fixed weights (white + 1.5e-14), no prior | 884 | 19 | 2.10 | −1.00 | 0.27 | 0.68 | −2.0 ± 1.4 |
| nav_simple_weights | b > 5 Rsun, fixed weights (white + 1.5e-14), no prior | 884 | 19 | 2.06 | −0.19 | 0.27 | 0.53 | −0.4 ± 1.1 |

The systematic error quoted in §1 is the rms deviation from the main value of: `combo_plasma_template`, `bin_300s`, `bin_1800s`, `noise_white`, `noise_ar1_no_offsets`, `trajectory_no_acceleration`, `trajectory_acceleration_free`, `track_offsets_2e-14`, `b_gt_8`, `elevation_20`, `tropo_dsn`, `awvr2`.

Remarks:
- The trajectory correction can be swapped (state only, free acceleration), and so can the troposphere calibration (radiometer unit 1 or 2, or the DSN's model), the ionosphere calibration (it cancels in the combination), the elevation cut, the bin size and the noise model, without moving λ beyond the quoted errors.
- Polynomials in the line-of-sight velocity instead of the state correction cannot represent the daily signature of a position error (χ²/dof drops to 0.5 because the noise model swallows it); they pull λ − 1 to −5.9 ± 4.8 and −8.8 ± 6.9 × 10⁻⁴, within their own errors.
- Starting from Horizons' trajectory instead of the navigation team's gives −2.6 ± 1.6 × 10⁻⁴ over June 6 – July 4 and −1.8 ± 3.1 × 10⁻⁴ over the main span (§6). Explicit per-track offsets with a 2 × 10⁻¹⁴ prior instead of the fitted random offsets give −2.2 ± 1.3 × 10⁻⁴; these two are the largest changes among the reasonable alternatives.
- Removing the tracks closest to the Sun (b < 5, 8 or 12 R☉) moves λ by less than its error. The information is in the tracks at 7–25 R☉, three quarters of it after conjunction: June 23 (b ≈ 7 R☉) and June 25 (b ≈ 13 R☉) carry 55% of it between them, and June 16 (b ≈ 14.5 R☉) another 14% (§2). The cuts at 8 and 12 R☉ remove June 23, and the jackknife error grows by 50–60%.
- Free per-track offsets throw away the day-to-day change of the signal and leave only its shape within a pass; λ is then determined only to ~1 × 10⁻³ (−19 ± 9 × 10⁻⁴).
- Fits to one side of conjunction only are weak (the trajectory correction absorbs most of a one-sided signal) but agree.

## 8. Why the precision is lower than Bertotti et al.'s

1. **Most of the gap is our own analysis.** Anderson, Lau & Giampieri (2004) fitted the same kind of data, closed-loop Doppler from JPL's navigation system combined as 1.0575901 y_KaKa − 0.0575901 y_XX like ours, with JPL's Orbit Determination Program, and got γ − 1 = (−1.3 ± 5.2) × 10⁻⁵ (random error only). That is 2.3 times Bertotti et al.'s ± 2.3 × 10⁻⁵, so the third link, which only the open-loop receivers recorded, explains at most a factor 2.3 of the factor 11 between our total error and theirs; the remaining factor 5 is ours. With independent errors weighted pass by pass, as Anderson et al. weighted theirs, our formal error would be ± 0.61 × 10⁻⁴ in γ (variant `noise_white`), close to theirs on a similar amount of data (161 h of 600-s points against their 131 h of 300-s points); but that fit's own jackknife error is ± 1.4 × 10⁻⁴, so independent errors understate the scatter. The correlated noise and day-to-day offsets of our noise model give ± 0.91 × 10⁻⁴, and the rest of the error budget (§1) multiplies that by 2.9. Counted in variance, the total of ± 2.6 × 10⁻⁴ is 12% formal error, 29% the bootstrap's excess over it (± 1.4 × 10⁻⁴), 24% the spread of the twelve variants (± 1.3 × 10⁻⁴) and 35% the trajectory term (± 1.6 × 10⁻⁴).
2. **No three-link plasma calibration.** Bertotti et al. combined X/X, X/Ka and Ka/Ka, which cancels the plasma exactly. The archive's closed-loop data have no X/Ka (§2), so we use the two-link combination, which leaves the up/down difference, and near the Sun inherits the closed-loop X/X data's excess noise. The passes within 5 R☉ (June 21–22), where the Shapiro signal is largest (3–4 × 10⁻¹⁰), therefore carry only about 1% of the information here (§2). Bertotti et al. had little there either: they used no pass on June 21 and only a handful of points on June 22, their noisiest pass by far (Supplementary Fig. S1), and Ashby & Bertotti (2010) give their minimum impact parameter as ≈ 6 R☉. The third link matters at 7–15 R☉. The three passes there that carry most of our information (June 16, 23 and 25: 69%) are 1.3–5 times noisier than our passes beyond 20 R☉ (600-s residual rms 2.9, 7.7 and 2.1 × 10⁻¹⁴, against 1.6 × 10⁻¹⁴), whereas Bertotti et al.'s three-link data were about as quiet on those days as far from the Sun (1–2.5 × 10⁻¹⁴ at 1000 s, Supplementary Fig. S1). The X/Ka recordings exist in the open-loop RSR files; processing them (the RSR headers give the receiver's frequency model; the carrier's residual frequency must be tracked in the 1 kHz complex samples) is the way to the full calibration.
3. **Excess noise of up to 1–3 × 10⁻¹⁴ per 600-s point** limits the far tracks (signal 0.3–1.2 × 10⁻¹⁰): 6 of the 16 tracks beyond 12 R☉ carry that much above their white noise, correlated from point to point. Its likely sources are errors of the troposphere calibration (on the far passes the radiometer and DSN calibrations differ by 1–3.5 × 10⁻¹⁴ from point to point, though by only ~0.5 × 10⁻¹⁴ in their daily means), and ocean-tide loading (~1 cm at Goldstone) and sub-daily Earth-orientation variations, which are not modelled here. The day-to-day part is smaller: the daily mean residuals scatter by 0.5 × 10⁻¹⁴ (§6). Bertotti et al. and Anderson et al. both fitted with JPL's Orbit Determination Program ("the most complete and sophisticated model available for fitting DSN data", in Anderson et al.'s 2003 preprint gr-qc/0308010, since withdrawn); we have the radiometers' zenith delays, our own station model (§3) and a nine-parameter correction to a reconstructed trajectory.

Not a difference: the DSS 25 Ka-band transmitter failed on June 18–20, the days of the largest signal before conjunction, and the published analyses lost those days too ("no data could be acquired for three days just before conjunction owing to a failure of the ground transmitter", Nature Fig. 2 caption; Anderson et al. 2004 count four scheduled passes). It limits every analysis of this experiment equally.

## 9. Caveats

- The trajectory we start from was fitted by JPL navigation with γ = 1, to these very data (plasma-calibrated with X and Ka). This matters only if its model had freedom that a Shapiro-shaped signal could hide in; it used smooth forces (state, RTG, solar pressure), which our nine-parameter correction spans to a few × 10⁻¹⁵ in y. Remainders of that size can still move λ by ~10⁻⁴ under data-driven weights, as the Horizons comparison shows (§6), which is why the trajectory term is in the error budget.
- The Ka uplink frequencies of June 14 and 17 (and 12 minutes on June 28) are estimated from X/X; those passes carry free biases, so they inform λ only through their shape within the pass.
- The closed-loop X/X data near the Sun are not consistent with a pure 1/f² plasma; the noise model deweights them rather than modelling them.
- Not modelled (each ~10⁻¹⁴, mostly diurnal): ocean-tide loading, sub-daily Earth orientation, antenna thermal deformation, the cm-level difference between ITRF93 and the frame of the Earth-orientation series, the planets' Shapiro delays (~10⁻¹⁵).
- The 2PN terms are in the model at their GR values (dropping them moves λ by −0.19 × 10⁻⁴; variant `no_2pn`).

## 10. Files and how to rebuild

- Code: `cassini_odf.py` (ODF reader), `cassini_tdf.py` (TDF reader), `cassini_spk.py` (SPK transfer-format reader, type 1), `cassini_media.py` (troposphere, ionosphere, Niell mapping, radiometer files), `cassini_doppler.py` (data, model, residuals, normal points, inventory; `run_cassini()` rebuilds everything), `cassini_plasma.py` (link combinations, plasma template), `cassini_gamma.py` (fits, jackknife, bootstrap, variants, outputs).
- Downloads (`fetch_data.py --cassini`): the ODFs and calibration files; added for this analysis the Sun's SSB vectors (Horizons), the X-uplink TDFs (~1.2 GB: the Ka/Ka of June 14 and 17), the radiometer path delays (pd1/pd2), the navigation SPK and the archive documents.
- Outputs in `data/cache/`: `cassini_fit.json` (main result with its noise parameters and nuisance values, every variant, the summary), `cassini_measured_series.parquet` (per 600-s point: date_utc, t_tdb, b_rsun, b_up_rsun, station, link, track, y_measured, sigma, y_gr, residual, dataset — the main combination plus X/X-only and Ka/Ka-only fits with and without the plasma template), `cassini_inventory.csv`, the 1-s observables and residuals and the 60-s normal points per model configuration.
- For the figure: y_measured (data minus fitted trajectory correction and biases) against date or b_rsun, with y_gr as GR's curve; average to an hour for display. Signs: y > 0 when the delay grows (before conjunction).
