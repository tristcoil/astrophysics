# Astrophysics Project Atlas

A master plan for astrophysics projects that are **scientifically accurate**, use **only real, freely available data**, and end in **visualizations worth keeping**. This document is the base we build from: each project below lists the science question, exactly where the data lives, how to access it, the analysis steps, the plots to make, and honest caveats.

---

## Ground rules

1. **Real data only.** Every figure is derived from mission archives or published catalogs. No synthetic stand-ins. (If we ever run a simulation — e.g. an N-body integration — it is labeled as a simulation, and its initial conditions still come from real data.)
2. **Honest analysis.** Uncertainties, quality cuts, and selection effects get documented in each project's own README. A pretty plot that hides its caveats is not the goal.
3. **Reproducible.** One folder per project. A fetch script downloads and caches raw data; figures are regenerable from scratch with one command.
4. **Credit the data.** Each archive has citation guidelines (Gaia, SDSS, GWOSC, etc.). Each project README notes its data provenance.

## Repository layout

```
astrophysics/
├── README.md              ← this file: the master plan
├── 01-exoplanet-transit/  ← one folder per project, numbered as we build them
│   ├── README.md          ← science question, data provenance, results, caveats
│   ├── fetch_data.py      ← downloads + caches raw data into data/
│   ├── analysis.py        ← (or a notebook) produces every figure from cached data
│   ├── data/              ← raw/cached data — never committed if we git init (gitignore it)
│   └── figures/           ← output plots
├── 02-.../
└── ...
```

## Toolbox

Python 3.10+ in a virtual environment. This machine has Python 3.14 but no global pip — a venv bootstraps its own:

```bash
cd ~/Desktop/astrophysics
python3 -m venv .venv
source .venv/bin/activate
pip install numpy scipy pandas matplotlib plotly astropy astroquery
```

Domain packages, installed per-project as needed:

| Package        | What for                                                    |
|----------------|-------------------------------------------------------------|
| `astropy`      | Units, coordinates, FITS files, cosmology models, time series — the backbone of everything |
| `astroquery`   | Programmatic access to Gaia, SDSS, JPL Horizons, MAST, VizieR, SkyView |
| `lightkurve`   | Kepler/K2/TESS light curves from the MAST archive            |
| `gwpy`, `gwosc`| LIGO/Virgo strain data and event catalogs                    |
| `healpy`       | HEALPix all-sky maps (Planck CMB)                            |
| `psrqpy`       | ATNF pulsar catalogue queries                                |
| `sunpy`        | Solar physics data (SDO imagery, GOES X-ray flux)            |
| `plotly`       | Interactive HTML figures (3D star maps, zoomable light curves) |
| `rebound`      | N-body orbital integrations (extension projects only)        |

> **Note on Python 3.14:** it is very new; if a compiled package (`healpy`, `gwpy`) has no wheel for it yet, the fallback is a venv on an older interpreter or a miniforge/conda environment. We'll deal with it per-project if it comes up.

---

## Project index

Difficulty: ★☆☆ = an afternoon · ★★☆ = a weekend · ★★★ = multi-session.

| #  | Project                                   | Field          | Difficulty | Download    | Payoff |
|----|-------------------------------------------|----------------|-----------|-------------|--------|
| A1 | Detect an exoplanet transit               | Exoplanets     | ★☆☆       | ~10–50 MB   | Find a planet in raw spacecraft data |
| A2 | Exoplanet population atlas                | Exoplanets     | ★☆☆       | ~5 MB       | The demographics of 5000+ worlds |
| B1 | Hertzsprung–Russell diagram from Gaia     | Stars          | ★☆☆       | ~50–200 MB  | All of stellar evolution in one plot |
| B2 | 3D map of the solar neighborhood          | Stars          | ★★☆       | ~100 MB     | Fly through real nearby stars |
| B3 | Hunt star clusters by their motion        | Stars          | ★★☆       | ~50 MB      | Clusters pop out of the noise |
| C1 | Gravitational waves: the GW150914 chirp   | Compact objects| ★★☆       | ~100 MB     | See + hear two black holes merge |
| C2 | Pulsar P–Ṗ diagram                        | Compact objects| ★☆☆       | <5 MB       | The whole neutron star zoo |
| C3 | Fast radio bursts (CHIME)                 | Compact objects| ★☆☆       | <5 MB       | The newest mystery in astronomy |
| D1 | Supernova Hubble diagram → dark energy    | Cosmology      | ★★☆       | ~10 MB      | Reproduce the 2011 Nobel Prize |
| D2 | The cosmic web from SDSS                  | Cosmology      | ★★☆       | ~50 MB      | Filaments and voids from real galaxies |
| D3 | Planck CMB map + acoustic peaks           | Cosmology      | ★★☆       | 0.3–1 GB*   | Baby picture of the universe |
| D4 | Dark matter from rotation curves (SPARC)  | Cosmology      | ★☆☆       | <5 MB       | Rubin's flat curves, reproduced |
| E1 | Kirkwood gaps in the asteroid belt        | Solar system   | ★★☆       | ~200 MB     | Resonances carve visible gaps |
| E2 | Solar system orbit visualizer (Horizons)  | Solar system   | ★★☆       | API, small  | Planets, ʻOumuamua, Voyager |
| E3 | Near-Earth object close approaches        | Solar system   | ★☆☆       | API, small  | What's passing by this month |
| F1 | 275 years of solar cycles + butterfly diagram | The Sun    | ★☆☆       | <10 MB      | The Sun's heartbeat since 1749 |
| F2 | The Sun today, in extreme UV (SDO)        | The Sun        | ★★☆       | ~100 MB     | Actual current images of the Sun |
| G1 | Color composites of real galaxies         | Imaging        | ★★☆       | varies      | Your own pretty pictures, honestly made |
| G2 | Meteorite landings world map              | Fun/easy       | ★☆☆       | ~10 MB      | Quick, satisfying, geographic |

\* D3 is small (~2 MB) if we only do the power spectrum; the full-resolution sky map is the heavy part.

**Suggested first arc:** A1 (transit) → B1 (HR diagram) → C1 (gravitational waves) → D1 (dark energy). That sequence ramps difficulty gently and touches four different subfields, four different archives, and four different kinds of visualization.

---

# A. Exoplanets

## A1. Detect an exoplanet transit ★☆☆

**Science question.** Can we find the periodic dimming of a star caused by an orbiting planet, directly in NASA spacecraft photometry — and measure the planet's size and orbital period from it?

**Background.** When a planet crosses (transits) its star's disk, the star dims by a fraction `δ ≈ (Rp/R★)²`. For a Jupiter crossing a Sun-like star that's ~1%; for an Earth, ~0.008%. Kepler (2009–2018) and TESS (2018–) measured stellar brightness precisely enough to see this. The dips repeat every orbital period, so folding the light curve on the right period stacks them into one sharp transit.

**Data.** Kepler/K2/TESS light curves, hosted at the MAST archive (STScI). Public domain, fetched automatically by `lightkurve`.

**Access.**
```python
import lightkurve as lk
search = lk.search_lightcurve("Kepler-10", mission="Kepler")
lc = search.download_all().stitch()   # all quarters, stitched
```

**Method.**
1. Download and stitch light curves for a known transiting system.
2. `lc.flatten()` — remove slow stellar variability, keep the transits.
3. Run a Box Least Squares (BLS) periodogram (`lc.to_periodogram(method="bls")`) to find the transit period — *we* find it, we don't look it up.
4. Phase-fold on the detected period; bin the folded curve.
5. Measure depth → planet radius via `Rp = R★ · √δ` (stellar radius from the NASA Exoplanet Archive); measure duration → sanity-check the orbit.
6. Compare our period/radius against published values — our numbers should land within a few percent.

**Good targets.** Start with a deep, unambiguous transit: HAT-P-7 b or TrES-2 b (giant planets, huge signal). Then Kepler-10 b (the first rocky exoplanet, ~0.015% depth — a real test). Ambitious extensions: the TRAPPIST-1 seven-planet system (K2 data, faint star), π Mensae c (TESS).

**Visualizations.**
- Raw multi-year light curve → zoomed to a single transit (the "there it is" moment).
- BLS periodogram with the planet's period spiking out.
- Phase-folded, binned transit with the fitted box/trapezoid model.
- To-scale diagram of star + silhouetted planet derived from our measured radius ratio.

**Watch out for.** Detrending can eat the transit if the window is too aggressive; quarters have instrumental offsets (stitching handles most of it); depth gives the radius *ratio* only — the absolute radius inherits the stellar radius uncertainty.

**Extensions.** Recover all 8 planets of Kepler-90 with iterative BLS masking; look for transit timing variations; try a star with *no* known planet and set an upper limit.

---

## A2. Exoplanet population atlas ★☆☆

**Science question.** What kinds of planets exist, and around what kinds of stars? What are the observational biases shaping what we see?

**Data.** NASA Exoplanet Archive (Caltech/IPAC) — the authoritative table of all confirmed exoplanets (5000+), with masses, radii, periods, discovery methods, and host-star properties. Free TAP/API access.

**Access.**
```python
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
planets = NasaExoplanetArchive.query_criteria(table="pscomppars", select="*")
```

**Method + visualizations.**
- **Radius vs. orbital period**, colored by discovery method — shows hot Jupiters, super-Earths, and *why* each method finds what it finds (transits favor short periods; radial velocity favors massive planets).
- **The radius valley**: histogram of planet radii for short-period planets shows the observed dip near ~1.8 Earth radii (Fulton et al. 2017) separating super-Earths from mini-Neptunes — real, published structure we can reproduce.
- **Mass–radius diagram** with density curves for iron/rock/water/H-He compositions — which worlds are rocky?
- **Discoveries per year, stacked by method** — the history of the field in one chart.
- **Habitable-zone plot**: stellar temperature vs. received stellar flux, with the conservative/optimistic HZ boundaries (Kopparapu et al. 2013), marking Earth, TRAPPIST-1 e/f/g, Proxima b.

**Watch out for.** This table is heterogeneous — parameters come from different papers with different methods. Plot error bars where they matter; never treat detection counts as occurrence rates without saying so (selection bias is the whole story here).

---

# B. Stars & the Milky Way (Gaia)

ESA's Gaia mission measured positions, parallaxes (→ distances), proper motions, and photometry for ~1.8 billion stars — the best map of the galaxy ever made. Data Release 3 is fully public via the ESA archive; `astroquery.gaia` submits ADQL (SQL-like) queries. The key relation: distance in parsecs = 1/parallax in arcseconds.

## B1. The Hertzsprung–Russell diagram ★☆☆

**Science question.** If you plot stellar color against intrinsic brightness for real stars, does the structure of stellar evolution appear — and can we identify each feature?

**Background.** The HR diagram is the single most important plot in stellar astrophysics. Stars are not scattered randomly: they trace the **main sequence** (hydrogen fusion, where stars spend ~90% of their lives), the **red giant branch**, the **red clump** (core helium fusion), and the **white dwarf sequence** (dead cores). With Gaia parallaxes we compute true absolute magnitudes: `M_G = G + 5·log10(parallax_mas) − 10`.

**Access.**
```python
from astroquery.gaia import Gaia
job = Gaia.launch_job_async("""
    SELECT phot_g_mean_mag, bp_rp, parallax
    FROM gaiadr3.gaia_source
    WHERE parallax > 10                -- within 100 pc
      AND parallax_over_error > 10    -- reliable distance
      AND ruwe < 1.4                  -- good astrometric fit
""")
```

**Method.** Query a clean nearby sample → compute M_G → 2D density plot (hexbin or log-scaled histogram; a scatter of 300k points is a blob) of M_G vs. BP−RP color → annotate every feature: main sequence, giant branch, red clump, white dwarfs, the Sun's position (M_G ≈ 4.67, BP−RP ≈ 0.82).

**Visualizations.**
- The annotated HR diagram (this is the keeper — worth polishing).
- Same diagram colored by tangential velocity — fast-moving halo stars trace a distinct, older population.
- Bonus quick win from the same archive: an **all-sky star-count density map** in galactic coordinates — the Milky Way's disk, dust lanes, and the Magellanic Clouds appear from nothing but star counts.

**Watch out for.** Naive distance = 1/parallax is only safe when parallax_over_error is high (our cut ensures this); no reddening correction is fine within ~100 pc, but say so; quality cuts remove faint white dwarfs — note the selection.

**Extensions.** Overlay theoretical isochrones (MIST/PARSEC, both free) to read off stellar ages; plot the white dwarf cooling sequence alone.

## B2. 3D map of the solar neighborhood ★★☆

**Science question.** What does our corner of the galaxy actually look like in three dimensions?

**Method.** Gaia query for stars within 50–100 pc with good parallaxes → convert (RA, Dec, parallax) to XYZ coordinates with `astropy.coordinates` → interactive 3D scatter in plotly, sized by brightness, colored by BP−RP (real color mapping: blue-white hot stars to red dwarfs) → label famous stars (Sirius, Vega, Alpha Centauri, Barnard's Star, Betelgeuse is too far — that itself is a teachable label).

**Visualizations.** A rotatable, zoomable HTML 3D map — the flagship interactive piece. Also: a 2D projection onto the galactic plane; a histogram of stellar types by distance showing that red dwarfs utterly dominate by number (and that even Gaia misses the faintest ones nearby — selection effects made visible).

**Extensions.** Add proper-motion vectors and animate the sky 100,000 years forward/backward (constellations dissolve); mark known exoplanet hosts from A2 in the 3D map — cross-project payoff.

## B3. Hunt star clusters by their motion ★★☆

**Science question.** Stars born together move together. Can we pick a real open cluster out of the field by kinematics alone?

**Method.** Query a few degrees around the Pleiades or Hyades → plot proper motion (pmra vs. pmdec): the cluster appears as a dense knot against the smeared field → select the knot, plot those stars on the sky and in an HR diagram → the cluster's HR diagram is razor-thin (single age!) compared to B1's field-star diagram — visual proof that clusters are single-age populations, and the basis of isochrone dating.

**Visualizations.** Four-panel story: sky field → proper-motion space (knot circled) → members on sky → members' HR diagram vs. field. Extension: fit an isochrone and *measure the cluster's age* (Pleiades ≈ 100 Myr).

---

# C. Compact objects & transients

## C1. Gravitational waves: the GW150914 chirp ★★☆

**Science question.** Can we pull the first gravitational wave ever detected out of the raw LIGO strain data, see the chirp, and hear it?

**Background.** On 2015-09-14, LIGO's two detectors recorded the spacetime ripple from two black holes (~36 and ~29 solar masses) merging into one of ~62 — the missing ~3 solar masses radiated as gravitational waves. The raw strain (fractional length change ~10⁻²¹) is buried in seismic and instrumental noise; band-passing and whitening reveal it. All LIGO/Virgo data is public at the Gravitational Wave Open Science Center (GWOSC).

**Access.**
```python
from gwpy.timeseries import TimeSeries
strain = TimeSeries.fetch_open_data("H1", 1126259446, 1126259478)  # GPS seconds around the event
```

**Method.**
1. Fetch 32 s of strain from both detectors (H1 Hanford, L1 Livingston).
2. Band-pass ~35–350 Hz + notch the power-line harmonics; whiten.
3. Zoom to the final 0.2 s: the chirp — amplitude and frequency both rising as the black holes spiral in.
4. Q-transform spectrogram: the upward-sweeping track.
5. Overlay H1 and L1 (shifted by the ~7 ms light-travel time between sites, one inverted) — the same signal in two detectors 3000 km apart is the point.
6. Write the whitened chirp out as a WAV file, plus a version frequency-shifted up a few hundred Hz (the famous audible "thump").

**Visualizations.** Raw vs. filtered strain; the two-detector overlay; the spectrogram chirp; and audio files.

**Watch out for.** Whitening/filter edge artifacts — trim the segment edges; we're doing signal *visualization*, not parameter estimation — quote published masses rather than deriving them (deriving them properly is a genuinely hard extension).

**Extensions.** GW170817 — the neutron star merger seen simultaneously in gravitational waves and light; the "stellar graveyard": masses of all cataloged mergers (GWTC catalogs, downloadable via the `gwosc` package) plotted as the famous ascending-dots chart; estimate the chirp mass from the frequency sweep rate `f(t)` — one formula, real astrophysics.

## C2. The pulsar P–Ṗ diagram ★☆☆

**Science question.** Plot every known pulsar's spin period against its slow-down rate — does the population sort itself into physically distinct families?

**Background.** Pulsars are rotating neutron stars. Period P and spin-down rate Ṗ yield the characteristic age `τ = P/(2Ṗ)` and surface magnetic field `B ≈ 3.2×10¹⁹ √(P·Ṗ)` gauss. On a log-log P–Ṗ plot: ordinary pulsars cluster in the middle; **millisecond pulsars** (old, "recycled" by accretion, P ~ few ms, tiny Ṗ) sit lower-left; **magnetars** (B ~ 10¹⁴–10¹⁵ G) upper-right; below the "death line," radio emission shuts off.

**Data.** ATNF Pulsar Catalogue (~3500 pulsars), queried with `psrqpy` — small and clean.

**Visualizations.** The P–Ṗ diagram with constant-age and constant-B diagonal lines, death line, binaries marked (they cluster with the millisecond pulsars — that's the recycling evidence); sky map in galactic coordinates (pulsars hug the galactic plane; millisecond pulsars, being old, don't); histogram of spin periods — the fastest known spins ~700 times per second, a city-sized object rotating at a good fraction of breakup.

## C3. Fast radio bursts ★☆☆

**Science question.** What do we actually know about FRBs — millisecond radio flashes, most from other galaxies, cause still debated?

**Data.** CHIME/FRB public catalog (Catalog 1: 536 bursts, 2021; later releases add thousands) — downloadable CSV.

**Method + visualizations.** All-sky map of burst positions (uniform sky ↔ extragalactic origin); dispersion measure (DM) histogram vs. the maximum the Milky Way can contribute along each sightline (the excess is the extragalactic smoking gun — explain DM as the frequency-dependent delay from free electrons en route); repeaters vs. one-offs (duration/bandwidth differences); for the handful of host-localized FRBs, DM vs. redshift — the Macquart relation, which weighs the invisible ionized baryons between galaxies.

---

# D. Cosmology

## D1. Supernova Hubble diagram → dark energy ★★☆

**Science question.** Using published Type Ia supernova distances, is the universe's expansion decelerating (as gravity alone predicts) or accelerating?

**Background.** Type Ia supernovae are standardizable candles: their calibrated peak luminosity is known, so their apparent brightness gives distance (distance modulus `μ = m − M = 5·log10(d_L/10 pc)`). Plot μ against redshift z and compare cosmological models. In 1998 two teams found distant supernovae *fainter* than a coasting universe predicts — expansion is accelerating. Dark energy. Nobel Prize 2011.

**Data.** Pantheon+ compilation — 1701 light curves of ~1550 SNe Ia, z = 0.001 to 2.26 — public on GitHub (`PantheonPlusSH0ES/DataRelease`): plain-text table with redshifts, distance moduli, uncertainties, plus the covariance matrix.

**Method.**
1. Load the table; plot the Hubble diagram (μ vs. z, log-x, error bars).
2. Compute model curves with `astropy.cosmology`: flat ΛCDM (Ωm≈0.3, ΩΛ≈0.7), matter-only Einstein–de Sitter (Ωm=1), and an empty coasting universe.
3. **Residual plot** relative to the coasting model — the data arc *above* it at high z; matter-only arcs below. This one panel is the discovery.
4. χ² fit of (H0, Ωm) in flat ΛCDM; draw confidence contours. Expect H0 ≈ 73 km/s/Mpc from this dataset (SH0ES calibration).
5. Note the **Hubble tension**: this local H0 (~73) vs. Planck's CMB-inferred value (~67.4) disagree at >4σ — one of the liveliest open problems in physics, and our two projects D1 and D3 sit on opposite sides of it.

**Watch out for.** For honest error bars use the provided covariance matrix (a diagonal-only fit is fine for the plot, but say which we did); low-z points calibrate, high-z points discriminate models.

**Extensions.** Recreate Hubble's original 1929 diagram from his published table of 24 galaxies (public domain — and his H0 was ~7× too high, a great story about systematics); fit w (dark energy equation of state).

## D2. The cosmic web from SDSS ★★☆

**Science question.** Map hundreds of thousands of real galaxies in 3D — is matter distributed randomly, or in filaments, walls, and voids?

**Data.** Sloan Digital Sky Survey spectroscopic catalog (millions of galaxy redshifts), via SQL query at SkyServer or `astroquery.sdss`. Free.

**Method.** Select main-sample galaxies in a thin declination slice (a few degrees) with z < ~0.15 → polar "wedge plot": angle = RA, radius = redshift (≈ distance), one dot per galaxy. The cosmic web appears: filaments, voids, and the **Sloan Great Wall** (~z ≈ 0.08, one of the largest known structures). Color dots by galaxy color (u−r or g−r): red ellipticals crowd the dense knots; blue star-formers fill the field — environment shapes galaxies, visible by eye.

**Watch out for.** It's a flux-limited survey: galaxies thin out at high z because faint ones drop out, not because space empties — annotate this on the figure. "Fingers of God" (radial smearing in clusters, from peculiar velocities) will be visible; label one — it's a feature, not a bug.

**Extensions.** Galaxy color bimodality (red sequence / blue cloud in a color–magnitude diagram); compare the wedge with a random catalog of the same selection function to show the clustering is real.

## D3. The cosmic microwave background ★★☆

**Science question.** Render the actual map of the universe at age 380,000 years, and extract the acoustic peaks that fix the universe's composition and geometry.

**Data.** Planck Legacy Archive (ESA): the SMICA CMB temperature map (HEALPix FITS; full resolution runs to hundreds of MB — lower-resolution versions exist at NASA's LAMBDA archive) and the published TT power spectrum (a ~2 MB text file). A power-spectrum-only version of this project skips the big download entirely.

**Method + visualizations.**
- `healpy` Mollweide projection of the temperature map, ±300 µK color scale — the fluctuations are 1 part in 100,000 of the 2.725 K background (worth stating on the figure).
- Published TT power spectrum: Dℓ vs. multipole ℓ, acoustic peaks annotated — first peak at ℓ ≈ 220 (~1° scale) ⇒ spatially flat universe; peak height ratios ⇒ baryon and dark matter densities. We plot real Planck points over the ΛCDM best-fit curve and annotate what each feature *means*.
- Optional: compute the power spectrum ourselves from the map with `healpy.anafast` and overlay it on Planck's published one — they'll agree at low/mid ℓ, and the high-ℓ divergence honestly shows why real analyses need masks and noise modeling.

**Watch out for.** Galactic-plane residuals in the map (mask exists for a reason); `anafast` without a mask is a demonstration, not a measurement — label it as such.

## D4. Dark matter from galaxy rotation curves ★☆☆

**Science question.** Stars orbiting far from a galaxy's center should slow down like outer planets do (v ∝ 1/√r) if visible matter is all there is. Do they?

**Data.** SPARC database (Case Western): high-quality rotation curves of 175 galaxies plus photometry-based visible-mass models. Small text files, free.

**Method + visualizations.** For a handful of galaxies, plot observed rotation velocity vs. radius against the Keplerian expectation from visible matter alone: observed curves stay **flat** far beyond where they should fall — the classic Rubin & Ford result and a pillar of dark matter. Overlay v_gas, v_disk, v_bulge components (provided) vs. observed; the gap is the dark halo. Small data, big physics — arguably the best payoff-per-megabyte in this whole document.

---

# E. Solar system

## E1. Kirkwood gaps in the asteroid belt ★★☆

**Science question.** Do Jupiter's orbital resonances really carve gaps in the asteroid belt — visible in a plain histogram of 1M+ real orbits?

**Data.** Minor Planet Center's `MPCORB.DAT` — orbital elements for every numbered/provisional minor planet (~1.4 million; ~90 MB gzipped, a few hundred MB unpacked). Fixed-width text, format documented on the MPC site. Alternative: JPL Small-Body Database query API.

**Method + visualizations.**
- Parse semi-major axes → histogram over 1.8–5.5 AU → sharp gaps at the mean-motion resonances with Jupiter: 3:1 (2.50 AU), 5:2 (2.82 AU), 7:3 (2.96 AU), 2:1 (3.28 AU) — annotate each. Populations *at* resonances too: Hildas (3:2, ~3.97 AU), Jupiter Trojans (1:1, 5.2 AU).
- Scatter a vs. eccentricity and a vs. inclination: asteroid **families** (Vesta, Eos, Koronis — collision debris) appear as clumps.
- Top-down snapshot of the belt: compute each asteroid's current position from its elements (astropy/rebound) — reveals the non-uniform structure, Trojan clouds leading/trailing Jupiter.

**Watch out for.** Discovery bias (small asteroids are found when near Earth); parsing fixed-width fields carefully (the MPC packed-date format is a rite of passage).

**Extension (labeled simulation).** `rebound` N-body integration of test particles + Jupiter for ~1 Myr: watch the 3:1 gap *form*. Mechanism, demonstrated.

## E2. Solar system orbit visualizer ★★☆

**Science question.** Where is everything, right now, and where is it going?

**Data.** JPL Horizons ephemerides via `astroquery.jplhorizons` — the same system JPL uses for actual mission navigation. Positions/velocities for planets, comets, asteroids, and spacecraft.

**Method + visualizations.** Fetch heliocentric positions over time → interactive 3D plotly figure and/or animated top-down view. Showcase trajectories: the inner planets' clockwork; a long-period comet's plunge; **ʻOumuamua and 2I/Borisov** (hyperbolic interstellar visitors — e > 1, in real coordinates); **Voyager 1 and 2's** grand-tour paths out of the system. Bonus: "night sky tonight" — compute current planet positions in the local sky with astropy.

## E3. Near-Earth object close approaches ★☆☆

**Data.** JPL's CNEOS close-approach API (`ssd-api.jpl.nasa.gov/cad.api`) and NASA NeoWs — free, no key needed for the former.

**Method + visualizations.** Timeline of upcoming close approaches: distance (in lunar distances) vs. date, dot size = estimated diameter; histogram of miss distances vs. size — plenty of small rocks pass inside the Moon's orbit every year, big ones are rare (honest risk communication, with the Torino scale for context); Apophis's well-measured 2029 approach (~38,000 km — inside geostationary orbit, zero impact risk, and we can plot exactly why we know).

---

# F. The Sun

## F1. 275 years of solar cycles ★☆☆

**Science question.** The Sun's activity cycles every ~11 years. What does the longest-running experiment in astrophysics — counting sunspots since before the French Revolution — actually show?

**Data.** SILSO (Royal Observatory of Belgium): official international sunspot number — monthly means since 1749, daily since 1818. Small CSVs, free. For the butterfly diagram: the RGO/USAF-NOAA sunspot-group positions dataset (1874–present), available via NOAA/solar-cycle archives.

**Method + visualizations.**
- Full 1749–present sunspot time series, cycles numbered 1–25; annotate the Dalton and (pre-1749, mention only) Maunder minima, and the current Cycle 25.
- The **butterfly diagram**: sunspot latitude vs. time — each cycle's spots start at mid-latitudes and drift equatorward, drawing butterfly wings. One of the most beautiful plots in science, from a plain scatter of real positions.
- 13-month smoothed series with cycle maxima marked; Cycle 25 vs. the official prediction.

## F2. The Sun today ★★☆

**Data.** SDO (Solar Dynamics Observatory) imagery via `sunpy` — actual near-real-time images of the Sun in extreme UV; GOES X-ray flux for flare monitoring.

**Method + visualizations.** Multi-wavelength grid of today's Sun (AIA 171/193/304 Å — different temperatures show corona, loops, prominences); overlay GOES X-ray flux to catch flares; if a big active region is present, track it across the disk over days (solar rotation, live).

---

# G. Imaging & fun

## G1. Color composites of real galaxies ★★☆

**Data.** Survey cutout services: DESI Legacy Survey viewer API, SDSS cutouts, `astroquery.skyview`; for the ambitious, raw Hubble/JWST public FITS frames from MAST.

**Method.** Fetch calibrated FITS images of M51, M31, M104 etc. in 3+ filters → align (`reproject`) → asinh-stretch RGB composite (`astropy.visualization.make_lupton_rgb` — the published algorithm used for SDSS's own images). Honest imaging: we document filter→channel mapping and stretch, i.e. what "true color" does and doesn't mean.

## G2. Meteorite landings map ★☆☆

**Data.** The Meteoritical Society's database of ~45,000 classified meteorites (via NASA's open data portal, CSV with coordinates, mass, class, year).

**Method + visualizations.** World map of falls vs. finds — instantly shows collection bias (deserts and Antarctica light up, oceans and jungles are blank: nobody finds black rocks in forests); mass distribution; discoveries over time. A relaxed one-evening project with a genuinely interesting selection-effects story.

---

## Output formats

Every project produces at least one **publication-quality static figure** (matplotlib, saved as PNG + SVG into `figures/`). Where interactivity earns its keep (B2's 3D star map, E2's orbits, zoomable light curves), we add a **self-contained interactive HTML** file (plotly) — openable in any browser, shareable as-is. Interactive pages must embed their (downsampled) data inline; keep embedded data well under ~5 MB for smooth loading.

## Practical notes

- **Be polite to archives**: cache every download in `data/` and never re-fetch what we have; use `TOP N` in ADQL while iterating on a query, remove it for the final run.
- **Big queries**: Gaia and SDSS can return millions of rows; we design the cut *first* (sky region, magnitude, quality flags), then fetch once.
- **If we `git init`** this folder: `.gitignore` gets `data/`, `.venv/`, `__pycache__/`. Figures and code are the deliverables; raw data is re-fetchable by script.
- **Each finished project's README** records: the science question, data source + retrieval date, exact quality cuts, headline result (with numbers), and known caveats. That's what makes this a portfolio rather than a pile of plots.

## Status

| Project | Status |
|---------|--------|
| A1 Exoplanet transit | not started |
| A2 Exoplanet atlas | not started |
| B1 HR diagram | not started |
| B2 3D solar neighborhood | not started |
| B3 Cluster hunting | not started |
| C1 GW150914 chirp | not started |
| C2 Pulsar P–Ṗ | not started |
| C3 Fast radio bursts | not started |
| D1 SN Hubble diagram | not started |
| D2 SDSS cosmic web | not started |
| D3 Planck CMB | not started |
| D4 Rotation curves | not started |
| E1 Kirkwood gaps | not started |
| E2 Orbit visualizer | not started |
| E3 NEO approaches | not started |
| F1 Solar cycles | not started |
| F2 The Sun today | not started |
| G1 Galaxy composites | not started |
| G2 Meteorite map | not started |
