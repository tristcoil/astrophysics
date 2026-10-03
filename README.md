# Astrophysics Project Atlas

> **Note: everything in this repository is AI-generated.** The plan, code, analyses, figures and text were all produced by AI (Anthropic's Claude, working in Claude Code) at the direction of a hobbyist. Results are compared with published values throughout, but nothing here has been reviewed by a professional astronomer.

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
| `emcee`        | MCMC sampling for posteriors (orbit fits, cosmology)         |
| `sncosmo`      | SALT supernova light-curve fits (D6)                         |
| `galpy`, `gala`| Galactic potentials, orbits and stream frames (H3, J6, J7)   |
| `spectral-cube`| Radio data cubes — 21-cm hydrogen (H2, J1)                   |
| `ehtim`        | Event Horizon Telescope interferometric imaging (I5, optional) |
| `lenstronomy`  | Strong gravitational lens modelling (J5 extension)           |
| `bilby`, `lalsuite` | Gravitational-wave parameter estimation and waveforms (I7, I9) |
| `qnm`          | Kerr black hole ringdown frequencies (I7)                    |
| `enterprise-pulsar`, `pint-pulsar`, `ceffyl` | Pulsar-timing noise models and gravitational-wave background fits (I8) |
| `pycs3`        | Time delays between lensed quasar images (D9)                |
| `colossus`     | Halo mass functions and cosmology helpers (D11, K2)          |
| `celerite2`    | Fast Gaussian-process fits to light curves (I17)             |
| MESA + MESA SDK | Stellar evolution code (Fortran, not a Python package; installs without root) (B11) |
| `skyfield`     | Positions of the Sun, Moon, planets and stars from JPL ephemerides (E8, B12) |
| `orbitize`     | Orbits of directly imaged planets (A6)                       |
| `dustmaps`     | 3D Galactic dust maps (H6)                                   |

> **Note on Python 3.14:** it is very new; if a compiled package (`healpy`, `gwpy`) has no wheel for it yet, the fallback is a venv on an older interpreter or a miniforge/conda environment. We'll deal with it per-project if it comes up.

---

## Project index

Difficulty: ★☆☆ = an afternoon · ★★☆ = a weekend · ★★★ = multi-session.

| #  | Project                                   | Field          | Difficulty | Download    | Payoff |
|----|-------------------------------------------|----------------|-----------|-------------|--------|
| A1 | Detect an exoplanet transit               | Exoplanets     | ★☆☆       | ~10–50 MB   | Find a planet in raw spacecraft data |
| A2 | Exoplanet population atlas                | Exoplanets     | ★☆☆       | ~5 MB       | The demographics of 5000+ worlds |
| A3 | Radial velocities: weigh a planet by its wobble | Exoplanets     | ★☆☆       | <20 MB      | 51 Peg b, the 2019 Nobel planet |
| A4 | The Earth as an exoplanet (DSCOVR)        | Exoplanets     | ★★☆       | ~0.5–5 GB   | Our planet as a single dot of light |
| A5 | Planets being born (ALMA discs)           | Exoplanets     | ★★☆       | ~1.3 GB     | Rings and gaps carved by young planets |
| A6 | Watch planets orbit another star          | Exoplanets     | ★★☆       | <1 MB       | HR 8799's four planets in motion |
| B1 | Hertzsprung–Russell diagram from Gaia     | Stars          | ★☆☆       | ~50–200 MB  | All of stellar evolution in one plot |
| B2 | 3D map of the solar neighborhood          | Stars          | ★★☆       | ~70 MB      | Fly through real nearby stars |
| B3 | Hunt star clusters by their motion        | Stars          | ★★☆       | ~50 MB      | Clusters pop out of the noise |
| B4 | Cepheids and the Leavitt law              | Stars          | ★☆☆       | ~5–200 MB   | The first step of the distance ladder |
| B5 | Asteroseismology: weigh stars by their sound | Stars          | ★★☆       | ~0.1–1 GB   | Hear red giants ring |
| B6 | How old are the oldest stars?             | Stars          | ★★☆       | ~100–300 MB | Can stars be older than the universe? |
| B7 | Reading starlight: spectroscopy           | Stars/galaxies | ★★☆       | ~50–500 MB  | Spectral classes, redshifts, the BPT diagram |
| B8 | A time-lapse of stellar evolution in clusters | Stars      | ★★☆       | ~0.3–1 GB   | Watch stars age across 7,200 clusters |
| B9 | Eclipsing binaries: weigh and measure stars | Stars        | ★★☆       | <50 MB      | Masses and radii to 2% |
| B10 | White dwarfs: the cooling embers         | Stars          | ★★☆       | ~0.1–0.5 GB | Dead stars freezing into crystals |
| B11 | Run a star to core collapse (MESA, simulation) | Stars    | ★★★       | a few GB (code) | The onion inside a star before it explodes |
| B12 | Constellations in motion                  | Stars          | ★☆☆       | ~0.1 GB     | The sky ±100,000 years, and stars that brush past the Sun |
| C1 | Gravitational waves: the GW150914 chirp   | Compact objects| ★★☆       | ~265 MB     | See + hear two black holes merge |
| C2 | Pulsar P–Ṗ diagram                        | Compact objects| ★☆☆       | <5 MB       | The whole neutron star zoo |
| C3 | Fast radio bursts (CHIME)                 | Compact objects| ★☆☆       | <5 MB       | The newest mystery in astronomy |
| C4 | GW170817: the neutron star merger         | Compact objects | ★★★       | ~0.5–1 GB   | Gravity travels at the speed of light |
| C5 | The black hole population (GWTC-4)        | Compact objects | ★☆☆       | <50 MB      | 218 mergers and the mass gaps |
| C6 | Gamma-ray bursts (Fermi GBM)              | Compact objects | ★☆☆       | <20 MB      | Two kinds of cosmic explosion |
| C7 | Neutrino astronomy: SN 1987A (+ IceCube)  | Compact objects | ★☆☆/★★★   | <1 MB/~100 MB | 24 particles weigh a supernova |
| C8 | What is inside a neutron star?            | Compact objects | ★★☆       | ~0.1–1 GB   | Matter denser than an atomic nucleus |
| C9 | Binary pulsars                            | Compact objects | ★★☆       | <5 MB       | Gravitational waves before LIGO |
| D1 | Supernova Hubble diagram → dark energy    | Cosmology      | ★★☆       | ~80 MB      | Reproduce the 2011 Nobel Prize |
| D2 | The cosmic web from SDSS                  | Cosmology      | ★★☆       | ~50 MB      | Filaments and voids from real galaxies |
| D3 | Planck CMB map + acoustic peaks           | Cosmology      | ★★☆       | 0.3–1 GB*   | Baby picture of the universe |
| D4 | Dark matter from rotation curves (SPARC)  | Cosmology      | ★☆☆       | <5 MB       | Rubin's flat curves, reproduced |
| D5 | Is dark energy changing? (DESI BAO)       | Cosmology      | ★★☆/★★★   | <5 MB/~1–5 GB | The 2025 hint, reproduced |
| D6 | Fresh supernovae onto our Hubble diagram  | Cosmology      | ★★☆       | ~50 MB      | Test D1 on supernovae it never saw |
| D7 | Is the universe really expanding?         | Cosmology      | ★☆☆       | ~0.1–0.5 GB | Three tests that tired light fails |
| D8 | The universe's speedometer: H(z)          | Cosmology      | ★☆☆       | <1 MB       | When the expansion began to speed up |
| D9 | H0 from lensed quasar time delays         | Cosmology      | ★★☆       | <50 MB      | A distance ladder with no rungs |
| D10 | Dark energy in the CMB (ISW effect)      | Cosmology      | ★★☆       | ~0.2–2 GB*  | Dark energy warms photons crossing superclusters |
| D11 | Counting galaxy clusters (eROSITA)       | Cosmology      | ★★☆       | <50 MB      | How clumpy is the universe? |
| E1 | Kirkwood gaps in the asteroid belt        | Solar system   | ★★☆       | ~200 MB     | Resonances carve visible gaps |
| E2 | Solar system orbit visualizer (Horizons)  | Solar system   | ★★☆       | API, small  | Planets, ʻOumuamua, Voyager |
| E3 | Near-Earth object close approaches        | Solar system   | ★☆☆       | API, small  | What's passing by this month |
| E4 | Meteor showers and their parent comets    | Solar system   | ★★☆       | ~0.2–1 GB   | Every shower traced to its comet |
| E5 | Dating surfaces by counting craters       | Solar system   | ★★☆       | ~300 MB     | Ages of the Moon and Mars from craters |
| E6 | Saturn's rings in radio light             | Solar system   | ★★☆       | ~0.1–1 GB   | Ringlets, gaps and the moons that carve them |
| E7 | Jupiter's weather from Hubble maps        | Solar system   | ★★☆       | ~1–3 GB     | Winds and the shrinking Great Red Spot |
| E8 | Total eclipses: draw the Moon's shadow    | Solar system   | ★★☆       | ~0.1–1 GB   | The 2026 and 2027 paths, and Baily's beads |
| E9 | Asteroids in 3D                           | Solar system   | ★★☆       | <0.5 GB     | Shapes, the spin barrier, and DART |
| E10 | Comet tails sculpted by sunlight         | Solar system   | ★★☆       | ~0.5 GB     | Comet Tsuchinshan–ATLAS in SOHO's view |
| F1 | 275 years of solar cycles + butterfly diagram | The Sun    | ★☆☆       | ~25 MB      | The Sun's heartbeat since 1749 |
| F2 | The Sun today, in extreme UV (SDO)        | The Sun        | ★★☆       | ~190 MB     | Actual current images of the Sun |
| F3 | The Gannon storm of May 2024              | The Sun        | ★★☆       | <0.5 GB     | From solar flare to aurora |
| F4 | The Sun's spectrum as one image           | The Sun        | ★☆☆       | ~40 MB      | Thousands of dark lines in one poster |
| G1 | Color composites of real galaxies         | Imaging        | ★★☆       | ~510 MB     | Your own pretty pictures, honestly made |
| G2 | Meteorite landings world map              | Fun/easy       | ★☆☆       | ~10 MB      | Quick, satisfying, geographic |
| H1 | Weigh the black hole at the Galactic Centre | Milky Way      | ★★☆       | <1 MB       | The 2020 Nobel, from S2's orbit |
| H2 | Our own rotation curve from 21-cm hydrogen | Milky Way      | ★★☆       | ~0.1–4 GB   | Our galaxy's rotation and spiral arms |
| H3 | Galactic archaeology with Gaia            | Milky Way      | ★★★       | ~0.2–1 GB   | Stellar streams and the galaxy we ate |
| H4 | The Sun's acceleration from quasars       | Milky Way      | ★★☆       | ~0.1–0.3 GB | Feel the Galaxy's pull in 1.6 million quasars |
| H5 | The Milky Way across the spectrum         | Milky Way      | ★☆☆       | <100 MB     | Ten all-sky maps, radio to gamma rays |
| H6 | Fly through the dust around the Sun       | Milky Way      | ★★☆       | ~3.3 GB     | The Local Bubble and the Radcliffe Wave in 3D |
| I1 | The M–σ relation                          | Black holes    | ★☆☆       | <1 MB       | Black holes grow with their galaxies |
| I2 | Light echoes: reverberation mapping       | Black holes    | ★★☆       | ~10–100 MB  | Weigh a black hole with light-travel time |
| I3 | Quasars across cosmic time                | Black holes    | ★★☆       | ~1–2 GB     | 750,000 feeding black holes |
| I4 | Stellar-mass black holes in X-ray binaries | Black holes    | ★★☆       | <50 MB      | The first black holes ever found |
| I5 | The shadow of M87* from EHT data          | Black holes    | ★★★       | ~50 MB      | The black hole image, from raw visibilities |
| I6 | Ray-trace a black hole (simulation)       | Black holes    | ★★☆       | none        | What the EHT should see |
| I7 | Hawking's area theorem: a black hole rings | Black holes    | ★★★       | ~0.3–1 GB   | The loudest merger tests Hawking |
| I8 | The nanohertz hum (pulsar timing arrays)  | Black holes    | ★★☆/★★★   | 5 MB/~0.6 GB | Supermassive black hole binaries everywhere |
| I9 | Dark sirens: H0 from black hole mergers   | Black holes    | ★★★       | ~1–5 GB     | The expansion rate from 235 mergers |
| I10 | Tidal disruption events                  | Black holes    | ★★☆       | <50 MB      | Black holes caught eating stars |
| I11 | Dormant black holes in Gaia orbits       | Black holes    | ★★☆       | ~0.1–0.3 GB | The nearest black holes, found by a wobble |
| I12 | A lone black hole found by microlensing  | Black holes    | ★★☆       | <10 MB      | Weigh a black hole that emits nothing |
| I13 | An intermediate-mass black hole in ω Cen? | Black holes   | ★★☆       | ~0.5–1 GB   | The missing middle, in seven fast stars |
| I14 | Megamasers: weigh a black hole, measure H0 | Black holes   | ★★★       | <5 MB       | One disk gives a mass and a distance |
| I15 | Quasars as backlights: the end of reionization | Black holes | ★★☆     | ~0.1–1 GB   | When the universe became transparent |
| I16 | Faster than light? Superluminal jets      | Black holes    | ★★☆       | <100 MB     | Jet blobs that seem to outrun light |
| I17 | The flickering of quasars                 | Black holes    | ★☆☆       | <200 MB     | 9,000 quasars on a random walk |
| I18 | Giant radio galaxies (LOFAR)              | Black holes    | ★★☆       | ~0.1–2 GB   | Jets 7 Mpc long |
| J1 | Rotation curves from raw velocity fields  | Galaxies       | ★★☆       | ~50–500 MB  | Measure galaxy rotation ourselves |
| J2 | Galaxy morphology with Galaxy Zoo         | Galaxies       | ★★☆       | ~200 MB     | Hubble's tuning fork, 300,000 galaxies |
| J3 | The galaxy life cycle                     | Galaxies       | ★★☆       | ~0.2–1 GB   | When the universe made its stars |
| J4 | Laniakea: the flows of galaxies           | Galaxies       | ★★☆       | ~20 MB      | Where the Milky Way is falling |
| J5 | Weighing galaxies with Einstein rings     | Galaxies       | ★★☆       | <5 MB       | Lensing weighs dark matter |
| J6 | Dwarf galaxies of the Local Group         | Galaxies       | ★☆☆       | <5 MB       | The darkest galaxies known |
| J7 | Andromeda is coming — or is it?           | Galaxies       | ★★☆       | <1 MB       | A coin flip for our galaxy's fate |
| J8 | The first galaxies with JWST              | Galaxies       | ★★★       | ~0.1–1 GB   | Galaxies 300 Myr after the Big Bang |
| J9 | Fly into the Hubble Ultra Deep Field      | Galaxies       | ★★☆       | ~0.5–2 GB   | 2,221 galaxies placed at their true distances |
| K1 | Why dark matter can't be atoms (BBN)      | Dark matter    | ★★☆       | <1 MB       | The first minutes weigh ordinary matter |
| K2 | Weigh a galaxy cluster three ways         | Dark matter    | ★★☆       | ~0.3 GB     | Zwicky's missing mass, measured |
| K3 | The Bullet Cluster                        | Dark matter    | ★★☆       | ~1 MB/~1 GB | Dark matter that passed straight through |
| K4 | Are black holes the dark matter?          | Dark matter    | ★★☆       | <10 MB      | 16 events rule out a halo of LIGO-sized black holes |
| K5 | Weigh the Milky Way's dark halo           | Dark matter    | ★★★       | ~0.1–1 GB   | How heavy is our galaxy, really? |
| K6 | Hunting dark matter particles underground | Dark matter    | ★★☆       | <5 MB       | WIMP limits, the neutrino fog and DAMA |
| L1 | What makes a supernova shine?             | Supernovae     | ★★☆       | <100 MB     | Radioactive nickel powers the glow |
| L2 | A census of exploding stars (ZTF)         | Supernovae     | ★★☆       | <100 MB     | The supernova zoo, and the missing big stars |
| L3 | Rewind a supernova remnant (Crab, Cas A)  | Supernovae     | ★★☆       | ~0.1–1 GB   | Date the explosions from the debris |
| L4 | Where your atoms were made                | Supernovae     | ★★☆       | ~0.8 GB     | The birthplaces of the periodic table |
| L5 | Runaway survivors of Type Ia supernovae   | Supernovae     | ★★☆       | <100 MB     | The fastest stars in the Galaxy |

\* D3 is small (~2 MB) if we only do the power spectrum; the full-resolution sky map is the heavy part. D10 reuses D3's maps.

**Suggested first arc:** A1 (transit) → B1 (HR diagram) → C1 (gravitational waves) → D1 (dark energy). That sequence ramps difficulty gently and touches four different subfields, four different archives, and four different kinds of visualization.

**Second-wave arcs** (A3–J8, added 2026-09-30):
- *The Hubble tension:* B4 (Cepheids) → C4 (standard siren) → D5 (BAO + DESI) → B6 (oldest stars): four independent routes into the question D1 (H0 = 73.5) and D3 (66.7) left open.
- *Black holes, small to large:* I4 → C5 → H1 → I5 + I6 → I1 → I2 + I3.
- *Galaxies:* J2 → J1 → J5 → J3 → J4 → J7.
- *Quick wins (★☆☆):* A3, B4, C5, C6, C7, I1, J6.

**Third-wave arcs** (I7–I14, D7–D11 and the new section K, added 2026-10-03):
- *Merging black holes:* C1 → I7 (a black hole rings) → C5 → I9 (mergers measure H0) → I8 (supermassive binaries).
- *Finding black holes that don't shine:* I11 (Gaia orbits) → I12 (microlensing) → I13 (ω Cen) → I10 (stars torn apart) → I14 (masers).
- *The expanding universe and dark energy:* D7 (is it expanding?) → D8 (the expansion rate over time) → D9 (H0 from lensing delays) → D10 (dark energy in the CMB) → D11 (dark energy slows structure).
- *The Hubble tension, continued:* I9, I14 and D9 add three routes to H0 that need no Cepheid ladder.
- *Dark matter:* K1 (not atoms) → K4 (not black holes) → K2 + K3 (clusters) → K5 (our own halo) → K6 (catching a particle).
- *Quick wins (★☆☆):* D7, D8.

**Fourth-wave arcs** (B8–B11, C8–C9, H4, I15–I17 and the new section L, added 2026-10-03):
- *A star's life, start to finish:* B9 (weigh stars) → B8 (watch them age) → B11 (inside, up to collapse) → L2 + L1 (the explosion) → L3 (the debris) → C8 or B10 (what remains).
- *Supernovae:* L1 → L2 → L3 → L4 → L5 (see section L).
- *Neutron stars:* C2 → C9 (weigh them) → C8 (what is inside) → L3 (birth kicks) → C4 (mergers) → L4 (gold).
- *Quasars:* I3 → I17 (flicker) → I16 (jets) → I15 (backlights) → H4 (fixed beacons).
- *Quick win (★☆☆):* I17.

**Fifth wave: the visual showcase** (A4–A6, B12, E6–E10, F3–F4, H5–H6, I18, J9, added 2026-10-03). These are chosen for what they look like, each with a "keeper" picture or animation:
- *Posters:* H5 (the Milky Way in ten wavelengths), F4 (the Sun's spectrum), A5 (newborn planetary systems), E6 (Saturn's rings), E7 (a decade of Jupiter).
- *Animations:* A4 (the Earth turning, the Moon crossing it), A6 (planets orbiting another star), B12 (constellations dissolving), E8 (the Moon's shadow sweeping the Earth), E10 (a comet passing the Sun), F3 (a solar storm from Sun to Earth).
- *3D:* H6 (dust clouds around the Sun), J9 (into the Hubble Ultra Deep Field), E9 (asteroid shapes), with B2 and D2 from earlier waves.
- *Quick wins (★☆☆):* B12, F4, H5.

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

## A3. Radial velocities: weigh a planet by its star's wobble ★☆☆

**Science question.** Can we detect a planet we cannot see, from the tiny back-and-forth Doppler shift it causes in its star, and weigh it?

**Background.** A star and its planet both orbit their common centre of mass, so the star's velocity along our line of sight swings with a semi-amplitude K. For 51 Pegasi b, K ≈ 56 m/s (walking pace). Jupiter moves the Sun by 12.5 m/s, and Earth moves it by 0.09 m/s. From the period P, K and the star's mass: `m sin i ≈ K √(1−e²) (P / 2πG)^(1/3) M★^(2/3)`. We only get *m sin i* because the orbit's tilt *i* is unknown, so the result is a minimum mass. 51 Peg b (Mayor & Queloz 1995) was the first planet found around a Sun-like star: Nobel Prize 2019.

**Data.** Two large public velocity archives: the Keck/HIRES LCES release (Butler et al. 2017; 60,949 velocities of 1,624 stars; VizieR `J/AJ/153/208`) and the HARPS RVBank (Trifonov et al. 2020, since updated to ~250,000 velocities of ~5,200 stars; VizieR `J/A+A/636/A74` and the MPIA project page). The RVBank also carries activity indicators for each spectrum. Stellar masses come from the NASA Exoplanet Archive (A2).

**Method + visualizations.**
- Raw velocities of 51 Peg (HD 217014) against time. Then a generalized Lomb–Scargle periodogram, in which the 4.23-day peak stands out, alongside its 1-day alias near 1.3 days.
- Keplerian fit (circular, then eccentric, with a stellar "jitter" term) → P, K, e → m sin i. Phase-folded velocity curve with the model.
- A multi-planet system next (e.g. 55 Cnc from Keck, or HD 10180 from HARPS): fit planets one at a time, subtracting each and re-running the periodogram.
- Hard mode: Proxima Cen b in HARPS (K ≈ 1.4 m/s).
- Stellar activity can fake a planet. Compare the velocity periodogram with the periodogram of the activity index.
- Our planets on A2's mass–period diagram.

**Check against.** 51 Peg b: P = 4.2308 d, K ≈ 56 m/s, m sin i ≈ 0.47 M_Jup.

**Watch out for.** Aliases from nightly sampling. Instrument offsets (e.g. HIRES's 2004 detector upgrade; fit them as free parameters). Jitter must be fitted, or the error bars lie. *m sin i* is a lower limit, so say so on every figure.

**Builds on.** A1 (a transiting planet with radial velocities gives both mass and radius, and so density), A2.

## A4. The Earth as an exoplanet ★★☆

**Science question.** From 1.5 million km away, NASA's DSCOVR spacecraft photographs the whole sunlit Earth every one to two hours. If Earth were a single dot of light around another star, what could we still learn about it: its day length, its oceans and clouds, its seasons? And can we make the movie of our planet turning, and of the Moon crossing its face?

**Background.** DSCOVR orbits the Sun–Earth L1 point, about 1.5 million km sunward of Earth. Its EPIC camera (Earth Polychromatic Imaging Camera) takes 2048 × 2048 images in ten narrow bands from 317 to 780 nm. Telescopes will see Earth-like exoplanets only as unresolved points, so the question is what a point's colour and brightness reveal. Summing each EPIC image into one number per band turns a decade of pictures into an exoplanet-style light curve. Its autocorrelation recovers the 24-hour rotation period, and clouds, continents and oceans modulate it (Jiang et al. 2018). On 16 July 2015 EPIC caught the Moon passing between the spacecraft and Earth, showing the lunar far side, never seen from Earth, in full sunlight.

**Data.** The EPIC archive and API (`epic.gsfc.nasa.gov`; no key needed): natural-colour and enhanced images with metadata (spacecraft, Sun and Moon positions) for 3,639 days since 13 June 2015, typically 13–22 images a day. The calibrated ten-band images (level 1B) are archived at NASA's Atmospheric Science Data Center.

**Method + visualizations.**
- The keeper: a year of Earth in one movie, with an "Earth as a dot" panel beside it showing the dot's colour and brightness.
- Light curves in each band; autocorrelation and periodogram → the rotation period. How few images per day would still recover it?
- Colour against time: is the Pacific-facing hemisphere darker and bluer than the Africa-facing one? How much does cloud cover change the answer?
- Seasons: the North Pole tips into view in June, the South Pole in December.
- The Moon transit of 16 July 2015 as an animation, and the Moon's brightness against Earth's (the Moon reflects ~12% of sunlight, Earth ~30%).

**Check against.** Jiang et al. 2018 (AJ 156, 26): the 24-hour period recovered by autocorrelation.

**Watch out for.** The Sun–Earth–spacecraft angle swings between ~4° and ~15° as DSCOVR loops around L1, which changes the brightness; correct for it. The natural-colour PNGs are processed for display, so use level-1B data for photometry. The archive has gaps, e.g. a safe mode from mid-2019 to early 2020.

**Builds on.** A1 and A2 (exoplanets as points of light), E2 (orbits), F2 (another planet-sized disk, imaged every day).

## A5. Planets being born: ALMA's protoplanetary discs ★★☆

**Science question.** Young stars are surrounded by discs of gas and dust from which planets form. ALMA's sharpest images show these discs carved into rings, gaps and spirals. Are we seeing young planets clearing their paths, and how massive would they have to be?

**Background.** A planet embedded in a disc pushes gas away from its orbit and opens a gap, and dust grains collect in rings where the gas pressure peaks. A gap's width and depth grow with the planet's mass relative to the star, and depend on the disc's viscosity and temperature. In 2014 ALMA's first long-baseline image showed concentric rings around HL Tau, a star about a million years old (ALMA Partnership 2015). DSHARP (Disk Substructures at High Angular Resolution Project; Andrews et al. 2018) imaged 20 nearby discs at 35 mas (~5 au) resolution in 1.25 mm continuum and CO. Rings and gaps are everywhere; a few discs have spirals. From the gap widths, Zhang et al. (2018) inferred planets from about Neptune's to Jupiter's mass.

**Data.** The DSHARP data release (`almascience.eso.org/almadata/lp/DSHARP/`): fiducial continuum images of all 20 discs as FITS (~64 MB each), CO line cubes (~2 GB each, optional) and the calibrated measurement sets. HL Tau's ALMA long-baseline images (ALMA science verification data). Gaia distances to each star.

**Method + visualizations.**
- The keeper: a gallery of the 20 discs on one physical scale, with the orbits of the Solar System's planets drawn for comparison.
- Deproject each disc (inclination, position angle) and make a radial brightness profile; find the gaps and rings.
- Gap widths → planet masses with Zhang et al.'s scaling; plot these unseen planets on A2's mass–separation diagram of known exoplanets.
- CO velocity map of one disc (e.g. HD 163296): Keplerian rotation → the star's mass, as with J1's spider diagrams. Look for the "kink" in the CO channels that Pinte et al. (2018) attributed to a planet.
- Spiral arms (Elias 2-27, WaOph 6): their pitch angles.

**Check against.** Huang et al. 2018 (DSHARP II): ring and gap positions; Zhang et al. 2018 (DSHARP VII): planet masses.

**Watch out for.** A gap does not prove a planet: snow lines, magnetic effects and dust growth can also make rings. Gap-to-mass conversions depend on the assumed viscosity and temperature, by factors of several. Continuum images trace millimetre-sized dust, not gas. The images are reconstructions from interferometer data (I5's caveat).

**Builds on.** A2 (exoplanet demographics), I5 (interferometry), J1 (velocity fields), B1 (young stars on the HR diagram).

## A6. Watch planets orbit another star ★★☆

**Science question.** Four giant planets circle the young star HR 8799. Astronomers have imaged them directly since 2008, and found them in archival Hubble images from 1998, so we can watch them move along their orbits. Can we fit those orbits ourselves, make the movie, and test whether such a crowded system can be stable?

**Background.** Direct imaging suppresses the star's light (with a coronagraph and image processing) to reveal planets tens of thousands of times fainter right beside it. HR 8799, 41 pc away and ~40 Myr old, has four planets of several Jupiter masses at ~15–70 au. The outermost takes about 450 years per orbit and the innermost about 45. Over a decade their motion is plain to see. Measuring each planet's position relative to the star (astrometry, to ~1 mas with the Gemini Planet Imager) constrains Keplerian orbits. Wang et al. (2018) found that the system is stable only if the four orbits are nearly coplanar and close to a 1:2:4:8 chain of period resonances. β Pictoris b, on a ~20-year orbit seen nearly edge-on, is a second example; its data ship with the `orbitize!` package.

**Data.** Published astrometry of the four HR 8799 planets compiled by Wang et al. (2018, AJ 156, 192; GPI points with ~1 mas errors), plus earlier Keck measurements (e.g. Konopacky et al. 2016). β Pic b astrometry from `orbitize!`'s example data. Gaia DR3 for the star's distance, and the Hipparcos–Gaia Catalog of Accelerations (Brandt 2021) for a dynamical mass constraint.

**Method + visualizations.**
- Fit each planet's orbit with `orbitize!` (the OFTI method for short arcs, MCMC otherwise) → posterior orbits, periods and inclinations.
- The keeper: an animation of the four planets moving along their fitted orbits from 1998 to 2050, with the measured positions as dots.
- Test coplanarity and the 1:2:4:8 chain; run posterior draws through N-body integrations (REBOUND, as in E1) — which survive 40 Myr?
- β Pic b's edge-on orbit, and its closest passages in front of the star.

**Check against.** Wang et al. 2018: orbit fits and stability results.

**Watch out for.** Short arcs leave orbital elements strongly correlated: quote posteriors, not single best fits. Different instruments have different plate scales and north angles; fit offsets between them. Planet masses come from evolutionary models and the age unless the dynamics constrain them.

**Builds on.** A3 (orbits from motion), E1 (REBOUND), I11 (astrometric orbits), A5 (planets still forming).

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

*B4–B7 go beyond Gaia: pulsating stars, starquakes, the ages of the oldest clusters, and the spectra that hold most of what we know about stars.*

## B4. Cepheids and the Leavitt law ★☆☆

**Science question.** Some pulsating stars reveal their true brightness through their period. Can we reproduce Henrietta Leavitt's 1912 discovery with ~9,600 stars and use it to measure the distance to another galaxy?

**Background.** Classical Cepheids are yellow supergiants that pulsate with periods of 1–100 days. The longer the period, the more luminous the star: the **Leavitt law**, or period–luminosity relation. Leavitt found it among 25 variables in the Small Magellanic Cloud. They all lie at about the same distance, so their apparent brightness tracks their true brightness. Calibrate the law's zero point with nearby Cepheids of known parallax and it measures distance. This is the first step of the ladder behind D1's H0 = 73.5. Dust dims and reddens stars. The **Wesenheit index** `W_I = I − 1.55 (V − I)` cancels that effect by construction and shrinks the scatter dramatically.

**Data.** The OGLE-IV Collection of Variable Stars, classical Cepheids of the Magellanic System (Soszyński et al. 2015–2019: 9,649 stars). It provides periods, mean I and V magnitudes, pulsation modes, and full I/V light curves, via `ftp.astrouw.edu.pl/ogle/ogle4/OCVS/` (`lmc/cep`, `smc/cep`). Milky Way Cepheids with Gaia DR3 parallaxes (`gaiadr3.vari_cepheid` joined to `gaia_source`) set the zero point. Leavitt's original 25 stars (Harvard College Observatory Circular 173, 1912) are public domain.

**Method + visualizations.**
- A few folded light curves (the classic sawtooth rise and slow decline).
- Leavitt's 25 stars next to OGLE's thousands: the same law, nearly 400× the sample.
- Period–luminosity in V, I and W for the LMC. Fundamental-mode and first-overtone pulsators form two parallel sequences. The scatter drops from ~0.2 mag in V to ~0.08 mag in W.
- SMC vs LMC offset → relative distance. The SMC's larger scatter shows its depth along the line of sight.
- Residuals from the law mapped across the LMC show which side of its tilted disk is nearer.
- Zero point from Gaia parallaxes of Milky Way Cepheids → absolute LMC distance.

**Check against.** LMC distance 49.59 ± 0.54 kpc (distance modulus 18.477; Pietrzyński et al. 2019, from eclipsing binaries, a method independent of ours). Wesenheit slope ≈ −3.3 (Soszyński et al. 2015).

**Watch out for.** Separate the pulsation modes before fitting. Gaia parallaxes have a zero-point offset (Lindegren et al. 2021 recipe; bright Cepheids may need more, as Riess et al. 2021 found). Milky Way Cepheids are more metal-rich than LMC ones, so state the metallicity term we assume. Blending in crowded fields makes Cepheids look too bright.

**Builds on.** D1 (this is its first ladder step), C7 (the LMC is SN 1987A's distance), J4 (the zero point of every galaxy distance).

**Extensions.** RR Lyrae stars (also in OGLE) as a second standard candle for old populations; Type II Cepheids.

## B5. Asteroseismology: weigh stars by their sound ★★☆

**Science question.** Stars ring like bells. Can we pick up the oscillations of red giants in Kepler photometry, and use them to weigh and measure the stars?

**Background.** Turbulent convection excites sound waves that make the whole star pulsate. Its power spectrum shows a comb of peaks centred on a frequency **νmax**, with a regular spacing **Δν** (the "large separation"). Δν scales with the square root of the mean density; νmax scales with surface gravity. With solar reference values (νmax☉ = 3090 µHz, Δν☉ = 135.1 µHz, Teff☉ = 5777 K):
`R/R☉ ≈ (νmax/νmax☉)(Δν/Δν☉)⁻²(Teff/Teff☉)^½` and `M/M☉ ≈ (νmax/νmax☉)³(Δν/Δν☉)⁻⁴(Teff/Teff☉)^1.5`.
Red giants oscillate over hours (νmax ~ 10–200 µHz), within reach of Kepler's 30-minute cadence (Nyquist 283 µHz).

**Data.** Kepler long-cadence light curves via `lightkurve` (the A1 pipeline). Comparison catalogues: Yu et al. (2018), seismic parameters for 16,000 Kepler red giants; APOKASC-3 (Pinsonneault et al. 2025), seismic masses, radii and evolutionary states with APOGEE spectroscopy.

**Method + visualizations.**
- Stitch ~4 years of photometry and high-pass filter it. Power density spectrum → fit the granulation background plus a Gaussian envelope → νmax.
- Δν from the autocorrelation of the spectrum. **Échelle diagram**: the spectrum cut into Δν-wide strips and stacked, in which the l = 0, 1, 2 modes line up as ridges. This is the keeper.
- Scaling relations → mass and radius for ~100 giants. Compare with the catalogues. The Δν ∝ νmax^0.75 relation across the sample.
- Seismic stars placed on B1's HR diagram.

**Check against.** Yu et al. (2018) and APOKASC-3 values of νmax, Δν, M and R, star by star.

**Watch out for.** The scaling relations need percent-level corrections (APOKASC-3 applies them; say whether we do). Long cadence cannot see νmax above ~250 µHz. Frequency resolution needs months of continuous data. Granulation noise dominates at low frequency.

**Builds on.** A1 (same archive and tools), B1 (the red clump, seen from the inside).

**Extensions.** Red clump vs red giant branch from mixed-mode period spacings (Bedding et al. 2011): stars of identical surface appearance, burning helium or not. Solar-like oscillations of main-sequence stars in Kepler short cadence (e.g. 16 Cyg A/B).

## B6. How old are the oldest stars? ★★☆

**Science question.** Globular clusters hold the Milky Way's oldest stars. How old are they? And are they younger than the universe that D1 and D3 measured?

**Background.** A cluster's stars formed together. As a cluster ages, the **main-sequence turnoff** (where stars exhaust core hydrogen) moves to fainter, redder stars. Fitting theoretical isochrones to the turnoff reads off the age. In the 1990s globular cluster ages (~15 Gyr) seemed older than a universe with a high H0 allowed: the "age crisis". It was resolved by Hipparcos distances and dark energy. Our own projects restage it. D1's supernovae give 12.5 Gyr with H0 = 73.5; D3's CMB fit gives 13.825 ± 0.033 Gyr.

**Data.** The HST UV Legacy Survey of Galactic Globular Clusters (HUGS; Piotto et al. 2015, Nardiello et al. 2018): astro-photometric catalogues of 56 globular clusters in five bands (F275W, F336W, F438W, F606W, F814W) with membership probabilities, a MAST High Level Science Product. MIST isochrones (Choi et al. 2016), with HST bolometric corrections. Distances from Baumgardt & Vasiliev (2021, Gaia + HST); metallicity and reddening from the Harris catalogue (2010 edition).

**Method + visualizations.**
- Clean colour–magnitude diagrams (quality flags, membership) of M92, NGC 6397, M15 and 47 Tuc.
- Fix distance and reddening from independent sources, then fit the turnoff and subgiant branch on an age grid. Monte Carlo over distance, reddening, metallicity and [α/Fe] → an age distribution for each cluster.
- Age vs metallicity: the metal-rich 47 Tuc is younger than the metal-poor clusters.
- The ultraviolet colours split each cluster into multiple populations (the "chromosome map"). The survey was designed for this.
- The ages on a timeline against D1's and D3's ages of the universe.

**Check against.** M92: 13.80 ± 0.75 Gyr (Ying et al. 2023).

**Watch out for.** A 0.1 mag distance error shifts the age by roughly 1 Gyr. Distance is the dominant term, so propagate it honestly. Isochrone physics (convection, diffusion, helium) differs between model sets, so repeat the fit with a second set (BaSTI or PARSEC) and quote the spread. Stars need a few hundred million years to form after the Big Bang, so an age of *nearly* 13.8 Gyr is not a contradiction.

**Builds on.** B1, B3 (the isochrone-fitting extension), D1, D3.

## B7. Reading starlight: spectroscopy ★★☆

**Science question.** Almost everything we know about the composition, temperature and motion of stars and galaxies comes from spectra. What can we read directly from real SDSS spectra?

**Background.** Cool gas in front of a hot source absorbs at specific wavelengths; hot, thin gas emits at them. Stellar spectral classes (O B A F G K M) are a temperature sequence. Hydrogen's Balmer lines peak in A stars near 9,500 K; metal lines and then molecular bands (TiO) take over in cooler stars. For galaxies, emission-line ratios such as [O III]/Hβ vs [N II]/Hα (the **BPT diagram**; Baldwin, Phillips & Terlevich 1981) tell whether the gas is lit by young stars or by an accreting black hole.

**Data.** SDSS DR17 optical spectra (FITS per spectrum, via `astroquery.sdss` or the SAS). SEGUE stellar parameters (SSPP: Teff, log g, [Fe/H]). The MPA-JHU emission-line measurements for SDSS galaxies (tables `galSpecLine`, `galSpecInfo` via SkyServer/CasJobs). The DR16 quasar catalogue for quasar spectra.

**Method + visualizations.**
- **The stellar sequence**: SDSS stars from hot to cool, normalized and stacked, with the Balmer lines, Ca II H&K, Mg b, Na D and TiO bands labelled. Balmer line strength vs Teff, peaking near 9,500 K. The Boltzmann and Saha equations explain why.
- **Measure redshifts ourselves**: cross-correlate a galaxy spectrum with a rest-frame template. For a quasar, identify Lyα, C IV, Mg II and Hβ by hand, compute z, and compare with the SDSS pipeline.
- **The BPT diagram** for all galaxies with all four lines detected (S/N > 3), with the Kewley et al. (2001) and Kauffmann et al. (2003) dividing lines. Star-forming, composite, Seyfert and LINER classes, and how the AGN fraction climbs with stellar mass.
- The 4000 Å break vs colour: old stars announce themselves.

**Check against.** Our redshifts vs the SDSS pipeline's; BPT class fractions vs Kauffmann et al. (2003).

**Watch out for.** The 3″ fibre sees only a galaxy's centre. Hβ emission sits inside stellar Hβ absorption (use the MPA-JHU corrected fluxes). The S/N cuts select for strong-line galaxies. Telluric (atmospheric) absorption bands are not stellar.

**Builds on.** D2 (same survey). Feeds I1 (velocity dispersions), I3 (quasar spectra), J3 (star-formation rates), J5 (lens spectra).

*B8–B11 follow a star's life: how massive stars are, how they age, what they leave behind, and what goes on inside them up to the moment of collapse. Section L picks up the story at the explosion.*

## B8. A time-lapse of stellar evolution from star clusters ★★☆

**Science question.** A star's life takes millions to billions of years, far too slow to watch. But a star cluster is a snapshot of stars that share one age, so lining clusters up by age makes a time-lapse. Can we watch the main sequence burn down, giants appear and white dwarfs pile up, using thousands of Gaia clusters?

**Background.** A star's main-sequence lifetime falls steeply with mass, roughly `t ≈ 10 Gyr × (M/M☉)^(−2.5)`: a 10 M☉ star lasts ~30 Myr, the Sun ~10 Gyr. So the older a cluster, the lower the mass at its **main-sequence turnoff**. Stars just past the turnoff cross the Hertzsprung gap quickly, climb the red giant branch, settle in the red clump (helium burning), and end as white dwarfs (below ~8 M☉) or in supernovae (above). Hunt & Reffert (2023) found 7,200 clusters in Gaia DR3 (4,114 of them highly reliable), with member lists, ages, distances and extinctions. Theoretical isochrones (MIST, PARSEC, BaSTI) predict the colour–magnitude diagram at each age.

**Data.** The Hunt & Reffert 2023 catalogue and membership lists (VizieR `J/A+A/673/A114`), with Gaia DR3 photometry of the members. MIST isochrones and tracks, now at `mist.science` (the old `waps.cfa.harvard.edu/MIST` site no longer responds), and BaSTI as a second model set. Well-studied anchor clusters: the Pleiades (B3), the Hyades, Praesepe, M67 and NGC 188.

**Method + visualizations.**
- Dereddened colour–magnitude diagrams of clusters binned by age, from ~10 Myr to several Gyr, as an animation: the turnoff slides down the main sequence, giants appear, the red clump forms.
- Turnoff mass against cluster age from isochrone fits → the mass–lifetime relation, compared with t ∝ M^(−2.5).
- Our own isochrone fits (age, distance, extinction, metallicity) for a dozen clusters, compared with Hunt & Reffert's neural-network ages.
- The Hertzsprung gap: the ratio of stars crossing it to stars in the red clump gives the ratio of time spent in each phase.
- White dwarfs in the nearest clusters (Hyades, Praesepe): the youngest clusters host the most massive white dwarfs, because their progenitors were heavier → the initial–final mass relation (B10).

**Check against.** Hunt & Reffert 2023 ages; literature ages of the Pleiades, Hyades, M67 and NGC 188; the initial–final mass relation of Cummings et al. (2018).

**Watch out for.** Differential reddening smears young clusters' diagrams. Unresolved binaries form a second sequence above the main sequence (B1's twin sequence). Membership lists include some field stars. Isochrone sets disagree most after the main sequence, and rotation broadens young clusters' turnoffs. Hunt & Reffert warn that a few thousand of their "clusters" are more likely unbound moving groups.

**Builds on.** B1 (HR diagram), B3 (cluster membership), B5 (giants from the inside), B6 (globular cluster ages).

## B9. Weigh and measure stars: eclipsing binaries ★★☆

**Science question.** Most stellar masses and radii are inferred from models. In eclipsing binaries we measure them directly: the eclipses give the radii, and the orbital motion gives the masses. Can we fit one system ourselves, then use hundreds of others to draw the mass–luminosity and mass–radius relations that decide how long stars live?

**Background.** When two stars eclipse each other, the eclipse depths and durations give their radii relative to the orbit and their surface-brightness ratio. Radial velocities of both stars give semi-amplitudes K₁ and K₂. With the inclination from the light curve (near 90°), Kepler's third law gives each mass: `M₁ sin³i = P K₂ (K₁ + K₂)² (1 − e²)^(3/2) / (2πG)`, and the same for M₂ with K₁. Masses and radii come out to 1–2%. Luminosity rises steeply with mass (roughly L ∝ M^3.5–4 for Sun-like stars), so the main-sequence lifetime, t ∝ M/L, falls as about M^(−2.5). That is why B8's clusters lose their massive stars first. DEBCat (Southworth 2015) lists 388 detached eclipsing binaries with masses and radii measured to ~2%.

**Data.** DEBCat (`astro.keele.ac.uk/jkt/debcat/`): a plain-text table of 388 systems, last updated September 2026. TESS light curves via `lightkurve` (A1's pipeline). Radial velocities for the chosen system from its published paper (VizieR). Gaia DR3's ~2.2 million eclipsing-binary candidates (`gaiadr3.vari_eclipsing_binary`) for the population view.

**Method + visualizations.**
- One system end to end: a TESS light curve fitted with a simple model of our own (two limb-darkened spheres), plus the published radial velocities → both masses and radii, compared with DEBCat.
- Mass–luminosity and mass–radius diagrams from DEBCat, with MIST main-sequence models at several ages and the Sun marked.
- Main-sequence lifetime from M/L → the t ∝ M^(−2.5) relation that sets B8's turnoffs.
- M dwarfs that are larger than models predict (the "radius inflation" problem, a few percent to ~10%).
- The sky map and period distribution of Gaia's eclipsing-binary candidates.

**Check against.** DEBCat's values for the fitted system; the mass–luminosity relations of Eker et al. (2015, 2018).

**Watch out for.** Light from a third star in TESS's 21″ pixels makes eclipses shallower. Starspots distort cool stars' light curves. Eccentric orbits need e and ω from eclipse timing and velocities together. DEBCat favours bright, well-separated systems.

**Builds on.** A1 (TESS light curves), A3 (radial velocities), B1 (HR diagram), B8 (lifetimes).

## B10. White dwarfs: the cooling embers of dead stars ★★☆

**Science question.** About 97% of stars, the Sun included, will end as white dwarfs: Earth-sized cores about as massive as the Sun, held up by electron degeneracy pressure, cooling for billions of years. Can we see the physics in Gaia's ~359,000 white dwarfs: the mass limit, the cooling clock, and the moment their cores freeze into crystals?

**Background.** A white dwarf has no fuel left and shines only by cooling. Degeneracy pressure makes heavier white dwarfs smaller (R ∝ M^(−1/3) at low mass), up to the **Chandrasekhar limit** of ~1.4 M☉. Above it a white dwarf cannot support itself; near it, one in a binary can explode as a Type Ia supernova (L5). Most white dwarfs weigh ~0.6 M☉. Cooling is slow and calculable, so the faintest white dwarfs date the Galactic disk (Winget et al. 1987). As the carbon–oxygen core cools it crystallizes, releasing latent heat that delays cooling by up to ~1 Gyr. Gaia revealed the resulting pile-up (Tremblay et al. 2019). Massive white dwarfs show an even longer delay, the "Q branch" (Cheng, Cummings & Ménard 2019), attributed to ²²Ne sinking in the crystallizing core (Blouin et al. 2021).

**Data.** Gentile Fusillo et al. 2021: the Gaia EDR3 white dwarf catalogue, with ~359,000 high-confidence candidates and temperatures, surface gravities and masses from model fits (VizieR `J/MNRAS/508/3877`). SDSS and DESI spectra to separate hydrogen (DA) and helium (DB) atmospheres. Sirius B as the benchmark: 1.018 ± 0.011 M☉ from its orbit (Bond et al. 2017).

**Method + visualizations.**
- The white dwarf sequence in Gaia's HR diagram, with cooling tracks for 0.6, 1.0 and 1.2 M☉ and the onset of crystallization marked; the pile-up and the Q branch.
- The mass distribution: a peak near 0.6 M☉ and nothing above ~1.4 M☉.
- Mass–radius relation from the photometric fits (radius from luminosity and temperature), over Chandrasekhar's curve, with Sirius B on it.
- The luminosity function of a volume-limited sample (within 100 pc): its faint-end cutoff → the age of the local disk.
- Cooling age against kinematics: older white dwarfs move faster (the disk heats its stars over time).

**Check against.** Tremblay et al. 2019 (position of the pile-up); Gentile Fusillo et al. 2021 (mass distribution); Sirius B (Bond et al. 2017).

**Watch out for.** Photometric masses depend on the assumed atmosphere (hydrogen or helium); assume the wrong one and the mass is wrong. Unresolved double white dwarfs look overluminous and therefore too light. Completeness drops for the faintest, coolest white dwarfs, exactly where the age signal sits. Turning the cutoff into a disk age needs an assumed star-formation history; state it.

**Builds on.** B1 (the white dwarf sequence), B8 (initial–final masses), B6 (stellar ages), L5 (white dwarfs that explode).

## B11. Run a star to core collapse with MESA (labeled simulation) ★★★

**Science question.** What happens inside a massive star in its final million years, and why does it end with an iron core collapsing in under a second? Can we evolve stars of 1, 12 and 20 M☉ ourselves with a professional stellar-evolution code, and compare each stage with real stars?

**Background.** Stars fuse hydrogen to helium, then helium to carbon and oxygen. A star above ~8 M☉ goes on to burn carbon, neon, oxygen and silicon, each stage faster than the last: hydrogen burning lasts millions of years, silicon burning only days (Woosley, Heger & Weaver 2002). The result is an onion of shells of different composition around an iron core. Fusing iron releases no energy, so when the iron core reaches ~1.3–1.8 M☉ it collapses into a neutron star; section L follows what happens next. MESA (Modules for Experiments in Stellar Astrophysics; Paxton et al. 2011–2019; Jermyn et al. 2023) is free and open-source, and its test suite includes ready-made runs from before the main sequence to core collapse (12 and 20 M☉), to a white dwarf (1 M☉), and a Type II-P explosion with its light curve. **Per ground rule 1 this is a simulation**: it is labelled as such, and every stage is compared with observations.

**Data.** Nothing to download but the code: MESA 26.4.1 (April 2026) and the MESA SDK (compilers and libraries for Linux, installed without root), a few GB in all. For comparison: B8's clusters, B1's HR diagram, the luminosities of red supergiants such as Betelgeuse, L2's supernova progenitors and L1's light curves.

**Method + visualizations.**
- Run the test cases `1M_pre_ms_to_wd`, `12M_pre_ms_to_core_collapse` and `20M_pre_ms_to_core_collapse` (hours to a day of CPU each).
- Evolutionary tracks over B1's Gaia HR diagram: the Sun's future path through the red giant tip to a white dwarf; massive stars crossing to red supergiants.
- **Kippenhahn diagrams**: convective zones and burning regions against log(time left before collapse), showing each burning stage shorter than the last.
- The onion: composition against mass coordinate at collapse, and the central temperature–density track toward degeneracy.
- Neutrino losses overtaking light in the late stages, which is why the last stages are so short.
- Extension: the `ccsn_IIp` test case explodes a model and computes its light curve; compare its plateau and nickel tail with L1's observed supernovae.

**Check against.** Figures in the MESA instrument papers; iron-core masses at collapse of ~1.3–1.8 M☉ (e.g. Sukhbold et al. 2016); B8's turnoffs; Betelgeuse's luminosity.

**Watch out for.** Results depend on choices: convective overshoot, mass-loss rates, the size of the nuclear network and the resolution. Run each case with two settings and show the spread. One-dimensional models cannot capture the 3D convection of the final shells. The explosion in `ccsn_IIp` is put in by hand as an energy injection, not computed from neutrino physics.

**Builds on.** B1 and B8 (tracks against real stars), B10 (the white dwarf ending), L1 and L2 (the explosion), C7 (neutrinos).

## B12. Constellations in motion, and stars that brush past the Sun ★☆☆

**Science question.** The constellations look fixed, but every star is moving. Can we run the sky forward and backward a hundred thousand years and watch familiar shapes dissolve, and find the stars that have passed, or will pass, closest to the Sun?

**Background.** A proper motion of a fraction of an arcsecond per year adds up to degrees over 100,000 years, so constellations change shape. In the Big Dipper, five stars (Merak, Phecda, Megrez, Alioth and Mizar) belong to the Ursa Major moving group and travel together, while Dubhe and Alkaid move differently, so the dipper bends. Distance and radial velocity give full 3D motions, and so close approaches. Of 33 million Gaia DR3 stars with full 3D motions, 61 pass within 1 pc of the Sun within a few million years (42 of them with high confidence). The closest is the K dwarf Gliese 710, which will pass 0.064 pc (~13,000 au) away in 1.3 Myr, through the Oort cloud (Bailer-Jones 2022). Scholz's Star passed about 0.25 pc away ~70,000 years ago (Mamajek et al. 2015).

**Data.** B2's Gaia DR3 + Hipparcos sample (Hipparcos for the brightest stars, which saturate Gaia). The Bailer-Jones (2022) encounter table (VizieR `J/ApJ/935/L9`). Constellation lines from Stellarium's open data. Radial velocities from Gaia DR3 and SIMBAD.

**Method + visualizations.**
- The keeper: animations of the Big Dipper, Orion and the Southern Cross from −100,000 to +100,000 years, using full 3D motions rather than straight-line proper motions.
- The whole night sky 1 Myr ago, today, and 1 Myr from now.
- An encounter timeline: distance from the Sun against time for the 61 stars; Gliese 710's path through the Oort cloud.
- The Sun's nearest neighbour through time: which star was closest at each epoch?

**Check against.** Bailer-Jones 2022: Gliese 710 at 0.0636 pc (90%: 0.0595–0.0678 pc) in 1.3 Myr.

**Watch out for.** Straight-line motion fails over long spans: integrate the 3D orbits, in a Galactic potential beyond ~1 Myr. Radial velocities for bright stars come from SIMBAD and vary in quality. The orbital motion of binary stars contaminates their proper motions, the main source of spurious encounters.

**Builds on.** B2 (the 3D neighbourhood), H3 (orbits in a potential), E4 (comets from the Oort cloud).

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

## C4. GW170817: the neutron star merger ★★★

**Science question.** On 17 August 2017 LIGO and Virgo recorded two neutron stars merging. 1.7 seconds later, Fermi saw a gamma-ray burst from the same patch of sky. What can we measure by combining the two?

**Background.** Neutron stars are lighter than GW150914's black holes, so the inspiral is slower. The signal sweeps through LIGO's band for ~100 s (about 3,000 cycles) instead of 0.2 s, which pins down the chirp mass very precisely. The source was 40 (+8/−14) Mpc away, in galaxy NGC 4993, identified from its kilonova glow. Three results followed:
- The gamma rays arrived 1.74 ± 0.05 s after the merger. After travelling 130 million years, this bounds the speed of gravity to within −3×10⁻¹⁵ and +7×10⁻¹⁶ of the speed of light.
- The signal amplitude gives distance directly, with no distance ladder. Combined with the host's recession velocity, it yields H0 = 70 (+12/−8) km/s/Mpc: a **standard siren**.
- Short gamma-ray bursts come from neutron star mergers.

**Data.** GWOSC: 4096 s of H1, L1 and V1 strain around GPS 1187008882.4. For L1, use the **glitch-subtracted** frame (`L-L1_CLEANED_HOFT_C02_T1700406_v3`): a loud instrumental glitch sat on top of the signal. GWOSC also provides the GWTC-1 parameter-estimation samples for GW170817 (distance, inclination). Fermi GBM data for trigger `bn170817529` (TTE/CTIME files, HEASARC FTP). The Hubble-flow velocity of NGC 4993 is from the standard-siren paper (Abbott et al. 2017, *Nature* 551, 85).

**Method + visualizations.**
1. Q-transform of each detector: a minute-long rising track. Show L1 before and after glitch removal.
2. Matched filter with our own post-Newtonian (TaylorF2) template bank, extending C1's code → chirp mass and merger time in each detector.
3. GBM light curve (NaI detectors, 50–300 keV) → burst onset relative to merger. Turn the 1.74 s delay into a Δv/c bound, assuming, as the paper did, an intrinsic emission delay of 0–10 s and a conservative 26 Mpc distance.
4. Standard siren: GWOSC distance posterior × NGC 4993's velocity → H0 posterior, drawn over D1's 73.5 and D3's 66.7 bands.

**Check against.** Network SNR 32.4; detector-frame chirp mass ≈ 1.198 M☉; delay 1.74 ± 0.05 s; H0 = 70.0 (+12.0/−8.0).

**Watch out for.** Virgo's weak signal (SNR ~2) is real information: the source sat near Virgo's blind spot. A low-order template gives only an approximate chirp mass; the published value needs 3.5PN phasing. The H0 posterior is wide because distance and inclination are degenerate. NGC 4993's own peculiar velocity (~300 km/s) is ~10% of its recession velocity.

**Builds on.** C1 (whitening, matched filter), C6 (gamma-ray bursts), D1 and D3 (the Hubble tension).

## C5. The black hole population ★☆☆

**Science question.** After ~218 detected mergers, what does the population of merging black holes look like, and where are the gaps?

**Background.** Each gravitational-wave event yields two masses, an effective spin (χ_eff) and a distance. The population has structure:
- a peak near 10 M☉ and a bump near 35 M☉;
- a predicted **pair-instability gap** (roughly 50–130 M☉), where massive stars should blow themselves apart completely and leave nothing behind;
- a lower mass gap between the heaviest neutron stars and the lightest black holes.

Events that fall inside the gaps are the interesting ones. GW190521 (85 + 66 M☉) sits in the upper gap. GW230529's primary (2.5–4.5 M☉) sits in the lower one. GW231123 (~225 M☉ total) is the heaviest so far.

**Data.** The GWOSC event API (`gwosc.org/eventapi/`), cumulative GWTC-4.0: 218 candidates with p_astro ≥ 0.5, with median and 90% intervals for masses, spins and distances. The full posterior samples (Zenodo) are optional and large.

**Method + visualizations.**
- The "stellar graveyard": every merger's masses with error bars, beside black holes and neutron stars known from light (I4, C2).
- m1–m2 plane; χ_eff distribution (mostly near zero, slightly positive); chirp mass vs distance, which shows the detection horizon.
- **Selection correction**: an inspiral is heard out to a distance ∝ M_chirp^(5/6), so the surveyed volume grows as M_chirp^(5/2). Weight each event by the inverse and compare the observed and approximate intrinsic mass distributions.

**Check against.** The LVK population papers (GWTC-3: Abbott et al. 2023, *PRX*; GWTC-4.0 population, 2025): the ~10 and ~35 M☉ features.

**Watch out for.** Medians hide broad, correlated posteriors. Source-frame masses assume a cosmology. Our selection weighting is crude (sensitivity changed between runs, and heavy systems merge outside the band), so label it as approximate.

**Builds on.** C1; pairs with I4 (black holes weighed by light) and C2 (neutron stars).

## C6. Gamma-ray bursts ★☆☆

**Science question.** Brief flashes of gamma rays arrive about once a day from random directions. Where do they come from, and why are there two kinds?

**Background.** Gamma-ray bursts were discovered in 1967 by the Vela satellites, which were watching for nuclear tests. BATSE (1991–2000) found them spread evenly over the sky. A Galactic population would crowd toward the Milky Way's plane or centre, so the bursts had to be very distant. BeppoSAX afterglows and redshifts confirmed this in 1997. Durations are **bimodal** (Kouveliotou et al. 1993):
- short bursts (under ~2 s) come from neutron star mergers, as C4's GRB 170817A proved;
- long bursts come from the collapse of massive stars and arrive with supernovae.

Redshifts reach z ≈ 8.2 (GRB 090423).

**Data.** The Fermi GBM Burst Catalog (HEASARC table `fermigbrst`, >4,000 bursts since 2008: T90 duration, fluence, peak flux, spectral fits, positions), via `astroquery.heasarc`. The Swift GRB table for redshifts. The BATSE 4B catalogue (HEASARC) for the historical comparison.

**Method + visualizations.**
- log T90 histogram → two-Gaussian mixture fit and the crossover duration.
- Hardness vs duration: short bursts are spectrally harder.
- Sky map in Galactic coordinates + BATSE's isotropy statistics (dipole ⟨cos θ⟩ toward the Galactic centre, quadrupole ⟨sin² b − ⅓⟩), corrected for GBM's exposure.
- Redshift distributions of long vs short bursts (Swift). GRB 221009A, the brightest burst on record, marked as the outlier it is.

**Check against.** Short/long crossover near ~2 s; isotropy statistics consistent with zero (as BATSE found).

**Watch out for.** GBM positions are uncertain by several degrees, which is fine for isotropy and useless for host galaxies. T90 depends on the instrument's energy band and sensitivity, so the crossover shifts between missions. Redshifts exist only for bursts with bright afterglows.

**Builds on.** C3 (the same isotropy argument used for FRBs), C4.

## C7. Neutrino astronomy ★☆☆ (IceCube extension ★★★)

**Science question.** Neutrinos barely interact with anything, yet on 23 February 1987 about two dozen of them from a supernova in the Large Magellanic Cloud were caught in three detectors. What can 24 particles tell us? And what does IceCube see today?

**Background.** When a massive star's core collapses into a neutron star, ~99% of the ~3×10⁵³ erg released escapes as neutrinos within ~10 s. Kamiokande-II (11 events), IMB (8) and Baksan (5) recorded SN 1987A's burst about three hours before its light was seen. If neutrinos have mass, lower-energy ones travel slightly slower and arrive later. The arrival spread therefore bounds the neutrino mass. This is the same logic as C3's dispersion measure, where lower frequencies arrive later.

**Data.** The published event lists (times, energies, uncertainties): Hirata et al. 1987/1988, Bionta et al. 1987, Alexeyev et al. 1988, compiled in Loredo & Lamb (2002, astro-ph/0107260) and Vissani (2015, arXiv:1409.4710). For the extension: the IceCube 10-year point-source data release (track events from 2008-04-06 to 2018-07-08, with effective areas and resolution functions), at HEASARC and on the IceCube data-release page.

**Method + visualizations.**
- Event times and energies from the three detectors on one plot (with each detector's clock uncertainty).
- Fit a thermal (Fermi–Dirac) spectrum → temperature (~4 MeV). Detector mass and cross-section → total energy emitted. Compare with the binding energy of a newborn neutron star (~3GM²/5R): proof that one formed.
- Neutrino-mass limit from the energy–time spread.
- IceCube extension: an all-sky map of a million track events; an unbinned likelihood point-source test at NGC 1068 and TXS 0506+056 using the release's response functions.

**Check against.** Total energy ~few × 10⁵³ erg; neutrino mass < 5.7 eV (Loredo & Lamb 2002, 95%); IceCube's 10-year result, where NGC 1068 is the hottest northern source at 2.9σ post-trial (Aartsen et al. 2020).

**Watch out for.** Tiny numbers mean Poisson statistics throughout. Kamiokande-II's absolute clock was uncertain by ±1 minute, and some of its low-energy events may be background. For IceCube, atmospheric neutrinos dominate, and trials factors matter. NGC 1068's later 4.2σ evidence (2022) used a reconstruction not included in this release.

**Builds on.** C3 (dispersion), B4 (the LMC's distance sets the energy scale), C2 (the neutron star that formed).

## C8. What is inside a neutron star? ★★☆

**Science question.** A neutron star packs more than the Sun's mass into a ball the size of a city, denser than an atomic nucleus. Nobody knows exactly how matter behaves at that density. Can we combine NICER's radius measurements, the heaviest known pulsars and GW170817 to test the candidate equations of state, and use pulsar glitches to probe the superfluid inside?

**Background.** Given an equation of state (pressure against density), the Tolman–Oppenheimer–Volkoff (TOV) equations, which describe hydrostatic equilibrium in general relativity, give a unique mass–radius curve and a maximum mass. "Stiff" matter makes large stars and a high maximum mass; "soft" matter makes small ones. Three kinds of measurement constrain it:
- **Maximum mass.** Radio timing gives PSR J0740+6620 2.08 ± 0.07 M☉ (Fonseca et al. 2021), which rules out equations of state that cannot hold 2 M☉.
- **Radius.** NICER maps hot spots on spinning pulsars. Light bending makes the pulse shape depend on M/R. PSR J0030+0451 (~1.4 M☉): R ≈ 12.7 (+1.3/−1.0) km in an updated analysis. PSR J0740+6620: 12.4 km (Riley et al. 2021) or 13.7 km (Miller et al. 2021), from two teams' models of the same data. PSR J0437−4715 (1.418 ± 0.044 M☉): 11.4 km (Choudhury et al. 2024), or 11.9–15.5 km (68%) in a 2025 reanalysis (arXiv:2512.08790).
- **Tidal deformability.** GW170817's waveform (C4) limits how much the two stars deformed each other before merging: R = 11.9 ± 1.4 km (LVK 2018).

Glitches are sudden spin-ups when superfluid neutrons inside the crust hand angular momentum to the rest of the star. Vela's regular glitches require at least ~1.4% of the moment of inertia to be in that superfluid (Link, Epstein & Lattimer 1999).

**Data.**
- NICER posterior samples on Zenodo: PSR J0030+0451 (Riley et al. 2019, record 3386449, a 0.9 GB archive) and PSR J0740+6620 (Riley et al. 2021, record 4697625; Miller et al. 2021, record 4670689, with files of 29–390 MB).
- GW170817 posterior samples with tidal deformability (GWOSC, as in C4).
- Equation-of-state tables from CompOSE (`compose.obspm.fr`), or the piecewise polytropes of Read et al. (2009).
- The Jodrell Bank Glitch Catalogue: 728 glitches in 223 pulsars (Basu et al. 2022; `jb.man.ac.uk/pulsar/glitches.html`).

**Method + visualizations.**
- Integrate the TOV equations for a dozen equations of state → mass–radius curves. Start with Oppenheimer & Volkoff's 1939 free neutron gas, whose maximum mass of ~0.7 M☉ is far too low: the historical lesson that nuclear forces matter.
- All constraints on one plot: the 2 M☉ line, NICER's mass–radius contours and GW170817's radius band. Which equations of state survive?
- A simple Bayesian fit of a piecewise polytrope to all three.
- Glitches: the size distribution (two populations), Vela's glitches every ~3 years, and glitch activity against spin-down rate → the superfluid fraction.

**Check against.** Riley et al. 2021 and Miller et al. 2021 radii for J0740+6620; LVK 2018: R = 11.9 ± 1.4 km; Fonseca et al. 2021: 2.08 ± 0.07 M☉; Link et al. 1999: ≥1.4%.

**Watch out for.** NICER radii depend on the assumed hot-spot shapes; the two teams' models give different radii from the same data, so show both. Equation-of-state tables treat the crust differently. Glitch catalogues miss small glitches. Tidal deformability depends on both stars' masses.

**Builds on.** C2 (pulsars), C4 (GW170817), C9 (masses from binary pulsars), I8 (millisecond pulsars), L3 (neutron stars at birth).

## C9. Binary pulsars: weighing neutron stars and the first evidence of gravitational waves ★★☆

**Science question.** In 1974 Hulse and Taylor found a pulsar orbiting another neutron star every 7.75 hours. Over decades its orbit shrank exactly as general relativity predicts for energy carried off by gravitational waves (Nobel Prize 1993). Can we repeat the test for every binary pulsar in the ATNF catalogue, and weigh the neutron stars along the way?

**Background.** A pulsar in a binary is a clock on an orbit. Timing it gives the Keplerian orbit plus **post-Keplerian** effects that, in general relativity, depend only on the two masses: periastron advance ω̇, time dilation γ, Shapiro delay (range r and shape s = sin i), and orbital decay Ṗ_b. Any two give both masses; each extra one tests the theory. The quadrupole formula (Peters & Mathews 1963) predicts `Ṗ_b = −(192π/5) (G/c³)^(5/3) (P_b/2π)^(−5/3) f(e) m_p m_c (m_p + m_c)^(−1/3)`, with `f(e) = (1 + 73e²/24 + 37e⁴/96)(1 − e²)^(−7/2)`. For PSR B1913+16 the observed decay is 0.9983 ± 0.0016 times the prediction (Weisberg & Huang 2016). The double pulsar confirms it to 1.3 × 10⁻⁴ (Kramer et al. 2021). The observed Ṗ_b must first be corrected for the Galaxy's acceleration and for the Shklovskii effect (the pulsar's sideways motion makes its distance grow). Neutron stars in double neutron star binaries cluster near 1.35 M☉, while recycled pulsars with white dwarf companions reach 2 M☉ (Özel & Freire 2016).

**Data.** The ATNF Pulsar Catalogue via `psrqpy` (as in C2): orbital parameters, PBDOT, OMDOT, GAMMA, Shapiro-delay M2 and SINI, proper motions and distances. Published timing solutions for the showcase systems (Weisberg & Huang 2016; Kramer et al. 2021). A Galactic potential from `galpy` (H3).

**Method + visualizations.**
- For each binary with a measured Ṗ_b: GR's prediction from the masses against the observed value after removing the Galactic and Shklovskii terms, as a ratio plot. Hulse–Taylor's famous parabola of the accumulated shift in periastron time, recomputed from the published parameters.
- Mass–mass diagrams for B1913+16 and the double pulsar: each post-Keplerian effect is a band, and all bands cross at one point. That crossing is the test of GR.
- The neutron star mass distribution from every system with two post-Keplerian parameters or a Shapiro delay: double neutron stars against pulsar–white dwarf binaries, and the 2 M☉ pulsars that matter for C8.
- Time to merger for each double neutron star (Peters 1964): which will become LIGO sources like C4's within the age of the universe?

**Check against.** Weisberg & Huang 2016: 0.9983 ± 0.0016; Kramer et al. 2021: agreement to 1.3 × 10⁻⁴ (95% confidence).

**Watch out for.** For nearby or wide binaries the Galactic and Shklovskii corrections can exceed the GR term, and distances from dispersion measure are often 20–50% uncertain. Catalogue values come from different papers and timing models (DD, ELL1, DDGR), so use consistent ones. Shapiro-delay masses need nearly edge-on orbits.

**Builds on.** C2 (pulsars, `psrqpy`), C1 and C4 (gravitational waves), A3 (Keplerian orbits), H3 (`galpy`).

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

## D5. Is dark energy changing? ★★☆ (own BAO measurement ★★★)

**Science question.** In 2024–25 the DESI survey reported a hint that dark energy is weakening over time. Can we reproduce that hint by adding DESI's measurements to D1's supernova likelihood? And how much does it depend on which supernova sample we use?

**Background.** Sound waves in the plasma of the early universe froze at recombination and left a preferred separation between galaxies. This is the **sound horizon** r_d ≈ 147 Mpc, the same physics that makes D3's CMB acoustic peaks. Measuring that ruler across the sky and along the line of sight at several redshifts gives D_M/r_d and D_H/r_d. This is **BAO**: baryon acoustic oscillations. Supernovae (D1) trace the *shape* of the expansion history, BAO its scale in units of r_d, and the CMB fixes r_d. Dark energy's equation of state is often written `w(a) = w0 + wa (1 − a)`; a cosmological constant is w0 = −1, wa = 0. DESI DR2 (2025) found that DESI + CMB + supernovae prefer w0 > −1, wa < 0. The significance depends on the supernova sample: 2.8σ with Pantheon+, 3.8σ with Union3, 4.2σ with DES-Y5.

**Data.**
- The DESI DR2 BAO table (DESI Collaboration 2025, arXiv:2503.14738): D_M/r_d and D_H/r_d with correlations at z = 0.3–2.33. The numbers are small enough to type in.
- D1's full-covariance Pantheon+ likelihood, which already matches the published chains.
- A compressed CMB prior on (θ*, Ω_b h², Ω_c h²), either from D3's own fit or the one given in the DESI paper.
- Optional: DES-Y5 supernovae (public on GitHub) to test the sample dependence.
- For our own BAO measurement: BOSS DR12 CMASS/LOWZ large-scale-structure catalogues with randoms (SDSS SAS), or DESI DR1 LSS catalogues (`data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/`).

**Method + visualizations.**
1. BAO alone in ΛCDM → Ωm and r_d·h. DESI's Ωm (≈0.30) vs D1's flat Ωm (0.332 ± 0.018): the mild disagreement that drives the result.
2. Add the CMB prior, then fit w0–wa (MCMC with `emcee`) for BAO + CMB + Pantheon+. Draw the w0–wa contours for each probe and for the combination. Compute Δχ² against ΛCDM → the significance.
3. Swap in DES-Y5 and watch the significance move.
4. Stretch: the correlation function ξ(s) of BOSS CMASS or DESI DR1 LRGs (Landy–Szalay estimator with randoms and weights, reusing D2's pair-count code) → the BAO bump near 100 h⁻¹ Mpc (Eisenstein et al. 2005) → our own D_V/r_d.
- Also: a Hubble diagram of the BAO distances with the ΛCDM and best w0wa curves, and the implied w(z) band.

**Check against.** DESI DR2 ΛCDM Ωm = 0.2975 ± 0.0086; 2.8σ preference for w0wa with Pantheon+.

**Watch out for.** w0–wa is a convenient parametrization, not a physical model. An apparent "phantom crossing" (w < −1 in the past) may be an artefact of it. The significance depends on supernova calibration (see the debate over DES-Y5's low-redshift sample, e.g. Efstathiou 2024). Quote the result as "Δχ² = X for 2 extra parameters", not as a discovery.

**Builds on.** D1 (likelihood), D2 (clustering code), D3 (CMB prior, sound horizon). Together with B4 and C4, a fourth independent angle on the Hubble tension.

## D6. Fresh supernovae onto our Hubble diagram ★★☆

**Science question.** Pantheon+ was frozen in 2022. Can we fit light curves of Type Ia supernovae that it never saw, and check whether they land on D1's Hubble diagram?

**Background.** Type Ia supernovae are *standardizable*, not standard. Brighter ones fade more slowly (larger "stretch" x1), and redder ones look dimmer (colour c). The SALT light-curve model fits (peak brightness, x1, c) to multi-band photometry. The Tripp relation `μ = m_B − M + α x1 − β c` then gives the distance.

**Data.** The ZTF Type Ia supernova Data Release 2 (Rigault et al. 2025): ~3,600 spectroscopically classified SNe Ia at z < 0.3 from 2018–2020, with public g/r/i light curves, a sample not in Pantheon+. For the live element: recent classifications from the Transient Name Server (TNS; free account for API use), the ZTF Bright Transient Survey explorer, and ZTF photometry from the ALeRCE broker API (public).

**Method + visualizations.**
- SALT3 fits with `sncosmo` (Milky Way dust from the SFD map), with fit quality cuts.
- Standardize with D1's α, β and M → distance moduli → overplot on D1's Hubble diagram and residual plot.
- Hubble-residual scatter before and after standardization: the "standardizable" in action.
- Supernovae classified in the last month, fitted the same way: the Hubble diagram updated with events a few weeks old.

**Check against.** Hubble-residual scatter ~0.15 mag after standardization; ZTF DR2 papers' own results.

**Watch out for.** ZTF's g and r bands constrain colour less well than Pantheon+'s surveys. At z < 0.02, peculiar velocities of ~300 km/s add ~0.1 mag of scatter. Cross-survey calibration offsets mean we should *not* fit H0 from this sample. Watch for misclassified peculiar types (Iax, 91bg-like).

**Builds on.** D1.

*D7–D11 ask whether the universe really expands, how its expansion rate has changed, and what dark energy does to the CMB and to the growth of structure.*

## D7. Is the universe really expanding? ★☆☆

**Science question.** Redshift grows with distance. Is that because space expands, or could light simply lose energy on its long journey ("tired light", Zwicky 1929)? Three tests separate the two. In an expanding universe distant events run in slow motion, the CMB was hotter in the past, and distant galaxies look far dimmer per unit area than their distance alone explains.

**Background.** Expansion stretches every interval by (1 + z), so a supernova at z = 1 takes twice as long to rise and fade. The CMB cools as space expands, so at redshift z it should be `T = 2.725 K × (1 + z)`. A galaxy's surface brightness falls as (1 + z)⁻⁴ (Tolman 1930). Tired light predicts no time dilation, an unchanging CMB temperature, and dimming of only (1 + z)⁻¹.

**Data.**
- The DES 5-year supernova release (GitHub `des-science/DES-SN5YR`, folder `0_DATA`): light curves of the 1,504 Type Ia supernovae at 0.1 < z < 1.2 that White et al. (2024) used to measure time dilation.
- Published CMB temperatures at high redshift: CO and C I absorption in quasar spectra (Noterdaeme et al. 2011 and later), Sunyaev–Zel'dovich spectra of Planck clusters, and water absorption at z = 6.34 (Riechers et al. 2022: 16.4–30.2 K, against ~20 K expected).
- The Tolman test: HST surface brightnesses of elliptical galaxies in three clusters at z ≈ 0.76–0.92 (Lubin & Sandage 2001).

**Method + visualizations.**
- DES light curves on an observed time axis, then with time divided by (1 + z): the spread collapses. Fit width ∝ (1 + z)^b → b.
- T_CMB vs z with the (1 + z) line and tired light's flat line; fit `T = T₀ (1 + z)^(1−β)`.
- Surface brightness vs z against (1 + z)⁻⁴ and (1 + z)⁻¹, allowing for younger galaxies being brighter.
- A summary panel: the three tests side by side, with tired light failing each.

**Check against.** White et al. 2024: b = 1.003 ± 0.005 (stat) ± 0.010 (sys). Noterdaeme et al. 2011: β = −0.007 ± 0.027.

**Watch out for.** Supernova light-curve widths also depend on stretch, which correlates with host type and, through selection, with redshift; compare like with like. Molecular thermometers need corrections for collisions and local radiation. The Tolman test is the weakest of the three because galaxies were brighter in the past, so it constrains expansion and evolution together.

**Builds on.** D1 and D6 (supernovae), D3 (FIRAS's 2.725 K).

## D8. The universe's speedometer ★☆☆

**Science question.** D1 inferred the acceleration from distances, which add up the whole expansion history. Can we measure the expansion rate H(z) itself, epoch by epoch, watch deceleration turn into acceleration, and extrapolate to the far future?

**Background.** **Cosmic chronometers** (Jimenez & Loeb 2002) use `H(z) = −1/(1 + z) · dz/dt`. The most massive passive galaxies formed their stars early and fast. The difference in their stellar ages between two nearby redshifts therefore measures dt directly (from B7's 4000 Å break or full spectral fits), and dz is just the redshift difference. No cosmological model goes in, though stellar-population models do. Two galaxies a fixed comoving distance apart separate at a rate proportional to ȧ = H(z)/(1 + z). ȧ falls while gravity decelerates the expansion and rises once dark energy dominates, so its minimum marks the switch. D1 put the switch at z = 0.59.

**Data.** 32 chronometer measurements of H(z) at 0.07 ≤ z ≤ 1.965, compiled in Moresco et al. (2022, *Living Reviews in Relativity*). 15 of them come with a full systematic covariance in Moresco's `CCcovariance` repository on GitLab. For comparison: D5's DESI D_H/r_d converted to H(z), and D1's fit.

**Method + visualizations.**
- H(z) against ΛCDM (D1/D3 parameters), a matter-only universe and an empty one.
- ȧ(z) and its dip. A smooth fit or a Gaussian process gives the transition redshift and its uncertainty.
- The deceleration parameter `q(z) = (1 + z) H′/H − 1`, and H0 from the chronometers alone.
- The future: a(t) from the Big Bang to +100 Gyr for ΛCDM (exponential growth), for D5's w0–wa fit, and for w < −1. With w < −1 a "Big Rip" tears everything apart in finite time (Caldwell, Kamionkowski & Weinberg 2003).

**Check against.** D1's transition redshift of 0.59; D5's DESI H(z) points.

**Watch out for.** Stellar-population systematics (initial mass function, metallicity, star-formation history) dominate the error budget, so a diagonal-only fit overstates the precision. Progenitor bias: galaxies selected as passive at different redshifts must belong to the same population. Extrapolating w0–wa into the future is speculation; label it as such.

**Builds on.** D1, D5, B7 (the 4000 Å break), J3 (passive galaxies).

## D9. H0 from the flicker of lensed quasars ★★☆

**Science question.** When a foreground galaxy splits a background quasar into several images, the light of each image travels a different path and arrives at a different time. When the quasar flickers, the images flicker in turn. Can we measure those delays ourselves from 15 years of monitoring, and turn them into H0 without any distance ladder?

**Background.** The delay between two images is `Δt = (D_Δt / c) ΔΦ`. ΔΦ, the difference in "Fermat potential", is set by the lens galaxy's mass distribution. The **time-delay distance** `D_Δt = (1 + z_l) D_l D_s / D_ls` is proportional to 1/H0 (Refsdal 1964). Delays run from days to months. Stars in the lens galaxy add slow microlensing wiggles to each image independently. The lens mass profile is the main systematic. Adding a uniform sheet of mass and rescaling leaves the images unchanged but changes the delays: the **mass-sheet degeneracy** (Falco, Gorenstein & Shapiro 1985). H0LiCOW found 73.3 (+1.7/−1.8) km/s/Mpc from six lenses with simple mass profiles (Wong et al. 2020). TDCOSMO 2025 used stellar kinematics to break the degeneracy and found 71.6 (+3.9/−3.3) from eight lenses.

**Data.** COSMOGRAIL (Millon et al. 2020): R-band light curves of 23 lensed quasars from 15 years with the 1.2 m Euler telescope (VizieR `J/A+A/640/A105`). Published lens-model products: H0LiCOW's time-delay distance posteriors (GitHub `shsuyu/H0LiCOW-public`) and the TDCOSMO 2025 results (arXiv:2506.03023).

**Method + visualizations.**
- Light curves of RX J1131−1231 and HE 0435−1223, one colour per image. Shifted by the delays, they line up.
- Delays measured our own way (cross-correlation plus a spline for microlensing, as in I2), then with `PyCS3`, COSMOGRAIL's own tool.
- Our delays × the published Fermat potentials → D_Δt → H0 for each lens. The combination drawn over D1, D3, C4 and I9.
- The mass-sheet degeneracy made visible: lens models with identical images but different delays.

**Check against.** Millon et al. 2020's delays for the same lenses; Wong et al. 2020: 73.3 (+1.7/−1.8); TDCOSMO 2025: 71.6 (+3.9/−3.3).

**Watch out for.** Microlensing variability is the main systematic in the delays: fit it, don't ignore it. Seasonal gaps alias delays near whole years. The lens model, not the delay, dominates the error budget, and published models carry their assumptions with them; say which ones we inherit. Mass along the line of sight also biases D_Δt.

**Builds on.** I2 (cross-correlation lags), J5 (lens models), D1, D3, C4 and I9 (H0).

## D10. Dark energy in the CMB: the integrated Sachs–Wolfe effect ★★☆

**Science question.** Can we see dark energy acting on the CMB? Photons crossing a huge supercluster or void on their way to us should pick up a tiny temperature change, but only if the gravitational wells and hills are decaying. In a flat universe, only dark energy makes them decay.

**Background.** A photon gains energy falling into a gravitational well and loses it climbing out. In a flat, matter-only universe the wells stay constant and the two effects cancel. If dark energy makes a well shallower while the photon is inside, the photon leaves with a net gain. Superclusters then appear as slightly hot spots and supervoids as cold spots: the **integrated Sachs–Wolfe (ISW) effect** (Sachs & Wolfe 1967; Crittenden & Turok 1996). The signal is a few µK against ±100 µK primary fluctuations, so it appears only after stacking many structures or cross-correlating with galaxy maps. Granett et al. (2008) stacked 50 supervoids and 50 superclusters from SDSS luminous red galaxies and found 9.6 ± 2.2 µK, larger than ΛCDM predicts. With BOSS voids and a matched filter, Nadathur & Crittenden (2016) measured 1.64 ± 0.53 times the ΛCDM amplitude (3.1σ).

**Data.** The Planck SMICA map and masks (from D3). Void and supercluster catalogues: BOSS DR12 voids (Mao et al. 2017, VizieR `J/ApJ/835/161`) and SDSS DR7 voids and superclusters (Nadathur & Hotchkiss 2014, `J/MNRAS/440/1248`). The page that hosted Granett et al.'s list is offline (try the Internet Archive). For cross-correlations: D2's galaxies or an all-sky galaxy map.

**Method + visualizations.**
- Stacked CMB patches centred on voids and on superclusters (compensated top-hat filter), with errors from random positions.
- The ΛCDM prediction from linear theory with D3's parameters, and the zero prediction of a matter-only universe.
- A galaxy–CMB cross-correlation, with errors from rotated maps.
- The CMB Cold Spot and the Eridanus supervoid: how much of the Cold Spot can a void explain?

**Check against.** Granett et al. 2008: 9.6 ± 2.2 µK; Nadathur & Crittenden 2016: A_ISW = 1.64 ± 0.53; Planck 2015 results XXI: ISW detected at ~4σ with all tracers combined.

**Watch out for.** Choosing the filter size or structure count after looking at the result inflates the significance, so fix them first. The primary CMB is the noise here, and better data cannot reduce it; only more structures can. Galactic foregrounds and the mask matter. Voids found in different tracers have different expected signals, so predict the signal per catalogue.

**Builds on.** D3 (maps, masks), D2 (galaxy catalogues), D1 and D5 (dark energy).

## D11. Counting galaxy clusters: how clumpy is the universe? ★★☆

**Science question.** Galaxy clusters are the largest objects gravity has had time to build. How many exist at each mass and epoch depends on how much matter there is and how fast structure grows, and dark energy slows that growth. Can we count eROSITA's X-ray clusters and measure how clumpy the universe is?

**Background.** The **halo mass function** (Press & Schechter 1974; Tinker et al. 2008) gives the number of dark matter halos per volume above mass M at redshift z. Its high-mass end is exponentially sensitive to Ωm and to σ8, the typical density fluctuation on 8 h⁻¹ Mpc scales. Surveys usually quote `S8 = σ8 (Ωm / 0.3)^0.5`. Several low-redshift measurements came out 2–3σ below the CMB's S8: the **S8 tension**. eROSITA's first all-sky survey (eRASS1) counted 5,259 clusters over 12,791 deg² and found S8 = 0.86 ± 0.01, in line with the CMB (Ghirardini et al. 2024). As of 2026, KiDS-Legacy cosmic shear agrees with the CMB, while DES Y6 sits 2.4–2.7σ low.

**Data.** The eRASS1 cluster catalogues (Bulbul et al. 2024): the primary catalogue of 12,247 optically confirmed clusters and groups (48 MB) and the cosmology subsample, on the eROSITA-DE DR1 site (`erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/BulbulE_DR1/`). Planck's PSZ2 catalogue (1,653 SZ detections) for comparison. Halo mass functions from `colossus`.

**Method + visualizations.**
- Sky map of the eRASS1 clusters (the western Galactic hemisphere) and their redshift distribution.
- Counts in bins of X-ray luminosity and redshift against predictions on an (Ωm, σ8) grid, using the published mass–observable relation and selection function.
- Ωm–σ8 contours, and S8 beside Planck (D3), KiDS-Legacy and DES Y6.
- Dark energy's brake: counts predicted with and without Λ, normalized to the same CMB fluctuations at z ≈ 1100.

**Check against.** Ghirardini et al. 2024: Ωm = 0.29 (+0.01/−0.02), σ8 = 0.88 ± 0.02, S8 = 0.86 ± 0.01.

**Watch out for.** The mass calibration (from weak lensing of the clusters) and the selection function decide everything: a 10% mass bias moves σ8 substantially. Using the published scaling relation tests consistency, not independence. Planck's SZ counts needed a hydrostatic mass bias (1 − b) of about 0.6 to agree with the CMB; K2 measures that bias directly.

**Builds on.** D3 (σ8 from the CMB), D2 (clustering), K2 (cluster masses), D5 (dark energy).

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

## E4. Meteor showers and their parent comets ★★☆

**Science question.** A meteor shower happens when Earth crosses a stream of debris. Can we find the showers blind in a database of camera-measured meteor orbits, and trace each one to the comet or asteroid that shed it?

**Background.** Several cameras filming the same meteor let us triangulate its path through the atmosphere and reconstruct the particle's orbit around the Sun before impact. Meteors of one shower travel on parallel paths, so by perspective they appear to radiate from one point on the sky, the **radiant**. Parent bodies:
- Perseids ← comet 109P/Swift–Tuttle
- Leonids ← 55P/Tempel–Tuttle
- Orionids and η-Aquariids ← 1P/Halley
- Quadrantids ← 2003 EH1
- Geminids ← the *asteroid* (3200) Phaethon

The **D-criterion** (Southworth & Hawkins 1963) measures how similar two orbits are.

**Data.** Global Meteor Network (GMN) trajectory summary files: a volunteer camera network running since December 2018, with millions of meteor orbits under CC BY 4.0 (`globalmeteornetwork.org/data/`; column definitions in the GMN orbit-data PDF). Parent orbits from the JPL Small-Body Database API (as in E1/E3). Shower names from the IAU Meteor Data Center list, used only *after* blind detection.

**Method + visualizations.**
- **Radiant map** in Sun-centred ecliptic coordinates: showers appear as knots on top of the sporadic background. That background has its own structure (the helion, antihelion, apex and toroidal sources), which is worth explaining.
- Blind clustering (radiant + entry speed + solar longitude, or D-criterion linking) → shower list → match to IAU names.
- D-criterion between each shower's mean orbit and every comet and near-Earth asteroid → parents found, not assumed.
- Entry-speed histogram: Leonids ~71 km/s (head-on, retrograde) vs Geminids ~35 km/s. Activity profile vs date; radiant drift day by day.
- A shower's orbits and its parent's in E2's 3D viewer.

**Check against.** The known shower–parent pairs above; published GMN shower parameters.

**Watch out for.** Most cameras are in the northern hemisphere, so southern showers are undercounted. Deceleration corrections matter for the orbits of short, faint meteors. The D-criterion threshold is a choice: show how the results change with it. Sporadic meteors contaminate weak showers.

**Builds on.** E1 (orbital elements), E2 (3D orbits), E3 (fireballs), G2 (the meteorites that survive).

## E5. Dating surfaces by counting craters ★★☆

**Science question.** Older surfaces have more craters. Can we read the ages of the Moon's dark "seas" and bright highlands, and of Mars's two hemispheres, from crater counts alone?

**Background.** The cumulative number of craters larger than a given diameter falls as a power law with size. The density of craters ≥1 km per km² grows with exposure time. A **chronology function** turns density into age. It is calibrated against radiometric ages of Apollo and Luna samples (Neukum et al. 2001): `N(≥1 km) = 5.44×10⁻¹⁴ (e^(6.93 T) − 1) + 8.38×10⁻⁴ T` per km², with T in Gyr. The exponential term is the early heavy bombardment; the linear term is today's steady rain. The oldest highlands are saturated: new craters erase old ones as fast as they form.

**Data.** Robbins (2019) Lunar Crater Database: ~1.3 million craters, complete above ~1–2 km. Robbins & Hynek (2012) Mars crater database: 384,343 craters ≥1 km. Both at USGS Astrogeology. Unit boundaries: the USGS Unified Geologic Map of the Moon (Fortezzo et al. 2020) and the Mars global geologic map (Tanaka et al. 2014).

**Method + visualizations.**
- Global crater-density maps. The maria stand out as smooth; Mars's north–south dichotomy is obvious.
- Size–frequency distributions and "R-plots" (a standard way to show crater densities relative to a reference slope) for maria (e.g. Imbrium, Tranquillitatis) vs highlands. Saturation shows at small sizes.
- Model ages from the chronology function → compare with Apollo landing sites' sample ages. Mars: ages of the northern plains, the southern highlands and the Tharsis volcanoes.
- Cross-check with E3: does today's measured impact rate match the chronology's linear term? The comparison needs a crater-to-impactor size scaling, which is itself worth explaining.

**Check against.** Apollo 11 mare basalts ~3.6–3.9 Gyr and Apollo 12 ~3.1–3.3 Gyr; crater-count mare ages (Hiesinger et al. 2011) mostly 3.3–3.8 Gyr, with a few maria much younger.

**Watch out for.** Circularity: the chronology is calibrated on these same Apollo sites, so agreement there is a consistency check, not an independent test. Secondary craters (thrown out by larger impacts) inflate small-crater counts. Measure areas on the sphere, not in map projection. Completeness depends on image resolution and lighting.

**Builds on.** E1 (the asteroid belt supplies the impactors), E3 (today's impact rate).

*E6–E10 are chosen for their visuals: Saturn's rings, Jupiter's storms, eclipse shadows, asteroid shapes and comet tails.*

## E6. Saturn's rings in radio light ★★☆

**Science question.** From Earth Saturn's rings look smooth, but they are countless orbiting ice particles arranged in thousands of ringlets, gaps and waves. Cassini shone radio waves through the rings to Earth. Can we read their structure, find the moons that sculpt them, and weigh the rings from the waves those moons raise?

**Background.** During a radio occultation, Cassini transmitted at 0.94, 3.6 and 13 cm wavelengths while passing behind the rings as seen from Earth. The dimming of the signal gives the optical depth along a track across the rings, reconstructed to ~1 km radial resolution. Moons sculpt the rings through resonances, just as Jupiter carves E1's Kirkwood gaps. Ring particles near 117,560 km orbit twice for every orbit of Mimas, and that 2:1 resonance holds the sharp outer edge of the B ring, the inner edge of the Cassini Division. Small moons clear the Encke gap (Pan) and the Keeler gap (Daphnis). At weaker resonances, moons launch spiral **density waves**. A wave's wavelength shrinks with distance from the resonance at a rate set by the ring's surface density, so the waves weigh the rings: tens of grams per cm² in the A ring. All the rings together weigh ~1.5 × 10¹⁹ kg, about 0.4 Mimas masses (Iess et al. 2019).

**Data.** The PDS Ring-Moon Systems Node's Cassini RSS ring occultation archive (`pds-rings.seti.org`, volumes `CORSS_8xxx`, data set CO-SR-RSS-4/5-OCC-V2.0): radial optical-depth profiles from occultations between 2005 and 2010. Cassini UVIS stellar occultations, also at the Rings Node, for finer resolution. Saturn's gravity field (J2, J4) and the moons' orbits from JPL for resonance locations.

**Method + visualizations.**
- The keeper: the whole ring system as one optical-depth profile from the C ring to the F ring, every named gap and ringlet labelled, and the same profile wrapped into a face-on image of the rings.
- Compute moon resonances (Mimas, Janus/Epimetheus, Prometheus, Pandora…) and mark where they fall: which edges and waves line up?
- Zoom on one A-ring density wave and fit its wavelength against radius → surface density.
- Compare the three radio wavelengths: differences in optical depth reveal how many centimetre-sized particles there are.
- Ingress against egress of one occultation, and different occultations against each other: some ringlets are eccentric.

**Check against.** The Mimas 2:1 resonance at the B ring's outer edge; published A-ring surface densities from density waves (e.g. Tiscareno et al. 2007); the ring mass of Iess et al. 2019.

**Watch out for.** Optical depth depends on the viewing geometry: convert to normal optical depth using the ring opening angle. Resolution varies between occultations. The densest parts of the B ring block the signal completely. Radius scales differ by a few km between occultations; align them on sharp edges.

**Builds on.** E1 (resonances), E2 (orbits), C1 (signal processing).

## E7. Jupiter's weather from a decade of Hubble maps ★★☆

**Science question.** Jupiter's clouds race around the planet in alternating bands, and the Great Red Spot, a storm wider than Earth, has been shrinking for over a century. Hubble has mapped the whole planet every year since 2014. Can we measure the winds ourselves and watch the Red Spot change?

**Background.** Jupiter's atmosphere is organized into east–west jets that alternate with latitude; the fastest, near 24°N, blows at ~150 m/s. Jupiter rotates in under 10 hours, so Hubble maps the whole planet on two consecutive rotations. Cross-correlating the two maps shows how far each cloud moved in ~10 hours, giving the wind at every latitude. The Great Red Spot was about 40,000 km long in the 1880s, 23,300 km when Voyager passed in 1979, and ~16,500 km in 2015, while also becoming rounder (Simon et al. 2015, 2018). OPAL (Outer Planet Atmospheres Legacy) has mapped Jupiter every year since 2014, and Saturn, Uranus and Neptune too.

**Data.** OPAL global maps of Jupiter for every year (MAST High Level Science Product `opal`): cylindrical maps from two rotations (A and B) in several filters, from the ultraviolet (F275W) to the near-infrared. JunoCam images (public) for optional close-ups.

**Method + visualizations.**
- The keeper: a decade of Jupiter maps as an animation, and the same maps wrapped onto spinning globes.
- The wind profile: cross-correlate each year's two maps in latitude strips → wind speed against latitude; compare the years.
- The Great Red Spot: fit an ellipse each year → length, width and drift in longitude; extend the series back in time with published historical sizes.
- Colour changes in the belts and zones, and new storms appearing.
- The same tools on OPAL's Saturn maps.

**Check against.** OPAL wind profiles (Simon et al. 2015) and Great Red Spot sizes (Simon et al. 2018).

**Watch out for.** The map projection stretches high latitudes; mask the distorted limb. Cross-correlation fails where features change shape or have no contrast. Longitudes are given in System III (tied to Jupiter's magnetic field), not to the clouds.

**Builds on.** E2 (planets), F2 (rotation from cross-correlation), G1 (colour images).

## E8. Total eclipses: draw the Moon's shadow ★★☆

**Science question.** On 12 August 2026 the Moon's shadow crossed Greenland, Iceland and Spain. On 2 August 2027 it will cross southern Spain, North Africa and Egypt, giving up to 6 minutes 23 seconds of totality near Luxor. Can we compute these paths ourselves from the planetary ephemeris, animate the shadow sweeping across Earth, and predict the Baily's beads from the mountains on the Moon's edge?

**Background.** A total eclipse happens where the Moon's umbra, its cone of full shadow, touches Earth. Its track follows from the positions of the Sun and Moon and from Earth's rotation. Bessel's method boils each eclipse down to a few polynomials (the "Besselian elements"). The Moon's distance varies, so sometimes the umbra falls short and the eclipse is annular. Eclipses repeat in Saros cycles of 18 years, 11 days and 8 hours: the 2027 eclipse belongs to Saros 136, which also gave the long eclipses of 1991 and 2009. The Moon's edge is not smooth. Its mountains and valleys break the last sliver of Sun into Baily's beads and shift contact times by seconds. NASA's maps do not include that, but the laser altimeter on the Lunar Reconnaissance Orbiter (LRO) has mapped the edge.

**Data.** JPL's DE440 ephemeris (via `skyfield`); Earth's shape (WGS84) and rotation (IERS); LRO LOLA topography (PDS Geosciences Node) for the lunar limb profile. NASA's eclipse predictions (Espenak, `eclipse.gsfc.nasa.gov`) and the Five Millennium Canon for checking.

**Method + visualizations.**
- Compute the umbra's track for the 2026 and 2027 eclipses: the central line, path limits and duration.
- The keeper: an animation of the shadow sweeping the globe, with the partial phases fading around it.
- Duration along each path; why the 2026 shadow, falling at a grazing angle, is stretched near Iceland.
- The lunar limb seen from one site, built from LOLA heights for that moment's libration → Baily's beads and corrected contact times.
- The whole Saros 136 family on one world map.

**Check against.** NASA's greatest durations: 2 min 18.2 s (2026) and 6 min 23.2 s (2027); Espenak & Meeus's Five Millennium Canon.

**Watch out for.** Time scales: mixing TT and UT (ΔT ≈ 69 s) shifts the path by tens of kilometres. The Moon's centre of figure is offset from its centre of mass by about half a kilometre. Refraction changes the observed contact times near the horizon.

**Builds on.** E2 (ephemerides), E5 (lunar topography), F2 (the Sun).

## E9. Asteroids in 3D ★★☆

**Science question.** Asteroids are not round. Their brightness changes as they spin, and from those changes we can recover their shapes. Spacecraft have seen some up close. Can we build a 3D gallery, find the spin rate at which rubble-pile asteroids would fly apart, and measure how much NASA's DART impact changed a small moon's orbit?

**Background.** An elongated asteroid brightens and fades twice per rotation. Light curves taken from many viewing angles can be inverted into a convex 3D shape and spin axis (Kaasalainen & Torppa 2001). DAMIT, a database run by Charles University in Prague, holds 16,098 such models of 10,758 asteroids. Plotting rotation period against size shows a sharp limit: almost no asteroid larger than a few hundred metres spins faster than once every ~2.2 hours, the rate at which a loose rubble pile would fly apart (Pravec & Harris 2000). Spacecraft confirmed that Itokawa, Ryugu and Bennu are rubble piles. In September 2022 NASA's DART spacecraft hit Dimorphos, the ~150 m moon of Didymos. The impact shortened its 11.9-hour orbit by 33.0 ± 1.0 minutes, far more than the ~7 minutes from the spacecraft's momentum alone: the ejected debris pushed too (Thomas et al. 2023).

**Data.** DAMIT (`astro.troja.mff.cuni.cz/projects/damit/`): downloadable shape models and spin states. The Asteroid Lightcurve Database (LCDB, Warner et al. 2009; `minplanobs.org`) for periods and sizes. Spacecraft shape models from the PDS Small Bodies Node (Eros, Itokawa, Bennu, Ryugu, Didymos and Dimorphos).

**Method + visualizations.**
- The keeper: a 3D gallery of asteroid shapes at true relative scale, from DAMIT and spacecraft, as interactive meshes.
- Make a light curve from a shape: render a model under sunlight as it rotates, and compare with a measured light curve.
- The spin barrier: period against diameter for tens of thousands of asteroids from the LCDB, with the 2.2-hour line. The fast spinners below ~200 m must be held together by strength.
- DART: Dimorphos's orbit before and after the impact, and the momentum enhancement factor β from the period change (Cheng et al. 2023).

**Check against.** Thomas et al. 2023: −33.0 ± 1.0 min (3σ); the ~2.2-hour spin barrier (Pravec & Harris 2000).

**Watch out for.** Light-curve inversion gives only convex shapes (concavities are invisible), and the pole has a mirror ambiguity. LCDB periods carry quality codes; use only reliable ones. Small, fast, bright asteroids are easier to measure, which biases the sample.

**Builds on.** E1 and E3 (asteroids), E5 (impacts), G1 (rendering).

## E10. Comet tails sculpted by sunlight ★★☆

**Science question.** A comet's dust tail curves away from the Sun because sunlight pushes on every grain. In October 2024 comet C/2023 A3 (Tsuchinshan–ATLAS) swept through the field of view of SOHO's coronagraph as the second-brightest comet SOHO had ever recorded. Can we reproduce the shape of its tail from physics, and read off when the dust left the nucleus and how big the grains are?

**Background.** Radiation pressure pushes dust outward with a force that, like gravity, falls off as 1/r². Each grain therefore moves as if the Sun were lighter by a factor (1 − β), where β, the ratio of radiation force to gravity, is larger for smaller grains: β ≈ 0.57/(ρ s) for density ρ in g/cm³ and radius s in µm. Grains released at the same moment but with different β lie along a **synchrone**; grains with the same β released at different times lie along a **syndyne** (Finson & Probstein 1968). Laid over an image, the two families read the tail like a clock and a sieve at once. C/2023 A3 passed 3° from the Sun on 9–10 October 2024 at a phase angle of 173°. Forward scattering by its dust made it shine at about magnitude −4 in SOHO's LASCO C3 coronagraph; only C/2006 P1 (McNaught, −5.5) had been brighter in SOHO's images.

**Data.** SOHO LASCO C3 images from 7–13 October 2024 (and of other comets, e.g. C/2011 W3 Lovejoy), via `sunpy` (the Virtual Solar Observatory) or the Helioviewer API. The comet's orbit from JPL Horizons (E2). Optionally STEREO's heliospheric imagers.

**Method + visualizations.**
- Register the LASCO C3 frames, subtract a model of the corona, and make a movie of the comet's passage.
- Compute syndynes and synchrones from the comet's orbit as seen from SOHO, and project them onto the images.
- The keeper: the tail with the best-fitting syndyne–synchrone grid overlaid → when the bright dust left the nucleus, and the grain sizes.
- The anti-tail: why dust appeared to point toward the Sun when Earth crossed the comet's orbital plane in mid-October.
- Brightness against phase angle: the forward-scattering surge.

**Check against.** SOHO's reported peak of about magnitude −4; Finson & Probstein's original analysis of comet Arend–Roland.

**Watch out for.** LASCO images are dominated by the corona and stray light, so subtract a background model. Compute the geometry for SOHO's position, not Earth's. Gas (ion) tails follow the solar wind, not syndynes: separate them by colour or shape. Near the Sun brightness is boosted by forward scattering, so don't read it as dust production.

**Builds on.** E2 (orbits), E4 (comets and meteor streams), F2 (solar spacecraft data), G1 (image processing).

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

## F3. The Gannon storm: the May 2024 superstorm from Sun to aurora ★★☆

**Science question.** On 10–11 May 2024 auroras were seen as far south as Mexico and southern Europe during the strongest geomagnetic storm in more than 20 years. Can we follow the event from the flares and coronal mass ejections on the Sun, through the solar wind at L1, to the shaking of Earth's magnetic field measured on the ground?

**Background.** In early May 2024, active region AR 13664 produced a series of X-class flares and fast coronal mass ejections (CMEs). Several CMEs merged on the way and reached Earth on 10 May. At L1 the solar wind jumped to 750–950 km/s and its magnetic field reached 73 nT, with a southward component of −50 nT. A southward field is the orientation that lets solar-wind energy into Earth's magnetosphere. On the ground the Kp index reached 9, the top of its scale. Dst, which tracks the ring of current around Earth, fell to −412 nT: the strongest storm since the 2003 Halloween storms, and the sixth-largest since 1957. The 1989 storm that blacked out Québec reached −589 nT; the 1859 Carrington event is estimated at between about −850 and −1,750 nT.

**Data.** SDO AIA images and GOES X-ray flux (F2's pipeline); SOHO LASCO coronagraph images and the CDAW CME catalogue (speeds); 1-minute solar wind data at L1 from OMNIWeb (NASA SPDF); Kp from GFZ Potsdam and provisional Dst from the Kyoto WDC; 1-minute magnetometer data from INTERMAGNET observatories at many latitudes (e.g. via the BGS GIN web service, no registration).

**Method + visualizations.**
- The keeper: one multi-panel animation on a shared clock: the Sun in AIA 131 and 171 Å, the CMEs in LASCO, the solar wind's speed, density and Bz at L1, and ground magnetometers from the Arctic to the equator.
- CME transit times: departure speed (CDAW) against arrival time at L1 → how much the solar wind slowed them.
- A chain of magnetometers: how far toward the equator the auroral currents reached, from where the largest field swings were.
- Dst against time, with the 1989 and 2003 storms and the Carrington estimates on one scale.
- Predict Dst from the solar wind with Burton et al.'s (1975) empirical model, and compare.

**Check against.** Kp 9; provisional Dst −412 nT; solar wind 750–950 km/s and Bz −50 nT at L1.

**Watch out for.** Provisional Dst is revised later; quote the version used. OMNI data are time-shifted to Earth's bow shock, not raw L1 times. Remove each magnetometer's quiet-day baseline. Merged CMEs make one-to-one associations between eruptions and arrivals ambiguous.

**Builds on.** F1 (the solar cycle), F2 (flares, SDO data), C1 (time series).

## F4. The Sun's spectrum as one image ★☆☆

**Science question.** Spread sunlight through a fine enough spectrograph and the rainbow is crossed by thousands of dark lines, each the fingerprint of an element in the Sun's atmosphere. Can we turn one of the most detailed solar spectra ever measured into a single poster, and identify the elements behind the lines?

**Background.** Fraunhofer mapped hundreds of dark lines in sunlight in 1814; in 1859–60 Kirchhoff and Bunsen matched them to elements heated in the laboratory. In 1868 a yellow line seen during an eclipse matched no known element: helium was found in the Sun before it was found on Earth. Each line's depth and width depend on how much of the element there is and on the temperature and pressure where the line forms, which is how solar abundances are measured (Asplund et al. 2021). Fourier-transform spectrometers at Kitt Peak resolved the solar spectrum in extraordinary detail. The famous NOAO poster of the visible spectrum, wrapped into rows, was made from such scans.

**Data.** Kurucz's 2005 solar flux atlas, built from spliced Kitt Peak FTS scans (`kurucz.harvard.edu/sun/fluxatlas2005/`: 13 files, ~40 MB in all). Line lists from the NIST Atomic Spectra Database or VALD. Optionally the IAG solar atlas (Reiners et al. 2016) for comparison.

**Method + visualizations.**
- The keeper: 380–700 nm wrapped into ~50 rows, each pixel coloured by its wavelength and darkened by the measured flux: the classic poster, rebuilt from the data.
- Identify the strongest lines (hydrogen Balmer, Na D, Mg b, Ca II H & K, iron) from NIST/VALD and colour-code the poster by element.
- Which dark bands come from Earth's own atmosphere (O₂, H₂O) rather than the Sun?
- A curve of growth from Fe I line strengths → the classic estimate of the Sun's iron abundance.
- Extension: the Sun's gravitational redshift (~0.6 km/s) hides under convective blueshifts of similar size; González Hernández et al. (2020) measured it with a laser frequency comb.

**Check against.** Line identifications in Moore, Minnaert & Houtgast (1966); the solar iron abundance log ε(Fe) = 7.46 (Asplund et al. 2021).

**Watch out for.** Read the atlas documentation on how the scans were spliced and normalized. Lines overlap heavily, especially in the blue. A curve-of-growth abundance is a teaching estimate (±0.2 dex), not a modern 3D non-LTE result.

**Builds on.** B7 (stellar spectra), F2 (the Sun), K1 (abundances).

---

# G. Imaging & fun

## G1. Color composites of real galaxies ★★☆

**Data.** Survey cutout services: DESI Legacy Survey viewer API, SDSS cutouts, `astroquery.skyview`; for the ambitious, raw Hubble/JWST public FITS frames from MAST.

**Method.** Fetch calibrated FITS images of M51, M31, M104 etc. in 3+ filters → align (`reproject`) → asinh-stretch RGB composite (`astropy.visualization.make_lupton_rgb` — the published algorithm used for SDSS's own images). Honest imaging: we document filter→channel mapping and stretch, i.e. what "true color" does and doesn't mean.

## G2. Meteorite landings map ★☆☆

**Data.** The Meteoritical Society's database of ~45,000 classified meteorites (via NASA's open data portal, CSV with coordinates, mass, class, year).

**Method + visualizations.** World map of falls vs. finds — instantly shows collection bias (deserts and Antarctica light up, oceans and jungles are blank: nobody finds black rocks in forests); mass distribution; discoveries over time. A relaxed one-evening project with a genuinely interesting selection-effects story.

---

# H. The Milky Way as a galaxy

Section B maps the stars around us. This section treats the Milky Way as a whole galaxy: what it weighs, how it rotates, what sits at its centre, and what it has eaten.

## H1. Weigh the black hole at the Galactic Centre ★★☆

**Science question.** A star called S2 circles the centre of the Milky Way every 16 years. From its orbit alone, how massive is the object it orbits, how far away is it, and must it be a black hole?

**Background.** Kepler's third law weighs whatever sits at the focus of an orbit: `M = 4π² a³ / (G P²)`. Astrometry measures the orbit's size in *angle*. Radial velocities measure it in *km/s*. Combining the two gives the distance to the Galactic Centre, R0, geometrically, with no ladder. S2's orbit has P ≈ 16.0 yr and e ≈ 0.88. At pericentre (May 2018) it passed ~120 AU from the centre at ~7,700 km/s, 2.5% of the speed of light. Four million solar masses inside a region that small can only be a black hole. Nobel Prize 2020 (Genzel and Ghez, with Penrose).

**Data.**
- The UCLA Galactic Center Group's S0-2 (= S2) astrometry and radial velocities for 1995–2018, published with Do et al. (2019, *Science*): `astro.ucla.edu/~ghezgroup/gc/gr/data.zip`. That page has been offline since about 2023; the Internet Archive's 1 September 2019 copy of the zip works (see `17-galactic-centre/fetch_data.py`).
- VLT + Keck S2 positions: Gillessen et al. 2009, VizieR `J/ApJ/707/L114`.
- Keck S0-2 radial velocities: Chu et al. 2018, VizieR `J/ApJ/854/12`.

**Method + visualizations.**
- S2's positions on the sky, coloured by epoch. Sgr A* sits at the ellipse's focus, not its centre.
- Fit a Keplerian orbit + M + R0 + reference-frame drift (13 parameters, as Do et al. did) with MCMC → M and R0, and the banana-shaped M–R0 correlation.
- Radial-velocity curve with the 2018 pericentre swing. Extension: test the ~200 km/s relativistic redshift at pericentre (fit a factor f: f = 0 is Newton, f = 1 is Einstein).
- A scale diagram: S2's pericentre vs the Solar System vs the event horizon. Density bound → why it must be a black hole.
- Prediction for I5: the shadow size of Sgr A* implied by our M and R0.

**Check against.** M = 4.28 × 10⁶ M☉ and R0 = 8.32 kpc (Gillessen et al. 2017); GRAVITY (2022): 4.297 × 10⁶ M☉ and 8.277 kpc.

**Watch out for.** Reference-frame alignment (via radio masers) is the main systematic. The Keck and VLT frames differ, so fit an offset rather than mixing them naively. Some epochs are confused by nearby stars; use the data flags. The 12′-per-orbit relativistic precession, which GRAVITY detected in 2020, is below what these data can see.

**Builds on.** E2 (Kepler's third law), D4 (weighing with orbits). Leads to I5 (the shadow) and I1 (where the Milky Way sits on M–σ).

## H2. Our own rotation curve from 21-cm hydrogen ★★☆

**Science question.** D4 measured the rotation curves of 175 other galaxies. Can we measure our own, from inside the disk, using the radio glow of hydrogen gas?

**Background.** Neutral hydrogen emits at 1420.4 MHz (21 cm) when its electron's spin flips, and the line's Doppler shift gives the gas velocity along the line of sight. Looking inward at Galactic longitude l (0° < l < 90°), the fastest gas on each sightline sits at the **tangent point**. There the orbit runs parallel to our view, at radius R = R0 sin l, and the circular speed is `V(R) = v_terminal + Θ0 sin l`. This "terminal velocity" method maps the inner Galaxy's rotation. Spiral arms show up as ridges in the longitude–velocity diagram.

**Data.** HI4PI (HI4PI Collaboration 2016, A&A 594, A116): the full-sky 21-cm survey, 16.2′ resolution, ~1.3 km/s channels, on CDS (`J/A+A/594/A116`). The full cubes are tens of GB. We need only the Galactic plane (|b| ≤ ~2°), so fetch plane tiles or per-position spectra through the VO rather than the whole sky. Comparison curves: Clemens (1985, CO), McClure-Griffiths & Dickey (2016, HI terminal velocities), Eilers et al. (2019, Gaia + APOGEE, outer Galaxy).

**Method + visualizations.**
- The **longitude–velocity diagram** at b = 0°, the keeper, with the spiral arms labelled (Perseus, Outer, Scutum–Centaurus, the Local arm and the "3-kpc arms").
- Terminal velocity vs longitude in the first quadrant (l = 15°–90°) and the fourth (l = 270°–345°) → V(R) inside the solar circle. The known north–south asymmetry.
- The Milky Way's curve plotted among D4's SPARC galaxies; M(<R0) ≈ V²R/G ≈ 10¹¹ M☉.
- Extension: a face-on map of the *outer* Galaxy's hydrogen (outside the solar circle each velocity maps to a single distance). The spiral arms, and the warp: the gas disk bends up on one side and down on the other (Levine, Blitz & Heiles 2006).

**Check against.** Θ0 ≈ 236 ± 7 km/s at R0 ≈ 8.15 kpc (Reid et al. 2019); the Clemens (1985) and McClure-Griffiths & Dickey (2016) curves.

**Watch out for.** Exclude |l| < 15°: the central bar drives non-circular motions there. Streaming motions near spiral arms bias terminal velocities by ~10 km/s. Inside the solar circle each velocity maps to two distances, so no face-on map there. R0 and Θ0 are assumed inputs, so propagate them.

**Builds on.** D4, B2. The same physics as J1, from inside.

## H3. Galactic archaeology with Gaia ★★★

**Science question.** The Milky Way grew by swallowing smaller galaxies and star clusters. Can we find their remains, as torn-up stellar streams and the debris of a giant ancient merger, in Gaia's six-dimensional data?

**Background.** When the Milky Way's tides pull apart a satellite, its stars spread along its orbit as a **stellar stream**. GD-1 is a thin stream from a long-dissolved cluster, stretching over 80° of sky, and it has gaps that may have been punched by dark matter clumps (Price-Whelan & Bonaca 2018). **Gaia–Sausage–Enceladus** (Belokurov et al. 2018; Helmi et al. 2018) was a galaxy that merged ~10 Gyr ago. Its stars still move on very elongated orbits, which form a "sausage" in velocity space. The **phase spiral** (Antoja et al. 2018) is a snail-shell pattern in the vertical positions and velocities of stars near the Sun: the disk is still ringing from a recent disturbance.

**Data.** Gaia DR3 via ADQL:
- the GD-1 field (a ~60° × 10° strip);
- the 33.8 million stars with Gaia radial velocities (for full 3D motions);
- Bailer-Jones et al. (2021) distances.

Optional: APOGEE DR17 chemistry ([Fe/H], [α/Fe]); accreted stars are chemically distinct.

**Method + visualizations.**
- GD-1: rotate into the stream's own frame (Koposov et al. 2010, built into `gala`), select by proper motion along the stream and with an old, metal-poor isochrone filter. The stream appears out of nothing: a before/after pair. Density along the stream → gaps and the "spur".
- Gaia–Sausage–Enceladus: Galactocentric velocities (astropy) → v_R vs v_φ for halo stars (the sausage). Energy vs angular momentum in a Milky Way potential (`galpy`/`gala`) → clumps of accreted debris (GSE, Sequoia, the Helmi streams).
- Phase spiral: z vs v_z for stars in a narrow ring around the Sun, coloured by v_φ.

**Check against.** The published figures: GD-1's track and gaps (Price-Whelan & Bonaca 2018), the sausage (Belokurov et al. 2018), the spiral (Antoja et al. 2018).

**Watch out for.** Gaia's archive queue is slow, so design each query tightly and cache it. Radial velocities exist only for brighter stars (G_RVS < 14), so the 6D sample is local. Selecting the halo by velocity shapes what we find, which is a circular argument if we are careless. Energy and angular momentum depend on the assumed potential.

**Builds on.** B1–B3 (Gaia skills), H2 (the rotation curve constrains the potential), B6 (the oldest stars).

## H4. Feel the Galaxy's pull: the Sun's acceleration from quasars ★★☆

**Science question.** The Sun circles the Milky Way at ~230 km/s, so it is always accelerating toward the Galactic Centre. The acceleration is tiny, ~0.2 nanometres per second per second, but Gaia can measure it: it makes every distant quasar appear to drift across the sky in a pattern that points at the Galactic Centre. Can we find this pattern in 1.6 million quasars?

**Background.** Light from any source reaches us tilted by our own velocity (aberration). As our velocity changes, the tilt changes, so all extragalactic sources appear to move in a dipole pattern of proper motions converging on the direction of our acceleration: **secular aberration drift**. For ~2.3 × 10⁻¹⁰ m/s² the amplitude is ~5 µas per year, the width of a hair seen from 4,000 km. Gaia EDR3 measured (2.32 ± 0.16) × 10⁻¹⁰ m/s², or 5.05 ± 0.35 µas/yr, pointing toward the Galactic Centre (Gaia Collaboration, Klioner et al. 2021). Circular motion predicts V₀²/R₀ ≈ (233 km/s)² / 8.2 kpc ≈ 2.1 × 10⁻¹⁰ m/s².

**Data.** Gaia DR3's quasar-like sources of the Gaia-CRF3 (`gaiadr3.agn_cross_id`, ~1.6 million), with their astrometry from `gaia_source`. R₀ and V₀ from H1 and H2.

**Method + visualizations.**
- Quasar proper motions should be zero: their distribution, centred on zero, with widths set by Gaia's errors.
- Fit vector spherical harmonics (the l = 1 "glide" is the acceleration; the l = 1 rotation is a check on the reference frame), with error weights, outlier rejection and a magnitude- and colour-dependent error model.
- The keeper: an all-sky map of binned proper-motion vectors converging on the Galactic Centre.
- Our acceleration against V₀²/R₀ from H1 and H2, and its direction against the Galactic Centre's.

**Check against.** Klioner et al. 2021: (2.32 ± 0.16) × 10⁻¹⁰ m/s², i.e. 5.05 ± 0.35 µas/yr.

**Watch out for.** The signal is hundreds of times smaller than a single quasar's proper-motion error, so only the average over a million sources reveals it. Gaia's own systematics also have large-scale patterns; the fit must separate the glide from frame rotation and magnitude-dependent effects. Some "quasars" are stars: reject sources with significant parallax or proper motion.

**Builds on.** H1 (R₀), H2 (V₀), I3 (quasars), B1–B3 (Gaia astrometry).

## H5. The Milky Way across the spectrum ★☆☆

**Science question.** Our galaxy looks completely different at each wavelength. Radio shows fast electrons spiralling in magnetic fields; 21 cm shows atomic hydrogen; far-infrared shows warm dust; near-infrared shows old stars that visible light cannot reach through the dust; X-rays show million-degree gas; gamma rays show cosmic rays striking gas. Can we put a dozen all-sky maps on one grid and see what each one reveals?

**Background.** NASA's "Multiwavelength Milky Way" posters stacked the Galactic plane in many bands. Each band traces something different:
- synchrotron radio emission (Haslam et al. 1982, 408 MHz);
- the CMB and Galactic foregrounds (WMAP, Planck);
- atomic hydrogen (HI4PI) and molecular gas through CO (Dame et al. 2001, updated 2022);
- dust (IRAS, Planck 857 GHz);
- stars (2MASS near-infrared, optical);
- ionized gas (Hα; Finkbeiner 2003);
- hot gas (ROSAT X-rays);
- gamma rays (Fermi LAT), where the Fermi bubbles reach 50° above and below the Galactic Centre (Su, Slatyer & Finkbeiner 2010).

**Data.** CDS's hips2fits service (`alasky.cds.unistra.fr/hips-image-services/hips2fits`, also wrapped in `astroquery`). It cuts any HiPS survey into an all-sky Mollweide FITS image in Galactic coordinates in a single request. The surveys: Haslam 408 MHz, WMAP K–V, Planck HFI 857 GHz and colour, HI4PI, Dame CO (2022), IRAS/IRIS, 2MASS, Mellinger optical, Finkbeiner Hα, GALEX ultraviolet, ROSAT RASS and Fermi LAT. A few MB each.

**Method + visualizations.**
- The keeper: strips of the Galactic plane (|b| < 10°) in every band, stacked and aligned in longitude: the classic poster, made from data.
- All-sky maps with a slider or a fade animation through the wavelengths.
- Dust against gas: IRAS or Planck dust against HI4PI + CO column density → the gas-to-dust ratio, and "dark gas" that neither hydrogen nor CO traces (Grenier et al. 2005).
- Structures that appear only in some bands: Loop I and the North Polar Spur (radio, X-ray), the Fermi bubbles (gamma rays), Orion and Taurus (CO, dust), the Magellanic Stream (hydrogen).
- RGB composites of three bands (e.g. radio, Hα, X-ray), captioned honestly (G1's rules).

**Check against.** Each survey's own published maps; NASA's multiwavelength Milky Way poster.

**Watch out for.** HiPS maps are resampled; for numbers, use the original HEALPix maps (LAMBDA, the Planck Legacy Archive, HI4PI). Every survey has its own resolution, units and artefacts (zodiacal residuals in Planck, stripes in ROSAT, gaps in IRAS). Colour scales are choices: state the stretch of each panel.

**Builds on.** D3 (Planck maps), H2 (HI4PI), B1 (Gaia star counts), G1 (honest colour), C6 (Fermi).

## H6. Fly through the dust clouds around the Sun ★★☆

**Science question.** Gaia measured the distances and the dimming of 54 million stars near the Sun, enough to map the dust between them in three dimensions. Can we fly through the clouds around us, see the Local Bubble we live in, and trace the Radcliffe Wave, a 2.7 kpc ripple of star-forming clouds?

**Background.** Dust dims and reddens starlight. Knowing the distance and the extinction of many stars along each line of sight places the dust in 3D. Edenhofer et al. (2024) built such a map out to 1.25 kpc: Nside 256 (~14′) on the sky, 516 distance bins, parsec-scale in distance. The Sun sits inside the **Local Bubble**, a cavity a few hundred parsecs across blown out by supernovae over the past ~14 Myr. Most young stars nearby formed on its surface (Zucker et al. 2022). Many nearby clouds (Orion, Perseus, Taurus, Cepheus…) line up along the **Radcliffe Wave**, 2.7 kpc long with an amplitude of ~160 pc (Alves et al. 2020), which is also oscillating (Konietzka et al. 2024).

**Data.** The Edenhofer et al. (2024) map on Zenodo (DOI 10.5281/zenodo.8187942): the mean map and its uncertainty, a 3.3 GB HEALPix FITS file. The record also holds 12 posterior samples (123 GB in all; not needed). The map is also available through the `dustmaps` package. Cloud distances and the Radcliffe Wave model (Zucker et al. 2020; Alves et al. 2020). Young stars from Gaia DR3 (B8).

**Method + visualizations.**
- The keeper: a volume-rendered 3D fly-through (plotly or a WebGL viewer) starting at the Sun and passing Taurus and Orion.
- Face-on and edge-on maps of the dust density, with the walls of the Local Bubble.
- The Radcliffe Wave: fit a damped sine wave to the cloud positions and overlay it on the map.
- Young stars (younger than ~20 Myr, from B8) over the dust: star formation on the bubble's surface.
- Integrated extinction to any distance → a dereddening tool for B1 and B8.

**Check against.** Alves et al. 2020 (length 2.7 kpc, amplitude ~160 pc); Zucker et al. 2022 (the bubble's history); a second map (Bayestar19 or Vergely et al. 2022).

**Watch out for.** Uncertain stellar distances smear dust along lines of sight. The map starts at 69 pc, so the very nearest clouds are missing. Maps from different teams disagree in detail; compare two. Rendering choices (opacity, colour) can create apparent structures.

**Builds on.** B1–B3 (Gaia), B8 (young clusters), H2 (gas), H5 (dust at all wavelengths), L3 (the supernovae that blew the bubble).

---

# I. Black holes

Black holes already run through the atlas: C1 hears two merge, C5 maps the merging population, H1 weighs the one at our Galaxy's centre. This section covers the rest, from ~10 M☉ black holes in X-ray binaries to billion-solar-mass ones powering quasars, to the first images of a black hole's shadow. Two numbers recur. The Schwarzschild radius is `R_s = 2GM/c² ≈ 3 km × (M/M☉)`. The Eddington luminosity, where radiation pressure balances gravity on infalling gas, is `L_Edd ≈ 1.26×10³⁸ (M/M☉) erg/s`.

**Suggested black hole arc:** I4 (stellar) → C5 (merging) → H1 (Sgr A*) → I5 + I6 (the shadow) → I1 (black holes and galaxies) → I2 + I3 (quasars).

## I1. The M–σ relation: black holes grow with their galaxies ★☆☆

**Science question.** A galaxy's central black hole weighs only ~0.5% as much as its central bulge, and its gravity dominates only a tiny region. Yet its mass tracks the random speeds σ of bulge stars thousands of times farther out. How tight is that relation, and what does it imply?

**Background.** The M–σ relation (Ferrarese & Merritt 2000; Gebhardt et al. 2000) is roughly M ∝ σ^4–5 with ~0.3 dex scatter. The black hole's gravity matters only inside its "sphere of influence", `r ≈ GM/σ²`: tens of parsecs, against kiloparsec-sized bulges. The relation therefore implies that black holes and galaxies regulated each other's growth, probably through energy fed back from the black hole during its quasar phase.

**Data.** van den Bosch (2016) compilation of ~230 dynamically measured black hole masses with σ, luminosities and sizes (GitHub `remcovandenbosch/Black-Hole-Mass-compilation`; VizieR `J/ApJ/831/134`). The Kormendy & Ho (2013) tables (ARA&A supplement). The AGN Black Hole Mass Database (Bentz & Katz 2015) for masses from reverberation mapping (I2).

**Method + visualizations.**
- log M vs log σ with errors in both, fitted with intrinsic scatter (MCMC). Points coloured by measurement method (stellar dynamics, gas dynamics, masers, reverberation).
- Classical bulges and ellipticals vs pseudobulges (Kormendy & Ho found that pseudobulges do not follow the relation).
- Our anchors marked: the Milky Way (H1, which sits *below* the relation) and M87 (I5).
- M vs bulge mass. The ~0.5% ratio.

**Check against.** Kormendy & Ho (2013): slope 4.38 ± 0.29, scatter 0.29 dex. McConnell & Ma (2013): slope 5.64. The difference comes from sample choices, so show both.

**Watch out for.** Resolution bias: masses can be measured only where the sphere of influence is resolved, which may push the relation high (Shankar et al. 2016). Masses scale with the assumed distance. σ is defined over different apertures in different papers.

**Builds on.** H1, I5 (anchor points), B7 (σ from spectra). Feeds I3 (black hole demographics).

## I2. Light echoes: weigh a black hole by reverberation mapping ★★☆

**Science question.** The broad emission lines of active galaxies come from gas orbiting close to the black hole, far too small to resolve with any telescope. Can we measure its size from light-travel-time echoes, and so weigh the black hole?

**Background.** The accretion disk's light flickers, and the surrounding broad-line gas brightens in response after a delay τ, the light-travel time. So the gas sits at R = cτ. The line's width gives its orbital speed ΔV, and `M = f · R ΔV² / G`, where f is a geometry factor of order unity calibrated on M–σ. Across many objects, size grows with luminosity as R ∝ L^~0.5 (the **R–L relation**). That relation is what lets a single spectrum estimate the mass of any of 750,000 quasars (I3).

**Data.**
- AGN STORM monitoring of NGC 5548: continuum and Hβ light curves (Pei et al. 2017, ApJ 837, 131, machine-readable tables with the article).
- SDSS Reverberation Mapping final release (Shen et al. 2024, ApJS 272, 26): 849 quasars, 11-year photometric and 7-year spectroscopic light curves, measured lags.
- The AGN Black Hole Mass Database (Bentz & Katz 2015) for the compiled lags and masses.

**Method + visualizations.**
- NGC 5548: continuum and Hβ light curves overlaid, then shifted by the lag (the echo made visible).
- Our own interpolated cross-correlation function (Gaskell & Peterson 1987), with flux-randomization and random-subset-selection resampling for errors (Peterson et al. 1998) → lag → R → M.
- R–L relation from the database and SDSS-RM lags → slope and scatter, which set the precision of every single-epoch mass in I3.
- Our NGC 5548 mass placed on I1's M–σ.

**Check against.** NGC 5548 Hβ lag 4.17 ± 0.36 days (Pei et al. 2017); R–L slope ≈ 0.53 (Bentz et al. 2013).

**Watch out for.** Seasonal gaps and uneven sampling (SDSS-RM's one-year alias). Host-galaxy starlight inflates the measured luminosity. The f factor is uncertain by a factor of ~2 for any single object. During STORM, NGC 5548's Hβ briefly stopped echoing its continuum (the "BLR holiday", Goad et al. 2016). Find it; it is part of the story.

**Builds on.** I1; feeds I3.

## I3. Quasars across cosmic time ★★☆

**Science question.** Quasars are supermassive black holes swallowing gas so fast they outshine their galaxies. When were they most common, how massive did they get, and how did some reach billions of solar masses within a billion years of the Big Bang?

**Background.** Single-epoch black hole masses come from one spectrum: line width + luminosity through I2's R–L relation. The **Eddington ratio** L/L_Edd measures how hard each black hole is feeding. Growing at the Eddington limit with 10% radiative efficiency, a black hole's mass e-folds every ~50 Myr. A seed of 100 M☉ therefore needs ~0.8 Gyr to reach 10⁹ M☉. Yet J0313−1806 at z = 7.64 already holds 1.6 × 10⁹ M☉ only 670 Myr after the Big Bang (Wang et al. 2021). Either the seeds were heavy or growth outran the Eddington limit. JWST's "little red dots" (2023–) have reopened the question.

**Data.** SDSS DR16 quasar catalogue (Lyke et al. 2020; 750,414 quasars). The Wu & Shen (2022, ApJS) catalogue of their spectral properties, black hole masses and Eddington ratios (GitHub `QiaoyaWu/sdss4_dr16q_tutorial`). Published samples of z > 6 quasars (e.g. the Fan, Bañados & Simcoe 2023 review).

**Method + visualizations.**
- Redshift distribution and sky map; luminosity vs redshift, with the survey's flux limit drawn.
- Black hole mass vs redshift with lines of constant Eddington ratio. Quasars pile up below the Eddington limit.
- Number density vs cosmic time for a uniformly selected, luminosity-limited subsample. The quasar era peaks at z ~ 2–3, when the universe was ~2–3 Gyr old, the same era as J3's peak of star formation.
- **The growth problem**: for each z > 6 quasar, the e-folds needed vs the age of the universe at that redshift (D3's cosmology) → the minimum seed mass.
- A composite quasar spectrum (stack thousands, as Vanden Berk et al. 2001 did), linking to B7.

**Check against.** The published quasar luminosity function (e.g. Richards et al. 2006); J0313−1806's numbers.

**Watch out for.** DR16Q's selection is heterogeneous, so number densities need a uniformly selected subsample. Single-epoch masses are uncertain by ~0.4 dex. Different lines serve at different redshifts (Hβ → Mg II → C IV), and outflow-driven shifts in C IV bias its masses. The flux limit means we see only the most luminous quasars at high z.

**Builds on.** I2, B7, D3; pairs with J3.

## I4. Stellar-mass black holes in X-ray binaries ★★☆

**Science question.** The first black holes ever identified sit in binary star systems that shine in X-rays. How are they weighed, and why do they look lighter than the black holes LIGO hears merging?

**Background.** The visible companion star's radial-velocity curve (period P, semi-amplitude K) gives the **mass function** `f(M) = P K³ / (2πG) = M_BH³ sin³i / (M_BH + M★)²`. This is a strict *lower limit* on the unseen mass. If f(M) exceeds ~3 M☉, above any possible neutron star, the object is a black hole whatever the geometry. V404 Cyg has f(M) = 6.08 ± 0.06 M☉. Cygnus X-1 (1972), the first black hole identified, was revised to 21.2 ± 2.2 M☉ in 2021 (Miller-Jones et al.). Known masses show a gap between ~2 and ~5 M☉ separating neutron stars from black holes. Gaia is now finding *dormant* black holes through the wobble of their companions, e.g. Gaia BH1 at ~9.6 M☉ from DR3's orbit solutions.

**Data.** BlackCAT (Corral-Santana et al. 2016; VizieR `J/A+A/587/A61`, with an updated online version): X-ray transients with orbital periods, K, mass functions and masses. Outburst light curves from MAXI (`maxi.riken.jp`) and the Swift/BAT transient monitor. The Gaia DR3 non-single-star catalogue for Gaia BH1. Neutron star masses for comparison (e.g. Özel & Freire 2016).

**Method + visualizations.**
- Mass functions from published P and K → histogram, then masses after marginalizing over inclination. The 2–5 M☉ gap against neutron star masses.
- The combined stellar graveyard: black holes weighed by light (here) vs by gravitational waves (C5). The GW ones are heavier, and part of the reason is that heavier mergers are louder.
- An outburst: the MAXI light curve and hardness–intensity diagram of MAXI J1820+070 or GX 339−4 as the disk switches between hard and soft states.
- Gaia BH1: mass function from its astrometric orbit.

**Check against.** BlackCAT masses; Cyg X-1 21.2 ± 2.2 M☉; V404 Cyg's mass function.

**Watch out for.** Inclination is the dominant systematic, and some argue the mass gap is partly an artefact of it (Kreidberg et al. 2012). We find these black holes only when they are X-ray bright or wobble visibly, which is a strong selection.

**Builds on.** A3 (radial-velocity fitting), C2 (neutron stars), C5 (GW masses).

## I5. The shadow of M87* from Event Horizon Telescope data ★★★

**Science question.** In 2019 the Event Horizon Telescope showed the first image of the glowing ring around a black hole. Can we measure the ring's size ourselves, straight from the calibrated interferometer data, and turn it into a mass?

**Background.** The EHT links radio dishes across Earth at 1.3 mm wavelength, giving a resolution of λ/D ≈ 25 microarcseconds. An interferometer does not take pictures. Each pair of telescopes measures one Fourier component (a "visibility") of the image, at a spatial frequency set by their separation as seen from the source. A thin ring has a characteristic visibility pattern: amplitudes that fall to a first null at a baseline length set by the ring's diameter. The ring encircles the black hole's shadow, whose size is set by the photon orbit: ~10.4 GM/c²D for a non-rotating hole, nearly independent of spin. M87 lies 16.8 Mpc away.

**Data.** EHT 2017 M87 calibrated data (EHT Collaboration 2019, Paper III): uvfits/csv files for four days and two bands, on GitHub (`eventhorizontelescope/2019-D01-01`) and CyVerse. The Sgr A* calibrated data (2022) are also on CyVerse.

**Method + visualizations.**
- The (u, v) coverage: how few Fourier samples the famous image rests on.
- Visibility amplitude vs baseline length: the ring's null, visible in raw data.
- Fit geometric models (thin ring, thick ring, crescent, each blurred by a Gaussian) to visibility amplitudes *and* closure phases. Closure phases sum phases around triangles of telescopes, cancelling each station's phase errors. → ring diameter.
- Diameter → mass, using the EHT's simulation-calibrated ring-to-mass factor (α ≈ 11) and the pure shadow factor (10.4) side by side, to show the systematic.
- Which earlier mass does the ring support: stellar dynamics (6.6 × 10⁹ M☉, Gebhardt et al. 2011) or gas dynamics (3.5 × 10⁹ M☉, Walsh et al. 2013)?
- Sgr A* extension: H1's mass and distance predict the ring size; test that prediction against the EHT's Sgr A* data.
- Optional: image the data with `ehtim` (regularized maximum likelihood).

**Check against.** Ring diameter 42 ± 3 µas; M = (6.5 ± 0.7) × 10⁹ M☉ (EHT 2019); Sgr A* ring 51.8 ± 2.3 µas (EHT 2022).

**Watch out for.** With coverage this sparse, images are model-dependent, which is why geometric fits come first. The bright ring is not the shadow edge itself, and the conversion factor comes from simulations; say so. Station gains are uncertain, so lean on closure quantities.

**Builds on.** C1 (signal processing), H1 (the Sgr A* prediction), I6 (what the ring should look like), I1 (M87 on M–σ).

## I6. Ray-trace a black hole (labeled simulation) ★★☆

**Science question.** What *should* a black hole look like? Can we compute how light bends around one, straight from general relativity, and compare the prediction with the EHT rings?

**Background.** Near a black hole, light can orbit at the **photon sphere** (r = 3GM/c²). Rays passing closer than an impact parameter of √27 GM/c² fall in, which carves out the shadow. A thin accretion disk seen at an angle appears warped: the far side's light is lensed over the top and under the bottom (Luminet's 1979 hand-drawn image, and the look of the film *Interstellar*). Gas moving toward us is Doppler-boosted, so one side of the ring is brighter. **Per ground rule 1 this is a simulation**: nothing is downloaded, and every input (mass, distance, viewing angle) is a measured value from H1 and I5.

**Method + visualizations.**
- Integrate photon paths in the Schwarzschild metric. The orbit equation `d²u/dφ² = −u + 3GMu²/c²` (u = 1/r) is enough. Trace one ray back from each image pixel and record where it hits a thin disk extending from the innermost stable orbit (6GM/c²) outward.
- Apply gravitational redshift and Doppler boosting (bolometric intensity ∝ g⁴) → a Luminet-style image, plus secondary and higher-order images.
- A sweep of viewing angles (animated), and a ray diagram showing the photon sphere.
- Blur to the EHT's ~20 µas resolution at M87's mass, distance and 17° jet inclination → compare ring diameter and brightness asymmetry with I5.
- Extension: spinning (Kerr) black holes from published shadow formulas. The shadow's size changes by only a few percent with spin, which is why it is such a clean mass measurement.

**Check against.** Shadow diameter √27 × 2GM/c² (analytic); Luminet (1979); I5's measured ring.

**Watch out for.** M87's emission comes from a hot, thin, magnetized flow, not a thin disk (the EHT used GRMHD simulations), so our image is illustrative and must be labelled as such. Numerical accuracy degrades for rays skimming the photon sphere.

**Builds on.** I5, H1, C1 (general relativity in action).

*I7–I14 follow black holes as they merge, from stellar pairs heard by LIGO to supermassive pairs heard by pulsars, and hunt for black holes that do not shine.*

## I7. Hawking's area theorem: hear a black hole ring ★★★

**Science question.** GW250114 is the loudest black hole merger recorded so far. After the merger, the new black hole rang like a struck bell for a few milliseconds. Can we pick out its tones in the public data, read its mass and spin from them, and test Hawking's 1971 theorem that the total area of black hole horizons can never shrink?

**Background.** A disturbed black hole sheds its distortions as **quasinormal modes**: damped sinusoids whose frequencies and decay times depend only on the final mass and spin (the "no-hair" property of Kerr black holes). The loudest is the fundamental (ℓ, m, n) = (2, 2, 0) tone; its overtones (n = 1, 2, …) die away faster. A Kerr black hole's horizon area is `A = 8π (GM/c²)² (1 + √(1 − χ²))`, where χ is the dimensionless spin. The area theorem says the remnant's area must exceed the sum of the two initial areas. GW250114 (14 January 2025; 33.6 + 32.2 M☉ → ~63 M☉, ~400 Mpc away) reached a network SNR of 80, three times GW150914's 26. The LVK confirmed the area increase with ~99.999% credibility, found the first overtone, and constrained the tones to within ±30% of the Kerr spectrum (LVK 2025, arXiv:2509.08054). GW150914's ringdown alone had given 97% (Isi et al. 2021).

**Data.** GWOSC event `GW250114_082203`: LIGO Hanford and Livingston strain (32 s and 4096 s, at 4 and 16 kHz). Parameter-estimation samples from the GWTC-5.0 release on Zenodo (linked from the event page) give the inspiral masses and spins. Kerr mode frequencies come from the `qnm` package (Stein 2019).

**Method + visualizations.**
- Whitened strain and Q-transform (C1's code): the chirp, then a brief ring.
- Ringdown-only fits started at a series of times after the peak: one damped sinusoid (220), then 220 + 221. Fit in the time domain with the noise covariance from the PSD (Isi et al. 2019). A short segment that starts abruptly breaks FFT-based methods. MCMC (`emcee`) → frequency and damping time → final mass and spin.
- The "black hole spectroscopy" plane: measured frequency and damping time over the Kerr prediction.
- Areas: A₁ + A₂ from the inspiral samples, the final area from the ringdown. A histogram of the fractional change, against the predicted ~65% increase.
- Audio of the ringdown, shifted up to an audible pitch.

**Check against.** LVK (arXiv:2509.08054): area increase at ~99.999% credibility; mode frequencies within ±30% of Kerr. GWOSC (GWTC-5.0): final mass 62.9 (+1.0/−1.1) M☉, network SNR 78.6.

**Watch out for.** The ringdown's start time is the main systematic. Start too early and the nonlinear merger leaks in, where overtones can soak it up and look significant; start too late and the signal is gone. Show every result against start time. Detector-frame masses are (1 + z) ≈ 1.09 times the source-frame ones. The `ringdown` package needs Python < 3.13, so on this machine we write the time-domain likelihood ourselves; it is short. The initial areas come from GR waveforms, so this is a consistency test of GR, not an independent proof; say so.

**Builds on.** C1 (whitening, Q-transform), C5 (the population), I6 (the photon sphere: roughly, the tones are light trapped there leaking away).

## I8. The nanohertz hum: supermassive black hole binaries in pulsar timing ★★☆ (own noise analysis ★★★)

**Science question.** In 2023 pulsar-timing arrays reported evidence for gravitational waves with periods of years. They most likely come from pairs of supermassive black holes spiralling together across the universe. Can we see the evidence ourselves: the spectrum of this hum, and the Hellings–Downs pattern that marks it as gravitational waves?

**Background.** Millisecond pulsars (C2) keep time to ~100 ns over decades. A passing gravitational wave changes the light-travel time between each pulsar and Earth, which shifts the pulses' arrival times. The part of the shift caused at Earth correlates pairs of pulsars in a way that depends only on their angle ζ on the sky. This is the **Hellings–Downs curve** (1983): `Γ(ζ) = ½ − ¼x + (3/2) x ln x`, with `x = (1 − cos ζ)/2`. A clock error would correlate all pulsars equally (a monopole), and an error in the Solar System model gives a dipole; only gravitational waves give this quadrupolar shape. Binaries of 10⁸–10¹⁰ M☉ black holes left by galaxy mergers produce a background with characteristic strain `h_c(f) = A (f / 1 yr⁻¹)^(−2/3)`, i.e. timing-residual power ∝ f^(−13/3). NANOGrav's 15-year data (67 pulsars) prefer Hellings–Downs correlations over uncorrelated noise with Bayes factors of 200–1000 (~3–4σ), with A = 2.4 (+0.7/−0.6) × 10⁻¹⁵ (Agazie et al. 2023). The European, Indian, Parkes, Chinese and MeerKAT arrays reported consistent signals. The International PTA is combining them into DR3 (~130 pulsars).

**Data.**
- The NANOGrav 15-year data set v2.1.0 (Zenodo, DOI 10.5281/zenodo.16051178, 639 MB): arrival times and timing models for 68 pulsars.
- Kernel-density representations of NANOGrav's free-spectrum posteriors (Zenodo record 10344086, 5 MB), which refit any spectral model in seconds with `ceffyl`.

**Method + visualizations.**
- Quick route: power-law fits to the free spectrum, with γ free and with γ = 13/3 → A. The spectrum with the supermassive-binary expectation drawn on.
- Full route (★★★): load the par/tim files with PINT (`pint-pulsar`). Model each pulsar's white noise, red noise and dispersion-measure variations with `enterprise`, then fit the common process. Then compute the **optimal statistic** (Anholm et al. 2009; Chamberlin et al. 2015): the 2,211 pulsar-pair correlations, binned by angle, over the Hellings–Downs curve, with monopole and dipole curves for contrast.
- A sky map of the 67 pulsars; timing residuals of a well-timed pulsar such as PSR J1909−3744.
- What A means: the number and masses of binaries it requires (Phinney 2001's energy argument), compared with I1 and I3's black hole demographics.

**Check against.** Agazie et al. 2023 (ApJL 951, L8): A = 2.4 (+0.7/−0.6) × 10⁻¹⁵ at f = 1/yr for γ = 13/3; Hellings–Downs correlations favoured with Bayes factors of 200–1000.

**Watch out for.** Pulsar noise models decide the answer: one pulsar's red noise can masquerade as a common signal. The detection claim rests on the correlations, not on the common spectrum. Pulsar distances are poorly known, so the optimal statistic uses only the Earth term. On Python 3.14, `enterprise-pulsar` and PINT are pure Python but untested here; `PTArcade` needs Python < 3.13.

**Builds on.** C2 (millisecond pulsars), C1 (gravitational waves), I1 and I3 (supermassive black hole demographics), J7 (galaxy mergers).

## I9. Dark sirens: the expansion rate from black hole mergers ★★★

**Science question.** C4 measured H0 from the one neutron star merger whose host galaxy was seen. Black hole mergers give off no light, but their signals still tell us how far away they are. Can hundreds of them, together with a catalogue of galaxies, measure the expansion rate of the universe?

**Background.** A merger's signal amplitude gives its luminosity distance directly, with no distance ladder: a **standard siren** (Schutz 1986). For the redshift, a "dark siren" analysis treats every galaxy in the event's 3D localization volume as a possible host. Over many events the true hosts' redshifts add up while the wrong ones average out. A "spectral siren" uses the black hole mass distribution instead. Masses reach us redshifted by (1 + z), so matching the observed peaks to C5's source-frame features near 10 and 35 M☉ gives z. GWTC-5.0 (May 2026) combined 235 dark sirens with GW170817 to get H0 = 71.0 (+9.0/−7.1) km/s/Mpc. For the first time the dark sirens alone constrained H0 better than the single bright siren.

**Data.**
- Sky maps and parameter-estimation samples from GWTC-2.1, GWTC-3, GWTC-4.0 and GWTC-5.0 (Zenodo, linked from gwosc.org), plus the O4a + O4b search-sensitivity injections for selection effects.
- GLADE+ (Dálya et al. 2022): 22.5 million galaxies and 750,000 quasars. The full text file is 6.4 GB. VizieR (`VII/291`) serves cone searches, so fetch only the galaxies inside each event's 90% sky area.

**Method + visualizations.**
1. 3D localization volumes for the best-localized events (GW190814 and GW170817 first), with the GLADE+ galaxies inside them.
2. For each event, an H0 likelihood that sums over candidate hosts weighted by luminosity, plus a term for the galaxies the catalogue misses. GLADE+ is complete only to ~44 Mpc (in B-band luminosity).
3. Selection effects: at a higher H0, the same detectable distances hold more galaxies. Normalize with the sensitivity injections (Mandel, Farr & Gair 2019).
4. Multiply the events' posteriors and show the result narrowing as events are added. Add GW170817's bright siren from C4, and draw D1's 73.5 and D3's 66.7.
5. A spectral-siren cross-check with a fixed mass model.

**Check against.** LVK GWTC-5.0 cosmology (arXiv:2605.27227): H0 = 71.0 (+9.0/−7.1) km/s/Mpc. GWTC-4.0 (arXiv:2509.04348): 75.4 (+12.8/−9.1) from 141 dark sirens plus GW170817. The LVK's `gwcosmo` code (on PyPI) for a like-for-like run.

**Watch out for.** Most events lie hundreds of Mpc away, where GLADE+ is very incomplete; without a completeness correction the answer is biased. Even catalogue-based dark sirens depend on the assumed mass distribution. GLADE+'s photometric redshifts have large errors. Single-event posteriors are nearly flat, so the information is in the product; check that no single event drives it.

**Builds on.** C4 (bright siren), C5 (mass distribution), D1 and D3 (the Hubble tension), D2 (galaxy catalogues).

## I10. Tidal disruption events: black holes caught eating stars ★★☆

**Science question.** Most supermassive black holes are dark and quiet because they have no gas to eat. Now and then a star strays too close and is torn apart, and the black hole flares for months. Can we find these flares in public survey data, test the predicted t^(−5/3) fade, and weigh black holes that are otherwise invisible?

**Background.** A star is pulled apart where the black hole's tides beat its own gravity, at `r_t ≈ R★ (M_BH / M★)^(1/3)`. Above ~10⁸ M☉ (the Hills mass, for a Sun-like star) r_t lies inside the horizon, and the star is swallowed whole without a flare. About half the debris stays bound and falls back at a rate ∝ t^(−5/3) (Rees 1988; Phinney 1989). The optical flares are blue and hot (a few × 10⁴ K), last months, and sit in galaxy nuclei. The Zwicky Transient Facility found 33 in a flux-limited sample over three years (Yao et al. 2023), in black holes of 10^5.1–10^8.2 M☉. In galaxies of ~10¹⁰ M☉ the rate is ≈3.2 × 10⁻⁵ per galaxy per year, about once every 30,000 years. It is highest in "green valley" galaxies (J3) and drops above ~10^7.5 M☉, as direct capture predicts. Many TDEs later settle onto a plateau whose brightness tracks the black hole mass (Mummery et al. 2024, 63 TDEs).

**Data.** ZTF photometry from the ALeRCE broker API (`api.alerce.online`; no account needed). Classifications from the Transient Name Server and the published samples (Yao et al. 2023; Hammerstein et al. 2023). Host galaxies from SDSS or the Legacy Survey, with velocity dispersions from B7 for I1's M–σ masses.

**Method + visualizations.**
- Light curves of Yao et al.'s 33 events aligned at peak. Fit the decline with (t − t₀)^(−p): how many give p ≈ 5/3?
- Colour and blackbody temperature against time: TDEs stay hot while supernovae cool and redden. That colour is how surveys tell them apart.
- Peak luminosity and decay time against black hole mass from M–σ, and the cutoff near the Hills mass.
- Late-time plateaus → black hole masses (Mummery et al. 2024's relation), compared with the M–σ masses.
- The rate per galaxy from the sample's survey volume (V_max). Then a live TDE from the last month's TNS classifications, fitted the same way.

**Check against.** Yao et al. 2023: ≈3.2 × 10⁻⁵ galaxy⁻¹ yr⁻¹ at M_gal ~ 10¹⁰ M☉; (3.7, 7.4, 1.6) × 10⁻⁵ in red, green and blue galaxies; black hole masses 10^5.1–10^8.2 M☉.

**Watch out for.** The t^(−5/3) law describes the fallback rate, not necessarily the optical light, which may come from shocks or reprocessing; fitted indices often differ. Nuclear supernovae and flaring AGN contaminate the sample. M–σ masses carry ~0.3–0.4 dex scatter. A flux limit favours luminous flares, which may favour certain black hole masses.

**Builds on.** I1 (masses), I3 (feeding black holes), D6 (ZTF and ALeRCE), B7 (TDE spectra: broad hydrogen and helium lines), J3 (green-valley galaxies).

## I11. Dormant black holes in Gaia's binary orbits ★★☆

**Science question.** The Milky Way should hold ~100 million stellar-mass black holes, and nearly all are dark: nothing feeds them, so they emit no X-rays. Gaia has found three by the wobble they give a companion star. Can we repeat the hunt in Gaia DR3's orbit catalogue, and be ready for DR4 on 2 December 2026?

**Background.** Gaia sees a star trace a tiny ellipse on the sky as it orbits an unseen companion. If the companion is dark, the ellipse's angular size a₀, the parallax ϖ and the period P give the **astrometric mass function** `(a₀/ϖ)³ / P² = M₂³ / (M₁ + M₂)²`, in AU, years and M☉. With the star's own mass M₁ from its colour and luminosity, solve for M₂. A companion of more than ~3 M☉ that gives off no light can only be a black hole. Three are known:
- Gaia BH1 (El-Badry et al. 2023): 9.6 M☉ in a 186-day orbit around a Sun-like star, 480 pc away;
- Gaia BH2: 8.9 M☉ in a 1,277-day orbit around a red giant, 1.16 kpc away;
- Gaia BH3 (Gaia Collaboration 2024): 32.7 M☉ in an 11.6-year orbit, 590 pc away, found in pre-release DR4 data. Its companion is a very metal-poor giant from a disrupted star cluster.

The same method has turned up ~20 neutron star candidates in wide orbits (El-Badry et al. 2024).

**Data.** Gaia DR3 `gaiadr3.nss_two_body_orbit`: 443,205 solutions, ~169,000 of them astrometric (`Orbital`, `AstroSpectroSB1` and variants). Join them to `gaia_source` and the astrophysical parameters (FLAME masses). The data and code released with Gaia BH3 ("Gaia Black Hole 3 tools", ESA). From 2 December 2026, DR4's epoch astrometry for every source.

**Method + visualizations.**
- Astrometric orbits with significant a₀ and good parallaxes → mass function → M₂ for each, given M₁. The orbits come as Thiele–Innes constants; convert them to a₀ first.
- M₂ against period for all solutions: ordinary binaries, white dwarf companions and the outliers. Recover Gaia BH1 and BH2 blind.
- Vet candidates as El-Badry et al. did. Is there light from the companion in the spectral energy distribution? Could a hidden inner binary (a triple) explain it? Is the solution poor (quality flags, RUWE)?
- Gaia BH1's orbit drawn on the sky, to scale.
- Why DR3 could not see BH3: its 11.6-year period is four times DR3's 34-month baseline.
- The stellar graveyard with C5 and I4: BH3 sits among LIGO's black holes, not among the X-ray binaries.

**Check against.** El-Badry et al. 2023a,b: BH1 9.62 ± 0.18 M☉, BH2 8.94 ± 0.34 M☉. Gaia Collaboration (Panuzzo et al.) 2024: BH3 32.70 ± 0.82 M☉.

**Watch out for.** Most DR3 orbits with "massive dark companions" are spurious, so check every candidate against the quality flags and the published validation. A dark companion can be a close pair of faint stars or white dwarfs, which inflates M₂. DR3's 34 months limit periods to under ~3 years. M₁ comes from stellar models, so propagate its uncertainty.

**Builds on.** I4 (mass functions, X-ray binaries), A3 (Keplerian orbits), B1–B3 (Gaia queries), H3 (BH3's companion came from the disrupted cluster ED-2).

## I12. A lone black hole found by microlensing ★★☆

**Science question.** A black hole with no companion emits nothing and orbits nothing. In 2011 one drifted in front of a star in the Galactic bulge, magnifying it for months and nudging its apparent position by a couple of milliarcseconds. Can we fit the event ourselves and weigh the only isolated black hole found so far?

**Background.** When a dark lens passes close to our line of sight to a star, the star brightens and fades symmetrically (the Paczyński curve). The magnification is `A(u) = (u² + 2) / (u √(u² + 4))`, where u is the lens–star separation in units of the **Einstein radius** θ_E (J5's lensing, on the scale of a single star). The event's timescale t_E = θ_E / μ_rel mixes up mass, distance and speed. Two extra measurements break the degeneracy. The **microlens parallax** π_E comes from the light curve's distortion as Earth orbits the Sun. θ_E itself comes from the shift in the star's apparent position, which is largest, θ_E/√8, at u = √2. Then `M = θ_E / (κ π_E)` with κ = 8.144 mas/M☉. For OGLE-2011-BLG-0462, HST measured the shift over 11 years: 7.15 ± 0.83 M☉ at 1.52 kpc, with no light from the lens (Sahu et al. 2025). Independent analyses at first disagreed: Lam et al. (2022) found 1.6–4.4 M☉, which allowed a neutron star. Systematic errors in the photometry and astrometry explained the gap.

**Data.** OGLE Early Warning System photometry of OGLE-2011-BLG-0462 (a public tarball on the EWS event page, 0.85 MB). The HST WFC3 images are public in MAST (DOI 10.17909/v97n-kh49); Sahu et al. (2022, 2025) plot the per-epoch deflections but do not tabulate them. Thousands of other EWS bulge events per season for comparison. The Gaia DR3 microlensing catalogue (`gaiadr3.vari_microlensing`; Wyrzykowski et al. 2023) for long events across the whole sky.

**Method + visualizations.**
- Fit a point-lens model to the OGLE light curve (t₀, u₀, t_E, blending), then add annual parallax (π_E,N, π_E,E). The long event becomes asymmetric.
- The t_E distribution for one season of EWS events. Our event sits in the long tail, where heavy lenses live (t_E ∝ √M).
- The astrometric track: the star's apparent position looping on the sky during the event, with the HST measurements.
- A mass–distance diagram, where the π_E and θ_E constraints cross at the lens.
- Why a 7 M☉ *star* at 1.5 kpc would have been easy to see, and was not.

**Check against.** Sahu et al. 2025: 7.15 ± 0.83 M☉, 1.52 ± 0.15 kpc, moving at 51 ± 8 km/s relative to its surroundings.

**Watch out for.** Parallax from ground-based photometry is subtle and easily biased by systematics, which caused the early disagreement. The source has a bright neighbour 0.4″ away, so its astrometry needs careful PSF subtraction. Blended light dilutes the magnification. Parallax fits have a sign degeneracy (u₀ > 0 or < 0), so fit both.

**Builds on.** J5 (Einstein radius), I4 and I11 (black holes in binaries), K4 (the same technique turned on dark matter).

## I13. An intermediate-mass black hole in ω Centauri? ★★☆

**Science question.** Black holes come small (5–150 M☉) and supermassive (over 10⁵ M☉), but the middle of the range is nearly empty. In 2024, seven stars racing through the core of ω Centauri pointed to a black hole of at least 8,200 M☉. Can we find those stars in the public Hubble catalogue, and test whether something else could hold them in?

**Background.** ω Cen is the Milky Way's most massive globular cluster, probably the stripped nucleus of a dwarf galaxy the Milky Way swallowed. Its central escape velocity is ~62 km/s. A star moving faster would leave the core within about a thousand years unless something massive and compact holds it. Häberle et al. (2024, *Nature*) found seven such stars within 3″ (0.08 pc) of the centre in 20 years of HST astrometry; the fastest crosses the sky at 113 km/s. They imply a black hole of at least ~8,200 M☉, and models allow up to ~50,000 M☉. A 2025 analysis that adds millisecond-pulsar accelerations (Bañares-Hernández et al.) prefers an extended cluster of dark stellar remnants (~2–3 × 10⁵ M☉) instead, with an IMBH below 6,000 M☉ (3σ). The question is open.

**Data.** oMEGACat II (Häberle et al. 2024, ApJ 970, 192): HST proper motions of 1,395,781 stars within ω Cen's half-light radius over a ~20-year baseline (Zenodo, DOI 10.5281/zenodo.11104046: astrometry 448 MB, plus seven photometric catalogues). The seven fast stars are listed in the *Nature* paper's Extended Data Table 1. oMEGACat I (Nitschai et al. 2023): MUSE line-of-sight velocities. Distance 5.43 kpc (Baumgardt & Vasiliev 2021).

**Method + visualizations.**
- Proper motions in the central 20″, converted to km/s (v = 4.74 µ d). Find the fast stars blind: speed against radius, compared with the 62 km/s escape speed, keeping only well-measured stars on the cluster's main sequence.
- A minimum enclosed mass for each star: a bound star at projected radius R with 2D speed v needs `M > v² R / 2G`. Both projections make this a lower limit.
- The velocity-dispersion profile from the whole catalogue, and Jeans models with and without a central point mass.
- How the answer moves with the adopted cluster centre; published centres differ by arcseconds.

**Check against.** Häberle et al. 2024: seven stars above 62 km/s; star A at 4.41 ± 0.08 mas/yr (113.0 ± 1.1 km/s); M ≥ ~8,200 M☉. Bañares-Hernández et al. 2025: below 6,000 M☉ (3σ).

**Watch out for.** Crowding: blends in the dense core can fake fast motions, so require many epochs and clean fits. The escape velocity depends on the unseen population of stellar-mass black holes and neutron stars, which is the 2025 counter-argument. Exclude foreground stars with the colour–magnitude diagram. Most stars have only 2D motions.

**Builds on.** H1 (weighing a central mass with stars), B3 and B6 (globular clusters, membership), I1 (where an IMBH would sit on M–σ), H3 (ω Cen as a swallowed galaxy's nucleus).

## I14. Megamasers: weigh a black hole, measure the universe ★★★

**Science question.** In a few galaxies, water molecules in a thin disk around the central black hole act as natural microwave lasers (masers). Radio interferometers can follow individual maser spots for years. Can we use them to weigh the black hole in NGC 4258 to a few percent, and measure that galaxy's distance by geometry alone?

**Background.** Masers at the disk's edges, seen side-on, move toward and away from us at up to ~1,000 km/s. Their speeds fall with radius as r^(−1/2), following Kepler's law, so ~4 × 10⁷ M☉ is packed inside 0.13 pc (Miyoshi et al. 1995; Humphreys et al. 2013). After Sgr A* (H1), this is the strongest case for a supermassive black hole. Masers on the disk's near side, in front of the black hole, speed up along our line of sight by about 9 km/s per year (centripetal acceleration, a = v²/r) and drift across the sky. The acceleration gives the disk radius in physical units and the VLBI map gives it as an angle, so their ratio is the distance: `D ≈ v² / (a θ)`. NGC 4258 is 7.576 ± 0.082 ± 0.076 Mpc away (Reid, Pesce & Riess 2019) and anchors the Cepheid ladder (B4). For galaxies 50–150 Mpc away, the same method measures H0 directly. The Megamaser Cosmology Project's six galaxies give H0 = 73.9 ± 3.0 km/s/Mpc (Pesce et al. 2020).

**Data.** For NGC 4258: VLBI positions and velocities from 18 epochs (Argon et al. 2007, VizieR `J/ApJ/659/1040`) and accelerations from 10 years of monitoring (Humphreys et al. 2008). Megamaser Cosmology Project maps and accelerations, e.g. NGC 5765b (Gao et al. 2016, `J/ApJ/817/128`) and CGCG 074-064 (Pesce et al. 2020, `J/ApJ/890/118`). The six distances from Pesce et al. 2020.

**Method + visualizations.**
- NGC 4258's position–velocity diagram: two Keplerian wings and the straight line of the near-side masers → the black hole mass and the enclosed-mass profile.
- Maser spots on the sky coloured by velocity: a warped thin disk seen edge-on.
- The distance: the simple estimate D ≈ v²/(aθ) from the near-side masers, then a warped-disk fit to positions, velocities and accelerations.
- H0 from the six MCP galaxies on a Hubble diagram with D1's supernovae, and the ladder re-anchored on our own NGC 4258 distance.
- The black hole on I1's M–σ relation.

**Check against.** Reid et al. 2019: 7.576 ± 0.082 (stat) ± 0.076 (sys) Mpc. Humphreys et al. 2013: 7.60 ± 0.17 ± 0.15 Mpc. Pesce et al. 2020: H0 = 73.9 ± 3.0. Gao et al. 2016: NGC 5765b's enclosed mass (4.55 ± 0.40) × 10⁷ M☉.

**Watch out for.** The disk is warped, and the orbits may be slightly eccentric (Humphreys et al. 2013 fit both). Matching spots across epochs is a tracking problem. MCP galaxies have peculiar velocities of a few hundred km/s, a few percent of their recession speeds; correct them with a flow model (J4). Six galaxies is a small sample, so the H0 error is mostly statistical.

**Builds on.** H1 (Keplerian black hole masses), B4 (distance ladder), I1 (M–σ), J4 (peculiar velocities), D1 (H0).

*I15–I17 return to quasars: as backlights for the young universe, as jets that seem to outrun light, and as flickering accretion discs.*

## I15. Quasars as backlights: the Lyman-α forest and the end of reionization ★★☆

**Science question.** Every cloud of hydrogen between us and a distant quasar absorbs a little of its light at the Lyman-α wavelength, shifted by the cloud's own redshift. The result is a "forest" of absorption lines that maps the cosmic web along the line of sight. Toward the most distant quasars, the forest turns into near-total darkness. When did the universe's hydrogen become transparent?

**Background.** Neutral hydrogen absorbs at 121.6 nm in its rest frame, so gas at redshift z absorbs at 121.6 (1 + z) nm. Blueward of the quasar's own Lyα emission line, the spectrum records the gas along the whole sightline. Even at z ≈ 3 only about one hydrogen atom in 100,000 is neutral, yet that is enough to make the forest. Gunn & Peterson (1965) pointed out that fully neutral gas would absorb everything, leaving a dark "trough". Quasars at z > 6 show long dark troughs (Becker et al. 2001; Fan et al. 2006). The XQR-30 spectra show large sightline-to-sightline scatter at z ≥ 5.4: reionization was still patchy then and ended by z ≈ 5.3, ~1.1 Gyr after the Big Bang (Bosman et al. 2022). The CMB (Planck's optical depth τ ≈ 0.054) puts the midpoint near z ≈ 7.7.

**Data.** The XQR-30 reduced spectra (D'Odorico et al. 2023; GitHub `XQR-30/Spectra`): 42 high-quality X-shooter spectra of quasars at z ≈ 5.8–6.6. SDSS DR16/17 quasar spectra (I3) for the forest at z ≈ 2–4.

**Method + visualizations.**
- A gallery of quasars from z = 2 to 6.5 aligned at Lyα: the forest thickens and turns black.
- Effective optical depth τ_eff = −ln⟨F⟩ in 50 Mpc/h chunks against z, from SDSS (z ≈ 2–4) and XQR-30 (z ≈ 5–6): the steep rise above z ≈ 5.5, and the scatter between sightlines.
- Dark gaps: the longest stretches with no transmitted light, sightline by sightline.
- The 1D flux power spectrum at z ≈ 3 from SDSS spectra: the cosmic web in absorption.
- Proximity zones: the size of the ionized bubble each z > 6 quasar carved around itself.

**Check against.** Bosman et al. 2022: models with a uniform ultraviolet background fit the data up to z ≈ 5.2 and are ruled out at z ≥ 5.4; reionization ends by z = 5.3. Becker et al. 2013: τ_eff(z) at z = 2–5.

**Watch out for.** The quasar's unabsorbed continuum is hidden under the forest and must be estimated. At low z, continuum fitting is the main systematic; at high z the continuum is predicted from the red side of the spectrum. Metal lines and damped Lyα absorbers contaminate. In the darkest troughs, sky-subtraction residuals dominate. Exclude the proximity zone from forest statistics.

**Builds on.** I3 (quasars), B7 (spectra), D2 (the cosmic web), J8 (the galaxies that did the reionizing), D3 (τ from the CMB).

## I16. Faster than light? Superluminal jets ★★☆

**Science question.** Radio images of quasar jets show blobs that appear to cross the sky at up to ~50 times the speed of light. Nothing moves faster than light, so what is going on, and what does it tell us about the jets launched near black holes?

**Background.** A blob moving at speed βc at a small angle θ to our line of sight nearly keeps pace with its own light, so the images we receive are compressed in time. Its apparent sideways speed is `β_app = β sin θ / (1 − β cos θ)`. This exceeds 1 for β close to 1, reaching a maximum of βγ at cos θ = β (predicted by Rees 1966). The same geometry boosts the jet's brightness (Doppler beaming), so jets aimed at us dominate bright samples (blazars). MOJAVE (Monitoring Of Jets in AGN with VLBA Experiments) has imaged hundreds of jets at 15 GHz since 1994 and tracked 1,744 features in 382 jets over at least five epochs each. Apparent speeds peak below ~10c and reach ~50c in PKS 0805−07 (Lister et al. 2016, 2019).

**Data.** MOJAVE kinematics tables (Lister et al. 2019, VizieR `J/ApJ/874/43`; Lister et al. 2021, `J/ApJ/923/30`), plus the project's online velocity table and VLBA images (FITS) at `cv.nrao.edu/MOJAVE`. Redshifts from the same tables.

**Method + visualizations.**
- An animated jet: MOJAVE images of 3C 279 or 3C 273 over two decades, with blobs moving outward.
- Our own fit of one feature's position against time → angular speed → β_app, using its redshift and D3's cosmology.
- The β_app distribution for all jets. Simulate a population with a range of Lorentz factors and random orientations, Doppler-boosted and flux-limited, and fit the observed distribution → the intrinsic Lorentz factors.
- β_app against brightness: are the fastest jets the most boosted?
- Accelerating features: how many jets speed up or slow down along their length?

**Check against.** Lister et al. 2019: speeds peaking below ~10c with a tail to ~50c, Lorentz factors in the parent population up to ~50, and 59% of the best-sampled features accelerating.

**Watch out for.** Features are fitted Gaussian blobs, and matching the same blob across epochs is a choice. Apparent speeds scale with the distance assumed for each redshift. Flux-limited samples favour boosted jets, so the observed distribution is not the intrinsic one; that is the point of the simulation. Some features are standing shocks that do not move.

**Builds on.** I3 (quasars), I5 (M87's jet base at the EHT), C6 (relativistic beaming in gamma-ray bursts), D3 (distances).

## I17. The flickering of quasars ★☆☆

**Science question.** Quasars brighten and fade irregularly by tens of percent over months and years. Is this flicker just noise, or does it carry information about the black hole and its accretion disc? What do ~9,000 quasars watched for 10 years tell us?

**Background.** Quasar light curves are well described by a **damped random walk**: a random process that wanders but is pulled back toward its mean on a characteristic timescale τ (Kelly et al. 2009; MacLeod et al. 2010), typically hundreds of days. Its structure function rises with time lag and flattens beyond τ. τ grows with black hole mass across the whole range from dwarf galaxies to the largest quasars, and matches the thermal timescale of the accretion disc. Accreting white dwarfs fall on the same relation (Burke et al. 2021). Some quasars change dramatically: "changing-look" quasars lose or gain their broad emission lines within years (LaMassa et al. 2015).

**Data.** SDSS Stripe 82 quasar light curves: 9,258 quasars with ~10 years of ugriz photometry, plus black hole masses (MacLeod et al. 2010; `faculty.washington.edu/ivezic/macleod/qso_dr7/`). ZTF light curves (IRSA or ALeRCE) to extend the baseline past 2018. I3's DR16Q catalogue.

**Method + visualizations.**
- Light curves of a few quasars across 20 years (SDSS then ZTF).
- A damped-random-walk fit for each quasar (Gaussian-process likelihood with `celerite2`) → τ and amplitude.
- The ensemble structure function, and how τ and amplitude depend on luminosity, wavelength and black hole mass; the Burke et al. relation.
- Bluer when brighter: colour against brightness.
- A hunt for changing-look candidates: quasars whose brightness changed by more than a magnitude between SDSS and ZTF.

**Check against.** MacLeod et al. 2010 (τ distribution); Burke et al. 2021 (τ against black hole mass).

**Watch out for.** τ is biased low unless the light curve spans at least ~10τ (Kozłowski 2017). Seasonal gaps. SDSS and ZTF photometry differ in calibration and colour terms. Host-galaxy light dilutes the variability of faint AGN.

**Builds on.** I2 (reverberation: the same flicker, echoed), I3 (black hole masses), D6 (ZTF data), I10 (nuclear transients).

## I18. Giant radio galaxies: the largest structures black holes build ★★☆

**Science question.** Some supermassive black holes launch jets that run far beyond their galaxy and inflate radio lobes millions of light-years across. LOFAR's sky survey found that such "giants" are common, and in 2024 one of them, Porphyrion, turned out to span 7 Mpc, the longest jets known. Can we find giants in the radio maps, measure them, and see them over images of their host galaxies?

**Background.** At low radio frequencies the old, extended lobes of radio galaxies are brightest. LOFAR's LoTSS DR2 imaged 27% of the northern sky at 144 MHz with 6″ resolution and catalogued 4,396,228 sources (Shimwell et al. 2022). A radio galaxy counts as a "giant" if it spans more than 0.7 Mpc. Machine learning plus citizen science raised the number of known giants above 11,000 (Oei et al. 2024). Porphyrion (Oei et al. 2024, *Nature*) spans 7 Mpc, from a galaxy seen when the universe was ~6.3 Gyr old. The lengths are projections, so true sizes are larger. The population tells us how long jets stay switched on and how the cosmic web confines them.

**Data.** LoTSS DR2 images and source catalogue (`lofar-surveys.org`: mosaics as FITS, plus a cutout service). The giant radio galaxy catalogues of Oei et al. (2023, 2024). Host galaxies and redshifts from the Legacy Survey and SDSS (G1, B7). NVSS at 1.4 GHz for spectral indices.

**Method + visualizations.**
- The keeper: a gallery of giants as radio contours over Legacy Survey colour images, with the Milky Way and Andromeda drawn to scale beside them.
- Angular sizes → projected lengths with D3's cosmology; the length distribution against Oei et al.'s.
- Spectral index from LoTSS's in-band images, or against NVSS → the lobes' ages from synchrotron ageing.
- Where the hosts sit on I1 and I3's black hole mass scale, and where the giants sit in D2's cosmic web.
- Porphyrion in full, with its host galaxy marked.

**Check against.** The length distributions of Oei et al. (2023, 2024); Porphyrion's 7 Mpc.

**Watch out for.** Matching separate lobes to one host is hard, which is why the catalogues needed citizen scientists. LoTSS's resolution and surface-brightness limit miss the faintest, largest lobes. Projected lengths underestimate true sizes.

**Builds on.** I16 (jets at parsec scale), I3 (quasars), D2 (the cosmic web), G1 (colour images).

---

# J. Galaxies

D2 mapped where galaxies are and D4 showed what holds them together. This section asks what galaxies *are*: how they rotate, what shapes they take, how they live and die, where they are falling, and what the first ones looked like.

**Suggested galaxy arc:** J2 (shapes) → J1 (rotation) → J5 (weighing with lensing) → J3 (life cycle) → J4 (flows) → J7 (our own fate).

## J1. Rotation curves from raw velocity fields ★★☆

**Science question.** D4 took its rotation curves from the SPARC tables. Can we measure them ourselves, from maps of how gas moves across a galaxy's face, and see what the circular-rotation picture hides?

**Background.** A tilted rotating disk has a characteristic velocity map, the "spider diagram": gas on one side approaches, gas on the other recedes. `v_obs = v_sys + v_rot(R) sin i cos θ`, where θ is the angle around the disk. A **tilted-ring model** fits rings of centre, systemic velocity, position angle, inclination and rotation speed, and so turns the map into a rotation curve. Radio 21-cm maps reach far beyond the stars; that is how flat rotation curves became undeniable (Bosma 1978). Optical integral-field spectroscopy maps stars and ionized gas in the inner parts. Stars rotate slower than gas because their random motions contribute part of the support against gravity.

**Data.** THINGS (Walter et al. 2008): VLA 21-cm data for 34 nearby galaxies, including NGC 3198, NGC 2403, NGC 2903 and DDO 154, with public cubes and velocity-field maps (`mpia.de/THINGS`). MaNGA DR17 (10,010 galaxies): stellar and Hα velocity maps from the Data Analysis Pipeline (MAPS files or the `marvin` package).

**Method + visualizations.**
- NGC 3198's velocity field (the spider diagram) and our own tilted-ring fit (scipy least squares): first with fixed geometry, then free per ring (outer warps appear) → rotation curve, approaching and receding sides separately.
- The residual velocity field: spiral-arm streaming, and NGC 2903's bar.
- Our curves over de Blok et al. (2008) and SPARC (D4). DDO 154's inner slope, the cusp–core problem, measured directly.
- MaNGA: stellar vs gas rotation for a few spirals (asymmetric drift), and a search for galaxies whose gas rotates *against* their stars.

**Check against.** de Blok et al. (2008) THINGS rotation curves; D4's SPARC curves for the same galaxies.

**Watch out for.** Simple intensity-weighted velocity maps are biased where line profiles are skewed; de Blok et al. used Gauss–Hermite fits to the cubes. Beam smearing flattens inner gradients, which matters most for dwarfs and the cusp–core question. v_rot ∝ 1/sin i, so nearly face-on galaxies are unreliable.

**Builds on.** D4, H2.

## J2. Galaxy morphology with Galaxy Zoo ★★☆

**Science question.** Hubble sorted galaxies into spirals, ellipticals and lenticulars in 1926. With 16 million volunteer classifications of 300,000 galaxies, how do shape, colour, mass and environment relate?

**Background.** Shape correlates with star formation. Spirals are mostly blue and forming stars; ellipticals are mostly red and quiet. The exceptions are informative: red spirals (Masters et al. 2010) and blue ellipticals (Schawinski et al. 2009). Ellipticals crowd into cluster cores: the **morphology–density relation** (Dressler 1980). Bars are more common in redder disk galaxies (Masters et al. 2011).

**Data.** Galaxy Zoo 2 (Willett et al. 2013): 304,122 SDSS galaxies with debiased vote fractions for every question (`data.galaxyzoo.org`). Galaxy Zoo DESI (Walmsley et al. 2023): machine-learning morphologies for 8.7 million galaxies. SDSS photometry and MPA-JHU stellar masses; D2's local-density estimates; image cutouts from the Legacy Survey (as in G1).

**Method + visualizations.**
- Clean samples (debiased fraction > 0.8): the colour–magnitude diagram split by shape, and how many galaxies are the "wrong" colour.
- Bar fraction vs colour and mass; number of spiral arms vs mass.
- Morphology–density with D2's density measure: the elliptical fraction climbing into clusters. Separate this from the colour–density trend D2 already found.
- A tuning-fork poster of real cutouts, chosen by vote fraction.
- ML extension: train a small CNN on cutouts using Galaxy Zoo labels and compare it with the published Zoobot predictions.

**Check against.** Willett et al. (2013) class fractions; Masters et al. (2010, 2011).

**Watch out for.** Distant galaxies look smoother to volunteers, so use debiased fractions and a redshift-limited sample. A vote fraction is not a probability. Keep training and test galaxies separate in the ML extension.

**Builds on.** D2 (density, colour), G1 (images); pairs with J3.

## J3. The galaxy life cycle ★★☆

**Science question.** Galaxies are either busily forming stars or have mostly stopped. How do they cross from one group to the other, and when did the universe as a whole form most of its stars?

**Background.** Star-forming galaxies lie on a tight **main sequence**: star-formation rate rises almost in proportion to stellar mass, with ~0.3 dex scatter. "Quenched" galaxies sit far below it, and the sparse region in between is the **green valley**. Quenching depends on mass (probably black hole feedback, I1) and on environment (D2's red fraction), and the two effects appear separable (Peng et al. 2010). Across cosmic time, the universe's star-formation rate peaked at z ≈ 1.9, ~3.5 Gyr after the Big Bang, and has declined since with an e-folding time of 3.9 Gyr (Madau & Dickinson 2014).

**Data.** The MPA-JHU value-added catalogue for SDSS DR8: star-formation rates (Brinchmann et al. 2004) and stellar masses (SkyServer tables `galSpecExtra`, `galSpecInfo`). GSWLC (Salim et al. 2016, 2018): UV-to-IR SED-fit masses and star-formation rates for ~700,000 galaxies. The Madau & Dickinson (2014) compilation of star-formation-rate density measurements (their Table 1).

**Method + visualizations.**
- The SFR–mass plane (density map): the main sequence fitted, the green valley and the quenched cloud.
- Quenched fraction vs mass and vs D2's local density → mass and environment quenching as separate effects.
- B7's BPT classes across the plane: AGN crowd the green valley.
- **Cosmic star-formation history**: the Madau & Dickinson points with their fit `ψ(z) = 0.015 (1+z)^2.7 / [1 + ((1+z)/2.9)^5.6]` M☉ yr⁻¹ Mpc⁻³. Integrate for the build-up of stellar mass; convert z to cosmic time with D3's cosmology.
- Our own z ≈ 0.1 anchor point from SDSS (total SFR per volume in a volume-limited sample) against their z = 0 value.

**Check against.** Peak at z ≈ 1.9; half of today's stars formed before z ≈ 1.3 (Madau & Dickinson 2014).

**Watch out for.** Star-formation indicators (Hα, UV + IR, SED fitting) disagree at the ~0.2–0.3 dex level. SDSS fibres see only galaxy centres (MPA-JHU applies aperture corrections). Assumptions about the stellar mass spectrum (the "IMF") shift everything: Madau & Dickinson assume Salpeter, while most catalogues use Chabrier (~0.2 dex offset).

**Builds on.** D2, B7, D3; pairs with I3 (black hole growth peaks at the same epoch) and J8 (the high-z end).

## J4. Laniakea: the flows of galaxies ★★☆

**Science question.** Galaxies don't only recede with the Hubble flow; they also fall toward concentrations of mass. Can we map those extra motions, find the Great Attractor, and explain why we move at ~370 km/s relative to the CMB?

**Background.** A galaxy's **peculiar velocity** is `v_pec = cz − H0 d`. It needs a distance measured without redshift: Tully–Fisher for spirals (rotation speed ↔ luminosity, the relation D4 reproduced), the Fundamental Plane for ellipticals, supernovae, Cepheids or the tip of the red giant branch. Each distance is uncertain by ~20%, so maps average many galaxies. The CMB dipole says the Sun moves at 369.8 km/s toward (l, b) = (264°, 48°), and the Local Group at ~620 km/s. Tully et al. (2014) defined our home supercluster, **Laniakea**, as the basin of galaxies flowing toward the Great Attractor region.

**Data.** Cosmicflows-4 (Tully et al. 2023, ApJ 944, 94): distances to 55,877 galaxies in 38,065 groups, on VizieR (`J/ApJ/944/94`) and the Extragalactic Distance Database (`edd.ifa.hawaii.edu`). The 2MASS Redshift Survey (Huchra et al. 2012; 44,599 galaxies over the whole sky) for the density field. The Planck CMB dipole (D3).

**Method + visualizations.**
- Cosmicflows-4 Hubble diagram → H0, then note that it is set by the Cepheid/TRGB zero point, i.e. B4's ladder.
- Peculiar velocities with an unbiased estimator (Watkins & Feldman 2015). Sky maps in supergalactic coordinates coloured by infall and outflow: toward Virgo, the Great Attractor and Shapley, and away from the Local Void.
- Bulk flow in spheres of growing radius vs the ΛCDM expectation.
- The direction of the gravitational pull (the dipole of 2MRS galaxies) vs the CMB dipole's direction.

**Check against.** Tully et al. (2023) H0 ≈ 75 km/s/Mpc; the CMB dipole direction and amplitude.

**Watch out for.** Malmquist bias: distance errors scatter galaxies in and out of volume, so work in log-distance ratios. The Milky Way's disk hides part of the sky (the "zone of avoidance"), including part of the Great Attractor. Cosmicflows-4 mixes methods with different zero points. H0 and the zero point are degenerate.

**Builds on.** D1 (H0), B4 (zero point), D3 (CMB dipole), D4 (Tully–Fisher), D2 (structure).

## J5. Weighing galaxies with Einstein rings ★★☆

**Science question.** When a massive galaxy sits almost exactly in front of a more distant one, gravity bends the background light into arcs or a ring. The ring's size weighs the foreground galaxy directly. How much of that mass is dark?

**Background.** The **Einstein radius** θ_E measures the mass inside it: `θ_E = √(4GM(<θ_E)/c² · D_ls/(D_l D_s))`. For a "singular isothermal sphere", the simplest model of a galaxy's total mass, `θ_E = 4π (σ/c)² D_ls/D_s`. The Sloan Lens ACS Survey (SLACS) found lenses by spotting two redshifts in one SDSS spectrum (a background galaxy's emission lines inside an elliptical's spectrum) and imaged them with Hubble. It confirmed 85 grade-A lenses.

**Data.** SLACS tables (Bolton et al. 2008; Auger et al. 2009): lens and source redshifts, Einstein radii, velocity dispersions, effective radii and stellar masses. HST ACS images of the lenses from MAST for display. The SDSS spectra themselves (B7).

**Method + visualizations.**
- Lensing σ (from θ_E and redshifts, using D3's cosmology) vs the stars' measured σ. The ratio is ≈1: ellipticals' total mass behaves as isothermal even though neither stars nor dark matter alone does.
- Projected mass inside θ_E vs stellar mass → the dark fraction, and how it depends on the IMF assumed for the stars.
- A gallery of HST rings with simple lens-model overlays, and one SDSS spectrum showing the background galaxy's [O II] line.
- Extension: model one lens ourselves (singular isothermal ellipsoid + background source, `lenstronomy`).

**Check against.** Bolton et al. (2008) σ ratio near unity; Auger et al. (2010) and Treu et al. (2010) dark fractions and IMF preference.

**Watch out for.** Velocity dispersions need aperture corrections (3″ fibre). The stellar mass, and so the dark fraction, depends on the IMF. The isothermal model is an idealization. The lens search favours massive, high-σ galaxies.

**Builds on.** D4 (dark matter by other means), B7 (spectra), G1 (HST images), I1 (σ), D3 (distances).

## J6. Dwarf galaxies of the Local Group ★☆☆

**Science question.** The faintest galaxies are the most dark-matter-dominated objects known. What does the census of the Milky Way's and Andromeda's small companions tell us about dark matter?

**Background.** In a dwarf spheroidal galaxy, stars move randomly at only 3–10 km/s. The mass inside the half-light radius follows from `M_½ ≈ 4σ²R_e/G` (Wolf et al. 2010). For the faintest "ultra-faint" dwarfs this gives mass-to-light ratios in the hundreds to thousands. Dozens of faint satellites found by SDSS, DES and others have eased the "missing satellites" problem. Some satellites orbit in thin planes (around the Milky Way and M31, Ibata et al. 2013), a puzzle for simple models. Several arrived as companions of the Large Magellanic Cloud.

**Data.** McConnachie (2012) catalogue of Local Group and nearby dwarf galaxies, updated online version with Gaia proper motions (`NearbyGalaxies_Jan2021_PUBLIC.fits`, at the CADC "Local Group and Nearby Dwarf Galaxies" page). The Local Volume Database (Pace 2024; GitHub `apace7/local_volume_database`). Gaia EDR3 systemic proper motions (Pace et al. 2022).

**Method + visualizations.**
- Luminosity vs size for galaxies and globular clusters: where does a tiny galaxy end and a star cluster begin? (Dark matter and a spread in metallicity are the tests.)
- M_½ vs luminosity: mass-to-light climbing to thousands at the faint end.
- 3D map of the satellites around the Milky Way and M31, with the planes.
- Orbits from Gaia proper motions in a Milky Way potential (`galpy`): which satellites are on first infall, and which came with the LMC.

**Check against.** McConnachie (2012) and Wolf et al. (2010) mass estimates; Strigari et al. (2008), ~10⁷ M☉ within 300 pc for all dwarfs.

**Watch out for.** Binary stars inflate the tiny velocity dispersions of ultra-faint dwarfs. Many have only tens of member stars with velocities. Some dwarfs are being tidally torn apart, so equilibrium masses are suspect there. Faint satellites are found only where deep surveys have looked.

**Builds on.** D4 (dark matter), H3 (dwarfs that were already destroyed); feeds J7.

## J7. Andromeda is coming — or is it? (labeled simulation) ★★☆

**Science question.** Andromeda approaches us at ~110 km/s. Textbooks say the two galaxies will merge in ~4.5 Gyr, but a 2025 study using Hubble and Gaia put the odds at about a coin flip within 10 Gyr. Can we reproduce the odds from the measured motions?

**Background.** In 1959 Kahn & Woltjer argued that the Milky Way and M31 must once have been flying apart and have since turned around. Treating them as a two-body orbit that started at the Big Bang gives their total mass. The result, several × 10¹² M☉, was early evidence for dark matter from just two numbers. The collision question hinges on Andromeda's *sideways* motion, which Hubble (van der Marel et al. 2012) and Gaia (van der Marel et al. 2019; Salomon et al. 2021) measured with large uncertainties. Sawala et al. (2025, *Nature Astronomy*) included M33 and the LMC and sampled the uncertainties. They found a ~50% chance of a merger within 10 Gyr and <2% for a head-on hit; the LMC merges with the Milky Way first.

**Data.** Measured inputs only: M31's distance and radial velocity; its HST and Gaia proper motions; M33's and the LMC's positions and velocities (J6's catalogue); published mass estimates with uncertainties.

**Method + visualizations.**
- The timing argument first: a radial Kepler orbit with an age of 13.8 Gyr (D3) → Local Group mass.
- Monte Carlo: draw from every measured uncertainty → integrate MW, M31, M33 and LMC as extended masses with dynamical friction for 10 Gyr → merger probability vs time. Then drop the LMC or M33, or use HST-only vs Gaia proper motions, and see what moves the answer.
- Animated orbits; probability-vs-time curves.
- Optional: a restricted N-body run (test-particle disks, Toomre & Toomre 1972 style) for what a merger would look like.
- **Labeled as a simulation per ground rule 1**, with every initial condition from data.

**Check against.** Sawala et al. (2025): ~50% within 10 Gyr; the timing-argument mass (e.g. Li & White 2008, ~5 × 10¹² M☉).

**Watch out for.** Galaxy masses dominate the uncertainty. Dynamical friction is approximate. The answer is a probability distribution, not a date.

**Builds on.** J6, H2 (Milky Way mass), E1/E2 (orbit integration), D3 (cosmic age).

## J8. The first galaxies with JWST ★★★

**Science question.** JWST now sees galaxies from when the universe was under 300 million years old. How far back does the record go, how do we know their distances, and are there "too many, too bright" early galaxies?

**Background.** Hydrogen absorbs light blueward of the Lyman-α line (121.6 nm). At z ~ 10–14 this **Lyman break** shifts into the near-infrared, so a galaxy vanishes from JWST's bluer filters: a "dropout". Photometric candidates need spectra. A famous z ~ 16 candidate turned out to be a dusty galaxy at z = 4.9. Spectroscopic records: JADES-GS-z14-0 at z ≈ 14.2–14.3 (2024) and MoM-z14 at z = 14.44 (2025). There are more bright galaxies at z > 10 than pre-JWST models predicted. The candidate explanations (bursty star formation, little dust, efficient early star formation) do not require overturning ΛCDM.

**Data.** JADES (the JWST Advanced Deep Extragalactic Survey) on MAST as High Level Science Products (`archive.stsci.edu/hlsp/jades`): NIRCam mosaics and photometric catalogues (DR5) and NIRSpec spectra with redshifts (DR3 onward).

**Method + visualizations.**
- Colour–colour dropout selection at z ~ 9–15 from the catalogue, compared with JADES's photometric redshifts.
- Spectra of the record holders: the Lyman break, plus GN-z11's surprising nitrogen lines.
- Redshift → cosmic time with D3's cosmology.
- A rough UV luminosity function at z ~ 10–14, drawn over published ones and pre-JWST predictions.
- Colour images of the deep field (G1's methods).

**Check against.** JADES spectroscopic redshifts; published z > 10 UV luminosity functions (e.g. Donnan et al. 2024; Finkelstein et al. 2024).

**Watch out for.** Photometric redshifts are contaminated by dusty z ~ 4–5 galaxies and brown dwarfs. The effective survey volume (completeness) is the hard part, so show the published luminosity functions as the primary result. Deep fields are small, so cosmic variance is large.

**Builds on.** G1 (imaging), B7 (spectra), J3 (the start of cosmic star formation), D3 (cosmology).

## J9. Fly into the Hubble Ultra Deep Field ★★☆

**Science question.** The Hubble Ultra Deep Field shows thousands of galaxies in a patch of sky about the size of a grain of sand held at arm's length. They are spread across most of cosmic history, but the picture flattens them onto one plane. Can we use measured redshifts to put each galaxy back at its true distance and fly into the field in 3D?

**Background.** Hubble's eXtreme Deep Field (XDF; Illingworth et al. 2013) combines about two million seconds of exposure on this patch. MUSE on ESO's Very Large Telescope then took a spectrum at every pixel. Its second data release measured secure redshifts for 2,221 sources (Bacon et al. 2023): 25 nearby galaxies, 677 [O II] emitters at z = 0.25–1.5, 201 galaxies at z = 1.5–2.8, and 1,308 Lyman-α emitters at z = 2.8–6.7. Many of the Lyman-α emitters are invisible in Hubble's images. Redshift turns into distance through D3's cosmology, giving a 3D cone that reaches back to ~800 million years after the Big Bang.

**Data.** The XDF images (MAST High Level Science Product `xdf`: ACS and WFC3 mosaics at 30 or 60 mas). The MUSE HUDF DR2 catalogue (Bacon et al. 2023; VizieR `J/A+A/670/A4`, full products at AMUSED, `amused.univ-lyon1.fr`). Optionally JADES photometric redshifts (J8) for galaxies without a MUSE redshift.

**Method + visualizations.**
- Cut a stamp of every galaxy with a MUSE redshift from an XDF colour image (G1's methods).
- Place each stamp at its comoving distance in a 3D scene (plotly or a WebGL viewer), scaled to its true physical size.
- The keeper: a fly-through movie from the Milky Way out to z ≈ 6.7, with a clock showing lookback time.
- The redshift histogram, with MUSE's harder range (z = 1.5–2.8, where few strong emission lines fall in its band) and the Lyman-α emitters invisible to Hubble marked.
- Galaxy sizes and colours against cosmic time.

**Check against.** Bacon et al. 2023's redshift counts; STScI's published visualizations of the field.

**Watch out for.** Comoving distance and light-travel distance differ: pick one and label it. The MUSE sample is selected by emission lines, so it is not a complete census. At different redshifts the same filter shows different rest-frame light.

**Builds on.** G1 (colour images), J8 (the JWST deep field), D3 (distances), J3 (galaxy evolution).

---

# K. Dark matter

D4 found dark matter in rotation curves, J5 in gravitational lensing, J6 in dwarf galaxies and H3 in the gaps of a stellar stream. All of those detect its gravity. This section asks what dark matter is *not* (atoms, black holes), weighs it where most of it sits (galaxy clusters and our own halo), and follows the experiments trying to catch a particle of it.

**Suggested dark matter arc:** K1 (not atoms) → K4 (not black holes) → K2 + K3 (clusters) → K5 (our own halo) → K6 (catching a particle).

## K1. Why dark matter cannot be atoms: Big Bang nucleosynthesis ★★☆

**Science question.** In its first ~20 minutes the universe fused hydrogen into deuterium, helium and a trace of lithium. How much of each formed depends mainly on one number, the density of ordinary matter. Can we measure that density from the lightest elements, compare it with the CMB's value from 380,000 years later, and show it is only about a sixth of all the matter there is?

**Background.** At about one second, neutrons and protons stop converting into each other, at a ratio n/p ≈ 1/6. Free neutrons then decay (lifetime ~880 s) until deuterium can survive the bath of photons, after ~3 minutes. By then n/p ≈ 1/7. Almost every remaining neutron ends up in helium-4, so the helium mass fraction is `Y_p ≈ 2(n/p)/(1 + n/p) ≈ 0.25`. Deuterium is the leftover that did not burn: the denser the baryons, the less survives. That makes D/H a sensitive **baryometer**, read as η (baryons per photon, ~6 × 10⁻¹⁰) or as Ω_b h². The measurements:
- D/H = (2.527 ± 0.030) × 10⁻⁵ in near-pristine gas clouds seen against quasars (Cooke, Pettini & Steidel 2018);
- Y_p = 0.2458 ± 0.0013 from metal-poor star-forming galaxies (LBT Y_p Project, 2026);
- lithium in old halo stars, ~3 times below the prediction (the "lithium problem").

The CMB's acoustic peaks give Ω_b h² independently: D3 found 0.02221 ± 0.00019. All matter together has Ωm h² ≈ 0.14, so baryons are ~16% of the matter. The rest is not made of atoms.

**Data.** Small published tables: Cooke et al. 2018 (seven D/H systems); Aver et al. 2021 and the LBT Y_p Project (2026) for helium; Spite-plateau lithium. BBN codes with current nuclear rates: PRyMordial (GitHub `vallima/PRyMordial`) and LINX (GitHub `cgiovanetti/LINX`, JAX-based); neither is on PyPI. D3's Ω_b h² and Ω_c h².

**Method + visualizations.**
- A toy calculation first: weak-interaction freeze-out, neutron decay and the deuterium bottleneck → Y_p(η). It should land within a few percent of the full codes.
- The full network (PRyMordial or LINX): D, ³He, ⁴He and ⁷Li against η, the classic "Schramm plot", with the measured bands across it.
- η from deuterium alone → Ω_b h², beside D3's CMB value. Two epochs, 3 minutes and 380,000 years after the Big Bang, give one answer.
- Ω_b / Ωm → the fraction of matter that is not baryons. Helium as a count of neutrino species (N_eff).
- The lithium problem, shown honestly.

**Check against.** Cooke et al. 2018: BBN Ω_b h² within 2σ of Planck. Planck 2018: Ω_b h² = 0.02237 ± 0.00015. LUNA's measured deuterium-burning rate (Mossa et al. 2020): Ω_b h² = 0.02233 ± 0.00036.

**Watch out for.** The d(p,γ)³He reaction rate shifts the predicted D/H by about as much as the measurement error. Codes with different rate compilations disagree at ~1σ, so say which rates we use. Helium abundances need corrections for temperature, density and ionization; the EMPRESS survey's lower value (~0.237) shows how much those matter. Lithium may have been destroyed inside the stars rather than missing from the Big Bang.

**Builds on.** D3 (Ω_b h² from the acoustic peaks), C3 (the FRB census found the baryons where BBN says they should be), B7 (spectra), C7 (neutrinos).

## K2. Weigh a galaxy cluster three ways ★★☆

**Science question.** In 1933 Zwicky found that the galaxies of the Coma cluster move far too fast to be held by their visible mass: the first evidence for dark matter. Today we can weigh clusters by their galaxies' motions, by the hot gas that fills them, and by lensing. Do the methods agree? And does the share of a cluster's mass in ordinary matter reveal the matter density of the whole universe?

**Background.**
- **Galaxy motions.** The virial theorem gives M ~ 3σ²R/G for line-of-sight velocity dispersion σ and cluster radius R.
- **Hot gas.** Gravity holds in the X-ray-emitting gas (10⁷–10⁸ K). In hydrostatic equilibrium, `M(<r) = −(k T r / G μ m_p) (d ln n / d ln r + d ln T / d ln r)`.
- **The baryon budget.** Clusters are big enough to hold a fair sample of the universe, so their (gas + stars) / total mass should equal Ω_b/Ωm. With Ω_b from K1 this gives Ωm ≈ 0.3, not the then-favoured 1 (White et al. 1993, "a challenge to cosmological orthodoxy").

**Data.** X-COP (Eckert et al. 2019; Ettori et al. 2019; Ghirardini et al. 2019): XMM-Newton + Planck SZ profiles of gas density, temperature, pressure, hydrostatic mass, gas fraction and stellar mass for 12 massive nearby clusters, out to R₂₀₀. The full release is one 315 MB archive on SWITCHdrive, linked from the X-COP data page (`dominiqueeckert.wixsite.com/xcop/data`); the older per-file links on astro.unige.ch now redirect elsewhere. D2's SDSS spectroscopic galaxies for velocity dispersions (D2 measured Coma's: 963 ± 28 km/s). Planck's PSZ2 catalogue for SZ-based masses.

**Method + visualizations.**
- Zwicky's argument redone with D2's Coma dispersion and today's H0 → a mass-to-light ratio in the hundreds.
- Virial masses for the X-COP clusters inside the SDSS footprint (members chosen by redshift with iterative clipping), against the X-COP hydrostatic and Planck SZ masses: the three-way mass comparison.
- Gas fraction profiles for all 12 clusters, plus stars → the baryon fraction at R₅₀₀ and R₂₀₀, rising toward the cosmic Ω_b/Ωm (K1, D3). Then Ωm = Ω_b / f_b.
- The hydrostatic mass bias: if the baryon fraction must equal the cosmic value, how much non-thermal pressure is allowed? Compare with the (1 − b) ≈ 0.6 that Planck's cluster counts needed (D11).

**Check against.** Eckert et al. 2019: non-thermal pressure ~6% at R₅₀₀ and ~10% at R₂₀₀; with Planck's mass calibration, 1 − b = 0.85 ± 0.05. Ettori et al. 2019: hydrostatic masses to 3% (statistical) at R₅₀₀.

**Watch out for.** Interlopers and substructure inflate virial masses, and mass scales as σ², so a 10% error in σ is a 20% error in mass. Hydrostatic masses miss the pressure from turbulence and bulk motions. Stars and intracluster light are a small but uncertain term. Gas fractions rise with radius, so compare clusters at the same overdensity.

**Builds on.** D2 (SDSS, Coma), D3 and K1 (Ω_b/Ωm), D11 (cluster mass bias), J5 (lensing masses).

## K3. The Bullet Cluster: dark matter that passed straight through ★★☆

**Science question.** About 150 million years ago two galaxy clusters collided at ~4,700 km/s. Their hot gas, which holds most of their ordinary matter, collided and stalled in the middle. Their galaxies and most of their mass, as mapped by lensing, passed straight through. Can we rebuild the image often called the most direct evidence for dark matter, and use it to limit whether dark matter particles collide with each other?

**Background.** In a relaxed cluster, gas and dark matter sit together, so lensing and X-rays cannot tell extra gravity apart from unseen mass. In a collision they separate. Clowe et al. (2006) found the lensing mass peaks on the galaxies, 8σ away from the X-ray gas that holds most of the baryons. Modified gravity without dark matter would put the mass peaks on the gas. The smaller cluster's gas, the "bullet", drives a bow shock with a Mach number of ~3. If dark matter particles scattered off one another, the dark matter would lag behind the galaxies. The absence of a lag limits the self-interaction cross-section to σ/m ≲ 1 cm²/g (Markevitch et al. 2004; Randall et al. 2008). Combining 72 cluster collisions gives σ/m < 0.47 cm²/g (Harvey et al. 2015).

**Data.** The Clowe et al. 2006 public release: the weak-lensing convergence (κ) map and the gas surface-mass map, as FITS files. The project page is gone, but the Internet Archive serves the files (`flamingos.astro.ufl.edu/1e0657/data/1e0657.release1.kappa.fits` and `1e0657_central_gasmass_mod.fits`). The ~500 ks of Chandra ACIS-I observations of 1E 0657−56 (Chandra Data Archive), binned from the event files with `astropy`. Optical images from the Legacy Survey or HST, as in G1.

**Method + visualizations.**
- The famous overlay, built from data: optical image, exposure-corrected X-ray gas and κ contours.
- Positions of the κ, gas and galaxy-light peaks, with bootstrap errors, and their offsets.
- A mass budget in apertures around each peak: lensing mass against gas mass against stellar mass. The gas dominates the baryons, yet the mass follows the galaxies.
- The bow shock: X-ray surface brightness across the front → density jump → Mach number (Rankine–Hugoniot).
- A simple self-interaction bound: the bullet's dark matter survived the passage, so its scattering optical depth (σ/m) × Σ must be ≲ 1.

**Check against.** Clowe et al. 2006: 8σ offset. Markevitch et al. 2004 and Randall et al. 2008: σ/m ≲ 1 cm²/g. Markevitch 2006: Mach ≈ 3.

**Watch out for.** κ maps are noisy, with edge effects, so compare peak significances rather than pixel values. The gas-mass map depends on a model of the gas. A single system does not rule out every alternative theory of gravity (some add massive neutrinos to cope), so state exactly what it does rule out. We see a 3D collision in projection.

**Builds on.** J5 (lensing), K2 (cluster masses), G1 (images), D2 (clusters).

## K4. Are black holes the dark matter? ★★☆

**Science question.** If the Milky Way's dark halo were made of compact objects (dead stars, planets, or black holes born in the Big Bang), some would pass in front of stars in the Large Magellanic Cloud and briefly magnify them. Twenty years of OGLE monitoring of 78.7 million LMC stars turned up 16 such events. What does that rule out, including the idea that LIGO's black holes are the dark matter?

**Background.** The probability that a given star is being lensed at any moment is the **optical depth** τ. For a halo made entirely of compact objects it is ~5 × 10⁻⁷ toward the LMC, whatever their mass (Paczyński 1986). The mass sets how long events last (t_E ∝ √M), and 10–100 M☉ objects would cause events lasting years, which earlier surveys could not detect. The MACHO project (Alcock et al. 2000) saw 13–17 events and suggested ~20% of the halo in ~0.5 M☉ objects. EROS-2 and OGLE found far fewer, consistent with ordinary stars in the Milky Way disk and in the LMC itself. After LIGO's first detections, ~30 M☉ primordial black holes were proposed as the dark matter (Bird et al. 2016). OGLE's 20-year search found no event longer than a year: 10 M☉ black holes can be at most 1.2% of the dark matter, and 100 M☉ ones at most 3% (Mróz et al. 2024, *Nature*).

**Data.** OGLE's LMC optical-depth release (Mróz et al. 2024, *ApJS* and *Nature*): event light curves, detection efficiencies against timescale for each field, source-star counts, posterior samples and the published limits. It is at `www.astrouw.edu.pl/ogle/ogle4/LMC_OPTICAL_DEPTH/` and on Zenodo (DOI 10.5281/zenodo.10879576). OGLE EWS bulge events (I12) for practice. The MACHO and EROS-2 published results for comparison.

**Method + visualizations.**
- Paczyński fits to the 13 events found by the pipeline (and the 3 found by eye) → their timescales.
- A model halo (cored isothermal or NFW) of compact objects of a single mass M, plus lensing by Milky Way disk stars and by the LMC's own stars → the expected event rate against t_E for each field.
- Multiply by OGLE's efficiencies and observing time → the expected number of events for a halo fraction f. Compare with what was seen using Poisson statistics → an upper limit on f at each M.
- The keeper: f_max against M from 10⁻⁵ to 10³ M☉, with MACHO's claim, EROS-2's limits and the LIGO-motivated window marked.

**Check against.** Mróz et al. 2024: τ = (0.121 ± 0.037) × 10⁻⁷; f < 1% for 1.8 × 10⁻⁴ to 6.3 M☉; at most 1.2% at 10 M☉, 3.0% at 100 M☉ and 11% at 1,000 M☉.

**Watch out for.** The halo model (density profile, velocity distribution, the LMC's distance and motion) moves the limits by tens of percent. A 2025 critique ("Eppur non si trovano", arXiv:2510.20005) argues exactly this; read it alongside. Self-lensing by LMC stars must be modelled, not ignored. A spread of black hole masses changes the limits. All events are fitted as single point lenses.

**Builds on.** I12 (microlensing fits), C5 (LIGO's black holes), B4 (the LMC's distance), K5 (the halo model).

## K5. Weigh the Milky Way's dark halo ★★★

**Science question.** How much dark matter does our galaxy hold, and how much is passing through the Earth right now? Estimates of the Milky Way's total mass span a factor of five, and a 2023 claim of a sharply falling rotation curve would make the halo surprisingly light. Can we weigh it several ways with Gaia and see which answers agree?

**Background.** Four routes:
- **The rotation curve** to ~25 kpc, from Jeans modelling of red giants (Eilers et al. 2019: 229 km/s at the Sun, falling slowly). Some Gaia DR3 analyses find a near-Keplerian drop beyond ~19 kpc and a total mass of only ~2 × 10¹¹ M☉ (Jiao et al. 2023). Others blame disequilibrium and a lopsided disk.
- **The escape speed.** The fastest stars still bound near the Sun mark v_esc ≈ 500–550 km/s, which gives M₂₀₀ ≈ 10¹² M☉ (Deason et al. 2019).
- **Distant tracers.** Globular clusters and satellite galaxies with Gaia proper motions, and their orbits.
- **The local density.** The up-and-down motions of disk stars measure the total density near the Sun. Subtracting stars and gas leaves ρ_DM ≈ 0.01 M☉/pc³ ≈ 0.4 GeV/cm³, the number every dark matter detector assumes (K6).

**Data.** Gaia DR3's 33.8 million stars with radial velocities (H3's queries) and Bailer-Jones distances. APOGEE DR17 + Gaia giants for the outer curve, and published curves (Eilers et al. 2019; Jiao et al. 2023; Ou et al. 2024). Globular cluster mean motions (Vasiliev & Baumgardt 2021, on Baumgardt's catalogue pages), J6's satellites and H2's 21-cm curve. A 2026 review of the debate (arXiv:2608.10189).

**Method + visualizations.**
- A rotation curve from Jeans modelling out to 25 kpc, with an error budget, over H2's and the published curves. Where does it fall, and by how much?
- The escape speed: fit the high-velocity tail of halo stars with f(v) ∝ (v_esc − v)^k (Leonard & Tremaine 1990) → v_esc → halo mass.
- Joint halo fits (NFW, Einasto, cored) to the rotation curve, the escape speed and the cluster and satellite motions → M₂₀₀ and concentration. Can a light halo fit everything?
- The local density: a vertical Jeans analysis of stars within ~1 kpc of the plane → ρ_DM in M☉/pc³ and GeV/cm³.

**Check against.** Eilers et al. 2019: 229.0 ± 0.2 km/s at R₀, slope −1.7 km/s/kpc. Deason et al. 2019: M₂₀₀ = 1.00 (+0.31/−0.24) × 10¹² M☉. Jiao et al. 2023: ~2 × 10¹¹ M☉. Local density ~0.3–0.6 GeV/cm³ (de Salas & Widmark 2021 review).

**Watch out for.** Jeans models assume a steady, axisymmetric galaxy, but the disk is warped and still ringing (H3's phase spiral). Distances to faint giants beyond 20 kpc are uncertain. The escape-speed fit depends on the assumed shape of the tail and on unbound stars (hypervelocity stars, LMC debris). The LMC holds ~10–20% of the Milky Way's mass and drags the outer halo; ignoring it biases the outer tracers.

**Builds on.** H2 (inner rotation curve), H3 (Gaia 6D, potentials), J6 and J7 (satellites, Local Group mass), D4 (other galaxies' curves), K6 (local density).

## K6. Hunting dark matter particles underground ★★☆

**Science question.** If dark matter is a sea of particles, the Earth moves through it at over 200 km/s and should occasionally knock a nucleus in a detector. Tonne-scale xenon experiments have seen nothing. One sodium-iodide experiment, DAMA/LIBRA, has reported a yearly signal for over two decades. Can we compute what a dark matter particle would do in a detector, rebuild a published limit from public data, and test DAMA's claim?

**Background.** The expected recoil rate depends on the local dark matter density (K5), its speed distribution, the target nucleus, and the particle's mass and cross-section. The "Standard Halo Model" is a Maxwellian with v₀ ≈ 238 km/s, cut off at the escape speed of 544 km/s, with ρ = 0.3 GeV/cm³ (Baxter et al. 2021). Recoils carry a few to tens of keV. Earth's orbit adds to the Sun's motion through the halo in June and subtracts in December, so the rate should peak around 2 June (Drukier, Freese & Spergel 1986). DAMA reports such a modulation: (0.01014 ± 0.00074) counts/day/kg/keV in 2–6 keV, at 13.7σ. Experiments with the same crystals see none: COSINE-100 and ANAIS-112 together exclude DAMA's signal at 4.7σ (1–6 keV, 2025). LZ, with 4.2 tonne-years, limits the WIMP–nucleon cross-section to below 2.2 × 10⁻⁴⁸ cm² at 40 GeV. Xenon detectors now see solar ⁸B neutrinos bumping nuclei, the start of the "neutrino fog" that will eventually limit such searches.

**Data.** LZ's 2024 results on HEPData (record 155182: limits, sensitivity bands, efficiencies and analysis data points). XENONnT and PandaX-4T published limits. DAMA's and COSINE-100/ANAIS-112's published modulation amplitudes and phases; residual-rate plots digitized only where needed, and labelled as such.

**Method + visualizations.**
- Recoil spectra dR/dE for xenon, argon, germanium and sodium iodide at 10, 50 and 500 GeV (Lewin & Smith 1996 formulae, Helm form factor).
- Annual modulation: the rate through the year, its amplitude (a few percent) and its phase.
- An approximate 90% limit from LZ's efficiency, exposure and observed events (Poisson, or Yellin's maximum gap), against the official curve, with the neutrino fog drawn on.
- DAMA against COSINE-100 and ANAIS-112: modulation amplitudes in 1–6 and 2–6 keV, with error bars, on one plot.
- The WIMP mass–cross-section plane, with limits over the decades.
- Extension (indirect detection): stack Fermi-LAT gamma-ray photons on J6's dwarf galaxies.

**Check against.** LZ 2024: 2.2 × 10⁻⁴⁸ cm² at 40 GeV (90% CL). XENONnT 2025: 1.7 × 10⁻⁴⁷ cm² at 30 GeV. COSINE-100 + ANAIS-112: 0.0005 ± 0.0019 counts/day/kg/keV in 1–6 keV.

**Watch out for.** Our simple limit will differ from the published profile-likelihood limit (backgrounds, signal acceptance), so quote the ratio. The halo parameters (ρ, v₀, v_esc) shift limits by tens of percent, more at low mass. Sodium-iodide experiments differ in energy calibration (quenching factors), so compare like with like. No WIMP signal does not mean no dark matter: axions and other candidates need different experiments.

**Builds on.** K5 (local density and speeds), J6 (dwarfs for gamma rays), C7 (neutrino detectors).

---

# L. Supernovae: how stars explode

D1 and D6 use one kind of supernova as a cosmic yardstick, and C7 caught the neutrinos of SN 1987A. This section asks what supernovae *are*. A star can blow up in two ways.

**Core collapse** (Types II, Ib and Ic). A star born with more than ~8 M☉ builds an iron core inside shells of lighter elements (B11). Fusing iron releases no energy, so nothing holds the core up for long. When it reaches ~1.4 M☉ it collapses in a fraction of a second, from a few thousand kilometres across to a ~30 km ball of nuclear-density matter. The core bounces and sends a shock wave outward, which stalls. Neutrinos pouring out of the hot newborn neutron star reheat the gas behind the shock and relaunch it; they carry off ~99% of the ~3 × 10⁵³ erg released (C7). Hours to a day later the shock breaks out through the surface, and the envelope flies off at ~10,000 km/s with ~10⁵¹ erg of kinetic energy. A neutron star (C8) or a black hole is left behind.

**Thermonuclear explosion** (Type Ia). A carbon–oxygen white dwarf (B10) in a binary ignites carbon, either by growing toward ~1.4 M☉ or through a merger or a "double detonation" (L5). A degenerate star cannot expand to cool itself, so the burning runs away and destroys the whole star. Nothing is left behind.

Both kinds make radioactive nickel-56, whose decay powers most of the light we see (L1). Supernovae are classified by their spectra: Type II show hydrogen and Type I do not. Within Type I, Ia show silicon, Ib show helium, and Ic show neither.

**Suggested supernova arc:** L1 (what makes them shine) → L2 (the zoo and the progenitors) → L3 (the debris) → L4 (the elements) → L5 (Type Ia survivors). B11 is the before, C8 and B10 the after.

## L1. What makes a supernova shine? Radioactive nickel ★★☆

**Science question.** A supernova stays bright for months, long after the explosion itself is over. Where does that light come from? Can we read the radioactive decay of nickel-56 → cobalt-56 → iron-56 in real light curves, and weigh the freshly made nickel?

**Background.** Explosive burning makes ⁵⁶Ni (half-life 6.1 days), which decays to ⁵⁶Co (half-life 77.2 days) and then to stable ⁵⁶Fe. The gamma rays from these decays are trapped in the dense ejecta and turned into optical light. At peak, the luminosity roughly equals the radioactive heating at that moment ("Arnett's rule"), so the peak luminosity gives the nickel mass: ~0.6 M☉ for a typical Type Ia, a few hundredths of a solar mass for a Type II. Later the light curve fades at cobalt's decay rate, 0.98 mag per 100 days, while the ejecta still trap the gamma rays, and faster once they start to escape. For SN 2014J in M82, INTEGRAL detected the cobalt lines at 847 and 1238 keV directly: 0.6 ± 0.1 M☉ of nickel (Churazov et al. 2014). SN 1987A's radioactive tail gives ~0.07 M☉.

**Data.** Light curves from the Open Supernova Catalog's data repositories on GitHub (`astrocatalogs`, one JSON file per supernova). ZTF Bright Transient Survey light curves via ALeRCE (as in D6 and I10). Published bolometric light curves of SN 1987A and SN 2011fe. INTEGRAL line fluxes for SN 2014J (Churazov et al. 2014).

**Method + visualizations.**
- Bolometric light curves, built from multi-band photometry with a blackbody or bolometric correction, for a Type Ia (SN 2011fe, SN 2014J) and Type II supernovae (SN 1987A, SN 2023ixf).
- Overplot the decay energy of 0.6 and 0.07 M☉ of ⁵⁶Ni: the tails follow cobalt's slope.
- Arnett's rule for dozens of ZTF Type Ia: peak luminosity and rise time → nickel mass. Brighter Type Ia fade more slowly (the Phillips relation behind D1's standardization) because they made more nickel.
- SN 2014J's gamma-ray lines → nickel mass, compared with the optical estimate.
- Type II-P: the plateau (hydrogen recombining as the ejecta cool), then the drop onto the cobalt tail, whose height gives the nickel mass.

**Check against.** Churazov et al. 2014: SN 2014J, 0.6 ± 0.1 M☉; SN 1987A: ~0.07 M☉.

**Watch out for.** A bolometric luminosity needs ultraviolet and infrared coverage that most supernovae lack, so corrections are model-dependent. Dust in the host galaxy dims and reddens. Arnett's rule holds to ~10–20%. Late-time gamma-ray escape must be modelled. Luminosities scale with distance squared.

**Builds on.** D1 and D6 (Type Ia light curves), C7 (SN 1987A), B11 (MESA's `ccsn_IIp` light curve), L2.

## L2. A census of exploding stars ★★☆

**Science question.** Supernovae come in many types. How common is each type, how bright, in which galaxies, and from which stars? The ZTF Bright Transient Survey has classified almost every supernova brighter than 18.5 mag in the northern sky since 2018. Can we build the supernova zoo, and test the puzzle of the missing massive progenitors?

**Background.** Type Ia supernovae occur in all galaxies, including old ellipticals. Core-collapse supernovae occur only where stars are forming, because their progenitors live less than ~50 Myr. In a magnitude-limited survey most supernovae are Type Ia, because they are the brightest. In a fixed volume, core-collapse events outnumber them (Li et al. 2011). Archival images taken before the explosion reveal the progenitors of many nearby Type II-P supernovae: red supergiants of ~8–17 M☉. Smartt (2009, 2015) noticed that none was more luminous than log L/L☉ ≈ 5.1 (~18 M☉), although red supergiants exist up to ~25 M☉: about 13 such progenitors were expected and about one was found. This is the **red supergiant problem**. Perhaps the heaviest stars collapse straight into black holes without a bright explosion. N6946-BH1, a red supergiant that brightened briefly and then vanished (Gerke et al. 2015; Adams et al. 2017), is the best candidate for such a "failed supernova".

**Data.** The ZTF BTS explorer (`sites.astro.caltech.edu/ztf/bts/explorer.php`): classifications, redshifts, peak magnitudes and host galaxies as CSV. The survey is spectroscopically 93% complete to 18.5 mag (Perley et al. 2020). Light curves via ALeRCE. Host galaxy colours and masses from SDSS and the Legacy Survey (J3). Progenitor detections and limits from Smartt (2015)'s compilation.

**Method + visualizations.**
- Type fractions against survey depth, and peak absolute magnitude distributions by type (Ia, II, IIn, Ib/c, superluminous): the luminosity functions.
- Volume-corrected rates (1/V_max), and type fractions in a volume-limited subset compared with Li et al. 2011.
- Host galaxies on J3's star-formation-rate–mass plane: Type Ia everywhere, core-collapse only in star-forming galaxies.
- Type II-P plateau lengths and brightnesses → ejecta masses and explosion energies (Popov 1993 scaling).
- Progenitor masses against the initial mass function: how many progenitors above 18 M☉ should we have seen? The red supergiant problem as a probability.

**Check against.** Perley et al. 2020: Type Ia peak absolute magnitudes from about −16.5 to −20.5, and a Type Ia rate of ~2.35 × 10⁴ Gpc⁻³ yr⁻¹. Smartt 2015: ~13 high-mass progenitors expected, about one found.

**Watch out for.** Classifications from low-resolution spectra mix subtypes. Dust in their host discs dims core-collapse supernovae. Progenitor masses depend on stellar models and on bolometric corrections for cool supergiants. Dust around a progenitor can hide its luminosity, one proposed answer to the red supergiant problem.

**Builds on.** D6 (ZTF), I10 (TDEs, a contaminant), J3 (host galaxies), B11 (progenitor models), L1.

## L3. Rewind a supernova remnant: the Crab and Cassiopeia A ★★☆

**Science question.** Chinese astronomers recorded a "guest star" in 1054 where the Crab Nebula now sits. Cassiopeia A exploded around 1670 with no confirmed sighting. Can we measure how fast their debris flies outward, run it backward to the explosion date, and find the neutron stars they left behind?

**Background.** Debris filaments move outward at thousands of km/s, so at the Crab's 2 kpc they shift by ~0.1–0.2″ per year. Comparing images taken decades apart gives each filament's proper motion. Extrapolating backward to a common point gives a "convergence date". For the Crab it comes out between 1105 and 1130, later than 1054. The filaments have therefore sped up since the explosion, pushed by the nebula of particles and fields that the Crab pulsar's spin-down powers (C2). Cas A's optical knots, imaged over 50 years, converge on 1671.3 ± 0.9 (Thorstensen, Fesen & van den Bergh 2001). Light echoes, the supernova's light scattered off dust and arriving centuries late, gave spectra of the original explosion: Cas A was a Type IIb (Krause et al. 2008). Neutron stars receive kicks at birth. Ordinary pulsars move at 246 ± 22 km/s on average across the sky (Hobbs et al. 2005), and Cas A's neutron star sits away from the remnant's expansion centre.

**Data.**
- Crab: Digitized Sky Survey plates (POSS-I from the 1950s, POSS-II from the 1990s) from STScI's DSS server, and modern Pan-STARRS or Legacy Survey images (as in G1): a ~60-year baseline. The 2025 CFHT study (arXiv:2502.19632; 19,974 proper-motion vectors) for comparison.
- Cas A: Chandra ACIS observations from 2000 to the 2020s (Chandra Data Archive), or HST images of the optical knots (MAST).
- Green's catalogue of 310 Galactic supernova remnants (2024 October version, `mrao.cam.ac.uk/surveys/snrs/`). ATNF pulsar proper motions (C2) for kicks.

**Method + visualizations.**
- Register old and new images on Gaia stars and blink them: an animation of 60 years of expansion. Measure filament motions by cross-correlating small patches.
- Motions against positions → the expansion centre and convergence date; the Crab's acceleration from the gap to 1054.
- Cas A's expansion from Chandra: radial profiles of the forward shock at each epoch → expansion rate and age.
- The Galaxy's remnants from Green's catalogue on the sky; sizes against ages for remnants with known dates, against the Sedov–Taylor blast-wave law `R ∝ (E t² / ρ)^(1/5)`.
- Kicks: pulsar speeds across the sky from ATNF, and pulsars traced back to the centres of their remnants.

**Check against.** Nugent 1998: Crab convergence ~1130; the 2025 CFHT measurement: 1105.5 ± 0.5. Thorstensen et al. 2001: Cas A 1671.3 ± 0.9. Hobbs et al. 2005: 246 ± 22 km/s.

**Watch out for.** Photographic plates respond nonlinearly and used different emulsions and filters from modern cameras: measure positions, not brightnesses. Filaments change shape over decades. Plate astrometry needs a good reference frame: Gaia stars, with their own proper motions wound back to the plate epoch. Pulsar distances from dispersion measure are uncertain by tens of percent.

**Builds on.** G1 (images), C2 (the Crab pulsar, proper motions), B2 and H3 (Gaia astrometry), L2.

## L4. Where your atoms were made ★★☆

**Science question.** Apart from hydrogen and most of the helium (K1), every atom in your body was made in stars. Oxygen and calcium came mostly from massive stars that exploded quickly; much of the iron from Type Ia supernovae that took up to billions of years longer; europium and gold at least partly from neutron star mergers. Can we read this history in the chemistry of 900,000 stars, and weigh the heavy elements made in GW170817's kilonova?

**Background.** Core-collapse supernovae enrich gas with α elements (O, Mg, Si, Ca, Ti) within a few Myr of star formation, giving a high α/Fe ratio. Type Ia supernovae add mostly iron, after delays from ~0.1 Gyr to many Gyr. So as a galaxy ages, [α/Fe] stays high until Type Ia start to contribute, then falls. The resulting "knee" in the [α/Fe]–[Fe/H] diagram marks how quickly that galaxy formed stars. Elements heavier than iron are built by neutron captures: slow ones (the s-process, in asymptotic-giant-branch stars, making e.g. barium) and rapid ones (the r-process: europium, gold). The kilonova of GW170817 (AT2017gfo) is best fitted by ~0.08 M☉ of ejecta in three components, one rich in lanthanides (Villar et al. 2017).

**Data.** GALAH DR4 (Buder et al. 2025): up to 32 elements for 917,588 stars with Gaia astrometry (the `galah_dr4_allstar` catalogue, 723 MB, from the GALAH public data on Data Central). APOGEE DR17 as a cross-check. H3's Gaia–Sausage–Enceladus stars, debris of another galaxy with its own knee. Villar et al. 2017's unified kilonova data set (647 measurements from 46 instruments) and their public model code.

**Method + visualizations.**
- [Mg/Fe] against [Fe/H] for the disc: the high- and low-α sequences and the knee. The same for H3's Gaia–Sausage–Enceladus stars, whose knee sits at lower [Fe/H] because a small galaxy formed stars slowly.
- [Eu/Fe] (r-process) and [Ba/Fe] (s-process) against [Fe/H]. Large scatter at low metallicity points to rare, prolific sources such as mergers.
- A "periodic table of origins" in the style of Johnson (2019), coloured by what our stars show.
- The kilonova: fit a blue + red radioactive model to AT2017gfo's light curve → ejecta and lanthanide masses. Multiplied by the neutron star merger rate (C5), is that enough to explain the Milky Way's europium?

**Check against.** Villar et al. 2017: blue ~0.020 M☉ at 0.27c, purple ~0.047 M☉, red (lanthanide-rich) ~0.011 M☉. GALAH DR4's published abundance trends (Buder et al. 2025).

**Watch out for.** Abundance systematics (non-LTE effects, especially for O, Na and Al; blended europium lines) depend on stellar type, so compare like with like (giants with giants). The survey's footprint and selection shape the [Fe/H] distribution. Kilonova fits depend on assumed opacities, so ejecta masses carry factor-of-two systematics.

**Builds on.** B7 (spectroscopy), H3 (accreted galaxies), K1 (the Big Bang's share), C4 (GW170817), L1 and L2.

## L5. Runaway survivors of Type Ia supernovae ★★☆

**Science question.** What exactly explodes in a Type Ia supernova is still debated. One scenario leaves a smoking gun: when one white dwarf in a tight pair explodes, its companion flies off at its orbital speed, over 1,000 km/s, faster than any normal star. Can we find these survivors with Gaia?

**Background.** In the "D⁶" scenario (dynamically driven double-degenerate double-detonation), a white dwarf pulls helium from a close white dwarf companion. The helium shell detonates and triggers a carbon detonation in the core, giving a Type Ia below the Chandrasekhar mass. The companion survives and is released at its orbital speed. Shen et al. (2018) found three such runaways in Gaia DR2 moving at 1,000–3,000 km/s, one of them traceable to a known supernova remnant. El-Badry et al. (2023) found four more, with radial velocities up to −2,285 km/s and birth speeds of ~2,200–2,500 km/s: the fastest stars known in the Milky Way. They are hot, small and hydrogen-free, with carbon- or helium-dominated atmospheres. Slower runaways (≳600 km/s) may be partly burnt leftovers of weaker Type Iax explosions.

**Data.** Gaia DR3 astrometry and photometry; published spectroscopy and velocities (Shen et al. 2018; El-Badry et al. 2023); Green's remnant catalogue (L3); B10's white dwarf catalogue.

**Method + visualizations.**
- Repeat the search: Gaia sources with tangential velocities above ~600 km/s (from proper motion and parallax) that sit between the white dwarfs and the main sequence in the HR diagram. How many candidates are real, and how many are parallax errors?
- The known runaways in the HR diagram and on the sky, with their 3D velocities from the published radial velocities.
- Orbits traced back in a Milky Way potential (`galpy`) → birthplaces. Which pass near known remnants within their flight times?
- The orbital-speed argument: the ejection speed against donor mass for white dwarf pairs at contact → which donors give 2,500 km/s.
- Survivor numbers against the Type Ia rate (L2) → what fraction of Type Ia supernovae could be D⁶?

**Check against.** Shen et al. 2018 (three D⁶ stars); El-Badry et al. 2023 (radial velocities of −1,694 and −2,285 km/s; birth speeds ~2,200–2,500 km/s).

**Watch out for.** Most high-velocity "white dwarfs" in Gaia are spurious, because bad parallaxes inflate velocities; require significant parallaxes and clean astrometry. Without spectra, candidates stay candidates. Orbits traced far back are dominated by the distance uncertainty.

**Builds on.** B10 (white dwarfs), D1 (Type Ia), H3 (`galpy` orbits), L3 (remnants), I11 (Gaia quality cuts).

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
| A1 Exoplanet transit | **done** — [01-exoplanet-transit/](01-exoplanet-transit/). HAT-P-7 b and Kepler-10 b recovered from Kepler photometry; both periods correct to ~4 parts in 10⁶, radii bracketing the published values. Also recovers HAT-P-7 b's orbital phase curve and 48 ppm secondary eclipse. |
| A2 Exoplanet atlas | **done** — [02-exoplanet-population/](02-exoplanet-population/). 6,372 planets from the NASA Exoplanet Archive, with the archive's formula-filled masses (47% of the table) and radii excluded. Radius valley measured at 1.82 R⊕ (68%: 1.80–1.85), tilting with period at slope −0.10 vs Van Eylen et al.'s −0.09; 96% of small planets below it fit bare rock, 90% above it are too big for rock. Of 275 habitable-zone planets, 17 are rocky-sized — 15 around red dwarfs, none around a Sun-like star. |
| A3 Radial velocities | **done** — [20-radial-velocities/](20-radial-velocities/). 2,328 velocities from ELODIE, Lick, Keck and HARPS (1989–2019), fitted with our own multi-instrument periodogram, Keplerian fitter and a noise model that treats each night's exposures as one measurement. 51 Peg b: P = 4.2307860 ± 0.0000064 d, K = 55.36 ± 0.32 m/s, m sin i = 0.461 ± 0.015 M_Jup (true mass 0.463 ± 0.023 with the planet's own CRIRES orbital speed); the 21 archived ELODIE spectra from before the 1995 announcement alone give 4.2305 d. 55 Cnc: five planets found blind; e's true 0.7365-d period beats its 2.817-d alias (none of 500 simulations with the alias as truth reproduce the data), and the velocities predict the TESS transits 6,700 orbits later to 66 ± 26 min: 7.75 ± 0.37 M⊕, 1.87 ± 0.05 R⊕, 6.5 g/cm³, pure rock. Proxima b: K = 1.22 ± 0.16 m/s, m sin i = 1.05 ± 0.14 M⊕, absent from all seven activity indicators. GJ 581 "d" found blind at 66.9 d, then traced to half the star's 137-d rotation in H-alpha. Plus an offline fold-it-yourself page. |
| A4 Earth as an exoplanet | not started |
| A5 ALMA protoplanetary discs | not started |
| A6 Directly imaged orbits | not started |
| B1 HR diagram | **done** — [02-hr-diagram/](02-hr-diagram/). 215,844 Gaia DR3 stars within 100 pc. Twin-binary sequence 0.75 ± 0.02 mag above the main sequence (0.753 predicted); red clump M_G 0.51 ± 0.02; the Jao gap confirmed at 5.6σ; the white dwarf A/B split; Parenago's discontinuity at BP−RP 0.84. Plus an all-sky star-count map of 1.8 billion sources. |
| B2 3D solar neighborhood | **done** — [06-solar-neighborhood/](06-solar-neighborhood/). 300,338 stars within 100 pc: 299,525 from Gaia DR3 (72,701 spurious sources removed with the Rybizki et al. fidelity classifier, 71% of them toward the crowded inner Galaxy) plus 813 from Hipparcos that are too bright or too binary for Gaia. Interactive 3D map in true star colours, with SIMBAD names for every star within 20 pc. Red dwarfs are 74% of the stars within 20 pc and none of the 205 visible to the naked eye; only 46% of the naked-eye sky lies inside 100 pc. The Hyades found blind as the densest clump, at 46.8 pc (literature 47.5), and Coma Berenices at 85.5 pc (85.9); both beat an isotropic null (142 and 48 neighbours vs 25). Photometric classes match SIMBAD spectral letters for 96% of dwarfs; 99.3% of the Gaia Catalogue of Nearby Stars recovered. Of 477 objects known within 10 pc, 330 found: every missed star is in a binary, and 70 of the 79 missed brown dwarfs are too faint. 2.0 ± 0.5% more stars below the Sun than above (4.4σ), the sign expected for a Sun ~20 pc above the disc. |
| B3 Cluster hunting | **done** — [02-cluster-hunting/](02-cluster-hunting/). The Pleiades, Praesepe and M67 found blind as knots in Gaia DR3 proper motions (1,358 / 953 / 1,208 members); distances within 0.2% of Hunt & Reffert (2023), Pleiades at 135.0 ± 0.2 pc; members' main sequences 2.3–6.5× thinner than the field's; perspective and depth effects of a shared 3D velocity confirmed. Isochrone age fitting is the open extension. |
| B4 Cepheids / Leavitt law | not started |
| B5 Asteroseismology | not started |
| B6 Oldest stars | not started |
| B7 Spectroscopy | not started |
| B8 Cluster time-lapse | not started |
| B9 Eclipsing binaries | not started |
| B10 White dwarfs | not started |
| B11 MESA to core collapse | not started |
| B12 Constellations in motion | not started |
| C1 GW150914 chirp | **done** — [02-gravitational-waves/](02-gravitational-waves/). An hour of H1/L1 strain from GWOSC, processed from raw with our own whitening, matched filter and Q-transform (plus WAV audio). The discovery paper's Figure 1 reproduced at correlation 0.995 (H1) / 0.983 (L1); network SNR 23.8 (published 23.6–25.2); Hanford 7.10 ± 0.19 ms after Livingston (published 6.9), inverted as the arm geometry predicts. Zero-crossing chirp mass 35.2 M☉ (the LIGO authors' same method: ~37); calibrated by 461 injections into real noise, 31.8 M☉ (68%: 28.3–36.0) vs GWTC-1's 31.2 detector-frame — which leaves only black holes compact enough. Louder than all 4,015 time slides (187 days of background): FAR < 2/yr from one hour of data. |
| C2 Pulsar P–Ṗ | **done** — [03-pulsar-ppdot/](03-pulsar-ppdot/). 2,764 neutron stars from ATNF catalogue v2.8.1 (plus an interactive version). Recycling: 81% of millisecond pulsars have a companion vs 2.7% of ordinary ones, and spin tracks companion type (He white dwarf 4.1 ms → CO white dwarf 19.5 ms → neutron star 62.5 ms). 47 radio pulsars shine beyond the death valley. The 716 Hz record holder implies density > 7×10¹³ g/cm³ and radius < 14.4 km. Crab's characteristic age is 34% too old against SN 1054, implying a 17–20 ms birth spin. Pulsar heights grow from 29 to 374 pc above the plane with age. 40% of cluster pulsars appear to spin *up*. The sky map is partly survey footprint: one FAST survey found 32% of field millisecond pulsars. |
| C3 Fast radio bursts | **done** — [04-fast-radio-bursts/](04-fast-radio-bursts/). The Second CHIME/FRB Catalog (Jan 2026; 4,539 bursts from 3,641 sources) plus 80 host-localized FRBs, with ATNF pulsars as a control sample inside the Galaxy. Dispersion law measured on single bursts: delay ∝ ν^−2 to within 0.0003 in the index, timing residuals 0.05–0.1 ms over 4–10 s sweeps. The median FRB has 10× the DM that the whole Milky Way disk holds in its direction; 99.4% exceed disk plus a generous halo. FRBs mildly *avoid* the Galactic plane (12% within 10° vs 17% isotropic; pulsars 73%). Repeater bursts are 2.2× wider with half the bandwidth (p ~ 10⁻⁸). Macquart relation: diffuse baryon fraction f_d = 1.04 (+0.16/−0.12), so the "missing" baryons are there; run with Macquart et al.'s own set-up, the same data give Ω_b h₇₀ = 0.051, matching the published values. |
| C4 GW170817 | not started |
| C5 Black hole population | not started |
| C6 Gamma-ray bursts | not started |
| C7 Neutrinos | not started |
| C8 Neutron star interior | not started |
| C9 Binary pulsars | not started |
| D1 SN Hubble diagram | **done** — [05-supernova-cosmology/](05-supernova-cosmology/). All 1,701 Pantheon+ light curves (z = 0.001–2.26) fitted with our own full-covariance likelihood, which matches the Pantheon+ team's published chains to within 0.05σ on every parameter. Dark energy and present-day acceleration both at 7.0σ (Δχ² = 49; the 1998 papers had 3–4σ): Ωm = 0.295 ± 0.054, ΩΛ = 0.609 ± 0.081 with curvature free, flat Ωm = 0.332 ± 0.018, w = −0.91 ± 0.15 (consistent with a cosmological constant). Cepheid-calibrated H0 = 73.55 ± 1.02 km/s/Mpc (73.53 ± 1.05 from SH0ES's 277 supernovae with no cosmological model), 5.4σ from Planck's 67.4. Acceleration began at z = 0.59, 5.3 Gyr ago; implied age 12.5 Gyr. Also found that the release's per-supernova error column is ~1.4× the covariance diagonal (χ² 698 vs 1403 for 1588 dof); the figures use the covariance. |
| D2 SDSS cosmic web | **done** — [07-cosmic-web/](07-cosmic-web/). All 656,795 galaxies of the SDSS DR17 Main Galaxy Sample. The classic 4°-slice wedge map (33,899 galaxies) with the Sloan Great Wall outlined (densest connected stretch 254 Mpc end to end, RA 168°–208°, z = 0.081; Gott et al. traced 450 Mpc by eye), a 49 Mpc void, and a finger-of-God zoom on Abell 1750. The web tested against randomness on a volume-limited sample (134,975 galaxies with Mr − 5 log h < −20, verified complete to z = 0.106): ~5 billion pair counts give a projected correlation function within 96–100% of Zehavi et al. (2011) point by point at 0.4–10 h⁻¹ Mpc, with r₀ = 5.55 ± 0.24 h⁻¹ Mpc, γ = 1.85 (published 5.54 ± 0.14, 1.83); a galaxy is ~50× more likely than chance to have a neighbour 1 Mpc away. ξ(r_p, π) shows both fingers of God and Kaiser squashing. The red fraction rises from 51% to 81% with local density, at fixed luminosity too. 711 clusters found by friends-of-friends and named from the Abell catalogue; Coma's velocity dispersion 963 ± 28 km/s (Colless & Dunn: 1015 ± 59). Plus an offline interactive 3D map with one-slab-at-a-time views. |
| D3 Planck CMB | **done** — [09-planck-cmb/](09-planck-cmb/). COBE/FIRAS: the CMB is a 2.725 K blackbody (χ² 45 for 42 dof; Compton y = (−0.1 ± 2.8)×10⁻⁶). Planck 2018 SMICA half-mission maps (1.6 GB): the sky map, a pure-noise half-difference map, and our own noise-free cross-spectrum, which matches Planck's published spectrum to 0.05% over ℓ = 500–1000 once two traps were fixed: mask mode-coupling (+11% at the first trough with the usual area correction, removed exactly via correlation functions, verified on a simulated sky) and a pixel window that the beam file's note says is included but is not. Above ℓ ≈ 1300 a point-source excess, D ∝ ℓ² with 56 ± 2 µK² at ℓ = 3000 (Planck's 143 GHz fit: 51). Seven acoustic peaks measured (first at ℓ = 220.9 ± 2.3, 5736 ± 69 µK²), evenly spaced with phase φ = 0.277 ± 0.008: adiabatic, not isocurvature. CAMB fit to the TT spectrum: Ω_b h² = 0.02221 ± 0.00019, Ω_c h² = 0.1214 ± 0.0021, n_s = 0.963 ± 0.005, age 13.825 ± 0.033 Gyr, all within 0.7σ of Planck's own TT+lowE results. H₀ = 66.7 ± 0.9 km/s/Mpc: 5.1σ from D1's supernova value. |
| D4 Rotation curves | **done** — [08-rotation-curves/](08-rotation-curves/). All 175 SPARC galaxies; stars fixed at the stellar-population mass-to-light ratio (0.5 at 3.6 µm), nothing tuned. In every showcase galaxy the orbits need 4–9× the pull the visible matter supplies at the edge (NGC 3198: beyond 7 kpc it needs more mass than all its stars and gas combined, 6.8× by 44 kpc). Maximum-disc test: in 131 of 153 galaxies no mass-to-light ratio, however extreme, explains both centre and edge; the 5 a plausible one explains all have short curves. Radial acceleration relation reproduced point by point (2,692 of 2,693 published points within 0.015 dex): g† = 1.11 × 10⁻¹⁰ m/s² vs published 1.20 ± 0.24 (their own points through our fit: 1.12), scatter 0.133 dex. Baryonic Tully–Fisher slope 3.87 ± 0.09 (published 3.85 ± 0.09). NFW haloes match ΛCDM concentrations for the big spirals but not the dwarf DDO 154 (−3.1σ): the cusp–core problem. Plus an offline interactive explorer of all 175 galaxies with a mass-to-light slider. |
| D5 Evolving dark energy | not started |
| D6 Fresh supernovae | not started |
| D7 Is the universe expanding? | not started |
| D8 Cosmic chronometers | not started |
| D9 Lensed quasar time delays | not started |
| D10 ISW effect | not started |
| D11 Cluster counts | not started |
| E1 Kirkwood gaps | **done** — [10-kirkwood-gaps/](10-kirkwood-gaps/). 1,568,901 orbits from the MPC's MPCORB (29 Sep 2026) plus AstDyS proper elements and families. A blind search for dips with a wall on each side finds the gaps; only then are they named: the six clean interior gaps (3:1, 8:3, 5:2, 7:3, 9:4, 11:5) all sit 0.6–2.6 mAU (≤ 0.1%) from where Kepler's third law puts the resonances with Jupiter — chance alignment of the strongest four ~6 × 10⁻⁷. The 4:1 and 2:1 are the belt's edges; a weak dip at 2.415 AU matches Mars' 1:2 resonance, backed by a 9.6× pile-up of captured asteroids in proper elements. The 3:1 core is 92–94% empty at every size from > 9 km to < 1 km, so not a discovery artefact (Kirkwood's 91 asteroids of 1866 already had none where 4.5 were expected). Hildas 97.9% locked in the 3:2 (control 50.7%), their triangle's corners measured at −58°, +58°, +172° from Jupiter; Trojans L4/L5 = 1.47 ± 0.03 for H < 15 vs the bias-corrected 1.45 ± 0.05 (raw 1.65, inflated by discovery bias). Families 1.7–14.5× tighter in proper than osculating inclination (Koronis 14.5×). Labelled N-body simulation (3,000 real-asteroid-shaped particles, Jupiter + Saturn from Horizons, 1 Myr, REBOUND): 72% cleared from the 3:1 (median 24,000 yr) and 49% from the 5:2 vs 7% elsewhere; the 2:1 does not clear in 1 Myr. |
| E2 Orbit visualizer | **done** — [12-orbit-visualizer/](12-orbit-visualizer/). 20 objects from JPL Horizons (83 queries): planets, Pluto, Ceres, Halley, Hale-Bopp, the three interstellar objects and the five escaping spacecraft, plus an interactive 3D viewer for 1972–2050. Kepler's third law from the orbits alone: slope 1.49980, the Sun's GM to +0.4 ± 0.5 ppm of JPL's, and Jupiter's mass (Sun/Jupiter 1048 ± 16, JPL 1047.35); the ppm-level residuals of every planet match a Newtonian ring model, and the outer planets orbit the Sun plus Jupiter and Saturn (+0.14%). Halley sweeps equal areas to 0.3% while its speed changes 60-fold. Interstellar v∞ 26.39, 32.28, 57.97 km/s, radiants within 1′ of published (ʻOumuamua's 17′ shift is its non-gravitational push). Voyager 1 reaches one light-day from Earth at 2026-11-18 10:16:07.0 UTC — NASA's time to 5 ms — and never comes back inside. |
| E3 NEO approaches | **done** — [11-neo-close-approaches/](11-neo-close-approaches/). JPL close-approach, Sentry, fireball and Horizons data. In the 30 days to 29 Sep 2026, 165 known asteroids passed within 0.05 au but only 28 are predicted for the next 60: 146 of the 165 were unknown a month earlier. Near misses become impact rates through the d² law, which holds for large asteroids (k = 2.16, 1.97) and fails for small ones (k = 0.62), measuring discovery bias: one ≥1 km impact per 0.88 Myr (published 0.6 Myr), ≥140 m per 33,000 yr; against the fireball record, surveys catch ~27% of 10 m and ~3% of 3 m rocks passing inside the Moon's orbit. All 2,206 Sentry objects rate Torino 0 (recomputed, matching all 46,129 dated ratings). Apophis 2029: 38,011 km ± 3 km, path bent 27°, Aten → Apollo (323- to 423-day year). |
| E4 Meteor showers | not started |
| E5 Crater dating | not started |
| E6 Saturn's rings | not started |
| E7 Jupiter's weather | not started |
| E8 Eclipse shadows | not started |
| E9 Asteroids in 3D | not started |
| E10 Comet tails | not started |
| F1 Solar cycles | **done** — [15-solar-cycles/](15-solar-cycles/). The SILSO sunspot number 1749–2026 plus all 260,240 daily records of 45,376 sunspot groups since 1874 (Greenwich, then USAF/NOAA). Our own smoothing and cycle finder reproduce SILSO's official table: all 25 maxima exact, 23 of 25 minima (the other two are ties). Cycles last 9.0–13.6 yr (mean 11.03 ± 1.18) and peak at 81–285; Gnevyshev–Ohl holds in 9 of 11 pairs, broken exactly where the literature says (4→5, 22→23). Butterfly: daily spot areas identical to Hathaway's derived file on all 55,509 days; the zones drift as 28.0° × exp(−t/97 months) (Hathaway 2011: 28°, 90), and stronger cycles start higher (r = 0.69, p = 0.007, as Solanki et al. 2008 found, against Hathaway's "no dependence"). Differential rotation from 93,496 day-to-day moves: Ω = 14.554 − 2.80 sin²φ °/day (Balthasar et al. 1986, same plates: 14.551 − 2.87). Cycle 25 peaked at 160.9 in Oct 2024 (SILSO's value, reproduced), 40% above the 2019 forecast and 9 months early, yet only 15th of 25; its early rise had ruled out the forecast by mid-2023, and the late peak came from the south, 1.5 yr after the north. Maunder minimum: no spot seen in 34 of 71 years despite a median 822 observations each; reconstructions differ 20× on its depth. Plus an offline zoomable explorer. |
| F2 The Sun today | **done** — [16-sun-today/](16-sun-today/). SDO's near-real-time images of 29 Sep 2026 20:00 UTC in all ten AIA channels plus an HMI magnetogram, and the week before (57 magnetograms, 57 AIA frames, GOES-18 X-rays). 35 flares, 15 of them C-class, none above C3.1. The week's most energetic, a C2.4 on the west limb from NOAA 4536, reached 10.6 MK 4 min before its X-ray peak (emission measure 1.1×10⁴⁸ cm⁻³ a minute after), about 1.7×10²⁹ erg; its loops then cooled through AIA's channels in temperature order (131 Å +2 min, 94 Å +9, 335 Å +87, the cooler ones still rising at +104). Differential rotation from magnetogram cross-correlation: ω = 14.47 − 1.81 sin²φ − 2.05 sin⁴φ °/day (Komm et al. 1993: 14.42 − 2.00 − 2.09), 2.8% faster than the surface gas; the equator turns in 24.9 days (26.7 seen from Earth). Hale's law: 10 of 10 regions lead with Cycle 25's polarity. Limb darkening within 2.4% of Allen's law, giving 4,990 K at τ = 0.1 to 6,010 K at τ = 1. Coronal scale heights give 0.9/1.2/1.7 MK in 171/193/211 Å. Plus a week-long time-lapse and an offline explorer page; `run.py --epoch now` redoes it all for today. |
| F3 Gannon storm | not started |
| F4 Solar spectrum poster | not started |
| G1 Galaxy composites | **done** — [13-galaxy-composites/](13-galaxy-composites/). Twelve galaxies (M51, M101, M74, M63, M81, M64, M104, NGC 5866, NGC 4631, M82, the Mice, M87) built from calibrated SDSS g, r, i with one physical recipe: Lupton et al.'s arcsinh stretch fixed in surface-brightness units, and white = the Sun's measured colour (the survey viewer's default whitens a colour no star has). Checked first: the survey's own JPEGs reproduced (99.995% of 7.2 M values within ±1); 278 stars within 0.004 mag of Gaia's synthetic photometry; whole galaxies within 0.006 mag of Brown et al. 2014 (colours 0.003) once foreground stars are filled rather than dropped. True colour: a g, r, i → CIE XYZ matrix fitted to 129 galaxy spectra (leave-one-out error Δu′v′ 0.003, one just-noticeable difference) makes every galaxy a shade of white, 4,570 K (M104, M82) to 6,330 K (M101), matching their spectra to ≤ 0.007. Too faint for colour anyway: the brightest core (M104) reaches 0.08 cd/m², and at most 1.5% of any galaxy's area passes the 0.005 cd/m² cone threshold; no telescope can raise surface brightness. M51 on one grid from seven telescopes (FUV to 21 cm, ≤ 0.24″ from Gaia once stray SIP keywords are removed); Hubble's mosaic registered (0.84″ offset) and cross-calibrated against SDSS (−2/+9/+4% of the spectrum's prediction). Legacy Surveys DR9-north keeps only 36% of M51's z light. Plus an offline make-your-own-composite explorer. |
| G2 Meteorite map | **done** — [14-meteorite-landings/](14-meteorite-landings/). All 45,716 meteorites of the Meteoritical Society catalogue (NASA's 2013 snapshot; the old Socrata link is dead, the file moved to `data.nasa.gov/docs/legacy/`), on our own equal-area Equal Earth maps. Repaired first: 13,748 without usable coordinates (incl. 6,214 at a dummy 0°,0° and 223 Yamato meteorites whose zeroed latitude would put Antarctica in Kenya), one found on Mars, and every year before 1583 one too small (a Gregorian 1 January read back in the Julian calendar — all 8 match the arithmetic). 72% of located finds come from Antarctica (8% of the land); Oman yields 9,900 per million km², 150× the rest of the land. Against camera-network flux (58 falls ≥100 g per million km² per year), only 1 fall in ~1,500 is recorded, and recovery ∝ population density^0.59 (95%: 0.44–0.75) — 1 in 170 above 300 people/km², ~1 in 4,000 below 10; yet falls held at ~6/yr since 1900 while population grew. Evatt et al.'s iron table reproduced from the independent snapshot (falls 60/1,107 vs 60/1,097; Antarctica 0.74% vs 0.68%), and extended: individual finds are 40% metal vs 5% of falls; desert meteorites bought through dealers are 18% achondrites vs 4% where the find site was recorded. Plus an offline interactive map of all 31,893 located meteorites with filters. |
| H1 Sgr A* mass | **done** — [17-galactic-centre/](17-galactic-centre/). S2's orbit from both teams' raw measurements: 45 Keck and 145 NTT/VLT positions (1992–2018) plus 116 velocities — UCLA's Do et al. 2019 release (its page is gone; fetched from the Internet Archive) and Gillessen et al. 2017 via VizieR — fitted by MCMC with each team's reference frame free (one shared frame: χ² worse by 1,537). Each team's own analysis reproduced first: all 16 of Do et al.'s Newtonian posteriors within 0.15σ, their GR mass and distance within 0.04σ; Gillessen et al.'s VLT fit within 0.3σ on every number. Both teams together: M = 3.96 ± 0.04 ± 0.02 million M☉, R₀ = 7.93 ± 0.04 ± 0.03 kpc (Kepler by hand: 126.9 mas = 1,006 AU, P = 16.045 yr → 3.96 million suns). At GRAVITY's distance our banana gives 4.289 ± 0.010 million M☉ (GRAVITY: 4.297); their 4% larger R₀ is the known Keck–VLTI tension. Einstein's redshift detected: f = 0.82 ± 0.16 on Keck positions (Do et al.: 0.88), 0.70 ± 0.15 adding the VLT's, Newton excluded at 5σ; a Newtonian fit overstates M and R₀ by 5% and fakes Gillessen et al.'s 28 km/s "recession" (1.5 GM/ac = 17.5 km/s). A free precession fit claims 5× Einstein at 5σ — a false signal from 1990s speckle positions; adaptive-optics data alone cannot measure it. Pericentre 115 AU (16 light-hours, 1,470 horizon radii) at 7,590 km/s: > 5.5 × 10¹⁵ M☉/pc³, where the longest-lived dark cluster (Maoz 1998, reproduced) lasts 8 × 10⁵ yr. For I5: predicted shadow 51.2 ± 0.3 µas (EHT ring 51.8 ± 2.3). Next pericentre 5 June 2034. Plus an offline orbit explorer. |
| H2 21-cm rotation curve | **done** — [18-milky-way-rotation/](18-milky-way-rotation/). The HI4PI 21-cm survey: 25 cubes, of which only the \|v\| ≤ 350 km/s channels were fetched by HTTP byte range (3.8 of 6.6 GB). Terminal velocities from McClure-Griffiths & Dickey's two-error-function edge fit match the VLA survey at all 748 first-quadrant longitudes (median −0.24 km/s, rms 0.57), although HI4PI has a 16′ beam against the VLA's 1′. Their linear fits reproduced (v_circ 172.7 / 178.1 vs 172.5 / 177.7), as is the north–south asymmetry (−5.13 vs −5.2 km/s at \|l\| = 55–67°); an outward drift of the LSR explains none of it. Rotation curve from the tangent points with the LSRK convention handled explicitly (the Sun's angular speed, 30.32 km/s/kpc, from Sgr A*): 207–217 km/s at 2.6–3.6 kpc rising to 236–245 at 6–7.6 kpc, with ~10 km/s arm bumps; 63 of 1,442 longitudes with ambiguous edges dropped. The 3D orbits of 144 masers, computed from Reid et al.'s astrometry, lag rotation by 4.7 km/s and fall inward at 4.7 (published 4.3 ± 5.6, 6.0 ± 1.4); Perseus kinematic distances are 1.17× the parallaxes (up to 4.7×). Mass 1.06 × 10¹¹ M☉ inside the Sun's orbit; the Milky Way rotates faster than 89% of SPARC galaxies. Outer disc mapped face-on to 30 kpc after re-deriving LBH06's elliptical-orbit correction (with their constants α = 7.9, β = −1.06 vs published 8.67, −1.08; anticentre mismatch χ² 1137 → 0.6). HI scale length 4.05 kpc with the IAU constants (Kalberla & Dedes: 4.0 for co-rotating gas), 3.2 kpc with ours. Warp tilt W1 = 1.62 kpc at 20.5 kpc (LBH06: 1.74), maximum at φ ≈ 80–90°. The disc flares from 0.6 to 3.3 kpc thick. The plane cubes alone would have cut off the northern warp's top in 20–44% of cells beyond 16 kpc, so 7 cubes above the plane were added. |
| H3 Galactic archaeology | **done** — [19-galactic-archaeology/](19-galactic-archaeology/). Gaia DR3 only: 7.5 million rows from ESA's archive as small sky-patch queries, each all-sky sample checked against a server count. **GD-1** found blind among 3.8 million stars, with the ancient cluster M92 (same age, metallicity and distance) as the colour–magnitude filter and the stream's motion found in the data (up to 89σ): 3,464 members over 90° of sky; distance 7.8 → 11.6 kpc, gradient 0.0098 ± 0.0012 mag/deg (Koposov 2010: ~0.012), confirmed by the members' mean parallax (−0.004 ± 0.008 mas); gap G-20 78% deep (5.6σ), a 65% gap at −6.5°, the densest point at −14.5° (PWB2018's progenitor candidate: −13.5°); the spur 33 stars vs 2 in the mirror box; 65 of 66 spectroscopic members (Huang+2019) recovered. **Gaia-Sausage-Enceladus**: 74,938 halo stars within 2.5 kpc (0.41%); with no velocity cut at all, using Gaia's own RVS metallicities of cool giants, the halo's anisotropy jumps from β = 0.26–0.34 below [M/H] ≈ −2.2 to 0.74–0.83 above −2.0 (Belokurov+2018: 0.2–0.4 → ~0.9); the usual velocity cut finds only ~half of the round metal-poor halo. GSE median [M/H] −1.16 (Feuillet+2021: −1.17), 0.037 ± 0.010 dex poorer in α than the thick disc; its stars form the halo's blue sequence. Found on the way: half of GSP-Spec's "[M/H] < −1.66" stars are hot disc stars with wrong metallicities. **Phase spiral**: 2.5 million stars (Antoja+2018: 0.94 million); the snail shell in counts, v_R and v_φ once each star is compared with a phase-mixed disc; the disturbance dated to 250–560 Myr (median 360) in the McMillan 2017 potential, whose local disc matches the measured surface density (76 vs 68 ± 4 M☉/pc²); Antoja: 300–900 Myr. Plus an offline halo explorer. |
| H4 Sun's acceleration from quasars | not started |
| H5 Multiwavelength Milky Way | not started |
| H6 3D dust fly-through | not started |
| I1 M–σ relation | not started |
| I2 Reverberation mapping | not started |
| I3 Quasars | not started |
| I4 X-ray binaries | not started |
| I5 M87* shadow | not started |
| I6 Black hole ray tracing | not started |
| I7 Area theorem (GW250114) | not started |
| I8 Pulsar-timing background | not started |
| I9 Dark sirens | not started |
| I10 Tidal disruption events | not started |
| I11 Gaia dormant black holes | not started |
| I12 Microlensing black hole | not started |
| I13 ω Cen black hole | not started |
| I14 Megamasers | not started |
| I15 Reionization from quasar spectra | not started |
| I16 Superluminal jets | not started |
| I17 Quasar variability | not started |
| I18 Giant radio galaxies | not started |
| J1 Velocity fields | not started |
| J2 Galaxy morphology | not started |
| J3 Galaxy life cycle | not started |
| J4 Laniakea flows | not started |
| J5 Einstein rings | not started |
| J6 Local Group dwarfs | not started |
| J7 Andromeda collision | not started |
| J8 First galaxies (JWST) | not started |
| J9 Hubble Ultra Deep Field in 3D | not started |
| K1 Big Bang nucleosynthesis | not started |
| K2 Cluster masses | not started |
| K3 Bullet Cluster | not started |
| K4 MACHOs and primordial black holes | not started |
| K5 Milky Way dark halo | not started |
| K6 Direct detection | not started |
| L1 Supernova nickel | not started |
| L2 Supernova census | not started |
| L3 Remnant expansion | not started |
| L4 Origin of the elements | not started |
| L5 Type Ia survivors | not started |
